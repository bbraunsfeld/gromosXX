/**
 * @file perturbed_nonbonded_set.cc
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction_types.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_parameter.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"
#include "../../../interaction/nonbonded/pairlist/pairlist_algorithm.h"

#include "../../../math/periodicity.h"

#include "../../../interaction/nonbonded/interaction/storage.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_outerloop.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_outerloop.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_term.h"

#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_pair.h"

#include "../../../math/volume.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_set.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_set.h"

#include "../../../util/debug.h"

#include "nonbonded_set_interface.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::Perturbed_Nonbonded_Set
::Perturbed_Nonbonded_Set(Pairlist_Algorithm & pairlist_alg, Nonbonded_Parameter & param, 
			  int rank, int num_threads)
  : Nonbonded_Set(pairlist_alg, param, rank, num_threads),
    m_perturbed_outerloop(param),
    m_perturbed_pair(param)
{
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::Perturbed_Nonbonded_Set
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim)
{
  DEBUG(4, "Perturbed_Nonbonded_Set::calculate_interactions");

  // zero forces, energies, virial...
  m_storage.zero();

  // need to update pairlist?
  const bool pairlist_update = !(sim.steps() % sim.param().pairlist.skip_step);
  if(pairlist_update){
    DEBUG(6, "\tdoing longrange...");
    
    //====================
    // create a pairlist
    //====================
    
    // zero the longrange forces, energies, virial
    m_longrange_storage.zero();

    // parallelisation using STRIDE:
    // chargegroup based pairlist can only use this one!!!!
    // TODO:
    // move decision to pairlist!!!
    DEBUG(6, "create a pairlist");
    
    if (topo.perturbed_solute().atoms().size() > 0){
      m_pairlist_alg.update_perturbed(topo, conf, sim, 
				      pairlist(), perturbed_pairlist(),
				      m_rank, topo.num_atoms(), m_num_threads);
      DEBUG(10,"m_pairlist_alg.update_perturbed");
    }
    else{
      // assure it's empty
      assert(perturbed_pairlist().size() == topo.num_atoms());
      perturbed_pairlist().clear();
      
      
      m_pairlist_alg.update(topo, conf, sim, 
			    pairlist(),
			    m_rank, topo.num_atoms(), m_num_threads);
            DEBUG(10,"m_pairlist_alg.update");
    }

    // Print pairlist
    /*
    std::cout << pairlist().solute_long << std::endl;
    std::cout << pairlist().solute_short << std::endl;
    std::cout << pairlist().solvent_long << std::endl;
    std::cout << pairlist().solvent_short << std::endl;
    std::cout << perturbed_pairlist().solute_long << std::endl;
    std::cout << perturbed_pairlist().solute_short << std::endl;
    std::cout << perturbed_pairlist().solvent_long << std::endl;
    std::cout << perturbed_pairlist().solvent_short << std::endl;
   */  
  }

  if (sim.param().polarise.cos) {
    //===============
    // polarisation
    //===============

    // calculate explicit polarisation of the molecules
    DEBUG(6, "\texplicit polarisation");
    start_subtimer("explicit polarisation");
    if (topo.perturbed_solute().atoms().size() > 0) {
      m_perturbed_outerloop.perturbed_electric_field_outerloop(topo, conf, sim,
                                       m_pairlist, m_perturbed_pairlist,
				       m_storage, m_longrange_storage, m_rank);
    } else {
      m_outerloop.electric_field_outerloop(topo, conf, sim, m_pairlist, 
				       m_storage, m_longrange_storage, m_rank);
    }
    stop_subtimer("explicit polarisation");
  }  

  if (pairlist_update){
    start_subtimer("longrange");
    m_outerloop.lj_crf_outerloop(topo, conf, sim,
            m_pairlist.solute_long, m_pairlist.solvent_long,
            m_longrange_storage, true /*longrange!*/, m_pairlist_alg.timer(),
            m_rank == 0);
    if (topo.perturbed_solute().atoms().size() > 0){
      DEBUG(6, "\tperturbed long range");
      m_perturbed_outerloop.perturbed_lj_crf_outerloop(topo, conf, sim,
              m_perturbed_pairlist.solute_long,
              m_longrange_storage);
    }
    stop_subtimer("longrange");
  }

  // calculate forces / energies
  DEBUG(6, "\tshort range interactions");
  start_subtimer("shortrange");

  m_outerloop.lj_crf_outerloop(topo, conf, sim,
          m_pairlist.solute_short, m_pairlist.solvent_short,
          m_storage, false, m_pairlist_alg.timer(),
          m_rank == 0);

  if (topo.perturbed_solute().atoms().size() > 0){
    DEBUG(6, "\tperturbed short range");
    m_perturbed_outerloop.perturbed_lj_crf_outerloop(topo, conf, sim, 
						     m_perturbed_pairlist.solute_short,
						     m_storage);
  }

  stop_subtimer("shortrange");

  start_subtimer("1,4 interaction");
  m_outerloop.one_four_outerloop(topo, conf, sim, m_storage, m_rank, m_num_threads);
  stop_subtimer("1,4 interaction");
  start_subtimer("LJ exceptions");
  m_outerloop.lj_exception_outerloop(topo, conf, sim, m_storage,  m_rank, m_num_threads);
  stop_subtimer("LJ exceptions");

  if (sim.param().nonbonded.rf_excluded) {
    start_subtimer("RF excluded");
    DEBUG(6, "\tRF excluded interactions and self term");
    m_outerloop.RF_excluded_outerloop(topo, conf, sim, m_storage,  m_rank, m_num_threads);
    stop_subtimer("RF excluded");
  }
  
  // add 1,4 - interactions
  if (m_rank == 0){
    if (sim.param().polarise.cos) {
      start_subtimer("polarisation self-energy");
      if (topo.perturbed_solute().atoms().size()) {
        m_perturbed_outerloop.perturbed_self_energy_outerloop(topo, conf, sim, m_storage);
      } else {
        m_outerloop.self_energy_outerloop(topo, conf, sim, m_storage);
      }
      stop_subtimer("polarisation self-energy");
    } 
    
    DEBUG(6, "\t1,4 - interactions");
    if (topo.perturbed_solute().atoms().size() > 0){
      start_subtimer("1,4 interaction");
      DEBUG(6, "\tperturbed 1,4 - interactions");
      m_perturbed_outerloop.perturbed_one_four_outerloop(topo, conf, sim, 
							 m_storage);
      stop_subtimer("1,4 interaction");
    }
   
    
    // possibly do the RF contributions due to excluded atoms
    if(sim.param().nonbonded.rf_excluded){
      start_subtimer("RF excluded");
      DEBUG(6, "\tRF excluded interactions and self term");
      if (topo.perturbed_solute().atoms().size() > 0){
	DEBUG(6, "\tperturbed RF excluded interactions and self term");
	m_perturbed_outerloop.perturbed_RF_excluded_outerloop(topo, conf, sim,
							      m_storage);
      }
      stop_subtimer("RF excluded");
    }

    DEBUG(6, "\tperturbed pairs");
    start_subtimer("perturbed pairs");
    m_perturbed_pair.perturbed_pair_outerloop(topo, conf, sim, m_storage);
    stop_subtimer("perturbed pairs");
  }
  
  // add long-range force
  DEBUG(6, "\t(set) add long range forces");
  m_storage.force += m_longrange_storage.force;

  // DEBUG
  // std::cout << "total nb" << std::endl;
  // for(int i=0; i<topo.num_atoms(); ++i){
  // std::cout << "force " << i << " = " << math::v2s(m_storage.force(i)) << "\n";
  // }
  // DEBUG
  
  // and long-range energies
  DEBUG(6, "\t(set) add long range energies");
  const unsigned int lj_e_size = unsigned(m_storage.energies.lj_energy.size());
  
  for(unsigned int i = 0; i < lj_e_size; ++i){
    for(unsigned int j = 0; j < lj_e_size; ++j){
      m_storage.energies.lj_energy[i][j] += 
	m_longrange_storage.energies.lj_energy[i][j];
      m_storage.energies.crf_energy[i][j] += 
	m_longrange_storage.energies.crf_energy[i][j];
    }
  }

  //ANITA
  if (sim.param().precalclam.nr_lambdas){
    for(unsigned int i = 0; i < sim.param().precalclam.nr_lambdas; ++i){
      for(unsigned int j = 0; j < lj_e_size; ++j){
        for(unsigned int k = 0; k < lj_e_size; ++k){
          m_storage.energies.A_lj_energy[i][j][k] +=
            m_longrange_storage.energies.A_lj_energy[i][j][k];
          m_storage.energies.B_lj_energy[i][j][k] +=
            m_longrange_storage.energies.B_lj_energy[i][j][k];
          m_storage.energies.A_crf_energy[i][j][k] +=
            m_longrange_storage.energies.A_crf_energy[i][j][k];
          m_storage.energies.B_crf_energy[i][j][k] +=
            m_longrange_storage.energies.B_crf_energy[i][j][k];
        }
      }
    }
  }//ANITA

  // add longrange virial
  if (sim.param().pcouple.virial){
    DEBUG(7, "\t(set) add long range virial");

	m_storage.virial_tensor += m_longrange_storage.virial_tensor;
  }
  
  // and long-range energy lambda-derivatives
  DEBUG(7, "(set) add long-range lambda-derivatives");
  
  const unsigned int lj_size 
    = unsigned(m_storage.perturbed_energy_derivatives.lj_energy.size());
  
  for(unsigned int i = 0; i < lj_size; ++i){
    for(unsigned int j = 0; j < lj_size; ++j){
      
      assert(m_storage.perturbed_energy_derivatives.
	     lj_energy.size() > i);
      assert(m_storage.perturbed_energy_derivatives.
	     lj_energy[i].size() > j);
      assert(m_storage.perturbed_energy_derivatives.
	     lj_energy.size() > j);
      assert(m_storage.perturbed_energy_derivatives.
	     lj_energy[j].size() > i);
      assert(m_longrange_storage.perturbed_energy_derivatives.
	     lj_energy.size() > i);
      assert(m_longrange_storage.perturbed_energy_derivatives.
	     lj_energy[i].size() > j);
      assert(m_longrange_storage.perturbed_energy_derivatives.
	     lj_energy.size() > j);
      assert(m_longrange_storage.perturbed_energy_derivatives.
	     lj_energy[j].size() > i);
      
      m_storage.perturbed_energy_derivatives.lj_energy[i][j] += 
	m_longrange_storage.perturbed_energy_derivatives.lj_energy[i][j];
      
      m_storage.perturbed_energy_derivatives.crf_energy[i][j] += 
	m_longrange_storage.perturbed_energy_derivatives.crf_energy[i][j];
      
    }
  }
  
  //ANITA
  if (sim.param().precalclam.nr_lambdas){
    for(unsigned int i = 0; i < sim.param().precalclam.nr_lambdas; ++i){
      for(unsigned int j = 0; j < lj_e_size; ++j){
        for(unsigned int k = 0; k < lj_e_size; ++k){
          m_storage.perturbed_energy_derivatives.A_lj_energy[i][j][k] +=
            m_longrange_storage.perturbed_energy_derivatives.A_lj_energy[i][j][k];
          m_storage.perturbed_energy_derivatives.B_lj_energy[i][j][k] +=
            m_longrange_storage.perturbed_energy_derivatives.B_lj_energy[i][j][k];
          m_storage.perturbed_energy_derivatives.A_crf_energy[i][j][k] +=
            m_longrange_storage.perturbed_energy_derivatives.A_crf_energy[i][j][k];
          m_storage.perturbed_energy_derivatives.B_crf_energy[i][j][k] +=
            m_longrange_storage.perturbed_energy_derivatives.B_crf_energy[i][j][k];
        }
      }
    }
  }//ANITA


  DEBUG(7, "(set) calculate interactions done!");

  return 0;
}

int interaction::Perturbed_Nonbonded_Set::update_configuration
(
 topology::Topology const & topo,
 configuration::Configuration & conf,
 simulation::Simulation const & sim)
{
  const int ljs = conf.current().energies.lj_energy.size();
  configuration::Energy & pe = conf.current().perturbed_energy_derivatives;

  Nonbonded_Set::update_configuration(topo, conf, sim);
  
  for(int i = 0; i < ljs; ++i){
    for(int j = 0; j < ljs; ++j){
      
      assert(pe.lj_energy.size() > unsigned(i));
      assert(pe.lj_energy[i].size() > unsigned(j));

      assert(m_storage.perturbed_energy_derivatives.
	     lj_energy.size() > unsigned(i));
      assert(m_storage.perturbed_energy_derivatives.
	     lj_energy[i].size() > (unsigned(j)));
	  
      pe.lj_energy[i][j] += 
	m_storage.perturbed_energy_derivatives.
	lj_energy[i][j];
      pe.crf_energy[i][j] += 
	m_storage.perturbed_energy_derivatives.
	crf_energy[i][j];
    }
    pe.self_energy[i] += m_storage.perturbed_energy_derivatives.self_energy[i];
  }

  // ANITA
  const unsigned int nr_lambdas = unsigned(m_storage.energies.A_lj_total.size());

  for(unsigned int i=0; i < nr_lambdas; ++i) {
    for(int j=0; j < ljs; ++j) {
      for(int k=0; k < ljs; ++k) {
        pe.A_lj_energy[i][j][k] +=
            m_storage.perturbed_energy_derivatives.A_lj_energy[i][j][k];
        pe.B_lj_energy[i][j][k] +=
            m_storage.perturbed_energy_derivatives.B_lj_energy[i][j][k];
        pe.A_crf_energy[i][j][k] +=
            m_storage.perturbed_energy_derivatives.A_crf_energy[i][j][k];
        pe.B_crf_energy[i][j][k] +=
            m_storage.perturbed_energy_derivatives.B_crf_energy[i][j][k];
      }
    }
  } //ANITA

  
  return 0;
}

/**
 * calculate the hessian for a given atom.
 * this will be VERY SLOW !
 */
int interaction::Perturbed_Nonbonded_Set
::calculate_hessian(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    unsigned int atom_i, unsigned int atom_j,
		    math::Matrix & hessian){
  
  if (topo.is_perturbed(atom_i) ||
      topo.is_perturbed(atom_j)){
    assert(false);
    return -1;
  }

  return Nonbonded_Set::calculate_hessian(topo, conf, sim,
					  atom_i, atom_j, hessian);
}

int interaction::Perturbed_Nonbonded_Set
::init(topology::Topology const & topo,
       configuration::Configuration const & conf,
       simulation::Simulation const & sim,
       std::ostream & os,
       bool quiet)
{
  DEBUG(7, "Perturbed Nonbonded Set :: init");
  DEBUG(7, "sizes: " 
	<< unsigned(conf.current().perturbed_energy_derivatives.bond_energy.size())
	<< " - " 
	<< unsigned(conf.current().perturbed_energy_derivatives.kinetic_energy.size()));

  Nonbonded_Set::init(topo, conf, sim, os, quiet);
  
  m_storage.perturbed_energy_derivatives.resize
    (unsigned(conf.current().perturbed_energy_derivatives.bond_energy.size()),
     unsigned(conf.current().perturbed_energy_derivatives.kinetic_energy.size()),
     unsigned(sim.param().precalclam.nr_lambdas)); //ANITA

  m_longrange_storage.perturbed_energy_derivatives.resize
    (unsigned(conf.current().perturbed_energy_derivatives.bond_energy.size()),
     unsigned(conf.current().perturbed_energy_derivatives.kinetic_energy.size()),
     unsigned(sim.param().precalclam.nr_lambdas));  //ANITA

  perturbed_pairlist().resize(topo.num_atoms());

  return 0;
}

