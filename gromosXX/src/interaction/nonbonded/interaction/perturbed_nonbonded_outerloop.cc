/**
 * @file perturbed_nonbonded_outerloop.cc
 * (template) methods of Perturbed_Nonbonded_Outerloop.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction_types.h>
#include <math/periodicity.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/interaction/perturbed_nonbonded_term.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_innerloop.h>

#include <interaction/nonbonded/interaction/perturbed_nonbonded_outerloop.h>

#include <interaction/nonbonded/interaction_spec.h>

#include <util/debug.h>
#include <util/template_split.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::Perturbed_Nonbonded_Outerloop
::Perturbed_Nonbonded_Outerloop(Nonbonded_Parameter & nbp)
  : m_param(nbp)
{
}

//==================================================
// interaction loops
//==================================================

//==================================================
// the perturbed interaction (outer) loops
//==================================================

/**
 * helper function to calculate perturbed forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
void interaction::Perturbed_Nonbonded_Outerloop
::perturbed_lj_crf_outerloop(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     Pairlist const & pairlist,
			     Storage & storage)
{
  SPLIT_PERTURBATION(_perturbed_lj_crf_outerloop,
		     topo, conf, sim,
		     pairlist, storage);
}

template<typename t_perturbation_details>
void interaction::Perturbed_Nonbonded_Outerloop
::_perturbed_lj_crf_outerloop(topology::Topology & topo,
			      configuration::Configuration & conf,
			      simulation::Simulation & sim,
			      Pairlist const & pairlist,
			      Storage & storage)
{
  SPLIT_PERT_INNERLOOP(_split_perturbed_lj_crf_outerloop,
		       topo, conf, sim,
		       pairlist, storage);
}


template<typename t_interaction_spec, typename  t_perturbation_details>
void interaction::Perturbed_Nonbonded_Outerloop
::_split_perturbed_lj_crf_outerloop(topology::Topology & topo,
				    configuration::Configuration & conf,
				    simulation::Simulation & sim,
				    Pairlist const & pairlist,
				    Storage & storage)
{  
  DEBUG(7, "\tcalculate perturbed interactions");  

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Perturbed_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details> innerloop(m_param);
  innerloop.init(sim);
  innerloop.set_lambda(topo.lambda(), topo.lambda_exp());

  std::vector<unsigned int>::const_iterator j_it, j_to;
  unsigned int i;
  unsigned int size_i = unsigned(pairlist.size());

  if (t_interaction_spec::do_bekker){

    periodicity.recalc_shift_vectors();
    
    int pc;
    unsigned int j;

    // translate the atom j

    for(i=0; i < size_i; ++i){
      
      for(j_it = pairlist[i].begin(),
	    j_to = pairlist[i].end();
	  j_it != j_to;
	  ++j_it){
      
	pc = (*j_it >> 26);
	j = (*j_it & 67108863);
      
	DEBUG(10, "\tperturbed nonbonded_interaction: i " << i << " j " << j
	      << " pc " << pc);
      
	innerloop.perturbed_lj_crf_innerloop(topo, conf, i, j, storage, periodicity, pc);
      }
    }
  }
  else{

    DEBUG(9, "perturbed nonbonded_interaction: no grid based pairlist");
    for(i=0; i < size_i; ++i){
    
      for(j_it = pairlist[i].begin(),
	    j_to = pairlist[i].end();
	  j_it != j_to;
	  ++j_it){

	DEBUG(10, "\tperturbed nonbonded_interaction: i "
	      << i << " j " << *j_it);
	// printf("nb pair %d - %d\n", i, *j_it);
	
	// shortrange, therefore store in simulation.system()
	innerloop.perturbed_lj_crf_innerloop(topo, conf, i, *j_it, storage, periodicity);
      }
      
    }
  }

  DEBUG(7, "end of function perturbed nonbonded interaction");  
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
void interaction::Perturbed_Nonbonded_Outerloop
::perturbed_one_four_outerloop(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       Storage & storage)
{
  SPLIT_PERTURBATION(_perturbed_one_four_outerloop,
		     topo, conf, sim, storage);
}

template<typename t_perturbation_details>
void interaction::Perturbed_Nonbonded_Outerloop
::_perturbed_one_four_outerloop(topology::Topology & topo,
				configuration::Configuration & conf,
				simulation::Simulation & sim,
				Storage & storage)
{
  SPLIT_PERT_INNERLOOP(_split_perturbed_one_four_outerloop,
		       topo, conf, sim, storage);
}

template<typename t_interaction_spec, typename  t_perturbation_details>
void interaction::Perturbed_Nonbonded_Outerloop
::_split_perturbed_one_four_outerloop(topology::Topology & topo,
				      configuration::Configuration & conf,
				      simulation::Simulation & sim,
				      Storage & storage)
{
  DEBUG(7, "\tcalculate perturbed 1,4-interactions");
  
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Perturbed_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details> innerloop(m_param);
  innerloop.init(sim);
  innerloop.set_lambda(topo.lambda(), topo.lambda_exp());
  
  std::set<int>::const_iterator it, to;
  std::map<unsigned int, topology::Perturbed_Atom>::const_iterator 
    mit=topo.perturbed_solute().atoms().begin(), 
    mto=topo.perturbed_solute().atoms().end();
  
  for(; mit!=mto; ++mit){
    it = mit->second.one_four_pair().begin();
    to = mit->second.one_four_pair().end();
    
    for( ; it != to; ++it){

      innerloop.perturbed_one_four_interaction_innerloop
	(topo, conf, mit->second.sequence_number(), *it, periodicity);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}  


/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
void interaction::Perturbed_Nonbonded_Outerloop
::perturbed_RF_excluded_outerloop(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  Storage & storage)
{
  SPLIT_PERTURBATION(_perturbed_RF_excluded_outerloop,
		     topo, conf, sim, storage);
}

template<typename  t_perturbation_details>
void interaction::Perturbed_Nonbonded_Outerloop
::_perturbed_RF_excluded_outerloop(topology::Topology & topo,
				   configuration::Configuration & conf,
				   simulation::Simulation & sim,
				   Storage & storage)
{
  SPLIT_PERT_INNERLOOP(_split_perturbed_RF_excluded_outerloop,
		       topo, conf, sim, storage);
}

template<typename t_interaction_spec, typename  t_perturbation_details>
void interaction::Perturbed_Nonbonded_Outerloop
::_split_perturbed_RF_excluded_outerloop(topology::Topology & topo,
					 configuration::Configuration & conf,
					 simulation::Simulation & sim,
					 Storage & storage)
{

  DEBUG(7, "\tcalculate perturbed excluded RF interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Perturbed_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details> innerloop(m_param);
  innerloop.init(sim);
  innerloop.set_lambda(topo.lambda(), topo.lambda_exp());

  std::map<unsigned int, topology::Perturbed_Atom>::const_iterator
    mit=topo.perturbed_solute().atoms().begin(),
    mto=topo.perturbed_solute().atoms().end();

  DEBUG(9, "\tSize of perturbed atoms " 
	<< unsigned(topo.perturbed_solute().atoms().size()));
  
  for(; mit!=mto; ++mit){
    innerloop.perturbed_RF_excluded_interaction_innerloop(topo, conf, mit, periodicity);
  }
}
