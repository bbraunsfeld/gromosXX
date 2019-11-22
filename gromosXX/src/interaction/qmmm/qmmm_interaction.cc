/**
 * @file qmmm_interaction.cc
 * Implements QMMM interaction
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>

#include <interaction/qmmm/qm_zone.h>
#include <interaction/qmmm/qm_worker.h>
#include <interaction/qmmm/nonbonded/qmmm_nonbonded_set.h>

#include <interaction/interaction.h>
#include <interaction/qmmm/qmmm_interaction.h>

#include <util/debug.h>
#include <util/error.h>

#include <math/boundary_checks.h>
/*
#ifdef OMP
#include <omp.h>
#endif
*/
#ifdef XXMPI
#include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::QMMM_Interaction * interaction::QMMM_Interaction::qmmm_ptr = nullptr;

interaction::QMMM_Interaction::QMMM_Interaction() : Interaction("QMMM Interaction")
                                                  , m_parameter()
                                                  , m_set_size(1)
                                                  , m_rank(0)
                                                  , m_size(1)
                                                  , m_worker(nullptr)
                                                  , m_qm_zone(nullptr)
  {
#ifdef XXMPI
    m_rank = MPI::COMM_WORLD.Get_rank();
    m_size = MPI::COMM_WORLD.Get_size();
#endif
  qmmm_ptr = this;
}

interaction::QMMM_Interaction::~QMMM_Interaction() {
  for (unsigned int i = 0; i < m_qmmm_nonbonded_set.size(); ++i) {
    DEBUG(12, "deleting set " << i);
    delete m_qmmm_nonbonded_set[i];
    m_qmmm_nonbonded_set[i] = nullptr;
  }
  if (m_worker != nullptr)
    DEBUG(12, "deleting QM Worker ");
    delete m_worker;
  if (m_qm_zone != nullptr)
    DEBUG(12, "deleting QM Zone ");
    delete m_qm_zone;
}


// Create QM here and provide pointer for polarisation?
  /** We will create:
   * 1. QMMM interaction, that will run QM worker, provide QMMM pairlist and
   *    electric field for polarisation
   *    Will be run only on rank 0, but multiple CPUs
   * 2. QMMM LJ interaction, that will use QMMM pairlist to calculate QMMM classically
   *    Will be run on all MPI nodes except rank 0
   */





/**
 * This should be separate class AddRemove : public QM_Link
 * 
configuration::Configuration interaction::QMMM_Interaction::AddRemove(topology::Topology &topo,
                                                                      configuration::Configuration &conf,
                                                                      simulation::Simulation &sim) {
    // create a conf for capped system. Required for the evaluation of the forces on the external
    //point charges in the MOPAC worker
    //
    configuration::Configuration qmmm_conf = conf;
    qmmm_conf.current().force=0.0;
    math::Vec posCap,dR;
    unsigned int m1,q1;
    const double rch=0.109; //Carbon-Hydrogen bond length
    const bool verbose = false;
    unsigned int pi=0;
        for (std::vector< std::pair<unsigned int,unsigned int> >::const_iterator
                     it = topo.qm_mm_pair().begin(); it != topo.qm_mm_pair().end(); ++it,++pi  )
        {
            q1=it->first;
            m1=it->second;
            dR = conf.current().pos(m1) - conf.current().pos(q1);
            posCap = (rch/abs(dR)    )    * dR + conf.current().pos(q1) ;
            qmmm_conf.current().pos(m1) = posCap;
        }
    return qmmm_conf;
}
*/

/**
 * This should be separate class AddRemove2 : public QM_Link
 *
int interaction::QMMM_Interaction::AddRemove2(topology::Topology &topo,
                                                                      configuration::Configuration &conf,
                                                                      simulation::Simulation &sim) {
    math::Vec posCap,dQM_MM,unit_MM_QM;
    math::Vec FQM,FMM;
    double d_quot=0;
    double d_cart=0;
    unsigned int m1,q1;
    const double rch=0.109; //Carbon-Hydrogen bond length
    const bool verbose = false;
    unsigned int pi=0;
    for (std::vector< std::pair<unsigned int,unsigned int> >::const_iterator
                 it = topo.qm_mm_pair().begin(); it != topo.qm_mm_pair().end(); ++it,++pi  )
    {
        q1=it->first;
        m1=it->second;

       // std::cout << "size qmpair" << topo.qm_mm_pair().size()<<  std::endl;
        dQM_MM = conf.current().pos(m1) - conf.current().pos(q1);
        unit_MM_QM=dQM_MM/abs(dQM_MM);
        posCap = (rch/abs(dQM_MM)    )    * dQM_MM + conf.current().pos(q1) ;
        d_quot=rch/abs(dQM_MM);
        //std::cout << "pi" << pi << "  abs(unit_mm)" << abs(unit_MM_QM) << std::endl;
        for (unsigned int i = 0; i < 3; i++) {
            d_cart=dQM_MM[i];
            FQM[i]=storage.force(q1)[i]+stor_link.force(pi)[i]*((1.-d_quot)+(d_cart*unit_MM_QM[i]*d_quot/abs(dQM_MM)));
            FMM[i]=stor_link.force(pi)[i]*((d_quot)-(d_cart*unit_MM_QM[i]*d_quot/abs(dQM_MM)));
            //FMM=conf.current().force(m1)[i]+stor_link.force(pi)[i]*((d_quot)-(d_cart*unit_MM_QM[i]*d_quot/abs(dQM_MM)));
        }
       // std::cout << "FQM " << stor_link.force(pi)[0] << "  " <<  stor_link.force(pi)[1] << "  "
        //          << stor_link.force(pi)[2] << std::endl;
        //  std::cout << "FQM " << FQM[0] << "  " << FQM[1] << "  " << FQM[2] << std::endl;
     //   std::cout << "FMM " << FMM[0] << "  " << FMM[1] << "  " << FMM[2] << std::endl;
 //       conf.current().force(q1)+=FQM;
//        conf.current().force(m1)+=FMM;
        storage.force(q1)=FQM; //otherwise doublecounting when QM-forces are added to conf.current().force
        storage.force(m1)=FMM; //otherwise doublecounting when QM-forces are added to conf.current().force
    }

}

*/
int interaction::QMMM_Interaction::scf_step(topology::Topology& topo,
                                            configuration::Configuration& conf,
                                            simulation::Simulation& sim) {
  // Update only COS in QM zone
  m_timer.start();
  DEBUG(15,"Updating COS in QM zone");
  m_timer.start("QM zone update");
  m_qm_zone->update_cos(topo, conf, sim);
  m_timer.stop("QM zone update");

  int err = 0;
  m_timer.start(m_worker->name());
  err = m_worker->run_QM(topo, conf, sim, *m_qm_zone);
  m_timer.stop(m_worker->name());
  if (err) return err;
  m_timer.stop();
  return 0;
}

void interaction::QMMM_Interaction::write_qm_data(topology::Topology& topo,
                                                  configuration::Configuration& conf,
                                                  const simulation::Simulation& sim) {
  // If geometry minimisation within QM is requested
  /*if (minimise) {
    //qm_zone->write_pos(conf.current().pos);
  }
  */
  DEBUG(15,"Writing QM data");
  m_timer.start("writing QM results");
  m_qm_zone->write_force(conf.current().force);
  if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical)
    m_qm_zone->write_charge(topo.charge());
  conf.current().energies.qm_total = m_qm_zone->QM_energy();
  m_timer.stop("writing QM results");
}


int interaction::QMMM_Interaction::calculate_interactions(topology::Topology& topo,
                                                          configuration::Configuration& conf,
                                                          simulation::Simulation& sim) {
  DEBUG(4, "QMMM_Interaction::calculate_interactions");
  int err = 0;

  // Do QMMM only on master (MPI)
  if (m_rank == 0) {
    m_timer.start();
    // Update QM Zone
    DEBUG(15,"Updating QM zone");
    m_timer.start("QM zone update");
    err = m_qm_zone->update(topo, conf, sim);
    m_timer.stop("QM zone update");
    if (err) return err;
    
    m_timer.start(m_worker->name());
    err = m_worker->run_QM(topo, conf, sim, *m_qm_zone);
    m_timer.stop(m_worker->name());
    if (err) return err;
    
    if (sim.param().qmmm.qmmm != simulation::qmmm_polarisable) {
      // in polarisable embedding, we will write the data after electric field evaluation
      this->write_qm_data(topo, conf, sim);
    }

    if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical
        || sim.param().qmmm.qm_lj) {
      this->calculate_nonbonded(topo, conf, sim);
    }
    m_timer.stop();
  }
  return 0;
}

int interaction::QMMM_Interaction::calculate_nonbonded(topology::Topology& topo,
                                                       configuration::Configuration& conf,
                                                       simulation::Simulation& sim)
  {
  DEBUG(4, "QMMM_Interaction::calculate_nonbonded");

  // allow multistep - calculate nonbonded only every nth step
  //int steps = sim.param().multistep.steps;
  //if (steps == 0) steps = 1;

  //if ((sim.steps() % steps) == 0) {
    for (unsigned i = 0; i < m_set_size; ++i) {
      if(m_qmmm_nonbonded_set[i]->calculate_interactions(topo, conf, sim)) {
        m_timer.stop();
	      return 1;
      }
    }
  //}
  DEBUG(6, "sets are done, adding things up...");
  this->store_set_data(topo, conf, sim);
  
  if (sim.param().pairlist.print &&
      (!(sim.steps() % sim.param().pairlist.skip_step))) {
    DEBUG(7, "print QM-MM pairlist...");
    this->print_pairlist(topo);
  }

  DEBUG(6, "QMMM_Interaction::calculate_nonbonded done");
  return 0;
}

void interaction::QMMM_Interaction::get_electric_field(const simulation::Simulation& sim
                                                     , math::VArray & electric_field) {
  // Write electric field
  m_timer.start();
  m_timer.start("Electric field writing");
  m_qm_zone->electric_field(sim, electric_field);
  m_timer.stop("Electric field writing");
  m_timer.stop();
}

int interaction::QMMM_Interaction::init(topology::Topology& topo,
            configuration::Configuration& conf,
            simulation::Simulation& sim,
            std::ostream& os,
            bool quiet) {
  if (!quiet)
    os << "QMMM INTERACTION\n";
    // Initial checks
  if (sim.param().qmmm.cutoff > 0.0 && 
          !math::boundary_check_cutoff(conf.current().box, sim.param().boundary.boundary,
          sim.param().qmmm.cutoff)) {
    io::messages.add("The RCUTQ cutoff is too large for the size of the "
            "computational box.", "QMMM_Interaction", io::message::error);
    return 1;
  }

  if (m_rank == 0) {
    DEBUG(15,"Creating QM Worker");
    m_worker = interaction::QM_Worker::get_instance(sim);
    if (m_worker == nullptr ||
        m_worker->init(sim))
      return 1;
    DEBUG(15,"QM Worker initialized");
    
    // Create QM_Zone
    m_qm_zone = new interaction::QM_Zone;
    DEBUG(15,"QM Zone created");

    if (m_qm_zone->init(topo, conf, sim)) return 1;
    DEBUG(15,"QM Zone built");
  }
  if (!quiet) {
    switch (sim.param().qmmm.qmmm) {
      case simulation::qmmm_mechanical:
        os << "\tmechanical";
        break;
      case simulation::qmmm_electrostatic:
        os << "\telectrostatic";
        break;
      case simulation::qmmm_polarisable:
        os << "\tpolarisable";
        break;
      default:
        os << "\tunknown";
    }
    os << " embedding scheme" << std::endl;
    os << "\tusing external ";
    switch (sim.param().qmmm.software) {
      case simulation::qm_mndo:
        os << "MNDO";
        break;
      case simulation::qm_dftb:
        os << "DFTB";
        break;
      case simulation::qm_mopac:
        os << "MOPAC";
        break;
      case simulation::qm_turbomole:
        os << "Turbomole";
        break;
      default:
        os << "unknown";
        break;
    }
    os << " program package" << std::endl;
    os.precision(3);
    if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical) {
      os << "\tQM-MM interactions will be treated classically" << std::endl
         << "\tusing standard cutoffs (RCUTP, RCUTL)" << std::endl
         << "\tRCUTQM value (" << sim.param().qmmm.cutoff << ") ignored";
    }
    else if (sim.param().qmmm.cutoff == 0.0) {
      os << "\tincluding all MM atoms in QM-MM interaction calculation";
    }
    else {
      if (sim.param().qmmm.atomic_cutoff) os << "\tatom-based cutoff";
      else os << "\tchargegroup-based cutoff";
    os << ": " << std::fabs(sim.param().qmmm.cutoff) << " nm";
    }
    os << std::endl;

    if (topo.qmmm_link().size()) {
      os << "\tusing link-atom scheme with capping atom" << std::endl
         << "\tdistance between QM link atom and capping atom: "
         << sim.param().qmmm.cap_length << std::endl
         << "\tQM-MM links:" << std::endl;
      for (std::set< std::pair<unsigned,unsigned> >::const_iterator
            it = topo.qmmm_link().begin()
          , to = topo.qmmm_link().end()
          ; it != to; ++it) {
        os << "\t\t" << it->first << "-" << it->second << std::endl;
      }
    }

    if (sim.param().qmmm.qm_lj)
      os << "\tLJ interactions between QM atoms activated" << std::endl;
    else
      os << "\tno LJ interactions between QM atoms" << std::endl;

    if (sim.param().qmmm.mm_scale > 0.0) {
      os << "\tMM point charges will be scaled using 2/pi atan(s*|R|) with s = " <<
        sim.param().qmmm.mm_scale << std::endl;
      os << "\t\t|R| is distance to closest QM atom" << std::endl;
    }

    os << "\tworker: " << m_worker->name() << std::endl
    // Some information on QM worker?
      << "\tQM zone: " << std::endl
      << "\t\tnumber of QM atoms: \t" << m_qm_zone->qm.size() << std::endl
      << "\t\tnumber of QM-MM links: \t" << m_qm_zone->link.size() << std::endl;
  }

  // Remove relevant bonded terms from topo
  DEBUG(15, "Removing bonded terms");
  this->remove_bonded_terms(topo, os, quiet);

  // Remove exclusions containing QM-QM and QM-MM interactions
  DEBUG(15, "Removing QM-QM and QM-MM exclusions");
  this->remove_exclusions(topo, sim, os, quiet);

  // Create nonbonded set for LJ interactions
  if (sim.param().force.nonbonded_vdw
      && (sim.param().qmmm.qmmm > simulation::qmmm_mechanical // Will do LJ between QM and MM
        || sim.param().qmmm.qm_lj)) // Will do LJ between QM atoms
    {
    this->init_nonbonded(topo, conf, sim, os, quiet);
  }

  if (!quiet)
    os << "END" << std::endl;
  DEBUG(9, "QMMM init done");
  return 0;
}

int interaction::QMMM_Interaction::init_nonbonded(topology::Topology& topo,
                                                  configuration::Configuration& conf,
                                                  simulation::Simulation& sim,
                                                  std::ostream& os,
                                                  bool quiet) {
  DEBUG(9, "QMMM nonbonded init");
  if (!quiet)
    os << "\tQMMM nonbonded interaction:\n";

  if (sim.param().multicell.multicell) {
    io::messages.add("MULTICELL is not allowed with QMMM"
                    , "QMMM_Interaction", io::message::error);
    return 1;
  }

  m_qmmm_nonbonded_set.clear();

// OpenMP and MPI parallelization if necessary
/*
#ifdef OMP
  m_set_size *= omp_get_num_threads();
#endif
#ifdef XXMPI
  m_set_size *= MPI::COMM_WORLD.Get_size();*/
  for (unsigned i = 0; i < m_set_size; ++i) {
    m_qmmm_nonbonded_set.push_back(
          new QMMM_Nonbonded_Set(*(this->m_qm_zone), this->m_timer
                                , m_parameter, i, m_set_size));
  }

  if (!quiet)
    os << "\t\tcreated " << m_qmmm_nonbonded_set.size() << " set"
       <<  (m_qmmm_nonbonded_set.size() > 1 ? "s" : "") << std::endl;

  std::vector<QMMM_Nonbonded_Set *>::iterator
    it = m_qmmm_nonbonded_set.begin(),
    to = m_qmmm_nonbonded_set.end();

  bool q = quiet;
  for (; it != to; ++it) {
      (*it)->init(topo, conf, sim, os, q);
    // only print first time...
    q = true;
  }
  DEBUG(9, "QMMM nonbonded init done");
  return 0;
}

void interaction::QMMM_Interaction::remove_bonded_terms(
                                topology::Topology& topo
                              , std::ostream& os
                              , bool quiet)
  {
  // Remove bonded terms that will be defined by QM
  // Definitions to simplify code
  typedef std::vector<topology::two_body_term_struct> twoBodyVec;
  typedef std::vector<topology::three_body_term_struct> threeBodyVec;
  typedef std::vector<topology::four_body_term_struct> fourBodyVec;
  typedef std::vector<topology::eight_body_term_struct> eightBodyVec;
  twoBodyVec& bonds = topo.solute().bonds();
  twoBodyVec& cgbonds = topo.solute().cgbonds();
  threeBodyVec& angles = topo.solute().angles();
  fourBodyVec& improper_dihedrals = topo.solute().improper_dihedrals();
  fourBodyVec& dihedrals = topo.solute().dihedrals();
  eightBodyVec& crossdihedrals = topo.solute().crossdihedrals();

  // Simple lambda function to check if the QM-MM link was defined
  auto are_linked = [&](unsigned qmi, unsigned mmi)->bool {
    return topo.qmmm_link().count(std::make_pair( qmi, mmi ));
  };
  if (!quiet)
    os << "\tterms removed from topology:\n";
  // Counter for neat printing
  unsigned count = 0;
  if (!quiet)
    os << "\t\tbonds:\n";
  // Delete QM-QM bonded terms and check, if QM-MM links were properly defined
  for (twoBodyVec::iterator b_it = bonds.begin(); b_it != bonds.end(); ) {
    // If QM-QM
    if (topo.is_qm(b_it->i) && topo.is_qm(b_it->j)) {
      if (!quiet)
        if (count == 0) os << "\t\t";
        os << (b_it->i + 1) << "-" << (b_it->j + 1) << " ";
        if (++count == 8) {
          os << "\n";
          count = 0;
        }
      DEBUG(4,"Erased bond: " << b_it->i + 1 << " - " << b_it->j + 1);
      b_it = bonds.erase(b_it);
    }
    // If QM-MM or MM-QM and not in QM-MM link
    else if ((topo.is_qm(b_it->i)
            && !are_linked( b_it->i, b_it->j ))
          ||
            (topo.is_qm(b_it->j)
            && !are_linked( b_it->j, b_it->i ))
          )
      {
      std::ostringstream msg;
      msg << "Bonded interaction across QM-MM boundary, but no QMMM link defined for"
          << " atoms " << b_it->i << " and " << b_it->j;
      io::messages.add(msg.str(), "QMMM_interaction", io::message::warning);
      ++b_it;
    }
    else {
      ++b_it;
    }
  }

  if (!quiet && cgbonds.size())
    os << "\t\tcoarse-grained bonds:\n";
  // Delete coarse-grained bonds between QM-QM
  for (twoBodyVec::iterator b_it = cgbonds.begin(); b_it != cgbonds.end(); ) {
    // If QM-QM
    if (topo.is_qm(b_it->i) && topo.is_qm(b_it->j)) {
      if (!quiet)
        os << "\t\t" << (b_it->i + 1) << " " << (b_it->j + 1) << "\n";
      DEBUG(4,"Erased coarse-grained bond: " << b_it->i + 1 << "-" << b_it->j + 1);
      b_it = cgbonds.erase(b_it);
      continue;
    } else ++b_it;
  }

  if (!quiet) {
    os << "\n\t\tbond angles:\n";
    count = 0;
  }
  // Delete (QM/MM)-QM-(QM/MM) bond-angle terms and check, if QM-MM links were properly defined
  for (threeBodyVec::iterator a_it = angles.begin(); a_it != angles.end(); ) {
    const unsigned i = a_it->i
                 , j = a_it->j
                 , k = a_it->k;
    const unsigned qm_count = topo.is_qm(i) + topo.is_qm(j) + topo.is_qm(k);
    if (qm_count == 0) { ++a_it; continue; }// MM-MM-MM
    if (topo.is_qm(j)) {
      a_it = angles.erase(a_it);
      if (!quiet) {
        if (count == 0) os << "\t\t";
        os << (i+1) << "-" << (j+1) << "-" << (k+1) << " ";
        if (++count == 5) {
          os << "\n";
          count = 0;
        }
      }
      DEBUG(4,"Erased bond-angle bending term: " << i+1 << "-" << j+1 << "-" << k+1);
      if (qm_count == 3) continue; // QM-QM-QM, nothing to check
    } else ++a_it;
    // Loop over neighbouring pairs and check, if QM-MM links were properly defined
    const std::vector<unsigned> indices = {i,j,k};
    for (std::vector<unsigned>::const_iterator
        i_it = indices.begin(), j_it = i_it + 1, i_to = indices.end() - 1
        ; i_it != i_to; ++i_it, ++j_it) {
      if ((   topo.is_qm(*i_it)
          && !topo.is_qm(*j_it)
          && !are_linked(*i_it, *j_it))
        ||
          (  !topo.is_qm(*i_it)
          &&  topo.is_qm(*j_it)
          && !are_linked(*j_it, *i_it)
        ))
        {
        std::ostringstream msg;
        msg << "Bond-angle bending interaction across QM-MM boundary, "
          << "but no QMMM link defined between atoms "
          << *i_it + 1 << " and " << *j_it + 1;
        io::messages.add(msg.str(), "QMMM_interaction", io::message::warning);
      }
    }
  }

  if (!quiet) {
    os << "\n\t\tdihedral angles:\n";
    count = 0;
  }
  // Delete (MM/QM)-QM-QM-(MM/QM) terms - dihedrals
  for (fourBodyVec::iterator d_it = dihedrals.begin(); d_it != dihedrals.end(); ) {
    const unsigned i = d_it->i
                 , j = d_it->j
                 , k = d_it->k
                 , l = d_it->l;
    const unsigned qm_count = topo.is_qm(i) + topo.is_qm(j) + topo.is_qm(k) + topo.is_qm(l);
    if (qm_count == 0) { ++d_it; continue; }
    if ((topo.is_qm(j) && topo.is_qm(k))) {
      d_it = dihedrals.erase(d_it);
      if (!quiet) {
        if (count == 0) os << "\t\t";
        os << (i+1) << "-" << (j+1) << "-" << (k+1) << "-" << (l+1) << " ";
        if (++count == 4) {
          os << "\n";
          count = 0;
        }
      }
      DEBUG(4,"Erased dihedral angle torsion term: "
            << i+1 << "-" << j+1 << "-" << k+1 << "-" << l+1);
      if (qm_count == 4) continue;
    } else ++d_it;
    // loop over neighbouring pairs and check, if QM-MM links were properly defined
    const std::vector<unsigned> indices = {i,j,k,l};
    for (std::vector<unsigned>::const_iterator
        i_it = indices.begin(), j_it = i_it + 1, to = indices.end()
        ; j_it != to; ++i_it, ++j_it) {
      if ((   topo.is_qm(*i_it)
          && !topo.is_qm(*j_it)
          && !are_linked(*i_it, *j_it))
        ||
          (  !topo.is_qm(*i_it)
          &&  topo.is_qm(*j_it)
          && !are_linked(*j_it, *i_it)
        ))
        {
        std::ostringstream msg;
        msg << "Dihedral angle torsion interaction across QM-MM boundary, "
            << "but no QMMM link defined between atoms "
            << *i_it + 1 << " and " << *j_it + 1;
        io::messages.add(msg.str(), "QMMM_interaction", io::message::warning);
      }
    }
  }
  DEBUG(15,"Dihedral angle torsion terms done");
  
  if (!quiet) {
    os << "\n\t\timproper dihedral angles:\n";
    count = 0;
  }
  // Delete QM-(MM/QM)-(MM/QM)-(MM/QM) improper dihedral
  for (fourBodyVec::iterator
      d_it = improper_dihedrals.begin(), d_to; d_it != improper_dihedrals.end(); ) {
    const unsigned i = d_it->i
                 , j = d_it->j
                 , k = d_it->k
                 , l = d_it->l;
    const unsigned qm_count = topo.is_qm(i) + topo.is_qm(j) + topo.is_qm(k) + topo.is_qm(l);
    if (qm_count == 0) { ++d_it; continue; }
    if ( topo.is_qm(d_it->i)
      /*&& topo.is_qm(d_it->j)
      && topo.is_qm(d_it->k)
      && topo.is_qm(d_it->l)*/ )
      {
      DEBUG(4,"Erased improper dihedral angle bending term: "
            << i + 1 << "-" << j + 1 << "-"
            << k + 1 << "-" << l + 1);
      d_it = improper_dihedrals.erase(d_it);
      if (!quiet) {
        if (count == 0) os << "\t\t";
        os << (i+1) << "-" << (j+1) << "-" << (k+1) << "-" << (l+1) << " ";
        if (++count == 4) {
          os << "\n";
          count = 0;
        }
      }
      if (qm_count == 4) continue;
    } else ++d_it;

    // loop over neighbouring pairs and check, if QM-MM links were properly defined
    const std::vector<unsigned> indices = {j,k,l};
    for (std::vector<unsigned>::const_iterator
        j_it = indices.begin(), to = indices.end()
        ; j_it != to; ++j_it) {
      if ((   topo.is_qm(i)
          && !topo.is_qm(*j_it)
          && !are_linked(i, *j_it))
        ||
          (  !topo.is_qm(i)
          &&  topo.is_qm(*j_it)
          && !are_linked(*j_it, i)
        ))
        {
        std::ostringstream msg;
        msg << "Improper dihedral angle bending interaction across QM-MM boundary, "
            << "but no QMMM link defined between atoms "
            << i + 1 << " and " << *j_it + 1;
        io::messages.add(msg.str(), "QMMM_interaction", io::message::warning);
      }
    }
  }

  // Delete All-QM crossdihedrals (where are they used?)
  for (eightBodyVec::iterator d_it = crossdihedrals.begin(); d_it != crossdihedrals.end(); ) {
    if (   (topo.is_qm(d_it->a))
        && (topo.is_qm(d_it->b))
        && (topo.is_qm(d_it->c))
        && (topo.is_qm(d_it->d))
        && (topo.is_qm(d_it->e))
        && (topo.is_qm(d_it->f))
        && (topo.is_qm(d_it->g))
        && (topo.is_qm(d_it->h))
      )
      {
      os << "\t\tremoved crossdihedral bending term: "
         << d_it->a + 1 << " " << d_it->b + 1 << " "
         << d_it->c + 1 << " " << d_it->d + 1 << " "
         << d_it->e + 1 << " " << d_it->f + 1 << " "
         << d_it->g + 1 << " " << d_it->h + 1 << std::endl;
      d_it = crossdihedrals.erase(d_it);
    } else ++d_it;
  }
}

void interaction::QMMM_Interaction::remove_exclusions(
                                topology::Topology& topo
                              , const simulation::Simulation& sim
                              , std::ostream& os
                              , bool quiet)
    /*** MECHANICAL EMBEDDING */
    /** We need to remove:
     *    1. QM-QM one-four pairs - should be 0
     *    2. QM-QM exclusions - to avoid RF correction term
     *                        - they are excluded anyway within pairlist
     *                        - keep copy and calculate LJ, if QMLJ is on
     *    3. We need to do QM-QM LJ exceptions here (if QMLJ is ON)
     * * QM-MM pairs are done classically with full charges, exclusions and LJs
     */

    /** ELECTROSTATIC EMBEDDING */
    /** The same as mechanical embedding, but also remove and make a copy of:
     *    1. QM-MM one-four pairs - should be run here with LJ only
     *                            - this maybe makes sense only for
     *                              link-atom schemes
     *    2. QM-MM exclusions - to avoid RF and LS correction term
     *    3. QM-MM exceptions - should be done here
     */
  {
  DEBUG(4, "Removing exclusions of QM atoms");
  topo.qm_all_exclusion().clear();
  topo.qm_all_exclusion().resize(topo.num_atoms());
  // Remove QM-QM and QM-MM exclusions
  for (unsigned int i = 0; i < topo.num_solute_atoms(); ++i) {
    const bool i_is_qm = topo.is_qm(i);
    if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical
        && !i_is_qm) continue;
        
    // One-four pairs - they use CS6 and CS12
    for (topology::excl_cont_t::value_type::const_iterator
          it = topo.one_four_pair(i).begin()
        ; it != topo.one_four_pair(i).end(); ) {
      switch (i_is_qm + topo.is_qm(*it)) {
        case 0 : { // MM-MM
          ++it;
          break;
        }
        case 1 : { // QM-MM
          if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
            // Make a copy, if we want to calculate them
            DEBUG(9, "Making copy of 1,4: " << i << " - " << *it);
            topo.qm_one_four_pair(i).insert(*it);
            topo.qm_all_exclusion(i).insert(*it);
            // Erase
            DEBUG(9, "Removing 1,4 pair: " << i << " - " << *it);
            it = topo.one_four_pair(i).erase(it);
          } else ++it;
          break;
        }
        case 2 : { // QM-QM
          if (sim.param().qmmm.qm_lj) {
            DEBUG(9, "Making copy of 1,4: " << i << " - " << *it);
            topo.qm_one_four_pair(i).insert(*it);
            topo.qm_all_exclusion(i).insert(*it);
          }
          DEBUG(9, "Removing 1,4 pair: " << i << " - " << *it);
          it = topo.one_four_pair(i).erase(it);
        }
      }
    }

    // Exclusions
    for (topology::excl_cont_t::value_type::const_iterator
          it = topo.exclusion(i).begin()
        ; it != topo.exclusion(i).end(); ) {
      switch (i_is_qm + topo.is_qm(*it)) {
        case 0 : { // MM-MM
          ++it;
          break;
        }
        case 1 : { // QM-MM
          if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
            // Make a copy
            DEBUG(9, "Making copy of exclusion: " << i << " - " << *it);
            topo.qm_exclusion(i).insert(*it);
            topo.qm_all_exclusion(i).insert(*it);
            // Erase
            DEBUG(9, "Removing exclusion: " << i << " - " << *it);
            it = topo.exclusion(i).erase(it);
          } else ++it;
          break;
        }
        case 2 : { // QM-QM
          if (sim.param().qmmm.qm_lj) {
            // Make a copy
            DEBUG(9, "Making copy of exclusion: " << i << " - " << *it);
            topo.qm_exclusion(i).insert(*it);
            topo.qm_all_exclusion(i).insert(*it);
          }
          DEBUG(9, "Removing exclusion: " << i << " - " << *it);
          it = topo.exclusion(i).erase(it);
        }
      }
    }
  }

  // LJ Exceptions - they use LJEX
  for (std::vector<topology::lj_exception_struct>::const_iterator
        it = topo.lj_exceptions().begin()
      ; it != topo.lj_exceptions().end(); ) {
    const unsigned i = it->i;
    const unsigned j = it->j;
    assert(i < j);
    switch (topo.is_qm(i) + topo.is_qm(j)) {
      case 0 : { // MM-MM
        ++it;
        break;
      }
      case 1 : { // QM-MM
        if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
          // Make a copy
          DEBUG(9, "Making copy of LJ exception: " << i << " - " << j);
          topo.qm_lj_exceptions().push_back(*it);
          topo.qm_all_exclusion(i).insert(j);
          // Erase
          DEBUG(9, "Removing LJ exception: " << i << " - " << j);
          it = topo.lj_exceptions().erase(it);
        } else ++it;
        break;
      }
      case 2 : { // QM-QM
        if (sim.param().qmmm.qm_lj) {
          DEBUG(9, "Making copy of LJ exception: " << i << " - " << j);
          topo.qm_lj_exceptions().push_back(*it);
          topo.qm_all_exclusion(i).insert(j);
        }
        // Erase
        DEBUG(9, "Removing LJ exception: " << i << " - " << j);
        it = topo.lj_exceptions().erase(it);
        break;
      }
    }
  }
  topo.update_all_exclusion();
}

void interaction::QMMM_Interaction::store_set_data(
        const topology::Topology& topo,
        configuration::Configuration& conf,
        const simulation::Simulation& sim
        ) {
  m_timer.start("set data summation");
  std::vector<QMMM_Nonbonded_Set *>::iterator
    it = m_qmmm_nonbonded_set.begin(),
    to = m_qmmm_nonbonded_set.end();

  // add the forces, energies, virial...
  for (; it != to; ++it) {
    DEBUG(7, "adding forces from set " << it - m_qmmm_nonbonded_set.begin());
    (*it)->update_configuration(topo, conf, sim);
  }
  m_timer.stop("set data summation");
}

int interaction::QMMM_Interaction::print_pairlist(const topology::Topology& topo
                                                , std::ostream & os) {
  DEBUG(4, "QMMM_Interaction::print_pairlist");

  Pairlist temp_solute, temp_solvent;
  temp_solute.resize(topo.num_atoms());
  temp_solvent.resize(topo.num_atoms());

  for (unsigned int atom_i = 0; atom_i < topo.num_atoms(); ++atom_i) {

    for (unsigned i = 0; i < m_set_size; ++i) {

      assert(m_qmmm_nonbonded_set.size() > unsigned(i));
      assert(m_qmmm_nonbonded_set[i]->pairlist().solute_short.size() > atom_i);
      assert(m_qmmm_nonbonded_set[i]->pairlist().solvent_short.size() > atom_i);

      for (unsigned int atom_j = 0;
              atom_j < m_qmmm_nonbonded_set[i]->pairlist().solute_short[atom_i].size();
              ++atom_j) {

        assert(temp_solute.size() > atom_i);
        assert(temp_solute.size() > m_qmmm_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j]);

        if (m_qmmm_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j] < atom_i)
          temp_solute[m_qmmm_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j]].push_back(atom_i);
        else
          temp_solute[atom_i].push_back(m_qmmm_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j]);
      }
      for (unsigned int atom_j = 0;
              atom_j < m_qmmm_nonbonded_set[i]->pairlist().solvent_short[atom_i].size();
              ++atom_j) {

        assert(temp_solvent.size() > atom_i);
        assert(temp_solvent.size() > m_qmmm_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j]);

        if (m_qmmm_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j] < atom_i)
          temp_solvent[m_qmmm_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j]].push_back(atom_i);
        else
          temp_solvent[atom_i].push_back(m_qmmm_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j]);
      }
    }
  }
  os << temp_solute << std::endl
     << temp_solvent << std::endl;
  return 0;
}

void interaction::QMMM_Interaction::print_timing(std::ostream & os)
{
  m_timer.print(os);
  m_worker->timer().print(os);
}