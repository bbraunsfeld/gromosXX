/**
 * @file nonbonded_innerloop.cc
 * template methods of Nonbonded_Innerloop
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>
::lj_crf_innerloop
(topology::Topology & topo, configuration::Configuration & conf,
 unsigned int i, unsigned int j,
 Storage & storage,
 Periodicity_type const & periodicity,
 int pc)
{
  // storage.energies.bond_energy[0] += j;

  DEBUG(8, "\tpair\t" << i << "\t" << j);
  
  math::Vec r;
  double f;
  double e_lj, e_crf;
  
  if (t_nonbonded_spec::do_bekker){
    r = conf.current().pos(i) + periodicity.shift(pc).pos
      - conf.current().pos(j);
    DEBUG(10, "\tpc=" << pc << " shift = " << periodicity.shift(pc).pos(0)
	  << " / " << periodicity.shift(pc).pos(1) 
	  << " / " << periodicity.shift(pc).pos(2));
    DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
  }
  else{
    periodicity.nearest_image(conf.current().pos(i), 
			      conf.current().pos(j), r);
    DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
    
  }
    
  const lj_parameter_struct &lj = 
    m_param->lj_parameter(topo.iac(i),
			  topo.iac(j));
  
  DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
  DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
  
  lj_crf_interaction(r, lj.c6, lj.c12,
		     topo.charge()(i) * 
		     topo.charge()(j),
		     f, e_lj, e_crf);
  
  // most common case
  if (t_nonbonded_spec::do_virial == math::molecular_virial){
	  math::Vec rf = f * r;
	  storage.force(i) += rf;
	  storage.force(j) -= rf;

	  for(int b=0; b<3; ++b){
		  const double rr = r(b) - conf.special().rel_mol_com_pos(i)(b) + conf.special().rel_mol_com_pos(j)(b);
		  for(int a=0; a<3; ++a){
			  storage.virial_tensor(b, a) += rr * rf(a);
		  }
	  }
  }
  else{
	for (int a=0; a<3; ++a){
    
		const double term = f * r(a);
		storage.force(i)(a) += term;
		storage.force(j)(a) -= term;

	    if (t_nonbonded_spec::do_virial == math::atomic_virial){
			for(int b=0; b<3; ++b){
				storage.virial_tensor(b, a) += 
				r(b) * term;
			}
		}
	}
  }

  // energy
  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
	<< " j " << topo.atom_energy_group(j));

  assert(storage.energies.lj_energy.size() > 
	 topo.atom_energy_group(i));
  assert(storage.energies.lj_energy.size() >
	 topo.atom_energy_group(j));
  
  storage.energies.lj_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_lj;
  
  storage.energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_crf;
  
}


template<typename t_nonbonded_spec>
void interaction::Nonbonded_Innerloop<t_nonbonded_spec>
::one_four_interaction_innerloop
(topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i, unsigned int j,
 Periodicity_type const & periodicity)
{
  DEBUG(8, "\t1,4-pair\t" << i << "\t" << j);
  
  math::Vec r;
  double f, e_lj, e_crf;
  
  periodicity.nearest_image(conf.current().pos(i), 
			    conf.current().pos(j), r);

  const lj_parameter_struct &lj = 
    m_param->lj_parameter(topo.iac(i),
			  topo.iac(j));
  
  DEBUG(11, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);
  
  lj_crf_interaction(r, lj.cs6, lj.cs12,
		     topo.charge()(i) * 
		     topo.charge()(j),
		     f, e_lj, e_crf);
  
  for (int a=0; a<3; ++a){

    const double term = f * r(a);
    // storage.force(i)(a) += term;
    // storage.force(j)(a) -= term;
    conf.current().force(i)(a) += term;
    conf.current().force(j)(a) -= term;
    
    if (t_nonbonded_spec::do_virial == math::atomic_virial){
      for(int b=0; b<3; ++b){
	conf.current().virial_tensor(b, a) += 
	  r(b) * term;
      }
    }
  }
  
  // energy
  conf.current().energies.lj_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_lj;
    
  conf.current().energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_crf;
    
  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
	<< " j " << topo.atom_energy_group(j));
    
}

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>
::RF_excluded_interaction_innerloop
(topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i,
 Periodicity_type const & periodicity)
{
  math::Vec r, f;
  double e_crf;
  
  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;

  std::set<int>::const_iterator it, to;
  it = topo.exclusion(i).begin();
  to = topo.exclusion(i).end();
  
  DEBUG(8, "\tself-term " << i );
  r=0;
  
  // this will only contribute in the energy, the force should be zero.
  rf_interaction(r,topo.charge()(i) * topo.charge()(i),
		 f, e_crf);
  conf.current().energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(i)] += 0.5 * e_crf;
  DEBUG(11, "\tcontribution " << 0.5*e_crf);
  
  for( ; it != to; ++it){
    
    DEBUG(11, "\texcluded pair " << i << " - " << *it);
    
    periodicity.nearest_image(pos(i), pos(*it), r);
    
    
    rf_interaction(r, topo.charge()(i) * 
		   topo.charge()(*it),
		   f, e_crf);
    
    force(i) += f;
    force(*it) -= f;
    
    if (t_nonbonded_spec::do_virial == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int b=0; b<3; ++b)
	  conf.current().virial_tensor(a, b) += 
	    r(a) * f(b);
      DEBUG(11, "\tatomic virial done");
    }

    // energy
    conf.current().energies.crf_energy[topo.atom_energy_group(i)]
      [topo.atom_energy_group(*it)] += e_crf;
    
    DEBUG(11, "\tcontribution " << e_crf);
    
  } // loop over excluded pairs
  
}

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>
::RF_solvent_interaction_innerloop
(topology::Topology & topo,
 configuration::Configuration & conf,
 topology::Chargegroup_Iterator const & cg_it,
 Periodicity_type const & periodicity)
{
  math::Vec r;
  double e_crf;
  
  math::VArray &pos   = conf.current().pos;

  // loop over the atoms
  topology::Atom_Iterator at_it = cg_it.begin(),
    at_to = cg_it.end();
  
  for ( ; at_it != at_to; ++at_it){
    DEBUG(11, "\tsolvent self term " << *at_it);
    // no solvent self term. The distance dependent part and the forces
    // are zero. The distance independent part should add up to zero 
    // for the energies and is left out.
    
    for(topology::Atom_Iterator at2_it=at_it+1; at2_it!=at_to; ++at2_it){
      
      DEBUG(11, "\tsolvent " << *at_it << " - " << *at2_it);
      periodicity.nearest_image(pos(*at_it), 
				pos(*at2_it), r);
      
      // for solvent, we don't calculate internal forces (rigid molecules)
      // and the distance independent parts should go to zero
      e_crf = -topo.charge()(*at_it) * 
	topo.charge()(*at2_it) * 
	math::four_pi_eps_i * 
	crf_2cut3i() * abs2(r);
      
      // energy
      conf.current().energies.crf_energy
	[topo.atom_energy_group(*at_it) ]
	[topo.atom_energy_group(*at2_it)] += e_crf;
    } // loop over at2_it
  } // loop over at_it
  
}