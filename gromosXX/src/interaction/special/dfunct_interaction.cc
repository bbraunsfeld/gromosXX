/**
 * @file dfunct_interaction.cc
 * dfunct
 * 
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"
#include "../../math/gmath.h"
//copied input output begin
#include <io/instream.h>
#include <io/blockinput.h>
#include <io/parameter/in_parameter.h>
#include <iostream>
#include <fstream>
//copied input output end
// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/dfunct_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

template<math::boundary_enum B>
static int _calculate_dfunct_substitution_form(topology::Topology& topo, 
																					     configuration::Configuration& conf, 
																					     simulation::Simulation& sim) {
	// shorten the code 
	int atom_i = sim.param().dfunct.atom_i; 
	int atom_j = sim.param().dfunct.atom_j;
	int atom_k = sim.param().dfunct.atom_k;
	int atom_l = sim.param().dfunct.atom_l;
	double r_0 = sim.param().dfunct.r_0;
	int d = sim.param().dfunct.d;
	double force = sim.param().dfunct.force;
	DEBUG(10, "DFUNCT atom_i " << atom_i << math::v2s(conf.current().pos(atom_i)));
	DEBUG(10, "DFUNCT atom_j " << atom_j << math::v2s(conf.current().pos(atom_j)));
	DEBUG(10, "DFUNCT atom_k " << atom_k << math::v2s(conf.current().pos(atom_k)));
	DEBUG(10, "DFUNCT atom_l " << atom_l << math::v2s(conf.current().pos(atom_l)));
	// atomic distances expressed as vectors
	// in GROMOS vector r_ji is defined as the vector from point i to point j (r_j - r_i)
	math::Vec dist_vec_ji, dist_vec_lk;
	// find nearest periodic copies
	math::Periodicity<B> periodicity(conf.current().box);
	periodicity.nearest_image(conf.current().pos(atom_j), conf.current().pos(atom_i), dist_vec_ji);
	periodicity.nearest_image(conf.current().pos(atom_l), conf.current().pos(atom_k), dist_vec_lk);
	double dist_ji   = math::abs(dist_vec_ji);
	double dist_lk   = math::abs(dist_vec_lk);
	DEBUG(10, "DFUNCT dist_vec_ji " << math::v2s(dist_vec_ji));
	DEBUG(10, "DFUNCT dist_vec_lk " << math::v2s(dist_vec_lk));
	DEBUG(10, "DFUNCT dist_ji " << dist_ji);
	DEBUG(10, "DFUNCT dist_lk " << dist_lk);
	// compute forces on atoms i, j, k, l and the combined biasing potential
	math::Vec prefactor_i =     force * (dist_vec_ji / dist_ji);
	math::Vec prefactor_k = d * force * (dist_vec_lk / dist_lk);
	double force_term = (dist_ji + d * dist_lk - r_0);
	math::Vec force_i =  prefactor_i  * force_term;
	math::Vec force_j = -force_i;
	math::Vec force_k =  prefactor_k * force_term;
	math::Vec force_l = -force_k;
	double V_bias = 0.5 * force * (dist_ji + d * dist_lk - r_0) * (dist_ji + d * dist_lk - r_0);
	DEBUG(10, "DFUNCT Prefactor i " << math::v2s(prefactor_i));
	DEBUG(10, "DFUNCT Prefactor k " << math::v2s(prefactor_k));
	DEBUG(10, "DFUNCT Force on i " << math::v2s(force_i));
  DEBUG(10, "DFUNCT Force on j " << math::v2s(force_j));
  DEBUG(10, "DFUNCT Force on k " << math::v2s(force_k));
  DEBUG(10, "DFUNCT Force on l " << math::v2s(force_l));
	DEBUG(10, "DFUNCT V_bias " << V_bias);
	// store forces and biasing potential
	conf.current().force(atom_i) += force_i;
	conf.current().force(atom_j) += force_j;
	conf.current().force(atom_k) += force_k;
	conf.current().force(atom_l) += force_l;
	conf.current().energies.distanceres_total += V_bias;
	return 0;
}

template<math::boundary_enum B>
static int _calculate_dfunct_diels_alder_form(topology::Topology& topo, 
																					     configuration::Configuration& conf, 
																					     simulation::Simulation& sim) {
	/**
	 * @brief to be implemented
	 * 
	 */
	// shorten the code 
	int    atom_i   = sim.param().dfunct.atom_i; 
	int    atom_j   = sim.param().dfunct.atom_j;
	int    atom_k   = sim.param().dfunct.atom_k;
	int    atom_l   = sim.param().dfunct.atom_l;
	double r_0   = sim.param().dfunct.r_0;
	int    d        = sim.param().dfunct.d;
	double force    = sim.param().dfunct.force;
	DEBUG(10, "DFUNCT atom_i " << atom_i);
	DEBUG(10, "DFUNCT atom_k " << atom_j);
	DEBUG(10, "DFUNCT atom_k " << atom_k);
	DEBUG(10, "DFUNCT atom_l " << atom_l);
	// atomic distances expressed as vectors
	math::Vec dist_vec_ij, dist_vec_kl, dist_vec_ijkl;
	// find nearest periodic copies
	math::Periodicity<B> periodicity(conf.current().box);
	periodicity.nearest_image(conf.current().pos(atom_i), -conf.current().pos(atom_j), dist_vec_ij);
	periodicity.nearest_image(conf.current().pos(atom_k), -conf.current().pos(atom_l), dist_vec_kl);
	double dist_ij   = math::abs(dist_vec_ij);
	double dist_kl   = math::abs(dist_vec_kl);
	DEBUG(10, "DFUNCT dist_vec_ij " << math::v2s(dist_vec_ij));
	DEBUG(10, "DFUNCT dist_vec_kl " << math::v2s(dist_vec_kl));
	DEBUG(10, "DFUNCT dist_ij " << dist_ij);
	DEBUG(10, "DFUNCT dist_kl " << dist_kl);
	// scale distances
	math::Vec dist_vec_ij_halfs = 0.5 * dist_vec_ij;
	math::Vec dist_vec_kl_halfs = 0.5 * dist_vec_kl;
	periodicity.nearest_image(dist_vec_ij, dist_vec_kl, dist_vec_ijkl);
	double dist_ijkl = math::abs(dist_vec_ijkl);
	DEBUG(10, "DFUNCT dist_vec_ijkl " << math::v2s(dist_vec_ijkl));
	DEBUG(10, "DFUNCT dist_ijkl " << dist_ijkl);
	// compute forces on atoms i, j, k, l and the combined biasing potential
	// the force is identical on all atoms - only the sign changes
	// force with respect to r_ij: k * (unit_vec) * (r_ij + d * r_kl - r_0)
	// force with respect to r_kl: d * k * (unit_vec) * (r_ij + d * r_kl - r_0) 
	// https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.0c01112 
	// math::Vec force_i = force * (dist_vec_ijkl / dist_ijkl) * (dist_ij + d * dist_kl - r_0);
	// math::Vec force_j = force_i;
	// math::Vec force_k = d * force * (dist_vec_ijkl / dist_ijkl) * (dist_ij + d * dist_kl - r_0);
	// math::Vec force_l = force_k;
	double V_bias = 0.5 * force * (dist_ij + d * dist_kl - r_0) * (dist_ij + d * dist_kl - r_0);
	math::Vec force_i = force * dist_vec_ijkl / dist_ijkl * (dist_ij + d * dist_kl - r_0);
	math::Vec force_j = force_i;
	math::Vec force_k = d * force_i;
	math::Vec force_l = force_k;
	DEBUG(10, "DFUNCT Force on i " << math::v2s(force_i));
  DEBUG(10, "DFUNCT Force on j " << math::v2s(force_j));
  DEBUG(10, "DFUNCT Force on k " << math::v2s(force_k));
  DEBUG(10, "DFUNCT Force on l " << math::v2s(force_l));
	// store forces and biasing potential
	conf.current().force(atom_i) += force_i;
	conf.current().force(atom_j) += force_j;
	conf.current().force(atom_k) += force_k;
	conf.current().force(atom_l) += force_l;
	conf.current().energies.distanceres_total += V_bias;
	return 0;
}

int interaction::DFunct_Interaction::init(topology::Topology& topo,
		     																  configuration::Configuration& conf,
		     																  simulation::Simulation& sim,
		     																  std::ostream& os,
		     																  bool quiet) {
  
	return 0;
}

int interaction::DFunct_Interaction::calculate_interactions(topology::Topology & topo,
				                                                    configuration::Configuration& conf,
				                                                    simulation::Simulation& sim) {
  switch (sim.param().dfunct.dfunct) {
		case simulation::dfunct_substitution:
			SPLIT_BOUNDARY(_calculate_dfunct_substitution_form, topo, conf, sim);
			break;

		case simulation::dfunct_diels_alder:
			SPLIT_BOUNDARY(_calculate_dfunct_diels_alder_form, topo, conf, sim);
			break;

		default:
      io::messages.add("DFunct functional not implemented", "DFunct_Interaction", io::message::critical);
      break;
		}
	
	return 0;
}