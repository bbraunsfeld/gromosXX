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
static int _calculate_dfunct_interactions(topology::Topology& topo, 
																					configuration::Configuration& conf, 
																					simulation::Simulation& sim) {
	// shorten the code 
	int    atom_i   = sim.param().dfunct.atom_i; 
	int    atom_j   = sim.param().dfunct.atom_j;
	int    atom_k   = sim.param().dfunct.atom_k;
	int    atom_l   = sim.param().dfunct.atom_l;
	double target   = sim.param().dfunct.target;
	int    d        = sim.param().dfunct.d;
	double force    = sim.param().dfunct.force;
	DEBUG(10, "DFUNCT atom_i " << math::v2s(atom_i));
	DEBUG(10, "DFUNCT atom_k " << math::v2s(atom_j));
	DEBUG(10, "DFUNCT atom_k " << math::v2s(atom_k));
	DEBUG(10, "DFUNCT atom_l " << math::v2s(atom_l));
	// atomic distances expressed as vectors
	math::Vec dist_vec_ij, dist_vec_kl, dist_vec_ijkl;
	// find nearest periodic copies
	math::Periodicity<B> periodicity(conf.current().box);
	periodicity.nearest_image(conf.current().pos(atom_i), -conf.current().pos(atom_j), dist_vec_ij);
	periodicity.nearest_image(conf.current().pos(atom_k), -conf.current().pos(atom_l), dist_vec_kl);
	double dist_ij   = math::abs(dist_vec_ij);
	double dist_kl   = math::abs(dist_vec_kl);
	DEBUG(10, "DFUNCT dist_vec_ij " << math::v2s(dist_ij));
	DEBUG(10, "DFUNCT dist_vec_kl " << math::v2s(dist_kl));
	DEBUG(10, "DFUNCT dist_ij " << math::v2s(dist_ij));
	DEBUG(10, "DFUNCT dist_kl " << math::v2s(dist_kl));
	// scale distances
	dist_vec_ij = 0.5 * dist_vec_ij;
	dist_vec_kl = 0.5 * dist_vec_kl;
	periodicity.nearest_image(dist_vec_ij, dist_vec_kl, dist_vec_ijkl);
	double dist_ijkl = math::abs(dist_vec_ijkl);
	DEBUG(10, "DFUNCT dist_vec_ijkl " << math::v2s(dist_vec_ijkl));
	DEBUG(10, "DFUNCT dist_ijkl " << math::v2s(dist_ijkl));
	// compute forces on atoms 1...4 and combined biasing potential
	// the force is identical on all atoms - just the sign changes
	// force: -k * (unit_vec) * (r_ij + d * r_kl - R_0) 
	// https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.0c01112 
	math::Vec force_on_atoms = -force * dist_vec_ijkl / dist_ijkl * (dist_ij + d * dist_kl - target);
	math::Vec force_i =  force_on_atoms;
	math::Vec force_j =  force_on_atoms;
	math::Vec force_k = -force_on_atoms;
	math::Vec force_l = -force_on_atoms;
	double V_bias = 0.5 * force * (dist_ij + d * dist_kl - target) * (dist_ij + d * dist_kl - target);
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
  SPLIT_BOUNDARY(_calculate_dfunct_interactions, topo, conf, sim);
	return 0;
}