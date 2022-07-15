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
	int atom_1 = sim.param().dfunct.atom_1; 
	int atom_2 = sim.param().dfunct.atom_2;
	int atom_3 = sim.param().dfunct.atom_3;
	int atom_4 = sim.param().dfunct.atom_4;
	double target = sim.param().dfunct.target;
	int d = sim.param().dfunct.d;
	double force = sim.param().dfunct.force;
	math::Periodicity<B> periodicity(conf.current().box);
	math::Vec dist_vec_12, dist_vec_34;
	math::Vec dist_vec_3412;
	periodicity.nearest_image(conf.current().pos(atom_1), -conf.current().pos(atom_2), dist_vec_12);
	periodicity.nearest_image(conf.current().pos(atom_3), -conf.current().pos(atom_4), dist_vec_34);
	dist_vec_12 = 0.5 * dist_vec_12;
	dist_vec_34 = 0.5 * dist_vec_34;
	periodicity.nearest_image(dist_vec_12, dist_vec_34, dist_vec_3412);
	double dist_12 = math::abs(dist_vec_12);
	double dist_34 = math::abs(dist_vec_34);
	double dist_3412 = math::abs(dist_vec_3412);
	math::Vec force_1, force_2, force_3, force_4;
	force_1 = -force * dist_vec_3412 / dist_3412 * (dist_3412 - target);
	force_2 =  force_1;
	force_3 = -force_1;
	force_4 = -force_1;
	double V_bias = 0.5 * force * (dist_12 + d * dist_34) * (dist_12 + d * dist_34);
	conf.current().force(atom_1) += force_1;
	conf.current().force(atom_2) += force_2;
	conf.current().force(atom_3) += force_3;
	conf.current().force(atom_4) += force_4;
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