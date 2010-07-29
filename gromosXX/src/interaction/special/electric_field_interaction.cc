/**
 * @file electric_field_interaction.cc
 * template external Electric_Field_Interaction
 */


#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>

#include <math/periodicity.h>

// special interactions
#include <interaction/interaction_types.h>

#include <interaction/special/electric_field_interaction.h>

#include <util/template_split.h>
#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate electric field interactions
 */

int interaction::Electric_Field_Interaction::
calculate_interactions(topology::Topology& topo,
                       configuration::Configuration& conf,
                       simulation::Simulation& sim)
{
  // loop over the atoms
  // unit E = e . nm^-2
  math::Vec E(sim.param().electric.Ef_x, sim.param().electric.Ef_y, sim.param().electric.Ef_z);
  double scale = math::four_pi_eps_i;
  

  for (unsigned int i = 0; i < topo.num_atoms(); ++i){
    double q = scale * topo.charge(i);
    conf.current().force(i) += q*E;
  }

  //The electric energies were not added to the total energy
  /* F = qE
   * U = - integral(Fdr) => U = -qEr
   * The energy depends on the (absolute) position of the particle
   * Maybe it doesn't make sense to use this energy within
   * periodic boundary conditions
   */

  

  return 0;
}
