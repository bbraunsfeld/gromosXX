/**
 * @file vgrid.h
 * nonbonded grid
 */

#ifndef INCLUDED_VGRID_H
#define INCLUDED_VGRID_H

namespace topology
{
  class Topology;
}
namespace configuration
{
  class Configuration;
}
namespace simulation
{
  class Simulation;
}

namespace interaction
{
  class Nonbonded_Parameter;
  
  void grid(topology::Topology const & topo,
	    configuration::Configuration & conf,
	    simulation::Simulation const & sim,
	    Nonbonded_Parameter & param);
  
} // interaction

#endif