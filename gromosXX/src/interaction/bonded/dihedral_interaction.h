/**
 * @file dihedral_interaction.h
 * dihedral interaction.
 */

#ifndef INCLUDED_DIHEDRAL_INTERACTION_H
#define INCLUDED_DIHEDRAL_INTERACTION_H

namespace configuration{
	class Configuration;
}
namespace topology{
	class Topology;
}
namespace simulation{
	class Simulation;
}

namespace interaction
{
  /**
   * @class Dihedral_Interaction
   * calculates the dihedral interactions.
   */
  class Dihedral_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Dihedral_Interaction() : Interaction("Dihedral") {}
    /**
     * Destructor.
     */
    virtual ~Dihedral_Interaction() {}

    /**
     * init
     */
 
   virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
                     bool quiet = false);
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  protected:

    /**
     * calculate nearest minimum
     */
    // double _calculate_nearest_minimum(double phi, int m, double pd);
    
  };
  
} // interaction

#endif
