/**
 * @file perturbed_soft_bond_interaction.h
 * perturbed soft bond interaction.
 */

#ifndef INCLUDED_PERTURBED_SOFT_BOND_INTERACTION
#define INCLUDED_PERTURBED_SOFT_BOND_INTERACTION

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
   * @class Perturbed_Soft_Bond_Interaction
   * calculates the perturbed bond interactions (soft harmonic).
   * Wang et al (2017) J. Chem. Theory Comput., 13 (1), 42–54
   */
  class Perturbed_Soft_Bond_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Soft_Bond_Interaction()
      : Interaction("PerturbedSoftBond") {}

    
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Soft_Bond_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) ;

    /**
     * calculate the interactions (force and energy, lambda derivative)
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  };
  
} // interaction

#endif
