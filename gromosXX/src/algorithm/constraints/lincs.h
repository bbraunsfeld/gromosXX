/**
 * @file lincs.h
 * the LINCS algorithm.
 */

#ifndef INCLUDED_LINCS_H
#define INCLUDED_LINCS_H

namespace algorithm {

  /**
   * @class Lincs
   * implements the lincs algorithm.
   */
  class Lincs : public Algorithm {
  public:
    /**
     * Constructor.
     */
    Lincs();

    /**
     * Destructor.
     */
    virtual ~Lincs();

    virtual int apply(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);


    /**
     * accessor to the constrained atoms
     */
    std::set<unsigned int> & constrained_atoms() {
      return m_constrained_atoms;
    }

    /**
     * accessor to the constrained atoms
     */
    const std::set<unsigned int> & constrained_atoms() const {
      return m_constrained_atoms;
    }

    /**
     * initialize startup positions and velocities
     * if required.
     */
    int init(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            std::ostream & os = std::cout,
            bool quiet = false);

  protected:

    /**
     * the atoms that are involved in the contraints
     */
    std::set<unsigned int> m_constrained_atoms;

  };

} //algorithm

#endif
