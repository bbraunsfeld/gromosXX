/**
 * @file turbomole_worker.h
 * The worker class for the Turbomole QM software
 */
#ifndef INCLUDED_TURBOMOLE_WORKER_H
#define	INCLUDED_TURBOMOLE_WORKER_H

namespace interaction {
  class QM_Worker;
  /**
   * @class Turbomole_Worker
   * a worker class which calls the Turbomole software
   */
  class Turbomole_Worker : public QM_Worker {
  public:
    /**
     * Constructor
     */
    Turbomole_Worker() : QM_Worker("Turbomole Worker") {
    //this->get_new_qmID();

    }
    /**
     * Destructor
     */
    virtual ~Turbomole_Worker();
    /**
     * initialise the QM worker
     * @return 0 if successful, non-zero on failure
     */
    virtual int init(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);
    /**
     * run a QM job in Turbomole
     * @param qm_pos a vector containing the QM atom positions
     * @param mm_atoms the MM atoms to include
     * @param storage the energies, forces, charges obtained
     * @return 0 if successful, non-zero if not.
     */
    virtual int run_QM(const topology::Topology & topo,
                       const configuration::Configuration & conf,
                       const simulation::Simulation & sim,
                       interaction::QM_Zone & qm_zone);
    private:
    /**
     * file name for TM working directory
     */
    std::string working_directory;
    std::vector<std::string> tool_chain;
  };
}

#endif	/* TURBOMOLE_WORKER_H */
