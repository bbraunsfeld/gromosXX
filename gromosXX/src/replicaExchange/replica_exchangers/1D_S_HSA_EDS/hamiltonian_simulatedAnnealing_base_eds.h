/* 
 * File:   hamiltonian_simulatedAnnealing_base_eds.h
 * Author: bschroed
 *
 * Created on April 18, 2018, 3:38 PM
 * Modified June 18, 2021 - bschroed, srieder
 */

#ifndef hamiltonian_simulatedAnnealing_base_eds_H
#define hamiltonian_simulatedAnnealing_base_eds_H

#include <replicaExchange/replica_exchangers/replica_exchange_base_interface.h>
#include <replicaExchange/replica_data.h>
#include <replicaExchange/repex_mpi.h>

//const

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <io/argument.h>
#include <util/usage.h>
#include <util/error.h>

#include <io/read_input.h>
#include <io/print_block.h>

#include <time.h>
#include <unistd.h>

#include <io/configuration/out_configuration.h>


#include <string>
#include <math/random.h>

#ifdef XXMPI
#include <mpi.h>
#endif

namespace re
{
     class hamiltonian_simulatedAnnealing_base_eds : public virtual replica_exchange_base_interface {
    public:
        hamiltonian_simulatedAnnealing_base_eds(io::Argument _args, 
                unsigned int cont, 
                unsigned int globalThreadID, 
                replica_graph_control & replicaGraphMPIControl,
                simulation::MpiControl & replica_mpi_control);
        /**
         * inits replica_reeds
         */
        void init() override;
        /** 
         * Initialize file and data for eds_stat output.
         */
        void init_eds_stat();
        
        
    protected:       
        virtual ~hamiltonian_simulatedAnnealing_base_eds();

       /*
       * Attributes
       */

        /**
         * for exchanging params easily
         */
        simulation::Parameter::reeds_struct& reedsParam;

       /**
       * stat information of all replicas for exchange optimisation 
       */
       std::map<ID_t, re::reeds_replica_stat_data > replicaStatData;
       /**
       *  output file stream for output file
       */
       std::map<ID_t, std::ofstream *> eds_stat_out;
       /**
        *  Exchange parameters
        */
       simulation::Parameter::eds_struct eds_para;
        /**
        * contains original forces which have to be reset after RE-EDS exchange energy calculation
        */
       math::VArray force_orig;
       /**
        * contains original virial which have to be reset after RE-EDS exchange energy calculation
        */
       math::Matrix virial_tensor_orig;
    
       
       /*
        * Functions
        */
       /**
        * Set RE - params
        */
       void setParams() override;

       void set_s();
       
        // RE-Exchange functions
        /**
        * Sets eds_struct() parameters to original value of replica
        */
        void reset_eds();
        /**
         * Sets  eds_struct() parameters to value of partner (for energy calculations)
         */
        void change_eds(const unsigned int partner);
        
        //Build exchange probabilities
        void determine_switch_probabilities() ;
        double calc_probability(const unsigned int partnerReplicaID) override;
        void swap_replicas_1D(const unsigned int partnerReplicaID);

         //EXECUTE SWAP:
         virtual void execute_swap(const unsigned int partnerReplicaID);

         /*
        * energy calculation for statistical purposes of eds_stat() in replica_exchange_base.cc
        * for given configuration with new smoothing parameter s.
        */
        double calc_energy_eds_stat(double s);
        double calculate_energy_core();
        double calculate_energy(const unsigned int partner);
    };
}
#endif /* hamiltonian_simulatedAnnealing_base_eds_H */
