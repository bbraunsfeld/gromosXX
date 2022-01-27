/* 
 * File:   replica_exchange_master_2d_l_T_HREMD.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:06 PM
 */

#include <replicaExchange/replica_exchangers/replica_exchange_slave_interface.h>
#include <replicaExchange/replica_exchangers/2D_T_lambda_REPEX/replica_exchange_slave_2d_l_T_HREMD.h>

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <io/argument.h>
#include <time.h>
#include <unistd.h>

#include <io/configuration/out_configuration.h>


#ifdef XXMPI
#include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE re
#define SUBMODULE replica_exchanger


re::replica_exchange_slave_2d_l_T_HREMD::replica_exchange_slave_2d_l_T_HREMD(io::Argument & _args,
                                                                             unsigned int cont,
                                                                             unsigned int globalThreadID,
                                                                             replica_graph_control & replicaGraphMPIControl,
                                                                             simulation::MpiControl & replica_mpi_control) :
        replica_exchange_base_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
        replica_exchange_base_2d_l_T_HREMD(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
        replica_exchange_slave_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control) {
        DEBUG(4, "replica_exchange_slave_2d_l_T_HREMD "<< globalThreadID <<":Constructor:\t START");

        DEBUG(4, "replica_exchange_slave_2d_l_T_HREMD "<< globalThreadID <<":Constructor:\t DONE");
}

re::replica_exchange_slave_2d_l_T_HREMD::~replica_exchange_slave_2d_l_T_HREMD() {

}