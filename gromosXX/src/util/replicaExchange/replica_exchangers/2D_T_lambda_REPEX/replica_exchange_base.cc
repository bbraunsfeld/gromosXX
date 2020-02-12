/* 
 * File:   replica_exchange_base.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:06 PM
 */

#include "replicaExchange/replica_mpi_tools.h"
#include "replicaExchange/replica/replica.h"
#include "replicaExchange/replica/replica_MPI_master.h"
#include "replicaExchange/replica/replica_MPI_slave.h"
#include "replicaExchange/replica_exchangers/2D_T_lambda_REPEX/replica_exchange_base.h"

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
#include <math/random.h>
#include <math/volume.h>
#include <string>

#ifdef XXMPI
#include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

util::replica_exchange_base::replica_exchange_base(io::Argument _args,
                                                   unsigned int cont, 
                                                   unsigned int globalThreadID,
                                                   replica_graph_mpi_control replicaGraphMPIControl,
                                                   simulation::mpi_control_struct replica_mpi_control) : 
        replica_exchange_base_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control){
#ifdef XXMPI
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":Constructor:\t START ");
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":Constructor:\t SIMULATIONID:  "<< simulationID);

  //construct replica obj
  DEBUG(4,"replica_exchange_base "<< globalThreadID <<":Constructor:\t  createReplica");
  createReplicas(cont, globalThreadID, replica_mpi_control);
  DEBUG(4,"replica_exchange_base "<< globalThreadID <<":Constructor:\t createdReplica T");

  //RE-Vars
  // set some variables
  stepsPerRun = replica->sim.param().step.number_of_steps;
  run = 0;
  total_runs = replica->sim.param().replica.trials + replica->sim.param().replica.equilibrate;
  partnerReplicaID = simulationID;
  time = replica->sim.time();
  steps = 0;
  switched = 0;
  
  replica->curentStepNumber=0;
  replica->totalStepNumber = total_runs*stepsPerRun;
  replica->stepsPerRun= stepsPerRun;

  const int numT = replica->sim.param().replica.num_T;

  T = replica->sim.param().replica.temperature[simulationID % numT];
  l = replica->sim.param().replica.lambda[simulationID / numT];
  dt = replica->sim.param().replica.dt[simulationID / numT];

  set_lambda();
  set_temp();
  
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":Constructor:\t Constructor \t DONE");
  #else
    throw "Cannot use send_to_master from replica_exchange_slave_eds without MPI!"; 
  #endif
}


util::replica_exchange_base::~replica_exchange_base() {
    delete replica;
}