/**
 * @file cudaKernel.h
 * contains struct definitions
 * contains function declarations used in GROMOSXX MD++
 */

#define CUKERNEL_TYPES_ONLY
#include "lib/types.h"
#include "parameter.h"
#include "pairlist.h"
#define gpu_status void
#undef CUKERNEL_TYPES_ONLY

namespace cudakernel {
  extern "C" {
    /**
     * initializes the MD-paramters allocates and initializes the memory on device
     * @param[inout] device_number the CUDA device number. If -1 is given, the driver will determine the value and return it in this parameter.
     * @param[in] num_atoms the number of atoms
     * @param[in] cutoff_short the short-range cutoff
     * @param[in] cutoff_long the long-range cutoff
     * @param[in] box the edge of the box
     * @param[in] num_atoms_per_mol number of atoms per solvent molecule
     * @param[in] estimated_neighbors_short the estimated number of atoms in the short-range pairlist
     * @param[in] estimated_neighbors_long the estimated number of atoms in the long-range pairlist
     * @param[in] crf_2cut3i a reaction field constant
     * @param[in] crf_cut a rection field constant
     * @param[in] crf_cut3i a reaction field constant
     * @param[in] lj_crf_params the nonbonded interaction parameters.
     */
    gpu_status * cudaInit(int & device_number,
            unsigned int num_atoms,
            double cutoff_short,
            double cutoff_long,
            double box_x,
            double box_y,
            double box_z, 
            unsigned int num_atoms_per_mol,
            unsigned int estimated_neighbors_short,
            unsigned int estimated_neighbors_long,
            double crf_2cut3i, double crf_cut, double crf_cut3i,
            lj_crf_parameter * lj_crf_params,
            unsigned int num_of_gpus,
            unsigned int gpu_id,
            int * error);

    /**
     * initialize the constraints solver
     * @param[in] constr the cosntraints
     * @param[in] mass the masses
     * @param[in] tol the tolerance for itertative algorithms like SHAKE
     */
    //void cudaInitConstraints(constraint * constr, double* mass, double tol, gpu_status * gpu_stat);
    /**
     * copies box parameters to GPU
     * @param[in] box the box edge length of the box.
     */
    int cudaCopyBox(gpu_status * gpu_stat, double box_x, double box_y, double box_z);

    /**
     * copies the positions CPU->GPU
     * @param pos the positions
     */
    int cudaCopyPositions(double * pos, gpu_status * gpu_stat);
    /**
     * calculates/updates the long/shortrange pairlists
     */
    void cudaCalcPairlist(gpu_status * gpu_stat);
    /**
     * calculates the long- or shortrange forces, energies and virial
     * adds the contributions of solvent-solvent interactions
     * @param[inout] forces the forces are added to this
     * @param[inout] virial the 3x3 virial tensor is added to this
     * @param[inout] lj_energy the Lennard-Jones energy is added to this
     * @param[inout] crf_energy the electrostatic energy is added to this
     * @param[in] longrange true: longrange interaction calculated, false: shortrange interactions calculated
     * @return 0 if successful or sum error codes otherwise.
     */
    int cudaCalcForces(double * forces, double * virial, double * lj_energy, double * crf_energy,
            bool longrange, gpu_status * gpu_stat);

    /**
     * Initializes the GPU for the constraints calculation
     */
    gpu_status * cudaInitConstraints(unsigned int num_of_gpus, unsigned int gpu_id, unsigned int num_atoms,
                                  unsigned int num_solvent_mol);
    /**
     * applyes to constraints to the positions. First call cudaCopyPositions if you haven't done so yet.
     * @param[inout] newpos the new positions to be constrained.
     * @param[in] oldpos the old positions
     * @param[out] shake_fail_mol -1 on success or a molecule index if the solving of constraints failed for this molecule.
     */
    int cudaConstraints(double * newpos, double * oldpos,
                    int & shake_fail_mol, gpu_status * gpu_stat);

    /**
     *  Initializes the GPU for the constraints calculation
     * @param[inout] device the CUDA device to use. If -1, the driver will determine and return the result
     * @param[in] constr pointer to the constraint lengths squared
     * @param[in] factor pointer to the factor matrix
     * @param[in] mass pointer to an array with the masses of the atoms in a molecule
     * @param[in] tol the tolerance
     * @param[in] num_gpus the number of GPUs
     * @param[in] gpu_id the device number of the gpu
     * @param[in] num_atoms The number of Atoms
     * @param[in] num_solvent_mol The number of solvent molecules
     */
    gpu_status * cudaInitGPU_Shake(int & device, double * constr,
                                  double * factor, double * mass, 
                                  double tol, unsigned int num_gpus, unsigned int gpu_id, unsigned int num_atoms,
                                  unsigned int num_solvent_mol, int * error);
    /**
     * applies constraints to the positions with M_SHAKE
     * @param[inout] newpos the newpositions, will be modified
     * @param[in] oldpos the old positions
     * @param[inout] shake_fail_mol at which molecule the shake algorithm fails
     * @param[in] gpu_stat all important pointers fo the memory handling on the GPU
     */
    int cudaGPU_Shake(double * newpos, double * oldpos,
                    int & shake_fail_mol, gpu_status * gpu_stat);
    /**
     * frees the memory allocated 
     */
    int CleanUp(gpu_status * gpu_stat);
    /**
     * a dummy function for linker checking
     */
    void test();
    
  }
}

