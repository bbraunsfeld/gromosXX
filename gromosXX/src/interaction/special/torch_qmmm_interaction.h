/**
 * @file torch_qmmm_interaction.h
 * torch_qmmm_interaction
 *
 */

#ifndef INCLUDED_TORCH_QMMM_INTERACTION_H
#define INCLUDED_TORCH_QMMM_INTERACTION_H

#include "../../stdheader.h"
#include "../../interaction/qmmm/qm_zone.h"
#include "../../interaction/qmmm/qmmm_interaction.h"
#include "../../interaction/special/torch_interaction.h"
#include "../../simulation/parameter.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include <torch/csrc/autograd/function.h>
#include <torch/csrc/autograd/variable.h>
#include <torch/torch.h>

namespace interaction {

/**
 * A specialized version of Torch_Interaction that sends the QM/MM zone to a
 * Torch model, typename T indicates numerical precision (half, single, double)
 */
template <typename T>
class Torch_QMMM_Interaction : public Torch_Interaction<T> {

public:
  /**
   * Initializes the interaction
   */
  Torch_QMMM_Interaction(const simulation::torch_model &model)
      : Torch_Interaction<T>(model, "Torch QM/MM Interface") {}

  /**
   * Deallocates resources
   */
  virtual ~Torch_QMMM_Interaction() = default;

  /**
   * Initializes everything necessary
   */
  virtual int init(topology::Topology &topo, configuration::Configuration &conf,
                   simulation::Simulation &sim, std::ostream &os = std::cout,
                   bool quiet = false) override;

private:
  /**
   * Sets-up a pointer to the QM zone of the simulation
   */
  virtual int init_qm_zone();

  /**
   * Prepares the coordinates according to the atom selection scheme selected
   */
  virtual int prepare_input(const simulation::Simulation &sim) override;

  /**
   * Sets-up the vector with QM atom types
   */
  virtual int init_qm_atom_numbers();

  /**
   * Sets-up the vector with QM coordinates - will also cast floating point precisions
   */
  virtual int prepare_qm_atoms();

  /**
   * Sets-up the vectors with MM atom types, charges, and coordinates - will
   * also cast floating point precision
   */
  virtual int prepare_mm_atoms();

  /**
   * Initializes tensors ready to go into the model
   */
  virtual int build_tensors(const simulation::Simulation &sim) override;

  /**
   * Forward pass of the model loaded
   */
  virtual int forward() override;

  /**
   * Backward pass of the model loaded
   */
  virtual int backward() override;

  /**
   * Gets energy from Torch and updates Gromos
   */
  virtual int update_energy(topology::Topology &topo,
                            configuration::Configuration &conf,
                            const simulation::Simulation &sim) override;

  /**
   * Gets forces from Torch and updates Gromos
   */
  virtual int update_forces(topology::Topology &topo,
                            configuration::Configuration &conf,
                            const simulation::Simulation &sim) override;

  /**
   * Saves Torch input data
   */
  virtual void save_torch_input(const unsigned int step
                              , const simulation::Simulation& sim) override;

  /**
   * Saves Torch output data
   */
  virtual void save_torch_output(const unsigned int step
                               , const simulation::Simulation& sim) override;
  
  /**
   * Saves the QM coordinates sent to Torch
   */
  virtual void save_input_coord(std::ofstream& ifs
                              , const unsigned int step);

  /**
   * Saves the MM coordinates and charges sent to Torch
   */
  virtual void save_input_point_charges(std::ofstream& ifs
                                      , const unsigned int step
                                      , const unsigned int ncharges);

  /**
   * Saves the gradients from backwards call on Torch model
   */
  virtual void save_output_gradients(std::ofstream& ifs
                                   , const unsigned int step);

  /**
   * Saves the point charge gradients from backwards call on Torch model
   */
  virtual void save_output_pc_gradients(std::ofstream& ifs
                                      , const unsigned int step);

  /**
   * Saves the charges calculated by Torch model (not implemented)
   */
  virtual void save_output_charges(std::ofstream& ifs
                                 , const unsigned int step);

  /**
   * Helper function to write the current step size
   */
  virtual void write_step_size(std::ofstream& ifs, 
                               const unsigned int step) const;

  /**
   * Helper function to write the header in coordinate files
   */
  virtual void write_coordinate_header(std::ofstream& ifs) const;

  /**
   * Helper function to write the footer in coordinate files
   */
  virtual void write_coordinate_footer(std::ofstream& ifs) const;

  /**
   * Helper function to write a single gradient to a file
   */
  virtual void write_gradient(const math::Vec& gradient, 
                              std::ofstream& inputfile_stream) const;

  /**
   * Helper function to write a single QM atom to a file
   */
  virtual void write_qm_atom(std::ofstream& inputfile_stream
                           , const int atomic_number
                           , const math::Vec& pos) const;

  /**
   * Helper function to write a single MM atom to a file
   */
  virtual void write_mm_atom(std::ofstream& inputfile_stream
                           , const int atomic_number
                           , const math::Vec& pos
                           , const double charge) const;

  /**
   * Helper function to write a single charge to a file (not completely implemented)
   */
  virtual void write_charge(std::ofstream& inputfile_stream
                          , const int atomic_number
                          , const double charge) const;

  /**
   * Helper function to open a file
   */
  virtual int open_input(std::ofstream& inputfile_stream, const std::string& input_file) const;
  
  /**
   * Computes the number of charges like in QM_Worker.cc
   * TODO: combine with QM_Worker.cc
   */
  virtual int get_num_charges(const simulation::Simulation &sim) const;

  /**
   * A copy of the QM zone, is refreshed every step
  */
  QM_Zone qm_zone;

  /**
   * A (non_owning) pointer to the QMMM interaction
   */
  QMMM_Interaction *qmmm_ptr;

  /**
   * The size of the QM zone (QM atoms + QM link atoms)
   */
  int natoms;

  /**
   * Number of point charges in the MM zone
   */
  int ncharges;

  /**
   * How many batches are sent to Torch each iteration
  */
  unsigned batch_size;

  /**
   * How many Cartesian coordinates
  */
  unsigned dimensions;

  /**
   * Handle to the input coordinate trajectory (Torch/QM) 
   */
  mutable std::ofstream input_coordinate_stream;

  /**
   * Handle to the input point charge trajectory (Torch/QM) 
   */
  mutable std::ofstream input_point_charge_stream;

  /**
   * Handle to the output gradient trajectory (Torch/QM) 
   */
  mutable std::ofstream output_gradient_stream;

  /**
   * Handle to the output point charge gradient trajectory (Torch/QM) 
   */
  mutable std::ofstream output_point_charge_gradient_stream;

  /**
   * Handle to the output charges trajectory (Torch/QM)
   * TODO: not implemented 
   */
  mutable std::ofstream output_charges_stream;

  /**
   * Atomic numbers of the QM zone as C style array
   */
  std::vector<int64_t> qm_atomic_numbers;

  /**
   * Atomic positions of the QM zone as C style array
   */
  std::vector<T> qm_positions;

  /**
   * Atomic numbers of the MM zone as C style array
   */
  std::vector<unsigned> mm_atomic_numbers;

  /**
   * Charges of the MM zone as C style array
   */
  std::vector<T> mm_charges;

  /**
   * Atomic positions of the MM zone as C style array
   */
  std::vector<T> mm_positions;

  /**
   * A (non-owning) tensor to hold QM atomic numbers
   */
  torch::Tensor qm_atomic_numbers_tensor;

  /**
   * A (non-owning) tensor to hold QM positions
   */
  torch::Tensor qm_positions_tensor;

  /**
   * A (non-owning) tensor to hold MM atomic positions
   */
  torch::Tensor mm_atomic_numbers_tensor;

  /**
   * A (non-owning) tensor to hold MM charges
   */
  torch::Tensor mm_charges_tensor;

  /**
   * A (non-owning) tensor to hold MM positions
   */
  torch::Tensor mm_positions_tensor;

  /**
   * A (non-owning) tensor to hold energies calculated
   */
  torch::Tensor energy_tensor;

  /**
   * A (non-owning) tensor to hold QM gradients calculated
   */
  torch::Tensor qm_gradient_tensor;

  /**
   * A (non-owning) tensor to hold MM gradients calculated
   */
  torch::Tensor mm_gradient_tensor;
};

} // namespace interaction

#endif /* INCLUDED_TORCH_QMMM_INTERACTION_H */