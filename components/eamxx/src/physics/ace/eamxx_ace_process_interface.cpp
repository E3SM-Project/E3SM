#include "eamxx_ace_process_interface.hpp"

#include <fstream>
#include <iostream>
#include <string>

namespace scream {

ACE::ACE(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  m_forward_steps           = params.get<int>("forward_steps");
  m_forward_steps_in_memory = params.get<int>("forward_steps_in_memory");
  m_checkpoint_path         = params.get<std::string>("checkpoint_path");
}

void ACE::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  m_grid = grids_manager->get_grid("Physics");

  constexpr auto Pa     = ekat::units::Pa;
  constexpr auto kg     = ekat::units::kg;
  constexpr auto m      = ekat::units::m;
  constexpr auto s      = ekat::units::s;
  constexpr auto K      = ekat::units::K;
  constexpr auto m2     = m * m;
  constexpr auto nondim = ekat::units::nondimensional();

  const auto horiz_wind_layout =
      m_grid->get_2d_vector_layout(2);  // 2D horizontal wind layout
  const auto layout_3d = m_grid->get_3d_scalar_layout(true);
  const auto layout_2d = m_grid->get_2d_scalar_layout();

  add_field<Updated>("T_mid", layout_3d, K, m_grid->name());
  add_field<Updated>("horiz_winds", horiz_wind_layout, m / s, m_grid->name());
  add_field<Updated>("ocnfrac", scalar2d, nondim, grid_name);
  add_field<Updated>("landfrac", scalar2d, nondim, grid_name);
  add_field<Updated>("icefrac", scalar2d, nondim, grid_name);
}

void ACE::initialize_impl(const RunType /* run_type */) {
  // Output configuration values
  std::cout << "Number of Forward Steps: " << m_forward_steps << std::endl;
  std::cout << "Forward Steps in Memory: " << m_forward_steps_in_memory
            << std::endl;
  std::cout << "Checkpoint Path: " << m_checkpoint_path << std::endl;

  // Check if the checkpoint file exists
  std::ifstream cp_file(m_checkpoint_path, std::ios::binary);
  if(!cp_file.good()) {
    std::cerr << "Checkpoint file not found: " << m_checkpoint_path << "\n";
    return;
  }
  cp_file.close();

  std::cout << "Loading model checkpoint from: " << m_checkpoint_path
            << std::endl;
  try {
    m_module = torch::jit::load(m_checkpoint_path);
  } catch(const c10::Error &e) {
    std::cerr << "Error loading the model checkpoint: " << e.what() << "\n";
    std::cerr << "Verify that the checkpoint was properly exported and "
                 "includes 'constants.pkl'.\n";
    return;
  }
}

void ACE::run_impl(const double dt) {
  auto device = (*m_module.parameters().begin()).device();
  // TODO: Replace with actual input data
  // For now, we are using random data for testing purposes
  // The input tensor shape is (batch_size, num_channels, height, width)
  // Where batch_size is 1, num_channels is 39, height is 180, and width is 360
  // See comment on bottom for mechanics
  torch::Tensor input = torch::rand({1, 39, 180, 360});
  input               = input.to(device);
  std::vector<torch::jit::IValue> inputs{input};
  std::cout << "Input tensor shape: " << input.sizes() << std::endl;
  std::cout << "Performing inference for " << m_forward_steps << " steps..."
            << std::endl;
  torch::Tensor output = m_module.forward(inputs).toTensor();
  std::cout << "Inference output shape: " << output.sizes() << std::endl;
  std::cout << "Inference completed. Results stored in: " << m_checkpoint_path
            << std::endl;
}

void ACE::finalize_impl() {}

}  // namespace scream

/*

 Additional comments:

 In the tensor, the 39 "channels" correspond to the following in_names
 When the emulator runs, it returns 44 "channels" corresponding to the out_names
 We thus need to get these from the native model state (e.g., T_mid is
 temperature) Remap onto the Gaussian 180x360 grid and remap onto the 8-level
 vertical grid

  in_names:
  - land_fraction
  - ocean_fraction
  - sea_ice_fraction
  - DSWRFtoa
  - HGTsfc
  - PRESsfc
  - surface_temperature
  - air_temperature_0 # _0 denotes the top most layer of the atmosphere
  - air_temperature_1
  - air_temperature_2
  - air_temperature_3
  - air_temperature_4
  - air_temperature_5
  - air_temperature_6
  - air_temperature_7
  - specific_total_water_0
  - specific_total_water_1
  - specific_total_water_2
  - specific_total_water_3
  - specific_total_water_4
  - specific_total_water_5
  - specific_total_water_6
  - specific_total_water_7
  - eastward_wind_0
  - eastward_wind_1
  - eastward_wind_2
  - eastward_wind_3
  - eastward_wind_4
  - eastward_wind_5
  - eastward_wind_6
  - eastward_wind_7
  - northward_wind_0
  - northward_wind_1
  - northward_wind_2
  - northward_wind_3
  - northward_wind_4
  - northward_wind_5
  - northward_wind_6
  - northward_wind_7

  out_names:
  - PRESsfc
  - surface_temperature
  - air_temperature_0
  - air_temperature_1
  - air_temperature_2
  - air_temperature_3
  - air_temperature_4
  - air_temperature_5
  - air_temperature_6
  - air_temperature_7
  - specific_total_water_0
  - specific_total_water_1
  - specific_total_water_2
  - specific_total_water_3
  - specific_total_water_4
  - specific_total_water_5
  - specific_total_water_6
  - specific_total_water_7
  - eastward_wind_0
  - eastward_wind_1
  - eastward_wind_2
  - eastward_wind_3
  - eastward_wind_4
  - eastward_wind_5
  - eastward_wind_6
  - eastward_wind_7
  - northward_wind_0
  - northward_wind_1
  - northward_wind_2
  - northward_wind_3
  - northward_wind_4
  - northward_wind_5
  - northward_wind_6
  - northward_wind_7
  - LHTFLsfc
  - SHTFLsfc
  - PRATEsfc
  - ULWRFsfc
  - ULWRFtoa
  - DLWRFsfc
  - DSWRFsfc
  - USWRFsfc
  - USWRFtoa
  - tendency_of_total_water_path_due_to_advection

  Additionally, the emulator calculates some diagnostics by default
*/
