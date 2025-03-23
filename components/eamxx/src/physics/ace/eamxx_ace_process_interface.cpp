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

  constexpr auto m      = ekat::units::m;
  constexpr auto s      = ekat::units::s;
  constexpr auto K      = ekat::units::K;
  constexpr auto nondim = ekat::units::Units::nondimensional();

  const auto horiz_wind_layout = m_grid->get_2d_vector_layout(2);
  const auto layout_3d         = m_grid->get_3d_scalar_layout(true);
  const auto layout_2d         = m_grid->get_2d_scalar_layout();

  add_field<Updated>("T_mid", layout_3d, K, m_grid->name());
  add_field<Updated>("horiz_winds", horiz_wind_layout, m / s, m_grid->name());
  add_field<Updated>("ocnfrac", layout_2d, nondim, m_grid->name());
  add_field<Required>("landfrac", layout_2d, nondim, m_grid->name());
  add_field<Updated>("icefrac", layout_2d, nondim, m_grid->name());

  // declare a layouts of 7x180x360
  using namespace ShortFieldTagsNames;
  FieldLayout scalar3d_39x180x360_mid{{CMP, CMP, CMP, CMP}, {1, 39, 180, 360}};
  FieldLayout scalar3d_44x180x360_mid{{CMP, CMP, CMP, CMP}, {1, 44, 180, 360}};

  add_field<Required>("InField", scalar3d_39x180x360_mid, nondim,
                      m_grid->name());
  add_field<Computed>("OutField", scalar3d_44x180x360_mid, nondim,
                      m_grid->name());
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
  // Get a view from an input field
  auto in_field_view = get_field_in("InField").get_view<const Real ****>();
  // Get a view from an output field
  auto out_field = get_field_out("OutField");
  out_field.deep_copy(0.0);
  auto out_field_view = out_field.get_view<Real ****>();

  auto device = (*m_module.parameters().begin()).device();

  // Create Torch tensors from Kokkos views
  torch::Tensor in_tensor = torch::from_blob((void *)in_field_view.data(),
                                             {1, 39, 180, 360}, torch::kFloat32)
                                .to(device);
  torch::Tensor out_tensor =
      torch::from_blob((void *)out_field_view.data(), {1, 44, 180, 360},
                       torch::kFloat32)
          .to(device);

  std::cout << "in_tensor shape: " << in_tensor.sizes() << std::endl;
  std::cout << "in_tensor device: " << in_tensor.device() << std::endl;
  std::cout << "in_tensor dtype: " << in_tensor.dtype() << std::endl;

  // Perform inference
  std::vector<torch::jit::IValue> inputs{in_tensor};
  std::cout << "Input tensor shape: " << in_tensor.sizes() << std::endl;
  std::cout << "Performing inference for " << m_forward_steps << " steps..."
            << std::endl;
  // torch::Tensor temp_out;
  try {
    torch::Tensor temp_out = m_module.forward(inputs).toTensor();
    std::memcpy(out_tensor.data_ptr<float>(), temp_out.data_ptr<float>(),
                temp_out.numel() * sizeof(float));
  } catch(const c10::Error &e) {
    std::cerr << "Error during inference: " << e.what() << std::endl;
    return;
  }

  // verify that the values in OutField are the same as in temp_out?
  out_field.sync_to_host();
  auto out_field_view_2 =
      get_field_out("OutField").get_view<const Real ****, Host>();
  auto out_tensor_cpu = out_tensor.to(torch::kCPU);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      for(int k = 0; k < 3; k++) {
        std::cout << "OutField[0][" << i << "][" << j << "][" << k
                  << "] = " << out_field_view_2(0, i, j, k) << std::endl;
        std::cout << "out_tensor[0][" << i << "][" << j << "][" << k
                  << "] = " << out_tensor[0][i][j][k].item<float>()
                  << std::endl;
      }
    }
  }
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
