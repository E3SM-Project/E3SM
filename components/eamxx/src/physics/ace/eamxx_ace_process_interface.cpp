#include "eamxx_ace_process_interface.hpp"

#include <fstream>
#include <iostream>
#include <string>

namespace scream {

ACE::ACE(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  m_checkpoint_path         = params.get<std::string>("checkpoint_path");
}

void ACE::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  m_ll_grid = grids_manager->get_grid("AceLL");
  m_pt_grid = grids_manager->get_grid("AcePG");
  m_sc_grid = grids_manager->get_grid("Physics");
  using namespace ShortFieldTagsNames;

  constexpr auto m      = ekat::units::m;
  constexpr auto s      = ekat::units::s;
  constexpr auto K      = ekat::units::K;
  constexpr auto nondim = ekat::units::Units::nondimensional();
  constexpr auto Pa    = ekat::units::Pa;
  constexpr auto W     = ekat::units::W;
  constexpr auto kg    = ekat::units::kg;
  constexpr auto N = ekat::units::N;
  auto m2 = m*m;
  auto s2 = s*s;
  
  // Get grid dimensions from the grid
  // Get the layout of the grid to determine dimensions

  auto sc_layout_2d_scalar = m_sc_grid->get_2d_scalar_layout();
  auto sc_layout_3d_scalar = m_sc_grid->get_3d_scalar_layout(true); // true = midpoints
  auto sc_layout_3d_scalar_int = m_sc_grid->get_3d_scalar_layout(false); // false = interfaces
  auto sc_layout_3d_vector = m_sc_grid->get_3d_vector_layout(true,2);
  auto sc_layout_2d_vector = m_sc_grid->get_2d_vector_layout(2);

  // Stuff sc_export expects from us
  add_field<Computed>("p_int",                sc_layout_3d_scalar_int,  Pa,     m_sc_grid->name());
  add_field<Computed>("pseudo_density",       sc_layout_3d_scalar,  Pa,     m_sc_grid->name());
  add_field<Computed>("phis",                 sc_layout_2d_scalar,      m2/s2,  m_sc_grid->name());
  add_field<Computed>("p_mid",                sc_layout_3d_scalar,  Pa,     m_sc_grid->name());
  add_field<Computed>("T_mid",                sc_layout_3d_scalar,  K,      m_sc_grid->name());
  add_field<Computed>("qv",                   sc_layout_3d_scalar,  kg/kg, m_sc_grid->name());
  add_field<Computed>("horiz_winds",          sc_layout_3d_vector,      m/s,    m_sc_grid->name());
  add_field<Computed>("sfc_flux_dir_nir",     sc_layout_2d_scalar,      W/m2,   m_sc_grid->name());
  add_field<Computed>("sfc_flux_dir_vis",     sc_layout_2d_scalar,      W/m2,   m_sc_grid->name());
  add_field<Computed>("sfc_flux_dif_nir",     sc_layout_2d_scalar,      W/m2,   m_sc_grid->name());
  add_field<Computed>("sfc_flux_dif_vis",     sc_layout_2d_scalar,      W/m2,   m_sc_grid->name());
  add_field<Computed>("sfc_flux_sw_net" ,     sc_layout_2d_scalar,      W/m2,   m_sc_grid->name());
  add_field<Computed>("sfc_flux_lw_dn"  ,     sc_layout_2d_scalar,      W/m2,   m_sc_grid->name());
  add_field<Computed>("precip_liq_surf_mass", sc_layout_2d_scalar,      kg/m2,  m_sc_grid->name());
  add_field<Computed>("precip_ice_surf_mass", sc_layout_2d_scalar,      kg/m2,  m_sc_grid->name());

  // Stuff we can (but don't have to) use from the import layer
  add_field<Required>("sfc_alb_dir_vis",  sc_layout_2d_scalar, nondim,  m_sc_grid->name());
  add_field<Required>("sfc_alb_dir_nir",  sc_layout_2d_scalar, nondim,  m_sc_grid->name());
  add_field<Required>("sfc_alb_dif_vis",  sc_layout_2d_scalar, nondim,  m_sc_grid->name());
  add_field<Required>("sfc_alb_dif_nir",  sc_layout_2d_scalar, nondim,  m_sc_grid->name());
  add_field<Required>("surf_lw_flux_up",  sc_layout_2d_scalar, W/m2,    m_sc_grid->name());
  add_field<Required>("surf_sens_flux",   sc_layout_2d_scalar, W/m2,    m_sc_grid->name());
  add_field<Required>("surf_evap",        sc_layout_2d_scalar, kg/m2/s, m_sc_grid->name());
  add_field<Required>("surf_mom_flux",    sc_layout_2d_vector, N/m2,    m_sc_grid->name());
  add_field<Required>("surf_radiative_T", sc_layout_2d_scalar, K,       m_sc_grid->name());
  add_field<Required>("T_2m",             sc_layout_2d_scalar, K,       m_sc_grid->name());
  add_field<Required>("qv_2m",            sc_layout_2d_scalar, kg/kg,   m_sc_grid->name());
  add_field<Required>("wind_speed_10m",   sc_layout_2d_scalar, m/s,     m_sc_grid->name());
  add_field<Required>("snow_depth_land",  sc_layout_2d_scalar, m,       m_sc_grid->name());
  add_field<Required>("ocnfrac",          sc_layout_2d_scalar, nondim,  m_sc_grid->name());
  add_field<Required>("landfrac",         sc_layout_2d_scalar, nondim,  m_sc_grid->name());
  add_field<Required>("icefrac",          sc_layout_2d_scalar, nondim,  m_sc_grid->name());

  FieldLayout scalar3d_39x180x360_mid{{CMP, CMP, LAT, LON}, {m_batch, m_in_ch, m_height, m_width}};
  FieldLayout scalar3d_44x180x360_mid{{CMP, CMP, LAT, LON}, {m_batch, m_out_ch, m_height, m_width}};

  add_field<Required>("InField", scalar3d_39x180x360_mid, nondim,
                      m_ll_grid->name());
  add_field<Computed>("OutField", scalar3d_44x180x360_mid, nondim,
                      m_ll_grid->name());
}

void ACE::initialize_impl(const RunType /* run_type */) {
  // Output configuration values
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

void ACE::preprocess () {

  /*
    MAY make use of these provided from the import layer:
  */


}
void ACE::postprocess () {

  /* 
    MUST provide these to export layer:
  */


}

void ACE::run_impl(const double dt) {
  auto device = (*m_module.parameters().begin()).device();

  // TODO: this is just to test mechanics of stuff, but will be replaced later
  // InField will follow the mechanics in the comment after the code (bottom)
  // Get a view from an input field
  auto in_field_view = get_field_in("InField").get_view<const Real ****>();

  torch::Tensor m_torch_input = torch::from_blob(
      (void *)in_field_view.data(), {m_batch, m_in_ch, m_height, m_width},
      torch::TensorOptions()
          .dtype(torch::kFloat32)
          .device(device)
          .memory_format(torch::MemoryFormat::Contiguous));

  // Debug prints
  std::cout << "in_tensor shape: " << m_torch_input.sizes() << std::endl;
  std::cout << "in_tensor device: " << m_torch_input.device() << std::endl;
  std::cout << "in_tensor dtype: " << m_torch_input.dtype() << std::endl;

  // Prepare for inference
  std::vector<torch::jit::IValue> inputs{m_torch_input};
  std::cout << "Input tensor shape: " << m_torch_input.sizes() << std::endl;

  torch::Tensor m_torch_output;
  try {
    m_torch_output = m_module.forward(inputs).toTensor().to(device);
  } catch(const c10::Error &e) {
    std::cerr << "Error during inference: " << e.what() << std::endl;
    return;
  }

  // Careful here, need to ensure LayoutRight! Torch::Tensor is LR by default!
  Kokkos::View<float ****, Kokkos::LayoutRight, scream::DefaultDevice,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      kokkos_output(m_torch_output.data_ptr<float>(), m_batch, m_out_ch, m_height, m_width);

  // TODO: this is just to test mechanics of stuff, but will be replaced later
  // OutField will follow the mechanics in the comment after the code (bottom)
  // Get a view from an output field
  auto out_field_view = get_field_out("OutField").get_view<Real ****>();
  Kokkos::deep_copy(out_field_view, kokkos_output);

  // verify that the values in OutField are the same as in m_torch_output?
  get_field_out("OutField").sync_to_host();
  auto out_field_view_2 =
      get_field_out("OutField").get_view<const Real ****, Host>();
  auto out_tensor_cpu = m_torch_output.to(torch::kCPU);

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      for(int k = 0; k < 3; k++) {
        std::cout << "out_fied[0][" << i << "][" << j << "][" << k
                  << "] = " << out_field_view_2(0, i, j, k) << std::endl;
        std::cout << "out_tens[0][" << i << "][" << j << "][" << k
                  << "] = " << out_tensor_cpu[0][i][j][k].item<float>()
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

*/
