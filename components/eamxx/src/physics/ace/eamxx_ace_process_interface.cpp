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
  m_ll_grid = grids_manager->get_grid("LatLonPhysics");
  m_pt_grid = grids_manager->get_grid("Physics");
  using namespace ShortFieldTagsNames;

  constexpr auto m      = ekat::units::m;
  constexpr auto s      = ekat::units::s;
  constexpr auto K      = ekat::units::K;
  constexpr auto nondim = ekat::units::Units::nondimensional();
  
  // Get grid dimensions from the grid
  // Get the layout of the grid to determine dimensions
  auto layout_2d = m_ll_grid->get_2d_scalar_layout();
  auto layout_3d = m_ll_grid->get_3d_scalar_layout(true); // true = midpoints
  // m_ll_grid->get_3d_vector_layout();
  // Extract dimensions from layouts
  int m_height = layout_2d.dim(0);  // First dimension is LAT
  int m_width = layout_2d.dim(1);  // Second dimension is LON
  int m_levels = layout_3d.dim(2);  // Third dimension in 3D layout is LEV
  
  // Define layouts
  const auto layout_3d_vector = m_ll_grid->get_3d_vector_layout(true,2);
  const auto horiz_wind_layout = FieldLayout({CMP, LAT, LON}, {2, m_height, m_width});
  const auto layout_3d_scalar = FieldLayout({LAT, LON, LEV}, {m_height, m_width, m_levels});
  const auto layout_2d_scalar = FieldLayout({LAT, LON}, {m_height, m_width});

  // Stuff we will use from the import layer
  add_field<Updated>("T_mid", layout_3d_scalar, K, m_ll_grid->name());
  add_field<Updated>("horiz_winds", horiz_wind_layout, m / s, m_ll_grid->name());
  add_field<Updated>("ocnfrac", layout_2d_scalar, nondim, m_ll_grid->name());
  add_field<Required>("landfrac", layout_2d_scalar, nondim, m_ll_grid->name());
  add_field<Updated>("icefrac", layout_2d_scalar, nondim, m_ll_grid->name());

  // Stuff we will provide to the export layer
  
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

  // add_field<Computed>("sfc_alb_dir_vis",  scalar2d, nondim,  grid_name);
  // add_field<Computed>("sfc_alb_dir_nir",  scalar2d, nondim,  grid_name);
  // add_field<Computed>("sfc_alb_dif_vis",  scalar2d, nondim,  grid_name);
  // add_field<Computed>("sfc_alb_dif_nir",  scalar2d, nondim,  grid_name);
  // add_field<Computed>("surf_lw_flux_up",  scalar2d, W/m2,    grid_name);
  // add_field<Computed>("surf_sens_flux",   scalar2d, W/m2,    grid_name);
  // add_field<Computed>("surf_evap",        scalar2d, kg/m2/s, grid_name);
  // add_field<Computed>("surf_mom_flux",    vector2d, N/m2,    grid_name);
  // add_field<Computed>("surf_radiative_T", scalar2d, K,       grid_name);
  // add_field<Computed>("T_2m",             scalar2d, K,       grid_name);
  // add_field<Computed>("qv_2m",            scalar2d, kg/kg,   grid_name);
  // add_field<Computed>("wind_speed_10m",   scalar2d, m/s,     grid_name);
  // add_field<Computed>("snow_depth_land",  scalar2d, m,       grid_name);
  // add_field<Computed>("ocnfrac",          scalar2d, nondim,  grid_name);
  // add_field<Computed>("landfrac",         scalar2d, nondim,  grid_name);
  // add_field<Computed>("icefrac",          scalar2d, nondim,  grid_name);
  // // Friction velocity [m/s]
  // add_field<Computed>("fv",               scalar2d, m/s,     grid_name);
  // // Aerodynamical resistance
  // add_field<Computed>("ram1",             scalar2d, s/m,     grid_name);
  // // Sea surface temperature [K]
  // add_field<Computed>("sst",              scalar2d, K,       grid_name);
  // //dust fluxes [kg/m^2/s]: Four flux values for eacch column
  // add_field<Computed>("dstflx",           vector4d, kg/m2/s, grid_name);

}
void ACE::postprocess () {

  /* 
    MUST provide these to export layer:
  */

  // // These fields are required for computation/exports
  // add_field<Required>("p_int",                scalar3d_layout_int,  Pa,     grid_name);
  // add_field<Required>("pseudo_density",       scalar3d_layout_mid,  Pa,     grid_name, ps);
  // add_field<Required>("phis",                 scalar2d_layout,      m2/s2,  grid_name);
  // add_field<Required>("p_mid",                scalar3d_layout_mid,  Pa,     grid_name, ps);
  // add_field<Required>("T_mid",                scalar3d_layout_mid,  K,      grid_name, ps);
  // add_tracer<Required>("qv", m_grid,  kg/kg, ps);
  // // TODO: Switch horiz_winds to using U and V, note right now there is an issue with when the subfields are created, so can't switch yet.
  // add_field<Required>("horiz_winds",          vector3d_layout,      m/s,    grid_name);
  // add_field<Required>("sfc_flux_dir_nir",     scalar2d_layout,      W/m2,   grid_name);
  // add_field<Required>("sfc_flux_dir_vis",     scalar2d_layout,      W/m2,   grid_name);
  // add_field<Required>("sfc_flux_dif_nir",     scalar2d_layout,      W/m2,   grid_name);
  // add_field<Required>("sfc_flux_dif_vis",     scalar2d_layout,      W/m2,   grid_name);
  // add_field<Required>("sfc_flux_sw_net" ,     scalar2d_layout,      W/m2,   grid_name);
  // add_field<Required>("sfc_flux_lw_dn"  ,     scalar2d_layout,      W/m2,   grid_name);
  // add_field<Required>("precip_liq_surf_mass", scalar2d_layout,      kg/m2,  grid_name);
  // add_field<Required>("precip_ice_surf_mass", scalar2d_layout,      kg/m2,  grid_name);

}

void ACE::run_impl(const double dt) {
  auto device = (*m_module.parameters().begin()).device();

  // TODO: this is just to test mechanics of stuff, but will be replaced later
  // InField will follow the mechanics in the comment after the code (bottom)
  // Get a view from an input field
  auto in_field_view = get_field_in("InField").get_view<const Real ****>();

  // TODO move to namelist? or make them members?
  const int batch = 1, in_ch = 39, out_ch = 44, height = 180, width = 360;

  torch::Tensor m_torch_input = torch::from_blob(
      (void *)in_field_view.data(), {m_batch, m_in_ch, height, width},
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
      kokkos_output(m_torch_output.data_ptr<float>(), m_batch, m_out_ch, height, width);

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
