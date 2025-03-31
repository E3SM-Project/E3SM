#include "eamxx_ace_process_interface.hpp"

#include <fstream>
#include <iostream>
#include <string>

#include "physics/share/physics_constants.hpp"
#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"

namespace scream {

ACE::ACE(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  m_checkpoint_path = params.get<std::string>("checkpoint_path");
  m_sc2ace_map = m_params.get<std::string>("sc2ace_mapfile");
  m_ace2sc_map = m_params.get<std::string>("ace2sc_mapfile");
  m_ace_print_verification = m_params.get<bool>("ace_print_verification", true);
}

void ACE::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  m_ll_grid = grids_manager->get_grid("AceLL");
  m_pt_grid = grids_manager->get_grid("AcePG");
  m_sc_grid = grids_manager->get_grid("Physics PG2");

  using namespace ShortFieldTagsNames;

  constexpr auto m = ekat::units::m;
  constexpr auto s = ekat::units::s;
  constexpr auto K = ekat::units::K;
  constexpr auto nondim = ekat::units::Units::nondimensional();
  constexpr auto Pa = ekat::units::Pa;
  constexpr auto W = ekat::units::W;
  constexpr auto kg = ekat::units::kg;
  constexpr auto N = ekat::units::N;
  auto m2 = m * m;
  auto s2 = s * s;

  // Get grid dimensions from the sc grid
  // Get the layout of the grid to determine dimensions

  auto sc_layout_2d_scalar = m_sc_grid->get_2d_scalar_layout();
  auto sc_layout_3d_scalar = m_sc_grid->get_3d_scalar_layout(true);  // true = midpoints
  auto sc_layout_3d_scalar_int = m_sc_grid->get_3d_scalar_layout(false);  // false = interfaces
  auto sc_layout_3d_vector = m_sc_grid->get_3d_vector_layout(true, 2);
  auto sc_layout_2d_vector = m_sc_grid->get_2d_vector_layout(2);

  // Stuff sc_export EXPECTS from ACE
  // NOTE: "Updated" seems to be necessary for the ll-sc syncing in initialize_impl
  // TODO: ask Luca to take a look into this?
  add_field<Updated>("p_int", sc_layout_3d_scalar_int, Pa, m_sc_grid->name());
  add_field<Updated>("pseudo_density", sc_layout_3d_scalar, Pa, m_sc_grid->name());
  add_field<Updated>("phis", sc_layout_2d_scalar, m2 / s2, m_sc_grid->name());
  add_field<Updated>("p_mid", sc_layout_3d_scalar, Pa, m_sc_grid->name());
  add_field<Updated>("T_mid", sc_layout_3d_scalar, K, m_sc_grid->name());
  add_field<Updated>("qv", sc_layout_3d_scalar, kg / kg, m_sc_grid->name());
  add_field<Updated>("horiz_winds", sc_layout_3d_vector, m / s, m_sc_grid->name());
  add_field<Updated>("sfc_flux_dir_nir", sc_layout_2d_scalar, W / m2, m_sc_grid->name());
  add_field<Updated>("sfc_flux_dir_vis", sc_layout_2d_scalar, W / m2, m_sc_grid->name());
  add_field<Updated>("sfc_flux_dif_nir", sc_layout_2d_scalar, W / m2, m_sc_grid->name());
  add_field<Updated>("sfc_flux_dif_vis", sc_layout_2d_scalar, W / m2, m_sc_grid->name());
  add_field<Updated>("sfc_flux_sw_net", sc_layout_2d_scalar, W / m2, m_sc_grid->name());
  add_field<Updated>("sfc_flux_lw_dn", sc_layout_2d_scalar, W / m2, m_sc_grid->name());
  add_field<Updated>("precip_liq_surf_mass", sc_layout_2d_scalar, kg / m2, m_sc_grid->name());
  add_field<Updated>("precip_ice_surf_mass", sc_layout_2d_scalar, kg / m2, m_sc_grid->name());

  // Stuff ACE can (but doesn't have to) use from the sc_import
  // NOTE: "Updated" seems to be necessary for the ll-sc syncing in initialize_impl
  // TODO: ask Luca to take a look into this?
  add_field<Updated>("sfc_alb_dir_vis", sc_layout_2d_scalar, nondim, m_sc_grid->name());
  add_field<Updated>("sfc_alb_dir_nir", sc_layout_2d_scalar, nondim, m_sc_grid->name());
  add_field<Updated>("sfc_alb_dif_vis", sc_layout_2d_scalar, nondim, m_sc_grid->name());
  add_field<Updated>("sfc_alb_dif_nir", sc_layout_2d_scalar, nondim, m_sc_grid->name());
  add_field<Updated>("surf_lw_flux_up", sc_layout_2d_scalar, W / m2, m_sc_grid->name());
  add_field<Updated>("surf_sens_flux", sc_layout_2d_scalar, W / m2, m_sc_grid->name());
  add_field<Updated>("surf_evap", sc_layout_2d_scalar, kg / m2 / s, m_sc_grid->name());
  add_field<Updated>("surf_mom_flux", sc_layout_2d_vector, N / m2, m_sc_grid->name());
  add_field<Updated>("surf_radiative_T", sc_layout_2d_scalar, K, m_sc_grid->name());
  add_field<Updated>("T_2m", sc_layout_2d_scalar, K, m_sc_grid->name());
  add_field<Updated>("qv_2m", sc_layout_2d_scalar, kg / kg, m_sc_grid->name());
  add_field<Updated>("wind_speed_10m", sc_layout_2d_scalar, m / s,m_sc_grid->name());
  add_field<Updated>("snow_depth_land", sc_layout_2d_scalar, m,m_sc_grid->name());
  add_field<Updated>("ocnfrac", sc_layout_2d_scalar, nondim, m_sc_grid->name());
  add_field<Updated>("landfrac", sc_layout_2d_scalar, nondim,m_sc_grid->name());
  add_field<Updated>("icefrac", sc_layout_2d_scalar, nondim, m_sc_grid->name());

  // Create two big fields for input and output for ease (can remove later)
  FieldLayout scalar3d_39x180x360_mid{{CMP, CMP, LAT, LON}, {m_batch, m_in_ch, m_height, m_width}};
  FieldLayout scalar3d_44x180x360_mid{{CMP, CMP, LAT, LON}, {m_batch, m_out_ch, m_height, m_width}};
  m_input_field = Field(FieldIdentifier("InField", scalar3d_39x180x360_mid, nondim, m_ll_grid->name()));
  m_output_field = Field(FieldIdentifier("OutField", scalar3d_44x180x360_mid, nondim, m_ll_grid->name()));
  m_input_field.allocate_view();
  m_output_field.allocate_view();
}

void ACE::initialize_impl(const RunType /* run_type */) {
  // Output configuration values
  std::cout << "ACE::initialize_impl: Checkpoint Path: " << m_checkpoint_path << std::endl;

  // Check if the checkpoint file exists
  std::ifstream cp_file(m_checkpoint_path, std::ios::binary);
  if(!cp_file.good()) {
    std::cerr << "ACE::initialize_impl: Checkpoint file not found: " << m_checkpoint_path << "\n";
    return;
  }
  cp_file.close();

  std::cout << "ACE::initialize_impl: Loading model checkpoint from: " << m_checkpoint_path << std::endl;
  try {
    m_module = torch::jit::load(m_checkpoint_path);
  } catch(const c10::Error &e) {
    std::cerr << "ACE::initialize_impl: Error loading the model checkpoint: " << e.what() << "\n";
    std::cerr << "ACE::initialize_impl: Verify that the checkpoint was properly exported and "
                 "includes 'constants.pkl'.\n";
    return;
  }

  // Initialize m_input_tensor, m_output_tensor here
  // TODO: in the future share underlying memory of fields instead of allocating tensors
  // TODO: maybe helpful to explicitly set device too
  m_input_tensor = torch::zeros(
      {m_batch, m_in_ch, m_height, m_width},
      torch::TensorOptions().dtype(torch::kFloat32).requires_grad(false));
  m_output_tensor = torch::zeros(
      {m_batch, m_out_ch, m_height, m_width},
      torch::TensorOptions().dtype(torch::kFloat32).requires_grad(false));

  // A lot of workarounds follow to get remapping going
  std::cout << "ACE::initialize_impl: setting up remappers" << std::endl;
  
  // Same as the set_grids, but here for pt grid instead of sc grid above
  auto pt_layout_2d_scalar = m_pt_grid->get_2d_scalar_layout();
  auto pt_layout_3d_scalar = m_pt_grid->get_3d_scalar_layout(true);  // true = midpoints
  auto pt_layout_3d_scalar_int = m_pt_grid->get_3d_scalar_layout(false);  // false = interfaces
  auto pt_layout_3d_vector = m_pt_grid->get_3d_vector_layout(true, 2);
  auto pt_layout_2d_vector = m_pt_grid->get_2d_vector_layout(2);

  // NOTE: both remappers are coarsening by construction
  // TODO: EAMxx remappers don't support remap_bwd() (could've made life easier...)
  m_sc2ace_remapper = std::make_shared<CoarseningRemapper>(m_sc_grid, m_sc2ace_map);
  m_ace2sc_remapper = std::make_shared<CoarseningRemapper>(m_pt_grid, m_ace2sc_map);

  // Register the fields for remapping, creating helper fields
  m_sc2ace_remapper->registration_begins();
  m_ace2sc_remapper->registration_begins();
  
  create_helper_field("p_int", pt_layout_3d_scalar_int, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("p_int"), get_helper_field("p_int"));
  m_ace2sc_remapper->register_field(get_helper_field("p_int"), get_field_out("p_int"));
  
  create_helper_field("pseudo_density", pt_layout_3d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("pseudo_density"), get_helper_field("pseudo_density"));
  m_ace2sc_remapper->register_field(get_helper_field("pseudo_density"), get_field_out("pseudo_density"));
  
  create_helper_field("phis", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("phis"), get_helper_field("phis"));
  m_ace2sc_remapper->register_field(get_helper_field("phis"), get_field_out("phis"));

  create_helper_field("p_mid", pt_layout_3d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("p_mid"), get_helper_field("p_mid"));
  m_ace2sc_remapper->register_field(get_helper_field("p_mid"), get_field_out("p_mid"));

  create_helper_field("T_mid", pt_layout_3d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("T_mid"), get_helper_field("T_mid"));
  m_ace2sc_remapper->register_field(get_helper_field("T_mid"), get_field_out("T_mid"));
  
  create_helper_field("qv", pt_layout_3d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("qv"), get_helper_field("qv"));
  m_ace2sc_remapper->register_field(get_helper_field("qv"), get_field_out("qv"));

  create_helper_field("horiz_winds", pt_layout_3d_vector, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("horiz_winds"), get_helper_field("horiz_winds"));
  m_ace2sc_remapper->register_field(get_helper_field("horiz_winds"), get_field_out("horiz_winds"));

  create_helper_field("sfc_flux_dir_nir", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("sfc_flux_dir_nir"), get_helper_field("sfc_flux_dir_nir"));
  m_ace2sc_remapper->register_field(get_helper_field("sfc_flux_dir_nir"), get_field_out("sfc_flux_dir_nir"));

  create_helper_field("sfc_flux_dir_vis", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("sfc_flux_dir_vis"), get_helper_field("sfc_flux_dir_vis"));
  m_ace2sc_remapper->register_field(get_helper_field("sfc_flux_dir_vis"), get_field_out("sfc_flux_dir_vis"));

  create_helper_field("sfc_flux_dif_nir", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("sfc_flux_dif_nir"), get_helper_field("sfc_flux_dif_nir"));
  m_ace2sc_remapper->register_field(get_helper_field("sfc_flux_dif_nir"), get_field_out("sfc_flux_dif_nir"));

  create_helper_field("sfc_flux_dif_vis", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("sfc_flux_dif_vis"), get_helper_field("sfc_flux_dif_vis"));

  create_helper_field("sfc_flux_sw_net", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("sfc_flux_sw_net"), get_helper_field("sfc_flux_sw_net"));
  m_ace2sc_remapper->register_field(get_helper_field("sfc_flux_sw_net"), get_field_out("sfc_flux_sw_net"));

  create_helper_field("sfc_flux_lw_dn", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("sfc_flux_lw_dn"), get_helper_field("sfc_flux_lw_dn"));
  m_ace2sc_remapper->register_field(get_helper_field("sfc_flux_lw_dn"), get_field_out("sfc_flux_lw_dn"));

  create_helper_field("precip_liq_surf_mass", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("precip_liq_surf_mass"), get_helper_field("precip_liq_surf_mass"));
  m_ace2sc_remapper->register_field(get_helper_field("precip_liq_surf_mass"), get_field_out("precip_liq_surf_mass"));

  create_helper_field("precip_ice_surf_mass", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("precip_ice_surf_mass"), get_helper_field("precip_ice_surf_mass"));
  m_ace2sc_remapper->register_field(get_helper_field("precip_ice_surf_mass"), get_field_out("precip_ice_surf_mass"));

  create_helper_field("sfc_alb_dir_vis", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("sfc_alb_dir_vis"), get_helper_field("sfc_alb_dir_vis"));
  m_ace2sc_remapper->register_field(get_helper_field("sfc_alb_dir_vis"), get_field_out("sfc_alb_dir_vis"));

  create_helper_field("sfc_alb_dir_nir", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("sfc_alb_dir_nir"), get_helper_field("sfc_alb_dir_nir"));
  m_ace2sc_remapper->register_field(get_helper_field("sfc_alb_dir_nir"), get_field_out("sfc_alb_dir_nir"));

  create_helper_field("sfc_alb_dif_vis", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("sfc_alb_dif_vis"), get_helper_field("sfc_alb_dif_vis"));
  m_ace2sc_remapper->register_field(get_helper_field("sfc_alb_dif_vis"), get_field_out("sfc_alb_dif_vis"));

  create_helper_field("sfc_alb_dif_nir", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("sfc_alb_dif_nir"), get_helper_field("sfc_alb_dif_nir"));
  m_ace2sc_remapper->register_field(get_helper_field("sfc_alb_dif_nir"), get_field_out("sfc_alb_dif_nir"));

  create_helper_field("surf_lw_flux_up", pt_layout_2d_scalar, m_pt_grid->name()); 
  m_sc2ace_remapper->register_field(get_field_out("surf_lw_flux_up"), get_helper_field("surf_lw_flux_up"));
  m_ace2sc_remapper->register_field(get_helper_field("surf_lw_flux_up"), get_field_out("surf_lw_flux_up"));

  create_helper_field("surf_sens_flux", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("surf_sens_flux"), get_helper_field("surf_sens_flux"));
  m_ace2sc_remapper->register_field(get_helper_field("surf_sens_flux"), get_field_out("surf_sens_flux"));

  create_helper_field("surf_evap", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("surf_evap"), get_helper_field("surf_evap"));
  m_ace2sc_remapper->register_field(get_helper_field("surf_evap"), get_field_out("surf_evap"));

  create_helper_field("surf_mom_flux", pt_layout_2d_vector, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("surf_mom_flux"), get_helper_field("surf_mom_flux"));
  m_ace2sc_remapper->register_field(get_helper_field("surf_mom_flux"), get_field_out("surf_mom_flux"));

  create_helper_field("surf_radiative_T", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("surf_radiative_T"), get_helper_field("surf_radiative_T"));
  m_ace2sc_remapper->register_field(get_helper_field("surf_radiative_T"), get_field_out("surf_radiative_T"));
  
  create_helper_field("T_2m", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("T_2m"), get_helper_field("T_2m"));
  m_ace2sc_remapper->register_field(get_helper_field("T_2m"), get_field_out("T_2m"));
  
  create_helper_field("qv_2m", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("qv_2m"), get_helper_field("qv_2m"));
  m_ace2sc_remapper->register_field(get_helper_field("qv_2m"), get_field_out("qv_2m"));

  create_helper_field("wind_speed_10m", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("wind_speed_10m"), get_helper_field("wind_speed_10m"));
  m_ace2sc_remapper->register_field(get_helper_field("wind_speed_10m"), get_field_out("wind_speed_10m"));

  create_helper_field("snow_depth_land", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("snow_depth_land"), get_helper_field("snow_depth_land"));
  m_ace2sc_remapper->register_field(get_helper_field("snow_depth_land"), get_field_out("snow_depth_land"));
  
  create_helper_field("ocnfrac", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("ocnfrac"), get_helper_field("ocnfrac"));
  m_ace2sc_remapper->register_field(get_helper_field("ocnfrac"), get_field_out("ocnfrac"));
  
  create_helper_field("landfrac", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("landfrac"), get_helper_field("landfrac"));
  m_ace2sc_remapper->register_field(get_helper_field("landfrac"), get_field_out("landfrac"));
  
  create_helper_field("icefrac", pt_layout_2d_scalar, m_pt_grid->name());
  m_sc2ace_remapper->register_field(get_field_out("icefrac"), get_helper_field("icefrac"));
  m_ace2sc_remapper->register_field(get_helper_field("icefrac"), get_field_out("icefrac"));

  m_sc2ace_remapper->registration_ends();
  m_ace2sc_remapper->registration_ends();

  // The follow is to ensure both the ll (lat-lon grid) and pt (a point-grid version thereof)
  // share the underlying data of the fields once remapped...

  // Same units as set_grids
  constexpr auto m = ekat::units::m;
  constexpr auto s = ekat::units::s;
  constexpr auto K = ekat::units::K;
  constexpr auto nondim = ekat::units::Units::nondimensional();
  constexpr auto Pa = ekat::units::Pa;
  constexpr auto W = ekat::units::W;
  constexpr auto kg = ekat::units::kg;
  auto m2 = m * m;
  auto s2 = s * s;

  // Same layouts as before, but for ll grid
  auto ll_layout_2d_scalar = m_ll_grid->get_2d_scalar_layout();
  auto ll_layout_3d_scalar = m_ll_grid->get_3d_scalar_layout(true);  // true = midpoints
  auto ll_layout_3d_scalar_int = m_ll_grid->get_3d_scalar_layout(false);  // false = interfaces
  auto ll_layout_3d_vector = m_ll_grid->get_3d_vector_layout(true, 2);

  // Get views of all helper fields created
  auto land_fraction_pt_view = get_helper_field("landfrac").get_view<Real *>();
  auto ocean_fraction_pt_view = get_helper_field("ocnfrac").get_view<Real *>();
  auto sea_ice_fraction_pt_view = get_helper_field("icefrac").get_view<Real *>();
  auto surf_radiative_T_pt_view = get_helper_field("surf_radiative_T").get_view<Real *>();
  auto T_mid_pt_view = get_helper_field("T_mid").get_view<Real **>();
  auto qv_pt_view = get_helper_field("qv").get_view<Real **>();
  auto wind_pt_view = get_helper_field("horiz_winds").get_view<Real ***>();
  auto p_mid_view = get_helper_field("p_mid").get_view<Real **>();
  auto p_int_view = get_helper_field("p_int").get_view<Real **>();
  auto dp_mid_view = get_helper_field("pseudo_density").get_view<Real **>();
  auto phis_view = get_helper_field("phis").get_view<Real *>();
  auto sfc_flux_dir_nir_view = get_helper_field("sfc_flux_dir_nir").get_view<Real *>();
  auto sfc_flux_dir_vis_view = get_helper_field("sfc_flux_dir_vis").get_view<Real *>();
  auto sfc_flux_dif_nir_view = get_helper_field("sfc_flux_dif_nir").get_view<Real *>();
  auto sfc_flux_dif_vis_view = get_helper_field("sfc_flux_dif_vis").get_view<Real *>();
  auto sfc_flux_sw_net_view = get_helper_field("sfc_flux_sw_net").get_view<Real *>();
  auto sfc_flux_lw_dn_view = get_helper_field("sfc_flux_lw_dn").get_view<Real *>();
  auto precip_liq_surf_mass_view = get_helper_field("precip_liq_surf_mass").get_view<Real *>();
  auto precip_ice_surf_mass_view = get_helper_field("precip_ice_surf_mass").get_view<Real *>();

  // Created FID for all fields with correct layouts and info
  FieldIdentifier land_fraction_ll_fid("landfrac_ll", ll_layout_2d_scalar, nondim, m_ll_grid->name());
  FieldIdentifier ocean_fraction_ll_fid("ocnfrac_ll", ll_layout_2d_scalar, nondim, m_ll_grid->name());
  FieldIdentifier sea_ice_fraction_ll_fid("icefrac_ll", ll_layout_2d_scalar, nondim, m_ll_grid->name());
  FieldIdentifier surf_radiative_T_ll_fid("surf_radiative_T_ll", ll_layout_2d_scalar, K, m_ll_grid->name());
  FieldIdentifier T_mid_ll_fid("T_mid_ll", ll_layout_3d_scalar, K, m_ll_grid->name());
  FieldIdentifier qv_ll_fid("qv_ll", ll_layout_3d_scalar, kg / kg, m_ll_grid->name());
  FieldIdentifier wind_ll_fid("horiz_winds_ll", ll_layout_3d_vector, m / s, m_ll_grid->name());
  FieldIdentifier p_mid_ll_fid("p_mid_ll", ll_layout_3d_scalar, Pa, m_ll_grid->name());
  FieldIdentifier p_int_ll_fid("p_int_ll", ll_layout_3d_scalar_int, Pa, m_ll_grid->name());
  FieldIdentifier dp_mid_ll_fid("pseudo_density_ll", ll_layout_3d_scalar, Pa, m_ll_grid->name());
  FieldIdentifier phis_ll_fid("phis_ll", ll_layout_2d_scalar, m2 / s2, m_ll_grid->name());
  FieldIdentifier sfc_flux_dir_nir_fid("sfc_flux_dir_nir_ll", ll_layout_2d_scalar, W / m2, m_ll_grid->name());
  FieldIdentifier sfc_flux_dir_vis_fid("sfc_flux_dir_vis_ll", ll_layout_2d_scalar, W / m2, m_ll_grid->name());
  FieldIdentifier sfc_flux_dif_nir_fid("sfc_flux_dif_nir_ll", ll_layout_2d_scalar, W / m2, m_ll_grid->name());
  FieldIdentifier sfc_flux_dif_vis_fid("sfc_flux_dif_vis_ll", ll_layout_2d_scalar, W / m2, m_ll_grid->name());
  FieldIdentifier sfc_flux_sw_net_fid("sfc_flux_sw_net_ll", ll_layout_2d_scalar, W / m2, m_ll_grid->name());
  FieldIdentifier sfc_flux_lw_dn_fid("sfc_flux_lw_dn_ll", ll_layout_2d_scalar, W / m2, m_ll_grid->name());
  FieldIdentifier precip_liq_surf_mass_fid("precip_liq_surf_mass_ll", ll_layout_2d_scalar, kg / m2, m_ll_grid->name());
  FieldIdentifier precip_ice_surf_mass_fid("precip_ice_surf_mass_ll", ll_layout_2d_scalar, kg / m2, m_ll_grid->name());

  // Properly create the view with correct dimension ordering
  // noting how how views increasing by +1 dimension (so, COLxLEV --> LATxLONxLEV, etc.)
  // TODO: does this construction necessitate the above views be
  // TODO: ask Luca why this is so and maybe fix?
  Field::view_dev_t<Real **> land_fraction_ll_v(land_fraction_pt_view.data(), m_height, m_width);
  Field::view_dev_t<Real **> ocean_fraction_ll_v(ocean_fraction_pt_view.data(), m_height, m_width);
  Field::view_dev_t<Real **> sea_ice_fraction_ll_v(sea_ice_fraction_pt_view.data(), m_height, m_width);
  Field::view_dev_t<Real **> surf_radiative_T_ll_v(surf_radiative_T_pt_view.data(), m_height, m_width);
  Field::view_dev_t<Real ***> T_mid_ll_v(T_mid_pt_view.data(), m_height, m_width, m_levels);
  Field::view_dev_t<Real ***> qv_ll_v(qv_pt_view.data(), m_height, m_width, m_levels);
  Field::view_dev_t<Real ****> wind_ll_v(wind_pt_view.data(), m_height, m_width, 2, m_levels);
  Field::view_dev_t<Real ***> p_mid_ll_v(p_mid_view.data(), m_height, m_width, m_levels);
  Field::view_dev_t<Real ***> p_int_ll_v(p_int_view.data(), m_height, m_width, m_levels + 1);
  Field::view_dev_t<Real ***> dp_mid_ll_v(dp_mid_view.data(), m_height, m_width, m_levels);
  Field::view_dev_t<Real **> phis_ll_v(phis_view.data(), m_height, m_width);
  Field::view_dev_t<Real **> sfc_flux_dir_nir_v(sfc_flux_dir_nir_view.data(), m_height, m_width);
  Field::view_dev_t<Real **> sfc_flux_dir_vis_v(sfc_flux_dir_vis_view.data(), m_height, m_width);
  Field::view_dev_t<Real **> sfc_flux_dif_nir_v(sfc_flux_dif_nir_view.data(), m_height, m_width);
  Field::view_dev_t<Real **> sfc_flux_dif_vis_v(sfc_flux_dif_vis_view.data(), m_height, m_width);
  Field::view_dev_t<Real **> sfc_flux_sw_net_v(sfc_flux_sw_net_view.data(), m_height, m_width);
  Field::view_dev_t<Real **> sfc_flux_lw_dn_v(sfc_flux_lw_dn_view.data(), m_height, m_width);
  Field::view_dev_t<Real **> precip_liq_surf_mass_v(precip_liq_surf_mass_view.data(), m_height, m_width);
  Field::view_dev_t<Real **> precip_ice_surf_mass_v(precip_ice_surf_mass_view.data(), m_height, m_width);

  // Okay, ready to roll
  m_land_fraction_ll = Field(land_fraction_ll_fid, land_fraction_ll_v);
  m_ocean_fraction_ll = Field(ocean_fraction_ll_fid, ocean_fraction_ll_v);
  m_sea_ice_fraction_ll = Field(sea_ice_fraction_ll_fid, sea_ice_fraction_ll_v);
  m_surf_radiative_T_ll = Field(surf_radiative_T_ll_fid, surf_radiative_T_ll_v);
  m_T_mid_ll = Field(T_mid_ll_fid, T_mid_ll_v);
  m_qv_ll = Field(qv_ll_fid, qv_ll_v);
  m_wind_ll = Field(wind_ll_fid, wind_ll_v);
  m_p_mid_ll = Field(p_mid_ll_fid, p_mid_ll_v);
  m_p_int_ll = Field(p_int_ll_fid, p_int_ll_v);
  m_dp_mid_ll = Field(dp_mid_ll_fid, dp_mid_ll_v);
  m_phis_ll = Field(phis_ll_fid, phis_ll_v);
  m_sfc_flux_dir_nir = Field(sfc_flux_dir_nir_fid, sfc_flux_dir_nir_v);
  m_sfc_flux_dir_vis = Field(sfc_flux_dir_vis_fid, sfc_flux_dir_vis_v);
  m_sfc_flux_dif_nir = Field(sfc_flux_dif_nir_fid, sfc_flux_dif_nir_v);
  m_sfc_flux_dif_vis = Field(sfc_flux_dif_vis_fid, sfc_flux_dif_vis_v);
  m_sfc_flux_sw_net = Field(sfc_flux_sw_net_fid, sfc_flux_sw_net_v);
  m_sfc_flux_lw_dn = Field(sfc_flux_lw_dn_fid, sfc_flux_lw_dn_v);
  m_precip_liq_surf_mass = Field(precip_liq_surf_mass_fid, precip_liq_surf_mass_v);
  m_precip_ice_surf_mass = Field(precip_ice_surf_mass_fid, precip_ice_surf_mass_v);
}

void ACE::preprocess() {
  /*

    The idea is to preprocess the fields into an organized data structure
    so that it can be input into the emulator, whose design is such that
    we pass to it 39 channels comprised of 39 fields, as described below.
    First, the sc_import layer provides ACE with some fields it may use.

    -----------------------------------------------------
    MAY make use of provided fields from sc_import layer:
    -----------------------------------------------------
    add_field<Computed>("sfc_alb_dir_vis",  scalar2d, nondim,  grid_name);
    add_field<Computed>("sfc_alb_dir_nir",  scalar2d, nondim,  grid_name);
    add_field<Computed>("sfc_alb_dif_vis",  scalar2d, nondim,  grid_name);
    add_field<Computed>("sfc_alb_dif_nir",  scalar2d, nondim,  grid_name);
    add_field<Computed>("surf_lw_flux_up",  scalar2d, W/m2,    grid_name);
    add_field<Computed>("surf_sens_flux",   scalar2d, W/m2,    grid_name);
    add_field<Computed>("surf_evap",        scalar2d, kg/m2/s, grid_name);
    add_field<Computed>("surf_mom_flux",    vector2d, N/m2,    grid_name);
    add_field<Computed>("surf_radiative_T", scalar2d, K,       grid_name);
    add_field<Computed>("T_2m",             scalar2d, K,       grid_name);
    add_field<Computed>("qv_2m",            scalar2d, kg/kg,   grid_name);
    add_field<Computed>("wind_speed_10m",   scalar2d, m/s,     grid_name);
    add_field<Computed>("snow_depth_land",  scalar2d, m,       grid_name);
    add_field<Computed>("ocnfrac",          scalar2d, nondim,  grid_name);
    add_field<Computed>("landfrac",         scalar2d, nondim,  grid_name);
    add_field<Computed>("icefrac",          scalar2d, nondim,  grid_name);
    // Friction velocity [m/s]
    add_field<Computed>("fv",               scalar2d, m/s,     grid_name);
    // Aerodynamical resistance
    add_field<Computed>("ram1",             scalar2d, s/m,     grid_name);
    // Sea surface temperature [K]
    add_field<Computed>("sst",              scalar2d, K,       grid_name);
    //dust fluxes [kg/m^2/s]: Four flux values for eacch column
    add_field<Computed>("dstflx",           vector4d, kg/m2/s, grid_name);

    ----------------------------
    ACE emulator requires these:
    ----------------------------
    - land_fraction
    - ocean_fraction
    - sea_ice_fraction
    - DSWRFtoa # solin
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
  */

  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  auto input_view = m_input_field.get_view<Real ****>();

  auto land_fraction = m_land_fraction_ll.get_view<const Real **>();
  auto ocean_fraction = m_ocean_fraction_ll.get_view<const Real **>();
  auto sea_ice_fraction = m_sea_ice_fraction_ll.get_view<const Real **>();
  auto surf_radiative_T = m_surf_radiative_T_ll.get_view<const Real **>();
  auto T_mid = m_T_mid_ll.get_view<const Real ***>();
  auto qv = m_qv_ll.get_view<const Real ***>();
  auto wind = m_wind_ll.get_view<const Real ****>();
  auto p_int = m_p_int_ll.get_view<const Real ***>();

  const auto num_lat = m_height;
  const auto policy = ESU::get_default_team_policy(m_width, m_height);

  // Populate the view of input data
  // TODO: bypass this and just fill in the tensor directly?
  // ALSO: note hacks!
  Kokkos::parallel_for(
      "ACE::preprocess", policy, KOKKOS_LAMBDA(const MT &team) {
        const int ilon = team.league_rank();
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, num_lat), [&](const int &ilat) {
              input_view(0, 0, ilat, ilon) = land_fraction(ilat, ilon);
              input_view(0, 1, ilat, ilon) = ocean_fraction(ilat, ilon);
              input_view(0, 2, ilat, ilon) = sea_ice_fraction(ilat, ilon);
              // DSWRFtoa (SOLIN): TODD (BEN???)
              input_view(0, 3, ilat, ilon) = 500.0;
              // HGTsfc (surface height): TODO
              input_view(0, 4, ilat, ilon) = 10.0;
              input_view(0, 5, ilat, ilon) = p_int(ilat, ilon, 8);
              input_view(0, 6, ilat, ilon) = surf_radiative_T(ilat, ilon);
              input_view(0, 7, ilat, ilon) = T_mid(ilat, ilon, 0);
              input_view(0, 8, ilat, ilon) = T_mid(ilat, ilon, 1);
              input_view(0, 9, ilat, ilon) = T_mid(ilat, ilon, 2);
              input_view(0, 10, ilat, ilon) = T_mid(ilat, ilon, 3);
              input_view(0, 11, ilat, ilon) = T_mid(ilat, ilon, 4);
              input_view(0, 12, ilat, ilon) = T_mid(ilat, ilon, 5);
              input_view(0, 13, ilat, ilon) = T_mid(ilat, ilon, 6);
              input_view(0, 14, ilat, ilon) = T_mid(ilat, ilon, 7);
              input_view(0, 15, ilat, ilon) = qv(ilat, ilon, 0);
              input_view(0, 16, ilat, ilon) = qv(ilat, ilon, 1);
              input_view(0, 17, ilat, ilon) = qv(ilat, ilon, 2);
              input_view(0, 18, ilat, ilon) = qv(ilat, ilon, 3);
              input_view(0, 19, ilat, ilon) = qv(ilat, ilon, 4);
              input_view(0, 20, ilat, ilon) = qv(ilat, ilon, 5);
              input_view(0, 21, ilat, ilon) = qv(ilat, ilon, 6);
              input_view(0, 22, ilat, ilon) = qv(ilat, ilon, 7);
              input_view(0, 23, ilat, ilon) = wind(ilat, ilon, 0, 0);
              input_view(0, 24, ilat, ilon) = wind(ilat, ilon, 0, 1);
              input_view(0, 25, ilat, ilon) = wind(ilat, ilon, 0, 2);
              input_view(0, 26, ilat, ilon) = wind(ilat, ilon, 0, 3);
              input_view(0, 27, ilat, ilon) = wind(ilat, ilon, 0, 4);
              input_view(0, 28, ilat, ilon) = wind(ilat, ilon, 0, 5);
              input_view(0, 29, ilat, ilon) = wind(ilat, ilon, 0, 6);
              input_view(0, 30, ilat, ilon) = wind(ilat, ilon, 0, 7);
              input_view(0, 31, ilat, ilon) = wind(ilat, ilon, 1, 0);
              input_view(0, 32, ilat, ilon) = wind(ilat, ilon, 1, 1);
              input_view(0, 33, ilat, ilon) = wind(ilat, ilon, 1, 2);
              input_view(0, 34, ilat, ilon) = wind(ilat, ilon, 1, 3);
              input_view(0, 35, ilat, ilon) = wind(ilat, ilon, 1, 4);
              input_view(0, 36, ilat, ilon) = wind(ilat, ilon, 1, 5);
              input_view(0, 37, ilat, ilon) = wind(ilat, ilon, 1, 6);
              input_view(0, 38, ilat, ilon) = wind(ilat, ilon, 1, 7);
            });
      });
}

void ACE::run_impl(const double dt) {
  // The main ACE run loop takes place here
  std::cout << "ACE::run_impl: sc2ace remapping" << std::endl;
  m_sc2ace_remapper->remap_fwd();
  std::cout << "ACE::run_impl: ACE preprocessing" << std::endl;
  preprocess();

  // Detect model device and ensure tensor is on the same device
  torch::Device model_device = torch::kCPU;
  for (const auto& param : m_module.parameters()) {
    model_device = param.device();
    break;  // We only need to check the first parameter
  }
  
  std::cout << "ACE::run_impl: Model is on " 
            << (model_device.is_cuda() ? "CUDA" : "CPU") 
            << " device" << std::endl;

  auto input_view = m_input_field.get_view<Real ****>();
  
  // we probably don't want to support dtype flexibility ... 
  // let's force everyone to use kfloat32 and be done with it
  // we should maybe even use lower precision...
  // anyway, this will likely fail below (in the torch call, due to conflict with saved )
  auto torch_dtype = std::is_same<Real, double>::value ? torch::kFloat64 : torch::kFloat32;

  // the tensor uses the same pointer as the big field view
  m_input_tensor = torch::from_blob(
    (void *)input_view.data(), {m_batch, m_in_ch, m_height, m_width},
    torch::TensorOptions()
        .dtype(torch_dtype)
        .device(model_device)
        .memory_format(torch::MemoryFormat::Contiguous));

  // likely unneeded, but just in case...
  m_input_tensor = m_input_tensor.to(model_device);
  
  std::cout << "ACE::run_impl: ACE inference" << std::endl;
  std::vector<torch::jit::IValue> inputs{m_input_tensor};
  try { 
    // TODO: this call creates and returns a tensor
    // TODO: investigate how to make it recycle a view from above   
    m_output_tensor = m_module.forward(inputs).toTensor();
  } catch(const c10::Error &e) {
    std::cerr << "ACE::run_impl: Error during inference: " << e.what() << std::endl;
    return;
  }

  // TODO: is there a way to avoid the forthcoming deep_copy?
  // Careful here, need to ensure LayoutRight! Torch::Tensor is LR by default!
  Kokkos::View<Real ****, Kokkos::LayoutRight, scream::DefaultDevice,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>>
       kokkos_output(m_output_tensor.data_ptr<Real>(), m_batch, m_out_ch, m_height, m_width);
  auto output_view = m_output_field.get_view<Real ****>();
  Kokkos::deep_copy(output_view, kokkos_output);

  std::cout << "ACE::run_impl: ACE postprocessing" << std::endl;
  postprocess();

  if (m_ace_print_verification) {
    std::cout << "ACE::run_impl: ACE verification" << std::endl;
    m_output_field.sync_to_host();
    auto out_field_view_2 = m_output_field.get_view<const Real ****, Host>();
    auto out_tensor_cpu   = m_output_tensor.to(torch::kCPU);

    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        for(int k = 0; k < 3; k++) {
          std::cout << "out_fied[0][" << i << "][" << j << "][" << k
                    << "] = " << out_field_view_2(0, i, j, k) << std::endl;
          std::cout << "out_tens[0][" << i << "][" << j << "][" << k
                    << "] = " << out_tensor_cpu[0][i][j][k].item<double>()
                    << std::endl;
        }
      }
    }
  }

  std::cout << "ACE::run_impl: ace2sc remapping" << std::endl;
  m_ace2sc_remapper->remap_fwd();

}

void ACE::postprocess() {
  /*

    The idea is to postprocess the fields into an organized data structure
    after it comes out of the emulator, whose design is such that
    we get 44 channels comprised of 44 fields, as described below.
    First, the sc_export layer EXPECTS ACE to provide it some fileds.

    -------------------------------------------------
    MUST provide required fields for sc_export layer:
    -------------------------------------------------
    add_field<Required>("p_int", scalar3d_layout_int, Pa, grid_name);
    add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa, grid_name, ps);
    add_field<Required>("phis", scalar2d_layout, m2/s2, grid_name);
    add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name, ps);
    add_field<Required>("T_mid", scalar3d_layout_mid, K, grid_name, ps);
    add_tracer<Required>("qv", m_grid, kg/kg, ps);
    add_field<Required>("horiz_winds", vector3d_layout, m/s, grid_name);
    add_field<Required>("sfc_flux_dir_nir", scalar2d_layout, W/m2, grid_name);
    add_field<Required>("sfc_flux_dir_vis", scalar2d_layout, W/m2, grid_name);
    add_field<Required>("sfc_flux_dif_nir", scalar2d_layout, W/m2, grid_name);
    add_field<Required>("sfc_flux_dif_vis", scalar2d_layout, W/m2, grid_name);
    add_field<Required>("sfc_flux_sw_net", scalar2d_layout, W/m2, grid_name);
    add_field<Required>("sfc_flux_lw_dn", scalar2d_layout, W/m2, grid_name);
    add_field<Required>("precip_liq_surf_mass", scalar2d_layout, kg/m2, grid_name);
    add_field<Required>("precip_ice_surf_mass", scalar2d_layout, kg/m2, grid_name);

    ---------------------------
    ACE emulator outputs these:
    ---------------------------
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
    - LHTFLsfc // latent_heat_flux
    - SHTFLsfc // sensible_heat_flux
    - PRATEsfc // precipitation_rate
    - ULWRFsfc // sfc_up_lw_radiative_flux
    - ULWRFtoa // toa_up_lw_radiative_flux
    - DLWRFsfc // sfc_down_lw_radiative_flux
    - DSWRFsfc // sfc_down_sw_radiative_flux
    - USWRFsfc // sfc_up_sw_radiative_flux
    - USWRFtoa // toa_up_sw_radiative_flux
    - tendency_of_total_water_path_due_to_advection
  */

  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  auto output_view = m_output_field.get_view<const Real ****>();

  auto T_mid = m_T_mid_ll.get_view<Real ***>();
  auto qv = m_qv_ll.get_view<Real ***>();
  auto wind = m_wind_ll.get_view<Real ****>();
  auto p_mid = m_p_mid_ll.get_view<Real ***>();
  auto p_int = m_p_int_ll.get_view<Real ***>();
  auto dp_mid = m_dp_mid_ll.get_view<Real ***>();
  auto phis = m_phis_ll.get_view<Real **>();
  auto sfc_flux_dir_nir = m_sfc_flux_dir_nir.get_view<Real **>();
  auto sfc_flux_dir_vis = m_sfc_flux_dir_vis.get_view<Real **>();
  auto sfc_flux_dif_nir = m_sfc_flux_dif_nir.get_view<Real **>();
  auto sfc_flux_dif_vis = m_sfc_flux_dif_vis.get_view<Real **>();
  auto sfc_flux_sw_net = m_sfc_flux_sw_net.get_view<Real **>();
  auto sfc_flux_lw_dn = m_sfc_flux_lw_dn.get_view<Real **>();
  auto precip_liq_surf_mass = m_precip_liq_surf_mass.get_view<Real **>();
  auto precip_ice_surf_mass = m_precip_ice_surf_mass.get_view<Real **>();

  const auto num_lat = m_height;
  const auto policy = ESU::get_default_team_policy(m_width, m_height);

  // Populate the field views
  // TODO: bypass this and just fill in the tensor directly?
  // ALSO: note hackery and trickery!!!
  Kokkos::parallel_for(
      "ACE::postprocess", policy, KOKKOS_LAMBDA(const MT &team) {
        const int ilon = team.league_rank();
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, num_lat), [&](const int &ilat) {
              // TODO: is this even legit conceptually? Taking dp and p_int etc.
              // TODO: second, copied these ak and bk values from an ai2/ace file
              // TODO: formula is: p_i = ak_i + bk_i * PRESsfc
              // std::map<int, std::pair<Real, Real>> m_ak_bk_values = {
              //   {0, {64.247, 0.0}},
              //   {1, {5167.14603, 0.0}},
              //   {2, {12905.42546, 0.01755}},
              //   {3, {13982.4677, 0.11746}},
              //   {4, {12165.28766, 0.2896}},
              //   {5, {8910.07678, 0.49806}},
              //   {6, {4955.72632, 0.72625}},
              //   {7, {2155.78385, 0.88192}}
              // };
              p_int(ilat, ilon, 8) = output_view(0, 0, ilat, ilon);
              p_mid(ilat, ilon, 0) = p_int(ilat, ilon, 8) * 0.0 + 64.247;
              p_int(ilat, ilon, 0) = (0.0 + p_mid(ilat, ilon, 0)) / 2;
              dp_mid(ilat, ilon, 0) = p_mid(ilat, ilon, 0) - p_int(ilat, ilon, 0);
              p_mid(ilat, ilon, 1) = p_int(ilat, ilon, 8) * 0.0 + 5167.14603;
              p_int(ilat, ilon, 1) = (p_mid(ilat, ilon, 0) + p_mid(ilat, ilon, 1)) / 2;
              dp_mid(ilat, ilon, 1) = p_mid(ilat, ilon, 1) - p_int(ilat, ilon, 1);
              p_mid(ilat, ilon, 2) = p_int(ilat, ilon, 8) * 0.01755 + 12905.42546;
              p_int(ilat, ilon, 2) = (p_mid(ilat, ilon, 1) + p_mid(ilat, ilon, 2)) / 2;
              dp_mid(ilat, ilon, 2) = p_mid(ilat, ilon, 2) - p_int(ilat, ilon, 2);
              p_mid(ilat, ilon, 3) = p_int(ilat, ilon, 8) * 0.11746 + 13982.4677;
              p_int(ilat, ilon, 3) = (p_mid(ilat, ilon, 2) + p_mid(ilat, ilon, 3)) / 2;
              dp_mid(ilat, ilon, 3) = p_mid(ilat, ilon, 3) - p_int(ilat, ilon, 3);
              p_mid(ilat, ilon, 4) = p_int(ilat, ilon, 8) * 0.2896 + 12165.28766;
              p_int(ilat, ilon, 4) = (p_mid(ilat, ilon, 3) + p_mid(ilat, ilon, 4)) / 2;
              dp_mid(ilat, ilon, 4) = p_mid(ilat, ilon, 4) - p_int(ilat, ilon, 4);
              p_mid(ilat, ilon, 5) = p_int(ilat, ilon, 8) * 0.49806 + 8910.07678;
              p_int(ilat, ilon, 5) = (p_mid(ilat, ilon, 4) + p_mid(ilat, ilon, 5)) / 2;
              dp_mid(ilat, ilon, 5) = p_mid(ilat, ilon, 5) - p_int(ilat, ilon, 5);
              p_mid(ilat, ilon, 6) = p_int(ilat, ilon, 8) * 0.72625 + 4955.72632;
              p_int(ilat, ilon, 6) = (p_mid(ilat, ilon, 5) + p_mid(ilat, ilon, 6)) / 2;
              dp_mid(ilat, ilon, 6) = p_mid(ilat, ilon, 6) - p_int(ilat, ilon, 6);
              p_mid(ilat, ilon, 7) = p_int(ilat, ilon, 8) * 0.88192 + 2155.78385;
              p_int(ilat, ilon, 7) = (p_mid(ilat, ilon, 6) + p_mid(ilat, ilon, 7)) / 2;
              dp_mid(ilat, ilon, 7) = p_mid(ilat, ilon, 7) - p_int(ilat, ilon, 7);
              phis(ilat, ilon) = 1.0;
              // Surface temperature ???
              // surface_temperature = output_view(0, 1, ilat, ilon);
              T_mid(ilat, ilon, 0) = output_view(0, 2, ilat, ilon);
              T_mid(ilat, ilon, 1) = output_view(0, 3, ilat, ilon);
              T_mid(ilat, ilon, 2) = output_view(0, 4, ilat, ilon);
              T_mid(ilat, ilon, 3) = output_view(0, 5, ilat, ilon);
              T_mid(ilat, ilon, 4) = output_view(0, 6, ilat, ilon);
              T_mid(ilat, ilon, 5) = output_view(0, 7, ilat, ilon);
              T_mid(ilat, ilon, 6) = output_view(0, 8, ilat, ilon);
              T_mid(ilat, ilon, 7) = output_view(0, 9, ilat, ilon);
              qv(ilat, ilon, 0) = output_view(0, 10, ilat, ilon);
              qv(ilat, ilon, 1) = output_view(0, 11, ilat, ilon);
              qv(ilat, ilon, 2) = output_view(0, 12, ilat, ilon);
              qv(ilat, ilon, 3) = output_view(0, 13, ilat, ilon);
              qv(ilat, ilon, 4) = output_view(0, 14, ilat, ilon);
              qv(ilat, ilon, 5) = output_view(0, 15, ilat, ilon);
              qv(ilat, ilon, 6) = output_view(0, 16, ilat, ilon);
              qv(ilat, ilon, 7) = output_view(0, 17, ilat, ilon);
              wind(ilat, ilon, 0, 0) = output_view(0, 18, ilat, ilon);
              wind(ilat, ilon, 0, 1) = output_view(0, 19, ilat, ilon);
              wind(ilat, ilon, 0, 2) = output_view(0, 20, ilat, ilon);
              wind(ilat, ilon, 0, 3) = output_view(0, 21, ilat, ilon);
              wind(ilat, ilon, 0, 4) = output_view(0, 22, ilat, ilon);
              wind(ilat, ilon, 0, 5) = output_view(0, 23, ilat, ilon);
              wind(ilat, ilon, 0, 6) = output_view(0, 24, ilat, ilon);
              wind(ilat, ilon, 0, 7) = output_view(0, 25, ilat, ilon);
              wind(ilat, ilon, 1, 0) = output_view(0, 26, ilat, ilon);
              wind(ilat, ilon, 1, 1) = output_view(0, 27, ilat, ilon);
              wind(ilat, ilon, 1, 2) = output_view(0, 28, ilat, ilon);
              wind(ilat, ilon, 1, 3) = output_view(0, 29, ilat, ilon);
              wind(ilat, ilon, 1, 4) = output_view(0, 30, ilat, ilon);
              wind(ilat, ilon, 1, 5) = output_view(0, 31, ilat, ilon);
              wind(ilat, ilon, 1, 6) = output_view(0, 32, ilat, ilon);
              wind(ilat, ilon, 1, 7) = output_view(0, 33, ilat, ilon);
              // LHTFLsfc // latent_heat_flux
              // LHTFLsfc            = output_view(0, 34, ilat, ilon);
              // SHTFLsfc // sensible_heat_flux
              // SHTFLsfc            = output_view(0, 35, ilat, ilon);
              // PRATEsfc // precipitation_rate
              precip_liq_surf_mass(ilat, ilon) = output_view(0, 36, ilat, ilon);
              precip_ice_surf_mass(ilat, ilon) = 0.0;
              // ULWRFsfc // sfc_up_lw_radiative_flux
              // ULWRFsfc            = output_view(0, 37, ilat, ilon);
              // ULWRFtoa // toa_up_lw_radiative_flux
              // ULWRFtoa            = output_view(0, 38, ilat, ilon);
              // DLWRFsfc // sfc_down_lw_radiative_flux
              sfc_flux_lw_dn(ilat, ilon) = output_view(0, 39, ilat, ilon);
              // DSWRFsfc // sfc_down_sw_radiative_flux
              // DSWRFsfc            = output_view(0, 40, ilat, ilon);
              // USWRFsfc // sfc_up_sw_radiative_flux
              // USWRFsfc            = output_view(0, 41, ilat, ilon);
              sfc_flux_sw_net(ilat, ilon) = output_view(0, 40, ilat, ilon) - output_view(0, 41, ilat, ilon);
              // USWRFtoa // toa_up_sw_radiative_flux
              // USWRFtoa            = output_view(0, 42, ilat, ilon);
              // ???
              // tendency_of_total_water_path_due_to_advection = output_view(0, 43, ilat, ilon);
              // For the lack of a better alternative, fill
              // dir_vis, etc. with random fractions of net surface balances
              // TODO: HACK!
              sfc_flux_dir_vis(ilat, ilon) = 0.4 * sfc_flux_sw_net(ilat, ilon);
              sfc_flux_dir_nir(ilat, ilon) = 0.3 * sfc_flux_sw_net(ilat, ilon);
              sfc_flux_dif_vis(ilat, ilon) = 0.2 * sfc_flux_sw_net(ilat, ilon);
              sfc_flux_dif_nir(ilat, ilon) = 0.1 * sfc_flux_sw_net(ilat, ilon);
            });
      });
}

void ACE::finalize_impl() {}

Field ACE::create_helper_field(const std::string &name,
                               const FieldLayout &layout,
                               const std::string &grid_name) {
  using namespace ekat::units;

  // For helper fields we don't bother w/ units
  // so we set them to non-dimensional
  // HACK: to get it to name them different names...
  auto my_name = name + "_helper";
  FieldIdentifier id(my_name, layout, Units::nondimensional(), grid_name);

  // Create the field.
  // Init with NaN's, so we spot instances of uninited memory usage
  Field f(id);
  f.get_header().get_alloc_properties().request_allocation(1);
  f.allocate_view();
  f.deep_copy(ekat::ScalarTraits<Real>::invalid());

  m_helper_fields[my_name] = f;
  return m_helper_fields[my_name];
}

}  // namespace scream
