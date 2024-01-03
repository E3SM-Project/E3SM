#include <physics/mam/eamxx_mam_optics_process_interface.hpp>
#include <share/property_checks/field_lower_bound_check.hpp>
#include <share/property_checks/field_within_interval_check.hpp>

#include "scream_config.h" // for SCREAM_CIME_BUILD
#include <ekat/ekat_assert.hpp>

#include "share/io/scorpio_input.hpp"
#include "share/grid/point_grid.hpp"

namespace scream
{

MAMOptics::MAMOptics(
    const ekat::Comm& comm,
    const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params),
    aero_config_() {
}

AtmosphereProcessType MAMOptics::type() const {
  return AtmosphereProcessType::Physics;
}

std::string MAMOptics::name() const {
  return "mam4_optics";
}

void MAMOptics::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;
  constexpr int nvars = mam4::ndrop::nvars;

  grid_ = grids_manager->get_grid("Physics");
  const auto& grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();      // number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels(); // number of levels per column
  nswbands_ = mam4::modal_aer_opt::nswbands; // number of shortwave bands
  nlwbands_ = mam4::modal_aer_opt::nlwbands; // number of longwave bands

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Define aerosol optics fields computed by this process.
  auto nondim = Units::nondimensional();
  //FieldLayout scalar3d_swband_layout { {COL, LEV, SWBND}, {ncol_, nlev_, nswbands_} };
  FieldLayout scalar3d_swbandp_layout { {COL, ILEV, SWBND}, {ncol_, nlev_+1, nswbands_} };
  FieldLayout scalar3d_lwband_layout { {COL, LEV, LWBND}, {ncol_, nlev_, nlwbands_} };
  FieldLayout scalar3d_layout_int{ {COL, ILEV}, {ncol_, nlev_+1} };

  // layout for 3D (2d horiz X 1d vertical) variables
  FieldLayout scalar3d_layout_mid{ {COL, LEV}, {ncol_, nlev_} };
  // FIXME: I switch the order of dimension.
  FieldLayout scalar_state_q_layout{ {COL, LEV, NGAS}, {ncol_, nlev_,nvars} };

  add_field<Required>("T_mid", scalar3d_layout_mid, K, grid_name); // Temperature
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name); // total pressure

  add_field<Required>("z_int", scalar3d_layout_int, m, grid_name); // vertical position at interface
  add_field<Required>("z_mid", scalar3d_layout_mid, m, grid_name); // vertical position pressure
  add_field<Required>("p_int", scalar3d_layout_int, Pa, grid_name); // total pressure
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa, grid_name);
  add_field<Required>("pseudo_density_dry", scalar3d_layout_mid, Pa, grid_name);
  // FIXME: dimenstion of state_q
  add_field<Updated>("state_q", scalar_state_q_layout, kg/kg, grid_name);
  add_field<Updated>("qqcw", scalar_state_q_layout, nondim, grid_name);
  add_field<Required>("cldn", scalar3d_layout_mid, nondim, grid_name); // could fraction

  // shortwave aerosol scattering asymmetry parameter [-]
  add_field<Computed>("aero_g_sw",   scalar3d_swbandp_layout, nondim, grid_name);
  // shortwave aerosol single-scattering albedo [-]
  add_field<Computed>("aero_ssa_sw", scalar3d_swbandp_layout, nondim, grid_name);
  // shortwave aerosol optical depth [-]
  add_field<Computed>("aero_tau_sw", scalar3d_swbandp_layout, nondim, grid_name);
  // longwave aerosol optical depth [-]
  add_field<Computed>("aero_tau_lw", scalar3d_lwband_layout, nondim, grid_name);
  // // aerosol extinction optical depth
  add_field<Computed>("aero_tau_forward", scalar3d_swbandp_layout, nondim, grid_name);

    // FIXME: this field doesn't belong here, but this is a convenient place to
  // FIXME: put it for now.
  // number mixing ratio for CCN
  using Spack      = ekat::Pack<Real,SCREAM_SMALL_PACK_SIZE>;
  using Pack       = ekat::Pack<Real,Spack::n>;
  constexpr int ps = Pack::n;
  // FieldLayout scalar3d_layout_mid { {COL, LEV}, {ncol_, nlev_} };
  add_field<Computed>("nccn", scalar3d_layout_mid, 1/kg, grid_name, ps);

add_field<Computed>("extinct",   scalar3d_layout_mid, 1/m, grid_name);
add_field<Computed>("absorb",   scalar3d_layout_mid, 1/m, grid_name);
//Aerosol optical depth 850 nm
FieldLayout scalar2d_layout_mid{ {COL}, {ncol_} };
add_field<Computed>("aodnir", scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("aoduv",  scalar2d_layout_mid, nondim, grid_name);

constexpr int ntot_amode = mam4::AeroConfig::num_modes();
//FIXME: I add NMODES to field_tag
FieldLayout scalar3d_layout_nmodes{ {COL, NMODES}, {ncol_, ntot_amode} };

add_field<Computed>("dustaodmode",  scalar3d_layout_nmodes, nondim, grid_name);
add_field<Computed>("aodmode",  scalar3d_layout_nmodes, nondim, grid_name);
add_field<Computed>("burdenmode",  scalar3d_layout_nmodes, nondim, grid_name);

add_field<Computed>("aodabsbc",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("aodvis",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("aodall",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("ssavis",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("aodabs",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("burdendust",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("burdenso4",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("burdenbc",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("burdenpom",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("burdensoa",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("burdenseasalt",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("burdenmom",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("momaod",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("dustaod",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("so4aod",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("pomaod",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("soaaod",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("bcaod",  scalar2d_layout_mid, nondim, grid_name);
add_field<Computed>("seasaltaod",  scalar2d_layout_mid, nondim, grid_name);

}

void MAMOptics::initialize_impl(const RunType run_type) {

    // FIXME: we have nvars in several process.
  constexpr int ntot_amode = mam4::AeroConfig::num_modes();

  dry_atm_.T_mid     = get_field_in("T_mid").get_view<const Real**>();
  dry_atm_.p_mid     = get_field_in("p_mid").get_view<const Real**>();
  // FIXME, there are two version of p_int in the nc file: p_dry_int and p_int
  // change to const Real
  p_int_     = get_field_in("p_int").get_view<const Real**>();

  dry_atm_.cldfrac   = get_field_in("cldn").get_view<const Real**>();
  // dry_atm_.pblh      = get_field_in("pbl_height").get_view<const Real*>();
  // FIXME: use const Real; why are using buffer in microphysics
  z_mid_     = get_field_in("z_mid").get_view<const Real**>();
  z_iface_   = get_field_in("z_int").get_view<const Real**>();

  p_del_     = get_field_in("pseudo_density").get_view<const Real**>();
  // FIXME: In the nc file, there is also pseudo_density_dry
  dry_atm_.p_del     = get_field_in("pseudo_density_dry").get_view<const Real**>();


  // FIXME: We need to get ssa_cmip6_sw_, af_cmip6_sw_, ext_cmip6_sw_, ext_cmip6_lw_ from a nc file.
  // aer_rad_props_sw inputs that are prescribed, i.e., we need a netcdf file.
  ssa_cmip6_sw_ = mam_coupling::view_3d ("ssa_cmip6_sw", ncol_, nlev_, nswbands_);
  af_cmip6_sw_ = mam_coupling::view_3d ("af_cmip6_sw", ncol_, nlev_, nswbands_);
  ext_cmip6_sw_= mam_coupling::view_3d ("ext_cmip6_sw", ncol_, nswbands_, nlev_);
  ext_cmip6_lw_ = mam_coupling::view_3d("ext_cmip6_lw_", ncol_, nlev_, nlwbands_);

  const int work_len = mam4::modal_aer_opt::get_work_len_aerosol_optics();
  work_ = mam_coupling::view_2d("work", ncol_, work_len);

  Kokkos::deep_copy(ssa_cmip6_sw_, 0.0);
  Kokkos::deep_copy(af_cmip6_sw_, 0.0);
  Kokkos::deep_copy(ext_cmip6_sw_, 0.0);
  Kokkos::deep_copy(ext_cmip6_lw_, 0.0);
 // read table info
{
  using namespace ShortFieldTagsNames;

  using view_1d_host = typename KT::view_1d<Real>::HostMirror;
  // Set up input structure to read data from file.
  // using strvec_t = std::vector<std::string>;
    // Views in aerosol_optics_device_data_ are allocated in the following fuctions.
  // Note: these function do not set values for aerosol_optics_device_data_.
  mam4::modal_aer_opt::set_complex_views_modal_aero(aerosol_optics_device_data_);
  mam4::modal_aer_opt::set_aerosol_optics_data_for_modal_aero_sw_views(aerosol_optics_device_data_);
  mam4::modal_aer_opt::set_aerosol_optics_data_for_modal_aero_lw_views(aerosol_optics_device_data_);

  mam_coupling::AerosolOpticsHostData aerosol_optics_host_data;

  std::map<std::string,FieldLayout>  layouts;
  std::map<std::string,view_1d_host> host_views;
  ekat::ParameterList rrtmg_params;

  mam_coupling::set_parameters_table(aerosol_optics_host_data,
                                     rrtmg_params,
                                     layouts,
                                     host_views );

  //FIXME: this name need to be pass in the input file.
  std::vector<std::string> name_table_modes=
  {
  "mam4_mode1_rrtmg_aeronetdust_c141106.nc",
  "mam4_mode2_rrtmg_c130628.nc",
  "mam4_mode3_rrtmg_aeronetdust_c141106.nc",
  "mam4_mode4_rrtmg_c130628.nc"
  };

  for (int imode = 0; imode < ntot_amode; imode++)
  {
    mam_coupling::read_rrtmg_table(name_table_modes[imode],
                                 imode,// mode No
                                 rrtmg_params, grid_,
                                 host_views,
                                 layouts,
                                 aerosol_optics_host_data,
                                 aerosol_optics_device_data_);
  }


// FIXME: we need to get this name from the yaml file.
std::string table_name_water = "water_refindex_rrtmg_c080910.nc";
// it will syn data to device.
mam_coupling::read_water_refindex(table_name_water, grid_,
                                  aerosol_optics_device_data_.crefwlw,
                                  aerosol_optics_device_data_.crefwsw);
//
{
 // make a list of host views
  std::map<std::string,view_1d_host> host_views_aero;
  // defines layouts
  std::map<std::string,FieldLayout>  layouts_aero;
  ekat::ParameterList params_aero;
  std::string surname_aero="aer";
  mam_coupling::set_refindex(surname_aero, params_aero, host_views_aero, layouts_aero );

  // read data
  /*
  soa:s-organic: /compyfs/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc
  dst:dust:      /compyfs/inputdata/atm/cam/physprops/dust_aeronet_rrtmg_c141106.nc
  ncl:seasalt:   /compyfs/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc
  so4:sulfate:   /compyfs/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc
  pom:p-organic: /compyfs/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709.nc
  bc :black-c:   /compyfs/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc
  mom:m-organic: /compyfs/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc */
  std::vector<std::string> name_table_aerosols=
  {
  "ocphi_rrtmg_c100508.nc", // soa:s-organic
  "dust_aeronet_rrtmg_c141106.nc",//dst:dust:
  "ssam_rrtmg_c100508.nc", // ncl:seasalt
  "sulfate_rrtmg_c080918.nc", // so4:sulfate
  "ocpho_rrtmg_c130709.nc", // pom:p-organic
  "bcpho_rrtmg_c100508.nc", // bc :black-c
  "poly_rrtmg_c130816.nc" // mom:m-organic
  };

  //FIXME: make a function that return a index given the species name
// specname_amode(ntot_aspectype) = (/ 'sulfate (0)   ',
//  'ammonium (1) ', 'nitrate (2)   ', &
      //  'p-organic (3) ', 's-organic (4) ', 'black-c (5)  ', &
      //  'seasalt (6)  ', 'dust  (7)    ', &
      //  'm-organic (8)' /)
  // FIXME: move this info to a conf file
  std::vector<int> species_ids=
  {
  4, // soa:s-organic
  7,//dst:dust:
  6, // ncl:seasalt
  0, // so4:sulfate
  3, // pom:p-organic
  5, // bc :black-c
  8 // mom:m-organic
  };


  constexpr int maxd_aspectype = mam4::ndrop::maxd_aspectype;
  auto specrefndxsw_host = mam_coupling::complex_view_2d::HostMirror(
        "specrefndxsw_host", nswbands_, maxd_aspectype);

  auto specrefndxlw_host = mam_coupling::complex_view_2d::HostMirror(
        "specrefndxlw_host", nlwbands_, maxd_aspectype);

  const int n_spec = size(name_table_aerosols);
  for (int ispec = 0; ispec < n_spec; ispec++)
  {
    // read data
    auto table_name = name_table_aerosols[ispec];
    // need to update table name
    params_aero.set("Filename",table_name);
    AtmosphereInput refindex_aerosol(params_aero, grid_,
    host_views_aero, layouts_aero);
    refindex_aerosol.read_variables();
    refindex_aerosol.finalize();

    // copy data to device
    int species_id =species_ids[ispec]; // FIXME: I need the species index for soa
    mam_coupling::set_refindex_aerosol(species_id,
                                      host_views_aero,
                                      specrefndxsw_host, // complex refractive index for water visible
                                      specrefndxlw_host);
  } // done ispec
  // reshape specrefndxsw_host and copy it to device
  mam4::modal_aer_opt::set_device_specrefindex(aerosol_optics_device_data_.specrefindex_sw, "short_wave", specrefndxsw_host);
  mam4::modal_aer_opt::set_device_specrefindex(aerosol_optics_device_data_.specrefindex_lw, "long_wave", specrefndxlw_host);
}

#if 0
{
  constexpr int coef_number = mam4::modal_aer_opt::coef_number;
    constexpr int refindex_real = mam4::modal_aer_opt::refindex_real;
  constexpr int refindex_im = mam4::modal_aer_opt::refindex_im;
auto host_1d = view_1d_host(refindex_real_lw_host.data(), refindex_real_lw_host.size());
scorpio::register_file(tables_filename_mode4,scorpio::Read);
std::string var_name = "refindex_real_lw";
std::vector<std::string> vec_of_dims;
vec_of_dims.push_back("lw_band");
vec_of_dims.push_back("refindex_real");
scorpio::register_dimension(tables_filename_mode4, vec_of_dims[0], vec_of_dims[0], nlwbands, 0);
scorpio::register_dimension(tables_filename_mode4, vec_of_dims[1], vec_of_dims[1], refindex_real, 0);
std::reverse(vec_of_dims.begin(),vec_of_dims.end());
const auto& fp_precision = "real";
std::string io_decomp_tag ="dt=real,grid-idx=0,layout=lw_band16-refindex_real7";
scorpio::register_variable(tables_filename_mode4, var_name, var_name,
                              vec_of_dims, fp_precision, io_decomp_tag);
std::vector<scorpio::offset_t> var_dof(host_1d.size());
std::iota(var_dof.begin(),var_dof.end(),0);
scorpio::set_dof(tables_filename_mode4,var_name,var_dof.size(),var_dof.data());
scorpio::set_decomp(tables_filename_mode4);
scorpio::grid_read_data_array(tables_filename_mode4,var_name,-1000,host_1d.data(), host_1d.size());
scorpio::eam_pio_closefile(tables_filename_mode4);}
#endif
}

}
void MAMOptics::run_impl(const double dt) {

  constexpr Real zero =0.0;

  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);

  // const auto policy = haero::ThreadTeamPolicy(1u, 1);
  /// outputs
  auto aero_g_sw   = get_field_out("aero_g_sw").get_view<Real***>();
  auto aero_ssa_sw = get_field_out("aero_ssa_sw").get_view<Real***>();
  auto aero_tau_sw = get_field_out("aero_tau_sw").get_view<Real***>();
  auto aero_tau_lw = get_field_out("aero_tau_lw").get_view<Real***>();

  auto aero_tau_forward = get_field_out("aero_tau_forward").get_view<Real***>();

  // NOTE: we do not compute this variable in aersol_optics
  auto aero_nccn   = get_field_out("nccn").get_view<Real**>(); // FIXME: get rid of this

  // they get modified in mam4xx
  state_q_ = get_field_out("state_q").get_view<Real***>();
  qqcw_ = get_field_out("qqcw").get_view<Real***>();

  const Real t = 0.0;

  if (false) { // remove when ready to do actual calculations
    // populate these fields with reasonable representative values
    // Kokkos::deep_copy(aero_g_sw, 0.5);
    // Kokkos::deep_copy(aero_ssa_sw, 0.7);
    // Kokkos::deep_copy(aero_tau_sw, 0.0);
    // Kokkos::deep_copy(aero_tau_lw, 0.0);
    // Kokkos::deep_copy(aero_nccn, 50.0);
  } else {

    // Compute optical properties on all local columns.
    // (Strictly speaking, we don't need this parallel_for here yet, but we leave
    //  it in anticipation of column-specific aerosol optics to come.)
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ThreadTeam& team) {
      const Int icol = team.league_rank(); // column index
      auto odap_aer_icol = ekat::subview(aero_tau_lw, icol);

      // FIXME: Get rid of this
      // auto nccn = ekat::subview(aero_nccn, icol);
      auto pmid =  ekat::subview(dry_atm_.p_mid, icol);
      auto temperature =  ekat::subview(dry_atm_.T_mid, icol);
      auto cldn = ekat::subview(dry_atm_.cldfrac, icol);

      // FIXME: interface pressure [Pa]
      auto pint =  ekat::subview(p_int_, icol);
      auto zm =  ekat::subview(z_mid_, icol);
      // FIXME: dry mass pressure interval [Pa]
      // FIXME:
      auto zi= ekat::subview(z_iface_, icol);
      auto pdel = ekat::subview(p_del_, icol);
      auto pdeldry = ekat::subview(dry_atm_.p_del, icol);
   auto state_q_icol = ekat::subview(state_q_, icol);
   auto qqcw_icol = ekat::subview(qqcw_, icol);

  auto ssa_cmip6_sw_icol = ekat::subview(ssa_cmip6_sw_, icol);
  auto af_cmip6_sw_icol = ekat::subview(af_cmip6_sw_, icol);
  auto ext_cmip6_sw_icol = ekat::subview(ext_cmip6_sw_, icol);
  auto ext_cmip6_lw_icol = ekat::subview(ext_cmip6_lw_, icol);

  // FIXME: check if this correct: Note that these variables have pver+1 levels
  // tau_w =>  aero_ssa_sw  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
   auto tau_w_icol = ekat::subview(aero_ssa_sw, icol);
  // tau_w_g => "aero_g_sw" (pcols,0:pver,nswbands) ! aerosol assymetry parameter * tau * w
  auto tau_w_g_icol = ekat::subview(aero_g_sw, icol);
  // tau_w_f(pcols,0:pver,nswbands) => aero_tau_forward  ? ! aerosol forward scattered fraction * tau * w
  auto tau_w_f_icol = ekat::subview(aero_tau_forward, icol);
  // tau  => aero_tau_sw (?)   (pcols,0:pver,nswbands) ! aerosol extinction optical depth
  auto tau_icol = ekat::subview(aero_tau_sw, icol);

  // work
  auto work_icol = ekat::subview(work_, icol);

    mam4::aer_rad_props::aer_rad_props_sw(team,  dt, zi, pmid,
    pint, temperature,
    zm, state_q_icol, qqcw_icol,
    pdel, pdeldry,
    cldn, ssa_cmip6_sw_icol,
    af_cmip6_sw_icol, ext_cmip6_sw_icol,
    tau_icol, tau_w_icol, tau_w_g_icol,
    tau_w_f_icol,
    aerosol_optics_device_data_, work_icol);

team.team_barrier();

mam4::aer_rad_props::aer_rad_props_lw(team,
    dt, pmid, pint,
    temperature, zm, zi,
    state_q_icol, qqcw_icol, pdel, pdeldry,
    cldn, ext_cmip6_lw_icol,
    aerosol_optics_device_data_,
    odap_aer_icol);


    });

  }
  printf("Done  with aerosol_optics \n");
}

void MAMOptics::finalize_impl()
{
}

} // namespace scream
