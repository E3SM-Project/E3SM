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

  grid_ = grids_manager->get_grid("Physics");
  const auto& grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();      // number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels(); // number of levels per column
  nswbands_ = mam4::modal_aer_opt::nswbands;//14;                           // number of shortwave bands
  nlwbands_ = mam4::modal_aer_opt::nlwbands;//16;                           // number of longwave bands

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Define aerosol optics fields computed by this process.
  auto nondim = Units::nondimensional();
  FieldLayout scalar3d_swband_layout { {COL, LEV, SWBND}, {ncol_, nlev_, nswbands_} };
  FieldLayout scalar3d_lwband_layout { {COL, LEV, LWBND}, {ncol_, nlev_, nlwbands_} };

  // layout for 3D (2d horiz X 1d vertical) variables
  FieldLayout scalar3d_layout_mid{ {COL, LEV}, {ncol_, nlev_} };
  // FIXME: I switch the order of dimension.
  FieldLayout scalar3d_layout_int{ {COL, ILEV}, {ncol_, nlev_+1} };

  add_field<Required>("T_mid", scalar3d_layout_mid, K, grid_name); // Temperature
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name); // total pressure
  add_field<Required>("cldfrac_tot", scalar3d_layout_mid, nondim, grid_name); // total pressure
  add_field<Required>("z_int", scalar3d_layout_int, m, grid_name); // vertical position at interface
  add_field<Required>("z_mid", scalar3d_layout_mid, m, grid_name); // vertical position pressure
  add_field<Required>("p_int", scalar3d_layout_int, Pa, grid_name); // total pressure
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa, grid_name);
  add_field<Required>("pseudo_density_dry", scalar3d_layout_mid, Pa, grid_name);

#if 1
  // shortwave aerosol scattering asymmetry parameter [-]
  add_field<Computed>("aero_g_sw",   scalar3d_swband_layout, nondim, grid_name);
  // shortwave aerosol single-scattering albedo [-]
  add_field<Computed>("aero_ssa_sw", scalar3d_swband_layout, nondim, grid_name);
  // shortwave aerosol optical depth [-]
  add_field<Computed>("aero_tau_sw", scalar3d_swband_layout, nondim, grid_name);
  // longwave aerosol optical depth [-]
  add_field<Computed>("aero_tau_lw", scalar3d_lwband_layout, nondim, grid_name);

  // // aerosol extinction optical depth
  add_field<Computed>("aero_tau_forward", scalar3d_swband_layout, nondim, grid_name);

    // FIXME: this field doesn't belong here, but this is a convenient place to
  // FIXME: put it for now.
  // number mixing ratio for CCN
  using Spack      = ekat::Pack<Real,SCREAM_SMALL_PACK_SIZE>;
  using Pack       = ekat::Pack<Real,Spack::n>;
  constexpr int ps = Pack::n;
  // FieldLayout scalar3d_layout_mid { {COL, LEV}, {ncol_, nlev_} };
  add_field<Computed>("nccn", scalar3d_layout_mid, 1/kg, grid_name, ps);

#endif



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

  // FieldLayout scalar3d_refrtablw_layout { {NMODES,LWBND, LEV},
  //                                    {mam4::AeroConfig::num_modes,
  //                                     mam4::modal_aer_opt::nlwbands,
  //                                     mam4::modal_aer_opt::refindex_real} };



  // Construct the grid needed for input:
  // auto grid = std::make_shared<PointGrid>("grid",num_local_cols,num_global_cols,source_data_nlevs,comm);
  // Kokkos::deep_copy(grid->get_dofs_gids().template get_view<gid_type*>(),unique_src_dofs);
  // grid->get_dofs_gids().sync_to_host();






}

void MAMOptics::initialize_impl(const RunType run_type) {

#if 1
  dry_atm_.T_mid     = get_field_in("T_mid").get_view<const Real**>();
  dry_atm_.p_mid     = get_field_in("p_mid").get_view<const Real**>();
  // FIXME, there are two version of p_int in the nc file: p_dry_int and p_int
  // change to const Real
  p_int_     = get_field_in("p_int").get_view<const Real**>();

  dry_atm_.cldfrac   = get_field_in("cldfrac_tot").get_view<const Real**>(); // FIXME: tot or liq?
  // dry_atm_.pblh      = get_field_in("pbl_height").get_view<const Real*>();
  // FIXME: use const Real; why are using buffer in microphysics
  z_mid_     = get_field_in("z_mid").get_view<const Real**>();
  z_iface_   = get_field_in("z_int").get_view<const Real**>();

  p_del_     = get_field_in("pseudo_density").get_view<const Real**>();
  // FIXME: In the nc file, there is also pseudo_density_dry
  dry_atm_.p_del     = get_field_in("pseudo_density_dry").get_view<const Real**>();

  // FIXME: we have nvars in several process.
  constexpr int nvars = mam4::ndrop::nvars;
  constexpr int nlwbands = mam4::modal_aer_opt::nlwbands;
  constexpr int nswbands = mam4::modal_aer_opt::nswbands;
  constexpr int maxd_aspectype = mam4::ndrop::maxd_aspectype;
  constexpr int ntot_amode = mam4::AeroConfig::num_modes();

  state_q_ = mam_coupling::view_3d("state_q_", ncol_, nlev_, nvars);
  Kokkos::deep_copy(state_q_,10);
  qqcw_ = mam_coupling::view_3d("qqcw_", ncol_, nlev_, nvars);
  Kokkos::deep_copy(qqcw_,10);

  // FIXME: work arrays
  mass_ = mam_coupling::view_2d("mass", ncol_,  nlev_);
  cheb_ = mam_coupling::view_3d("cheb", ncol_, mam4::modal_aer_opt::ncoef, nlev_);

  dgnumwet_m_ = mam_coupling::view_3d("dgnumwet_m", ncol_, nlev_, ntot_amode);
  dgnumdry_m_ = mam_coupling::view_3d("dgnumdry_m", ncol_, nlev_, ntot_amode);

  radsurf_ = mam_coupling::view_2d("radsurf",ncol_, nlev_);
  logradsurf_ = mam_coupling::view_2d("logradsurf",ncol_,nlev_);

  air_density_ = mam_coupling::view_2d("air_density", ncol_,  nlev_);
  ext_cmip6_sw_inv_m_ = mam_coupling::view_3d ("ext_cmip6_sw_inv_m", ncol_, nswbands, nlev_);


  const int nwbands = nlwbands > nswbands ? nlwbands: nswbands;

  specrefindex_=mam_coupling::complex_view_3d("specrefindex", ncol_,
  mam4::modal_aer_opt::max_nspec, nwbands);
  qaerwat_m_ = mam_coupling::view_3d ("qaerwat_m", ncol_, nlev_, ntot_amode);
  ext_cmip6_lw_inv_m_ = mam_coupling::view_3d ("ext_cmip6_lw_inv_m", ncol_, nlev_, nlwbands);


  // aer_rad_props_sw inputs that are prescribed, i.e., we need a netcdf file.
  ssa_cmip6_sw_ = mam_coupling::view_3d ("ssa_cmip6_sw", ncol_, nlev_, nswbands);
  af_cmip6_sw_ = mam_coupling::view_3d ("af_cmip6_sw", ncol_, nlev_, nswbands);
  ext_cmip6_sw_= mam_coupling::view_3d ("ext_cmip6_sw", ncol_, nlev_, nswbands);
  Kokkos::deep_copy(ssa_cmip6_sw_, 1.0);
  Kokkos::deep_copy(af_cmip6_sw_, 1.0);
  Kokkos::deep_copy(ext_cmip6_sw_, 1.0);

  ext_cmip6_lw_ = mam_coupling::view_2d("ext_cmip6_lw_", nlev_, nlwbands);
  // odap_aer_ = mam_coupling::view_2d("odap_aer_", nlev_, nlwbands);

  Kokkos::deep_copy(ext_cmip6_lw_, 1.0);
  // Kokkos::deep_copy(odap_aer_, 1.0);


#endif
 // read table info
#if 1
{
  using namespace ShortFieldTagsNames;

  constexpr int refindex_real = mam4::modal_aer_opt::refindex_real;
  constexpr int refindex_im = mam4::modal_aer_opt::refindex_im;
  constexpr int coef_number = mam4::modal_aer_opt::coef_number;

  using view_1d_host = typename KT::view_1d<Real>::HostMirror;
  // Set up input structure to read data from file.
  using strvec_t = std::vector<std::string>;

  mam_coupling::AerosolOpticsDeviceData aerosol_optics_device_data;
// FIXME: move to a function
for (int d1 = 0; d1 < ntot_amode; ++d1)
{
  for (int d5 = 0; d5 < nlwbands; ++d5) {
        absplw_[d1][d5] =
            mam_coupling::view_3d("absplw_", coef_number,
                                             refindex_real,
                                             refindex_im);
        // FIXME
        aerosol_optics_device_data.absplw[d1][d5]=absplw_[d1][d5];
  }

   for (int d3 = 0; d3 < nlwbands; ++d3) {
        refrtablw_[d1][d3] = mam_coupling::view_1d("refrtablw",
       refindex_real);
        refitablw_[d1][d3] = mam_coupling::view_1d("refitablw",
       refindex_im);
         // FIXME
         aerosol_optics_device_data.refrtablw[d1][d3]=refrtablw_[d1][d3];
         aerosol_optics_device_data.refitablw[d1][d3]=refitablw_[d1][d3];
  } // d3

  for (int d3 = 0; d3 < nswbands; ++d3) {
    refrtabsw_[d1][d3] = mam_coupling::view_1d ("refrtabsw", refindex_real);
    refitabsw_[d1][d3] = mam_coupling::view_1d("refitabsw", refindex_im);
    //FIXME
    aerosol_optics_device_data.refrtabsw[d1][d3]=refrtabsw_[d1][d3];
    aerosol_optics_device_data.refitabsw[d1][d3]=refitabsw_[d1][d3];
    } // d3

      // allocate memory for tables: abspsw,extpsw,asmpsw

     for (int d5 = 0; d5 < nswbands; ++d5) {
        abspsw_[d1][d5] =
            mam_coupling::view_3d("abspsw", coef_number, refindex_real, refindex_im);
        extpsw_[d1][d5] =
            mam_coupling::view_3d("extpsw", coef_number, refindex_real, refindex_im);
        asmpsw_[d1][d5] =
            mam_coupling::view_3d("asmpsw", coef_number, refindex_real, refindex_im);
     //FIXME
    aerosol_optics_device_data.abspsw[d1][d5]=abspsw_[d1][d5];
    aerosol_optics_device_data.extpsw[d1][d5]=extpsw_[d1][d5];
    aerosol_optics_device_data.asmpsw[d1][d5]=asmpsw_[d1][d5];
     } // d5
  }

#if 1

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
                                 aerosol_optics_device_data);
  }

#endif

// FIXME: we need to get this name from the yaml file.
std::string table_name = "water_refindex_rrtmg_c080910.nc";
mam_coupling::read_water_refindex(table_name, grid_, crefwlw_,crefwsw_);

//
{
   // FIXME: I need to get these variables from a netCDF file.
  specrefndxsw_ = mam_coupling::complex_view_2d("specrefndxlw_",nswbands, maxd_aspectype );
  // Kokkos::deep_copy(specrefndxsw_, 1.0);
  specrefndxlw_  = mam_coupling::complex_view_2d("specrefndxlw_",nlwbands, maxd_aspectype );
  // Kokkos::deep_copy(specrefndxlw_, 1.0);

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
    mam_coupling::copy_refindex_aerosol_to_device(species_id,
                                             host_views_aero,
                                             specrefndxsw_, // complex refractive index for water visible
                                             specrefndxlw_);

  } // done ispec

}



#if 0
{
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
#endif
}
void MAMOptics::run_impl(const double dt) {

  constexpr Real zero =0.0;
#if 1
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);
#if 1
  // get the aerosol optics fields
  auto aero_g_sw   = get_field_out("aero_g_sw").get_view<Real***>();
  auto aero_ssa_sw = get_field_out("aero_ssa_sw").get_view<Real***>();
  auto aero_tau_sw = get_field_out("aero_tau_sw").get_view<Real***>();
  auto aero_tau_lw = get_field_out("aero_tau_lw").get_view<Real***>();

  auto aero_tau_forward = get_field_out("aero_tau_forward").get_view<Real***>();
  // diagnostics

  // aerosol extinction [1/m]
  auto extinct = get_field_out("extinct").get_view<Real**>();
  // aerosol absorption [1/m]
  auto absorb = get_field_out("absorb").get_view<Real**>();

  auto aodnir = get_field_out("aodnir").get_view<Real*>();
  auto aoduv = get_field_out("aoduv").get_view<Real*>();
  // dustaodmode[ntot_amode],
  auto dustaodmode = get_field_out("dustaodmode").get_view<Real**>();
  // aodmode[ntot_amode]
  auto aodmode = get_field_out("aodmode").get_view<Real**>();
  // burdenmode[ntot_amode]
  auto burdenmode = get_field_out("burdenmode").get_view<Real**>();

  auto aodabsbc = get_field_out("aodabsbc").get_view<Real*>();
  auto aodvis = get_field_out("aodvis").get_view<Real*>();
  auto aodall = get_field_out("aodall").get_view<Real*>();
  auto  ssavis = get_field_out("ssavis").get_view<Real*>();
  auto  aodabs = get_field_out("aodabs").get_view<Real*>();
  auto  burdendust= get_field_out("burdendust").get_view<Real*>();

  auto  burdenso4= get_field_out("burdenso4").get_view<Real*>();
  auto  burdenbc= get_field_out("burdenbc").get_view<Real*>();
  auto  burdenpom= get_field_out("burdenpom").get_view<Real*>();
  auto  burdensoa= get_field_out("burdensoa").get_view<Real*>();

  auto  burdenseasalt= get_field_out("burdenseasalt").get_view<Real*>();
  auto  burdenmom= get_field_out("burdenmom").get_view<Real*>();
  auto  momaod= get_field_out("momaod").get_view<Real*>();
  auto  dustaod= get_field_out("dustaod").get_view<Real*>();

  auto  so4aod= get_field_out("so4aod").get_view<Real*>();
  auto  pomaod= get_field_out("pomaod").get_view<Real*>();
  auto  soaaod= get_field_out("soaaod").get_view<Real*>();
  auto  bcaod= get_field_out("bcaod").get_view<Real*>();

  auto  seasaltaod= get_field_out("seasaltaod").get_view<Real*>();

  Kokkos::deep_copy(bcaod,zero);



  printf("dt %e\n",dt);

  auto aero_nccn   = get_field_out("nccn").get_view<Real**>(); // FIXME: get rid of this
#endif
  // constexpr int pver = mam4::nlev;
  constexpr int ntot_amode=mam4::AeroConfig::num_modes();
  constexpr int maxd_aspectype = mam4::ndrop::maxd_aspectype;
  constexpr int nspec_max = mam4::ndrop::nspec_max;
  constexpr int num_aerosol_ids = mam4::AeroConfig::num_aerosol_ids();

  const Real t = 0.0;


      mam4::AeroId specname_amode[9] = {mam4::AeroId::SO4,  // sulfate
                                      mam4::AeroId::None, // ammonium
                                      mam4::AeroId::None, // nitrate
                                      mam4::AeroId::POM,  // p-organic
                                      mam4::AeroId::SOA,  // s-organic
                                      mam4::AeroId::BC,   // black-c
                                      mam4::AeroId::NaCl, // seasalt
                                      mam4::AeroId::DST,  // dust
                                      mam4::AeroId::MOM}; // m-organic


  if (false) { // remove when ready to do actual calculations
    // populate these fields with reasonable representative values
    // Kokkos::deep_copy(aero_g_sw, 0.5);
    // Kokkos::deep_copy(aero_ssa_sw, 0.7);
    // Kokkos::deep_copy(aero_tau_sw, 0.0);
    // Kokkos::deep_copy(aero_tau_lw, 0.0);
    // Kokkos::deep_copy(aero_nccn, 50.0);
  } else {


    const bool do_adjust = true;
    const bool do_aitacc_transfer = true;
    const bool update_mmr = false;

    // Compute optical properties on all local columns.
    // (Strictly speaking, we don't need this parallel_for here yet, but we leave
    //  it in anticipation of column-specific aerosol optics to come.)
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ThreadTeam& team) {
      const Int icol = team.league_rank(); // column index
      // auto g_sw = ekat::subview(aero_g_sw, icol);
      // auto ssa_sw = ekat::subview(aero_ssa_sw, icol);
      // auto odap_aer_icol = ekat::subview(aero_tau_sw, icol);
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
 #if 0
      printf("temperature %e\n",temperature(0));
      printf("pmid %e\n",pmid(0));
      printf("cldn %e\n",cldn(0));
      printf("state_q_ %e\n",state_q_(0,0));
      printf("zm %e\n",zm(0));
      printf("pdeldry %e\n",pdeldry(0));
      printf("pdel %e\n",pdel(0));
#endif
      int nspec_amode[ntot_amode];

          int lspectype_amode[maxd_aspectype][ntot_amode];
          int lmassptr_amode[maxd_aspectype][ntot_amode];
          Real specdens_amode[maxd_aspectype];
          Real spechygro[maxd_aspectype];
          int numptr_amode[ntot_amode];
          int mam_idx[ntot_amode][nspec_max];
          int mam_cnst_idx[ntot_amode][nspec_max];

          mam4::ndrop::get_e3sm_parameters(nspec_amode, lspectype_amode, lmassptr_amode,
                              numptr_amode, specdens_amode, spechygro, mam_idx,
                              mam_cnst_idx);

          Real sigmag_amode[ntot_amode] = {0.18000000000000000e+001,
                                           0.16000000000000001e+001,
                                           0.18000000000000000e+001,
                                           0.16000000238418579e+001};

#if 1
          Real inv_density[ntot_amode]
                          [num_aerosol_ids] = {};
          Real num2vol_ratio_min[ntot_amode] = {};
          Real num2vol_ratio_max[ntot_amode] = {};
          Real num2vol_ratio_max_nmodes[ntot_amode] = {};
          Real num2vol_ratio_min_nmodes[ntot_amode] = {};
          Real num2vol_ratio_nom_nmodes[ntot_amode] = {};
          Real dgnmin_nmodes[ntot_amode] = {};
          Real dgnmax_nmodes[ntot_amode] = {};
          Real dgnnom_nmodes[ntot_amode] = {};
          Real mean_std_dev_nmodes[ntot_amode] = {};
          // outputs
          bool noxf_acc2ait[num_aerosol_ids] = {};
          int n_common_species_ait_accum = {};
          int ait_spec_in_acc[num_aerosol_ids] = {};
          int acc_spec_in_ait[num_aerosol_ids] = {};

           mam4::modal_aero_calcsize::init_calcsize(
              inv_density, num2vol_ratio_min, num2vol_ratio_max,
              num2vol_ratio_max_nmodes, num2vol_ratio_min_nmodes,
              num2vol_ratio_nom_nmodes, dgnmin_nmodes, dgnmax_nmodes,
              dgnnom_nmodes, mean_std_dev_nmodes,
              // outputs
              noxf_acc2ait, n_common_species_ait_accum, ait_spec_in_acc,
              acc_spec_in_ait);

          team.team_barrier();
          {

          for (int kk = mam4::ndrop::top_lev; kk < nlev_; ++kk) {

            const auto state_q_k = Kokkos::subview(state_q_, icol, kk, Kokkos::ALL());
            const auto qqcw_k = Kokkos::subview(qqcw_, icol, kk, Kokkos::ALL());
            auto dgncur_i = Kokkos::subview(dgnumdry_m_, icol, kk, Kokkos::ALL());
            Real dgncur_c[ntot_amode] = {};
            mam4::modal_aero_calcsize::modal_aero_calcsize_sub(
                state_q_k.data(), // in
                qqcw_k.data(),    // in
                dt, do_adjust, do_aitacc_transfer, update_mmr, lmassptr_amode,
                numptr_amode,
                inv_density, // in
                num2vol_ratio_min, num2vol_ratio_max, num2vol_ratio_max_nmodes,
                num2vol_ratio_min_nmodes, num2vol_ratio_nom_nmodes,
                dgnmin_nmodes, dgnmax_nmodes, dgnnom_nmodes,
                mean_std_dev_nmodes,
                // outputs
                noxf_acc2ait, n_common_species_ait_accum, ait_spec_in_acc,
                acc_spec_in_ait, dgncur_i.data(), dgncur_c);
          } // k
          }
          team.team_barrier();

#endif


  #if 1

   auto mass_icol = ekat::subview(mass_, icol);
   auto radsurf_icol = ekat::subview(radsurf_, icol);
   auto logradsurf_icol = ekat::subview(logradsurf_, icol);

   auto cheb_icol = ekat::subview(cheb_, icol);
   auto dgnumwet_m_icol = ekat::subview(dgnumwet_m_, icol);
   auto dgnumdry_m_icol = ekat::subview(dgnumdry_m_, icol);
   auto specrefindex_icol = ekat::subview(specrefindex_, icol);
   auto qaerwat_m_icol = ekat::subview(qaerwat_m_, icol);
   auto ext_cmip6_lw_inv_m_icol = ekat::subview(ext_cmip6_lw_inv_m_, icol);
   auto state_q_icol = ekat::subview(state_q_, icol);

#endif

#if 1

  auto ssa_cmip6_sw_icol = ekat::subview(ssa_cmip6_sw_, icol);
  auto af_cmip6_sw_icol = ekat::subview(af_cmip6_sw_, icol);
  auto ext_cmip6_sw_icol = ekat::subview(ext_cmip6_sw_, icol);

  // FIXME: check if this correct: Note that these variables have pver+1 levels
  // tau_w =>  aero_ssa_sw  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
   auto tau_w_icol = ekat::subview(aero_ssa_sw, icol);
  // tau_w_g => "aero_g_sw" (pcols,0:pver,nswbands) ! aerosol assymetry parameter * tau * w
  auto tau_w_g_icol = ekat::subview(aero_g_sw, icol);
  // tau_w_f(pcols,0:pver,nswbands) => aero_tau_forward  ? ! aerosol forward scattered fraction * tau * w
  auto tau_w_f_icol = ekat::subview(aero_tau_forward, icol);
  // tau  => aero_tau_sw (?)   (pcols,0:pver,nswbands) ! aerosol extinction optical depth
  auto tau_icol = ekat::subview(aero_tau_sw, icol);

  auto air_density_icol = ekat::subview(air_density_, icol);

  auto ext_cmip6_sw_inv_m_icol = ekat::subview(ext_cmip6_sw_inv_m_, icol);

  auto extinct_icol = ekat::subview(extinct, icol);
  // aerosol absorption [1/m]
  auto absorb_icol = ekat::subview(absorb, icol);
  // dustaodmode[ntot_amode],
  auto dustaodmode_icol = ekat::subview(dustaodmode, icol);
  // aodmode[ntot_amode]
  auto aodmode_icol = ekat::subview(aodmode, icol);
  // burdenmode[ntot_amode]
  auto burdenmode_icol = ekat::subview(burdenmode, icol);

#if 1
    // FIXME: following variables, wihout 2, are producing nans.
    Real aodabsbc2 = zero;
    Real bcaod2 = zero;
    Real dustaod2 = zero;
    Real momaod2 = zero;
    Real pomaod2 = zero;
    Real seasaltaod2 =zero;
    Real so4aod2 = zero;
    Real soaaod2 = zero;

    team.team_barrier();

     mam4::aer_rad_props::aer_rad_props_sw( dt, zi, pmid,
    pint, temperature,
    zm, state_q_icol,
    pdel, pdeldry,
    cldn, ssa_cmip6_sw_icol,
    af_cmip6_sw_icol, ext_cmip6_sw_icol,
    // nnite, idxnite,
    // is_cmip6_volc,
    // const ColumnView qqcw_fld[pcnst],
    tau_icol, tau_w_icol, tau_w_g_icol,
    tau_w_f_icol,
    nspec_amode,
    sigmag_amode,
    lmassptr_amode,
    spechygro, specdens_amode,
    lspectype_amode,
    specrefndxsw_, // specrefndxsw( nswbands, maxd_aspectype )
    crefwlw_,
    crefwsw_,
    // FIXME
    specname_amode,
    extpsw_,
    abspsw_,
    asmpsw_,
    refrtabsw_,
    refitabsw_,
    // diagnostic
    extinct_icol, //        ! aerosol extinction [1/m]
    absorb_icol,  //         ! aerosol absorption [1/m]
    aodnir[icol], aoduv[icol], dustaodmode_icol.data(),
    aodmode_icol.data(), burdenmode_icol.data(), aodabsbc2,
    aodvis[icol], aodall[icol], ssavis[icol], aodabs[icol], burdendust[icol],
    burdenso4[icol], burdenbc[icol], burdenpom[icol], burdensoa[icol],
    burdenseasalt[icol], burdenmom[icol], momaod2, dustaod2,
    so4aod2, // total species AOD
    pomaod2, soaaod2, bcaod2, seasaltaod2,
    // work views
    mass_icol, air_density_icol, cheb_icol,
    dgnumwet_m_icol, dgnumdry_m_icol,
    radsurf_icol, logradsurf_icol,
    specrefindex_icol, qaerwat_m_icol,
    ext_cmip6_sw_inv_m_icol);

    // printf("aodabsbc2 %e \n",aodabsbc2);
#endif
#endif


          team.team_barrier();
          {
          for (int kk = mam4::ndrop::top_lev; kk < nlev_; ++kk) {

            const auto state_q_k = Kokkos::subview(state_q_, icol, kk, Kokkos::ALL());
            const auto qqcw_k = Kokkos::subview(qqcw_, icol, kk, Kokkos::ALL());
            auto dgncur_i = Kokkos::subview(dgnumdry_m_, icol, kk, Kokkos::ALL());
            Real dgncur_c[ntot_amode] = {};
            mam4::modal_aero_calcsize::modal_aero_calcsize_sub(
                state_q_k.data(), // in
                qqcw_k.data(),    // in
                dt, do_adjust, do_aitacc_transfer, update_mmr, lmassptr_amode,
                numptr_amode,
                inv_density, // in
                num2vol_ratio_min, num2vol_ratio_max, num2vol_ratio_max_nmodes,
                num2vol_ratio_min_nmodes, num2vol_ratio_nom_nmodes,
                dgnmin_nmodes, dgnmax_nmodes, dgnnom_nmodes,
                mean_std_dev_nmodes,
                // outputs
                noxf_acc2ait, n_common_species_ait_accum, ait_spec_in_acc,
                acc_spec_in_ait, dgncur_i.data(), dgncur_c);
          } // k
          }

mam4::aer_rad_props::aer_rad_props_lw(
    dt, pmid, pint,
    temperature, zm, zi,
    state_q_icol, pdel, pdeldry,
    cldn, ext_cmip6_lw_,
    // const ColumnView qqcw_fld[pcnst],
    odap_aer_icol,
    //
    nspec_amode, sigmag_amode,
    lmassptr_amode,
    spechygro, specdens_amode,
    lspectype_amode,
    specrefndxlw_,
    crefwlw_,
    crefwsw_,
    absplw_,
    refrtablw_,
    refitablw_,
    // work views
    mass_icol, cheb_icol, dgnumwet_m_icol,
    dgnumdry_m_icol, radsurf_icol,
    logradsurf_icol, specrefindex_icol,
    qaerwat_m_icol, ext_cmip6_lw_inv_m_icol);


    });

// #if 0
//     Kokkos::deep_copy(aero_g_sw, 0.5);
//     Kokkos::deep_copy(aero_ssa_sw, 0.7);
//     Kokkos::deep_copy(aero_tau_sw, 0.0);
//     Kokkos::deep_copy(aero_tau_lw, 0.0);
//     Kokkos::deep_copy(aero_nccn, 50.0);
// #endif

  }
  printf("Done  with aerosol_optics \n");
#endif
}

void MAMOptics::finalize_impl()
{
}

} // namespace scream
