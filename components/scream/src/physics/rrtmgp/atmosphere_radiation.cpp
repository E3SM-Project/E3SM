#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "physics/rrtmgp/atmosphere_radiation.hpp"
#include "physics/rrtmgp/rrtmgp_utils.hpp"
#include "physics/rrtmgp/share/shr_orb_mod_c2f.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "share/util/scream_column_ops.hpp"
#include "cpp/rrtmgp/mo_gas_concentrations.h"
#include "YAKL.h"
#include "ekat/ekat_assert.hpp"

namespace scream {

  using KT = KokkosTypes<DefaultDevice>;
  using ExeSpace = KT::ExeSpace;
  using MemberType = KT::MemberType;

RRTMGPRadiation::
RRTMGPRadiation (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}  // RRTMGPRadiation::RRTMGPRadiation

// =========================================================================================
void RRTMGPRadiation::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {

  using namespace ekat::units;

  // Gather the active gases from the rrtmgp parameter list and assign to the m_gas_names vector.
  auto active_gases = m_params.get<std::vector<std::string>>("active_gases");
  for (auto& it : active_gases) {
    // Make sure only unique names are added
    if (std::find(m_gas_names.begin(), m_gas_names.end(), it) == m_gas_names.end()) {
      m_gas_names.push_back(it);
    }
  }
  m_ngas = m_gas_names.size();

  // Declare the set of fields used by rrtmgp
  auto kgkg = kg/kg;
  kgkg.set_string("kg/kg");
  auto m3 = m * m * m;
  auto Wm2 = W / m / m;
  Wm2.set_string("W/m2");
  auto nondim = m/m;  // dummy unit for non-dimensional fields
  auto micron = m / 1000000;

  using namespace ShortFieldTagsNames;

  const auto& grid_name = m_params.get<std::string>("Grid");
  m_grid = grids_manager->get_grid(grid_name);
  m_ncol = m_grid->get_num_local_dofs();
  m_nlay = m_grid->get_num_vertical_levels();
  m_lat  = m_grid->get_geometry_data("lat");
  m_lon  = m_grid->get_geometry_data("lon");

  // Set up dimension layouts
  FieldLayout scalar2d_layout     { {COL   }, {m_ncol    } };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_ncol,m_nlay} };
  FieldLayout scalar3d_layout_int { {COL,ILEV}, {m_ncol,m_nlay+1} };
  FieldLayout scalar3d_swband_layout { {COL,SWBND,LEV}, {m_ncol, m_nswbands, m_nlay} };
  FieldLayout scalar3d_lwband_layout { {COL,LWBND,LEV}, {m_ncol, m_nlwbands, m_nlay} };

  constexpr int ps = SCREAM_SMALL_PACK_SIZE;
  // Set required (input) fields here
  add_field<Required>("p_mid" , scalar3d_layout_mid, Pa, grid_name, ps);
  add_field<Required>("p_int", scalar3d_layout_int, Pa, grid_name, ps);
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa, grid_name, ps);
  // WARNING! We are switch surface albedo variables to "Updated" to allow rrtmgp
  // property checks to repair them as a precondition check.  TODO: Change back
  // to "Required" when surface coupling is it's own process and can handle repair
  // locally.
  add_field<Updated>("sfc_alb_dir_vis", scalar2d_layout, nondim, grid_name);
  add_field<Updated>("sfc_alb_dir_nir", scalar2d_layout, nondim, grid_name);
  add_field<Updated>("sfc_alb_dif_vis", scalar2d_layout, nondim, grid_name);
  add_field<Updated>("sfc_alb_dif_nir", scalar2d_layout, nondim, grid_name);
  // End WARNING Message
  add_field<Required>("qc", scalar3d_layout_mid, kgkg, grid_name, ps);
  add_field<Required>("qi", scalar3d_layout_mid, kgkg, grid_name, ps);
  add_field<Required>("cldfrac_tot", scalar3d_layout_mid, nondim, grid_name, ps);
  add_field<Required>("eff_radius_qc", scalar3d_layout_mid, micron, grid_name, ps);
  add_field<Required>("eff_radius_qi", scalar3d_layout_mid, micron, grid_name, ps);
  add_field<Required>("qv",scalar3d_layout_mid,kgkg,grid_name, ps);
  add_field<Required>("surf_lw_flux_up",scalar2d_layout,W/(m*m),grid_name);
  // Set of required gas concentration fields
  for (auto& it : m_gas_names) {
    if (it == "h2o") { /* Special case where water vapor is called h2o in radiation */
      // do nothing, qv has already been added.
    } else {
      add_field<Required>(it,scalar3d_layout_mid,kgkg,grid_name, ps);
    }
  }
  // Required aerosol optical properties from SPA
  m_do_aerosol_rad = m_params.get<bool>("do_aerosol_rad",true);
  if (m_do_aerosol_rad) {
    add_field<Required>("aero_tau_sw", scalar3d_swband_layout, nondim, grid_name, ps);
    add_field<Required>("aero_ssa_sw", scalar3d_swband_layout, nondim, grid_name, ps);
    add_field<Required>("aero_g_sw"  , scalar3d_swband_layout, nondim, grid_name, ps);
    add_field<Required>("aero_tau_lw", scalar3d_lwband_layout, nondim, grid_name, ps);
  }

  // Set computed (output) fields
  add_field<Updated >("T_mid"     , scalar3d_layout_mid, K  , grid_name, ps);
  add_field<Computed>("SW_flux_dn", scalar3d_layout_int, Wm2, grid_name, ps);
  add_field<Computed>("SW_flux_up", scalar3d_layout_int, Wm2, grid_name, ps);
  add_field<Computed>("SW_flux_dn_dir", scalar3d_layout_int, Wm2, grid_name, ps);
  add_field<Computed>("LW_flux_up", scalar3d_layout_int, Wm2, grid_name, ps);
  add_field<Computed>("LW_flux_dn", scalar3d_layout_int, Wm2, grid_name, ps);
  add_field<Computed>("rad_heating_pdel", scalar3d_layout_mid, Pa*K/s, grid_name, ps);

  // Translation of variables from EAM
  // --------------------------------------------------------------
  // EAM name | EAMXX name       | Description
  // --------------------------------------------------------------
  // soll       sfc_flux_dir_nir   solar near-IR direct flux
  // sols       sfc_flux_dir_vis   solar UV/visible direct flux
  // solld      sfc_flux_dif_nir   solar near-ID diffuse flux
  // solsd      sfc_flux_dif_vis   solar UV/visible diffuse flux
  // netsw      sfc_flux_sw_net    net (down - up) SW flux at surface
  // flwds      sfc_flux_lw_dn     downwelling LW flux at surface
  // --------------------------------------------------------------
  add_field<Computed>("sfc_flux_dir_nir", scalar2d_layout, Wm2, grid_name);
  add_field<Computed>("sfc_flux_dir_vis", scalar2d_layout, Wm2, grid_name);
  add_field<Computed>("sfc_flux_dif_nir", scalar2d_layout, Wm2, grid_name);
  add_field<Computed>("sfc_flux_dif_vis", scalar2d_layout, Wm2, grid_name);
  add_field<Computed>("sfc_flux_sw_net" , scalar2d_layout, Wm2, grid_name);
  add_field<Computed>("sfc_flux_lw_dn"  , scalar2d_layout, Wm2, grid_name);
}  // RRTMGPRadiation::set_grids

size_t RRTMGPRadiation::requested_buffer_size_in_bytes() const
{
  const size_t interface_request = Buffer::num_1d_ncol*m_ncol*sizeof(Real) +
                                   Buffer::num_2d_nlay*m_ncol*m_nlay*sizeof(Real) +
                                   Buffer::num_2d_nlay_p1*m_ncol*(m_nlay+1)*sizeof(Real) +
                                   Buffer::num_2d_nswbands*m_ncol*m_nswbands*sizeof(Real) +
                                   Buffer::num_3d_nlev_nswbands*m_ncol*(m_nlay+1)*m_nswbands*sizeof(Real) +
                                   Buffer::num_3d_nlev_nlwbands*m_ncol*(m_nlay+1)*m_nlwbands*sizeof(Real) +
                                   Buffer::num_3d_nlay_nswbands*m_ncol*(m_nlay)*m_nswbands*sizeof(Real) +
                                   Buffer::num_3d_nlay_nlwbands*m_ncol*(m_nlay)*m_nlwbands*sizeof(Real);

  return interface_request;
} // RRTMGPRadiation::requested_buffer_size
// =========================================================================================

void RRTMGPRadiation::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d array
  m_buffer.mu0 = decltype(m_buffer.mu0)("mu0", mem, m_ncol);
  mem += m_buffer.mu0.totElems();
  m_buffer.sfc_alb_dir_vis = decltype(m_buffer.sfc_alb_dir_vis)("sfc_alb_dir_vis", mem, m_ncol);
  mem += m_buffer.sfc_alb_dir_vis.totElems();
  m_buffer.sfc_alb_dir_nir = decltype(m_buffer.sfc_alb_dir_nir)("sfc_alb_dir_nir", mem, m_ncol);
  mem += m_buffer.sfc_alb_dir_nir.totElems();
  m_buffer.sfc_alb_dif_vis = decltype(m_buffer.sfc_alb_dif_vis)("sfc_alb_dif_vis", mem, m_ncol);
  mem += m_buffer.sfc_alb_dif_vis.totElems();
  m_buffer.sfc_alb_dif_nir = decltype(m_buffer.sfc_alb_dif_nir)("sfc_alb_dif_nir", mem, m_ncol);
  mem += m_buffer.sfc_alb_dif_nir.totElems();
  m_buffer.cosine_zenith = decltype(m_buffer.cosine_zenith)(mem, m_ncol);
  mem += m_buffer.cosine_zenith.size();
  m_buffer.sfc_flux_dir_vis = decltype(m_buffer.sfc_flux_dir_vis)("sfc_flux_dir_vis", mem, m_ncol);
  mem += m_buffer.sfc_flux_dir_vis.totElems();
  m_buffer.sfc_flux_dir_nir = decltype(m_buffer.sfc_flux_dir_nir)("sfc_flux_dir_nir", mem, m_ncol);
  mem += m_buffer.sfc_flux_dir_nir.totElems();
  m_buffer.sfc_flux_dif_vis = decltype(m_buffer.sfc_flux_dif_vis)("sfc_flux_dif_vis", mem, m_ncol);
  mem += m_buffer.sfc_flux_dif_vis.totElems();
  m_buffer.sfc_flux_dif_nir = decltype(m_buffer.sfc_flux_dif_nir)("sfc_flux_dif_nir", mem, m_ncol);
  mem += m_buffer.sfc_flux_dif_nir.totElems();

  // 2d arrays
  m_buffer.p_lay = decltype(m_buffer.p_lay)("p_lay", mem, m_ncol, m_nlay);
  mem += m_buffer.p_lay.totElems();
  m_buffer.t_lay = decltype(m_buffer.t_lay)("t_lay", mem, m_ncol, m_nlay);
  mem += m_buffer.t_lay.totElems();
  m_buffer.p_del = decltype(m_buffer.p_del)("p_del", mem, m_ncol, m_nlay);
  mem += m_buffer.p_del.totElems();
  m_buffer.qc = decltype(m_buffer.qc)("qc", mem, m_ncol, m_nlay);
  mem += m_buffer.qc.totElems();
  m_buffer.qi = decltype(m_buffer.qi)("qi", mem, m_ncol, m_nlay);
  mem += m_buffer.qi.totElems();
  m_buffer.cldfrac_tot = decltype(m_buffer.cldfrac_tot)("cldfrac_tot", mem, m_ncol, m_nlay);
  mem += m_buffer.cldfrac_tot.totElems();
  m_buffer.eff_radius_qc = decltype(m_buffer.eff_radius_qc)("eff_radius_qc", mem, m_ncol, m_nlay);
  mem += m_buffer.eff_radius_qc.totElems();
  m_buffer.eff_radius_qi = decltype(m_buffer.eff_radius_qi)("eff_radius_qi", mem, m_ncol, m_nlay);
  mem += m_buffer.eff_radius_qi.totElems();
  m_buffer.tmp2d = decltype(m_buffer.tmp2d)("tmp2d", mem, m_ncol, m_nlay);
  mem += m_buffer.tmp2d.totElems();
  m_buffer.lwp = decltype(m_buffer.lwp)("lwp", mem, m_ncol, m_nlay);
  mem += m_buffer.lwp.totElems();
  m_buffer.iwp = decltype(m_buffer.iwp)("iwp", mem, m_ncol, m_nlay);
  mem += m_buffer.iwp.totElems();
  m_buffer.sw_heating = decltype(m_buffer.sw_heating)("sw_heating", mem, m_ncol, m_nlay);
  mem += m_buffer.sw_heating.totElems();
  m_buffer.lw_heating = decltype(m_buffer.lw_heating)("lw_heating", mem, m_ncol, m_nlay);
  mem += m_buffer.lw_heating.totElems();
  m_buffer.rad_heating = decltype(m_buffer.rad_heating)("rad_heating", mem, m_ncol, m_nlay);
  mem += m_buffer.rad_heating.totElems();
  // 3d arrays
  m_buffer.p_lev = decltype(m_buffer.p_lev)("p_lev", mem, m_ncol, m_nlay+1);
  mem += m_buffer.p_lev.totElems();
  m_buffer.t_lev = decltype(m_buffer.t_lev)("t_lev", mem, m_ncol, m_nlay+1);
  mem += m_buffer.t_lev.totElems();
  m_buffer.sw_flux_up = decltype(m_buffer.sw_flux_up)("sw_flux_up", mem, m_ncol, m_nlay+1);
  mem += m_buffer.sw_flux_up.totElems();
  m_buffer.sw_flux_dn = decltype(m_buffer.sw_flux_dn)("sw_flux_dn", mem, m_ncol, m_nlay+1);
  mem += m_buffer.sw_flux_dn.totElems();
  m_buffer.sw_flux_dn_dir = decltype(m_buffer.sw_flux_dn_dir)("sw_flux_dn_dir", mem, m_ncol, m_nlay+1);
  mem += m_buffer.sw_flux_dn_dir.totElems();
  m_buffer.lw_flux_up = decltype(m_buffer.lw_flux_up)("lw_flux_up", mem, m_ncol, m_nlay+1);
  mem += m_buffer.lw_flux_up.totElems();
  m_buffer.lw_flux_dn = decltype(m_buffer.lw_flux_dn)("lw_flux_dn", mem, m_ncol, m_nlay+1);
  mem += m_buffer.lw_flux_dn.totElems();
  // 3d nswbands
  m_buffer.sw_bnd_flux_up = decltype(m_buffer.sw_bnd_flux_up)("sw_bnd_flux_up", mem, m_ncol, m_nlay+1, m_nswbands);
  mem += m_buffer.sw_bnd_flux_up.totElems();
  m_buffer.sw_bnd_flux_dn = decltype(m_buffer.sw_bnd_flux_dn)("sw_bnd_flux_dn", mem, m_ncol, m_nlay+1, m_nswbands);
  mem += m_buffer.sw_bnd_flux_dn.totElems();
  m_buffer.sw_bnd_flux_dir = decltype(m_buffer.sw_bnd_flux_dir)("sw_bnd_flux_dir", mem, m_ncol, m_nlay+1, m_nswbands);
  mem += m_buffer.sw_bnd_flux_dir.totElems();
  m_buffer.sw_bnd_flux_dif = decltype(m_buffer.sw_bnd_flux_dif)("sw_bnd_flux_dif", mem, m_ncol, m_nlay+1, m_nswbands);
  mem += m_buffer.sw_bnd_flux_dif.totElems();
  // 3d nlwbands
  m_buffer.lw_bnd_flux_up = decltype(m_buffer.lw_bnd_flux_up)("lw_bnd_flux_up", mem, m_ncol, m_nlay+1, m_nlwbands);
  mem += m_buffer.lw_bnd_flux_up.totElems();
  m_buffer.lw_bnd_flux_dn = decltype(m_buffer.lw_bnd_flux_dn)("lw_bnd_flux_dn", mem, m_ncol, m_nlay+1, m_nlwbands);
  mem += m_buffer.lw_bnd_flux_dn.totElems();
  // Surface albedos
  m_buffer.sfc_alb_dir = decltype(m_buffer.sfc_alb_dir)("sfc_alb_dir", mem, m_ncol, m_nswbands);
  mem += m_buffer.sfc_alb_dir.totElems();
  m_buffer.sfc_alb_dif = decltype(m_buffer.sfc_alb_dif)("sfc_alb_dif", mem, m_ncol, m_nswbands);
  mem += m_buffer.sfc_alb_dif.totElems();

  m_buffer.aero_tau_sw = decltype(m_buffer.aero_tau_sw)("aero_tau_sw", mem, m_ncol, m_nlay, m_nswbands);
  mem += m_buffer.aero_tau_sw.totElems();
  m_buffer.aero_ssa_sw = decltype(m_buffer.aero_ssa_sw)("aero_ssa_sw", mem, m_ncol, m_nlay, m_nswbands);
  mem += m_buffer.aero_ssa_sw.totElems();
  m_buffer.aero_g_sw   = decltype(m_buffer.aero_g_sw  )("aero_g_sw"  , mem, m_ncol, m_nlay, m_nswbands);
  mem += m_buffer.aero_g_sw.totElems();
  m_buffer.aero_tau_lw = decltype(m_buffer.aero_tau_lw)("aero_tau_lw", mem, m_ncol, m_nlay, m_nlwbands);
  mem += m_buffer.aero_tau_lw.totElems();

  size_t used_mem = (reinterpret_cast<Real*>(mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for RRTMGPRadiation.");
} // RRTMGPRadiation::init_buffers

void RRTMGPRadiation::initialize_impl(const RunType /* run_type */) {
  using PC = scream::physics::Constants<Real>;

  // Determine rad timestep, specified as number of atm steps
  m_rad_freq_in_steps = m_params.get<Int>("rad_frequency", 1);

  // Determine orbital year. If Orbital Year is negative, use current year
  // from timestamp for orbital year; if positive, use provided orbital year
  // for duration of simulation.
  m_orbital_year = m_params.get<Int>("Orbital Year",-9999);
  // Get orbital parameters from yaml file
  m_orbital_eccen = m_params.get<Int>("Orbital Eccentricity",-9999);
  m_orbital_obliq = m_params.get<Int>("Orbital Obliquity"   ,-9999);
  m_orbital_mvelp = m_params.get<Int>("Orbital MVELP"       ,-9999);

  // Determine whether or not we are using a fixed solar zenith angle (positive value)
  m_fixed_solar_zenith_angle = m_params.get<Real>("Fixed Solar Zenith Angle", -9999);

  // Whether or not to do MCICA subcolumn sampling
  m_do_subcol_sampling = m_params.get<bool>("do_subcol_sampling",false);
  EKAT_REQUIRE_MSG(not m_do_subcol_sampling, "Error! RRTMGP does not yet support do_subcol_sampling = true");

  // Initialize yakl
  if(!yakl::isInitialized()) { yakl::init(); }

  // Names of active gases
  auto gas_names_yakl_offset = string1d("gas_names",m_ngas);
  m_gas_mol_weights          = view_1d_real("gas_mol_weights",m_ngas);
  /* the lookup function for getting the gas mol weights doesn't work on device. */
  auto gas_mol_w_host = Kokkos::create_mirror_view(m_gas_mol_weights);
  for (int igas = 0; igas < m_ngas; igas++) {  
    /* Note: YAKL starts the index from 1 */
    gas_names_yakl_offset(igas+1)   = m_gas_names[igas];
    gas_mol_w_host[igas]            = PC::get_gas_mol_weight(m_gas_names[igas]);
  }
  Kokkos::deep_copy(m_gas_mol_weights,gas_mol_w_host);

  // Initialize GasConcs object to pass to RRTMGP initializer;
  gas_concs.init(gas_names_yakl_offset,m_ncol,m_nlay);
  rrtmgp::rrtmgp_initialize(gas_concs, m_atm_logger);

  // Set property checks for fields in this process

  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("T_mid"),m_grid,140.0, 500.0,false);
  add_precondition_check<FieldWithinIntervalCheck>(get_field_out("sfc_alb_dir_vis"),m_grid,0.0,1.0,true);
  add_precondition_check<FieldWithinIntervalCheck>(get_field_out("sfc_alb_dir_nir"),m_grid,0.0,1.0,true);
  add_precondition_check<FieldWithinIntervalCheck>(get_field_out("sfc_alb_dif_vis"),m_grid,0.0,1.0,true);
  add_precondition_check<FieldWithinIntervalCheck>(get_field_out("sfc_alb_dif_nir"),m_grid,0.0,1.0,true);
}

// =========================================================================================

void RRTMGPRadiation::run_impl (const int dt) {
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using PC = scream::physics::Constants<Real>;
  using CO = scream::ColumnOps<DefaultDevice,Real>;

  // get a host copy of lat/lon 
  auto h_lat  = Kokkos::create_mirror_view(m_lat);
  auto h_lon  = Kokkos::create_mirror_view(m_lon);
  Kokkos::deep_copy(h_lat,m_lat);
  Kokkos::deep_copy(h_lon,m_lon);

  // Get data from the FieldManager
  auto d_pmid = get_field_in("p_mid").get_view<const Real**>();
  auto d_pint = get_field_in("p_int").get_view<const Real**>();
  auto d_pdel = get_field_in("pseudo_density").get_view<const Real**>();
  auto d_sfc_alb_dir_vis = get_field_in("sfc_alb_dir_vis").get_view<const Real*>();
  auto d_sfc_alb_dir_nir = get_field_in("sfc_alb_dir_nir").get_view<const Real*>();
  auto d_sfc_alb_dif_vis = get_field_in("sfc_alb_dif_vis").get_view<const Real*>();
  auto d_sfc_alb_dif_nir = get_field_in("sfc_alb_dif_nir").get_view<const Real*>();
  auto d_qv = get_field_in("qv").get_view<const Real**>();
  auto d_qc = get_field_in("qc").get_view<const Real**>();
  auto d_qi = get_field_in("qi").get_view<const Real**>();
  auto d_cldfrac_tot = get_field_in("cldfrac_tot").get_view<const Real**>();
  auto d_rel = get_field_in("eff_radius_qc").get_view<const Real**>();
  auto d_rei = get_field_in("eff_radius_qi").get_view<const Real**>();
  auto d_surf_lw_flux_up = get_field_in("surf_lw_flux_up").get_view<const Real*>();
  // Output fields
  auto d_tmid = get_field_out("T_mid").get_view<Real**>();
  using SmallPack = ekat::Pack<Real,SCREAM_SMALL_PACK_SIZE>;
  const int n_lay_w_pack = SCREAM_SMALL_PACK_SIZE*ekat::npack<SmallPack>(m_nlay);
  view_3d_real d_aero_tau_sw("aero_tau_sw",m_ncol,m_nswbands,n_lay_w_pack);
  view_3d_real d_aero_ssa_sw("aero_ssa_sw",m_ncol,m_nswbands,n_lay_w_pack);
  view_3d_real d_aero_g_sw  ("aero_g_sw"  ,m_ncol,m_nswbands,n_lay_w_pack);
  view_3d_real d_aero_tau_lw("aero_tau_lw",m_ncol,m_nlwbands,n_lay_w_pack);
  if (m_do_aerosol_rad) {
    Kokkos::deep_copy(d_aero_tau_sw,get_field_in("aero_tau_sw").get_view<const Real***>());
    Kokkos::deep_copy(d_aero_ssa_sw,get_field_in("aero_ssa_sw").get_view<const Real***>());
    Kokkos::deep_copy(d_aero_g_sw  ,get_field_in("aero_g_sw"  ).get_view<const Real***>());
    Kokkos::deep_copy(d_aero_tau_lw,get_field_in("aero_tau_lw").get_view<const Real***>());
  } else {
    Kokkos::deep_copy(d_aero_tau_sw,0.0);
    Kokkos::deep_copy(d_aero_ssa_sw,0.0);
    Kokkos::deep_copy(d_aero_g_sw  ,0.0);
    Kokkos::deep_copy(d_aero_tau_lw,0.0);
    
  } 
  auto d_sw_flux_up = get_field_out("SW_flux_up").get_view<Real**>();
  auto d_sw_flux_dn = get_field_out("SW_flux_dn").get_view<Real**>();
  auto d_sw_flux_dn_dir = get_field_out("SW_flux_dn_dir").get_view<Real**>();
  auto d_lw_flux_up = get_field_out("LW_flux_up").get_view<Real**>();
  auto d_lw_flux_dn = get_field_out("LW_flux_dn").get_view<Real**>();
  auto d_rad_heating_pdel = get_field_out("rad_heating_pdel").get_view<Real**>();
  auto d_sfc_flux_dir_vis = get_field_out("sfc_flux_dir_vis").get_view<Real*>();
  auto d_sfc_flux_dir_nir = get_field_out("sfc_flux_dir_nir").get_view<Real*>();
  auto d_sfc_flux_dif_vis = get_field_out("sfc_flux_dif_vis").get_view<Real*>();
  auto d_sfc_flux_dif_nir = get_field_out("sfc_flux_dif_nir").get_view<Real*>();
  auto d_sfc_flux_sw_net = get_field_out("sfc_flux_sw_net").get_view<Real*>();
  auto d_sfc_flux_lw_dn  = get_field_out("sfc_flux_lw_dn").get_view<Real*>();

  // Create YAKL arrays. RRTMGP expects YAKL arrays with styleFortran, i.e., data has ncol
  // as the fastest index. For this reason we must copy the data.
  auto p_lay           = m_buffer.p_lay;
  auto t_lay           = m_buffer.t_lay;
  auto p_lev           = m_buffer.p_lev;
  auto p_del           = m_buffer.p_del;
  auto t_lev           = m_buffer.t_lev;
  auto mu0             = m_buffer.mu0;
  auto sfc_alb_dir     = m_buffer.sfc_alb_dir;
  auto sfc_alb_dif     = m_buffer.sfc_alb_dif;
  auto sfc_alb_dir_vis = m_buffer.sfc_alb_dir_vis;
  auto sfc_alb_dir_nir = m_buffer.sfc_alb_dir_nir;
  auto sfc_alb_dif_vis = m_buffer.sfc_alb_dif_vis;
  auto sfc_alb_dif_nir = m_buffer.sfc_alb_dif_nir;
  auto qc              = m_buffer.qc;
  auto qi              = m_buffer.qi;
  auto cldfrac_tot     = m_buffer.cldfrac_tot;
  auto rel             = m_buffer.eff_radius_qc;
  auto rei             = m_buffer.eff_radius_qi;
  auto sw_flux_up      = m_buffer.sw_flux_up;
  auto sw_flux_dn      = m_buffer.sw_flux_dn;
  auto sw_flux_dn_dir  = m_buffer.sw_flux_dn_dir;
  auto lw_flux_up      = m_buffer.lw_flux_up;
  auto lw_flux_dn      = m_buffer.lw_flux_dn;
  auto sw_bnd_flux_up  = m_buffer.sw_bnd_flux_up;
  auto sw_bnd_flux_dn  = m_buffer.sw_bnd_flux_dn;
  auto sw_bnd_flux_dir = m_buffer.sw_bnd_flux_dir;
  auto sw_bnd_flux_dif = m_buffer.sw_bnd_flux_dif;
  auto lw_bnd_flux_up  = m_buffer.lw_bnd_flux_up;
  auto lw_bnd_flux_dn  = m_buffer.lw_bnd_flux_dn;
  auto sfc_flux_dir_vis = m_buffer.sfc_flux_dir_vis;
  auto sfc_flux_dir_nir = m_buffer.sfc_flux_dir_nir;
  auto sfc_flux_dif_vis = m_buffer.sfc_flux_dif_vis;
  auto sfc_flux_dif_nir = m_buffer.sfc_flux_dif_nir;
  auto aero_tau_sw     = m_buffer.aero_tau_sw;
  auto aero_ssa_sw     = m_buffer.aero_ssa_sw;
  auto aero_g_sw       = m_buffer.aero_g_sw;
  auto aero_tau_lw     = m_buffer.aero_tau_lw;

  constexpr auto stebol = PC::stebol;
  const auto ncol = m_ncol;
  const auto nlay = m_nlay;
  const auto nlwbands = m_nlwbands;
  const auto nswbands = m_nswbands;

  // Compute orbital parameters; these are used both for computing
  // the solar zenith angle and also for computing total solar
  // irradiance scaling (tsi_scaling).
  double obliqr, lambm0, mvelpp;
  auto ts = timestamp();
  auto orbital_year = m_orbital_year;
  auto eccen = m_orbital_eccen;
  auto obliq = m_orbital_obliq;
  auto mvelp = m_orbital_mvelp;
  if (eccen >= 0 && obliq >= 0 && mvelp >= 0) {
    // use fixed oribal parameters; to force this, we need to set
    // orbital_year to SHR_ORB_UNDEF_INT, which is exposed through
    // our c2f bridge as shr_orb_undef_int_c2f
    orbital_year = shr_orb_undef_int_c2f;
  } else if (orbital_year < 0) {
    // compute orbital parameters based on current year
    orbital_year = ts.get_year();
  }
  shr_orb_params_c2f(&orbital_year, &eccen, &obliq, &mvelp, 
                     &obliqr, &lambm0, &mvelpp);
  // Use the orbital parameters to calculate the solar declination and eccentricity factor
  double delta, eccf;
  auto calday = ts.frac_of_year_in_days();
  shr_orb_decl_c2f(calday, eccen, mvelpp, lambm0,
                   obliqr, &delta, &eccf);

  // Are we going to update fluxes and heating this step?
  auto update_rad = scream::rrtmgp::radiation_do(m_rad_freq_in_steps, ts.get_num_steps());

  // Copy data from the FieldManager to the YAKL arrays
  {
    // Determine the cosine zenith angle
    // NOTE: Since we are bridging to F90 arrays this must be done on HOST and then
    //       deep copied to a device view.
    auto d_mu0 = m_buffer.cosine_zenith;
    auto h_mu0 = Kokkos::create_mirror_view(d_mu0);
    if (m_fixed_solar_zenith_angle > 0) {
      for (int i=0; i<m_ncol; i++) {
        h_mu0(i) = m_fixed_solar_zenith_angle;
      }
    } else {
      // Now use solar declination to calculate zenith angle for all points
      for (int i=0;i<m_ncol;i++) {
        double lat = h_lat(i)*PC::Pi/180.0;  // Convert lat/lon to radians
        double lon = h_lon(i)*PC::Pi/180.0;
        h_mu0(i) = shr_orb_cosz_c2f(calday, lat, lon, delta, dt);
      }
    }
    Kokkos::deep_copy(d_mu0,h_mu0);

    // dz and T_int will need to be computed
    view_2d_real d_tint("T_int", m_ncol, m_nlay+1);
    view_2d_real d_dz  ("dz",    m_ncol, m_nlay);

    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_ncol, m_nlay);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();

      // Calculate dz
      const auto pseudo_density = ekat::subview(d_pdel, i);
      const auto p_mid          = ekat::subview(d_pmid, i);
      const auto T_mid          = ekat::subview(d_tmid, i);
      const auto qv             = ekat::subview(d_qv,   i);
      const auto dz             = ekat::subview(d_dz,   i);
      PF::calculate_dz<Real>(team, pseudo_density, p_mid, T_mid, qv, dz);
      team.team_barrier();

      // Calculate T_int from longwave flux up from the surface, assuming
      // blackbody emission with emissivity of 1.
      // TODO: Does land model assume something other than emissivity of 1? If so
      // we should use that here rather than assuming perfect blackbody emission.
      // NOTE: RRTMGP can accept vertical ordering surface to toa, or toa to
      // surface. The input data for the standalone test is ordered surface to
      // toa, but SCREAM in general assumes data is toa to surface. We account
      // for this here by swapping bc_top and bc_bot in the case that the input
      // data is ordered surface to toa.
      const auto T_int = ekat::subview(d_tint, i);
      const auto P_mid = ekat::subview(d_pmid, i);
      const int itop = (P_mid(0) < P_mid(nlay-1)) ? 0 : nlay-1;
      const Real bc_top = T_mid(itop);
      const Real bc_bot = sqrt(sqrt(d_surf_lw_flux_up(i)/stebol));
      if (itop == 0) {
          CO::compute_interface_values_linear(team, nlay, T_mid, dz, bc_top, bc_bot, T_int);
      } else {
          CO::compute_interface_values_linear(team, nlay, T_mid, dz, bc_bot, bc_top, T_int);
      }
      team.team_barrier();

      mu0(i+1) = d_mu0(i);
      sfc_alb_dir_vis(i+1) = d_sfc_alb_dir_vis(i);
      sfc_alb_dir_nir(i+1) = d_sfc_alb_dir_nir(i);
      sfc_alb_dif_vis(i+1) = d_sfc_alb_dif_vis(i);
      sfc_alb_dif_nir(i+1) = d_sfc_alb_dif_nir(i);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlay), [&] (const int& k) {
        p_lay(i+1,k+1)       = d_pmid(i,k);
        t_lay(i+1,k+1)       = d_tmid(i,k);
        p_del(i+1,k+1)       = d_pdel(i,k);
        qc(i+1,k+1)          = d_qc(i,k);
        qi(i+1,k+1)          = d_qi(i,k);
        rel(i+1,k+1)         = d_rel(i,k);
        rei(i+1,k+1)         = d_rei(i,k);
        p_lev(i+1,k+1)       = d_pint(i,k);
        t_lev(i+1,k+1)       = d_tint(i,k);
      });

      p_lev(i+1,nlay+1) = d_pint(i,nlay);
      t_lev(i+1,nlay+1) = d_tint(i,nlay);

      // Note that RRTMGP expects ordering (col,lay,bnd) but the FM keeps things in (col,bnd,lay) order
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nswbands*nlay), [&] (const int&idx) {
          auto b = idx / nlay;
          auto k = idx % nlay;
          aero_tau_sw(i+1,k+1,b+1) = d_aero_tau_sw(i,b,k);
          aero_ssa_sw(i+1,k+1,b+1) = d_aero_ssa_sw(i,b,k);
          aero_g_sw  (i+1,k+1,b+1) = d_aero_g_sw  (i,b,k);
      });
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlwbands*nlay), [&] (const int&idx) {
          auto b = idx / nlay;
          auto k = idx % nlay;
          aero_tau_lw(i+1,k+1,b+1) = d_aero_tau_lw(i,b,k);
      });
    });
  }
  Kokkos::fence();

  // Populate GasConcs object to pass to RRTMGP driver
  auto tmp2d = m_buffer.tmp2d;
  for (int igas = 0; igas < m_ngas; igas++) {
    auto name = m_gas_names[igas];
    auto fm_name = name=="h2o" ? "qv" : name;
    auto d_temp  = get_field_in(fm_name).get_view<const Real**>();
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_nlay, m_ncol);
    const auto gas_mol_weights = m_gas_mol_weights;

    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int k = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ncol), [&] (const int& i) {
        tmp2d(i+1,k+1) = PF::calculate_vmr_from_mmr(gas_mol_weights[igas],d_qv(i,k),d_temp(i,k)); // Note that for YAKL arrays i and k start with index 1
      });
    });
    Kokkos::fence();

    gas_concs.set_vmr(name, tmp2d);
  }

  // Set layer cloud fraction.
  //
  // If not doing subcolumn sampling for mcica, we want to make sure we use grid-mean
  // condensate for computing cloud optical properties, because we are assuming the
  // entire column is completely clear or cloudy. Thus, in this case we want to set
  // cloud fraction to 0 or 1. Note that we could choose an alternative threshold
  // criteria here, like qc + qi > 1e-5 or something.
  //
  // If we *are* doing subcolumn sampling for MCICA, then keep cloud fraction as input
  // from cloud fraction parameterization, wherever that is computed.
  auto do_subcol_sampling = m_do_subcol_sampling;
  if (not do_subcol_sampling) {
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_ncol, m_nlay);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlay), [&] (const int& k) {
        if (d_cldfrac_tot(i,k) > 0) {
          cldfrac_tot(i+1,k+1) = 1;
        } else {
          cldfrac_tot(i+1,k+1) = 0;
        }
      });
    });
  } else {
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_ncol, m_nlay);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlay), [&] (const int& k) {
        cldfrac_tot(i+1,k+1) = d_cldfrac_tot(i,k);
      });
    });
  }
  Kokkos::fence();

  // Compute layer cloud mass (per unit area)
  auto lwp = m_buffer.lwp;
  auto iwp = m_buffer.iwp;
  scream::rrtmgp::mixing_ratio_to_cloud_mass(qc, cldfrac_tot, p_del, lwp);
  scream::rrtmgp::mixing_ratio_to_cloud_mass(qi, cldfrac_tot, p_del, iwp);
  // Convert to g/m2 (needed by RRTMGP)
  {
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_nlay, m_ncol);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int k = team.league_rank()+1; // Note that for YAKL arrays i and k start with index 1
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ncol), [&] (const int& icol) {
      int i = icol+1;
      lwp(i,k) *= 1e3;
      iwp(i,k) *= 1e3;
    });
  });
  }
  Kokkos::fence();

  // Compute band-by-band surface_albedos. This is needed since
  // the AD passes broadband albedos, but rrtmgp require band-by-band.
  rrtmgp::compute_band_by_band_surface_albedos(
    m_ncol, m_nswbands,
    sfc_alb_dir_vis, sfc_alb_dir_nir,
    sfc_alb_dif_vis, sfc_alb_dif_nir,
    sfc_alb_dir, sfc_alb_dif);

  // Run RRTMGP driver
  if (update_rad) {
    rrtmgp::rrtmgp_main(
      m_ncol, m_nlay,
      p_lay, t_lay, p_lev, t_lev,
      gas_concs,
      sfc_alb_dir, sfc_alb_dif, mu0,
      lwp, iwp, rel, rei,
      aero_tau_sw, aero_ssa_sw, aero_g_sw,
      aero_tau_lw,
      sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
      lw_flux_up, lw_flux_dn, 
      sw_bnd_flux_up, sw_bnd_flux_dn, sw_bnd_flux_dir,
      lw_bnd_flux_up, lw_bnd_flux_dn, 
      eccf, m_atm_logger
    );
  }

  // Compute and apply heating tendency
  auto sw_heating  = m_buffer.sw_heating;
  auto lw_heating  = m_buffer.lw_heating;
  auto rad_heating = m_buffer.rad_heating;
  if (update_rad) {
    rrtmgp::compute_heating_rate(
      sw_flux_up, sw_flux_dn, p_del, sw_heating
    );
    rrtmgp::compute_heating_rate(
      lw_flux_up, lw_flux_dn, p_del, lw_heating
    );
    {
      const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_nlay, m_ncol);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int k = team.league_rank()+1; // Note that for YAKL arrays i and k start with index 1
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ncol), [&] (const int& icol) {
          int i = icol+1;
          // Combine SW and LW heating into a net heating tendency
          rad_heating(i,k) = sw_heating(i,k) + lw_heating(i,k);
          // Apply heating tendency to temperature
          t_lay(i,k) = t_lay(i,k) + rad_heating(i,k) * dt;
          // Compute pdel * rad_heating to conserve energy across timesteps in case
          // levels change before next heating update.
          d_rad_heating_pdel(icol,k-1) = d_pdel(icol,k-1) * rad_heating(i,k);
        });
      });
    }
    Kokkos::fence();
  } else {
    {
      const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_nlay, m_ncol);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int k = team.league_rank()+1; // Note that for YAKL arrays i and k start with index 1
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ncol), [&] (const int& icol) {
          int i = icol+1;
          // Back out net heating tendency from the pdel weighted quantity we keep in the FM
          rad_heating(i,k) = d_rad_heating_pdel(icol,k-1) / d_pdel(icol,k-1);
          // Apply heating tendency to temperature
          t_lay(i,k) = t_lay(i,k) + rad_heating(i,k) * dt;
        });
      });
    }
    Kokkos::fence();
  }

  // Index to surface (bottom of model); used to get surface fluxes below
  const int kbot = nlay+1;

  // Compute diffuse flux as difference between total and direct; use YAKL parallel_for here because these are YAKL objects
  parallel_for(Bounds<3>(m_nswbands,m_nlay+1,m_ncol), YAKL_LAMBDA(int ibnd, int ilev, int icol) {
    sw_bnd_flux_dif(icol,ilev,ibnd) = sw_bnd_flux_dn(icol,ilev,ibnd) - sw_bnd_flux_dir(icol,ilev,ibnd);
  });

  // Compute surface fluxes
  rrtmgp::compute_broadband_surface_fluxes(
      m_ncol, kbot, m_nswbands,
      sw_bnd_flux_dir, sw_bnd_flux_dif, 
      sfc_flux_dir_vis, sfc_flux_dir_nir, 
      sfc_flux_dif_vis, sfc_flux_dif_nir
  );

  // Copy output data back to FieldManager
  if (update_rad) {
    {
      const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_ncol, m_nlay);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int i = team.league_rank();
        d_sfc_flux_dir_nir(i) = sfc_flux_dir_nir(i+1);
        d_sfc_flux_dir_vis(i) = sfc_flux_dir_vis(i+1);
        d_sfc_flux_dif_nir(i) = sfc_flux_dif_nir(i+1);
        d_sfc_flux_dif_vis(i) = sfc_flux_dif_vis(i+1);
        d_sfc_flux_sw_net(i)  = sw_flux_dn(i+1,kbot) - sw_flux_up(i+1,kbot);
        d_sfc_flux_lw_dn(i)   = lw_flux_dn(i+1,kbot);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlay+1), [&] (const int& k) {
          d_sw_flux_up(i,k)     = sw_flux_up(i+1,k+1);
          d_sw_flux_dn(i,k)     = sw_flux_dn(i+1,k+1);
          d_sw_flux_dn_dir(i,k) = sw_flux_dn_dir(i+1,k+1);
          d_lw_flux_up(i,k)     = lw_flux_up(i+1,k+1);
          d_lw_flux_dn(i,k)     = lw_flux_dn(i+1,k+1);
        });
      });
    }
  }
  // Temperature is always updated
  {
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_ncol, m_nlay);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlay), [&] (const int& k) {
        d_tmid(i,k) = t_lay(i+1,k+1);
      });
    });
  }
   
}
// =========================================================================================

void RRTMGPRadiation::finalize_impl  () {
  gas_concs.reset();
  rrtmgp::rrtmgp_finalize();

  // Finalize YAKL
  yakl::finalize();
}
// =========================================================================================


}  // namespace scream
