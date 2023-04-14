#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "physics/rrtmgp/atmosphere_radiation.hpp"
#include "physics/rrtmgp/rrtmgp_utils.hpp"
#include "physics/rrtmgp/shr_orb_mod_c2f.hpp"
#include "physics/share/scream_trcmix.hpp"
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
  auto m2 = m * m;
  auto Wm2 = W / m / m;
  Wm2.set_string("W/m2");
  auto nondim = m/m;  // dummy unit for non-dimensional fields
  auto micron = m / 1000000;
  auto molmol = mol/mol;
  molmol.set_string("mol/mol");

  using namespace ShortFieldTagsNames;

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  m_ncol = m_grid->get_num_local_dofs();
  m_nlay = m_grid->get_num_vertical_levels();
  m_lat  = m_grid->get_geometry_data("lat");
  m_lon  = m_grid->get_geometry_data("lon");

  // Figure out radiation column chunks stats
  m_col_chunk_size = std::min(m_params.get("column_chunk_size", m_ncol),m_ncol);
  m_num_col_chunks = (m_ncol+m_col_chunk_size-1) / m_col_chunk_size;
  m_col_chunk_beg.resize(m_num_col_chunks+1,0);
  for (int i=0; i<m_num_col_chunks; ++i) {
    m_col_chunk_beg[i+1] = std::min(m_ncol,m_col_chunk_beg[i] + m_col_chunk_size);
  }
  this->log(LogLevel::debug,
            "[RRTMGP::set_grids] Col chunking stats:\n"
            "  - Chunk size: " + std::to_string(m_col_chunk_size) + "\n"
            "  - Number of chunks: " + std::to_string(m_num_col_chunks) + "\n");

  // Set up dimension layouts
  m_nswgpts = m_params.get<int>("nswgpts",112);
  m_nlwgpts = m_params.get<int>("nlwgpts",128);
  FieldLayout scalar2d_layout     { {COL   }, {m_ncol    } };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_ncol,m_nlay} };
  FieldLayout scalar3d_layout_int { {COL,ILEV}, {m_ncol,m_nlay+1} };
  FieldLayout scalar3d_swband_layout { {COL,SWBND,LEV}, {m_ncol, m_nswbands, m_nlay} };
  FieldLayout scalar3d_lwband_layout { {COL,LWBND,LEV}, {m_ncol, m_nlwbands, m_nlay} };
  FieldLayout scalar3d_swgpts_layout { {COL,SWGPT,LEV}, {m_ncol, m_nswgpts, m_nlay} };
  FieldLayout scalar3d_lwgpts_layout { {COL,LWGPT,LEV}, {m_ncol, m_nlwgpts, m_nlay} };

  constexpr int ps = SCREAM_SMALL_PACK_SIZE;
  // Set required (input) fields here
  add_field<Required>("p_mid" , scalar3d_layout_mid, Pa, grid_name, ps);
  add_field<Required>("p_int", scalar3d_layout_int, Pa, grid_name, ps);
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa, grid_name, ps);
  add_field<Required>("sfc_alb_dir_vis", scalar2d_layout, nondim, grid_name);
  add_field<Required>("sfc_alb_dir_nir", scalar2d_layout, nondim, grid_name);
  add_field<Required>("sfc_alb_dif_vis", scalar2d_layout, nondim, grid_name);
  add_field<Required>("sfc_alb_dif_nir", scalar2d_layout, nondim, grid_name);
  add_field<Required>("qc", scalar3d_layout_mid, kgkg, grid_name, ps);
  add_field<Required>("qi", scalar3d_layout_mid, kgkg, grid_name, ps);
  add_field<Required>("cldfrac_tot", scalar3d_layout_mid, nondim, grid_name, ps);
  add_field<Required>("eff_radius_qc", scalar3d_layout_mid, micron, grid_name, ps);
  add_field<Required>("eff_radius_qi", scalar3d_layout_mid, micron, grid_name, ps);
  add_field<Required>("qv",scalar3d_layout_mid,kgkg,grid_name, ps);
  add_field<Required>("surf_lw_flux_up",scalar2d_layout,W/(m*m),grid_name);
  // Set of required gas concentration fields
  for (auto& it : m_gas_names) {
    // Add gas VOLUME mixing ratios (moles of gas / moles of air; what actually gets input to RRTMGP)
    if (it == "o3") {
      // o3 is read from file, or computed by chemistry
      add_field<Updated >(it + "_volume_mix_ratio", scalar3d_layout_mid, molmol, grid_name, ps);
    } else {
      // the rest are computed from prescribed surface values
      add_field<Computed>(it + "_volume_mix_ratio", scalar3d_layout_mid, molmol, grid_name, ps);
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
  add_field<Computed>("SW_clrsky_flux_dn", scalar3d_layout_int, Wm2, grid_name, ps);
  add_field<Computed>("SW_clrsky_flux_up", scalar3d_layout_int, Wm2, grid_name, ps);
  add_field<Computed>("SW_clrsky_flux_dn_dir", scalar3d_layout_int, Wm2, grid_name, ps);
  add_field<Computed>("LW_clrsky_flux_up", scalar3d_layout_int, Wm2, grid_name, ps);
  add_field<Computed>("LW_clrsky_flux_dn", scalar3d_layout_int, Wm2, grid_name, ps);
  add_field<Computed>("rad_heating_pdel", scalar3d_layout_mid, Pa*K/s, grid_name, ps);
  // Cloud properties added as computed fields for diagnostic purposes
  add_field<Computed>("cldlow"        , scalar2d_layout, nondim, grid_name);
  add_field<Computed>("cldmed"        , scalar2d_layout, nondim, grid_name);
  add_field<Computed>("cldhgh"        , scalar2d_layout, nondim, grid_name);
  add_field<Computed>("cldtot"        , scalar2d_layout, nondim, grid_name);

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

  // Boundary flux fields for energy and mass conservation checks
  if (has_column_conservation_check()) {
    add_field<Computed>("vapor_flux", scalar2d_layout, kg/m2/s, grid_name);
    add_field<Computed>("water_flux", scalar2d_layout, m/s,     grid_name);
    add_field<Computed>("ice_flux",   scalar2d_layout, m/s,     grid_name);
    add_field<Computed>("heat_flux",  scalar2d_layout, W/m2,    grid_name);
  }
}  // RRTMGPRadiation::set_grids

size_t RRTMGPRadiation::requested_buffer_size_in_bytes() const
{
  const size_t interface_request =
    Buffer::num_1d_ncol*m_col_chunk_size +
    Buffer::num_2d_nlay*m_col_chunk_size*m_nlay +
    Buffer::num_2d_nlay_p1*m_col_chunk_size*(m_nlay+1) +
    Buffer::num_2d_nswbands*m_col_chunk_size*m_nswbands +
    Buffer::num_3d_nlev_nswbands*m_col_chunk_size*(m_nlay+1)*m_nswbands +
    Buffer::num_3d_nlev_nlwbands*m_col_chunk_size*(m_nlay+1)*m_nlwbands +
    Buffer::num_3d_nlay_nswbands*m_col_chunk_size*(m_nlay)*m_nswbands +
    Buffer::num_3d_nlay_nlwbands*m_col_chunk_size*(m_nlay)*m_nlwbands +
    Buffer::num_3d_nlay_nswgpts*m_col_chunk_size*(m_nlay)*m_nswgpts +
    Buffer::num_3d_nlay_nlwgpts*m_col_chunk_size*(m_nlay)*m_nlwgpts;

  return interface_request * sizeof(Real);
} // RRTMGPRadiation::requested_buffer_size
// =========================================================================================

void RRTMGPRadiation::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d arrays
  m_buffer.mu0 = decltype(m_buffer.mu0)("mu0", mem, m_col_chunk_size);
  mem += m_buffer.mu0.totElems();
  m_buffer.sfc_alb_dir_vis = decltype(m_buffer.sfc_alb_dir_vis)("sfc_alb_dir_vis", mem, m_col_chunk_size);
  mem += m_buffer.sfc_alb_dir_vis.totElems();
  m_buffer.sfc_alb_dir_nir = decltype(m_buffer.sfc_alb_dir_nir)("sfc_alb_dir_nir", mem, m_col_chunk_size);
  mem += m_buffer.sfc_alb_dir_nir.totElems();
  m_buffer.sfc_alb_dif_vis = decltype(m_buffer.sfc_alb_dif_vis)("sfc_alb_dif_vis", mem, m_col_chunk_size);
  mem += m_buffer.sfc_alb_dif_vis.totElems();
  m_buffer.sfc_alb_dif_nir = decltype(m_buffer.sfc_alb_dif_nir)("sfc_alb_dif_nir", mem, m_col_chunk_size);
  mem += m_buffer.sfc_alb_dif_nir.totElems();
  m_buffer.sfc_flux_dir_vis = decltype(m_buffer.sfc_flux_dir_vis)("sfc_flux_dir_vis", mem, m_col_chunk_size);
  mem += m_buffer.sfc_flux_dir_vis.totElems();
  m_buffer.sfc_flux_dir_nir = decltype(m_buffer.sfc_flux_dir_nir)("sfc_flux_dir_nir", mem, m_col_chunk_size);
  mem += m_buffer.sfc_flux_dir_nir.totElems();
  m_buffer.sfc_flux_dif_vis = decltype(m_buffer.sfc_flux_dif_vis)("sfc_flux_dif_vis", mem, m_col_chunk_size);
  mem += m_buffer.sfc_flux_dif_vis.totElems();
  m_buffer.sfc_flux_dif_nir = decltype(m_buffer.sfc_flux_dif_nir)("sfc_flux_dif_nir", mem, m_col_chunk_size);
  mem += m_buffer.sfc_flux_dif_nir.totElems();
  m_buffer.cosine_zenith = decltype(m_buffer.cosine_zenith)(mem, m_col_chunk_size);
  mem += m_buffer.cosine_zenith.size();

  // 2d arrays
  m_buffer.p_lay = decltype(m_buffer.p_lay)("p_lay", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.p_lay.totElems();
  m_buffer.t_lay = decltype(m_buffer.t_lay)("t_lay", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.t_lay.totElems();
  m_buffer.p_del = decltype(m_buffer.p_del)("p_del", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.p_del.totElems();
  m_buffer.qc = decltype(m_buffer.qc)("qc", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.qc.totElems();
  m_buffer.qi = decltype(m_buffer.qi)("qi", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.qi.totElems();
  m_buffer.cldfrac_tot = decltype(m_buffer.cldfrac_tot)("cldfrac_tot", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.cldfrac_tot.totElems();
  m_buffer.eff_radius_qc = decltype(m_buffer.eff_radius_qc)("eff_radius_qc", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.eff_radius_qc.totElems();
  m_buffer.eff_radius_qi = decltype(m_buffer.eff_radius_qi)("eff_radius_qi", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.eff_radius_qi.totElems();
  m_buffer.tmp2d = decltype(m_buffer.tmp2d)("tmp2d", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.tmp2d.totElems();
  m_buffer.lwp = decltype(m_buffer.lwp)("lwp", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.lwp.totElems();
  m_buffer.iwp = decltype(m_buffer.iwp)("iwp", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.iwp.totElems();
  m_buffer.sw_heating = decltype(m_buffer.sw_heating)("sw_heating", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.sw_heating.totElems();
  m_buffer.lw_heating = decltype(m_buffer.lw_heating)("lw_heating", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.lw_heating.totElems();
  // 3d arrays
  m_buffer.p_lev = decltype(m_buffer.p_lev)("p_lev", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.p_lev.totElems();
  m_buffer.t_lev = decltype(m_buffer.t_lev)("t_lev", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.t_lev.totElems();
  m_buffer.sw_flux_up = decltype(m_buffer.sw_flux_up)("sw_flux_up", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_flux_up.totElems();
  m_buffer.sw_flux_dn = decltype(m_buffer.sw_flux_dn)("sw_flux_dn", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_flux_dn.totElems();
  m_buffer.sw_flux_dn_dir = decltype(m_buffer.sw_flux_dn_dir)("sw_flux_dn_dir", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_flux_dn_dir.totElems();
  m_buffer.lw_flux_up = decltype(m_buffer.lw_flux_up)("lw_flux_up", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_flux_up.totElems();
  m_buffer.lw_flux_dn = decltype(m_buffer.lw_flux_dn)("lw_flux_dn", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_flux_dn.totElems();
  m_buffer.sw_clrsky_flux_up = decltype(m_buffer.sw_clrsky_flux_up)("sw_clrsky_flux_up", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clrsky_flux_up.totElems();
  m_buffer.sw_clrsky_flux_dn = decltype(m_buffer.sw_clrsky_flux_dn)("sw_clrsky_flux_dn", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clrsky_flux_dn.totElems();
  m_buffer.sw_clrsky_flux_dn_dir = decltype(m_buffer.sw_clrsky_flux_dn_dir)("sw_clrsky_flux_dn_dir", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clrsky_flux_dn_dir.totElems();
  m_buffer.lw_clrsky_flux_up = decltype(m_buffer.lw_clrsky_flux_up)("lw_clrsky_flux_up", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clrsky_flux_up.totElems();
  m_buffer.lw_clrsky_flux_dn = decltype(m_buffer.lw_clrsky_flux_dn)("lw_clrsky_flux_dn", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clrsky_flux_dn.totElems();
  // 3d arrays with nswbands dimension (shortwave fluxes by band)
  m_buffer.sw_bnd_flux_up = decltype(m_buffer.sw_bnd_flux_up)("sw_bnd_flux_up", mem, m_col_chunk_size, m_nlay+1, m_nswbands);
  mem += m_buffer.sw_bnd_flux_up.totElems();
  m_buffer.sw_bnd_flux_dn = decltype(m_buffer.sw_bnd_flux_dn)("sw_bnd_flux_dn", mem, m_col_chunk_size, m_nlay+1, m_nswbands);
  mem += m_buffer.sw_bnd_flux_dn.totElems();
  m_buffer.sw_bnd_flux_dir = decltype(m_buffer.sw_bnd_flux_dir)("sw_bnd_flux_dir", mem, m_col_chunk_size, m_nlay+1, m_nswbands);
  mem += m_buffer.sw_bnd_flux_dir.totElems();
  m_buffer.sw_bnd_flux_dif = decltype(m_buffer.sw_bnd_flux_dif)("sw_bnd_flux_dif", mem, m_col_chunk_size, m_nlay+1, m_nswbands);
  mem += m_buffer.sw_bnd_flux_dif.totElems();
  // 3d arrays with nlwbands dimension (longwave fluxes by band)
  m_buffer.lw_bnd_flux_up = decltype(m_buffer.lw_bnd_flux_up)("lw_bnd_flux_up", mem, m_col_chunk_size, m_nlay+1, m_nlwbands);
  mem += m_buffer.lw_bnd_flux_up.totElems();
  m_buffer.lw_bnd_flux_dn = decltype(m_buffer.lw_bnd_flux_dn)("lw_bnd_flux_dn", mem, m_col_chunk_size, m_nlay+1, m_nlwbands);
  mem += m_buffer.lw_bnd_flux_dn.totElems();
  // 2d arrays with extra nswbands dimension (surface albedos by band)
  m_buffer.sfc_alb_dir = decltype(m_buffer.sfc_alb_dir)("sfc_alb_dir", mem, m_col_chunk_size, m_nswbands);
  mem += m_buffer.sfc_alb_dir.totElems();
  m_buffer.sfc_alb_dif = decltype(m_buffer.sfc_alb_dif)("sfc_alb_dif", mem, m_col_chunk_size, m_nswbands);
  mem += m_buffer.sfc_alb_dif.totElems();
  // 3d arrays with extra band dimension (aerosol optics by band)
  m_buffer.aero_tau_sw = decltype(m_buffer.aero_tau_sw)("aero_tau_sw", mem, m_col_chunk_size, m_nlay, m_nswbands);
  mem += m_buffer.aero_tau_sw.totElems();
  m_buffer.aero_ssa_sw = decltype(m_buffer.aero_ssa_sw)("aero_ssa_sw", mem, m_col_chunk_size, m_nlay, m_nswbands);
  mem += m_buffer.aero_ssa_sw.totElems();
  m_buffer.aero_g_sw   = decltype(m_buffer.aero_g_sw  )("aero_g_sw"  , mem, m_col_chunk_size, m_nlay, m_nswbands);
  mem += m_buffer.aero_g_sw.totElems();
  m_buffer.aero_tau_lw = decltype(m_buffer.aero_tau_lw)("aero_tau_lw", mem, m_col_chunk_size, m_nlay, m_nlwbands);
  mem += m_buffer.aero_tau_lw.totElems();
  // 3d arrays with extra ngpt dimension (cloud optics by gpoint; primarily for debugging)
  m_buffer.cld_tau_sw_gpt = decltype(m_buffer.cld_tau_sw_gpt)("cld_tau_sw_gpt", mem, m_col_chunk_size, m_nlay, m_nswgpts);
  mem += m_buffer.cld_tau_sw_gpt.totElems();
  m_buffer.cld_tau_lw_gpt = decltype(m_buffer.cld_tau_lw_gpt)("cld_tau_lw_gpt", mem, m_col_chunk_size, m_nlay, m_nlwgpts);
  mem += m_buffer.cld_tau_lw_gpt.totElems();

  size_t used_mem = (reinterpret_cast<Real*>(mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for RRTMGPRadiation.");
} // RRTMGPRadiation::init_buffers

void RRTMGPRadiation::initialize_impl(const RunType /* run_type */) {
  using PC = scream::physics::Constants<Real>;

  // Determine rad timestep, specified as number of atm steps
  m_rad_freq_in_steps = m_params.get<Int>("rad_frequency", 1);

  // Determine orbital year. If orbital_year is negative, use current year
  // from timestamp for orbital year; if positive, use provided orbital year
  // for duration of simulation.
  m_orbital_year = m_params.get<Int>("orbital_year",-9999);
  // Get orbital parameters from yaml file
  m_orbital_eccen = m_params.get<Int>("orbital_eccentricity",-9999);
  m_orbital_obliq = m_params.get<Int>("orbital_obliquity"   ,-9999);
  m_orbital_mvelp = m_params.get<Int>("orbital_mvelp"       ,-9999);

  // Determine whether or not we are using a fixed solar zenith angle (positive value)
  m_fixed_solar_zenith_angle = m_params.get<Real>("Fixed Solar Zenith Angle", -9999);

  // Get prescribed surface values of greenhouse gases
  m_co2vmr     = m_params.get<Real>("co2vmr", 388.717e-6);
  m_n2ovmr     = m_params.get<Real>("n2ovmr", 323.141e-9);
  m_ch4vmr     = m_params.get<Real>("ch4vmr", 1807.851e-9);
  m_f11vmr     = m_params.get<Real>("f11vmr", 768.7644e-12);
  m_f12vmr     = m_params.get<Real>("f12vmr", 531.2820e-12);
  m_n2vmr      = m_params.get<Real>("n2vmr", 0.7906);
  m_covmr      = m_params.get<Real>("covmr", 1.0e-7);

  // Whether or not to do MCICA subcolumn sampling
  m_do_subcol_sampling = m_params.get<bool>("do_subcol_sampling",true);

  // Initialize yakl
  yakl_init();

  // Names of active gases
  auto gas_names_yakl_offset = string1d("gas_names",m_ngas);
  m_gas_mol_weights          = view_1d_real("gas_mol_weights",m_ngas);
  // the lookup function for getting the gas mol weights doesn't work on device
  auto gas_mol_w_host = Kokkos::create_mirror_view(m_gas_mol_weights);
  for (int igas = 0; igas < m_ngas; igas++) {
    const auto& gas_name = m_gas_names[igas];

    /* Note: YAKL starts the index from 1 */
    gas_names_yakl_offset(igas+1)   = gas_name;
    gas_mol_w_host[igas]            = PC::get_gas_mol_weight(gas_name);

  }
  Kokkos::deep_copy(m_gas_mol_weights,gas_mol_w_host);

  // Initialize GasConcs object to pass to RRTMGP initializer;
  std::string coefficients_file_sw = m_params.get<std::string>("rrtmgp_coefficients_file_sw");
  std::string coefficients_file_lw = m_params.get<std::string>("rrtmgp_coefficients_file_lw");
  std::string cloud_optics_file_sw = m_params.get<std::string>("rrtmgp_cloud_optics_file_sw");
  std::string cloud_optics_file_lw = m_params.get<std::string>("rrtmgp_cloud_optics_file_lw");
  m_gas_concs.init(gas_names_yakl_offset,m_col_chunk_size,m_nlay);
  rrtmgp::rrtmgp_initialize(
          m_gas_concs,
          coefficients_file_sw, coefficients_file_lw,
          cloud_optics_file_sw, cloud_optics_file_lw,
          m_atm_logger
  );

  // Set property checks for fields in this process
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("T_mid"),m_grid,100.0, 500.0,false);
}

// =========================================================================================

void RRTMGPRadiation::run_impl (const double dt) {
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using PC = scream::physics::Constants<Real>;
  using CO = scream::ColumnOps<DefaultDevice,Real>;

  // get a host copy of lat/lon
  auto h_lat  = m_lat.get_view<const Real*,Host>();
  auto h_lon  = m_lon.get_view<const Real*,Host>();

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
  auto d_sw_clrsky_flux_up = get_field_out("SW_clrsky_flux_up").get_view<Real**>();
  auto d_sw_clrsky_flux_dn = get_field_out("SW_clrsky_flux_dn").get_view<Real**>();
  auto d_sw_clrsky_flux_dn_dir = get_field_out("SW_clrsky_flux_dn_dir").get_view<Real**>();
  auto d_lw_clrsky_flux_up = get_field_out("LW_clrsky_flux_up").get_view<Real**>();
  auto d_lw_clrsky_flux_dn = get_field_out("LW_clrsky_flux_dn").get_view<Real**>();
  auto d_rad_heating_pdel = get_field_out("rad_heating_pdel").get_view<Real**>();
  auto d_sfc_flux_dir_vis = get_field_out("sfc_flux_dir_vis").get_view<Real*>();
  auto d_sfc_flux_dir_nir = get_field_out("sfc_flux_dir_nir").get_view<Real*>();
  auto d_sfc_flux_dif_vis = get_field_out("sfc_flux_dif_vis").get_view<Real*>();
  auto d_sfc_flux_dif_nir = get_field_out("sfc_flux_dif_nir").get_view<Real*>();
  auto d_sfc_flux_sw_net = get_field_out("sfc_flux_sw_net").get_view<Real*>();
  auto d_sfc_flux_lw_dn  = get_field_out("sfc_flux_lw_dn").get_view<Real*>();
  auto d_cldlow = get_field_out("cldlow").get_view<Real*>();
  auto d_cldmed = get_field_out("cldmed").get_view<Real*>();
  auto d_cldhgh = get_field_out("cldhgh").get_view<Real*>();
  auto d_cldtot = get_field_out("cldtot").get_view<Real*>();

  constexpr auto stebol = PC::stebol;
  const auto nlay = m_nlay;
  const auto nlwbands = m_nlwbands;
  const auto nswbands = m_nswbands;
  const auto nlwgpts = m_nlwgpts;

  // Are we going to update fluxes and heating this step?
  auto ts = timestamp();
  auto update_rad = scream::rrtmgp::radiation_do(m_rad_freq_in_steps, ts.get_num_steps());

  if (update_rad) {
    // On each chunk, we internally "reset" the GasConcs object to subview the concs 3d array
    // with the correct ncol dimension. So let's keep a copy of the original (ref-counted)
    // array, to restore at the end inside the m_gast_concs object.
    auto gas_concs = m_gas_concs.concs;

    // Compute orbital parameters; these are used both for computing
    // the solar zenith angle and also for computing total solar
    // irradiance scaling (tsi_scaling).
    double obliqr, lambm0, mvelpp;
    auto orbital_year = m_orbital_year;
    auto eccen = m_orbital_eccen;
    auto obliq = m_orbital_obliq;
    auto mvelp = m_orbital_mvelp;
    if (eccen >= 0 && obliq >= 0 && mvelp >= 0) {
      // use fixed orbital parameters; to force this, we need to set
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
    auto calday = ts.frac_of_year_in_days() + 1;  // Want day + fraction; calday 1 == Jan 1 0Z
    shr_orb_decl_c2f(calday, eccen, mvelpp, lambm0,
                     obliqr, &delta, &eccf);

    // Loop over each chunk of columns
    for (int ic=0; ic<m_num_col_chunks; ++ic) {
      const int beg  = m_col_chunk_beg[ic];
      const int ncol = m_col_chunk_beg[ic+1] - beg;
      this->log(LogLevel::debug,
                "[RRTMGP::run_impl] Col chunk beg,end: " + std::to_string(beg) + ", " + std::to_string(beg+ncol) + "\n");


      // Create YAKL arrays. RRTMGP expects YAKL arrays with styleFortran, i.e., data has ncol
      // as the fastest index. For this reason we must copy the data.
      auto subview_1d = [&](const real1d v) -> real1d {
        return real1d(v.label(),v.myData,ncol);
      };
      auto subview_2d = [&](const real2d v) -> real2d {
        return real2d(v.label(),v.myData,ncol,v.dimension[1]);
      };
      auto subview_3d = [&](const real3d v) -> real3d {
        return real3d(v.label(),v.myData,ncol,v.dimension[1],v.dimension[2]);
      };

      auto p_lay           = subview_2d(m_buffer.p_lay);
      auto t_lay           = subview_2d(m_buffer.t_lay);
      auto p_lev           = subview_2d(m_buffer.p_lev);
      auto p_del           = subview_2d(m_buffer.p_del);
      auto t_lev           = subview_2d(m_buffer.t_lev);
      auto mu0             = subview_1d(m_buffer.mu0);
      auto sfc_alb_dir     = subview_2d(m_buffer.sfc_alb_dir);
      auto sfc_alb_dif     = subview_2d(m_buffer.sfc_alb_dif);
      auto sfc_alb_dir_vis = subview_1d(m_buffer.sfc_alb_dir_vis);
      auto sfc_alb_dir_nir = subview_1d(m_buffer.sfc_alb_dir_nir);
      auto sfc_alb_dif_vis = subview_1d(m_buffer.sfc_alb_dif_vis);
      auto sfc_alb_dif_nir = subview_1d(m_buffer.sfc_alb_dif_nir);
      auto qc              = subview_2d(m_buffer.qc);
      auto qi              = subview_2d(m_buffer.qi);
      auto cldfrac_tot     = subview_2d(m_buffer.cldfrac_tot);
      auto rel             = subview_2d(m_buffer.eff_radius_qc);
      auto rei             = subview_2d(m_buffer.eff_radius_qi);
      auto sw_flux_up      = subview_2d(m_buffer.sw_flux_up);
      auto sw_flux_dn      = subview_2d(m_buffer.sw_flux_dn);
      auto sw_flux_dn_dir  = subview_2d(m_buffer.sw_flux_dn_dir);
      auto lw_flux_up      = subview_2d(m_buffer.lw_flux_up);
      auto lw_flux_dn      = subview_2d(m_buffer.lw_flux_dn);
      auto sw_clrsky_flux_up      = subview_2d(m_buffer.sw_clrsky_flux_up);
      auto sw_clrsky_flux_dn      = subview_2d(m_buffer.sw_clrsky_flux_dn);
      auto sw_clrsky_flux_dn_dir  = subview_2d(m_buffer.sw_clrsky_flux_dn_dir);
      auto lw_clrsky_flux_up      = subview_2d(m_buffer.lw_clrsky_flux_up);
      auto lw_clrsky_flux_dn      = subview_2d(m_buffer.lw_clrsky_flux_dn);
      auto sw_bnd_flux_up  = subview_3d(m_buffer.sw_bnd_flux_up);
      auto sw_bnd_flux_dn  = subview_3d(m_buffer.sw_bnd_flux_dn);
      auto sw_bnd_flux_dir = subview_3d(m_buffer.sw_bnd_flux_dir);
      auto sw_bnd_flux_dif = subview_3d(m_buffer.sw_bnd_flux_dif);
      auto lw_bnd_flux_up  = subview_3d(m_buffer.lw_bnd_flux_up);
      auto lw_bnd_flux_dn  = subview_3d(m_buffer.lw_bnd_flux_dn);
      auto sfc_flux_dir_vis = subview_1d(m_buffer.sfc_flux_dir_vis);
      auto sfc_flux_dir_nir = subview_1d(m_buffer.sfc_flux_dir_nir);
      auto sfc_flux_dif_vis = subview_1d(m_buffer.sfc_flux_dif_vis);
      auto sfc_flux_dif_nir = subview_1d(m_buffer.sfc_flux_dif_nir);
      auto aero_tau_sw     = subview_3d(m_buffer.aero_tau_sw);
      auto aero_ssa_sw     = subview_3d(m_buffer.aero_ssa_sw);
      auto aero_g_sw       = subview_3d(m_buffer.aero_g_sw);
      auto aero_tau_lw     = subview_3d(m_buffer.aero_tau_lw);
      auto cld_tau_sw_gpt  = subview_3d(m_buffer.cld_tau_sw_gpt);
      auto cld_tau_lw_gpt  = subview_3d(m_buffer.cld_tau_lw_gpt);

      // Set gas concs to "view" only the first ncol columns
      m_gas_concs.ncol = ncol;
      m_gas_concs.concs = subview_3d(gas_concs);

      // Copy data from the FieldManager to the YAKL arrays
      {
        // Determine the cosine zenith angle
        // NOTE: Since we are bridging to F90 arrays this must be done on HOST and then
        //       deep copied to a device view.
        auto d_mu0 = m_buffer.cosine_zenith;
        auto h_mu0 = Kokkos::create_mirror_view(d_mu0);
        if (m_fixed_solar_zenith_angle > 0) {
          for (int i=0; i<ncol; i++) {
            h_mu0(i) = m_fixed_solar_zenith_angle;
          }
        } else {
          // Now use solar declination to calculate zenith angle for all points
          for (int i=0;i<ncol;i++) {
            double lat = h_lat(i+beg)*PC::Pi/180.0;  // Convert lat/lon to radians
            double lon = h_lon(i+beg)*PC::Pi/180.0;
            h_mu0(i) = shr_orb_cosz_c2f(calday, lat, lon, delta, m_rad_freq_in_steps * dt);
          }
        }
        Kokkos::deep_copy(d_mu0,h_mu0);

        // dz and T_int will need to be computed
        view_2d_real d_tint("T_int", ncol, m_nlay+1);
        view_2d_real d_dz  ("dz",    ncol, m_nlay);

        const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
          const int i = team.league_rank();
          const int icol = i+beg;

          // Calculate dz
          const auto pseudo_density = ekat::subview(d_pdel, icol);
          const auto p_mid          = ekat::subview(d_pmid, icol);
          const auto T_mid          = ekat::subview(d_tmid, icol);
          const auto qv             = ekat::subview(d_qv,   icol);
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
          const auto P_mid = ekat::subview(d_pmid, icol);
          const int itop = (P_mid(0) < P_mid(nlay-1)) ? 0 : nlay-1;
          const Real bc_top = T_mid(itop);
          const Real bc_bot = sqrt(sqrt(d_surf_lw_flux_up(icol)/stebol));
          if (itop == 0) {
              CO::compute_interface_values_linear(team, nlay, T_mid, dz, bc_top, bc_bot, T_int);
          } else {
              CO::compute_interface_values_linear(team, nlay, T_mid, dz, bc_bot, bc_top, T_int);
          }
          team.team_barrier();

          mu0(i+1) = d_mu0(i);
          sfc_alb_dir_vis(i+1) = d_sfc_alb_dir_vis(icol);
          sfc_alb_dir_nir(i+1) = d_sfc_alb_dir_nir(icol);
          sfc_alb_dif_vis(i+1) = d_sfc_alb_dif_vis(icol);
          sfc_alb_dif_nir(i+1) = d_sfc_alb_dif_nir(icol);

          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
            p_lay(i+1,k+1)       = d_pmid(icol,k);
            t_lay(i+1,k+1)       = d_tmid(icol,k);
            p_del(i+1,k+1)       = d_pdel(icol,k);
            qc(i+1,k+1)          = d_qc(icol,k);
            qi(i+1,k+1)          = d_qi(icol,k);
            rel(i+1,k+1)         = d_rel(icol,k);
            rei(i+1,k+1)         = d_rei(icol,k);
            p_lev(i+1,k+1)       = d_pint(icol,k);
            t_lev(i+1,k+1)       = d_tint(i,k);
          });

          p_lev(i+1,nlay+1) = d_pint(icol,nlay);
          t_lev(i+1,nlay+1) = d_tint(i,nlay);

          // Note that RRTMGP expects ordering (col,lay,bnd) but the FM keeps things in (col,bnd,lay) order
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nswbands*nlay), [&] (const int&idx) {
              auto b = idx / nlay;
              auto k = idx % nlay;
              aero_tau_sw(i+1,k+1,b+1) = d_aero_tau_sw(icol,b,k);
              aero_ssa_sw(i+1,k+1,b+1) = d_aero_ssa_sw(icol,b,k);
              aero_g_sw  (i+1,k+1,b+1) = d_aero_g_sw  (icol,b,k);
          });
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlwbands*nlay), [&] (const int&idx) {
              auto b = idx / nlay;
              auto k = idx % nlay;
              aero_tau_lw(i+1,k+1,b+1) = d_aero_tau_lw(icol,b,k);
          });
        });
      }
      Kokkos::fence();

      // Populate GasConcs object to pass to RRTMGP driver
      // set_vmr requires the input array size to have the correct size,
      // and the last chunk may have less columns, so create a temp of
      // correct size that uses m_buffer.tmp2d's pointer
      //
      // h2o is taken from qv and requies no initialization here;
      // o3 is computed elsewhere (either read from file or computed by chemistry);
      // n2 and co are set to constants and are not handled by trcmix;
      // the rest are handled by trcmix
      real2d tmp2d = subview_2d(m_buffer.tmp2d);
      const auto gas_mol_weights = m_gas_mol_weights;
      for (int igas = 0; igas < m_ngas; igas++) {
        auto name = m_gas_names[igas];
        auto d_vmr = get_field_out(name + "_volume_mix_ratio").get_view<Real**>();
        if (name == "h2o") {
          // h2o is (wet) mass mixing ratio in FM, otherwise known as "qv", which we've already read in above
          // Convert to vmr
          const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
            const int i = team.league_rank();
            const int icol = i + beg;
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
              d_vmr(icol,k) = PF::calculate_vmr_from_mmr(gas_mol_weights[igas],d_qv(icol,k),d_qv(icol,k));
            });
          });
          Kokkos::fence();
        } else if (name == "o3") {
          // We read o3 in as a vmr already
        } else if (name == "n2") {
          // n2 prescribed as a constant value
          Kokkos::deep_copy(d_vmr, m_params.get<Real>("n2vmr", 0.7906));
        } else if (name == "co") {
          // co prescribed as a constant value
          Kokkos::deep_copy(d_vmr, m_params.get<Real>("covmr", 1.0e-7));
        } else {
          // This gives (dry) mass mixing ratios
          scream::physics::trcmix(
            name, m_lat.get_view<const Real*>(), d_pmid, d_vmr,
            m_co2vmr, m_n2ovmr, m_ch4vmr, m_f11vmr, m_f12vmr
          );
          // Back out volume mixing ratios
          const auto air_mol_weight = PC::MWdry;
          const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_ncol, m_nlay);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
            const int i = team.league_rank();
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
              d_vmr(i,k) = air_mol_weight / gas_mol_weights[igas] * d_vmr(i,k);
            });
          });
        }

        // Copy to YAKL
        const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
          const int i = team.league_rank();
          const int icol = i + beg;
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
            tmp2d(i+1,k+1) = d_vmr(icol,k); // Note that for YAKL arrays i and k start with index 1
          });
        });
        Kokkos::fence();

        // Populate GasConcs object
        m_gas_concs.set_vmr(name, tmp2d);
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
      auto lwp = m_buffer.lwp;
      auto iwp = m_buffer.iwp;
      if (not do_subcol_sampling) {
        const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
          const int i = team.league_rank();
          const int icol = i + beg;
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
            if (d_cldfrac_tot(icol,k) > 0) {
              cldfrac_tot(i+1,k+1) = 1;
            } else {
              cldfrac_tot(i+1,k+1) = 0;
            }
          });
        });
      } else {
        const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
          const int i = team.league_rank();
          const int icol = i + beg;
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
            cldfrac_tot(i+1,k+1) = d_cldfrac_tot(icol,k);
          });
        });
      }
      Kokkos::fence();
      // Compute layer cloud mass (per unit area)
      scream::rrtmgp::mixing_ratio_to_cloud_mass(qc, cldfrac_tot, p_del, lwp);
      scream::rrtmgp::mixing_ratio_to_cloud_mass(qi, cldfrac_tot, p_del, iwp);
      // Convert to g/m2 (needed by RRTMGP)
      {
      const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int i = team.league_rank();
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
          // Note that for YAKL arrays i and k start with index 1
          lwp(i+1,k+1) *= 1e3;
          iwp(i+1,k+1) *= 1e3;
        });
      });
      }
      Kokkos::fence();

      // Compute band-by-band surface_albedos. This is needed since
      // the AD passes broadband albedos, but rrtmgp require band-by-band.
      rrtmgp::compute_band_by_band_surface_albedos(
        ncol, nswbands,
        sfc_alb_dir_vis, sfc_alb_dir_nir,
        sfc_alb_dif_vis, sfc_alb_dif_nir,
        sfc_alb_dir, sfc_alb_dif);

      // Compute cloud optical properties here?

      // Run RRTMGP driver
      rrtmgp::rrtmgp_main(
        ncol, m_nlay,
        p_lay, t_lay, p_lev, t_lev,
        m_gas_concs,
        sfc_alb_dir, sfc_alb_dif, mu0,
        lwp, iwp, rel, rei, cldfrac_tot,
        aero_tau_sw, aero_ssa_sw, aero_g_sw, aero_tau_lw,
        cld_tau_sw_gpt, cld_tau_lw_gpt,
        sw_flux_up       , sw_flux_dn       , sw_flux_dn_dir       , lw_flux_up       , lw_flux_dn,
        sw_clrsky_flux_up, sw_clrsky_flux_dn, sw_clrsky_flux_dn_dir, lw_clrsky_flux_up, lw_clrsky_flux_dn,
        sw_bnd_flux_up   , sw_bnd_flux_dn   , sw_bnd_flux_dir      , lw_bnd_flux_up   , lw_bnd_flux_dn,
        eccf, m_atm_logger
      );

      // Update heating tendency
      auto sw_heating  = m_buffer.sw_heating;
      auto lw_heating  = m_buffer.lw_heating;
      rrtmgp::compute_heating_rate(
        sw_flux_up, sw_flux_dn, p_del, sw_heating
      );
      rrtmgp::compute_heating_rate(
        lw_flux_up, lw_flux_dn, p_del, lw_heating
      );
      {
        const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
          const int idx = team.league_rank();
          const int icol = idx+beg;
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& ilay) {
            // Combine SW and LW heating into a net heating tendency; use d_rad_heating_pdel temporarily
            // Note that for YAKL arrays i and k start with index 1
            d_rad_heating_pdel(icol,ilay) = sw_heating(idx+1,ilay+1) + lw_heating(idx+1,ilay+1);
          });
        });
      }
      Kokkos::fence();

      // Index to surface (bottom of model); used to get surface fluxes below
      const int kbot = nlay+1;

      // Compute diffuse flux as difference between total and direct
      Kokkos::parallel_for(Kokkos::RangePolicy<ExeSpace>(0,nswbands*(nlay+1)*ncol),
                           KOKKOS_LAMBDA (const int idx) {
        // CAREFUL: these are YAKL arrays, with "LayoutLeft". So make the indices stride accordingly, and add 1.
        const int ibnd = (idx / ncol) / (nlay+1) + 1;
        const int ilev = (idx / ncol) % (nlay+1) + 1;
        const int icol =  idx % ncol + 1;
        sw_bnd_flux_dif(icol,ilev,ibnd) = sw_bnd_flux_dn(icol,ilev,ibnd) - sw_bnd_flux_dir(icol,ilev,ibnd);
      });
      // Compute surface fluxes
      rrtmgp::compute_broadband_surface_fluxes(
          ncol, kbot, nswbands,
          sw_bnd_flux_dir, sw_bnd_flux_dif,
          sfc_flux_dir_vis, sfc_flux_dir_nir,
          sfc_flux_dif_vis, sfc_flux_dif_nir
      );

      // Compute diagnostic total cloud area (vertically-projected cloud cover)
      auto cldlow = real1d("cldlow", ncol);
      auto cldmed = real1d("cldmed", ncol);
      auto cldhgh = real1d("cldhgh", ncol);
      auto cldtot = real1d("cldtot", ncol);
      // NOTE: limits for low, mid, and high clouds are mostly taken from EAM F90 source, with the
      // exception that I removed the restriction on low clouds to be above (numerically lower pressures)
      // 1200 hPa, and on high clouds to be below (numerically high pressures) 50 hPa. This probably
      // does not matter in practice, as clouds probably should not be produced above 50 hPa and we
      // should not be encountering surface pressure above 1200 hPa, but in the event that things go off
      // the rails we might want to look at these still.
      rrtmgp::compute_cloud_area(ncol, nlay, nlwgpts, 700e2, std::numeric_limits<Real>::max(), p_lay, cld_tau_lw_gpt, cldlow);
      rrtmgp::compute_cloud_area(ncol, nlay, nlwgpts, 400e2,                            700e2, p_lay, cld_tau_lw_gpt, cldmed);
      rrtmgp::compute_cloud_area(ncol, nlay, nlwgpts,     0,                            400e2, p_lay, cld_tau_lw_gpt, cldhgh);
      rrtmgp::compute_cloud_area(ncol, nlay, nlwgpts,     0, std::numeric_limits<Real>::max(), p_lay, cld_tau_lw_gpt, cldtot);

      // Copy output data back to FieldManager
      const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int i = team.league_rank();
        const int icol = i + beg;
        d_sfc_flux_dir_nir(icol) = sfc_flux_dir_nir(i+1);
        d_sfc_flux_dir_vis(icol) = sfc_flux_dir_vis(i+1);
        d_sfc_flux_dif_nir(icol) = sfc_flux_dif_nir(i+1);
        d_sfc_flux_dif_vis(icol) = sfc_flux_dif_vis(i+1);
        d_sfc_flux_sw_net(icol)  = sw_flux_dn(i+1,kbot) - sw_flux_up(i+1,kbot);
        d_sfc_flux_lw_dn(icol)   = lw_flux_dn(i+1,kbot);
        d_cldlow(icol) = cldlow(i+1);
        d_cldmed(icol) = cldmed(i+1);
        d_cldhgh(icol) = cldhgh(i+1);
        d_cldtot(icol) = cldtot(i+1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay+1), [&] (const int& k) {
          d_sw_flux_up(icol,k)            = sw_flux_up(i+1,k+1);
          d_sw_flux_dn(icol,k)            = sw_flux_dn(i+1,k+1);
          d_sw_flux_dn_dir(icol,k)        = sw_flux_dn_dir(i+1,k+1);
          d_lw_flux_up(icol,k)            = lw_flux_up(i+1,k+1);
          d_lw_flux_dn(icol,k)            = lw_flux_dn(i+1,k+1);
          d_sw_clrsky_flux_up(icol,k)     = sw_clrsky_flux_up(i+1,k+1);
          d_sw_clrsky_flux_dn(icol,k)     = sw_clrsky_flux_dn(i+1,k+1);
          d_sw_clrsky_flux_dn_dir(icol,k) = sw_clrsky_flux_dn_dir(i+1,k+1);
          d_lw_clrsky_flux_up(icol,k)     = lw_clrsky_flux_up(i+1,k+1);
          d_lw_clrsky_flux_dn(icol,k)     = lw_clrsky_flux_dn(i+1,k+1);
        });
      });
    } // loop over chunk

    // Restore the refCounted array.
    m_gas_concs.concs = gas_concs;

  } // update_rad

  // Apply temperature tendency; if we updated radiation this timestep, then d_rad_heating_pdel should
  // contain actual heating rate, not pdel scaled heating rate. Otherwise, if we have NOT updated the
  // radiative heating, then we need to back out the heating from the rad_heating*pdel term that we carry
  // across timesteps to conserve energy.
  const int ncols = m_ncol;
  const int nlays = m_nlay;
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncols, nlays);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlays), [&] (const int& k) {
      if (update_rad) {
        d_tmid(i,k) = d_tmid(i,k) + d_rad_heating_pdel(i,k) * dt;
        d_rad_heating_pdel(i,k) = d_pdel(i,k) * d_rad_heating_pdel(i,k);
      } else {
        auto rad_heat = d_rad_heating_pdel(i,k) / d_pdel(i,k);
        d_tmid(i,k) = d_tmid(i,k) + rad_heat * dt;
      }
    });
  });

  // If necessary, set appropriate boundary fluxes for energy and mass conservation checks.
  // Any boundary fluxes not included in radiation interface are set to 0.
  if (has_column_conservation_check()) {
    auto vapor_flux = get_field_out("vapor_flux").get_view<Real*>();
    auto water_flux = get_field_out("water_flux").get_view<Real*>();
    auto ice_flux   = get_field_out("ice_flux").get_view<Real*>();
    auto heat_flux  = get_field_out("heat_flux").get_view<Real*>();

    const int ncols = m_ncol;
    const int nlays = m_nlay;
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncols, nlays);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int icol = team.league_rank();

      vapor_flux(icol) = 0;
      water_flux(icol) = 0;
      ice_flux(icol)   = 0;

      const auto fsns = d_sw_flux_dn(icol, nlays) - d_sw_flux_up(icol, nlays);
      const auto fsnt = d_sw_flux_dn(icol, 0)     - d_sw_flux_up(icol, 0);
      const auto flns = d_lw_flux_up(icol, nlays) - d_lw_flux_dn(icol, nlays);
      const auto flnt = d_lw_flux_up(icol, 0)     - d_lw_flux_dn(icol, 0);

      heat_flux(icol) = (fsnt - fsns) - (flnt - flns);
    });
  }

}
// =========================================================================================

void RRTMGPRadiation::finalize_impl  () {
  m_gas_concs.reset();
  rrtmgp::rrtmgp_finalize();

  // Finalize YAKL
  yakl_finalize();
}
// =========================================================================================

}  // namespace scream
