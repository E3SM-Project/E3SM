#include "physics/rrtmgp/eamxx_rrtmgp_interface.hpp"
#include "physics/rrtmgp/eamxx_rrtmgp_process_interface.hpp"
#include "physics/rrtmgp/rrtmgp_utils.hpp"
#include "physics/rrtmgp/shr_orb_mod_c2f.hpp"
#include "physics/share/eamxx_trcmix.hpp"

#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_fv_phys_rrtmgp_active_gases_workaround.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"
#include "share/util/eamxx_column_ops.hpp"

#include "ekat/ekat_assert.hpp"

#include "cpp/rrtmgp/mo_gas_concentrations.h"
#ifdef RRTMGP_ENABLE_YAKL
#include "YAKL.h"
#endif

namespace scream {

using KT = KokkosTypes<DefaultDevice>;
using ExeSpace = KT::ExeSpace;
using MemberType = KT::MemberType;

namespace {

struct ConvertToRrtmgpSubview
{
  int beg;
  int ncol;

  template <typename View>
  View subview1d(const View& v) const {
    return View(v, std::make_pair(beg, beg+ncol));
  }

  template <typename View>
  View subview2d_impl(const View& v, const int inner_dim) const {
    return View(v, std::make_pair(beg, beg+ncol), std::make_pair(0, inner_dim));
  }

#ifdef RRTMGP_LAYOUT_LEFT
  template <typename FieldView, typename BufferView>
  BufferView subview2d(const FieldView&, const BufferView& buffer_view, const int inner_dim) const {
    return BufferView(buffer_view, std::make_pair(0, ncol), Kokkos::ALL);
  }
#else
  // Be sure to trim the excess
  // items from the field manager views due to simd packs. If we don't trim, then
  // check_range_k will fail due to looking at unused values. Once rrtmgp can handle
  // packs, this won't be necessary.
  template <typename FieldView, typename BufferView>
  FieldView subview2d(const FieldView& field_view, const BufferView&, const int inner_dim) const {
    return subview2d_impl(field_view, inner_dim);
  }
#endif

  template <typename View>
  View subview3d(const View& v) const {
    // The range assumes these are buffer views, not fields
    return View(v, std::make_pair(0, ncol), Kokkos::ALL, Kokkos::ALL);
  }
};

}

RRTMGPRadiation::
RRTMGPRadiation (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Gather the active gases from the rrtmgp parameter list and assign to the m_gas_names vector.
  const auto& active_gases = m_params.get<std::vector<std::string>>("active_gases");
  for (auto& it : active_gases) {
    // Make sure only unique names are added
    if (not ekat::contains(m_gas_names,it)) {
      m_gas_names.push_back(it);
      if (it=="o3") {
        TraceGasesWorkaround::singleton().add_active_gas(it + "_volume_mix_ratio");
      }
    }
  }

  m_ngas = m_gas_names.size();
}

void RRTMGPRadiation::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {

  using namespace ekat::units;
  using namespace ekat::prefixes;

  // Declare the set of fields used by rrtmgp
  Units m2(m*m,"m2");
  auto nondim = Units::nondimensional();
  auto micron = micro*m;

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  m_ncol = m_grid->get_num_local_dofs();
  m_nlay = m_grid->get_num_vertical_levels();

  if (m_iop_data_manager) {
    // For IOP runs, we need to use the lat/lon from the
    // IOP files instead of the geometry data. Deep copy
    // to device and sync to host since both will be used.
    m_lat = m_grid->get_geometry_data("lat").clone();
    m_lat.deep_copy(m_iop_data_manager->get_params().get<Real>("target_latitude"));
    m_lat.sync_to_host();

    m_lon = m_grid->get_geometry_data("lon").clone();
    m_lon.deep_copy(m_iop_data_manager->get_params().get<Real>("target_longitude"));
    m_lon.sync_to_host();
  } else {
    m_lat = m_grid->get_geometry_data("lat");
    m_lon = m_grid->get_geometry_data("lon");
  }

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
  FieldLayout scalar2d = m_grid->get_2d_scalar_layout();
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(true);
  FieldLayout scalar3d_int = m_grid->get_3d_scalar_layout(false);
  FieldLayout scalar3d_swband = m_grid->get_3d_vector_layout(true,m_nswbands,"swband");
  FieldLayout scalar3d_lwband = m_grid->get_3d_vector_layout(true,m_nlwbands,"lwband");
  FieldLayout scalar3d_swgpts = m_grid->get_3d_vector_layout(true,m_nswgpts,"swgpt");
  FieldLayout scalar3d_lwgpts = m_grid->get_3d_vector_layout(true,m_nlwgpts,"lwgpt");

  // Set required (input) fields here
  add_field<Required>("p_mid" , scalar3d_mid, Pa, grid_name);
  add_field<Required>("p_int", scalar3d_int, Pa, grid_name);
  add_field<Required>("pseudo_density", scalar3d_mid, Pa, grid_name);
  add_field<Required>("sfc_alb_dir_vis", scalar2d, nondim, grid_name);
  add_field<Required>("sfc_alb_dir_nir", scalar2d, nondim, grid_name);
  add_field<Required>("sfc_alb_dif_vis", scalar2d, nondim, grid_name);
  add_field<Required>("sfc_alb_dif_nir", scalar2d, nondim, grid_name);
  add_field<Required>("qc", scalar3d_mid, kg/kg, grid_name);
  add_field<Required>("nc", scalar3d_mid, 1/kg, grid_name);
  add_field<Required>("qi", scalar3d_mid, kg/kg, grid_name);
  add_field<Required>("cldfrac_tot", scalar3d_mid, nondim, grid_name);
  add_field<Required>("eff_radius_qc", scalar3d_mid, micron, grid_name);
  add_field<Required>("eff_radius_qi", scalar3d_mid, micron, grid_name);
  add_field<Required>("qv",scalar3d_mid,kg/kg,grid_name);
  add_field<Required>("surf_lw_flux_up",scalar2d,W/(m*m),grid_name);
  // Set of required gas concentration fields
  for (auto& it : m_gas_names) {
    // Add gas VOLUME mixing ratios (moles of gas / moles of air; what actually gets input to RRTMGP)
    if (it == "o3") {
      // o3 is read from file, or computed by chemistry
      add_field<Required>(it + "_volume_mix_ratio", scalar3d_mid, mol/mol, grid_name);
    } else {
      // the rest are computed by RRTMGP from prescribed surface values
      // NOTE: this may change at some point
      add_field<Computed>(it + "_volume_mix_ratio", scalar3d_mid, mol/mol, grid_name);
    }
  }
  // Required aerosol optical properties from SPA
  m_do_aerosol_rad = m_params.get<bool>("do_aerosol_rad",true);
  if (m_do_aerosol_rad) {
    add_field<Required>("aero_tau_sw", scalar3d_swband, nondim, grid_name);
    add_field<Required>("aero_ssa_sw", scalar3d_swband, nondim, grid_name);
    add_field<Required>("aero_g_sw"  , scalar3d_swband, nondim, grid_name);
    add_field<Required>("aero_tau_lw", scalar3d_lwband, nondim, grid_name);
  }

  // Whether we do extra clean/clear sky calculations
  m_extra_clnclrsky_diag = m_params.get<bool>("extra_clnclrsky_diag", false);
  m_extra_clnsky_diag    = m_params.get<bool>("extra_clnsky_diag", false);

  // Set computed (output) fields
  add_field<Updated >("T_mid"     , scalar3d_mid, K  , grid_name);
  add_field<Computed>("SW_flux_dn", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("SW_flux_up", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("SW_flux_dn_dir", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("LW_flux_up", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("LW_flux_dn", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("SW_clnclrsky_flux_dn", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("SW_clnclrsky_flux_up", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("SW_clnclrsky_flux_dn_dir", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("SW_clrsky_flux_dn", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("SW_clrsky_flux_up", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("SW_clrsky_flux_dn_dir", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("SW_clnsky_flux_dn", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("SW_clnsky_flux_up", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("SW_clnsky_flux_dn_dir", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("LW_clnclrsky_flux_up", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("LW_clnclrsky_flux_dn", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("LW_clrsky_flux_up", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("LW_clrsky_flux_dn", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("LW_clnsky_flux_up", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("LW_clnsky_flux_dn", scalar3d_int, W/m2, grid_name);
  add_field<Computed>("rad_heating_pdel", scalar3d_mid, Pa*K/s, grid_name);
  // Cloud properties added as computed fields for diagnostic purposes
  add_field<Computed>("cldlow"        , scalar2d, nondim, grid_name);
  add_field<Computed>("cldmed"        , scalar2d, nondim, grid_name);
  add_field<Computed>("cldhgh"        , scalar2d, nondim, grid_name);
  add_field<Computed>("cldtot"        , scalar2d, nondim, grid_name);
  // 0.67 micron and 10.5 micron optical depth (needed for COSP)
  add_field<Computed>("dtau067"       , scalar3d_mid, nondim, grid_name);
  add_field<Computed>("dtau105"       , scalar3d_mid, nondim, grid_name);
  add_field<Computed>("sunlit"        , scalar2d    , nondim, grid_name);
  add_field<Computed>("cldfrac_rad"   , scalar3d_mid, nondim, grid_name);
  // Cloud-top diagnostics following AeroCom recommendation
  add_field<Computed>("T_mid_at_cldtop", scalar2d, K, grid_name);
  add_field<Computed>("p_mid_at_cldtop", scalar2d, Pa, grid_name);
  add_field<Computed>("cldfrac_ice_at_cldtop", scalar2d, nondim, grid_name);
  add_field<Computed>("cldfrac_liq_at_cldtop", scalar2d, nondim, grid_name);
  add_field<Computed>("cldfrac_tot_at_cldtop", scalar2d, nondim, grid_name);
  add_field<Computed>("cdnc_at_cldtop", scalar2d, 1 / (m * m * m), grid_name);
  add_field<Computed>("eff_radius_qc_at_cldtop", scalar2d, micron, grid_name);
  add_field<Computed>("eff_radius_qi_at_cldtop", scalar2d, micron, grid_name);

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
  add_field<Computed>("sfc_flux_dir_nir", scalar2d, W/m2, grid_name);
  add_field<Computed>("sfc_flux_dir_vis", scalar2d, W/m2, grid_name);
  add_field<Computed>("sfc_flux_dif_nir", scalar2d, W/m2, grid_name);
  add_field<Computed>("sfc_flux_dif_vis", scalar2d, W/m2, grid_name);
  add_field<Computed>("sfc_flux_sw_net" , scalar2d, W/m2, grid_name);
  add_field<Computed>("sfc_flux_lw_dn"  , scalar2d, W/m2, grid_name);

  // Boundary flux fields for energy and mass conservation checks
  if (has_column_conservation_check()) {
    add_field<Computed>("vapor_flux", scalar2d, kg/m2/s, grid_name);
    add_field<Computed>("water_flux", scalar2d, m/s,     grid_name);
    add_field<Computed>("ice_flux",   scalar2d, m/s,     grid_name);
    add_field<Computed>("heat_flux",  scalar2d, W/m2,    grid_name);
  }

  // Load bands bounds from coefficients files and compute the band centerpoint.
  // Store both in the grid (if not already present)
  const auto cm = centi*m;
  for (std::string prefix : {"sw", "lw"} ) {
    int nbands = prefix == "sw" ? m_nswbands : m_nlwbands;

    if (not m_grid->has_geometry_data(prefix + "band_bounds")) {
      using namespace ShortFieldTagsNames;

      // NOTE: use append, so we get builtin name for (CMP,2) dim, without hard-coding it here
      FieldLayout layout({CMP},{nbands},{prefix+"band"});
      layout.append_dim(CMP,2);
      Field bounds (FieldIdentifier(prefix + "band_bounds", layout, 1/cm, grid_name));
      bounds.allocate_view();

      std::string fname = m_params.get<std::string>("rrtmgp_coefficients_file_" + prefix);
      scorpio::register_file(fname,scorpio::FileMode::Read);
      scorpio::read_var(fname,"bnd_limits_wavenumber",bounds.get_view<Real**,Host>().data());
      scorpio::release_file(fname);

      bounds.sync_to_dev();
      m_grid->set_geometry_data(bounds);
    }

    // If no bounds were in the grid, the bands centerpoint likely wouldn't either. Still, let's check...
    if (not m_grid->has_geometry_data(prefix + "bands")) {
      auto bounds = m_grid->get_geometry_data(prefix + "band_bounds");
      auto bounds_h = bounds.get_view<const Real**,Host>();

      auto bands = bounds.subfield(1,0).clone(prefix + "band");
      auto bands_h = bands.get_view<Real*,Host>();
      for (int i=0; i<nbands; ++i) {
        bands_h(i) = (bounds_h(i,0) + bounds_h(i,1)) / 2;
      }
      bands.sync_to_dev();
      m_grid->set_geometry_data(bands);
    }
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
    Buffer::num_3d_nlay_nlwgpts*m_col_chunk_size*(m_nlay)*m_nlwgpts *
#if defined(RRTMGP_ENABLE_YAKL) && defined(RRTMGP_ENABLE_KOKKOS)
    2;
#else
    1;
#endif

  return interface_request * sizeof(Real);
} // RRTMGPRadiation::requested_buffer_size
// =========================================================================================

void RRTMGPRadiation::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

#ifdef RRTMGP_ENABLE_YAKL
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
  m_buffer.z_del = decltype(m_buffer.z_del)("z_del", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.z_del.totElems();
  m_buffer.p_del = decltype(m_buffer.p_del)("p_del", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.p_del.totElems();
  m_buffer.qc = decltype(m_buffer.qc)("qc", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.qc.totElems();
  m_buffer.nc = decltype(m_buffer.nc)("nc", mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.nc.totElems();
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
  m_buffer.p_lev = decltype(m_buffer.p_lev)("p_lev", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.p_lev.totElems();
  m_buffer.t_lev = decltype(m_buffer.t_lev)("t_lev", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.t_lev.totElems();
  m_buffer.d_tint = decltype(m_buffer.d_tint)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.d_tint.size();
  m_buffer.d_dz  = decltype(m_buffer.d_dz )(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.d_dz.size();
  // 3d arrays
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
  m_buffer.sw_clnclrsky_flux_up = decltype(m_buffer.sw_clnclrsky_flux_up)("sw_clnclrsky_flux_up", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clnclrsky_flux_up.totElems();
  m_buffer.sw_clnclrsky_flux_dn = decltype(m_buffer.sw_clnclrsky_flux_dn)("sw_clnclrsky_flux_dn", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clnclrsky_flux_dn.totElems();
  m_buffer.sw_clnclrsky_flux_dn_dir = decltype(m_buffer.sw_clnclrsky_flux_dn_dir)("sw_clnclrsky_flux_dn_dir", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clnclrsky_flux_dn_dir.totElems();
  m_buffer.sw_clrsky_flux_up = decltype(m_buffer.sw_clrsky_flux_up)("sw_clrsky_flux_up", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clrsky_flux_up.totElems();
  m_buffer.sw_clrsky_flux_dn = decltype(m_buffer.sw_clrsky_flux_dn)("sw_clrsky_flux_dn", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clrsky_flux_dn.totElems();
  m_buffer.sw_clrsky_flux_dn_dir = decltype(m_buffer.sw_clrsky_flux_dn_dir)("sw_clrsky_flux_dn_dir", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clrsky_flux_dn_dir.totElems();
  m_buffer.sw_clnsky_flux_up = decltype(m_buffer.sw_clnsky_flux_up)("sw_clnsky_flux_up", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clnsky_flux_up.totElems();
  m_buffer.sw_clnsky_flux_dn = decltype(m_buffer.sw_clnsky_flux_dn)("sw_clnsky_flux_dn", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clnsky_flux_dn.totElems();
  m_buffer.sw_clnsky_flux_dn_dir = decltype(m_buffer.sw_clnsky_flux_dn_dir)("sw_clnsky_flux_dn_dir", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clnsky_flux_dn_dir.totElems();
  m_buffer.lw_clnclrsky_flux_up = decltype(m_buffer.lw_clnclrsky_flux_up)("lw_clnclrsky_flux_up", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clnclrsky_flux_up.totElems();
  m_buffer.lw_clnclrsky_flux_dn = decltype(m_buffer.lw_clnclrsky_flux_dn)("lw_clnclrsky_flux_dn", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clnclrsky_flux_dn.totElems();
  m_buffer.lw_clrsky_flux_up = decltype(m_buffer.lw_clrsky_flux_up)("lw_clrsky_flux_up", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clrsky_flux_up.totElems();
  m_buffer.lw_clrsky_flux_dn = decltype(m_buffer.lw_clrsky_flux_dn)("lw_clrsky_flux_dn", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clrsky_flux_dn.totElems();
  m_buffer.lw_clnsky_flux_up = decltype(m_buffer.lw_clnsky_flux_up)("lw_clnsky_flux_up", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clnsky_flux_up.totElems();
  m_buffer.lw_clnsky_flux_dn = decltype(m_buffer.lw_clnsky_flux_dn)("lw_clnsky_flux_dn", mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clnsky_flux_dn.totElems();
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
  m_buffer.cld_tau_sw_bnd = decltype(m_buffer.cld_tau_sw_bnd)("cld_tau_sw_bnd", mem, m_col_chunk_size, m_nlay, m_nswbands);
  mem += m_buffer.cld_tau_sw_bnd.totElems();
  m_buffer.cld_tau_lw_bnd = decltype(m_buffer.cld_tau_lw_bnd)("cld_tau_lw_bnd", mem, m_col_chunk_size, m_nlay, m_nlwbands);
  mem += m_buffer.cld_tau_lw_bnd.totElems();
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
  // 1d arrays
  m_buffer.mu0_k = decltype(m_buffer.mu0_k)(mem, m_col_chunk_size);
  mem += m_buffer.mu0_k.size();
  m_buffer.sfc_alb_dir_vis_k = decltype(m_buffer.sfc_alb_dir_vis_k)(mem, m_col_chunk_size);
  mem += m_buffer.sfc_alb_dir_vis_k.size();
  m_buffer.sfc_alb_dir_nir_k = decltype(m_buffer.sfc_alb_dir_nir_k)(mem, m_col_chunk_size);
  mem += m_buffer.sfc_alb_dir_nir_k.size();
  m_buffer.sfc_alb_dif_vis_k = decltype(m_buffer.sfc_alb_dif_vis_k)(mem, m_col_chunk_size);
  mem += m_buffer.sfc_alb_dif_vis_k.size();
  m_buffer.sfc_alb_dif_nir_k = decltype(m_buffer.sfc_alb_dif_nir_k)(mem, m_col_chunk_size);
  mem += m_buffer.sfc_alb_dif_nir_k.size();
  m_buffer.sfc_flux_dir_vis_k = decltype(m_buffer.sfc_flux_dir_vis_k)(mem, m_col_chunk_size);
  mem += m_buffer.sfc_flux_dir_vis_k.size();
  m_buffer.sfc_flux_dir_nir_k = decltype(m_buffer.sfc_flux_dir_nir_k)(mem, m_col_chunk_size);
  mem += m_buffer.sfc_flux_dir_nir_k.size();
  m_buffer.sfc_flux_dif_vis_k = decltype(m_buffer.sfc_flux_dif_vis_k)(mem, m_col_chunk_size);
  mem += m_buffer.sfc_flux_dif_vis_k.size();
  m_buffer.sfc_flux_dif_nir_k = decltype(m_buffer.sfc_flux_dif_nir_k)(mem, m_col_chunk_size);
  mem += m_buffer.sfc_flux_dif_nir_k.size();
  m_buffer.cosine_zenith = decltype(m_buffer.cosine_zenith)(mem, m_col_chunk_size);
  mem += m_buffer.cosine_zenith.size();

  // 2d arrays
  m_buffer.p_lay_k = decltype(m_buffer.p_lay_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.p_lay_k.size();
  m_buffer.t_lay_k = decltype(m_buffer.t_lay_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.t_lay_k.size();
  m_buffer.z_del_k = decltype(m_buffer.z_del_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.z_del_k.size();
  m_buffer.p_del_k = decltype(m_buffer.p_del_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.p_del_k.size();
  m_buffer.qc_k = decltype(m_buffer.qc_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.qc_k.size();
  m_buffer.nc_k = decltype(m_buffer.nc_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.nc_k.size();
  m_buffer.qi_k = decltype(m_buffer.qi_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.qi_k.size();
  m_buffer.cldfrac_tot_k = decltype(m_buffer.cldfrac_tot_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.cldfrac_tot_k.size();
  m_buffer.eff_radius_qc_k = decltype(m_buffer.eff_radius_qc_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.eff_radius_qc_k.size();
  m_buffer.eff_radius_qi_k = decltype(m_buffer.eff_radius_qi_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.eff_radius_qi_k.size();
  m_buffer.tmp2d_k = decltype(m_buffer.tmp2d_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.tmp2d_k.size();
  m_buffer.lwp_k = decltype(m_buffer.lwp_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.lwp_k.size();
  m_buffer.iwp_k = decltype(m_buffer.iwp_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.iwp_k.size();
  m_buffer.sw_heating_k = decltype(m_buffer.sw_heating_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.sw_heating_k.size();
  m_buffer.lw_heating_k = decltype(m_buffer.lw_heating_k)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.lw_heating_k.size();
  m_buffer.p_lev_k = decltype(m_buffer.p_lev_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.p_lev_k.size();
  m_buffer.t_lev_k = decltype(m_buffer.t_lev_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.t_lev_k.size();
  m_buffer.d_tint = decltype(m_buffer.d_tint)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.d_tint.size();
  m_buffer.d_dz  = decltype(m_buffer.d_dz)(mem, m_col_chunk_size, m_nlay);
  mem += m_buffer.d_dz.size();
  // 3d arrays
  m_buffer.sw_flux_up_k = decltype(m_buffer.sw_flux_up_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_flux_up_k.size();
  m_buffer.sw_flux_dn_k = decltype(m_buffer.sw_flux_dn_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_flux_dn_k.size();
  m_buffer.sw_flux_dn_dir_k = decltype(m_buffer.sw_flux_dn_dir_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_flux_dn_dir_k.size();
  m_buffer.lw_flux_up_k = decltype(m_buffer.lw_flux_up_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_flux_up_k.size();
  m_buffer.lw_flux_dn_k = decltype(m_buffer.lw_flux_dn_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_flux_dn_k.size();
  m_buffer.sw_clnclrsky_flux_up_k = decltype(m_buffer.sw_clnclrsky_flux_up_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clnclrsky_flux_up_k.size();
  m_buffer.sw_clnclrsky_flux_dn_k = decltype(m_buffer.sw_clnclrsky_flux_dn_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clnclrsky_flux_dn_k.size();
  m_buffer.sw_clnclrsky_flux_dn_dir_k = decltype(m_buffer.sw_clnclrsky_flux_dn_dir_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clnclrsky_flux_dn_dir_k.size();
  m_buffer.sw_clrsky_flux_up_k = decltype(m_buffer.sw_clrsky_flux_up_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clrsky_flux_up_k.size();
  m_buffer.sw_clrsky_flux_dn_k = decltype(m_buffer.sw_clrsky_flux_dn_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clrsky_flux_dn_k.size();
  m_buffer.sw_clrsky_flux_dn_dir_k = decltype(m_buffer.sw_clrsky_flux_dn_dir_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clrsky_flux_dn_dir_k.size();
  m_buffer.sw_clnsky_flux_up_k = decltype(m_buffer.sw_clnsky_flux_up_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clnsky_flux_up_k.size();
  m_buffer.sw_clnsky_flux_dn_k = decltype(m_buffer.sw_clnsky_flux_dn_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clnsky_flux_dn_k.size();
  m_buffer.sw_clnsky_flux_dn_dir_k = decltype(m_buffer.sw_clnsky_flux_dn_dir_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.sw_clnsky_flux_dn_dir_k.size();
  m_buffer.lw_clnclrsky_flux_up_k = decltype(m_buffer.lw_clnclrsky_flux_up_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clnclrsky_flux_up_k.size();
  m_buffer.lw_clnclrsky_flux_dn_k = decltype(m_buffer.lw_clnclrsky_flux_dn_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clnclrsky_flux_dn_k.size();
  m_buffer.lw_clrsky_flux_up_k = decltype(m_buffer.lw_clrsky_flux_up_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clrsky_flux_up_k.size();
  m_buffer.lw_clrsky_flux_dn_k = decltype(m_buffer.lw_clrsky_flux_dn_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clrsky_flux_dn_k.size();
  m_buffer.lw_clnsky_flux_up_k = decltype(m_buffer.lw_clnsky_flux_up_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clnsky_flux_up_k.size();
  m_buffer.lw_clnsky_flux_dn_k = decltype(m_buffer.lw_clnsky_flux_dn_k)(mem, m_col_chunk_size, m_nlay+1);
  mem += m_buffer.lw_clnsky_flux_dn_k.size();
  // 3d arrays with nswbands dimension (shortwave fluxes by band)
  m_buffer.sw_bnd_flux_up_k = decltype(m_buffer.sw_bnd_flux_up_k)(mem, m_col_chunk_size, m_nlay+1, m_nswbands);
  mem += m_buffer.sw_bnd_flux_up_k.size();
  m_buffer.sw_bnd_flux_dn_k = decltype(m_buffer.sw_bnd_flux_dn_k)(mem, m_col_chunk_size, m_nlay+1, m_nswbands);
  mem += m_buffer.sw_bnd_flux_dn_k.size();
  m_buffer.sw_bnd_flux_dir_k = decltype(m_buffer.sw_bnd_flux_dir_k)(mem, m_col_chunk_size, m_nlay+1, m_nswbands);
  mem += m_buffer.sw_bnd_flux_dir_k.size();
  m_buffer.sw_bnd_flux_dif_k = decltype(m_buffer.sw_bnd_flux_dif_k)(mem, m_col_chunk_size, m_nlay+1, m_nswbands);
  mem += m_buffer.sw_bnd_flux_dif_k.size();
  // 3d arrays with nlwbands dimension (longwave fluxes by band)
  m_buffer.lw_bnd_flux_up_k = decltype(m_buffer.lw_bnd_flux_up_k)(mem, m_col_chunk_size, m_nlay+1, m_nlwbands);
  mem += m_buffer.lw_bnd_flux_up_k.size();
  m_buffer.lw_bnd_flux_dn_k = decltype(m_buffer.lw_bnd_flux_dn_k)(mem, m_col_chunk_size, m_nlay+1, m_nlwbands);
  mem += m_buffer.lw_bnd_flux_dn_k.size();
  // 2d arrays with extra nswbands dimension (surface albedos by band)
  m_buffer.sfc_alb_dir_k = decltype(m_buffer.sfc_alb_dir_k)(mem, m_col_chunk_size, m_nswbands);
  mem += m_buffer.sfc_alb_dir_k.size();
  m_buffer.sfc_alb_dif_k = decltype(m_buffer.sfc_alb_dif_k)(mem, m_col_chunk_size, m_nswbands);
  mem += m_buffer.sfc_alb_dif_k.size();
  // 3d arrays with extra band dimension (aerosol optics by band)
  m_buffer.aero_tau_sw_k = decltype(m_buffer.aero_tau_sw_k)(mem, m_col_chunk_size, m_nlay, m_nswbands);
  mem += m_buffer.aero_tau_sw_k.size();
  m_buffer.aero_ssa_sw_k = decltype(m_buffer.aero_ssa_sw_k)(mem, m_col_chunk_size, m_nlay, m_nswbands);
  mem += m_buffer.aero_ssa_sw_k.size();
  m_buffer.aero_g_sw_k   = decltype(m_buffer.aero_g_sw_k  )(mem, m_col_chunk_size, m_nlay, m_nswbands);
  mem += m_buffer.aero_g_sw_k.size();
  m_buffer.aero_tau_lw_k = decltype(m_buffer.aero_tau_lw_k)(mem, m_col_chunk_size, m_nlay, m_nlwbands);
  mem += m_buffer.aero_tau_lw_k.size();
  // 3d arrays with extra ngpt dimension (cloud optics by gpoint; primarily for debugging)
  m_buffer.cld_tau_sw_gpt_k = decltype(m_buffer.cld_tau_sw_gpt_k)(mem, m_col_chunk_size, m_nlay, m_nswgpts);
  mem += m_buffer.cld_tau_sw_gpt_k.size();
  m_buffer.cld_tau_lw_gpt_k = decltype(m_buffer.cld_tau_lw_gpt_k)(mem, m_col_chunk_size, m_nlay, m_nlwgpts);
  mem += m_buffer.cld_tau_lw_gpt_k.size();
  m_buffer.cld_tau_sw_bnd_k = decltype(m_buffer.cld_tau_sw_bnd_k)(mem, m_col_chunk_size, m_nlay, m_nswbands);
  mem += m_buffer.cld_tau_sw_bnd_k.size();
  m_buffer.cld_tau_lw_bnd_k = decltype(m_buffer.cld_tau_lw_bnd_k)(mem, m_col_chunk_size, m_nlay, m_nlwbands);
  mem += m_buffer.cld_tau_lw_bnd_k.size();
#endif

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
  m_orbital_eccen = m_params.get<double>("orbital_eccentricity",-9999);
  m_orbital_obliq = m_params.get<double>("orbital_obliquity"   ,-9999);
  m_orbital_mvelp = m_params.get<double>("orbital_mvelp"       ,-9999);

  // Value for prescribing an invariant solar constant (i.e. total solar irradiance at
  // TOA).  Used for idealized experiments such as RCE. Disabled when value is less than 0.
  m_fixed_total_solar_irradiance = m_params.get<double>("fixed_total_solar_irradiance", -9999);

  // Determine whether or not we are using a fixed solar zenith angle (positive value)
  m_fixed_solar_zenith_angle = m_params.get<double>("Fixed Solar Zenith Angle", -9999);

  // Get prescribed surface values of greenhouse gases
  m_co2vmr     = m_params.get<double>("co2vmr", 388.717e-6);
  m_n2ovmr     = m_params.get<double>("n2ovmr", 323.141e-9);
  m_ch4vmr     = m_params.get<double>("ch4vmr", 1807.851e-9);
  m_f11vmr     = m_params.get<double>("f11vmr", 768.7644e-12);
  m_f12vmr     = m_params.get<double>("f12vmr", 531.2820e-12);
  m_n2vmr      = m_params.get<double>("n2vmr", 0.7906);
  m_covmr      = m_params.get<double>("covmr", 1.0e-7);

  // Whether or not to do MCICA subcolumn sampling
  m_do_subcol_sampling = m_params.get<bool>("do_subcol_sampling",true);

  // Initialize yakl
  init_kls();

  // Names of active gases
  auto gas_names_yakl_offset = string1dv(m_ngas);
  m_gas_mol_weights          = real1dk("gas_mol_weights",m_ngas);
  // the lookup function for getting the gas mol weights doesn't work on device
  auto gas_mol_w_host = Kokkos::create_mirror_view(m_gas_mol_weights);
  for (int igas = 0; igas < m_ngas; igas++) {
    const auto& gas_name = m_gas_names[igas];

    gas_names_yakl_offset[igas] = gas_name;
    gas_mol_w_host[igas]        = PC::get_gas_mol_weight(gas_name);

  }
  Kokkos::deep_copy(m_gas_mol_weights,gas_mol_w_host);

  // Initialize GasConcs object to pass to RRTMGP initializer;
  std::string coefficients_file_sw = m_params.get<std::string>("rrtmgp_coefficients_file_sw");
  std::string coefficients_file_lw = m_params.get<std::string>("rrtmgp_coefficients_file_lw");
  std::string cloud_optics_file_sw = m_params.get<std::string>("rrtmgp_cloud_optics_file_sw");
  std::string cloud_optics_file_lw = m_params.get<std::string>("rrtmgp_cloud_optics_file_lw");
#ifdef RRTMGP_ENABLE_YAKL
  m_gas_concs.init(gas_names_yakl_offset,m_col_chunk_size,m_nlay);
  rrtmgp::rrtmgp_initialize(
          m_gas_concs,
          coefficients_file_sw, coefficients_file_lw,
          cloud_optics_file_sw, cloud_optics_file_lw,
          m_atm_logger
  );
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
  const double multiplier = m_params.get<double>("pool_size_multiplier", 1.0);

  m_gas_concs_k.init(gas_names_yakl_offset,m_col_chunk_size,m_nlay);
  interface_t::rrtmgp_initialize(
          m_gas_concs_k,
          coefficients_file_sw, coefficients_file_lw,
          cloud_optics_file_sw, cloud_optics_file_lw,
          m_atm_logger,
          multiplier
  );
  VALIDATE_KOKKOS(m_gas_concs, m_gas_concs_k);
  VALIDATE_KOKKOS(rrtmgp::k_dist_sw, *interface_t::k_dist_sw_k);
  VALIDATE_KOKKOS(rrtmgp::k_dist_lw, *interface_t::k_dist_lw_k);
  VALIDATE_KOKKOS(rrtmgp::cloud_optics_sw, *interface_t::cloud_optics_sw_k);
  VALIDATE_KOKKOS(rrtmgp::cloud_optics_lw, *interface_t::cloud_optics_lw_k);
#endif

  // Set property checks for fields in this process
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("T_mid"),m_grid,100.0, 500.0,false);

  // VMR of n2 and co is currently prescribed as a constant value, read from file
  if (has_computed_field("n2_volume_mix_ratio",m_grid->name())) {
    auto n2_vmr = get_field_out("n2_volume_mix_ratio").get_view<Real**>();
    Kokkos::deep_copy(n2_vmr, m_params.get<double>("n2vmr", 0.7906));
  }
  if (has_computed_field("co_volume_mix_ratio",m_grid->name())) {
    auto co_vmr = get_field_out("co_volume_mix_ratio").get_view<Real**>();
    Kokkos::deep_copy(co_vmr, m_params.get<double>("covmr", 1.0e-7));
  }
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
  auto d_nc = get_field_in("nc").get_view<const Real**>();
  auto d_qi = get_field_in("qi").get_view<const Real**>();
  auto d_cldfrac_tot = get_field_in("cldfrac_tot").get_view<const Real**>();
  auto d_rel = get_field_in("eff_radius_qc").get_view<const Real**>();
  auto d_rei = get_field_in("eff_radius_qi").get_view<const Real**>();
  auto d_surf_lw_flux_up = get_field_in("surf_lw_flux_up").get_view<const Real*>();
  // Output fields
  auto d_tmid = get_field_out("T_mid").get_view<Real**>();
  auto d_cldfrac_rad = get_field_out("cldfrac_rad").get_view<Real**>();

  // Aerosol optics only exist if m_do_aerosol_rad is true, so declare views and copy from FM if so
  using view_3d = Field::view_dev_t<const Real***>;
  view_3d d_aero_tau_sw;
  view_3d d_aero_ssa_sw;
  view_3d d_aero_g_sw;
  view_3d d_aero_tau_lw;
  if (m_do_aerosol_rad) {
    d_aero_tau_sw = get_field_in("aero_tau_sw").get_view<const Real***>();
    d_aero_ssa_sw = get_field_in("aero_ssa_sw").get_view<const Real***>();
    d_aero_g_sw   = get_field_in("aero_g_sw"  ).get_view<const Real***>();
    d_aero_tau_lw = get_field_in("aero_tau_lw").get_view<const Real***>();
  }
  auto d_sw_flux_up = get_field_out("SW_flux_up").get_view<Real**>();
  auto d_sw_flux_dn = get_field_out("SW_flux_dn").get_view<Real**>();
  auto d_sw_flux_dn_dir = get_field_out("SW_flux_dn_dir").get_view<Real**>();
  auto d_lw_flux_up = get_field_out("LW_flux_up").get_view<Real**>();
  auto d_lw_flux_dn = get_field_out("LW_flux_dn").get_view<Real**>();
  auto d_sw_clnclrsky_flux_up = get_field_out("SW_clnclrsky_flux_up").get_view<Real**>();
  auto d_sw_clnclrsky_flux_dn = get_field_out("SW_clnclrsky_flux_dn").get_view<Real**>();
  auto d_sw_clnclrsky_flux_dn_dir = get_field_out("SW_clnclrsky_flux_dn_dir").get_view<Real**>();
  auto d_sw_clrsky_flux_up = get_field_out("SW_clrsky_flux_up").get_view<Real**>();
  auto d_sw_clrsky_flux_dn = get_field_out("SW_clrsky_flux_dn").get_view<Real**>();
  auto d_sw_clrsky_flux_dn_dir = get_field_out("SW_clrsky_flux_dn_dir").get_view<Real**>();
  auto d_sw_clnsky_flux_up = get_field_out("SW_clnsky_flux_up").get_view<Real**>();
  auto d_sw_clnsky_flux_dn = get_field_out("SW_clnsky_flux_dn").get_view<Real**>();
  auto d_sw_clnsky_flux_dn_dir = get_field_out("SW_clnsky_flux_dn_dir").get_view<Real**>();
  auto d_lw_clnclrsky_flux_up = get_field_out("LW_clnclrsky_flux_up").get_view<Real**>();
  auto d_lw_clnclrsky_flux_dn = get_field_out("LW_clnclrsky_flux_dn").get_view<Real**>();
  auto d_lw_clrsky_flux_up = get_field_out("LW_clrsky_flux_up").get_view<Real**>();
  auto d_lw_clrsky_flux_dn = get_field_out("LW_clrsky_flux_dn").get_view<Real**>();
  auto d_lw_clnsky_flux_up = get_field_out("LW_clnsky_flux_up").get_view<Real**>();
  auto d_lw_clnsky_flux_dn = get_field_out("LW_clnsky_flux_dn").get_view<Real**>();
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
  // Outputs for COSP
  auto d_dtau067 = get_field_out("dtau067").get_view<Real**>();
  auto d_dtau105 = get_field_out("dtau105").get_view<Real**>();
  auto d_sunlit = get_field_out("sunlit").get_view<Real*>();

  Kokkos::deep_copy(d_dtau067,0.0);
  Kokkos::deep_copy(d_dtau105,0.0);
  // Outputs for AeroCom cloud-top diagnostics
  auto d_T_mid_at_cldtop = get_field_out("T_mid_at_cldtop").get_view<Real *>();
  auto d_p_mid_at_cldtop = get_field_out("p_mid_at_cldtop").get_view<Real *>();
  auto d_cldfrac_ice_at_cldtop =
      get_field_out("cldfrac_ice_at_cldtop").get_view<Real *>();
  auto d_cldfrac_liq_at_cldtop =
      get_field_out("cldfrac_liq_at_cldtop").get_view<Real *>();
  auto d_cldfrac_tot_at_cldtop =
      get_field_out("cldfrac_tot_at_cldtop").get_view<Real *>();
  auto d_cdnc_at_cldtop = get_field_out("cdnc_at_cldtop").get_view<Real *>();
  auto d_eff_radius_qc_at_cldtop =
      get_field_out("eff_radius_qc_at_cldtop").get_view<Real *>();
  auto d_eff_radius_qi_at_cldtop =
      get_field_out("eff_radius_qi_at_cldtop").get_view<Real *>();

  constexpr auto stebol = PC::stebol;
  const auto nlay = m_nlay;
  const auto nlwbands = m_nlwbands;
  const auto nswbands = m_nswbands;
  const auto nlwgpts = m_nlwgpts;
  const auto do_aerosol_rad = m_do_aerosol_rad;

  // Are we going to update fluxes and heating this step?
  auto ts = start_of_step_ts();
  auto update_rad = scream::rrtmgp::radiation_do(m_rad_freq_in_steps, ts.get_num_steps());

  if (update_rad) {
    // On each chunk, we internally "reset" the GasConcs object to subview the concs 3d array
    // with the correct ncol dimension. So let's keep a copy of the original (ref-counted)
    // array, to restore at the end inside the m_gast_concs object.
#ifdef RRTMGP_ENABLE_YAKL
    auto gas_concs = m_gas_concs.concs;
    auto orig_ncol = m_gas_concs.ncol;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    auto gas_concs_k = m_gas_concs_k.concs;
    auto orig_ncol_k = m_gas_concs_k.ncol;
#endif

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

    // Overwrite eccf if using a fixed solar constant.
    auto fixed_total_solar_irradiance = m_fixed_total_solar_irradiance;
    if (fixed_total_solar_irradiance >= 0){
       eccf = fixed_total_solar_irradiance/1360.9;
    }

    // Precompute VMR for all gases, on all cols, before starting the chunks loop
    //
    // h2o is taken from qv
    // o3 is computed elsewhere (either read from file or computed by chemistry);
    // n2 and co are set to constants and are not handled by trcmix;
    // the rest are handled by trcmix
    const auto gas_mol_weights = m_gas_mol_weights;
    for (int igas = 0; igas < m_ngas; igas++) {
      auto name = m_gas_names[igas];

      // We read o3 in as a vmr already. Also, n2 and co are currently set
      // as a constant value, read from file during init. Skip these.
      if (name=="o3" or name == "n2" or name == "co") continue;

      auto d_vmr = get_field_out(name + "_volume_mix_ratio").get_view<Real**>();
      if (name == "h2o") {
        // h2o is (wet) mass mixing ratio in FM, otherwise known as "qv", which we've already read in above
        // Convert to vmr
        const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_ncol, m_nlay);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
          const int icol = team.league_rank();
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
            d_vmr(icol,k) = PF::calculate_vmr_from_mmr(gas_mol_weights[igas],d_qv(icol,k),d_qv(icol,k));
          });
        });
        Kokkos::fence();
      } else {
        // This gives (dry) mass mixing ratios
        scream::physics::trcmix(
          name, m_nlay, m_lat.get_view<const Real*>(), d_pmid, d_vmr,
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
    }

    // Loop over each chunk of columns
    for (int ic=0; ic<m_num_col_chunks; ++ic) {
      const int beg  = m_col_chunk_beg[ic];
      const int ncol = m_col_chunk_beg[ic+1] - beg;
      this->log(LogLevel::debug,
                "[RRTMGP::run_impl] Col chunk beg,end: " + std::to_string(beg) + ", " + std::to_string(beg+ncol) + "\n");

      // d_tint and d_dz are used in eamxx calls and therefore
      // must be layout right
      ulrreal2dk d_tint = ulrreal2dk(m_buffer.d_tint.data(), m_col_chunk_size, m_nlay+1);
      ulrreal2dk d_dz   = ulrreal2dk(m_buffer.d_dz.data(), m_col_chunk_size, m_nlay);
      auto d_mu0 = m_buffer.cosine_zenith;
#ifdef RRTMGP_ENABLE_YAKL
      TIMED_INLINE_KERNEL(init_views,
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
      auto z_del           = subview_2d(m_buffer.z_del);
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
      auto nc              = subview_2d(m_buffer.nc);
      auto qi              = subview_2d(m_buffer.qi);
      auto cldfrac_tot     = subview_2d(m_buffer.cldfrac_tot);
      auto rel             = subview_2d(m_buffer.eff_radius_qc);
      auto rei             = subview_2d(m_buffer.eff_radius_qi);
      auto sw_flux_up      = subview_2d(m_buffer.sw_flux_up);
      auto sw_flux_dn      = subview_2d(m_buffer.sw_flux_dn);
      auto sw_flux_dn_dir  = subview_2d(m_buffer.sw_flux_dn_dir);
      auto lw_flux_up      = subview_2d(m_buffer.lw_flux_up);
      auto lw_flux_dn      = subview_2d(m_buffer.lw_flux_dn);
      auto sw_clnclrsky_flux_up      = subview_2d(m_buffer.sw_clnclrsky_flux_up);
      auto sw_clnclrsky_flux_dn      = subview_2d(m_buffer.sw_clnclrsky_flux_dn);
      auto sw_clnclrsky_flux_dn_dir  = subview_2d(m_buffer.sw_clnclrsky_flux_dn_dir);
      auto sw_clrsky_flux_up      = subview_2d(m_buffer.sw_clrsky_flux_up);
      auto sw_clrsky_flux_dn      = subview_2d(m_buffer.sw_clrsky_flux_dn);
      auto sw_clrsky_flux_dn_dir  = subview_2d(m_buffer.sw_clrsky_flux_dn_dir);
      auto sw_clnsky_flux_up      = subview_2d(m_buffer.sw_clnsky_flux_up);
      auto sw_clnsky_flux_dn      = subview_2d(m_buffer.sw_clnsky_flux_dn);
      auto sw_clnsky_flux_dn_dir  = subview_2d(m_buffer.sw_clnsky_flux_dn_dir);
      auto lw_clnclrsky_flux_up      = subview_2d(m_buffer.lw_clnclrsky_flux_up);
      auto lw_clnclrsky_flux_dn      = subview_2d(m_buffer.lw_clnclrsky_flux_dn);
      auto lw_clrsky_flux_up      = subview_2d(m_buffer.lw_clrsky_flux_up);
      auto lw_clrsky_flux_dn      = subview_2d(m_buffer.lw_clrsky_flux_dn);
      auto lw_clnsky_flux_up      = subview_2d(m_buffer.lw_clnsky_flux_up);
      auto lw_clnsky_flux_dn      = subview_2d(m_buffer.lw_clnsky_flux_dn);
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
      auto cld_tau_sw_bnd  = subview_3d(m_buffer.cld_tau_sw_bnd);
      auto cld_tau_lw_bnd  = subview_3d(m_buffer.cld_tau_lw_bnd);
      auto cld_tau_sw_gpt  = subview_3d(m_buffer.cld_tau_sw_gpt);
      auto cld_tau_lw_gpt  = subview_3d(m_buffer.cld_tau_lw_gpt);
                   );
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      ConvertToRrtmgpSubview conv = {beg, ncol};
      TIMED_INLINE_KERNEL(init_views,

      // Note, ncol will not necessary be m_col_chunk_size because the number of cols
      // will not always be evenly divided by m_col_chunk_size. In most cases, the
      // extra space will not cause any problems, but it does sometimes.
      auto p_lay_k           = conv.subview2d(d_pmid, m_buffer.p_lay_k, m_nlay);
      auto t_lay_k           = conv.subview2d(d_tmid, m_buffer.t_lay_k, m_nlay);
      auto p_lev_k           = conv.subview2d(d_pint, m_buffer.p_lev_k, m_nlay+1);
      auto z_del_k           = conv.subview2d(d_dz, m_buffer.z_del_k, m_nlay);
      auto p_del_k           = conv.subview2d(d_pdel, m_buffer.p_del_k, m_nlay);
      auto t_lev_k           = conv.subview2d(d_tint, m_buffer.t_lev_k, m_nlay+1);
      auto sfc_alb_dir_k     = m_buffer.sfc_alb_dir_k;
      auto sfc_alb_dif_k     = m_buffer.sfc_alb_dif_k;
      auto sfc_alb_dir_vis_k = conv.subview1d(d_sfc_alb_dir_vis);
      auto sfc_alb_dir_nir_k = conv.subview1d(d_sfc_alb_dir_nir);
      auto sfc_alb_dif_vis_k = conv.subview1d(d_sfc_alb_dif_vis);
      auto sfc_alb_dif_nir_k = conv.subview1d(d_sfc_alb_dif_nir);
      auto qc_k              = conv.subview2d(d_qc, m_buffer.qc_k, m_nlay);
      auto nc_k              = conv.subview2d(d_nc, m_buffer.nc_k, m_nlay);
      auto qi_k              = conv.subview2d(d_qi, m_buffer.qi_k, m_nlay);
      auto cldfrac_tot_k     = m_buffer.cldfrac_tot_k;
      auto rel_k             = conv.subview2d(d_rel, m_buffer.eff_radius_qc_k, m_nlay);
      auto rei_k             = conv.subview2d(d_rei, m_buffer.eff_radius_qi_k, m_nlay);
      auto sw_flux_up_k      = conv.subview2d(d_sw_flux_up, m_buffer.sw_flux_up_k, m_nlay+1);
      auto sw_flux_dn_k      = conv.subview2d(d_sw_flux_dn, m_buffer.sw_flux_dn_k, m_nlay+1);
      auto sw_flux_dn_dir_k  = conv.subview2d(d_sw_flux_dn_dir, m_buffer.sw_flux_dn_dir_k, m_nlay+1);
      auto lw_flux_up_k      = conv.subview2d(d_lw_flux_up, m_buffer.lw_flux_up_k, m_nlay+1);
      auto lw_flux_dn_k      = conv.subview2d(d_lw_flux_dn, m_buffer.lw_flux_dn_k, m_nlay+1);
      auto sw_clnclrsky_flux_up_k      = conv.subview2d(d_sw_clnclrsky_flux_up, m_buffer.sw_clnclrsky_flux_up_k, m_nlay+1);
      auto sw_clnclrsky_flux_dn_k      = conv.subview2d(d_sw_clnclrsky_flux_dn, m_buffer.sw_clnclrsky_flux_dn_k, m_nlay+1);
      auto sw_clnclrsky_flux_dn_dir_k  = conv.subview2d(d_sw_clnclrsky_flux_dn_dir, m_buffer.sw_clnclrsky_flux_dn_dir_k, m_nlay+1);
      auto sw_clrsky_flux_up_k      = conv.subview2d(d_sw_clrsky_flux_up, m_buffer.sw_clrsky_flux_up_k, m_nlay+1);
      auto sw_clrsky_flux_dn_k      = conv.subview2d(d_sw_clrsky_flux_dn, m_buffer.sw_clrsky_flux_dn_k, m_nlay+1);
      auto sw_clrsky_flux_dn_dir_k  = conv.subview2d(d_sw_clrsky_flux_dn_dir, m_buffer.sw_clrsky_flux_dn_dir_k, m_nlay+1);
      auto sw_clnsky_flux_up_k      = conv.subview2d(d_sw_clnsky_flux_up, m_buffer.sw_clnsky_flux_up_k, m_nlay+1);
      auto sw_clnsky_flux_dn_k      = conv.subview2d(d_sw_clnsky_flux_dn, m_buffer.sw_clnsky_flux_dn_k, m_nlay+1);
      auto sw_clnsky_flux_dn_dir_k  = conv.subview2d(d_sw_clnsky_flux_dn_dir, m_buffer.sw_clnsky_flux_dn_dir_k, m_nlay+1);
      auto lw_clnclrsky_flux_up_k   = conv.subview2d(d_lw_clnclrsky_flux_up, m_buffer.lw_clnclrsky_flux_up_k, m_nlay+1);
      auto lw_clnclrsky_flux_dn_k   = conv.subview2d(d_lw_clnclrsky_flux_dn, m_buffer.lw_clnclrsky_flux_dn_k, m_nlay+1);
      auto lw_clrsky_flux_up_k      = conv.subview2d(d_lw_clrsky_flux_up, m_buffer.lw_clrsky_flux_up_k, m_nlay+1);
      auto lw_clrsky_flux_dn_k      = conv.subview2d(d_lw_clrsky_flux_dn, m_buffer.lw_clrsky_flux_dn_k, m_nlay+1);
      auto lw_clnsky_flux_up_k      = conv.subview2d(d_lw_clnsky_flux_up, m_buffer.lw_clnsky_flux_up_k, m_nlay+1);
      auto lw_clnsky_flux_dn_k      = conv.subview2d(d_lw_clnsky_flux_dn, m_buffer.lw_clnsky_flux_dn_k, m_nlay+1);
      auto sw_bnd_flux_up_k  = m_buffer.sw_bnd_flux_up_k;
      auto sw_bnd_flux_dn_k  = m_buffer.sw_bnd_flux_dn_k;
      auto sw_bnd_flux_dir_k = m_buffer.sw_bnd_flux_dir_k;
      auto sw_bnd_flux_dif_k = m_buffer.sw_bnd_flux_dif_k;
      auto lw_bnd_flux_up_k  = m_buffer.lw_bnd_flux_up_k;
      auto lw_bnd_flux_dn_k  = m_buffer.lw_bnd_flux_dn_k;
      auto sfc_flux_dir_vis_k = conv.subview1d(d_sfc_flux_dir_vis);
      auto sfc_flux_dir_nir_k = conv.subview1d(d_sfc_flux_dir_nir);
      auto sfc_flux_dif_vis_k = conv.subview1d(d_sfc_flux_dif_vis);
      auto sfc_flux_dif_nir_k = conv.subview1d(d_sfc_flux_dif_nir);
      auto aero_tau_sw_k     = m_buffer.aero_tau_sw_k;
      auto aero_ssa_sw_k     = m_buffer.aero_ssa_sw_k;
      auto aero_g_sw_k       = m_buffer.aero_g_sw_k;
      auto aero_tau_lw_k     = m_buffer.aero_tau_lw_k;
      auto cld_tau_sw_bnd_k  = conv.subview3d(m_buffer.cld_tau_sw_bnd_k);
      auto cld_tau_lw_bnd_k  = conv.subview3d(m_buffer.cld_tau_lw_bnd_k);
      auto cld_tau_sw_gpt_k  = conv.subview3d(m_buffer.cld_tau_sw_gpt_k);
      auto cld_tau_lw_gpt_k  = conv.subview3d(m_buffer.cld_tau_lw_gpt_k);
                   );
#endif

      // Set gas concs to "view" only the first ncol columns
#ifdef RRTMGP_ENABLE_YAKL
      m_gas_concs.ncol = ncol;
      m_gas_concs.concs = subview_3d(gas_concs);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      m_gas_concs_k.ncol = ncol;
      m_gas_concs_k.concs = conv.subview3d(gas_concs_k);
#endif

      // Copy data from the FieldManager to the YAKL arrays
      {
        // Determine the cosine zenith angle
        // NOTE: Since we are bridging to F90 arrays this must be done on HOST and then
        //       deep copied to a device view.
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

        const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
        TIMED_KERNEL(
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
          const int itop = (p_mid(0) < p_mid(nlay-1)) ? 0 : nlay-1;
          const Real bc_top = T_mid(itop);
          const Real bc_bot = sqrt(sqrt(d_surf_lw_flux_up(icol)/stebol));
          if (itop == 0) {
            CO::compute_interface_values_linear(team, nlay, T_mid, dz, bc_top, bc_bot, T_int);
          } else {
            CO::compute_interface_values_linear(team, nlay, T_mid, dz, bc_bot, bc_top, T_int);
          }
          team.team_barrier();

#ifdef RRTMGP_ENABLE_YAKL
          mu0(i+1) = d_mu0(i);
          sfc_alb_dir_vis(i+1) = d_sfc_alb_dir_vis(icol);
          sfc_alb_dir_nir(i+1) = d_sfc_alb_dir_nir(icol);
          sfc_alb_dif_vis(i+1) = d_sfc_alb_dif_vis(icol);
          sfc_alb_dif_nir(i+1) = d_sfc_alb_dif_nir(icol);

          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
            p_lay(i+1,k+1)       = d_pmid(icol,k);
            t_lay(i+1,k+1)       = d_tmid(icol,k);
            z_del(i+1,k+1)       = d_dz(i,k);
            p_del(i+1,k+1)       = d_pdel(icol,k);
            qc(i+1,k+1)          = d_qc(icol,k);
            nc(i+1,k+1)          = d_nc(icol,k);
            qi(i+1,k+1)          = d_qi(icol,k);
            rel(i+1,k+1)         = d_rel(icol,k);
            rei(i+1,k+1)         = d_rei(icol,k);
            p_lev(i+1,k+1)       = d_pint(icol,k);
            t_lev(i+1,k+1)       = d_tint(i,k);
          });

          p_lev(i+1,nlay+1) = d_pint(icol,nlay);
          t_lev(i+1,nlay+1) = d_tint(i,nlay);

          // Note that RRTMGP expects ordering (col,lay,bnd) but the FM keeps things in (col,bnd,lay) order
          if (do_aerosol_rad) {
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
          } else {
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nswbands*nlay), [&] (const int&idx) {
                auto b = idx / nlay;
                auto k = idx % nlay;
                aero_tau_sw(i+1,k+1,b+1) = 0;
                aero_ssa_sw(i+1,k+1,b+1) = 0;
                aero_g_sw  (i+1,k+1,b+1) = 0;
            });
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlwbands*nlay), [&] (const int&idx) {
                auto b = idx / nlay;
                auto k = idx % nlay;
                aero_tau_lw(i+1,k+1,b+1) = 0;
            });
          }
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
#ifdef RRTMGP_LAYOUT_LEFT
          // Copy to layout left buffer views
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
            p_lay_k(i,k)       = d_pmid(icol,k);
            t_lay_k(i,k)       = d_tmid(icol,k);
            z_del_k(i,k)       = d_dz(i,k);
            p_del_k(i,k)       = d_pdel(icol,k);
            qc_k(i,k)          = d_qc(icol,k);
            nc_k(i,k)          = d_nc(icol,k);
            qi_k(i,k)          = d_qi(icol,k);
            rel_k(i,k)         = d_rel(icol,k);
            rei_k(i,k)         = d_rei(icol,k);
            p_lev_k(i,k)       = d_pint(icol,k);
            t_lev_k(i,k)       = d_tint(i,k);
          });

          p_lev_k(i,nlay) = d_pint(icol,nlay);
          t_lev_k(i,nlay) = d_tint(i,nlay);
#endif

          // Note that RRTMGP expects ordering (col,lay,bnd) but the FM keeps things in (col,bnd,lay) order
          if (do_aerosol_rad) {
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nswbands*nlay), [&] (const int&idx) {
              auto b = idx / nlay;
              auto k = idx % nlay;
              aero_tau_sw_k(i,k,b) = d_aero_tau_sw(icol,b,k);
              aero_ssa_sw_k(i,k,b) = d_aero_ssa_sw(icol,b,k);
              aero_g_sw_k  (i,k,b) = d_aero_g_sw  (icol,b,k);
            });
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlwbands*nlay), [&] (const int&idx) {
              auto b = idx / nlay;
              auto k = idx % nlay;
              aero_tau_lw_k(i,k,b) = d_aero_tau_lw(icol,b,k);
            });
          } else {
            // cuda complains (in warning only) about these being not allowed...
            // Kokkos::deep_copy(aero_tau_sw_k, 0);
            // Kokkos::deep_copy(aero_ssa_sw_k, 0);
            // Kokkos::deep_copy(aero_g_sw_k  , 0);
            // Kokkos::deep_copy(aero_tau_lw_k, 0);
            // So, do the manual labor instead:
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nswbands*nlay), [&] (const int&idx) {
              auto b = idx / nlay;
              auto k = idx % nlay;
              aero_tau_sw_k(i,k,b) = 0.0;
              aero_ssa_sw_k(i,k,b) = 0.0;
              aero_g_sw_k  (i,k,b) = 0.0;
            });
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlwbands*nlay), [&] (const int&idx) {
              auto b = idx / nlay;
              auto k = idx % nlay;
              aero_tau_lw_k(i,k,b) = 0.0;
            });
          }
#endif
        });
                     );
      }
      Kokkos::fence();
#ifdef RRTMGP_ENABLE_KOKKOS
      COMPARE_ALL_WRAP(std::vector<real3d>({aero_tau_sw, aero_ssa_sw, aero_g_sw, aero_tau_lw}),
                       std::vector<real3dk>({aero_tau_sw_k, aero_ssa_sw_k, aero_g_sw_k, aero_tau_lw_k}));
#endif


      // Populate GasConcs object to pass to RRTMGP driver
      // set_vmr requires the input array size to have the correct size,
      // and the last chunk may have less columns, so create a temp of
      // correct size that uses m_buffer.tmp2d's pointer
#ifdef RRTMGP_ENABLE_YAKL
      real2d tmp2d = subview_2d(m_buffer.tmp2d);
#endif
      for (int igas = 0; igas < m_ngas; igas++) {
        auto name = m_gas_names[igas];
        auto full_name = name + "_volume_mix_ratio";

        // 'o3' is marked as 'Required' rather than 'Computed', so we need to get the proper field
        auto f = name=="o3" ? get_field_in(full_name) : get_field_out(full_name);
        auto d_vmr = f.get_view<const Real**>();
#ifdef RRTMGP_ENABLE_KOKKOS
        auto tmp2d_k = conv.subview2d_impl(d_vmr, m_nlay);
#endif

#ifdef RRTMGP_ENABLE_YAKL
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
#endif

        // Populate GasConcs object
#ifdef RRTMGP_ENABLE_YAKL
        m_gas_concs.set_vmr(name, tmp2d);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
        COMPARE_WRAP(tmp2d, tmp2d_k);
        m_gas_concs_k.set_vmr(name, tmp2d_k);
#endif
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
#ifdef RRTMGP_ENABLE_YAKL
      auto lwp = m_buffer.lwp;
      auto iwp = m_buffer.iwp;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      auto lwp_k = m_buffer.lwp_k;
      auto iwp_k = m_buffer.iwp_k;
#endif
      if (not do_subcol_sampling) {
        const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
          const int i = team.league_rank();
          const int icol = i + beg;
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
#ifdef RRTMGP_ENABLE_YAKL
            if (d_cldfrac_tot(icol,k) > 0) {
              cldfrac_tot(i+1,k+1) = 1;
            } else {
              cldfrac_tot(i+1,k+1) = 0;
            }
            d_cldfrac_rad(icol,k) = cldfrac_tot(i+1,k+1);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
            if (d_cldfrac_tot(icol,k) > 0) {
              cldfrac_tot_k(i,k) = 1;
            } else {
              cldfrac_tot_k(i,k) = 0;
            }
            d_cldfrac_rad(icol,k) = cldfrac_tot_k(i,k);
#endif
          });
        });
      } else {
        const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
          const int i = team.league_rank();
          const int icol = i + beg;
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
#ifdef RRTMGP_ENABLE_YAKL
            cldfrac_tot(i+1,k+1) = d_cldfrac_tot(icol,k);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
            cldfrac_tot_k(i,k) = d_cldfrac_tot(icol,k);
#endif
            d_cldfrac_rad(icol,k) = d_cldfrac_tot(icol,k);
          });
        });
      }
      Kokkos::fence();
#ifdef RRTMGP_ENABLE_KOKKOS
      COMPARE_WRAP(cldfrac_tot, cldfrac_tot_k);
#endif

      // Compute layer cloud mass (per unit area)
#ifdef RRTMGP_ENABLE_YAKL
      rrtmgp::mixing_ratio_to_cloud_mass(qc, cldfrac_tot, p_del, lwp);
      rrtmgp::mixing_ratio_to_cloud_mass(qi, cldfrac_tot, p_del, iwp);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      interface_t::mixing_ratio_to_cloud_mass(qc_k, cldfrac_tot_k, p_del_k, lwp_k);
      interface_t::mixing_ratio_to_cloud_mass(qi_k, cldfrac_tot_k, p_del_k, iwp_k);
      COMPARE_ALL_WRAP(std::vector<real2d>({lwp, iwp}),
                       std::vector<real2dk>({lwp_k, iwp_k}));
#endif
      // Convert to g/m2 (needed by RRTMGP)
      {
      const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int i = team.league_rank();
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
          // Note that for YAKL arrays i and k start with index 1
#ifdef RRTMGP_ENABLE_YAKL
          lwp(i+1,k+1) *= 1e3;
          iwp(i+1,k+1) *= 1e3;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
          lwp_k(i,k) *= 1e3;
          iwp_k(i,k) *= 1e3;
#endif
        });
      });
      }
      Kokkos::fence();

      // Compute band-by-band surface_albedos. This is needed since
      // the AD passes broadband albedos, but rrtmgp require band-by-band.
#ifdef RRTMGP_ENABLE_YAKL
      TIMED_KERNEL(
      rrtmgp::compute_band_by_band_surface_albedos(
        ncol, nswbands,
        sfc_alb_dir_vis, sfc_alb_dir_nir,
        sfc_alb_dif_vis, sfc_alb_dif_nir,
        sfc_alb_dir, sfc_alb_dif);
                   );
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      TIMED_KERNEL(
      interface_t::compute_band_by_band_surface_albedos(
        ncol, nswbands,
        sfc_alb_dir_vis_k, sfc_alb_dir_nir_k,
        sfc_alb_dif_vis_k, sfc_alb_dif_nir_k,
        sfc_alb_dir_k, sfc_alb_dif_k);
                   );
      COMPARE_ALL_WRAP(std::vector<real2d>({sfc_alb_dir, sfc_alb_dif}),
                       std::vector<real2dk>({sfc_alb_dir_k, sfc_alb_dif_k}));
#endif
      // Compute cloud optical properties here?

      // Run RRTMGP driver
#ifdef RRTMGP_ENABLE_YAKL
      TIMED_KERNEL(
      rrtmgp::rrtmgp_main(
        ncol, m_nlay,
        p_lay, t_lay, p_lev, t_lev,
        m_gas_concs,
        sfc_alb_dir, sfc_alb_dif, mu0,
        lwp, iwp, rel, rei, cldfrac_tot,
        aero_tau_sw, aero_ssa_sw, aero_g_sw, aero_tau_lw,
        cld_tau_sw_bnd, cld_tau_lw_bnd,
        cld_tau_sw_gpt, cld_tau_lw_gpt,
        sw_flux_up       , sw_flux_dn       , sw_flux_dn_dir       , lw_flux_up       , lw_flux_dn,
        sw_clnclrsky_flux_up, sw_clnclrsky_flux_dn, sw_clnclrsky_flux_dn_dir,
        sw_clrsky_flux_up, sw_clrsky_flux_dn, sw_clrsky_flux_dn_dir,
        sw_clnsky_flux_up, sw_clnsky_flux_dn, sw_clnsky_flux_dn_dir,
        lw_clnclrsky_flux_up, lw_clnclrsky_flux_dn,
        lw_clrsky_flux_up, lw_clrsky_flux_dn,
        lw_clnsky_flux_up, lw_clnsky_flux_dn,
        sw_bnd_flux_up   , sw_bnd_flux_dn   , sw_bnd_flux_dir      , lw_bnd_flux_up   , lw_bnd_flux_dn,
        eccf, m_atm_logger,
        m_extra_clnclrsky_diag, m_extra_clnsky_diag
      );
                   );
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      TIMED_KERNEL(
      interface_t::rrtmgp_main(
        ncol, m_nlay,
        p_lay_k, t_lay_k, p_lev_k, t_lev_k,
        m_gas_concs_k,
        sfc_alb_dir_k, sfc_alb_dif_k, d_mu0,
        lwp_k, iwp_k, rel_k, rei_k, cldfrac_tot_k,
        aero_tau_sw_k, aero_ssa_sw_k, aero_g_sw_k, aero_tau_lw_k,
        cld_tau_sw_bnd_k, cld_tau_lw_bnd_k,
        cld_tau_sw_gpt_k, cld_tau_lw_gpt_k,
        sw_flux_up_k, sw_flux_dn_k, sw_flux_dn_dir_k, lw_flux_up_k, lw_flux_dn_k,
        sw_clnclrsky_flux_up_k, sw_clnclrsky_flux_dn_k, sw_clnclrsky_flux_dn_dir_k,
        sw_clrsky_flux_up_k, sw_clrsky_flux_dn_k, sw_clrsky_flux_dn_dir_k,
        sw_clnsky_flux_up_k, sw_clnsky_flux_dn_k, sw_clnsky_flux_dn_dir_k,
        lw_clnclrsky_flux_up_k, lw_clnclrsky_flux_dn_k,
        lw_clrsky_flux_up_k, lw_clrsky_flux_dn_k,
        lw_clnsky_flux_up_k, lw_clnsky_flux_dn_k,
        sw_bnd_flux_up_k, sw_bnd_flux_dn_k, sw_bnd_flux_dir_k, lw_bnd_flux_up_k, lw_bnd_flux_dn_k,
        eccf, m_atm_logger,
        m_extra_clnclrsky_diag, m_extra_clnsky_diag
      );
                   );
      COMPARE_ALL_WRAP(std::vector<real2d>({
        sw_flux_up, sw_flux_dn, sw_flux_dn_dir, lw_flux_up, lw_flux_dn,
        sw_clnclrsky_flux_up, sw_clnclrsky_flux_dn, sw_clnclrsky_flux_dn_dir,
        sw_clrsky_flux_up, sw_clrsky_flux_dn, sw_clrsky_flux_dn_dir,
        sw_clnsky_flux_up, sw_clnsky_flux_dn, sw_clnsky_flux_dn_dir,
        lw_clnclrsky_flux_up, lw_clnclrsky_flux_dn,
        lw_clrsky_flux_up, lw_clrsky_flux_dn,
        lw_clnsky_flux_up, lw_clnsky_flux_dn}),
                       std::vector<real2dk>({
        sw_flux_up_k, sw_flux_dn_k, sw_flux_dn_dir_k, lw_flux_up_k, lw_flux_dn_k,
        sw_clnclrsky_flux_up_k, sw_clnclrsky_flux_dn_k, sw_clnclrsky_flux_dn_dir_k,
        sw_clrsky_flux_up_k, sw_clrsky_flux_dn_k, sw_clrsky_flux_dn_dir_k,
        sw_clnsky_flux_up_k, sw_clnsky_flux_dn_k, sw_clnsky_flux_dn_dir_k,
        lw_clnclrsky_flux_up_k, lw_clnclrsky_flux_dn_k,
        lw_clrsky_flux_up_k, lw_clrsky_flux_dn_k,
        lw_clnsky_flux_up_k, lw_clnsky_flux_dn_k}));
      COMPARE_ALL_WRAP(std::vector<real3d>({sw_bnd_flux_up, sw_bnd_flux_dn, sw_bnd_flux_dir, lw_bnd_flux_up, lw_bnd_flux_dn}),
                       std::vector<real3dk>({sw_bnd_flux_up_k, sw_bnd_flux_dn_k, sw_bnd_flux_dir_k, lw_bnd_flux_up_k, lw_bnd_flux_dn_k}));
#endif

      // Update heating tendency
#ifdef RRTMGP_ENABLE_YAKL
      TIMED_INLINE_KERNEL(heating_tendency,
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
                   );
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      TIMED_INLINE_KERNEL(heating_tendency,
      auto sw_heating_k  = m_buffer.sw_heating_k;
      auto lw_heating_k  = m_buffer.lw_heating_k;
      rrtmgp::compute_heating_rate(
        sw_flux_up_k, sw_flux_dn_k, p_del_k, sw_heating_k
      );
      rrtmgp::compute_heating_rate(
        lw_flux_up_k, lw_flux_dn_k, p_del_k, lw_heating_k
      );
      {
        const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
          const int idx = team.league_rank();
          const int icol = idx+beg;
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& ilay) {
            // Combine SW and LW heating into a net heating tendency; use d_rad_heating_pdel temporarily
            // Note that for YAKL arrays i and k start with index 1
            d_rad_heating_pdel(icol,ilay) = sw_heating_k(idx,ilay) + lw_heating_k(idx,ilay);
          });
        });
      }
      Kokkos::fence();
                   );
      COMPARE_ALL_WRAP(std::vector<real2d>({sw_heating, lw_heating}),
                       std::vector<real2dk>({sw_heating_k, lw_heating_k}));
#endif

      // Index to surface (bottom of model); used to get surface fluxes below
#ifdef RRTMGP_ENABLE_YAKL
      const int kbot = nlay+1;

      TIMED_KERNEL(
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
                   );
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      const int kbot_k = nlay;

      TIMED_KERNEL(
      // Compute diffuse flux as difference between total and direct
      Kokkos::parallel_for(Kokkos::RangePolicy<ExeSpace>(0,nswbands*(nlay+1)*ncol),
                           KOKKOS_LAMBDA (const int idx) {
        // CAREFUL: these are YAKL arrays, with "LayoutLeft". So make the indices stride accordingly, and add 1.
        const int ibnd = (idx / ncol) / (nlay+1);
        const int ilev = (idx / ncol) % (nlay+1);
        const int icol =  idx % ncol;
        sw_bnd_flux_dif_k(icol,ilev,ibnd) = sw_bnd_flux_dn_k(icol,ilev,ibnd) - sw_bnd_flux_dir_k(icol,ilev,ibnd);
      });
      // Compute surface fluxes
      interface_t::compute_broadband_surface_fluxes(
          ncol, kbot_k, nswbands,
          sw_bnd_flux_dir_k, sw_bnd_flux_dif_k,
          sfc_flux_dir_vis_k, sfc_flux_dir_nir_k,
          sfc_flux_dif_vis_k, sfc_flux_dif_nir_k
      );
                   );
      COMPARE_ALL_WRAP(std::vector<real1d>({sfc_flux_dir_vis, sfc_flux_dir_nir, sfc_flux_dif_vis, sfc_flux_dif_nir}),
                       std::vector<real1dk>({sfc_flux_dir_vis_k, sfc_flux_dir_nir_k, sfc_flux_dif_vis_k, sfc_flux_dif_nir_k}));
#endif

      // Compute diagnostic total cloud area (vertically-projected cloud cover)
#ifdef RRTMGP_ENABLE_YAKL
      TIMED_KERNEL(
      real1d cldlow ("cldlow", d_cldlow.data() + m_col_chunk_beg[ic], ncol);
      real1d cldmed ("cldmed", d_cldmed.data() + m_col_chunk_beg[ic], ncol);
      real1d cldhgh ("cldhgh", d_cldhgh.data() + m_col_chunk_beg[ic], ncol);
      real1d cldtot ("cldtot", d_cldtot.data() + m_col_chunk_beg[ic], ncol);
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
                   );
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      TIMED_KERNEL(
      real1dk cldlow_k (d_cldlow.data() + m_col_chunk_beg[ic], ncol);
      real1dk cldmed_k (d_cldmed.data() + m_col_chunk_beg[ic], ncol);
      real1dk cldhgh_k (d_cldhgh.data() + m_col_chunk_beg[ic], ncol);
      real1dk cldtot_k (d_cldtot.data() + m_col_chunk_beg[ic], ncol);
      // NOTE: limits for low, mid, and high clouds are mostly taken from EAM F90 source, with the
      // exception that I removed the restriction on low clouds to be above (numerically lower pressures)
      // 1200 hPa, and on high clouds to be below (numerically high pressures) 50 hPa. This probably
      // does not matter in practice, as clouds probably should not be produced above 50 hPa and we
      // should not be encountering surface pressure above 1200 hPa, but in the event that things go off
      // the rails we might want to look at these still.
      interface_t::compute_cloud_area(ncol, nlay, nlwgpts, 700e2, std::numeric_limits<Real>::max(), p_lay_k, cld_tau_lw_gpt_k, cldlow_k);
      interface_t::compute_cloud_area(ncol, nlay, nlwgpts, 400e2,                            700e2, p_lay_k, cld_tau_lw_gpt_k, cldmed_k);
      interface_t::compute_cloud_area(ncol, nlay, nlwgpts,     0,                            400e2, p_lay_k, cld_tau_lw_gpt_k, cldhgh_k);
      interface_t::compute_cloud_area(ncol, nlay, nlwgpts,     0, std::numeric_limits<Real>::max(), p_lay_k, cld_tau_lw_gpt_k, cldtot_k);
                   );
      COMPARE_ALL_WRAP(std::vector<real1d>({cldlow, cldmed, cldhgh, cldtot}),
                       std::vector<real1dk>({cldlow_k, cldmed_k, cldhgh_k, cldtot_k}));
#endif

      // Compute cloud-top diagnostics following AeroCOM recommendation
#ifdef RRTMGP_ENABLE_YAKL
      TIMED_INLINE_KERNEL(cloud_top,

      // Get visible 0.67 micron band for COSP
      auto idx_067 = rrtmgp::get_wavelength_index_sw(0.67e-6);
      // Get IR 10.5 micron band for COSP
      auto idx_105 = rrtmgp::get_wavelength_index_lw(10.5e-6);

      // Compute cloud-top diagnostics following AeroCom recommendation
      real1d T_mid_at_cldtop ("T_mid_at_cldtop", d_T_mid_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1d p_mid_at_cldtop ("p_mid_at_cldtop", d_p_mid_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1d cldfrac_ice_at_cldtop ("cldfrac_ice_at_cldtop", d_cldfrac_ice_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1d cldfrac_liq_at_cldtop ("cldfrac_liq_at_cldtop", d_cldfrac_liq_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1d cldfrac_tot_at_cldtop ("cldfrac_tot_at_cldtop", d_cldfrac_tot_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1d cdnc_at_cldtop ("cdnc_at_cldtop", d_cdnc_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1d eff_radius_qc_at_cldtop ("eff_radius_qc_at_cldtop", d_eff_radius_qc_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1d eff_radius_qi_at_cldtop ("eff_radius_qi_at_cldtop", d_eff_radius_qi_at_cldtop.data() + m_col_chunk_beg[ic], ncol);

      rrtmgp::compute_aerocom_cloudtop(
          ncol, nlay, t_lay, p_lay, p_del, z_del, qc, qi, rel, rei, cldfrac_tot,
          nc, T_mid_at_cldtop, p_mid_at_cldtop, cldfrac_ice_at_cldtop,
          cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
          eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop);
                          );
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      TIMED_INLINE_KERNEL(cloud_top,
      // Get visible 0.67 micron band for COSP
      auto idx_067_k = interface_t::get_wavelength_index_sw_k(0.67e-6);
      // Get IR 10.5 micron band for COSP
      auto idx_105_k = interface_t::get_wavelength_index_lw_k(10.5e-6);

      real1dk T_mid_at_cldtop_k (d_T_mid_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1dk p_mid_at_cldtop_k (d_p_mid_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1dk cldfrac_ice_at_cldtop_k (d_cldfrac_ice_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1dk cldfrac_liq_at_cldtop_k (d_cldfrac_liq_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1dk cldfrac_tot_at_cldtop_k (d_cldfrac_tot_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1dk cdnc_at_cldtop_k (d_cdnc_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1dk eff_radius_qc_at_cldtop_k (d_eff_radius_qc_at_cldtop.data() + m_col_chunk_beg[ic], ncol);
      real1dk eff_radius_qi_at_cldtop_k (d_eff_radius_qi_at_cldtop.data() + m_col_chunk_beg[ic], ncol);

      interface_t::compute_aerocom_cloudtop(
          ncol, nlay, t_lay_k, p_lay_k, p_del_k, z_del_k, qc_k, qi_k, rel_k, rei_k, cldfrac_tot_k,
          nc_k, T_mid_at_cldtop_k, p_mid_at_cldtop_k, cldfrac_ice_at_cldtop_k,
          cldfrac_liq_at_cldtop_k, cldfrac_tot_at_cldtop_k, cdnc_at_cldtop_k,
          eff_radius_qc_at_cldtop_k, eff_radius_qi_at_cldtop_k);
                          );
      COMPARE_ALL_WRAP(std::vector<real1d>({
            T_mid_at_cldtop, p_mid_at_cldtop, cldfrac_ice_at_cldtop,
            cldfrac_liq_at_cldtop, cldfrac_tot_at_cldtop, cdnc_at_cldtop,
            eff_radius_qc_at_cldtop, eff_radius_qi_at_cldtop}),
                       std::vector<real1dk>({
            T_mid_at_cldtop_k, p_mid_at_cldtop_k, cldfrac_ice_at_cldtop_k,
            cldfrac_liq_at_cldtop_k, cldfrac_tot_at_cldtop_k, cdnc_at_cldtop_k,
            eff_radius_qc_at_cldtop_k, eff_radius_qi_at_cldtop_k}));
#endif

      // Copy output data back to FieldManager
      const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_nlay);
#ifdef RRTMGP_ENABLE_YAKL
      TIMED_KERNEL(
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int i = team.league_rank();
        const int icol = i + beg;
        d_sfc_flux_dir_nir(icol) = sfc_flux_dir_nir(i+1);
        d_sfc_flux_dir_vis(icol) = sfc_flux_dir_vis(i+1);
        d_sfc_flux_dif_nir(icol) = sfc_flux_dif_nir(i+1);
        d_sfc_flux_dif_vis(icol) = sfc_flux_dif_vis(i+1);
        d_sfc_flux_sw_net(icol)  = sw_flux_dn(i+1,kbot) - sw_flux_up(i+1,kbot);
        d_sfc_flux_lw_dn(icol)   = lw_flux_dn(i+1,kbot);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay+1), [&] (const int& k) {
          d_sw_flux_up(icol,k)            = sw_flux_up(i+1,k+1);
          d_sw_flux_dn(icol,k)            = sw_flux_dn(i+1,k+1);
          d_sw_flux_dn_dir(icol,k)        = sw_flux_dn_dir(i+1,k+1);
          d_lw_flux_up(icol,k)            = lw_flux_up(i+1,k+1);
          d_lw_flux_dn(icol,k)            = lw_flux_dn(i+1,k+1);
          d_sw_clnclrsky_flux_up(icol,k)     = sw_clnclrsky_flux_up(i+1,k+1);
          d_sw_clnclrsky_flux_dn(icol,k)     = sw_clnclrsky_flux_dn(i+1,k+1);
          d_sw_clnclrsky_flux_dn_dir(icol,k) = sw_clnclrsky_flux_dn_dir(i+1,k+1);
          d_sw_clrsky_flux_up(icol,k)     = sw_clrsky_flux_up(i+1,k+1);
          d_sw_clrsky_flux_dn(icol,k)     = sw_clrsky_flux_dn(i+1,k+1);
          d_sw_clrsky_flux_dn_dir(icol,k) = sw_clrsky_flux_dn_dir(i+1,k+1);
          d_sw_clnsky_flux_up(icol,k)     = sw_clnsky_flux_up(i+1,k+1);
          d_sw_clnsky_flux_dn(icol,k)     = sw_clnsky_flux_dn(i+1,k+1);
          d_sw_clnsky_flux_dn_dir(icol,k) = sw_clnsky_flux_dn_dir(i+1,k+1);
          d_lw_clnclrsky_flux_up(icol,k)     = lw_clnclrsky_flux_up(i+1,k+1);
          d_lw_clnclrsky_flux_dn(icol,k)     = lw_clnclrsky_flux_dn(i+1,k+1);
          d_lw_clrsky_flux_up(icol,k)     = lw_clrsky_flux_up(i+1,k+1);
          d_lw_clrsky_flux_dn(icol,k)     = lw_clrsky_flux_dn(i+1,k+1);
          d_lw_clnsky_flux_up(icol,k)     = lw_clnsky_flux_up(i+1,k+1);
          d_lw_clnsky_flux_dn(icol,k)     = lw_clnsky_flux_dn(i+1,k+1);
        });
        // Extract optical properties for COSP
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
          d_dtau067(icol,k) = cld_tau_sw_bnd(i+1,k+1,idx_067);
          d_dtau105(icol,k) = cld_tau_lw_bnd(i+1,k+1,idx_105);
        });
        if (d_sw_clrsky_flux_dn(icol,0) > 0) {
            d_sunlit(icol) = 1.0;
        } else {
            d_sunlit(icol) = 0.0;
        }
      });
                   );
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
      TIMED_KERNEL(
      Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int i = team.league_rank();
        const int icol = i + beg;
        d_sfc_flux_sw_net(icol)  = sw_flux_dn_k(i,kbot_k) - sw_flux_up_k(i,kbot_k);
        d_sfc_flux_lw_dn(icol)   = lw_flux_dn_k(i,kbot_k);
#ifdef RRTMGP_LAYOUT_LEFT
        // Copy from layout left buffer views back to layout right fields
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay+1), [&] (const int& k) {
          d_sw_flux_up(icol,k)            = sw_flux_up_k(i,k);
          d_sw_flux_dn(icol,k)            = sw_flux_dn_k(i,k);
          d_sw_flux_dn_dir(icol,k)        = sw_flux_dn_dir_k(i,k);
          d_lw_flux_up(icol,k)            = lw_flux_up_k(i,k);
          d_lw_flux_dn(icol,k)            = lw_flux_dn_k(i,k);
          d_sw_clnclrsky_flux_up(icol,k)     = sw_clnclrsky_flux_up_k(i,k);
          d_sw_clnclrsky_flux_dn(icol,k)     = sw_clnclrsky_flux_dn_k(i,k);
          d_sw_clnclrsky_flux_dn_dir(icol,k) = sw_clnclrsky_flux_dn_dir_k(i,k);
          d_sw_clrsky_flux_up(icol,k)     = sw_clrsky_flux_up_k(i,k);
          d_sw_clrsky_flux_dn(icol,k)     = sw_clrsky_flux_dn_k(i,k);
          d_sw_clrsky_flux_dn_dir(icol,k) = sw_clrsky_flux_dn_dir_k(i,k);
          d_sw_clnsky_flux_up(icol,k)     = sw_clnsky_flux_up_k(i,k);
          d_sw_clnsky_flux_dn(icol,k)     = sw_clnsky_flux_dn_k(i,k);
          d_sw_clnsky_flux_dn_dir(icol,k) = sw_clnsky_flux_dn_dir_k(i,k);
          d_lw_clnclrsky_flux_up(icol,k)     = lw_clnclrsky_flux_up_k(i,k);
          d_lw_clnclrsky_flux_dn(icol,k)     = lw_clnclrsky_flux_dn_k(i,k);
          d_lw_clrsky_flux_up(icol,k)     = lw_clrsky_flux_up_k(i,k);
          d_lw_clrsky_flux_dn(icol,k)     = lw_clrsky_flux_dn_k(i,k);
          d_lw_clnsky_flux_up(icol,k)     = lw_clnsky_flux_up_k(i,k);
          d_lw_clnsky_flux_dn(icol,k)     = lw_clnsky_flux_dn_k(i,k);
        });
#endif
        // Extract optical properties for COSP
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
          d_dtau067(icol,k) = cld_tau_sw_bnd_k(i,k,idx_067_k);
          d_dtau105(icol,k) = cld_tau_lw_bnd_k(i,k,idx_105_k);
        });
        if (d_sw_clrsky_flux_dn(icol,0) > 0) {
            d_sunlit(icol) = 1.0;
        } else {
            d_sunlit(icol) = 0.0;
        }
      });
                   );
#ifdef RRTMGP_ENABLE_YAKL
      // Sync back to gas_concs_k
      real3dk temp(gas_concs_k, std::make_pair(0, ncol), Kokkos::ALL, Kokkos::ALL);
      Kokkos::deep_copy(temp, m_gas_concs_k.concs);
#endif
#endif
    } // loop over chunk

    // Restore the refCounted array.
#ifdef RRTMGP_ENABLE_YAKL
    m_gas_concs.concs = gas_concs;
    m_gas_concs.ncol = orig_ncol;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
    m_gas_concs_k.concs = gas_concs_k;
    m_gas_concs_k.ncol = orig_ncol_k;
#endif
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
#ifdef RRTMGP_ENABLE_YAKL
  m_gas_concs.reset();
  rrtmgp::rrtmgp_finalize();
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
  m_gas_concs_k.reset();
  // Finalize the interface, passing a bool for rank 0
  // to print info about memory stats on that rank
  interface_t::rrtmgp_finalize(m_comm.am_i_root());
#endif

  finalize_kls();
}
// =========================================================================================

}  // namespace scream
