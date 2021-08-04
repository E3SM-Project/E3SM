#include "atmosphere_prescribed_aerosol.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{
  using namespace spa;
// =========================================================================================
SPA::SPA (const ekat::Comm& comm, const ekat::ParameterList& params)
 : m_spa_comm (comm)
 , m_spa_params (params)
{
  // Nothing to do here
}

// =========================================================================================
void SPA::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg/kg;
  Q.set_string("kg/kg");
  Units nondim(0,0,0,0,0,0,0);

  const auto& grid_name = m_spa_params.get<std::string>("Grid");
  auto grid  = grids_manager->get_grid(grid_name);
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces 
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  FieldLayout scalar2d_layout     { {COL},     {m_num_cols} };
  FieldLayout scalar1d_layout_mid { {LEV},     {m_num_levs} };

  // Set of fields used strictly as input
  constexpr int ps = Pack::n;
  add_field<Required>("PS"         , scalar2d_layout,     Pa,     grid_name, ps);
  add_field<Required>("PS_beg"     , scalar2d_layout,     Pa,     grid_name, ps);
  add_field<Required>("PS_end"     , scalar2d_layout,     Pa,     grid_name, ps);
  add_field<Required>("hyam"       , scalar1d_layout_mid, nondim, grid_name, ps);
  add_field<Required>("hybm"       , scalar1d_layout_mid, nondim, grid_name, ps);
  add_field<Updated>("CCN3_beg"   , scalar3d_layout_mid, 1/kg,   grid_name, ps); // TODO:  This and p_mid should actually be "required"
  add_field<Updated>("CCN3_end"   , scalar3d_layout_mid, 1/kg,   grid_name, ps); // TODO:  This and p_mid should actually be "required"
  add_field<Updated>("p_mid"      , scalar3d_layout_mid, 1/kg,   grid_name, ps);

  // Set of fields used strictly as output
  add_field<Computed>("nc_activated", scalar3d_layout_mid, 1/kg,   grid_name,ps);

  // Set of fields used as input and output
  // - There are no fields used as both input and output.
}
// =========================================================================================
int SPA::requested_buffer_size_in_bytes() const
{
  const Int nk_pack    = ekat::npack<Spack>(m_num_levs);
  const Int nk_pack_p1 = ekat::npack<Spack>(m_num_levs+1);

  // Number of Reals needed by local views in the interface
  const int interface_request =
      // 1d view scalar, size (ncol)
      Buffer::num_1d_scalar*m_num_cols*sizeof(Real) +
      // 2d view packed, size (ncol, nlev_packs)
      Buffer::num_2d_vector*m_num_cols*nk_pack*sizeof(Spack) +
      Buffer::num_2dp1_vector*m_num_cols*nk_pack_p1*sizeof(Spack);

  // Number of Reals needed by the WorkspaceManager
  const auto policy       = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nk_pack);
  const int wsm_request   = WSM::get_total_bytes_needed(nk_pack, 3, policy);

  return interface_request + wsm_request;
}

// =========================================================================================
void SPA::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d scalar views
  m_buffer.ps_src = decltype(m_buffer.ps_src)(mem, m_num_cols);
  mem += m_buffer.ps_src.size();

  // 2d scalar views
//  m_buffer.col_location = decltype(m_buffer.col_location)(mem, m_num_cols, 3);
//  mem += m_buffer.col_location.size();

  Spack* s_mem = reinterpret_cast<Spack*>(mem);

  // 2d packed views
  const Int nk_pack    = ekat::npack<Spack>(m_num_levs);
  const Int nk_pack_p1 = ekat::npack<Spack>(m_num_levs+1);

  m_buffer.p_mid_src = decltype(m_buffer.p_mid_src)(s_mem, m_num_cols, nk_pack);
  s_mem += m_buffer.p_mid_src.size();
  m_buffer.ccn3_src = decltype(m_buffer.ccn3_src)(s_mem, m_num_cols, nk_pack);
  s_mem += m_buffer.ccn3_src.size();

  // WSM data
  m_buffer.wsm_data = s_mem;

  // Compute workspace manager size to check used memory
  // vs. requested memory
  const auto policy  = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nk_pack);
  const int wsm_size = WSM::get_total_bytes_needed(nk_pack, 3, policy)/sizeof(Spack);
  s_mem += wsm_size;

  int used_mem = (reinterpret_cast<Real*>(s_mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for SPA.");
}

// =========================================================================================
void SPA::initialize_impl (const util::TimeStamp& /* t0 */)
{
}

// =========================================================================================
void SPA::run_impl (const Real dt)
{
  // Calculate ice cloud fraction and total cloud fraction given the liquid cloud fraction
  // and the ice mass mixing ratio. 
  const auto& hyam         = m_spa_fields_in["hyam"].get_reshaped_view<const Pack*>();
  const auto& hybm         = m_spa_fields_in["hybm"].get_reshaped_view<const Pack*>();
  const auto& PS           = m_spa_fields_in["PS"].get_reshaped_view<const Real*>();
  const auto& ps_beg       = m_spa_fields_in["PS_beg"].get_reshaped_view<const Real*>();
  const auto& ps_end       = m_spa_fields_in["PS_end"].get_reshaped_view<const Real*>();
  const auto& p_mid        = m_spa_fields_out["p_mid"].get_reshaped_view<Pack**>(); // TODO this and ccn need to be const when the EKAT API is fixed
  const auto& ccn3_beg     = m_spa_fields_out["CCN3_beg"].get_reshaped_view<Pack**>();  // TODO: when time interp is put in build ccn_src from two timesteps of data.
  const auto& ccn3_end     = m_spa_fields_out["CCN3_end"].get_reshaped_view<Pack**>();  // TODO: when time interp is put in build ccn_src from two timesteps of data.
  const auto& nc_activated = m_spa_fields_out["nc_activated"].get_reshaped_view<Pack**>();
  auto p_mid_src    = m_buffer.p_mid_src;
  auto ps_src       = m_buffer.ps_src;
  auto ccn3_src     = m_buffer.ccn3_src;

  // First Step: Horizontal Interpolation if needed - Skip for Now

  // Second Step: Temporal Interpolation
  /* Gather time information for interpolation */
  auto ts = timestamp();
  Real t_now = ts.get_julian_day();
  Real t_beg_month = ts.get_julian_day(ts.get_years(),ts.get_months(),0,0);
  Real days_this_month = ts.get_dpm();
  
  /* First determine PS for the source data at this time */
  {
  using ExeSpace = typename KT::ExeSpace;
  Kokkos::parallel_for("time-interp-PS",
    m_num_cols,
    KOKKOS_LAMBDA(const int& i) {
      ps_src(i) = ps_beg(i) + (t_now - t_beg_month) * ( ps_end(i) - ps_beg(i) )/ days_this_month;
    });
  }
  /* Next, interpolate all of the ccn data */
  {
  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(m_num_levs);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_num_cols, nk_pack);
  Kokkos::parallel_for(
    "time-interp-aero",
    policy,
    KOKKOS_LAMBDA(const KT::MemberType& team) {

      const Int i = team.league_rank();
      const auto& ccn3_src_sub = ekat::subview(ccn3_src,i);
      const auto& ccn3_beg_sub = ekat::subview(ccn3_beg,i);
      const auto& ccn3_end_sub = ekat::subview(ccn3_end,i);
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {
          ccn3_src_sub(k) = ccn3_beg_sub(k) + (t_now - t_beg_month) * ( ccn3_end_sub(k) - ccn3_beg_sub(k) )/ days_this_month;
      });
    });
  }

  // Final Step: Vertical Interpolation of aerosol data at all columns
  /* First construct the source pressure interface levels */
  SPAFunc::reconstruct_pressure_profile(m_num_cols,m_num_levs,hyam,hybm,ps_src,p_mid_src);
  /* Next, use EKAT interpolation to remap CCN data onto current set of pressure levels */
  SPAFunc::aero_vertical_remap(m_num_cols,m_num_levs,m_num_levs,p_mid_src,p_mid,ccn3_src,nc_activated);  // switch to ccn3_src when ready.
 
  // TODO, when all three interpolations are finished in run_impl, create a "SPA_main" routine that calls all three in order and simplify run_impl. 

  // Advance current timestamp.
  ts += dt;
  for (auto& f : m_spa_fields_out) {
    f.second.get_header().get_tracking().update_time_stamp(ts);
  }
}

// =========================================================================================
void SPA::finalize_impl()
{
  // Do nothing
}

void SPA::set_required_field_impl (const Field<const Real>& f) {

  const auto& name = f.get_header().get_identifier().name();
  m_spa_fields_in.emplace(name,f);

  // Add myself as customer to the field
  add_me_as_customer(f);
}

void SPA::set_computed_field_impl (const Field<      Real>& f) {

  const auto& name = f.get_header().get_identifier().name();
  m_spa_fields_out.emplace(name,f);

  // Add myself as provider for the field
  add_me_as_provider(f);
}

} // namespace scream
