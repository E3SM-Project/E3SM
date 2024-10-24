#include "eamxx_tms_process_interface.hpp"

#include "physics/tms/tms_functions.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{

// =========================================================================================
TurbulentMountainStress::TurbulentMountainStress (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Do nothing
}

// =========================================================================================
void TurbulentMountainStress::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // Define some useful units. The units of mixing ratio
  // Q are technically non-dimensional. Nevertheless,
  // for output reasons, we like to see 'kg/kg'.
  const auto nondim = Units::nondimensional();
  const auto m2 = pow(m,2);

  // Initialize grid from grids manager
  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  EKAT_REQUIRE_MSG(grid_name=="Physics PG2",
                   "Error! TMS process can only be used with \"Physics PG2\" physics grid. "
                   "Current physics grid is "+grid_name+".\n");

  m_ncols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_nlevs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Add required/computed fields
  auto scalar2d     = m_grid->get_2d_scalar_layout();
  auto vector2d     = m_grid->get_2d_vector_layout(2);
  auto scalar3d_mid = m_grid->get_3d_scalar_layout(true);
  auto vector3d_mid = m_grid->get_3d_vector_layout(true,2);

  constexpr int ps = Spack::n;
  add_field<Required>("horiz_winds",    vector3d_mid, m/s,    grid_name,            ps);
  add_field<Required>("T_mid",          scalar3d_mid, K,      grid_name,            ps);
  add_field<Required>("p_mid",          scalar3d_mid, Pa,     grid_name,            ps);
  add_field<Required>("pseudo_density", scalar3d_mid, Pa,     grid_name,            ps);
  add_field<Required>("sgh30",          scalar2d    , m,      grid_name);
  add_field<Required>("landfrac",       scalar2d    , nondim, grid_name);
  add_tracer<Required>("qv", m_grid, kg/kg, ps);

  add_field<Computed>("surf_drag_coeff_tms", scalar2d, kg/(m2*s), grid_name);
  add_field<Computed>("wind_stress_tms",     vector2d, N/m2,      grid_name);
}

// =========================================================================================
void TurbulentMountainStress::initialize_impl (const RunType /* run_type */)
{
  // Do nothing
}

// =========================================================================================
void TurbulentMountainStress::run_impl (const double /* dt */)
{
  // Helper views
  const auto pseudo_density = get_field_in("pseudo_density").get_view<const Spack**>();
  const auto qv             = get_field_in("qv").get_view<const Spack**>();
  const auto dz             = m_buffer.dz;
  const auto z_int          = m_buffer.z_int;

  // Input views
  const auto horiz_winds = get_field_in("horiz_winds").get_view<const Spack***>();
  const auto T_mid       = get_field_in("T_mid").get_view<const Spack**>();
  const auto p_mid       = get_field_in("p_mid").get_view<const Spack**>();
  const auto sgh30       = get_field_in("sgh30").get_view<const Real*>();
  const auto landfrac    = get_field_in("landfrac").get_view<const Real*>();
  const auto exner       = m_buffer.exner;
  const auto z_mid       = m_buffer.z_mid;

  // Output views
  const auto surf_drag_coeff_tms = get_field_out("surf_drag_coeff_tms").get_view<Real*>();
  const auto wind_stress_tms     = get_field_out("wind_stress_tms").get_view<Real**>();

  // Preprocess inputs
  const int ncols = m_ncols;
  const int nlevs = m_nlevs;
  const int nlev_packs = ekat::npack<Spack>(nlevs);
  // calculate_z_int contains a team-level parallel_scan, which requires a special policy
  const auto scan_policy = ekat::ExeSpaceUtils<TMSFunctions::KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncols, nlev_packs);
  Kokkos::parallel_for(scan_policy, KOKKOS_LAMBDA (const TMSFunctions::KT::MemberType& team) {
    const int i = team.league_rank();

    const auto p_mid_i = ekat::subview(p_mid, i);
    const auto exner_i = ekat::subview(exner, i);
    const auto pseudo_density_i = ekat::subview(pseudo_density, i);
    const auto T_mid_i = ekat::subview(T_mid, i);
    const auto qv_i = ekat::subview(qv, i);
    const auto dz_i = ekat::subview(dz, i);
    const auto z_int_i = ekat::subview(z_int, i);
    const auto z_mid_i = ekat::subview(z_mid, i);

    // Calculate exner
    PF::exner_function<Spack>(team, p_mid_i, exner_i);

    // Calculate z_mid
    PF::calculate_dz(team, pseudo_density_i, p_mid_i, T_mid_i, qv_i, dz_i);
    const Real z_surf = 0.0; // For now, set z_int(i,nlevs) = z_surf = 0
    team.team_barrier();
    PF::calculate_z_int(team, nlevs, dz_i, z_surf, z_int_i);
    team.team_barrier();
    PF::calculate_z_mid(team, nlevs, z_int_i, z_mid_i);
  });

  // Compute TMS
  TMSFunctions::compute_tms(ncols, nlevs,
                            ekat::scalarize(horiz_winds),
                            ekat::scalarize(T_mid),
                            ekat::scalarize(p_mid),
                            ekat::scalarize(exner),
                            ekat::scalarize(z_mid),
                            sgh30, landfrac,
                            surf_drag_coeff_tms, wind_stress_tms);
}

// =========================================================================================
void TurbulentMountainStress::finalize_impl()
{
  // Do nothing
}
// =========================================================================================
size_t TurbulentMountainStress::requested_buffer_size_in_bytes() const
{
  const int nlev_packs  = ekat::npack<Spack>(m_nlevs);
  const int nlevi_packs = ekat::npack<Spack>(m_nlevs+1);
  return Buffer::num_2d_midpoint_views*m_ncols*nlev_packs*sizeof(Spack) +
         Buffer::num_2d_interface_views*m_ncols*nlevi_packs*sizeof(Spack);
}
// =========================================================================================
void TurbulentMountainStress::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
                   "Error! Buffers size not sufficient.\n");

  Spack* mem = reinterpret_cast<Spack*>(buffer_manager.get_memory());
  const int nlev_packs  = ekat::npack<Spack>(m_nlevs);
  const int nlevi_packs = ekat::npack<Spack>(m_nlevs+1);

  uview_2d* buffer_mid_view_ptrs[Buffer::num_2d_midpoint_views] = {
    &m_buffer.exner, &m_buffer.dz, &m_buffer.z_mid
  };

  for (int i=0; i<Buffer::num_2d_midpoint_views; ++i) {
    *buffer_mid_view_ptrs[i] = uview_2d(mem, m_ncols, nlev_packs);
    mem += buffer_mid_view_ptrs[i]->size();
  }

  uview_2d* buffer_int_view_ptrs[Buffer::num_2d_interface_views] = {
    &m_buffer.z_int
  };

  for (int i=0; i<Buffer::num_2d_interface_views; ++i) {
    *buffer_int_view_ptrs[i] = uview_2d(mem, m_ncols, nlevi_packs);
    mem += buffer_int_view_ptrs[i]->size();
  }

  size_t used_mem = (reinterpret_cast<Real*>(mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem == requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for TurbulentMountainStress.");
}
} // namespace scream
