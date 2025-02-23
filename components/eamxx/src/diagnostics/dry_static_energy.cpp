#include "diagnostics/dry_static_energy.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

namespace scream
{

// =========================================================================================
DryStaticEnergyDiagnostic::
DryStaticEnergyDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void DryStaticEnergyDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  auto m2  = pow(m,2);
  auto s2  = pow(s,2);

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  auto scalar2d = grid->get_2d_scalar_layout();
  auto scalar3d = grid->get_3d_scalar_layout(true);

  // The fields required for this diagnostic to be computed
  add_field<Required>("T_mid",          scalar3d, K,      grid_name);
  add_field<Required>("pseudo_density", scalar3d, Pa,     grid_name);
  add_field<Required>("p_mid",          scalar3d, Pa,     grid_name);
  add_field<Required>("qv",             scalar3d, kg/kg,  grid_name);
  add_field<Required>("phis",           scalar2d, m2/s2,  grid_name);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar3d, m2/s2, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();

  // Initialize a 2d view of dz to be used in compute_diagnostic
  m_tmp_mid = view_2d("",m_num_cols,m_num_levs);
  m_tmp_int = view_2d("",m_num_cols,m_num_levs+1);
}
// =========================================================================================
void DryStaticEnergyDiagnostic::compute_diagnostic_impl()
{
  using MemberType = typename KT::MemberType;
  using PF = PhysicsFunctions<DefaultDevice>;

  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(m_num_cols, m_num_levs);

  const auto& dse                = m_diagnostic_output.get_view<Real**>();
  const auto& T_mid              = get_field_in("T_mid").get_view<const Real**>();
  const auto& p_mid              = get_field_in("p_mid").get_view<const Real**>();
  const auto& qv_mid             = get_field_in("qv").get_view<const Real**>();
  const auto& pseudo_density_mid = get_field_in("pseudo_density").get_view<const Real**>();
  const auto& phis               = get_field_in("phis").get_view<const Real*>();

  // Set surface geopotential for this diagnostic
  const Real surf_geopotential = 0.0;

  const int num_levs = m_num_levs;
  auto      tmp_mid  = m_tmp_mid;
  auto      tmp_int  = m_tmp_int;

  Kokkos::parallel_for("DryStaticEnergyDiagnostic",
                       default_policy,
                       KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank();
    const auto& dz_s    = ekat::subview(tmp_mid,icol);
    const auto& z_int_s = ekat::subview(tmp_int,icol);
    const auto& z_mid_s = dz_s; // Reuse the memory for z_mid, but set a new variable for code readability.
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, num_levs), [&] (const Int& ilev) {
      dz_s(ilev) = PF::calculate_dz(pseudo_density_mid(icol,ilev), p_mid(icol,ilev), T_mid(icol,ilev), qv_mid(icol,ilev));
    });
    team.team_barrier();

    PF::calculate_z_int(team,num_levs,dz_s,surf_geopotential,z_int_s);
    team.team_barrier();

    PF::calculate_z_mid(team,num_levs,z_int_s,z_mid_s);
    team.team_barrier();

    const auto& dse_s = ekat::subview(dse,icol);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, num_levs), [&] (const Int& ilev) {
      dse_s(ilev) = PF::calculate_dse(T_mid(icol,ilev),z_mid_s(ilev),phis(icol));
    });
    team.team_barrier();
  });
  Kokkos::fence();
}

} //namespace scream
