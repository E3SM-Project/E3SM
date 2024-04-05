#include "diagnostics/aodvis.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

namespace scream {

AODVis::AODVis(const ekat::Comm &comm,
                                     const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  // Nothing to do here
}

void AODVis::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto grid             = grids_manager->get_grid("Physics");
  const auto &grid_name = grid->name();

  const auto nondim = Units::nondimensional();

  m_ncols = grid->get_num_local_dofs();
  m_nlevs = grid->get_num_vertical_levels();

  // Define layouts we need (both inputs and outputs)
  // TODO: define this properly (dunno state of things)
  auto scalar3d_swband_layout = grid->scalar3d_swband_layout(true);
  auto scalar2d = grid->get_scalar2d(true);

  // The fields required for this diagnostic to be computed
  add_field<Required>("aero_tau_sw", scalar3d_swband_layout, nondim, grid_name);

  // Construct and allocate the aodvis field
  // We are going to assume we have nondim units here for ease
  FieldIdentifier fid("VisibileAerosolOpticalDepth", scalar2d, nondim, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void AODVis::compute_diagnostic_impl()
{
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  const auto aod = m_diagnostic_output.get_view<Real*>();
  const auto tau = get_field_in("aero_tau_sw").get_view<const Real***>();

  const auto num_levs = m_num_levs;
  const auto policy = ESU::get_default_team_policy(m_num_cols, m_num_levs);
  Kokkos::parallel_for("Compute " + name(), policy,
                       KOKKOS_LAMBDA(const MT& team) {
    const int icol = team.league_rank();
    // TODO: take a slice of tau at 9th band (vis band)
    auto tau_vis = tau(:, 9, :) 
    auto tau_icol = ekat::subview(tau_vis,icol);
    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, num_levs),
                            [&] (const int& ilev, Real& lsum) {
      lsum += tau_icol(ilev);
    },aod(icol));
    team.team_barrier();
  });
}

}  // namespace scream
