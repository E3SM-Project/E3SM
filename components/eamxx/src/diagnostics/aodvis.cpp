#include "diagnostics/aodvis.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

namespace scream {

AODVis::AODVis(const ekat::Comm &comm, const ekat::ParameterList &params)
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
  // TODO: don't hardcode this!
  m_swbnd = 14;

  // Define layouts we need (both inputs and outputs)
  FieldLayout scalar3d_swband_layout{{COL, SWBND, LEV},
                                     {m_ncols, m_swbnd, m_nlevs}};
  FieldLayout scalar1d_layout{{COL}, {m_ncols}};

  // The fields required for this diagnostic to be computed
  add_field<Required>("aero_tau_sw", scalar3d_swband_layout, nondim, grid_name);

  // Construct and allocate the aodvis field
  // We are going to assume we have nondim units here for ease
  FieldIdentifier fid("AerosolOpticalDepth550nm", scalar1d_layout, nondim,
                      grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void AODVis::compute_diagnostic_impl() {
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  const auto aod = m_diagnostic_output.get_view<Real *>();
  // TODO: don't hardcode swbnd 10
  // Get slice of tau at swbnd 10
  const auto tau_vis =
      get_field_in("aero_tau_sw").subfield(1, 10).get_view<Real **>();

  const auto num_levs = m_nlevs;
  const auto policy   = ESU::get_default_team_policy(m_ncols, m_nlevs);
  Kokkos::parallel_for(
      "Compute " + name(), policy, KOKKOS_LAMBDA(const MT &team) {
        const int icol = team.league_rank();
        auto tau_icol  = ekat::subview(tau_vis, icol);
        Kokkos::parallel_reduce(
            Kokkos::TeamVectorRange(team, num_levs),
            [&](const int &ilev, Real &lsum) { lsum += tau_icol(ilev); },
            aod(icol));
        team.team_barrier();
      });
}

}  // namespace scream
