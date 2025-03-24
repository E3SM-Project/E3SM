#include "diagnostics/aodvis.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

#include "share/util/eamxx_universal_constants.hpp"

namespace scream {

AODVis::AODVis(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  // Nothing to do here
}

void AODVis::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  auto grid             = grids_manager->get_grid("Physics");
  const auto &grid_name = grid->name();

  const auto nondim = Units::nondimensional();

  m_ncols = grid->get_num_local_dofs();
  m_nlevs = grid->get_num_vertical_levels();

  // Define layouts we need (both inputs and outputs)
  auto vector3d = grid->get_3d_vector_layout(true, m_swbands, "swband");
  auto scalar2d = grid->get_2d_scalar_layout();

  // The fields required for this diagnostic to be computed
  add_field<Required>("aero_tau_sw", vector3d, nondim, grid_name);
  add_field<Required>("sunlit",      scalar2d, nondim, grid_name);

  // Construct and allocate the aodvis field
  FieldIdentifier fid(name(), scalar2d, nondim, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void AODVis::compute_diagnostic_impl() {
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  Real var_fill_value = constants::DefaultFillValue<Real>().value;

  const auto aod     = m_diagnostic_output.get_view<Real *>();
  const auto tau_vis = get_field_in("aero_tau_sw")
                           .subfield(1, m_vis_bnd)
                           .get_view<const Real **>();
  const auto sunlit = get_field_in("sunlit").get_view<const Real *>();

  const auto num_levs = m_nlevs;
  const auto policy   = ESU::get_default_team_policy(m_ncols, m_nlevs);
  Kokkos::parallel_for(
      "Compute " + m_diagnostic_output.name(), policy, KOKKOS_LAMBDA(const MT &team) {
        const int icol = team.league_rank();
        if(sunlit(icol) == 0.0) {
          aod(icol) = var_fill_value;
        } else {
          auto tau_icol = ekat::subview(tau_vis, icol);
          aod(icol)     = ESU::view_reduction(team, 0, num_levs, tau_icol);
        }
      });
}

}  // namespace scream
