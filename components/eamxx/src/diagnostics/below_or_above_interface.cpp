#include "diagnostics/below_or_above_interface.hpp"

#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_team_policy_utils.hpp>
#include <ekat_reduction_utils.hpp>

namespace scream {

BelowOrAboveInterface::BelowOrAboveInterface(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  m_name = m_params.get<std::string>("input_field");
  m_type = m_params.get<std::string>("above_or_below");
}

void BelowOrAboveInterface::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  const auto &gname = m_params.get<std::string>("grid_name");
  add_field<Required>(m_name, gname);
}

void AODVis::initialize_impl(const RunType /*run_type*/) {
  const auto &f   = get_field_in(m_name);
  const auto &fid = f.get_header().get_identifier();
  const auto &gn  = fid.get_grid_name();
  const auto &layout = fid.get_layout();
  
  // ensure last dim in layout is ILEV
  using namespace ShortFieldTagsNames;
  EKAT_REQUIRE_MSG(layout.dim( layout.rank()-1) == ILEV,
      "Error! BelowOrAboveInterface diagnostic requires the last dimension of the input field to be ILEV.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n"); 
  auto scalar3d_mid = grid->get_3d_scalar_layout(true);

  // Construct and allocate the aodvis field
  FieldIdentifier fid(m_name + "_" + m_type, scalar3d_mid, fid.get_units(), gn);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void AODVis::compute_diagnostic_impl() {
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using TPF = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
  using RU  = ekat::ReductionUtils<typename KT::ExeSpace>;

  const auto f = get_field_in(m_name).get_view<const Real **>();
  const auto d = m_diagnostic_output.get_view<Real **>();
  const auto num_cols = get_field_in(m_name).get_header().get_num_local_dofs();
  const auto num_levs = get_field_in(m_name).get_header().get_num_vertical_levels();
  // note -1 to loop over LEVs, not ILEVs
  const auto policy   = TPF::get_default_team_policy(num_cols, num_levs-1);

  auto d_type = m_type;
  Kokkos::parallel_for("BelowOrAboveInterface", policy, KOKKOS_LAMBDA(const MT& team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_levels), [&] (const int ilev) {
      if (d_type == "above") {
        // interface 0 is above midpoint 0
        d(icol, ilev) = f(icol, ilev);
      } else {
        // interface 1 is below midpoint 0
        d(icol, ilev) = f(icol, ilev+1);
      }
    });
  });
}

}  // namespace scream
