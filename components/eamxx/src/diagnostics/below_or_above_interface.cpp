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

void BelowOrAboveInterface::initialize_impl(const RunType /*run_type*/) {
  const auto &f   = get_field_in(m_name);
  const auto &fid = f.get_header().get_identifier();
  const auto &gn  = fid.get_grid_name();
  const auto &layout = fid.get_layout();
  
  // // ensure last dim in layout is ILEV
  // using namespace ShortFieldTagsNames;
  // EKAT_REQUIRE_MSG(layout.dim( layout.rank()-1) == ILEV,
  //     "Error! BelowOrAboveInterface diagnostic requires the last dimension of the input field to be ILEV.\n"
  //     " - field name: " + fid.name() + "\n"
  //     " - field layout: " + layout.to_string() + "\n"); 
  
  using namespace ShortFieldTagsNames;
  FieldLayout diag_layout = {{COL, LEV}, {layout.dims()[0], layout.dims()[1]-1}};

  // Construct and allocate the output field
  FieldIdentifier fid_out(m_name + "_" + m_type, diag_layout, fid.get_units(), gn);
  m_diagnostic_output = Field(fid_out);
  m_diagnostic_output.allocate_view();
}

void BelowOrAboveInterface::compute_diagnostic_impl() {
  using KT  = KokkosTypes<DefaultDevice>;
  using MT  = typename KT::MemberType;
  using TPF = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
  using RU  = ekat::ReductionUtils<typename KT::ExeSpace>;

  const auto f = get_field_in(m_name).get_view<const Real **>();
  const auto d = m_diagnostic_output.get_view<Real **>();
  const auto& layout = get_field_in(m_name).get_header().get_identifier().get_layout();
  const auto num_cols = layout.dims()[0];
  const auto num_levs = layout.dims()[1];
  // note -1 to loop over LEVs, not ILEVs
  const auto policy   = TPF::get_default_team_policy(num_cols, num_levs-1);

  // diag_offset is 1 if m_type is below, otherwise it is 0
  auto diag_offset = (m_type == "below")? 1 : 0;
  Kokkos::parallel_for("BelowOrAboveInterface", policy, KOKKOS_LAMBDA(const MT& team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_levs-1), [&] (const int ilev) {
      d(icol, ilev) = f(icol, ilev+diag_offset);
    });
  });
}

}  // namespace scream
