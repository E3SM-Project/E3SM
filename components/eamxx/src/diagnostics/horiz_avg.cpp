#include "diagnostics/horiz_avg.hpp"

#include "share/field/field_utils.hpp"

namespace scream {

HorizAvgDiag::HorizAvgDiag(const ekat::Comm &comm,
                           const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  const auto &fname = m_params.get<std::string>("field_name");
  m_diag_name       = fname + "_horiz_avg";
}

void HorizAvgDiag::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  const auto &fn = m_params.get<std::string>("field_name");
  const auto &gn = m_params.get<std::string>("grid_name");
  const auto g   = grids_manager->get_grid("Physics");

  add_field<Required>(fn, gn);

  // first clone the area unscaled, we will scale it later in initialize_impl
  m_scaled_area = g->get_geometry_data("area").clone();
}

void HorizAvgDiag::initialize_impl(const RunType /*run_type*/) {
  using namespace ShortFieldTagsNames;
  const auto &f      = get_fields_in().front();
  const auto &fid    = f.get_header().get_identifier();
  const auto &layout = fid.get_layout();

  EKAT_REQUIRE_MSG(layout.rank() >= 1 && layout.rank() <= 3,
                   "Error! Field rank not supported by HorizAvgDiag.\n"
                   " - field name: " +
                       fid.name() +
                       "\n"
                       " - field layout: " +
                       layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(layout.tags()[0] == COL,
                   "Error! HorizAvgDiag diagnostic expects a layout starting "
                   "with the 'COL' tag.\n"
                   " - field name  : " +
                       fid.name() +
                       "\n"
                       " - field layout: " +
                       layout.to_string() + "\n");

  FieldIdentifier d_fid(m_diag_name, layout.clone().strip_dim(COL),
                        fid.get_units(), fid.get_grid_name());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  // scale the area field
  auto total_area = field_sum<Real>(m_scaled_area, &m_comm);
  m_scaled_area.scale(sp(1.0) / total_area);
}

void HorizAvgDiag::compute_diagnostic_impl() {
  const auto &f = get_fields_in().front();
  const auto &d = m_diagnostic_output;
  // Call the horiz_contraction impl that will take care of everything
  horiz_contraction<Real>(d, f, m_scaled_area, &m_comm);
}

}  // namespace scream
