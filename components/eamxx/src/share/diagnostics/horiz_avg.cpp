#include "horiz_avg.hpp"

#include "share/field/field_utils.hpp"

namespace scream {

HorizAvgDiag::
HorizAvgDiag(const ekat::Comm &comm,
             const ekat::ParameterList &params)
 : AtmosphereDiagnostic(comm, params)
{
  EKAT_REQUIRE_MSG (m_params.isParameter("field_name"),
      "[HorizAvgDiag] Error! Missing required param 'field_name'\n");
}

void HorizAvgDiag::create_requests()
{
  const auto &fn = m_params.get<std::string>("field_name");
  const auto &gn = m_params.get<std::string>("grid_name");
  const auto g   = m_grids_manager->get_grid("physics");

  add_field<Required>(fn, gn);

  m_area = g->get_geometry_data("area");
}

void HorizAvgDiag::initialize_impl(const RunType /*run_type*/)
{
  using namespace ShortFieldTagsNames;

  const auto &f      = get_fields_in().front();
  const auto &fid    = f.get_header().get_identifier();
  const auto &layout = fid.get_layout();

  EKAT_REQUIRE_MSG(layout.rank() >= 1 && layout.rank() <= 3,
      "Error! Field rank not supported by HorizAvgDiag.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(layout.tags()[0] == COL,
      "Error! HorizAvgDiag diagnostic expects a layout starting with the 'COL' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  const auto &fname = m_params.get<std::string>("field_name");
  const auto dname  = fname + "_horiz_avg";

  auto d_fid = fid.clone(dname).reset_layout(layout.clone().strip_dim(COL));
  m_diagnostic_output = Field(d_fid,true);

  if (f.has_valid_mask()) {
    // Output valid_mask: 1 where den > 0 (i.e., at least one valid column)
    m_diagnostic_output.create_valid_mask(Field::MaskInit::None);

    // To use at runtime to compute sum(1*w) via horiz_avg (with f_in=1)
    m_ones = f.clone();
    m_ones.deep_copy(1);
    m_ones.set_valid_mask(f.get_valid_mask());

    m_denom = m_diagnostic_output.clone("denom");
  } else {
    // Since area is constant and there is no masking, we can pre-compute
    // the area sum, and set weight = area/sum(area), so that we avoid
    // one loop at runtime.
    // NOTE: clone m_area, since the field from the grid is read-only
    m_area = m_area.clone("scaled_area");

    m_denom = m_area.subfield(COL,0).clone(); // This should create a 0d field
    m_ones = m_area.clone("ones");
    m_ones.deep_copy(1);

    horiz_contraction(m_denom, m_ones, m_area, m_comm);
    m_denom.sync_to_host();
    auto d = m_denom.get_internal_view_data<Real,Host>()[0];
    m_area.scale(1/d);

    // Free up the memory
    m_ones = Field();
    m_denom = Field();
  }
}

void HorizAvgDiag::compute_diagnostic_impl()
{
  const auto &f = get_fields_in().front();

  // sum(w * f) possibly masked
  horiz_contraction(m_diagnostic_output, f, m_area, m_comm);

  if (f.has_valid_mask()) {
    // Denominator: sum(weight) masked
    horiz_contraction(m_denom, m_ones, m_area, m_comm);

    auto& nonzero_denom = m_diagnostic_output.get_valid_mask();
    compute_mask(m_denom,0,Comparison::NE,nonzero_denom);

    m_diagnostic_output.scale_inv(m_denom,nonzero_denom);
  }
}

}  // namespace scream
