#include "horiz_avg.hpp"

#include "share/field/field_utils.hpp"
#include "share/util/eamxx_universal_constants.hpp"

namespace scream {

HorizAvg::
HorizAvg(const ekat::Comm &comm,
             const ekat::ParameterList &params,
             const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm, params, grid)
{
  EKAT_REQUIRE_MSG (m_params.isParameter("field_name"),
      "[HorizAvg] Error! Missing required param 'field_name'\n");

  using namespace ShortFieldTagsNames;

  const auto &fn = m_params.get<std::string>("field_name");
  m_field_in_names.push_back(fn);

  // The weight: the grid's area (default), 1, or a field (as in X.weighted(W))
  m_weight_name = m_params.get<std::string>("weight_name","");
  if (not m_weight_name.empty()) {
    EKAT_REQUIRE_MSG (not m_params.isParameter("weighting_method"),
        "[HorizAvg] Error! Use either weighting_method or weight_name, not both.\n");
  }
  const auto method = m_params.get<std::string>("weighting_method","area");
  EKAT_REQUIRE_MSG (method=="area" or method=="none",
      "[HorizAvg] Error! Invalid weighting_method.\n"
      " - weighting_method: " + method + "\n"
      " - valid methods: 'area', 'none'\n");

  if (m_weight_name=="area" or (m_weight_name.empty() and method=="area")) {
    // The cell area is grid geometry, not a field
    m_area = m_grid->get_geometry_data("area");
  } else if (m_weight_name.empty()) {
    FieldIdentifier ones_fid("ones",m_grid->get_2d_scalar_layout(),ekat::units::none,m_grid->name());
    m_area = Field(ones_fid,true);
    m_area.deep_copy(1);
  } else {
    m_weight_is_field = true;
    // X.weighted(X) is allowed; the field is then an input once
    if (m_weight_name!=fn) {
      m_field_in_names.push_back(m_weight_name);
    }
  }
}

void HorizAvg::initialize_impl()
{
  using namespace ShortFieldTagsNames;

  const auto &f      = m_fields_in.at(m_field_in_names.front());
  const auto &fid    = f.get_header().get_identifier();
  const auto &layout = fid.get_layout();

  EKAT_REQUIRE_MSG(layout.rank() >= 1 && layout.rank() <= 3,
      "Error! Field rank not supported by HorizAvg.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(layout.tags()[0] == COL,
      "Error! HorizAvg diagnostic expects a layout starting with the 'COL' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  if (m_weight_is_field) {
    const auto& w   = m_fields_in.at(m_weight_name);
    const auto& wid = w.get_header().get_identifier();
    EKAT_REQUIRE_MSG (wid.get_layout().rank()==1 and wid.get_layout().tags()[0]==COL and
                      w.data_type()==f.data_type(),
        "Error! HorizAvg weight must be a (COL) field with the data type of the averaged field.\n"
        " - field name   : " + fid.name() + "\n"
        " - weight name  : " + wid.name() + "\n"
        " - weight layout: " + wid.get_layout().to_string() + "\n");
    EKAT_REQUIRE_MSG (not w.has_valid_mask(),
        "Error! HorizAvg cannot use a masked field as weight.\n"
        "  As in xarray, weights cannot be missing; mask the averaged field instead.\n"
        " - field name : " + fid.name() + "\n"
        " - weight name: " + wid.name() + "\n");
    // Read in place: the weight is another field, updated before this diag runs
    m_area = w;
  }

  const auto &fname = m_params.get<std::string>("field_name");
  const auto dname  = fname + "_horiz_avg";

  auto d_fid = fid.clone(dname).reset_layout(layout.clone().strip_dim(COL));
  m_diagnostic_output = Field(d_fid,true);

  if (f.has_valid_mask() or m_weight_is_field) {
    // Output valid_mask: 1 where den != 0 (e.g., at least one valid column)
    m_diagnostic_output.create_valid_mask(Field::MaskInit::None);
    m_diagnostic_output.get_header().set_may_be_filled(true);

    // To use at runtime to compute sum(1*w) via horiz_avg (with f_in=1)
    m_ones = f.clone();
    m_ones.deep_copy(1);
    if (f.has_valid_mask()) {
      m_ones.set_valid_mask(f.get_valid_mask());
    }

    m_denom = m_diagnostic_output.clone("denom");
  } else {
    // Since the weight is constant and there is no masking, we can pre-compute
    // the area sum, and set weight = area/sum(area), so that we avoid
    // one loop at runtime.
    // NOTE: clone m_area, since the field from the grid is read-only
    m_area = m_area.clone("scaled_area", CloneFlags::CopyData);

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

void HorizAvg::compute_impl()
{
  const auto &f = m_fields_in.at(m_field_in_names.front());

  // sum(w * f) possibly masked
  horiz_contraction(m_diagnostic_output, f, m_area, m_comm);

  if (m_denom.is_allocated()) {
    // Denominator: sum(weight) masked
    horiz_contraction(m_denom, m_ones, m_area, m_comm);

    auto& nonzero_denom = m_diagnostic_output.get_valid_mask();
    compute_mask(m_denom,0,Comparison::NE,nonzero_denom);

    m_diagnostic_output.scale_inv(m_denom,nonzero_denom);
    // IO relies on fill_value for masked-out entries
    m_diagnostic_output.deep_copy(constants::fill_value<Real>,nonzero_denom,true);
  }
}

}  // namespace scream
