#include "field_prev.hpp"

#include "share/util/eamxx_universal_constants.hpp"

namespace scream {

FieldPrevDiag::
FieldPrevDiag(const ekat::Comm &comm,
              const ekat::ParameterList &params)
 : AtmosphereDiagnostic(comm, params)
{
  EKAT_REQUIRE_MSG(params.isParameter("field_name"),
                   "Error! FieldPrevDiag requires 'field_name' in its "
                   "input parameters.\n");

  m_name = m_params.get<std::string>("field_name");
}

void FieldPrevDiag::
create_requests()
{
  const auto &gname = m_params.get<std::string>("grid_name");
  add_field<Required>(m_name, gname);
}

void FieldPrevDiag::initialize_impl(const RunType /*run_type*/) {
  const auto &f   = get_field_in(m_name);
  const auto &fid = f.get_header().get_identifier();
  const auto &gn  = fid.get_grid_name();

  EKAT_REQUIRE_MSG(
      f.data_type() == DataType::RealType,
      "Error! FieldPrevDiag only supports Real data type fields.\n"
      " - field name: " +
          fid.name() +
          "\n"
          " - field data type: " +
          e2str(f.data_type()) + "\n");

  const auto &layout = fid.get_layout();

  // Output has the same layout and units as the input field
  FieldIdentifier d_fid(m_name + "_prev", layout.clone(), fid.get_units(), gn);
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  // Storage for the field value at the start of the timestep
  FieldIdentifier prev_fid(m_name + "_prev_store", layout.clone(), fid.get_units(), gn);
  m_f_prev = Field(prev_fid);
  m_f_prev.allocate_view();
}

void FieldPrevDiag::init_timestep(const util::TimeStamp &start_of_step) {
  const auto &f_curr = get_field_in(m_name);
  m_f_prev.deep_copy(f_curr);
  m_f_prev.get_header().get_tracking().update_time_stamp(start_of_step);
}

void FieldPrevDiag::compute_diagnostic_impl() {
  const auto &prev_ts = m_f_prev.get_header().get_tracking().get_time_stamp();

  if (prev_ts.is_valid()) {
    // We have a stored value from the start of the previous timestep
    m_diagnostic_output.deep_copy(m_f_prev);
  } else {
    // init_timestep has not been called yet; fill with invalid sentinel
    m_diagnostic_output.deep_copy(constants::fill_value<Real>);
  }
}

}  // namespace scream
