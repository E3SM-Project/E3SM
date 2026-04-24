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

void FieldPrevDiag::initialize_impl() {
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
  const auto &f_curr_ts = f_curr.get_header().get_tracking().get_time_stamp();

  // Only capture f_curr if it has been initialized (valid timestamp).
  // If the source is a computed diagnostic that has not been evaluated yet
  // (e.g., on the very first timestep before any computation has run),
  // its data is Kokkos-zero-initialized — not a meaningful "t=0" value.
  // Leaving m_f_prev with an invalid timestamp causes compute_diagnostic_impl
  // to return fill_value, which is safer than silently propagating zeros into
  // downstream arithmetic (e.g., X_minus_X_prev_times_something can overflow
  // if X_prev is zero and X is fill_value from another uninitialized diag).
  if (f_curr_ts.is_valid()) {
    m_f_prev.deep_copy(f_curr);
    m_f_prev.get_header().get_tracking().update_time_stamp(start_of_step);
  }
}

void FieldPrevDiag::compute_diagnostic_impl() {
  const auto &prev_ts = m_f_prev.get_header().get_tracking().get_time_stamp();

  if (prev_ts.is_valid()) {
    // Normal case: return the value captured at the start of this timestep.
    m_diagnostic_output.deep_copy(m_f_prev);
  } else {
    // m_f_prev was never captured because the source field had no valid
    // timestamp when init_timestep ran (e.g., a derived diagnostic on the
    // very first step).  By the time compute_diagnostic_impl is called the
    // source has always been computed (the base-class timestamp guard
    // ensures this), so we use the current source value as a stand-in for
    // "X at t=0".  This gives X − X_prev = 0 on the first step and avoids
    // fill_value propagation into downstream arithmetic that could overflow
    // (fill_value × fill_value → Inf).
    m_diagnostic_output.deep_copy(get_field_in(m_name));
  }
}

}  // namespace scream
