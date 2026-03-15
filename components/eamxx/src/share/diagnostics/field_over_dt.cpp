#include "field_over_dt.hpp"

#include "share/util/eamxx_universal_constants.hpp"

namespace scream {

FieldOverDtDiag::
FieldOverDtDiag(const ekat::Comm &comm,
                const ekat::ParameterList &params)
 : AtmosphereDiagnostic(comm, params)
{
  EKAT_REQUIRE_MSG(params.isParameter("field_name"),
                   "Error! FieldOverDtDiag requires 'field_name' in its "
                   "input parameters.\n");

  m_name = m_params.get<std::string>("field_name");
}

void FieldOverDtDiag::
create_requests()
{
  const auto &gname = m_params.get<std::string>("grid_name");
  add_field<Required>(m_name, gname);
}

void FieldOverDtDiag::initialize_impl(const RunType /*run_type*/) {
  const auto &f   = get_field_in(m_name);
  const auto &fid = f.get_header().get_identifier();
  const auto &gn  = fid.get_grid_name();

  EKAT_REQUIRE_MSG(
      f.data_type() == DataType::RealType,
      "Error! FieldOverDtDiag only supports Real data type fields.\n"
      " - field name: " + fid.name() + "\n"
      " - field data type: " + e2str(f.data_type()) + "\n");

  const auto &layout = fid.get_layout();

  using namespace ekat::units;
  auto diag_units = fid.get_units() / s;

  FieldIdentifier d_fid(m_name + "_over_dt", layout.clone(), diag_units, gn);
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();
}

void FieldOverDtDiag::init_timestep(const util::TimeStamp &start_of_step) {
  m_start_ts = start_of_step;
}

void FieldOverDtDiag::compute_diagnostic_impl() {
  const auto &f = get_field_in(m_name);

  if (!m_start_ts.is_valid()) {
    // init_timestep has not been called yet; fill with invalid sentinel
    m_diagnostic_output.deep_copy(constants::fill_value<Real>);
    return;
  }

  const auto &curr_ts = f.get_header().get_tracking().get_time_stamp();
  const std::int64_t dt = curr_ts - m_start_ts;

  EKAT_REQUIRE_MSG(dt > 0,
      "Error! FieldOverDtDiag: dt must be positive.\n"
      " - field name: " + m_name + "\n"
      " - start timestamp: " + m_start_ts.to_string() + "\n"
      " - curr timestamp:  " + curr_ts.to_string() + "\n");

  m_diagnostic_output.deep_copy(f);
  m_diagnostic_output.scale(1.0 / dt);
}

}  // namespace scream
