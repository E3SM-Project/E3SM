#include "field_over_dt.hpp"

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
  m_diagnostic_output = Field(d_fid,true);
  if (f.has_valid_mask()) {
    m_diagnostic_output.set_valid_mask(f.get_valid_mask());
    m_diagnostic_output.get_header().set_may_be_filled(true);
  }
}

void FieldOverDtDiag::init_timestep(const util::TimeStamp &start_of_step) {
  m_start_ts = start_of_step;
  EKAT_REQUIRE_MSG (m_start_ts.is_valid(),
      "Error! Initializing FieldOverDtDiag timestep with an invalid time stamp.\n"
      " - diag field name: " + m_diagnostic_output.name() + "\n");
}

void FieldOverDtDiag::compute_diagnostic_impl()
{
  const auto &f = get_field_in(m_name);

  const auto &curr_ts = f.get_header().get_tracking().get_time_stamp();
  EKAT_REQUIRE_MSG (curr_ts.is_valid() and m_start_ts.is_valid(),
      "Error! FieldOverDtDiag does not work if you don't call init_timestep first.\n"
      " - diag field name: " + m_diagnostic_output.name() + "\n");

  const std::int64_t dt = curr_ts.seconds_from(m_start_ts);

  EKAT_REQUIRE_MSG(dt > 0,
      "Error! FieldOverDtDiag: dt must be positive.\n"
      " - field name: " + m_name + "\n"
      " - start timestamp: " + m_start_ts.to_string() + "\n"
      " - curr timestamp:  " + curr_ts.to_string() + "\n");

  // diag = 0*diag + 1/dt*f
  const auto dt_inv = 1 / static_cast<Real>(dt);
  if (f.has_valid_mask()) {
    m_diagnostic_output.update(f,dt_inv,0,f.get_valid_mask());
    // TODO: remove when IO handles fill value internally
    m_diagnostic_output.deep_copy(constants::fill_value<Real>,f.get_valid_mask(),true);
  } else {
    m_diagnostic_output.update(f,dt_inv,0);
  }
}

}  // namespace scream
