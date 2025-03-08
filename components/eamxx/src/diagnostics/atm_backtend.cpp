#include "diagnostics/atm_backtend.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

#include "share/util/eamxx_universal_constants.hpp"

namespace scream {

AtmBackTendDiag::
AtmBackTendDiag(const ekat::Comm &comm,
                const ekat::ParameterList &params)
 : AtmosphereDiagnostic(comm, params)
{
  EKAT_REQUIRE_MSG(params.isParameter("Tendency Name"),
                   "Error! AtmBackTendDiag requires 'Tendency Name' in its "
                   "input parameters.\n");

  m_name = m_params.get<std::string>("Tendency Name");
}

void AtmBackTendDiag::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  const auto &gname = m_params.get<std::string>("grid_name");
  add_field<Required>(m_name, gname);
}

void AtmBackTendDiag::initialize_impl(const RunType /*run_type*/) {
  const auto &f   = get_field_in(m_name);
  const auto &fid = f.get_header().get_identifier();
  const auto &gn  = fid.get_grid_name();

  // Sanity checks
  const auto &layout = fid.get_layout();
  EKAT_REQUIRE_MSG(
      f.data_type() == DataType::RealType,
      "Error! AtmBackTendDiag only supports Real data type field.\n"
      " - field name: " +
          fid.name() +
          "\n"
          " - field data type: " +
          e2str(f.data_type()) + "\n");

  using namespace ekat::units;
  // The units are the same except per second
  auto diag_units = fid.get_units() / s;

  // All good, create the diag output
  FieldIdentifier d_fid(m_name + "_atm_backtend", layout.clone(), diag_units, gn);
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  // Let's also create the previous field
  FieldIdentifier prev_fid(m_name + "_atm_backtend_prev", layout.clone(), diag_units, gn);
  m_f_prev = Field(prev_fid);
  m_f_prev.allocate_view();
}

void AtmBackTendDiag::init_timestep(const util::TimeStamp &start_of_step) {
  const auto &f_curr = get_field_in(m_name);
  m_f_prev.deep_copy(f_curr);
  m_f_prev.get_header().get_tracking().update_time_stamp(start_of_step);
}

void AtmBackTendDiag::compute_diagnostic_impl() {
  Real var_fill_value = constants::DefaultFillValue<Real>().value;
  std::int64_t dt;

  const auto &f       = get_field_in(m_name);
  const auto &curr_ts = f.get_header().get_tracking().get_time_stamp();
  const auto &prev_ts = m_f_prev.get_header().get_tracking().get_time_stamp();

  if(prev_ts.is_valid()) {
    // This diag was called before, so we have a valid value for m_f_prev,
    // and can compute the tendency
    dt = curr_ts - prev_ts;
    m_f_prev.update(f, 1.0 / dt, -1.0 / dt);
    m_diagnostic_output.deep_copy(m_f_prev);
  } else {
    // This is the first time we evaluate this diag. We cannot compute a tend
    // yet, so fill with an invalid value
    m_diagnostic_output.deep_copy(var_fill_value);
  }
}

}  // namespace scream
