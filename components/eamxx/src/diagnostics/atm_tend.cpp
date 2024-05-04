#include "diagnostics/atm_tend.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

#include "share/util/scream_universal_constants.hpp"

namespace scream {

AtmTendDiag::AtmTendDiag(const ekat::Comm &comm,
                         const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  EKAT_REQUIRE_MSG(params.isParameter("Field Name"),
                   "Error! AtmTendDiag requires 'Field Name' in its "
                   "input parameters.\n");

  m_name = m_params.get<std::string>("Field Name");
}

std::string AtmTendDiag::name() const { return m_name + "_atm_tend"; }

void AtmTendDiag::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  const auto &gname = m_params.get<std::string>("grid_name");
  add_field<Required>(m_name, gname);
}

void AtmTendDiag::initialize_impl(const RunType /*run_type*/) {
  const auto &f   = get_field_in(m_name);
  const auto &fid = f.get_header().get_identifier();
  const auto &gn  = fid.get_grid_name();

  // Sanity checks
  using namespace ShortFieldTagsNames;
  const auto &layout = fid.get_layout();
  EKAT_REQUIRE_MSG(f.data_type() == DataType::RealType,
                   "Error! AtmTendDiag only supports Real data type field.\n"
                   " - field name: " +
                       fid.name() +
                       "\n"
                       " - field data type: " +
                       e2str(f.data_type()) + "\n");

  using namespace ekat::units;
  // The units are the same except per second
  auto diag_units = fid.get_units() / s;
  // TODO: set the units string correctly by appending "/s"

  // All good, create the diag output
  FieldIdentifier d_fid(name(), layout.clone(), diag_units, gn);
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  // Let's also create the previous field
  FieldIdentifier prev_fid(name() + "_prev", layout.clone(), diag_units, gn);
  m_f_prev = Field(prev_fid);
  m_f_prev.allocate_view();
}
void AtmTendDiag::compute_diagnostic_impl() {
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
  m_f_prev.deep_copy(f);
  m_f_prev.get_header().get_tracking().update_time_stamp(curr_ts);
}

}  // namespace scream
