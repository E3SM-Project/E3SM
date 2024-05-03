#include "diagnostics/atm_tend.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

#include "share/util/scream_universal_constants.hpp"

namespace scream {

AtmTendDiag::AtmTendDiag(const ekat::Comm &comm,
                         const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  EKAT_REQUIRE_MSG(params.isParameter("Tend Name"),
                   "Error! AtmTendDiag requires 'Tend Name' in its "
                   "input parameters.\n");

  m_name = m_params.get<std::string>("Tend Name");
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

  // Sanity checks
  using namespace ShortFieldTagsNames;
  const auto &layout = fid.get_layout();
  EKAT_REQUIRE_MSG(f.data_type() == DataType::RealType,
                   "Error! FieldAtHeight only supports Real data type field.\n"
                   " - field name: " +
                       fid.name() +
                       "\n"
                       " - field data type: " +
                       e2str(f.data_type()) + "\n");

  // All good, create the diag output
  FieldIdentifier d_fid(name(), layout.clone(), fid.get_units(),
                        fid.get_grid_name());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  // Let's also create the previous field
  FieldIdentifier prev_fid(name() + "_prev", layout.clone(), fid.get_units(),
                           fid.get_grid_name());
  m_field_prev = Field(prev_fid);
  m_field_prev.allocate_view();
}
void AtmTendDiag::compute_diagnostic_impl() {
  Real var_fill_value = constants::DefaultFillValue<Real>().value;
  std::int64_t dt;
  auto tts = m_diagnostic_output.get_header().get_tracking().get_time_stamp();

  const auto &f = get_field_in(m_name);

  if(m_ts.is_valid()) {
    dt       = tts - m_ts;
    auto ddt = static_cast<Real>(dt);
    m_ts     = tts;
    m_field_prev.update(f, 1 / ddt, -1 / ddt);
    m_diagnostic_output.deep_copy(m_field_prev);
  } else {
    m_diagnostic_output.deep_copy(var_fill_value);
    m_ts = tts;
  }
  m_field_prev.deep_copy(f);
}

}  // namespace scream
