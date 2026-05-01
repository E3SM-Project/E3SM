#include "field_prev.hpp"

#include "share/util/eamxx_universal_constants.hpp"

namespace scream {

FieldPrev::
FieldPrev(const ekat::Comm &comm,
              const ekat::ParameterList &params,
              const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm, params, grid)
{
  EKAT_REQUIRE_MSG(params.isParameter("field_name"),
                   "Error! FieldPrev requires 'field_name' in its "
                   "input parameters.\n");

  m_name = m_params.get<std::string>("field_name");
  m_field_in_names.push_back(m_name);
}

void FieldPrev::initialize_impl() {
  const auto &f   = m_fields_in.at(m_name);
  const auto &fid = f.get_header().get_identifier();
  const auto &gn  = fid.get_grid_name();

  EKAT_REQUIRE_MSG(
      f.data_type() == DataType::RealType,
      "Error! FieldPrev only supports Real data type fields.\n"
      " - field name: " +
          fid.name() +
          "\n"
          " - field data type: " +
          e2str(f.data_type()) + "\n");

  const auto &layout = fid.get_layout();

  // Output has the same layout and units as the input field
  FieldIdentifier d_fid(m_name + "_prev", layout.clone(), fid.get_units(), gn);
  m_diagnostic_output = Field(d_fid,true);

  // NOTE: even if the input field is masked, we must use a different one, since
  // at the first timestep, the field may have not been computed yet (but it's not masked per se)
  auto mask = m_diagnostic_output.create_valid_mask();
  m_diagnostic_output.get_header().set_may_be_filled(true);
  mask.deep_copy(0);

  // TODO: remove when IO stops relying on mask=0 entries being already set to FillValue
  m_diagnostic_output.deep_copy(constants::fill_value<Real>);
}

void FieldPrev::init_timestep(const util::TimeStamp &start_of_step)
{
  // At the timestep begin, the field is at the "prev" timestep.
  const auto& f_prev = m_fields_in.at(m_name);
  const auto& prev_ts = f_prev.get_header().get_tracking().get_time_stamp();

  if (prev_ts.is_valid()) {
    // Field was inited/computed at least once, so it makes sense to copy it
    if (f_prev.has_valid_mask()) {
      m_diagnostic_output.deep_copy(f_prev,f_prev.get_valid_mask());
      m_diagnostic_output.get_valid_mask().deep_copy(f_prev.get_valid_mask());
    } else {
      m_diagnostic_output.deep_copy(f_prev);
      m_diagnostic_output.get_valid_mask().deep_copy(1);
    }
  } else {
    m_diagnostic_output.get_valid_mask().deep_copy(0);
  }
}

void FieldPrev::compute_impl() {
  // Nothing to do here
}

}  // namespace scream
