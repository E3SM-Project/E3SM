#include "share/diagnostics/field_broadcast.hpp"

namespace scream {

FieldBroadcast::
FieldBroadcast (const ekat::Comm &comm,
                const ekat::ParameterList &params,
                const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm, params, grid)
{
  EKAT_REQUIRE_MSG (params.isParameter("field_name") and params.isParameter("target_name"),
      "Error! FieldBroadcast requires 'field_name' and 'target_name' in its input parameters.\n");

  m_name     = m_params.get<std::string>("field_name");
  m_tgt_name = m_params.get<std::string>("target_name");
  m_field_in_names.push_back(m_name);
  if (m_tgt_name!=m_name) {
    m_field_in_names.push_back(m_tgt_name);
  }
}

void FieldBroadcast::initialize_impl ()
{
  const auto& f   = m_fields_in.at(m_name);
  const auto& tgt = m_fields_in.at(m_tgt_name);
  const auto& fid = f.get_header().get_identifier();
  const auto& tl  = tgt.get_header().get_identifier().get_layout();

  EKAT_REQUIRE_MSG (fid.get_grid_name()==tgt.get_header().get_identifier().get_grid_name(),
      "Error! FieldBroadcast requires the field and the target to be on the same grid.\n"
      " - field name: " + m_name + "\n"
      " - target name: " + m_tgt_name + "\n");

  // NOTE: this throws if the layout of f is not an ordered subset of the target one
  m_broadcasted = f.broadcast_to(tl);

  // NOTE: most diags ASSUME `get_view` works, but for a broadcasted field it does NOT work.
  //       Hence, we must deep_copy at run time. If we ever make other diags check whether
  //       get_view or get_strided_view can be used (or if we switch to ALWAYS strided view)
  //       we can change m_diagnostic_output to be the same as m_broadcasted
  FieldIdentifier d_fid(m_name + "_broadcast_to_" + m_tgt_name, tl, fid.get_units(),
                        fid.get_grid_name(), fid.data_type());
  m_diagnostic_output = Field(d_fid,true);
  if (f.has_valid_mask()) {
    m_diagnostic_output.create_valid_mask();
    m_diagnostic_output.get_header().set_may_be_filled(true);
  }
}

void FieldBroadcast::compute_impl ()
{
  // See comment above: remove when m_diagnostic_output aliases m_broadcasted
  if (m_broadcasted.has_valid_mask()) {
    const auto& mask = m_broadcasted.get_valid_mask();
    m_diagnostic_output.get_valid_mask().deep_copy(mask);
    m_diagnostic_output.deep_copy(m_broadcasted,mask);
  } else {
    m_diagnostic_output.deep_copy(m_broadcasted);
  }
}

}  // namespace scream
