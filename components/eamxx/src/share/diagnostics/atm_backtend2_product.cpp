// ============================================================
//  atm_backtend2_product.cpp
// ============================================================
#include "atm_backtend2_product.hpp"
#include <ekat/ekat_assert.hpp>

namespace scream {

AtmBacktend2Product::AtmBacktend2Product(const ekat::Comm& comm,
                                        const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm, params)
{
  EKAT_REQUIRE_MSG(params.isParameter("tendency_name"),
    "Error! AtmBacktend2Product requires 'tendency_name' in its input parameters.\n");

  m_name      = params.get<std::string>("tendency_name");
  m_diag_name = m_name + "_atm_backtend2_product";
}

void AtmBacktend2Product::create_requests()
{
  const auto &gname = m_params.get<std::string>("grid_name");
  add_field<Required>(m_name, gname);
}

void AtmBacktend2Product::initialize_impl(const RunType /* run_type */)
{
  const auto& f_in   = get_field_in(m_name);
  const auto& fid    = f_in.get_header().get_identifier();
  const auto& layout = fid.get_layout();
  const auto& gn     = fid.get_grid_name();

  EKAT_REQUIRE_MSG(
      f_in.data_type() == DataType::RealType,
      "Error! AtmBacktend2Product only supports Real data type field.\n"
      " - field name: " + fid.name() + "\n"
      " - field data type: " + e2str(f_in.data_type()) + "\n");

  // Output units are [input_units]^2.
  FieldIdentifier d_fid(m_diag_name, layout.clone(), fid.get_units() * fid.get_units(), gn);
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  m_f_pre  = f_in.clone(m_name + "_sbp_f_pre");
  m_bt_pre = f_in.clone(m_name + "_sbp_bt_pre");

  m_f_pre.deep_copy(f_in);
  m_bt_pre.deep_copy(0.0);
}

void AtmBacktend2Product::compute_diagnostic_impl()
{
  const auto& f_in = get_field_in(m_name);

  // diag = f_cur - f_pre = bt_cur
  m_diagnostic_output.deep_copy(f_in);
  m_diagnostic_output.update(m_f_pre, -1.0, 1.0);

  // diag = bt_cur * old_bt_pre
  m_diagnostic_output.scale(m_bt_pre);

  // Advance carry-forward state: bt_pre = bt_cur = f_cur - old_f_pre
  m_bt_pre.deep_copy(f_in);
  m_bt_pre.update(m_f_pre, -1.0, 1.0);

  // Advance carry-forward state: f_pre = f_cur
  m_f_pre.deep_copy(f_in);
}

} // namespace scream
