#include "binary_ops.hpp"

#include "share/physics/physics_constants.hpp"

namespace {
  constexpr int OP_PLUS = 0;
  constexpr int OP_MINUS = 1;
  constexpr int OP_TIMES = 2;
  constexpr int OP_OVER = 3;
}

namespace scream {

// parse string to get the operator code
int get_binary_operator_code(const std::string& op) {
  if (op == "plus") return OP_PLUS; // addition
  if (op == "minus") return OP_MINUS; // subtraction
  if (op == "times") return OP_TIMES; // multiplication
  if (op == "over") return OP_OVER; // division
  return -1; // invalid operator
}
// apply binary operation on two input units
ekat::units::Units apply_binary_op(const ekat::units::Units& a, const ekat::units::Units& b, const int op_code) {
  switch (op_code) {
    case OP_PLUS: [[fallthrough]];
    case OP_MINUS:
      EKAT_REQUIRE_MSG(a == b, "Error! Addition/subtraction requires compatible units.\n");
      return a;
    case OP_TIMES:
      return a * b;
    case OP_OVER:
      return a / b;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized/unsupported operation code (" + std::to_string(op_code) + ").\n");
  }
  return ekat::units::Units::invalid();
}
// apply binary operation on two input fields
void apply_binary_op(Field& d, const Field& a, const Real b, const int op_code){
  switch (op_code) {
    case OP_PLUS:
      d.deep_copy(b);
      d.update(a,1,1);
      break;
    case OP_MINUS:
      d.deep_copy(b);
      d.update(a,1,-1);
      break;
    case OP_TIMES:
      d.deep_copy(a);
      d.scale(b);
      break;
    case OP_OVER:
      d.deep_copy(a);
      d.scale(1 / b);
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized/unsupported binary op code (" + std::to_string(op_code) + ")\n");
  }
}
void apply_binary_op(Field& d, const Real a, const Field& b, const int op_code){
  switch (op_code) {
    case OP_PLUS:
      d.deep_copy(a);
      d.update(b,1,1);
      break;
    case OP_MINUS:
      d.deep_copy(a);
      d.update(b,-1,1);
      break;
    case OP_TIMES:
      d.deep_copy(a);
      d.scale(b);
      break;
    case OP_OVER:
      d.deep_copy(a);
      d.scale_inv(b);
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized/unsupported binary op code (" + std::to_string(op_code) + ")\n");
  }
}
void apply_binary_op(Field& d, const Field& a, const Field& b, const int op_code){
  d.deep_copy(a);
  switch (op_code) {
    case OP_PLUS:
      d.update(b, 1, 1); // addition
      break;
    case OP_MINUS:
      d.update(b, -1, 1); // subtraction
      break;
    case OP_TIMES:
      d.scale(b); // multiplication
      break;
    case OP_OVER:
      d.scale_inv(b); // division
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized/unsupported binary op code (" + std::to_string(op_code) + ")\n");
  }
}
void apply_binary_op(Field& d, const Real& a, const Real& b, const int op_code){
  switch (op_code) {
    case OP_PLUS:
      d.deep_copy(a+b);
      break;
    case OP_MINUS:
      d.deep_copy(a-b);
      break;
    case OP_TIMES:
      d.deep_copy(a*b);
      break;
    case OP_OVER:
      d.deep_copy(a/b);
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized/unsupported binary op code (" + std::to_string(op_code) + ")\n");
  }
}

BinaryOpsDiag::
BinaryOpsDiag(const ekat::Comm &comm,
                const ekat::ParameterList &params)
 : AtmosphereDiagnostic(comm, params)
{
  m_arg1_name = m_params.get<std::string>("arg1");
  m_arg2_name = m_params.get<std::string>("arg2");
  m_binary_op_str = m_params.get<std::string>("binary_op");
  m_binary_op_code= get_binary_operator_code(m_binary_op_str);
  
  // Validate operator
  EKAT_REQUIRE_MSG(m_binary_op_code >= 0,
                   "Error! Invalid binary operator: '" + m_binary_op_str + "'\n"
                   "Valid operators are: plus, minus, times, over\n");
}

void BinaryOpsDiag::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  const auto &gname = m_params.get<std::string>("grid_name");
  const auto& pc_dict = physics::Constants<Real>::dictionary();
  m_arg1_is_field = pc_dict.count(m_arg1_name)==0;
  m_arg2_is_field = pc_dict.count(m_arg2_name)==0;
  if (m_arg1_is_field)
    add_field<Required>(m_arg1_name, gname);
  if (m_arg2_is_field)
    add_field<Required>(m_arg2_name, gname);
}

void BinaryOpsDiag::initialize_impl(const RunType /*run_type*/)
{
  // get the inputs specs
  
  const auto& dict = physics::Constants<Real>::dictionary();

  if (m_arg1_is_field and m_arg2_is_field) {
    const auto& fid1   = get_field_in(m_arg1_name).get_header().get_identifier();
    const auto& fid2   = get_field_in(m_arg2_name).get_header().get_identifier();

    // Must be on same layout, same datatype, same grid
    EKAT_REQUIRE_MSG (fid1.get_layout() == fid2.get_layout(),
      "Error! BinaryOpsDiag requires both input fields to have the same layout.\n"
      " - field 1 name: " + fid1.name() + "\n"
      " - field 1 layout: " + fid1.get_layout().to_string() + "\n"
      " - field 2 name: " + fid2.name() + "\n"
      " - field 2 layout: " + fid2.get_layout().to_string() + "\n");
    EKAT_REQUIRE_MSG (fid1.data_type() == fid2.data_type(),
      "Error! BinaryOpsDiag requires both input fields to have the same data type.\n"
      " - field 1 name: " + fid1.name() + "\n"
      " - field 1 data type: " + e2str(fid1.data_type()) + "\n"
      " - field 2 name: " + fid2.name() + "\n"
      " - field 2 data type: " + e2str(fid2.data_type()) + "\n");
    EKAT_REQUIRE_MSG (fid1.get_grid_name() == fid1.get_grid_name(),
      "Error! BinaryOpsDiag requires both input fields to be on the same grid.\n"
      " - field 1 name: " + fid1.name() + "\n"
      " - field 1 grid name: " + fid1.get_grid_name() + "\n"
      " - field 2 name: " + fid2.name() + "\n"
      " - field 2 grid name: " + fid2.get_grid_name() + "\n");
  }

  // All good, create the diag output
  auto dl = m_arg1_is_field ? get_field_in(m_arg1_name).get_header().get_identifier().get_layout()
                            : (m_arg2_is_field ? get_field_in(m_arg2_name).get_header().get_identifier().get_layout()
                                               : FieldLayout({},{}));
  const auto& u1 = m_arg1_is_field ? get_field_in(m_arg1_name).get_header().get_identifier().get_units()
                                   : dict.at(m_arg1_name).units;
  const auto& u2 = m_arg2_is_field ? get_field_in(m_arg2_name).get_header().get_identifier().get_units()
                                   : dict.at(m_arg2_name).units;
  auto diag_units = apply_binary_op(u1, u2, m_binary_op_code);
  auto gn = m_params.get<std::string>("grid_name");
  auto diag_name = m_arg1_name + "_" + m_binary_op_str + "_" + m_arg2_name;
  FieldIdentifier d_fid(diag_name, dl, diag_units, gn);
  m_diagnostic_output = Field(d_fid,true);

  if (not m_arg1_is_field and not m_arg2_is_field) {
    // We can pre-compute the diag now
    compute_diagnostic_impl();
    m_diagnostic_output.get_header().get_tracking().update_time_stamp(start_of_step_ts());
  }
}

void BinaryOpsDiag::compute_diagnostic_impl()
{
  const auto& dict = physics::Constants<Real>::dictionary();
  if (m_arg1_is_field and m_arg2_is_field) {
    const auto& f1 = get_field_in(m_arg1_name);
    const auto& f2 = get_field_in(m_arg2_name);
    apply_binary_op(m_diagnostic_output, f1, f2, m_binary_op_code);
  } else if (m_arg1_is_field) {
    const auto& f1 = get_field_in(m_arg1_name);
    const auto  c2 = dict.at(m_arg2_name).value;
    apply_binary_op(m_diagnostic_output, f1, c2, m_binary_op_code);
  } else if (m_arg2_is_field) {
    // We can do scale/scale_inv for * and /, but for + and - we must deep copy arg2 first
    const auto  c1 = physics::Constants<Real>::dictionary().at(m_arg1_name).value;
    const auto& f2 = get_field_in(m_arg2_name);
    apply_binary_op(m_diagnostic_output, c1, f2, m_binary_op_code);
  } else {
    const auto  c1 = dict.at(m_arg1_name).value;
    const auto  c2 = dict.at(m_arg2_name).value;
    apply_binary_op(m_diagnostic_output, c1, c2, m_binary_op_code);
  }
}

}  // namespace scream
