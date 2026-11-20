#include "share/diagnostics/binary_op.hpp"

#include "share/physics/physics_constants.hpp"
#include "share/util/eamxx_quantity.hpp"

#include <cmath>
#include <cstdlib>

namespace {
  constexpr int OP_PLUS = 0;
  constexpr int OP_MINUS = 1;
  constexpr int OP_TIMES = 2;
  constexpr int OP_OVER = 3;
}

namespace scream {

namespace {
// Key, in the header of a rank-0 field made of constants alone, of its Quantity<Real>
constexpr auto constant_key = "constant_value";

// A non-field operand is a named physical constant (Rgas, P0, ...) or a numeric
// literal. Only expressions can write a literal -- (qc+qv)*2 -- but accept both.
// NOTE: strtod because from_chars for floating point is missing on some
//       compilers we support (e.g. compy).
bool str2real (const std::string& s, Real& value) {
  if (s.empty()) return false;

  char* endptr;
  const char* startptr = s.c_str();
  const double d = std::strtod(startptr,&endptr);

  if (endptr==startptr or *endptr!='\0') return false;
  // strtod also takes "nan"/"inf"; those are field names to us, not numbers.
  if (not std::isfinite(d)) return false;

  value = static_cast<Real>(d);
  return true;
}

// Resolves 's' as a scalar operand, if it is one. A numeric literal has no units
// of its own, so it gets invalid units (see apply_binary_op). The dictionary is
// checked first: a lookup is cheaper than a strtod, and named constants are the
// common case. Called once per operand, from the constructor.
std::optional<Quantity<Real>> scalar_operand (const std::string& s) {
  const auto& dict = physics::Constants<Real>::dictionary();
  auto it = dict.find(s);
  if (it!=dict.end()) {
    return it->second;
  }
  Real value;
  if (str2real(s,value)) {
    return Quantity<Real>(value);
  }
  return std::nullopt;
}
} // anonymous namespace

// parse string to get the operator code
int get_binary_operator_code(const std::string& op) {
  if (op == "plus") return OP_PLUS; // addition
  if (op == "minus") return OP_MINUS; // subtraction
  if (op == "times") return OP_TIMES; // multiplication
  if (op == "over") return OP_OVER; // division
  return -1; // invalid operator
}
// apply binary operation on two input units. A numeric literal has invalid units: in a sum it
// takes the other's, in a product it is a pure number. A sum of literals alone is a literal.
ekat::units::Units apply_binary_op(ekat::units::Units a, ekat::units::Units b, const int op_code) {
  const auto invalid = ekat::units::Units::invalid();

  if (a==invalid and b==invalid) {
    return invalid;
  }

  if (op_code==OP_PLUS or op_code==OP_MINUS) {
    if (a==invalid) return b;
    if (b==invalid) return a;
    EKAT_REQUIRE_MSG(a == b,
        "Error! Addition/subtraction requires compatible units.\n"
        " - lhs units: " + a.get_si_string() + "\n"
        " - rhs units: " + b.get_si_string() + "\n");
    return a;
  }

  EKAT_REQUIRE_MSG (op_code==OP_TIMES or op_code==OP_OVER,
      "Error! Unrecognized/unsupported operation code (" + std::to_string(op_code) + ").\n");
  if (a==invalid) a = ekat::units::none;
  if (b==invalid) b = ekat::units::none;
  return op_code==OP_TIMES ? a*b : a/b;
}
// apply binary operation on two input fields
void apply_binary_op(const Field& d, const Field& a, const Real b, const int op_code){
  bool masked = d.has_valid_mask();
  switch (op_code) {
    case OP_PLUS:
      if (masked) {
        d.update(a,1,0,b,d.get_valid_mask());
      } else {
        d.update(a,1,0,b);
      }
      break;
    case OP_MINUS:
      if (masked) {
        d.update(a,1,0,-b,d.get_valid_mask());
      } else {
        d.update(a,1,0,-b);
      }
      break;
    case OP_TIMES:
      if (masked) {
        d.update(a,b,0,0,d.get_valid_mask());
      } else {
        d.update(a,b,0,0);
      }
      break;
    case OP_OVER:
      if (masked) {
        d.update(a,1/b,0,0,d.get_valid_mask());
      } else {
        d.update(a,1/b,0,0);
      }
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized/unsupported binary op code (" + std::to_string(op_code) + ")\n");
  }
}

void apply_binary_op(const Field& d, const Real a, const Field& b, const int op_code)
{
  bool masked = d.has_valid_mask();
  switch (op_code) {
    case OP_PLUS:
      if (masked) {
        d.update(b,1,0,a,d.get_valid_mask());
      } else {
        d.update(b,1,0,a);
      }
      break;
    case OP_MINUS:
      if (masked) {
        d.update(b,-1,0,a,d.get_valid_mask());
      } else {
        d.update(b,-1,0,a);
      }
      break;
    case OP_TIMES:
      if (masked) {
        d.update(b,a,0,0,d.get_valid_mask());
      } else {
        d.update(b,a,0,0);
      }
      break;
    case OP_OVER:
      if (masked) {
        auto mask = d.get_valid_mask();
        d.deep_copy(a,mask);
        d.scale_inv(b,mask);
      } else {
        d.deep_copy(a);
        d.scale_inv(b);
      }
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized/unsupported binary op code (" + std::to_string(op_code) + ")\n");
  }
}

void apply_binary_op(const Field& d, const Field& a, const Field& b, const int op_code)
{
  bool masked = d.has_valid_mask();
  if (masked)
    d.deep_copy(a,d.get_valid_mask());
  else
    d.deep_copy(a);

  switch (op_code) {
    case OP_PLUS:
      if (masked)
        d.update(b, 1, 1, d.get_valid_mask());
      else
        d.update(b, 1, 1);
      break;
    case OP_MINUS:
      if (masked)
        d.update(b, -1, 1, d.get_valid_mask());
      else
        d.update(b, -1, 1);
      break;
    case OP_TIMES:
      if (masked)
        d.scale(b, d.get_valid_mask());
      else
        d.scale(b);
      break;
    case OP_OVER:
      if (masked)
        d.scale_inv(b, d.get_valid_mask());
      else
        d.scale_inv(b);
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized/unsupported binary op code (" + std::to_string(op_code) + ")\n");
  }
}
// apply binary operation on two scalars, with their units
Quantity<Real> scalar_op(const Quantity<Real>& a, const Quantity<Real>& b, const int op_code){
  const auto units = apply_binary_op(a.units,b.units,op_code);
  switch (op_code) {
    case OP_PLUS:  return Quantity<Real>(a.value+b.value,units);
    case OP_MINUS: return Quantity<Real>(a.value-b.value,units);
    case OP_TIMES: return Quantity<Real>(a.value*b.value,units);
    case OP_OVER:
      EKAT_REQUIRE_MSG (b.value!=0,
          "Error! Division by zero in an operation on constants.\n"
          " - lhs: " + std::to_string(a.value) + "\n");
      return Quantity<Real>(a.value/b.value,units);
    default:
      EKAT_ERROR_MSG("Error! Unrecognized/unsupported binary op code (" + std::to_string(op_code) + ")\n");
  }
  return Quantity<Real>(); // Unreachable, but some compilers complain
}

BinaryOp::
BinaryOp(const ekat::Comm &comm,
         const ekat::ParameterList &params,
         const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm, params, grid)
{
  m_arg1_name = m_params.get<std::string>("arg1");
  m_arg2_name = m_params.get<std::string>("arg2");
  m_binary_op_str = m_params.get<std::string>("binary_op");
  m_binary_op_code= get_binary_operator_code(m_binary_op_str);
  
  // Validate operator
  EKAT_REQUIRE_MSG(m_binary_op_code >= 0,
                   "Error! Invalid binary operator: '" + m_binary_op_str + "'\n"
                   "Valid operators are: plus, minus, times, over\n");

  // A non-field operand's value/units never change, so resolve them here rather
  // than on every compute_impl call.
  m_arg1_scalar = scalar_operand(m_arg1_name);
  m_arg2_scalar = scalar_operand(m_arg2_name);
  if (not m_arg1_scalar)
    m_field_in_names.push_back(m_arg1_name);
  if (not m_arg2_scalar)
    m_field_in_names.push_back(m_arg2_name);
}

void BinaryOp::initialize_impl()
{
  // A rank-0 field made of constants alone is a scalar (it is still an input, so computed first)
  auto as_constant = [&](std::optional<Quantity<Real>>& scalar, const std::string& name) {
    if (not scalar) {
      const auto& h = m_fields_in.at(name).get_header();
      if (h.has_extra_data(constant_key)) {
        scalar = h.get_extra_data<Quantity<Real>>(constant_key);
      }
    }
  };
  as_constant(m_arg1_scalar,m_arg1_name);
  as_constant(m_arg2_scalar,m_arg2_name);

  if (not m_arg1_scalar and not m_arg2_scalar) {
    const auto& fid1   = m_fields_in.at(m_arg1_name).get_header().get_identifier();
    const auto& fid2   = m_fields_in.at(m_arg2_name).get_header().get_identifier();

    // Must be on same layout, same datatype, same grid
    const bool rank0_vs_not = (fid1.get_layout().rank()==0) != (fid2.get_layout().rank()==0);
    EKAT_REQUIRE_MSG (fid1.get_layout() == fid2.get_layout(),
      "Error! BinaryOp requires both input fields to have the same layout.\n"
      " - field 1 name: " + fid1.name() + "\n"
      " - field 1 layout: " + fid1.get_layout().to_string() + "\n"
      " - field 2 name: " + fid2.name() + "\n"
      " - field 2 layout: " + fid2.get_layout().to_string() + "\n" +
      (rank0_vs_not
        ? " - note: a rank-0 field is a scalar operand only if it is made of constants alone\n"
          "   (e.g., Rair/Cpair). Any other must be broadcast explicitly, as in X.broadcast_like(Y).\n"
        : ""));
    EKAT_REQUIRE_MSG (fid1.data_type() == fid2.data_type(),
      "Error! BinaryOp requires both input fields to have the same data type.\n"
      " - field 1 name: " + fid1.name() + "\n"
      " - field 1 data type: " + e2str(fid1.data_type()) + "\n"
      " - field 2 name: " + fid2.name() + "\n"
      " - field 2 data type: " + e2str(fid2.data_type()) + "\n");
    EKAT_REQUIRE_MSG (fid1.get_grid_name() == fid2.get_grid_name(),
      "Error! BinaryOp requires both input fields to be on the same grid.\n"
      " - field 1 name: " + fid1.name() + "\n"
      " - field 1 grid name: " + fid1.get_grid_name() + "\n"
      " - field 2 name: " + fid2.name() + "\n"
      " - field 2 grid name: " + fid2.get_grid_name() + "\n");
  }

  // All good, create the diag output
  auto fid = [&](const std::string& name) -> const FieldIdentifier& {
    return m_fields_in.at(name).get_header().get_identifier();
  };
  auto dl = m_arg1_scalar ? (m_arg2_scalar ? FieldLayout({},{}) : fid(m_arg2_name).get_layout())
                          : fid(m_arg1_name).get_layout();
  auto u1 = m_arg1_scalar ? m_arg1_scalar->units : fid(m_arg1_name).get_units();
  auto u2 = m_arg2_scalar ? m_arg2_scalar->units : fid(m_arg2_name).get_units();
  // The units of an operation on literals alone are invalid, but a field needs units
  auto diag_units = apply_binary_op(u1, u2, m_binary_op_code);
  if (diag_units==ekat::units::Units::invalid()) {
    diag_units = ekat::units::none;
  }
  auto gn = m_params.get<std::string>("grid_name");
  auto diag_name = m_arg1_name + "_" + m_binary_op_str + "_" + m_arg2_name;
  FieldIdentifier d_fid(diag_name, dl, diag_units, gn);
  add_output(d_fid);

  if (m_arg1_scalar and m_arg2_scalar) {
    // We can pre-compute the diag now. It is a constant, which others use as a scalar.
    const auto q = scalar_op(*m_arg1_scalar, *m_arg2_scalar, m_binary_op_code);
    get().deep_copy(q.value);
    get().get_header().set_extra_data(constant_key,q);
  } else {
    m_arg1_has_mask = not m_arg1_scalar and m_fields_in.at(m_arg1_name).has_valid_mask();
    m_arg2_has_mask = not m_arg2_scalar and m_fields_in.at(m_arg2_name).has_valid_mask();
    if (m_arg1_has_mask or m_arg2_has_mask) {
      get().create_valid_mask();
      get().get_header().set_may_be_filled(true);
    }
  }
}

void BinaryOp::compute_impl()
{
  // First, if one (or both) the args are masked, compute the diag mask
  if (m_arg1_has_mask) {
    get().get_valid_mask().deep_copy(m_fields_in.at(m_arg1_name).get_valid_mask());
    if (m_arg2_has_mask)
      get().get_valid_mask().scale(m_fields_in.at(m_arg2_name).get_valid_mask());
  } else if (m_arg2_has_mask) {
    get().get_valid_mask().deep_copy(m_fields_in.at(m_arg2_name).get_valid_mask());
  }

  // Note the case where both args are scalars was handled in the initialize method
  if (not m_arg1_scalar and not m_arg2_scalar) {
    const auto& f1 = m_fields_in.at(m_arg1_name);
    const auto& f2 = m_fields_in.at(m_arg2_name);
    apply_binary_op(get(), f1, f2, m_binary_op_code);
  } else if (not m_arg1_scalar) {
    const auto& f1 = m_fields_in.at(m_arg1_name);
    apply_binary_op(get(), f1, m_arg2_scalar->value, m_binary_op_code);
  } else if (not m_arg2_scalar) {
    // We can do scale/scale_inv for * and /, but for + and - we must deep copy arg2 first
    const auto& f2 = m_fields_in.at(m_arg2_name);
    apply_binary_op(get(), m_arg1_scalar->value, f2, m_binary_op_code);
  }

  // Until IO is ready to fully rely on valid_mask fields, we must set diag=fill_value where invalid
  if (get().has_valid_mask()) {
    constexpr auto fv = constants::fill_value<Real>;
    get().deep_copy(fv,get().get_valid_mask(),true);
  }
}

}  // namespace scream
