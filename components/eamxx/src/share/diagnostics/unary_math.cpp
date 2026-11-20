#include "share/diagnostics/unary_math.hpp"

#include "share/field/field_pfor.hpp"

#include <ekat_units.hpp>

#include <cmath>
#include <cstdlib>

namespace scream {

namespace {

using Op = UnaryMath::Op;

// Whether x is in the domain of the op (only log, sqrt and pow can be outside it)
KOKKOS_INLINE_FUNCTION
bool in_domain (const Op op, const Real x, const Real p, const bool p_is_int) {
  switch (op) {
    case Op::Log:  return x>0;
    case Op::Sqrt: return x>=0;
    case Op::Pow:  return (x>0) or (x==0 and p>=0) or (x<0 and p_is_int);
    default:       return true;
  }
}

KOKKOS_INLINE_FUNCTION
Real apply (const Op op, const Real x, const Real p) {
  switch (op) {
    case Op::Log:  return Kokkos::log(x);
    case Op::Exp:  return Kokkos::exp(x);
    case Op::Sqrt: return Kokkos::sqrt(x);
    case Op::Abs:  return Kokkos::abs(x);
    case Op::Pow:  return Kokkos::pow(x,p);
  }
  return 0; // Unreachable, but some compilers complain
}

// NOTE: kernels that go through field_pfor must be functors with a templated
//       operator() (not generic lambdas), since extended lambdas cannot be generic.

// No output mask: compute everywhere, so an invalid x gives NaN (or -inf)
struct UnaryMathKernel {
  Op   op;
  Real p;

  template<typename YView, typename XView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const YView& y, const XView& x, const Idx&... i) const {
    y(i...) = apply(op,x(i...),p);
  }
};

// Output mask, already holding X's mask (or 1 if X has none). Entries already
// masked are skipped. If mask_invalid=true, entries outside the op's domain are
// masked; otherwise they are computed anyways (giving NaN or -inf).
struct UnaryMathMaskedKernel {
  Op   op;
  Real p;
  bool p_is_int;
  bool mask_invalid;

  template<typename YView, typename MaskView, typename XView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const YView& y, const MaskView& mask, const XView& x, const Idx&... i) const {
    if (mask(i...)==0) {
      return;
    }
    const auto xi = x(i...);
    if (mask_invalid and not in_domain(op,xi,p,p_is_int)) {
      mask(i...) = 0;
    } else {
      y(i...) = apply(op,xi,p);
    }
  }
};

// The exponent as n/d, with 1<=d<=12, if it is one
bool as_rational (const Real p, ekat::RationalConstant::iType& n, ekat::RationalConstant::iType& d) {
  for (d=1; d<=12; ++d) {
    const auto pd = static_cast<double>(p)*d;
    const auto r  = std::round(pd);
    if (std::abs(pd-r) <= 1e-12*std::max(1.0,std::abs(pd))) {
      n = static_cast<ekat::RationalConstant::iType>(r);
      return true;
    }
  }
  return false;
}

} // anonymous namespace

UnaryMath::
UnaryMath (const ekat::Comm &comm,
           const ekat::ParameterList &params,
           const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm, params, grid)
{
  EKAT_REQUIRE_MSG (m_params.isParameter("field_name"),
      "Error! UnaryMath requires 'field_name' in its input parameters.\n");
  EKAT_REQUIRE_MSG (m_params.isParameter("op"),
      "Error! UnaryMath requires 'op' in its input parameters.\n");

  m_field_name = m_params.get<std::string>("field_name");
  m_op_str     = m_params.get<std::string>("op");

  if      (m_op_str=="log")  m_op = Op::Log;
  else if (m_op_str=="exp")  m_op = Op::Exp;
  else if (m_op_str=="sqrt") m_op = Op::Sqrt;
  else if (m_op_str=="abs")  m_op = Op::Abs;
  else if (m_op_str=="pow")  m_op = Op::Pow;
  else {
    EKAT_ERROR_MSG (
        "Error! Invalid op for UnaryMath.\n"
        " - input: '" + m_op_str + "'\n"
        " - valid ops: log, exp, sqrt, abs, pow\n");
  }

  const bool is_pow = m_op==Op::Pow;
  EKAT_REQUIRE_MSG (is_pow == m_params.isParameter("exponent"),
      "Error! UnaryMath takes 'exponent' if and only if op is 'pow'.\n"
      " - op: " + m_op_str + "\n");
  if (is_pow) {
    // NOTE: strtod because from_chars for floating point is missing on some
    //       compilers we support (e.g. compy).
    const auto& s = m_params.get<std::string>("exponent");
    char* end;
    const double p = std::strtod(s.c_str(),&end);
    EKAT_REQUIRE_MSG (not s.empty() and *end=='\0' and std::isfinite(p),
        "Error! UnaryMath's 'exponent' must be a finite number.\n"
        " - input: '" + s + "'\n");
    m_exponent = p;
  }

  const auto invalid = m_params.get<std::string>("invalid","nan");
  EKAT_REQUIRE_MSG (invalid=="nan" or invalid=="mask",
      "Error! Invalid value for UnaryMath's 'invalid'.\n"
      " - input: '" + invalid + "'\n"
      " - valid values: 'nan', 'mask'\n");
  EKAT_REQUIRE_MSG (invalid=="nan" or m_op==Op::Log or m_op==Op::Sqrt or is_pow,
      "Error! 'invalid' only applies to log, sqrt and pow, which have invalid inputs.\n"
      " - op: " + m_op_str + "\n");
  m_mask_invalid = invalid=="mask";

  const auto units = m_params.get<std::string>("units","check");
  EKAT_REQUIRE_MSG (units=="check" or units=="ignore",
      "Error! Invalid value for UnaryMath's 'units'.\n"
      " - input: '" + units + "'\n"
      " - valid values: 'check', 'ignore'\n");
  EKAT_REQUIRE_MSG (units=="check" or m_op==Op::Log or m_op==Op::Exp or is_pow,
      "Error! 'units' only applies to log, exp and pow; sqrt and abs always have valid units.\n"
      " - op: " + m_op_str + "\n");
  m_ignore_units = units=="ignore";

  m_field_in_names.push_back(m_field_name);
}

void UnaryMath::initialize_impl()
{
  using namespace ekat::units;

  const auto& f   = m_fields_in.at(m_field_name);
  const auto& fid = f.get_header().get_identifier();

  EKAT_REQUIRE_MSG (f.data_type()==DataType::RealType,
      "Error! UnaryMath only supports Real data type fields.\n"
      " - field name: " + fid.name() + "\n"
      " - field data type: " + e2str(f.data_type()) + "\n");

  const auto& u = fid.get_units();
  Units diag_units = none;
  if (not m_ignore_units) {
    switch (m_op) {
      case Op::Log: [[fallthrough]];
      case Op::Exp:
        EKAT_REQUIRE_MSG (u.is_dimensionless(),
            "Error! UnaryMath: " + m_op_str + " requires a dimensionless input.\n"
            " - field name: " + fid.name() + "\n"
            " - field units: " + u.to_string() + "\n"
            " - to apply it to the field's values regardless, use units='ignore'.\n");
        break;
      case Op::Sqrt:
        diag_units = sqrt(u);
        break;
      case Op::Abs:
        diag_units = u;
        break;
      case Op::Pow:
      {
        // Dimensionless and unscaled (e.g., kg/kg) stays so for any exponent
        if (not u.is_dimensionless() or u.is_scaled()) {
          ekat::RationalConstant::iType n, d;
          EKAT_REQUIRE_MSG (as_rational(m_exponent,n,d),
              "Error! UnaryMath: the units of pow(X,p) are u^p, which needs p=n/d with d<=12.\n"
              " - field name: " + fid.name() + "\n"
              " - field units: " + u.to_string() + "\n"
              " - exponent: " + m_params.get<std::string>("exponent") + "\n"
              " - to apply it to the field's values regardless, use units='ignore'.\n");
          diag_units = pow(u,ekat::RationalConstant(n,d));
        }
        break;
      }
    }
  }

  std::string diag_name = m_op_str + "_of_" + m_field_name;
  if (m_op==Op::Pow) {
    diag_name += "_to_" + m_params.get<std::string>("exponent");
  }
  FieldIdentifier d_fid(diag_name, fid.get_layout().clone(), diag_units, fid.get_grid_name());
  add_output(d_fid);

  m_out_masked = f.has_valid_mask() or m_mask_invalid;
  if (m_out_masked) {
    get().create_valid_mask();
    get().get_header().set_may_be_filled(true);
  }
}

void UnaryMath::compute_impl()
{
  using details::field_pfor;
  using details::view_as;

  const auto& f = m_fields_in.at(m_field_name);
  const auto  d = get();

  if (not m_out_masked) {
    field_pfor(UnaryMathKernel{m_op,m_exponent},
               view_as<Real>(d), view_as<const Real>(f));
  } else {
    const auto& mask = d.get_valid_mask();
    if (f.has_valid_mask()) {
      mask.deep_copy(f.get_valid_mask());
    } else {
      mask.deep_copy(1);
    }
    const bool p_is_int = m_exponent==std::round(m_exponent);
    field_pfor(UnaryMathMaskedKernel{m_op,m_exponent,p_is_int,m_mask_invalid},
               view_as<Real>(d), view_as<int>(mask), view_as<const Real>(f));

    // Until IO is ready to fully rely on valid_mask fields, we must set diag=fill_value where invalid
    d.deep_copy(constants::fill_value<Real>,mask,true);
  }
  Kokkos::fence();
}

} // namespace scream
