#include "share/diagnostics/unary_math.hpp"

#include "share/field/field_pfor.hpp"

#include <ekat_units.hpp>

namespace scream {

namespace {

using Op = UnaryMath::Op;

// Whether x is in the domain of the op (only log and sqrt have one)
KOKKOS_INLINE_FUNCTION
bool in_domain (const Op op, const Real x) {
  switch (op) {
    case Op::Log:  return x>0;
    case Op::Sqrt: return x>=0;
    default:       return true;
  }
}

KOKKOS_INLINE_FUNCTION
Real apply (const Op op, const Real x) {
  switch (op) {
    case Op::Log:  return Kokkos::log(x);
    case Op::Exp:  return Kokkos::exp(x);
    case Op::Sqrt: return Kokkos::sqrt(x);
    case Op::Abs:  return Kokkos::abs(x);
  }
  return 0; // Unreachable, but some compilers complain
}

// NOTE: kernels that go through field_pfor must be functors with a templated
//       operator() (not generic lambdas), since extended lambdas cannot be generic.
// NOTE: the op is a runtime member rather than a template parameter. field_pfor
//       instantiates a kernel for every rank and layout combination, so one
//       kernel per op would multiply the compile time. All threads take the same
//       branch, so the switch costs little.

// Output not masked: compute everywhere, so an invalid x gives NaN (or -inf)
struct UnaryMathKernel {
  Op op;

  template<typename YView, typename XView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const YView& y, const XView& x, const Idx&... i) const {
    y(i...) = apply(op,x(i...));
  }
};

// Output mask (for invalid=mask) from the domain alone, for an unmasked x.
// NOTE: only compares x with 0, so it is safe where x holds a fill value
struct DomainMaskKernel {
  Op op;

  template<typename MaskView, typename XView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const MaskView& mask, const XView& x, const Idx&... i) const {
    mask(i...) = in_domain(op,x(i...)) ? 1 : 0;
  }
};

// Same, for a masked x: the output mask is the input one, minus the entries out of domain
struct DomainMaskAndKernel {
  Op op;

  template<typename MaskView, typename InMaskView, typename XView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const MaskView& mask, const InMaskView& in_mask, const XView& x, const Idx&... i) const {
    mask(i...) = (in_mask(i...)!=0 and in_domain(op,x(i...))) ? 1 : 0;
  }
};

// Output masked: compute only where the mask is set, and set the fill value
// elsewhere (until IO relies on valid_mask fields alone)
struct MaskedKernel {
  Op   op;
  Real fill_value;

  template<typename YView, typename MaskView, typename XView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const YView& y, const MaskView& mask, const XView& x, const Idx&... i) const {
    y(i...) = mask(i...)!=0 ? apply(op,x(i...)) : fill_value;
  }
};

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
  else {
    EKAT_ERROR_MSG (
        "Error! Invalid op for UnaryMath.\n"
        " - input: '" + m_op_str + "'\n"
        " - valid ops: log, exp, sqrt, abs\n");
  }

  const auto invalid = m_params.get<std::string>("invalid","nan");
  EKAT_REQUIRE_MSG (invalid=="nan" or invalid=="mask",
      "Error! Invalid value for UnaryMath's 'invalid'.\n"
      " - input: '" + invalid + "'\n"
      " - valid values: 'nan', 'mask'\n");
  EKAT_REQUIRE_MSG (invalid=="nan" or m_op==Op::Log or m_op==Op::Sqrt,
      "Error! 'invalid' only applies to log and sqrt, which have invalid inputs.\n"
      " - op: " + m_op_str + "\n");
  m_mask_invalid = invalid=="mask";

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

  // NOTE: do not use get(name,default): it would store the default in the params, which would
  //       change the diag hash (computed to share the diag among streams) after the fact
  const auto units = m_params.isParameter("units") ? m_params.get<std::string>("units")
                                                   : std::string("check");
  EKAT_REQUIRE_MSG (units=="check" or units=="ignore",
      "Error! Invalid value for UnaryMath's 'units'.\n"
      " - input: '" + units + "'\n"
      " - valid values: 'check', 'ignore'\n");
  EKAT_REQUIRE_MSG (units=="check" or m_op==Op::Log or m_op==Op::Exp,
      "Error! 'units' only applies to log and exp; sqrt and abs always have valid units.\n"
      " - op: " + m_op_str + "\n");

  const auto& u = fid.get_units();
  Units diag_units = none;
  if (m_op==Op::Log or m_op==Op::Exp) {
    // NOTE: is_dimensionless() ignores the scale: g/kg has no dimension, but its
    //       values are 1000 times the quantity, so log/exp of them is not log/exp of it.
    EKAT_REQUIRE_MSG (units=="ignore" or u==none,
        "Error! UnaryMath: " + m_op_str + " requires a dimensionless, unscaled input.\n"
        " - field name: " + fid.name() + "\n"
        " - field units: " + u.to_string() + "\n"
        " - kg/kg is fine, but g/kg or K are not.\n"
        " - to apply it to the field's values regardless, use units='ignore'.\n");
  } else if (m_op==Op::Sqrt) {
    diag_units = sqrt(u);
  } else {
    diag_units = u;
  }

  auto& diag = add_output(FieldIdentifier(m_op_str + "_of_" + m_field_name,
                                          fid.get_layout().clone(), diag_units,
                                          fid.get_grid_name()));

  m_out_masked = f.has_valid_mask() or m_mask_invalid;
  if (m_mask_invalid) {
    diag.create_valid_mask();
  } else if (f.has_valid_mask()) {
    // The output is invalid exactly where the input is
    diag.set_valid_mask(f.get_valid_mask());
  }
}

void UnaryMath::compute_impl()
{
  using details::field_pfor;
  using details::view_as;

  const auto& f = m_fields_in.at(m_field_name);
  const auto  d = get();

  if (not m_out_masked) {
    field_pfor(UnaryMathKernel{m_op},
               view_as<Real>(d), view_as<const Real>(f));
  } else {
    const auto& mask = d.get_valid_mask();
    if (m_mask_invalid) {
      if (f.has_valid_mask()) {
        field_pfor(DomainMaskAndKernel{m_op},
                   view_as<int>(mask), view_as<const int>(f.get_valid_mask()),
                   view_as<const Real>(f));
      } else {
        field_pfor(DomainMaskKernel{m_op},
                   view_as<int>(mask), view_as<const Real>(f));
      }
    }
    // Kernels on the same execution space run in order: no fence needed between them
    field_pfor(MaskedKernel{m_op,constants::fill_value<Real>},
               view_as<Real>(d), view_as<const int>(mask), view_as<const Real>(f));
  }

  // field_pfor is asynchronous, but compute_impl must be blocking
  Kokkos::fence();
}

} // namespace scream
