#include "share/diagnostics/field_pow.hpp"

#include "share/field/field_pfor.hpp"

#include <ekat_units.hpp>

#include <algorithm>
#include <charconv>
#include <cmath>

namespace scream {

namespace {

using Kind = FieldPow::Kind;

// What the kernels need to know about the exponent
struct PowSpec {
  Kind kind;
  Real p;
  bool p_is_int;
};

// Whether x^p is defined. 0^0 is not, but p=0 is rejected on construction.
KOKKOS_INLINE_FUNCTION
bool in_domain (const PowSpec& s, const Real x) {
  return x>0 or (x==0 and s.p>0) or (x<0 and s.p_is_int);
}

KOKKOS_INLINE_FUNCTION
Real apply (const PowSpec& s, const Real x) {
  switch (s.kind) {
    case Kind::Square:  return x*x;
    case Kind::Cube:    return x*x*x;
    case Kind::Inverse: return 1/x;
    default:            return Kokkos::pow(x,s.p);
  }
}

// NOTE: kernels that go through field_pfor must be functors with a templated
//       operator() (not generic lambdas), since extended lambdas cannot be generic.
// NOTE: the kind of power is a runtime member rather than a template parameter,
//       since field_pfor instantiates a kernel for every rank and layout combination.

// Output not masked: compute everywhere, so an invalid x gives NaN (or inf)
struct PowKernel {
  PowSpec spec;

  template<typename YView, typename XView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const YView& y, const XView& x, const Idx&... i) const {
    y(i...) = apply(spec,x(i...));
  }
};

// Output mask (for invalid=mask) from the domain alone, for an unmasked x.
// NOTE: only compares x with 0, so it is safe where x holds a fill value
struct DomainMaskKernel {
  PowSpec spec;

  template<typename MaskView, typename XView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const MaskView& mask, const XView& x, const Idx&... i) const {
    mask(i...) = in_domain(spec,x(i...)) ? 1 : 0;
  }
};

// Same, for a masked x: the output mask is the input one, minus the entries out of domain
struct DomainMaskAndKernel {
  PowSpec spec;

  template<typename MaskView, typename InMaskView, typename XView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const MaskView& mask, const InMaskView& in_mask, const XView& x, const Idx&... i) const {
    mask(i...) = (in_mask(i...)!=0 and in_domain(spec,x(i...))) ? 1 : 0;
  }
};

// Output masked: compute only where the mask is set, and set the fill value
// elsewhere (until IO relies on valid_mask fields alone)
struct MaskedPowKernel {
  PowSpec spec;
  Real    fill_value;

  template<typename YView, typename MaskView, typename XView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const YView& y, const MaskView& mask, const XView& x, const Idx&... i) const {
    y(i...) = mask(i...)!=0 ? apply(spec,x(i...)) : fill_value;
  }
};

// Parses all of s as a number of type T (int or floating point).
// NOTE: from_chars rather than strto*/sto*: it does not depend on the locale, it
//       does not skip whitespace or take a '+', and it reports an out-of-range
//       value, where strtoll saturates (which would silently be another exponent).
//       dexpr already needs its floating point overload, to read literals.
template<typename T>
bool str2num (const std::string& s, T& value)
{
  const char* first = s.data();
  const char* last  = first + s.size();
  const auto [ptr,ec] = std::from_chars(first,last,value);
  return ec==std::errc{} and ptr==last;
}

// A unit is raised to a rational power, which a decimal exponent is not.
// So look for n/d that matches it, trying denominators 1,2,..,max_den in turn.
// (The match is to 1e-12, relative, so 0.5 and -1.5 do, but 0.2857 does not.)
// An exponent written as 'n/d' does not go through this.
constexpr int max_den = 12;
bool decimal_to_fraction (const double p, long long& n, long long& d)
{
  // A numerator past this cannot be a long long (and the cast would be undefined)
  constexpr double max_num = 9.0e18;
  for (d=1; d<=max_den; ++d) {
    const double pd = p*d;
    const double r  = std::round(pd);
    if (std::abs(r)<max_num and std::abs(pd-r) <= 1e-12*std::max(1.0,std::abs(pd))) {
      n = static_cast<long long>(r);
      return true;
    }
  }
  return false;
}

} // anonymous namespace

FieldPow::
FieldPow (const ekat::Comm &comm,
          const ekat::ParameterList &params,
          const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm, params, grid)
{
  EKAT_REQUIRE_MSG (m_params.isParameter("field_name"),
      "Error! FieldPow requires 'field_name' in its input parameters.\n");
  EKAT_REQUIRE_MSG (m_params.isParameter("exponent"),
      "Error! FieldPow requires 'exponent' in its input parameters.\n");

  m_field_name   = m_params.get<std::string>("field_name");
  m_exponent_str = m_params.get<std::string>("exponent");

  // An exponent is a fraction 'n/d', or a decimal number.
  // NOTE: everything about the exponent (zero, integer, which fast path, the fraction
  //       for the units) is decided here in double, from the string, so that narrowing
  //       to Real (which can round 2+2^-24 to 2, in single precision) changes none of it.
  const auto& s = m_exponent_str;
  const auto slash = s.find('/');
  double p;
  bool is_int;
  long long int_value = 0;
  if (slash!=std::string::npos) {
    EKAT_REQUIRE_MSG (str2num(s.substr(0,slash),m_num) and str2num(s.substr(slash+1),m_den) and m_den>0,
        "Error! FieldPow's 'exponent' must be a fraction of two integers (with a positive denominator) or a number.\n"
        " - input: '" + s + "'\n");
    p = static_cast<double>(m_num)/static_cast<double>(m_den);
    m_has_fraction = true;
    is_int = m_num % m_den == 0;
    int_value = m_num / m_den;
  } else {
    EKAT_REQUIRE_MSG (str2num(s,p) and std::isfinite(p),
        "Error! FieldPow's 'exponent' must be a fraction of two integers or a finite number.\n"
        " - input: '" + s + "'\n");
    m_has_fraction = decimal_to_fraction(p,m_num,m_den);
    is_int = p==std::round(p);
    if (is_int and std::abs(p)<9.0e18) {
      int_value = static_cast<long long>(p);
    }
  }
  m_exponent        = static_cast<Real>(p);
  m_exponent_is_int = is_int;

  EKAT_REQUIRE_MSG (p!=0 and m_exponent!=0,
      "Error! FieldPow's 'exponent' cannot be 0 (or so small that it is, as a Real): "
      "x^0 is 1, and 0^0 is not defined.\n"
      " - input: '" + s + "'\n");

  // The fast paths are for these exact integers only
  if      (is_int and int_value==2)  m_kind = Kind::Square;
  else if (is_int and int_value==3)  m_kind = Kind::Cube;
  else if (is_int and int_value==-1) m_kind = Kind::Inverse;
  else                               m_kind = Kind::General;

  const auto invalid = m_params.get<std::string>("invalid","nan");
  EKAT_REQUIRE_MSG (invalid=="nan" or invalid=="mask",
      "Error! Invalid value for FieldPow's 'invalid'.\n"
      " - input: '" + invalid + "'\n"
      " - valid values: 'nan', 'mask'\n");
  m_mask_invalid = invalid=="mask";

  m_field_in_names.push_back(m_field_name);
}

void FieldPow::initialize_impl()
{
  using namespace ekat::units;

  const auto& f   = m_fields_in.at(m_field_name);
  const auto& fid = f.get_header().get_identifier();

  EKAT_REQUIRE_MSG (f.data_type()==DataType::RealType,
      "Error! FieldPow only supports Real data type fields.\n"
      " - field name: " + fid.name() + "\n"
      " - field data type: " + e2str(f.data_type()) + "\n");

  const auto units = m_params.get<std::string>("units","check");
  EKAT_REQUIRE_MSG (units=="check" or units=="ignore",
      "Error! Invalid value for FieldPow's 'units'.\n"
      " - input: '" + units + "'\n"
      " - valid values: 'check', 'ignore'\n");

  const auto& u = fid.get_units();
  Units diag_units = none;
  // Dimensionless and unscaled (e.g., kg/kg) stays so, for any exponent
  if (units=="check" and (not u.is_dimensionless() or u.is_scaled())) {
    EKAT_REQUIRE_MSG (m_has_fraction,
        "Error! FieldPow: the units of X^p are u^p, which needs p to be a fraction n/d.\n"
        " - field name: " + fid.name() + "\n"
        " - field units: " + u.to_string() + "\n"
        " - exponent: " + m_exponent_str + "\n"
        " - write the exponent as 'n/d' (or as a decimal that is one, with d<=" +
          std::to_string(max_den) + "), or use units='ignore'.\n");
    diag_units = pow(u,ekat::RationalConstant(m_num,m_den));
  }

  // '/' is not welcome in a name
  auto exp_name = m_exponent_str;
  std::replace(exp_name.begin(),exp_name.end(),'/','_');

  auto& diag = add_output(FieldIdentifier(m_field_name + "_pow_" + exp_name,
                                          fid.get_layout().clone(), diag_units,
                                          fid.get_grid_name()));

  m_out_masked = f.has_valid_mask() or m_mask_invalid;
  if (m_mask_invalid) {
    diag.create_valid_mask();
  } else if (f.has_valid_mask()) {
    // The output is invalid exactly where the input is
    diag.set_valid_mask(f.get_valid_mask());
  }
  if (m_out_masked) {
    diag.get_header().set_may_be_filled(true);
  }
}

void FieldPow::compute_impl()
{
  using details::field_pfor;
  using details::view_as;

  const auto& f = m_fields_in.at(m_field_name);
  const auto  d = get();

  const PowSpec spec{m_kind,m_exponent,m_exponent_is_int};

  if (not m_out_masked) {
    field_pfor(PowKernel{spec},
               view_as<Real>(d), view_as<const Real>(f));
  } else {
    const auto& mask = d.get_valid_mask();
    if (m_mask_invalid) {
      if (f.has_valid_mask()) {
        field_pfor(DomainMaskAndKernel{spec},
                   view_as<int>(mask), view_as<const int>(f.get_valid_mask()),
                   view_as<const Real>(f));
      } else {
        field_pfor(DomainMaskKernel{spec},
                   view_as<int>(mask), view_as<const Real>(f));
      }
    }
    // Kernels on the same execution space run in order: no fence needed between them
    field_pfor(MaskedPowKernel{spec,constants::fill_value<Real>},
               view_as<Real>(d), view_as<const int>(mask), view_as<const Real>(f));
  }

  // field_pfor is asynchronous, but compute_impl must be blocking
  Kokkos::fence();
}

} // namespace scream
