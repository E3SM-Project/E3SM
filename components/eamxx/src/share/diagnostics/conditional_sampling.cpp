#include "conditional_sampling.hpp"

#include "share/field/field_utils.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_team_policy_utils.hpp>

#include <charconv>
#include <cstdlib>
#include <string>

namespace scream {

namespace {
template<int N>
using MDRange = Kokkos::MDRangePolicy<
                  typename KokkosTypes<DefaultDevice>::ExeSpace,
                  Kokkos::Rank<N,Kokkos::Iterate::Right,Kokkos::Iterate::Right>
                >;
}

std::pair<int,bool> str2int (const std::string& s) {
  // Check if rhs is a value
  auto beg = s.data();
  auto end = beg + s.size();
  int i;
  auto [ptr,ec] = std::from_chars(beg,end,i);
  if (ec == std::errc{} && ptr == end) {
    return std::make_pair(i,true);
  }
  return std::make_pair(0,false);
}

// NOTE: the overload of from_chars for floating point was not immediately added
//       by compilers supporting C++17. Some compilers (e.g., on compy) do not
//       support it. Hence, stick with strtod.
std::pair<double,bool> str2dbl (const std::string& s) {
  if (s.empty()) return {0.0, false};

  char* endptr;
  const char* startptr = s.c_str();
  double d = std::strtod(startptr, &endptr);

  // Check if conversion happened and reached the end of the string
  if (endptr != startptr && *endptr == '\0') {
    return {d, true};
  }

  return {0.0, false};
}

// This is to avoid using original diag string in output field name,
// since user may have used >= and similar, which are not great in an nc var name
std::string cmp2str (const Comparison cmp)
{
  std::string s;
  switch (cmp) {
    case Comparison::EQ: s = "eq"; break;
    case Comparison::NE: s = "ne"; break;
    case Comparison::GT: s = "gt"; break;
    case Comparison::GE: s = "ge"; break;
    case Comparison::LT: s = "lt"; break;
    case Comparison::LE: s = "le"; break;
    default:
      EKAT_ERROR_MSG ("Unsupported/unrecognized Comparison enum value.\n");
  }
  return s;
}

Comparison str2cmp (const std::string& cmp_str)
{
  Comparison cmp = Comparison::EQ;
  if (cmp_str == "eq" || cmp_str == "==") {
    cmp = Comparison::EQ;
  } else if (cmp_str == "ne" || cmp_str == "!=") {
    cmp = Comparison::NE;
  } else if (cmp_str == "gt" || cmp_str == ">") {
    cmp = Comparison::GT;
  } else if (cmp_str == "ge" || cmp_str == ">=") {
    cmp = Comparison::GE;
  } else if (cmp_str == "lt" || cmp_str == "<") {
    cmp = Comparison::LT;
  } else if (cmp_str == "le" || cmp_str == "<=") {
    cmp = Comparison::LE;
  } else {
    EKAT_ERROR_MSG ("Error! Invalid comparison operator string.\n"
        " - cmp string   : " + cmp_str + "\n"
        " - valid choices: eq, ==, ne, !=, gt, >, ge, >=, lt, <, le, <=\n");
  }
  return cmp;
}

ConditionalSampling::
ConditionalSampling(const ekat::Comm &comm, const ekat::ParameterList &params)
 : AtmosphereDiagnostic(comm, params)
{
  m_input_f = m_params.get<std::string>("field_name");
  m_condition_lhs = m_params.get<std::string>("condition_lhs");
  m_condition_rhs = m_params.get<std::string>("condition_rhs");
  auto cmp = m_params.get<std::string>("condition_cmp");
  m_condition_cmp = str2cmp(cmp);

  m_lhs_is_lev = m_condition_lhs=="lev";

  auto rhs_as_double = str2dbl(m_condition_rhs);
  auto rhs_as_int    = str2int(m_condition_rhs);

  EKAT_REQUIRE_MSG (not m_lhs_is_lev or (rhs_as_int.second and rhs_as_int.first>=0),
      "Error! Conditional sampling of the form X_where_lev_CMP_Z requires Z to be a positive integer.\n"
      " - X  : " + m_input_f + "\n"
      " - CMP: " + cmp + "\n"
      " - Z  : " + m_condition_rhs + "\n");

  m_rhs_is_field = not (rhs_as_double.second or rhs_as_int.second);
  if (rhs_as_int.second)
    m_rhs_value = rhs_as_int.first;
  else if (rhs_as_double.second)
    m_rhs_value = rhs_as_double.first;

  m_diag_is_mask = m_input_f == "mask";

  m_diag_name = m_input_f + "_where_" + m_condition_lhs + "_" + cmp2str(m_condition_cmp) + "_" + m_condition_rhs;
}

void ConditionalSampling::create_requests()
{
  using namespace ShortFieldTagsNames;
  const auto &gn = m_params.get<std::string>("grid_name");
  const auto g   = m_grids_manager->get_grid("physics");

  // Special case: if the input field is "mask", we're just computing the mask where the condition holds
  if (not m_diag_is_mask) {
    add_field<Required>(m_input_f, gn);
  }

  // Special case: if lhs field is "lev", we don't need a lhs field.
  // We can actually precompute the output mask for a 1d col
  if (m_lhs_is_lev) {
    FieldIdentifier lev_fid("lev_mask",g->get_vertical_layout(LEV),ekat::units::none,g->name(),DataType::IntType);
    m_lev_mask = Field(lev_fid,true);
    auto lev_idx = m_lev_mask.clone("lev_idx");
    auto lev_idx_h = lev_idx.get_view<int*,Host>();
    for (int k=0; k<g->get_num_vertical_levels(); ++k)
      lev_idx_h(k) = k;
    lev_idx.sync_to_dev();
    compute_mask(lev_idx,m_rhs_value,m_condition_cmp,m_lev_mask);
  } else {
    add_field<Required>(m_condition_lhs, gn);
  }

  if (m_rhs_is_field) {
    add_field<Required>(m_condition_rhs, gn);
  }
}

void ConditionalSampling::initialize_impl(const RunType /*run_type*/)
{
  if (m_diag_is_mask) {
    if (m_lhs_is_lev) {
      m_diagnostic_output = m_lev_mask.clone(m_diag_name);
    } else {
      auto dfid = get_field_in(m_condition_lhs).get_header().get_identifier().clone(m_diag_name);
      dfid.reset_dtype(DataType::IntType).reset_units(ekat::units::none);
      m_diagnostic_output = Field(dfid,true);
    }
  } else {
    auto xfid = get_field_in(m_input_f).get_header().get_identifier();
    m_diagnostic_output = Field(xfid.clone(m_diag_name),true);
    m_diagnostic_output.create_valid_mask();
  }

  // If lhs is "lev" but diag is not "mask", we can still precompute the diag mask by
  // broadcasting m_lev_mask
  if (m_lhs_is_lev and not m_diag_is_mask) {
    auto mask = m_diagnostic_output.get_valid_mask();
    const auto& dims = mask.get_header().get_identifier().get_layout().dims();
    auto lev_mask_v = m_lev_mask.get_view<const int*>();
    switch (mask.rank()) {
      case 1:
        mask.deep_copy(m_lev_mask);
        break;
      case 2:
        {
          auto mask_v = mask.get_view<int**>();
          auto set_idx = KOKKOS_LAMBDA(int i, int k) {
            mask_v(i,k) = lev_mask_v(k);
          };
          MDRange<2> p({0,0},{dims[0],dims[1]});
          Kokkos::parallel_for(p,set_idx);
        } break;
      case 3:
        {
          auto mask_v = mask.get_view<int***>();
          auto set_idx = KOKKOS_LAMBDA(int i, int j, int k) {
            mask_v(i,j,k) = lev_mask_v(k);
          };
          MDRange<3> p({0,0,0},{dims[0],dims[1],dims[2]});
          Kokkos::parallel_for(p,set_idx);
        } break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported rank  in ConditionalSampling initialization.\n"
            " - diag name: " + m_diag_name + "\n"
            " - diag rank: " + std::to_string(mask.rank()) + "\n");
    }
  }
}

void ConditionalSampling::compute_diagnostic_impl()
{
  // Compute the mask (unless lhs is "lev", in which case we already did)
  auto& mask = m_diag_is_mask ? m_diagnostic_output : m_diagnostic_output.get_valid_mask();
  if (not m_lhs_is_lev) {
    const auto& lhs = get_field_in(m_condition_lhs);
    if (m_rhs_is_field) {
      const auto& rhs = get_field_in(m_condition_rhs);
      compute_mask(lhs,rhs,m_condition_cmp,mask);
      if (rhs.has_valid_mask())
        mask.scale(rhs.get_valid_mask());
    } else {
      compute_mask(lhs,m_rhs_value,m_condition_cmp,mask);
    }
    if (lhs.has_valid_mask())
      mask.scale(lhs.get_valid_mask());
  }

  if (not m_diag_is_mask) {
    const auto& f = get_field_in(m_input_f);
    m_diagnostic_output.deep_copy(f,mask);
    constexpr auto fv = constants::fill_value<Real>;
    m_diagnostic_output.deep_copy(fv,mask,true);
  }
}

} // namespace scream
