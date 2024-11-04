#ifndef SCREAM_FIELD_UTILS_HPP
#define SCREAM_FIELD_UTILS_HPP

#include "share/field/field_utils_impl.hpp"

namespace scream {

// Check that two fields store the same entries.
// NOTE: if the field is padded, padding entries are NOT checked.
inline bool views_are_equal(const Field& f1, const Field& f2, const ekat::Comm* comm = nullptr) {
  EKAT_REQUIRE_MSG (f1.data_type()==f2.data_type(),
      "Error! Views have different data types.\n");

  switch (f1.data_type()) {
    case DataType::IntType:
      return impl::views_are_equal<const int>(f1,f2,comm);
    case DataType::FloatType:
      return impl::views_are_equal<const float>(f1,f2,comm);
    case DataType::DoubleType:
      return impl::views_are_equal<const double>(f1,f2,comm);
    default:
      EKAT_ERROR_MSG ("Error! Unrecognized field data type.\n");
  }
}

template<typename Engine, typename PDF>
void randomize (const Field& f, Engine& engine, PDF&& pdf)
{
  EKAT_REQUIRE_MSG(f.is_allocated(),
      "Error! Cannot randomize the values of a field not yet allocated.\n");

  // Deduce scalar type from pdf
  using ST = decltype(pdf(engine));

  // Check compatibility between PDF and field data type
  const auto data_type = f.data_type();
  EKAT_REQUIRE_MSG (
      (std::is_same<ST,int>::value && data_type==DataType::IntType) ||
      (std::is_same<ST,float>::value && data_type==DataType::FloatType) ||
      (std::is_same<ST,double>::value && data_type==DataType::DoubleType),
      "Error! Field data type incompatible with input PDF.\n");

  impl::randomize<ST>(f,engine,pdf);
}

// Compute a random perturbation of a field for all view entries whose
// level index satisfies the mask.
// Input:
//   - f:                  Field to perturbed. Required to have level midpoint
//                         tag as last dimension.
//   - engine:             Random number engine.
//   - pdf:                Random number distribution where a random value (say,
//                         pertval) is taken s.t.
//                           field_view(i0,...,iN) *= (1 + pertval)
//   - base_seed:          Seed used for creating the engine input.
//   - level_mask:         Mask (size of the level dimension of f) where f(i0,...,k) is
//                         perturbed if level_mask(k)=true
//   - dof_gids:           Field containing global DoF IDs for columns of f (if applicable)
template<typename Engine, typename PDF, typename MaskType>
void perturb (const Field& f,
              Engine& engine,
              PDF&& pdf,
              const int base_seed,
              const MaskType& level_mask,
              const Field& dof_gids = Field())
{
  EKAT_REQUIRE_MSG(f.is_allocated(),
       	           "Error! Cannot perturb the values of a field not yet allocated.\n");

  // Deduce scalar type from pdf
  using ST = decltype(pdf(engine));

  // Check compatibility between PDF and field data type
  const auto data_type = f.data_type();
  EKAT_REQUIRE_MSG((std::is_same_v<ST,int> && data_type==DataType::IntType) or
                   (std::is_same_v<ST,float> && data_type==DataType::FloatType) or
                   (std::is_same_v<ST,double> && data_type==DataType::DoubleType),
                   "Error! Field data type incompatible with input PDF.\n");

  using namespace ShortFieldTagsNames;
  const auto& fl = f.get_header().get_identifier().get_layout();

  // Field we are perturbing should have a level dimension,
  // and it is required to be the last dimension
  EKAT_REQUIRE_MSG(fl.rank()>0 &&
                   (fl.tags().back() == LEV || fl.tags().back() == ILEV),
                   "Error! Trying to perturb field \""+f.name()+"\", but field "
	                 "does not have LEV or ILEV as last dimension.\n"
                   "  - field name: " + f.name() + "\n"
                   "  - field layout: " + fl.to_string() + "\n");

  if (fl.has_tag(COL)) {
    // If field has a column dimension, it should be the first dimension
    EKAT_REQUIRE_MSG(fl.tag(0) == COL,
                     "Error! Trying to perturb field \""+f.name()+"\", but field "
	                   "does not have COL as first dimension.\n"
                     "  - field name: " + f.name() + "\n"
                     "  - field layout: " + fl.to_string() + "\n");

    const auto& dof_gids_fl = dof_gids.get_header().get_identifier().get_layout();
    EKAT_REQUIRE_MSG(dof_gids_fl.dim(0) == fl.dim(COL),
                     "Error! Field of DoF GIDs should have the same size as "
                     "perturbed field's column dimension.\n"
                     "  - dof_gids dim: " + std::to_string(dof_gids_fl.dim(0)) + "\n"
                     "  - field name: " + f.name() + "\n"
                     "  - field layout: " + fl.to_string() + "\n");
    EKAT_REQUIRE_MSG(dof_gids.data_type() == DataType::IntType,
                     "Error! DoF GIDs field must have \"int\" as data type.\n");
  }

  impl::perturb<ST>(f, engine, pdf, base_seed, level_mask, dof_gids);
}

template<typename ST>
ST frobenius_norm(const Field& f, const ekat::Comm* comm = nullptr)
{
  // Check compatibility between ST and field data type
  const auto data_type = f.data_type();
  EKAT_REQUIRE_MSG (data_type==DataType::FloatType || data_type==DataType::DoubleType,
      "Error! Frobenius norm only allowed for floating-point field value types.\n");

  EKAT_REQUIRE_MSG (
      (std::is_same<ST,float>::value && data_type==DataType::FloatType) ||
      (std::is_same<ST,double>::value && data_type==DataType::DoubleType),
      "Error! Field data type incompatible with template argument.\n");

  return impl::frobenius_norm<ST>(f,comm);
}

template<typename ST>
ST field_sum(const Field& f, const ekat::Comm* comm = nullptr)
{
  // Check compatibility between ST and field data type
  const auto data_type = f.get_header().get_identifier().data_type();

  EKAT_REQUIRE_MSG (
      (std::is_same<ST,int>::value && data_type==DataType::IntType) ||
      (std::is_same<ST,float>::value && data_type==DataType::FloatType) ||
      (std::is_same<ST,double>::value && data_type==DataType::DoubleType),
      "Error! Field data type incompatible with template argument.\n");

  return impl::field_sum<ST>(f,comm);
}

template<typename ST>
ST field_max(const Field& f, const ekat::Comm* comm = nullptr)
{
  // Check compatibility between ST and field data type
  const auto data_type = f.data_type();

  EKAT_REQUIRE_MSG (
      (std::is_same<ST,int>::value && data_type==DataType::IntType) ||
      (std::is_same<ST,float>::value && data_type==DataType::FloatType) ||
      (std::is_same<ST,double>::value && data_type==DataType::DoubleType),
      "Error! Field data type incompatible with template argument.\n");

  return impl::field_max<ST>(f,comm);
}

template<typename ST>
ST field_min(const Field& f, const ekat::Comm* comm = nullptr)
{
  // Check compatibility between ST and field data type
  const auto data_type = f.data_type();

  EKAT_REQUIRE_MSG (
      (std::is_same<ST,int>::value && data_type==DataType::IntType) ||
      (std::is_same<ST,float>::value && data_type==DataType::FloatType) ||
      (std::is_same<ST,double>::value && data_type==DataType::DoubleType),
      "Error! Field data type incompatible with template argument.\n");

  return impl::field_min<ST>(f,comm);
}

// Prints the value of a field at a certain location, specified by tags and indices.
// If the field layout contains all the location tags, we will slice the field along
// those tags, and print it. E.g., f might be a <COL,LEV> field, and the tags/indices
// refer to a single column, in which case we'll print a whole column worth of data.
inline void
print_field_hyperslab (const Field& f,
                       const std::vector<FieldTag>& tags = {},
                       const std::vector<int>& indices = {},
                       std::ostream& out = std::cout)
{
  const auto dt = f.data_type();
  const auto rank = f.rank();

  EKAT_REQUIRE_MSG (rank>=static_cast<int>(tags.size()),
      "Error! Requested location incompatible with field rank.\n"
      "  - field name: " + f.name() + "\n"
      "  - field rank: " + std::to_string(rank) + "\n"
      "  - requested indices: (" + ekat::join(indices,",") + "\n");

  switch (dt) {
    case DataType::IntType:
      impl::print_field_hyperslab<int>(f,tags,indices,out,rank,0);
      break;
    case DataType::FloatType:
      impl::print_field_hyperslab<float>(f,tags,indices,out,rank,0);
      break;
    case DataType::DoubleType:
      impl::print_field_hyperslab<double>(f,tags,indices,out,rank,0);
      break;
    default:
      EKAT_ERROR_MSG ("[print_field_hyperslab] Error! Invalid/unsupported data type.\n"
          " - field name: " + f.name() + "\n");
  }
}

} // namespace scream

#endif // SCREAM_FIELD_UTILS_HPP
