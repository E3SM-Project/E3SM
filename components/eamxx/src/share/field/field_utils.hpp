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
                       const std::vector<FieldTag>& tags,
                       const std::vector<int>& indices,
                       std::ostream& out)
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
  }
}

} // namespace scream

#endif // SCREAM_FIELD_UTILS_HPP
