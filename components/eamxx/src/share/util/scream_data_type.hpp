#ifndef SCREAM_DATA_TYPE_HPP
#define SCREAM_DATA_TYPE_HPP

#include "share/util/scream_utils.hpp"

#include <ekat/ekat_assert.hpp>

#include <string>

namespace scream
{

// An enum for specifying fields data type
// WARNING: if you add new data types, ensure that the underlying
//          int value is such that, if i<j, then there is NO narrowing
//          conversion going from i to j. If that's not possible, then
//          you need to change the impl of is_narrowing_conversion below.
enum class DataType {
  Invalid    = 0,
  IntType    = 1,
  FloatType  = 2,
  DoubleType = 3,
#ifdef SCREAM_DOUBLE_PRECISION
  RealType   = DoubleType
#else
  RealType   = FloatType
#endif
};

template<typename ST>
DataType get_data_type () {
  // if statements are compiled out
  if (std::is_same<ST,int>::value) {
    return DataType::IntType;
  } else if (std::is_same<ST,float>::value) {
    return DataType::FloatType;
  } else if (std::is_same<ST,double>::value) {
    return DataType::DoubleType;
  } else {
    EKAT_ERROR_MSG ("Error! Unsupported data type.\n"
        " - typeid(ST): " + std::string(typeid(ST).name()) + "\n");
  }
  return DataType::Invalid;
}

inline bool is_narrowing_conversion (const DataType from, const DataType to) {
  auto from_i = etoi(from);
  auto to_i = etoi(to);
  return from_i>to_i;
}

inline std::string e2str (const DataType data_type) {
  switch (data_type) {
    case DataType::IntType:    return "int";
    case DataType::FloatType:  return "float";
    case DataType::DoubleType: return "double";
    default:
      EKAT_ERROR_MSG("Error! Unsupported DataType value.\n");
  }
}

inline int get_type_size (const DataType data_type) {
  switch (data_type) {
    case DataType::IntType:    return sizeof(int);
    case DataType::FloatType:  return sizeof(float);
    case DataType::DoubleType: return sizeof(double);
    default:
      EKAT_ERROR_MSG("Error! Unsupported DataType value.\n");
  }
}

} // namespace scream

#endif // SCREAM_DATA_TYPE_HPP
