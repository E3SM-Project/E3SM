/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_ERRORDEFS_HPP
#define HOMMEXX_ERRORDEFS_HPP

#ifndef NDEBUG
#define DEBUG_PRINT(...) \
  do {                   \
    printf(__VA_ARGS__); \
  } while(false)
// This macro always evaluates eval, but
// This enables us to define variables specifically for use
// in asserts Note this can still cause issues
#define DEBUG_EXPECT(eval, expected) \
do {                                 \
  auto v = eval;                     \
  assert(v == expected);             \
} while(false)
#else
#define DEBUG_PRINT(...) \
do {                     \
} while(false)
#define DEBUG_EXPECT(eval, expected) \
do {                                 \
  eval;                              \
} while(false)
#endif

#ifdef DEBUG_TRACE
#define TRACE_PRINT(...) \
do {                     \
  printf(__VA_ARGS__);   \
} while(false)
#else
#define TRACE_PRINT(...) \
do {                     \
} while(false)
#endif

#include <HommexxEnums.hpp>
#include <string>
#include <sstream>

namespace Homme {
namespace Errors {

void runtime_check(bool cond, const std::string& message, int code);
void runtime_abort(const std::string& message, int code);

static constexpr int err_unknown_option               = 11;
static constexpr int err_not_implemented              = 12;
static constexpr int err_invalid_options_combination  = 13;
static constexpr int err_negative_layer_thickness     = 101;

template<typename T>
void option_error(const std::string& location,
                  const std::string& option,
                  const T& value)
{
  std::stringstream msg;
  msg << "Error in " << location << ": " << "unsupported value '"
      << value << "' for input parameter '" << option << "'.";

  runtime_abort(msg.str(),err_not_implemented);
}

template<typename T>
void check_option(const std::string& location,
                  const std::string& option,
                  const T& actual_value,
                  const std::initializer_list<T>& admissible_values)
{
  bool bad_value = true;
  for (const auto& value : admissible_values) {
    if (value==actual_value) {
      bad_value = false;
    }
  }

  if (bad_value) {
    option_error(location,option,actual_value);
  }
}

template<typename T>
void check_option (const std::string& location,
                   const std::string& option,
                   const T& value, const T& ref_value,
                   const ComparisonOp& relation)
{
  bool bad_inputs = false;
  switch (relation) {
    case ComparisonOp::EQ: if (value!=ref_value) { bad_inputs = true; } break;
    case ComparisonOp::NE: if (value==ref_value) { bad_inputs = true; } break;
    case ComparisonOp::GT: if (value<=ref_value) { bad_inputs = true; } break;
    case ComparisonOp::LT: if (value>=ref_value) { bad_inputs = true; } break;
    case ComparisonOp::GE: if (value<ref_value)  { bad_inputs = true; } break;
    case ComparisonOp::LE: if (value>ref_value)  { bad_inputs = true; } break;
  }

  if (bad_inputs) {
    const std::string cmp_str[6] {"==", "!=", ">", "<", ">=", "<="};
    std::stringstream msg;
    msg << "Error in " << location << ": " << "unsupported value '"
        << value << "' for input parameter '" << option << "'.\n"
        << " The value should satisfy " << option << " " << cmp_str[etoi(relation)]
        << " " << ref_value << ".";

    runtime_abort(msg.str(),err_invalid_options_combination);
  }
}

template<typename T1, typename T2>
void check_options_relation(const std::string& location,
                            const std::string& option1,
                            const std::string& option2,
                            const T1& value1, const T2& value2,
                            const ComparisonOp& relation)
{
  bool bad_inputs = false;
  switch (relation) {
    case ComparisonOp::EQ: if (value1!=value2) { bad_inputs = true; } break;
    case ComparisonOp::NE: if (value1==value2) { bad_inputs = true; } break;
    case ComparisonOp::GT: if (value1<=value2) { bad_inputs = true; } break;
    case ComparisonOp::LT: if (value1>=value2) { bad_inputs = true; } break;
    case ComparisonOp::GE: if (value1<value2)  { bad_inputs = true; } break;
    case ComparisonOp::LE: if (value1>value2)  { bad_inputs = true; } break;
  }

  if (bad_inputs) {
    const std::string cmp_str[6] {"==", "!=", ">", "<", ">=", "<="};
    std::stringstream msg;
    msg << "Error in " << location << ": " << "unsupported combination for input parameters '"
        << option1 << "' (" << value1 << ") and '" << option2 << "' (" << value2 << ").\n"
        << " The two should satisfy " << option1 << " " << cmp_str[etoi(relation)]
        << " " << option2 << ".";

    runtime_abort(msg.str(),err_invalid_options_combination);
  }
}

} // namespace Errors
} // namespace Homme

#endif // HOMMEXX_ERRORDEFS_HPP
