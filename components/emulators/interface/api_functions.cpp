#include "api_functions.hpp"
#include "harness.hpp"

#include <field_types.hpp>

#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace e3sm::emulator {

namespace {
std::vector<std::string> split_string(const std::string& original,
                                      const char delimiter) {
  std::vector<std::string> pieces;
  std::stringstream ss(original);
  std::string item;
  while (std::getline(ss, item, delimiter)) {
    pieces.push_back(item);
  }
  return pieces;
}

ScalarType dtype_from_fortran(int code) {
  switch (code) {
  case 0:
    return ScalarType::Real64;
  case 1:
    return ScalarType::Int32;
  }
}

} // namespace

extern "C" {
void* harness_create() {
  try {
    return new Harness{};
  } catch (...) {
    return nullptr;
  }
}

void harness_destroy(void* handle) { delete static_cast<Harness*>(handle); }

void harness_register_field(void* handle, const char* name, int extent,
                            int dtype, void* data) {
  Harness* harness = static_cast<Harness*>(handle);
  FieldView view =
      FieldView{make_field_view(dtype_from_fortran(dtype), extent, data)};
  harness->register_field(std::string{name}, view);
}
} // extern "c"
} // namespace e3sm::emulator
