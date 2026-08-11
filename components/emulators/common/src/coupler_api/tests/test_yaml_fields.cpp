#include <coupler.hpp>
#include <iostream>
#include <string>

using namespace e3sm::coupler;
int main() {
  std::cout << "Parsing " << TEST_FIELDS_YAML << std::endl;

  Coupler coupler;
  const auto fields =
      coupler.read_coupling_fields_from_yaml(std::string(TEST_FIELDS_YAML));

  for (const auto& field : fields) {
    std::cout << field;
  }
  std::cout << std::endl;

  return 0;
}
