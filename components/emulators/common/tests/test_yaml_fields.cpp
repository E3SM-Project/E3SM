#include <ekat_yaml.hpp>

#include <iostream>
#include <string>

int main() {
  std::cout << "Parsing " << TEST_FIELDS_YAML << std::endl;
  const ekat::ParameterList params = ekat::parse_yaml_file(TEST_FIELDS_YAML);


  return 0;
}
