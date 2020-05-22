#ifndef SCREAM_PARSE_YAML_FILE_HPP
#define SCREAM_PARSE_YAML_FILE_HPP

#include "ekat/scream_parameter_list.hpp"
#include <string>

namespace scream {

ParameterList parse_yaml_file (const std::string& fname);
void parse_yaml_file (const std::string& fname, ParameterList& params);

} // namespace scream

#endif // SCREAM_PARSE_YAML_FILE_HPP
