#include "c_api_types.hpp"
#include <vector>
#include <sstream>
#include <string>

namespace {

std::vector<std::string> split_string(const std::string &original, const char delimiter) {
  std::vector<std::string> pieces;
  std::stringstream ss(original);
  std::string item;
  while (std::getline(ss, item, delimiter)) {
    pieces.push_back(item);
  }
  return pieces;
}

}


extern "C"{

}
