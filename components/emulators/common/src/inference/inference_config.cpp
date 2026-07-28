/**
 * @file inference_config.cpp
 * @brief InferenceConfig implementation.
 */

#include "inference_config.hpp"

#include "inference_error.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>

namespace emulator {
namespace inference {

namespace {

std::string trim(const std::string &s) {
  const auto first = s.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    return "";
  }
  return s.substr(first, s.find_last_not_of(" \t\r\n") - first + 1);
}

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

bool parse_bool(const std::string &value, const std::string &key) {
  const std::string v = lower(trim(value));
  if (v == "true" || v == "1" || v == "yes" || v == "on" || v == ".true.") {
    return true;
  }
  if (v == "false" || v == "0" || v == "no" || v == "off" || v == ".false.") {
    return false;
  }
  EMULATOR_INFER_REQUIRE(false, "Could not read '" << value << "' for '" << key
                                                   << "' as a boolean.");
  return false;
}

int parse_int(const std::string &value, const std::string &key) {
  try {
    return std::stoi(trim(value));
  } catch (const std::exception &) {
    EMULATOR_INFER_REQUIRE(false, "Could not read '"
                                      << value << "' for '" << key
                                      << "' as an integer.");
  }
  return 0;
}

/// Strip an optional pair of matching quotes (namelists often carry them).
std::string unquote(const std::string &s) {
  if (s.size() >= 2 && (s.front() == '"' || s.front() == '\'') &&
      s.back() == s.front()) {
    return s.substr(1, s.size() - 2);
  }
  return s;
}

} // namespace

void InferenceConfig::apply(const std::string &key, const std::string &raw) {
  const std::string k = lower(trim(key));
  const std::string value = unquote(trim(raw));

  if (k == "backend") {
    backend = lower(value);
  } else if (k == "model_path" || k == "model") {
    model_path = value;
  } else if (k == "input") {
    inputs.push_back(value);
  } else if (k == "output") {
    outputs.push_back(value);
  } else if (k == "input_channels") {
    input_channels = parse_int(value, k);
  } else if (k == "output_channels") {
    output_channels = parse_int(value, k);
  } else if (k == "verbose") {
    verbose = parse_bool(value, k);
  } else if (k.rfind("option.", 0) == 0) {
    options[k.substr(7)] = value;
  } else {
    // Unknown keys are options rather than errors: the Python side owns its
    // own settings and should not need a C++ change to add one.
    options[k] = value;
  }
}

std::string InferenceConfig::get(const std::string &key,
                                 const std::string &fallback) const {
  const auto it = options.find(key);
  return it == options.end() ? fallback : it->second;
}

int InferenceConfig::get_int(const std::string &key, int fallback) const {
  const auto it = options.find(key);
  return it == options.end() ? fallback : parse_int(it->second, key);
}

bool InferenceConfig::get_bool(const std::string &key, bool fallback) const {
  const auto it = options.find(key);
  return it == options.end() ? fallback : parse_bool(it->second, key);
}

InferenceConfig InferenceConfig::from_file(const std::string &path,
                                           const std::string &prefix) {
  InferenceConfig config;
  std::ifstream ifs(path);
  if (!ifs) {
    return config; // no file, no settings: the defaults are a valid run
  }

  std::string line;
  while (std::getline(ifs, line)) {
    const auto comment = line.find_first_of("#!");
    if (comment != std::string::npos) {
      line = line.substr(0, comment);
    }

    const auto colon = line.find(':');
    if (colon == std::string::npos) {
      continue; // not a setting; the file may hold other things
    }
    std::string key = trim(line.substr(0, colon));
    const std::string value = line.substr(colon + 1);

    if (!prefix.empty()) {
      if (lower(key).rfind(lower(prefix), 0) != 0) {
        continue; // belongs to the component, not to us
      }
      key = key.substr(prefix.size());
    }
    if (trim(key).empty()) {
      continue;
    }
    config.apply(key, value);
  }
  return config;
}

} // namespace inference
} // namespace emulator
