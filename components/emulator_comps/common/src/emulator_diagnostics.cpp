/**
 * @file emulator_diagnostics.cpp
 * @brief Implementation of diagnostic configuration utilities.
 */

#include "emulator_diagnostics.hpp"
#include <algorithm>
#include <cctype>

namespace emulator {

namespace {
// Helper to convert string to lowercase
std::string to_lower(const std::string &s) {
  std::string result = s;
  std::transform(result.begin(), result.end(), result.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return result;
}
} // namespace

FrequencyUnit str_to_freq_unit(const std::string &s) {
  std::string lower = to_lower(s);

  if (lower == "nsteps" || lower == "nstep") {
    return FrequencyUnit::NSTEPS;
  } else if (lower == "nsecs" || lower == "nseconds" || lower == "nsecond") {
    return FrequencyUnit::NSECS;
  } else if (lower == "nmins" || lower == "nminutes" || lower == "nminute") {
    return FrequencyUnit::NMINS;
  } else if (lower == "nhours" || lower == "nhour") {
    return FrequencyUnit::NHOURS;
  } else if (lower == "ndays" || lower == "nday") {
    return FrequencyUnit::NDAYS;
  } else if (lower == "nmonths" || lower == "nmonth") {
    return FrequencyUnit::NMONTHS;
  } else if (lower == "nyears" || lower == "nyear") {
    return FrequencyUnit::NYEARS;
  } else if (lower == "none" || lower == "never") {
    return FrequencyUnit::NONE;
  }

  return FrequencyUnit::NONE;
}

std::string freq_unit_to_str(FrequencyUnit u) {
  switch (u) {
  case FrequencyUnit::NSTEPS:
    return "nsteps";
  case FrequencyUnit::NSECS:
    return "nsecs";
  case FrequencyUnit::NMINS:
    return "nmins";
  case FrequencyUnit::NHOURS:
    return "nhours";
  case FrequencyUnit::NDAYS:
    return "ndays";
  case FrequencyUnit::NMONTHS:
    return "nmonths";
  case FrequencyUnit::NYEARS:
    return "nyears";
  case FrequencyUnit::NONE:
  default:
    return "none";
  }
}

OutputAvgType str_to_avg_type(const std::string &s) {
  std::string lower = to_lower(s);

  if (lower == "instant" || lower == "instantaneous") {
    return OutputAvgType::INSTANT;
  } else if (lower == "average" || lower == "avg" || lower == "mean") {
    return OutputAvgType::AVERAGE;
  } else if (lower == "min" || lower == "minimum") {
    return OutputAvgType::MIN;
  } else if (lower == "max" || lower == "maximum") {
    return OutputAvgType::MAX;
  } else if (lower == "std" || lower == "stddev" || lower == "stdev") {
    return OutputAvgType::STD;
  } else if (lower == "sum" || lower == "total" || lower == "accumulated") {
    return OutputAvgType::SUM;
  }

  return OutputAvgType::INSTANT;
}

std::string avg_type_to_str(OutputAvgType t) {
  switch (t) {
  case OutputAvgType::INSTANT:
    return "instant";
  case OutputAvgType::AVERAGE:
    return "average";
  case OutputAvgType::MIN:
    return "min";
  case OutputAvgType::MAX:
    return "max";
  case OutputAvgType::STD:
    return "std";
  case OutputAvgType::SUM:
    return "sum";
  default:
    return "instant";
  }
}

OutputPrecision str_to_precision(const std::string &s) {
  std::string lower = to_lower(s);

  if (lower == "float64" || lower == "double" || lower == "f64") {
    return OutputPrecision::FLOAT64;
  }

  return OutputPrecision::FLOAT32;
}

std::string precision_to_str(OutputPrecision p) {
  switch (p) {
  case OutputPrecision::FLOAT64:
    return "float64";
  case OutputPrecision::FLOAT32:
  default:
    return "float32";
  }
}

std::string file_type_suffix(FileType t) {
  switch (t) {
  case FileType::RESTART:
    return ".atm.r.";
  case FileType::HISTORY_RESTART:
    return ".atm.rh.";
  case FileType::HISTORY:
  default:
    return ".atm.h.";
  }
}

} // namespace emulator
