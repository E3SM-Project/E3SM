/**
 * @file emulator_config.cpp
 * @brief Implementation of configuration parsing functions.
 *
 * Parses YAML configuration files using yaml-cpp to populate
 * EmulatorConfig structures.
 *
 * @see emulator_config.hpp for structure definitions
 */

#include "emulator_config.hpp"
#include "emulator_diagnostics.hpp"
#include <fstream>
#include <yaml-cpp/yaml.h>

namespace emulator {

/**
 * @brief Parse configuration from a YAML file.
 *
 * Reads the YAML file and extracts configuration for the specified
 * component section (e.g., "eatm"). All subsections (build, runtime,
 * model_io, coupling) are parsed if present.
 */
EmulatorConfig parse_emulator_config(const std::string &yaml_file,
                                     const std::string &section_name) {
  EmulatorConfig config;

  YAML::Node root = YAML::LoadFile(yaml_file);

  if (!root[section_name]) {
    throw std::runtime_error("YAML config missing '" + section_name +
                             "' section: " + yaml_file);
  }

  YAML::Node section = root[section_name];

  // Parse build section (optional - may be set by buildnml)
  if (section["build"]) {
    YAML::Node build = section["build"];
    if (build["grid_name"]) {
      config.build.grid_name = build["grid_name"].as<std::string>();
    }
    if (build["inference_backend"]) {
      config.build.inference_backend =
          build["inference_backend"].as<std::string>();
    }
  }

  // Parse runtime section
  if (section["runtime"]) {
    YAML::Node runtime = section["runtime"];
    if (runtime["model_path"]) {
      config.runtime.model_path = runtime["model_path"].as<std::string>();
    }
    if (runtime["ic_file"]) {
      config.runtime.ic_file = runtime["ic_file"].as<std::string>();
    }
    if (runtime["enabled"]) {
      config.runtime.enabled = runtime["enabled"].as<bool>();
    }
  }

  // Parse model I/O section
  if (section["model_io"]) {
    YAML::Node model_io = section["model_io"];
    if (model_io["input_variables"]) {
      for (const auto &var : model_io["input_variables"]) {
        config.model_io.input_variables.push_back(var.as<std::string>());
      }
    }
    if (model_io["output_variables"]) {
      for (const auto &var : model_io["output_variables"]) {
        config.model_io.output_variables.push_back(var.as<std::string>());
      }
    }
    if (model_io["spatial_mode"]) {
      config.model_io.spatial_mode = model_io["spatial_mode"].as<bool>();
    }
  }

  // Parse coupling section
  if (section["coupling"]) {
    YAML::Node coupling = section["coupling"];
    if (coupling["zero_init_exports"]) {
      config.coupling.zero_init_exports =
          coupling["zero_init_exports"].as<bool>();
    }
    if (coupling["debug"]) {
      config.coupling.debug = coupling["debug"].as<bool>();
    }
  }

  // Parse diagnostics section
  if (section["diagnostics"]) {
    YAML::Node diag = section["diagnostics"];

    // Parse history streams
    if (diag["history"]) {
      for (const auto &stream_node : diag["history"]) {
        OutputStreamConfig stream_cfg;

        if (stream_node["stream_name"]) {
          stream_cfg.stream_name = stream_node["stream_name"].as<std::string>();
        }
        if (stream_node["filename_prefix"]) {
          stream_cfg.filename_prefix =
              stream_node["filename_prefix"].as<std::string>();
        }
        if (stream_node["frequency"]) {
          stream_cfg.frequency = stream_node["frequency"].as<int>();
        }
        if (stream_node["frequency_unit"]) {
          stream_cfg.frequency_unit =
              str_to_freq_unit(stream_node["frequency_unit"].as<std::string>());
        }
        if (stream_node["averaging"]) {
          stream_cfg.avg_type =
              str_to_avg_type(stream_node["averaging"].as<std::string>());
        }
        if (stream_node["precision"]) {
          stream_cfg.precision =
              str_to_precision(stream_node["precision"].as<std::string>());
        }
        if (stream_node["max_snapshots_per_file"]) {
          stream_cfg.max_snapshots_per_file =
              stream_node["max_snapshots_per_file"].as<int>();
        }
        if (stream_node["fields"]) {
          for (const auto &field : stream_node["fields"]) {
            stream_cfg.fields.push_back(field.as<std::string>());
          }
        }

        config.diagnostics.history_streams.push_back(stream_cfg);
      }
    }

    // Parse restart config
    if (diag["restart"]) {
      YAML::Node restart = diag["restart"];
      if (restart["enabled"]) {
        config.diagnostics.restart.enabled = restart["enabled"].as<bool>();
      }
      if (restart["frequency"]) {
        config.diagnostics.restart.frequency = restart["frequency"].as<int>();
      }
      if (restart["frequency_unit"]) {
        config.diagnostics.restart.frequency_unit =
            str_to_freq_unit(restart["frequency_unit"].as<std::string>());
      }
      if (restart["filename_prefix"]) {
        config.diagnostics.restart.filename_prefix =
            restart["filename_prefix"].as<std::string>();
      }
    }

    // Parse history restart config
    if (diag["history_restart"]) {
      YAML::Node hrst = diag["history_restart"];
      if (hrst["enabled"]) {
        config.diagnostics.history_restart.enabled = hrst["enabled"].as<bool>();
      }
      if (hrst["filename_prefix"]) {
        config.diagnostics.history_restart.filename_prefix =
            hrst["filename_prefix"].as<std::string>();
      }
    }
  }

  // Legacy: grid_file at top level
  if (section["grid_file"]) {
    config.grid_file = section["grid_file"].as<std::string>();
  }

  return config;
}

/**
 * @brief Parse configuration with fallback to defaults.
 *
 * Gracefully handles missing or malformed config files by returning
 * default configuration values.
 */
EmulatorConfig
parse_emulator_config_with_defaults(const std::string &yaml_file,
                                    const std::string &section_name,
                                    bool verbose) {
  (void)verbose; // Currently unused (removed logging)

  EmulatorConfig config;

  // Check if file exists
  std::ifstream test(yaml_file);
  if (!test.is_open()) {
    return config;
  }
  test.close();

  try {
    config = parse_emulator_config(yaml_file, section_name);
  } catch (...) {
    // Silently use defaults on parse error
  }

  return config;
}

} // namespace emulator
