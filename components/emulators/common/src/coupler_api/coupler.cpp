#include "field_registry.hpp"
#include <coupler.hpp>
#include <ekat_yaml.hpp>
#include <filesystem>
#include <vector>

namespace e3sm::coupler {

/**
 * @brief Parse a YAML file containing the coupling field configuration.
 *
 * The YAML file must contain a top-level `fields` mapping. Each field is a
 * named sublist describing its coupling kind, merge type, sources,
 * destinations, and metadata. The file is parsed using ekat::YAMLParser
 * into an ekat::ParameterList.
 *
 * Expected YAML format:
 * fields:
 *   z:
 *     kind: state
 *     name: z
 *     merge_type: direct
 *     attributes:
 *       attname: sa_z
 *       longname: height at the lowest model level
 *       stdname: height
 *       units: m
 *     sources:
 *     - atm
 *     destinations:
 *     - lnd
 *     - ice
 *
 * @param filename: Path to the coupling configuration YAML file.
 * @return A vector of coupling field descriptions parsed from the file.
 *
 * @throws std::runtime_error If the configuration file does not exist or
 *         cannot be parsed.
 */
std::vector<CoupledFieldEntry>
read_coupling_fields_from_yaml(const std::string& filename) {

  if (!std::filesystem::exists(filename)) {
    throw std::runtime_error("Coupling configuration file not found: " +
                             filename);
  }

  const ekat::ParameterList params = ekat::parse_yaml_file(filename);
  // fields are the named coupling fields
  const auto& fields = params.sublist("fields");
  const auto& field_names = fields.sublist_names();

  std::vector<CoupledFieldEntry> entries;
  entries.reserve(field_names.size());

  for (const auto& field_name : field_names) {
    // get the sub list for the field
    const auto& field = fields.sublist(field_name);

    const auto kind = field.get<std::string>("kind");

    const auto merge_type_str = field.get<std::string>("merge_type");

    MergeType merge_type;
    if (merge_type_str == "direct") {
      merge_type = MergeType::Direct;
    } else if (merge_type_str == "scale") {
      merge_type = MergeType::ScaledByFraction;
    } else {
      std::runtime_error("Unknow merged type specificed for " + field_name +
                         ":" + merge_type_str);
    }

    const auto& attributes = field.sublist("attributes");
    const auto attname = attributes.get<std::string>("attname");
    const auto longname = attributes.get<std::string>("longname");
    const auto stdname = attributes.get<std::string>("stdname");
    const auto units = attributes.get<std::string>("units");

    const auto sources = field.get<std::vector<std::string>>("sources");

    const auto destinations =
        field.get<std::vector<std::string>>("destinations");

    entries.emplace_back(CoupledFieldEntry{
        .id = field_name,
        .merge_type = merge_type,
        .attributes = {.name = field_name,
                       .long_name = longname,
                       .standard_name = stdname,
                       .units = units},
        .sources = std::move(sources),
        .destinations = std::move(destinations),
    });
  }
}
} // namespace e3sm::coupler
