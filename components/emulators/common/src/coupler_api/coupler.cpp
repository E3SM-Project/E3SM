#include "field_registry.hpp"
#include <algorithm>
#include <coupler.hpp>
#include <ekat_yaml.hpp>
#include <filesystem>
#include <stdexcept>
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
Coupler::read_coupling_fields_from_yaml(const std::filesystem::path& filename) {

  if (!std::filesystem::exists(filename)) {
    throw std::runtime_error("Coupling configuration file not found: " +
                             filename.string());
  }

  const ekat::ParameterList params = ekat::parse_yaml_file(filename.string());
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
  return entries;
}

std::string to_string(const CoupledFieldEntry& entry) {
  std::ostringstream out;

  out << entry.id << ":\n";
  out << "  merge_type: " << to_string(entry.merge_type) << '\n';

  out << "  sources:";
  for (const auto& source : entry.sources) {
    out << ' ' << source;
  }
  out << '\n';

  out << "  destinations:";
  for (const auto& destination : entry.destinations) {
    out << ' ' << destination;
  }
  out << '\n';

  out << "  attributes:\n";
  out << "    attname: " << entry.attributes.name << '\n';
  out << "    longname: " << entry.attributes.long_name << '\n';
  out << "    stdname: " << entry.attributes.standard_name << '\n';
  out << "    units: " << entry.attributes.units << '\n';

  return out.str();
}

void Coupler::build_routes(const std::filesystem::path& filename) {
  auto fields = read_coupling_fields_from_yaml(filename);

  routes_.clear();
  routes_.reserve(fields.size());

  for (const auto& field : fields) {
    CouplingRoute route{.merge_type = field.merge_type};
    for (const auto& component : field.sources) {
      const auto id = registry_.get_id(component, field.id);
      route.sources.push_back(id);
      coupling_plans_[component].export_ids.push_back(id);
    }
    for (const auto& component : field.destinations) {
      const auto id = registry_.get_id(component, field.id);
      route.destinations.push_back(id);
      coupling_plans_[component].import_ids.push_back(id);
    }
    routes_.emplace_back(std::move(route));
  }
  build_buffers();
}

const ActiveCouplingFields&
Coupler::coupling_plan(const std::string& component_name) const {
  return coupling_plans_.at(component_name);
}

/**
 * @brief: Create the buffers to allow transfer between components.
 *
 * Currently only MergeType::Direct is supported, will implement others
 * once MOAB is better integrated
 */
void Coupler::build_buffers() {
  buffers_.clear();

  for (auto& route : routes_) {
    if (route.merge_type != MergeType::Direct) {
      throw std::runtime_error("Only direct merge type is currently supported");
    }

    if (route.sources.size() != 1) {
      throw std::runtime_error("Direct merge requires 1 source");
    }
    const auto source_id = route.sources.front();
    const auto& source_field = registry_.get(source_id);
    route.buffer.resize(source_field.size);

    auto& source_buffers = buffers_[source_field.component];
    source_buffers.exports.push_back({
        .id = source_id,
        .data = std::span<double>{route.buffer},
    });

    for (const auto dest_id : route.destinations) {
      const auto& dest_field = registry_.get(dest_id);
      auto& dest_buffer = buffers_[dest_field.component];

      dest_buffer.imports.push_back({
          .id = dest_id,
          .data = std::span<double>{route.buffer},
      });
    }
  } // routes loop
}

void Coupler::print_field_buffer(std::ostream& os, FieldID id,
                                 const std::span<const double> buffer,
                                 std::size_t count) const {
  const auto& field = registry_.get(id);

  os << "      " << field.component << "::" << field.attributes.name
     << " [id=" << id << ", size=" << field.size
     << ", buffer_size=" << buffer.size()
     << ", buffer=" << static_cast<const void*>(buffer.data()) << "]\n";

  if (count > 0) {
    os << "       [ ";
    for (std::size_t i = 0; i < std::min(count, buffer.size()); i++) {
      os << buffer[i] << ", ";
    }
    os << "]\n";
  }
}

void Coupler::print_coupling_plan(std::ostream& os) const {
  os << "Coupling routes:\n";

  for (std::size_t i = 0; i < routes_.size(); ++i) {
    const auto& route = routes_[i];

    os << "  Route " << i << '\n';
    os << "    merge: " << to_string(route.merge_type) << '\n';

    os << "    sources:\n";
    for (const auto id : route.sources) {
      print_field_buffer(os, id, route.buffer);
    }

    os << "    destinations:\n";
    for (const auto id : route.destinations) {
      print_field_buffer(os, id, route.buffer);
    }
  }

  os << "Buffers By Component\n";
  for (const auto& [component, bufs] : buffers_) {
    os << "  Component " << component << "\n";
    os << "    Exports: \n";
    for (const auto& buffer : bufs.exports) {
      print_field_buffer(os, buffer.id, buffer.data);
    }
    os << "    Imports: \n";
    for (const auto& buffer : bufs.imports) {
      print_field_buffer(os, buffer.id, buffer.data);
    }
  }

  os << "\nComponent plans:\n";

  for (const auto& [component, plan] : coupling_plans_) {
    os << "  " << component << '\n';

    os << "    exports:\n";
    for (const auto id : plan.export_ids) {
      const auto& field = registry_.get(id);

      os << "      " << field.attributes.name << " [id=" << id << "]\n";
    }

    os << "    imports:\n";
    for (const auto id : plan.import_ids) {
      const auto& field = registry_.get(id);

      os << "      " << field.attributes.name << " [id=" << id << "]\n";
    }
  }
}

} // namespace e3sm::coupler
