#include "coupler.hpp"
#include "field_registry.hpp"

#include <ekat_yaml.hpp>

#include <algorithm>
#include <filesystem>
#include <functional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace e3sm::emulator::coupler {

namespace {
template <typename Map> std::string join_keys(const Map& map) {
  std::string joined;
  for (const auto& [key, value] : map) {
    joined += (joined.empty() ? "" : ", ") + key;
  }
  return joined.empty() ? "none" : joined;
}

void require_unique_components(const CoupledFieldEntry& field,
                               const std::vector<std::string>& components,
                               const std::string& label) {
  for (auto it = components.begin(); it != components.end(); ++it) {
    if (std::find(components.begin(), it, *it) != it) {
      throw std::runtime_error("Route '" + field.id + "': component '" + *it +
                               "' is listed more than once in " + label);
    }
  }
}
} // namespace

/**
 * @brief Parse a YAML file containing the coupling field configuration.
 *
 * The YAML file must contain a top-level `fields` mapping. Each key is a route
 * ID, and its sublist names the field the route carries and describes its
 * coupling kind, merge type, sources, and destinations. Several routes can
 * carry the same field, such as a prognostic field leaving and returning to a
 * component. The file is parsed using ekat::YAMLParser into an
 * ekat::ParameterList.
 *
 * Expected YAML format:
 * fields:
 *   z:
 *     kind: state
 *     name: z
 *     merge_type: direct
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
    throw std::runtime_error("Coupling plan file not found: " +
                             filename.string());
  }

  try {
    const ekat::ParameterList params = ekat::parse_yaml_file(filename.string());
    // fields are the named coupling fields
    const auto& fields = params.sublist("fields");
    const auto& route_ids = fields.sublist_names();

    std::vector<CoupledFieldEntry> entries;
    entries.reserve(route_ids.size());

    for (const auto& route_id : route_ids) {
      // get the sub list for the route
      const auto& field = fields.sublist(route_id);
      const auto kind = field.get<std::string>("kind");
      auto field_name = field.get<std::string>("name");

      const auto merge_type_str = field.get<std::string>("merge_type");
      MergeType merge_type;
      if (merge_type_str == "direct") {
        merge_type = MergeType::Direct;
      } else if (merge_type_str == "scale") {
        merge_type = MergeType::ScaledByFraction;
      } else {
        throw std::runtime_error("Unknown merge type specified for " +
                                 route_id + ": '" + merge_type_str +
                                 "' (expected 'direct' or 'scale')");
      }

      const auto sources = field.get<std::vector<std::string>>("sources");
      const auto destinations =
          field.get<std::vector<std::string>>("destinations");

      entries.emplace_back(CoupledFieldEntry{
          .id = route_id,
          .field_name = std::move(field_name),
          .merge_type = merge_type,
          .sources = std::move(sources),
          .destinations = std::move(destinations),
      });
    }
    return entries;
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to read coupling plan '" +
                             filename.string() + "': " + e.what());
  }
}

std::string to_string(const CoupledFieldEntry& entry) {
  std::ostringstream out;

  out << entry.id << ":\n";
  out << "  name: " << entry.field_name << '\n';
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
  return out.str();
}

void Coupler::build_routes(const std::filesystem::path& filename) {
  auto fields = read_coupling_fields_from_yaml(filename);

  routes_.clear();
  routes_.reserve(fields.size());

  // A second route into the same input would silently overwrite the first
  std::unordered_map<FieldID, std::string> delivered_by;

  for (const auto& field : fields) {
    if (field.destinations.empty()) {
      throw std::runtime_error("Route '" + field.id + "' has no destinations");
    }
    require_unique_components(field, field.sources, "sources");
    require_unique_components(field, field.destinations, "destinations");

    CouplingRoute route{.merge_type = field.merge_type,
                        .name = field.id,
                        .field_name = field.field_name};
    for (const auto& component : field.sources) {
      if (!registry_.contains(component, field.field_name, FieldRole::Export) &&
          registry_.contains(component, field.field_name, FieldRole::Import)) {
        throw std::runtime_error("Route '" + field.id + "': '" + component +
                                 "' is listed as a source of '" +
                                 field.field_name +
                                 "', but it only imports this field. Sources "
                                 "must export the field.");
      }
      const auto id =
          registry_.get_id(component, field.field_name, FieldRole::Export);
      route.sources.push_back(id);
      coupling_plans_[component].export_ids.push_back(id);
    }
    for (const auto& component : field.destinations) {
      if (!registry_.contains(component, field.field_name, FieldRole::Import) &&
          registry_.contains(component, field.field_name, FieldRole::Export)) {
        throw std::runtime_error("Route '" + field.id + "': '" + component +
                                 "' is listed as a destination of '" +
                                 field.field_name +
                                 "', but it only exports this field. "
                                 "Destinations must import the field.");
      }
      const auto id =
          registry_.get_id(component, field.field_name, FieldRole::Import);
      const auto [delivered, is_new] = delivered_by.try_emplace(id, field.id);
      if (!is_new) {
        throw std::runtime_error("Routes '" + delivered->second + "' and '" +
                                 field.id + "' both deliver '" +
                                 field.field_name + "' to '" + component + "'");
      }
      route.destinations.push_back(id);
      coupling_plans_[component].import_ids.push_back(id);
    }
    routes_.emplace_back(std::move(route));
  }
}

const ActiveCouplingFields&
Coupler::coupling_plan(const std::string& component_name) const {
  const auto plan = coupling_plans_.find(component_name);
  if (plan == coupling_plans_.end()) {
    throw std::runtime_error("Component '" + component_name +
                             "' is not in the coupling plan. Components in "
                             "the plan: " +
                             join_keys(coupling_plans_));
  }
  return plan->second;
}

ComponentBuffers& Coupler::component_buffers(std::string_view component_name) {
  const auto buffers = buffers_.find(std::string(component_name));
  if (buffers == buffers_.end()) {
    throw std::runtime_error("Component '" + std::string(component_name) +
                             "' has no coupling buffers. Components with "
                             "buffers: " +
                             join_keys(buffers_));
  }
  return buffers->second;
}

/**
 * @brief: Create the buffers to allow transfer between components.
 *  If a route has no external_buffer, internal storage is allocated and used.
 */
void Coupler::build_buffers(RouteBufferMap external_buffers,
                            const ComponentIntervals& intervals) {
  buffers_.clear();

  for (auto& route : routes_) {
    if (route.merge_type != MergeType::Direct) {
      throw std::runtime_error("Route '" + route.name + "': merge type " +
                               to_string(route.merge_type) +
                               " is not supported yet; use 'direct'");
    }

    if (route.sources.size() != 1) {
      throw std::runtime_error("Route '" + route.name +
                               "': a direct merge requires exactly 1 source, "
                               "but " +
                               std::to_string(route.sources.size()) +
                               " are listed");
    }
    const auto source_id = route.sources.front();
    const auto& source_field = registry_.get(source_id);
    RouteBufferMap::iterator it = external_buffers.find(route.name);
    if (it == external_buffers.end()) {
      route.owned_buffer.resize(source_field.size);
      route.buffer = std::span(route.owned_buffer);
    } else {
      route.buffer = it->second;
      if (route.buffer.size() != source_field.size) {
        throw std::runtime_error(
            "Route '" + route.name + "': the external buffer has " +
            std::to_string(route.buffer.size()) + " value(s), but source '" +
            source_field.component + "' has " +
            std::to_string(source_field.size));
      }
    }

    route.source_component = source_field.component;
    route.destination_components.clear();
    route.averages = false;
    const auto source_interval = intervals.find(route.source_component);
    std::size_t n_slower_destinations = 0;
    for (const auto dest_id : route.destinations) {

      const auto& component = registry_.get(dest_id).component;
      route.destination_components.push_back(component);

      const auto dest_interval = intervals.find(component);
      if (source_interval != intervals.end() &&
          dest_interval != intervals.end() &&
          dest_interval->second > source_interval->second) {
        ++n_slower_destinations;
      }
    }
    if (n_slower_destinations > 0) {
      if (n_slower_destinations != route.destinations.size()) {
        throw std::runtime_error(
            "Route '" + route.name +
            "': some destinations run less often than source '" +
            route.source_component +
            "' and others don't, so they'd need different values. Use separate "
            "routes.");
      }
      route.averages = true;
      route.accumulator.assign(route.buffer.size(), 0.0);
      route.mean.assign(route.buffer.size(), 0.0);
      route.n_samples = 0;
    }

    auto& source_buffers = buffers_[source_field.component];
    source_buffers.exports.push_back({
        .id = source_id,
        .data = route.buffer,
    });

    for (const auto dest_id : route.destinations) {
      const auto& dest_field = registry_.get(dest_id);
      if (dest_field.size != route.buffer.size()) {
        throw std::runtime_error(
            "Route '" + route.name + "': destination '" + dest_field.component +
            "' has " + std::to_string(dest_field.size) +
            " value(s), but source '" + source_field.component + "' has " +
            std::to_string(route.buffer.size()));
      }
      auto& dest_buffer = buffers_[dest_field.component];

      dest_buffer.imports.push_back({
          .id = dest_id,
          .data = route.averages ? std::span<const double>(route.mean)
                                 : std::span<const double>(route.buffer),
      });
    }
  } // routes loop
}

void Coupler::accumulate_exports(std::string_view component_name) {
  for (auto& route : routes_) {
    if (!route.averages || route.source_component != component_name) {
      continue;
    }
    std::ranges::transform(route.accumulator, route.buffer,
                           route.accumulator.begin(), std::plus<>{});
    ++route.n_samples;
  }
}

void Coupler::prepare_imports(std::string_view component_name) {
  for (auto& route : routes_) {
    if (!route.averages ||
        std::ranges::find(route.destination_components, component_name) ==
            route.destination_components.end()) {
      continue;
    }
    if (route.n_samples == 0) {
      throw std::runtime_error(
          "Route '" + route.name + "': '" + std::string(component_name) +
          "' is importing a time mean, but source '" + route.source_component +
          "' has exported nothing since the last import");
    }
    const auto n_samples = static_cast<double>(route.n_samples);
    std::ranges::transform(route.accumulator, route.mean.begin(),
                           [n_samples](double sum) { return sum / n_samples; });
    std::ranges::fill(route.accumulator, 0.0);
    route.n_samples = 0;
  }
}

void Coupler::print_field_buffer(std::ostream& os, FieldID id,
                                 const std::span<const double> buffer,
                                 std::size_t count) const {
  const auto& field = registry_.get(id);

  os << "      " << field.component << "::" << field.name << " [id=" << id
     << ", size=" << field.size << ", buffer_size=" << buffer.size()
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

      os << "      " << field.name << " [id=" << id << "]\n";
    }

    os << "    imports:\n";
    for (const auto id : plan.import_ids) {
      const auto& field = registry_.get(id);

      os << "      " << field.name << " [id=" << id << "]\n";
    }
  }
}

} // namespace e3sm::emulator::coupler
