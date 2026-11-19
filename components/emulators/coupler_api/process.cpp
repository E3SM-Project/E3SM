#include "process.hpp"
#include "field_registry.hpp"
#include "field_types.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>

namespace e3sm::emulator::coupler {

namespace {
void check_buffer_count(const char* kind, std::size_t expected,
                        std::size_t actual) {
  if (expected != actual) {
    throw std::runtime_error(std::string("Incorrect number of ") + kind +
                             " buffers: expected " + std::to_string(expected) +
                             ", got " + std::to_string(actual));
  }
}

void check_buffer_match(const char* kind, const ExportBuffer& local, FieldID id,
                        std::size_t size) {
  if (local.id != id) {
    throw std::runtime_error(
        std::string(kind) + " field ID mismatch: expected " +
        std::to_string(local.id) + ", got " + std::to_string(id));
  }
  if (local.data.size() != size) {
    throw std::runtime_error(std::string(kind) +
                             " size mismatch for field ID " +
                             std::to_string(id) + ": local field has " +
                             std::to_string(local.data.size()) +
                             " value(s), buffer has " + std::to_string(size));
  }
}
} // namespace

void CouplingAdapter::register_coupled_field(FieldRegistry& registry,
                                             RegisteredField metadata,
                                             FieldView data) {
  const auto id = registry.register_field(std::move(metadata));

  fields.push_back({
      .id = id,
      .data = data,

  });
}

void CouplingAdapter::configure(const ActiveCouplingFields& plan) {
  active_exports = resolve(plan.export_ids);
  active_imports = resolve(plan.import_ids);
}

std::span<const int> CouplingAdapter::get_local_gids() const {
  return std::span<const int>{local_gids};
}
void CouplingAdapter::set_local_gids(std::vector<int> gids) {
  local_gids = std::move(gids);
}

void CouplingAdapter::export_fields(std::span<ExportBuffer> buffers) const {
  check_buffer_count("export", active_exports.size(), buffers.size());

  for (std::size_t i = 0; i < buffers.size(); ++i) {
    const auto& local = active_exports[i];
    auto& remote = buffers[i];

    bool copy = remote.data.copy_from(local.data);
    if (copy) {
      copies.push_back({.id = remote.id,
                        .bytes = remote.data.size_bytes(),
                        .role = FieldRole::Export});
    }
  }
}

void CouplingAdapter::import_fields(std::span<const ImportBuffer> buffers) {
  check_buffer_count("import", active_imports.size(), buffers.size());

  for (std::size_t i = 0; i < buffers.size(); ++i) {
    auto& local = active_imports[i];
    const auto& remote = buffers[i];

    bool copy = local.data.copy_from(remote.data);
    if (copy) {
      copies.push_back({.id = remote.id,
                        .bytes = remote.data.size_bytes(),
                        .role = FieldRole::Export});
    }
  }
}

/**
 * @brief Helper to resolve the FieldBinding corresponding to a FieldID
 * */
FieldBufferList CouplingAdapter::resolve(std::span<const FieldID> ids) const {

  FieldBufferList result;
  result.reserve(ids.size());

  for (FieldID id : ids) {
    auto it = std::ranges::find_if(
        fields, [id](const ExportBuffer& field) { return field.id == id; });

    if (it == fields.end()) {
      throw std::runtime_error("Coupling plan references field ID " +
                               std::to_string(id) +
                               ", which this component did not register");
    }

    result.push_back(*it);
  }

  return result;
}

} // namespace e3sm::emulator::coupler
