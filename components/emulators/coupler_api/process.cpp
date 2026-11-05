#include <process.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>

namespace e3sm::coupler {

namespace {
void check_buffer_count(const char* kind, std::size_t expected,
                        std::size_t actual) {
  if (expected != actual) {
    throw std::runtime_error(std::string("Incorrect number of ") + kind +
                             " buffers: expected " + std::to_string(expected) +
                             ", got " + std::to_string(actual));
  }
}

void check_buffer_match(const char* kind, const FieldBinding& local, FieldID id,
                        std::size_t size) {
  if (local.id != id) {
    throw std::runtime_error(std::string(kind) + " field ID mismatch: expected " +
                             std::to_string(local.id) + ", got " +
                             std::to_string(id));
  }
  if (local.data.size() != size) {
    throw std::runtime_error(std::string(kind) + " size mismatch for field ID " +
                             std::to_string(id) + ": local field has " +
                             std::to_string(local.data.size()) +
                             " value(s), buffer has " + std::to_string(size));
  }
}
} // namespace

void CouplingAdapter::register_coupled_field(FieldRegistry& registry,
                                             RegisteredField metadata,
                                             std::span<double> data) {
  const auto id = registry.register_field(std::move(metadata));

  fields_.push_back({
      .id = id,
      .data = data,
  });
}

void CouplingAdapter::configure(const ActiveCouplingFields& plan) {
  active_exports_ = resolve(plan.export_ids);
  active_imports_ = resolve(plan.import_ids);
}

std::span<const int> CouplingAdapter::get_local_gids() const {
  return std::span<const int>{local_gids_};
}
void CouplingAdapter::set_local_gids(std::vector<int> gids) {
  local_gids_ = std::move(gids);
}

void CouplingAdapter::export_fields(std::span<ExportBuffer> buffers) const {
  check_buffer_count("export", active_exports_.size(), buffers.size());

  for (std::size_t i = 0; i < buffers.size(); ++i) {
    const auto& local = active_exports_[i];
    auto& remote = buffers[i];

    if (local.data.data() == remote.data.data()) {
      continue; // Route shares this field's storage
    }

    check_buffer_match("Export", local, remote.id, remote.data.size());
    std::ranges::copy(local.data, remote.data.begin());
  }
}

void CouplingAdapter::import_fields(std::span<const ImportBuffer> buffers) {
  check_buffer_count("import", active_imports_.size(), buffers.size());

  for (std::size_t i = 0; i < buffers.size(); ++i) {
    auto& local = active_imports_[i];
    const auto& remote = buffers[i];

    if (local.data.data() == remote.data.data()) {
      continue; // Route shares this field's storage
    }

    check_buffer_match("Import", local, remote.id, remote.data.size());
    std::ranges::copy(remote.data, local.data.begin());
  }
}

/**
 * @brief Helper to resolve the FieldBinding corresponding to a FieldID
 * */
std::vector<FieldBinding>
CouplingAdapter::resolve(std::span<const FieldID> ids) const {
  std::vector<FieldBinding> result;
  result.reserve(ids.size());

  for (FieldID id : ids) {
    auto it = std::ranges::find_if(
        fields_, [id](const FieldBinding& field) { return field.id == id; });

    if (it == fields_.end()) {
      throw std::runtime_error(
          "Coupling plan references field ID " + std::to_string(id) +
          ", which this component did not register");
    }

    result.push_back(*it);
  }

  return result;
}

} // namespace e3sm::coupler
