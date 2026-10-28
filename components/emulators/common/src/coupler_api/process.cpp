#include <process.hpp>

#include <algorithm>

namespace e3sm::coupler {

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

void CouplingAdapter::export_fields(std::span<ExportBuffer> buffers) const {
  if (buffers.size() != active_exports_.size()) {
    throw std::runtime_error("Incorrect number of export buffers");
  }

  for (std::size_t i = 0; i < buffers.size(); ++i) {
    const auto& local = active_exports_[i];
    auto& remote = buffers[i];

    if (local.id != remote.id) {
      throw std::runtime_error("Export field ID mismatch");
    }

    if (local.data.size() != remote.data.size()) {
      throw std::runtime_error("Export field size mismatch");
    }

    std::ranges::copy(local.data, remote.data.begin());
  }
}

void CouplingAdapter::import_fields(std::span<const ImportBuffer> buffers) {
  if (buffers.size() != active_imports_.size()) {
    throw std::runtime_error("Incorrect number of import buffers");
  }

  for (std::size_t i = 0; i < buffers.size(); ++i) {
    auto& local = active_imports_[i];
    const auto& remote = buffers[i];

    if (local.id != remote.id) {
      throw std::runtime_error("Import field ID mismatch");
    }

    if (local.data.size() != remote.data.size()) {
      throw std::runtime_error("Import field size mismatch");
    }

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
          "Coupling configuration references unknown field");
    }

    result.push_back(*it);
  }

  return result;
}

} // namespace e3sm::coupler
