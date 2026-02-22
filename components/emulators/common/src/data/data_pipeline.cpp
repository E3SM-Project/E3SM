/**
 * @file data_pipeline.cpp
 * @brief Implementation of DataPipeline.
 */

#include "data_pipeline.hpp"

#include <stdexcept>
#include <utility>

namespace emulator {

// ── Construction ──────────────────────────────────────────────────────

DataPipeline::DataPipeline(std::string name, int npoints, DataLayout layout)
    : m_name(std::move(name)), m_npoints(npoints), m_layout(layout) {}

DataPipeline::~DataPipeline() = default;
DataPipeline::DataPipeline(DataPipeline &&) noexcept = default;
DataPipeline &DataPipeline::operator=(DataPipeline &&) noexcept = default;

// ── Assembly ──────────────────────────────────────────────────────────

void DataPipeline::set_import_source(std::unique_ptr<DataSource> source) {
  m_import_source = std::move(source);
}

void DataPipeline::set_export_source(std::unique_ptr<DataSource> source) {
  m_export_source = std::move(source);
}

void DataPipeline::set_tensor_adapter(std::unique_ptr<TensorAdapter> adapter) {
  m_tensor_adapter = std::move(adapter);
}

void DataPipeline::set_io_adapter(std::unique_ptr<IOAdapter> io) {
  m_io_adapter = std::move(io);
}

void DataPipeline::add_import_field(const FieldSpec &spec) {
  if (m_initialized) {
    throw std::runtime_error(
        "DataPipeline::add_import_field(): already initialized");
  }
  m_import_fields.push_back(spec);
}

void DataPipeline::add_export_field(const FieldSpec &spec) {
  if (m_initialized) {
    throw std::runtime_error(
        "DataPipeline::add_export_field(): already initialized");
  }
  m_export_fields.push_back(spec);
}

void DataPipeline::auto_register_fields() {
  if (m_initialized) {
    throw std::runtime_error(
        "DataPipeline::auto_register_fields(): already initialized");
  }
  if (m_import_source) {
    auto fields = m_import_source->available_fields();
    for (auto &f : fields) {
      m_import_fields.push_back(std::move(f));
    }
  }
  if (m_export_source) {
    auto fields = m_export_source->available_fields();
    for (auto &f : fields) {
      m_export_fields.push_back(std::move(f));
    }
  }
}

void DataPipeline::initialize() {
  if (m_initialized) {
    throw std::runtime_error("DataPipeline::initialize(): already initialized");
  }

  // Build import view
  if (!m_import_fields.empty()) {
    m_import_view =
        std::make_unique<DataView>(m_name + "_import", m_npoints, m_layout);
    m_import_view->add_fields(m_import_fields);
    m_import_view->allocate();
  }

  // Build export view
  if (!m_export_fields.empty()) {
    m_export_view =
        std::make_unique<DataView>(m_name + "_export", m_npoints, m_layout);
    m_export_view->add_fields(m_export_fields);
    m_export_view->allocate();
  }

  if (!m_import_view && !m_export_view) {
    throw std::runtime_error(
        "DataPipeline::initialize(): no fields registered on either view");
  }

  m_initialized = true;
}

// ── Runtime ───────────────────────────────────────────────────────────

void DataPipeline::ingest() {
  if (!m_import_source) {
    throw std::runtime_error("DataPipeline::ingest(): no import source set");
  }
  if (!m_import_view) {
    throw std::runtime_error(
        "DataPipeline::ingest(): no import view (no import fields)");
  }
  m_import_source->import_data(*m_import_view);
}

void DataPipeline::egest() {
  if (!m_export_source) {
    throw std::runtime_error("DataPipeline::egest(): no export source set");
  }
  if (!m_export_view) {
    throw std::runtime_error(
        "DataPipeline::egest(): no export view (no export fields)");
  }
  m_export_source->export_data(*m_export_view);
}

// ── Tensor access ─────────────────────────────────────────────────────

const double *DataPipeline::tensor_data() const {
  if (!m_tensor_adapter) {
    throw std::runtime_error(
        "DataPipeline::tensor_data(): no tensor adapter set");
  }
  if (!m_import_view) {
    throw std::runtime_error("DataPipeline::tensor_data(): no import view");
  }
  return m_tensor_adapter->data_ptr(*m_import_view);
}

std::vector<int64_t> DataPipeline::tensor_shape() const {
  if (!m_tensor_adapter) {
    throw std::runtime_error(
        "DataPipeline::tensor_shape(): no tensor adapter set");
  }
  if (!m_import_view) {
    throw std::runtime_error("DataPipeline::tensor_shape(): no import view");
  }
  return m_tensor_adapter->tensor_shape(*m_import_view);
}

void DataPipeline::write_tensor(const void *ml_output) {
  if (!m_tensor_adapter) {
    throw std::runtime_error(
        "DataPipeline::write_tensor(): no tensor adapter set");
  }
  if (!m_export_view) {
    throw std::runtime_error("DataPipeline::write_tensor(): no export view");
  }
  m_tensor_adapter->from_tensor(ml_output, *m_export_view);
}

// ── IO access ─────────────────────────────────────────────────────────

void DataPipeline::write_io(const std::string &filename) {
  if (!m_io_adapter) {
    throw std::runtime_error("DataPipeline::write_io(): no IO adapter set");
  }
  if (!m_import_view) {
    throw std::runtime_error("DataPipeline::write_io(): no import view");
  }
  m_io_adapter->write(*m_import_view, filename);
}

void DataPipeline::read_io(const std::string &filename) {
  if (!m_io_adapter) {
    throw std::runtime_error("DataPipeline::read_io(): no IO adapter set");
  }
  if (!m_import_view) {
    throw std::runtime_error("DataPipeline::read_io(): no import view");
  }
  m_io_adapter->read(*m_import_view, filename);
}

// ── Direct access ─────────────────────────────────────────────────────

DataView &DataPipeline::import_view() {
  if (!m_import_view) {
    throw std::runtime_error("DataPipeline::import_view(): no import view");
  }
  return *m_import_view;
}

const DataView &DataPipeline::import_view() const {
  if (!m_import_view) {
    throw std::runtime_error("DataPipeline::import_view(): no import view");
  }
  return *m_import_view;
}

DataView &DataPipeline::export_view() {
  if (!m_export_view) {
    throw std::runtime_error("DataPipeline::export_view(): no export view");
  }
  return *m_export_view;
}

const DataView &DataPipeline::export_view() const {
  if (!m_export_view) {
    throw std::runtime_error("DataPipeline::export_view(): no export view");
  }
  return *m_export_view;
}

} // namespace emulator
