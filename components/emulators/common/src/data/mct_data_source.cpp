/**
 * @file mct_data_source.cpp
 * @brief Implementation of MctDataSource.
 */

#include "mct_data_source.hpp"

#include <stdexcept>

namespace emulator {

MctDataSource::MctDataSource(const double *import_buf, int nfields, int npoints,
                             DataLayout layout, std::vector<FieldSpec> fields,
                             std::vector<std::pair<int, int>> field_map)
    : m_import_buf(import_buf), m_nfields(nfields), m_npoints(npoints),
      m_layout(layout), m_fields(std::move(fields)),
      m_field_map(std::move(field_map)) {}

void MctDataSource::set_export_buffer(double *buf) { m_export_buf = buf; }

void MctDataSource::update_import_buffer(const double *buf) {
  m_import_buf = buf;
}

std::vector<FieldSpec> MctDataSource::available_fields() const {
  return m_fields;
}

void MctDataSource::import_data(DataView &view) {
  if (!m_import_buf) {
    throw std::runtime_error(
        "MctDataSource::import_data(): import buffer is null");
  }
  view.import_from(m_import_buf, m_nfields, m_npoints, m_layout, m_field_map);
}

void MctDataSource::export_data(const DataView &view) {
  if (!m_export_buf) {
    return; // No export buffer set — no-op (read-only source)
  }
  view.export_to(m_export_buf, m_nfields, m_npoints, m_layout, m_field_map);
}

} // namespace emulator
