/**
 * @file data_view.cpp
 * @brief Implementation of the DataView class.
 */

#include "data_view.hpp"
#include "decomp_info.hpp"

#include <algorithm>
#include <cstring>
#include <stdexcept>

namespace emulator {

// ── Construction ──────────────────────────────────────────────────────

DataView::DataView(std::string name, int npoints, DataLayout layout)
    : m_name(std::move(name)), m_npoints(npoints), m_layout(layout) {
  if (m_npoints < 0) {
    throw std::runtime_error("DataView: npoints must be non-negative");
  }
}

DataView::~DataView() = default;
DataView::DataView(DataView &&) noexcept = default;
DataView &DataView::operator=(DataView &&) noexcept = default;

// ── Field management ──────────────────────────────────────────────────

void DataView::add_field(const FieldSpec &spec) {
  if (m_allocated) {
    throw std::runtime_error(
        "DataView::add_field(): cannot add fields after allocation");
  }
  if (spec.name.empty()) {
    throw std::runtime_error("DataView::add_field(): field name is empty");
  }
  if (m_field_index.count(spec.name)) {
    throw std::runtime_error("DataView::add_field(): duplicate field name '" +
                             spec.name + "'");
  }
  int idx = static_cast<int>(m_fields.size());
  m_fields.push_back(spec);
  m_field_index[spec.name] = idx;
}

void DataView::add_fields(const std::vector<FieldSpec> &specs) {
  for (const auto &s : specs) {
    add_field(s);
  }
}

void DataView::allocate() {
  if (m_allocated) {
    throw std::runtime_error("DataView::allocate(): already allocated");
  }
  if (m_fields.empty()) {
    throw std::runtime_error("DataView::allocate(): no fields registered");
  }
  std::size_t total = static_cast<std::size_t>(num_fields()) *
                      static_cast<std::size_t>(num_points());
  m_data.assign(total, 0.0);
  m_allocated = true;
}

int DataView::field_index(const std::string &name) const {
  auto it = m_field_index.find(name);
  if (it == m_field_index.end()) {
    throw std::runtime_error("DataView::field_index(): unknown field '" + name +
                             "'");
  }
  return it->second;
}

const FieldSpec &DataView::field_spec(int idx) const {
  if (idx < 0 || idx >= num_fields()) {
    throw std::runtime_error("DataView::field_spec(): index out of range");
  }
  return m_fields[static_cast<std::size_t>(idx)];
}

const FieldSpec &DataView::field_spec(const std::string &name) const {
  return field_spec(field_index(name));
}

// ── Data access ───────────────────────────────────────────────────────

double *DataView::data() {
  if (!m_allocated) {
    throw std::runtime_error("DataView::data(): not allocated");
  }
  return m_data.data();
}

const double *DataView::data() const {
  if (!m_allocated) {
    throw std::runtime_error("DataView::data(): not allocated");
  }
  return m_data.data();
}

std::size_t DataView::size() const { return m_data.size(); }

double *DataView::field_data(int field_idx) {
  if (!m_allocated) {
    throw std::runtime_error("DataView::field_data(): not allocated");
  }
  if (field_idx < 0 || field_idx >= num_fields()) {
    throw std::runtime_error(
        "DataView::field_data(): field index out of range");
  }
  // Return pointer to first element of this field
  return &m_data[offset(field_idx, 0)];
}

const double *DataView::field_data(int field_idx) const {
  if (!m_allocated) {
    throw std::runtime_error("DataView::field_data(): not allocated");
  }
  if (field_idx < 0 || field_idx >= num_fields()) {
    throw std::runtime_error(
        "DataView::field_data(): field index out of range");
  }
  return &m_data[offset(field_idx, 0)];
}

double *DataView::field_data(const std::string &name) {
  return field_data(field_index(name));
}

const double *DataView::field_data(const std::string &name) const {
  return field_data(field_index(name));
}

double &DataView::operator()(int field_idx, int point_idx) {
  return m_data[offset(field_idx, point_idx)];
}

double DataView::operator()(int field_idx, int point_idx) const {
  return m_data[offset(field_idx, point_idx)];
}

std::pair<int, int> DataView::shape() const {
  if (m_layout == DataLayout::FIELD_MAJOR) {
    return {num_fields(), num_points()};
  }
  return {num_points(), num_fields()};
}

// ── Import / Export ───────────────────────────────────────────────────

void DataView::import_from(const double *src, int src_nfields, int src_npoints,
                           DataLayout src_layout,
                           const std::vector<std::pair<int, int>> &field_map) {
  if (!m_allocated) {
    throw std::runtime_error("DataView::import_from(): not allocated");
  }

  // Build the mapping: vector of (src_field, dst_field) pairs
  std::vector<std::pair<int, int>> mapping = field_map;
  if (mapping.empty()) {
    int n = std::min(src_nfields, num_fields());
    mapping.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      mapping.emplace_back(i, i);
    }
  }

  int npts = std::min(src_npoints, num_points());

  for (const auto &[sf, df] : mapping) {
    for (int p = 0; p < npts; ++p) {
      // Source offset
      std::size_t src_off;
      if (src_layout == DataLayout::FIELD_MAJOR) {
        src_off = static_cast<std::size_t>(sf) *
                      static_cast<std::size_t>(src_npoints) +
                  static_cast<std::size_t>(p);
      } else {
        src_off = static_cast<std::size_t>(p) *
                      static_cast<std::size_t>(src_nfields) +
                  static_cast<std::size_t>(sf);
      }
      // Destination offset (in our layout)
      m_data[offset(df, p)] = src[src_off];
    }
  }
}

void DataView::export_to(
    double *dst, int dst_nfields, int dst_npoints, DataLayout dst_layout,
    const std::vector<std::pair<int, int>> &field_map) const {
  if (!m_allocated) {
    throw std::runtime_error("DataView::export_to(): not allocated");
  }

  // Build the mapping: vector of (src_field (this), dst_field) pairs
  std::vector<std::pair<int, int>> mapping = field_map;
  if (mapping.empty()) {
    int n = std::min(num_fields(), dst_nfields);
    mapping.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      mapping.emplace_back(i, i);
    }
  }

  int npts = std::min(num_points(), dst_npoints);

  for (const auto &[sf, df] : mapping) {
    for (int p = 0; p < npts; ++p) {
      // Destination offset
      std::size_t dst_off;
      if (dst_layout == DataLayout::FIELD_MAJOR) {
        dst_off = static_cast<std::size_t>(df) *
                      static_cast<std::size_t>(dst_npoints) +
                  static_cast<std::size_t>(p);
      } else {
        dst_off = static_cast<std::size_t>(p) *
                      static_cast<std::size_t>(dst_nfields) +
                  static_cast<std::size_t>(df);
      }
      // Source offset (from our layout)
      dst[dst_off] = m_data[offset(sf, p)];
    }
  }
}

// ── Utilities ─────────────────────────────────────────────────────────

void DataView::zero() { std::fill(m_data.begin(), m_data.end(), 0.0); }

void DataView::fill(double value) {
  std::fill(m_data.begin(), m_data.end(), value);
}

void DataView::copy_from(const DataView &other) {
  if (!m_allocated || !other.m_allocated) {
    throw std::runtime_error(
        "DataView::copy_from(): both views must be allocated");
  }
  if (m_data.size() != other.m_data.size()) {
    throw std::runtime_error("DataView::copy_from(): size mismatch");
  }
  if (m_layout == other.m_layout) {
    std::memcpy(m_data.data(), other.m_data.data(),
                m_data.size() * sizeof(double));
  } else {
    // Transpose copy
    for (int f = 0; f < num_fields(); ++f) {
      for (int p = 0; p < num_points(); ++p) {
        m_data[offset(f, p)] = other.m_data[other.offset(f, p)];
      }
    }
  }
}

// ── Decomposition (optional) ──────────────────────────────────────────

void DataView::set_decomp(DecompInfo decomp) {
  decomp.validate();
  if (decomp.npoints_local != m_npoints) {
    throw std::runtime_error("DataView::set_decomp(): decomp npoints_local (" +
                             std::to_string(decomp.npoints_local) +
                             ") != view npoints (" + std::to_string(m_npoints) +
                             ")");
  }
  m_decomp = std::make_unique<DecompInfo>(std::move(decomp));
}

bool DataView::has_decomp() const { return m_decomp != nullptr; }

const DecompInfo &DataView::decomp() const {
  if (!m_decomp) {
    throw std::runtime_error(
        "DataView::decomp(): no decomposition info attached");
  }
  return *m_decomp;
}

// ── Private helpers ───────────────────────────────────────────────────

std::size_t DataView::offset(int field_idx, int point_idx) const {
  if (m_layout == DataLayout::FIELD_MAJOR) {
    return static_cast<std::size_t>(field_idx) *
               static_cast<std::size_t>(num_points()) +
           static_cast<std::size_t>(point_idx);
  }
  // POINT_MAJOR
  return static_cast<std::size_t>(point_idx) *
             static_cast<std::size_t>(num_fields()) +
         static_cast<std::size_t>(field_idx);
}

} // namespace emulator
