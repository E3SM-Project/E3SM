/**
 * @file data_view.hpp
 * @brief Non-owning view abstraction for zero-copy data access.
 *
 * Provides a lightweight, non-owning view over contiguous data that can
 * wrap raw pointers, std::vector, or (in the future) Kokkos views without
 * copying data.
 */

#ifndef EMULATOR_DATA_VIEW_HPP
#define EMULATOR_DATA_VIEW_HPP

#include <cstddef>
#include <vector>

namespace emulator {

/**
 * @brief Memory layout for multi-dimensional data.
 */
enum class DataLayout {
  ROW_MAJOR,    ///< C-style: last dimension varies fastest
  COLUMN_MAJOR, ///< Fortran-style: first dimension varies fastest
  UNKNOWN       ///< Layout not specified
};

/**
 * @brief Non-owning view over contiguous data.
 *
 * Provides zero-copy access to data stored elsewhere. Does not manage
 * lifetime - the caller must ensure the underlying data remains valid.
 *
 * ## Usage
 *
 * ```cpp
 * std::vector<double> data = {1.0, 2.0, 3.0};
 * DataView<double> view(data);
 * // or from raw pointer:
 * DataView<double> view(ptr, size);
 * ```
 *
 * ## Future Kokkos Integration
 *
 * When implementing the LAPIS backend, DataView can be constructed
 * directly from Kokkos::View without copying, enabling zero-copy
 * data paths between E3SM components and ML inference.
 *
 * @tparam T Element type (typically double or float)
 */
template <typename T> class DataView {
public:
  /// Default constructor: empty view
  DataView() : m_data(nullptr), m_size(0), m_layout(DataLayout::UNKNOWN) {}

  /// Construct from raw pointer and size
  DataView(T *data, std::size_t size, DataLayout layout = DataLayout::ROW_MAJOR)
      : m_data(data), m_size(size), m_layout(layout) {}

  /// Construct from raw pointer and size (const version)
  DataView(const T *data, std::size_t size,
           DataLayout layout = DataLayout::ROW_MAJOR)
      : m_data(const_cast<T *>(data)), m_size(size), m_layout(layout) {}

  /// Construct from std::vector (non-owning)
  explicit DataView(std::vector<T> &vec,
                    DataLayout layout = DataLayout::ROW_MAJOR)
      : m_data(vec.data()), m_size(vec.size()), m_layout(layout) {}

  /// Construct from const std::vector (non-owning)
  explicit DataView(const std::vector<T> &vec,
                    DataLayout layout = DataLayout::ROW_MAJOR)
      : m_data(const_cast<T *>(vec.data())), m_size(vec.size()),
        m_layout(layout) {}

  /// Get raw pointer to data
  T *data() { return m_data; }
  const T *data() const { return m_data; }

  /// Get number of elements
  std::size_t size() const { return m_size; }

  /// Check if view is empty
  bool empty() const { return m_size == 0 || m_data == nullptr; }

  /// Get memory layout
  DataLayout layout() const { return m_layout; }

  /// Element access (bounds checking in debug builds)
  T &operator[](std::size_t idx) { return m_data[idx]; }
  const T &operator[](std::size_t idx) const { return m_data[idx]; }

  /// Iterator support
  T *begin() { return m_data; }
  T *end() { return m_data + m_size; }
  const T *begin() const { return m_data; }
  const T *end() const { return m_data + m_size; }

private:
  T *m_data;
  std::size_t m_size;
  DataLayout m_layout;
};

/**
 * @brief Convenience type aliases.
 */
using DoubleView = DataView<double>;
using FloatView = DataView<float>;

} // namespace emulator

#endif // EMULATOR_DATA_VIEW_HPP
