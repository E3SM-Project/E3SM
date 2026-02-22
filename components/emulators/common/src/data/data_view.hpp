/**
 * @file data_view.hpp
 * @brief DataView — owning, layout-aware data container for emulator fields.
 *
 * DataView is the central abstraction for bridging coupler data,
 * scorpio-based I/O, and ML/AI tensor backends. It owns a contiguous
 * buffer of doubles organized as (nfields x npoints) or (npoints x nfields)
 * depending on the chosen DataLayout, supports copy-in/copy-out with
 * external buffers (e.g. MCT attribute vectors), and exposes raw pointers
 * for zero-overhead handoff to tensor adapters and IO adapters.
 */

#ifndef E3SM_EMULATOR_DATA_VIEW_HPP
#define E3SM_EMULATOR_DATA_VIEW_HPP

#include "field_spec.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Forward-declare DecompInfo so callers that don't need it
// don't have to include decomp_info.hpp.
namespace emulator {
struct DecompInfo;
}

namespace emulator {

/**
 * @brief Memory layout of the 2D data buffer inside a DataView.
 *
 * The MCT attribute vector uses FIELD_MAJOR (rAttr(nfields, npoints) in
 * Fortran column-major, which is row-major nfields x npoints in C).
 * ML backends typically prefer POINT_MAJOR (each row is a point with
 * nfields features).
 */
enum class DataLayout {
  FIELD_MAJOR, ///< data[field * npoints + point] — contiguous per-field
  POINT_MAJOR  ///< data[point * nfields + field] — contiguous per-point
};

/**
 * @brief Owning, layout-aware data container for emulator fields.
 *
 * ## Lifecycle
 * 1. Construct with a name and npoints (optionally attach DecompInfo later).
 * 2. Register fields via add_field() / add_fields().
 * 3. Call allocate() to reserve the contiguous buffer.
 * 4. Use import_from() to copy data in from an external source (MCT).
 * 5. Access data via data(), field_data(), or operator() for processing.
 * 6. Use export_to() to copy data back to an external destination (MCT).
 *
 * ## Memory layout
 * - FIELD_MAJOR: buffer[field_idx * npoints + point_idx]
 *   — each field is a contiguous slice; good for per-field IO.
 * - POINT_MAJOR: buffer[point_idx * nfields + field_idx]
 *   — each point's features are contiguous; good for ML inference.
 */
class DataView {
public:
  /**
   * @brief Construct a DataView.
   *
   * @param name    Identifier for this view (e.g. "atm_import")
   * @param npoints Number of local grid points
   * @param layout  Memory layout (default: POINT_MAJOR for ML use)
   */
  DataView(std::string name, int npoints,
           DataLayout layout = DataLayout::POINT_MAJOR);

  ~DataView();

  // Non-copyable, movable
  DataView(const DataView &) = delete;
  DataView &operator=(const DataView &) = delete;
  DataView(DataView &&) noexcept;
  DataView &operator=(DataView &&) noexcept;

  // ── Field management ──────────────────────────────────────────────

  /**
   * @brief Register a field. Must be called before allocate().
   * @throws std::runtime_error if already allocated or name collides
   */
  void add_field(const FieldSpec &spec);

  /**
   * @brief Register multiple fields at once.
   */
  void add_fields(const std::vector<FieldSpec> &specs);

  /**
   * @brief Allocate the contiguous data buffer for all registered fields.
   *
   * Must be called exactly once, after all fields are registered.
   * @throws std::runtime_error if no fields or already allocated
   */
  void allocate();

  /** @brief Number of registered fields. */
  int num_fields() const { return static_cast<int>(m_fields.size()); }

  /** @brief Local points on this rank. */
  int num_points() const { return m_npoints; }

  /** @brief Whether allocate() has been called. */
  bool is_allocated() const { return m_allocated; }

  /**
   * @brief Look up a field index by name.
   * @return Field index (0-based)
   * @throws std::runtime_error if name not found
   */
  int field_index(const std::string &name) const;

  /** @brief Get field spec by index. */
  const FieldSpec &field_spec(int idx) const;

  /** @brief Get field spec by name. */
  const FieldSpec &field_spec(const std::string &name) const;

  // ── Data access ───────────────────────────────────────────────────

  /**
   * @brief Raw pointer to the entire contiguous buffer.
   * @throws std::runtime_error if not allocated
   */
  double *data();
  const double *data() const;

  /**
   * @brief Total number of doubles in the buffer (nfields * npoints).
   */
  std::size_t size() const;

  /**
   * @brief Pointer to the start of a single field's data.
   *
   * For FIELD_MAJOR, this is a contiguous slice of npoints doubles.
   * For POINT_MAJOR, the data is strided (stride = nfields).
   *
   * @param field_idx 0-based field index
   * @throws std::runtime_error if not allocated or index out of range
   */
  double *field_data(int field_idx);
  const double *field_data(int field_idx) const;

  /** @brief Overload accepting a field name. */
  double *field_data(const std::string &name);
  const double *field_data(const std::string &name) const;

  /**
   * @brief Element access.
   *
   * @param field_idx 0-based field index
   * @param point_idx 0-based point index
   */
  double &operator()(int field_idx, int point_idx);
  double operator()(int field_idx, int point_idx) const;

  /**
   * @brief Logical shape of the 2D view.
   *
   * @return {dim0, dim1} where:
   *   FIELD_MAJOR → {nfields, npoints}
   *   POINT_MAJOR → {npoints, nfields}
   */
  std::pair<int, int> shape() const;

  // ── Import / Export (copy-in / copy-out) ──────────────────────────

  /**
   * @brief Copy data from an external buffer into this DataView.
   *
   * Handles layout transposition if src_layout differs from this
   * view's layout.
   *
   * @param src         Pointer to the source buffer
   * @param src_nfields Number of fields in the source buffer
   * @param src_npoints Number of points in the source buffer
   * @param src_layout  Memory layout of the source buffer
   * @param field_map   Mapping: {src_field_idx, dst_field_idx} pairs.
   *                    If empty, copies min(src_nfields, nfields)
   *                    fields in order.
   */
  void import_from(const double *src, int src_nfields, int src_npoints,
                   DataLayout src_layout,
                   const std::vector<std::pair<int, int>> &field_map = {});

  /**
   * @brief Copy data from this DataView into an external buffer.
   *
   * Handles layout transposition if dst_layout differs from this
   * view's layout.
   *
   * @param dst         Pointer to the destination buffer
   * @param dst_nfields Number of fields in the destination buffer
   * @param dst_npoints Number of points in the destination buffer
   * @param dst_layout  Memory layout of the destination buffer
   * @param field_map   Mapping: {src_field_idx (this view), dst_field_idx}
   *                    pairs. If empty, copies min(nfields, dst_nfields)
   *                    fields in order.
   */
  void export_to(double *dst, int dst_nfields, int dst_npoints,
                 DataLayout dst_layout,
                 const std::vector<std::pair<int, int>> &field_map = {}) const;

  // ── Utilities ─────────────────────────────────────────────────────

  /** @brief Zero the entire buffer. */
  void zero();

  /** @brief Fill the entire buffer with a constant value. */
  void fill(double value);

  /** @brief Copy all field data from another DataView of the same shape. */
  void copy_from(const DataView &other);

  // ── Decomposition (optional) ──────────────────────────────────────

  /**
   * @brief Attach domain decomposition metadata for parallel IO.
   *
   * Not needed for standalone ML use. Required before using IOAdapter.
   */
  void set_decomp(DecompInfo decomp);

  /** @brief Whether decomposition info has been attached. */
  bool has_decomp() const;

  /**
   * @brief Get decomposition info.
   * @throws std::runtime_error if no decomp attached
   */
  const DecompInfo &decomp() const;

  // ── Accessors ─────────────────────────────────────────────────────

  const std::string &name() const { return m_name; }
  DataLayout layout() const { return m_layout; }

private:
  /** @brief Compute buffer offset for (field_idx, point_idx). */
  std::size_t offset(int field_idx, int point_idx) const;

  std::string m_name;
  int m_npoints;
  DataLayout m_layout;

  std::vector<FieldSpec> m_fields;
  std::unordered_map<std::string, int> m_field_index; ///< name → index

  std::vector<double> m_data;
  bool m_allocated = false;

  std::unique_ptr<DecompInfo> m_decomp; ///< Optional decomposition
};

} // namespace emulator

#endif // E3SM_EMULATOR_DATA_VIEW_HPP
