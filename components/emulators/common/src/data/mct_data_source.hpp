/**
 * @file mct_data_source.hpp
 * @brief Concrete DataSource wrapping MCT attribute vector buffers.
 *
 * MctDataSource is the bridge between E3SM's MCT coupler layer and the
 * emulator data pipeline. It wraps a flat double* buffer (MCT rAttr)
 * and delegates to DataView::import_from() / export_to().
 */

#ifndef E3SM_EMULATOR_MCT_DATA_SOURCE_HPP
#define E3SM_EMULATOR_MCT_DATA_SOURCE_HPP

#include "data_source.hpp"

#include <utility>
#include <vector>

namespace emulator {

/**
 * @brief DataSource that wraps an MCT attribute vector (double*).
 *
 * ## Usage
 * ```cpp
 * // During init: configure source with MCT buffer info
 * MctDataSource src(rAttr_ptr, nfields, npoints,
 *                   DataLayout::FIELD_MAJOR, field_specs);
 *
 * // During run: import MCT data into DataView
 * src.import_data(import_view);
 *
 * // After processing: export results back to MCT
 * src.set_export_buffer(export_rAttr_ptr);
 * src.export_data(export_view);
 * ```
 */
class MctDataSource : public DataSource {
public:
  /**
   * @brief Construct an MctDataSource.
   *
   * @param import_buf  Pointer to the MCT import buffer (read-only)
   * @param nfields     Number of fields in the MCT buffer
   * @param npoints     Number of points in the MCT buffer
   * @param layout      Memory layout of the MCT buffer
   * @param fields      Field specifications for the MCT fields
   * @param field_map   Optional field mapping: {src_idx, dst_idx} pairs.
   *                    If empty, maps fields 1:1 in order.
   */
  MctDataSource(const double *import_buf, int nfields, int npoints,
                DataLayout layout, std::vector<FieldSpec> fields,
                std::vector<std::pair<int, int>> field_map = {});

  /**
   * @brief Set the export buffer for writing data back to MCT.
   *
   * Must be called before export_data() if bidirectional use is needed.
   * @param buf Pointer to the MCT export buffer (writable)
   */
  void set_export_buffer(double *buf);

  /**
   * @brief Update the import buffer pointer (e.g. for next timestep).
   * @param buf New pointer to the MCT import buffer
   */
  void update_import_buffer(const double *buf);

  // ── DataSource interface ────────────────────────────────────────

  std::vector<FieldSpec> available_fields() const override;
  void import_data(DataView &view) override;
  void export_data(const DataView &view) override;

private:
  const double *m_import_buf;
  double *m_export_buf = nullptr;
  int m_nfields;
  int m_npoints;
  DataLayout m_layout;
  std::vector<FieldSpec> m_fields;
  std::vector<std::pair<int, int>> m_field_map;
};

} // namespace emulator

#endif // E3SM_EMULATOR_MCT_DATA_SOURCE_HPP
