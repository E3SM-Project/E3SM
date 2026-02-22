/**
 * @file data_source.hpp
 * @brief Abstract interface for populating and draining DataViews.
 *
 * DataSource is the input/output abstraction that allows DataView to
 * receive data from arbitrary external sources (MCT coupler buffers,
 * NetCDF files, tensor outputs, test fixtures, etc.) without knowing
 * the specifics of each source.
 */

#ifndef E3SM_EMULATOR_DATA_SOURCE_HPP
#define E3SM_EMULATOR_DATA_SOURCE_HPP

#include "data_view.hpp"
#include "field_spec.hpp"

#include <vector>

namespace emulator {

/**
 * @brief Abstract interface for anything that can fill/drain a DataView.
 *
 * ## Concrete implementations
 * - MctDataSource   : wraps MCT attribute vector (double* + layout)
 * - (future) FileDataSource   : wraps a NetCDF/HDF5 file path
 * - (future) TensorDataSource : wraps an existing tensor buffer
 *
 * ## Usage
 * ```cpp
 * auto source = std::make_unique<MctDataSource>(...);
 * source->import_data(view);   // external → DataView
 * // ... process ...
 * source->export_data(view);   // DataView → external
 * ```
 */
class DataSource {
public:
  virtual ~DataSource() = default;

  /**
   * @brief Report which fields this source provides.
   *
   * Used by DataPipeline to auto-register fields on a DataView.
   */
  virtual std::vector<FieldSpec> available_fields() const = 0;

  /**
   * @brief Copy data from the external source into the DataView.
   *
   * The DataView must already be allocated with compatible fields.
   * @param view Destination DataView (modified in-place)
   */
  virtual void import_data(DataView &view) = 0;

  /**
   * @brief Copy data from the DataView back to the external source.
   *
   * Default implementation is a no-op — not all sources are writable.
   * @param view Source DataView (read-only)
   */
  virtual void export_data(const DataView &view) { (void)view; }
};

} // namespace emulator

#endif // E3SM_EMULATOR_DATA_SOURCE_HPP
