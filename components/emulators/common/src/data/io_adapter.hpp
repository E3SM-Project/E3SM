/**
 * @file io_adapter.hpp
 * @brief Abstract interface for DataView ↔ parallel I/O.
 *
 * IOAdapter encapsulates the mapping from DataView fields and decomposition
 * metadata to a parallel I/O library (e.g. scorpio/PIO). The abstract
 * interface keeps emulator_common free of scorpio compile-time dependencies;
 * a concrete ScorpioIOAdapter is provided by a higher layer or a specific
 * emulator component that links against scorpio.
 */

#ifndef E3SM_EMULATOR_DATA_IO_ADAPTER_HPP
#define E3SM_EMULATOR_DATA_IO_ADAPTER_HPP

#include "data_view.hpp"

#include <string>

namespace emulator {

/**
 * @brief Abstract interface for reading/writing DataView data to files.
 *
 * ## Concrete implementation guide (e.g. ScorpioIOAdapter)
 *
 * A concrete adapter would:
 * 1. In init_io(): call pio_initdecomp() using view.decomp().dof
 *    and view.decomp().global_dims.
 * 2. In define_variables(): call pio_def_var() for each field using
 *    view.field_spec(i).name, .units, .long_name.
 * 3. In write(): call pio_write_darray() per field using
 *    view.field_data(i) as the source buffer.
 * 4. In read(): call pio_read_darray() per field into
 *    view.field_data(i).
 *
 * ## Lifecycle
 * 1. Construct with backend-specific config (PIO system, IO type, etc.)
 * 2. Call init_io(view) to set up decomposition.
 * 3. Call write(view, filename) or read(view, filename) as needed.
 * 4. Destructor or explicit close cleans up IO resources.
 */
class IOAdapter {
public:
  virtual ~IOAdapter() = default;

  /**
   * @brief Initialize I/O decomposition from the DataView's metadata.
   *
   * Sets up whatever internal bookkeeping the backend needs (e.g.
   * pio_initdecomp, creating IO descriptors).
   *
   * @param view DataView whose decomp() and field metadata are used.
   */
  virtual void init_io(const DataView &view) = 0;

  /**
   * @brief Define netCDF variables for all fields in the DataView.
   *
   * Called between file creation and the first write. Maps each
   * FieldSpec to a netCDF variable with appropriate dimensions,
   * units, and long_name attributes.
   *
   * @param view DataView whose field specs are used.
   */
  virtual void define_variables(const DataView &view) = 0;

  /**
   * @brief Write all fields from the DataView to a file.
   *
   * @param view     Source DataView (field data is read)
   * @param filename Path to the output file
   */
  virtual void write(const DataView &view,
                     const std::string &filename) = 0;

  /**
   * @brief Read field data from a file into the DataView.
   *
   * @param view     Destination DataView (field data is overwritten)
   * @param filename Path to the input file
   */
  virtual void read(DataView &view,
                    const std::string &filename) = 0;

  /**
   * @brief Close any open file handles and release IO resources.
   */
  virtual void close() = 0;
};

} // namespace emulator

#endif // E3SM_EMULATOR_DATA_IO_ADAPTER_HPP
