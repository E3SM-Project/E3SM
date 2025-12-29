/**
 * @file emulator_io.hpp
 * @brief I/O wrapper using SCORPIO/PIO C interface.
 *
 * Provides a C++ interface to the SCORPIO (PIO) library for parallel
 * NetCDF I/O operations. Similar to EAMxx's eamxx_scorpio_interface.
 */

#ifndef EMULATOR_IO_HPP
#define EMULATOR_IO_HPP

#include <mpi.h>
#include <string>
#include <vector>

// PIO C interface from SCORPIO
extern "C" {
#include <pio.h>
}

namespace emulator {

/**
 * @brief Static class providing parallel I/O using SCORPIO/PIO.
 *
 * Wraps the SCORPIO (Scalable Parallel I/O) C library to provide
 * convenient C++ methods for reading and writing NetCDF files in
 * parallel. All methods are static as only one PIO system is active.
 *
 * ## Initialization
 * Call `initialize()` once per component with the MPI communicator.
 * Call `finalize()` before MPI_Finalize.
 *
 * ## File Operations
 * - `open_file()` / `create_file()` - Open or create NetCDF files
 * - `close_file()` - Close an open file
 * - `sync_file()` - Flush writes to disk
 *
 * ## Variable I/O
 * - `read_var_1d()`, `read_var_2d()` - Read array data
 * - `write_var_1d()`, `write_var_2d()` - Write array data
 * - `has_var()` - Check if variable exists
 *
 * @note Currently uses serial NetCDF (PIO_IOTYPE_NETCDF) for simplicity.
 *       For large-scale runs, consider PIO_IOTYPE_PNETCDF.
 */
class EmulatorIO {
public:
  // =========================================================================
  // Initialization
  // =========================================================================

  /**
   * @brief Initialize the PIO subsystem.
   *
   * Must be called once before any file operations. Typically called
   * from EmulatorComp::create_instance().
   *
   * @param comm MPI communicator for parallel I/O
   * @param comp_name Component name (for logging)
   */
  static void initialize(MPI_Comm comm, const std::string &comp_name);

  /**
   * @brief Finalize the PIO subsystem.
   *
   * Should be called before MPI_Finalize to release resources.
   */
  static void finalize();

  /**
   * @brief Check if PIO is initialized.
   * @return true if initialize() has been called successfully
   */
  static bool is_initialized() { return s_initialized; }

  // =========================================================================
  // File Operations
  // =========================================================================

  /**
   * @brief Open an existing NetCDF file for reading.
   * @param filename Path to the file
   * @return NetCDF file ID (ncid), or -1 on failure
   */
  static int open_file(const std::string &filename);

  /**
   * @brief Create a new NetCDF file for writing.
   * @param filename Path to the new file
   * @return NetCDF file ID (ncid), or -1 on failure
   */
  static int create_file(const std::string &filename);

  /**
   * @brief Close an open file.
   * @param ncid File ID from open_file() or create_file()
   */
  static void close_file(int ncid);

  /**
   * @brief Synchronize (flush) file contents to disk.
   * @param ncid File ID
   */
  static void sync_file(int ncid);

  // =========================================================================
  // Variable I/O
  // =========================================================================

  /**
   * @brief Read a 1D double array from a NetCDF variable.
   * @param ncid File ID
   * @param varname Variable name
   * @param data Output buffer (must be pre-allocated)
   * @param size Number of elements to read
   * @return true on success, false if variable not found or read failed
   */
  static bool read_var_1d(int ncid, const std::string &varname, double *data,
                          int size);

  /**
   * @brief Read a 2D double array from a NetCDF variable.
   * @param ncid File ID
   * @param varname Variable name
   * @param data Output buffer (must be pre-allocated, row-major)
   * @param nx Size in x dimension
   * @param ny Size in y dimension
   * @return true on success
   */
  static bool read_var_2d(int ncid, const std::string &varname, double *data,
                          int nx, int ny);

  /**
   * @brief Read a 2D slice from a 3D variable (time, lat, lon).
   * @param ncid File ID
   * @param varname Variable name
   * @param data Output buffer
   * @param nx Size in x dimension
   * @param ny Size in y dimension
   * @param time_idx Time index to read
   * @return true on success
   */
  static bool read_var_3d_slice(int ncid, const std::string &varname,
                                double *data, int nx, int ny, int time_idx);

  /**
   * @brief Read a 1D integer array from a NetCDF variable.
   * @param ncid File ID
   * @param varname Variable name
   * @param data Output buffer
   * @param size Number of elements to read
   * @return true on success
   */
  static bool read_var_1d_int(int ncid, const std::string &varname, int *data,
                              int size);

  /**
   * @brief Write a 1D double array to a NetCDF variable.
   * @param ncid File ID
   * @param varname Variable name (must already exist)
   * @param data Input data buffer
   * @param size Number of elements to write
   * @return true on success
   */
  static bool write_var_1d(int ncid, const std::string &varname,
                           const double *data, int size);

  /**
   * @brief Write a 2D double array to a NetCDF variable.
   * @param ncid File ID
   * @param varname Variable name (must already exist)
   * @param data Input data buffer (row-major)
   * @param nx Size in x dimension
   * @param ny Size in y dimension
   * @return true on success
   */
  static bool write_var_2d(int ncid, const std::string &varname,
                           const double *data, int nx, int ny);

  // =========================================================================
  // Dimension/Variable Operations
  // =========================================================================

  /**
   * @brief Define a new dimension in a NetCDF file.
   * @param ncid File ID (must be in define mode)
   * @param dimname Dimension name
   * @param length Dimension length
   * @return Dimension ID, or -1 on failure
   */
  static int define_dim(int ncid, const std::string &dimname, int length);

  /**
   * @brief Get the size of a dimension.
   * @param ncid File ID
   * @param dimname Dimension name
   * @return Dimension length, or -1 if not found
   */
  static int get_dim_size(int ncid, const std::string &dimname);

  /**
   * @brief Check if a variable exists in the file.
   * @param ncid File ID
   * @param varname Variable name
   * @return true if variable exists
   */
  static bool has_var(int ncid, const std::string &varname);

  /**
   * @brief Define a new variable in a NetCDF file.
   * @param ncid File ID (must be in define mode)
   * @param varname Variable name
   * @param nctype NetCDF type (e.g., NC_DOUBLE)
   * @param dimids Vector of dimension IDs
   * @return Variable ID, or -1 on failure
   */
  static int define_var(int ncid, const std::string &varname, int nctype,
                        const std::vector<int> &dimids);

  /**
   * @brief End define mode and switch to data mode.
   *
   * Must be called after defining dimensions and variables, before
   * writing data.
   * @param ncid File ID
   * @return true on success
   */
  static bool end_def(int ncid);

private:
  static MPI_Comm s_comm;    ///< MPI communicator
  static int s_iosysid;      ///< PIO I/O system ID
  static bool s_initialized; ///< Initialization flag
  static int s_rank;         ///< MPI rank (for root-only logging)
};

/** @brief Default fill value for missing data. */
constexpr double FILLVALUE = 1.0e20;

} // namespace emulator

#endif // EMULATOR_IO_HPP
