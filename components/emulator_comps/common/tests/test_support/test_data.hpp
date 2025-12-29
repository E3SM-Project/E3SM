/**
 * @file test_data.hpp
 * @brief Synthetic test data generators for emulator component tests.
 *
 * Provides utilities for creating test grids, field data, and comparison
 * functions for unit and integration tests.
 */

#ifndef EMULATOR_TEST_DATA_HPP
#define EMULATOR_TEST_DATA_HPP

#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <random>
#include <vector>

namespace emulator {
namespace testing {

//------------------------------------------------------------------------------
// Test Grid Structure
//------------------------------------------------------------------------------

/**
 * @brief Test grid data structure.
 *
 * Contains all data needed to set up a test grid for emulator components.
 */
struct TestGrid {
  int ncol;        ///< Local number of columns
  int ncol_global; ///< Global number of columns
  int nlev;        ///< Number of vertical levels
  int nx;          ///< Grid dimension in x
  int ny;          ///< Grid dimension in y

  std::vector<double> lat;  ///< Latitude (degrees)
  std::vector<double> lon;  ///< Longitude (degrees)
  std::vector<double> area; ///< Cell area (steradians)
  std::vector<int> gids;    ///< Global column IDs (1-indexed)
};

//------------------------------------------------------------------------------
// Grid Creation Functions
//------------------------------------------------------------------------------

/**
 * @brief Create a uniform lat-lon test grid.
 *
 * Creates a simple uniform grid for testing. Grid points are distributed
 * uniformly in latitude and longitude.
 *
 * @param ncol Number of columns
 * @param nlev Number of vertical levels (default 72)
 * @return TestGrid with synthetic data
 */
inline TestGrid create_test_grid(int ncol, int nlev = 72) {
  TestGrid grid;
  grid.ncol = ncol;
  grid.ncol_global = ncol;
  grid.nlev = nlev;

  // Compute grid dimensions (approximate square)
  grid.ny = static_cast<int>(std::sqrt(static_cast<double>(ncol)));
  if (grid.ny < 1)
    grid.ny = 1;
  grid.nx = (ncol + grid.ny - 1) / grid.ny;

  // Allocate arrays
  grid.lat.resize(ncol);
  grid.lon.resize(ncol);
  grid.area.resize(ncol);
  grid.gids.resize(ncol);

  // Create uniform lat-lon distribution
  double dlat = 180.0 / grid.ny;
  double dlon = 360.0 / grid.nx;
  double cell_area = (4.0 * M_PI) / ncol; // Approximate equal area

  for (int i = 0; i < ncol; ++i) {
    int iy = i / grid.nx;
    int ix = i % grid.nx;

    // Latitude: -90 to +90
    grid.lat[i] = -90.0 + (iy + 0.5) * dlat;

    // Longitude: 0 to 360
    grid.lon[i] = (ix + 0.5) * dlon;

    // Area (approximate)
    grid.area[i] = cell_area;

    // Global ID (1-indexed as expected by E3SM)
    grid.gids[i] = i + 1;
  }

  return grid;
}

/**
 * @brief Create a partitioned test grid for MPI tests.
 *
 * Creates a grid that is partitioned across MPI ranks using block distribution.
 *
 * @param ncol_global Global number of columns
 * @param rank MPI rank
 * @param nprocs Number of MPI processes
 * @param nlev Number of vertical levels (default 72)
 * @return TestGrid with local partition
 */
inline TestGrid create_partitioned_grid(int ncol_global, int rank, int nprocs,
                                        int nlev = 72) {
  // Block distribution with remainder
  int base_cols = ncol_global / nprocs;
  int remainder = ncol_global % nprocs;

  // Ranks 0..remainder-1 get one extra column
  int ncol_local = base_cols + (rank < remainder ? 1 : 0);
  int start_gid = rank * base_cols + std::min(rank, remainder);

  // Create global grid to extract from
  TestGrid global = create_test_grid(ncol_global, nlev);

  // Create local grid
  TestGrid grid;
  grid.ncol = ncol_local;
  grid.ncol_global = ncol_global;
  grid.nlev = nlev;
  grid.nx = global.nx;
  grid.ny = global.ny;

  // Allocate and fill local arrays
  grid.lat.resize(ncol_local);
  grid.lon.resize(ncol_local);
  grid.area.resize(ncol_local);
  grid.gids.resize(ncol_local);

  for (int i = 0; i < ncol_local; ++i) {
    int global_idx = start_gid + i;
    grid.lat[i] = global.lat[global_idx];
    grid.lon[i] = global.lon[global_idx];
    grid.area[i] = global.area[global_idx];
    grid.gids[i] = global_idx + 1; // 1-indexed
  }

  return grid;
}

//------------------------------------------------------------------------------
// Field Data Generators
//------------------------------------------------------------------------------

/**
 * @brief Create a field with random values.
 *
 * @param size Number of elements
 * @param min Minimum value (default 0.0)
 * @param max Maximum value (default 1.0)
 * @param seed Random seed (default 42)
 * @return Vector of random values
 */
inline std::vector<double> create_random_field(int size, double min = 0.0,
                                               double max = 1.0,
                                               int seed = 42) {
  std::vector<double> data(size);
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(min, max);

  for (int i = 0; i < size; ++i) {
    data[i] = dist(rng);
  }
  return data;
}

/**
 * @brief Create a field with constant value.
 *
 * @param size Number of elements
 * @param value Constant value
 * @return Vector of constant values
 */
inline std::vector<double> create_constant_field(int size, double value) {
  return std::vector<double>(size, value);
}

/**
 * @brief Create a field with linear gradient from 0 to 1.
 *
 * @param size Number of elements
 * @return Vector with values from 0 to 1
 */
inline std::vector<double> create_gradient_field(int size) {
  std::vector<double> data(size);
  for (int i = 0; i < size; ++i) {
    data[i] = static_cast<double>(i) / static_cast<double>(size - 1);
  }
  return data;
}

/**
 * @brief Create temperature-like field data.
 *
 * Creates a field with realistic temperature values (250-320 K).
 *
 * @param size Number of elements
 * @param seed Random seed
 * @return Vector of temperature-like values
 */
inline std::vector<double> create_temperature_field(int size, int seed = 42) {
  return create_random_field(size, 250.0, 320.0, seed);
}

/**
 * @brief Create pressure-like field data.
 *
 * Creates a field with realistic surface pressure values (95000-105000 Pa).
 *
 * @param size Number of elements
 * @param seed Random seed
 * @return Vector of pressure-like values
 */
inline std::vector<double> create_pressure_field(int size, int seed = 42) {
  return create_random_field(size, 95000.0, 105000.0, seed);
}

//------------------------------------------------------------------------------
// Comparison Utilities
//------------------------------------------------------------------------------

/**
 * @brief Compare two floating point values with tolerance.
 *
 * @param a First value
 * @param b Second value
 * @param rtol Relative tolerance (default 1e-10)
 * @param atol Absolute tolerance (default 1e-14)
 * @return true if values are approximately equal
 */
inline bool approx_equal(double a, double b, double rtol = 1e-10,
                         double atol = 1e-14) {
  double diff = std::abs(a - b);
  double threshold = atol + rtol * std::max(std::abs(a), std::abs(b));
  return diff <= threshold;
}

/**
 * @brief Compare two arrays with tolerance.
 *
 * @param a First array
 * @param b Second array
 * @param n Number of elements
 * @param rtol Relative tolerance
 * @param atol Absolute tolerance
 * @return true if all elements are approximately equal
 */
inline bool arrays_equal(const double *a, const double *b, int n,
                         double rtol = 1e-10, double atol = 1e-14) {
  for (int i = 0; i < n; ++i) {
    if (!approx_equal(a[i], b[i], rtol, atol)) {
      return false;
    }
  }
  return true;
}

/**
 * @brief Compare two vectors with tolerance.
 */
inline bool vectors_equal(const std::vector<double> &a,
                          const std::vector<double> &b, double rtol = 1e-10,
                          double atol = 1e-14) {
  if (a.size() != b.size())
    return false;
  return arrays_equal(a.data(), b.data(), static_cast<int>(a.size()), rtol,
                      atol);
}

/**
 * @brief Find the maximum absolute difference between two arrays.
 */
inline double max_diff(const double *a, const double *b, int n) {
  double max_d = 0.0;
  for (int i = 0; i < n; ++i) {
    max_d = std::max(max_d, std::abs(a[i] - b[i]));
  }
  return max_d;
}

} // namespace testing
} // namespace emulator

#endif // EMULATOR_TEST_DATA_HPP
