/**
 * @file test_config.hpp
 * @brief Test configuration utilities for emulator component tests.
 *
 * Provides paths and configuration for test execution, including
 * test data directory access and MPI rank utilities.
 */

#ifndef EMULATOR_TEST_CONFIG_HPP
#define EMULATOR_TEST_CONFIG_HPP

#include <cstdlib>
#include <iostream>
#include <string>

namespace emulator {
namespace testing {

/**
 * @brief Get path to test data directory.
 *
 * Returns the directory containing test data files.
 * Set at compile time via EMULATOR_TEST_DATA_DIR define.
 */
inline std::string get_test_data_dir() {
#ifdef EMULATOR_TEST_DATA_DIR
  return EMULATOR_TEST_DATA_DIR;
#else
  return ".";
#endif
}

/**
 * @brief Get path to a specific test data file.
 *
 * @param filename Name of the file (relative to test data directory)
 * @return Full path to the file
 */
inline std::string get_test_data_file(const std::string &filename) {
  auto path = get_test_data_dir() + "/" + filename;
  return path;
}

/**
 * @brief Check if running in CI environment.
 */
inline bool is_ci_environment() {
  const char *ci = std::getenv("CI");
  return ci != nullptr && std::string(ci) == "true";
}

/**
 * @brief Get the number of MPI ranks for the current test.
 */
inline int get_test_nprocs() {
  int nprocs = 1;
#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif
  return nprocs;
}

/**
 * @brief Get the current MPI rank.
 */
inline int get_test_rank() {
  int rank = 0;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  return rank;
}

} // namespace testing
} // namespace emulator

#endif // EMULATOR_TEST_CONFIG_HPP
