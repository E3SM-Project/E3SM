/**
 * @file decomp_info.hpp
 * @brief Domain decomposition metadata for parallel DataViews.
 */

#ifndef E3SM_EMULATOR_DATA_DECOMP_INFO_HPP
#define E3SM_EMULATOR_DATA_DECOMP_INFO_HPP

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace emulator {

/**
 * @brief Describes how a global grid is decomposed onto this MPI rank.
 *
 * Carries the local-to-global mapping (DOF array) needed by scorpio
 * for parallel I/O decomposition, plus the local/global point counts.
 *
 * The DOF array uses 1-based indexing to match Fortran/PIO conventions.
 */
struct DecompInfo {
  int npoints_local = 0;   ///< Number of grid points owned by this rank
  int npoints_global = 0;  ///< Total global grid points

  /**
   * @brief Global DOF indices for local points (1-based, PIO convention).
   *
   * dof[i] is the 1-based global index of local point i.
   * Size must equal npoints_local after initialization.
   */
  std::vector<int64_t> dof;

  /**
   * @brief Global grid dimensions (optional).
   *
   * For unstructured grids: {ncol}.
   * For structured grids: {nlon, nlat}.
   * For 3D data: {ncol, nlev} or {nlon, nlat, nlev}.
   * Used by IOAdapter for netCDF dimension definitions.
   */
  std::vector<int> global_dims;

  DecompInfo() = default;

  DecompInfo(int nlocal, int nglobal,
             std::vector<int64_t> dof_ = {},
             std::vector<int> dims = {})
      : npoints_local(nlocal), npoints_global(nglobal),
        dof(std::move(dof_)), global_dims(std::move(dims)) {}

  /**
   * @brief Validate internal consistency.
   * @throws std::runtime_error if dof size doesn't match npoints_local
   */
  void validate() const {
    if (!dof.empty() &&
        static_cast<int>(dof.size()) != npoints_local) {
      throw std::runtime_error(
          "DecompInfo: dof size (" + std::to_string(dof.size()) +
          ") != npoints_local (" + std::to_string(npoints_local) + ")");
    }
  }
};

} // namespace emulator

#endif // E3SM_EMULATOR_DATA_DECOMP_INFO_HPP
