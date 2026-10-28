#ifndef E3SM_EMULATOR_GRID_DESC_HPP
#define E3SM_EMULATOR_GRID_DESC_HPP

#include <mpi.h>
#include <string>
#include <vector>

namespace e3sm::emulator {
/**
 * @brief C++-version EmulatorDimensionDesc
  -std::string name;
  -std::size_t local_extent;
  -std::size_t global_extent;
 */
struct DimensionDesc {
  std::string name;
  std::size_t local_extent;
  std::size_t global_extent;
};
std::string to_string(const DimensionDesc& desc);

/**
 * @brief Struct for holding the grid info for this process
 *  and it's possible dimentions;
   @fields:
   - int grid_type;
   - int nx;
   - int ny;
   - std::vector<DimensionDesc> dimensions;

   - std::vector<int> col_gids;
   - std::vector<double> lat;
   - std::vector<double> lon;
   - std::vector<double> area;
 *
 */
struct GridDesc {
  int grid_type;
  int nx;
  int ny;
  int num_global_columns;
  int num_local_columns;
  MPI_Comm comm;

  std::vector<DimensionDesc> dimensions;
  std::vector<int> col_gids;
  std::vector<double> lat;
  std::vector<double> lon;
  std::vector<double> area;
};
} // namespace e3sm::emulator
#endif
