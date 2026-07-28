/**
 * @file inference_context.hpp
 * @brief The resources a component hands the inference layer.
 *
 * Everything a model needs to know about *where* it is running and *which
 * part of the world it owns*, and nothing about what the model is.  Both
 * halves come from the MCT coupling layer: the component communicator, and
 * the horizontal decomposition the coupler already assigned to this rank.
 *
 * The rank and size are the *component's*, not the job's, which is the whole
 * point of passing them explicitly.  A distributed model that discovers its
 * rank from `SLURM_PROCID` / `SLURM_NTASKS` sees the entire coupled job, so
 * its process group waits forever for ocean and land ranks that will never
 * join it.
 */

#ifndef E3SM_EMULATOR_INFERENCE_CONTEXT_HPP
#define E3SM_EMULATOR_INFERENCE_CONTEXT_HPP

#include <string>
#include <vector>

namespace emulator {
namespace inference {

/**
 * @brief Ranks and horizontal decomposition.
 *
 * Default-constructed, it describes a serial run with no grid, which is what
 * tools and unit tests want.
 */
struct InferenceContext {
  // --- parallel resources, from the component communicator ---------------
  int rank = 0;          ///< Rank within the component communicator
  int size = 1;          ///< Size of the component communicator
  std::string node_name; ///< Hostname, for logging

  // --- horizontal decomposition, from the coupler ------------------------
  int nx = 0;                ///< Global longitude points (0 if unstructured)
  int ny = 0;                ///< Global latitude points
  int num_global_cols = 0;   ///< Global column count
  std::vector<int> col_gids; ///< 1-based global ids of this rank's columns
  std::vector<double> lat;   ///< Latitude of each local column [degrees]
  std::vector<double> lon;   ///< Longitude of each local column [degrees]

  int num_local_cols() const { return static_cast<int>(col_gids.size()); }

  /// @brief True if this rank is the one that speaks for the component.
  bool is_root() const { return rank == 0; }

  /// @brief Store the decomposition the coupler gave this rank.
  void set_grid(int nx_, int ny_, int num_global_cols_, const int *gids,
                const double *lat_, const double *lon_, int num_local_cols_);

  /// @brief One-line summary, for logs.
  std::string to_string() const;
};

/**
 * @brief Build a context from a Fortran MPI communicator handle.
 *
 * Built without MPI -- or under MPI but never launched, as a unit test is --
 * this returns a serial context, so the same call site works either way.
 *
 * @param fcomm Fortran communicator handle, as passed through the MCT layer.
 */
InferenceContext make_context(int fcomm);

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_INFERENCE_CONTEXT_HPP
