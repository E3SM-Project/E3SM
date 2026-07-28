/**
 * @file inference_context.cpp
 * @brief InferenceContext implementation.
 *
 * The only translation unit in the inference layer that includes <mpi.h>.
 */

#include "inference_context.hpp"

#include "inference_error.hpp"

#include <sstream>

#ifdef EMULATOR_HAVE_MPI
#include <mpi.h>
#include <unistd.h>
#endif

namespace emulator {
namespace inference {

void InferenceContext::set_grid(int nx_, int ny_, int num_global_cols_,
                                const int *gids, const double *lat_,
                                const double *lon_, int num_local_cols_) {
  EMULATOR_INFER_REQUIRE(num_local_cols_ >= 0,
                         "Negative local column count " << num_local_cols_
                                                        << ".");
  nx = nx_;
  ny = ny_;
  num_global_cols = num_global_cols_;
  const auto n = static_cast<std::size_t>(num_local_cols_);
  col_gids.assign(gids, gids + n);
  lat.assign(lat_, lat_ + n);
  lon.assign(lon_, lon_ + n);
}

std::string InferenceContext::to_string() const {
  std::ostringstream oss;
  oss << "rank " << rank << "/" << size << " on "
      << (node_name.empty() ? "?" : node_name) << ", grid " << nx << "x" << ny
      << ", " << num_local_cols() << " of " << num_global_cols << " columns";
  return oss.str();
}

#ifdef EMULATOR_HAVE_MPI

namespace {

std::string hostname() {
  char buf[256];
  if (::gethostname(buf, sizeof(buf)) != 0) {
    return "";
  }
  buf[sizeof(buf) - 1] = '\0';
  return std::string(buf);
}

} // namespace

InferenceContext make_context(int fcomm) {
  InferenceContext context;

  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized == 0) {
    return context; // MPI linked but never started (tools, unit tests)
  }

  MPI_Comm comm = MPI_Comm_f2c(static_cast<MPI_Fint>(fcomm));
  if (comm == MPI_COMM_NULL) {
    return context;
  }

  MPI_Comm_rank(comm, &context.rank);
  MPI_Comm_size(comm, &context.size);
  context.node_name = hostname();

  return context;
}

#else // !EMULATOR_HAVE_MPI

InferenceContext make_context(int fcomm) {
  (void)fcomm;
  return InferenceContext();
}

#endif // EMULATOR_HAVE_MPI

} // namespace inference
} // namespace emulator
