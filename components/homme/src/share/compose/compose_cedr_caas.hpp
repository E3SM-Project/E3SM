#ifndef INCLUDE_COMPOSE_CEDR_CAAS_HPP
#define INCLUDE_COMPOSE_CEDR_CAAS_HPP

#include "cedr_caas.hpp"

#include "compose.hpp"
#ifndef COMPOSE_PORT

namespace homme {
namespace compose {

// This class is for the case of running the F90 code with HORIZ_OPENMP. We
// explicitly use Kokkos::Serial here so we can run the Kokkos kernels in the
// super class w/o triggering an expecution-space initialization error in
// Kokkos. Then the _horiz_omp versions of the CAAS impls parallelizes using the
// HORIZ_OPENMP framework instead of Kokkos.
struct CAAS : public cedr::caas::CAAS<Kokkos::Serial> {
  typedef cedr::caas::CAAS<Kokkos::Serial> Super;

  CAAS (const cedr::mpi::Parallel::Ptr& p, const cedr::Int nlclcells,
        const typename Super::UserAllReducer::Ptr& uar)
    : Super(p, nlclcells, uar)
  {}

  void run () override { run_horiz_omp(); }

private:
  void run_horiz_omp();
  void reduce_locally_horiz_omp();
  void finish_locally_horiz_omp();
};

} // namespace compose
} // namespace homme
#endif

#endif
