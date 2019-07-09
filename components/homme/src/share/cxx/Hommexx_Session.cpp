/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#include "Hommexx_Session.hpp"
#include "ExecSpaceDefs.hpp"
#include "profiling.hpp"
#include "mpi/Comm.hpp"

// If possible, avoid including MpiContext,
// which forces to compile all the mpi folder
#ifndef HOMMEXX_NO_MPI_CONTEXT
  #include "mpi/MpiContext.hpp"
#endif

// If possible, avoid including Context,
// which forces to compile all the functors
#ifndef HOMMEXX_NO_CONTEXT
  #include "Context.hpp"
#endif

#include "vector/vector_pragmas.hpp"

#include <iostream>

namespace Homme
{

std::string active_avx_string () {
  std::string s;
#if defined __AVX512F__
  s += " - AVX512F";
#endif
#if defined __AVX2__
  s += " - AVX2";
#endif
#if defined __AVX__
  s += " - AVX";
#endif
  return s;
}

void initialize_hommexx_session ()
{
  /* Make certain profiling is only done for code we're working on */
  profiling_pause();

  /* Set Environment variables to control how many
   * threads/processors Kokkos uses */
  initialize_kokkos();

#ifdef HOMMEXX_NO_MPI_CONTEXT
  // In the unit tests, we may not need the MpiContext object.
  // Forcing to use MpiContext, also forces to compile the src code for
  // all the object that context stores (or else, undef ref linking errors).
  // If a test *does* need the context, that test can add the MpiContext and
  // all its sources to the set of files to compile.
  // For the initialization purpose, however, the whole MPI_COMM_WORLD is fine.
  const Comm comm(MPI_COMM_WORLD);
#else
  const auto comm = MpiContext::singleton().get_comm();
#endif
  if (comm.root()) {
    ExecSpace::print_configuration(std::cout, true);
    // Print configure-time settings.
#ifdef HOMMEXX_SHA1
    std::cout << "HOMMEXX SHA1: " << HOMMEXX_SHA1 << "\n";
#endif
    std::cout << "HOMMEXX AVX_VERSION: " << HOMMEXX_AVX_VERSION << "\n";
    std::cout << "HOMMEXX VECTOR_SIZE: " << VECTOR_SIZE << "\n";
    std::cout << "HOMMEXX vector tag: " << Scalar::label() << "\n";
    std::cout << "HOMMEXX active AVX set:" << active_avx_string() << "\n";
    std::cout << "HOMMEXX MPI_ON_DEVICE: " << HOMMEXX_MPI_ON_DEVICE << "\n";
    std::cout << "HOMMEXX CUDA_(MIN/MAX)_WARP_PER_TEAM: " << HOMMEXX_CUDA_MIN_WARP_PER_TEAM
              << " / " << HOMMEXX_CUDA_MAX_WARP_PER_TEAM << "\n";
#ifndef HOMMEXX_NO_VECTOR_PRAGMAS
    std::cout << "HOMMEXX has vector pragmas\n";
#else
    std::cout << "HOMMEXX doesn't have vector pragmas\n";
#endif

#ifdef HOMMEXX_CONFIG_IS_CMAKE
    std::cout << "HOMMEXX configured with CMake\n";
# ifdef HAVE_CONFIG_H
    std::cout << "HOMMEXX has config.h.c\n";
# endif
#else
    std::cout << "HOMMEXX provided best default values in Config.hpp\n";
#endif
  }
}

void finalize_hommexx_session ()
{
#ifndef HOMMEXX_NO_MPI_CONTEXT
  MpiContext::finalize_singleton();
#endif
#ifndef HOMMEXX_NO_CONTEXT
  Context::finalize_singleton();
#endif

  Kokkos::finalize();
}

} // namespace Homme
