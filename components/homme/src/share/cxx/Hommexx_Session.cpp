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

#include "Context.hpp"

#include "vector/vector_pragmas.hpp"

#include <iostream>

namespace Homme
{

// Default settings for homme's session
bool Session::m_inited = false;
bool Session::m_handle_kokkos = true;
bool Session::m_throw_instead_of_abort = false;

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

void print_homme_config_settings () {
  // Print configure-time settings.
#ifdef HOMMEXX_SHA1
  std::cout << "HOMMEXX SHA1: " << HOMMEXX_SHA1 << "\n";
#endif
  std::cout << "HOMMEXX VECTOR_SIZE: " << VECTOR_SIZE << "\n";
  std::cout << "HOMMEXX vector tag: " << Scalar::label() << "\n";
  std::cout << "HOMMEXX active AVX set:" << active_avx_string() << "\n";
  std::cout << "HOMMEXX MPI_ON_DEVICE: " << HOMMEXX_MPI_ON_DEVICE << "\n";
#ifdef HOMMEXX_CUDA_SHARE_BUFFER
  std::cout << "HOMMEXX CUDA_SHARE_BUFFER: on\n";
#else
  std::cout << "HOMMEXX CUDA_SHARE_BUFFER: off\n";
#endif
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

void initialize_hommexx_session ()
{
  // If hommexx session is not currently inited, then init it.
  if (!Session::m_inited) {
    /* Make certain profiling is only done for code we're working on */
    profiling_pause();

    /* Set Environment variables to control how many
     * threads/processors Kokkos uses */
    if (Session::m_handle_kokkos) {
      initialize_kokkos();
    }

    // Note: at this point, the Comm *should* already be created.
    const auto& comm = Context::singleton().get<Comm>();
    if (comm.root()) {
      ExecSpace::print_configuration(std::cout, true);
      print_homme_config_settings ();
    }

    Session::m_inited = true;
  }
}

void finalize_hommexx_session ()
{
  if (Session::m_inited) {
    Context::finalize_singleton();

    if (Session::m_handle_kokkos) {
      Kokkos::finalize();
    }
  }

  Session::m_inited = false;
}

} // namespace Homme
