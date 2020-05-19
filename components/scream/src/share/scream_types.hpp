#ifndef SCREAM_TYPES_HPP
#define SCREAM_TYPES_HPP

#include "scream_config.hpp"
#include "ekat/ekat_kokkos_helpers.hpp"

/*
 * Header contains globally useful types for Scream.
 * The global Int and Real types are defined here along
 * with a type dictionary for accessing commonly-used
 * Kokkos types.
 */

namespace scream {

#ifdef SCREAM_DOUBLE_PRECISION
using Real = double;
#else
using Real = float;
#endif

typedef int Int;

// The default device we use.
using DefaultDevice = Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;

// A device type to force host execution
using HostDevice = Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultHostExecutionSpace::memory_space>;

// Inject ekat's kokkos types and meta utils into scream namespace
using ekat::KokkosTypes;
using ekat::Unmanaged;
using ekat::MemoryManaged;
using ekat::MemoryUnmanaged;

// An enum to be used with object that have 'repository'-like behavior
enum class RepoState {
  Clean,
  Open,
  Closed
};

} // namespace scream

#endif // SCREAM_TYPES_HPP
