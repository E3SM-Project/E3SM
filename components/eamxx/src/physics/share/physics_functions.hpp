#ifndef PHYSICS_FUNCTIONS_HPP
#define PHYSICS_FUNCTIONS_HPP

#include "physics_constants.hpp"

#include "share/eamxx_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream {
namespace physics {

/*
 * Functions is a stateless struct used to encapsulate a
 * number of functions for shared physics. We use the ETI pattern for
 * these functions.
 *
 * Assumptions:
 *  - Kokkos team policies have a vector length of 1
 */

template <typename ScalarT, typename DeviceT>
struct Functions
{

  enum SaturationFcn { Polysvp1 = 0, MurphyKoop = 1};

  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  template <typename S>
  using BigPack = ekat::Pack<Scalar,SCREAM_PACK_SIZE>;
  template <typename S>
  using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

  using IntSmallPack = SmallPack<Int>;
  using Pack         = BigPack<Scalar>;
  using Spack        = SmallPack<Scalar>;

  template <typename S>
  using Mask = ekat::Mask<Pack::n>;

  template <typename S>
  using SmallMask = ekat::Mask<SmallPack<S>::n>;

  using Smask = SmallMask<Scalar>;

  using KT = KokkosTypes<Device>;

  using C = scream::physics::Constants<Scalar>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;
  template <typename S>
  using uview_2d = typename ekat::template Unmanaged<view_2d<S> >;

  using MemberType = typename KT::MemberType;

  using Workspace = typename ekat::WorkspaceManager<Spack, Device>::Workspace;

  //
  // --------- Functions ---------
  //

  // For all functions, range_mask represents a mask on the Spack variables which
  // determines which elements are active (not masked).

  //  compute saturation vapor pressure
  //  polysvp1 returned in units of pa.
  //  t is input in units of k.
  //  ice refers to saturation with respect to liquid (false) or ice (true)
  KOKKOS_FUNCTION
  static Spack polysvp1(const Spack& t, const bool ice, const Smask& range_mask, const char* caller=nullptr);

  //  compute saturation vapor pressure using Murphy and Koop(2005) formulation
  //  MurphyKoop_svp returned in units of pa.
  //  t is input in units of k.
  //  ice refers to saturation with respect to liquid (false) or ice (true)
  KOKKOS_FUNCTION
  static Spack MurphyKoop_svp(const Spack& t, const bool ice, const Smask& range_mask, const char* caller=nullptr);

  // Calls a function to obtain the saturation vapor pressure, and then computes
  // and returns the dry saturation mixing ratio, with respect to either liquid or ice,
  // depending on value of 'ice'
  KOKKOS_FUNCTION
  static Spack qv_sat_dry(const Spack& t_atm, const Spack& p_atm, const bool ice, const Smask& range_mask, const SaturationFcn func_idx = MurphyKoop, const char* caller=nullptr);

  // Calls qv_sat_dry and converts it to wet mixing ratio
  KOKKOS_FUNCTION
  static Spack qv_sat_wet(const Spack& t_atm, const Spack& p_atm, const bool ice, const Smask& range_mask, const Spack& dp_wet, const Spack& dp_dry, 
                          const SaturationFcn func_idx = MurphyKoop, const char* caller=nullptr);

  //checks temperature for negatives and NaNs
  KOKKOS_FUNCTION
  static void check_temperature(const Spack& t_atm, const char* caller, const Smask& range_mask);

};

} // namespace physics
} // namespace scream

// If a GPU build, make all code available to the translation unit; otherwise,
// ETI is used.
#ifdef EAMXX_ENABLE_GPU
# include "physics_saturation_impl.hpp"
#endif

#endif // PHYSICS_FUNCTIONS_HPP
