#ifndef SCREAM_LIN_INTERP_HPP
#define SCREAM_LIN_INTERP_HPP

#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_assert.hpp"
#include "share/util/scream_kokkos_utils.hpp"
#include "share/scream_pack.hpp"
#include "share/scream_pack_kokkos.hpp"

namespace scream {
namespace util {

template <typename ScalarT, typename DeviceT=DefaultDevice>
struct LinInterp
{
  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  using KT = KokkosTypes<Device>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  using ExeSpace    = typename KT::ExeSpace;
  using MemberType  = typename KT::MemberType;
  using TeamPolicy  = typename KT::TeamPolicy;

  // Testing has shown that li runs better on SKX with pk=1
  static constexpr int LI_PACKN = SCREAM_POSSIBLY_NO_PACK_SIZE;

  using Pack    = scream::pack::Pack<Scalar, LI_PACKN>;
  using IntPack = scream::pack::Pack<int, LI_PACKN>;

  //
  // ------ public API -------
  //

  LinInterp(int ncol, int km1, int km2, Scalar minthresh);

  int km1_pack() const { return m_km1_pack; }
  int km2_pack() const { return m_km2_pack; }

  template<typename V>
  KOKKOS_INLINE_FUNCTION
  void setup(const MemberType& team,
             const V& x1,
             const V& x2) const;

  // Linearly interpolate y(x1) onto coordinates x2
  template <typename V>
  KOKKOS_INLINE_FUNCTION
  void lin_interp(const MemberType& team,
                  const V& x1,
                  const V& x2,
                  const V& y1,
                  const V& y2) const;

  //
  // -------- Internal API, data ------
  //

  KOKKOS_INLINE_FUNCTION
  static void setup_impl(
    const MemberType& team, const LinInterp& liv, const view_1d<const Pack>& x1, const view_1d<const Pack>& x2);

  KOKKOS_INLINE_FUNCTION
  static void lin_interp_impl(
    const MemberType& team,
    const LinInterp& liv,
    const view_1d<const Pack>& x1, const view_1d<const Pack>& x2, const view_1d<const Pack>& y1,
    const view_1d<Pack>& y2);


  int m_ncol;
  int m_km1;
  int m_km2;
  int m_km1_pack;
  int m_km2_pack;
  Scalar m_minthresh;
  TeamPolicy m_policy;
  view_2d<IntPack> m_indx_map; // [x2-idx] -> x1-idx
};

#include "scream_lin_interp_impl.hpp"

} //namespace util
} //namespace scream

#endif
