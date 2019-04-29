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

  LinInterp(int ncol, int km1, int km2, Scalar minthresh) :
    m_ncol(ncol),
    m_km1(km1),
    m_km2(km2),
    m_km1_pack(scream::pack::npack<Pack>(km1)),
    m_km2_pack(scream::pack::npack<Pack>(km2)),
    m_minthresh(minthresh),
    m_policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_km2_pack)),
    m_indx_map("m_indx_map", ncol, scream::pack::npack<IntPack>(km2))
#ifndef NDEBUG
    , m_indx_map_dbg("m_indx_map_dbg", ncol, scream::pack::npack<IntPack>(km2))
#endif
  {}

  int km1_pack() const { return m_km1_pack; }
  int km2_pack() const { return m_km2_pack; }

  template<typename V>
  KOKKOS_INLINE_FUNCTION
  void setup(const MemberType& team,
             const V& x1,
             const V& x2) const
  {
    setup_nlogn(team, *this, scream::pack::repack<Pack::n>(x1), scream::pack::repack<Pack::n>(x2));
#ifndef NDEBUG
    setup_n2(team, *this, scream::pack::repack<Pack::n>(x1), scream::pack::repack<Pack::n>(x2));
#endif
  }

  // Linearly interpolate y(x1) onto coordinates x2
  template <typename V>
  KOKKOS_INLINE_FUNCTION
  void lin_interp(const MemberType& team,
                  const V& x1,
                  const V& x2,
                  const V& y1,
                  const V& y2) const
  {
    lin_interp_impl(team,
                    *this,
                    scream::pack::repack<Pack::n>(x1),
                    scream::pack::repack<Pack::n>(x2),
                    scream::pack::repack<Pack::n>(y1),
                    scream::pack::repack<Pack::n>(y2));
  }

  KOKKOS_INLINE_FUNCTION
  static void lin_interp_impl(const MemberType& team,
                              const LinInterp& liv,
                              const view_1d<const Pack>& x1, const view_1d<const Pack>& x2, const view_1d<const Pack>& y1,
                              const view_1d<Pack>& y2)
  {
    auto x1s = scream::pack::scalarize(x1);
    auto y1s = scream::pack::scalarize(y1);

    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, liv.m_km2_pack), [&] (Int k2) {
      const auto indx_pk = liv.m_indx_map(i, k2);
#ifndef NDEBUG
      const auto indx_pk_dbg = liv.m_indx_map_dbg(i, k2);
      for (int s = 0; s < Pack::n; ++s) {
        if (k2*Pack::n + s < liv.m_km2) {
          scream_kassert(indx_pk[s] == indx_pk_dbg[s]);
        }
      }
#endif
      const auto end_mask = indx_pk == liv.m_km1 - 1;
      if (end_mask.any()) {
        const auto not_end = !end_mask;
        scream_masked_loop(end_mask, s) {
          int k1 = indx_pk[s];
          y2(k2)[s] = y1s(k1) + (y1s(k1)-y1s(k1-1))*(x2(k2)[s]-x1s(k1))/(x1s(k1)-x1s(k1-1));
        }
        scream_masked_loop(not_end, s) {
          int k1 = indx_pk[s];
          y2(k2)[s] = y1s(k1) + (y1s(k1+1)-y1s(k1))*(x2(k2)[s]-x1s(k1))/(x1s(k1+1)-x1s(k1));
        }
      }
      else {
        Pack x1p, x1p1, y1p, y1p1;
        scream::pack::index_and_shift<1>(x1s, indx_pk, x1p, x1p1);
        scream::pack::index_and_shift<1>(y1s, indx_pk, y1p, y1p1);
        const auto& x2p = x2(k2);

        y2(k2) = y1p + (y1p1-y1p)*(x2p-x1p)/(x1p1-x1p);
      }

      y2(k2).set(y2(k2) < liv.m_minthresh, liv.m_minthresh);
    });
  }

#ifndef NDEBUG
  KOKKOS_INLINE_FUNCTION
  static void setup_n2(const MemberType& team, const LinInterp& liv, const view_1d<const Pack>& x1, const view_1d<const Pack>& x2)
  {
    auto x1s = scream::pack::scalarize(x1);
    auto idxs = liv.m_indx_map_dbg;

    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, liv.m_km2_pack), [&] (Int k2) {
      for (int s = 0; s < Pack::n; ++s) {
        if( x2(k2)[s] <= x1s(0) ) { // x2[k2] comes before x1[0]
          idxs(i, k2)[s] = 0;
        }
        else if( x2(k2)[s] >= x1s(liv.m_km1-1) ) { // x2[k2] comes after x1[-1]
          idxs(i, k2)[s] = liv.m_km1-1;
        }
        else {
          for (int k1 = 1; k1 < liv.m_km1; ++k1) { // scan over x1
            if( (x2(k2)[s]>=x1s(k1-1)) && (x2(k2)[s]<x1s(k1)) ) { // check if x2[k2] lies within x1[k1-1] and x1[k1]
              idxs(i, k2)[s] = k1-1;
            }
          }
        }
      }
    });
  }
#endif

  KOKKOS_INLINE_FUNCTION
  static void setup_nlogn(const MemberType& team, const LinInterp& liv, const view_1d<const Pack>& x1, const view_1d<const Pack>& x2)
  {
    auto x1s = scream::pack::scalarize(x1);

    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, liv.m_km2_pack), [&] (Int k2) {
      for (int s = 0; s < Pack::n; ++s) {
        const Scalar x1_indv = x2(k2)[s];
        auto begin = x1s.data();
        auto upper = begin + liv.m_km1;

        auto ub = util::upper_bound(begin, upper, x1_indv);
        int x1_idx = ub - begin;
        if (x1_idx > 0) {
          --x1_idx;
        }
        liv.m_indx_map(i, k2)[s] = x1_idx;
      }
    });
  }

  int m_ncol;
  int m_km1;
  int m_km2;
  int m_km1_pack;
  int m_km2_pack;
  Scalar m_minthresh;
  TeamPolicy m_policy;
  view_2d<IntPack> m_indx_map; // [x2-idx] -> x1-idx
#ifndef NDEBUG
  view_2d<IntPack> m_indx_map_dbg; // [x2-idx] -> x1-idx
#endif
};

} //namespace util
} //namespace scream

#endif
