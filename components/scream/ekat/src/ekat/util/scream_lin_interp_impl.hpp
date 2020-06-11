// Never include this header directly, only scream_lin_interp.hpp should include it

template <typename ScalarT, typename DeviceT>
LinInterp<ScalarT, DeviceT>::LinInterp(int ncol, int km1, int km2, Scalar minthresh) :
  m_ncol(ncol),
  m_km1(km1),
  m_km2(km2),
  m_km1_pack(scream::pack::npack<Pack>(km1)),
  m_km2_pack(scream::pack::npack<Pack>(km2)),
  m_minthresh(minthresh),
  m_policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_km2_pack)),
  m_indx_map("m_indx_map", ncol, scream::pack::npack<IntPack>(km2))
{}

template<typename ScalarT, typename DeviceT>
template<typename V>
KOKKOS_INLINE_FUNCTION
void LinInterp<ScalarT, DeviceT>::setup(const MemberType& team,
                                        const V& x1,
                                        const V& x2) const
{
  setup_impl(team, *this, scream::pack::repack<Pack::n>(x1), scream::pack::repack<Pack::n>(x2));
}

// Linearly interpolate y(x1) onto coordinates x2
template <typename ScalarT, typename DeviceT>
template <typename V>
KOKKOS_INLINE_FUNCTION
void LinInterp<ScalarT, DeviceT>::lin_interp(const MemberType& team,
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

template <typename ScalarT, typename DeviceT>
KOKKOS_INLINE_FUNCTION
void LinInterp<ScalarT, DeviceT>::lin_interp_impl(
  const MemberType& team,
  const LinInterp& liv,
  const view_1d<const Pack>& x1, const view_1d<const Pack>& x2, const view_1d<const Pack>& y1,
  const view_1d<Pack>& y2)
{
  auto x1s = scream::pack::scalarize(x1);
  auto y1s = scream::pack::scalarize(y1);

  const int i = team.league_rank();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, liv.m_km2_pack), [&] (Int k2) {
      const auto indx_pk = liv.m_indx_map(i, k2);
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

template <typename ScalarT, typename DeviceT>
KOKKOS_INLINE_FUNCTION
void LinInterp<ScalarT, DeviceT>::setup_impl(
  const MemberType& team, const LinInterp& liv, const view_1d<const Pack>& x1, const view_1d<const Pack>& x2)
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
