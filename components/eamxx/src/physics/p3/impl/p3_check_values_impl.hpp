
#ifndef P3_CHECK_VALUES_IMPL_HPP
#define P3_CHECK_VALUES_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*!------------------------------------------------------------------------------------
  Checks current values of prognotic variables for reasonable values and
  stops and prints values if they are out of specified allowable ranges.

  'check_consistency' means include trap for inconsistency in moments;
  otherwise, only trap for Q, T_atm, and negative Qx, etc.  This option is here
  to allow for Q<qsmall.and.N>nsmall or Q>qsmall.and.N<small which can be produced
  at the leading edges due to sedimentation and whose values are accpetable
  since lambda limiters are later imposed after SEDI (so one does not necessarily
  want to trap for inconsistency after sedimentation has been called).

  The value 'source_ind' indicates the approximate location in 'p3_main'
  from where 'check_values' was called before it resulted in a trap.
  -----------------------------------------------------------------------------------
*/
template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::check_values(const uview_1d<const Spack>& qv, const uview_1d<const Spack>& temp, const Int& ktop, const Int& kbot,
               const Int& timestepcount, const bool& force_abort, const Int& source_ind, const MemberType& team,
               const uview_1d<const Scalar>& col_loc)
{
  constexpr Scalar T_low  = 173.;
  constexpr Scalar T_high = 323.;
  constexpr Scalar Q_high = 40.e-3;
  constexpr Scalar Q_low  = 0.;

  Int kmin, kmax;
  ekat::impl::set_min_max(ktop, kbot, kmin, kmax, Spack::n);

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, kmax-kmin+1), [&] (int pk_) {

    bool trap{false};

    const int pk = kmin + pk_;

    const auto t_gt_low_bound   = temp(pk) > T_low;
    const auto t_lt_high_bound  = temp(pk) < T_high;

    const auto qv_gt_low_bound  = qv(pk) > Q_low;
    const auto qv_lt_high_bound = qv(pk) < Q_high;

    const auto t_out_bounds  = !(t_gt_low_bound && t_lt_high_bound);
    const auto qv_out_bounds = !(qv_gt_low_bound && qv_lt_high_bound);

    if (t_out_bounds.any()) {
      for (int s=0; s<Spack::n; ++s) {
        trap = true;
        //printf ("** WARNING IN P3_MAIN -- src, gcol, lon, lat, lvl, tstep, T: %d, %d, %13.6f, %13.6f, %d, %d, %13.6f\n"
        //,source_ind,static_cast<int>(col_loc(0)),col_loc(1),col_loc(2),pk,timestepcount,temp(pk)[s]);
      }
    }

    if (qv_out_bounds.any()) {
      for (int s=0; s<Spack::n; ++s) {
        // trap = .true.  !note, tentatively no trap, since Qv could be negative passed in to mp
        //printf ("** WARNING IN P3_MAIN -- src, gcol, lon, lat, lvl, tstep, Qv: %d, %d, %13.6f, %13.6f, %d, %d, %13.6f\n"
        //        ,source_ind,static_cast<int>(col_loc(0)),col_loc(1),col_loc(2),pk,timestepcount,qv(pk)[s]);
      }
    }

    if (trap && force_abort) {
      // printf ("**********************************************************\n");
      // printf ("** DEBUG TRAP IN P3_MAIN, s/r CHECK_VALUES -- source: %d\n ",source_ind);
      // printf ("**********************************************************\n");
      EKAT_KERNEL_REQUIRE( source_ind == 100 );
    }
  });
}

} // namespace p3
} // namespace scream

#endif
