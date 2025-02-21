#ifndef SCREAM_COLUMN_OPS_HPP
#define SCREAM_COLUMN_OPS_HPP

#include "share/util/eamxx_combine_ops.hpp"
#include "share/eamxx_types.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_math.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/util//ekat_arch.hpp"
#include "ekat/ekat_pack.hpp"

#include <type_traits>

namespace scream {

/*
 *  ColumnOps: a series of utility kernels that operate on a single column
 *
 *  This class is responsible of implementing common kernels used in
 *  scream to compute quantities at level midpoints and level interfaces.
 *  For instance, compute interface quantities from midpoints
 *  ones, or integrate over a column, or compute increments of interface
 *  quantities (which will be defined at midpoints).
 *  The kernels are meant to be launched from within a parallel region, with
 *  team policy. More precisely, they are meant to be called from the outer most
 *  parallel region. In other words, you should *not* be inside a TeamVectorRange
 *  parallel loop when calling these kernels, since these kernels will attempt
 *  to create such loops. Furthermore, these kernels *assume* that the team policy
 *  vector length (on CUDA) is 1. We have no way of checking this (the vector length
 *  is stored in the policy, but not in the team member), so you must make
 *  sure that this is the case.
 *
 *  In the compute_* methods, InputProvider can either be a functor (e.g., a lambda)
 *  or a 1d view. The only requirement is that operator()(int) is defined,
 *  and returns a Pack<ScalarType,PackSize>.
 *  For instance, one could use a lambda to compute the midpoint average of
 *  the product of two interface quantities, like this:
 *
 *    using col_ops = ColumnOps<DefaultDevice,Real,N>;
 *    using pack_type = typename col_ops::pack_type;
 *
 *    auto prod = [&](const int k)->pack_type { return x(k)*y(k); }
 *    col_ops::compute_midpoint_values(team,nlevs,prod,output);
 *
 *  Notes:
 *   - all methods have a different impl, depending on whether PackSize>1.
 *     The if statement is evaluated at compile-time, so there is no run-time
 *     penalization. The only requirement is that both branches must compile.
 *   - some methods accept a non-type template argument of type CombineMode.
 *     This argument can be used to specify how the result of the calculation
 *     should be written in the output view. E.g., if CM=Update, the output
 *     view y will be updated as y = beta*y + alpha*f(x). The values alpha
 *     and beta are used only if CM needs them, and an error is thrown if the
 *     user specifies non-trivial alpha/beta when they are not needed.
 *     See eamxx_combine_ops.hpp for more details.
 *
 *  RECALL: k=0 is the model top, while k=num_mid_levels+1 is the surface!
 */

template<typename DeviceType, typename ScalarType>
class ColumnOps {
public:
  // Expose template params
  using device_type = DeviceType;
  using scalar_type = ScalarType;

  template<int PackSize>
  using pack_type = ekat::Pack<scalar_type,PackSize>;

  template<typename ScalarT>
  static constexpr bool is_simd () {
    return ekat::ScalarTraits<ScalarT>::is_simd;
  }

  template<typename ScalarT>
  KOKKOS_FUNCTION
  static constexpr int pack_size () {
    return sizeof(ScalarT) / sizeof(scalar_type);
  }

  // Kokkos types
  using exec_space = typename device_type::execution_space;
  using KT = ekat::KokkosTypes<device_type>;

  using TeamMember = typename KT::MemberType;

  template<typename ScalarT,typename MT = Kokkos::MemoryManaged>
  using view_1d = typename KT::template view_1d<ScalarT,MT>;

  KOKKOS_INLINE_FUNCTION
  static constexpr scalar_type one  () { return scalar_type(1); }
  KOKKOS_INLINE_FUNCTION
  static constexpr scalar_type zero () { return scalar_type(0); }

  template<typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static void debug_checks (const int num_levels, const view_1d<ScalarT,MT>& x) {

    // Mini function to check that InputProvider supports op()(int)->pack_type,
    // and that the number of levels is compatible with pack_type and x's size.

    EKAT_KERNEL_ASSERT_MSG (num_levels>=0 && pack_size<ScalarT>()*x.extent_int(0)>=num_levels,
        "Error! Number of levels out of bounds.\n");

    using ret_type = decltype(std::declval<InputProvider>()(0));
    using raw_ret_type = typename std::remove_const<typename std::remove_reference<ret_type>::type>::type;

    static_assert(std::is_same<raw_ret_type,ScalarT>::value,
      "Error! InputProvider should expose op()(int), returning a ScalarT.\n");

    static_assert(!ekat::OnGpu<exec_space>::value || pack_size<ScalarT>()==1,
                  "Error! Do not use PackSize>1 on GPU.\n");
  }

  // Compute X at level midpoints, given X at level interfaces
  template<typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static void
  compute_midpoint_values (const TeamMember& team,
                           const int num_mid_levels,
                           const InputProvider& x_i,
                           const view_1d<ScalarT,MT>& x_m)
  {
    compute_midpoint_values<CombineMode::Replace>(team,num_mid_levels,x_i,x_m,1,0);
  }
  // Compute X at level midpoints, given X at level interfaces
  template<CombineMode CM, typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static void
  compute_midpoint_values (const TeamMember& team,
                           const int num_mid_levels,
                           const InputProvider& x_i,
                           const view_1d<ScalarT,MT>& x_m,
                           const scalar_type alpha = one(),
                           const scalar_type beta = zero())
  {
    // Sanity checks
    debug_checks<InputProvider>(num_mid_levels,x_m);

    compute_midpoint_values_impl<CM,InputProvider,ScalarT,MT>(team,num_mid_levels,x_i,x_m,alpha,beta);
  }

  // Compute X at level interfaces, given X at level midpoints and top and bot bc.
  // Note: with proper bc, and with constant dz, then x_int(x_mid(x_int))==x_int.
  template<typename InputProvider1, typename InputProvider2, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static void
  compute_interface_values_linear (const TeamMember&          team,
                                   const int                  num_mid_levels,
                                   const InputProvider1&      x_m,
                                   const InputProvider2&      dz,
                                   const scalar_type&         bc_top,
                                   const scalar_type&         bc_bot,
                                   const view_1d<ScalarT,MT>& x_i)
  {
    // Sanity checks
    debug_checks<InputProvider1>(num_mid_levels+1,x_i);

    compute_interface_values_linear_impl(team,num_mid_levels,x_m,dz,bc_top,bc_bot,x_i);
  }

  // Compute X at level interfaces, given X at level midpoints and top or bot bc.
  // Notes:
  //  - If FixTop=true, the bc value is imposed at the top, and a scan sum from the
  //    top is used to retrieve the other interface values. If FixTop=false, the
  //    bc is imposed at the bottom, and a scan sum from the bottom is used to
  //    retrieve the other interface values.
  //  - with proper bc, then x_int(x_mid(x_int))==x_int.
  //  - CAVEAT: this implementation is not monotonic, in the sense that it can
  //            create maxima/minima larger than the input x_m, even at the
  //            "interior" interfaces. E.g., x_m==1, bc=0, yield x_i=[0,2,0,2,...].
  template<bool FixTop, typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static void
  compute_interface_values_compatible (const TeamMember& team,
                                       const int num_mid_levels,
                                       const InputProvider& x_m,
                                       const scalar_type& bc,
                                       const view_1d<ScalarT,MT>& x_i)
  {
    // Sanity checks
    debug_checks<InputProvider>(num_mid_levels+1,x_i);

    compute_interface_values_compatible_impl<FixTop>(team,num_mid_levels,x_m,bc,x_i);
  }

  // Given X at level interfaces, compute dX at level midpoints.
  template<typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static void
  compute_midpoint_delta (const TeamMember& team,
                          const int num_mid_levels,
                          const InputProvider& x_i,
                          const view_1d<ScalarT,MT>& dx_m,
                          const scalar_type alpha = one(),
                          const scalar_type beta = zero())
  {
    compute_midpoint_delta<CombineMode::Replace>(team,num_mid_levels,x_i,dx_m,alpha,beta);
  }

  // Given X at level interfaces, compute dX at level midpoints.
  template<CombineMode CM, typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static void
  compute_midpoint_delta (const TeamMember& team,
                          const int num_mid_levels,
                          const InputProvider& x_i,
                          const view_1d<ScalarT,MT>& dx_m,
                          const scalar_type alpha = one(),
                          const scalar_type beta = zero())
  {
    // Sanity checks
    debug_checks<InputProvider>(num_mid_levels,dx_m);

    compute_midpoint_delta_impl<CM>(team,num_mid_levels,x_i,dx_m,alpha,beta);
  }

  // Scan sum of a quantity defined at midpoints, to retrieve its integral at interfaces.
  // This function is the logical inverse of the one above.
  // Notes:
  //  - FromTop: true means we scan over [0,num_mid_levels), while false means the scan
  //             is over (num_mid_levels,0]. Recall that ilev=0 is the top of the model.
  //  - InputProvider: must provide an input al all mid levels
  //  - s0: used as bc value at k=0 (FromTop=true) or k=num_mid_levels (FromTop=false)
  template<bool FromTop, typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static void
  column_scan (const TeamMember& team,
               const int num_mid_levels,
               const InputProvider& dx_m,
               const view_1d<ScalarT,MT>& x_i,
               const scalar_type& s0 = zero())
  {
    // Sanity checks
    debug_checks<InputProvider>(num_mid_levels+1,x_i);

    // Scan's impl is quite lengthy, so split impl into two fcns, depending on pack size.
    column_scan_impl<FromTop>(team,num_mid_levels,dx_m,x_i,s0);
  }

protected:

  // ------------ Impls of midpoint_value ------------- //

  template<CombineMode CM, typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<(pack_size<ScalarT>()==1)>::type
  compute_midpoint_values_impl (const TeamMember& team,
                                const int num_mid_levels,
                                const InputProvider& x_i,
                                const view_1d<ScalarT,MT>& x_m,
                                const scalar_type alpha,
                                const scalar_type beta)
  {
    // For GPU (or any build with pack size 1), things are simpler
    team_parallel_for(team,num_mid_levels,
                      [&](const int& k) {
      auto tmp = ( x_i(k) + x_i(k+1) ) / 2.0;
      combine<CM>(tmp, x_m(k), alpha, beta);
    });
  }

  template<CombineMode CM, typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<(pack_size<ScalarT>()>1)>::type
  compute_midpoint_values_impl (const TeamMember& team,
                                const int num_mid_levels,
                                const InputProvider& x_i,
                                const view_1d<ScalarT,MT>& x_m,
                                const scalar_type alpha,
                                const scalar_type beta)
  {
    using pack_type = ScalarT;
    using pack_info = ekat::PackInfo<pack_size<ScalarT>()>;

    const auto NUM_MID_PACKS = pack_info::num_packs(num_mid_levels);
    const auto NUM_INT_PACKS = pack_info::num_packs(num_mid_levels+1);

    // It is easier to read if we check whether #int_packs==#mid_packs.
    // This lambda can be used in both cases to process packs that have a next pack
    auto shift_and_avg = [&](const int k){
      // Shift's first arg is the scalar to put in the "empty" spot at the end
      pack_type tmp = ekat::shift_left(x_i(k+1)[0], x_i(k));
      tmp += x_i(k);
      tmp /= 2.0;
      combine<CM>(tmp, x_m(k), alpha, beta);
    };
    if (NUM_MID_PACKS==NUM_INT_PACKS) {
      const auto LAST_PACK =  NUM_MID_PACKS-1;

      // Use SIMD operations only on NUM_MID_PACKS-1, since mid pack
      // does not have a 'next' one
      team_parallel_for(team,NUM_MID_PACKS-1,shift_and_avg);

      team_single (team, [&]() {
        // Last level pack treated separately, since int pack k+1 does not exist.
        // Shift's first arg is the scalar to put in the "empty" spot at the end.
        // In this case, we don't need it, so just uze zero.
        const auto& xi_last = x_i(LAST_PACK);
        pack_type tmp = ekat::shift_left(zero(), xi_last);

        tmp += xi_last;
        tmp /= 2.0;
        combine<CM>(tmp, x_m(LAST_PACK), alpha, beta);
      });
    } else {
      // We can use SIMD operations on all NUM_MID_PACKS mid packs,
      // since x_i(k+1) is *always* fine
      team_parallel_for(team,NUM_MID_PACKS,shift_and_avg);
    }
  }

  // ------------ Impls of midpoint_delta ------------- //

  template<CombineMode CM, typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<(pack_size<ScalarT>()==1)>::type
  compute_midpoint_delta_impl (const TeamMember& team,
                               const int num_mid_levels,
                               const InputProvider& x_i,
                               const view_1d<ScalarT,MT>& dx_m,
                               const scalar_type alpha,
                               const scalar_type beta)
  {
    // For GPU (or any build with pack size 1), things are simpler
    team_parallel_for(team,num_mid_levels,
                      [&](const int& k) {
      auto tmp = x_i(k+1)-x_i(k);
      combine<CM>(tmp,dx_m(k),alpha,beta);
    });
  }

  template<CombineMode CM, typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<(pack_size<ScalarT>()>1)>::type
  compute_midpoint_delta_impl (const TeamMember& team,
                               const int num_mid_levels,
                               const InputProvider& x_i,
                               const view_1d<ScalarT,MT>& dx_m,
                               const scalar_type alpha,
                               const scalar_type beta)
  {
    using pack_info = ekat::PackInfo<pack_size<ScalarT>()>;

    const auto NUM_MID_PACKS     = pack_info::num_packs(num_mid_levels);
    const auto NUM_INT_PACKS     = pack_info::num_packs(num_mid_levels+1);

    // It is easier to read if we check whether #int_packs==#mid_packs.
    // This lambda can be used in both cases to process packs that have a next pack
    auto shift_and_subtract = [&](const int k){
      auto tmp = ekat::shift_left(x_i(k+1)[0], x_i(k));
      combine<CM>(tmp - x_i(k),dx_m(k),alpha,beta);
    };
    if (NUM_MID_PACKS==NUM_INT_PACKS) {
      const auto LAST_PACK = NUM_MID_PACKS - 1;

      // Use SIMD operations only on NUM_MID_PACKS-1, since mid pack
      // does not have a 'next' one
      team_parallel_for(team,NUM_MID_PACKS-1,shift_and_subtract);

      // Last pack does not have a next one, so needs to be treated separately and serially.
      team_single(team, [&](){
        // Shift's first arg is the scalar to put in the "empty" spot at the end.
        // In this case, we don't need it, so just uze zero.
        const auto& xi_last = x_i(LAST_PACK);
        auto tmp = ekat::shift_left(zero(),xi_last);
        combine<CM>(tmp - xi_last,dx_m(LAST_PACK),alpha,beta);
      });
    } else {
      // We can use SIMD operations on all NUM_MID_PACKS mid packs,
      // since x_i(k+1) is *always* fine
      team_parallel_for(team,NUM_MID_PACKS,shift_and_subtract);
    }
  }

  // ------------ Impls of column_scan ------------- //

  template<bool FromTop,typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<(pack_size<ScalarT>()==1)>::type
  column_scan_impl (const TeamMember& team,
                    const int num_mid_levels,
                    const InputProvider& dx_m,
                    const view_1d<ScalarT,MT>& x_i,
                    const scalar_type& s0 = zero())
  {
    // If statement is evaluated at compile time, and compiled away.
    if (FromTop) {
      team_single(team,[&](){
        x_i(0) = s0;
      });
      // No need for a barrier here

      team_parallel_scan(team,num_mid_levels,
                         [&](const int k, ScalarT& accumulator, const bool last) {
        accumulator += dx_m(k);
        if (last) {
          x_i(k+1) = s0 + accumulator;
        }
      });
    } else {
      team_single(team,[&](){
        x_i(num_mid_levels) = s0;
      });
      // No need for a barrier here

      team_parallel_scan(team,num_mid_levels,
                         [&](const int k, ScalarT& accumulator, const bool last) {
        const auto k_bwd = num_mid_levels - k - 1;
        accumulator += dx_m(k_bwd);
        if (last) {
          x_i(k_bwd) = s0 + accumulator;
        }
      });
    }
  }

  template<bool FromTop, typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<(pack_size<ScalarT>()>1)>::type
  column_scan_impl (const TeamMember& team,
                    const int num_mid_levels,
                    const InputProvider& dx_m,
                    const view_1d<ScalarT,MT>& x_i,
                    const scalar_type& s0 = zero())
  {
    EKAT_KERNEL_ASSERT_MSG(pack_size<ScalarT>() <= num_mid_levels,
                           "Error! Currently, column_scan_impl() is not implemented for pack_size > num_mid_levels.");

    using pack_type = ScalarT;
    constexpr int PackLength = pack_size<ScalarT>();
    using pack_info = ekat::PackInfo<PackLength>;

    const int NUM_MID_PACKS = pack_info::num_packs(num_mid_levels);
    const int NUM_INT_PACKS = pack_info::num_packs(num_mid_levels+1);
    const int LAST_INT_PACK = NUM_INT_PACKS - 1;

    // If statement is evaluated at compile time, and compiled away.
    if (FromTop) {
      // Strategy:
      //  1. Do a packed reduction of x_m, to get x_i(k) = dx_m(0)+...+dx_m(k-1)
      //  2. Let s = s0 + reduce(x_i(k)). The rhs is the sum of all dx_m
      //     on all "physical" levels on all previous packs (plus bc).
      //  3. Do the scan over the current pack: x_i(k)[i] = s + dx_m(k)[0,...,i-1]

      // It is easier to read if we check whether #int_packs==#mid_packs.
      // This lambda can be used so long as there is a 'next' pack x_i(k+1);
      auto packed_scan_from_top = [&](const int& k, pack_type& accumulator, const bool last) {
        accumulator += dx_m(k);
        if (last) {
          x_i(k+1) = accumulator;
        }
      };

      if (NUM_MID_PACKS==NUM_INT_PACKS) {
        // Compute sum of previous packs (hence, stop at second-to-last pack)
        team_parallel_scan(team,NUM_MID_PACKS-1,packed_scan_from_top);
        team.team_barrier();

        // On each pack, reduce the result of the previous step, add the bc s0,
        // then do the scan sum within the current pack.
        team_parallel_for(team,NUM_INT_PACKS,
                          [&](const int k) {
          scalar_type s = s0;
          // If k==0, x_i(k) does not contain any scan sum (and may contain garbage), so ignore it
          if (k!=0) {
            ekat::reduce_sum<false>(x_i(k),s);
          }
          x_i(k)[0] = s;

          const auto this_pack_end = pack_info::vec_end(num_mid_levels+1,k);
          for (int i=1; i<this_pack_end; ++i) {
            x_i(k)[i] = x_i(k)[i-1] + dx_m(k)[i-1];
          }
        });
      } else {
        // Compute sum of previous packs
        team_parallel_scan(team,NUM_MID_PACKS,packed_scan_from_top);
        team.team_barrier();

        // On each pack, reduce the result of the previous step, add the bc s0,
        // then do the scan sum within the current pack.
        team_parallel_for(team,NUM_INT_PACKS,
                          [&](const int k) {
          scalar_type s = s0;
          // If k==0, x_i(k) does not contain any scan sum (and may contain garbage), so ignore it
          if (k!=0) {
            ekat::reduce_sum<false>(x_i(k),s);
          }
          x_i(k)[0] = s;

          // Note: for the last interface, this_pack_end==1, so we will *not* access
          //       dx_m(LAST_INT_PACK) (which would be OOB).
          const auto this_pack_end = pack_info::vec_end(num_mid_levels+1,k);
          for (int i=1; i<this_pack_end; ++i) {
            x_i(k)[i] = x_i(k)[i-1] + dx_m(k)[i-1];
          }
        });
      }
    } else {
      // Strategy:
      //  1. Do a packed reduction of x_m, to get x_i(k) = dx_m(k)+...+dx_m(num_mid_levs-1)
      //  2. Let s = s0 + reduce(x_i(k)). The rhs is the sum of all dx_m
      //     on all "physical" levels on all subsequent packs (plus bc).
      //  3. Do the scan over the current pack: x_i(k)[i] = s + dx_m(k)[0,...,i-1]

      // It is easier to read if we check whether #int_packs==#mid_packs.
      if (NUM_MID_PACKS==NUM_INT_PACKS) {
        // The easiest thing to do is to do the scan sum in the last pack,
        // where dx_m contains junk, then call this routine again, but
        // for num_mid_levels=(NUM_MID_PACKS-1)*PackLength.
        team_single(team,[&]() {
          auto& xi_last = x_i(NUM_INT_PACKS-1);
          const auto& dxm_last = dx_m(NUM_MID_PACKS-1);
          auto LAST_INT_VEC_END = pack_info::last_vec_end(num_mid_levels+1);
          xi_last[LAST_INT_VEC_END-1] = s0;
          for (int i=LAST_INT_VEC_END-2; i>=0; --i) {
            xi_last[i] = xi_last[i+1] + dxm_last[i];
          }
        });
        team.team_barrier();
        column_scan_impl<FromTop>(team,(NUM_MID_PACKS-1)*PackLength,dx_m,x_i,x_i(NUM_INT_PACKS-1)[0]);
      } else {
        // In this case, all packs of dx_m are full of meaningful values.
        auto packed_scan_from_bot = [&](const int& k, pack_type& accumulator, const bool last) {
          const auto k_bwd = NUM_MID_PACKS - k - 1;
          accumulator += dx_m(k_bwd);
          if (last) {
            x_i(k_bwd-1) = accumulator;
          }
        };
        team_parallel_scan(team,NUM_MID_PACKS-1,packed_scan_from_bot);

        // Need to wait for the packed scan to be done before we move fwd
        team.team_barrier();

        // Now let s = s0 + reduce_sum(x_i(k)). The second term is the sum of all dx_m
        // on all "physical" levels on all packs above current one.
        // Then, do the scan over the current pack: x_i(k)[i] = s + dx_m(k)[i,...,pack_end]

        // The last int pack only needs s0, and the second to last had no "following"
        // midpoints pack, and we didn't write anything in it during the scan sum,
        // so fill it with 0's
        team_single(team,[&]() {
          x_i(LAST_INT_PACK)[0] = s0;
          x_i(LAST_INT_PACK-1) = 0;
        });

        team_parallel_for(team,NUM_MID_PACKS,
                          [&](const int k) {
          const auto k_bwd = NUM_MID_PACKS - k - 1;

          scalar_type s = s0;
          if (k_bwd<NUM_MID_PACKS) {
            ekat::reduce_sum(x_i(k_bwd),s);
          }

          auto& xi_kbwd = x_i(k_bwd);
          const auto& dxm_kbwd = dx_m(k_bwd);
          xi_kbwd[PackLength-1] = s + dxm_kbwd[PackLength-1];
          for (int i=PackLength-2; i>=0; --i) {
            xi_kbwd[i] = xi_kbwd[i+1]+dxm_kbwd[i];
          }
        });
      }
    }
  }

  // ------------ Impls of compute_interface_values ------------- //

  template<typename InputProvider1, typename InputProvider2, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<(pack_size<ScalarT>()==1)>::type
  compute_interface_values_linear_impl (const TeamMember&          team,
                                        const int                  num_mid_levels,
                                        const InputProvider1&      x_m,
                                        const InputProvider2&      dz,
                                        const scalar_type&         bc_top,
                                        const scalar_type&         bc_bot,
                                        const view_1d<ScalarT,MT>& x_i)
  {
    // Pack size 1 yields a simple impl
    team_parallel_for(team,num_mid_levels+1, [&](const int k) {
      if (k==0)                   x_i(k) = bc_top;
      else if (k==num_mid_levels) x_i(k) = bc_bot;
      else                        x_i(k) = (x_m(k)*dz(k-1) + x_m(k-1)*dz(k))/(dz(k-1) + dz(k));
    });
  }

  template<typename InputProvider1, typename InputProvider2, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<(pack_size<ScalarT>()>1)>::type
  compute_interface_values_linear_impl (const TeamMember&          team,
                                        const int                  num_mid_levels,
                                        const InputProvider1&      x_m,
                                        const InputProvider2&      dz,
                                        const scalar_type&         bc_top,
                                        const scalar_type&         bc_bot,
                                        const view_1d<ScalarT,MT>& x_i)
  {
    using PackType           = ScalarT;
    constexpr int PackLength = pack_size<PackType>();
    using IntPackType        = ekat::Pack<Int,PackLength>;
    using PackInfo           = ekat::PackInfo<PackLength>;
    const auto num_int_packs = PackInfo::num_packs(num_mid_levels+1);

    const auto s_x_m = ekat::scalarize(x_m);
    const auto s_dz  = ekat::scalarize(dz);

    team_parallel_for(team,num_int_packs, [&](const int k) {
      // Setup Masks for boundary conditions
      const auto range        = ekat::range<IntPackType>(k*PackType::n);
      const auto at_top       = (range==0);
      const auto at_bottom    = (range==num_mid_levels);
      const auto in_interior  = (range>0 && range<num_mid_levels);

      // Calculate shift. Mask out 0 so shift does not
      // attempt to access index -1 or num_mid_levels.
      auto range_no_boundary = range;
      range_no_boundary.set(range<1 || range>=num_mid_levels, 1);
      PackType x_m_k, x_m_km1, dz_k, dz_km1;
      ekat::index_and_shift<-1>(s_x_m, range_no_boundary, x_m_k, x_m_km1);
      ekat::index_and_shift<-1>(s_dz,  range_no_boundary, dz_k,  dz_km1);

      // Calculate interface values
      x_i(k).set(at_top,      bc_top);
      x_i(k).set(at_bottom,   bc_bot);
      x_i(k).set(in_interior, (x_m_k*dz_km1 + x_m_km1*dz_k)/(dz_km1 + dz_k));
    });
  }

  template<bool FixTop, typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<(pack_size<ScalarT>()==1)>::type
  compute_interface_values_compatible_impl (const TeamMember& team,
                                            const int num_mid_levels,
                                            const InputProvider& x_m,
                                            const scalar_type& bc,
                                            const view_1d<ScalarT,MT>& x_i)
  {
    // Helper function that returns (-1)^k
    auto m1_pow_k = [](const int k) -> scalar_type {
      return 1 - 2*(k%2);
    };

    // The expression of x_i is (N=num_mid_levels)
    //   x_i(k+1) = (-1)^k [ -x_i(0) + 2\Sum_{n=0}^k (-1)^n x_m(n) ]
    //   x_i(k)   = (-1)^k [ (-1)^N x_i(N) + 2\Sum_{n=k}^{N-1} (-1)^n x_m(n) ]
    // for the cases where we fix top and bot respectively. In both cases,
    // we do a scan sum of (-1)^n x_m(n), using the alt_sign impl of column_scan.
    auto scan_input = [&](const int k) -> ScalarT {
      // The first term is -1 for k odd and 1 for k even.
      return 2*m1_pow_k(k) * x_m(k);
    };

    column_scan_impl<FixTop>(team,num_mid_levels,scan_input,x_i,0);

    if (FixTop) {

      // Need to add -x_i(0), adn multiply everything by (-1)^k
      team.team_barrier();
      team_parallel_for(team,num_mid_levels+1,
                        [&](const int k) {
        x_i(k) -= bc;
        x_i(k) *= m1_pow_k(k+1);
      });
    } else {
      // Need to add (-1)^N x_i(N), adn multiply everything by (-1)^k
      team.team_barrier();
      team_parallel_for(team,num_mid_levels+1,
                        [&](const int k) {
        x_i(k) += m1_pow_k(num_mid_levels) * bc;
        x_i(k) *= m1_pow_k(k);
      });
    }
  }

  template<bool FixTop, typename InputProvider, typename ScalarT, typename MT>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<(pack_size<ScalarT>()>1)>::type
  compute_interface_values_compatible_impl (const TeamMember& team,
                                            const int num_mid_levels,
                                            const InputProvider& x_m,
                                            const scalar_type& bc,
                                            const view_1d<ScalarT,MT>& x_i)
  {
    using pack_type = ScalarT;
    constexpr int PackLength = pack_size<ScalarT>();
    using pack_info = ekat::PackInfo<PackLength>;

    auto m1_pow_k = [](const int k)->int {
      return 1 - 2*(k%2);
    };
    // Store a pack of (-1)^k
    pack_type sign = 0;
    vector_simd
    for (int i=0; i<PackLength; ++i) {
      sign[i] = m1_pow_k(i);
    }

    // Scanned quantity is (-1)^n x_m(n)
    auto lambda = [&](const int k)->pack_type {
      return sign*x_m(k);
    };

    // Do a scan sum with 0 bc.
    // Note: the 2nd template arge tells column_scan_impl to perform the final
    // reduction on a single pack using 'interleaved_reduce_sum'
    column_scan_impl<FixTop>(team,num_mid_levels,lambda,x_i);
    team.team_barrier();

    const auto NUM_INT_PACKS = pack_info::num_packs(num_mid_levels+1);
    if (FixTop) {
      // Final formula:
      //   x_i(k+1) = (-1)^k [ -x_i(0) + 2\Sum_{n=0}^k (-1)^n x_m(n) ]
      // At this stage, x_i(k) contains the part within the \Sum

      team_parallel_for(team,NUM_INT_PACKS,
                        [&](const int k) {
        x_i(k) = sign*(bc - 2.0*x_i(k));
      });
    } else {
      // Final formula:
      //   x_i(k) = (-1)^k [ (-1)^N x_i(N) + 2\Sum_{n=k}^{N-1} (-1)^n x_m(n) ]
      // At this stage, x_i(k) contains the part within the \Sum

      team_parallel_for(team,NUM_INT_PACKS,
                        [&](const int k) {
        x_i(k) = sign*(bc*m1_pow_k(num_mid_levels) + 2.0*x_i(k));
      });
    }
  }

  // +--------------------------------------------+
  // |        Kokkos::parallel_X wrappers         |
  // +--------------------------------------------+

  // Runs the input lambda with a TeamVectorRange parallel for over [0,count) range
  template<typename Lambda>
  KOKKOS_INLINE_FUNCTION
  static void team_parallel_for (const TeamMember& team,
                                 const int count,
                                 const Lambda& f)
  {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,count),f);
  }

  // Runs the input lambda with a TeamVectorRange parallel for over [start,end) range
  template<typename Lambda>
  KOKKOS_INLINE_FUNCTION
  static void team_parallel_for (const TeamMember& team,
                                 const int start, const int end,
                                 const Lambda& f)
  {
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,start,end),f);
  }

  // Runs the input lambda with a TeamVectorRange parallel scan over [0,count) range
  template<typename Lambda>
  KOKKOS_INLINE_FUNCTION
  static void team_parallel_scan (const TeamMember& team,
                                  const int count,
                                  const Lambda& f)
  {
    team_parallel_scan(team,0,count,f);
  }

  // Runs the input lambda with a TeamVectorRange parallel scan over [start,end) range
  template<typename Lambda>
  KOKKOS_INLINE_FUNCTION
  static void team_parallel_scan (const TeamMember& team,
                                  const int start, const int end,
                                  const Lambda& f)
  {
    auto is_pow_of_2 = [](const int n)->bool {
      // This seems funky, but write down a pow of 2 and a non-pow of 2 in binary (both positive),
      // and you'll see why it works
      return n>0 && (n & (n-1))==0;
    };
    EKAT_KERNEL_REQUIRE_MSG (!ekat::OnGpu<typename device_type::execution_space>::value ||
                             is_pow_of_2(team.team_size()),
      "Error! Team-level parallel_scan on CUDA only works for team sizes that are power of 2.\n"
      "       You could try to reduce the team size to the previous pow of 2.\n");
    Kokkos::parallel_scan(Kokkos::TeamThreadRange(team,start,end),f);
  }

  // Runs the input lambda only for one team thread
  template<typename Lambda>
  KOKKOS_INLINE_FUNCTION
  static void team_single (const TeamMember& team,
                           const Lambda& f)
  {
    Kokkos::single(Kokkos::PerTeam(team),f);
  }

};

} // namespace Homme

#endif // SCREAM_COLUMN_OPS_HPP
