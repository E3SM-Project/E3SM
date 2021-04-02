#ifndef SCREAM_COLUMN_OPS_HPP
#define SCREAM_COLUMN_OPS_HPP

#include <type_traits>
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_math.hpp"
#include "share/util/scream_combine_ops.hpp"
#include "share/scream_types.hpp"

#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/util//ekat_arch.hpp"
#include "ekat/ekat_pack.hpp"

namespace scream {

/*
 *  ColumnOps: a series of utility kernels that operate on a single column
 *
 *  This class is responsible of implementing common kernels used in
 *  scream to compute quantities at level midpoints and level interfaces.
 *  For instance, compute interface quantities from midpoints
 *  ones, or integrate over a column, or compute increments of midpoint
 *  quantities (which will be defined at interfaces).
 *  The kernels are meant to be launched from within a parallel region, with
 *  team policy. More precisely, they are meant to be called from the outer most
 *  parallel region. In other words, you should *not* be inside a TeamThreadRange
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
 *    col_ops::compute_midpoint_values(team,prod,output);
 *
 *  Note: all methods have a different impl, depending on whether PackSize>1.
 *        The if statement is evaluated at compile-time, so there is no run-time
 *        penalization. The only requirement is that both branches must compile.
 *        Also, the impl for PackSize>1 is not thread safe (no ||for or single),
 *        so it would *not* work on GPU. Hence, the whole class has a static_assert
 *        to make sure we have PackSize==1 if we are in a GPU build.
 */

template<typename DeviceType, typename ScalarType, int PackSize>
class ColumnOps {
public:
  // Expose template params
  using device_type = DeviceType;
  using scalar_type = ScalarType;
  using pack_type   = ekat::Pack<scalar_type,PackSize>;

  // Pack info
  enum : int {
    pack_size = PackSize,
  };
  using pack_info = ekat::PackInfo<pack_size>;

  // Kokkos types
  using exec_space = typename device_type::execution_space;
  using KT = ekat::KokkosTypes<device_type>;

  using TeamMember = typename KT::MemberType;

  template<typename T>
  using view_1d = typename KT::template view_1d<T>;


  static constexpr scalar_type one  () { return scalar_type(1); }
  static constexpr scalar_type zero () { return scalar_type(0); }

  // All functions have an 'input' provider template parameter. This can
  // either be a lambda (to allow input calculation on the fly) or a 1d view.
  // By default, it is a view.
  using DefaultProvider = view_1d<const pack_type>;

  // Runs the input lambda with a TeamThreadRange parallel for over [0,count) range
  template<typename Lambda>
  KOKKOS_INLINE_FUNCTION
  static void team_parallel_for (const TeamMember& team,
                                 const int count,
                                const Lambda& f) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,count),f);
  }

  // Runs the input lambda with a TeamThreadRange parallel for over [start,end) range
  template<typename Lambda>
  KOKKOS_INLINE_FUNCTION
  static void team_parallel_for (const TeamMember& team,
                                 const int start,
                                 const int end,
                                 const Lambda& f) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,start,end),f);
  }

  // Runs the input lambda with a TeamThreadRange parallel scan over [0,count) range
  template<typename Lambda>
  KOKKOS_INLINE_FUNCTION
  static void team_parallel_scan (const TeamMember& team,
                                  const int count,
                                  const Lambda& f) {
    Kokkos::parallel_scan(Kokkos::TeamThreadRange(team,count),f);
  }
  template<typename Lambda>
  KOKKOS_INLINE_FUNCTION
  static void team_parallel_scan (const TeamMember& team,
                                  const int start,
                                  const int end,
                                  const Lambda& f) {
    Kokkos::parallel_scan(Kokkos::TeamThreadRange(team,start,end),f);
  }

  // Runs the input lambda only for one team thread
  template<typename Lambda>
  KOKKOS_INLINE_FUNCTION
  static void team_single (const TeamMember& team,
                           const Lambda& f) {
    Kokkos::single(Kokkos::PerTeam(team),f);
  }

  template<typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  static void debug_checks (const int num_levels, const view_1d<pack_type>& x) {

    // Mini function to check that InputProvider supports op()(int)->pack_type,
    // and that the number of levels is compatible with pack_type and x's size.

    const auto npacks_min = pack_info::num_packs(num_levels);
    EKAT_KERNEL_ASSERT_MSG (num_levels>=0 && x.extent_int(0)>=npacks_min,
        "Error! Number of levels out of bounds.\n");

    using ret_type = decltype(std::declval<InputProvider>()(0));
    using raw_ret_type = typename std::remove_const<typename std::remove_reference<ret_type>::type>::type;

    static_assert(std::is_same<raw_ret_type,pack_type>::value,
      "Error! InputProvider should expose op()(int), returning a pack type.\n");
  }

  // Safety checks
  static_assert(!ekat::OnGpu<exec_space>::value || pack_size==1,
                "Error! ColumnOps impl would be buggy on gpu if VECTOR_SIZE>1.\n");

  // Compute X at level midpoints, given X at level interfaces
  template<CombineMode CM = CombineMode::Replace,
           typename InputProvider = DefaultProvider>
  KOKKOS_INLINE_FUNCTION
  static void
  compute_midpoint_values (const TeamMember& team,
                           const int num_mid_levels,
                           const InputProvider& x_i,
                           const view_1d<pack_type>& x_m,
                           const scalar_type alpha = one(),
                           const scalar_type beta = zero())
  {
    // Sanity checks
    debug_checks<InputProvider>(num_mid_levels,x_m);

    // For GPU (or any build with pack size 1), things are simpler
    if (pack_size==1) {
      team_parallel_for(team,num_mid_levels,
                        [&](const int& k) {
        auto tmp = ( x_i(k) + x_i(k+1) ) / 2.0;
        combine<CM>(tmp, x_m(k), alpha, beta);
      });
    } else {

      const auto NUM_MID_PACKS     = pack_info::num_packs(num_mid_levels);
      const auto LAST_MID_PACK     = pack_info::last_pack_idx(num_mid_levels);
      const auto LAST_INT_PACK     = pack_info::last_pack_idx(num_mid_levels+1);
      const auto LAST_INT_PACK_END = pack_info::last_vec_end(num_mid_levels+1);

      // Try to use SIMD operations as much as possible.
      team_parallel_for(team,NUM_MID_PACKS-1,
                        [&](const int k){
        // Shift's first arg is the scalar to put in the "empty" spot at the end
        pack_type tmp = ekat::shift_left(x_i(k+1)[0], x_i(k));
        tmp += x_i(k);
        tmp /= 2.0;
        combine<CM>(tmp, x_m(k), alpha, beta);
      });

      team_single (team, [&]() {
        // Last level pack treated separately, since int pack k+1 may not exist,
        // depending on whether NUM_MID_PACKS==NUM_INT_PACKS.
        // The shift left might not even "need" the 2nd arg (if num int packs == num mid packs),
        // but passing the last x_int value takes care of both scenarios.
        pack_type tmp = ekat::shift_left(x_i(LAST_INT_PACK)[LAST_INT_PACK_END-1], x_i(LAST_MID_PACK));

        tmp += x_i(LAST_MID_PACK);
        tmp /= 2.0;
        combine<CM>(tmp, x_m(LAST_MID_PACK), alpha, beta);
      });
    }
  }

  // Compute X at level interfaces, given X at level midpoints and top/bot bc values
  // Note: with proper bc, then x_int(x_mid(x_int))==x_int. However, this version
  //       weighs neighboring cells equally, even if one is larger than the other.
  // RECALL: k=0 is the model top, while k=num_mid_levels+1 is the surface!
  template<bool FixTop,
           CombineMode CM = CombineMode::Replace,
           typename InputProvider = DefaultProvider>
  KOKKOS_INLINE_FUNCTION
  static void
  compute_interface_values (const TeamMember& team,
                            const int num_mid_levels,
                            const InputProvider& x_m,
                            const scalar_type& bc,
                            const view_1d<pack_type>& x_i)
  {
    compute_interface_values_impl<pack_size,FixTop>(team,num_mid_levels,x_m,bc,x_i);
  }

  // Given X at level interfaces, compute dX at level midpoints.
  template<CombineMode CM = CombineMode::Replace,
           typename InputProvider = DefaultProvider>
  KOKKOS_INLINE_FUNCTION
  static void
  compute_midpoint_delta (const TeamMember& team,
                          const int num_mid_levels,
                          const InputProvider& x_i,
                          const view_1d<pack_type>& dx_m,
                          const scalar_type alpha = one(),
                          const scalar_type beta = zero())
  {
    // Sanity checks
    debug_checks<InputProvider>(num_mid_levels,dx_m);

    // For GPU (or any build with pack size 1), things are simpler
    if (pack_size==1) {
      team_parallel_for(team,num_mid_levels,
                        [&](const int& k) {
        auto tmp = x_i(k+1)-x_i(k);
        combine<CM>(tmp,dx_m(k),alpha,beta);
      });
    } else {
      const auto NUM_MID_PACKS     = pack_info::num_packs(num_mid_levels);
      const auto LAST_MID_PACK     = NUM_MID_PACKS - 1;
      const auto LAST_INT_PACK     = pack_info::last_pack_idx(num_mid_levels+1);
      const auto LAST_INT_PACK_END = pack_info::last_vec_end(num_mid_levels+1);

      // Try to use SIMD operations as much as possible.
      team_parallel_for(team,NUM_MID_PACKS-1,
                        [&](const int k){
        auto tmp = ekat::shift_left(x_i(k+1)[0], x_i(k));
        combine<CM>(tmp - x_i(k),dx_m(k),alpha,beta);
      });

      team_single(team, [&](){
        // Last pack does not necessarily have a next pack, so needs to be treated apart.
        auto tmp = ekat::shift_left(x_i(LAST_INT_PACK)[LAST_INT_PACK_END-1],x_i(LAST_MID_PACK));
        combine<CM>(tmp - x_i(LAST_MID_PACK),dx_m(LAST_MID_PACK),alpha,beta);
      });
    }
  }

  // Scan sum of a quantity defined at midpoints, to retrieve its integral at interfaces.
  // This function is the logical inverse of the one above
  // Notes:
  //  - FromTop: true means we scan over [0,num_mid_levels), while false is the opposite.
  //    RECALL: k=0 is the model top, while k=num_mid_levels+1 is the surface!
  //  - InputProvider: must provide an input al all mid levels
  //  - s0: used as bc value at k=0 (Forward) or k=num_mid_levels (Backward)
  template<bool FromTop, typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  static void
  column_scan (const TeamMember& team,
               const int num_mid_levels,
               const InputProvider& dx_m,
               const view_1d<pack_type>& x_i,
               const scalar_type& s0 = zero())
  {
    // Sanity checks
    debug_checks<InputProvider>(num_mid_levels+1,x_i);

    // Scan's impl is quite lengthy, so split impl into two fcns, depending on pack size.
    column_scan_impl<pack_size,false,FromTop>(team,num_mid_levels,dx_m,x_i,s0);
  }

protected:

  // Reduce a pack in two steps: first, reduce odd and even entries
  // separately, then combine the results. Can help with accuracy if
  // the pack contains alternating + and - entries of similar magnitude.
  KOKKOS_INLINE_FUNCTION
  static scalar_type interleaved_reduce_sum (const pack_type& p) {
    scalar_type s_even = zero();
    scalar_type s_odd  = zero();
    vector_simd
    for (int i=0; i<pack_size/2; ++i) {
      s_even += p[2*i];
      s_odd  += p[2*i+1];
    }
    return s_even+s_odd;
  }

  template<int PackLength,bool InterleavedReduction,
           bool FromTop,typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<PackLength==1>::type
  column_scan_impl (const TeamMember& team,
                    const int num_mid_levels,
                    const InputProvider& dx_m,
                    const view_1d<pack_type>& x_i,
                    const scalar_type& s0 = zero())
  {
    // If statement is evaluated at compile time, and compiled away.
    if (FromTop) {
      team_single(team,[&](){
        x_i(0)[0] = s0;
      });
      // No need for a barrier here

      team_parallel_scan(team,num_mid_levels,
                         [&](const int k, scalar_type& accumulator, const bool last) {
        accumulator += dx_m(k)[0];
        if (last) {
          x_i(k+1)[0] = s0 + accumulator;
        }
      });
    } else {
      team_single(team,[&](){
        x_i(num_mid_levels)[0] = s0;
      });
      // No need for a barrier here

      team_parallel_scan(team,num_mid_levels,
                         [&](const int k, scalar_type& accumulator, const bool last) {
        const auto k_bwd = num_mid_levels - k - 1;
        accumulator += dx_m(k_bwd)[0];
        if (last) {
          x_i(k_bwd)[0] = s0 + accumulator;
        }
      });
    }
  }

  template<int PackLength,bool InterleavedReduction,
           bool FromTop,typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<(PackLength>1)>::type
  column_scan_impl (const TeamMember& team,
                    const int num_mid_levels,
                    const InputProvider& dx_m,
                    const view_1d<pack_type>& x_i,
                    const scalar_type& s0 = zero())
  {
    // If statement is evaluated at compile time, and compiled away.
    const int NUM_MID_PACKS = pack_info::num_packs(num_mid_levels);
    const int NUM_INT_PACKS = pack_info::num_packs(num_mid_levels+1);
    const int LAST_INT_PACK = NUM_INT_PACKS - 1;
    const int LAST_MID_PACK = NUM_MID_PACKS - 1;

    if (FromTop) {
      // Te last interface pack is problematic, and needs to be pre-filled with zeros.
      team_single(team,[&](){
        x_i(LAST_INT_PACK) = zero();
      });
      team.team_barrier();

      // First, do a packed reduction of x_m, to get x_i(k) = dx_m(0)+...+dx_m(k-1)
      // Note: stop at NUM_INT_PACKS-2, to avoid OOB access.
      team_parallel_scan(team,NUM_INT_PACKS-1,
                         [&](const int& k, pack_type& accumulator, const bool last) {
        accumulator += dx_m(k);
        if (last) {
          x_i(k+1) = accumulator;
        }
      });
      team.team_barrier();

      // Now let s = s0 + reduce_sum(x_i(k)). The second term is the sum of all dx_m
      // on all "physical" levels on all previous packs.
      // Then, do the scan over the current pack: x_i(k)[i] = s + dx_m(k)[0,...,i-1]
      team_parallel_for(team,NUM_INT_PACKS,
                        [&](const int k) {
        scalar_type s = s0;
        // If k==0, x_i(k) does not contain any scan sum (and may contain garbage), so ignore it
        if (k!=0) {
          if (InterleavedReduction) {
            s += interleaved_reduce_sum (x_i(k));
          } else {
            ekat::reduce_sum<false>(x_i(k),s);
          }
        }
        x_i(k)[0] = s;

        // Note: if NUM_INT_PACKS>NUM_MID_PACKS, this_pack_end==1 on last int pack,
        // so we will *not* access dx_m(LAST_INT_PACK) (which would be OOB).
        const auto this_pack_end = pack_info::vec_end(num_mid_levels+1,k);
        for (int i=1; i<this_pack_end; ++i) {
          x_i(k)[i] = x_i(k)[i-1] + dx_m(k)[i-1];
        }
      });
    } else {
      // First, do a packed reduction of x_m, to get x_i(k) = dx_m(k)+...+dx_m(last_mid_pack)

      // The last interface pack as well as possibly the 2nd last (depending on whether
      // NUM_INT_PACKS==NUM_MID_PACKS or not) are problematic, and need to be pre-filled with zeros.
      team_single(team,[&](){
        x_i(LAST_INT_PACK) = zero();
        x_i(LAST_MID_PACK) = zero();
      });
      team.team_barrier();

      const auto LAST_MID_PACK_END = pack_info::last_vec_end(num_mid_levels);
      team_parallel_scan(team,1,NUM_MID_PACKS,
                         [&](const int& k, pack_type& accumulator, const bool last) {
        const auto k_bwd = NUM_MID_PACKS - k;
        accumulator += dx_m(k_bwd);
        if (k==0 && LAST_MID_PACK_END!=pack_size) {
          // Note: dx_m(LAST_MID_PACK) might contain junk, so we need to zero out
          // the trailing part of the accumulator (if any)
          for (int i=LAST_MID_PACK_END+1; i<pack_size; ++i) {
            accumulator = zero();
          }
        }

        if (last) {
          x_i(k_bwd-1) = accumulator;
        }
      });

      // Need to wait for the packed scan to be done before we move fwd
      team.team_barrier();

      // Now let s = s0 + reduce_sum(x_i(k)). The second term is the sum of all dx_m
      // on all "physical" levels on all packs above current one.
      // Then, do the scan over the current pack: x_i(k)[i] = s + dx_m(k)[i,...,pack_end]

      // If NUM_INT_PACKS>NUM_MID_PACKS, the last int pack only needs s0
      if (NUM_INT_PACKS>NUM_MID_PACKS) {
        team_single(team,[&]() {
          x_i(LAST_INT_PACK)[0] = s0;
        });
      }
      team_parallel_for(team,NUM_MID_PACKS,
                        [&](const int k) {
        const auto k_bwd = NUM_MID_PACKS - k - 1;
        const auto this_pack_end = pack_info::vec_end(num_mid_levels+1,k_bwd);

        scalar_type s = s0;
        if (k_bwd<NUM_MID_PACKS) {
          ekat::reduce_sum(x_i(k_bwd),s);
        }

        auto tally = s;
        auto& x_i_kbwd = x_i(k_bwd);
        for (int i=this_pack_end-1; i>=0; --i) {
          tally += dx_m(k_bwd)[i];
          x_i_kbwd[i] = tally;
        }
      });
    }
  }

  template<int PackLength, bool FixTop,
           CombineMode CM = CombineMode::Replace,
           typename InputProvider = DefaultProvider>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<PackLength==1>::type
  compute_interface_values_impl (const TeamMember& team,
                                 const int num_mid_levels,
                                 const InputProvider& x_m,
                                 const scalar_type& bc,
                                 const view_1d<pack_type>& x_i)
  {
    // Sanity checks
    debug_checks<InputProvider>(num_mid_levels+1,x_i);

    auto sign = [](const int k) -> scalar_type {
      return 1 - 2*(k%2);
    };

    auto lambda_odd = [&](const int k)->pack_type {
      return k%2==1 ? 2*x_m(k) : 0;
    };
    auto lambda_even = [&](const int k)->pack_type {
      return k%2==0 ? 2*x_m(k) : 0;
    };

    // Do a scan sum with 0 bc.
    if (FixTop) {
      team_single(team,[&](){
        x_i(0)[0] = zero();
      });
      // No need for a barrier here

      // To avoid adding + and - quantities (which may have similar magnitude),
      // split the scan sum of even and odd entries.
      team_parallel_scan(team,(num_mid_levels+1)/2,
                         [&](const int k, scalar_type& accumulator, const bool last) {
        accumulator += lambda_even(2*k)[0];
        if (last) {
          x_i(2*k+1)[0] = accumulator;
        }
      });
      team_parallel_scan(team,num_mid_levels/2,
                         [&](const int k, scalar_type& accumulator, const bool last) {
        accumulator += lambda_odd(2*k+1)[0];
        if (last) {
          x_i(2*k+2)[0] = accumulator;
        }
      });

      // Now, even entries of x_i contain scan_sum of previous even midpoints,
      // and similartly for odd entries. To obtain a complete scan sum, all
      // we need is to subtract the previous inteface from the current.
      // CAREFUL: do this starting from last interface!
      team.team_barrier();
      team_parallel_for(team,num_mid_levels,
                        [&](const int k) {
        const int k_bwd = num_mid_levels - k - 1;
        x_i(k_bwd+1) -= x_i(k_bwd);
        x_i(k_bwd+1) -= sign(k_bwd)*bc;
      });
      team_single(team,[&](){
        x_i(0)[0] = bc;
      });
    } else {
      team_single(team,[&](){
        x_i(num_mid_levels)[0] = 0;
      });
      // No need for a barrier here

      // To avoid adding + and - quantities (which may have similar magnitude),
      // split the scan sum of even and odd entries.
      team_parallel_scan(team,(num_mid_levels+1)/2,
                         [&](const int k, scalar_type& accumulator, const bool last) {
        const int k_bwd = (num_mid_levels+1)/2 - k - 1;
        accumulator += lambda_even(2*k_bwd)[0];
        if (last) {
          x_i(2*k_bwd)[0] = accumulator;
        }
      });
      team_parallel_scan(team,num_mid_levels/2,
                         [&](const int k, scalar_type& accumulator, const bool last) {
        const int k_bwd = num_mid_levels/2 - k - 1;
        accumulator += lambda_odd(2*k_bwd+1)[0];
        if (last) {
          x_i(2*k_bwd+1)[0] = accumulator;
        }
      });

      // Now, even entries of x_i contain scan_sum of previous even midpoints,
      // and similartly for odd entries. To obtain a complete scan sum, all
      // we need is to subtract the previous inteface from the current.
      team.team_barrier();
      team_parallel_for(team,num_mid_levels+1,
                        [&](const int k) {
        x_i(k-1) -= x_i(k);
        x_i(k-1) -= sign(k)*bc;
      });
      team_single(team,[&](){
        x_i(num_mid_levels)[0] = bc;
      });
    }
  }

  template<int PackLength, bool FixTop,
           CombineMode CM = CombineMode::Replace,
           typename InputProvider = DefaultProvider>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<(PackLength>1)>::type
  compute_interface_values_impl (const TeamMember& team,
                                 const int num_mid_levels,
                                 const InputProvider& x_m,
                                 const scalar_type& bc,
                                 const view_1d<pack_type>& x_i)
  {
    // Sanity checks
    debug_checks<InputProvider>(num_mid_levels+1,x_i);

    auto m1_pow_k = [](const int k)->int {
      return 1 - 2*(k%2);
    };
    // Store a pack of (-1)^k
    pack_type sign = 0;
    vector_simd
    for (int i=0; i<pack_size; ++i) {
      sign[i] = m1_pow_k(i);
    }

    // Scanned quantity is (-1)^n x_m(n)
    auto lambda = [&](const int k)->pack_type {
      return sign*x_m(k);
    };

    // Do a scan sum with 0 bc.
    // Note: the 2nd template arge tells column_scan_impl to perform the final
    // reduction on a single pack using 'interleaved_reduce_sum'
    column_scan_impl<pack_size,true,FixTop>(team,num_mid_levels,lambda,x_i);
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

};

} // namespace Homme

#endif // SCREAM_COLUMN_OPS_HPP
