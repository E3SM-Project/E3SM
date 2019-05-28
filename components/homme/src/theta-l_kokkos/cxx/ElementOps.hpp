#ifndef HOMMEXX_ELEMENT_OPS_HPP
#define HOMMEXX_ELEMENT_OPS_HPP

#include "KernelVariables.hpp"
#include "HommexxEnums.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/VectorUtils.hpp"

namespace Homme {

/*
 *  ElementOps: a series of utility kernels inside an element
 *
 *  The name is mimicing the f90 module name, but ColumnOps would have
 *  been a more appropriate name, probably. In fact, this class is responsible
 *  of implementing common kernels used in the Theta model to compute quantities
 *  at level midpoints and level interfaces. For instance, compute interface
 *  quantities from midpoints ones, or integrate over a column, or compute
 *  increments of midpoint quantities (which will be defined at interfaces).
 *  The kernels are meant to be launched from within a parallel region, with
 *  team policy. More precisely, they are meant to be called from a parallel
 *  region dispatched over the number of thread in a single team. In other words,
 *  you should be inside a TeamThreadRange parallel loop before calling these
 *  kernels, but you should *not* be inside a ThreadVectorRange loop, since these
 *  kernels will attempt create such loops.
 */

class ElementOps {
public:

  using MIDPOINTS = ColInfo<NUM_PHYSICAL_LEV>;
  using INTERFACES = ColInfo<NUM_INTERFACE_LEV>;

  ElementOps () = default;

  KOKKOS_INLINE_FUNCTION
  void compute_midpoint_values (const KernelVariables& kv,
                                ExecViewUnmanaged<const Scalar [NUM_LEV_P]> x_i,
                                ExecViewUnmanaged<      Scalar [NUM_LEV]  > x_m) const
  {
    // Compute midpoint quanitiy. Note: the if statement is evaluated at compile time, so no penalization. Only requirement is both branches must compile.
    if (OnGpu<ExecSpace>::value) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_PHYSICAL_LEV),
                           [=](const int& ilev) {
        x_m(ilev) = (x_i(ilev) + x_i(ilev+1))/2.0;
      });
    } else {
      // Try to use SIMD operations as much as possible.
      for (int ilev=0; ilev<MIDPOINTS::LastPack; ++ilev) {
        Scalar tmp = x_i(ilev);
        tmp.shift_left(1);
        tmp[VECTOR_END] = x_i(ilev+1)[0];
        x_m(ilev) = (x_i(ilev) + tmp) / 2;
      }

      // Last level pack treated separately, since ilev+1 may throw depending if NUM_LEV=NUM_LEV_P
      Scalar tmp = x_i(MIDPOINTS::LastPack);
      tmp.shift_left(1);
      tmp[MIDPOINTS::LastVecEnd] = x_i(INTERFACES::LastPack)[INTERFACES::LastVecEnd];
      x_m(MIDPOINTS::LastPack) = (x_i(MIDPOINTS::LastPack) + tmp) / 2;
    }
  }

  // Computes the average of x_i*y_i and adds it to xy_m
  KOKKOS_INLINE_FUNCTION
  void update_midpoint_values_with_product (const KernelVariables& kv, const Real coeff,
                                            ExecViewUnmanaged<const Scalar [NUM_LEV_P]> x_i,
                                            ExecViewUnmanaged<const Scalar [NUM_LEV_P]> y_i,
                                            ExecViewUnmanaged<      Scalar [NUM_LEV]  > xy_m) const
  {
    // Compute midpoint quanitiy. Note: the if statement is evaluated at compile time, so no penalization. Only requirement is both branches must compile.
    if (OnGpu<ExecSpace>::value) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_PHYSICAL_LEV),
                           [=](const int& ilev) {
        xy_m(ilev) += coeff * (x_i(ilev)*y_i(ilev) + x_i(ilev+1)*y_i(ilev+1))/2.0;
      });
    } else {
      // Try to use SIMD operations as much as possible: the first NUM_LEV-1 packs can be vectorized
      for (int ilev=0; ilev<MIDPOINTS::LastPack; ++ilev) {
        Scalar tmp = x_i(ilev);
        tmp *= y_i(ilev);
        tmp.shift_left(1);
        tmp[VECTOR_END] = x_i(ilev+1)[0]*y_i(ilev+1)[0];
        tmp += x_i(ilev)*y_i(ilev);
        tmp *= coeff/2.0;
        xy_m(ilev) += tmp;
      }

      // Last level pack treated separately, since ilev+1 may throw depending if NUM_LEV=NUM_LEV_P
      Scalar tmp = x_i(MIDPOINTS::LastPack);
      tmp *= y_i(MIDPOINTS::LastPack);
      tmp.shift_left(1);
      tmp[VECTOR_END] = x_i(INTERFACES::LastPack)[0]*y_i(INTERFACES::LastPack)[0];
      tmp += x_i(MIDPOINTS::LastPack)*y_i(MIDPOINTS::LastPack);
      tmp *= coeff/2.0;
      xy_m(MIDPOINTS::LastPack) += tmp;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void compute_interface_values (const KernelVariables& kv,
                                 ExecViewUnmanaged<const Scalar [NUM_LEV]  > x_m,
                                 ExecViewUnmanaged<      Scalar [NUM_LEV_P]> x_i) const
  {
    // Compute interface quanitiy.
    if (OnGpu<ExecSpace>::value) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,1,NUM_PHYSICAL_LEV),
                           [=](const int& ilev) {
        x_i(ilev) = (x_m(ilev) + x_m(ilev-1)) / 2.0;
      });
      // Fix the top/bottom
      Kokkos::single(Kokkos::PerThread(kv.team),[&](){
        x_i(0) = x_m(0);
        x_i(NUM_INTERFACE_LEV-1) = x_m(NUM_PHYSICAL_LEV-1);
      });
    } else {
      // Try to use SIMD operations as much as possible: the last NUM_LEV-1 packs are treated uniformly, and can be vectorized
      for (int ilev=1; ilev<NUM_LEV; ++ilev) {
        Scalar tmp = x_m(ilev);
        tmp.shift_right(1);
        tmp[0] = x_m(ilev-1)[VECTOR_END];
        x_i(ilev) = (x_m(ilev) + tmp) / 2.0;
      }

      // First pack does not have a previous pack, and the extrapolation of the 1st interface is x_i = x_m.
      // Luckily, shift_right inserts leading 0's, so the formula is almost the same
      Scalar tmp = x_m(0);
      tmp.shift_right(1);
      x_i(0) = (x_m(0) + tmp) / 2.0;
      x_i(0)[0] = x_m(0)[0];

      // The last interface is x_i=x_m.
      x_i(INTERFACES::LastPack)[INTERFACES::LastVecEnd] = x_m(MIDPOINTS::LastPack)[MIDPOINTS::LastVecEnd];
    }
  }

  // Similar to the above, but uses layers thicknesses as averaging weights
  KOKKOS_INLINE_FUNCTION
  void compute_interface_values (const KernelVariables& kv,
                                 ExecViewUnmanaged<const Scalar [NUM_LEV]  > dp_m,
                                 ExecViewUnmanaged<const Scalar [NUM_LEV_P]> dp_i,
                                 ExecViewUnmanaged<const Scalar [NUM_LEV]  > x_m,
                                 ExecViewUnmanaged<      Scalar [NUM_LEV_P]> x_i) const
  {
    // Compute interface quanitiy.
    if (OnGpu<ExecSpace>::value) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,1,NUM_PHYSICAL_LEV),
                           [=](const int& ilev) {
        x_i(ilev) = (x_m(ilev)*dp_m(ilev) + x_m(ilev-1)*dp_m(ilev-1)) / (2.0*dp_i(ilev));
      });
      // Fix the top/bottom
      Kokkos::single(Kokkos::PerThread(kv.team),[&](){
        x_i(0) = x_m(0);
        x_i(NUM_INTERFACE_LEV-1) = x_m(NUM_PHYSICAL_LEV-1);
      });
    } else {
      // Try to use SIMD operations as much as possible: the last NUM_LEV-1 packs are treated uniformly, and can be vectorized
      for (int ilev=1; ilev<NUM_LEV; ++ilev) {
        Scalar tmp = x_m(ilev)*dp_m(ilev);
        tmp.shift_right(1);
        tmp[0] = x_m(ilev-1)[VECTOR_END]*dp_m(ilev-1)[VECTOR_END];
        x_i(ilev) = (x_m(ilev)*dp_m(ilev) + tmp) / (2.0*dp_i(ilev));
      }

      // First pack does not have a previous pack, and the extrapolation of the 1st interface is x_i = x_m.
      // Luckily, dp_i(0) = dp_m(0), and shift_right inserts leading 0's, so the formula is almost the same
      Scalar tmp = x_m(0)*dp_m(0);
      tmp.shift_right(1);
      x_i(0) = (x_m(0)*dp_m(0) + tmp) / (2.0*dp_i(0));
      x_i(0)[0] = x_m(0)[0];

      // The last interface is p_i=p_m.
      x_i(INTERFACES::LastPack)[INTERFACES::LastVecEnd] = x_m(MIDPOINTS::LastPack)[MIDPOINTS::LastVecEnd];
    }
  }

  KOKKOS_INLINE_FUNCTION
  void compute_midpoint_delta (const KernelVariables& kv,
                               ExecViewUnmanaged<const Scalar [NUM_LEV_P]> x_i,
                               ExecViewUnmanaged<      Scalar [NUM_LEV]  > dx_m) const
  {
    // Compute increment of interface values at midpoints.
    if (OnGpu<ExecSpace>::value) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,0,NUM_PHYSICAL_LEV),
                           [=](const int& ilev) {
        dx_m(ilev) = x_i(ilev+1)-x_i(ilev);
      });
    } else {
      // Try to use SIMD operations as much as possible. First NUM_LEV-1 packs can be treated the same
      for (int ilev=0; ilev<MIDPOINTS::LastPack; ++ilev) {
        Scalar tmp = x_i(ilev);
        tmp.shift_left(1);
        tmp[VECTOR_END] = x_i(ilev+1)[0];
        dx_m(ilev) = tmp - x_i(ilev);
      }

      // Last pack does not necessarily have a next pack, so needs to be treated a part.
      Scalar tmp = x_i(MIDPOINTS::LastPack);
      tmp.shift_left(1);
      tmp[MIDPOINTS::LastVecEnd] = x_i(INTERFACES::LastPack)[INTERFACES::LastVecEnd];
      dx_m(MIDPOINTS::LastPack) = tmp - x_i(MIDPOINTS::LastPack);
    }
  }

  template<BCType bcType = BCType::Zero>
  KOKKOS_INLINE_FUNCTION
  void compute_interface_delta (const KernelVariables& kv,
                                ExecViewUnmanaged<const Scalar [NUM_LEV]  > x_m,
                                ExecViewUnmanaged<      Scalar [NUM_LEV_P]> dx_i) const
  {
    static_assert (bcType==BCType::Zero || bcType == BCType::DoNothing,
                   "Error! Invalid bcType for interface delta calculation.\n");

    // Compute increment of midpoint values at interfaces. Top and bottom interfaces are set to 0.
    if (OnGpu<ExecSpace>::value) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,1,NUM_PHYSICAL_LEV),
                           [=](const int& ilev) {
        dx_i(ilev) = x_m(ilev)-x_m(ilev-1);
      });

      if (bcType==BCType::Zero) {
        // Fix the top/bottom
        Kokkos::single(Kokkos::PerThread(kv.team),[&](){
          dx_i(0) = dx_i(NUM_INTERFACE_LEV-1) = 0.0;
        });
      }
    } else {
      // Try to use SIMD operations as much as possible
      for (int ilev=1; ilev<NUM_LEV; ++ilev) {
        Scalar tmp = x_m(ilev);
        tmp.shift_right(1);
        tmp[0] = x_m(ilev-1)[VECTOR_END];
        dx_i(ilev) = x_m(ilev) - tmp;
      }

      // First pack does not have a previous pack. Luckily, shift_right inserts leading 0's, so the formula is the same, without the tmp[0] modification
      Scalar tmp = x_m(0);
      tmp.shift_right(1);
      dx_i(0) = x_m(0) - tmp;

      if (bcType==BCType::Zero) {
        // Fix the top/bottom levels
        dx_i(0)[0] = dx_i(INTERFACES::LastPack)[INTERFACES::LastVecEnd] = 0.0;
      }
    }
  }

  // Note: Forward=true means from k=0 to k=NUM_INTERFACE_LEV, false is the other way around
  // Note: the first value of sum_i (at 0 or NUM_INTERFACE_LEV, depending on Forward), is
  //       assumed to be VALID. In other words, the boundary condition of the integral must
  //       be set from OUTSIDE this kernel
  template<bool Forward,bool Inclusive,int LENGTH,typename Lambda>
  KOKKOS_INLINE_FUNCTION
  void column_scan (const KernelVariables& kv,
                    const Lambda& input_provider,
                    const ExecViewUnmanaged<Scalar [ColInfo<LENGTH>::NumPacks]>& sum) const
  {
    column_scan_impl<ExecSpace,Forward,Inclusive,LENGTH>(kv,input_provider,sum);
  }

  template<typename ExecSpaceType,bool Forward,bool Inclusive,int LENGTH,typename Lambda>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<!OnGpu<ExecSpaceType>::value>::type
  column_scan_impl (const KernelVariables& /* kv */,
                    const Lambda& input_provider,
                    const ExecViewUnmanaged<Scalar [ColInfo<LENGTH>::NumPacks]>& sum,
                    const Real s0 = 0.0) const
  {
    // It is easier to write two loops for Forward true/false. There's no runtime penalty,
    // since the if is evaluated at compile time, so no big deal.
    constexpr int lastPack = ColInfo<LENGTH>::NumPacks - 1;
    if (Forward) {
      // Running integral
      Real integration = s0;

      for (int ilev = 0; ilev<ColInfo<LENGTH>::NumPacks; ++ilev) {
        // In all but the last level pack, the loop is over the whole vector
        const int vec_end = ilev == lastPack ? ColInfo<LENGTH>::LastVecEnd
                                             : VECTOR_END;

        auto input = input_provider(ilev);
        // Integrate
        auto& sum_val = sum(ilev);
        sum_val[0] = integration + (Inclusive ? input[0] : 0.0);
        for (int iv = 1; iv <= vec_end; ++iv) {
          sum_val[iv] = sum_val[iv - 1] + (Inclusive ? input[iv] : input[iv-1]);
        }

        // Update running integral
        integration = sum_val[vec_end] + (Inclusive ? 0.0 : input[vec_end]);;
      }
    } else {
      // Running integral
      Real integration = s0;

      for (int ilev = lastPack; ilev >= 0; --ilev) {
        // In all but the last level pack, the loop is over the whole vector
        const int vec_start = ilev == lastPack ? ColInfo<LENGTH>::LastVecEnd
                                               : VECTOR_END;

        auto input = input_provider(ilev);
        // Integrate
        auto& sum_val = sum(ilev);
        sum_val[vec_start] = integration + (Inclusive ? input[vec_start] : 0.0);
        for (int iv = vec_start - 1; iv >= 0; --iv) {
          sum_val[iv] = sum_val[iv + 1] + (Inclusive ? input[iv] : input[iv+1]);
        }

        // Update running integral
        integration = sum_val[0] + (Inclusive ? 0.0 : input[0]);
      }
    }
  }

  template<typename ExecSpaceType,bool Forward,bool Inclusive,int LENGTH,typename Lambda>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<OnGpu<ExecSpaceType>::value>::type
  column_scan_impl (const KernelVariables& kv,
                    const Lambda& input_provider,
                    const ExecViewUnmanaged<Scalar [ColInfo<LENGTH>::NumPacks]>& sum,
                    const Real s0 = 0.0) const
  {
    // On GPU we rely on the fact that Scalar is basically double[1].
    static_assert (!OnGpu<ExecSpaceType>::value || ColInfo<LENGTH>::NumPacks==LENGTH, "Error! In a GPU build we expect VECTOR_SIZE=1.\n");

    if (Forward) {
      // accumulate input in [0,LENGTH].
      Dispatch<ExecSpaceType>::parallel_scan(kv.team, LENGTH,
                                            [&](const int k, Real& accumulator, const bool last) {
        accumulator += input_provider(k)[0];
        if (k==0) {
          accumulator += s0;
        }

        constexpr int last_idx = Inclusive ? LENGTH-1 : LENGTH-2;
        constexpr int offset = Inclusive ? 0 : 1;
        if (last && k<=last_idx) {
          sum(k+offset) = accumulator;
        }
      });
    } else {
      // accumulate input in [LENGTH,0].
      Dispatch<ExecSpaceType>::parallel_scan(kv.team, LENGTH,
                                            [&](const int k, Real& accumulator, const bool last) {
        // level must range in (LENGTH,0], while k ranges in [0, LENGTH).
        const int k_bwd = LENGTH-k-1;

        accumulator += input_provider(k_bwd)[0];
        if (k==0) {
          accumulator += s0;
        }

        constexpr int last_idx = Inclusive ? 0 : 1;
        constexpr int offset = Inclusive ? 0 : 1;
        if (last && k_bwd>=last_idx) {
          sum(k_bwd-offset) = accumulator;
        }
      });
    }
  }

  // Special case where input is on midpoints, but output is on interfaces.
  // In this case (for forward case), we perform sum(k+1) = sum(k) + provider(k),
  // for k=0,NUM_LEV. This can be done with an exclusive sum, using sum(0) as
  // initial value. Similarly for backward sum
  // Note: we are *assuming* that the first (or last, for bwd) entry of  sum
  //       contains the desired initial value
  template<bool Forward>
  KOKKOS_INLINE_FUNCTION
  void column_scan_mid_to_int (const KernelVariables& kv,
                               const ExecViewUnmanaged<const Scalar[NUM_LEV]>& input,
                               const ExecViewUnmanaged<Scalar [NUM_LEV_P]>& sum) const
  {
    auto input_provider = [&](const int ilev)->Scalar {
      return input(ilev);
    };

    using Specs = ColInfo<NUM_INTERFACE_LEV>;
    const Real s0 = Forward ? sum(0)[0] : sum(Specs::LastPack)[Specs::LastVecEnd];
    column_scan_impl<ExecSpace,Forward,false,NUM_INTERFACE_LEV>(kv,input_provider,sum,s0);
  }

  // Special case where input is on midpoints, but output is on interfaces.
  // In this case, we perform sum(k+1) = sum(k) + provider(k), k=0,NUM_LEV.
  // For this reason, we do an exclusive sum, with initial value sum(0).
  template<bool Forward,typename Lambda>
  KOKKOS_INLINE_FUNCTION
  void column_scan_mid_to_int (const KernelVariables& kv,
                               const Lambda& input_provider,
                               const ExecViewUnmanaged<Scalar [NUM_LEV_P]>& sum) const
  {
    if (Forward) {
      // It's safe to pass the output as it is, and claim is Exclusive over NUM_INTERFACE_LEV
      column_scan_impl<ExecSpace,true,false,NUM_INTERFACE_LEV>(kv,input_provider,sum,sum(0)[0]);
    } else {
      // Tricky: the input is not defined at NUM_INTEFACE_LEV-1. So we cast this scan sum
      //         into an inclusive sum over NUM_PHYSICAL_LEV, with output cropped to NUM_LEV packs.
      // Note: we also need to init sum at NUM_PHYSICAL_LEV-1
      ExecViewUnmanaged<Scalar[NUM_LEV]> sum_cropped(sum.data());
      const Real s0 = sum(INTERFACES::LastPack)[INTERFACES::LastVecEnd];
      sum_cropped(MIDPOINTS::LastPack)[MIDPOINTS::LastVecEnd] = s0;
      column_scan_impl<ExecSpace,false,true,NUM_PHYSICAL_LEV>(kv,input_provider,sum_cropped,s0);
    }
  }
};

} // namespace Homme

#endif // HOMMEXX_ELEMENT_OPS_HPP
