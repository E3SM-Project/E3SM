#ifndef HOMMEXX_COLUMN_OPS_HPP
#define HOMMEXX_COLUMN_OPS_HPP

#include "KernelVariables.hpp"
#include "ErrorDefs.hpp"
#include "HommexxEnums.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/VectorUtils.hpp"

namespace Homme {

// Small helper function to combine a new value with an old one.
// The template argument help reducing the number of operations
// performed (the if is resolved at compile time). In the most
// complete form, the function performs
//    result = beta*result + alpha*newVal
// This routine should have no overhead compared to a manual
// update (assuming you call it with the proper CM)
template<CombineMode CM, typename Scalar1, typename Scalar2>
KOKKOS_FORCEINLINE_FUNCTION
void combine (const Scalar1& newVal, Scalar2& result,
              const Real alpha = 1.0, const Real beta = 1.0){
  switch (CM) {
    case CombineMode::Replace:
      result = newVal;
      break;
    case CombineMode::Scale:
      result = alpha*newVal;
      break;
    case CombineMode::Update:
      result *= beta;
      result += newVal;
      break;
    case CombineMode::ScaleUpdate:
      result *= beta;
      result += alpha*newVal;
      break;
    case CombineMode::ScaleAdd:
      result += alpha*newVal;
      break;
    case CombineMode::Add:
      result += newVal;
      break;
    case CombineMode::ProdUpdate:
      result *= newVal;
      break;
  }
}

/*
 *  ColumnOps: a series of utility kernels inside an element
 * 
 *  This class is responsible of implementing common kernels used in the
 *  preqx and theta models to compute quantities at level midpoints and
 *  level interfaces. For instance, compute interface quantities from midpoints
 *  ones, or integrate over a column, or compute increments of midpoint
 *  quantities (which will be defined at interfaces).
 *  The kernels are meant to be launched from within a parallel region, with
 *  team policy. More precisely, they are meant to be called from a parallel
 *  region dispatched over the number of thread in a single team. In other words,
 *  you should be inside a TeamThreadRange parallel loop before calling these
 *  kernels, but you should *not* be inside a ThreadVectorRange loop, since these
 *  kernels will attempt to create such loops.
 * 
 *  In the compute_* methods, InputProvider can either be a functor (e.g., a lambda) or a 1d view.
 *  The only requirement is that operator()(int)->Scalar is defined.
 *  For instance, one could use a lambda to compute the midpoint average of the product of two
 *  interface quantities, like this:
 *  
 *    auto prod = [&](const int ilev)->Scalar { return x(ilev)*y(ilev); }
 *    column_ops.compute_midpoint_values(kv,prod,output);
 */
  
class ColumnOps {
public:
  using MIDPOINTS  = ColInfo<NUM_PHYSICAL_LEV>;
  using INTERFACES = ColInfo<NUM_INTERFACE_LEV>;

  using DefaultMidProvider = ExecViewUnmanaged<const Scalar [NUM_LEV]>;
  using DefaultIntProvider = ExecViewUnmanaged<const Scalar [NUM_LEV_P]>;

  template<CombineMode CM>
  KOKKOS_INLINE_FUNCTION
  static constexpr bool needsAlpha () {
    return CM==CombineMode::Scale || CM==CombineMode::ScaleAdd || CM==CombineMode::ScaleUpdate;
  }

  template<CombineMode CM>
  KOKKOS_INLINE_FUNCTION
  static constexpr bool needsBeta () {
    return CM==CombineMode::Update || CM==CombineMode::ScaleUpdate;
  }

#ifndef NDEBUG
  template<CombineMode CM>
  KOKKOS_INLINE_FUNCTION
  static void sanity_check (const Real alpha, const Real beta) {
    assert ((needsAlpha<CM>() || alpha==1.0) &&
            "Error! Input alpha would be discarded by the requested combine mode.\n");
    assert ((needsBeta<CM>() || beta==0.0) &&
            "Error! Input beta would be discarded by the requested combine mode.\n");
  }
#endif

  template<CombineMode CM = CombineMode::Replace, typename InputProvider = DefaultIntProvider>
  KOKKOS_INLINE_FUNCTION
  static void compute_midpoint_values (const KernelVariables& kv,
                                const InputProvider& x_i,
                                const ExecViewUnmanaged<Scalar [NUM_LEV]>& x_m,
                                const Real alpha = 1.0, const Real beta = 0.0)
  {
#ifndef NDEBUG
    sanity_check<CM>(alpha,beta);
#endif

    // Compute midpoint quanitiy. Note: the if statement is evaluated at compile time, so no penalization. Only requirement is both branches must compile.
    if (OnGpu<ExecSpace>::value) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_PHYSICAL_LEV),
                           [&](const int& ilev) {
        Scalar tmp = (x_i(ilev) + x_i(ilev+1))/2.0;
        combine<CM>(tmp, x_m(ilev), alpha, beta);
      });
    } else {
      constexpr int LAST_MID_PACK     = MIDPOINTS::LastPack;
      constexpr int LAST_MID_PACK_END = MIDPOINTS::LastPackEnd;
      constexpr int LAST_INT_PACK     = INTERFACES::LastPack;
      constexpr int LAST_INT_PACK_END = INTERFACES::LastPackEnd;

      // Try to use SIMD operations as much as possible.
      for (int ilev=0; ilev<LAST_MID_PACK; ++ilev) {
        Scalar tmp = x_i(ilev);
        tmp.shift_left(1);
        tmp[VECTOR_END] = x_i(ilev+1)[0];
        tmp += x_i(ilev);
        tmp /= 2.0;
        combine<CM>(tmp, x_m(ilev), alpha, beta);
      }

      // Last level pack treated separately, since ilev+1 may throw depending if NUM_LEV=NUM_LEV_P
      Scalar tmp = x_i(LAST_MID_PACK);
      tmp.shift_left(1);
      tmp[LAST_MID_PACK_END] = x_i(LAST_INT_PACK)[LAST_INT_PACK_END];
      tmp += x_i(LAST_MID_PACK);
      tmp /= 2.0;
      combine<CM>(tmp, x_m(LAST_MID_PACK), alpha, beta);
    }
  }

  template<CombineMode CM = CombineMode::Replace, typename InputProvider = DefaultMidProvider>
  KOKKOS_INLINE_FUNCTION
  static void compute_interface_values (const KernelVariables& kv,
                                 const InputProvider& x_m,
                                 const ExecViewUnmanaged<Scalar [NUM_LEV_P]>& x_i,
                                 const Real alpha = 1.0, const Real beta = 0.0)
  {
#ifndef NDEBUG
    sanity_check<CM>(alpha,beta);
#endif

    // Compute interface quanitiy.
    if (OnGpu<ExecSpace>::value) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,1,NUM_PHYSICAL_LEV),
                           [=](const int& ilev) {
        Scalar tmp = (x_m(ilev) + x_m(ilev-1)) / 2.0;
        combine<CM>(tmp, x_i(ilev), alpha, beta);
      });
      // Fix the top/bottom
      Kokkos::single(Kokkos::PerThread(kv.team),[&](){
        combine<CM>(x_m(0), x_i(0), alpha, beta);
        combine<CM>(x_m(NUM_PHYSICAL_LEV-1), x_i(NUM_INTERFACE_LEV-1), alpha, beta);
      });
    } else {
      constexpr int LAST_MID_PACK     = MIDPOINTS::LastPack;
      constexpr int LAST_MID_PACK_END = MIDPOINTS::LastPackEnd;
      constexpr int LAST_INT_PACK     = INTERFACES::LastPack;
      constexpr int LAST_INT_PACK_END = INTERFACES::LastPackEnd;

      // Try to use SIMD operations as much as possible: the last NUM_LEV-1 packs are treated uniformly, and can be vectorized
      for (int ilev=1; ilev<NUM_LEV; ++ilev) {
        Scalar tmp = x_m(ilev);
        tmp.shift_right(1);
        tmp[0] = x_m(ilev-1)[VECTOR_END];
        tmp += x_m(ilev);
        tmp /= 2.0;
        combine<CM>(tmp, x_i(ilev), alpha, beta);
      }

      // First pack does not have a previous pack, and the extrapolation of the 1st interface is x_i = x_m.
      // Luckily, shift_right inserts leading 0's, so the formula is almost the same
      Scalar tmp = x_m(0);
      tmp.shift_right(1);
      tmp += x_m(0);
      tmp /= 2.0;
      combine<CM>(tmp, x_i(0), alpha, beta);

      // Fix top/bottom
      combine<CM>(x_m(0)[0], x_i(0)[0], alpha, beta);
      combine<CM>(x_m(LAST_MID_PACK)[LAST_MID_PACK_END], x_i(LAST_INT_PACK)[LAST_INT_PACK_END], alpha, beta);
    }
  }

  // Similar to the above, but uses midpoints/interface weights when computing the average
  template<CombineMode CM = CombineMode::Replace,
           typename WeightsMidProvider = DefaultMidProvider,
           typename WeightsIntProvider = DefaultIntProvider,
           typename InputProvider = DefaultMidProvider>
  KOKKOS_INLINE_FUNCTION
  static void compute_interface_values (const KernelVariables& kv,
                                 const WeightsMidProvider& weights_m,
                                 const WeightsIntProvider& weights_i,
                                 const InputProvider& x_m,
                                 const ExecViewUnmanaged<Scalar [NUM_LEV_P]>& x_i,
                                 const Real alpha = 1.0, const Real beta = 0.0)
  {
#ifndef NDEBUG
    sanity_check<CM>(alpha,beta);
#endif

    // Compute interface quanitiy.
    if (OnGpu<ExecSpace>::value) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,1,NUM_PHYSICAL_LEV),
                           [=](const int& ilev) {
        Scalar tmp = (x_m(ilev)*weights_m(ilev) + x_m(ilev-1)*weights_m(ilev-1)) / (2.0*weights_i(ilev));
        combine<CM>(tmp,x_i(ilev),alpha,beta);
      });
      // Fix the top/bottom
      Kokkos::single(Kokkos::PerThread(kv.team),[&](){
        combine<CM>(x_m(0),x_i(0),alpha,beta);
        combine<CM>(x_m(NUM_PHYSICAL_LEV-1),x_i(NUM_INTERFACE_LEV-1),alpha,beta);
      });
    } else {
      constexpr int LAST_MID_PACK     = MIDPOINTS::LastPack;
      constexpr int LAST_MID_PACK_END = MIDPOINTS::LastPackEnd;
      constexpr int LAST_INT_PACK     = INTERFACES::LastPack;
      constexpr int LAST_INT_PACK_END = INTERFACES::LastPackEnd;

      // Try to use SIMD operations as much as possible: the last NUM_LEV-1 packs are treated uniformly, and can be vectorized
      for (int ilev=1; ilev<NUM_LEV; ++ilev) {
        Scalar tmp = x_m(ilev)*weights_m(ilev);
        tmp.shift_right(1);
        tmp[0] = x_m(ilev-1)[VECTOR_END]*weights_m(ilev-1)[VECTOR_END];
        tmp += x_m(ilev)*weights_m(ilev);
        tmp /= 2.0*weights_i(ilev);
        combine<CM>(tmp,x_i(ilev),alpha,beta);
      }

      // First pack does not have a previous pack, and the extrapolation of the 1st interface is x_i = x_m.
      // Luckily, dp_i(0) = dp_m(0), and shift_right inserts leading 0's, so the formula is almost the same
      Scalar tmp = x_m(0)*weights_m(0);
      tmp.shift_right(1);
      tmp += x_m(0)*weights_m(0);
      tmp /= 2.0*weights_i(0);
      combine<CM>(tmp, x_i(0), alpha, beta);

      // Fix top/bottom
      combine<CM>(x_m(0)[0], x_i(0)[0], alpha, beta);
      combine<CM>(x_m(LAST_MID_PACK)[LAST_MID_PACK_END],x_i(LAST_INT_PACK)[LAST_INT_PACK_END],alpha,beta);
    }
  }

  template<CombineMode CM = CombineMode::Replace,
           typename InputProvider = DefaultIntProvider>
  KOKKOS_INLINE_FUNCTION
  static void compute_midpoint_delta (const KernelVariables& kv,
                               const InputProvider& x_i,
                               const ExecViewUnmanaged<Scalar [NUM_LEV]>& dx_m,
                               const Real alpha = 1.0, const Real beta = 0.0)
  {
    // Compute increment of interface values at midpoints.
    if (OnGpu<ExecSpace>::value) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,0,NUM_PHYSICAL_LEV),
                           [=](const int& ilev) {
        Scalar tmp = x_i(ilev+1)-x_i(ilev);
        combine<CM>(tmp,dx_m(ilev),alpha,beta);
      });
    } else {
#ifndef NDEBUG
      sanity_check<CM>(alpha,beta);
#endif
      constexpr int LAST_MID_PACK     = MIDPOINTS::LastPack;
      constexpr int LAST_MID_PACK_END = MIDPOINTS::LastPackEnd;
      constexpr int LAST_INT_PACK     = INTERFACES::LastPack;
      constexpr int LAST_INT_PACK_END = INTERFACES::LastPackEnd;

      // Try to use SIMD operations as much as possible. First NUM_LEV-1 packs can be treated the same
      for (int ilev=0; ilev<LAST_MID_PACK; ++ilev) {
        Scalar tmp = x_i(ilev);
        tmp.shift_left(1);
        tmp[VECTOR_END] = x_i(ilev+1)[0];
        combine<CM>(tmp - x_i(ilev),dx_m(ilev),alpha,beta);
      }

      // Last pack does not necessarily have a next pack, so needs to be treated a part.
      Scalar tmp = x_i(LAST_MID_PACK);
      tmp.shift_left(1);
      tmp[LAST_MID_PACK_END] = x_i(LAST_INT_PACK)[LAST_INT_PACK_END];
      combine<CM>(tmp - x_i(LAST_MID_PACK),dx_m(LAST_MID_PACK),alpha,beta);
    }
  }

  template<CombineMode CM = CombineMode::Replace,
           BCType bcType = BCType::Zero,
           typename InputProvider = DefaultMidProvider>
  KOKKOS_INLINE_FUNCTION
  static void compute_interface_delta (const KernelVariables& kv,
                                const InputProvider& x_m,
                                const ExecViewUnmanaged<Scalar [NUM_LEV_P]> dx_i,
                                const Real bcVal = 0.0,
                                const Real alpha = 1.0, const Real beta = 0.0)
  {
#ifndef NDEBUG
    sanity_check<CM>(alpha,beta);
#endif

    static_assert (bcType==BCType::Zero || bcType==BCType::Value || bcType == BCType::DoNothing,
                   "Error! Invalid bcType for interface delta calculation.\n");

    // Compute increment of midpoint values at interfaces. Top and bottom interfaces are set to 0.
    if (OnGpu<ExecSpace>::value) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,1,NUM_PHYSICAL_LEV),
                           [=](const int& ilev) {
        combine<CM>(x_m(ilev)-x_m(ilev-1),dx_i(ilev),alpha,beta);
      });

      if (bcType==BCType::Zero) {
        // Fix the top/bottom
        Kokkos::single(Kokkos::PerThread(kv.team),[&](){
          combine<CM>(0.0, dx_i(0)[0], alpha, beta);
          combine<CM>(0.0, dx_i(NUM_INTERFACE_LEV-1)[0], alpha, beta);
        });
      } else if (bcType==BCType::Value) {
        // Fix the top/bottom
        Kokkos::single(Kokkos::PerThread(kv.team),[&](){
          combine<CM>(bcVal, dx_i(0)[0], alpha, beta);
          combine<CM>(bcVal, dx_i(NUM_INTERFACE_LEV-1)[0], alpha, beta);
        });
      }
    } else {
      // Try to use SIMD operations as much as possible
      for (int ilev=1; ilev<NUM_LEV; ++ilev) {
        Scalar tmp = x_m(ilev);
        tmp.shift_right(1);
        tmp[0] = x_m(ilev-1)[VECTOR_END];
        combine<CM>(x_m(ilev) - tmp, dx_i(ilev), alpha, beta);
      }

      // First pack does not have a previous pack. Luckily, shift_right inserts leading 0's, so the formula is the same, without the tmp[0] modification
      Scalar tmp = x_m(0);
      tmp.shift_right(1);
      combine<CM>(x_m(0) - tmp, dx_i(0), alpha, beta);

      constexpr int LAST_PACK     = ColInfo<NUM_INTERFACE_LEV>::LastPack;
      constexpr int LAST_PACK_END = ColInfo<NUM_INTERFACE_LEV>::LastPackEnd;

      if (bcType==BCType::Zero) {
        // Fix the top/bottom levels
        combine<CM>(0.0, dx_i(0)[0], alpha, beta);
        combine<CM>(0.0, dx_i(LAST_PACK)[LAST_PACK_END], alpha, beta);
      } else if (bcType==BCType::Value) {
        // Fix the top/bottom levels
        combine<CM>(bcVal, dx_i(0)[0], alpha, beta);
        combine<CM>(bcVal, dx_i(LAST_PACK)[LAST_PACK_END], alpha, beta);
      }
    }
  }

  // Note: Forward=true means from k=0 to k=NUM_INTERFACE_LEV, false is the other way around
  // Note: the first value of sum_i (at 0 or NUM_INTERFACE_LEV, depending on Forward), is
  //       assumed to be VALID. In other words, the boundary condition of the integral must
  //       be set from OUTSIDE this kernel
  // Note: InputProvider could be a lambda or a 1d view.
  template<bool Forward, bool Inclusive, int LENGTH, typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  static void column_scan (const KernelVariables& kv,
                    const InputProvider& input_provider,
                    const ExecViewUnmanaged<Scalar [ColInfo<LENGTH>::NumPacks]>& sum,
                    const Real s0 = 0.0)
  {
    column_scan_impl<ExecSpace,Forward,Inclusive,LENGTH>(kv,input_provider,sum,s0);
  }

  template<typename ExecSpaceType,bool Forward,bool Inclusive,int LENGTH,typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<!OnGpu<ExecSpaceType>::value>::type
  column_scan_impl (const KernelVariables& /* kv */,
                    const InputProvider& input_provider,
                    const ExecViewUnmanaged<Scalar [ColInfo<LENGTH>::NumPacks]>& sum,
                    const Real s0 = 0.0)
  {
    constexpr int OFFSET             = Inclusive ? 0 : 1;
    constexpr int LOOP_RAW_SIZE      = LENGTH - OFFSET;
    constexpr int LOOP_LAST_PACK     = ColInfo<LOOP_RAW_SIZE>::LastPack;
    constexpr int LOOP_LAST_PACK_END = ColInfo<LOOP_RAW_SIZE>::LastPackEnd;
    constexpr int LOOP_SIZE          = ColInfo<LOOP_RAW_SIZE>::NumPacks;
    constexpr int LAST_PACK          = ColInfo<LENGTH>::LastPack;
    constexpr int LAST_PACK_END      = ColInfo<LENGTH>::LastPackEnd;

    // It is easier to write two loops for Forward true/false. There's no runtime penalty,
    // since the if is evaluated at compile time, so no big deal.
    if (Forward) {

      // Running integral
      Real integration = s0;

      for (int ilev = 0; ilev<LOOP_SIZE; ++ilev) {
        // In all but the last level pack, the loop is over the whole pack
        const int vec_end = (ilev == LOOP_LAST_PACK ? LOOP_LAST_PACK_END : VECTOR_END);

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

      // In an exclusive sum, the procedure above failed to update the last level
      if (!Inclusive) {
        sum(LAST_PACK)[LAST_PACK_END] = integration;
      }
    } else {
      // Running integral
      Real integration = s0;

      // In an exclusive sum, the procedure below would fail to add 
      // the input's last level to the output's second-to-last level
      if (!Inclusive) {
        integration += input_provider(LAST_PACK)[LAST_PACK_END];
      }

      for (int ipack=0; ipack<LOOP_SIZE; ++ipack) {
        const int ilev = LOOP_LAST_PACK-ipack;

        // In all but the last level pack, the loop is over the whole vector
        const int vec_start = (ilev == LOOP_LAST_PACK ? LOOP_LAST_PACK_END : VECTOR_END);

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

  template<typename ExecSpaceType,bool Forward,bool Inclusive,int LENGTH,typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if<OnGpu<ExecSpaceType>::value>::type
  column_scan_impl (const KernelVariables& kv,
                    const InputProvider& input_provider,
                    const ExecViewUnmanaged<Scalar [ColInfo<LENGTH>::NumPacks]>& sum,
                    const Real s0 = 0.0)
  {
    // On GPU we rely on the fact that Scalar is basically double[1].
    static_assert (!OnGpu<ExecSpaceType>::value || ColInfo<LENGTH>::NumPacks==LENGTH, "Error! In a GPU build we expect VECTOR_SIZE=1.\n");

    if (Forward) {
      // If exclusive, no need to go to access last input level
      constexpr int offset = Inclusive ? 0 : 1;
      constexpr int loop_size = LENGTH - offset;
      Dispatch<ExecSpaceType>::parallel_scan(kv.team, loop_size,
                                            [&](const int k, Real& accumulator, const bool last) {
        accumulator += input_provider(k)[0];
        if (k==0) {
          // First entry from the bottom: add initial value
          accumulator += s0;
        }

        if (last && k<loop_size) {
          sum(k+offset) = accumulator;
        }
      });
    } else {
      // If exclusive, no need to go to access first input level
      constexpr int offset = Inclusive ? 0 : 1;
      constexpr int loop_size = LENGTH - offset;
      Dispatch<ExecSpaceType>::parallel_scan(kv.team, loop_size,
                                            [&](const int k, Real& accumulator, const bool last) {
        // level must range in (LENGTH,0], while k ranges in [0, LENGTH).
        const int k_bwd = LENGTH - k - 1;

        accumulator += input_provider(k_bwd)[0];
        if (k==0) {
          // First entry from the top: add initial value
          accumulator += s0;
        }

        if (last) {
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
  template<bool Forward,typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  static void column_scan_mid_to_int (const KernelVariables& kv,
                               const InputProvider& input_provider,
                               const ExecViewUnmanaged<Scalar [NUM_LEV_P]>& sum)
  {
    if (Forward) {
      // It's safe to pass the output as it is, and claim is Exclusive over NUM_INTERFACE_LEV
      column_scan_impl<ExecSpace,true,false,NUM_INTERFACE_LEV>(kv,input_provider,sum,sum(0)[0]);
    } else {
      // Tricky: likely, the provider does not provide input at NUM_INTEFACE_LEV-1. So we cast this scan sum
      //         into an inclusive sum over NUM_PHYSICAL_LEV, with output cropped to NUM_LEV packs.
      // Note: we also need to init sum at NUM_PHYSICAL_LEV-1

      constexpr int LAST_MID_PACK     = MIDPOINTS::LastPack;
      constexpr int LAST_MID_PACK_END = MIDPOINTS::LastPackEnd;
      constexpr int LAST_INT_PACK     = INTERFACES::LastPack;
      constexpr int LAST_INT_PACK_END = INTERFACES::LastPackEnd;

      ExecViewUnmanaged<Scalar[NUM_LEV]> sum_cropped(sum.data());
      const Real s0 = sum(LAST_INT_PACK)[LAST_INT_PACK_END];
      sum_cropped(LAST_MID_PACK)[LAST_MID_PACK_END] = s0;
      column_scan_impl<ExecSpace,false,true,NUM_PHYSICAL_LEV>(kv,input_provider,sum_cropped,s0);
    }
  }

  template<int LENGTH, typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  static void column_reduction (const KernelVariables& kv,
                                const InputProvider& input,
                                Real& sum)
  {
#ifdef HOMMEXX_BFB_TESTING
    Dispatch<>::parallel_reduce(kv.team,Kokkos::ThreadVectorRange(kv.team,LENGTH),
                                [&](const int k, Real& accumulator) {
      if (!OnGpu<ExecSpace>::value) {
        const int ilev = k / VECTOR_SIZE;
        const int ivec = k % VECTOR_SIZE;
        accumulator += input(ilev)[ivec];
      } else {
        accumulator += input(k)[0];
      }
    },sum);
#else
    // To squeeze some perf out on CPU, do reduction at Scalar level,
    // then reduce the Scalar at the end. Be CAREFUL: the last pack may
    // contain some garbage if VECTOR_SIZE does not divide LENGTH!

    constexpr int  NUM_PACKS     = ColInfo<LENGTH>::NumPacks;
    constexpr int  LAST_PACK     = ColInfo<LENGTH>::LastPack;
    constexpr int  LAST_PACK_LEN = ColInfo<LENGTH>::LastPackLen;
    constexpr bool HAS_GARBAGE   = LAST_PACK_LEN!=VECTOR_SIZE;
    constexpr int  LOOP_LENGTH   = HAS_GARBAGE ? std::max(0,NUM_PACKS-1) : NUM_PACKS;

    Scalar packed_sum;
    Dispatch<>::parallel_reduce(kv.team,Kokkos::ThreadVectorRange(kv.team,LOOP_LENGTH),
                                [&](const int ilev, Scalar& accumulator){
      accumulator += input(ilev);
    },packed_sum);
    if (HAS_GARBAGE) {
      // Safety check: we assume on GPU VECTOR_SIZE=1.
      assert (!OnGpu<ExecSpace>::value);
      for (int k=0; k<LAST_PACK_LEN; ++k) {
        packed_sum[k] += input(LAST_PACK)[k];
      }
    }
    sum = packed_sum.reduce_add();
#endif
  }
};

} // namespace Homme

#endif // HOMMEXX_COLUMN_OPS_HPP
