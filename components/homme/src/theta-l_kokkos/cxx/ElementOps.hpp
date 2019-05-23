#ifndef HOMMEXX_ELEMENT_OPS_HPP
#define HOMMEXX_ELEMENT_OPS_HPP

#include "KernelVariables.hpp"
#include "PhysicalConstants.hpp"
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

  // Note: the last midpoint is in the pack NUM_LEV, at position (NUM_PHYSICAL_LEV-1)%VECTOR_SIZE;
  //       the last interface is in the pack NUM_LEV_P, one position after the last midpoint (mod VECTOR_SIZE).
  //       This is true regardless of whether NUM_LEV_P>NUM_LEV or not, so one formula works for all cases:
  //         - if NUM_LEV=NUM_LEV_P, then last_midpoint_vec_idx<VECTOR_SIZE-1, so the formula clearly works
  //         - otherwise, last_midpoint_vec_idx=VECTOR_SIZE-1, so last_interface_vec_idx=VECTOR_SIZE%VECTOR_SIZE=0, which is correct.
  // Note: these are just for convenience, not speedup. Begin compile time constants, they would be resolved at compile time anyways.
  static constexpr int last_midpoint_vec_idx  = (NUM_PHYSICAL_LEV + VECTOR_SIZE - 1) % VECTOR_SIZE;
  static constexpr int last_interface_vec_idx = (last_midpoint_vec_idx+1) % VECTOR_SIZE;
  static constexpr int LAST_LEV = NUM_LEV-1;
  static constexpr int LAST_LEV_P = NUM_LEV_P-1;
  static constexpr int VECTOR_END = VECTOR_SIZE-1;

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
      for (int ilev=0; ilev<LAST_LEV; ++ilev) {
        Scalar tmp = x_i(ilev);
        tmp.shift_left(1);
        tmp[VECTOR_END] = x_i(ilev+1)[0];
        x_m(ilev) = (x_i(ilev) + tmp) / 2;
      }

      // Last level pack treated separately, since ilev+1 may throw depending if NUM_LEV=NUM_LEV_P
      Scalar tmp = x_i(LAST_LEV);
      tmp.shift_left(1);
      tmp[last_midpoint_vec_idx] = x_i(LAST_LEV_P)[last_interface_vec_idx];
      x_m(LAST_LEV) = (x_i(LAST_LEV) + tmp) / 2;
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
      for (int ilev=0; ilev<LAST_LEV; ++ilev) {
        Scalar tmp = x_i(ilev);
        tmp *= y_i(ilev);
        tmp.shift_left(1);
        tmp[VECTOR_END] = x_i(ilev+1)[0]*y_i(ilev+1)[0];
        tmp += x_i(ilev)*y_i(ilev);
        tmp *= coeff/2.0;
        xy_m(ilev) += tmp;
      }

      // Last level pack treated separately, since ilev+1 may throw depending if NUM_LEV=NUM_LEV_P
      Scalar tmp = x_i(LAST_LEV);
      tmp *= y_i(LAST_LEV);
      tmp.shift_left(1);
      tmp[VECTOR_END] = x_i(LAST_LEV_P)[0]*y_i(LAST_LEV_P)[0];
      tmp += x_i(LAST_LEV)*y_i(LAST_LEV);
      tmp *= coeff/2.0;
      xy_m(LAST_LEV) += tmp;
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
      x_i(LAST_LEV_P)[last_interface_vec_idx] = x_m(LAST_LEV)[last_midpoint_vec_idx];
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
      x_i(LAST_LEV_P)[last_interface_vec_idx] = x_m(LAST_LEV)[last_midpoint_vec_idx];
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
      for (int ilev=0; ilev<LAST_LEV; ++ilev) {
        Scalar tmp = x_i(ilev);
        tmp.shift_left(1);
        tmp[VECTOR_END] = x_i(ilev+1)[0];
        dx_m(ilev) = tmp - x_i(ilev);
      }

      // Last pack does not necessarily have a next pack, so needs to be treated a part.
      Scalar tmp = x_i(LAST_LEV);
      tmp.shift_left(1);
      tmp[last_midpoint_vec_idx] = x_i(LAST_LEV_P)[last_interface_vec_idx];
      dx_m(LAST_LEV) = tmp - x_i(LAST_LEV);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void compute_interface_delta (const KernelVariables& kv,
                                ExecViewUnmanaged<const Scalar [NUM_LEV]  > x_m,
                                ExecViewUnmanaged<      Scalar [NUM_LEV_P]> dx_i) const
  {
    // Compute increment of midpoint values at interfaces. Top and bottom interfaces are set to 0.
    if (OnGpu<ExecSpace>::value) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,1,NUM_PHYSICAL_LEV),
                           [=](const int& ilev) {
        dx_i(ilev) = x_m(ilev)-x_m(ilev-1);
      });
      // Fix the top/bottom
      Kokkos::single(Kokkos::PerThread(kv.team),[&](){
        dx_i(0) = dx_i(NUM_INTERFACE_LEV-1) = 0.0;
      });
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

      // Fix the top/bottom levels
      dx_i(0)[0] = dx_i(LAST_LEV_P)[last_interface_vec_idx] = 0.0;
    }
  }

  // Note: Forward=true means from k=0 to k=NUM_INTERFACE_LEV, false is the other way around
  // Note: the first value of sum_i (at 0 or NUM_INTERFACE_LEV, depending on Forward), is
  //       assumed to be VALID. In other words, the boundary condition of the integral must
  //       be set from OUTSIDE this kernel
  template<bool Forward>
  KOKKOS_INLINE_FUNCTION
  void column_scan_sum (const KernelVariables& kv,
                        const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& x_m,
                        const ExecViewUnmanaged<      Scalar [NUM_LEV_P]>& sum_i) const
  {
    column_scan_sum<Forward,ExecSpace>(kv,x_m,sum_i);
  }

  template<bool Forward, typename ExecSpaceType>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<!OnGpu<ExecSpaceType>::value>::type
  column_scan_sum (const KernelVariables& /* kv */,
                   const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& x_m,
                   const ExecViewUnmanaged<      Scalar [NUM_LEV_P]>& sum_i) const
  {
    // It is easier to write two loops for Forward true/false. There's no runtime penalty,
    // since the if is evaluated at compile time, so no big deal.
    if (Forward) {
      // Running integral
      Real integration = sum_i(LAST_LEV_P)[last_interface_vec_idx];

      for (int ilev = 0; ilev<NUM_LEV; ++ilev) {
        // In all but the last level pack, the loop is over the whole vector
        const int vec_end = (ilev == (NUM_LEV - 1) ? last_midpoint_vec_idx
                                                   : VECTOR_SIZE - 1);

        // Integrate
        auto& sum_i_val = sum_i(ilev);
        sum_i_val[0] = integration;
        for (int iv = 1; iv <= vec_end; ++iv)
          sum_i_val[iv] = sum_i_val[iv - 1] + x_m(ilev)[iv - 1];

        // Update running integral
        integration = sum_i_val[vec_end] + x_m(ilev)[vec_end];
      }
      // Last interface
      sum_i(LAST_LEV_P)[last_interface_vec_idx] = integration;
    } else {
      // Running integral
      Real integration = sum_i(LAST_LEV_P)[last_interface_vec_idx];

      for (int ilev = NUM_LEV - 1; ilev >= 0; --ilev) {
        // In all but the last level pack, the loop is over the whole vector
        const int vec_start = (ilev == (NUM_LEV - 1) ? last_midpoint_vec_idx
                                                     : VECTOR_SIZE - 1);

        // Integrate
        auto& sum_i_val = sum_i(ilev);
        sum_i_val[vec_start] = integration;
        for (int iv = vec_start - 1; iv >= 0; --iv)
          sum_i_val[iv] = sum_i_val[iv + 1] + x_m(ilev)[iv + 1];

        // Update running integral
        integration = sum_i_val[0] + x_m(ilev)[0];
      }
    }
  }

  template<bool Forward, typename ExecSpaceType>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<OnGpu<ExecSpaceType>::value>::type
  column_scan_sum (const KernelVariables& kv,
                   const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& x_m,
                   const ExecViewUnmanaged<      Scalar [NUM_LEV_P]>& sum_i) const
  {
    // On GPU we rely on the fact that Scalar is basically double[1].
    static_assert (!OnGpu<ExecSpaceType>::value || NUM_LEV==NUM_PHYSICAL_LEV, "Error! In a GPU build we expect VECTOR_SIZE=1.\n");

    if (Forward) {
      // accumulate x_m in [0,NUM_PHYSICAL_LEV].
      Dispatch<ExecSpaceType>::parallel_scan(kv.team, NUM_INTERFACE_LEV-1,
                                            [&](const int k, Real& accumulator, const bool last) {
        accumulator += x_m(k)[0];

        if (last) {
          sum_i(k) = accumulator;
        }
      });
    } else {
      // accumulate x_m in [NUM_PHYSICAL_LEV,0].
      Dispatch<ExecSpaceType>::parallel_scan(kv.team, NUM_INTERFACE_LEV-1,
                                            [&](const int k, Real& accumulator, const bool last) {
        // level must range in [NUM_PHYSICAL_LEV-1,0], while k ranges in [0, NUM_PHYSICAL_LEV-2].
        const int level = NUM_INTERFACE_LEV-k-1;

        accumulator += x_m(level)[0];

        if (last) {
          sum_i(level) = accumulator;
        }
      });
    }
  }

  KOKKOS_INLINE_FUNCTION
  void compute_phi_i (const KernelVariables& kv,
                      const ExecViewUnmanaged<const Scalar [NP][NP][NUM_LEV]  >& vtheta_dp,
                      const ExecViewUnmanaged<const Scalar [NP][NP][NUM_LEV]  >& exner,
                      const ExecViewUnmanaged<const Scalar [NP][NP][NUM_LEV]  >& pnh,
                      const ExecViewUnmanaged<      Scalar [NP][NP][NUM_LEV_P]>& phi_i) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP; 
      const int jgp = idx % NP;
      compute_phi_i_impl<ExecSpace>(kv,Homme::subview(vtheta_dp,igp,jgp),
                                       Homme::subview(exner,igp,jgp),
                                       Homme::subview(pnh,igp,jgp),
                                       Homme::subview(phi_i,igp,jgp));
    });
  }

  KOKKOS_INLINE_FUNCTION
  void compute_phi_i (const KernelVariables& kv,
                      const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& vtheta_dp,
                      const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& exner,
                      const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& pnh,
                      const ExecViewUnmanaged<      Scalar [NUM_LEV_P]>& phi_i) const
  {
    compute_phi_i_impl<ExecSpace>(kv,vtheta_dp,exner,pnh,phi_i);
  }

  // KOKKOS_INLINE_FUNCTION
  // void compute_phinh_i (const KernelVariables& kv, const Real phis,
  //                       const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& vtheta_dp,
  //                       const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& p,
  //                       const ExecViewUnmanaged<      Scalar [NUM_LEV_P]>& phi_i) const {
  //   compute_phinh_i_impl<ExecSpace>(kv,phis,vtheta_dp,p,phi_i);
  // }

  template<typename ExecSpaceType>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<!OnGpu<ExecSpaceType>::value>::type
  compute_phi_i_impl (const KernelVariables& /* kv */,
                      const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& vtheta_dp,
                      const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& exner,
                      const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& pnh,
                      const ExecViewUnmanaged<      Scalar [NUM_LEV_P]>& phi_i) const
  {
    // Running integral
    // Note: phi_i(NUM_PHYSICAL_LEV) is ASSUMED to be init-ed outside this function
    Real integration = phi_i(LAST_LEV_P)[last_interface_vec_idx];

    for (int ilev = NUM_LEV - 1; ilev >= 0; --ilev) {
      // In all but the last level pack, the loop is over the whole vector
      const int vec_start = (ilev == (NUM_LEV - 1) ? last_midpoint_vec_idx
                                                   : VECTOR_SIZE - 1);

      // Precompute this product as a SIMD operation
      auto temp = PhysicalConstants::Rgas*vtheta_dp(ilev)*exner(ilev)/pnh(ilev);

      // Integrate
      auto& phi_i_val = phi_i(ilev);
      phi_i_val[vec_start] = integration;
      for (int iv = vec_start - 1; iv >= 0; --iv)
        phi_i_val[iv] = phi_i_val[iv + 1] + temp[iv + 1];

      // Update running integral
      integration = phi_i_val[0] + temp[0];
    }
  }

  template<typename ExecSpaceType>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<OnGpu<ExecSpaceType>::value>::type
  compute_phi_i_impl (const KernelVariables& kv,
                      const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& vtheta_dp,
                      const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& exner,
                      const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& pnh,
                      const ExecViewUnmanaged<      Scalar [NUM_LEV_P]>& phi_i)
  {
    // On GPU we rely on the fact that Scalar is basically double[1].
    static_assert (!OnGpu<ExecSpaceType>::value || NUM_LEV==NUM_PHYSICAL_LEV, "Error! In a GPU build we expect VECTOR_SIZE=1.\n");

    // Note: phi_i(NUM_PHYSICAL_LEV) is ASSUMED to be init-ed outside this function

    // accumulate Rgas*vtheta_dp*exner/pnh in [NUM_PHYSICAL_LEV,1].
    Dispatch<ExecSpaceType>::parallel_scan(kv.team, NUM_INTERFACE_LEV-1,
                                          [&](const int k, Real& accumulator, const bool last) {
      // level must range in [NUM_PHYSICAL_LEV-1,1], while k ranges in [0, NUM_PHYSICAL_LEV-2].
      const int level = NUM_INTERFACE_LEV-k-1;

      accumulator += PhysicalConstants::Rgas*vtheta_dp(level)[0]*exner(level)[0]/pnh(level)[0]; 

      if (last) {
        phi_i(level) = accumulator;
      }
    });
  }

  // template<typename ExecSpaceType>
  // KOKKOS_INLINE_FUNCTION
  // typename std::enable_if<!OnGpu<ExecSpaceType>::value>::type
  // compute_phinh_i_impl (const KernelVariables& [> kv <], const Real phis,
  //                       const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& vtheta_dp,
  //                       const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& p,
  //                       const ExecViewUnmanaged<      Scalar [NUM_LEV_P]>& phinh_i) const {
  //   // Running integral
  //   Real integration = phinh_i(LAST_LEV_P)[last_interface_vec_idx] = phis;

  //   for (int ilev = NUM_LEV - 1; ilev >= 0; --ilev) {
  //     // In all but the last level pack, the loop is over the whole vector
  //     const int vec_start = (ilev == (NUM_LEV - 1) ? last_midpoint_idx
  //                                                  : VECTOR_SIZE - 1);

  //     // Precompute this product as a SIMD operation
  //     auto temp = p(ilev);
  //     temp /= PhysicalConstants::p0;
  //     pow_update(temp,PhysicalConstants::kappa-1.0);
  //     temp *= vtheta_dp(ilev);
  //     temp *= PhysicalConstants::Rgas/PhysicalConstants::p0;

  //     // Integrate
  //     auto& phinh_i_val = phinh_i(ilev);
  //     phinh_i_val[vec_start] = integration;
  //     for (int iv = vec_start - 1; iv >= 0; --iv)
  //       phinh_i_val[iv] = phinh_i_val[iv + 1] + temp[iv + 1];

  //     // Update running integral
  //     integration = phinh_i_val[0] + temp[0];
  //   }
  // }

  // template<typename ExecSpaceType>
  // KOKKOS_INLINE_FUNCTION
  // typename std::enable_if<OnGpu<ExecSpaceType>::value>::type
  // compute_phinh_i_impl (const KernelVariables& kv, const Real phis,
  //                       const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& vtheta_dp,
  //                       const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& p,
  //                       const ExecViewUnmanaged<      Scalar [NUM_LEV_P]>& phinh_i) const {

  //   // On GPU we rely on the fact that Scalar is basically double[1].
  //   static_assert (!OnGpu<ExecSpaceType>::value || NUM_LEV==NUM_PHYSICAL_LEV, "Error! In a GPU build we expect VECTOR_SIZE=1.\n");

  //   // Initialize phinh_i at the surface
  //   phinh_i(NUM_INTERFACE_LEV-1) = phis;

  //   // accumulate Rgas*vtheta_dp*(p/p0)^(k-1)/p0 in [NUM_PHYSICAL_LEV,1].
  //   Dispatch<ExecSpaceType>::parallel_scan(kv.team, NUM_INTERFACE_LEV-1,
  //                                         [&](const int k, Real& accumulator, const bool last) {
  //     // level must range in [NUM_PHYSICAL_LEV-1,1], while k ranges in [0, NUM_PHYSICAL_LEV-2].
  //     const int level = NUM_INTERFACE_LEV-k-1;

  //     accumulator += PhysicalConstants::Rgas*vtheta_dp(level)[0] *
  //                    std::pow(p(level)[0]/PhysicalConstants::p0,PhysicalConstants::kappa-1.0) /
  //                    PhysicalConstants::p0; 

  //     if (last) {
  //       phinh_i(level) = accumulator;
  //     }
  //   });
  // }

};

} // namespace Homme

#endif // HOMMEXX_ELEMENT_OPS_HPP
