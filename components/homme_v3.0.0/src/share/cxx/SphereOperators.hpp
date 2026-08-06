/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_SPHERE_OPERATORS_HPP
#define HOMMEXX_SPHERE_OPERATORS_HPP

#include "Types.hpp"
#include "ElementsGeometry.hpp"
#include "CombineOps.hpp"
#include "ReferenceElement.hpp"
#include "Dimensions.hpp"
#include "KernelVariables.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/ViewUtils.hpp"

#include <Kokkos_Core.hpp>

namespace Homme {

class SphereOperators
{
  static constexpr int MAX_NUM_LEV = NUM_LEV_P;
  static constexpr int NUM_2D_VECTOR_BUFFERS = 2;
  static constexpr int NUM_3D_SCALAR_BUFFERS = 3;
  static constexpr int NUM_3D_VECTOR_BUFFERS = 3;

  // These two short names will be used to extract subviews from the 3d buffers,
  // since we can no longer use the auto keyword. Before we used to do
  //   const auto& gv = Homme::subview(vector_buf_ml,kv.team_idx,0);
  // However, now vector_buf_ml has last timension NUM_LEV_P, but we may
  // need gv to be smaller. Hence, we simply grab the pointer from the subview,
  // but then we need to explicitly tell the compiler the type of the result,
  // which can no longer be deduced. Like this:
  //   vector_buf<NUM_LEV gv(Homme::subview(vector_buf_ml,kv.team_idx,0).data());

  template<int NL>
  using scalar_buf = ExecViewUnmanaged<Scalar[NP][NP][NL]>;

  template<int NL>
  using vector_buf = ExecViewUnmanaged<Scalar[2][NP][NP][NL]>;
  
  // std::min is constexpr only from c++14 on.
  template<int M, int N>
  struct Min {
    static constexpr int value = M > N ? N : M;
  };

  template<int NUM_LEVELS>
  using DefaultProvider = ExecViewUnmanaged<const Scalar [NP][NP][NUM_LEVELS]>;
public:


  SphereOperators () = default;

  SphereOperators (const ElementsGeometry& geometry,
                   const ReferenceElement& ref_FE)
  {
    setup (geometry, ref_FE);
  }

  //only for unit tests
  SphereOperators (const Real scale_factor, const Real laplacian_rigid_factor)
  {
    m_scale_factor_inv = 1/scale_factor;
    m_laplacian_rigid_factor = laplacian_rigid_factor;
  }

  void setup (const ElementsGeometry& geometry,
              const ReferenceElement& ref_FE) {
    // Get reference element stuff
    dvv = ref_FE.get_deriv();
    m_mp = ref_FE.get_mass();

    // Get all needed 2d fields from elements
    m_d        = geometry.m_d;
    m_dinv     = geometry.m_dinv;
    m_metdet   = geometry.m_metdet;
    m_metinv   = geometry.m_metinv;
    m_spheremp = geometry.m_spheremp;
    m_scale_factor_inv = 1/geometry.m_scale_factor;
    m_laplacian_rigid_factor = geometry.m_laplacian_rigid_factor;
  }

  template<typename... Tags>
  void allocate_buffers (const Kokkos::TeamPolicy<ExecSpace,Tags...>& team_policy)
  {
    const int num_parallel_iterations = team_policy.league_size();
    const int alloc_dim = OnGpu<ExecSpace>::value ?
                          num_parallel_iterations : std::min(get_num_concurrent_teams(team_policy),num_parallel_iterations);

    if (vector_buf_ml.extent_int(0)<alloc_dim) {
      vector_buf_sl = decltype(vector_buf_sl)("",alloc_dim);
      scalar_buf_ml = decltype(scalar_buf_ml)("",alloc_dim);
      vector_buf_ml = decltype(vector_buf_ml)("",alloc_dim);
    }
  }

  void allocate_buffers (const TeamUtils<ExecSpace>& tu)
  {
    const int alloc_dim = tu.get_num_ws_slots();

    if (vector_buf_ml.extent_int(0)<alloc_dim) {
      vector_buf_sl = decltype(vector_buf_sl)("",alloc_dim);
      scalar_buf_ml = decltype(scalar_buf_ml)("",alloc_dim);
      vector_buf_ml = decltype(vector_buf_ml)("",alloc_dim);
    }
  }

  // This one is used in the unit tests
  void set_views (const ExecViewManaged<const Real         [NP][NP]>  dvv_in,
                  const ExecViewManaged<const Real * [2][2][NP][NP]>  d,
                  const ExecViewManaged<const Real * [2][2][NP][NP]>  dinv,
                  const ExecViewManaged<const Real * [2][2][NP][NP]>  metinv,
                  const ExecViewManaged<const Real *       [NP][NP]>  metdet,
                  const ExecViewManaged<const Real *       [NP][NP]>  spheremp,
                  const ExecViewManaged<const Real         [NP][NP]>  mp)
  {
    dvv = dvv_in;
    m_d = d;
    m_dinv = dinv;
    m_metinv = metinv;
    m_metdet = metdet;
    m_spheremp = spheremp;
    m_mp = mp;
  }

  // ================ SINGLE-LEVEL IMPLEMENTATION =========================== //

  KOKKOS_INLINE_FUNCTION void
  gradient_sphere_sl (const KernelVariables &kv,
                      const ExecViewUnmanaged<const Real    [NP][NP]>& scalar,
                      const ExecViewUnmanaged<      Real [2][NP][NP]>& grad_s) const
  {
    // Make sure the buffers have been created
    assert (vector_buf_sl.size()>0);

    const auto& D_inv = Homme::subview(m_dinv,kv.ie);
    const auto& temp_v_buf = Homme::subview(vector_buf_sl,kv.team_idx,0);
    constexpr int np_squared = NP * NP;
    // TODO: Use scratch space for this
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int j = loop_idx / NP;
      const int l = loop_idx % NP;
      Real dsdx(0), dsdy(0);
      for (int i = 0; i < NP; ++i) {
        dsdx += dvv(l, i) * scalar(j, i);
        dsdy += dvv(l, i) * scalar(i, j);
      }
      temp_v_buf(0, j, l) = dsdx * m_scale_factor_inv;
      temp_v_buf(1, l, j) = dsdy * m_scale_factor_inv;
    });
    kv.team_barrier();

    constexpr int grad_iters = 2 * NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, grad_iters),
                         [&](const int loop_idx) {
      const int h = (loop_idx / NP) / NP;
      const int i = (loop_idx / NP) % NP;
      const int j = loop_idx % NP;
      grad_s(h, j, i) = D_inv(h, 0, j, i) * temp_v_buf(0, j, i) +
                        D_inv(h, 1, j, i) * temp_v_buf(1, j, i);
    });
    kv.team_barrier();
  }

  KOKKOS_INLINE_FUNCTION void
  gradient_sphere_update_sl (const KernelVariables &kv,
                             const ExecViewUnmanaged<const Real    [NP][NP]>& scalar,
                             const ExecViewUnmanaged<      Real [2][NP][NP]>& grad_s) const
  {
    // Make sure the buffers have been created
    assert (vector_buf_sl.size()>0);

    constexpr int np_squared = NP * NP;
    const auto& D_inv = Homme::subview(m_dinv,kv.ie);
    const auto& temp_v_buf = Homme::subview(vector_buf_sl,kv.team_idx,0);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int j = loop_idx / NP;
      const int l = loop_idx % NP;
      Real dsdx(0), dsdy(0);
      for (int i = 0; i < NP; ++i) {
        dsdx += dvv(l, i) * scalar(j, i);
        dsdy += dvv(l, i) * scalar(i, j);
      }
      temp_v_buf(0, j, l) = dsdx * m_scale_factor_inv;
      temp_v_buf(1, l, j) = dsdy * m_scale_factor_inv;
    });
    kv.team_barrier();

    constexpr int grad_iters = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, grad_iters),
                         [&](const int loop_idx) {
      const int i = loop_idx / NP;
      const int j = loop_idx % NP;
      const auto& tmp0 = temp_v_buf(0,i,j);
      const auto& tmp1 = temp_v_buf(1,i,j);
      grad_s(0,i,j) += D_inv(0,0,i,j) * tmp0 + D_inv(0,1,i,j) * tmp1;
      grad_s(1,i,j) += D_inv(1,0,i,j) * tmp0 + D_inv(1,1,i,j) * tmp1;
    });
    kv.team_barrier();
  }

  KOKKOS_INLINE_FUNCTION void
  divergence_sphere_sl (const KernelVariables &kv,
                        const ExecViewUnmanaged<const Real [2][NP][NP]>& v,
                        const ExecViewUnmanaged<      Real    [NP][NP]>& div_v) const
  {
    // Make sure the buffers have been created
    assert (vector_buf_sl.size()>0);

    const auto& metdet = Homme::subview(m_metdet,kv.ie);
    const auto& D_inv = Homme::subview(m_dinv,kv.ie);
    const auto& gv_buf = Homme::subview(vector_buf_sl,kv.team_idx,0);
    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      const auto& v0 = v(0,igp,jgp);
      const auto& v1 = v(1,igp,jgp);
      gv_buf(0,igp,jgp) = (D_inv(0,0,igp,jgp) * v0 + D_inv(1,0,igp,jgp) * v1) * metdet(igp,jgp);
      gv_buf(1,igp,jgp) = (D_inv(0,1,igp,jgp) * v0 + D_inv(1,1,igp,jgp) * v1) * metdet(igp,jgp);
    });
    kv.team_barrier();

    constexpr int div_iters = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, div_iters),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Real dudx = 0.0, dvdy = 0.0;
      for (int kgp = 0; kgp < NP; ++kgp) {
        dudx += dvv(jgp, kgp) * gv_buf(0, igp, kgp);
        dvdy += dvv(igp, kgp) * gv_buf(1, kgp, jgp);
      }
      div_v(igp,jgp) = (dudx + dvdy) * ((1.0 / metdet(igp,jgp)) *
                                         m_scale_factor_inv);
    });
    kv.team_barrier();
  }

  KOKKOS_INLINE_FUNCTION void
  divergence_sphere_wk_sl (const KernelVariables &kv,
                           const ExecViewUnmanaged<const Real [2][NP][NP]>& v,
                           const ExecViewUnmanaged<      Real    [NP][NP]>& div_v) const
  {
    // Make sure the buffers have been created
    assert (vector_buf_sl.size()>0);

    const auto& D_inv = Homme::subview(m_dinv,kv.ie);
    const auto& spheremp = Homme::subview(m_spheremp,kv.ie);
    const auto& gv_buf = Homme::subview(vector_buf_sl,kv.team_idx,0);

    // copied from strong divergence as is but without metdet
    // conversion to contravariant
    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      const auto& v0 = v(0,igp,jgp);
      const auto& v1 = v(1,igp,jgp);
      gv_buf(0,igp,jgp) = D_inv(0,0,igp,jgp) * v0 + D_inv(1,0,igp,jgp) * v1;
      gv_buf(1,igp,jgp) = D_inv(0,1,igp,jgp) * v0 + D_inv(1,1,igp,jgp) * v1;
    });
    kv.team_barrier();

    // in strong div
    // kgp = i in strong code, jgp=j, igp=l
    // in weak div, n is like j in strong div,
    // n(weak)=j(strong)=jgp
    // m(weak)=l(strong)=igp
    // j(weak)=i(strong)=kgp
    constexpr int div_iters = NP * NP;
    // keeping indices' names as in F
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, div_iters),
                         [&](const int loop_idx) {
      // Note: for this one time, it is better if m strides faster, due to
      //       the way the views are accessed.
      const int mgp = loop_idx % NP;
      const int ngp = loop_idx / NP;
      Real dd = 0.0;
      for (int jgp = 0; jgp < NP; ++jgp) {
        dd -= (spheremp(ngp, jgp) * gv_buf(0, ngp, jgp) * dvv(jgp, mgp) +
               spheremp(jgp, mgp) * gv_buf(1, jgp, mgp) * dvv(jgp, ngp)) *
              m_scale_factor_inv;
      }
      div_v(ngp, mgp) = dd;
    });
    kv.team_barrier();

  } // end of divergence_sphere_wk_sl

  // Note that divergence_sphere requires scratch space of 3 x NP x NP Reals
  // This must be called from the device space
  KOKKOS_INLINE_FUNCTION void
  vorticity_sphere_sl (const KernelVariables &kv,
                       const ExecViewUnmanaged<const Real [NP][NP]>& u,
                       const ExecViewUnmanaged<const Real [NP][NP]>& v,
                       const ExecViewUnmanaged<      Real [NP][NP]>& vort) const
  {
    // Make sure the buffers have been created
    assert (vector_buf_sl.size()>0);

    const auto& D = Homme::subview(m_d,kv.ie);
    const auto& metdet = Homme::subview(m_metdet,kv.ie);
    const auto& vcov_buf = Homme::subview(vector_buf_sl,kv.team_idx,0);

    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      const auto& u_ij = u(igp,jgp);
      const auto& v_ij = v(igp,jgp);
      vcov_buf(0,igp,jgp) = D(0,0,igp,jgp) * u_ij + D(0,1,igp,jgp) * v_ij;
      vcov_buf(1,igp,jgp) = D(1,0,igp,jgp) * u_ij + D(1,1,igp,jgp) * v_ij;
    });
    kv.team_barrier();

    constexpr int vort_iters = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, vort_iters),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Real dudy = 0.0;
      Real dvdx = 0.0;
      for (int kgp = 0; kgp < NP; ++kgp) {
        dvdx += dvv(jgp, kgp) * vcov_buf(1, igp, kgp);
        dudy += dvv(igp, kgp) * vcov_buf(0, kgp, jgp);
      }

      vort(igp, jgp) = (dvdx - dudy) * ((1.0 / metdet(igp, jgp)) *
                                        m_scale_factor_inv);
    });
    kv.team_barrier();
  }

  // analog of fortran's laplace_wk_sphere
  // Single level implementation
  KOKKOS_INLINE_FUNCTION void
  laplace_wk_sl (const KernelVariables &kv,
                 const ExecViewUnmanaged<const Real [NP][NP]>& field,
                 const ExecViewUnmanaged<      Real [NP][NP]>& laplace) const
  {
    // Make sure the buffers have been created
    assert (vector_buf_sl.size()>0);

    const auto& grad_s = Homme::subview(vector_buf_sl, kv.team_idx, 1);
    gradient_sphere_sl(kv, field, grad_s);
    divergence_sphere_wk_sl(kv, grad_s, laplace);
  } // end of laplace_wk_sl

  // ================ MULTI-LEVEL IMPLEMENTATION =========================== //

  // Note: if you are puzzled by the use/need of ViewConst, don't worry, you're not alone.
  //       To be clear, using `const typename ExecViewUnmanaged<Scalar[...]>::const_type`
  //       would also work. I prefer ViewConst cause it does not change the 'template
  //       signature of the type', while Kokkos::View<...>::const_type has a pre-defined
  //       templated signature (with data type, layout and memory space).
  //       Anyhow, back to the point, the key fact is that (quoting cppreference),
  //
  //        "The nested-name-specifier (everything to the left of the scope resolution
  //         operator ::) of a type that was specified using a qualified-id" [do not
  //         participate in template argument deduction].
  //
  //       So, without the `typename Blah::type` trick, the type of scalar is used to deduce
  //       the template arguments, but if the input view had data type non-const, this
  //       makes the template deduction fail. On the other hand, *with* the `typename Blah::type`
  //       trick, only the type of the output view is used to deduce the template args. This
  //       succeeds, and then the copy constructor of View is used to produce a View<const T>
  //       from a View<T>. Magic.

  template<int NUM_LEV_OUT, typename InputProvider>
  KOKKOS_INLINE_FUNCTION void
  gradient_sphere (const KernelVariables &kv,
                   const InputProvider& scalar,
                   const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_OUT]>& grad_s,
                   const int NUM_LEV_REQUEST) const
  {
    assert(NUM_LEV_REQUEST>=0);
    assert(NUM_LEV_REQUEST<=NUM_LEV_OUT);

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    const auto& D_inv = Homme::subview(m_dinv, kv.ie);

    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        Scalar v0, v1;
        for (int kgp = 0; kgp < NP; ++kgp) {
          v0 += dvv(jgp, kgp) * scalar(igp, kgp, ilev);
          v1 += dvv(igp, kgp) * scalar(kgp, jgp, ilev);
        }
        v0 *= m_scale_factor_inv;
        v1 *= m_scale_factor_inv;
        grad_s(0,igp,jgp,ilev) = D_inv(0,0,igp,jgp) * v0 + D_inv(0,1,igp,jgp) * v1;
        grad_s(1,igp,jgp,ilev) = D_inv(1,0,igp,jgp) * v0 + D_inv(1,1,igp,jgp) * v1;
      });
    });
    kv.team_barrier();
  }

  template<int NUM_LEV_OUT, typename InputProvider, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  gradient_sphere (const KernelVariables &kv,
                   const InputProvider& scalar,
                   const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_OUT]>& grad_s) const
  {
    static_assert(NUM_LEV_REQUEST>=0, "Error! Invalid value for NUM_LEV_REQUEST.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_OUT, "Error! Output view does not have enough levels.\n");

    gradient_sphere<NUM_LEV_OUT, InputProvider>(kv, scalar, grad_s, NUM_LEV_REQUEST);
  }

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  gradient_sphere_update (const KernelVariables &kv,
                          const typename ViewConst<ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_IN]>>::type& scalar,
                          const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_OUT]>& grad_s) const
  {
    static_assert(NUM_LEV_REQUEST>=0, "Error! Invalid value for NUM_LEV_REQUEST.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_IN, "Error! Input view does not have enough levels.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_OUT, "Error! Output view does not have enough levels.\n");

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    const auto& D_inv = Homme::subview(m_dinv, kv.ie);
    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        Scalar dsdx, dsdy;
        for (int kgp = 0; kgp < NP; ++kgp) {
          dsdx += dvv(jgp, kgp) * scalar(igp, kgp, ilev);
          dsdy += dvv(igp, kgp) * scalar(kgp, jgp, ilev);
        }
        dsdx *= m_scale_factor_inv;
        dsdy *= m_scale_factor_inv;
        grad_s(0,igp,jgp,ilev) += D_inv(0,0,igp,jgp) * dsdx + D_inv(0,1,igp,jgp) * dsdy;
        grad_s(1,igp,jgp,ilev) += D_inv(1,0,igp,jgp) * dsdx + D_inv(1,1,igp,jgp) * dsdy;
      });
    });
    kv.team_barrier();
  }

  template<int NUM_LEV_OUT, typename InputProvider, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  divergence_sphere (const KernelVariables &kv,
                     const InputProvider& v,
                     const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_OUT]>& div_v,
                     const Real alpha = 1.0, const Real beta = 0.0) const
  {
    divergence_sphere_cm<CombineMode::Replace,
                         InputProvider,
                         NUM_LEV_OUT,
                         NUM_LEV_REQUEST>(kv,v,div_v,alpha,beta);

  }

  template<int NUM_LEV_OUT, typename InputProvider>
  KOKKOS_INLINE_FUNCTION void
  divergence_sphere_nlev (const KernelVariables &kv,
                          const InputProvider& v,
                          const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_OUT]>& div_v,
                          const int NUM_LEV_REQUEST,
                          const Real alpha = 1.0, const Real beta = 0.0) const
  {
    divergence_sphere_cm<CombineMode::Replace,
                         InputProvider,
                         NUM_LEV_OUT>(kv,v,div_v,alpha,beta,NUM_LEV_REQUEST);

  }

  template<CombineMode CM, typename InputProvider, int NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  divergence_sphere_cm (const KernelVariables &kv,
                        const InputProvider& v,
                        const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_OUT]>& div_v,
                        const Real alpha, const Real beta,
                        const int NUM_LEV_REQUEST) const
  {
    assert(NUM_LEV_REQUEST>=0);
    assert(NUM_LEV_REQUEST<=NUM_LEV_OUT);

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    const auto& D_inv = Homme::subview(m_dinv, kv.ie);
    const auto& metdet = Homme::subview(m_metdet, kv.ie);
    vector_buf<NUM_LEV_OUT> gv_buf(Homme::subview(vector_buf_ml,kv.team_idx, 0).data());
    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        const auto& v0 = v(0, igp, jgp, ilev);
        const auto& v1 = v(1, igp, jgp, ilev);
        gv_buf(0,igp,jgp,ilev) = (D_inv(0,0,igp,jgp) * v0 + D_inv(1,0,igp,jgp) * v1) * metdet(igp,jgp);
        gv_buf(1,igp,jgp,ilev) = (D_inv(0,1,igp,jgp) * v0 + D_inv(1,1,igp,jgp) * v1) * metdet(igp,jgp);
      });
    });
    kv.team_barrier();

    // j, l, i -> i, j, k
    constexpr int div_iters = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, div_iters),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        Scalar dudx, dvdy;
        for (int kgp = 0; kgp < NP; ++kgp) {
          dudx += dvv(jgp, kgp) * gv_buf(0, igp, kgp, ilev);
          dvdy += dvv(igp, kgp) * gv_buf(1, kgp, jgp, ilev);
        }
        combine<CM>((dudx + dvdy) * (1.0 / metdet(igp, jgp) * m_scale_factor_inv),
                     div_v(igp, jgp, ilev), alpha, beta);
      });
    });
    kv.team_barrier();
  }

  template<CombineMode CM, typename InputProvider,
           int NUM_LEV_OUT, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  divergence_sphere_cm (const KernelVariables &kv,
                        const InputProvider& v,
                        const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_OUT]>& div_v,
                        const Real alpha = 1.0, const Real beta = 0.0) const
  {
    static_assert(NUM_LEV_REQUEST>=0, "Error! Invalid value for NUM_LEV_REQUEST.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_OUT, "Error! Output view does not have enough levels.\n");

    divergence_sphere_cm<CM, InputProvider, NUM_LEV_OUT>(kv, v, div_v, alpha, beta, NUM_LEV_REQUEST);
  }

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  divergence_sphere_update (const KernelVariables &kv,
                            const Real alpha, const bool add_hyperviscosity,
                            const typename ViewConst<ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_IN]>>::type& vstar,
                            const typename ViewConst<ExecViewUnmanaged<Scalar    [NP][NP][NUM_LEV_IN]>>::type& qdp,
                            // On input, qtens_biharmonic if add_hyperviscosity, undefined
                            // if not; on output, qtens.
                            const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_OUT]>& qtens) const
  {
    static_assert(NUM_LEV_REQUEST>=0, "Error! Invalid value for NUM_LEV_REQUEST.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_IN, "Error! Input view does not have enough levels.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_OUT, "Error! Output view does not have enough levels.\n");

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    const auto& D_inv = Homme::subview(m_dinv, kv.ie);
    const auto& metdet = Homme::subview(m_metdet, kv.ie);
    vector_buf<NUM_LEV_REQUEST> gv(Homme::subview(vector_buf_ml,kv.team_idx,0).data());
    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        const auto& qdpijk = qdp(igp, jgp, ilev);
        const auto v0 = vstar(0, igp, jgp, ilev) * qdpijk;
        const auto v1 = vstar(1, igp, jgp, ilev) * qdpijk;
        gv(0,igp,jgp,ilev) = (D_inv(0,0,igp,jgp) * v0 + D_inv(1,0,igp,jgp) * v1) * metdet(igp,jgp);
        gv(1,igp,jgp,ilev) = (D_inv(0,1,igp,jgp) * v0 + D_inv(1,1,igp,jgp) * v1) * metdet(igp,jgp);
      });
    });
    kv.team_barrier();

    constexpr int div_iters = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, div_iters),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        Scalar dudx, dvdy;
        for (int kgp = 0; kgp < NP; ++kgp) {
          dudx += dvv(jgp, kgp) * gv(0, igp, kgp, ilev);
          dvdy += dvv(igp, kgp) * gv(1, kgp, jgp, ilev);
        }
        const Scalar qtensijk0 = add_hyperviscosity ? qtens(igp,jgp,ilev) : 0;
        qtens(igp,jgp,ilev) = (qdp(igp,jgp,ilev) +
                               alpha*((dudx + dvdy) * (1.0 / metdet(igp,jgp) * m_scale_factor_inv)) +
                               qtensijk0);
      });
    });
    kv.team_barrier();
  }

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  vorticity_sphere (const KernelVariables &kv,
                    const typename ViewConst<ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_IN]>>::type& u,
                    const typename ViewConst<ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_IN]>>::type& v,
                    const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_OUT]>& vort) const
  {
    static_assert(NUM_LEV_REQUEST>=0, "Error! Invalid value for NUM_LEV_REQUEST.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_IN, "Error! Input view does not have enough levels.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_OUT, "Error! Output view does not have enough levels.\n");

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    const auto& D = Homme::subview(m_d, kv.ie);
    const auto& metdet = Homme::subview(m_metdet, kv.ie);
    vector_buf<NUM_LEV_REQUEST> vcov_buf(Homme::subview(vector_buf_ml,kv.team_idx,0).data());
    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        const auto& u_ijk = u(igp, jgp, ilev);
        const auto& v_ijk = v(igp, jgp, ilev);
        vcov_buf(0,igp,jgp,ilev) = D(0,0,igp,jgp) * u_ijk + D(0,1,igp,jgp) * v_ijk;
        vcov_buf(1,igp,jgp,ilev) = D(1,0,igp,jgp) * u_ijk + D(1,1,igp,jgp) * v_ijk;
      });
    });
    kv.team_barrier();

    constexpr int vort_iters = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, vort_iters),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        Scalar dudy, dvdx;
        for (int kgp = 0; kgp < NP; ++kgp) {
          dvdx += dvv(jgp, kgp) * vcov_buf(1, igp, kgp, ilev);
          dudy += dvv(igp, kgp) * vcov_buf(0, kgp, jgp, ilev);
        }
        vort(igp, jgp, ilev) = (dvdx - dudy) * (1.0 / metdet(igp, jgp) *
                                                m_scale_factor_inv);
      });
    });
    kv.team_barrier();
  }

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  vorticity_sphere (const KernelVariables &kv,
                    const typename ViewConst<ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_IN]>>::type& v,
                    const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_OUT]>& vort,
                    const int NUM_LEV_REQUEST) const
  {
    assert(NUM_LEV_REQUEST>=0);
    assert(NUM_LEV_REQUEST<=NUM_LEV_IN);
    assert(NUM_LEV_REQUEST<=NUM_LEV_OUT);

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    const auto& D = Homme::subview(m_d, kv.ie);
    const auto& metdet = Homme::subview(m_metdet, kv.ie);
    vector_buf<NUM_LEV_OUT> sphere_buf(Homme::subview(vector_buf_ml,kv.team_idx,0).data());
    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        const auto& v0 = v(0,igp,jgp,ilev);
        const auto& v1 = v(1,igp,jgp,ilev);
        sphere_buf(0,igp,jgp,ilev) = D(0,0,igp,jgp) * v0 + D(0,1,igp,jgp) * v1;
        sphere_buf(1,igp,jgp,ilev) = D(1,0,igp,jgp) * v0 + D(1,1,igp,jgp) * v1;
      });
    });
    kv.team_barrier();

    constexpr int vort_iters = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, vort_iters),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        Scalar dudy, dvdx;
        for (int kgp = 0; kgp < NP; ++kgp) {
          dvdx += dvv(jgp, kgp) * sphere_buf(1, igp, kgp, ilev);
          dudy += dvv(igp, kgp) * sphere_buf(0, kgp, jgp, ilev);
        }
        vort(igp, jgp, ilev) = (dvdx - dudy) * (1.0 / metdet(igp, jgp) *
                                                m_scale_factor_inv);
      });
    });
    kv.team_barrier();
  }

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  vorticity_sphere (const KernelVariables &kv,
                    const typename ViewConst<ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_IN]>>::type& v,
                    const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_OUT]>& vort) const
  {
    static_assert(NUM_LEV_REQUEST>=0, "Error! Invalid value for NUM_LEV_REQUEST.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_IN, "Error! Input view does not have enough levels.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_OUT, "Error! Output view does not have enough levels.\n");

    vorticity_sphere<NUM_LEV_OUT,NUM_LEV_IN>(kv, v, vort, NUM_LEV_REQUEST);
  }

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  divergence_sphere_wk (const KernelVariables &kv,
                        // On input, a field whose divergence is sought; on
                        // output, the view's data are invalid.
                        const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_IN]>& v,
                        const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_OUT]>& div_v,
                        const int NUM_LEV_REQUEST) const
  {
    assert(NUM_LEV_REQUEST>=0);
    assert(NUM_LEV_REQUEST<=NUM_LEV_IN);
    assert(NUM_LEV_REQUEST<=NUM_LEV_OUT);

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    const auto& D_inv = Homme::subview(m_dinv, kv.ie);
    const auto& spheremp = Homme::subview(m_spheremp, kv.ie);
    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        const auto v0 = v(0,igp,jgp,ilev);
        const auto v1 = v(1,igp,jgp,ilev);
        v(0,igp,jgp,ilev) = D_inv(0, 0, igp, jgp) * v0 + D_inv(1, 0, igp, jgp) * v1;
        v(1,igp,jgp,ilev) = D_inv(0, 1, igp, jgp) * v0 + D_inv(1, 1, igp, jgp) * v1;
      });
    });
    kv.team_barrier();

    constexpr int div_iters = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, div_iters),
                         [&](const int loop_idx) {
      // Note: for this one time, it is better if m strides faster, due to
      //       the way the views are accessed.
      const int mgp = loop_idx % NP;
      const int ngp = loop_idx / NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        Scalar dd;
        // TODO: move multiplication by scale_factor_inv outside the loop
        for (int jgp = 0; jgp < NP; ++jgp) {
          // Here, v is the temporary buffer, aliased on the input v.
          dd -= (spheremp(ngp, jgp) * v(0, ngp, jgp, ilev) * dvv(jgp, mgp) +
                 spheremp(jgp, mgp) * v(1, jgp, mgp, ilev) * dvv(jgp, ngp)) *
                m_scale_factor_inv;
        }
        div_v(ngp, mgp, ilev) = dd;
      });
    });
    kv.team_barrier();

  }//end of divergence_sphere_wk

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  divergence_sphere_wk (const KernelVariables &kv,
                        // On input, a field whose divergence is sought; on
                        // output, the view's data are invalid.
                        const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_IN]>& v,
                        const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_OUT]>& div_v) const
  {
    divergence_sphere_wk<NUM_LEV_OUT, NUM_LEV_IN>(kv, v, div_v, NUM_LEV_REQUEST);
  }//end of divergence_sphere_wk

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  laplace_simple (const KernelVariables &kv,
                  const typename ViewConst<ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_IN]>>::type& field,
                  const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_OUT]>& laplace,
                  const int NUM_LEV_REQUEST) const
  {
    assert(NUM_LEV_REQUEST<=NUM_LEV_IN);
    assert(NUM_LEV_REQUEST<=NUM_LEV_OUT);

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    vector_buf<NUM_LEV_OUT> grad_s(Homme::subview(vector_buf_ml, kv.team_idx, 0).data());
    gradient_sphere<NUM_LEV_OUT,decltype(field)>(kv, field, grad_s, NUM_LEV_REQUEST);
    divergence_sphere_wk<NUM_LEV_OUT,NUM_LEV_OUT>(kv, grad_s, laplace, NUM_LEV_REQUEST);
  }//end of laplace_simple

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  laplace_simple (const KernelVariables &kv,
                  const typename ViewConst<ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_IN]>>::type& field,
                  const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_OUT]>& laplace) const
  {
    static_assert(NUM_LEV_REQUEST>=0, "Error! Invalid value for NUM_LEV_REQUEST.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_IN, "Error! Input view does not have enough levels.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_OUT, "Error! Output view does not have enough levels.\n");

    laplace_simple<NUM_LEV_OUT,NUM_LEV_IN>(kv, field, laplace, NUM_LEV_REQUEST);
  }//end of laplace_simple

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  laplace_tensor(const KernelVariables &kv,
                 const ExecViewUnmanaged<const Real   [2][2][NP][NP]>&              tensorVisc,
                 const typename ViewConst<ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_IN]>>::type&  field,         // input
                 const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_OUT]>& laplace) const
  {
    static_assert(NUM_LEV_REQUEST>=0, "Error! Invalid value for NUM_LEV_REQUEST.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_IN, "Error! Input view does not have enough levels.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_OUT, "Error! Output view does not have enough levels.\n");

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    vector_buf<NUM_LEV_REQUEST> grad_s(Homme::subview(vector_buf_ml, kv.team_idx, 1).data());
    vector_buf<NUM_LEV_REQUEST> sphere_buf(Homme::subview(vector_buf_ml, kv.team_idx, 2).data());

    gradient_sphere<NUM_LEV_REQUEST,decltype(field),NUM_LEV_REQUEST>(kv, field, grad_s);
    //now multiply tensorVisc(:,:,i,j)*grad_s(i,j) (matrix*vector, independent of i,j )
    //but it requires a temp var to store a result. the result is then placed to grad_s,
    //or should it be an extra temp var instead of an extra loop?
    constexpr int num_iters = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, num_iters),
                       [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        const auto& grad_s0 = grad_s(0,igp,jgp,ilev);
        const auto& grad_s1 = grad_s(1,igp,jgp,ilev);
        sphere_buf(0,igp,jgp,ilev) = tensorVisc(0,0,igp,jgp) * grad_s0 + tensorVisc(1,0,igp,jgp) * grad_s1;
        sphere_buf(1,igp,jgp,ilev) = tensorVisc(0,1,igp,jgp) * grad_s0 + tensorVisc(1,1,igp,jgp) * grad_s1;
      });
    });
    kv.team_barrier();

    divergence_sphere_wk<NUM_LEV_OUT,NUM_LEV_REQUEST,NUM_LEV_REQUEST>(kv, sphere_buf, laplace);
  }//end of laplace_tensor

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  curl_sphere_wk_testcov (const KernelVariables &kv,
                          const typename ViewConst<ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_IN]>>::type& scalar,
                          const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_OUT]>& curls,
                          const int NUM_LEV_REQUEST) const
  {
    assert(NUM_LEV_REQUEST>=0);
    assert(NUM_LEV_REQUEST<=NUM_LEV_IN);
    assert(NUM_LEV_REQUEST<=NUM_LEV_OUT);

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    const auto& D = Homme::subview(m_d, kv.ie);
    vector_buf<NUM_LEV_OUT> sphere_buf(Homme::subview(vector_buf_ml,kv.team_idx,0).data());
    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared), [&](const int loop_idx) {
      const int ngp = loop_idx / NP;
      const int mgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        auto& sb0 = sphere_buf(0, ngp, mgp, ilev);
        auto& sb1 = sphere_buf(1, ngp, mgp, ilev);
        sb0 = 0;
        sb1 = 0;
        for (int jgp = 0; jgp < NP; ++jgp) {
          sb0 -= m_mp(jgp,mgp)*scalar(jgp,mgp,ilev)*dvv(jgp,ngp);
          sb1 += m_mp(ngp,jgp)*scalar(ngp,jgp,ilev)*dvv(jgp,mgp);
        }
      });
    });
    kv.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP; //slowest
      const int jgp = loop_idx % NP; //fastest
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        const auto& sb0 = sphere_buf(0, igp, jgp, ilev);
        const auto& sb1 = sphere_buf(1, igp, jgp, ilev);
        curls(0,igp,jgp,ilev) = (D(0,0,igp,jgp) * sb0 + D(1,0,igp,jgp) * sb1) * m_scale_factor_inv;
        curls(1,igp,jgp,ilev) = (D(0,1,igp,jgp) * sb0 + D(1,1,igp,jgp) * sb1) * m_scale_factor_inv;
      });
    });
    kv.team_barrier();
  }

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  curl_sphere_wk_testcov (const KernelVariables &kv,
                          const typename ViewConst<ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_IN]>>::type& scalar,
                          const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_OUT]>& curls) const
  {
    static_assert(NUM_LEV_REQUEST>=0, "Error! Invalid value for NUM_LEV_REQUEST.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_IN, "Error! Input view does not have enough levels.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_OUT, "Error! Output view does not have enough levels.\n");

    curl_sphere_wk_testcov<NUM_LEV_OUT,NUM_LEV_IN>(kv, scalar, curls, NUM_LEV_REQUEST);
  }

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  curl_sphere_wk_testcov_update (const KernelVariables &kv, const Real alpha, const Real beta,
                                 const typename ViewConst<ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_IN]>>::type& scalar,
                                 const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_OUT]>& curls,
                                 const int NUM_LEV_REQUEST) const
  {
    assert(NUM_LEV_REQUEST>=0);
    assert(NUM_LEV_REQUEST<=NUM_LEV_IN);
    assert(NUM_LEV_REQUEST<=NUM_LEV_OUT);

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    const auto& D = Homme::subview(m_d, kv.ie);
    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared), [&](const int loop_idx) {
      const int ngp = loop_idx / NP;
      const int mgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        Scalar sb0, sb1;
        for (int jgp = 0; jgp < NP; ++jgp) {
          sb0 -= m_mp(jgp,mgp)*scalar(jgp,mgp,ilev)*dvv(jgp,ngp);
          sb1 += m_mp(ngp,jgp)*scalar(ngp,jgp,ilev)*dvv(jgp,mgp);
        }
        curls(0,ngp,mgp,ilev) = beta*curls(0,ngp,mgp,ilev) + alpha *
                                ( D(0,0,ngp,mgp)*sb0 + D(1,0,ngp,mgp)*sb1 )
                                * m_scale_factor_inv;
        curls(1,ngp,mgp,ilev) = beta*curls(1,ngp,mgp,ilev) + alpha *
                                ( D(0,1,ngp,mgp)*sb0 + D(1,1,ngp,mgp)*sb1 )
                              * m_scale_factor_inv;
      });
    });
    kv.team_barrier();
  }

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  curl_sphere_wk_testcov_update (const KernelVariables &kv, const Real alpha, const Real beta,
                                 const typename ViewConst<ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_IN]>>::type& scalar,
                                 const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_OUT]>& curls) const
  {
    static_assert(NUM_LEV_REQUEST>=0, "Error! Invalid value for NUM_LEV_REQUEST.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_IN, "Error! Input view does not have enough levels.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_OUT, "Error! Output view does not have enough levels.\n");

    curl_sphere_wk_testcov_update<NUM_LEV_OUT,NUM_LEV_IN>(kv, alpha, beta, scalar, curls, NUM_LEV_REQUEST);
  }

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  grad_sphere_wk_testcov (const KernelVariables &kv,
                          const typename ViewConst<ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_IN]>>::type& scalar,
                          const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_OUT]>& grads,
                          const int NUM_LEV_REQUEST) const
  {
    assert(NUM_LEV_REQUEST>=0);
    assert(NUM_LEV_REQUEST<=NUM_LEV_IN);
    assert(NUM_LEV_REQUEST<=NUM_LEV_OUT);

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    const auto& D = Homme::subview(m_d, kv.ie);
    const auto& metinv = Homme::subview(m_metinv, kv.ie);
    const auto& metdet = Homme::subview(m_metdet, kv.ie);
    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared), [&](const int loop_idx) {
      const int ngp = loop_idx / NP;
      const int mgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        Scalar b0, b1;
        for (int jgp = 0; jgp < NP; ++jgp) {
          const auto& mpnj = m_mp(ngp,jgp);
          const auto& mpjm = m_mp(jgp,mgp);
          const auto& md = metdet(ngp,mgp);
          const auto& snj = scalar(ngp,jgp,ilev);
          const auto& sjm = scalar(jgp,mgp,ilev);
          const auto& djm = dvv(jgp,mgp);
          const auto& djn = dvv(jgp,ngp);
          b0 -= (mpnj * metinv(0,0,ngp,mgp) * md * snj * djm +
                 mpjm * metinv(0,1,ngp,mgp) * md * sjm * djn);
          b1 -= (mpnj * metinv(1,0,ngp,mgp) * md * snj * djm +
                 mpjm * metinv(1,1,ngp,mgp) * md * sjm * djn);
        }
        grads(0,ngp,mgp,ilev) = (D(0,0,ngp,mgp) * b0 + D(1,0,ngp,mgp) * b1) * m_scale_factor_inv;
        grads(1,ngp,mgp,ilev) = (D(0,1,ngp,mgp) * b0 + D(1,1,ngp,mgp) * b1) * m_scale_factor_inv;
      });
    });
    kv.team_barrier();
  }

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  grad_sphere_wk_testcov (const KernelVariables &kv,
                          const typename ViewConst<ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV_IN]>>::type& scalar,
                          const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_OUT]>& grads) const
  {
    static_assert(NUM_LEV_REQUEST>=0, "Error! Invalid value for NUM_LEV_REQUEST.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_IN, "Error! Input view does not have enough levels.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_OUT, "Error! Output view does not have enough levels.\n");

    grad_sphere_wk_testcov<NUM_LEV_OUT,NUM_LEV_IN>(kv, scalar, grads, NUM_LEV_REQUEST);
  }

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  vlaplace_sphere_wk_cartesian (const KernelVariables &kv,
                                const ExecViewUnmanaged<const Real [2][2][NP][NP]>&          tensorVisc,
                                const ExecViewUnmanaged<const Real [2][3][NP][NP]>&          vec_sph2cart,
                                const typename ViewConst<ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_IN]>>::type& vector,
                                const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_OUT]>& laplace) const
  {
    static_assert(NUM_LEV_REQUEST>=0, "Error! Invalid value for NUM_LEV_REQUEST.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_IN, "Error! Input view does not have enough levels.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_OUT, "Error! Output view does not have enough levels.\n");

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    const auto& spheremp = Homme::subview(m_spheremp, kv.ie);
    scalar_buf<NUM_LEV_REQUEST> laplace0(Homme::subview(scalar_buf_ml,kv.team_idx,0).data());
    scalar_buf<NUM_LEV_REQUEST> laplace1(Homme::subview(scalar_buf_ml,kv.team_idx,1).data());
    scalar_buf<NUM_LEV_REQUEST> laplace2(Homme::subview(scalar_buf_ml,kv.team_idx,2).data());
    constexpr int np_squared = NP * NP;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP; //slowest
      const int jgp = loop_idx % NP; //fastest
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        const auto& v0 = vector(0,igp,jgp,ilev);
        const auto& v1 = vector(1,igp,jgp,ilev);
        laplace0(igp,jgp,ilev) = vec_sph2cart(0,0,igp,jgp)*v0 + vec_sph2cart(1,0,igp,jgp)*v1 ;
        laplace1(igp,jgp,ilev) = vec_sph2cart(0,1,igp,jgp)*v0 + vec_sph2cart(1,1,igp,jgp)*v1 ;
        laplace2(igp,jgp,ilev) = vec_sph2cart(0,2,igp,jgp)*v0 + vec_sph2cart(1,2,igp,jgp)*v1 ;
      });
    });
    kv.team_barrier();

    // Use laplace* as input, and then overwrite it with the output (saves temporaries)
    laplace_tensor<NUM_LEV_REQUEST,NUM_LEV_REQUEST,NUM_LEV_REQUEST>(kv,tensorVisc,laplace0,laplace0);
    laplace_tensor<NUM_LEV_REQUEST,NUM_LEV_REQUEST,NUM_LEV_REQUEST>(kv,tensorVisc,laplace1,laplace1);
    laplace_tensor<NUM_LEV_REQUEST,NUM_LEV_REQUEST,NUM_LEV_REQUEST>(kv,tensorVisc,laplace2,laplace2);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP; //slowest
      const int jgp = loop_idx % NP; //fastest
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
#define UNDAMPRRCART
#ifdef UNDAMPRRCART
        laplace(0,igp,jgp,ilev) = vec_sph2cart(0,0,igp,jgp)*laplace0(igp,jgp,ilev)
                                + vec_sph2cart(0,1,igp,jgp)*laplace1(igp,jgp,ilev)
                                + vec_sph2cart(0,2,igp,jgp)*laplace2(igp,jgp,ilev)
                                + 2.0*spheremp(igp,jgp)*vector(0,igp,jgp,ilev)
                                        *(m_laplacian_rigid_factor)*(m_laplacian_rigid_factor);

        laplace(1,igp,jgp,ilev) = vec_sph2cart(1,0,igp,jgp)*laplace0(igp,jgp,ilev)
                                + vec_sph2cart(1,1,igp,jgp)*laplace1(igp,jgp,ilev)
                                + vec_sph2cart(1,2,igp,jgp)*laplace2(igp,jgp,ilev)
                                + 2.0*spheremp(igp,jgp)*vector(1,igp,jgp,ilev)
                                        *(m_laplacian_rigid_factor)*(m_laplacian_rigid_factor);
#else
        laplace(0,igp,jgp,ilev) = vec_sph2cart(0,0,igp,jgp)*laplace0(igp,jgp,ilev)
                                + vec_sph2cart(0,1,igp,jgp)*laplace1(igp,jgp,ilev)
                                + vec_sph2cart(0,2,igp,jgp)*laplace2(igp,jgp,ilev);
        laplace(1,igp,jgp,ilev) = vec_sph2cart(1,0,igp,jgp)*laplace0(igp,jgp,ilev)
                                + vec_sph2cart(1,1,igp,jgp)*laplace1(igp,jgp,ilev)
                                + vec_sph2cart(1,2,igp,jgp)*laplace2(igp,jgp,ilev);
#endif
      });
    });
    kv.team_barrier();
  } // end of vlaplace_sphere_wk_cartesian

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  vlaplace_sphere_wk_contra (const KernelVariables &kv, const Real nu_ratio,
                             const typename ViewConst<ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_IN]>>::type& vector,
                             const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_OUT]>& laplace,
                             const int NUM_LEV_REQUEST) const
  {
    assert(NUM_LEV_REQUEST>=0);
    assert(NUM_LEV_REQUEST<=NUM_LEV_IN);
    assert(NUM_LEV_REQUEST<=NUM_LEV_OUT);

    // Make sure the buffers have been created
    assert (vector_buf_ml.size()>0);

    const auto& spheremp = Homme::subview(m_spheremp, kv.ie);
    scalar_buf<NUM_LEV_OUT> div (Homme::subview(scalar_buf_ml,kv.team_idx,0).data());
    scalar_buf<NUM_LEV_OUT> vort(Homme::subview(scalar_buf_ml,kv.team_idx,0).data());
    vector_buf<NUM_LEV_OUT> grad_curl_cov(Homme::subview(vector_buf_ml,kv.team_idx,1).data());
    constexpr int np_squared = NP * NP;

    // grad(div(v))
    divergence_sphere_nlev(kv,vector,div,NUM_LEV_REQUEST);
    if (nu_ratio>0 && nu_ratio!=1.0) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                          [&](const int loop_idx) {
        const int igp = loop_idx / NP; //slow
        const int jgp = loop_idx % NP; //fast
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
          div(igp,jgp,ilev) *= nu_ratio;
        });
      });
      kv.team_barrier();
    }
    grad_sphere_wk_testcov<NUM_LEV_OUT,NUM_LEV_OUT>(kv,div,grad_curl_cov,NUM_LEV_REQUEST);

    // curl(curl(v))
    vorticity_sphere<NUM_LEV_OUT,NUM_LEV_IN>(kv,vector,vort,NUM_LEV_REQUEST);
    curl_sphere_wk_testcov_update<NUM_LEV_OUT,NUM_LEV_OUT>(kv,-1.0,1.0,vort,grad_curl_cov,NUM_LEV_REQUEST);

    const auto re2 = m_laplacian_rigid_factor*m_laplacian_rigid_factor;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, np_squared),
                        [&](const int loop_idx) {
      const int igp = loop_idx / NP; //slow
      const int jgp = loop_idx % NP; //fast
      const auto f = 2.0*spheremp(igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV_REQUEST), [&] (const int& ilev) {
        laplace(0,igp,jgp,ilev) = f*vector(0,igp,jgp,ilev)*re2 + grad_curl_cov(0,igp,jgp,ilev);
        laplace(1,igp,jgp,ilev) = f*vector(1,igp,jgp,ilev)*re2 + grad_curl_cov(1,igp,jgp,ilev);
      });
     });
     kv.team_barrier();
  }//end of vlaplace_sphere_wk_contra

  template<int NUM_LEV_OUT, int NUM_LEV_IN = NUM_LEV_OUT, int NUM_LEV_REQUEST = NUM_LEV_OUT>
  KOKKOS_INLINE_FUNCTION void
  vlaplace_sphere_wk_contra (const KernelVariables &kv, const Real nu_ratio,
                             const typename ViewConst<ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_IN]>>::type& vector,
                             const ExecViewUnmanaged<Scalar [2][NP][NP][NUM_LEV_OUT]>& laplace) const
  {
    static_assert(NUM_LEV_REQUEST>=0, "Error! Invalid value for NUM_LEV_REQUEST.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_IN, "Error! Input view does not have enough levels.\n");
    static_assert(NUM_LEV_REQUEST<=NUM_LEV_OUT, "Error! Output view does not have enough levels.\n");

    vlaplace_sphere_wk_contra<NUM_LEV_OUT,NUM_LEV_IN>(kv, nu_ratio, vector, laplace, NUM_LEV_REQUEST);
  }//end of vlaplace_sphere_wk_contra

  // The buffers should be enough to handle any single call to any
  // single sphere operator.
  // One might prefer them to be private, but they are handy for
  // users of SphereOperators for other things, and using the same
  // memory buffers for multiple things gives performance in mem
  // b/w-limited computations.

  ExecViewManaged<Real   * [NUM_2D_VECTOR_BUFFERS][2][NP][NP]>                vector_buf_sl;
  ExecViewManaged<Scalar * [NUM_3D_SCALAR_BUFFERS][NP][NP][MAX_NUM_LEV]>      scalar_buf_ml;
  ExecViewManaged<Scalar * [NUM_3D_VECTOR_BUFFERS][2][NP][NP][MAX_NUM_LEV]>   vector_buf_ml;


  ExecViewManaged<const Real [NP][NP]>          dvv;
  ExecViewManaged<const Real [NP][NP]>          m_mp;

  ExecViewManaged<const Real * [NP][NP]>        m_spheremp;
  ExecViewManaged<const Real * [NP][NP]>        m_rspheremp;
  ExecViewManaged<const Real * [2][2][NP][NP]>  m_metinv;
  ExecViewManaged<const Real * [NP][NP]>        m_metdet;
  ExecViewManaged<const Real * [2][2][NP][NP]>  m_d;
  ExecViewManaged<const Real * [2][2][NP][NP]>  m_dinv;

  Real m_scale_factor_inv, m_laplacian_rigid_factor;
};

} // namespace Homme

#endif // HOMMEXX_SPHERE_OPERATORS_HPP
