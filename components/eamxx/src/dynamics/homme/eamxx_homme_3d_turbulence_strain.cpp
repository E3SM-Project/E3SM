#include "eamxx_homme_process_interface.hpp"

// HOMMEXX includes
#include "Context.hpp"
#include "ElementsGeometry.hpp"
#include "ElementsState.hpp"
#include "ReferenceElement.hpp"
#include "TimeLevel.hpp"
#include "Types.hpp"

// Scream includes
#include "dynamics/homme/homme_dimensions.hpp"

// EKAT includes
#include <ekat_assert.hpp>

namespace scream
{

namespace {

// Build one Cartesian velocity component from the local horizontal wind only.
// vec_sph2cart layout is assumed to be:
//   (ie, ihoriz, icart, igp, jgp)
// where ihoriz=0,1 are the two local horizontal velocity components
// and icart=0,1,2 are x,y,z Cartesian components.
template <typename VDynViewType, typename VecSph2CartViewType>
KOKKOS_INLINE_FUNCTION
auto horiz_wind_to_cart_component(
    const VDynViewType& v_dyn,
    const VecSph2CartViewType& vec_sph2cart,
    const int ie,
    const int tl,
    const int icart,
    const int igp,
    const int jgp,
    const int ilev)
{
  const auto u = v_dyn(ie,tl,0,igp,jgp,ilev);
  const auto v = v_dyn(ie,tl,1,igp,jgp,ilev);

  return vec_sph2cart(ie,0,icart,igp,jgp) * u
       + vec_sph2cart(ie,1,icart,igp,jgp) * v;
}

} // anonymous namespace

void HommeDynamics::compute_horizontal_derivs_of_car_velocity ()
{
  using namespace Homme;

  constexpr int NGP = HOMMEXX_NP;
  constexpr int NVL = HOMMEXX_NUM_LEV;

  const auto& c        = Context::singleton();
  const auto& state    = c.get<ElementsState>();
  const auto& geom     = c.get<ElementsGeometry>();
  const auto& ref_fe   = c.get<ReferenceElement>();
  const auto& tl       = c.get<TimeLevel>();

  const int nelem = m_dyn_grid->get_num_local_dofs() / (NGP*NGP);
  const int n0    = tl.n0;

  // Output views:
  //   grad_Ux_dyn(ie,dir,igp,jgp,ilev), dir=0,1
  //   grad_Uy_dyn(ie,dir,igp,jgp,ilev)
  //   grad_Uz_dyn(ie,dir,igp,jgp,ilev)

  auto grad_Ux_dyn = m_helper_fields.at("grad_Ux_dyn").template get_view<Homme::Scalar*[2][NGP][NGP][NVL]>();
  auto grad_Uy_dyn = m_helper_fields.at("grad_Uy_dyn").template get_view<Homme::Scalar*[2][NGP][NGP][NVL]>();
  auto grad_Uz_dyn = m_helper_fields.at("grad_Uz_dyn").template get_view<Homme::Scalar*[2][NGP][NGP][NVL]>();

  const auto dvv             = ref_fe.get_deriv();
  const auto dinv            = geom.m_dinv;
  const auto vec_sph2cart    = geom.m_vec_sph2cart;
  const Real scale_factor_inv = 1.0 / geom.m_scale_factor;

  // Flatten over element, GLL point, and level.
  Kokkos::parallel_for(
      "compute_horizontal_derivs_of_car_velocity",
      Kokkos::RangePolicy<KT::ExeSpace>(0, nelem*NGP*NGP*NVL),
      KOKKOS_LAMBDA (const int idx) {

    const int ie   =  idx / (NGP*NGP*NVL);
    const int rem1 =  idx % (NGP*NGP*NVL);
    const int igp  =  rem1 / (NGP*NVL);
    const int rem2 =  rem1 % (NGP*NVL);
    const int jgp  =  rem2 / NVL;
    const int ilev =  rem2 % NVL;

    // ------------------------------------------------------------------
    // Compute the two intermediate directional derivatives for each of
    // the three Cartesian velocity components.
    //
    // This mirrors HOMME's gradient_sphere algebra:
    //   temp0(igp,jgp) = sum_k dvv(jgp,k) * scalar(igp,k)
    //   temp1(igp,jgp) = sum_k dvv(igp,k) * scalar(k,jgp)
    // then
    //   grad(h,igp,jgp) = dinv(h,0,igp,jgp)*temp0 + dinv(h,1,igp,jgp)*temp1
    //
    // Here scalar is one of Ux, Uy, Uz.
    // ------------------------------------------------------------------

    Scalar dsdx_Ux = 0.0;
    Scalar dsdy_Ux = 0.0;
    Scalar dsdx_Uy = 0.0;
    Scalar dsdy_Uy = 0.0;
    Scalar dsdx_Uz = 0.0;
    Scalar dsdy_Uz = 0.0;

    for (int kgp = 0; kgp < NGP; ++kgp) {
      const Scalar Ux_row = horiz_wind_to_cart_component(
          state.m_v, vec_sph2cart, ie, n0, 0, igp, kgp, ilev);
      const Scalar Ux_col = horiz_wind_to_cart_component(
          state.m_v, vec_sph2cart, ie, n0, 0, kgp, jgp, ilev);

      const Scalar Uy_row = horiz_wind_to_cart_component(
          state.m_v, vec_sph2cart, ie, n0, 1, igp, kgp, ilev);
      const Scalar Uy_col = horiz_wind_to_cart_component(
          state.m_v, vec_sph2cart, ie, n0, 1, kgp, jgp, ilev);

      const Scalar Uz_row = horiz_wind_to_cart_component(
          state.m_v, vec_sph2cart, ie, n0, 2, igp, kgp, ilev);
      const Scalar Uz_col = horiz_wind_to_cart_component(
          state.m_v, vec_sph2cart, ie, n0, 2, kgp, jgp, ilev);

      dsdx_Ux += dvv(jgp,kgp) * Ux_row;
      dsdy_Ux += dvv(igp,kgp) * Ux_col;

      dsdx_Uy += dvv(jgp,kgp) * Uy_row;
      dsdy_Uy += dvv(igp,kgp) * Uy_col;

      dsdx_Uz += dvv(jgp,kgp) * Uz_row;
      dsdy_Uz += dvv(igp,kgp) * Uz_col;
    }

    dsdx_Ux *= scale_factor_inv;
    dsdy_Ux *= scale_factor_inv;
    dsdx_Uy *= scale_factor_inv;
    dsdy_Uy *= scale_factor_inv;
    dsdx_Uz *= scale_factor_inv;
    dsdy_Uz *= scale_factor_inv;

    // Final transformed horizontal gradients in the local 2-dir basis.
    for (int h = 0; h < 2; ++h) {
      grad_Ux_dyn(ie,h,igp,jgp,ilev) =
          dinv(ie,h,0,igp,jgp) * dsdx_Ux
        + dinv(ie,h,1,igp,jgp) * dsdy_Ux;

      grad_Uy_dyn(ie,h,igp,jgp,ilev) =
          dinv(ie,h,0,igp,jgp) * dsdx_Uy
        + dinv(ie,h,1,igp,jgp) * dsdy_Uy;

      grad_Uz_dyn(ie,h,igp,jgp,ilev) =
          dinv(ie,h,0,igp,jgp) * dsdx_Uz
        + dinv(ie,h,1,igp,jgp) * dsdy_Uz;
    }
  });

  Kokkos::fence();
}

} // namespace scream
