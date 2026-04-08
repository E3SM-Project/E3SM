#include "eamxx_homme_process_interface.hpp"

// HOMMEXX includes
#include "Elements.hpp"
#include "ElementsDerivedState.hpp"
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
    const int ilev_pack)
{
  const auto u = v_dyn(ie,tl,0,igp,jgp,ilev_pack);
  const auto v = v_dyn(ie,tl,1,igp,jgp,ilev_pack);

  return vec_sph2cart(ie,0,icart,igp,jgp) * u
       + vec_sph2cart(ie,1,icart,igp,jgp) * v;
}

template <typename GradViewType, typename VecSph2CartViewType>
KOKKOS_INLINE_FUNCTION
Real cart_grad_to_local_component(
    const GradViewType& grad_Ux_dyn,
    const GradViewType& grad_Uy_dyn,
    const GradViewType& grad_Uz_dyn,
    const VecSph2CartViewType& vec_sph2cart,
    const int ie,
    const int idir,   // local velocity component index: 0,1
    const int jdir,   // local derivative direction index: 0,1
    const int igp,
    const int jgp,
    const int ilev)
{
  // grad_U[xyz]_dyn(ie,jdir,igp,jgp,ilev) are the Cartesian velocity-component
  // gradients with respect to local horizontal direction jdir.
  //
  // Project that Cartesian gradient vector back onto the local basis idir.
  return vec_sph2cart(ie,idir,0,igp,jgp) * grad_Ux_dyn(ie,jdir,igp,jgp,ilev)
       + vec_sph2cart(ie,idir,1,igp,jgp) * grad_Uy_dyn(ie,jdir,igp,jgp,ilev)
       + vec_sph2cart(ie,idir,2,igp,jgp) * grad_Uz_dyn(ie,jdir,igp,jgp,ilev);
}

} // anonymous namespace

void HommeDynamics::compute_horizontal_derivs_of_car_velocity ()
{
  using namespace Homme;

  constexpr int NGP  = HOMMEXX_NP;
  constexpr int VLEN = VECTOR_SIZE;

  const auto& c      = Context::singleton();
  const auto& state  = c.get<ElementsState>();
  const auto& geom   = c.get<ElementsGeometry>();
  const auto& ref_fe = c.get<ReferenceElement>();
  const auto& tl     = c.get<TimeLevel>();

  const int nelem        = m_dyn_grid->get_num_local_dofs() / (NGP*NGP);
  const int n0           = tl.n0;
  const int nlev_pack    = state.m_v.extent_int(5);
  const int nlev_scalar  = m_helper_fields.at("grad_Ux_dyn")
                             .template get_view<Real*****>().extent_int(4);
  const int nlevi_scalar = state.m_w_i.extent_int(4) * VLEN;

  const auto w_int_dyn = state.m_w_i;

  auto grad_Ux_dyn = m_helper_fields.at("grad_Ux_dyn").template get_view<Real*****>();
  auto grad_Uy_dyn = m_helper_fields.at("grad_Uy_dyn").template get_view<Real*****>();
  auto grad_Uz_dyn = m_helper_fields.at("grad_Uz_dyn").template get_view<Real*****>();

  const auto dvv              = ref_fe.get_deriv();
  const auto dinv             = geom.m_dinv;
  const auto vec_sph2cart     = geom.m_vec_sph2cart;
  const Real scale_factor_inv = 1.0 / geom.m_scale_factor;

  using TeamPolicy = Kokkos::TeamPolicy<KT::ExeSpace>;
  using MemberType = typename TeamPolicy::member_type;

  Kokkos::parallel_for(
      "compute_horizontal_derivs_of_car_velocity",
      TeamPolicy(nelem, Kokkos::AUTO()),
      KOKKOS_LAMBDA (const MemberType& team) {

    const int ie = team.league_rank();

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, NGP*NGP),
        [&] (const int idx) {

      const int igp = idx / NGP;
      const int jgp = idx % NGP;

      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(team, nlev_pack),
          [&] (const int ilev_pack) {

        Scalar dsdx_Ux = 0.0;
        Scalar dsdy_Ux = 0.0;
        Scalar dsdx_Uy = 0.0;
        Scalar dsdy_Uy = 0.0;
        Scalar dsdx_Uz = 0.0;
        Scalar dsdy_Uz = 0.0;

        for (int kgp = 0; kgp < NGP; ++kgp) {

          Scalar Ux_row = 0.0;
          Scalar Ux_col = 0.0;
          Scalar Uy_row = 0.0;
          Scalar Uy_col = 0.0;
          Scalar Uz_row = 0.0;
          Scalar Uz_col = 0.0;

          for (int s = 0; s < VLEN; ++s) {
            const int ilev = ilev_pack*VLEN + s;

            if (ilev < nlev_scalar) {

              // Horizontal wind contribution projected into Cartesian components.
              const Real Ux_row_h = horiz_wind_to_cart_component(
                  state.m_v, vec_sph2cart, ie, n0, 0, igp, kgp, ilev_pack)[s];

              const Real Ux_col_h = horiz_wind_to_cart_component(
                  state.m_v, vec_sph2cart, ie, n0, 0, kgp, jgp, ilev_pack)[s];

              const Real Uy_row_h = horiz_wind_to_cart_component(
                  state.m_v, vec_sph2cart, ie, n0, 1, igp, kgp, ilev_pack)[s];

              const Real Uy_col_h = horiz_wind_to_cart_component(
                  state.m_v, vec_sph2cart, ie, n0, 1, kgp, jgp, ilev_pack)[s];

              const Real Uz_row_h = horiz_wind_to_cart_component(
                  state.m_v, vec_sph2cart, ie, n0, 2, igp, kgp, ilev_pack)[s];

              const Real Uz_col_h = horiz_wind_to_cart_component(
                  state.m_v, vec_sph2cart, ie, n0, 2, kgp, jgp, ilev_pack)[s];

              // Vertical velocity contribution.
              // This is a hack, replace with ColOps later
              const auto w_row_pack = w_int_dyn(ie,n0,igp,kgp,ilev_pack);
              const auto w_col_pack = w_int_dyn(ie,n0,kgp,jgp,ilev_pack);

              const Real w_row = w_row_pack[s];
              const Real w_col = w_col_pack[s];

              // Project vertical velocity into Cartesian components.
              const Real wx_row = vec_sph2cart(ie,0,0,igp,kgp) * w_row;
              const Real wy_row = vec_sph2cart(ie,0,1,igp,kgp) * w_row;
              const Real wz_row = vec_sph2cart(ie,0,2,igp,kgp) * w_row;

              const Real wx_col = vec_sph2cart(ie,1,0,kgp,jgp) * w_col;
              const Real wy_col = vec_sph2cart(ie,1,1,kgp,jgp) * w_col;
              const Real wz_col = vec_sph2cart(ie,1,2,kgp,jgp) * w_col;

              Ux_row[s] = Ux_row_h + wx_row;
              Ux_col[s] = Ux_col_h + wx_col;

              Uy_row[s] = Uy_row_h + wy_row;
              Uy_col[s] = Uy_col_h + wy_col;

              Uz_row[s] = Uz_row_h + wz_row;
              Uz_col[s] = Uz_col_h + wz_col;
            }
          }

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

        for (int h = 0; h < 2; ++h) {
          const Scalar gx =
              dinv(ie,h,0,igp,jgp) * dsdx_Ux
            + dinv(ie,h,1,igp,jgp) * dsdy_Ux;

          const Scalar gy =
              dinv(ie,h,0,igp,jgp) * dsdx_Uy
            + dinv(ie,h,1,igp,jgp) * dsdy_Uy;

          const Scalar gz =
              dinv(ie,h,0,igp,jgp) * dsdx_Uz
            + dinv(ie,h,1,igp,jgp) * dsdy_Uz;

          for (int s = 0; s < VLEN; ++s) {
            const int ilev = ilev_pack*VLEN + s;
            if (ilev < nlev_scalar) {
              grad_Ux_dyn(ie,h,igp,jgp,ilev) = gx[s];
              grad_Uy_dyn(ie,h,igp,jgp,ilev) = gy[s];
              grad_Uz_dyn(ie,h,igp,jgp,ilev) = gz[s];
            }
          }
        }

      }); // ThreadVectorRange
    });   // TeamThreadRange
  });     // TeamPolicy

  Kokkos::fence();
}

void HommeDynamics::contract_to_local_strain2 ()
{
  using namespace Homme;

  constexpr int NGP  = HOMMEXX_NP;
  constexpr int VLEN = VECTOR_SIZE;

  const auto& c     = Context::singleton();
  const auto& geom  = c.get<ElementsGeometry>();
  const auto& elems = c.get<Elements>();

  const int nelem = m_dyn_grid->get_num_local_dofs() / (NGP*NGP);

  auto grad_Ux_dyn = m_helper_fields.at("grad_Ux_dyn").template get_view<Real*****>();
  auto grad_Uy_dyn = m_helper_fields.at("grad_Uy_dyn").template get_view<Real*****>();
  auto grad_Uz_dyn = m_helper_fields.at("grad_Uz_dyn").template get_view<Real*****>();

  auto strain2_dyn = m_helper_fields.at("strain2_dyn").template get_view<Real****>();

  const auto turb_strain2 = elems.m_derived.m_turb_strain2;
  const auto vec_sph2cart = geom.m_vec_sph2cart;

  const int nlev_scalar = grad_Ux_dyn.extent_int(4);
  const int nlev_pack   = turb_strain2.extent_int(3);

  using TeamPolicy = Kokkos::TeamPolicy<KT::ExeSpace>;
  using MemberType = typename TeamPolicy::member_type;

  Kokkos::parallel_for(
      "contract_to_local_strain2",
      TeamPolicy(nelem, Kokkos::AUTO()),
      KOKKOS_LAMBDA (const MemberType& team) {

    const int ie = team.league_rank();

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, NGP*NGP),
        [&] (const int idx) {

      const int igp = idx / NGP;
      const int jgp = idx % NGP;

      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(team, nlev_pack),
          [&] (const int ilev_pack) {

        Scalar strain2_pack(0.0);

        for (int s = 0; s < VLEN; ++s) {
          const int ilev = ilev_pack*VLEN + s;

          if (ilev < nlev_scalar) {
            const Real A00 = cart_grad_to_local_component(
                grad_Ux_dyn, grad_Uy_dyn, grad_Uz_dyn,
                vec_sph2cart, ie, 0, 0, igp, jgp, ilev);

            const Real A01 = cart_grad_to_local_component(
                grad_Ux_dyn, grad_Uy_dyn, grad_Uz_dyn,
                vec_sph2cart, ie, 0, 1, igp, jgp, ilev);

            const Real A10 = cart_grad_to_local_component(
                grad_Ux_dyn, grad_Uy_dyn, grad_Uz_dyn,
                vec_sph2cart, ie, 1, 0, igp, jgp, ilev);

            const Real A11 = cart_grad_to_local_component(
                grad_Ux_dyn, grad_Uy_dyn, grad_Uz_dyn,
                vec_sph2cart, ie, 1, 1, igp, jgp, ilev);

            const Real S00 = A00;
            const Real S11 = A11;
            const Real S01 = 0.5 * (A01 + A10);

            const Real strain2_val =
                2.0 * (S00*S00 + 2.0*S01*S01 + S11*S11);

            strain2_pack[s] = strain2_val;
            strain2_dyn(ie,igp,jgp,ilev) = strain2_val;
          } else {
            strain2_pack[s] = 0.0;
          }
        }

        turb_strain2(ie,igp,jgp,ilev_pack) = strain2_pack;
      }); // ThreadVectorRange
    });   // TeamThreadRange
  });     // TeamPolicy

  Kokkos::fence();
}

} // namespace scream
