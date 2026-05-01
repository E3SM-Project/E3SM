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
#include "PhysicalConstants.hpp"

// Scream includes
#include "dynamics/homme/homme_dimensions.hpp"
#include "ColumnOps.hpp"

// EKAT includes
#include <ekat_assert.hpp>

namespace scream
{

namespace {

// Build one Cartesian velocity component from the local horizontal wind only.
// vec_sph2cart layout is assumed to be:
//   (ie, isph, icart, igp, jgp)
// where isph=0,1,2 are the local spherical basis directions
// (zonal, meridional, radial) and icart=0,1,2 are x,y,z Cartesian components.
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
  const auto& elems  = c.get<Elements>();

  const int nelem       = m_dyn_grid->get_num_local_dofs() / (NGP*NGP);
  const int n0          = tl.n0;
  const int nlev_pack   = state.m_v.extent_int(5);
  const int nlev_scalar = m_helper_fields.at("grad_Ux_dyn")
                            .template get_view<Real*****>().extent_int(4);

  const auto w_int_dyn = state.m_w_i;

  using MidColumn = decltype(Homme::subview(elems.m_derived.m_turb_shear_strain3d, 0, 0, 0));
  using IntColumn = decltype(Homme::subview(state.m_w_i, 0, 0, 0, 0));

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

    // Construct the KernelVariables object needed by ColumnOps.
    KernelVariables kv(team, ie);

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, NGP*NGP),
        [&] (const int idx) {

      const int igp = idx / NGP;
      const int jgp = idx % NGP;

      // One packed vertical accumulator per midpoint level.
      Scalar dsdx_Ux[NUM_LEV];
      Scalar dsdy_Ux[NUM_LEV];
      Scalar dsdx_Uy[NUM_LEV];
      Scalar dsdy_Uy[NUM_LEV];
      Scalar dsdx_Uz[NUM_LEV];
      Scalar dsdy_Uz[NUM_LEV];

      for (int k = 0; k < nlev_pack; ++k) {
        dsdx_Ux[k] = 0.0;
        dsdy_Ux[k] = 0.0;
        dsdx_Uy[k] = 0.0;
        dsdy_Uy[k] = 0.0;
        dsdx_Uz[k] = 0.0;
        dsdy_Uz[k] = 0.0;
      }

      for (int kgp = 0; kgp < NGP; ++kgp) {
        // Sweep the element stencil to build local horizontal derivatives of the
        // full 3D Cartesian velocity field at each GLL point.

        // Full interface columns for w at the two stencil points.
        const auto w_int_row = Homme::subview(w_int_dyn, ie, n0, igp, kgp);
        const auto w_int_col = Homme::subview(w_int_dyn, ie, n0, kgp, jgp);

        // Temporary midpoint columns.
        Scalar w_mid_row_buf[NUM_LEV];
        Scalar w_mid_col_buf[NUM_LEV];

        MidColumn w_mid_row(w_mid_row_buf);
        MidColumn w_mid_col(w_mid_col_buf);

        // Put vertical velocity on the midpoint grid
        ColumnOps::compute_midpoint_values(kv, w_int_row, w_mid_row);
        ColumnOps::compute_midpoint_values(kv, w_int_col, w_mid_col);

        team.team_barrier();

        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(team, nlev_pack),
            [&] (const int ilev_pack) {

          Scalar Ux_row = 0.0;
          Scalar Ux_col = 0.0;
          Scalar Uy_row = 0.0;
          Scalar Uy_col = 0.0;
          Scalar Uz_row = 0.0;
          Scalar Uz_col = 0.0;

          for (int s = 0; s < VLEN; ++s) {
            const int ilev = ilev_pack*VLEN + s;

            if (ilev < nlev_scalar) {

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

              // Add the vertical motion to the horizontal wind so the gradients
              // are taken on the full resolved 3D velocity, not just (u,v).
              const Real w_row = w_mid_row(ilev_pack)[s];
              const Real w_col = w_mid_col(ilev_pack)[s];

              const Real wx_row = vec_sph2cart(ie,2,0,igp,kgp) * w_row;
              const Real wy_row = vec_sph2cart(ie,2,1,igp,kgp) * w_row;
              const Real wz_row = vec_sph2cart(ie,2,2,igp,kgp) * w_row;

              const Real wx_col = vec_sph2cart(ie,2,0,kgp,jgp) * w_col;
              const Real wy_col = vec_sph2cart(ie,2,1,kgp,jgp) * w_col;
              const Real wz_col = vec_sph2cart(ie,2,2,kgp,jgp) * w_col;

              Ux_row[s] = Ux_row_h + wx_row;
              Ux_col[s] = Ux_col_h + wx_col;

              Uy_row[s] = Uy_row_h + wy_row;
              Uy_col[s] = Uy_col_h + wy_col;

              Uz_row[s] = Uz_row_h + wz_row;
              Uz_col[s] = Uz_col_h + wz_col;
            }
          }

          dsdx_Ux[ilev_pack] += dvv(jgp,kgp) * Ux_row;
          dsdy_Ux[ilev_pack] += dvv(igp,kgp) * Ux_col;

          dsdx_Uy[ilev_pack] += dvv(jgp,kgp) * Uy_row;
          dsdy_Uy[ilev_pack] += dvv(igp,kgp) * Uy_col;

          dsdx_Uz[ilev_pack] += dvv(jgp,kgp) * Uz_row;
          dsdy_Uz[ilev_pack] += dvv(igp,kgp) * Uz_col;
        });

        team.team_barrier();
      }

      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(team, nlev_pack),
          [&] (const int ilev_pack) {

        // Convert reference-element derivatives into physical horizontal
        // gradients using the local metric tensor on the curved element.
        Scalar gx0 = (dinv(ie,0,0,igp,jgp) * dsdx_Ux[ilev_pack]
                    + dinv(ie,0,1,igp,jgp) * dsdy_Ux[ilev_pack]) * scale_factor_inv;

        Scalar gy0 = (dinv(ie,0,0,igp,jgp) * dsdx_Uy[ilev_pack]
                    + dinv(ie,0,1,igp,jgp) * dsdy_Uy[ilev_pack]) * scale_factor_inv;

        Scalar gz0 = (dinv(ie,0,0,igp,jgp) * dsdx_Uz[ilev_pack]
                    + dinv(ie,0,1,igp,jgp) * dsdy_Uz[ilev_pack]) * scale_factor_inv;

        Scalar gx1 = (dinv(ie,1,0,igp,jgp) * dsdx_Ux[ilev_pack]
                    + dinv(ie,1,1,igp,jgp) * dsdy_Ux[ilev_pack]) * scale_factor_inv;

        Scalar gy1 = (dinv(ie,1,0,igp,jgp) * dsdx_Uy[ilev_pack]
                    + dinv(ie,1,1,igp,jgp) * dsdy_Uy[ilev_pack]) * scale_factor_inv;

        Scalar gz1 = (dinv(ie,1,0,igp,jgp) * dsdx_Uz[ilev_pack]
                    + dinv(ie,1,1,igp,jgp) * dsdy_Uz[ilev_pack]) * scale_factor_inv;

        for (int s = 0; s < VLEN; ++s) {
          const int ilev = ilev_pack*VLEN + s;
          if (ilev < nlev_scalar) {
            grad_Ux_dyn(ie,0,igp,jgp,ilev) = gx0[s];
            grad_Uy_dyn(ie,0,igp,jgp,ilev) = gy0[s];
            grad_Uz_dyn(ie,0,igp,jgp,ilev) = gz0[s];

            grad_Ux_dyn(ie,1,igp,jgp,ilev) = gx1[s];
            grad_Uy_dyn(ie,1,igp,jgp,ilev) = gy1[s];
            grad_Uz_dyn(ie,1,igp,jgp,ilev) = gz1[s];
          }
        }
      });

    });
  });

  Kokkos::fence();
}

void HommeDynamics::compute_vertical_derivs ()
{
  using namespace Homme;

  constexpr int NGP  = HOMMEXX_NP;
  constexpr int VLEN = VECTOR_SIZE;

  const auto& c      = Context::singleton();
  const auto& state  = c.get<ElementsState>();
  const auto& tl     = c.get<TimeLevel>();
  const auto& elems  = c.get<Elements>();

  const int nelem       = m_dyn_grid->get_num_local_dofs() / (NGP*NGP);
  const int n0          = tl.n0;
  const int nlev_pack   = state.m_v.extent_int(5);
  const int nlev_scalar = m_helper_fields.at("grad_dz_dyn")
                            .template get_view<Real*****>().extent_int(4);

  const auto v_dyn        = state.m_v;
  const auto w_int_dyn    = state.m_w_i;
  const auto phinh_i_dyn  = state.m_phinh_i;

  using MidColumn = decltype(Homme::subview(elems.m_derived.m_turb_shear_strain3d, 0, 0, 0));
  using IntColumn = decltype(Homme::subview(state.m_w_i, 0, 0, 0, 0));

  auto grad_vertical = m_helper_fields.at("grad_dz_dyn").template get_view<Real*****>();

  using TeamPolicy = Kokkos::TeamPolicy<KT::ExeSpace>;
  using MemberType = typename TeamPolicy::member_type;

  Kokkos::parallel_for(
      "compute_vertical_derivs",
      TeamPolicy(nelem, Kokkos::AUTO()),
      KOKKOS_LAMBDA (const MemberType& team) {

    const int ie = team.league_rank();

    KernelVariables kv(team, ie);

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, NGP*NGP),
        [&] (const int idx) {

      const int igp = idx / NGP;
      const int jgp = idx % NGP;

      // Column views
      const auto u_mid = Homme::subview(v_dyn,     ie, n0, 0, igp, jgp);
      const auto v_mid = Homme::subview(v_dyn,     ie, n0, 1, igp, jgp);
      const auto w_int = Homme::subview(w_int_dyn, ie, n0,    igp, jgp);
      const auto phi_int = Homme::subview(phinh_i_dyn, ie, n0, igp, jgp);

      // Temporary midpoint and interface columns.
      Scalar w_mid_buf[NUM_LEV];
      Scalar grad_u_int_buf[NUM_LEV_P];
      Scalar grad_v_int_buf[NUM_LEV_P];
      Scalar grad_w_int_buf[NUM_LEV_P];
      Scalar grad_u_mid_buf[NUM_LEV];
      Scalar grad_v_mid_buf[NUM_LEV];
      Scalar grad_w_mid_buf[NUM_LEV];

      MidColumn w_mid(w_mid_buf);
      IntColumn grad_u_int(grad_u_int_buf);
      IntColumn grad_v_int(grad_v_int_buf);
      IntColumn grad_w_int(grad_w_int_buf);
      MidColumn grad_u_mid(grad_u_mid_buf);
      MidColumn grad_v_mid(grad_v_mid_buf);
      MidColumn grad_w_mid(grad_w_mid_buf);

      // ------------------------------------------
      // 1) Put w on midpoint levels
      // ------------------------------------------
      ColumnOps::compute_midpoint_values(kv, w_int, w_mid);

      team.team_barrier();

      // ------------------------------------------
      // 2) Compute vertical gradients on interface levels
      // ------------------------------------------
      //
      // Use HOMME's interface geopotential to define the distance between
      // adjacent midpoint levels, matching SHOC's interface-gradient staggering.
      constexpr Real gravit = PhysicalConstants::g;

      for (int ilev = 0; ilev <= nlev_scalar; ++ilev) {
        const int ilev_pack = ilev / VLEN;
        const int s = ilev % VLEN;

        if (ilev == 0 || ilev == nlev_scalar) {
          grad_u_int(ilev_pack)[s] = 0.0;
          grad_v_int(ilev_pack)[s] = 0.0;
          grad_w_int(ilev_pack)[s] = 0.0;
        } else {
          const int ilev_above = ilev - 1;
          const int pack_above = ilev_above / VLEN;
          const int s_above = ilev_above % VLEN;

          const Real phi_int_above = phi_int(pack_above)[s_above];
          const Real phi_int_below = (s < VLEN-1) ? phi_int(ilev_pack)[s+1] : phi_int(ilev_pack+1)[0];
          const Real dz_int = (phi_int_above - phi_int_below) / (2.0*gravit);
          const Real inv_dz_int = 1.0 / dz_int;

          // HOMME levels are indexed top to bottom, so the numerator is
          // above-minus-below for derivatives with respect to increasing height.
          grad_u_int(ilev_pack)[s] =
            (u_mid(pack_above)[s_above] - u_mid(ilev_pack)[s]) * inv_dz_int;

          grad_v_int(ilev_pack)[s] =
            (v_mid(pack_above)[s_above] - v_mid(ilev_pack)[s]) * inv_dz_int;

          grad_w_int(ilev_pack)[s] =
            (w_mid(pack_above)[s_above] - w_mid(ilev_pack)[s]) * inv_dz_int;
        }
      }

      team.team_barrier();

      // ------------------------------------------
      // 3) Interpolate interface gradients back to midpoints
      // ------------------------------------------
      ColumnOps::compute_midpoint_values(kv, grad_u_int, grad_u_mid);
      ColumnOps::compute_midpoint_values(kv, grad_v_int, grad_v_mid);
      ColumnOps::compute_midpoint_values(kv, grad_w_int, grad_w_mid);

      team.team_barrier();

      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(team, nlev_pack),
          [&] (const int ilev_pack) {

        for (int s = 0; s < VLEN; ++s) {
          const int ilev = ilev_pack*VLEN + s;

          if (ilev < nlev_scalar) {
            grad_vertical(ie,0,igp,jgp,ilev) = grad_u_mid(ilev_pack)[s];
            grad_vertical(ie,1,igp,jgp,ilev) = grad_v_mid(ilev_pack)[s];
            grad_vertical(ie,2,igp,jgp,ilev) = grad_w_mid(ilev_pack)[s];
          }
        }
      });

    });
  });

  Kokkos::fence();
}

void HommeDynamics::contract_to_local_strain3d ()
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

  auto grad_vertical = m_helper_fields.at("grad_dz_dyn").template get_view<Real*****>();

  auto shear_strain3d_dyn = m_helper_fields.at("shear_strain3d_dyn").template get_view<Real****>();

  const auto turb_shear_strain3d = elems.m_derived.m_turb_shear_strain3d;
  const auto vec_sph2cart = geom.m_vec_sph2cart;

  const int nlev_scalar = grad_Ux_dyn.extent_int(4);
  const int nlev_pack   = turb_shear_strain3d.extent_int(3);

  using TeamPolicy = Kokkos::TeamPolicy<KT::ExeSpace>;
  using MemberType = typename TeamPolicy::member_type;

  Kokkos::parallel_for(
      "contract_to_local_strain3d",
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

        Scalar shear_strain3d_pack(0.0);

        for (int s = 0; s < VLEN; ++s) {
          const int ilev = ilev_pack*VLEN + s;

          if (ilev < nlev_scalar) {
            // Project the Cartesian gradient tensor back into the local
            // element basis so SHOC receives a physically meaningful shear rate.
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

            const Real A20 = cart_grad_to_local_component(
                grad_Ux_dyn, grad_Uy_dyn, grad_Uz_dyn,
                vec_sph2cart, ie, 2, 0, igp, jgp, ilev);

            const Real A21 = cart_grad_to_local_component(
                grad_Ux_dyn, grad_Uy_dyn, grad_Uz_dyn,
                vec_sph2cart, ie, 2, 1, igp, jgp, ilev);

            const Real A02 = grad_vertical(ie,0,igp,jgp,ilev); // du/dz
            const Real A12 = grad_vertical(ie,1,igp,jgp,ilev); // dv/dz
            const Real A22 = grad_vertical(ie,2,igp,jgp,ilev); // dw/dz

            const Real S00 = A00;
            const Real S11 = A11;
            const Real S22 = A22;

            const Real S01 = 0.5 * (A01 + A10);
            const Real S02 = 0.5 * (A02 + A20);
            const Real S12 = 0.5 * (A12 + A21);

            // Form 2*S_ij*S_ij, the shear-strain invariant that enters the
            // turbulence production term passed from dynamics to physics.
            const Real shear_strain3d_val =
                2.0 * (S00*S00 + S11*S11 + S22*S22
                     + 2.0*S01*S01
                     + 2.0*S02*S02
                     + 2.0*S12*S12);

            shear_strain3d_pack[s] = shear_strain3d_val;
            shear_strain3d_dyn(ie,igp,jgp,ilev) = shear_strain3d_val;
          } else {
            shear_strain3d_pack[s] = 0.0;
          }
        }

        turb_shear_strain3d(ie,igp,jgp,ilev_pack) = shear_strain3d_pack;
      });
    });
  });

  Kokkos::fence();
}

} // namespace scream
