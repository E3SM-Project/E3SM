#include "eamxx_homme_process_interface.hpp"

// HOMMEXX includes
#include "Context.hpp"
#include "ElementsGeometry.hpp"
#include "ElementsState.hpp"
#include "ReferenceElement.hpp"
#include "TimeLevel.hpp"
#include "Types.hpp"
#include "utilities/ViewUtils.hpp"

// Scream includes
#include "dynamics/homme/homme_dimensions.hpp"
#include "share/util/eamxx_column_ops.hpp"

// EKAT includes
#include <ekat_assert.hpp>

namespace scream
{

namespace {

// Project a local velocity vector (u,v,w) into one Cartesian component.
template <typename BasisViewType>
KOKKOS_INLINE_FUNCTION
Real local_to_cart_component(
    const BasisViewType& basis_sph2cart,
    const Real u,
    const Real v,
    const Real w)
{
  return basis_sph2cart(0) * u
       + basis_sph2cart(1) * v
       + basis_sph2cart(2) * w;
}

} // anonymous namespace

void HommeDynamics::compute_horizontal_derivs_of_car_velocity ()
{
  using namespace Homme;
  using ColOps = scream::ColumnOps<DefaultDevice, Real>;

  constexpr int NGP  = HOMMEXX_NP;

  const auto& c      = Context::singleton();
  const auto& state  = c.get<ElementsState>();
  const auto& geom   = c.get<ElementsGeometry>();
  const auto& ref_fe = c.get<ReferenceElement>();
  const auto& tl     = c.get<TimeLevel>();

  const int nelem       = m_dyn_grid->get_num_local_dofs() / (NGP*NGP);
  const int n0          = tl.n0;
  const int nlev_scalar = m_helper_fields.at("grad_Ux_dyn")
                            .template get_view<Real*****>().extent_int(4);

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
  const int ncols = nelem*NGP*NGP;
  const TeamPolicy policy(ncols, Kokkos::AUTO());
  const auto w_mid_row_all = m_w_mid_row_all;
  const auto w_mid_col_all = m_w_mid_col_all;
  const auto dsdx_Ux_all = m_dsdx_Ux_all;
  const auto dsdy_Ux_all = m_dsdy_Ux_all;
  const auto dsdx_Uy_all = m_dsdx_Uy_all;
  const auto dsdy_Uy_all = m_dsdy_Uy_all;
  const auto dsdx_Uz_all = m_dsdx_Uz_all;
  const auto dsdy_Uz_all = m_dsdy_Uz_all;

  Kokkos::parallel_for(
      "compute_horizontal_derivs_of_car_velocity",
      policy,
      KOKKOS_LAMBDA (const MemberType& team) {

    const int ie  =  team.league_rank() / (NGP*NGP);
    const int igp = (team.league_rank() / NGP) % NGP;
    const int jgp =  team.league_rank() % NGP;
    const int icol = team.league_rank();

    // Grab the scratch storage associated with this (ie,igp,jgp) column.
    const auto w_mid_row = Kokkos::subview(w_mid_row_all, icol, Kokkos::ALL());
    const auto w_mid_col = Kokkos::subview(w_mid_col_all, icol, Kokkos::ALL());
    using Unmanaged = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
    using Real1dU = typename ColOps::template view_1d<Real, Unmanaged>;
    const Real1dU w_mid_row_1d(w_mid_row.data(), nlev_scalar);
    const Real1dU w_mid_col_1d(w_mid_col.data(), nlev_scalar);
    const auto dsdx_Ux = Kokkos::subview(dsdx_Ux_all, icol, Kokkos::ALL());
    const auto dsdy_Ux = Kokkos::subview(dsdy_Ux_all, icol, Kokkos::ALL());
    const auto dsdx_Uy = Kokkos::subview(dsdx_Uy_all, icol, Kokkos::ALL());
    const auto dsdy_Uy = Kokkos::subview(dsdy_Uy_all, icol, Kokkos::ALL());
    const auto dsdx_Uz = Kokkos::subview(dsdx_Uz_all, icol, Kokkos::ALL());
    const auto dsdy_Uz = Kokkos::subview(dsdy_Uz_all, icol, Kokkos::ALL());

    // Accumulate reference-element derivatives in the two local horizontal directions.
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_scalar), [&] (const int ilev) {
      dsdx_Ux(ilev) = 0;
      dsdy_Ux(ilev) = 0;
      dsdx_Uy(ilev) = 0;
      dsdy_Uy(ilev) = 0;
      dsdx_Uz(ilev) = 0;
      dsdy_Uz(ilev) = 0;
    });
    team.team_barrier();

    for (int kgp = 0; kgp < NGP; ++kgp) {
      // The horizontal stencil uses interface w, so first put the two stencil
      // columns of vertical velocity onto midpoint levels.
      const auto w_row_i = Kokkos::subview(w_int_dyn, ie, n0, igp, kgp, Kokkos::ALL());
      const auto w_col_i = Kokkos::subview(w_int_dyn, ie, n0, kgp, jgp, Kokkos::ALL());
      const auto w_row_i_real = Homme::viewAsReal(w_row_i);
      const auto w_col_i_real = Homme::viewAsReal(w_col_i);

      ColOps::compute_midpoint_values(team, nlev_scalar, w_row_i_real, w_mid_row_1d);
      team.team_barrier();
      ColOps::compute_midpoint_values(team, nlev_scalar, w_col_i_real, w_mid_col_1d);
      team.team_barrier();

      const auto row_x = Kokkos::subview(vec_sph2cart, ie, Kokkos::ALL(), 0, igp, kgp);
      const auto row_y = Kokkos::subview(vec_sph2cart, ie, Kokkos::ALL(), 1, igp, kgp);
      const auto row_z = Kokkos::subview(vec_sph2cart, ie, Kokkos::ALL(), 2, igp, kgp);
      const auto col_x = Kokkos::subview(vec_sph2cart, ie, Kokkos::ALL(), 0, kgp, jgp);
      const auto col_y = Kokkos::subview(vec_sph2cart, ie, Kokkos::ALL(), 1, kgp, jgp);
      const auto col_z = Kokkos::subview(vec_sph2cart, ie, Kokkos::ALL(), 2, kgp, jgp);
      const auto u_row_view =
          Homme::viewAsReal(Kokkos::subview(state.m_v, ie, n0, 0, igp, kgp, Kokkos::ALL()));
      const auto v_row_view =
          Homme::viewAsReal(Kokkos::subview(state.m_v, ie, n0, 1, igp, kgp, Kokkos::ALL()));
      const auto u_col_view =
          Homme::viewAsReal(Kokkos::subview(state.m_v, ie, n0, 0, kgp, jgp, Kokkos::ALL()));
      const auto v_col_view =
          Homme::viewAsReal(Kokkos::subview(state.m_v, ie, n0, 1, kgp, jgp, Kokkos::ALL()));

      // Build the full Cartesian velocity on the two stencil lines and apply
      // the derivative matrix weights along each local direction.
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_scalar), [&] (const int ilev) {
        const Real u_row = u_row_view(ilev);
        const Real v_row = v_row_view(ilev);
        const Real w_row = w_mid_row(ilev);

        const Real u_col = u_col_view(ilev);
        const Real v_col = v_col_view(ilev);
        const Real w_col = w_mid_col(ilev);

        dsdx_Ux(ilev) += dvv(jgp,kgp) * local_to_cart_component(row_x, u_row, v_row, w_row);
        dsdy_Ux(ilev) += dvv(igp,kgp) * local_to_cart_component(col_x, u_col, v_col, w_col);

        dsdx_Uy(ilev) += dvv(jgp,kgp) * local_to_cart_component(row_y, u_row, v_row, w_row);
        dsdy_Uy(ilev) += dvv(igp,kgp) * local_to_cart_component(col_y, u_col, v_col, w_col);

        dsdx_Uz(ilev) += dvv(jgp,kgp) * local_to_cart_component(row_z, u_row, v_row, w_row);
        dsdy_Uz(ilev) += dvv(igp,kgp) * local_to_cart_component(col_z, u_col, v_col, w_col);
      });
      team.team_barrier();
    }

    // Convert the reference-element derivatives into physical horizontal
    // gradients using the inverse metric tensor on this curved element.
    const auto dinv_ij = Kokkos::subview(dinv, ie, Kokkos::ALL(), Kokkos::ALL(), igp, jgp);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_scalar), [&] (const int ilev) {
      grad_Ux_dyn(ie,0,igp,jgp,ilev) = (dinv_ij(0,0) * dsdx_Ux(ilev) + dinv_ij(0,1) * dsdy_Ux(ilev)) * scale_factor_inv;
      grad_Uy_dyn(ie,0,igp,jgp,ilev) = (dinv_ij(0,0) * dsdx_Uy(ilev) + dinv_ij(0,1) * dsdy_Uy(ilev)) * scale_factor_inv;
      grad_Uz_dyn(ie,0,igp,jgp,ilev) = (dinv_ij(0,0) * dsdx_Uz(ilev) + dinv_ij(0,1) * dsdy_Uz(ilev)) * scale_factor_inv;

      grad_Ux_dyn(ie,1,igp,jgp,ilev) = (dinv_ij(1,0) * dsdx_Ux(ilev) + dinv_ij(1,1) * dsdy_Ux(ilev)) * scale_factor_inv;
      grad_Uy_dyn(ie,1,igp,jgp,ilev) = (dinv_ij(1,0) * dsdx_Uy(ilev) + dinv_ij(1,1) * dsdy_Uy(ilev)) * scale_factor_inv;
      grad_Uz_dyn(ie,1,igp,jgp,ilev) = (dinv_ij(1,0) * dsdx_Uz(ilev) + dinv_ij(1,1) * dsdy_Uz(ilev)) * scale_factor_inv;
    });
  });

  Kokkos::fence();
}

void HommeDynamics::compute_local_strain_components3d ()
{
  using namespace Homme;

  constexpr int NGP  = HOMMEXX_NP;

  const auto& c    = Context::singleton();
  const auto& geom = c.get<ElementsGeometry>();

  const int nelem = m_dyn_grid->get_num_local_dofs() / (NGP*NGP);

  auto grad_Ux_dyn = m_helper_fields.at("grad_Ux_dyn").template get_view<Real*****>();
  auto grad_Uy_dyn = m_helper_fields.at("grad_Uy_dyn").template get_view<Real*****>();
  auto grad_Uz_dyn = m_helper_fields.at("grad_Uz_dyn").template get_view<Real*****>();

  auto shear_components_dyn = m_helper_fields.at("shear_strain3d_components_dyn").template get_view<Real*****>();

  const auto vec_sph2cart = geom.m_vec_sph2cart;

  const int nlev_scalar = grad_Ux_dyn.extent_int(4);

  using Policy = Kokkos::MDRangePolicy<KT::ExeSpace, Kokkos::Rank<4>>;
  const Policy policy({0, 0, 0, 0}, {nelem, NGP, NGP, nlev_scalar});

  Kokkos::parallel_for(
      "compute_local_strain_components3d",
      policy,
      KOKKOS_LAMBDA (const int ie, const int igp, const int jgp, const int ilev) {

    // The stored gradients are Cartesian components differentiated along the
    // two local horizontal directions. Project them back into the local basis
    // so SHOC receives the six local shear-tensor components it expects.
    const Real gx0 = grad_Ux_dyn(ie,0,igp,jgp,ilev);
    const Real gy0 = grad_Uy_dyn(ie,0,igp,jgp,ilev);
    const Real gz0 = grad_Uz_dyn(ie,0,igp,jgp,ilev);

    const Real gx1 = grad_Ux_dyn(ie,1,igp,jgp,ilev);
    const Real gy1 = grad_Uy_dyn(ie,1,igp,jgp,ilev);
    const Real gz1 = grad_Uz_dyn(ie,1,igp,jgp,ilev);

    const Real b0_0 = vec_sph2cart(ie, 0, 0, igp, jgp);
    const Real b0_1 = vec_sph2cart(ie, 0, 1, igp, jgp);
    const Real b0_2 = vec_sph2cart(ie, 0, 2, igp, jgp);

    const Real b1_0 = vec_sph2cart(ie, 1, 0, igp, jgp);
    const Real b1_1 = vec_sph2cart(ie, 1, 1, igp, jgp);
    const Real b1_2 = vec_sph2cart(ie, 1, 2, igp, jgp);

    const Real b2_0 = vec_sph2cart(ie, 2, 0, igp, jgp);
    const Real b2_1 = vec_sph2cart(ie, 2, 1, igp, jgp);
    const Real b2_2 = vec_sph2cart(ie, 2, 2, igp, jgp);

    shear_components_dyn(ie,0,igp,jgp,ilev) = b0_0 * gx0 + b0_1 * gy0 + b0_2 * gz0;
    shear_components_dyn(ie,1,igp,jgp,ilev) = b0_0 * gx1 + b0_1 * gy1 + b0_2 * gz1;
    shear_components_dyn(ie,2,igp,jgp,ilev) = b1_0 * gx0 + b1_1 * gy0 + b1_2 * gz0;
    shear_components_dyn(ie,3,igp,jgp,ilev) = b1_0 * gx1 + b1_1 * gy1 + b1_2 * gz1;
    shear_components_dyn(ie,4,igp,jgp,ilev) = b2_0 * gx0 + b2_1 * gy0 + b2_2 * gz0;
    shear_components_dyn(ie,5,igp,jgp,ilev) = b2_0 * gx1 + b2_1 * gy1 + b2_2 * gz1;
  });

  Kokkos::fence();
}

} // namespace scream
