#include "eamxx_homme_process_interface.hpp"

// HOMMEXX includes
#include "Context.hpp"
#include "ElementsGeometry.hpp"
#include "ElementsState.hpp"
#include "HybridVCoord.hpp"
#include "ReferenceElement.hpp"
#include "TimeLevel.hpp"
#include "Tracers.hpp"
#include "Types.hpp"
#include "utilities/ViewUtils.hpp"

// Scream includes
#include "dynamics/homme/homme_dimensions.hpp"
#include "share/util/eamxx_column_ops.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"

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

void HommeDynamics::compute_horizontal_derivs_for_3d_turbulence_and_leonard ()
{
  using namespace Homme;
  using PF = PhysicsFunctions<DefaultDevice>;

  constexpr int NGP  = HOMMEXX_NP;

  const auto& c      = Context::singleton();
  const auto& state  = c.get<ElementsState>();
  const auto& geom   = c.get<ElementsGeometry>();
  const auto& hvcoord = c.get<HybridVCoord>();
  const auto& ref_fe = c.get<ReferenceElement>();
  const auto& tl     = c.get<TimeLevel>();

  const int nelem       = m_dyn_grid->get_num_local_dofs() / (NGP*NGP);
  const int n0          = tl.n0;
  const int n0_qdp      = tl.n0_qdp;
  const int qc_idx      = m_qc_idx;
  const Real ptop       = hvcoord.ps0 * hvcoord.hybrid_ai0;
  const auto& grad_Ux_field = m_helper_fields.at("grad_Ux_dyn");
  const auto& grad_Ux_layout = grad_Ux_field.get_header().get_identifier().get_layout();
  const int nlev_scalar = grad_Ux_layout.dims().back();

  const auto w_int_dyn = state.m_w_i;
  const auto vtheta_dp_dyn = state.m_vtheta_dp;
  const auto dp_dyn = state.m_dp3d;
  const auto& tracers = c.get<Tracers>();
  const auto qdp_dyn = tracers.qdp;

  auto grad_Ux_dyn = m_helper_fields.at("grad_Ux_dyn").template get_view<Real*****>();
  auto grad_Uy_dyn = m_helper_fields.at("grad_Uy_dyn").template get_view<Real*****>();
  auto grad_Uz_dyn = m_helper_fields.at("grad_Uz_dyn").template get_view<Real*****>();
  auto wthl_leonard_base_dyn = m_helper_fields.at("wthl_leonard_base_dyn").template get_view<Real****>();
  auto wqt_leonard_base_dyn = m_helper_fields.at("wqt_leonard_base_dyn").template get_view<Real****>();
  auto uw_leonard_base_dyn = m_helper_fields.at("uw_leonard_base_dyn").template get_view<Real****>();
  auto vw_leonard_base_dyn = m_helper_fields.at("vw_leonard_base_dyn").template get_view<Real****>();

  const auto dvv              = ref_fe.get_deriv();
  const auto dinv             = geom.m_dinv;
  const auto vec_sph2cart     = geom.m_vec_sph2cart;
  const Real scale_factor_inv = 1.0 / geom.m_scale_factor;

  using TeamPolicy = Kokkos::TeamPolicy<KT::ExeSpace>;
  using MemberType = typename TeamPolicy::member_type;
  const int ncols = nelem*NGP*NGP;
  const TeamPolicy policy(ncols, Kokkos::AUTO());

  Kokkos::parallel_for(
      "compute_horizontal_derivs_for_3d_turbulence_and_leonard",
      policy,
      KOKKOS_LAMBDA (const MemberType& team) {

    const int ie  =  team.league_rank() / (NGP*NGP);
    const int igp = (team.league_rank() / NGP) % NGP;
    const int jgp =  team.league_rank() % NGP;

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_scalar), [&] (const int ilev) {
      Real dsdx_ux = 0;
      Real dsdy_ux = 0;
      Real dsdx_uy = 0;
      Real dsdy_uy = 0;
      Real dsdx_uz = 0;
      Real dsdy_uz = 0;
      Real dsdx_theta_l = 0;
      Real dsdy_theta_l = 0;
      Real dsdx_qt = 0;
      Real dsdy_qt = 0;
      Real dsdx_u = 0;
      Real dsdy_u = 0;
      Real dsdx_v = 0;
      Real dsdy_v = 0;

      for (int kgp = 0; kgp < NGP; ++kgp) {
        // The horizontal stencil uses interface w, so average the two
        // adjacent interface values onto midpoint levels on the fly.
        const auto w_row_i = Kokkos::subview(w_int_dyn, ie, n0, igp, kgp, Kokkos::ALL());
        const auto w_col_i = Kokkos::subview(w_int_dyn, ie, n0, kgp, jgp, Kokkos::ALL());
        const auto w_row_i_real = Homme::viewAsReal(w_row_i);
        const auto w_col_i_real = Homme::viewAsReal(w_col_i);

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
        const auto vtheta_dp_row =
            Homme::viewAsReal(Kokkos::subview(vtheta_dp_dyn, ie, n0, igp, kgp, Kokkos::ALL()));
        const auto dp_row =
            Homme::viewAsReal(Kokkos::subview(dp_dyn, ie, n0, igp, kgp, Kokkos::ALL()));
        const auto qvdp_row =
            Homme::viewAsReal(Kokkos::subview(qdp_dyn, ie, n0_qdp, 0, igp, kgp, Kokkos::ALL()));
        const auto vtheta_dp_col =
            Homme::viewAsReal(Kokkos::subview(vtheta_dp_dyn, ie, n0, kgp, jgp, Kokkos::ALL()));
        const auto dp_col =
            Homme::viewAsReal(Kokkos::subview(dp_dyn, ie, n0, kgp, jgp, Kokkos::ALL()));
        const auto qvdp_col =
            Homme::viewAsReal(Kokkos::subview(qdp_dyn, ie, n0_qdp, 0, kgp, jgp, Kokkos::ALL()));

        const Real u_row = u_row_view(ilev);
        const Real v_row = v_row_view(ilev);
        const Real w_row = 0.5 * (w_row_i_real(ilev) + w_row_i_real(ilev + 1));

        const Real u_col = u_col_view(ilev);
        const Real v_col = v_col_view(ilev);
        const Real w_col = 0.5 * (w_col_i_real(ilev) + w_col_i_real(ilev + 1));
        const Real qv_row = qvdp_row(ilev) / dp_row(ilev);
        const Real qv_col = qvdp_col(ilev) / dp_col(ilev);
        Real qc_row = 0;
        Real qc_col = 0;
        if (qc_idx >= 0) {
          const auto qcdp_row =
              Homme::viewAsReal(Kokkos::subview(qdp_dyn, ie, n0_qdp, qc_idx, igp, kgp, Kokkos::ALL()));
          const auto qcdp_col =
              Homme::viewAsReal(Kokkos::subview(qdp_dyn, ie, n0_qdp, qc_idx, kgp, jgp, Kokkos::ALL()));
          qc_row = qcdp_row(ilev) / dp_row(ilev);
          qc_col = qcdp_col(ilev) / dp_col(ilev);
        }
        const Real theta_row =
            PF::calculate_temperature_from_virtual_temperature(vtheta_dp_row(ilev) / dp_row(ilev), qv_row);
        const Real theta_col =
            PF::calculate_temperature_from_virtual_temperature(vtheta_dp_col(ilev) / dp_col(ilev), qv_col);
        Real p_int_top_row = ptop;
        Real p_int_top_col = ptop;
        for (int k = 0; k < ilev; ++k) {
          p_int_top_row += dp_row(k);
          p_int_top_col += dp_col(k);
        }
        const Real T_row = PF::calculate_T_from_theta(theta_row, p_int_top_row + 0.5*dp_row(ilev));
        const Real T_col = PF::calculate_T_from_theta(theta_col, p_int_top_col + 0.5*dp_col(ilev));
        const Real theta_l_row = PF::calculate_thetal_from_theta(theta_row, T_row, qc_row);
        const Real theta_l_col = PF::calculate_thetal_from_theta(theta_col, T_col, qc_col);
        const Real qt_row = qv_row + qc_row;
        const Real qt_col = qv_col + qc_col;

        dsdx_ux += dvv(jgp,kgp) * local_to_cart_component(row_x, u_row, v_row, w_row);
        dsdy_ux += dvv(igp,kgp) * local_to_cart_component(col_x, u_col, v_col, w_col);
        dsdx_u += dvv(jgp,kgp) * u_row;
        dsdy_u += dvv(igp,kgp) * u_col;

        dsdx_uy += dvv(jgp,kgp) * local_to_cart_component(row_y, u_row, v_row, w_row);
        dsdy_uy += dvv(igp,kgp) * local_to_cart_component(col_y, u_col, v_col, w_col);
        dsdx_v += dvv(jgp,kgp) * v_row;
        dsdy_v += dvv(igp,kgp) * v_col;

        dsdx_uz += dvv(jgp,kgp) * local_to_cart_component(row_z, u_row, v_row, w_row);
        dsdy_uz += dvv(igp,kgp) * local_to_cart_component(col_z, u_col, v_col, w_col);
        dsdx_theta_l += dvv(jgp,kgp) * theta_l_row;
        dsdy_theta_l += dvv(igp,kgp) * theta_l_col;
        dsdx_qt += dvv(jgp,kgp) * qt_row;
        dsdy_qt += dvv(igp,kgp) * qt_col;
      }

      const auto dinv_ij = Kokkos::subview(dinv, ie, Kokkos::ALL(), Kokkos::ALL(), igp, jgp);
      grad_Ux_dyn(ie,0,igp,jgp,ilev) = (dinv_ij(0,0) * dsdx_ux + dinv_ij(0,1) * dsdy_ux) * scale_factor_inv;
      grad_Uy_dyn(ie,0,igp,jgp,ilev) = (dinv_ij(0,0) * dsdx_uy + dinv_ij(0,1) * dsdy_uy) * scale_factor_inv;
      grad_Uz_dyn(ie,0,igp,jgp,ilev) = (dinv_ij(0,0) * dsdx_uz + dinv_ij(0,1) * dsdy_uz) * scale_factor_inv;
      const Real grad_thl_0 = (dinv_ij(0,0) * dsdx_theta_l + dinv_ij(0,1) * dsdy_theta_l) * scale_factor_inv;
      const Real grad_qt_0 = (dinv_ij(0,0) * dsdx_qt + dinv_ij(0,1) * dsdy_qt) * scale_factor_inv;
      const Real grad_u_0 = (dinv_ij(0,0) * dsdx_u + dinv_ij(0,1) * dsdy_u) * scale_factor_inv;
      const Real grad_v_0 = (dinv_ij(0,0) * dsdx_v + dinv_ij(0,1) * dsdy_v) * scale_factor_inv;

      grad_Ux_dyn(ie,1,igp,jgp,ilev) = (dinv_ij(1,0) * dsdx_ux + dinv_ij(1,1) * dsdy_ux) * scale_factor_inv;
      grad_Uy_dyn(ie,1,igp,jgp,ilev) = (dinv_ij(1,0) * dsdx_uy + dinv_ij(1,1) * dsdy_uy) * scale_factor_inv;
      grad_Uz_dyn(ie,1,igp,jgp,ilev) = (dinv_ij(1,0) * dsdx_uz + dinv_ij(1,1) * dsdy_uz) * scale_factor_inv;
      const Real grad_thl_1 = (dinv_ij(1,0) * dsdx_theta_l + dinv_ij(1,1) * dsdy_theta_l) * scale_factor_inv;
      const Real grad_qt_1 = (dinv_ij(1,0) * dsdx_qt + dinv_ij(1,1) * dsdy_qt) * scale_factor_inv;
      const Real grad_u_1 = (dinv_ij(1,0) * dsdx_u + dinv_ij(1,1) * dsdy_u) * scale_factor_inv;
      const Real grad_v_1 = (dinv_ij(1,0) * dsdx_v + dinv_ij(1,1) * dsdy_v) * scale_factor_inv;

      const Real b2_0 = vec_sph2cart(ie, 2, 0, igp, jgp);
      const Real b2_1 = vec_sph2cart(ie, 2, 1, igp, jgp);
      const Real b2_2 = vec_sph2cart(ie, 2, 2, igp, jgp);
      const Real dw_dloc0 = b2_0 * grad_Ux_dyn(ie,0,igp,jgp,ilev)
                          + b2_1 * grad_Uy_dyn(ie,0,igp,jgp,ilev)
                          + b2_2 * grad_Uz_dyn(ie,0,igp,jgp,ilev);
      const Real dw_dloc1 = b2_0 * grad_Ux_dyn(ie,1,igp,jgp,ilev)
                          + b2_1 * grad_Uy_dyn(ie,1,igp,jgp,ilev)
                          + b2_2 * grad_Uz_dyn(ie,1,igp,jgp,ilev);
      wthl_leonard_base_dyn(ie,igp,jgp,ilev) = dw_dloc0 * grad_thl_0 + dw_dloc1 * grad_thl_1;
      wqt_leonard_base_dyn(ie,igp,jgp,ilev) = dw_dloc0 * grad_qt_0 + dw_dloc1 * grad_qt_1;
      uw_leonard_base_dyn(ie,igp,jgp,ilev) = dw_dloc0 * grad_u_0 + dw_dloc1 * grad_u_1;
      vw_leonard_base_dyn(ie,igp,jgp,ilev) = dw_dloc0 * grad_v_0 + dw_dloc1 * grad_v_1;
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

  auto& shear_components_field = m_helper_fields.at("shear_strain3d_components_dyn");
  const auto& shear_components_layout = shear_components_field.get_header().get_identifier().get_layout();
  auto shear_components_dyn = shear_components_field.template get_view<Real*****>();

  const auto vec_sph2cart = geom.m_vec_sph2cart;

  const int nlev_scalar = shear_components_layout.dims().back();

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
