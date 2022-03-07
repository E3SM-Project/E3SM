#include "control/surface_coupling.hpp"

#include "share/util/scream_common_physics_functions.hpp"

namespace scream {
namespace control {

void SurfaceCoupling::do_export (const bool init_phase)
{
  if (m_num_scream_exports==0) {
    return;
  }

  using KT = KokkosTypes<device_type>;
  using policy_type = KT::RangePolicy;
  using PF = PhysicsFunctions<device_type>;
  using C = scream::physics::Constants<Real>;

  // For each export fields that is not trivially exist in the field
  // manager (see Case 2 in register_export()), calculate correct data
  // values.
  const bool scream_ad_run =
      (m_field_mgr->has_field("qv") && m_field_mgr->has_field("T_mid") &&
       m_field_mgr->has_field("p_mid") && m_field_mgr->has_field("pseudo_density") &&
       m_field_mgr->has_field("precip_liq_surf") && m_field_mgr->has_field("precip_ice_surf"));
  if (scream_ad_run) {
    const int last_entry = m_num_levs-1;
    const auto& qv              = m_field_mgr->get_field("qv").get_view<const Real**>();
    const auto& T_mid           = m_field_mgr->get_field("T_mid").get_view<const Real**>();
    const auto& p_mid           = m_field_mgr->get_field("p_mid").get_view<const Real**>();
    const auto& pseudo_density  = m_field_mgr->get_field("pseudo_density").get_view<const Real**>();
    const auto& precip_liq_surf = m_field_mgr->get_field("precip_liq_surf").get_view<const Real*>();
    const auto& precip_ice_surf = m_field_mgr->get_field("precip_ice_surf").get_view<const Real*>();
    const auto l_dz             = dz;
    const auto l_z_int          = z_int;
    const auto l_z_mid          = z_mid;
    const auto l_Sa_z           = Sa_z;
    const auto l_Sa_ptem        = Sa_ptem;
    const auto l_Sa_dens        = Sa_dens;
    const auto l_Faxa_rainl     = Faxa_rainl;
    const auto l_Faxa_snowl     = Faxa_snowl;

    // Local copy, to deal with CUDA's handling of *this.
    const int num_levs = m_num_levs;

    const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(m_num_cols, num_levs);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) {
      const int i = team.league_rank();

      const auto qv_i              = ekat::subview(qv, i);
      const auto T_mid_i           = ekat::subview(T_mid, i);
      const auto p_mid_i           = ekat::subview(p_mid, i);
      const auto pseudo_density_i  = ekat::subview(pseudo_density, i);
      const auto dz_i              = ekat::subview(l_dz, i);
      const auto z_int_i           = ekat::subview(l_z_int, i);
      const auto z_mid_i           = ekat::subview(l_z_mid, i);

      // Compute vertical layer thickness
      PF::calculate_dz(team, pseudo_density_i, p_mid_i, T_mid_i, qv_i, dz_i);
      team.team_barrier();

      // Compute vertical layer heights. Use z_int(nlevs) = z_surf = 0.0.
      const Real z_surf = 0.0;
      PF::calculate_z_int(team, num_levs, dz_i, z_surf, z_int_i);
      team.team_barrier();
      PF::calculate_z_mid(team, num_levs, z_int_i, z_mid_i);
      team.team_barrier();

      l_Sa_z(i)       = z_mid_i(last_entry);
      l_Sa_ptem(i)    = PF::calculate_theta_from_T(T_mid_i(last_entry), p_mid_i(last_entry));
      l_Sa_dens(i)    = PF::calculate_density(pseudo_density_i(last_entry), dz_i(last_entry));
      l_Faxa_rainl(i) = precip_liq_surf(i)*C::RHO_H2O;
      l_Faxa_snowl(i) = precip_ice_surf(i)*C::RHO_H2O; //p3_ice_sed_impl.hpp uses INV_RHO_H2O
    });
  }

  // Local copies, to deal with CUDA's handling of *this.
  const auto scream_exports = m_scream_exports_dev;
  const auto cpl_exports_view_d = m_cpl_exports_view_d;
  const int num_cols = m_num_cols;

  // Pack the fields
  auto pack_policy   = policy_type (0,m_num_scream_exports*num_cols);
  Kokkos::parallel_for(pack_policy, KOKKOS_LAMBDA(const int& i) {
    const int ifield = i / num_cols;
    const int icol   = i % num_cols;
    const auto& info = scream_exports(ifield);
    const auto offset = icol*info.col_stride + info.col_offset;

    // during the initial export, some fields may need to be skipped
    // since values have not been computed inside SCREAM at the time
    bool do_export = (!init_phase || info.do_initial_export);
    if (do_export) {
      cpl_exports_view_d(icol,info.cpl_idx) = info.data[offset];
    }
  });

  // Deep copy fields from device to cpl host array
  Kokkos::deep_copy(m_cpl_exports_view_h,m_cpl_exports_view_d);
}

} // namespace control
} // namespace scream
