#include "gw_functions.hpp"
#include "eamxx_gw_process_interface.hpp"

#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>

#include <array>

namespace scream
{

/*------------------------------------------------------------------------------------------------*/
GWDrag::GWDrag(const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params) {
  // Nothing to do here
}
/*------------------------------------------------------------------------------------------------*/
void GWDrag::create_requests() {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  constexpr int pack_size = Spack::n;

  // retrieve runtime options
  GWF::s_common_init.load_runtime_options(m_params);
  GWF::s_convect_init.load_runtime_options(m_params);
  GWF::s_front_init.load_runtime_options(m_params);

  m_grid = m_grids_manager->get_grid("physics");

  const auto& grid_name = m_grid->name();
  const auto layout     = m_grid->get_3d_scalar_layout(true);

  // retrieve local grid parameters
  m_ncol = m_grid->get_num_local_dofs();
  m_nlev = m_grid->get_num_vertical_levels();

  const auto nondim = Units::nondimensional();
  const auto m2     = pow(m,2);
  const auto s2     = pow(s,2);
  const auto K2     = pow(K,2);

  FieldLayout scalar2d     = m_grid->get_2d_scalar_layout();        // 2D variables
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(true);    // 3D variables at mid-levels
  FieldLayout scalar3d_int = m_grid->get_3d_scalar_layout(false);   // 3D variables at interfaces
  FieldLayout vector3d_mid = m_grid->get_3d_vector_layout(true,2);  // horiz_wind field

  // Input variables
  add_field<Required>("p_mid",                scalar3d_mid, Pa,     grid_name, pack_size);
  add_field<Required>("p_int",                scalar3d_int, Pa,     grid_name, pack_size);
  add_field<Required>("pseudo_density",       scalar3d_mid, Pa,     grid_name, pack_size);
  add_field<Required>("phis",                 scalar2d    , m2/s2,  grid_name);
  add_field<Required>("omega",                scalar3d_mid, Pa/s,   grid_name, pack_size);
  add_field<Required>("landfrac",             scalar2d    , nondim, grid_name);
  add_field<Required>("sgh",                  scalar2d    , nondim, grid_name);

  // Input/Output variables
  add_field <Updated>("T_mid",                scalar3d_mid, K,      grid_name, pack_size);
  add_tracer<Updated>("qv",                   m_grid,       kg/kg,             pack_size);
  add_tracer<Updated>("qc",                   m_grid,       kg/kg,             pack_size);
  add_tracer<Updated>("qi",                   m_grid,       kg/kg,             pack_size);
  add_field <Updated>("horiz_winds",          vector3d_mid, m/s,    grid_name, pack_size);

  // // Output variables
  // add_field <Updated>("precip_liq_surf_mass", scalar2d,     kg/m2,  grid_name, "ACCUMULATED");
  // add_field <Updated>("precip_ice_surf_mass", scalar2d,     kg/m2,  grid_name, "ACCUMULATED");

  // Diagnostic Outputs
  add_field<Computed>("gw_activity",          scalar2d,     nondim, grid_name);
  add_field<Computed>("gw_T_mid_tend",        scalar3d_mid, K/s,    grid_name, pack_size);
  add_field<Computed>("gw_qv_tend",           scalar3d_mid, kg/kg/s,grid_name, pack_size);
  add_field<Computed>("gw_u_tend",            scalar3d_mid, m/s/s,  grid_name, pack_size);
  add_field<Computed>("gw_v_tend",            scalar3d_mid, m/s/s,  grid_name, pack_size);
}

/*------------------------------------------------------------------------------------------------*/
void GWDrag::initialize_impl (const RunType) {
  // Set property checks for fields in this process
  add_postcondition_check<Interval>(get_field_out("T_mid"),       m_grid,100.0,400.0,false);
  add_postcondition_check<Interval>(get_field_out("horiz_winds"), m_grid,-200.0, 200.0,false);

  // // Set phase speeds
  // cref = (/ (dc * l, l = -pgwv, pgwv) /)

  // GWF::gw_common_init( m_nlev,
  //                      pgwv,
  //                      dc,
  //                      cref_in,
  //                      orographic_only_in,
  //                      do_molec_diff_in,
  //                      tau_0_ubc_in,
  //                      nbot_molec_in,
  //                      ktop_in,
  //                      kbotbg_in,
  //                      fcrit2_in,
  //                      kwv_in,
  //                      alpha_in );
}

/*------------------------------------------------------------------------------------------------*/
void GWDrag::run_impl (const double dt) {
  const int nlev_mid_packs = ekat::npack<Spack>(m_nlev);
  //----------------------------------------------------------------------------
  // get fields

  // variables not updated by GWD
  const auto& phis        = get_field_in("phis")          .get_view<const Real*>();
  const auto& p_mid       = get_field_in("p_mid")         .get_view<const Spack**>();
  const auto& p_int       = get_field_in("p_int")         .get_view<const Spack**>();
  const auto& p_del       = get_field_in("pseudo_density").get_view<const Spack**>();
  const auto& omega       = get_field_in("omega")         .get_view<const Spack**>();
  const auto& landfrac    = get_field_in("landfrac")      .get_view<const Real*>();
  const auto& sgh         = get_field_in("sgh")           .get_view<const Real*>();

  // variables updated by GWD
  const auto& T_mid       = get_field_out("T_mid")        .get_view<Spack**>();
  const auto& qv          = get_field_out("qv")           .get_view<Spack**>();
  const auto& qc          = get_field_out("qc")           .get_view<Spack**>();
  const auto& qi          = get_field_out("qi")           .get_view<Spack**>();
  const auto& hwinds_fld  = get_field_out("horiz_winds");
  const auto& uwind       = hwinds_fld.get_component(0)   .get_view<Spack**>();
  const auto& vwind       = hwinds_fld.get_component(1)   .get_view<Spack**>();
  //----------------------------------------------------------------------------
  // // calculate altitude on interfaces (z_int) and mid-points (z_mid)

  // // create temporaries to avoid "Implicit capture" warning
  // const auto loc_p_mid = p_mid;
  // const auto loc_p_del = p_del;
  // const auto loc_T_mid = T_mid;
  // const auto loc_qv    = qv;
  // auto loc_z_mid = m_buffer.z_mid;
  // auto loc_z_del = m_buffer.z_del;
  // auto loc_z_int = m_buffer.z_int;
  // auto loc_nlev = m_nlev;

  // Kokkos::parallel_for(scan_policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
  //   const int i = team.league_rank();
  //   const auto p_mid_i = ekat::subview(loc_p_mid, i);
  //   const auto p_del_i = ekat::subview(loc_p_del, i);
  //   const auto T_mid_i = ekat::subview(loc_T_mid, i);
  //   const auto qv_i    = ekat::subview(loc_qv,    i);
  //   auto z_mid_i = ekat::subview(loc_z_mid, i);
  //   auto z_del_i = ekat::subview(loc_z_del, i);
  //   auto z_int_i = ekat::subview(loc_z_int, i);
  //   auto z_surf = 0.0; // z_mid & z_int are altitude above the surface
  //   PF::calculate_dz(team, p_del_i, p_mid_i, T_mid_i, qv_i, z_del_i);
  //   team.team_barrier();
  //   PF::calculate_z_int(team, loc_nlev, z_del_i, z_surf, z_int_i);
  //   team.team_barrier();
  //   PF::calculate_z_mid(team, loc_nlev, z_int_i, z_mid_i);
  //   team.team_barrier();
  // });
  //----------------------------------------------------------------------------

  // // Calculate local molecular diffusivity
  // if (do_molec_diff) {
  // }

  const auto loc_common_init = GWF::s_common_init;
  const auto loc_convect_init = GWF::s_convect_init;
  const auto loc_front_init = GWF::s_front_init;

  // Convective gravity waves (Beres scheme)
  if (loc_common_init.use_gw_convect) {

    // // Determine wave sources
    // GWF::gw_beres_src();

    // // Solve for the drag profile
    // GWF::gw_drag_prof();

    // add the diffusion coefficients
    // do k = 0, pver
    //   egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
    // end do

    // Store constituents tendencies
    // do m=1, pcnst
    //    do k = 1, pver
    //       ptend%q(:ncol,k,m) = qtgw(:,k,m)
    //    end do
    // end do

    // ! add the momentum tendencies to the output tendency arrays
    // do k = 1, pver
    //    ptend%u(:ncol,k) = utgw(:,k)
    //    ptend%v(:ncol,k) = vtgw(:,k)
    //    ptend%s(:ncol,k) = ttgw(:,k)
    // end do

    // // Momentum & energy conservation
    // GWF::momentum_energy_conservation();

  }

  // Frontally generated gravity waves
  if (loc_common_init.use_gw_frontal) {
    // GWF::gw_cm_src();
    // GWF::gw_drag_prof();
    // GWF::momentum_energy_conservation();
  }

  // Orographic stationary gravity waves
  if (loc_common_init.use_gw_orographic) {
    
    // // Determine the orographic wave source
    // GWF::gw_oro_src(team,
    //                 loc_common_init,
    //                 m_ncol,
    //                 uwind,
    //                 vwind,
    //                 T_mid,
    //                 sgh,
    //                 p_mid,
    //                 p_int,
    //                 p_del,
    //                 z_mid,
    //                 nm,
    //                 src_level,
    //                 tend_level,
    //                 tau,
    //                 ubm,
    //                 ubi,
    //                 xv,
    //                 yv,
    //                 c );

    // Solve for the drag profile with orographic sources
    // GWF::gw_drag_prof();

    // // GW energy fixer
    // do k = 1, pver
    //    utgw(:,k) = utgw(:,k) * cam_in%landfrac(:ncol)
    //    ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
    //    vtgw(:,k) = vtgw(:,k) * cam_in%landfrac(:ncol)
    //    ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
    //    ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
    // enddo

    // dE = 0.0
    // do k = 1, pver
    //    dE(:ncol) = dE(:ncol) &
    //              - dpm(:ncol,k)*(ptend%u(:ncol,k) * (u(:ncol,k)+ptend%u(:ncol,k)*0.5_r8*dt) &
    //                             +ptend%v(:ncol,k) * (v(:ncol,k)+ptend%v(:ncol,k)*0.5_r8*dt) &
    //                             +ptend%s(:ncol,k) )
    // enddo
    // dE(:ncol)=dE(:ncol) / (pint(:ncol,pver+1) - pint(:ncol,1))

    // do k = 1, pver
    //    ptend%s(:ncol,k) = ptend%s(:ncol,k) + dE(:ncol)
    //    ttgw(:ncol,k) = ( ttgw(:ncol,k) + dE(:ncol) ) / cpairv(:ncol, k, lchnk)
    // enddo

    // // add orographic constituent tendencies to total constituent tendencies
    // do m = 1, pcnst
    //     do k = 1, pver
    //        ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
    //     end do
    // end do

  }

  // Convert the tendencies for the dry constituents to dry air basis.
  // do m = 1, pcnst
  //    if (cnst_type(m).eq.'dry') then
  //       do k = 1, pver
  //          do i = 1, ncol
  //             ptend%q(i,k,m) = ptend%q(i,k,m)*state1%pdel(i,k)/state1%pdeldry(i,k)
  //          end do
  //       end do
  //    end if
  // end do

}
/*------------------------------------------------------------------------------------------------*/
size_t GWDrag::requested_buffer_size_in_bytes() const
{
  const int nlev_mid_packs = ekat::npack<Spack>(m_nlev);
  const int nlev_int_packs = ekat::npack<Spack>(m_nlev+1);
  size_t gw_buffer_size = 0;

  gw_buffer_size += Buffer::num_2d_mid_views*m_ncols*nlev_mid_packs*sizeof(Pack);
  gw_buffer_size += Buffer::num_2d_int_views*m_ncols*nlev_int_packs*sizeof(Pack);

  return gw_buffer_size;
}
/*------------------------------------------------------------------------------------------------*/
void GWDrag::init_buffers(const ATMBufferManager &buffer_manager)
{
  auto buffer_chk = ( buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(buffer_chk,"Error! Buffers size not sufficient.\n");
  //----------------------------------------------------------------------------
  Pack* mem = reinterpret_cast<Pack*>(buffer_manager.get_memory());
  const int nlev_mid_packs = ekat::npack<Spack>(m_nlev);
  const int nlev_int_packs = ekat::npack<Spack>(m_nlev+1);
  //----------------------------------------------------------------------------
  uview_2d* buffer_mid_view_ptrs[Buffer::num_2d_midpoint_views] = {
    &m_buffer.z_del,
    &m_buffer.z_mid
  };
  for (int i=0; i<Buffer::num_2d_midpoint_views; ++i) {
    *buffer_mid_view_ptrs[i] = uview_2d(mem, m_ncols, nlev_packs);
    mem += buffer_mid_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  uview_2d* buffer_int_view_ptrs[Buffer::num_2d_interface_views] = {
    &m_buffer.z_int
  };
  for (int i=0; i<Buffer::num_2d_interface_views; ++i) {
    *buffer_int_view_ptrs[i] = uview_2d(mem, m_ncols, nlevi_packs);
    mem += buffer_int_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  size_t used_mem = (reinterpret_cast<Real*>(mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem == requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for TurbulentMountainStress.");
}
/*------------------------------------------------------------------------------------------------*/
void GWDrag::finalize_impl ()
{
  // placeholder for final cleanup
}
/*------------------------------------------------------------------------------------------------*/
} // namespace scream
