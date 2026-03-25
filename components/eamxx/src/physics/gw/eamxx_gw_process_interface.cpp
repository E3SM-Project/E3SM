#include "gw_functions.hpp"
#include "eamxx_gw_process_interface.hpp"

#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

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
  constexpr int pack_size = Pack::n;

  m_grid = m_grids_manager->get_grid("physics");

  const auto& grid_name = m_grid->name();
  const auto layout     = m_grid->get_3d_scalar_layout(true);

  // retrieve local grid parameters
  m_ncol = m_grid->get_num_local_dofs();
  m_nlev = m_grid->get_num_vertical_levels();
  m_npgw = m_params.get<int>("pgwv") * 2 + 1;

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

  // defaults for gw_common_init()
  bool do_molec_diff_default = false; // Flag for molecular diffusion
  int nbot_molec_default = 0;         // bottom level for molecular diffusion
  int ktop_default = 0;               // Top level for gravity waves.
  Real kwv_default = 6.28e-5;         // Effective horizontal wave number (100 km wavelength)

  // calculate interface reference pressures
  const auto hyai = m_grid->get_geometry_data("hyai").get_view<const Real*>();
  const auto hybi = m_grid->get_geometry_data("hybi").get_view<const Real*>();
  Kokkos::View<Real*, Kokkos::HostSpace> pref_int("pref_int", hyai.size());
  // Kokkos::parallel_for(Kokkos::RangePolicy<>(0, m_nlev+1), KOKKOS_LAMBDA (const int k) {
  Kokkos::parallel_for("calculate_pref_int", 
    Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, hyai.size()), 
    KOKKOS_LAMBDA (const int k) {
      pref_int(k) = PC::P0.value * hyai(k) + PC::P0.value * hybi(k);
  });
  Kokkos::fence();

  GWF::gw_common_init( m_params, 
                       m_nlev,
                       pref_int,
                       do_molec_diff_default,
                       nbot_molec_default,
                       ktop_default,
                       kwv_default);

  std::string gw_drag_file = m_params.get<std::string>("gw_drag_file");
  scorpio::register_file(gw_drag_file,scorpio::FileMode::Read);
  const int PS_dim_size = scorpio::get_dimlen(gw_drag_file, "PS"); // Phase Speed [m/s]
  const int MW_dim_size = scorpio::get_dimlen(gw_drag_file, "MW"); // Mean Wind in Heating [m/s]
  const int HD_dim_size = scorpio::get_dimlen(gw_drag_file, "HD"); // Heating Depth [km]
  Kokkos::View<Real***, Kokkos::HostSpace> mfcc_host("mfcc_host", HD_dim_size, MW_dim_size, PS_dim_size);
  scorpio::read_var(gw_drag_file,"mfcc",mfcc_host.data());
  scorpio::release_file(gw_drag_file);

  GWF::gw_convect_init( m_params, mfcc_host );

  GWF::gw_front_init( m_params, pref_int );

  // Set property checks for fields in this process
  using Interval = FieldWithinIntervalCheck;
  using LowerBound = FieldLowerBoundCheck;
  add_postcondition_check<Interval>(get_field_out("T_mid"),       m_grid,100.0,400.0,false);
  add_postcondition_check<Interval>(get_field_out("horiz_winds"), m_grid,-200.0, 200.0,false);
}

/*------------------------------------------------------------------------------------------------*/
void GWDrag::run_impl (const double dt) {
  using PC = scream::physics::Constants<Real>;
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
  const int nlev_mid_packs = ekat::npack<Pack>(m_nlev);
  const int nlev_int_packs = ekat::npack<Pack>(m_nlev+1);
  const auto team_policy = TPF::get_default_team_policy(m_ncol, nlev_mid_packs);
  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(m_ncol, nlev_mid_packs);
  //----------------------------------------------------------------------------
  // get fields

  // variables not updated by GWD
  const auto& phis        = get_field_in("phis")          .get_view<const Real*>();
  const auto& p_mid       = get_field_in("p_mid")         .get_view<const Pack**>();
  const auto& p_int       = get_field_in("p_int")         .get_view<const Pack**>();
  const auto& p_del       = get_field_in("pseudo_density").get_view<const Pack**>();
  const auto& omega       = get_field_in("omega")         .get_view<const Pack**>();
  const auto& landfrac    = get_field_in("landfrac")      .get_view<const Real*>();
  const auto& sgh         = get_field_in("sgh")           .get_view<const Real*>();

  // variables updated by GWD
  const auto& T_mid       = get_field_out("T_mid")        .get_view<Pack**>();
  const auto& qv          = get_field_out("qv")           .get_view<Pack**>();
  const auto& qc          = get_field_out("qc")           .get_view<Pack**>();
  const auto& qi          = get_field_out("qi")           .get_view<Pack**>();
  const auto& hwinds_fld  = get_field_out("horiz_winds");
  const auto& uwind       = hwinds_fld.get_component(0)   .get_view<Pack**>();
  const auto& vwind       = hwinds_fld.get_component(1)   .get_view<Pack**>();
  //----------------------------------------------------------------------------
  // create local temporaries to avoid "Implicit capture" warning
  const auto loc_common_init = GWF::s_common_init;
  const auto loc_convect_init = GWF::s_convect_init;
  const auto loc_front_init = GWF::s_front_init;

  const auto loc_p_mid = p_mid;
  const auto loc_p_int = p_int;
  const auto loc_p_del = p_del;
  const auto loc_T_mid = T_mid;
  const auto loc_qv    = qv;

  auto loc_z_mid       = m_buffer.z_mid;
  auto loc_z_del       = m_buffer.z_del;
  auto loc_z_int       = m_buffer.z_int;
  auto loc_T_int       = m_buffer.T_int;
  auto loc_N_mid       = m_buffer.N_mid;
  auto loc_N_int       = m_buffer.N_int;
  auto loc_rho_int     = m_buffer.rho_int;
  //----------------------------------------------------------------------------
  // calculate altitude on interfaces (z_int) and mid-points (z_mid)
  Kokkos::parallel_for(scan_policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();
    const auto p_mid_i = ekat::subview(loc_p_mid, i);
    const auto p_del_i = ekat::subview(loc_p_del, i);
    const auto T_mid_i = ekat::subview(loc_T_mid, i);
    const auto qv_i = ekat::subview(loc_qv,    i);
    auto z_mid_i = ekat::subview(loc_z_mid, i);
    auto z_del_i = ekat::subview(loc_z_del, i);
    auto z_int_i = ekat::subview(loc_z_int, i);
    auto z_surf = 0.0; // z_mid & z_int are altitude above the surface
    PF::calculate_dz(team, p_del_i, p_mid_i, T_mid_i, qv_i, z_del_i);
    team.team_barrier();
    PF::calculate_z_int(team, m_nlev, z_del_i, z_surf, z_int_i);
    team.team_barrier();
    PF::calculate_z_mid(team, m_nlev, z_int_i, z_mid_i);
    team.team_barrier();
  });
  //----------------------------------------------------------------------------

  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
    const Int i = team.league_rank();
    // Get single-column subviews of all inputs
    // Note: scalarize includes pack padding, so we need to slice to actual sizes
    const auto T_mid_s   = ekat::scalarize(ekat::subview(loc_T_mid,  i));
    const auto p_mid_s   = ekat::scalarize(ekat::subview(loc_p_mid,  i));
    const auto p_int_s   = ekat::scalarize(ekat::subview(loc_p_int,  i));
    const auto T_int_s   = ekat::scalarize(ekat::subview(loc_T_int,  i));
    const auto N_mid_s   = ekat::scalarize(ekat::subview(loc_N_mid,  i));
    const auto N_int_s   = ekat::scalarize(ekat::subview(loc_N_int,  i));
    const auto rho_int_s = ekat::scalarize(ekat::subview(loc_rho_int,i));

    // Slice to exclude padding: mid levels = [0, m_nlev), int levels = [0, m_nlev+1)
    const auto T_mid_i   = ekat::subview(T_mid_s,   Kokkos::pair<int, int>{0, m_nlev});
    const auto p_mid_i   = ekat::subview(p_mid_s,   Kokkos::pair<int, int>{0, m_nlev});
    const auto p_int_i   = ekat::subview(p_int_s,   Kokkos::pair<int, int>{0, m_nlev+1});
    const auto T_int_i   = ekat::subview(T_int_s,   Kokkos::pair<int, int>{0, m_nlev}); // note: interface data w/ m_nlev dimension
    const auto N_mid_i   = ekat::subview(N_mid_s,   Kokkos::pair<int, int>{0, m_nlev});
    const auto N_int_i   = ekat::subview(N_int_s,   Kokkos::pair<int, int>{0, m_nlev+1});
    const auto rho_int_i = ekat::subview(rho_int_s, Kokkos::pair<int, int>{0, m_nlev+1});

    GWF::gw_prof(team, m_nlev, PC::Cpair.value,
                 T_mid_i, p_mid_i, p_int_i,
                 rho_int_i, T_int_i, N_mid_i, N_int_i);

  });

  // // Calculate local molecular diffusivity
  // if (do_molec_diff) {
  //   ???
  // }


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

    Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
      const Int i = team.league_rank();

      // Get single-column subviews of all inputs
      const auto uwind_i = ekat::scalarize(ekat::subview(uwind, i));
      const auto vwind_i = ekat::scalarize(ekat::subview(vwind, i));
      const auto T_mid_i = ekat::scalarize(ekat::subview(T_mid, i));
      const auto p_mid_i = ekat::scalarize(ekat::subview(p_mid, i));
      const auto p_int_i = ekat::scalarize(ekat::subview(p_int, i));
      const auto p_del_i = ekat::scalarize(ekat::subview(p_del, i));
      const auto z_mid_i = ekat::scalarize(ekat::subview(m_buffer.z_mid, i));
      const auto N_mid_i = ekat::scalarize(ekat::subview(m_buffer.N_mid, i));
      const auto tau_i   = ekat::scalarize(ekat::subview(m_buffer.tau, i));
      const auto ubm_i   = ekat::scalarize(ekat::subview(m_buffer.ubm, i));
      const auto ubi_i   = ekat::scalarize(ekat::subview(m_buffer.ubi, i));
      const auto c_i     = ekat::scalarize(ekat::subview(m_buffer.c, i));

      Int src_lev; // level index of gravity wave source
      Int tnd_lev; // lowest level index where tendencies are allowed
      Real xv;      // zonal unit vector of source wind
      Real yv;      // meridional unit vector of source wind

      // Determine the orographic wave source
      Kokkos::printf("DEBUG: calling gw_oro_src\n");
      GWF::gw_oro_src(team, GWF::s_common_init, m_nlev, m_npgw,
                      uwind_i, vwind_i, T_mid_i, sgh(i),
                      p_mid_i, p_int_i, p_del_i, z_mid_i, N_mid_i,
                      src_lev, tnd_lev,
                      tau_i, ubm_i, ubi_i, xv, yv, c_i );
      Kokkos::printf("DEBUG: done\n");

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

    });
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
  const int nlev_mid_packs = ekat::npack<Pack>(m_nlev);
  const int nlev_int_packs = ekat::npack<Pack>(m_nlev+1);
  size_t gw_buffer_size = 0;

  gw_buffer_size += Buffer::num_3d_mid_views*m_ncol*m_npgw*nlev_mid_packs*sizeof(Pack);
  gw_buffer_size += Buffer::num_2d_mid_views*m_ncol*nlev_mid_packs*sizeof(Pack);
  gw_buffer_size += Buffer::num_2d_int_views*m_ncol*nlev_int_packs*sizeof(Pack);
  gw_buffer_size += Buffer::num_2d_pgw_views*m_ncol*m_npgw*sizeof(Pack);

  return gw_buffer_size;
}
/*------------------------------------------------------------------------------------------------*/
void GWDrag::init_buffers(const ATMBufferManager &buffer_manager)
{
  auto buffer_chk = ( buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(buffer_chk,"Error! Buffers size not sufficient.\n");
  //----------------------------------------------------------------------------
  Pack* mem = reinterpret_cast<Pack*>(buffer_manager.get_memory());
  const int nlev_mid_packs = ekat::npack<Pack>(m_nlev);
  const int nlev_int_packs = ekat::npack<Pack>(m_nlev+1);
  //----------------------------------------------------------------------------
  uview_3d* buffer_3d_mid_view_ptrs[Buffer::num_3d_mid_views] = {
    &m_buffer.tau
  };
  for (int i=0; i<Buffer::num_3d_mid_views; ++i) {
    *buffer_3d_mid_view_ptrs[i] = uview_3d(mem, m_ncol, m_npgw, nlev_mid_packs);
    mem += buffer_3d_mid_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  uview_2d* buffer_2d_mid_view_ptrs[Buffer::num_2d_mid_views] = {
    &m_buffer.z_del,
    &m_buffer.z_mid,
    &m_buffer.T_int,
    &m_buffer.N_mid,
    &m_buffer.ubm
  };
  for (int i=0; i<Buffer::num_2d_mid_views; ++i) {
    *buffer_2d_mid_view_ptrs[i] = uview_2d(mem, m_ncol, nlev_mid_packs);
    mem += buffer_2d_mid_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  uview_2d* buffer_2d_int_view_ptrs[Buffer::num_2d_int_views] = {
    &m_buffer.z_int,
    // &m_buffer.T_int,
    &m_buffer.N_int,
    &m_buffer.rho_int,
    &m_buffer.ubi
  };
  for (int i=0; i<Buffer::num_2d_int_views; ++i) {
    *buffer_2d_int_view_ptrs[i] = uview_2d(mem, m_ncol, nlev_int_packs);
    mem += buffer_2d_int_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  uview_2d* buffer_2d_pgw_view_ptrs[Buffer::num_2d_pgw_views] = {
    &m_buffer.c
  };
  for (int i=0; i<Buffer::num_2d_pgw_views; ++i) {
    *buffer_2d_pgw_view_ptrs[i] = uview_2d(mem, m_ncol, m_npgw);
    mem += buffer_2d_pgw_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  size_t used_mem = (reinterpret_cast<Real*>(mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem == requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for GWDrag.");
}
/*------------------------------------------------------------------------------------------------*/
void GWDrag::finalize_impl ()
{
  // placeholder for final cleanup
}
/*------------------------------------------------------------------------------------------------*/
} // namespace scream
