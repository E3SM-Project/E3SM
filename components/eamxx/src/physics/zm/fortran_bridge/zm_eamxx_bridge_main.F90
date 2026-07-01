module zm_eamxx_bridge_main
  !-----------------------------------------------------------------------------
  ! Purpose: This replicates the functionality of zm_conv_intr.F90 to provide
  ! an interface between EAMxx and EAM's fortran Zhang-McFarlane Deeep Cu scheme
  !-----------------------------------------------------------------------------
  use iso_c_binding
  use cam_logfile,   only: iulog
  use shr_sys_mod,   only: shr_sys_flush
  use zm_eamxx_bridge_params, only: masterproc, r8, pver, pverp, top_lev
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! public methods
  public :: zm_eamxx_bridge_init_c
  public :: zm_eamxx_bridge_run_c

!===================================================================================================
#include "eamxx_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif
!===================================================================================================
contains
!===================================================================================================

subroutine zm_eamxx_bridge_init_c( pver_in, limcnv_in, &
                                  trig_dcape_in, trig_ull_in, &
                                  clos_dyn_adj_in, mcsp_enabled_in, &
                                  mcsp_t_coeff_in, mcsp_q_coeff_in, &
                                  mcsp_mom_coeff_in, mcsp_use_full_shear_in ) bind(C)
  use mpi
  use zm_conv,       only: zm_const, zm_param
  use zm_conv_types, only: zm_const_set_for_testing, zm_param_set_for_testing
  use zm_conv_types, only: zm_param_mpi_broadcast, zm_param_print
  use zm_eamxx_bridge_wv_saturation, only: wv_sat_init
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int), value, intent(in) :: pver_in
  integer(kind=c_int), value, intent(in) :: limcnv_in
  logical(kind=c_bool),value, intent(in) :: trig_dcape_in
  logical(kind=c_bool),value, intent(in) :: trig_ull_in
  logical(kind=c_bool),value, intent(in) :: clos_dyn_adj_in
  logical(kind=c_bool),value, intent(in) :: mcsp_enabled_in
  real(kind=c_real),   value, intent(in) :: mcsp_t_coeff_in
  real(kind=c_real),   value, intent(in) :: mcsp_q_coeff_in
  real(kind=c_real),   value, intent(in) :: mcsp_mom_coeff_in
  logical(kind=c_bool),value, intent(in) :: mcsp_use_full_shear_in
  !-----------------------------------------------------------------------------
  ! Local variables
  integer :: mpi_rank, ierror
  !-----------------------------------------------------------------------------
  pver  = pver_in
  pverp = pver+1
  top_lev = 1
  !-----------------------------------------------------------------------------
  ! obtain master process ID
  call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, ierror)
  masterproc = .false.
  if (mpi_rank==0) masterproc = .true.
  !-----------------------------------------------------------------------------
  ! set ZM constants and parameters
  call zm_const_set_for_testing(zm_const)
  call zm_param_set_for_testing(zm_param)
  zm_param%limcnv = limcnv_in ! override testing value when running the fortran bridge
  call zm_param_mpi_broadcast(zm_param)
  ! override some settings
  zm_param%zm_microp       = .false.
  zm_param%old_snow        = .true.
  zm_param%trig_dcape      = trig_dcape_in
  zm_param%trig_ull        = trig_ull_in
  zm_param%clos_dyn_adj    = clos_dyn_adj_in
  zm_param%mcsp_enabled    = mcsp_enabled_in
  zm_param%mcsp_t_coeff       = mcsp_t_coeff_in
  zm_param%mcsp_q_coeff       = mcsp_q_coeff_in
  zm_param%mcsp_mom_coeff     = mcsp_mom_coeff_in
  zm_param%mcsp_use_full_shear = mcsp_use_full_shear_in
  call zm_param_print(zm_param)
  !-----------------------------------------------------------------------------
  call wv_sat_init()
  !-----------------------------------------------------------------------------
  return
end subroutine zm_eamxx_bridge_init_c

!===================================================================================================

subroutine zm_eamxx_bridge_run_c( ncol, dtime, is_first_step, &
                                  state_phis, state_zm, state_zi, &
                                  state_p_mid, state_p_int, state_p_del, &
                                  state_t, state_qv, state_u, state_v, &
                                  state_omega, state_cldfrac, state_pblh, tpert, landfrac, &
                                  t_star_in, q_star_in, &
                                  output_prec, output_snow, output_cape, output_dcape, output_activity, &
                                  output_tend_t, output_tend_q, output_tend_u, output_tend_v, &
                                  output_rain_prod, output_snow_prod, &
                                  output_prec_flux, output_snow_flux, output_mass_flux, &
                                  output_dlf, mcsp_freq, mcsp_shear, zm_depth, &
                                  mcsp_ds_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out, &
                                  evap_ds_out, evap_dq_out ) bind(C)
  use zm_conv,                  only: zm_const, zm_param
  use zm_aero_type,             only: zm_aero_t
  use zm_microphysics_state,    only: zm_microp_st
  use zm_eamxx_bridge_methods,  only: zm_tend_init, zm_physics_update
  use zm_conv,                  only: zm_conv_main, zm_conv_evap
  use zm_conv_mcsp,             only: zm_conv_mcsp_tend
  use zm_transport,             only: zm_transport_momentum
  ! use zm_transport,             only: zm_transport_tracer
  use zm_conv_types, only: zm_param_print, zm_const_print
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int),               value, intent(in   ) :: ncol               ! 01 number of columns on rank
  real(kind=c_real),                 value, intent(in   ) :: dtime              ! 02 time step
  logical(kind=c_bool),              value, intent(in   ) :: is_first_step      ! 03 flag for first step
  real(kind=c_real),  dimension(ncol),      intent(in   ) :: state_phis         ! 04 input state surface geopotential height
  real(kind=c_real),  dimension(ncol,pver), intent(in   ) :: state_zm           ! 05 input state altitude at mid-levels
  real(kind=c_real),  dimension(ncol,pverp),intent(in   ) :: state_zi           ! 06 input state altitude at interfaces
  real(kind=c_real),  dimension(ncol,pver), intent(in   ) :: state_p_mid        ! 07 input state mid-point pressure
  real(kind=c_real),  dimension(ncol,pverp),intent(in   ) :: state_p_int        ! 08 input state interface pressure
  real(kind=c_real),  dimension(ncol,pver), intent(in   ) :: state_p_del        ! 09 input state pressure thickness
  real(kind=c_real),  dimension(ncol,pver), intent(in   ) :: state_t            ! 10 input state temperature
  real(kind=c_real),  dimension(ncol,pver), intent(in   ) :: state_qv           ! 11 input state water vapor
  real(kind=c_real),  dimension(ncol,pver), intent(in   ) :: state_u            ! 12 input state zonal wind
  real(kind=c_real),  dimension(ncol,pver), intent(in   ) :: state_v            ! 13 input state meridional wind
  real(kind=c_real),  dimension(ncol,pver), intent(in   ) :: state_omega        ! 14 input state vertical pressure velocity
  real(kind=c_real),  dimension(ncol,pver), intent(in   ) :: state_cldfrac      ! 15 input state cloud fraction              (cld)
  real(kind=c_real),  dimension(ncol),      intent(in   ) :: state_pblh         ! 16 input planetary boundary layer height   (pblh)
  real(kind=c_real),  dimension(ncol),      intent(in   ) :: tpert              ! 17 input parcel temperature perturbation
  real(kind=c_real),  dimension(ncol),      intent(in   ) :: landfrac           ! 18 land fraction
  real(kind=c_real),  dimension(ncol,pver), intent(inout) :: t_star_in          ! 19 DCAPE T from time step n-1
  real(kind=c_real),  dimension(ncol,pver), intent(inout) :: q_star_in          ! 20 DCAPE q from time step n-1
  real(kind=c_real),  dimension(ncol),      intent(  out) :: output_prec        ! 21 output total precipitation              (prec)
  real(kind=c_real),  dimension(ncol),      intent(  out) :: output_snow        ! 22 output frozen precipitation             (snow)
  real(kind=c_real),  dimension(ncol),      intent(  out) :: output_cape        ! 23 output convective avail. pot. energy    (cape)
  real(kind=c_real),  dimension(ncol),      intent(  out) :: output_dcape       ! 24 output dynamic cape
  integer(kind=c_int),dimension(ncol),      intent(  out) :: output_activity    ! 25 integer deep convection activity flag   (ideep)
  real(kind=c_real),  dimension(ncol,pver), intent(  out) :: output_tend_t      ! 26 output tendency of temperature          (ptend_loc_s)
  real(kind=c_real),  dimension(ncol,pver), intent(  out) :: output_tend_q      ! 27 output tendency of water vapor          (ptend_loc_q)
  real(kind=c_real),  dimension(ncol,pver), intent(  out) :: output_tend_u      ! 28 output tendency of zonal wind           (ptend_loc_u)
  real(kind=c_real),  dimension(ncol,pver), intent(  out) :: output_tend_v      ! 29 output tendency of meridional wind      (ptend_loc_v)
  real(kind=c_real),  dimension(ncol,pver), intent(  out) :: output_rain_prod   ! 30 rain production rate                    (rprd)
  real(kind=c_real),  dimension(ncol,pver), intent(  out) :: output_snow_prod   ! 31 snow production rate                    (sprd)
  real(kind=c_real),  dimension(ncol,pverp),intent(  out) :: output_prec_flux   ! 32 output precip flux at each mid-levels   (flxprec/pflx)
  real(kind=c_real),  dimension(ncol,pverp),intent(  out) :: output_snow_flux   ! 33 output precip flux at each mid-levels   (flxsnow)
  real(kind=c_real),  dimension(ncol,pverp),intent(  out) :: output_mass_flux   ! 34 output convective mass flux--m sub c    (mcon)
  real(kind=c_real),  dimension(ncol,pver), intent(  out) :: output_dlf         ! 35 detrained convective cloud water        (dlf)
  real(kind=c_real),  dimension(ncol),      intent(  out) :: mcsp_freq          ! 36 MCSP diagnostic output
  real(kind=c_real),  dimension(ncol),      intent(  out) :: mcsp_shear         ! 37 MCSP diagnostic output
  real(kind=c_real),  dimension(ncol),      intent(  out) :: zm_depth           ! 38 MCSP diagnostic output
  real(kind=c_real),  dimension(ncol,pver), intent(  out) :: mcsp_ds_out        ! 39 MCSP tendency
  real(kind=c_real),  dimension(ncol,pver), intent(  out) :: mcsp_dq_out        ! 40 MCSP tendency
  real(kind=c_real),  dimension(ncol,pver), intent(  out) :: mcsp_du_out        ! 41 MCSP tendency
  real(kind=c_real),  dimension(ncol,pver), intent(  out) :: mcsp_dv_out        ! 42 MCSP tendency
  real(kind=c_real),  dimension(ncol,pver), intent(  out) :: evap_ds_out        ! 43 zm_conv_evap tendency
  real(kind=c_real),  dimension(ncol,pver), intent(  out) :: evap_dq_out        ! 44 zm_conv_evap tendency
  !-----------------------------------------------------------------------------
  ! Local variables
  integer :: i,k

  logical(kind=c_bool) :: loc_is_first_step

  ! arguments for zm_conv_main - order somewhat consistent with current interface
  integer,  dimension(ncol)      :: jctop          ! output top-of-deep-convection indices
  integer,  dimension(ncol)      :: jcbot          ! output bot-of-deep-convection indices
  real(r8), dimension(ncol,pver) :: zdu            ! detraining mass flux
  real(r8), dimension(ncol,pver) :: mu             ! updraft cloud mass flux
  real(r8), dimension(ncol,pver) :: md             ! downdraft cloud mass flux
  real(r8), dimension(ncol,pver) :: du             ! detrainment in updraft
  real(r8), dimension(ncol,pver) :: eu             ! entrainment in updraft
  real(r8), dimension(ncol,pver) :: ed             ! entrainment in downdraft
  real(r8), dimension(ncol,pver) :: dp             ! layer thickness [mb]
  real(r8), dimension(ncol)      :: dsubcld        ! sub-cloud layer thickness
  integer,  dimension(ncol)      :: jt             ! top level index of convection
  integer,  dimension(ncol)      :: maxg           ! gathered values of maxi
  integer,  dimension(ncol)      :: ideep          ! flag to indicate ZM is active
  integer                        :: lengath        ! number of gathered columns per chunk
  real(r8), dimension(ncol)      :: rliq           ! reserved liquid (not yet in cldliq) for energy integrals
  real(r8), dimension(ncol,pver), target :: local_t_star ! DCAPE T from time step n-1
  real(r8), dimension(ncol,pver), target :: local_q_star ! DCAPE q from time step n-1
  ! real(r8), dimension(ncol)      :: dcape          ! DCAPE cape change
  type(zm_aero_t)                :: aero           ! derived type for aerosol information
  type(zm_microp_st)             :: microp_st      ! ZM microphysics data structure

  real(r8), dimension(ncol,pver) :: state_s        ! dry static energy
  real(r8), dimension(ncol,pver) :: zm_qc          ! convective in-cloud liquid water

  ! local copy of state variables for calling zm_conv_evap()
  real(r8), dimension(ncol,pver) :: local_state_t
  real(r8), dimension(ncol,pver) :: local_state_qv
  real(r8), dimension(ncol,pver) :: local_state_zm
  real(r8), dimension(ncol,pverp):: local_state_zi

  real(r8), dimension(ncol,pver) :: output_tend_s       ! dry static energy tendency used to set output_tend_t

  ! temporary local tendencies for calling zm_conv_evap()
  real(r8), dimension(ncol,pver) :: local_tend_s        ! temporary tendency of dry static energy   (ptend_loc_s)
  real(r8), dimension(ncol,pver) :: local_tend_q        ! temporary tendency of water vapor         (ptend_loc_q)
  real(r8), dimension(ncol,pver) :: local_tend_u        ! temporary tendency of zonal wind
  real(r8), dimension(ncol,pver) :: local_tend_v        ! temporary tendency of meridional wind

  real(r8), dimension(ncol,pver) :: tend_s_snwprd       ! DSE tend from snow production
  real(r8), dimension(ncol,pver) :: tend_s_snwevmlt     ! DSE tend from snow evap/melt
  real(r8), dimension(ncol,pver) :: ntprprd             ! net precip production in layer
  real(r8), dimension(ncol,pver) :: ntsnprd             ! net snow production in layer

  ! used in momentum transport calculations
  real(r8), dimension(ncol,pver,2) :: tx_winds
  real(r8), dimension(ncol,pver,2) :: tx_wind_tend
  real(r8), dimension(ncol,pver,2) :: tx_pguall
  real(r8), dimension(ncol,pver,2) :: tx_pgdall
  real(r8), dimension(ncol,pver,2) :: tx_icwu
  real(r8), dimension(ncol,pver,2) :: tx_icwd

  !-----------------------------------------------------------------------------
  ! initialize various thing

  loc_is_first_step = is_first_step

  !-----------------------------------------------------------------------------
  ! initialize output tendencies - normally done by physics_ptend_init()

  do i = 1,ncol
    output_prec(i) = 0
    output_snow(i) = 0
    output_cape(i) = 0
    output_activity(i) = 0
    do k = 1,pver
      output_tend_t(i,k) = 0
      output_tend_q(i,k) = 0
      output_tend_u(i,k) = 0
      output_tend_v(i,k) = 0
      output_rain_prod(i,k) = 0
      output_snow_prod(i,k) = 0
      output_prec_flux(i,k) = 0
      output_snow_flux(i,k) = 0
      output_mass_flux(i,k) = 0
      output_dlf(i,k) = 0
    end do
  end do

  !-----------------------------------------------------------------------------
  ! populate local copies of state variables for zm_conv_evap()

  do i = 1,ncol
    do k = 1,pver
      local_state_t (i,k) = state_t (i,k)
      local_state_qv(i,k) = state_qv(i,k)
      local_state_zm(i,k) = state_zm(i,k)
      local_state_zi(i,k) = state_zi(i,k)
      local_t_star  (i,k) = t_star_in(i,k)
      local_q_star  (i,k) = q_star_in(i,k)
    end do
  end do

  !-----------------------------------------------------------------------------
  ! Call the primary Zhang-McFarlane convection parameterization

  call zm_conv_main(ncol, ncol, pver, pverp, loc_is_first_step, dtime, &
                    state_t, state_qv, state_omega, &
                    state_p_mid, state_p_int, state_p_del, &
                    state_phis, state_zm, state_zi, state_pblh, &
                    tpert, landfrac, local_t_star, local_q_star, &
                    lengath, ideep, maxg, jctop, jcbot, jt, &
                    output_prec, output_tend_s, output_tend_q, &
                    output_cape, output_dcape, output_mass_flux, output_prec_flux, &
                    zdu, mu, eu, du, md, ed, dp, dsubcld, &
                    zm_qc, rliq, output_rain_prod, output_dlf, &
                    aero, microp_st )

  !-----------------------------------------------------------------------------
  ! mesoscale coherent structure parameterization (MCSP)- modifies tendencies from zm_conv_main() prior to updating the state

  if (zm_param%mcsp_enabled) then

    ! initialize local output tendencies for MCSP
    call zm_tend_init( ncol, pver, local_tend_s, local_tend_q, local_tend_u, local_tend_v )

    do i = 1,ncol
      do k = 1,pver
        state_s(i,k) = state_t(i,k)*zm_const%cpair
      end do
    end do

    ! perform the MCSP calculations
    call zm_conv_mcsp_tend( ncol, ncol, pver, pverp, &
                            dtime, jctop, zm_const, zm_param, &
                            state_p_mid, state_p_int, state_p_del, &
                            state_s, state_qv, state_u, state_v, &
                            output_tend_s, output_tend_q, &
                            local_tend_s, local_tend_q, &
                            local_tend_u, local_tend_v, &
                            mcsp_ds_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out, &
                            mcsp_freq, mcsp_shear, zm_depth )

    ! add MCSP tendencies to ZM convective tendencies
    do i = 1,ncol
      do k = 1,pver
        output_tend_s(i,k) = output_tend_s(i,k) + local_tend_s(i,k)
        output_tend_q(i,k) = output_tend_q(i,k) + local_tend_q(i,k)
        output_tend_u(i,k) = output_tend_u(i,k) + local_tend_u(i,k)
        output_tend_v(i,k) = output_tend_v(i,k) + local_tend_v(i,k)
      end do
    end do

  end if

  !-----------------------------------------------------------------------------
  ! apply tendencies from zm_conv_main() & MCSP to local copy of state variables

  call zm_physics_update( ncol, dtime, &
                          state_phis, local_state_zm, local_state_zi, &
                          state_p_mid, state_p_int, state_p_del, &
                          local_state_t, local_state_qv, &
                          output_tend_s, output_tend_q)

  !-----------------------------------------------------------------------------
  ! Compute the precipitation, rain evaporation, and snow formation/melting
  ! Note - this routine expects an updated state following zm_conv_main() (+MCSP)

  ! initialize local output tendencies for zm_conv_evap()
  call zm_tend_init( ncol, pver, local_tend_s, local_tend_q, local_tend_u, local_tend_v )

  ! perform the convective evaporation calculations
  call zm_conv_evap(ncol, ncol, pver, pverp, dtime, &
                    state_p_mid, state_p_del, &
                    local_state_t, local_state_qv, &
                    output_rain_prod, state_cldfrac, &
                    local_tend_s, local_tend_q, &
                    tend_s_snwprd, tend_s_snwevmlt, &
                    output_prec, output_snow, ntprprd, ntsnprd, &
                    output_prec_flux, output_snow_flux, microp_st)

  ! add tendencies from zm_conv_evap() to output tendencies
  do i = 1,ncol
    do k = 1,pver
      output_tend_s(i,k) = output_tend_s(i,k) + local_tend_s(i,k)
      output_tend_q(i,k) = output_tend_q(i,k) + local_tend_q(i,k)
      evap_ds_out(i,k)   = local_tend_s(i,k)
      evap_dq_out(i,k)   = local_tend_q(i,k)
    end do
  end do

  !-----------------------------------------------------------------------------
  ! convective momentum transport

  ! initialize local output tendencies for zm_transport_momentum()
  call zm_tend_init( ncol, pver, local_tend_s, local_tend_q, tx_wind_tend(:,:,1), tx_wind_tend(:,:,2) )

  do i = 1,ncol
    do k = 1,pver
      ! zm_transport_momentum expects winds that may have been modified by MCSP,
      ! but normally this is disabled, so U/V tendencies will be zero
      tx_winds(i,k,1) = state_u(i,k) + output_tend_u(i,k)*dtime
      tx_winds(i,k,2) = state_v(i,k) + output_tend_v(i,k)*dtime
    end do
  end do

  call zm_transport_momentum( ncol, ncol, pver, pverp, tx_winds, 2, &
                              mu, md, du, eu, ed, dp, &
                              jt, maxg, ideep, 1, lengath, &
                              tx_wind_tend, tx_pguall, tx_pgdall, &
                              tx_icwu, tx_icwd, dtime, local_tend_s )

  ! add tendencies from zm_transport_momentum() to output tendencies
  do i = 1,ncol
    do k = 1,pver
      output_tend_s(i,k) = output_tend_s(i,k) + local_tend_s(i,k)
      output_tend_u(i,k) = output_tend_u(i,k) + tx_wind_tend(i,k,1)
      output_tend_v(i,k) = output_tend_v(i,k) + tx_wind_tend(i,k,2)
    end do
  end do

  !-----------------------------------------------------------------------------
  ! convective tracer transport

  ! this is just a placeholder for now
  ! call zm_transport_tracer(...)

  !-----------------------------------------------------------------------------
  ! populate deep convection activity flag

  if (lengath.gt.0) then
    do i=1,lengath
      output_activity(ideep(i)) = 1
    end do
  end if

  !-----------------------------------------------------------------------------
  ! convert dry static energy tendency to temperature tendency

  do i = 1,ncol
    do k = 1,pver
      output_tend_t(i,k) = output_tend_s(i,k)/zm_const%cpair
    end do
  end do

  !-----------------------------------------------------------------------------
  return
end subroutine zm_eamxx_bridge_run_c

!===================================================================================================

end module zm_eamxx_bridge_main
