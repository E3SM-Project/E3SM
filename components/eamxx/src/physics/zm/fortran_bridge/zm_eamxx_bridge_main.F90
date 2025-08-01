module zm_eamxx_bridge_main
  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  use iso_c_binding
  use cam_logfile,   only: iulog
  use shr_sys_mod,   only: shr_sys_flush
  use zm_eamxx_bridge_params, only: masterproc, r8, pcols, pver, pverp, top_lev
  !-----------------------------------------------------------------------------
  implicit none

  private
  !-----------------------------------------------------------------------------
  ! public methods for bridging
  public :: zm_eamxx_bridge_init_c
  public :: zm_eamxx_bridge_run_c
  !-----------------------------------------------------------------------------
  ! public variables

  !-----------------------------------------------------------------------------
  ! private variables

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

subroutine zm_eamxx_bridge_init_c( pcol_in, pver_in ) bind(C)
  use mpi
  use zm_conv,       only: zm_const, zm_param
  use zm_conv_types, only: zm_const_set_for_testing, zm_param_set_for_testing
  use zm_conv_types, only: zm_param_mpi_broadcast, zm_param_print
  use zm_eamxx_bridge_wv_saturation, only: wv_sat_init
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int), value, intent(in) :: pcol_in
  integer(kind=c_int), value, intent(in) :: pver_in
  !-----------------------------------------------------------------------------
  ! Local variables
  integer :: mpi_rank, ierror
  !-----------------------------------------------------------------------------
  pcols = pcol_in
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
  call zm_param_mpi_broadcast(zm_param)
  if (masterproc) call zm_param_print(zm_param)
  !-----------------------------------------------------------------------------
  ! make sure we are turning off the extra stuff
  zm_param%zm_microp       = .false.
  zm_param%trig_dcape      = .false.
  zm_param%trig_ull        = .false.
  zm_param%clos_dyn_adj    = .false.
  zm_param%mcsp_enabled    = .true.
  !-----------------------------------------------------------------------------
  call wv_sat_init()
  !-----------------------------------------------------------------------------
  return
end subroutine zm_eamxx_bridge_init_c

!===================================================================================================

subroutine zm_eamxx_bridge_run_c( ncol, dtime, is_first_step, &
                                  state_phis, state_zm, state_zi, &
                                  state_p_mid, state_p_int, state_p_del, &
                                  state_t, state_qv, state_qc, state_u, state_v, &
                                  state_omega, state_cldfrac, state_pblh, landfrac, &
                                  output_prec, output_cape, &
                                  output_tend_s, output_tend_q, &
                                  output_tend_u, output_tend_v, &
                                  output_prec_flux, output_mass_flux ) bind(C)
  use zm_conv,                  only: zm_const, zm_param
  use zm_aero_type,             only: zm_aero_t
  use zm_microphysics_state,    only: zm_microp_st
  use zm_eamxx_bridge_methods,  only: zm_tend_init, zm_physics_update
  use zm_conv,                  only: zm_convr, zm_conv_evap
  use zm_conv_mcsp,             only: zm_conv_mcsp_tend
  use zm_transport,             only: zm_transport_momentum
  ! use zm_transport,             only: zm_transport_tracer
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int),               value, intent(in   ) :: ncol
  real(kind=c_real),                 value, intent(in   ) :: dtime
  logical(kind=c_bool),              value, intent(in   ) :: is_first_step
  real(kind=c_real), dimension(pcols),      intent(in   ) :: state_phis         ! input state surface geopotential height
  real(kind=c_real), dimension(pcols,pver), intent(in   ) :: state_zm           ! input state altitude at mid-levels
  real(kind=c_real), dimension(pcols,pverp),intent(in   ) :: state_zi           ! input state altitude at interfaces
  real(kind=c_real), dimension(pcols,pver), intent(in   ) :: state_p_mid        ! input state mid-point pressure
  real(kind=c_real), dimension(pcols,pverp),intent(in   ) :: state_p_int        ! input state interface pressure
  real(kind=c_real), dimension(pcols,pver), intent(in   ) :: state_p_del        ! input state pressure thickness
  real(kind=c_real), dimension(pcols,pver), intent(in   ) :: state_t            ! input state temperature
  real(kind=c_real), dimension(pcols,pver), intent(in   ) :: state_qv           ! input state water vapor
  real(kind=c_real), dimension(pcols,pver), intent(in   ) :: state_qc           ! input state cloud liquid water          (ql)
  real(kind=c_real), dimension(pcols,pver), intent(in   ) :: state_u            ! input state zonal wind
  real(kind=c_real), dimension(pcols,pver), intent(in   ) :: state_v            ! input state meridional wind
  real(kind=c_real), dimension(pcols,pver), intent(in   ) :: state_omega        ! input state vertical pressure velocity
  real(kind=c_real), dimension(pcols,pver), intent(in   ) :: state_cldfrac      ! input state cloud fraction              (cld)
  real(kind=c_real), dimension(pcols),      intent(in   ) :: state_pblh         ! input planetary boundary layer height   (pblh)
  real(kind=c_real), dimension(pcols),      intent(in   ) :: landfrac           ! land fraction
  real(kind=c_real), dimension(pcols),      intent(  out) :: output_prec        ! output total precipitation              (prec)
  real(kind=c_real), dimension(pcols),      intent(  out) :: output_cape        ! output convective avail. pot. energy    (cape)
  real(kind=c_real), dimension(pcols,pver), intent(  out) :: output_tend_s      ! output tendency of dry static energy    (ptend_loc_s)
  real(kind=c_real), dimension(pcols,pver), intent(  out) :: output_tend_q      ! output tendency of water vapor          (ptend_loc_q)
  real(kind=c_real), dimension(pcols,pver), intent(  out) :: output_tend_u      ! output tendency of zonal wind           (ptend_loc_u)
  real(kind=c_real), dimension(pcols,pver), intent(  out) :: output_tend_v      ! output tendency of meridional wind      (ptend_loc_v)
  real(kind=c_real), dimension(pcols,pverp),intent(  out) :: output_prec_flux   ! output precip flux at each mid-levels   (pflx)
  real(kind=c_real), dimension(pcols,pverp),intent(  out) :: output_mass_flux   ! output convective mass flux--m sub c    (mcon)
  !-----------------------------------------------------------------------------
  ! Local variables
  integer :: i,k

  integer :: diag_kt
  integer :: diag_kb

  logical :: loc_is_first_step

  ! arguments for zm_convr - order consistent with current interface
  integer :: lchnk = 0

  integer,  dimension(pcols)      :: jctop          ! output top-of-deep-convection indices
  integer,  dimension(pcols)      :: jcbot          ! output bot-of-deep-convection indices
  ! real(r8), dimension(pcols,pverp):: mcon           ! convective mass flux--m sub c
  real(r8), dimension(pcols,pver) :: cme            ! condensation - evaporation
  ! real(r8), dimension(pcols)      :: cape           ! convective available potential energy
  real(r8), dimension(pcols)      :: tpert          ! thermal temperature excess
  real(r8), dimension(pcols,pver) :: dlf            ! detrained convective cloud water mixing ratio
  real(r8), dimension(pcols,pverp):: pflx           ! precip flux at each level
  real(r8), dimension(pcols,pver) :: zdu            ! detraining mass flux
  real(r8), dimension(pcols,pver) :: rprd           ! rain production rate
  real(r8), dimension(pcols,pver) :: mu             ! upward cloud mass flux
  real(r8), dimension(pcols,pver) :: md             ! entrainment in updraft
  real(r8), dimension(pcols,pver) :: du             ! detrainment in updraft
  real(r8), dimension(pcols,pver) :: eu             ! downward cloud mass flux
  real(r8), dimension(pcols,pver) :: ed             ! entrainment in downdraft
  real(r8), dimension(pcols,pver) :: dp             ! layer thickness [mb]
  real(r8), dimension(pcols)      :: dsubcld        ! sub-cloud layer thickness
  integer,  dimension(pcols)      :: jt             ! top level index of convection
  integer,  dimension(pcols)      :: maxg           ! gathered values of maxi
  integer,  dimension(pcols)      :: ideep          ! flag to indicate ZM is active
  integer                         :: lengath        ! number of gathered columns per chunk
  real(r8), dimension(pcols)      :: rliq           ! reserved liquid (not yet in cldliq) for energy integrals
  real(r8), dimension(pcols,pver), target :: t_star ! DCAPE T from time step n-1
  real(r8), dimension(pcols,pver), target :: q_star ! DCAPE q from time step n-1
  real(r8), dimension(pcols)      :: dcape          ! DCAPE cape change
  real(r8), dimension(pcols,pver) :: qi             ! grid slice of cloud ice
  real(r8), dimension(pcols,pver) :: dif            ! detrained convective cloud ice mixing ratio
  real(r8), dimension(pcols,pver) :: dnlf           ! detrained convective cloud water num concen
  real(r8), dimension(pcols,pver) :: dnif           ! detrained convective cloud ice num concen
  real(r8), dimension(pcols,pver) :: dsf            ! detrained convective snow mixing ratio
  real(r8), dimension(pcols,pver) :: dnsf           ! detrained convective snow num concen
  real(r8), dimension(pcols,pver) :: sprd           ! snow production rate
  real(r8), dimension(pcols)      :: rice           ! reserved ice (not yet in cldice) for energy integrals
  real(r8), dimension(pcols,pver) :: frz            ! freezing rate
  real(r8), dimension(pcols,pver) :: mudpcu         ! width parameter of droplet size distr
  real(r8), dimension(pcols,pver) :: lambdadpcu     ! slope of cloud liquid size distr
  type(zm_aero_t)                 :: aero           ! derived type for aerosol information
  type(zm_microp_st)              :: microp_st      ! ZM microphysics data structure
  real(r8), dimension(pcols,pver) :: wuc            ! pbuf variable for in-cloud vertical velocity

  real(r8), dimension(pcols,pver) :: state_s

  ! local copy of state variables for calling zm_conv_evap()
  real(r8), dimension(pcols,pver) :: local_state_t
  real(r8), dimension(pcols,pver) :: local_state_qv
  real(r8), dimension(pcols,pver) :: local_state_zm
  real(r8), dimension(pcols,pverp):: local_state_zi

  ! temporary local tendencies for calling zm_conv_evap()
  real(r8), dimension(pcols,pver) :: local_tend_s        ! output tendency of dry static energy   (ptend_loc_s)
  real(r8), dimension(pcols,pver) :: local_tend_q        ! output tendency of water vapor         (ptend_loc_q)
  real(r8), dimension(pcols,pver) :: local_tend_u        ! output tendency of zonal wind
  real(r8), dimension(pcols,pver) :: local_tend_v        ! output tendency of meridional wind

  real(r8), dimension(pcols,pver) :: tend_s_snwprd
  real(r8), dimension(pcols,pver) :: tend_s_snwevmlt
  real(r8), dimension(pcols,pver) :: snow
  real(r8), dimension(pcols,pver) :: ntprprd
  real(r8), dimension(pcols,pver) :: ntsnprd
  ! real(r8), dimension(pcols,pverp):: flxprec
  real(r8), dimension(pcols,pverp):: flxsnow

  ! used in momentum transport calculations
   real(r8), dimension(pcols,pver,2) :: tx_winds
   real(r8), dimension(pcols,pver,2) :: tx_wind_tend
   real(r8), dimension(pcols,pver,2) :: tx_pguall
   real(r8), dimension(pcols,pver,2) :: tx_pgdall
   real(r8), dimension(pcols,pver,2) :: tx_icwu
   real(r8), dimension(pcols,pver,2) :: tx_icwd
  
  logical :: old_snow  = .true. ! flag to use snow production from zm_conv_evap - set false to use zm microphysics
  !-----------------------------------------------------------------------------
  loc_is_first_step = is_first_step
  diag_kt = 1
  diag_kb = pver
  !-----------------------------------------------------------------------------
  ! if (is_first_step) then
  !   write(iulog,*) 'zm_eamxx_bridge_run_c - pcols: ',pcols
  !   write(iulog,*) 'zm_eamxx_bridge_run_c - ncol: ',ncol
  !   call shr_sys_flush(iulog)
  ! end if
  !-----------------------------------------------------------------------------
  ! ! assign fake values for checking data made it back to C++
  ! do i = 1,ncol
  !   output_prec(i) = 1
  !   output_cape(i) = 2
  !   do k = 1,pver
  !     output_tend_s(i,k) = 10000 + i*100.0 + k*1.0
  !     output_tend_q(i,k) = 20000 + i*100.0 + k*1.0
  !   end do
  ! end do
  !-----------------------------------------------------------------------------
  ! initialize output tendencies - normally done by physics_ptend_init()
  do i = 1,ncol
    output_prec(i) = 0
    do k = 1,pver
      output_tend_s(i,k) = 0
      output_tend_q(i,k) = 0
      output_tend_u(i,k) = 0
      output_tend_v(i,k) = 0
    end do
  end do
  !-----------------------------------------------------------------------------
  ! update local copy of state variables for zm_conv_evap()
  do i = 1,ncol
    do k = 1,pver
      local_state_t (i,k) = state_t (i,k)
      local_state_qv(i,k) = state_qv(i,k)
      local_state_zm(i,k) = state_zm(i,k)
      local_state_zi(i,k) = state_zi(i,k)
    end do
  end do
  !-----------------------------------------------------------------------------
  ! do i = 1,ncol
  !   write(iulog,*) 'zm_eamxx_bridge_run_c 0 - prec(',i,') : ',output_prec(i)
  !   do k = diag_kt,diag_kb
  !     write(iulog,*) 'zm_eamxx_bridge_run_c 0 - (',i,',',k,') pmid / tend_s / tend_q : ',state_p_mid(i,k),' / ',output_tend_s(i,k),' / ',output_tend_q(i,k)
  !   end do
  ! end do
  !-----------------------------------------------------------------------------
  ! Call the primary Zhang-McFarlane convection parameterization
  call zm_convr( lchnk, ncol, loc_is_first_step, &
                 state_t, state_qv, &
                 output_prec, &
                 jctop, jcbot, &
                 state_pblh, &
                 state_zm, state_phis, state_zi, &
                 output_tend_q, output_tend_s, &
                 state_p_mid, state_p_int, state_p_del, state_omega, &
                 0.5*dtime, &
                 output_mass_flux, &
                 cme, &
                 output_cape, &
                 tpert, &
                 dlf, &
                 pflx, &
                 zdu, &
                 rprd, &
                 mu, md, du, eu, ed, dp, &
                 dsubcld, &
                 jt, &
                 maxg, ideep, lengath, &
                 state_qc, &
                 rliq, &
                 landfrac, &
                 t_star, q_star, dcape, &
                 aero, &
                 qi, dif, dnlf, dnif, dsf, dnsf, sprd, rice, frz, &
                 mudpcu, lambdadpcu, &
                 microp_st, &
                 wuc )
  !-----------------------------------------------------------------------------
  ! do i = 1,ncol
  !   write(iulog,*) 'zm_eamxx_bridge_run_c 1 - (',i,') prec / cape : ',output_prec(i),' / ',output_cape(i)
  !   ! do k = diag_kt,diag_kb
  !   !   write(iulog,*) 'zm_eamxx_bridge_run_c 1 - (',i,',',k,') pmid / tend_s / tend_q : ',state_p_mid(i,k),' / ',output_tend_s(i,k),' / ',output_tend_q(i,k)
  !   ! end do
  ! end do

  !----------------------------------------------------------------------------
  ! mesoscale coherent structure parameterization (MCSP)
  ! This modifies the tendencies from zm_convr() prior to updating the state
  if (zm_param%mcsp_enabled) then

    ! initialize local output tendencies for MCSP
    call zm_tend_init( ncol, pver, local_tend_s, local_tend_q, local_tend_u, local_tend_v )

    do i = 1,ncol
      do k = 1,pver
        state_s(i,k) = state_t(i,k)*zm_const%cpair
      end do
    end do

    ! perform the MCSP calculations
    call zm_conv_mcsp_tend( lchnk, pcols, ncol, pver, pverp, &
                            dtime, jctop, zm_const, zm_param, &
                            state_p_mid, state_p_int, state_p_del, &
                            state_s, state_qv, state_u, state_v, &
                            output_tend_s, output_tend_q, &
                            local_tend_s, local_tend_q, &
                            local_tend_u, local_tend_v )

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

  !----------------------------------------------------------------------------
  ! apply tendencies from zm_convr() & MCSP to local copy of state variables
  call zm_physics_update( ncol, dtime, state_phis, local_state_zm, local_state_zi, &
                          state_p_mid, state_p_int, state_p_del, &
                          local_state_t, local_state_qv, &
                          output_tend_s, output_tend_q)

  !-----------------------------------------------------------------------------
  ! Compute the precipitation, rain evaporation, and snow formation/melting
  ! Note - this routine expects an updated state following zm_convr() (+MCSP)

  ! initialize local output tendencies for zm_conv_evap()
  call zm_tend_init( ncol, pver, local_tend_s, local_tend_q, local_tend_u, local_tend_v )

  ! perform the convective evaporation calculations
  call zm_conv_evap(ncol, lchnk, &
                    local_state_t, state_p_mid, state_p_del, local_state_qv, &
                    local_tend_s, &
                    tend_s_snwprd, &
                    tend_s_snwevmlt, &
                    local_tend_q, &
                    rprd, state_cldfrac, &
                    dtime, &
                    output_prec, &
                    snow, &
                    ntprprd, &
                    ntsnprd, &
                    output_prec_flux, &
                    flxsnow, &
                    sprd, old_snow)

  ! add tendencies from zm_conv_evap() to output tendencies
  do i = 1,ncol
    do k = 1,pver
      output_tend_s(i,k) = output_tend_s(i,k) + local_tend_s(i,k)
      output_tend_q(i,k) = output_tend_q(i,k) + local_tend_q(i,k)
    end do
  end do

  ! apply tendencies from zm_conv_evap() to local copy of state variables
  call zm_physics_update( ncol, dtime, state_phis, local_state_zm, local_state_zi, &
                          state_p_mid, state_p_int, state_p_del, &
                          local_state_t, local_state_qv, &
                          local_tend_s, local_tend_q)

  !----------------------------------------------------------------------------
  ! convective momentum transport

  ! initialize local output tendencies for zm_conv_evap()
  call zm_tend_init( ncol, pver, local_tend_s, local_tend_q, tx_wind_tend(:,:,1), tx_wind_tend(:,:,2) )

  do i = 1,ncol
    do k = 1,pver
      tx_winds(i,k,1) = state_u(i,k)
      tx_winds(i,k,2) = state_v(i,k)
    end do
  end do

  call zm_transport_momentum( ncol, tx_winds, 2, &
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

  !----------------------------------------------------------------------------
  ! convective tracer transport

  !-----------------------------------------------------------------------------
  ! do i = 1,ncol
  !   write(iulog,*) 'zm_eamxx_bridge_run_c 2 - (',i,') prec / cape : ',output_prec(i),' / ',output_cape(i)
  !   ! do k = diag_kt,diag_kb
  !   !   write(iulog,*) 'zm_eamxx_bridge_run_c 2 - (',i,',',k,') pmid / tend_s / tend_q : ',state_p_mid(i,k),' / ',output_tend_s(i,k),' / ',output_tend_q(i,k)
  !   !   ! write(iulog,*) 'zm_eamxx_bridge_run_c 2 - (',i,',',k,') output_prec_flux : ',output_prec_flux(i,k)
  !   !   ! write(iulog,*) 'zm_eamxx_bridge_run_c 2 - (',i,',',k,') flxprec          : ',flxprec(i,k)
  !   ! end do
  ! end do
  !-----------------------------------------------------------------------------
  return
end subroutine zm_eamxx_bridge_run_c

!===================================================================================================

end module zm_eamxx_bridge_main
