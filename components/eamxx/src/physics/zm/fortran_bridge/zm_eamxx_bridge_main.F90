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
  use zm_conv_types,   only: zm_const_t, zm_param_t
  use zm_conv,         only: zm_const, zm_param
  use zm_conv_types,   only: zm_const_set_for_testing, zm_param_set_for_testing
  use zm_conv_types,   only: zm_param_mpi_broadcast, zm_param_print
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
  zm_param%mcsp_enabled    = .false.
  zm_param%trig_dcape      = .false.
  zm_param%trig_ull        = .false.
  zm_param%clos_dyn_adj    = .false.
  !-----------------------------------------------------------------------------
  call wv_sat_init()
  !-----------------------------------------------------------------------------
  return
end subroutine zm_eamxx_bridge_init_c

!===================================================================================================

subroutine zm_eamxx_bridge_run_c( ncol, dtime, is_first_step, &
                                  state_phis, &
                                  state_p_mid, state_p_int, state_p_del, &
                                  state_t, state_qv, state_qc, &
                                  state_omega, state_pblh, &
                                  output_prec, output_tend_s, output_tend_q, &
                                  output_prec_flux, output_mass_flux ) bind(C)
  use zm_aero_type,          only: zm_aero_t
  use zm_microphysics_state, only: zm_microp_st
  use zm_conv,               only: zm_convr
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int),               value, intent(in ) :: ncol
  logical(kind=c_real),              value, intent(in ) :: dtime
  logical(kind=c_bool),              value, intent(in ) :: is_first_step
  real(kind=c_real), dimension(pcols),      intent(in ) :: state_phis        ! input state surface geopotential height
  real(kind=c_real), dimension(pcols,pver), intent(in ) :: state_p_mid       ! input state mid-point pressure
  real(kind=c_real), dimension(pcols,pver), intent(in ) :: state_p_int       ! input state interface pressure
  real(kind=c_real), dimension(pcols,pver), intent(in ) :: state_p_del       ! input state pressure thickness
  real(kind=c_real), dimension(pcols,pver), intent(in ) :: state_t           ! input state temperature
  real(kind=c_real), dimension(pcols,pver), intent(in ) :: state_qv          ! input state water vapor
  real(kind=c_real), dimension(pcols,pver), intent(in ) :: state_qc          ! input state cloud liquid water
  real(kind=c_real), dimension(pcols,pver), intent(in ) :: state_omega       ! input state vertical pressure velocity
  real(kind=c_real), dimension(pcols),      intent(in ) :: state_pblh        ! input planetary boundary layer height
  real(kind=c_real), dimension(pcols),      intent(out) :: output_prec       ! output total precipitation            (prec)
  real(kind=c_real), dimension(pcols,pver), intent(out) :: output_tend_s     ! output tendency of dry static energy  (ptend_loc_s)
  real(kind=c_real), dimension(pcols,pver), intent(out) :: output_tend_q     ! output tendency of water vapor        (ptend_loc_q)
  real(kind=c_real), dimension(pcols,pverp),intent(out) :: output_prec_flux  ! output precip flux at each mid-levels (pflx)
  real(kind=c_real), dimension(pcols,pverp),intent(out) :: output_mass_flux  ! output convective mass flux--m sub c  (mcon)

  !-----------------------------------------------------------------------------
  ! Local variables
  integer :: i,j,k

  ! arguments for zm_convr - order consistent with current interface
  integer :: lchnk = 0

  integer,  dimension(pcols)      :: jctop        ! output top-of-deep-convection indices
  integer,  dimension(pcols)      :: jcbot        ! output bot-of-deep-convection indices
  ! real(r8), dimension(pcols,pver) :: state_zm     ! input state altitude at mid-levels
  ! real(r8), dimension(pcols,pverp):: state_zi     ! input state altitude at interfaces
  real(r8), dimension(pcols,pverp):: mcon         ! convective mass flux--m sub c
  real(r8), dimension(pcols,pver) :: cme          ! condensation - evaporation
  real(r8), dimension(pcols)      :: cape         ! convective available potential energy
  ! real(r8), dimension(pcols)      :: tpert        ! thermal temperature excess
  ! real(r8), dimension(pcols,pver) :: dlf          ! detrained convective cloud water mixing ratio
  ! real(r8), dimension(pcols,pverp):: pflx         ! precip flux at each level
  ! real(r8), dimension(pcols,pver) :: zdu          ! detraining mass flux
  ! real(r8), dimension(pcols,pver) :: rprd         ! rain production rate
  ! real(r8), dimension(pcols,pver) :: mu           ! upward cloud mass flux
  ! real(r8), dimension(pcols,pver) :: md           ! entrainment in updraft
  ! real(r8), dimension(pcols,pver) :: du           ! detrainment in updraft
  ! real(r8), dimension(pcols,pver) :: eu           ! downward cloud mass flux
  ! real(r8), dimension(pcols,pver) :: ed           ! entrainment in downdraft
  ! real(r8), dimension(pcols,pver) :: dp           ! layer thickness [mb]
  ! real(r8), dimension(pcols)      :: dsubcld      ! sub-cloud layer thickness
  ! integer,  dimension(pcols)      :: jt           ! top level index of convection
  ! integer,  dimension(pcols)      :: maxg         ! gathered values of maxi
  ! integer,  dimension(pcols)      :: ideep        ! flag to indicate ZM is active
  ! integer                         :: lengath      ! number of gathered columns per chunk
  ! real(r8), dimension(pcols,pver) :: ql           ! grid slice of cloud liquid water
  ! real(r8), dimension(pcols)      :: rliq         ! reserved liquid (not yet in cldliq) for energy integrals
  ! real(r8), dimension(pcols)      :: landfrac     ! land fraction
  ! real(r8), dimension(pcols,pver), target :: t_star       ! DCAPE T from time step n-1
  ! real(r8), dimension(pcols,pver), target :: q_star       ! DCAPE q from time step n-1
  ! real(r8), dimension(pcols)      :: dcape        ! DCAPE cape change
  ! type(zm_aero_t), allocatable    :: aero(:)      ! derived type for aerosol information
  ! real(r8), dimension(pcols,pver) :: qi           ! grid slice of cloud ice
  ! real(r8), dimension(pcols,pver) :: dif          ! detrained convective cloud ice mixing ratio
  ! real(r8), dimension(pcols,pver) :: dnlf         ! detrained convective cloud water num concen
  ! real(r8), dimension(pcols,pver) :: dnif         ! detrained convective cloud ice num concen
  ! real(r8), dimension(pcols,pver) :: dsf          ! detrained convective snow mixing ratio
  ! real(r8), dimension(pcols,pver) :: dnsf         ! detrained convective snow num concen
  ! real(r8), dimension(pcols,pver) :: sprd         ! ???
  ! real(r8), dimension(pcols)      :: rice         ! reserved ice (not yet in cldice) for energy integrals
  ! real(r8), dimension(pcols,pver) :: frz          ! ???
  ! real(r8), dimension(pcols,pver) :: mudpcu       ! width parameter of droplet size distr
  ! real(r8), dimension(pcols,pver) :: lambdadpcu   ! slope of cloud liquid size distr
  ! type(zm_microp_st)              :: microp_st    ! ZM microphysics data structure
  ! real(r8), dimension(pcols,pver) :: wuc          ! pbuf variable for in-cloud vertical velocity
  !-----------------------------------------------------------------------------
  ! if (is_first_step) then
  !   write(iulog,*) 'zm_eamxx_bridge_run_c - pcols: ',pcols
  !   write(iulog,*) 'zm_eamxx_bridge_run_c - ncol: ',ncol
  !   call shr_sys_flush(iulog)
  ! end if
  !-----------------------------------------------------------------------------
  ! assign fake values for checking data made it back to C++
  do i = 1,ncol
    output_prec(i) = 1
    do k = 1,pver
      output_tend_s(i,k) = 10000 + i*100.0 + k*1.0
      output_tend_q(i,k) = 20000 + i*100.0 + k*1.0
      ! output_tend_s(i,k) = 2.0
      ! output_tend_q(i,k) = 3.0
    end do
  end do
  !-----------------------------------------------------------------------------
  ! do i = 1,ncol
  !   ! write(iulog,*) 'zm_eamxx_bridge_run_c - prec(',i,') : ',output_prec(i)
  !   ! write(iulog,*) 'zm_eamxx_bridge_run_c - phis(',i,') : ',state_phis(i)
  !   do k = 1,6
  !     write(iulog,*) 'zm_eamxx_bridge_run_c - (',i,',',k,') pmid / tend_s / tend_q : ',state_p_mid(i,k),' / ',output_tend_s(i,k),' / ',output_tend_q(i,k)
  !     ! write(iulog,*) 'zm_eamxx_bridge_run_c - tend_s(',i,',',k,') : ',output_tend_s(i,k)
  !     ! write(iulog,*) 'zm_eamxx_bridge_run_c - pmid(',i,',',k,') : ',state_p_mid(i,k)
  !   end do
  ! end do
  !-----------------------------------------------------------------------------
  ! ! Call the primary Zhang-McFarlane convection parameterization
  ! call zm_convr( lchnk, ncol, is_first_step, &
  !                state_t, state_q, &
  !                prec, &
  !                jctop, jcbot, &
  !                pblh, &
  !                state_zm, state_phis, state_zi, &
  !                output_tend_q, output_tend_s, &
  !                state_p_mid, state_p_int, state_p_del, state_omega, &
  !                0.5*dtime, &
  !                mcon, &
  !                cme, &
  !                cape, &
  !                tpert, &
  !                dlf, &
  !                pflx, &
  !                zdu, &
  !                rprd, &
  !                mu, md, du, eu, ed, dp, &
  !                dsubcld, &
  !                jt, &
  !                maxg, ideep, lengath, &
  !                ql, &
  !                rliq, &
  !                landfrac, &
  !                t_star, q_star, dcape, &
  !                aero(lchnk), &
  !                qi, dif, dnlf, dnif, dsf, dnsf, sprd, rice, frz, &
  !                mudpcu, &
  !                lambdadpcu, &
  !                microp_st, &
  !                wuc )
  !-----------------------------------------------------------------------------
  ! call zm_conv_evap(state1%ncol, state1%lchnk, &
  !                   state1%t, state1%pmid, state1%pdel, &
  !                   state1%q(1:pcols,1:pver,1), &
  !                   ptend_loc%s, &
  !                   tend_s_snwprd, tend_s_snwevmlt, &
  !                   ptend_loc%q(:pcols,:pver,1), &
  !                   rprd, cld, ztodt, prec, snow, &
  !                   ntprprd, ntsnprd, &
  !                   flxprec, flxsnow, &
  !                   sprd, old_snow)
  !-----------------------------------------------------------------------------
  return
end subroutine zm_eamxx_bridge_run_c

!===================================================================================================


end module zm_eamxx_bridge_main
