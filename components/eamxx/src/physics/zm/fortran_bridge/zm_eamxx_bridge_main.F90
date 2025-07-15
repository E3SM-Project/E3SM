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

subroutine zm_eamxx_bridge_init_c( pcols_in, pver_in ) bind(C)
  use mpi
  use zm_conv_types,   only: zm_const_t, zm_param_t
  use zm_conv,         only: zm_const, zm_param
  use zm_conv_types,   only: zm_const_set_for_testing, zm_param_set_for_testing
  use zm_conv_types,   only: zm_param_mpi_broadcast, zm_param_print
  use zm_eamxx_bridge_wv_saturation, only: wv_sat_init
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int), value, intent(in) :: pcols_in
  integer(kind=c_int), value, intent(in) :: pver_in
  !-----------------------------------------------------------------------------
  ! Local variables
  integer :: mpi_rank, ierror
  !-----------------------------------------------------------------------------
  pcols = pcols_in
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

subroutine zm_eamxx_bridge_run_c( ncol, is_first_step, state_phis, state_t, state_q ) bind(C)
  use zm_aero_type,          only: zm_aero_t
  use zm_microphysics_state, only: zm_microp_st
  use zm_conv,               only: zm_convr
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int),  value,              intent(in) :: ncol
  logical(kind=c_bool), value,              intent(in) :: is_first_step
  real(kind=c_real), dimension(pcols)     , intent(in) :: state_phis   ! input state surface geopotential height
  real(kind=c_real), dimension(pcols,pver), intent(in) :: state_t      ! input state temperature
  real(kind=c_real), dimension(pcols,pver), intent(in) :: state_q      ! input state water vapor
  !-----------------------------------------------------------------------------
  ! Local variables
  integer :: i
  ! arguments for zm_convr - order consistent with current interface
  ! integer  :: lchnk = 0
  ! 
  ! logical  :: is_first_step
  ! real(r8), dimension(pcols,pver) :: state_t      ! input state temperature
  ! real(r8), dimension(pcols,pver) :: state_q      ! input state water vapor
  ! real(r8), dimension(pcols)      :: prec         ! output precipitation
  ! integer,  dimension(pcols)      :: jctop        ! output top-of-deep-convection indices
  ! integer,  dimension(pcols)      :: jcbot        ! output bot-of-deep-convection indices
  ! real(r8), dimension(pcols)      :: pblh         ! input planetary boundary layer height
  ! real(r8), dimension(pcols,pver) :: state_zm     ! input state altitude at mid-levels
  
  ! real(r8), dimension(pcols,pverp):: state_zi     ! input state altitude at interfaces
  ! real(r8), dimension(pcols,pver) :: ptend_loc_q  ! output tendency of water vapor
  ! real(r8), dimension(pcols,pver) :: ptend_loc_s  ! output tendency of dry statis energy
  ! real(r8), dimension(pcols,pver) :: state_pmid   ! input state pressure at mid-levels
  ! real(r8), dimension(pcols,pverp):: state_pint   ! input state pressure at interfaces
  ! real(r8), dimension(pcols,pver) :: state_pdel   ! input state pressure thickness
  ! real(r8), dimension(pcols,pver) :: state_omega  ! input state vertical pressure velocity
  ! real(r8)                        :: ztodt        ! model time increment x2
  ! real(r8), dimension(pcols,pverp):: mcon         ! convective mass flux--m sub c
  ! real(r8), dimension(pcols,pver) :: cme          ! condensation - evaporation
  ! real(r8), dimension(pcols)      :: cape         ! convective available potential energy
  ! real(r8), dimension(pcols)      :: tpert        ! thermal temperature excess
  ! real(r8), dimension(pcols,pver) :: dlf          ! detrained convective cloud water mixing ratio
  ! real(r8), dimension(pcols,pver) :: pflx         ! precip flux at each level
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
  ! write(iulog,*) 'zm_eamxx_bridge_run_c - 00 - is_first_step: ',is_first_step
  call shr_sys_flush(iulog)
  if (is_first_step) then
    write(iulog,*) 'zm_eamxx_bridge_run_c - 00 - pcols: ',pcols
    call shr_sys_flush(iulog)
    write(iulog,*) 'zm_eamxx_bridge_run_c - 00 - ncol: ',ncol
    call shr_sys_flush(iulog)
  end if
  !-----------------------------------------------------------------------------
  write(iulog,*) 'zm_eamxx_bridge_run_c - 01'
  call shr_sys_flush(iulog)
  do i = 1,ncol
    write(iulog,*) 'zm_eamxx_bridge_run_c - 01 - phis(',i,')   : ',state_phis(i)
    write(iulog,*) 'zm_eamxx_bridge_run_c - 01 - t(',i,',pver) : ',state_t(i,pver)
    ! write(iulog,*) 'zm_eamxx_bridge_run_c - 01 - t(',i,',1)    : ',state_t(i,1)
    write(iulog,*) 'zm_eamxx_bridge_run_c - 01 - q(',i,',pver) : ',state_q(i,pver)
    ! write(iulog,*) 'zm_eamxx_bridge_run_c - 01 - q(',i,',1)    : ',state_q(i,1)
    call shr_sys_flush(iulog)
  end do
  !-----------------------------------------------------------------------------
  ! ! Call the primary Zhang-McFarlane convection parameterization
  ! call zm_convr( lchnk, ncol, is_first_step, &
  !                state_t, state_q, &
  !                prec, &
  !                jctop, jcbot, &
  !                pblh, &
  !                state_zm, state_phis, state_zi, &
  !                ptend_loc_q, ptend_loc_s, &
  !                state_pmid, state_pint, state_pdel, state_omega, &
  !                0.5*ztodt, &
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
  return
end subroutine zm_eamxx_bridge_run_c

!===================================================================================================


end module zm_eamxx_bridge_main
