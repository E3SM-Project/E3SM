module check_energy

!---------------------------------------------------------------------------------
! Purpose:
!
! Module to check 
!   1. vertically integrated total energy and water conservation for each
!      column within the physical parameterizations
!
!   2. global mean total energy conservation between the physics output state
!      and the input state on the next time step.
!
!   3. add a globally uniform heating term to account for any change of total energy in 2.
!
! Author: Byron Boville  Oct 31, 2002
!         
! Modifications:
!   2003-03  Boville  Add global energy check and fixer.        
!   2016-08  Kai Zhang (kai.zhang@pnnl.gov) & Phil Rasch 
!               1. Modifications to check_energy_chng (water conservation check part) for
!                  MG2 (now includes prognostic rain and snow)    
!               2. Better printout information for energy and water conservation check 
!               3. Additional water conservation check utilities 
!---------------------------------------------------------------------------------
!
!   2020-01  O. Guba Correct energy density function
!
!---------------------------------------------------------------------------------
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver, begchunk, endchunk
  use spmd_utils,      only: masterproc
  
  use phys_gmean,      only: gmean
  use physconst,       only: gravit, latvap, latice, cpair, cpairv
  use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_ptend_init
  use constituents,    only: cnst_get_ind, pcnst, cnst_name, cnst_get_type_byind, &
                             icldliq, icldice, irain, isnow
  use time_manager,    only: is_first_step
  use cam_logfile,     only: iulog
  use cam_abortutils,  only: endrun 
  use phys_control,    only: ieflx_opt

  implicit none
  private

! Public types:
  public check_tracers_data

! Public methods
  public :: check_energy_defaultopts ! set default namelist values
  public :: check_energy_setopts   ! set namelist values
  public :: check_energy_register  ! register fields in physics buffer
  public :: check_energy_get_integrals ! get energy integrals computed in check_energy_gmean
  public :: check_energy_init      ! initialization of module
  public :: check_energy_timestep_init  ! timestep initialization of energy integrals and cumulative boundary fluxes
  public :: check_energy_chng      ! check changes in integrals against cumulative boundary fluxes
  public :: check_energy_gmean     ! global means of physics input and output total energy
  public :: check_energy_fix       ! add global mean energy difference as a heating
  public :: check_tracers_init      ! initialize tracer integrals and cumulative boundary fluxes
  public :: check_tracers_chng      ! check changes in integrals against cumulative boundary fluxes
  public :: check_tracers_fini      ! free memory associated with check_tracers_data type variable

  public :: qflx_gmean              ! calculate global mean of qflx for water conservation check 
  public :: check_qflx              ! output qflx at certain locations for water conservation check  
  public :: check_prect             ! output prect at certain locations for water conservation check  
  public :: check_water             ! output water path at certain locations for water conservation check  

  public :: ieflx_gmean             ! calculate global mean of ieflx 
  public :: check_ieflx_fix         ! add ieflx to sensible heat flux 

  public :: energy_helper_eam_def

! Private module data

  logical  :: print_energy_errors = .false.
  character(len=16) :: microp_scheme 

  real(r8) :: teout_glob           ! global mean energy of output state
  real(r8) :: teinp_glob           ! global mean energy of input state
  real(r8) :: tedif_glob           ! global mean energy difference
  real(r8) :: psurf_glob           ! global mean surface pressure
  real(r8) :: ptopb_glob           ! global mean top boundary pressure
  real(r8) :: heat_glob            ! global mean heating rate
  real(r8) :: ieflx_glob           ! global mean implied internal energy flux 

! Physics buffer indices
  
  integer  :: teout_idx  = 0       ! teout index in physics buffer 
  integer  :: dtcore_idx = 0       ! dtcore index in physics buffer 

  type check_tracers_data
     real(r8), allocatable :: tracer(:,:) ! initial vertically integrated total (kinetic + static) energy
     real(r8), allocatable :: tracer_tnd(:,:) ! cumulative boundary flux of total energy
     integer :: count(pcnst)               ! count of values with significant imbalances
  end type check_tracers_data


!===============================================================================
contains
!===============================================================================

subroutine check_energy_defaultopts( &
   print_energy_errors_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   logical,          intent(out), optional :: print_energy_errors_out
!-----------------------------------------------------------------------

   if ( present(print_energy_errors_out) ) then
      print_energy_errors_out = print_energy_errors
   endif

end subroutine check_energy_defaultopts

!================================================================================================

subroutine check_energy_setopts( &
   print_energy_errors_in)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   logical,          intent(in), optional :: print_energy_errors_in
!-----------------------------------------------------------------------

   if ( present(print_energy_errors_in) ) then
      print_energy_errors = print_energy_errors_in
   endif

end subroutine check_energy_setopts

!================================================================================================

  subroutine check_energy_register()
!
! Register fields in the physics buffer.
! 
!-----------------------------------------------------------------------
    
    use physics_buffer, only : pbuf_add_field, dtype_r8, dyn_time_lvls
    use physics_buffer, only : pbuf_register_subcol
    use subcol_utils,   only : is_subcol_on

!-----------------------------------------------------------------------

! Request physics buffer space for fields that persist across timesteps.

    call pbuf_add_field('TEOUT', 'global',dtype_r8 , (/pcols,dyn_time_lvls/),      teout_idx)
    call pbuf_add_field('DTCORE','global',dtype_r8,  (/pcols,pver,dyn_time_lvls/),dtcore_idx)

    if(is_subcol_on()) then
      call pbuf_register_subcol('TEOUT', 'phys_register', teout_idx)
      call pbuf_register_subcol('DTCORE', 'phys_register', dtcore_idx)
    end if

  end subroutine check_energy_register

!===============================================================================

subroutine check_energy_get_integrals( tedif_glob_out, heat_glob_out )

!----------------------------------------------------------------------- 
! Purpose: Return energy integrals
!-----------------------------------------------------------------------

     real(r8), intent(out), optional :: tedif_glob_out
     real(r8), intent(out), optional :: heat_glob_out

!-----------------------------------------------------------------------

   if ( present(tedif_glob_out) ) then
      tedif_glob_out = tedif_glob
   endif
   if ( present(heat_glob_out) ) then
      heat_glob_out = heat_glob
   endif

end subroutine check_energy_get_integrals

!================================================================================================

  subroutine check_energy_init()
!
! Initialize the energy conservation module
! 
!-----------------------------------------------------------------------
    use cam_history,       only: addfld, horiz_only, add_default
    use phys_control,      only: phys_getopts

    implicit none

    logical          :: history_budget
    integer          :: history_budget_histfile_num ! output history file number for budget fields

!-----------------------------------------------------------------------

    call phys_getopts( history_budget_out = history_budget, &
                       microp_scheme_out  = microp_scheme,   &
                       history_budget_histfile_num_out = history_budget_histfile_num)

! register history variables
    call addfld('TEINP', horiz_only,    'A', 'W/m2', 'Total energy of physics input')
    call addfld('TEOUT', horiz_only,    'A', 'W/m2', 'Total energy of physics output')
    call addfld('TEFIX', horiz_only,    'A', 'W/m2', 'Total energy after fixer')
    call addfld('EFIX',   horiz_only,  'A', 'W/m2', 'Effective sensible heat flux due to energy fixer')
    call addfld('DTCORE', (/ 'lev' /), 'A', 'K/s' , 'T tendency due to dynamical core')

    call addfld('BC01Q', horiz_only,   'I', 'kg/m2', 'q after process')
    call addfld('BC01QL', horiz_only,  'I', 'kg/m2', 'ql after process')
    call addfld('BC01QI', horiz_only,  'I', 'kg/m2', 'qi after process')
    call addfld('BC01QR', horiz_only,  'I', 'kg/m2', 'qr after process')
    call addfld('BC01QS', horiz_only,  'I', 'kg/m2', 'qs after process')
    call addfld('BC01TW', horiz_only,  'I', 'kg/m2', 'total water after process')

    call addfld('BC02Q', horiz_only,   'I', 'kg/m2', 'q after process')
    call addfld('BC02QL', horiz_only,  'I', 'kg/m2', 'ql after process')
    call addfld('BC02QI', horiz_only,  'I', 'kg/m2', 'qi after process')
    call addfld('BC02QR', horiz_only,  'I', 'kg/m2', 'qr after process')
    call addfld('BC02QS', horiz_only,  'I', 'kg/m2', 'qs after process')
    call addfld('BC02TW', horiz_only,  'I', 'kg/m2', 'total water after process')

    call addfld('AC01QFLX', horiz_only,    'A', 'kg/m2/s', 'total water change due to water flux ')
    call addfld('AC02QFLX', horiz_only,    'A', 'kg/m2/s', 'total water change due to water flux ')
    call addfld('BC01QFLX', horiz_only,    'A', 'kg/m2/s', 'total water change due to water flux ')


    if (masterproc) then
       write (iulog,*) ' print_energy_errors is set', print_energy_errors
    endif

    if ( history_budget ) then
       call add_default ('DTCORE   '  , history_budget_histfile_num, ' ')
    end if

    if(ieflx_opt>0) then 
       call addfld('IEFLX',    horiz_only, 'A', 'W/m2', 'Implied internal energy flux')
       call addfld('SHFLXORI', horiz_only, 'A', 'W/m2', 'SHFLX before adding IEFLX')
       call add_default ('IEFLX', 1, ' ') 
    end if 

  end subroutine check_energy_init

!===============================================================================

  subroutine check_energy_timestep_init(state, tend, pbuf, col_type)
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field
!-----------------------------------------------------------------------
! Compute initial values of energy and water integrals, 
! zero cumulative tendencies
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state),   intent(inout)    :: state
    type(physics_tend ),   intent(inout)    :: tend
    type(physics_buffer_desc), pointer      :: pbuf(:)
    integer, optional                       :: col_type  ! Flag inidicating whether using grid or subcolumns
!---------------------------Local storage-------------------------------

    real(r8) :: ke(state%ncol)                     ! vertical integral of kinetic energy
    real(r8) :: se(state%ncol)                     ! vertical integral of static energy
    real(r8) :: wv(state%ncol)                     ! vertical integral of water (vapor)
    real(r8) :: wl(state%ncol)                     ! vertical integral of water (liquid)
    real(r8) :: wi(state%ncol)                     ! vertical integral of water (ice)
    real(r8) :: te(state%ncol)       
    real(r8) :: tw(state%ncol)

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    real(r8) :: wr(state%ncol)                     ! vertical integral of rain
    real(r8) :: ws(state%ncol)                     ! vertical integral of snow
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol

    call energy_helper_eam_def(state%u,state%v,state%T,state%q,state%ps,state%pdel,state%phis, &
                                   ke,se,wv,wl,wi,wr,ws,te,tw, &
                                   ncol)

    state%te_ini(:ncol) = te(:ncol)
    state%tw_ini(:ncol) = tw(:ncol)

    state%te_cur(:ncol) = state%te_ini(:ncol)
    state%tw_cur(:ncol) = state%tw_ini(:ncol)


! zero cummulative boundary fluxes 
    tend%te_tnd(:ncol) = 0._r8
    tend%tw_tnd(:ncol) = 0._r8

    state%count = 0

! initialize physics buffer
    if (is_first_step()) then
       call pbuf_set_field(pbuf, teout_idx, state%te_ini, col_type=col_type)
    end if

  end subroutine check_energy_timestep_init

!===============================================================================

  subroutine check_energy_chng(state, tend, name, nstep, ztodt,        &
       flx_vap, flx_cnd, flx_ice, flx_sen)

!-----------------------------------------------------------------------
! Check that the energy and water change matches the boundary fluxes
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    use cam_history,       only: outfld

    type(physics_state)    , intent(inout) :: state
    type(physics_tend )    , intent(inout) :: tend
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: flx_vap(:)          ! (pcols) - boundary flux of vapor         (kg/m2/s)
    real(r8), intent(in   ) :: flx_cnd(:)          ! (pcols) -boundary flux of liquid+ice    (m/s) (precip?)
    real(r8), intent(in   ) :: flx_ice(:)          ! (pcols) -boundary flux of ice           (m/s) (snow?)
    real(r8), intent(in   ) :: flx_sen(:)          ! (pcols) -boundary flux of sensible heat (w/m2)

!******************** BAB ******************************************************
!******* Note that the precip and ice fluxes are in precip units (m/s). ********
!******* I would prefer to have kg/m2/s.                                ********
!******* I would also prefer liquid (not total) and ice fluxes          ********
!*******************************************************************************

!---------------------------Local storage-------------------------------

    real(r8) :: te_xpd(state%ncol)                 ! expected value (f0 + dt*boundary_flux)
    real(r8) :: te_dif(state%ncol)                 ! energy of input state - original energy
    real(r8) :: te_tnd(state%ncol)                 ! tendency from last process
    real(r8) :: te_rer(state%ncol)                 ! relative error in energy column
    real(r8) :: te_err(state%ncol)                 ! absolute error in energy column

    real(r8) :: tw_xpd(state%ncol)                 ! expected value (w0 + dt*boundary_flux)
    real(r8) :: tw_dif(state%ncol)                 ! tw_inp - original water
    real(r8) :: tw_tnd(state%ncol)                 ! tendency from last process
    real(r8) :: tw_rer(state%ncol)                 ! relative error in water column
    real(r8) :: tw_err(state%ncol)                 ! absolute error in water column

    real(r8) :: ke(state%ncol)                     ! vertical integral of kinetic energy
    real(r8) :: se(state%ncol)                     ! vertical integral of static energy
    real(r8) :: wv(state%ncol)                     ! vertical integral of water (vapor)
    real(r8) :: wl(state%ncol)                     ! vertical integral of water (liquid)
    real(r8) :: wi(state%ncol)                     ! vertical integral of water (ice)

    real(r8) :: te(state%ncol)                     ! vertical integral of total energy
    real(r8) :: tw(state%ncol)                     ! vertical integral of total water

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    real(r8) :: wr(state%ncol)                     ! vertical integral of rain
    real(r8) :: ws(state%ncol)                     ! vertical integral of snow
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol

    call energy_helper_eam_def(state%u,state%v,state%T,state%q,state%ps,state%pdel,state%phis, &
                                   ke,se,wv,wl,wi,wr,ws,te,tw, &
                                   ncol)

    ! compute expected values and tendencies
    do i = 1, ncol
       ! change in static energy and total water
       te_dif(i) = te(i) - state%te_cur(i)
       tw_dif(i) = tw(i) - state%tw_cur(i)

       ! expected tendencies from boundary fluxes for last process
       te_tnd(i) = flx_vap(i)*(latvap+latice) - (flx_cnd(i) - flx_ice(i))*1000._r8*latice + flx_sen(i)
       tw_tnd(i) = flx_vap(i) - flx_cnd(i) *1000._r8

       ! cummulative tendencies from boundary fluxes
       tend%te_tnd(i) = tend%te_tnd(i) + te_tnd(i)
       tend%tw_tnd(i) = tend%tw_tnd(i) + tw_tnd(i)

       ! expected new values from previous state plus boundary fluxes
       te_xpd(i) = state%te_cur(i) + te_tnd(i)*ztodt
       tw_xpd(i) = state%tw_cur(i) + tw_tnd(i)*ztodt

       ! absolute error, expected value - input state / previous state 
       te_err(i) = te_xpd(i) - te(i)

       ! relative error, expected value - input state / previous state 
       te_rer(i) = (te_xpd(i) - te(i)) / state%te_cur(i)
    end do

    ! absolute error for total water (allow for dry atmosphere)
    tw_err = 0._r8
    tw_err(:ncol) = tw_xpd(:ncol) - tw(:ncol)

    ! relative error for total water (allow for dry atmosphere)
    tw_rer = 0._r8
    where (state%tw_cur(:ncol) > 0._r8) 
       tw_rer(:ncol) = (tw_xpd(:ncol) - tw(:ncol)) / state%tw_cur(:ncol)
    end where

    ! error checking
    ! the relative error threshold for the water budget has been reduced to 1.e-10
    ! to avoid messages generated by QNEG3 calls
    ! PJR- change to identify if error in energy or water 
    if (print_energy_errors) then
       if (any(abs(te_rer(1:ncol)) > 1.E-14_r8 )) then
          write(iulog,"(a,a,i8,i5)") "significant en_conservation errors from process, nstep, chunk ", &
                name, nstep, lchnk
          write(iulog,"(a5,2a10,6a15)") ' i', 'lat', 'lon', 'en', 'en_from_flux', 'diff', 'exptd diff', 'rerr', 'cum_diff'
	  do i = 1,ncol
             if (abs(te_rer(i)) > 1.E-14_r8 ) then 
                state%count = state%count + 1
                write(iulog,"(i5,2f10.2,6e15.7)") i, state%lat(i), state%lon(i), te(i),te_xpd(i),te_dif(i),  &
                      te_tnd(i)*ztodt,te_rer(i), tend%te_tnd(i)*ztodt
             endif
          end do
       end if

       if (any(abs(tw_rer(1:ncol)) > 1.E-10_r8)) then
          write(iulog,"(a,a,i8,i5)") "significant w_conservation errors from process, nstep, chunk ", &
               name, nstep, lchnk
          write(iulog,"(a5,2a10,6a15)") ' i', 'lat', 'lon', 'tw', 'tw_from_flux', 'diff', 'exptd diff', 'rerr', 'cum_diff'
          do i = 1, ncol
             if ( abs(tw_rer(i)) > 1.E-10_r8) then
                state%count = state%count + 1
                write(iulog,"(i5,2f10.2,6e15.7)") i, state%lat(i), state%lon(i),tw(i),tw_xpd(i),tw_dif(i), &
                     tw_tnd(i)*ztodt, tw_rer(i), tend%tw_tnd(i)*ztodt
             end if
          end do
       end if
    end if

    ! copy new value to state
    do i = 1, ncol
       state%te_cur(i) = te(i)
       state%tw_cur(i) = tw(i)
    end do

  end subroutine check_energy_chng


!===============================================================================
  subroutine check_energy_gmean(state, pbuf2d, dtime, nstep)

    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    
!-----------------------------------------------------------------------
! Compute global mean total energy of physics input and output states
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
    type(physics_buffer_desc),    pointer    :: pbuf2d(:,:)

    real(r8), intent(in) :: dtime        ! physics time step
    integer , intent(in) :: nstep        ! current timestep number

!---------------------------Local storage-------------------------------
    integer :: ncol                      ! number of active columns
    integer :: lchnk                     ! chunk index

    real(r8) :: te(pcols,begchunk:endchunk,3)   
                                         ! total energy of input/output states (copy)
    real(r8) :: te_glob(3)               ! global means of total energy
    real(r8), pointer :: teout(:)
!-----------------------------------------------------------------------

    ! Copy total energy out of input and output states
#ifdef CPRCRAY
!DIR$ CONCURRENT
#endif
    do lchnk = begchunk, endchunk
       ncol = state(lchnk)%ncol
       ! input energy
       te(:ncol,lchnk,1) = state(lchnk)%te_ini(:ncol)
       ! output energy
       call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk),teout_idx, teout)

       te(:ncol,lchnk,2) = teout(1:ncol)
       ! surface pressure for heating rate
       te(:ncol,lchnk,3) = state(lchnk)%pint(:ncol,pver+1)
    end do

    ! Compute global means of input and output energies and of
    ! surface pressure for heating rate (assume uniform ptop)
    call gmean(te, te_glob, 3)

    if (begchunk .le. endchunk) then
       teinp_glob = te_glob(1)
       teout_glob = te_glob(2)
       psurf_glob = te_glob(3)
       ptopb_glob = state(begchunk)%pint(1,1)

       ! Global mean total energy difference
       tedif_glob =  teinp_glob - teout_glob
       heat_glob  = -tedif_glob/dtime * gravit / (psurf_glob - ptopb_glob)

       if (masterproc) then
          write(iulog,'(1x,a9,1x,i8,4(1x,e25.17))') "nstep, te", nstep, teinp_glob, teout_glob, heat_glob, psurf_glob
       end if
    else
       heat_glob = 0._r8
    end if  !  (begchunk .le. endchunk)
    
  end subroutine check_energy_gmean

!===============================================================================

subroutine ieflx_gmean(state, tend, pbuf2d, cam_in, cam_out, nstep)

!!...................................................................
!! Calculate global mean of the implied internal energy flux 
!! 
!! This subroutien is called only when ieflx_opt > 0 
!! 
!! Author: Kai Zhang 
!!...................................................................

    use camsrfexch,       only: cam_out_t, cam_in_t
    use physics_buffer,   only: physics_buffer_desc, pbuf_get_field, pbuf_get_chunk, pbuf_set_field 
    use cam_history,      only: outfld
    use phys_control,     only: ieflx_opt

    integer , intent(in) :: nstep        ! current timestep number
    type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
    type(physics_tend ), intent(in   ), dimension(begchunk:endchunk) :: tend
    type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
    type(cam_out_t),                    dimension(begchunk:endchunk) :: cam_out
    type(physics_buffer_desc),          pointer :: pbuf2d(:,:)

    real(r8), parameter :: cpsw = 3.996e3  !! cp of sea water used in MPAS [J/kg/K]
    real(r8), parameter :: rhow = 1.e3     !! density of water [kg/m3]
 
    integer :: ncol                      ! number of active columns
    integer :: lchnk                     ! chunk index

    real(r8) :: ieflx(pcols) 

    real(r8) :: qflx(pcols,begchunk:endchunk) !qflx [kg/m2/s]
    real(r8) :: rain(pcols,begchunk:endchunk) !rain [m/s] 
    real(r8) :: snow(pcols,begchunk:endchunk) !snow [m/s] 
    real(r8) :: ienet(pcols,begchunk:endchunk) !ieflx net [W/m2] or [J/m2/s]

!- 
    ieflx_glob = 0._r8

    qflx = 0._r8 
    rain = 0._r8 
    snow = 0._r8 
    ienet = 0._r8 

#ifdef CPRCRAY
!DIR$ CONCURRENT
#endif
    do lchnk = begchunk, endchunk

       ncol = state(lchnk)%ncol
       qflx(:ncol,lchnk) = cam_in(lchnk)%cflx(:ncol,1)
       snow(:ncol,lchnk) = cam_out(lchnk)%precsc(:ncol) + cam_out(lchnk)%precsl(:ncol)
       rain(:ncol,lchnk) = cam_out(lchnk)%precc(:ncol)  + cam_out(lchnk)%precl(:ncol) - snow(:ncol,lchnk) 

       select case (ieflx_opt) 

       !!..................................................................................... 
       !! Calculate the internal energy flux at surface (imitate what is considered in the ocean model)   
       !! 
       !! ieflx_opt = 1 : air temperature in the lowest model layer will be used 
       !! ieflx_opt = 2 : skin temperature (from lnd/ocn/ice components) will be used  
       !! 
       !! ieflx_opt = 2 is recommended for now. 
       !! 
       !! (rhow*) converts the unit of precipitation from m/s to kg/m2/s 
       !!..................................................................................... 

       case(1) 
          ienet(:ncol,lchnk) = cpsw * qflx(:ncol,lchnk) * cam_in(lchnk)%ts(:ncol) - & 
                               cpsw * rhow * ( rain(:ncol,lchnk) + snow(:ncol,lchnk) ) * cam_out(lchnk)%tbot(:ncol)
       case(2) 
          ienet(:ncol,lchnk) = cpsw * qflx(:ncol,lchnk) * cam_in(lchnk)%ts(:ncol) - & 
                               cpsw * rhow * ( rain(:ncol,lchnk) + snow(:ncol,lchnk) ) * cam_in(lchnk)%ts(:ncol)
       case(3) 
          ienet(:ncol,lchnk) = cam_in(lchnk)%h2otemp(:ncol)
       case default 
          call endrun('*** incorrect ieflx_opt ***')
       end select 


    end do

    call gmean(ienet, ieflx_glob)

#ifdef CPRCRAY
!DIR$ CONCURRENT
#endif
    do lchnk = begchunk, endchunk

       ieflx = ieflx_glob

       call outfld('IEFLX', ieflx, pcols, lchnk)

    end do

!!!    if (begchunk .le. endchunk) then
!!!       if (masterproc) then
!!!          write(iulog,'(1x,a12,1x,i8,4(1x,e25.17))') "nstep, ieflx, ieup, iedn ", nstep, ieflx_glob, ieup_glob, iedn_glob 
!!!       end if
!!!    end if

  end subroutine ieflx_gmean


!===============================================================================
  subroutine check_ieflx_fix(lchnk, ncol, nstep, shflx)

!!
!! Add the global mean internal energy flux to the sensible heat flux 
!!
!! This subroutien is called only when ieflx_opt > 0 
!! 
!! Called by typhsac 
!! 

    use cam_history,       only: outfld

    integer, intent(in   ) :: nstep          ! time step number
    integer, intent(in   ) :: lchnk  
    integer, intent(in   ) :: ncol
    real(r8),intent(inout) :: shflx(pcols) 

    integer :: i

    call outfld('SHFLXORI', shflx, pcols, lchnk)

    if(nstep>1) then 
       do i = 1, ncol
          shflx(i) = shflx(i) + ieflx_glob 
       end do
    end if 

    return
  end subroutine check_ieflx_fix

!===============================================================================

subroutine qflx_gmean(state, tend, cam_in, dtime, nstep)

!!...................................................................
!! Calculate global mean of qflx for water conservation check 
!! 
!! Author: Kai Zhang 
!!...................................................................

    use camsrfexch,       only: cam_out_t, cam_in_t
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk

    ! Compute global mean qflx

    type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
    type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
    type(physics_buffer_desc),    pointer    :: pbuf2d(:,:)
    type(physics_tend ), intent(in   ), dimension(begchunk:endchunk) :: tend

    real(r8), intent(in) :: dtime        ! physics time step
    integer , intent(in) :: nstep        ! current timestep number

    integer :: ncol                      ! number of active columns
    integer :: lchnk                     ! chunk index

    real(r8) :: qflx(pcols,begchunk:endchunk)
    real(r8) :: qflx_glob

    qflx_glob = 0._r8 

#ifdef CPRCRAY
!DIR$ CONCURRENT
#endif
    do lchnk = begchunk, endchunk
       ncol = state(lchnk)%ncol
       qflx(:ncol,lchnk) = cam_in(lchnk)%cflx(:ncol,1)
    end do

    call gmean(qflx, qflx_glob)

    if (begchunk .le. endchunk) then
       if (masterproc) then
          write(iulog,'(1x,a12,1x,i8,4(1x,e25.17))') "nstep, qflx ", nstep, qflx_glob
       end if
    end if

  end subroutine qflx_gmean


!===============================================================================
  subroutine check_energy_fix(state, ptend, nstep, eshflx)

!-----------------------------------------------------------------------
! Add heating rate required for global mean total energy conservation
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    use cam_history, only: outfld
    use scamMod, only: heat_glob_scm, single_column, use_replay

    type(physics_state), intent(in   ) :: state
    type(physics_ptend), intent(out)   :: ptend

    integer , intent(in   ) :: nstep          ! time step number
    real(r8), intent(out  ) :: eshflx(pcols)  ! effective sensible heat flux
    real(r8) :: heat_out(pcols)

!---------------------------Local storage-------------------------------
    integer  :: i                        ! column
    integer  :: ncol                     ! number of atmospheric columns in chunk
    integer  :: lchnk
!-----------------------------------------------------------------------
    ncol = state%ncol
    lchnk = state%lchnk

    call physics_ptend_init(ptend, state%psetcols, 'chkenergyfix', ls=.true.)

#if ( defined OFFLINE_DYN )
    ! disable the energy fix for offline driver
    heat_glob = 0._r8
#endif
! add (-) global mean total energy difference as heating
    if (single_column .and. use_replay) then
      heat_glob = heat_glob_scm(1)
    endif
    
    ! In single column model we do NOT want to take into
    !   consideration the dynamics energy fixer.  Since only
    !   one column of dynamics is active, this data will 
    !   essentially be garbage. 
    if (single_column .and. .not. use_replay) then
      heat_glob = 0._r8
    endif
    
    ptend%s(:ncol,:pver) = heat_glob
!!$    write(iulog,*) "chk_fix: heat", state%lchnk, ncol, heat_glob

#if ( defined E3SM_SCM_REPLAY )
    if (nstep > 0) then
      heat_out(:ncol) = heat_glob
      call outfld('heat_glob',  heat_out(:ncol), pcols, lchnk)
    endif
#endif

! compute effective sensible heat flux
    do i = 1, ncol
       eshflx(i) = heat_glob * (state%pint(i,pver+1) - state%pint(i,1)) / gravit
    end do
!!!    if (nstep > 0) write(iulog,*) "heat", heat_glob, eshflx(1)

    return
  end subroutine check_energy_fix


!===============================================================================
  subroutine check_tracers_init(state, tracerint)

!-----------------------------------------------------------------------
! Compute initial values of tracers integrals, 
! zero cumulative tendencies
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------

    type(physics_state),   intent(in)    :: state
    type(check_tracers_data), intent(out)   :: tracerint

!---------------------------Local storage-------------------------------

    real(r8) :: tr(pcols)                          ! vertical integral of tracer
    real(r8) :: trpdel(pcols, pver)                ! pdel for tracer

    integer ncol                                   ! number of atmospheric columns
    integer ierror                                 ! allocate status return
    integer  i,k,m                                 ! column, level,constituent indices

!-----------------------------------------------------------------------

    ncol  = state%ncol
    allocate (tracerint%tracer(pcols,pcnst), stat=ierror)
    if ( ierror /= 0 ) call endrun('CHECK_TRACERS_INIT error: allocation error tracer')

    allocate (tracerint%tracer_tnd(pcols,pcnst), stat=ierror)
    if ( ierror /= 0 ) call endrun('CHECK_TRACERS_INIT error: allocation error tracer_tnd')

    do m = 1,pcnst

       if ( any(m == (/ 1, icldliq, icldice, &
                           irain,   isnow    /)) ) exit   ! dont process water substances
                                                            ! they are checked in check_energy
       if (cnst_get_type_byind(m).eq.'dry') then
          trpdel(:ncol,:) = state%pdeldry(:ncol,:)
       else
          trpdel(:ncol,:) = state%pdel(:ncol,:)
       endif

       ! Compute vertical integrals of tracer
       tr = 0._r8
       do k = 1, pver
          do i = 1, ncol
             tr(i) = tr(i) + state%q(i,k,m)*trpdel(i,k)/gravit
          end do
       end do

       ! Compute vertical integrals of frozen static tracers and total water.
       do i = 1, ncol
          tracerint%tracer(i,m) = tr(i)
       end do

       ! zero cummulative boundary fluxes 
       tracerint%tracer_tnd(:ncol,m) = 0._r8

       tracerint%count(m) = 0

    end do

    return
  end subroutine check_tracers_init

!===============================================================================
  subroutine check_tracers_fini(tracerint)

!-----------------------------------------------------------------------
! Deallocate storage assoicated with check_tracers_data type variable
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------

    type(check_tracers_data), intent(inout)   :: tracerint

!-----------------------------------------------------------------------

    deallocate(tracerint%tracer)
    deallocate(tracerint%tracer_tnd)

    return
  end subroutine check_tracers_fini

!===============================================================================
  subroutine check_tracers_chng(state, tracerint, name, nstep, ztodt, cflx)

!-----------------------------------------------------------------------
! Check that the tracers and water change matches the boundary fluxes
! these checks are not save when there are tracers transformations, as 
! they only check to see whether a mass change in the column is
! associated with a flux
!-----------------------------------------------------------------------

    implicit none

!------------------------------Arguments--------------------------------

    type(physics_state)    , intent(in   ) :: state
    type(check_tracers_data), intent(inout) :: tracerint! tracers integrals and boundary fluxes
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: cflx(pcols,pcnst)       ! boundary flux of tracers       (kg/m2/s)

!---------------------------Local storage-------------------------------

    real(r8) :: tracer_inp(pcols,pcnst)                   ! total tracer of new (input) state
    real(r8) :: tracer_xpd(pcols,pcnst)                   ! expected value (w0 + dt*boundary_flux)
    real(r8) :: tracer_dif(pcols,pcnst)                   ! tracer_inp - original tracer
    real(r8) :: tracer_tnd(pcols,pcnst)                   ! tendency from last process
    real(r8) :: tracer_rer(pcols,pcnst)                   ! relative error in tracer column

    real(r8) :: tr(pcols)                           ! vertical integral of tracer
    real(r8) :: trpdel(pcols, pver)                       ! pdel for tracer

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: m                            ! tracer index
    character(len=8) :: tracname   ! tracername
!-----------------------------------------------------------------------
!!$    if (.true.) return

    lchnk = state%lchnk
    ncol  = state%ncol

    do m = 1,pcnst

       if ( any(m == (/ 1, icldliq, icldice, &
                           irain,   isnow    /)) ) exit   ! dont process water substances
                                                            ! they are checked in check_energy

       tracname = cnst_name(m)
       if (cnst_get_type_byind(m).eq.'dry') then
          trpdel(:ncol,:) = state%pdeldry(:ncol,:)
       else
          trpdel(:ncol,:) = state%pdel(:ncol,:)
       endif

       ! Compute vertical integrals tracers
       tr = 0._r8
       do k = 1, pver
          do i = 1, ncol
             tr(i) = tr(i) + state%q(i,k,m)*trpdel(i,k)/gravit
          end do
       end do

       ! Compute vertical integrals of tracer
       do i = 1, ncol
          tracer_inp(i,m) = tr(i)
       end do

       ! compute expected values and tendencies
       do i = 1, ncol
          ! change in tracers 
          tracer_dif(i,m) = tracer_inp(i,m) - tracerint%tracer(i,m)

          ! expected tendencies from boundary fluxes for last process
          tracer_tnd(i,m) = cflx(i,m)

          ! cummulative tendencies from boundary fluxes
          tracerint%tracer_tnd(i,m) = tracerint%tracer_tnd(i,m) + tracer_tnd(i,m)

          ! expected new values from original values plus boundary fluxes
          tracer_xpd(i,m) = tracerint%tracer(i,m) + tracerint%tracer_tnd(i,m)*ztodt

          ! relative error, expected value - input value / original 
          tracer_rer(i,m) = (tracer_xpd(i,m) - tracer_inp(i,m)) / tracerint%tracer(i,m)
       end do

!! final loop for error checking
!    do i = 1, ncol

!! error messages
!       if (abs(enrgy_rer(i)) > 1.E-14 .or. abs(water_rer(i)) > 1.E-14) then
!          tracerint%count = tracerint%count + 1
!          write(iulog,*) "significant conservations error after ", name,        &
!               " count", tracerint%count, " nstep", nstep, "chunk", lchnk, "col", i
!          write(iulog,*) enrgy_inp(i),enrgy_xpd(i),enrgy_dif(i),tracerint%enrgy_tnd(i)*ztodt,  &
!               enrgy_tnd(i)*ztodt,enrgy_rer(i)
!          write(iulog,*) water_inp(i),water_xpd(i),water_dif(i),tracerint%water_tnd(i)*ztodt,  &
!               water_tnd(i)*ztodt,water_rer(i)
!       end if
!    end do


       ! final loop for error checking
       if ( maxval(tracer_rer) > 1.E-14_r8 ) then
          write(iulog,*) "CHECK_TRACERS TRACER large rel error"
          write(iulog,*) tracer_rer
       endif

       do i = 1, ncol
          ! error messages
          if (abs(tracer_rer(i,m)) > 1.E-14_r8 ) then
             tracerint%count = tracerint%count + 1
             write(iulog,*) "CHECK_TRACERS TRACER significant conservation error after ", name,        &
                  " count", tracerint%count, " nstep", nstep, "chunk", lchnk, "col",i
             write(iulog,*)' process name, tracname, index ',  name, tracname, m
             write(iulog,*)" input integral              ",tracer_inp(i,m)
             write(iulog,*)" expected integral           ", tracer_xpd(i,m)
             write(iulog,*)" input - inital integral     ",tracer_dif(i,m)
             write(iulog,*)" cumulative tend      ",tracerint%tracer_tnd(i,m)*ztodt
             write(iulog,*)" process tend         ",tracer_tnd(i,m)*ztodt
             write(iulog,*)" relative error       ",tracer_rer(i,m)
             call endrun()
          end if
       end do
    end do

    return
  end subroutine check_tracers_chng


  subroutine check_water(state, tend, name, nstep, ztodt)

!!...................................................................
!! Check water path at certain locations for water conservation check 
!! 
!! Water species include water vapor, droplet, ice, rain, and snow 
!!
!! Author: Kai Zhang 
!!...................................................................

!! Arguments 
!!...................................................................
    use cam_history,       only: outfld

    type(physics_state)    , intent(in) :: state
    type(physics_tend )    , intent(in) :: tend
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)

!! Local 
!!...................................................................

    real(r8) :: ke(state%ncol) 
    real(r8) :: se(state%ncol) 
    real(r8) :: te(state%ncol) 
    real(r8) :: wv(state%ncol)                     ! vertical integral of water (vapor)
    real(r8) :: wl(state%ncol)                     ! vertical integral of water (liquid)
    real(r8) :: wi(state%ncol)                     ! vertical integral of water (ice)
    real(r8) :: tw(state%ncol)                     ! vertical integral of total water

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    real(r8) :: wr(state%ncol)                     ! vertical integral of rain
    real(r8) :: ws(state%ncol)                     ! vertical integral of snow
!!...................................................................

    lchnk = state%lchnk
    ncol  = state%ncol

    call energy_helper_eam_def(state%u,state%v,state%T,state%q,state%ps,state%pdel,state%phis, &
                                   ke,se,wv,wl,wi,wr,ws,te,tw, &
                                   ncol)

    if(name.eq.'PHYBC01') then 
       call outfld('BC01Q',           wv,pcols   ,lchnk   )
       call outfld('BC01QL',          wl,pcols   ,lchnk   )
       call outfld('BC01QI',          wi,pcols   ,lchnk   )
       call outfld('BC01QR',          wr,pcols   ,lchnk   )
       call outfld('BC01QS',          ws,pcols   ,lchnk   )
       call outfld('BC01TW',          tw,pcols   ,lchnk   )
    end if

    if(name.eq.'PHYBC02') then
       call outfld('BC02Q',           wv,pcols   ,lchnk   )
       call outfld('BC02QL',          wl,pcols   ,lchnk   )
       call outfld('BC02QI',          wi,pcols   ,lchnk   )
       call outfld('BC02QR',          wr,pcols   ,lchnk   )
       call outfld('BC02QS',          ws,pcols   ,lchnk   )
       call outfld('BC02TW',          tw,pcols   ,lchnk   )
    end if

  end subroutine check_water 




  subroutine check_qflx(state, tend, name, nstep, ztodt, qflx)

!!...................................................................
!! Output qflx at certain locations for water conservation check 
!!...................................................................

    use cam_history,       only: outfld
    
    type(physics_state)    , intent(in) :: state
    type(physics_tend )    , intent(in) :: tend
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: qflx(:)             ! (pcols) - boundary flux of vapor (kg/m2/s)

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    
!!...................................................................

    lchnk = state%lchnk
    ncol  = state%ncol
    
    if(name.eq.'PHYAC01') then
       call outfld('AC01QFLX', qflx, pcols, lchnk)
    end if

    if(name.eq.'PHYAC02') then
       call outfld('AC02QFLX', qflx, pcols, lchnk)
    end if

    if(name.eq.'PHYBC01') then
       call outfld('BC01QFLX', qflx, pcols, lchnk)
    end if

  end subroutine check_qflx




  subroutine check_prect(state, tend, name, nstep, ztodt, prect)

!!...................................................................
!! Output precipitation at certain locations for water conservation check 
!!...................................................................
    use cam_history,       only: outfld

    type(physics_state)    , intent(in) :: state
    type(physics_tend )    , intent(in) :: tend
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: prect(:)            ! (pcols) - precipitation

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices

!!...................................................................

    lchnk = state%lchnk
    ncol  = state%ncol

    if(name.eq.'PHYBC01') then
       call outfld('BC01PR', prect, pcols, lchnk)
    end if

    if(name.eq.'PHYBC02') then
       call outfld('BC02PR', prect, pcols, lchnk)
    end if

  end subroutine check_prect

!====================================================================

  subroutine energy_helper_eam_def(u,v,T,q,ps,pdel,phis, &
                                   ke,se,wv,wl,wi,wr,ws,te,tw, &     
                                   ncol,teloc,psterm)

!state vars are of size psetcols,pver, so, not exactly correct
    real(r8), intent(in) :: u(pcols,pver) 
    real(r8), intent(in) :: v(pcols,pver) 
    real(r8), intent(in) :: T(pcols,pver) 
    real(r8), intent(in) :: q(pcols,pver,pcnst) 
    real(r8), intent(in) :: ps(pcols) 
    real(r8), intent(in) :: pdel(pcols,pver) 
    real(r8), intent(in) :: phis(pcols) 


    real(r8), intent(inout) :: ke(ncol)     ! vertical integral of kinetic energy
    real(r8), intent(inout) :: se(ncol)     ! vertical integral of static energy
    real(r8), intent(inout) :: wv(ncol)     ! vertical integral of water (vapor)
    real(r8), intent(inout) :: wl(ncol)     ! vertical integral of water (liquid)
    real(r8), intent(inout) :: wi(ncol)     ! vertical integral of water (ice)
    real(r8), intent(inout) :: te(ncol)     ! vertical integral of total energy
    real(r8), intent(inout) :: tw(ncol)     ! vertical integral of total water
    real(r8), intent(inout) :: wr(ncol)     ! vertical integral of rain
    real(r8), intent(inout) :: ws(ncol)     ! vertical integral of snow

! do not use in this version
    real(r8), intent(inout), optional :: teloc(pcols,pver) 
    real(r8), intent(inout), optional :: psterm(pcols) 

    integer, intent(in) :: ncol                   
    integer :: i,k                            

    !if (icldliq > 1  .and.  icldice > 1 .and. irain > 1 .and. isnow > 1) then
       do i = 1, ncol
          call energy_helper_eam_def_column(u(i,:),v(i,:),T(i,:),q(i,1:pver,1:pcnst),&
                                   ps(i),pdel(i,:),phis(i), &
                                   ke(i),se(i),wv(i),wl(i),wi(i),wr(i),ws(i),te(i),tw(i) )                             
       enddo
    !else
    !   call endrun('energy_helper...column is not implemented if water forms do not exist')
    !endif

  end subroutine energy_helper_eam_def

  subroutine energy_helper_eam_def_column(u,v,T,q,ps,pdel,phis, &
                                   ke,se,wv,wl,wi,wr,ws,te,tw, &
                                   teloc,psterm)

!state vars are of size psetcols,pver, so, not exactly correct
    real(r8), intent(in) :: u(pver)
    real(r8), intent(in) :: v(pver)
    real(r8), intent(in) :: T(pver)
    real(r8), intent(in) :: q(pver,pcnst)
    real(r8), intent(in) :: ps
    real(r8), intent(in) :: pdel(pver)
    real(r8), intent(in) :: phis

    real(r8), intent(inout) :: ke     ! vertical integral of kinetic energy
    real(r8), intent(inout) :: se     ! vertical integral of static energy
    real(r8), intent(inout) :: wv     ! vertical integral of water (vapor)
    real(r8), intent(inout) :: wl     ! vertical integral of water (liquid)
    real(r8), intent(inout) :: wi     ! vertical integral of water (ice)
    real(r8), intent(inout) :: te     ! vertical integral of total energy
    real(r8), intent(inout) :: tw     ! vertical integral of total water
    real(r8), intent(inout) :: wr     ! vertical integral of rain
    real(r8), intent(inout) :: ws     ! vertical integral of snow

    real(r8), intent(inout), optional :: teloc(pver)
    real(r8), intent(inout), optional :: psterm

    integer :: i,k

    ! Compute vertical integrals of dry static energy and water (vapor, liquid, ice)
    ke = 0._r8
    se = 0._r8
    wv = 0._r8
    wl = 0._r8
    wi = 0._r8
    wr = 0._r8
    ws = 0._r8

    !keep it bfb and fast
    if (present(teloc) .and. present(psterm))then
       teloc = 0.0; psterm = 0.0
       do k = 1, pver
          teloc(k) = 0.5_r8*(u(k)**2 + v(k)**2)*pdel(k)/gravit &
                   + t(k)*cpair*pdel(k)/gravit &
                   + (latvap+latice)*q(k,1       )*pdel(k)/gravit
       end do
       if (icldliq > 1 .and. irain > 1) then
          do k = 1, pver
             teloc(k) = teloc(k) &
                      + latice*(q(k,icldliq) + q(k,irain))*pdel(k)/gravit
          end do
       end if
       psterm = phis*ps/gravit
    endif

    do k = 1, pver
       ke = ke + 0.5_r8*(u(k)**2 + v(k)**2)*pdel(k)/gravit
       se = se +         t(k)*cpair*pdel(k)/gravit
       wv = wv + q(k,1      )*pdel(k)/gravit
    end do
    se = se + phis*ps/gravit

    if (icldliq > 1 .and. icldice > 1) then
       do k = 1, pver
          wl = wl + q(k,icldliq)*pdel(k)/gravit
          wi = wi + q(k,icldice)*pdel(k)/gravit
       end do
    end if

    if (microp_scheme == 'P3') then
       ! In the case where the micro-physics scheme is P3 there is no snow
       ! constituent and isnow = -1.  So we still calculate the rest of the
       ! consituents and ensure that ws = 0.0.  NOTE! This change will likely
       ! lead to conflicts with upstream E3SM any time check_energy is
       ! changed...  The most important thing is to avoid any SNOW calculations
       ! when P3 is the microphysics scheme.
       ws = 0.0_r8

       if (irain > 1) then
          do k = 1, pver
             wr = wr + q(k,irain)*pdel(k)/gravit
          end do
       end if
    else if (irain > 1 .and. isnow > 1) then
       do k = 1, pver
          wr = wr + q(k,irain)*pdel(k)/gravit
          ws = ws + q(k,isnow)*pdel(k)/gravit
       end do
    end if

    ! Compute vertical integrals of frozen static energy and total water.
    te = se + ke + (latvap+latice)*wv + latice*( wl + wr )
    tw = wv + wl + wi + wr + ws

  end subroutine energy_helper_eam_def_column



end module check_energy

