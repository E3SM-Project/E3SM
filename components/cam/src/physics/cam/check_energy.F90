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

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver, begchunk, endchunk
  use spmd_utils,      only: masterproc
  
  use phys_gmean,      only: gmean
  use physconst,       only: gravit, latvap, latice, cpair, cpairv
  use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_ptend_init
  use constituents,    only: cnst_get_ind, pcnst, cnst_name, cnst_get_type_byind
  use time_manager,    only: is_first_step
  use cam_logfile,     only: iulog
  use cam_abortutils,  only: endrun 
  use phys_control,    only: ieflx_opt !!l_ieflx_fix

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

  public :: qflx_gmean              ! calculate global mean of qflx for water conservation check 
  public :: check_qflx              ! output qflx at certain locations for water conservation check  
  public :: check_prect             ! output prect at certain locations for water conservation check  
  public :: check_water             ! output water path at certain locations for water conservation check  

  public :: ieflx_gmean             ! calculate global mean of ieflx 
  public :: check_ieflx_fix         ! add ieflx to sensible heat flux 


! Private module data

  logical  :: print_energy_errors = .false.

  real(r8) :: teout_glob           ! global mean energy of output state
  real(r8) :: teinp_glob           ! global mean energy of input state
  real(r8) :: tedif_glob           ! global mean energy difference
  real(r8) :: psurf_glob           ! global mean surface pressure
  real(r8) :: ptopb_glob           ! global mean top boundary pressure
  real(r8) :: heat_glob            ! global mean heating rate
  real(r8) :: ieflx_glob           ! global mean implied internal energy flux 

! Physics buffer indices
  
  integer  :: ieflx_idx  = 0       ! teout index in physics buffer 
  integer  :: teout_idx  = 0       ! teout index in physics buffer 
  integer  :: dtcore_idx = 0       ! dtcore index in physics buffer 

  type check_tracers_data
     real(r8) :: tracer(pcols,pcnst)       ! initial vertically integrated total (kinetic + static) energy
     real(r8) :: tracer_tnd(pcols,pcnst)   ! cumulative boundary flux of total energy
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
    call pbuf_add_field('IEFLX', 'global',dtype_r8 , (/pcols,dyn_time_lvls/),      ieflx_idx)
    call pbuf_add_field('DTCORE','global',dtype_r8,  (/pcols,pver,dyn_time_lvls/),dtcore_idx)

    if(is_subcol_on()) then
      call pbuf_register_subcol('TEOUT', 'phys_register', teout_idx)
      call pbuf_register_subcol('IEFLX', 'phys_register', ieflx_idx)
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
       call addfld('SHFLXFIX', horiz_only, 'A', 'W/m2', 'SHFLX after adding IEFLX')
       call add_default ('SHFLXFIX', 1, ' ') 
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

    real(r8) :: ieflx(state%ncol)                     ! vertical integral of kinetic energy

    real(r8),allocatable :: cpairv_loc(:,:,:)

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    real(r8) :: wr(state%ncol)                     ! vertical integral of rain
    real(r8) :: ws(state%ncol)                     ! vertical integral of snow
    integer :: ixrain
    integer :: ixsnow
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain, abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow, abort=.false.)

    ! cpairv_loc needs to be allocated to a size which matches state and ptend
    ! If psetcols == pcols, cpairv is the correct size and just copy into cpairv_loc
    ! If psetcols > pcols and all cpairv match cpair, then assign the constant cpair

    if (state%psetcols == pcols) then
       allocate (cpairv_loc(state%psetcols,pver,begchunk:endchunk))
       cpairv_loc(:,:,:) = cpairv(:,:,:)
    else if (state%psetcols > pcols .and. all(cpairv(:,:,:) == cpair)) then
       allocate(cpairv_loc(state%psetcols,pver,begchunk:endchunk))
       cpairv_loc(:,:,:) = cpair
    else
       call endrun('check_energy_timestep_init: cpairv is not allowed to vary when subcolumns are turned on')
    end if

! Compute vertical integrals of dry static energy and water (vapor, liquid, ice)
    ke = 0._r8
    se = 0._r8
    wv = 0._r8
    wl = 0._r8
    wi = 0._r8
    wr = 0._r8
    ws = 0._r8

    ieflx = 0._r8

    do k = 1, pver
       do i = 1, ncol
          ke(i) = ke(i) + 0.5_r8*(state%u(i,k)**2 + state%v(i,k)**2)*state%pdel(i,k)/gravit
          se(i) = se(i) + state%s(i,k         )*state%pdel(i,k)/gravit
!!! cam6  se(i) = se(i) +         state%t(i,k)*cpairv_loc(i,k,lchnk)*state%pdel(i,k)/gravit
          wv(i) = wv(i) + state%q(i,k,1       )*state%pdel(i,k)/gravit
       end do
    end do
!!! cam6    do i = 1, ncol
!!! cam6       se(i) = se(i) + state%phis(i)*state%ps(i)/gravit
!!! cam6    end do

    ! Don't require cloud liq/ice to be present.  Allows for adiabatic/ideal phys.
    if (ixcldliq > 1  .and.  ixcldice > 1) then
       do k = 1, pver
          do i = 1, ncol
             wl(i) = wl(i) + state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit
             wi(i) = wi(i) + state%q(i,k,ixcldice)*state%pdel(i,k)/gravit
          end do
       end do
    end if

    if (ixrain   > 1  .and.  ixsnow   > 1 ) then
       do k = 1, pver
          do i = 1, ncol
             wr(i) = wr(i) + state%q(i,k,ixrain)*state%pdel(i,k)/gravit
             ws(i) = ws(i) + state%q(i,k,ixsnow)*state%pdel(i,k)/gravit
          end do
       end do
    end if


! Compute vertical integrals of frozen static energy and total water.
    do i = 1, ncol
!!     state%te_ini(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
       state%te_ini(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*( wl(i) + wr(i) ) 
       state%tw_ini(i) = wv(i) + wl(i) + wi(i) + wr(i) + ws(i) 

       state%te_cur(i) = state%te_ini(i)
       state%tw_cur(i) = state%tw_ini(i)
    end do

! zero cummulative boundary fluxes 
    tend%te_tnd(:ncol) = 0._r8
    tend%tw_tnd(:ncol) = 0._r8

    state%count = 0

! initialize physics buffer
    if (is_first_step()) then
       call pbuf_set_field(pbuf, teout_idx, state%te_ini, col_type=col_type)
       call pbuf_set_field(pbuf, ieflx_idx, ieflx,        col_type=col_type)
    end if

    deallocate(cpairv_loc)

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

    real(r8),allocatable :: cpairv_loc(:,:,:)

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    real(r8) :: wr(state%ncol)                     ! vertical integral of rain
    real(r8) :: ws(state%ncol)                     ! vertical integral of snow
    integer :: ixrain
    integer :: ixsnow
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain, abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow, abort=.false.)

    ! cpairv_loc needs to be allocated to a size which matches state and ptend
    ! If psetcols == pcols, cpairv is the correct size and just copy into cpairv_loc
    ! If psetcols > pcols and all cpairv match cpair, then assign the constant cpair

    if (state%psetcols == pcols) then
       allocate (cpairv_loc(state%psetcols,pver,begchunk:endchunk))
       cpairv_loc(:,:,:) = cpairv(:,:,:)
    else if (state%psetcols > pcols .and. all(cpairv(:,:,:) == cpair)) then
       allocate(cpairv_loc(state%psetcols,pver,begchunk:endchunk))
       cpairv_loc(:,:,:) = cpair
    else
       call endrun('check_energy_chng: cpairv is not allowed to vary when subcolumns are turned on')
    end if

    ! Compute vertical integrals of dry static energy and water (vapor, liquid, ice)
    ke = 0._r8
    se = 0._r8
    wv = 0._r8
    wl = 0._r8
    wi = 0._r8
    wr = 0._r8
    ws = 0._r8

    do k = 1, pver
       do i = 1, ncol
          ke(i) = ke(i) + 0.5_r8*(state%u(i,k)**2 + state%v(i,k)**2)*state%pdel(i,k)/gravit
          se(i) = se(i) + state%s(i,k         )*state%pdel(i,k)/gravit
!!!cam6   se(i) = se(i) + state%t(i,k)*cpairv_loc(i,k,lchnk)*state%pdel(i,k)/gravit
          wv(i) = wv(i) + state%q(i,k,1       )*state%pdel(i,k)/gravit
       end do
    end do
!!!cam6    do i = 1, ncol
!!!cam6       se(i) = se(i) + state%phis(i)*state%ps(i)/gravit
!!!cam6    end do

    ! Don't require cloud liq/ice to be present.  Allows for adiabatic/ideal phys.
    if (ixcldliq > 1  .and.  ixcldice > 1) then
       do k = 1, pver
          do i = 1, ncol
             wl(i) = wl(i) + state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit
             wi(i) = wi(i) + state%q(i,k,ixcldice)*state%pdel(i,k)/gravit
          end do
       end do
    end if

    if (ixrain   > 1  .and.  ixsnow   > 1 ) then
       do k = 1, pver
          do i = 1, ncol
             wr(i) = wr(i) + state%q(i,k,ixrain)*state%pdel(i,k)/gravit
             ws(i) = ws(i) + state%q(i,k,ixsnow)*state%pdel(i,k)/gravit
          end do
       end do
    end if

    ! Compute vertical integrals of frozen static energy and total water.
    do i = 1, ncol
!!     te(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
       te(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*( wl(i) + wr(i) )
       tw(i) = wv(i) + wl(i) + wi(i) + wr(i) + ws(i)
    end do

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

    deallocate(cpairv_loc)

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
!DIR$ CONCURRENT
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
!! Author: Kai Zhang 
!!...................................................................

    use camsrfexch,       only: cam_out_t, cam_in_t
    use physics_buffer,   only: physics_buffer_desc, pbuf_get_field, pbuf_get_chunk, pbuf_set_field 
    use cam_history,      only: outfld
    use phys_control,     only: ieflx_opt

    ! Compute global mean qflx

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

    real(r8), pointer :: ieflx(:) 

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

!DIR$ CONCURRENT
    do lchnk = begchunk, endchunk

       ncol = state(lchnk)%ncol
       qflx(:ncol,lchnk) = cam_in(lchnk)%cflx(:ncol,1)
       snow(:ncol,lchnk) = cam_out(lchnk)%precsc(:ncol) + cam_out(lchnk)%precsl(:ncol)
       rain(:ncol,lchnk) = cam_out(lchnk)%precc(:ncol)  + cam_out(lchnk)%precl(:ncol) - snow(:ncol,lchnk) 

       !! the calculation below (rhow*) converts the unit of precipitation from m/s to kg/m2/s 

       select case (ieflx_opt) 

       case(1) 
          ienet(:ncol,lchnk) = cpsw * qflx(:ncol,lchnk) * cam_in(lchnk)%ts(:ncol) - & 
                               cpsw * rhow * ( rain(:ncol,lchnk) + snow(:ncol,lchnk) ) * cam_out(lchnk)%tbot(:ncol)
       case(2) 
          ienet(:ncol,lchnk) = cpsw * qflx(:ncol,lchnk) * cam_in(lchnk)%ts(:ncol) - & 
                               cpsw * rhow * ( rain(:ncol,lchnk) + snow(:ncol,lchnk) ) * cam_in(lchnk)%ts(:ncol)
       case default 
          call endrun('*** incorrect ieflx_opt ***')
       end select 

       !! put it to pbuf for more comprehensive treatment in the future 

       call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk),ieflx_idx, ieflx)

       ieflx(:ncol) = ienet(:ncol,lchnk) 

       call outfld('IEFLX', ieflx(:ncol), pcols, lchnk)

    end do

    call gmean(ienet, ieflx_glob)

!!!    if (begchunk .le. endchunk) then
!!!       if (masterproc) then
!!!          write(iulog,'(1x,a12,1x,i8,4(1x,e25.17))') "nstep, ieflx, ieup, iedn ", nstep, ieflx_glob, ieup_glob, iedn_glob 
!!!       end if
!!!    end if

  end subroutine ieflx_gmean


!===============================================================================
  subroutine check_ieflx_fix(lchnk, ncol, nstep, shflx)

!!
!! Add implied internal energy flux to the sensible heat flux 
!!
!! Called by typhsac 
!! 

    use cam_history,       only: outfld

    integer, intent(in   ) :: nstep          ! time step number
    integer, intent(in   ) :: lchnk  
    integer, intent(in   ) :: ncol
    real(r8),intent(inout) :: shflx(pcols) 

    integer :: i

    if(nstep>1) then 
       do i = 1, ncol
          shflx(i) = shflx(i) + ieflx_glob 
       end do
    end if 

    call outfld('SHFLXFIX', shflx, pcols, lchnk)

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

!DIR$ CONCURRENT
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

    type(physics_state), intent(in   ) :: state
    type(physics_ptend), intent(out)   :: ptend

    integer , intent(in   ) :: nstep          ! time step number
    real(r8), intent(out  ) :: eshflx(pcols)  ! effective sensible heat flux

!---------------------------Local storage-------------------------------
    integer  :: i                        ! column
    integer  :: ncol                     ! number of atmospheric columns in chunk
!-----------------------------------------------------------------------
    ncol = state%ncol

    call physics_ptend_init(ptend, state%psetcols, 'chkenergyfix', ls=.true.)

#if ( defined OFFLINE_DYN )
    ! disable the energy fix for offline driver
    heat_glob = 0._r8
#endif
! add (-) global mean total energy difference as heating
    ptend%s(:ncol,:pver) = heat_glob
!!$    write(iulog,*) "chk_fix: heat", state%lchnk, ncol, heat_glob

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
    integer  i,k,m                                 ! column, level,constituent indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ and tracer indices
    integer :: ixrain, ixsnow                      ! RAINQM and SNOWQM indices

!-----------------------------------------------------------------------

    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain,   abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow,   abort=.false.)

    do m = 1,pcnst

       if ( any(m == (/ 1, ixcldliq, ixcldice, &
                           ixrain,   ixsnow    /)) ) exit   ! dont process water substances
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
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: ixrain, ixsnow                      ! RAINQM and SNOWQM indices
    integer :: m                            ! tracer index
    character(len=8) :: tracname   ! tracername
!-----------------------------------------------------------------------
!!$    if (.true.) return

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain,   abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow,   abort=.false.)

    do m = 1,pcnst

       if ( any(m == (/ 1, ixcldliq, ixcldice, &
                           ixrain,   ixsnow    /)) ) exit   ! dont process water substances
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

    real(r8) :: wv(state%ncol)                     ! vertical integral of water (vapor)
    real(r8) :: wl(state%ncol)                     ! vertical integral of water (liquid)
    real(r8) :: wi(state%ncol)                     ! vertical integral of water (ice)
    real(r8) :: tw(state%ncol)                     ! vertical integral of total water

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    real(r8) :: wr(state%ncol)                     ! vertical integral of rain
    real(r8) :: ws(state%ncol)                     ! vertical integral of snow
    integer :: ixrain
    integer :: ixsnow
!!...................................................................

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain, abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow, abort=.false.)


!! Compute vertical integrals of all water species (vapor, liquid, ice, rain, snow)
!!...................................................................
    wv = 0._r8
    wl = 0._r8
    wi = 0._r8
    wr = 0._r8
    ws = 0._r8

    do k = 1, pver
       do i = 1, ncol
          wv(i) = wv(i) + state%q(i,k,1       )*state%pdel(i,k)/gravit
       end do
    end do

    ! Don't require cloud liq/ice to be present.  Allows for adiabatic/ideal phys.
    if (ixcldliq > 1  .and.  ixcldice > 1) then
       do k = 1, pver
          do i = 1, ncol
             wl(i) = wl(i) + state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit
             wi(i) = wi(i) + state%q(i,k,ixcldice)*state%pdel(i,k)/gravit
          end do
       end do
    end if

    if (ixrain   > 1  .and.  ixsnow   > 1 ) then
       do k = 1, pver
          do i = 1, ncol
             wr(i) = wr(i) + state%q(i,k,ixrain)*state%pdel(i,k)/gravit
             ws(i) = ws(i) + state%q(i,k,ixsnow)*state%pdel(i,k)/gravit
          end do
       end do
    end if

!! Total water path
!!...................................................................
    do i = 1, ncol
       tw(i) = wv(i) + wl(i) + wi(i) + wr(i) + ws(i)
    end do

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


end module check_energy

