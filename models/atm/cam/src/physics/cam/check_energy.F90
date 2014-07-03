
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
!   03.03.29  Boville  Add global energy check and fixer.        
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver, pverp, begchunk, endchunk
  use spmd_utils,      only: masterproc
  
  use phys_gmean,      only: gmean
  use physconst,       only: gravit, latvap, latice
  use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_ptend_init
  use constituents,    only: cnst_get_ind, pcnst, cnst_name, cnst_get_type_byind
  use time_manager,    only: is_first_step
  use cam_logfile,     only: iulog

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


! Private module data

  logical  :: print_energy_errors = .false.

  real(r8) :: teout_glob           ! global mean energy of output state
  real(r8) :: teinp_glob           ! global mean energy of input state
  real(r8) :: tedif_glob           ! global mean energy difference
  real(r8) :: psurf_glob           ! global mean surface pressure
  real(r8) :: ptopb_glob           ! global mean top boundary pressure
  real(r8) :: heat_glob            ! global mean heating rate

! Physics buffer indices
  
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
    
    use physics_buffer, only : pbuf_add_field, dtype_r8, pbuf_times

!-----------------------------------------------------------------------

! Request physics buffer space for fields that persist across timesteps.

    call pbuf_add_field('TEOUT', 'global',dtype_r8 , (/pcols,pbuf_times/),     teout_idx)
    call pbuf_add_field('DTCORE','global',dtype_r8,  (/pcols,pver,pbuf_times/),dtcore_idx )

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
    use cam_history,       only: addfld, phys_decomp

    implicit none
!-----------------------------------------------------------------------

! register history variables
    call addfld('TEINP   ', 'W/m2', 1,    'A', 'Total energy of physics input',    phys_decomp)
    call addfld('TEOUT   ', 'W/m2', 1,    'A', 'Total energy of physics output',   phys_decomp)
    call addfld('TEFIX   ', 'W/m2', 1,    'A', 'Total energy after fixer',         phys_decomp)
    call addfld('DTCORE'  , 'K/s' , pver, 'I', 'T tendency due to dynamical core', phys_decomp)

    if (masterproc) then
       write (iulog,*) ' print_energy_errors is set', print_energy_errors
    endif

  end subroutine check_energy_init

!===============================================================================

  subroutine check_energy_timestep_init(state, tend, pbuf)
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field, pbuf_times
!-----------------------------------------------------------------------
! Compute initial values of energy and water integrals, 
! zero cumulative tendencies
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state),   intent(inout)    :: state
    type(physics_tend ),   intent(inout)    :: tend
    type(physics_buffer_desc), pointer      :: pbuf(:)
!---------------------------Local storage-------------------------------

    real(r8) :: ke(pcols)                          ! vertical integral of kinetic energy
    real(r8) :: se(pcols)                          ! vertical integral of static energy
    real(r8) :: wv(pcols)                          ! vertical integral of water (vapor)
    real(r8) :: wl(pcols)                          ! vertical integral of water (liquid)
    real(r8) :: wi(pcols)                          ! vertical integral of water (ice)
    real(r8) :: tr(pcols)                          ! vertical integral of tracer
    real(r8) :: trpdel(pcols, pver)                ! pdel for tracer

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: itim                                ! pbuf time index
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)

! Compute vertical integrals of dry static energy and water (vapor, liquid, ice)
    ke = 0._r8
    se = 0._r8
    wv = 0._r8
    wl = 0._r8
    wi = 0._r8
    do k = 1, pver
       do i = 1, ncol
          ke(i) = ke(i) + 0.5_r8*(state%u(i,k)**2 + state%v(i,k)**2)*state%pdel(i,k)/gravit
          se(i) = se(i) + state%s(i,k         )*state%pdel(i,k)/gravit
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

! Compute vertical integrals of frozen static energy and total water.
    do i = 1, ncol
       state%te_ini(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
       state%tw_ini(i) = wv(i) + wl(i) + wi(i)

       state%te_cur(i) = state%te_ini(i)
       state%tw_cur(i) = state%tw_ini(i)
    end do

! zero cummulative boundary fluxes 
    tend%te_tnd(:ncol) = 0._r8
    tend%tw_tnd(:ncol) = 0._r8

    state%count = 0

! initialize physics buffer
    if (is_first_step()) then
       call pbuf_set_field(pbuf, teout_idx, state%te_ini)
    end if

  end subroutine check_energy_timestep_init

!===============================================================================

  subroutine check_energy_chng(state, tend, name, nstep, ztodt,        &
       flx_vap, flx_cnd, flx_ice, flx_sen)

!-----------------------------------------------------------------------
! Check that the energy and water change matches the boundary fluxes
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state)    , intent(inout) :: state
    type(physics_tend )    , intent(inout) :: tend
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: flx_vap(pcols)      ! boundary flux of vapor         (kg/m2/s)
    real(r8), intent(in   ) :: flx_cnd(pcols)      ! boundary flux of liquid+ice    (m/s) (precip?)
    real(r8), intent(in   ) :: flx_ice(pcols)      ! boundary flux of ice           (m/s) (snow?)
    real(r8), intent(in   ) :: flx_sen(pcols)      ! boundary flux of sensible heat (w/m2)

!******************** BAB ******************************************************
!******* Note that the precip and ice fluxes are in precip units (m/s). ********
!******* I would prefer to have kg/m2/s.                                ********
!******* I would also prefer liquid (not total) and ice fluxes          ********
!*******************************************************************************

!---------------------------Local storage-------------------------------

    real(r8) :: te_xpd(pcols)                   ! expected value (f0 + dt*boundary_flux)
    real(r8) :: te_dif(pcols)                   ! energy of input state - original energy
    real(r8) :: te_tnd(pcols)                   ! tendency from last process
    real(r8) :: te_rer(pcols)                   ! relative error in energy column

    real(r8) :: tw_xpd(pcols)                   ! expected value (w0 + dt*boundary_flux)
    real(r8) :: tw_dif(pcols)                   ! tw_inp - original water
    real(r8) :: tw_tnd(pcols)                   ! tendency from last process
    real(r8) :: tw_rer(pcols)                   ! relative error in water column

    real(r8) :: ke(pcols)                          ! vertical integral of kinetic energy
    real(r8) :: se(pcols)                          ! vertical integral of static energy
    real(r8) :: wv(pcols)                          ! vertical integral of water (vapor)
    real(r8) :: wl(pcols)                          ! vertical integral of water (liquid)
    real(r8) :: wi(pcols)                          ! vertical integral of water (ice)

    real(r8) :: te(pcols)                          ! vertical integral of total energy
    real(r8) :: tw(pcols)                          ! vertical integral of total water

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    logical :: flag
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)

    ! Compute vertical integrals of dry static energy and water (vapor, liquid, ice)
    ke = 0._r8
    se = 0._r8
    wv = 0._r8
    wl = 0._r8
    wi = 0._r8
    do k = 1, pver
       do i = 1, ncol
          ke(i) = ke(i) + 0.5_r8*(state%u(i,k)**2 + state%v(i,k)**2)*state%pdel(i,k)/gravit
          se(i) = se(i) + state%s(i,k         )*state%pdel(i,k)/gravit
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

    ! Compute vertical integrals of frozen static energy and total water.
    do i = 1, ncol
       te(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
       tw(i) = wv(i) + wl(i) + wi(i)
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

       ! relative error, expected value - input state / previous state 
       te_rer(i) = (te_xpd(i) - te(i)) / state%te_cur(i)
    end do

    ! relative error for total water (allow for dry atmosphere)
    tw_rer = 0._r8
    where (state%tw_cur(:ncol) > 0._r8) 
       tw_rer(:ncol) = (tw_xpd(:ncol) - tw(:ncol)) / state%tw_cur(:ncol)
    end where

    ! error checking
    if (print_energy_errors) then
       if (any(abs(te_rer(1:ncol)) > 1.E-14_r8 .or. abs(tw_rer(1:ncol)) > 1.E-10_r8)) then
          do i = 1, ncol
             ! the relative error threshold for the water budget has been reduced to 1.e-10
             ! to avoid messages generated by QNEG3 calls
             ! PJR- change to identify if error in energy or water 
             if (abs(te_rer(i)) > 1.E-14_r8 ) then 
                state%count = state%count + 1
                write(iulog,*) "significant energy conservation error after ", name,        &
                      " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
                write(iulog,*) te(i),te_xpd(i),te_dif(i),tend%te_tnd(i)*ztodt,  &
                      te_tnd(i)*ztodt,te_rer(i)
             endif
             if ( abs(tw_rer(i)) > 1.E-10_r8) then
                state%count = state%count + 1
                write(iulog,*) "significant water conservation error after ", name,        &
                      " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
                write(iulog,*) tw(i),tw_xpd(i),tw_dif(i),tend%tw_tnd(i)*ztodt,  &
                      tw_tnd(i)*ztodt,tw_rer(i)
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

    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk, pbuf_old_tim_idx
    
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
    integer :: itim                      ! time  index in pbuf
    integer :: lchnk                     ! chunk index

    real(r8) :: te(pcols,begchunk:endchunk,3)   
                                         ! total energy of input/output states (copy)
    real(r8) :: te_glob(3)               ! global means of total energy
    real(r8), pointer :: teout(:)
!-----------------------------------------------------------------------

    ! Find previous total energy in buffer
    itim = pbuf_old_tim_idx()

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
    integer  :: i,k                      ! column, level indexes
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

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k,m                                 ! column, level,constituent indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ and tracer indices

!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)

    do m = 1,pcnst

       if (m.eq.1.or.m.eq.ixcldice.or.m.eq.ixcldliq) exit   ! dont process water substances
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

    use abortutils, only: endrun 


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
    integer :: m                            ! tracer index
    character(len=8) :: tracname   ! tracername
!-----------------------------------------------------------------------
!!$    if (.true.) return

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)

    do m = 1,pcnst
       
       if (m.eq.1.or.m.eq.ixcldice.or.m.eq.ixcldliq) exit   ! dont process water substances
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


end module check_energy
