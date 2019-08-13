!=======================================================================
!
!BOP
!
! !MODULE: CICE_InitMod - performs CICE initialization
!
! !DESCRIPTION:
!
!  This module contains the CICE initialization routine that sets model
!  parameters and initializes the grid and CICE state variables.
!
! !REVISION HISTORY:
!  SVN:$Id: CICE_InitMod.F90 52 2007-01-30 18:04:24Z eclare $
!
!  authors Elizabeth C. Hunke, LANL
!          William H. Lipscomb, LANL
!          Philip W. Jones, LANL
!
! 2006: Converted to free form source (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module CICE_InitMod
!
! !USES:
!
      use ice_aerosol
      use ice_age
      use ice_FY
      use ice_calendar
      use ice_communicate
      use ice_diagnostics
      use ice_domain
      use ice_dyn_evp
      use ice_exit
      use ice_fileunits
      use ice_flux
      use ice_grid
      use ice_history
      use ice_restart
      use ice_init
      use ice_itd
      use ice_kinds_mod
      use ice_lvl
      use ice_mechred
      use ice_meltpond
      use ice_shortwave
      use ice_state, only: tr_aero, tr_iage, tr_FY, tr_pond, tr_lvl
      use ice_therm_itd
      use ice_therm_vertical
      use ice_timers
      use ice_transport_driver
      use ice_transport_remap
      use ice_work

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: CICE_Init

!
!EOP
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: CICE_Init - initialize CICE model
!
! !DESCRIPTION:
!
!  Initialize CICE model.
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine CICE_Init( mpicom_ice )
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         mpicom_ice ! communicator for sequential ccsm
!EOP
!
!     local temporary variables

      character(len=char_len_long) :: fname
      integer  (kind=int_kind)     :: ltmp

      call init_communicate( mpicom_ice ) ! initial setup for message passing
      call init_fileunits       ! set unit numbers (including nu_diag)
      call input_data           ! namelist variables
      call init_work            ! work arrays

      call init_domain_blocks   ! set up block decomposition
      call init_grid1           ! domain distribution
      call init_ice_timers      ! initialize all timers
      call ice_timer_start(timer_total)   ! start timing entire run
      call init_grid2           ! grid variables

      call init_calendar        ! initialize some calendar stuff
      call init_hist (dt)       ! initialize output history file
      call init_evp (dt_dyn)    ! define evp dynamics parameters, variables
      call init_coupler_flux    ! initialize fluxes exchanged with coupler
      call init_thermo_vertical ! initialize vertical thermodynamics
      call init_itd             ! initialize ice thickness distribution
      call init_state           ! initialize the ice state
      call init_transport       ! initialize horizontal transport

      if (runtype == 'continue') then ! start from core restart file
         ! given by pointer in ice_in
         if (restart_format == 'bin') then
            call restartfile()        
         else
            call restartfile(usepio=.true.)
         end if
         call calendar(time)          ! For continuation runs.
      else if (restart) then          
         ! ice_ic = core restart file or 'default' or 'none'
         ltmp = len_trim(ice_ic)
         if (ice_ic(ltmp-2:ltmp) /= '.nc') then
            call restartfile (ice_ic)      
         else
            call restartfile (usepio=.true., ice_ic=ice_ic)
         end if
      endif

      ! tracers
      if (tr_iage) call init_age        ! ice age tracer
      if (tr_FY)   call init_FY         ! FY area tracer
      if (tr_lvl)  call init_lvl        ! level ice tracer
      if (tr_pond) call init_meltponds  ! melt ponds
      if (tr_aero) call init_aerosol    ! ice aerosol MH

      call init_diags           ! initialize diagnostic output points
      call init_history_therm   ! initialize thermo history variables
      call init_history_dyn     ! initialize dynamic history variables

      ! Initialize shortwave components using swdn from previous timestep
      ! if restarting. These components will be scaled to current forcing
      ! in prep_radiation.
      if (runtype == 'continue' .or. restart) &
         call init_shortwave    ! initialize radiative transfer

#if (defined _NOIO)
! Not enought memory on BGL to write a history file yet!
!      if(write_ic .and. .not.prescribed_ice) call ice_write_hist(dt)
#else
      if(write_ic .and. .not.prescribed_ice) call ice_write_hist(dt)
#endif
      write_ic = .false.

      end subroutine CICE_Init

!=======================================================================

      end module CICE_InitMod

!=======================================================================
