!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_varying_forcing
!
!> \brief Module to handle time-varying forcing for the ocean core
!> \author Steven Brus, Xylar Asay-Davis, LANL
!> \date 2/3/2019 
!> \details
!>  Manages the wind and pressure forcing for the ocean.
!>  Based very closely on the seaice forcing module and
!>  framework written by Adrian Turner.
!
!-----------------------------------------------------------------------

module ocn_time_varying_forcing

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_forcing
  use mpas_stream_manager
  use ocn_framework_forcing
  use ocn_constants
  use mpas_log, only: mpas_log_write

  implicit none

  private
  public :: &
       ocn_time_varying_forcing_init, &
       ocn_time_varying_forcing_get, &
       ocn_time_varying_forcing_write_restart_times

contains

!-----------------------------------------------------------------------
! initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_varying_forcing_init
!
!> \brief Driver for ocean core forcing initialization
!> \author Steven Brus, Xylar Asay-Davis, LANL
!> \date 2/3/2019
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine ocn_time_varying_forcing_init(domain)!{{{

    type (domain_type) :: domain

    logical, pointer :: &
         config_use_time_varying_atmospheric_forcing, &
         config_use_time_varying_land_ice_forcing

    call MPAS_pool_get_config(domain % configs, "config_use_time_varying_atmospheric_forcing", config_use_time_varying_atmospheric_forcing)

    ! init the atmospheric forcing
    if (config_use_time_varying_atmospheric_forcing) then
       call init_atmospheric_forcing(domain)
    endif

    call MPAS_pool_get_config(domain % configs, "config_use_time_varying_land_ice_forcing", config_use_time_varying_land_ice_forcing)

    ! init other forcings, to be added later as needed
    if (config_use_time_varying_land_ice_forcing) then
      call init_land_ice_forcing(domain)
    endif

    ! init other forcings, to be added later as needed
    !if (config_use_time_varying_****_forcing) then
    !  call init_****_forcing(domain)
    !endif

  end subroutine ocn_time_varying_forcing_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_atmospheric_forcing
!
!> \brief
!> \author Steven Brus, LANL
!> \date 1/8/2019
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_atmospheric_forcing(domain)!{{{

    type (domain_type) :: domain

    character(len=strKIND), pointer :: &
         config_atmospheric_forcing_type

    call MPAS_pool_get_config(domain % configs, "config_time_varying_atmospheric_forcing_type", config_atmospheric_forcing_type)

    select case (trim(config_atmospheric_forcing_type))
      case ("WINDPRES")
         call init_atmospheric_forcing_WINDPRES(domain)
      case default
         call mpas_log_write("Atmospheric forcing type unknown: "//trim(config_atmospheric_forcing_type), MPAS_LOG_CRIT)
    end select

  end subroutine init_atmospheric_forcing!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_atmospheric_forcing_WINDPRES
!
!> \brief Initialize the forcing objects for wind and pressure
!> \author Steven Brus, LANL
!> \date 1/8/2019
!> \details
!>  This routine calls the MPAS_forcing module subroutines that initializes
!>  the forcings type
!
!-----------------------------------------------------------------------

  subroutine init_atmospheric_forcing_WINDPRES(domain)!{{{

    type(domain_type) :: domain

    character(len=strKIND), pointer :: &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration, &
         forcingInterval, &
         forcingReferenceTime

    logical, pointer :: &
         config_do_restart

    ! get atmospheric forcing configuration options
    call MPAS_pool_get_config(&
         domain % configs, &
         "config_time_varying_atmospheric_forcing_start_time", &
         config_forcing_start_time)

    call MPAS_pool_get_config(&
         domain % configs, &
         "config_time_varying_atmospheric_forcing_cycle_start", &
         config_forcing_cycle_start)

    call MPAS_pool_get_config(&
         domain % configs, &
         "config_time_varying_atmospheric_forcing_cycle_duration", &
         config_forcing_cycle_duration)

    call MPAS_pool_get_config(&
         domain % configs, &
         "config_time_varying_atmospheric_forcing_interval", &
         forcingInterval)

    call MPAS_pool_get_config(&
         domain % configs, &
         "config_time_varying_atmospheric_forcing_reference_time", &
         forcingReferenceTime)

    call MPAS_pool_get_config(&
         domain % configs, &
         "config_do_restart", &
         config_do_restart)

    ! create the atmospheric forcing group
    call MPAS_forcing_init_group(&
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_atmospheric_forcing", &            ! forcingGroupName
         domain, &                               ! domain
         config_forcing_start_time, &            ! startTimeStr
         config_forcing_cycle_start, &           ! forcingCycleStart
         config_forcing_cycle_duration, &        ! forcingCycleDuration
         config_do_restart)!, &                  ! restart
         !.false.)                               ! forcingCycleStartInclusive

    ! windSpeedU init
    call MPAS_forcing_init_field(&
         domain % streamManager, &               ! streamManager
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_atmospheric_forcing", &            ! forcingGroupName
         "windSpeedU", &                         ! forcingName
         "atmospheric_forcing", &                ! forcingStreamID
         "timeVaryingForcing", &                 ! poolname
         "windSpeedU", &                         ! fieldname
         "linear", &                             ! interpolationType
         forcingReferenceTime, &                 ! forcingReferenceTimeStr
         forcingInterval)!,                      ! forcingIntervalStr
         !"next")                                ! forcingInitializationType

    ! windSpeedV init
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         forcingGroupHead, &
         "ocn_atmospheric_forcing", &
         "windSpeedV", &
         "atmospheric_forcing", &
         "timeVaryingForcing", &
         "windSpeedV", &
         "linear", &
         forcingReferenceTime, &
         forcingInterval)!, &
         !"next")

    ! atmosphericPressure init
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         forcingGroupHead, &
         "ocn_atmospheric_forcing", &
         "atmosPressure", &
         "atmospheric_forcing", &
         "timeVaryingForcing", &
         "atmosPressure", &
         "linear", &
         forcingReferenceTime, &
         forcingInterval)!, &
         !"next")

    call MPAS_forcing_init_field_data(&
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_atmospheric_forcing", &            ! forcingGroupName
         domain % streamManager, &               ! streamManager
         config_do_restart, &                    ! restart
         .false.)                                ! interpolateAtInit

  end subroutine init_atmospheric_forcing_WINDPRES!}}}

!-----------------------------------------------------------------------
! runtime
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_varying_forcing_get
!
!> \brief Retrieve forcing data during time stepping
!> \author Steven Brus, Xylar Asay-Davis, LANL
!> \date 2/3/2019
!> \details
!>  This routine calls the MPAS_forcing routine that will perform the
!>  forcing data aquisition and interpolation during timestepping
!
!-----------------------------------------------------------------------

  subroutine ocn_time_varying_forcing_get(streamManager, domain, simulationClock)!{{{

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    type (MPAS_clock_type) :: simulationClock

    logical, pointer :: &
         config_use_time_varying_atmospheric_forcing, &
         config_use_time_varying_land_ice_forcing

    call MPAS_pool_get_config(domain % configs, "config_use_time_varying_atmospheric_forcing", config_use_time_varying_atmospheric_forcing)

    if (config_use_time_varying_atmospheric_forcing) then

        call atmospheric_forcing(&
             streamManager, &
             domain, &
             simulationClock)

    endif

    call MPAS_pool_get_config(domain % configs, "config_use_time_varying_land_ice_forcing", &
                              config_use_time_varying_land_ice_forcing)

    if (config_use_time_varying_land_ice_forcing) then

        call land_ice_forcing(&
             streamManager, &
             domain, &
             simulationClock)

    endif

    ! call other forcings, to be added if needed
    ! if (config_use_****_forcing) then
    !   call ****_forcing(&
    !        streamManager, &
    !        domain, &
    !        simulationClock)
    ! endif


  end subroutine ocn_time_varying_forcing_get!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  atmospheric_forcing
!
!> \brief
!> \author Steven Brus, LANL
!> \date 1/8/2019
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine atmospheric_forcing(streamManager, domain, simulationClock)!{{{

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    type (MPAS_clock_type) :: simulationClock

    type (MPAS_time_type) :: currentForcingTime

    type (block_type), pointer :: block

    character(len=StrKIND), pointer :: config_dt

    character(len=strKIND), pointer :: config_atmospheric_forcing_type

    type (MPAS_timeInterval_type) :: timeStepESMF

    character(len=StrKIND) :: timeStamp

    real (kind=RKIND) :: dtSim

    integer :: err

    ! configurations
    call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)
    call mpas_pool_get_config(domain % configs, 'config_time_varying_atmospheric_forcing_type', config_atmospheric_forcing_type)

    ! convert config_dt to real
    call mpas_set_timeInterval(timeStepESMF, timeString=config_dt,ierr=err)
    call mpas_get_timeInterval(timeStepESMF, dt=dtSim)

    ! use the forcing layer to get data
    call MPAS_forcing_get_forcing(&
         forcingGroupHead, &               ! forcingGroupHead
         "ocn_atmospheric_forcing", &      ! forcingGroupName
         streamManager, &                  ! streamManager
         dtSim)                            ! dt

    call MPAS_forcing_get_forcing_time(&
         forcingGroupHead, &               ! forcingGroupHead
         "ocn_atmospheric_forcing", &      ! forcingGroupName
         currentForcingTime)               ! forcingTime
    !call mpas_get_time(curr_time=currentForcingTime, dateTimeString=timeStamp, ierr=err)
    !call mpas_log_write('Forcing time ' // trim(timeStamp))

    ! perform post forcing
    block => domain % blocklist
    do while (associated(block))

       ! apply drag law to compute wind stress from wind speed
       if (trim(config_atmospheric_forcing_type) == "WINDPRES") then
          call post_atmospheric_forcing(block)
       endif

       block => block % next
    end do

  end subroutine atmospheric_forcing!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  post_atmospheric_forcing
!
!> \brief Calculate wind stress from wind speed and apply ramp
!> \author Steven Brus, LANL
!> \date 1/8/2019
!> \details
!>  Uses drag law from Garratt, J.R., 1977. Review of drag coefficients
!>  over oceans and continents. Monthly weather review, 105(7), pp.915-929.
!>
!
!-----------------------------------------------------------------------

  subroutine post_atmospheric_forcing(block)!{{{

    type (block_type), pointer :: block

    type (MPAS_time_type) :: currentForcingTime

    type (mpas_pool_type), pointer :: &
         mesh, &
         forcingPool, &
         timeVaryingForcingPool, &
         diagnosticsPool

    real(kind=RKIND), dimension(:), pointer :: &
         windSpeedU, &
         windSpeedV, &
         windSpeedMagnitude, &
         windStressZonal, &
         windStressMeridional, &
         atmosPressure, &
         atmosphericPressure

    real(kind=RKIND) :: &
         windStressCoefficient, &
         rhoAir, &
         ramp, &
         currentTime, &
         windStressCoefficientLimit

    real(kind=RKIND), pointer:: &
        config_forcing_ramp, &
        daysSinceStartOfSim

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    call MPAS_pool_get_config(block % configs, "config_time_varying_atmospheric_forcing_ramp", config_forcing_ramp)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "forcing", forcingPool)
    call MPAS_pool_get_subpool(block % structs, "timeVaryingForcing", timeVaryingForcingPool)
    call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnosticsPool)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(diagnosticsPool, "daysSinceStartOfSim", daysSinceStartOfSim)
    call MPAS_pool_get_array(timeVaryingForcingPool, "windSpeedU", windSpeedU)
    call MPAS_pool_get_array(timeVaryingForcingPool, "windSpeedV", windSpeedV)
    call MPAS_pool_get_array(timeVaryingForcingPool, "windSpeedMagnitude", windSpeedMagnitude)
    call MPAS_pool_get_array(forcingPool, "windStressZonal", windStressZonal)
    call MPAS_pool_get_array(forcingPool, "windStressMeridional", windStressMeridional)

    call MPAS_pool_get_array(timeVaryingForcingPool, "atmosPressure", atmosPressure)
    call MPAS_pool_get_array(forcingPool, "atmosphericPressure", atmosphericPressure)

    rhoAir = 1.225_RKIND
    windStressCoefficientLimit = 0.0035_RKIND

    ramp = tanh((2.0_RKIND*daysSinceStartOfSim)/config_forcing_ramp)

    !$omp do schedule(runtime) private(iCell, windStressCoefficient)
    do iCell = 1, nCells

       windSpeedU(iCell) = ramp*windSpeedU(iCell)
       windSpeedV(iCell) = ramp*windSpeedV(iCell)
       windSpeedMagnitude(iCell) = sqrt(windSpeedU(iCell)**2 + windSpeedV(iCell)**2)
       windStressCoefficient = (0.75_RKIND + 0.067_RKIND * windSpeedMagnitude(iCell)) * 0.001_RKIND   ! Garratt 1977
       if (windStressCoefficient > windStressCoefficientLimit) then
         windStressCoefficient = windStressCoefficientLimit
       end if

       windStressZonal(iCell) =  rhoAir * windSpeedU(iCell) * windSpeedMagnitude(iCell) * windStressCoefficient
       windStressMeridional(iCell) = rhoAir * windSpeedV(iCell) * windSpeedMagnitude(iCell) * windStressCoefficient

       atmosphericPressure(iCell) = ((1.0_RKIND-ramp)*atm_ref_pressure + ramp*atmosPressure(iCell)) - atm_ref_pressure

    enddo
    !$omp end do

  end subroutine post_atmospheric_forcing!}}}


!-----------------------------------------------------------------------
! restart
!-----------------------------------------------------------------------
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_varying_forcing_write_restart_times
!
!> \brief
!> \author Steven Brus, LANL
!> \date 1/8/2019
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine ocn_time_varying_forcing_write_restart_times(domain)!{{{

    type(domain_type) :: domain

    call MPAS_forcing_write_restart_times(forcingGroupHead)

  end subroutine ocn_time_varying_forcing_write_restart_times!}}}

!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_land_ice_forcing
!
!> \brief
!> \author Xylar Asay-Davis, LANL
!> \date 2/3/2019
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_land_ice_forcing(domain)!{{{

    type (domain_type) :: domain

    character(len=strKIND), pointer :: &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration, &
         forcingInterval, &
         forcingReferenceTime

    logical, pointer :: &
         config_do_restart


    ! get land-ice forcing configuration options
    call mpas_pool_get_config(domain % configs, "config_time_varying_land_ice_forcing_start_time", &
         config_forcing_start_time)

    call mpas_pool_get_config(domain % configs, "config_time_varying_land_ice_forcing_cycle_start", &
         config_forcing_cycle_start)

    call mpas_pool_get_config(domain % configs, "config_time_varying_land_ice_forcing_cycle_duration", &
         config_forcing_cycle_duration)

    call mpas_pool_get_config(domain % configs, "config_time_varying_land_ice_forcing_interval", &
         forcingInterval)

    call mpas_pool_get_config(domain % configs, "config_time_varying_land_ice_forcing_reference_time", &
         forcingReferenceTime)

    call mpas_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    ! create the land-ice forcing group
    call mpas_forcing_init_group(forcingGroupHead, "ocn_land_ice_forcing", domain, config_forcing_start_time, &
         config_forcing_cycle_start, config_forcing_cycle_duration, config_do_restart)

    call mpas_forcing_init_field(domain % streamManager, forcingGroupHead, "ocn_land_ice_forcing", &
         "landIcePressure", "land_ice_forcing", "timeVaryingForcing", "landIcePressureForcing", "linear", &
         forcingReferenceTime, forcingInterval)

    call mpas_forcing_init_field(domain % streamManager, forcingGroupHead, "ocn_land_ice_forcing", &
         "landIceFraction", "land_ice_forcing", "timeVaryingForcing", "landIceFractionForcing", "linear", &
         forcingReferenceTime, forcingInterval)

    call mpas_forcing_init_field(domain % streamManager, forcingGroupHead, "ocn_land_ice_forcing", &
         "landIceDraft", "land_ice_forcing", "timeVaryingForcing", "landIceDraftForcing", "linear", &
         forcingReferenceTime, forcingInterval)

    call mpas_forcing_init_field_data(&
         forcingGroupHead, &                     ! forcingGroupHead
         "ocn_land_ice_forcing", &               ! forcingGroupName
         domain % streamManager, &               ! streamManager
         config_do_restart, &                    ! restart
         .true.)                                 ! interpolateAtInit

  end subroutine init_land_ice_forcing!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  land_ice_forcing
!
!> \brief Compute land-ice pressure
!> \author Xylar Asay-Davis, LANL
!> \date 2/3/2019
!> \details Interpolate the pressure from overlying land ice
!>
!
!-----------------------------------------------------------------------

  subroutine land_ice_forcing(streamManager, domain, simulationClock)!{{{

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    type (MPAS_clock_type) :: simulationClock

    type (MPAS_time_type) :: currentForcingTime

    type (block_type), pointer :: block

    character(len=StrKIND), pointer :: config_dt

    type (MPAS_timeInterval_type) :: timeStepESMF

    character(len=StrKIND) :: timeStamp

    real (kind=RKIND) :: dtSim

    integer :: err

    type (mpas_pool_type), pointer :: &
         mesh, &
         forcingPool, &
         timeVaryingForcingPool

    real(kind=RKIND), dimension(:), pointer :: &
         landIceFractionForcing, &
         landIceFraction, &
         landIcePressureForcing, &
         landIcePressure, &
         landIceDraftForcing, &
         landIceDraft

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    ! configurations
    call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)

    ! convert config_dt to real
    call mpas_set_timeInterval(timeStepESMF, timeString=config_dt,ierr=err)
    call mpas_get_timeInterval(timeStepESMF, dt=dtSim)

    ! use the forcing layer to get data
    call mpas_forcing_get_forcing(forcingGroupHead, "ocn_land_ice_forcing", streamManager, dtSim)

    call mpas_forcing_get_forcing_time(forcingGroupHead, "ocn_land_ice_forcing", currentForcingTime)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "forcing", forcingPool)
       call MPAS_pool_get_subpool(block % structs, "timeVaryingForcing", timeVaryingForcingPool)
   
       call MPAS_pool_get_dimension(mesh, "nCells", nCells)
   
       call MPAS_pool_get_array(timeVaryingForcingPool, "landIceFractionForcing", landIceFractionForcing)
       call MPAS_pool_get_array(timeVaryingForcingPool, "landIcePressureForcing", landIcePressureForcing)
       call MPAS_pool_get_array(timeVaryingForcingPool, "landIceDraftForcing", landIceDraftForcing)
       call MPAS_pool_get_array(forcingPool, "landIceFraction", landIceFraction)
       call MPAS_pool_get_array(forcingPool, "landIcePressure", landIcePressure)
       call MPAS_pool_get_array(forcingPool, "landIceDraft", landIceDraft)
   
       do iCell = 1, nCells
   
          landIceFraction(iCell) = landIceFractionForcing(iCell)
          landIcePressure(iCell) = landIcePressureForcing(iCell)
          landIceDraft(iCell) = landIceDraftForcing(iCell)
   
       enddo
       block => block % next
    end do

  end subroutine land_ice_forcing!}}}

end module ocn_time_varying_forcing
! vim: foldmethod=marker
