! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_surface_restoring
!
!> \brief MPAS ocean restoring
!> \author Todd Ringler
!> \date   06/08/2015
!> \details
!>  This module contains routines for computing the surface tracer flux due to restoring
!
!-----------------------------------------------------------------------

module ocn_tracer_surface_restoring

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_timekeeping
   use mpas_forcing
   use mpas_stream_manager
   use ocn_constants
   use ocn_config
   use ocn_framework_forcing
   use ocn_diagnostics_variables
   use ocn_mesh

   implicit none

   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_tracer_surface_restoring_compute,          &
             ocn_get_surfaceSalinityData,                   &
             ocn_salinity_restoring_forcing_write_restart,  &
             ocn_tracer_surface_restoring_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_surface_restoring_compute
!
!> \brief   computes a surface tracer flux due to surface restoring
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine computes a surface tracer flux due to surface restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_surface_restoring_compute(groupName, nTracers, nCells, tracers, pistonVelocity, &
                 tracersSurfaceRestoringValue, tracersSurfaceFlux, indexSalinity,  &
                 use_surface_salinity_monthly_restoring, salinity_restoring_constant_piston_velocity, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character (len=*), intent(in) :: groupName !< Input: Name of tracer group

      ! scalars
      integer, intent(in) ::  &
          nTracers,           &
          indexSalinity,      &
          nCells

      real (kind=RKIND), intent(in) :: salinity_restoring_constant_piston_velocity

      ! logicals
      logical, intent(in) ::  &
          use_surface_salinity_monthly_restoring

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers

      ! two dimensional ararys
      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         pistonVelocity,     &
         tracersSurfaceRestoringValue

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
        tracersSurfaceFlux

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, iLevel, iTracer

      err = 0

      !$omp parallel
      !$omp do schedule(runtime) private(iTracer, iLevel)
      do iCell=1,nCells
        iLevel = minLevelCell(iCell)  ! base surface flux restoring on tracer fields in the top layer

        do iTracer=1,nTracers

          ! For monthly salinity restoring, tracersSurfaceRestoringValue contains the zero-mean deltaS
          if (trim(groupName) == 'activeTracers' &
                .and. iTracer == indexSalinity   &
                .and. use_surface_salinity_monthly_restoring) then
          tracersSurfaceFlux(iTracer, iCell) =   tracersSurfaceFlux(iTracer, iCell) + &
                                                 salinity_restoring_constant_piston_velocity * &
                                                 tracersSurfaceRestoringValue(iTracer,iCell)
          else
          tracersSurfaceFlux(iTracer, iCell) =  tracersSurfaceFlux(iTracer, iCell) - &
                                                pistonVelocity(iTracer,iCell) *      &
                                                 (tracers(iTracer, iLevel, iCell) - tracersSurfaceRestoringValue(iTracer,iCell))
          endif
        enddo
      enddo
      !$omp end do
      !$omp end parallel

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_surface_restoring_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_surface_restoring_init
!
!> \brief   Initializes ocean surface restoring
!> \author  Todd Ringler
!> \date    06/09/2015
!> \details
!>  This routine initializes fields required for tracer surface flux restoring
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_surface_restoring_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_surface_restoring_init!}}}

!***********************************************************************

!***********************************************************************
!
!  routine get_surfaceSalinityData
!
!> \brief   retrieve data needed to compute surface salinity restoring using monthly climatology
!> \author  Mathew Maltrud
!> \date    09/23/16
!> \details
!>  This routine calls mpas_forcing routines to acquire needed surface salinity forcing data and interpolates
!>    between time levels.
!
!-----------------------------------------------------------------------

    subroutine ocn_get_surfaceSalinityData( streamManager, &
        domain, &
        simulationClock, &
        firstTimeStep) !{{{

        type (MPAS_streamManager_type), intent(inout) :: streamManager

        type (domain_type) :: domain
        type (MPAS_timeInterval_type) :: timeStepSurfaceSalinity
        type (MPAS_clock_type) :: simulationClock

        type (MPAS_Time_Type) :: currTime
        character(len=strKind) :: timeStamp

        logical, intent(in) :: firstTimeStep
        real(kind=RKIND) :: dt, sumAreaDeltaS, sumArea, avgDeltaS, deltaS, sumAreaDeltaSGlobal, sumAreaGlobal
        real(kind=RKIND) :: avgDeltaS1
        type (block_type), pointer :: block

        type (dm_info) :: dminfo

        type (mpas_pool_type), pointer :: forcingPool
        type (mpas_pool_type), pointer :: meshPool
        type (mpas_pool_type), pointer :: statePool
        type (mpas_pool_type), pointer :: tracersPool
        type (mpas_pool_type), pointer :: surfaceSalinityMonthlyForcing
        type (mpas_pool_type), pointer :: tracersSurfaceRestoringFieldsPool

        real (kind=RKIND), dimension(:), pointer :: &
           surfaceSalinityMonthlyClimatologyValue, iceFraction, areaCell

        real (kind=RKIND), dimension(:,:), pointer :: &
           activeTracersSurfaceRestoringValue

        real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

        integer, pointer :: nCells, nCellsSolve
        integer, pointer :: indexSalinityPtr, indexSalinitySurfaceRestoringValuePtr
        integer :: indexSalinity, indexSalinitySurfaceRestoringValue
        integer :: iCell, timeLevel
        integer, dimension(:), pointer :: landIceMask

        character(len=strKIND) :: &
           forcingIntervalMonthly,  &
           forcingReferenceTimeMonthly

        integer :: ierr

        integer, parameter :: nSums = 2
        real (kind=RKIND), dimension(nSums) :: reductions, sums

        ! initialize monthly forcing to be read from file

        if (firstTimeStep) then

           currTime = mpas_get_clock_time( simulationClock, MPAS_NOW, ierr)
           call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
           timeStamp = '0000'//trim(timeStamp(5:))

           forcingIntervalMonthly = "0000-01-00_00:00:00"
           forcingReferenceTimeMonthly = "0000-01-15_00:00:00"

           call MPAS_forcing_init_group( forcingGroupHead,  &
                "surfaceSalinityMonthlyClimatology", &
                domain, &
                timeStamp, &
                '0000-01-01_00:00:00', &
                '0001-00-00_00:00:00', &
                config_do_restart)

           call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'surfaceSalinityMonthlyClimatology', &  !  forcing group name
                'surfaceSalinityMonthlyClimatologyValue', &  !  array name
                'surface_salinity_monthly_data', &  !  stream name
                'surfaceSalinityMonthlyForcing',  & !  pool name
                'surfaceSalinityMonthlyClimatologyValue',  &   !  array name
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

           call MPAS_forcing_init_field_data( forcingGroupHead, &
                'surfaceSalinityMonthlyClimatology', &
                domain % streamManager, &
                config_do_restart, &
                .false.)

           return
        endif  !  first timestep

        call MPAS_set_timeInterval(timeStepSurfaceSalinity,timeString=config_surface_salinity_monthly_restoring_compute_interval)
        call MPAS_get_timeInterval(timeStepSurfaceSalinity,dt=dt)

        call MPAS_pool_get_subpool(domain % blocklist % structs, 'surfaceSalinityMonthlyForcing',  &
                                                                  surfaceSalinityMonthlyForcing)
        call MPAS_pool_get_array(surfaceSalinityMonthlyForcing, 'surfaceSalinityMonthlyClimatologyValue',   &
                                                                 surfaceSalinityMonthlyClimatologyValue)

        call MPAS_forcing_get_forcing(forcingGroupHead, &
             'surfaceSalinityMonthlyClimatology', streamManager, dt)

      sumAreaDeltaS = 0.0_RKIND
      sumArea = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields',tracersSurfaceRestoringFieldsPool)
        ! Use time level 1, which is always the new time level
        timeLevel = 1
        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)

        call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
        call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', &
                                                   indexSalinityPtr)
        call mpas_pool_get_dimension(tracersSurfaceRestoringFieldsPool, &
                                'index_salinitySurfaceRestoringValue',  &
                                 indexSalinitySurfaceRestoringValuePtr)
        indexSalinity = indexSalinityPtr
        indexSalinitySurfaceRestoringValue = &
                                 indexSalinitySurfaceRestoringValuePtr
        call mpas_pool_get_subpool(block % structs, 'surfaceSalinityMonthlyForcing',  &
                                                                  surfaceSalinityMonthlyForcing)
        call mpas_pool_get_array(surfaceSalinityMonthlyForcing, 'surfaceSalinityMonthlyClimatologyValue',   &
                                                                 surfaceSalinityMonthlyClimatologyValue)

        call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersSurfaceRestoringValue', &
           activeTracersSurfaceRestoringValue)

        call mpas_pool_get_array(forcingPool, 'iceFraction', iceFraction)
        call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)
        call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

        ! This is not in a threaded region, so no openMP pragmas are needed.
        if ( associated(landIceMask)) then

          if (config_salinity_restoring_under_sea_ice) then

            ! Simulation has landIceMask AND config_salinity_restoring_under_sea_ice=.true.
            do iCell = 1, nCells
              if (landIceMask(iCell)==1) then
                ! Turn off salinity restoring in this cell
                activeTracersSurfaceRestoringValue(indexSalinitySurfaceRestoringValue,iCell) = 0.0_RKIND
              else
                ! Turn on salinity restoring in this cell
                deltaS = surfaceSalinityMonthlyClimatologyValue(iCell) - activeTracers(indexSalinity,1,iCell)
                if (deltaS >  config_salinity_restoring_max_difference) deltaS =  config_salinity_restoring_max_difference
                if (deltaS < -config_salinity_restoring_max_difference) deltaS = -config_salinity_restoring_max_difference

                ! Salinity restoring below sea ice is always on, regardless of iceFraction value
                activeTracersSurfaceRestoringValue(indexSalinitySurfaceRestoringValue,iCell) = deltaS
              endif
            end do

          else ! config_salinity_restoring_under_sea_ice = .false.

            ! Simulation has landIceMask AND config_salinity_restoring_under_sea_ice=.false. (default)
            do iCell = 1, nCells
              if (landIceMask(iCell)==1) then
                ! Turn off salinity restoring in this cell
                activeTracersSurfaceRestoringValue(indexSalinitySurfaceRestoringValue,iCell) = 0.0_RKIND
              else
                ! Turn on salinity restoring in this cell
                deltaS = surfaceSalinityMonthlyClimatologyValue(iCell) - activeTracers(indexSalinity,1,iCell)
                if (deltaS >  config_salinity_restoring_max_difference) deltaS =  config_salinity_restoring_max_difference
                if (deltaS < -config_salinity_restoring_max_difference) deltaS = -config_salinity_restoring_max_difference

                ! Salinity restoring below sea ice tapers below partial sea ice
                ! coverage, from full in the open ocean to zero when iceFraction=1.0
                activeTracersSurfaceRestoringValue(indexSalinitySurfaceRestoringValue,iCell) = &
                   deltaS*(1.0_RKIND - iceFraction(iCell))
              endif
            end do

          endif

        else  ! associated(landIceMask)) = .false.

          if (config_salinity_restoring_under_sea_ice) then

            ! Simulation has NO landIceMask AND config_salinity_restoring_under_sea_ice=.true.
            do iCell = 1, nCells
              deltaS = surfaceSalinityMonthlyClimatologyValue(iCell) - activeTracers(indexSalinity,1,iCell)
              if (deltaS >  config_salinity_restoring_max_difference) deltaS =  config_salinity_restoring_max_difference
              if (deltaS < -config_salinity_restoring_max_difference) deltaS = -config_salinity_restoring_max_difference

              ! Salinity restoring below sea ice is always on, regardless of iceFraction value
              activeTracersSurfaceRestoringValue(indexSalinitySurfaceRestoringValue,iCell) = deltaS
            end do

          else ! config_salinity_restoring_under_sea_ice = .false.

            ! Simulation has NO landIceMask AND config_salinity_restoring_under_sea_ice=.false. (default)
            do iCell = 1, nCells
              deltaS = surfaceSalinityMonthlyClimatologyValue(iCell) - activeTracers(indexSalinity,1,iCell)
              if (deltaS >  config_salinity_restoring_max_difference) deltaS =  config_salinity_restoring_max_difference
              if (deltaS < -config_salinity_restoring_max_difference) deltaS = -config_salinity_restoring_max_difference

              ! Salinity restoring below sea ice tapers below partial sea ice
              ! coverage, from full in the open ocean to zero when iceFraction=1.0
              activeTracersSurfaceRestoringValue(indexSalinitySurfaceRestoringValue,iCell) = &
                 deltaS*(1.0_RKIND - iceFraction(iCell))
            end do

          endif

        endif

        do iCell=1,nCellsSolve
           deltaS = activeTracersSurfaceRestoringValue(indexSalinitySurfaceRestoringValue,iCell)
           if (deltaS .ne. 0.0_RKIND) then
              sumAreaDeltaS = sumAreaDeltaS + deltaS*areaCell(iCell)
              sumArea = sumArea + areaCell(iCell)
           endif
        enddo

        block => block % next
      end do

      ! Global sum to subtract global mean of deltaS
      dminfo = domain % dminfo
      sums(1) = sumAreaDeltaS
      sums(2) = sumArea
      call mpas_dmpar_sum_real_array(dminfo, nSums, sums(1:nSums), reductions(1:nSums))
      sumAreaDeltaSGlobal = reductions(1)
      sumAreaGlobal       = reductions(2)
      avgDeltaS1 = sumAreaDeltaSGlobal/(sumAreaGlobal + 1.e-20_RKIND)
      !LPV: There is no guarantee that a global sum of deltaS is BFB across
      !decompositions, the next line rounds each value to 10 decimal places to
      !ensure BFB results
      avgDeltaS = float (int(avgDeltaS1 * 1.0E10_RKIND + 0.5_RKIND)) / 1.0E10_RKIND

      block => domain % blocklist
      do while (associated(block))

        call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
        call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
        call mpas_pool_get_array(forcingPool, 'iceFraction', iceFraction)
        call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields',tracersSurfaceRestoringFieldsPool)
        call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersSurfaceRestoringValue', &
           activeTracersSurfaceRestoringValue)
        call mpas_pool_get_dimension(tracersSurfaceRestoringFieldsPool, &
                                'index_salinitySurfaceRestoringValue',  &
                                 indexSalinitySurfaceRestoringValuePtr)
        indexSalinitySurfaceRestoringValue = &
                                 indexSalinitySurfaceRestoringValuePtr

        do iCell = 1, nCells
          deltaS = activeTracersSurfaceRestoringValue(indexSalinitySurfaceRestoringValue,iCell)
          if (deltaS .ne. 0.0_RKIND) then
             activeTracersSurfaceRestoringValue(indexSalinitySurfaceRestoringValue,iCell) =  &
                deltaS - avgDeltaS
             salinitySurfaceRestoringTendency(iCell) = &
                activeTracersSurfaceRestoringValue(indexSalinitySurfaceRestoringValue,iCell) &
                * config_salinity_restoring_constant_piston_velocity
          else
             salinitySurfaceRestoringTendency(iCell) = 0.0_RKIND
          endif
        enddo

        block => block % next
      end do

    end subroutine ocn_get_surfaceSalinityData!}}}

!***********************************************************************
!
!  routine ocn_salinity_restoring_forcing_write_restart
!
!> \brief   writes restart timestamp for salinity restoring data to be read in on future restart
!> \author  Mathew Maltrud
!> \date    10/17/2016

!
!-----------------------------------------------------------------------

   subroutine ocn_salinity_restoring_forcing_write_restart(domain)!{{{

      type(domain_type) :: domain

      call MPAS_forcing_write_restart_times(forcingGroupHead)

    end subroutine ocn_salinity_restoring_forcing_write_restart!}}}

!***********************************************************************

end module ocn_tracer_surface_restoring

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
