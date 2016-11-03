! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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

   implicit none

   type(MPAS_forcing_group_type),pointer :: surfaceSalinityForcingGroup

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

   subroutine ocn_tracer_surface_restoring_compute(groupName, nTracers, nCellsSolve, tracers, pistonVelocity, &
                 tracersSurfaceRestoringValue, tracersSurfaceFlux, indexSalinity,  &
                 use_surface_salinity_monthly_restoring, err)!{{{

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
          nCellsSolve

      ! logicals
      logical, intent(in) ::  &
          use_surface_salinity_monthly_restoring

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers

      ! two dimensional ararys
      real (kind=RKIND), dimension(:,:), intent(in) :: &
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

      iLevel = 1  ! base surface flux restoring on tracer fields in the top layer

      !$omp do schedule(runtime) private(iTracer)
      do iCell=1,nCellsSolve
        do iTracer=1,nTracers

!maltrud for monthly salinity restoring, tracersSurfaceRestoringValue contains the zero-mean deltaS
          if (trim(groupName) == 'activeTracers' &
                .and. iTracer == indexSalinity   &
                .and. use_surface_salinity_monthly_restoring) then
!maltrud bug?
!         tracersSurfaceFlux(iTracer, iCell) =   tracersSurfaceFlux(iTracer, iCell) - &
          tracersSurfaceFlux(iTracer, iCell) =   tracersSurfaceFlux(iTracer, iCell) + &
                                                 pistonVelocity(iTracer,iCell) *      &
                                                 tracersSurfaceRestoringValue(iTracer,iCell)
          else
          tracersSurfaceFlux(iTracer, iCell) =   tracersSurfaceFlux(iTracer, iCell) - &
                                                 pistonVelocity(iTracer,iCell) *      &
                                                 (tracers(iTracer, iLevel, iCell) - tracersSurfaceRestoringValue(iTracer,iCell))
          endif
        enddo
      enddo
      !$omp end do

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

        logical, intent(in) :: firstTimeStep
        character(len=strKind), pointer :: config_dt
        real(kind=RKIND) :: dt, sumAreaDeltaS, sumArea, avgDeltaS, deltaS, sumAreaDeltaSGlobal, sumAreaGlobal

        type (block_type), pointer :: block

        type (dm_info) :: dminfo

        type (mpas_pool_type), pointer :: forcingPool
        type (mpas_pool_type), pointer :: meshPool
        type (mpas_pool_type), pointer :: statePool
        type (mpas_pool_type), pointer :: tracersPool
        type (mpas_pool_type), pointer :: surfaceSalinityMonthlyForcing
        type (mpas_pool_type), pointer :: tracersSurfaceRestoringFieldsPool

        real (kind=RKIND), pointer :: salinity_restoring_max_difference, salinity_restoring_constant_piston_velocity

        real (kind=RKIND), dimension(:), pointer :: &
           surfaceSalinityMonthlyClimatologyValue, iceFraction, areaCell

        real (kind=RKIND), dimension(:,:), pointer :: &
           activeTracersPistonVelocity,  &
           activeTracersSurfaceRestoringValue

        real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

        integer, pointer :: nCellsSolve, indexSalinity, indexSalinitySurfaceRestoringValue, indexSalinityPistonVelocity
        integer :: iCell, timeLevel

        character(len=strKIND) :: &
           forcingIntervalMonthly,  &
           forcingReferenceTimeMonthly

        character(len=strKind), pointer :: &
           config_forcing_restart_file
  
        logical, pointer :: &
           config_do_restart

        integer, parameter :: nSums = 2
        real (kind=RKIND), dimension(nSums) :: reductions, sums

! initialize monthly forcing to be read from file

        if (firstTimeStep) then

           forcingIntervalMonthly = "0000-01-00_00:00:00"
           forcingReferenceTimeMonthly = "0000-01-15_00:00:00"

           call MPAS_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
           call MPAS_pool_get_config(domain % configs, 'config_forcing_restart_file', config_forcing_restart_file)

           call MPAS_forcing_init_group( surfaceSalinityForcingGroup,  &
                "surfaceSalinityMonthlyClimatology", &
                domain, &
                '0000-01-01_00:00:00', &
                '0000-01-01_00:00:00', &
                '0001-00-00_00:00:00', &
                config_do_restart, &
                config_forcing_restart_file)
!        .false.)  ! OPTIONAL?

           call mpas_pool_get_subpool(domain % blocklist % structs, 'surfaceSalinityMonthlyForcing',  &
                                                                     surfaceSalinityMonthlyForcing)
           call mpas_pool_get_array(surfaceSalinityMonthlyForcing, 'surfaceSalinityMonthlyClimatologyValue',   &
                                                                    surfaceSalinityMonthlyClimatologyValue)
           call MPAS_forcing_init_field( domain % streamManager, &
                surfaceSalinityForcingGroup, &
                'surfaceSalinityMonthlyClimatology', &  !  forcing group name
                'surfaceSalinityMonthlyClimatologyValue', &  !  array name
                'surface_salinity_monthly_data', &  !  stream name
                'surfaceSalinityMonthlyForcing',  & !  pool name
                'surfaceSalinityMonthlyClimatologyValue',  &   !  array name
                'linear',  &
                forcingReferenceTimeMonthly,  &
                forcingIntervalMonthly)

           call MPAS_forcing_init_field_data( surfaceSalinityForcingGroup, &
                'surfaceSalinityMonthlyClimatology', &
                domain % streamManager)

        endif  !  first timestep

        call MPAS_pool_get_config(domain%configs, 'config_dt', config_dt)
        call MPAS_pool_get_config(domain%configs, 'config_salinity_restoring_max_difference',  &
                                                          salinity_restoring_max_difference)
        call MPAS_pool_get_config(domain%configs, 'config_salinity_restoring_constant_piston_velocity',  &
                                                          salinity_restoring_constant_piston_velocity)

        call mpas_set_timeInterval(timeStepSurfaceSalinity,timeString=config_dt)
        call mpas_get_timeInterval(timeStepSurfaceSalinity,dt=dt)

        call mpas_pool_get_subpool(domain % blocklist % structs, 'surfaceSalinityMonthlyForcing',  &
                                                                  surfaceSalinityMonthlyForcing)
        call mpas_pool_get_array(surfaceSalinityMonthlyForcing, 'surfaceSalinityMonthlyClimatologyValue',   &
                                                                 surfaceSalinityMonthlyClimatologyValue)

        call MPAS_forcing_get_forcing(surfaceSalinityForcingGroup, &
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
!maltrud make sure of timeLevel
        timeLevel = 1
        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)

        call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)
        call mpas_pool_get_dimension(tracersSurfaceRestoringFieldsPool, 'index_salinitySurfaceRestoringValue',  &
                                                                         indexSalinitySurfaceRestoringValue)
        call mpas_pool_get_dimension(tracersSurfaceRestoringFieldsPool, 'index_salinityPistonVelocity',  &
                                                                         indexSalinityPistonVelocity)

        call mpas_pool_get_subpool(block % structs, 'surfaceSalinityMonthlyForcing',  &
                                                                  surfaceSalinityMonthlyForcing)
        call mpas_pool_get_array(surfaceSalinityMonthlyForcing, 'surfaceSalinityMonthlyClimatologyValue',   &
                                                                 surfaceSalinityMonthlyClimatologyValue)

        call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersSurfaceRestoringValue', activeTracersSurfaceRestoringValue)
        call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersPistonVelocity', activeTracersPistonVelocity)

        call mpas_pool_get_array(forcingPool, 'iceFraction', iceFraction)

        call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

        do iCell = 1, nCellsSolve
           deltaS = surfaceSalinityMonthlyClimatologyValue(iCell) - activeTracers(indexSalinity,1,iCell)
           if (deltaS >  salinity_restoring_max_difference) deltaS =  salinity_restoring_max_difference
           if (deltaS < -salinity_restoring_max_difference) deltaS = -salinity_restoring_max_difference
!maltrud debug
!          if (iceFraction(iCell) > 1.e-10_RKIND) deltaS = 0.0_RKIND
           if (iceFraction(iCell) > 1.e-10_RKIND) then
              deltaS = 0.0_RKIND
!maltrud maybe say deltaS = (1-iceFrac)*deltaS for iceFrac < 0.1, but = 0 for > 0.1?

!maltrud debug
           else
           sumAreaDeltaS = sumAreaDeltaS + deltaS*areaCell(iCell)
           sumArea = sumArea + areaCell(iCell)
!maltrud debug
           endif

           activeTracersSurfaceRestoringValue(indexSalinitySurfaceRestoringValue,iCell) = deltaS

           activeTracersPistonVelocity(indexSalinityPistonVelocity,iCell) =  &
                 salinity_restoring_constant_piston_velocity

        enddo

         block => block % next
      end do

!       do global sum to subtract global mean of deltaS
        dminfo = domain % dminfo
        sums(1) = sumAreaDeltaS
        sums(2) = sumArea
        call mpas_dmpar_sum_real_array(dminfo, nSums, sums(1:nSums), reductions(1:nSums))
        sumAreaDeltaSGlobal = reductions(1)
        sumAreaGlobal       = reductions(2)
        avgDeltaS = sumAreaDeltaSGlobal/(sumAreaGlobal + 1.e-20_RKIND)

!maltrud debug
!write(*,*)' sumAreaDeltaS, sumAreaDeltaSGlobal, avgDeltaS = ',sumAreaDeltaS, sumAreaDeltaSGlobal, avgDeltaS

      block => domain % blocklist
      do while (associated(block))

        call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
        call mpas_pool_get_array(forcingPool, 'iceFraction', iceFraction)
        call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields',tracersSurfaceRestoringFieldsPool)
        call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersSurfaceRestoringValue', activeTracersSurfaceRestoringValue)
        call mpas_pool_get_dimension(tracersSurfaceRestoringFieldsPool, 'index_salinitySurfaceRestoringValue',  &
                                                                         indexSalinitySurfaceRestoringValue)

        do iCell = 1, nCellsSolve
!maltrud debug
!          if (iceFraction(iCell) <= 1.e-10_RKIND)  &
           if (iceFraction(iCell) <= 1.e-10_RKIND)  &
           activeTracersSurfaceRestoringValue(indexSalinitySurfaceRestoringValue,iCell) =  &
              activeTracersSurfaceRestoringValue(indexSalinitySurfaceRestoringValue,iCell) - avgDeltaS
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

      character(len=strKind), pointer :: config_forcing_restart_file

      call MPAS_pool_get_config(domain % configs, "config_forcing_restart_file", config_forcing_restart_file)

      call MPAS_forcing_write_restart_times(domain, surfaceSalinityForcingGroup, config_forcing_restart_file)

    end subroutine ocn_salinity_restoring_forcing_write_restart!}}}

!***********************************************************************

end module ocn_tracer_surface_restoring

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
