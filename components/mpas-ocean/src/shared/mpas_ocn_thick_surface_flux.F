! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_thick_surface_flux
!
!> \brief MPAS ocean surface fluxes for thickness
!> \author Doug Jacobsen
!> \date   12/17/12
!> \details
!>  This module contains the routine for computing
!>  tendencies for thickness from surface fluxes
!
!-----------------------------------------------------------------------

module ocn_thick_surface_flux

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_forcing
   use mpas_stream_manager
   use ocn_constants
   use ocn_config
   use ocn_framework_forcing
   use ocn_mesh

   use ocn_forcing

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

   public :: ocn_thick_surface_flux_tend, &
             ocn_thick_surface_flux_init, &
             ocn_get_hosingData,  &
             ocn_hosing_write_restart

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: surfaceThicknessFluxOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_thick_surface_flux_tend
!
!> \brief   Computes tendency term from horizontal advection of thickness
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \details
!>  This routine computes the horizontal advection tendency for
!>  thicknes based on current state and user choices of forcings.
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_surface_flux_tend(transmissionCoefficients, transmissionCoefficientsRunoff, &
      surfaceThicknessFlux, surfaceThicknessFluxRunoff, tend, err)!{{{
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         transmissionCoefficients,     &!< Input: Coefficients for the transmission of surface fluxes
         transmissionCoefficientsRunoff !< Input: Coefficients for the transmission of surface fluxes due to river runoff

      real (kind=RKIND), dimension(:), intent(in) :: &
         surfaceThicknessFlux,       &!< Input: surface flux of thickness
         surfaceThicknessFluxRunoff   !< Input: surface flux of thickness due to river runoff


      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: thickness tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, k

      real (kind=RKIND) :: remainingFlux, remainingFluxRunoff

      err = 0

      if (.not. surfaceThicknessFluxOn) return

      call mpas_timer_start("thick surface flux")

#ifdef MPAS_OPENACC
      !$acc enter data copyin(transmissionCoefficients, transmissionCoefficientsRunoff)

      !$acc parallel loop &
      !$acc     present(tend, surfaceThicknessFlux, surfaceThicknessFluxRunoff, &
      !$acc     transmissionCoefficients, transmissionCoefficientsRunoff, minLevelCell, maxLevelCell) &
      !$acc     private(k, remainingFlux, remainingFluxRunoff)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(remainingFlux, remainingFluxRunoff, k)
#endif
      do iCell = 1, nCellsOwned
        remainingFlux = 1.0_RKIND
        remainingFluxRunoff = 1.0_RKIND
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          remainingFlux = remainingFlux - transmissionCoefficients(k, iCell)
          remainingFluxRunoff = remainingFluxRunoff - transmissionCoefficientsRunoff(k, iCell)

          tend(k, iCell) = tend(k, iCell) + surfaceThicknessFlux(iCell) * transmissionCoefficients(k, iCell)  &
                                          + surfaceThicknessFluxRunoff(iCell) * transmissionCoefficientsRunoff(k, iCell)
        end do

        if(maxLevelCell(iCell) > 0 .and. remainingFlux > 0.0_RKIND) then
          tend(maxLevelCell(iCell), iCell) = tend(maxLevelCell(iCell), iCell) + remainingFlux * surfaceThicknessFlux(iCell)
        end if

        if(maxLevelCell(iCell) > 0 .and. remainingFluxRunoff > 0.0_RKIND) then
          tend(maxLevelCell(iCell), iCell) = tend(maxLevelCell(iCell), iCell) &
             + remainingFluxRunoff * surfaceThicknessFluxRunoff(iCell)
        end if
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(transmissionCoefficients, transmissionCoefficientsRunoff)
#endif

      call mpas_timer_stop("thick surface flux")

   !--------------------------------------------------------------------

   end subroutine ocn_thick_surface_flux_tend!}}}

!***********************************************************************
!
!  routine ocn_thick_surface_flux_init
!
!> \brief   Initializes ocean horizontal thickness surface fluxes
!> \author  Doug Jacobsen
!> \date    12/17/12
!> \details
!>  This routine initializes quantities related to thickness
!>  surface fluxes in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_surface_flux_init(domain,err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

!-----------------------------------------------------------------
! variables for freshwater hosing
!-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: forcingPool

      real (kind=RKIND), dimension(:), pointer :: &
         hosingFlux,     & ! input flux in forcing pool
         freshWaterHosing  ! input flux in stream

      character(len=strKIND) :: &
         forcingIntervalAnnual,  &
         forcingReferenceTimeAnnual

      type (MPAS_Time_Type) :: currTime
      integer :: ierr
      character(len=StrKIND)  :: timeStamp

      err = 0

      surfaceThicknessFluxOn = .true.

      if (config_disable_thick_sflux) then
         surfaceThicknessFluxOn = .false.
      end if

!-----------------------------------------------------------------
! freshwater hosing
!-----------------------------------------------------------------

      if (.not. config_use_hosing_flux) return

!-----------------------------------------------------------------------
! initialize annual forcing to be read from file
!-----------------------------------------------------------------------

      forcingIntervalAnnual = "0001-00-00_00:00:00"
      forcingReferenceTimeAnnual = "1850-07-01_00:00:00"

      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)

      call MPAS_forcing_init_group( forcingGroupHead,  &
                "hosingAnnualFreshwaterForcing", &
                domain, &
                trim(timeStamp), &
                '1850-07-01_00:00:00', &
                '0600-00-00_00:00:00', &
                config_do_restart)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
      call mpas_pool_get_array(forcingPool, 'freshWaterHosing', freshWaterHosing)

      call MPAS_forcing_init_field( domain % streamManager, &
                forcingGroupHead, &
                'hosingAnnualFreshwaterForcing', &
                'freshWaterHosing', &
                'hosing_annual_flux', &
                'forcing',  &
                'freshWaterHosing',  &
                'linear',  &
                forcingReferenceTimeAnnual,  &
                forcingIntervalAnnual)

      call MPAS_forcing_init_field_data( forcingGroupHead, &
        'hosingAnnualFreshwaterForcing',  &
         domain % streamManager,  &
         config_do_restart, &
         .false.)

   !--------------------------------------------------------------------

   end subroutine ocn_thick_surface_flux_init!}}}

!***********************************************************************
!***********************************************************************
!
!  routine ocn_get_hosingData
!
!> \brief   retrieve data needed for freshwater hosing surface fluxes
!> \author  Mathew Maltrud
!> \date    11/30/2022
!> \details
!>  This routine calls mpas_forcing routines to acquire needed hsoing data and interpolates
!>    between time levels.  directly copied from ocn_get_CFCData.
!
!-----------------------------------------------------------------------

    subroutine ocn_get_hosingData( streamManager, &
        domain, &
        simulationClock, &
        firstTimeStep) !{{{

        type (MPAS_streamManager_type), intent(inout) :: streamManager

        type (domain_type) :: domain
        type (MPAS_timeInterval_type) :: timeStepHosing
        type (MPAS_clock_type) :: simulationClock

        logical, intent(in) :: firstTimeStep
        real(kind=RKIND) :: dt

        type (mpas_pool_type), pointer :: forcingPool
        type (mpas_pool_type), pointer :: meshPool

        real (kind=RKIND), dimension(:), pointer :: &
         hosingFlux,         &
         freshWaterHosing

        integer, pointer :: nCells
        integer :: iCell

        call mpas_set_timeInterval(timeStepHosing,timeString=config_dt)
        call mpas_get_timeInterval(timeStepHosing,dt=dt)

        if (firstTimeStep .and. config_do_restart) then
          call MPAS_forcing_get_forcing(forcingGroupHead, &
             'hosingAnnualFreshwaterForcing', streamManager, 0.0_RKIND)
        else
          call MPAS_forcing_get_forcing(forcingGroupHead, &
             'hosingAnnualFreshwaterForcing', streamManager, dt)
        endif

        call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
        call mpas_pool_get_array(forcingPool, 'hosingFlux', hosingFlux)
        call mpas_pool_get_array(forcingPool, 'freshWaterHosing', freshWaterHosing)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        do iCell = 1, nCells
           hosingFlux(iCell) = freshWaterHosing(iCell)
        enddo
    end subroutine ocn_get_hosingData!}}}

!***********************************************************************
!
!  routine ocn_hosing_write_restart
!
!> \brief   writes restart timestamp for hosing data to be read in on future restart
!> \author  Mathew Maltrud
!> \date    11/30/2022

!
!-----------------------------------------------------------------------

   subroutine ocn_hosing_write_restart(domain)!{{{

      type(domain_type) :: domain

      call MPAS_forcing_write_restart_times(forcingGroupHead)

    end subroutine ocn_hosing_write_restart!}}}

!-----------------------------------------------------------------------

end module ocn_thick_surface_flux

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
