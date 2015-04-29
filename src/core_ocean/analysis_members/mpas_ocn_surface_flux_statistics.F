! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_surface_flux_statistics
!
!> \brief MPAS ocean analysis member: surface statistics
!> \author Todd Ringler
!> \date   April 24, 2015
!> \details
!>  MPAS ocean analysis member: surface statistics
!
!-----------------------------------------------------------------------

module ocn_surface_flux_statistics

   use mpas_grid_types
   use mpas_timer
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_diagnostics_routines

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

   public :: ocn_setup_packages_surface_flux_statistics, &
             ocn_init_surface_flux_statistics, &
             ocn_compute_surface_flux_statistics, &
             ocn_restart_surface_flux_statistics, &
             ocn_finalize_surface_flux_statistics

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type (timer_node), pointer :: amSurfaceStatisticsTimer

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_setup_packages_surface_flux_statistics
!
!> \brief   Set up packages for MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details 
!>  This routine is intended to configure the packages for this MPAS
!>   ocean analysis member
!
!-----------------------------------------------------------------------

   subroutine ocn_setup_packages_surface_flux_statistics(configPool, packagePool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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

      logical, pointer :: amSurfaceFluxStatisticsPkgActive

      err = 0

      call mpas_pool_get_package(packagePool, 'amSurfaceFluxStatisticsPkgActive', amSurfaceFluxStatisticsPkgActive)

      ! turn on package for this analysis member
      amSurfaceFluxStatisticsPkgActive = .true.

      write(6,*) ' active surface statistics'

   end subroutine ocn_setup_packages_surface_flux_statistics!}}}


!***********************************************************************
!
!  routine ocn_init_surface_flux_statistics
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_surface_flux_statistics(domain, err)!{{{

      use mpas_packages

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

      write(6,*) ' surf stat inited'

   end subroutine ocn_init_surface_flux_statistics!}}}

!***********************************************************************
!
!  routine ocn_compute_surface_flux_statistics
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_surface_flux_statistics(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: amSurfaceStatisticsPool
      type (mpas_pool_type), pointer :: amSurfaceStatisticsLocationPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: forcingPool

      real (kind=RKIND), dimension(:,:), pointer ::   minSurfaceStatistics, &
                                                      maxSurfaceStatistics, &
                                                      avgSurfaceStatistics

      ! pointers to data in pools to be analyzed
      real (kind=RKIND), dimension(:), pointer :: ssh
      real (kind=RKIND), dimension(:), pointer :: surfaceWindStressMagnitude
      real (kind=RKIND), dimension(:,:), pointer :: kineticEnergyCell

!latentHeatFlux packages="bulkForcing"
!sensibleHeatFlux packages="bulkForcing"
!longWaveHeatFluxUp packages="bulkForcing"
!longWaveHeatFluxDown packages="bulkForcing"
!seaIceHeatFlux packages="bulkForcing"
!shortWaveHeatFlux packages="bulkForcing"

!evaporationFlux packages="bulkForcing"
!seaIceFreshWaterFlux packages="bulkForcing"
!riverRunoffFlux packages="bulkForcing"
!iceRunoffFlux packages="bulkForcing"
!rainFlux packages="bulkForcing"
!snowFlux packages="bulkForcing"
!seaIceEnergy  (packages="frazilIce)

!surfaceThicknessFlux
!surfaceTracerFlux  (temp, salinity)
!penetrativeTemperatureFlux

!seaIceSalinityFlux
!surfaceWindStressMagnitude
!windStressZonal packages="bulkForcing"
!windStressMeridional packages="bulkForcing"
!ssh
!seaSurfacePressure

      ! pointers to data in mesh pool
      integer, pointer :: nCells, nCellsSolve, nSurfaceFluxStatisticsFields, nOceanRegions
      real (kind=RKIND), dimension(:), pointer ::  areaCell, lonCell, latCell

      ! scratch space
      type(field2DReal), pointer :: workArrayField
      real (kind=RKIND), dimension(:,:), pointer :: workArray
      type(field1DReal), pointer :: workMaskField, workMinField, workMaxField, workSumField
      real (kind=RKIND), dimension(:), pointer :: workMask, workMin, workMax, workSum

      ! local variables
      integer :: nDefinedDataFields
      integer :: iCell, iRegion, iTracer, err_tmp

      ! package flag
      logical, pointer :: amSurfaceFluxStatisticsPkgActive

      ! assume no error
      err = 0

      ! check to see if package is on
      call mpas_pool_get_package(ocnPackages, 'amSurfaceFluxStatisticsPkgActive', amSurfaceFluxStatisticsPkgActive)
      if ( .not. amSurfaceFluxStatisticsPkgActive ) return

      ! start timer
      call mpas_timer_start("compute_surface_flux_statistics", .false., amSurfaceStatisticsTimer)

      ! set highest level pointer
      dminfo = domain % dminfo

      ! loop over all ocean regions
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nOceanRegions', nOceanRegions)
      do iRegion=1,nOceanRegions

      ! get pointers to analysis member arrays
      call mpas_pool_get_subpool(domain % blocklist % structs, 'amSurfaceStatistics', amSurfaceStatisticsPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'amSurfaceStatisticsLocation', amSurfaceStatisticsLocationPool)
      call mpas_pool_get_array(amSurfaceStatisticsPool, 'minSurfaceStatistics', minSurfaceStatistics)
      call mpas_pool_get_array(amSurfaceStatisticsPool, 'maxSurfaceStatistics', maxSurfaceStatistics)
      call mpas_pool_get_array(amSurfaceStatisticsPool, 'avgSurfaceStatistics', avgSurfaceStatistics)

      ! get pointers to scratch variables
      call mpas_pool_get_subpool(domain % blocklist % structs, 'surfaceStatisticsScratch', scratchPool)
      call mpas_pool_get_field(scratchPool, 'workArray', workArrayField)
      call mpas_pool_get_field(scratchPool, 'workMask', workMaskField)
      call mpas_pool_get_field(scratchPool, 'workMin', workMinField)
      call mpas_pool_get_field(scratchPool, 'workMax', workMaxField)
      call mpas_pool_get_field(scratchPool, 'workSum', workSumField)
      call mpas_allocate_scratch_field(workArrayField, .true.)
      call mpas_allocate_scratch_field(workMaskField, .true.)
      call mpas_allocate_scratch_field(workMinField, .true.)
      call mpas_allocate_scratch_field(workMaxField, .true.)
      call mpas_allocate_scratch_field(workSumField, .true.)
      workArray => workArrayField % array
      workMask => workMaskField % array
      workMin => workMinField % array
      workMax => workMaxField % array
      workSum => workSumField % array

      ! loop over blocks
      ! NOTE: code is not valid for multiple blocks !
      block => domain % blocklist
      do while (associated(block))
         ! get pointers to pools
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

         ! get pointers to mesh
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nSurfaceFluxStatisticsFields', nSurfaceFluxStatisticsFields)
         call mpas_pool_get_dimension(block % dimensions, 'nOceanRegions', nOceanRegions)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
         call mpas_pool_get_array(meshPool, 'latCell', latCell)

         ! test to make sure the arrays are big enough
         nDefinedDataFields = size(avgSurfaceStatistics(:,1))
         if (nDefinedDataFields > nSurfaceFluxStatisticsFields) then
              ! call error/finalize here
         endif

    ! HERE -- how to get num_tracers
    !    call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)

         ! get pointers to data that will be analyzed
         call mpas_pool_get_array(statePool, 'ssh', ssh)  
         call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell)
         call mpas_pool_get_array(forcingPool, 'surfaceWindStressMagnitude', surfaceWindStressMagnitude)

         ! copy data into work array
         workArray(1,:) = areaCell(:)
         workArray(2,:) = ssh(:)
         workArray(3,:) = kineticEnergyCell(1,:)

         call compute_mask(nCells, nCellsSolve, iRegion, lonCell, latCell, workMask)
         call compute_statistics(nDefinedDataFields, nCells, workArray, workMask, workMin, workMax, workSum)

         block => block % next
      end do

      ! communication
      call mpas_dmpar_sum_real_array(dminfo, nDefinedDataFields, workSum(1:nDefinedDataFields), avgSurfaceStatistics(1:nDefinedDataFields,iRegion))
      call mpas_dmpar_min_real_array(dminfo, nDefinedDataFields, workMin(1:nDefinedDataFields), minSurfaceStatistics(1:nDefinedDataFields,iRegion))
      call mpas_dmpar_min_real_array(dminfo, nDefinedDataFields, workMax(1:nDefinedDataFields), maxSurfaceStatistics(1:nDefinedDataFields,iRegion))

      ! normalize averages
      avgSurfaceStatistics(:,iRegion) = avgSurfaceStatistics(:,iRegion) / avgSurfaceStatistics(1,iRegion)

      end do ! iRegion

      ! deallocate scratch fields
      call mpas_deallocate_scratch_field(workArrayField, .true.)
      call mpas_deallocate_scratch_field(workMaskField, .true.)
      call mpas_deallocate_scratch_field(workMinField, .true.)
      call mpas_deallocate_scratch_field(workMaxField, .true.)
      call mpas_deallocate_scratch_field(workSumField, .true.)

      ! stop timer
      call mpas_timer_stop("surface_flux_statistics", amSurfaceStatisticsTimer)

   contains

   subroutine compute_mask(nCells, nCellsSolve, iRegion, lonCell, latCell, workMask)
   ! this subroutines produces a 0/1 mask that is multiplied with workArray to 
   ! allow for min/max/avg to represent specific regions of the ocean domain
   integer, intent(in) :: nCells, nCellsSolve, iRegion
   real(kind=RKIND), dimension(:), intent(out) :: lonCell, latCell
   real(kind=RKIND), dimension(:), intent(out) :: workMask
   integer :: iCell
   real :: rtd

   rtd = 180.0_RKIND / (4.0*atan(1.0))
   workMask(:) = 0.0_RKIND
   do iCell=1,nCellsSolve
      workMask(iCell) = 1.0_RKIND
   enddo

   if (iRegion.eq.1) then
      ! Arctic
      do iCell=1,nCellsSolve
        if(latCell(iCell).lt. 60.0*rtd) workMask(iCell) = 0.0_RKIND
      enddo
   elseif (iRegion.eq.2) then
      ! Equatorial
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt. 15.0*rtd) workMask(iCell) = 0.0_RKIND
        if(latCell(iCell).lt.-15.0*rtd) workMask(iCell) = 0.0_RKIND
      enddo
   elseif (iRegion.eq.3) then
      ! Southern Ocean
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt.-50.0*rtd) workMask(iCell) = 0.0_RKIND
      enddo
   elseif (iRegion.eq.4) then
      ! Nino 3
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt.  5.0*rtd) workMask(iCell) = 0.0_RKIND
        if(latCell(iCell).lt. -5.0*rtd) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).lt.210.0*rtd) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).gt.270.0*rtd) workMask(iCell) = 0.0_RKIND
      enddo
   elseif (iRegion.eq.5) then
      ! Nino 4
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt.  5.0*rtd) workMask(iCell) = 0.0_RKIND
        if(latCell(iCell).lt. -5.0*rtd) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).lt.160.0*rtd) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).gt.210.0*rtd) workMask(iCell) = 0.0_RKIND
      enddo
   elseif (iRegion.eq.6) then
      ! Nino 3.4
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt.  5.0*rtd) workMask(iCell) = 0.0_RKIND
        if(latCell(iCell).lt. -5.0*rtd) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).lt.190.0*rtd) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).gt.240.0*rtd) workMask(iCell) = 0.0_RKIND
      enddo
   else
      ! global (do nothing!)
   endif

   end subroutine compute_mask


   subroutine compute_statistics(nDefinedDataFields, nCells, workArray, workMask, workMin, workMax, workSum)
   ! this subroutines does the actual summing, min, max, masking ect
   ! this hides the messy code from the high-level subroutine

   integer, intent(in) :: nDefinedDataFields, nCells
   real(kind=RKIND), dimension(:,:), intent(in) :: workArray
   real(kind=RKIND), dimension(:), intent(in) :: workMask
   real(kind=RKIND), dimension(:), intent(out) :: workMin, workMax, workSum
   integer :: iCell, iData

   workSum = 0.0
   do iCell=1,nCells
    workSum(1) = workSum(1) + workArray(1,iCell)
    workSum(2:nDefinedDataFields) = workSum(2:nDefinedDataFields) + workArray(1,iCell)*workArray(2:nDefinedDataFields,iCell)*workMask(iCell)
   enddo

   do iData=1,nDefinedDataFields
      workMin(iData) = minval(workArray(iData,:)*workMask(:))
      workMax(iData) = maxval(workArray(iData,:)*workMask(:))
   enddo

   end subroutine compute_statistics

   end subroutine ocn_compute_surface_flux_statistics!}}}

!***********************************************************************
!
!  routine ocn_restart_surface_flux_statistics
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_surface_flux_statistics(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

   end subroutine ocn_restart_surface_flux_statistics!}}}

!***********************************************************************
!
!  routine ocn_finalize_surface_flux_statistics
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_surface_flux_statistics(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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

      err = 0

   end subroutine ocn_finalize_surface_flux_statistics!}}}

end module ocn_surface_flux_statistics

! vim: foldmethod=marker
