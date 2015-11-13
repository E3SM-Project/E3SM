! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_regional_statistics
!
!> \brief MPAS sea ice analysis mode member: regional_statistics
!> \author Adrian K. Turner
!> \date   6th September 2015
!> \details
!>  MPAS sea ice analysis mode member: regional_statistics
!>  calculates regional statistics
!>
!-----------------------------------------------------------------------

module cice_regional_statistics

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

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

   public :: cice_bootstrap_regional_statistics, &
             cice_init_regional_statistics, &
             cice_precompute_regional_statistics, &
             cice_compute_regional_statistics, &
             cice_restart_regional_statistics, &
             cice_finalize_regional_statistics

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine cice_bootstrap_regional_statistics
!
!> \brief   Bootstrap MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    10th November 2015
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_bootstrap_regional_statistics(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

   end subroutine cice_bootstrap_regional_statistics!}}}

!***********************************************************************
!
!  routine cice_init_regional_statistics
!
!> \brief   Initialize MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    6th September 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_init_regional_statistics(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

      type (block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           regionsPool

      real(kind=RKIND), dimension(:), pointer :: &
           latCell
      
      integer, dimension(:,:), pointer :: &
           regionMask

      integer, pointer :: &
           nCells

      integer :: &
           iCell

      err = 0

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)
         
         call MPAS_pool_get_dimension(meshPool, "nCells", nCells)
         
         call MPAS_pool_get_array(meshPool, "latCell", latCell)
         
         call MPAS_pool_get_array(regionsPool, "regionMask", regionMask)

         ! initialize to zero everywhere
         regionMask = 0

         ! global
         regionMask(1,:) = 1

         ! northern hemisphere
         do iCell = 1, nCells
            if (latCell(iCell) >= 0.0_RKIND) &
                 regionMask(2,iCell) = 1
         enddo ! iCell

         ! southern hemisphere
         do iCell = 1, nCells
            if (latCell(iCell) < 0.0_RKIND) &   
                 regionMask(3,iCell) = 1
         enddo ! iCell

         block => block % next
      enddo

   end subroutine cice_init_regional_statistics!}}}

!***********************************************************************
!
!  routine cice_precompute_regional_statistics
!
!> \brief   Compute MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    9th September 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_precompute_regional_statistics(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

      err = 0

   end subroutine cice_precompute_regional_statistics!}}}

!***********************************************************************
!
!  routine cice_compute_regional_statistics
!
!> \brief   Compute MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    6th September 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_compute_regional_statistics(domain, instance, timeLevel, err)!{{{

      use cice_mesh, only: &
           cice_interpolate_vertex_to_cell

      use ice_constants_colpkg, only: &
           awtvdr, &
           awtidr, &
           awtvdf, &
           awtidf, &
           rhoi, &
           rhos

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

      type (dm_info) :: &
           dminfo

      type (block_type), pointer :: &
           block

      type (mpas_pool_type), pointer :: &
           meshPool, &
           regionsPool, &
           tracersAggregatePool, &
           shortwavePool, &
           velocitySolverPool, &
           boundaryPool, &
           regionalStatisticsAMPool

      integer, pointer :: &
           nRegions, &
           nCellsSolve

      integer :: &
           iRegion, &
           iCell

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           iceAreaCell, &
           iceVolumeCell, &
           snowVolumeCell, &
           solarZenithAngleCosine, &
           albedoVisibleDirectCell, &
           albedoIRDirectCell, &
           albedoVisibleDiffuseCell, &
           albedoIRDiffuseCell, &
           icePressure, &
           uVelocity, &
           vVelocity, &
           uVelocityCell, &
           vVelocityCell

      integer, dimension(:,:), pointer :: &
           regionMask

      real(kind=RKIND), dimension(:), pointer :: &
           totalIceArea, &
           totalIceExtent, &
           totalIceVolume, &
           totalSnowVolume, &
           totalKineticEnergy, &
           rmsIceSpeed, &
           averageAlbedo, &
           maximumIceVolume, &
           maximumIcePressure

      real(kind=RKIND), dimension(:), allocatable :: &
           globalSumsIn, &
           globalMaxsIn, &
           globalSumsOut, &
           globalMaxsOut
      
      integer :: &
           nSums, &
           nMaxs, &
           iSum, &
           iMax

      real(kind=RKIND), pointer :: &
           iceExtentLimit

      real(kind=RKIND), parameter :: &
           m2_to_km2 = 1.0e-6_RKIND, &
           m3_to_km3 = 1.0e-9_RKIND, &
           Nm_to_kNm = 1.0e-3_RKIND

      err = 0

      dminfo = domain % dminfo

      call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

      nSums = 8
      nMaxs = 2

      allocate(globalSumsIn(nRegions*nSums))
      allocate(globalMaxsIn(nRegions*nMaxs))
      globalSumsIn = 0.0_RKIND
      globalMaxsIn = -1.0e30_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_config(block % configs, "config_AM_regionalStatistics_ice_extent_limit", iceExtentLimit)

         call MPAS_pool_get_subpool(block % structs, 'mesh', meshPool)
         call MPAS_pool_get_subpool(block % structs, 'regions', regionsPool)
         call MPAS_pool_get_subpool(block % structs, 'tracers_aggregate', tracersAggregatePool)
         call MPAS_pool_get_subpool(block % structs, 'shortwave', shortwavePool)
         call MPAS_pool_get_subpool(block % structs, 'velocity_solver', velocitySolverPool)
         call MPAS_pool_get_subpool(block % structs, 'boundary', boundaryPool)
         call MPAS_pool_get_subpool(block % structs, 'regionalStatisticsAM', regionalStatisticsAMPool)

         call MPAS_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

         call MPAS_pool_get_array(meshPool, 'areaCell', areaCell)         

         call MPAS_pool_get_array(tracersAggregatePool, 'iceAreaCell', iceAreaCell)
         call MPAS_pool_get_array(tracersAggregatePool, 'iceVolumeCell', iceVolumeCell)
         call MPAS_pool_get_array(tracersAggregatePool, 'snowVolumeCell', snowVolumeCell)

         call MPAS_pool_get_array(shortwavePool, "solarZenithAngleCosine", solarZenithAngleCosine)
         call MPAS_pool_get_array(shortwavePool, "albedoVisibleDirectCell", albedoVisibleDirectCell)
         call MPAS_pool_get_array(shortwavePool, "albedoIRDirectCell", albedoIRDirectCell)
         call MPAS_pool_get_array(shortwavePool, "albedoVisibleDiffuseCell", albedoVisibleDiffuseCell)
         call MPAS_pool_get_array(shortwavePool, "albedoIRDiffuseCell", albedoIRDiffuseCell)

         call MPAS_pool_get_array(velocitySolverPool, "icePressure", icePressure)
         call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
         call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)

         call MPAS_pool_get_array(regionsPool, 'regionMask', regionMask)

         call MPAS_pool_get_array(regionalStatisticsAMPool, 'uVelocityCell', uVelocityCell)
         call MPAS_pool_get_array(regionalStatisticsAMPool, 'vVelocityCell', vVelocityCell)

         ! cell centre velocity for velocity statistics
         call cice_interpolate_vertex_to_cell(meshPool, boundaryPool, uVelocityCell, uVelocity)
         call cice_interpolate_vertex_to_cell(meshPool, boundaryPool, vVelocityCell, vVelocity)

         do iCell = 1, nCellsSolve
            do iRegion = 1, nRegions

               if (regionMask(iRegion,iCell) == 1) then

                  ! total ice area
                  iSum = (iRegion-1) * nSums + 1
                  globalSumsIn(iSum) = globalSumsIn(iSum) + iceAreaCell(iCell) * areaCell(iCell)

                  ! total ice extent
                  if (iceAreaCell(iCell) > iceExtentLimit) then

                     iSum = (iRegion-1) * nSums + 2
                     globalSumsIn(iSum) = globalSumsIn(iSum) + areaCell(iCell)

                  endif

                  ! total ice volume
                  iSum = (iRegion-1) * nSums + 3
                  globalSumsIn(iSum) = globalSumsIn(iSum) + iceVolumeCell(iCell) * areaCell(iCell)

                  ! total snow volume
                  iSum = (iRegion-1) * nSums + 4
                  globalSumsIn(iSum) = globalSumsIn(iSum) + snowVolumeCell(iCell) * areaCell(iCell)

                  ! kinetic energy
                  iSum = (iRegion-1) * nSums + 5
                  globalSumsIn(iSum) = globalSumsIn(iSum) + 0.5_RKIND * areaCell(iCell) * &
                       (snowVolumeCell(iCell) * rhos + iceVolumeCell(iCell) * rhoi) * &
                       (uVelocityCell(iCell)**2 + vVelocityCell(iCell)**2)

                  ! total mass (for RMS ice speed)
                  iSum = (iRegion-1) * nSums + 6
                  globalSumsIn(iSum) = globalSumsIn(iSum) + areaCell(iCell) * &
                       (snowVolumeCell(iCell) * rhos + iceVolumeCell(iCell) * rhoi)

                  ! average albedo
                  if (solarZenithAngleCosine(iCell) > 0.0_RKIND) then

                     iSum = (iRegion-1) * nSums + 7
                     globalSumsIn(iSum) = globalSumsIn(iSum) + &
                          areaCell(iCell) * &
                          (awtvdr * albedoVisibleDirectCell(iCell) + &
                           awtidr * albedoIRDirectCell(iCell) + &
                           awtvdf * albedoVisibleDiffuseCell(iCell) + &
                           awtidf * albedoIRDiffuseCell(iCell))

                     iSum = (iRegion-1) * nSums + 8
                     globalSumsIn(iSum) = globalSumsIn(iSum) + &
                          areaCell(iCell)

                  endif

                  ! maximum ice volume
                  iMax = (iRegion-1) * nMaxs + 1
                  globalMaxsIn(iMax) = max(globalMaxsIn(iMax), iceVolumeCell(iCell))

                  ! maximum ice pressure
                  iMax = (iRegion-1) * nMaxs + 2
                  globalMaxsIn(iMax) = max(globalMaxsIn(iMax), icePressure(iCell))

               endif ! regionMask == 1

            enddo ! iRegion
         enddo ! iCell

         block => block % next
      enddo

      ! MPI calls
      allocate(globalSumsOut(nRegions*nSums))
      allocate(globalMaxsOut(nRegions*nMaxs))

      call mpas_dmpar_sum_real_array(dminfo, nRegions*nSums, globalSumsIn, globalSumsOut)
      call mpas_dmpar_max_real_array(dminfo, nRegions*nMaxs, globalMaxsIn, globalMaxsOut)

      ! set Registry variables
      call MPAS_pool_get_subpool(domain % blocklist % structs, 'regionalStatisticsAM', regionalStatisticsAMPool)

      call MPAS_pool_get_array(regionalStatisticsAMPool, 'totalIceArea', totalIceArea)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'totalIceExtent', totalIceExtent)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'totalIceVolume', totalIceVolume)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'totalSnowVolume', totalSnowVolume)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'totalKineticEnergy', totalKineticEnergy)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'rmsIceSpeed', rmsIceSpeed)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'averageAlbedo', averageAlbedo)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'maximumIceVolume', maximumIceVolume)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'maximumIcePressure', maximumIcePressure)

      do iRegion = 1, nRegions

         totalIceArea(iRegion)       = globalSumsOut((iRegion-1) * nSums + 1)
         totalIceExtent(iRegion)     = globalSumsOut((iRegion-1) * nSums + 2)
         totalIceVolume(iRegion)     = globalSumsOut((iRegion-1) * nSums + 3)
         totalSnowVolume(iRegion)    = globalSumsOut((iRegion-1) * nSums + 4)
         totalKineticEnergy(iRegion) = globalSumsOut((iRegion-1) * nSums + 5)
         rmsIceSpeed(iRegion)        = globalSumsOut((iRegion-1) * nSums + 5) / &
                                       max(globalSumsOut((iRegion-1) * nSums + 6), 1e-11_RKIND)
         averageAlbedo(iRegion)      = globalSumsOut((iRegion-1) * nSums + 7) / &
                                       max(globalSumsOut((iRegion-1) * nSums + 8), 1e-11_RKIND)
         maximumIceVolume(iRegion)   = globalMaxsOut((iRegion-1) * nMaxs + 1)
         maximumIcePressure(iRegion) = globalMaxsOut((iRegion-1) * nMaxs + 2)

      enddo ! iRegion

      ! unit conversions
      totalIceArea    = totalIceArea    * m2_to_km2
      totalIceExtent  = totalIceExtent  * m2_to_km2
      totalIceVolume  = totalIceVolume  * m3_to_km3
      totalSnowVolume = totalSnowVolume * m3_to_km3

      maximumIcePressure = maximumIcePressure * Nm_to_kNm ! N/m -> kN/m

      deallocate(globalSumsOut)
      deallocate(globalMaxsOut)

      deallocate(globalSumsIn)
      deallocate(globalMaxsIn)

   end subroutine cice_compute_regional_statistics!}}}

!***********************************************************************
!
!  routine cice_restart_regional_statistics
!
!> \brief   Save restart for MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    6th September 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_restart_regional_statistics(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

   end subroutine cice_restart_regional_statistics!}}}

!***********************************************************************
!
!  routine cice_finalize_regional_statistics
!
!> \brief   Finalize MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    6th September 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_finalize_regional_statistics(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

   end subroutine cice_finalize_regional_statistics!}}}

!-----------------------------------------------------------------------

end module cice_regional_statistics

! vim: foldmethod=marker
