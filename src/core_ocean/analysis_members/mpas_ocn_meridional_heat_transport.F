! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_meridional_heat_transport
!
!> \brief MPAS ocean analysis core member: meridional_heat_transport
!> \author Mark Petersen
!> \date   March 2014
!> \details
!>  MPAS ocean analysis core member: meridional_heat_transport
!>  Compute zonal means of selected variables
!
!-----------------------------------------------------------------------

module ocn_meridional_heat_transport

   use mpas_derived_types
   use mpas_pool_routines
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

   public :: ocn_init_meridional_heat_transport, &
             ocn_compute_meridional_heat_transport, &
             ocn_restart_meridional_heat_transport, &
             ocn_finalize_meridional_heat_transport

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: nMerHeatTransBinsUsed

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_meridional_heat_transport
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2014
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_meridional_heat_transport(domain, err)!{{{

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

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meridionalHeatTransportAMPool
      type (mpas_pool_type), pointer :: meshPool

      integer ::  iBin
      integer, pointer ::  nMerHeatTransBins

      real (kind=RKIND) :: binWidth
      ! These are array size 1 because mpas_dmpar_min_real_array calls require arrays.
      real (kind=RKIND), dimension(1) :: minBin, maxBin, minBinDomain, maxBinDomain
      ! the variable used to discriminate cells into Bins (either the y-value or the latitude)
      real (kind=RKIND), dimension(:), pointer ::  binBoundaryMerHeatTrans, binVariable

      !number of latitude bins specified in the config
      integer, pointer :: config_AM_meridionalHeatTransport_num_bins
      !smallest and highest latitude specified in the config
      real (kind=RKIND), pointer :: config_AM_meridionalHeatTransport_min_bin, config_AM_meridionalHeatTransport_max_bin

      !determines if the simulation was run on a sphere or on a plane
      logical, pointer :: on_a_sphere

      !!!! Region variables
      !! region MHT calculation variables
      real (kind=RKIND) :: maskFactor
      integer :: curRegion, i, j, iCell

      !! region arrays/variables
      character (len=STRKIND), dimension(:), pointer :: regionGroupNames
      integer, dimension(:, :), pointer :: regionCellMasks, regionsInGroup
      integer, dimension(:), pointer ::  nRegionsInGroup
      integer, pointer :: nRegions, nRegionGroups, maxRegionsInGroup, nCellsSolve
      real (kind=RKIND), dimension(:), pointer :: minLatRegionLocal, maxLatRegionLocal, &
                                                  minLatRegionGlobal, maxLatRegionGlobal
      character (len=STRKIND), pointer :: additionalRegion

      !! region preliminary variables
      integer :: regionGroupNumber, regionsInAddGroup, regionGroupOffset

      !!region pool
      type (mpas_pool_type), pointer :: regionPool

      !! region dimensions
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegions', nRegions)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegionGroups', nRegionGroups)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'maxRegionsInGroup', maxRegionsInGroup)

      !! region config for MHT
      call mpas_pool_get_config(domain % configs, 'config_AM_meridionalHeatTransport_region_group', &
                                additionalRegion)

      !! get region values
      call mpas_pool_get_subpool(domain % blocklist % structs, 'regions', regionPool)
      call mpas_pool_get_array(regionPool, 'regionsInGroup', regionsInGroup)
      call mpas_pool_get_array(regionPool, 'nRegionsInGroup', nRegionsInGroup)
      call mpas_pool_get_array(regionPool, 'regionGroupNames', regionGroupNames)

      regionGroupOffset = 1
      regionGroupNumber = 0

      if(additionalRegion /= '') then
         !!! region preliminaries
         !!! figure out the region group number that matches the configured additional region's name
         do i = 1, nRegionGroups
            if (regionGroupNames(i) .eq. additionalRegion) then
               regionGroupNumber = i
               !!! determine offset to compensate for several region groups in the
               !!! regions file
               do j = 1, i - 1
                 regionGroupOffset = regionGroupOffset + nRegionsInGroup(j)
               end do
            end if
         end do
      end if

      if(additionalRegion /= '') then
         regionsInAddGroup = nRegionsInGroup(regionGroupNumber)
         allocate(minLatRegionLocal(maxRegionsInGroup))
         allocate(maxLatRegionLocal(maxRegionsInGroup))
         allocate(minLatRegionGlobal(maxRegionsInGroup))
         allocate(maxLatRegionGlobal(maxRegionsInGroup))
      end if
      !!!! END Region variables

      dminfo = domain % dminfo

      err = 0

      minBin =  1.0e34_RKIND
      maxBin = -1.0e34_RKIND

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nMerHeatTransBins', nMerHeatTransBins)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'meridionalHeatTransportAM', meridionalHeatTransportAMPool)

      call mpas_pool_get_config(domain % configs, 'config_AM_meridionalHeatTransport_num_bins', &
                                config_AM_meridionalHeatTransport_num_bins)
      call mpas_pool_get_config(domain % configs, 'config_AM_meridionalHeatTransport_min_bin', &
                                config_AM_meridionalHeatTransport_min_bin)
      call mpas_pool_get_config(domain % configs, 'config_AM_meridionalHeatTransport_max_bin', &
                                config_AM_meridionalHeatTransport_max_bin)

      nMerHeatTransBinsUsed = config_AM_meridionalHeatTransport_num_bins

      call mpas_pool_get_array(meridionalHeatTransportAMPool, 'binBoundaryMerHeatTrans', binBoundaryMerHeatTrans)

      ! Find min and max values of binning variable.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

         ! Get region-specific variables from pools
         if(additionalRegion /= '') then
            call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_array(regionPool, 'regionCellMasks', regionCellMasks)
         end if

         ! Bin by latitude on a sphere, by yCell otherwise.
         if (on_a_sphere) then
            call mpas_pool_get_array(meshPool, 'latCell', binVariable)
         else
            call mpas_pool_get_array(meshPool, 'yCell', binVariable)
         end if

         minBin = min(minBin, minval(binVariable) )
         maxBin = max(maxBin, maxval(binVariable) )

         ! If using regions, iterate through groups and add as needed
         if(additionalRegion /= '') then
            do i = 1, regionsInAddGroup
               curRegion = regionsInGroup(i, regionGroupNumber)
               do iCell = 1, nCellsSolve
                  if (regionCellMasks(curRegion, iCell) .eq. 1) then
                     minLatRegionLocal(i) = min(minLatRegionLocal(i), binVariable(iCell))
                     maxLatRegionLocal(i) = max(maxLatRegionLocal(i), binVariable(iCell))
                  end if
               end do
            end do
         end if
         block => block % next
      end do

      call mpas_dmpar_min_real_array(dminfo, 1, minBin, minBinDomain)
      call mpas_dmpar_max_real_array(dminfo, 1, maxBin, maxBinDomain)

      ! Set up bins.
      binBoundaryMerHeatTrans = -1.0e34_RKIND

      ! Change min and max bin bounds to configuration settings, if applicable.
      if (config_AM_meridionalHeatTransport_min_bin > -1.0e33_RKIND) then
         minBinDomain(1) = config_AM_meridionalHeatTransport_min_bin
      else
         ! use measured min value, but decrease slightly to include least value.
         minBinDomain(1) = minBinDomain(1) - 1.0e-10_RKIND * abs(minBinDomain(1))
      end if

      if (config_AM_meridionalHeatTransport_max_bin > -1.0e33_RKIND) then
         maxBinDomain(1) = config_AM_meridionalHeatTransport_max_bin
      else
         ! use measured max value, but increase slightly to include max value.
         maxBinDomain(1) = maxBinDomain(1) + 1.0e-10_RKIND * abs(maxBinDomain(1))
      end if

      binBoundaryMerHeatTrans(1) = minBinDomain(1)
      binWidth = (maxBinDomain(1) - minBinDomain(1)) / nMerHeatTransBinsUsed

      do iBin = 2, nMerHeatTransBinsUsed
         binBoundaryMerHeatTrans(iBin) = binBoundaryMerHeatTrans(iBin-1) + binWidth
      end do
      binBoundaryMerHeatTrans(nMerHeatTransBinsUsed+1) = binBoundaryMerHeatTrans(nMerHeatTransBinsUsed) + binWidth

   end subroutine ocn_init_meridional_heat_transport!}}}

!***********************************************************************
!
!  routine ocn_compute_meridional_heat_transport
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2014
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_meridional_heat_transport(domain, timeLevel, err)!{{{

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
      type (mpas_pool_type), pointer :: meridionalHeatTransportAMPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      ! REGION POOL
      type (mpas_pool_type), pointer :: regionPool

      integer :: iTracer, k, iCell, kMax, i, iEdge
      integer :: iBin, iField, nMerHeatTransVariables
      integer, pointer :: nCellsSolve, nVertLevels, nMerHeatTransBins, indexTemperature
      integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgeSignOnCell, cellsOnEdge, edgesOnCell

      real (kind=RKIND) :: div_huT
      real (kind=RKIND), dimension(:), pointer ::  areaCell, binVariable, binBoundaryMerHeatTrans, dvEdge
      real (kind=RKIND), dimension(:), pointer ::  merHeatTransLat
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge, normalTransportVelocity
      real (kind=RKIND), dimension(:,:), pointer :: merHeatTransLatZ
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      real (kind=RKIND), dimension(:,:), allocatable :: mht_meridional_integral
      real (kind=RKIND), dimension(:,:,:), allocatable :: sumMerHeatTrans, totalSumMerHeatTrans

      logical, pointer :: on_a_sphere

      !!!! REGION VARIABLES
      real (kind=RKIND) :: maskFactor
      integer :: curRegion, regionGroupOffset, j
      character (len=STRKIND) :: currentName
      real (kind=RKIND), dimension(:,:,:,:), allocatable :: sumMerHeatTransRegion, totalSumMerHeatTransRegion
      real (kind=RKIND), dimension(:,:,:), pointer :: merHeatTransLatZRegion
      real (kind=RKIND), dimension(:,:), pointer :: merHeatTransLatRegion
      character (len=STRKIND), dimension(:), pointer :: regionNames, regionGroupNames
      integer, dimension(:, :), pointer :: regionCellMasks, regionVertexMasks, regionsInGroup
      integer, dimension(:), pointer ::  nRegionsInGroup
      integer, pointer :: nRegions, nRegionGroups, maxRegionsInGroup
      character (len=STRKIND), pointer :: additionalRegion
      integer :: regionGroupNumber, regionsInAddGroup
      real (kind=RKIND), dimension(:,:,:), allocatable :: mht_meridional_integral_region
      !!!! END REGION VARIABLES

      err = 0
      dminfo = domain % dminfo

      call mpas_pool_get_subpool(domain % blocklist % structs, 'meridionalHeatTransportAM', meridionalHeatTransportAMPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)

      nMerHeatTransVariables = 1

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nMerHeatTransBins', nMerHeatTransBins)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meridionalHeatTransportAMPool, 'binBoundaryMerHeatTrans', binBoundaryMerHeatTrans)

      !!!! v
      !! region dimensions
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegions', nRegions)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegionGroups', nRegionGroups)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'maxRegionsInGroup', maxRegionsInGroup)

      !! region config for MHT
      call mpas_pool_get_config(domain % configs, 'config_AM_meridionalHeatTransport_region_group', &
                                additionalRegion)

      !! get region values
      call mpas_pool_get_subpool(domain % blocklist % structs, 'regions', regionPool)
      call mpas_pool_get_array(regionPool, 'regionsInGroup', regionsInGroup)
      call mpas_pool_get_array(regionPool, 'nRegionsInGroup', nRegionsInGroup)
      call mpas_pool_get_array(regionPool, 'regionNames', regionNames)
      call mpas_pool_get_array(regionPool, 'regionGroupNames', regionGroupNames)
      call mpas_pool_get_array(meridionalHeatTransportAMPool, 'merHeatTransLatRegion', merHeatTransLatRegion)
      call mpas_pool_get_array(meridionalHeatTransportAMPool, 'merHeatTransLatZRegion', merHeatTransLatZRegion)

      regionGroupOffset = 1
      regionGroupNumber = 0

      !! if a region is selected, set up region group offsets etc
      if(additionalRegion /= '') then
         !!! region preliminaries
         !!! figure out the region group number that matches the configured additional region's name
         do i = 1, nRegionGroups
            if (regionGroupNames(i) .eq. additionalRegion) then
               regionGroupNumber = i
               !!! determine offset to compensate for several region groups in the
               !!! regions file
               do j = 1, i - 1
                 regionGroupOffset = regionGroupOffset + nRegionsInGroup(j)
               end do
            end if
         end do
      end if

      if(additionalRegion /= '') then
         regionsInAddGroup = nRegionsInGroup(regionGroupNumber)
      end if
      !!!! end region initialization

      !! allocate MHT calculation arrays
      allocate(sumMerHeatTrans(nMerHeatTransVariables,nVertLevels,nMerHeatTransBinsUsed))
      allocate(totalSumMerHeatTrans(nMerHeatTransVariables,nVertLevels,nMerHeatTransBinsUsed))
      allocate(mht_meridional_integral(nVertLevels,nMerHeatTransBinsUsed))

      !! allocate region-specific arrays
      if(additionalRegion /= '') then
!         allocate(totalSumMerHeatTransLatZRegion(nVertLevels, nMerHeatTransBinsUsed + 1, maxRegionsInGroup))
         allocate(sumMerHeatTransRegion(nMerHeatTransVariables, nVertLevels, nMerHeatTransBinsUsed, maxRegionsInGroup))
         allocate(totalSumMerHeatTransRegion(nMerHeatTransVariables, nVertLevels, nMerHeatTransBinsUsed, maxRegionsInGroup))
         allocate(mht_meridional_integral_region(nVertLevels,nMerHeatTransBinsUsed, maxRegionsInGroup))
      end if

      sumMerHeatTrans = 0.0_RKIND
      if(additionalRegion /= '') then
         sumMerHeatTransRegion = 0.0_RKIND
         totalSumMerHeatTrans = 0.0_RKIND
      end if

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)

         call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

         if(additionalRegion /= '') then
            call mpas_pool_get_array(regionPool, 'regionCellMasks', regionCellMasks)
         end if

         ! Bin by latitude on a sphere, by yCell otherwise.
         if (on_a_sphere) then
            call mpas_pool_get_array(meshPool, 'latCell', binVariable)
         else
            call mpas_pool_get_array(meshPool, 'yCell', binVariable)
         end if

         do iCell = 1, nCellsSolve
            kMax = maxLevelCell(iCell)

            if (binVariable(iCell) .lt. binBoundaryMerHeatTrans(1)) cycle

            do iBin = 1, nMerHeatTransBinsUsed
               if (binVariable(iCell) .lt. binBoundaryMerHeatTrans(iBin+1)) then

                  do k = 1, kMax

                     ! Compute divergence of huT, i.e. layerThicknessEdge * normalTransportVelocity * temperature, at an edge
                     ! for meridional heat transport.  Here we use a centered difference to compute the temperature at
                     ! the edge, which is an approximation to the actual edge temperature used in the horizontal
                     ! advection scheme (for example, FCT).  We expect that the error in this approximation is small.
                     ! Here we do not divide by the area, as one normally does in a divergence calculation, so that
                     ! div_huT is weighted by area here.
                     iField = 1
                     div_huT = 0.0_RKIND
                     do i = 1, nEdgesOnCell(iCell)
                        iEdge = edgesOnCell(i, iCell)
                        div_huT = div_huT - layerThicknessEdge(k, iEdge) * normalTransportVelocity(k, iEdge) &
                             * 0.5_RKIND * (activeTracers(indexTemperature,k,cellsOnEdge(1,iEdge)) &
                             + activeTracers(indexTemperature,k,cellsOnEdge(2,iEdge))) &
                             * edgeSignOnCell(i, iCell) * dvEdge(iEdge)
                     end do
                     sumMerHeatTrans(iField,k,iBin) = sumMerHeatTrans(iField,k,iBin) + div_huT

                     !!!!! region-specific MHT
                     if(additionalRegion /= '') then
                        do i = 1, regionsInAddGroup
                           curRegion = regionsInGroup(i, regionGroupNumber)
                           sumMerHeatTransRegion(iField,k,iBin,curRegion) = sumMerHeatTrans(iField,k,iBin) * &
                                                                    regionCellMasks(curRegion, iCell)
                        end do
                        !!!!! end region-specific MHT
                     end if

                  end do
                  exit
               endif
            end do
         end do

         block => block % next
      end do

      ! mpi summation over all processors
      ! Note the input and output arrays are of the same dimension, so summation is
      ! over the domain decompositon (by processor), not over an array index.
      call mpas_dmpar_sum_real_array(dminfo, nVertLevels*nMerHeatTransBinsUsed*nMerHeatTransVariables, &
                                     sumMerHeatTrans, totalSumMerHeatTrans)

      !!!! Region version
      if(additionalRegion /= '') then
         do i = 1, regionsInAddGroup
            curRegion = regionsInGroup(i, regionGroupNumber)
            call mpas_dmpar_sum_real_array(dminfo, nVertLevels*nMerHeatTransBinsUsed*nMerHeatTransVariables, &
                                           sumMerHeatTransRegion(:,:,:,curRegion), &
                                        totalSumMerHeatTransRegion(:,:,:,curRegion))
         end do
      end if

      ! Even though these variables do not include an index that is decomposed amongst
      ! domain partitions, we assign them within a block loop so that all blocks have the
      ! correct values for writing output.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_dimension(block % dimensions, 'nMerHeatTransBins', nMerHeatTransBins)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

         call mpas_pool_get_subpool(block % structs, 'meridionalHeatTransportAM', meridionalHeatTransportAMPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)

         call mpas_pool_get_array(meridionalHeatTransportAMPool, 'meridionalHeatTransportLat', merHeatTransLat)
         call mpas_pool_get_array(meridionalHeatTransportAMPool, 'meridionalHeatTransportLatZ', merHeatTransLatZ)

         if(additionalRegion /= '') then
            call mpas_pool_get_array(meridionalHeatTransportAMPool, 'merHeatTransLatRegion', merHeatTransLatRegion)
            call mpas_pool_get_array(meridionalHeatTransportAMPool, 'merHeatTransLatZRegion', merHeatTransLatZRegion)
         end if

         do iBin = 1, nMerHeatTransBinsUsed
            do k = 1, nVertLevels

               ! MHT = sum ( div(huT) A ) * rho c_p, in PW
               ! where the sum is over each latitude bin
               ! Here we simply multiply by (rho c_p) and convert to PW:
               iField = 1
               mht_meridional_integral(k,iBin) = totalSumMerHeatTrans(iField,k,iBin)*rho_sw*cp_sw*1.0e-15_RKIND
               ! compute MHT per bin for regions
               if(additionalRegion /= '') then
                  do i = 1, regionsInAddGroup
                     curRegion = regionsInGroup(i, regionGroupNumber)
                     mht_meridional_integral_region(k,iBin,curRegion) = totalSumMerHeatTransRegion(iField,k,iBin,curRegion) &
                                * rho_sw*cp_sw*1.0e-15_RKIND
                  end do
               end if
            end do
         end do

         ! Compute integral of ( sum ( div(huT) A ) * rho c_p ) from southernmost latitude to bin boundary.
         ! Note that mht_meridional_integral is indexed by bin, spanning 1:nMerHeatTransBinsUsed, while
         ! merHeatTransLatZ (second dimension) is indexed by bin boundary, spanning 1:nMerHeatTransBinsUsed+1
         merHeatTransLatZ(:,1) = 0.0_RKIND
         do iBin = 2, nMerHeatTransBinsUsed+1
            merHeatTransLatZ(:,iBin) = merHeatTransLatZ(:,iBin-1) + mht_meridional_integral(:,iBin-1)
            ! integrate MHT for regions
            if(additionalRegion /= '') then
               do i = 1, regionsInAddGroup
                  curRegion = regionsInGroup(i, regionGroupNumber)
                  merHeatTransLatZRegion(:,iBin,curRegion) = merHeatTransLatZRegion(:,iBin-1,curRegion) +&
                            mht_meridional_integral_region(:,iBin-1,curRegion)
               end do
            end if
         end do

         ! merHeatTransLatZ is a function of depth.  Sum in vertical to get
         ! merHeatTransLat, a single value for each latitude bin boundary.
         ! merHeatTransLat is indexed by bin boundary, spanning 1:nMerHeatTransBinsUsed+1
         do iBin = 1, nMerHeatTransBinsUsed+1
            merHeatTransLat(iBin) = sum(merHeatTransLatZ(:,iBin))
            ! sum up MHT for regions
            if(additionalRegion /= '') then
               do i = 1, regionsInAddGroup
                  curRegion = regionsInGroup(i, regionGroupNumber)
                  merHeatTransLatRegion(iBin,curRegion) = sum(merHeatTransLatZRegion(:,iBin,curRegion))
               end do
            end if
         end do

         block => block % next
      end do

      call mpas_dmpar_sum_real_array(dminfo, nVertLevels*nMerHeatTransBinsUsed*nMerHeatTransVariables, &
                                     sumMerHeatTrans, totalSumMerHeatTrans)

      deallocate(sumMerHeatTrans)
      deallocate(totalSumMerHeatTrans)
      deallocate(mht_meridional_integral)

      !!!! region clean-up
      if(additionalRegion /= '') then
         deallocate(totalSumMerHeatTransRegion)
         deallocate(sumMerHeatTransRegion)
      end if
   end subroutine ocn_compute_meridional_heat_transport!}}}

!***********************************************************************
!
!  routine ocn_restart_meridional_heat_transport
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2014
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_meridional_heat_transport(domain, err)!{{{

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

   end subroutine ocn_restart_meridional_heat_transport!}}}

!***********************************************************************
!
!  routine ocn_finalize_meridional_heat_transport
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2014
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_meridional_heat_transport(domain, err)!{{{

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

   end subroutine ocn_finalize_meridional_heat_transport!}}}

end module ocn_meridional_heat_transport

! vim: foldmethod=marker
