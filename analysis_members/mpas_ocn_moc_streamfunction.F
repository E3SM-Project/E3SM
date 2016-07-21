! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_TEM_PLATE
!
!> \brief MPAS ocean analysis mode member: TEM_PLATE
!> \author FILL_IN_AUTHOR
!> \date   FILL_IN_DATE
!> \details
!>  MPAS ocean analysis mode member: TEM_PLATE
!>  In order to add a new analysis member, do the following:
!>  1. Copy these to your new analysis member name:
!>     cp mpas_ocn_TEMPLATE.F mpas_ocn_your_new_name.F
!>     cp Registry_ocn_TEMPLATE.xml Registry_ocn_your_new_name.xml
!>
!>  2. In those two new files, replace the following text:
!>     temPlate, TEM_PLATE, FILL_IN_AUTHOR, FILL_IN_DATE
!>     Typically temPlate uses camel case (variable names), like yourNewName,
!>     while TEM_PLATE uses underscores (subroutine names), like your_new_name.
!>     note: do not replace 'filename_template' in Registry_ocn_yourNewName.xml
!>
!>  3. Add a #include line for your registry to
!>     Registry_analysis_members.xml
!>
!>  4. In mpas_ocn_analysis_driver.F, add a use statement for your new analysis member.
!>     In addition, add lines for your analysis member, and replace TEM_PLATE
!>     and temPlate as described in step 2. There should be 5 places that need additions:
!>      - Adding the analysis member name to the analysis member list
!>      - Adding an init if test can subroutine call
!>      - Adding a compute if test can subroutine call
!>      - Adding a restart if test can subroutine call
!>      - Adding a finalize if test can subroutine call
!>
!>  5. In src/core_ocean/analysis_members/Makefile, add your
!>     new analysis member to the list of members. See another analysis member
!>     in that file for an example.
!>     NOTE: If your analysis member depends on other files, add a dependency
!>           line for the member and list them there. See okubo weiss for an example.
!>
!-----------------------------------------------------------------------

module ocn_moc_streamfunction

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

   public :: ocn_init_moc_streamfunction, &
             ocn_compute_moc_streamfunction, &
             ocn_restart_moc_streamfunction, &
             ocn_finalize_moc_streamfunction

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   
   integer :: nMocStreamfunctionBinsUsed

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_moc_streamfunction
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Nils H. Feige
!> \date    2016-04-08
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_moc_streamfunction(domain, err)!{{{

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
      type (mpas_pool_type), pointer :: mocStreamfunctionAMPool
      type (mpas_pool_type), pointer :: meshPool

      integer ::  iBin
     
      real (kind=RKIND) :: binWidth
      ! These are array size 1 because mpas_dmpar_min_real_array calls require arrays.
      real (kind=RKIND), dimension(1) :: minBin, maxBin, minBinDomain, maxBinDomain
      ! the variable used to discriminate cells into Bins (either the y-value or the latitude)
      real (kind=RKIND), dimension(:), pointer :: binVariable, binBoundaryMocStreamfunction

      !number of latitude bins specified in the config
      integer, pointer :: config_AM_mocStreamfunction_num_bins
      !smallest and highest latitude specified in the config
      real (kind=RKIND), pointer :: config_AM_mocStreamfunction_min_bin, config_AM_mocStreamfunction_max_bin

      !determines if the simulation was run on a sphere or on a plane
      logical, pointer :: on_a_sphere
      
      !!!! REGION STUFF
      !! region moc calculation variables
      integer :: currentRegion, i, iCell
                         
      !! region arrays/variables
      character (len=STRKIND), dimension(:), pointer :: regionGroupNames
      integer, dimension(:, :), pointer :: regionCellMasks, regionsInGroup
      integer, dimension(:), pointer ::  nRegionsInGroup
      integer, pointer :: nRegions, nRegionGroups, maxRegionsInGroup, nCellsSolve
      real (kind=RKIND), dimension(:,:), pointer :: minMaxLatRegion
      real (kind=RKIND), dimension(:)  , pointer :: minLatRegionLocal, maxLatRegionLocal, tminLatRegionLocal, tmaxLatRegionLocal
      character (len=STRKIND), pointer :: additionalRegion
      
      !! region preliminary variables
      integer :: regionGroupNumber, regionsInAddGroup
      
      !!region pool
      type (mpas_pool_type), pointer :: regionPool
      
      !! region dimensions
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegions', nRegions)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegionGroups', nRegionGroups)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'maxRegionsInGroup', maxRegionsInGroup)
      
      !! region config for moc
      call mpas_pool_get_config(domain % configs, 'config_AM_mocStreamfunction_region_group', &
                                additionalRegion)
      
      !! get region values
      call mpas_pool_get_subpool(domain % blocklist % structs, 'regions', regionPool)
      call mpas_pool_get_array(regionPool, 'regionsInGroup', regionsInGroup)
      call mpas_pool_get_array(regionPool, 'nRegionsInGroup', nRegionsInGroup)
      call mpas_pool_get_array(regionPool, 'regionGroupNames', regionGroupNames)
      
      !!! region preliminaries
      do i = 1, nRegionGroups
         if (regionGroupNames(i) .eq. additionalRegion) then
            regionGroupNumber = i
         end if
      end do
      
      regionsInAddGroup = nRegionsInGroup(regionGroupNumber)
      !!!! END REGION STUFF
      
      allocate(minLatRegionLocal(maxRegionsInGroup))
      allocate(maxLatRegionLocal(maxRegionsInGroup))
      allocate(tminLatRegionLocal(maxRegionsInGroup))
      allocate(tmaxLatRegionLocal(maxRegionsInGroup))
      
      dminfo = domain % dminfo

      err = 0

      minBin =  1.0e34_RKIND
      maxBin = -1.0e34_RKIND

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mocStreamfunctionAM', mocStreamfunctionAMPool)

      call mpas_pool_get_config(domain % configs, 'config_AM_mocStreamfunction_num_bins', &
                                config_AM_mocStreamfunction_num_bins)
      call mpas_pool_get_config(domain % configs, 'config_AM_mocStreamfunction_min_bin', &
                                config_AM_mocStreamfunction_min_bin)
      call mpas_pool_get_config(domain % configs, 'config_AM_mocStreamfunction_max_bin', &
                                config_AM_mocStreamfunction_max_bin)
      
      call mpas_pool_get_array(mocStreamfunctionAMPool, 'minMaxLatRegion', minMaxLatRegion)
      
      minLatRegionLocal(:) =  4.0_RKIND
      maxLatRegionLocal(:) = -4.0_RKIND                      
                             
      nMocStreamfunctionBinsUsed = config_AM_mocStreamfunction_num_bins

      call mpas_pool_get_array(mocStreamfunctionAMPool, 'binBoundaryMocStreamfunction', binBoundaryMocStreamfunction)

      ! Find min and max values of binning variable. For the whole domain as well as for each region
      ! in the current region group.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         
         call mpas_pool_get_array(regionPool, 'regionCellMasks', regionCellMasks)
         call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)	!!!! this could probably NOT be pulled before the loop.

         ! Bin by latitude on a sphere, by yCell otherwise.
         if (on_a_sphere) then
            call mpas_pool_get_array(meshPool, 'latCell', binVariable)
         else
            call mpas_pool_get_array(meshPool, 'yCell', binVariable)
         end if

         minBin = min(minBin, minval(binVariable) )
         maxBin = max(maxBin, maxval(binVariable) )
         
         do i = 1, regionsInAddGroup
            currentRegion = regionsInGroup(i, regionGroupNumber)
            do iCell = 1, nCellsSolve
               if (regionCellMasks(currentRegion, iCell) .eq. 1) then
                  minLatRegionLocal(i) = min(minLatRegionLocal(i), binVariable(iCell))
                  maxLatRegionLocal(i) = max(maxLatRegionLocal(i), binVariable(iCell))
               end if
            end do
         end do

         block => block % next
      end do

      call mpas_dmpar_min_real_array(dminfo, 1, minBin, minBinDomain)
      call mpas_dmpar_max_real_array(dminfo, 1, maxBin, maxBinDomain)
      
      call mpas_dmpar_min_real_array(dminfo, maxRegionsInGroup, minLatRegionLocal(:), tminLatRegionLocal(:))
      call mpas_dmpar_max_real_array(dminfo, maxRegionsInGroup, maxLatRegionLocal(:), tmaxLatRegionLocal(:))
      minMaxLatRegion(1, :) = tminLatRegionLocal(:)
      minMaxLatRegion(2, :) = tmaxLatRegionLocal(:)
      
      deallocate(minLatRegionLocal)
      deallocate(maxLatRegionLocal)
      deallocate(tminLatRegionLocal)
      deallocate(tmaxLatRegionLocal)

      !print *, 'mins:', minMaxLatRegion(1,:)
      !print *, 'maxs:', minMaxLatRegion(2,:)

      ! Set up bins.
      binBoundaryMocStreamfunction = -1.0e34_RKIND

      ! Change min and max bin bounds to configuration settings, if applicable.
      if (config_AM_mocStreamfunction_min_bin > -1.0e33_RKIND) then
         minBinDomain(1) = config_AM_mocStreamfunction_min_bin
      else
         ! use measured min value, but decrease slightly to include least value.
         minBinDomain(1) = minBinDomain(1) - 1.0e-10_RKIND * abs(minBinDomain(1))
      end if

      if (config_AM_mocStreamfunction_max_bin > -1.0e33_RKIND) then
         maxBinDomain(1) = config_AM_mocStreamfunction_max_bin
      else
         ! use measured max value, but increase slightly to include max value.
         maxBinDomain(1) = maxBinDomain(1) + 1.0e-10_RKIND * abs(maxBinDomain(1))
      end if

      binBoundaryMocStreamfunction(1) = minBinDomain(1)
      binWidth = (maxBinDomain(1) - minBinDomain(1)) / nMocStreamfunctionBinsUsed

      ! Use the same bin boundaries for the regions and the global MOC.
      do iBin = 2, nMocStreamfunctionBinsUsed + 1
         binBoundaryMocStreamfunction(iBin) = binBoundaryMocStreamfunction(iBin-1) + binWidth
      end do
      !!! integrate next line into loop by simply going from 2 to nMocStreamfunctionBinsUsed+1.
      !binBoundaryMocStreamfunction(nMocStreamfunctionBinsUsed+1) = binBoundaryMocStreamfunction(nMocStreamfunctionBinsUsed) + binWidth

   end subroutine ocn_init_moc_streamfunction!}}}

!***********************************************************************
!
!  routine ocn_compute_moc_streamfunction
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Nils H. Feige
!> \date    2016-04-08
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_moc_streamfunction(domain, timeLevel, err)!{{{
      implicit none
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

      type (mpas_pool_type), pointer :: mocStreamfunctionAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: regionPool
      type (mpas_pool_type), pointer :: transectPool
     
      integer, pointer :: nVertLevels, nCellsSolve
      integer :: iCell, iBin, genericCounter, k
      real (kind=RKIND) :: binWidth
      real (kind=RKIND), dimension(:,:), pointer :: mocStreamValLatAndDepthLocal
	  real (kind=RKIND), dimension(:), pointer :: latCell
      real (kind=RKIND), dimension(:), pointer ::  areaCell, binBoundaryMocStreamfunction
      real (kind=RKIND), dimension(:,:), pointer :: mocStreamvalLatAndDepth, mocStreamValLatAndDepthTotal
      real (kind=RKIND), dimension(:,:), pointer :: vertVelocityTop
      real (kind=RKIND), dimension(:,:), pointer :: sumVertBinVelocity
      character (len=STRKIND), pointer :: verticalVelocityArrayName, normalVelocityArrayName
      
      !!!! TRANSECT VARIABLES !!!!
      integer, pointer :: nEdgesSolve, num_tracers
      integer :: iEdge, iTransect, c1, c2, currentTransect
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: transectEdgeMasks, transectEdgeMaskSigns, cellsOnEdge, &
            transectsInGroup
 
      real (kind=RKIND) :: m3ps_to_Sv
      real (kind=RKIND), dimension(:), pointer ::  dvEdge, refLayerThickness
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, normalVelocity
      real (kind=RKIND), dimension(:,:), allocatable ::  sumTransport, totalSumTransport
      character (len=STRKIND), dimension(:), pointer :: transectNames, transectGroupNames
      integer, dimension(:), pointer ::  nTransectsInGroup
      integer, pointer :: nTransects, nTransectGroups, maxTransectsInGroup
      character (len=STRKIND), pointer :: additionalTransect
      integer :: transectGroupNumber, transectsInAddGroup
      !!!! END TRANSECT VARIABLES !!!!
      
      !!!! REGION VARIABLES
      real (kind=RKIND) :: maskFactor
      integer :: currentRegion, i
      real (kind=RKIND), dimension(:,:,:), pointer :: mocStreamValLatAndDepthRegionLocal, &
                         mocStreamvalLatAndDepthRegion, mocStreamValLatAndDepthRegionTotal, &
                         sumVertBinVelocityRegion
      integer, dimension(:, :), pointer :: regionCellMasks, regionVertexMasks, regionsInGroup
      character (len=STRKIND), dimension(:), pointer :: regionNames, regionGroupNames
      integer, dimension(:), pointer ::  nRegionsInGroup
      integer, pointer :: nRegions, nRegionGroups, maxRegionsInGroup
      character (len=STRKIND), pointer :: additionalRegion
      integer :: regionGroupNumber, regionsInAddGroup

      !!!! END REGION VARIABLES
      
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)
            
      !!!! REGION INITIALIZATION
      !! region dimensions
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegions', nRegions)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nRegionGroups', nRegionGroups)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'maxRegionsInGroup', maxRegionsInGroup)
      
      !! region config for moc
      call mpas_pool_get_config(domain % configs, 'config_AM_mocStreamfunction_region_group', &
                                additionalRegion)
      
      !! get region values
      call mpas_pool_get_subpool(domain % blocklist % structs, 'regions', regionPool)
      call mpas_pool_get_array(regionPool, 'regionsInGroup', regionsInGroup)
      call mpas_pool_get_array(regionPool, 'nRegionsInGroup', nRegionsInGroup)
      call mpas_pool_get_array(regionPool, 'regionNames', regionNames)
      call mpas_pool_get_array(regionPool, 'regionGroupNames', regionGroupNames)
      
         !print *, 'nRegions', nRegions
         !print *, 'nRegionGroups', nRegionGroups
         !print *, 'maxRegionsInGroup', maxRegionsInGroup
         !print *, 'regionsInGroup', regionsInGroup
         !print *, 'regionNames', regionNames
         !print *, 'regionGroupNames', regionGroupNames
      
      !!! region preliminaries
         !print *, 'additional Region: ', additionalRegion
      do i = 1, nRegionGroups
         if (regionGroupNames(i) .eq. additionalRegion) then
            regionGroupNumber = i
            !print *, 'found region with the same name:', regionGroupNames(i)
         end if
      end do
      
      regionsInAddGroup = nRegionsInGroup(regionGroupNumber)
         !print *, 'the region group has', regionsInAddGroup, ' regions'
      
         !print *, 'the numbers of the regions in the group are:'
         !do i = 1, regionsInAddGroup
            !print *, regionsInGroup(i, regionGroupNumber)
         !end do
      
         !print *, 'the names of the regions are:'
         !do i = 1, regionsInAddGroup
            !print *, regionNames(regionsInGroup(i, regionGroupNumber))
         !end do
      
      !! allocate regional moc calculation arrays
      allocate(mocStreamValLatAndDepthRegionLocal(nMocStreamfunctionBinsUsed + 1, nVertLevels, maxRegionsInGroup))
      allocate(sumVertBinVelocityRegion(nMocStreamfunctionBinsUsed + 1, nVertLevels, maxRegionsInGroup))
      allocate(mocStreamValLatAndDepthRegionTotal(nMocStreamfunctionBinsUsed + 1, nVertLevels, maxRegionsInGroup))
      
      mocStreamValLatAndDepthRegionLocal = 0.0_RKIND
      mocStreamValLatAndDepthRegionTotal = 0.0_RKIND
      sumVertBinVelocityRegion = 0.0_RKIND
      !!!! END REGION INITIALIZATION
      
      !!!! TRANSECT INITIALIZATION
      !! transect dimensions
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nTransects', nTransects)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nTransectGroups', nTransectGroups)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'maxTransectsInGroup', maxTransectsInGroup)
      
      !!!! Transect config for MOC
      call mpas_pool_get_config(domain % configs, 'config_AM_mocStreamfunction_transect_group', &
                                additionalTransect)
      print *, additionalTransect                          
      !! get transect values
      call mpas_pool_get_subpool(domain % blocklist % structs, 'transects', transectPool)
      call mpas_pool_get_array(transectPool, 'transectsInGroup', transectsInGroup)
      call mpas_pool_get_array(transectPool, 'nTransectsInGroup', nTransectsInGroup)
      call mpas_pool_get_array(transectPool, 'transectNames', transectNames)
      call mpas_pool_get_array(transectPool, 'transectGroupNames', transectGroupNames)
      do i = 1, nTransectGroups
         if (transectGroupNames(i) .eq. additionalTransect) then
            transectGroupNumber = i
         end if
      end do
      
      transectsInAddGroup = nTransectsInGroup(transectGroupNumber)
      
      allocate(sumTransport(nVertLevels,maxTransectsInGroup))
      allocate(totalSumTransport(nVertLevels,maxTransectsInGroup))
      
      m3ps_to_Sv = 1e-6
      !!!! END TRANSECT INITIALIZATION
      
      if (transectsInAddGroup .ne. regionsInAddGroup) then
         !!!FAIL HORRIBLY!!!
         print *, 'transectsInGroup count does not match regionsInGroup count!', transectsInAddGroup, regionsInAddGroup
      end if
      
      err = 0
      
      dminfo = domain % dminfo
      
      allocate(mocStreamValLatAndDepthLocal(nMocStreamfunctionBinsUsed + 1, nVertLevels))
      allocate(sumVertBinVelocity(nMocStreamfunctionBinsUsed + 1, nVertLevels))
      allocate(mocStreamValLatAndDepthTotal(nMocStreamfunctionBinsUsed + 1, nVertLevels))
      
      mocStreamValLatAndDepthLocal = 0.0_RKIND
      mocStreamValLatAndDepthTotal = 0.0_RKIND
      sumVertBinVelocity = 0.0_RKIND
      
      call mpas_pool_get_config(domain % configs, 'config_AM_mocStreamfunction_vertical_velocity_value', &
                                verticalVelocityArrayName)

      call mpas_pool_get_config(domain % configs, 'config_AM_mocStreamfunction_normal_velocity_value', &
                                normalVelocityArrayName)
      
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'mocStreamfunctionAM', mocStreamfunctionAMPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         call mpas_pool_get_array(mocStreamfunctionAMPool, 'binBoundaryMocStreamfunction', binBoundaryMocStreamfunction)

         binWidth = (binBoundaryMocStreamfunction(nMocStreamfunctionBinsUsed + 1) - binBoundaryMocStreamfunction(1)) & 
         / nMocStreamfunctionBinsUsed
         !!! NOT SURE ABOUT THIS BUT POSSIBLY BETTER
         !!! binWidth = binBoundaryMocStreamfunction(2) - binBoundaryMocStreamfunction(1)
         
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_array(diagnosticsPool, verticalVelocityArrayName, vertVelocityTop)
         call mpas_pool_get_array(meshPool, 'latCell', latCell)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         
         call mpas_pool_get_array(regionPool, 'regionCellMasks', regionCellMasks)
         
         !!!! TRANSECT DOMAINSPLIT VARIABLES
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_array(transectPool,'transectEdgeMaskSigns',transectEdgeMaskSigns)
         call mpas_pool_get_array(transectPool,'transectEdgeMasks',transectEdgeMasks)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_array(statePool, normalVelocityArrayName, normalVelocity)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness)
         
         
         !!!! TRANSECT CALCULATION
         !!!! NOTE: iTransect == iRegion !!!! This might have to change!
         !!!! Define transect group, find current transect group number, find
         !!!! nTransects in that group, find transects in group. THEN for each
         !!!! transect in the group, there should be a region in the region
         !!!! group so that the region number equals the transect number

         sumTransport = 0.0_RKIND
         do iTransect = 1,transectsInAddGroup
            currentTransect = transectsInGroup(iTransect, transectGroupNumber)
            do iEdge = 1,nEdgesSolve
               c1 = cellsOnEdge(1,iEdge)
               c2 = cellsOnEdge(2,iEdge)
               do k = 1, maxLevelEdgeTop(iEdge)
                  sumTransport(k,iTransect) = sumTransport(k,iTransect) + &
                     transectEdgeMaskSigns(currentTransect,iEdge) &
                     * transectEdgeMasks(currentTransect, iEdge) &
                     * normalVelocity(k,iEdge)*dvEdge(iEdge) &
                     * 0.5_RKIND*(layerThickness(k,c1) + layerThickness(k,c2))
               end do
            end do
            do k = 2, nVertLevels
               mocStreamValLatAndDepthRegionLocal(1, k, iTransect) = & 
                  mocStreamValLatAndDepthRegionLocal(1, k - 1, iTransect) & 
                  + sumTransport(k - 1, iTransect)
            end do
         end do
 
         !!!! END TRANSECT CALCULATION
         
         do iCell = 1,nCellsSolve
            !!!!do iBin = 1, nMocStreamfunctionBinsUsed
            !!!! useless inner loop. calculate bin number directly instead of going through this loop up to nMocStreamfunctionBinsUsed times
            !!!! BETTER: iBin = int(latCell(iCell) / binWidth) +1 !!! or plus 2, if the values of the first bin should be stored in the second column of the array
            !!!!if (latCell(iCell) .gt. binBoundaryMocStreamfunction(iBin) .and. latCell(iCell) .lt. binBoundaryMocStreamfunction(iBin + 1)) then
                
                iBin = MAX(int((latCell(iCell) - binBoundaryMocStreamfunction(1)) / binWidth), 2)
                
                !!!! TESTING CODE. THIS CODE ALTERS THE VERTICAL VELOCITY TO PRODUCE A PREDICTABLE 
                ! PATTERN IF THE CALCULATION IS CORRECT. TESTS SHOW THIS PREDICTABLE PATTERN SO THIS 
                ! IS NOT USED ANYMORE.  
                !
                !  if (latCell(iCell) .lt. -0.785398_RKIND .or. latCell(iCell) .gt. 0.785398_RKIND) then
                !     vertVelocityTop(:, iCell) = 0.0_RKIND
                !  else
                !     vertVelocityTop(1:25, iCell) = 0.0_RKIND
                !     vertVelocityTop(76:100, iCell) = 0.0_RKIND
                !        if (latCell(iCell) .lt. 0.0_RKIND) then
                !           do k = 26, 75
                !              vertVelocityTop(k, iCell) = 1.0_RKIND + k - 26
                !           end do
                !        else
                !           do k = 26, 75
                !              vertVelocityTop(k, iCell) = -1.0_RKIND - k + 26
                !           end do
                !        end if
                !  end if
                !!!! END TESTING CODE
                
                  !!!!! REGION STUFF

                  do i = 1, regionsInAddGroup
                     currentRegion = regionsInGroup(i, regionGroupNumber)
                     sumVertBinVelocityRegion(iBin, :, i) = sumVertBinVelocityRegion(iBin, :, i) + (vertVelocityTop(:, iCell) * &
                              areaCell(iCell) * regionCellMasks(currentRegion, iCell))
                  end do

                  !!!!! END REGION STUFF
                  sumVertBinVelocity(iBin, :) = sumVertBinVelocity(iBin, :) + (vertVelocityTop(:, iCell) * areaCell(iCell))
         !!!!     end if
         !!!!   end do
         end do
         
         do iBin = 2, nMocStreamfunctionBinsUsed + 1
            mocStreamValLatAndDepthLocal(iBin, :) = mocStreamValLatAndDepthLocal(iBin-1, :) &
               + sumVertBinVelocity(iBin, :)
            !do i = 1, regionsInAddGroup !!! should be commented in again at some point 
            !!!! CURRENTLY ONLY THE FIRST REGION WILL BE CALCULATED 
            !!! THIS IS CONSISTENT WITH TESTCASE WITH ONLY ONE TRANSECT
            mocStreamValLatAndDepthRegionLocal(iBin, :, :) = mocStreamValLatAndDepthRegionLocal(iBin-1, :, :) &
               + sumVertBinVelocityRegion(iBin, :, :)
            !end do
         end do  

         block => block % next
     end do
     
     call mpas_dmpar_sum_real_array(dminfo, nVertLevels * (nMocStreamfunctionBinsUsed + 1), mocStreamValLatAndDepthLocal, &
         mocStreamvalLatAndDepthTotal)
         
     call mpas_dmpar_sum_real_array(dminfo, nVertLevels * (nMocStreamfunctionBinsUsed + 1) * maxRegionsInGroup, mocStreamValLatAndDepthRegionLocal, &
         mocStreamvalLatAndDepthRegionTotal)
         
     call mpas_pool_get_subpool(domain % blocklist % structs, 'mocStreamfunctionAM', mocStreamfunctionAMPool)
     call mpas_pool_get_array(mocStreamfunctionAMPool, 'mocStreamvalLatAndDepth', mocStreamvalLatAndDepth)
     mocStreamvalLatAndDepth = mocStreamvalLatAndDepthTotal * m3ps_to_Sv
     
     call mpas_pool_get_array(mocStreamfunctionAMPool, 'mocStreamvalLatAndDepthRegion', mocStreamvalLatAndDepthRegion)
     mocStreamvalLatAndDepthRegion = mocStreamvalLatAndDepthRegionTotal * m3ps_to_Sv

     deallocate(mocStreamvalLatAndDepthTotal)
     deallocate(mocStreamvalLatAndDepthLocal)
     deallocate(sumVertBinVelocity)
     
     deallocate(mocStreamvalLatAndDepthRegionTotal)
     deallocate(mocStreamvalLatAndDepthRegionLocal)
     deallocate(sumVertBinVelocityRegion)
     
     !!!! TRANSECT CELANUP
     deallocate(sumTransport)
     deallocate(totalSumTransport)

   end subroutine ocn_compute_moc_streamfunction!}}}

!***********************************************************************
!
!  routine ocn_restart_moc_streamfunction
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Nils H. Feige
!> \date    2016-04-08
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_moc_streamfunction(domain, err)!{{{

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

   end subroutine ocn_restart_moc_streamfunction!}}}

!***********************************************************************
!
!  routine ocn_finalize_moc_streamfunction
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Nils H. Feige
!> \date    2016-04-08
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_moc_streamfunction(domain, err)!{{{

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

   end subroutine ocn_finalize_moc_streamfunction!}}}

end module ocn_moc_streamfunction
