! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_global_ocean
!
!> \brief MPAS ocean initialize case -- Global Ocean
!> \author Doug Jacobsen
!> \date   03/04/2014
!> \details
!>  This module contains the routines for initializing the 
!>  the global ocean test case
!
!-----------------------------------------------------------------------

module ocn_init_global_ocean

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_io
   use mpas_io_streams
   use mpas_dmpar

   use ocn_init_cell_markers

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

   public :: ocn_init_setup_global_ocean, &
             ocn_init_validate_global_ocean

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: nDepth
   integer :: nLatTracer, nLonTracer
   integer :: nLatWind, nLonWind
   integer :: nLatTopo, nLonTopo
   type (field1DReal) :: depthIC
   type (field1DReal) :: windLat, windLon
   type (field1DReal) :: topoLat, topoLon
   type (field1DReal) :: tracerLat, tracerLon
   type (field2DReal) :: topoIC, zonalWindIC, meridionalWindIC
   type (field3DReal) :: temperatureIC, salinityIC

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean
!
!> \brief   Setup for global ocean test case
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details 
!>  This routine sets up the initial conditions for the global ocean test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_global_ocean(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      type (mpas_pool_type), pointer :: meshPool
      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      logical, pointer :: config_global_ocean_cull_inland_seas

      logical, pointer :: on_a_sphere

      iErr = 0

      call mpas_pool_get_config(domain % configs, 'config_init_configuration', config_init_configuration)

      if (trim(config_init_configuration) /= "global_ocean") return

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      if ( .not. on_a_sphere ) call mpas_dmpar_global_abort('ERROR: The global ocean configuration can only be applied to a spherical mesh. Exiting...')

      call mpas_pool_get_config(domain % configs, 'config_global_ocean_cull_inland_seas', config_global_ocean_cull_inland_seas)

      write(stderrUnit,*) 'Reading depth levels.'
      call ocn_init_setup_global_ocean_read_depth_levels(domain, iErr)

      write(stderrUnit,*) 'Reading topography data.'
      call ocn_init_setup_global_ocean_read_topo(domain, iErr)
      write(stderrUnit,*) 'Interpolating topography data.'
      call ocn_init_setup_global_ocean_interpolate_topo(domain, iErr)
      write(stderrUnit,*) 'Cleaning up topography IC fields'
      call ocn_init_global_ocean_destroy_topo_fields()

      if (config_global_ocean_cull_inland_seas) then
         write(stderrUnit,*) 'Removing inland seas.'
         call ocn_init_setup_global_ocean_cull_inland_seas(domain, iErr)
      end if


      write(stderrUnit,*) 'Reading temperature IC.'
      call ocn_init_setup_global_ocean_read_temperature(domain, iErr)
      write(stderrUnit,*) 'Reading salinity IC.'
      call ocn_init_setup_global_ocean_read_salinity(domain, iErr)
      write(stderrUnit,*) 'Reading Lat/Lon tracer coordinates'
      call ocn_init_setup_global_ocean_read_tracer_lat_lon(domain, iErr)
      write(stderrUnit,*) 'Interpolating tracers'
      call ocn_init_setup_global_ocean_interpolate_tracers(domain, iErr)
      write(stderrUnit,*) 'Cleaning up tracer IC fields'
      call ocn_init_global_ocean_destroy_tracer_fields()

      write(stderrUnit,*) 'Reading windstress IC.'
      call ocn_init_setup_global_ocean_read_windstress(domain, iErr)
      write(stderrUnit,*) 'Interpolating windstress.'
      call ocn_init_setup_global_ocean_interpolate_windstress(domain, iErr)
      write(stderrUnit,*) 'Destroying windstress fields'
      call ocn_init_global_ocean_destroy_windstress_fields()

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_global_ocean!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_topo
!
!> \brief   Read the topography IC file
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details 
!>  This routine reads the topography IC file, including latitude and longitude
!>   information for topography data.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_topo(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (MPAS_Stream_type) :: topographyStream

       character (len=StrKIND), pointer :: config_global_ocean_topography_file, config_global_ocean_topography_lat_varname, &
                                           config_global_ocean_topography_nlat_dimname, config_global_ocean_topography_lon_varname, &
                                           config_global_ocean_topography_nlon_dimname, config_global_ocean_topography_varname

       logical, pointer :: config_global_ocean_topography_latlon_degrees

       integer :: iLat, iLon

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_file', config_global_ocean_topography_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_lat_varname', config_global_ocean_topography_lat_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_nlat_dimname', config_global_ocean_topography_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_lon_varname', config_global_ocean_topography_lon_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_nlon_dimname', config_global_ocean_topography_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_varname', config_global_ocean_topography_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_latlon_degrees', config_global_ocean_topography_latlon_degrees)

       ! Define stream for depth levels
       call MPAS_createStream(topographyStream, domain % iocontext, config_global_ocean_topography_file, MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)

       ! Setup topoLat, topoLon, and topoIC fields for stream to be read in
       topoLat % fieldName = trim(config_global_ocean_topography_lat_varname)
       topoLat % dimSizes(1) = nLatTopo
       topoLat % dimNames(1) = trim(config_global_ocean_topography_nlat_dimname)
       topoLat % isVarArray = .false.
       topoLat % isPersistent = .true.
       topoLat % isActive = .true.
       topoLat % hasTimeDimension = .false.
       topoLat % block => domain % blocklist
       allocate(topoLat % array(nLatTopo))

       topoLon % fieldName = trim(config_global_ocean_topography_lon_varname)
       topoLon % dimSizes(1) = nLonTopo
       topoLon % dimNames(1) = trim(config_global_ocean_topography_nlon_dimname)
       topoLon % isVarArray = .false.
       topoLon % isPersistent = .true.
       topoLon % isActive = .true.
       topoLon % hasTimeDimension = .false.
       topoLon % block => domain % blocklist
       allocate(topoLon % array(nLonTopo))

       topoIC % fieldName = trim(config_global_ocean_topography_varname)
       topoIC % dimSizes(1) = nLonTopo
       topoIC % dimSizes(2) = nLatTopo
       topoIC % dimNames(1) = trim(config_global_ocean_topography_nlon_dimname)
       topoIC % dimNames(2) = trim(config_global_ocean_topography_nlat_dimname)
       topoIC % isVarArray = .false.
       topoIC % isPersistent = .true.
       topoIC % isActive = .true.
       topoIC % hasTimeDimension = .false.
       topoIC % block => domain % blocklist
       allocate(topoIC % array(nLonTopo, nLatTopo))

       ! Add topoLat, topoLon, and topoIC fields to stream
       call MPAS_streamAddField(topographyStream, topoLat, iErr)
       call MPAS_streamAddField(topographyStream, topoLon, iErr)
       call MPAS_streamAddField(topographyStream, topoIC, iErr)

       ! Read stream
       call MPAS_readStream(topographyStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(topographyStream)

       if (config_global_ocean_topography_latlon_degrees) then
          topoLat % array(:) = topoLat % array(:) * pii / 180.0_RKIND 
          topoLon % array(:) = topoLon % array(:) * pii / 180.0_RKIND
       end if

       do iLon = 1, nLonTopo
          if (topoLon % array(iLon) < 0.0_RKIND) then
             topoLon % array(iLon) = 2.0_RKIND * pii + topoLon % array(iLon) 
          end if
       end do

    end subroutine ocn_init_setup_global_ocean_read_topo!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_interpolate_topo
!
!> \brief   Interpolate the topography IC to MPAS mesh
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details 
!>  This routine interpolates topography data to the MPAS mesh. Currently it
!>   uses a bilinear interpolation
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_interpolate_topo(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (mpas_pool_type), pointer :: meshPool, scratchPool, statePool, verticalMeshPool

       real (kind=RKIND) :: currentLat, currentLon
       real (kind=RKIND) :: dist, minDist
       real (kind=RKIND) :: alpha, beta, depthLat1, depthLat2, proposedDepth

       real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, bottomDepth, bottomDepthObserved, refBottomDepth
       real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness

       integer, pointer :: nCells, nCellsSolve, nVertLevels

       type (field1DInteger), pointer :: maxLevelCellField, smoothedLevelsField
       integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
       integer, dimension(:, :), pointer :: cellsOnCell

       integer :: latSearch, lonSearch, searchIdx
       integer :: iCell, coc, j, k, maxLevel
       integer :: minimum_levels

       logical, pointer :: config_global_ocean_smooth_topography
       real (kind=RKIND), pointer :: config_global_ocean_minimum_depth

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_minimum_depth', config_global_ocean_minimum_depth)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_smooth_topography', config_global_ocean_smooth_topography)

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
          call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

          call mpas_pool_get_array(meshPool, 'latCell', latCell)
          call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
          call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
          call mpas_pool_get_array(meshPool, 'bottomDepthObserved', bottomDepthObserved)
          call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

          do k = 1, nVertLevels
             if (refBottomDepth(k).gt.config_global_ocean_minimum_depth) then
                minimum_levels = k
                write (stdoutUnit,'(a,f8.2,2a,i5,a,f8.2,a)') 'config_global_ocean_minimum_depth=',config_global_ocean_minimum_depth,' m.  ', &
                   'Setting minimum layer index to ',minimum_levels, ' with a bottom depth of ', refBottomDepth(k), ' m.'
                exit
             end if
          end do

          do iCell = 1, nCells
             currentLat = latCell(iCell)
             currentLon = lonCell(iCell)

             lonSearch = 1
             minDist = 2.0_RKIND * pii
             do searchIdx = 1, nLonTopo 
                dist = abs(currentLon - topoLon % array(searchIdx))
                if (dist < minDist) then
                   minDist = dist
                   lonSearch = searchIdx
                end if
             end do

             latSearch = 1
             minDist = 2.0_RKIND * pii
             do searchIdx = 1, nLatTopo
                dist = abs(currentLat - topoLat % array(searchIdx))
                if (dist < minDist) then
                   minDist = dist
                   latSearch = searchIdx
                end if
             end do

             if (topoIC % array(lonSearch, latSearch) < 0.0_RKIND) then
                ! Record depth of the bottom of the ocean, before any alterations for modeling purposes.
                bottomDepthObserved(iCell) = abs(topoIC % array(lonSearch, latSearch))

                ! Enforce minimum depth
                bottomDepth(iCell) = max(bottomDepthObserved(iCell), refBottomDepth(minimum_levels))

                maxLevelCell(iCell) = -1
                do k = 1, nVertLevels
                   if (refBottomDepth(k) >= bottomDepth(iCell)) then
                      maxLevelCell(iCell) = k
                      exit
                   end if
                end do

                if (maxLevelCell(iCell) == -1) then
                   maxLevelCell(iCell) = nVertLevels
                   bottomDepth(iCell) = refBottomDepth( nVertLevels )
                end if

             else
                bottomDepth(iCell) = 0.0_RKIND
                maxLevelCell(iCell) = -1
             end if
          end do

          ! Smooth depth levels. Enforce different in maxLevelCell to only be a maximum
          ! of 1 vertical level between two neighboring cells.
          if (config_global_ocean_smooth_topography) then
             call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

             call mpas_pool_get_field(scratchPool, 'smoothedLevels', smoothedLevelsField)

             call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
             call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)

             call mpas_allocate_scratch_field(smoothedLevelsField, .true.)

             maxLevelCell(nCells+1) = -1
             smoothedLevelsField % array = maxLevelCell

             do iCell = 1, nCellsSolve
                maxLevel = 0
                do j = 1, nEdgesOnCell(iCell)
                   coc = cellsOnCell(j, iCell)
                   maxLevel = max(maxLevel, maxLevelCell(coc))
                end do

                if (maxLevel < maxLevelCell(iCell) ) then
                   smoothedLevelsField % array(iCell) = maxLevel + 1
                   bottomDepth(iCell) = refBottomDepth(maxLevel + 1)
                end if
             end do

             maxLevelCell(:) = smoothedLevelsField % array(:)

             call mpas_deallocate_scratch_field(smoothedLevelsField, .true.)
          end if

          block_ptr => block_ptr % next
       end do

       call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_field(meshPool, 'maxLevelCell', maxLevelCellField)
       call mpas_dmpar_exch_halo_field(maxLevelCellField)

       ! Set layerThickness based on refBottomDepth and bottomDepth
       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
          call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)

          call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

          call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
          call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

          call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

          call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

          do iCell = 1, nCellsSolve 
             if (maxLevelCell(iCell) > 0) then

                ! By going to maxLevelCell, this loop sets the layer Thickness as the full cell at the bottom.
                layerThickness(1, iCell) = refBottomDepth(1)
                do k = 2, maxLevelCell(iCell)
                   layerThickness(k, iCell) = refBottomDepth(k) - refBottomDepth(k-1)
                end do

                ! The following lines could be used for partial bottom cells, but only if the temperature is interpolated in the vertical as well.
                ! In version 3.0, one may alter the IC for partial bottom cells on start-up in MPAS.
                !k = maxLevelCell(iCell)
                !layerThickness(k, iCell) =  bottomDepth(iCell) - refBottomDepth(k-1)

                restingThickness(:, iCell) = layerThickness(:, iCell)
             end if
          end do

          block_ptr => block_ptr % next
       end do

    end subroutine ocn_init_setup_global_ocean_interpolate_topo!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_cull_inland_seas
!
!> \brief   Read the topography IC file
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details 
!>  This routine removes all inland seas. These are defined as isolated ocean cells.
!>   It uses a parallel version of an advancing front algorithm which might not be
!>   optimal for this purpose.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_cull_inland_seas(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (mpas_pool_type), pointer :: scratchPool, meshPool

       type (field1DInteger), pointer :: cullStackField, touchedCellField, oceanCellField

       real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, bottomDepth
       integer, dimension(:), pointer :: stack, oceanMask, touchMask
       integer, pointer :: stackSize

       real (kind=RKIND) :: currentLat, currentLon
       real (kind=RKIND) :: dist, minDist

       integer :: iCell
       integer :: localStackSize, globalStackSize
       integer :: j, coc
       integer :: touched

       integer, pointer :: nCells, nCellsSolve, nVertLevels
       integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
       integer, dimension(:, :), pointer :: cellsOnCell

       iErr = 0

       call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)

       call mpas_pool_get_field(scratchPool, 'cullStack', cullStackField)
       call mpas_pool_get_field(scratchPool, 'touchedCell', touchedCellField)
       call mpas_pool_get_field(scratchPool, 'oceanCell', oceanCellField)

       call mpas_allocate_scratch_field(cullStackField, .false.)
       call mpas_allocate_scratch_field(touchedCellField, .false.)
       call mpas_allocate_scratch_field(oceanCellField, .false.)

       ! Seed all deepest points for advancing front algorithm
       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

          call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
          call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

          call mpas_pool_get_array(meshPool, 'latCell', latCell)
          call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

          call mpas_pool_get_array(scratchPool, 'cullStack', stack)
          call mpas_pool_get_array(scratchPool, 'oceanCell', oceanMask)
          call mpas_pool_get_array(scratchPool, 'touchedCell', touchMask)
          call mpas_pool_get_array(scratchPool, 'cullStackSize', stackSize)

          stack(:) = 0
          oceanMask(:) = 0
          touchMask(:) = 0
          stackSize = 0

          ! Add all cells that have maxLevelCell == nVertLevels to stack
          do iCell = 1, nCellsSolve
             if (maxLevelCell(iCell) == nVertLevels) then
                stackSize = stackSize + 1
                stack(stackSize) = iCell
                touchMask(iCell) = 1
                oceanMask(iCell) = 1
             end if
          end do

          block_ptr => block_ptr % next
       end do

       ! Advancing front algorithm continues until all stacks on all processes are empty.
       globalStackSize = 1
       do while(globalStackSize /= 0)
          ! Advance front on each block with a non-zero stack until stack is empty.
          block_ptr => domain % blocklist
          do while(associated(block_ptr))
             call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
             call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

             call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
             call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
             call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

             call mpas_pool_get_array(scratchPool, 'cullStack', stack)
             call mpas_pool_get_array(scratchPool, 'oceanCell', oceanMask)
             call mpas_pool_get_array(scratchPool, 'touchedCell', touchMask)
             call mpas_pool_get_array(scratchPool, 'cullStackSize', stackSize)

             touched = 0
             do while(stackSize > 0)
                iCell = stack(stackSize)
                stackSize = stackSize - 1
                do j = 1, nEdgesOnCell(iCell)
                   coc = cellsOnCell(j, iCell)
                   if (touchMask(coc) == 0 .and. bottomDepth(coc) > 0.0_RKIND) then
                      oceanMask(coc) = 1
                      stackSize = stackSize + 1
                      stack(stackSize) = coc
                   end if
                   touchMask(coc) = 1
                   touched = touched + 1
                end do
             end do

             block_ptr => block_ptr % next
          end do

          ! Perform a halo exchange on oceanMask
          call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
          call mpas_pool_get_field(scratchPool, 'oceanCell', oceanCellField)
          call mpas_dmpar_exch_halo_field(oceanCellField)

          ! Check to see if any cells have been masked as ocean in the halo that have not been touched.
          ! If there are any, add them to the stack. Also, compute globalStackSize
          localStackSize = 0
          block_ptr => domain % blocklist
          do while(associated(block_ptr))
             call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
             call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

             call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
             call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

             call mpas_pool_get_array(scratchPool, 'cullStack', stack)
             call mpas_pool_get_array(scratchPool, 'oceanCell', oceanMask)
             call mpas_pool_get_array(scratchPool, 'touchedCell', touchMask)
             call mpas_pool_get_array(scratchPool, 'cullStackSize', stackSize)

             do iCell = nCellsSolve, nCells
                if (oceanMask(iCell) == 1 .and. touchMask(iCell) == 0) then
                   stackSize = stackSize + 1
                   stack(stackSize) = iCell
                   touchMask(iCell) = 1
                end if
             end do

             localStackSize = localStackSize + stackSize
             block_ptr => block_ptr % next
          end do

          call mpas_dmpar_sum_int(domain % dminfo, localStackSize, globalStackSize)
       end do

       ! Mark all cells that aren't ocean cells for removal
       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

          call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

          call mpas_pool_get_array(scratchPool, 'oceanCell', oceanMask)

          do iCell = 1, nCellsSolve
             if (oceanMask(iCell) == 0) then
                maxLevelCell(iCell) = -1
             end if
          end do
          block_ptr => block_ptr % next
       end do

       call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)

       call mpas_pool_get_field(scratchPool, 'cullStack', cullStackField)
       call mpas_pool_get_field(scratchPool, 'touchedCell', touchedCellField)
       call mpas_pool_get_field(scratchPool, 'oceanCell', oceanCellField)

       call mpas_deallocate_scratch_field(cullStackField, .false.)
       call mpas_deallocate_scratch_field(touchedCellField, .false.)
       call mpas_deallocate_scratch_field(oceanCellField, .false.)

       block_ptr => domain % blocklist
       do while (associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

          call ocn_mark_maxlevelcell(meshPool, iErr)
          block_ptr => block_ptr % next
       end do

    end subroutine ocn_init_setup_global_ocean_cull_inland_seas!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_depth_levels
!
!> \brief   Read depth levels for global ocean test case
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details 
!>  This routine reads the depth levels from the temperature IC file and sets
!>  refBottomDepth accordingly
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_depth_levels(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (MPAS_Stream_type) :: depthStream

       type (mpas_pool_type), pointer :: meshPool

       character (len=StrKIND), pointer :: config_global_ocean_depth_file, config_global_ocean_depth_varname, &
                                           config_global_ocean_depth_dimname

       real (kind=RKIND), pointer :: config_global_ocean_depth_conversion_factor

       integer :: k, iCell

       real (kind=RKIND), dimension(:), pointer :: refBottomDepth

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_file', config_global_ocean_depth_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_varname', config_global_ocean_depth_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_dimname', config_global_ocean_depth_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_conversion_factor', config_global_ocean_depth_conversion_factor)

       ! Define stream for depth levels
       call MPAS_createStream(depthStream, domain % iocontext, config_global_ocean_depth_file, MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)

       ! Setup depth field for stream to be read in
       depthIC % fieldName = trim(config_global_ocean_depth_varname)
       depthIC % dimSizes(1) = nDepth
       depthIC % dimNames(1) = trim(config_global_ocean_depth_dimname)
       depthIC % isVarArray = .false.
       depthIC % isPersistent = .true.
       depthIC % isActive = .true.
       depthIC % hasTimeDimension = .false.
       depthIC % block => domain % blocklist
       allocate(depthIC % array(nDepth))

       ! Add depth field to stream
       call MPAS_streamAddField(depthStream, depthIC, iErr)

       ! Read stream
       call MPAS_readStream(depthStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(depthStream)
       depthIC % array(:) = depthIC % array(:) * config_global_ocean_depth_conversion_factor

       ! Set refBottomDepth depending on depth levels. And convert appropriately
       block_ptr => domain % blocklist
       do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
         ! depthIC is the mid-depth of each layer.  Convert to bottom depth.
         refBottomDepth(1) = 2.0 * depthIC % array(1)
         do k=2,nDepth
            refBottomDepth(k) = refBottomDepth(k-1) + 2*(depthIC % array(k) - refBottomDepth(k-1))
         enddo

         block_ptr => block_ptr % next
       end do

    end subroutine ocn_init_setup_global_ocean_read_depth_levels!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_tracer_lat_lon
!
!> \brief   Read Lat/Lon for tracers in global ocean test case
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details 
!>  This routine reads the latitude and longitude coordinats for tracers from the temperature IC file.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_tracer_lat_lon(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (MPAS_Stream_type) :: tracerStream

       character (len=StrKIND), pointer :: config_global_ocean_temperature_file, config_global_ocean_tracer_lat_varname, &
                                           config_global_ocean_tracer_nlat_dimname, config_global_ocean_tracer_lon_varname, &
                                           config_global_ocean_tracer_nlon_dimname

       logical, pointer :: config_global_ocean_tracer_latlon_degrees

       integer :: iLat, iLon

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_temperature_file', config_global_ocean_temperature_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_lat_varname', config_global_ocean_tracer_lat_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_nlat_dimname', config_global_ocean_tracer_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_lon_varname', config_global_ocean_tracer_lon_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_nlon_dimname', config_global_ocean_tracer_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_latlon_degrees', config_global_ocean_tracer_latlon_degrees)

       ! Define stream for depth levels
       call MPAS_createStream(tracerStream, domain % iocontext, config_global_ocean_temperature_file, MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)

       ! Setup tracerLat and tracerLon fields for stream to be read in
       tracerLat % fieldName = trim(config_global_ocean_tracer_lat_varname)
       tracerLat % dimSizes(1) = nLatTracer
       tracerLat % dimNames(1) = trim(config_global_ocean_tracer_nlat_dimname)
       tracerLat % isVarArray = .false.
       tracerLat % isPersistent = .true.
       tracerLat % isActive = .true.
       tracerLat % hasTimeDimension = .false.
       tracerLat % block => domain % blocklist
       allocate(tracerLat % array(nLatTracer))

       tracerLon % fieldName = trim(config_global_ocean_tracer_lon_varname)
       tracerLon % dimSizes(1) = nLonTracer
       tracerLon % dimNames(1) = trim(config_global_ocean_tracer_nlon_dimname)
       tracerLon % isVarArray = .false.
       tracerLon % isPersistent = .true.
       tracerLon % isActive = .true.
       tracerLon % hasTimeDimension = .false.
       tracerLon % block => domain % blocklist
       allocate(tracerLon % array(nLonTracer))

       ! Add tracerLat and tracerLon fields to stream
       call MPAS_streamAddField(tracerStream, tracerLat, iErr)
       call MPAS_streamAddField(tracerStream, tracerLon, iErr)

       ! Read stream
       call MPAS_readStream(tracerStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(tracerStream)

       if (config_global_ocean_tracer_latlon_degrees) then
          do iLat = 1, nLatTracer
             tracerLat % array(iLat) = tracerLat % array(iLat) * pii / 180.0_RKIND
          end do

          do iLon = 1, nLonTracer
             tracerLon % array(iLon) = tracerLon % array(iLon) * pii / 180.0_RKIND
          end do
       end if

       do iLon = 1, nLonTracer
          if (tracerLon % array(iLon) < 0.0_RKIND) then
             tracerLon % array(iLon) = 2.0_RKIND * pii + tracerLon % array(iLon)
          end if
       end do

    end subroutine ocn_init_setup_global_ocean_read_tracer_lat_lon!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_temperature
!
!> \brief   Read temperature ICs for global ocean test case
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details 
!>  This routine reads the temperature field from the temperature IC file.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_temperature(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (MPAS_Stream_type) :: temperatureStream

       character (len=StrKIND), pointer :: config_global_ocean_temperature_file, config_global_ocean_temperature_varname, &
                                           config_global_ocean_tracer_nlon_dimname, config_global_ocean_tracer_nlat_dimname, &
                                           config_global_ocean_depth_dimname

       integer :: k

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_temperature_file', config_global_ocean_temperature_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_temperature_varname', config_global_ocean_temperature_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_nlon_dimname', config_global_ocean_tracer_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_nlat_dimname', config_global_ocean_tracer_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_dimname', config_global_ocean_depth_dimname)

       ! Define stream for temperature IC
       call MPAS_createStream(temperatureStream, domain % iocontext, config_global_ocean_temperature_file, MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)

       ! Setup temperature field for stream to be read in
       temperatureIC % fieldName = trim(config_global_ocean_temperature_varname)
       temperatureIC % dimSizes(1) = nLonTracer
       temperatureIC % dimSizes(2) = nLatTracer
       temperatureIC % dimSizes(3) = nDepth
       temperatureIC % dimNames(1) = trim(config_global_ocean_tracer_nlon_dimname)
       temperatureIC % dimNames(2) = trim(config_global_ocean_tracer_nlat_dimname)
       temperatureIC % dimNames(3) = trim(config_global_ocean_depth_dimname)
       temperatureIC % isVarArray = .false.
       temperatureIC % isPersistent = .true.
       temperatureIC % isActive = .true.
       temperatureIC % hasTimeDimension = .false.
       temperatureIC % block => domain % blocklist
       allocate(temperatureIC % array(nLonTracer, nLatTracer, nDepth))

       ! Add temperature field to stream
       call MPAS_streamAddField(temperatureStream, temperatureIC, iErr)

       ! Read stream
       call MPAS_readStream(temperatureStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(temperatureStream)

    end subroutine ocn_init_setup_global_ocean_read_temperature!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_salinity
!
!> \brief   Read salinity ICs for global ocean test case
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details 
!>  This routine reads the salinity field from the salinity IC file.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_salinity(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (MPAS_Stream_type) :: salinityStream

       character (len=StrKIND), pointer :: config_global_ocean_salinity_file, config_global_ocean_salinity_varname, &
                                           config_global_ocean_tracer_nlon_dimname, config_global_ocean_tracer_nlat_dimname, &
                                           config_global_ocean_depth_dimname

       integer :: k

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_salinity_file', config_global_ocean_salinity_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_salinity_varname', config_global_ocean_salinity_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_nlon_dimname', config_global_ocean_tracer_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_nlat_dimname', config_global_ocean_tracer_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_dimname', config_global_ocean_depth_dimname)

       ! Define stream for salinity IC
       call MPAS_createStream(salinityStream, domain % iocontext, config_global_ocean_salinity_file, MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)

       ! Setup salinity field for stream to be read in
       salinityIC % fieldName = trim(config_global_ocean_salinity_varname)
       salinityIC % dimSizes(1) = nLonTracer
       salinityIC % dimSizes(2) = nLatTracer
       salinityIC % dimSizes(3) = nDepth
       salinityIC % dimNames(1) = trim(config_global_ocean_tracer_nlon_dimname)
       salinityIC % dimNames(2) = trim(config_global_ocean_tracer_nlat_dimname)
       salinityIC % dimNames(3) = trim(config_global_ocean_depth_dimname)
       salinityIC % isVarArray = .false.
       salinityIC % isPersistent = .true.
       salinityIC % isActive = .true.
       salinityIC % hasTimeDimension = .false.
       salinityIC % block => domain % blocklist
       allocate(salinityIC % array(nLonTracer, nLatTracer, nDepth))

       ! Add salinity field to stream
       call MPAS_streamAddField(salinityStream, salinityIC, iErr)

       ! Read stream
       call MPAS_readStream(salinityStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(salinityStream)

    end subroutine ocn_init_setup_global_ocean_read_salinity!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_interoplate_tracers
!
!> \brief   Interpolate tracer quantities to MPAS grid
!> \author  Doug Jacobsen
!> \date    03/05/2014
!> \details 
!>  This routine interpolates the temperature/salinity data read in from the
!>  initial condition file to the MPAS grid. Currently it uses a nearest neighbor interpolation.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_interpolate_tracers(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr
       type (mpas_pool_type), pointer :: meshPool, statePool, scratchPool, tracersPool, forcingPool
       type (mpas_pool_type), pointer :: tracersSurfaceRestoringFieldsPool, tracersInteriorRestoringFieldsPool

       real (kind=RKIND) :: currentLat, currentLon, counter
       real (kind=RKIND) :: minDist, dist
       real (kind=RKIND) :: x, x1, x2, y, y1, y2, coef, coef11, coef12, coef21, coef22
       integer :: iLat, iLon, iSmooth, j, coc
       integer :: latSearch, lonSearch
       integer :: iCell, k
       integer :: xInd1, xInd2, yInd1, yInd2
       integer, pointer :: idxSalinity, idxTemperature, nCells, nCellsSolve, idxTracer1

       type (field2DReal), pointer :: smoothedTemperatureField, smoothedSalinityField
       type (field3DReal), pointer :: activeTracersField

       integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
       integer, dimension(:, :), pointer :: cellsOnCell

       real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
       real (kind=RKIND), dimension(:, :), pointer :: smoothedTemperature, smoothedSalinity
       real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, debugTracers
       real (kind=RKIND), dimension(:, :), pointer ::    activeTracersPistonVelocity, activeTracersSurfaceRestoringValue
       real (kind=RKIND), dimension(:, :, :), pointer :: activeTracersInteriorRestoringValue, activeTracersInteriorRestoringRate

       character (len=StrKIND), pointer :: config_global_ocean_tracer_method
       integer, pointer ::  config_global_ocean_smooth_TS_iterations
       real (kind=RKIND), pointer :: config_global_ocean_piston_velocity
       real (kind=RKIND), pointer :: config_global_ocean_interior_restore_rate

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_method', config_global_ocean_tracer_method)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_smooth_TS_iterations', config_global_ocean_smooth_TS_iterations)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_piston_velocity', config_global_ocean_piston_velocity)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_interior_restore_rate', config_global_ocean_interior_restore_rate)

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
          call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
          call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields', tracersSurfaceRestoringFieldsPool)
          call mpas_pool_get_subpool(forcingPool, 'tracersInteriorRestoringFields', tracersInteriorRestoringFieldsPool)

          call mpas_pool_get_dimension(tracersPool, 'index_temperature', idxTemperature)
          call mpas_pool_get_dimension(tracersPool, 'index_salinity', idxSalinity)
          call mpas_pool_get_dimension(tracersPool, 'index_tracer1', idxTracer1)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

          call mpas_pool_get_array(meshPool, 'latCell', latCell)
          call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

          call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
          call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)

          if (config_global_ocean_tracer_method .eq. "nearest_neighbor") then
             do iCell = 1, nCells
                currentLat = latCell(iCell)
                currentLon = lonCell(iCell)

                lonSearch = 1
                minDist = 2.0_RKIND * pii
                do iLon = 1, nLonTracer
                   dist = abs(currentLon - tracerLon % array(iLon))
                   if (dist < minDist) then
                      minDist = dist
                      lonSearch = iLon
                   end if
                end do

                latSearch = 1
                minDist = 2.0_RKIND * pii
                do iLat = 1, nLatTracer
                   dist = abs(currentLat - tracerLat % array(iLat))
                   if (dist < minDist) then
                      minDist = dist
                      latSearch = iLat
                   end if
                end do

                do k = 1, maxLevelCell(iCell)
                   if ( associated(activeTracers) ) then
                      activeTracers(idxTemperature, k, iCell) = temperatureIC % array(lonSearch, latSearch, k)
                      activeTracers(idxSalinity, k, iCell) = salinityIC % array(lonSearch, latSearch, k)
                   end if
                end do
             end do

          elseif (config_global_ocean_tracer_method .eq. "bilinear_interpolation") then

             do iCell = 1, nCells
                x = lonCell(iCell)
                y = latCell(iCell)

                ! Set up bilinear interpolation indices in longitude, watching for periodic boundary at 0 and 2 pi
                xInd1 = 0
                if (x .le. tracerLon % array(1)) then
                   xInd1 = nLonTracer
                   xInd2 = 1
                   x1 = tracerLon % array(xInd1) - 2.0*pii
                   x2 = tracerLon % array(xInd2)
                elseif (x .ge. tracerLon % array(nLonTracer)) then
                   xInd1 = nLonTracer
                   xInd2 = 1
                   x1 = tracerLon % array(xInd1)
                   x2 = tracerLon % array(xInd2) + 2.0*pii
                else
                   do iLon = 1, nLonTracer-1
                      if (x .le. tracerLon % array(iLon+1)) then
                         xInd1 = iLon
                         xInd2 = iLon+1
                         x1 = tracerLon % array(xInd1)
                         x2 = tracerLon % array(xInd2)
                         exit
                      end if
                   end do
                endif

                yInd1 = 0
                if (y .le. tracerLat % array(1)) then
                   ! if south of the southernmost data point, extrapolate as a constant in latitude
                   yInd1 = 1
                   yInd2 = 1
                   coef   = 1.0_RKIND/(x2-x1)
                   coef11 = 1.0_RKIND*(x2-x )
                   coef21 = 1.0_RKIND*(x -x1)
                   coef12 = 0.0_RKIND
                   coef22 = 0.0_RKIND
                elseif (y .ge. tracerLat % array(nLatTracer)) then
                   ! if north of the northernmost data point, extrapolate as a constant in latitude
                   yInd1 = nLatTracer
                   yInd2 = nLatTracer
                   coef   = 1.0_RKIND/(x2-x1)
                   coef11 = 1.0_RKIND*(x2-x )
                   coef21 = 1.0_RKIND*(x -x1)
                   coef12 = 0.0_RKIND
                   coef22 = 0.0_RKIND
                else
                   ! Set up bilinear interpolation coefficients in latitude
                   do iLat = 1, nLatTracer-1
                      if (y .le. tracerLat % array(iLat+1)) then
                         yInd1 = iLat
                         yInd2 = iLat+1
                         exit
                      end if
                   end do
                   y1 = tracerLat % array(yInd1)
                   y2 = tracerLat % array(yInd2)
                   coef   = 1.0_RKIND/(x2-x1)/(y2-y1)
                   coef11 = 1.0_RKIND*(x2-x )*(y2-y )
                   coef21 = 1.0_RKIND*(x -x1)*(y2-y )
                   coef12 = 1.0_RKIND*(x2-x )*(y -y1)
                   coef22 = 1.0_RKIND*(x -x1)*(y -y1)
                endif

                ! Assign T&S using bilinear interpolation
                ! formulas from http://en.wikipedia.org/wiki/Bilinear_interpolation
                do k = 1, maxLevelCell(iCell)

                   if ( associated(activeTracers) ) then
                      activeTracers(idxTemperature, k, iCell) = coef*( &
                             coef11* temperatureIC % array(xInd1,yInd1, k) &
                           + coef21* temperatureIC % array(xInd2,yInd1, k) &
                           + coef12* temperatureIC % array(xInd1,yInd2, k) &
                           + coef22* temperatureIC % array(xInd2,yInd2, k) )

                      activeTracers(idxSalinity, k, iCell) = coef*( &
                             coef11* salinityIC % array(xInd1,yInd1, k) &
                           + coef21* salinityIC % array(xInd2,yInd1, k) &
                           + coef12* salinityIC % array(xInd1,yInd2, k) &
                           + coef22* salinityIC % array(xInd2,yInd2, k) )
                   end if

                end do

             end do

          else
             write(stderrUnit,*) 'ERROR: Invalid choice of config_global_ocean_tracer_method.'
             iErr = 1
             call mpas_dmpar_finalize(domain % dminfo)
          endif


          block_ptr => block_ptr % next
       end do

       ! Smooth temperature and salinity.
       if (config_global_ocean_smooth_TS_iterations .gt. 0) then
          call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)

          call mpas_pool_get_field(scratchPool, 'smoothedTemperature', smoothedTemperatureField)
          call mpas_pool_get_field(scratchPool, 'smoothedSalinity', smoothedSalinityField)

          call mpas_allocate_scratch_field(smoothedTemperatureField, .false.)
          call mpas_allocate_scratch_field(smoothedSalinityField, .false.)

          do iSmooth = 1,config_global_ocean_smooth_TS_iterations

             block_ptr => domain % blocklist
             do while(associated(block_ptr))
                call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
                call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
                call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)
                call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
                call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields', tracersSurfaceRestoringFieldsPool)
                call mpas_pool_get_subpool(forcingPool, 'tracersInteriorRestoringFields', tracersInteriorRestoringFieldsPool)

                call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
                call mpas_pool_get_dimension(tracersPool, 'index_temperature', idxTemperature)
                call mpas_pool_get_dimension(tracersPool, 'index_salinity', idxSalinity)

                call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
                call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
                call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)

                call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)

                call mpas_pool_get_array(scratchPool, 'smoothedTemperature', smoothedTemperature)
                call mpas_pool_get_array(scratchPool, 'smoothedSalinity', smoothedSalinity)

                maxLevelCell(nCells+1) = -1

                do iCell = 1, nCells
                   do k = 1, maxLevelCell(iCell)
                      if ( associated(activeTracers) ) then
                         smoothedtemperature(k, iCell) = activeTracers(idxTemperature, k, iCell)
                         smoothedsalinity(k, iCell) = activeTracers(idxSalinity, k, iCell)
                      end if
                      counter = 1

                      do j = 1, nEdgesOnCell(iCell)
                         coc = cellsOnCell(j, iCell)
                         ! check if coc not 0 (or nCells+1)?
                         if (k .le. maxLevelCell(coc)) then

                            if ( associated(activeTracers) ) then
                               smoothedtemperature(k, iCell) = smoothedtemperature(k, iCell) + activeTracers (idxTemperature, k, coc)
                               smoothedsalinity(k, iCell) = smoothedsalinity(k, iCell) + activeTracers(idxSalinity, k, coc)
                            end if
                            counter = counter + 1

                         end if
                      end do ! edgesOnCell

                      smoothedtemperature(k, iCell) = smoothedtemperature(k, iCell) / counter
                      smoothedsalinity(k, iCell) = smoothedsalinity(k, iCell) / counter

                   end do ! k level

                end do ! iCell

                if ( associated(activeTracers) ) then
                   activeTracers(idxTemperature, :, :) = smoothedtemperature(:,:) 
                   activeTracers(idxSalinity, :, :) = smoothedsalinity(:,:) 
                end if

                block_ptr => block_ptr % next
             end do

             call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
             call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
             call mpas_pool_get_field(tracersPool, 'activeTracers', activeTracersField,1)

             call mpas_dmpar_exch_halo_field(activeTracersField)

          end do ! iSmooth

          call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
          call mpas_pool_get_field(scratchPool, 'smoothedTemperature', smoothedTemperatureField)
          call mpas_pool_get_field(scratchPool, 'smoothedSalinity', smoothedSalinityField)
          call mpas_deallocate_scratch_field(smoothedTemperatureField, .false.)
          call mpas_deallocate_scratch_field(smoothedSalinityField, .false.)
       endif

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
          call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

          call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
          call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)

          call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields', tracersSurfaceRestoringFieldsPool)
          call mpas_pool_get_subpool(forcingPool, 'tracersInteriorRestoringFields', tracersInteriorRestoringFieldsPool)
          call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersSurfaceRestoringValue', activeTracersSurfaceRestoringValue, 1)
          call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersPistonVelocity', activeTracersPistonVelocity, 1)
          call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, 'activeTracersInteriorRestoringValue', activeTracersInteriorRestoringValue, 1)
          call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, 'activeTracersInteriorRestoringRate', activeTracersInteriorRestoringRate, 1)

          ! set interior restoring values and rate
          if ( associated(activeTracersInteriorRestoringValue) .and. associated(activeTracers) ) then
             activeTracersInteriorRestoringValue(:, :, :) = activeTracers(:, :, :)
          end if

          if ( associated(activeTracersInteriorRestoringRate) ) then
             activeTracersInteriorRestoringRate(:, :, :) = config_global_ocean_interior_restore_rate
          end if

           ! set surface restoring values and rate
          if ( associated(activeTracersSurfaceRestoringValue) .and. associated(activeTracers) ) then
             activeTracersSurfaceRestoringValue(:, :) = activeTracers(:, 1, :)
          end if

          if ( associated(activeTracersPistonVelocity) ) then
             activeTracersPistonVelocity(:, :) = config_global_ocean_piston_velocity
          end if

          block_ptr => block_ptr % next
       end do


    end subroutine ocn_init_setup_global_ocean_interpolate_tracers!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_windstress
!
!> \brief   Read the windstress IC file
!> \author  Doug Jacobsen
!> \date    03/07/2014
!> \details 
!>  This routine reads the windstress IC file, including latitude and longitude
!>   information for windstress data.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_windstress(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (MPAS_Stream_type) :: windstressStream

       integer :: iLat, iLon

       character (len=StrKIND), pointer :: config_global_ocean_windstress_file, config_global_ocean_windstress_lat_varname, &
                                           config_global_ocean_windstress_nlat_dimname, config_global_ocean_windstress_lon_varname, &
                                           config_global_ocean_windstress_nlon_dimname, config_global_ocean_windstress_zonal_varname, &
                                           config_global_ocean_windstress_meridional_varname

       logical, pointer :: config_global_ocean_windstress_latlon_degrees

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_file', config_global_ocean_windstress_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_lat_varname', config_global_ocean_windstress_lat_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_nlat_dimname', config_global_ocean_windstress_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_lon_varname', config_global_ocean_windstress_lon_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_nlon_dimname', config_global_ocean_windstress_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_zonal_varname', config_global_ocean_windstress_zonal_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_meridional_varname', config_global_ocean_windstress_meridional_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_latlon_degrees', config_global_ocean_windstress_latlon_degrees)

       ! Define stream for depth levels
       call MPAS_createStream(windstressStream, domain % iocontext, config_global_ocean_windstress_file, MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)

       ! Setup windLat, windLon, and windIC fields for stream to be read in
       windLat % fieldName = trim(config_global_ocean_windstress_lat_varname)
       windLat % dimSizes(1) = nLatWind
       windLat % dimNames(1) = trim(config_global_ocean_windstress_nlat_dimname)
       windLat % isVarArray = .false.
       windLat % isPersistent = .true.
       windLat % isActive = .true.
       windLat % hasTimeDimension = .false.
       windLat % block => domain % blocklist
       allocate(windLat % array(nLatWind))

       windLon % fieldName = trim(config_global_ocean_windstress_lon_varname)
       windLon % dimSizes(1) = nLonWind
       windLon % dimNames(1) = trim(config_global_ocean_windstress_nlon_dimname)
       windLon % isVarArray = .false.
       windLon % isPersistent = .true.
       windLon % isActive = .true.
       windLon % hasTimeDimension = .false.
       windLon % block => domain % blocklist
       allocate(windLon % array(nLonWind))

       zonalWindIC % fieldName = trim(config_global_ocean_windstress_zonal_varname)
       zonalWindIC % dimSizes(1) = nLonWind
       zonalWindIC % dimSizes(2) = nLatWind
       zonalWindIC % dimNames(1) = trim(config_global_ocean_windstress_nlon_dimname)
       zonalWindIC % dimNames(2) = trim(config_global_ocean_windstress_nlat_dimname)
       zonalWindIC % isVarArray = .false.
       zonalWindIC % isPersistent = .true.
       zonalWindIC % isActive = .true.
       zonalWindIC % hasTimeDimension = .false.
       zonalWindIC % block => domain % blocklist
       allocate(zonalWindIC % array(nLonWind, nLatWind))

       meridionalWindIC % fieldName = trim(config_global_ocean_windstress_meridional_varname)
       meridionalWindIC % dimSizes(1) = nLonWind
       meridionalWindIC % dimSizes(2) = nLatWind
       meridionalWindIC % dimNames(1) = trim(config_global_ocean_windstress_nlon_dimname)
       meridionalWindIC % dimNames(2) = trim(config_global_ocean_windstress_nlat_dimname)
       meridionalWindIC % isVarArray = .false.
       meridionalWindIC % isPersistent = .true.
       meridionalWindIC % isActive = .true.
       meridionalWindIC % hasTimeDimension = .false.
       meridionalWindIC % block => domain % blocklist
       allocate(meridionalWindIC % array(nLonWind, nLatWind))

       ! Add windLat, windLon, and windIC fields to stream
       call MPAS_streamAddField(windstressStream, windLat, iErr)
       call MPAS_streamAddField(windstressStream, windLon, iErr)
       call MPAS_streamAddField(windstressStream, zonalWindIC, iErr)
       call MPAS_streamAddField(windstressStream, meridionalWindIC, iErr)

       ! Read stream
       call MPAS_readStream(windstressStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(windstressStream)

       if (config_global_ocean_windstress_latlon_degrees) then
          windLat % array(:) = windLat % array(:) * pii / 180.0_RKIND 
          windLon % array(:) = windLon % array(:) * pii / 180.0_RKIND
       end if

       do iLon = 1, nLonWind
          if (windLon % array(iLon) < 0.0_RKIND) then
             windLon % array(iLon) = 2.0_RKIND * pii + windLon % array(iLon) 
          end if
       end do

    end subroutine ocn_init_setup_global_ocean_read_windstress!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_interpolate_windstress
!
!> \brief   Interpolate the windstress IC to MPAS mesh
!> \author  Doug Jacobsen
!> \date    03/07/2014
!> \details 
!>  This routine interpolates windstress data to the MPAS mesh. Currently it
!>   uses a bilinear interpolation
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_interpolate_windstress(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (mpas_pool_type), pointer :: meshPool, forcingPool

       real (kind=RKIND) :: currentLat, currentLon
       real (kind=RKIND) :: zonalWind, meridionalWind
       real (kind=RKIND) :: angle
       real (kind=RKIND) :: dist, minDist
       real (kind=RKIND) :: x, x1, x2, y, y1, y2, coef, coef11, coef12, coef21, coef22

       integer :: ilat, iLon
       integer :: latSearch, lonSearch
       integer :: iEdge
       integer :: xInd1, xInd2, yInd1, yInd2

       real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge, angleEdge, surfaceStress

       integer, pointer :: nEdgesSolve, nEdges

       character (len=StrKIND), pointer :: config_global_ocean_windstress_method
       real (kind=RKIND), pointer :: config_global_ocean_windstress_conversion_factor

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_method', config_global_ocean_windstress_method)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_conversion_factor', config_global_ocean_windstress_conversion_factor)

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

          call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
          call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

          call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
          call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge)
          call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)

          call mpas_pool_get_array(forcingPool, 'surfaceStress', surfaceStress)

          if (config_global_ocean_windstress_method .eq. "nearest_neighbor") then
             do iEdge = 1, nEdgesSolve
                currentLat = latEdge(iEdge)
                currentLon = lonEdge(iEdge)
                angle = angleEdge(iEdge)

                minDist = 2.0_RKIND * pii
                lonSearch = 1
                do iLon = 1, nLonWind
                   dist = abs(currentLon - windLon % array(iLon))
                   if (dist < minDist) then
                      minDist = dist
                      lonSearch = iLon
                   end if
                end do

                minDist = 2.0_RKIND * pii
                latSearch = 1
                do iLat = 1, nLatWind
                   dist = abs(currentLat - windLat % array(iLat))
                   if (dist < minDist) then
                      minDist = dist
                      latSearch = iLat
                   end if
                end do

                zonalWind = zonalWindIC % array(lonSearch, latSearch) * config_global_ocean_windstress_conversion_factor
                meridionalWind = meridionalWindIC % array(lonSearch, latSearch) * config_global_ocean_windstress_conversion_factor

                surfaceStress(iEdge) = zonalWind * cos(angle) + meridionalWind * sin(angle)
             end do

          elseif (config_global_ocean_windstress_method .eq. "bilinear_interpolation") then

             do iEdge = 1, nEdges
                x = lonEdge(iEdge)
                y = latEdge(iEdge)
                angle = angleEdge(iEdge)

                ! Set up bilinear interpolation indices in longitude, watching for periodic boundary at 0 and 2 pi
                xInd1 = 0
                if (x .le. windLon % array(1)) then
                   xInd1 = nLonWind
                   xInd2 = 1
                   x1 = windLon % array(xInd1) - 2.0_RKIND*pii
                   x2 = windLon % array(xInd2)
                elseif (x .ge. windLon % array(nLonWind)) then
                   xInd1 = nLonWind
                   xInd2 = 1
                   x1 = windLon % array(xInd1)
                   x2 = windLon % array(xInd2) + 2.0_RKIND*pii
                else
                   do iLon = 1, nLonWind-1
                      if (x .le. windLon % array(iLon+1)) then
                         xInd1 = iLon
                         xInd2 = iLon+1
                         x1 = windLon % array(xInd1)
                         x2 = windLon % array(xInd2)
                         exit
                      end if
                   end do
                endif

                yInd1 = 0
                if (y .le. windLat % array(1)) then
                   ! if south of the southernmost data point, extrapolate as a constant in latitude
                   yInd1 = 1
                   yInd2 = 1
                   coef   = 1.0_RKIND/(x2-x1)
                   coef11 = 1.0_RKIND*(x2-x )
                   coef21 = 1.0_RKIND*(x -x1)
                   coef12 = 0.0_RKIND
                   coef22 = 0.0_RKIND
                elseif (y .ge. windLat % array(nLatWind)) then
                   ! if north of the northernmost data point, extrapolate as a constant in latitude
                   yInd1 = nLatWind
                   yInd2 = nLatWind
                   coef   = 1.0_RKIND/(x2-x1)
                   coef11 = 1.0_RKIND*(x2-x )
                   coef21 = 1.0_RKIND*(x -x1)
                   coef12 = 0.0_RKIND
                   coef22 = 0.0_RKIND
                else
                   ! Set up bilinear interpolation coefficients in latitude
                   do iLat = 1, nLatWind-1
                      if (y .le. windLat % array(iLat+1)) then
                         yInd1 = iLat
                         yInd2 = iLat+1
                         exit
                      end if
                   end do
                   y1 = windLat % array(yInd1)
                   y2 = windLat % array(yInd2)
                   coef   = 1.0_RKIND/(x2-x1)/(y2-y1)
                   coef11 = 1.0_RKIND*(x2-x )*(y2-y )
                   coef21 = 1.0_RKIND*(x -x1)*(y2-y )
                   coef12 = 1.0_RKIND*(x2-x )*(y -y1)
                   coef22 = 1.0_RKIND*(x -x1)*(y -y1)
                endif

                zonalWind = coef*config_global_ocean_windstress_conversion_factor*( &
                     coef11* zonalWindIC % array(xInd1, yInd1) &
                     + coef21* zonalWindIC % array(xInd2, yInd1) &
                     + coef12* zonalWindIC % array(xInd1, yInd2) &
                     + coef22* zonalWindIC % array(xInd2, yInd2) )

                meridionalWind = coef*config_global_ocean_windstress_conversion_factor*( &
                     coef11* meridionalWindIC % array(xInd1, yInd1) &
                     + coef21* meridionalWindIC % array(xInd2, yInd1) &
                     + coef12* meridionalWindIC % array(xInd1, yInd2) &
                     + coef22* meridionalWindIC % array(xInd2, yInd2) )

                surfaceStress(iEdge) = zonalWind * cos(angle) + meridionalWind * sin(angle)

             end do

          else
             write(stderrUnit,*) 'ERROR: Invalid choice of config_global_ocean_windstress_method.'
             iErr = 1
             call mpas_dmpar_finalize(domain % dminfo)
          endif

          block_ptr => block_ptr % next
       end do

    end subroutine ocn_init_setup_global_ocean_interpolate_windstress!}}}

!***********************************************************************
!
!  routine ocn_init_global_ocean_destroy_tracer_fields
!
!> \brief   Tracer field cleanup routine
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details 
!>  This routine destroys the fields that were created to hold tracer
!>  initial condition information
!
!-----------------------------------------------------------------------

    subroutine ocn_init_global_ocean_destroy_tracer_fields()!{{{
        deallocate(temperatureIC % array)
        deallocate(salinityIC % array)
        deallocate(tracerLat % array)
        deallocate(tracerLon % array)
    end subroutine ocn_init_global_ocean_destroy_tracer_fields!}}}

!***********************************************************************
!
!  routine ocn_init_global_ocean_destroy_topo_fields
!
!> \brief   Topography field cleanup routine
!> \author  Doug Jacobsen
!> \date    03/07/2014
!> \details 
!>  This routine destroys the fields that were created to hold topography
!>  initial condition information
!
!-----------------------------------------------------------------------

    subroutine ocn_init_global_ocean_destroy_topo_fields()!{{{
        deallocate(topoIC % array)
        deallocate(topoLat % array)
        deallocate(topoLon % array)
    end subroutine ocn_init_global_ocean_destroy_topo_fields!}}}

!***********************************************************************
!
!  routine ocn_init_global_ocean_destroy_windstress_fields
!
!> \brief   Windstress field cleanup routine
!> \author  Doug Jacobsen
!> \date    03/07/2014
!> \details 
!>  This routine destroys the fields that were created to hold windstress
!>  initial condition information
!
!-----------------------------------------------------------------------

    subroutine ocn_init_global_ocean_destroy_windstress_fields()!{{{
        deallocate(zonalWindIC % array)
        deallocate(meridionalWindIC % array)
        deallocate(windLat % array)
        deallocate(windLon % array)
    end subroutine ocn_init_global_ocean_destroy_windstress_fields!}}}

!***********************************************************************
!
!  routine ocn_init_validate_global_ocean
!
!> \brief   Validation for global ocean test case
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details 
!>  This routine validates the configuration options for the global ocean test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_global_ocean(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout), target :: iocontext

      type (mpas_io_context_type), pointer :: iocontext_ptr

      integer, intent(out) :: iErr
      type (MPAS_IO_Handle_type) :: inputFile

      character (len=StrKIND), pointer :: config_init_configuration, config_global_ocean_depth_file, &
                                          config_global_ocean_depth_dimname, config_global_ocean_temperature_file, &
                                          config_global_ocean_salinity_file, config_global_ocean_tracer_nlat_dimname, &
                                          config_global_ocean_tracer_nlon_dimname, config_global_ocean_topography_file, &
                                          config_global_ocean_topography_nlat_dimname, config_global_ocean_topography_nlon_dimname, &
                                          config_global_ocean_windstress_file, config_global_ocean_windstress_nlat_dimname, &
                                          config_global_ocean_windstress_nlon_dimname

      integer, pointer :: config_vert_levels

      iocontext_ptr => iocontext

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('global_ocean')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_global_ocean_depth_file', config_global_ocean_depth_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_depth_dimname', config_global_ocean_depth_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_temperature_file', config_global_ocean_temperature_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_salinity_file', config_global_ocean_salinity_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_nlat_dimname', config_global_ocean_tracer_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_nlon_dimname', config_global_ocean_tracer_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_file', config_global_ocean_topography_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_nlat_dimname', config_global_ocean_topography_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_nlon_dimname', config_global_ocean_topography_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_file', config_global_ocean_windstress_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_nlat_dimname', config_global_ocean_windstress_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_nlon_dimname', config_global_ocean_windstress_nlon_dimname)

      inputFile = MPAS_io_open(config_global_ocean_depth_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)

      call MPAS_io_inq_dim(inputFile, config_global_ocean_depth_dimname, nDepth, iErr)

      call MPAS_io_close(inputFile, iErr)

      inputFile = MPAS_io_open(config_global_ocean_temperature_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)

      call MPAS_io_inq_dim(inputFile, config_global_ocean_tracer_nlat_dimname, nLatTracer, iErr)
      call MPAS_io_inq_dim(inputFile, config_global_ocean_tracer_nlon_dimname, nLonTracer, iErr)

      call MPAS_io_close(inputFile, iErr)

      inputFile = MPAS_io_open(config_global_ocean_topography_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)

      call MPAS_io_inq_dim(inputFile, config_global_ocean_topography_nlat_dimname, nLatTopo, iErr)
      call MPAS_io_inq_dim(inputFile, config_global_ocean_topography_nlon_dimname, nLonTopo, iErr)

      call MPAS_io_close(inputFile, iErr)

      inputFile = MPAS_io_open(config_global_ocean_windstress_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)

      call MPAS_io_inq_dim(inputFile, config_global_ocean_windstress_nlat_dimname, nLatWind, iErr)
      call MPAS_io_inq_dim(inputFile, config_global_ocean_windstress_nlon_dimname, nLonWind, iErr)

      call MPAS_io_close(inputFile, iErr)

      if (config_vert_levels <= 0 .and. nDepth > 0) then
         config_vert_levels = nDepth
      else if(config_vert_levels <= 0) then
         write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Not given a usable value for vertical levels.'
         iErr = 1
      end if

      if (trim(config_global_ocean_temperature_file) == 'none') then
         write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Invalid filename for config_global_ocean_temperature_file'
         iErr = 1
      end if

      if (trim(config_global_ocean_salinity_file) == 'none') then
         write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Invalid filename for config_global_ocean_salinity_file'
         iErr = 1
      end if

      if (trim(config_global_ocean_depth_file) == 'none') then
         write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Invalid filename for config_global_ocean_depth_file'
         iErr = 1
      end if

      if (trim(config_global_ocean_topography_file) == 'none') then
         write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Invalid filename for config_global_ocean_topography_file'
         iErr = 1
      end if

      if (trim(config_global_ocean_windstress_file) == 'none') then
         write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Invalid filename for config_global_ocean_windstress_file'
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_global_ocean!}}}

!***********************************************************************

end module ocn_init_global_ocean

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
