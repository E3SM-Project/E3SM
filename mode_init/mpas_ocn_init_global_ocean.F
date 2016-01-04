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
   use mpas_stream_manager
   use mpas_timekeeping
   use mpas_dmpar

   use ocn_constants
   use ocn_init_cell_markers
   use ocn_init_vertical_grids
   use ocn_init_interpolation

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
   integer :: nLatTracer, nLonTracer, nLonSW, nLatSW
   integer :: nLatWind, nLonWind
   integer :: nLatTopo, nLonTopo
   type (field1DReal) :: depthIC
   type (field1DReal) :: windLat, windLon
   type (field1DReal) :: topoLat, topoLon
   type (field1DReal) :: tracerLat, tracerLon
   type (field1DReal) :: swDataLat, swDataLon
   type (field2DReal) :: topoIC, zonalWindIC, meridionalWindIC, chlorophyllIC, zenithAngleIC, clearSkyIC
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

      character (len=StrKIND), pointer :: config_init_configuration, config_sw_absorption_type
      logical, pointer :: config_global_ocean_cull_inland_seas

      logical, pointer :: on_a_sphere

      iErr = 0

      call mpas_pool_get_config(domain % configs, 'config_init_configuration', config_init_configuration)
      call mpas_pool_get_config(domain % configs, 'config_sw_absorption_type', config_sw_absorption_type)

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
      call ocn_init_setup_global_ocean_create_model_topo(domain, iErr)
      write(stderrUnit,*) 'Cleaning up topography IC fields'
      call ocn_init_global_ocean_destroy_topo_fields()

      if(trim(config_sw_absorption_type) == 'ohlmann00') then

      write(stderrUnit,*) 'Reading penetrating shortwave lat/lon data'
      call ocn_init_setup_global_ocean_read_swData_lat_lon(domain,iErr)
      write(stderrUnit,*) 'Interpolating penetrating shortwave data'
      call ocn_init_setup_global_ocean_interpolate_swData(domain,iErr)
      write(stderrUnit,*) 'Cleaning penetrating shortwave data'
      call ocn_init_global_ocean_destroy_swData_fields()
      
      endif
      
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
!  routine ocn_init_setup_global_ocean_create_model_topo
!
!> \brief   Interpolate the topography IC to MPAS mesh
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details 
!>  This routine interpolates topography data to the MPAS mesh. Currently it
!>   uses a bilinear interpolation
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_create_model_topo(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (mpas_pool_type), pointer :: meshPool, scratchPool, statePool, verticalMeshPool, diagnosticsPool

       real (kind=RKIND) :: currentLat, currentLon
       real (kind=RKIND) :: dist, minDist
       real (kind=RKIND) :: alpha, beta, depthLat1, depthLat2, proposedDepth

       real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, bottomDepth, bottomDepthObserved, &
            refBottomDepth, refLayerThickness, refZMid
       real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness, zMid

       integer, pointer :: nCells, nCellsSolve, nVertLevels

       type (field1DInteger), pointer :: maxLevelCellField, smoothedLevelsField
       integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
       integer, dimension(:, :), pointer :: cellsOnCell

       integer :: latSearch, lonSearch, searchIdx
       integer :: iCell, coc, j, k, maxLevel
       integer :: minimum_levels

       character (len=StrKIND), pointer :: config_global_ocean_topography_method
       logical, pointer :: config_global_ocean_smooth_topography
       real (kind=RKIND), pointer :: config_global_ocean_minimum_depth

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_method', config_global_ocean_topography_method)
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

          ! Record depth of the bottom of the ocean, before any alterations for modeling purposes.
          if (config_global_ocean_topography_method .eq. "nearest_neighbor") then
 
             call ocn_init_interpolation_nearest_horiz(topoLon % array, topoLat % array, &
                                                       topoIC % array, nLonTopo, nLatTopo, &
                                                       lonCell, latCell, bottomDepthObserved, nCells, &
                                                       inXPeriod = 2.0_RKIND * pii)
          elseif (config_global_ocean_topography_method .eq. "bilinear_interpolation") then
             call ocn_init_interpolation_bilinear_horiz(topoLon % array, topoLat % array, &
                                                        topoIC % array, nLonTopo, nLatTopo, &
                                                        lonCell, latCell, bottomDepthObserved, nCells, &
                                                        inXPeriod = 2.0_RKIND * pii)
          else
             write(stderrUnit,*) 'ERROR: Invalid choice of config_global_ocean_topography_method.'
             iErr = 1
             call mpas_dmpar_finalize(domain % dminfo)
          endif

          do iCell = 1, nCells
             ! Record depth of the bottom of the ocean, before any alterations for modeling purposes.
             ! Flip the sign to positive down.
             bottomDepthObserved(iCell) = -bottomDepthObserved(iCell)
             if (bottomDepthObserved(iCell) > 0.0_RKIND) then

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

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
          call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

          call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
          call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

          call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
          call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)
          call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)

          call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
          call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

          ! Compute refLayerThickness and refZMid
          call ocn_compute_layerThickness_zMid_from_bottomDepth(refLayerThickness,refZMid, &
               refBottomDepth,refBottomDepth(nVertLevels), &
               nVertLevels,nVertLevels,iErr)

          do iCell = 1, nCells 
             call ocn_alter_bottomDepth_for_pbcs(bottomDepth(iCell), refBottomDepth, maxLevelCell(iCell), iErr)

             ! Compute LayerThickness and zMid
             call ocn_compute_layerThickness_zMid_from_bottomDepth(layerThickness(:,iCell),zMid(:,iCell), &
                  refBottomDepth,bottomDepth(iCell), &
                  maxLevelCell(iCell),nVertLevels,iErr)

             ! Compute restingThickness
             restingThickness(:, iCell) = layerThickness(:, iCell)
          end do

          call ocn_compute_Haney_number(domain % dminfo, meshPool, diagnosticsPool, iErr)

          block_ptr => block_ptr % next
       end do

     end subroutine ocn_init_setup_global_ocean_create_model_topo!}}}

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
!  routine ocn_init_setup_global_ocean_read_swData_lat_lon
!
!> \brief   Read Lat/Lon for swData in global ocean test case
!> \author  Luke Van Roekel
!> \date    11/16/2015
!> \details 
!>  This routine reads the latitude and longitude coordinats for swData (chlorophyll, 
!           clearSkyRadiation, zenithangle) from the swData IC file.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_swData_lat_lon(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (MPAS_Stream_type) :: SWStream

       character (len=StrKIND), pointer :: config_global_ocean_swData_file, config_global_ocean_swData_lat_varname, &
                                           config_global_ocean_swData_nlat_dimname, config_global_ocean_swData_lon_varname, &
                                           config_global_ocean_swData_nlon_dimname

       logical, pointer :: config_global_ocean_swData_latlon_degrees

       integer :: iLat, iLon

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_file', config_global_ocean_swData_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_lat_varname', config_global_ocean_swData_lat_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_nlat_dimname', config_global_ocean_swData_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_lon_varname', config_global_ocean_swData_lon_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_nlon_dimname', config_global_ocean_swData_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_latlon_degrees', config_global_ocean_swData_latlon_degrees)

       ! Define stream for depth levels
       call MPAS_createStream(SWStream, domain % iocontext, config_global_ocean_swData_file, MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)

       ! Setup tracerLat and tracerLon fields for stream to be read in
       swDataLat % fieldName = trim(config_global_ocean_swData_lat_varname)
       swDataLat % dimSizes(1) = nLatSW
       swDataLat % dimNames(1) = trim(config_global_ocean_swData_nlat_dimname)
       swDataLat % isVarArray = .false.
       swDataLat % isPersistent = .true.
       swDataLat % isActive = .true.
       swDataLat % hasTimeDimension = .false.
       swDataLat % block => domain % blocklist
       allocate(swDataLat % array(nLatSW))

       swDataLon % fieldName = trim(config_global_ocean_swData_lon_varname)
       swDataLon % dimSizes(1) = nLonSW
       swDataLon % dimNames(1) = trim(config_global_ocean_swData_nlon_dimname)
       swDataLon % isVarArray = .false.
       swDataLon % isPersistent = .true.
       swDataLon % isActive = .true.
       swDataLon % hasTimeDimension = .false.
       swDataLon % block => domain % blocklist
       allocate(swDataLon % array(nLonSW))

       ! Add tracerLat and tracerLon fields to stream
       call MPAS_streamAddField(SWStream, swDataLat, iErr)
       call MPAS_streamAddField(SWStream, swDataLon, iErr)

       ! Read stream
       call MPAS_readStream(SWStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(SWStream)

       if (config_global_ocean_swData_latlon_degrees) then
          do iLat = 1, nLatSW
             swDataLat % array(iLat) = swDataLat % array(iLat) * pii / 180.0_RKIND
          end do

          do iLon = 1, nLonSW
             swDataLon % array(iLon) = swDataLon % array(iLon) * pii / 180.0_RKIND
          end do
       end if

       do iLon = 1, nLonSW
          if (swDataLon % array(iLon) < 0.0_RKIND) then
             swDataLon % array(iLon) = 2.0_RKIND * pii + swDataLon % array(iLon)
          end if
       end do

    end subroutine ocn_init_setup_global_ocean_read_swData_lat_lon!}}}
    
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
       real (kind=RKIND) :: zMidPBC
       integer :: iLat, iLon, iSmooth, j, coc
       integer :: latSearch, lonSearch
       integer :: iCell, k, km1
       integer :: xInd1, xInd2, yInd1, yInd2
       integer, pointer :: idxSalinity, idxTemperature, nCells, nVertLevels, nCellsSolve, idxTracer1

       type (field2DReal), pointer :: smoothedTemperatureField, smoothedSalinityField, interpTracerField
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

       ! These variables needed to interpolate tracers for partial bottom cells.
       ! Might remove once we interpolate to an arbitrary vertical grid.
       logical, pointer :: config_alter_ICs_for_pbcs
       character (len=StrKIND), pointer :: config_pbc_alteration_type
       type (mpas_pool_iterator_type) :: groupItr
       real (kind=RKIND), dimension(:), pointer :: bottomDepth, refBottomDepth
       real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroup
       real (kind=RKIND), dimension(:), allocatable :: zMidZLevel

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

          if(associated(debugTracers)) then
             debugTracers(idxTracer1,:,:) = 1.0_RKIND
          end if

          if ( associated(activeTracers) ) then
             call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)
             call mpas_pool_get_field(scratchPool, 'interpTracer', interpTracerField)
             call mpas_allocate_scratch_field(interpTracerField, .true.)
             if (config_global_ocean_tracer_method .eq. "nearest_neighbor") then
                call ocn_init_interpolation_nearest_horiz(tracerLon % array, tracerLat % array, &
                                                          temperatureIC % array, nLonTracer, nLatTracer, &
                                                          lonCell, latCell, interpTracerField % array, nCells, &
                                                          inXPeriod = 2.0_RKIND * pii)
                activeTracers(idxTemperature,:,:) = interpTracerField % array(:,:)

                call ocn_init_interpolation_nearest_horiz(tracerLon % array, tracerLat % array, &
                                                          salinityIC % array, nLonTracer, nLatTracer, &
                                                          lonCell, latCell, interpTracerField % array, nCells, &
                                                          inXPeriod = 2.0_RKIND * pii)
                activeTracers(idxSalinity,:,:) = interpTracerField % array(:,:)

             elseif (config_global_ocean_tracer_method .eq. "bilinear_interpolation") then
                call ocn_init_interpolation_bilinear_horiz(tracerLon % array, tracerLat % array, &
                                                          temperatureIC % array, nLonTracer, nLatTracer, &
                                                          lonCell, latCell, interpTracerField % array, nCells, &
                                                          inXPeriod = 2.0_RKIND * pii)
                activeTracers(idxTemperature,:,:) = interpTracerField % array(:,:)

                call ocn_init_interpolation_bilinear_horiz(tracerLon % array, tracerLat % array, &
                                                          salinityIC % array, nLonTracer, nLatTracer, &
                                                          lonCell, latCell, interpTracerField % array, nCells, &
                                                          inXPeriod = 2.0_RKIND * pii)
                activeTracers(idxSalinity,:,:) = interpTracerField % array(:,:)

             else
                write(stderrUnit,*) 'ERROR: Invalid choice of config_global_ocean_tracer_method.'
                iErr = 1
                call mpas_dmpar_finalize(domain % dminfo)
             endif
             call mpas_deallocate_scratch_field(interpTracerField, .true.)
          end if

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

             if ( activeTracersField % isActive ) then
                call mpas_dmpar_exch_halo_field(activeTracersField)
             end if

          end do ! iSmooth

          call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
          call mpas_pool_get_field(scratchPool, 'smoothedTemperature', smoothedTemperatureField)
          call mpas_pool_get_field(scratchPool, 'smoothedSalinity', smoothedSalinityField)
          call mpas_deallocate_scratch_field(smoothedTemperatureField, .false.)
          call mpas_deallocate_scratch_field(smoothedSalinityField, .false.)
       endif

       ! Interpolate tracers for partial bottom cells.
       ! This can be removed once we interpolate to an arbitrary vertical grid.
       call mpas_pool_get_config(domain % configs, 'config_alter_ICs_for_pbcs', config_alter_ICs_for_pbcs)
       call mpas_pool_get_config(domain % configs, 'config_pbc_alteration_type', config_pbc_alteration_type)
       if (config_alter_ICs_for_pbcs.and.config_pbc_alteration_type .eq. 'partial_cell') then

          block_ptr => domain % blocklist
          do while(associated(block_ptr))
             call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
             call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
             call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

             call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
             call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

             call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
             call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
             call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)

             call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)

             allocate(zMidZLevel(nVertLevels))
             zMidZLevel(1) = - 0.5*(refBottomDepth(1)) ! could add SSH here
             do k = 2, nVertLevels
                zMidZLevel(k) = - 0.5*(refBottomDepth(k) + refBottomDepth(k-1))
             end do

             call mpas_pool_begin_iteration(tracersPool)
             do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
                if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                   call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroup, 1)


                   if ( associated(tracersGroup) ) then
                      do iCell = 1, nCells
                         ! Linearly interpolate the initial tracers for new location of bottom cell for PBCs
                         k = maxLevelCell(iCell)
                         if (k>1) then
                            zMidPBC = -0.5_RKIND * (bottomDepth(iCell) + refBottomDepth(k-1))
                            km1 = max(k-1,1)
                            tracersGroup(:, k, iCell) = tracersGroup(:, k, iCell) &
                                 + (tracersGroup(:, km1, iCell) - tracersGroup(:, k, iCell)) &
                                 /(zMidZLevel(km1) - zMidZLevel(k) + 1.0e-16_RKIND) &
                                 *(zMidPBC - zMidZLevel(k))
                         endif
                      end do
                   end if
                end if
             end do

             deallocate(zMidZLevel)
             block_ptr => block_ptr % next
          end do
       endif
       ! end: Interpolate tracers for partial bottom cells.

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
!  routine ocn_init_setup_global_ocean_interoplate_swData
!
!> \brief   Interpolate penetrating shortwave radiation quantities to MPAS grid
!> \author  Luke Van Roekel
!> \date    11/11/2015
!> \details 
!>  This routine interpolates the penetrating swData data read in from the
!>  initial condition file to the MPAS grid. Currently it uses a nearest neighbor interpolation.
!
!-----------------------------------------------------------------------

subroutine ocn_init_setup_global_ocean_interpolate_swData(domain, iErr)!{{{
   type (domain_type), intent(inout) :: domain
   integer, intent(out) :: iErr

   type (block_type), pointer :: block_ptr
   type (MPAS_Stream_type) :: zenithStream, chlorophyllStream, clearSkyStream
   type (mpas_pool_type), pointer :: meshPool, statePool, scratchPool, shortwavePool, diagnosticsPool

   type(MPAS_TimeInterval_type) :: timeStep ! time step interval
   type(MPAS_Time_Type) :: currentTime
   character(len=STRKIND) :: currentTimeStamp
   
   real (kind=RKIND) :: currentLat, currentLon, counter
   real (kind=RKIND) :: minDist, dist, dt
   real (kind=RKIND) :: x, x1, x2, y, y1, y2, coef, coef11, coef12, coef21, coef22
   
   integer :: iLat, iLon, j, coc
   integer :: latSearch, lonSearch, timelen
   integer :: iCell, monIndex
   integer :: xInd1, xInd2, yInd1, yInd2
   integer, pointer :: nCells, nCellsSolve, maxLevelCell

   integer, dimension(12), parameter :: daysInMonth     = (/31,28,31,30,31,30,31,31,30,31,30,31/)
   integer, dimension(:), pointer :: nEdgesOnCell
   integer, dimension(:, :), pointer :: cellsOnCell

   real (kind=RKIND), dimension(:), pointer :: chlorophyllData, zenithAngle, clearSkyRadiation
   real (kind=RKIND), dimension(:), pointer :: latCell, lonCell

   character (len=StrKIND), pointer :: config_global_ocean_swData_method, xtime
   character (len=StrKIND), pointer :: config_global_ocean_swData_file, config_global_ocean_zenithAngle_varname, &
                                       config_global_ocean_swData_nlon_dimname, config_global_ocean_swData_nlat_dimname, &
                                       config_global_ocean_chlorophyll_varname, config_global_ocean_clearSky_varname

   iErr = 0

   call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_method', config_global_ocean_swData_method)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_file', config_global_ocean_swData_file)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_zenithAngle_varname', config_global_ocean_zenithAngle_varname)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_chlorophyll_varname', config_global_ocean_chlorophyll_varname)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_clearSky_varname', config_global_ocean_clearSky_varname)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_nlon_dimname', config_global_ocean_swData_nlon_dimname)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_nlat_dimname', config_global_ocean_swData_nlat_dimname)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_method', config_global_ocean_swData_method)
   ! Define stream for zenithAngle IC
   call MPAS_createStream(zenithStream, domain % iocontext, config_global_ocean_swData_file, MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)
   call MPAS_createStream(chlorophyllStream, domain % iocontext, config_global_ocean_swData_file, MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)
   call MPAS_createStream(clearSkyStream, domain % iocontext, config_global_ocean_swData_file, MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)
   

   ! Setup zenithAngle field for stream to be read in
   
   zenithAngleIC % fieldName = trim(config_global_ocean_zenithAngle_varname)
   zenithAngleIC % dimSizes(1) = nLonSW
   zenithAngleIC % dimSizes(2) = nLatSW
   zenithAngleIC % dimNames(1) = trim(config_global_ocean_swData_nlon_dimname)
   zenithAngleIC % dimNames(2) = trim(config_global_ocean_swData_nlat_dimname)
 
   zenithAngleIC % isVarArray = .false.
   zenithAngleIC % isPersistent = .true.
   zenithAngleIC % isActive = .true.
   zenithAngleIC % hasTimeDimension = .false.
   zenithAngleIC % block => domain % blocklist
   allocate(zenithAngleIC % array(nLonSW, nLatSW))
      
   ! Setup zenithAngle field for stream to be read in
   chlorophyllIC % fieldName = trim(config_global_ocean_chlorophyll_varname)
   chlorophyllIC % dimSizes(1) = nLonSW
   chlorophyllIC % dimSizes(2) = nLatSW
   chlorophyllIC % dimNames(1) = trim(config_global_ocean_swData_nlon_dimname)
   chlorophyllIC % dimNames(2) = trim(config_global_ocean_swData_nlat_dimname)
   chlorophyllIC % isVarArray = .false.
   chlorophyllIC % isPersistent = .true.
   chlorophyllIC % isActive = .true.
   chlorophyllIC % hasTimeDimension = .false.
   chlorophyllIC % block => domain % blocklist
   allocate(chlorophyllIC % array(nLonSW, nLatSW))

   ! Setup zenithAngle field for stream to be read in
   clearSKYIC % fieldName = trim(config_global_ocean_clearSky_varname)
   clearSKYIC % dimSizes(1) = nLonSW
   clearSKYIC % dimSizes(2) = nLatSW
   clearSKYIC % dimNames(1) = trim(config_global_ocean_swData_nlon_dimname)
   clearSKYIC % dimNames(2) = trim(config_global_ocean_swData_nlat_dimname)
   clearSKYIC % isVarArray = .false.
   clearSKYIC % isPersistent = .true.
   clearSKYIC % isActive = .true.
   clearSKYIC % hasTimeDimension = .false.
   clearSKYIC % block => domain % blocklist
   allocate(clearSKYIC % array(nLonSW, nLatSW))
   ! Add chlorophyll field to stream
      
   call MPAS_streamAddField(zenithStream, zenithAngleIC, iErr)
   call MPAS_streamAddField(chlorophyllStream, chlorophyllIC, iErr)
   call MPAS_streamAddField(clearSkyStream, clearSKYIC, iErr)
     
   do monIndex=1,12
   block_ptr => domain % blocklist
   do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
      call mpas_pool_get_subpool(block_ptr % structs, 'shortwave', shortwavePool)
      call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
      
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      call mpas_pool_get_array(shortwavePool, 'chlorophyllData', chlorophyllData)
      call mpas_pool_get_array(shortWavePool, 'zenithAngle', zenithAngle)
      call mpas_pool_get_array(shortWavePool, 'clearSkyRadiation', clearSkyRadiation)

   ! Read stream
      
      call MPAS_readStream(zenithStream, monIndex, iErr)
      call MPAS_readStream(chlorophyllStream, monIndex, iErr)
      call MPAS_readStream(clearSkyStream, monIndex, iErr)
     
      
      if (config_global_ocean_swData_method .eq. "nearest_neighbor") then
         call ocn_init_interpolation_nearest_horiz(swDataLon % array, swDataLat % array, &
                                                   chlorophyllIC % array, nLonSW, nLatSW, &
                                                   lonCell, latCell, chlorophyllData, nCells, &
                                                   inXPeriod = 2.0_RKIND * pii)

         call ocn_init_interpolation_nearest_horiz(swDataLon % array, swDataLat % array, &
                                                   zenithAngleIC % array, nLonSW, nLatSW, &
                                                   lonCell, latCell, zenithAngle, nCells, &
                                                   inXPeriod = 2.0_RKIND * pii)

         call ocn_init_interpolation_nearest_horiz(swDataLon % array, swDataLat % array, &
                                                   clearSKYIC % array, nLonSW, nLatSW, &
                                                   lonCell, latCell, clearSkyRadiation, nCells, &
                                                   inXPeriod = 2.0_RKIND * pii)

      elseif (config_global_ocean_swData_method .eq. "bilinear_interpolation") then
         call ocn_init_interpolation_bilinear_horiz(swDataLon % array, swDataLat % array, &
                                                    chlorophyllIC % array, nLonSW, nLatSW, &
                                                    lonCell, latCell, chlorophyllData, nCells, &
                                                    inXPeriod = 2.0_RKIND * pii)
 
         call ocn_init_interpolation_bilinear_horiz(swDataLon % array, swDataLat % array, &
                                                    zenithAngleIC % array, nLonSW, nLatSW, &
                                                    lonCell, latCell, zenithAngle, nCells, &
                                                    inXPeriod = 2.0_RKIND * pii)

         call ocn_init_interpolation_bilinear_horiz(swDataLon % array, swDataLat % array, &
                                                    clearSKYIC % array, nLonSW, nLatSW, &
                                                    lonCell, latCell, clearSkyRadiation, nCells, &
                                                    inXPeriod = 2.0_RKIND * pii)
     else
         write(stderrUnit,*) 'ERROR: Invalid choice of config_global_ocean_swData_method.'
         iErr = 1
         call mpas_dmpar_finalize(domain % dminfo)
      endif

      block_ptr => block_ptr % next
   end do !loop on blocks
 
  ! increment clock with month string

   currentTime = mpas_get_clock_time(domain % clock, MPAS_NOW, iErr)
   call mpas_get_time(currentTime, dateTimeString=currentTimeStamp)

   xtime=currentTimeStamp
   call mpas_stream_mgr_write(domain % streamManager, streamID='shortwave_forcing_data_init', &
                              forceWriteNow=.true., ierr=ierr)
   call mpas_set_timeInterval(timeStep, dt=real(daysInMonth(monIndex),RKIND)*86400.0_RKIND)
   call mpas_advance_clock(domain % clock, timeStep)
   
   enddo  !ends loop over months
   
   
   ! Close stream
   call MPAS_closeStream(zenithStream)
   call MPAS_closeStream(chlorophyllStream)
   call MPAS_closeStream(clearSkyStream)

! reset mpas clock for other streams and final write

   currentTime = mpas_get_clock_time(domain % clock, MPAS_START_TIME, iErr)
   call mpas_set_clock_time(domain%clock, currentTime  , MPAS_NOW,iErr)
   currentTime = mpas_get_clock_time(domain % clock, MPAS_NOW, iErr)
   call mpas_get_time(currentTime, dateTimeString=currentTimeStamp)

   xtime=currentTimeStamp
   call mpas_stream_mgr_reset_alarms(domain%streamManager, streamID='shortwave_forcing_data_init', &
                                     direction=MPAS_STREAM_OUTPUT, ierr=ierr)

end subroutine ocn_init_setup_global_ocean_interpolate_swData!}}}

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
       real (kind=RKIND) :: dist, minDist
       real (kind=RKIND) :: x, x1, x2, y, y1, y2, coef, coef11, coef12, coef21, coef22

       integer :: ilat, iLon
       integer :: latSearch, lonSearch
       integer :: iCell
       integer :: xInd1, xInd2, yInd1, yInd2

       real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, windStressZonal, windStressMeridional

       integer, pointer :: nCellsSolve, nCells

       logical, pointer :: config_use_bulk_wind_stress
       character (len=StrKIND), pointer :: config_global_ocean_windstress_method
       real (kind=RKIND), pointer :: config_global_ocean_windstress_conversion_factor

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_use_bulk_wind_stress', config_use_bulk_wind_stress)
       if (.not.config_use_bulk_wind_stress) then
          write(stderrUnit,'(A)') ' WARNING: wind stress not initialized because config_use_bulk_wind_stress = .false.'
          return
       endif

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_method', config_global_ocean_windstress_method)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_conversion_factor', config_global_ocean_windstress_conversion_factor)

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

          call mpas_pool_get_array(meshPool, 'latCell', latCell)
          call mpas_pool_get_array(meshPool, 'lonCell', lonCell)

          call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
          call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)

          if (config_global_ocean_windstress_method .eq. "nearest_neighbor") then
             call ocn_init_interpolation_nearest_horiz(windLon % array, windLat % array, &
                                                       zonalWindIC % array, nLonWind, nLatWind, &
                                                       lonCell, latCell, windStressZonal, nCells, &
                                                       inXPeriod = 2.0_RKIND * pii)

             call ocn_init_interpolation_nearest_horiz(windLon % array, windLat % array, &
                                                       meridionalWindIC % array, nLonWind, nLatWind, &
                                                       lonCell, latCell, windStressMeridional, nCells, &
                                                       inXPeriod = 2.0_RKIND * pii)

          elseif (config_global_ocean_windstress_method .eq. "bilinear_interpolation") then
             call ocn_init_interpolation_bilinear_horiz(windLon % array, windLat % array, &
                                                       zonalWindIC % array, nLonWind, nLatWind, &
                                                       lonCell, latCell, windStressZonal, nCells, &
                                                       inXPeriod = 2.0_RKIND*pii)

             call ocn_init_interpolation_bilinear_horiz(windLon % array, windLat % array, &
                                                       meridionalWindIC % array, nLonWind, nLatWind, &
                                                       lonCell, latCell, windStressMeridional, nCells, &
                                                       inXPeriod = 2.0_RKIND*pii)

          else
             write(stderrUnit,*) 'ERROR: Invalid choice of config_global_ocean_windstress_method.'
             iErr = 1
             call mpas_dmpar_finalize(domain % dminfo)
          endif

          windStressZonal(:) = windStressZonal(:) * config_global_ocean_windstress_conversion_factor
          windStressMeridional(:) = windStressMeridional(:) * config_global_ocean_windstress_conversion_factor

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
!  routine ocn_init_global_ocean_destroy_swData_fields
!
!> \brief   penetrating shortwave data fields cleanup routine
!> \author  Luke Van Roekel
!> \date    11/11/2015
!> \details 
!>  This routine destroys the fields that were created to hold penetrating sw radiation data
!>  initial condition information
!
!-----------------------------------------------------------------------

    subroutine ocn_init_global_ocean_destroy_swData_fields()!{{{
        deallocate(chlorophyllIC % array)
        deallocate(zenithAngleIC % array)
        deallocate(clearSKYIC % array)
    end subroutine ocn_init_global_ocean_destroy_swData_fields!}}}

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

      integer, intent(out) :: iErr
      type (MPAS_IO_Handle_type) :: inputFile

      character (len=StrKIND), pointer :: config_init_configuration, config_global_ocean_depth_file, &
                                          config_global_ocean_depth_dimname, config_global_ocean_temperature_file, &
                                          config_global_ocean_salinity_file, config_global_ocean_tracer_nlat_dimname, &
                                          config_global_ocean_tracer_nlon_dimname, config_global_ocean_topography_file, &
                                          config_global_ocean_topography_nlat_dimname, config_global_ocean_topography_nlon_dimname, &
                                          config_global_ocean_windstress_file, config_global_ocean_windstress_nlat_dimname, &
                                          config_global_ocean_windstress_nlon_dimname, config_global_ocean_swData_file, &
                                          config_global_ocean_swData_nlon_dimname, config_global_ocean_swData_nlat_dimname

      integer, pointer :: config_vert_levels
      type (mpas_io_context_type), pointer :: iocontext_ptr

      iocontext_ptr => iocontext

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
      call mpas_pool_get_config(configPool, 'config_global_ocean_swData_file', config_global_ocean_swData_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_swData_nlat_dimname', config_global_ocean_swData_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_swData_nlon_dimname', config_global_ocean_swData_nlon_dimname)


      inputFile = MPAS_io_open(config_global_ocean_depth_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)

      call MPAS_io_inq_dim(inputFile, config_global_ocean_depth_dimname, nDepth, iErr)

      call MPAS_io_close(inputFile, iErr)

      inputFile = MPAS_io_open(config_global_ocean_temperature_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)

      call MPAS_io_inq_dim(inputFile, config_global_ocean_tracer_nlat_dimname, nLatTracer, iErr)
      call MPAS_io_inq_dim(inputFile, config_global_ocean_tracer_nlon_dimname, nLonTracer, iErr)

      call MPAS_io_close(inputFile, iErr)

      inputFile = MPAS_io_open(config_global_ocean_swData_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)
 
      call MPAS_io_inq_dim(inputFile, config_global_ocean_swData_nlat_dimname, nLatSW, iErr)
      call MPAS_io_inq_dim(inputFile, config_global_ocean_swData_nlon_dimname, nLonSW, iErr)

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
