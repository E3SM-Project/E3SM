! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_isomip_plus
!
!> \brief MPAS ocean initialize case -- ISOMIP+
!> \author Xylar Asay-Davis
!> \date   07/21/2015
!> \details
!>  This module contains the routines for initializing the
!>  the ISOMIP+ experiments
!
!-----------------------------------------------------------------------

module ocn_init_isomip_plus

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_io
   use mpas_io_streams
   use mpas_dmpar

   use ocn_constants

   use ocn_init_cell_markers

   use ocn_init_vertical_grids
   use ocn_init_ssh_and_landIcePressure
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

   public :: ocn_init_setup_isomip_plus, &
             ocn_init_validate_isomip_plus

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: nyTopo, nxTopo
   type (field1DReal) :: topoY, topoX
   type (field2DReal) :: iceTopoIC, bedTopoIC, oceanFracTopoIC, &
                         floatingFracTopoIC, smoothedDraftMaskTopoIC

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_isomip_plus
!
!> \brief   Setup for ISOMIP+ experiments
!> \author  Xylar Asay-Davis
!> \date    07/21/2015
!> \details
!>  This routine sets up the initial conditions for the ISOMIP+ experiments
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_isomip_plus(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      type (mpas_pool_type), pointer :: meshPool, statePool, diagnosticsPool, &
                                        verticalMeshPool, forcingPool, &
                                        tracersPool
      type (mpas_pool_type), pointer :: tracersInteriorRestoringFieldsPool

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration, &
                                          config_isomip_plus_vertical_level_distribution, &
                                          config_land_ice_flux_mode

      real(kind=RKIND), pointer :: config_isomip_plus_max_bottom_depth, &
                                   config_isomip_plus_init_top_temp, &
                                   config_isomip_plus_init_bot_temp, &
                                   config_isomip_plus_init_top_sal, &
                                   config_isomip_plus_init_bot_sal, &
                                   config_isomip_plus_restore_top_temp, &
                                   config_isomip_plus_restore_bot_temp, &
                                   config_isomip_plus_restore_top_sal, &
                                   config_isomip_plus_restore_bot_sal, &
                                   config_isomip_plus_restore_rate, &
                                   config_isomip_plus_restore_evap_rate, &
                                   config_isomip_plus_restore_xMin, &
                                   config_isomip_plus_restore_xMax, &
                                   config_isomip_plus_coriolis_parameter, &
                                   config_isomip_plus_effective_density

      type (block_type), pointer :: block_ptr

      logical, pointer :: on_a_sphere

      integer, pointer :: nCells, nVertLevels
      integer, pointer :: index_temperature, index_salinity, index_tracer1

      integer, dimension(:), pointer :: maxLevelCell

      real(kind=RKIND), dimension(:), pointer :: refBottomDepth, &
                                                 fCell, fEdge, fVertex, effectiveDensityInLandIce, xCell, &
                                                 refZMid, refLayerThickness

      real(kind=RKIND), dimension(:,:), pointer :: zMid
      real(kind=RKIND), dimension(:,:), pointer :: layerThickness
      real(kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      real(kind=RKIND), dimension(:,:,:), pointer :: debugTracers

      real(kind=RKIND), dimension(:,:,:), pointer :: activeTracersInteriorRestoringValue, activeTracersInteriorRestoringRate
      real (kind=RKIND), dimension(:), pointer :: evaporationFlux, seaIceSalinityFlux, seaIceHeatFlux


      real(kind=RKIND), dimension(:), pointer :: columnThicknessFraction

      integer :: iCell, k
      real(kind=RKIND) :: frac, evapRate

      iErr = 0

      call mpas_pool_get_config(domain % configs, 'config_init_configuration', config_init_configuration)

      if (trim(config_init_configuration) /= "isomip_plus") return

      call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      if ( on_a_sphere ) call mpas_dmpar_global_abort('ERROR: ISOMIP+ can only be applied to a planar mesh.  Exiting...')

      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_vertical_level_distribution', &
         config_isomip_plus_vertical_level_distribution)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_max_bottom_depth', config_isomip_plus_max_bottom_depth)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_init_top_temp', config_isomip_plus_init_top_temp)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_init_bot_temp', config_isomip_plus_init_bot_temp)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_init_top_sal', config_isomip_plus_init_top_sal)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_init_bot_sal', config_isomip_plus_init_bot_sal)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_restore_top_temp', config_isomip_plus_restore_top_temp)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_restore_bot_temp', config_isomip_plus_restore_bot_temp)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_restore_top_sal', config_isomip_plus_restore_top_sal)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_restore_bot_sal', config_isomip_plus_restore_bot_sal)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_restore_rate', config_isomip_plus_restore_rate)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_restore_evap_rate', config_isomip_plus_restore_evap_rate)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_restore_xMin', config_isomip_plus_restore_xMin)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_restore_xMax', config_isomip_plus_restore_xMax)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_coriolis_parameter', config_isomip_plus_coriolis_parameter)
      call mpas_pool_get_config(domain % configs, 'config_isomip_plus_effective_density', config_isomip_plus_effective_density)


      allocate(columnThicknessFraction(nVertLevels))
      if(trim(config_isomip_plus_vertical_level_distribution) == "constant") then
         columnThicknessFraction(:) = 1.0_RKIND/nVertLevels
      end if

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)
        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)

        ! Set refBottomDepth
        refBottomDepth(1) =  columnThicknessFraction(1)*abs(config_isomip_plus_max_bottom_depth)
        do k = 2, nVertLevels
            refBottomDepth(k) = refBottomDepth(k-1) + columnThicknessFraction(k)*abs(config_isomip_plus_max_bottom_depth)
        end do

        ! Compute refLayerThickness and refZMid
        call ocn_compute_layerThickness_zMid_from_bottomDepth(refLayerThickness,refZMid, &
             refBottomDepth,refBottomDepth(nVertLevels), &
             nVertLevels,nVertLevels,iErr)

        block_ptr => block_ptr % next
      end do

      write(stderrUnit,*) 'Reading topography data.'
      call ocn_init_setup_isomip_plus_read_topo(domain, iErr)
      write(stderrUnit,*) 'Interpolating topography data.'
      call ocn_init_setup_isomip_plus_interpolate_topo(domain, iErr)
      write(stderrUnit,*) 'Cleaning up topography IC fields'
      call ocn_init_isomip_plus_destroy_topo_fields()

      ! compute the vertical grid (layerThickness, restingThickness, maxLevelCell, zMid) bottomDepth and refBottomDepth
      call ocn_init_ssh_and_landIcePressure_vertical_grid(domain, iErr)

      if(iErr .ne. 0) then
        write(stderrUnit,*) 'ERROR: ocn_init_ssh_and_landIcePressure_vertical_grid failed.'
        return
      end if

      block_ptr => domain % blocklist
      do while (associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)

        do iCell = 1, nCells
          do k = 1, maxLevelCell(iCell)
            frac = (0.0_RKIND-zMid(k, iCell))/(0.0_RKIND-config_isomip_plus_max_bottom_depth)
            activeTracers(index_temperature, k, iCell) = (1.0_RKIND - frac)*config_isomip_plus_init_top_temp &
                                                 + frac*config_isomip_plus_init_bot_temp
            activeTracers(index_salinity, k, iCell) = (1.0_RKIND - frac)*config_isomip_plus_init_top_sal &
                                              + frac*config_isomip_plus_init_bot_sal
          end do
        end do

        block_ptr => block_ptr % next
      end do


      ! compute or update the land-ice pressure (or possibly SSH), also computing density and bottomPressure along the way
      ! If this is the initial guess, the vertical grid and activeTracers may also be recomputed based on SSH
      call ocn_init_ssh_and_landIcePressure_balance(domain, iErr)

      if(iErr .ne. 0) then
        write(stderrUnit,*) 'ERROR: ocn_init_ssh_and_landIcePressure_balance failed.'
        return
      end if


      block_ptr => domain % blocklist
      do while (associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'fCell', fCell)
        call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
        call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

        call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
        call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
        call mpas_pool_get_dimension(tracersPool, 'index_tracer1', index_tracer1)

        call mpas_pool_get_subpool(forcingPool, 'tracersInteriorRestoringFields', tracersInteriorRestoringFieldsPool)
        call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, 'activeTracersInteriorRestoringRate', &
                activeTracersInteriorRestoringRate, 1)
        call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, 'activeTracersInteriorRestoringValue', &
                activeTracersInteriorRestoringValue, 1)

        call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
        call mpas_pool_get_array(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFlux)
        call mpas_pool_get_array(forcingPool, 'seaIceHeatFlux', seaIceHeatFlux)

        if(associated(activeTracersInteriorRestoringValue)) then
          activeTracersInteriorRestoringRate(:,:,:) = 0.0_RKIND
          activeTracersInteriorRestoringValue(:,:,:) = 0.0_RKIND
          do iCell = 1, nCells
            do k = 1, maxLevelCell(iCell)
              frac = (0.0_RKIND-zMid(k, iCell))/(0.0_RKIND-config_isomip_plus_max_bottom_depth)
              activeTracersInteriorRestoringValue(index_temperature, k, iCell) &
                = (1.0_RKIND - frac)*config_isomip_plus_restore_top_temp &
                + frac*config_isomip_plus_restore_bot_temp
              activeTracersInteriorRestoringValue(index_salinity, k, iCell) &
                = (1.0_RKIND - frac)*config_isomip_plus_restore_top_sal &
                + frac*config_isomip_plus_restore_bot_sal
            end do
            ! fraction of the way between southern and northern boundaries of restoring region
            frac = max((xCell(iCell)-config_isomip_plus_restore_xMin) &
              /(config_isomip_plus_restore_xMax-config_isomip_plus_restore_xMin), 0.0_RKIND)
            ! restoring rate converted from 1/days to 1/s
            activeTracersInteriorRestoringRate(:,:,iCell) = frac*config_isomip_plus_restore_rate/86400.0_RKIND
          end do
        end if

        if(associated(evaporationFlux) .and. associated(seaIceSalinityFlux) &
            .and. associated(seaIceHeatFlux)) then
          do iCell = 1, nCells
            if((xCell(iCell) >= config_isomip_plus_restore_xMin)  &
                .and. (xCell(iCell) <= config_isomip_plus_restore_xMax)) then
              ! convert to m/s, negative for evaporation rather than precipitation
              evapRate = -config_isomip_plus_restore_evap_rate/31536000.0_RKIND
              evaporationFlux(iCell) = rho_sw*evapRate
              ! associated salt flux in kg/m^2/s of salt
              seaIceSalinityFlux(iCell) = evapRate*config_isomip_plus_restore_top_sal/sflux_factor
              ! associated heat flux in W/m^2
              seaIceHeatFlux(iCell) = evapRate*config_isomip_plus_restore_top_temp/hflux_factor
            else
              evaporationFlux(iCell) = 0.0_RKIND
              seaIceSalinityFlux(iCell) = 0.0_RKIND
              seaIceHeatFlux(iCell) = 0.0_RKIND
            end if
          end do
        end if

        if(associated(debugTracers)) then
          debugTracers(index_tracer1,:,:) = 1.0_RKIND
        end if

        call ocn_mark_maxlevelcell(meshPool, iErr)

        ! Set Coriolis parameters
        fCell(:) = config_isomip_plus_coriolis_parameter
        fEdge(:) = config_isomip_plus_coriolis_parameter
        fVertex(:) = config_isomip_plus_coriolis_parameter

        if(config_land_ice_flux_mode == 'coupled') then
          call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', effectiveDensityInLandIce, 1)
          effectiveDensityInLandIce(:) = config_isomip_plus_effective_density
        end if

        block_ptr => block_ptr % next
      end do

      call ocn_compute_Haney_number(domain, iErr)

      if(iErr .ne. 0) then
        write(stderrUnit,*) 'ERROR: ocn_compute_Haney_number failed.'
        return
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_isomip_plus!}}}

!***********************************************************************
!
!  routine ocn_init_setup_isomip_plus_read_topo
!
!> \brief   Read the topography IC file
!> \author  Xylar Asay-Davis
!> \date    07/21/2015
!> \details
!>  This routine reads the topography IC file, including x and y coordinates.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_isomip_plus_read_topo(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (MPAS_Stream_type) :: topographyStream

       character (len=StrKIND), pointer :: config_isomip_plus_topography_file

       integer :: iY, iX

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_isomip_plus_topography_file', config_isomip_plus_topography_file)

       ! Define stream for depth levels
       call MPAS_createStream(topographyStream, domain % iocontext, config_isomip_plus_topography_file, &
          MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)

       ! Setup topoY, topoX, bedTopoIC and iceTopoIC fields for stream to be read in
       topoY % fieldName = 'y'
       topoY % dimSizes(1) = nyTopo
       topoY % dimNames(1) = 'y'
       topoY % isVarArray = .false.
       topoY % isPersistent = .true.
       topoY % isActive = .true.
       topoY % hasTimeDimension = .false.
       topoY % block => domain % blocklist
       allocate(topoY % attLists(1))
       allocate(topoY % array(nyTopo))

       topoX % fieldName = 'x'
       topoX % dimSizes(1) = nxTopo
       topoX % dimNames(1) = 'x'
       topoX % isVarArray = .false.
       topoX % isPersistent = .true.
       topoX % isActive = .true.
       topoX % hasTimeDimension = .false.
       topoX % block => domain % blocklist
       allocate(topoX % attLists(1))
       allocate(topoX % array(nxTopo))

       bedTopoIC % fieldName = 'Z_bed'
       bedTopoIC % dimSizes(1) = nxTopo
       bedTopoIC % dimSizes(2) = nyTopo
       bedTopoIC % dimNames(1) = 'x'
       bedTopoIC % dimNames(2) = 'y'
       bedTopoIC % isVarArray = .false.
       bedTopoIC % isPersistent = .true.
       bedTopoIC % isActive = .true.
       bedTopoIC % hasTimeDimension = .false.
       bedTopoIC % block => domain % blocklist
       allocate(bedTopoIC % attLists(1))
       allocate(bedTopoIC % array(nxTopo, nyTopo))

       iceTopoIC % fieldName = 'Z_ice_draft'
       iceTopoIC % dimSizes(1) = nxTopo
       iceTopoIC % dimSizes(2) = nyTopo
       iceTopoIC % dimNames(1) = 'x'
       iceTopoIC % dimNames(2) = 'y'
       iceTopoIC % isVarArray = .false.
       iceTopoIC % isPersistent = .true.
       iceTopoIC % isActive = .true.
       iceTopoIC % hasTimeDimension = .false.
       iceTopoIC % block => domain % blocklist
       allocate(iceTopoIC % attLists(1))
       allocate(iceTopoIC % array(nxTopo, nyTopo))

       ! ocean fraction is 1.0 - landFraction
       oceanFracTopoIC % fieldName = 'landFraction'
       oceanFracTopoIC % dimSizes(1) = nxTopo
       oceanFracTopoIC % dimSizes(2) = nyTopo
       oceanFracTopoIC % dimNames(1) = 'x'
       oceanFracTopoIC % dimNames(2) = 'y'
       oceanFracTopoIC % isVarArray = .false.
       oceanFracTopoIC % isPersistent = .true.
       oceanFracTopoIC % isActive = .true.
       oceanFracTopoIC % hasTimeDimension = .false.
       oceanFracTopoIC % block => domain % blocklist
       allocate(oceanFracTopoIC % attLists(1))
       allocate(oceanFracTopoIC % array(nxTopo, nyTopo))

       floatingFracTopoIC % fieldName = 'floatingIceFraction'
       floatingFracTopoIC % dimSizes(1) = nxTopo
       floatingFracTopoIC % dimSizes(2) = nyTopo
       floatingFracTopoIC % dimNames(1) = 'x'
       floatingFracTopoIC % dimNames(2) = 'y'
       floatingFracTopoIC % isVarArray = .false.
       floatingFracTopoIC % isPersistent = .true.
       floatingFracTopoIC % isActive = .true.
       floatingFracTopoIC % hasTimeDimension = .false.
       floatingFracTopoIC % block => domain % blocklist
       allocate(floatingFracTopoIC % attLists(1))
       allocate(floatingFracTopoIC % array(nxTopo, nyTopo))

       smoothedDraftMaskTopoIC % fieldName = 'smoothedDraftMask'
       smoothedDraftMaskTopoIC % dimSizes(1) = nxTopo
       smoothedDraftMaskTopoIC % dimSizes(2) = nyTopo
       smoothedDraftMaskTopoIC % dimNames(1) = 'x'
       smoothedDraftMaskTopoIC % dimNames(2) = 'y'
       smoothedDraftMaskTopoIC % isVarArray = .false.
       smoothedDraftMaskTopoIC % isPersistent = .true.
       smoothedDraftMaskTopoIC % isActive = .true.
       smoothedDraftMaskTopoIC % hasTimeDimension = .false.
       smoothedDraftMaskTopoIC % block => domain % blocklist
       allocate(smoothedDraftMaskTopoIC % attLists(1))
       allocate(smoothedDraftMaskTopoIC % array(nxTopo, nyTopo))

       ! Add topoY, topoX, bedTopoIC and iceTopoIC fields to stream
       call MPAS_streamAddField(topographyStream, topoY, iErr)
       call MPAS_streamAddField(topographyStream, topoX, iErr)
       call MPAS_streamAddField(topographyStream, bedTopoIC, iErr)
       call MPAS_streamAddField(topographyStream, iceTopoIC, iErr)
       call MPAS_streamAddField(topographyStream, oceanFracTopoIC, iErr)
       call MPAS_streamAddField(topographyStream, floatingFracTopoIC, iErr)
       call MPAS_streamAddField(topographyStream, smoothedDraftMaskTopoIC, iErr)

       ! Read stream
       call MPAS_readStream(topographyStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(topographyStream)

       ! we read in the land fraction, and oceanFrac = 1.0 - landFraction
       oceanFracTopoIC % array(:,:) = 1.0_RKIND - oceanFracTopoIC % array(:,:)

    end subroutine ocn_init_setup_isomip_plus_read_topo!}}}

!***********************************************************************
!
!  routine ocn_init_setup_isomip_plus_interpolate_topo
!
!> \brief   Interpolate the topography IC to MPAS mesh
!> \author  Xylar Asay-Davis
!> \date    07/21/2015
!> \details
!>  This routine interpolates topography data to the MPAS mesh. Currently it
!>   uses a bilinear interpolation
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_isomip_plus_interpolate_topo(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (mpas_pool_type), pointer :: meshPool, verticalMeshPool, &
                                         forcingPool, diagnosticsPool, &
                                         statePool

       real (kind=RKIND), dimension(:), pointer :: yCell, xCell, bottomDepth, refBottomDepth, &
                                                   ssh, landIceFraction, oceanFracObserved, &
                                                   bottomDepthObserved

       integer, pointer :: nCells, nVertLevels

       integer, dimension(:), pointer :: maxLevelCell, modifySSHMask

       integer :: iCell, k, maxLevel

       integer, pointer :: config_isomip_plus_minimum_levels
       real(kind=RKIND), pointer :: config_isomip_plus_min_column_thickness, config_isomip_plus_min_ocean_fraction

       real (kind=RKIND) :: depth

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_isomip_plus_minimum_levels', &
          config_isomip_plus_minimum_levels)
       call mpas_pool_get_config(domain % configs, 'config_isomip_plus_min_column_thickness', &
          config_isomip_plus_min_column_thickness)
       call mpas_pool_get_config(domain % configs, 'config_isomip_plus_min_ocean_fraction', config_isomip_plus_min_ocean_fraction)

       ! weight bed topo and ice draft by the fraction of ocean before interpolating
       bedTopoIC % array(:,:) = bedTopoIC % array(:,:) * oceanFracTopoIC % array(:,:)
       iceTopoIC % array(:,:) = iceTopoIC % array(:,:) * oceanFracTopoIC % array(:,:)

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

          call mpas_pool_get_array(meshPool, 'yCell', yCell)
          call mpas_pool_get_array(meshPool, 'xCell', xCell)
          call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
          call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
          call mpas_pool_get_array(meshPool, 'bottomDepthObserved', bottomDepthObserved)
          call mpas_pool_get_array(meshPool, 'oceanFracObserved', oceanFracObserved)

          call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
          call mpas_pool_get_array(diagnosticsPool, 'modifySSHMask', modifySSHMask)
          call mpas_pool_get_array(forcingPool, 'landIceFraction', landIceFraction)


          ! use landIceFraction as temporary storage here
          call ocn_init_interpolation_bilinear_horiz(topoX % array, topoY % array, &
                                                     smoothedDraftMaskTopoIC % array, nxTopo, nyTopo, &
                                                     xCell, yCell, landIceFraction, nCells)

          do iCell = 1,nCells
            if(landIceFraction(iCell) > 0.01_RKIND) then
              modifySSHMask(iCell) = 1
            else
              modifySSHMask(iCell) = 0
            end if
          end do

          call ocn_init_interpolation_bilinear_horiz(topoX % array, topoY % array, &
                                                     oceanFracTopoIC % array, nxTopo, nyTopo, &
                                                     xCell, yCell, oceanFracObserved, nCells)

          call ocn_init_interpolation_bilinear_horiz(topoX % array, topoY % array, &
                                                     floatingFracTopoIC % array, nxTopo, nyTopo, &
                                                     xCell, yCell, landIceFraction, nCells)

          call ocn_init_interpolation_bilinear_horiz(topoX % array, topoY % array, &
                                                     bedTopoIC % array, nxTopo, nyTopo, &
                                                     xCell, yCell, bottomDepthObserved, nCells)

          bottomDepthObserved(:) = -bottomDepthObserved(:)

          call ocn_init_interpolation_bilinear_horiz(topoX % array, topoY % array, &
                                                     iceTopoIC % array, nxTopo, nyTopo, &
                                                     xCell, yCell, ssh, nCells)

          do iCell = 1,nCells
            if(oceanFracObserved(iCell) > config_isomip_plus_min_ocean_fraction) then
              ! we multiplied by the fraction of ocean before; divide by the interpolated ocean fraction
              bottomDepthObserved(iCell) = bottomDepthObserved(iCell)/oceanFracObserved(iCell)
            else
              bottomDepthObserved(iCell) = 0.0_RKIND
            end if

            if(bottomDepthObserved(iCell) > 0.0_RKIND) then
              ! we multiplied by the fraction of ocean before; divide by the interpolated ocean fraction
              ssh(iCell) = ssh(iCell)/oceanFracObserved(iCell)
              ! bottomDepth is lowered (increased) if necessary to accommodate the minimum column thickness
              bottomDepth(iCell) = max(bottomDepthObserved(iCell),-ssh(iCell)+config_isomip_plus_min_column_thickness)
              ! landIceFraction is the fraction of the *ocean* that is covered by land ice, so divide
              ! by the fraction of the cell that is ocean
              landIceFraction(iCell) = landIceFraction(iCell)/oceanFracObserved(iCell)

              if(landIceFraction(iCell) < 0.5_RKIND) then
                 landIceFraction(iCell) = 0.0_RKIND
              end if

              maxLevelCell(iCell) = -1
              do k = 1, nVertLevels
                depth = refBottomDepth(k)

                if (depth >= bottomDepth(iCell)) then
                  maxLevelCell(iCell) = k
                  exit
                end if
              end do

              if (maxLevelCell(iCell) == -1) then
                maxLevelCell(iCell) = nVertLevels
                bottomDepth(iCell) = refBottomDepth( nVertLevels )
                ! this time we may need to raise the ssh if the column is too thin because we can't go any deeper
                ssh(iCell) = max(ssh(iCell), -bottomDepth(iCell)+config_isomip_plus_min_column_thickness)
              else if (maxLevelCell(iCell) <= config_isomip_plus_minimum_levels) then
                maxLevelCell(iCell) = config_isomip_plus_minimum_levels
                bottomDepth(iCell) = refBottomDepth( config_isomip_plus_minimum_levels )
              end if

            else ! not an ocean cell
              bottomDepth(iCell) = 0.0_RKIND
              maxLevelCell(iCell) = -1
              ssh(iCell) = 0.0_RKIND
              landIceFraction(iCell) = 0.0_RKIND
              modifySSHMask(iCell) = 0
            end if
          end do !iCell

          ! Enforce minimum number of layers in ocean cells.
          do iCell = 1, nCells
              if (maxLevelCell(iCell) > 0 .and. maxLevelCell(iCell) < config_isomip_plus_minimum_levels) then
                  maxLevelCell(iCell) = config_isomip_plus_minimum_levels
                  bottomDepth(iCell) = refBottomDepth(config_isomip_plus_minimum_levels)
              end if
          end do

          block_ptr => block_ptr % next
       end do

    end subroutine ocn_init_setup_isomip_plus_interpolate_topo!}}}

!***********************************************************************
!
!  routine ocn_init_isomip_plus_destroy_topo_fields
!
!> \brief   Topography field cleanup routine
!> \author  Xylar Asay-Davis
!> \date    07/21/2015
!> \details
!>  This routine destroys the fields that were created to hold topography
!>  initial condition information
!
!-----------------------------------------------------------------------

    subroutine ocn_init_isomip_plus_destroy_topo_fields()!{{{
        deallocate(bedTopoIC % array)
        deallocate(iceTopoIC % array)
        deallocate(oceanFracTopoIC % array)
        deallocate(floatingFracTopoIC % array)
        deallocate(topoY % array)
        deallocate(topoX % array)
    end subroutine ocn_init_isomip_plus_destroy_topo_fields!}}}

!***********************************************************************
!
!  routine ocn_init_validate_isomip_plus
!
!> \brief   Validation for ISOMIP+ experiments
!> \author  Xylar Asay-Davis
!> \date    07/21/2015
!> \details
!>  This routine validates the configuration options for the ISOMIP+ experiments.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_isomip_plus(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout), target :: iocontext
      integer, intent(out) :: iErr

      type (mpas_io_context_type), pointer :: iocontext_ptr

      type (MPAS_IO_Handle_type) :: inputFile

      character (len=StrKIND), pointer :: config_init_configuration, &
                                          config_isomip_plus_topography_file, &
                                          config_isomip_plus_vertical_level_distribution, &
                                          config_land_ice_flux_mode

      integer, pointer :: config_vert_levels, config_isomip_plus_vert_levels

      iErr = 0

      iocontext_ptr => iocontext

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('isomip_plus')) return


      call mpas_pool_get_config(configPool, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
      if(config_land_ice_flux_mode .ne. 'standalone' .and. config_land_ice_flux_mode .ne. 'coupled') then
         write(stderrUnit,*) 'ERROR: Validation failed for isomip_plus. config_land_ice_flux_mode must either',&
            ' be standalone or coupled.'
         iErr = 1
         return
      end if

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_isomip_plus_vert_levels', config_isomip_plus_vert_levels)

      if(config_vert_levels <= 0 .and. config_isomip_plus_vert_levels > 0) then
         config_vert_levels = config_isomip_plus_vert_levels
      else if (config_vert_levels <= 0) then
         write(stderrUnit,*) 'ERROR: Validation failed for isomip_plus. Not given a usable value for vertical levels.'
         iErr = 1
         return
      end if

      call mpas_pool_get_config(configPool, 'config_isomip_plus_vertical_level_distribution', &
         config_isomip_plus_vertical_level_distribution)

      if(trim(config_isomip_plus_vertical_level_distribution) .ne. "constant") then
         write(stderrUnit,*) 'ERROR: Validation failed for isomip_plus. Bad vertical level distribution.'
         iErr = 1
         return
      end if


      call mpas_pool_get_config(configPool, 'config_isomip_plus_topography_file', config_isomip_plus_topography_file)

      if (trim(config_isomip_plus_topography_file) == 'none') then
         write(stderrUnit,*) 'ERROR: Validation failed for ISOMIP+. Invalid filename for config_isomip_plus_topography_file'
         iErr = 1
         return
      end if

      inputFile = MPAS_io_open(config_isomip_plus_topography_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)
      if(iErr > 0) then
         write(stderrUnit,*) 'ERROR: Validation failed for isomip_plus. Could not open file:', &
            trim(config_isomip_plus_topography_file)
         return
      end if

      call MPAS_io_inq_dim(inputFile, 'x', nxTopo, iErr)
      if(iErr /= 0) then
         write(stderrUnit,*) 'ERROR: Validation failed for isomip_plus. Could not find dim x in file:', &
            trim(config_isomip_plus_topography_file)
         return
      end if

      call MPAS_io_inq_dim(inputFile, 'y', nyTopo, iErr)
      if(iErr /= 0) then
         write(stderrUnit,*) 'ERROR: Validation failed for isomip_plus. Could not find dim y in file:', &
            trim(config_isomip_plus_topography_file)
         return
      end if

      call MPAS_io_close(inputFile, iErr)
      if(iErr /= 0) then
         write(stderrUnit,*) 'ERROR: Validation failed for isomip_plus. Could not close file:', &
            trim(config_isomip_plus_topography_file)
         return
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_isomip_plus!}}}

!***********************************************************************

end module ocn_init_isomip_plus

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
