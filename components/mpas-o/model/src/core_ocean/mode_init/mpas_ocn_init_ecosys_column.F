! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_ecosys_column
!
!> \brief MPAS ocean initialize case -- BGC (ecosys + DMS + MacroMolecules)
!> \author Mathew Maltrud
!> \date   11/01/2015
!> \details
!>  This module contains the routines for initializing the
!>  the ecosys column test configuration. This in a
!>  single column configuration.
!
!-----------------------------------------------------------------------

module ocn_init_ecosys_column

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_io_streams

   use ocn_init_cell_markers
   use ocn_init_vertical_grids

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

   public :: ocn_init_setup_ecosys_column, &
             ocn_init_setup_ecosys_read_column, &
             ocn_init_validate_ecosys_column

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type (field2DReal) :: columnIC

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_ecosys_column
!
!> \brief   Setup for ecosys column test configuration
!> \author  Mathew Maltrud
!> \date    11/01/2015
!> \details
!>  This routine sets up the initial conditions for the ecosys column test configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_ecosys_column(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool, statePool
      type (mpas_pool_type), pointer :: diagnosticsPool, forcingPool
      type (mpas_pool_type), pointer :: ecosysAuxiliary  ! additional forcing fields

      type (mpas_pool_type), pointer :: tracersPool

      integer, pointer :: nVertLevels, nVertLevelsP1, nCellsSolve, index_dummy

      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, refZMid, vertCoordMovementWeights
      real (kind=RKIND), dimension(:), pointer :: bottomDepth
      real (kind=RKIND), dimension(:, :), pointer :: layerThickness, restingThickness
      real (kind=RKIND), dimension(:), pointer :: PH_PREV, PH_PREV_ALT_CO2
      real (kind=RKIND), dimension(:, :), pointer :: PH_PREV_3D, PH_PREV_ALT_CO2_3D
      real (kind=RKIND), dimension(:, :, :), pointer :: activeTracers, ecosysTracers, DMSTracers,  &
         MacroMoleculesTracers

      real (kind=RKIND), dimension(:), pointer :: interfaceLocations

      real (kind=RKIND), allocatable, dimension(:,:) :: ecoFieldColumn

      integer :: iCell, iEdge, iVertex, iField, k, numTracersTotal, nVertLevelsInputColumn

      integer, allocatable, dimension(:) :: indexField

      character (len=StrKIND) :: fieldName

      character (len=StrKIND), pointer :: config_init_configuration, &
                                          config_ecosys_column_TS_filename, &
                                          config_ecosys_column_ecosys_filename, &
                                          config_ecosys_column_vertical_grid

      integer, pointer ::                 config_ecosys_column_vert_levels

      real (kind=RKIND), pointer ::       config_ecosys_column_bottom_depth

      ! assume no error
      iErr = 0

      ! get and test if this is the configuration specified
      call mpas_pool_get_config(domain % configs, 'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('ecosys_column')) return

      ! build the vertical grid
      ! intent(out) is interfaceLocations. An array ranging from 0 to 1
      call mpas_pool_get_config(domain % configs, 'config_ecosys_column_vertical_grid', config_ecosys_column_vertical_grid)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevelsP1', nVertLevelsP1)
      allocate(interfaceLocations(nVertLevelsP1))
      call ocn_generate_vertical_grid(config_ecosys_column_vertical_grid, interfaceLocations)

      ! load the remaining configuration parameters
      call mpas_pool_get_config(domain % configs, 'config_ecosys_column_bottom_depth', config_ecosys_column_bottom_depth)
      call mpas_pool_get_config(domain % configs, 'config_ecosys_column_TS_filename', config_ecosys_column_TS_filename)
      call mpas_pool_get_config(domain % configs, 'config_ecosys_column_ecosys_filename', config_ecosys_column_ecosys_filename)
      call mpas_pool_get_config(domain % configs, 'config_ecosys_column_vert_levels', config_ecosys_column_vert_levels)

      nVertLevelsInputColumn = config_ecosys_column_vert_levels

      ! load data that required to initialize the ocean simulation
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
        call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, 1)
        call mpas_pool_get_array(tracersPool, 'DMSTracers', DMSTracers, 1)
        call mpas_pool_get_array(tracersPool, 'MacroMoleculesTracers', MacroMoleculesTracers, 1)
        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

        call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

        call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV', PH_PREV)
        call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV_ALT_CO2', PH_PREV_ALT_CO2)
        call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV_3D', PH_PREV_3D)
        call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV_ALT_CO2_3D', PH_PREV_ALT_CO2_3D)

        ! Set refBottomDepth and refBottomDepthTopOfCell
        do k = 1, nVertLevels
           refBottomDepth(k) = config_ecosys_column_bottom_depth * interfaceLocations(k+1)
           refZMid(k) = - 0.5_RKIND * config_ecosys_column_bottom_depth * (interfaceLocations(k) + interfaceLocations(k+1))
        end do

        ! Set vertCoordMovementWeights
        vertCoordMovementWeights(:) = 1.0_RKIND

        if (nVertLevelsInputColumn /= nVertLevels) return

        numTracersTotal = 32  ! T,S + 30 eco
        allocate(ecoFieldColumn(nVertLevelsInputColumn, numTracersTotal))
        allocate(indexField(numTracersTotal))

        if ( associated(activeTracers) ) then
           fieldName = 'temperature'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_TS_filename, &
                nVertLevelsInputColumn, 1, ecoFieldColumn, iErr)

           fieldName = 'salinity'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_TS_filename, &
                nVertLevelsInputColumn, 2, ecoFieldColumn, iErr)

           call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_dummy)
           indexField(1) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_dummy)
           indexField(2) = index_dummy
           do iCell = 1, nCellsSolve
              do k = 1, nVertLevels
                 activeTracers(indexField(1), k, iCell) = ecoFieldColumn(k,1)
                 activeTracers(indexField(2), k, iCell) = ecoFieldColumn(k,2)
              end do
           end do
        end if

        if ( associated(ecosysTracers) ) then

           fieldName = 'PO4'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 3, ecoFieldColumn, iErr)
           fieldName = 'NO3'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 4, ecoFieldColumn, iErr)
           fieldName = 'SiO3'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 5, ecoFieldColumn, iErr)
           fieldName = 'NH4'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 6, ecoFieldColumn, iErr)
           fieldName = 'Fe'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 7, ecoFieldColumn, iErr)
           fieldName = 'O2'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 8, ecoFieldColumn, iErr)
           fieldName = 'DIC'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 9, ecoFieldColumn, iErr)
           fieldName = 'DIC_ALT_CO2'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 10, ecoFieldColumn, iErr)
           fieldName = 'ALK'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 11, ecoFieldColumn, iErr)
           fieldName = 'DOC'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 12, ecoFieldColumn, iErr)
           fieldName = 'DON'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 13, ecoFieldColumn, iErr)
           fieldName = 'DOFe'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 14, ecoFieldColumn, iErr)
           fieldName = 'DOP'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 15, ecoFieldColumn, iErr)
           fieldName = 'DOPr'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 16, ecoFieldColumn, iErr)
           fieldName = 'DONr'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 17, ecoFieldColumn, iErr)
           fieldName = 'zooC'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 18, ecoFieldColumn, iErr)
           fieldName = 'spChl'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 19, ecoFieldColumn, iErr)
           fieldName = 'spC'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 20, ecoFieldColumn, iErr)
           fieldName = 'spFe'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 21, ecoFieldColumn, iErr)
           fieldName = 'spCaCO3'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 22, ecoFieldColumn, iErr)
           fieldName = 'diatChl'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 23, ecoFieldColumn, iErr)
           fieldName = 'diatC'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 24, ecoFieldColumn, iErr)
           fieldName = 'diatFe'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 25, ecoFieldColumn, iErr)
           fieldName = 'diatSi'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 26, ecoFieldColumn, iErr)
           fieldName = 'diazChl'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 27, ecoFieldColumn, iErr)
           fieldName = 'diazC'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 28, ecoFieldColumn, iErr)
           fieldName = 'diazFe'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 29, ecoFieldColumn, iErr)
           fieldName = 'phaeoChl'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 30, ecoFieldColumn, iErr)
           fieldName = 'phaeoC'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 31, ecoFieldColumn, iErr)
           fieldName = 'phaeoFe'
           call ocn_init_setup_ecosys_read_column(domain, fieldName, config_ecosys_column_ecosys_filename, &
                nVertLevelsInputColumn, 32, ecoFieldColumn, iErr)

           call mpas_pool_get_dimension(tracersPool, 'index_PO4',      index_dummy)
           indexField(3) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_NO3',      index_dummy)
           indexField(4) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_SiO3',     index_dummy)
           indexField(5) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_NH4',      index_dummy)
           indexField(6) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_Fe',       index_dummy)
           indexField(7) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_O2',       index_dummy)
           indexField(8) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_DIC',      index_dummy)
           indexField(9) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_DIC_ALT_CO2', index_dummy)
           indexField(10) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_ALK',      index_dummy)
           indexField(11) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_DOC',      index_dummy)
           indexField(12) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_DON',      index_dummy)
           indexField(13) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_DOFe',     index_dummy)
           indexField(14) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_DOP',      index_dummy)
           indexField(15) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_DOPr',     index_dummy)
           indexField(16) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_DONr',     index_dummy)
           indexField(17) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_zooC',     index_dummy)
           indexField(18) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_spChl',    index_dummy)
           indexField(19) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_spC',      index_dummy)
           indexField(20) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_spFe',     index_dummy)
           indexField(21) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_spCaCO3',  index_dummy)
           indexField(22) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_diatChl',  index_dummy)
           indexField(23) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_diatC',    index_dummy)
           indexField(24) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_diatFe',   index_dummy)
           indexField(25) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_diatSi',   index_dummy)
           indexField(26) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_diazChl',  index_dummy)
           indexField(27) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_diazC',    index_dummy)
           indexField(28) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_diazFe',   index_dummy)
           indexField(29) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_phaeoChl', index_dummy)
           indexField(30) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_phaeoC',   index_dummy)
           indexField(31) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_phaeoFe',  index_dummy)
           indexField(32) = index_dummy

           do iField = 3, numTracersTotal
              do iCell = 1, nCellsSolve
                 do k = 1, nVertLevels
                    ecosysTracers(indexField(iField), k, iCell) = ecoFieldColumn(k,iField)
                 end do
              end do
           end do

           do iCell = 1, nCellsSolve
              PH_PREV(iCell) = 8.0_RKIND
              PH_PREV_ALT_CO2(iCell) = 8.0_RKIND
              do k = 1, nVertLevels
                 PH_PREV_3D(k, iCell) = 8.0_RKIND
                 PH_PREV_ALT_CO2_3D(k, iCell) = 8.0_RKIND
              end do
           end do

        end if  !  associated(ecosysTracers)

        if ( associated(DMSTracers) ) then
           call mpas_pool_get_dimension(tracersPool, 'index_DMS', index_dummy)
           indexField(1) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_DMSP', index_dummy)
           indexField(2) = index_dummy
           do iCell = 1, nCellsSolve
              do k = 1, nVertLevels
                 DMSTracers(indexField(1), k, iCell) = 0.0_RKIND
                 DMSTracers(indexField(2), k, iCell) = 0.0_RKIND
              end do
           end do
        end if  !  associated(DMSTracers)

        if ( associated(MacroMoleculesTracers) ) then
           call mpas_pool_get_dimension(tracersPool, 'index_PROT', index_dummy)
           indexField(1) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_POLY', index_dummy)
           indexField(2) = index_dummy
           call mpas_pool_get_dimension(tracersPool, 'index_LIP', index_dummy)
           indexField(3) = index_dummy
           do iCell = 1, nCellsSolve
              do k = 1, nVertLevels
                 MacroMoleculesTracers(indexField(1), k, iCell) = 0.0_RKIND
                 MacroMoleculesTracers(indexField(2), k, iCell) = 0.0_RKIND
                 MacroMoleculesTracers(indexField(3), k, iCell) = 0.0_RKIND
              end do
           end do
        end if  !  associated(MacroMoleculesTracers)

        do iCell = 1, nCellsSolve
           ! Set layerThickness
           do k = 1, nVertLevels
              layerThickness(k, iCell) = config_ecosys_column_bottom_depth * (interfaceLocations(k+1) - interfaceLocations(k))
              restingThickness(k, iCell) = layerThickness(k, iCell)
           end do

           ! Set bottomDepth
           bottomDepth(iCell) = config_ecosys_column_bottom_depth

           ! Set maxLevelCell
           maxLevelCell(iCell) = nVertLevels
        end do

        block_ptr => block_ptr % next
      end do

      deallocate(interfaceLocations)
      deallocate(ecoFieldColumn, indexField)

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_ecosys_column!}}}

!***********************************************************************
!
!  routine ocn_init_validate_ecosys_column
!
!> \brief   Validation for ecosys column test case
!> \author  Mathew Maltrud
!> \date    11/01/2015
!> \details
!>  This routine validates the configuration options for the ecosys column test configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_ecosys_column(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_ecosys_column_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('ecosys_column')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_ecosys_column_vert_levels', config_ecosys_column_vert_levels)

      if(config_vert_levels <= 0 .and. config_ecosys_column_vert_levels > 0) then
         config_vert_levels = config_ecosys_column_vert_levels
      else if(config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for ecosys column test case. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_ecosys_column!}}}

!***********************************************************************

!
!  routine ocn_init_setup_ecosys_read_column
!
!> \brief   Read a column of a specified field from a given file
!> \author  Mathew Maltrud
!> \date    11/01/2014
!> \details
!>  This routine reads a column of a specified field from a given file
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_ecosys_read_column(domain, fieldName, fileName, &
              nVertLevelsInputColumn, iField, ecoFieldColumn, iErr)!{{{

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr
      integer, intent(in) :: nVertLevelsInputColumn, iField
      character (len=StrKIND), intent(in) :: fieldName, fileName
      real (kind=RKIND), dimension(:,:), intent(inout) :: ecoFieldColumn

       type (block_type), pointer :: block_ptr

       type (MPAS_Stream_type) :: columnStream

       character (len=StrKIND), pointer :: config_global_ocean_temperature_file, config_global_ocean_temperature_varname, &
                                           config_global_ocean_tracer_nlon_dimname, config_global_ocean_tracer_nlat_dimname, &
                                           config_global_ocean_depth_dimname

       integer :: k

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_nlon_dimname', &
                                 config_global_ocean_tracer_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_dimname', config_global_ocean_depth_dimname)

       ! Define stream for reading a column
!      call MPAS_createStream(columnStream, domain % iocontext, fileName, MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)
       call MPAS_createStream(columnStream, domain % iocontext, fileName, MPAS_IO_NETCDF, MPAS_IO_READ)

       ! Setup field for stream to be read in
       columnIC % fieldName = trim(fieldName)
       columnIC % dimSizes(1) = nVertLevelsInputColumn
       columnIC % dimSizes(2) = 1
       columnIC % dimNames(1) = 'nVertLevels'
       columnIC % dimNames(2) = 'nCells'
       columnIC % isVarArray = .false.
       columnIC % isPersistent = .true.
       columnIC % isActive = .true.
       columnIC % hasTimeDimension = .false.
       columnIC % block => domain % blocklist
       allocate(columnIC % array(nVertLevelsInputColumn, 1))

       ! Add column field to stream
       call MPAS_streamAddField(columnStream, columnIC, iErr)

       ! Read stream
       call MPAS_readStream(columnStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(columnStream)

       do k = 1, nVertLevelsInputColumn
          ecoFieldColumn(k,iField) = columnIC % array(k,1)
       end do

    end subroutine ocn_init_setup_ecosys_read_column

!***********************************************************************

end module ocn_init_ecosys_column

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
