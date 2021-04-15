! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_init
!
!> \brief MPAS-ocean initialization routines.
!> \author Mark Petersen
!> \date   December 2013
!> \details
!>  This module contains routines to initialize variables at the
!>    beginning of an MPAS-Ocean simulation, or when starting the
!>    ocean analysis core.
!
!-----------------------------------------------------------------------

module ocn_init_routines

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_dmpar
   use mpas_constants

   use mpas_rbf_interpolation
   use mpas_vector_operations
   use mpas_vector_reconstruction
   use mpas_tracer_advection_helpers

   use ocn_diagnostics
   use ocn_diagnostics_variables
   use ocn_gm
   use ocn_constants
   use ocn_config
   use ocn_mesh

   use ocn_surface_land_ice_fluxes
   use ocn_forcing

   interface ocn_init_add_tau_metadata
      module procedure ocn_init_add_tau_metadata_real
      module procedure ocn_init_add_tau_metadata_int
      module procedure ocn_init_add_tau_metadata_logical
      module procedure ocn_init_add_tau_metadata_character
   end interface

   private

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

   public :: &
      ocn_init_routines_compute_min_max_level, &
      ocn_init_routines_compute_mesh_scaling, &
      ocn_init_routines_setup_sign_and_index_fields, &
      ocn_init_routines_vert_coord, &
      ocn_init_routines_block, &
      ocn_init_metadata

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_routines_compute_mesh_scaling
!
!> \brief   set up mesh scaling variables
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine initializes meshScalingDel2 and meshScalingDel4
!
!-----------------------------------------------------------------------
   subroutine ocn_init_routines_compute_mesh_scaling(meshPool, & !{{{
                 config_hmix_scaleWithMesh, config_maxMeshDensity, &
                 config_hmix_use_ref_cell_width, config_hmix_ref_cell_width)

      type (mpas_pool_type), intent(inout) :: meshPool
      logical, intent(in) :: config_hmix_scaleWithMesh
      real (kind=RKIND), intent(in) :: config_maxMeshDensity
      logical, intent(in) :: config_hmix_use_ref_cell_width
      real (kind=RKIND), intent(in) :: config_hmix_ref_cell_width

      integer :: iEdge, cell1, cell2
      integer, pointer :: nEdges
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND) :: cellWidth
      real (kind=RKIND), dimension(:), pointer :: meshDensity, meshScalingDel2, meshScalingDel4, meshScaling
      real (kind=RKIND), dimension(:), pointer :: areaCell

      call mpas_pool_get_array(meshPool, 'meshDensity', meshDensity)
      call mpas_pool_get_array(meshPool, 'meshScalingDel2', meshScalingDel2)
      call mpas_pool_get_array(meshPool, 'meshScalingDel4', meshScalingDel4)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      if (config_hmix_scaleWithMesh) then
         if (config_hmix_use_ref_cell_width) then
            ! Mesh scaling is set by areaCell and config_hmix_ref_cell_width
            ! See description of config_hmix_ref_cell_width in Registry.xml for
            ! more detail.
            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               ! Effective cell width at edge iEdge, assuming neighboring cells
               ! are circles for this calculation.
               cellWidth = 2.0_RKIND * sqrt((areaCell(cell1) + areaCell(cell2) ) / 2.0_RKIND / pii)
               meshScalingDel2(iEdge) = cellWidth / config_hmix_ref_cell_width
               meshScalingDel4(iEdge) = ( cellWidth / config_hmix_ref_cell_width )**3
            end do
         else
            ! Mesh scaling is set by meshDensity. This is both confusing and
            ! inconvenient, as the flags like config_mom_del2 need to be reset
            ! for every resolution. It is kept for backwards compatibility, but
            ! should become defunct.
            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               meshScalingDel2(iEdge) = 1.0_RKIND / ( ((meshDensity(cell1) + meshDensity(cell2) ) / 2.0_RKIND) &
                                      / config_maxMeshDensity)**(1.0_RKIND / 4.0_RKIND)  ! goes as dc**1
               meshScalingDel4(iEdge) = 1.0_RKIND / ( ((meshDensity(cell1) + meshDensity(cell2) ) / 2.0_RKIND) &
                                      / config_maxMeshDensity)**(3.0_RKIND / 4.0_RKIND)  ! goes as dc**3
            end do
         end if
      else
         ! If config_hmix_scaleWithMesh is false, hmix coefficients do not vary with
         ! cell size but remain constant across the domain.
         meshScalingDel2(:) = 1.0_RKIND
         meshScalingDel4(:) = 1.0_RKIND
      end if

   end subroutine ocn_init_routines_compute_mesh_scaling!}}}

!***********************************************************************
!
!  routine ocn_init_routines_block
!
!> \brief   Initialize blocks within MPAS-Ocean core
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine calls all block-level initializations required to begin a
!>  simulation with MPAS-Ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_init_routines_block(block, dt, err)!{{{

      type (block_type), intent(inout) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool, tracersPool
      type (mpas_pool_type), pointer :: forcingPool, scratchPool
      integer :: i, iEdge, iCell, k
      integer :: err1

      integer, dimension(:,:), pointer :: highOrderAdvectionMaskTmp, boundaryCellTmp
      integer, dimension(:,:), pointer :: edgeSignOnCellTmp, edgeSignOnVertexTmp
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity

      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroup

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels

      type (mpas_pool_iterator_type) :: groupItr

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
      call mpas_pool_get_dimension(block % dimensions, 'nVertices', nVertices)
      call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(meshPool, 'highOrderAdvectionMask', highOrderAdvectionMaskTmp)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCellTmp)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

      call mpas_initialize_deriv_two(meshPool, derivTwo, err)
      call mpas_tracer_advection_coefficients(meshPool, &
          config_horiz_tracer_adv_order, derivTwo, advCoefs, &
          advCoefs3rd, nAdvCellsForEdge, advCellsForEdge, &
          err1, maxLevelCell, highOrderAdvectionMaskTmp, &
          boundaryCellTmp)
      err = ior(err, err1)

      ! Update ocn_mesh variables for highOrderAdvectionMask
      highOrderAdvectionMask = real(highOrderAdvectionMaskTmp, kind=RKIND)

      if (.not. config_do_restart) then
         do iCell=1,nCells
            boundaryLayerDepth(iCell) = layerThickness(minLevelCell(iCell), iCell) * 0.5_RKIND
            indMLD(iCell) = minLevelCell(iCell)
         end do
      end if

      call ocn_diagnostic_solve(dt,  statePool, forcingPool, meshPool, scratchPool, tracersPool)

      ! initialize velocities and active tracers on land to be zero.
      areaCell(nCells+1) = -1.0e34_RKIND

      layerThickness(:, nCells+1) = 0.0_RKIND


      do iEdge=1, nEdges
         normalVelocity(maxLevelEdgeTop(iEdge)+1:maxLevelEdgeBot(iEdge), iEdge) = 0.0_RKIND

         normalVelocity(maxLevelEdgeBot(iEdge)+1:nVertLevels,iEdge) = -1.0e34_RKIND
         
         normalVelocity(1:minLevelEdgeTop(iEdge)-1,iEdge) = -1.0e34_RKIND
      end do

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroup, 1)
            if ( associated(tracersGroup) ) then
               do iCell=1,nCells
                  tracersGroup(:, maxLevelCell(iCell)+1:nVertLevels,iCell) =  -1.0e34_RKIND
                  tracersGroup(:, 1:minLevelCell(iCell)-1,iCell) =  -1.0e34_RKIND
               end do
            end if
         end if
      end do

      ! ------------------------------------------------------------------
      ! Accumulating various parametrizations of the transport velocity
      ! ------------------------------------------------------------------
      do iEdge = 1, nEdges
         normalTransportVelocity(:, iEdge) = normalVelocity(:, iEdge)
      end do

      ! ------------------------------------------------------------------
      ! End: Accumulating various parametrizations of the transport velocity
      ! ------------------------------------------------------------------

      call ocn_init_routines_compute_mesh_scaling(meshPool, config_hmix_scaleWithMesh, config_maxMeshDensity, &
              config_hmix_use_ref_cell_width, config_hmix_ref_cell_width)

      call mpas_rbf_interp_initialize(meshPool)
      call mpas_initialize_tangent_vectors(meshPool, edgeTangentVectors)

      call mpas_init_reconstruct(meshPool, includeHalos=.true.)

      call mpas_reconstruct(meshPool, normalVelocity,        &
                       velocityX,            &
                       velocityY,            &
                       velocityZ,            &
                       velocityZonal,        &
                       velocityMeridional    &
                      )

      call mpas_pool_initialize_time_levels(statePool)

      ! compute land-ice fluxes for potential output at startup
      call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, forcingPool, err1, 1)
      err = ior(err, err1)
      call ocn_surface_land_ice_fluxes_build_arrays(meshPool, &
                                                    forcingPool, scratchPool, statePool, dt, err1)
      err = ior(err, err1)



   end subroutine ocn_init_routines_block!}}}


!***********************************************************************
!
!  routine ocn_init_metadata
!
!> \brief   Initialize any metadata for this processor
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine sets up any metadata for this MPI task and it's associated threads.
!>  The meta data could be related to performance data, or information about
!>  all the blocks on this processor.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_metadata(domain)!{{{
      type (domain_type), intent(inout) :: domain

      type (block_type), pointer :: block

      character (len=StrKIND) :: metaDataName

      integer :: iHalo
      integer :: numBlocks
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray, nVerticesArray

      numBlocks = 0

      block => domain % blocklist
      do while ( associated(block) )
         numBlocks = numBlocks + 1

         call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
         call mpas_pool_get_dimension(block % dimensions, 'nVerticesArray', nVerticesArray)

         call ocn_init_add_tau_metadata( 'numCellHalos', size(nCellsArray) )
         call ocn_init_add_tau_metadata( 'numEdgeHalos', size(nEdgesArray) )
         call ocn_init_add_tau_metadata( 'numVertexHalos', size(nVerticesArray) )

         do iHalo = 1, size(nCellsArray)
            write(metaDataName, '(a8, i2)' ) 'cellHalo', iHalo
            call ocn_init_add_tau_metadata( metaDataName, nCellsArray(iHalo) )
         end do

         do iHalo = 1, size(nEdgesArray)
            write(metaDataName, '(a8, i2)' ) 'edgeHalo', iHalo
            call ocn_init_add_tau_metadata( metaDataName, nEdgesArray(iHalo) )
         end do

         do iHalo = 1, size(nVerticesArray)
            write(metaDataName, '(a10, i2)' ) 'vertexHalo', iHalo
            call ocn_init_add_tau_metadata( metaDataName, nVerticesArray(iHalo) )
         end do

         block => block % next
      end do

      call ocn_init_add_tau_metadata( 'numBlocks', numBlocks )

   end subroutine ocn_init_metadata!}}}


!***********************************************************************
!
!  routine ocn_init_add_tau_metadata_real
!
!> \brief   Add tau real metadata
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine adds a real value metadata for TAU to this task.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_add_tau_metadata_real( dataName, dataValue )!{{{
      character (len=*), intent(in) :: dataName
      real (kind=RKIND), intent(in) :: dataValue

      character (len=StrKIND) :: dataString

#ifdef MPAS_TAU
      dataString = ''
      write( dataString, * ) dataValue
      call tau_metadata( trim(dataName), trim(dataString) )
#endif

   end subroutine ocn_init_add_tau_metadata_real!}}}

!***********************************************************************
!
!  routine ocn_init_add_tau_metadata_int
!
!> \brief   Add tau integer metadata
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine adds an integer value metadata for TAU to this task.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_add_tau_metadata_int( dataName, dataValue )!{{{
      character (len=*), intent(in) :: dataName
      integer, intent(in) :: dataValue

      character (len=StrKIND) :: dataString

#ifdef MPAS_TAU
      dataString = ''
      write( dataString, * ) dataValue
      call tau_metadata( trim(dataName), trim(dataString) )
#endif

   end subroutine ocn_init_add_tau_metadata_int!}}}

!***********************************************************************
!
!  routine ocn_init_add_tau_metadata_logical
!
!> \brief   Add tau logical metadata
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine adds a logical value metadata for TAU to this task.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_add_tau_metadata_logical( dataName, dataValue )!{{{
      character (len=*), intent(in) :: dataName
      logical, intent(in) :: dataValue

      character (len=StrKIND) :: dataString

#ifdef MPAS_TAU
      dataString = ''
      write( dataString, * ) dataValue
      call tau_metadata( trim(dataName), trim(dataString) )
#endif

   end subroutine ocn_init_add_tau_metadata_logical!}}}

!***********************************************************************
!
!  routine ocn_init_add_tau_metadata_character
!
!> \brief   Add tau character metadata
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine adds a character value metadata for TAU to this task.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_add_tau_metadata_character( dataName, dataValue )!{{{
      character (len=*), intent(in) :: dataName
      character (len=*), intent(in) :: dataValue

#ifdef MPAS_TAU
      call tau_metadata( trim(dataName), trim(dataValue) )
#endif

   end subroutine ocn_init_add_tau_metadata_character!}}}



end module ocn_init_routines

! vim: foldmethod=marker
