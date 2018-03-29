! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_setup
!
!> \brief MPAS land ice setup module
!> \author Matt Hoffman
!> \date   17 April 2011
!> \details
!>  This module contains various subroutines for
!>  setting up the land ice core.
!
!-----------------------------------------------------------------------
module li_setup

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_kind_types
   use mpas_dmpar
   use mpas_log

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   type (mpas_pool_type), public, pointer :: liConfigs  !< Public parameter: pool of config options
   type (mpas_pool_type), public, pointer :: liPackages


   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------
   public :: li_setup_config_options, &
             li_setup_vertical_grid, &
             li_setup_sign_and_index_fields, &
             li_setup_wachspress_vertex_to_cell_weights, &
             li_interpolate_vertex_to_cell_2d, &
             li_cells_to_vertices_1dfield_using_kiteAreas, &
             li_calculate_layerThickness, &
             li_compute_gradient_2d


   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------



!***********************************************************************

contains


!***********************************************************************
!
!  routine li_setup_config_options
!
!> \brief   Makes any setup changes needed based on chosen config options
!> \author  Matt Hoffman
!> \date    16 April 2014
!> \details
!>  This routine makes any adjustments as needed based on which
!>  config options were chosen.
!
!-----------------------------------------------------------------------

   subroutine li_setup_config_options( domain, err )

      use mpas_timekeeping

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err            !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      character(len=StrKIND), pointer :: config_thermal_solver
      character(len=StrKIND), pointer :: config_tracer_advection

      err = 0

      ! Make config pool publicly available in this module
      liConfigs => domain % configs
      liPackages => domain % packages

      ! ---
      ! Config-specific setup occurs here
      ! ---

      ! If thermal evolution is enabled, force tracer advection to be on
      call mpas_pool_get_config(liConfigs, 'config_thermal_solver', config_thermal_solver)
      call mpas_pool_get_config(liConfigs, 'config_tracer_advection', config_tracer_advection)
      if (( (trim(config_thermal_solver) == 'temperature') .or. &
            (trim(config_thermal_solver) == 'enthalpy') ) .and. &
            (trim(config_tracer_advection) /= 'fo') )then
         call mpas_log_write("Setting config_tracer_advection='fo' because a thermal solver has been selected. " // &
            "(config_tracer_advection was set to: " // trim(config_tracer_advection) // ")", MPAS_LOG_WARN)
         config_tracer_advection = 'fo'
      endif


      ! ---
      ! Print this run's configs to the log file
      ! ---
      call mpas_log_write("")
      call mpas_log_write("MPASLI is using the following configuration:")
      call mpas_log_write("============================================")
      call mpas_pool_print_summary(liConfigs, MPAS_POOL_CONFIG)
      call mpas_log_write("============================================")
      call mpas_log_write("")

   !--------------------------------------------------------------------
   end subroutine li_setup_config_options



!***********************************************************************
!
!  routine li_setup_vertical_grid
!
!> \brief   Initializes vertical coord system
!> \author  Matt Hoffman
!> \date    20 April 2012
!> \details
!>  This routine initializes the vertical coord system.
!
!-----------------------------------------------------------------------

   subroutine li_setup_vertical_grid(meshPool, geometryPool, err)

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
      type (mpas_pool_type), intent(inout) :: meshPool  !< Input/Output: mesh object
      type (mpas_pool_type), intent(inout) :: geometryPool  !< Input/Output: geometry object

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err            !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      ! Pool pointers
      integer, pointer :: nVertLevels ! Dimensions
      real (kind=RKIND), dimension(:), pointer :: layerThicknessFractions, layerCenterSigma, layerInterfaceSigma
      real (kind=RKIND), dimension(:), pointer :: thickness
      logical, pointer :: config_do_restart
      ! Truly locals
      integer :: k
      real (kind=RKIND) :: fractionTotal

      ! Get pool stuff
      call mpas_pool_get_config(liConfigs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      ! layerThicknessFractions is provided by input
      call mpas_pool_get_array(meshPool, 'layerThicknessFractions', layerThicknessFractions)
      call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)
      call mpas_pool_get_array(meshPool, 'layerInterfaceSigma', layerInterfaceSigma)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)

      ! Check that layerThicknessFractions are valid
      ! TODO - switch to having the user input the sigma levels instead???
      if (.not. config_do_restart) then   ! This would be applied an additional time on each restart,
                                          ! messing up the vertical levels (very, very slightly)
      fractionTotal = sum(layerThicknessFractions)
      if (fractionTotal /= 1.0_RKIND) then
         if (abs(fractionTotal - 1.0_RKIND) > 0.001_RKIND) then
            call mpas_log_write('The sum of layerThicknessFractions is different from 1.0 by more than 0.001.', MPAS_LOG_ERR)
            err = 1
         end if
         call mpas_log_write('Adjusting upper layerThicknessFrac by small amount because sum of ' // &
            'layerThicknessFractions is slightly different from 1.0.')
         layerThicknessFractions(1) = layerThicknessFractions(1) - (fractionTotal - 1.0_RKIND)
      endif
      endif

      ! layerCenterSigma is the fractional vertical position (0-1) of each layer center,
      ! with 0.0 at the ice surface and 1.0 at the ice bed
      ! layerInterfaceSigma is the fractional vertical position (0-1) of each layer interface,
      ! with 0.0 at the ice surface and 1.0 at the ice bed.
      ! Interface 1 is the surface, interface 2 is between layers 1 and 2, etc., and interface nVertLevels+1 is the bed.
      layerCenterSigma(1) = 0.5_RKIND * layerThicknessFractions(1)
      layerInterfaceSigma(1) = 0.0_RKIND
      do k = 2, nVertLevels
         layerCenterSigma(k) = layerCenterSigma(k-1) + 0.5_RKIND * layerThicknessFractions(k-1) &
            + 0.5_RKIND * layerThicknessFractions(k)
         layerInterfaceSigma(k) = layerInterfaceSigma(k-1) + layerThicknessFractions(k-1)
      end do
      layerInterfaceSigma(nVertLevels+1) = 1.0_RKIND

   !--------------------------------------------------------------------
   end subroutine li_setup_vertical_grid



!***********************************************************************
!
!  routine li_setup_sign_and_index_fields
!
!> \brief   Determines signs for various mesh items
!> \author  Matt Hoffman - based on code by Doug Jacobsen
!> \date    20 April 2012
!> \details
!>  This routine determines the sign for various mesh items.
!
!-----------------------------------------------------------------------
   subroutine li_setup_sign_and_index_fields(meshPool)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool  !< Input: mesh object

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

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      ! Pool pointers
      integer, pointer :: nCells !, nVertices, vertexDegree
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell, cellsOnEdge !, edgesOnVertex, cellsOnVertex, verticesOnCell, verticesOnEdge
      integer, dimension(:,:), pointer :: edgeSignOnCell !, edgeSignOnVertex, kiteIndexOnCell
      ! Truly locals
       integer :: iCell, iEdge, iVertex, i, j, k

      ! Get pool stuff
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

       edgeSignOnCell = 0.0_RKIND
       !edgeSignOnVertex = 0.0_RKIND
       !kiteIndexOnCell = 0.0_RKIND
       ! If needed, edgeSignOnVertex and kiteIndexOnCell can also be setup here.

       do iCell = 1, nCells
         do i = 1, nEdgesOnCell(iCell)
           iEdge = edgesOnCell(i, iCell)
           !iVertex = verticesOnCell(i, iCell)

           ! Vector points from cell 1 to cell 2
           if(iCell == cellsOnEdge(1, iEdge)) then
             edgeSignOnCell(i, iCell) = -1
           else
             edgeSignOnCell(i, iCell) =  1
           end if

           !do j = 1, vertexDegree
           !  if(cellsOnVertex(j, iVertex) == iCell) then
           !    kiteIndexOnCell(i, iCell) = j
           !  end if
           !end do
         end do
       end do

       !do iVertex = 1, nVertices
       !  do i = 1, vertexDegree
       !    iEdge = edgesOnVertex(i, iVertex)
       !
       !    ! Vector points from vertex 1 to vertex 2
       !    if(iVertex == verticesOnEdge(1, iEdge)) then
       !      edgeSignOnVertex(i, iVertex) = -1
       !    else
       !      edgeSignOnVertex(i, iVertex) =  1
       !    end if
       !  end do
       !end do

   !--------------------------------------------------------------------
   end subroutine li_setup_sign_and_index_fields


!***********************************************************************
!
!  routine li_setup_wachspress_vertex_to_cell_weights
!
!> \brief   Calculates weights to interpolate from vertices to cell centers
!> \author  Matt Hoffman
!> \date    29 Aug 2016
!> \details
!>  This routine determines the sign for various mesh items.
!
!-----------------------------------------------------------------------
   subroutine li_setup_wachspress_vertex_to_cell_weights(meshPool)
      use mpas_geometry_utils

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout), target :: meshPool  !< Input: mesh object

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! Pool pointers
      real(kind=RKIND), dimension(:,:), pointer :: wachspressWeightVertex
      real(kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real(kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell
      integer, pointer :: nCells, maxEdges
      integer :: iCell, iVertex, v
      integer :: nVerticesOnThisCell
      real(kind=RKIND), dimension(:,:), allocatable :: vertexCoordsOnCell
      type (mpas_pool_type), pointer :: meshPoolPointer

      ! Get pool stuff
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)
      call mpas_pool_get_array(meshPool, 'wachspressWeightVertex', wachspressWeightVertex)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)

      meshPoolPointer => meshPool ! mpas_wachspress_coordinates expected a pointer to the mesh pool instead of just the mesh pool itself

      allocate(vertexCoordsOnCell(3, maxEdges))

      do iCell = 1, nCells
         nVerticesOnThisCell = nEdgesOnCell(iCell)
         do v = 1, nVerticesOnThisCell
            iVertex = verticesOnCell(v, iCell)
            vertexCoordsOnCell(:, v) = (/ xVertex(iVertex), yVertex(iVertex), zVertex(iVertex) /)
         enddo
         wachspressWeightVertex(1:nVerticesOnThisCell, iCell) = mpas_wachspress_coordinates( &
               nVerticesOnThisCell, vertexCoordsOnCell(:, 1:nVerticesOnThisCell), &
               (/ xCell(iCell), yCell(iCell), zCell(iCell) /), meshPoolPointer)
      end do

      deallocate(vertexCoordsOnCell)
   !--------------------------------------------------------------------
   end subroutine li_setup_wachspress_vertex_to_cell_weights


!***********************************************************************
!
!  routine li_interpolate_vertex_to_cell_2d
!
!> \brief   Interpolates from vertices to cell centers using Wachspress functions
!> \author  Matt Hoffman
!> \date    29 Aug 2016
!> \details
!>  This routine interpolated from vertices to cell center values using
!>  Wachspress functions in MPAS operators.
!
!-----------------------------------------------------------------------
   subroutine li_interpolate_vertex_to_cell_2d(meshPool, vertexValue, cellValue)
      use mpas_geometry_utils

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool  !< Input: mesh object

      real(kind=RKIND), dimension(:,:), pointer, intent(in) :: vertexValue  !< value on vertices

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(:,:), pointer, intent(in) :: cellValue  !< value on cells

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! Pool pointers
      real(kind=RKIND), dimension(:,:), pointer :: wachspressWeightVertex
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell
      integer, pointer :: nCells, maxEdges
      integer :: iCell, iVertex, v
      integer :: nVerticesOnThisCell

      ! Get pool stuff
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)
      call mpas_pool_get_array(meshPool, 'wachspressWeightVertex', wachspressWeightVertex)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)

      cellValue(:,:) = 0.0_RKIND
      do iCell = 1, nCells
         nVerticesOnThisCell = nEdgesOnCell(iCell)
         do v = 1, nVerticesOnThisCell
            iVertex = verticesOnCell(v, iCell)
            cellValue(:, iCell) = cellValue(:, iCell) + wachspressWeightVertex(v, iCell) * vertexValue(:, iVertex)
         enddo
      end do

   !--------------------------------------------------------------------
   end subroutine li_interpolate_vertex_to_cell_2d


!***********************************************************************
!
!  subroutine li_cells_to_vertices_1dfield_using_kiteAreas
!
!> \brief   Converts a 1d scalar field from cells to vertices
!> \author  Matt Hoffman
!> \date    21 May 2012
!> \details
!>  This routine converts a 1d scalar field from cells to vertices.
!>  It will give garbage values on obtuse triangles!  But it does work
!>  on periodic meshes.
!>  TODO: It would be more efficient to calculate the weights once on init and then only
!>  perform the interp. in this routine.
!-----------------------------------------------------------------------
   subroutine li_cells_to_vertices_1dfield_using_kiteAreas(meshPool, fieldCells, fieldVertices)
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information
      real (kind=RKIND), dimension(:), intent(in) :: &
         fieldCells    !< Input: field on cells

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:), intent(out) :: &
         fieldVertices    !< Input: field on vertices

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
      integer, dimension(:,:), pointer :: cellsOnVertex
      integer, pointer :: nVertices, vertexDegree
      integer :: iCell, icell2, iVertex, cellIndex
      real (kind=RKIND) :: fVertexAccum, baryweight, weightAccum

      ! Get needed items from mesh pool
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)

      ! Calculate h on vertices using barycentric interpolation
      do iVertex = 1, nVertices  ! Loop over vertices
        fVertexAccum = 0.0_RKIND
        weightAccum = 0.0_RKIND
        ! Loop over cells on this vertex
        do iCell = 1, vertexDegree
          cellIndex = cellsOnVertex(iCell, iVertex)
          baryweight = 0.0_RKIND
          do iCell2 = 1, vertexDegree
            if (iCell2 /= icell) baryweight = baryweight + 0.5 * kiteAreasOnVertex(iCell2, iVertex)
          enddo
          fVertexAccum = fVertexAccum + baryweight * fieldCells(cellIndex)  ! add the contribution from this cell's kite
          weightAccum = weightAccum + kiteAreasOnVertex(iCell, iVertex)  ! This doesn't match areaTriangle for obtuse triangles!!!
        enddo
        fieldVertices(iVertex) = fVertexAccum / weightAccum  ! I assume this should never be 0...
      enddo

   end subroutine li_cells_to_vertices_1dfield_using_kiteAreas


!***********************************************************************
!
!  subroutine li_calculate_layerThickness
!
!> \brief   Calculates the thickness of each layer, given the total thickness
!> \author  William Lipscomb
!> \date    16 Feb 2016
!> \details
!> This routine calculates layerThickness in each cell and column,
!> given thickness in each cell. The calculation is very simple, but is
!> packaged in a subroutine to ensure that it is done the same way
!> in different parts of the code.
!-----------------------------------------------------------------------

   subroutine li_calculate_layerThickness(meshPool, thickness, layerThickness)

     !-----------------------------------------------------------------
     !
     ! input variables
     !
     !-----------------------------------------------------------------

     type (mpas_pool_type), intent(in) :: &
          meshPool                                     !< Input: mesh object

     real(kind=RKIND), dimension(:), intent(in) :: &
          thickness                                    !< Input: ice thickness

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

     real(kind=RKIND), dimension(:,:), intent(out) :: &
          layerThickness                               !< Output: thickness of each layer

     !-----------------------------------------------------------------
     !
     ! local variables
     !
     !-----------------------------------------------------------------

     integer, pointer :: &
          nCells,                    & ! number of cells
          nVertLevels                  ! number of vertical layers

     real (kind=RKIND), dimension(:), pointer :: &
          layerThicknessFractions      ! fractional thickness in each layer

     integer :: iCell, k

     call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
     call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

     call mpas_pool_get_array(meshPool, 'layerThicknessFractions', layerThicknessFractions)

     do iCell = 1, nCells
        do k = 1, nVertLevels
           layerThickness(k,iCell) = thickness(iCell) * layerThicknessFractions(k)
        enddo
     enddo

   end subroutine li_calculate_layerThickness


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine li_compute_gradient_2d
!
!> \brief  compute the gradient of a 2D scalar field
!> \author Matt Hoffman, William Lipscomb
!> \date   April 2015, Feb 2018
!> \details
!>  This routine computes the x and y coordinates of the gradient of a
!>  2D scalar field. The gradient is located at cell centers.
!>  Modified from a version in sea ice incremental remapping.
!
!-----------------------------------------------------------------------

  subroutine li_compute_gradient_2d(&
       meshPool,                     &
       field,                        &
       xGrad,  yGrad,                &
       err)

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------
    type (mpas_pool_type), intent(in) :: &
         meshPool                                     !< Input: mesh object

    real(kind=RKIND), dimension(:), intent(in) :: &
         field               !< Input: scalar field defined at cell centers

    !-----------------------------------------------------------------
    ! input/output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------
    real(kind=RKIND), dimension(:), intent(out) :: &
         xGrad, & !< Output: gradient in x-direction
         yGrad    !< Output: gradient in y-direction

    integer, intent(out) :: err            !< Output: error flag

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------
    integer, pointer :: &
         nCells,           & !< number of cells
         maxEdges            !< max number of edges per cell

    integer, dimension(:), pointer :: &
         nEdgesOnCell        !< number of edges per cell

    integer, dimension(:,:), pointer ::  &
         edgesOnCell,      & !< index for each edge of a given cell
         cellsOnCell,      & !< index for each cell neighbor of a given cell
         cellsOnEdge         !< index for each cell neighbor of a given edge

    real(kind=RKIND), dimension(:), pointer :: &
         dcEdge              !< distance between the 2 cell centers on each side of an edge

    real(kind=RKIND), dimension(:,:,:), pointer ::  &
         coeffs_reconstruct   !< coefficients for reconstructing the gradient at a cell center,
                             !<  given normal components on edges

    integer :: iCell, iEdge, iEdgeOnCell, iCellNeighbor

    real(kind=RKIND) :: &
         signGradient      ! = 1 or -1, depending on which direction is taken as positive at a given edge

    real(kind=RKIND), dimension(:), allocatable ::   &
         normalGrad  ! normal components of the gradient, defined on cell edges.  For one cell.

    err = 0

    ! get needed variables
    call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
    call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)
    call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
    call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
    call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
    call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
    call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
    call mpas_pool_get_array(meshPool, 'coeffs_reconstruct', coeffs_reconstruct)

    ! find dimensions and allocate arrays
    allocate(normalGrad(maxEdges))

    ! initialize the gradient
    xGrad(:) = 0.0_RKIND
    yGrad(:) = 0.0_RKIND

    ! loop over cells
    do iCell = 1, nCells

       ! initialize normal gradient components on edges of this cell
       normalGrad(:) = 0.0_RKIND

       ! loop over edges of this cell
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          iCellNeighbor = cellsOnCell(iEdgeOnCell, iCell)

          ! compute the normal component of the gradient on this edge
          if (iCellNeighbor >= 1 .and. iCellNeighbor <= nCells) then ! there is a cell neighbor on this edge
             iEdge = edgesOnCell(iEdgeOnCell,iCell)
             if (iCell == cellsOnEdge(1,iEdge)) then
                signGradient =  1.0_RKIND
             else
                signGradient = -1.0_RKIND
             end if
             normalGrad(iEdgeOnCell) = signGradient * (field(iCellNeighbor) - field(iCell)) / dcEdge(iEdge)
          else   ! there is no cell neighbor on this edge
             ! set gradient component = 0
             normalGrad(iEdgeOnCell) = 0.0_RKIND
          endif

          ! add the contribution of this normal component to the reconstructed
          ! gradient at the cell center (in global x/y/z coordinates)
          xGrad(iCell) = xGrad(iCell) + coeffs_reconstruct(1,iEdgeOnCell,iCell) * normalGrad(iEdgeOnCell)
          yGrad(iCell) = yGrad(iCell) + coeffs_reconstruct(2,iEdgeOnCell,iCell) * normalGrad(iEdgeOnCell)
       enddo  ! iEdgeOnCell

    enddo  ! iCell

    ! cleanup
    deallocate(normalGrad)

  end subroutine li_compute_gradient_2d


!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************



end module li_setup
