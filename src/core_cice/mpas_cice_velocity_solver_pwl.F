!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_velocity_solver_pwl
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_velocity_solver_pwl

  use mpas_grid_types

  implicit none

  private
  save

  public :: &
       cice_init_velocity_solver_pwl, &
       cice_internal_stress_pwl, &
       cice_strain_tensor_pwl, &
       cice_stress_divergence_pwl

contains

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_velocity_solver_pwl
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_init_velocity_solver_pwl(&
       mesh, &
       velocity_variational, &
       velocity_pwl, &
       boundary, &
       rotateCartesianGrid, &
       includeMetricTerms)!{{{

    use cice_mesh, only: &
         cice_cell_vertices_at_vertex

    use cice_velocity_solver_variational, only: &
         cice_calc_local_coords, &
         cice_calc_variational_metric_terms
    
    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    type(MPAS_pool_type), pointer :: &
         velocity_variational, & !< Input/Output: 
         velocity_pwl, &         !< Input/Output: 
         boundary                !< Input/Output: 

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input: 
         includeMetricTerms     !< Input: 

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex

    real(kind=RKIND), dimension(:), pointer :: &
         tanLatVertexRotatedOverRadius

    real(kind=RKIND), dimension(:,:), pointer :: &
         xLocal, &
         yLocal, &
         basisSubArea

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisGradientU, &
         basisGradientV

    real(kind=RKIND), dimension(:,:,:,:), pointer :: &
         basisIntegralsMetric

    call MPAS_pool_get_array(velocity_variational, "xLocal", xLocal)
    call MPAS_pool_get_array(velocity_variational, "yLocal", yLocal)
    call MPAS_pool_get_array(velocity_variational, "cellVerticesAtVertex", cellVerticesAtVertex)
    call MPAS_pool_get_array(velocity_variational, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)

    call MPAS_pool_get_array(velocity_pwl, "basisSubArea", basisSubArea)
    call MPAS_pool_get_array(velocity_pwl, "basisGradientU", basisGradientU)
    call MPAS_pool_get_array(velocity_pwl, "basisGradientV", basisGradientV)
    call MPAS_pool_get_array(velocity_pwl, "basisIntegralsMetric", basisIntegralsMetric)

    call cice_calc_local_coords(&
         mesh, &
         xLocal, &
         yLocal, &
         rotateCartesianGrid)

    call cice_calc_variational_metric_terms(&
         mesh, &
         tanLatVertexRotatedOverRadius, &
         rotateCartesianGrid, &
         includeMetricTerms)

    call cice_cell_vertices_at_vertex(&
         mesh, &
         cellVerticesAtVertex)

    call init_velocity_solver_pwl_basis(&
         mesh, &
         basisSubArea, &
         basisGradientU, &
         basisGradientV, &
         basisIntegralsMetric, &
         xLocal, &
         yLocal)

  end subroutine cice_init_velocity_solver_pwl!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_velocity_solver_pwl_basis
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_velocity_solver_pwl_basis(&
       mesh, &
       basisSubArea, &
       basisGradientU, &
       basisGradientV, &
       basisIntegralsMetric, &
       xLocal, &
       yLocal)!{{{

    use cice_numerics, only: &
         cice_solve_linear_basis_system

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         basisSubArea !< Output: 
         
    real(kind=RKIND), dimension(:,:,:), intent(out) :: &    
         basisGradientU, & !< Output: 
         basisGradientV    !< Output: 

    real(kind=RKIND), dimension(:,:,:,:), intent(out) :: & 
         basisIntegralsMetric !< Output:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         xLocal, & !< Input: 
         yLocal    !< Input: 

    real(kind=RKIND) :: &
         xPWLCentre, &
         yPWLCentre, &
         alphaPWL, &
         a, b, c, s

    integer :: &
         iCell, &
         iVertexOnCell, &
         iEdgeOnCell, &
         iEdge, &
         iVertexOnCell1, &
         iVertexOnCell2, &
         iSubCell, &
         jEdgeOnCell

    real(kind=RKIND), dimension(3,3) :: &
         leftMatrix

    real(kind=RKIND), dimension(3) :: &
         rightHandSide, &
         solutionVector

    integer, pointer :: &
         nCells, &
         maxEdges

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell, &
         edgesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         dvEdge

    real(kind=RKIND), dimension(:,:), allocatable :: &
         subBasisGradientU, &
         subBasisGradientV

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "maxEdges", maxEdges)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "dvEdge", dvEdge)

    allocate(subBasisGradientU(maxEdges,3))
    allocate(subBasisGradientV(maxEdges,3))

    ! loop over cells
    do iCell = 1, nCells
       
       alphaPWL = 1.0_RKIND / real(nEdgesOnCell(iCell),RKIND)

       ! determine cell centre for piecewise linear basis
       xPWLCentre = 0.0_RKIND
       yPWLCentre = 0.0_RKIND

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          xPWLCentre = xPWLCentre + alphaPWL * xLocal(iVertexOnCell,iCell)
          yPWLCentre = yPWLCentre + alphaPWL * yLocal(iVertexOnCell,iCell)

       enddo ! iVertexOnCell

       ! calculate the area of the subcells
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          iEdge = edgesOnCell(iEdgeOnCell,iCell)
          iVertexOnCell1 = iEdgeOnCell
          iVertexOnCell2 = wrapped_index(iEdgeOnCell + 1, nEdgesOnCell(iCell))

          c = dvEdge(iEdge)
          a = sqrt((xLocal(iVertexOnCell1,iCell) - xPWLCentre)**2 + &
                   (yLocal(iVertexOnCell1,iCell) - yPWLCentre)**2)
          b = sqrt((xLocal(iVertexOnCell2,iCell) - xPWLCentre)**2 + &
                   (yLocal(iVertexOnCell2,iCell) - yPWLCentre)**2)

          s = (a + b + c) * 0.5_RKIND

          ! Heron's formula 
          basisSubArea(iEdgeOnCell,iCell) = sqrt(s * (s-a) * (s-b) * (s-c))

       enddo ! iEdgeOnCell

       ! calculate the linear basis on the sub triangle

       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          iVertexOnCell1 = iEdgeOnCell
          iVertexOnCell2 = wrapped_index(iEdgeOnCell + 1, nEdgesOnCell(iCell))

          ! set up left hand matrix
          leftMatrix(1,1) = xLocal(iVertexOnCell1,iCell) - xPWLCentre
          leftMatrix(1,2) = yLocal(iVertexOnCell1,iCell) - yPWLCentre
          leftMatrix(1,3) = 1.0_RKIND
          
          leftMatrix(2,1) = xLocal(iVertexOnCell2,iCell) - xPWLCentre
          leftMatrix(2,2) = yLocal(iVertexOnCell2,iCell) - yPWLCentre
          leftMatrix(2,3) = 1.0_RKIND
          
          leftMatrix(3,1) = 0.0_RKIND
          leftMatrix(3,2) = 0.0_RKIND
          leftMatrix(3,3) = 1.0_RKIND

          ! first basis
          rightHandSide(1) = 1.0_RKIND
          rightHandSide(2) = 0.0_RKIND
          rightHandSide(3) = 0.0_RKIND

          call cice_solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

          subBasisGradientU(1,iEdgeOnCell) = solutionVector(1)
          subBasisGradientV(1,iEdgeOnCell) = solutionVector(2)

          ! second basis
          rightHandSide(1) = 0.0_RKIND
          rightHandSide(2) = 1.0_RKIND
          rightHandSide(3) = 0.0_RKIND

          call cice_solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

          subBasisGradientU(2,iEdgeOnCell) = solutionVector(1)
          subBasisGradientV(2,iEdgeOnCell) = solutionVector(2)

          ! third basis
          subBasisGradientU(3,iEdgeOnCell) = 1.0_RKIND - subBasisGradientU(1,iEdgeOnCell) - subBasisGradientU(2,iEdgeOnCell)
          subBasisGradientV(3,iEdgeOnCell) = 1.0_RKIND - subBasisGradientV(1,iEdgeOnCell) - subBasisGradientV(2,iEdgeOnCell)

       enddo ! iEdgeOnCell

       ! use the linear sub area basis to calculate the PWL basis
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          ! loop over subcells
          do iSubCell = 1, nEdgesOnCell(iCell)

             ! array (index of the basis vertex, subarea value)
             basisGradientU(iEdgeOnCell,iSubCell,iCell) = subBasisGradientU(3,iSubCell) * alphaPWL
             basisGradientV(iEdgeOnCell,iSubCell,iCell) = subBasisGradientV(3,iSubCell) * alphaPWL

             if (iSubCell == iEdgeOnCell) then

                basisGradientU(iEdgeOnCell,iSubCell,iCell) = basisGradientU(iEdgeOnCell,iSubCell,iCell) + subBasisGradientU(1,iSubCell)
                basisGradientV(iEdgeOnCell,iSubCell,iCell) = basisGradientV(iEdgeOnCell,iSubCell,iCell) + subBasisGradientV(1,iSubCell)

             else if (iSubCell == wrapped_index(iEdgeOnCell - 1, nEdgesOnCell(iCell))) then

                basisGradientU(iEdgeOnCell,iSubCell,iCell) = basisGradientU(iEdgeOnCell,iSubCell,iCell) + subBasisGradientU(2,iSubCell)
                basisGradientV(iEdgeOnCell,iSubCell,iCell) = basisGradientV(iEdgeOnCell,iSubCell,iCell) + subBasisGradientV(2,iSubCell)

             endif

          enddo ! iSubCell
          
       enddo ! iEdgeOnCell

       ! basis integrals for the metric terms
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)
          do jEdgeOnCell = 1, nEdgesOnCell(iCell)
             do iSubCell = 1, nEdgesOnCell(iCell)

                basisIntegralsMetric(iEdgeOnCell,jEdgeOnCell,iSubCell,iCell) = alphaPWL**2

                if ( (iSubCell == iEdgeOnCell .and. &
                      iSubCell == jEdgeOnCell) .or. &
                     (iSubCell == wrapped_index(iEdgeOnCell - 1, nEdgesOnCell(iCell)) .and. &
                      iSubCell == wrapped_index(jEdgeOnCell - 1, nEdgesOnCell(iCell))) ) then

                   basisIntegralsMetric(iEdgeOnCell,jEdgeOnCell,iSubCell,iCell) = 2.0_RKIND * (alphaPWL + 1)**2

                endif

                if ( (iSubCell == iEdgeOnCell .and. &
                      iSubCell == wrapped_index(jEdgeOnCell - 1, nEdgesOnCell(iCell))) .or. &
                     (iSubCell == wrapped_index(iEdgeOnCell - 1, nEdgesOnCell(iCell)) .and. &
                      iSubCell == jEdgeOnCell) ) then

                   basisIntegralsMetric(iEdgeOnCell,jEdgeOnCell,iSubCell,iCell) = 2.0_RKIND * alphaPWL**2 + 2.0_RKIND * alphaPWL + 1.0_RKIND

                endif

                basisIntegralsMetric(iEdgeOnCell,jEdgeOnCell,iSubCell,iCell) = basisIntegralsMetric(iEdgeOnCell,jEdgeOnCell,iSubCell,iCell) * &
                     basisSubArea(iSubCell,iCell) / 12.0_RKIND

             enddo ! iSubCell
          enddo ! jEdgeOnCell
       enddo ! iSubCell

    enddo ! iCell

    deallocate(subBasisGradientU)
    deallocate(subBasisGradientV)

  end subroutine init_velocity_solver_pwl_basis!}}}

!-----------------------------------------------------------------------
! time step
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_internal_stress_pwl
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_internal_stress_pwl(&
       mesh, &
       velocity_variational, &
       velocity_pwl, &
       uVelocity, &
       vVelocity, &
       icePressure, &
       replacementPressure, &
       dtElastic, &
       evpDamping, &
       stressDivergenceU, &
       stressDivergenceV, &
       solveStress, &
       solveVelocity)!{{{

    use cice_velocity_solver_variational, only: &
         cice_stress_tensor_variational

    type (MPAS_pool_type), pointer :: &
         mesh, &   !< Input/Output: 
         velocity_variational, & !< Input/Output: 
         velocity_pwl    !< Input/Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, & !< Input: 
         vVelocity    !< Input:   

    real(kind=RKIND), dimension(:), intent(inout) :: &
         icePressure, &      !< Input/Output: 
         replacementPressure !< Input/Output: 

    real(kind=RKIND), intent(in) :: &
         dtElastic !< Input: 

    logical, intent(in) :: &
         evpDamping !< Input: 

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, & !< Output: 
         stressDivergenceV    !< Output: 

    integer, dimension(:), intent(in) :: &
         solveStress, & !< Input: 
         solveVelocity  !< Input: 

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex

    real(kind=RKIND), dimension(:), pointer :: &
         tanLatVertexRotatedOverRadius

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11, &
         strain22, &
         strain12, &
         stress11, &
         stress22, &
         stress12, &
         basisSubArea

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisGradientUPWL, &
         basisGradientVPWL

    real(kind=RKIND), dimension(:,:,:,:), pointer :: &    
         basisIntegralsMetric

    call MPAS_pool_get_array(velocity_variational, "strain11", strain11)
    call MPAS_pool_get_array(velocity_variational, "strain22", strain22)
    call MPAS_pool_get_array(velocity_variational, "strain12", strain12)
    call MPAS_pool_get_array(velocity_variational, "stress11", stress11)
    call MPAS_pool_get_array(velocity_variational, "stress22", stress22)
    call MPAS_pool_get_array(velocity_variational, "stress12", stress12)
    call MPAS_pool_get_array(velocity_variational, "cellVerticesAtVertex", cellVerticesAtVertex)
    call MPAS_pool_get_array(velocity_variational, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)

    call MPAS_pool_get_array(velocity_pwl, "basisGradientU", basisGradientUPWL)
    call MPAS_pool_get_array(velocity_pwl, "basisGradientV", basisGradientVPWL)
    call MPAS_pool_get_array(velocity_pwl, "basisSubArea", basisSubArea)
    call MPAS_pool_get_array(velocity_pwl, "basisIntegralsMetric", basisIntegralsMetric)

    call cice_strain_tensor_pwl(mesh, &
         strain11, &
         strain22, &
         strain12, &
         uVelocity, &
         vVelocity, &
         basisGradientUPWL, &
         basisGradientVPWL, &
         tanLatVertexRotatedOverRadius, &
         solveStress)
    
    call cice_stress_tensor_variational(mesh, &
         stress11, &
         stress22, &
         stress12, &
         strain11, &
         strain22, &
         strain12, &
         icePressure, &
         replacementPressure, &
         solveStress, &
         dtElastic, evpDamping)
    
    call cice_stress_divergence_pwl(mesh, &
         stressDivergenceU, &
         stressDivergenceV, &
         stress11, &
         stress22, &
         stress12, &
         basisGradientUPWL, &
         basisGradientVPWL, &
         basisSubArea, &
         basisIntegralsMetric, &
         tanLatVertexRotatedOverRadius, &
         cellVerticesAtVertex, &
         solveVelocity)

  end subroutine cice_internal_stress_pwl!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_strain_tensor_pwl
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------  

  subroutine cice_strain_tensor_pwl(&
       mesh, &
       strain11, &
       strain22, &
       strain12, &
       uVelocity, &
       vVelocity, &
       basisGradientU, &
       basisGradientV, &
       tanLatVertexRotatedOverRadius, &
       solveStress)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         strain11, & !< Output: 
         strain22, & !< Output: 
         strain12    !< Output: 
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, & !< Input: 
         vVelocity, & !< Input: 
         tanLatVertexRotatedOverRadius !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         basisGradientU, & !< Input: 
         basisGradientV    !< Input: 
    
    integer, dimension(:), intent(in) :: &
         solveStress !< Input: 

    real(kind=RKIND) :: &
         alphaPWL

    integer :: &
         iCell, &
         jVertexOnCell, &
         iVertexOnCell, &
         iVertex, &
         jVertex, &
         iSubCell1, &
         iSubCell2

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)

    ! loop over cells
    do iCell = 1, nCells

       strain11(:,iCell) = 0.0_RKIND
       strain22(:,iCell) = 0.0_RKIND
       strain12(:,iCell) = 0.0_RKIND

       if (solveStress(iCell) == 1) then

          alphaPWL = 1.0_RKIND / real(nEdgesOnCell(iCell), RKIND)

          ! loop over velocity points surrounding cell - location of stress and derivative
          do jVertexOnCell = 1, nEdgesOnCell(iCell)

             jVertex = verticesOnCell(jVertexOnCell,iCell)

             iSubCell1 = jVertexOnCell
             iSubCell2 = wrapped_index(jVertexOnCell + 1, nEdgesOnCell(iCell))

             ! loop over basis functions
             do iVertexOnCell = 1, nEdgesOnCell(iCell)

                iVertex = verticesOnCell(iVertexOnCell,iCell)

                strain11(jVertexOnCell,iCell) = strain11(jVertexOnCell,iCell) + &
                     uVelocity(iVertex) * 0.5_RKIND * (basisGradientU(iVertexOnCell,iSubCell1,iCell) + basisGradientU(iVertexOnCell,iSubCell2,iCell))

                strain22(jVertexOnCell,iCell) = strain22(jVertexOnCell,iCell) + &
                     vVelocity(iVertex) * 0.5_RKIND * (basisGradientV(iVertexOnCell,iSubCell1,iCell) + basisGradientV(iVertexOnCell,iSubCell2,iCell))

                strain12(jVertexOnCell,iCell) = strain12(jVertexOnCell,iCell) + 0.5_RKIND * ( &
                     uVelocity(iVertex) * 0.5_RKIND * (basisGradientV(iVertexOnCell,iSubCell1,iCell) + basisGradientV(iVertexOnCell,iSubCell2,iCell)) + &
                     vVelocity(iVertex) * 0.5_RKIND * (basisGradientU(iVertexOnCell,iSubCell1,iCell) + basisGradientU(iVertexOnCell,iSubCell2,iCell)))

             enddo ! iVertexOnCell

             ! metric terms
             strain11(jVertexOnCell,iCell) = strain11(jVertexOnCell,iCell) - &
                  vVelocity(jVertex) * tanLatVertexRotatedOverRadius(jVertex)

             strain12(jVertexOnCell,iCell) = strain12(jVertexOnCell,iCell) + &
                  uVelocity(jVertex) * tanLatVertexRotatedOverRadius(jVertex) * 0.5_RKIND

          enddo ! jVertexOnCell

       endif ! solveStress

    enddo ! iCell

  end subroutine cice_strain_tensor_pwl!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_stress_divergence_pwl
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_stress_divergence_pwl(&
       mesh, &
       stressDivergenceU, &
       stressDivergenceV, &
       stress11, &
       stress22, &
       stress12, &
       basisGradientU, &
       basisGradientV, &
       basisSubArea, &
       basisIntegralsMetric, &
       tanLatVertexRotatedOverRadius, &
       cellVerticesAtVertex, &
       solveVelocity)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, & !< Output: 
         stressDivergenceV    !< Output: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &    
         stress11, & !< Input: 
         stress22, & !< Input: 
         stress12    !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &   
         basisGradientU, & !< Input: 
         basisGradientV    !< Input: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         basisSubArea !< Input: 

    real(kind=RKIND), dimension(:,:,:,:), intent(in) :: & 
         basisIntegralsMetric !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &     
         tanLatVertexRotatedOverRadius !< Input:

    integer, dimension(:,:), intent(in) :: &
         cellVerticesAtVertex !< Input: 

    integer, dimension(:), intent(in) :: &
         solveVelocity !< Input: 

    real(kind=RKIND) :: &
         stressDivergenceUCell, &
         stressDivergenceVCell

    integer :: &
         iVertex, &
         iCellOnVertex, &
         iCell, &
         iVertexOnCell, &
         jVertexOnCell, &
         iSubCell, &
         iVertex2

    integer, pointer :: &
         nVertices, &
         vertexDegree

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         cellsOnVertex, &
         verticesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         areaTriangle

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "areaTriangle", areaTriangle)

    ! loop over velocity positions
    do iVertex = 1, nVertices

       stressDivergenceU(iVertex) = 0.0_RKIND
       stressDivergenceV(iVertex) = 0.0_RKIND

       if (solveVelocity(iVertex) == 1) then

          ! loop over surrounding cells
          do iCellOnVertex = 1, vertexDegree
             
             ! get the cell number of this cell
             iCell = cellsOnVertex(iCellOnVertex, iVertex)
             
             ! get the vertexOnCell number of the iVertex velocity point from cell iCell
             jVertexOnCell = cellVerticesAtVertex(iCellOnVertex,iVertex)

             stressDivergenceUCell = 0.0_RKIND
             stressDivergenceVCell = 0.0_RKIND

             ! loop over the vertices of the surrounding cell
             do iVertexOnCell = 1, nEdgesOnCell(iCell)

                iVertex2 = verticesOnCell(iVertexOnCell,iCell)

                ! loop over subcells
                do iSubCell = 1, nEdgesOnCell(iCell)

                   ! normal terms
                   stressDivergenceUCell = stressDivergenceUCell + &
                        stress11(iVertexOnCell,iCell) * basisSubArea(iSubCell,iCell) * basisGradientU(iVertexOnCell,iSubCell,iCell) + &
                        stress12(iVertexOnCell,iCell) * basisSubArea(iSubCell,iCell) * basisGradientV(iVertexOnCell,iSubCell,iCell)
                   
                   stressDivergenceVCell = stressDivergenceVCell + &
                        stress22(iVertexOnCell,iCell) * basisSubArea(iSubCell,iCell) * basisGradientV(iVertexOnCell,iSubCell,iCell) + &
                        stress12(iVertexOnCell,iCell) * basisSubArea(iSubCell,iCell) * basisGradientU(iVertexOnCell,iSubCell,iCell)

                   ! metric terms
                   stressDivergenceUCell = stressDivergenceUCell + &
                        stress12(iVertexOnCell,iCell) * basisIntegralsMetric(iVertexOnCell,jVertexOnCell,iSubCell,iCell) * &
                        tanLatVertexRotatedOverRadius(iVertex2)
                   
                   stressDivergenceVCell = stressDivergenceVCell - &
                        stress11(iVertexOnCell,iCell) * basisIntegralsMetric(iVertexOnCell,jVertexOnCell,iSubCell,iCell) * &
                        tanLatVertexRotatedOverRadius(iVertex2)

                enddo ! iSubCell
                
             enddo ! jVertex

             stressDivergenceU(iVertex) = stressDivergenceU(iVertex) - stressDivergenceUCell
             stressDivergenceV(iVertex) = stressDivergenceV(iVertex) - stressDivergenceVCell
             
          enddo ! iCellOnVertex

          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) / areaTriangle(iVertex)
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) / areaTriangle(iVertex)

       endif ! solveVelocity

    enddo ! iVertex
 
  end subroutine cice_stress_divergence_pwl!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  pwl_basis_gradient
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function pwl_basis_gradient(&
       nEdgesOnCell, &
       basisGradient, &
       alphaPWL, &
       basisVertexOnCell, &
       iSubCell) &
       result(grad)!{{{
    
    integer, intent(in) :: &
         nEdgesOnCell, &      !< Input: 
         basisVertexOnCell, & !< Input: basis function vertex
         iSubCell             !< Input: subcell to calculate

    real(kind=RKIND), intent(in) :: &
         alphaPWL !< Input: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         basisGradient !< Input: 

    real(kind=RKIND) :: grad
    
    grad = basisGradient(1,iSubCell) * &
           (alphaPWL + merge(1.0_RKIND, 0.0_RKIND, iSubCell == basisVertexOnCell)) + &
           basisGradient(2,iSubCell) * &
           (alphaPWL + merge(1.0_RKIND, 0.0_RKIND, iSubCell == wrapped_index(basisVertexOnCell - 1, nEdgesOnCell)))

  end function pwl_basis_gradient!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wrapped_index
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function wrapped_index(&
       input, &
       nelements) &
       result(output)!{{{

    integer, intent(in) :: &
         input, &  !< Input: 
         nelements !< Input: 

    integer :: output

    output = mod(input - 1, nelements) + 1

  end function wrapped_index!}}}

!-----------------------------------------------------------------------

end module cice_velocity_solver_pwl
