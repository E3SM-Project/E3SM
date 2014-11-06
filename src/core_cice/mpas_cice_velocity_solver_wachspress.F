!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_velocity_solver_wachspress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_velocity_solver_wachspress

  use mpas_grid_types

  implicit none

  private
  save

  public :: &
       cice_init_velocity_solver_wachspress

contains

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_velocity_solver_wachspress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_init_velocity_solver_wachspress(&
       mesh, &
       velocity_variational, &
       boundary, &
       rotateCartesianGrid, &
       includeMetricTerms)!{{{

    use cice_mesh, only: &
         cice_cell_vertices_at_vertex

    use cice_velocity_solver_variational_shared, only: &
         cice_calc_local_coords, &
         cice_calc_variational_metric_terms

    use cice_debug, only: &
         cice_plot_cells_planar

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    type(MPAS_pool_type), pointer :: &
         velocity_variational, & !< Input/Output: 
         boundary                !< Input/Output: 

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input: 
         includeMetricTerms     !< Input: 

    integer :: &
         iCell, &
         iVertex

    integer, pointer :: &
         nCells, &
         maxEdges, i1, i2

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex

    real(kind=RKIND), dimension(:), pointer :: &    
         tanLatVertexRotatedOverRadius

    real(kind=RKIND), dimension(:,:), allocatable :: &
         xLocal, &
         yLocal

    real(kind=RKIND), dimension(:,:), allocatable :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND), dimension(:,:,:), allocatable :: &
         wachspressKappa

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisGradientU, &
         basisGradientV, &
         basisIntegralsU, &
         basisIntegralsV, &
         basisIntegralsMetric

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "maxEdges", maxEdges)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)

    call MPAS_pool_get_array(velocity_variational, "cellVerticesAtVertex", cellVerticesAtVertex)
    call MPAS_pool_get_array(velocity_variational, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)
    call MPAS_pool_get_array(velocity_variational, "basisGradientU", basisGradientU)
    call MPAS_pool_get_array(velocity_variational, "basisGradientV", basisGradientV)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsU", basisIntegralsU)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsV", basisIntegralsV)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsMetric", basisIntegralsMetric)

    allocate(xLocal(maxEdges,nCells))
    allocate(yLocal(maxEdges,nCells))
    allocate(wachspressKappa(maxEdges,maxEdges,nCells))
    allocate(wachspressA(maxEdges,nCells))
    allocate(wachspressB(maxEdges,nCells))

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

    call calc_wachspress_coefficients(&
         mesh, &
         wachspressKappa, &
         wachspressA, &
         wachspressB, &
         xLocal, &
         yLocal)

    call calculate_wachspress_derivatives(&
         mesh, &
         basisGradientU, &
         basisGradientV, &
         xLocal, &
         yLocal, &
         wachspressA, & 
         wachspressB, &
         wachspressKappa)

    call integrate_wachspress(&
         mesh, &
         basisIntegralsU, &
         basisIntegralsV, &
         basisIntegralsMetric, &
         xLocal, &        
         yLocal, &
         wachspressA, &   
         wachspressB, &
         wachspressKappa)

    call cice_cell_vertices_at_vertex(&
         mesh, &
         cellVerticesAtVertex)

    iCell = 1
    iVertex = 1

    call plot_wachpress(mesh, &
         iCell, &
         nEdgesOnCell(iCell), &
         iVertex, &
         wachspressKappa(:,:,iCell), &
         wachspressA(:,iCell), &
         wachspressB(:,iCell), &
         xLocal(:,iCell), &
         yLocal(:,iCell))

    stop

    deallocate(xLocal)
    deallocate(yLocal)
    deallocate(wachspressKappa)
    deallocate(wachspressA)
    deallocate(wachspressB)

    iCell = 4000
    do i1 = 1, maxEdges
       do i2 = 1, maxEdges
          write(70,*) i1, i2, &
               basisGradientU(i1,i2,iCell), basisGradientV(i1,i2,iCell), &
               basisIntegralsMetric(i1,i2,iCell), &
               basisIntegralsU(i1,i2,iCell), basisIntegralsV(i1,i2,iCell)
       enddo
    enddo

    call cice_plot_cells_planar(mesh, basisGradientU(3,3,:), "basisGradientU_wach")
    call cice_plot_cells_planar(mesh, basisGradientV(3,3,:), "basisGradientV_wach")
    call cice_plot_cells_planar(mesh, basisIntegralsMetric(3,3,:), "basisIntegralsMetric_wach")
    call cice_plot_cells_planar(mesh, basisIntegralsU(3,3,:), "basisIntegralsU_wach")
    call cice_plot_cells_planar(mesh, basisIntegralsV(3,3,:), "basisIntegralsV_wach")
    !stop

  end subroutine cice_init_velocity_solver_wachspress!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calc_wachspress_coefficients
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine calc_wachspress_coefficients(&
       mesh, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       xLocal, &
       yLocal)!{{{
    
    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         wachspressKappa !< Output: 

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         wachspressA, & !< Output: 
         wachspressB    !< Output: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         xLocal, & !< Input: 
         yLocal    !< Input: 

    integer :: &
         iCell, &
         iVertex, &
         i0, &
         i1, &
         i2, &
         jVertex

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)

    ! loop over cells
    do iCell = 1, nCells

       ! loop over vertices
       do iVertex = 1, nEdgesOnCell(iCell)
          
          ! end points of line segment
          i1 = iVertex - 1
          i2 = iVertex
          if (i1 < 1) i1 = i1 + nEdgesOnCell(iCell)
          
          ! solve for the line segment equation
          wachspressA(iVertex, iCell) = (yLocal(i2,iCell) - yLocal(i1,iCell)) / (xLocal(i1,iCell) * yLocal(i2,iCell) - xLocal(i2,iCell) * yLocal(i1,iCell))
          wachspressB(iVertex, iCell) = (xLocal(i1,iCell) - xLocal(i2,iCell)) / (xLocal(i1,iCell) * yLocal(i2,iCell) - xLocal(i2,iCell) * yLocal(i1,iCell))

       enddo ! iVertex

       ! loop over vertices
       do iVertex = 1, nEdgesOnCell(iCell)

          ! determine kappa
          wachspressKappa(1,iVertex,iCell) = 1.0_RKIND

          do jVertex = 2, nEdgesOnCell(iCell)
             
             ! previous, this and next vertex
             i0 = jVertex - 1
             i1 = jVertex
             i2 = jVertex + 1
             if (i2 > nEdgesOnCell(iCell)) i2 = i2 - nEdgesOnCell(iCell)
             
             wachspressKappa(jVertex,iVertex,iCell) = wachspressKappa(jVertex-1,iVertex,iCell) * &
                  (wachspressA(i2,iCell) * (xLocal(i0,iCell) - xLocal(i1,iCell)) + wachspressB(i2,iCell) * (yLocal(i0,iCell) - yLocal(i1,iCell))) / &
                  (wachspressA(i0,iCell) * (xLocal(i1,iCell) - xLocal(i0,iCell)) + wachspressB(i0,iCell) * (yLocal(i1,iCell) - yLocal(i0,iCell)))

          enddo ! jVertex
          
       enddo ! iVertex
       
    enddo ! iCell

  end subroutine calc_wachspress_coefficients!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wachspress_basis_function
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function wachspress_basis_function(&
       nEdgesOnCell, &
       iVertex, &
       x, &
       y, &
       wachspressKappa, &
       wachspressA, &
       wachspressB) &
       result(wachpress)!{{{

    integer, intent(in) :: &
         nEdgesOnCell, & !< Input: 
         iVertex         !< Input: 

    real(kind=RKIND), intent(in) :: &
         x, & !< Input: 
         y    !< Input: 
    
    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input: 
         wachspressB    !< Input: 

    real(kind=RKIND) :: &
         wachpress !< Output: 

    real(kind=RKIND) :: &
         numerator, &
         denominator, &
         numerator_ivertex

    integer :: &
         jVertex

    ! sum over numerators to get denominator
    denominator = 0.0_RKIND

    do jVertex = 1, nEdgesOnCell

       numerator = wachspress_numerator(nEdgesOnCell, jVertex, iVertex, x, y, wachspressKappa, wachspressA, wachspressB)

       denominator = denominator + numerator

       if (jvertex == iVertex) then
          numerator_ivertex = numerator
       endif

    enddo ! jVertex

    wachpress = numerator_ivertex / denominator

  end function wachspress_basis_function!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wachspress_basis_derivative
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function wachspress_basis_derivative(&
       nEdgesOnCell, &
       iVertex, &
       x, &
       y, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       iDerivativeType) &
       result(wachspress)!{{{

    integer, intent(in) :: &
         nEdgesOnCell, & !< Input: 
         iVertex, &      !< Input: 
         iDerivativeType !< Input: 

    real(kind=RKIND), intent(in) :: &
         x, & !< Input: 
         y    !< Input: 
    
    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input: 
         wachspressB    !< Input: 

    real(kind=RKIND) :: &
         wachspress !< Output: 

    real(kind=RKIND) :: &
         numerator, &
         derivative, &
         denominator, &
         sum_of_derivatives, &
         numerator_ivertex, &
         derivative_ivertex

    integer :: &
         jVertex

    ! sum over numerators to get denominator
    denominator = 0.0_RKIND
    sum_of_derivatives = 0.0_RKIND

    do jVertex = 1, nEdgesOnCell

       numerator  = wachspress_numerator(nEdgesOnCell, jVertex, iVertex, x, y, wachspressKappa, wachspressA, wachspressB)
       derivative = wachspress_numerator_derivative(nEdgesOnCell, jVertex, iVertex, x, y, wachspressKappa, wachspressA, wachspressB, iDerivativeType)

       denominator        = denominator        + numerator
       sum_of_derivatives = sum_of_derivatives + derivative

       if (jvertex == iVertex) then
          numerator_ivertex  = numerator
          derivative_ivertex = derivative
       endif

    enddo ! jVertex

    wachspress = derivative_ivertex / denominator - (numerator_ivertex / denominator**2) * sum_of_derivatives

  end function wachspress_basis_derivative!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wachspress_numerator
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function wachspress_numerator(&
       nEdgesOnCell, &
       jVertex, &
       iVertex, &
       x, &
       y, &
       wachspressKappa, &
       wachspressA, &
       wachspressB) &
       result(numerator)!{{{

    integer, intent(in) :: &
         nEdgesOnCell, & !< Input: 
         jVertex, &      !< Input: 
         iVertex         !< Input: 

    real(kind=RKIND), intent(in) :: &
         x, & !< Input: 
         y    !< Input: 
    
    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input: 
         wachspressB    !< Input: 

    real(kind=RKIND) :: &
         numerator !< Output: 

    integer :: &
         kVertex, &
         i1, &
         i2
    
    i1 = jVertex 
    i2 = jVertex + 1
    if (i2 > nEdgesOnCell) i2 = i2 - nEdgesOnCell

    numerator = 1.0_RKIND

    do kVertex = 1, nEdgesOnCell

       if (kVertex /= i1 .and. kVertex /= i2) then
          numerator = numerator * wachspress_edge_equation(kVertex, x, y, wachspressA, wachspressB)
       endif

    enddo ! jVertex

    numerator = numerator * wachspressKappa(jVertex,iVertex)

  end function wachspress_numerator!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wachspress_numerator_derivative
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function wachspress_numerator_derivative(&
       nEdgesOnCell, &
       jVertex, &
       iVertex, &
       x, &
       y, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       iDerivativeType) &
       result(derivative)!{{{

    integer, intent(in) :: &
         nEdgesOnCell, & !< Input: 
         jVertex, &      !< Input: 
         iVertex, &      !< Input: 
         iDerivativeType !< Input: 

    real(kind=RKIND), intent(in) :: &
         x, & !< Input: 
         y    !< Input: 
    
    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input: 
         wachspressB    !< Input: 

    real(kind=RKIND) :: &
         derivative !< Output: 

    real(kind=RKIND) :: &
         sum_of_products, &
         product

    integer :: &
         kVertex, &
         lVertex, &
         i1, &
         i2

    i1 = jVertex 
    i2 = jVertex + 1
    if (i2 > nEdgesOnCell) i2 = i2 - nEdgesOnCell

    sum_of_products = 0.0_RKIND

    do kVertex = 1, nEdgesOnCell

       if (kVertex /= i1 .and. kVertex /= i2) then
          
          product = 1.0_RKIND

          do lVertex = 1, nEdgesOnCell
             
             if (lVertex /= i1 .and. lVertex /= i2) then
                
                if (lVertex == kVertex) then

                   product = product * wachspress_edge_equation_derivative(lVertex, wachspressA, wachspressB, iDerivativeType)
                
                else

                   product = product * wachspress_edge_equation(lVertex, x, y, wachspressA, wachspressB)

                endif

             endif
             
          enddo ! jVertex

          sum_of_products = sum_of_products + product
          
       endif

    enddo ! jVertex

    derivative = sum_of_products * wachspressKappa(jVertex,iVertex)

  end function wachspress_numerator_derivative!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wachspress_edge_equation
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function wachspress_edge_equation(&
       iVertex, &
       x, &
       y, &
       wachspressA, &
       wachspressB) &
       result(edge_equation)!{{{

    integer, intent(in) :: &
         iVertex !< Input: 

    real(kind=RKIND), intent(in) :: &
         x, & !< Input: 
         y    !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input: 
         wachspressB    !< Input: 

    real(kind=RKIND) :: &
         edge_equation !< Output: 

    edge_equation = 1.0_RKIND - wachspressA(iVertex) * x - wachspressB(iVertex) * y

  end function wachspress_edge_equation!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  wachspress_edge_equation_derivative
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function wachspress_edge_equation_derivative(&
       iVertex, &
       wachspressA, &
       wachspressB, &
       iDerivativeType) &
       result(derivative)!{{{

    integer, intent(in) :: &
         iVertex, &      !< Input: 
         iDerivativeType !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input: 
         wachspressB    !< Input: 

    real(kind=RKIND) :: &
         derivative !< Output: 

    if (iDerivativeType == 1) then

       derivative = -wachspressA(iVertex)

    else if (iDerivativeType == 2) then

       derivative = -wachspressB(iVertex)

    endif

  end function wachspress_edge_equation_derivative!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calculate_wachspress_derivatives
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine calculate_wachspress_derivatives(&
       mesh, &
       basisGradientU, &
       basisGradientV, &
       xLocal, &
       yLocal, &
       wachspressA, &
       wachspressB, &
       wachspressKappa)!{{{
    
    ! basisGradientUV(jVertexOnCell,iVertexOnCell,iCell)
    ! iCell         : The cell the gradients are based in
    ! iVertexOnCell : The vertex basis function the gradient is calculated from
    ! jVertexOnCell : The vertex location the gradients are calculated at
    
    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         basisGradientU, & !< Output: 
         basisGradientV    !< Output: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &    
         wachspressA, & !< Input: 
         wachspressB, & !< Input: 
         xLocal,      & !< Input: 
         yLocal         !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &  
         wachspressKappa !< Input: 

    integer :: &
         iCell, &
         iBasisVertex, &
         iGradientVertex

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

       ! loop over vertices - basis function
       do iBasisVertex = 1, nEdgesOnCell(iCell)

          ! loop over vertices again - derivative position
          do iGradientVertex = 1, nEdgesOnCell(iCell)
          
             basisGradientU(iBasisVertex,iGradientVertex,iCell) = &
                  wachspress_basis_derivative(nEdgesOnCell(iCell), &
                                              iBasisVertex, &
                                              xLocal(iGradientVertex,iCell), &
                                              yLocal(iGradientVertex,iCell), &
                                              wachspressKappa(:,:,iCell), &
                                              wachspressA(:,iCell), &
                                              wachspressB(:,iCell), &
                                              1)

             basisGradientV(iBasisVertex,iGradientVertex,iCell) = &
                  wachspress_basis_derivative(nEdgesOnCell(iCell), &
                                              iBasisVertex, &
                                              xLocal(iGradientVertex,iCell), &
                                              yLocal(iGradientVertex,iCell), &
                                              wachspressKappa(:,:,iCell), &
                                              wachspressA(:,iCell), &
                                              wachspressB(:,iCell), &
                                              2)

           enddo ! iGradientVertex

       enddo ! iBasisVertex

    enddo ! iCell

  end subroutine calculate_wachspress_derivatives!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_triangle_mapping
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine get_triangle_mapping(&
       mapping, &
       x1, y1, &
       x2, y2, &
       u1, v1, &
       u2, v2)!{{{

    real(kind=RKIND), dimension(2,2), intent(out) :: &
         mapping !< Output: 

    real(kind=RKIND), intent(in) :: &
         x1, & !< Input: 
         y1, & !< Input:
         x2, & !< Input:
         y2, & !< Input:
         u1, & !< Input:
         v1, & !< Input:
         u2, & !< Input:
         v2    !< Input:

    mapping(1,1) = (u2*y1 - u1*y2) / (x2*y1 - x1*y2)
    mapping(1,2) = (u1*x2 - u2*x1) / (y1*x2 - y2*x1)

    mapping(2,1) = (v2*y1 - v1*y2) / (x2*y1 - x1*y2)
    mapping(2,2) = (v1*x2 - v2*x1) / (y1*x2 - y2*x1)

  end subroutine get_triangle_mapping!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  use_triangle_mapping
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine use_triangle_mapping(&
       u, v, &
       x, y, &
       mapping)!{{{

    real(kind=RKIND), intent(out) :: &
         u, v !< Output: 

    real(kind=RKIND), intent(in) :: &
         x, y !< Input:

    real(kind=RKIND), dimension(2,2), intent(in) :: &
         mapping !< Input:

    u = mapping(1,1) * x + mapping(1,2) * y
    v = mapping(2,1) * x + mapping(2,2) * y

  end subroutine use_triangle_mapping!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  integrate_wachspress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine integrate_wachspress(&
       mesh, &
       basisIntegralsU, &
       basisIntegralsV, &
       basisIntegralsMetric, &
       xLocal, &
       yLocal, &
       wachspressA, &
       wachspressB, &
       wachspressKappa)!{{{

    ! basisIntegralsUV (iStressVertex,iVelocityVertex,iCell)
    ! iCell         : cell integrals are performed on
    ! iStressVertex : vertex number of Wachspress function 
    ! iVelocityVertex : vertex number of Wachspress derivative function 
    ! Sij

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         basisIntegralsU, &   !< Output: 
         basisIntegralsV, &   !< Output: 
         basisIntegralsMetric !< Output: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         xLocal, &      !< Input:
         yLocal, &      !< Input:
         wachspressA, & !< Input:
         wachspressB    !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &    
         wachspressKappa !< Input:

    real(kind=RKIND) :: &
         integration

    integer :: &
         iCell, &
         iStressVertex, &
         iVelocityVertex, &
         iIntegrationType

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)

    do iCell = 1, nCells

       if (iCell == 1) then

       do iStressVertex = 1, nEdgesOnCell(iCell)

          do iVelocityVertex = 1, nEdgesOnCell(iCell)

             do iIntegrationType = 1, 3

                call integrate_wachspress_polygon(&
                     integration, &
                     nEdgesOnCell(iCell), &
                     iStressVertex, &
                     iVelocityVertex, &
                     xLocal(:,iCell), &
                     yLocal(:,iCell), &
                     wachspressA(:,iCell), &
                     wachspressB(:,iCell), &
                     wachspressKappa(:,:,iCell), &
                     iIntegrationType)

                if (iIntegrationType == 1) then

                   basisIntegralsU(iStressVertex,iVelocityVertex,iCell) = integration

                else if (iIntegrationType == 2) then

                   basisIntegralsV(iStressVertex,iVelocityVertex,iCell) = integration

                else if (iIntegrationType == 3) then

                   basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) = integration

                endif

             enddo ! iIntegrationType

          enddo ! jVertex

       enddo ! iVertex

       else

          basisIntegralsU(:,:,iCell)      = basisIntegralsU(:,:,1)
          basisIntegralsV(:,:,iCell)      = basisIntegralsV(:,:,1)
          basisIntegralsMetric(:,:,iCell) = basisIntegralsMetric(:,:,1)

       endif

    enddo ! iCell

  end subroutine integrate_wachspress!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  integrate_wachspress_polygon
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine integrate_wachspress_polygon(&
       integration, &
       nEdgesOnCell, &
       iStressVertex, &
       iVelocityVertex, &
       xLocal, &
       yLocal, &
       wachspressA, &
       wachspressB, &
       wachspressKappa, &
       iIntegrationType)!{{{

    integer, intent(in) :: &
         nEdgesOnCell, &  !< Input:
         iStressVertex, & !< Input:
         iVelocityVertex, & !< Input:
         iIntegrationType  !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xLocal, & !< Input:
         yLocal    !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB    !< Input:

    real(kind=RKIND), intent(out) :: &
         integration !< Output: 

    real(kind=RKIND) :: &
         integration_subtriangle

    real(kind=RKIND), dimension(2,2) :: &
         mapping

    integer :: &
         iSubTriangle, &
         i1, &
         i2

    integration = 0.0_RKIND

    do iSubTriangle = 1, nEdgesOnCell

       i1 = iSubTriangle
       i2 = iSubTriangle + 1
       if (i2 > nEdgesOnCell) i2 = i2 - nEdgesOnCell

       call get_triangle_mapping(&
            mapping, &
            1.0_RKIND, 0.0_RKIND, &
            0.0_RKIND, 1.0_RKIND, &
            xLocal(i1), yLocal(i1), &
            xLocal(i2), yLocal(i2))

       call integrate_wachspress_subtriangle(&
            integration_subtriangle, &
            nEdgesOnCell, &
            iStressVertex, &
            iVelocityVertex, &
            wachspressKappa(:,:), &
            wachspressA(:), &
            wachspressB(:), &
            iIntegrationType, &
            mapping)

       integration = integration + integration_subtriangle

    enddo ! iSubTriangle

  end subroutine integrate_wachspress_polygon!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  integrate_wachspress_subtriangle
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine integrate_wachspress_subtriangle(&
       integration, &
       nEdgesOnCell, &
       iStressVertex, &
       iVelocityVertex, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       iIntegrationType, &
       mapping)!{{{

    real(kind=RKIND), intent(out) :: &
         integration !< Output: 

    integer, intent(in) :: &
         nEdgesOnCell, &  !< Input:
         iStressVertex, & !< Input:
         iVelocityVertex, & !< Input:
         iIntegrationType  !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB    !< Input:

    real(kind=RKIND), dimension(2,2), intent(in) :: &
         mapping !< Input:

    real(kind=RKIND) :: &
         jacobian

    real(kind=RKIND) :: &
         scaling, &
         x, &
         y, &
         u, &
         v

    integer :: &
         i, j, k

    integer, parameter :: n = 10

    integration = 0.0_RKIND

    do i = 0, n
       do j = 0, n

          if (i<=n-j) then

             if (i==n .or. j==n .or. (i==0 .and. j==0)) then

                scaling = 1.0_RKIND

             else if ((j==0 .and. i/=0 .and. i/=n) .or. (i==0 .and. j/=0 .and. j/=n) .or. (i==n-j .and. i/=0 .and. j/=0)) then

                scaling = 3.0_RKIND

             else

                scaling = 6.0_RKIND 

             endif

             u = real(i,RKIND) / real(n,RKIND)
             v = real(j,RKIND) / real(n,RKIND)
             
             call use_triangle_mapping(&
                  x, y, &
                  u, v, &
                  mapping)
             
             jacobian = mapping(1,1) * mapping(2,2) - mapping(1,2) * mapping(2,1)
             
             ! area test
             !integration = integration + scaling * jacobian
             
             ! actual integration
             if (iIntegrationType == 1 .or. iIntegrationType == 2) then

                integration = integration + scaling * jacobian * &
                     wachspress_basis_function  (nEdgesOnCell, iStressVertex,   x, y, wachspressKappa, wachspressA, wachspressB) * &
                     wachspress_basis_derivative(nEdgesOnCell, iVelocityVertex, x, y, wachspressKappa, wachspressA, wachspressB, iIntegrationType)
                
             else if (iIntegrationType == 3) then

                integration = integration + scaling * jacobian * &
                     wachspress_basis_function(nEdgesOnCell, iStressVertex,   x, y, wachspressKappa, wachspressA, wachspressB) * &
                     wachspress_basis_function(nEdgesOnCell, iVelocityVertex, x, y, wachspressKappa, wachspressA, wachspressB)

             endif

          endif

       enddo ! i
    enddo ! j

    integration = integration / (6.0_RKIND * real(n,RKIND)**2)

  end subroutine integrate_wachspress_subtriangle!}}}

!-----------------------------------------------------------------------
! Hex specific plotting stuff
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_wachpress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine plot_wachpress(&
       mesh, &
       iCell, &
       nEdgesOnCell, &
       iVertex, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       xLocal, &
       yLocal)!{{{
    
    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    integer, intent(in) :: &
         iCell, &        !< Input: 
         nEdgesOnCell, & !< Input: 
         iVertex         !< Input: 

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         wachspressKappa !< Output: 

    real(kind=RKIND), dimension(:), intent(out) :: &
         wachspressA, & !< Output: 
         wachspressB    !< Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         xLocal, & !< Input: 
         yLocal    !< Input: 

    real(kind=RKIND) :: &
         xmin, xmax, ymin, ymax
    
    real(kind=RKIND), dimension(2,2) :: &
         mapping

    integer :: &
         iSubTriangle, &
         i1, &
         i2, &
         iVertexOnCell, &
         iVertex2, &
         iObject

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         xCell, &
         yCell

    ! init variables
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)

    xmin =  1e30
    xmax = -1e30
    ymin =  1e30
    ymax = -1e30

    iObject = 1

    do iVertexOnCell = 1, nEdgesOnCell

       iVertex2 = verticesOnCell(iVertexOnCell,iCell)

       xmin = min(xmin, xVertex(iVertex2)) - 1000.0_RKIND
       xmax = max(xmax, xVertex(iVertex2)) + 1000.0_RKIND
       ymin = min(ymin, yVertex(iVertex2)) - 1000.0_RKIND
       ymax = max(ymax, yVertex(iVertex2)) + 1000.0_RKIND
       
       open(55,file='test.txt')
       if (iVertexOnCell == iVertex) then
          write(55,*) xVertex(iVertex2), yVertex(iVertex2)
       endif
       close(55)

    enddo ! iVertex

    open(55,file="wachspress.txt")

    write(55,fmt='(a,f10.2,a,f10.2,a)') "set xrange [",xmin,":",xmax,"]"
    write(55,fmt='(a,f10.2,a,f10.2,a)') "set yrange [",ymin,":",ymax,"]"

    ! loop over subtriangles
    do iSubTriangle = 1, nEdgesOnCell

       i1 = iSubTriangle
       i2 = iSubTriangle + 1
       if (i2 > nEdgesOnCell) i2 = i2 - nEdgesOnCell

       ! get the triangle mapping
       call get_triangle_mapping(mapping, &
                                 1.0_RKIND,0.0_RKIND,&
                                 0.0_RKIND,1.0_RKIND,&
                                 xLocal(i1),yLocal(i1),&
                                 xLocal(i2),yLocal(i2))    
       
       ! plot from subtriangle
       call plot_subtriangle_2(nEdgesOnCell, iSubTriangle, iVertex, &
                             wachspressKappa, wachspressA, wachspressB, &
                             mapping, xCell(iCell), yCell(iCell), iObject)

    enddo ! iSubTriangle

    close(55)

  end subroutine plot_wachpress!}}}


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_subtriangle
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine plot_subtriangle(&
       nEdgesOnCell, &
       iSubTriangle, &
       iVertex, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       mapping, &
       x0, &
       y0, &
       iObject)!{{{

    integer, intent(in) :: &
         nEdgesOnCell, & !< Input: 
         iSubTriangle, & !< Input: 
         iVertex         !< Input: 

    integer, intent(inout) :: &
         iObject         !< Input: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input: 
         wachspressB    !< Input: 

    real(kind=RKIND), dimension(2,2), intent(in) :: &
         mapping !< Input: 

    real(kind=RKIND), intent(in) :: &
         x0, & !< Input: 
         y0    !< Input: 

    real(kind=RKIND) :: &
         u, v, &
         x, y, &
         wachspress, &
         wachspress1

    integer :: &
         i, j

    real(kind=RKIND) :: &
         x1, x2, x3, x4, x5, &
         y1, y2, y3, y4, y5, &
         d, &
         dx, dy

    integer, parameter :: n = 5

    logical, parameter :: &
         lfinitedifference = .false.

    d = 1.0_RKIND / real(n, RKIND)

    do i = 0, n
       do j = 0, n

          if (i<=n-j-1) then

             iObject = iObject + 1  

             u = real(i,RKIND) / real(n,RKIND)
             v = real(j,RKIND) / real(n,RKIND)

             call use_triangle_mapping(x, y, u, v, mapping)

             if (lfinitedifference) then
                wachspress = &
                     wachspress_basis_function(nEdgesOnCell, iVertex, x, y, wachspressKappa, wachspressA, wachspressB)
                
                dx = 0.0_RKIND
                dy = 100.0_RKIND
                
                x1 = x + dx
                y1 = y + dy
                
                wachspress1 = &
                     wachspress_basis_function(nEdgesOnCell, iVertex, x1, y1, wachspressKappa, wachspressA, wachspressB)             
                
                wachspress = (wachspress1 - wachspress) / dy
                
             else

                wachspress = &
                     wachspress_basis_function(nEdgesOnCell, iVertex, x, y, wachspressKappa, wachspressA, wachspressB)     
                     !wachspress_basis_derivative(nEdgesOnCell, iVertex, x, y, wachspressKappa, wachspressA, wachspressB, 2)

             endif

             if (i<n-j-1) then
                
                call use_triangle_mapping(x1, y1, u,   v,   mapping)
                call use_triangle_mapping(x2, y2, u,   v+d, mapping)
                call use_triangle_mapping(x3, y3, u+d, v+d, mapping)
                call use_triangle_mapping(x4, y4, u+d, v,   mapping)
                x5 = x1
                y5 = y1
                
                x1 = x1 + x0 ; y1 = y1 + y0
                x2 = x2 + x0 ; y2 = y2 + y0
                x3 = x3 + x0 ; y3 = y3 + y0
                x4 = x4 + x0 ; y4 = y4 + y0
                x5 = x5 + x0 ; y5 = y5 + y0
                
                write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') "set object ",iObject," polygon from ",&
                     x1,",",y1," to ",x2,",",y2," to ",x3,",",y3," to ",x4,",",y4," to ",x5,",",y5     

                write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', wachspress, ' fillstyle solid'

             else

                call use_triangle_mapping(x1, y1, u,   v,   mapping)
                call use_triangle_mapping(x2, y2, u,   v+d, mapping)
                call use_triangle_mapping(x4, y4, u+d, v,   mapping)
                x5 = x1
                y5 = y1

                x1 = x1 + x0 ; y1 = y1 + y0
                x2 = x2 + x0 ; y2 = y2 + y0
                x4 = x4 + x0 ; y4 = y4 + y0
                x5 = x5 + x0 ; y5 = y5 + y0
                
                write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') "set object ",iObject," polygon from ",&
                     x1,",",y1," to ",x2,",",y2," to ",x4,",",y4," to ",x5,",",y5    

                write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', wachspress, ' fillstyle solid'  

             endif

          endif

       enddo ! j
    enddo ! i

  end subroutine plot_subtriangle!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_subtriangle_2
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 31 October 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine plot_subtriangle_2(&
       nEdgesOnCell, &
       iSubTriangle, &
       iVertex, &
       wachspressKappa, &
       wachspressA, &
       wachspressB, &
       mapping, &
       x0, &
       y0, &
       iObject)!{{{

    integer, intent(in) :: &
         nEdgesOnCell, & !< Input: 
         iSubTriangle, & !< Input: 
         iVertex         !< Input:
 
    integer, intent(inout) :: &
         iObject         !< Input: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input: 
         wachspressB    !< Input: 

    real(kind=RKIND), dimension(2,2), intent(in) :: &
         mapping !< Input: 

    real(kind=RKIND), intent(in) :: &
         x0, & !< Input: 
         y0    !< Input: 

    integer :: &
         i, j, &
         ii, jj

    real(kind=RKIND) :: &
         u, v, &
         uu0, vv0, &
         uu1, vv1, &
         uu2, vv2, &
         xx0, yy0, &
         xx1, yy1, &
         xx2, yy2, &
         f0, f1, f2, &
         f_avg, &
         d, &
         u0, u1, u2, &
         v0, v1, v2, &
         uuu, vvv, &
         x, y

    integer, parameter :: n = 50
    integer, parameter :: nn = 1

    logical, parameter :: &
         lfinitedifference = .false.

    d = (1.0_RKIND / real(n,RKIND))

    ! loop over integration segments
    do i = 0, n
       do j = 0, n

          ! forward integration triangle
          if (i+j<=n-1) then

             ! subtriangle starting point
             u0 = real(i,RKIND) / real(n,RKIND)
             v0 = real(j,RKIND) / real(n,RKIND)
             call use_triangle_mapping(x, y, u0, v0, mapping)
             f0 = wachspress_basis_function(nEdgesOnCell, iVertex, x, y, wachspressKappa, wachspressA, wachspressB)

             u1 = real(i+1,RKIND) / real(n,RKIND)
             v1 = real(j  ,RKIND) / real(n,RKIND)
             call use_triangle_mapping(x, y, u1, v1, mapping)
             f1 = wachspress_basis_function(nEdgesOnCell, iVertex, x, y, wachspressKappa, wachspressA, wachspressB)

             u2 = real(i  ,RKIND) / real(n,RKIND)
             v2 = real(j+1,RKIND) / real(n,RKIND)
             call use_triangle_mapping(x, y, u2, v2, mapping)
             f2 = wachspress_basis_function(nEdgesOnCell, iVertex, x, y, wachspressKappa, wachspressA, wachspressB)

             u = u0
             v = v0

             ! break up the forward triangle for sub-plotting
             do ii = 0, nn
                do jj = 0, nn

                   ! forward triangle
                   if (ii+jj<=nn-1) then

                      uu0 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv0 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu1 = real(ii+1,RKIND) / real(nn,RKIND)
                      vv1 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu2 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv2 = real(jj+1,RKIND) / real(nn,RKIND)

                      ! subtriangle location in canonical coordinates
                      uuu = (uu0 + uu1 + uu2) / 3.0_RKIND
                      vvv = (vv0 + vv1 + vv2) / 3.0_RKIND

                      ! linear basis
                      f_avg = f0 * (-uuu-vvv+1) + &
                              f1 * uuu + &
                              f2 * vvv

                      uu0 = u + uu0 * d
                      vv0 = v + vv0 * d
                      call use_triangle_mapping(xx0, yy0, uu0, vv0, mapping)

                      uu1 = u + uu1 * d
                      vv1 = v + vv1 * d
                      call use_triangle_mapping(xx1, yy1, uu1, vv1, mapping)

                      uu2 = u + uu2 * d
                      vv2 = v + vv2 * d
                      call use_triangle_mapping(xx2, yy2, uu2, vv2, mapping)

                      xx0 = x0 + xx0 ; yy0 = y0 + yy0
                      xx1 = x0 + xx1 ; yy1 = y0 + yy1
                      xx2 = x0 + xx2 ; yy2 = y0 + yy2

                      write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') "set object ",iObject," polygon from ",&
                           xx0,",",yy0," to ",xx1,",",yy1," to ",xx2,",",yy2," to ",xx0,",",yy0   

                      write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', f_avg, ' fillstyle solid'

                      iObject = iObject + 1  

                   endif

                   ! backward triangle
                   if (ii+jj<=nn .and. ii/=0 .and. jj/=0) then

                      uu0 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv0 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu1 = real(ii-1,RKIND) / real(nn,RKIND)
                      vv1 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu2 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv2 = real(jj-1,RKIND) / real(nn,RKIND)

                      ! subtriangle location in canonical coordinates
                      uuu = (uu0 + uu1 + uu2) / 3.0_RKIND
                      vvv = (vv0 + vv1 + vv2) / 3.0_RKIND

                      ! linear basis
                      f_avg = f0 * (-uuu-vvv+1) + &
                              f1 * uuu + &
                              f2 * vvv

                      uu0 = u + uu0 * d
                      vv0 = v + vv0 * d
                      call use_triangle_mapping(xx0, yy0, uu0, vv0, mapping)

                      uu1 = u + uu1 * d
                      vv1 = v + vv1 * d
                      call use_triangle_mapping(xx1, yy1, uu1, vv1, mapping)

                      uu2 = u + uu2 * d
                      vv2 = v + vv2 * d
                      call use_triangle_mapping(xx2, yy2, uu2, vv2, mapping)

                      xx0 = x0 + xx0 ; yy0 = y0 + yy0
                      xx1 = x0 + xx1 ; yy1 = y0 + yy1
                      xx2 = x0 + xx2 ; yy2 = y0 + yy2

                      write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') "set object ",iObject," polygon from ",&
                           xx0,",",yy0," to ",xx1,",",yy1," to ",xx2,",",yy2," to ",xx0,",",yy0   

                      write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', f_avg, ' fillstyle solid'

                      iObject = iObject + 1  

                   endif

                enddo ! jj
             enddo ! ii

          endif

          ! backward integration triangle
          if (i+j<=n .and. i/=0 .and. j/=0) then

             ! subtriangle starting point
             u0 = real(i,RKIND) / real(n,RKIND)
             v0 = real(j,RKIND) / real(n,RKIND)
             call use_triangle_mapping(x, y, u0, v0, mapping)
             f0 = wachspress_basis_function(nEdgesOnCell, iVertex, x, y, wachspressKappa, wachspressA, wachspressB)

             u1 = real(i-1,RKIND) / real(n,RKIND)
             v1 = real(j  ,RKIND) / real(n,RKIND)
             call use_triangle_mapping(x, y, u1, v1, mapping)
             f1 = wachspress_basis_function(nEdgesOnCell, iVertex, x, y, wachspressKappa, wachspressA, wachspressB)

             u2 = real(i  ,RKIND) / real(n,RKIND)
             v2 = real(j-1,RKIND) / real(n,RKIND)
             call use_triangle_mapping(x, y, u2, v2, mapping)
             f2 = wachspress_basis_function(nEdgesOnCell, iVertex, x, y, wachspressKappa, wachspressA, wachspressB)

             u = u0 - d
             v = v0 - d

             ! break up the backward triangle for sub-plotting
             do ii = 0, nn
                do jj = 0, nn

                   ! forward triangle
                   if (ii+jj>=nn .and. ii/=nn .and. jj/=nn) then

                      uu0 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv0 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu1 = real(ii+1,RKIND) / real(nn,RKIND)
                      vv1 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu2 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv2 = real(jj+1,RKIND) / real(nn,RKIND)

                      ! subtriangle location in canonical coordinates
                      uuu = (uu0 + uu1 + uu2) / 3.0_RKIND
                      vvv = (vv0 + vv1 + vv2) / 3.0_RKIND

                      uuu = 1.0_RKIND - uuu
                      vvv = 1.0_RKIND - vvv

                      ! linear basis
                      f_avg = f0 * (-uuu-vvv+1) + &
                              f1 * uuu + &
                              f2 * vvv

                      uu0 = u + uu0 * d
                      vv0 = v + vv0 * d
                      call use_triangle_mapping(xx0, yy0, uu0, vv0, mapping)

                      uu1 = u + uu1 * d
                      vv1 = v + vv1 * d
                      call use_triangle_mapping(xx1, yy1, uu1, vv1, mapping)

                      uu2 = u + uu2 * d
                      vv2 = v + vv2 * d
                      call use_triangle_mapping(xx2, yy2, uu2, vv2, mapping)

                      xx0 = x0 + xx0 ; yy0 = y0 + yy0
                      xx1 = x0 + xx1 ; yy1 = y0 + yy1
                      xx2 = x0 + xx2 ; yy2 = y0 + yy2

                      write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') "set object ",iObject," polygon from ",&
                           xx0,",",yy0," to ",xx1,",",yy1," to ",xx2,",",yy2," to ",xx0,",",yy0   

                      write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', f_avg, ' fillstyle solid'

                      iObject = iObject + 1  

                   endif

                   ! backward triangle
                   if (ii+jj>=nn+1) then

                      uu0 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv0 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu1 = real(ii-1,RKIND) / real(nn,RKIND)
                      vv1 = real(jj  ,RKIND) / real(nn,RKIND)

                      uu2 = real(ii  ,RKIND) / real(nn,RKIND)
                      vv2 = real(jj-1,RKIND) / real(nn,RKIND)

                      ! subtriangle location in canonical coordinates
                      uuu = (uu0 + uu1 + uu2) / 3.0_RKIND
                      vvv = (vv0 + vv1 + vv2) / 3.0_RKIND

                      uuu = 1.0_RKIND - uuu
                      vvv = 1.0_RKIND - vvv

                      ! linear basis
                      f_avg = f0 * (-uuu-vvv+1) + &
                              f1 * uuu + &
                              f2 * vvv

                      uu0 = u + uu0 * d
                      vv0 = v + vv0 * d
                      call use_triangle_mapping(xx0, yy0, uu0, vv0, mapping)

                      uu1 = u + uu1 * d
                      vv1 = v + vv1 * d
                      call use_triangle_mapping(xx1, yy1, uu1, vv1, mapping)

                      uu2 = u + uu2 * d
                      vv2 = v + vv2 * d
                      call use_triangle_mapping(xx2, yy2, uu2, vv2, mapping)

                      xx0 = x0 + xx0 ; yy0 = y0 + yy0
                      xx1 = x0 + xx1 ; yy1 = y0 + yy1
                      xx2 = x0 + xx2 ; yy2 = y0 + yy2

                      write(55,fmt='(a,i5,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') "set object ",iObject," polygon from ",&
                           xx0,",",yy0," to ",xx1,",",yy1," to ",xx2,",",yy2," to ",xx0,",",yy0   

                      write(55,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', f_avg, ' fillstyle solid'

                      iObject = iObject + 1  

                   endif

                enddo ! jj
             enddo ! ii

          endif

       enddo ! j
    enddo ! i

  end subroutine plot_subtriangle_2!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_wachspress_square
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine plot_wachspress_square(&
       mesh, &
       wachspressKappa, &
       wachspressA, &
       wachspressB)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         wachspressKappa !< Input: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressA, & !< Input: 
         wachspressB    !< Input: 

    integer :: &
         iCell, &
         iVertex, &
         iVertexOnCell

    real(kind=RKIND) :: &
         x, y, x0, y0, dx, dy, wachspress

    integer :: i, j

    integer, parameter :: ni = 100
    integer, parameter :: nj = 100

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell
    
    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex

    ! init variables
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nedgesOnCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)

    iCell = 1

    x0 = 1e30
    y0 = 1e30
    
    do iVertexOnCell = 1, nEdgesOnCell(iCell)
       
       iVertex = verticesOnCell(iVertexOnCell,iCell)

       x0 = min(x0,xVertex(iVertex))
       y0 = min(y0,yVertex(iVertex))

    enddo ! 
    
    dx = 16000.0_RKIND / real(ni,RKIND)
    dy = 16000.0_RKIND / real(nj,RKIND)
       
    iVertexOnCell = 1


    do i = 1, ni
       do j = 1, nj

          x = x0 + (real(i,RKIND) - 0.5_RKIND) * dx - 24000.0_RKIND
          y = x0 + (real(j,RKIND) - 0.5_RKIND) * dy - 24000.0_RKIND

          wachspress =  &
               wachspress_basis_derivative(&
               nEdgesOnCell(iCell), &
               iVertexOnCell, &
               x, &
               y, &
               wachspressKappa(:,:,iCell), &
               wachspressA(:,iCell), &
               wachspressB(:,iCell),1)

          write(55,*) x+0.5_RKIND*dx,x-0.5_RKIND*dx,y+0.5_RKIND*dy,y-0.5_RKIND*dy, wachspress          

       enddo ! j
    enddo ! i

  end subroutine plot_wachspress_square!}}}

!-----------------------------------------------------------------------

end module cice_velocity_solver_wachspress
