module cice_dyn_hex_variational

  use mpas_grid_types

  implicit none

  private
  public :: init_dynamics_hex_variational, &
            strain_tensor_hex_variational, &
            stress_tensor_hex_variational, &
            stress_divergence_hex_variational

contains

  !-------------------------------------------------------------
  ! Initialization
  !-------------------------------------------------------------

  subroutine init_dynamics_hex_variational(mesh, hexvar, boundary)

    use cice_dynamics_shared, only: interior_vertices, &
                                    calc_cell_vertices_at_vertex

    type(mesh_type), intent(in) :: mesh

    type(hexvar_type),   pointer :: hexvar
    type(boundary_type), pointer :: boundary

    integer :: &
         iCell, &
         iVertex

    call interior_vertices(mesh, &
         boundary % interiorVertex % array)

    call calc_local_coords(mesh, &
         hexvar % xLocal % array, &
         hexvar % yLocal % array)

    call calc_wachspress_coefficients(mesh, &
         hexvar % wachspressKappa % array, &
         hexvar % wachspressA % array, &
         hexvar % wachspressB % array, &
         hexvar % xLocal % array, &
         hexvar % yLocal % array)

    call calculate_wachspress_derivatives(mesh, &
         hexvar % basisGradientU % array,  &
         hexvar % basisGradientV % array,  &
         hexvar % xLocal % array,          &
         hexvar % yLocal % array,          &
         hexvar % wachspressA % array,     & 
         hexvar % wachspressB % array,     &
         hexvar % wachspressKappa % array)

    call integrate_wachspress(mesh, &
         hexvar % basisIntegralsU % array, &
         hexvar % basisIntegralsV % array, &
         hexvar % xLocal % array,          &        
         hexvar % yLocal % array,          &
         hexvar % wachspressA % array,     &   
         hexvar % wachspressB % array,     &
         hexvar % wachspressKappa % array)

    call calc_cell_vertices_at_vertex(mesh, &
         hexvar % cellVerticesAtVertex % array)

    iCell = 1
    iVertex = 1

    call plot_wachpress(mesh, &
         iCell, &
         mesh % nEdgesOnCell % array(iCell), &
         iVertex, &
         hexvar % wachspressKappa % array(:,:,iCell), &
         hexvar % wachspressA % array(:,iCell), &
         hexvar % wachspressB % array(:,iCell), &
         hexvar % xLocal % array(:,iCell), &
         hexvar % yLocal % array(:,iCell))

  end subroutine init_dynamics_hex_variational
  
  !-------------------------------------------------------------

  subroutine calc_local_coords(mesh, xLocal, yLocal)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xLocal, &
         yLocal

    integer :: &
         iCell, &
         iVertex, &
         iVertexOnCell

    do iCell = 1, mesh % nCells

       do iVertex = 1, mesh % nEdgesOnCell % array(iCell)

          iVertexOnCell = mesh % verticesOnCell % array(iVertex, iCell)

          xLocal(iVertex,iCell) = mesh % xVertex % array(iVertexOnCell) - mesh % xCell % array(iCell)
          yLocal(iVertex,iCell) = mesh % yVertex % array(iVertexOnCell) - mesh % yCell % array(iCell)

       enddo ! iVertex

    enddo ! iCell

  end subroutine calc_local_coords

  !-------------------------------------------------------------

  subroutine calc_wachspress_coefficients(mesh, wachspressKappa, wachspressA, wachspressB, xLocal, yLocal)
    
    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         wachspressKappa

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         xLocal, &
         yLocal

    integer :: &
         iCell, &
         iVertex, &
         i0, &
         i1, &
         i2, &
         jVertex, &
         nEdgesOnCell

    ! loop over cells
    do iCell = 1, mesh % nCells

       ! number of edges
       nEdgesOnCell = mesh % nEdgesOnCell % array(iCell)
       
       ! loop over vertices
       do iVertex = 1, nEdgesOnCell
          
          ! end points of line segment
          i1 = iVertex - 1
          i2 = iVertex
          if (i1 < 1) i1 = i1 + nEdgesOnCell
          
          ! solve for the line segment equation
          wachspressA(iVertex, iCell) = (yLocal(i2,iCell) - yLocal(i1,iCell)) / (xLocal(i1,iCell) * yLocal(i2,iCell) - xLocal(i2,iCell) * yLocal(i1,iCell))
          wachspressB(iVertex, iCell) = (xLocal(i1,iCell) - xLocal(i2,iCell)) / (xLocal(i1,iCell) * yLocal(i2,iCell) - xLocal(i2,iCell) * yLocal(i1,iCell))

       enddo ! iVertex

       ! loop over vertices
       do iVertex = 1, nEdgesOnCell

          ! determine kappa
          wachspressKappa(1,iVertex,iCell) = 1.0_RKIND

          do jVertex = 2, nEdgesOnCell
             
             ! previous, this and next vertex
             i0 = jVertex - 1
             i1 = jVertex
             i2 = jVertex + 1
             if (i2 > nEdgesOnCell) i2 = i2 - nEdgesOnCell
             
             wachspressKappa(jVertex,iVertex,iCell) = wachspressKappa(jVertex-1,iVertex,iCell) * &
                  (wachspressA(i2,iCell) * (xLocal(i0,iCell) - xLocal(i1,iCell)) + wachspressB(i2,iCell) * (yLocal(i0,iCell) - yLocal(i1,iCell))) / &
                  (wachspressA(i0,iCell) * (xLocal(i1,iCell) - xLocal(i0,iCell)) + wachspressB(i0,iCell) * (yLocal(i1,iCell) - yLocal(i0,iCell)))

          enddo ! jVertex
          
       enddo ! iVertex
       
    enddo ! iCell

  end subroutine calc_wachspress_coefficients

  !-------------------------------------------------------------

  function wachspress_basis_function(nEdgesOnCell, iVertex, x, y, wachspressKappa, wachspressA, wachspressB) result(wachpress)

    integer, intent(in) :: &
         nEdgesOnCell, &
         iVertex

    real(kind=RKIND), intent(in) :: &
         x, &
         y
    
    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND) :: &
         wachpress

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

  end function wachspress_basis_function

  !-------------------------------------------------------------

  function wachspress_basis_derivative(nEdgesOnCell, iVertex, x, y, wachspressKappa, wachspressA, wachspressB, iDerivativeType) result(wachspress)

    integer, intent(in) :: &
         nEdgesOnCell, &
         iVertex, &
         iDerivativeType

    real(kind=RKIND), intent(in) :: &
         x, &
         y
    
    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND) :: &
         wachspress

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

  end function wachspress_basis_derivative

  !-------------------------------------------------------------

  function wachspress_numerator(nEdgesOnCell, jVertex, iVertex, x, y, wachspressKappa, wachspressA, wachspressB) result(numerator)

    integer, intent(in) :: &
         nEdgesOnCell, &
         jVertex, &
         iVertex

    real(kind=RKIND), intent(in) :: &
         x, &
         y
    
    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND) :: numerator

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

  end function wachspress_numerator

  !-------------------------------------------------------------

  function wachspress_numerator_derivative(nEdgesOnCell, jVertex, iVertex, x, y, wachspressKappa, wachspressA, wachspressB, iDerivativeType) result(derivative)

    integer, intent(in) :: &
         nEdgesOnCell, &
         jVertex, &
         iVertex, &
         iDerivativeType

    real(kind=RKIND), intent(in) :: &
         x, &
         y
    
    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND) :: &
         derivative

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

  end function wachspress_numerator_derivative

  !-------------------------------------------------------------

  function wachspress_edge_equation(iVertex, x, y, wachspressA, wachspressB) result(edge_equation)

    integer, intent(in) :: &
         iVertex

    real(kind=RKIND), intent(in) :: &
         x, &
         y

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND) :: edge_equation

    edge_equation = 1.0_RKIND - wachspressA(iVertex) * x - wachspressB(iVertex) * y

  end function wachspress_edge_equation

  !-------------------------------------------------------------

  function wachspress_edge_equation_derivative(iVertex, wachspressA, wachspressB, iDerivativeType) result(derivative)

    integer, intent(in) :: &
         iVertex, &
         iDerivativeType

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND) :: derivative

    if (iDerivativeType == 1) then

       derivative = -wachspressA(iVertex)

    else if (iDerivativeType == 2) then

       derivative = -wachspressB(iVertex)

    endif

  end function wachspress_edge_equation_derivative

  !-------------------------------------------------------------

  subroutine calculate_wachspress_derivatives(mesh, &
                                              basisGradientU, basisGradientV, &
                                              xLocal,         yLocal,         &
                                              wachspressA,    wachspressB,    &
                                              wachspressKappa)

    ! basisGradientUV(jVertexOnCell,iVertexOnCell,iCell)
    ! iCell         : The cell the gradients are based in
    ! iVertexOnCell : The vertex basis function the gradient is calculated from
    ! jVertexOnCell : The vertex location the gradients are calculated at
    
    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         basisGradientU, &
         basisGradientV

    real(kind=RKIND), dimension(:,:), intent(in) :: &    
         wachspressA, &
         wachspressB, &
         xLocal,      &
         yLocal

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &  
         wachspressKappa

    integer :: &
         nEdgesOnCell

    integer :: &
         iCell, &
         iVertexOnCell, &
         jVertexOnCell, &
         jVertex

    ! loop over cells
    do iCell = 1, mesh % nCells

       nEdgesOnCell = mesh % nEdgesOnCell % array(iCell)

       ! loop over vertices - basis function
       do iVertexOnCell = 1, nEdgesOnCell

          ! loop over vertices again - derivative position
          do jVertexOnCell = 1, nEdgesOnCell

             jVertex = mesh % verticesOnCell % array(jVertexOnCell,iCell)
          
             basisGradientU(jVertexOnCell,iVertexOnCell,iCell) = &
                  wachspress_basis_derivative(nEdgesOnCell, &
                                              iVertexOnCell, &
                                              xLocal(jVertexOnCell,iCell), &
                                              yLocal(jVertexOnCell,iCell), &
                                              wachspressKappa(:,:,iCell), &
                                              wachspressA(:,iCell), &
                                              wachspressB(:,iCell), &
                                              1)

             basisGradientV(jVertexOnCell,iVertexOnCell,iCell) = &
                  wachspress_basis_derivative(nEdgesOnCell, &
                                              iVertexOnCell, &
                                              xLocal(jVertexOnCell,iCell), &
                                              yLocal(jVertexOnCell,iCell), &
                                              wachspressKappa(:,:,iCell), &
                                              wachspressA(:,iCell), &
                                              wachspressB(:,iCell), &
                                              2)

           enddo ! jVertexOnCell

       enddo ! iVertexOnCell

    enddo ! iCell

  end subroutine calculate_wachspress_derivatives

  !-------------------------------------------------------------

  subroutine get_triangle_mapping(mapping, x1, y1, x2, y2, u1, v1, u2, v2)

    real(kind=RKIND), dimension(2,2), intent(out) :: &
         mapping

    real(kind=RKIND), intent(in) :: &
         x1, & ! input
         y1, &
         x2, & 
         y2, &
         u1, & ! output
         v1, &
         u2, &
         v2

    mapping(1,1) = (u2*y1 - u1*y2) / (x2*y1 - x1*y2)
    mapping(1,2) = (u1*x2 - u2*x1) / (y1*x2 - y2*x1)

    mapping(2,1) = (v2*y1 - v1*y2) / (x2*y1 - x1*y2)
    mapping(2,2) = (v1*x2 - v2*x1) / (y1*x2 - y2*x1)

  end subroutine get_triangle_mapping

  !-------------------------------------------------------------

  subroutine use_triangle_mapping(u, v, x, y, mapping)

    real(kind=RKIND), intent(out) :: &
         u, v

    real(kind=RKIND), intent(in) :: &
         x, y

    real(kind=RKIND), dimension(2,2), intent(in) :: &
         mapping

    u = mapping(1,1) * x + mapping(1,2) * y
    v = mapping(2,1) * x + mapping(2,2) * y

  end subroutine use_triangle_mapping

  !-------------------------------------------------------------

  subroutine integrate_wachspress(mesh, &
                                  basisIntegralsU, basisIntegralsV, &
                                  xLocal,          yLocal,          &
                                  wachspressA,     wachspressB,     &
                                  wachspressKappa)

    ! basisIntegralsUV (iVertexOnCell,jVertexOnCell,iCell)
    ! iCell         : cell integrals are performed on
    ! iVertexOnCell : vertex number of Wachspress function 
    ! jVertexOnCell : vertex number of Wachspress derivative function 
    ! Sij

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         basisIntegralsU, &
         basisIntegralsV

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         xLocal, &
         yLocal, &
         wachspressA, &
         wachspressB

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &    
         wachspressKappa

    type(hexvar_type), pointer :: hexvar

    real(kind=RKIND) :: &
         integration

    integer :: &
         iCell, &
         iVertexOnCell, &
         jVertexOnCell, &
         iDerivativeType, &
         nEdgesOnCell

    do iCell = 1, mesh % nCells

       if (iCell == 1) then

       nEdgesOnCell = mesh % nEdgesOnCell % array(iCell)

       do iVertexOnCell = 1, nEdgesOnCell

          do jVertexOnCell = 1, nEdgesOnCell

             do iDerivativeType = 1, 2

                call integrate_wachspress_polygon(&
                     integration,                nEdgesOnCell,         &
                     iVertexOnCell,              jvertexOnCell,        &
                     xLocal(:,iCell),            yLocal(:,iCell),      &
                     wachspressA(:,iCell),       wachspressB(:,iCell), &
                     wachspressKappa(:,:,iCell), iDerivativeType)

                if (iDerivativeType == 1) then

                   basisIntegralsU(iVertexOnCell,jVertexOnCell,iCell) = integration

                else if (iDerivativeType == 2) then

                   basisIntegralsV(iVertexOnCell,jVertexOnCell,iCell) = integration

                endif

             enddo ! iDerivativeType

          enddo ! jVertex

       enddo ! iVertex

       else

          basisIntegralsU(:,:,iCell) = basisIntegralsU(:,:,1)
          basisIntegralsV(:,:,iCell) = basisIntegralsV(:,:,1)

       endif

    enddo ! iCell

  end subroutine integrate_wachspress

  !-------------------------------------------------------------

  subroutine integrate_wachspress_polygon(integration,     nEdgesOnCell,   &
                                          iVertexOnCell,   jvertexOnCell,  &
                                          xLocal,          yLocal,         &
                                          wachspressA,     wachspressB,    &
                                          wachspressKappa, iDerivativeType)

    integer, intent(in) :: &
         nEdgesOnCell, &
         iVertexOnCell, &
         jvertexOnCell, &
         iDerivativeType

    real(kind=RKIND), dimension(:), intent(in) :: &
         xLocal, &
         yLocal

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND) :: &
         integration, &
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

       call get_triangle_mapping(mapping, &
                                 1.0_RKIND, 0.0_RKIND, &
                                 0.0_RKIND, 1.0_RKIND, &
                                 xLocal(i1), yLocal(i1), &
                                 xLocal(i2), yLocal(i2))

       call integrate_wachspress_subtriangle(integration_subtriangle, nEdgesOnCell, iVertexOnCell, jvertexOnCell, &
                                             wachspressKappa(:,:), wachspressA(:), wachspressB(:), iDerivativeType, mapping)

       integration = integration + integration_subtriangle

    enddo ! iSubTriangle

  end subroutine integrate_wachspress_polygon

  !-------------------------------------------------------------

  subroutine integrate_wachspress_subtriangle(integration,nEdgesOnCell, iVertexOnCell, jvertexOnCell, &
       wachspressKappa, wachspressA, wachspressB, iDerivativeType, mapping)

    real(kind=RKIND), intent(out) :: &
         integration

    integer, intent(in) :: &
         nEdgesOnCell, &
         iVertexOnCell, &
         jvertexOnCell, &
         iDerivativeType

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND), dimension(2,2), intent(in) :: &
         mapping

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
             
             call use_triangle_mapping(x, y, u, v, mapping)
             
             jacobian = mapping(1,1) * mapping(2,2) - mapping(1,2) * mapping(2,1)
             
             ! area test
             !integration = integration + scaling * jacobian
             
             ! actual integration
             integration = integration + scaling * jacobian * &
                  wachspress_basis_function(nEdgesOnCell, iVertexOnCell, x, y, wachspressKappa, wachspressA, wachspressB) * &
                  wachspress_basis_derivative(nEdgesOnCell, jVertexOnCell, x, y, wachspressKappa, wachspressA, wachspressB, iDerivativeType)

          endif

       enddo ! i
    enddo ! j

    integration = integration / (6.0_RKIND * real(n,RKIND)**2)

  end subroutine integrate_wachspress_subtriangle

  !-------------------------------------------------------------
  ! time step
  !-------------------------------------------------------------
  
  subroutine strain_tensor_hex_variational(mesh, &
                                           strain11,       strain22,       &
                                           strain12,                       &
                                           uVelocity,      vVelocity,      &
                                           basisGradientU, basisGradientV, &
                                           solveStress)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         strain11, &
         strain22, &
         strain12
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         basisGradientU, &
         basisGradientV
    
    integer, dimension(:), intent(in) :: &
         solveStress

    integer :: &
         iCell, &
         jVertexOnCell, &
         iVertexOnCell, &
         iVertex, &
         jVertex

    ! loop over cells
    do iCell = 1, mesh % nCells

       strain11(:,iCell) = 0.0_RKIND
       strain22(:,iCell) = 0.0_RKIND
       strain12(:,iCell) = 0.0_RKIND

       if (solveStress(iCell) == 1) then

          ! loop over velocity points surrounding cell - location of stress and derivative
          do jVertexOnCell = 1, mesh % nEdgesOnCell % array(iCell)

             jVertex = mesh % verticesOnCell % array(jVertexOnCell,iCell)

             ! loop over basis functions
             do iVertexOnCell = 1, mesh % nEdgesOnCell % array(iCell)

                iVertex = mesh % verticesOnCell % array(iVertexOnCell,iCell)

                strain11(jVertexOnCell,iCell) = strain11(jVertexOnCell,iCell) + &
                     uVelocity(iVertex) * basisGradientU(jVertexOnCell,iVertexOnCell,iCell)

                strain22(jVertexOnCell,iCell) = strain22(jVertexOnCell,iCell) + &
                     vVelocity(iVertex) * basisGradientV(jVertexOnCell,iVertexOnCell,iCell)

                strain12(jVertexOnCell,iCell) = strain12(jVertexOnCell,iCell) + 0.5_RKIND * (&
                     uVelocity(iVertex) * basisGradientV(jVertexOnCell,iVertexOnCell,iCell) + &
                     vVelocity(iVertex) * basisGradientU(jVertexOnCell,iVertexOnCell,iCell))

             enddo ! iVertexOnCell

          enddo ! jVertexOnCell

       endif ! solveStress

    enddo ! iCell

  end subroutine strain_tensor_hex_variational

  !-------------------------------------------------------------

  subroutine stress_tensor_hex_variational(mesh, &
                                           stress11, stress22, stress12,     &
                                           strain11, strain22, strain12,     &
                                           icePressure, replacementPressure, &
                                           solveStress,                      &
                                           dtElastic,   evpDamping)

    use cice_dynamics_shared, only: evp_constitutive_relation, &
                                    linear_constitutive_relation

    use cice_testing, only: iCellTest

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND), dimension(:), intent(out) :: &
         replacementPressure

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND), dimension(:), intent(in) :: &
         icePressure

    integer, dimension(:), intent(in) :: &
         solveStress

    real(kind=RKIND), intent(in) :: &
         dtElastic

    logical, intent(in) :: &
         evpDamping

    real(kind=RKIND) :: &
         replacementPressureCell

    integer :: &
         iCell, &
         iVertexOnCell

    do iCell = 1, mesh % nCells

       if (solveStress(iCell) == 1) then

          do iVertexOnCell = 1, mesh % nEdgesOnCell % array(iCell)

             call evp_constitutive_relation(stress11(iVertexOnCell,iCell), stress22(iVertexOnCell,iCell), stress12(iVertexOnCell,iCell), &
                  strain11(iVertexOnCell,iCell), strain22(iVertexOnCell,iCell), strain12(iVertexOnCell,iCell), &
                  icePressure(iCell), replacementPressureCell, &
                  mesh % areaCell % array(iCell),              &
                  dtElastic, evpDamping)

             if (iVertexOnCell == 1) then

                replacementPressure(iCell) = replacementPressureCell

             endif

             !call linear_constitutive_relation(stress11(iVertexOnCell,iCell), stress22(iVertexOnCell,iCell), stress12(iVertexOnCell,iCell), &
             !                                  strain11(iVertexOnCell,iCell), strain22(iVertexOnCell,iCell), strain12(iVertexOnCell,iCell))

          enddo ! iVertexOnCell

       else

          stress11(:,iCell) = 0.0_RKIND
          stress22(:,iCell) = 0.0_RKIND
          stress12(:,iCell) = 0.0_RKIND

       endif ! solveStress

    enddo ! iCell

  end subroutine stress_tensor_hex_variational

  !-------------------------------------------------------------

  subroutine stress_divergence_hex_variational(mesh, &
                                               stressDivergenceU,    stressDivergenceV, &
                                               stress11,             stress22,         &
                                               stress12,                                &
                                               basisIntegralsU,      basisIntegralsV,   &
                                               cellVerticesAtVertex, solveVelocity)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, &
         stressDivergenceV

    real(kind=RKIND), dimension(:,:), intent(in) :: &    
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &   
         basisIntegralsU, &
         basisIntegralsV

    integer, dimension(:,:), intent(in) :: &
         cellVerticesAtVertex

    integer, dimension(:), intent(in) :: &
         solveVelocity

    real(kind=RKIND) :: &
         stressDivergenceUCell, &
         stressDivergenceVCell

    integer :: &
         iVertex, &
         iCellOnVertex, &
         iCell, &
         iVertexOnCell, &
         jVertexOnCell

    ! loop over velocity positions
    do iVertex = 1, mesh % nVertices

       stressDivergenceU(iVertex) = 0.0_RKIND
       stressDivergenceV(iVertex) = 0.0_RKIND

       if (solveVelocity(iVertex) == 1) then

          ! loop over surrounding cells
          do iCellOnVertex = 1, mesh % vertexDegree
             
             ! get the cell number of this cell
             iCell = mesh % cellsOnVertex % array(iCellOnVertex, iVertex)
             
             ! get the vertexOnCell number of the iVertex velocity point from cell iCell
             jVertexOnCell = cellVerticesAtVertex(iCellOnVertex,iVertex)

             stressDivergenceUCell = 0.0_RKIND
             stressDivergenceVCell = 0.0_RKIND

             ! loop over the vertices of the surrounding cell
             do iVertexOnCell = 1, mesh % nEdgesOnCell % array(iCell)

                stressDivergenceUCell = stressDivergenceUCell + &
                     stress11(iVertexOnCell,iCell) * basisIntegralsU(iVertexOnCell,jVertexOnCell,iCell) + &
                     stress12(iVertexOnCell,iCell) * basisIntegralsV(iVertexOnCell,jVertexOnCell,iCell)
                
                stressDivergenceVCell = stressDivergenceVCell + &
                     stress22(iVertexOnCell,iCell) * basisIntegralsV(iVertexOnCell,jVertexOnCell,iCell) + &
                     stress12(iVertexOnCell,iCell) * basisIntegralsU(iVertexOnCell,jVertexOnCell,iCell)
                
             enddo ! jVertex

             stressDivergenceU(iVertex) = stressDivergenceU(iVertex) - stressDivergenceUCell
             stressDivergenceV(iVertex) = stressDivergenceV(iVertex) - stressDivergenceVCell
             
          enddo ! iCellOnVertex

          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) / mesh % areaTriangle % array(iVertex)
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) / mesh % areaTriangle % array(iVertex)

       endif ! solveVelocity

    enddo ! iVertex

  end subroutine stress_divergence_hex_variational

  !-------------------------------------------------------------
  ! Hex specific plotting stuff
  !-------------------------------------------------------------

  subroutine plot_wachpress(mesh, iCell, nEdgesOnCell, iVertex, wachspressKappa, wachspressA, wachspressB, xLocal, yLocal)
    
    type(mesh_type), intent(in) :: mesh

    integer, intent(in) :: &
         iCell, &
         nEdgesOnCell, &
         iVertex

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         wachspressKappa

    real(kind=RKIND), dimension(:), intent(out) :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND), dimension(:), intent(in) :: &
         xLocal, &
         yLocal

    real(kind=RKIND) :: &
         xmin, xmax, ymin, ymax
    
    real(kind=RKIND), dimension(2,2) :: &
         mapping

    integer :: &
         iSubTriangle, &
         i1, &
         i2, &
         iVertexOnCell, &
         iVertex2

    xmin =  1e30
    xmax = -1e30
    ymin =  1e30
    ymax = -1e30

    do iVertexOnCell = 1, nEdgesOnCell

       iVertex2 = mesh % verticesOnCell % array(iVertexOnCell,iCell)

       xmin = min(xmin,mesh % xVertex % array(iVertex2)) - 1000.0_RKIND
       xmax = max(xmax,mesh % xVertex % array(iVertex2)) + 1000.0_RKIND
       ymin = min(ymin,mesh % yVertex % array(iVertex2)) - 1000.0_RKIND
       ymax = max(ymax,mesh % yVertex % array(iVertex2)) + 1000.0_RKIND
       
       open(55,file='test.txt')
       if (iVertexOnCell == iVertex) then
          write(55,*) mesh % xVertex % array(iVertex2), mesh % yVertex % array(iVertex2)
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
       call plot_subtriangle(nEdgesOnCell, iSubTriangle, iVertex, &
                             wachspressKappa, wachspressA, wachspressB, &
                             mapping, mesh % xCell % array(iCell), mesh % yCell % array(iCell))

    enddo ! iSubTriangle

    close(55)
stop
  end subroutine plot_wachpress

  !-------------------------------------------------------------
  
  subroutine plot_subtriangle(nEdgesOnCell, iSubTriangle, iVertex, wachspressKappa, wachspressA, wachspressB, mapping, x0, y0)

    use cice_testing, only: matlab_jet

    integer, intent(in) :: &
         nEdgesOnCell, &
         iSubTriangle, &
         iVertex

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND), dimension(2,2), intent(in) :: &
         mapping

    real(kind=RKIND), intent(in) :: &
         x0, y0


    real(kind=RKIND) :: &
         u, v, &
         x, y, &
         wachspress, &
         wachspress1

    integer :: &
         i, j, iObject

    real(kind=RKIND) :: &
         x1, x2, x3, x4, x5, &
         y1, y2, y3, y4, y5, &
         d, minf, maxf, fValue, &
         dx, dy
    
    !real(kind=RKIND), parameter :: &
    !     fMin = -1.7759450831950581E-004_RKIND, &
    !     fMax = 5.5313696223098112E-005_RKIND
    real(kind=RKIND), parameter :: &
         fMin = 0.0_RKIND, &
         fMax = 1.0_RKIND

    integer, parameter :: n = 50

    character(len=7) :: &
         color

    logical, parameter :: &
         lfinitedifference = .false.

    iObject = ((iSubTriangle - 1) * (n+1) * (n+2)) / 2

    d = 1.0_RKIND / real(n, RKIND)
 
    minf = 1e30
    maxf = -1e30

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

             fValue = (wachspress - fMin) / (fMax - fMin)
                
             minf = min(minf,wachspress)
             maxf = max(maxf,wachspress)

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

                color = matlab_jet(max(min(fValue,1.0_RKIND),0.0_RKIND))
                write(55,fmt='(a,i5,a,a,a)') "set object ",iObject,' fc rgb "', color, '" fillstyle solid'  

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

                color = matlab_jet(max(min(fValue,1.0_RKIND),0.0_RKIND))
                write(55,fmt='(a,i5,a,a,a)') "set object ",iObject,' fc rgb "', color, '" fillstyle solid'   

             endif

          endif

       enddo ! j
    enddo ! i

    write(*,*) minf, maxf

  end subroutine plot_subtriangle

  !-------------------------------------------------------------

  subroutine plot_wachspress_square(mesh,wachspressKappa, wachspressA, wachspressB)
    type(mesh_type), intent(in) :: mesh
    real(kind=RKIND), dimension(:,:,:) :: wachspressKappa
    real(kind=RKIND), dimension(:,:) :: wachspressA, wachspressB

    integer :: &
         iCell, &
         nEdgesOnCell, &
         iVertex, &
         iVertexOnCell

    real(kind=RKIND) :: &
         x, y, x0, y0, dx, dy, wachspress

    integer :: i, j

    integer, parameter :: ni = 100
    integer, parameter :: nj = 100

    iCell = 1

    nEdgesOnCell = mesh % nEdgesOnCell % array(iCell)

    x0 = 1e30
    y0 = 1e30
    
    do iVertexOnCell = 1, nEdgesOnCell
       
       iVertex = mesh % verticesOnCell % array(iVertexOnCell,iCell)

       x0 = min(x0,mesh % xVertex % array(iVertex))
       y0 = min(y0,mesh % yVertex % array(iVertex))

    enddo ! 
    
    dx = 16000.0_RKIND / real(ni,RKIND)
    dy = 16000.0_RKIND / real(nj,RKIND)
       
    iVertexOnCell = 1


    do i = 1, ni
       do j = 1, nj

          x = x0 + (real(i,RKIND) - 0.5_RKIND) * dx - 24000.0_RKIND
          y = x0 + (real(j,RKIND) - 0.5_RKIND) * dy - 24000.0_RKIND

          wachspress =  &
               wachspress_basis_derivative(nEdgesOnCell, iVertexOnCell, x, y, wachspressKappa(:,:,iCell), wachspressA(:,iCell), wachspressB(:,iCell),1)

          write(55,*) x+0.5_RKIND*dx,x-0.5_RKIND*dx,y+0.5_RKIND*dy,y-0.5_RKIND*dy, wachspress          

       enddo ! j
    enddo ! i

  end subroutine plot_wachspress_square

  !-------------------------------------------------------------


end module cice_dyn_hex_variational
