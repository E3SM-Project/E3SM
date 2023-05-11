!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_velocity_solver_wachspress
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_velocity_solver_wachspress

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timer
  use mpas_log, only: mpas_log_write

  implicit none

  private
  save

  public :: &
       seaice_init_velocity_solver_wachspress

contains

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_velocity_solver_wachspress
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_velocity_solver_wachspress(domain)!{{{

    use mpas_timer

    use seaice_mesh, only: &
         seaice_calc_local_coords

    use seaice_wachspress_basis, only: &
         seaice_calc_wachspress_coefficients

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         blockPtr

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocityVariationalPool

    integer, pointer :: &
         nCells, &
         maxEdges

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         zVertex, &
         xCell, &
         yCell, &
         zCell

    real(kind=RKIND), dimension(:,:), allocatable :: &
         xLocal, &
         yLocal

    logical, pointer :: &
         config_rotate_cartesian_grid, &
         config_include_metric_terms, &
         on_a_sphere

    character(len=strKIND), pointer :: &
         config_wachspress_integration_type

    integer, pointer :: &
         config_wachspress_integration_order

    real(kind=RKIND), pointer :: &
         sphere_radius

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisGradientU, &
         basisGradientV, &
         basisIntegralsU, &
         basisIntegralsV, &
         basisIntegralsMetric

    real(kind=RKIND), dimension(:,:), allocatable :: &
         wachspressA, &
         wachspressB

    real(kind=RKIND), dimension(:,:,:), allocatable :: &
         wachspressKappa

    integer :: &
         iCell, &
         iStressVertex, &
         iVelocityVertex, &
         iSubCell, &
         iVertex

    call mpas_timer_start("Velocity solver Wachpress init")

    call MPAS_pool_get_config(domain % configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)
    call MPAS_pool_get_config(domain % configs, "config_include_metric_terms", config_include_metric_terms)
    call MPAS_pool_get_config(domain % configs, "config_wachspress_integration_type", config_wachspress_integration_type)
    call MPAS_pool_get_config(domain % configs, "config_wachspress_integration_order", config_wachspress_integration_order)

    blockPtr => domain % blocklist
    do while (associated(blockPtr))

       call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_variational", velocityVariationalPool)

       call MPAS_pool_get_config(meshPool, "on_a_sphere", on_a_sphere)
       call MPAS_pool_get_config(meshPool, "sphere_radius", sphere_radius)

       call MPAS_pool_get_dimension(meshPool, "nCells", nCells)
       call MPAS_pool_get_dimension(meshPool, "maxEdges", maxEdges)

       call MPAS_pool_get_array(meshPool, "nEdgesOnCell", nEdgesOnCell)
       call MPAS_pool_get_array(meshPool, "verticesOnCell", verticesOnCell)
       call MPAS_pool_get_array(meshPool, "xVertex", xVertex)
       call MPAS_pool_get_array(meshPool, "yVertex", yVertex)
       call MPAS_pool_get_array(meshPool, "zVertex", zVertex)
       call MPAS_pool_get_array(meshPool, "xCell", xCell)
       call MPAS_pool_get_array(meshPool, "yCell", yCell)
       call MPAS_pool_get_array(meshPool, "zCell", zCell)

       call MPAS_pool_get_array(velocityVariationalPool, "basisGradientU", basisGradientU)
       call MPAS_pool_get_array(velocityVariationalPool, "basisGradientV", basisGradientV)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsU", basisIntegralsU)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsV", basisIntegralsV)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsMetric", basisIntegralsMetric)

       call mpas_timer_start("variational calc_local_coords")
       allocate(xLocal(maxEdges,nCells))
       allocate(yLocal(maxEdges,nCells))

       call seaice_calc_local_coords(&
            xLocal, &
            yLocal, &
            nCells, &
            nEdgesOnCell, &
            verticesOnCell, &
            xVertex, &
            yVertex, &
            zVertex, &
            xCell, &
            yCell, &
            zCell, &
            config_rotate_cartesian_grid, &
            on_a_sphere)
       call mpas_timer_stop("variational calc_local_coords")

       allocate(wachspressKappa(maxEdges,maxEdges,nCells))
       allocate(wachspressA(maxEdges,nCells))
       allocate(wachspressB(maxEdges,nCells))

       call mpas_timer_start("wachpress calc_coefficients")
       call seaice_calc_wachspress_coefficients(&
            wachspressKappa, &
            wachspressA, &
            wachspressB, &
            nCells, &
            nEdgesOnCell, &
            xLocal, &
            yLocal)
       call mpas_timer_stop("wachpress calc_coefficients")

       call mpas_timer_start("wachpress calc_derivatives")
       call calculate_wachspress_derivatives(&
            basisGradientU, &
            basisGradientV, &
            nCells, &
            maxEdges, &
            nEdgesOnCell, &
            xLocal, &
            yLocal, &
            wachspressA, &
            wachspressB, &
            wachspressKappa)
       call mpas_timer_stop("wachpress calc_derivatives")

       call mpas_timer_start("wachpress integrate")
       call integrate_wachspress(&
            basisIntegralsU, &
            basisIntegralsV, &
            basisIntegralsMetric, &
            nCells, &
            nEdgesOnCell, &
            xLocal, &
            yLocal, &
            wachspressA, &
            wachspressB, &
            wachspressKappa, &
            config_wachspress_integration_type, &
            config_wachspress_integration_order)
       call mpas_timer_stop("wachpress integrate")

       deallocate(wachspressKappa)
       deallocate(wachspressA)
       deallocate(wachspressB)
       deallocate(xLocal)
       deallocate(yLocal)

       blockPtr => blockPtr % next
    enddo

    call mpas_timer_stop("Velocity solver Wachpress init")

  end subroutine seaice_init_velocity_solver_wachspress!}}}

!-----------------------------------------------------------------------
! Integration
!-----------------------------------------------------------------------

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
       basisIntegralsU, &
       basisIntegralsV, &
       basisIntegralsMetric, &
       nCells, &
       nEdgesOnCell, &
       xLocal, &
       yLocal, &
       wachspressA, &
       wachspressB, &
       wachspressKappa, &
       integrationType, &
       integrationOrder)!{{{

    use seaice_triangle_quadrature, only: &
         seaice_triangle_quadrature_rules

    ! basisIntegralsUV (iStressVertex,iVelocityVertex,iCell)
    ! iCell         : cell integrals are performed on
    ! iStressVertex : vertex number of Wachspress function
    ! iVelocityVertex : vertex number of Wachspress derivative function
    ! Sij

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         basisIntegralsU, &   !< Output:
         basisIntegralsV, &   !< Output:
         basisIntegralsMetric !< Output:

    integer, intent(in) :: &
         nCells !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         xLocal, &      !< Input:
         yLocal, &      !< Input:
         wachspressA, & !< Input:
         wachspressB    !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         wachspressKappa !< Input:

    character(len=strKIND), intent(in) :: &
         integrationType !< Input:

    integer, intent(in) :: &
         integrationOrder !< Input:

    integer :: &
         iCell, &
         iStressVertex, &
         iVelocityVertex

    integer :: &
         nIntegrationPoints

    real(kind=RKIND), dimension(:), allocatable :: &
         integrationU, &
         integrationV, &
         integrationWeights

    real(kind=RKIND) :: &
         normalizationFactor

    ! Quadrature rules
    call seaice_triangle_quadrature_rules(&
         integrationType, &
         integrationOrder, &
         nIntegrationPoints, &
         integrationU, &
         integrationV, &
         integrationWeights, &
         normalizationFactor)

    !$omp parallel do default(shared) private(iStressVertex, iVelocityVertex)
    do iCell = 1, nCells

       do iVelocityVertex = 1, nEdgesOnCell(iCell)

          do iStressVertex = 1, nEdgesOnCell(iCell)

             basisIntegralsU(iStressVertex,iVelocityVertex,iCell)      = 0.0_RKIND
             basisIntegralsV(iStressVertex,iVelocityVertex,iCell)      = 0.0_RKIND
             basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) = 0.0_RKIND

             call integrate_wachspress_polygon(&
                  basisIntegralsU(iStressVertex,iVelocityVertex,iCell), &
                  basisIntegralsV(iStressVertex,iVelocityVertex,iCell), &
                  basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell), &
                  nEdgesOnCell(iCell), &
                  iStressVertex, &
                  iVelocityVertex, &
                  xLocal(:,iCell), &
                  yLocal(:,iCell), &
                  wachspressA(:,iCell), &
                  wachspressB(:,iCell), &
                  wachspressKappa(:,:,iCell), &
                  nIntegrationPoints, &
                  integrationU, &
                  integrationV, &
                  integrationWeights, &
                  normalizationFactor)

          enddo ! jVertex

       enddo ! iVertex

    enddo ! iCell

    deallocate(integrationU)
    deallocate(integrationV)
    deallocate(integrationWeights)

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
       basisIntegralsU, &
       basisIntegralsV, &
       basisIntegralsMetric, &
       nEdgesOnCell, &
       iStressVertex, &
       iVelocityVertex, &
       xLocal, &
       yLocal, &
       wachspressA, &
       wachspressB, &
       wachspressKappa, &
       nIntegrationPoints, &
       integrationU, &
       integrationV, &
       integrationWeights, &
       normalizationFactor)!{{{

    use seaice_mesh, only: &
         seaice_wrapped_index

    use seaice_wachspress_basis, only: &
         seaice_wachspress_indexes, &
         seaice_wachspress_basis_function, &
         seaice_wachspress_basis_derivative

    real(kind=RKIND), intent(inout) :: &
         basisIntegralsU, &   !< Input/Output:
         basisIntegralsV, &   !< Input/Output:
         basisIntegralsMetric !< Input/Output:

    integer, intent(in) :: &
         nEdgesOnCell, &  !< Input:
         iStressVertex, & !< Input:
         iVelocityVertex  !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xLocal, & !< Input:
         yLocal    !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         wachspressA, & !< Input:
         wachspressB    !< Input:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         wachspressKappa !< Input:

    integer, intent(in) :: &
         nIntegrationPoints !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         integrationU, & !< Input:
         integrationV, & !< Input:
         integrationWeights !< Input:

    real(kind=RKIND), intent(in) :: &
         normalizationFactor !< Input:

    integer, dimension(nEdgesOnCell) :: &
         nEdgesOnCellSubset

    integer, dimension(nEdgesOnCell,nEdgesOnCell) :: &
         vertexIndexSubset

    real(kind=RKIND) :: &
         basisIntegralsSubTriangleTmp, &
         basisIntegralsUSubTriangle, &
         basisIntegralsVSubTriangle, &
         basisIntegralsMetricSubTriangle

    real(kind=RKIND), dimension(nIntegrationPoints) :: &
         x, &
         y, &
         stressBasisFunction, &
         velocityBasisFunction, &
         velocityBasisDerivativeU, &
         velocityBasisDerivativeV

    real(kind=RKIND), dimension(2,2) :: &
         mapping

    real(kind=RKIND), dimension(nEdgesOnCell) :: &
         jacobian

    integer :: &
         iIntegrationPoint, &
         iSubTriangle, &
         i1, &
         i2

    call seaice_wachspress_indexes(&
         nEdgesOnCell, &
         nEdgesOnCellSubset, &
         vertexIndexSubset)

    do iSubTriangle = 1, nEdgesOnCell

       i1 = iSubTriangle
       i2 = seaice_wrapped_index(iSubTriangle + 1, nEdgesOnCell)

       call get_triangle_mapping(&
            mapping, &
            jacobian(iSubTriangle), &
            1.0_RKIND, 0.0_RKIND, &
            0.0_RKIND, 1.0_RKIND, &
            xLocal(i1), yLocal(i1), &
            xLocal(i2), yLocal(i2))

       !in-lined use_triangle_mapping
       do iIntegrationPoint = 1, nIntegrationPoints

          x(iIntegrationPoint) = mapping(1,1) * integrationU(iIntegrationPoint) + &
                                 mapping(1,2) * integrationV(iIntegrationPoint)
          y(iIntegrationPoint) = mapping(2,1) * integrationU(iIntegrationPoint) + &
                                 mapping(2,2) * integrationV(iIntegrationPoint)

       enddo ! iIntegrationPoint

       call seaice_wachspress_basis_function(&
            nEdgesOnCell, iStressVertex, x, y, &
            wachspressKappa, wachspressA, wachspressB, &
            nEdgesOnCellSubset, vertexIndexSubset, &
            stressBasisFunction)

       call seaice_wachspress_basis_function(&
            nEdgesOnCell, iVelocityVertex, x, y, &
            wachspressKappa, wachspressA, wachspressB, &
            nEdgesOnCellSubset, vertexIndexSubset, &
            velocityBasisFunction)

       call seaice_wachspress_basis_derivative(&
            nEdgesOnCell, iVelocityVertex, x, y, &
            wachspressKappa, wachspressA, wachspressB, &
            nEdgesOnCellSubset, vertexIndexSubset, &
            velocityBasisDerivativeU, &
            velocityBasisDerivativeV)

       basisIntegralsUSubTriangle      = 0.0_RKIND
       basisIntegralsVSubTriangle      = 0.0_RKIND
       basisIntegralsMetricSubTriangle = 0.0_RKIND

       do iIntegrationPoint = 1, nIntegrationPoints

          basisIntegralsSubTriangleTmp = &
               jacobian(iSubTriangle) * &
               integrationWeights(iIntegrationPoint) * &
               stressBasisFunction(iIntegrationPoint)

          basisIntegralsUSubTriangle = basisIntegralsUSubTriangle + &
               basisIntegralsSubTriangleTmp * &
               velocityBasisDerivativeU(iIntegrationPoint)

          basisIntegralsVSubTriangle = basisIntegralsVSubTriangle + &
               basisIntegralsSubTriangleTmp * &
               velocityBasisDerivativeV(iIntegrationPoint)

          basisIntegralsMetricSubTriangle = basisIntegralsMetricSubTriangle + &
               basisIntegralsSubTriangleTmp * &
               velocityBasisFunction(iIntegrationPoint)

       enddo ! iIntegrationPoint

       basisIntegralsU      = basisIntegralsU      + basisIntegralsUSubTriangle      / normalizationFactor
       basisIntegralsV      = basisIntegralsV      + basisIntegralsVSubTriangle      / normalizationFactor
       basisIntegralsMetric = basisIntegralsMetric + basisIntegralsMetricSubTriangle / normalizationFactor

    enddo ! iSubTriangle

  end subroutine integrate_wachspress_polygon!}}}

!-----------------------------------------------------------------------
! Remapping
!-----------------------------------------------------------------------

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
       jacobian, &
       x1, y1, &
       x2, y2, &
       u1, v1, &
       u2, v2)!{{{

    real(kind=RKIND), dimension(2,2), intent(out) :: &
         mapping !< Output:

    real(kind=RKIND), intent(out) :: &
         jacobian !< Output:

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

    jacobian = mapping(1,1) * mapping(2,2) - mapping(1,2) * mapping(2,1)

  end subroutine get_triangle_mapping!}}}

!-----------------------------------------------------------------------
! Wachspress derivatives
!-----------------------------------------------------------------------

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
       basisGradientU, &
       basisGradientV, &
       nCells, &
       maxEdges, &
       nEdgesOnCell, &
       xLocal, &
       yLocal, &
       wachspressA, &
       wachspressB, &
       wachspressKappa)!{{{

    use seaice_mesh, only: &
         seaice_wrapped_index

    use seaice_wachspress_basis, only: &
         seaice_wachspress_indexes, &
         seaice_wachspress_basis_derivative

    ! basisGradientUV(jVertexOnCell,iVertexOnCell,iCell)
    ! iCell         : The cell the gradients are based in
    ! iVertexOnCell : The vertex basis function the gradient is calculated from
    ! jVertexOnCell : The vertex location the gradients are calculated at

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         basisGradientU, & !< Output:
         basisGradientV    !< Output:

    integer, intent(in) :: &
         nCells, & !< Input:
         maxEdges  !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell !< Input:

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

    integer, dimension(:), allocatable :: &
         nEdgesOnCellSubset

    integer, dimension(:,:), allocatable :: &
         vertexIndexSubset

    real(kind=RKIND), dimension(:), allocatable :: &
         x, y, derivativeU, derivativeV

    allocate(x(maxEdges))
    allocate(y(maxEdges))

    allocate(derivativeU(maxEdges))
    allocate(derivativeV(maxEdges))

    allocate(nEdgesOnCellSubset(maxEdges))
    allocate(vertexIndexSubset(maxEdges,maxEdges))

    ! loop over cells
    do iCell = 1, nCells

       call seaice_wachspress_indexes(&
            nEdgesOnCell(iCell), &
            nEdgesOnCellSubset(1:nEdgesOnCell(iCell)), &
            vertexIndexSubset(1:nEdgesOnCell(iCell),1:nEdgesOnCell(iCell)))

       ! loop over vertices again - derivative position
       do iGradientVertex = 1, nEdgesOnCell(iCell)

          x(iGradientVertex) = xLocal(iGradientVertex,iCell)
          y(iGradientVertex) = yLocal(iGradientVertex,iCell)

       enddo ! iGradientVertex

       ! loop over vertices - basis function
       do iBasisVertex = 1, nEdgesOnCell(iCell)

          call seaice_wachspress_basis_derivative(&
               nEdgesOnCell(iCell), &
               iBasisVertex, &
               x(1:nEdgesOnCell(iCell)), &
               y(1:nEdgesOnCell(iCell)), &
               wachspressKappa(:,:,iCell), &
               wachspressA(:,iCell), &
               wachspressB(:,iCell), &
               nEdgesOnCellSubset(1:nEdgesOnCell(iCell)), &
               vertexIndexSubset(1:nEdgesOnCell(iCell),1:nEdgesOnCell(iCell)), &
               derivativeU(1:nEdgesOnCell(iCell)), &
               derivativeV(1:nEdgesOnCell(iCell)))

          basisGradientU(iBasisVertex,:,iCell) = 0.0_RKIND
          basisGradientV(iBasisVertex,:,iCell) = 0.0_RKIND

          iGradientVertex = iBasisVertex
          basisGradientU(iBasisVertex,iGradientVertex,iCell) = derivativeU(iGradientVertex)
          basisGradientV(iBasisVertex,iGradientVertex,iCell) = derivativeV(iGradientVertex)

          iGradientVertex = seaice_wrapped_index(iBasisVertex - 1, nEdgesOnCell(iCell))
          basisGradientU(iBasisVertex,iGradientVertex,iCell) = derivativeU(iGradientVertex)
          basisGradientV(iBasisVertex,iGradientVertex,iCell) = derivativeV(iGradientVertex)

          iGradientVertex = seaice_wrapped_index(iBasisVertex + 1, nEdgesOnCell(iCell))
          basisGradientU(iBasisVertex,iGradientVertex,iCell) = derivativeU(iGradientVertex)
          basisGradientV(iBasisVertex,iGradientVertex,iCell) = derivativeV(iGradientVertex)

       enddo ! iBasisVertex

    enddo ! iCell

    deallocate(nEdgesOnCellSubset)
    deallocate(vertexIndexSubset)

    deallocate(x)
    deallocate(y)

    deallocate(derivativeU)
    deallocate(derivativeV)

  end subroutine calculate_wachspress_derivatives!}}}

!-----------------------------------------------------------------------

end module seaice_velocity_solver_wachspress
