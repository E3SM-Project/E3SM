!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_velocity_solver_weak
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module cice_velocity_solver_weak

  use mpas_derived_types
  use mpas_pool_routines

  implicit none

  private
  save

  public :: &
       cice_init_velocity_solver_weak, &
       cice_internal_stress_weak, &
       cice_strain_tensor_weak, &
       cice_stress_tensor_weak, &
       cice_stress_divergence_weak, &
       cice_final_divergence_shear_weak

contains

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_velocity_solver_weak
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_init_velocity_solver_weak(&
       mesh, &
       boundary, &
       velocity_weak, &
       rotateCartesianGrid)!{{{

    use cice_mesh, only: &
         cice_normal_vectors

    type(MPAS_pool_type), pointer, intent(inout) :: &
         mesh !< Input/Output:

    type(MPAS_pool_type), pointer :: &
         velocity_weak, & !< Input/Output:
         boundary   !< Input/Output:

    logical, intent(in) :: &
         rotateCartesianGrid !< Input:

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         normalVectorPolygon, &
         normalVectorTriangle

    real(kind=RKIND), dimension(:), pointer :: &
         latCellRotated, &
         latVertexRotated

    integer, dimension(:), pointer :: &
         interiorVertex

    call MPAS_pool_get_array(velocity_weak, "normalVectorPolygon", normalVectorPolygon)
    call MPAS_pool_get_array(velocity_weak, "normalVectorTriangle", normalVectorTriangle)
    call MPAS_pool_get_array(velocity_weak, "latCellRotated", latCellRotated)
    call MPAS_pool_get_array(velocity_weak, "latVertexRotated", latVertexRotated)
    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call cice_normal_vectors(&
         mesh, &
         normalVectorPolygon, &
         normalVectorTriangle, &
         interiorVertex, &
         rotateCartesianGrid, &
         .true., &
         latCellRotated, &
         latVertexRotated)

  end subroutine cice_init_velocity_solver_weak!}}}

!-----------------------------------------------------------------------
! Time step
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_internal_stress_weak
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_internal_stress_weak(domain)!{{{

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocityWeakPool, &
         velocitySolverPool

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         icePressure, &
         replacementPressure, &
         stressDivergenceU, &
         stressDivergenceV, &
         strain11, &
         strain22, &
         strain12, &
         stress11, &
         stress22, &
         stress12, &
         latCellRotated, &
         latVertexRotated

    real(kind=RKIND), pointer :: &
         elasticTimeStep

    logical, pointer :: &
         revisedEVP

    integer, dimension(:), pointer :: &
         solveStress, &
         solveVelocity

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         normalVectorPolygon, &
         normalVectorTriangle

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_revised_evp", revisedEVP)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocityWeakPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "icePressure", icePressure)
       call MPAS_pool_get_array(velocitySolverPool, "elasticTimeStep", elasticTimeStep)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)

       call MPAS_pool_get_array(velocityWeakPool, "normalVectorPolygon", normalVectorPolygon)
       call MPAS_pool_get_array(velocityWeakPool, "normalVectorTriangle", normalVectorTriangle)
       call MPAS_pool_get_array(velocityWeakPool, "latCellRotated", latCellRotated)
       call MPAS_pool_get_array(velocityWeakPool, "latVertexRotated", latVertexRotated)
       call MPAS_pool_get_array(velocityWeakPool, "strain11", strain11)
       call MPAS_pool_get_array(velocityWeakPool, "strain22", strain22)
       call MPAS_pool_get_array(velocityWeakPool, "strain12", strain12)
       call MPAS_pool_get_array(velocityWeakPool, "stress11", stress11)
       call MPAS_pool_get_array(velocityWeakPool, "stress22", stress22)
       call MPAS_pool_get_array(velocityWeakPool, "stress12", stress12)
       call MPAS_pool_get_array(velocityWeakPool, "replacementPressure", replacementPressure)

       call cice_strain_tensor_weak(&
            meshPool, &
            strain11, &
            strain22, &
            strain12, &
            uVelocity, &
            vVelocity, &
            normalVectorPolygon, &
            latCellRotated, &
            solveStress)

       call cice_stress_tensor_weak(&
            meshPool, &
            stress11, &
            stress22, &
            stress12, &
            strain11, &
            strain22, &
            strain12, &
            icePressure, &
            replacementPressure, &
            solveStress, &
            elasticTimeStep, &
            revisedEVP)

       call cice_stress_divergence_weak(&
            meshPool, &
            stressDivergenceU, &
            stressDivergenceV, &
            stress11, &
            stress22, &
            stress12, &
            normalVectorTriangle, &
            latVertexRotated, &
            solveVelocity)

       block => block % next
    end do

  end subroutine cice_internal_stress_weak!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_strain_tensor_weak
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_strain_tensor_weak(&
       mesh, &
       strain11, &
       strain22, &
       strain12, &
       uVelocity, &
       vVelocity, &
       normalVectorPolygon, &
       latCellRotated, &
       solveStress)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         strain11, & !< Output:
         strain22, & !< Output:
         strain12    !< Output:

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, &   !< Input:
         vVelocity, &   !< Input:
         latCellRotated !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorPolygon !< Input:

    integer, dimension(:), intent(in) :: &
         solveStress !< Input:

    integer :: &
         iCell, &
         iEdgeOnCell, &
         iEdge, &
         iVertexOnEdge, &
         iVertex

    real(kind=RKIND) :: &
         uVelocityEdge, &
         vVelocityEdge, &
         uCellCentre, &
         vCellCentre

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell, &
         edgesOnCell, &
         verticesOnEdge

    real(kind=RKIND), pointer :: &
         sphere_radius

    real(kind=RKIND) :: &
         sphereRadius

    real(kind=RKIND), dimension(:), pointer :: &
         dvEdge, &
         areaCell

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnEdge", verticesOnEdge)
    call MPAS_pool_get_array(mesh, "dvEdge", dvEdge)
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    ! planar cases with zero sphere radius
    sphereRadius = sphere_radius
    if (sphereRadius == 0.0_RKIND) sphereRadius = 1.0_RKIND

    do iCell = 1, nCells

       strain11(iCell) = 0.0_RKIND
       strain22(iCell) = 0.0_RKIND
       strain12(iCell) = 0.0_RKIND

       if (solveStress(iCell) == 1) then

          uCellCentre = 0.0_RKIND
          vCellCentre = 0.0_RKIND

          do iEdgeOnCell = 1, nEdgesOnCell(iCell)

             ! cell centre velocities
             iVertex = verticesOnCell(iEdgeOnCell,iCell)

             uCellCentre = uCellCentre + uVelocity(iVertex)
             vCellCentre = vCellCentre + vVelocity(iVertex)

             ! interpolated edge velocity
             iEdge = edgesOnCell(iEdgeOnCell,iCell)

             uVelocityEdge = 0.0_RKIND
             vVelocityEdge = 0.0_RKIND

             do iVertexOnEdge = 1, 2

                iVertex = verticesOnEdge(iVertexOnEdge,iEdge)

                uVelocityEdge = uVelocityEdge + uVelocity(iVertex)
                vVelocityEdge = vVelocityEdge + vVelocity(iVertex)

             enddo ! iVertexOnEdge

             uVelocityEdge = uVelocityEdge / 2.0_RKIND
             vVelocityEdge = vVelocityEdge / 2.0_RKIND

             ! summation over edges
             strain11(iCell) = strain11(iCell) + uVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) * dvEdge(iEdge)
             strain22(iCell) = strain22(iCell) + vVelocityEdge * normalVectorPolygon(2,iEdgeOnCell,iCell) * dvEdge(iEdge)
             strain12(iCell) = strain12(iCell) + 0.5_RKIND * ( &
                  uVelocityEdge * normalVectorPolygon(2,iEdgeOnCell,iCell) + &
                  vVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) ) * dvEdge(iEdge)

          enddo ! iEdgeOnCell

          uCellCentre = uCellCentre / real(nEdgesOnCell(iCell), RKIND)
          vCellCentre = vCellCentre / real(nEdgesOnCell(iCell), RKIND)

          strain11(iCell) = strain11(iCell) / areaCell(iCell)
          strain22(iCell) = strain22(iCell) / areaCell(iCell)
          strain12(iCell) = strain12(iCell) / areaCell(iCell)

          ! metric terms
          strain11(iCell) = strain11(iCell) - (vCellCentre * tan(latCellRotated(iCell))) / sphereRadius
          strain12(iCell) = strain12(iCell) + (uCellCentre * tan(latCellRotated(iCell)) * 0.5_RKIND) / sphereRadius

          !if (abs(strain11(iCell)) < 1e-10_RKIND) strain11(iCell) = 0.0_RKIND

       endif ! solveStress

    enddo ! iCell

  end subroutine cice_strain_tensor_weak!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_stress_tensor_weak
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_stress_tensor_weak(&
       mesh, &
       stress11, &
       stress22, &
       stress12, &
       strain11, &
       strain22, &
       strain12, &
       icePressure, &
       replacementPressure, &
       solveStress, &
       dtElastic, &
       revisedEVP)!{{{

    use cice_velocity_solver_constitutive_relation, only: &
         cice_evp_constitutive_relation, &
         cice_evp_constitutive_relation_revised

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(inout) :: &
         stress11, &         !< Input/Output:
         stress22, &         !< Input/Output:
         stress12, &         !< Input/Output:
         replacementPressure !< Input/Output:

    real(kind=RKIND), dimension(:), intent(inout) :: &
         strain11, & !< Input/Output:
         strain22, & !< Input/Output:
         strain12, & !< Input/Output:
         icePressure !< Input/Output:

    integer, dimension(:), intent(in) :: &
         solveStress !< Input:

    real(kind=RKIND), intent(in) :: &
         dtElastic !< Input:

    logical, intent(in) :: &
         revisedEVP !< Input:

    integer :: &
         iCell

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    if (.not. revisedEVP) then

       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             call cice_evp_constitutive_relation(&
                  stress11(iCell), &
                  stress22(iCell), &
                  stress12(iCell), &
                  strain11(iCell), &
                  strain22(iCell), &
                  strain12(iCell), &
                  icePressure(iCell), &
                  replacementPressure(iCell), &
                  areaCell(iCell), &
                  dtElastic)

          else

             stress11(iCell) = 0.0_RKIND
             stress22(iCell) = 0.0_RKIND
             stress12(iCell) = 0.0_RKIND

          endif ! solveStress

       end do ! iCell

    else

       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             call cice_evp_constitutive_relation_revised(&
                  stress11(iCell), &
                  stress22(iCell), &
                  stress12(iCell), &
                  strain11(iCell), &
                  strain22(iCell), &
                  strain12(iCell), &
                  icePressure(iCell), &
                  replacementPressure(iCell), &
                  areaCell(iCell))

          else

             stress11(iCell) = 0.0_RKIND
             stress22(iCell) = 0.0_RKIND
             stress12(iCell) = 0.0_RKIND

          endif ! solveStress

       end do ! iCell

    endif

  end subroutine cice_stress_tensor_weak!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_stress_divergence_weak
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_stress_divergence_weak(&
       mesh, &
       stressDivergenceU, &
       stressDivergenceV, &
       stress11, &
       stress22, &
       stress12, &
       normalVectorTriangle, &
       latVertexRotated,  &
       solveVelocity)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, & !< Output:
         stressDivergenceV    !< Output:

    real(kind=RKIND), dimension(:), intent(in) :: &
         stress11, &      !< Input:
         stress22, &      !< Input:
         stress12, &      !< Input:
         latVertexRotated !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorTriangle !< Input:

    integer, dimension(:), intent(in) :: &
         solveVelocity !< Input:

    real(kind=RKIND) :: &
         stress11Edge, &
         stress22Edge, &
         stress12Edge, &
         stress11Vertex, &
         stress22Vertex, &
         stress12Vertex

    integer :: &
         iVertex, &
         iVertexDegree, &
         iEdge, &
         iCellOnEdge, &
         iCell

    integer, pointer :: &
         nVerticesSolve, &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         cellsOnVertex, &
         edgesOnVertex, &
         cellsOnEdge

    real(kind=RKIND), pointer :: &
         sphere_radius

    real(kind=RKIND) :: &
         sphereRadius

    real(kind=RKIND), dimension(:), pointer :: &
         areaTriangle, &
         dcEdge

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nVertices", nVerticesSolve)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "edgesOnVertex", edgesOnVertex)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "areaTriangle", areaTriangle)
    call MPAS_pool_get_array(mesh, "dcEdge", dcEdge)

    ! planar cases with zero sphere radius
    sphereRadius = sphere_radius
    if (sphereRadius == 0.0_RKIND) sphereRadius = 1.0_RKIND

    do iVertex = 1, nVerticesSolve

       stressDivergenceU(iVertex) = 0.0_RKIND
       stressDivergenceV(iVertex) = 0.0_RKIND

       if (solveVelocity(iVertex) == 1) then

          stress11Vertex = 0.0_RKIND
          stress22Vertex = 0.0_RKIND
          stress12Vertex = 0.0_RKIND

          do iVertexDegree = 1, vertexDegree

             ! vertex stresses
             iCell = cellsOnVertex(iVertexDegree,iVertex)

             stress11Vertex = stress11Vertex + stress11(iCell)
             stress22Vertex = stress22Vertex + stress22(iCell)
             stress12Vertex = stress12Vertex + stress12(iCell)

             ! interpolated edge velocity
             iEdge = edgesOnVertex(iVertexDegree,iVertex)

             stress11Edge = 0.0_RKIND
             stress22Edge = 0.0_RKIND
             stress12Edge = 0.0_RKIND

             do iCellOnEdge = 1, 2

                iCell = cellsOnEdge(iCellOnEdge,iEdge)

                stress11Edge = stress11Edge + stress11(iCell)
                stress22Edge = stress22Edge + stress22(iCell)
                stress12Edge = stress12Edge + stress12(iCell)

             enddo ! iCellOnEdge

             stress11Edge = stress11Edge / 2.0_RKIND
             stress22Edge = stress22Edge / 2.0_RKIND
             stress12Edge = stress12Edge / 2.0_RKIND

             stressDivergenceU(iVertex) = stressDivergenceU(iVertex) + &
                  (stress11Edge * normalVectorTriangle(1,iVertexDegree,iVertex) + &
                   stress12Edge * normalVectorTriangle(2,iVertexDegree,iVertex)) * dcEdge(iEdge)

             stressDivergenceV(iVertex) = stressDivergenceV(iVertex) + &
                  (stress22Edge * normalVectorTriangle(2,iVertexDegree,iVertex) + &
                   stress12Edge * normalVectorTriangle(1,iVertexDegree,iVertex)) * dcEdge(iEdge)

          enddo ! iVertexDegree

          stress11Vertex = stress11Vertex / real(vertexDegree, RKIND)
          stress22Vertex = stress22Vertex / real(vertexDegree, RKIND)
          stress12Vertex = stress12Vertex / real(vertexDegree, RKIND)

          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) / areaTriangle(iVertex)
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) / areaTriangle(iVertex)

          ! metric terms
          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) - &
               (tan(latVertexRotated(iVertex)) * stress12Vertex * 2.0_RKIND) / sphereRadius
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) + &
               (tan(latVertexRotated(iVertex)) * (stress11Vertex - stress22Vertex)) / sphereRadius

       endif ! solveVelocity

    enddo ! iVertex

  end subroutine cice_stress_divergence_weak!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  final_divergence_shear_weak
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date July 9th 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_final_divergence_shear_weak(block)

    use cice_velocity_solver_constitutive_relation, only: &
         eccentricitySquared

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocityWeakPool, &
         velocitySolverPool, &
         ridgingPool

    integer, pointer :: &
         nCellsSolve

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         strain11, &
         strain22, &
         strain12, &
         divergence, &
         shear, &
         ridgeConvergence, &
         ridgeShear

    real(kind=RKIND), dimension(:), allocatable :: &
         Delta

    real(kind=RKIND) :: &
         strainDivergence, &
         strainTension, &
         strainShearing

    logical, pointer :: &
         config_use_column_package

    integer :: &
         iCell

    call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
    call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocityWeakPool)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

    call MPAS_pool_get_dimension(meshPool, "nCellsSolve", nCellsSolve)

    call MPAS_pool_get_array(meshPool, "nEdgesOnCell", nEdgesOnCell)

    call MPAS_pool_get_array(velocityWeakPool, "strain11", strain11)
    call MPAS_pool_get_array(velocityWeakPool, "strain22", strain22)
    call MPAS_pool_get_array(velocityWeakPool, "strain12", strain12)

    call MPAS_pool_get_array(velocitySolverPool, "divergence", divergence)
    call MPAS_pool_get_array(velocitySolverPool, "shear", shear)

    allocate(Delta(nCellsSolve))

    do iCell = 1, nCellsSolve

       strainDivergence = strain11(iCell) + strain22(iCell)
       strainTension    = strain11(iCell) - strain22(iCell)
       strainShearing   = strain12(iCell) * 2.0_RKIND

       Delta = sqrt(strainDivergence**2 + (strainTension**2 + strainShearing**2) / eccentricitySquared)

       divergence(iCell) = strainDivergence
       shear(iCell)      = sqrt(strainTension**2 + strainShearing**2)

    enddo ! iCell

    ! ridging parameters
    call MPAS_pool_get_config(block % configs, "config_use_column_package", config_use_column_package)

    if (config_use_column_package) then

       call MPAS_pool_get_subpool(block % structs, "ridging", ridgingPool)

       call MPAS_pool_get_array(ridgingPool, "ridgeConvergence", ridgeConvergence)
       call MPAS_pool_get_array(ridgingPool, "ridgeShear", ridgeShear)

       do iCell = 1, nCellsSolve

          ridgeConvergence(iCell) = -min(divergence(iCell),0.0_RKIND)
          ridgeShear(iCell)       = 0.5_RKIND * (Delta(iCell) - abs(divergence(iCell)))

       enddo ! iCell

    endif ! config_use_column_package

    ! cleanup
    deallocate(Delta)

  end subroutine cice_final_divergence_shear_weak

!-----------------------------------------------------------------------

end module cice_velocity_solver_weak
