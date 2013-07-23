module cice_dyn_hex_weak

  use mpas_grid_types

  implicit none

  private
  public :: init_dynamics_hex_weak, &
            strain_tensor_hex_weak, &
            stress_tensor_hex_weak, &
            stress_divergence_hex_weak, &
            stress_divergence_hex_weak_corners

  public :: spherical_test_strain_weak, &
            spherical_test_divergence_stress_weak, &
            spherical_test_divergence_stress_weak_polygon

contains

  !-------------------------------------------------------------
  ! Initialization
  !-------------------------------------------------------------

  subroutine init_dynamics_hex_weak(mesh, boundary, normal, gridCurvature)

    use cice_dynamics_shared, only: interior_vertices, &
                                    normal_vectors

    type(mesh_type), intent(inout) :: mesh

    type(normal_type),   pointer :: normal
    type(boundary_type), pointer :: boundary

    character(len=*), intent(in) :: gridCurvature

    mesh % xCell % array = mesh % xCell % array / mesh % sphere_radius
    mesh % yCell % array = mesh % yCell % array / mesh % sphere_radius
    mesh % zCell % array = mesh % zCell % array / mesh % sphere_radius

    mesh % xVertex % array = mesh % xVertex % array / mesh % sphere_radius
    mesh % yVertex % array = mesh % yVertex % array / mesh % sphere_radius
    mesh % zVertex % array = mesh % zVertex % array / mesh % sphere_radius

    mesh % xEdge % array = mesh % xEdge % array / mesh % sphere_radius
    mesh % yEdge % array = mesh % yEdge % array / mesh % sphere_radius
    mesh % zEdge % array = mesh % zEdge % array / mesh % sphere_radius

    mesh % dcEdge % array = mesh % dcEdge % array / mesh % sphere_radius
    mesh % dvEdge % array = mesh % dvEdge % array / mesh % sphere_radius

    mesh % areaCell % array = mesh % areaCell % array / (mesh % sphere_radius)**2
    mesh % areaTriangle % array = mesh % areaTriangle % array / (mesh % sphere_radius)**2

    mesh % sphere_radius = 1.0_RKIND

    call interior_vertices(mesh, &
         boundary % interiorVertex % array)

    call normal_vectors(mesh, &
         normal % normalVectorPolygon % array,         &
         normal % normalVectorTriangle % array,        &
         normal % normalVectorTriangleCorners % array, &
         normal % latCellRotated % array,              &
         normal % latVertexRotated % array,            &
         boundary % interiorVertex % array,            &
         gridCurvature)

  end subroutine init_dynamics_hex_weak

  !-------------------------------------------------------------
  ! Time step
  !-------------------------------------------------------------

  subroutine strain_tensor_hex_weak(mesh, &
                                    strain11, strain22,   &
                                    strain12,             &
                                    uVelocity, vVelocity, &
                                    normalVectorPolygon,  &
                                    latCellRotated,       &
                                    solveStress)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, &
         vVelocity, &
         latCellRotated

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorPolygon

    integer, dimension(:), intent(in) :: &
         solveStress

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

    real(kind=RKIND), parameter :: &
         rad_to_deg = 180.0_RKIND / 3.14159265359_RKIND
     
    do iCell = 1, mesh % nCells

       strain11(iCell) = 0.0_RKIND
       strain22(iCell) = 0.0_RKIND
       strain12(iCell) = 0.0_RKIND

       if (solveStress(iCell) == 1) then

          uCellCentre = 0.0_RKIND
          vCellCentre = 0.0_RKIND

          do iEdgeOnCell = 1, mesh % nEdgesOnCell % array(iCell)

             ! cell centre velocities
             iVertex = mesh % verticesOnCell % array(iEdgeOnCell,iCell)

             uCellCentre = uCellCentre + uVelocity(iVertex)
             vCellCentre = vCellCentre + vVelocity(iVertex)

             ! interpolated edge velocity
             iEdge = mesh % edgesOnCell % array(iEdgeOnCell,iCell)

             uVelocityEdge = 0.0_RKIND
             vVelocityEdge = 0.0_RKIND

             do iVertexOnEdge = 1, 2

                iVertex = mesh % verticesOnEdge % array(iVertexOnEdge,iEdge)

                uVelocityEdge = uVelocityEdge + uVelocity(iVertex)
                vVelocityEdge = vVelocityEdge + vVelocity(iVertex)

             enddo ! iVertexOnEdge

             uVelocityEdge = uVelocityEdge / 2.0_RKIND
             vVelocityEdge = vVelocityEdge / 2.0_RKIND

             ! summation over edges
             strain11(iCell) = strain11(iCell) + uVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) * mesh % dvEdge % array(iEdge)
             strain22(iCell) = strain22(iCell) + vVelocityEdge * normalVectorPolygon(2,iEdgeOnCell,iCell) * mesh % dvEdge % array(iEdge)
             strain12(iCell) = strain12(iCell) + 0.5_RKIND * ( &
                  uVelocityEdge * normalVectorPolygon(2,iEdgeOnCell,iCell) + &
                  vVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) ) * mesh % dvEdge % array(iEdge)

          enddo ! iEdgeOnCell

          uCellCentre = uCellCentre / real(mesh % nEdgesOnCell % array(iCell), RKIND)
          vCellCentre = vCellCentre / real(mesh % nEdgesOnCell % array(iCell), RKIND)

          strain11(iCell) = strain11(iCell) / mesh % areaCell % array(iCell)
          strain22(iCell) = strain22(iCell) / mesh % areaCell % array(iCell)
          strain12(iCell) = strain12(iCell) / mesh % areaCell % array(iCell)

          ! metric terms
          strain11(iCell) = strain11(iCell) - (vCellCentre * tan(latCellRotated(iCell))) / mesh % sphere_radius
          strain12(iCell) = strain12(iCell) + (uCellCentre * tan(latCellRotated(iCell)) * 0.5_RKIND) / mesh % sphere_radius

       endif ! solveStress

    enddo ! iCell

  end subroutine strain_tensor_hex_weak

  !-------------------------------------------------------------

  subroutine stress_tensor_hex_weak(mesh, &
                                    stress11, stress22, stress12,     &
                                    strain11, strain22, strain12,     &
                                    icePressure, replacementPressure, &
                                    solveStress,                      &
                                    dtElastic,   evpDamping)

    use cice_dynamics_shared, only: evp_constitutive_relation, &
                                    linear_constitutive_relation

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(inout) :: &
         stress11, &
         stress22, &
         stress12, &
         replacementPressure

    real(kind=RKIND), dimension(:), intent(in) :: &
         strain11, &
         strain22, &
         strain12, &
         icePressure

    integer, dimension(:), intent(in) :: &
         solveStress

    real(kind=RKIND), intent(in) :: &
         dtElastic

    logical, intent(in) :: &
         evpDamping

    integer :: iCell

    do iCell = 1, mesh % nCells

       if (solveStress(iCell) == 1) then

          call evp_constitutive_relation(stress11(iCell),    stress22(iCell),    stress12(iCell), &
                                         strain11(iCell),    strain22(iCell),    strain12(iCell), &
                                         icePressure(iCell), replacementPressure(iCell),          &
                                         mesh % areaCell % array(iCell),                          &
                                         dtElastic,          evpDamping)

          !call linear_constitutive_relation(stress11(iCell), stress22(iCell), stress12(iCell), &
          !                                  strain11(iCell), strain22(iCell), strain12(iCell))

       else

          stress11(iCell) = 0.0_RKIND
          stress22(iCell) = 0.0_RKIND
          stress12(iCell) = 0.0_RKIND

       endif ! solveStress

    end do ! iCell

  end subroutine stress_tensor_hex_weak
  
  !-------------------------------------------------------------

  subroutine stress_divergence_hex_weak(mesh, &
                                        stressDivergenceU,    stressDivergenceV, &
                                        stress11,             stress22,          &
                                        stress12,                                &
                                        normalVectorTriangle, latVertexRotated,  &
                                        solveVelocity)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, &
         stressDivergenceV

    real(kind=RKIND), dimension(:), intent(in) :: &
         stress11, &
         stress22, &
         stress12, &
         latVertexRotated

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorTriangle

    integer, dimension(:), intent(in) :: &
         solveVelocity

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

    do iVertex = 1, mesh % nVertices

       stressDivergenceU(iVertex) = 0.0_RKIND
       stressDivergenceV(iVertex) = 0.0_RKIND

       if (solveVelocity(iVertex) == 1) then

          stress11Vertex = 0.0_RKIND
          stress22Vertex = 0.0_RKIND
          stress12Vertex = 0.0_RKIND
          
          do iVertexDegree = 1, mesh % vertexDegree
             
             ! vertex stresses
             iCell = mesh % cellsOnVertex % array(iVertexDegree,iVertex)

             stress11Vertex = stress11Vertex + stress11(iCell)
             stress22Vertex = stress22Vertex + stress22(iCell)
             stress12Vertex = stress12Vertex + stress12(iCell)
             
             ! interpolated edge velocity
             iEdge = mesh % edgesOnVertex % array(iVertexDegree,iVertex)
             
             stress11Edge = 0.0_RKIND
             stress22Edge = 0.0_RKIND
             stress12Edge = 0.0_RKIND
             
             do iCellOnEdge = 1, 2
                
                iCell = mesh % cellsOnEdge % array(iCellOnEdge,iEdge)
                
                stress11Edge = stress11Edge + stress11(iCell)
                stress22Edge = stress22Edge + stress22(iCell)
                stress12Edge = stress12Edge + stress12(iCell)
                
             enddo ! iCellOnEdge
             
             stress11Edge = stress11Edge / 2.0_RKIND
             stress22Edge = stress22Edge / 2.0_RKIND
             stress12Edge = stress12Edge / 2.0_RKIND
             
             stressDivergenceU(iVertex) = stressDivergenceU(iVertex) + &
                  (stress11Edge * normalVectorTriangle(1,iVertexDegree,iVertex) + &
                   stress12Edge * normalVectorTriangle(2,iVertexDegree,iVertex)) * mesh % dcEdge % array(iEdge)
             
             stressDivergenceV(iVertex) = stressDivergenceV(iVertex) + &
                  (stress22Edge * normalVectorTriangle(2,iVertexDegree,iVertex) + &
                   stress12Edge * normalVectorTriangle(1,iVertexDegree,iVertex)) * mesh % dcEdge % array(iEdge)

          enddo ! iVertexDegree

          stress11Vertex = stress11Vertex / real(mesh % vertexDegree, RKIND)
          stress22Vertex = stress22Vertex / real(mesh % vertexDegree, RKIND)
          stress12Vertex = stress12Vertex / real(mesh % vertexDegree, RKIND)
          
          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) / mesh % areaTriangle % array(iVertex)
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) / mesh % areaTriangle % array(iVertex)

          ! metric terms
          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) - &
               (tan(latVertexRotated(iVertex)) * stress12Vertex * 2.0_RKIND) / mesh % sphere_radius
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) + &
               (tan(latVertexRotated(iVertex)) * (stress11Vertex - stress22Vertex)) / mesh % sphere_radius

       endif ! solveVelocity

    enddo ! iVertex

  end subroutine stress_divergence_hex_weak

  !-------------------------------------------------------------

  subroutine stress_divergence_hex_weak_corners(mesh, &
                                                stressDivergenceU,           stressDivergenceV, &
                                                stress11,                    stress22,          &
                                                stress12,                                       &
                                                normalVectorTriangleCorners, latVertexRotated,  &
                                                solveVelocity)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, &
         stressDivergenceV

    real(kind=RKIND), dimension(:), intent(in) :: &
         stress11, &
         stress22, &
         stress12, &
         latVertexRotated

    real(kind=RKIND), dimension(:,:,:,:), intent(in) :: &
         normalVectorTriangleCorners

    integer, dimension(:), intent(in) :: &
         solveVelocity

    real(kind=RKIND) :: &
         stress11Vertex, &
         stress22Vertex, &
         stress12Vertex

    integer :: &
         iVertex, &
         iVertexDegree, &
         iEdge, &
         iCellOnEdge, &
         iCell

    do iVertex = 1, mesh % nVertices

       stressDivergenceU(iVertex) = 0.0_RKIND
       stressDivergenceV(iVertex) = 0.0_RKIND

       if (solveVelocity(iVertex) == 1) then

          stress11Vertex = 0.0_RKIND
          stress22Vertex = 0.0_RKIND
          stress12Vertex = 0.0_RKIND
          
          do iVertexDegree = 1, mesh % vertexDegree
             
             ! vertex stresses
             iCell = mesh % cellsOnVertex % array(iVertexDegree,iVertex)

             stress11Vertex = stress11Vertex + stress11(iCell)
             stress22Vertex = stress22Vertex + stress22(iCell)
             stress12Vertex = stress12Vertex + stress12(iCell)
             
             ! interpolated edge velocity
             iEdge = mesh % edgesOnVertex % array(iVertexDegree,iVertex)
             
             do iCellOnEdge = 1, 2
                
                iCell = mesh % cellsOnEdge % array(iCellOnEdge,iEdge)
                
                stressDivergenceU(iVertex) = stressDivergenceU(iVertex) + &
                     (stress11(iCell) * normalVectorTriangleCorners(1,iCellOnEdge,iVertexDegree,iVertex) + &
                      stress12(iCell) * normalVectorTriangleCorners(2,iCellOnEdge,iVertexDegree,iVertex)) * mesh % dcEdge % array(iEdge)
             
                stressDivergenceV(iVertex) = stressDivergenceV(iVertex) + &
                     (stress22(iCell) * normalVectorTriangleCorners(2,iCellOnEdge,iVertexDegree,iVertex) + &
                      stress12(iCell) * normalVectorTriangleCorners(1,iCellOnEdge,iVertexDegree,iVertex)) * mesh % dcEdge % array(iEdge)

             enddo ! iCellOnEdge

          enddo ! iVertexDegree

          stress11Vertex = stress11Vertex / real(mesh % vertexDegree, RKIND)
          stress22Vertex = stress22Vertex / real(mesh % vertexDegree, RKIND)
          stress12Vertex = stress12Vertex / real(mesh % vertexDegree, RKIND)
          
          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) / (mesh % areaTriangle % array(iVertex) * 2.0_RKIND)
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) / (mesh % areaTriangle % array(iVertex) * 2.0_RKIND)

          ! metric terms
          !stressDivergenceU(iVertex) = stressDivergenceU(iVertex) - &
          !     (tan(latVertexRotated(iVertex)) * stress12Vertex * 2.0_RKIND) / mesh % sphere_radius
          !stressDivergenceV(iVertex) = stressDivergenceV(iVertex) + &
          !     (tan(latVertexRotated(iVertex)) * (stress11Vertex - stress22Vertex)) / mesh % sphere_radius

       endif ! solveVelocity

    enddo ! iVertex

  end subroutine stress_divergence_hex_weak_corners

  !-------------------------------------------------------------

  subroutine stress_divergence_hex_weak_polygon(mesh, &
                                        stressDivergenceU,   stressDivergenceV, &
                                        stress11,            stress22,          &
                                        stress12,                               &
                                        normalVectorPolygon, latCellRotated,    &
                                        solveVelocity)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, &
         stressDivergenceV

    real(kind=RKIND), dimension(:), intent(in) :: &
         stress11, &
         stress22, &
         stress12, &
         latCellRotated

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorPolygon

    integer, dimension(:), intent(in) :: &
         solveVelocity

    real(kind=RKIND) :: &
         stress11Edge, &
         stress22Edge, &
         stress12Edge, &
         stress11Cell, &
         stress22Cell, &
         stress12Cell

    integer :: &
         iVertex, &
         iEdgeOnCell, &
         iEdge, &
         iVertexOnEdge, &
         iCell

    do iCell = 1, mesh % nCells

       stressDivergenceU(iCell) = 0.0_RKIND
       stressDivergenceV(iCell) = 0.0_RKIND

       if (solveVelocity(iCell) == 1) then

          stress11Cell = 0.0_RKIND
          stress22Cell = 0.0_RKIND
          stress12Cell = 0.0_RKIND
          
          do iEdgeOnCell = 1, mesh % nEdgesOnCell % array(iCell)
             
             ! vertex stresses
             iVertex = mesh % verticesOnCell % array(iEdgeOnCell,iCell)

             stress11Cell = stress11Cell + stress11(iVertex)
             stress22Cell = stress22Cell + stress22(iVertex)
             stress12Cell = stress12Cell + stress12(iVertex)
             
             ! interpolated edge velocity
             iEdge = mesh % edgesOnCell % array(iEdgeOnCell,iCell)
             
             stress11Edge = 0.0_RKIND
             stress22Edge = 0.0_RKIND
             stress12Edge = 0.0_RKIND
             
             do iVertexOnEdge = 1, 2
                
                iVertex = mesh % verticesOnEdge % array(iVertexOnEdge,iEdge)
                
                stress11Edge = stress11Edge + stress11(iVertex)
                stress22Edge = stress22Edge + stress22(iVertex)
                stress12Edge = stress12Edge + stress12(iVertex)
                
             enddo ! iVertexOnEdge
             
             stress11Edge = stress11Edge / 2.0_RKIND
             stress22Edge = stress22Edge / 2.0_RKIND
             stress12Edge = stress12Edge / 2.0_RKIND
             
             stressDivergenceU(iCell) = stressDivergenceU(iCell) + &
                  (stress11Edge * normalVectorPolygon(1,iEdgeOnCell,iCell) + &
                   stress12Edge * normalVectorPolygon(2,iEdgeOnCell,iCell)) * mesh % dvEdge % array(iEdge)
             
             stressDivergenceV(iCell) = stressDivergenceV(iCell) + &
                  (stress22Edge * normalVectorPolygon(2,iEdgeOnCell,iCell) + &
                   stress12Edge * normalVectorPolygon(1,iEdgeOnCell,iCell)) * mesh % dvEdge % array(iEdge)

             !write(*,*) iCell, iEdgeOnCell, stress11Edge, normalVectorPolygon(1,iEdgeOnCell,iCell), &
             !      stress12Edge, normalVectorPolygon(2,iEdgeOnCell,iCell), mesh % dvEdge % array(iEdge)

          enddo ! iEdgeOnCell

          stress11Cell = stress11Cell / real(mesh % nEdgesOnCell % array(iCell), RKIND)
          stress22Cell = stress22Cell / real(mesh % nEdgesOnCell % array(iCell), RKIND)
          stress12Cell = stress12Cell / real(mesh % nEdgesOnCell % array(iCell), RKIND)
          
          !write(*,*) iCell, stressDivergenceU(iCell), stressDivergenceU(iCell) / mesh % areaCell % array(iCell), mesh % areaCell % array(iCell)

          stressDivergenceU(iCell) = stressDivergenceU(iCell) / mesh % areaCell % array(iCell)
          stressDivergenceV(iCell) = stressDivergenceV(iCell) / mesh % areaCell % array(iCell)

          ! metric terms
          !stressDivergenceU(iCell) = stressDivergenceU(iCell) - &
          !     (tan(latCellRotated(iCell)) * stress12Cell * 2.0_RKIND) / mesh % sphere_radius
          !stressDivergenceV(iCell) = stressDivergenceV(iCell) + &
          !     (tan(latCellRotated(iCell)) * (stress11Cell - stress22Cell)) / mesh % sphere_radius

       endif ! solveVelocity

    enddo ! iCell
!stop
  end subroutine stress_divergence_hex_weak_polygon

  !-------------------------------------------------------------
  ! test routines
  !-------------------------------------------------------------

  subroutine spherical_test_strain_weak(mesh, hexdyn, hexweak, normal)

    use cice_testing, only: spherical_test_strain, plot_pop_grid, rms_difference

    use cice_dynamics_shared, only: latlon_from_xyz, grid_rotation_forward, grid_rotation_forward_nochange
    
    type(mesh_type), intent(in) :: mesh
    
    type(normal_type),   pointer :: normal
    type(hexdyn_type),   pointer :: hexdyn
    type(hexweak_type),  pointer :: hexweak
    
    real(kind=RKIND), dimension(:), allocatable :: &
         strain11_test, &
         strain22_test, &
         strain12_test, &
         strain11_diff, &
         strain22_diff, &
         strain12_diff, &
         latVertexRotated, &
         lonVertexRotated, &
         latCellRotated, &
         lonCellRotated

    integer :: &
         iCell, &
         iCellOnCell, &
         iVertex

    character(len=200), parameter :: &
         test_type = "sinusoidal2"

    real(kind=RKIND) :: &
         up, vp, ut, vt, &
         xp, yp, zp, &
         rms_strain11, &
         rms_strain22, &
         rms_strain12

    real(kind=RKIND), parameter :: &
         pi = 3.14159265359_RKIND, &
         deg_to_rad = pi / 180.0_RKIND, &
         rad_to_deg = 180.0_RKIND / pi

    hexdyn % solveStress % array = 1

    allocate(strain11_test(mesh % nCells))
    allocate(strain22_test(mesh % nCells))
    allocate(strain12_test(mesh % nCells))
    allocate(strain11_diff(mesh % nCells))
    allocate(strain22_diff(mesh % nCells))
    allocate(strain12_diff(mesh % nCells))

    allocate(lonVertexRotated(mesh % nVertices))
    allocate(latVertexRotated(mesh % nVertices))
    allocate(lonCellRotated(mesh % nCells))
    allocate(latCellRotated(mesh % nCells))

    do iVertex = 1, mesh % nVertices
    
       call grid_rotation_forward(xp, yp, zp, &
            mesh % xVertex % array(iVertex), &
            mesh % yVertex % array(iVertex), &
            mesh % zVertex % array(iVertex))

       call latlon_from_xyz(latVertexRotated(iVertex), lonVertexRotated(iVertex), &
            xp, yp, zp, mesh % sphere_radius)

    enddo ! iVertex

    do iCell = 1, mesh % nCells

       call grid_rotation_forward(xp, yp, zp, &
            mesh % xCell % array(iCell), &
            mesh % yCell % array(iCell), &
            mesh % zCell % array(iCell))

       call latlon_from_xyz(latCellRotated(iCell), lonCellRotated(iCell), &
            xp, yp, zp, mesh % sphere_radius)

    enddo ! iCell
    
    ! set test velocities and expected outputs
    call spherical_test_strain(mesh, &
         hexdyn % uVelocity % array, &
         hexdyn % vVelocity % array, &
         lonVertexRotated,           &
         latVertexRotated,           &
         strain11_test,              &
         strain22_test,              &
         strain12_test,              &
         lonCellRotated,             &
         latCellRotated,             &
         test_type)

    ! calculate strain
    call strain_tensor_hex_weak(mesh, &
         hexweak % strain11 % array,           &
         hexweak % strain22 % array,           &
         hexweak % strain12 % array,           &
         hexdyn % uVelocity % array,           &
         hexdyn % vVelocity % array,           &
         normal % normalVectorPolygon % array, &
         normal % latCellRotated % array,      &
         hexdyn % solveStress % array)

    ! calculate difference
    strain11_diff = hexweak % strain11 % array(1:mesh % nCells) - &
                    strain11_test
    strain22_diff = hexweak % strain22 % array(1:mesh % nCells) - &
                    strain22_test
    strain12_diff = hexweak % strain12 % array(1:mesh % nCells) - &
                    strain12_test

    do iCell = 1, mesh % nCells

       write(*,*) iCell, hexdyn % solveStress % array(iCell), &
            hexweak % strain11 % array(iCell), strain11_test(iCell), &
            hexweak % strain22 % array(iCell), strain22_test(iCell), &
            hexweak % strain12 % array(iCell), strain12_test(iCell), &
            hexweak % strain12 % array(iCell) - strain12_test(iCell)

    enddo ! iCell

   ! set cells to solve
    do iCell = 1, mesh % nCells
       hexdyn % solveStress % array(iCell) = 1
       !do iCellOnCell = 1, mesh % nEdgesOnCell % array(iCell)
       !   ! cut off boundary cells
       !   if (mesh % cellsOnCell % array(iCellOnCell,iCell) > mesh % nCells) then
       !      hexdyn % solveStress % array(iCell) = 0
       !      exit
       !   endif
       !enddo ! iCellOnCell
       ! cut off low latitude cells
       if (mesh % latCell % array(iCell) < 30.0 * deg_to_rad .and. mesh % latCell % array(iCell) > -30.0 * deg_to_rad) then
          hexdyn % solveStress % array(iCell) = 0
       endif
    enddo ! iCell

    open(55,file='strain11.txt')
    open(56,file='strain22.txt')
    open(57,file='strain12.txt')
    do iCell = 1, mesh % nCells
       write(55,*) iCell, mesh % lonCell % array(iCell) * rad_to_deg, mesh % latCell % array(iCell) * rad_to_deg, &
            hexweak % strain11 % array(iCell), strain11_test(iCell), strain11_diff(iCell), &
            strain11_diff(iCell) * real(hexdyn % solveStress % array(iCell),RKIND)
       write(56,*) iCell, mesh % lonCell % array(iCell) * rad_to_deg, mesh % latCell % array(iCell) * rad_to_deg, &
            hexweak % strain22 % array(iCell), strain22_test(iCell), strain22_diff(iCell), &
            strain22_diff(iCell) * real(hexdyn % solveStress % array(iCell),RKIND)
       write(57,*) iCell, mesh % lonCell % array(iCell) * rad_to_deg, mesh % latCell % array(iCell) * rad_to_deg, &
            hexweak % strain12 % array(iCell), strain12_test(iCell), strain12_diff(iCell), &
            strain12_diff(iCell) * real(hexdyn % solveStress % array(iCell),RKIND)
    enddo ! iCell
    close(55)
    close(56)
    close(57)

    ! rms difference
    call rms_difference( &
         hexweak % strain11 % array(1:mesh % nCells),   &
         strain11_test,                                 &
         hexdyn % solveStress % array(1:mesh % nCells), &
         mesh % nCells,                                 &
         rms_strain11)

    call rms_difference( &
         hexweak % strain22 % array(1:mesh % nCells),   &
         strain22_test,                                 &
         hexdyn % solveStress % array(1:mesh % nCells), &
         mesh % nCells,                                 &
         rms_strain22)

    call rms_difference( &
         hexweak % strain12 % array(1:mesh % nCells),   &
         strain12_test,                                 &
         hexdyn % solveStress % array(1:mesh % nCells), &
         mesh % nCells,                                 &
         rms_strain12)

    write(*,*) "RMS:", rms_strain11, rms_strain22, rms_strain12

    deallocate(strain11_test)
    deallocate(strain22_test)
    deallocate(strain12_test)
    deallocate(strain11_diff)
    deallocate(strain22_diff)
    deallocate(strain12_diff)

    deallocate(lonVertexRotated)
    deallocate(latVertexRotated)
    deallocate(lonCellRotated)
    deallocate(latCellRotated)

  end subroutine spherical_test_strain_weak

  !-------------------------------------------------------------

  subroutine spherical_test_divergence_stress_weak(mesh, hexdyn, hexweak, normal)

    use cice_testing, only: spherical_test_divergence_stress, rms_difference

    use cice_dynamics_shared, only: grid_rotation_forward, latlon_from_xyz

    type(mesh_type), intent(in) :: mesh

    type(normal_type),   pointer :: normal
    type(hexdyn_type),   pointer :: hexdyn
    type(hexweak_type),  pointer :: hexweak

    real(kind=RKIND), dimension(:), allocatable :: &
         stressDivergenceU_test, &
         stressDivergenceV_test, &
         stressDivergenceU_diff, &
         stressDivergenceV_diff, &
         latVertexRotated, &
         lonVertexRotated, &
         latCellRotated, &
         lonCellRotated
    
    real(kind=RKIND) :: &
         xp, yp, zp, &
         rms_divergenceU, &
         rms_divergenceV

    integer :: &
         iCell, &
         iVertex

    real(kind=RKIND), parameter :: &
         pi = 3.14159265359_RKIND, &
         deg_to_rad = pi / 180.0_RKIND, &
         rad_to_deg = 180.0_RKIND / pi

    character(len=200), parameter :: &
         test_type = "test4"

    hexdyn % solveVelocity % array = 1

    allocate(stressDivergenceU_test(mesh % nVertices))
    allocate(stressDivergenceV_test(mesh % nVertices))
    allocate(stressDivergenceU_diff(mesh % nVertices))
    allocate(stressDivergenceV_diff(mesh % nVertices))

    allocate(lonVertexRotated(mesh % nVertices))
    allocate(latVertexRotated(mesh % nVertices))
    allocate(lonCellRotated(mesh % nCells))
    allocate(latCellRotated(mesh % nCells))

    do iVertex = 1, mesh % nVertices
    
       call grid_rotation_forward(xp, yp, zp, &
            mesh % xVertex % array(iVertex), &
            mesh % yVertex % array(iVertex), &
            mesh % zVertex % array(iVertex))

       call latlon_from_xyz(latVertexRotated(iVertex), lonVertexRotated(iVertex), &
            xp, yp, zp, mesh % sphere_radius)

    enddo ! iVertex

    do iCell = 1, mesh % nCells

       call grid_rotation_forward(xp, yp, zp, &
            mesh % xCell % array(iCell), &
            mesh % yCell % array(iCell), &
            mesh % zCell % array(iCell))

       call latlon_from_xyz(latCellRotated(iCell), lonCellRotated(iCell), &
            xp, yp, zp, mesh % sphere_radius)

    enddo ! iCell

    ! set test velocities and expected outputs
    call spherical_test_divergence_stress(mesh, &
         hexweak % stress11 % array, &
         hexweak % stress22 % array, &
         hexweak % stress12 % array, &
         lonCellRotated,             &
         latCellRotated,             &
         stressDivergenceU_test,     &
         stressDivergenceV_test,     &
         lonVertexRotated,           &
         latVertexRotated,           &
         test_type)

    ! calculate divergence of stress
    call stress_divergence_hex_weak(mesh, &
         hexdyn % stressDivergenceU % array,    &
         hexdyn % stressDivergenceV % array,    &
         hexweak % stress11 % array,            &
         hexweak % stress22 % array,            &
         hexweak % stress12 % array,            &
         normal % normalVectorTriangle % array, &
         normal % latVertexRotated % array,     &
         hexdyn % solveVelocity % array)

    ! calculate difference
    stressDivergenceU_diff = hexdyn % stressDivergenceU % array(1:mesh % nVertices) - &
                             stressDivergenceU_test
    stressDivergenceV_diff = hexdyn % stressDivergenceV % array(1:mesh % nVertices) - &
                             stressDivergenceV_test

    do iVertex = 1, mesh % nVertices

       write(*,*) iVertex, hexdyn % solveVelocity % array(iVertex), &
            hexdyn % stressDivergenceU % array(iVertex), stressDivergenceU_test(iVertex), &
            hexdyn % stressDivergenceV % array(iVertex), stressDivergenceV_test(iVertex), &
            hexdyn % stressDivergenceU % array(iVertex) - stressDivergenceU_test(iVertex), &
            hexdyn % stressDivergenceV % array(iVertex) - stressDivergenceV_test(iVertex)

    enddo ! iVertex

    ! set vertices to solve
    do iVertex = 1, mesh % nVertices
       hexdyn % solveVelocity % array(iVertex) = 1
       ! cut off low latitude cells
       if (mesh % latVertex % array(iVertex) < 30.0 * deg_to_rad .and. mesh % latVertex % array(iVertex) > -30.0 * deg_to_rad) then
          hexdyn % solveVelocity % array(iVertex) = 0
       endif
    enddo ! iVertex

    open(55,file='divergenceu.txt')
    open(56,file='divergencev.txt')
    do iVertex = 1, mesh % nVertices
       write(55,*) iVertex, mesh % lonVertex % array(iVertex) * rad_to_deg, mesh % latVertex % array(iVertex) * rad_to_deg, &
            hexdyn % stressDivergenceU % array(iVertex), stressDivergenceU_test(iVertex), stressDivergenceU_diff(iVertex), &
            stressDivergenceU_diff(iVertex) * real(hexdyn % solveVelocity % array(iVertex),RKIND)
       write(56,*) iVertex, mesh % lonVertex % array(iVertex) * rad_to_deg, mesh % latVertex % array(iVertex) * rad_to_deg, &
            hexdyn % stressDivergenceV % array(iVertex), stressDivergenceV_test(iVertex), stressDivergenceV_diff(iVertex), &
            stressDivergenceV_diff(iVertex) * real(hexdyn % solveVelocity % array(iVertex),RKIND)
    enddo ! iVertex
    close(55)
    close(56)

    ! rms difference
    call rms_difference( &
         hexdyn % stressDivergenceU % array(1:mesh % nVertices), &
         stressDivergenceU_test,                                  &
         hexdyn % solveVelocity % array(1:mesh % nVertices),      &
         mesh % nVertices,                                        &
         rms_divergenceU)

    call rms_difference( &
         hexdyn % stressDivergenceV % array(1:mesh % nVertices), &
         stressDivergenceV_test,                                  &
         hexdyn % solveVelocity % array(1:mesh % nVertices),      &
         mesh % nVertices,                                        &
         rms_divergenceV)

    write(*,*) "RMS:", rms_divergenceU, rms_divergenceV

    ! cleanup
    deallocate(stressDivergenceU_test)
    deallocate(stressDivergenceV_test)
    deallocate(stressDivergenceU_diff)
    deallocate(stressDivergenceV_diff)

    deallocate(lonVertexRotated)
    deallocate(latVertexRotated)
    deallocate(lonCellRotated)
    deallocate(latCellRotated)

  end subroutine spherical_test_divergence_stress_weak

  !-------------------------------------------------------------

  subroutine spherical_test_divergence_stress_weak_polygon(mesh, hexdyn, hexweak, normal)

    use cice_testing, only: spherical_test_divergence_stress, rms_difference

    use cice_dynamics_shared, only: grid_rotation_forward, latlon_from_xyz

    type(mesh_type), intent(in) :: mesh

    type(normal_type),   pointer :: normal
    type(hexdyn_type),   pointer :: hexdyn
    type(hexweak_type),  pointer :: hexweak

    real(kind=RKIND), dimension(:), allocatable :: &
         stressDivergenceU,      &
         stressDivergenceV,      &
         stressDivergenceU_test, &
         stressDivergenceV_test, &
         stressDivergenceU_diff, &
         stressDivergenceV_diff, &
         stress11,               &
         stress22,               &
         stress12,               &
         latVertexRotated,       &
         lonVertexRotated,       &
         latCellRotated,         &
         lonCellRotated

    integer, dimension(:), allocatable :: &
         solveVelocity
    
    real(kind=RKIND) :: &
         xp, yp, zp, &
         rms_divergenceU, &
         rms_divergenceV

    integer :: &
         iCell, &
         iVertex

    real(kind=RKIND), parameter :: &
         pi = 3.14159265359_RKIND, &
         deg_to_rad = pi / 180.0_RKIND, &
         rad_to_deg = 180.0_RKIND / pi

    character(len=200), parameter :: &
         test_type = "test1"


    allocate(solveVelocity(mesh % nCells))

    allocate(stressDivergenceU(mesh % nCells))
    allocate(stressDivergenceV(mesh % nCells))
    allocate(stressDivergenceU_test(mesh % nCells))
    allocate(stressDivergenceV_test(mesh % nCells))
    allocate(stressDivergenceU_diff(mesh % nCells))
    allocate(stressDivergenceV_diff(mesh % nCells))

    allocate(stress11(mesh % nVertices))
    allocate(stress22(mesh % nVertices))
    allocate(stress12(mesh % nVertices))

    allocate(lonVertexRotated(mesh % nVertices))
    allocate(latVertexRotated(mesh % nVertices))
    allocate(lonCellRotated(mesh % nCells))
    allocate(latCellRotated(mesh % nCells))


    solveVelocity = 1

    do iVertex = 1, mesh % nVertices
    
       call grid_rotation_forward(xp, yp, zp, &
            mesh % xVertex % array(iVertex), &
            mesh % yVertex % array(iVertex), &
            mesh % zVertex % array(iVertex))

       call latlon_from_xyz(latVertexRotated(iVertex), lonVertexRotated(iVertex), &
            xp, yp, zp, mesh % sphere_radius)

    enddo ! iVertex

    do iCell = 1, mesh % nCells

       call grid_rotation_forward(xp, yp, zp, &
            mesh % xCell % array(iCell), &
            mesh % yCell % array(iCell), &
            mesh % zCell % array(iCell))

       call latlon_from_xyz(latCellRotated(iCell), lonCellRotated(iCell), &
            xp, yp, zp, mesh % sphere_radius)

    enddo ! iCell

    ! set test velocities and expected outputs
    call spherical_test_divergence_stress(mesh, &
         stress11,               &
         stress22,               &
         stress12,               &
         lonVertexRotated,       &
         latVertexRotated,       &
         stressDivergenceU_test, &
         stressDivergenceV_test, &
         lonCellRotated,         &
         latCellRotated,         &
         test_type)

    ! calculate divergence of stress
    call stress_divergence_hex_weak_polygon(mesh, &
         stressDivergenceU,                    &
         stressDivergenceV,                    &
         stress11,                             &
         stress22,                             &
         stress12,                             &
         normal % normalVectorPolygon % array, &
         normal % latCellRotated % array,      &
         solveVelocity)

    ! calculate difference
    stressDivergenceU_diff = stressDivergenceU - &
                             stressDivergenceU_test
    stressDivergenceV_diff = stressDivergenceV - &
                             stressDivergenceV_test

    do iCell = 1, mesh % nCells

       write(*,*) iVertex, solveVelocity(iCell), &
            stressDivergenceU(iCell), stressDivergenceU_test(iCell), &
            stressDivergenceV(iCell), stressDivergenceV_test(iCell), &
            stressDivergenceU(iCell) - stressDivergenceU_test(iCell), &
            stressDivergenceV(iCell) - stressDivergenceV_test(iCell)

    enddo ! iCell

    ! set vertices to solve
    do iCell = 1, mesh % nCells
       solveVelocity(iCell) = 1
       ! cut off low latitude cells
       if (mesh % latCell % array(iCell) < 30.0 * deg_to_rad .and. mesh % latCell % array(iCell) > -30.0 * deg_to_rad) then
          solveVelocity(iCell) = 0
       endif
    enddo ! iVertex

    !open(55,file='divergenceu.txt')
    !open(56,file='divergencev.txt')
    !do iVertex = 1, mesh % nVertices
    !   write(55,*) iVertex, mesh % lonVertex % array(iVertex) * rad_to_deg, mesh % latVertex % array(iVertex) * rad_to_deg, &
    !        hexdyn % stressDivergenceU % array(iVertex), stressDivergenceU_test(iVertex), stressDivergenceU_diff(iVertex), &
    !        stressDivergenceU_diff(iVertex) * real(hexdyn % solveVelocity % array(iVertex),RKIND)
    !   write(56,*) iVertex, mesh % lonVertex % array(iVertex) * rad_to_deg, mesh % latVertex % array(iVertex) * rad_to_deg, &
    !        hexdyn % stressDivergenceV % array(iVertex), stressDivergenceV_test(iVertex), stressDivergenceV_diff(iVertex), &
    !        stressDivergenceV_diff(iVertex) * real(hexdyn % solveVelocity % array(iVertex),RKIND)
    !enddo ! iVertex
    !close(55)
    !close(56)

    ! rms difference
    call rms_difference( &
         stressDivergenceU,      &
         stressDivergenceU_test, &
         solveVelocity,          &
         mesh % nCells,          &
         rms_divergenceU)

    call rms_difference( &
         stressDivergenceV,      &
         stressDivergenceV_test, &
         solveVelocity,          &
         mesh % nCells,          &
         rms_divergenceV)

    write(*,*) "RMS:", rms_divergenceU, rms_divergenceV

    ! cleanup
    deallocate(stressDivergenceU)
    deallocate(stressDivergenceV)
    deallocate(stressDivergenceU_test)
    deallocate(stressDivergenceV_test)
    deallocate(stressDivergenceU_diff)
    deallocate(stressDivergenceV_diff)

    deallocate(stress11)
    deallocate(stress22)
    deallocate(stress12)

    deallocate(lonVertexRotated)
    deallocate(latVertexRotated)
    deallocate(lonCellRotated)
    deallocate(latCellRotated)

  end subroutine spherical_test_divergence_stress_weak_polygon

  !-------------------------------------------------------------

end module cice_dyn_hex_weak


