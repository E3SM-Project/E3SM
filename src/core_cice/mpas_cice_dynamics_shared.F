module cice_dynamics_shared

  use mpas_grid_types

  implicit none


  private

  ! Weak formulation
  public :: init_evp, &
            boundary_cells, &
            boundary_vertices, &
            boundary_edges, &
            interior_vertices, &
            calc_cell_vertices_at_vertex, &
            normal_vectors, &
            ice_strength, &
            evp_constitutive_relation, &
            linear_constitutive_relation, &
            air_stress, &
            ocean_stress, &
            ocean_stress_coefficient, &
            solve_velocity, &
            surface_tilt, &
            solve_linear_basis_system, &
            interpolate_cell_to_vertex, &
            principal_stresses, &
            grid_rotation_forward, &
            grid_rotation_forward_nochange, &
            latlon_vector_rotation_forward, &
            latlon_vector_rotation_backward, &
            latlon_from_xyz, &
            dot_product_3space

  real(kind=RKIND), parameter, public :: &
       sinOceanTurningAngle = 0.0_RKIND, & ! northern hemisphere
       cosOceanTurningAngle = 1.0_RKIND, & ! northern hemisphere
       areaMinimum = 0.001_RKIND, &
       massMinimum = 0.01_RKIND
  
  real(kind=RKIND), parameter :: &
       eccentricity = 2.0_RKIND, &
       dampingTimescaleParameter = 0.36_RKIND, &
       puny = 1.0e-11_RKIND

  real(kind=RKIND), public :: &
       dtDynamics, &
       dtElastic, &
       dampingTimescale, &
       elapsedTime, &
       evpDampingCriterion

  integer, parameter, public :: &
       nElasticSubcycle = 120 ! will be from namelist

contains

  !-------------------------------------------------------------
  ! Initialization
  !-------------------------------------------------------------

  subroutine init_evp(dt)

    real(kind=RKIND), intent(in) :: dt

    dtDynamics = dt

    dtElastic = dtDynamics / real(nElasticSubcycle,RKIND)

    dampingTimescale = dampingTimescaleParameter * dtDynamics

    evpDampingCriterion = (1230.0_RKIND * dampingTimescale) / dtElastic**2 

    elapsedTime = 0.0_RKIND

  end subroutine init_evp

  !-------------------------------------------------------------

  subroutine ice_strength(mesh,        icePressure,   &
                          iceAreaCell, iceVolumeCell, &
                          solveStress)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         icePressure

    real(kind=RKIND), dimension(:), intent(in) :: &
         iceAreaCell, &
         iceVolumeCell

    integer, dimension(:), intent(in) :: &
         solveStress

    integer :: &
         iCell

    real(kind=RKIND), parameter :: &
         Pstar = 2.75e4_RKIND, & ! constant in Hibler strength formula 
         Cstar = 20.0_RKIND      ! constant in Hibler strength formula 

    do iCell = 1, mesh % nCells

       if (solveStress(iCell) == 1) then

          icePressure(iCell) = Pstar * iceVolumeCell(iCell) * exp(-Cstar*(1.0_RKIND-iceAreaCell(iCell)))

       else

          icePressure(iCell) = 0.0_RKIND

       endif ! solveStress

    enddo ! iCell

  end subroutine ice_strength

  !-------------------------------------------------------------

  subroutine interpolate_cell_to_vertex(mesh, variableVertex, variableCell)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         variableVertex

    real(kind=RKIND), dimension(:), intent(in) :: &
         variableCell

    real(kind=RKIND) :: &
         totalArea

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell

    do iVertex = 1, mesh % nVertices
    
       variableVertex(iVertex) = 0.0_RKIND
       totalArea = 0.0_RKIND

       do iVertexDegree = 1, mesh % vertexDegree

          iCell = mesh % cellsOnVertex % array(iVertexDegree,iVertex) 

          !variableVertex(iVertex) = variableVertex(iVertex) + mesh % kiteAreasOnVertex % array(iVertexDegree, iVertex) * variableCell(iCell)

          variableVertex(iVertex) = variableVertex(iVertex) + mesh % areaCell % array(iCell) * variableCell(iCell)
          totalArea = totalArea + mesh % areaCell % array(iCell)

       enddo ! iVertexDegree

       variableVertex(iVertex) = variableVertex(iVertex) / totalArea

       !variableVertex(iVertex) = variableVertex(iVertex) / mesh % areaTriangle % array(iVertex)

    enddo ! iVertex    

  end subroutine interpolate_cell_to_vertex

  !-------------------------------------------------------------

  subroutine evp_constitutive_relation(stress11,    stress22,    stress12, &
                                       strain11,    strain22,    strain12, &
                                       icePressure, replacementPressure,   &
                                       areaCell,    dtElastic,   evpDamping)

    real(kind=RKIND), intent(inout) :: &
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND), intent(in) :: &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND), intent(in) :: &
         icePressure, & 
         dtElastic,   &
         areaCell

    real(kind=RKIND), intent(out) :: &
         replacementPressure

    logical, intent(in) :: &
         evpDamping

    real(kind=RKIND) :: &
         strainDivergence,   &
         strainTension,      &
         strainShearing,     &
         stress1,            &
         stress2,            &
         Delta,              &
         pressureCoefficient
    
    ! convert from stress11 to stress1 etc
    strainDivergence = strain11 + strain22
    strainTension    = strain11 - strain22
    strainShearing   = strain12 * 2.0_RKIND
    
    stress1 = stress11 + stress22
    stress2 = stress11 - stress22

    ! perform the constituitive relation
    Delta = sqrt(strainDivergence**2 + &
         (strainTension**2 + strainShearing**2) / eccentricity**2)

    if (evpDamping) then
       pressureCoefficient = icePressure  / max(Delta,4.0_RKIND*puny)
       replacementPressure = pressureCoefficient * Delta
       pressureCoefficient = min(pressureCoefficient,evpDampingCriterion*areaCell)
    else
       pressureCoefficient = icePressure / max(Delta,puny)
       replacementPressure = pressureCoefficient * Delta
    endif

    pressureCoefficient = (pressureCoefficient * dtElastic) / (2.0_RKIND * dampingTimescale)

    stress1  = (stress1  + pressureCoefficient * (strainDivergence - Delta)) / &
               (1.0_RKIND + (0.5_RKIND * dtElastic) / dampingTimescale)

    stress2  = (stress2  + pressureCoefficient * strainTension) / &
               (1.0_RKIND + (0.5_RKIND * dtElastic * eccentricity**2) / dampingTimescale)

    stress12 = (stress12 + pressureCoefficient * strainShearing * 0.5_RKIND) / &
               (1.0_RKIND + (0.5_RKIND * dtElastic * eccentricity**2) / dampingTimescale)

    write(*,*) pressureCoefficient, strainTension, (1.0_RKIND + (0.5_RKIND * dtElastic * eccentricity**2) / dampingTimescale)
    write(*,*) stress1, stress2, stress12 ; stop

    ! convert back
    stress11 = 0.5_RKIND * (stress1 + stress2)
    stress22 = 0.5_RKIND * (stress1 - stress2)

  end subroutine evp_constitutive_relation

  !-------------------------------------------------------------

  subroutine solve_velocity(nPoints,           solvePoints,       &
                            uVelocity,         vVelocity,         &
                            mass,              coriolisParameter, &
                            stressDivergenceU, stressDivergenceV, &
                            airStressU,        airStressV,        &
                            surfaceTiltForceU, surfaceTiltForceV, &
                            oceanStressU,      oceanStressV,      &
                            oceanStressCoeff,  dtElastic,mesh)

    integer, intent(in) :: &
         nPoints

    integer, dimension(:), intent(in) :: &
         solvePoints

    real(kind=RKIND), dimension(:), intent(out) :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:), intent(in) :: &
         mass,              &
         coriolisParameter, &
         stressDivergenceU, &
         stressDivergenceV, &
         airStressU,        &
         airStressV,        &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         oceanStressU,      &
         oceanStressV,      &
         oceanStressCoeff

    real(kind=RKIND), intent(in) :: &
         dtElastic
    type(mesh_type), intent(in) :: mesh
    real(kind=RKIND), dimension(2) :: &
         rightHandSide

    real(kind=RKIND), dimension(2,2) :: &
         leftMatrix

    real(kind=RKIND) :: &
         solutionDenominator

    integer :: &
         iPoint,i,j

    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then

          ! U equation
          leftMatrix(1,1) =  mass(iPoint) / dtElastic                 + oceanStressCoeff(iPoint) * cosOceanTurningAngle
          leftMatrix(1,2) = -mass(iPoint) * coriolisParameter(iPoint) - oceanStressCoeff(iPoint) * sinOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iPoint))

          ! V equation
          leftMatrix(2,1) =  mass(iPoint) * coriolisParameter(iPoint) + oceanStressCoeff(iPoint) * sinOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iPoint))
          leftMatrix(2,2) =  mass(iPoint) / dtElastic                 + oceanStressCoeff(iPoint) * cosOceanTurningAngle

          ! right hand side of matrix solve
          rightHandSide(1) = stressDivergenceU(iPoint) + airStressU(iPoint) + surfaceTiltForceU(iPoint) + &
                             oceanStressCoeff(iPoint) * oceanStressU(iPoint) + (mass(iPoint) * uVelocity(iPoint)) / dtElastic

          rightHandSide(2) = stressDivergenceV(iPoint) + airStressV(iPoint) + surfaceTiltForceV(iPoint) + &
                             oceanStressCoeff(iPoint) * oceanStressV(iPoint) + (mass(iPoint) * vVelocity(iPoint)) / dtElastic

          ! solve the equation
          solutionDenominator = leftMatrix(1,1) * leftMatrix(2,2) - leftMatrix(1,2) * leftMatrix(2,1)

          uVelocity(iPoint) = (leftMatrix(2,2) * rightHandSide(1) - leftMatrix(1,2) * rightHandSide(2)) / solutionDenominator
          vVelocity(iPoint) = (leftMatrix(1,1) * rightHandSide(2) - leftMatrix(2,1) * rightHandSide(1)) / solutionDenominator

          write(70,*) mesh % POPindxiv % array(iPoint), mesh % POPindxjv % array(iPoint), &
               uVelocity(iPoint), coriolisParameter(iPoint)

       else

          uVelocity(iPoint) = 0.0_RKIND
          vVelocity(iPoint) = 0.0_RKIND

       endif ! solvePoints

    enddo ! iCell

  end subroutine solve_velocity

  !-------------------------------------------------------------

  subroutine air_stress(airStressU, airStressV, uAirVelocity, vAirVelocity, airDensity, iceArea)

    real(kind=RKIND), dimension(:), intent(out) :: &
         airStressU, &
         airStressV

    real(kind=RKIND), dimension(:), intent(in) :: &
         uAirVelocity, &
         vAirVelocity, &
         airDensity, &
         iceArea

    real(kind=RKIND) :: &
         windSpeed

    integer :: &
         nPoints, &
         iPoint

    real(kind=RKIND), parameter :: &
         airStressCoeff = 0.0012_RKIND 
    
    nPoints = size(airStressU,1)

    do iPoint = 1, nPoints

       windSpeed = sqrt(uAirVelocity(iPoint)**2 + vAirVelocity(iPoint)**2)

       airStressU(iPoint) = airDensity(iPoint) * windSpeed * airStressCoeff * uAirVelocity(iPoint) * iceArea(iPoint)
       airStressV(iPoint) = airDensity(iPoint) * windSpeed * airStressCoeff * vAirVelocity(iPoint) * iceArea(iPoint)

    enddo ! iPoint

  end subroutine air_stress

  !-------------------------------------------------------------

  subroutine ocean_stress(oceanStressU,      oceanStressV,   &
                          uOceanVelocity,    vOceanVelocity, &
                          coriolisParameter, solvePoints)

    real(kind=RKIND), dimension(:), intent(out) :: &
         oceanStressU, &
         oceanStressV

    real(kind=RKIND), dimension(:), intent(in) :: &
         uOceanVelocity, &
         vOceanVelocity, &
         coriolisParameter

    integer, dimension(:), intent(in) :: &
         solvePoints

    integer :: &
         nPoints, &
         iPoint

    nPoints = size(oceanStressU,1)

    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then

          oceanStressU(iPoint) = uOceanVelocity(iPoint) * cosOceanTurningAngle - &
                                 vOceanVelocity(iPoint) * sinOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iPoint))
          oceanStressV(iPoint) = uOceanVelocity(iPoint) * sinOceanTurningAngle + &
                                 vOceanVelocity(iPoint) * cosOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iPoint))
          
       else

          oceanStressU(iPoint) = 0.0_RKIND
          oceanStressV(iPoint) = 0.0_RKIND

       endif ! solvePoints

    enddo ! iPoint

  end subroutine ocean_stress

  !-------------------------------------------------------------

  subroutine ocean_stress_coefficient(oceanStressCoeff, &
                                      uOceanVelocity, vOceanVelocity, &
                                      uVelocity, vVelocity, &
                                      iceArea, solvePoints) 

    real(kind=RKIND), dimension(:), intent(out) :: &
         oceanStressCoeff

    real(kind=RKIND), dimension(:), intent(in) :: &
         uOceanVelocity, &
         vOceanVelocity, &
         uVelocity, &
         vVelocity, &
         iceArea

    integer, dimension(:), intent(in) :: &
         solvePoints

    real(kind=RKIND), parameter :: &
         densityWater      = 1026.0_RKIND, &
         iceOceanDragCoeff = 0.00536_RKIND

    integer :: &
         nPoints, &
         iPoint

    nPoints = size(oceanStressCoeff,1)

    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then

          oceanStressCoeff(iPoint) = iceOceanDragCoeff * densityWater * iceArea(iPoint) * &
                                     sqrt((uOceanVelocity(iPoint) - uVelocity(iPoint))**2 + &
                                          (vOceanVelocity(iPoint) - vVelocity(iPoint))**2)
       else

          oceanStressCoeff(iPoint) = 0.0_RKIND

       endif

    enddo ! iPoint

  end subroutine ocean_stress_coefficient
  
  !-------------------------------------------------------------

  subroutine surface_tilt(surfaceTiltForceU, surfaceTiltForceV, &
                          uOceanVelocity,    vOceanVelocity,    &
                          totalMass,         coriolisParameter, &
                          solvePoints)

    real(kind=RKIND), dimension(:), intent(out) :: &
         surfaceTiltForceU, &
         surfaceTiltForceV

    real(kind=RKIND), dimension(:), intent(in) :: &
         uOceanVelocity,  &
         vOceanVelocity,  &
         totalMass,       &
         coriolisParameter

    integer, dimension(:), intent(in) :: &
         solvePoints

    integer :: &
         nPoints, &
         iPoint

    nPoints = size(surfaceTiltForceU,1)

    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then

          surfaceTiltForceU(iPoint) = -coriolisParameter(iPoint) * totalMass(iPoint) * vOceanVelocity(iPoint)
          surfaceTiltForceV(iPoint) =  coriolisParameter(iPoint) * totalMass(iPoint) * uOceanVelocity(iPoint)

       else

          surfaceTiltForceU(iPoint) = 0.0_RKIND
          surfaceTiltForceV(iPoint) = 0.0_RKIND

       endif ! solvePoints

    enddo ! iPoint

  end subroutine surface_tilt

  !-------------------------------------------------------------

  subroutine principal_stresses(principalStress1, principalStress2, &
                                stress11, stress22, stress12, &
                                replacementPressure)

    real(kind=RKIND), dimension(:), intent(out) :: &
         principalStress1, &
         principalStress2

    real(kind=RKIND), dimension(:), intent(in) :: &
         stress11, &
         stress22, &
         stress12, &
         replacementPressure

    real(kind=RKIND) :: &
         sqrtContents

    integer :: &
         iPoint, &
         nPoints

    nPoints = size(principalStress1)

    open(48,file='principal.txt')

    do iPoint = 1, nPoints

       if (replacementPressure(iPoint) > puny) then

          sqrtContents = (stress11(iPoint) + stress22(iPoint))**2 - &
                         4.0_RKIND * stress11(iPoint) * stress22(iPoint) + &
                         4.0_RKIND * stress12(iPoint)**2

          principalStress1(iPoint) = 0.5_RKIND * (stress11(iPoint) + stress22(iPoint)) + 0.5_RKIND * sqrt(sqrtContents)
          principalStress2(iPoint) = 0.5_RKIND * (stress11(iPoint) + stress22(iPoint)) - 0.5_RKIND * sqrt(sqrtContents)

          principalStress1(iPoint) = principalStress1(iPoint) / replacementPressure(iPoint)
          principalStress2(iPoint) = principalStress2(iPoint) / replacementPressure(iPoint)

       else
          
          principalStress1(iPoint) = 1.0e30_RKIND
          principalStress2(iPoint) = 1.0e30_RKIND

       endif

       write(48,*) principalStress1(iPoint), principalStress2(iPoint)

    enddo ! iPoint

    close(48)

  end subroutine principal_stresses

  !-------------------------------------------------------------
  ! Testing
  !-------------------------------------------------------------

  subroutine linear_constitutive_relation(stress11, stress22, stress12, &
                                          strain11, strain22, strain12)

    use cice_testing, only: lambda

    real(kind=RKIND), intent(out) :: &
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND), intent(in) :: &
         strain11, &
         strain22, &
         strain12

    stress11 = lambda * strain11
    stress22 = lambda * strain22
    stress12 = lambda * strain12

  end subroutine linear_constitutive_relation

  !-------------------------------------------------------------
  ! Legacy
  !-------------------------------------------------------------

  subroutine boundary_cells(mesh, boundaryCell, boundaryCell2, interiorCell) 

    type(mesh_type), intent(in) :: mesh

    integer, dimension(:), intent(out) :: &
         boundaryCell, &
         boundaryCell2, &
         interiorCell
         
    integer :: &
         iCell, &
         iCellOnCell, &
         iAdjacentCell, &
         nAdjacentCells, &
         nEdgesOnCell

    do iCell = 1, mesh % nCells

       boundaryCell(iCell) = 0
       interiorCell(iCell) = 0

       nEdgesOnCell = mesh % nEdgesOnCell % array(iCell)
       nAdjacentCells = 0

       do iCellOnCell = 1, nEdgesOnCell

          iAdjacentCell = mesh % cellsOnCell % array(iCellOnCell,iCell)

          if (iAdjacentCell <= mesh % nCells) then

             nAdjacentCells = nAdjacentCells + 1

          endif

          if (nAdjacentCells == nEdgesOnCell) interiorCell(iCell) = 1
          if (nAdjacentCells <  nEdgesOnCell) boundaryCell(iCell) = 1

       enddo ! iCellOnCell

    enddo ! iCell

    do iCell = 1, mesh % nCells

       boundaryCell2(iCell) = 0

       if (boundaryCell(iCell) == 0) then

          do iCellOnCell = 1, mesh % nEdgesOnCell % array(iCell)

             iAdjacentCell = mesh % cellsOnCell % array(iCellOnCell,iCell)

             if (boundaryCell(iAdjacentCell) == 1) then

                boundaryCell2(iCell) = 1

             endif

          enddo ! iCellOnCell

       endif

    enddo ! iCell

  end subroutine boundary_cells

  !-------------------------------------------------------------

  subroutine boundary_vertices(mesh, boundaryVertex, interiorVertex, interiorVertex2, boundaryCell)

    type(mesh_type), intent(in) :: mesh

    integer, dimension(:), intent(out) :: &
         boundaryVertex, & ! boundary vertices of boundary vertex triangles
         interiorVertex, & ! all normal shaped interior vertices
         interiorVertex2

    integer, dimension(:), intent(in) :: &
         boundaryCell

    integer :: &
         nInteriorAdjacentCells

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell

    logical :: &
         linterior

    ! boundary vertices
    do iVertex = 1, mesh % nVertices

       boundaryVertex(iVertex) = 0
       interiorVertex(iVertex) = 0
       
       nInteriorAdjacentCells = 0

       do iVertexDegree = 1, mesh % vertexDegree

          iCell = mesh % cellsOnVertex % array(iVertexDegree, iVertex)
          
          if (iCell >= 1 .and. iCell <= mesh % nCells) then
             
             nInteriorAdjacentCells = nInteriorAdjacentCells + 1
             
          endif
          
       enddo ! iVertexDegree

       ! boundary vertices that have triangles
       if (nInteriorAdjacentCells == 2) then

          boundaryVertex(iVertex) = 1
          
       endif

       ! vertex points we directly calculate on
       if (nInteriorAdjacentCells == 3) then

          interiorVertex(iVertex) = 1
          
       endif
       
    enddo ! iVertex


    ! second interior vertices
    do iVertex = 1, mesh % nVertices

       interiorVertex2(iVertex) = 0

       linterior = .true.

       do iVertexDegree = 1, mesh % vertexDegree

          iCell = mesh % cellsOnVertex % array(iVertexDegree, iVertex)
          
          if (boundaryCell(iCell) == 1) then
             
             linterior = .false.
             
          endif
          
       enddo ! iVertexDegree

       if (linterior) interiorVertex2(iVertex) = 1

    enddo ! iVertex

  end subroutine boundary_vertices

  !-------------------------------------------------------------

  subroutine boundary_edges(mesh, boundaryEdge)

    type(mesh_type), intent(in) :: mesh

    integer, dimension(:), intent(out) :: &
         boundaryEdge

    integer :: &
         nInteriorAdjacentCells

    integer :: &
         iCell, &
         iEdge, &
         iCellsOnEdge
    
    ! boundary edges
    do iEdge = 1, mesh % nEdges
       
       do iCellsOnEdge = 1, 2
          
          iCell = mesh % cellsOnEdge % array(iCellsOnEdge,iEdge)
          
          if (iCell > mesh % nCells) then
             boundaryEdge(iEdge) = 1
          else
             boundaryEdge(iEdge) = 0
          endif
          
       enddo ! iCellsOnEdge

    enddo ! iEdge
    
  end subroutine boundary_edges

  !-------------------------------------------------------------

  subroutine interior_vertices(mesh, interiorVertex)

    type(mesh_type), intent(in) :: mesh

    integer, dimension(:), intent(out) :: &
         interiorVertex

    integer :: &
         nInteriorAdjacentCells

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell

    ! boundary vertices
    do iVertex = 1, mesh % nVertices

       interiorVertex(iVertex) = 0
       
       nInteriorAdjacentCells = 0

       do iVertexDegree = 1, mesh % vertexDegree

          iCell = mesh % cellsOnVertex % array(iVertexDegree, iVertex)

          if (iCell >= 1 .and. iCell <= mesh % nCells) then
             
             nInteriorAdjacentCells = nInteriorAdjacentCells + 1
             
          endif
          
       enddo ! iVertexDegree

       ! vertex points we directly calculate on
       if (nInteriorAdjacentCells == mesh % vertexDegree) then

          interiorVertex(iVertex) = 1
          
       endif
       
    enddo ! iVertex

  end subroutine interior_vertices

  !-------------------------------------------------------------

   subroutine calc_cell_vertices_at_vertex(mesh, cellVerticesAtVertex)

    type(mesh_type), intent(in) :: mesh

    integer, dimension(:,:), intent(out) :: &
         cellVerticesAtVertex

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell, &
         iVertexOnCell, &
         jVertex

    do iVertex = 1, mesh % nVertices

       do iVertexDegree = 1, mesh % vertexDegree

          cellVerticesAtVertex(iVertexDegree,iVertex) = 0

          iCell = mesh % cellsOnVertex % array(iVertexDegree, iVertex)

          do iVertexOnCell = 1, mesh % nEdgesOnCell % array(iCell)

             jVertex = mesh % verticesOnCell % array(iVertexOnCell,iCell)

             if (iVertex == jVertex) then

                cellVerticesAtVertex(iVertexDegree,iVertex) = iVertexOnCell

             endif

          enddo ! iVertexOnCell

       enddo ! iVertexDegree

    enddo ! iVertex

   end subroutine calc_cell_vertices_at_vertex

  !-------------------------------------------------------------
  ! normal vectors
  !-------------------------------------------------------------

  subroutine normal_vectors(mesh, &
                            normalVectorPolygon, normalVectorTriangle, normalVectorTriangleCorners, &
                            latCellRotated,      latVertexRotated,     &
                            interiorVertex,      gridCurvature)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon, &
         normalVectorTriangle

    real(kind=RKIND), dimension(:,:,:,:), intent(out) :: &
         normalVectorTriangleCorners

    real(kind=RKIND), dimension(:), intent(out) :: &
         latCellRotated, &
         latVertexRotated

    integer, dimension(:), intent(in) :: &
         interiorVertex

    character(len=*), intent(in) :: &
         gridCurvature

    if (trim(gridCurvature) == "planar") then

       call normal_vectors_planar_polygon(mesh, normalVectorPolygon)
       call normal_vectors_planar_triangle(mesh, normalVectorTriangle, interiorVertex)

    else if (trim(gridCurvature) == "spherical") then

       call normal_vectors_spherical_polygon_metric(mesh, normalVectorPolygon, latCellRotated)
       call normal_vectors_spherical_triangle_metric(mesh, normalVectorTriangle, latVertexRotated, interiorVertex)

    endif

  end subroutine normal_vectors

  !-------------------------------------------------------------

  subroutine normal_vectors_planar_polygon(mesh, normalVectorPolygon)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon

    integer :: &
         iCell, &
         iEdgeOnCell, &
         iEdge

    real(kind=RKIND) :: &
         dx, dy

    do iCell = 1, mesh % nCells

       do iEdgeOnCell = 1, mesh % nEdgesOnCell % array(iCell)
          
          iEdge = mesh % edgesOnCell % array(iEdgeOnCell,iCell)

          dx = mesh % xEdge % array(iEdge) - mesh % xCell % array(iCell)
          dy = mesh % yEdge % array(iEdge) - mesh % yCell % array(iCell)

          normalVectorPolygon(1,iEdgeOnCell,iCell) = dx / sqrt(dx**2 + dy**2)
          normalVectorPolygon(2,iEdgeOnCell,iCell) = dy / sqrt(dx**2 + dy**2)

       enddo ! iEdgeOnCell

    enddo ! iCell

  end subroutine normal_vectors_planar_polygon

  !-------------------------------------------------------------

  subroutine normal_vectors_planar_triangle(mesh, normalVectorTriangle, interiorVertex)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorTriangle

    integer, dimension(:), intent(in) :: &
         interiorVertex

    integer :: &
         iVertex, &
         iVertexDegree, &
         iEdge

    real(kind=RKIND) :: &
         dx, dy

    do iVertex = 1, mesh % nVertices

       normalVectorTriangle(:,:,iVertex) = 0.0_RKIND

       if (interiorVertex(iVertex) == 1) then

          do iVertexDegree = 1, mesh % vertexDegree
             
             iEdge = mesh % edgesOnVertex % array(iVertexDegree,iVertex)
             
             dx = mesh % xEdge % array(iEdge) - mesh % xVertex % array(iVertex)
             dy = mesh % yEdge % array(iEdge) - mesh % yVertex % array(iVertex)
             
             normalVectorTriangle(1,iVertexDegree,iVertex) = dx / sqrt(dx**2 + dy**2)
             normalVectorTriangle(2,iVertexDegree,iVertex) = dy / sqrt(dx**2 + dy**2)
             
          enddo ! iVertexDegree

       endif ! interiorVertex

    enddo ! iVertex

  end subroutine normal_vectors_planar_triangle

  !-------------------------------------------------------------

  subroutine normal_vectors_spherical_polygon_metric(mesh, normalVectorPolygon, latCellRotated)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon

    real(kind=RKIND), dimension(:), intent(out) :: &
         latCellRotated
    
    real(kind=RKIND), dimension(3) :: &
         cellCentreRotated, &
         cellCentreEquator, &
         edgeRotated, &
         vertexRotated1, &
         vertexRotated2, &
         edgeEquator, &
         vertexEquator1, &
         vertexEquator2, &
         vertexVector, &
         normalGreatCircle, &
         eastwardsVector

    real(kind=RKIND) :: &
         lonCellRotated, &
         normalGreatCircleNorm

    real(kind=RKIND), dimension(3,3) :: &
         yRotationMatrix, &
         zRotationMatrix

    integer :: &
         iCell, &
         iEdgeOnCell, &
         iEdge, &
         iVertex1, &
         iVertex2

    real(kind=RKIND), parameter :: &
         rad_to_deg = 180.0_RKIND / 3.14159265359_RKIND

    yRotationMatrix = 0.0_RKIND
    zRotationMatrix = 0.0_RKIND

    yRotationMatrix(2,2) = 1.0_RKIND
    zRotationMatrix(3,3) = 1.0_RKIND

    do iCell = 1, mesh % nCells

       ! rotate the cell centre to rotated geographical grid
       call grid_rotation_forward(cellCentreRotated(1),        cellCentreRotated(2),        cellCentreRotated(3), &
                                  mesh % xCell % array(iCell), mesh % yCell % array(iCell), mesh % zCell % array(iCell))

       ! calculate lon and lat of the cell centre in the rotated coordinate system
       lonCellRotated        = atan2(cellCentreRotated(2), cellCentreRotated(1))
       latCellRotated(iCell) = asin(cellCentreRotated(3) / mesh % sphere_radius)

       ! rotate the cell centre onto the rotated equator
       yRotationMatrix(1,1) =  cos(latCellRotated(iCell))
       yRotationMatrix(1,3) =  sin(latCellRotated(iCell))
       yRotationMatrix(3,1) = -sin(latCellRotated(iCell))
       yRotationMatrix(3,3) =  cos(latCellRotated(iCell))

       zRotationMatrix(1,1) =  cos(-lonCellRotated)
       zRotationMatrix(1,2) = -sin(-lonCellRotated)
       zRotationMatrix(2,1) =  sin(-lonCellRotated)
       zRotationMatrix(2,2) =  cos(-lonCellRotated)

       ! this should be (r, 0, 0)
       !cellCentreEquator = matmul(yRotationMatrix, matmul(zRotationMatrix, cellCentreRotated))
       
       do iEdgeOnCell = 1, mesh % nEdgesOnCell % array(iCell)
          
          iEdge = mesh % edgesOnCell % array(iEdgeOnCell,iCell)

          ! vector in great circle plane
          iVertex1 = mesh % verticesOnEdge % array(1,iEdge)
          iVertex2 = mesh % verticesOnEdge % array(2,iEdge)

          ! perform grid rotations
          call grid_rotation_forward(edgeRotated(1),                   edgeRotated(2),                   edgeRotated(3), &
                                     mesh % xEdge % array(iEdge),      mesh % yEdge % array(iEdge),      mesh % zEdge % array(iEdge))
          
          call grid_rotation_forward(vertexRotated1(1),                vertexRotated1(2),                vertexRotated1(3), &
                                     mesh % xVertex % array(iVertex1), mesh % yVertex % array(iVertex1), mesh % zVertex % array(iVertex1))

          call grid_rotation_forward(vertexRotated2(1),                vertexRotated2(2),                vertexRotated2(3), &
                                     mesh % xVertex % array(iVertex2), mesh % yVertex % array(iVertex2), mesh % zVertex % array(iVertex2))

          ! rotate to equator in new coords
          edgeEquator    = matmul(yRotationMatrix, matmul(zRotationMatrix, edgeRotated))
          vertexEquator1 = matmul(yRotationMatrix, matmul(zRotationMatrix, vertexRotated1))
          vertexEquator2 = matmul(yRotationMatrix, matmul(zRotationMatrix, vertexRotated2))

          ! vector joining vertices
          vertexVector = vertexEquator2 - vertexEquator1

          ! form great circle plane normal vector
          call cross_product_3space(normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                    vertexVector(1),      vertexVector(2),      vertexVector(3), &
                                    edgeEquator(1),       edgeEquator(2),       edgeEquator(3))

          if (iCell == mesh % cellsOnEdge % array(2,iEdge)) then
             normalGreatCircle = -1.0_RKIND * normalGreatCircle
          endif

          ! normalize normal vector
          normalGreatCircleNorm = sqrt(normalGreatCircle(1)**2 + normalGreatCircle(2)**2 + normalGreatCircle(3)**2)
          normalGreatCircle(1) = normalGreatCircle(1) / normalGreatCircleNorm
          normalGreatCircle(2) = normalGreatCircle(2) / normalGreatCircleNorm
          normalGreatCircle(3) = normalGreatCircle(3) / normalGreatCircleNorm

          ! eastwards vector at edge
          eastwardsVector(1) = -edgeEquator(2)
          eastwardsVector(2) =  edgeEquator(1)
          eastwardsVector(3) = 0.0_RKIND
          eastwardsVector = eastwardsVector / &
                            sqrt(eastwardsVector(1)**2 + eastwardsVector(2)**2)

          call dot_product_3space(normalVectorPolygon(1,iEdgeOnCell,iCell), &
                                  normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                  eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
          normalVectorPolygon(2,iEdgeOnCell,iCell) = sign(1.0_RKIND,normalGreatCircle(3))&
                                                   * sqrt(1.0_RKIND - max(min(normalVectorPolygon(1,iEdgeOnCell,iCell),1.0_RKIND),-1.0_RKIND)**2)

       enddo ! iEdgeOnCell
       
    enddo ! iCell

  end subroutine normal_vectors_spherical_polygon_metric

  !-------------------------------------------------------------

  subroutine normal_vectors_spherical_polygon(mesh, normalVectorPolygon)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon

    real(kind=RKIND), dimension(3) :: &
         normalGreatCircle, &
         eastwardsVector

    real(kind=RKIND) :: &
         xEdge,    yEdge,    zEdge,    &
         xVertex1, yVertex1, zVertex1, &
         xVertex2, yVertex2, zVertex2

    real(kind=RKIND) :: &
         normalGreatCircleNorm

    integer :: &
         iCell, &
         iEdgeOnCell, &
         iEdge, &
         iVertex1, &
         iVertex2

    real(kind=RKIND), parameter :: &
         rad_to_deg = 180.0_RKIND / 3.14159265359_RKIND

    do iCell = 1, mesh % nCells

       do iEdgeOnCell = 1, mesh % nEdgesOnCell % array(iCell)
          
          iEdge = mesh % edgesOnCell % array(iEdgeOnCell,iCell)

          ! vector in great circle plane
          iVertex1 = mesh % verticesOnEdge % array(1,iEdge)
          iVertex2 = mesh % verticesOnEdge % array(2,iEdge)

          ! perform grid rotations
          call grid_rotation_forward_nochange(xEdge, yEdge, zEdge, &
                                     mesh % xEdge % array(iEdge), mesh % yEdge % array(iEdge), mesh % zEdge % array(iEdge))
          
          call grid_rotation_forward_nochange(xVertex1, yVertex1, zVertex1, &
                                     mesh % xVertex % array(iVertex1), mesh % yVertex % array(iVertex1), mesh % zVertex % array(iVertex1))

          call grid_rotation_forward_nochange(xVertex2, yVertex2, zVertex2, &
                                     mesh % xVertex % array(iVertex2), mesh % yVertex % array(iVertex2), mesh % zVertex % array(iVertex2))

          ! form great circle plane normal vector
          call cross_product_3space(normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                    xVertex2 - xVertex1,  yVertex2 - yVertex1,  zVertex2 - zVertex1, &
                                    xEdge,                yEdge,                zEdge)

          if (iCell == mesh % cellsOnEdge % array(2,iEdge)) then
             normalGreatCircle = -1.0_RKIND * normalGreatCircle
          endif

          ! normalize normal vector
          normalGreatCircleNorm = sqrt(normalGreatCircle(1)**2 + normalGreatCircle(2)**2 + normalGreatCircle(3)**2)
          normalGreatCircle(1) = normalGreatCircle(1) / normalGreatCircleNorm
          normalGreatCircle(2) = normalGreatCircle(2) / normalGreatCircleNorm
          normalGreatCircle(3) = normalGreatCircle(3) / normalGreatCircleNorm

          ! eastwards vector at edge
          eastwardsVector(1) = -yEdge
          eastwardsVector(2) =  xEdge
          eastwardsVector(3) = 0.0_RKIND
          eastwardsVector = eastwardsVector / &
                            sqrt(eastwardsVector(1)**2 + eastwardsVector(2)**2)

          call dot_product_3space(normalVectorPolygon(1,iEdgeOnCell,iCell), &
                                  normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                  eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
          normalVectorPolygon(2,iEdgeOnCell,iCell) = sign(1.0_RKIND,normalGreatCircle(3))&
                                                   * sqrt(1.0_RKIND - max(min(normalVectorPolygon(1,iEdgeOnCell,iCell),1.0_RKIND),-1.0_RKIND)**2)

       enddo ! iEdgeOnCell
       
    enddo ! iCell

  end subroutine normal_vectors_spherical_polygon

  !-------------------------------------------------------------

  subroutine normal_vectors_spherical_triangle_corners_metric(mesh, normalVectorTriangleCorners, latVertexRotated, interiorVertex)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:,:), intent(out) :: &
         normalVectorTriangleCorners

    real(kind=RKIND), dimension(:), intent(out) :: &
         latVertexRotated

    integer, dimension(:), intent(in) :: &
         interiorVertex

    real(kind=RKIND), dimension(3) :: &
         vertexRotated, &
         vertexEquator, &
         edgeRotated, &
         cellRotated1, &
         cellRotated2, &
         edgeEquator, &
         cellEquator1, &
         cellEquator2, &
         cellVector, &
         normalGreatCircle, &
         eastwardsVector

    real(kind=RKIND) :: &
         lonVertexRotated, &
         normalGreatCircleNorm

    real(kind=RKIND), dimension(3,3) :: &
         yRotationMatrix, &
         zRotationMatrix

    integer :: &
         iVertex, &
         iVertexDegree, &
         iEdge, &
         iCell1, &
         iCell2

    real(kind=RKIND), parameter :: &
         pi = 3.14159265359_RKIND, &
         deg_to_rad = pi / 180.0_RKIND, &
         rad_to_deg = 180.0_RKIND / pi

    yRotationMatrix = 0.0_RKIND
    zRotationMatrix = 0.0_RKIND

    yRotationMatrix(2,2) = 1.0_RKIND
    zRotationMatrix(3,3) = 1.0_RKIND

    do iVertex = 1, mesh % nVertices

       normalVectorTriangleCorners(:,:,:,iVertex) = 0.0_RKIND

       if (interiorVertex(iVertex) == 1) then

          ! rotate the cell centre to rotated geographical grid
          call grid_rotation_forward(vertexRotated(1),                vertexRotated(2),                vertexRotated(3), &
                                     mesh % xVertex % array(iVertex), mesh % yVertex % array(iVertex), mesh % zVertex % array(iVertex))
 
          ! calculate lon and lat of the cell centre in the rotated coordinate system
          lonVertexRotated          = atan2(vertexRotated(2), vertexRotated(1))
          latVertexRotated(iVertex) = asin(vertexRotated(3) / mesh % sphere_radius)
          
          ! rotate the cell centre onto the rotated equator
          yRotationMatrix(1,1) =  cos(latVertexRotated(iVertex))
          yRotationMatrix(1,3) =  sin(latVertexRotated(iVertex))
          yRotationMatrix(3,1) = -sin(latVertexRotated(iVertex))
          yRotationMatrix(3,3) =  cos(latVertexRotated(iVertex))
          
          zRotationMatrix(1,1) =  cos(-lonVertexRotated)
          zRotationMatrix(1,2) = -sin(-lonVertexRotated)
          zRotationMatrix(2,1) =  sin(-lonVertexRotated)
          zRotationMatrix(2,2) =  cos(-lonVertexRotated)

          ! this should be (r, 0, 0)
          !vertexEquator = matmul(yRotationMatrix, matmul(zRotationMatrix, vertexRotated))
          
          do iVertexDegree = 1, mesh % vertexDegree
             
             iEdge = mesh % edgesOnVertex % array(iVertexDegree,iVertex)

             ! vector in great circle plane
             iCell1 = mesh % cellsOnEdge % array(1,iEdge)
             iCell2 = mesh % cellsOnEdge % array(2,iEdge)

             ! perform grid rotations
             call grid_rotation_forward(cellRotated1(1),              cellRotated1(2),              cellRotated1(3), &
                                        mesh % xCell % array(iCell1), mesh % yCell % array(iCell1), mesh % zCell % array(iCell1))

             call grid_rotation_forward(cellRotated2(1),              cellRotated2(2),              cellRotated2(3), &
                                        mesh % xCell % array(iCell2), mesh % yCell % array(iCell2), mesh % zCell % array(iCell2))

             ! rotate to equator in new coords
             cellEquator1 = matmul(yRotationMatrix, matmul(zRotationMatrix, cellRotated1))
             cellEquator2 = matmul(yRotationMatrix, matmul(zRotationMatrix, cellRotated2))

             ! vector joining vertices
             cellVector = cellEquator2 - cellEquator1

             ! first cell 
             
             ! form great circle plane normal vector
             call cross_product_3space(normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                       cellVector(1),        cellVector(2),        cellVector(3), &
                                       cellEquator1(1),      cellEquator1(2),     cellEquator1(3))

             if (iVertex == mesh % verticesOnEdge % array(1,iEdge)) then
                normalGreatCircle = -1.0_RKIND * normalGreatCircle
             endif

             ! normalize normal vector
             normalGreatCircleNorm = sqrt(normalGreatCircle(1)**2 + normalGreatCircle(2)**2 + normalGreatCircle(3)**2)
             normalGreatCircle(1) = normalGreatCircle(1) / normalGreatCircleNorm
             normalGreatCircle(2) = normalGreatCircle(2) / normalGreatCircleNorm
             normalGreatCircle(3) = normalGreatCircle(3) / normalGreatCircleNorm
             
             ! eastwards vector at edge
             eastwardsVector(1) = -cellEquator1(2)
             eastwardsVector(2) =  cellEquator1(1)
             eastwardsVector(3) = 0.0_RKIND
             eastwardsVector = eastwardsVector / &
                               sqrt(eastwardsVector(1)**2 + eastwardsVector(2)**2)
             
             call dot_product_3space(normalVectorTriangleCorners(1,1,iVertexDegree,iVertex), &
                                     normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                     eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
             normalVectorTriangleCorners(2,1,iVertexDegree,iVertex) = sign(1.0_RKIND,normalGreatCircle(3)) &
                                                                    * sqrt(1.0_RKIND - max(min(normalVectorTriangleCorners(1,1,iVertexDegree,iVertex),1.0_RKIND),-1.0_RKIND)**2)

             ! second cell 
             
             ! form great circle plane normal vector
             call cross_product_3space(normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                       cellVector(1),        cellVector(2),        cellVector(3), &
                                       cellEquator2(1),      cellEquator2(2),     cellEquator2(3))

             if (iVertex == mesh % verticesOnEdge % array(1,iEdge)) then
                normalGreatCircle = -1.0_RKIND * normalGreatCircle
             endif

             ! normalize normal vector
             normalGreatCircleNorm = sqrt(normalGreatCircle(1)**2 + normalGreatCircle(2)**2 + normalGreatCircle(3)**2)
             normalGreatCircle(1) = normalGreatCircle(1) / normalGreatCircleNorm
             normalGreatCircle(2) = normalGreatCircle(2) / normalGreatCircleNorm
             normalGreatCircle(3) = normalGreatCircle(3) / normalGreatCircleNorm
             
             ! eastwards vector at edge
             eastwardsVector(1) = -cellEquator2(2)
             eastwardsVector(2) =  cellEquator2(1)
             eastwardsVector(3) = 0.0_RKIND
             eastwardsVector = eastwardsVector / &
                               sqrt(eastwardsVector(1)**2 + eastwardsVector(2)**2)
             
             call dot_product_3space(normalVectorTriangleCorners(1,2,iVertexDegree,iVertex), &
                                     normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                     eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
             normalVectorTriangleCorners(2,2,iVertexDegree,iVertex) = sign(1.0_RKIND,normalGreatCircle(3)) &
                                                                    * sqrt(1.0_RKIND - max(min(normalVectorTriangleCorners(1,2,iVertexDegree,iVertex),1.0_RKIND),-1.0_RKIND)**2)

          enddo ! iVertexDegree

       endif ! interiorVertex

    enddo ! iVertex

  end subroutine normal_vectors_spherical_triangle_corners_metric

  !-------------------------------------------------------------

  subroutine normal_vectors_spherical_triangle_metric(mesh, normalVectorTriangle, latVertexRotated, interiorVertex)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorTriangle

    real(kind=RKIND), dimension(:), intent(out) :: &
         latVertexRotated

    integer, dimension(:), intent(in) :: &
         interiorVertex

    real(kind=RKIND), dimension(3) :: &
         vertexRotated, &
         vertexEquator, &
         edgeRotated, &
         cellRotated1, &
         cellRotated2, &
         edgeEquator, &
         cellEquator1, &
         cellEquator2, &
         cellVector, &
         normalGreatCircle, &
         eastwardsVector

    real(kind=RKIND) :: &
         lonVertexRotated, &
         normalGreatCircleNorm

    real(kind=RKIND), dimension(3,3) :: &
         yRotationMatrix, &
         zRotationMatrix

    integer :: &
         iVertex, &
         iVertexDegree, &
         iEdge, &
         iCell1, &
         iCell2

    real(kind=RKIND), parameter :: &
         pi = 3.14159265359_RKIND, &
         deg_to_rad = pi / 180.0_RKIND, &
         rad_to_deg = 180.0_RKIND / pi

    yRotationMatrix = 0.0_RKIND
    zRotationMatrix = 0.0_RKIND

    yRotationMatrix(2,2) = 1.0_RKIND
    zRotationMatrix(3,3) = 1.0_RKIND

    do iVertex = 1, mesh % nVertices

       normalVectorTriangle(:,:,iVertex) = 0.0_RKIND

       if (interiorVertex(iVertex) == 1) then

          ! rotate the cell centre to rotated geographical grid
          call grid_rotation_forward(vertexRotated(1),                vertexRotated(2),                vertexRotated(3), &
                                     mesh % xVertex % array(iVertex), mesh % yVertex % array(iVertex), mesh % zVertex % array(iVertex))
 
          ! calculate lon and lat of the cell centre in the rotated coordinate system
          lonVertexRotated          = atan2(vertexRotated(2), vertexRotated(1))
          latVertexRotated(iVertex) = asin(vertexRotated(3) / mesh % sphere_radius)
          
          ! rotate the cell centre onto the rotated equator
          yRotationMatrix(1,1) =  cos(latVertexRotated(iVertex))
          yRotationMatrix(1,3) =  sin(latVertexRotated(iVertex))
          yRotationMatrix(3,1) = -sin(latVertexRotated(iVertex))
          yRotationMatrix(3,3) =  cos(latVertexRotated(iVertex))
          
          zRotationMatrix(1,1) =  cos(-lonVertexRotated)
          zRotationMatrix(1,2) = -sin(-lonVertexRotated)
          zRotationMatrix(2,1) =  sin(-lonVertexRotated)
          zRotationMatrix(2,2) =  cos(-lonVertexRotated)

          ! this should be (r, 0, 0)
          !vertexEquator = matmul(yRotationMatrix, matmul(zRotationMatrix, vertexRotated))
          
          do iVertexDegree = 1, mesh % vertexDegree
             
             iEdge = mesh % edgesOnVertex % array(iVertexDegree,iVertex)

             ! vector in great circle plane
             iCell1 = mesh % cellsOnEdge % array(1,iEdge)
             iCell2 = mesh % cellsOnEdge % array(2,iEdge)

             ! perform grid rotations
             call grid_rotation_forward(edgeRotated(1),               edgeRotated(2),               edgeRotated(3), &
                                        mesh % xEdge % array(iEdge),  mesh % yEdge % array(iEdge),  mesh % zEdge % array(iEdge))

             call grid_rotation_forward(cellRotated1(1),              cellRotated1(2),              cellRotated1(3), &
                                        mesh % xCell % array(iCell1), mesh % yCell % array(iCell1), mesh % zCell % array(iCell1))

             call grid_rotation_forward(cellRotated2(1),              cellRotated2(2),              cellRotated2(3), &
                                        mesh % xCell % array(iCell2), mesh % yCell % array(iCell2), mesh % zCell % array(iCell2))

             ! rotate to equator in new coords
             edgeEquator  = matmul(yRotationMatrix, matmul(zRotationMatrix, edgeRotated))
             cellEquator1 = matmul(yRotationMatrix, matmul(zRotationMatrix, cellRotated1))
             cellEquator2 = matmul(yRotationMatrix, matmul(zRotationMatrix, cellRotated2))

             ! vector joining vertices
             cellVector = cellEquator2 - cellEquator1

             ! form great circle plane normal vector
             call cross_product_3space(normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                       cellVector(1),        cellVector(2),        cellVector(3), &
                                       edgeEquator(1),       edgeEquator(2),       edgeEquator(3))

             if (iVertex == mesh % verticesOnEdge % array(1,iEdge)) then
                normalGreatCircle = -1.0_RKIND * normalGreatCircle
             endif

             ! normalize normal vector
             normalGreatCircleNorm = sqrt(normalGreatCircle(1)**2 + normalGreatCircle(2)**2 + normalGreatCircle(3)**2)
             normalGreatCircle(1) = normalGreatCircle(1) / normalGreatCircleNorm
             normalGreatCircle(2) = normalGreatCircle(2) / normalGreatCircleNorm
             normalGreatCircle(3) = normalGreatCircle(3) / normalGreatCircleNorm
             
             ! eastwards vector at edge
             eastwardsVector(1) = -edgeEquator(2)
             eastwardsVector(2) =  edgeEquator(1)
             eastwardsVector(3) = 0.0_RKIND
             eastwardsVector = eastwardsVector / &
                               sqrt(eastwardsVector(1)**2 + eastwardsVector(2)**2)
             
             call dot_product_3space(normalVectorTriangle(1,iVertexDegree,iVertex), &
                                     normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                     eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
             normalVectorTriangle(2,iVertexDegree,iVertex) = sign(1.0_RKIND,normalGreatCircle(3)) &
                                                           * sqrt(1.0_RKIND - max(min(normalVectorTriangle(1,iVertexDegree,iVertex),1.0_RKIND),-1.0_RKIND)**2)

          enddo ! iVertexDegree

       endif ! interiorVertex

    enddo ! iVertex

  end subroutine normal_vectors_spherical_triangle_metric

  !-------------------------------------------------------------

  subroutine normal_vectors_spherical_triangle(mesh, normalVectorTriangle, interiorVertex)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorTriangle

    integer, dimension(:), intent(in) :: &
         interiorVertex

    real(kind=RKIND), dimension(3) :: &
         normalGreatCircle, &
         eastwardsVector

    real(kind=RKIND) :: &
         xEdge,  yEdge,  zEdge,  &
         xCell1, yCell1, zCell1, &
         xCell2, yCell2, zCell2

    integer :: &
         iVertex, &
         iVertexDegree, &
         iEdge, &
         iCell1, &
         iCell2

    do iVertex = 1, mesh % nVertices

       normalVectorTriangle(:,:,iVertex) = 0.0_RKIND

       if (interiorVertex(iVertex) == 1) then

          do iVertexDegree = 1, mesh % vertexDegree
             
             iEdge = mesh % edgesOnVertex % array(iVertexDegree,iVertex)

             ! vector in great circle plane
             iCell1 = mesh % cellsOnEdge % array(1,iEdge)
             iCell2 = mesh % cellsOnEdge % array(2,iEdge)

             ! perform grid rotations
             call grid_rotation_forward(xEdge, yEdge, zEdge, &
                                        mesh % xEdge % array(iEdge), mesh % yEdge % array(iEdge), mesh % zEdge % array(iEdge))

             call grid_rotation_forward(xCell1, yCell1, zCell1, &
                                        mesh % xCell % array(iCell1), mesh % yCell % array(iCell1), mesh % zCell % array(iCell1))

             call grid_rotation_forward(xCell2, yCell2, zCell2, &
                                        mesh % xCell % array(iCell2), mesh % yCell % array(iCell2), mesh % zCell % array(iCell2))

             ! form great circle plane normal vector
             call cross_product_3space(normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                       xCell2 - xCell1,      yCell2 - yCell1,      zCell2 - zCell1,      &
                                       xEdge,                yEdge,                zEdge)

             if (iVertex == mesh % verticesOnEdge % array(2,iEdge)) then
                normalGreatCircle = -1.0_RKIND * normalGreatCircle
             endif
             
             ! eastwards vector at edge
             eastwardsVector(1) = -yEdge
             eastwardsVector(2) =  xEdge
             eastwardsVector(3) = 0.0_RKIND
             eastwardsVector = eastwardsVector / &
                               sqrt(eastwardsVector(1)**2 + eastwardsVector(2)**2)
             
             call dot_product_3space(normalVectorTriangle(1,iVertexDegree,iVertex), &
                                     normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                     eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
             normalVectorTriangle(2,iVertexDegree,iVertex) = sign(1.0_RKIND,normalGreatCircle(3))&
                                                           * sqrt(1.0_RKIND - min(normalVectorTriangle(1,iVertexDegree,iVertex),1.0_RKIND)**2)
             
          enddo ! iVertexDegree

       endif ! interiorVertex

    enddo ! iVertex

  end subroutine normal_vectors_spherical_triangle

  !-------------------------------------------------------------

  subroutine dot_product_3space(dot, x1, y1, z1, x2, y2, z2)

    real(kind=RKIND), intent(out) :: &
         dot

    real(kind=RKIND), intent(in) :: &
         x1, y1, z1, &
         x2, y2, z2

    dot = x1 * x2 + y1 * y2 + z1 * z2

  end subroutine dot_product_3space

  !-------------------------------------------------------------

  subroutine cross_product_3space(cp1, cp2, cp3, x1, y1, z1, x2, y2, z2)

    real(kind=RKIND), intent(out) :: &
         cp1, cp2, cp3

    real(kind=RKIND), intent(in) :: &
         x1, y1, z1, &
         x2, y2, z2

    cp1 = y1 * z2 - z1 * y2
    cp2 = z1 * x2 - x1 * z2
    cp3 = x1 * y2 - y1 * x2

  end subroutine cross_product_3space

  !-------------------------------------------------------------
  ! rotated lat-lon grid
  !-------------------------------------------------------------

  subroutine grid_rotation_forward(xp, yp, zp, x, y, z)

    ! rotate xyz coordinates from geographical grid to rotated grid with poles on real equator

    real(kind=RKIND), intent(out) :: &
         xp, & ! x position of point on rotated grid
         yp, & ! y position of point on rotated grid
         zp    ! z position of point on rotated grid

    real(kind=RKIND), intent(in) :: &
         x,  & ! x position of point on geographical grid
         y,  & ! y position of point on geographical grid
         z     ! z position of point on geographical grid

    call grid_rotation_forward_nochange(xp, yp, zp, x, y, z)

    !xp = -z
    !yp = y
    !zp = x

  end subroutine grid_rotation_forward

  !-------------------------------------------------------------

  subroutine grid_rotation_forward_nochange(xp, yp, zp, x, y, z)

    ! leave the grid unchanged - for testing

    real(kind=RKIND), intent(out) :: &
         xp, & ! x position of point on rotated grid
         yp, & ! y position of point on rotated grid
         zp    ! z position of point on rotated grid

    real(kind=RKIND), intent(in) :: &
         x,  & ! x position of point on geographical grid
         y,  & ! y position of point on geographical grid
         z     ! z position of point on geographical grid

    xp = x
    yp = y
    zp = z

  end subroutine grid_rotation_forward_nochange

  !-------------------------------------------------------------

  subroutine grid_rotation_backward(x, y, z, xp, yp, zp)

    ! rotate xyz coordinates from rotated grid with poles on real equator to geographical grid

    real(kind=RKIND), intent(out) :: &
         x,  & ! x position of point on geographical grid
         y,  & ! y position of point on geographical grid
         z     ! z position of point on geographical grid

    real(kind=RKIND), intent(in) :: &
         xp, & ! x position of point on rotated grid
         yp, & ! y position of point on rotated grid
         zp    ! z position of point on rotated grid

    call grid_rotation_backward_nochange(x, y, z, xp, yp, zp)

    !x = zp
    !y = yp
    !z = -xp

  end subroutine grid_rotation_backward

  !-------------------------------------------------------------

  subroutine grid_rotation_backward_nochange(x, y, z, xp, yp, zp)

    ! rotate xyz coordinates from rotated grid with poles on real equator to geographical grid

    real(kind=RKIND), intent(out) :: &
         x,  & ! x position of point on geographical grid
         y,  & ! y position of point on geographical grid
         z     ! z position of point on geographical grid

    real(kind=RKIND), intent(in) :: &
         xp, & ! x position of point on rotated grid
         yp, & ! y position of point on rotated grid
         zp    ! z position of point on rotated grid

    x = xp
    y = yp
    z = zp

  end subroutine grid_rotation_backward_nochange

  !-------------------------------------------------------------

  subroutine xyz_vector_rotation_forward(vxp, vyp, vzp, vx, vy, vz)

    ! rotate a xyz vector from geographical grid to rotated grid with poles on real equator

    real(kind=RKIND), intent(out) :: &
         vxp, & ! x component of velocity vector on rotated grid
         vyp, & ! y component of velocity vector on rotated grid
         vzp    ! z component of velocity vector on rotated grid

    real(kind=RKIND), intent(in) :: &
         vx,  & ! x component of velocity vector on geographical grid
         vy,  & ! y component of velocity vector on geographical grid
         vz     ! z component of velocity vector on geographical grid
    
    call xyz_vector_rotation_forward_nochange(vxp, vyp, vzp, vx, vy, vz) 

    !vxp = -vz
    !vyp =  vy
    !vzp =  vx

  end subroutine xyz_vector_rotation_forward

  !-------------------------------------------------------------

  subroutine xyz_vector_rotation_forward_nochange(vxp, vyp, vzp, vx, vy, vz)

    ! rotate a xyz vector from geographical grid to rotated grid with poles on real equator

    real(kind=RKIND), intent(out) :: &
         vxp, & ! x component of velocity vector on rotated grid
         vyp, & ! y component of velocity vector on rotated grid
         vzp    ! z component of velocity vector on rotated grid

    real(kind=RKIND), intent(in) :: &
         vx,  & ! x component of velocity vector on geographical grid
         vy,  & ! y component of velocity vector on geographical grid
         vz     ! z component of velocity vector on geographical grid
    
    vxp = vx
    vyp = vy
    vzp = vz

  end subroutine xyz_vector_rotation_forward_nochange

  !-------------------------------------------------------------

  subroutine xyz_vector_rotation_backward(vx, vy, vz, vxp, vyp, vzp)

    ! rotate xyz vector from rotated grid with poles on real equator to geographical grid

    real(kind=RKIND), intent(out) :: &
         vx,  & ! x component of velocity vector on geographical grid
         vy,  & ! y component of velocity vector on geographical grid
         vz     ! z component of velocity vector on geographical grid

    real(kind=RKIND), intent(in) :: &
         vxp, & ! x component of velocity vector on rotated grid
         vyp, & ! y component of velocity vector on rotated grid
         vzp    ! z component of velocity vector on rotated grid

    call xyz_vector_rotation_backward_nochange(vx, vy, vz, vxp, vyp, vzp)

    !vx =  vzp
    !vy =  vyp
    !vz = -vxp

  end subroutine xyz_vector_rotation_backward

  !-------------------------------------------------------------

  subroutine xyz_vector_rotation_backward_nochange(vx, vy, vz, vxp, vyp, vzp)

    ! rotate xyz vector from rotated grid with poles on real equator to geographical grid

    real(kind=RKIND), intent(out) :: &
         vx,  & ! x component of velocity vector on geographical grid
         vy,  & ! y component of velocity vector on geographical grid
         vz     ! z component of velocity vector on geographical grid

    real(kind=RKIND), intent(in) :: &
         vxp, & ! x component of velocity vector on rotated grid
         vyp, & ! y component of velocity vector on rotated grid
         vzp    ! z component of velocity vector on rotated grid

    vx = vxp
    vy = vyp
    vz = vzp

  end subroutine xyz_vector_rotation_backward_nochange

  !-------------------------------------------------------------

  subroutine latlon_vector_rotation_forward(up, vp, u, v, lat, lon, x, y, z, r)

    ! rotate a latlon vector from geographical grid to rotated grid with poles on real equator

    real(kind=RKIND), intent(out) :: &
         up,   & ! u component of velocity in rotated grid
         vp      ! v component of velocity in rotated grid

    real(kind=RKIND), intent(in) :: &
         u,    & ! u component of velocity on geographical grid
         v,    & ! v component of velocity on geographical grid
         lat,  & ! latitude of point on geographical grid
         lon,  & ! longitude of point on geographical grid
         x,    & ! x position of point on geographical grid
         y,    & ! y position of point on geographical grid
         z,    & ! z position of point on geographical grid
         r       ! radius of the earth 

    real(kind=RKIND) :: &
         xp,   & ! x position of point on rotated grid
         yp,   & ! y position of point on rotated grid
         zp,   & ! z position of point on rotated grid
         latp, & ! latitude of point on rotated grid
         lonp, & ! longitude of point on rotated grid
         vx,   & ! x component of velocity vector on geographical grid
         vy,   & ! y component of velocity vector on geographical grid
         vz,   & ! z component of velocity vector on geographical grid
         vxp,  & ! x component of velocity vector on rotated grid
         vyp,  & ! y component of velocity vector on rotated grid
         vzp     ! z component of velocity vector on rotated grid

    ! perform rotation of the point from geographical grid to rotated grid
    call grid_rotation_forward(xp, yp, zp, x, y, z)

    ! calculate latitude and longitude of the point in the rotated grid
    call latlon_from_xyz(latp, lonp, xp, yp, zp, r)

    ! convert lat lon vector to xyz vector on geographical grid
    call latlon_vector_to_xyz_vector(vx, vy, vz, u, v, lat, lon)

    ! perform rotation of geographical xyz vector to rotated grid
    call xyz_vector_rotation_forward(vxp, vyp, vzp, vx, vy, vz)

    ! convert xyz vector to lat lon vector on rotated grid
    call xyz_vector_to_latlon_vector(up, vp, vxp, vyp, vzp, latp, lonp)

  end subroutine latlon_vector_rotation_forward

  !-------------------------------------------------------------

  subroutine latlon_vector_rotation_backward(u, v, up, vp, lat, lon, x, y, z, r)

    ! rotate latlon vector from rotated grid with poles on real equator to geographical grid

    real(kind=RKIND), intent(out) :: &
         u,    & ! u component of velocity on geographical grid
         v       ! v component of velocity on geographical grid

    real(kind=RKIND), intent(in) :: &
         up,   & ! u component of velocity in rotated grid
         vp,   & ! v component of velocity in rotated grid
         lat,  & ! latitude of point on geographical grid
         lon,  & ! longitude of point on geographical grid
         x,    & ! x position of point on geographical grid
         y,    & ! y position of point on geographical grid
         z,    & ! z position of point on geographical grid
         r       ! radius of the earth 

    real(kind=RKIND) :: &
         xp,   & ! x position of point on rotated grid
         yp,   & ! y position of point on rotated grid
         zp,   & ! z position of point on rotated grid
         latp, & ! latitude of point on rotated grid
         lonp, & ! longitude of point on rotated grid
         vx,   & ! x component of velocity vector on geographical grid
         vy,   & ! y component of velocity vector on geographical grid
         vz,   & ! z component of velocity vector on geographical grid
         vxp,  & ! x component of velocity vector on rotated grid
         vyp,  & ! y component of velocity vector on rotated grid
         vzp     ! z component of velocity vector on rotated grid

    ! perform rotation of the point from geographical grid to rotated grid
    call grid_rotation_forward(xp, yp, zp, x, y, z)

    ! calculate latitude and longitude of the point in the rotated grid
    call latlon_from_xyz(latp, lonp, xp, yp, zp, r)

    ! convert lat lon vector to xyz vector on rotated grid
    call latlon_vector_to_xyz_vector(vxp, vyp, vzp, up, vp, latp, lonp)

    ! perform rotation of rotated xyz vector to geographical grid
    call xyz_vector_rotation_backward(vx, vy, vz, vxp, vyp, vzp)

    ! convert xyz vector to lat lon vector on geographical grid
    call xyz_vector_to_latlon_vector(u, v, vx, vy, vz, lat, lon)
     
  end subroutine latlon_vector_rotation_backward

  !-------------------------------------------------------------

  subroutine latlon_vector_to_xyz_vector(vx, vy, vz, u, v, lat, lon)

    ! convert a latlon vector to a xyz vector

    real(kind=RKIND), intent(out) :: &
         vx,  & ! x component of velocity vector
         vy,  & ! y component of velocity vector
         vz     ! z component of velocity vector

    real(kind=RKIND), intent(in) :: &
         u,   & ! u component of velocity
         v,   & ! v component of velocity
         lat, & ! latitude of point 
         lon    ! longitude of point 

    vx = (-u) * sin(lon) - v * sin(lat) * cos(lon)
    vy =   u  * cos(lon) - v * sin(lat) * sin(lon)
    vz =                   v * cos(lat)

  end subroutine latlon_vector_to_xyz_vector

  !-------------------------------------------------------------

  subroutine xyz_vector_to_latlon_vector(u, v, vx, vy, vz, lat, lon)

    ! convert a xyz vector vector to a latlon vector

    real(kind=RKIND), intent(out) :: &
         u,   & ! u component of velocity
         v      ! v component of velocity

    real(kind=RKIND), intent(in) :: &
         vx,  & ! x component of velocity vector
         vy,  & ! y component of velocity vector
         vz,  & ! z component of velocity vector
         lat, & ! latitude of point 
         lon    ! longitude of point 

    u = (-sin(lon)) * vx + &
          cos(lon)  * vy
    
    !write(*,*) 1, (-sin(lat)) * cos(lon), vx, (-sin(lat)) * cos(lon) * vx, (-sin(lat)) * cos(lon),  vy, (-sin(lat)) * cos(lon) * vy, cos(lat),  vz, cos(lat)  * vz

    v = (-sin(lat)) * cos(lon) * vx + &
        (-sin(lat)) * sin(lon) * vy + &
          cos(lat)  * vz

  end subroutine xyz_vector_to_latlon_vector

  !-------------------------------------------------------------

  subroutine latlon_from_xyz(lat, lon, x, y, z, r)

    ! given xyz coordinates determine the latitude and longitude

    real(kind=RKIND), intent(out) :: &
         lat, & ! latitude of point 
         lon    ! longitude of point 

    real(kind=RKIND), intent(in) :: &
         x,   & ! x position of point
         y,   & ! y position of point
         z,   & ! z position of point
         r      ! radius of earth

    lon = atan2(y, x)
    lat = asin(z/r)

  end subroutine latlon_from_xyz

  !-------------------------------------------------------------
  ! LU decomposition
  !-------------------------------------------------------------

  subroutine solve_linear_basis_system(leftMatrix, rightHandSide, solutionVector)

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         leftMatrix
    real(kind=RKIND), dimension(:), intent(in) :: &
         rightHandSide
    real(kind=RKIND), dimension(:), intent(out) :: &
         solutionVector

    real(kind=RKIND), allocatable, dimension(:,:) :: a
    real(kind=RKIND), allocatable, dimension(:) :: b
    integer, allocatable, dimension(:) :: indx
    real(kind=RKIND) :: d

    integer :: n

    n = size(rightHandSide)

    allocate(a(n,n))
    allocate(b(n))
    allocate(indx(n))

    a(:,:) = leftMatrix(:,:)
    b(:) = rightHandSide(:)

    call ludcmp(a,indx,d)

    call lubksb(a,indx,b)

    solutionVector(:) = b(:)

    deallocate(a)
    deallocate(b)
    deallocate(indx)

  end subroutine solve_linear_basis_system

  !-------------------------------------------------------------

  SUBROUTINE ludcmp(a,indx,d)

    IMPLICIT NONE
    REAL(kind=RKIND), DIMENSION(:,:), INTENT(INOUT) :: a
    INTEGER, DIMENSION(:), INTENT(OUT) :: indx
    REAL(kind=RKIND), INTENT(OUT) :: d
    REAL(kind=RKIND), DIMENSION(size(a,1)) :: vv
    REAL(kind=RKIND), PARAMETER :: TINY=1.0e-20_RKIND
    INTEGER :: j,n,imax
    n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
    d=1.0
    vv=maxval(abs(a),dim=2)
    if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')
    vv=1.0_RKIND/vv
    do j=1,n
       imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
       if (j /= imax) then
          call swap(a(imax,:),a(j,:))
          d=-d
          vv(imax)=vv(j)
       end if
       indx(j)=imax
       if (a(j,j) == 0.0) a(j,j)=TINY
       a(j+1:n,j)=a(j+1:n,j)/a(j,j)
       a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
    end do
  END SUBROUTINE ludcmp

  !-------------------------------------------------------------

  SUBROUTINE lubksb(a,indx,b)

    IMPLICIT NONE
    REAL(kind=RKIND), DIMENSION(:,:), INTENT(IN) :: a
    INTEGER, DIMENSION(:), INTENT(IN) :: indx
    REAL(kind=RKIND), DIMENSION(:), INTENT(INOUT) :: b
    INTEGER :: i,n,ii,ll
    REAL(kind=RKIND) :: summ
    n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
    ii=0
    do i=1,n
       ll=indx(i)
       summ=b(ll)
       b(ll)=b(i)
       if (ii /= 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
       else if (summ /= 0.0) then
          ii=i
       end if
       b(i)=summ
    end do
    do i=n,1,-1
       b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    end do
  END SUBROUTINE lubksb

  !-------------------------------------------------------------

  FUNCTION assert_eq(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq
    if (n1 == n2 .and. n2 == n3) then
       assert_eq=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq

  !-------------------------------------------------------------

  FUNCTION imaxloc(arr)
    REAL(kind=RKIND), DIMENSION(:), INTENT(IN) :: arr
    INTEGER :: imaxloc
    INTEGER, DIMENSION(1) :: imax
    imax=maxloc(arr(:))
    imaxloc=imax(1)
  END FUNCTION imaxloc

  !-------------------------------------------------------------

  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror

  !-------------------------------------------------------------

  FUNCTION outerprod(a,b)
    REAL(kind=RKIND), DIMENSION(:), INTENT(IN) :: a,b
    REAL(kind=RKIND), DIMENSION(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod

  !-------------------------------------------------------------

  SUBROUTINE swap(a,b)
    REAL(kind=RKIND), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(kind=RKIND), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap

  !-------------------------------------------------------------

end module cice_dynamics_shared
