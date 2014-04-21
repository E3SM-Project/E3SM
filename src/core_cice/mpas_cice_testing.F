module cice_testing

  use mpas_grid_types

  implicit none

  private
  public :: &
       divergence_stress_test_velocity_set, &
       divergence_stress_test_stress_set_hex, &
       divergence_stress_test_stress_set_tri, &
       divergence_stress_test_stress_set_weak, &
       init_square_test_case_hex, &
       init_spherical_test_case

  real(kind=RKIND), parameter :: &
       rhoi = 917.0_RKIND

  ! parameter in constant stress strain relation
  real(kind=RKIND), parameter, public :: lambda = 1.0_RKIND

contains

  !-------------------------------------------------------------
  ! Spherical test case
  !-------------------------------------------------------------
  
  subroutine init_spherical_test_case(&
       mesh, &
       iceAreaCell, &
       iceVolumeCell, &
       totalMassCell, &
       surfaceTemperature, &
       airTemperature, &
       uVelocity, &
       vVelocity, &
       interiorVertex)
    
    use mpas_constants, only: omega
    use cice_diagnostics, only: init_ijpop_from_ivertex

    type(MPAS_pool_type), intent(inout) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         iceAreaCell, &
         iceVolumeCell, &
         totalMassCell, &
         surfaceTemperature, &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:), intent(in) :: &
         airTemperature

    integer, dimension(:), intent(in) :: &
         interiorVertex

    real(kind=RKIND), parameter :: &
         iceHeightCell = 1.0_RKIND

    integer :: &
         iCell, &
         iVertex

    real(kind=RKIND), parameter :: &  
         pi = 3.14159265359_RKIND, &
         degreesToRadians = pi / 180.0_RKIND, &
         initialIceEdgeLatitudeNorthernHemisphere =  70.0_RKIND * degreesToRadians, &
         initialIceEdgeLatitudeSouthernHemisphere = -60.0_RKIND * degreesToRadians, &
         circle_radius = 1.0e6_RKIND

    real(kind=RKIND) :: &
         x, y, z

    character(len=200), parameter :: &
         test_type = "icecaps"
         !test_type = "circleofice"

    integer, pointer :: &
         nCells, &
         nVertices

    real(kind=RKIND), dimension(:), pointer :: &
         fVertex, &
         latVertex, &
         latCell, &
         xCell, &
         yCell, &
         zCell

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "fVertex", fVertex)
    call MPAS_pool_get_array(mesh, "latVertex", latVertex)
    call MPAS_pool_get_array(mesh, "latCell", latCell)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    call init_ijpop_from_ivertex(mesh)

    fVertex = 2.0_RKIND * omega * sin(latVertex)

    if (trim(test_type) == "icecaps") then

       do iCell = 1, nCells

          if (latCell(iCell) > initialIceEdgeLatitudeNorthernHemisphere .or. &
              latCell(iCell) < initialIceEdgeLatitudeSouthernHemisphere) then
          !if (latCell(iCell) < initialIceEdgeLatitudeSouthernHemisphere) then

             ! ice present
             iceAreaCell(iCell)        = 1.0_RKIND
             iceVolumeCell(iCell)      = iceAreaCell(iCell) * iceHeightCell
             totalMassCell(iCell)      = iceVolumeCell(iCell) * rhoi
             surfaceTemperature(iCell) = 0.0_RKIND

          else

             ! no ice
             iceAreaCell(iCell)   = 0.0_RKIND
             iceVolumeCell(iCell) = 0.0_RKIND
             totalMassCell(iCell) = 0.0_RKIND
             surfaceTemperature(iCell) = 0.0_RKIND

          endif

       end do ! iCell

       do iVertex = 1, nVertices

          if (interiorVertex(iVertex) == 1) then

             uVelocity(iVertex) = 0.0_RKIND
             vVelocity(iVertex) = 1.0_RKIND

          else

             uVelocity(iVertex) = 0.0_RKIND
             vVelocity(iVertex) = 0.0_RKIND

          endif

       enddo ! iVertex

    else if (trim(test_type) == "circleofice") then

       ! circle of ice at equator
       do iCell = 1, nCells

          x = xCell(iCell)
          y = yCell(iCell)
          z = zCell(iCell)

          if (sqrt(x**2+z**2) < circle_radius .and. y > 0.0_RKIND) then

             iceAreaCell(iCell)        = 1.0_RKIND
             iceVolumeCell(iCell)      = iceAreaCell(iCell) * iceHeightCell
             totalMassCell(iCell)      = iceVolumeCell(iCell) * rhoi
             surfaceTemperature(iCell) = 0.0_RKIND

          else

             ! no ice
             iceAreaCell(iCell)   = 0.0_RKIND
             iceVolumeCell(iCell) = 0.0_RKIND
             totalMassCell(iCell) = 0.0_RKIND
             surfaceTemperature(iCell) = 0.0_RKIND

          endif

       enddo ! iCell

       do iVertex = 1, nVertices

          if (interiorVertex(iVertex) == 1) then

             uVelocity(iVertex) = 1.0_RKIND
             vVelocity(iVertex) = 0.0_RKIND

          else

             uVelocity(iVertex) = 0.0_RKIND
             vVelocity(iVertex) = 0.0_RKIND

          endif

       enddo ! iVertex

    endif ! test_type

  end subroutine init_spherical_test_case

  !-------------------------------------------------------------
  ! Square test case
  !-------------------------------------------------------------

  subroutine init_square_test_case_hex(mesh, tracers, icestate, hexdyn, hexfor, boundary, atmos_forcing)

    type(MPAS_pool_type), pointer, intent(inout) :: mesh
    type(MPAS_pool_type), pointer :: tracers
    type(MPAS_pool_type), pointer :: icestate
    type(MPAS_pool_type), pointer :: hexdyn
    type(MPAS_pool_type), pointer :: hexfor
    type(MPAS_pool_type), pointer :: boundary
    type(MPAS_pool_type), pointer :: atmos_forcing

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         xVertex, &
         yVertex, &
         uAirVelocity, &
         vAirVelocity, &
         uOceanVelocity, &
         vOceanVelocity

    ! init variables
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)

    call MPAS_pool_get_array(atmos_forcing, "uAirVelocity", uAirVelocity)
    call MPAS_pool_get_array(atmos_forcing, "vAirVelocity", vAirVelocity)

    call MPAS_pool_get_array(hexfor, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(hexfor, "vOceanVelocity", vOceanVelocity)

    call square_test_correct_positions(mesh)

    call init_atmos_velocity(uAirVelocity, &
                             vAirVelocity, &
                             xCell,        &
                             yCell,        &
                             0.0_RKIND)

    call init_ocean_velocity(uOceanVelocity, &
                             vOceanVelocity, &
                             xVertex,        &
                             yVertex)

    call init_ice_state(mesh, tracers, icestate, hexdyn, hexfor, boundary)

  end subroutine init_square_test_case_hex

  !-------------------------------------------------------------

  subroutine init_ocean_velocity(uOceanVelocity, vOceanVelocity, x, y)

    real(kind=RKIND), dimension(:), intent(out) :: &
         uOceanVelocity, &
         vOceanVelocity

    real(kind=RKIND), dimension(:), intent(in) :: &
         x, y

    real(kind=RKIND), parameter :: a = 0.1_RKIND  
   
    real(kind=RKIND), parameter :: &
         Lx = 1.28e6_RKIND, &
         Ly = 1.28e6_RKIND

    integer :: iPoint, nPoints

    nPoints = size(uOceanVelocity,1)

    do iPoint = 1, nPoints
       
       uOceanVelocity(iPoint) =  a * ((2.0_RKIND * y(iPoint) - Ly) / Ly)

       vOceanVelocity(iPoint) = -a * ((2.0_RKIND * x(iPoint) - Lx) / Lx)
       
    enddo ! iCell

  end subroutine init_ocean_velocity

  !-------------------------------------------------------------

  subroutine init_atmos_velocity(uAirVelocity, vAirVelocity, xin, yin, time)

    real(kind=RKIND), dimension(:), intent(out) :: &
         uAirVelocity, &
         vAirVelocity

    real(kind=RKIND), dimension(:), intent(in) :: &
         xin, yin

    real(kind=RKIND), intent(in) :: &
         time

    real(kind=RKIND), parameter :: pi = &
         3.14159265358979323846_RKIND
         !3.14159265359_RKIND

    real(kind=RKIND), parameter :: a = 5.0_RKIND
    real(kind=RKIND), parameter :: b = 3.0_RKIND

    real(kind=RKIND), parameter :: theta = 4.0_RKIND * 24.0_RKIND * 3600.0_RKIND

    real(kind=RKIND), parameter :: &
         Lx = 1.28e6_RKIND, &
         Ly = 1.28e6_RKIND

    real(kind=RKIND) :: &
         x, y, &
         xmin, xmax, ymin, ymax

    integer :: iPoint, nPoints

    xmin = minval(xin)
    xmax = maxval(xin)
    ymin = minval(yin)
    ymax = maxval(yin)

    nPoints = size(uAirVelocity,1)

    do iPoint = 1, nPoints

       x = xin(iPoint)
       y = yin(iPoint)

       uAirVelocity(iPoint) = &!a * (y / Ly)
            a + (sin((2.0_RKIND * pi * time) / theta) - b) * sin(2.0_RKIND * pi * (x / Lx)) * sin(pi * (y / Ly))
       vAirVelocity(iPoint) = &!0.0_RKIND!&
            a + (sin((2.0_RKIND * pi * time) / theta) - b) * sin(2.0_RKIND * pi * (y / Ly)) * sin(pi * (x / Lx))

       !uAirVelocity(iPoint) = a
       !vAirVelocity(iPoint) = a

       !write(*,*) iPoint, x, y, uAirVelocity(iPoint), vAirVelocity(iPoint), time

    enddo ! iPoint
!stop
  end subroutine init_atmos_velocity

  !-------------------------------------------------------------

  subroutine init_ice_state(mesh, tracers, icestate, hexdyn, hexfor, boundary)

    type(MPAS_pool_type), pointer, intent(in) :: mesh
    type(MPAS_pool_type), pointer :: tracers
    type(MPAS_pool_type), pointer :: icestate
    type(MPAS_pool_type), pointer :: hexdyn
    type(MPAS_pool_type), pointer :: hexfor
    type(MPAS_pool_type), pointer :: boundary

    real(kind=RKIND) :: iceThickness

    real(kind=RKIND) :: x    

    real(kind=RKIND), parameter :: &
         Lx = 1.28e6_RKIND

    integer :: &
         iCell, &
         iVertex

    integer, pointer :: &
         nCells, &
         nVertices

    integer, dimension(:), pointer :: &
         interiorVertex

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         iceAreaCell, &
         iceVolumeCell, &
         totalMassCell, &
         uVelocity, &
         vVelocity, &
         uOceanVelocity, &
         vOceanVelocity
    
    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "xCell", xCell)

    call MPAS_pool_get_array(tracers, "iceAreaCell", iceAreaCell)
    call MPAS_pool_get_array(tracers, "iceVolumeCell", iceVolumeCell)
    call MPAS_pool_get_array(icestate, "totalMassCell", totalMassCell)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call MPAS_pool_get_array(hexdyn, "uVelocity", uVelocity)
    call MPAS_pool_get_array(hexdyn, "vVelocity", vVelocity)

    call MPAS_pool_get_array(hexfor, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(hexfor, "vOceanVelocity", vOceanVelocity)

    iceThickness = 2.0_RKIND

    do iCell = 1, nCells

       x = xCell(iCell)

       iceAreaCell(iCell) = &!0.95_RKIND!&
            max(min(x / Lx, 1.0_RKIND), 0.0_RKIND)

       iceVolumeCell(iCell) = &
            iceThickness * iceAreaCell(iCell)

       totalMassCell(iCell) = iceVolumeCell(iCell) * rhoi
          
    enddo ! iCell

    do iVertex = 1, nVertices
    
       if (interiorVertex(iVertex) == 1) then

          uVelocity(iVertex) = uOceanVelocity(iVertex)
          vVelocity(iVertex) = vOceanVelocity(iVertex)

       else

          uVelocity(iVertex) = 0.0_RKIND
          vVelocity(iVertex) = 0.0_RKIND

       endif

    enddo ! iVertex

  end subroutine init_ice_state

  !-------------------------------------------------------------

  subroutine square_test_correct_positions(mesh)

    type(MPAS_pool_type), intent(inout) :: mesh

    ! periodic quad - use with ocean
    real(kind=RKIND), parameter :: &
         dx = -16000.0_RKIND, &
         dy = -16000.0_RKIND

    ! periodic hex - use with ocean82x94.nc
    !real(kind=RKIND), parameter :: &
    !     dx = -16000.0_RKIND * (5.0_RKIND / 4.0_RKIND), &
    !     dy = -16000.0_RKIND * (sqrt(3.0_RKIND)/2.0_RKIND) - 16000.0_RKIND * (2.0_RKIND / sqrt(3.0_RKIND)) * 0.25_RKIND

    ! longitude/latitude parameters - Barrow AK
    real(kind=RKIND), parameter :: &
         longitudeSquare = -156.5_RKIND, &
         latitudeSquare  = 71.35_RKIND, &
         radius          = 6.37e6_RKIND, &
         pi              = 3.14159265358979323846_RKIND, &
         rad_to_deg      = 180.0_RKIND / pi, &
         distanceToAngle = rad_to_deg / radius, &
         omega           = 7.292e-5_RKIND
    
    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         lonCell, &
         latCell, &
         xVertex, &
         yVertex, &
         lonVertex, &
         latVertex, &
         xEdge, &
         yEdge, &
         lonEdge, &
         latEdge, &
         fCell, &
         fVertex, &
         fEdge

    ! init variables
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "lonCell", lonCell)
    call MPAS_pool_get_array(mesh, "latCell", latCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "lonVertex", lonVertex)
    call MPAS_pool_get_array(mesh, "latVertex", latVertex)
    call MPAS_pool_get_array(mesh, "xEdge", xEdge)
    call MPAS_pool_get_array(mesh, "yEdge", yEdge)
    call MPAS_pool_get_array(mesh, "lonEdge", lonEdge)
    call MPAS_pool_get_array(mesh, "latEdge", latEdge)
    call MPAS_pool_get_array(mesh, "fCell", fCell)
    call MPAS_pool_get_array(mesh, "fVertex", fVertex)
    call MPAS_pool_get_array(mesh, "fEdge", fEdge)

    ! Cell
    xCell = xCell + dx
    yCell = yCell + dy

    lonCell = xCell * distanceToAngle + longitudeSquare
    latCell = yCell * distanceToAngle + latitudeSquare

    lonCell = lonCell / rad_to_deg
    latCell = latCell / rad_to_deg

    ! Vertex
    xVertex = xVertex + dx
    yVertex = yVertex + dy

    lonVertex = xVertex * distanceToAngle + longitudeSquare
    latVertex = yVertex * distanceToAngle + latitudeSquare

    lonVertex = lonVertex / rad_to_deg
    latVertex = latVertex / rad_to_deg

    ! Edge
    xEdge = xEdge + dx
    yEdge = yEdge + dy

    lonEdge = xEdge * distanceToAngle + longitudeSquare
    latEdge = yEdge * distanceToAngle + latitudeSquare

    lonEdge = lonEdge / rad_to_deg
    latEdge = latEdge / rad_to_deg

    ! fvalues
    fCell   = 2.0_RKIND * omega * sin(latCell)
    fVertex = 2.0_RKIND * omega * sin(latVertex)
    fEdge   = 2.0_RKIND * omega * sin(latEdge)

  end subroutine square_test_correct_positions

  !-------------------------------------------------------------
  ! stress divergence operator test velocities
  !-------------------------------------------------------------

  subroutine divergence_stress_test_velocity_set(uVelocity, vVelocity, x, y, type)
    
    real(kind=RKIND), dimension(:), intent(out) :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:), intent(in) :: &
         x, y
    
    character(len=*), intent(in) :: &
         type

    integer :: iPoint, nPoints

    real(kind=RKIND), parameter :: &
         velocityConstantU = 112.87654_RKIND, &
         velocityConstantV = -34.5678_RKIND, &
         velocityScale = 1.0_RKIND

    nPoints = size(uVelocity,1)

    select case (type)
    case ("zero")
       write(*,*) "zero velocities"

       ! zero velocities
       do iPoint = 1, nPoints
          uVelocity(iPoint) = 0.0_RKIND
          vVelocity(iPoint) = 0.0_RKIND
       enddo ! iPoint

    case ("constant")
       write(*,*) "constant velocities"

       ! constant velocities
       do iPoint = 1, nPoints
          uVelocity(iPoint) = velocityConstantU
          vVelocity(iPoint) = velocityConstantV
       enddo ! iPoint

    case ("linearx")
       write(*,*) "linearx velocities"

       ! constant velocities
       do iPoint = 1, nPoints
          uVelocity(iPoint) = x(iPoint)
          vVelocity(iPoint) = 0.0_RKIND
       enddo ! iPoint

    case ("lineary")
       write(*,*) "lineary velocities"

       ! constant velocities
       do iPoint = 1, nPoints
          uVelocity(iPoint) = 0.0_RKIND
          vVelocity(iPoint) = y(iPoint)
       enddo ! iPoint

    case ("constantsig12")
       write(*,*) "constant sigma_12"

       ! constant velocities
       do iPoint = 1, nPoints
          uVelocity(iPoint) = y(iPoint)
          vVelocity(iPoint) = x(iPoint)
       enddo ! iPoint

    case ("div1")
       write(*,*) "div1 velocities"

       ! velocity field that gives constant divergence
       do iPoint = 1, nPoints
          uVelocity(iPoint) = x(iPoint)**2 / (4.0_RKIND * lambda) + (x(iPoint) * y(iPoint)) / lambda
          vVelocity(iPoint) = y(iPoint)**2 / (4.0_RKIND * lambda) + (x(iPoint) * y(iPoint)) / lambda
       enddo ! iPoint

    case ("divx1")
       write(*,*) "divx1 velocities"

       ! velocity field that gives constant divergence
       do iPoint = 1, nPoints
          uVelocity(iPoint) = 0.5_RKIND * x(iPoint)**2
          vVelocity(iPoint) = 0.0_RKIND
       enddo ! iPoint

    case ("divy1")
       write(*,*) "divy1 velocities"

       ! velocity field that gives constant divergence
       do iPoint = 1, nPoints
          uVelocity(iPoint) = 0.0_RKIND
          vVelocity(iPoint) = 0.5_RKIND * y(iPoint)**2
       enddo ! iPoint

    case ("s12")
       write(*,*) "s12 velocities"

       ! velocity field that gives constant divergence
       do iPoint = 1, nPoints
          uVelocity(iPoint) = 0.5_RKIND * y(iPoint)**2
          vVelocity(iPoint) = 0.5_RKIND * x(iPoint)**2
       enddo ! iPoint

    end select

  end subroutine divergence_stress_test_velocity_set

  !-------------------------------------------------------------

  subroutine divergence_stress_test_stress_set_hex(mesh, stress1, stress2, stress12)

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         stress1, &
         stress2, &
         stress12

    real(kind=RKIND) :: &
         xpy, &
         x, y

    integer :: &
         iCell, &
         iVertexOnCell, &
         iVertex

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)

    do iCell = 1, nCells

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          iVertex = verticesOnCell(iVertexOnCell,iCell)

          x = xVertex(iVertex)
          y = yVertex(iVertex)

          xpy = x + y

          ! divu = 1 ; divv = 1
          stress1(iVertexOnCell,iCell)  =  1.5_RKIND * xpy
          stress2(iVertexOnCell,iCell)  = -0.5_RKIND * xpy
          stress12(iVertexOnCell,iCell) =  0.5_RKIND * xpy

          ! divu = 1 ; divv = 0
          !stress1(iVertexOnCell,iCell)  =  0.5_RKIND * x
          !stress2(iVertexOnCell,iCell)  =  0.5_RKIND * x
          !stress12(iVertexOnCell,iCell) =  0.5_RKIND * y

          ! divu = 0 ; divv = 1
          !stress1(iVertexOnCell,iCell)  =  0.5_RKIND * y
          !stress2(iVertexOnCell,iCell)  = -0.5_RKIND * y
          !stress12(iVertexOnCell,iCell) =  0.5_RKIND * x

          ! others
          !stress1(iVertexOnCell,iCell)  =  0.5_RKIND * x
          !stress2(iVertexOnCell,iCell)  =  0.5_RKIND * x
          !stress12(iVertexOnCell,iCell) =  0.0_RKIND

          !stress1(iVertexOnCell,iCell)  =  0.0_RKIND
          !stress2(iVertexOnCell,iCell)  =  0.0_RKIND
          !stress12(iVertexOnCell,iCell) =  xpy

          !stress1(iVertexOnCell,iCell)  =  y
          !stress2(iVertexOnCell,iCell)  =  -y
          !stress12(iVertexOnCell,iCell) =  0.0_RKIND

          stress1(iVertexOnCell,iCell)  =  0.5_RKIND * x + y
          stress2(iVertexOnCell,iCell)  =  0.5_RKIND * y + x
          stress12(iVertexOnCell,iCell) =  0.5_RKIND * (x + y)

          stress1(iVertexOnCell,iCell)  =  x
          stress2(iVertexOnCell,iCell)  =  0.0_RKIND
          stress12(iVertexOnCell,iCell) =  0.0_RKIND


       enddo ! iVertexOnCell

    enddo ! iCell

  end subroutine divergence_stress_test_stress_set_hex

  !-------------------------------------------------------------

  subroutine divergence_stress_test_stress_set_tri(mesh, stress1, stress2, stress12)

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         stress1, &
         stress2, &
         stress12

    real(kind=RKIND) :: &
         xpy, &
         x, y

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell

    integer, pointer :: &
         nVertices, &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell

    ! init variables
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)

    do iVertex = 1, nVertices

       do iVertexDegree = 1, vertexDegree

          iCell = cellsOnVertex(iVertexDegree,iVertex)

          x = xCell(iCell)
          y = yCell(iCell)

          xpy = x + y

          ! divu = 1 ; divv = 1
          !stress1(iVertexDegree,iVertex)  =  1.5_RKIND * xpy
          !stress2(iVertexDegree,iVertex)  = -0.5_RKIND * xpy
          !stress12(iVertexDegree,iVertex) =  0.5_RKIND * xpy

          ! divu = 1 ; divv = 0
          !stress1(iVertexDegree,iVertex)  =  0.5_RKIND * x
          !stress2(iVertexDegree,iVertex)  =  0.5_RKIND * x
          !stress12(iVertexDegree,iVertex) =  0.5_RKIND * y

          ! divu = 0 ; divv = 1
          !stress1(iVertexDegree,iVertex)  =  0.5_RKIND * y
          !stress2(iVertexDegree,iVertex)  = -0.5_RKIND * y
          !stress12(iVertexDegree,iVertex) =  0.5_RKIND * x

          ! others
          !stress1(iVertexDegree,iVertex)  =  0.5_RKIND * x
          !stress2(iVertexDegree,iVertex)  =  0.5_RKIND * x
          !stress12(iVertexDegree,iVertex) =  0.0_RKIND

          !stress1(iVertexDegree,iVertex)  =  0.0_RKIND
          !stress2(iVertexDegree,iVertex)  =  0.0_RKIND
          !stress12(iVertexDegree,iVertex) =  xpy

          !stress1(iVertexDegree,iVertex)  =  y
          !stress2(iVertexDegree,iVertex)  =  -y
          !stress12(iVertexDegree,iVertex) =  0.0_RKIND

          !stress1(iVertexDegree,iVertex)  =  0.5_RKIND * x + y
          !stress2(iVertexDegree,iVertex)  =  0.5_RKIND * y + x
          !stress12(iVertexDegree,iVertex) =  0.5_RKIND * (x + y)

          stress1(iVertexDegree,iVertex)  =  x
          stress2(iVertexDegree,iVertex)  =  0.0_RKIND
          stress12(iVertexDegree,iVertex) =  0.0_RKIND


       enddo ! iVertexOnCell

    enddo ! iCell

  end subroutine divergence_stress_test_stress_set_tri

  !-------------------------------------------------------------

  subroutine divergence_stress_test_stress_set_weak(mesh, stress11, stress22, stress12)

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(inout) :: &
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND) :: &
         xpy, &
         x, y

    integer :: &
         iCell

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell
    
    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)

    do iCell = 1, nCells

       x = xCell(iCell)
       y = yCell(iCell)

       xpy = x + y

       ! divu = 1 ; divv = 0
       !stress11(iCell) =  x
       !stress22(iCell) =  0.0_RKIND
       !stress12(iCell) =  0.0_RKIND

       ! divu = 0 ; divv = 1
       !stress11(iCell) =  0.0_RKIND
       !stress22(iCell) =  y
       !stress12(iCell) =  0.0_RKIND

       ! divu = 1 ; divv = 1
       !stress11(iCell) =  x
       !stress22(iCell) =  y
       !stress12(iCell) =  0.0_RKIND

       ! divu = 1 ; divv = 1
       !stress11(iCell) =  0.5_RKIND * xpy
       !stress22(iCell) =  0.5_RKIND * xpy
       !stress12(iCell) =  0.5_RKIND * xpy

       ! divu = 1 ; divv = 1
       stress11(iCell) =  100.0_RKIND
       stress22(iCell) =  -1000.0_RKIND
       stress12(iCell) =  1.0_RKIND * xpy

    enddo ! iCell

  end subroutine divergence_stress_test_stress_set_weak

end module cice_testing

