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
       spherical_test_strain, &
       spherical_test_divergence_stress, &
       init_spherical_test_case

  real(kind=RKIND), parameter :: &
       pi = 3.14159265359_RKIND, &
       deg_to_rad = pi / 180.0_RKIND, &
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

    type(mesh_type), intent(inout) :: mesh

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
         initialIceEdgeLatitudeSouthernHemisphere = -60.0_RKIND * degreesToRadians

    call init_ijpop_from_ivertex(mesh)

    mesh % fVertex % array = 2.0_RKIND * omega * sin(mesh % latVertex % array)

    do iCell = 1, mesh % nCells

       !if (mesh % latCell % array(iCell) > initialIceEdgeLatitudeNorthernHemisphere .or. &
       !    mesh % latCell % array(iCell) < initialIceEdgeLatitudeSouthernHemisphere) then
       if (mesh % latCell % array(iCell) < initialIceEdgeLatitudeSouthernHemisphere) then

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

    do iVertex = 1, mesh % nVertices

       if (interiorVertex(iVertex) == 1) then

          uVelocity(iVertex) = 0.0_RKIND
          vVelocity(iVertex) = 1.0_RKIND
          
       else

          uVelocity(iVertex) = 0.0_RKIND
          vVelocity(iVertex) = 0.0_RKIND

       endif

    enddo ! iVertex
    
  end subroutine init_spherical_test_case

  !-------------------------------------------------------------
  ! Square test case
  !-------------------------------------------------------------

  subroutine init_square_test_case_hex(mesh, icestate, hexdyn, hexfor, boundary, atmos_forcing)

    type(mesh_type), intent(inout) :: mesh

    type(icestate_type),      pointer :: icestate
    type(hexdyn_type),        pointer :: hexdyn
    type(hexfor_type),        pointer :: hexfor
    type(boundary_type),      pointer :: boundary
    type(atmos_forcing_type), pointer :: atmos_forcing

    call square_test_correct_positions(mesh)

    call init_atmos_velocity(atmos_forcing % uAirVelocity % array, &
                             atmos_forcing % vAirVelocity % array, &
                             mesh % xCell % array,                 &
                             mesh % yCell % array,                 &
                             0.0_RKIND)

    call init_ocean_velocity(hexfor % uOceanVelocity % array, &
                             hexfor % vOceanVelocity % array, &
                             mesh % xVertex % array,          &
                             mesh % yVertex % array)

    call init_ice_state(mesh, icestate, hexdyn, hexfor, boundary)

  end subroutine init_square_test_case_hex

  !-------------------------------------------------------------

  !subroutine init_square_test_case_tri(mesh, icestate, forcing)

  !  type(mesh_type), intent(in) :: mesh

  !  type(icestate_type), pointer :: icestate
  !  type(forcing_type),  pointer :: forcing

  !  call init_atmos_velocity(forcing % uAirVelocity % array, &
  !                           forcing % vAirVelocity % array, &
  !                           mesh % xCell % array,           &
  !                           mesh % yCell % array,           &
  !                           0.0_RKIND)

  !  call init_ocean_velocity(forcing % uOceanVelocity % array, &
  !                           forcing % vOceanVelocity % array, &
  !                           mesh % xCell % array,             &
  !                           mesh % yCell % array)

  !  !call init_ice_state(mesh, icestate)

  !end subroutine init_square_test_case_tri

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

  subroutine init_ice_state(mesh, icestate, hexdyn, hexfor, boundary)

    type(mesh_type), intent(in) :: mesh

    type(icestate_type), pointer :: icestate
    type(hexdyn_type),   pointer :: hexdyn
    type(hexfor_type),   pointer :: hexfor
    type(boundary_type), pointer :: boundary

    real(kind=RKIND) :: iceThickness

    real(kind=RKIND) :: x    

    real(kind=RKIND), parameter :: &
         Lx = 1.28e6_RKIND

    integer :: &
         iCell, &
         iVertex

    iceThickness = 2.0_RKIND

    do iCell = 1, mesh % nCells

       x = mesh % xCell % array(iCell)

       icestate % tracers % array(icestate % index_iceAreaCell, 1, iCell) = &!0.95_RKIND!&
            max(min(x / Lx, 1.0_RKIND), 0.0_RKIND)

       icestate % tracers % array(icestate % index_iceVolumeCell, 1, :) = &
            iceThickness * icestate % tracers % array(icestate % index_iceAreaCell, 1, iCell)

       icestate % totalMassCell % array(iCell) = icestate % tracers % array(icestate % index_iceVolumeCell, 1, iCell) * rhoi
          
    enddo ! iCell

    do iVertex = 1, mesh % nVertices
    
       if (boundary % interiorVertex % array(iVertex) == 1) then

          hexdyn % uVelocity % array(iVertex) = hexfor % uOceanVelocity % array(iVertex)
          hexdyn % vVelocity % array(iVertex) = hexfor % vOceanVelocity % array(iVertex)

       else

          hexdyn % uVelocity % array(iVertex) = 0.0_RKIND
          hexdyn % vVelocity % array(iVertex) = 0.0_RKIND

       endif

    enddo ! iVertex

  end subroutine init_ice_state

  !-------------------------------------------------------------

  subroutine square_test_correct_positions(mesh)

    type(mesh_type), intent(inout) :: mesh

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

    ! Cell
    mesh % xCell % array = mesh % xCell % array + dx
    mesh % yCell % array = mesh % yCell % array + dy

    mesh % lonCell % array = mesh % xCell % array * distanceToAngle + longitudeSquare
    mesh % latCell % array = mesh % yCell % array * distanceToAngle + latitudeSquare

    mesh % lonCell % array = mesh % lonCell % array / rad_to_deg
    mesh % lonCell % array = mesh % lonCell % array / rad_to_deg

    ! Vertex
    mesh % xVertex % array = mesh % xVertex % array + dx
    mesh % yVertex % array = mesh % yVertex % array + dy

    mesh % lonVertex % array = mesh % xVertex % array * distanceToAngle + longitudeSquare
    mesh % latVertex % array = mesh % yVertex % array * distanceToAngle + latitudeSquare

    mesh % lonVertex % array = mesh % lonVertex % array / rad_to_deg
    mesh % latVertex % array = mesh % latVertex % array / rad_to_deg

    ! Edge
    mesh % xEdge % array = mesh % xEdge % array + dx
    mesh % yEdge % array = mesh % yEdge % array + dy

    mesh % lonEdge % array = mesh % xEdge % array * distanceToAngle + longitudeSquare
    mesh % latEdge % array = mesh % yEdge % array * distanceToAngle + latitudeSquare

    mesh % lonEdge % array = mesh % lonEdge % array / rad_to_deg
    mesh % latEdge % array = mesh % latEdge % array / rad_to_deg

    ! fvalues
    mesh % fCell % array   = 2.0_RKIND * omega * sin(mesh % latCell % array)
    mesh % fVertex % array = 2.0_RKIND * omega * sin(mesh % latVertex % array)
    mesh % fEdge % array   = 2.0_RKIND * omega * sin(mesh % latEdge % array)

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

    type(mesh_type), intent(in) :: mesh

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

    do iCell = 1, mesh % nCells

       do iVertexOnCell = 1, mesh % nEdgesOnCell % array(iCell)

          iVertex = mesh % verticesOnCell % array(iVertexOnCell,iCell)

          x = mesh % xVertex % array(iVertex)
          y = mesh % yVertex % array(iVertex)

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

    type(mesh_type), intent(in) :: mesh

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

    do iVertex = 1, mesh % nVertices

       do iVertexDegree = 1, mesh % vertexDegree

          iCell = mesh % cellsOnVertex % array(iVertexDegree,iVertex)

          x = mesh % xCell % array(iCell)
          y = mesh % yCell % array(iCell)

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

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(inout) :: &
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND) :: &
         xpy, &
         x, y

    integer :: &
         iCell

    do iCell = 1, mesh % nCells

       x = mesh % xCell % array(iCell)
       y = mesh % yCell % array(iCell)

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

  !-------------------------------------------------------------
  ! spherical tests - strain
  !-------------------------------------------------------------

  subroutine spherical_test_strain(mesh, &
       uVelocity,         &
       vVelocity,         &
       longitudeVelocity, &
       latitudeVelocity,  &
       strain11,          &
       strain22,          &
       strain12,          &
       longitudeStrain,   &
       latitudeStrain,    &
       test_type)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         uVelocity, &
         vVelocity, &
         strain11,  &
         strain22,  &
         strain12

    real(kind=RKIND), dimension(:), intent(in) :: &
         longitudeVelocity, &
         latitudeVelocity,  &
         longitudeStrain,   &
         latitudeStrain

    character(len=*), intent(in) :: &
         test_type

    real(kind=RKIND) :: &
         uVelocity_test, &
         vVelocity_test, &
         du_dlon, &
         du_dlat, &
         dv_dlon, &
         dv_dlat, &
         radius

    integer :: &
         nPointsVelocity, &
         nPointsStrain, &
         iPoint

    nPointsVelocity = size(uVelocity)
    nPointsStrain   = size(strain11)

    radius = mesh % sphere_radius

    ! set velocity points
    do iPoint = 1, nPointsVelocity

       call spherical_test_strain_velocities( &
            uVelocity(iPoint),         &
            vVelocity(iPoint),         &
            du_dlon,                   &
            du_dlat,                   &
            dv_dlon,                   &
            dv_dlat,                   &
            longitudeVelocity(iPoint), &
            latitudeVelocity(iPoint),  &
            test_type)

    enddo ! iPoint

    !write(*,*) "Velocities:", trim(test_type)
    !write(*,*) minval(uVelocity), maxval(uVelocity)
    !write(*,*) minval(vVelocity), maxval(vVelocity)

    ! set strain points
    do iPoint = 1, nPointsStrain

       call spherical_test_strain_velocities( &
            uVelocity_test,          &
            vVelocity_test,          &
            du_dlon,                 &
            du_dlat,                 &
            dv_dlon,                 &
            dv_dlat,                 &
            longitudeStrain(iPoint), &
            latitudeStrain(iPoint),  &
            test_type)

       strain11(iPoint) = strain11_component(uVelocity_test, vVelocity_test, &
                                             du_dlon, du_dlat, dv_dlon, dv_dlat, &
                                             radius, longitudeStrain(iPoint), latitudeStrain(iPoint))
       strain22(iPoint) = strain22_component(uVelocity_test, vVelocity_test, &
                                             du_dlon, du_dlat, dv_dlon, dv_dlat, &
                                             radius, longitudeStrain(iPoint), latitudeStrain(iPoint))
       strain12(iPoint) = strain12_component(uVelocity_test, vVelocity_test, &
                                             du_dlon, du_dlat, dv_dlon, dv_dlat, &
                                             radius, longitudeStrain(iPoint), latitudeStrain(iPoint))

    enddo ! iPoint

    !write(*,*) "Strains:", trim(test_type)
    !write(*,*) minval(strain11), maxval(strain11)
    !write(*,*) minval(strain22), maxval(strain22)
    !write(*,*) minval(strain12), maxval(strain12)

  end subroutine spherical_test_strain

  !-------------------------------------------------------------
  
  subroutine spherical_test_strain_velocities( &
       u,       &
       v,       &
       du_dlon, &
       du_dlat, &
       dv_dlon, &
       dv_dlat, &
       lon,     &
       lat,     &
       test_type)

    real(kind=RKIND), intent(out) :: &
         u,       &
         v,       &
         du_dlon, &
         du_dlat, &
         dv_dlon, &
         dv_dlat
    
    real(kind=RKIND), intent(in) :: &
         lon, &
         lat
    
    character(len=*), intent(in) :: &
         test_type

    if (trim(test_type) == "zero") then

       u = 0.0_RKIND
       v = 0.0_RKIND

       du_dlon = 0.0_RKIND
       du_dlat = 0.0_RKIND
       dv_dlon = 0.0_RKIND
       dv_dlat = 0.0_RKIND

    else if (trim(test_type) == "zonal") then

       u = 1.0_RKIND
       v = 0.0_RKIND

       du_dlon = 0.0_RKIND
       du_dlat = 0.0_RKIND
       dv_dlon = 0.0_RKIND
       dv_dlat = 0.0_RKIND

    else if (trim(test_type) == "meridonal") then

       u = 0.0_RKIND
       v = 1.0_RKIND

       du_dlon = 0.0_RKIND
       du_dlat = 0.0_RKIND
       dv_dlon = 0.0_RKIND
       dv_dlat = 0.0_RKIND

    else if (trim(test_type) == "solid_body") then

       u = cos(lat)
       v = 0.0_RKIND

       du_dlon = 0.0_RKIND
       du_dlat = -sin(lat)
       dv_dlon = 0.0_RKIND
       dv_dlat = 0.0_RKIND

    else if (trim(test_type) == "sinusoidal1") then

       u = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       v = 0.0_RKIND

       du_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       du_dlat = -2.0_RKIND * cos(lon) * sin(2.0_RKIND * lat)
       dv_dlon = 0.0_RKIND
       dv_dlat = 0.0_RKIND

    else if (trim(test_type) == "sinusoidal2") then

       u = 0.0_RKIND
       v = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))

       du_dlon = 0.0_RKIND
       du_dlat = 0.0_RKIND
       dv_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dv_dlat = -2.0_RKIND * cos(lon) * sin(2.0_RKIND * lat)

    endif

  end subroutine spherical_test_strain_velocities

  !-------------------------------------------------------------

  function strain11_component(u, v, du_dlon, du_dlat, dv_dlon, dv_dlat, r, lon, lat) result(e11)

    real(kind=RKIND), intent(in) :: &
         u,       &
         v,       &
         du_dlon, &
         du_dlat, &
         dv_dlon, &
         dv_dlat, &                                  
         r,       &
         lon,     &
         lat

    real(kind=RKIND) :: e11

    e11 = (1.0_RKIND / (r * cos(lat))) * (du_dlon - v * sin(lat)) 

  end function strain11_component

  !-------------------------------------------------------------

  function strain22_component(u, v, du_dlon, du_dlat, dv_dlon, dv_dlat, r, lon, lat) result(e22)

    real(kind=RKIND), intent(in) :: &
         u,       &
         v,       &
         du_dlon, &
         du_dlat, &
         dv_dlon, &
         dv_dlat, &                                  
         r,       &
         lon,     &
         lat

    real(kind=RKIND) :: e22

    e22 = dv_dlat / r

  end function strain22_component

  !-------------------------------------------------------------

  function strain12_component(u, v, du_dlon, du_dlat, dv_dlon, dv_dlat, r, lon, lat) result(e12)

    real(kind=RKIND), intent(in) :: &
         u,       &
         v,       &
         du_dlon, &
         du_dlat, &
         dv_dlon, &
         dv_dlat, &                                  
         r,       &
         lon,     &
         lat

    real(kind=RKIND) :: e12

    e12 = (0.5_RKIND / r) * (du_dlat + u * tan(lat) + dv_dlon / cos(lat))

  end function strain12_component

  !-------------------------------------------------------------
  ! spherical tests - strain
  !-------------------------------------------------------------

  subroutine spherical_test_divergence_stress(mesh, &
       stress11,            &
       stress22,            &
       stress12,            &
       longitudeStress,     &
       latitudeStress,      &
       stressDivergenceU,   &
       stressDivergenceV,   &
       longitudeDivergence, &
       latitudeDivergence,  &
       test_type)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, &
         stressDivergenceV, &
         stress11,          &
         stress22,          &
         stress12
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         longitudeStress,     &
         latitudeStress,      &
         longitudeDivergence, &
         latitudeDivergence

    character(len=*), intent(in) :: &
         test_type

    real(kind=RKIND) :: &
         stress11_test,  &
         stress22_test,  &
         stress12_test,  &
         dstress11_dlon, &
         dstress11_dlat, &
         dstress22_dlon, &
         dstress22_dlat, &
         dstress12_dlon, &
         dstress12_dlat, &
         radius

    integer :: &
         nPointsStress, &
         nPointsDivergence, &
         iPoint

    radius = mesh % sphere_radius

    nPointsStress     = size(stress11)
    nPointsDivergence = size(stressDivergenceU)

    ! set stress points
    do iPoint = 1, nPointsStress

       call spherical_test_divergence_stress_stresses( &
         stress11(iPoint),        &
         stress22(iPoint),        &
         stress12(iPoint),        &
         dstress11_dlon,          &
         dstress11_dlat,          &
         dstress22_dlon,          &
         dstress22_dlat,          &
         dstress12_dlon,          &
         dstress12_dlat,          &
         longitudeStress(iPoint), &
         latitudeStress(iPoint),  &
         test_type)

    enddo ! iPoint

    ! set divergence poinys
    do iPoint = 1, nPointsDivergence

       call spherical_test_divergence_stress_stresses( &
         stress11_test,               &
         stress22_test,               &
         stress12_test,               &
         dstress11_dlon,              &
         dstress11_dlat,              &
         dstress22_dlon,              &
         dstress22_dlat,              &
         dstress12_dlon,              &
         dstress12_dlat,              &
         longitudeDivergence(iPoint), &
         latitudeDivergence(iPoint),  &
         test_type)

       stressDivergenceU(iPoint) = divergenceStressU(stress11_test, stress22_test, stress12_test, &
                                                     dstress11_dlon, dstress11_dlat, dstress22_dlon, dstress22_dlat, dstress12_dlon, dstress12_dlat, &
                                                     radius, longitudeDivergence(iPoint), latitudeDivergence(iPoint))
       stressDivergenceV(iPoint) = divergenceStressV(stress11_test, stress22_test, stress12_test, &
                                                     dstress11_dlon, dstress11_dlat, dstress22_dlon, dstress22_dlat, dstress12_dlon, dstress12_dlat, &
                                                     radius, longitudeDivergence(iPoint), latitudeDivergence(iPoint))

       write(55,*) iPoint, longitudeDivergence(iPoint), latitudeDivergence(iPoint), stressDivergenceU(iPoint), stressDivergenceV(iPoint), radius

    enddo ! iPoint

  end subroutine spherical_test_divergence_stress

  !-------------------------------------------------------------

  subroutine spherical_test_divergence_stress_stresses( &
       stress11,       &
       stress22,       &
       stress12,       &
       dstress11_dlon, &
       dstress11_dlat, &
       dstress22_dlon, &
       dstress22_dlat, &
       dstress12_dlon, &
       dstress12_dlat, &
       lon,            &
       lat,            &
       test_type)
    
    real(kind=RKIND), intent(out) :: &
         stress11,       &
         stress22,       &
         stress12,       &
         dstress11_dlon, &
         dstress11_dlat, &
         dstress22_dlon, &
         dstress22_dlat, &
         dstress12_dlon, &
         dstress12_dlat

    real(kind=RKIND), intent(in) :: &
         lon, &
         lat

    character(len=*), intent(in) :: &
         test_type

    if (trim(test_type) == "zero") then
       
       stress11 = 0.0_RKIND
       stress22 = 0.0_RKIND
       stress12 = 0.0_RKIND
       
       dstress11_dlon = 0.0_RKIND
       dstress11_dlat = 0.0_RKIND
       dstress22_dlon = 0.0_RKIND
       dstress22_dlat = 0.0_RKIND
       dstress12_dlon = 0.0_RKIND
       dstress12_dlat = 0.0_RKIND

    else if (trim(test_type) == "const11") then

       stress11 = 1.0_RKIND
       stress22 = 0.0_RKIND
       stress12 = 0.0_RKIND
       
       dstress11_dlon = 0.0_RKIND
       dstress11_dlat = 0.0_RKIND
       dstress22_dlon = 0.0_RKIND
       dstress22_dlat = 0.0_RKIND
       dstress12_dlon = 0.0_RKIND
       dstress12_dlat = 0.0_RKIND

    else if (trim(test_type) == "const22") then

       stress11 = 0.0_RKIND
       stress22 = 1.0_RKIND
       stress12 = 0.0_RKIND
       
       dstress11_dlon = 0.0_RKIND
       dstress11_dlat = 0.0_RKIND
       dstress22_dlon = 0.0_RKIND
       dstress22_dlat = 0.0_RKIND
       dstress12_dlon = 0.0_RKIND
       dstress12_dlat = 0.0_RKIND

    else if (trim(test_type) == "const12") then

       stress11 = 0.0_RKIND
       stress22 = 0.0_RKIND
       stress12 = 1.0_RKIND
       
       dstress11_dlon = 0.0_RKIND
       dstress11_dlat = 0.0_RKIND
       dstress22_dlon = 0.0_RKIND
       dstress22_dlat = 0.0_RKIND
       dstress12_dlon = 0.0_RKIND
       dstress12_dlat = 0.0_RKIND

    else if (trim(test_type) == "test1") then

       stress11 = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       stress22 = 0.0_RKIND
       stress12 = 0.0_RKIND
       
       dstress11_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dstress11_dlat = -cos(lon) * 2.0_RKIND * sin(2.0_RKIND * lat)
       dstress22_dlon = 0.0_RKIND
       dstress22_dlat = 0.0_RKIND
       dstress12_dlon = 0.0_RKIND
       dstress12_dlat = 0.0_RKIND

    else if (trim(test_type) == "test2") then

       stress11 = 0.0_RKIND
       stress22 = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       stress12 = 0.0_RKIND
       
       dstress11_dlon = 0.0_RKIND
       dstress11_dlat = 0.0_RKIND
       dstress22_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dstress22_dlat = -cos(lon) * 2.0_RKIND * sin(2.0_RKIND * lat)
       dstress12_dlon = 0.0_RKIND
       dstress12_dlat = 0.0_RKIND

    else if (trim(test_type) == "test3") then

       stress11 = 0.0_RKIND
       stress22 = 0.0_RKIND
       stress12 = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       
       dstress11_dlon = 0.0_RKIND
       dstress11_dlat = 0.0_RKIND
       dstress22_dlon = 0.0_RKIND
       dstress22_dlat = 0.0_RKIND
       dstress12_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dstress12_dlat = -cos(lon) * 2.0_RKIND * sin(2.0_RKIND * lat)

    else if (trim(test_type) == "test4") then

       stress11 = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       stress22 = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       stress12 = cos(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       
       dstress11_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dstress11_dlat = -cos(lon) * 2.0_RKIND * sin(2.0_RKIND * lat)
       dstress22_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dstress22_dlat = -cos(lon) * 2.0_RKIND * sin(2.0_RKIND * lat)
       dstress12_dlon = -sin(lon) * (1.0_RKIND + cos(2.0_RKIND * lat))
       dstress12_dlat = -cos(lon) * 2.0_RKIND * sin(2.0_RKIND * lat)

    else
       
       write(*,*) "Unknown test case"
       stop

    endif

  end subroutine spherical_test_divergence_stress_stresses

  !-------------------------------------------------------------

  function divergenceStressU(stress11, stress22, stress12, &
                             dstress11_dlon, dstress11_dlat, dstress22_dlon, dstress22_dlat, dstress12_dlon, dstress12_dlat, &
                             r, lon, lat) result(divu)

    real(kind=RKIND), intent(in) :: &
         stress11,       &
         stress22,       &
         stress12,       &
         dstress11_dlon, &
         dstress11_dlat, &
         dstress22_dlon, &
         dstress22_dlat, &
         dstress12_dlon, &
         dstress12_dlat, &
         r,              &
         lon,            &
         lat

    real(kind=RKIND) :: divu

    divu = (1.0_RKIND / (r * cos(lat))) * dstress11_dlon + &
           (1.0_RKIND / r)              * dstress12_dlat - &
           (2.0_RKIND / r) * tan(lat)   * stress12

  end function divergenceStressU

  !-------------------------------------------------------------

  function divergenceStressV(stress11, stress22, stress12, &
                             dstress11_dlon, dstress11_dlat, dstress22_dlon, dstress22_dlat, dstress12_dlon, dstress12_dlat, &
                             r, lon, lat) result(divv)

    real(kind=RKIND), intent(in) :: &
         stress11,       &
         stress22,       &
         stress12,       &
         dstress11_dlon, &
         dstress11_dlat, &
         dstress22_dlon, &
         dstress22_dlat, &
         dstress12_dlon, &
         dstress12_dlat, &
         r,              &
         lon,            &
         lat

    real(kind=RKIND) :: divv

    divv = (1.0_RKIND / (r * cos(lat))) * dstress12_dlon + &
           (1.0_RKIND / r)              * dstress22_dlat + &
           (1.0_RKIND / r) * tan(lat)   * stress11       - &
           (1.0_RKIND / r) * tan(lat)   * stress22

  end function divergenceStressV

  !-------------------------------------------------------------

end module cice_testing

