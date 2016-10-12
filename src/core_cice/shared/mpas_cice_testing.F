!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_testing
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_testing

  use mpas_derived_types
  use mpas_pool_routines

  implicit none

  private
  save

  public :: &
       cice_divergence_stress_test_velocity_set, &
       cice_divergence_stress_test_stress_set_hex, &
       cice_divergence_stress_test_stress_set_tri, &
       cice_divergence_stress_test_stress_set_weak, &
       cice_init_square_test_case_hex, &
       cice_init_square_point_test_case_hex, &
       cice_init_spherical_test_case

contains

!-----------------------------------------------------------------------
! Spherical test case
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_spherical_test_case
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_init_spherical_test_case(&
       mesh, &
       iceAreaCell, &
       iceVolumeCell, &
       totalMassCell, &
       surfaceTemperature, &
       airTemperature, &
       uVelocity, &
       vVelocity, &
       interiorVertex)

    use cice_constants, only: omega
    use cice_constants, only: ciceDegreesToRadians
    use cice_debug, only: cice_init_ijpop_from_ivertex

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
         initialIceEdgeLatitudeNorthernHemisphere =  70.0_RKIND * ciceDegreesToRadians, &
         initialIceEdgeLatitudeSouthernHemisphere = -60.0_RKIND * ciceDegreesToRadians, &
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
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(mesh, "fVertex", fVertex)
    call MPAS_pool_get_array(mesh, "latVertex", latVertex)
    call MPAS_pool_get_array(mesh, "latCell", latCell)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    call cice_init_ijpop_from_ivertex(mesh)

    fVertex = 2.0_RKIND * omega * sin(latVertex)

    if (trim(test_type) == "icecaps") then

       do iCell = 1, nCells

          if (latCell(iCell) > initialIceEdgeLatitudeNorthernHemisphere .or. &
              latCell(iCell) < initialIceEdgeLatitudeSouthernHemisphere) then
          !if (latCell(iCell) < initialIceEdgeLatitudeSouthernHemisphere) then

             ! ice present
             iceAreaCell(iCell)        = 1.0_RKIND
             iceVolumeCell(iCell)      = iceAreaCell(iCell) * iceHeightCell
             surfaceTemperature(iCell) = 0.0_RKIND

          else

             ! no ice
             iceAreaCell(iCell)   = 0.0_RKIND
             iceVolumeCell(iCell) = 0.0_RKIND
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

  end subroutine cice_init_spherical_test_case

!-----------------------------------------------------------------------
! Square test case
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_square_test_case_hex
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_init_square_test_case_hex(&
       block, &
       configs)

    type(block_type), intent(inout) :: &
         block

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         velocity_solver, &
         boundary, &
         atmos_forcing, &
         ocean_coupling

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         uAirVelocity, &
         vAirVelocity, &
         airDensity, &
         uOceanVelocity, &
         vOceanVelocity

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmos_forcing)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)

    call MPAS_pool_get_array(atmos_forcing, "uAirVelocity", uAirVelocity)
    call MPAS_pool_get_array(atmos_forcing, "vAirVelocity", vAirVelocity)
    call MPAS_pool_get_array(atmos_forcing, "airDensity", airDensity)

    call MPAS_pool_get_array(ocean_coupling, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(ocean_coupling, "vOceanVelocity", vOceanVelocity)

    call square_test_correct_positions(mesh)

    call init_square_test_case_atmos(&
         uAirVelocity, &
         vAirVelocity, &
         airDensity,   &
         xCell,        &
         yCell,        &
         0.0_RKIND)

    call init_square_test_case_ocean(&
         uOceanVelocity, &
         vOceanVelocity, &
         xCell,          &
         yCell)

    call init_square_test_case_state(&
         mesh,          &
         tracers,       &
         velocity_solver,        &
         ocean_coupling, &
         boundary)

  end subroutine cice_init_square_test_case_hex

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_square_test_case_ocean
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_square_test_case_ocean(&
       uOceanVelocity, &
       vOceanVelocity, &
       x, &
       y)

    real(kind=RKIND), dimension(:), intent(out) :: &
         uOceanVelocity, &
         vOceanVelocity

    real(kind=RKIND), dimension(:), intent(in) :: &
         x, &
         y

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

  end subroutine init_square_test_case_ocean

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_square_test_case_atmos
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_square_test_case_atmos(&
       uAirVelocity, &
       vAirVelocity, &
       airDensity, &
       xin, &
       yin, &
       time)

    use cice_constants, only: pii

    real(kind=RKIND), dimension(:), intent(out) :: &
         uAirVelocity, &
         vAirVelocity, &
         airDensity

    real(kind=RKIND), dimension(:), intent(in) :: &
         xin, &
         yin

    real(kind=RKIND), intent(in) :: &
         time

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

       ! velocities
       uAirVelocity(iPoint) = &!a * (y / Ly)
            a + (sin((2.0_RKIND * pii * time) / theta) - b) * sin(2.0_RKIND * pii * (x / Lx)) * sin(pii * (y / Ly))
       vAirVelocity(iPoint) = &!0.0_RKIND!&
            a + (sin((2.0_RKIND * pii * time) / theta) - b) * sin(2.0_RKIND * pii * (y / Ly)) * sin(pii * (x / Lx))

       !uAirVelocity(iPoint) = a
       !vAirVelocity(iPoint) = a

       !write(*,*) iPoint, x, y, uAirVelocity(iPoint), vAirVelocity(iPoint), time

       ! air density
       airDensity(iPoint) = 1.3_RKIND

    enddo ! iPoint

  end subroutine init_square_test_case_atmos

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_square_test_case_state
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_square_test_case_state(&
       mesh, &
       tracers, &
       velocity_solver, &
       ocean_coupling, &
       boundary)

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh

    type(MPAS_pool_type), pointer :: &
         tracers, &
         velocity_solver, &
         ocean_coupling, &
         boundary

    real(kind=RKIND) :: &
         iceThickness

    real(kind=RKIND) :: &
         x

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
         uVelocity, &
         vVelocity, &
         uOceanVelocity, &
         vOceanVelocity, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(mesh, "xCell", xCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)

    call MPAS_pool_get_array(velocity_solver, "uOceanVelocityVertex", uOceanVelocityVertex)
    call MPAS_pool_get_array(velocity_solver, "vOceanVelocityVertex", vOceanVelocityVertex)

    call MPAS_pool_get_array(ocean_coupling, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(ocean_coupling, "vOceanVelocity", vOceanVelocity)

    ! thickness and area

    iceThickness = 2.0_RKIND

    do iCell = 1, nCells

       x = xCell(iCell)

       iceAreaCategory(1,:,iCell) = &!0.95_RKIND!&
            max(min(x / Lx, 1.0_RKIND), 0.0_RKIND)

       iceVolumeCategory(1,:,iCell) = &
            iceThickness * iceAreaCategory(1,:,iCell)

       snowVolumeCategory(1,:,iCell) = 0.0_RKIND

    enddo ! iCell

  end subroutine init_square_test_case_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  square_test_correct_positions
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine square_test_correct_positions(mesh)

    use cice_constants, only: &
         ciceRadiansToDegrees

    type(MPAS_pool_type), intent(inout) :: &
         mesh

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
         distanceToAngle = ciceRadiansToDegrees / radius, &
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
         fVertex

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
    call MPAS_pool_get_array(mesh, "fVertex", fVertex)

    ! Cell
    xCell = xCell + dx
    yCell = yCell + dy

    lonCell = xCell * distanceToAngle + longitudeSquare
    latCell = yCell * distanceToAngle + latitudeSquare

    lonCell = lonCell / ciceRadiansToDegrees
    latCell = latCell / ciceRadiansToDegrees

    ! Vertex
    xVertex = xVertex + dx
    yVertex = yVertex + dy

    lonVertex = xVertex * distanceToAngle + longitudeSquare
    latVertex = yVertex * distanceToAngle + latitudeSquare

    lonVertex = lonVertex / ciceRadiansToDegrees
    latVertex = latVertex / ciceRadiansToDegrees

    ! Edge
    xEdge = xEdge + dx
    yEdge = yEdge + dy

    lonEdge = xEdge * distanceToAngle + longitudeSquare
    latEdge = yEdge * distanceToAngle + latitudeSquare

    lonEdge = lonEdge / ciceRadiansToDegrees
    latEdge = latEdge / ciceRadiansToDegrees

    ! fvalues
    fVertex = 2.0_RKIND * omega * sin(latVertex)

  end subroutine square_test_correct_positions

!-----------------------------------------------------------------------
! point square test case
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_square_point_test_case_hex
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_init_square_point_test_case_hex(&
       block, &
       configs)

    type(block_type), intent(inout) :: &
         block

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs

    type(MPAS_pool_type), pointer :: &
         tracers, &
         mesh

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         surfaceTemperature, &
         iceVolumeCategory, &
         snowVolumeCategory

    integer, parameter :: &
         iCellPoint = 81*40   ! global index for first cell with ice
!!         iCellPoint = 420   ! global index for first cell with ice

    !WHL - added the remaining variables for parallel runs
    integer, dimension(:),   pointer :: indexToCellID
    integer, pointer :: nCellsSolve
    integer :: iCell

    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call mpas_pool_get_array(mesh, 'indexToCellID', indexToCellID)
    call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)

    !WHL - iCellPoint is a global cell index.  For parallel runs, we need to check whether
    !      the cell with this global index is owned by the local processor.

    !WHL - old code commented out
!    iceAreaCategory(:,:,iCellPoint:iCellPoint+9) = 1.0_RKIND
!    surfaceTemperature(:,:,iCellPoint:iCellPoint+9) = 1.0_RKIND
!    iceVolumeCategory(:,:,iCellPoint:iCellPoint+9) = 1.0_RKIND
!    snowVolumeCategory(:,:,iCellPoint:iCellPoint+9) = 1.0_RKIND

!!    write(stderrUnit,*) 'Initial cells with ice (localID, globalID):'

    do iCell = 1, nCellsSolve
       if (indexToCellID(iCell) >= iCellPoint .and. indexToCellID(iCell) <= iCellPoint+9) then
          iceAreaCategory(:,:,iCell) = 1.0_RKIND
          surfaceTemperature(:,:,iCell) = 1.0_RKIND
          iceVolumeCategory(:,:,iCell) = 1.0_RKIND
          snowVolumeCategory(:,:,iCell) = 1.0_RKIND
!!          write(stderrUnit,*) iCell, indexToCellID(iCell)
       endif
    enddo
!!    write(stderrUnit,*) ' '

  end subroutine cice_init_square_point_test_case_hex

!-----------------------------------------------------------------------
! stress divergence operator test velocities
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_divergence_stress_test_velocity_set
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_divergence_stress_test_velocity_set(&
       uVelocity, &
       vVelocity, &
       x, &
       y, &
       type)

    real(kind=RKIND), dimension(:), intent(out) :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:), intent(in) :: &
         x, &
         y

    character(len=*), intent(in) :: &
         type

    integer :: &
         iPoint, &
         nPoints

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
          uVelocity(iPoint) = x(iPoint)**2 / 4.0_RKIND + x(iPoint) * y(iPoint)
          vVelocity(iPoint) = y(iPoint)**2 / 4.0_RKIND + x(iPoint) * y(iPoint)
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

  end subroutine cice_divergence_stress_test_velocity_set

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_divergence_stress_test_stress_set_hex
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_divergence_stress_test_stress_set_hex(&
       mesh, &
       stress1, &
       stress2, &
       stress12)

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh

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
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

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

  end subroutine cice_divergence_stress_test_stress_set_hex

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_divergence_stress_test_stress_set_tri
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_divergence_stress_test_stress_set_tri(&
       mesh, &
       stress1, &
       stress2, &
       stress12)

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh

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
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)

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

  end subroutine cice_divergence_stress_test_stress_set_tri

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_divergence_stress_test_stress_set_weak
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_divergence_stress_test_stress_set_weak(&
       mesh, &
       stress11, &
       stress22, &
       stress12)

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh

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
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

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

  end subroutine cice_divergence_stress_test_stress_set_weak

!-----------------------------------------------------------------------

end module cice_testing
