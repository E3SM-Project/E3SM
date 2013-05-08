module cice_testing

  use mpas_grid_types

  implicit none

  private
  public :: divergence_stress_test_velocity_set, &
            divergence_stress_test_stress_set_hex, &
            divergence_stress_test_stress_set_tri, &
            divergence_stress_test_stress_set_weak, &
            gnuplot_cell, &
            gnuplot_triangle, &
            matlab_jet, &
            plot_boundary_triangles, &
            boundary_locations, &
            labels_verticesOnCell, &
            labels_verticesDegreeOnVertex, &
            writeout_minmax, &
            writeout_on_cell_real, &
            writeout_on_vertex_real, &
            writeout_on_vertex_on_cell_real, &
            writeout_array, &
            writeout_edgecell, &
            writeout_edgeedgecell, &
            rotate_ninety, &
            reverse_grids, &
            gnuplot_vertexvector, &
            find_nearest_cell, &
            init_square_test_case_hex, &
            plot_cells_spherical, &
            plot_pop_grid

  integer, public :: iCellTest

  integer :: &
       iObject = 1, &
       iLabel  = 1, &
       iArrow  = 1

  real(kind=RKIND) :: &
       xmin =  1.0e30_RKIND, &
       xmax = -1.0e30_RKIND, &
       ymin =  1.0e30_RKIND, &
       ymax = -1.0e30_RKIND

    real(kind=RKIND), parameter :: &
         pi = 3.14159265359_RKIND, &
         deg_to_rad = pi / 180.0_RKIND

  ! parameter in constant stress strain relation
  real(kind=RKIND), parameter, public :: lambda = 1.0_RKIND

contains

  !-------------------------------------------------------------
  ! Square test case
  !-------------------------------------------------------------

  subroutine init_square_test_case_hex(mesh, icestate, hexdyn, hexfor, boundary)

    type(mesh_type), intent(inout) :: mesh

    type(icestate_type), pointer :: icestate
    type(hexdyn_type),   pointer :: hexdyn
    type(hexfor_type),   pointer :: hexfor
    type(boundary_type), pointer :: boundary

    call square_test_correct_positions(mesh)

    call init_atmos_velocity(hexfor % uAirVelocity % array, &
                             hexfor % vAirVelocity % array, &
                             mesh % xCell % array,        &
                             mesh % yCell % array,        &
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
         Lx = 1.28e6_RKIND, &
         rhoi = 917.0_RKIND

    integer :: &
         iCell, &
         iVertex

    iceThickness = 2.0_RKIND

    do iCell = 1, mesh % nCells

       x = mesh % xCell % array(iCell)

       icestate % iceAreaCell % array(iCell) = &!0.95_RKIND!&
            max(min(x / Lx, 1.0_RKIND), 0.0_RKIND)

       icestate % iceVolumeCell % array(iCell) = &
            iceThickness * icestate % iceAreaCell % array(iCell)

       icestate % totalMassCell % array(iCell) = icestate % iceVolumeCell % array(iCell) * rhoi
          
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
  ! general plotting
  !-------------------------------------------------------------
  
  subroutine gnuplot_cell(mesh, cellArray, filename, append, nofill) 
    
    type(mesh_type), intent(in) :: mesh
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         cellArray
    
    character(len=*), intent(in) :: &
         filename
    
    logical, optional, intent(in) :: &
         append, &
         nofill

    real(kind=RKIND) :: &
         x, y, x0, y0, xc, yc

    real(kind=RKIND) :: &
         fMin, fMax, &
         fValue

    integer :: &
         iCell, &
         iVertexOnCell, &
         iVertex

    character(len=400) :: &
         stroutvertex, &
         strout, &
         stroutlabel, &
         stroutint

    integer :: &
         red, &
         green, &
         blue

    character(len=7) :: &
         color

    ! optional arguments 
    logical :: &
         lappend, &
         lnofill

    ! optional arguments
    lappend = .false.
    if (present(append)) lappend = append

    lnofill = .false.
    if (present(nofill)) lnofill = nofill

    fMin = minval(cellArray(1:mesh % nCells))
    fMax = maxval(cellArray(1:mesh % nCells))
    !fMin = 1
    !fMax = mesh % nCells

    if (lappend) then
       open(55,file=filename,action="write",position="append")
    else
       open(55,file=filename,action="write")
    end if

    open(56,file="border.txt",action="write")

    write(56,*) minval(mesh % cellsOnCell % array), maxval(mesh % cellsOnCell % array)
    write(56,*) minloc(mesh % cellsOnCell % array), maxloc(mesh % cellsOnCell % array)
    write(56,*) size(mesh % cellsOnCell % array,1), size(mesh % cellsOnCell % array,2)

    do iCell = 1, mesh % nCells

       xc = mesh % xCell % array(iCell)
       yc = mesh % yCell % array(iCell)

       write(*,*) iCell, xc, yc

       xmin = min(xmin,xc) ; xmax = max(xmax,xc)
       ymin = min(ymin,yc) ; ymax = max(ymax,yc)

       !write(stroutint,fmt='(i5)') iCell
       !write(stroutlabel,fmt='(a,i5,a,a,a,f14.2,a,f14.2,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at ',xc, ",",yc, " center"
       !iLabel = iLabel + 1
       !write(55,*) trim(stroutlabel)

       fValue = (cellArray(iCell) - fMin) / (fMax - fMin) 
       !fValue = (real(iCell,RKIND) - fMin) / (fMax - fMin) 

       write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "

       do iVertexOnCell = 1, mesh % nEdgesOnCell % array(iCell)

          !write(*,*) iCell, iVertexOnCell

          write(56,*) iCell, mesh % nCells, iVertexOnCell, mesh % cellsOnCell % array(iVertexOnCell, iCell)

          iVertex = mesh % verticesOnCell % array(iVertexOnCell, iCell)

          x = mesh % xVertex % array(iVertex)
          y = mesh % yVertex % array(iVertex)

          xmin = min(xmin,x) ; xmax = max(xmax,x)
          ymin = min(ymin,y) ; ymax = max(ymax,y)

          if (iVertexOnCell == 1) then
             x0 = x
             y0 = y
          endif

          write(stroutvertex,fmt='(f14.2,a,f14.2,a)') x, ",", y, " to "

          strout = trim(strout)//trim(stroutvertex)

       enddo ! iVertex

       write(stroutvertex,fmt='(f14.2,a,f14.2)') x0, ",", y0

       strout = trim(strout)//trim(stroutvertex)

       write(55,*) trim(strout)

       color = matlab_jet(fValue)
       !if (iCell == 3272) color = "#000000"

       if (lnofill) then
          write(strout,fmt='(a,i5,a,a,a)') "set object ",iObject,' fillstyle empty border lt -1'
       else
          write(strout,fmt='(a,i5,a,a,a)') "set object ",iObject,' fc rgb "', color, '" fillstyle solid'
       endif

       write(55,*) trim(strout)

       iObject = iObject + 1

    enddo ! iCell

    write(stroutint,fmt='(e20.8)') fMax
    write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.9'
    iLabel = iLabel + 1
    write(55,*) trim(stroutlabel)


    write(stroutint,fmt='(e20.8)') fMin
    write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.8'
    iLabel = iLabel + 1
    write(55,*) trim(stroutlabel)

    close(55)
    close(56)

    !call gnuplot_triangle(mesh, mesh % xVertex % array, filename, .true., mesh % nCells, .true.)
!stop
  end subroutine gnuplot_cell

  !-------------------------------------------------------------

  subroutine gnuplot_triangle(mesh, vertexArray, filename, append, nofill)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(in) :: &
         vertexArray

    character(len=*), intent(in) :: &
         filename

    logical, optional, intent(in) :: &
         append, &
         nofill

    real(kind=RKIND) :: &
         x, y, x0, y0, xv, yv

    real(kind=RKIND) :: &
         fMin, fMax, &
         fValue

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell

    character(len=400) :: &
         stroutvertex, &
         strout, &
         stroutlabel, &
         stroutint

    integer :: &
         red, &
         green, &
         blue

    character(len=7) :: &
         color

    logical :: atBoundary

    ! optional arguments 
    logical :: &
         lappend, &
         lnofill

    ! optional arguments
    lappend = .false.
    if (present(append)) lappend = append

    lnofill = .false.
    if (present(nofill)) lnofill = nofill

    ! min/maxes
    fMin = 1e30
    fMax = -1e30

    ! loop over triangles
    do iVertex = 1, mesh % nVertices
       atBoundary = .false.

       ! loop over triangle vertices
       do iVertexDegree = 1, mesh % vertexDegree

          iCell = mesh % cellsOnVertex % array(iVertexDegree,iVertex)

          if (iCell > mesh % nCells) then
             atBoundary = .true.
          end if
       enddo ! iVertexDegree

       if (.not. atBoundary) then
          fMin = min(fMin,vertexArray(iVertex))
          fMax = max(fMax,vertexArray(iVertex))
       endif
    
    enddo ! iVertex
    write(*,*) "gt:", fMin, fMax! ; stop

    !fMin = minval(vertexArray)
    !fMax = maxval(vertexArray)
    !fMin = 1
    !fMax = mesh % nVertices

    if (lappend) then
       open(55,file=filename,action="write",position="append")
    else
       open(55,file=filename,action="write")
    end if

    ! loop over triangles
    do iVertex = 1, mesh % nVertices

       xv = mesh % xVertex % array(iVertex)
       yv = mesh % yVertex % array(iVertex)

       !write(*,*) iVertex, xv, yv

       xmin = min(xmin,xv) ; xmax = max(xmax,xv)
       ymin = min(ymin,yv) ; ymax = max(ymax,yv)

       !write(stroutint,fmt='(i5)') iVertex
       !write(stroutlabel,fmt='(a,i5,a,a,a,f14.2,a,f14.2,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at ',xv, ",",yv, " center"
       !iLabel = iLabel + 1
       !write(55,*) trim(stroutlabel)

       atBoundary = .false.

       fValue = (vertexArray(iVertex) - fMin) / (fMax - fMin) 
       !fValue = (real(iVertex,RKIND) - fMin) / (fMax - fMin) 
       !write(*,*) iVertex, vertexArray(iVertex), fMin, fMax

       write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "

       ! loop over triangle vertices
       do iVertexDegree = 1, mesh % vertexDegree

          iCell = mesh % cellsOnVertex % array(iVertexDegree,iVertex)

          if (iCell > mesh % nCells) then
             atBoundary = .true.
          end if

          x = mesh % xCell % array(iCell)
          y = mesh % yCell % array(iCell)

          if (.not. atBoundary) then
             xmin = min(xmin,x) ; xmax = max(xmax,x)
             ymin = min(ymin,y) ; ymax = max(ymax,y)
          endif

          if (iVertexDegree == 1) then
             x0 = x
             y0 = y
          endif

          write(stroutvertex,fmt='(f14.2,a,f14.2,a)') x, ",", y, " to "

          strout = trim(strout)//trim(stroutvertex)

       enddo ! iVertexDegree

       write(stroutvertex,fmt='(f14.2,a,f14.2)') x0, ",", y0

       strout = trim(strout)//trim(stroutvertex)

       if (.not. atBoundary) write(55,*) trim(strout)

       color = matlab_jet(fValue)

       if (lnofill) then
          write(strout,fmt='(a,i5,a,a,a)') "set object ",iObject,' fillstyle empty border lt -1'
       else
          write(strout,fmt='(a,i5,a,a,a)') "set object ",iObject,' fc rgb "', color, '" fillstyle solid'
       endif

       if (.not. atBoundary) iObject = iObject + 1

       if (.not. atBoundary) write(55,*) trim(strout)

    enddo ! iVertex

    write(stroutint,fmt='(e20.8)') fMax
    write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.9 rotate right'
    iLabel = iLabel + 1
    write(55,*) trim(stroutlabel)

    write(stroutint,fmt='(e20.8)') fMin
    write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.6 rotate right'
    iLabel = iLabel + 1
    write(55,*) trim(stroutlabel)

    close(55)

  end subroutine gnuplot_triangle

  !-------------------------------------------------------------

  subroutine gnuplot_vertexvector(mesh, vertexArrayU, vertexArrayV, filename)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(in) :: &
         vertexArrayU, &
         vertexArrayV

    character(len=*), intent(in) :: &
         filename

    integer, parameter :: plot_every = 10

    real(kind=RKIND), parameter :: &
         plotLength = 100000.0_RKIND

    character(len=200) :: strout

    real(kind=RKIND) :: &
         x1, x2, y1, y2, &
         n1, n2, &
         arrowLength, &
         arrowLengthMax, &
         arrowLengthPlot, &
         random

    integer :: iVertex

    ! find maximum vector
    arrowLengthMax = -1e30

    do iVertex = 1, mesh % nVertices

       n1 = vertexArrayU(iVertex)
       n2 = vertexArrayV(iVertex)
       arrowLength = sqrt(n1**2 + n2**2)

       arrowLengthMax = Max(arrowLengthMax,arrowLength)

    enddo ! iVertex

    write(*,*) arrowLengthMax


    open(55,file=filename)

    do iVertex = 1, mesh % nVertices
    
       random = rand()

       if (random <= 1.0_RKIND / real(plot_every,RKIND)) then


          n1 = vertexArrayU(iVertex)
          n2 = vertexArrayV(iVertex)
          arrowLength = sqrt(n1**2 + n2**2)
          n1 = n1 / arrowLength
          n2 = n2 / arrowLength

          arrowLengthPlot = (arrowLength / arrowLengthMax) * plotLength

          x1 = mesh % xVertex % array(iVertex)
          y1 = mesh % yVertex % array(iVertex)

          x2 = x1 + n1 * arrowLengthPlot
          y2 = y1 + n2 * arrowLengthPlot

          write(*,*) arrowLengthPlot, n1, n2, x1, y1, x2, y2

          write(strout,fmt='(a,i5,a,f14.2,a,f14.2,a,f14.2,a,f14.2)') "set arrow ",iArrow," from ",x1,",",y1," to ",x2,",",y2
          iArrow = iArrow + 1

          write(55,*) trim(strout)

       endif

    enddo ! iVertex

    close(55)

  end subroutine gnuplot_vertexvector

  !-------------------------------------------------------------

  function hexstring(red,green,blue) result(hex)

    integer, intent(in) :: red,green,blue
    character(len=7) :: hex

    integer :: hex1, hex2

    hex = "#"//dectohex(red)//dectohex(green)//dectohex(blue)

  end function hexstring

  !-------------------------------------------------------------

  function dectohex(dec) result(hex)

    integer, intent(in) :: dec
    character(len=2) :: hex

    integer :: hex1, hex2

    character(len=1), dimension(16), parameter :: hexchars = &
         (/'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'/)

    hex1 = floor(real(dec, RKIND) / 16.0_RKIND)
    hex2 = dec - hex1 * 16

    hex = hexchars(hex1+1)//hexchars(hex2+1)

  end function dectohex

  !-------------------------------------------------------------

  function matlab_jet(f) result(color)

    real(kind=RKIND), intent(in) :: f
    character(len=7) :: color

    real(kind=RKIND) :: r_red, r_green, r_blue
    integer :: red, green, blue

    r_blue  = max(min(min(4.0_RKIND * f + 0.5_RKIND, -4.0_RKIND * f + 2.5_RKIND),1.0_RKIND),0.0_RKIND)

    r_green = max(min(min(4.0_RKIND * f - 0.5_RKIND, -4.0_RKIND * f + 3.5_RKIND),1.0_RKIND),0.0_RKIND)

    r_red   = max(min(min(4.0_RKIND * f - 1.5_RKIND, -4.0_RKIND * f + 4.5_RKIND),1.0_RKIND),0.0_RKIND)

    red   = nint(r_red   * 255)
    green = nint(r_green * 255)
    blue  = nint(r_blue  * 255)

    color = hexstring(red,green,blue)

  end function matlab_jet

  !-------------------------------------------------------------

  subroutine boundary_locations(mesh, boundaryCell, boundaryVertex, interiorVertex, boundaryEdge)

    type(mesh_type), intent(in) :: mesh

    integer, dimension(:), intent(in) :: &
         boundaryCell, &
         boundaryVertex, &
         interiorVertex, &
         boundaryEdge

    integer :: &
         iCell, &
         iVertex, &
         iEdge

    ! identify border cells
    open(55,file="border_cell.txt",action="write")

    do iCell = 1, mesh % nCells

       if (boundaryCell(iCell) == 1) then

          write(55,*) mesh % xCell % array(iCell), mesh % yCell % array(iCell)

       endif

    enddo ! iCell

    close(55)

    ! identify border vertices
    open(55,file="border_vertex.txt",action="write")

    do iVertex = 1, mesh % nVertices

       if (boundaryVertex(iVertex) == 1) then

          write(55,*) mesh % xVertex % array(iVertex), mesh % yVertex % array(iVertex)

       endif

    enddo ! iVertex

    close(55)

    ! identify interior vertices
    open(55,file="interior_vertex.txt",action="write")

    do iVertex = 1, mesh % nVertices

       if (interiorVertex(iVertex) == 1) then

          write(55,*) mesh % xVertex % array(iVertex), mesh % yVertex % array(iVertex)

       endif

    enddo ! iVertex

    close(55)

    ! identify border edges
    open(55,file="border_edge.txt",action="write")

    do iEdge = 1, mesh % nEdges

       if (boundaryEdge(iEdge) == 1) then

          write(55,*) mesh % xEdge % array(iEdge), mesh % yEdge % array(iEdge)

       endif

    enddo ! iVertex

    close(55)

  end subroutine boundary_locations

  !-------------------------------------------------------------

  subroutine plot_boundary_triangles(mesh, boundaryVertex, interiorVertex, boundaryEdge)

    type(mesh_type), intent(in) :: mesh

    integer, dimension(:), intent(in) :: &
         boundaryVertex, &
         interiorVertex, &
         boundaryEdge

    real(kind=RKIND) :: &
         x, y, &
         x0, y0

    character(len=200) :: &
         stroutvertex, &
         strout

    integer :: &
         iEdgeBoundary, &
         iEdge, &
         iCellOnEdge, &
         iCell, &
         iVertexOnEdge, &
         iVertex, &
         iVertexDegree

    integer :: start_index

    start_index = mesh % nCells + mesh % nVertices + 10

    open(55,file="boundary_triangles_edge.txt",action="write")

    ! edge triangles
    iEdgeBoundary = 0

    do iEdge = 1, mesh % nEdges

       write(0,*) iEdge, boundaryEdge(iEdge)

       if (boundaryEdge(iEdge) == 1) then
          ! am on boundary
          iEdgeBoundary = iEdgeBoundary + 1

          write(strout,fmt='(a,i5,a)') "set object ",iEdgeBoundary+start_index," polygon from "

          ! find the interior cell
          do iCellOnEdge = 1, 2

             iCell = mesh % cellsOnEdge %array(iCellOnEdge,iEdge)

             if (iCell >= 1 .and. iCell <= mesh % nCells) then

                ! this cell is interior to domain
                x = mesh % xCell % array(iCell)
                y = mesh % yCell % array(iCell)               
                x0 = x
                y0 = y

                write(stroutvertex,fmt='(f10.2,a,f10.2,a)') x, ",", y, " to "

                strout = trim(strout)//trim(stroutvertex)

             endif

          enddo ! iCellOnEdge

          ! find the adjacent vertices
          do iVertexOnEdge = 1, 2

             iVertex = mesh % verticesOnEdge % array(iVertexOnEdge,iEdge)

             x = mesh % xVertex % array(iVertex)
             y = mesh % yVertex % array(iVertex)

             write(stroutvertex,fmt='(f10.2,a,f10.2,a)') x, ",", y, " to "

             strout = trim(strout)//trim(stroutvertex)

          enddo ! iVertexOnEdge

          write(stroutvertex,fmt='(f10.2,a,f10.2)') x0, ",", y0

          strout = trim(strout)//trim(stroutvertex)

          write(55,*) trim(strout)

          write(strout,fmt='(a,i5,a,a,a)') "set object ",iEdgeBoundary+start_index,' fillstyle empty border lt -1'

          write(55,*) trim(strout)

       endif

    enddo ! iEdge

    close(55)

    open(55,file="boundary_triangles_vertex.txt",action="write")

    ! loop over vertices
    do iVertex = 1, mesh % nVertices

       ! boundary vertex adjacent to two interior cells
       if (boundaryVertex(iVertex) == 1) then

          iEdgeBoundary = iEdgeBoundary + 1

          write(strout,fmt='(a,i5,a)') "set object ",iEdgeBoundary+start_index," polygon from "

          ! starting vertex point
          x = mesh % xVertex % array(iVertex)
          y = mesh % yVertex % array(iVertex)
          x0 = x
          y0 = y

          write(stroutvertex,fmt='(f10.2,a,f10.2,a)') x, ",", y, " to "

          strout = trim(strout)//trim(stroutvertex)

          ! the two interior cells
          do iVertexDegree = 1, mesh % vertexDegree

             ! adjacent cell
             iCell = mesh % cellsOnVertex % array(iVertexDegree, iVertex)

             if (iCell >= 1 .and. iCell <= mesh % nCells) then
                ! interior adjacent cell

                x = mesh % xCell % array(iCell)
                y = mesh % yCell % array(iCell)

                write(stroutvertex,fmt='(f10.2,a,f10.2,a)') x, ",", y, " to "

                strout = trim(strout)//trim(stroutvertex)

             endif

          enddo ! iVertexDegree

          write(stroutvertex,fmt='(f10.2,a,f10.2)') x0, ",", y0

          strout = trim(strout)//trim(stroutvertex)

          write(55,*) trim(strout)

          write(strout,fmt='(a,i5,a,a,a)') "set object ",iEdgeBoundary+start_index,' fillstyle empty border lt -1'

          write(55,*) trim(strout)

       endif

    enddo ! iVertex

    close(55)

  end subroutine plot_boundary_triangles

  !-------------------------------------------------------------

  subroutine labels_verticesOnCell(mesh, verticesOnCell, filename) 

    type(mesh_type), intent(in) :: mesh

    integer, dimension(:,:), intent(in) :: &
         verticesOnCell

    character(len=*), intent(in) :: &
         filename

    integer :: &
         iCell

    do iCell = 1, mesh % nCells

       call label_verticesOnCell(mesh, iCell, verticesOnCell(:,iCell), filename, .true.) 

    enddo ! iCell

  end subroutine labels_verticesOnCell

  !-------------------------------------------------------------

  subroutine label_verticesOnCell(mesh, iCell, verticesOnCell, filename, append) 

    type(mesh_type), intent(in) :: mesh

    integer, intent(in) :: &
         iCell

    integer, dimension(:), intent(in) :: &
         verticesOnCell

    character(len=*), intent(in) :: &
         filename

    logical, intent(in), optional :: &
         append

    real(kind=RKIND) :: &
         x, y, x0, y0, xc, yc, xl, yl

    real(kind=RKIND) :: &
         fMin, fMax

    character(len=200) :: &
         stroutvertex, &
         strout, &
         stroutlabel, &
         stroutint

    integer :: &
         iVertexOnCell, &
         iVertex

    logical :: &
         lappend

    lappend = .false.
    if (present(append)) lappend = .true.

    if (lappend) then
       open(55,file=filename,position="append")
    else
       open(55,file=filename)
    endif

    write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "

    xc = mesh % xCell % array(iCell)
    yc = mesh % yCell % array(iCell)

    ! draw cell
    do iVertexOnCell = 1, mesh % nEdgesOnCell % array(iCell)

       iVertex = mesh % verticesOnCell % array(iVertexOnCell,iCell)

       x = mesh % xVertex % array(iVertex)
       y = mesh % yVertex % array(iVertex)

       if (iVertexOnCell == 1) then
          x0 = x
          y0 = y
       endif

       write(stroutvertex,fmt='(f10.2,a,f10.2,a)') x, ",", y, " to "

       strout = trim(strout)//trim(stroutvertex)

       xl = 0.3_RKIND * xc + 0.7_RKIND * x
       yl = 0.3_RKIND * yc + 0.7_RKIND * y

       !write(stroutint,fmt='(i5)') iVertexOnCell
       write(stroutint,fmt='(i5)') verticesOnCell(iVertexOnCell)
       write(stroutlabel,fmt='(a,i5,a,a,a,f10.2,a,f10.2,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at ',xl, ",",yl, ' center font "Helvetica,8"'

       write(55,*) trim(stroutlabel)

       iLabel = iLabel + 1

    enddo ! iVertexOnCell

    write(stroutvertex,fmt='(f10.2,a,f10.2)') x0, ",", y0

    strout = trim(strout)//trim(stroutvertex)

    write(55,*) trim(strout)

    write(strout,fmt='(a,i5,a,a,a)') "set object ",iObject,' fillstyle empty border lt -1'

    write(55,*) trim(strout)

    iObject = iObject + 1

    close(55)

  end subroutine label_verticesOnCell

  !-------------------------------------------------------------

  subroutine labels_verticesDegreeOnVertex(mesh, verticesDegreeOnVertex, interiorVertex, filename) 

    type(mesh_type), intent(in) :: mesh

    integer, dimension(:,:), intent(in) :: &
         verticesDegreeOnVertex

    integer, dimension(:), intent(in) :: &
         interiorVertex

    character(len=*), intent(in) :: &
         filename

    integer :: &
         iVertex

    do iVertex = 1, mesh % nVertices

       if (interiorVertex(iVertex) == 1) then

          call label_verticesDegreeOnVertex(mesh, iVertex, verticesDegreeOnVertex(:,iVertex), filename, .true.) 

       endif

    enddo ! iVertex

  end subroutine labels_verticesDegreeOnVertex

  !-------------------------------------------------------------

  subroutine label_verticesDegreeOnVertex(mesh, iVertex, verticesDegreeOnVertex, filename, append) 

    type(mesh_type), intent(in) :: mesh

    integer, intent(in) :: &
         iVertex

    integer, dimension(:), intent(in) :: &
         verticesDegreeOnVertex

    character(len=*), intent(in) :: &
         filename

    logical, intent(in), optional :: &
         append

    real(kind=RKIND) :: &
         x, y, x0, y0, xv, yv, xl, yl

    real(kind=RKIND) :: &
         fMin, fMax

    character(len=200) :: &
         stroutvertex, &
         strout, &
         stroutlabel, &
         stroutint

    integer :: &
         iVertexDegree, &
         iCell

    logical :: &
         lappend

    lappend = .false.
    if (present(append)) lappend = .true.

    if (lappend) then
       open(55,file=filename,position="append")
    else
       open(55,file=filename)
    endif

    write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "

    xv = mesh % xVertex % array(iVertex)
    yv = mesh % yVertex % array(iVertex)

    ! draw cell
    do iVertexDegree = 1, mesh % vertexDegree

       iCell = mesh % cellsOnVertex % array(iVertexDegree,iVertex)

       x = mesh % xCell % array(iCell)
       y = mesh % yCell % array(iCell)

       if (iVertexDegree == 1) then
          x0 = x
          y0 = y
       endif

       write(stroutvertex,fmt='(f10.2,a,f10.2,a)') x, ",", y, " to "

       strout = trim(strout)//trim(stroutvertex)

       xl = 0.3_RKIND * xv + 0.7_RKIND * x
       yl = 0.3_RKIND * yv + 0.7_RKIND * y

       !write(stroutint,fmt='(i5)') iVertexDegree
       write(stroutint,fmt='(i5)') verticesDegreeOnVertex(iVertexDegree)
       write(stroutlabel,fmt='(a,i5,a,a,a,f10.2,a,f10.2,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at ',xl, ",",yl, ' center font "Helvetica,8"'
       iLabel = iLabel + 1

       write(55,*) trim(stroutlabel)

    enddo ! iVertexOnCell

    write(stroutvertex,fmt='(f10.2,a,f10.2)') x0, ",", y0

    strout = trim(strout)//trim(stroutvertex)

    write(55,*) trim(strout)

    write(strout,fmt='(a,i5,a,a,a)') "set object ",iObject,' fillstyle empty border lt -1'

    write(55,*) trim(strout)

    iObject = iObject + 1

    close(55)

  end subroutine label_verticesDegreeOnVertex

  !-------------------------------------------------------------

  subroutine writeout_minmax()

    character(len=200) :: strout
    real(kind=RKIND), parameter :: border = 7000.0_RKIND

    real(kind=RKIND) :: xming, xmaxg, yming, ymaxg

    xming = min(xmin,ymin)-border ; yming = xming
    xmaxg = max(xmax,ymax)+border ; ymaxg = xmaxg

    write(strout,fmt='(a,f10.2,a,f10.2,a,f10.2,a,f10.2,a)') &
         "set size square ; set xrange [", xming, ":", xmaxg, "] ; set yrange [", yming, ":", ymaxg, "]"
    write(*,*) trim(strout)

  end subroutine writeout_minmax

  !-------------------------------------------------------------
  ! write out variables routines
  !-------------------------------------------------------------

  subroutine writeout_on_cell_real(mesh, arrayOnCell, filename)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(in) :: &
         arrayOnCell
    
    character(len=*), intent(in) :: &
         filename

    integer :: iCell

    open(55,file=filename)

    do iCell = 1, mesh % nCells

       write(55,*) iCell, arrayOnCell(iCell), &
            mesh % xCell % array(iCell), mesh % yCell % array(iCell)

    enddo ! iCell

    close(55)

  end subroutine writeout_on_cell_real

  !-------------------------------------------------------------

  subroutine writeout_on_vertex_real(mesh, arrayOnVertex, filename)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(in) :: &
         arrayOnVertex
    
    character(len=*), intent(in) :: &
         filename

    integer :: iVertex

    open(55,file=filename)

    do iVertex = 1, mesh % nVertices

       write(55,*) iVertex, arrayOnVertex(iVertex), &
            mesh % xVertex % array(iVertex), mesh % yVertex % array(iVertex)

    enddo ! iVertex

    close(55)

  end subroutine writeout_on_vertex_real

  !-------------------------------------------------------------

  subroutine writeout_on_vertex_on_cell_real(mesh, arrayVertexOnCell, filename)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         arrayVertexOnCell
    
    character(len=*), intent(in) :: &
         filename

    integer :: &
         iCell, &
         iVertexOnCell, &
         iVertex

    open(55,file=filename)

    do iCell = 1, mesh % nCells

       do iVertexOnCell = 1, mesh % nEdgesOnCell % array(iCell)

          iVertex = mesh % verticesOnCell % array(iVertexOnCell,iCell)

          write(55,*) iCell, iVertexOnCell, arrayVertexOnCell(iVertexOnCell,iCell), &
               mesh % xVertex % array(iVertex), mesh % yVertex % array(iVertex), &
               arrayVertexOnCell(iVertexOnCell,iCell) / (mesh % xVertex % array(iVertex) + mesh % yVertex % array(iVertex))

       enddo ! iVertexOnCell

    enddo ! iCell

    close(55)

  end subroutine writeout_on_vertex_on_cell_real

  !-------------------------------------------------------------

  subroutine writeout_array(n, arr, filename)

    integer, intent(in) :: n

    real(kind=RKIND), dimension(n), intent(in) :: &
         arr

    character(len=*), intent(in) :: &
         filename

    integer :: i

    open(55,file=filename)

    do i = 1, n

       write(55,*) i, arr(i)

    enddo ! i

    close(55)

  end subroutine writeout_array

  !-------------------------------------------------------------
    
  subroutine writeout_edgecell(mesh, arrayEdgeCell, filename)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         arrayEdgeCell
    
    character(len=*), intent(in) :: &
         filename

    integer :: &
         iCell, &
         iEdge

    open(55,file=filename)

    do iCell = 1, mesh % nCells
       
       do iEdge = 1, mesh % nEdgesOnCell % array(iCell)

          write(55,*) iCell, iEdge, arrayEdgeCell(iEdge,iCell)

       enddo ! iEdge

    enddo ! iCell

    close(55)

  end subroutine writeout_edgecell

  !-------------------------------------------------------------

  subroutine writeout_edgeedgecell(mesh, arrayEdgeEdgeCell, filename)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         arrayEdgeEdgeCell
    
    character(len=*), intent(in) :: &
         filename

    integer :: &
         iCell, &
         iEdge, &
         jEdge

    open(55,file=filename)

    do iCell = 1, mesh % nCells
       
       do iEdge = 1, mesh % nEdgesOnCell % array(iCell)

          do jEdge = 1, mesh % nEdgesOnCell % array(iCell) 

             write(55,*) iCell, iEdge, jEdge, arrayEdgeEdgeCell(jEdge,iEdge,iCell)

          enddo ! jEdge

       enddo ! iEdge

    enddo ! iCell

    close(55)

  end subroutine writeout_edgeedgecell

  !-------------------------------------------------------------

  function find_nearest_cell(mesh,x,y) result(iNear)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), intent(in) :: &
         x, y

    integer :: iNear

    real(kind=RKIND) :: &
         distance, &
         min_distance

    integer :: &
         iCell

    iNear = -1
    min_distance = 1e30_RKIND

    do iCell = 1, mesh % nCells

       distance = sqrt((x - mesh % xCell % array(iCell))**2 + (y - mesh % yCell % array(iCell))**2)

       if (distance < min_distance) then

          min_distance = distance

          iNear = iCell

       endif

    enddo ! iCell    

  end function find_nearest_cell

  !-------------------------------------------------------------

  subroutine rotate_ninety(mesh)

    type(mesh_type), intent(in) :: mesh

    integer :: &
         iCell, &
         iVertex, &
         iEdge

    real(kind=RKIND) :: x, y

    do iCell = 1, mesh % nCells

       x = mesh % xCell % array(iCell)
       y = mesh % yCell % array(iCell)

       mesh % xCell % array(iCell) = y
       mesh % yCell % array(iCell) = -x

    enddo ! iCell

    do iVertex = 1, mesh % nVertices

       x = mesh % xVertex % array(iVertex)
       y = mesh % yVertex % array(iVertex)

       mesh % xVertex % array(iVertex) = y
       mesh % yVertex % array(iVertex) = -x

    enddo ! iVertex

    do iEdge = 1, mesh % nEdges

       x = mesh % xEdge % array(iEdge)
       y = mesh % yEdge % array(iEdge)

       mesh % xEdge % array(iEdge) = y
       mesh % yEdge % array(iEdge) = -x

    enddo ! iVertex

  end subroutine rotate_ninety

  !-------------------------------------------------------------
  
  subroutine reverse_grids(mesh)

    type(mesh_type), intent(inout) :: mesh

    real(kind=RKIND), dimension(:), allocatable :: tmp_ncells, tmp_nvertices
    integer, dimension(:,:), allocatable :: tmp_maxedgesncells, tmp_vertexdegreenvertices


    integer :: nCells0, nVertices0, vertexDegree0, maxEdges0

    !write(*,*) lbound(mesh % xCell % array), ubound(mesh % xCell % array), size(mesh % xCell % array,1), mesh % nCells
    !stop


    nCells0 = mesh % nCells
    nVertices0 = mesh % nVertices
    vertexDegree0 = mesh % vertexDegree
    maxEdges0 = mesh % maxEdges

    allocate(tmp_ncells(1:mesh % nCells+1))
    allocate(tmp_nvertices(1:mesh % nVertices+1))
    allocate(tmp_maxedgesncells(1:mesh % maxEdges,1:mesh % nCells+1))
    allocate(tmp_vertexdegreenvertices(1:mesh % vertexDegree,1:mesh % nVertices+1))

    ! x position
    tmp_ncells    = mesh % xCell % array
    tmp_nvertices = mesh % xVertex % array
    deallocate(mesh % xCell % array)
    deallocate(mesh % xVertex % array)
    allocate(mesh % xCell % array(1:nVertices0+1))
    allocate(mesh % xVertex % array(1:nCells0+1))
    mesh % xCell % array = tmp_nvertices
    mesh % xVertex % array = tmp_ncells

    ! x position
    tmp_ncells    = mesh % yCell % array
    tmp_nvertices = mesh % yVertex % array
    deallocate(mesh % yCell % array)
    deallocate(mesh % yVertex % array)
    allocate(mesh % yCell % array(1:nVertices0+1))
    allocate(mesh % yVertex % array(1:nCells0+1))
    mesh % yCell % array = tmp_nvertices
    mesh % yVertex % array = tmp_ncells

    ! z position
    tmp_ncells    = mesh % zCell % array
    tmp_nvertices = mesh % zVertex % array
    deallocate(mesh % zCell % array)
    deallocate(mesh % zVertex % array)
    allocate(mesh % zCell % array(1:nVertices0+1))
    allocate(mesh % zVertex % array(1:nCells0+1))
    mesh % zCell % array = tmp_nvertices
    mesh % zVertex % array = tmp_ncells

    ! nEdgesOnCell
    deallocate(mesh % nEdgesOnCell % array)
    allocate(mesh % nEdgesOnCell % array(nVertices0+1))
    mesh % nEdgesOnCell % array = 3


    ! verticesOnCell / cellsOnVertex
    tmp_maxedgesncells = mesh % verticesOnCell % array
    tmp_vertexdegreenvertices = mesh % cellsOnVertex % array
    deallocate(mesh % verticesOnCell % array)
    deallocate(mesh % cellsOnVertex % array)
    allocate(mesh % verticesOnCell % array(vertexDegree0,nVertices0+1))
    allocate(mesh % cellsOnVertex % array(maxEdges0,nCells0+1))
    mesh % verticesOnCell % array = tmp_vertexdegreenvertices
    mesh % cellsOnVertex % array = tmp_maxedgesncells

    ! area
    deallocate(mesh % areaCell % array)
    allocate(mesh % areaCell % array(nVertices0+1))
    mesh % areaCell % array = mesh % areaTriangle % array
    
    ! sizes
    mesh % nCells = nVertices0
    mesh % nVertices = nCells0
    mesh % vertexDegree = maxEdges0
    mesh % maxEdges = vertexDegree0

    deallocate(tmp_ncells)
    deallocate(tmp_nvertices)
    deallocate(tmp_maxedgesncells)
    deallocate(tmp_vertexdegreenvertices)

    !stop
  end subroutine reverse_grids

  !-------------------------------------------------------------
  ! Spherical plotting
  !-------------------------------------------------------------

  subroutine plot_cells_spherical(mesh, cellArray)

    type(mesh_type), intent(in) :: mesh
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         cellArray

    real(kind=RKIND) :: &
         xc, yc, zc, &
         xv, yv, zv, &
         xv0, yv0, zv0

    real(kind=RKIND) :: &
         x0, y0, z0, &
         x1, y1, z1, &
         x2, y2, z2, &
         x3, y3, z3, &
         x4, y4, z4

    real(kind=RKIND) :: &
         lat0, lon0, &
         lat1, lon1, &
         lat2, lon2, &
         lat3, lon3, &
         lat4, lon4

    logical :: lnofill

    real(kind=RKIND) :: &
         fValue, &
         fMin, &
         fMax

    integer :: &
         iCell, &
         iVertexOnCell, &
         iVertex

    integer, parameter :: &
         nLat = 20, &
         nLon = 40

    integer :: &
         iLat, iLon

    real(kind=RKIND) :: &
         dLat, &
         dLon

    character(len=7) :: &
         color

    character(len=400) :: &
         stroutvertex, &
         strout, &
         stroutlabel, &
         stroutint

    lnofill = .false.

    fMin = 250.0_RKIND!minval(cellArray(1:mesh % nCells))
    fMax = 300.0_RKIND!maxval(cellArray(1:mesh % nCells))

    !fMin = real(1,RKIND)
    !fMax = real(mesh % nCells,RKIND)

    open(55,file="spherical.txt",action="write")
    open(56,file="clearPlot.gnuplot",action="write")

    do iCell = 1, mesh % nCells

       xc = mesh % xCell % array(iCell)
       yc = mesh % yCell % array(iCell)
       zc = mesh % zCell % array(iCell)

       call perform_rotation(xc, yc, zc)

       if (xc > 0.0_RKIND) then

          fValue = (cellArray(iCell) - fMin) / (fMax - fMin)
          !fValue = (real(iCell,RKIND) - fMin) / (fMax - fMin)

          write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "

          do iVertexOnCell = 1, mesh % nEdgesOnCell % array(iCell)

             iVertex = mesh % verticesOnCell % array(iVertexOnCell, iCell)

             xv = mesh % xVertex % array(iVertex)
             yv = mesh % yVertex % array(iVertex)
             zv = mesh % zVertex % array(iVertex)

             call perform_rotation(xv, yv, zv)

             if (iVertexOnCell == 1) then
                xv0 = xv
                yv0 = yv
                zv0 = zv
             endif

             write(stroutvertex,fmt='(f14.2,a,f14.2,a)') yv, ",", zv, " to "
             strout = trim(strout)//trim(stroutvertex)

          enddo ! iVertexOnCell

          write(stroutvertex,fmt='(f14.2,a,f14.2)') yv0, ",", zv0
          strout = trim(strout)//trim(stroutvertex)
          write(55,*) trim(strout)

          color = matlab_jet(fValue)
          if (lnofill) then
             write(strout,fmt='(a,i5,a,a,a)') "set object ",iObject,' fillstyle empty border lt -1'
          else
             write(strout,fmt='(a,i5,a,a,a)') "set object ",iObject,' fc rgb "', color, '" fillstyle solid'
          endif

          write(55,*) trim(strout)

          write(56,*) "unset object ", iObject
          iObject = iObject + 1

       endif

    enddo ! iCell

    ! now write on latitude/longitude lines
    dLat = 180.0_RKIND / real(nLat,RKIND)
    dLon = 360.0_RKIND / real(nLon,RKIND)

    do iLat = 1, nLat-2
       do iLon = 1, nLon

          lat0 = (real(iLat,RKIND) + 0.5_RKIND) * dLat - 90.0_RKIND
          lon0 = (real(iLon,RKIND) + 0.5_RKIND) * dLat

          lat1 = (real(iLat,RKIND) + 0.0_RKIND) * dLat - 90.0_RKIND
          lon1 = (real(iLon,RKIND) + 0.0_RKIND) * dLat

          lat2 = (real(iLat,RKIND) + 1.0_RKIND) * dLat - 90.0_RKIND
          lon2 = (real(iLon,RKIND) + 0.0_RKIND) * dLat

          lat3 = (real(iLat,RKIND) + 1.0_RKIND) * dLat - 90.0_RKIND
          lon3 = (real(iLon,RKIND) + 1.0_RKIND) * dLat

          lat4 = (real(iLat,RKIND) + 0.0_RKIND) * dLat - 90.0_RKIND
          lon4 = (real(iLon,RKIND) + 1.0_RKIND) * dLat

          call latlon_to_xyz(mesh, x0, y0, z0, lat0, lon0)
          call perform_rotation(x0, y0, z0)

          if (x0 > 0.0_RKIND) then

             call latlon_to_xyz(mesh, x1, y1, z1, lat1, lon1) ; call perform_rotation(x1, y1, z1)
             call latlon_to_xyz(mesh, x2, y2, z2, lat2, lon2) ; call perform_rotation(x2, y2, z2)
             call latlon_to_xyz(mesh, x3, y3, z3, lat3, lon3) ; call perform_rotation(x3, y3, z3)
             call latlon_to_xyz(mesh, x4, y4, z4, lat4, lon4) ; call perform_rotation(x4, y4, z4)
             
             write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "
             
             write(stroutvertex,fmt='(f14.2,a,f14.2,a)') y1, ",", z1, " to " ; strout = trim(strout)//trim(stroutvertex)
             write(stroutvertex,fmt='(f14.2,a,f14.2,a)') y2, ",", z2, " to " ; strout = trim(strout)//trim(stroutvertex)
             write(stroutvertex,fmt='(f14.2,a,f14.2,a)') y3, ",", z3, " to " ; strout = trim(strout)//trim(stroutvertex)
             write(stroutvertex,fmt='(f14.2,a,f14.2,a)') y4, ",", z4, " to " ; strout = trim(strout)//trim(stroutvertex)
             write(stroutvertex,fmt='(f14.2,a,f14.2)')   y1, ",", z1         ; strout = trim(strout)//trim(stroutvertex)

             write(55,*) trim(strout)

             write(strout,fmt='(a,i5,a,a,a)') "set object ",iObject,' fillstyle empty border rgb "gray" lw 0.2'

             write(55,*) trim(strout)

             write(56,*) "unset object ", iObject
             iObject = iObject + 1

          endif

       enddo ! iLon

    enddo ! iLat

    ! prime meridian
    do iLat = 1, nLat

       lat1 = real(iLat-1,RKIND) * dLat - 90.0_RKIND
       lon1 = 0.0_RKIND
       
       lat2 = real(iLat,RKIND) * dLat - 90.0_RKIND
       lon2 = 0.0_RKIND

       call latlon_to_xyz(mesh, x1, y1, z1, lat1, lon1) ; call perform_rotation(x1, y1, z1)
       call latlon_to_xyz(mesh, x2, y2, z2, lat2, lon2) ; call perform_rotation(x2, y2, z2)

       if (0.5_RKIND * (x1 + x2) > 0.0_RKIND) then

          write(strout,fmt='(a,i5,a,f14.2,a,f14.2,a,f14.2,a,f14.2,a)') "set arrow ",iArrow," from ", y1, ", ", z1, " to ", y2, ", ", z2, " nohead lt 4"
          write(55,*) trim(strout)

          write(56,*) "unset arrow ", iArrow
          iArrow = iArrow + 1

       endif

    enddo ! iLat

    ! equator
    do iLon = 1, nLon

       lat1 = 0.0_RKIND
       lon1 = real(iLon-1,RKIND) * dLat
       
       lat2 = 0.0_RKIND
       lon2 = real(iLon,RKIND) * dLat

       call latlon_to_xyz(mesh, x1, y1, z1, lat1, lon1) ; call perform_rotation(x1, y1, z1)
       call latlon_to_xyz(mesh, x2, y2, z2, lat2, lon2) ; call perform_rotation(x2, y2, z2)

       if (0.5_RKIND * (x1 + x2) > 0.0_RKIND) then

          write(strout,fmt='(a,i5,a,f14.2,a,f14.2,a,f14.2,a,f14.2,a)') "set arrow ",iArrow," from ", y1, ", ", z1, " to ", y2, ", ", z2, " nohead lt 3"
          write(55,*) trim(strout)

          write(56,*) "unset arrow ", iArrow
          iArrow = iArrow + 1

       endif

    enddo ! iLon

    ! min/max labels
    write(stroutint,fmt='(e20.8)') fMax
    write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.9'
    write(55,*) trim(stroutlabel)
    write(56,*) "unset label ", iLabel
    iLabel = iLabel + 1

    write(stroutint,fmt='(e20.8)') fMin
    write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.8'
    write(55,*) trim(stroutlabel)
    write(56,*) "unset label ", iLabel
    iLabel = iLabel + 1

    close(55)
    close(56)

    write(strout,fmt='(a,f15.3,a,f15.3,a,f15.3,a,f15.3,a)') "set xrange [",-mesh % sphere_radius*1.1_RKIND,":",mesh % sphere_radius*1.1_RKIND,&
                       "] ; set yrange [",-mesh % sphere_radius*1.1_RKIND,":",mesh % sphere_radius*1.1_RKIND,"]"
    write(*,*) trim(strout)

  end subroutine plot_cells_spherical

  !-------------------------------------------------------------

  subroutine latlon_to_xyz(mesh, x, y, z, lat, lon)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), intent(in) :: &
         lat, lon

    real(kind=RKIND), intent(out) :: &
         x, y, z

    z = mesh % sphere_radius * sin(lat*deg_to_rad)
    x = mesh % sphere_radius * cos(lon*deg_to_rad) * cos(lat*deg_to_rad)
    y = mesh % sphere_radius * sin(lon*deg_to_rad) * cos(lat*deg_to_rad)

  end subroutine latlon_to_xyz

  !-------------------------------------------------------------

  subroutine perform_rotation(x, y, z)

    real(kind=RKIND), intent(inout) :: &
         x, y, z

    real(kind=RKIND), dimension(3,3) :: &
         Rx, Ry, Rz

    real(kind=RKIND), dimension(3) :: &
         v1, v2

    real(kind=RKIND) :: &
         thetax, & ! left-right - leave at zero
         thetay, & ! up-down (positive look at north pole)
         thetaz    ! rotate about axis (positive look west)

    thetax = 0.0_RKIND
    thetay = 50.0_RKIND
    thetaz = 70.0_RKIND

    v1(1) = x ; v1(2) = y ; v1(3) = z

    thetax = thetax * deg_to_rad
    thetay = thetay * deg_to_rad
    thetaz = thetaz * deg_to_rad

    Rx(1,1) = 1.0_RKIND    ; Rx(1,2) = 0.0_RKIND    ; Rx(1,3) = 0.0_RKIND
    Rx(2,1) = 0.0_RKIND    ; Rx(2,2) =  cos(thetax) ; Rx(2,3) = -sin(thetax)
    Rx(3,1) = 0.0_RKIND    ; Rx(3,2) =  sin(thetax) ; Rx(3,3) =  cos(thetax)

    Ry(1,1) =  cos(thetay) ; Ry(1,2) = 0.0_RKIND    ; Ry(1,3) =  sin(thetay)
    Ry(2,1) = 0.0_RKIND    ; Ry(2,2) = 1.0_RKIND    ; Ry(2,3) = 0.0_RKIND
    Ry(3,1) = -sin(thetay) ; Ry(3,2) = 0.0_RKIND    ; Ry(3,3) =  cos(thetay)

    Rz(1,1) =  cos(thetaz) ; Rz(1,2) = -sin(thetaz) ; Rz(1,3) = 0.0_RKIND
    Rz(2,1) =  sin(thetaz) ; Rz(2,2) =  cos(thetaz) ; Rz(2,3) = 0.0_RKIND
    Rz(3,1) = 0.0_RKIND    ; Rz(3,2) = 0.0_RKIND    ; Rz(3,3) = 1.0_RKIND

    v2 = matmul(Rx,matmul(Ry,matmul(Rz,v1)))

    x = v2(1) ; y = v2(2) ; z = v2(3)

  end subroutine perform_rotation

  !-------------------------------------------------------------

  subroutine plot_pop_grid(mesh, cellArray)

    type(mesh_type), intent(in) :: mesh
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         cellArray

    integer :: &
         i, j

    open(55,file='pop_grid.txt',action='write')

    do i = 1, mesh % POP_nx
       do j = 1, mesh % POP_ny

          write(55,*) real(i,RKIND)-0.5_RKIND, real(i,RKIND)+0.5_RKIND, &
                     real(j,RKIND)-0.5_RKIND, real(j,RKIND)+0.5_RKIND, &
                     cellArray(mesh % cellID_ij % array(i,j))

       enddo ! j
    enddo ! i

    close(55)

  end subroutine plot_pop_grid

  !-------------------------------------------------------------

end module cice_testing

