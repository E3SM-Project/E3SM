!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_debug
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module cice_debug

  use mpas_derived_types
  use mpas_pool_routines

  implicit none

  private
  save

  public :: &
       cice_prepare_output_variables, &
       cice_test_case_output, &
       cice_rms_difference, &
       cice_init_ijpop_from_ivertex, &
       cice_plot_spherical_latlon, &
       cice_plot_pop_grid_cells, &
       cice_plot_pop_grid_vertices, &
       cice_time_output_variable_real, &
       cice_time_output_variable_real_2D, &
       cice_time_output_variable_int, &
       cice_time_output_variable_int_2D, &
       cice_plot_vertices_planar, &
       cice_plot_cells_planar, &
       cice_debug_output_array, &
       cice_debug_output_block, &
       cice_debug_output_domain, &
       cice_debug_write_out_fields, &
       cice_writeout_configs, &
       cice_debug_write_out_field, &
       cice_debug_rotate_winds, &
       cice_debug_processor_number, &
       cice_debug_forcing_restart_test

  integer, private :: &
       iObject = 2, &
       iLabel  = 1, &
       iArrow  = 1

  integer, public :: &
       istep = 0

contains

!-----------------------------------------------------------------------
! write out variables routines
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  find_nearest_cell
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  function find_nearest_cell(mesh,x,y) result(iNear)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), intent(in) :: &
         x, & !< Input:
         y    !< Input:

    integer :: iNear

    real(kind=RKIND) :: &
         distance, &
         min_distance

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

    iNear = -1
    min_distance = 1e30_RKIND

    do iCell = 1, nCells

       distance = sqrt((x - xCell(iCell))**2 + (y - yCell(iCell))**2)

       if (distance < min_distance) then

          min_distance = distance

          iNear = iCell

       endif

    enddo ! iCell

  end function find_nearest_cell!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  rotate_ninety
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine rotate_ninety(mesh)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    integer :: &
         iCell, &
         iVertex, &
         iEdge

    real(kind=RKIND) :: x, y

    integer, pointer :: &
         nCells, &
         nVertices, &
         nEdges

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         xVertex, &
         yVertex, &
         xEdge, &
         yEdge

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "nEdges", nEdges)

    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "xEdge", xEdge)
    call MPAS_pool_get_array(mesh, "yEdge", yEdge)

    do iCell = 1, nCells

       x = xCell(iCell)
       y = yCell(iCell)

       xCell(iCell) = y
       yCell(iCell) = -x

    enddo ! iCell

    do iVertex = 1, nVertices

       x = xVertex(iVertex)
       y = yVertex(iVertex)

       xVertex(iVertex) = y
       yVertex(iVertex) = -x

    enddo ! iVertex

    do iEdge = 1, nEdges

       x = xEdge(iEdge)
       y = yEdge(iEdge)

       xEdge(iEdge) = y
       yEdge(iEdge) = -x

    enddo ! iVertex

  end subroutine rotate_ninety!}}}

!-----------------------------------------------------------------------
! Spherical plotting
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_vertices_spherical
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine plot_vertices_spherical(mesh, cellArray, interiorVertex, filename, tag)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         cellArray !< Input:

    integer, dimension(:), intent(in) :: &
         interiorVertex !< Input:

    character(len=*), optional, intent(in) :: &
         filename !< Input:

    integer, intent(in), optional :: &
         tag !< Input:

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

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell

    integer, parameter :: &
         nLat = 20, &
         nLon = 40

    integer :: &
         iLat, iLon

    real(kind=RKIND) :: &
         dLat, &
         dLon

    character(len=400) :: &
         stroutvertex, &
         strout, &
         stroutlabel, &
         stroutint

    integer, pointer :: &
         nVertices, &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    real(kind=RKIND), pointer :: &
         sphere_radius

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         zVertex, &
         xCell, &
         yCell, &
         zCell

    logical, parameter :: &
         lgrid = .false., &
         lequator = .false., &
         lprimemeridian = .false.

    real(kind=RKIND) :: &
         minArray, &
         maxArray

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    iObject = 1

    lnofill = .false.

    call open_filename_tag(55, trim(filename), tag)

    write(55,*) "set xrange [-1:1]"
    write(55,*) "set yrange [-1:1]"
    !write(55,*) "set cbrange [",minval(cellArray(1:nVertices)),":",maxval(cellArray(1:nVertices)),"]"
    !write(55,*) "set cbrange [-0.3:0.3]"

    write(55,*) "set pm3d"
    write(55,*) "set size square"
    write(55,*) "unset key"
    write(55,*) "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0,"//&
         " 4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0 )"

    minArray = 1e30
    maxArray = -1e30

    do iVertex = 1, nVertices

       if (interiorVertex(iVertex) == 1) then

          xc = xVertex(iVertex) / sphere_radius
          yc = yVertex(iVertex) / sphere_radius
          zc = zVertex(iVertex) / sphere_radius

          call perform_rotation(xc, yc, zc)

          if (xc > 0.0_RKIND) then

             minArray = min(minArray,cellArray(iVertex))
             maxArray = max(maxArray,cellArray(iVertex))

             write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "

             do iVertexDegree = 1, vertexDegree

                iCell = cellsOnVertex(iVertexDegree, iVertex)

                xv = xCell(iCell) / sphere_radius
                yv = yCell(iCell) / sphere_radius
                zv = zCell(iCell) / sphere_radius

                call perform_rotation(xv, yv, zv)

                if (iVertexDegree == 1) then
                   xv0 = xv
                   yv0 = yv
                   zv0 = zv
                endif

                write(stroutvertex,fmt='(e14.6,a,e14.6,a)') yv, ", ", zv, " to "
                strout = trim(strout)//trim(stroutvertex)

             enddo ! iVertexDegree

             write(stroutvertex,fmt='(e14.6,a,e14.6)') yv0, ", ", zv0
             strout = trim(strout)//trim(stroutvertex)
             write(55,*) trim(strout)

             if (lnofill) then
                write(strout,fmt='(a,i5,a,a,a)') "set object ",iObject,' fillstyle empty border lt -1'
             else
                write(strout,fmt='(a,i5,a,e14.6,a)') &
                     "set object ",iObject,' fc palette cb ', cellArray(iVertex), ' fillstyle solid'
             endif

             write(55,*) trim(strout)

             iObject = iObject + 1

          endif

       endif ! interiorVertex

    enddo ! iVertex

    write(55,*) "set cbrange [",minArray,":",maxArray,"]"

    ! now write on latitude/longitude lines
    if (lgrid) then
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

                x1 = x1 / sphere_radius ; y1 = y1 / sphere_radius ; z1 = z1 / sphere_radius
                x2 = x2 / sphere_radius ; y2 = y2 / sphere_radius ; z2 = z2 / sphere_radius
                x3 = x3 / sphere_radius ; y3 = y3 / sphere_radius ; z3 = z3 / sphere_radius
                x4 = x4 / sphere_radius ; y4 = y4 / sphere_radius ; z4 = z4 / sphere_radius

                write(stroutvertex,fmt='(e14.6,a,e14.6,a)') y1, ",", z1, " to " ; strout = trim(strout)//trim(stroutvertex)
                write(stroutvertex,fmt='(e14.6,a,e14.6,a)') y2, ",", z2, " to " ; strout = trim(strout)//trim(stroutvertex)
                write(stroutvertex,fmt='(e14.6,a,e14.6,a)') y3, ",", z3, " to " ; strout = trim(strout)//trim(stroutvertex)
                write(stroutvertex,fmt='(e14.6,a,e14.6,a)') y4, ",", z4, " to " ; strout = trim(strout)//trim(stroutvertex)
                write(stroutvertex,fmt='(e14.6,a,e14.6)')   y1, ",", z1         ; strout = trim(strout)//trim(stroutvertex)

                write(55,*) trim(strout)

                write(strout,fmt='(a,i5,a,a,a)') "set object ",iObject,' fillstyle empty border rgb "gray" lw 0.2'

                write(55,*) trim(strout)

                iObject = iObject + 1

             endif

          enddo ! iLon

       enddo ! iLat
    endif ! lgrid

    ! prime meridian
    if (lprimemeridian) then
       do iLat = 1, nLat

          lat1 = real(iLat-1,RKIND) * dLat - 90.0_RKIND
          lon1 = 0.0_RKIND

          lat2 = real(iLat,RKIND) * dLat - 90.0_RKIND
          lon2 = 0.0_RKIND

          call latlon_to_xyz(mesh, x1, y1, z1, lat1, lon1) ; call perform_rotation(x1, y1, z1)
          call latlon_to_xyz(mesh, x2, y2, z2, lat2, lon2) ; call perform_rotation(x2, y2, z2)

          x1 = x1 / sphere_radius ; y1 = y1 / sphere_radius ; z1 = z1 / sphere_radius
          x2 = x2 / sphere_radius ; y2 = y2 / sphere_radius ; z2 = z2 / sphere_radius

          if (0.5_RKIND * (x1 + x2) > 0.0_RKIND) then

             write(strout,fmt='(a,i5,a,e14.6,a,e14.6,a,e14.6,a,e14.6,a)') &
                  "set arrow ",iArrow," from ", y1, ", ", z1, " to ", y2, ", ", z2, " nohead lt 4"
             write(55,*) trim(strout)

             iArrow = iArrow + 1

          endif

       enddo ! iLat
    endif ! lprimemeridian

    ! equator
    if (lequator) then
       do iLon = 1, nLon

          lat1 = 0.0_RKIND
          lon1 = real(iLon-1,RKIND) * dLat

          lat2 = 0.0_RKIND
          lon2 = real(iLon,RKIND) * dLat

          call latlon_to_xyz(mesh, x1, y1, z1, lat1, lon1) ; call perform_rotation(x1, y1, z1)
          call latlon_to_xyz(mesh, x2, y2, z2, lat2, lon2) ; call perform_rotation(x2, y2, z2)

          x1 = x1 / sphere_radius ; y1 = y1 / sphere_radius ; z1 = z1 / sphere_radius
          x2 = x2 / sphere_radius ; y2 = y2 / sphere_radius ; z2 = z2 / sphere_radius

          if (0.5_RKIND * (x1 + x2) > 0.0_RKIND) then

             write(strout,fmt='(a,i5,a,e14.6,a,e14.6,a,e14.6,a,e14.6,a)') &
                  "set arrow ",iArrow," from ", y1, ", ", z1, " to ", y2, ", ", z2, " nohead lt 3"
             write(55,*) trim(strout)

             iArrow = iArrow + 1

          endif

       enddo ! iLon
    endif ! lequator

    ! min/max labels
    !write(stroutint,fmt='(e20.8)') fMax_use
    !write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.9'
    !write(55,*) trim(stroutlabel)
    !iLabel = iLabel + 1

    !write(stroutint,fmt='(e20.8)') fMin_use
    !write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.8'
    !write(55,*) trim(stroutlabel)
    !iLabel = iLabel + 1

    write(55,*) "plot -10"

    close(55)

  end subroutine plot_vertices_spherical!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_cells_spherical
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine plot_cells_spherical(mesh, cellArray, filename, tag)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         cellArray !< Input:

    character(len=*), optional, intent(in) :: &
         filename !< Input:

    integer, intent(in), optional :: &
         tag !< Input:

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

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    real(kind=RKIND), pointer :: &
         sphere_radius

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         zCell, &
         xVertex, &
         yVertex, &
         zVertex

    logical, parameter :: &
         lgrid = .false., &
         lequator = .false., &
         lprimemeridian = .false.

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)

    iObject = 1

    lnofill = .false.

    call open_filename_tag(55, trim(filename), tag)

    write(55,*) "set xrange [-1:1]"
    write(55,*) "set yrange [-1:1]"
    write(55,*) "set cbrange [",minval(cellArray(1:nCells)),":",maxval(cellArray(1:nCells)),"]"
    !write(55,*) "set cbrange [-0.3:0.3]"

    write(55,*) "set pm3d"
    write(55,*) "set size square"
    write(55,*) "unset key"
    write(55,*) "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0, "//&
         "4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0 )"

    do iCell = 1, nCells

       xc = xCell(iCell) / sphere_radius
       yc = yCell(iCell) / sphere_radius
       zc = zCell(iCell) / sphere_radius

       call perform_rotation(xc, yc, zc)

       if (xc > 0.0_RKIND) then

          write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "

          do iVertexOnCell = 1, nEdgesOnCell(iCell)

             iVertex = verticesOnCell(iVertexOnCell, iCell)

             xv = xVertex(iVertex) / sphere_radius
             yv = yVertex(iVertex) / sphere_radius
             zv = zVertex(iVertex) / sphere_radius

             call perform_rotation(xv, yv, zv)

             if (iVertexOnCell == 1) then
                xv0 = xv
                yv0 = yv
                zv0 = zv
             endif

             write(stroutvertex,fmt='(e14.6,a,e14.6,a)') yv, ",", zv, " to "
             strout = trim(strout)//trim(stroutvertex)

          enddo ! iVertexOnCell

          write(stroutvertex,fmt='(e14.6,a,e14.6)') yv0, ",", zv0
          strout = trim(strout)//trim(stroutvertex)
          write(55,*) trim(strout)

          if (lnofill) then
             write(strout,fmt='(a,i5,a,a,a)') "set object ",iObject,' fillstyle empty border lt -1'
          else
             write(strout,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', cellArray(iCell), ' fillstyle solid'
          endif

          write(55,*) trim(strout)

          iObject = iObject + 1

       endif

    enddo ! iCell

    ! now write on latitude/longitude lines
    if (lgrid) then
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

                x1 = x1 / sphere_radius ; y1 = y1 / sphere_radius ; z1 = z1 / sphere_radius
                x2 = x2 / sphere_radius ; y2 = y2 / sphere_radius ; z2 = z2 / sphere_radius
                x3 = x3 / sphere_radius ; y3 = y3 / sphere_radius ; z3 = z3 / sphere_radius
                x4 = x4 / sphere_radius ; y4 = y4 / sphere_radius ; z4 = z4 / sphere_radius

                write(stroutvertex,fmt='(f14.2,a,f14.2,a)') y1, ",", z1, " to " ; strout = trim(strout)//trim(stroutvertex)
                write(stroutvertex,fmt='(f14.2,a,f14.2,a)') y2, ",", z2, " to " ; strout = trim(strout)//trim(stroutvertex)
                write(stroutvertex,fmt='(f14.2,a,f14.2,a)') y3, ",", z3, " to " ; strout = trim(strout)//trim(stroutvertex)
                write(stroutvertex,fmt='(f14.2,a,f14.2,a)') y4, ",", z4, " to " ; strout = trim(strout)//trim(stroutvertex)
                write(stroutvertex,fmt='(f14.2,a,f14.2)')   y1, ",", z1         ; strout = trim(strout)//trim(stroutvertex)

                write(55,*) trim(strout)

                write(strout,fmt='(a,i5,a,a,a)') "set object ",iObject,' fillstyle empty border rgb "gray" lw 0.2'

                write(55,*) trim(strout)

                iObject = iObject + 1

             endif

          enddo ! iLon

       enddo ! iLat
    endif ! lgrid

    ! prime meridian
    if (lprimemeridian) then
       do iLat = 1, nLat

          lat1 = real(iLat-1,RKIND) * dLat - 90.0_RKIND
          lon1 = 0.0_RKIND

          lat2 = real(iLat,RKIND) * dLat - 90.0_RKIND
          lon2 = 0.0_RKIND

          call latlon_to_xyz(mesh, x1, y1, z1, lat1, lon1) ; call perform_rotation(x1, y1, z1)
          call latlon_to_xyz(mesh, x2, y2, z2, lat2, lon2) ; call perform_rotation(x2, y2, z2)

          x1 = x1 / sphere_radius ; y1 = y1 / sphere_radius ; z1 = z1 / sphere_radius
          x2 = x2 / sphere_radius ; y2 = y2 / sphere_radius ; z2 = z2 / sphere_radius

          if (0.5_RKIND * (x1 + x2) > 0.0_RKIND) then

             write(strout,fmt='(a,i5,a,f14.2,a,f14.2,a,f14.2,a,f14.2,a)') &
                  "set arrow ",iArrow," from ", y1, ", ", z1, " to ", y2, ", ", z2, " nohead lt 4"
             write(55,*) trim(strout)

             iArrow = iArrow + 1

          endif

       enddo ! iLat
    endif ! lprimemeridian

    ! equator
    if (lequator) then
       do iLon = 1, nLon

          lat1 = 0.0_RKIND
          lon1 = real(iLon-1,RKIND) * dLat

          lat2 = 0.0_RKIND
          lon2 = real(iLon,RKIND) * dLat

          call latlon_to_xyz(mesh, x1, y1, z1, lat1, lon1) ; call perform_rotation(x1, y1, z1)
          call latlon_to_xyz(mesh, x2, y2, z2, lat2, lon2) ; call perform_rotation(x2, y2, z2)

          x1 = x1 / sphere_radius ; y1 = y1 / sphere_radius ; z1 = z1 / sphere_radius
          x2 = x2 / sphere_radius ; y2 = y2 / sphere_radius ; z2 = z2 / sphere_radius

          if (0.5_RKIND * (x1 + x2) > 0.0_RKIND) then

             write(strout,fmt='(a,i5,a,f14.2,a,f14.2,a,f14.2,a,f14.2,a)') &
                  "set arrow ",iArrow," from ", y1, ", ", z1, " to ", y2, ", ", z2, " nohead lt 3"
             write(55,*) trim(strout)

             iArrow = iArrow + 1

          endif

       enddo ! iLon
    endif ! lequator

    ! min/max labels
    !write(stroutint,fmt='(e20.8)') fMax_use
    !write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.9'
    !write(55,*) trim(stroutlabel)
    !write(56,*) "unset label ", iLabel
    !iLabel = iLabel + 1

    !write(stroutint,fmt='(e20.8)') fMin_use
    !write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.8'
    !write(55,*) trim(stroutlabel)
    !write(56,*) "unset label ", iLabel
    !iLabel = iLabel + 1

    write(55,*) "plot -10"

    close(55)

  end subroutine plot_cells_spherical!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  latlon_to_xyz
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine latlon_to_xyz(mesh, x, y, z, lat, lon)!{{{

    use cice_constants, only: ciceDegreesToRadians

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), intent(in) :: &
         lat, & !< Input:
         lon    !< Input:

    real(kind=RKIND), intent(out) :: &
         x, & !< Output:
         y, & !< Output:
         z    !< Output:

    real(kind=RKIND), pointer :: &
         sphere_radius

    ! init variables
    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    z = sphere_radius * sin(lat*ciceDegreesToRadians)
    x = sphere_radius * cos(lon*ciceDegreesToRadians) * cos(lat*ciceDegreesToRadians)
    y = sphere_radius * sin(lon*ciceDegreesToRadians) * cos(lat*ciceDegreesToRadians)

  end subroutine latlon_to_xyz!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  perform_rotation
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine perform_rotation(x, y, z)!{{{

    use cice_constants, only: &
         ciceDegreesToRadians

    real(kind=RKIND), intent(inout) :: &
         x, & !< Input/Output:
         y, & !< Input/Output:
         z    !< Input/Output:

    real(kind=RKIND), dimension(3,3) :: &
         Rx, Ry, Rz

    real(kind=RKIND), dimension(3) :: &
         v1, v2

    real(kind=RKIND) :: &
         thetax, & ! left-right - leave at zero
         thetay, & ! up-down (positive look at north pole)
         thetaz    ! rotate about axis (positive look west)

    thetax = 0.0_RKIND
    thetay = -90.0_RKIND
    thetaz = 90.0_RKIND

    v1(1) = x ; v1(2) = y ; v1(3) = z

    thetax = thetax * ciceDegreesToRadians
    thetay = thetay * ciceDegreesToRadians
    thetaz = thetaz * ciceDegreesToRadians

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

  end subroutine perform_rotation!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_pop_grid_cells
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_plot_pop_grid_cells(mesh, cellArray, filename, tag)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         cellArray !< Input:

    character(len=*), optional, intent(in) :: &
         filename !< Input:

    integer, optional, intent(in) :: &
         tag !< Input:

    integer :: &
         i, j, iCell

    logical :: lfound

    character(len=5) :: &
         strtag

    character(len=200) :: &
         filename_use

    integer, pointer :: &
         POP_nx, &
         POP_ny, &
         nCells

    integer, dimension(:), pointer :: &
         POPindxi, &
         POPindxj

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "POP_nx", POP_nx)
    call MPAS_pool_get_dimension(mesh, "POP_ny", POP_ny)

    call MPAS_pool_get_array(mesh, "POPindxi", POPindxi)
    call MPAS_pool_get_array(mesh, "POPindxj", POPindxj)

    if (present(filename)) then

       if (present(tag)) then
          write(strtag,fmt='(i5.5)') tag
          filename_use = filename(1:len(trim(filename))-4)//"_"//trim(strtag)//filename(len(trim(filename))-3:)
          open(55,file='history/'//trim(filename_use),action='write')
       else
          open(55,file='history/'//trim(filename),action='write')
       endif

    else
       open(55,file='history/'//'pop_grid.txt',action='write')
    endif

    do i = 1, POP_nx
       do j = 1, POP_ny

          lfound = .false.

          do iCell = 1, nCells

             if (POPindxi(iCell) == i .and. POPindxj(iCell) == j) then

                lfound = .true.

                write(55,*) real(i,RKIND)-0.5_RKIND, real(i,RKIND)+0.5_RKIND, &
                            real(j,RKIND)-0.5_RKIND, real(j,RKIND)+0.5_RKIND, &
                            cellArray(iCell)

                exit

             endif

          enddo ! iCell

          !if (.not. lfound) write(55,*) real(i,RKIND)-0.5_RKIND, real(i,RKIND)+0.5_RKIND, &
          !                              real(j,RKIND)-0.5_RKIND, real(j,RKIND)+0.5_RKIND, &
          !                              -99.0_RKIND

       enddo ! j
    enddo ! i

    close(55)

    !do iCell = 1, nCells

    !   write(*,*) iCell, POPindxi(iCell), POPindxj(iCell)

    !enddo ! iCell

    !stop

  end subroutine cice_plot_pop_grid_cells!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_pop_grid_vertices
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_plot_pop_grid_vertices(mesh, vertexArray, interiorVertex, filename, tag)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         vertexArray !< Input:

    integer, dimension(:), intent(in) :: &
         interiorVertex !< Input:

    character(len=*), optional, intent(in) :: &
         filename !< Input:

    integer, optional, intent(in) :: &
         tag !< Input:

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell, &
         i, &
         j

    character(len=5) :: &
         strtag

    character(len=200) :: &
         filename_use

    real(kind=RKIND), dimension(:,:,:), allocatable :: &
         arrayOrdered

    integer, pointer :: &
         POP_nx, &
         POP_ny, &
         nVertices

    integer, dimension(:), pointer :: &
         POPindxiv, &
         POPindxjv

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "POP_nx", POP_nx)
    call MPAS_pool_get_dimension(mesh, "POP_ny", POP_ny)

    call MPAS_pool_get_array(mesh, "POPindxiv", POPindxiv)
    call MPAS_pool_get_array(mesh, "POPindxjv", POPindxjv)

    allocate(arrayOrdered(POP_nx,POP_ny,2))
    arrayOrdered = 0.0_RKIND

    if (present(filename)) then

       if (present(tag)) then
          write(strtag,fmt='(i5.5)') tag
          filename_use = filename(1:len(trim(filename))-4)//"_"//trim(strtag)//filename(len(trim(filename))-3:)
          open(55,file='history/'//trim(filename_use),action='write')
       else
          open(55,file='history/'//trim(filename),action='write')
       endif

    else
       open(55,file='history/'//'pop_grid.txt',action='write')
    endif

    do iVertex = 1, nVertices

       if (interiorVertex(iVertex) == 1) then

          i = POPindxiv(iVertex)
          j = POPindxjv(iVertex)

          if (i > 0 .and. i <= size(arrayOrdered,1) .and. j > 0 .and. j <= size(arrayOrdered,2)) then

             arrayOrdered(i,j,1) = 1.0_RKIND
             arrayOrdered(i,j,2) = vertexArray(iVertex)

          endif

       endif ! interiorVertex

    enddo ! iVertex

    do i = 1, POP_nx
       do j = 1, POP_ny

          if (arrayOrdered(i,j,1) == 1.0_RKIND) then

             write(55,*) &
                  real(i,RKIND)-0.5_RKIND, real(i,RKIND)+0.5_RKIND, &
                  real(j,RKIND)-0.5_RKIND, real(j,RKIND)+0.5_RKIND, &
                  arrayOrdered(i,j,2)

          endif

       enddo! j
    enddo ! i

    close(55)

    deallocate(arrayOrdered)

  end subroutine cice_plot_pop_grid_vertices!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_ijpop_from_ivertex
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_init_ijpop_from_ivertex(mesh)!{{{

    type(MPAS_pool_type), intent(inout) :: &
         mesh !< Input/Output:

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell, &
         i, j, &
         imin, jmin, &
         imax, jmax, &
         POP_nx_2

    logical :: l_boundary

    integer, pointer :: &
         nCells, &
         nVertices, &
         vertexDegree, &
         POP_nx

    integer, dimension(:), pointer :: &
         POPindxi, &
         POPindxj, &
         POPindxiv, &
         POPindxjv

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    real(kind=RKIND), dimension(:), pointer :: &
         areaTriangle, &
         areaCell, &
         latVertex

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_dimension(mesh, "POP_nx", POP_nx)

    call MPAS_pool_get_array(mesh, "POPindxi", POPindxi)
    call MPAS_pool_get_array(mesh, "POPindxj", POPindxj)
    call MPAS_pool_get_array(mesh, "POPindxiv", POPindxiv)
    call MPAS_pool_get_array(mesh, "POPindxjv", POPindxjv)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "areaTriangle", areaTriangle)
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)
    call MPAS_pool_get_array(mesh, "latVertex", latVertex)

    do iVertex = 1, nVertices

       imin = 1000000000
       jmin = 1000000000

       imax = -1000000000
       jmax = -1000000000

       ! first we find the minimum and maximum POP i value of a cell surrounding the vertex point
       do iVertexDegree = 1, vertexDegree

          iCell = cellsOnVertex(iVertexDegree,iVertex)

          i = POPindxi(iCell)
          j = POPindxj(iCell)

          imin = min(imin,i)
          jmin = min(jmin,j)

          imax = max(imax,i)
          jmax = max(jmax,j)

       enddo ! iVertexDegree

       ! decide if at border
       POP_nx_2 = nint(real(POP_nx,RKIND) / 2.0_RKIND)

       l_boundary = .false.
       if (imin < POP_nx_2 .and. imax > POP_nx_2) l_boundary = .true.

       do iVertexDegree = 1, vertexDegree

          iCell = cellsOnVertex(iVertexDegree,iVertex)

          i = POPindxi(iCell)
          j = POPindxj(iCell)

          if (l_boundary .and. i > POP_nx_2) i = i - POP_nx
          imin = min(imin,i)

          jmin = min(jmin,j)

       enddo ! iVertexDegree

       if (l_boundary .and. imin < 1) imin = imin + POP_nx

       POPindxiv(iVertex) = imin
       POPindxjv(iVertex) = jmin

    enddo ! iVertex

    !open(11,file="vertexareas.txt")
    !do iVertex = 1, nVertices
    !   write(11,*) POPindxiv(iVertex), POPindxjv(iVertex), &
    !        areaTriangle(iVertex)
    !enddo
    !close(11)

    !open(11,file="cellareas.txt")
    !do iCell = 1, nCells
    !   write(11,*) POPindxi(iCell), POPindxj(iCell), &
    !        areaCell(iCell)
    !enddo
    !close(11)

    !open(11,file="vertexlat.txt")
    !do iVertex = 1, nVertices
    !   write(11,*) POPindxiv(iVertex), POPindxjv(iVertex), &
    !        latVertex(iVertex)
    !enddo
    !close(11)

  end subroutine cice_init_ijpop_from_ivertex!}}}

!-----------------------------------------------------------------------
! RMS comparison
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_rms_difference
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_rms_difference(array1, array2, mask, n, rms)!{{{

    real(kind=RKIND), dimension(:), intent(in) :: &
         array1, & !< Input:
         array2    !< Input:

    integer, dimension(:), intent(in) :: &
         mask !< Input:

    integer, intent(in) :: &
         n !< Input:

    real(kind=RKIND), intent(out) :: &
         rms !< Output:

    integer :: &
         i, &
         num

    rms = 0.0_RKIND
    num = 0

    do i = 1, n

       if (mask(i) == 1) then

          rms = rms + (array1(i) - array2(i))**2
          num = num + 1

       endif

    enddo ! i

    rms = sqrt(rms / real(num, RKIND))

  end subroutine cice_rms_difference!}}}

!-----------------------------------------------------------------------
! planar plots
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_plot_vertices_planar
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_plot_vertices_planar(mesh, array, interiorVertex, filename)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         array !< Input:

    integer, dimension(:), intent(in) :: &
         interiorVertex !< Input:

    character(len=*), intent(in) :: &
         filename !< Input:

    integer, pointer :: &
         nVertices, &
         vertexDegree

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    integer :: &
         iCell, &
         iVertex, &
         iVertexDegree, &
         iObject

    real(kind=RKIND) :: &
         xv, yv, &
         xv0, yv0

    character(len=1000) :: &
         strout, &
         stroutvertex

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)

    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)

    iObject = 1

    open(55,file=trim(filename))

    write(55,*) "#", minval(array), maxval(array)
    write(55,*) "set cbrange [",minval(array),":",maxval(array),"]"
    write(55,*) "set xrange  [",minval(xCell),":",maxval(xCell),"]"
    write(55,*) "set yrange  [",minval(yCell),":",maxval(yCell),"]"

    write(55,*) "set pm3d"
    write(55,*) "set size square"
    write(55,*) "unset key"
    write(55,*) "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0, "//&
         "4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0 )"

    do iVertex = 1, nVertices

       if (interiorVertex(iVertex) == 1) then

          write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "

          do iVertexDegree = 1, vertexDegree

             iCell = cellsOnVertex(iVertexDegree, iVertex)

             xv = xCell(iCell)
             yv = yCell(iCell)

             if (iVertexDegree == 1) then
                xv0 = xv
                yv0 = yv
             endif

             write(stroutvertex,fmt='(f14.2,a,f14.2,a)') xv, ",", yv, " to "
             strout = trim(strout)//trim(stroutvertex)

          enddo ! iVertexDegree

          write(stroutvertex,fmt='(f14.2,a,f14.2)') xv0, ",", yv0
          strout = trim(strout)//trim(stroutvertex)
          write(55,*) trim(strout)

          write(strout,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', array(iVertex), ' fillstyle solid'
          write(55,*) trim(strout)

          iObject = iObject + 1

       endif

    enddo ! iVertex

    write(55,*) "plot ", minval(yCell)

    close(55)

  end subroutine cice_plot_vertices_planar!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_plot_cells_planar
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_plot_cells_planar(mesh, array, filename, tag)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         array !< Input:

    character(len=*), intent(in) :: &
         filename !< Input:

    integer, intent(in), optional :: &
         tag !< Input:

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    integer :: &
         iCell, &
         iVertex, &
         iVertexOnCell, &
         iObject

    real(kind=RKIND) :: &
         xv, yv, &
         xv0, yv0

    character(len=1000) :: &
         strout, &
         stroutvertex

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)

    iObject = 1

    call open_filename_tag(55, filename, tag)

    write(55,*) "#", minval(array), maxval(array)
    write(55,*) "set cbrange [",minval(array),":",maxval(array),"]"
    write(55,*) "set xrange  [",minval(xVertex),":",maxval(xVertex),"]"
    write(55,*) "set yrange  [",minval(yVertex),":",maxval(yVertex),"]"

    write(55,*) "set pm3d"
    write(55,*) "set size square"
    write(55,*) "unset key"
    write(55,*) "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0, "//&
         "4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0 )"

    do iCell = 1, nCells

       write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          iVertex = verticesOnCell(iVertexOnCell, iCell)

          xv = xVertex(iVertex)
          yv = yVertex(iVertex)

          if (iVertexOnCell == 1) then
             xv0 = xv
             yv0 = yv
          endif

          write(stroutvertex,fmt='(f14.2,a,f14.2,a)') xv, ",", yv, " to "
          strout = trim(strout)//trim(stroutvertex)

       enddo ! iVertexOnCell

       write(stroutvertex,fmt='(f14.2,a,f14.2)') xv0, ",", yv0
       strout = trim(strout)//trim(stroutvertex)
       write(55,*) trim(strout)

       write(strout,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', array(iCell), ' fillstyle solid'
       write(55,*) trim(strout)

       iObject = iObject + 1

    enddo ! iVertex

    write(55,*) "plot ", minval(yVertex)

    close(55)

  end subroutine cice_plot_cells_planar!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  open_filename_tag
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine open_filename_tag(unit, filename, tag)!{{{

    integer, intent(in) :: &
         unit !< Input:

    character(len=*), intent(in) :: &
         filename !< Input:

    integer, optional, intent(in) :: &
         tag !< Input:

    character(len=2000) :: &
         strtag, &
         filename_use

    if (present(tag)) then
       write(strtag,fmt='(i5.5)') tag
       filename_use = filename(1:len(trim(filename))-4)//"_"//trim(strtag)//filename(len(trim(filename))-3:)
       open(unit,file=trim(filename_use),action='write')
    else
       open(unit,file=trim(filename),action='write')
    endif

  end subroutine open_filename_tag!}}}

!-----------------------------------------------------------------------
! plot spherical on lat lon
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_plot_spherical_latlon
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_plot_spherical_latlon(mesh, array, filename, tag)!{{{

    use cice_constants, only: &
         pii

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         array !< Input:

    character(len=*), intent(in) :: &
         filename !< Input:

    integer, optional, intent(in) :: &
         tag !< Input:

    integer, pointer :: &
         nCells, &
         nVertices, &
         vertexDegree

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         pointVertexOnPoint

    real(kind=RKIND), dimension(:), pointer :: &
         lonPointVertex, &
         latPointVertex, &
         lonPoint, &
         latPoint

    integer, dimension(:), allocatable :: &
         nPointVerticesOnPoint, &
         plotPoint

    integer :: &
         iPoint, &
         iPointVertexOnPoint, &
         iPointVertex, &
         nPoints, &
         array_size

    real(kind=RKIND) :: &
         plotLonPointVertex, &
         plotLatPointVertex, &
         plotLonPointVertex0, &
         plotLatPointVertex0

    character(len=2000) :: &
         strout, &
         stroutvertex

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    if (size(array) == nCells+1 .or. size(array) == nCells) then

       array_size = nCells

       nPoints = nCells
       call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
       allocate(nPointVerticesOnPoint(nCells))
       nPointVerticesOnPoint(:) = nEdgesOnCell(1:nCells)
       call MPAS_pool_get_array(mesh, "verticesOnCell", pointVertexOnPoint)
       call MPAS_pool_get_array(mesh, "lonVertex", lonPointVertex)
       call MPAS_pool_get_array(mesh, "latVertex", latPointVertex)
       call MPAS_pool_get_array(mesh, "lonCell", lonPoint)
       call MPAS_pool_get_array(mesh, "latCell", latPoint)
       allocate(plotPoint(nCells))
       plotPoint = 1

    else if (size(array) == nVertices+1 .or. size(array) == nVertices) then

       array_size = nVertices

       nPoints = nVertices
       call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)
       allocate(nPointVerticesOnPoint(nVertices))
       nPointVerticesOnPoint(:) = vertexDegree
       call MPAS_pool_get_array(mesh, "cellsOnVertex", pointVertexOnPoint)
       call MPAS_pool_get_array(mesh, "lonCell", lonPointVertex)
       call MPAS_pool_get_array(mesh, "latCell", latPointVertex)
       call MPAS_pool_get_array(mesh, "lonVertex", lonPoint)
       call MPAS_pool_get_array(mesh, "latVertex", latPoint)
       allocate(plotPoint(nVertices))
       do iPoint = 1, nVertices
          plotPoint(iPoint) = 1
          do iPointVertexOnPoint = 1, vertexDegree
             if (pointVertexOnPoint(iPointVertexOnPoint,iPoint) > nCells) then
                plotPoint(iPoint) = 0
             endif
          enddo ! iPointVertexOnPoint
       enddo ! iPoint

    else
       write(*,*) "cice_plot_spherical_latlon: size not supported: ", size(array)
    endif

    call open_filename_tag(55, filename, tag)

    write(55,*) "set xrange [0:6.283]"
    write(55,*) "set yrange [-1.571:1.571]"
    write(55,*) "set cbrange [",minval(array(1:array_size)),":",maxval(array(1:array_size)),"]"

    write(55,*) "set pm3d"
    write(55,*) "set size square"
    write(55,*) "unset key"
    write(55,*) "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0, "//&
         "4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0 )"

    iObject = 1

    do iPoint = 1, nPoints

       if (plotPoint(iPoint) == 1) then

          write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "

          do iPointVertexOnPoint = 1, nPointVerticesOnPoint(iPoint)

             iPointVertex = pointVertexOnPoint(iPointVertexOnPoint,iPoint)

             plotLonPointVertex = lonPointVertex(iPointVertex)
             plotLatPointVertex = latPointVertex(iPointVertex)

             if (lonPoint(iPoint) <= pii) then

                if (plotLonPointVertex >= 0.9_RKIND * 2.0_RKIND * pii) then
                   plotLonPointVertex = plotLonPointVertex - 2.0_RKIND * pii
                endif

             endif

             if (lonPoint(iPoint) > pii) then

                if (plotLonPointVertex <= 0.1_RKIND * 2.0_RKIND * pii) then
                   plotLonPointVertex = plotLonPointVertex + 2.0_RKIND * pii
                endif

             endif

             if (iPointVertexOnPoint == 1) then

                plotLonPointVertex0 = plotLonPointVertex
                plotLatPointVertex0 = plotLatPointVertex

             endif

             write(stroutvertex,fmt='(f14.2,a,f14.2,a)') plotLonPointVertex, ",", plotLatPointVertex, " to "
             strout = trim(strout)//trim(stroutvertex)

          enddo ! iPointVertexOnPoint

          write(stroutvertex,fmt='(f14.2,a,f14.2)') plotLonPointVertex0, ",", plotLatPointVertex0
          strout = trim(strout)//trim(stroutvertex)
          write(55,*) trim(strout)

          write(strout,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', array(iPoint), ' fillstyle solid'
          write(55,*) trim(strout)

          iObject = iObject + 1

       endif ! plotPoint

    enddo ! iPoint

    deallocate(nPointVerticesOnPoint)
    deallocate(plotPoint)

    write(55,*) "plot -10"

    close(55)

  end subroutine cice_plot_spherical_latlon!}}}

!-----------------------------------------------------------------------
! strip_plot
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  strip_plot
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine strip_plot(array, filename, tag)!{{{

    real(kind=RKIND), dimension(:), intent(in) :: &
         array !< Input:

    character(len=*), intent(in) :: &
         filename !< Input:

    integer, intent(in) :: &
         tag !< Input:

    integer :: &
         iCell

    integer, parameter :: &
         iCellPoint = 81*40

    real(kind=RKIND) :: &
         diff

    call open_filename_tag(55, filename, tag)

    do iCell = iCellPoint - 5, iCellPoint + 45

       diff = array(iCell) - 1.0_RKIND
       if (array(iCell) == 0.0_RKIND) diff = 0

       write(55,*) iCell, array(iCell), diff

    enddo ! iCell

    close(55)

  end subroutine strip_plot!}}}

!-----------------------------------------------------------------------
! test case output
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_test_case_output
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_test_case_output(domain, itimestep)!{{{

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    integer, intent(in) :: &
         itimestep !< Input:

    type(block_type), pointer :: &
         block

    type (MPAS_pool_type), pointer :: &
         configs

    logical, pointer :: &
         config_test_case_diag

    character(len=strKIND), pointer :: &
         config_test_case_diag_type

    block => domain % blocklist
    do while (associated(block))

       configs => block % configs

       call MPAS_pool_get_config(configs, "config_test_case_diag", config_test_case_diag)

       if (config_test_case_diag) then

          call MPAS_pool_get_config(configs, "config_test_case_diag_type", config_test_case_diag_type)

          select case (trim(config_test_case_diag_type))
          case ("global_velocity")
             call test_case_output_global_velocity(block, itimestep)
          case ("square_advection")
             call test_case_output_square_advection(block, itimestep)
          case ("advection")
             !call test_case_output_advection(block, itimestep)
          case default
             write(*,*) "Test case diagnostics not implemented!", trim(config_test_case_diag_type)
             stop
          end select

       endif

       block => block % next
    end do

  end subroutine cice_test_case_output!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  test_case_output_global_velocity
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine test_case_output_global_velocity(block, itimestep)!{{{

    use cice_mesh, only: &
         cice_latlon_vector_rotation_backward

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    integer, intent(in) :: &
         itimestep !< Input:

    type (MPAS_pool_type), pointer :: &
         configs, &
         mesh, &
         velocity_solver, &
         boundary, &
         atmos_forcing, &
         tracers_aggregate

    integer, dimension(:), pointer :: &
         interiorVertex

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         airStressVertexU, &
         airStressVertexV, &
         airStressCellU, &
         airStressCellV, &
         totalMassVertexfVertex, &
         uAirVelocity, &
         vAirVelocity, &
         stressDivergenceU, &
         stressDivergenceV, &
         latVertex, &
         lonVertex, &
         xVertex, &
         yVertex, &
         zVertex, &
         latCell, &
         lonCell, &
         xCell, &
         yCell, &
         zCell, &
         iceAreaCell

    real(kind=RKIND), dimension(:), allocatable :: &
         uVelocityGeographical, &
         vVelocityGeographical, &
         airStressVertexUGeographical, &
         airStressVertexVGeographical, &
         airStressCellUGeographical, &
         airStressCellVGeographical, &
         uAirVelocityGeographical, &
         vAirVelocityGeographical, &
         stressDivergenceUGeographical, &
         stressDivergenceVGeographical

    integer, pointer :: &
         nVertices, &
         nCells

    logical, pointer :: &
         config_rotate_cartesian_grid

    real(kind=RKIND), pointer :: &
         sphere_radius

    integer :: &
         iVertex, &
         iCell

    configs => block % configs
    call MPAS_pool_get_config(configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmos_forcing)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)

    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "latVertex", latVertex)
    call MPAS_pool_get_array(mesh, "lonVertex", lonVertex)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)

    call MPAS_pool_get_array(mesh, "latCell", latCell)
    call MPAS_pool_get_array(mesh, "lonCell", lonCell)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)
    call MPAS_pool_get_array(velocity_solver, "airStressVertexU", airStressVertexU)
    call MPAS_pool_get_array(velocity_solver, "airStressVertexV", airStressVertexV)
    call MPAS_pool_get_array(velocity_solver, "airStressCellU", airStressCellU)
    call MPAS_pool_get_array(velocity_solver, "airStressCellV", airStressCellV)
    call MPAS_pool_get_array(velocity_solver, "totalMassVertexfVertex", totalMassVertexfVertex)
    call MPAS_pool_get_array(velocity_solver, "stressDivergenceU", stressDivergenceU)
    call MPAS_pool_get_array(velocity_solver, "stressDivergenceV", stressDivergenceV)

    call MPAS_pool_get_array(atmos_forcing, "uAirVelocity", uAirVelocity)
    call MPAS_pool_get_array(atmos_forcing, "vAirVelocity", vAirVelocity)

    call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)

    allocate(uVelocityGeographical(nVertices))
    allocate(vVelocityGeographical(nVertices))
    allocate(airStressVertexUGeographical(nVertices))
    allocate(airStressVertexVGeographical(nVertices))
    allocate(airStressCellUGeographical(nCells))
    allocate(airStressCellVGeographical(nCells))
    allocate(uAirVelocityGeographical(nCells))
    allocate(vAirVelocityGeographical(nCells))
    allocate(stressDivergenceUGeographical(nVertices))
    allocate(stressDivergenceVGeographical(nVertices))

    do iVertex = 1, nVertices
       call cice_latlon_vector_rotation_backward(&
            uVelocityGeographical(iVertex), &
            vVelocityGeographical(iVertex), &
            uVelocity(iVertex), &
            vVelocity(iVertex), &
            latVertex(iVertex), &
            lonVertex(iVertex), &
            xVertex(iVertex), &
            yVertex(iVertex), &
            zVertex(iVertex), &
            sphere_radius, &
            config_rotate_cartesian_grid)
       call cice_latlon_vector_rotation_backward(&
            airStressVertexUGeographical(iVertex), &
            airStressVertexVGeographical(iVertex), &
            airStressVertexU(iVertex), &
            airStressVertexV(iVertex), &
            latVertex(iVertex), &
            lonVertex(iVertex), &
            xVertex(iVertex), &
            yVertex(iVertex), &
            zVertex(iVertex), &
            sphere_radius, &
            config_rotate_cartesian_grid)
       call cice_latlon_vector_rotation_backward(&
            stressDivergenceUGeographical(iVertex), &
            stressDivergenceVGeographical(iVertex), &
            stressDivergenceU(iVertex), &
            stressDivergenceV(iVertex), &
            latVertex(iVertex), &
            lonVertex(iVertex), &
            xVertex(iVertex), &
            yVertex(iVertex), &
            zVertex(iVertex), &
            sphere_radius, &
            config_rotate_cartesian_grid)
    enddo ! iVertex

    do iCell = 1, nCells
       call cice_latlon_vector_rotation_backward(&
            uAirVelocityGeographical(iCell), &
            vAirVelocityGeographical(iCell), &
            uAirVelocity(iCell), &
            vAirVelocity(iCell), &
            latCell(iCell), &
            lonCell(iCell), &
            xCell(iCell), &
            yCell(iCell), &
            zCell(iCell), &
            sphere_radius, &
            config_rotate_cartesian_grid)
       call cice_latlon_vector_rotation_backward(&
            airStressCellUGeographical(iCell), &
            airStressCellVGeographical(iCell), &
            airStressCellU(iCell), &
            airStressCellV(iCell), &
            latCell(iCell), &
            lonCell(iCell), &
            xCell(iCell), &
            yCell(iCell), &
            zCell(iCell), &
            sphere_radius, &
            config_rotate_cartesian_grid)
    enddo ! iCell
    !uAirVelocityGeographical = uAirVelocity(1:nCells)
    !vAirVelocityGeographical = vAirVelocity(1:nCells)

    ! air velocity
    call plot_cells_spherical(mesh, uAirVelocityGeographical, "uAirVelocity.txt", itimestep)
    call plot_cells_spherical(mesh, vAirVelocityGeographical, "vAirVelocity.txt", itimestep)

    write(*,*) "uAirVelocity: ", sum(abs(uAirVelocityGeographical(1:nCells))), nCells
    write(*,*) "vAirVelocity: ", sum(abs(vAirVelocityGeographical(1:nCells))), nCells

    ! air stress on cells
    call plot_cells_spherical(mesh, airStressCellUGeographical, "airStressCellU.txt", itimestep)
    call plot_cells_spherical(mesh, airStressCellVGeographical, "airStressCellV.txt", itimestep)

    write(*,*) "airStressCellU: ", sum(abs(airStressCellUGeographical(1:nCells))), nCells
    write(*,*) "airStressCellV: ", sum(abs(airStressCellVGeographical(1:nCells))), nCells

    ! area on cell
    call plot_cells_spherical(mesh, iceAreaCell, "iceAreaCell.txt", itimestep)

    write(*,*) "iceAreaCell: ", sum(abs(iceAreaCell(1:nCells))), nCells

    ! air stress on vertices
    call plot_vertices_spherical(mesh, airStressVertexUGeographical, interiorVertex, "airStressVertexU.txt", itimestep)
    call plot_vertices_spherical(mesh, airStressVertexVGeographical, interiorVertex, "airStressVertexV.txt", itimestep)

    write(*,*) "airStressVertexU: ", sum(abs(airStressVertexUGeographical(1:nVertices))), nVertices, sum(interiorVertex)
    write(*,*) "airStressVertexV: ", sum(abs(airStressVertexVGeographical(1:nVertices))), nVertices, sum(interiorVertex)

    ! coriolis force coefficient
    call plot_vertices_spherical(mesh, totalMassVertexfVertex, interiorVertex, "totalMassVertexfVertex.txt", itimestep)

    write(*,*) "totalMassVertexfVertex: ", sum(abs(totalMassVertexfVertex(1:nVertices))), nVertices, sum(interiorVertex)

    ! ice internal stress
    call plot_vertices_spherical(mesh, stressDivergenceUGeographical, interiorVertex, "stressDivergenceU.txt", itimestep)
    call plot_vertices_spherical(mesh, stressDivergenceVGeographical, interiorVertex, "stressDivergenceV.txt", itimestep)

    write(*,*) "stressDivergenceU: ", sum(abs(stressDivergenceUGeographical(1:nVertices)), &
         mask=interiorVertex(1:nVertices)==1), nVertices, sum(interiorVertex)
    write(*,*) "stressDivergenceV: ", sum(abs(stressDivergenceVGeographical(1:nVertices)), &
         mask=interiorVertex(1:nVertices)==1), nVertices, sum(interiorVertex)

    ! ice velocity
    call plot_vertices_spherical(mesh, uVelocityGeographical, interiorVertex, "uVelocity.txt", itimestep)
    call plot_vertices_spherical(mesh, vVelocityGeographical, interiorVertex, "vVelocity.txt", itimestep)

    write(*,*) "uVelocity: ", sum(abs(uVelocityGeographical(1:nVertices)), &
         mask=interiorVertex(1:nVertices)==1), nVertices, sum(interiorVertex)
    write(*,*) "vVelocity: ", sum(abs(vVelocityGeographical(1:nVertices)), &
         mask=interiorVertex(1:nVertices)==1), nVertices, sum(interiorVertex)

    deallocate(uVelocityGeographical)
    deallocate(vVelocityGeographical)
    deallocate(airStressVertexUGeographical)
    deallocate(airStressVertexVGeographical)
    deallocate(airStressCellUGeographical)
    deallocate(airStressCellVGeographical)
    deallocate(uAirVelocityGeographical)
    deallocate(vAirVelocityGeographical)
    deallocate(stressDivergenceUGeographical)
    deallocate(stressDivergenceVGeographical)

  end subroutine test_case_output_global_velocity!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  test_case_output_square_advection
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine test_case_output_square_advection(block, itimestep)!{{{

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    integer, intent(in) :: &
         itimestep !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    call strip_plot(iceAreaCategory(1,1,:), "iceAreaCategory_strip.txt", itimestep)
    call strip_plot(surfaceTemperature(1,1,:), "surfaceTemperature_strip.txt", itimestep)
    call strip_plot(iceVolumeCategory(1,1,:), "iceVolumeCategory_strip.txt", itimestep)
    call strip_plot(snowVolumeCategory(1,1,:), "snowVolumeCategory_strip.txt", itimestep)

  end subroutine test_case_output_square_advection!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  test_case_output_advection
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine test_case_output_advection(block, itimestep)!{{{

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    integer, intent(in) :: &
         itimestep !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature

    if (mod(itimestep, 24) == 0) then

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
       call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

       call cice_plot_spherical_latlon(mesh, iceAreaCategory(1,1,:), "iceAreaCategory.txt", itimestep)
       call cice_plot_spherical_latlon(mesh, iceVolumeCategory(1,1,:), "iceVolumeCategory.txt", itimestep)
       call cice_plot_spherical_latlon(mesh, snowVolumeCategory(1,1,:), "snowVolumeCategory.txt", itimestep)
       call cice_plot_spherical_latlon(mesh, surfaceTemperature(1,1,:), "surfaceTemperature.txt", itimestep)

    endif

  end subroutine test_case_output_advection!}}}

!-----------------------------------------------------------------------
! Standard output
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_prepare_output_variables
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date July, 14th 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_prepare_output_variables(block)

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         plotting, &
         biogeochemistry

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCategoryPlot, &
         iceVolumeCategoryPlot, &
         snowVolumeCategoryPlot, &
         surfaceTemperaturePlot, &
         iceThicknessPlot, &
         snowThicknessPlot, &
         netBrineHeight, &
         totalVerticalSalinity, &
         totalSkeletalAlgae

    real(kind=RKIND), dimension(:,:), pointer :: &
         totalVerticalBiologyIce, &
         totalVerticalBiologySnow

    integer :: &
         iCell

    integer, pointer :: &
         nCells

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "plotting", plotting)
    call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    call MPAS_pool_get_array(plotting, "iceAreaCategoryPlot", iceAreaCategoryPlot)
    call MPAS_pool_get_array(plotting, "iceVolumeCategoryPlot", iceVolumeCategoryPlot)
    call MPAS_pool_get_array(plotting, "snowVolumeCategoryPlot", snowVolumeCategoryPlot)
    call MPAS_pool_get_array(plotting, "surfaceTemperaturePlot", surfaceTemperaturePlot)
    call MPAS_pool_get_array(plotting, "iceThicknessPlot", iceThicknessPlot)
    call MPAS_pool_get_array(plotting, "snowThicknessPlot", snowThicknessPlot)

    call MPAS_pool_get_array(biogeochemistry, "netBrineHeight", netBrineHeight)
    call MPAS_pool_get_array(biogeochemistry, "totalSkeletalAlgae", totalSkeletalAlgae)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalBiologyIce", totalVerticalBiologyIce)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalBiologySnow", totalVerticalBiologySnow)
    call MPAS_pool_get_array(biogeochemistry, "totalVerticalSalinity", totalVerticalSalinity)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    iceAreaCategoryPlot(:)    = iceAreaCategory(1,1,:)
    iceVolumeCategoryPlot(:)  = iceVolumeCategory(1,1,:)
    snowVolumeCategoryPlot(:) = snowVolumeCategory(1,1,:)
    surfaceTemperaturePlot(:) = surfaceTemperature(1,1,:)

    do iCell = 1, nCells

       if (iceAreaCategoryPlot(iCell) > 0.0_RKIND) then
          iceThicknessPlot(iCell)  = iceVolumeCategoryPlot(iCell)  / iceAreaCategoryPlot(iCell)
          snowThicknessPlot(iCell) = snowVolumeCategoryPlot(iCell) / iceAreaCategoryPlot(iCell)
       else
          netBrineHeight(iCell) = 0.0_RKIND
          totalVerticalSalinity(iCell) = 0.0_RKIND
          totalVerticalBiologyIce(:,iCell) = 0.0_RKIND
          totalVerticalBiologySnow(:,iCell) = 0.0_RKIND
          iceThicknessPlot(iCell)  = 0.0_RKIND
          snowThicknessPlot(iCell) = 0.0_RKIND
       endif
    enddo ! iCell

  end subroutine cice_prepare_output_variables

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_time_output_variable_real
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date August 14th 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_time_output_variable_real(mesh, array, clock, name)

    use mpas_timekeeping

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh

    real(kind=RKIND), dimension(:), intent(in) :: &
         array

    type(MPAS_Clock_type), intent(in) :: &
         clock

    character(len=*), intent(in) :: &
         name

    integer :: &
         iCell

    type(MPAS_Time_type) :: &
         curr_time

    character(len=strKIND) :: &
         dateTimeString, &
         filename

    curr_time = mpas_get_clock_time(clock, MPAS_NOW)

    call mpas_get_time(curr_time, dateTimeString=dateTimeString)

    write(filename,fmt='(a,a,a,a)') trim(name), "_", trim(dateTimeString), '.txt'

    open(22,file=trim(filename))

    do iCell = 1, size(array)

       write(22,*) iCell, array(iCell)

    enddo ! iCell

    close(22)

    open(22,file="diffs.txt",position="append")

    write(22,*) "diff ", trim(filename), " ../restart_1/", trim(filename)

    close(22)

  end subroutine cice_time_output_variable_real

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_time_output_variable_real_2D
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date August 10th 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_time_output_variable_real_2D(mesh, array, clock, name)

    use mpas_timekeeping

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         array

    type(MPAS_Clock_type), intent(in) :: &
         clock

    character(len=*), intent(in) :: &
         name

    integer :: &
         iCell1, iCell2

    type(MPAS_Time_type) :: &
         curr_time

    character(len=strKIND) :: &
         dateTimeString, &
         filename

    curr_time = mpas_get_clock_time(clock, MPAS_NOW)

    call mpas_get_time(curr_time, dateTimeString=dateTimeString)

    write(filename,fmt='(a,a,a,a)') trim(name), "_", trim(dateTimeString), '.txt'

    open(22,file=trim(filename))

    do iCell1 = 1, size(array,1)
    do iCell2 = 1, size(array,2)

       write(22,*) iCell1, iCell2, array(iCell1,iCell2)

    enddo ! iCell2
    enddo ! iCell1

    close(22)

    open(22,file="diffs.txt",position="append")

    write(22,*) "diff ", trim(filename), " ../restart_1/", trim(filename)

    close(22)

  end subroutine cice_time_output_variable_real_2D

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_time_output_variable_int
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date August 14th 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_time_output_variable_int(mesh, array, clock, name)

    use mpas_timekeeping

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh

    integer, dimension(:), intent(in) :: &
         array

    type(MPAS_Clock_type), intent(in) :: &
         clock

    character(len=*), intent(in) :: &
         name

    integer :: &
         iCell

    type(MPAS_Time_type) :: &
         curr_time

    character(len=strKIND) :: &
         dateTimeString, &
         filename

    curr_time = mpas_get_clock_time(clock, MPAS_NOW)

    call mpas_get_time(curr_time, dateTimeString=dateTimeString)

    write(filename,fmt='(a,a,a,a)') trim(name), "_", trim(dateTimeString), '.txt'

    open(22,file=trim(filename))

    do iCell = 1, size(array)

       write(22,*) iCell, array(iCell)

    enddo ! iCell

    close(22)

    open(22,file="diffs.txt",position="append")

    write(22,*) "diff ", trim(filename), " ../restart_1/", trim(filename)

    close(22)

  end subroutine cice_time_output_variable_int

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_time_output_variable_int_2D
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date August 19th 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_time_output_variable_int_2D(mesh, array, clock, name)

    use mpas_timekeeping

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh

    integer, dimension(:,:), intent(in) :: &
         array

    type(MPAS_Clock_type), intent(in) :: &
         clock

    character(len=*), intent(in) :: &
         name

    integer :: &
         iCell1, iCell2

    type(MPAS_Time_type) :: &
         curr_time

    character(len=strKIND) :: &
         dateTimeString, &
         filename

    curr_time = mpas_get_clock_time(clock, MPAS_NOW)

    call mpas_get_time(curr_time, dateTimeString=dateTimeString)

    write(filename,fmt='(a,a,a,a)') trim(name), "_", trim(dateTimeString), '.txt'

    open(22,file=trim(filename))

    do iCell1 = 1, size(array,1)
    do iCell2 = 1, size(array,2)

       write(22,*) iCell1, iCell2, array(iCell1,iCell2)

    enddo ! iCell1
    enddo ! iCell2

    close(22)

    open(22,file="diffs.txt",position="append")

    write(22,*) "diff ", trim(filename), " ../restart_1/", trim(filename)

    close(22)

    stop

  end subroutine cice_time_output_variable_int_2D

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_debug_output_array
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date Feburary 11th 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_debug_output_array(array,fileout)

    real(kind=RKIND), dimension(:), intent(in) :: &
         array

    character(len=*), optional :: fileout

    integer :: recl

    recl = size(array) * 8

    if (present(fileout)) then
       open(60,file=trim(fileout),form="unformatted",access="direct",recl=recl)
    else
       open(60,file="mpas_array_comparison.dat",form="unformatted",access="direct",recl=recl)
    endif

    write(60,rec=1) array

    close(60)

    !write(stderrUnit,*) "cice_debug_output_array"

  end subroutine cice_debug_output_array

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_debug_output_block
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date Feburary 11th 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_debug_output_block(block)

    type(block_type) :: block

    type(mpas_pool_type), pointer :: &
         pool

    real(kind=RKIND), dimension(:), pointer :: &
         array

    call MPAS_pool_get_subpool(block % structs, "", pool)

    call MPAS_pool_get_array(pool, "", array)

    call cice_debug_output_array(array)

    stop

  end subroutine cice_debug_output_block

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_debug_output_block
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date Feburary 11th 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_debug_output_domain(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(mpas_pool_type), pointer :: &
         pool

    !real(kind=RKIND), dimension(:,:,:), pointer :: &
    real(kind=RKIND), dimension(:), pointer :: &
         array

    !character(len=200), parameter :: &
    !     poolName  = "tracers", &
    !     fieldName = "iceVolumeCategory"
    character(len=200), parameter :: &
         poolName  = "atmos_coupling", &
         fieldName = "airPotentialTemperature"

    write(stderrUnit,*) trim(poolName), " ", trim(fieldName)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, trim(poolName), pool)

       call MPAS_pool_get_array(pool, trim(fieldName), array, 1)

       !call cice_debug_output_array(array(1,1,:))
       call cice_debug_output_array(array(:))

       block => block % next
    end do

    stop

  end subroutine cice_debug_output_domain

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
! cice_debug_write_out_field
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date April 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_debug_write_out_field(domain, poolname, fieldname, cicename)

    use cice_constants, only: &
         pii

    type(domain_type) :: domain

    character(len=*), intent(in) :: &
         poolname, &
         fieldname, &
         cicename

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         pool

    type(MPAS_pool_field_info_type) :: fieldInfo

    type(field1DReal), pointer :: array1D
    type(field2DReal), pointer :: array2D
    type(field3DReal), pointer :: array3D
    type(field4DReal), pointer :: array4D

    character(len=1000) :: fileout

    integer :: iDim1, iDim2, iDim3

    block => domain % blocklist

    call MPAS_pool_get_subpool(block % structs, trim(poolname), pool)

    call MPAS_pool_get_field_info(pool, trim(fieldname), fieldInfo)

    write(stderrUnit,*) "cice_debug_write_out_field: ", fieldInfo % nDims, " ", &
         trim(poolname), " ", trim(fieldname), " ", trim(cicename)

    select case (fieldInfo % nDims)
    case(1)

       call MPAS_pool_get_field(pool, trim(fieldname), array1D)

       if (trim(fieldname) == "lonCell") where(array1D % array > pii) array1D % array = array1D % array - 2.0_RKIND * pii

       if (array1D % isActive) then
          fileout = "comparison/pop_mpas_array_comparison_"//trim(cicename)//".dat"
          call cice_debug_output_array(array1D % array, trim(fileout))
       endif

    case(2)

       call MPAS_pool_get_field(pool, trim(fieldname), array2D)
       if (array2D % isActive) then
          do iDim1 = 1, size(array2D % array,1)
             write(fileout,fmt='(a,i1,a)') "comparison/pop_mpas_array_comparison_"//trim(cicename)//"_",iDim1,".dat"
             call cice_debug_output_array(array2D % array(iDim1,:), trim(fileout))
          enddo
       endif

    case(3)

       call MPAS_pool_get_field(pool, trim(fieldname), array3D)
       if (array3D % isActive) then
          do iDim1 = 1, size(array3D % array,1)
             do iDim2 = 1, size(array3D % array,2)
                if (size(array3D % array,1) == 1 .and. trim(fieldname) /= "snowEnthalpy" .and. &
                     trim(fieldname) /= "absorbedShortwaveSnowLayer") then
                   write(fileout,fmt='(a,i1,a)') "comparison/pop_mpas_array_comparison_"//trim(cicename)//"_",iDim2,".dat"
                else
                   write(fileout,fmt='(a,i1,a,i1,a)') &
                        "comparison/pop_mpas_array_comparison_"//trim(cicename)//"_",iDim1,"_",iDim2,".dat"
                endif
                call cice_debug_output_array(array3D % array(iDim1,iDim2,:), trim(fileout))
             enddo
          enddo
       endif

    case(4)

       call MPAS_pool_get_field(pool, trim(fieldname), array4D)
       if (array4D % isActive) then
          do iDim1 = 1, size(array4D % array,1)
             do iDim2 = 1, size(array4D % array,2)
                do iDim3 = 1, size(array4D % array,3)
                   write(fileout,fmt='(a,i1,a,i1,a,i1,a)') &
                        "comparison/pop_mpas_array_comparison_"//trim(cicename)//"_",iDim1,"_",iDim2,"_",iDim3,".dat"
                   call cice_debug_output_array(array4D % array(iDim1,iDim2,iDim3,:), trim(fileout))
                enddo
             enddo
          enddo
       endif

    end select

  end subroutine cice_debug_write_out_field

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
! cice_debug_write_out_fields
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date April 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_debug_write_out_fields(domain, callingSubroutine, location)

    type(domain_type) :: domain

    character(len=*), intent(in) :: callingSubroutine

    integer, intent(in) :: location
#if 0
    ! parameters where we want to writeout output
    character(len=strKIND), parameter :: callingSubroutineDesired = "cice_timestep"
    integer, parameter :: locationDesired = 2
    integer, parameter :: iStepDesired = 24

! cice_init_column_shortwave           1           0
! column_radiation           1           0
! column_radiation           2           0
! cice_init_column_shortwave           2           0
! column_prep_radiation           1           1
! column_prep_radiation           2           1
! column_vertical_thermodynamics           1           1
! column_vertical_thermodynamics           2           1
! column_itd_thermodynamics           1           1
! column_itd_thermodynamics           2           1
! column_ridging           1           1
! column_ridging           2           1
! column_radiation           1           1
! column_radiation           2           1
! cice_timestep
! cice_run_velocity_solver

    write(stderrUnit,*) trim(callingSubroutine), trim(callingSubroutineDesired), location, locationDesired, iStep, iStepDesired

    if (trim(callingSubroutine) == trim(callingSubroutineDesired) .and. &
        location == locationDesired .and. &
        iStep == iStepDesired) then

       write(stderrUnit,*) "cice_debug_write_out_fields: ", trim(callingSubroutine), location, iStep

       ! tracer
       call cice_debug_write_out_field(domain, "tracers", "iceAreaCategory", "aicen")
       call cice_debug_write_out_field(domain, "tracers", "iceVolumeCategory", "vicen")
       call cice_debug_write_out_field(domain, "tracers", "snowVolumeCategory", "vsnon")
       call cice_debug_write_out_field(domain, "tracers", "surfaceTemperature", "Tsfcn")
       call cice_debug_write_out_field(domain, "tracers", "iceEnthalpy", "qicen")
       call cice_debug_write_out_field(domain, "tracers", "iceSalinity", "sicen")
       call cice_debug_write_out_field(domain, "tracers", "snowEnthalpy", "qsnon")
       call cice_debug_write_out_field(domain, "tracers", "iceAge", "iagen")
       call cice_debug_write_out_field(domain, "tracers", "firstYearIceArea", "FYn")
       call cice_debug_write_out_field(domain, "tracers", "levelIceArea", "alvln")
       call cice_debug_write_out_field(domain, "tracers", "levelIceVolume", "vlvln")
       call cice_debug_write_out_field(domain, "tracers", "pondArea", "apndn")
       call cice_debug_write_out_field(domain, "tracers", "pondDepth", "hpndn")
       call cice_debug_write_out_field(domain, "tracers", "pondLidThickness", "ipndn")
       call cice_debug_write_out_field(domain, "tracers", "skeletalAlgaeConc", "bgc_N_skn")
       call cice_debug_write_out_field(domain, "tracers", "skeletalNitrateConc", "bgc_Nit_skn")
       call cice_debug_write_out_field(domain, "tracers", "verticalAlgaeConc", "bgc_Nn")
       call cice_debug_write_out_field(domain, "tracers", "verticalNitrateConc", "bgc_Nitn")
       call cice_debug_write_out_field(domain, "tracers", "brinefraction", "fbrin")
       call cice_debug_write_out_field(domain, "tracers", "verticalSalinity", "zSaln")

       call cice_debug_write_out_field(domain, "tracers_aggregate", "iceAreaCell", "aice")
       call cice_debug_write_out_field(domain, "tracers_aggregate", "iceVolumeCell", "vice")
       call cice_debug_write_out_field(domain, "tracers_aggregate", "snowVolumeCell", "vsno")
       call cice_debug_write_out_field(domain, "tracers_aggregate", "surfaceTemperatureCell", "Tsfc")
       call cice_debug_write_out_field(domain, "tracers_aggregate", "iceEnthalpyCell", "qice")
       call cice_debug_write_out_field(domain, "tracers_aggregate", "iceSalinityCell", "sice")
       call cice_debug_write_out_field(domain, "tracers_aggregate", "snowEnthalpyCell", "qsno")
       call cice_debug_write_out_field(domain, "tracers_aggregate", "iceAgeCell", "iage")
       call cice_debug_write_out_field(domain, "tracers_aggregate", "firstYearIceAreaCell", "FY")
       call cice_debug_write_out_field(domain, "tracers_aggregate", "levelIceAreaCell", "alvl")
       call cice_debug_write_out_field(domain, "tracers_aggregate", "levelIceVolumeCell", "vlvl")
       call cice_debug_write_out_field(domain, "tracers_aggregate", "pondAreaCell", "apnd")
       call cice_debug_write_out_field(domain, "tracers_aggregate", "pondDepthCell", "hpnd")
       call cice_debug_write_out_field(domain, "tracers_aggregate", "pondLidThicknessCell", "ipnd")

       call cice_debug_write_out_field(domain, "icestate", "iceAreaCellInitial", "aice_init")
       call cice_debug_write_out_field(domain, "icestate", "iceAreaCategoryInitial", "aicen_init")
       call cice_debug_write_out_field(domain, "icestate", "iceVolumeCategoryInitial", "vicen_init")
       !call cice_debug_write_out_field(domain, "icestate", "snowVolumeCategoryInitial", "vsnon_init")
       call cice_debug_write_out_field(domain, "icestate", "openWaterArea", "aice0")

       call cice_debug_write_out_field(domain, "atmos_coupling", "airPotentialTemperature", "potT")
       call cice_debug_write_out_field(domain, "atmos_coupling", "airTemperature", "Tair")
       call cice_debug_write_out_field(domain, "atmos_coupling", "airSpecificHumidity", "Qa")
       call cice_debug_write_out_field(domain, "atmos_coupling", "airDensity", "rhoa")
       call cice_debug_write_out_field(domain, "atmos_coupling", "shortwaveVisibleDirectDown", "swvdr")
       call cice_debug_write_out_field(domain, "atmos_coupling", "shortwaveVisibleDiffuseDown", "swvdf")
       call cice_debug_write_out_field(domain, "atmos_coupling", "shortwaveIRDirectDown", "swidr")
       call cice_debug_write_out_field(domain, "atmos_coupling", "shortwaveIRDiffuseDown", "swidf")
       call cice_debug_write_out_field(domain, "atmos_coupling", "longwaveDown", "flw")
       call cice_debug_write_out_field(domain, "atmos_coupling", "rainfallRate", "frain")
       call cice_debug_write_out_field(domain, "atmos_coupling", "snowfallRate", "fsnow")
       call cice_debug_write_out_field(domain, "atmos_coupling", "uAirVelocity", "uatm")
       call cice_debug_write_out_field(domain, "atmos_coupling", "vAirVelocity", "vatm")
       call cice_debug_write_out_field(domain, "atmos_coupling", "atmosReferenceTemperature2m", "")
       call cice_debug_write_out_field(domain, "atmos_coupling", "atmosReferenceHumidity2m", "")
       call cice_debug_write_out_field(domain, "atmos_coupling", "atmosReferenceSpeed10m", "")

       call cice_debug_write_out_field(domain, "atmos_forcing", "windSpeed", "wind")
       call cice_debug_write_out_field(domain, "atmos_forcing", "shortwaveDown", "fsw")
       !call cice_debug_write_out_field(domain, "atmos_forcing", "sensibleTransferCoefficient", "shcoef")
       !call cice_debug_write_out_field(domain, "atmos_forcing", "latentTransferCoefficient", "lhcoef")

       call cice_debug_write_out_field(domain, "ocean_coupling", "seaSurfaceTemperature", "sst")
       call cice_debug_write_out_field(domain, "ocean_coupling", "seaSurfaceSalinity", "sss")
       call cice_debug_write_out_field(domain, "ocean_coupling", "seaFreezingTemperature", "Tf")
       call cice_debug_write_out_field(domain, "ocean_coupling", "freezingMeltingPotential", "frzmlt")
       !call cice_debug_write_out_field(domain, "ocean_coupling", "uOceanVelocity", "uocn")
       !call cice_debug_write_out_field(domain, "ocean_coupling", "vOceanVelocity", "vocn")

       !call cice_debug_write_out_field(domain, "velocity_solver", "uOceanVelocityVertex", "uocn")
       !call cice_debug_write_out_field(domain, "velocity_solver", "vOceanVelocityVertex", "vocn")
       call cice_debug_write_out_field(domain, "velocity_geographical", "uOceanVelocityVertexGeo", "uocn")
       call cice_debug_write_out_field(domain, "velocity_geographical", "vOceanVelocityVertexGeo", "vocn")

       call cice_debug_write_out_field(domain, "ridging", "ridgeConvergence", "rdg_conv")
       call cice_debug_write_out_field(domain, "ridging", "ridgeShear", "rdg_shear")
       call cice_debug_write_out_field(domain, "ridging", "areaLossRidge", "dardg1dt")
       call cice_debug_write_out_field(domain, "ridging", "areaGainRidge", "dardg2dt")
       call cice_debug_write_out_field(domain, "ridging", "iceVolumeRidged", "dvirdgdt")
       call cice_debug_write_out_field(domain, "ridging", "openingRateRidge", "opening")
       call cice_debug_write_out_field(domain, "ridging", "ridgeParticipationFunction", "aparticn")
       call cice_debug_write_out_field(domain, "ridging", "ratioRidgeThicknessToIce", "krdgn")
       call cice_debug_write_out_field(domain, "ridging", "fractionNewRidgeArea", "aredistn")
       call cice_debug_write_out_field(domain, "ridging", "fractionNewRidgeVolume", "vredistn")
       call cice_debug_write_out_field(domain, "ridging", "areaLossRidgeCategory", "dardg1ndt")
       call cice_debug_write_out_field(domain, "ridging", "areaGainRidgeCategory", "dardg2ndt")
       call cice_debug_write_out_field(domain, "ridging", "iceVolumeRidgedCategory", "dvirdgndt")
       call cice_debug_write_out_field(domain, "ridging", "raftingIceArea", "araftn")
       call cice_debug_write_out_field(domain, "ridging", "raftingIceVolume", "vraftn")

       call cice_debug_write_out_field(domain, "melt_growth_rates", "lateralIceMeltFraction", "rside")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "surfaceIceMelt", "meltt")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "surfaceIceMeltCategory", "melttn")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "basalIceMelt", "meltb")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "basalIceMeltCategory", "meltbn")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "lateralIceMelt", "meltl")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "snowMelt", "melts")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "snowMeltCategory", "meltsn")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "congelation", "congel")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "congelationCategory", "congeln")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "snowiceFormation", "snoice")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "snowiceFormationCategory", "snoicen")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "snowThicknessChange", "dsnow")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "snowThicknessChangeCategory", "dsnown")
       call cice_debug_write_out_field(domain, "melt_growth_rates", "frazilFormation", "frazil")

       call cice_debug_write_out_field(domain, "atmos_fluxes", "surfaceHeatFlux", "fsurf")
       call cice_debug_write_out_field(domain, "atmos_fluxes", "surfaceHeatFluxCategory", "fsurfn")
       call cice_debug_write_out_field(domain, "atmos_fluxes", "surfaceConductiveFlux", "fcondtop")
       call cice_debug_write_out_field(domain, "atmos_fluxes", "surfaceConductiveFluxCategory", "fcondtopn")
       call cice_debug_write_out_field(domain, "atmos_fluxes", "longwaveUp", "flwout")
       call cice_debug_write_out_field(domain, "atmos_fluxes", "sensibleHeatFlux", "fsens")
       call cice_debug_write_out_field(domain, "atmos_fluxes", "latentHeatFlux", "flat")
       call cice_debug_write_out_field(domain, "atmos_fluxes", "latentHeatFluxCategory", "flatn")
       call cice_debug_write_out_field(domain, "atmos_fluxes", "evaporativeWaterFlux", "evap")

       call cice_debug_write_out_field(domain, "ocean_fluxes", "oceanFreshWaterFlux", "fresh")
       call cice_debug_write_out_field(domain, "ocean_fluxes", "oceanSaltFlux", "fsalt")
       call cice_debug_write_out_field(domain, "ocean_fluxes", "oceanHeatFlux", "fhocn")
       call cice_debug_write_out_field(domain, "ocean_fluxes", "oceanShortwaveFlux", "fswthru")

       call cice_debug_write_out_field(domain, "shortwave", "solarZenithAngleCosine", "coszen")
       call cice_debug_write_out_field(domain, "shortwave", "albedoVisibleDirectCategory", "alvdrn")
       call cice_debug_write_out_field(domain, "shortwave", "albedoVisibleDiffuseCategory", "alvdfn")
       call cice_debug_write_out_field(domain, "shortwave", "albedoIRDirectCategory", "alidrn")
       call cice_debug_write_out_field(domain, "shortwave", "albedoIRDiffuseCategory", "alidfn")
       call cice_debug_write_out_field(domain, "shortwave", "albedoVisibleDirectCell", "alvdr")
       call cice_debug_write_out_field(domain, "shortwave", "albedoVisibleDiffuseCell", "alvdf")
       call cice_debug_write_out_field(domain, "shortwave", "albedoIRDirectCell", "alidr")
       call cice_debug_write_out_field(domain, "shortwave", "albedoIRDiffuseCell", "alidf")
       call cice_debug_write_out_field(domain, "shortwave", "albedoVisibleDirectArea", "alvdr_ai")
       call cice_debug_write_out_field(domain, "shortwave", "albedoVisibleDiffuseArea", "alvdf_ai")
       call cice_debug_write_out_field(domain, "shortwave", "albedoIRDirectArea", "alidr_ai")
       call cice_debug_write_out_field(domain, "shortwave", "albedoIRDiffuseArea", "alidf_ai")
       call cice_debug_write_out_field(domain, "shortwave", "shortwaveScalingFactor", "scale_factor")
       call cice_debug_write_out_field(domain, "shortwave", "surfaceShortwaveFlux", "fswsfcn")
       call cice_debug_write_out_field(domain, "shortwave", "interiorShortwaveFlux", "fswintn")
       call cice_debug_write_out_field(domain, "shortwave", "penetratingShortwaveFlux", "fswthrun")
       call cice_debug_write_out_field(domain, "shortwave", "shortwaveLayerPenetration", "fswpenln")
       call cice_debug_write_out_field(domain, "shortwave", "absorbedShortwaveFlux", "fswabs")
       call cice_debug_write_out_field(domain, "shortwave", "absorbedShortwaveSnowLayer", "Sswabsn")
       call cice_debug_write_out_field(domain, "shortwave", "absorbedShortwaveIceLayer", "Iswabsn")
       call cice_debug_write_out_field(domain, "shortwave", "bareIceAlbedoCategory", "albicen")
       call cice_debug_write_out_field(domain, "shortwave", "snowAlbedoCategory", "albsnon")
       call cice_debug_write_out_field(domain, "shortwave", "pondAlbedoCategory", "albpndn")
       call cice_debug_write_out_field(domain, "shortwave", "bareIceAlbedoCell", "albice")
       call cice_debug_write_out_field(domain, "shortwave", "snowAlbedoCell", "albsno")
       call cice_debug_write_out_field(domain, "shortwave", "pondAlbedoCell", "albpnd")
       call cice_debug_write_out_field(domain, "shortwave", "effectivePondAreaCategory", "apeffn")
       call cice_debug_write_out_field(domain, "shortwave", "effectivePondAreaCell", "apeff_ai")

       call cice_debug_write_out_field(domain, "drag", "airOceanDragCoefficientRatio", "Cdn_atm_ocn")
       call cice_debug_write_out_field(domain, "drag", "oceanDragCoefficient", "Cdn_ocn")
       call cice_debug_write_out_field(domain, "drag", "oceanDragCoefficientSkin", "Cdn_ocn_skin")
       call cice_debug_write_out_field(domain, "drag", "oceanDragCoefficientFloe", "Cdn_ocn_floe")
       call cice_debug_write_out_field(domain, "drag", "oceanDragCoefficientKeel", "Cdn_ocn_keel")
       call cice_debug_write_out_field(domain, "drag", "airDragCoefficient", "Cdn_atm")
       call cice_debug_write_out_field(domain, "drag", "airDragCoefficientSkin", "Cdn_atm_skin")
       call cice_debug_write_out_field(domain, "drag", "airDragCoefficientFloe", "Cdn_atm_floe")
       call cice_debug_write_out_field(domain, "drag", "airDragCoefficientPond", "Cdn_atm_pond")
       call cice_debug_write_out_field(domain, "drag", "airDragCoefficientRidge", "Cdn_atm_rdg")
       call cice_debug_write_out_field(domain, "drag", "dragFreeboard", "hfreebd")
       call cice_debug_write_out_field(domain, "drag", "dragIceSnowDraft", "hdraft")
       call cice_debug_write_out_field(domain, "drag", "dragRidgeHeight", "hridge")
       call cice_debug_write_out_field(domain, "drag", "dragRidgeSeparation", "distrdg")
       call cice_debug_write_out_field(domain, "drag", "dragKeelDepth", "hkeel")
       call cice_debug_write_out_field(domain, "drag", "dragKeelSeparation", "dkeel")
       call cice_debug_write_out_field(domain, "drag", "dragFloeLength", "lfloe")
       call cice_debug_write_out_field(domain, "drag", "dragFloeSeparation", "dfloe")

       call cice_debug_write_out_field(domain, "ponds", "pondFreshWaterFlux", "fpond")
       call cice_debug_write_out_field(domain, "ponds", "pondSnowDepthDifference", "dhsn")
       call cice_debug_write_out_field(domain, "ponds", "pondLidMeltFluxFraction", "ffracn")

       call cice_debug_write_out_field(domain, "aerosols", "atmosAerosolFlux", "faero_atm")
       call cice_debug_write_out_field(domain, "aerosols", "oceanAerosolFlux", "faero_ocn")
       !call cice_debug_write_out_field(domain, "aerosols", "aerosolMassExtinctionCrossSection", "kaer_tab")
       !call cice_debug_write_out_field(domain, "aerosols", "aerosolSingleScatterAlbedo", "waer_tab")
       !call cice_debug_write_out_field(domain, "aerosols", "aerosolAsymmetryParameter", "gaer_tab")

       !call cice_debug_write_out_field(domain, "biogeochemistry", "newlyFormedIce", "first_ice")
       call cice_debug_write_out_field(domain, "biogeochemistry", "oceanBioFluxes", "flux_bio")
       call cice_debug_write_out_field(domain, "biogeochemistry", "oceanBioConcentrations", "ocean_bio")

       call cice_debug_write_out_field(domain, "initial", "initialSalinityProfile", "salinz")
       call cice_debug_write_out_field(domain, "initial", "initialMeltingTemperatureProfile", "Tmltz")

       call cice_debug_write_out_field(domain, "mesh", "lonCell", "TLON")
       call cice_debug_write_out_field(domain, "mesh", "latCell", "TLAT")

       ! velocity related stuff
       call cice_debug_write_out_field(domain, "mesh", "lonVertex", "ULON")
       call cice_debug_write_out_field(domain, "mesh", "latVertex", "ULAT")

       !call cice_debug_write_out_field(domain, "velocity_solver", "uVelocity", "uvel")
       !call cice_debug_write_out_field(domain, "velocity_solver", "vVelocity", "vvel")
       !call cice_debug_write_out_field(domain, "velocity_solver", "stressDivergenceU", "strintx")
       !call cice_debug_write_out_field(domain, "velocity_solver", "stressDivergenceV", "strinty")
       !call cice_debug_write_out_field(domain, "velocity_solver", "airStressVertexU", "strairx")
       !call cice_debug_write_out_field(domain, "velocity_solver", "airStressVertexV", "strairy")
       !call cice_debug_write_out_field(domain, "velocity_solver", "oceanStressU", "strocnx")
       !call cice_debug_write_out_field(domain, "velocity_solver", "oceanStressV", "strocny")
       !call cice_debug_write_out_field(domain, "velocity_solver", "surfaceTiltForceU", "strtltx")
       !call cice_debug_write_out_field(domain, "velocity_solver", "surfaceTiltForceV", "strtlty")

       call cice_debug_write_out_field(domain, "velocity_solver", "oceanStressCellU", "strocnxT")
       call cice_debug_write_out_field(domain, "velocity_solver", "oceanStressCellV", "strocnyT")

       call cice_debug_write_out_field(domain, "velocity_solver", "airStressCellU", "strairxT")
       call cice_debug_write_out_field(domain, "velocity_solver", "airStressCellV", "strairyT")

       call cice_debug_write_out_field(domain, "velocity_solver", "icePressure", "strength")

       call cice_debug_write_out_field(domain, "velocity_geographical", "uVelocityGeo", "uvel")
       call cice_debug_write_out_field(domain, "velocity_geographical", "vVelocityGeo", "vvel")
       call cice_debug_write_out_field(domain, "velocity_geographical", "stressDivergenceUGeo", "strintx")
       call cice_debug_write_out_field(domain, "velocity_geographical", "stressDivergenceVGeo", "strinty")
       call cice_debug_write_out_field(domain, "velocity_geographical", "airStressVertexUGeo", "strairx")
       call cice_debug_write_out_field(domain, "velocity_geographical", "airStressVertexVGeo", "strairy")
       call cice_debug_write_out_field(domain, "velocity_geographical", "oceanStressUGeo", "strocnx")
       call cice_debug_write_out_field(domain, "velocity_geographical", "oceanStressVGeo", "strocny")
       call cice_debug_write_out_field(domain, "velocity_geographical", "surfaceTiltForceUGeo", "strtltx")
       call cice_debug_write_out_field(domain, "velocity_geographical", "surfaceTiltForceVGeo", "strtlty")

       !call cice_debug_write_out_field(domain, "mesh_test", "vectorCellU", "vectorCellU")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorCellV", "vectorCellV")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexU", "vectorVertexU")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexV", "vectorVertexV")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorCellURotate", "vectorCellURotate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorCellVRotate", "vectorCellVRotate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexURotate", "vectorVertexURotate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexVRotate", "vectorVertexVRotate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexUInterpolate", "vectorVertexUInterpolate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexVInterpolate", "vectorVertexVInterpolate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexURotateInterpolate", "vectorVertexURotateInterpolate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexVRotateInterpolate", "vectorVertexVRotateInterpolate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorCellURotateRotate", "vectorCellURotateRotate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorCellVRotateRotate", "vectorCellVRotateRotate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexURotateRotate", "vectorVertexURotateRotate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexVRotateRotate", "vectorVertexVRotateRotate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexURotateInterpolateRotate", &
       !                                                     "vectorVertexURotateInterpolateRotate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexVRotateInterpolateRotate", &
       !                                                     "vectorVertexVRotateInterpolateRotate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorCellUDiff", "vectorCellUDiff")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorCellVDiff", "vectorCellVDiff")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexUDiff", "vectorVertexUDiff")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexVDiff", "vectorVertexVDiff")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexUInterpolateDiff", "vectorVertexUInterpolateDiff")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexVInterpolateDiff", "vectorVertexVInterpolateDiff")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexUDiffRotate", "vectorVertexUDiffRotate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexVDiffRotate", "vectorVertexVDiffRotate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexUDiffNoRotate", "vectorVertexUDiffNoRotate")
       !call cice_debug_write_out_field(domain, "mesh_test", "vectorVertexVDiffNoRotate", "vectorVertexVDiffNoRotate")


       stop

    endif
#endif
  end subroutine cice_debug_write_out_fields

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_writeout_config
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date April 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_writeout_config(domain, unit, mpasConfig, popConfig, type)

    type(domain_type) :: domain

    integer, intent(in) :: unit

    character(len=*), intent(in) :: mpasConfig, popConfig, type

    character(len=20) :: popConfigFixed

    logical, pointer :: configLogical
    integer, pointer :: configInteger
    character(len=strKIND), pointer :: configCharacter
    real(kind=RKIND), pointer :: configReal
    integer :: configIntegerInt

    write(popConfigFixed,*) trim(popConfig)

    select case (trim(type))
    case ("logical")
       call MPAS_pool_get_config(domain % blocklist % configs, trim(mpasConfig), configLogical)
       write(unit,*) popConfigFixed, configLogical
    case ("integer")
       call MPAS_pool_get_config(domain % blocklist % configs, trim(mpasConfig), configInteger)
       write(unit,*) popConfigFixed, configInteger
    case ("character")
       call MPAS_pool_get_config(domain % blocklist % configs, trim(mpasConfig), configCharacter)
       if ( trim(mpasConfig) == "config_thermodynamics_type" .or. &
            trim(mpasConfig) == "config_itd_conversion_type" .or. &
            trim(mpasConfig) == "config_category_bounds_type" .or. &
            trim(mpasConfig) == "config_ice_strength_formulation" .or. &
            trim(mpasConfig) == "config_ridging_participation_function" .or. &
            trim(mpasConfig) == "config_ridging_redistribution_function") then
          configIntegerInt = config_cice_int(trim(mpasConfig), configCharacter)
          write(unit,*) popConfigFixed, configIntegerInt
       else
          write(unit,*) popConfigFixed, " ", trim(configCharacter)
       endif
    case ("real")
       call MPAS_pool_get_config(domain % blocklist % configs, trim(mpasConfig), configReal)
       write(unit,*) popConfigFixed, configReal
    end select

  end subroutine cice_writeout_config

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  config_cice_int
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date April 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  function config_cice_int(configName, configValue) result(configValueCice)

    character(len=*), intent(in) :: &
         configName, &
         configValue

    integer :: configValueCice

    select case (trim(configName))

    ! ktherm
    case ("config_thermodynamics_type")

       select case (trim(configValue))
       case ("zero layer")
          configValueCice = 0
       case ("BL99")
          configValueCice = 1
       case ("mushy")
          configValueCice = 2
       end select

    ! kitd
    case ("config_itd_conversion_type")

       select case (trim(configValue))
       case ("delta function")
          configValueCice = 0
       case ("linear remap")
          configValueCice = 1
       end select

    ! kcatbound
    case ("config_category_bounds_type")

       select case (trim(configValue))
       case ("single category")
          configValueCice = -1
       case ("original")
          configValueCice = 0
       case ("new")
          configValueCice = 1
       case ("WMO")
          configValueCice = 2
       case ("asymptotic")
          configValueCice = 3
       end select

    ! kstrength
    case ("config_ice_strength_formulation")

       select case (trim(configValue))
       case ("Hibler79")
          configValueCice = 0
       case ("Rothrock75")
          configValueCice = 1
       end select

    ! krdg_partic
    case ("config_ridging_participation_function")

       select case (trim(configValue))
       case ("Thorndike75")
          configValueCice = 0
       case ("exponential")
          configValueCice = 1
       end select

    ! krdg_redist
    case ("config_ridging_redistribution_function")

       select case (trim(configValue))
       case ("Hibler80")
          configValueCice = 0
       case ("exponential")
          configValueCice = 1
       end select

    end select

  end function config_cice_int

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_writeout_configs
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date April 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_writeout_configs(domain)

    type(domain_type) :: domain

    open(55,file="pop_mpas_config_comparison_mpas.txt",action="write")

    call cice_writeout_config(domain, 55, "config_use_ice_age", "tr_iage", "logical")
    call cice_writeout_config(domain, 55, "config_use_first_year_ice", "tr_FY", "logical")
    call cice_writeout_config(domain, 55, "config_use_level_ice", "tr_lvl", "logical")
    call cice_writeout_config(domain, 55, "config_use_meltponds", "tr_pond", "logical")
    call cice_writeout_config(domain, 55, "config_use_cesm_meltponds", "tr_pond_cesm", "logical")
    call cice_writeout_config(domain, 55, "config_use_level_meltponds", "tr_pond_lvl", "logical")
    call cice_writeout_config(domain, 55, "config_use_topo_meltponds", "tr_pond_topo", "logical")
    call cice_writeout_config(domain, 55, "config_use_aerosols", "tr_aero", "logical")
    call cice_writeout_config(domain, 55, "config_use_brine", "tr_brine", "logical")
    call cice_writeout_config(domain, 55, "config_use_vertical_zsalinity", "solve_zsal", "logical")
    call cice_writeout_config(domain, 55, "config_use_vertical_biochemistry", "solve_zbgc", "logical")
    call cice_writeout_config(domain, 55, "config_use_shortwave_bioabsorption", "dEdd_algae", "logical")
    call cice_writeout_config(domain, 55, "config_use_vertical_tracers", "z_tracers", "logical")
    call cice_writeout_config(domain, 55, "config_use_skeletal_biochemistry", "skl_bgc", "logical")
    call cice_writeout_config(domain, 55, "config_use_nitrate", "tr_bgc_Nit", "logical")
    call cice_writeout_config(domain, 55, "config_use_carbon", "tr_bgc_C", "logical")
    call cice_writeout_config(domain, 55, "config_use_chlorophyll", "tr_bgc_chl", "logical")
    call cice_writeout_config(domain, 55, "config_use_ammonium", "tr_bgc_Am", "logical")
    call cice_writeout_config(domain, 55, "config_use_silicate", "tr_bgc_Sil", "logical")
    call cice_writeout_config(domain, 55, "config_use_DMS", "tr_bgc_DMS", "logical")
    call cice_writeout_config(domain, 55, "config_use_nonreactive", "tr_bgc_PON", "logical")
    call cice_writeout_config(domain, 55, "config_use_humics", "tr_bgc_hum", "logical")
    call cice_writeout_config(domain, 55, "config_use_DON", "tr_bgc_DON", "logical")
    call cice_writeout_config(domain, 55, "config_use_iron", "tr_bgc_iron", "logical")
    call cice_writeout_config(domain, 55, "config_use_modal_aerosols", "modal_aero", "logical")
    call cice_writeout_config(domain, 55, "config_use_zaerosols", "tr_zaero", "logical")

    call cice_writeout_config(domain, 55, "config_shortwave_type", "shortwave", "character")
    call cice_writeout_config(domain, 55, "config_albedo_type", "albedo_type", "character")
    call cice_writeout_config(domain, 55, "config_visible_ice_albedo", "albicev", "real")
    call cice_writeout_config(domain, 55, "config_infrared_ice_albedo", "albicei", "real")
    call cice_writeout_config(domain, 55, "config_visible_snow_albedo", "albsnowv", "real")
    call cice_writeout_config(domain, 55, "config_infrared_snow_albedo", "albsnowi", "real")
    call cice_writeout_config(domain, 55, "config_variable_albedo_thickness_limit", "ahmax", "real")
    call cice_writeout_config(domain, 55, "config_ice_shortwave_tuning_parameter", "R_ice", "real")
    call cice_writeout_config(domain, 55, "config_pond_shortwave_tuning_parameter", "R_pnd", "real")
    call cice_writeout_config(domain, 55, "config_snow_shortwave_tuning_parameter", "R_snw", "real")
    call cice_writeout_config(domain, 55, "config_temp_change_snow_grain_radius_change", "dT_mlt", "real")
    call cice_writeout_config(domain, 55, "config_max_melting_snow_grain_radius", "rsnw_mlt", "real")

    call cice_writeout_config(domain, 55, "config_snow_to_ice_transition_depth", "hs0", "real")
    call cice_writeout_config(domain, 55, "config_adjust_shortwave_meltwater_infiltration", "snowinfil", "logical")
    call cice_writeout_config(domain, 55, "config_pond_refreezing_type", "frzpnd", "character")
    call cice_writeout_config(domain, 55, "config_pond_flushing_timescale", "dpscale", "real")
    call cice_writeout_config(domain, 55, "config_min_meltwater_retained_fraction", "rfracmin", "real")
    call cice_writeout_config(domain, 55, "config_max_meltwater_retained_fraction", "rfracmax", "real")
    call cice_writeout_config(domain, 55, "config_pond_depth_to_fraction_ratio", "pndaspect", "real")
    call cice_writeout_config(domain, 55, "config_snow_on_pond_ice_tapering_parameter", "hs1", "real")
    call cice_writeout_config(domain, 55, "config_critical_pond_ice_thickness", "hp1", "real")

    call cice_writeout_config(domain, 55, "config_thermodynamics_type", "ktherm", "character")
    call cice_writeout_config(domain, 55, "config_heat_conductivity_type", "conduct", "character")
    call cice_writeout_config(domain, 55, "config_rapid_mode_channel_radius", "a_rapid_mode", "real")
    call cice_writeout_config(domain, 55, "config_rapid_model_critical_Ra", "Rac_rapid_mode", "real")
    call cice_writeout_config(domain, 55, "config_rapid_mode_aspect_ratio", "aspect_rapid_mode", "real")
    call cice_writeout_config(domain, 55, "config_slow_mode_drainage_strength", "dSdt_slow_mode", "real")
    call cice_writeout_config(domain, 55, "config_slow_mode_critical_porosity", "phi_c_slow_mode", "real")
    call cice_writeout_config(domain, 55, "config_congelation_ice_porosity", "phi_i_mushy", "real")

    call cice_writeout_config(domain, 55, "config_itd_conversion_type", "kitd", "character")
    call cice_writeout_config(domain, 55, "config_category_bounds_type", "kcatbound", "character")

    call cice_writeout_config(domain, 55, "config_ice_strength_formulation", "kstrength", "character")
    call cice_writeout_config(domain, 55, "config_ridging_participation_function", "krdg_partic", "character")
    call cice_writeout_config(domain, 55, "config_ridging_redistribution_function", "krdg_redist", "character")
    call cice_writeout_config(domain, 55, "config_ridiging_efolding_scale", "mu_rdg", "real")

    call cice_writeout_config(domain, 55, "config_atmos_boundary_method", "atmbndy", "character")
    call cice_writeout_config(domain, 55, "config_calc_surface_stresses", "calc_strair", "logical")
    call cice_writeout_config(domain, 55, "config_calc_surface_temperature", "calc_Tsfc", "logical")
    call cice_writeout_config(domain, 55, "config_use_form_drag", "formdrag", "logical")
    call cice_writeout_config(domain, 55, "config_use_high_frequency_coupling", "highfreq", "logical")
    call cice_writeout_config(domain, 55, "config_boundary_layer_iteration_number", "natmiter", "integer")

    call cice_writeout_config(domain, 55, "config_use_ocean_mixed_layer", "oceanmixed_ice", "logical")
    call cice_writeout_config(domain, 55, "config_min_friction_velocity", "ustar_min", "real")

    close(55)

  end subroutine cice_writeout_configs

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_debug_rotate_winds
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date April 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_debug_rotate_winds(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         atmos_coupling, &
         ocean_coupling

    real(kind=RKIND), dimension(:), pointer :: &
         uAirVelocity, &
         vAirVelocity, &
         uOceanVelocity, &
         vOceanVelocity

    real(kind=RKIND), dimension(:,:), allocatable :: &
         ULAT, &
         ULON, &
         HTN, &
         HTE, &
         HUS, &
         HUW, &
         ANGLE, &
         ANGLET

    integer, parameter :: &
         nx = 320, &
         ny = 384

    integer, dimension(:), allocatable :: &
         iCell_map, i_map, j_map

    integer, pointer :: &
         nCells

    integer :: &
         recl, &
         iCell

    real(kind=RKIND) :: &
         angle_0, &
         angle_w, &
         angle_s, &
         angle_sw, &
         workx, &
         worky

    real(kind=RKIND), parameter :: &
         pi  = 3.14159265358979323846_RKIND, &
         pi2 = 2.0_RKIND * pi

    block => domain % blocklist
    do while (associated(block))

    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

    call MPAS_pool_get_array(atmos_coupling, "uAirVelocity", uAirVelocity)
    call MPAS_pool_get_array(atmos_coupling, "vAirVelocity", vAirVelocity)

    call MPAS_pool_get_array(ocean_coupling, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(ocean_coupling, "vOceanVelocity", vOceanVelocity)

    call MPAS_pool_get_dimension(atmos_coupling, "nCells", nCells)

    allocate(ULAT(0:nx+1,0:ny+1))
    allocate(ULON(0:nx+1,0:ny+1))
    allocate(HTN(0:nx+1,0:ny+1))
    allocate(HTE(0:nx+1,0:ny+1))
    allocate(HUS(0:nx+1,0:ny+1))
    allocate(HUW(0:nx+1,0:ny+1))
    allocate(ANGLE(0:nx+1,0:ny+1))
    allocate(ANGLET(0:nx+1,0:ny+1))

    ULAT = 0.0_RKIND
    ULON = 0.0_RKIND
    HTN = 0.0_RKIND
    HTE = 0.0_RKIND
    HUS = 0.0_RKIND
    HUW = 0.0_RKIND
    ANGLE = 0.0_RKIND
    ANGLET = 0.0_RKIND

    recl = 2 * nx * ny
    open(55, file="/usr/projects/climate/akt/cice_mpas_comparison/cice/column_pkg/rundir/grid", &
         form="unformatted", access="direct", recl=recl)
    read(55,rec=1) ULAT(1:nx,1:ny)
    read(55,rec=2) ULON(1:nx,1:ny)
    read(55,rec=3) HTN(1:nx,1:ny)
    read(55,rec=4) HTE(1:nx,1:ny)
    read(55,rec=5) HUS(1:nx,1:ny)
    read(55,rec=6) HUW(1:nx,1:ny)
    read(55,rec=7) ANGLE(1:nx,1:ny)
    close(55)

    ! halo
    ANGLE(0,:)    = ANGLE(nx,:)
    ANGLE(nx+1,:) = ANGLE(1,:)
    ANGLE(:,0)    = ANGLE(:,ny)
    ANGLE(:,ny+1) = ANGLE(:,1)

    allocate(iCell_map(nCells))
    allocate(i_map(nCells))
    allocate(j_map(nCells))

    open(55,file="/usr/projects/climate/akt/cice_mpas_comparison/cice/column_pkg/rundir/mapping_pop_to_mpas_gx1.txt")
    do iCell = 1, nCells
       read(55,*) iCell_map(iCell), i_map(iCell), j_map(iCell)
    enddo ! iCell
    close(55)

    do iCell = 1, nCells

       angle_0  = ANGLE(i_map(iCell)  ,j_map(iCell)  ) !   w----0
       angle_w  = ANGLE(i_map(iCell)-1,j_map(iCell)  ) !   |    |
       angle_s  = ANGLE(i_map(iCell),  j_map(iCell)-1) !   |    |
       angle_sw = ANGLE(i_map(iCell)-1,j_map(iCell)-1) !   sw---s

       if ( angle_0 < 0.0_RKIND ) then
          if ( abs(angle_w - angle_0) > pi) &
                   angle_w = angle_w  - pi2
          if ( abs(angle_s - angle_0) > pi) &
                   angle_s = angle_s  - pi2
          if ( abs(angle_sw - angle_0) > pi) &
                   angle_sw = angle_sw - pi2
       endif

       ANGLET(i_map(iCell),j_map(iCell)) = angle_0 * 0.25_RKIND + angle_w * 0.25_RKIND &
                                         + angle_s * 0.25_RKIND + angle_sw* 0.25_RKIND

    enddo ! iCell

    do iCell = 1, nCells

       workx      = uAirVelocity(iCell) ! wind velocity, m/s
       worky      = vAirVelocity(iCell)
       uAirVelocity(iCell) = workx*cos(ANGLET(i_map(iCell),j_map(iCell))) & ! convert to POP grid
                           + worky*sin(ANGLET(i_map(iCell),j_map(iCell)))   ! note uatm, vatm, wind
       vAirVelocity(iCell) = worky*cos(ANGLET(i_map(iCell),j_map(iCell))) & !  are on the T-grid here
                           - workx*sin(ANGLET(i_map(iCell),j_map(iCell)))

       !workx      = uOceanVelocity(iCell)
       !worky      = vOceanVelocity(iCell)
       !uOceanVelocity(iCell) = workx*cos(ANGLET(i_map(iCell),j_map(iCell))) & ! convert to POP grid
       !                      + worky*sin(ANGLET(i_map(iCell),j_map(iCell)))   ! note uatm, vatm, wind
       !vOceanVelocity(iCell) = worky*cos(ANGLET(i_map(iCell),j_map(iCell))) & !  are on the T-grid here
       !                      - workx*sin(ANGLET(i_map(iCell),j_map(iCell)))

    enddo ! iCell

    deallocate(ULAT)
    deallocate(ULON)
    deallocate(HTN)
    deallocate(HTE)
    deallocate(HUS)
    deallocate(HUW)
    deallocate(ANGLE)
    deallocate(ANGLET)

    deallocate(iCell_map)
    deallocate(i_map)
    deallocate(j_map)

    block => block % next
    enddo

  end subroutine cice_debug_rotate_winds

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_debug_processor_number
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date April 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_debug_processor_number(domain)

    type(domain_type) :: domain

    real(kind=RKIND) :: proc_id

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: diagnosticsPool

    real(kind=RKIND), dimension(:), pointer :: processorNumber

    proc_id = real(domain % dminfo % my_proc_id, RKIND)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnosticsPool)

       call MPAS_pool_get_array(diagnosticsPool, "processorNumber", processorNumber)

       processorNumber = proc_id

       block => block % next
    enddo

  end subroutine cice_debug_processor_number

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_debug_forcing_restart_test
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date April 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_debug_forcing_restart_test(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         atmosCouplingPool, &
         atmosCouplingRestartTestPool, &
         oceanCouplingPool, &
         oceanCouplingRestartTestPool

    real(kind=RKIND), dimension(:), pointer :: &
         airLevelHeight, &
         airPotentialTemperature, &
         airTemperature, &
         airSpecificHumidity, &
         airDensity, &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown, &
         longwaveDown, &
         rainfallRate, &
         snowfallRate, &
         uAirVelocity, &
         vAirVelocity, &
         atmosReferenceSpeed10m, &
         atmosReferenceTemperature2m, &
         atmosReferenceHumidity2m, &
         airLevelHeightRestartTest, &
         airPotentialTemperatureRestartTest, &
         airTemperatureRestartTest, &
         airSpecificHumidityRestartTest, &
         airDensityRestartTest, &
         shortwaveVisibleDirectDownRestartTest, &
         shortwaveVisibleDiffuseDownRestartTest, &
         shortwaveIRDirectDownRestartTest, &
         shortwaveIRDiffuseDownRestartTest, &
         longwaveDownRestartTest, &
         rainfallRateRestartTest, &
         snowfallRateRestartTest, &
         uAirVelocityRestartTest, &
         vAirVelocityRestartTest, &
         atmosReferenceSpeed10mRestartTest, &
         atmosReferenceTemperature2mRestartTest, &
         atmosReferenceHumidity2mRestartTest

    real(kind=RKIND), dimension(:), pointer :: &
         seaSurfaceTemperature, &
         seaSurfaceSalinity, &
         seaFreezingTemperature, &
         freezingMeltingPotential, &
         uOceanVelocity, &
         vOceanVelocity, &
         seaSurfaceTiltU, &
         seaSurfaceTiltV, &
         oceanMixedLayerDepth, &
         oceanHeatFluxConvergence, &
         seaSurfaceTemperatureRestartTest, &
         seaSurfaceSalinityRestartTest, &
         seaFreezingTemperatureRestartTest, &
         freezingMeltingPotentialRestartTest, &
         uOceanVelocityRestartTest, &
         vOceanVelocityRestartTest, &
         seaSurfaceTiltURestartTest, &
         seaSurfaceTiltVRestartTest, &
         oceanMixedLayerDepthRestartTest, &
         oceanHeatFluxConvergenceRestartTest

    block => domain % blocklist
    do while (associated(block))

       ! atmosphere coupling
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCouplingPool)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling_restart_test", atmosCouplingRestartTestPool)

       call MPAS_pool_get_array(atmosCouplingPool, "airLevelHeight", airLevelHeight)
       call MPAS_pool_get_array(atmosCouplingPool, "airPotentialTemperature", airPotentialTemperature)
       call MPAS_pool_get_array(atmosCouplingPool, "airTemperature", airTemperature)
       call MPAS_pool_get_array(atmosCouplingPool, "airSpecificHumidity", airSpecificHumidity)
       call MPAS_pool_get_array(atmosCouplingPool, "airDensity", airDensity)
       call MPAS_pool_get_array(atmosCouplingPool, "shortwaveVisibleDirectDown", shortwaveVisibleDirectDown)
       call MPAS_pool_get_array(atmosCouplingPool, "shortwaveVisibleDiffuseDown", shortwaveVisibleDiffuseDown)
       call MPAS_pool_get_array(atmosCouplingPool, "shortwaveIRDirectDown", shortwaveIRDirectDown)
       call MPAS_pool_get_array(atmosCouplingPool, "shortwaveIRDiffuseDown", shortwaveIRDiffuseDown)
       call MPAS_pool_get_array(atmosCouplingPool, "longwaveDown", longwaveDown)
       call MPAS_pool_get_array(atmosCouplingPool, "rainfallRate", rainfallRate)
       call MPAS_pool_get_array(atmosCouplingPool, "snowfallRate", snowfallRate)
       call MPAS_pool_get_array(atmosCouplingPool, "uAirVelocity", uAirVelocity)
       call MPAS_pool_get_array(atmosCouplingPool, "vAirVelocity", vAirVelocity)
       call MPAS_pool_get_array(atmosCouplingPool, "atmosReferenceSpeed10m", atmosReferenceSpeed10m)
       call MPAS_pool_get_array(atmosCouplingPool, "atmosReferenceTemperature2m", atmosReferenceTemperature2m)
       call MPAS_pool_get_array(atmosCouplingPool, "atmosReferenceHumidity2m", atmosReferenceHumidity2m)

       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "airLevelHeightRestartTest", airLevelHeightRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "airPotentialTemperatureRestartTest", &
                                                               airPotentialTemperatureRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "airTemperatureRestartTest", airTemperatureRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "airSpecificHumidityRestartTest", airSpecificHumidityRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "airDensityRestartTest", airDensityRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "shortwaveVisibleDirectDownRestartTest", &
                                                               shortwaveVisibleDirectDownRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "shortwaveVisibleDiffuseDownRestartTest", &
                                                               shortwaveVisibleDiffuseDownRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "shortwaveIRDirectDownRestartTest", shortwaveIRDirectDownRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "shortwaveIRDiffuseDownRestartTest", &
                                                               shortwaveIRDiffuseDownRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "longwaveDownRestartTest", longwaveDownRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "rainfallRateRestartTest", rainfallRateRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "snowfallRateRestartTest", snowfallRateRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "uAirVelocityRestartTest", uAirVelocityRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "vAirVelocityRestartTest", vAirVelocityRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "atmosReferenceSpeed10mRestartTest", &
                                                               atmosReferenceSpeed10mRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "atmosReferenceTemperature2mRestartTest", &
                                                               atmosReferenceTemperature2mRestartTest)
       call MPAS_pool_get_array(atmosCouplingRestartTestPool, "atmosReferenceHumidity2mRestartTest", &
                                                               atmosReferenceHumidity2mRestartTest)

       airLevelHeightRestartTest              = airLevelHeightRestartTest              + airLevelHeight
       airPotentialTemperatureRestartTest     = airPotentialTemperatureRestartTest     + airPotentialTemperature
       airTemperatureRestartTest              = airTemperatureRestartTest              + airTemperature
       airSpecificHumidityRestartTest         = airSpecificHumidityRestartTest         + airSpecificHumidity
       airDensityRestartTest                  = airDensityRestartTest                  + airDensity
       shortwaveVisibleDirectDownRestartTest  = shortwaveVisibleDirectDownRestartTest  + shortwaveVisibleDirectDown
       shortwaveVisibleDiffuseDownRestartTest = shortwaveVisibleDiffuseDownRestartTest + shortwaveVisibleDiffuseDown
       shortwaveIRDirectDownRestartTest       = shortwaveIRDirectDownRestartTest       + shortwaveIRDirectDown
       shortwaveIRDiffuseDownRestartTest      = shortwaveIRDiffuseDownRestartTest      + shortwaveIRDiffuseDown
       longwaveDownRestartTest                = longwaveDownRestartTest                + longwaveDown
       rainfallRateRestartTest                = rainfallRateRestartTest                + rainfallRate
       snowfallRateRestartTest                = snowfallRateRestartTest                + snowfallRate
       uAirVelocityRestartTest                = uAirVelocityRestartTest                + uAirVelocity
       vAirVelocityRestartTest                = vAirVelocityRestartTest                + vAirVelocity
       atmosReferenceSpeed10mRestartTest      = atmosReferenceSpeed10mRestartTest      + atmosReferenceSpeed10m
       atmosReferenceTemperature2mRestartTest = atmosReferenceTemperature2mRestartTest + atmosReferenceTemperature2m
       atmosReferenceHumidity2mRestartTest    = atmosReferenceHumidity2mRestartTest    + atmosReferenceHumidity2m

       ! atmosphere coupling
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling_restart_test", oceanCouplingRestartTestPool)

       call MPAS_pool_get_array(oceanCouplingPool, "seaSurfaceTemperature", seaSurfaceTemperature)
       call MPAS_pool_get_array(oceanCouplingPool, "seaSurfaceSalinity", seaSurfaceSalinity)
       call MPAS_pool_get_array(oceanCouplingPool, "seaFreezingTemperature", seaFreezingTemperature)
       call MPAS_pool_get_array(oceanCouplingPool, "freezingMeltingPotential", freezingMeltingPotential)
       call MPAS_pool_get_array(oceanCouplingPool, "uOceanVelocity", uOceanVelocity)
       call MPAS_pool_get_array(oceanCouplingPool, "vOceanVelocity", vOceanVelocity)
       call MPAS_pool_get_array(oceanCouplingPool, "seaSurfaceTiltU", seaSurfaceTiltU)
       call MPAS_pool_get_array(oceanCouplingPool, "seaSurfaceTiltV", seaSurfaceTiltV)
       call MPAS_pool_get_array(oceanCouplingPool, "oceanMixedLayerDepth", oceanMixedLayerDepth)
       call MPAS_pool_get_array(oceanCouplingPool, "oceanHeatFluxConvergence", oceanHeatFluxConvergence)

       call MPAS_pool_get_array(oceanCouplingRestartTestPool, "seaSurfaceTemperatureRestartTest", seaSurfaceTemperatureRestartTest)
       call MPAS_pool_get_array(oceanCouplingRestartTestPool, "seaSurfaceSalinityRestartTest", seaSurfaceSalinityRestartTest)
       call MPAS_pool_get_array(oceanCouplingRestartTestPool, "seaFreezingTemperatureRestartTest", &
                                                               seaFreezingTemperatureRestartTest)
       call MPAS_pool_get_array(oceanCouplingRestartTestPool, "freezingMeltingPotentialRestartTest", &
                                                               freezingMeltingPotentialRestartTest)
       call MPAS_pool_get_array(oceanCouplingRestartTestPool, "uOceanVelocityRestartTest", uOceanVelocityRestartTest)
       call MPAS_pool_get_array(oceanCouplingRestartTestPool, "vOceanVelocityRestartTest", vOceanVelocityRestartTest)
       call MPAS_pool_get_array(oceanCouplingRestartTestPool, "seaSurfaceTiltURestartTest", seaSurfaceTiltURestartTest)
       call MPAS_pool_get_array(oceanCouplingRestartTestPool, "seaSurfaceTiltVRestartTest", seaSurfaceTiltVRestartTest)
       call MPAS_pool_get_array(oceanCouplingRestartTestPool, "oceanMixedLayerDepthRestartTest", oceanMixedLayerDepthRestartTest)
       call MPAS_pool_get_array(oceanCouplingRestartTestPool, "oceanHeatFluxConvergenceRestartTest", &
                                                               oceanHeatFluxConvergenceRestartTest)

       seaSurfaceTemperatureRestartTest    = seaSurfaceTemperatureRestartTest    + seaSurfaceTemperature
       seaSurfaceSalinityRestartTest       = seaSurfaceSalinityRestartTest       + seaSurfaceSalinity
       seaFreezingTemperatureRestartTest   = seaFreezingTemperatureRestartTest   + seaFreezingTemperature
       freezingMeltingPotentialRestartTest = freezingMeltingPotentialRestartTest + freezingMeltingPotential
       uOceanVelocityRestartTest           = uOceanVelocityRestartTest           + uOceanVelocity
       vOceanVelocityRestartTest           = vOceanVelocityRestartTest           + vOceanVelocity
       seaSurfaceTiltURestartTest          = seaSurfaceTiltURestartTest          + seaSurfaceTiltU
       seaSurfaceTiltVRestartTest          = seaSurfaceTiltVRestartTest          + seaSurfaceTiltV
       oceanMixedLayerDepthRestartTest     = oceanMixedLayerDepthRestartTest     + oceanMixedLayerDepth
       oceanHeatFluxConvergenceRestartTest = oceanHeatFluxConvergenceRestartTest + oceanHeatFluxConvergence

       block => block % next
    enddo

  end subroutine cice_debug_forcing_restart_test

!-----------------------------------------------------------------------

end module cice_debug
