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
       cice_debug_output_array

  integer, private :: &
       iObject = 2, &
       iLabel  = 1, &
       iArrow  = 1

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
    write(55,*) "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0, 4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0 )"

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
                write(strout,fmt='(a,i5,a,e14.6,a)') "set object ",iObject,' fc palette cb ', cellArray(iVertex), ' fillstyle solid'
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

             write(strout,fmt='(a,i5,a,e14.6,a,e14.6,a,e14.6,a,e14.6,a)') "set arrow ",iArrow," from ", y1, ", ", z1, " to ", y2, ", ", z2, " nohead lt 4"
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

             write(strout,fmt='(a,i5,a,e14.6,a,e14.6,a,e14.6,a,e14.6,a)') "set arrow ",iArrow," from ", y1, ", ", z1, " to ", y2, ", ", z2, " nohead lt 3"
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
    write(55,*) "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0, 4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0 )"

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

             write(strout,fmt='(a,i5,a,f14.2,a,f14.2,a,f14.2,a,f14.2,a)') "set arrow ",iArrow," from ", y1, ", ", z1, " to ", y2, ", ", z2, " nohead lt 4"
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

             write(strout,fmt='(a,i5,a,f14.2,a,f14.2,a,f14.2,a,f14.2,a)') "set arrow ",iArrow," from ", y1, ", ", z1, " to ", y2, ", ", z2, " nohead lt 3"
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
       POP_nx_2 = nint(real(POP_nx) / 2.0)

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
    write(55,*) "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0, 4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0 )"

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
    write(55,*) "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0, 4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0 )"

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

    use mpas_constants, only: &
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
    write(55,*) "set palette defined (0 0.0 0.0 0.5, 1 0.0 0.0 1.0, 2 0.0 0.5 1.0, 3 0.0 1.0 1.0, 4 0.5 1.0 0.5, 5 1.0 1.0 0.0, 6 1.0 0.5 0.0, 7 1.0 0.0 0.0, 8 0.5 0.0 0.0 )"

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

  subroutine cice_test_case_output(block, itimestep)!{{{

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    integer, intent(in) :: &
         itimestep !< Input: 

    type (MPAS_pool_type), pointer :: &
         configs

    logical, pointer :: &
         config_test_case_diag

    character(len=strKIND), pointer :: &
         config_test_case_diag_type

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
         icestate

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
    call MPAS_pool_get_subpool(block % structs, "icestate", icestate)

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

    call MPAS_pool_get_array(icestate, "iceAreaCell", iceAreaCell)    

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

    write(*,*) "stressDivergenceU: ", sum(abs(stressDivergenceUGeographical(1:nVertices)), mask=interiorVertex(1:nVertices)==1), nVertices, sum(interiorVertex)
    write(*,*) "stressDivergenceV: ", sum(abs(stressDivergenceVGeographical(1:nVertices)), mask=interiorVertex(1:nVertices)==1), nVertices, sum(interiorVertex)

    ! ice velocity
    call plot_vertices_spherical(mesh, uVelocityGeographical, interiorVertex, "uVelocity.txt", itimestep)
    call plot_vertices_spherical(mesh, vVelocityGeographical, interiorVertex, "vVelocity.txt", itimestep)

    write(*,*) "uVelocity: ", sum(abs(uVelocityGeographical(1:nVertices)), mask=interiorVertex(1:nVertices)==1), nVertices, sum(interiorVertex)
    write(*,*) "vVelocity: ", sum(abs(vVelocityGeographical(1:nVertices)), mask=interiorVertex(1:nVertices)==1), nVertices, sum(interiorVertex)

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
         plotting

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
         iceThickness, &
         snowThickness

    integer :: &
         iCell

    integer, pointer :: &
         nCells

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "plotting", plotting)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    call MPAS_pool_get_array(plotting, "iceAreaCategoryPlot", iceAreaCategoryPlot)
    call MPAS_pool_get_array(plotting, "iceVolumeCategoryPlot", iceVolumeCategoryPlot)
    call MPAS_pool_get_array(plotting, "snowVolumeCategoryPlot", snowVolumeCategoryPlot)
    call MPAS_pool_get_array(plotting, "surfaceTemperaturePlot", surfaceTemperaturePlot)
    call MPAS_pool_get_array(plotting, "iceThickness", iceThickness)
    call MPAS_pool_get_array(plotting, "snowThickness", snowThickness)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    iceAreaCategoryPlot(:)    = iceAreaCategory(1,1,:)
    iceVolumeCategoryPlot(:)  = iceVolumeCategory(1,1,:)
    snowVolumeCategoryPlot(:) = snowVolumeCategory(1,1,:)
    surfaceTemperaturePlot(:) = surfaceTemperature(1,1,:)

    do iCell = 1, nCells

       if (iceAreaCategoryPlot(iCell) > 0.0_RKIND) then
          iceThickness(iCell)  = iceVolumeCategoryPlot(iCell)  / iceAreaCategoryPlot(iCell)
          snowThickness(iCell) = snowVolumeCategoryPlot(iCell) / iceAreaCategoryPlot(iCell)
       else
          iceThickness(iCell)  = 0.0_RKIND
          snowThickness(iCell) = 0.0_RKIND
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

  subroutine cice_debug_output_array(array)

    real(kind=RKIND), dimension(:), intent(in) :: &
         array

    integer :: recl

    recl = size(array,1) * 8

    open(60,file="mpas_array_comparison.dat",form="unformatted",access="direct",recl=recl)

    write(60,rec=1) array

    close(60)

  end subroutine cice_debug_output_array

!-----------------------------------------------------------------------

end module cice_debug

