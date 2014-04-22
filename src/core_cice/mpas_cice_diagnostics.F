module cice_diagnostics

  use mpas_grid_types

  implicit none

  private
  public :: &
       writeout_state, &
       gnuplot_cell, &
       matlab_jet, &
       plot_pop_grid_cells, &
       rms_difference, &
       plot_pop_grid_vertices, &
       gnuplot_triangle, &
       writeout_minmax, &
       plot_cells_spherical, &
       plot_vertices_spherical, &
       plot_state, &
       init_ijpop_from_ivertex, &
       output_meridonal_strip

  integer, public :: iCellTest
  
  integer :: &
       iObject = 2, &
       iLabel  = 1, &
       iArrow  = 1

  real(kind=RKIND) :: &
       xmin =  1.0e30_RKIND, &
       xmax = -1.0e30_RKIND, &
       ymin =  1.0e30_RKIND, &
       ymax = -1.0e30_RKIND

contains

  !-------------------------------------------------------------
  ! general plotting
  !-------------------------------------------------------------
  
  subroutine writeout_state(block, nstep)

    type(block_type), intent(in) :: block

    integer, intent(in) :: &
         nstep

    type (MPAS_pool_type), pointer :: &
         boundary, &
         icestate, &
         mesh, &
         hexdyn

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         iceAreaCell, &
         iceVolumeCell         

    integer, dimension(:), pointer :: &
         interiorVertex

    integer, parameter :: &
         nstep_skip = 100

    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "hexdyn", hexdyn)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)

    call MPAS_pool_get_array(icestate, "iceAreaCell", iceAreaCell)
    call MPAS_pool_get_array(icestate, "iceVolumeCell", iceVolumeCell)

    call MPAS_pool_get_array(hexdyn, "uVelocity", uVelocity)
    call MPAS_pool_get_array(hexdyn, "vVelocity", vVelocity)    
    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)    

    if (nstep == 0 .or. nstep == 1 .or. mod(nstep,nstep_skip) == 0) then

       ! cell centre quantities
       call plot_pop_grid_cells(mesh, iceAreaCell,   "iceAreaCell.txt",   nstep)
       call plot_pop_grid_cells(mesh, iceVolumeCell, "iceVolumeCell.txt", nstep)
       
       ! vertex quantities
       call plot_pop_grid_vertices(mesh, uVelocity, interiorVertex, "uVelocity.txt", nstep)
       call plot_pop_grid_vertices(mesh, vVelocity, interiorVertex, "vVelocity.txt", nstep)

    endif

  end subroutine writeout_state
  
  !-------------------------------------------------------------

  subroutine gnuplot_cell(mesh, cellArray, filename, append, nofill) 
    
    type(MPAS_pool_type), pointer, intent(in) :: mesh
    
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

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         cellsOnCell, &
         verticesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         xVertex, &
         yVertex

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "cellsOnCell", cellsOnCell)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)

    ! optional arguments
    lappend = .false.
    if (present(append)) lappend = append

    lnofill = .false.
    if (present(nofill)) lnofill = nofill

    fMin = minval(cellArray(1:nCells))
    fMax = maxval(cellArray(1:nCells))
    !fMin = 1
    !fMax = nCells

    if (lappend) then
       open(55,file=filename,action="write",position="append")
    else
       open(55,file=filename,action="write")
    end if

    open(56,file="border.txt",action="write")

    write(56,*) minval(cellsOnCell), maxval(cellsOnCell)
    write(56,*) minloc(cellsOnCell), maxloc(cellsOnCell)
    write(56,*) size(cellsOnCell,1), size(cellsOnCell,2)

    do iCell = 1, nCells

       xc = xCell(iCell)
       yc = yCell(iCell)

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

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          !write(*,*) iCell, iVertexOnCell

          write(56,*) iCell, nCells, iVertexOnCell, cellsOnCell(iVertexOnCell, iCell)

          iVertex = verticesOnCell(iVertexOnCell, iCell)

          x = xVertex(iVertex)
          y = yVertex(iVertex)

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

  end subroutine gnuplot_cell

  !-------------------------------------------------------------

  subroutine gnuplot_triangle(mesh, vertexArray, filename, append, nofill)

    type(MPAS_pool_type), pointer, intent(in) :: mesh

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

    integer, pointer :: &
         nCells, &
         nVertices, &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         xCell, &
         yCell

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)

    ! optional arguments
    lappend = .false.
    if (present(append)) lappend = append

    lnofill = .false.
    if (present(nofill)) lnofill = nofill

    ! min/maxes
    fMin = 1e30
    fMax = -1e30

    ! loop over triangles
    do iVertex = 1, nVertices
       atBoundary = .false.

       ! loop over triangle vertices
       do iVertexDegree = 1, vertexDegree

          iCell = cellsOnVertex(iVertexDegree,iVertex)

          if (iCell > nCells) then
             atBoundary = .true.
          end if
       enddo ! iVertexDegree

       if (.not. atBoundary) then
          fMin = min(fMin,vertexArray(iVertex))
          fMax = max(fMax,vertexArray(iVertex))
       endif
    
    enddo ! iVertex
    write(*,*) "gt:", fMin, fMax! ; stop
    !fMin = -0.212_RKIND
    !fMax = 0.033_RKIND

    !fMin = minval(vertexArray)
    !fMax = maxval(vertexArray)
    !fMin = 1
    !fMax = nVertices

    if (lappend) then
       open(55,file=filename,action="write",position="append")
    else
       open(55,file=filename,action="write")
    end if

    ! loop over triangles
    do iVertex = 1, nVertices

       xv = xVertex(iVertex)
       yv = yVertex(iVertex)

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
       do iVertexDegree = 1, vertexDegree

          iCell = cellsOnVertex(iVertexDegree,iVertex)

          if (iCell > nCells) then
             atBoundary = .true.
          end if

          x = xCell(iCell)
          y = yCell(iCell)

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

    type(MPAS_pool_type), pointer, intent(in) :: mesh

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

    integer, pointer :: &
         nVertices

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex

    ! init variables
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)

    ! find maximum vector
    arrowLengthMax = -1e30

    do iVertex = 1, nVertices

       n1 = vertexArrayU(iVertex)
       n2 = vertexArrayV(iVertex)
       arrowLength = sqrt(n1**2 + n2**2)

       arrowLengthMax = Max(arrowLengthMax,arrowLength)

    enddo ! iVertex

    write(*,*) arrowLengthMax


    open(55,file=filename)

    do iVertex = 1, nVertices
    
       random = rand()

       if (random <= 1.0_RKIND / real(plot_every,RKIND)) then


          n1 = vertexArrayU(iVertex)
          n2 = vertexArrayV(iVertex)
          arrowLength = sqrt(n1**2 + n2**2)
          n1 = n1 / arrowLength
          n2 = n2 / arrowLength

          arrowLengthPlot = (arrowLength / arrowLengthMax) * plotLength

          x1 = xVertex(iVertex)
          y1 = yVertex(iVertex)

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

  subroutine plot_boundary_triangles(mesh, boundaryVertex, interiorVertex, boundaryEdge)

    type(MPAS_pool_type), pointer, intent(in) :: mesh

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

    integer, pointer :: &
         nCells, &
         nVertices, &
         nEdges, &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         cellsOnEdge, &
         verticesOnEdge, &
         cellsOnVertex

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         xVertex, &
         yVertex

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "nEdges", nEdges)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "verticesOnEdge", verticesOnEdge)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)

    start_index = nCells + nVertices + 10

    open(55,file="boundary_triangles_edge.txt",action="write")

    ! edge triangles
    iEdgeBoundary = 0

    do iEdge = 1, nEdges

       write(0,*) iEdge, boundaryEdge(iEdge)

       if (boundaryEdge(iEdge) == 1) then
          ! am on boundary
          iEdgeBoundary = iEdgeBoundary + 1

          write(strout,fmt='(a,i5,a)') "set object ",iEdgeBoundary+start_index," polygon from "

          ! find the interior cell
          do iCellOnEdge = 1, 2

             iCell = cellsOnEdge(iCellOnEdge,iEdge)

             if (iCell >= 1 .and. iCell <= nCells) then

                ! this cell is interior to domain
                x = xCell(iCell)
                y = yCell(iCell)               
                x0 = x
                y0 = y

                write(stroutvertex,fmt='(f10.2,a,f10.2,a)') x, ",", y, " to "

                strout = trim(strout)//trim(stroutvertex)

             endif

          enddo ! iCellOnEdge

          ! find the adjacent vertices
          do iVertexOnEdge = 1, 2

             iVertex = verticesOnEdge(iVertexOnEdge,iEdge)

             x = xVertex(iVertex)
             y = yVertex(iVertex)

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
    do iVertex = 1, nVertices

       ! boundary vertex adjacent to two interior cells
       if (boundaryVertex(iVertex) == 1) then

          iEdgeBoundary = iEdgeBoundary + 1

          write(strout,fmt='(a,i5,a)') "set object ",iEdgeBoundary+start_index," polygon from "

          ! starting vertex point
          x = xVertex(iVertex)
          y = yVertex(iVertex)
          x0 = x
          y0 = y

          write(stroutvertex,fmt='(f10.2,a,f10.2,a)') x, ",", y, " to "

          strout = trim(strout)//trim(stroutvertex)

          ! the two interior cells
          do iVertexDegree = 1, vertexDegree

             ! adjacent cell
             iCell = cellsOnVertex(iVertexDegree, iVertex)

             if (iCell >= 1 .and. iCell <= nCells) then
                ! interior adjacent cell

                x = xCell(iCell)
                y = yCell(iCell)

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

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    integer, dimension(:,:), intent(in) :: &
         verticesOnCell

    character(len=*), intent(in) :: &
         filename

    integer :: &
         iCell

    integer, pointer :: &
         nCells

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)

    do iCell = 1, nCells

       call label_verticesOnCell(mesh, iCell, verticesOnCell(:,iCell), filename, .true.) 

    enddo ! iCell

  end subroutine labels_verticesOnCell

  !-------------------------------------------------------------

  subroutine label_verticesOnCell(mesh, iCell, verticesOnCell0, filename, append) 

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    integer, intent(in) :: &
         iCell

    integer, dimension(:), intent(in) :: &
         verticesOnCell0

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

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         xVertex, &
         yVertex

    ! init variables
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)

    lappend = .false.
    if (present(append)) lappend = .true.

    if (lappend) then
       open(55,file=filename,position="append")
    else
       open(55,file=filename)
    endif

    write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "

    xc = xCell(iCell)
    yc = yCell(iCell)

    ! draw cell
    do iVertexOnCell = 1, nEdgesOnCell(iCell)

       iVertex = verticesOnCell(iVertexOnCell,iCell)

       x = xVertex(iVertex)
       y = yVertex(iVertex)

       if (iVertexOnCell == 1) then
          x0 = x
          y0 = y
       endif

       write(stroutvertex,fmt='(f10.2,a,f10.2,a)') x, ",", y, " to "

       strout = trim(strout)//trim(stroutvertex)

       xl = 0.3_RKIND * xc + 0.7_RKIND * x
       yl = 0.3_RKIND * yc + 0.7_RKIND * y

       !write(stroutint,fmt='(i5)') iVertexOnCell
       write(stroutint,fmt='(i5)') verticesOnCell0(iVertexOnCell)
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

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    integer, dimension(:,:), intent(in) :: &
         verticesDegreeOnVertex

    integer, dimension(:), intent(in) :: &
         interiorVertex

    character(len=*), intent(in) :: &
         filename

    integer :: &
         iVertex

    integer, pointer :: &
         nVertices

    ! init variables
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)

    do iVertex = 1, nVertices

       if (interiorVertex(iVertex) == 1) then

          call label_verticesDegreeOnVertex(mesh, iVertex, verticesDegreeOnVertex(:,iVertex), filename, .true.) 

       endif

    enddo ! iVertex

  end subroutine labels_verticesDegreeOnVertex

  !-------------------------------------------------------------

  subroutine label_verticesDegreeOnVertex(mesh, iVertex, verticesDegreeOnVertex, filename, append) 

    type(MPAS_pool_type), pointer, intent(in) :: mesh

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

    integer, pointer :: &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         cellsOnVertex
 
    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         xCell, &
         yCell

    ! init variables
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)

    lappend = .false.
    if (present(append)) lappend = .true.

    if (lappend) then
       open(55,file=filename,position="append")
    else
       open(55,file=filename)
    endif

    write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "

    xv = xVertex(iVertex)
    yv = yVertex(iVertex)

    ! draw cell
    do iVertexDegree = 1, vertexDegree

       iCell = cellsOnVertex(iVertexDegree,iVertex)

       x = xCell(iCell)
       y = yCell(iCell)

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

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(in) :: &
         arrayOnCell
    
    character(len=*), intent(in) :: &
         filename

    integer :: iCell

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)

    open(55,file=filename)

    do iCell = 1, nCells

       write(55,*) iCell, arrayOnCell(iCell), &
            xCell(iCell), yCell(iCell)

    enddo ! iCell

    close(55)

  end subroutine writeout_on_cell_real

  !-------------------------------------------------------------

  subroutine writeout_on_vertex_real(mesh, arrayOnVertex, filename)

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(in) :: &
         arrayOnVertex
    
    character(len=*), intent(in) :: &
         filename

    integer :: iVertex

    integer, pointer :: &
         nVertices

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex

    open(55,file=filename)

    do iVertex = 1, nVertices

       write(55,*) iVertex, arrayOnVertex(iVertex), &
            xVertex(iVertex), yVertex(iVertex)

    enddo ! iVertex

    close(55)

  end subroutine writeout_on_vertex_real

  !-------------------------------------------------------------

  subroutine writeout_on_vertex_on_cell_real(mesh, arrayVertexOnCell, filename)

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         arrayVertexOnCell
    
    character(len=*), intent(in) :: &
         filename

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

    open(55,file=filename)

    do iCell = 1, nCells

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          iVertex = verticesOnCell(iVertexOnCell,iCell)

          write(55,*) iCell, iVertexOnCell, arrayVertexOnCell(iVertexOnCell,iCell), &
               xVertex(iVertex), yVertex(iVertex), &
               arrayVertexOnCell(iVertexOnCell,iCell) / (xVertex(iVertex) + yVertex(iVertex))

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

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         arrayEdgeCell
    
    character(len=*), intent(in) :: &
         filename

    integer :: &
         iCell, &
         iEdge

    integer, pointer :: &
         nCells
    
    integer, dimension(:), pointer :: &
         nEdgesOnCell

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)

    open(55,file=filename)

    do iCell = 1, nCells
       
       do iEdge = 1, nEdgesOnCell(iCell)

          write(55,*) iCell, iEdge, arrayEdgeCell(iEdge,iCell)

       enddo ! iEdge

    enddo ! iCell

    close(55)

  end subroutine writeout_edgecell

  !-------------------------------------------------------------

  subroutine writeout_edgeedgecell(mesh, arrayEdgeEdgeCell, filename)

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         arrayEdgeEdgeCell
    
    character(len=*), intent(in) :: &
         filename

    integer :: &
         iCell, &
         iEdge, &
         jEdge

    integer, pointer :: &
         nCells
    
    integer, dimension(:), pointer :: &
         nEdgesOnCell

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)

    open(55,file=filename)

    do iCell = 1, nCells
       
       do iEdge = 1, nEdgesOnCell(iCell)

          do jEdge = 1, nEdgesOnCell(iCell) 

             write(55,*) iCell, iEdge, jEdge, arrayEdgeEdgeCell(jEdge,iEdge,iCell)

          enddo ! jEdge

       enddo ! iEdge

    enddo ! iCell

    close(55)

  end subroutine writeout_edgeedgecell

  !-------------------------------------------------------------

  function find_nearest_cell(mesh,x,y) result(iNear)

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), intent(in) :: &
         x, y

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
    call MPAS_pool_get_array(mesh, "nCells", nCells)
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

  end function find_nearest_cell

  !-------------------------------------------------------------

  subroutine rotate_ninety(mesh)

    type(MPAS_pool_type), pointer, intent(in) :: mesh

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
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "nEdges", nEdges)
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

  end subroutine rotate_ninety

  !-------------------------------------------------------------
  ! Spherical plotting
  !-------------------------------------------------------------

  subroutine plot_state(mesh, icestate, atmos_forcing, time, dt)

    type(MPAS_pool_type), pointer, intent(in) :: mesh
    type(MPAS_pool_type), pointer :: icestate
    type(MPAS_pool_type), pointer :: atmos_forcing

    real(kind=RKIND), intent(in) :: &
         time, &
         dt

    character(len=200) :: &
         suffix

    integer :: &
         iCell

    integer, parameter :: &
         output_freq = 24

    real(kind=RKIND), dimension(:), allocatable :: &
         iceThickness, &
         snowThickness

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell

    real(kind=RKIND), dimension(:), pointer :: &
         surfaceTemperature, &
         airPotentialTemperature, &
         longwaveDown

    ! init variables
    call MPAS_pool_get_array(icestate, "iceAreaCell", iceAreaCell)
    call MPAS_pool_get_array(icestate, "iceVolumeCell", iceVolumeCell)
    call MPAS_pool_get_array(icestate, "snowVolumeCell", snowVolumeCell)
    call MPAS_pool_get_array(atmos_forcing, "surfaceTemperature", surfaceTemperature)
    call MPAS_pool_get_array(atmos_forcing, "airPotentialTemperature", airPotentialTemperature)
    call MPAS_pool_get_array(atmos_forcing, "longwaveDown", longwaveDown)

    if (mod(nint(time/dt),output_freq) == 0) then

       allocate(iceThickness(nCells))
       allocate(snowThickness(nCells))

       do iCell = 1, nCells

          if (iceAreaCell(iCell) > 0.0_RKIND) then
             iceThickness(iCell)  = iceVolumeCell(iCell) / &
                                    iceAreaCell(iCell)
             snowThickness(iCell) = snowVolumeCell(iCell) / &
                                    iceAreaCell(iCell)
          else
             iceThickness(iCell)  = 0.0_RKIND
             snowThickness(iCell) = 0.0_RKIND
          endif

          if (iCell == 7755) write(*,*) time, iceThickness(iCell), snowThickness(iCell), surfaceTemperature(iCell)

       enddo ! iCell

       !return

       write(suffix,fmt='(a,i10.10,a)') "_", nint(time), ".txt"
       
       call plot_cells_spherical(mesh, airPotentialTemperature, "movie/data/airTemperature"//trim(suffix), -50.0_RKIND, 50.0_RKIND)
       call plot_cells_spherical(mesh, longwaveDown, "movie/data/longwave"//trim(suffix), 0.0_RKIND, 1000.0_RKIND)

       call plot_cells_spherical(mesh, iceAreaCell, "movie/data/iceAreaCell"//trim(suffix),          0.0_RKIND, 1.0_RKIND)
       call plot_cells_spherical(mesh, iceThickness,       "movie/data/iceThickness"//trim(suffix),         0.0_RKIND, 2.0_RKIND)
       call plot_cells_spherical(mesh, snowThickness,      "movie/data/snowThickness"//trim(suffix),        0.0_RKIND, 2.0_RKIND)
       call plot_cells_spherical(mesh, surfaceTemperature, "movie/data/surfaceTemperature"//trim(suffix), -60.0_RKIND, 0.1_RKIND)

       deallocate(iceThickness)
       deallocate(snowThickness)

       !stop

    endif

  end subroutine plot_state

  !-------------------------------------------------------------

  subroutine plot_vertices_spherical(mesh, cellArray, interiorVertex, filename, fmin, fmax)

    type(MPAS_pool_type), pointer, intent(in) :: mesh
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         cellArray

    integer, dimension(:), intent(in) :: &
         interiorVertex

    character(len=*), optional, intent(in) :: &
         filename
    
    real(kind=RKIND), optional, intent(in) :: &
         fmin, &
         fmax

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
         fMin_use, &
         fMax_use

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

    character(len=7) :: &
         color

    character(len=400) :: &
         stroutvertex, &
         strout, &
         stroutlabel, &
         stroutint

    character(len=200) :: &
         filename_use

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
    
    ! init variables
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "sphere_radius", sphere_radius)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    iObject = 1

    lnofill = .false.

    if (present(fmin) .and. present(fmax)) then
       fMin_use = fmin
       fMax_use = fmax
    else
       fMin_use = minval(cellArray(1:nVertices))
       fMax_use = maxval(cellArray(1:nVertices))
    endif

    filename_use = "spherical_vertices.txt"
    if (present(filename)) filename_use = trim(filename)

    open(55,file=trim(filename_use),action="write")
    open(56,file="clearPlot.gnuplot",action="write")

    do iVertex = 1, nVertices

       if (interiorVertex(iVertex) == 1) then

          xc = xVertex(iVertex)
          yc = yVertex(iVertex)
          zc = zVertex(iVertex)
          
          call perform_rotation(xc, yc, zc)
          
          if (xc > 0.0_RKIND) then
             
             fValue = (cellArray(iVertex) - fMin_use) / (fMax_use - fMin_use)
             !fValue = (real(iCell,RKIND) - fMin_use) / (fMax_use - fMin_use)
             
             write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "
             
             do iVertexDegree = 1, vertexDegree
                
                iCell = cellsOnVertex(iVertexDegree, iVertex)
                
                xv = xCell(iCell)
                yv = yCell(iCell)
                zv = zCell(iCell)
                
                call perform_rotation(xv, yv, zv)
                
                if (iVertexDegree == 1) then
                   xv0 = xv
                   yv0 = yv
                   zv0 = zv
                endif
                
                write(stroutvertex,fmt='(f14.2,a,f14.2,a)') yv, ",", zv, " to "
                strout = trim(strout)//trim(stroutvertex)
                
             enddo ! iVertexDegree
             
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
          
       endif ! interiorVertex
       
    enddo ! iVertex
       
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
    write(stroutint,fmt='(e20.8)') fMax_use
    write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.9'
    write(55,*) trim(stroutlabel)
    write(56,*) "unset label ", iLabel
    iLabel = iLabel + 1

    write(stroutint,fmt='(e20.8)') fMin_use
    write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.8'
    write(55,*) trim(stroutlabel)
    write(56,*) "unset label ", iLabel
    iLabel = iLabel + 1

    close(55)
    close(56)

    write(strout,fmt='(a,f15.3,a,f15.3,a,f15.3,a,f15.3,a)') &
         "set xrange [",-sphere_radius*1.1_RKIND,":",sphere_radius*1.1_RKIND,&
         "] ; set yrange [",-sphere_radius*1.1_RKIND,":",sphere_radius*1.1_RKIND,"]"
    !write(*,*) trim(strout)

  end subroutine plot_vertices_spherical

  !-------------------------------------------------------------

  subroutine plot_cells_spherical(mesh, cellArray, filename, fmin, fmax)

    type(MPAS_pool_type), pointer, intent(in) :: mesh
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         cellArray

    character(len=*), optional, intent(in) :: &
         filename
    
    real(kind=RKIND), optional, intent(in) :: &
         fmin, &
         fmax

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
         fMin_use, &
         fMax_use

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

    character(len=200) :: &
         filename_use

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

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "sphere_radius", sphere_radius)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)

    iObject = 1

    lnofill = .false.

    if (present(fmin) .and. present(fmax)) then
       fMin_use = fmin
       fMax_use = fmax
    else
       fMin_use = minval(cellArray(1:nCells))
       fMax_use = maxval(cellArray(1:nCells))
    endif

    filename_use = "spherical_cells.txt"
    if (present(filename)) filename_use = trim(filename)

    open(55,file=trim(filename_use),action="write")
    open(56,file="clearPlot.gnuplot",action="write")

    do iCell = 1, nCells

       xc = xCell(iCell)
       yc = yCell(iCell)
       zc = zCell(iCell)

       call perform_rotation(xc, yc, zc)

       if (xc > 0.0_RKIND) then

          fValue = (cellArray(iCell) - fMin_use) / (fMax_use - fMin_use)
          !fValue = (real(iCell,RKIND) - fMin_use) / (fMax_use - fMin_use)

          !write(*,*) iCell, iObject, fValue, cellArray(iCell)

          write(strout,fmt='(a,i5,a)') "set object ",iObject," polygon from "

          do iVertexOnCell = 1, nEdgesOnCell(iCell)

             iVertex = verticesOnCell(iVertexOnCell, iCell)

             xv = xVertex(iVertex)
             yv = yVertex(iVertex)
             zv = zVertex(iVertex)

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
    write(stroutint,fmt='(e20.8)') fMax_use
    write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.9'
    write(55,*) trim(stroutlabel)
    write(56,*) "unset label ", iLabel
    iLabel = iLabel + 1

    write(stroutint,fmt='(e20.8)') fMin_use
    write(stroutlabel,fmt='(a,i5,a,a,a)') "set label ",iLabel,' "',trim(adjustl(stroutint)),'" at graph 1.02, 0.8'
    write(55,*) trim(stroutlabel)
    write(56,*) "unset label ", iLabel
    iLabel = iLabel + 1

    close(55)
    close(56)

    write(strout,fmt='(a,f15.3,a,f15.3,a,f15.3,a,f15.3,a)') &
         "set xrange [",-sphere_radius*1.1_RKIND,":",sphere_radius*1.1_RKIND,&
         "] ; set yrange [",-sphere_radius*1.1_RKIND,":",sphere_radius*1.1_RKIND,"]"
    !write(*,*) trim(strout)

  end subroutine plot_cells_spherical

  !-------------------------------------------------------------

  subroutine latlon_to_xyz(mesh, x, y, z, lat, lon)

    use cice_constants, only: deg_to_rad

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), intent(in) :: &
         lat, lon

    real(kind=RKIND), intent(out) :: &
         x, y, z

    real(kind=RKIND), pointer :: &
         sphere_radius

    ! init variables
    call MPAS_pool_get_array(mesh, "sphere_radius", sphere_radius)

    z = sphere_radius * sin(lat*deg_to_rad)
    x = sphere_radius * cos(lon*deg_to_rad) * cos(lat*deg_to_rad)
    y = sphere_radius * sin(lon*deg_to_rad) * cos(lat*deg_to_rad)

  end subroutine latlon_to_xyz

  !-------------------------------------------------------------

  subroutine perform_rotation(x, y, z)

    use cice_constants, only: deg_to_rad

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
    thetay = 60.0_RKIND
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

  subroutine plot_pop_grid_cells(mesh, cellArray, filename, tag)

    type(MPAS_pool_type), pointer, intent(in) :: mesh
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         cellArray

    character(len=*), optional, intent(in) :: &
         filename

    integer, optional, intent(in) :: &
         tag

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
    call MPAS_pool_get_array(mesh, "POP_nx", POP_nx)
    call MPAS_pool_get_array(mesh, "POP_ny", POP_ny)
    call MPAS_pool_get_array(mesh, "nCells", nCells)
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

  end subroutine plot_pop_grid_cells

  !-------------------------------------------------------------

  subroutine plot_pop_grid_vertices(mesh, vertexArray, interiorVertex, filename, tag)

    type(MPAS_pool_type), pointer, intent(in) :: mesh
    
    real(kind=RKIND), dimension(:), intent(in) :: &
         vertexArray

    integer, dimension(:), intent(in) :: &
         interiorVertex

    character(len=*), optional, intent(in) :: &
         filename

    integer, optional, intent(in) :: &
         tag

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
    call MPAS_pool_get_array(mesh, "POP_nx", POP_nx)
    call MPAS_pool_get_array(mesh, "POP_ny", POP_ny)
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
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

  end subroutine plot_pop_grid_vertices

  !-------------------------------------------------------------

  subroutine init_ijpop_from_ivertex(mesh)

    type(MPAS_pool_type), intent(inout) :: mesh

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
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "POP_nx", POP_nx)
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
    
    open(11,file="vertexareas.txt")
    do iVertex = 1, nVertices
       write(11,*) POPindxiv(iVertex), POPindxjv(iVertex), &
            areaTriangle(iVertex)
    enddo
    close(11)
    
    open(11,file="cellareas.txt")
    do iCell = 1, nCells
       write(11,*) POPindxi(iCell), POPindxj(iCell), &
            areaCell(iCell)
    enddo
    close(11)

    open(11,file="vertexlat.txt")
    do iVertex = 1, nVertices
       write(11,*) POPindxiv(iVertex), POPindxjv(iVertex), &
            latVertex(iVertex)
    enddo
    close(11)


  end subroutine init_ijpop_from_ivertex

  !-------------------------------------------------------------

  subroutine output_meridonal_strip(mesh, array, filename, imeridonal, istep, step_freq)

    type(MPAS_pool_type), intent(inout) :: mesh

    real(kind=RKIND), dimension(:), intent(in) :: &
         array

    character(len=*), intent(in) :: &
         filename

    integer, intent(in) :: &
         imeridonal, &
         istep, &
         step_freq

    integer :: &
         iCell, &
         j, &
         ny

    real(kind=RKIND), allocatable, dimension(:) :: &
         stripval

    integer, allocatable, dimension(:) :: &
         stripmask, &
         stripcell

    integer, pointer :: &
         nCells
    
    integer, dimension(:), pointer :: &
         POPindxi, &
         POPindxj

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "POPindxi", POPindxi)
    call MPAS_pool_get_array(mesh, "POPindxj", POPindxj)

    if (mod(istep,step_freq) == 0) then
       
       ny = maxval(POPindxj)

       allocate(stripval(ny))
       allocate(stripmask(ny))
       allocate(stripcell(ny))

       stripval  = -999.0_RKIND
       stripmask = 0

       do iCell = 1, nCells

          if (POPindxi(iCell) == imeridonal) then

             j = POPindxj(iCell)

             stripval(j) = array(iCell)
             stripmask(j) = 1
             stripcell(j) = iCell

          endif

       enddo ! i

       open(55,file=trim(filename), position="append")

       do j = 1, ny

          if (stripmask(j) == 1) then

             write(55,*) j, stripcell(j), stripval(j)

          endif

       enddo ! j

       write(55,*)

       close(55)

       deallocate(stripval)
       deallocate(stripmask)
       deallocate(stripcell)

    endif

  end subroutine output_meridonal_strip

  !-------------------------------------------------------------
  ! RMS comparison
  !-------------------------------------------------------------

  subroutine rms_difference(array1, array2, mask, n, rms)

    real(kind=RKIND), dimension(:), intent(in) :: &
         array1, &
         array2

    integer, dimension(:), intent(in) :: &
         mask

    integer, intent(in) :: &
         n

    real(kind=RKIND), intent(out) :: &
         rms

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

  end subroutine rms_difference

  !-------------------------------------------------------------
  ! global statistics
  !-------------------------------------------------------------

  subroutine global_statistics(block, nstep)

    type(block_type), intent(inout) :: block

    integer, intent(in) :: &
         nstep

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell

    type (MPAS_pool_type), pointer :: &
         mesh, &
         icestate

    integer :: &
         iCell

    ! (1): Global, (2): Northern hemipshere, (3): Southern hemisphere
    real(kind=RKIND), dimension(3) :: &
         iceAreaCellMin,       iceAreaCellMax, &
         iceVolumeCellMin,     iceVolumeCellMax, &
         snowVolumeCellMin,    snowVolumeCellMax, &
         totalIceVolume

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell, &
         latCell

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "icestate", icestate)

    call MPAS_pool_get_array(icestate, "iceAreaCell", iceAreaCell)
    call MPAS_pool_get_array(icestate, "iceVolumeCell", iceVolumeCell)
    call MPAS_pool_get_array(icestate, "snowVolumeCell", snowVolumeCell)

    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)
    call MPAS_pool_get_array(mesh, "latCell", latCell)

    ! init min/max
    iceAreaCellMin(:)       =  1.0e30_RKIND
    iceAreaCellMax(:)       = -1.0e30_RKIND
    iceVolumeCellMin(:)     =  1.0e30_RKIND
    iceVolumeCellMax(:)     = -1.0e30_RKIND
    snowVolumeCellMin(:)    =  1.0e30_RKIND
    snowVolumeCellMax(:)    = -1.0e30_RKIND

    ! init totals
    totalIceVolume(:) = 0.0_RKIND

    ! determine aggregate quantities
    do iCell = 1, nCells

       ! whole globe

       ! min/max area volume
       iceAreaCellMin(1)    = min(iceAreaCellMin(1),    iceAreaCell(iCell))
       iceAreaCellMax(1)    = max(iceAreaCellMin(1),    iceAreaCell(iCell))
       iceVolumeCellMin(1)  = min(iceVolumeCellMin(1),  iceVolumeCell(iCell))
       iceVolumeCellMax(1)  = max(iceVolumeCellMin(1),  iceVolumeCell(iCell))
       snowVolumeCellMin(1) = min(snowVolumeCellMin(1), snowVolumeCell(iCell))
       snowVolumeCellMax(1) = max(snowVolumeCellMin(1), snowVolumeCell(iCell))

       ! total ice volume
       totalIceVolume(1) = totalIceVolume(1) + areaCell(iCell) * iceVolumeCell(iCell)

       if (latCell(iCell) > 0.0_RKIND) then

          ! northern hemisphere

          ! min/max area volume
          iceAreaCellMin(2)    = min(iceAreaCellMin(2),    iceAreaCell(iCell))
          iceAreaCellMax(2)    = max(iceAreaCellMin(2),    iceAreaCell(iCell))
          iceVolumeCellMin(2)  = min(iceVolumeCellMin(2),  iceVolumeCell(iCell))
          iceVolumeCellMax(2)  = max(iceVolumeCellMin(2),  iceVolumeCell(iCell))
          snowVolumeCellMin(2) = min(snowVolumeCellMin(2), snowVolumeCell(iCell))
          snowVolumeCellMax(2) = max(snowVolumeCellMin(2), snowVolumeCell(iCell))

          ! total ice volume
          totalIceVolume(2) = totalIceVolume(2) + areaCell(iCell) * iceVolumeCell(iCell)

       else if (latCell(iCell) < 0.0_RKIND) then

          ! southern hemisphere

          ! min/max area volume
          iceAreaCellMin(3)    = min(iceAreaCellMin(3),    iceAreaCell(iCell))
          iceAreaCellMax(3)    = max(iceAreaCellMin(3),    iceAreaCell(iCell))
          iceVolumeCellMin(3)  = min(iceVolumeCellMin(3),  iceVolumeCell(iCell))
          iceVolumeCellMax(3)  = max(iceVolumeCellMin(3),  iceVolumeCell(iCell))
          snowVolumeCellMin(3) = min(snowVolumeCellMin(3), snowVolumeCell(iCell))
          snowVolumeCellMax(3) = max(snowVolumeCellMin(3), snowVolumeCell(iCell))

          ! total ice volume
          totalIceVolume(3) = totalIceVolume(3) + areaCell(iCell) * iceVolumeCell(iCell)

       endif

    enddo ! iCell

    open(55,file="global_statistics.txt")

    write(55,*) nstep, &
         iceAreaCellMin(1),    iceAreaCellMax(1), &
         iceVolumeCellMin(1),  iceVolumeCellMax(1), &
         snowVolumeCellMin(1), snowVolumeCellMax(1), &
         totalIceVolume(1)

    close(55)

    open(55,file="global_statistics_nh.txt")

    write(55,*) nstep, &
         iceAreaCellMin(2),    iceAreaCellMax(2), &
         iceVolumeCellMin(2),  iceVolumeCellMax(2), &
         snowVolumeCellMin(2), snowVolumeCellMax(2), &
         totalIceVolume(2)

    close(55)

    open(55,file="global_statistics_sh.txt")

    write(55,*) nstep, &
         iceAreaCellMin(3),    iceAreaCellMax(3), &
         iceVolumeCellMin(3),  iceVolumeCellMax(3), &
         snowVolumeCellMin(3), snowVolumeCellMax(3), &
         totalIceVolume(3)

    close(55)

  end subroutine global_statistics

  !-------------------------------------------------------------

end module cice_diagnostics

