#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module zoltan_mod
  use kinds, only : iulog, real_kind
  use parallel_mod, only : abortmp
  use params_mod,             only : SFCURVE, ZOLTAN2RCB, ZOLTAN2MJ, &
                                       ZOLTAN2RIB, ZOLTAN2HSFC, ZOLTAN2PATOH, &
                                       ZOLTAN2PHG, ZOLTAN2METIS, &
                                       ZOLTAN2PARMETIS, ZOLTAN2SCOTCH, ZOLTAN2PTSCOTCH, &
                                       ZOLTAN2BLOCK, ZOLTAN2CYCLIC, ZOLTAN2RANDOM, &
                                       ZOLTAN2ZOLTAN, ZOLTAN2ND, ZOLTAN2PARMA, &
                                       ZOLTAN2MJRCB, ZOLTAN2_1PHASEMAP,  &
                                       Z2_NO_TASK_MAPPING, Z2_TASK_MAPPING, &
                                       Z2_OPTIMIZED_TASK_MAPPING
  implicit none

  private 
  integer, parameter :: VertexWeight = 1
  integer, parameter :: EdgeWeight = 1

  public :: genzoltanpart, getfixmeshcoordinates, printMetrics, is_zoltan_partition, is_zoltan_task_mapping

contains

  function is_zoltan_partition(partmethod) result (zm)
  integer :: partmethod
  logical :: zm

  zm=.false.
  if (partmethod .eq. ZOLTAN2RCB .OR. &
       partmethod .eq. ZOLTAN2MJ .OR.  &
       partmethod .eq. ZOLTAN2RIB .OR. &
       partmethod .eq. ZOLTAN2HSFC .OR. &
       partmethod .eq. ZOLTAN2PATOH .OR. &
       partmethod .eq. ZOLTAN2PHG .OR. &
       partmethod .eq. ZOLTAN2METIS .OR. &
       partmethod .eq. ZOLTAN2PARMETIS .OR. &
       partmethod .eq. ZOLTAN2PARMA .OR. &
       partmethod .eq. ZOLTAN2SCOTCH .OR. &
       partmethod .eq. ZOLTAN2PTSCOTCH .OR. &
       partmethod .eq. ZOLTAN2BLOCK .OR. &
       partmethod .eq. ZOLTAN2CYCLIC .OR. &
       partmethod .eq. ZOLTAN2RANDOM .OR. &
       partmethod .eq. ZOLTAN2ZOLTAN .OR. &
       partmethod .eq. ZOLTAN2MJRCB .OR. &
       partmethod .eq. ZOLTAN2_1PHASEMAP .OR. &
       partmethod .eq. ZOLTAN2ND) zm=.true.
  end function is_zoltan_partition

  function is_zoltan_task_mapping(z2_map_method) result (zm)
  integer :: z2_map_method
  logical :: zm

  zm=.false.
  if ( z2_map_method .eq. Z2_TASK_MAPPING .OR. &
       z2_map_method .eq. Z2_OPTIMIZED_TASK_MAPPING ) zm=.true.
  end function is_zoltan_task_mapping



  subroutine getfixmeshcoordinates(GridVertex, coord_dim1, coord_dim2, coord_dim3, coord_dimension) !result(cartResult)

    use gridgraph_mod,          only : GridVertex_t
    use control_mod,            only : coord_transform_method
    use cube_mod,               only : cube_xstart, cube_xend, cube_ystart, cube_yend,convert_gbl_index
    use coordinate_systems_mod, only : cartesian2D_t, cartesian3D_t, cubedsphere2cart
    use dimensions_mod,         only : ne
    use params_mod,             only : SPHERE_COORDS, CUBE_COORDS, FACE_2D_LB_COORDS

    type (GridVertex_t), intent(inout) :: GridVertex(:)
    real (kind=real_kind),allocatable, intent(inout) :: coord_dim1(:)
    real (kind=real_kind),allocatable, intent(inout) :: coord_dim2(:)
    real (kind=real_kind),allocatable, intent(inout) :: coord_dim3(:)
    integer, intent(inout) :: coord_dimension
    !type(cartesian3D_t)     ,  allocatable   :: cartResult(:)
    !type (real(kind=real_kind)),  allocatable   :: coord_dim1(:)
    !type (real(kind=real_kind)),  allocatable   :: coord_dim2(:)
    !type (real(kind=real_kind)),  allocatable   :: coord_dim3(:)


    integer ie,je,face_no, nelem, i
    real (kind=real_kind)  :: dx,dy, startx, starty
    real (kind=real_kind)  :: x,y,z, centerx, centery
    type (cartesian2D_t) :: corner1, corner2, corner3, corner4
    type(cartesian3D_t)  :: cart1, cart2, cart3, cart4

    nelem = SIZE(GridVertex)

    allocate(coord_dim1(nelem))
    allocate(coord_dim2(nelem))
    allocate(coord_dim3(nelem))

    if (coord_transform_method == SPHERE_COORDS) then
        coord_dimension = 3
    elseif (coord_transform_method == CUBE_COORDS) then
        coord_dimension = 3
    elseif (coord_transform_method == FACE_2D_LB_COORDS) then
        coord_dimension = 2
    else
        coord_dimension = 3
    endif

    do i=1,nelem
      call convert_gbl_index(GridVertex(i)%number,ie,je,face_no)
      dx = (cube_xend-cube_xstart)/ne
      dy = (cube_yend-cube_ystart)/ne
      startx = cube_xstart+ie*dx
      starty = cube_ystart+je*dy


      corner1%x = startx
      corner1%y = starty
      corner2%x = startx+dx
      corner2%y = starty
      corner3%x = startx+dx
      corner3%y = starty+dy
      corner4%x = startx
      corner4%y = starty+dy
      centerx = (cube_xend-cube_xstart)/2
      centery = (cube_yend-cube_ystart)/2

      !!!!!!!!SPHERE COORDS
      if (coord_transform_method == SPHERE_COORDS) then

        !coord_dimension = 3
        cart1 = cubedsphere2cart(corner1, face_no)
        cart2 = cubedsphere2cart(corner2, face_no)
        cart3 = cubedsphere2cart(corner3, face_no)
        cart4 = cubedsphere2cart(corner4, face_no)
        x = (cart1%x + cart2%x + cart3%x + cart4%x) / 4.0
        y = (cart1%y + cart2%y + cart3%y + cart4%y) / 4.0
        z = (cart1%z + cart2%z + cart3%z + cart4%z) / 4.0

      elseif (coord_transform_method == CUBE_COORDS) then
        !!!!CUBE COORDS
        !coord_dimension = 3
        if (face_no .EQ. 1) then
            x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) * face_no)
            y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
            z = 1
        elseif (face_no .EQ. 2) then
            x = 1
            y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
            z = -(corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) * face_no)
        elseif(face_no .EQ. 3) then
            x = -(corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) * face_no)
            y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
            z = -1
        elseif(face_no .EQ. 4) then
            x = -1
            y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
            z = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) * face_no)
        elseif(face_no .EQ. 5) then
            x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) * face_no)
            y = 1
            z = (corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
        elseif(face_no .EQ. 6) then
            x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) * face_no)
            y = -1
            z = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
        endif

      elseif (coord_transform_method == FACE_2D_LB_COORDS) then
        !FACE COORDINATES each for face has an half up or down
        !coord_dimension = 2

        if (face_no .EQ. 1) then
            x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 + ((cube_xend-cube_xstart) * face_no)
            y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
            z = 1
        elseif (face_no .EQ. 2) then
            !x = 1
            y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
            !z = -(corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) * face_no)
            z = 1
            x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 + ((cube_xend-cube_xstart) * face_no)
        elseif(face_no .EQ. 3) then
            !x = -(corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) * face_no)
            !y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
            !z = -1
            x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 + ((cube_xend-cube_xstart) * face_no)
            y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
            z = 1
        elseif(face_no .EQ. 4) then
            !x = -1
            !y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
            !z = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) * face_no)
            x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 + ((cube_xend-cube_xstart) * face_no)
            y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
            z = 1
        elseif(face_no .EQ. 5) then
            x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) )
            !y = 1
            !z = (corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
            z = 1
            y = (corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart))

            if (y > centery + cube_ystart) then
                x = x + (cube_xend-cube_xstart)
                y = -y + (cube_yend-cube_ystart)
            else
                x = -x + 3 * (cube_xend-cube_xstart)
                y = y - cube_ystart + cube_yend
            endif

        elseif(face_no .EQ. 6) then
            x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) )
            !y = -1
            !z = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
            z = 1
            y = (corner1%y + corner2%y + corner3%y + corner4%y) / 4.0! - ((cube_yend-cube_ystart))


            if (x > centerx + cube_xstart) then
                z = x
                x = y
                y = z
                z = 1
                y = y - ((cube_yend-cube_ystart))
                x = x + 2 * ((cube_xend-cube_xstart) )
            else
                z = x
                x = y
                y = z
                z = 1
                x = -x + 4 * (cube_xend-cube_xstart)
                y = -y - ((cube_yend-cube_ystart))
            endif

        endif
      else
    !      if (face_no .EQ. 1) then
    !        x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 + ((cube_xend-cube_xstart) * face_no)
    !        y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 + ((cube_yend-cube_ystart) * face_no)
    !        z = 1
    !      elseif (face_no .EQ. 2) then
    !        !x = 1
    !        y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
    !        !z = -(corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) * face_no)
    !        z = 1
    !        x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 + ((cube_xend-cube_xstart) * face_no)
    !      elseif(face_no .EQ. 3) then
    !        !x = -(corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) * face_no)
    !        !y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
    !        !z = -1
    !        x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 + ((cube_xend-cube_xstart) * face_no)
    !        y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
    !        z = 1
    !      elseif(face_no .EQ. 4) then
    !        !x = -1
    !        !y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
    !        !z = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) * face_no)
    !        x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 + ((cube_xend-cube_xstart) * face_no)
    !        y = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
    !        z = 1
    !      elseif(face_no .EQ. 5) then
    !        x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) )
    !        !y = 1
    !        !z = (corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
    !        z = -1
    !        y = (corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart))
    !
    !        if (y > centery + cube_ystart) then
    !            x = x + (cube_xend-cube_xstart)
    !            y = -y + (cube_yend-cube_ystart)
    !        else
    !            x = -x + 3 * (cube_xend-cube_xstart)
    !            y = y - cube_ystart + cube_yend
    !        endif
    !
    !      elseif(face_no .EQ. 6) then
    !        x = (corner1%x + corner2%x + corner3%x + corner4%x) / 4.0 !+ ((cube_xend-cube_xstart) )
    !        !y = -1
    !        !z = -(corner1%y + corner2%y + corner3%y + corner4%y) / 4.0 !+ ((cube_yend-cube_ystart) * face_no)
    !        z = 2
    !        y = (corner1%y + corner2%y + corner3%y + corner4%y) / 4.0! - ((cube_yend-cube_ystart))
    !
    !
    !        if (Y < centerY + cube_ystart) then
    !            y = -y - ((cube_yend-cube_ystart))
    !            x = x +  ((cube_xend-cube_xstart) )
    !        else
    !            x = -x + 3 * (cube_xend-cube_xstart)
    !            y = +y - ((cube_yend-cube_ystart))
    !        endif
    !
    !      endif
      endif

      coord_dim1(i) = x
      coord_dim2(i) = y
      coord_dim3(i) = z
      !write(iulog,*) i, x,y,z
      !coord_dim1(i) = x
      !coord_dim2(i) = y
      !coord_dim3(i) = z
    end do
  end subroutine getfixmeshcoordinates


  subroutine printMetrics(GridEdge,GridVertex, comm)
    use gridgraph_mod, only : GridVertex_t, GridEdge_t
    use dimensions_mod , only : npart

    implicit none
    type (GridVertex_t), intent(inout) :: GridVertex(:)
    type (GridEdge_t),   intent(inout) :: GridEdge(:)
    integer,             intent(in) :: comm

    integer , target, allocatable :: xadj(:),adjncy(:)    ! Adjacency structure for METIS
    real(kind=REAL_KIND),  target, allocatable :: vwgt(:),adjwgt(:)    ! Weights for the adj struct for METIS
    integer                       :: nelem_edge,nelem

    nelem_edge=SIZE(GridEdge)
    nelem = SIZE(GridVertex)

    allocate(xadj(nelem+1))
    allocate(vwgt(nelem))
    allocate(adjncy(nelem_edge))
    allocate(adjwgt(nelem_edge))
    vwgt(:)=VertexWeight
    call CreateMeshGraph(GridVertex,xadj,adjncy,adjwgt)
#if TRILINOS_HAVE_ZOLTAN2
    CALL Z2PRINTMETRICS(nelem,xadj,adjncy,adjwgt,vwgt,npart, comm, GridVertex%processor_number)
#else
    call abortmp("ERROR: Zoltan partition option not available")
#endif
  end subroutine printMetrics

  subroutine genzoltanpart(GridEdge,GridVertex, comm, coord_dim1, coord_dim2, coord_dim3, coord_dimension)
    use gridgraph_mod, only : GridVertex_t, GridEdge_t, freegraph, createsubgridgraph, printgridvertex
    use kinds, only : int_kind
    use dimensions_mod , only : nmpi_per_node, npart, nnodes, nelem
    !use control_mod, only:  partmethod
    !use params_mod, only : wrecursive
    use, intrinsic :: iso_c_binding, only : C_CHAR, C_NULL_CHAR
    use control_mod, only : partmethod, z2_map_method

    implicit none 
    type (GridVertex_t), intent(inout) :: GridVertex(:)
    type (GridEdge_t),   intent(inout) :: GridEdge(:)
    integer,             intent(in) :: comm
    real (kind=real_kind),intent(in) :: coord_dim1(:)
    real (kind=real_kind),intent(in) :: coord_dim2(:)
    real (kind=real_kind),intent(in) :: coord_dim3(:)
    integer, intent(inout) :: coord_dimension


    integer , target, allocatable :: xadj(:),adjncy(:)    ! Adjacency structure for METIS
    real(kind=REAL_KIND),  target, allocatable :: vwgt(:),adjwgt(:)    ! Weights for the adj struct for METIS

    integer , target, allocatable :: xadj_nl(:),adjncy_nl(:)
    integer , target, allocatable :: vwgt_nl(:),adjwgt_nl(:)

    type (GridVertex_t), allocatable   :: SubVertex(:)

    !integer(kind=int_kind), allocatable          :: part(:)
    integer, allocatable          :: part_nl(:),local2global_nl(:)
    integer, allocatable          :: part_fl(:),local2global_fl(:)


    integer, allocatable          :: cnt(:),newnum(:),oldnum(:)

    integer                       :: nelem_edge,numflag,edgecut,wgtflag
    integer                       :: head_part,tail_part
    integer                       :: options(5)
    integer                       :: i,j,ii,ig,in,ip,if
    integer                       :: nelem_nl,nelem_fl,newPartition
    integer                       :: partitionmethod,numpartitions
    integer                       :: nodes_per_frame
    logical , parameter           :: Debug = .true.
    nelem_edge = SIZE(GridEdge) 

    !print *, "nelem = ", nelem
    !print *, "nelem_edge = ", nelem_edge
    !print *, "npart = ", npart

    !allocate(tpwgts(npart))
    !allocate(part(nelem))
    allocate(xadj(nelem+1))
    allocate(vwgt(nelem))
    allocate(adjncy(nelem_edge))
    allocate(adjwgt(nelem_edge))

    call CreateMeshGraph(GridVertex,xadj,adjncy,adjwgt)
    vwgt(:)=VertexWeight
#if TRILINOS_HAVE_ZOLTAN2
    CALL ZOLTANPART(nelem,xadj,adjncy,adjwgt,vwgt, npart, comm, coord_dim1, coord_dim2, coord_dim3,coord_dimension,  GridVertex%processor_number, partmethod, z2_map_method)
#else
    call abortmp("ERROR: Zoltan partition option not available")
#endif
  end subroutine genzoltanpart


  subroutine CreateMeshGraph(GridVertex,xadj,adjncy,adjwgt)
    use gridgraph_mod, only : GridVertex_t, num_neighbors
    use kinds, only : int_kind
    type (GridVertex_t), intent(in) :: GridVertex(:)
    integer,intent(out)           :: xadj(:),adjncy(:)
    real(kind=REAL_KIND), intent(out)           :: adjwgt(:)

    integer                         :: i,j,k,ii,jj
    integer                         :: degree,nelem
    !integer(kind=int_kind),allocatable  :: neigh_list(:),sort_indices(:)
    !real(kind=REAL_KIND),allocatable    :: neigh_wgt(:)
    integer(kind=int_kind) :: neigh_list(num_neighbors), &
                              sort_indices(num_neighbors)
    real(kind=REAL_KIND) :: neigh_wgt(num_neighbors)
    integer                         :: max_neigh

    integer :: start, cnt

    nelem = SIZE(GridVertex)

    degree = 0
    ii     = 1
    ii = 0
    max_neigh = num_neighbors

    do i=1,nelem
       !print *, "i = ", i
       xadj(i)      = ii
       degree = 0
       neigh_list=0


       do j=1,num_neighbors
          cnt = GridVertex(i)%nbrs_ptr(j+1) -  GridVertex(i)%nbrs_ptr(j) 
          start =  GridVertex(i)%nbrs_ptr(j) 
          !print *, "j,cnt,start = ", j, cnt, start
          do k=0, cnt-1
             !print *, "k,wgt = ", k, GridVertex(i)%nbrs_wgt(start+k)
             if(GridVertex(i)%nbrs_wgt(start+k) .gt. 0) then 
                degree = degree + 1
                !adjncy(ii+degree-1)   = GridVertex(i)%nbrs(start+k)
                !neigh_list(degree)    = GridVertex(i)%nbrs(start+k) - 1
                adjncy(ii+degree) = GridVertex(i)%nbrs(start+k) - 1
                !adjwgt(ii+degree-1)   = GridVertex(i)%nbrs_wgt(start+k)*EdgeWeight
                !neigh_wgt(degree)     = GridVertex(i)%nbrs_wgt(start+k)*EdgeWeight
                adjwgt(ii+degree) = GridVertex(i)%nbrs_wgt(start+k)*EdgeWeight
             endif
          enddo
       enddo
       !if (degree > max_neigh) then
       ! write(iulog,*) "degree", degree, "max_neigh", max_neigh
       ! call abortmp( "number of neighbors found exceeds expected max")
       !endif

       !call sort(degree,neigh_list,sort_indices)
       !degree       = COUNT(GridVertex(i)%nbrs_wgt(:) .gt. 0) 
       ! Copy the sorted adjncy list in

       !do j=1,degree
          !adjncy(ii+j-1) = neigh_list(sort_indices(j))
          !adjncy(ii+j) = neigh_list(sort_indices(j))
          !adjncy(ii+j) = neigh_list(j)
          !adjwgt(ii+j-1) = neigh_wgt(sort_indices(j))
          !adjwgt(ii+j) = neigh_wgt(j)
       !enddo
       !print *, "degree = ", degree
       !print *, "wgts = ", GridVertex(i)%nbrs_wgt(:) 
       ii           = ii + degree
    enddo
    xadj(nelem+1)     = ii
    !open(unit=11, file="csr.txt")
    !write(11,*) xadj
    !write(11,*) adjncy
    !write(11,*) adjwgt
    !close(11)

  end subroutine CreateMeshGraph

 subroutine sort(n,list,index)
    use kinds, only : int_kind
    implicit none
    integer, intent(in) :: n
    integer(kind=int_kind), intent(in) :: list(n)
    integer(kind=int_kind), intent(inout) :: index(n)

    ! Local variables
    integer :: i,iloc,nct,ii,lmax,lmin
    logical :: msk(n)
    logical,parameter  :: Debug =.FALSE.

    msk=.TRUE.
    if(Debug) write(iulog,*)'sort: point #1'
    nct=0
    do i=1,n
       if(list(i) .eq. 0) then 
          msk(i)=.FALSE.
       else
          nct = nct + 1
       endif
    enddo
    if(Debug) write(iulog,*)'sort: point #2   list: ',list
    if(Debug) write(iulog,*)'sort: point #2.1 msk: ',msk

    lmax=maxval(list)
    do i=1,nct
       if(Debug) write(iulog,*)'sort: point #3'

       !     pgf90 Rel 3.1-4i:  minloc() with mask is buggy   
       !     iloc = minloc(list,dim=1,mask=msk)
       lmin=lmax
       iloc=-1      
       do ii=1,n
          if (msk(ii) .and. list(ii)<= lmin) iloc=ii
       enddo
       if (iloc==-1) call abortmp( "sort() error")

       index(i)=iloc
       if(Debug) write(iulog,*)'sort: point #4'
       msk(iloc)=.FALSE.
       if(Debug) write(iulog,*)'sort: point #5'
       if(Debug) write(iulog,*)'sort: i, msk',i,msk
    enddo
    if(Debug) write(iulog,*)'sort: point #6'
    !DBG   write(iulog,*)'sort: list is:',list
    !DBG   write(iulog,*)'sort: index is:',index
    !DBG   stop

  end subroutine sort
end module zoltan_mod
