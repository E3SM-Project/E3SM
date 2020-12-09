#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module mesh_mod

 use kinds, only : real_kind, long_kind
  use physical_constants, only : DD_PI
  use control_mod, only : MAX_FILE_LEN

#ifndef HOMME_WITHOUT_PIOLIBRARY
  use netcdf ! _EXTERNAL
#endif

  implicit none
  logical, public           :: MeshUseMeshFile = .false.
#ifndef HOMME_WITHOUT_PIOLIBRARY
  public  :: MeshOpen           ! Must be called first


  integer, parameter :: MXSTLN = 32
#endif
  ! ===============================
  ! Public methods for mesh_mod
  ! ===============================

  public  :: MeshCubeEdgeCount  ! called anytime afer MeshOpen
  public  :: MeshCubeElemCount  ! called anytime afer MeshOpen
  public  :: MeshCubeTopology   ! called afer MeshOpen
  public  :: MeshCubeTopologyCoords   ! called afer MeshOpen

#ifndef HOMME_WITHOUT_PIOLIBRARY
  public  :: MeshSetCoordinates ! called after MeshCubeTopology
  public  :: MeshPrint          ! show the contents of the Mesh after it has been loaded into the module
  public  :: MeshClose
  ! ===============================
  ! Private members
  ! ===============================

  integer,private,parameter :: nfaces          = 6 ! number of faces on the cube
  integer,private,parameter :: nInnerElemEdge  = 8 ! number of edges for an interior element

  character (len=MAX_FILE_LEN), private              :: p_mesh_file_name
  integer                     , private              :: p_ncid
  integer                     , private              :: p_number_elements
  integer                     , private              :: p_number_elements_per_face
  integer                     , private              :: p_number_blocks
  integer                     , private              :: p_number_nodes
  integer                     , private              :: p_number_dimensions
  integer                     , private              :: p_number_neighbor_edges
   real(kind=real_kind)       , private, allocatable :: p_node_coordinates(:,:)
  integer                     , private, allocatable :: p_connectivity(:,:)


  ! Not used will eliminate later
  integer                     , private ::  p_elem_block_ids

  ! ===============================
  ! Private methods
  ! ===============================

  private :: create_index_table
  private :: find_side_neighbors
  private :: find_corner_neighbors
  private :: get_node_coordinates
  private :: get_2D_sub_coordinate_indexes
  private :: mesh_connectivity
  private :: cube_face_element_centroids
  private :: smallest_diameter_element
  private :: cube_to_cube_coordinates
  private :: sphere_to_cube_coordinates
  private :: initialize_space_filling_curve

  private :: handle_error
  private :: open_mesh_file
  private :: close_mesh_file
  private :: get_number_of_elements
  private :: get_number_of_dimensions
  private :: get_number_of_elements_per_face
  private :: get_number_of_nodes
  private :: get_number_of_element_blocks
  private :: get_node_multiplicity
  private :: get_face_connectivity

  ! Not used will eliminate later
  private :: get_block_ids

  ! Efficient SFC.
  private :: mesh_init_sfc_from_centroids_orig, mesh_init_sfc_from_centroids, mesh_theta2idx, gen_mesh_data
  public :: test_mesh_init_sfc_from_centroids

#endif

contains

#ifndef HOMME_WITHOUT_PIOLIBRARY
!======================================================================
!  subroutine handle_error
!======================================================================
  subroutine handle_error (status, file, line)
    use parallel_mod, only : abortmp
    implicit none
    integer,            intent(in) :: status
    character (len=*),  intent(in) :: file
    integer,            intent(in) :: line
    print *, file,':', line, ': ', trim(nf90_strerror(status))
    call abortmp("Terminating program due to netcdf error while obtaining mesh information, please see message in standard output.")
  end subroutine handle_error

!======================================================================
!  open_mesh_file()
!
!> Open the netcdf file containing the mesh.
!! Assign the holder to the file to p_ncid so everyone else knows
!! how to use it without passing the argument around.
!======================================================================
  subroutine open_mesh_file()
    implicit none
    integer                        :: status

    status = nf90_open(p_mesh_file_name, NF90_NOWRITE, p_ncid)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

    MeshUseMeshFile = .true.

  end subroutine open_mesh_file

!======================================================================
! subroutine close_mesh_file()
!======================================================================

  subroutine close_mesh_file()
    implicit none
    integer              :: status

    status = nf90_close(p_ncid)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

  end subroutine close_mesh_file

!======================================================================
! function get_number_of_dimensions()
!======================================================================

  function get_number_of_dimensions() result(number_dimensions)
    implicit none
    integer              :: number_dimensions

     ! local variables
     integer              :: status, number_of_dim_id

     ! Get the id of 'num_elem', if such dimension is not there panic and quit :P
    status = nf90_inq_dimid(p_ncid, "num_dim", number_of_dim_id)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

    ! How many values for 'num_elem' are there?
    status = nf90_inquire_dimension(p_ncid, number_of_dim_id, len = number_dimensions)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

  end function get_number_of_dimensions

!======================================================================
! function get_number_of_elements()
!======================================================================

  function get_number_of_elements() result(number_elements)
    implicit none
    integer              :: number_elements
    ! local variables
    integer              :: status, number_of_elements_id

    ! Get the id of 'num_elem', if such dimension is not there panic and quit :P
    status = nf90_inq_dimid(p_ncid, "num_elem", number_of_elements_id)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

    ! How many values for 'num_elem' are there?
    status = nf90_inquire_dimension(p_ncid, number_of_elements_id, len = number_elements)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

  end function get_number_of_elements

!======================================================================
!  function get_number_of_nodes()
!======================================================================
  function get_number_of_nodes() result(number_nodes)
    implicit none
    integer              :: number_nodes
    ! local variables
    integer              :: status, number_of_nodes_id

    ! Get the id of 'num_nodes', if such dimension is not there panic and quit :P
    status = nf90_inq_dimid(p_ncid, "num_nodes", number_of_nodes_id)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

    ! How many values for 'num_nodes' are there?
    status = nf90_inquire_dimension(p_ncid, number_of_nodes_id, len = number_nodes)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

  end function get_number_of_nodes


!======================================================================
!  function get_number_of_element_blocks()
!======================================================================
  function get_number_of_element_blocks() result(number_element_blocks)
    use parallel_mod, only : abortmp
    implicit none
    integer              :: number_element_blocks
    ! local variables
    integer              :: status, number_of_element_blocks_id

    ! Get the id of 'num_el_blk', if such dimension is not there panic and quit :P
    status = nf90_inq_dimid(p_ncid, "num_el_blk", number_of_element_blocks_id)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

    ! How many values for 'num_el_blk' are there?
    status = nf90_inquire_dimension(p_ncid, number_of_element_blocks_id, len = number_element_blocks)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)

     if (number_element_blocks /= 1) then
        if (number_element_blocks /= 6  ) then
            call abortmp('Reading cube-sphere from input file is not supported')
        else
           call abortmp('Number of elements blocks not exactly 1 (sphere) or 6 (cube)')
        endif
     endif

  end function get_number_of_element_blocks

!======================================================================
!  function get_number_of_elements_per_face()
!======================================================================
    function get_number_of_elements_per_face() result(number_elements_per_face)
    use parallel_mod, only : abortmp
    implicit none
    integer             :: number_elements_per_face

    integer               :: face_num ! For each of the face, we get the information
    character(len=MXSTLN) :: element_type ! Each face is composed of elements of certain type
    integer               :: number_elements_in_face ! How many elements in this face
    integer               :: num_nodes_per_elem ! How many nodes in each element
    integer               :: number_of_attributes ! How many attributes in the face

    integer               :: status, dimension_id

    if (p_number_blocks == 0)  then
       call abortmp('get_number_of_elements_per_face called before MeshOpen')
    else if (p_number_blocks == 1) then ! we are in the presence of a sphere
       ! First we get sure the number of nodes per element is four
       status = nf90_inq_dimid(p_ncid, "num_nod_per_el1", dimension_id)
       if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
       status = nf90_inquire_dimension(p_ncid, dimension_id, len =  num_nodes_per_elem)
       if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
       if (num_nodes_per_elem /= 4)  call abortmp('Number of nodes per element is not four')
       ! now we check how many elements there are in the face
       status = nf90_inq_dimid(p_ncid, "num_el_in_blk1", dimension_id)
       if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
       status = nf90_inquire_dimension(p_ncid, dimension_id, len = number_elements_in_face)
       if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
       number_elements_per_face =  number_elements_in_face
    else if (p_number_blocks == 6) then ! we are in the presence of a cube-sphere
       call abortmp('Reading a mesh for a cube-sphere is not supported')
    else
       call abortmp('Number of elements blocks not exactly 1 (sphere) or 6 (cube)')
    end if

  end function get_number_of_elements_per_face

!======================================================================
! This function is used to set the value of p_elem_block_ids  but such variable is never used
!======================================================================
  function get_block_ids(idexo) result(block_ids)
    use parallel_mod, only : abortmp
    implicit none
    integer(kind=long_kind), intent(in)  :: idexo
    integer(kind=long_kind)              :: block_ids(p_number_blocks)

    block_ids = 0

  end function get_block_ids



!======================================================================
! subroutine get_face_connectivity
!======================================================================
  subroutine get_face_connectivity()
    use parallel_mod, only : abortmp
    implicit none

    integer              :: var_id, status

    status = nf90_inq_varid(p_ncid, "connect1", var_id)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
    status = nf90_get_var(p_ncid, var_id, p_connectivity)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
  end subroutine get_face_connectivity

!======================================================================
! subroutine get_node_multiplicity
!======================================================================
  subroutine get_node_multiplicity(node_multiplicity)
    use parallel_mod, only : abortmp
    use dimensions_mod, only : max_elements_attached_to_node
    implicit none
    integer, intent(out) :: node_multiplicity(:)
    integer              :: node_num(4)

    integer              :: k, number_nodes

    node_multiplicity(:) = 0
    number_nodes = SIZE(node_multiplicity)
    ! check this external buffer was allocated correctly
    if (number_nodes /= p_number_nodes) call abortmp('Number of nodes does not matches size of node multiplicity array')
    ! for each node, we have for four other nodes

    if (minval(p_connectivity) < 1 .or. number_nodes < maxval(p_connectivity)) then
       call abortmp('get_node_multiplicity: Node number less than 1 or greater than max.')
    end if

    do k=1,p_number_elements_per_face
       node_num = p_connectivity(:,k)
       node_multiplicity(node_num) = node_multiplicity(node_num) + 1
    enddo

    if (minval(node_multiplicity) < 3 .or. max_elements_attached_to_node < maxval(node_multiplicity)) then
      print *, 'minval(node_multiplicity)', minval(node_multiplicity)
      print *, 'maxval(node_multiplicity)', maxval(node_multiplicity),&
           ' and max_elements_attached_to_node ',max_elements_attached_to_node
      call abortmp('get_node_multiplicity: Number of elements attached to node less than 3 or greater than maximum.')
    endif

  end subroutine get_node_multiplicity

!======================================================================
!  subroutine get_node_coordinates ()
!======================================================================
  subroutine get_node_coordinates ()
    use coordinate_systems_mod, only : cartesian3D_t
    use parallel_mod, only : abortmp

    implicit none
    integer              :: var_id, status

    status = nf90_inq_varid(p_ncid, "coord", var_id)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
    status = nf90_get_var(p_ncid, var_id, p_node_coordinates)
    if(status /= nf90_NoErr) call handle_error(status, __FILE__, __LINE__)
  end subroutine get_node_coordinates

  ! ================================================================================
  !
  ! -----------------Internal private routines that do not use netCDF IO -----------
  !
  ! ================================================================================

!======================================================================
! subroutine get_2D_sub_coordinate_indexes
!======================================================================
  subroutine get_2D_sub_coordinate_indexes(x, y, sgnx, sgny, face_no)
     implicit none
    integer, intent(in)              :: face_no
    integer, intent(out)             :: x,y
    integer, intent(out)             :: sgnx, sgny
    if (face_no == 1 .or. face_no == 3) then
       x = 2
       y = 3
    else if (face_no == 2 .or. face_no == 4) then
       x = 1
       y = 3
    else
       x = 2
       y = 1
    endif
    if (face_no == 1 .or. face_no == 4 .or. face_no == 5) then
       sgnx =  1
       sgny =  1
    else if (face_no == 2 .or. face_no == 3) then
       sgnx = -1
       sgny =  1
    else
       sgnx =  1
       sgny = -1
    endif
  end subroutine get_2D_sub_coordinate_indexes



!======================================================================
! subroutine mesh_connectivity(connect)
!
! puts the transpose of p_connectivity into connect
!======================================================================

  subroutine  mesh_connectivity (connect)
    use parallel_mod, only : abortmp
    implicit none
    integer,  intent(out) :: connect(p_number_elements,4)

    integer :: k, j

    if (0 == p_number_blocks)  call abortmp('mesh_connectivity called before MeshOpen')
    j=0
    do k=1, p_number_elements_per_face
       j=j+1
       connect(j,:) = p_connectivity(:,k)
    enddo

    if (j /= p_number_elements) call abortmp('mesh_connectivity: Number of elements in side sets not equal to total elements')

    if (minval(connect) < 1 .or. maxval(connect) > p_number_nodes) then
       call abortmp('mesh_connectivity: Node number out of bounds')
    end if

  end subroutine mesh_connectivity
!======================================================================
! subroutine create_index_table()
!
! this is needed to detremine side and corner neighbors
!======================================================================

  subroutine create_index_table(index_table, element_nodes)

    use dimensions_mod, only : max_elements_attached_to_node
    use parallel_mod, only : abortmp

    integer, allocatable, intent(inout)  :: index_table(:,:)
    integer             ,  intent(in)    :: element_nodes(p_number_elements, 4)
    integer                              :: cnt, cnt_index, node
    integer                              :: k, ll

    !Create an index table so that we can find neighbors on O(n)
    ! so for each node, we want to know which elements it is part of
    allocate(index_table(p_number_nodes, max_elements_attached_to_node + 1))

    !the last column in the index table is a count of the number of elements
    index_table = 0

    cnt_index =  max_elements_attached_to_node + 1

    do k=1,p_number_elements
        do ll=1,4
           node = element_nodes(k, ll) !the node
           cnt = index_table(node, cnt_index)  !how many elements for that node already in table
           cnt = cnt + 1 !increment since we are adding an element
           if (cnt >  max_elements_attached_to_node) then
              call abortmp('Found a node in too many elements.')
           endif
           index_table(node, cnt_index) = cnt
           index_table(node, cnt) = k !put the element in the indextable
        enddo
    enddo

  end subroutine create_index_table

!======================================================================
! subroutine find_side_neighbors()
!
! find the element neighbors to the n,s,e,w and put them in GridVertex_t
!  (only 1 neighbor to the n,s,e,w)
!======================================================================
  subroutine find_side_neighbors (GridVertex, normal_to_homme_ordering, element_nodes, edge_wgt, index_table)
    use coordinate_systems_mod, only : cartesian3D_t
    use gridgraph_mod, only   : GridVertex_t
    use parallel_mod, only : abortmp
    use dimensions_mod, only : max_elements_attached_to_node

    implicit none
    integer             ,  intent(in)    :: normal_to_homme_ordering(8)
    integer             ,  intent(in)    :: element_nodes(p_number_elements, 4)
    integer             ,  intent(in)    :: edge_wgt
    integer             ,  intent(in)     :: index_table(:,:)
    type (GridVertex_t) ,  intent(inout) :: GridVertex(:)

    integer                              :: i_node(2), my_node(2)
    integer                              :: neighbor, direction, init_size
    integer                              :: j,k,ll,i, m
    integer                              :: i_elem, jump, end_i
    integer                              :: loc, node, cnt_index, cnt, a_count(2)
    logical                              :: found
    if (0 == p_number_blocks)  call abortmp('find_side_neighbors called before MeshOpen')


    !the last column in the index table is a count of the number of elements
    cnt_index =  max_elements_attached_to_node + 1

    !use index table to find neighbors
    do k=1,p_number_elements  ! for each element k
       !set the side weights
       GridVertex(k)%nbrs_wgt(1:4) = edge_wgt
       do ll=1,4 !loop through the four sides

          jump = normal_to_homme_ordering(ll)
          loc = GridVertex(k)%nbrs_ptr(jump)

          if (GridVertex(k)%nbrs(loc) == 0) then  !if side is not set yet, then
             !look for side element
             found = .false.
             neighbor = 0

             my_node(1) = element_nodes(k, ll)
             a_count(1) = index_table(my_node(1), cnt_index)
             my_node(2) = element_nodes(k, mod(ll,4)+1)
             a_count(2) = index_table(my_node(2), cnt_index)

             !loop through the elements that are in the index table for each node
             !and find the element number and direction of the side neighbor
             do m = 1,2
                if (found) exit
                end_i = a_count(m)
                do i = 1, end_i
                   if (found) exit
                   i_elem = index_table(my_node(m),i)
                   if (i_elem /= k) then !k is the element we are setting sides for
                      do j=1,4 !loop through each of i_elem's four sides
                         i_node(1) = element_nodes(i_elem, j)
                         i_node(2) = element_nodes(i_elem, mod(j,4)+1)
                         if ( (i_node(1) == my_node(2) .and. i_node(2) == my_node(1)) .or. &
                              (i_node(1) == my_node(1) .and. i_node(2) == my_node(2)) ) then
                            neighbor = i_elem
                            direction = j
                            found = .true.
                            !found a match
                            exit
                         end if
                      end do ! j loop
                   end if
                enddo ! i loop
             enddo !m loop

             if (neighbor == 0) call abortmp('find_side_neighbor: Neighbor not found! Every side should have a neighbor.')

             GridVertex(k)%nbrs(loc) = neighbor
             jump = normal_to_homme_ordering(direction)
             loc = GridVertex(neighbor)%nbrs_ptr(jump)
             GridVertex(neighbor)%nbrs(loc)= k
          endif
       enddo !  ll loop => 4 sides
    enddo ! k loop: each element

    do k=1,p_number_elements
       do ll=1,4
          if ( 0 == GridVertex(k)%nbrs(ll)) then
             call abortmp('Found one side of one element witout a neighbor.  Bummer!')
          end if
       end do
    end do

  end subroutine find_side_neighbors

!======================================================================
! function smallest_diameter_element
!======================================================================

  function smallest_diameter_element(element_nodes) result(min_diameter)
    use parallel_mod, only   : abortmp

    implicit none
    integer             ,intent(in)  :: element_nodes(:,:)

    integer                          :: i, j
    integer                          :: node_numbers(4)
    real(kind=real_kind)             :: coordinates (4,3)
    real                             :: x(3), y(3), r(3), d, min_diameter

    if (SIZE(element_nodes,dim=1) /= p_number_elements) then
       call abortmp('smallest_diameter_element:Element count check failed in &
            &exodus_mesh. Connectivity array length not equal to number of elements.')
    end if
    if ( p_number_elements_per_face /= p_number_elements) then
       call abortmp('smallest_diameter_element: Element count check failed in &
            &exodus_mesh. Element array length not equal to sum of face.')
    end if

    min_diameter = 9999999.
    do i=1, p_number_elements
       node_numbers = element_nodes(i,:)
       coordinates = p_node_coordinates(node_numbers,:)
       ! smallest side length
       do j=1,4
          x = coordinates(j         ,:)
          y = coordinates(1+MOD(j,4),:)
          r = x-y
          d  = dot_product(r,r)
          if (d < min_diameter ) then
             min_diameter = d
          end if
       end do
       ! smallest diameter length
       do j=1,2
          x = coordinates(j         ,:)
          y = coordinates(2+MOD(j,4),:)
          r = x-y
          d  = dot_product(r,r)
          if (d < min_diameter ) then
             min_diameter = d
          end if
       end do
    enddo
    min_diameter = SQRT(min_diameter)
  end function smallest_diameter_element

!======================================================================
!  subroutine cube_to_cube_coordinates
!======================================================================

  subroutine cube_to_cube_coordinates (cube_coor, node_coor, face_number)
    use parallel_mod,    only   : abortmp
    use physical_constants, only : dd_pi
    implicit none
    real(kind=real_kind),    intent(in)  :: node_coor(4,3)
    integer,                 intent(in)  :: face_number
    real(kind=real_kind),    intent(out) :: cube_coor(4,2)
    real(kind=real_kind)                 :: test_coor(4,2)

    integer                              :: i,j,x_index, y_index, sgnx, sgny
    call get_2D_sub_coordinate_indexes(x_index, y_index, sgnx, sgny, face_number)
    cube_coor(:,1) = sgnx*node_coor(:,x_index)
    cube_coor(:,2) = sgny*node_coor(:,y_index)
  end subroutine cube_to_cube_coordinates


!======================================================================
!  subroutine sphere_to_cube_coordinates
!======================================================================

  subroutine sphere_to_cube_coordinates (cube_coor, node_coor, face_number)
    use coordinate_systems_mod, only   : cartesian3D_t, cartesian2d_t, spherical_polar_t, &
         change_coordinates, sphere2cubedsphere
    implicit none
    real(kind=real_kind),    intent(in)   :: node_coor(4,3)
    integer,                 intent(in)   :: face_number
    real(kind=real_kind),    intent(out)  :: cube_coor(4,2)
    integer                               :: i, l
    type(cartesian2d_t)                   :: cart(4)

    do i=1,4
       cart(i) = sphere2cubedsphere(change_coordinates(node_coor(i,:)), face_number)
    end do
    cube_coor(:,1) = cart(:)%x
    cube_coor(:,2) = cart(:)%y
  end subroutine sphere_to_cube_coordinates


!======================================================================
!  subroutine cube_face_element_centroids
!======================================================================

  subroutine cube_face_element_centroids(centroids, face_numbers, element_nodes)
    use parallel_mod,    only   : abortmp
    implicit none
    integer            , intent(in)  :: element_nodes(:,:)
    integer,             intent(in)  :: face_numbers    (p_number_elements)
    real(kind=real_kind),intent(out) :: centroids       (p_number_elements,2)
    real(kind=real_kind)             :: coordinates(4,3)
    real(kind=real_kind)             :: cube_coor  (4,2)
    integer                          :: i, j, node_numbers(4)

    if (0 == p_number_blocks)  call abortmp('cube_face_element_centroids called before MeshOpen')
    if (SIZE(element_nodes,dim=1) /= p_number_elements) then
       call abortmp('cube_face_element_centroids:Element count check failed in &
            &exodus_mesh. Connectivity array length not equal to number of elements.')
    end if
    if ( p_number_elements_per_face /= p_number_elements ) then
       call abortmp('cube_face_element_centroids: Element count check failed in &
            &exodus_mesh. Element array length not equal to sum of face.')
    end if

    do i=1, p_number_elements
       node_numbers = element_nodes(i,:)
       coordinates = p_node_coordinates(node_numbers,:)
       if (6 == p_number_blocks) then
          call cube_to_cube_coordinates   (cube_coor, coordinates, face_numbers(i))
       else
          call sphere_to_cube_coordinates (cube_coor, coordinates, face_numbers(i))
       end if
       centroids(i,:) = SUM(cube_coor,dim=1)/4.0
    enddo
  end subroutine cube_face_element_centroids

!======================================================================
! subroutine initialize_space_filling_curve
!======================================================================
  subroutine initialize_space_filling_curve(GridVertex, element_nodes)
    use gridgraph_mod, only   : GridVertex_t
    use parallel_mod,  only   : abortmp

    implicit none
    type (GridVertex_t), intent(inout) :: GridVertex(:)
    integer            , intent(in)    :: element_nodes(:,:)

    real(kind=real_kind)               :: centroids(p_number_elements,2)
    integer                            :: face_numbers(p_number_elements)
    real(kind=real_kind)               :: h
    integer                            :: i, sfcidx(p_number_elements)

    if (SIZE(GridVertex) /= p_number_elements) then
       call abortmp('initialize_space_filling_curve:Element count check failed &
            &in exodus_mesh. Vertex array length not equal to number of elements.')
    end if
    if (SIZE(element_nodes,dim=1) /= p_number_elements) then
       call abortmp('initialize_space_filling_curve:Element count check failed &
            &in exodus_mesh. Connectivity array length not equal to number of elements.')
    end if

    face_numbers(:) = GridVertex(:)%face_number
    h = smallest_diameter_element    (                         element_nodes)

    call cube_face_element_centroids (centroids, face_numbers, element_nodes)

    if (h<.00001) call abortmp('initialize_space_filling_curve: Unreasonably small element found. less than .00001')

    call mesh_init_sfc_from_centroids(p_number_elements, face_numbers, centroids, h, sfcidx)

    do i = 1, p_number_elements
       GridVertex(i)%SpaceCurve = sfcidx(i)
       if (-1==GridVertex(i)%SpaceCurve) then
          write (*,*) " Error in projecting element ",i," to space filling curve."
          write (*,*) " Face:",face_numbers(i)
          write (*,*) " Centroid:",centroids(i,:)
       end if
    end do
  end subroutine initialize_space_filling_curve

  subroutine mesh_init_sfc_from_centroids_orig(nelem, face_numbers, centroids, h, sfcidx)
    use spacecurve_mod, only  : GenspaceCurve

    integer, intent(in) :: nelem, face_numbers(nelem)
    real(kind=real_kind), intent(in) :: centroids(nelem,2), h
    integer, intent(out) :: sfcidx(nelem)
    integer, allocatable :: Mesh2(:,:), Mesh2_map(:,:), sfcij(:,:)
    integer :: ne, ne2, sfc_index, i, j, i2, j2, face
    real(kind=real_kind) :: x, y
    real(kind=real_kind), parameter :: pi = 3.141592653589793116d0
    integer, parameter :: nfaces = 6

    ne = ceiling(pi/h);
    ! find the smallest ne2 which is a power of 2 and ne2>ne
    ne2=2**ceiling(log(real(ne))/log(2d0))

    allocate(Mesh2(ne2,ne2))
    allocate(Mesh2_map(ne2,ne2))
    allocate(sfcij(0:ne2*ne2,2))
    ! create a reverse index array for Mesh2
    ! j = Mesh2(i,j)
    ! (i,j) = (sfcij(j,1),sfci(j,2))
    call GenspaceCurve(Mesh2)  ! SFC partition for ne2
    do j2=1,ne2
       do i2=1,ne2
          j=Mesh2(i2,j2)
          sfcij(j,1)=i2
          sfcij(j,2)=j2
       enddo
    enddo
    sfc_index = 0
    do face = 1,nfaces
       ! associate every element on the ne x ne mesh (Mesh) with its closest
       ! element on the ne2 x ne2 mesh (Mesh2) Store this as a map from Mesh2 ->
       ! Mesh in Mesh2_map.  elements in Mesh2 which are not mapped get assigned
       ! a value of 0
       Mesh2_map = 0
       do i=1,nelem
          if (face_numbers(i) == face ) then
             x = centroids(i,1)
             y = centroids(i,2)
             ! map this element to an (i2,j2) element
             ! [ -pi/4, pi/4 ]  -> [ 0, ne2 ]
             i2 = nint( (0.5 + 2.0*x/pi)*ne2 + .5 )
             j2 = nint( (0.5 + 2.0*y/pi)*ne2 + .5 )
             if (face == 4 .or. face == 6 )               i2 = ne2-i2+1
             if (face == 1 .or. face == 2 .or. face == 6) j2 = ne2-j2+1
             if (i2<1  ) i2=1
             if (i2>ne2) i2=ne2
             if (j2<1  ) j2=1
             if (j2>ne2) j2=ne2
             Mesh2_map(i2,j2)=i
          end if
       end do
       ! generate a SFC for Mesh with the same ordering as the elements in Mesh2
       ! which map to Mesh.
       do j = 0, ne2*ne2-1
          i2 = sfcij(j,1)
          j2 = sfcij(j,2)
          i = Mesh2_map(i2,j2)
          if (i /= 0) then
             ! (i2,j2) element maps to element
             sfcidx(i) = sfc_index
             sfc_index = sfc_index+1
          endif
       enddo
    enddo
    deallocate(Mesh2)
    deallocate(Mesh2_map)
    deallocate(sfcij)
  end subroutine mesh_init_sfc_from_centroids_orig

  function mesh_theta2idx(x, ne) result(i)
    real(kind=real_kind), intent(in) :: x
    integer, intent(in) :: ne
    real(kind=real_kind), parameter :: pi = 3.141592653589793116d0, half = 0.5d0
    integer :: i

    i = nint((half + 2*x/pi)*ne + half)
    if (i < 1) i = 1
    if (i > ne) i = ne
  end function mesh_theta2idx

  subroutine mesh_init_sfc_from_centroids(nelem, face_numbers, centroids, h, sfcidx)
    use spacecurve_mod, only : sfcmap_t, sfcmap_init, sfcmap_pos2i, sfcmap_finalize
    use sort_mod, only : sortints

    integer, intent(in) :: nelem, face_numbers(nelem)
    real(kind=real_kind), intent(in) :: centroids(nelem,2), h
    integer, intent(out) :: sfcidx(nelem)
    real(kind=real_kind), parameter :: pi = 3.141592653589793116d0
    integer, parameter :: nfaces = 6
    type (sfcmap_t) :: sfcmap
    integer :: ne, ne2, i, face, i2, j2, pos(2), o, sfc, ne2sq
    real(kind=real_kind) :: x, y
    integer, allocatable :: inds(:,:)

    ne = ceiling(pi/h);
    ! find the smallest ne2 which is a power of 2 and ne2>ne. could instead use
    ! find_next_factorable for a smaller ne2.
    ne2 = 2**ceiling(log(real(ne))/log(2d0))
    ne2sq = ne2*ne2

    o = sfcmap_init(ne2, sfcmap)
    allocate(inds(2,nelem))
    do i = 1, nelem
       face = face_numbers(i)
       x = centroids(i,1)
       y = centroids(i,2)
       i2 = mesh_theta2idx(x, ne2)
       j2 = mesh_theta2idx(y, ne2)
       if (face == 4 .or. face == 6)                i2 = ne2-i2+1
       if (face == 1 .or. face == 2 .or. face == 6) j2 = ne2-j2+1
       pos = (/ i2, j2 /)
       inds(2,i) = i
       o = sfcmap_pos2i(sfcmap, pos, sfc)
       inds(1,i) = ne2*ne2*(face-1) + sfc
    end do
    call sortints(inds)
    do i = 1, nelem
       sfcidx(inds(2,i)) = i-1
    end do
    deallocate(inds)
    call sfcmap_finalize(sfcmap)
  end subroutine mesh_init_sfc_from_centroids

  subroutine gen_mesh_data(nx, ny, centroids, h)
    integer, intent(in) :: nx, ny
    real(kind=real_kind), intent(out) :: centroids(ny,nx,2), h
    integer :: iy, ix, d, i, j, iy1, ix1
    real(kind=real_kind) :: dx, dy, xo, yo, x, y, r(2), pmin(2), pmax(2), x0, y0, x1, y1
    real(kind=real_kind), parameter :: pi = 3.141592653589793116d0, delta = 0.01d0

    dx = 1.0d0 / (nx-1)
    dy = 1.0d0 / (ny-1)

    pmin = (/ 1, 1 /)
    pmax = (/ 0, 0 /)

    ! Generate an interesting mesh.
    do ix = 1, nx
       do iy = 1, ny
          call random_number(r)
          r = 2*(r - 0.5d0)
          xo = (ix-1)*dx
          yo = (iy-1)*dy
          x = xo + 0.2*dx*r(1) + 0.1d0*cos(pi*yo)
          y = yo + 0.2*dy*r(2) + 0.1d0*sin(pi*xo)
          centroids(iy,ix,:) = (/ x, y /)
          do d = 1, 2
             pmin(d) = min(pmin(d), centroids(iy,ix,d))
             pmax(d) = max(pmax(d), centroids(iy,ix,d))
          end do
       end do
    end do
    ! Normalize it to fit in (0,1)x(0,1)
    do ix = 1, nx
       do iy = 1, ny
          do d = 1, 2
             centroids(iy,ix,d) = delta + (centroids(iy,ix,d) - pmin(d)) / &
                  (2*delta + (pmax(d) - pmin(d)))
          end do
       end do
    end do
    ! Get the smallest element dimension.
    h = 1
    do ix = 1, nx
       do iy = 1, ny
          x0 = centroids(iy,ix,1)
          y0 = centroids(iy,ix,2)
          do i = -1, 1
             do j = -1, 1
                if (i == 0 .and. j == 0) cycle
                ix1 = modulo(ix + i - 1, nx) + 1
                iy1 = modulo(iy + j - 1, ny) + 1
                x1 = centroids(iy1,ix1,1)
                y1 = centroids(iy1,ix1,2)
                h = min(h, sqrt((x1 - x0)**2 + (y1 - y0)**2))
             end do
          end do
       end do
    end do
    ! Scale to [-pi/4, pi/4]
    do ix = 1, nx
       do iy = 1, ny
          do d = 1, 2
             centroids(iy,ix,d) = 0.5d0*pi*(centroids(iy,ix,d) - 0.5d0)
          end do
       end do
    end do
  end subroutine gen_mesh_data

  subroutine test_mesh_init_sfc_from_centroids(nx, ny)
    integer, intent(in) :: nx, ny
    real(kind=real_kind), allocatable :: centroids(:,:)
    integer, allocatable :: face_numbers(:), sfcidxo(:), sfcidxn(:)
    real(kind=real_kind) :: start, finish
    real(kind=real_kind) :: h
    integer :: i, nelem

    nelem = nx*ny
    allocate(centroids(nelem,2), face_numbers(nelem), &
         sfcidxo(nelem), sfcidxn(nelem))
    do i = 1, nelem
       face_numbers(i) = modulo(i, 6) + 1
    end do
    call gen_mesh_data(nx, ny, centroids, h)
    !call mesh_print(nx, ny, centroids, h, 'centroids.hy')
    call cpu_time(start)
    call mesh_init_sfc_from_centroids_orig(nelem, face_numbers, centroids, h, sfcidxo)
    call cpu_time(finish)
    print *, 'orig',finish-start
    call cpu_time(start)
    call mesh_init_sfc_from_centroids(nelem, face_numbers, centroids, h, sfcidxn)
    call cpu_time(finish)
    print *, 'new ',finish-start
    do i = 1, nelem
       if (sfcidxo(i) /= sfcidxn(i)) then
          print '(a,i6,i6,i6,i3,i6,i6)','ERROR',i,nx,ny,face_numbers(i),sfcidxo(i),sfcidxn(i)
          exit
       end if
    end do
    deallocate(centroids, face_numbers, sfcidxo, sfcidxn)
  end subroutine test_mesh_init_sfc_from_centroids

!======================================================================
! subroutine find_corner_neighbors
!======================================================================

  subroutine find_corner_neighbors  (GridVertex, normal_to_homme_ordering, element_nodes, corner_wgt, index_table)
    use parallel_mod,           only : abortmp
    use gridgraph_mod,          only : GridVertex_t
    use dimensions_mod,         only : max_elements_attached_to_node, max_corner_elem
    use control_mod, only: north, south, east, west, neast,seast, nwest,swest
    implicit none

    type (GridVertex_t), intent(inout) :: GridVertex(:)
    integer            , intent(in)    :: normal_to_homme_ordering(8)
    integer            , intent(in)    :: element_nodes(p_number_elements, 4)
    integer            , intent(in)    :: corner_wgt
    integer            , intent(in)    :: index_table(:,:)

    integer                          :: node_elements (2*max_elements_attached_to_node)
    integer                          :: elem_neighbor (4*max_elements_attached_to_node)
    integer                          :: nbr_cnt(4)
    integer                          :: elem_nbr_start, start
    integer                          :: i, j, k, ll, jj, kk
    integer                          :: node, loc, cnt, cnt_index
    integer                          :: corner_array(max_corner_elem), orig_pos(max_corner_elem)
    integer                          :: face_array(max_corner_elem), a_corner_elems(max_corner_elem)
    integer                          :: corner_sides(2)
    integer                          :: side_elem, corner_elem, tmp_s

    !the last column in the index table is a count of the number of elements
    cnt_index =  max_elements_attached_to_node + 1

    do i=1, p_number_elements  !loop through all elements
       node_elements(:) = 0
       elem_neighbor(:) = 0
       elem_nbr_start = 0
       nbr_cnt(:) = 0

       do j=1,4  !check each of the 4 nodes at the element corners
          node = element_nodes(i,j)
          cnt = index_table(node, cnt_index)
          if (cnt < 3 .or. max_elements_attached_to_node < cnt) then
             call abortmp('find_corner_neighbors: Number of elements attached to node less than 3 or greater than maximum.')
          endif

          node_elements(1:cnt) = index_table(node, 1:cnt)

          !now node_elements contains the element neighbors to that node - so grab the
          ! corner neighbors - these are the ones that are not already side neighbors (or myself)
          k = 0
          do ll=1,cnt
             if ( i /= node_elements(ll) .and. & !not me
                  GridVertex(i)%nbrs(1) /= node_elements(ll) .and. & !not side 1
                  GridVertex(i)%nbrs(2) /= node_elements(ll) .and. & ! etc ...
                  GridVertex(i)%nbrs(3) /= node_elements(ll) .and. &
                  GridVertex(i)%nbrs(4) /= node_elements(ll)) then
                k = k + 1
                elem_neighbor(elem_nbr_start + k) = node_elements(ll)
             end if
          end do ! end of ll loop for multiplicity

          !keep track of where we are starting in elem_neighbor for each corner j
          elem_nbr_start = elem_nbr_start + k
          nbr_cnt(j) = k !how many neighbors in this corner
       end do ! end of j loop through 4 nodes


       ! now that we have done the 4 corners we can populate nbrs and nbrs_ptr
       ! with the corners in the proper order (clockwise) in neighbors
       ! also we can add the corner weight

       do j=5,8 !loop through 4 corners
          elem_nbr_start = 1
          !easiest to do the corner in ascending order - find loc
          do jj = 5,8
             ll = normal_to_homme_ordering(jj)
             if (j == ll) then
                loc = jj
                exit
             end if
             elem_nbr_start = elem_nbr_start + nbr_cnt(jj-4)
          end do

          start =  GridVertex(i)%nbrs_ptr(j)
          cnt = nbr_cnt(loc - 4)
          GridVertex(i)%nbrs_ptr(j+1) = start + cnt

          if (cnt > 0) then
             GridVertex(i)%nbrs(start : start + cnt-1) = &
                  elem_neighbor(elem_nbr_start : elem_nbr_start + cnt -1)
             GridVertex(i)%nbrs_face(start : start + cnt - 1) = &
                  GridVertex(elem_neighbor(elem_nbr_start : elem_nbr_start + cnt -1))%face_number
             GridVertex(i)%nbrs_wgt(start : start + cnt-1)  = corner_wgt

          end if

          ! within each corner neighbor, lets list the corners in clockwise order
          if (cnt > 1) then !cnt is the number of neighbors in this corner j
                            !there can be at most max_corner element of these

             a_corner_elems = 0
             a_corner_elems(1:cnt) = elem_neighbor(elem_nbr_start : elem_nbr_start + cnt -1)
             !corner-sides(2) is clockwise of corner_side(1)
             corner_array= 0
             orig_pos = 0
             select case (j)
             case(neast)
                corner_sides(1) = north
                corner_sides(2) = east
             case(seast)
                corner_sides(1) = east
                corner_sides(2) = south
             case(swest)
                corner_sides(1) = south
                corner_sides(2) = west
             case(nwest)
                corner_sides(1) = west
                corner_sides(2) = north
             end select

             !so the first element to list touches  corner_sides(1) element
             side_elem = GridVertex(i)%nbrs(corner_sides(1))

             !loop though the corner elements and see if any have a side neighbor
             !that = side_elem
             do k = 1,cnt !number of corner elements
                corner_elem = a_corner_elems(k)
                do kk = 1,4 !number of sides to check
                   loc = GridVertex(corner_elem)%nbrs_ptr(kk)
                   tmp_s = GridVertex(corner_elem)%nbrs(loc)
                   if (tmp_s == side_elem) then
                      corner_array(1) = corner_elem
                      orig_pos(1) = k
                      exit
                   endif
                enddo
                if  (corner_array(1)> 0) exit
             enddo
             if (corner_array(1)==0) then
                print *, i, cnt
                call abortmp('find_corner_neighbors (1) : mistake finding corner neighbor order')
             endif

             !if cnt == 2, we are done (we know the order of neighbors)
             if (cnt ==2) then
                if (corner_array(1) ==  a_corner_elems(1)) then
                   corner_array(2) =  a_corner_elems(2)
                   orig_pos(2) = 2
                else
                   corner_array(2) =  a_corner_elems(1)
                   orig_pos(2) = 1
                end if
             else !cnt = 3 or 4
                !find which corner element borders corner_sides(2)
                side_elem = GridVertex(i)%nbrs(corner_sides(2))
                do k = 1,cnt
                   corner_elem = a_corner_elems(k)
                   do kk = 1,4
                      loc = GridVertex(corner_elem)%nbrs_ptr(kk)
                      tmp_s = GridVertex(corner_elem)%nbrs(loc)
                      if (tmp_s == side_elem) then
                         corner_array(4) = corner_elem
                         orig_pos(4) = k
                         exit
                      endif
                   enddo
                   if  (corner_array(4)> 0) exit
                enddo
                if (corner_array(4)==0 .or. corner_array(4) == corner_array(1)) then
                   print *, i, cnt
                   call abortmp('find_corner_neighbors (2) : mistake finding corner neighbor order')
                endif

                !now if cnt = 3 then we are done
                if (cnt ==3) then
                   corner_array(3) = corner_array(4)
                   orig_pos(3) = orig_pos(4)

                   do k = 1,cnt !find the "middle" element
                      if (k /= orig_pos(1) .and. k /= orig_pos(3)) then
                         orig_pos(2) = k
                         corner_array(2) = a_corner_elems(k)
                         exit
                      endif
                   enddo
                else  !cnt = 4
                   !which of the two unassigned elements borders the element in
                   !corner_array(1) => put in corner_array(2)
                   side_elem = corner_array(1)

                   do k = 1,cnt
                      corner_elem = a_corner_elems(k)
                      if (corner_elem == corner_array(4) .or. corner_elem == corner_array(1)) then
                         cycle
                      else
                         do kk = 1,4 !check each side
                            loc = GridVertex(corner_elem)%nbrs_ptr(kk)
                            tmp_s = GridVertex(corner_elem)%nbrs(loc)
                            if (tmp_s == side_elem) then
                               corner_array(2) = corner_elem
                               orig_pos(2) = k
                               exit
                            endif
                         enddo
                      endif
                      if  (corner_array(2)> 0) exit
                   enddo
                   !now put the remaining one in pos 3
                   do k = 1,cnt
                      corner_elem = a_corner_elems(k)
                      if (corner_elem /= corner_array(4) .and. corner_elem /= &
                           corner_array(2) .and. corner_elem /= corner_array(1)) then
                         corner_array(3) = corner_elem
                         orig_pos(3) = k
                         exit
                      endif
                   enddo
                endif ! end of cnt=4
             endif! end of not cnt=2

             !now re-set the elements in this corner
             GridVertex(i)%nbrs(start : start + cnt-1) = corner_array(1:cnt)
             !nbrs_wgt are the same - nothing to do
             !fix neighbors face
             do k = 1,cnt
                face_array(k) = GridVertex(i)%nbrs_face(start + orig_pos(k) - 1)
             end do
             GridVertex(i)%nbrs_face(start : start + cnt - 1) = face_array(1:cnt)
          endif !end of cnt > 1 loop for corners

       end do !j loop through each corner

    end do ! end of i loop through elements
  end subroutine find_corner_neighbors

  ! ================================================================================
  !
  ! -------------------------------Public Methods-----------------------------------
  !
  ! ================================================================================

!======================================================================
! subroutine MeshOpen
!======================================================================

  subroutine MeshOpen(mesh_file_name, par)
    use parallel_mod, only : abortmp, parallel_t
    use kinds, only : real_kind, iulog

    implicit none
    character (len=*), intent(in) :: mesh_file_name
    type (parallel_t), intent(in) :: par

    integer              :: ncid
    integer, allocatable :: node_multiplicity(:)
    integer              :: k

    p_mesh_file_name    = mesh_file_name
    call open_mesh_file ()

    p_number_elements   = get_number_of_elements       ()
    p_number_nodes      = get_number_of_nodes          ()
    p_number_blocks     = get_number_of_element_blocks ()
    p_number_dimensions = get_number_of_dimensions     ()

    if (p_number_dimensions /= 3) then
       call abortmp('The number of dimensions must be 3, otherwise the mesh algorithms will not work')
    endif

    ! Only spheres are allowed in input files.
    if (par%masterproc) then
       if (p_number_blocks == 1) then
          write(iulog,*) "Since the mesh file has only one block, it is assumed to be a sphere."
       endif
    end if

    if (p_number_blocks /= 1) then
       call abortmp('Number of elements blocks not exactly 1 (sphere)')
    end if

    p_number_elements_per_face = get_number_of_elements_per_face()
    ! Because all elements are in one face, this value must match  p_number_elements
    if ( p_number_elements /= p_number_elements_per_face) then
       call abortmp('The value of the total number of elements does not match all the elements found in face 1')
    end if

    allocate( p_connectivity(4,p_number_elements_per_face) )
    p_connectivity(:,:)=0
    ! extract the connectivity from the netcdf file
    call get_face_connectivity()

    allocate(node_multiplicity(p_number_nodes))
    call get_node_multiplicity(node_multiplicity)

    ! tricky:  For each node with multiplicity n, there are n(n-1) neighbor links
    ! created.  But this counts each edge twice, so:  n(n-1) -n
    ! Should be the same as SUM(SIZE(GridVertex(i)%nbrs(j)%n),i=1:p_number_elements,j=1:8)
    ! p_number_neighbor_edges = dot_product(mult,mult) - 2*sum(mult)
    p_number_neighbor_edges = 0
    do k=1,p_number_nodes
      p_number_neighbor_edges = p_number_neighbor_edges + node_multiplicity(k)*(node_multiplicity(k)-2)
    end do

    deallocate(node_multiplicity)

    ! allocate the space for the coordinates, this is used in many functions
    allocate(p_node_coordinates(p_number_nodes, p_number_dimensions))
    call get_node_coordinates()

    if (p_number_elements_per_face /= p_number_elements) then
       call abortmp('MeshOpen: Total number of elements not equal to the number of elements on face 1!')
    end if

  end subroutine MeshOpen

!======================================================================
! subroutine MeshClose
!
! This routine acts as a destructor cleaning the memory allocated in MeshOpen
! which acts as a constructor allocated dynamical memory for the nodes coordinates.
!======================================================================

  subroutine MeshClose

    ! release memory
    deallocate(p_node_coordinates)
    deallocate(p_connectivity)
    ! let the file go
    call close_mesh_file ()

  end subroutine MeshClose


!======================================================================
! subroutine MeshPrint
!======================================================================


  subroutine MeshPrint(par)
    use parallel_mod, only : abortmp, parallel_t
    implicit none
    type (parallel_t), intent(in) :: par
    if (par%masterproc) then
       print *, 'This are the values for file ', trim(p_mesh_file_name)
       print *, 'The value for the number of dimensions (num_dim) is ', p_number_dimensions
       print *, 'The number of elements in the mesh file is ', p_number_elements
       print *, 'The number of nodes in the mesh file is ', p_number_nodes
       print *, 'The number of blocks in the mesh file is ',  p_number_blocks
       print *, 'The number of elements in the face 1 (sphere) is ',  p_number_elements_per_face
       if ( p_number_elements == p_number_elements) then
          print *, 'The value of the total number of elements does match all the elements found in face 1 (the only face)'
       else
          print *, 'The value of the total number of elements does not match all the elements found in face 1'
          print *, 'This message should not be appearing, there is something wrong in the code'
       endif
       print *, 'The number of neighbor edges ', p_number_neighbor_edges
       !print *, 'The node connectivity are (compare with ncdump -v connect1) ', p_connectivity
       !print *, ' ========================================================='
       !print *, 'The node coordinates are (compare with ncdump -v coord) ', p_node_coordinates
    endif

  end subroutine MeshPrint

!======================================================================
! subroutine MeshCubeTopologyCoords
!======================================================================
#endif

   subroutine MeshCubeTopologyCoords(GridEdge, GridVertex, coord_dim1, coord_dim2, coord_dim3, coord_dimension)
    use parallel_mod,           only : abortmp
    use dimensions_mod,         only : np,  max_elements_attached_to_node
    use coordinate_systems_mod, only : cartesian2D_t, cartesian3D_t, cube_face_number_from_cart, &
                                       cube_face_number_from_sphere, cart2cubedsphere
    use gridgraph_mod,          only : GridVertex_t
    use gridgraph_mod,          only : GridEdge_t
    use cube_mod,               only : CubeSetupEdgeIndex
    use gridgraph_mod,          only : initgridedge, num_neighbors
    use control_mod,            only : north, south, east, west, neast, seast, swest, nwest, partmethod, &
                                       z2_map_method, coord_transform_method
    use params_mod,             only : SFCURVE, SPHERE_COORDS, CUBE_COORDS, FACE_2D_LB_COORDS

    use kinds, only : iulog, real_kind
    use zoltan_mod, only : is_zoltan_partition, is_zoltan_task_mapping

    implicit none
    type (GridEdge_t),   intent(inout) :: GridEdge(:)
    type (GridVertex_t), intent(inout) :: GridVertex(:)
    real(kind=real_kind),allocatable, intent(inout) :: coord_dim1(:)
    real(kind=real_kind),allocatable, intent(inout) :: coord_dim2(:)
    real(kind=real_kind),allocatable, intent(inout) :: coord_dim3(:)
    integer, intent(inout) :: coord_dimension

#ifndef HOMME_WITHOUT_PIOLIBRARY
    real (kind=real_kind)            :: x,y,z, rangex
    real(kind=real_kind)             :: coordinates(4,3)
    real(kind=real_kind)             :: centroid(3)
    type (cartesian3D_t)             :: face_center
    type (cartesian2D_t)             :: face_center_2d

    integer                          :: i, j, k, ll, m, loc, face_no
    integer                          :: element_nodes(p_number_elements, 4)
    integer                          :: EdgeWgtP,CornerWgt
    integer                          :: normal_to_homme_ordering(8)
    integer                          :: node_numbers(4)
    integer, allocatable             :: index_table(:,:)


    normal_to_homme_ordering(1) = south
    normal_to_homme_ordering(2) =  east
    normal_to_homme_ordering(3) = north
    normal_to_homme_ordering(4) =  west
    normal_to_homme_ordering(5) = swest
    normal_to_homme_ordering(6) = seast
    normal_to_homme_ordering(7) = neast
    normal_to_homme_ordering(8) = nwest

    coord_dimension = 3


    if (SIZE(GridVertex) /= p_number_elements) then
       call abortmp('MeshCubeTopology: Element count check failed in exodus_mesh. &
            &Vertex array length not equal to number of elements.')
    end if
    if (p_number_elements_per_face /= p_number_elements) then
       call abortmp('MeshCubeTopology: Element count check failed in exodus_mesh. &
            &Element array length not equal to sum of face.')
    end if

    EdgeWgtP = np
    CornerWgt = 1


    call mesh_connectivity (element_nodes)


    do i=1, p_number_elements
       GridVertex(i)%number           = i
       GridVertex(i)%face_number      = 0
       GridVertex(i)%processor_number = 0
       GridVertex(i)%SpaceCurve       = 0

       GridVertex(i)%nbrs(:) = 0
       GridVertex(i)%nbrs_face(:) = 0
       GridVertex(i)%nbrs_wgt(:) = 0
       GridVertex(i)%nbrs_wgt_ghost(:) = 1

       !each elements has one side neighbor (first 4)
       GridVertex(i)%nbrs_ptr(1) = 1
       GridVertex(i)%nbrs_ptr(2) = 2
       GridVertex(i)%nbrs_ptr(3) = 3
       GridVertex(i)%nbrs_ptr(4) = 4
       !don't know about corners yet
       GridVertex(i)%nbrs_ptr(5:num_neighbors+1) = 5

    end do


    !create index table to find neighbors
    call create_index_table(index_table, element_nodes)


    ! side neighbors
    call find_side_neighbors(GridVertex, normal_to_homme_ordering, element_nodes, EdgeWgtP, index_table)
    if (is_zoltan_partition(partmethod) .or. is_zoltan_task_mapping(z2_map_method)) then
        allocate(coord_dim1(p_number_elements))
        allocate(coord_dim2(p_number_elements))
        allocate(coord_dim3(p_number_elements))
    endif

    ! set vertex faces
    do i=1, p_number_elements
       node_numbers = element_nodes(i,:)
       coordinates = p_node_coordinates(node_numbers,:)
       centroid = SUM(coordinates, dim=1)/4.0
       face_center%x = centroid(1)
       face_center%y = centroid(2)
       face_center%z = centroid(3)
       GridVertex(i)%face_number = cube_face_number_from_cart(face_center)

       if (is_zoltan_partition(partmethod) .or. is_zoltan_task_mapping(z2_map_method)) then
                rangex = 1.7 ! ~pi/2
                if (coord_transform_method == SPHERE_COORDS) then
                    coord_dim1(i) = face_center%x
                    coord_dim2(i) = face_center%y
                    coord_dim3(i) = face_center%z
                elseif (coord_transform_method == CUBE_COORDS) then
                    face_no = GridVertex(i)%face_number
                    face_center_2d = cart2cubedsphere(face_center, face_no)

                    if (face_no .EQ. 1) then
                        x = face_center_2d%x
                        y = -(face_center_2d%y)  !+ ((cube_yend-cube_ystart) * face_no)
                        z = 1
                    elseif (face_no .EQ. 2) then
                        x = 1
                        y = -(face_center_2d%y)  !+ ((cube_yend-cube_ystart) * face_no)
                        z = -(face_center_2d%x)  !+ ((cube_xend-cube_xstart) * face_no)
                    elseif(face_no .EQ. 3) then
                        x = -(face_center_2d%x)  !+ ((cube_xend-cube_xstart) * face_no)
                        y = -(face_center_2d%y)  !+ ((cube_yend-cube_ystart) * face_no)
                        z = -1
                    elseif(face_no .EQ. 4) then
                        x = -1
                        y = -(face_center_2d%y) !+ ((cube_yend-cube_ystart) * face_no)
                        z = (face_center_2d%x)  !+ ((cube_xend-cube_xstart) * face_no)
                    elseif(face_no .EQ. 5) then
                        x = (face_center_2d%x) !+ ((cube_xend-cube_xstart) * face_no)
                        y = 1
                        z = (face_center_2d%y) !+ ((cube_yend-cube_ystart) * face_no)
                    elseif(face_no .EQ. 6) then
                        x = (face_center_2d%x) !+ ((cube_xend-cube_xstart) * face_no)
                        y = -1
                        z = -(face_center_2d%y) !+ ((cube_yend-cube_ystart) * face_no)
                    endif

                elseif (coord_transform_method == FACE_2D_LB_COORDS) then
                    face_no = GridVertex(i)%face_number
                    face_center_2d = cart2cubedsphere(face_center, face_no)
                    if (face_no .EQ. 1) then
                        x = (face_center_2d%x)  + rangex * face_no !+ ((cube_xend-cube_xstart) * face_no)
                        y = -(face_center_2d%y)
                        z = 1
                    elseif (face_no .EQ. 2) then
                        y = -(face_center_2d%y)
                        z = 1
                        x = (face_center_2d%x)  + rangex * face_no  !+ ((cube_xend-cube_xstart) * face_no)
                    elseif(face_no .EQ. 3) then
                        x = (face_center_2d%x)   + rangex * face_no !+ ((cube_xend-cube_xstart) * face_no)
                        y = -(face_center_2d%y)
                        z = 1
                    elseif(face_no .EQ. 4) then
                        x = (face_center_2d%x)   + rangex * face_no !+ ((cube_xend-cube_xstart) * face_no)
                        y = -(face_center_2d%y)
                        z = 1
                    elseif(face_no .EQ. 5) then
                        x = (face_center_2d%x)
                        z = 1
                        y = (face_center_2d%y)

                        if (y > 0) then
                            x = x + rangex
                            y = -y + rangex
                        else
                            x = -x + 3 * rangex
                            y = y + rangex
                        endif


                        !if (y > centery + cube_ystart) then
                            !x = x + (cube_xend-cube_xstart)
                            !y = -y + (cube_yend-cube_ystart)
                        !else
                            !x = -x + 3 * (cube_xend-cube_xstart)
                            !y = y - cube_ystart + cube_yend
                        !endif

                    elseif(face_no .EQ. 6) then
                        x = (face_center_2d%x)
                        z = 1
                        y = (face_center_2d%y)
                        if (x > 0) then
                            z = x
                            x = y
                            y = z
                            z = 1
                            y = y - rangex
                            x = x + 2 * rangex
                        else
                            z = x
                            x = y
                            y = z
                            z = 1
                            x = -x + 4 * rangex
                            y = -y - rangex
                        endif

                        !if (x > centerx + cube_xstart) then
                        !    z = x
                        !    x = y
                        !    y = z
                        !    z = 1
                        !    y = y - ((cube_yend-cube_ystart))
                        !    x = x + 2 * ((cube_xend-cube_xstart) )
                        !else
                        !    z = x
                        !    x = y
                        !    y = z
                        !    z = 1
                        !    x = -x + 4 * (cube_xend-cube_xstart)
                        !    y = -y - ((cube_yend-cube_ystart))
                        !endif

                    endif
                else
                    coord_dim1(i) = face_center%x
                    coord_dim2(i) = face_center%y
                    coord_dim3(i) = face_center%z
                endif

                coord_dim1(i) = x
                coord_dim2(i) = y
                coord_dim3(i) = z
        endif
    end do

    if (coord_transform_method == SPHERE_COORDS) then
        coord_dimension = 3
    elseif (coord_transform_method == CUBE_COORDS) then
        coord_dimension = 3
    elseif (coord_transform_method == FACE_2D_LB_COORDS) then
        coord_dimension = 2
    else
        coord_dimension = 3
    endif

    ! set side neighbor faces
    do i=1, p_number_elements
       do j=1,4 !look at each side
          k = normal_to_homme_ordering(j)
          loc =  GridVertex(i)%nbrs_ptr(k)
          ll = GridVertex(i)%nbrs(loc)
          GridVertex(i)%nbrs_face(loc) = GridVertex(ll)%face_number
       end do
    end do

    ! find corner neighbor and faces (weights added also)
    call find_corner_neighbors  (GridVertex, normal_to_homme_ordering, element_nodes, CornerWgt, index_table)


    !done with the index table
    deallocate(index_table)


    call initgridedge(GridEdge,GridVertex)
    do i=1,SIZE(GridEdge)
       call CubeSetupEdgeIndex(GridEdge(i))
    enddo

    if(partmethod .eq. SFCURVE) then
        call initialize_space_filling_curve(GridVertex, element_nodes)
    endif
#endif
  end subroutine MeshCubeTopologyCoords

!======================================================================
! subroutine MeshCubeTopology
!======================================================================
   subroutine MeshCubeTopology(GridEdge, GridVertex)
    use parallel_mod,           only : abortmp
    use dimensions_mod,         only : np,  max_elements_attached_to_node
    use coordinate_systems_mod, only : cartesian3D_t, cube_face_number_from_cart, cube_face_number_from_sphere
    use gridgraph_mod,          only : GridVertex_t
    use gridgraph_mod,          only : GridEdge_t
    use cube_mod,               only : CubeSetupEdgeIndex
    use gridgraph_mod,          only : initgridedge, num_neighbors
    use control_mod,            only : north, south, east, west, neast, seast, swest, nwest

    implicit none
    type (GridEdge_t),   intent(inout) :: GridEdge(:)
    type (GridVertex_t), intent(inout) :: GridVertex(:)

#ifndef HOMME_WITHOUT_PIOLIBRARY
    real(kind=real_kind)             :: coordinates(4,3)
    real(kind=real_kind)             :: centroid(3)
    type (cartesian3D_t)             :: face_center

    integer                          :: i, j, k, ll, m, loc
    integer                          :: element_nodes(p_number_elements, 4)
    integer                          :: EdgeWgtP,CornerWgt
    integer                          :: normal_to_homme_ordering(8)
    integer                          :: node_numbers(4)
    integer, allocatable             :: index_table(:,:)

    normal_to_homme_ordering(1) = south
    normal_to_homme_ordering(2) =  east
    normal_to_homme_ordering(3) = north
    normal_to_homme_ordering(4) =  west
    normal_to_homme_ordering(5) = swest
    normal_to_homme_ordering(6) = seast
    normal_to_homme_ordering(7) = neast
    normal_to_homme_ordering(8) = nwest

    if (SIZE(GridVertex) /= p_number_elements) then
       call abortmp('MeshCubeTopology: Element count check failed in exodus_mesh. &
            &Vertex array length not equal to number of elements.')
    end if
    if (p_number_elements_per_face /= p_number_elements) then
       call abortmp('MeshCubeTopology: Element count check failed in exodus_mesh. &
            &Element array length not equal to sum of face.')
    end if

    EdgeWgtP = np
    CornerWgt = 1


    call mesh_connectivity (element_nodes)

    do i=1, p_number_elements
       GridVertex(i)%number           = i
       GridVertex(i)%face_number      = 0
       GridVertex(i)%processor_number = 0
       GridVertex(i)%SpaceCurve       = 0

       GridVertex(i)%nbrs(:) = 0
       GridVertex(i)%nbrs_face(:) = 0
       GridVertex(i)%nbrs_wgt(:) = 0
       GridVertex(i)%nbrs_wgt_ghost(:) = 1

       !each elements has one side neighbor (first 4)
       GridVertex(i)%nbrs_ptr(1) = 1
       GridVertex(i)%nbrs_ptr(2) = 2
       GridVertex(i)%nbrs_ptr(3) = 3
       GridVertex(i)%nbrs_ptr(4) = 4
       !don't know about corners yet
       GridVertex(i)%nbrs_ptr(5:num_neighbors+1) = 5

    end do

    !create index table to find neighbors
    call create_index_table(index_table, element_nodes)

    ! side neighbors
    call find_side_neighbors(GridVertex, normal_to_homme_ordering, element_nodes, EdgeWgtP, index_table)

    ! set vertex faces
    do i=1, p_number_elements
       node_numbers = element_nodes(i,:)
       coordinates = p_node_coordinates(node_numbers,:)
       centroid = SUM(coordinates, dim=1)/4.0
       face_center%x = centroid(1)
       face_center%y = centroid(2)
       face_center%z = centroid(3)
       GridVertex(i)%face_number = cube_face_number_from_cart(face_center)
    end do

    ! set side neighbor faces
    do i=1, p_number_elements
       do j=1,4 !look at each side
          k = normal_to_homme_ordering(j)
          loc =  GridVertex(i)%nbrs_ptr(k)
          ll = GridVertex(i)%nbrs(loc)
          GridVertex(i)%nbrs_face(loc) = GridVertex(ll)%face_number
       end do
    end do

    ! find corner neighbor and faces (weights added also)
    call find_corner_neighbors  (GridVertex, normal_to_homme_ordering, element_nodes, CornerWgt, index_table)

    !done with the index table
    deallocate(index_table)


    call initgridedge(GridEdge,GridVertex)
    do i=1,SIZE(GridEdge)
       call CubeSetupEdgeIndex(GridEdge(i))
    enddo

    call initialize_space_filling_curve(GridVertex, element_nodes)
#endif
  end subroutine MeshCubeTopology

!======================================================================
! subroutine MeshSetCoordinates(elem)
!======================================================================

  subroutine MeshSetCoordinates(elem)
    use element_mod, only : element_t
    use parallel_mod, only : abortmp
    use coordinate_systems_mod, only   : cartesian3D_t, cartesian2d_t, spherical_polar_t, &
                                         change_coordinates, sphere2cubedsphere
    implicit none

    type (element_t), intent(inout)  :: elem(:)
#ifndef HOMME_WITHOUT_PIOLIBRARY
    integer                          :: connectivity(p_number_elements,4)
    integer                          :: node_multiplicity(p_number_nodes)
    integer                          :: face_no, i, k, l
    integer                          :: number
    integer                          :: node_num(4)
    real(kind=real_kind)             :: coordinates(4,3)
    real(kind=real_kind)             :: cube_coor  (4,2)

    real(kind=real_kind)             :: x_double
    real                             :: x_real
    type (cartesian2d_t)             :: cart2

    connectivity     =0
    node_multiplicity=0

    call mesh_connectivity (connectivity)

    do k=1,p_number_elements
       node_num = connectivity(k,:)
       node_multiplicity(node_num(:)) = node_multiplicity(node_num(:)) + 1
    end do

    do k=1,SIZE(elem)
       number  = elem(k)%vertex%number
       face_no = elem(k)%vertex%face_number
       node_num = connectivity(number,:)
       coordinates = p_node_coordinates(node_num,:)

      if (6 == p_number_blocks) then
         call cube_to_cube_coordinates   (cube_coor, coordinates, face_no)
      else
         call sphere_to_cube_coordinates (cube_coor, coordinates, face_no)
      end if
!      elem(k)%node_numbers         = node_num
!      elem(k)%node_multiplicity(:) = node_multiplicity(node_num(:))
      elem(k)%corners(:)%x         = cube_coor(:,1)
      elem(k)%corners(:)%y         = cube_coor(:,2)
    end do
#endif
  end subroutine MeshSetCoordinates

!======================================================================
!function MeshCubeEdgeCount()
!======================================================================
  function MeshCubeEdgeCount()  result(nedge)
    use parallel_mod, only : abortmp
    implicit none
    integer                     :: nedge
#ifndef HOMME_WITHOUT_PIOLIBRARY
    if (0 == p_number_blocks)  call abortmp('MeshCubeEdgeCount called before MeshOpenMesh')
    if (MeshUseMeshFile) then
      ! should be the same as SUM(SIZE(GridVertex(i)%nbrs(j)%n),i=1:p_number_elements,j=1:nInnerElemEdge)
      ! the total number of neighbors.
      nedge = p_number_neighbor_edges
    else
      call abortmp('Error in MeshCubeEdgeCount: Should not call for non-exodus mesh file.')
    endif
#endif
  end function MeshCubeEdgeCount

  function MeshCubeElemCount()  result(nelem)
    use parallel_mod, only : abortmp
    implicit none
    integer                     :: nelem
#ifndef HOMME_WITHOUT_PIOLIBRARY
    if (0 == p_number_blocks)  call abortmp('MeshCubeElemCount called before MeshOpenMesh')
    if (MeshUseMeshFile) then
      nelem = p_number_elements
    else
      call abortmp('Error in MeshCubeElemCount: Should not call for non-exodus mesh file.')
    end if
#endif
  end function MeshCubeElemCount

#ifndef HOMME_WITHOUT_PIOLIBRARY
  subroutine test_private_methods
    implicit none
    integer                          :: element_nodes(p_number_elements, 4)
    call mesh_connectivity (element_nodes)
  end subroutine test_private_methods
#endif

end module mesh_mod
