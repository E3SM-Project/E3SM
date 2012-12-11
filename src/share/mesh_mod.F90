#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module mesh_mod

 use kinds, only : real_kind, long_kind
  use physical_constants, only : DD_PI
  use control_mod, only : MAX_FILE_LEN

  use netcdf ! _EXTERNAL

  implicit none
  logical, public           :: MeshUseMeshFile = .false.
  public  :: MeshOpen           ! Must be called first

  
  integer, parameter :: MXSTLN = 32

  ! ===============================
  ! Public methods for mesh_mod
  ! ===============================
  
  public  :: MeshCubeEdgeCount  ! called anytime afer MeshOpen
  public  :: MeshCubeElemCount  ! called anytime afer MeshOpen
  public  :: MeshCubeTopology   ! called afer MeshOpen
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


  private :: find_side_neighbor_information
  private :: find_side_neighbors
  private :: find_corner_neighbors
  private :: find_nodal_neighbors
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


contains

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
      print *, 'maxval(node_multiplicity)', maxval(node_multiplicity),' and max_elements_attached_to_node ',max_elements_attached_to_node
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
!  function find_side_neighbor_information
!
!   return the element number and direction of the side neighbor (n,s,e,w)
!   that shares side specified in nodes
!======================================================================
  function find_side_neighbor_information(me, nodes, nel, element_nodes, direction) result(neighbor)
    use parallel_mod, only : abortmp
    integer, intent(in)       :: me !my element
    integer, intent(in)       :: nodes(2) !nodes defining one of the sides
    integer, intent(in)       :: nel !number of elements
    integer, intent(out)      :: direction
    integer, intent(in)       :: element_nodes(nel, 4) 
    integer                   :: i,j
    integer                   :: neighbor
    integer                   :: side_nodes(2)
    neighbor = 0
    do i=1,nel !loop through ALL elements
       if (i /= me) then
          do j=1,4 !loop through each elements four sides
             side_nodes(1) = element_nodes(i,j)
             side_nodes(2) = element_nodes(i,mod(j,4)+1)
             if ( (side_nodes(1) == nodes(2) .and. side_nodes(2) == nodes(1)) .or. &
                  (side_nodes(1) == nodes(1) .and. side_nodes(2) == nodes(2)) ) then
                neighbor = i
                direction = j
                return  !found a match - we can return
             end if
          end do
       end if
    end do
    if (neighbor == 0) call abortmp('find_side_neighbor_information: Neighbor not found! Every side should have a neighbor.') 
  end function find_side_neighbor_information
  
!======================================================================
! function find_nodal_neighbors
!======================================================================
 
  function find_nodal_neighbors(node, element_nodes, elements) result(multiplicity)
    use parallel_mod, only : abortmp
    use dimensions_mod, only : max_elements_attached_to_node
    implicit none
    ! parameters 
    integer, intent(in)    :: node
    integer, intent(in)    :: element_nodes(p_number_elements, 4)
    integer, intent(out)   :: elements      (2*max_elements_attached_to_node)
    ! return value
    integer                :: multiplicity
    ! local variables
    integer                :: i,j

    multiplicity = 0
    do i=1,p_number_elements  !loop through all elements
       do j=1,4  !loop through each node
          if ( element_nodes(i,j) == node ) then
             multiplicity = multiplicity + 1
             elements(multiplicity) = i !keep track of matches
          end if
       end do
    end do
    if (multiplicity < 3 .or. max_elements_attached_to_node < multiplicity) then
      call abortmp('find_nodal_neighbors: Number of elements attached to node less than 3 or greater than maximum.')
      !at least 3 elements attached to each node - no more than 7 for mesh (in dimensions_mod)
    endif
    
  end function find_nodal_neighbors

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
! subroutine find_side_neighbors()
!
! find the element neighbors to the n,s,e,w and put them in GridVertex_t
!  (only 1 neighbor to the n,s,e,w
!======================================================================
  subroutine find_side_neighbors (GridVertex, normal_to_homme_ordering, element_nodes, edge_wgt)
    use coordinate_systems_mod, only : cartesian3D_t
    use gridgraph_mod, only   : GridVertex_t
    use parallel_mod, only : abortmp
    implicit none
    integer             ,  intent(in)    :: normal_to_homme_ordering(8)
    integer             ,  intent(in)    :: element_nodes(p_number_elements, 4)
    integer             ,  intent(in)    :: edge_wgt
    type (GridVertex_t) ,  intent(inout) :: GridVertex(:)
    
    integer                              :: side_nodes(2)
    integer                              :: neighbor, direction, init_size
    integer                              :: j,k,ll
    integer                              :: loc

    if (0 == p_number_blocks)  call abortmp('find_side_neighbors called before MeshOpen')

    
    do k=1,p_number_elements  ! for each element

       !set the side weights
       GridVertex(k)%nbrs_wgt(1:4) = edge_wgt

       do ll=1,4 !loop through the four sides

          j = normal_to_homme_ordering(ll)
          loc = GridVertex(k)%nbrs_ptr(j)

          if (GridVertex(k)%nbrs(loc) == 0) then  !if not set yet, then

             side_nodes(1) = element_nodes(k, ll)
             side_nodes(2) = element_nodes(k, mod(ll,4)+1)
             neighbor = find_side_neighbor_information(k, side_nodes,     &
                  p_number_elements, &
                  element_nodes,     &
                  direction)

             GridVertex(k)%nbrs(loc) = neighbor         
             
             j = normal_to_homme_ordering(direction)
             loc = GridVertex(neighbor)%nbrs_ptr(j)
             GridVertex(neighbor)%nbrs(loc)= k         

          endif
       enddo
    enddo
    
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
       call abortmp('smallest_diameter_element:Element count check failed in exodus_mesh. Connectivity array length not equal to number of elements.')
    end if
    if ( p_number_elements_per_face /= p_number_elements) then
       call abortmp('smallest_diameter_element: Element count check failed in exodus_mesh. Element array length not equal to sum of face.')
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
    real,                intent(out) :: centroids       (p_number_elements,2)
    real(kind=real_kind)             :: coordinates(4,3) 
    real(kind=real_kind)             :: cube_coor  (4,2) 
    integer                          :: i, j, node_numbers(4)
    
    if (0 == p_number_blocks)  call abortmp('cube_face_element_centroids called before MeshOpen')
    if (SIZE(element_nodes,dim=1) /= p_number_elements) then
       call abortmp('cube_face_element_centroids:Element count check failed in exodus_mesh. Connectivity array length not equal to number of elements.')
    end if
    if ( p_number_elements_per_face /= p_number_elements ) then
       call abortmp('cube_face_element_centroids: Element count check failed in exodus_mesh. Element array length not equal to sum of face.')
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
    use spacecurve_mod, only  : GenspaceCurve
    
    implicit none
    type (GridVertex_t), intent(inout) :: GridVertex(:)
    integer            , intent(in)    :: element_nodes(:,:)
    
    integer,allocatable                :: Mesh2(:,:),Mesh2_map(:,:),sfcij(:,:)
    
    real                               :: centroids(p_number_elements,2)
    integer                            :: face_numbers(p_number_elements)
    real                               :: x, y, h
    integer                            :: i, j, i2, j2, ne, ne2
    integer                            :: sfc_index, face, nelem
    
    if (SIZE(GridVertex) /= p_number_elements) then
       call abortmp('initialize_space_filling_curve:Element count check failed in exodus_mesh. Vertex array length not equal to number of elements.')
    end if
    if (SIZE(element_nodes,dim=1) /= p_number_elements) then
       call abortmp('initialize_space_filling_curve:Element count check failed in exodus_mesh. Connectivity array length not equal to number of elements.')
    end if
    
    face_numbers(:) = GridVertex(:)%face_number
    h = smallest_diameter_element    (                         element_nodes)

    call cube_face_element_centroids (centroids, face_numbers, element_nodes)
    
    if (h<.00001) call abortmp('initialize_space_filling_curve: Unreasonably small element found. less than .00001')
    
    ne = CEILING(0.5*DD_PI/(h/2));
    
    ! find the smallest ne2 which is a power of 2 and ne2>ne
    ne2=2**ceiling( log(real(ne))/log(2d0) )
    if (ne2<ne) call abortmp('initialize_space_filling_curve: Fatel SFC error')
    
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
    
    
    GridVertex(:)%SpaceCurve=-1
    sfc_index   = 0
    do face = 1,nfaces
       ! associate every element on the ne x ne mesh (Mesh)
       ! with its closest element on the ne2 x ne2 mesh (Mesh2)
       ! Store this as a map from Mesh2 -> Mesh in Mesh2_map.
       ! elements in Mesh2 which are not mapped get assigned a value of 0
       Mesh2_map=0
       do i=1,p_number_elements
          if (face_numbers(i) == face ) then
             x = centroids(i,1)
             y = centroids(i,2)
             ! map this element to an (i2,j2) element
             ! [ -DD_PI/4, DD_PI/4 ]  -> [ 0, ne2 ]
             i2=nint( (0.5 + 2.0*x/DD_PI)*ne2 + .5 )
             j2=nint( (0.5 + 2.0*y/DD_PI)*ne2 + .5 )
             if (face == 4 .or. face == 6 )               i2 = ne2-i2+1
             if (face == 1 .or. face == 2 .or. face == 6) j2 = ne2-j2+1
             if (i2<1  ) i2=1
             if (i2>ne2) i2=ne2
             if (j2<1  ) j2=1
             if (j2>ne2) j2=ne2
             Mesh2_map(i2,j2)=i
          end if
       end do
       
       ! generate a SFC for Mesh with the same ordering as the 
       ! elements in Mesh2 which map to Mesh.
       do j=0,ne2*ne2-1
          i2=sfcij(j,1)
          j2=sfcij(j,2)
          i=Mesh2_map(i2,j2)
          if (i/=0) then
             ! (i2,j2) element maps to element
             GridVertex(i)%SpaceCurve=sfc_index
             sfc_index=sfc_index+1
          endif
       enddo
    enddo
    deallocate(Mesh2)
    deallocate(Mesh2_map)
    deallocate(sfcij)
    
    if (minval(GridVertex(:)%SpaceCurve) == -1) then
       do i=1,p_number_elements
          if (-1==GridVertex(i)%SpaceCurve) then
             write (*,*) " Error in projecting element ",i," to space filling curve."
             write (*,*) " Face:",face_numbers(i)
             write (*,*) " Centroid:",centroids(i,:)
          end if
       end do
       call abortmp('initialize_space_filling_curve: Vertex not on SpaceCurve')
    end if
    
  end subroutine initialize_space_filling_curve
  
!======================================================================
! subroutine find_corner_neighbors
!======================================================================

  subroutine find_corner_neighbors  (GridVertex, normal_to_homme_ordering, element_nodes, corner_wgt)
    use parallel_mod,           only : abortmp
    use gridgraph_mod,          only : GridVertex_t
    use dimensions_mod,         only : max_elements_attached_to_node
    implicit none
    
    type (GridVertex_t), intent(inout) :: GridVertex(:)
    integer            , intent(in)  :: normal_to_homme_ordering(8)
    integer            , intent(in)  :: element_nodes(p_number_elements, 4)
    integer            , intent(in) :: corner_wgt

    integer                          :: node_elements (2*max_elements_attached_to_node)
    integer                          :: elem_neighbor (4*max_elements_attached_to_node)
    integer                          :: nbr_cnt(4)
    integer                          :: multiplicity, elem_nbr_start, start
    integer                          :: i, j, k, ll, loc, cnt, jj
    
    do i=1, p_number_elements  !loop through all elements
       node_elements(:) = 0
       elem_neighbor(:) = 0
       elem_nbr_start = 0
       nbr_cnt(:) = 0

       do j=1,4  !check each of the 4 nodes at the element corners
          multiplicity = find_nodal_neighbors(element_nodes(i,j), element_nodes, node_elements)
          !returns all of the element neighbors to that node (in node_elements)
          k = 0
          do ll=1,multiplicity !node_elements that are not already side neighbors (or myself)
                              !are corner neighbors
             if (                         i /= node_elements(ll) .and. &
                  GridVertex(i)%nbrs(1) /= node_elements(ll) .and. &
                  GridVertex(i)%nbrs(2) /= node_elements(ll) .and. &
                  GridVertex(i)%nbrs(3) /= node_elements(ll) .and. &
                  GridVertex(i)%nbrs(4) /= node_elements(ll)) then   
                k = k + 1
                elem_neighbor(elem_nbr_start + k) = node_elements(ll)
             end if
          end do ! end of ll loop for multiplicity

          !keep track of where we are starting in elem_neighbor for each j
          elem_nbr_start = elem_nbr_start + k
          nbr_cnt(j) = k !how many neighbors in this corner
       end do ! end of j loop through 4 nodes


       !now that we have done the 4 corners we can populate nbrs and nbrs_ptr 
       ! (in cons. order)
       !also we can add the corner weight

       do j=5,8 !loop through 4 corners
          elem_nbr_start = 1
          !easiest to do these in ascending order - find loc
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

       end do
       
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
    type (GridEdge_t),   intent(out) :: GridEdge(:)
    type (GridVertex_t), intent(out) :: GridVertex(:)

    real(kind=real_kind)             :: coordinates(4,3) 
    real(kind=real_kind)             :: centroid(3)
    type (cartesian3D_t)             :: face_center

    integer                          :: i, j, k, ll, m, loc
    integer                          :: element_nodes(p_number_elements, 4)
    integer                          :: EdgeWgtP,CornerWgt
    integer                          :: normal_to_homme_ordering(8)
    integer                          :: node_numbers(4)
    integer                          :: max_size

    normal_to_homme_ordering(1) = south
    normal_to_homme_ordering(2) =  east
    normal_to_homme_ordering(3) = north 
    normal_to_homme_ordering(4) =  west
    normal_to_homme_ordering(5) = swest
    normal_to_homme_ordering(6) = seast
    normal_to_homme_ordering(7) = neast
    normal_to_homme_ordering(8) = nwest

    if (SIZE(GridVertex) /= p_number_elements) then
       call abortmp('MeshCubeTopology: Element count check failed in exodus_mesh. Vertex array length not equal to number of elements.')
    end if
    if (p_number_elements_per_face /= p_number_elements) then
       call abortmp('MeshCubeTopology: Element count check failed in exodus_mesh. Element array length not equal to sum of face.')
    end if

    EdgeWgtP = np
    CornerWgt = 1


    call mesh_connectivity (element_nodes)


    !INITIALIZE the Grid Vertex for each element
    !max number of neighbors:
    !set init_size  (4 side neighbors + 4 corners (-3 is for 2 sides and element itself)
    max_size = 4 + 4*(max_elements_attached_to_node -3)    

    !TODO: can either reduce storage at the end if have alloced too much
    ! or can alloc less and then increase if needed

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


    ! side neighbors 
    call find_side_neighbors  (GridVertex, normal_to_homme_ordering, element_nodes, EdgeWgtP)

   
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
    call find_corner_neighbors  (GridVertex, normal_to_homme_ordering, element_nodes, CornerWgt)

  
    call initgridedge(GridEdge,GridVertex) 
    do i=1,SIZE(GridEdge)
       call CubeSetupEdgeIndex(GridEdge(i)) 
    enddo

    call initialize_space_filling_curve(GridVertex, element_nodes)
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
      elem(k)%node_numbers         = node_num 
      elem(k)%node_multiplicity(:) = node_multiplicity(node_num(:))
      elem(k)%corners(:)%x         = cube_coor(:,1)
      elem(k)%corners(:)%y         = cube_coor(:,2)
    end do
  end subroutine MeshSetCoordinates

!======================================================================
!function MeshCubeEdgeCount()
!======================================================================
  function MeshCubeEdgeCount()  result(nedge)
    use parallel_mod, only : abortmp
    implicit none
    integer                     :: nedge
    if (0 == p_number_blocks)  call abortmp('MeshCubeEdgeCount called before MeshOpenMesh')
    if (MeshUseMeshFile) then
      ! should be the same as SUM(SIZE(GridVertex(i)%nbrs(j)%n),i=1:p_number_elements,j=1:nInnerElemEdge)
      ! the total number of neighbors.
      nedge = p_number_neighbor_edges
    else
      call abortmp('Error in MeshCubeEdgeCount: Should not call for non-exodus mesh file.')
    endif

  end function MeshCubeEdgeCount

  function MeshCubeElemCount()  result(nelem)
    use parallel_mod, only : abortmp
    implicit none
    integer                     :: nelem
    if (0 == p_number_blocks)  call abortmp('MeshCubeElemCount called before MeshOpenMesh')
    if (MeshUseMeshFile) then
      nelem = p_number_elements
    else
      call abortmp('Error in MeshCubeElemCount: Should not call for non-exodus mesh file.')
    end if
  end function MeshCubeElemCount

  subroutine test_private_methods
    implicit none
    integer                          :: element_nodes(p_number_elements, 4)
    call mesh_connectivity (element_nodes)
  end subroutine test_private_methods


end module mesh_mod


