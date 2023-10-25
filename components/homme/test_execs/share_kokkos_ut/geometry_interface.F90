module geometry_interface_mod

  use element_mod,    only : element_t
  use gridgraph_mod,  only : GridVertex_t, GridEdge_t
  use hybrid_mod,     only : hybrid_t
  use kinds,          only : real_kind
  use metagraph_mod,  only : MetaVertex_t
  use parallel_mod,   only : parallel_t, abortmp

  implicit none

  type (GridVertex_t), dimension(:), allocatable :: GridVertex
  type (GridEdge_t),   dimension(:), allocatable :: GridEdge

  type (MetaVertex_t) :: MetaVertex
  type (parallel_t)   :: par
  type (hybrid_t)     :: hybrid

  type (element_t), dimension(:), pointer, public :: elem

  integer :: nelemd

  public :: initmp_f90
  public :: init_cube_geometry_f90
  public :: init_connectivity_f90
  public :: init_c_connectivity_f90
  public :: cleanup_geometry_f90

#include <mpif.h>

contains

  subroutine initmp_f90 () bind(c)
    use parallel_mod, only : initmp
    use hybrid_mod,   only : hybrid_create

    par = initmp()

    ! No horizontal openmp, for now
    hybrid = hybrid_create(par,0,1)
  end subroutine initmp_f90


  subroutine init_cube_geometry_f90 (ne_in) bind(c)
    use iso_c_binding,   only : c_int
    use cube_mod,        only: CubeTopology, CubeElemCount, CubeEdgeCount
    use dimensions_mod,  only: nelem, nelemd, npart, np
    use dimensions_mod,  only: ne
    use gridgraph_mod,   only: allocate_gridvertex_nbrs
    use spacecurve_mod,  only: genspacepart
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: ne_in
    !
    ! Locals
    !
    integer :: ie, num_edges

    ne = ne_in
    nelem = CubeElemCount()
    num_edges = CubeEdgeCount()

    if (nelem < par%nprocs) then
       call abortmp('Error: too many MPI tasks. Run the test with less MPI ranks, or increase ne.')
    end if

    ! Allocate
    allocate (GridVertex(nelem))
    allocate (GridEdge(num_edges))

    do ie=1,nelem
      call allocate_gridvertex_nbrs(GridVertex(ie))
    enddo

    ! Generate mesh connectivity
    call CubeTopology(GridEdge, GridVertex)

    ! Set number of partitions before partitioning mesh
    npart = par%nprocs

    ! Partition mesh among processes
    call genspacepart(GridEdge, GridVertex)
  end subroutine init_cube_geometry_f90

  subroutine init_connectivity_f90 () bind(c)
    use dimensions_mod, only : nelem, nelemd, nlev, nlevp
    use element_mod,    only : allocate_element_desc
    use metagraph_mod,  only : initMetaGraph, LocalElemCount
    use parallel_mod,   only : iam, MPI_MAX, MPIinteger_t
    use schedtype_mod,  only : Schedule
    use schedule_mod,   only : genEdgeSched
    !
    ! Locals
    !
    integer :: Global2Local(nelem)
    integer, pointer :: map_ptr(:)
    integer :: nelemd_max, ie, ip, ierr
    integer, target :: my_num_elems
    integer, target, allocatable :: num_elems(:)
    integer, pointer :: my_num_elems_ptr, num_elems_ptr(:)

    call initMetaGraph(iam,MetaVertex,GridVertex,GridEdge)

    nelemd = LocalElemCount(MetaVertex)

    if (nelemd .le. 0) then
      call abortmp ("Something went wrong during domain partitioning, &
                     and a process is left with no elements. Please, &
                     run the test with more elements, or less MPI ranks.")
    endif

    allocate (num_elems(par%nprocs))
    my_num_elems_ptr => my_num_elems
    num_elems_ptr => num_elems
    my_num_elems = nelemd

    allocate (elem(nelemd))
    call allocate_element_desc(elem)

    allocate (Schedule(1))
    call genEdgeSched(elem,iam,Schedule(1),MetaVertex)

    ! Defaults all local ids to 0 (meaning not on this process)
    Global2Local = 0
    do ie=1,SIZE(MetaVertex%members)
      Global2Local(MetaVertex%members(ie)%number) = ie
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,Global2Local,nelem,MPIinteger_t,MPI_MAX,par%comm,ierr)

    ! Pass info to C
    call init_c_connectivity_f90 (nelemd,Global2Local)

  end subroutine init_connectivity_f90

  subroutine init_c_connectivity_f90 (nelemd,Global2Local)
    use gridgraph_mod,  only : GridEdge_t
    use dimensions_mod, only : nelem
    interface
      subroutine init_connectivity (num_local_elems, max_corner_elems) bind (c)
        use iso_c_binding, only : c_int
        !
        ! Inputs
        !
        integer (kind=c_int), intent(in) :: num_local_elems, max_corner_elems
      end subroutine init_connectivity
      subroutine finalize_connectivity () bind(c)
      end subroutine finalize_connectivity

      subroutine add_connection (first_lid,  first_gid,  first_pos,  first_pid, &
                                 second_lid, second_gid, second_pos, second_pid) bind(c)
        use iso_c_binding, only : c_int
        !
        ! Inputs
        !
        integer (kind=c_int), intent(in) :: first_lid,  first_gid,  first_pos,  first_pid
        integer (kind=c_int), intent(in) :: second_lid, second_gid, second_pos, second_pid
      end subroutine add_connection
    end interface

    !
    ! Inputs
    !
    integer, dimension(nelem), intent(in) :: Global2Local
    integer, intent(in) :: nelemd
    !
    ! Locals
    !
    integer :: ie, num_edges, ierr, max_corner_elems
    type(GridEdge_t) :: e

    max_corner_elems = 1 ! always structured cubed-sphere in unit tests
    call init_connectivity(nelemd, max_corner_elems)

    num_edges = SIZE(GridEdge)
    do ie=1,num_edges
      e = GridEdge(ie)
      call add_connection(Global2Local(e%head%number),e%head%number,e%head_dir,e%head%processor_number, &
                          Global2Local(e%tail%number),e%tail%number,e%tail_dir,e%tail%processor_number)
    enddo

    call finalize_connectivity()
  end subroutine init_c_connectivity_f90

  subroutine cleanup_geometry_f90 () bind(c)
    use schedtype_mod,  only : Schedule
    use gridgraph_mod,  only : deallocate_gridvertex_nbrs
    use dimensions_mod, only : nelem
    !
    ! Locals
    !
    integer :: ie

    do ie=1,nelem
      call deallocate_gridvertex_nbrs(GridVertex(ie))
    enddo
    deallocate(elem)
    deallocate(GridVertex)
    deallocate(GridEdge)

    deallocate(Schedule(1)%SendCycle)
    deallocate(Schedule(1)%RecvCycle)
    deallocate(Schedule(1)%MoveCycle)
    deallocate(Schedule(1)%pIndx)
    deallocate(Schedule(1)%gIndx)
    deallocate(Schedule(1)%Local2Global)
    deallocate(Schedule)
  end subroutine cleanup_geometry_f90

end module geometry_interface_mod
