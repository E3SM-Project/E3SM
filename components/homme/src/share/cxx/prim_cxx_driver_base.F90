#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
module prim_cxx_driver_base

  use prim_driver_base, only: deriv1

  implicit none


  private

  public :: prim_init1
  public :: prim_finalize

  private :: generate_global_to_local
  public  :: init_cxx_connectivity
  public  :: setup_element_pointers

#include <mpif.h>

  contains

  subroutine prim_init1(elem, par, dom_mt, tl)
    use iso_c_binding,    only : c_int, c_loc
    use derivative_mod,   only : derivinit
    use dimensions_mod,   only : nelemd, np
    use domain_mod,       only : domain1d_t
    use element_mod,      only : element_t
    use kinds,            only : iulog, real_kind
    use parallel_mod,     only : parallel_t
    use time_mod,         only : TimeLevel_t, TimeLevel_init
    use prim_driver_base, only : prim_init1_geometry, prim_init1_elem_arrays, &
                                 prim_init1_cleanup, prim_init1_buffers,      &
                                 MetaVertex, GridEdge, deriv1
#ifndef CAM
    use prim_driver_base, only : prim_init1_no_cam
#endif

    interface
      subroutine reset_cxx_comm (f_comm) bind(c)
        use iso_c_binding, only: c_int
        !
        ! Inputs
        !
        integer(kind=c_int), intent(in) :: f_comm
      end subroutine reset_cxx_comm
      subroutine initialize_hommexx_session() bind(c)
      end subroutine initialize_hommexx_session

    end interface
    !
    ! Inputs
    !
    type (element_t),   pointer     :: elem(:)
    type (parallel_t),  intent(in)  :: par
    type (domain1d_t),  pointer     :: dom_mt(:)
    type (timelevel_t), intent(out) :: tl

    ! Initialize MPI comm in C++
    call reset_cxx_comm (INT(par%comm,c_int))

    ! Initialize kokkos before any environment changes from the Fortran
    call initialize_hommexx_session()

#ifndef CAM
    call prim_init1_no_cam(par)
#endif

    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call derivinit(deriv1)

    ! ==================================
    ! Initialize and partition the geometry
    ! ==================================
    call prim_init1_geometry(elem,par,dom_mt)

    ! ==================================
    ! Initialize C++ mpi communication structures
    ! ==================================
    call init_cxx_connectivity(nelemd,GridEdge,MetaVertex,par)

    ! Cleanup the tmp stuff used in prim_init1_geometry
    call prim_init1_cleanup()

    ! ==================================
    ! Initialize the buffers for exchanges
    ! ==================================
    call prim_init1_buffers(elem,par)

    ! ==================================
    ! Initialize element pointers
    ! ==================================
    call setup_element_pointers(elem)

    ! ==================================
    ! Initialize element arrays (fluxes and state)
    ! ==================================
    call prim_init1_elem_arrays(elem,par)

    ! Initialize the time levels
    call TimeLevel_init(tl)

    if(par%masterproc) write(iulog,*) 'end of prim_init1'
  end subroutine prim_init1

  subroutine prim_finalize ()
    use prim_driver_base, only: prim_finalize_base=>prim_finalize
    interface
      subroutine finalize_hommexx_session() bind(c)
      end subroutine finalize_hommexx_session
    end interface

    call prim_finalize_base()

    ! This call lets C++ destroy the singleton containing all the views,
    ! and finalize the Kokkos execution space.
    call finalize_hommexx_session()
  end subroutine prim_finalize

  subroutine init_cxx_connectivity (nelemd, GridEdge, MetaVertex, par)
    use dimensions_mod, only : nelem
    use gridgraph_mod,  only : GridEdge_t
    use metagraph_mod,  only : MetaVertex_t
    use parallel_mod,   only : parallel_t
    !
    ! Interfaces
    !
    interface
      subroutine init_connectivity (num_local_elems) bind (c)
        use iso_c_binding, only : c_int
        !
        ! Inputs
        !
        integer (kind=c_int), intent(in) :: num_local_elems
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
    integer, intent(in) :: nelemd
    type(GridEdge_t),   intent(in) :: GridEdge(:)
    type (parallel_t),  intent(in) :: par
    type(MetaVertex_t), intent(in) :: MetaVertex
    !
    ! Locals
    !
    integer :: Global2Local(nelem)
    integer :: ie, num_edges
    type(GridEdge_t) :: e

    ! Generate a global-to-local map of the meta vertices
    call generate_global_to_local(MetaVertex,Global2Local,par)

    ! Initialize C++ connectivity structure
    call init_connectivity(nelemd)

    ! Add all connections to the C++ structure
    num_edges = SIZE(GridEdge)
    do ie=1,num_edges
      e = GridEdge(ie)
      call add_connection(Global2Local(e%head%number),e%head%number,e%head_dir,e%head%processor_number, &
                          Global2Local(e%tail%number),e%tail%number,e%tail_dir,e%tail%processor_number)
    enddo

    call finalize_connectivity()
  end subroutine init_cxx_connectivity

  subroutine setup_element_pointers (elem)
    use element_mod,    only : element_t
    use element_state,  only : allocate_element_arrays, setup_element_pointers_ie
    use dimensions_mod, only : nelemd
    !
    ! Inputs
    !
    type (element_t), intent(inout) :: elem(:)
    !
    ! Locals
    !
    integer :: ie

    call allocate_element_arrays(nelemd)

    do ie=1,nelemd
      call setup_element_pointers_ie(ie, elem(ie)%state, elem(ie)%derived, elem(ie)%accum)
    enddo

  end subroutine setup_element_pointers

!!!!!!!!!!!!!!!!!!!!!!! PRIVATE SUBROUTINES BELOW THIS LINE !!!!!!!!!!!!!!!!!!!!!!!

  subroutine generate_global_to_local (MetaVertex, Global2Local, par)
    use dimensions_mod, only : nelem
    use metagraph_mod,  only : MetaVertex_t
    use parallel_mod,   only : parallel_t, MPI_MAX, MPIinteger_t
    !
    ! Inputs
    !
    type (parallel_t),  intent(in) :: par
    type(MetaVertex_t), intent(in) :: MetaVertex
    integer, intent(out) :: Global2Local(nelem)
    !
    ! Locals
    !
    integer :: ie, ierr

    ! Defaults all local ids to 0 (meaning not on this process)
    Global2Local = 0
    do ie=1,SIZE(MetaVertex%members)
      Global2Local(MetaVertex%members(ie)%number) = ie
    enddo

    call MPI_Allreduce(MPI_IN_PLACE,Global2Local,nelem,MPIinteger_t,MPI_MAX,par%comm,ierr)

  end subroutine generate_global_to_local

end module
