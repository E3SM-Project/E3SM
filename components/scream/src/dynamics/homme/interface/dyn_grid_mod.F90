module dyn_grid_mod

  use kinds,              only: iulog
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use dimensions_mod,     only: nelem, nelemd, nelemdmax, np
  use homme_context_mod,  only: par, npes, iam, masterproc

  implicit none
  private

! We need MPI in here, so include it
#include <mpif.h>

  public :: dyn_grid_init, cleanup_grid_init_data

contains

  subroutine dyn_grid_init ()
    use prim_driver_base,     only: prim_init1_geometry, MetaVertex, GridEdge
    use prim_cxx_driver_base, only: init_cxx_connectivity
    use dimensions_mod,       only: nelemd
    use parallel_mod,         only: abortmp
    use homme_context_mod,    only: is_parallel_inited, elem, par, dom_mt

    if (.not. is_parallel_inited) then
      call abortmp ("Error! 'homme_init_parallel_f90' must be called *before* init_dyn_grid_f90.\n")
    endif

    if (masterproc) write(iulog,*) "Initializing dynamics grid..."

    call prim_init1_geometry(elem,par,dom_mt)

    call init_cxx_connectivity(nelemd,GridEdge,MetaVertex,par)

  end subroutine dyn_grid_init

  subroutine cleanup_grid_init_data ()
    use prim_driver_base, only: prim_init1_cleanup

    ! Cleanup the tmp stuff used in prim_init1_geometry
    call prim_init1_cleanup()
  end subroutine cleanup_grid_init_data

end module dyn_grid_mod
