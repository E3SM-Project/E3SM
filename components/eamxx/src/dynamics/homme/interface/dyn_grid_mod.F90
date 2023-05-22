module dyn_grid_mod

  use kinds,              only: iulog
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use dimensions_mod,     only: nelem, nelemd, nelemdmax, np
  use edgetype_mod,       only: EdgeBuffer_t

  implicit none
  private

! We need MPI in here, so include it
#include <mpif.h>

  public :: dyn_grid_init, get_my_dyn_data, cleanup_grid_init_data

  type (EdgeBuffer_t) :: edge

contains

  subroutine dyn_grid_init ()
    use prim_driver_base,     only: prim_init1_geometry, MetaVertex, GridEdge
    use prim_cxx_driver_base, only: init_cxx_connectivity
    use dimensions_mod,       only: nelemd
    use parallel_mod,         only: abortmp
    use homme_context_mod,    only: is_parallel_inited, elem, par, dom_mt, masterproc
    use edge_mod_base,        only: initEdgeBuffer

    if (.not. is_parallel_inited) then
      call abortmp ("Error! 'homme_init_parallel_f90' must be called *before* init_dyn_grid_f90.\n")
    endif

    if (masterproc) then
      write(iulog,*) "Initializing dynamics grid..."
    endif

    ! Initialize geometric structures in homme
    call prim_init1_geometry(elem,par,dom_mt)
    call init_cxx_connectivity(nelemd,GridEdge,MetaVertex,par)

    ! Init internal edge buffer, used in get_my_dyn_data
    call initEdgeBuffer(par,edge,elem,1)

  end subroutine dyn_grid_init

  subroutine get_my_dyn_data (dg_gids, cg_gids, elgpgp, lat, lon)
    use iso_c_binding,     only: c_int, c_double
    use dimensions_mod,    only: nelemd, np
    use homme_context_mod, only: elem, par
    use shr_const_mod,     only: pi=>SHR_CONST_PI
    use bndry_mod_base,    only: bndry_exchangeV
    use edge_mod_base,     only: edgeVpack_nlyr, edgeVunpack_nlyr
    use kinds,             only: real_kind, int_kind
    use dof_mod,           only: genLocalDof
    !
    ! Inputs
    !
    real(kind=c_double), intent(out) :: lat (:,:,:), lon(:,:,:)
    integer(kind=c_int), intent(out) :: cg_gids (:), dg_gids(:), elgpgp(:,:)
    !
    ! Local(s)
    !
    real(kind=real_kind), allocatable :: el_cg_gids (:,:,:)  ! Homme's bex stuff only works with reals
    integer(kind=int_kind), allocatable :: el_dg_gids (:,:,:)  ! Homme's getLocalDof might not work with c_int
    integer :: idof, ip,jp, ie, icol

    ! Get the gids
    allocate(el_cg_gids(np,np,nelemd))
    allocate(el_dg_gids(np,np,nelemd))
    el_dg_gids = 0
    el_cg_gids = 0
    do ie=1,nelemd
      call genLocalDof(elem(ie)%vertex%number,np,el_dg_gids(:,:,ie))
      do icol=1,elem(ie)%idxP%NumUniquePts
        ip = elem(ie)%idxP%ia(icol)
        jp = elem(ie)%idxP%ja(icol)
        el_cg_gids(ip,jp,ie) = elem(ie)%idxP%UniquePtOffset + icol - 1
      enddo
      call edgeVpack_nlyr(edge,elem(ie)%desc,el_cg_gids(:,:,ie),1,0,1)
    enddo
    call bndry_exchangeV(par,edge)
    do ie=1,nelemd
      call edgeVunpack_nlyr(edge,elem(ie)%desc,el_cg_gids(:,:,ie),1,0,1)
      do ip=1,np
        do jp=1,np
          idof = (ie-1)*16+(jp-1)*4+ip
          cg_gids(idof) = INT(el_cg_gids(ip,jp,ie),kind=c_int)
          dg_gids(idof) = INT(el_dg_gids(ip,jp,ie),kind=c_int)
          lat(ip,jp,ie)  = elem(ie)%spherep(ip,jp)%lat * 180.0_c_double/pi
          lon(ip,jp,ie)  = elem(ie)%spherep(ip,jp)%lon * 180.0_c_double/pi
          elgpgp(1,idof) = ie-1
          elgpgp(2,idof) = jp-1
          elgpgp(3,idof) = ip-1
        enddo
      enddo
    enddo
  end subroutine get_my_dyn_data

  subroutine cleanup_grid_init_data ()
    use edge_mod_base,    only: FreeEdgeBuffer

    ! Cleanup edge used in get_my_dyn_data
    ! NOTE: do not call prim_init1_cleanup, since compose,
    ! which has not been inited yet, needs some grid connectivity info
    call FreeEdgeBuffer(edge)
  end subroutine cleanup_grid_init_data

end module dyn_grid_mod
