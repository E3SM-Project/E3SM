#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module moabmesh_mod    
!  use, intrinsic :: ISO_C_BINDING

  use kinds, only : real_kind, int_kind, longdouble_kind, iulog
!  use edge_mod, only : ghostbuffertr_t, initghostbufferTR, freeghostbuffertr, &
!       ghostVpack, ghostVunpack,  edgebuffer_t, initEdgebuffer

  use dimensions_mod, only: nelem, nelemd, nelemdmax, nlev, ne, nc, nhe, nlev, ntrac, np, ntrac_d,ns, nhr, nhc
!  use time_mod, only : timelevel_t
  use element_mod, only : element_t, timelevels
!  use fvm_control_volume_mod, only: fvm_struct
!  use fvm_mod, only: cellghostbuf, edgeveloc
  use hybrid_mod, only : hybrid_t
  use parallel_mod, only : parallel_t
  use coordinate_systems_mod, only: spherical_polar_t, cartesian3D_t, spherical_to_cart

!  use control_mod, only: tracer_transport_type, TRACERTRANSPORT_MOAB , TRACERTRANSPORT_MOAB_LINEAR

  implicit none

#include "moab/MOABConfig.h"
!#include "moab/iMOAB.h"

!  public :: moabmesh

contains

  SUBROUTINE errorout(ierr, message)
  integer ierr
  character*(*) message
  if (ierr.ne.0) then
    print *, message
    call exit (1)
  end if
  return
  end subroutine

  subroutine create_moab_mesh(hybrid, elem, nets, nete)

    use ISO_C_BINDING

    type (element_t), intent(inout) :: elem(:)

    type (hybrid_t)      , intent(in) :: hybrid

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete

    integer ierr, i, j, ie, pid
    integer filenum

    real(kind=real_kind), allocatable, target :: moab_vert_coords(:)
    INTEGER (C_INT) ,allocatable , target :: moab_corner_quads(:)
    integer moab_dims_fvc, moab_dim_cquads, ix, stride, idx ! used for indexing in loops
    integer moab_dims_cell
    type(cartesian3D_t) :: tmppt ! used to convert to cartesian from spherical
    integer nelemd
    character*10 appname
    integer , external :: iMOAB_InitializeFortran, iMOAB_RegisterFortranApplication


 !  hybrid%par

     nelemd = nete-nets+1
 !    moab_dims_fvc = (nc+1)*(nc+1)*(nete-nets+1)*3
     moab_dim_cquads = (nete-nets+1)*4

!     allocate(moab_fvm_coords(moab_dims_fvc))
     allocate(moab_corner_quads(moab_dim_cquads))
     do ie=nets,nete
       ix = ie-nets
!       do j=1,4
!         moab_corner_quads(ix*4+j) = elem(ie)%node_numbers(j)
!       enddo
!      write(iulog,*)( elem(ie)%node_numbers(j), j=1,4)
       moab_corner_quads(ix*4+1) = elem(ie)%gdofP(1,1)
       moab_corner_quads(ix*4+2) = elem(ie)%gdofP(np,1)
       moab_corner_quads(ix*4+3) = elem(ie)%gdofP(np,np)
       moab_corner_quads(ix*4+4) = elem(ie)%gdofP(1,np)
       ! some coordinates?
     enddo


!     call create_mesh(moabdata%imeshInstance, moabdata%fineSet, moabdata%depSet, &
!        moabdata%intxSet, c_loc(moab_fvm_coords ) , &
!        c_loc(moab_corner_quads), nc, nelemd, tracer_transport_type, moabdata%par%comm, ierr)

      ierr = iMOAB_InitializeFortran()
      call errorout(ierr, 'fail to initialize iMOAB')

      appname = 'HM_COARSE'

      ierr = iMOAB_RegisterFortranApplication(appname, hybrid%par, pid)
      call errorout(ierr, 'fail to register homme spectral mesh')

     ! deallocate
     deallocate(moab_corner_quads)

  end subroutine create_moab_mesh
  

end module moabmesh_mod
