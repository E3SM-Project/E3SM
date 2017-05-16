#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dimensions_mod
#ifdef CAM
  use constituents, only : qsize_d=>pcnst ! _EXTERNAL
#endif
  implicit none
  private

! set MAX number of tracers.  actual number of tracers is a run time argument  
#ifdef CAM
#else
#ifdef QSIZE_D
  integer, parameter         :: qsize_d=QSIZE_D
#else
  integer, parameter         :: qsize_d=4
#endif
#endif


  integer, parameter, public :: np = NP

  integer         :: qsize = 0

  integer, parameter, public :: npsq = np*np
  integer, parameter, public :: nlev=PLEV
  integer, parameter, public :: nlevp=nlev+1

!  params for a mesh 
!  integer, public, parameter :: max_elements_attached_to_node = 7
!  integer, public, parameter :: s_nv = 2*max_elements_attached_to_node 
  !default for non-refined mesh (note that these are *not* parameters now)
  integer, public  :: max_elements_attached_to_node = 4
  integer, public  :: s_nv = 6
  integer, public  :: max_corner_elem               = 1 !max_elements_attached_to_node-3
  integer, public  :: max_neigh_edges               = 8 !4 + 4*max_corner_elem

  public :: qsize,qsize_d

  integer, public :: ne
  integer, public :: nelem       ! total number of elements
  integer, public :: nelemd      ! number of elements per MPI task
  integer, public :: nelemdmax   ! max number of elements on any MPI task
  integer, public :: nnodes,npart,nmpi_per_node
  integer, public :: GlobalUniqueCols

  public :: set_mesh_dimensions

contains

  subroutine set_mesh_dimensions()

    ! new "params"
    max_elements_attached_to_node = 7  ! variable resolution
    s_nv = 2*max_elements_attached_to_node 

    !recalculate these
    max_corner_elem               = max_elements_attached_to_node-3
    max_neigh_edges               = 4 + 4*max_corner_elem

  end subroutine set_mesh_dimensions

end module dimensions_mod

