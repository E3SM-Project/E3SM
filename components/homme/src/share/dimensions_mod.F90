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

  integer, parameter, public  :: max_elements_attached_to_node = 7 ! RRM meshes
! defaults for cubed sphere grids:
  integer, public  :: max_corner_elem               = 1  !  max_elements_attached_to_node-3
  integer, public  :: max_neigh_edges               = 8  !  4 + 4*max_corner_elem

  public :: qsize,qsize_d

  integer, public :: ne
  integer, public :: ne_x,ne_y   ! used for planar topology- number of elements in each direction
  integer, public :: nelem       ! total number of elements
  integer, public :: nelemd      ! number of elements per MPI task
  integer, public :: nelemdmax   ! max number of elements on any MPI task
  integer, public :: nnodes,npart,nmpi_per_node
  integer, public :: GlobalUniqueCols

  public :: set_mesh_dimensions

contains

  subroutine set_mesh_dimensions()

    !recalculate these for RRM meshes
    max_corner_elem               = max_elements_attached_to_node-3
    max_neigh_edges               = 4 + 4*max_corner_elem

  end subroutine set_mesh_dimensions
end module dimensions_mod
