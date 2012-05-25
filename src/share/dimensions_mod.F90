#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dimensions_mod
#ifdef CAM
  use constituents, only : qsize_d=>pcnst ! _EXTERNAL
  use constituents, only : ntrac_d=>pcnst ! _EXTERNAL
#endif
  implicit none
  private

! set MAX number of tracers.  actual number of tracers is a run time argument  
#ifndef CAM
  integer, parameter         :: qsize_d=4           ! SE tracers  
  integer, parameter         :: ntrac_d=4          ! CSLAM tracers
#endif
  integer, parameter, public :: nvar = 3 ! FI # dependent variables 


#ifdef _CSLAM 
  ! settings used by 2D CSLAM test code: 
  integer, parameter, public :: np = NC+1 ! this can be delted once CSLAM has it own I/O
  ! nc  : nc*nc is the number of finite volume cells on each element (we have 
  !       6*ne*ne*nc*nc cells on the sphere),   
  integer, parameter, public :: nc=NC       

  ! set the number of tracers
  integer         :: ntrac =NTRAC
  integer         :: qsize=0
#else   
  integer, parameter, public :: np = NP
  integer, parameter, public :: nc = 4
  integer         :: ntrac = 0
  integer         :: qsize=qsize_d
#endif

  ! CSLAM dimensions:
  integer, parameter, public :: ngpc=4    !number of Gausspoints for the CSLAM integral approximation  
  ! nhe ... halo/depth of the extended element (CFL number), now only nhe=1 is tested !
  !         this number sets where we have to calculate the reconstruction in the halo!
  integer, parameter, public :: nhe=1     !number/depth of the extended element (CFL number)
  ! nhc ... halo/depth for the tracer values, only cubic reconstruction is supported
  !         now, therefore nhc=nhe+3    
  integer, parameter, public :: nhc=nhe+3 
  integer, public :: npdg = 0  ! dg degree for hybrid cg/dg element  0=disabled 

  integer, parameter, public :: npsq = np*np
  integer, parameter, public :: nlev=PLEV
  integer, parameter, public :: nlevp=nlev+1
#ifdef MESH
  integer, parameter, public :: max_elements_attached_to_node = 7  ! variable resolution
#else
  integer, parameter, public :: max_elements_attached_to_node = 4
#endif
  integer, parameter, public :: max_corner_elem               = max_elements_attached_to_node-3
  integer, parameter, public :: max_neigh_edges               = 4 + 4*max_corner_elem
  public :: qsize,qsize_d,ntrac_d,ntrac

  integer, public :: ne
  integer, public :: nelem       ! total number of elements
  integer, public :: nelemd      ! number of elements per MPI task
  integer, public :: nelemdmax   ! max number of elements on any MPI task
  integer, public :: nPhysProc                          ! This is the number of physics processors/ per dynamics processor
  integer, public :: nnodes,npart,nmpi_per_node
  integer, public :: GlobalUniqueCols

end module dimensions_mod

