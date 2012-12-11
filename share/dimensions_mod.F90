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
  integer, parameter         :: qsize_d=4          ! SE tracers  
  integer, parameter         :: ntrac_d=4          ! fvm tracers
#endif
  integer, parameter, public :: nvar = 3 ! FI # dependent variables 

 
  integer, parameter, public :: np = NP
  integer, parameter, public :: nc = NC
  integer         :: ntrac = 0
  integer         :: qsize = 0

  ! fvm dimensions:
  integer, parameter, public :: ngpc=4    !number of Gausspoints for the fvm integral approximation  
  ! nhe ... halo/depth of the extended element (CFL number), now only nhe=1 is tested !
  !         this number sets where we have to calculate the reconstruction in the halo!
  integer, parameter, public :: nhe=1     !number/depth of the extended element (CFL number)
  ! nhc ... halo/depth for the tracer values, only cubic reconstruction is supported
  !         now, therefore nhc=nhe+3    
  integer, parameter, public :: nhc=nhe+3 
  
  ! constants for SPELT
  integer, parameter, public :: nip=3     !number of interpolation values, works only for this
  integer, parameter, public :: nipm=nip-1
  integer, parameter, public :: nep=nipm*nc+1      ! number of points in an element  
  
  
  integer, public :: npdg = 0  ! dg degree for hybrid cg/dg element  0=disabled 

  integer, parameter, public :: npsq = np*np
  integer, parameter, public :: nlev=PLEV
  integer, parameter, public :: nlevp=nlev+1


#ifdef MESH
  integer, public, parameter :: max_elements_attached_to_node = 7
  integer, public, parameter :: s_nv = 2*max_elements_attached_to_node 
#else
  integer, public,parameter :: max_elements_attached_to_node = 4
  integer, public, parameter :: s_nv = 6
#endif

  integer, public, parameter :: max_corner_elem               = max_elements_attached_to_node-3
  integer, public, parameter :: max_neigh_edges               = 4 + 4*max_corner_elem

  public :: qsize,qsize_d,ntrac_d,ntrac

  integer, public :: ne
  integer, public :: nelem       ! total number of elements
  integer, public :: nelemd      ! number of elements per MPI task
  integer, public :: nelemdmax   ! max number of elements on any MPI task
  integer, public :: nPhysProc                          ! This is the number of physics processors/ per dynamics processor
  integer, public :: nnodes,npart,nmpi_per_node
  integer, public :: GlobalUniqueCols



  public :: set_mesh_dimensions

contains

  subroutine set_mesh_dimensions()

    !This will be used when we removed the mesh compile flag

!    max_elements_attached_to_node = 7  ! variable resolution
!    s_nv = 2*max_elements_attached_to_node 

!    max_corner_elem               = max_elements_attached_to_node-3
!    max_neigh_edges               = 4 + 4*max_corner_elem


  end subroutine set_mesh_dimensions


end module dimensions_mod

