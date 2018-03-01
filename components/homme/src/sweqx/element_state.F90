#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module element_state

  use kinds,                  only: real_kind, long_kind, int_kind
  use dimensions_mod,         only: np, npsq, nlev, nlevp, qsize_d

  implicit none
  private
  integer, public, parameter :: timelevels = 3


  type, public :: elem_state_t

    ! prognostic variables for shallow-water solver
     real (kind=real_kind) :: p(np,np,nlev,timelevels)
     real (kind=real_kind) :: ps(np,np)                               ! surface geopotential
     real (kind=real_kind) :: gradps(np,np,2)                         ! gradient of surface geopotential
     real (kind=real_kind) :: v(np,np,2,nlev,timelevels)              ! contravarient comp

  end type elem_state_t

  !___________________________________________________________________
  type, public :: derived_state_t
  end type derived_state_t

  type, public :: elem_accum_t
  end type elem_accum_t


contains
end module 
