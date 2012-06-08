#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module physics_types_mod
  use kinds, only : real_kind, int_kind
  use dimensions_mod, only : np

  ! JPE: This parameter must match the number of variables in the state
  ! structure all of which are assumed to be of kind=real_kind.  This is a
  ! requirement for restart I/O.

  integer(kind=int_kind), parameter, public :: PhysStateComponents=1

  type, public :: physics_state_t
    sequence
    real  (kind=real_kind) :: CBMF(np,np)
  end type physics_state_t

  type, public :: physics_accum_t
    sequence
    real  (kind=real_kind) :: Precip(np,np)
  end type physics_accum_t

  type, public :: physics_surfc_t
    sequence
    real  (kind=real_kind) :: Wd(np,np)
    real  (kind=real_kind) :: Tprime(np,np)
    real  (kind=real_kind) :: Qprime(np,np)
    real  (kind=real_kind) :: Precip(np,np)
  end type physics_surfc_t

  type, public :: physics_t
    type(physics_state_t) :: state
    type(physics_surfc_t) :: surfc
    type(physics_accum_t) :: accum
  end type physics_t

  type(physics_t), public, allocatable, target :: pelem(:)

end module physics_types_mod
