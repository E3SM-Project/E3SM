#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module types_mod
  
  use kinds, only : real_kind

  integer, parameter, public :: RKMaxStages = 100

  type, public :: Rk_t
     integer                                          :: Stages
     real (kind=real_kind) :: alpha0(RKMaxStages)
     real (kind=real_kind) :: alpha(RKMaxStages)
     real (kind=real_kind) :: beta(RKMaxStages)
     ! CFL
     real (kind=real_kind)                            :: RKCFL
  end type

end module types_mod
