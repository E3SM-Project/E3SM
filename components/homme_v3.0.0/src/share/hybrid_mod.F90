#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

! ===========================================
! Module to support hybrid programming model
! hybrid_t is assumed to be a private struct
! ===========================================

module hybrid_mod
  use parallel_mod, only : parallel_t
implicit none
private

  type, public :: hybrid_t
     type (parallel_t) :: par
     integer           :: ithr
     integer           :: hthreads
     integer           :: vthreads
     logical           :: masterthread
  end type

  public :: hybrid_create

contains

  function hybrid_create(par,ithr,hthreads) result(hybrid)
      type (parallel_t), intent(in) :: par
      integer          , intent(in) :: ithr
      integer          , intent(in) :: hthreads
      type (hybrid_t)               :: hybrid

      hybrid%par      = par      ! relies on parallel_mod copy constructor
      hybrid%ithr     = ithr
      hybrid%hthreads = hthreads
      hybrid%masterthread = (par%masterproc .and. ithr==0)

  end function hybrid_create 

end module hybrid_mod
