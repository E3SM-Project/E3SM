module elm_finalizeMod

  !-----------------------------------------------------------------------
  ! Performs land model cleanup
  !
  !
  implicit none
  save
  public   ! By default everything is public
  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------
  !
  public :: final
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine final( )
    !
    ! !DESCRIPTION:
    ! Finalize land surface model
    !
#ifdef USE_PETSC_LIB
#include <petsc/finclude/petsc.h>
#endif
    ! !USES:
#ifdef USE_PETSC_LIB
    use petscsys
#endif
    !
    ! !ARGUMENTS
    implicit none
    !

#ifdef USE_PETSC_LIB
    PetscErrorCode        :: ierr

#endif

  end subroutine final

end module elm_finalizeMod
