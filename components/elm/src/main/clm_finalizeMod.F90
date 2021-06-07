module clm_finalizeMod

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
    use clm_varctl             , only : use_vsfm
#ifdef USE_PETSC_LIB
    use petscsys
#endif
    !
    ! !ARGUMENTS
    implicit none
    !

#ifdef USE_PETSC_LIB
    PetscErrorCode        :: ierr

    if (use_vsfm) then
       call PetscFinalize(ierr)
    endif
#endif

  end subroutine final

end module clm_finalizeMod
