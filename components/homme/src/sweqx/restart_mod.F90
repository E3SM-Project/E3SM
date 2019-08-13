#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module restart_mod 

   !------------------
   use kinds , only : int_kind, real_kind
   !------------------
   use dimensions_mod, only : nelemd
   !------------------
   use parallel_mod, only : parallel_t, MPIreal_t
   !------------------
   !use time_mod
   !------------------
   use element_state, only : elem_state_t
   !------------------
   !use control_mod
   !------------------
   use restart_io_mod, only : nwordsrestartbuffer_t, restartbuffer, File_elem_t, statedesc_t, &
        createstatedescriptor, ConstructElementFile, PrintStateDescriptor, addstatefield
   !------------------
   !use timer_mod
   !------------------
   implicit none

private 
    
   integer,parameter              :: RestartVersion=1
   ! =========================================
   !  Some variables used by all MPI routines
   ! =========================================
   integer                         :: errorcode,errorlen,ierr
   character(len=80)               :: errorstring
#if 0
#ifdef _MPI
   integer(kind=MPI_OFFSET_KIND)   :: offset,nbytes
#else
   integer(kind=int_kind)          :: offset,nbytes
#endif
#endif
   ! ====================================================
   !  Routines for Restart files
   ! ====================================================
   public :: initRestartFile

contains 
! =========================================================
! initRestartFile:
!
!  Initalizes MPI I-O  by seting up some MPI datastructures
! =========================================================
   subroutine initRestartFile(state,par,File)
    type (elem_state_t) :: state
    type (parallel_t),intent(in)    :: par
    type (File_elem_t),intent(out) :: File

    integer                      :: ie,ig,ierr

    integer                      :: count
    integer,allocatable          :: blklen(:),disp(:),oldtype(:)
    integer                      :: len
    type (StateDesc_t)           :: RestDesc
    integer                      :: NumComponents
    integer                      :: type
    integer                      :: i
 
    !=========================================
    !  Copy over the parallel_t datastructure
    !=========================================
    File%par = par

    !=========================================
    !  Allocate restart buffer
    !=========================================
    allocate(RestartBuffer(nelemd))

    !================================================================
    !   Construct the descriptor of the state variable for MPI I/O
    !================================================================
    NumComponents = 13     !  number of variable components in the state buffer

    RestDesc = CreateStateDescriptor(NumComponents)

    type = MPIReal_t       !  All the types are Real
    !=========================================
    ! Add all the fields in the State variable
    !=========================================

    len = SIZE(state%p)
    call AddStateField(RestDesc,len,type)

    len = SIZE(state%v)
    call AddStateField(RestDesc,len,type)

    len = SIZE(state%ps)
    call AddStateField(RestDesc,len,type)

    len = SIZE(state%gradps)
    call AddStateField(RestDesc,len,type)

#if defined(_MPI) && defined(_PRESTART)
    call PrintStateDescriptor(RestDesc)
    call ConstructElementFile(RestDesc,File,ierr)
#endif

    nwordsRestartBuffer_t=RestDesc%nwords

    print *,'initRestartFile: the count for CalcStateLength is ',nwordsRestartBuffer_t


    end subroutine initRestartFile

end module restart_mod
