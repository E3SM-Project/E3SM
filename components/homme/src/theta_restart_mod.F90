#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_restart_mod 
   !------------------
   use kinds , only : real_kind, int_kind
   !------------------
   use dimensions_mod, only : nelemd
   !------------------
   use parallel_mod, only : parallel_t, MPIreal_t, abortmp
   !------------------
   use element_state, only : elem_state_t
   !------------------
   use restart_io_mod, only : nwordsRestartBuffer_t, RestartBuffer,  File_elem_t, &
        StateDesc_t, createstatedescriptor, AddStateField, constructelementfile, &
        collective_io_read, collective_io_write, printstatedescriptor
   !------------------
   implicit none

private 
   ! =========================================
   !  Some variables used by all MPI routines
   ! =========================================
   integer                         :: errorcode,errorlen,ierr
   character(len=80)               :: errorstring
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
   subroutine initRestartFile(state, par,File)
    type (elem_state_t), intent(in) :: state
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
    logical,parameter            :: Debug = .FALSE.

    !=========================================
    !  Copy over the parallel_t datastructure
    !=========================================
    File%par = par

    !=========================================
    !  Allocate restart buffer
    !=========================================
    if(Debug) print *,'initRestartFile: nelemd: ',nelemd

    collective_io_read = .true.
    collective_io_write = .true.

    allocate(RestartBuffer(nelemd))

    !================================================================
    !   Construct the descriptor of the state variable for MPI I/O
    !================================================================
    NumComponents = 9   ! THIS NUMBER MUST MATCH NUMBER OF TIMES AddStateField is called below
    RestDesc = CreateStateDescriptor(NumComponents)

    type = MPIReal_t       !  All the types are Real
    !=========================================
    ! for PRESTART, must add *all* the fields in the State variable
    !=========================================
    len = SIZE(state%v)
    call AddStateField(RestDesc,len,type)
#ifdef MODEL_THETA_L
    len = SIZE(state%w_i)
#else
    len = SIZE(state%w)
#endif
    call AddStateField(RestDesc,len,type)
#ifdef MODEL_THETA_L
    len = SIZE(state%vtheta_dp)
#else
    len = SIZE(state%theta_dp_cp)
#endif
    call AddStateField(RestDesc,len,type)

    len = SIZE(state%ps_v)
    call AddStateField(RestDesc,len,type)

#ifdef MODEL_THETA_L
    len = SIZE(state%phinh_i)
#else
    len = SIZE(state%phinh)
#endif
    call AddStateField(RestDesc,len,type)

    len = SIZE(state%phis)
    call AddStateField(RestDesc,len,type)

    len = SIZE(state%Q)
    call AddStateField(RestDesc,len,type)
    len = SIZE(state%Qdp)
    call AddStateField(RestDesc,len,type)

    len = SIZE(state%dp3d)
    call AddStateField(RestDesc,len,type)

#if defined(_MPI) && defined(_PRESTART)
    if(Debug) call PrintStateDescriptor(RestDesc)
    call ConstructElementFile(RestDesc,File,ierr)
#endif
    nwordsRestartBuffer_t=RestDesc%nwords

    end subroutine initRestartFile
end module prim_restart_mod
