! Trimmed-down version of glc_communicate including just what is needed for cism unit tests,
! in order to avoid dependencies

module glc_communicate

   use glc_kinds_mod
   use shr_sys_mod, only : shr_sys_abort

   implicit none
   public
   save

   integer(int_kind), parameter :: my_task = 0
   integer(int_kind), parameter :: master_task = 0

contains
   
   subroutine exit_message_environment(ierr)
      integer (int_kind), intent(out) :: ierr

      return
   end subroutine exit_message_environment

   subroutine abort_message_environment(ierr)
      integer (int_kind), intent(out) :: ierr

      call shr_sys_abort('glc_communicate.F90: abort_message_environment')
   end subroutine abort_message_environment
      
end module glc_communicate
