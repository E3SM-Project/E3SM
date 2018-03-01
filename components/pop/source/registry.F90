!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module registry 

!BOP
! !MODULE: registry
!
! !DESCRIPTION:
!  This module provides a means for registering, checking, and
!     recording events that have occurred in CCSM POP
!
! !REVISION HISTORY:
!  SVN:$Id: registry.F90 12674 2008-10-31 22:21:32Z njn01 $

! !USES:

   use kinds_mod
   use exit_mod
   use io_tools
 
   implicit none
   private
   save
 
! !PUBLIC MEMBER FUNCTIONS:
   public ::                 &
      init_registry,         &
      registry_match,        &
      register_string,       &
      registry_err_check,    &
      trap_registry_failure

!EOP
!BOC

   integer (int_kind), parameter ::  &
      max_registry_size = 200        ! maximum size of registry
 
   integer (int_kind) ::  &
      registry_failure_count, &
      registry_size
 
   character (char_len), dimension (max_registry_size) ::  &
      registry_storage
 
!EOC
!***********************************************************************

   contains

!***********************************************************************
 subroutine init_registry 

!-----------------------------------------------------------------------
!
!  This routine initializes the registry storage array and
!  sets the failure counter to zero
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      n                     ! dummy loop index
 
 
   call reset_registry_failure_count
   registry_size  = 0

   do n=1,max_registry_size
     registry_storage(n) = ' '
   end do 
 
 end subroutine init_registry

 
 function registry_match (string)

!-----------------------------------------------------------------------
!  This function checks to see if a string has already been registered
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     input variables
!-----------------------------------------------------------------------
 
   character (*), intent(in) :: string    
 
!-----------------------------------------------------------------------
!     output variables
!-----------------------------------------------------------------------
 
   logical (log_kind) :: registry_match    !  T ==> string is registered    
 
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------

   integer (int_kind) ::  n                 ! dummy loop index

!-----------------------------------------------------------------------
!
!     search to determine if string has already been registered
!
!-----------------------------------------------------------------------

   registry_match = .false.
 
   string_search: do n=1,max_registry_size
     if ( registry_storage(n)  == string) then
       registry_match = .true.
       exit string_search
     endif
   end do string_search
    
 end function registry_match  

 
 subroutine reset_registry_failure_count
     registry_failure_count = 0
 end subroutine reset_registry_failure_count
 
 
 subroutine register_string (string)

!-----------------------------------------------------------------------
!     this routine registers a character string in registry_storage
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     input variables
!-----------------------------------------------------------------------
 
   character (*), intent(in) :: string    ! string added to registry
 

    
!-----------------------------------------------------------------------
!     if string is not already defined, add string to registry
!-----------------------------------------------------------------------
 
   if (.not. registry_match(string) ) then   
     registry_size = registry_size + 1
     registry_storage(registry_size) = string
   endif
 
 end subroutine register_string
 
 
 subroutine registry_err_check (string,string_present,caller)

!-----------------------------------------------------------------------
!    This routine complains if a string is in the registry but
!    should not be, or is not in the registry but should be.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!   input variables
!-----------------------------------------------------------------------
 
   character (*), intent(in) ::   &
      string,                     & ! test string
      caller                        ! calling routine name
 
   logical (log_kind),intent(in) ::  &
      string_present                  ! T ==> want string to be IN the registry

!-----------------------------------------------------------------------
!   local variables
!-----------------------------------------------------------------------
    
   character (char_len) :: message   ! error message
 
 
!-----------------------------------------------------------------------
!   check for error conditions; if error exits, print message and
!   increment registry_failure_count
!-----------------------------------------------------------------------
 
 
   if ((registry_match(string) .neqv. string_present)) then
 
      if (string_present) then
        write(message,1100) 'registry_error:', trim(string), &
             'has NOT been registered -- calling routine is ', &
              trim(caller)
      else
        write(message,1100) 'registry_error:', trim(string), &
             'has ALREADY been registered -- calling routine is ', &
              trim(caller)
      endif
    
      1100 format(1x, 4a)
 
      call document('registry_err_check',message)
 
      registry_failure_count = registry_failure_count + 1
   else
 
   endif ! registry_match
    
 end subroutine registry_err_check

 
 subroutine trap_registry_failure
!-----------------------------------------------------------------------
!
!    This subroutine checks to see if there have been any registry
!    failures.  If any have occurred, then the model will stop.
!-----------------------------------------------------------------------

   if (registry_failure_count /= 0) then
    call exit_POP (sigAbort, &
      'Registry failure count > 0 ; search output for "registry_error" for info')  
   endif
 
 end subroutine trap_registry_failure
 
 
 end module registry 

