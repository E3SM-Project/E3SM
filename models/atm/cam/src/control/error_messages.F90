module error_messages

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! General purpose routines for issuing error messages.
   ! 
   ! Author: B. Eaton
   ! 
   !----------------------------------------------------------------------- 
   use abortutils,  only: endrun
   use cam_logfile, only: iulog

   implicit none
   save
   private
   public :: &
      alloc_err,      &! Issue error message after non-zero return from an allocate statement.
      handle_err,     &! Issue error message after non-zero return from anything
      handle_ncerr     ! Handle error returns from netCDF library procedures.

   ! If an error message string is not empty, abort with that string as the
   ! error message.
   public :: handle_errmsg

!##############################################################################
contains
!##############################################################################

   subroutine alloc_err( istat, routine, name, nelem )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Issue error message after non-zero return from an allocate statement.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      integer, intent(in) ::&
         istat           ! status from allocate statement
      character(len=*), intent(in) ::&
         routine,       &! routine that called allocate
         name            ! name of array
      integer, intent(in) ::&
         nelem           ! number of elements attempted to allocate
      !-----------------------------------------------------------------------

      if ( istat .ne. 0 ) then
         write(iulog,*)'ERROR trying to allocate memory in routine: ' &
                   //trim(routine)
         write(iulog,*)'  Variable name: '//trim(name)
         write(iulog,*)'  Number of elements: ',nelem
         call endrun ('ALLOC_ERR')
      end if

      return

   end subroutine alloc_err

!##############################################################################

   subroutine handle_err( istat, msg )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Issue error message after non-zero return from anything.
      !
      ! Author: T. Henderson
      !----------------------------------------------------------------------- 

      integer,          intent(in) :: istat  ! status, zero = "no error"
      character(len=*), intent(in) :: msg    ! error message to print
      !-----------------------------------------------------------------------

      if ( istat .ne. 0 ) then
         call endrun (trim(msg))
      end if

      return

   end subroutine handle_err

!##############################################################################

   subroutine handle_ncerr( ret, mes, line )
      
      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Check netCDF library function return code.  If error detected 
      ! issue error message then abort.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

!-----------------------------------------------------------------------
     use netcdf
!-----------------------------------------------------------------------

      integer, intent(in) ::&
         ret                 ! return code from netCDF library routine
      character(len=*), intent(in) ::&
         mes                 ! message to be printed if error detected
      integer, intent(in), optional :: line
      !-----------------------------------------------------------------------

      if ( ret .ne. NF90_NOERR ) then
         if(present(line)) then
            write(iulog,*) mes, line
         else	
            write(iulog,*) mes
         end if
         write(iulog,*) nf90_strerror( ret )
         call endrun ('HANDLE_NCERR')
      endif

      return

   end subroutine handle_ncerr

!##############################################################################

   subroutine handle_errmsg(errmsg, subname, extra_msg)

     ! String that is asserted to be null.
     character(len=*), intent(in)           :: errmsg
     ! Name of procedure generating the message.
     character(len=*), intent(in), optional :: subname
     ! Additional message from the procedure calling this one.
     character(len=*), intent(in), optional :: extra_msg

     if (trim(errmsg) /= "") then

        if (present(extra_msg)) &
             write(iulog,*) "handle_errmsg: &
             &Message from caller: ",trim(extra_msg)

        if (present(subname)) then
           call endrun("ERROR: handle_errmsg: "// &
                trim(subname)//": "//trim(errmsg))
        else
           call endrun("ERROR: handle_errmsg: "// &
                "Error message received from routine: "//trim(errmsg))
        end if

     end if

   end subroutine handle_errmsg

!##############################################################################

end module error_messages
