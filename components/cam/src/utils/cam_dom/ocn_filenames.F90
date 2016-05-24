module ocn_filenames

!----------------------------------------------------------------------- 
!
! DESCRIPTION
! Module and methods to handle filenames needed for the model. This 
! includes input filenames, and most output filenames that the model
! uses. All filenames that the model uses will use methods or data
! constructed by this module. In some cases (such as the cam_history module)
! other modules or routines will store the actual filenames used, but
! this module is used to determine the names.
!
!----------------------------------------------------------------------- 

   use shr_kind_mod,     only: shr_kind_cs, shr_kind_cl
   use shr_sys_mod,      only: shr_sys_abort
   use filenames,        only: caseid    !CAM specific use

   implicit none

   integer, parameter :: nlen = shr_kind_cl                ! String length

CONTAINS

  character(len=nlen) function restart_filename( year, month, day, sec, case )

  !----------------------------------------------------------------------- 
  !
  ! !DESCRIPTION: Create a filename from a filename specifyer. The 
  ! filename specifyer includes codes for setting things such as the
  ! year, month, day, seconds in day, caseid, and tape number. This
  ! routine is private to filenames.F90
  !
  ! Interpret filename specifyer string with: 
  !
  !      %c for case, 
  !      %t for optional number argument sent into function
  !      %y for year
  !      %m for month
  !      %d for day
  !      %s for second
  !      %% for the "%" character
  !
  ! If the filename specifyer has spaces " ", they will be trimmed out
  ! of the resulting filename.
  !
  !----------------------------------------------------------------------- 
  !
  ! Arguments
  !
  integer         , intent(in) :: year           ! Simulation year
  integer         , intent(in) :: month          ! Simulation month
  integer         , intent(in) :: day            ! Simulation day
  integer         , intent(in) :: sec            ! Seconds into current simulation day
  character(len=*), intent(in), optional :: case ! Optional casename
  !
  ! Local variables
  !
  character(len=nlen) :: string    ! Temporary character string 
  character(len=nlen) :: format    ! Format character string 
  integer :: i, n                  ! Loop variables
  logical :: done
  character(len=nlen) :: filename_spec = '%c.camdom.r.%y-%m-%d-%s' ! ocn restarts
  !----------------------------------------------------------------------- 

  if ( len_trim(filename_spec) == 0 )then
     call shr_sys_abort ('RESTART_FILENAME: filename specifier is empty')
  end if
  if ( index(trim(filename_spec)," ") /= 0 )then
     call shr_sys_abort ('RESTART_FILENAME: filename specifier can not contain a space:'//trim(filename_spec))
  end if
  !
  ! Go through each character in the filename specifyer and interpret if special string
  !
  i = 1
  restart_filename = ''
  do while ( i <= len_trim(filename_spec) )
     !
     ! If following is an expansion string
     !
     if ( filename_spec(i:i) == "%" )then
        i = i + 1
        select case( filename_spec(i:i) )
        case( 'c' )   ! caseid
           if ( present(case) )then
              string = trim(case)
           else
              string = trim(caseid)
           end if
        case( 'y' )   ! year
           if ( year > 99999   ) then
              format = '(i6.6)'
           else if ( year > 9999    ) then
              format = '(i5.5)'
           else
              format = '(i4.4)'
           end if
           write(string,format) year
        case( 'm' )   ! month
           write(string,'(i2.2)') month
        case( 'd' )   ! day
           write(string,'(i2.2)') day
        case( 's' )   ! second
           write(string,'(i5.5)') sec
        case( '%' )   ! percent character
           string = "%"
        case default
           call shr_sys_abort ('RESTART_FILENAME: Invalid expansion character: '//filename_spec(i:i))
        end select
        !
        ! Otherwise take normal text up to the next "%" character
        !
     else
        n = index( filename_spec(i:), "%" )
        if ( n == 0 ) n = len_trim( filename_spec(i:) ) + 1
        if ( n == 0 ) exit 
        string = filename_spec(i:n+i-2)
        i = n + i - 2
     end if
     if ( len_trim(restart_filename) == 0 )then
        restart_filename = trim(string)
     else
        if ( (len_trim(restart_filename)+len_trim(string)) >= nlen )then
           call shr_sys_abort ('RESTART_FILENAME: Resultant filename too long')
        end if
        restart_filename = trim(restart_filename) // trim(string)
     end if
     i = i + 1
  end do
  if ( len_trim(restart_filename) == 0 )then
     call shr_sys_abort('RESTART_FILENAME: Resulting filename is empty')
  end if

end function restart_filename

end module ocn_filenames
