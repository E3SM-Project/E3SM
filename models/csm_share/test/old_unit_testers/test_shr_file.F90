program test_shr_file
use shr_sys_mod,  only: shr_sys_abort, shr_sys_system
use shr_file_mod, only: shr_file_getUnit, shr_file_freeUnit, &
                        shr_file_chDir, shr_file_chStdIn, shr_file_chStdOut
!
! unit test of the shr_file_mod module
! 
write(6,*) 'Test file get/put: '
call test_getput()

write(6,*) 'Test units: '
call test_unit()

! Test the stdio series of subroutines
write(6,*) 'Test stdio: '
call test_stdio()

stop "Tests Pass"

contains

subroutine test_stdio()
use shr_sys_mod, only: shr_sys_getenv, shr_sys_chdir
integer, parameter :: nModels = 5
character(len=3), parameter :: models(nmodels) = (/"atm", "lnd", "ice", "ocn", "cpl"/)
character(len=256) :: nlfile
character(len=256) :: pwd, cwd
character(len=256), parameter :: logfile = "test_shr_file.log"
integer :: i, unit, j
integer :: rcode
logical :: exists
namelist /atm_inparm/ j
namelist /lnd_inparm/ j
namelist /ocn_inparm/ j
namelist /ice_inparm/ j
namelist /cpl_inparm/ j

call shr_sys_getenv( "pwd", pwd, rcode )
call shr_sys_system( "/bin/rm "//trim(logfile), rcode )
do i = 1, nModels
   call shr_sys_system( "/bin/rm "//models(i)//".log", rcode )
   if ( i == 1 )then
     open(6,file=logfile,status="new")
   else
     open(6,file=logfile,status="old", position="append")
   end if
   write(6,*) "test model: ", models(i)
   write(6,*) "test chdir: "
   call shr_file_chDir(models(i),rcodeOut=rcode)
   if ( rcode /= 0 )then
      call shr_sys_abort( "error: chDir returns error code" )
   end if
   call shr_sys_getenv( "pwd", cwd, rcode )
   !if ( trim(pwd)//"/nl" /= cwd )then
   !   write(6,*) 'pwd = ', trim(pwd)
   !   write(6,*) 'cwd = ', trim(cwd)
   !   call shr_sys_abort( "error: chDir did not go to correct directory" )
   !end if
   write(6,*) "test chstdin: "
   if ( (models(i) == "atm") .or. (models(i) == "lnd") )then
      call shr_file_chStdIn(models(i), NLFilename=nlfile,rcodeOut=rcode)
      unit = shr_file_getUnit()
      inquire(file=nlfile,exist=exists)
      if ( .not. exists )then
         call shr_sys_abort( "error: nlfilename does NOT exist: "//trim(nlfile) )
      end if
      open(unit,file=trim(nlfile),status="old")
   else
      call shr_file_chStdIn(models(i),rcodeOut=rcode)
      unit = 5
   end if
   if ( rcode /= 0 )then
      call shr_sys_abort( "error: chstdin returns error code" )
   end if
   if (      models(i) == "atm" )then
      read(unit,nml=atm_inparm,iostat=rcode)
   else if ( models(i) == "lnd" )then
      read(unit,nml=lnd_inparm,iostat=rcode)
   else if ( models(i) == "ocn" )then
      read(unit,nml=ocn_inparm,iostat=rcode)
   else if ( models(i) == "ice" )then
      read(unit,nml=ice_inparm,iostat=rcode)
   else if ( models(i) == "cpl" )then
      read(unit,nml=cpl_inparm,iostat=rcode)
   end if
   close(unit)
   if ( rcode /= 0 )then
      call shr_sys_abort( "error: reading namelist returns error code" )
   end if
   write(6,*) "test chstdout: "
   call shr_file_chStdOut(models(i),rcodeOut=rcode)
   if ( rcode /= 0 )then
      call shr_sys_abort( "error: chstdout returns error code" )
   end if
   write(6,*) "<<<<<<<<This should go to the "//models(i)//" stdout file>>>>>>>>>>"
   call shr_sys_chdir("..",rcode)
   close(6)
end do

end subroutine test_stdio

subroutine is_prefix( filename, expPrefix, nExpPrefix )
use shr_file_mod, only: shr_file_queryPrefix, shr_file_noPrefix
character(*), intent(IN) :: filename
character(*), intent(IN) :: ExpPrefix
integer,      intent(IN) :: nExpPrefix

integer :: nPrefix
character(256) :: Prefix

nPrefix = shr_file_queryPrefix( filename, prefix=prefix )
if ( nPrefix /= nExpPrefix .or. trim(prefix) /= trim(ExpPrefix) )then
   write(6,*) 'Prefix = ', trim(prefix), 'Expected = ', trim(ExpPrefix), " End"
   write(6,*) 'N-Prefix = ', nPrefix, 'N-Expected = ', nExpPrefix
   call shr_sys_abort( "error: wrong prefix type or wrong returned prefix length" )
end if

end subroutine is_prefix

subroutine test_getput()
use shr_file_mod, only: shr_file_queryPrefix, shr_file_get, shr_file_put, shr_file_noPrefix, &
                        shr_file_nullPrefix, shr_file_cpPrefix, shr_file_mssPrefix, &
                        shr_file_hpssPrefix
character(256) :: filename
character(256) :: prefix
integer :: nprefix


filename = "/long:directory_d/sub-directory::/file:with_colon.txt"
call is_prefix( filename, "", shr_file_noPrefix )
filename = "cp:/longdirectory_d/sub-directory::/file:with_colon.txt"
call is_prefix( filename, "cp:", shr_file_cpPrefix )
filename = "null:/long:directory_d/sub-directory::/file:with_colon.txt"
call is_prefix( filename, "null:", shr_file_nullPrefix )
filename = "mss:/long:directory_d/sub-directory::/file:with_colon.txt"
call is_prefix( filename, "mss:", shr_file_mssPrefix )
filename = "hpss:/long:directory_d/sub-directory::/file:with_colon.txt"
call is_prefix( filename, "hpss:", shr_file_hpssPrefix )
filename = "file:with_colon.txt"
call is_prefix( filename, "", shr_file_noPrefix )

end subroutine test_getput

subroutine test_unit()
integer, parameter :: mxUnits = 89
integer :: unit(mxUnits)
integer, parameter :: mxRandom = 5
integer, parameter :: Random(mxRandom) = (/ 4, 36, 91, 92, 95 /)
integer, parameter :: mxTaken = 30
integer, parameter :: taken(mxTaken) = (/  3,  9, 11, 21, 23, 25, 28, 30, 33, 35, &
                                          37, 39, 40, 42, 43, 45, 49, 52, 53, 55, &
                                          60, 61, 63, 64, 65, 66, 67, 69, 80, 82  /)
integer :: i, j
logical :: opened

! Test the get unit number routine
do k = 1, 2   ! Loop through this series twice to make sure things ok
   ! Open some random unit numbers
   do i = 1, mxRandom
      call open_file(random(i))
   end do
   ! First take a bunch of units with explicit unit numbers
   do i = 1, mxTaken
      j = shr_file_getUnit( taken(i) )
      call open_file(taken(i))
      if ( j /= taken(i) )then
         call shr_sys_abort( "error: get unit did NOT grab the correct unit" )
      end if
   end do
   ! Now loop through and take all other unit numbers
   do i = 1, mxUnits-mxTaken-mxRandom
      unit(i) = shr_file_getUnit()
      inquire(unit(i), opened=opened )
      if ( opened )then
         call shr_sys_abort( "error: get unit got a unit already opened" )
      end if
      call open_file(unit(i))
      do j = 1, mxTaken
         if ( unit(i) == taken(j) )then
            call shr_sys_abort( "error: get unit got a unit already taken" )
         end if
      end do
      do j = 1, i-1
         if ( unit(i) == unit(j) )then
            call shr_sys_abort( "error: get unit got a unit already taken" )
         end if
      end do
   end do
   ! Free units taken
   do i = 1, mxUnits-mxTaken-mxRandom
      call close_file(unit(i) )
      call shr_file_freeUnit( unit(i) )
   end do
   do i = 1, mxTaken
      call close_file(taken(i) )
      call shr_file_freeUnit( taken(i) )
   end do
   do i = 1, mxRandom
      call close_file(random(i))
   end do
end do
end subroutine test_unit

subroutine open_file(unit)
integer :: unit
character(len=256) :: tmp

write(6,*) "take unit", unit
write(tmp,"('tmp',i3.3,'.dat')") unit
open(unit, file=tmp, status="new")
end subroutine open_file

subroutine close_file(unit)
integer :: unit
close(unit,status="delete")
write(6,*) "free unit", unit
end subroutine close_file


end program test_shr_file
