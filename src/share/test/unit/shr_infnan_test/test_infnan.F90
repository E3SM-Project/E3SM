program test_infnan

!
! This is a test for the shr_infnan_mod module. It was created using the
! pre-CTest system, with minimal changes to keep it working. So it may not
! be a great example of a CTest test now.
!

use shr_kind_mod, only: r8 => shr_kind_r8
use shr_kind_mod, only: r4 => shr_kind_r4
use shr_kind_mod, only: i8 => shr_kind_i8
use shr_kind_mod, only: i4 => shr_kind_i4
use shr_infnan_mod

implicit none

real(r8) :: x, zero
real(r4) :: y
real(r8) :: r8array(100), r82Darray(10,10), r83Darray(4,4,4)
real(r8) :: r84Darray(3,3,3,3), r85Darray(2,2,2,2,2)
real(r8) :: inf
real(r8) :: nan
real(r8) :: nans
real(r4) :: spnan
real(r4) :: spnans
integer(i8), parameter :: dpinfpat = int(O'0777600000000000000000',i8)
integer(i8), parameter :: dpnanpat = int(O'0777700000000000000000',i8)
integer(i8), parameter :: dpnanspat = int(O'0777610000000000000000',i8)
integer(i4), parameter :: spnanpat =  int(Z'7FC00000',i4)
integer(i4), parameter :: spnanspat = int(Z'7FC10000',i4)
intrinsic :: count

inf = transfer(dpinfpat,inf)
nan = transfer(dpnanpat,nan)
nans = transfer(dpnanspat,nans)
spnan = transfer( spnanpat,spnan)
spnans = transfer( spnanspat,spnans)

x    = 0.0
zero = 0.0

call assert(       shr_infnan_isnan( nan ),    "Test that value set to nan is nan" )
call assert(       shr_infnan_isnan( nans ),   "Test that value set to nans is nan" )
call assert(       shr_infnan_isnan( spnan ),  "Test that value set to sp nan is nan" )
call assert(       shr_infnan_isnan( spnans ), "Test that value set to sp nans is nan" )
call assert( .not. shr_infnan_isnan( 1.0_r8 ), "Test that value set to one is NOT nan" )
call assert( .not. shr_infnan_isnan( 1.0_r4 ), "Test that value set to SP one is NOT nan" )
call assert( .not. shr_infnan_isnan( huge(x) ), "Test that value set to huge is NOT nan" )
x    = 1.0/zero
call assert( .not. shr_infnan_isnan( x ),     "Test that 1/0 is NOT nan" )
x    = -1.0/zero
call assert( .not. shr_infnan_isnan( x ),    "Test that -1/0 is NOT nan" )

r8array(:)      = 1.0d00
r8array(10)     = nan
r8array(15)     = nan
r82Darray(:,:)  = 1.0d00
r82Darray(5,5)  = nan
r82Darray(10,7) = nan
r82Darray(7,9)  = nan
r83Darray(:,:,:) = 1.0d00
r83Darray(4,2,2) = nan
r83Darray(3,1,2) = nan
r83Darray(1,1,1) = nan
r83Darray(1,1,4) = nan
r84Darray(:,:,:,:) = 1.0d00
r84Darray(3,2,2,1) = nan
r84Darray(3,1,2,1) = nan
r84Darray(1,1,1,1) = nan
r84Darray(1,1,3,1) = nan
r84Darray(1,2,3,1) = nan
r85Darray(:,:,:,:,:) = 1.0d00
r85Darray(1,2,2,1,1) = nan
r85Darray(1,1,2,1,2) = nan
r85Darray(1,1,1,2,1) = nan
r85Darray(1,2,2,2,1) = nan
r85Darray(1,2,1,1,2) = nan
r85Darray(1,1,1,1,1) = nan
call assert(       any(shr_infnan_isnan( r8array )),      "Test that array with 2 nans is nan" )
call assert(       count(shr_infnan_isnan( r8array )) == 2,    "Test that there are 2 nans in that array" )
call assert(       any(shr_infnan_isnan( r82Darray )),    "Test that 2D array with 3 nans is nan" )
call assert(       count(shr_infnan_isnan( r82Darray )) == 3,  "Test that there are 3 nans in that array" )
call assert(       any(shr_infnan_isnan( r83Darray )),    "Test that 3D array with 4 nans is nan" )
call assert(       count(shr_infnan_isnan( r83Darray )) == 4,  "Test that there are 4 nans in that array" )
call assert(       any(shr_infnan_isnan( r84Darray )),    "Test that 4D array with 5 nans is nan" )
call assert(       count(shr_infnan_isnan( r84Darray )) == 5,  "Test that there are 5 nans in that array" )
call assert(       any(shr_infnan_isnan( r85Darray )),    "Test that 5D array with 6 nans is nan" )
call assert(       count(shr_infnan_isnan( r85Darray )) == 6,  "Test that there are 6 nans in that array" )
call assert(       shr_infnan_isposinf( inf ),   "Test that value set to inf is inf" )
call assert( .not. shr_infnan_isposinf( 1.0_r8 ), "Test that value set to one is NOT inf" )
call assert( .not. shr_infnan_isposinf( 1.0_r4 ), "Test that value set to SP one is NOT inf" )
call assert(       shr_infnan_isneginf( -inf ),   "Test that value set to -inf is -inf" )
call assert( .not. shr_infnan_isneginf( 1.0_r8 ), "Test that value set to one is NOT -inf" )
call assert( .not. shr_infnan_isneginf( 1.0_r4 ), "Test that value set to SP one is NOT -inf" )
x = 1.0/zero
call assert(       shr_infnan_isposinf( x ),      "Test that 1/0 is inf" )
x = -1.0/zero
call assert(       shr_infnan_isneginf( x ),     "Test that -1/0 is -inf" )

x = -1.0
call assert(       shr_infnan_isnan( sqrt(x) ), "Test that sqrt-1 is nan" )
call assert(       shr_infnan_isnan( log(x) ),  "Test that log-1 is nan" )

x = shr_infnan_nan
call assert(       shr_infnan_isnan( x ),  "Test that shr_infnan_nan sets r8 to nan" )
y = shr_infnan_nan
call assert(       shr_infnan_isnan( y ),  "Test that shr_infnan_nan sets r4 to nan" )

x = shr_infnan_inf
call assert(       shr_infnan_isinf( x ),  "Test that shr_infnan_inf sets r8 to inf" )
y = shr_infnan_inf
call assert(       shr_infnan_isinf( y ),  "Test that shr_infnan_inf sets r4 to inf" )

x = shr_infnan_posinf
call assert(       shr_infnan_isposinf( x ),  "Test that shr_infnan_posinf sets r8 to +inf" )
y = shr_infnan_posinf
call assert(       shr_infnan_isposinf( y ),  "Test that shr_infnan_posinf sets r4 to +inf" )

x = shr_infnan_neginf
call assert(       shr_infnan_isneginf( x ),  "Test that shr_infnan_neginf sets r8 to -inf" )
y = shr_infnan_neginf
call assert(       shr_infnan_isneginf( y ),  "Test that shr_infnan_neginf sets r4 to -inf" )

x = shr_infnan_to_r8(shr_infnan_qnan)
call assert(       shr_infnan_isnan( x ),  "Test that shr_infnan_to_r8(shr_infnan_qnan) sets r8 to nan" )
y = shr_infnan_to_r4(shr_infnan_qnan)
call assert(       shr_infnan_isnan( y ),  "Test that shr_infnan_to_r4(shr_infnan_qnan) sets r4 to nan" )

x = shr_infnan_to_r8(shr_infnan_snan)
call assert(       shr_infnan_isnan( x ),  "Test that shr_infnan_to_r8(shr_infnan_snan) sets r8 to nan" )
y = shr_infnan_to_r4(shr_infnan_snan)
call assert(       shr_infnan_isnan( y ),  "Test that shr_infnan_to_r4(shr_infnan_snan) sets r4 to nan" )

x = shr_infnan_to_r8(shr_infnan_posinf)
call assert(       shr_infnan_isposinf( x ),  "Test that shr_infnan_to_r8(shr_infnan_posinf) sets r8 to +inf" )
y = shr_infnan_to_r4(shr_infnan_posinf)
call assert(       shr_infnan_isposinf( y ),  "Test that shr_infnan_to_r4(shr_infnan_posinf) sets r4 to +inf" )

x = shr_infnan_to_r8(shr_infnan_neginf)
call assert(       shr_infnan_isneginf( x ),  "Test that shr_infnan_to_r8(shr_infnan_neginf) sets r8 to -inf" )
y = shr_infnan_to_r4(shr_infnan_neginf)
call assert(       shr_infnan_isneginf( y ),  "Test that shr_infnan_to_r4(shr_infnan_neginf) sets r4 to -inf" )

contains

  subroutine assert(val, msg)
    logical, intent(in) :: val
    character(len=*), intent(in) :: msg

    if (.not. val) then
       print *, msg
       stop 1
    end if

  end subroutine assert

end program test_infnan
