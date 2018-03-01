module infnan

  ! Relabel shr_infnan public members, just for the convenience
  ! of having shorter names.

use shr_infnan_mod, only: &
     isnan => shr_infnan_isnan, &
     isinf => shr_infnan_isinf, &
     isposinf => shr_infnan_isposinf, &
     isneginf => shr_infnan_isneginf, &
     nan => shr_infnan_nan, &
     inf => shr_infnan_inf, &
     posinf => shr_infnan_posinf, &
     neginf => shr_infnan_neginf, &
     assignment(=) 

implicit none
private
save

! Weird thing here: if you make the module public,
! ifort 12 has an ICE. But ifort 13 has been fixed.
public :: isnan
public :: isinf
public :: isposinf
public :: isneginf
public :: nan
public :: inf
public :: posinf
public :: neginf
public :: assignment(=)

end module infnan
