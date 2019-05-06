module radconstants

! provide stubs to allow building with no radiation scheme active

use shr_kind_mod,   only: r8 => shr_kind_r8
use cam_abortutils, only: endrun

implicit none
private
save

integer, parameter, public :: nswbands = 1
integer, parameter, public :: nlwbands = 1
integer, parameter, public :: idx_sw_diag = 1
integer, parameter, public :: idx_lw_diag = 1
integer, parameter, public :: idx_nir_diag = 1
integer, parameter, public :: idx_uv_diag = 1
integer, parameter, public :: nrh = 1
integer, parameter, public :: ot_length = 32

public :: rad_gas_index

integer, public, parameter :: gasnamelength = 1
integer, public, parameter :: nradgas = 1
character(len=gasnamelength), public, parameter :: gaslist(nradgas) &
   = (/' '/)

!### M. Branson additions so that radae.F90 will compile successfully
! Index of volc. cnt. abs. 0500--0650 cm-1
!integer, public, parameter :: idx_LW_0500_0650=3
integer, public, parameter :: idx_LW_0500_0650=1
! Index of volc. cnt. abs. 0650--0800 cm-1
!integer, public, parameter :: idx_LW_0650_0800=4
integer, public, parameter :: idx_LW_0650_0800=1
! Index of volc. cnt. abs. 0800--1000 cm-1
!integer, public, parameter :: idx_LW_0800_1000=5
integer, public, parameter :: idx_LW_0800_1000=1
! Index of volc. cnt. abs. 1000--1200 cm-1
!integer, public, parameter :: idx_LW_1000_1200=6
integer, public, parameter :: idx_LW_1000_1200=1
! Index of volc. cnt. abs. 1200--2000 cm-1
!integer, public, parameter :: idx_LW_1200_2000=7
integer, public, parameter :: idx_LW_1200_2000=1

!========================================================================================
contains
!========================================================================================

integer function rad_gas_index(gasname)

   character(len=*),intent(in) :: gasname

   call endrun('rad_gas_index: ERROR: this is a stub')

end function rad_gas_index

end module radconstants
