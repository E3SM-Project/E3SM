module mkvarpar

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varpar
!
! !DESCRIPTION:
! Module containing CLM parameters
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_REARTH
!
! !PUBLIC TYPES:
  implicit none
  save
!
  integer, parameter :: nlevsoi = 10  ! number of soil layers
  integer, parameter :: numstdpft = 16! number of standard PFT types
  integer, parameter :: noveg   = 0   ! value for non-vegetated pft
  integer, parameter :: nlevurb = 5   ! number of urban layers
  integer, parameter :: numsolar = 2  ! number of solar types (Direct,Diffuse)
  integer, parameter :: numrad = 2    ! number of solar bands (VIS,NIR)
  real(r8),parameter :: elev_thresh  = 2600._r8   ! elevation threshold for screening urban areas
  real(r8),parameter :: re = SHR_CONST_REARTH*0.001

!
!EOP
!-----------------------------------------------------------------------

end module mkvarpar
