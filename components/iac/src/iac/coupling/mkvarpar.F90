module mkvarpar

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkvarpar
!
! !DESCRIPTION:
! Module containing CLM parameters
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
  integer, parameter :: nlevsoi = 10  ! number of soil layers
  integer, parameter :: numpft  = 16  ! number of plant types
  integer, parameter :: noveg   = 0   ! value for non-vegetated pft
  integer, parameter :: nlevurb = 15  ! number of urban layers
  integer, parameter :: numsolar = 2  ! number of solar types (Direct,Diffuse)
  integer, parameter :: numrad = 2    ! number of solar bands (VIS,NIR)
  real(r8),parameter :: elev_thresh  = 2040._r8   ! elevation threshold for screening urban areas
  real(R8),parameter :: SHR_CONST_REARTH  = 6.37122e6_R8    ! radius of earth ~ m
  real(r8),parameter :: re = SHR_CONST_REARTH*0.001  ! radius of earth (km)
!
!EOP
!-----------------------------------------------------------------------

end module mkvarpar
