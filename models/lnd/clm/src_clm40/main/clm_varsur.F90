
module clm_varsur

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varsur
!
! !DESCRIPTION:
! Module containing 2-d surface boundary data information
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! land model grid - moved to domainMod
!
! surface boundary data, these are all "gdc" local 
!
  integer , allocatable :: vegxy(:,:) ! vegetation type
  real(r8), allocatable,target :: wtxy(:,:)  ! subgrid weights

  real(r8),allocatable :: pctspec(:)         ! percent of spec lunits wrt gcell

  real(r8), allocatable,target :: topoxy(:,:)  ! subgrid glacier_mec sfc elevation
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2005-11-01 Moved grid to domainMod, T Craig
!
!EOP
!-----------------------------------------------------------------------

end module clm_varsur
