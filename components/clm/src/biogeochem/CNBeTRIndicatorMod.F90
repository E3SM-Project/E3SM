module CNBeTRIndicatorMod
  !
  ! DESCRIPTION
  ! code to switch on/off of gap mortality and phenology relevant processes.

  use shr_kind_mod        , only : r8 => shr_kind_r8
implicit none
public
integer, parameter :: pid_leafn_to_litter = 1
integer, parameter :: pid_frootn_to_litter = 2
integer, parameter :: pid_livestemn_to_litter =3

integer, parameter :: gid_m_leafn_to_litter  = 1
integer, parameter :: gid_m_frootn_to_litter = 2
integer, parameter :: gid_m_livestemn_to_litter = 3
integer, parameter :: gid_m_deadstemn_to_litter = 4
integer, parameter :: gid_m_livecrootn_to_litter = 5
integer, parameter :: gid_m_deadcrootn_to_litter = 6
integer, parameter :: gid_m_retransn_to_litter = 7
integer, parameter :: gid_m_leafn_storage_to_litter = 8
integer, parameter :: gid_m_frootn_storage_to_litter = 9
integer, parameter :: gid_m_livestemn_storage_to_litter = 10
integer, parameter :: gid_m_deadstemn_storage_to_litter = 11
integer, parameter :: gid_m_livecrootn_storage_to_litter = 12
integer, parameter :: gid_m_deadcrootn_storage_to_litter = 13
integer, parameter :: gid_m_leafn_xfer_to_litter = 14
integer, parameter :: gid_m_frootn_xfer_to_litter = 15
integer, parameter :: gid_m_livestemn_xfer_to_litter = 16
integer, parameter :: gid_m_deadstemn_xfer_to_litter = 17
integer, parameter :: gid_m_livecrootn_xfer_to_litter = 18
integer, parameter :: gid_m_deadcrootn_xfer_to_litter = 19

real(r8) :: pheno_indicator(3)
real(r8) :: gap_indicator(19)

contains
  subroutine set_pheno_indicators
  implicit none

  pheno_indicator(:) = 1._r8
  return
  pheno_indicator(pid_leafn_to_litter) = 0._r8
  pheno_indicator(pid_frootn_to_litter) = 0._r8
  pheno_indicator(pid_livestemn_to_litter) = 0._r8

  end subroutine set_pheno_indicators

!----------------------------------------------------------------------
  subroutine set_gap_indicators
  implicit none

  gap_indicator(:) = 1._r8
  return
  gap_indicator(gid_m_leafn_to_litter)  = 0._r8
  gap_indicator(gid_m_frootn_to_litter) = 0._r8
  gap_indicator(gid_m_livestemn_to_litter) = 0._r8
  gap_indicator(gid_m_deadstemn_to_litter) = 0._r8
  gap_indicator(gid_m_livecrootn_to_litter) = 0._r8
  gap_indicator(gid_m_deadcrootn_to_litter) = 0._r8
  gap_indicator(gid_m_retransn_to_litter) = 0._r8
  gap_indicator(gid_m_leafn_storage_to_litter) = 0._r8
  gap_indicator(gid_m_frootn_storage_to_litter) = 0._r8
  gap_indicator(gid_m_livestemn_storage_to_litter) = 0._r8
  gap_indicator(gid_m_deadstemn_storage_to_litter) = 0._r8
  gap_indicator(gid_m_livecrootn_storage_to_litter) = 0._r8
  gap_indicator(gid_m_deadcrootn_storage_to_litter) = 0._r8
  gap_indicator(gid_m_leafn_xfer_to_litter) = 0._r8
  gap_indicator(gid_m_frootn_xfer_to_litter) = 0._r8
  gap_indicator(gid_m_livestemn_xfer_to_litter) = 0._r8
  gap_indicator(gid_m_deadstemn_xfer_to_litter) = 0._r8
  gap_indicator(gid_m_livecrootn_xfer_to_litter) = 0._r8
  gap_indicator(gid_m_deadcrootn_xfer_to_litter) = 0._r8

  end subroutine set_gap_indicators


end module CNBeTRIndicatorMod
