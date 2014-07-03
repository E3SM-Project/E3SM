module CNDVInitMod

  !-----------------------------------------------------------------------
  ! !MODULE: initCNDVMod
  !
  ! !DESCRIPTION:
  ! Contains cold start initial values and time constant (and flux / diagnostic vars) 
  ! for CNDV scheme.
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: initColdCNDV
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initColdCNDV(bounds)
    !
    ! !DESCRIPTION:
    ! CNDV related initializations
    !
    ! !USES:
    use clmtype       , only : pdgvs, pcs
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use shr_const_mod , only : SHR_CONST_TKFRZ
    use decompMod     , only : bounds_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: p           ! pft index
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       pdgvs%present(p)   = .false.
       pdgvs%crownarea(p) = 0._r8
       pdgvs%nind(p)      = 0._r8
       pdgvs%t_mo_min(p)  = 1.0e+36_r8
       pdgvs%agdd20(p)    = 0._r8
       pdgvs%tmomin20(p)  = SHR_CONST_TKFRZ - 5._r8 !initialize this way for Phenology code
       pcs%leafcmax(p)    = 0._r8
    end do

  end subroutine initColdCNDV

end module CNDVInitMod
