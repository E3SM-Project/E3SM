
module CNDVEcosystemDyniniMod

#if (defined CNDV)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNDVEcosystemDyniniMod
!
! !DESCRIPTION:
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public  :: CNDVEcosystemDynini ! CNDV related initializations
!
! !REVISION HISTORY:
! Created by Sam Levis following DGVMEcosystemDynMod by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNDVEcosystemDynini
!
! !INTERFACE:
  subroutine CNDVEcosystemDynini()
!
! !DESCRIPTION:
! CNDV related initializations
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use decompMod    , only : get_proc_bounds, get_proc_global
    use shr_const_mod, only : SHR_CONST_PI, SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from LPJ initialization subroutines)
!         Sam Levis (adapted for CNDV coupling; eliminated redunant parameters)
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: g,p,n           ! indices
    integer  :: begp, endp      ! per-proc beginning and ending pft indices
    integer  :: begc, endc      !              "                column indices
    integer  :: begl, endl      !              "                landunit indices
    integer  :: begg, endg      !              "                gridcell indices
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    pptr => clm3%g%l%c%p

    ! ---------------------------------------------------------------
    ! Some of the following came from LPJ subroutine initgrid
    ! ---------------------------------------------------------------

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    do p = begp,endp
       pptr%pdgvs%present(p)   = .false.
       pptr%pdgvs%crownarea(p) = 0._r8
       pptr%pdgvs%nind(p)      = 0._r8
       pptr%pcs%leafcmax(p)    = 0._r8
       pptr%pdgvs%t_mo_min(p)  = 1.0e+36_r8
    end do

    do g = begg,endg
       gptr%gdgvs%agdd20(g)   = 0._r8
       gptr%gdgvs%tmomin20(g) = SHR_CONST_TKFRZ - 5._r8 !initialize this way for Phenology code
    end do

  end subroutine CNDVEcosystemDynini

#endif

end module CNDVEcosystemDyniniMod
