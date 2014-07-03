module FracWetMod
!-----------------------------------------------------------------------
! !DESCRIPTION:
! Determine fraction of vegetated surfaces which are wet and
! fraction of elai which is dry.
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: FracWet
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
  subroutine FracWet(numf, filter)
!
! !DESCRIPTION:
! Determine fraction of vegetated surfaces which are wet and
! fraction of elai which is dry. The variable ``fwet'' is the
! fraction of all vegetation surfaces which are wet including
! stem area which contribute to evaporation. The variable ``fdry''
! is the fraction of elai which is dry because only leaves
! can transpire.  Adjusted for stem area which does not transpire.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: numf                  ! number of filter non-lake points
    integer, intent(in) :: filter(numf)          ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine Hydrology1 in module Hydrology1Mod
!
! !REVISION HISTORY:
! Created by Keith Oleson and M. Vertenstein
! 03/08/29 Mariana Vertenstein : Migrated to vectorized code
!
! !LOCAL VARIABLES:
    integer  :: fp,p             ! indices
    real(r8) :: vegt             ! frac_veg_nosno*lsai
    real(r8) :: dewmxi           ! inverse of maximum allowed dew [1/mm]
!-----------------------------------------------------------------------


   associate(& 
   frac_veg_nosno     => pps%frac_veg_nosno , & ! Input:  [integer (:)]  fraction of veg not covered by snow (0/1 now) [-]
   dewmx              => pps%dewmx          , & ! Input:  [real(r8) (:)]  Maximum allowed dew [mm]                
   elai               => pps%elai           , & ! Input:  [real(r8) (:)]  one-sided leaf area index with burying by snow
   esai               => pps%esai           , & ! Input:  [real(r8) (:)]  one-sided stem area index with burying by snow
   h2ocan             => pws%h2ocan         , & ! Input:  [real(r8) (:)]  total canopy water (mm H2O)             
   fwet               => pps%fwet           , & ! Output: [real(r8) (:)]  fraction of canopy that is wet (0 to 1) 
   fdry               => pps%fdry             & ! Output: [real(r8) (:)]  fraction of foliage that is green and dry [-] (new)
   )

    ! Compute fraction of canopy that is wet and dry

    do fp = 1,numf
       p = filter(fp)
       if (frac_veg_nosno(p) == 1) then
          if (h2ocan(p) > 0._r8) then
             vegt    = frac_veg_nosno(p)*(elai(p) + esai(p))
             dewmxi  = 1.0_r8/dewmx(p)
             fwet(p) = ((dewmxi/vegt)*h2ocan(p))**0.666666666666_r8
             fwet(p) = min (fwet(p),1.0_r8)   ! Check for maximum limit of fwet
          else
             fwet(p) = 0._r8
          end if
          fdry(p) = (1._r8-fwet(p))*elai(p)/(elai(p)+esai(p))
       else
          fwet(p) = 0._r8
          fdry(p) = 0._r8
       end if
    end do

    end associate 
   end subroutine FracWet

end module FracWetMod
