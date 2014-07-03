module H2OSfcMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate surface water hydrology
  !
  ! !USES:
  use shr_log_mod   , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: FracH2oSfc     ! Calculate fraction of land surface that is wet
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine FracH2oSfc(bounds, num_h2osfc, filter_h2osfc, frac_h2osfc, no_update)
    !
    ! !DESCRIPTION:
    ! Determine fraction of land surfaces which are submerged  
    ! based on surface microtopography and surface water storage.
    !
    ! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use shr_const_mod, only : shr_const_pi
    use shr_spfn_mod , only : erf => shr_spfn_erf
    use clm_varcon   , only : istsoil, istcrop
    use decompMod    , only : bounds_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in)  :: bounds                      ! bounds
    integer , intent(in)           :: num_h2osfc                  ! number of column points in column filter
    integer , intent(in)           :: filter_h2osfc(:)            ! column filter 
    real(r8), intent(inout)        :: frac_h2osfc( bounds%begc: ) ! fractional surface water (mm) [col]
    integer , intent(in), optional :: no_update                   ! flag to make calculation w/o updating variables
    !
    ! !LOCAL VARIABLES:
    integer :: c,f,l          ! indices
    real(r8):: d,fd,dfdd      ! temporary variable for frac_h2oscs iteration
    real(r8):: sigma          ! microtopography pdf sigma in mm
    real(r8):: min_h2osfc
    !-----------------------------------------------------------------------

   SHR_ASSERT_ALL((ubound(frac_h2osfc) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))

   associate(& 
   h2osfc        =>  cws%h2osfc        , & ! Input:  [real(r8) (:)]  surface water (mm)                                
   h2osno        =>  cws%h2osno        , & ! Input:  [real(r8) (:)]  snow water (mm H2O)                               
   micro_sigma   =>  cps%micro_sigma   , & ! Input:  [real(r8) (:)]  microtopography pdf sigma (m)                     
   h2osoi_liq    =>  cws%h2osoi_liq    , & ! Output: [real(r8) (:,:)]  liquid water (col,lyr) [kg/m2]                  
   frac_sno      =>  cps%frac_sno      , & ! Output: [real(r8) (:)]  fraction of ground covered by snow (0 to 1)       
   frac_sno_eff  =>  cps%frac_sno_eff    & ! Output: [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)  
   )
 
    ! arbitrary lower limit on h2osfc for safer numerics...
    min_h2osfc=1.e-8_r8

    do f = 1, num_h2osfc
       c = filter_h2osfc(f)
       l = col%landunit(c)
       ! h2osfc only calculated for soil vegetated land units
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          !  Use newton-raphson method to iteratively determine frac_h20sfc
          !  based on amount of surface water storage (h2osfc) and 
          !  microtopography variability (micro_sigma)
          if (h2osfc(c) > min_h2osfc) then
             ! a cutoff is needed for numerical reasons...(nonconvergence after 5 iterations)
             d=0.0

             sigma=1.0e3*micro_sigma(c) ! convert to mm
             do l=1,10
                fd = 0.5*d*(1.0_r8+erf(d/(sigma*sqrt(2.0)))) &
                        +sigma/sqrt(2.0*shr_const_pi)*exp(-d**2/(2.0*sigma**2)) &
                        -h2osfc(c)
                dfdd = 0.5*(1.0_r8+erf(d/(sigma*sqrt(2.0))))
                
                d = d - fd/dfdd
             enddo
             !--  update the submerged areal fraction using the new d value
             frac_h2osfc(c) = 0.5*(1.0_r8+erf(d/(sigma*sqrt(2.0))))

          else
             frac_h2osfc(c) = 0._r8
             h2osoi_liq(c,1) = h2osoi_liq(c,1) + h2osfc(c)
             h2osfc(c)=0._r8
          endif

          if (.not. present(no_update)) then

             ! adjust fh2o, fsno when sum is greater than zero
             if (frac_sno(c) > (1._r8 - frac_h2osfc(c)) .and. h2osno(c) > 0) then

                if (frac_h2osfc(c) > 0.01_r8) then             
                   frac_h2osfc(c) = max(1.0_r8 - frac_sno(c),0.01_r8)
                   frac_sno(c) = 1.0_r8 - frac_h2osfc(c)
                else
                   frac_sno(c) = 1.0_r8 - frac_h2osfc(c)
                endif
                frac_sno_eff(c)=frac_sno(c)
             
             endif

          endif ! end of no_update construct

       else !if landunit not istsoil/istcrop, set frac_h2osfc to zero
          frac_h2osfc(c) = 0._r8
       endif

    end do
       
    end associate 
  end subroutine FracH2oSfc

end module H2OSfcMod
