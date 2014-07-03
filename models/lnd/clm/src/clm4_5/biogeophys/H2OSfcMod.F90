module H2OSfcMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: H2OSfcMod
!
! !DESCRIPTION:
! Calculate surface water hydrology
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: FracH2oSfc     ! Calculate fraction of land surface that is wet
!
! !REVISION HISTORY:
! Created by 09/15/07 Sean Swenson
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: FracH2oSfc
!
! !INTERFACE:
  subroutine FracH2oSfc(lbc, ubc, num_h2osfc, filter_h2osfc,frac_h2osfc,no_update)
!
! !DESCRIPTION:
! Determine fraction of land surfaces which are submerged  
! based on surface microtopography and surface water storage.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use shr_const_mod       , only : shr_const_pi
    use clm_varcon          , only : istsoil, istcrop
    use clm_varctl,   only: iulog
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbc, ubc                  ! column bounds
    integer , intent(in) :: num_h2osfc                ! number of column points in column filter
    integer , intent(in) :: filter_h2osfc(ubc-lbc+1)  ! column filter 
    real(r8), intent(inout) :: frac_h2osfc(lbc:ubc)   ! fractional surface water (mm)
    integer , intent(in), optional :: no_update       ! flag to make calculation w/o updating variables
!
! !CALLED FROM:
! subroutine Hydrology1 in module Hydrology1Mod
!
! !REVISION HISTORY:
! 09/24/07 Created by S. Swenson 
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    real(r8), pointer :: h2osfc(:)         ! surface water (mm)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: micro_sigma(:)    ! microtopography pdf sigma (m)
    real(r8), pointer :: frac_sno(:)       ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: frac_sno_eff(:)       ! eff. fraction of ground covered by snow (0 to 1)
    integer , pointer :: snl(:)            ! minus number of snow layers
    real(r8), pointer :: h2osno(:)         ! snow water (mm H2O)
    real(r8), pointer :: h2osoi_liq(:,:)   ! liquid water (col,lyr) [kg/m2]
    real(r8), pointer :: topo_slope(:)     ! topographic slope
    real(r8), pointer :: topo_ndx(:)       ! topographic slope
    integer , pointer :: ltype(:)          ! landunit type
    integer , pointer :: clandunit(:)      ! columns's landunit

    !intrinsic :: derf
    real(r8)  :: derf
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: c,f,l         ! indices
    real(r8):: d,fd,dfdd      ! temporary variable for frac_h2oscs iteration
    real(r8):: sigma          ! microtopography pdf sigma in mm
    real(r8):: min_h2osfc,minslope,maxslope,temp_norm,slopemax

!-----------------------------------------------------------------------

! Assign local pointers to derived subtypes components (column-level)

    h2osoi_liq          => clm3%g%l%c%cws%h2osoi_liq
    h2osfc              => clm3%g%l%c%cws%h2osfc
    micro_sigma         => clm3%g%l%c%cps%micro_sigma
    topo_slope          => clm3%g%l%c%cps%topo_slope
    topo_ndx            => clm3%g%l%c%cps%topo_ndx
    ltype               => clm3%g%l%itype
    clandunit           => clm3%g%l%c%landunit

    frac_sno            => clm3%g%l%c%cps%frac_sno 
    frac_sno_eff        => clm3%g%l%c%cps%frac_sno_eff
    snl                 => clm3%g%l%c%cps%snl
    h2osno              => clm3%g%l%c%cws%h2osno
 
    ! arbitrary lower limit on h2osfc for safer numerics...
    min_h2osfc=1.e-8_r8

    do f = 1, num_h2osfc
       c = filter_h2osfc(f)
       l = clandunit(c)
       ! h2osfc only calculated for soil vegetated land units
       if (ltype(l) == istsoil .or. ltype(l) == istcrop) then

          !  Use newton-raphson method to iteratively determine frac_h20sfc
          !  based on amount of surface water storage (h2osfc) and 
          !  microtopography variability (micro_sigma)
          if (h2osfc(c) > min_h2osfc) then
             ! a cutoff is needed for numerical reasons...(nonconvergence after 5 iterations)
             d=0.0

             sigma=1.0e3*micro_sigma(c) ! convert to mm
             do l=1,10
                fd = 0.5*d*(1.0_r8+derf(d/(sigma*sqrt(2.0)))) &
                        +sigma/sqrt(2.0*shr_const_pi)*exp(-d**2/(2.0*sigma**2)) &
                        -h2osfc(c)
                dfdd = 0.5*(1.0_r8+derf(d/(sigma*sqrt(2.0))))
                
                d = d - fd/dfdd
             enddo
             !--  update the submerged areal fraction using the new d value
             frac_h2osfc(c) = 0.5*(1.0_r8+derf(d/(sigma*sqrt(2.0))))

          else
             frac_h2osfc(c) = 0._r8
             h2osoi_liq(c,1) = h2osoi_liq(c,1) + h2osfc(c)
             h2osfc(c)=0._r8
          endif

          if (.not. present(no_update)) then

             ! adjust fh2o, fsno when sum is greater than zero
             ! energy balance error when h2osno > 0 and snl = 0
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
       
  end subroutine FracH2oSfc

end module H2OSfcMod
