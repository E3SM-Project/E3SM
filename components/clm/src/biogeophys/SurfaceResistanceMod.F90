module SurfaceResistanceMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines for calculation of surface resistances of the different tracers
  ! transported with BeTR. The surface here refers to water and soil, not including canopy
  !
  ! !USES:
  use shr_kind_mod  , only: r8 => shr_kind_r8
  use shr_const_mod , only: SHR_CONST_TKFRZ
  use clm_varctl    , only: iulog
  use SoilStateType , only: soilstate_type
  use WaterStateType, only: waterstate_type 
  
  implicit none
  save
  private
  integer :: soil_stress_method   !choose the method for soil resistance calculation
  
  integer, parameter :: leepielke_1992 = 0 !
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: calc_soilevap_stress
  public :: do_soilevap_beta
  public :: init_soil_stress
  !
  ! !REVISION HISTORY:
  ! 6/25/2013 Created by Jinyun Tang
  !-----------------------------------------------------------------------
  
contains

  !-----------------------------------------------------------------------
  subroutine init_soil_stress()
   !
   !DESCRIPTIONS
   ! initialize method for soil stress calculation
   implicit none
   
   soil_stress_method = leepielke_1992
   
   end subroutine init_soil_stress
   
   !------------------------------------------------------------------------------   
   subroutine calc_soilevap_stress(bounds, num_nolakec, filter_nolakec, &
        soilstate_vars, waterstate_vars)
     !
     ! DESCRIPTIONS
     ! compute the stress factor for soil evaporation calculation
     !
     use shr_kind_mod  , only : r8 => shr_kind_r8     
     use shr_const_mod , only : SHR_CONST_PI  
     use decompMod     , only : bounds_type
     use ColumnType    , only : col
     use LandunitType  , only : lun
     use abortutils    , only : endrun      
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type)     , intent(in)    :: bounds    ! bounds   
     integer               , intent(in)    :: num_nolakec
     integer               , intent(in)    :: filter_nolakec(:)
     type(soilstate_type)  , intent(inout) :: soilstate_vars
     type(waterstate_type) , intent(in)    :: waterstate_vars

     character(len=32) :: subname = 'calc_soilevap_stress'  ! subroutine name
     associate(                &
          soilbeta =>  soilstate_vars%soilbeta_col  & ! Output: [real(r8) (:)] factor that reduces ground evaporation
          )
   
       !select the right method and do the calculation
       select case (soil_stress_method)

       case (leepielke_1992)
          call calc_beta_leepielke1992(bounds, num_nolakec, filter_nolakec, &
               soilstate_vars, waterstate_vars, soilbeta(bounds%begc:bounds%endc))

       case default
          call endrun(subname // ':: a soilevap stress function must be specified!')     
       end select

     end associate

   end subroutine calc_soilevap_stress
   
   !------------------------------------------------------------------------------   
   subroutine calc_beta_leepielke1992(bounds, num_nolakec, filter_nolakec, &
        soilstate_vars, waterstate_vars, soilbeta)
     !
     ! DESCRIPTION
     ! compute the lee-pielke beta factor to scal actual soil evaporation from potential evaporation
     !
     ! USES
     use shr_kind_mod    , only : r8 => shr_kind_r8     
     use shr_const_mod   , only : SHR_CONST_PI
     use shr_log_mod     , only : errMsg => shr_log_errMsg   
     use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
     use decompMod       , only : bounds_type
     use clm_varcon      , only : denh2o, denice
     use landunit_varcon , only : istice, istice_mec, istwet, istsoil, istcrop
     use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
     use column_varcon   , only : icol_road_imperv, icol_road_perv
     use ColumnType      , only : col
     use LandunitType    , only : lun
     !
     implicit none
     type(bounds_type)     , intent(in)    :: bounds    ! bounds   
     integer               , intent(in)    :: num_nolakec
     integer               , intent(in)    :: filter_nolakec(:)
     type(soilstate_type)  , intent(in)    :: soilstate_vars
     type(waterstate_type) , intent(in)    :: waterstate_vars
     real(r8)              , intent(inout) :: soilbeta(bounds%begc:bounds%endc)

     !local variables
     real(r8) :: fac, fac_fc, wx      !temporary variables
     integer  :: c, l, fc     !indices

     SHR_ASSERT_ALL((ubound(soilbeta)    == (/bounds%endc/)), errMsg(__FILE__, __LINE__))

     associate(                                              &
          watsat      =>    soilstate_vars%watsat_col      , & ! Input:  [real(r8) (:,:)] volumetric soil water at saturation (porosity)
          watfc       =>    soilstate_vars%watfc_col       , & ! Input:  [real(r8) (:,:)] volumetric soil water at field capacity
          
          h2osoi_ice  =>    waterstate_vars%h2osoi_ice_col , & ! Input:  [real(r8) (:,:)] ice lens (kg/m2)                       
          h2osoi_liq  =>    waterstate_vars%h2osoi_liq_col , & ! Input:  [real(r8) (:,:)] liquid water (kg/m2)                   
          frac_sno    =>    waterstate_vars%frac_sno_col   , & ! Input:  [real(r8) (:)] fraction of ground covered by snow (0 to 1)
          frac_h2osfc =>    waterstate_vars%frac_h2osfc_col  & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
          )

       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = col%landunit(c)   
          if (lun%itype(l)/=istwet .AND. lun%itype(l)/=istice  &
               .AND. lun%itype(l)/=istice_mec) then
             if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
                wx   = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/col%dz(c,1)
                fac  = min(1._r8, wx/watsat(c,1))
                fac  = max( fac, 0.01_r8 )
                !! Lee and Pielke 1992 beta, added by K.Sakaguchi
                if (wx < watfc(c,1) ) then  !when water content of ths top layer is less than that at F.C.
                   fac_fc  = min(1._r8, wx/watfc(c,1))  !eqn5.66 but divided by theta at field capacity
                   fac_fc  = max( fac_fc, 0.01_r8 )
                   ! modify soil beta by snow cover. soilbeta for snow surface is one
                   soilbeta(c) = (1._r8-frac_sno(c)-frac_h2osfc(c)) &
                        *0.25_r8*(1._r8 - cos(SHR_CONST_PI*fac_fc))**2._r8 &
                        + frac_sno(c)+ frac_h2osfc(c)
                else   !when water content of ths top layer is more than that at F.C.
                   soilbeta(c) = 1._r8
                end if
             else if (col%itype(c) == icol_road_perv) then
                soilbeta(c) = 0._r8
             else if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall) then
                soilbeta(c) = 0._r8          
             else if (col%itype(c) == icol_roof .or. col%itype(c) == icol_road_imperv) then
                soilbeta(c) = 0._r8
             endif
          else
             soilbeta(c) =   1._r8
          endif
       enddo

     end associate

   end subroutine calc_beta_leepielke1992
   
   !------------------------------------------------------------------------------   
   function do_soilevap_beta()result(lres)
     !
     !DESCRIPTION
     ! return true if the moisture stress for soil evaporation is computed as beta factor
     ! otherwise false
     implicit none
     logical :: lres

     if(soil_stress_method==leepielke_1992)then
        lres=.true.
     else
        lres=.false.
     endif
     return

   end function do_soilevap_beta

end module SurfaceResistanceMod
