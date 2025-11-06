module SurfaceResistanceMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines for calculation of surface resistances of the different tracers
  ! transported with BeTR. The surface here refers to water and soil, not including canopy
  !
  ! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_const_mod , only : SHR_CONST_TKFRZ
  use elm_varctl    , only : iulog
  use SoilStateType , only : soilstate_type
  use WaterStateType, only : waterstate_type
  use ColumnDataType, only : col_ws
  use abortutils        , only : endrun


  implicit none
  save
  private
  integer :: soil_stress_method   !choose the method for soil resistance calculation
  integer, parameter :: leepielke_1992 = 0 !
  !$acc declare create(soil_stress_method)

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: calc_soilevap_stress
  public :: do_soilevap_beta
  public :: init_soil_stress
  public :: getlblcef
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
        soilstate_vars)
     !
     ! DESCRIPTIONS
     ! compute the stress factor for soil evaporation calculation
     !
      !$acc routine seq
     use shr_kind_mod  , only : r8 => shr_kind_r8
     use shr_const_mod , only : SHR_CONST_PI
     use decompMod     , only : bounds_type
     use ColumnType    , only : col_pp
     use LandunitType  , only : lun_pp
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type)     , intent(in)    :: bounds    ! bounds
     integer               , intent(in)    :: num_nolakec
     integer               , intent(in)    :: filter_nolakec(:)
     type(soilstate_type)  , intent(inout) :: soilstate_vars

     !character(len=32) :: subname = 'calc_soilevap_stress'  ! subroutine name
     associate(                &
          soilbeta =>  soilstate_vars%soilbeta_col,  & ! Output: [real(r8) (:)] factor that reduces ground evaporation
          dsl      =>  soilstate_vars%dsl_col     ,  & ! Output: [real(r8) (:)] soil dry surface layer thickness
          soilresis=>  soilstate_vars%soilresis_col  & ! Output: [real(r8) (:)] soil evaporation resistance
          )

       !select the right method and do the calculation
       select case (soil_stress_method)

       case (leepielke_1992)
          call calc_beta_leepielke1992(bounds, num_nolakec, filter_nolakec, &
               soilstate_vars, soilbeta(bounds%begc:bounds%endc))

       case (sl_14)
          call cal_soil_resistance_sl14(bounds, num_nolakec, filter_nolakec, &
               soilstate_vars, waterstate_vars, temperature_vars, &
               dsl(bounds%begc:bounds%endc), soilresis(bounds%begc,bounds%endc))

       case (sz_09)
          call calc_soil_resistance_sz09(bounds, num_nolakec, filter_nolakec, &
               soilstate_vars, waterstate_vars, &
               dsl(bounds%begc:bounds%endc), soilresis(bounds%begc:bounds%endc))

       case default
#ifndef _OPENACC
          call endrun('calc_soilevap_stress' //':: a soilevap stress function must be specified!')
#endif
       end select

     end associate

   end subroutine calc_soilevap_stress

   !------------------------------------------------------------------------------
   subroutine calc_beta_leepielke1992(bounds, num_nolakec, filter_nolakec, &
        soilstate_vars, soilbeta)
     !
     ! DESCRIPTION
     ! compute the lee-pielke beta factor to scal actual soil evaporation from potential evaporation
     !
     ! USES
      !$acc routine seq
     use shr_kind_mod    , only : r8 => shr_kind_r8
     use shr_const_mod   , only : SHR_CONST_PI
     use decompMod       , only : bounds_type
     use elm_varcon      , only : denh2o, denice
     use landunit_varcon , only : istice, istice_mec, istwet, istsoil, istcrop
     use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
     use column_varcon   , only : icol_road_imperv, icol_road_perv
     use ColumnType      , only : col_pp
     use LandunitType    , only : lun_pp
     use elm_varctl      , only : use_vsfm
     !
     implicit none
     type(bounds_type)     , intent(in)    :: bounds    ! bounds
     integer               , intent(in)    :: num_nolakec
     integer               , intent(in)    :: filter_nolakec(:)
     type(soilstate_type)  , intent(in)    :: soilstate_vars
     real(r8)              , intent(inout) :: soilbeta(bounds%begc:bounds%endc)

     !local variables
     real(r8) :: fac, fac_fc, wx      !temporary variables
     integer  :: c, l, fc     !indices


     associate(                                              &
          watsat      =>    soilstate_vars%watsat_col      , & ! Input:  [real(r8) (:,:)] volumetric soil water at saturation (porosity)
          watfc       =>    soilstate_vars%watfc_col       , & ! Input:  [real(r8) (:,:)] volumetric soil water at field capacity
          watmin      =>    soilstate_vars%watmin_col      , & ! Input:  [real(r8) (:,:)] min volumetric soil water
          sucmin      =>    soilstate_vars%sucmin_col      , & ! Input:  [real(r8) (:,:)] min volumetric soil water
          soilp_col   =>    col_ws%soilp      , & ! Input:  [real(r8) (:,:)] soil water pressure (Pa)

          h2osoi_ice  =>    col_ws%h2osoi_ice , & ! Input:  [real(r8) (:,:)] ice lens (kg/m2)
          h2osoi_liq  =>    col_ws%h2osoi_liq , & ! Input:  [real(r8) (:,:)] liquid water (kg/m2)
          frac_sno    =>    col_ws%frac_sno   , & ! Input:  [real(r8) (:)] fraction of ground covered by snow (0 to 1)
          frac_h2osfc =>    col_ws%frac_h2osfc  & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
          )

       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = col_pp%landunit(c)
          if (lun_pp%itype(l)/=istwet .AND. lun_pp%itype(l)/=istice  &
               .AND. lun_pp%itype(l)/=istice_mec) then
             if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
                wx   = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/col_pp%dz(c,1)
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
                if ( use_vsfm ) then
                   if ((wx < watmin(c,1)) .or. (soilp_col(c,1) < sucmin(c,1))) then
                      soilbeta(c) = 0._r8
                   end if
                end if
             else if (col_pp%itype(c) == icol_road_perv) then
                if (.not. use_vsfm) then
                   soilbeta(c) = 0._r8
                else
                   wx   = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/col_pp%dz(c,1)
                   fac  = min(1._r8, wx/watsat(c,1))
                   fac  = max( fac, 0.01_r8 )
                   if (wx < watfc(c,1) ) then  !when water content of ths top layer is less than that at F.C.
                      if (wx >= watmin(c,1) .and. soilp_col(c,1) >= sucmin(c,1) ) then
                         fac_fc  = min(1._r8, wx/watfc(c,1))  !eqn5.66 but divided by theta at field capacity
                         fac_fc  = max( fac_fc, 0.01_r8 )
                         ! modify soil beta by snow cover. soilbeta for snow surface is one
                         soilbeta(c) = (1._r8-frac_sno(c)-frac_h2osfc(c)) &
                              *0.25_r8*(1._r8 - cos(SHR_CONST_PI*fac_fc))**2._r8 &
                              + frac_sno(c)+ frac_h2osfc(c)
                      else
                         soilbeta(c) = 0._r8
                      endif
                   else   !when water content of ths top layer is more than that at F.C.
                      soilbeta(c) = 1._r8
                   end if
                endif
             else if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall) then
                soilbeta(c) = 0._r8
             else if (col_pp%itype(c) == icol_roof .or. col_pp%itype(c) == icol_road_imperv) then
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
     !$acc routine seq
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

    !------------------------------------------------------------------------------

   function getlblcef(rho,temp)result(cc)
     !$acc routine seq
     !compute the scaling paramter for laminar boundary resistance
     !the laminar boundary layer resistance is formulated as
     !Rb=2/(k*ustar)*(Sci/Pr)^(2/3)
     !cc = Rb*ustar
     !   = 2/k*(Sci/Pr)^(2/3)
     !   Pr=0.72, Prandtl number
     !   Sci = v/Di, Di is diffusivity of gas i
     !   v : kinetic viscosity

     use elm_varcon         , only :  vkc

     real(r8), intent(in) :: rho
     real(r8), intent(in) :: temp              ! air temperature

     real(r8), parameter  :: C = 120._r8       ! K
     real(r8), parameter  :: T0 = 291.25_r8    ! K
     real(r8), parameter  :: mu0 = 18.27e-6_r8 ! Pa s
     real(r8), parameter  :: prandtl = 0.72    !
     real(r8)             :: mu, diffh2o, sc
     real(r8)             :: cc                ! unitless scaling factor

     !compute the kinetic viscosity
     mu      = mu0 * (T0+C)/(temp+C) * (temp/T0)**(1.5)/rho !m^2 s^-1
     diffh2o = 0.229e-4_r8*(temp/273.15_r8)**1.75_r8        !m^2 s^-1
     sc      = mu/diffh2o                                   !schmidt number

     cc      = 2._r8/vkc*(Sc/Prandtl)**(2._r8/3._r8)

     return
   end function getlblcef

!------------------------------------------------------------------------------   
   subroutine calc_soil_resistance_sl14(bounds, num_nolakec, filter_nolakec, &
        soilstate_inst, waterstate_inst, temperature_inst, dsl, soilresis)
     !
     ! DESCRIPTION
     ! compute the lee-pielke beta factor to scal actual soil evaporation from
     ! potential evaporation
     ! xueyanz moved sl_14 (Swenson and Lawrence, 2014) from CLM5
     !
     ! USES
     use shr_kind_mod    , only : r8 => shr_kind_r8     
     use shr_const_mod   , only : SHR_CONST_PI
     use shr_log_mod     , only : errMsg => shr_log_errMsg   
     use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
     use decompMod       , only : bounds_type
     use clm_varcon      , only : denh2o, denice
     use landunit_varcon , only : istice_mec, istwet, istsoil, istcrop
     use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
     use column_varcon   , only : icol_road_imperv, icol_road_perv
     use ColumnType      , only : col_pp
     use LandunitType    , only : lun_pp
     !
     implicit none
     type(bounds_type)     , intent(in)    :: bounds    ! bounds   
     integer               , intent(in)    :: num_nolakec
     integer               , intent(in)    :: filter_nolakec(:)
     type(soilstate_type)  , intent(in)    :: soilstate_inst
     type(waterstate_type) , intent(in)    :: waterstate_inst
     type(temperature_type), intent(in)    :: temperature_inst
     real(r8)              , intent(inout) :: dsl(bounds%begc:bounds%endc)
     real(r8)              , intent(inout) :: soilresis(bounds%begc:bounds%endc)

   !local variables
     real(r8) :: aird, eps, dg, d0, vwc_liq
     real(r8) :: eff_por_top
     integer  :: c, l, fc     !indices
     
     SHR_ASSERT_ALL((ubound(dsl)    == (/bounds%endc/)), errMsg(sourcefile,
__LINE__))
     SHR_ASSERT_ALL((ubound(soilresis)    == (/bounds%endc/)),
errMsg(sourcefile, __LINE__))

     associate(                                              &
          dz                =>    col_pp%dz                             , & !
Input:  [real(r8) (:,:) ]  layer thickness (m)                             
          watsat            =>    soilstate_inst%watsat_col      , & ! Input:
[real(r8) (:,:)] volumetric soil water at saturation (porosity)
          bsw               =>    soilstate_inst%bsw_col             , & !
Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
          sucsat            =>    soilstate_inst%sucsat_col          , & !
Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
!          eff_porosity      =>    soilstate_inst%eff_porosity_col    , & !
!          Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice         
          t_soisno          =>    temperature_inst%t_soisno_col      ,  & !
Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                       
         
          h2osoi_ice        =>    waterstate_inst%h2osoi_ice_col , & ! Input:
[real(r8) (:,:)] ice lens (kg/m2)                       
          h2osoi_liq        =>    waterstate_inst%h2osoi_liq_col  & ! Input:
[real(r8) (:,:)] liquid water (kg/m2)                   
          )

!xueyanz
!open(12,file='t.txt')
!open(13,file='eff_por_top.txt')
!open(14,file='vwc_liq.txt')
!open(15,file='aird.txt')
!open(16,file='sucsat.txt')
!open(17,file='bsw.txt')
!open(18,file='watsat.txt')
   do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      l = col_pp%landunit(c)  
      if (lun_pp%itype(l)/=istwet .AND. lun_pp%itype(l)/=istice_mec) then
         if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
            vwc_liq = max(h2osoi_liq(c,1),1.0e-6_r8)/(dz(c,1)*denh2o)
! eff_porosity not calculated til SoilHydrology
             eff_por_top = max(0.01_r8,watsat(c,1)-min(watsat(c,1),
h2osoi_ice(c,1)/(dz(c,1)*denice)))

! calculate diffusivity and air free pore space
            aird = watsat(c,1)*(sucsat(c,1)/1.e7_r8)**(1./bsw(c,1))
            d0 = 2.12e-5*(t_soisno(c,1)/273.15)**1.75 ![Bitelli et al., JH, 08]
            eps = watsat(c,1) - aird
            dg = eps*d0*(eps/watsat(c,1))**(3._r8/max(3._r8,bsw(c,1)))
            
!      dsl(c) = dzmm(c,1)*max(0.001_r8,(0.8*eff_porosity(c,1) - vwc_liq)) &
! try arbitrary scaling (not top layer thickness)
!            dsl(c) = 15._r8*max(0.001_r8,(0.8*eff_porosity(c,1) - vwc_liq)) &
            dsl(c) = 15._r8*max(0.001_r8,(0.8*eff_por_top - vwc_liq)) &
                 !           /max(0.001_r8,(watsat(c,1)- aird))
                 /max(0.001_r8,(0.8*watsat(c,1)- aird))
            
            dsl(c)=max(dsl(c),0._r8)
            dsl(c)=min(dsl(c),200._r8)
            
            soilresis(c) = dsl(c)/(dg*eps*1.e3) + 20._r8
            soilresis(c) = min(1.e6_r8,soilresis(c))
!            write(12,*) dg*eps/d0 !xueyanz
!            write(13,*) eff_por_top
!            write(14,*) vwc_liq
!            write(15,*) aird
!            write(16,*) sucsat(c,1)
!            write(17,*) bsw(c,1)
!            write(18,*) watsat(c,1)
         else if (col_pp%itype(c) == icol_road_perv) then
            soilresis(c) = 1.e6_r8
         else if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) ==
icol_shadewall) then
            soilresis(c) = 1.e6_r8          
         else if (col_pp%itype(c) == icol_roof .or. col_pp%itype(c) ==
icol_road_imperv) then
            soilresis(c) = 1.e6_r8
         endif   
      else
         soilresis(c) =   0._r8
      endif
   enddo   
   end associate
   end subroutine calc_soil_resistance_sl14

!------------------------------------------------------------------------------   
   function do_soil_resistance_sl14()result(lres)
     !
     !DESCRIPTION
     ! return true if the soil evaporative resistance is computed using a DSL
     ! otherwise false
     implicit none
     logical :: lres

     if(soil_stress_method==sl_14)then
        lres=.true.
     else
        lres=.false.
     endif
     return

   end function do_soil_resistance_sl14

!------------------------------------------------------------------------------   
   subroutine calc_soil_resistance_sz09(bounds, num_nolakec, filter_nolakec, &
        soilstate_inst, waterstate_inst,dsl, soilresis)
     !
     ! DESCRIPTION
     ! compute the lee-pielke beta factor to scal actual soil evaporation from
     ! potential evaporation
     ! xueyanz implemented sz_09 (Sakaguchi and Zeng, 2009) DSL scheme here 
     !
     ! USES
     use shr_kind_mod    , only : r8 => shr_kind_r8     
     use shr_const_mod   , only : SHR_CONST_PI
     use shr_log_mod     , only : errMsg => shr_log_errMsg   
     use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
     use decompMod       , only : bounds_type
     use clm_varcon      , only : denh2o, denice
     use landunit_varcon , only : istice_mec, istwet, istsoil, istcrop
     use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
     use column_varcon   , only : icol_road_imperv, icol_road_perv
     use ColumnType      , only : col_pp
     use LandunitType    , only : lun_pp
     !
     implicit none
     type(bounds_type)     , intent(in)    :: bounds    ! bounds   
     integer               , intent(in)    :: num_nolakec
     integer               , intent(in)    :: filter_nolakec(:)
     type(soilstate_type)  , intent(in)    :: soilstate_inst
     type(waterstate_type) , intent(in)    :: waterstate_inst
!     type(temperature_type), intent(in)    :: temperature_inst
     real(r8)              , intent(inout) :: dsl(bounds%begc:bounds%endc)
     real(r8)              , intent(inout) :: soilresis(bounds%begc:bounds%endc)

   !local variables
     real(r8) :: dg, vwc_liq
     integer  :: c, l, fc     !indices
     
     SHR_ASSERT_ALL((ubound(dsl)    == (/bounds%endc/)), errMsg(sourcefile,
__LINE__))
     SHR_ASSERT_ALL((ubound(soilresis)    == (/bounds%endc/)),
errMsg(sourcefile, __LINE__))

     associate(                                              &
          dz                =>    col_pp%dz                             , & !
Input:  [real(r8) (:,:) ]  layer thickness (m)                             
          watsat            =>    soilstate_inst%watsat_col      , & ! Input:
[real(r8) (:,:)] volumetric soil water at saturation (porosity)
          watmin            =>    soilstate_inst%watmin_col      , & ! col
minimum volumetric soil water (nlevsoi)
          bsw               =>    soilstate_inst%bsw_col             , & !
Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
          h2osoi_ice        =>    waterstate_inst%h2osoi_ice_col , & ! Input:
[real(r8) (:,:)] ice lens (kg/m2)                       
          h2osoi_liq        =>    waterstate_inst%h2osoi_liq_col  & ! Input:
[real(r8) (:,:)] liquid water (kg/m2)                   
          )

!xueyanz
!open(14,file='vwc_liq.txt')
!open(17,file='bsw.txt')
!open(18,file='watsat.txt')
!open(19,file='watmin.txt')
   do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      l = col_pp%landunit(c)  
      if (lun_pp%itype(l)/=istwet .AND. lun_pp%itype(l)/=istice_mec) then
         if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
            vwc_liq = max(h2osoi_liq(c,1),1.0e-6_r8)/(dz(c,1)*denh2o)

! calculate diffusivity and air free pore space

            dg =
2.2e-5*((watsat(c,1))**2.)*((1.-watmin(c,1)/watsat(c,1))**(2.+3*bsw(c,1)))
            
!      dsl(c) = dzmm(c,1)*max(0.001_r8,(0.8*eff_porosity(c,1) - vwc_liq)) &
! try arbitrary scaling (not top layer thickness)
!            dsl(c) = 15._r8* (exp((1-min(1.,vwc_liq/watsat(c,1)))**5._r8)-1.) &
!                 /(2.71828-1.)
             dsl(c) = 50._r8* (exp((1-min(1.,vwc_liq/watsat(c,1)))**2.5_r8)-1.)
&
                      /(2.71828-1.)          !xueyanz 
            dsl(c)=max(dsl(c),0._r8)
            dsl(c)=min(dsl(c),200._r8)
            
            soilresis(c) = dsl(c)/(dg*1.e3)
            soilresis(c) = min(1.e6_r8,soilresis(c))
!            write(14,*) vwc_liq
!            write(17,*) bsw(c,1)
!            write(18,*) watsat(c,1)
!            write(19,*) watmin(c,1)
         else if (col_pp%itype(c) == icol_road_perv) then
            soilresis(c) = 1.e6_r8
         else if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) ==
icol_shadewall) then
            soilresis(c) = 1.e6_r8          
         else if (col_pp%itype(c) == icol_roof .or. col_pp%itype(c) ==
icol_road_imperv) then
            soilresis(c) = 1.e6_r8
         endif   
      else
         soilresis(c) =   0._r8
      endif
   enddo   
   end associate
   end subroutine calc_soil_resistance_sz09
!------------------------------------------------------------------------------  
 
   function do_soil_resistance_sz09()result(lres)
     !
     !DESCRIPTION
     ! return true if the soil evaporative resistance is computed using a DSL
     ! otherwise false
     implicit none
     logical :: lres

     if(soil_stress_method==sz_09)then
        lres=.true.
     else
        lres=.false.
     endif
     return

   end function do_soil_resistance_sz09
!------------------------------------------------------------------------------ 

end module SurfaceResistanceMod
