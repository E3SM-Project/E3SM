module TracerParamSetMod
#include "bshr_assert.h"

  use bshr_kind_mod            , only : r8 => shr_kind_r8
  use BeTR_decompMod           , only : bounds_type  => betr_bounds_type
 implicit none
 private
 public :: get_lgsorb_KL_Xsat, get_lnsorb_Kd
 public :: get_henrycef, get_equilibrium_scal
 public :: get_tauliq, get_taugas
 public :: get_aqueous_diffusivity, get_gas_diffusivity
 public :: get_diffusivity_ratio_gas2h2o
 public :: rhosat
contains

  subroutine get_lgsorb_KL_Xsat(tracerfamilyname, isoilorder, KL, Xsat)

  implicit none
  character(len=*), intent(in) :: tracerfamilyname
  integer , intent(in) :: isoilorder
  real(r8), intent(out) :: KL
  real(r8), intent(out) :: Xsat

  KL =1._r8
  Xsat = 0._r8


  end subroutine get_lgsorb_KL_Xsat



!-------------------------------------------------------------------------------
  function get_lnsorb_Kd(tracerfamily)result(Kd)

  !
  ! DESCRIPTION
  ! get the linear distribution parameter for modeling sorption
  implicit none
  character(len=*), intent(in) :: tracerfamily
  real(r8) :: Kd

  select case(trim(tracerfamily))
  case ('DOM')
    Kd = 1._r8
  case ('DOC')
    Kd = 1._r8
  case ('P_SOL')
    Kd = 1._r8
  case ('NH3x')
    Kd = 1._r8
  case default
    Kd = 0._r8
  end select
  end function get_lnsorb_Kd




!-------------------------------------------------------------------------------
   function get_diffusivity_ratio_gas2h2o(trcid, temp, betrtracer_vars)result(ratio)
   !
   ! DESCRIPTIONS:
   ! Compute the ratio of gas phase diffusivities for different volatile species in air with respect to that of water vapor
   !
   ! USES
   !
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   integer               , intent(in) :: trcid
   real(r8)              , intent(in) :: temp
   type(betrtracer_type) , intent(in) :: betrtracer_vars                                          ! betr configuration information

   real(r8) :: diff, ratio, diffh2o
   character(len=255) :: subname = 'get_diffusivity_ratio_gas2h2o'

   diffh2o=0.229e-4_r8*(temp/273.15_r8)**1.75_r8

   if(trcid==betrtracer_vars%id_trc_n2)then
      ratio=1.93e-5_r8*(temp/273.0_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_o2)then
      ratio=1.8e-5_r8*(temp/273.0_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_ar)then
      ratio=1.6e-5_r8*(temp/273.0_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_co2x)then
      ratio=1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_air_co2x)then
      ratio=1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_arrt_co2x)then
      ratio=1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_hrsoi_co2x)then
      ratio=1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_ch4)then
      ratio=1.9e-5_r8*(temp/298.0_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_nh3x)then
      ratio=0.211e-4_r8*(temp/273.15_r8)**1.75_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_no)then
      ratio=0.199e-4_r8*(temp/293.15)**1.5_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_n2o)then
      ratio=0.159e-4_r8*(temp/293.15)**1.5_r8/diffh2o
!isotopes
   elseif(trcid==betrtracer_vars%id_trc_o18_h2o)then
      ratio = 0.9691_r8            !0.9723_r8  from Merlivat, 1978, 0.9691 from Cappa et al, 2003
   elseif(trcid==betrtracer_vars%id_trc_d_h2o)then
      ratio = 0.9839_r8            !0.9755 from Merlivat, 1978, 0.9839 from Cappa et al., 2003
   elseif(trcid==betrtracer_vars%id_trc_c13_co2x)then
      ratio = 0.9957_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_c14_co2x)then
      ratio = 0.9913_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_o18_co2x)then
      ratio = 0.9913_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   elseif(trcid==betrtracer_vars%id_trc_o17_co2x)then
      ratio = 0.9957_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8/diffh2o
   endif
   return
   end function get_diffusivity_ratio_gas2h2o

!-------------------------------------------------------------------------------

   function get_gas_diffusivity(trcid, temp, betrtracer_vars)result(diff)
   !
   ! !DESCRIPTIONS
   !
   !compute gaseous diffusivity for volatile species
   !
   ! USES
   !
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   integer               , intent(in) :: trcid
   real(r8)              , intent(in) :: temp  !kelvin
   type(betrtracer_type) , intent(in) :: betrtracer_vars                                          ! betr configuration information

   !local variable
   real(r8) :: diff
   character(len=255) :: subname = 'get_gas_diffusivity'

   if(trcid==betrtracer_vars%id_trc_n2)then
      diff=1.93e-5_r8*(temp/273.0_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_ar)then
      diff=1.61e-5_r8*(temp/273.0_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_o2)then
      diff=1.8e-5_r8*(temp/273.0_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_co2x)then
      diff=1.47e-5_r8*(temp/273.15_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_air_co2x)then
      diff=1.47e-5_r8*(temp/273.15_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_arrt_co2x)then
      diff=1.47e-5_r8*(temp/273.15_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_hrsoi_co2x)then
      diff=1.47e-5_r8*(temp/273.15_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_ch4)then
      diff=1.9e-5_r8*(temp/298.0_r8)**1.82_r8
   elseif(trcid==betrtracer_vars%id_trc_nh3x)then
      diff=0.211e-4_r8*(temp/273.15_r8)**1.75_r8
   elseif(trcid==betrtracer_vars%id_trc_no)then
      diff=0.199e-4_r8*(temp/293.15)**1.5_r8
   elseif(trcid==betrtracer_vars%id_trc_n2o)then
      diff=0.159e-4_r8*(temp/293.15)**1.5_r8
!isotopes
   !use the kinetic theory of gases, assuming the collison diameters are same
   !between the light and heavy isotopmer, assuming the bath gas is dry air
   elseif(trcid==betrtracer_vars%id_trc_o18_h2o)then
      diff = 0.9723_r8*0.226e-4_r8*(temp/273.15_r8)**1.75_r8     !from Merlivat, 1978
   elseif(trcid==betrtracer_vars%id_trc_d_h2o)then
      diff = 0.9755_r8*0.226e-4_r8*(temp/273.15_r8)**1.75_r8     !from Merlivat, 1978
   elseif(trcid==betrtracer_vars%id_trc_c13_co2x)then
      diff = 0.9958_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8      ! from Trudinger, 1997
   elseif(trcid==betrtracer_vars%id_trc_o18_co2x)then
      diff = 0.9913_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8      !e.g. Wingate et al., 2009
   elseif(trcid==betrtracer_vars%id_trc_c14_co2x)then
      diff = 0.9918_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8      ! from Trudinger, 1997
   elseif(trcid==betrtracer_vars%id_trc_o17_co2x)then
      diff = 0.9957_r8*1.47e-5_r8*(temp/273.15_r8)**1.82_r8
   endif
   return
   end function get_gas_diffusivity

!-------------------------------------------------------------------------------

   function get_aqueous_diffusivity(trcid, temp, betrtracer_vars, is_dom)result(diff)
   !
   ! Descriptions:
   ! Compute aqueous diffusivity for volatile species
   !
   ! USES
   !
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   integer               ,           intent(in) :: trcid
   real(r8)              ,           intent(in) :: temp
   type(betrtracer_type) ,           intent(in) :: betrtracer_vars                                          ! betr configuration information
   logical               , optional, intent(in) :: is_dom

   !local variable
   real(r8) :: diff       ! m2/s
   character(len=255) :: subname ='get_aqueous_diffusivity'
   if(present(is_dom))then
     if(is_dom)then
       diff=2.6e-9_r8
       return
     endif
   endif
   if(trcid==betrtracer_vars%id_trc_n2)then
      diff=2.57e-9_r8*(temp/273.0_r8)
   elseif(trcid==betrtracer_vars%id_trc_o2)then
      diff=2.4e-9_r8*temp/298.0_r8
   elseif(trcid==betrtracer_vars%id_trc_ar)then
      diff=2.15e-9_r8*temp/298.0_r8
   elseif(trcid==betrtracer_vars%id_trc_co2x)then
      diff=1.81e-6_r8*exp(-2032.6/temp)
   elseif(trcid==betrtracer_vars%id_trc_air_co2x)then
      diff=1.81e-6_r8*exp(-2032.6/temp)
   elseif(trcid==betrtracer_vars%id_trc_arrt_co2x)then
      diff=1.81e-6_r8*exp(-2032.6/temp)
   elseif(trcid==betrtracer_vars%id_trc_hrsoi_co2x)then
      diff=1.81e-6_r8*exp(-2032.6/temp)
   elseif(trcid==betrtracer_vars%id_trc_ch4)then
      diff=1.5e-9_r8*temp/298.0_r8
   elseif(trcid==betrtracer_vars%id_trc_no3x)then
      diff=2.6e-9_r8*temp/298.15_r8
   elseif(trcid==betrtracer_vars%id_trc_no2x)then
      diff=2.6e-9_r8*temp/298.15_r8  !considers revision
   elseif(trcid==betrtracer_vars%id_trc_nh3x)then
      diff=1.64e-5_r8*temp/298.15_r8
!isotopes
   elseif(trcid==betrtracer_vars%id_trc_d_h2o)then
      diff=0.9833_r8*1e-9_r8*exp(-(535400._r8/temp-1393.3_r8)/temp+2.1876_r8)
   elseif(trcid==betrtracer_vars%id_trc_o18_h2o)then
      diff=0.9669_r8*1e-9_r8*exp(-(535400._r8/temp-1393.3_r8)/temp+2.1876_r8)
   elseif(trcid==betrtracer_vars%id_trc_o18_co2x)then
      !theoretical calculations based on molecular dynamics indicate the fractionation
      !between carbonate and bicarbonate due to diffusion is less than 1 per mil.
      !so set it to the diffusivity of the base value
      diff=1.81e-6_r8*exp(-2032.6/temp)
   elseif(trcid==betrtracer_vars%id_trc_c13_co2x)then
      diff=1.81e-6_r8*exp(-2032.6/temp)
   elseif(trcid==betrtracer_vars%id_trc_c14_co2x)then
      diff=1.81e-6_r8*exp(-2032.6/temp)
   else
      diff=2.6e-9_r8
   endif

   end function get_aqueous_diffusivity



!-------------------------------------------------------------------------------
   function get_taugas(eff_por, airvol, bsw)result(taugas)
   !
   !Descriptions
   !compute the tortuosity for the gas diffusion
   !Reference:
   !Millington and Quirk, 1961; Maggi et al., 2008
   !Moldrup et al, 2003
   !
   ! USES
   !

   implicit none
   real(r8),           intent(in) :: eff_por   ! effective porosity
   real(r8),           intent(in) :: airvol    ! air filled volume
   real(r8), optional, intent(in) :: bsw       ! clapp-hornber shape parameter

   real(r8) :: taugas
   character(len=255) :: subname ='get_taugas'

   if(eff_por == 0._r8)then
      taugas = 0._r8
   else
      if(present(bsw))then
         !modified from Eq.(5) in Moldrup et al., 2003
         taugas = (airvol/eff_por)**(3._r8/bsw)*airvol
      else
         taugas= eff_por**(1._r8/3._r8)*(airvol/eff_por)**(7._r8/3._r8)
      endif
   endif
   end function get_taugas


!-------------------------------------------------------------------------------
   function get_tauliq(eff_por, liqvol, bsw)result(tauliq)
   !
   !DESCRIPTION:
   !compute tortuosity for solute diffusion
   !Reference:
   !Millington and Quirk, 1961; Maggi et al., 2008
   !Moldrup et al, 2003
   !
   ! USES
   !
   implicit none
   real(r8),           intent(in) :: eff_por   !effective porosity
   real(r8),           intent(in) :: liqvol    !liquid water filled volume
   real(r8), optional, intent(in) :: bsw       !clapp-hornberg shape parameter

   real(r8) :: tauliq
   character(len=255) :: subname ='get_tauliq'

   if(eff_por == 0._r8)then
      tauliq = 0._r8
   else
      if(present(bsw))then
         !Modified from Eq.(3) in Moldrup et al., 2003.
         tauliq=(min(liqvol/eff_por,1._r8))**(bsw/3._r8-1._r8)*liqvol
      else
         tauliq=eff_por**(1._r8/3._r8)*(min(liqvol/eff_por,1._r8))**(7._r8/3._r8)
      endif
   endif
   end function get_tauliq


!-------------------------------------------------------------------------------
   function get_equilibrium_scal(temp, pH, tracer, betrtracer_vars)result(rscal)

   !DESCRIPTION:
   !
   !obtain the equilibrium scaling factor for species that
   !can exist in multiple aqueous forms through hydrolysis.
   !obtain mole fractions stored as carbonate and bicarbonate
   !dlogK/dT=dH/(2.303RT^2)
   !logK=logK(T_0)+dH*(1/(2.303R*T_0)-1/(2.303R*T))
   !
   ! Fang and Moncrieff, 1999, A model for soil CO2 production and transport 1:
   ! through the scale, the bulk aqueous tracer is, rscal * aqueous(netural)
   ! !USES
   !
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   real(r8)              , intent(in) :: temp, pH
   character(len=*)      , intent(in) :: tracer
   type(betrtracer_type) , intent(in) :: betrtracer_vars                                          ! betr configuration information

   !local variables
   !log value of dissociation constants, obtained under molar concentration (mol/dm3)
   real(r8), parameter :: Tref=298.15 ! Kelvin
   real(r8), parameter :: co2reflogK1=-6.352_r8   !25 Celcius
   real(r8), parameter :: co2reflogK2=-10.33_r8   !25 Celcisum
   real(r8), parameter :: co2dH1=-2.0e3_r8 ! J/mol
   real(r8), parameter :: co2dH2=-3.5e3_r8 ! J/mol
   real(r8), parameter :: nh3logK=9.24_r8  ! from Maggi et al. (2008)
   real(r8), parameter :: no3logK=1.30_r8  ! from Maggi et al. (2008)
   real(r8), parameter :: R=8.3144
   real(r8)            :: co2logK1, co2logK2, rscal
   character(len=255)  :: subname ='get_equilibrium_scal'

   if(trim(tracer)=='CO2x')then
      !H2CO3  <--> H(+)+HCO3(-)    K1
      !HCO3(-)<--> H(+)+CO3(2-)    K2
      !1.e3_r8 converts from mol/dm3 to mol/m3
      co2logK1 = co2reflogK1+co2dH1*(1._r8/(2.303_r8*R*Tref)-1._r8/(2.303_r8*R*temp))
      co2logK2 = co2reflogK2+co2dH2*(1._r8/(2.303_r8*R*Tref)-1._r8/(2.303_r8*R*temp))
      rscal = 1._r8+10._r8**(co2logK1)*10._r8**(-pH)*(1._r8+10._r8**(co2logK2)*10._r8**(-pH))*1.e3_r8

   elseif(trim(tracer)=='NH3x')then
      !NH3H2O <--> NH4(+) + OH(-)
      !10._r8**(-pH) gives mol / dm3
      rscal = 1._r8+10._r8**(nh3logK)*10._r8**(-pH)*1.e3_r8
   elseif(trim(tracer)=='NO3x')then
      !HNO3 <--> NO3(-) + H(+)
      rscal = 1._r8+10._r8**(no3logK)*10._r8**(-pH)*1.e3_r8
   else
      rscal = 1._r8  ! no rescal for other tracers
   endif
   return
   end function get_equilibrium_scal


!-------------------------------------------------------------------------------
   function get_henrycef(temp, trcid, betrtracer_vars)result(henry)
   !
   ! !DESCRIPTION
   !Compute the henry's law coefficient
   !There are unconsidered isotopic effect on the henry's law constant.
   !Some theoretical comments on such effect can be found in Jancso, 2002
   !REFERENCE: compilation of henry's law constants for inorganic and organic
   ! species of potential importance in environmental chemistry, Rolf Sander

   !USES
   use BeTRTracerType     , only : betrtracer_type
   implicit none
   real(r8)              , intent(in) :: temp
   integer               , intent(in) :: trcid
   type(betrtracer_type) , intent(in) :: betrtracer_vars                                          ! betr configuration information

   !local variable
   real(r8) :: henry    !unit[M/atm]=[mol_aq/dm3_aq]/[atm]
   real(r8) :: es
   character(len=255) :: subname ='get_henrycef'

   if(trcid == betrtracer_vars%id_trc_ar)then
      henry=1.4e-3_r8*exp(-1500._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_o2)then
      henry=1.3e-3_r8*exp(-1500._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid == betrtracer_vars%id_trc_nh3x)then
      henry=5.6e1_r8*exp(-4100._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_n2)then
      henry=6.1e-4_r8*exp(-1300._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_n2o)then
      henry=2.5e-2_r8*exp(-2600._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_no)then
      henry=1.9e-3*exp(-1500._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_c13_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_c14_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_o17_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_o18_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_air_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_arrt_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid==betrtracer_vars%id_trc_hrsoi_co2x)then
      henry=3.4e-2_r8*exp(-2400._r8*(1._r8/temp-1._r8/298.15_r8))
   else if(trcid == betrtracer_vars%id_trc_ch4)then
      henry=1.3e-3_r8*exp(-1700._r8*(1._r8/temp-1._r8/298.15_r8))  !mol dm^-3 amt^-1
   endif
   return
   end function get_henrycef


  !-----------------------------------------------------------------------
  subroutine rhoSat(T, rho, rhodT)
    ! compute the saturated vapor pressure density and its derivative against the temperature
    ! jyt
    use betr_varcon,    only: rwat => brwat
    use bshr_const_mod, only: SHR_CONST_TKFRZ

    implicit none
    real(r8)           , intent(in)  :: T
    real(r8)           , intent(out) :: rho
    real(r8), optional , intent(out) :: rhodT


    !------------------
    !local parameters
    ! For water vapor (temperature range 0C-100C)
    real(r8), parameter :: a0 =  6.11213476_r8
    real(r8), parameter :: a1 =  0.444007856_r8
    real(r8), parameter :: a2 =  0.143064234e-01_r8
    real(r8), parameter :: a3 =  0.264461437e-03_r8
    real(r8), parameter :: a4 =  0.305903558e-05_r8
    real(r8), parameter :: a5 =  0.196237241e-07_r8
    real(r8), parameter :: a6 =  0.892344772e-10_r8
    real(r8), parameter :: a7 = -0.373208410e-12_r8
    real(r8), parameter :: a8 =  0.209339997e-15_r8
    ! For derivative:water vapor
    real(r8), parameter :: b0 =  0.444017302_r8
    real(r8), parameter :: b1 =  0.286064092e-01_r8
    real(r8), parameter :: b2 =  0.794683137e-03_r8
    real(r8), parameter :: b3 =  0.121211669e-04_r8
    real(r8), parameter :: b4 =  0.103354611e-06_r8
    real(r8), parameter :: b5 =  0.404125005e-09_r8
    real(r8), parameter :: b6 = -0.788037859e-12_r8
    real(r8), parameter :: b7 = -0.114596802e-13_r8
    real(r8), parameter :: b8 =  0.381294516e-16_r8
    ! For ice (temperature range -75C-0C)
    real(r8), parameter :: c0 =  6.11123516_r8
    real(r8), parameter :: c1 =  0.503109514_r8
    real(r8), parameter :: c2 =  0.188369801e-01_r8
    real(r8), parameter :: c3 =  0.420547422e-03_r8
    real(r8), parameter :: c4 =  0.614396778e-05_r8
    real(r8), parameter :: c5 =  0.602780717e-07_r8
    real(r8), parameter :: c6 =  0.387940929e-09_r8
    real(r8), parameter :: c7 =  0.149436277e-11_r8
    real(r8), parameter :: c8 =  0.262655803e-14_r8
    ! For derivative:ice
    real(r8), parameter :: d0 =  0.503277922_r8
    real(r8), parameter :: d1 =  0.377289173e-01_r8
    real(r8), parameter :: d2 =  0.126801703e-02_r8
    real(r8), parameter :: d3 =  0.249468427e-04_r8
    real(r8), parameter :: d4 =  0.313703411e-06_r8
    real(r8), parameter :: d5 =  0.257180651e-08_r8
    real(r8), parameter :: d6 =  0.133268878e-10_r8
    real(r8), parameter :: d7 =  0.394116744e-13_r8
    real(r8), parameter :: d8 =  0.498070196e-16_r8

    real(r8) :: T_limit
    real(r8) :: td, es, esdT

    T_limit = T - SHR_CONST_TKFRZ
    if (T_limit > 100.0_r8) T_limit=100.0_r8
    if (T_limit < -75.0_r8) T_limit=-75.0_r8

    td       = T_limit
    if (td >= 0.0_r8) then
       es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
               + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
       esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
               + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
       es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
               + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
       esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
               + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    endif

    es    = es    * 100._r8            ! pa
    rho   = es/(rwat*T)                !kg  m^-3

    if(present(rhodT)) rhodT= esdT/(rwat*T)-rho/T         !kg  m^-3 K^-1

  end subroutine rhoSat


end module TracerParamSetMod
