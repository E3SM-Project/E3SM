!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_ConstantsMod

!BOP
! !MODULE: POP_ConstantsMod
!
! !DESCRIPTION:
!  This module defines a variety of physical and numerical constants
!  used throughout the Parallel Ocean Program.  This are initialized
!  either to internal constants or, if coupled to CCSM, to CCSM shr 
!  constants.
!
! !REVISION HISTORY:
!  SVN:$Id: $
!  CVS:$Name: POP_2_0_1 $

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_IOUnitsMod
   use netcdf
#ifdef CCSMCOUPLED
   use shr_const_mod
#endif

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_ConstantsInit

!  !PUBLIC DATA MEMBERS:

   ! numbers

   real (POP_r8), public :: &
      POP_pi,               &! pi
      POP_twoPi,            &! 2*pi
      POP_halfPi,           &! pi/2
      POP_tiny,             &! a very small number
      POP_big,              &! a very big number
      POP_undefinedR8        ! a value to use for undefined vars

   real (POP_r4), public :: &
      POP_undefinedR4        ! a value to use for undefined vars

   integer (POP_i4), public :: &
      POP_undefinedI4        ! a value to use for undefined vars

   character (POP_charLength), public ::  &
      POP_charBlank          ! empty character string

   ! physical constants
   ! note that most internal ocean constants are in cgs units
   !  while atmosphere and surface flux constants are sometimes
   !  in MKS units

   real (POP_r8), public ::   &
      POP_grav               ,&! gravit. accel. (cm/s^2)
      POP_omega              ,&! angular vel. of Earth 1/s
      POP_radiusEarth        ,&! radius of Earth (cm)
      POP_rhoAir             ,&! ambient air density (kg/m^3)
      POP_rhoFW              ,&! density of fresh water (g/cm^3)
      POP_rhoSW              ,&! density of salt water (g/cm^3)
      POP_CpSW               ,&! specific heat salt water
      POP_CpAir              ,&! heat capacity of air (J/kg/K)
      POP_speedSound         ,&! speed of sound (cm/s)
      POP_vonKarman          ,&! von Karman constant
      POP_emissivity         ,&!
      POP_stefanBoltzmann    ,&! W/m^2/K^4
      POP_latentHeatVapor    ,&! lat heat of vaporization (erg/g)
      POP_latentHeatVaporMKS ,&! lat heat of vaporization (J/kg)
      POP_latentHeatFusion   ,&! lat heat of fusion (erg/g)
      POP_latentHeatFusionMKS,&! lat heat of fusion (J/kg)
      POP_seaIceSalinity     ,&! salinity of sea ice formed (psu)
      POP_ocnRefSalinity       ! ocean reference salinity (psu)

   !  conversion factors

   real (POP_r8), public ::   &
      POP_degreeToRadian,     &! degree-radian conversion
      POP_radianToDegree,     &! degree-radian conversion
      POP_T0Kelvin,           &! zero point for Celcius
      POP_MeterPerCM,         &! meters per cm
      POP_cmPerMeter,         &! cm per meter
      POP_saltToPpt,          &! salt (g/g) to ppt
      POP_pptToSalt,          &! salt ppt to g/g
      POP_massToSv,           &! mass flux to Sverdrups
      POP_heatToPW,           &! heat flux to Petawatts
      POP_saltToSvppt,        &! salt flux to Sv*ppt
      POP_saltToMmPerDay,     &! salt to water (mm/day)
      POP_momentumFactor,     &! wind stress (N/m^2) to vel flux (cm^2/s^2)
      POP_heatFluxFactor,     &! heat flux (W/m^2) to temp flux (C*cm/s)
      POP_FWFluxFactor,       &! fw flux (kg/m^2/s) to salt flux (msu*cm/s)
      POP_SaltFluxFactor,     &! salt flux (kg/m^2/s) to salt flux (msu*cm/s)
      POP_FWMassToFWFlux       ! fw flux (kg/m^2/s) to fw flux (cm/s)

#ifdef CCSMCOUPLED
   real(POP_r8), public ::  &
      T0_Kelvin                ! zero point for Celsius
#endif
!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_ConstantsInit
! !INTERFACE:

 subroutine POP_ConstantsInit(errorCode)

! !DESCRIPTION:
!  This subroutine initializes standard constants using either
!  internal values or other accepted values.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode             ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------

   integer (POP_i4) :: n   ! dummy loop index

   character (POP_charLength) :: &
      outFormat            ! output format for writing constants

!-----------------------------------------------------------------------
!
!  define numbers and character constants
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   POP_pi          = 4.0_POP_r8*atan(1.0_POP_r8)
#ifdef COUPLEDCCSM
   POP_pi          = SHR_CONST_PI
#endif
   POP_twoPi       = 2.0_POP_r8*POP_pi
   POP_halfPi      = 0.5_POP_r8*POP_pi
   POP_tiny        = 1.0e-20_POP_r8
   POP_big         = 1.0e+30_POP_r8
   POP_undefinedR8 = nf90_fill_double
   POP_undefinedR4 = nf90_fill_real
   POP_undefinedI4 = nf90_fill_int

   do n=1,POP_charLength
     POP_charBlank(n:n) = ' '
   end do


!-----------------------------------------------------------------------
!
!  physical constants
!  note that most internal ocean constants are in cgs units
!   while atmosphere and surface flux constants are sometimes
!   in MKS units
!
!  some of these constants may be over-ridden by CSM-defined
!  constants if the CSM shared constants are available
!
!-----------------------------------------------------------------------

   POP_grav                = 980.6_POP_r8 ! gravit. accel. (cm/s^2)
   POP_omega               = 7.292123625e-5_POP_r8  ! angular vel. Earth 1/s
   POP_radiusEarth         = 6370.0e5_POP_r8        ! radius of Earth (cm)
   POP_rhoAir              = 1.2_POP_r8   ! ambient air density (kg/m^3)
   POP_rhoFW               = 1.0_POP_r8   ! avg. water density (g/cm^3)
   POP_rhoSW               = 4.1_POP_r8/3.996_POP_r8  ! density of salt water (g/cm^3)
   POP_CpSW                = 3.996e7_POP_r8  ! specific heat salt water
   POP_CpAir               = 1005.0_POP_r8   ! heat capacity of air (J/kg/K)
   POP_speedSound          = 1.5e5_POP_r8    ! speed of sound (cm/s)
   POP_vonKarman           = 0.4_POP_r8      ! von Karman constant
   POP_emissivity          = 1.0_POP_r8      !
   POP_stefanBoltzmann     = 567.0e-10_POP_r8 !  W/m^2/K^4
   POP_latentHeatVaporMKS  = 2.5e6_POP_r8    ! lat heat of vaporization (J/kg)
   POP_latentHeatFusion    = 3.34e9_POP_r8   ! lat heat of fusion (erg/g)
   POP_latentHeatFusionMKS = 3.34e5_POP_r8   ! lat heat of fusion (J/kg)
   POP_seaIceSalinity      = 4.0_POP_r8      ! (psu)
   POP_ocnRefSalinity      = 34.7_POP_r8     ! (psu)

#ifdef CCSMCOUPLED
   POP_grav                = SHR_CONST_G*100.0_POP_r8  ! cm/s^2
   POP_omega               = SHR_CONST_OMEGA         ! rad/s
   POP_radiusEarth         = SHR_CONST_REARTH*100.0_POP_r8 ! cm
   POP_rhoAir              = SHR_CONST_RHODAIR       ! kg/m^3
   POP_rhoFW               = SHR_CONST_RHOFW*0.001_POP_r8  ! g/cm^3
   POP_rhoSW               = SHR_CONST_RHOSW*0.001_POP_r8  ! g/cm^3
   POP_CpSW                = SHR_CONST_CPSW*10000.0_POP_r8 ! erg/g/K
   POP_CpAir               = SHR_CONST_CPDAIR        ! J/kg/K
   POP_vonKarman           = SHR_CONST_KARMAN
   POP_stefanBoltzmann     = SHR_CONST_STEBOL        ! W/m^2/K^4
   POP_latentHeatVaporMKS  = SHR_CONST_LATVAP        ! J/kg
   POP_latentHeatFusion    = SHR_CONST_LATICE*10000.0_POP_r8 ! erg/g
   POP_latentHeatFusionMKS = SHR_CONST_LATICE        ! J/kg 
   POP_seaIceSalinity      = SHR_CONST_ICE_REF_SAL   ! psu
   POP_ocnRefSalinity      = SHR_CONST_OCN_REF_SAL   ! psu
#endif

#ifdef ZERO_SEA_ICE_REF_SAL
    POP_seaIceSalinity       = 0.0_POP_r8
#endif


!-----------------------------------------------------------------------
!
!  conversion factors
!
!-----------------------------------------------------------------------

   POP_degreeToRadian = POP_pi/180.0_POP_r8
   POP_radianToDegree = 180.0_POP_r8/POP_pi
   POP_T0Kelvin       = 273.16_POP_r8       ! zero point for Celsius
#ifdef CCSMCOUPLED
   T0_Kelvin          = SHR_CONST_TKFRZ     ! zero point for Celsius
#endif
   POP_MeterPerCM     = .01_POP_r8          ! meters per cm
   POP_cmPerMeter     = 100._POP_r8         ! cm per meter
   POP_saltToPpt      = 1000._POP_r8        ! salt (g/g) to ppt
   POP_pptToSalt      = 1.e-3_POP_r8        ! salt ppt to g/g
   POP_massToSv       = 1.0e-12_POP_r8      ! mass flux to Sverdrups
   POP_heatToPW       = 4.186e-15_POP_r8    ! heat flux to Petawatts
   POP_saltToSvppt    = 1.0e-9_POP_r8       ! salt flux to Sv*ppt
   POP_saltToMmPerDay = 3.1536e+5_POP_r8    ! salt to water (mm/day)

!-----------------------------------------------------------------------
!
!  convert windstress (N/m^2) to velocity flux (cm^2/s^2):
!  -------------------------------------------------------
!    windstress in (N/m^2) = (kg/s^2/m) = 10(g/s^2/cm) = 10(dyn/cm^2)
!    assume here that density of seawater rho = 1 (g/cm^3)
!
!    vel_flux   = windstress / rho
!    vel_flux (cm^2/s^2) = windstress (N/m^2)*10 (g/s^2/cm)/(N/m^2)
!                          / [1 (g/cm^3)]
!                        = windstress (N/m^2)
!                          * momentum_factor ((cm^2/s^2)/N/m^2)
!    ==>  momentum_factor = 10
!
!-----------------------------------------------------------------------

   POP_momentumFactor = 10.0_POP_r8

!-----------------------------------------------------------------------
!
!  convert heat, solar flux (W/m^2) to temperature flux (C*cm/s):
!  --------------------------------------------------------------
!    heat_flux in (W/m^2) = (J/s/m^2) = 1000(g/s^3)
!    density of seawater rho_sw in (g/cm^3)
!    specific heat of seawater cp_sw in (erg/g/C) = (cm^2/s^2/C)
!
!    temp_flux          = heat_flux / (rho_sw*cp_sw)
!    temp_flux (C*cm/s) = heat_flux (W/m^2)
!                         * 1000 (g/s^3)/(W/m^2)
!                         / [(rho_sw*cp_sw) (g/cm/s^2/C)]
!
!                       = heat_flux (W/m^2)
!                         * hflux_factor (C*cm/s)/(W/m^2)
!
!    ==>  hflux_factor = 1000/(rho_sw*cp_sw)
!
!-----------------------------------------------------------------------

   POP_heatFluxFactor = 1000.0_POP_r8/(POP_rhoSW*POP_CpSW)

!-----------------------------------------------------------------------
!
!  convert fresh water flux (kg/m^2/s) to virtual salt flux (msu*cm/s):
!  --------------------------------------------------------------------
!    ocean reference salinity in (o/oo=psu)
!    density of freshwater rho_fw = 1.0 (g/cm^3)
!    h2o_flux in (kg/m^2/s) = 0.1 (g/cm^2/s)
!
!    salt_flux            = - h2o_flux * ocn_ref_salinity / rho_fw
!    salt_flux (msu*cm/s) = - h2o_flux (kg/m^2/s)
!                           * ocn_ref_salinity (psu)
!                           * 1.e-3 (msu/psu)
!                           * 0.1 (g/cm^2/s)/(kg/m^2/s)
!                           / 1.0 (g/cm^3)
!                         = - h2o_flux (kg/m^2/s)
!                           * ocn_ref_salinity (psu)
!                           * fwflux_factor (cm/s)(msu/psu)/(kg/m^2/s)
!
!    ==>  fwflux_factor = 1.e-4
!
!    salt_flux(msu*cm/s) = h2oflux(kg/m^2/s) * salinity_factor
!
!    ==> salinity_factor = - ocn_ref_salinity(psu) * fwflux_factor
!
!-----------------------------------------------------------------------

   POP_FWFluxFactor   = -POP_ocnRefSalinity*1.e-4_POP_r8

!-----------------------------------------------------------------------
!
!  convert salt flux (kg/m^2/s) to salt flux (msu*cm/s):
!  -----------------------------------------------------
!    density of freshwater rho_fw = 1.0 (g/cm^3)
!    salt_flux_kg in (kg/m^2/s) = 0.1 (g/cm^2/s)
!
!    salt_flux            = - h2o_flux * ocn_ref_salinity / rho_fw
!    salt_flux (msu*cm/s) = salt_flux_kg (kg/m^2/s)
!                           * 0.1 (g/cm^2/s)/(kg/m^2/s)
!                           / 1.0 (g/cm^3)
!                         = salt_flux_kg (kg/m^2/s)
!                           * sflux_factor (msu*cm/s)/(kg/m^2/s)
!
!    ==>  sflux_factor = 0.1
!
!-----------------------------------------------------------------------

   POP_SaltFluxFactor = 0.1_POP_r8

!-----------------------------------------------------------------------
!
!  convert fresh water mass flux (kg/m^2/s) to fresh water flux (cm/s):
!  --------------------------------------------------------------------
!    density of freshwater rho_fw = 1.0 (g/cm^3)
!    h2o_flux in (kg/m^2/s) = 0.1 (g/cm^2/s)
!
!    fw_flux  = h2o_flux / rho_fw
!    fw_flux (cm/s) = h2o_flux (kg/m^2/s)
!                     * 0.1 (g/cm^2/s)/(kg/m^2/s)
!                     / 1.0 (g/cm^3)
!                   = h2o_flux (kg/m^2/s)
!                     * fwmass_to_fwflux (cm/s)/(kg/m^2/s)
!
!    ==>  fwmass_to_fwflux = 0.1
!
!-----------------------------------------------------------------------

   POP_FWMassToFWFlux = 0.1_POP_r8

!-----------------------------------------------------------------------
!
!  Document important constants
!
!-----------------------------------------------------------------------

   if (POP_myTask == POP_masterTask) then
      write(POP_stdout,POP_blankFormat)
      write(POP_stdout,POP_delimFormat)
      write(POP_stdout,'(a30)') 'Physical constant values used:'
      write(POP_stdout,POP_delimFormat)
      write(POP_stdout,POP_blankFormat)

      outFormat = '(a34,1pe22.15)'

      write(POP_stdout,outFormat) ' Pi                             = ', &
                                  POP_pi
      write(POP_stdout,outFormat) ' Gravity                        = ', &
                                  POP_grav
      write(POP_stdout,outFormat) ' Omega                          = ', &
                                  POP_omega
      write(POP_stdout,outFormat) ' Gravity                        = ', &
                                  POP_grav
      write(POP_stdout,outFormat) ' Omega                          = ', &
                                  POP_omega
      write(POP_stdout,outFormat) ' Earth radius                   = ', &
                                  POP_radiusEarth
      write(POP_stdout,outFormat) ' Density air                    = ', &
                                  POP_rhoAir
      write(POP_stdout,outFormat) ' Density fresh water            = ', &
                                  POP_rhoFW
      write(POP_stdout,outFormat) ' Density salt  water            = ', &
                                  POP_rhoSW
      write(POP_stdout,outFormat) ' Spec. heat salt water          = ', &
                                  POP_CpSW
      write(POP_stdout,outFormat) ' Spec. heat air                 = ', &
                                  POP_CpAir
      write(POP_stdout,outFormat) ' Sound speed                    = ', &
                                  POP_speedSound
      write(POP_stdout,outFormat) ' von Karman constant            = ', &
                                  POP_vonKarman
      write(POP_stdout,outFormat) ' Emissivity                     = ', &
                                  POP_emissivity
      write(POP_stdout,outFormat) ' Stefan Boltzmann               = ', &
                                  POP_stefanBoltzmann
      write(POP_stdout,outFormat) ' Latent heat vaporization (MKS) = ', &
                                  POP_latentHeatVaporMKS
      write(POP_stdout,outFormat) ' Latent heat fusion             = ', &
                                  POP_latentHeatFusion
      write(POP_stdout,outFormat) ' Latent heat fusion (MKS)       = ', &
                                  POP_latentHeatFusionMKS
      write(POP_stdout,outFormat) ' Sea ice salinity               = ', &
                                  POP_seaIceSalinity
      write(POP_stdout,outFormat) ' Ocean ref. salinity            = ', &
                                  POP_ocnRefSalinity
   endif

!EOC
!-----------------------------------------------------------------------

 end subroutine POP_ConstantsInit

!***********************************************************************

 end module POP_ConstantsMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
