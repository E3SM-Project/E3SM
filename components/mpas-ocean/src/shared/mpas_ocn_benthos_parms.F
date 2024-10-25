! copyright (c) 2013,  los alamos national security, llc (lans)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_benthos_parms
!
!> \brief MPAS ocean benthos parameters
!> \author Nicole Jeffery
!> \date   07/21/2020
!> \details
!> This module contains parameters and type definitions for the
!> ocean benthos biogeochemical model.
!
!-----------------------------------------------------------------------

module ocn_benthos_parms

  implicit none

  private
  save

  !--------------------------------------------------------------------
  !
  ! Public member functions
  !
  !--------------------------------------------------------------------

  public :: benthos_parms_init

  !-----------------------------------------------------------------------------
  !   kinds definitions
  !-----------------------------------------------------------------------------

  integer, parameter, public ::               &
      benthos_char           = 256                    ,&
      benthos_log            = kind(.true.)           ,&
      benthos_i4             = selected_int_kind(6)   ,&
      benthos_i8             = selected_int_kind(13)  ,&
      benthos_r4             = selected_real_kind(6)  ,&
      benthos_r8             = selected_real_kind(13)

  real (KIND=benthos_r8), public :: T0_Kelvin_benthos
  !-----------------------------------------------------------------------------
  !   floating point constants used across ecosystem module
  !-----------------------------------------------------------------------------

  real (KIND=benthos_r8), parameter, public :: &
       c0_benthos  =   0.0_benthos_r8,    &
       c1_benthos  =   1.0_benthos_r8,    &
       c2_benthos  =   2.0_benthos_r8,    &
       c10_benthos = 10.0_benthos_r8,     &
       p5_benthos  =  0.5_benthos_r8,     &
       p1_benthos  =  0.1_benthos_r8

  real(KIND=benthos_r8), parameter, public :: &
       puny          = 1.0e-15_benthos_r8, & ! small number
       puny16        = 1.0e-16_benthos_r8, & ! small number
       sec_per_hour  = 3600.0_benthos_r8,  & ! number of secondes per hour
       sec_per_day   = 86400.0_benthos_r8,  & ! number of seconds per day
       sec_per_year  = 31536000.0_benthos_r8, & ! number of seconds per year
       days_per_sec  = c1_benthos / sec_per_day,            & ! number of days in a second
       years_per_sec = c1_benthos / sec_per_year, & ! number of years in a second
       m2percm2      = 0.0001_benthos_r8, &      ! m2 per cm2
       mM_umolperL   = 1000.0_benthos_r8, &  ! mM to umol/L
       fluxSed_o     = 1.22_benthos_r8 / sec_per_year, & ! Prescribed sediment flux (kg/m2/s)
       deepFlux_h2s  = -50.0_benthos_r8 / sec_per_year, & ! upward flux for h2s (mmol/m2/s)
       pocFluxRate_o = 6000.0_benthos_r8 / sec_per_year, &   ! POC flux rate (mmol/m2/s)
       bottom_water_density = 1027.0_benthos_r8, & ! kg/m3
       sediment_density = 2650.0_benthos_r8, & ! (kg/m3)
       mass_to_vol   = 1.0e3_benthos_r8 * bottom_water_density, & ! mol/kg -> mmol/m3
!       mass_to_vol   = 1.0e6_benthos_r8 * sediment_density, &
       vol_to_mass   = c1_benthos / mass_to_vol, &   ! mmol/m3 -> mol/kg
       mmol_per_kg_water = 5.55556e4_benthos_r8, &     ! mmol of water/ kg of water
       degreeC_to_K = 273.15_benthos_r8  , &
       kg_per_mg = 1.0e-6_benthos_r8, & ! miligrams to kilograms conversion
       gSi_per_molSi = 28.1_benthos_r8, & !g/mol silica molecular weight
       gC_per_molC = 12.0107_benthos_r8, & !g/mol carbon molecular weight
       massOpal_per_massbSi = 2.1_benthos_r8, & ! kg Opal per kg biogenic silica; Fahl and Nothig, 2007
       mass_per_molPIC = 100.1_benthos_r8, & !g/mol CaCO3 molecular weight
       mass_per_massPOC = 2.0_benthos_r8, & ! mass per mass of particulate organic carbon; Fahl and Nothig, 2007
       deg_to_rad = 3.14159265359_benthos_r8/180.0_benthos_r8

  !-----------------------------------------------------------------------------
  !   Fixed biogeochemical parameters
  !-----------------------------------------------------------------------------

  integer (KIND=benthos_i4), parameter, public :: &
       nPrimaryReactions   = 6,  & ! number of primary reactions
       nSecondaryReactions = 19, & ! number of secondary reactions
       nCarbonateReactions = 3     ! number of carbonate dissolution reactions

  integer (KIND=benthos_i4), dimension(:), allocatable, public :: &
       transport_index  ! contains mole ratios of secondary reaction sources

  real(KIND=benthos_r8), dimension(:,:), allocatable, public :: &
       secondarySourceStoich, & ! contains mole ratios of secondary reaction sources
       secondarySinkStoich,   & ! contains mole ratios of secondary reaction sinks
       primarySourceStoich,   & ! contains mole ratios of primary reaction sources
       primarySinkStoich,     & ! contains mole ratios of primary reaction sinks
       carbonateSourceStoich, & ! contains mole ratios of carbonate reaction sources
       carbonateSinkStoich,    & ! contains mole ratios of carbonate reaction sinks
       elementRatios                  ! converts mol tracer to mol elements

  real(KIND=benthos_r8), parameter, public :: &
       benthosDepth     = 0.3_benthos_r8,   & ! (m) depth of the benthos active layer
       benthos_salt_min = 0.1_benthos_r8,   & ! minimum salinity (ppt)
       benthos_dic_min  = benthos_salt_min/35.0_benthos_r8 * 1944.0_benthos_r8, & ! minium dic (mmol/m3)
       benthos_alk_min  = benthos_salt_min/35.0_benthos_r8 * 2225.0_benthos_r8, & ! minium alkalinity (mmolEq/m3)
       oceanMagnesium = 53.0_benthos_r8/sediment_density, & ! (mol/kg)  concentration of Mg in the ocean
       CtoP = 106.0_benthos_r8, & ! carbon to phosphorus ratio of sinking organic matter  (117 in BGC)
       NtoP = 16.0_benthos_r8,  & ! carbon to nitrogen ratio of organic sinking matter)
       CtoP_r = 290.0_benthos_r8, & ! carbon to phosphorus ratio of sinking organic matter  (117 in BGC)
       NtoP_r = 29.0_benthos_r8,  & ! carbon to nitrogen ratio of organic sinking matter)
       ironBoundPFraction = 0.175_benthos_r8, & ! fraction of bound p to iron
       fracHighlyReactive = 0.2_benthos_r8, & ! fraction of POC flux that is highly reactive
       fracLessReactive = 0.24_benthos_r8, & ! Fraction of POC flux that is less reactive
       fracRefractory = 0.56_benthos_r8, &   ! Fraction of POC flux that is refractory
       fracCalcite = 0.50_benthos_r8 , &     ! Fraction of caco3 flux that is calcite
       fracAragonite = 0.35_benthos_r8, &    ! Fraction of caco3 flux that is aragonite
       fracMgCalcite = 0.15_benthos_r8, &    ! Fraction of caco3 flux that is 15% Mg calcite
       sedFlux1 = 0.87478367_benthos_r8, & ! Middelburg et al (1997) exponent parameter w1 for sediment flux
       sedFlux2 = 0.00042512_benthos_r8, &   ! as a function of ocean depth.  exponent parameter w2
       sedFlux3 = 0.033_benthos_r8*years_per_sec ! scale parameter wo ==> result is sediment flux in units of m/y
                                                 ! to m/s

  !-----------------------------------------------------------------------------
  !  Primary Rate Constants from Reed et al 2011

  real (KIND=benthos_r8), parameter, public :: &
      k_alpha = 1.62_benthos_r8/sec_per_year, &   ! per year  Labile rate constant for primary redox of OM^a
      k_beta = 0.0086_benthos_r8/sec_per_year, & ! per year  Semi-labile rate constant for primary redox of OM^b
      k_ref = c0_benthos, &   ! refractory
      k_o2 = 20.0_benthos_r8, & ! mmol/m3  half sat/inhibition for O2
      k_no3 = 4.0_benthos_r8, & ! mmol/m3 half sat/inhibition for NO3
      k_mno2a = 4.0_benthos_r8, & ! mmol/kg half sat/inhibition for MnO2 (umol/g)
      k_feoh3a = 65.0_benthos_r8, & ! mmol/kg half sat/inhibition for fe(OH)3
      k_so4 = 1600.0_benthos_r8, & ! mmol/m3 half sat/inhibition for SO4
      reductionFractionSO4 = 0.075_benthos_r8, & ! reduces redox rates for P5-P6
      epsilon_Phosphate = 1.0  ! accelerated phosphate remineralization (2.0 for Arkona Basin)
  
  ! Secondary Rate Constants from Reed et al 2011
  real(KIND=benthos_r8), parameter, public :: &       
      k_s1 = 10000.0_benthos_r8/sec_per_year/mM_umolperL, & !/mmol m3 per s  from /mmol L per year
      k_s2 = 20000.0_benthos_r8/sec_per_year/mM_umolperL, & !/mmol m3 per s  from per year  800-20,000 is the range
      k_s3 = 140000.0_benthos_r8/sec_per_year/mM_umolperL, & !/mmol m3 per s  from per year
      k_s4 = 300.0_benthos_r8/sec_per_year/mM_umolperL, & !/mmol m3 per s from /mmol L per year
      k_s5 = c1_benthos/sec_per_year/mM_umolperL, &  !/mmol m3 per s from /mmol L per year
      k_s6 = 160.0_benthos_r8/sec_per_year/mM_umolperL, & !/mmol m3 per s from /mmol L per year
      k_s7 = 10000000.0_benthos_r8/sec_per_year/mM_umolperL, & !/mmol m3 per s from /mmol L per year
      k_s8 = c2_benthos/sec_per_year/mM_umolperL, &
      k_s9 = 20.0_benthos_r8/sec_per_year/mM_umolperL, &
      k_s10 = 8.0_benthos_r8/sec_per_year/mM_umolperL, &
      k_s11 = 100.0_benthos_r8/sec_per_year/mM_umolperL, &
      k_s12 = c10_benthos/sec_per_year/mM_umolperL, &
      k_s13 = 3.0_benthos_r8/sec_per_year, &  !per year
      k_s14 = 7.0_benthos_r8/sec_per_year/mM_umolperL, &   !  why is this in per (mmol/m3 s)  Now g/umol/s
      k_s15 = 0.6_benthos_r8/sec_per_year, & !per year
      k_s16 = 1.8_benthos_r8/sec_per_year ! per year
  
  ! Carbonate Rate constants
  ! Constants from Walter and Morse (1985) reduced by inhibition factor of 0.1 as in Krumins et al

  real(KIND=benthos_r8), parameter, public :: &
       inhibition = 0.1_benthos_r8, & ! inhibition factor from Krumins et al
       kp1 = 14.296e-6_benthos_r8, & ! mmoles/mg/hr !14.296 umoles/g/hr (Walter and Morse (1985) table 4)barnacle Balanus
       kp2 = 22.198e-6_benthos_r8, & ! mmol/mg/hr !22.198 umoles/g/hr (Walter and Morse (1985) table 4)green algae Halimeda
       kp3 = 16.78e-6_benthos_r8,  & ! mmol/mg/hr !16.78 umoles/g/hr (Walter and Morse (1985table 4) foraminifera Peneroplis
       n1 = 2.74_benthos_r8,     & ! Walter and morse 1985 table 4; power law for barnacle Balanus
       n2 = 2.43_benthos_r8,     & ! Walter and morse 1985 table 4; power law for green algae Halimeda
       n3 = 3.51_benthos_r8,     & !  Walter and morse 1985 table 4; 3.51 power law for foraminifera Peneroplis
! old spin-up
!       kp1 = 457.1_benthos_r8,   & ! 14.296E-6 mmoles/mg/hr !14.296 umoles/g/hr (Walter and Morse (1985) table 4)barnacle Balanus
!       kp2 = 1258.9_benthos_r8,  & ! 22.198E-6 mmol/mg/hr !22.198 umoles/g/hr (Walter and Morse (1985) table 4)green algae Halimeda
!       kp3 = 660.69_benthos_r8,  & ! 16.78E-6 mmol/mg/hr !16.78 umoles/g/hr (Walter and Morse (1985table 4) foraminifera Peneroplis
!       n1 = 2.74_benthos_r8,     & ! Walter and morse 1985 table 4; power law for barnacle Balanus
!       n2 = 2.43_benthos_r8,     & ! Walter and morse 1985 table 4; power law for green algae Halimeda
!       n3 = 3.61_benthos_r8,     & !  Walter and morse 1985 table 4; 3.51 power law for foraminifera Peneroplis
       mm_calcite = 100.0869_benthos_r8,   & ! g/mol molar mass of calcite
       mm_aragonite = 100.0869_benthos_r8, & ! g/mol molar mass of aragonite
       mm_mgcalcite = 97.7201_benthos_r8     ! g/mol molar mass of 15% Mg + 85% Ca + CO3
  !-----------------------------------------------------------------------------
  !   functional group
  !-----------------------------------------------------------------------------

  type, public :: benthos_indices_type
   integer (KIND=benthos_i4) :: &
      poca_ind,          & ! particulate labile organic carbon
      pocb_ind,          & ! particulate semi-labile organic carbon
      pocc_ind,          & ! particulate refractory organic carbon
      pona_ind,          & ! particulate labile organic nitrogen
      ponb_ind,          & ! particulate semi-labile organic nitrogen
      ponc_ind,          & ! particulate refractory organic nitrogen
      popa_ind,          & ! particulate labile organic phophorus
      popb_ind,          & ! particulate semi-labile organic phophorus
      popc_ind,          & ! particulate refractory organic phophorus
      o2_ind,            & ! dissolved oxygen
      nh4_ind,           & ! dissolved ammonia
      h2po4_ind,         & ! dissolved inorganic phophate
      co2_ind,           & ! dissolved carbon dioxide
      no3_ind,           & ! dissolved nitrate
      mno2a_ind,         & ! particulate amorphous manganese oxyhydroxide
      mno2b_ind,         & ! particulate crystalline manganese oxyhydroxide
      mn_ind,            & ! dissolved manganese(II)
      feoh3a_ind,        & ! particulate amorphous iron oxyhydroxide
      feoh3b_ind,        & ! particulate crystalline iron oxyhydroxide
      fe_ind,            & ! dissolved iron(II)
      fepa_ind,          & ! particulate iron-bound phosphate amorphous
      fepb_ind,          & ! particulate iron-bound phosphate crystalline
      so4_ind,           & ! dissolved sulfate
      h2s_ind,           & ! dissolved hydrogen sulfide
      ch4_ind,           & ! dissolved methane
      hco3_ind,          & ! dissolved bicarbonate
      fes_ind,           & ! particulate iron monosulfide
      fes2_ind,          & ! particulate pyrite
      s_ind,             & ! particulate elemental sulfur
      caco3a_ind,        & ! particulate calcite
      caco3b_ind,        & ! particulate aragonite
      co3_ind,           & ! dissolved carbonate
      camgco3_ind,       & ! partiuclate 15% magnesium calcite
      alk_ind,           & ! alkalinity
      dic_ind              ! dissolved inorganic carbon

   character (benthos_char), allocatable, dimension(:) ::  &
      short_name,      & ! short name of variable
      long_name,       & ! long name of variable
      units              ! units of variable

  end type benthos_indices_type

  type, public :: benthos_element_indices_type
   integer (KIND=benthos_i4) :: &
      carbon_ind,        & ! total carbon
      oxygen_ind,        & ! total oxygen
      nitrogen_ind,      & ! total nitrogen
      phosphorus_ind,    & ! total phosphorus
      sulfur_ind,        & ! total sulfur
      manganese_ind,     & ! total manganese
      iron_ind             ! total iron

   character (benthos_char), allocatable, dimension(:) ::  &
      short_name,      & ! short name of variable
      long_name,       & ! long name of variable
      units              ! units of variable

  end type benthos_element_indices_type

  type, public :: ocean_bgc_indices_type
   integer (benthos_i4) :: &
      po4_ind,          & ! dissolved inorganic phosphate
      no3_ind,          & ! dissolved inorganic nitrate
      sio3_ind,         & ! dissolved inorganic silicate
      nh4_ind,          & ! dissolved ammonia
      fe_ind,           & ! dissolved inorganic iron
      lig_ind,          & ! ligands
      o2_ind,           & ! dissolved oxygen
      dic_ind,          & ! dissolved inorganic carbon
      dic_alt_co2_ind,  & ! dissolved inorganic carbon with alternative CO2
      alk_ind,          & ! alkalinity
      alk_alt_co2_ind,  & ! alkalinity with alternateive CO2
      doc_ind,          & ! dissolved organic carbon
      don_ind,          & ! dissolved organic nitrogen
      dop_ind,          & ! dissolved organic phosphorus
      dopr_ind,         & ! refractory DOP
      donr_ind,         & ! refractory DON
      docr_ind,         & ! refractory DOC
      zooC_ind,         & ! zooplankton carbon
      spChl_ind,        & ! small phyto chlorophyll
      spC_ind,          & ! small phyto carbon
      spFe_ind,         & ! small phyto iron
      spP_ind,          & ! small phyto phosphorus
      spCaCO3_ind,      & ! small phyto calcium carbonate
      diatChl_ind,      & ! diatom chlorophyll
      diatC_ind,        & ! diatom carbon
      diatFe_ind,       & ! diatom iron
      diatP_ind,        & ! diatom phosphorus
      diatSi_ind,       & ! diatom silicon
      !phaeoChl_ind,     & ! Phaeocystis chlorophyll
      !phaeoC_ind,       & ! Phaeocystis carbon
      !phaeoFe_ind,      & ! Phaeocystis iron
      diazChl_ind,      & ! diazotroph chlorophyll
      diazC_ind,        & ! diazotroph carbon
      diazFe_ind,       & ! diazotroph iron
      diazP_ind           ! diazotroph phosphorus

   integer (benthos_i4) :: &
      sp_ind,          & ! autotroph index for small phyto
      diat_ind,        & ! autotroph index for diatoms
      diaz_ind           ! autotroph index for diazotrophs
   !   phaeo_ind          ! autotroph index for Phaeocystis

   !character (benthos_char), allocatable, dimension(:) ::  &
   !   short_name,      & ! short name of variable
   !   long_name,       & ! long name of variable
   !   units              ! units of variable

  end type ocean_bgc_indices_type

  type, public :: benthos_input_type

   real (KIND=benthos_r8), allocatable, dimension(:,:,:) :: &
!      tracerConc, &
      benthosTracerBulk
   real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      PH_PREV_3D, &
      deepStorage
    integer (KIND=benthos_i4), allocatable, dimension(:) :: &
      tracerType, &
      bottom_level
    real (KIND=benthos_r8), allocatable, dimension(:) :: &
      oceanBottomDepth, &
      oceanBottomTemperature, &
      oceanBottomSalinity, &
      oceanBottomSilicate, &   ! and any other bottom values that are needed!!!
      oceanBottomDensity, &
      oceanBottomPressure, &
      oceanBottomPhosphate, &
      oceanBottomAmmonium, &
      oceanBottomNitrate, &
      oceanBottomOxygen, &
      oceanBottomIron, &
      oceanBottomDIC, &
      oceanBottomAlkalinity, &
      oceanBottomManganese, &
      oceanPOCFlux, &
      oceanPOPFlux, &
      oceanPONFlux, &
      oceanParticulateIronFlux, &
      oceanCalciteFlux, &
      oceanBiogenicSilicaFlux, &
      oceanDustFlux, &
      lithogenicFlux, &
      terraPOCFlux

  end type benthos_input_type

  type, public :: benthos_forcing_type   ! not sure I need this?
   real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      deepStoragePrescribedFlux,        &
      benthosOceanSedimentFlux

   real (KIND=benthos_r8), allocatable, dimension(:) :: &
      sediment_FLUX_IN,          &
      iceFraction,               &
      iceKeelDepth,              &
      fastIceFraction

   real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      depositionFlux, &
      riverFlux, &
      gasFlux, &
      seaIceFlux, &
      netFlux

  end type benthos_forcing_type

  type, public :: benthos_output_type
   real (KIND=benthos_r8), allocatable, dimension(:,:,:) :: &
      benthosTendencies, &
      benthosReactionTendencies, &
      benthosTransportTendencies, &
      elementTendencies, &
      elementReactionTendencies, &
      elementTransportTendencies

   real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      elementDeepStorage, &
      PH_PREV_3D, &
      deepStorage
   
    real (KIND=benthos_r8), allocatable, dimension(:) :: &
      oceanBottomDepth, &        !fields that are actually used in the benthos routine
      oceanBottomTemperature, &
      oceanBottomSalinity, &
      oceanBottomDensity, &
      oceanBottomPressure, &
      oceanBottomSilicate, &
      oceanBottomPhosphate, &
      oceanBottomAmmonium, &
      oceanBottomNitrate, &
      oceanBottomOxygen, &
      oceanBottomIron, &
      oceanBottomDIC, &
      oceanBottomAlkalinity, &
      oceanBottomCarbonDioxide, &
      oceanBottomManganese, &
      oceanBottomHydrogenSulfide, &
      oceanBottomMethane, &
      oceanBottomBicarbonate, &
      oceanBottomCarbonate, &
      oceanBottomSulfate, &
      oceanPOCaFlux, & !
      oceanPOCbFlux, & !
      oceanPOCcFlux, & !
      oceanPOPaFlux, & !
      oceanPOPbFlux, & !
      oceanPOPcFlux, & !
      oceanPONaFlux, & !
      oceanPONbFlux, & !
      oceanPONcFlux, & !
      oceanCalciteFlux, & !
      oceanAragoniteFlux, & !
      oceanMgCalciteFlux, & !
      oceanManganeseOxideaFlux, & !
      oceanManganeseOxidebFlux, & !
      oceanIronOxyhydroxideaFlux, & !
      oceanIronOxyhydroxidebFlux, & !
      oceanIronBoundPhosphorusaFlux, & !
      oceanIronBoundPhosphorusbFlux, & !
      oceanPrecipManganeseOxidea, & !
      oceanPrecipManganeseOxideb, &  !
      oceanPrecipIronOxyhydroxidea, & !
      oceanPrecipIronOxyhydroxideb, & !
      oceanPrecipIronBoundPhosphorusa, & !
      oceanPrecipIronBoundPhosphorusb, &
      oceanTotalSedimentMassFlux

  end type benthos_output_type

  type, public :: benthos_flux_diagnostics_type
   real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      deepStorageDiffVelFlux,           &
      deepStorageFlux,                  &
      benthosOceanDiffVelFlux,          &
      benthosOceanFlux,                 &
      elementSedimentation,             &
      elementBurialExchange,            &
      elementOceanExchange

  end type benthos_flux_diagnostics_type

  type, public :: benthos_diagnostics_type
! 2D stuff
    real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      diag_primarySourceTend,   &
      diag_primarySinkTend,     &
      diag_secondarySourceTend, &
      diag_secondarySinkTend,   &
      diag_carbonateSourceTend, &
      diag_carbonateSinkTend,   &
      relErrorTransport, &
      relErrorReactions, &   ! not using
      relErrorElements, &
      relErrorElementsReactions

! more 2D stuff
    real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      diag_netElements,             &
      diag_netBenthosTracersCell,   &
      diag_netBenthosReactionTend,  &
      diag_netElementsReactionTend,  &
      diag_netBenthosTransportTend, &
      diag_netElementsTransportTend, &
      diag_netBenthosTend,          &
      diag_netElementsTend

  end type benthos_diagnostics_type

  !----------------------------------------------------------------------------
  !   ecosystem parameters accessible via namelist input
  !----------------------------------------------------------------------------

   integer(KIND=benthos_i4), public :: &
       benthosTestCase

   real (KIND=benthos_i8), public :: &
       fluxCorrect  ! 1.0 if useFluxCorrection = 1, 0 otherwise

   logical(KIND=benthos_log), public :: &
       useCarbonateSaturation,     & ! Turns on carbonate calculations
       useSecondaryReactions,      & ! Compute secondary reactionset
       useDepthDependentPorosity,  & ! fix Porosity to 0.5 or use  data estimate of porosity distribution with depth
       useNonZeroDiffusivity,      & ! set diffusivity to zero or use biodiffusion and depth dependent molecular diffusion
       useSedimentation,           & ! Turn on or off sedimentation rates
       usePositiveSedimentation,   & ! Set sign of sedimentation rates
       useFastSedimentation,       & ! use fast sedimentation rates (for testing advection)
       useBGCSinkingFlux,          & ! Turn on and off  biochemical sediment fluxes
       useStepInitialProfiles,     & ! use Step function Initial profiles
       useBenthicReactions,        & ! Turn on or off reaction terms
       useFluxCorrection,          & ! correct flux in check_conservation_FCT.m
       useOceanConc,               & ! Turn on dissolved biogeochemical tracers in the ocean
       useDeepSource,              & ! Turn on deep source from below the active benthos layers
       useConstantDiffusivity,     & ! Turn off the depth dependence of diffusivty except no flux
       ! bottom boundary and no flux top boundary for solids
       useOceanBottomState,        & ! Use the ocean ecosys model to force the benthos
       useBenthosMonthlyForcing      ! Use data file to force the benthos
  !---------------------------------------------------------------------
  !     Vertical transport parameters
  !---------------------------------------------------------------------

   real(KIND=benthos_r8), allocatable, dimension(:,:), public :: &
        benthosSoluteDiffusivity, & ! 1/s  normalized solute diffusivity D/h^2
        benthosSoluteDiffusivityI   ! 1/s  normalized solute diffusivity on interface grid
        
   real(KIND=benthos_r8), allocatable, dimension(:), public :: &
        vertBenthosGridThickI,   & ! normalized thickness of interface vertical grid levels
                                   ! used in column summing tracer fields
        vertBenthosGridThick ,   & ! normalized thickness of vertical grid levels for interior points
        benthosMidPointI,        & ! (m) depth of mid-point of interface vertical Grid
        benthosGridPointsI,      & ! location of interface grid values (positive down) n+1
        benthosGridPoints,       & ! location of grid values n
        benthosPorosity,         & ! liquid fraction of benthos grid levels
        benthosPorosityI,        & ! liquid fraction at benthos interface levels
        benthosSolidPorosity,    & ! solid fraction of benthos grid levels
        benthosSolidPorosityI,   & ! solid fraction at benthos interface levels
        benthosDiffusivity,      & ! solute diffusivity at benthos grid levels
        benthosDiffusivityI,     & ! solute diffusivity at benthos interface grid levels
        benthosSolidDiffusivity, & ! solid diffusivity at benthos grid levels
        benthosSolidDiffusivityI, &! solid diffusivity at benthos interface grid levels
        benthosMolecularDiffusivity ! solute molecular diffusivity
        
  real(KIND=benthos_r8), parameter, public :: &
        porosity_inf   = 0.877_benthos_r8, & ! porosity at the bottom boundary
        porosity_o     = 0.943_benthos_r8, & ! porosity at the benthos surface
        porosity_gamma = 0.05301_benthos_r8, & ! (m) Porosity e-folding distance
        fauna_biomass_max = 1000.0_benthos_r8, & ! (g/m2) saturation biomass
        fauna_biomass     = 1000.0_benthos_r8, & ! (g/m2) surface biomass use in biodiffusion
        molecular_diff    = 0.035_benthos_r8, & ! (m2/y) molecular diffusivity (Hensen et al 1998 for nitrate)
        max_bio_diff      = 5.411_benthos_r8, & ! (cm2/y) maximum biodiffusivity
        x_biodiffusion    = 0.0075_benthos_r8, & ! 0.005 but use 0.01 to start (m) biodiffusion e-folding length scale
        T_max_biodiffusion = 286.85_benthos_r8  ! (K) maximum biodiffusion temperature dependence
                                           ! Reed et al 2011

  !---------------------------------------------------------------------
  !     Temperature parameters
  !---------------------------------------------------------------------

!  real(KIND=benthos_r8), parameter :: &
!       Tref = 30.0_BGC_r8, & ! reference temperature (C)
!       Q_10 = 1.5_BGC_r8     ! factor for temperature dependence (non-dim)

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  subroutine benthos_parms_init

! !INPUT/OUTPUT parameterS:

    !---------------------------------------------------------------------------
    !   default namelist settings:  NONE YET
    !---------------------------------------------------------------------------

    ! turn of carbonate saturation subroutines
    useCarbonateSaturation = .true.

    ! turn of carbonate saturation subroutines
    useSecondaryReactions = .true.

    ! fix Porosity to 0.5 or use  data estimate of porosity distribution with depth
    useDepthDependentPorosity = .true.

    ! set diffusivity to zero or use biodiffusion and depth dependent molecular diffusion
    useNonZeroDiffusivity = .true.

    ! Turn on or off sedimentation rates
    useSedimentation = .true.

    ! Set sign of sedimentation rates
    usePositiveSedimentation = .true.

    ! use fast sedimentation rates (for testing advection)
    useFastSedimentation = .false.

    ! Turn on and off  biochemical sediment fluxes
    useBGCSinkingFlux = .true.

    ! use Step function Initial profiles
    useStepInitialProfiles = .false.

    ! Turn on or off reaction terms
    useBenthicReactions = .true.

    ! correct flux in check_conservation_FCT.m
    useFluxCorrection = .false.

    ! Turn on dissolved biogeochemical tracers in the ocean
    useOceanConc = .true.

    ! Turn on deep source from below the active benthos layers
    useDeepSource = .true.

    ! Turn off the depth dependence of diffusivty except no flux bottom boundary and no flux top boundary for solids
    useConstantDiffusivity = .false.

    ! Set to default configuration
    benthosTestCase = 0

  END SUBROUTINE benthos_parms_init

  !*****************************************************************************

END MODULE ocn_benthos_parms
