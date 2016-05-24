module CanopyFluxesMod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CanopyFluxesMod
!
! !DESCRIPTION:
! Calculates the leaf temperature and the leaf fluxes,
! transpiration, photosynthesis and  updates the dew
! accumulation due to evaporation.
!
! !USES:
   use abortutils,   only: endrun
   use clm_varctl,   only: iulog, use_c13, use_cn, use_cndv
   use shr_sys_mod,  only: shr_sys_flush
!
! !PUBLIC TYPES:
   implicit none
   save
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: CanopyFluxes !Calculates the leaf temperature and leaf fluxes
!
! !PRIVATE MEMBER FUNCTIONS:
   private :: Stomata     !Leaf stomatal resistance and leaf photosynthesis
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 4/25/05, Peter Thornton: replaced old Stomata subroutine with what
!   used to be called StomataCN, as part of migration to new sun/shade
!   algorithms. 
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CanopyFluxes 
!
! !INTERFACE:
  subroutine CanopyFluxes(lbg, ubg, lbc, ubc, lbp, ubp, &
                          num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! 1. Calculates the leaf temperature:
! 2. Calculates the leaf fluxes, transpiration, photosynthesis and
!    updates the dew accumulation due to evaporation.
!
! Method:
! Use the Newton-Raphson iteration to solve for the foliage
! temperature that balances the surface energy budget:
!
! f(t_veg) = Net radiation - Sensible - Latent = 0
! f(t_veg) + d(f)/d(t_veg) * dt_veg = 0     (*)
!
! Note:
! (1) In solving for t_veg, t_grnd is given from the previous timestep.
! (2) The partial derivatives of aerodynamical resistances, which cannot
!     be determined analytically, are ignored for d(H)/dT and d(LE)/dT
! (3) The weighted stomatal resistance of sunlit and shaded foliage is used
! (4) Canopy air temperature and humidity are derived from => Hc + Hg = Ha
!                                                          => Ec + Eg = Ea
! (5) Energy loss is due to: numerical truncation of energy budget equation
!     (*); and "ecidif" (see the code) which is dropped into the sensible
!     heat
! (6) The convergence criteria: the difference, del = t_veg(n+1)-t_veg(n)
!     and del2 = t_veg(n)-t_veg(n-1) less than 0.01 K, and the difference
!     of water flux from the leaf between the iteration step (n+1) and (n)
!     less than 0.1 W/m2; or the iterative steps over 40.
!
! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use clmtype
    use clm_atmlnd         , only : clm_a2l
    use clm_time_manager   , only : get_step_size, get_prev_date
    use clm_varpar         , only : nlevgrnd, nlevsno
    use clm_varcon         , only : sb, cpair, hvap, vkc, grav, denice, &
                                    denh2o, tfrz, csoilc, tlsai_crit, alpha_aero, &
                                    isecspday, degpsec
    use shr_const_mod      , only : SHR_CONST_TKFRZ
    use pftvarcon          , only : nirrig
    use QSatMod            , only : QSat
    use FrictionVelocityMod, only : FrictionVelocity, MoninObukIni
    use spmdMod            , only : masterproc
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbg, ubg                    ! gridcell bounds
    integer, intent(in) :: lbc, ubc                    ! column bounds
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: num_nolakep                 ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)   ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 12/19/01, Peter Thornton
! Changed tg to t_grnd for consistency with other routines
! 1/29/02, Peter Thornton
! Migrate to new data structures, new calling protocol. For now co2 and
! o2 partial pressures are hardwired, but they should be coming in from
! forc_pco2 and forc_po2. Keeping the same hardwired values as in CLM2 to
! assure bit-for-bit results in the first comparisons.
! 27 February 2008: Keith Oleson; Sparse/dense aerodynamic parameters from
! X. Zeng
! 6 March 2009: Peter Thornton; Daylength control on Vcmax, from Bill Bauerle
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in variables
!
   integer , pointer :: frac_veg_nosno(:) ! frac of veg not covered by snow (0 OR 1 now) [-]
   integer , pointer :: ivt(:)         ! pft vegetation type
   integer , pointer :: pcolumn(:)     ! pft's column index
   integer , pointer :: plandunit(:)   ! pft's landunit index
   integer , pointer :: pgridcell(:)   ! pft's gridcell index
   real(r8), pointer :: forc_th(:)     ! atmospheric potential temperature (Kelvin)
   real(r8), pointer :: t_grnd(:)      ! ground surface temperature [K]
   real(r8), pointer :: thm(:)         ! intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
   real(r8), pointer :: qg(:)          ! specific humidity at ground surface [kg/kg]
   real(r8), pointer :: thv(:)         ! virtual potential temperature (kelvin)
   real(r8), pointer :: z0mv(:)        ! roughness length over vegetation, momentum [m]
   real(r8), pointer :: z0hv(:)        ! roughness length over vegetation, sensible heat [m]
   real(r8), pointer :: z0qv(:)        ! roughness length over vegetation, latent heat [m]
   real(r8), pointer :: z0mg(:)        ! roughness length of ground, momentum [m]
   real(r8), pointer :: dqgdT(:)       ! temperature derivative of "qg"
   real(r8), pointer :: htvp(:)        ! latent heat of evaporation (/sublimation) [J/kg]
   real(r8), pointer :: emv(:)         ! ground emissivity
   real(r8), pointer :: emg(:)         ! vegetation emissivity
   real(r8), pointer :: forc_pbot(:)   ! atmospheric pressure (Pa)
   real(r8), pointer :: forc_pco2(:)   ! partial pressure co2 (Pa)

   ! 4/14/05: PET
   ! Adding isotope code
   real(r8), pointer :: forc_pc13o2(:) ! partial pressure c13o2 (Pa)
   
   real(r8), pointer :: forc_po2(:)    ! partial pressure o2 (Pa)
   real(r8), pointer :: forc_q(:)      ! atmospheric specific humidity (kg/kg)
   real(r8), pointer :: forc_u(:)      ! atmospheric wind speed in east direction (m/s)
   real(r8), pointer :: forc_v(:)      ! atmospheric wind speed in north direction (m/s)
   real(r8), pointer :: forc_hgt_u_pft(:) !observational height of wind at pft level [m]
   real(r8), pointer :: forc_rho(:)    ! density (kg/m**3)
   real(r8), pointer :: forc_lwrad(:)  ! downward infrared (longwave) radiation (W/m**2)
   real(r8), pointer :: displa(:)      ! displacement height (m)
   real(r8), pointer :: elai(:)        ! one-sided leaf area index with burying by snow
   real(r8), pointer :: esai(:)        ! one-sided stem area index with burying by snow
   real(r8), pointer :: fdry(:)        ! fraction of foliage that is green and dry [-]
   real(r8), pointer :: fwet(:)        ! fraction of canopy that is wet (0 to 1)
   real(r8), pointer :: laisun(:)      ! sunlit leaf area
   real(r8), pointer :: laisha(:)      ! shaded leaf area
   real(r8), pointer :: sabv(:)        ! solar radiation absorbed by vegetation (W/m**2)
   real(r8), pointer :: watsat(:,:)    ! volumetric soil water at saturation (porosity)
   real(r8), pointer :: watdry(:,:)    ! btran parameter for btran=0
   real(r8), pointer :: watopt(:,:)    ! btran parameter for btran = 1
   real(r8), pointer :: h2osoi_ice(:,:)! ice lens (kg/m2)
   real(r8), pointer :: h2osoi_liq(:,:)! liquid water (kg/m2)
   real(r8), pointer :: dz(:,:)        ! layer depth (m)
   real(r8), pointer :: t_soisno(:,:)  ! soil temperature (Kelvin)
   real(r8), pointer :: sucsat(:,:)    ! minimum soil suction (mm)
   real(r8), pointer :: bsw(:,:)       ! Clapp and Hornberger "b"
   real(r8), pointer :: rootfr(:,:)    ! fraction of roots in each soil layer
   real(r8), pointer :: dleaf(:)       ! characteristic leaf dimension (m)
   real(r8), pointer :: smpso(:)       ! soil water potential at full stomatal opening (mm)
   real(r8), pointer :: smpsc(:)       ! soil water potential at full stomatal closure (mm)
   real(r8), pointer :: frac_sno(:)    ! fraction of ground covered by snow (0 to 1)
   real(r8), pointer :: htop(:)        ! canopy top(m)
   real(r8), pointer :: snowdp(:)      ! snow height (m)
   real(r8), pointer :: soilbeta(:)    ! soil wetness relative to field capacity
   real(r8), pointer :: lat(:)         ! latitude (radians)
   real(r8), pointer :: decl(:)        ! declination angle (radians)
   real(r8), pointer :: max_dayl(:)    !maximum daylength for this column (s)
   real(r8), pointer :: londeg(:)      ! longitude (degrees) (for calculation of local time)
   
!
! local pointers to implicit inout arguments
!
   real(r8), pointer :: cgrnds(:)      ! deriv. of soil sensible heat flux wrt soil temp [w/m2/k]
   real(r8), pointer :: cgrndl(:)      ! deriv. of soil latent heat flux wrt soil temp [w/m**2/k]
   real(r8), pointer :: t_veg(:)       ! vegetation temperature (Kelvin)
   real(r8), pointer :: t_ref2m(:)     ! 2 m height surface air temperature (Kelvin)
   real(r8), pointer :: q_ref2m(:)     ! 2 m height surface specific humidity (kg/kg)
   real(r8), pointer :: t_ref2m_r(:)   ! Rural 2 m height surface air temperature (Kelvin)
   real(r8), pointer :: rh_ref2m(:)    ! 2 m height surface relative humidity (%)
   real(r8), pointer :: rh_ref2m_r(:)  ! Rural 2 m height surface relative humidity (%)
   real(r8), pointer :: h2ocan(:)      ! canopy water (mm H2O)
   real(r8), pointer :: cisun(:)       !sunlit intracellular CO2 (Pa)
   real(r8), pointer :: cisha(:)       !shaded intracellular CO2 (Pa)
!
! local pointers to implicit out arguments
!
   real(r8), pointer :: rb1(:)             ! boundary layer resistance (s/m)
   real(r8), pointer :: cgrnd(:)           ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
   real(r8), pointer :: dlrad(:)           ! downward longwave radiation below the canopy [W/m2]
   real(r8), pointer :: ulrad(:)           ! upward longwave radiation above the canopy [W/m2]
   real(r8), pointer :: ram1(:)            ! aerodynamical resistance (s/m)
   real(r8), pointer :: btran(:)           ! transpiration wetness factor (0 to 1)
   real(r8), pointer :: rssun(:)           ! sunlit stomatal resistance (s/m)
   real(r8), pointer :: rssha(:)           ! shaded stomatal resistance (s/m)
   real(r8), pointer :: psnsun(:)          ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: psnsha(:)          ! shaded leaf photosynthesis (umol CO2 /m**2/ s)

   ! 4/14/05: PET
   ! Adding isotope code
   real(r8), pointer :: c13_psnsun(:)      ! sunlit leaf photosynthesis (umol 13CO2 /m**2/ s)
   real(r8), pointer :: c13_psnsha(:)      ! shaded leaf photosynthesis (umol 13CO2 /m**2/ s)
   ! 4/21/05: PET
   ! Adding isotope code
   real(r8), pointer :: rc13_canair(:)     !C13O2/C12O2 in canopy air
   real(r8), pointer :: rc13_psnsun(:)     !C13O2/C12O2 in sunlit canopy psn flux
   real(r8), pointer :: rc13_psnsha(:)     !C13O2/C12O2 in shaded canopy psn flux
   real(r8), pointer :: alphapsnsun(:)     !fractionation factor in sunlit canopy psn flux
   real(r8), pointer :: alphapsnsha(:)     !fractionation factor in shaded canopy psn flux
   
   real(r8), pointer :: qflx_tran_veg(:)   ! vegetation transpiration (mm H2O/s) (+ = to atm)
   real(r8), pointer :: dt_veg(:)          ! change in t_veg, last iteration (Kelvin)
   real(r8), pointer :: qflx_evap_veg(:)   ! vegetation evaporation (mm H2O/s) (+ = to atm)
   real(r8), pointer :: eflx_sh_veg(:)     ! sensible heat flux from leaves (W/m**2) [+ to atm]
   real(r8), pointer :: taux(:)            ! wind (shear) stress: e-w (kg/m/s**2)
   real(r8), pointer :: tauy(:)            ! wind (shear) stress: n-s (kg/m/s**2)
   real(r8), pointer :: eflx_sh_grnd(:)    ! sensible heat flux from ground (W/m**2) [+ to atm]
   real(r8), pointer :: qflx_evap_soi(:)   ! soil evaporation (mm H2O/s) (+ = to atm)
   real(r8), pointer :: fpsn(:)            ! photosynthesis (umol CO2 /m**2 /s)
   real(r8), pointer :: rootr(:,:)         ! effective fraction of roots in each soil layer
   real(r8), pointer :: rresis(:,:)        ! root resistance by layer (0-1)  (nlevgrnd)	
   real(r8), pointer :: irrig_rate(:)      ! current irrigation rate [mm/s]
   integer, pointer  :: n_irrig_steps_left(:) ! # of time steps for which we still need to irrigate today
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
   real(r8), parameter :: btran0 = 0.0_r8  ! initial value
   real(r8), parameter :: zii = 1000.0_r8  ! convective boundary layer height [m]
   real(r8), parameter :: beta = 1.0_r8    ! coefficient of conective velocity [-]
   real(r8), parameter :: delmax = 1.0_r8  ! maxchange in  leaf temperature [K]
   real(r8), parameter :: dlemin = 0.1_r8  ! max limit for energy flux convergence [w/m2]
   real(r8), parameter :: dtmin = 0.01_r8  ! max limit for temperature convergence [K]
   integer , parameter :: itmax = 40       ! maximum number of iteration [-]
   integer , parameter :: itmin = 2        ! minimum number of iteration [-]
   real(r8), parameter :: irrig_min_lai = 0.0_r8           ! Minimum LAI for irrigation
   real(r8), parameter :: irrig_btran_thresh = 0.999999_r8 ! Irrigate when btran falls below 0.999999 rather than 1 to allow for round-off error
   integer , parameter :: irrig_start_time = isecspday/4   ! Time of day to check whether we need irrigation, seconds (0 = midnight). 
                                                           ! We start applying the irrigation in the time step FOLLOWING this time, 
                                                           ! since we won't begin irrigating until the next call to Hydrology1
   integer , parameter :: irrig_length = isecspday/6       ! Desired amount of time to irrigate per day (sec). Actual time may 
                                                           ! differ if this is not a multiple of dtime. Irrigation won't work properly 
                                                           ! if dtime > secsperday
   real(r8), parameter :: irrig_factor = 0.7_r8            ! Determines target soil moisture level for irrigation. If h2osoi_liq_so 
                                                           ! is the soil moisture level at which stomata are fully open and 
                                                           ! h2osoi_liq_sat is the soil moisture level at saturation (eff_porosity), 
                                                           ! then the target soil moisture level is 
                                                           !     (h2osoi_liq_so + irrig_factor*(h2osoi_liq_sat - h2osoi_liq_so)). 
                                                           ! A value of 0 means that the target soil moisture level is h2osoi_liq_so. 
                                                           ! A value of 1 means that the target soil moisture level is h2osoi_liq_sat

   !added by K.Sakaguchi for litter resistance
   real(r8), parameter :: lai_dl = 0.5_r8  ! placeholder for (dry) plant litter area index (m2/m2)
   real(r8), parameter :: z_dl = 0.05_r8   ! placeholder for (dry) litter layer thickness (m)
   !added by K.Sakaguchi for stability formulation
   real(r8), parameter :: ria  = 0.5_r8    ! free parameter for stable formulation (currently = 0.5, "gamma" in Sakaguchi&Zeng,2008)
   real(r8) :: dtime                 ! land model time step (sec)
   real(r8) :: zldis(lbp:ubp)        ! reference height "minus" zero displacement height [m]
   real(r8) :: zeta                  ! dimensionless height used in Monin-Obukhov theory
   real(r8) :: wc                    ! convective velocity [m/s]
   real(r8) :: dth(lbp:ubp)          ! diff of virtual temp. between ref. height and surface
   real(r8) :: dthv(lbp:ubp)         ! diff of vir. poten. temp. between ref. height and surface
   real(r8) :: dqh(lbp:ubp)          ! diff of humidity between ref. height and surface
   real(r8) :: obu(lbp:ubp)          ! Monin-Obukhov length (m)
   real(r8) :: um(lbp:ubp)           ! wind speed including the stablity effect [m/s]
   real(r8) :: ur(lbp:ubp)           ! wind speed at reference height [m/s]
   real(r8) :: uaf(lbp:ubp)          ! velocity of air within foliage [m/s]
   real(r8) :: temp1(lbp:ubp)        ! relation for potential temperature profile
   real(r8) :: temp12m(lbp:ubp)      ! relation for potential temperature profile applied at 2-m
   real(r8) :: temp2(lbp:ubp)        ! relation for specific humidity profile
   real(r8) :: temp22m(lbp:ubp)      ! relation for specific humidity profile applied at 2-m
   real(r8) :: ustar(lbp:ubp)        ! friction velocity [m/s]
   real(r8) :: tstar                 ! temperature scaling parameter
   real(r8) :: qstar                 ! moisture scaling parameter
   real(r8) :: thvstar               ! virtual potential temperature scaling parameter
   real(r8) :: taf(lbp:ubp)          ! air temperature within canopy space [K]
   real(r8) :: qaf(lbp:ubp)          ! humidity of canopy air [kg/kg]
   real(r8) :: rpp                   ! fraction of potential evaporation from leaf [-]
   real(r8) :: rppdry                ! fraction of potential evaporation through transp [-]
   real(r8) :: cf                    ! heat transfer coefficient from leaves [-]
   real(r8) :: rb(lbp:ubp)           ! leaf boundary layer resistance [s/m]
   real(r8) :: rah(lbp:ubp,2)        ! thermal resistance [s/m]
   real(r8) :: raw(lbp:ubp,2)        ! moisture resistance [s/m]
   real(r8) :: wta                   ! heat conductance for air [m/s]
   real(r8) :: wtg(lbp:ubp)          ! heat conductance for ground [m/s]
   real(r8) :: wtl                   ! heat conductance for leaf [m/s]
   real(r8) :: wta0(lbp:ubp)         ! normalized heat conductance for air [-]
   real(r8) :: wtl0(lbp:ubp)         ! normalized heat conductance for leaf [-]
   real(r8) :: wtg0                  ! normalized heat conductance for ground [-]
   real(r8) :: wtal(lbp:ubp)         ! normalized heat conductance for air and leaf [-]
   real(r8) :: wtga                  ! normalized heat cond. for air and ground  [-]
   real(r8) :: wtaq                  ! latent heat conductance for air [m/s]
   real(r8) :: wtlq                  ! latent heat conductance for leaf [m/s]
   real(r8) :: wtgq(lbp:ubp)         ! latent heat conductance for ground [m/s]
   real(r8) :: wtaq0(lbp:ubp)        ! normalized latent heat conductance for air [-]
   real(r8) :: wtlq0(lbp:ubp)        ! normalized latent heat conductance for leaf [-]
   real(r8) :: wtgq0                 ! normalized heat conductance for ground [-]
   real(r8) :: wtalq(lbp:ubp)        ! normalized latent heat cond. for air and leaf [-]
   real(r8) :: wtgaq                 ! normalized latent heat cond. for air and ground [-]
   real(r8) :: el(lbp:ubp)           ! vapor pressure on leaf surface [pa]
   real(r8) :: deldT                 ! derivative of "el" on "t_veg" [pa/K]
   real(r8) :: qsatl(lbp:ubp)        ! leaf specific humidity [kg/kg]
   real(r8) :: qsatldT(lbp:ubp)      ! derivative of "qsatl" on "t_veg"
   real(r8) :: e_ref2m               ! 2 m height surface saturated vapor pressure [Pa]
   real(r8) :: de2mdT                ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
   real(r8) :: qsat_ref2m            ! 2 m height surface saturated specific humidity [kg/kg]
   real(r8) :: dqsat2mdT             ! derivative of 2 m height surface saturated specific humidity on t_ref2m
   real(r8) :: air(lbp:ubp),bir(lbp:ubp),cir(lbp:ubp)  ! atmos. radiation temporay set
   real(r8) :: dc1,dc2               ! derivative of energy flux [W/m2/K]
   real(r8) :: delt                  ! temporary
   real(r8) :: delq(lbp:ubp)         ! temporary
   real(r8) :: del(lbp:ubp)          ! absolute change in leaf temp in current iteration [K]
   real(r8) :: del2(lbp:ubp)         ! change in leaf temperature in previous iteration [K]
   real(r8) :: dele(lbp:ubp)         ! change in latent heat flux from leaf [K]
   real(r8) :: dels                  ! change in leaf temperature in current iteration [K]
   real(r8) :: det(lbp:ubp)          ! maximum leaf temp. change in two consecutive iter [K]
   real(r8) :: efeb(lbp:ubp)         ! latent heat flux from leaf (previous iter) [mm/s]
   real(r8) :: efeold                ! latent heat flux from leaf (previous iter) [mm/s]
   real(r8) :: efpot                 ! potential latent energy flux [kg/m2/s]
   real(r8) :: efe(lbp:ubp)          ! water flux from leaf [mm/s]
   real(r8) :: efsh                  ! sensible heat from leaf [mm/s]
   real(r8) :: obuold(lbp:ubp)       ! monin-obukhov length from previous iteration
   real(r8) :: tlbef(lbp:ubp)        ! leaf temperature from previous iteration [K]
   real(r8) :: ecidif                ! excess energies [W/m2]
   real(r8) :: err(lbp:ubp)          ! balance error
   real(r8) :: erre                  ! balance error
   real(r8) :: co2(lbp:ubp)          ! atmospheric co2 partial pressure (pa)

   ! 4/14/05: PET
   ! Adding isotope code
   real(r8) :: c13o2(lbp:ubp)        ! atmospheric c13o2 partial pressure (pa)
   
   real(r8) :: o2(lbp:ubp)           ! atmospheric o2 partial pressure (pa)
   real(r8) :: svpts(lbp:ubp)        ! saturation vapor pressure at t_veg (pa)
   real(r8) :: eah(lbp:ubp)          ! canopy air vapor pressure (pa)
   real(r8) :: s_node                ! vol_liq/eff_porosity
   real(r8) :: smp_node              ! matrix potential
   real(r8) :: vol_ice               ! partial volume of ice lens in layer
   real(r8) :: eff_porosity          ! effective porosity in layer
   real(r8) :: vol_liq               ! partial volume of liquid water in layer
   integer  :: itlef                 ! counter for leaf temperature iteration [-]
   integer  :: nmozsgn(lbp:ubp)      ! number of times stability changes sign
   real(r8) :: w                     ! exp(-LSAI)
   real(r8) :: csoilcn               ! interpolated csoilc for less than dense canopies
   real(r8) :: fm(lbp:ubp)           ! needed for BGC only to diagnose 10m wind speed
   real(r8) :: wtshi                 ! sensible heat resistance for air, grnd and leaf [-]
   real(r8) :: wtsqi                 ! latent heat resistance for air, grnd and leaf [-]
   integer  :: j                     ! soil/snow level index
   integer  :: p                     ! pft index
   integer  :: c                     ! column index
   integer  :: l                     ! landunit index
   integer  :: g                     ! gridcell index
   integer  :: fp                    ! lake filter pft index
   integer  :: fn                    ! number of values in pft filter
   integer  :: fnorig                ! number of values in pft filter copy
   integer  :: fnold                 ! temporary copy of pft count
   integer  :: f                     ! filter index
   integer  :: filterp(ubp-lbp+1)    ! temporary filter
   integer  :: fporig(ubp-lbp+1)     ! temporary filter
   real(r8) :: displa_loc(lbp:ubp)   ! temporary copy
   real(r8) :: z0mv_loc(lbp:ubp)     ! temporary copy
   real(r8) :: z0hv_loc(lbp:ubp)     ! temporary copy
   real(r8) :: z0qv_loc(lbp:ubp)     ! temporary copy
   logical  :: found                 ! error flag for canopy above forcing hgt
   integer  :: index                 ! pft index for error
   real(r8) :: egvf                  ! effective green vegetation fraction
   real(r8) :: lt                    ! elai+esai
   real(r8) :: ri                    ! stability parameter for under canopy air (unitless)
   real(r8) :: csoilb                ! turbulent transfer coefficient over bare soil (unitless)
   real(r8) :: ricsoilc              ! modified transfer coefficient under dense canopy (unitless)
   real(r8) :: snowdp_c              ! critical snow depth to cover plant litter (m)
   real(r8) :: rdl                   ! dry litter layer resistance for water vapor  (s/m)
   real(r8) :: elai_dl               ! exposed (dry) plant litter area index
   real(r8) :: fsno_dl               ! effective snow cover over plant litter
   real(r8) :: dayl                  ! daylength (s)
   real(r8) :: temp                  ! temporary, for daylength calculation
   real(r8) :: dayl_factor(lbp:ubp)  ! scalar (0-1) for daylength effect on Vcmax
   integer  :: yr                       ! year at start of time step
   integer  :: mon                      ! month at start of time step
   integer  :: day                      ! day at start of time step
   integer  :: time                     ! time at start of time step (seconds after 0Z)
   integer  :: local_time               ! local time at start of time step (seconds after solar midnight)
   integer  :: irrig_nsteps_per_day     ! number of time steps per day in which we irrigate
   logical  :: check_for_irrig(lbp:ubp) ! where do we need to check soil moisture to see if we need to irrigate?
   logical  :: frozen_soil(lbc:ubc)     ! set to true if we have encountered a frozen soil layer
   real(r8) :: vol_liq_so               ! partial volume of liquid water in layer for which smp_node = smpso
   real(r8) :: h2osoi_liq_so            ! liquid water corresponding to vol_liq_so for this layer [kg/m2]
   real(r8) :: h2osoi_liq_sat           ! liquid water corresponding to eff_porosity for this layer [kg/m2]
   real(r8) :: deficit                  ! difference between desired soil moisture level for this layer and current soil moisture level [kg/m2]
!------------------------------------------------------------------------------

   ! Assign local pointers to derived type members (gridcell-level)

   forc_lwrad     => clm_a2l%forc_lwrad
   forc_pco2      => clm_a2l%forc_pco2

   forc_pc13o2    => clm_a2l%forc_pc13o2

   forc_po2       => clm_a2l%forc_po2
   forc_q         => clm_a2l%forc_q
   forc_pbot      => clm_a2l%forc_pbot
   forc_u         => clm_a2l%forc_u
   forc_v         => clm_a2l%forc_v
   forc_th        => clm_a2l%forc_th
   forc_rho       => clm_a2l%forc_rho
   lat            => grc%lat
   londeg         => grc%londeg

   ! Assign local pointers to derived type members (column-level)

   t_soisno       => ces%t_soisno
   watsat         => cps%watsat
   watdry         => cps%watdry 
   watopt         => cps%watopt 
   h2osoi_ice     => cws%h2osoi_ice
   dz             => cps%dz
   h2osoi_liq     => cws%h2osoi_liq
   sucsat         => cps%sucsat
   bsw            => cps%bsw
   emg            => cps%emg
   t_grnd         => ces%t_grnd
   qg             => cws%qg
   thv            => ces%thv
   dqgdT          => cws%dqgdT
   htvp           => cps%htvp
   z0mg           => cps%z0mg
   frac_sno       => cps%frac_sno
   snowdp         => cps%snowdp
   soilbeta       => cws%soilbeta
   decl           => cps%decl
   max_dayl       => cps%max_dayl

   ! Assign local pointers to derived type members (pft-level)

   rb1            => pps%rb1
   ivt            => pft%itype
   pcolumn        => pft%column
   plandunit      => pft%landunit
   pgridcell      => pft%gridcell
   frac_veg_nosno => pps%frac_veg_nosno
   btran          => pps%btran
   rootfr         => pps%rootfr
   rootr          => pps%rootr
   rresis         => pps%rresis
   emv            => pps%emv
   t_veg          => pes%t_veg
   displa         => pps%displa
   z0mv           => pps%z0mv
   z0hv           => pps%z0hv
   z0qv           => pps%z0qv
   ram1           => pps%ram1
   htop           => pps%htop
   rssun          => pps%rssun
   rssha          => pps%rssha
   cisun          => pps%cisun
   cisha          => pps%cisha
   psnsun         => pcf%psnsun
   psnsha         => pcf%psnsha

   ! 4/14/05: PET
   ! Adding isotope code
   c13_psnsun     => pc13f%psnsun
   c13_psnsha     => pc13f%psnsha
   ! 4/21/05: PET
   ! Adding isotope code
   rc13_canair    => pepv%rc13_canair
   rc13_psnsun    => pepv%rc13_psnsun
   rc13_psnsha    => pepv%rc13_psnsha
   alphapsnsun    => pps%alphapsnsun
   alphapsnsha    => pps%alphapsnsha
   
   elai           => pps%elai
   esai           => pps%esai
   fdry           => pps%fdry
   laisun         => pps%laisun
   laisha         => pps%laisha
   qflx_tran_veg  => pwf%qflx_tran_veg
   fwet           => pps%fwet
   h2ocan         => pws%h2ocan
   dt_veg         => pps%dt_veg
   sabv           => pef%sabv
   qflx_evap_veg  => pwf%qflx_evap_veg
   eflx_sh_veg    => pef%eflx_sh_veg
   taux           => pmf%taux
   tauy           => pmf%tauy
   eflx_sh_grnd   => pef%eflx_sh_grnd
   qflx_evap_soi  => pwf%qflx_evap_soi
   t_ref2m        => pes%t_ref2m
   q_ref2m        => pes%q_ref2m
   t_ref2m_r      => pes%t_ref2m_r
   rh_ref2m_r     => pes%rh_ref2m_r
   rh_ref2m       => pes%rh_ref2m
   dlrad          => pef%dlrad
   ulrad          => pef%ulrad
   cgrnds         => pef%cgrnds
   cgrndl         => pef%cgrndl
   cgrnd          => pef%cgrnd
   fpsn           => pcf%fpsn
   forc_hgt_u_pft => pps%forc_hgt_u_pft
   thm            => pes%thm
   irrig_rate         => cps%irrig_rate
   n_irrig_steps_left => cps%n_irrig_steps_left
      
   ! Assign local pointers to derived type members (ecophysiological)

   dleaf          => pftcon%dleaf
   smpso          => pftcon%smpso
   smpsc          => pftcon%smpsc

   ! Determine step size

   dtime = get_step_size()
   irrig_nsteps_per_day = ((irrig_length + (dtime - 1))/dtime)  ! round up

   ! Filter pfts where frac_veg_nosno is non-zero

   fn = 0
   do fp = 1,num_nolakep
      p = filter_nolakep(fp)
      if (frac_veg_nosno(p) /= 0) then
         fn = fn + 1
         filterp(fn) = p
      end if
   end do

   ! Initialize

   do f = 1, fn
      p = filterp(f)
      del(p)    = 0._r8  ! change in leaf temperature from previous iteration
      efeb(p)   = 0._r8  ! latent head flux from leaf for previous iteration
      wtlq0(p)  = 0._r8
      wtalq(p)  = 0._r8
      wtgq(p)   = 0._r8
      wtaq0(p)  = 0._r8
      obuold(p) = 0._r8
      btran(p)  = btran0
   end do
   
   ! calculate daylength control for Vcmax
   do f = 1, fn
      p=filterp(f)
      c=pcolumn(p)
      g=pgridcell(p)
      ! calculate daylength
      temp = -(sin(lat(g))*sin(decl(c)))/(cos(lat(g)) * cos(decl(c)))
      temp = min(1._r8,max(-1._r8,temp))
      dayl = 2.0_r8 * 13750.9871_r8 * acos(temp)
      ! calculate dayl_factor as the ratio of (current:max dayl)^2
      ! set a minimum of 0.01 (1%) for the dayl_factor
      dayl_factor(p)=min(1._r8,max(0.01_r8,(dayl*dayl)/(max_dayl(c)*max_dayl(c))))
   end do

   rb1(lbp:ubp) = 0._r8

   ! Effective porosity of soil, partial volume of ice and liquid (needed for btran)
   ! and root resistance factors

   do j = 1,nlevgrnd
      do f = 1, fn
         p = filterp(f)
         c = pcolumn(p)
         l = plandunit(p)

         ! Root resistance factors

         vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
         eff_porosity = watsat(c,j)-vol_ice
         vol_liq = min(eff_porosity, h2osoi_liq(c,j)/(dz(c,j)*denh2o))
         if (vol_liq .le. 0._r8 .or. t_soisno(c,j) .le. tfrz-2._r8) then
            rootr(p,j) = 0._r8
         else
            s_node = max(vol_liq/eff_porosity,0.01_r8)
            smp_node = max(smpsc(ivt(p)), -sucsat(c,j)*s_node**(-bsw(c,j)))

            rresis(p,j) = min( (eff_porosity/watsat(c,j))* &
                          (smp_node - smpsc(ivt(p))) / (smpso(ivt(p)) - smpsc(ivt(p))), 1._r8)
            rootr(p,j) = rootfr(p,j)*rresis(p,j)
            btran(p) = btran(p) + rootr(p,j)
         endif 
      end do
   end do

   ! Normalize root resistances to get layer contribution to ET

   do j = 1,nlevgrnd
      do f = 1, fn
         p = filterp(f)
         if (btran(p) .gt. btran0) then
           rootr(p,j) = rootr(p,j)/btran(p)
         else
           rootr(p,j) = 0._r8
         end if
      end do
   end do

   ! Determine if irrigation is needed (over irrigated soil columns)

   ! First, determine in what grid cells we need to bother 'measuring' soil water, to see if we need irrigation
   ! Also set n_irrig_steps_left for these grid cells
   ! n_irrig_steps_left(p) > 0 is ok even if irrig_rate(p) ends up = 0
   ! in this case, we'll irrigate by 0 for the given number of time steps
   call get_prev_date(yr, mon, day, time)  ! get time as of beginning of time step
   do f = 1, fn
      p = filterp(f)
      c = pcolumn(p)
      g = pgridcell(p)
      if (ivt(p) == nirrig .and. elai(p) > irrig_min_lai .and. btran(p) < irrig_btran_thresh) then
         ! see if it's the right time of day to start irrigating:
         local_time = modulo(time + nint(londeg(g)/degpsec), isecspday)
         if (modulo(local_time - irrig_start_time, isecspday) < dtime) then
            ! it's time to start irrigating
            check_for_irrig(p)    = .true.
            n_irrig_steps_left(c) = irrig_nsteps_per_day
            irrig_rate(c)         = 0._r8  ! reset; we'll add to this later
         else
            check_for_irrig(p)    = .false.
         end if
      else  ! non-irrig pft or elai<=irrig_min_lai or btran>irrig_btran_thresh
         check_for_irrig(p)       = .false.
      end if

   end do


   ! Now 'measure' soil water for the grid cells identified above and see if the soil is dry enough to warrant irrigation
   frozen_soil(:) = .false.
   do j = 1,nlevgrnd
      do f = 1, fn
         p = filterp(f)
         c = pcolumn(p)
         if (check_for_irrig(p) .and. .not. frozen_soil(c)) then
            ! if level L was frozen, then we don't look at any levels below L
            if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
               frozen_soil(c) = .true.
            else if (rootfr(p,j) > 0._r8) then
               ! determine soil water deficit in this layer:

               ! Calculate vol_liq_so - i.e., vol_liq at which smp_node = smpso - by inverting the above equations 
               ! for the root resistance factors
               vol_ice      = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))  ! this duplicates the above equation for vol_ice
               eff_porosity = watsat(c,j)-vol_ice  ! this duplicates the above equation for eff_porosity
               vol_liq_so   = eff_porosity * (-smpso(ivt(p))/sucsat(c,j))**(-1/bsw(c,j))

               ! Translate vol_liq_so and eff_porosity into h2osoi_liq_so and h2osoi_liq_sat and calculate deficit
               h2osoi_liq_so  = vol_liq_so * denh2o * dz(c,j)
               h2osoi_liq_sat = eff_porosity * denh2o * dz(c,j)
               deficit        = max((h2osoi_liq_so + irrig_factor*(h2osoi_liq_sat - h2osoi_liq_so)) - h2osoi_liq(c,j), 0._r8)

               ! Add deficit to irrig_rate, converting units from mm to mm/sec
               irrig_rate(c)  = irrig_rate(c) + deficit/(dtime*irrig_nsteps_per_day)

            end if  ! else if (rootfr(p,j) .gt. 0)
         end if     ! if (check_for_irrig(p) .and. .not. frozen_soil(c))
      end do        ! do f
   end do           ! do j

 ! Modify aerodynamic parameters for sparse/dense canopy (X. Zeng)
   do f = 1, fn
      p = filterp(f)
      c = pcolumn(p)

      lt = min(elai(p)+esai(p), tlsai_crit)
      egvf =(1._r8 - alpha_aero * exp(-lt)) / (1._r8 - alpha_aero * exp(-tlsai_crit))
      displa(p) = egvf * displa(p)
      z0mv(p)   = exp(egvf * log(z0mv(p)) + (1._r8 - egvf) * log(z0mg(c)))
      z0hv(p)   = z0mv(p)
      z0qv(p)   = z0mv(p)

  end do

   found = .false.
   do f = 1, fn
      p = filterp(f)
      c = pcolumn(p)
      g = pgridcell(p)

      ! Net absorbed longwave radiation by canopy and ground
      ! =air+bir*t_veg**4+cir*t_grnd(c)**4

      air(p) =   emv(p) * (1._r8+(1._r8-emv(p))*(1._r8-emg(c))) * forc_lwrad(g)
      bir(p) = - (2._r8-emv(p)*(1._r8-emg(c))) * emv(p) * sb
      cir(p) =   emv(p)*emg(c)*sb

      ! Saturated vapor pressure, specific humidity, and their derivatives
      ! at the leaf surface

      call QSat (t_veg(p), forc_pbot(g), el(p), deldT, qsatl(p), qsatldT(p))

      ! Determine atmospheric co2 and o2

      co2(p) = forc_pco2(g)
      o2(p)  = forc_po2(g)
      
      if (use_c13) then
         ! 4/14/05: PET
         ! Adding isotope code
         c13o2(p) = forc_pc13o2(g)
      end if
      
      ! Initialize flux profile

      nmozsgn(p) = 0

      taf(p) = (t_grnd(c) + thm(p))/2._r8
      qaf(p) = (forc_q(g)+qg(c))/2._r8

      ur(p) = max(1.0_r8,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
      dth(p) = thm(p)-taf(p)
      dqh(p) = forc_q(g)-qaf(p)
      delq(p) = qg(c) - qaf(p)
      dthv(p) = dth(p)*(1._r8+0.61_r8*forc_q(g))+0.61_r8*forc_th(g)*dqh(p)
      zldis(p) = forc_hgt_u_pft(p) - displa(p)

      ! Check to see if the forcing height is below the canopy height
      if (zldis(p) < 0._r8) then
         found = .true.
         index = p
      end if

   end do

   if (found) then
      write(iulog,*)'Error: Forcing height is below canopy height for pft index ',index
      call endrun()
   end if

   do f = 1, fn
      p = filterp(f)
      c = pcolumn(p)

      ! Initialize Monin-Obukhov length and wind speed

      call MoninObukIni(ur(p), thv(c), dthv(p), zldis(p), z0mv(p), um(p), obu(p))

   end do

   ! Set counter for leaf temperature iteration (itlef)

   itlef = 0    
   fnorig = fn
   fporig(1:fn) = filterp(1:fn)

   ! Make copies so that array sections are not passed in function calls to friction velocity
   
   do f = 1, fn
      p = filterp(f)
      displa_loc(p) = displa(p)
      z0mv_loc(p) = z0mv(p)
      z0hv_loc(p) = z0hv(p)
      z0qv_loc(p) = z0qv(p)
   end do

   ! Begin stability iteration

   ITERATION : do while (itlef <= itmax .and. fn > 0)

      ! Determine friction velocity, and potential temperature and humidity
      ! profiles of the surface boundary layer

      call FrictionVelocity (lbp, ubp, fn, filterp, &
                             displa_loc, z0mv_loc, z0hv_loc, z0qv_loc, &
                             obu, itlef+1, ur, um, ustar, &
                             temp1, temp2, temp12m, temp22m, fm)

      do f = 1, fn
         p = filterp(f)
         c = pcolumn(p)
         g = pgridcell(p)

         tlbef(p) = t_veg(p)
         del2(p) = del(p)

         ! Determine aerodynamic resistances

         ram1(p)  = 1._r8/(ustar(p)*ustar(p)/um(p))
         rah(p,1) = 1._r8/(temp1(p)*ustar(p))
         raw(p,1) = 1._r8/(temp2(p)*ustar(p))

         ! Bulk boundary layer resistance of leaves

         uaf(p) = um(p)*sqrt( 1._r8/(ram1(p)*um(p)) )
         cf  = 0.01_r8/(sqrt(uaf(p))*sqrt(dleaf(ivt(p))))
         rb(p)  = 1._r8/(cf*uaf(p))
         rb1(p) = rb(p)
  
         ! Parameterization for variation of csoilc with canopy density from
         ! X. Zeng, University of Arizona

         w = exp(-(elai(p)+esai(p)))

         ! changed by K.Sakaguchi from here
         ! transfer coefficient over bare soil is changed to a local variable
         ! just for readability of the code (from line 680)
         csoilb = (vkc/(0.13_r8*(z0mg(c)*uaf(p)/1.5e-5_r8)**0.45_r8))

         !compute the stability parameter for ricsoilc  ("S" in Sakaguchi&Zeng,2008)

         ri = ( grav*htop(p) * (taf(p) - t_grnd(c)) ) / (taf(p) * uaf(p) **2.00_r8)

         !! modify csoilc value (0.004) if the under-canopy is in stable condition

         if ( (taf(p) - t_grnd(c) ) > 0._r8) then
               ! decrease the value of csoilc by dividing it with (1+gamma*min(S, 10.0))
               ! ria ("gmanna" in Sakaguchi&Zeng, 2008) is a constant (=0.5)
               ricsoilc = csoilc / (1.00_r8 + ria*min( ri, 10.0_r8) )
               csoilcn = csoilb*w + ricsoilc*(1._r8-w)
         else
              csoilcn = csoilb*w + csoilc*(1._r8-w)
         end if

         !! Sakaguchi changes for stability formulation ends here

         rah(p,2) = 1._r8/(csoilcn*uaf(p))
         raw(p,2) = rah(p,2)

         ! Stomatal resistances for sunlit and shaded fractions of canopy.
         ! Done each iteration to account for differences in eah, tv.

         svpts(p) = el(p)                         ! pa
         eah(p) = forc_pbot(g) * qaf(p) / 0.622_r8   ! pa
      end do

      ! 4/25/05, PET: Now calling the sun/shade version of Stomata by default
      call Stomata (fn, filterp, lbp, ubp, svpts, eah, o2, co2, rb, dayl_factor, phase='sun')
      call Stomata (fn, filterp, lbp, ubp, svpts, eah, o2, co2, rb, dayl_factor, phase='sha')

      do f = 1, fn
         p = filterp(f)
         c = pcolumn(p)
         g = pgridcell(p)

         ! Sensible heat conductance for air, leaf and ground
         ! Moved the original subroutine in-line...

         wta    = 1._r8/rah(p,1)             ! air
         wtl    = (elai(p)+esai(p))/rb(p)    ! leaf
         wtg(p) = 1._r8/rah(p,2)             ! ground
         wtshi  = 1._r8/(wta+wtl+wtg(p))

         wtl0(p) = wtl*wtshi         ! leaf
         wtg0    = wtg(p)*wtshi      ! ground
         wta0(p) = wta*wtshi         ! air

         wtga    = wta0(p)+wtg0      ! ground + air
         wtal(p) = wta0(p)+wtl0(p)   ! air + leaf

         ! Fraction of potential evaporation from leaf

         if (fdry(p) .gt. 0._r8) then
            rppdry  = fdry(p)*rb(p)*(laisun(p)/(rb(p)+rssun(p)) + &
                                     laisha(p)/(rb(p)+rssha(p)))/elai(p)
         else
            rppdry = 0._r8
         end if
         efpot = forc_rho(g)*wtl*(qsatl(p)-qaf(p))

         if (efpot > 0._r8) then
            if (btran(p) > btran0) then
               qflx_tran_veg(p) = efpot*rppdry
               rpp = rppdry + fwet(p)
            else
               !No transpiration if btran below 1.e-10
               rpp = fwet(p)
               qflx_tran_veg(p) = 0._r8
            end if
            !Check total evapotranspiration from leaves
            rpp = min(rpp, (qflx_tran_veg(p)+h2ocan(p)/dtime)/efpot)
         else
            !No transpiration if potential evaporation less than zero
            rpp = 1._r8
            qflx_tran_veg(p) = 0._r8
         end if

         ! Update conductances for changes in rpp
         ! Latent heat conductances for ground and leaf.
         ! Air has same conductance for both sensible and latent heat.
         ! Moved the original subroutine in-line...

         wtaq    = frac_veg_nosno(p)/raw(p,1)                        ! air
         wtlq    = frac_veg_nosno(p)*(elai(p)+esai(p))/rb(p) * rpp   ! leaf

         !Litter layer resistance. Added by K.Sakaguchi
         snowdp_c = z_dl ! critical depth for 100% litter burial by snow (=litter thickness)
         fsno_dl = snowdp(c)/snowdp_c    ! effective snow cover for (dry)plant litter
         elai_dl = lai_dl*(1._r8 - min(fsno_dl,1._r8)) ! exposed (dry)litter area index
         rdl = ( 1._r8 - exp(-elai_dl) ) / ( 0.004_r8*uaf(p)) ! dry litter layer resistance

         ! add litter resistance and Lee and Pielke 1992 beta
         if (delq(p) .lt. 0._r8) then  !dew. Do not apply beta for negative flux (follow old rsoil)
            wtgq(p) = frac_veg_nosno(p)/(raw(p,2)+rdl)
         else
            wtgq(p) = soilbeta(c)*frac_veg_nosno(p)/(raw(p,2)+rdl)
         end if

         wtsqi   = 1._r8/(wtaq+wtlq+wtgq(p))

         wtgq0    = wtgq(p)*wtsqi      ! ground
         wtlq0(p) = wtlq*wtsqi         ! leaf
         wtaq0(p) = wtaq*wtsqi         ! air

         wtgaq    = wtaq0(p)+wtgq0     ! air + ground
         wtalq(p) = wtaq0(p)+wtlq0(p)  ! air + leaf

         dc1 = forc_rho(g)*cpair*wtl
         dc2 = hvap*forc_rho(g)*wtlq

         efsh   = dc1*(wtga*t_veg(p)-wtg0*t_grnd(c)-wta0(p)*thm(p))
         efe(p) = dc2*(wtgaq*qsatl(p)-wtgq0*qg(c)-wtaq0(p)*forc_q(g))

         ! Evaporation flux from foliage

         erre = 0._r8
         if (efe(p)*efeb(p) < 0._r8) then
            efeold = efe(p)
            efe(p)  = 0.1_r8*efeold
            erre = efe(p) - efeold
         end if
         dt_veg(p) = (sabv(p) + air(p) + bir(p)*t_veg(p)**4 + &
              cir(p)*t_grnd(c)**4 - efsh - efe(p)) / &
              (- 4._r8*bir(p)*t_veg(p)**3 +dc1*wtga +dc2*wtgaq*qsatldT(p))
         t_veg(p) = tlbef(p) + dt_veg(p)
         dels = dt_veg(p)
         del(p)  = abs(dels)
         err(p) = 0._r8
         if (del(p) > delmax) then
            dt_veg(p) = delmax*dels/del(p)
            t_veg(p) = tlbef(p) + dt_veg(p)
            err(p) = sabv(p) + air(p) + bir(p)*tlbef(p)**3*(tlbef(p) + &
                 4._r8*dt_veg(p)) + cir(p)*t_grnd(c)**4 - &
                 (efsh + dc1*wtga*dt_veg(p)) - (efe(p) + &
                 dc2*wtgaq*qsatldT(p)*dt_veg(p))
         end if

         ! Fluxes from leaves to canopy space
         ! "efe" was limited as its sign changes frequently.  This limit may
         ! result in an imbalance in "hvap*qflx_evap_veg" and
         ! "efe + dc2*wtgaq*qsatdt_veg"

         efpot = forc_rho(g)*wtl*(wtgaq*(qsatl(p)+qsatldT(p)*dt_veg(p)) &
            -wtgq0*qg(c)-wtaq0(p)*forc_q(g))
         qflx_evap_veg(p) = rpp*efpot
         
         ! Calculation of evaporative potentials (efpot) and
         ! interception losses; flux in kg m**-2 s-1.  ecidif
         ! holds the excess energy if all intercepted water is evaporated
         ! during the timestep.  This energy is later added to the
         ! sensible heat flux.

         ecidif = 0._r8
         if (efpot > 0._r8 .and. btran(p) > btran0) then
            qflx_tran_veg(p) = efpot*rppdry
         else
            qflx_tran_veg(p) = 0._r8
         end if
         ecidif = max(0._r8, qflx_evap_veg(p)-qflx_tran_veg(p)-h2ocan(p)/dtime)
         qflx_evap_veg(p) = min(qflx_evap_veg(p),qflx_tran_veg(p)+h2ocan(p)/dtime)

         ! The energy loss due to above two limits is added to
         ! the sensible heat flux.

         eflx_sh_veg(p) = efsh + dc1*wtga*dt_veg(p) + err(p) + erre + hvap*ecidif

         ! Re-calculate saturated vapor pressure, specific humidity, and their
         ! derivatives at the leaf surface

         call QSat(t_veg(p), forc_pbot(g), el(p), deldT, qsatl(p), qsatldT(p))

         ! Update vegetation/ground surface temperature, canopy air
         ! temperature, canopy vapor pressure, aerodynamic temperature, and
         ! Monin-Obukhov stability parameter for next iteration.

         taf(p) = wtg0*t_grnd(c) + wta0(p)*thm(p) + wtl0(p)*t_veg(p)
         qaf(p) = wtlq0(p)*qsatl(p) + wtgq0*qg(c) + forc_q(g)*wtaq0(p)

         ! Update Monin-Obukhov length and wind speed including the
         ! stability effect

         dth(p) = thm(p)-taf(p)
         dqh(p) = forc_q(g)-qaf(p)
         delq(p) = wtalq(p)*qg(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(g)

         tstar = temp1(p)*dth(p)
         qstar = temp2(p)*dqh(p)

         thvstar = tstar*(1._r8+0.61_r8*forc_q(g)) + 0.61_r8*forc_th(g)*qstar
         zeta = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*thv(c))

         if (zeta >= 0._r8) then     !stable
            zeta = min(2._r8,max(zeta,0.01_r8))
            um(p) = max(ur(p),0.1_r8)
         else                     !unstable
            zeta = max(-100._r8,min(zeta,-0.01_r8))
            wc = beta*(-grav*ustar(p)*thvstar*zii/thv(c))**0.333_r8
            um(p) = sqrt(ur(p)*ur(p)+wc*wc)
         end if
         obu(p) = zldis(p)/zeta

         if (obuold(p)*obu(p) < 0._r8) nmozsgn(p) = nmozsgn(p)+1
         if (nmozsgn(p) >= 4) obu(p) = zldis(p)/(-0.01_r8)
         obuold(p) = obu(p)

      end do   ! end of filtered pft loop

      ! Test for convergence

      itlef = itlef+1
      if (itlef > itmin) then
         do f = 1, fn
            p = filterp(f)
            dele(p) = abs(efe(p)-efeb(p))
            efeb(p) = efe(p)
            det(p)  = max(del(p),del2(p))
         end do
         fnold = fn
         fn = 0
         do f = 1, fnold
            p = filterp(f)
            if (.not. (det(p) < dtmin .and. dele(p) < dlemin)) then
               fn = fn + 1
               filterp(fn) = p
            end if
         end do
      end if

   end do ITERATION     ! End stability iteration

   fn = fnorig
   filterp(1:fn) = fporig(1:fn)

   do f = 1, fn
      p = filterp(f)
      c = pcolumn(p)
      g = pgridcell(p)

      ! Energy balance check in canopy

      err(p) = sabv(p) + air(p) + bir(p)*tlbef(p)**3*(tlbef(p) + 4._r8*dt_veg(p)) &
         + cir(p)*t_grnd(c)**4 - eflx_sh_veg(p) - hvap*qflx_evap_veg(p)

      ! Fluxes from ground to canopy space

      delt    = wtal(p)*t_grnd(c)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)
      taux(p) = -forc_rho(g)*forc_u(g)/ram1(p)
      tauy(p) = -forc_rho(g)*forc_v(g)/ram1(p)
      eflx_sh_grnd(p) = cpair*forc_rho(g)*wtg(p)*delt
      qflx_evap_soi(p) = forc_rho(g)*wtgq(p)*delq(p)

      ! 2 m height air temperature

      t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))
      t_ref2m_r(p) = t_ref2m(p)

      ! 2 m height specific humidity

      q_ref2m(p) = forc_q(g) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

      ! 2 m height relative humidity

      call QSat(t_ref2m(p), forc_pbot(g), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)
      rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)
      rh_ref2m_r(p) = rh_ref2m(p)

      ! Downward longwave radiation below the canopy

      dlrad(p) = (1._r8-emv(p))*emg(c)*forc_lwrad(g) + &
         emv(p)*emg(c)*sb*tlbef(p)**3*(tlbef(p) + 4._r8*dt_veg(p))

      ! Upward longwave radiation above the canopy

      ulrad(p) = ((1._r8-emg(c))*(1._r8-emv(p))*(1._r8-emv(p))*forc_lwrad(g) &
         + emv(p)*(1._r8+(1._r8-emg(c))*(1._r8-emv(p)))*sb*tlbef(p)**3*(tlbef(p) + &
         4._r8*dt_veg(p)) + emg(c)*(1._r8-emv(p))*sb*t_grnd(c)**4)

      ! Derivative of soil energy flux with respect to soil temperature

      cgrnds(p) = cgrnds(p) + cpair*forc_rho(g)*wtg(p)*wtal(p)
      cgrndl(p) = cgrndl(p) + forc_rho(g)*wtgq(p)*wtalq(p)*dqgdT(c)
      cgrnd(p)  = cgrnds(p) + cgrndl(p)*htvp(c)

      ! Update dew accumulation (kg/m2)
      
      h2ocan(p) = max(0._r8,h2ocan(p)+(qflx_tran_veg(p)-qflx_evap_veg(p))*dtime)
      
      ! total photosynthesis

      fpsn(p) = psnsun(p)*laisun(p) + psnsha(p)*laisha(p)
      
      if (use_cn) then
         if (use_c13) then
            ! 4/14/05: PET
            ! Adding isotope code
            rc13_canair(p) = c13o2(p)/(co2(p)-c13o2(p))
            rc13_psnsun(p) = rc13_canair(p)/alphapsnsun(p)
            rc13_psnsha(p) = rc13_canair(p)/alphapsnsha(p)
            c13_psnsun(p) = psnsun(p) * (rc13_psnsun(p)/(1._r8+rc13_psnsun(p)))
            c13_psnsha(p) = psnsha(p) * (rc13_psnsha(p)/(1._r8+rc13_psnsha(p)))
            !write(iulog,*) p,ivt(p),btran(p),psnsun(p),psnsha(p),alphapsnsun(p),alphapsnsha(p)
         end if
      end if
      
   end do

   ! Filter out pfts which have small energy balance errors; report others

   fnold = fn
   fn = 0
   do f = 1, fnold
      p = filterp(f)
      if (abs(err(p)) > 0.1_r8) then
         fn = fn + 1
         filterp(fn) = p
      end if
   end do

   do f = 1, fn
      p = filterp(f)
      write(iulog,*) 'energy balance in canopy ',p,', err=',err(p)
   end do

   end subroutine CanopyFluxes

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Stomata
!
! !INTERFACE:
   subroutine Stomata (fn, filterp, lbp, ubp, ei, ea, o2, co2, rb, dayl_factor, phase)
!
! !DESCRIPTION: 
! Leaf stomatal resistance and leaf photosynthesis. Modifications for CN code.

! !REVISION HISTORY:
! 22 January 2004: Created by Peter Thornton
! 4/14/05: Peter Thornton: Converted Ci from local variable to pps struct member
!    now returns cisun or cisha per pft as implicit output argument.
!    Also sets alphapsnsun and alphapsnsha.
! 4/25/05, Peter Thornton: Adopted as the default code for CLM, together with
!   modifications for sun/shade canopy.  Renamed from StomataCN to Stomata,
!   and eliminating the older Stomata subroutine
! 3/6/09: Peter Thornton; added dayl_factor control on Vcmax, from Bill Bauerle

! !USES:
     use shr_kind_mod , only : r8 => shr_kind_r8
     use shr_const_mod, only : SHR_CONST_TKFRZ, SHR_CONST_RGAS
     use clmtype
     use clm_atmlnd   , only : clm_a2l
     use spmdMod      , only: masterproc
     use pftvarcon    , only : nbrdlf_dcd_tmp_shrub
     use pftvarcon    , only : nsoybean, npcropmin
!
! !ARGUMENTS:
     implicit none
     integer , intent(in)    :: fn                 ! size of pft filter
     integer , intent(in)    :: filterp(fn)        ! pft filter
     integer , intent(in)    :: lbp, ubp           ! pft bounds
     real(r8), intent(in)    :: ei(lbp:ubp)        ! vapor pressure inside leaf (sat vapor press at tl) (pa)
     real(r8), intent(in)    :: ea(lbp:ubp)        ! vapor pressure of canopy air (pa)
     real(r8), intent(in)    :: o2(lbp:ubp)        ! atmospheric o2 concentration (pa)
     real(r8), intent(in)    :: co2(lbp:ubp)       ! atmospheric co2 concentration (pa)
     real(r8), intent(inout) :: rb(lbp:ubp)        ! boundary layer resistance (s/m)
     real(r8), intent(in)    :: dayl_factor(lbp:ubp) ! scalar (0-1) for daylength
     character(len=*), intent(in) :: phase         ! 'sun' or 'sha'
!
! !CALLED FROM:
! subroutine CanopyFluxes in this module
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in variables
! new ecophys variables (leafcn, flnr) added 1/26/04
!
     integer , pointer :: pcolumn(:)     ! pft's column index
     integer , pointer :: pgridcell(:)   ! pft's gridcell index
     integer , pointer :: ivt(:)         ! pft vegetation type
     real(r8), pointer :: qe25(:)        ! quantum efficiency at 25C (umol CO2 / umol photon)
     real(r8), pointer :: c3psn(:)       ! photosynthetic pathway: 0. = c4, 1. = c3
     real(r8), pointer :: mp(:)          ! slope of conductance-to-photosynthesis relationship
     real(r8), pointer :: tgcm(:)        ! air temperature at agcm reference height (kelvin)
     real(r8), pointer :: forc_pbot(:)   ! atmospheric pressure (Pa)
     real(r8), pointer :: tl(:)          ! leaf temperature (Kelvin)
     real(r8), pointer :: btran(:)       ! soil water transpiration factor (0 to 1)
     real(r8), pointer :: apar(:)        ! par absorbed per unit lai (w/m**2)
     real(r8), pointer :: leafcn(:)      ! leaf C:N (gC/gN)
     real(r8), pointer :: flnr(:)        ! fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
     real(r8), pointer :: sla(:)         ! specific leaf area, projected area basis (m^2/gC)
     real(r8), pointer :: fnitr(:)       ! foliage nitrogen limitation factor (-)
!
! local pointers to implicit inout variables
!
     real(r8), pointer :: rs(:)          ! leaf stomatal resistance (s/m)
     real(r8), pointer :: psn(:)         ! foliage photosynthesis (umol co2 /m**2/ s) [always +]
     real(r8), pointer :: ci(:)          ! intracellular leaf CO2 (Pa)
     real(r8), pointer :: alphapsn(:)    ! 13C fractionation factor for PSN ()
!
! local pointers to implicit out variables
!
     real(r8), pointer :: lnc(:)         ! leaf N concentration per unit projected LAI (gN leaf/m^2)
     real(r8), pointer :: vcmx(:)        ! maximum rate of carboxylation (umol co2/m**2/s)
!
!
! !LOCAL VARIABLES:
!EOP
!
     real(r8), parameter :: mpe = 1.e-6_r8   ! prevents overflow error if division by zero
     integer , parameter :: niter = 3     ! number of iterations
     integer  :: f,p,c,g ! indices
     integer  :: iter    ! iteration index
     real(r8) :: ab      ! used in statement functions
     real(r8) :: bc      ! used in statement functions
     real(r8) :: f1      ! generic temperature response (statement function)
     real(r8) :: f2      ! generic temperature inhibition (statement function)
     real(r8) :: tc      ! leaf temperature (degree celsius)
     real(r8) :: cs      ! co2 concentration at leaf surface (pa)
     real(r8) :: kc      ! co2 michaelis-menten constant (pa)
     real(r8) :: ko      ! o2 michaelis-menten constant (pa)
     real(r8) :: atmp    ! intermediate calculations for rs
     real(r8) :: btmp    ! intermediate calculations for rs
     real(r8) :: ctmp    ! intermediate calculations for rs
     real(r8) :: q       ! intermediate calculations for rs
     real(r8) :: r1,r2   ! roots for rs
     real(r8) :: ppf     ! absorb photosynthetic photon flux (umol photons/m**2/s)
     real(r8) :: wc      ! rubisco limited photosynthesis (umol co2/m**2/s)
     real(r8) :: wj      ! light limited photosynthesis (umol co2/m**2/s)
     real(r8) :: we      ! export limited photosynthesis (umol co2/m**2/s)
     real(r8) :: cp      ! co2 compensation point (pa)
     real(r8) :: awc     ! intermediate calcuation for wc
     real(r8) :: j       ! electron transport (umol co2/m**2/s)
     real(r8) :: cea     ! constrain ea or else model blows up
     real(r8) :: cf      ! s m**2/umol -> s/m
     real(r8) :: rsmax0  ! maximum stomatal resistance [s/m]
     real(r8) :: kc25    ! co2 michaelis-menten constant at 25c (pa)
     real(r8) :: akc     ! q10 for kc25
     real(r8) :: ko25    ! o2 michaelis-menten constant at 25c (pa)
     real(r8) :: ako     ! q10 for ko25
     real(r8) :: bp      ! minimum leaf conductance (umol/m**2/s)
     ! additional variables for new treatment of Vcmax, Peter Thornton, 1/26/04
     real(r8) :: act25   ! (umol/mgRubisco/min) Rubisco activity at 25 C
     real(r8) :: act     ! (umol/mgRubisco/min) Rubisco activity
     real(r8) :: q10act  ! (DIM) Q_10 for Rubisco activity
     real(r8) :: fnr     ! (gRubisco/gN in Rubisco)
!------------------------------------------------------------------------------

     ! Set statement functions

     f1(ab,bc) = ab**((bc-25._r8)/10._r8)
     f2(ab) = 1._r8 + exp((-2.2e05_r8+710._r8*(ab+SHR_CONST_TKFRZ))/(SHR_CONST_RGAS*0.001_r8*(ab+SHR_CONST_TKFRZ)))

     ! Assign local pointers to derived type members (pft-level)

     pcolumn   => pft%column
     pgridcell => pft%gridcell
     ivt       => pft%itype
     tl        => pes%t_veg
     btran     => pps%btran
     if (phase == 'sun') then
        apar   => pef%parsun
        rs     => pps%rssun
        psn    => pcf%psnsun
        ci     => pps%cisun

        alphapsn  => pps%alphapsnsun

        sla    => pps%slasun
        lnc    => pps%lncsun   
        vcmx   => pps%vcmxsun   
     else if (phase == 'sha') then
        apar   => pef%parsha
        rs     => pps%rssha
        psn    => pcf%psnsha
        ci     => pps%cisha
        sla    => pps%slasha   

        alphapsn  => pps%alphapsnsha

        lnc    => pps%lncsha   
        vcmx   => pps%vcmxsha
     end if

     ! Assign local pointers to derived type members (gridcell-level)

     forc_pbot => clm_a2l%forc_pbot

     ! Assign local pointers to derived type members (column-level)

     tgcm        => pes%thm

     ! Assign local pointers to pft constants
     ! new ecophys constants added 1/26/04

     qe25      => pftcon%qe25
     c3psn     => pftcon%c3psn
     mp        => pftcon%mp
     leafcn    => pftcon%leafcn
     flnr      => pftcon%flnr
     fnitr     => pftcon%fnitr

     ! Set constant values

     kc25  = 30._r8
     akc   = 2.1_r8
     ko25  = 30000._r8
     ako   = 1.2_r8
     bp    = 2000._r8

     ! New constants for CN code, added 1/26/04

     act25 = 3.6_r8
     q10act = 2.4_r8
     fnr = 7.16_r8
     
     ! Convert rubisco activity units from umol/mgRubisco/min -> umol/gRubisco/s

     act25 = act25 * 1000.0_r8 / 60.0_r8

     do f = 1, fn
        p = filterp(f)
        c = pcolumn(p)
        g = pgridcell(p)

        ! Initialize rs=rsmax and psn=0 because calculations are performed only
        ! when apar > 0, in which case rs <= rsmax and psn >= 0
        ! Set constants

        rsmax0 = 2.e4_r8
        cf = forc_pbot(g)/(SHR_CONST_RGAS*0.001_r8*tgcm(p))*1.e06_r8
        if (apar(p) <= 0._r8) then          ! night time
           rs(p) = min(rsmax0, 1._r8/bp * cf)
           psn(p) = 0._r8
           lnc(p) = 0._r8
           vcmx(p) = 0._r8
           if (use_c13) then
              alphapsn(p) = 1._r8
           end if
        else                             ! day time
           tc = tl(p) - SHR_CONST_TKFRZ
           ppf = 4.6_r8 * apar(p)                  
           j = ppf * qe25(ivt(p))
           kc = kc25 * f1(akc,tc)       
           ko = ko25 * f1(ako,tc)
           awc = kc * (1._r8+o2(p)/ko)
           cp = 0.5_r8*kc/ko*o2(p)*0.21_r8

           ! Modification for shrubs proposed by X.D.Z
           ! Why does he prefer this line here instead of in subr.
           ! CanopyFluxes? (slevis)
           ! Equivalent modification for soy following AgroIBIS
           if (use_cndv) then
              if (ivt(p) == nbrdlf_dcd_tmp_shrub) btran(p) = min(1._r8, btran(p) * 3.33_r8)
           end if
           if (ivt(p) == nsoybean) btran(p) = min(1._r8, btran(p) * 1.25_r8)
           
           ! new calculations for vcmax, 1/26/04
           lnc(p) = 1._r8 / (sla(p) * leafcn(ivt(p)))
           act = act25 * f1(q10act,tc)
           if (use_cn) then
              if ( ivt(p) < npcropmin )then
                 vcmx(p) = lnc(p) * flnr(ivt(p)) * fnr * act / f2(tc) * btran(p) * &
                      dayl_factor(p)
              else
                 vcmx(p) = 101._r8 * f1(q10act,tc) / f2(tc) * btran(p) * dayl_factor(p)
              end if
           else
              vcmx(p) = lnc(p) * flnr(ivt(p)) * fnr * act / f2(tc) * btran(p) * &
                   dayl_factor(p) * fnitr(ivt(p))
           end if
           
           ! First guess ci

           ci(p) = 0.7_r8*co2(p)*c3psn(ivt(p)) + 0.4_r8*co2(p)*(1._r8-c3psn(ivt(p)))

           ! rb: s/m -> s m**2 / umol

           rb(p) = rb(p)/cf 

           ! Constrain ea

           cea = max(0.25_r8*ei(p)*c3psn(ivt(p))+0.40_r8*ei(p)*(1._r8-c3psn(ivt(p))), min(ea(p),ei(p)) ) 

           ! ci iteration for 'actual' photosynthesis

           do iter = 1, niter
              wj = max(ci(p)-cp,0._r8)*j/(ci(p)+2._r8*cp)*c3psn(ivt(p)) + j*(1._r8-c3psn(ivt(p)))
              wc = max(ci(p)-cp,0._r8)*vcmx(p)/(ci(p)+awc)*c3psn(ivt(p)) + vcmx(p)*(1._r8-c3psn(ivt(p)))
              we = 0.5_r8*vcmx(p)*c3psn(ivt(p)) + 4000._r8*vcmx(p)*ci(p)/forc_pbot(g)*(1._r8-c3psn(ivt(p))) 
              psn(p) = min(wj,wc,we) 
              cs = max( co2(p)-1.37_r8*rb(p)*forc_pbot(g)*psn(p), mpe )
              atmp = mp(ivt(p))*psn(p)*forc_pbot(g)*cea / (cs*ei(p)) + bp
              btmp = ( mp(ivt(p))*psn(p)*forc_pbot(g)/cs + bp ) * rb(p) - 1._r8
              ctmp = -rb(p)
              if (btmp >= 0._r8) then
                 q = -0.5_r8*( btmp + sqrt(btmp*btmp-4._r8*atmp*ctmp) )
              else
                 q = -0.5_r8*( btmp - sqrt(btmp*btmp-4._r8*atmp*ctmp) )
              end if
              r1 = q/atmp
              r2 = ctmp/q
              rs(p) = max(r1,r2)
              ci(p) = max( cs-psn(p)*forc_pbot(g)*1.65_r8*rs(p), 0._r8 )
           end do

           ! rs, rb:  s m**2 / umol -> s/m 

           rs(p) = min(rsmax0, rs(p)*cf)
           rb(p) = rb(p) * cf 
           
           if (use_c13) then
              ! 4/14/05: PET
              ! Adding isotope code
              alphapsn(p) = 1._r8 + (((c3psn(ivt(p)) * (4.4_r8 + (22.6_r8*(ci(p)/co2(p))))) + &
                   ((1._r8 - c3psn(ivt(p))) * 4.4_r8))/1000._r8)
              !alphapsn(p) = 1._r8
              !write(iulog,*) 'in StomataCN ',p,ivt(p),c3psn(ivt(p)),ci(p),co2(p),alphapsn(p)
           end if
           
        end if

     end do

  end subroutine Stomata

end module CanopyFluxesMod
