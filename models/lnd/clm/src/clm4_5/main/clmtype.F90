module clmtype

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clmtype
!
! !DESCRIPTION: 
! Define derived type hierarchy. Includes declaration of
! the clm derived type and 1d mapping arrays. 
!
! -------------------------------------------------------- 
! gridcell types can have values of 
! -------------------------------------------------------- 
!   1 => default
! -------------------------------------------------------- 
! landunits types can have values of (see clm_varcon.F90)
! -------------------------------------------------------- 
!   1  => (istsoil)    soil (vegetated or bare soil landunit)
!   2  => (istice)     land ice
!   3  => (istdlak)    deep lake
!   4  => (istslak) shall lake (not currently implemented; SLake implementation has variable depth)
!   5  => (istwet)     wetland
!   6  => (isturb)     urban
!   7  => (istice_mec) land ice (multiple elevation classes) 
!   8  => (istcrop)    crop (only for crop configuration)
! -------------------------------------------------------- 
! column types can have values of
! -------------------------------------------------------- 
!   1  => (istsoil)          soil (vegetated or bare soil)
!   2  => (istice)           land ice
!   3  => (istdlak)          deep lake
!   4  => (istslak)          shallow lake 
!   5  => (istwet)           wetland
!   7  => (istice_mec)       land ice (multiple elevation classes)   
!   61 => (icol_roof)        urban roof
!   62 => (icol_sunwall)     urban sunwall
!   63 => (icol_shadewall)   urban shadewall
!   64 => (icol_road_imperv) urban impervious road
!   65 => (icol_road_perv)   urban pervious road
! -------------------------------------------------------- 
! pft types can have values of
! -------------------------------------------------------- 
!   0  => not vegetated
!   1  => needleleaf evergreen temperate tree
!   2  => needleleaf evergreen boreal tree
!   3  => needleleaf deciduous boreal tree
!   4  => broadleaf evergreen tropical tree
!   5  => broadleaf evergreen temperate tree
!   6  => broadleaf deciduous tropical tree
!   7  => broadleaf deciduous temperate tree
!   8  => broadleaf deciduous boreal tree
!   9  => broadleaf evergreen shrub
!   10 => broadleaf deciduous temperate shrub
!   11 => broadleaf deciduous boreal shrub
!   12 => c3 arctic grass
!   13 => c3 non-arctic grass
!   14 => c4 grass
!   15 => c3_crop
!   16 => c3_irrigated
!   17 => corn
!   18 => irrigated corn
!   19 => spring temperate cereal
!   20 => irrigated spring temperate cereal
!   21 => winter temperate cereal
!   22 => irrigated winter temperate cereal
!   23 => soybean
!   24 => irrigated soybean
! -------------------------------------------------------- 
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use domainMod   , only: domain_type
!
! !PUBLIC TYPES:
  implicit none

  private
!                              
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!F. Li and S. Levis (11/06/12)
!*******************************************************************************
!----------------------------------------------------
! Begin definition of conservation check structures
!----------------------------------------------------
! energy balance structure
!----------------------------------------------------
type, public :: energy_balance_type
   real(r8), pointer :: errsoi(:)        !soil/lake energy conservation error (W/m**2)
   real(r8), pointer :: errseb(:)        !surface energy conservation error (W/m**2)
   real(r8), pointer :: errsol(:)        !solar radiation conservation error (W/m**2)
   real(r8), pointer :: errlon(:)        !longwave radiation conservation error (W/m**2)
end type energy_balance_type

!----------------------------------------------------
! water balance structure
!----------------------------------------------------
type, public :: water_balance_type
   real(r8), pointer :: begwb(:)         !water mass begining of the time step
   real(r8), pointer :: endwb(:)         !water mass end of the time step
   real(r8), pointer :: errh2o(:)        !water conservation error (mm H2O)
end type water_balance_type

!----------------------------------------------------
! carbon balance structure
!----------------------------------------------------
type, public :: carbon_balance_type
   real(r8), pointer :: begcb(:)         !carbon mass, beginning of time step (gC/m**2)
   real(r8), pointer :: endcb(:)         !carbon mass, end of time step (gC/m**2)
   real(r8), pointer :: errcb(:)         !carbon balance error for the timestep (gC/m**2)
end type carbon_balance_type

!----------------------------------------------------
! nitrogen balance structure
!----------------------------------------------------
type, public :: nitrogen_balance_type
   real(r8), pointer :: begnb(:)         !nitrogen mass, beginning of time step (gN/m**2)
   real(r8), pointer :: endnb(:)         !nitrogen mass, end of time step (gN/m**2)
   real(r8), pointer :: errnb(:)         !nitrogen balance error for the timestep (gN/m**2)
end type nitrogen_balance_type

!----------------------------------------------------
! End definition of conservation check structures
!----------------------------------------------------
!*******************************************************************************

!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the pft_type level
!----------------------------------------------------
! pft physical state variables structure
!----------------------------------------------------
type, public :: pft_pstate_type
   real(r8), pointer :: prec60(:)              !60-day running mean of tot. precipitation (mm/s) added by F. Li and S. Levis
   real(r8), pointer :: prec10(:)              !10-day running mean of tot. precipitation (mm/s) added by F. Li and S. Levis
   integer , pointer :: frac_veg_nosno(:)       !fraction of vegetation not covered by snow (0 OR 1) [-] 
   integer , pointer :: frac_veg_nosno_alb(:)   !fraction of vegetation not covered by snow (0 OR 1) [-] 
   real(r8), pointer :: emv(:)                  !vegetation emissivity
   real(r8), pointer :: z0mv(:)                 !roughness length over vegetation, momentum [m]
   real(r8), pointer :: z0hv(:)                 !roughness length over vegetation, sensible heat [m]
   real(r8), pointer :: z0qv(:)                 !roughness length over vegetation, latent heat [m]
   real(r8), pointer :: rootfr(:,:)             !fraction of roots in each soil layer  (nlevgrnd)
   real(r8), pointer :: rootr(:,:)              !effective fraction of roots in each soil layer  (nlevgrnd)
   real(r8), pointer :: rresis(:,:)             !root resistance by layer (0-1)  (nlevgrnd)
   real(r8), pointer :: dewmx(:)                !Maximum allowed dew [mm]
   real(r8), pointer :: rssun(:)                !sunlit stomatal resistance (s/m)
   real(r8), pointer :: rhal(:)                 !relative humidity of the canopy air vs leaf
   real(r8), pointer :: vpdal(:)                !vpd of the canopy air vs leaf
   real(r8), pointer :: rssha(:)                !shaded stomatal resistance (s/m)
   real(r8), pointer :: rssun_z(:,:)            !canopy layer: sunlit leaf stomatal resistance (s/m)
   real(r8), pointer :: rssha_z(:,:)            !canopy layer: shaded leaf stomatal resistance (s/m)
   real(r8), pointer :: laisun(:)               !sunlit projected leaf area index
   real(r8), pointer :: laisha(:)               !shaded projected leaf area index
   real(r8), pointer :: laisun_z(:,:)           !sunlit leaf area for canopy layer
   real(r8), pointer :: laisha_z(:,:)           !shaded leaf area for canopy layer
   real(r8), pointer :: btran(:)                !transpiration wetness factor (0 to 1)
   real(r8), pointer :: btran2(:)               ! root zone soil wetness factor (0 to 1) added by F. Li and S. Levis
   real(r8), pointer :: fsun(:)                 !sunlit fraction of canopy
   real(r8), pointer :: tlai(:)                 !one-sided leaf area index, no burying by snow
   real(r8), pointer :: tsai(:)                 !one-sided stem area index, no burying by snow
   real(r8), pointer :: elai(:)                 !one-sided leaf area index with burying by snow
   real(r8), pointer :: esai(:)                 !one-sided stem area index with burying by snow
   real(r8), pointer :: fwet(:)                 !fraction of canopy that is wet (0 to 1)
   real(r8), pointer :: fdry(:)                 !fraction of foliage that is green and dry [-] (new)
   real(r8), pointer :: dt_veg(:)               !change in t_veg, last iteration (Kelvin)
   real(r8), pointer :: htop(:)                 !canopy top (m)
   real(r8), pointer :: hbot(:)                 !canopy bottom (m)
   real(r8), pointer :: z0m(:)                  !momentum roughness length (m)
   real(r8), pointer :: displa(:)               !displacement height (m)
   real(r8), pointer :: albd(:,:)               !surface albedo (direct)                              (numrad)
   real(r8), pointer :: albi(:,:)               !surface albedo (diffuse)                             (numrad)
   real(r8), pointer :: fabd(:,:)               !flux absorbed by canopy per unit direct flux         (numrad)
   real(r8), pointer :: fabd_sun(:,:)           !flux absorbed by sunlit canopy per unit direct flux  (numrad)
   real(r8), pointer :: fabd_sha(:,:)           !flux absorbed by shaded canopy per unit direct flux  (numrad)
   real(r8), pointer :: fabi(:,:)               !flux absorbed by canopy per unit diffuse flux        (numrad)
   real(r8), pointer :: fabi_sun(:,:)           !flux absorbed by sunlit canopy per unit diffuse flux (numrad)
   real(r8), pointer :: fabi_sha(:,:)           !flux absorbed by shaded canopy per unit diffuse flux (numrad)
   real(r8), pointer :: ftdd(:,:)               !down direct flux below canopy per unit direct flx    (numrad)
   real(r8), pointer :: ftid(:,:)               !down diffuse flux below canopy per unit direct flx   (numrad)
   real(r8), pointer :: ftii(:,:)               !down diffuse flux below canopy per unit diffuse flx  (numrad)
   real(r8), pointer :: vcmaxcintsun(:)         ! leaf to canopy scaling coefficient, sunlit leaf vcmax
   real(r8), pointer :: vcmaxcintsha(:)         ! leaf to canopy scaling coefficient, shaded leaf vcmax
   integer , pointer :: ncan(:)                 !number of canopy layers
   integer , pointer :: nrad(:)                 !number of canopy layers, above snow for radiative transfer
   real(r8), pointer :: fabd_sun_z(:,:)         !absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
   real(r8), pointer :: fabd_sha_z(:,:)         !absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
   real(r8), pointer :: fabi_sun_z(:,:)         !absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
   real(r8), pointer :: fabi_sha_z(:,:)         !absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
   real(r8), pointer :: fsun_z(:,:)             !sunlit fraction of canopy layer
   real(r8), pointer :: tlai_z(:,:)             !tlai increment for canopy layer
   real(r8), pointer :: tsai_z(:,:)             !tsai increment for canopy layer
   real(r8), pointer :: u10(:)                  !10-m wind (m/s) (for dust model)
   real(r8), pointer :: u10_clm(:)              !10-m wind (m/s)
   real(r8), pointer :: va(:)                   !atmospheric wind speed plus convective velocity (m/s)
   real(r8), pointer :: ram1(:)                 !aerodynamical resistance (s/m)
   real(r8), pointer :: ram1_lake(:)            !aerodynamical resistance (s/m)
   integer , pointer :: burndate(:)             ! crop burn date
   real(r8), pointer :: rh_leaf(:)              !fractional humidity at leaf surface (dimensionless)
   real(r8), pointer :: rhaf(:)                 !fractional humidity of canopy air (dimensionless)
   real(r8), pointer :: fv(:)                   !friction velocity (m/s) (for dust model)
   real(r8), pointer :: forc_hgt_u_pft(:)       !wind forcing height (10m+z0m+d) (m)
   real(r8), pointer :: forc_hgt_t_pft(:)       !temperature forcing height (10m+z0m+d) (m)
   real(r8), pointer :: forc_hgt_q_pft(:)       !specific humidity forcing height (10m+z0m+d) (m)
   real(r8), pointer :: lfpftd(:)               ! decrease of pft weight (0-1) on the column for the timestep added by F. Li and S. Levis
   ! Variables for prognostic crop model
   real(r8), pointer :: hdidx(:)                ! cold hardening index?
   real(r8), pointer :: cumvd(:)                ! cumulative vernalization d?ependence?
   real(r8), pointer :: htmx(:)                 ! max hgt attained by a crop during yr (m)
   real(r8), pointer :: vf(:)                   ! vernalization factor for cereal
   real(r8), pointer :: gddmaturity(:)          ! growing degree days (gdd) needed to harvest (ddays)
   real(r8), pointer :: gdd0(:)                 ! growing degree-days base  0C from planting  (ddays)
   real(r8), pointer :: gdd8(:)                 ! growing degree-days base  8C from planting  (ddays)
   real(r8), pointer :: gdd10(:)                ! growing degree-days base 10C from planting  (ddays)
   real(r8), pointer :: gdd020(:)               ! 20-year average of gdd0                     (ddays)
   real(r8), pointer :: gdd820(:)               ! 20-year average of gdd8                     (ddays)
   real(r8), pointer :: gdd1020(:)              ! 20-year average of gdd10                    (ddays)
   real(r8), pointer :: gddplant(:)             ! accum gdd past planting date for crop       (ddays)
   real(r8), pointer :: gddtsoi(:)              ! growing degree-days from planting (top two soil layers) (ddays)
   real(r8), pointer :: huileaf(:)              ! heat unit index needed from planting to leaf emergence
   real(r8), pointer :: huigrain(:)             ! heat unit index needed to reach vegetative maturity
   real(r8), pointer :: aleafi(:)               ! saved leaf allocation coefficient from phase 2
   real(r8), pointer :: astemi(:)               ! saved stem allocation coefficient from phase 2
   real(r8), pointer :: aleaf(:)                ! leaf allocation coefficient
   real(r8), pointer :: astem(:)                ! stem allocation coefficient
   logical , pointer :: croplive(:)             ! Flag, true if planted, not harvested
   logical , pointer :: cropplant(:)            ! Flag, true if planted
   integer , pointer :: harvdate(:)             ! harvest date
                                                ! cropplant and harvdate could be 2D to facilitate rotation
   integer , pointer :: idop(:)                 ! date of planting
   integer , pointer :: peaklai(:)              ! 1: max allowed lai; 0: not at max
   real(r8), pointer :: vds(:)                  !deposition velocity term (m/s) (for dry dep SO4, NH4NO3)
   real(r8), pointer :: alphapsnsun(:) !sunlit 13c fractionation ([])
   real(r8), pointer :: alphapsnsha(:) !shaded 13c fractionation ([])
   real(r8), pointer :: sandfrac(:)    ! sand fraction
   real(r8), pointer :: clayfrac(:)    ! clay fraction
   ! for dry deposition of chemical tracers
   real(r8), pointer :: mlaidiff(:)    ! difference between lai month one and month two
   real(r8), pointer :: rb1(:)         ! aerodynamical resistance (s/m)
   real(r8), pointer :: annlai(:,:)    ! 12 months of monthly lai from input data set  

   ! New variable for methane code
#ifdef LCH4
   real(r8), pointer :: grnd_ch4_cond(:)    !tracer conductance for boundary layer [m/s]
   real(r8), pointer :: canopy_cond(:)    !tracer conductance for canopy [m/s]
#endif
   ! and vertical profiles for calculating fluxes
   real(r8), pointer :: leaf_prof(:,:)          ! (1/m) profile of leaves
   real(r8), pointer :: froot_prof(:,:)         ! (1/m) profile of fine roots
   real(r8), pointer :: croot_prof(:,:)         ! (1/m) profile of coarse roots
   real(r8), pointer :: stem_prof(:,:)          ! (1/m) profile of stems
end type pft_pstate_type

type, public :: pft_psynstate_type
   logical, pointer :: c3flag(:)                ! true if C3 and false if C4
   real(r8), pointer :: ac(:,:)                  ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
   real(r8), pointer :: aj(:,:)                  ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
   real(r8), pointer :: ap(:,:)                  ! product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
   real(r8), pointer :: ag(:,:)                  ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
   real(r8), pointer :: an(:,:)                  ! net leaf photosynthesis (umol CO2/m**2/s)   
   real(r8), pointer :: vcmax_z(:,:)             ! maximum rate of carboxylation (umol co2/m**2/s)
   real(r8), pointer :: cp(:)                    ! CO2 compensation point (Pa)
   real(r8), pointer :: kc(:)                    ! Michaelis-Menten constant for CO2 (Pa)
   real(r8), pointer :: ko(:)                    ! Michaelis-Menten constant for O2 (Pa)
   real(r8), pointer :: qe(:)                    ! quantum efficiency, used only for C4 (mol CO2 / mol photons)
   real(r8), pointer :: tpu_z(:,:)               ! triose phosphate utilization rate (umol CO2/m**2/s)
   real(r8), pointer :: kp_z(:,:)                ! initial slope of CO2 response curve (C4 plants)
   real(r8), pointer :: theta_cj(:)              ! empirical curvature parameter for ac, aj photosynthesis co-limitation
   real(r8), pointer :: bbb(:)                   ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
   real(r8), pointer :: mbb(:)                   ! Ball-Berry slope of conductance-photosynthesis relationship
   real(r8), pointer :: gs_mol(:,:)              ! leaf stomatal conductance (umol H2O/m**2/s)
   real(r8), pointer :: gb_mol(:)                ! leaf boundary layer conductance (umol H2O/m**2/s)
end type pft_psynstate_type
!----------------------------------------------------
! pft ecophysiological constants structure
!----------------------------------------------------
type, public :: pft_epc_type
   integer , pointer :: noveg(:)                !value for not vegetated
   integer , pointer :: tree(:)                 !tree or not?
   real(r8), pointer :: smpso(:)                !soil water potential at full stomatal opening (mm)
   real(r8), pointer :: smpsc(:)                !soil water potential at full stomatal closure (mm)
   real(r8), pointer :: fnitr(:)                !foliage nitrogen limitation factor (-)
   real(r8), pointer :: foln(:)                 !foliage nitrogen (%)
   real(r8), pointer :: dleaf(:)                !characteristic leaf dimension (m)
   real(r8), pointer :: c3psn(:)                !photosynthetic pathway: 0. = c4, 1. = c3
   real(r8), pointer :: xl(:)                   !leaf/stem orientation index
   real(r8), pointer :: rhol(:,:)               !leaf reflectance: 1=vis, 2=nir   (numrad)
   real(r8), pointer :: rhos(:,:)               !stem reflectance: 1=vis, 2=nir   (numrad)
   real(r8), pointer :: taul(:,:)               !leaf transmittance: 1=vis, 2=nir (numrad)
   real(r8), pointer :: taus(:,:)               !stem transmittance: 1=vis, 2=nir (numrad)
   real(r8), pointer :: z0mr(:)                 !ratio of momentum roughness length to canopy top height (-)
   real(r8), pointer :: displar(:)              !ratio of displacement height to canopy top height (-)
   real(r8), pointer :: roota_par(:)            !CLM rooting distribution parameter [1/m]
   real(r8), pointer :: rootb_par(:)            !CLM rooting distribution parameter [1/m]
   ! new variables for CN code
   real(r8), pointer :: dwood(:)           !wood density (gC/m3)
   real(r8), pointer :: slatop(:)    !specific leaf area at top of canopy, projected area basis [m^2/gC]
   real(r8), pointer :: dsladlai(:)  !dSLA/dLAI, projected area basis [m^2/gC]
   real(r8), pointer :: leafcn(:)    !leaf C:N (gC/gN)
   real(r8), pointer :: flnr(:)      !fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
   real(r8), pointer :: woody(:)     !binary flag for woody lifeform (1=woody, 0=not woody)
   real(r8), pointer :: lflitcn(:)      !leaf litter C:N (gC/gN)
   real(r8), pointer :: frootcn(:)      !fine root C:N (gC/gN)
   real(r8), pointer :: livewdcn(:)     !live wood (phloem and ray parenchyma) C:N (gC/gN)
   real(r8), pointer :: deadwdcn(:)     !dead wood (xylem and heartwood) C:N (gC/gN)
   real(r8), pointer :: graincn(:)      !grain C:N (gC/gN) for prognostic crop model
   real(r8), pointer :: froot_leaf(:)   !allocation parameter: new fine root C per new leaf C (gC/gC)
   real(r8), pointer :: stem_leaf(:)    !allocation parameter: new stem c per new leaf C (gC/gC)
   real(r8), pointer :: croot_stem(:)   !allocation parameter: new coarse root C per new stem C (gC/gC)
   real(r8), pointer :: flivewd(:)      !allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
   real(r8), pointer :: fcur(:)         !allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
   real(r8), pointer :: lf_flab(:)      !leaf litter labile fraction
   real(r8), pointer :: lf_fcel(:)      !leaf litter cellulose fraction
   real(r8), pointer :: lf_flig(:)      !leaf litter lignin fraction
   real(r8), pointer :: fr_flab(:)      !fine root litter labile fraction
   real(r8), pointer :: fr_fcel(:)      !fine root litter cellulose fraction
   real(r8), pointer :: fr_flig(:)      !fine root litter lignin fraction
   real(r8), pointer :: leaf_long(:)    !leaf longevity (yrs)
   real(r8), pointer :: evergreen(:)    !binary flag for evergreen leaf habit (0 or 1)
   real(r8), pointer :: stress_decid(:) !binary flag for stress-deciduous leaf habit (0 or 1)
   real(r8), pointer :: season_decid(:) !binary flag for seasonal-deciduous leaf habit (0 or 1)

!  fire variables added by F. Li and S. Levis
    
  ! combustion completeness factors (0 to 1)
   real(r8), pointer :: cc_leaf(:)        !combustion completeness factor for leaf
   real(r8), pointer :: cc_lstem(:)       !combustion completeness factor for live stem
   real(r8), pointer :: cc_dstem(:)       !combustion completeness factor for dead stem
   real(r8), pointer :: cc_other(:)       !combustion completeness factor for other plant tissues
 !  mortality factors (0 to 1)
   real(r8), pointer :: fm_leaf(:)        !fire-related mortality factor for leaf
   real(r8), pointer :: fm_lstem(:)       !fire-related mortality factor for live stem
   real(r8), pointer :: fm_dstem(:)       !fire-related mortality factor for dead stem
   real(r8), pointer :: fm_other(:)       !fire-related mortality factor for other plant tissues
   real(r8), pointer :: fm_root(:)        !fire-related mortality factor for fine roots 
   real(r8), pointer :: fm_lroot(:)       !fire-related mortality factor for live roots
   real(r8), pointer :: fm_droot(:)       !fire-related mortality factor for dead roots

   real(r8), pointer :: rootprof_beta(:)            !CLM rooting distribution parameter for C and N inputs [unitless]
   ! new variables for crop code
   real(r8), pointer :: fertnitro(:)    ! fertilizer applied
   real(r8), pointer :: fleafcn(:)      ! C:N during grain fill; leaf
   real(r8), pointer :: ffrootcn(:)     ! C:N during grain fill; froot
   real(r8), pointer :: fstemcn(:)      ! C:N during grain fill; stem
end type pft_epc_type

type, public :: decomp_cascade_type
   !-- properties of each pathway along decomposition cascade 
   character(len=8), pointer :: cascade_step_name(:)      ! name of transition
   integer,  pointer :: cascade_donor_pool(:)             ! which pool is C taken from for a given decomposition step
   integer,  pointer :: cascade_receiver_pool(:)          ! which pool is C added to for a given decomposition step
   !-- properties of each decomposing pool
   logical,  pointer :: floating_cn_ratio_decomp_pools(:) ! TRUE => pool has fixed C:N ratio
   character(len=8), pointer :: decomp_pool_name_restart(:)       ! name of pool for restart files
   character(len=8), pointer :: decomp_pool_name_history(:)       ! name of pool for history files
   character(len=20), pointer :: decomp_pool_name_long(:)       ! name of pool for netcdf long names
   character(len=8), pointer :: decomp_pool_name_short(:)       ! name of pool for netcdf short names
   logical, pointer :: is_litter(:)                       ! TRUE => pool is a litter pool
   logical, pointer :: is_soil(:)                         ! TRUE => pool is a soil pool
   logical, pointer :: is_cwd(:)                          ! TRUE => pool is a cwd pool
   real(r8), pointer :: initial_cn_ratio(:)               ! c:n ratio for initialization of pools
   real(r8), pointer :: initial_stock(:)                  ! initial concentration for seeding at spinup
   logical, pointer :: is_metabolic(:)                       ! TRUE => pool is metabolic material
   logical, pointer :: is_cellulose(:)                       ! TRUE => pool is cellulose
   logical, pointer :: is_lignin(:)                          ! TRUE => pool is lignin
   real(r8), pointer :: spinup_factor(:)                  ! factor by which to scale AD and relevant processes by

end type decomp_cascade_type

#if (defined CNDV)
!----------------------------------------------------
! pft DGVM-specific ecophysiological constants structure
!----------------------------------------------------
type, public :: pft_dgvepc_type
   real(r8), pointer :: crownarea_max(:)   !tree maximum crown area [m2]
   real(r8), pointer :: tcmin(:)           !minimum coldest monthly mean temperature [units?]
   real(r8), pointer :: tcmax(:)           !maximum coldest monthly mean temperature [units?]
   real(r8), pointer :: gddmin(:)          !minimum growing degree days (at or above 5 C)
   real(r8), pointer :: twmax(:)           !upper limit of temperature of the warmest month [units?]
   real(r8), pointer :: reinickerp(:)      !parameter in allometric equation
   real(r8), pointer :: allom1(:)          !parameter in allometric
   real(r8), pointer :: allom2(:)          !parameter in allometric
   real(r8), pointer :: allom3(:)          !parameter in allometric
end type pft_dgvepc_type
#endif

!----------------------------------------------------
! pft ecophysiological variables structure
!----------------------------------------------------
type, public :: pft_epv_type
   real(r8), pointer :: dormant_flag(:)         !dormancy flag
   real(r8), pointer :: days_active(:)          !number of days since last dormancy
   real(r8), pointer :: onset_flag(:)           !onset flag
   real(r8), pointer :: onset_counter(:)        !onset days counter
   real(r8), pointer :: onset_gddflag(:)        !onset flag for growing degree day sum
   real(r8), pointer :: onset_fdd(:)            !onset freezing degree days counter
   real(r8), pointer :: onset_gdd(:)            !onset growing degree days
   real(r8), pointer :: onset_swi(:)            !onset soil water index
   real(r8), pointer :: offset_flag(:)          !offset flag
   real(r8), pointer :: offset_counter(:)       !offset days counter
   real(r8), pointer :: offset_fdd(:)           !offset freezing degree days counter
   real(r8), pointer :: offset_swi(:)           !offset soil water index
   real(r8), pointer :: fert_counter(:)         !>0 fertilize; <=0 not
   real(r8), pointer :: grain_flag(:)           !1: grain fill stage; 0: not
   real(r8), pointer :: lgsf(:)                 !long growing season factor [0-1]
   real(r8), pointer :: bglfr(:)                !background litterfall rate (1/s)
   real(r8), pointer :: bgtr(:)                 !background transfer growth rate (1/s)
   real(r8), pointer :: dayl(:)                 !daylength (seconds)
   real(r8), pointer :: prev_dayl(:)            !daylength from previous timestep (seconds)
   real(r8), pointer :: annavg_t2m(:)           !annual average 2m air temperature (K)
   real(r8), pointer :: tempavg_t2m(:)          !temporary average 2m air temperature (K)
   real(r8), pointer :: gpp(:)                  !GPP flux before downregulation (gC/m2/s)
   real(r8), pointer :: availc(:)               !C flux available for allocation (gC/m2/s)
   real(r8), pointer :: xsmrpool_recover(:)     !C flux assigned to recovery of negative cpool (gC/m2/s)
   real(r8), pointer :: xsmrpool_c13ratio(:)    !C13/C(12+13) ratio for xsmrpool (proportion)
   real(r8), pointer :: alloc_pnow(:)           !fraction of current allocation to display as new growth (DIM)
   real(r8), pointer :: c_allometry(:)          !C allocation index (DIM)
   real(r8), pointer :: n_allometry(:)          !N allocation index (DIM)
   real(r8), pointer :: plant_ndemand(:)        !N flux required to support initial GPP (gN/m2/s)
   real(r8), pointer :: tempsum_potential_gpp(:)!temporary annual sum of potential GPP
   real(r8), pointer :: annsum_potential_gpp(:) !annual sum of potential GPP
   real(r8), pointer :: tempmax_retransn(:)     !temporary annual max of retranslocated N pool (gN/m2)
   real(r8), pointer :: annmax_retransn(:)      !annual max of retranslocated N pool (gN/m2)
   real(r8), pointer :: avail_retransn(:)       !N flux available from retranslocation pool (gN/m2/s)
   real(r8), pointer :: plant_nalloc(:)         !total allocated N flux (gN/m2/s)
   real(r8), pointer :: plant_calloc(:)         !total allocated C flux (gC/m2/s)
   real(r8), pointer :: excess_cflux(:)         !C flux not allocated due to downregulation (gC/m2/s)
   real(r8), pointer :: downreg(:)              !fractional reduction in GPP due to N limitation (DIM)
   real(r8), pointer :: prev_leafc_to_litter(:) !previous timestep leaf C litterfall flux (gC/m2/s)
   real(r8), pointer :: prev_frootc_to_litter(:)!previous timestep froot C litterfall flux (gC/m2/s)
   real(r8), pointer :: tempsum_npp(:)          !temporary annual sum of NPP (gC/m2/yr)
   real(r8), pointer :: annsum_npp(:)           !annual sum of NPP (gC/m2/yr)
#if (defined CNDV)
   real(r8), pointer :: tempsum_litfall(:)      !temporary annual sum of litfall (gC/m2/yr)
   real(r8), pointer :: annsum_litfall(:)       !annual sum of litfall (gC/m2/yr)
#endif
   real(r8), pointer :: rc13_canair(:)          !C13O2/C12O2 in canopy air
   real(r8), pointer :: rc13_psnsun(:)          !C13O2/C12O2 in sunlit canopy psn flux
   real(r8), pointer :: rc13_psnsha(:)          !C13O2/C12O2 in shaded canopy psn flux
   real(r8), pointer :: rc14_atm(:)             !C14O2/C12O2 in atmosphere
end type pft_epv_type                        

!----------------------------------------------------
! pft energy state variables structure
!----------------------------------------------------
type, public :: pft_estate_type
   real(r8), pointer :: t_ref2m(:)            !2 m height surface air temperature (Kelvin)
   real(r8), pointer :: t_ref2m_min(:)        !daily minimum of average 2 m height surface air temperature (K)
   real(r8), pointer :: t_ref2m_max(:)        !daily maximum of average 2 m height surface air temperature (K)
   real(r8), pointer :: t_ref2m_min_inst(:)   !instantaneous daily min of average 2 m height surface air temp (K)
   real(r8), pointer :: t_ref2m_max_inst(:)   !instantaneous daily max of average 2 m height surface air temp (K)
   real(r8), pointer :: q_ref2m(:)            !2 m height surface specific humidity (kg/kg)
   real(r8), pointer :: t_ref2m_u(:)          !Urban 2 m height surface air temperature (Kelvin)
   real(r8), pointer :: t_ref2m_r(:)          !Rural 2 m height surface air temperature (Kelvin)
   real(r8), pointer :: t_ref2m_min_u(:)      !Urban daily minimum of average 2 m height surface air temperature (K)
   real(r8), pointer :: t_ref2m_min_r(:)      !Rural daily minimum of average 2 m height surface air temperature (K)
   real(r8), pointer :: t_ref2m_max_u(:)      !Urban daily maximum of average 2 m height surface air temperature (K)
   real(r8), pointer :: t_ref2m_max_r(:)      !Rural daily maximum of average 2 m height surface air temperature (K)
   real(r8), pointer :: t_ref2m_min_inst_u(:) !Urban instantaneous daily min of average 2 m height surface air temp (K)
   real(r8), pointer :: t_ref2m_min_inst_r(:) !Rural instantaneous daily min of average 2 m height surface air temp (K)
  real(r8), pointer :: t_ref2m_max_inst_u(:) !Urban instantaneous daily max of average 2 m height surface air temp (K)
   real(r8), pointer :: t_ref2m_max_inst_r(:) !Rural instantaneous daily max of average 2 m height surface air temp (K)
   real(r8), pointer :: a10tmin(:)            ! 10-day running mean of min 2-m temperature
   real(r8), pointer :: a5tmin(:)             ! 5-day running mean of min 2-m temperature
   real(r8), pointer :: t10(:)                !10-day running mean of the 2 m temperature (K)
   real(r8), pointer :: rh_ref2m(:)           !2 m height surface relative humidity (%)
   real(r8), pointer :: rh_ref2m_u(:)         !Urban 2 m height surface relative humidity (%)
   real(r8), pointer :: rh_ref2m_r(:)         !Rural 2 m height surface relative humidity (%)
   real(r8), pointer :: t_veg(:)              !vegetation temperature (Kelvin)
   real(r8), pointer :: thm(:)                !intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
end type pft_estate_type

!----------------------------------------------------
! pft water state variables structure
!----------------------------------------------------
type, public :: pft_wstate_type
   real(r8), pointer :: h2ocan(:)         !canopy water (mm H2O)
end type pft_wstate_type

!----------------------------------------------------
! pft carbon state variables structure
!----------------------------------------------------
type, public :: pft_cstate_type
   real(r8), pointer :: leafcmax(:)           ! (gC/m2) ann max leaf C
   ! variables for prognostic crop model
   real(r8), pointer :: grainc(:)             ! (gC/m2) grain C
   real(r8), pointer :: grainc_storage(:)     ! (gC/m2) grain C storage
   real(r8), pointer :: grainc_xfer(:)        ! (gC/m2) grain C transfer
   !
   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: xsmrpool(:)           ! (gC/m2) abstract C pool to meet excess MR demand
   real(r8), pointer :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation
   ! summary (diagnostic) state variables, not involved in mass balance
   real(r8), pointer :: dispvegc(:)           ! (gC/m2) displayed veg carbon, excluding storage and cpool
   real(r8), pointer :: storvegc(:)           ! (gC/m2) stored vegetation carbon, excluding cpool
   real(r8), pointer :: totvegc(:)            ! (gC/m2) total vegetation carbon, excluding cpool
   real(r8), pointer :: totpftc(:)            ! (gC/m2) total pft-level carbon, including cpool
#if (defined CN)
   real(r8), pointer :: woodc(:)              ! (gC/m2) wood C
#endif
end type pft_cstate_type

!----------------------------------------------------
! pft nitrogen state variables structure
!----------------------------------------------------
type, public :: pft_nstate_type
   ! variables for prognostic crop model
   real(r8), pointer :: grainn(:)             ! (gN/m2) grain N 
   real(r8), pointer :: grainn_storage(:)     ! (gN/m2) grain N storage
   real(r8), pointer :: grainn_xfer(:)        ! (gN/m2) grain N transfer
   !
   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N 
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
   real(r8), pointer :: npool(:)              ! (gN/m2) temporary plant N pool
   real(r8), pointer :: pft_ntrunc(:)         ! (gN/m2) pft-level sink for N truncation
   ! summary (diagnostic) state variables, not involved in mass balance
   real(r8), pointer :: dispvegn(:)           ! (gN/m2) displayed veg nitrogen, excluding storage
   real(r8), pointer :: storvegn(:)           ! (gN/m2) stored vegetation nitrogen
   real(r8), pointer :: totvegn(:)            ! (gN/m2) total vegetation nitrogen
   real(r8), pointer :: totpftn(:)            ! (gN/m2) total pft-level nitrogen
end type pft_nstate_type

!----------------------------------------------------
! pft VOC state variables structure
!----------------------------------------------------
type, public :: pft_vstate_type
   real(r8), pointer :: t_veg24(:)             ! 24hr average vegetation temperature (K)
   real(r8), pointer :: t_veg240(:)            ! 240hr average vegetation temperature (Kelvin)
   real(r8), pointer :: fsd24(:)               ! 24hr average of direct beam radiation 
   real(r8), pointer :: fsd240(:)              ! 240hr average of direct beam radiation 
   real(r8), pointer :: fsi24(:)               ! 24hr average of diffuse beam radiation 
   real(r8), pointer :: fsi240(:)              ! 240hr average of diffuse beam radiation 
   real(r8), pointer :: fsun24(:)              ! 24hr average of sunlit fraction of canopy 
   real(r8), pointer :: fsun240(:)             ! 240hr average of sunlit fraction of canopy
   real(r8), pointer :: elai_p(:)              ! leaf area index average over timestep 
end type pft_vstate_type

#if (defined CNDV)
!----------------------------------------------------
! pft DGVM state variables structure
!----------------------------------------------------
type, public :: pft_dgvstate_type
   real(r8), pointer :: agddtw(:)              !accumulated growing degree days above twmax
   real(r8), pointer :: agdd(:)                !accumulated growing degree days above 5
   real(r8), pointer :: t_mo(:)                !30-day average temperature (Kelvin)
   real(r8), pointer :: t_mo_min(:)            !annual min of t_mo (Kelvin)
   real(r8), pointer :: prec365(:)             !365-day running mean of tot. precipitation
   logical , pointer :: present(:)             !whether PFT present in patch
   logical , pointer :: pftmayexist(:)         !if .false. then exclude seasonal decid pfts from tropics
   real(r8), pointer :: nind(:)                !number of individuals (#/m**2)
   real(r8), pointer :: lm_ind(:)              !individual leaf mass
   real(r8), pointer :: lai_ind(:)             !LAI per individual
   real(r8), pointer :: fpcinc(:)              !foliar projective cover increment (fraction) 
   real(r8), pointer :: fpcgrid(:)             !foliar projective cover on gridcell (fraction)
   real(r8), pointer :: fpcgridold(:)          !last yr's fpcgrid
   real(r8), pointer :: crownarea(:)           !area that each individual tree takes up (m^2)
   real(r8), pointer :: greffic(:)
   real(r8), pointer :: heatstress(:)
end type pft_dgvstate_type
#endif

!----------------------------------------------------
! pft energy flux variables structure
!----------------------------------------------------
type, public :: pft_eflux_type
   real(r8), pointer :: sabg_soil(:)         !solar radiation absorbed by soil (W/m**2)
   real(r8), pointer :: sabg_snow(:)         !solar radiation absorbed by snow (W/m**2)
   real(r8), pointer :: sabg_chk(:)          !fsno weighted sum (needed by balancecheck, because fsno changes midway)
   real(r8), pointer :: sabg(:)              !solar radiation absorbed by ground (W/m**2)
   real(r8), pointer :: sabv(:)              !solar radiation absorbed by vegetation (W/m**2)
   real(r8), pointer :: fsa(:)               !solar radiation absorbed (total) (W/m**2)
   real(r8), pointer :: fsa_u(:)             !urban solar radiation absorbed (total) (W/m**2)
   real(r8), pointer :: fsa_r(:)             !rural solar radiation absorbed (total) (W/m**2)
   real(r8), pointer :: fsr(:)               !solar radiation reflected (W/m**2)
   real(r8), pointer :: parsun_z(:,:)        !absorbed PAR for sunlit leaves in canopy layer (W/m**2)
   real(r8), pointer :: parsha_z(:,:)        !absorbed PAR for shaded leaves in canopy layer (W/m**2)
   real(r8), pointer :: dlrad(:)             !downward longwave radiation below the canopy [W/m2]
   real(r8), pointer :: ulrad(:)             !upward longwave radiation above the canopy [W/m2]
   real(r8), pointer :: eflx_lh_tot(:)       !total latent heat flux (W/m**2)  [+ to atm]
   real(r8), pointer :: eflx_lh_tot_u(:)     !urban total latent heat flux (W/m**2)  [+ to atm]
   real(r8), pointer :: eflx_lh_tot_r(:)     !rural total latent heat flux (W/m**2)  [+ to atm]
   real(r8), pointer :: eflx_lh_grnd(:)      !ground evaporation heat flux (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_soil_grnd(:)    !soil heat flux (W/m**2) [+ = into soil]
   real(r8), pointer :: eflx_soil_grnd_u(:)  !urban soil heat flux (W/m**2) [+ = into soil]
   real(r8), pointer :: eflx_soil_grnd_r(:)  !rural soil heat flux (W/m**2) [+ = into soil]
   real(r8), pointer :: eflx_sh_tot(:)       !total sensible heat flux (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_sh_tot_u(:)     !urban total sensible heat flux (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_sh_tot_r(:)     !rural total sensible heat flux (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_sh_grnd(:)      !sensible heat flux from ground (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_sh_snow(:)      !sensible heat flux from snow (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_sh_soil(:)      !sensible heat flux from soil (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_sh_h2osfc(:)    !sensible heat flux from surface water (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_sh_veg(:)       !sensible heat flux from leaves (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_lh_vege(:)      !veg evaporation heat flux (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_lh_vegt(:)      !veg transpiration heat flux (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_wasteheat_pft(:) !sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
   real(r8), pointer :: eflx_heat_from_ac_pft(:) !sensible heat flux put back into canyon due to removal by AC (W/m**2)
   real(r8), pointer :: eflx_traffic_pft(:)      !traffic sensible heat flux (W/m**2)
   real(r8), pointer :: eflx_anthro(:)           !total anthropogenic heat flux (W/m**2)
   real(r8), pointer :: cgrnd(:)             !deriv. of soil energy flux wrt to soil temp [w/m2/k]
   real(r8), pointer :: cgrndl(:)            !deriv. of soil latent heat flux wrt soil temp [w/m**2/k]
   real(r8), pointer :: cgrnds(:)            !deriv. of soil sensible heat flux wrt soil temp [w/m2/k]
   real(r8), pointer :: eflx_gnet(:)         !net heat flux into ground (W/m**2)
! New lake field
   real(r8), pointer :: eflx_grnd_lake(:)    !net heat flux into lake / snow surface, excluding light transmission (W/m**2)
   real(r8), pointer :: dgnetdT(:)           !derivative of net ground heat flux wrt soil temp (W/m**2 K)
   real(r8), pointer :: eflx_lwrad_out(:)    !emitted infrared (longwave) radiation (W/m**2)
   real(r8), pointer :: eflx_lwrad_net(:)    !net infrared (longwave) rad (W/m**2) [+ = to atm]
   real(r8), pointer :: eflx_lwrad_net_u(:)  !urban net infrared (longwave) rad (W/m**2) [+ = to atm]
   real(r8), pointer :: eflx_lwrad_net_r(:)  !rural net infrared (longwave) rad (W/m**2) [+ = to atm]
   real(r8), pointer :: netrad(:)            !net radiation (W/m**2) [+ = to sfc]
   real(r8), pointer :: fsds_vis_d(:)        !incident direct beam vis solar radiation (W/m**2)
   real(r8), pointer :: fsds_nir_d(:)        !incident direct beam nir solar radiation (W/m**2)
   real(r8), pointer :: fsds_vis_i(:)        !incident diffuse vis solar radiation (W/m**2)
   real(r8), pointer :: fsds_nir_i(:)        !incident diffuse nir solar radiation (W/m**2)
   real(r8), pointer :: fsr_vis_d(:)         !reflected direct beam vis solar radiation (W/m**2)
   real(r8), pointer :: fsr_nir_d(:)         !reflected direct beam nir solar radiation (W/m**2)
   real(r8), pointer :: fsr_vis_i(:)         !reflected diffuse vis solar radiation (W/m**2)
   real(r8), pointer :: fsr_nir_i(:)         !reflected diffuse nir solar radiation (W/m**2)
   real(r8), pointer :: fsds_vis_d_ln(:)     !incident direct beam vis solar radiation at local noon (W/m**2)
   real(r8), pointer :: fsds_vis_i_ln(:)     !incident diffuse beam vis solar radiation at local noon (W/m**2)
   real(r8), pointer :: parveg_ln(:)         !absorbed par by vegetation at local noon (W/m**2)
   real(r8), pointer :: fsds_nir_d_ln(:)     !incident direct beam nir solar radiation at local noon (W/m**2)
   real(r8), pointer :: fsr_vis_d_ln(:)      !reflected direct beam vis solar radiation at local noon (W/m**2)
   real(r8), pointer :: fsr_nir_d_ln(:)      !reflected direct beam nir solar radiation at local noon (W/m**2)
   real(r8), pointer :: sabg_lyr(:,:)     ! absorbed radiation in each snow layer and top soil layer (pft,lyr) [W/m2]
   real(r8), pointer :: sabg_pen(:)       ! (rural) shortwave radiation penetrating top soisno layer [W/m2]
   real(r8), pointer :: sfc_frc_aer(:)    ! surface forcing of snow with all aerosols (pft) [W/m2]
   real(r8), pointer :: sfc_frc_bc(:)     ! surface forcing of snow with BC (pft) [W/m2]
   real(r8), pointer :: sfc_frc_oc(:)     ! surface forcing of snow with OC (pft) [W/m2]
   real(r8), pointer :: sfc_frc_dst(:)    ! surface forcing of snow with dust (pft) [W/m2]
   real(r8), pointer :: sfc_frc_aer_sno(:)! surface forcing of snow with all aerosols, averaged only when snow is present (pft) [W/m2]
   real(r8), pointer :: sfc_frc_bc_sno(:) ! surface forcing of snow with BC, averaged only when snow is present (pft) [W/m2]
   real(r8), pointer :: sfc_frc_oc_sno(:) ! surface forcing of snow with OC, averaged only when snow is present (pft) [W/m2]
   real(r8), pointer :: sfc_frc_dst_sno(:)! surface forcing of snow with dust, averaged only when snow is present (pft) [W/m2]
   real(r8), pointer :: fsr_sno_vd(:)     ! reflected direct beam vis solar radiation from snow (W/m**2)
   real(r8), pointer :: fsr_sno_nd(:)     ! reflected direct beam NIR solar radiation from snow (W/m**2)
   real(r8), pointer :: fsr_sno_vi(:)     ! reflected diffuse vis solar radiation from snow (W/m**2)
   real(r8), pointer :: fsr_sno_ni(:)     ! reflected diffuse NIR solar radiation from snow (W/m**2)
   real(r8), pointer :: fsds_sno_vd(:)    ! incident visible, direct radiation on snow (for history files)  [W/m2]
   real(r8), pointer :: fsds_sno_nd(:)    ! incident near-IR, direct radiation on snow (for history files)  [W/m2]
   real(r8), pointer :: fsds_sno_vi(:)    ! incident visible, diffuse radiation on snow (for history files) [W/m2]
   real(r8), pointer :: fsds_sno_ni(:)    ! incident near-IR, diffuse radiation on snow (for history files) [W/m2]
end type pft_eflux_type

!----------------------------------------------------
! pft momentum flux variables structure
!----------------------------------------------------
type, public :: pft_mflux_type
   real(r8),pointer ::  taux(:)           !wind (shear) stress: e-w (kg/m/s**2)
   real(r8),pointer ::  tauy(:)           !wind (shear) stress: n-s (kg/m/s**2)
end type pft_mflux_type

!----------------------------------------------------
! pft water flux variables structure
!----------------------------------------------------
type, public :: pft_wflux_type
   real(r8), pointer :: qflx_prec_intr(:) !interception of precipitation [mm/s]
   real(r8), pointer :: qflx_prec_grnd(:) !water onto ground including canopy runoff [kg/(m2 s)]
   real(r8), pointer :: qflx_rain_grnd(:) !rain on ground after interception (mm H2O/s) [+]
   real(r8), pointer :: qflx_snow_grnd(:) !snow on ground after interception (mm H2O/s) [+]
   real(r8), pointer :: qflx_snwcp_ice(:) !excess snowfall due to snow capping (mm H2O /s) [+]
   real(r8), pointer :: qflx_snwcp_liq(:) !excess rainfall due to snow capping (mm H2O /s) [+]
   real(r8), pointer :: qflx_evap_veg(:)  !vegetation evaporation (mm H2O/s) (+ = to atm)
   real(r8), pointer :: qflx_tran_veg(:)  !vegetation transpiration (mm H2O/s) (+ = to atm)
   real(r8), pointer :: qflx_evap_can(:)  !evaporation from leaves and stems 
   real(r8), pointer :: qflx_evap_soi(:)  !soil evaporation (mm H2O/s) (+ = to atm)
   real(r8), pointer :: qflx_evap_tot(:)  !qflx_evap_soi + qflx_evap_can + qflx_tran_veg
   real(r8), pointer :: qflx_evap_grnd(:) !ground surface evaporation rate (mm H2O/s) [+]
   real(r8), pointer :: qflx_dew_grnd(:)  !ground surface dew formation (mm H2O /s) [+]
   real(r8), pointer :: qflx_sub_snow(:)  !sublimation rate from snow pack (mm H2O /s) [+]
   real(r8), pointer :: qflx_dew_snow(:)  !surface dew added to snow pack (mm H2O /s) [+]
   real(r8), pointer :: qflx_ev_snow(:)   !snow evaporation (mm H2O/s) (+ = to atm)
   real(r8), pointer :: qflx_ev_soil(:)   !soil evaporation (mm H2O/s) (+ = to atm)
   real(r8), pointer :: qflx_ev_h2osfc(:) !h2osfc evaporation (mm H2O/s) (+ = to atm)
end type pft_wflux_type

!----------------------------------------------------
! pft carbon flux variables structure
!----------------------------------------------------
type, public :: pft_cflux_type
   real(r8), pointer :: psnsun(:)         !sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: psnsha(:)         !shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: psnsun_z(:,:)     !canopy layer: sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: psnsha_z(:,:)     !canopy layer: shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: cisun_z(:,:)      !intracellular sunlit leaf CO2 (Pa)
   real(r8), pointer :: cisha_z(:,:)      !intracellular shaded leaf CO2 (Pa)
   real(r8), pointer :: lmrsun(:)         !sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
   real(r8), pointer :: lmrsha(:)         !shaded leaf maintenance respiration rate (umol CO2/m**2/s)
   real(r8), pointer :: lmrsun_z(:,:)     !canopy layer: sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
   real(r8), pointer :: lmrsha_z(:,:)     !canopy layer: shaded leaf maintenance respiration rate (umol CO2/m**2/s)
   real(r8), pointer :: fpsn(:)           !photosynthesis (umol CO2 /m**2 /s)
   real(r8), pointer :: fco2(:)           !net CO2 flux (umol CO2 /m**2 /s) [+ = to atm]
   real(r8), pointer :: psnsun_wc(:)      !Rubsico-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: psnsha_wc(:)      !Rubsico-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: fpsn_wc(:)        !Rubisco-limited photosynthesis (umol CO2 /m**2 /s)
   real(r8), pointer :: psnsun_wj(:)      !RuBP-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: psnsha_wj(:)      !RuBP-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: fpsn_wj(:)        !RuBP-limited photosynthesis (umol CO2 /m**2 /s)
   real(r8), pointer :: psnsun_wp(:)      !product-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: psnsha_wp(:)      !product-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: fpsn_wp(:)        !product-limited photosynthesis (umol CO2 /m**2 /s)
   ! new variables for CN code
   ! gap mortality fluxes
   real(r8), pointer :: m_leafc_to_litter(:)                 ! leaf C mortality (gC/m2/s)
   real(r8), pointer :: m_leafc_storage_to_litter(:)         ! leaf C storage mortality (gC/m2/s)
   real(r8), pointer :: m_leafc_xfer_to_litter(:)            ! leaf C transfer mortality (gC/m2/s)
   real(r8), pointer :: m_frootc_to_litter(:)                ! fine root C mortality (gC/m2/s)
   real(r8), pointer :: m_frootc_storage_to_litter(:)        ! fine root C storage mortality (gC/m2/s)
   real(r8), pointer :: m_frootc_xfer_to_litter(:)           ! fine root C transfer mortality (gC/m2/s)
   real(r8), pointer :: m_livestemc_to_litter(:)             ! live stem C mortality (gC/m2/s)
   real(r8), pointer :: m_livestemc_storage_to_litter(:)     ! live stem C storage mortality (gC/m2/s)
   real(r8), pointer :: m_livestemc_xfer_to_litter(:)        ! live stem C transfer mortality (gC/m2/s)
   real(r8), pointer :: m_deadstemc_to_litter(:)             ! dead stem C mortality (gC/m2/s)
   real(r8), pointer :: m_deadstemc_storage_to_litter(:)     ! dead stem C storage mortality (gC/m2/s)
   real(r8), pointer :: m_deadstemc_xfer_to_litter(:)        ! dead stem C transfer mortality (gC/m2/s)
   real(r8), pointer :: m_livecrootc_to_litter(:)            ! live coarse root C mortality (gC/m2/s)
   real(r8), pointer :: m_livecrootc_storage_to_litter(:)    ! live coarse root C storage mortality (gC/m2/s)
   real(r8), pointer :: m_livecrootc_xfer_to_litter(:)       ! live coarse root C transfer mortality (gC/m2/s)
   real(r8), pointer :: m_deadcrootc_to_litter(:)            ! dead coarse root C mortality (gC/m2/s)
   real(r8), pointer :: m_deadcrootc_storage_to_litter(:)    ! dead coarse root C storage mortality (gC/m2/s)
   real(r8), pointer :: m_deadcrootc_xfer_to_litter(:)       ! dead coarse root C transfer mortality (gC/m2/s)
   real(r8), pointer :: m_gresp_storage_to_litter(:)         ! growth respiration storage mortality (gC/m2/s)
   real(r8), pointer :: m_gresp_xfer_to_litter(:)            ! growth respiration transfer mortality (gC/m2/s)
   ! harvest mortality fluxes
   real(r8), pointer :: hrv_leafc_to_litter(:)               ! leaf C harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_leafc_storage_to_litter(:)       ! leaf C storage harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_leafc_xfer_to_litter(:)          ! leaf C transfer harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_frootc_to_litter(:)              ! fine root C harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_frootc_storage_to_litter(:)      ! fine root C storage harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_frootc_xfer_to_litter(:)         ! fine root C transfer harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_livestemc_to_litter(:)           ! live stem C harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_livestemc_storage_to_litter(:)   ! live stem C storage harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_livestemc_xfer_to_litter(:)      ! live stem C transfer harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_deadstemc_to_prod10c(:)          ! dead stem C harvest to 10-year product pool (gC/m2/s)
   real(r8), pointer :: hrv_deadstemc_to_prod100c(:)         ! dead stem C harvest to 100-year product pool (gC/m2/s)
   real(r8), pointer :: hrv_deadstemc_storage_to_litter(:)   ! dead stem C storage harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_deadstemc_xfer_to_litter(:)      ! dead stem C transfer harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_livecrootc_to_litter(:)          ! live coarse root C harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_livecrootc_storage_to_litter(:)  ! live coarse root C storage harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_livecrootc_xfer_to_litter(:)     ! live coarse root C transfer harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_deadcrootc_to_litter(:)          ! dead coarse root C harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_deadcrootc_storage_to_litter(:)  ! dead coarse root C storage harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_deadcrootc_xfer_to_litter(:)     ! dead coarse root C transfer harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_gresp_storage_to_litter(:)       ! growth respiration storage harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_gresp_xfer_to_litter(:)          ! growth respiration transfer harvest mortality (gC/m2/s)
   real(r8), pointer :: hrv_xsmrpool_to_atm(:)               ! excess MR pool harvest mortality (gC/m2/s)

   ! PFT-level fire C fluxes added by F. Li and S. Levis
   real(r8), pointer :: m_leafc_to_fire(:)             ! (gC/m2/s) fire C emissions from leafc 
   real(r8), pointer :: m_leafc_storage_to_fire(:)     ! (gC/m2/s) fire C emissions from leafc_storage             
   real(r8), pointer :: m_leafc_xfer_to_fire(:)        ! (gC/m2/s) fire C emissions from leafc_xfer
   real(r8), pointer :: m_livestemc_to_fire(:)         ! (gC/m2/s) fire C emissions from livestemc
   real(r8), pointer :: m_livestemc_storage_to_fire(:) ! (gC/m2/s) fire C emissions from livestemc_storage       
   real(r8), pointer :: m_livestemc_xfer_to_fire(:)    ! (gC/m2/s) fire C emissions from livestemc_xfer
   real(r8), pointer :: m_deadstemc_to_fire(:)         ! (gC/m2/s) fire C emissions from deadstemc_xfer
   real(r8), pointer :: m_deadstemc_storage_to_fire(:) ! (gC/m2/s) fire C emissions from deadstemc_storage         
   real(r8), pointer :: m_deadstemc_xfer_to_fire(:)    ! (gC/m2/s) fire C emissions from deadstemc_xfer
   real(r8), pointer :: m_frootc_to_fire(:)            ! (gC/m2/s) fire C emissions from frootc
   real(r8), pointer :: m_frootc_storage_to_fire(:)    ! (gC/m2/s) fire C emissions from frootc_storage
   real(r8), pointer :: m_frootc_xfer_to_fire(:)       ! (gC/m2/s) fire C emissions from frootc_xfer
   real(r8), pointer :: m_livecrootc_to_fire(:)        ! (gC/m2/s) fire C emissions from livecrootc
   real(r8), pointer :: m_livecrootc_storage_to_fire(:)! (gC/m2/s) fire C emissions from livecrootc_storage     
   real(r8), pointer :: m_livecrootc_xfer_to_fire(:)   ! (gC/m2/s) fire C emissions from livecrootc_xfer
   real(r8), pointer :: m_deadcrootc_to_fire(:)        ! (gC/m2/s) fire C emissions from deadcrootc
   real(r8), pointer :: m_deadcrootc_storage_to_fire(:)! (gC/m2/s) fire C emissions from deadcrootc_storage 
   real(r8), pointer :: m_deadcrootc_xfer_to_fire(:)   ! (gC/m2/s) fire C emissions from deadcrootc_xfer
   real(r8), pointer :: m_gresp_storage_to_fire(:)     ! (gC/m2/s) fire C emissions from gresp_storage 
   real(r8), pointer :: m_gresp_xfer_to_fire(:)        ! (gC/m2/s) fire C emissions from gresp_xfer
   real(r8), pointer :: m_leafc_to_litter_fire(:)      ! (gC/m2/s) from leafc to litter c due to fire
   real(r8), pointer :: m_leafc_storage_to_litter_fire(:)  ! (gC/m2/s) from leafc_storage to litter C  due to fire               
   real(r8), pointer :: m_leafc_xfer_to_litter_fire(:)     ! (gC/m2/s) from leafc_xfer to litter C  due to fire               
   real(r8), pointer :: m_livestemc_to_litter_fire(:)      ! (gC/m2/s) from livestemc to litter C  due to fire               
   real(r8), pointer :: m_livestemc_storage_to_litter_fire(:)  ! (gC/m2/s) from livestemc_storage to litter C due to fire      
   real(r8), pointer :: m_livestemc_xfer_to_litter_fire(:) !(gC/m2/s) from livestemc_xfer to litter C due to fire      
   real(r8), pointer :: m_livestemc_to_deadstemc_fire(:)   !(gC/m2/s) from livestemc to deadstemc due to fire       
   real(r8), pointer :: m_deadstemc_to_litter_fire(:)      !(gC/m2/s) from deadstemc to litter C due to fire      
   real(r8), pointer :: m_deadstemc_storage_to_litter_fire(:)  !(gC/m2/s) from deadstemc_storage to litter C due to fire               
   real(r8), pointer :: m_deadstemc_xfer_to_litter_fire(:)     !(gC/m2/s) from deadstemc_xfer to litter C due to fire               
   real(r8), pointer :: m_frootc_to_litter_fire(:)             !(gC/m2/s) from frootc to litter C due to fire               
   real(r8), pointer :: m_frootc_storage_to_litter_fire(:)     !(gC/m2/s) from frootc_storage to litter C due to fire               
   real(r8), pointer :: m_frootc_xfer_to_litter_fire(:)        !(gC/m2/s) from frootc_xfer to litter C due to fire               
   real(r8), pointer :: m_livecrootc_to_litter_fire(:)         !(gC/m2/s) from livecrootc to litter C due to fire                     
   real(r8), pointer :: m_livecrootc_storage_to_litter_fire(:) !(gC/m2/s) from livecrootc_storage to litter C due to fire                     
   real(r8), pointer :: m_livecrootc_xfer_to_litter_fire(:)    !(gC/m2/s) from livecrootc_xfer to litter C due to fire                     
   real(r8), pointer :: m_livecrootc_to_deadcrootc_fire(:)     !(gC/m2/s) from livecrootc to deadstemc due to fire        
   real(r8), pointer :: m_deadcrootc_to_litter_fire(:)         !(gC/m2/s) from deadcrootc to litter C due to fire                       
   real(r8), pointer :: m_deadcrootc_storage_to_litter_fire(:) !(gC/m2/s) from deadcrootc_storage to litter C due to fire                       
   real(r8), pointer :: m_deadcrootc_xfer_to_litter_fire(:)    !(gC/m2/s) from deadcrootc_xfer to litter C due to fire                       
   real(r8), pointer :: m_gresp_storage_to_litter_fire(:)      !(gC/m2/s) from gresp_storage to litter C due to fire                       
   real(r8), pointer :: m_gresp_xfer_to_litter_fire(:)          !(gC/m2/s) from gresp_xfer to litter C due to fire                       
 
   ! phenology fluxes from transfer pools                     
   real(r8), pointer :: grainc_xfer_to_grainc(:)             ! grain C growth from storage for prognostic crop(gC/m2/s)
   real(r8), pointer :: leafc_xfer_to_leafc(:)               ! leaf C growth from storage (gC/m2/s)
   real(r8), pointer :: frootc_xfer_to_frootc(:)             ! fine root C growth from storage (gC/m2/s)
   real(r8), pointer :: livestemc_xfer_to_livestemc(:)       ! live stem C growth from storage (gC/m2/s)
   real(r8), pointer :: deadstemc_xfer_to_deadstemc(:)       ! dead stem C growth from storage (gC/m2/s)
   real(r8), pointer :: livecrootc_xfer_to_livecrootc(:)     ! live coarse root C growth from storage (gC/m2/s)
   real(r8), pointer :: deadcrootc_xfer_to_deadcrootc(:)     ! dead coarse root C growth from storage (gC/m2/s)
   ! leaf and fine root litterfall                           
   real(r8), pointer :: leafc_to_litter(:)                   ! leaf C litterfall (gC/m2/s)
   real(r8), pointer :: frootc_to_litter(:)                  ! fine root C litterfall (gC/m2/s)
   real(r8), pointer :: livestemc_to_litter(:)               ! live stem C litterfall (gC/m2/s)
   real(r8), pointer :: grainc_to_food(:)                    ! grain C to food for prognostic crop(gC/m2/s)
   ! maintenance respiration fluxes                          
   real(r8), pointer :: leaf_mr(:)                           ! leaf maintenance respiration (gC/m2/s)
   real(r8), pointer :: froot_mr(:)                          ! fine root maintenance respiration (gC/m2/s)
   real(r8), pointer :: livestem_mr(:)                       ! live stem maintenance respiration (gC/m2/s)
   real(r8), pointer :: livecroot_mr(:)                      ! live coarse root maintenance respiration (gC/m2/s)
   real(r8), pointer :: grain_mr(:)                          ! crop grain or organs maint. respiration (gC/m2/s)
   real(r8), pointer :: leaf_curmr(:)                        ! leaf maintenance respiration from current GPP (gC/m2/s)
   real(r8), pointer :: froot_curmr(:)                       ! fine root maintenance respiration from current GPP (gC/m2/s)
   real(r8), pointer :: livestem_curmr(:)                    ! live stem maintenance respiration from current GPP (gC/m2/s)
   real(r8), pointer :: livecroot_curmr(:)                   ! live coarse root maintenance respiration from current GPP (gC/m2/s)
   real(r8), pointer :: grain_curmr(:)                       ! crop grain or organs maint. respiration from current GPP (gC/m2/s)
   real(r8), pointer :: leaf_xsmr(:)                         ! leaf maintenance respiration from storage (gC/m2/s)
   real(r8), pointer :: froot_xsmr(:)                        ! fine root maintenance respiration from storage (gC/m2/s)
   real(r8), pointer :: livestem_xsmr(:)                     ! live stem maintenance respiration from storage (gC/m2/s)
   real(r8), pointer :: livecroot_xsmr(:)                    ! live coarse root maintenance respiration from storage (gC/m2/s)
   real(r8), pointer :: grain_xsmr(:)                        ! crop grain or organs maint. respiration from storage (gC/m2/s)
   ! photosynthesis fluxes                                   
   real(r8), pointer :: psnsun_to_cpool(:)                   ! C fixation from sunlit canopy (gC/m2/s)
   real(r8), pointer :: psnshade_to_cpool(:)                 ! C fixation from shaded canopy (gC/m2/s)
   ! allocation fluxes, from current GPP                     
   real(r8), pointer :: cpool_to_xsmrpool(:)                 ! allocation to maintenance respiration storage pool (gC/m2/s)
   real(r8), pointer :: cpool_to_grainc(:)                   ! allocation to grain C for prognostic crop(gC/m2/s)
   real(r8), pointer :: cpool_to_grainc_storage(:)           ! allocation to grain C storage for prognostic crop(gC/m2/s)
   real(r8), pointer :: cpool_to_leafc(:)                    ! allocation to leaf C (gC/m2/s)
   real(r8), pointer :: cpool_to_leafc_storage(:)            ! allocation to leaf C storage (gC/m2/s)
   real(r8), pointer :: cpool_to_frootc(:)                   ! allocation to fine root C (gC/m2/s)
   real(r8), pointer :: cpool_to_frootc_storage(:)           ! allocation to fine root C storage (gC/m2/s)
   real(r8), pointer :: cpool_to_livestemc(:)                ! allocation to live stem C (gC/m2/s)
   real(r8), pointer :: cpool_to_livestemc_storage(:)        ! allocation to live stem C storage (gC/m2/s)
   real(r8), pointer :: cpool_to_deadstemc(:)                ! allocation to dead stem C (gC/m2/s)
   real(r8), pointer :: cpool_to_deadstemc_storage(:)        ! allocation to dead stem C storage (gC/m2/s)
   real(r8), pointer :: cpool_to_livecrootc(:)               ! allocation to live coarse root C (gC/m2/s)
   real(r8), pointer :: cpool_to_livecrootc_storage(:)       ! allocation to live coarse root C storage (gC/m2/s)
   real(r8), pointer :: cpool_to_deadcrootc(:)               ! allocation to dead coarse root C (gC/m2/s)
   real(r8), pointer :: cpool_to_deadcrootc_storage(:)       ! allocation to dead coarse root C storage (gC/m2/s)
   real(r8), pointer :: cpool_to_gresp_storage(:)            ! allocation to growth respiration storage (gC/m2/s)
   ! growth respiration fluxes                               
   real(r8), pointer :: xsmrpool_to_atm(:)                   ! excess MR pool harvest mortality (gC/m2/s)
   real(r8), pointer :: cpool_leaf_gr(:)                     ! leaf growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_leaf_storage_gr(:)             ! leaf growth respiration to storage (gC/m2/s)
   real(r8), pointer :: transfer_leaf_gr(:)                  ! leaf growth respiration from storage (gC/m2/s)
   real(r8), pointer :: cpool_froot_gr(:)                    ! fine root growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_froot_storage_gr(:)            ! fine root  growth respiration to storage (gC/m2/s)
   real(r8), pointer :: transfer_froot_gr(:)                 ! fine root  growth respiration from storage (gC/m2/s)
   real(r8), pointer :: cpool_livestem_gr(:)                 ! live stem growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_livestem_storage_gr(:)         ! live stem growth respiration to storage (gC/m2/s)
   real(r8), pointer :: transfer_livestem_gr(:)              ! live stem growth respiration from storage (gC/m2/s)
   real(r8), pointer :: cpool_deadstem_gr(:)                 ! dead stem growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_deadstem_storage_gr(:)         ! dead stem growth respiration to storage (gC/m2/s)
   real(r8), pointer :: transfer_deadstem_gr(:)              ! dead stem growth respiration from storage (gC/m2/s)
   real(r8), pointer :: cpool_livecroot_gr(:)                ! live coarse root growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_livecroot_storage_gr(:)        ! live coarse root growth respiration to storage (gC/m2/s)
   real(r8), pointer :: transfer_livecroot_gr(:)             ! live coarse root growth respiration from storage (gC/m2/s)
   real(r8), pointer :: cpool_deadcroot_gr(:)                ! dead coarse root growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_deadcroot_storage_gr(:)        ! dead coarse root growth respiration to storage (gC/m2/s)
   real(r8), pointer :: transfer_deadcroot_gr(:)             ! dead coarse root growth respiration from storage (gC/m2/s)
   ! growth respiration for prognostic crop model
   real(r8), pointer :: cpool_grain_gr(:)                    ! grain growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_grain_storage_gr(:)            ! grain growth respiration to storage (gC/m2/s)
   real(r8), pointer :: transfer_grain_gr(:)                 ! grain growth respiration from storage (gC/m2/s)
   ! annual turnover of storage to transfer pools            
   real(r8), pointer :: grainc_storage_to_xfer(:)            ! grain C shift storage to transfer for prognostic crop model (gC/m2/s)
   real(r8), pointer :: leafc_storage_to_xfer(:)             ! leaf C shift storage to transfer (gC/m2/s)
   real(r8), pointer :: frootc_storage_to_xfer(:)            ! fine root C shift storage to transfer (gC/m2/s)
   real(r8), pointer :: livestemc_storage_to_xfer(:)         ! live stem C shift storage to transfer (gC/m2/s)
   real(r8), pointer :: deadstemc_storage_to_xfer(:)         ! dead stem C shift storage to transfer (gC/m2/s)
   real(r8), pointer :: livecrootc_storage_to_xfer(:)        ! live coarse root C shift storage to transfer (gC/m2/s)
   real(r8), pointer :: deadcrootc_storage_to_xfer(:)        ! dead coarse root C shift storage to transfer (gC/m2/s)
   real(r8), pointer :: gresp_storage_to_xfer(:)             ! growth respiration shift storage to transfer (gC/m2/s)
   ! turnover of livewood to deadwood
   real(r8), pointer :: livestemc_to_deadstemc(:)            ! live stem C turnover (gC/m2/s)
   real(r8), pointer :: livecrootc_to_deadcrootc(:)          ! live coarse root C turnover (gC/m2/s)
   ! summary (diagnostic) flux variables, not involved in mass balance
   real(r8), pointer :: gpp(:)            ! (gC/m2/s) gross primary production 
   real(r8), pointer :: mr(:)             ! (gC/m2/s) maintenance respiration
   real(r8), pointer :: current_gr(:)     ! (gC/m2/s) growth resp for new growth displayed in this timestep
   real(r8), pointer :: transfer_gr(:)    ! (gC/m2/s) growth resp for transfer growth displayed in this timestep
   real(r8), pointer :: storage_gr(:)     ! (gC/m2/s) growth resp for growth sent to storage for later display
   real(r8), pointer :: gr(:)             ! (gC/m2/s) total growth respiration
   real(r8), pointer :: ar(:)             ! (gC/m2/s) autotrophic respiration (MR + GR)
   real(r8), pointer :: rr(:)             ! (gC/m2/s) root respiration (fine root MR + total root GR)
   real(r8), pointer :: npp(:)            ! (gC/m2/s) net primary production
   real(r8), pointer :: agnpp(:)          ! (gC/m2/s) aboveground NPP
   real(r8), pointer :: bgnpp(:)          ! (gC/m2/s) belowground NPP
   real(r8), pointer :: litfall(:)        ! (gC/m2/s) litterfall (leaves and fine roots)
   real(r8), pointer :: vegfire(:)        ! (gC/m2/s) pft-level fire loss (obsolete, mark for removal)
   real(r8), pointer :: wood_harvestc(:)  ! (gC/m2/s) pft-level wood harvest (to product pools)
   real(r8), pointer :: pft_cinputs(:)    ! (gC/m2/s) pft-level carbon inputs (for balance checking)
   real(r8), pointer :: pft_coutputs(:)   ! (gC/m2/s) pft-level carbon outputs (for balance checking)
#if (defined CN)
   ! CLAMP summary (diagnostic) variables, not involved in mass balance
   real(r8), pointer :: frootc_alloc(:)   ! (gC/m2/s) pft-level fine root C alloc
   real(r8), pointer :: frootc_loss(:)    ! (gC/m2/s) pft-level fine root C loss
   real(r8), pointer :: leafc_alloc(:)    ! (gC/m2/s) pft-level leaf C alloc
   real(r8), pointer :: leafc_loss(:)     ! (gC/m2/s) pft-level leaf C loss
   real(r8), pointer :: woodc_alloc(:)    ! (gC/m2/s) pft-level wood C alloc
   real(r8), pointer :: woodc_loss(:)     ! (gC/m2/s) pft-level wood C loss
#endif
   ! new variables for fire code
   real(r8), pointer :: pft_fire_closs(:) ! (gC/m2/s) total pft-level fire C loss 
! For aerenchyma calculations in CH4 code
#if (defined LCH4)
   real(r8), pointer :: annavg_agnpp(:)          ! (gC/m2/s) annual average aboveground NPP
   real(r8), pointer :: annavg_bgnpp(:)          ! (gC/m2/s) annual average belowground NPP
   real(r8), pointer :: tempavg_agnpp(:)          ! (gC/m2/s) temp. average aboveground NPP
   real(r8), pointer :: tempavg_bgnpp(:)          ! (gC/m2/s) temp. average belowground NPP
#endif
end type pft_cflux_type

!----------------------------------------------------
! pft nitrogen flux variables structure
!----------------------------------------------------
type, public :: pft_nflux_type
   ! new variables for CN code
   ! gap mortality fluxes
   real(r8), pointer :: m_leafn_to_litter(:)                ! leaf N mortality (gN/m2/s)
   real(r8), pointer :: m_frootn_to_litter(:)               ! fine root N mortality (gN/m2/s)
   real(r8), pointer :: m_leafn_storage_to_litter(:)        ! leaf N storage mortality (gN/m2/s)
   real(r8), pointer :: m_frootn_storage_to_litter(:)       ! fine root N storage mortality (gN/m2/s)
   real(r8), pointer :: m_livestemn_storage_to_litter(:)    ! live stem N storage mortality (gN/m2/s)
   real(r8), pointer :: m_deadstemn_storage_to_litter(:)    ! dead stem N storage mortality (gN/m2/s)
   real(r8), pointer :: m_livecrootn_storage_to_litter(:)   ! live coarse root N storage mortality (gN/m2/s)
   real(r8), pointer :: m_deadcrootn_storage_to_litter(:)   ! dead coarse root N storage mortality (gN/m2/s)
   real(r8), pointer :: m_leafn_xfer_to_litter(:)           ! leaf N transfer mortality (gN/m2/s)
   real(r8), pointer :: m_frootn_xfer_to_litter(:)          ! fine root N transfer mortality (gN/m2/s)
   real(r8), pointer :: m_livestemn_xfer_to_litter(:)       ! live stem N transfer mortality (gN/m2/s)
   real(r8), pointer :: m_deadstemn_xfer_to_litter(:)       ! dead stem N transfer mortality (gN/m2/s)
   real(r8), pointer :: m_livecrootn_xfer_to_litter(:)      ! live coarse root N transfer mortality (gN/m2/s)
   real(r8), pointer :: m_deadcrootn_xfer_to_litter(:)      ! dead coarse root N transfer mortality (gN/m2/s)
   real(r8), pointer :: m_livestemn_to_litter(:)            ! live stem N mortality (gN/m2/s)
   real(r8), pointer :: m_deadstemn_to_litter(:)            ! dead stem N mortality (gN/m2/s)
   real(r8), pointer :: m_livecrootn_to_litter(:)           ! live coarse root N mortality (gN/m2/s)
   real(r8), pointer :: m_deadcrootn_to_litter(:)           ! dead coarse root N mortality (gN/m2/s)
   real(r8), pointer :: m_retransn_to_litter(:)             ! retranslocated N pool mortality (gN/m2/s)
   ! harvest mortality fluxes
   real(r8), pointer :: hrv_leafn_to_litter(:)                ! leaf N harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_frootn_to_litter(:)               ! fine root N harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_leafn_storage_to_litter(:)        ! leaf N storage harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_frootn_storage_to_litter(:)       ! fine root N storage harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_livestemn_storage_to_litter(:)    ! live stem N storage harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_deadstemn_storage_to_litter(:)    ! dead stem N storage harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_livecrootn_storage_to_litter(:)   ! live coarse root N storage harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_deadcrootn_storage_to_litter(:)   ! dead coarse root N storage harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_leafn_xfer_to_litter(:)           ! leaf N transfer harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_frootn_xfer_to_litter(:)          ! fine root N transfer harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_livestemn_xfer_to_litter(:)       ! live stem N transfer harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_deadstemn_xfer_to_litter(:)       ! dead stem N transfer harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_livecrootn_xfer_to_litter(:)      ! live coarse root N transfer harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_deadcrootn_xfer_to_litter(:)      ! dead coarse root N transfer harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_livestemn_to_litter(:)            ! live stem N harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_deadstemn_to_prod10n(:)           ! dead stem N harvest to 10-year product pool (gN/m2/s)
   real(r8), pointer :: hrv_deadstemn_to_prod100n(:)          ! dead stem N harvest to 100-year product pool (gN/m2/s)
   real(r8), pointer :: hrv_livecrootn_to_litter(:)           ! live coarse root N harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_deadcrootn_to_litter(:)           ! dead coarse root N harvest mortality (gN/m2/s)
   real(r8), pointer :: hrv_retransn_to_litter(:)             ! retranslocated N pool harvest mortality (gN/m2/s)

  ! PFT-level fire N fluxes added by F. Li and S. Levis
    real(r8), pointer :: m_leafn_to_fire(:)             ! (gN/m2/s) fire N emissions from leafn 
   real(r8), pointer :: m_leafn_storage_to_fire(:)     ! (gN/m2/s) fire N emissions from leafn_storage            
   real(r8), pointer :: m_leafn_xfer_to_fire(:)        ! (gN/m2/s) fire N emissions from leafn_xfer     
   real(r8), pointer :: m_livestemn_to_fire(:)         ! (gN/m2/s) fire N emissions from livestemn 
   real(r8), pointer :: m_livestemn_storage_to_fire(:) ! (gN/m2/s) fire N emissions from livestemn_storage      
   real(r8), pointer :: m_livestemn_xfer_to_fire(:)    ! (gN/m2/s) fire N emissions from livestemn_xfer
   real(r8), pointer :: m_deadstemn_to_fire(:)         ! (gN/m2/s) fire N emissions from deadstemn
   real(r8), pointer :: m_deadstemn_storage_to_fire(:) ! (gN/m2/s) fire N emissions from deadstemn_storage         
   real(r8), pointer :: m_deadstemn_xfer_to_fire(:)    ! (gN/m2/s) fire N emissions from deadstemn_xfer
   real(r8), pointer :: m_frootn_to_fire(:)            ! (gN/m2/s) fire N emissions from frootn
   real(r8), pointer :: m_frootn_storage_to_fire(:)    ! (gN/m2/s) fire N emissions from frootn_storage
   real(r8), pointer :: m_frootn_xfer_to_fire(:)       ! (gN/m2/s) fire N emissions from frootn_xfer
   real(r8), pointer :: m_livecrootn_to_fire(:)        ! (gN/m2/s) fire N emissions from m_livecrootn_to_fire
   real(r8), pointer :: m_livecrootn_storage_to_fire(:)! (gN/m2/s) fire N emissions from livecrootn_storage     
   real(r8), pointer :: m_livecrootn_xfer_to_fire(:)   ! (gN/m2/s) fire N emissions from livecrootn_xfer
   real(r8), pointer :: m_deadcrootn_to_fire(:)        ! (gN/m2/s) fire N emissions from deadcrootn
   real(r8), pointer :: m_deadcrootn_storage_to_fire(:)! (gN/m2/s) fire N emissions from deadcrootn_storage  
   real(r8), pointer :: m_deadcrootn_xfer_to_fire(:)   ! (gN/m2/s) fire N emissions from deadcrootn_xfer
   real(r8), pointer :: m_retransn_to_fire(:)          ! (gN/m2/s) fire N emissions from retransn
  real(r8), pointer :: m_leafn_to_litter_fire(:)       ! (gN/m2/s) from leafn to litter N  due to fire               
   real(r8), pointer :: m_leafn_storage_to_litter_fire(:) ! (gN/m2/s) from leafn_storage to litter N  due to fire                              
   real(r8), pointer :: m_leafn_xfer_to_litter_fire(:)    ! (gN/m2/s) from leafn_xfer to litter N  due to fire                              
   real(r8), pointer :: m_livestemn_to_litter_fire(:)      ! (gN/m2/s) from livestemn to litter N  due to fire                              
   real(r8), pointer :: m_livestemn_storage_to_litter_fire(:) ! (gN/m2/s) from livestemn_storage to litter N  due to fire                                     
   real(r8), pointer :: m_livestemn_xfer_to_litter_fire(:)    ! (gN/m2/s) from livestemn_xfer to litter N  due to fire                                     
   real(r8), pointer :: m_livestemn_to_deadstemn_fire(:)      ! (gN/m2/s) from livestemn to deadstemn N  due to fire                                     
   real(r8), pointer :: m_deadstemn_to_litter_fire(:)         ! (gN/m2/s) from deadstemn to litter N  due to fire                                     
   real(r8), pointer :: m_deadstemn_storage_to_litter_fire(:) ! (gN/m2/s) from deadstemn_storage to litter N  due to fire                                               
   real(r8), pointer :: m_deadstemn_xfer_to_litter_fire(:)    ! (gN/m2/s) from deadstemn_xfer to litter N  due to fire                                               
   real(r8), pointer :: m_frootn_to_litter_fire(:)            ! (gN/m2/s) from frootn to litter N  due to fire                                               
   real(r8), pointer :: m_frootn_storage_to_litter_fire(:)    ! (gN/m2/s) from frootn_storage to litter N  due to fire                                               
   real(r8), pointer :: m_frootn_xfer_to_litter_fire(:)       ! (gN/m2/s) from frootn_xfer to litter N  due to fire                                               
   real(r8), pointer :: m_livecrootn_to_litter_fire(:)        ! (gN/m2/s) from livecrootn to litter N  due to fire                                               
   real(r8), pointer :: m_livecrootn_storage_to_litter_fire(:)! (gN/m2/s) from livecrootn_storage to litter N  due to fire                                                     
   real(r8), pointer :: m_livecrootn_xfer_to_litter_fire(:)   ! (gN/m2/s) from livecrootn_xfer to litter N  due to fire                                                     
   real(r8), pointer :: m_livecrootn_to_deadcrootn_fire(:)    ! (gN/m2/s) from livecrootn_xfer to deadcrootn due to fire                                                     
   real(r8), pointer :: m_deadcrootn_to_litter_fire(:)        ! (gN/m2/s) from deadcrootn to deadcrootn due to fire                                                       
   real(r8), pointer :: m_deadcrootn_storage_to_litter_fire(:)! (gN/m2/s) from deadcrootn_storage to deadcrootn due to fire                                                        
   real(r8), pointer :: m_deadcrootn_xfer_to_litter_fire(:)   ! (gN/m2/s) from deadcrootn_xfer to deadcrootn due to fire                                                         
   real(r8), pointer :: m_retransn_to_litter_fire(:)          ! (gN/m2/s) from retransn to deadcrootn due to fire                                                         
   

   ! phenology fluxes from transfer pool
   real(r8), pointer :: grainn_xfer_to_grainn(:)            ! grain N growth from storage for prognostic crop model (gN/m2/s)
   real(r8), pointer :: leafn_xfer_to_leafn(:)              ! leaf N growth from storage (gN/m2/s)
   real(r8), pointer :: frootn_xfer_to_frootn(:)            ! fine root N growth from storage (gN/m2/s)
   real(r8), pointer :: livestemn_xfer_to_livestemn(:)      ! live stem N growth from storage (gN/m2/s)
   real(r8), pointer :: deadstemn_xfer_to_deadstemn(:)      ! dead stem N growth from storage (gN/m2/s)
   real(r8), pointer :: livecrootn_xfer_to_livecrootn(:)    ! live coarse root N growth from storage (gN/m2/s)
   real(r8), pointer :: deadcrootn_xfer_to_deadcrootn(:)    ! dead coarse root N growth from storage (gN/m2/s)
   ! litterfall fluxes
   real(r8), pointer :: livestemn_to_litter(:)              ! livestem N to litter (gN/m2/s)
   real(r8), pointer :: grainn_to_food(:)                   ! grain N to food for prognostic crop (gN/m2/s)
   real(r8), pointer :: leafn_to_litter(:)                  ! leaf N litterfall (gN/m2/s)
   real(r8), pointer :: leafn_to_retransn(:)                ! leaf N to retranslocated N pool (gN/m2/s)
   real(r8), pointer :: frootn_to_retransn(:)               ! fine root N to retranslocated N pool (gN/m2/s)
   real(r8), pointer :: frootn_to_litter(:)                 ! fine root N litterfall (gN/m2/s)
   ! allocation fluxes
   real(r8), pointer :: retransn_to_npool(:)                ! deployment of retranslocated N (gN/m2/s)       
   real(r8), pointer :: sminn_to_npool(:)                   ! deployment of soil mineral N uptake (gN/m2/s)
   real(r8), pointer :: npool_to_grainn(:)                  ! allocation to grain N for prognostic crop (gN/m2/s)
   real(r8), pointer :: npool_to_grainn_storage(:)          ! allocation to grain N storage for prognostic crop (gN/m2/s)
   real(r8), pointer :: npool_to_leafn(:)                   ! allocation to leaf N (gN/m2/s)
   real(r8), pointer :: npool_to_leafn_storage(:)           ! allocation to leaf N storage (gN/m2/s)
   real(r8), pointer :: npool_to_frootn(:)                  ! allocation to fine root N (gN/m2/s)
   real(r8), pointer :: npool_to_frootn_storage(:)          ! allocation to fine root N storage (gN/m2/s)
   real(r8), pointer :: npool_to_livestemn(:)               ! allocation to live stem N (gN/m2/s)
   real(r8), pointer :: npool_to_livestemn_storage(:)       ! allocation to live stem N storage (gN/m2/s)
   real(r8), pointer :: npool_to_deadstemn(:)               ! allocation to dead stem N (gN/m2/s)
   real(r8), pointer :: npool_to_deadstemn_storage(:)       ! allocation to dead stem N storage (gN/m2/s)
   real(r8), pointer :: npool_to_livecrootn(:)              ! allocation to live coarse root N (gN/m2/s)
   real(r8), pointer :: npool_to_livecrootn_storage(:)      ! allocation to live coarse root N storage (gN/m2/s)
   real(r8), pointer :: npool_to_deadcrootn(:)              ! allocation to dead coarse root N (gN/m2/s)
   real(r8), pointer :: npool_to_deadcrootn_storage(:)      ! allocation to dead coarse root N storage (gN/m2/s)
   ! annual turnover of storage to transfer pools           
   real(r8), pointer :: grainn_storage_to_xfer(:)           ! grain N shift storage to transfer for prognostic crop (gN/m2/s)
   real(r8), pointer :: leafn_storage_to_xfer(:)            ! leaf N shift storage to transfer (gN/m2/s)
   real(r8), pointer :: frootn_storage_to_xfer(:)           ! fine root N shift storage to transfer (gN/m2/s)
   real(r8), pointer :: livestemn_storage_to_xfer(:)        ! live stem N shift storage to transfer (gN/m2/s)
   real(r8), pointer :: deadstemn_storage_to_xfer(:)        ! dead stem N shift storage to transfer (gN/m2/s)
   real(r8), pointer :: livecrootn_storage_to_xfer(:)       ! live coarse root N shift storage to transfer (gN/m2/s)
   real(r8), pointer :: deadcrootn_storage_to_xfer(:)       ! dead coarse root N shift storage to transfer (gN/m2/s)
   real(r8), pointer :: fert(:)                             ! applied fertilizer (gN/m2/s)
   real(r8), pointer :: soyfixn(:)                          ! soybean fixed N (gN/m2/s)
   ! turnover of livewood to deadwood, with retranslocation 
   real(r8), pointer :: livestemn_to_deadstemn(:)           ! live stem N turnover (gN/m2/s)
   real(r8), pointer :: livestemn_to_retransn(:)            ! live stem N to retranslocated N pool (gN/m2/s)
   real(r8), pointer :: livecrootn_to_deadcrootn(:)         ! live coarse root N turnover (gN/m2/s)
   real(r8), pointer :: livecrootn_to_retransn(:)           ! live coarse root N to retranslocated N pool (gN/m2/s)
   ! summary (diagnostic) flux variables, not involved in mass balance
   real(r8), pointer :: ndeploy(:)                          ! total N deployed to growth and storage (gN/m2/s)
   real(r8), pointer :: pft_ninputs(:)                      ! total N inputs to pft-level (gN/m2/s)
   real(r8), pointer :: pft_noutputs(:)                     ! total N outputs from pft-level (gN/m2/s)
   real(r8), pointer :: wood_harvestn(:)                    ! total N losses to wood product pools (gN/m2/s)
   ! new variables for fire code 
   real(r8), pointer :: pft_fire_nloss(:)                   ! total pft-level fire N loss (gN/m2/s) 
end type pft_nflux_type

!----------------------------------------------------
! pft VOC fluxes structure for history output
!----------------------------------------------------
type, public :: megan_out_type
   real(r8), pointer :: flux_out(:)   !(n_megan_comps) MEGAN flux [ug C m-2 h-1]
endtype megan_out_type

!----------------------------------------------------
! pft VOC flux variables structure
!----------------------------------------------------
type, public :: pft_vflux_type
   real(r8), pointer :: vocflx_tot(:)     !total VOC flux into atmosphere [moles/m2/sec]
   real(r8), pointer :: vocflx(:,:)       !(num_mech_comps) MEGAN flux [moles/m2/sec]
   real(r8), pointer :: Eopt_out(:)       !Eopt coefficient
   real(r8), pointer :: topt_out(:)       !topt coefficient
   real(r8), pointer :: alpha_out(:)      !alpha coefficient
   real(r8), pointer :: cp_out(:)         !cp coefficient
   real(r8), pointer :: paru_out(:)
   real(r8), pointer :: par24u_out(:)
   real(r8), pointer :: par240u_out(:)
   real(r8), pointer :: para_out(:)
   real(r8), pointer :: par24a_out(:)
   real(r8), pointer :: par240a_out(:)
   real(r8), pointer :: gamma_out(:)
   real(r8), pointer :: gammaL_out(:)
   real(r8), pointer :: gammaT_out(:)
   real(r8), pointer :: gammaP_out(:)
   real(r8), pointer :: gammaA_out(:)
   real(r8), pointer :: gammaS_out(:)
   real(r8), pointer :: gammaC_out(:)
   type(megan_out_type), pointer :: meg(:) ! points to output fluxes
end type pft_vflux_type

!----------------------------------------------------
! pft dry dep velocity variables structure
!----------------------------------------------------
type, public :: pft_depvd_type
   real(r8), pointer :: drydepvel(:,:)
end type pft_depvd_type

!----------------------------------------------------
! pft dust flux variables structure
!----------------------------------------------------
type, public :: pft_dflux_type
   real(r8), pointer :: flx_mss_vrt_dst(:,:)    !(ndst)  !surface dust emission (kg/m**2/s) [ + = to atm]
   real(r8), pointer :: flx_mss_vrt_dst_tot(:)  !total dust flux into atmosphere
   real(r8), pointer :: vlc_trb(:,:)            !(ndst) turbulent deposition velocity (m/s)
   real(r8), pointer :: vlc_trb_1(:)            !turbulent deposition velocity 1(m/s)
   real(r8), pointer :: vlc_trb_2(:)            !turbulent deposition velocity 2(m/s)
   real(r8), pointer :: vlc_trb_3(:)            !turbulent deposition velocity 3(m/s)
   real(r8), pointer :: vlc_trb_4(:)            !turbulent deposition velocity 4(m/s)
end type pft_dflux_type

!----------------------------------------------------
! End definition of structures defined at the pft_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the column_type level
!----------------------------------------------------
! column physical state variables structure
!----------------------------------------------------
type, public :: column_pstate_type
   type(pft_pstate_type) :: pps_a            !pft-level pstate variables averaged to the column
   integer , pointer :: snl(:)                !number of snow layers
   integer , pointer :: isoicol(:)            !soil color class

   !F. Li and S. Levis
   real(r8), pointer :: gdp_lf(:)             ! global real gdp data (k US$/capita)
   real(r8), pointer :: peatf_lf(:)           ! global peatland fraction data (0-1)
   integer, pointer :: abm_lf(:)             ! global peak month of crop fire emissions 
   real(r8), pointer :: lgdp_col(:)           !gdp limitation factor for fire occurrence (0-1)
   real(r8), pointer :: lgdp1_col(:)          !gdp limitation factor for fire spreading (0-1)
   real(r8), pointer :: lpop_col(:)           !pop limitation factor for fire spreading (0-1)

   real(r8), pointer :: bsw(:,:)              !Clapp and Hornberger "b" (nlevgrnd)  
   real(r8), pointer :: watsat(:,:)           !volumetric soil water at saturation (porosity) (nlevgrnd) 
   real(r8), pointer :: watdry(:,:)           !btran parameter for btran=0
   real(r8), pointer :: watopt(:,:)           !btran parameter for btran = 1
   real(r8), pointer :: hksat(:,:)            !hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd) 
   real(r8), pointer :: hksat_min(:,:)        !mineral hksat
   real(r8), pointer :: tk_hist(:,:)          !thermal conductivity
   real(r8), pointer :: cv_hist(:,:)          !heat capacity
   real(r8), pointer :: sucsat(:,:)           !minimum soil suction (mm) (nlevgrnd) 
   real(r8), pointer :: hkdepth(:)            !decay factor (m)
   real(r8), pointer :: wtfact(:)             !maximum saturated fraction for a gridcell
   real(r8), pointer :: fracice(:,:)          !fractional impermeability (-)
   real(r8), pointer :: csol(:,:)             !heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd) 
   real(r8), pointer :: tkmg(:,:)             !thermal conductivity, soil minerals  [W/m-K] (new) (nlevgrnd) 
   real(r8), pointer :: tkdry(:,:)            !thermal conductivity, dry soil (W/m/Kelvin) (nlevgrnd) 
   real(r8), pointer :: tksatu(:,:)           !thermal conductivity, saturated soil [W/m-K] (new) (nlevgrnd) 
   real(r8), pointer :: smpmin(:)             !restriction for min of soil potential (mm) (new)
   real(r8), pointer :: gwc_thr(:)            !threshold soil moisture based on clay content
   real(r8), pointer :: mss_frc_cly_vld(:)    ![frc] Mass fraction clay limited to 0.20
   real(r8), pointer :: mbl_bsn_fct(:)        !basin factor
   logical , pointer :: do_capsnow(:)         !true => do snow capping
   real(r8), pointer :: snow_depth(:)             !snow height of snow covered area (m)
   real(r8), pointer :: snowdp(:)             ! gridcell averaged snow height (m)
   real(r8), pointer :: frac_sno(:)           !fraction of ground covered by snow (0 to 1)
   real(r8), pointer :: frac_sno_eff(:)       !fraction of ground covered by snow (0 to 1)
   real(r8), pointer :: frac_h2osfc(:)        !fractional area with surface water greater than zero
   real(r8), pointer :: frac_h2osfc_temp(:)   !temporay fractional area with surface water greater than zero
   real(r8), pointer :: topo_std(:)           !gridcell topographic standard deviation (m)
   real(r8), pointer :: topo_ndx(:)           !gridcell topographic index
   real(r8), pointer :: topo_slope(:)         !gridcell topographic slope
   real(r8), pointer :: micro_sigma(:)        ! microtopography pdf sigma (m)
   real(r8), pointer :: h2osfc_thresh(:)      ! level at which h2osfc "percolates"
   real(r8), pointer :: n_melt(:)             ! SCA shape parameter
   real(r8), pointer :: zi(:,:)               !interface level below a "z" level (m) (-nlevsno+0:nlevgrnd) 
   real(r8), pointer :: dz(:,:)               !layer thickness (m)  (-nlevsno+1:nlevgrnd) 
   real(r8), pointer :: z(:,:)                !layer depth (m) (-nlevsno+1:nlevgrnd) 
   real(r8), pointer :: frac_iceold(:,:)      !fraction of ice relative to the tot water (new) (-nlevsno+1:nlevgrnd) 
   integer , pointer :: imelt(:,:)            !flag for melting (=1), freezing (=2), Not=0 (new) (-nlevsno+1:nlevgrnd) 
   real(r8), pointer :: eff_porosity(:,:)     !effective porosity = porosity - vol_ice (nlevgrnd) 
   real(r8), pointer :: emg(:)                !ground emissivity
   real(r8), pointer :: z0mg(:)               !roughness length over ground, momentum [m]
   real(r8), pointer :: z0hg(:)               !roughness length over ground, sensible heat [m]
   real(r8), pointer :: z0qg(:)               !roughness length over ground, latent heat [m]
   real(r8), pointer :: htvp(:)               !latent heat of vapor of water (or sublimation) [j/kg]
   real(r8), pointer :: beta(:)               !coefficient of convective velocity [-]
   real(r8), pointer :: zii(:)                !convective boundary height [m]
   real(r8), pointer :: albgrd(:,:)           !ground albedo (direct) (numrad)
   real(r8), pointer :: albgri(:,:)           !ground albedo (diffuse) (numrad)
   real(r8), pointer :: rootr_column(:,:)     !effective fraction of roots in each soil layer (nlevgrnd)  
   real(r8), pointer :: rootfr_road_perv(:,:) !fraction of roots in each soil layer for urban pervious road
   real(r8), pointer :: rootr_road_perv(:,:)  !effective fraction of roots in each soil layer of urban pervious road

   real(r8), pointer :: wf(:)                 !soil water as frac. of whc for top 0.05 m (0-1) (only comment changed by F. Li and S. Levis)
   real(r8), pointer :: wf2(:)                !soil water as frac. of whc for top 0.17 m (0-1) added by F. Li and S. Levis

!  real(r8), pointer :: xirrig(:)             !irrigation rate
   real(r8), pointer :: max_dayl(:)           !maximum daylength for this column (s)
#if (defined VICHYDRO)
   real(r8), pointer :: b_infil(:)          !b infiltration parameter
   real(r8), pointer :: ds(:)               !fracton of Dsmax where non-linear baseflow begins
   real(r8), pointer :: dsmax(:)            !max. velocity of baseflow (mm/day)
   real(r8), pointer :: Wsvic(:)            !fraction of maximum soil moisutre where non-liear base flow occurs
   real(r8), pointer :: c_param(:)          !baseflow exponent (Qb)
   real(r8), pointer :: expt(:,:)           !pore-size distribution related paramter(Q12)
   real(r8), pointer :: ksat(:,:)           !Saturated hydrologic conductivity
   real(r8), pointer :: phi_s(:,:)          !soil moisture dissusion parameter
   real(r8), pointer :: depth(:,:)          !layer depth of upper layer 
   real(r8), pointer :: porosity(:,:)       !soil porisity (1-bulk_density/soil_density)
   real(r8), pointer :: max_moist(:,:)      !max layer moist + ice (mm)
   real(r8), pointer :: vic_clm_fract(:,:,:)!fraction of VIC layers in CLM layers
#endif
   ! new variables for CN code
   real(r8), pointer :: decl(:)               ! solar declination angle (radians)
   real(r8), pointer :: coszen(:)             ! cosine of solar zenith angle
   real(r8), pointer :: soilpsi(:,:)          ! soil water potential in each soil layer (MPa)
   real(r8), pointer :: bd(:,:)               ! bulk density of dry soil material [kg/m^3]
   real(r8), pointer :: fpi_vr(:,:)           ! fraction of potential immobilization (no units)
   real(r8), pointer :: fpi(:)                ! fraction of potential immobilization (no units)
   real(r8), pointer :: rf_decomp_cascade(:,:,:)              ! respired fraction in decomposition step (frac)
   real(r8), pointer :: pathfrac_decomp_cascade(:,:,:)        ! what fraction of C leaving a given pool passes through a given transition (frac)
   real(r8), pointer :: nfixation_prof(:,:)                   ! (1/m) profile for N fixation additions
   real(r8), pointer :: ndep_prof(:,:)                        ! (1/m) profile for N fixation additions
   real(r8), pointer :: alt(:)                                ! current depth of thaw
   real(r8), pointer :: altmax(:)                             ! maximum annual depth of thaw
   real(r8), pointer :: altmax_lastyear(:)                    ! prior year maximum annual depth of thaw
   integer, pointer :: alt_indx(:)                            ! current depth of thaw
   integer, pointer :: altmax_indx(:)                         ! maximum annual depth of thaw
   integer, pointer :: altmax_lastyear_indx(:)                ! prior year maximum annual depth of thaw
   real(r8), pointer :: som_adv_coef(:,:)                     ! SOM advective flux (m/s)
   real(r8), pointer :: som_diffus_coef(:,:)                  ! SOM diffusivity due to bio/cryo-turbation (m2/s)

#ifdef NITRIF_DENITRIF
   real(r8), pointer :: tmean_monthly_max_vr(:,:)     ! maximumn monthly-mean soil temperature
   real(r8), pointer :: tmean_monthly_vr(:,:)         ! monthly-mean soil temperature
#endif
   real(r8), pointer :: fpg(:)           !fraction of potential gpp (no units)
   real(r8), pointer :: annsum_counter(:) !seconds since last annual accumulator turnover
   real(r8), pointer :: cannsum_npp(:)    !annual sum of NPP, averaged from pft-level (gC/m2/yr)
   real(r8), pointer :: col_lag_npp(:)        ! (gC/m2/s) lagged net primary production
   real(r8), pointer :: cannavg_t2m(:)    !annual average of 2m air temperature, averaged from pft-level (K)
   real(r8), pointer :: watfc(:,:)        !volumetric soil water at field capacity (nlevsoi)

   ! F. Li and S. Levis
   real(r8), pointer :: nfire(:)        ! fire counts (count/km2/timestep), valid only in Reg. C
   real(r8), pointer :: fsr_pft(:)      ! fire spread rate in pft level (m/s)
   real(r8), pointer :: fsr_col(:)      ! fire spread rate at column level (m/s)
   real(r8), pointer :: fd_col(:)       ! fire duration at column level (hr)
   real(r8), pointer :: fd_pft(:)       ! fire duration in pft level    (hr)
   real(r8), pointer :: prec60_col(:)              !60-day running mean of tot. precipitation (mm/s)
   real(r8), pointer :: prec10_col(:)              !10-day running mean of tot. precipitation (mm/s)
   real(r8), pointer :: lfc(:)         ! conversion area fraction of BET and BDT that haven't burned before (0-1)
   real(r8), pointer :: lfc2(:)        ! conversion area fraction of BET and BDT that burned in this timestep ((timestep)-1)
   real(r8), pointer :: dtrotr_col(:)   ! annual decreased fraction coverage of BET on the gridcell (0-1)
   real(r8), pointer :: trotr1_col(:)   ! pft weight of BET and BDT on the gridcell(0-1)
   real(r8), pointer :: trotr2_col(:)   ! pft weight of BDT on the gridcell (0-1)
   real(r8), pointer :: cropf_col(:)    ! crop fraction in veg column (0-1)
   real(r8), pointer :: baf_crop(:)     ! baf for cropland per time step(0-1)
   real(r8), pointer :: baf_peatf(:)    ! baf for peatland per time step (0-1)
   real(r8), pointer :: fbac(:)         ! total burned area out of conversion (0-1)
   real(r8), pointer :: fbac1(:)        ! burned area out of conversion region due to land use fire (0-1)
   real(r8), pointer :: btran_col(:)    ! btran2 at column level (0-1)
   real(r8), pointer :: wtlf(:)         ! fractional coverage of non-crop PFTs (0-1)
   real(r8), pointer :: lfwt(:)         ! fractional coverage of non-crop and non-bare-soil PFTs (0-1)
   real(r8), pointer :: farea_burned(:)       !timestep fractional area burned (0-1) 


   real(r8), pointer :: albsnd_hst(:,:)       ! snow albedo, direct, for history files (col,bnd) [frc]
   real(r8), pointer :: albsni_hst(:,:)       ! snow albedo, diffuse, for history files (col,bnd) [frc]
   real(r8), pointer :: albsod(:,:)           ! soil albedo: direct (col,bnd) [frc]
   real(r8), pointer :: albsoi(:,:)           ! soil albedo: diffuse (col,bnd) [frc]
   real(r8), pointer :: flx_absdv(:,:)        ! absorbed flux per unit incident direct flux: VIS (col,lyr) [frc]
   real(r8), pointer :: flx_absdn(:,:)        ! absorbed flux per unit incident direct flux: NIR (col,lyr) [frc]
   real(r8), pointer :: flx_absiv(:,:)        ! absorbed flux per unit incident diffuse flux: VIS (col,lyr) [frc]
   real(r8), pointer :: flx_absin(:,:)        ! absorbed flux per unit incident diffuse flux: NIR (col,lyr) [frc]
   real(r8), pointer :: snw_rds(:,:)          ! snow grain radius (col,lyr) [m^-6, microns]
   real(r8), pointer :: snw_rds_top(:)        ! snow grain radius, top layer (col) [m^-6, microns]
   real(r8), pointer :: sno_liq_top(:)        ! snow liquid water fraction (mass), top layer (col) [fraction]
   real(r8), pointer :: mss_bcpho(:,:)        ! mass of hydrophobic BC in snow (col,lyr) [kg]
   real(r8), pointer :: mss_bcphi(:,:)        ! mass of hydrophillic BC in snow (col,lyr) [kg]
   real(r8), pointer :: mss_bctot(:,:)        ! total mass of BC in snow (pho+phi) (col,lyr) [kg]
   real(r8), pointer :: mss_bc_col(:)         ! column-integrated mass of total BC (col) [kg]
   real(r8), pointer :: mss_bc_top(:)         ! top-layer mass of total BC (col) [kg]
   real(r8), pointer :: mss_ocpho(:,:)        ! mass of hydrophobic OC in snow (col,lyr) [kg]
   real(r8), pointer :: mss_ocphi(:,:)        ! mass of hydrophillic OC in snow (col,lyr) [kg]
   real(r8), pointer :: mss_octot(:,:)        ! total mass of OC in snow (pho+phi) (col,lyr) [kg]
   real(r8), pointer :: mss_oc_col(:)         ! column-integrated mass of total OC (col) [kg]
   real(r8), pointer :: mss_oc_top(:)         ! top-layer mass of total OC (col) [kg]
   real(r8), pointer :: mss_dst1(:,:)         ! mass of dust species 1 in snow (col,lyr) [kg]
   real(r8), pointer :: mss_dst2(:,:)         ! mass of dust species 2 in snow (col,lyr) [kg]
   real(r8), pointer :: mss_dst3(:,:)         ! mass of dust species 3 in snow (col,lyr) [kg]
   real(r8), pointer :: mss_dst4(:,:)         ! mass of dust species 4 in snow (col,lyr) [kg]
   real(r8), pointer :: mss_dsttot(:,:)       ! total mass of dust in snow (col,lyr) [kg]
   real(r8), pointer :: mss_dst_col(:)        ! column-integrated mass of dust in snow (col) [kg]
   real(r8), pointer :: mss_dst_top(:)        ! top-layer mass of dust in snow (col) [kg]
   real(r8), pointer :: h2osno_top(:)         ! top-layer mass of snow (col) [kg]
   real(r8), pointer :: mss_cnc_bcphi(:,:)    ! mass concentration of hydrophilic BC in snow (col,lyr) [kg/kg]
   real(r8), pointer :: mss_cnc_bcpho(:,:)    ! mass concentration of hydrophilic BC in snow (col,lyr) [kg/kg]
   real(r8), pointer :: mss_cnc_ocphi(:,:)    ! mass concentration of hydrophilic OC in snow (col,lyr) [kg/kg]
   real(r8), pointer :: mss_cnc_ocpho(:,:)    ! mass concentration of hydrophilic OC in snow (col,lyr) [kg/kg]
   real(r8), pointer :: mss_cnc_dst1(:,:)     ! mass concentration of dust species 1 in snow (col,lyr) [kg/kg]
   real(r8), pointer :: mss_cnc_dst2(:,:)     ! mass concentration of dust species 2 in snow (col,lyr) [kg/kg]
   real(r8), pointer :: mss_cnc_dst3(:,:)     ! mass concentration of dust species 3 in snow (col,lyr) [kg/kg]
   real(r8), pointer :: mss_cnc_dst4(:,:)     ! mass concentration of dust species 4 in snow (col,lyr) [kg/kg]
   real(r8), pointer :: albgrd_pur(:,:)       ! pure snow ground direct albedo (numrad)
   real(r8), pointer :: albgri_pur(:,:)       ! pure snow ground diffuse albedo (numrad)
   real(r8), pointer :: albgrd_bc(:,:)        ! ground direct albedo without BC  (numrad)
   real(r8), pointer :: albgri_bc(:,:)        ! ground diffuse albedo without BC (numrad)
   real(r8), pointer :: albgrd_oc(:,:)        ! ground direct albedo without OC  (numrad)
   real(r8), pointer :: albgri_oc(:,:)        ! ground diffuse albedo without OC (numrad)
   real(r8), pointer :: albgrd_dst(:,:)       ! ground direct albedo without dust  (numrad)
   real(r8), pointer :: albgri_dst(:,:)       ! ground diffuse albedo without dust (numrad)
   real(r8), pointer :: dTdz_top(:)           ! temperature gradient in top layer  [K m-1]
   real(r8), pointer :: snot_top(:)           ! temperature of top snow layer [K]
   ! new variables for S Lake code
   real(r8), pointer :: ws(:)              ! surface friction velocity (m/s)
   real(r8), pointer :: ks(:)              ! coefficient for calculation of decay of eddy diffusivity with depth
   real(r8), pointer :: dz_lake(:,:)       ! lake layer thickness (m)  (1:nlevlak)
   real(r8), pointer :: z_lake(:,:)        ! layer depth for lake (m)
   real(r8), pointer :: savedtke1(:)       ! top level eddy conductivity from previous timestep (W/mK)
   real(r8), pointer :: cellsand(:,:)      ! sand value for gridcell containing column (1:nlevsoi)
   real(r8), pointer :: cellclay(:,:)      ! clay value for gridcell containing column (1:nlevsoi)
   real(r8), pointer :: cellorg(:,:)       ! organic matter for gridcell containing column (1:nlevsoi)
   real(r8), pointer :: lakedepth(:)       ! variable lake depth (m)
   real(r8), pointer :: etal(:)            ! lake extinction coefficient from surface data (1/m)
   real(r8), pointer :: lakefetch(:)       ! lake fetch from surface data (m)
   real(r8), pointer :: ust_lake(:)        ! friction velocity (m/s)
   ! end new variables for S Lake code
   ! New variables for finundated in methane code
#ifdef LCH4
   real(r8), pointer :: zwt0(:)             ! coefficient for determining finundated (m)
   real(r8), pointer :: f0(:)               ! maximum inundated fraction for a gridcell (for methane code)
   real(r8), pointer :: p3(:)               ! coefficient for determining finundated (m)
   ! added by Lei Meng for pH effects of methane production
   real(r8), pointer :: pH(:)               ! pH values
#endif
   ! End New variables for methane code
   real(r8), pointer :: irrig_rate(:)         ! current irrigation rate [mm/s]
   integer, pointer  :: n_irrig_steps_left(:) ! number of time steps for which we still need to irrigate today (if 0, ignore irrig_rate)
   real(r8), pointer :: forc_pbot(:)          ! surface atm pressure, downscaled to column (Pa)
   real(r8), pointer :: forc_rho(:)           ! surface air density, downscaled to column (kg/m^3)
   real(r8), pointer :: glc_frac(:)           ! ice fractional area
   real(r8), pointer :: glc_topo(:)           ! surface elevation (m)
end type column_pstate_type

!----------------------------------------------------
! column energy state variables structure
!----------------------------------------------------
type, public :: column_estate_type
   type(pft_estate_type):: pes_a              !pft-level energy state variables averaged to the column
   real(r8), pointer :: t_grnd(:)             !ground temperature (Kelvin)
   real(r8), pointer :: t_grnd_u(:)           !Urban ground temperature (Kelvin)
   real(r8), pointer :: t_grnd_r(:)           !Rural ground temperature (Kelvin)
   real(r8), pointer :: dt_grnd(:)            !change in t_grnd, last iteration (Kelvin)
   real(r8), pointer :: t_soisno(:,:)         !soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd) 
   real(r8), pointer :: t_soi_10cm(:)         !soil temperature in top 10cm of soil (Kelvin)
   real(r8), pointer :: tsoi17(:)            !soil temperature in top 17cm of soil (Kelvin) by F. Li and S. Levis
   real(r8), pointer :: t_lake(:,:)           !lake temperature (Kelvin)  (1:nlevlak)          
   real(r8), pointer :: tssbef(:,:)           !soil/snow temperature before update (-nlevsno+1:nlevgrnd) 
   real(r8), pointer :: thv(:)                !virtual potential temperature (kelvin)
   real(r8), pointer :: hc_soi(:)             !soil heat content (MJ/m2)
   real(r8), pointer :: hc_soisno(:)          !soil plus snow heat content (MJ/m2)
   real(r8), pointer :: forc_t(:)             !atm temperature, downscaled to column (Kelvin)
   real(r8), pointer :: forc_th(:)            !atm potl temperature, downscaled to column (Kelvin)
   real(r8), pointer :: t_h2osfc(:) 	      !surface water temperature
   real(r8), pointer :: t_h2osfc_bef(:)       !surface water temperature from time-step before
end type column_estate_type

!----------------------------------------------------
! column water state variables structure
!----------------------------------------------------
type, public :: column_wstate_type
   type(pft_wstate_type):: pws_a              !pft-level water state variables averaged to the column
   real(r8), pointer :: h2osfc(:)             !surface water (mm H2O)
   real(r8), pointer :: qg_snow(:)            !ground specific humidity [kg/kg]
   real(r8), pointer :: qg_soil(:)            !ground specific humidity [kg/kg]
   real(r8), pointer :: qg_h2osfc(:)          !ground specific humidity [kg/kg]
   real(r8), pointer :: swe_old(:,:)          !initial snow water
   real(r8), pointer :: h2osno(:)             !snow water (mm H2O)
   real(r8), pointer :: errh2osno(:)          !imbalance in snow water (mm H2O)
   real(r8), pointer :: snow_sources(:)       !snow sources (mm H2O/s)
   real(r8), pointer :: snow_sinks(:)         !snow sinks (mm H2O/s)
   real(r8), pointer :: h2osoi_liq(:,:)       !liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
   real(r8), pointer :: h2osoi_ice(:,:)       !ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
   real(r8), pointer :: h2osoi_liqice_10cm(:) !liquid water + ice lens in top 10cm of soil (kg/m2)
   real(r8), pointer :: h2osoi_vol(:,:)       !volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)  
   real(r8), pointer :: h2osno_old(:)         !snow mass for previous time step (kg/m2) (new)
   real(r8), pointer :: qg(:)                 !ground specific humidity [kg/kg]
   real(r8), pointer :: dqgdT(:)              !d(qg)/dT
   real(r8), pointer :: snowice(:)            !average snow ice lens
   real(r8), pointer :: snowliq(:)            !average snow liquid water
   real(r8) ,pointer :: soilalpha(:)          !factor that reduces ground saturated specific humidity (-)
   real(r8), pointer :: soilbeta(:)           !factor that reduces ground evaporation L&P1992(-)
   real(r8) ,pointer :: soilalpha_u(:)        !urban factor that reduces ground saturated specific humidity (-)
   real(r8), pointer :: zwt(:)                !water table depth
   real(r8), pointer :: frost_table(:)        !frost table depth
   real(r8), pointer :: zwt_perched(:)        !perched water table depth
   real(r8), pointer :: int_snow(:)           !integrated snowfall (mm H2O)
   real(r8), pointer :: fcov(:)               !fractional impermeable area
   real(r8), pointer :: wa(:)                 !water in the unconfined aquifer (mm)
   real(r8), pointer :: qcharge(:)            !aquifer recharge rate (mm/s)
   real(r8), pointer :: smp_l(:,:)            !soil matric potential (mm)
   real(r8), pointer :: hk_l(:,:)             !hydraulic conductivity (mm/s)
   real(r8), pointer :: fsat(:)               !fractional area with water table at surface
   real(r8), pointer :: forc_q(:)             !atm specific humidity, downscaled to column (kg/kg)
#if (defined VICHYDRO)
   real(r8), pointer :: moist(:,:)            !soil moisture (kg/m2) for VIC soil layers
   real(r8), pointer :: ice(:,:)              !soil ice (kg/m2) for VIC soil layers
   real(r8), pointer :: moist_vol(:,:)        !volumetric soil moisture for VIC soil layers
   real(r8), pointer :: max_infil(:)          !maximum infiltration rate calculated by VIC
   real(r8), pointer :: i_0(:)                !average saturation in top soil layers in VIC
#endif
#ifdef LCH4
   real(r8), pointer :: finundated(:)         !fractional inundated area (excluding dedicated wetland columns)
#endif
   ! new variables for S Lake code
   real(r8), pointer :: lake_icefrac(:,:)     ! mass fraction of lake layer that is frozen
   real(r8), pointer :: lake_icethick(:)      ! ice thickness (m) (integrated if lakepuddling)
   ! end new variables for S Lake code
end type column_wstate_type

!----------------------------------------------------
! column carbon state variables structure
!----------------------------------------------------
type, public :: column_cstate_type
   type(pft_cstate_type):: pcs_a              !pft-level carbon state variables averaged to the column
   ! NOTE: the soilc variable is used by the original CLM C-cycle code,
   ! and is not used by the CN code
   real(r8), pointer :: soilc(:)              !soil carbon (kg C /m**2)
   ! all c pools involved in decomposition
   real(r8), pointer :: decomp_cpools_vr(:,:,:)         ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
   real(r8), pointer :: col_ctrunc_vr(:,:)              ! (gC/m3) vertically-resolved column-level sink for C truncation

   !fire-related variables added by F. Li and S. Levis
   real(r8), pointer :: rootc_col(:)     !root carbon at column level (gC/m2)
   real(r8), pointer :: totvegc_col(:)   !column-level totvegc (gC/m2)  
   real(r8), pointer :: leafc_col(:)     !column-level leafc (gC/m2)       
   real(r8), pointer :: fuelc(:)         ! fuel avalability factor for Reg.C (0-1)
   real(r8), pointer :: fuelc_crop(:)    ! fuel avalability factor for Reg.A (0-1)

  ! pools for dynamic landcover
   real(r8), pointer :: seedc(:)                        ! (gC/m2) column-level pool for seeding new PFTs
   real(r8), pointer :: prod10c(:)                      ! (gC/m2) wood product C pool, 10-year lifespan
   real(r8), pointer :: prod100c(:)                     ! (gC/m2) wood product C pool, 100-year lifespan
   real(r8), pointer :: totprodc(:)                     ! (gC/m2) total wood product C
   ! summary (diagnostic) state variables, not involved in mass balance
   real(r8), pointer :: decomp_cpools(:,:)              ! (gC/m2)  decomposing (litter, cwd, soil) c pools
   real(r8), pointer :: decomp_cpools_1m(:,:)           ! (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
   real(r8), pointer :: cwdc(:)                         ! (gC/m2) Diagnostic: coarse woody debris C
   real(r8), pointer :: col_ctrunc(:)                   ! (gC/m2) column-level sink for C truncation
   real(r8), pointer :: totlitc(:)                      ! (gC/m2) total litter carbon
   real(r8), pointer :: totsomc(:)                      ! (gC/m2) total soil organic matter carbon
   real(r8), pointer :: totlitc_1m(:)                   ! (gC/m2) total litter carbon to 1 meter
   real(r8), pointer :: totsomc_1m(:)                   ! (gC/m2) total soil organic matter carbon to 1 meter
   real(r8), pointer :: totecosysc(:)                   ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
   real(r8), pointer :: totcolc(:)                      ! (gC/m2) total column carbon, incl veg and cpool
end type column_cstate_type

!----------------------------------------------------
! column methane variables structure
!----------------------------------------------------
#ifdef LCH4
type, public :: column_ch4_type
   ! new variables for CH4 code
   ! column-level methane fluxes
   real(r8), pointer :: ch4_prod_depth_sat(:,:) ! CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_prod_depth_unsat(:,:) ! CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_prod_depth_lake(:,:)! CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_oxid_depth_sat(:,:) ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_oxid_depth_unsat(:,:) !CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_oxid_depth_lake(:,:) ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_aere_depth_sat(:,:) ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_aere_depth_unsat(:,:) ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_tran_depth_sat(:,:) ! CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_tran_depth_unsat(:,:) ! CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_ebul_depth_sat(:,:) ! CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_ebul_depth_unsat(:,:) ! CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: ch4_ebul_total_sat(:)   ! Total column CH4 ebullition (mol/m2/s)
   real(r8), pointer :: ch4_ebul_total_unsat(:)   ! Total column CH4 ebullition (mol/m2/s)
   real(r8), pointer :: ch4_surf_aere_sat(:)   ! CH4 aerenchyma flux to atmosphere (after oxidation) (mol/m2/s)
   real(r8), pointer :: ch4_surf_aere_unsat(:)   ! CH4 aerenchyma flux to atmosphere (after oxidation) (mol/m2/s)
   real(r8), pointer :: ch4_surf_ebul_sat(:)   ! CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
   real(r8), pointer :: ch4_surf_ebul_unsat(:)   ! CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
   real(r8), pointer :: ch4_surf_ebul_lake(:)   ! CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
   real(r8), pointer :: co2_aere_depth_sat(:,:) ! CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: co2_aere_depth_unsat(:,:) ! CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: o2_oxid_depth_sat(:,:)  ! O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: o2_oxid_depth_unsat(:,:)  ! O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: o2_aere_depth_sat(:,:)  ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: o2_aere_depth_unsat(:,:)  ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: o2_decomp_depth_sat(:,:) !O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   real(r8), pointer :: o2_decomp_depth_unsat(:,:)!O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   real(r8), pointer :: co2_decomp_depth_sat(:,:)  ! CO2 production during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   real(r8), pointer :: co2_decomp_depth_unsat(:,:)  ! CO2 production during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   real(r8), pointer :: co2_oxid_depth_sat(:,:) ! CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: co2_oxid_depth_unsat(:,:) ! CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
   real(r8), pointer :: conc_o2_sat(:,:)        ! O2 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: conc_o2_unsat(:,:)        ! O2 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: conc_o2_lake(:,:)        ! O2 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: conc_ch4_sat(:,:)       ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: conc_ch4_unsat(:,:)       ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: conc_ch4_lake(:,:)      ! CH4 conc in each soil layer (mol/m3) (nlevsoi)
   real(r8), pointer :: ch4_surf_diff_sat(:)    ! CH4 surface flux (mol/m2/s)
   real(r8), pointer :: ch4_surf_diff_unsat(:)    ! CH4 surface flux (mol/m2/s)
   real(r8), pointer :: ch4_surf_diff_lake(:)   ! CH4 surface flux (mol/m2/s)
   real(r8), pointer :: ch4_dfsat_flux(:)       ! CH4 flux to atm due to decreasing fsat (kg C/m^2/s) [+]
   ! Other variables
   real(r8), pointer :: zwt_ch4_unsat(:) ! depth of water table for unsaturated fraction (m)
   real(r8), pointer :: fsat_bef(:)    !fsat from previous timestep
   real(r8), pointer :: lake_soilc(:,:)      ! total soil organic matter found in level (g C / m^3) (nlevsoi)
   real(r8), pointer :: lake_raw(:)        !aerodynamic resistance for moisture (s/m)
   real(r8), pointer :: totcolch4(:)       ! total methane found in soil column (g C / m^2)
   real(r8), pointer :: fphr(:,:)         ! fraction of potential heterotrophic respiration
   real(r8), pointer :: annsum_counter(:)     ! seconds since last annual accumulator turnover
   real(r8), pointer :: tempavg_somhr(:)         ! temporary average SOM heterotrophic resp. (gC/m2/s)
   real(r8), pointer :: annavg_somhr(:)          ! annual average SOM heterotrophic resp. (gC/m2/s)
   real(r8), pointer :: tempavg_finrw(:)    ! respiration-weighted annual average of finundated
   real(r8), pointer :: annavg_finrw(:)    ! respiration-weighted annual average of finundated
   real(r8), pointer :: sif(:) ! (unitless) ratio applied to sat. prod. to account for seasonal inundation
   real(r8), pointer :: o2stress_unsat(:,:) ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
   real(r8), pointer :: o2stress_sat(:,:) ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
   real(r8), pointer :: ch4stress_unsat(:,:) ! Ratio of methane available to the total per-timestep methane sinks (nlevsoi)
   real(r8), pointer :: ch4stress_sat(:,:) ! Ratio of methane available to the total per-timestep methane sinks (nlevsoi)
   real(r8), pointer :: qflx_surf_lag(:)	! time-lagged surface runoff (mm H2O /s)
   real(r8), pointer :: finundated_lag(:)       ! time-lagged fractional inundated area
   real(r8), pointer :: layer_sat_lag(:,:) ! Lagged saturation status of soil layer in the unsaturated zone (1 = sat)


end type column_ch4_type
#endif


!----------------------------------------------------
! column nitrogen state variables structure
!----------------------------------------------------
type, public :: column_nstate_type
   type(pft_nstate_type):: pns_a              !pft-level nitrogen state variables averaged to the column
   ! all n pools involved in decomposition
   real(r8), pointer :: decomp_npools_vr(:,:,:)    ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   real(r8), pointer :: sminn_vr(:,:)              ! (gN/m3) vertically-resolved soil mineral N
   real(r8), pointer :: col_ntrunc_vr(:,:)         ! (gN/m3) vertically-resolved column-level sink for N truncation
#ifdef NITRIF_DENITRIF
   real(r8), pointer :: smin_no3_vr(:,:)           ! (gN/m3) vertically-resolved soil mineral NO3
   real(r8), pointer :: smin_no3(:)                ! (gN/m2) soil mineral NO3 pool
   real(r8), pointer :: smin_nh4_vr(:,:)           ! (gN/m3) vertically-resolved soil mineral NH4
   real(r8), pointer :: smin_nh4(:)                ! (gN/m2) soil mineral NH4 pool
#endif
   ! wood product pools, for dynamic landcover
   real(r8), pointer :: seedn(:)                   ! (gN/m2) column-level pool for seeding new PFTs
   real(r8), pointer :: prod10n(:)                 ! (gN/m2) wood product N pool, 10-year lifespan
   real(r8), pointer :: prod100n(:)                ! (gN/m2) wood product N pool, 100-year lifespan
   real(r8), pointer :: totprodn(:)                ! (gN/m2) total wood product N
   ! summary (diagnostic) state variables, not involved in mass balance
   real(r8), pointer :: decomp_npools(:,:)         ! (gN/m2)  decomposing (litter, cwd, soil) N pools
   real(r8), pointer :: decomp_npools_1m(:,:)           ! (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
   real(r8), pointer :: sminn(:)                   ! (gN/m2) soil mineral N
   real(r8), pointer :: col_ntrunc(:)              ! (gN/m2) column-level sink for N truncation
   real(r8), pointer :: cwdn(:)                    ! (gN/m2) Diagnostic: coarse woody debris N
   real(r8), pointer :: totlitn(:)                 ! (gN/m2) total litter nitrogen
   real(r8), pointer :: totsomn(:)                 ! (gN/m2) total soil organic matter nitrogen
   real(r8), pointer :: totlitn_1m(:)              ! (gN/m2) total litter nitrogen to 1 meter
   real(r8), pointer :: totsomn_1m(:)              ! (gN/m2) total soil organic matter nitrogen to 1 meter
   real(r8), pointer :: totecosysn(:)              ! (gN/m2) total ecosystem nitrogen, incl veg 
   real(r8), pointer :: totcoln(:)                 ! (gN/m2) total column nitrogen, incl veg
end type column_nstate_type

!----------------------------------------------------
! column VOC state variables structure
!----------------------------------------------------
type, public :: column_vstate_type
   type(pft_vstate_type):: pvs_a              !pft-level VOC state variables averaged to the column
end type column_vstate_type

#if (defined CNDV)
!----------------------------------------------------
! column DGVM state variables structure
!----------------------------------------------------
type, public :: column_dgvstate_type
   type(pft_dgvstate_type):: pdgvs_a
end type column_dgvstate_type
#endif

!----------------------------------------------------
! column dust state variables structure
!----------------------------------------------------
type, public :: column_dstate_type
   real(r8), pointer :: dummy_entry(:)
end type column_dstate_type

!----------------------------------------------------
! column energy flux variables structure
!----------------------------------------------------
type, public :: column_eflux_type
   type(pft_eflux_type):: pef_a	              ! pft-level energy flux variables averaged to the column
   real(r8), pointer :: eflx_snomelt(:)       ! snow melt heat flux (W/m**2)
   real(r8), pointer :: eflx_snomelt_u(:)     ! urban snow melt heat flux (W/m**2)
   real(r8), pointer :: eflx_snomelt_r(:)     ! rural snow melt heat flux (W/m**2)
   real(r8), pointer :: eflx_impsoil(:)	      ! implicit evaporation for soil temperature equation
   real(r8), pointer :: eflx_fgr12(:)         ! ground heat flux between soil layers 1 and 2 (W/m2)
   real(r8), pointer :: eflx_fgr(:,:)         ! (rural) soil downward heat flux (W/m2) (1:nlevgrnd)
   ! Urban variable
   real(r8), pointer :: eflx_building_heat(:) ! heat flux from urban building interior to urban walls, roof (W/m**2)
   real(r8), pointer :: eflx_urban_ac(:)      ! urban air conditioning flux (W/m**2)
   real(r8), pointer :: eflx_urban_heat(:)    ! urban heating flux (W/m**2)
   real(r8), pointer :: eflx_bot(:)           ! heat flux from beneath the soil or ice column (W/m**2)
                                              ! positive upward; usually eflx_bot >= 0
end type column_eflux_type

!----------------------------------------------------
! column momentum flux variables structure
!----------------------------------------------------
type, public :: column_mflux_type
   type(pft_mflux_type)::  pmf_a        ! pft-level momentum flux variables averaged to the column
end type column_mflux_type

!----------------------------------------------------
! column water flux variables structure
!----------------------------------------------------
type, public :: column_wflux_type
   type(pft_wflux_type):: pwf_a	        ! pft-level water flux variables averaged to the column
   real(r8), pointer :: qflx_infl(:)	! infiltration (mm H2O /s)
   real(r8), pointer :: qflx_surf(:)	! surface runoff (mm H2O /s)
   real(r8), pointer :: qflx_drain(:) 	! sub-surface runoff (mm H2O /s)
   real(r8), pointer :: qflx_top_soil(:)! net water input into soil from top (mm/s)
   real(r8), pointer :: qflx_h2osfc_to_ice(:) ! conversion of h2osfc to ice
   real(r8), pointer :: qflx_h2osfc_surf(:)   !surface water runoff
   real(r8), pointer :: qflx_snow_h2osfc(:)   !snow falling on surface water
   real(r8), pointer :: qflx_drain_perched(:) ! sub-surface runoff from perched wt (mm H2O /s)
   real(r8), pointer :: qflx_floodc(:) 	      ! flood water flux at column level
   real(r8), pointer :: qflx_sl_top_soil(:)   ! liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
   real(r8), pointer :: qflx_snomelt(:)       ! snow melt (mm H2O /s)
   real(r8), pointer :: qflx_snow_melt(:)     ! snow melt (net)
   real(r8), pointer :: qflx_qrgwl(:) 	    ! qflx_surf at glaciers, wetlands, lakes
   real(r8), pointer :: qflx_runoff(:) 	    ! total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
   real(r8), pointer :: qflx_runoff_u(:)    ! Urban total runoff (qflx_drain+qflx_surf) (mm H2O /s)
   real(r8), pointer :: qflx_runoff_r(:)    ! Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
   real(r8), pointer :: qmelt(:) 	    ! snow melt [mm/s]
   real(r8), pointer :: h2ocan_loss(:)      ! mass balance correction term for dynamic weights
   real(r8), pointer :: qflx_rsub_sat(:)    ! soil saturation excess [mm/s]
   real(r8), pointer :: flx_bc_dep_dry(:)   ! dry (BCPHO+BCPHI) BC deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_bc_dep_wet(:)   ! wet (BCPHI) BC deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_bc_dep_pho(:)   ! hydrophobic BC deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_bc_dep_phi(:)   ! hydrophillic BC deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_bc_dep(:)       ! total (dry+wet) BC deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_oc_dep_dry(:)   ! dry (OCPHO+OCPHI) OC deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_oc_dep_wet(:)   ! wet (OCPHI) OC deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_oc_dep_pho(:)   ! hydrophobic OC deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_oc_dep_phi(:)   ! hydrophillic OC deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_oc_dep(:)       ! total (dry+wet) OC deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_dst_dep_dry1(:) ! dust species 1 dry deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_dst_dep_wet1(:) ! dust species 1 wet deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_dst_dep_dry2(:) ! dust species 2 dry deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_dst_dep_wet2(:) ! dust species 2 wet deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_dst_dep_dry3(:) ! dust species 3 dry deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_dst_dep_wet3(:) ! dust species 3 wet deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_dst_dep_dry4(:) ! dust species 4 dry deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_dst_dep_wet4(:) ! dust species 4 wet deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: flx_dst_dep(:)      ! total (dry+wet) dust deposition on ground (positive definite) (col) [kg/s]
   real(r8), pointer :: qflx_snofrz_lyr(:,:)! snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
   real(r8), pointer :: qflx_snofrz_col(:)  ! column-integrated snow freezing rate (positive definite) (col) [kg m-2 s-1]
   real(r8), pointer :: qflx_irrig(:)     !irrigation flux (mm H2O/s)
   real(r8), pointer :: qflx_glcice(:)      ! net flux of new glacial ice (growth - melt) (mm H2O/s), passed to GLC
   real(r8), pointer :: qflx_glcice_frz(:)  ! ice growth (positive definite) (mm H2O/s)
   real(r8), pointer :: qflx_glcice_melt(:) ! ice melt (positive definite) (mm H2O/s)
   real(r8), pointer :: glc_rofi(:)         ! ice runoff passed from GLC to CLM (mm H2O /s)
   real(r8), pointer :: glc_rofl(:)         ! liquid runoff passed from GLC to CLM (mm H2O /s)
end type column_wflux_type

!----------------------------------------------------
! column carbon flux variables structure
!----------------------------------------------------
type, public :: column_cflux_type
   type(pft_cflux_type):: pcf_a                           ! pft-level carbon flux variables averaged to the column
   ! phenology: litterfall and crop fluxes
   real(r8), pointer :: phenology_c_to_litr_met_c(:,:)             ! C fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gC/m3/s)
   real(r8), pointer :: phenology_c_to_litr_cel_c(:,:)             ! C fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gC/m3/s)
   real(r8), pointer :: phenology_c_to_litr_lig_c(:,:)             ! C fluxes associated with phenology (litterfall and crop) to litter lignin pool (gC/m3/s)
   ! gap mortality
   real(r8), pointer :: gap_mortality_c_to_litr_met_c(:,:)         ! C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
   real(r8), pointer :: gap_mortality_c_to_litr_cel_c(:,:)         ! C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
   real(r8), pointer :: gap_mortality_c_to_litr_lig_c(:,:)         ! C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
   real(r8), pointer :: gap_mortality_c_to_cwdc(:,:)               ! C fluxes associated with gap mortality to CWD pool (gC/m3/s)
   ! fire
   real(r8), pointer :: fire_mortality_c_to_cwdc(:,:)              ! C fluxes associated with fire mortality to CWD pool (gC/m3/s)
   ! harvest
   real(r8), pointer :: harvest_c_to_litr_met_c(:,:)               ! C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
   real(r8), pointer :: harvest_c_to_litr_cel_c(:,:)               ! C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
   real(r8), pointer :: harvest_c_to_litr_lig_c(:,:)               ! C fluxes associated with harvest to litter lignin pool (gC/m3/s)
   real(r8), pointer :: harvest_c_to_cwdc(:,:)                     ! C fluxes associated with harvest to CWD pool (gC/m3/s)
   ! new variables for CN code
   real(r8), pointer :: hrv_deadstemc_to_prod10c(:)                ! dead stem C harvest mortality to 10-year product pool (gC/m2/s)        
   real(r8), pointer :: hrv_deadstemc_to_prod100c(:)               ! dead stem C harvest mortality to 100-year product pool (gC/m2/s)        
   ! column-level fire fluxes
   real(r8), pointer :: m_decomp_cpools_to_fire_vr(:,:,:)          ! vertically-resolved decomposing C fire loss (gC/m3/s)
   real(r8), pointer :: m_decomp_cpools_to_fire(:,:)               ! vertically-integrated (diagnostic) decomposing C fire loss (gC/m2/s)
   real(r8), pointer :: m_c_to_litr_met_fire(:,:)        ! C from leaf, froot, xfer and storage C to litter labile C by fire (gC/m3/s) 
   real(r8), pointer :: m_c_to_litr_cel_fire(:,:)        ! C from leaf, froot, xfer and storage C to litter cellulose C by fire (gC/m3/s) 
   real(r8), pointer :: m_c_to_litr_lig_fire(:,:)        ! C from leaf, froot, xfer and storage C to litter lignin C by fire (gC/m3/s) 
   real(r8), pointer :: lf_conv_cflux(:)                ! (gC/m2/s) conversion C flux due to BET and BDT area decreasing (immediate loss to atm)
   real(r8), pointer :: somc_fire(:)                     ! (gC/m2/s) carbon emissions due to peat burning

   ! decomposition fluxes
   real(r8), pointer :: decomp_cascade_hr_vr(:,:,:)            ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
   real(r8), pointer :: decomp_cascade_hr(:,:)                 ! vertically-integrated (diagnostic) het. resp. from decomposing C pools (gC/m2/s)
   real(r8), pointer :: decomp_cascade_ctransfer_vr(:,:,:)     ! vertically-resolved C transferred along deomposition cascade (gC/m3/s)
   real(r8), pointer :: decomp_cascade_ctransfer(:,:)          ! vertically-integrated (diagnostic) C transferred along deomposition cascade (gC/m2/s)
   real(r8), pointer :: decomp_cpools_sourcesink(:,:,:)        ! (gC/m3/timestep)  change in decomposing c pools.  Used to update concentrations concurrently with vertical transport
   real(r8), pointer :: decomp_k(:,:,:)                        ! rate constant for decomposition (1./sec)
   real(r8), pointer :: hr_vr(:,:)                             ! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
   real(r8), pointer :: o_scalar(:,:)                          ! fraction by which decomposition is limited by anoxia
   real(r8), pointer :: w_scalar(:,:)                          ! fraction by which decomposition is limited by moisture availability
   real(r8), pointer :: t_scalar(:,:)                          ! fraction by which decomposition is limited by temperature
   real(r8), pointer :: som_c_leached(:)                           ! total SOM C loss from vertical transport (gC/m^2/s)
   real(r8), pointer :: decomp_cpools_leached(:,:)                 ! C loss from vertical transport from each decomposing C pool (gC/m^2/s)
   real(r8), pointer :: decomp_cpools_transport_tendency(:,:,:)    ! C tendency due to vertical transport in decomposing C pools (gC/m^3/s)
#ifdef NITRIF_DENITRIF
   real(r8), pointer :: phr_vr(:,:)                            ! potential hr (not N-limited) (gC/m3/s)
#endif
   ! dynamic landcover fluxes
#ifdef CN
   real(r8), pointer :: dwt_seedc_to_leaf(:)      ! (gC/m2/s) seed source to PFT-level
   real(r8), pointer :: dwt_seedc_to_deadstem(:)  ! (gC/m2/s) seed source to PFT-level
   real(r8), pointer :: dwt_conv_cflux(:)         ! (gC/m2/s) conversion C flux (immediate loss to atm)
   real(r8), pointer :: dwt_prod10c_gain(:)       ! (gC/m2/s) addition to 10-yr wood product pool
   real(r8), pointer :: dwt_prod100c_gain(:)      ! (gC/m2/s) addition to 100-yr wood product pool
   real(r8), pointer :: dwt_frootc_to_litr_met_c(:,:)   ! (gC/m3/s) fine root to litter due to landcover change
   real(r8), pointer :: dwt_frootc_to_litr_cel_c(:,:)   ! (gC/m3/s) fine root to litter due to landcover change
   real(r8), pointer :: dwt_frootc_to_litr_lig_c(:,:)   ! (gC/m3/s) fine root to litter due to landcover change
   real(r8), pointer :: dwt_livecrootc_to_cwdc(:,:) ! (gC/m3/s) live coarse root to CWD due to landcover change
   real(r8), pointer :: dwt_deadcrootc_to_cwdc(:,:) ! (gC/m3/s) dead coarse root to CWD due to landcover change
   real(r8), pointer :: dwt_closs(:)              ! (gC/m2/s) total carbon loss from product pools and conversion
   real(r8), pointer :: landuseflux(:)            ! (gC/m2/s) dwt_closs+product_closs
   real(r8), pointer :: landuptake(:)             ! (gC/m2/s) nee-landuseflux
   ! wood product pool loss fluxes
   real(r8), pointer :: prod10c_loss(:)           ! (gC/m2/s) decomposition loss from 10-yr wood product pool
   real(r8), pointer :: prod100c_loss(:)          ! (gC/m2/s) decomposition loss from 100-yr wood product pool
   real(r8), pointer :: product_closs(:)          ! (gC/m2/s) total wood product carbon loss
#endif
   ! summary (diagnostic) flux variables, not involved in mass balance
   real(r8), pointer :: lithr(:)         ! (gC/m2/s) litter heterotrophic respiration 
   real(r8), pointer :: somhr(:)         ! (gC/m2/s) soil organic matter heterotrophic respiration
   real(r8), pointer :: hr(:)            ! (gC/m2/s) total heterotrophic respiration
   real(r8), pointer :: sr(:)            ! (gC/m2/s) total soil respiration (HR + root resp)
   real(r8), pointer :: er(:)            ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
   real(r8), pointer :: litfire(:)       ! (gC/m2/s) litter fire losses
   real(r8), pointer :: somfire(:)       ! (gC/m2/s) soil organic matter fire losses
   real(r8), pointer :: totfire(:)       ! (gC/m2/s) total ecosystem fire losses
   real(r8), pointer :: nep(:)           ! (gC/m2/s) net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink
   real(r8), pointer :: nbp(:)           ! (gC/m2/s) net biome production, includes fire, landuse, and harvest flux, positive for sink
   real(r8), pointer :: nee(:)           ! (gC/m2/s) net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source
   real(r8), pointer :: col_cinputs(:)   ! (gC/m2/s) total column-level carbon inputs (for balance check)
   real(r8), pointer :: col_coutputs(:)  ! (gC/m2/s) total column-level carbon outputs (for balance check) 

#if (defined CN)
   ! CLAMP summary (diagnostic) flux variables, not involved in mass balance
   real(r8), pointer :: cwdc_hr(:)       ! (gC/m2/s) col-level coarse woody debris C heterotrophic respiration
   real(r8), pointer :: cwdc_loss(:)     ! (gC/m2/s) col-level coarse woody debris C loss
   real(r8), pointer :: litterc_loss(:)  ! (gC/m2/s) col-level litter C loss
#endif

   ! new variables for fire
   real(r8), pointer :: col_fire_closs(:) ! (gC/m2/s) total column-level fire C loss
end type column_cflux_type

!----------------------------------------------------
! column nitrogen flux variables structure
!----------------------------------------------------
type, public :: column_nflux_type
   type(pft_nflux_type):: pnf_a        !pft-level nitrogen flux variables averaged to the column
   ! new variables for CN code
   ! deposition fluxes
   real(r8), pointer :: ndep_to_sminn(:)                   ! atmospheric N deposition to soil mineral N (gN/m2/s)
   real(r8), pointer :: nfix_to_sminn(:)                   ! symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s) 
   real(r8), pointer :: fert_to_sminn(:)                   ! fertilizer N to soil mineral N (gN/m2/s)
   real(r8), pointer :: soyfixn_to_sminn(:)                ! soybean fixation to soil mineral N (gN/m2/s)
   ! phenology: litterfall and crop fluxes
   real(r8), pointer :: phenology_n_to_litr_met_n(:,:)             ! N fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gN/m3/s)
   real(r8), pointer :: phenology_n_to_litr_cel_n(:,:)             ! N fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gN/m3/s)
   real(r8), pointer :: phenology_n_to_litr_lig_n(:,:)             ! N fluxes associated with phenology (litterfall and crop) to litter lignin pool (gN/m3/s)
   ! gap mortality
   real(r8), pointer :: gap_mortality_n_to_litr_met_n(:,:)         ! N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
   real(r8), pointer :: gap_mortality_n_to_litr_cel_n(:,:)         ! N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
   real(r8), pointer :: gap_mortality_n_to_litr_lig_n(:,:)         ! N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
   real(r8), pointer :: gap_mortality_n_to_cwdn(:,:)               ! N fluxes associated with gap mortality to CWD pool (gN/m3/s)
   ! fire
   real(r8), pointer :: fire_mortality_n_to_cwdn(:,:)              ! N fluxes associated with fire mortality to CWD pool (gN/m3/s)
   ! harvest
   real(r8), pointer :: harvest_n_to_litr_met_n(:,:)               ! N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
   real(r8), pointer :: harvest_n_to_litr_cel_n(:,:)               ! N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
   real(r8), pointer :: harvest_n_to_litr_lig_n(:,:)               ! N fluxes associated with harvest to litter lignin pool (gN/m3/s)
   real(r8), pointer :: harvest_n_to_cwdn(:,:)                     ! N fluxes associated with harvest to CWD pool (gN/m3/s)
   !
   real(r8), pointer :: hrv_deadstemn_to_prod10n(:)                ! dead stem N harvest mortality to 10-year product pool (gN/m2/s)
   real(r8), pointer :: hrv_deadstemn_to_prod100n(:)               ! dead stem N harvest mortality to 100-year product pool (gN/m2/s)
   real(r8), pointer :: m_decomp_npools_to_fire_vr(:,:,:)          ! vertically-resolved decomposing N fire loss (gN/m3/s)
   real(r8), pointer :: m_decomp_npools_to_fire(:,:)               ! vertically-integrated (diagnostic) decomposing N fire loss (gN/m2/s)
   ! column-level fire N fluxes added by F. Li and S. Levis
   real(r8), pointer :: m_n_to_litr_met_fire(:,:)        ! N from leaf, froot, xfer and storage N to litter labile N by fire (gN/m3/s) 
   real(r8), pointer :: m_n_to_litr_cel_fire(:,:)        ! N from leaf, froot, xfer and storage N to litter cellulose N by fire (gN/m3/s) 
   real(r8), pointer :: m_n_to_litr_lig_fire(:,:)        ! N from leaf, froot, xfer and storage N to litter lignin N by fire (gN/m3/s) 

   ! decomposition fluxes
   real(r8), pointer :: decomp_cascade_ntransfer_vr(:,:,:)    ! vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
   real(r8), pointer :: decomp_cascade_ntransfer(:,:)         ! vert-int (diagnostic) transfer of N from donor to receiver pool along decomp. cascade (gN/m2/s)
   real(r8), pointer :: decomp_cascade_sminn_flux_vr(:,:,:)   ! vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
   real(r8), pointer :: decomp_cascade_sminn_flux(:,:)        ! vert-int (diagnostic) mineral N flux for transition along decomposition cascade (gN/m2/s)
   real(r8), pointer :: decomp_npools_sourcesink(:,:,:)       ! (gN/m3)  change in decomposing n pools (sum of all additions and subtractions from stateupdate1).  Used to update concentrations concurrently with vertical transport
   ! vertically-resolved immobilization fluxes
   real(r8), pointer :: potential_immob_vr(:,:)       ! vertically-resolved potential N immobilization (gN/m3/s) at each level
   real(r8), pointer :: potential_immob(:)            ! vert-int (diagnostic) potential N immobilization (gN/m2/s)
   real(r8), pointer :: actual_immob_vr(:,:)          ! vertically-resolved actual N immobilization (gN/m3/s) at each level
   real(r8), pointer :: actual_immob(:)               ! vert-int (diagnostic) actual N immobilization (gN/m2/s)
   real(r8), pointer :: sminn_to_plant_vr(:,:)        ! vertically-resolved plant uptake of soil mineral N (gN/m3/s)
   real(r8), pointer :: sminn_to_plant(:)             ! vert-int (diagnostic) plant uptake of soil mineral N (gN/m2/s)
   real(r8), pointer :: supplement_to_sminn_vr(:,:)   ! vertically-resolved supplemental N supply (gN/m3/s)
   real(r8), pointer :: supplement_to_sminn(:)        ! vert-int (diagnostic) supplemental N supply (gN/m2/s)
   real(r8), pointer :: gross_nmin_vr(:,:)            ! vertically-resolved gross rate of N mineralization (gN/m3/s)
   real(r8), pointer :: gross_nmin(:)                 ! vert-int (diagnostic) gross rate of N mineralization (gN/m2/s)
   real(r8), pointer :: net_nmin_vr(:,:)              ! vertically-resolved net rate of N mineralization (gN/m3/s)
   real(r8), pointer :: net_nmin(:)                   ! vert-int (diagnostic) net rate of N mineralization (gN/m2/s)
#ifdef NITRIF_DENITRIF
   ! nitrification / denitrification fluxes
   real(r8), pointer :: f_nit_vr(:,:)              ! (gN/m3/s) soil nitrification flux
   real(r8), pointer :: f_denit_vr(:,:)            ! (gN/m3/s) soil denitrification flux
   real(r8), pointer :: f_nit(:)                   ! (gN/m2/s) soil nitrification flux
   real(r8), pointer :: f_denit(:)                 ! (gN/m2/s) soil denitrification flux

   real(r8), pointer :: pot_f_nit_vr(:,:)          ! (gN/m3/s) potential soil nitrification flux
   real(r8), pointer :: pot_f_denit_vr(:,:)        ! (gN/m3/s) potential soil denitrification flux
   real(r8), pointer :: pot_f_nit(:)               ! (gN/m2/s) potential soil nitrification flux
   real(r8), pointer :: pot_f_denit(:)             ! (gN/m2/s) potential soil denitrification flux
   real(r8), pointer :: n2_n2o_ratio_denit_vr(:,:)    ! ratio of N2 to N2O production by denitrification [gN/gN]
   real(r8), pointer :: f_n2o_denit_vr(:,:)        ! flux of N2o from denitrification [gN/m^3/s]
   real(r8), pointer :: f_n2o_denit(:)             ! flux of N2o from denitrification [gN/m^2/s]
   real(r8), pointer :: f_n2o_nit_vr(:,:)          ! flux of N2o from nitrification [gN/m^3/s]
   real(r8), pointer :: f_n2o_nit(:)               ! flux of N2o from nitrification [gN/m^2/s]

   ! immobilization / uptake fluxes
   real(r8), pointer :: actual_immob_no3_vr(:,:)   ! vertically-resolved actual immobilization of NO3 (gN/m3/s)
   real(r8), pointer :: actual_immob_nh4_vr(:,:)   ! vertically-resolved actual immobilization of NH4 (gN/m3/s)
   real(r8), pointer :: smin_no3_to_plant_vr(:,:)  ! vertically-resolved plant uptake of soil NO3 (gN/m3/s)
   real(r8), pointer :: smin_nh4_to_plant_vr(:,:)  ! vertically-resolved plant uptake of soil NH4 (gN/m3/s)
   real(r8), pointer :: actual_immob_no3(:)        ! actual immobilization of NO3 (gN/m2/s)
   real(r8), pointer :: actual_immob_nh4(:)        ! actual immobilization of NH4 (gN/m2/s)
   real(r8), pointer :: smin_no3_to_plant(:)       ! plant uptake of soil NO3 (gN/m2/s)
   real(r8), pointer :: smin_nh4_to_plant(:)       ! plant uptake of soil Nh4 (gN/m2/s)
   ! leaching fluxes
   real(r8), pointer :: smin_no3_leached_vr(:,:)    ! vertically-resolved soil mineral NO3 loss to leaching (gN/m3/s)
   real(r8), pointer :: smin_no3_leached(:)         ! soil mineral NO3 pool loss to leaching (gN/m2/s)
   real(r8), pointer :: smin_no3_runoff_vr(:,:)     ! vertically-resolved rate of mineral NO3 loss with runoff (gN/m3/s)
   real(r8), pointer :: smin_no3_runoff(:)          ! soil mineral NO3 pool loss to runoff (gN/m2/s)

   ! NITRIF_DENITRIF diagnostic quantities
   real(r8), pointer :: smin_no3_massdens_vr(:,:)      ! (ugN / g soil) soil nitrate concentration
   real(r8), pointer :: soil_bulkdensity(:,:)          ! (kg soil / m3) bulk density of soil
   real(r8), pointer :: k_nitr_t_vr(:,:)
   real(r8), pointer :: k_nitr_ph_vr(:,:)
   real(r8), pointer :: k_nitr_h2o_vr(:,:)
   real(r8), pointer :: k_nitr_vr(:,:)
   real(r8), pointer :: wfps_vr(:,:)
   real(r8), pointer :: fmax_denit_carbonsubstrate_vr(:,:)
   real(r8), pointer :: fmax_denit_nitrate_vr(:,:)
   real(r8), pointer :: f_denit_base_vr(:,:)   ! nitrification and denitrification fluxes
   real(r8), pointer :: diffus(:,:) !diffusivity (m2/s)
   real(r8), pointer :: ratio_k1(:,:)
   real(r8), pointer :: ratio_no3_co2(:,:)
   real(r8), pointer :: soil_co2_prod(:,:)
   real(r8), pointer :: fr_WFPS(:,:)

   real(r8), pointer :: r_psi(:,:)
   real(r8), pointer :: anaerobic_frac(:,:)
#else
   ! denitrification fluxes
   real(r8), pointer :: sminn_to_denit_decomp_cascade_vr(:,:,:)   ! vertically-resolved denitrification along decomp cascade (gN/m3/s) 
   real(r8), pointer :: sminn_to_denit_decomp_cascade(:,:)        ! vertically-integrated (diagnostic) denitrification along decomp cascade (gN/m2/s) 
   real(r8), pointer :: sminn_to_denit_excess_vr(:,:)             ! vertically-resolved denitrification from excess mineral N pool (gN/m3/s)
   real(r8), pointer :: sminn_to_denit_excess(:)                  ! vertically-integrated (diagnostic) denitrification from excess mineral N pool (gN/m2/s)
   ! leaching fluxes
   real(r8), pointer :: sminn_leached_vr(:,:)                     ! vertically-resolved soil mineral N pool loss to leaching (gN/m3/s)
   real(r8), pointer :: sminn_leached(:)                          ! soil mineral N pool loss to leaching (gN/m2/s)
#endif
   ! dynamic landcover fluxes
   real(r8), pointer :: dwt_seedn_to_leaf(:)            ! (gN/m2/s) seed source to PFT-level
   real(r8), pointer :: dwt_seedn_to_deadstem(:)        ! (gN/m2/s) seed source to PFT-level
   real(r8), pointer :: dwt_conv_nflux(:)               ! (gN/m2/s) conversion N flux (immediate loss to atm)
   real(r8), pointer :: dwt_prod10n_gain(:)             ! (gN/m2/s) addition to 10-yr wood product pool
   real(r8), pointer :: dwt_prod100n_gain(:)            ! (gN/m2/s) addition to 100-yr wood product pool
   real(r8), pointer :: dwt_frootn_to_litr_met_n(:,:)   ! (gN/m3/s) fine root to litter due to landcover change
   real(r8), pointer :: dwt_frootn_to_litr_cel_n(:,:)   ! (gN/m3/s) fine root to litter due to landcover change
   real(r8), pointer :: dwt_frootn_to_litr_lig_n(:,:)   ! (gN/m3/s) fine root to litter due to landcover change
   real(r8), pointer :: dwt_livecrootn_to_cwdn(:,:)     ! (gN/m3/s) live coarse root to CWD due to landcover change
   real(r8), pointer :: dwt_deadcrootn_to_cwdn(:,:)     ! (gN/m3/s) dead coarse root to CWD due to landcover change
   real(r8), pointer :: dwt_nloss(:)                    ! (gN/m2/s) total nitrogen loss from product pools and conversion
   ! wood product pool loss fluxes
   real(r8), pointer :: prod10n_loss(:)                 ! (gN/m2/s) decomposition loss from 10-yr wood product pool
   real(r8), pointer :: prod100n_loss(:)                ! (gN/m2/s) decomposition loss from 100-yr wood product pool
   real(r8), pointer :: product_nloss(:)                ! (gN/m2/s) total wood product nitrogen loss
   ! summary (diagnostic) flux variables, not involved in mass balance
   real(r8), pointer :: denit(:)                        ! total rate of denitrification (gN/m2/s)
   real(r8), pointer :: col_ninputs(:)                  ! column-level N inputs (gN/m2/s)
   real(r8), pointer :: col_noutputs(:)                 ! column-level N outputs (gN/m2/s)
   ! new variables for fire
   real(r8), pointer :: col_fire_nloss(:)               ! total column-level fire N loss (gN/m2/s)

   real(r8), pointer :: som_n_leached(:)                           ! total SOM N loss from vertical transport (gN/m^2/s)
   real(r8), pointer :: decomp_npools_leached(:,:)                 ! N loss from vertical transport from each decomposing N pool (gN/m^2/s)
   real(r8), pointer :: decomp_npools_transport_tendency(:,:,:)    ! N tendency due to vertical transport in decomposing N pools (gN/m^3/s)
end type column_nflux_type

!----------------------------------------------------
! column dust flux variables structure
!----------------------------------------------------
type, public :: column_dflux_type
   type(pft_dflux_type):: pdf_a         !pft-level dust flux variables averaged to the column
end type column_dflux_type

!----------------------------------------------------
! End definition of structures defined at the column_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the landunit_type level
!----------------------------------------------------
! landunit physical state variables structure
! note - landunit type can be vegetated (includes bare soil), deep lake, 
! shallow lake, wetland, glacier or urban
!----------------------------------------------------
type, public :: landunit_pstate_type
   type(column_pstate_type):: cps_a            !column-level physical state variables averaged to landunit
   ! Urban variables
   real(r8), pointer :: t_building(:)         ! internal building temperature (K)
   real(r8), pointer :: t_building_max(:)     ! maximum internal building temperature (K)
   real(r8), pointer :: t_building_min(:)     ! minimum internal building temperature (K)
   real(r8), pointer :: tk_wall(:,:)          ! thermal conductivity of urban wall (W/m/K)
   real(r8), pointer :: tk_roof(:,:)          ! thermal conductivity of urban roof (W/m/K)
   real(r8), pointer :: tk_improad(:,:)       ! thermal conductivity of urban impervious road (W/m/K)
   real(r8), pointer :: cv_wall(:,:)          ! heat capacity of urban wall (J/m^3/K)
   real(r8), pointer :: cv_roof(:,:)          ! heat capacity of urban roof (J/m^3/K)
   real(r8), pointer :: cv_improad(:,:)       ! heat capacity of urban impervious road (J/m^3/K)
   real(r8), pointer :: thick_wall(:)         ! total thickness of urban wall (m)
   real(r8), pointer :: thick_roof(:)         ! total thickness of urban roof (m)
   integer, pointer :: nlev_improad(:)        ! number of impervious road layers (-)
   real(r8), pointer :: vf_sr(:)              ! view factor of sky for road
   real(r8), pointer :: vf_wr(:)              ! view factor of one wall for road
   real(r8), pointer :: vf_sw(:)              ! view factor of sky for one wall
   real(r8), pointer :: vf_rw(:)              ! view factor of road for one wall
   real(r8), pointer :: vf_ww(:)              ! view factor of opposing wall for one wall
   real(r8), pointer :: taf(:)                ! urban canopy air temperature (K)
   real(r8), pointer :: qaf(:)                ! urban canopy air specific humidity (kg/kg)
   real(r8), pointer :: sabs_roof_dir(:,:)       ! direct solar absorbed by roof per unit ground area per unit incident flux
   real(r8), pointer :: sabs_roof_dif(:,:)       ! diffuse solar absorbed by roof per unit ground area per unit incident flux
   real(r8), pointer :: sabs_sunwall_dir(:,:)    ! direct  solar absorbed by sunwall per unit wall area per unit incident flux
   real(r8), pointer :: sabs_sunwall_dif(:,:)    ! diffuse solar absorbed by sunwall per unit wall area per unit incident flux
   real(r8), pointer :: sabs_shadewall_dir(:,:)  ! direct  solar absorbed by shadewall per unit wall area per unit incident flux
   real(r8), pointer :: sabs_shadewall_dif(:,:)  ! diffuse solar absorbed by shadewall per unit wall area per unit incident flux
   real(r8), pointer :: sabs_improad_dir(:,:)    ! direct  solar absorbed by impervious road per unit ground area per unit incident flux
   real(r8), pointer :: sabs_improad_dif(:,:)    ! diffuse solar absorbed by impervious road per unit ground area per unit incident flux
   real(r8), pointer :: sabs_perroad_dir(:,:)    ! direct  solar absorbed by pervious road per unit ground area per unit incident flux
   real(r8), pointer :: sabs_perroad_dif(:,:)    ! diffuse solar absorbed by pervious road per unit ground area per unit incident flux
end type landunit_pstate_type

!----------------------------------------------------
! landunit energy flux variables structure
!----------------------------------------------------
type, public :: landunit_eflux_type
   type(column_eflux_type)::	cef_a		! column-level energy flux variables averaged to landunit
   ! Urban variables
   real(r8), pointer :: eflx_traffic_factor(:)  ! multiplicative traffic factor for sensible heat flux from urban traffic (-)
   real(r8), pointer :: eflx_traffic(:)         ! traffic sensible heat flux (W/m**2)
   real(r8), pointer :: eflx_wasteheat(:)       ! sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
   real(r8), pointer :: eflx_heat_from_ac(:)    ! sensible heat flux to be put back into canyon due to removal by AC (W/m**2)
end type landunit_eflux_type

!----------------------------------------------------
! End definition of structures defined at the landunit_type level
!----------------------------------------------------
!*******************************************************************************

!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the gridcell_type level
!----------------------------------------------------
! gridcell physical state variables structure
!----------------------------------------------------
type, public :: gridcell_pstate_type
   type(column_pstate_type):: cps_a   !column-level physical state variables averaged to gridcell
end type gridcell_pstate_type

!----------------------------------------------------
! gridcell energy state variables structure
!----------------------------------------------------
type, public :: gridcell_estate_type
   type(column_estate_type):: ces_a            !column-level energy state variables averaged to gridcell
   real(r8), pointer :: gc_heat1(:)            ! initial gridcell total heat content
   real(r8), pointer :: gc_heat2(:)            ! post land cover change total heat content
end type gridcell_estate_type

!----------------------------------------------------
! gridcell water state variables structure
!----------------------------------------------------
type, public :: gridcell_wstate_type
   type(column_wstate_type):: cws_a            !column-level water state variables averaged to gridcell
   real(r8), pointer :: gc_liq1(:)             ! initial gridcell total h2o liq content
   real(r8), pointer :: gc_liq2(:)             ! post land cover change total liq content
   real(r8), pointer :: gc_ice1(:)             ! initial gridcell total h2o liq content
   real(r8), pointer :: gc_ice2(:)             ! post land cover change total ice content
end type gridcell_wstate_type

!----------------------------------------------------
! gridcell carbon state variables structure
!----------------------------------------------------
type, public :: gridcell_cstate_type
   type(column_cstate_type):: ccs_a            !column-level carbon state variables averaged to gridcell
end type gridcell_cstate_type

!----------------------------------------------------
! gridcell CH4 flux variables structure
!----------------------------------------------------
#ifdef LCH4
type, public :: gridcell_ch4_type
   real(r8), pointer :: c_atm(:,:)            ! Atmospheric conc of CH4, O2, CO2 (mol/m3)
   real(r8), pointer :: ch4co2f(:)     ! gridcell CO2 production from CH4 oxidation (g C/m**2/s)
   real(r8), pointer :: ch4prodg(:)    !gridcell average CH4 production (g C/m^2/s)
   real(r8), pointer :: nem(:)    !gridcell average net methane correction to CO2 flux (g C/m^2/s)
end type gridcell_ch4_type
#endif

!----------------------------------------------------
! gridcell nitrogen state variables structure
!----------------------------------------------------
type, public :: gridcell_nstate_type
   type(column_nstate_type):: cns_a            !column-level nitrogen state variables averaged to gridcell
end type gridcell_nstate_type

!----------------------------------------------------
! gridcell VOC state variables structure
!----------------------------------------------------
type, public :: gridcell_vstate_type
   type(column_vstate_type):: cvs_a            !column-level VOC state variables averaged to gridcell
end type gridcell_vstate_type

!----------------------------------------------------
! gridcell VOC emission factor variables structure (heald)
!----------------------------------------------------
type, public :: gridcell_efstate_type
   real(r8), pointer      :: efisop(:,:)    ! isoprene emission factors
end type gridcell_efstate_type

!----------------------------------------------------
! gridcell dust state variables structure
!----------------------------------------------------
type, public :: gridcell_dstate_type
   type(column_dstate_type):: cds_a            !column-level dust state variables averaged to gridcell
end type gridcell_dstate_type

#if (defined CNDV)
!----------------------------------------------------
! gridcell DGVM state variables structure
!----------------------------------------------------
type, public :: gridcell_dgvstate_type
   real(r8), pointer :: agdd20(:)      !20-yr running mean of agdd
   real(r8), pointer :: tmomin20(:)    !20-yr running mean of tmomin
   real(r8), pointer :: t10min(:)      !ann minimum of 10-day running mean (K)
end type gridcell_dgvstate_type
#endif

!----------------------------------------------------
! gridcell energy flux variables structure
!----------------------------------------------------
type, public :: gridcell_eflux_type
   type(column_eflux_type):: cef_a             !column-level energy flux variables averaged to gridcell
   real(r8), pointer :: eflx_sh_totg(:)   ! total grid-level sensible heat flux
   real(r8), pointer :: eflx_dynbal(:)    ! dynamic land cover change conversion energy flux
end type gridcell_eflux_type

!----------------------------------------------------
! gridcell momentum flux variables structure
!-- -------------------------------------------------
type, public :: gridcell_mflux_type
   type(pft_mflux_type):: pmf_a                !pft-level momentum flux variables averaged to gridcell
end type gridcell_mflux_type

!----------------------------------------------------
! gridcell water flux variables structure
!----------------------------------------------------
type, public :: gridcell_wflux_type
   type(column_wflux_type):: cwf_a           ! column-level water flux variables averaged to gridcell
   real(r8), pointer :: qflx_runoffg(:)      ! total grid-level liq runoff
   real(r8), pointer :: qflx_snwcp_iceg(:)   ! total grid-level ice runoff
   real(r8), pointer :: qflx_liq_dynbal(:)   ! liq dynamic land cover change conversion runoff flux
   real(r8), pointer :: qflx_ice_dynbal(:)   ! ice dynamic land cover change conversion runoff flux
   real(r8), pointer :: qflx_floodg(:)       ! total grid-level flood water flux
end type gridcell_wflux_type

!----------------------------------------------------
! gridcell carbon flux variables structure
!----------------------------------------------------
type, public :: gridcell_cflux_type
   type(column_cflux_type):: ccf_a             !column-level carbon flux variables averaged to gridcell
end type gridcell_cflux_type

!----------------------------------------------------
! gridcell nitrogen flux variables structure
!----------------------------------------------------
type, public :: gridcell_nflux_type
   type(column_nflux_type):: cnf_a             !column-level nitrogen flux variables averaged to gridcell
end type gridcell_nflux_type

!----------------------------------------------------
! gridcell dust flux variables structure
!----------------------------------------------------
type, public :: gridcell_dflux_type
   type(pft_dflux_type):: pdf_a                !pft-level dust flux variables averaged to gridcell
end type gridcell_dflux_type

!----------------------------------------------------
! End definition of structures defined at the gridcell_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of spatial scaling hierarchy
!----------------------------------------------------

!----------------------------------------------------
! define the pft structure
!----------------------------------------------------

type, public :: pft_type

   ! g/l/c/p hierarchy, local g/l/c/p cells only
   integer, pointer :: column(:)        !index into column level quantities
   real(r8), pointer :: wtcol(:)        !weight (relative to column) 
   integer, pointer :: landunit(:)      !index into landunit level quantities
   real(r8), pointer :: wtlunit(:)      !weight (relative to landunit) 
   integer, pointer :: gridcell(:)      !index into gridcell level quantities
   real(r8), pointer :: wtgcell(:)      !weight (relative to gridcell) 

   ! topological mapping functionality
   integer , pointer :: itype(:)        !pft vegetation 
   integer , pointer :: mxy(:)          !m index for laixy(i,j,m),etc.
   logical , pointer :: active(:)       !true=>do computations on this pft (see reweightMod for details)

   ! conservation check structures for the pft level
   type(energy_balance_type)   :: pebal !energy balance structure
   type(water_balance_type)    :: pwbal !water balance structure
   type(carbon_balance_type)   :: pcbal !carbon balance structure
   type(nitrogen_balance_type) :: pnbal !nitrogen balance structure
   
#if (defined CNDV)
   ! DGVM state variables
   type(pft_dgvstate_type) :: pdgvs     !pft DGVM state variables
#endif
   
   ! CN ecophysiological variables
   type(pft_epv_type)    :: pepv        !pft ecophysiological variables
   
   ! state variables defined at the pft level
   type(pft_pstate_type) :: pps         !physical state variables
   type(pft_psynstate_type)::ppsyns     !photosynthesis relevant variables
   type(pft_estate_type) :: pes         !pft energy state
   type(pft_wstate_type) :: pws         !pft water state
   type(pft_cstate_type) :: pcs         !pft carbon state
   type(pft_nstate_type) :: pns         !pft nitrogen state
   type(pft_vstate_type) :: pvs         !pft VOC state

   ! flux variables defined at the pft level
   type(pft_eflux_type)  :: pef         !pft energy flux
   type(pft_mflux_type)  :: pmf         !pft momentum flux
   type(pft_wflux_type)  :: pwf         !pft water flux
   type(pft_cflux_type)  :: pcf         !pft carbon flux
   type(pft_nflux_type)  :: pnf         !pft nitrogen flux
   type(pft_vflux_type)  :: pvf         !pft VOC flux
   type(pft_dflux_type)  :: pdf         !pft dust flux
   type(pft_depvd_type)  :: pdd         !dry dep velocity
   
   type(pft_cstate_type) :: pc13s       !pft carbon-13 state
   type(pft_cflux_type)  :: pc13f       !pft carbon-13 flux

   type(pft_cstate_type) :: pc14s       !pft carbon-14 state
   type(pft_cflux_type)  :: pc14f       !pft carbon-14 flux
   
end type pft_type

!----------------------------------------------------
! define the column structure
!----------------------------------------------------

type, public :: column_type

   type(pft_type)   :: p       !plant functional type (pft) data structure 

   ! g/l/c/p hierarchy, local g/l/c/p cells only
   integer , pointer :: landunit(:)     !index into landunit level quantities
   real(r8), pointer :: wtlunit(:)      !weight (relative to landunit)
   integer , pointer :: gridcell(:)     !index into gridcell level quantities
   real(r8), pointer :: wtgcell(:)      !weight (relative to gridcell)
   integer , pointer :: pfti(:)         !beginning pft index for each column
   integer , pointer :: pftf(:)         !ending pft index for each column
   integer , pointer :: npfts(:)        !number of pfts for each column
   
   ! topological mapping functionality
   integer , pointer :: itype(:)        !column type
   logical , pointer :: active(:)       !true=>do computations on this column (see reweightMod for details)

   ! conservation check structures for the column level
   type(energy_balance_type)   :: cebal !energy balance structure
   type(water_balance_type)    :: cwbal !water balance structure
   type(carbon_balance_type)   :: ccbal !carbon balance structure
   type(nitrogen_balance_type) :: cnbal !nitrogen balance structure
   
   ! state variables defined at the column level
   type(column_pstate_type) :: cps      !column physical state variables
   type(column_estate_type) :: ces      !column energy state
   type(column_wstate_type) :: cws      !column water state
   type(column_cstate_type) :: ccs      !column carbon state
   type(column_nstate_type) :: cns      !column nitrogen state
   type(column_dstate_type) :: cds      !column dust state
   
   ! flux variables defined at the column level
   type(column_eflux_type) :: cef       !column energy flux
   type(column_mflux_type) :: cmf       !column momentum flux
   type(column_wflux_type) :: cwf       !column water flux
   type(column_cflux_type) :: ccf       !column carbon flux
#ifdef LCH4
   type(column_ch4_type)   :: cch4      !column CH4 variables
#endif
   type(column_nflux_type) :: cnf       !column nitrogen flux
   type(column_dflux_type) :: cdf       !column dust flux

#if (defined CNDV)
   ! dgvm variables defined at the column level
   type (column_dgvstate_type) :: cdgvs !column DGVM structure
#endif
   
   type(column_cstate_type) :: cc13s    !column carbon-13 state
   type(column_cflux_type)  :: cc13f    !column carbon-13 flux

   type(column_cstate_type) :: cc14s    !column carbon-14 state
   type(column_cflux_type)  :: cc14f    !column carbon-14 flux
   
end type column_type

!----------------------------------------------------
! define the geomorphological land unit structure
!----------------------------------------------------

type, public :: landunit_type
   type(column_type) :: c                 !column data structure (soil/snow/canopy columns)

   ! g/l/c/p hierarchy, local g/l/c/p cells only
   integer , pointer :: gridcell(:)       !index into gridcell level quantities
   real(r8), pointer :: wtgcell(:)        !weight (relative to gridcell)
   integer , pointer :: coli(:)           !beginning column index per landunit
   integer , pointer :: colf(:)           !ending column index for each landunit
   integer , pointer :: ncolumns(:)       !number of columns for each landunit
   integer , pointer :: pfti(:)           !beginning pft index for each landunit
   integer , pointer :: pftf(:)           !ending pft index for each landunit
   integer , pointer :: npfts(:)          !number of pfts for each landunit

   ! Urban canyon related properties
   real(r8), pointer :: canyon_hwr(:)     ! urban landunit canyon height to width ratio (-)   
   real(r8), pointer :: wtroad_perv(:)    ! urban landunit weight of pervious road column to total road (-)
   real(r8), pointer :: wtlunit_roof(:)   ! weight of roof with respect to urban landunit (-)

   ! Urban related info MV - this should be moved to land physical state - MV
   real(r8), pointer :: ht_roof(:)        ! height of urban roof (m)
   real(r8), pointer :: wind_hgt_canyon(:)! height above road at which wind in canyon is to be computed (m)
   real(r8), pointer :: z_0_town(:)       ! urban landunit momentum roughness length (m)
   real(r8), pointer :: z_d_town(:)       ! urban landunit displacement height (m)
   
   ! topological mapping functionality
   integer , pointer :: itype(:)        !landunit type
   logical , pointer :: ifspecial(:)    !true=>landunit is not vegetated
   logical , pointer :: lakpoi(:)       !true=>lake point
   logical , pointer :: urbpoi(:)       !true=>urban point
   logical , pointer :: glcmecpoi(:)    !true=>glacier_mec point
   integer , pointer :: udenstype(:)    !urban density type
   logical , pointer :: active(:)       !true=>do computations on this landunit (see reweightMod for details)

   ! state variables defined at the land unit level
   type(landunit_pstate_type) :: lps    !land unit physical state variables
   
   ! flux variables defined at the landunit level
   type(landunit_eflux_type) :: lef     !average of energy fluxes all columns
end type landunit_type

!----------------------------------------------------
! define the gridcell structure
!----------------------------------------------------

type, public :: gridcell_type

   type(landunit_type) :: l             !geomorphological landunits

   ! g/l/c/p hierarchy, local g/l/c/p cells only
   integer, pointer :: luni(:)          !beginning landunit index 
   integer, pointer :: lunf(:)          !ending landunit index 
   integer, pointer :: nlandunits(:)    !number of landunit for each gridcell
   integer, pointer :: coli(:)          !beginning column index
   integer, pointer :: colf(:)          !ending column index
   integer, pointer :: ncolumns(:)      !number of columns for each gridcell
   integer, pointer :: pfti(:)          !beginning pft index
   integer, pointer :: pftf(:)          !ending pft index
   integer, pointer :: npfts(:)         !number of pfts for each gridcell

   ! topological mapping functionality, local 1d gdc arrays
   integer , pointer :: gindex(:)       !global index
   real(r8), pointer :: area(:)         !total land area, gridcell (km^2)
   real(r8), pointer :: lat(:)          !latitude (radians)
   real(r8), pointer :: lon(:)          !longitude (radians)
   real(r8), pointer :: latdeg(:)       !latitude (degrees)
   real(r8), pointer :: londeg(:)       !longitude (degrees)
   integer , pointer :: gindex_a(:)     !"atm" global index
   real(r8), pointer :: lat_a(:) 	!"atm" latitude (radians) for albedo
   real(r8), pointer :: lon_a(:)        !"atm" longitude (radians) for albedo
   real(r8), pointer :: latdeg_a(:)     !"atm" latitude (degrees) for albedo
   real(r8), pointer :: londeg_a(:)     !"atm" longitude (degrees) for albedo

   real(r8), pointer :: gris_mask(:)    !Greenland ice sheet mask 
   real(r8), pointer :: gris_area(:)    !Greenland ice-covered area per gridcell (km^2)
   real(r8), pointer :: aais_mask(:)    !Antarctic ice sheet mask 
   real(r8), pointer :: aais_area(:)    !Antarctic ice-covered area per gridcell (km^2)
   real(r8), pointer :: tws(:)          !total water storage (mm H2O)

#if (defined CNDV)
   ! dgvm variables defined at the gridcell level
   type(gridcell_dgvstate_type):: gdgvs !gridcell DGVM structure
#endif

   ! state variables defined at the gridcell level
   type(gridcell_estate_type) :: ges    !average of energy states all landunits
   type(gridcell_wstate_type) :: gws    !average of water states all landunits
   type(gridcell_efstate_type):: gve	!gridcell VOC emission factors
   
   ! flux variables defined at the gridcell level
   type(gridcell_eflux_type) :: gef     !average of energy fluxes all landunits
   type(gridcell_wflux_type) :: gwf     !average of water fluxes all landunits
#ifdef LCH4
   type(gridcell_ch4_type)   :: gch4    !average of CH4 fluxes all landunits
#endif

end type gridcell_type

!----------------------------------------------------
! define the top-level (model) structure 
!----------------------------------------------------

type, public :: model_type
   ! lower level in hierarch
   type(gridcell_type) :: g    !gridicell data structure
   integer  :: ngridcells      !number of gridcells for this process
   real(r8) :: area            !total land area for all gridcells (km^2)

end type model_type

!----------------------------------------------------
! End definition of spatial scaling hierarchy
!----------------------------------------------------
!*******************************************************************************

!*******************************************************************************
!----------------------------------------------------
! Declare single instance of clmtype
!----------------------------------------------------
type(model_type)    , public, target     , save :: clm3

!----------------------------------------------------
! Declare single instance of array of ecophysiological constant types
!----------------------------------------------------
type(pft_epc_type), public, target, save :: pftcon


!----------------------------------------------------
! Declare single instance of array of decomposition cascade constant types
!----------------------------------------------------
type(decomp_cascade_type), public, target, save :: decomp_cascade_con

#if (defined CNDV)
!----------------------------------------------------
! Declare single instance of array of dgvm ecophysiological constant types
!----------------------------------------------------
type(pft_dgvepc_type), public, target, save :: dgv_pftcon
#endif

character(len=16), parameter, public :: grlnd  = 'lndgrid'      ! name of lndgrid
character(len=16), parameter, public :: namea  = 'gridcellatm'  ! name of atmgrid
character(len=16), parameter, public :: nameg  = 'gridcell'     ! name of gridcells
character(len=16), parameter, public :: namel  = 'landunit'     ! name of landunits
character(len=16), parameter, public :: namec  = 'column'       ! name of columns
character(len=16), parameter, public :: namep  = 'pft'          ! name of pfts

!
!EOP
!----------------------------------------------------------------------- 
end module clmtype  
