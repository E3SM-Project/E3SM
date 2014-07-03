module SLakeHydrologyMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SLakeHydrologyMod
!
! !DESCRIPTION:
! Calculation of Lake Hydrology. Full hydrology, aerosol deposition, etc. of snow layers is
! done. However, there is no infiltration, and the water budget is balanced with 
! qflx_qrgwl. Lake water mass is kept constant. The soil is simply maintained at
! volumetric saturation if ice melting frees up pore space. Likewise, if the water
! portion alone at some point exceeds pore capacity, it is reduced. This is consistent
! with the possibility of initializing the soil layer with excess ice.
! 
! If snow layers are present over an unfrozen lake, and the top layer of the lake
! is capable of absorbing the latent heat without going below freezing, 
! the snow-water is runoff and the latent heat is subtracted from the lake.
!
! Minimum snow layer thickness for lakes has been increased to avoid instabilities with 30 min timestep.
! Also frost / dew is prevented from being added to top snow layers that have already melted during the phase change step.
!
! !PUBLIC TYPES:
  implicit none
  save
  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SLakeHydrology        ! Calculates soil/snow hydrology
!
! !REVISION HISTORY:
! Created by Zack Subin, 2009
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SLakeHydrology
!
! !INTERFACE:
  subroutine SLakeHydrology(lbc, ubc, lbp, ubp, num_lakec, filter_lakec, &
                               num_lakep, filter_lakep &
  ! Snow filter for lakes is not returned to driver.  That's okay, because it looks like it is only
  ! needed for the call to SnowAge_grain, which will be done at the bottom of this module.
                        )
!
! !DESCRIPTION:
!
! WARNING: This subroutine assumes lake columns have one and only one pft.
!
! Sequence is:
!  SLakeHydrology:
!    Do needed tasks from Hydrology1, Biogeophysics2, & top of Hydrology2.
!    -> SnowWater:             change of snow mass and snow water onto soil
!    -> SnowCompaction:        compaction of snow layers
!    -> CombineSnowLayers:     combine snow layers that are thinner than minimum
!    -> DivideSnowLayers:      subdivide snow layers that are thicker than maximum
!    Add water to soil if melting has left it with open pore space.
!    If snow layers are found above a lake with unfrozen top layer, whose top
!    layer has enough heat to melt all the snow ice without freezing, do so
!    and eliminate the snow layers.
!    Cleanup and do water balance.
!    Do SNICAR stuff and diagnostics.
!    Call SnowAge_grain (it must be done here because the snow filters at the driver level are non-lakec only.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_atmlnd      , only : clm_a2l
    use clm_varcon      , only : denh2o, denice, spval, hfus, tfrz, cpliq, cpice
    use SLakeCon        , only : lsadz
    use clm_varpar      , only : nlevsno, nlevgrnd, nlevsoi
    use clm_varctl      , only : iulog
    use SnowHydrologyMod, only : SnowCompaction, CombineSnowLayers, &
                                 SnowWater, BuildSnowFilter
    use SnowHydrologyMod, only : DivideSnowLayers_Lake
    use clm_time_manager, only : get_step_size, is_perpetual
    use SNICARMod           , only : SnowAge_grain, snw_rds_min
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                  ! column bounds
    integer, intent(in) :: lbp, ubp                  ! pft bounds
    integer, intent(in) :: num_lakec                 ! number of column lake points in column filter
    integer, intent(in) :: filter_lakec(ubc-lbc+1)   ! column filter for lake points
    integer, intent(in) :: num_lakep                 ! number of pft lake points in column filter
    integer, intent(in) :: filter_lakep(ubp-lbp+1)   ! pft filter for lake points
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Created by Zack Subin
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    real(r8), pointer :: frac_sno_eff(:)  ! needed for snicar code
    real(r8), pointer :: qflx_floodg(:)   ! gridcell flux of flood water from RTM
    real(r8), pointer :: qflx_floodc(:)   ! column flux of flood water from RTM
    real(r8), pointer :: frost_table(:)   ! frost table depth (m)
    real(r8), pointer :: zwt_perched(:)   ! perched water table depth (m)
    real(r8), pointer :: qflx_drain_perched(:)! perched wt sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_h2osfc_surf(:)  ! surface water runoff (mm H2O /s)
    real(r8), pointer :: qflx_snow_melt(:)! net snow melt
    real(r8), pointer :: qflx_rsub_sat(:) !soil saturation excess [mm h2o/s]
    integer , pointer :: pcolumn(:)       ! pft's column index
    integer , pointer :: pgridcell(:)     ! pft's gridcell index
    integer , pointer :: cgridcell(:)     ! column's gridcell
    integer , pointer :: clandunit(:)     ! column's landunit
    real(r8), pointer :: watsat(:,:)      ! volumetric soil water at saturation (porosity)
    real(r8), pointer :: z(:,:)           ! layer depth  (m)
    real(r8), pointer :: dz_lake(:,:)     ! layer thickness for lake (m)
    real(r8), pointer :: forc_rain(:)     ! rain rate [mm/s]
    real(r8), pointer :: forc_snow(:)     ! snow rate [mm/s]
    real(r8), pointer :: begwb(:)         ! water mass begining of the time step
    real(r8), pointer :: qflx_evap_tot(:) ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(r8), pointer :: forc_t(:)        ! atmospheric temperature (Kelvin)
    logical , pointer :: do_capsnow(:)    ! true => do snow capping
    real(r8), pointer :: t_grnd(:)        ! ground temperature (Kelvin)
    real(r8), pointer :: qflx_evap_soi(:) ! soil evaporation (mm H2O/s) (+ = to atm)
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: dz(:,:)          ! layer thickness depth (m)
    real(r8), pointer :: zi(:,:)          ! interface depth (m)
    integer , pointer :: snl(:)           ! number of snow layers
    real(r8), pointer :: h2osno(:)        ! snow water (mm H2O)
    real(r8), pointer :: snow_depth(:)        ! snow height (m)
    real(r8), pointer :: lake_icefrac(:,:)! mass fraction of lake layer that is frozen
    real(r8), pointer :: t_lake(:,:)      ! lake temperature (Kelvin)
    real(r8), pointer :: qflx_snomelt(:)  ! snow melt (mm H2O /s)
    real(r8), pointer :: eflx_snomelt(:)  ! snow melt heat flux (W/m**2)
    real(r8), pointer :: eflx_sh_tot(:)   ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_grnd(:)  ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_soil_grnd(:)! heat flux into snow / lake (W/m**2) [+ = into soil]
                                          ! Here this includes the whole lake radiation absorbed.
    real(r8), pointer :: eflx_gnet(:)     ! net heat flux into ground (W/m**2)
    real(r8), pointer :: eflx_grnd_lake(:)! net heat flux into lake / snow surface, excluding light transmission (W/m**2)


!
! local pointers to implicit out arguments
!
    real(r8), pointer :: endwb(:)         ! water mass end of the time step
    real(r8), pointer :: snowice(:)       ! average snow ice lens
    real(r8), pointer :: snowliq(:)       ! average snow liquid water
    real(r8), pointer :: t_soisno(:,:)    ! snow temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)  ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_vol(:,:)  ! volumetric soil water [m3/m3]
    real(r8), pointer :: qflx_drain(:)    ! sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_surf(:)     ! surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_infl(:)     ! infiltration (mm H2O /s)
    real(r8), pointer :: qflx_qrgwl(:)    ! qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: qflx_runoff(:)   ! total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(r8), pointer :: qcharge(:)       ! aquifer recharge rate (mm/s)
    real(r8), pointer :: qflx_top_soil(:)      ! net water input into soil from top (mm/s)
    real(r8), pointer :: qflx_sl_top_soil(:)   ! liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
    real(r8), pointer :: qflx_prec_grnd(:)     ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_prec_grnd_col(:) ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_snow_grnd_pft(:) ! snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snow_grnd_col(:) ! snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_rain_grnd(:)     ! rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: frac_iceold(:,:)      ! fraction of ice relative to the tot water
    real(r8), pointer :: qflx_evap_tot_col(:)  ! pft quantity averaged to the column (assuming one pft)
    real(r8) ,pointer :: soilalpha(:)       ! factor that reduces ground saturated specific humidity (-)
    real(r8), pointer :: zwt(:)             ! water table depth
    real(r8), pointer :: fcov(:)            ! fractional area with water table at surface
    real(r8), pointer :: fsat(:)            ! fractional area with water table at surface
    real(r8), pointer :: rootr_column(:,:)  ! effective fraction of roots in each soil layer
    real(r8), pointer :: qflx_evap_grnd(:)  ! ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_evap_grnd_col(:) ! ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_sub_snow(:)   ! sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow_col(:) ! sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_snow(:)   ! surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_snow_col(:) ! surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_grnd(:)   ! ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_grnd_col(:)    ! ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: qflx_rain_grnd_col(:)   ! rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snwcp_ice_col(:)   ! excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snwcp_ice(:)       ! excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snwcp_liq_col(:)   ! excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snwcp_liq(:)       ! excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_irrig(:)           ! irrigation flux (mm H2O /s)
    !New SNICAR variables from Hydrology1
    real(r8), pointer :: snw_rds(:,:)          ! effective snow grain radius (col,lyr) [microns, m^-6]
    real(r8), pointer :: mss_bcpho(:,:)        ! mass of hydrophobic BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bcphi(:,:)        ! mass of hydrophilic BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bctot(:,:)        ! total mass of BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bc_col(:)         ! total column mass of BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bc_top(:)         ! total top-layer mass of BC (col,lyr) [kg]
    real(r8), pointer :: mss_ocpho(:,:)        ! mass of hydrophobic OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocphi(:,:)        ! mass of hydrophilic OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_octot(:,:)        ! total mass of OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_oc_col(:)         ! total column mass of OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_oc_top(:)         ! total top-layer mass of OC (col,lyr) [kg]
    real(r8), pointer :: mss_dst1(:,:)         ! mass of dust species 1 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst2(:,:)         ! mass of dust species 2 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst3(:,:)         ! mass of dust species 3 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst4(:,:)         ! mass of dust species 4 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dsttot(:,:)       ! total mass of dust in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst_col(:)        ! total column mass of dust in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst_top(:)        ! total top-layer mass of dust in snow (col,lyr) [kg]
    !Additional SNICAR variables from Hydrology2
    real(r8), pointer :: mss_cnc_bcphi(:,:) ! mass concentration of BC species 1 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_bcpho(:,:) ! mass concentration of BC species 2 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_ocphi(:,:) ! mass concentration of OC species 1 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_ocpho(:,:) ! mass concentration of OC species 2 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst1(:,:)  ! mass concentration of dust species 1 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst2(:,:)  ! mass concentration of dust species 2 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst3(:,:)  ! mass concentration of dust species 3 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst4(:,:)  ! mass concentration of dust species 4 (col,lyr) [kg/kg]

    
    ! New Diagnostics
    real(r8), pointer :: snot_top(:)        ! snow temperature in top layer (col) [K]
    real(r8), pointer :: dTdz_top(:)        ! temperature gradient in top layer (col) [K m-1]
    real(r8), pointer :: snw_rds_top(:)     ! effective snow grain size, top layer(col) [microns]
    real(r8), pointer :: sno_liq_top(:)     ! liquid water fraction in top snow layer (col) [frc]
    real(r8), pointer :: h2osno_top(:)      ! mass of snow in top layer (col) [kg]


!!!!!!
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: p,fp,g,l,c,j,fc,jtop             ! indices
    integer  :: num_shlakesnowc                  ! number of column snow points
    integer  :: filter_shlakesnowc(ubc-lbc+1)    ! column filter for snow points
    integer  :: num_shlakenosnowc                ! number of column non-snow points
    integer  :: filter_shlakenosnowc(ubc-lbc+1)  ! column filter for non-snow points
    real(r8) :: dtime                        ! land model time step (sec)
    integer  :: newnode                      ! flag when new snow node is set, (1=yes, 0=no)
    real(r8) :: dz_snowf                     ! layer thickness rate change due to precipitation [mm/s]
    real(r8) :: bifall                       ! bulk density of newly fallen dry snow [kg/m3]
    real(r8) :: fracsnow(lbp:ubp)            ! frac of precipitation that is snow
    real(r8) :: fracrain(lbp:ubp)            ! frac of precipitation that is rain
    real(r8) :: qflx_prec_grnd_snow(lbp:ubp) ! snow precipitation incident on ground [mm/s]
    real(r8) :: qflx_prec_grnd_rain(lbp:ubp) ! rain precipitation incident on ground [mm/s]
    real(r8) :: qflx_evap_soi_lim            ! temporary evap_soi limited by top snow layer content [mm/s]
    real(r8) :: h2osno_temp                  ! temporary h2osno [kg/m^2]
    real(r8) :: sumsnowice(lbc:ubc)          ! sum of snow ice if snow layers found above unfrozen lake [kg/m&2]
    logical  :: unfrozen(lbc:ubc)            ! true if top lake layer is unfrozen with snow layers above
    real(r8) :: heatrem                      ! used in case above [J/m^2]
    real(r8) :: heatsum(lbc:ubc)             ! used in case above [J/m^2]
    real(r8) :: snowmass                     ! liquid+ice snow mass in a layer [kg/m2]
    real(r8) :: snowcap_scl_fct              ! temporary factor used to correct for snow capping
    real(r8), parameter :: snow_bd = 250._r8 ! assumed snow bulk density (for lakes w/out resolved snow layers) [kg/m^3]
                                             ! Should only be used for frost below.
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (gridcell-level)

    forc_rain => clm_a2l%forc_rain
    forc_snow => clm_a2l%forc_snow
    forc_t    => clm_a2l%forc_t

    ! Assign local pointers to derived subtypes components (column-level)

    frac_sno_eff  => clm3%g%l%c%cps%frac_sno_eff 
    cgridcell     => clm3%g%l%c%gridcell
    clandunit     => clm3%g%l%c%landunit
    snl           => clm3%g%l%c%cps%snl
    t_grnd        => clm3%g%l%c%ces%t_grnd
    h2osno        => clm3%g%l%c%cws%h2osno
    snowice       => clm3%g%l%c%cws%snowice
    snowliq       => clm3%g%l%c%cws%snowliq
    zwt           => clm3%g%l%c%cws%zwt
    fcov          => clm3%g%l%c%cws%fcov
    fsat           => clm3%g%l%c%cws%fsat
    qcharge       => clm3%g%l%c%cws%qcharge
    qflx_top_soil => clm3%g%l%c%cwf%qflx_top_soil
    qflx_prec_grnd_col => clm3%g%l%c%cwf%pwf_a%qflx_prec_grnd
    qflx_evap_grnd_col => clm3%g%l%c%cwf%pwf_a%qflx_evap_grnd
    qflx_dew_grnd_col  => clm3%g%l%c%cwf%pwf_a%qflx_dew_grnd
    qflx_dew_snow_col  => clm3%g%l%c%cwf%pwf_a%qflx_dew_snow
    qflx_sub_snow_col  => clm3%g%l%c%cwf%pwf_a%qflx_sub_snow
    qflx_snwcp_ice_col => clm3%g%l%c%cwf%pwf_a%qflx_snwcp_ice
    watsat        => clm3%g%l%c%cps%watsat
    z             => clm3%g%l%c%cps%z
    dz            => clm3%g%l%c%cps%dz
    zi            => clm3%g%l%c%cps%zi
    t_soisno      => clm3%g%l%c%ces%t_soisno
    h2osoi_ice    => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq    => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_vol    => clm3%g%l%c%cws%h2osoi_vol
    qflx_drain    => clm3%g%l%c%cwf%qflx_drain
    qflx_surf     => clm3%g%l%c%cwf%qflx_surf
    qflx_infl     => clm3%g%l%c%cwf%qflx_infl
    qflx_qrgwl    => clm3%g%l%c%cwf%qflx_qrgwl
    qflx_runoff    => clm3%g%l%c%cwf%qflx_runoff
    qflx_irrig    =>  clm3%g%l%c%cwf%qflx_irrig
    endwb         => clm3%g%l%c%cwbal%endwb
    begwb         => clm3%g%l%c%cwbal%begwb
    dz_lake       => clm3%g%l%c%cps%dz_lake
    t_lake         => clm3%g%l%c%ces%t_lake
    lake_icefrac   => clm3%g%l%c%cws%lake_icefrac
    do_capsnow    => clm3%g%l%c%cps%do_capsnow
    snow_depth        => clm3%g%l%c%cps%snow_depth
    qflx_snow_grnd_col => clm3%g%l%c%cwf%pwf_a%qflx_snow_grnd
    frac_iceold   => clm3%g%l%c%cps%frac_iceold
    qflx_evap_tot_col => clm3%g%l%c%cwf%pwf_a%qflx_evap_tot
    soilalpha    => clm3%g%l%c%cws%soilalpha
    zwt          => clm3%g%l%c%cws%zwt
    rootr_column => clm3%g%l%c%cps%rootr_column
    qflx_rain_grnd_col => clm3%g%l%c%cwf%pwf_a%qflx_rain_grnd
    qflx_snomelt => clm3%g%l%c%cwf%qflx_snomelt
    eflx_snomelt => clm3%g%l%c%cef%eflx_snomelt
    ! Use column variables here.
    qflx_snwcp_ice_col => clm3%g%l%c%cwf%pwf_a%qflx_snwcp_ice
    qflx_snwcp_liq_col => clm3%g%l%c%cwf%pwf_a%qflx_snwcp_liq
    !SNICAR variables from Hydrology1
    snw_rds            => clm3%g%l%c%cps%snw_rds
    mss_bcpho          => clm3%g%l%c%cps%mss_bcpho
    mss_bcphi          => clm3%g%l%c%cps%mss_bcphi
    mss_bctot          => clm3%g%l%c%cps%mss_bctot
    mss_bc_col         => clm3%g%l%c%cps%mss_bc_col
    mss_bc_top         => clm3%g%l%c%cps%mss_bc_top
    mss_ocpho          => clm3%g%l%c%cps%mss_ocpho
    mss_ocphi          => clm3%g%l%c%cps%mss_ocphi
    mss_octot          => clm3%g%l%c%cps%mss_octot
    mss_oc_col         => clm3%g%l%c%cps%mss_oc_col
    mss_oc_top         => clm3%g%l%c%cps%mss_oc_top
    mss_dst1           => clm3%g%l%c%cps%mss_dst1
    mss_dst2           => clm3%g%l%c%cps%mss_dst2
    mss_dst3           => clm3%g%l%c%cps%mss_dst3
    mss_dst4           => clm3%g%l%c%cps%mss_dst4
    mss_dsttot         => clm3%g%l%c%cps%mss_dsttot
    mss_dst_col        => clm3%g%l%c%cps%mss_dst_col
    mss_dst_top        => clm3%g%l%c%cps%mss_dst_top
    ! Diagnostics
    snot_top          => clm3%g%l%c%cps%snot_top
    dTdz_top          => clm3%g%l%c%cps%dTdz_top
    snw_rds_top       => clm3%g%l%c%cps%snw_rds_top
    sno_liq_top       => clm3%g%l%c%cps%sno_liq_top
    h2osno_top        => clm3%g%l%c%cps%h2osno_top
    ! SNICAR variables from Hydrology2
    mss_cnc_bcphi     => clm3%g%l%c%cps%mss_cnc_bcphi
    mss_cnc_bcpho     => clm3%g%l%c%cps%mss_cnc_bcpho
    mss_cnc_ocphi     => clm3%g%l%c%cps%mss_cnc_ocphi
    mss_cnc_ocpho     => clm3%g%l%c%cps%mss_cnc_ocpho
    mss_cnc_dst1      => clm3%g%l%c%cps%mss_cnc_dst1
    mss_cnc_dst2      => clm3%g%l%c%cps%mss_cnc_dst2
    mss_cnc_dst3      => clm3%g%l%c%cps%mss_cnc_dst3
    mss_cnc_dst4      => clm3%g%l%c%cps%mss_cnc_dst4
    ! Flooding terms
    qflx_floodg       => clm_a2l%forc_flood
    qflx_floodc       => clm3%g%l%c%cwf%qflx_floodc
    frost_table       => clm3%g%l%c%cws%frost_table
    zwt_perched       => clm3%g%l%c%cws%zwt_perched
    qflx_drain_perched=> clm3%g%l%c%cwf%qflx_drain_perched
    qflx_h2osfc_surf  => clm3%g%l%c%cwf%qflx_h2osfc_surf
    qflx_snow_melt    => clm3%g%l%c%cwf%qflx_snow_melt
    qflx_rsub_sat     => clm3%g%l%c%cwf%qflx_rsub_sat
    qflx_top_soil     => clm3%g%l%c%cwf%qflx_top_soil
    qflx_sl_top_soil  => clm3%g%l%c%cwf%qflx_sl_top_soil



    ! Assign local pointers to derived type members (pft-level)

    pcolumn       => clm3%g%l%c%p%column
    pgridcell      => clm3%g%l%c%p%gridcell
    qflx_sub_snow  => clm3%g%l%c%p%pwf%qflx_sub_snow
    qflx_evap_grnd => clm3%g%l%c%p%pwf%qflx_evap_grnd
    qflx_dew_snow  => clm3%g%l%c%p%pwf%qflx_dew_snow
    qflx_dew_grnd  => clm3%g%l%c%p%pwf%qflx_dew_grnd
    qflx_prec_grnd     => clm3%g%l%c%p%pwf%qflx_prec_grnd
    qflx_snow_grnd_pft => clm3%g%l%c%p%pwf%qflx_snow_grnd
    qflx_rain_grnd     => clm3%g%l%c%p%pwf%qflx_rain_grnd
    qflx_evap_tot => clm3%g%l%c%p%pwf%qflx_evap_tot
    qflx_evap_soi  => clm3%g%l%c%p%pwf%qflx_evap_soi
    qflx_snwcp_ice => clm3%g%l%c%p%pwf%qflx_snwcp_ice
    qflx_snwcp_liq => clm3%g%l%c%p%pwf%qflx_snwcp_liq
    eflx_sh_tot    => clm3%g%l%c%p%pef%eflx_sh_tot
    eflx_sh_grnd   => clm3%g%l%c%p%pef%eflx_sh_grnd
    eflx_soil_grnd => clm3%g%l%c%p%pef%eflx_soil_grnd
    eflx_gnet      => clm3%g%l%c%p%pef%eflx_gnet
    eflx_grnd_lake => clm3%g%l%c%p%pef%eflx_grnd_lake


    ! Determine step size

    dtime = get_step_size()

    ! Add soil water to water balance.
    do j = 1, nlevgrnd
      do fc = 1, num_lakec
         c = filter_lakec(fc)
         begwb(c) = begwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
      end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Do precipitation onto ground, etc., from Hydrology1.

    do fp = 1, num_lakep
       p = filter_lakep(fp)
       g = pgridcell(p)
       c = pcolumn(p)

       qflx_prec_grnd_snow(p) = forc_snow(g)
       qflx_prec_grnd_rain(p) = forc_rain(g)
       qflx_prec_grnd(p) = qflx_prec_grnd_snow(p) + qflx_prec_grnd_rain(p)

       if (do_capsnow(c)) then
          qflx_snwcp_ice(p) = qflx_prec_grnd_snow(p)
          qflx_snwcp_liq(p) = qflx_prec_grnd_rain(p)
          qflx_snow_grnd_pft(p) = 0._r8
          qflx_rain_grnd(p) = 0._r8
       else
          qflx_snwcp_ice(p) = 0._r8
          qflx_snwcp_liq(p) = 0._r8
          qflx_snow_grnd_pft(p) = qflx_prec_grnd_snow(p)           ! ice onto ground (mm/s)
          qflx_rain_grnd(p)     = qflx_prec_grnd_rain(p)           ! liquid water onto ground (mm/s)
       end if
       ! Assuming one PFT; needed for below
       qflx_snow_grnd_col(c) = qflx_snow_grnd_pft(p)
       qflx_rain_grnd_col(c) = qflx_rain_grnd(p)

    end do ! (end pft loop)

    ! Determine snow height and snow water

    do fc = 1, num_lakec
       c = filter_lakec(fc)
       g = cgridcell(c)

       ! Use Alta relationship, Anderson(1976); LaChapelle(1961),
       ! U.S.Department of Agriculture Forest Service, Project F,
       ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

       if (do_capsnow(c)) then
          dz_snowf = 0._r8
       else
          if (forc_t(g) > tfrz + 2._r8) then
             bifall=50._r8 + 1.7_r8*(17.0_r8)**1.5_r8
          else if (forc_t(g) > tfrz - 15._r8) then
             bifall=50._r8 + 1.7_r8*(forc_t(g) - tfrz + 15._r8)**1.5_r8
          else
             bifall=50._r8
          end if
          dz_snowf = qflx_snow_grnd_col(c)/bifall
          snow_depth(c) = snow_depth(c) + dz_snowf*dtime
          h2osno(c) = h2osno(c) + qflx_snow_grnd_col(c)*dtime  ! snow water equivalent (mm)
       end if

       ! When the snow accumulation exceeds 40 mm, initialize snow layer
       ! Currently, the water temperature for the precipitation is simply set
       ! as the surface air temperature

       newnode = 0    ! flag for when snow node will be initialized
       if (snl(c) == 0 .and. qflx_snow_grnd_col(c) > 0.0_r8 .and. snow_depth(c) >= 0.01_r8 + lsadz) then
          newnode = 1
          snl(c) = -1
          dz(c,0) = snow_depth(c)                       ! meter
          z(c,0) = -0.5_r8*dz(c,0)
          zi(c,-1) = -dz(c,0)
          t_soisno(c,0) = min(tfrz, forc_t(g))      ! K
          h2osoi_ice(c,0) = h2osno(c)               ! kg/m2
          h2osoi_liq(c,0) = 0._r8                   ! kg/m2
          frac_iceold(c,0) = 1._r8

          ! intitialize SNICAR variables for fresh snow:
          snw_rds(c,0)    = snw_rds_min

          mss_bcpho(c,:)  = 0._r8
          mss_bcphi(c,:)  = 0._r8
          mss_bctot(c,:)  = 0._r8
          mss_bc_col(c)   = 0._r8
          mss_bc_top(c)   = 0._r8

          mss_ocpho(c,:)  = 0._r8
          mss_ocphi(c,:)  = 0._r8
          mss_octot(c,:)  = 0._r8
          mss_oc_col(c)   = 0._r8
          mss_oc_top(c)   = 0._r8

          mss_dst1(c,:)   = 0._r8
          mss_dst2(c,:)   = 0._r8
          mss_dst3(c,:)   = 0._r8
          mss_dst4(c,:)   = 0._r8
          mss_dsttot(c,:) = 0._r8
          mss_dst_col(c)  = 0._r8
          mss_dst_top(c)  = 0._r8
       end if

       ! The change of ice partial density of surface node due to precipitation.
       ! Only ice part of snowfall is added here, the liquid part will be added
       ! later.

       if (snl(c) < 0 .and. newnode == 0) then
          h2osoi_ice(c,snl(c)+1) = h2osoi_ice(c,snl(c)+1)+dtime*qflx_snow_grnd_col(c)
          dz(c,snl(c)+1) = dz(c,snl(c)+1)+dz_snowf*dtime
       end if

    end do

    ! Calculate sublimation and dew, adapted from HydrologyLake and Biogeophysics2.

    do fp = 1,num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)
       jtop = snl(c)+1

       qflx_evap_grnd(p) = 0._r8
       qflx_sub_snow(p) = 0._r8
       qflx_dew_snow(p) = 0._r8
       qflx_dew_grnd(p) = 0._r8

       if (jtop <= 0) then ! snow layers
          j = jtop
          ! Assign ground evaporation to sublimation from soil ice or to dew
          ! on snow or ground

          if (qflx_evap_soi(p) >= 0._r8) then
          ! for evaporation partitioning between liquid evap and ice sublimation, 
          ! use the ratio of liquid to (liquid+ice) in the top layer to determine split
          ! Since we're not limiting evap over lakes, but still can't remove more from top
          ! snow layer than there is there, create temp. limited evap_soi.
             qflx_evap_soi_lim = min(qflx_evap_soi(p), (h2osoi_liq(c,j)+h2osoi_ice(c,j))/dtime)
             if ((h2osoi_liq(c,j)+h2osoi_ice(c,j)) > 0._r8) then
                qflx_evap_grnd(p) = max(qflx_evap_soi_lim*(h2osoi_liq(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))), 0._r8)
             else
                qflx_evap_grnd(p) = 0._r8
             end if
             qflx_sub_snow(p) = qflx_evap_soi_lim - qflx_evap_grnd(p)     
          else
             ! if (t_grnd(c) < tfrz) then
             ! Causes rare blowup when thin snow layer should completely melt and has a high temp after thermal physics,
             ! but then is not eliminated in SnowHydrology because of this added frost. Also see below removal of
             ! completely melted single snow layer.
             if (t_grnd(c) < tfrz .and. t_soisno(c,j) < tfrz) then
                qflx_dew_snow(p) = abs(qflx_evap_soi(p))
             ! If top layer is only snow layer, SnowHydrology won't eliminate it if dew is added.
             else if (j < 0 .or. (t_grnd(c) == tfrz .and. t_soisno(c,j) == tfrz)) then
                qflx_dew_grnd(p) = abs(qflx_evap_soi(p))
             end if
          end if
          ! Update the pft-level qflx_snowcap
          ! This was moved in from Hydrology2 to keep all pft-level
          ! calculations out of Hydrology2
          if (do_capsnow(c)) then
              qflx_snwcp_ice(p) = qflx_snwcp_ice(p) + qflx_dew_snow(p) 
              qflx_snwcp_liq(p) = qflx_snwcp_liq(p) + qflx_dew_grnd(p)
          end if

       else ! No snow layers: do as in HydrologyLake but with actual clmtype variables
          if (qflx_evap_soi(p) >= 0._r8) then
             ! Sublimation: do not allow for more sublimation than there is snow
             ! after melt.  Remaining surface evaporation used for infiltration.
             qflx_sub_snow(p) = min(qflx_evap_soi(p), h2osno(c)/dtime)
             qflx_evap_grnd(p) = qflx_evap_soi(p) - qflx_sub_snow(p)
          else
             if (t_grnd(c) < tfrz-0.1_r8) then
                qflx_dew_snow(p) = abs(qflx_evap_soi(p))
             else
                qflx_dew_grnd(p) = abs(qflx_evap_soi(p))
             end if
          end if

          ! Update snow pack for dew & sub.

          h2osno_temp = h2osno(c)
          if (do_capsnow(c)) then
             h2osno(c) = h2osno(c) - qflx_sub_snow(p)*dtime
             qflx_snwcp_ice(p) = qflx_snwcp_ice(p) + qflx_dew_snow(p) 
             qflx_snwcp_liq(p) = qflx_snwcp_liq(p) + qflx_dew_grnd(p)
          else
             h2osno(c) = h2osno(c) + (-qflx_sub_snow(p)+qflx_dew_snow(p))*dtime
          end if
          if (h2osno_temp > 0._r8) then
             snow_depth(c) = snow_depth(c) * h2osno(c) / h2osno_temp
          else
             snow_depth(c) = h2osno(c)/snow_bd !Assume a constant snow bulk density = 250.
          end if

#if (defined PERGRO)
          if (abs(h2osno(c)) < 1.e-10_r8) h2osno(c) = 0._r8
          if (h2osno(c) == 0._r8) snow_depth(c) = 0._r8
#else
          h2osno(c) = max(h2osno(c), 0._r8)
#endif

       end if

       qflx_snwcp_ice_col(c) = qflx_snwcp_ice(p)
       qflx_snwcp_liq_col(c) = qflx_snwcp_liq(p)


    end do

    ! pft averages must be done here -- BEFORE SNOW CALCULATIONS AS THEY USE IT.
    ! for output to history tape and other uses
    ! (note that pft2col is called before SLakeHydrology, so we can't use that routine
    ! to do these column -> pft averages)
    do fp = 1,num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)
       qflx_evap_tot_col(c) = qflx_evap_tot(p)
       qflx_prec_grnd_col(c) = qflx_prec_grnd(p)
       qflx_evap_grnd_col(c) = qflx_evap_grnd(p)
       qflx_dew_grnd_col(c) = qflx_dew_grnd(p)
       qflx_dew_snow_col(c) = qflx_dew_snow(p)
       qflx_sub_snow_col(c) = qflx_sub_snow(p)

    enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Determine initial snow/no-snow filters (will be modified possibly by
    ! routines CombineSnowLayers and DivideSnowLayers below)

    call BuildSnowFilter(lbc, ubc, num_lakec, filter_lakec, &
         num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc)

    ! specify snow fraction
    do fc = 1, num_lakec
       c = filter_lakec(fc)
       if (h2osno(c) > 0.0_r8) then
          frac_sno_eff(c)     = 1._r8 
       else
          frac_sno_eff(c)     = 0._r8 
       endif
    enddo

    ! Determine the change of snow mass and the snow water onto soil

    call SnowWater(lbc, ubc, num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc)


    ! Determine soil hydrology
    ! Here this consists only of making sure that soil is saturated even as it melts and
    ! pore space opens up. Conversely, if excess ice is melting and the liquid water exceeds the
    ! saturation value, then remove water.

    do j = 1,nlevsoi  !nlevgrnd
                      ! changed to nlevsoi on 8/11/10 to make consistent with non-lake bedrock
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
          ! Could have changed during phase change! (Added 8/11/10)

          if (h2osoi_vol(c,j) < watsat(c,j)) then
             h2osoi_liq(c,j) = (watsat(c,j)*dz(c,j) - h2osoi_ice(c,j)/denice)*denh2o
          ! h2osoi_vol will be updated below, and this water addition will come from qflx_qrgwl
          else if (h2osoi_liq(c,j) > watsat(c,j)*denh2o*dz(c,j)) then
             h2osoi_liq(c,j) = watsat(c,j)*denh2o*dz(c,j)
          ! Another way to do this would be: if h2osoi_vol > watsat then remove min(h2osoi_liq,
            !(h2osoi_vol-watsat)*dz*denh2o) from h2osoi_liq.  The question is whether the excess ice
            ! melts first or last (or simultaneously) to the pore ice.  Because excess ice is often in chunks,
            ! requiring greater convergence of heat to melt, assume it melts last.
            ! This will also improve the initialization behavior or in an occasionally warm year, the excess ice
            ! won't start going away if a layer is briefly at freezing.

          ! Allow up to 10% excess ice over watsat in refreezing soil,
          ! e.g. heaving soil.  (As with > 10% excess ice modeling, and for the lake water,
          ! the thermal conductivity will be adjusted down to compensate for the fact that the nominal dz is smaller
          ! than the real soil volume.)  The current solution is consistent but perhaps unrealistic in real soils,
          ! where slow drainage may occur during freezing; drainage is only assumed to occur here when >10% excess
          ! ice melts. The latter is more likely to be permanent rather than seasonal anyway. Attempting to remove the
          ! ice volume after some has already frozen during the timestep would not conserve energy unless this were
          ! incorporated into the ice stream.

          end if

       end do
    end do
!!!!!!!!!!

    if (.not. is_perpetual()) then

       ! Natural compaction and metamorphosis.

       call SnowCompaction(lbc, ubc, num_shlakesnowc, filter_shlakesnowc)

       ! Combine thin snow elements

       call CombineSnowLayers(lbc, ubc, num_shlakesnowc, filter_shlakesnowc)

       ! Divide thick snow elements

       call DivideSnowLayers_Lake(lbc, ubc, num_shlakesnowc, filter_shlakesnowc)

    else

       do fc = 1, num_shlakesnowc
          c = filter_shlakesnowc(fc)
          h2osno(c) = 0._r8
       end do
       do j = -nlevsno+1,0
          do fc = 1, num_shlakesnowc
             c = filter_shlakesnowc(fc)
             if (j >= snl(c)+1) then
                h2osno(c) = h2osno(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
             end if
          end do
       end do

    end if

    ! Check for single completely unfrozen snow layer over lake.  Modeling this ponding is unnecessary and
    ! can cause instability after the timestep when melt is completed, as the temperature after melt can be
    ! excessive because the fluxes were calculated with a fixed ground temperature of freezing, but the 
    ! phase change was unable to restore the temperature to freezing.
    do fp = 1, num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)

       j = 0

       if (snl(c) == -1 .and. h2osoi_ice(c,j) == 0._r8) then
          ! Remove layer
          ! Take extra heat of layer and release to sensible heat in order to maintain energy conservation.
          heatrem = cpliq*h2osoi_liq(c,j)*(t_soisno(c,j) - tfrz)
          eflx_sh_tot(p) = eflx_sh_tot(p) + heatrem/dtime
          eflx_sh_grnd(p) = eflx_sh_grnd(p) + heatrem/dtime  ! Added this line 7/22/11 for consistency.
          eflx_soil_grnd(p) = eflx_soil_grnd(p) - heatrem/dtime
          eflx_gnet(p) = eflx_gnet(p) - heatrem/dtime
          eflx_grnd_lake(p) = eflx_grnd_lake(p) - heatrem/dtime
          qflx_sl_top_soil(c) = qflx_sl_top_soil(c) + h2osno(c)
          snl(c) = 0
          h2osno(c) = 0._r8
          snow_depth(c) = 0._r8
          ! Rest of snow layer book-keeping will be done below.
       end if
    end do


    ! Check for snow layers above lake with unfrozen top layer.  Mechanically,
    ! the snow will fall into the lake and melt or turn to ice.  If the top layer has
    ! sufficient heat to melt the snow without freezing, then that will be done.
    ! Otherwise, the top layer will undergo freezing, but only if the top layer will
    ! not freeze completely.  Otherwise, let the snow layers persist and melt by diffusion.
    do fc = 1, num_lakec
       c = filter_lakec(fc)

       if (t_lake(c,1) > tfrz .and. lake_icefrac(c,1) == 0._r8 .and. snl(c) < 0) then
          unfrozen(c) = .true.
       else
          unfrozen(c) = .false.
       end if
    end do

    do j = -nlevsno+1,0
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          if (unfrozen(c)) then
             if (j == -nlevsno+1) then
                sumsnowice(c) = 0._r8
                heatsum(c) = 0._r8
             end if
             if (j >= snl(c)+1) then
                sumsnowice(c) = sumsnowice(c) + h2osoi_ice(c,j)
                heatsum(c) = heatsum(c) + h2osoi_ice(c,j)*cpice*(tfrz - t_soisno(c,j)) &
                           + h2osoi_liq(c,j)*cpliq*(tfrz - t_soisno(c,j))
             end if
          end if
       end do
    end do

    do fc = 1, num_lakec
       c = filter_lakec(fc)

       if (unfrozen(c)) then
          heatsum(c) = heatsum(c) + sumsnowice(c)*hfus
          heatrem = (t_lake(c,1) - tfrz)*cpliq*denh2o*dz_lake(c,1) - heatsum(c)

          if (heatrem + denh2o*dz_lake(c,1)*hfus > 0._r8) then            
             ! Remove snow and subtract the latent heat from the top layer.
             qflx_snomelt(c) = qflx_snomelt(c) + h2osno(c)/dtime
             eflx_snomelt(c) = eflx_snomelt(c) + h2osno(c)*hfus/dtime
             ! update snow melt for this case
             qflx_snow_melt(c)     = qflx_snow_melt(c)  + qflx_snomelt(c)

             qflx_sl_top_soil(c) = qflx_sl_top_soil(c) + h2osno(c)

             h2osno(c) = 0._r8
             snow_depth(c) = 0._r8
             snl(c) = 0
             ! The rest of the bookkeeping for the removed snow will be done below.
             if (heatrem > 0._r8) then ! simply subtract the heat from the layer
                t_lake(c,1) = t_lake(c,1) - heatrem/(cpliq*denh2o*dz_lake(c,1))
             else !freeze part of the layer
                t_lake(c,1) = tfrz
                lake_icefrac(c,1) = -heatrem/(denh2o*dz_lake(c,1)*hfus)
             end if
          end if
       end if
    end do
!!!!!!!!!!!!

    ! Set empty snow layers to zero

    do j = -nlevsno+1,0
       do fc = 1, num_shlakesnowc
          c = filter_shlakesnowc(fc)
          if (j <= snl(c) .and. snl(c) > -nlevsno) then
             h2osoi_ice(c,j) = 0._r8
             h2osoi_liq(c,j) = 0._r8
             t_soisno(c,j) = 0._r8
             dz(c,j) = 0._r8
             z(c,j) = 0._r8
             zi(c,j-1) = 0._r8
          end if
       end do
    end do

    ! Build new snow filter

    call BuildSnowFilter(lbc, ubc, num_lakec, filter_lakec, &
         num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc)

    ! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice
    ! over all snow layers for history output

    do fc = 1, num_lakec
       c = filter_lakec(fc)
       snowice(c) = 0._r8
       snowliq(c) = 0._r8
    end do

    do j = -nlevsno+1, 0
       do fc = 1, num_shlakesnowc
          c = filter_shlakesnowc(fc)
          if (j >= snl(c)+1) then
             snowice(c) = snowice(c) + h2osoi_ice(c,j)
             snowliq(c) = snowliq(c) + h2osoi_liq(c,j)
          end if
       end do
    end do

    ! Determine ending water balance and volumetric soil water

    do fc = 1, num_lakec
       
       c = filter_lakec(fc)
       endwb(c) = h2osno(c)
    end do

    do j = 1, nlevgrnd
       do fc = 1, num_lakec
          c = filter_lakec(fc)
          endwb(c) = endwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
          h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
       end do
    end do

!!!!!!!!!!!!!
    ! Do history variables and set special landunit runoff (adapted from end of HydrologyLake)
    do fp = 1,num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)
       g = pgridcell(p)

       zwt_perched(c)    = spval
       frost_table(c)    = spval
       qflx_drain_perched(c)= 0._r8
       qflx_h2osfc_surf(c)  = 0._r8
       qflx_rsub_sat(c)     = 0._r8
       qflx_infl(c)      = 0._r8
       qflx_surf(c)      = 0._r8
       qflx_drain(c)     = 0._r8
       qflx_irrig(c)     = 0._r8
       rootr_column(c,:) = spval
       soilalpha(c)      = spval
       zwt(c)            = spval
       fcov(c)           = spval
       fsat(c)           = spval
       qcharge(c)        = spval

       ! Insure water balance using qflx_qrgwl
       qflx_qrgwl(c)     = forc_rain(g) + forc_snow(g) - qflx_evap_tot(p) - qflx_snwcp_ice(p) &
                         - (endwb(c)-begwb(c))/dtime + qflx_floodg(g)
       qflx_floodc(c)    = qflx_floodg(g)
       qflx_runoff(c)    = qflx_drain(c) + qflx_surf(c) + qflx_qrgwl(c)
       qflx_top_soil(c)  = qflx_prec_grnd_rain(p) + qflx_snomelt(c)

    enddo

    !  SNICAR Code and diagnostics

    !  Calculate column-integrated aerosol masses, and
    !  mass concentrations for radiative calculations and output
    !  (based on new snow level state, after SnowFilter is rebuilt.
    !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there 
    !  can be zero snow layers but an active column in filter)

    do fc = 1, num_shlakesnowc
       c = filter_shlakesnowc(fc)

       ! Zero column-integrated aerosol mass before summation
       mss_bc_col(c)  = 0._r8
       mss_oc_col(c)  = 0._r8
       mss_dst_col(c) = 0._r8

       do j = -nlevsno+1, 0

          ! layer mass of snow:
          snowmass = h2osoi_ice(c,j)+h2osoi_liq(c,j)

          ! Correct the top layer aerosol mass to account for snow capping. 
          ! This approach conserves the aerosol mass concentration
          ! (but not the aerosol amss) when snow-capping is invoked

          if (j == snl(c)+1) then
             if (do_capsnow(c)) then
                snowcap_scl_fct = snowmass / (snowmass+(qflx_snwcp_ice_col(c)*dtime))
                                                        ! Make sure column variable here
                mss_bcpho(c,j) = mss_bcpho(c,j)*snowcap_scl_fct
                mss_bcphi(c,j) = mss_bcphi(c,j)*snowcap_scl_fct
                mss_ocpho(c,j) = mss_ocpho(c,j)*snowcap_scl_fct
                mss_ocphi(c,j) = mss_ocphi(c,j)*snowcap_scl_fct
                
                mss_dst1(c,j)  = mss_dst1(c,j)*snowcap_scl_fct
                mss_dst2(c,j)  = mss_dst2(c,j)*snowcap_scl_fct
                mss_dst3(c,j)  = mss_dst3(c,j)*snowcap_scl_fct
                mss_dst4(c,j)  = mss_dst4(c,j)*snowcap_scl_fct 
             endif
          endif

          if (j >= snl(c)+1) then
             mss_bctot(c,j)     = mss_bcpho(c,j) + mss_bcphi(c,j)
             mss_bc_col(c)      = mss_bc_col(c)  + mss_bctot(c,j)
             mss_cnc_bcphi(c,j) = mss_bcphi(c,j) / snowmass
             mss_cnc_bcpho(c,j) = mss_bcpho(c,j) / snowmass

             mss_octot(c,j)     = mss_ocpho(c,j) + mss_ocphi(c,j)
             mss_oc_col(c)      = mss_oc_col(c)  + mss_octot(c,j)
             mss_cnc_ocphi(c,j) = mss_ocphi(c,j) / snowmass
             mss_cnc_ocpho(c,j) = mss_ocpho(c,j) / snowmass
             
             mss_dsttot(c,j)    = mss_dst1(c,j)  + mss_dst2(c,j) + mss_dst3(c,j) + mss_dst4(c,j)
             mss_dst_col(c)     = mss_dst_col(c) + mss_dsttot(c,j)
             mss_cnc_dst1(c,j)  = mss_dst1(c,j)  / snowmass
             mss_cnc_dst2(c,j)  = mss_dst2(c,j)  / snowmass
             mss_cnc_dst3(c,j)  = mss_dst3(c,j)  / snowmass
             mss_cnc_dst4(c,j)  = mss_dst4(c,j)  / snowmass
         
          else
             !set variables of empty snow layers to zero
             snw_rds(c,j)       = 0._r8

             mss_bcpho(c,j)     = 0._r8
             mss_bcphi(c,j)     = 0._r8
             mss_bctot(c,j)     = 0._r8
             mss_cnc_bcphi(c,j) = 0._r8
             mss_cnc_bcpho(c,j) = 0._r8

             mss_ocpho(c,j)     = 0._r8
             mss_ocphi(c,j)     = 0._r8
             mss_octot(c,j)     = 0._r8
             mss_cnc_ocphi(c,j) = 0._r8
             mss_cnc_ocpho(c,j) = 0._r8

             mss_dst1(c,j)      = 0._r8
             mss_dst2(c,j)      = 0._r8
             mss_dst3(c,j)      = 0._r8
             mss_dst4(c,j)      = 0._r8
             mss_dsttot(c,j)    = 0._r8
             mss_cnc_dst1(c,j)  = 0._r8
             mss_cnc_dst2(c,j)  = 0._r8
             mss_cnc_dst3(c,j)  = 0._r8
             mss_cnc_dst4(c,j)  = 0._r8
          endif
       enddo
       
       ! top-layer diagnostics
       h2osno_top(c)  = h2osoi_ice(c,snl(c)+1) + h2osoi_liq(c,snl(c)+1)
       mss_bc_top(c)  = mss_bctot(c,snl(c)+1)
       mss_oc_top(c)  = mss_octot(c,snl(c)+1)
       mss_dst_top(c) = mss_dsttot(c,snl(c)+1)
    enddo
    
    ! Zero mass variables in columns without snow
    do fc = 1, num_shlakenosnowc
       c = filter_shlakenosnowc(fc)
            
       h2osno_top(c)      = 0._r8
       snw_rds(c,:)       = 0._r8

       mss_bc_top(c)      = 0._r8
       mss_bc_col(c)      = 0._r8    
       mss_bcpho(c,:)     = 0._r8
       mss_bcphi(c,:)     = 0._r8
       mss_bctot(c,:)     = 0._r8
       mss_cnc_bcphi(c,:) = 0._r8
       mss_cnc_bcpho(c,:) = 0._r8

       mss_oc_top(c)      = 0._r8
       mss_oc_col(c)      = 0._r8    
       mss_ocpho(c,:)     = 0._r8
       mss_ocphi(c,:)     = 0._r8
       mss_octot(c,:)     = 0._r8
       mss_cnc_ocphi(c,:) = 0._r8
       mss_cnc_ocpho(c,:) = 0._r8

       mss_dst_top(c)     = 0._r8
       mss_dst_col(c)     = 0._r8
       mss_dst1(c,:)      = 0._r8
       mss_dst2(c,:)      = 0._r8
       mss_dst3(c,:)      = 0._r8
       mss_dst4(c,:)      = 0._r8
       mss_dsttot(c,:)    = 0._r8
       mss_cnc_dst1(c,:)  = 0._r8
       mss_cnc_dst2(c,:)  = 0._r8
       mss_cnc_dst3(c,:)  = 0._r8
       mss_cnc_dst4(c,:)  = 0._r8

       ! top-layer diagnostics (spval is not averaged when computing history fields)
       snot_top(c)        = spval
       dTdz_top(c)        = spval
       snw_rds_top(c)     = spval
       sno_liq_top(c)     = spval
    enddo

    !Must be done here because the snow filter used in Hydrology2 & the Driver are for non-lake columns.
    call SnowAge_grain(lbc, ubc, num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc)

  end subroutine SLakeHydrology

end module SLakeHydrologyMod
