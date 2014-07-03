module Hydrology2Mod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Hydrology2Mod
!
! !DESCRIPTION:
! Calculation of soil/snow hydrology.
!
! !USES:
   use shr_kind_mod, only : r8 => shr_kind_r8
   use clm_varctl,   only : iulog
   use abortutils,   only : endrun

! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Hydrology2        ! Calculates soil/snow hydrology
!
! !REVISION HISTORY:
! 2/28/02 Peter Thornton: Migrated to new data structures.
! 7/12/03 Forrest Hoffman ,Mariana Vertenstein : Migrated to vector code
! 11/05/03 Peter Thornton: Added calculation of soil water potential
!   for use in CN phenology code.
! 04/25/07 Keith Oleson: CLM3.5 Hydrology
!F. Li and S. Levis (11/06/12) for wf2 and tsoi17
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Hydrology2
!
! !INTERFACE:
  subroutine Hydrology2(lbc, ubc, lbp, ubp, &
                        num_nolakec, filter_nolakec, &
                        num_hydrologyc, filter_hydrologyc, &
                        num_urbanc, filter_urbanc, &
                        num_snowc, filter_snowc, &
                        num_nosnowc, filter_nosnowc)
!
! !DESCRIPTION:
! This is the main subroutine to execute the calculation of soil/snow
! hydrology
! Calling sequence is:
!  Hydrology2:                 surface hydrology driver
!    -> SnowWater:             change of snow mass and snow water onto soil
!    -> SurfaceRunoff:         surface runoff
!    -> Infiltration:          infiltration into surface soil layer
!    -> SoilWater:             soil water movement between layers
!          -> Tridiagonal      tridiagonal matrix solution
!    -> Drainage:              subsurface runoff
!    -> SnowCompaction:        compaction of snow layers
!    -> CombineSnowLayers:     combine snow layers that are thinner than minimum
!    -> DivideSnowLayers:      subdivide snow layers that are thicker than maximum
!
! !USES:
    use clmtype
    use clm_atmlnd      , only : clm_a2l
    use clm_varcon      , only : denh2o, denice, istice, istwet, istsoil, isturb, istice_mec, spval, &
                                 icol_roof, icol_road_imperv, icol_road_perv, icol_sunwall, &
                                 icol_shadewall, istdlak, &
                                 tfrz, hfus, grav
    use clm_varcon      , only : istcrop
    use clm_varctl      , only : glc_dyntopo
    use clm_varpar      , only : nlevgrnd, nlevsno, nlevsoi, nlevurb
    use SnowHydrologyMod, only : SnowCompaction, CombineSnowLayers, DivideSnowLayers, &
                                 SnowWater, BuildSnowFilter
    use SoilHydrologyMod, only : Infiltration, SoilWater, Drainage, SurfaceRunoff
    use clm_time_manager, only : get_step_size, get_nstep, is_perpetual
#if (defined VICHYDRO)
    use CLMVICMapMod    , only : CLMVICMap
#endif

!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                    ! column bounds
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer, intent(in) :: num_hydrologyc              ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(ubc-lbc+1)! column filter for soil points
    integer, intent(in) :: num_urbanc                  ! number of column urban points in column filter
    integer, intent(in) :: filter_urbanc(ubc-lbc+1)    ! column filter for urban points
    integer  :: num_snowc                  ! number of column snow points
    integer  :: filter_snowc(ubc-lbc+1)    ! column filter for snow points
    integer  :: num_nosnowc                ! number of column non-snow points
    integer  :: filter_nosnowc(ubc-lbc+1)  ! column filter for non-snow points
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    real(r8), pointer :: snow_depth(:)      !snow height of snow covered area (m)
    real(r8), pointer :: snowdp(:)             ! gridcell averaged snow height (m)
    real(r8), pointer :: frac_sno_eff(:)  !eff.  snow cover fraction (col) [frc]
    real(r8), pointer :: qflx_evap_soi(:) ! soil evaporation
    real(r8), pointer :: h2osfc(:)        ! surface water (mm)
    real(r8), pointer :: frac_h2osfc(:)   ! fraction of ground covered by surface water (0 to 1)
    real(r8), pointer :: t_h2osfc(:) 	  ! surface water temperature
    real(r8), pointer :: qflx_drain_perched(:)    ! sub-surface runoff from perched zwt (mm H2O /s)
    real(r8), pointer :: qflx_floodg(:)   ! gridcell flux of flood water from RTM
    real(r8), pointer :: qflx_h2osfc_surf(:)!surface water runoff (mm/s)
    logical , pointer :: cactive(:)       ! true=>do computations on this column (see reweightMod for details)
    integer , pointer :: cgridcell(:)     ! column's gridcell
    integer , pointer :: clandunit(:)     ! column's landunit
    integer , pointer :: ityplun(:)       ! landunit type
    integer , pointer :: ctype(:)         ! column type
    integer , pointer :: snl(:)           ! number of snow layers
    real(r8), pointer :: h2ocan(:)        ! canopy water (mm H2O)
    real(r8), pointer :: h2osno(:)        ! snow water (mm H2O)
    real(r8), pointer :: watsat(:,:)      ! volumetric soil water at saturation (porosity)
    real(r8), pointer :: sucsat(:,:)      ! minimum soil suction (mm)
    real(r8), pointer :: bsw(:,:)         ! Clapp and Hornberger "b"
    real(r8), pointer :: z(:,:)           ! layer depth  (m)
    real(r8), pointer :: forc_rain(:)     ! rain rate [mm/s]
    real(r8), pointer :: forc_snow(:)     ! snow rate [mm/s]
    real(r8), pointer :: begwb(:)         ! water mass begining of the time step
    real(r8), pointer :: qflx_evap_tot(:) ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(r8), pointer :: smpmin(:)        ! restriction for min of soil potential (mm)
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: dz(:,:)          ! layer thickness depth (m)
    real(r8), pointer :: zi(:,:)          ! interface depth (m)
    real(r8), pointer :: zwt(:)           ! water table depth (m)
    real(r8), pointer :: fcov(:)          ! fractional impermeable area
    real(r8), pointer :: fsat(:)          ! fractional area with water table at surface
    real(r8), pointer :: wa(:)            ! water in the unconfined aquifer (mm)
    real(r8), pointer :: qcharge(:)       ! aquifer recharge rate (mm/s)
    real(r8), pointer :: smp_l(:,:)       ! soil matrix potential [mm]
    real(r8), pointer :: hk_l(:,:)        ! hydraulic conductivity (mm/s)
    real(r8), pointer :: qflx_rsub_sat(:) ! soil saturation excess [mm h2o/s]
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: endwb(:)         ! water mass end of the time step
    real(r8), pointer :: wf(:)            ! soil water as frac. of whc for top 0.05 m  !F. Li and S. Levis
    real(r8), pointer :: wf2(:)           ! soil water as frac. of whc for top 0.17 m added by F. Li and S. Levis 
    real(r8), pointer :: snowice(:)       ! average snow ice lens
    real(r8), pointer :: snowliq(:)       ! average snow liquid water
    real(r8), pointer :: t_grnd(:)        ! ground temperature (Kelvin)
    real(r8), pointer :: t_soisno(:,:)    ! soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)  ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(r8), pointer :: t_soi_10cm(:)         ! soil temperature in top 10cm of soil (Kelvin)
    real(r8), pointer :: tsoi17(:)         ! soil temperature in top 17cm of soil (Kelvin) added by F. Li and S. Levis
   real(r8), pointer :: h2osoi_liqice_10cm(:) ! liquid water + ice lens in top 10cm of soil (kg/m2)
    real(r8), pointer :: h2osoi_vol(:,:)  ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(r8), pointer :: qflx_drain(:)    ! sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_surf(:)     ! surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_infl(:)     ! infiltration (mm H2O /s)
    real(r8), pointer :: qflx_qrgwl(:)    ! qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: qflx_irrig(:)    ! irrigation flux (mm H2O /s)
    real(r8), pointer :: qflx_runoff(:)   ! total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(r8), pointer :: qflx_runoff_u(:) ! Urban total runoff (qflx_drain+qflx_surf) (mm H2O /s)
    real(r8), pointer :: qflx_runoff_r(:) ! Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(r8), pointer :: t_grnd_u(:)      ! Urban ground temperature (Kelvin)
    real(r8), pointer :: t_grnd_r(:)      ! Rural ground temperature (Kelvin)
    real(r8), pointer :: qflx_snwcp_ice(:)! excess snowfall due to snow capping (mm H2O /s) [+]`
    real(r8), pointer :: soilpsi(:,:)     ! soil water potential in each soil layer (MPa)

    real(r8), pointer :: snot_top(:)        ! snow temperature in top layer (col) [K]
    real(r8), pointer :: dTdz_top(:)        ! temperature gradient in top layer (col) [K m-1]
    real(r8), pointer :: snw_rds(:,:)       ! effective snow grain radius (col,lyr) [microns, m^-6]
    real(r8), pointer :: snw_rds_top(:)     ! effective snow grain size, top layer(col) [microns]
    real(r8), pointer :: sno_liq_top(:)     ! liquid water fraction in top snow layer (col) [frc]
    real(r8), pointer :: frac_sno(:)        ! snow cover fraction (col) [frc]
    real(r8), pointer :: h2osno_top(:)      ! mass of snow in top layer (col) [kg]

    real(r8), pointer :: mss_bcpho(:,:)     ! mass of hydrophobic BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bcphi(:,:)     ! mass of hydrophillic BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bctot(:,:)     ! total mass of BC (pho+phi) (col,lyr) [kg]
    real(r8), pointer :: mss_bc_col(:)      ! total mass of BC in snow column (col) [kg]
    real(r8), pointer :: mss_bc_top(:)      ! total mass of BC in top snow layer (col) [kg]
    real(r8), pointer :: mss_cnc_bcphi(:,:) ! mass concentration of BC species 1 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_bcpho(:,:) ! mass concentration of BC species 2 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_ocpho(:,:)     ! mass of hydrophobic OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocphi(:,:)     ! mass of hydrophillic OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_octot(:,:)     ! total mass of OC (pho+phi) (col,lyr) [kg]
    real(r8), pointer :: mss_oc_col(:)      ! total mass of OC in snow column (col) [kg]
    real(r8), pointer :: mss_oc_top(:)      ! total mass of OC in top snow layer (col) [kg]
    real(r8), pointer :: mss_cnc_ocphi(:,:) ! mass concentration of OC species 1 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_ocpho(:,:) ! mass concentration of OC species 2 (col,lyr) [kg/kg]

    real(r8), pointer :: mss_dst1(:,:)      ! mass of dust species 1 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst2(:,:)      ! mass of dust species 2 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst3(:,:)      ! mass of dust species 3 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst4(:,:)      ! mass of dust species 4 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dsttot(:,:)    ! total mass of dust in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst_col(:)     ! total mass of dust in snow column (col) [kg]
    real(r8), pointer :: mss_dst_top(:)     ! total mass of dust in top snow layer (col) [kg]
    real(r8), pointer :: mss_cnc_dst1(:,:)  ! mass concentration of dust species 1 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst2(:,:)  ! mass concentration of dust species 2 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst3(:,:)  ! mass concentration of dust species 3 (col,lyr) [kg/kg]
    real(r8), pointer :: mss_cnc_dst4(:,:)  ! mass concentration of dust species 4 (col,lyr) [kg/kg]
    logical , pointer :: do_capsnow(:)      ! true => do snow capping
    real(r8), pointer :: qflx_glcice(:)     ! flux of new glacier ice (mm H2O /s)
    real(r8), pointer :: qflx_glcice_frz(:) ! ice growth (positive definite) (mm H2O/s)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: g,l,c,j,fc                 ! indices
    integer  :: nstep                      ! time step number
    real(r8) :: dtime                      ! land model time step (sec)
    real(r8) :: vol_liq(lbc:ubc,1:nlevgrnd)! partial volume of liquid water in layer
    real(r8) :: dwat(lbc:ubc,1:nlevgrnd)   ! change in soil water
    real(r8) :: hk(lbc:ubc,1:nlevgrnd)     ! hydraulic conductivity (mm h2o/s)
    real(r8) :: dhkdw(lbc:ubc,1:nlevgrnd)  ! d(hk)/d(vol_liq)
    real(r8) :: psi,vwc,fsattmp,psifrz     ! temporary variables for soilpsi calculation
#if (defined CN) 
    real(r8) :: watdry                     ! temporary
    real(r8) :: rwat(lbc:ubc)              ! soil water wgted by depth to maximum depth of 0.5 m
    real(r8) :: swat(lbc:ubc)              ! same as rwat but at saturation
    real(r8) :: rz(lbc:ubc)                ! thickness of soil layers contributing to rwat (m)
    real(r8) :: tsw                        ! volumetric soil water to 0.5 m
    real(r8) :: stsw                       ! volumetric soil water to 0.5 m at saturation
#endif
    real(r8) :: snowmass                   ! liquid+ice snow mass in a layer [kg/m2]
    real(r8) :: snowcap_scl_fct            ! temporary factor used to correct for snow capping
    real(r8) :: fracl                      ! fraction of soil layer contributing to 10cm total soil water
    real(r8) :: s_node                     ! soil wetness (-)
    real(r8) :: icefrac(lbc:ubc,1:nlevsoi)

!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (gridcell-level)

    forc_rain => clm_a2l%forc_rain
    forc_snow => clm_a2l%forc_snow

    ! Assign local pointers to derived subtypes components (landunit-level)

    ityplun => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    snow_depth      => clm3%g%l%c%cps%snow_depth
    snowdp             => clm3%g%l%c%cps%snowdp
    frac_sno_eff      => clm3%g%l%c%cps%frac_sno_eff 
    qflx_evap_soi     => clm3%g%l%c%cwf%pwf_a%qflx_evap_soi
    h2osfc            => clm3%g%l%c%cws%h2osfc
    frac_h2osfc       => clm3%g%l%c%cps%frac_h2osfc
    t_h2osfc          => clm3%g%l%c%ces%t_h2osfc
    qflx_drain_perched=> clm3%g%l%c%cwf%qflx_drain_perched
    qflx_floodg       => clm_a2l%forc_flood
    qflx_h2osfc_surf  => clm3%g%l%c%cwf%qflx_h2osfc_surf
    cactive           => clm3%g%l%c%active
    cgridcell         => clm3%g%l%c%gridcell
    clandunit         => clm3%g%l%c%landunit
    ctype             => clm3%g%l%c%itype
    snl               => clm3%g%l%c%cps%snl
    t_grnd            => clm3%g%l%c%ces%t_grnd
    h2ocan            => clm3%g%l%c%cws%pws_a%h2ocan
    h2osno            => clm3%g%l%c%cws%h2osno
    wf                => clm3%g%l%c%cps%wf
    wf2               => clm3%g%l%c%cps%wf2
    snowice           => clm3%g%l%c%cws%snowice
    snowliq           => clm3%g%l%c%cws%snowliq
    zwt               => clm3%g%l%c%cws%zwt
    fcov              => clm3%g%l%c%cws%fcov
    fsat              => clm3%g%l%c%cws%fsat
    wa                => clm3%g%l%c%cws%wa
    qcharge           => clm3%g%l%c%cws%qcharge
    watsat            => clm3%g%l%c%cps%watsat
    sucsat            => clm3%g%l%c%cps%sucsat
    bsw               => clm3%g%l%c%cps%bsw
    z                 => clm3%g%l%c%cps%z
    dz                => clm3%g%l%c%cps%dz
    zi                => clm3%g%l%c%cps%zi
    t_soisno          => clm3%g%l%c%ces%t_soisno
    h2osoi_ice        => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq        => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_vol        => clm3%g%l%c%cws%h2osoi_vol
    t_soi_10cm         => clm3%g%l%c%ces%t_soi_10cm
    tsoi17             => clm3%g%l%c%ces%tsoi17
    h2osoi_liqice_10cm => clm3%g%l%c%cws%h2osoi_liqice_10cm
    qflx_evap_tot     => clm3%g%l%c%cwf%pwf_a%qflx_evap_tot
    qflx_drain        => clm3%g%l%c%cwf%qflx_drain
    qflx_surf         => clm3%g%l%c%cwf%qflx_surf
    qflx_infl         => clm3%g%l%c%cwf%qflx_infl
    qflx_qrgwl        => clm3%g%l%c%cwf%qflx_qrgwl
    qflx_irrig        => clm3%g%l%c%cwf%qflx_irrig
    endwb             => clm3%g%l%c%cwbal%endwb
    begwb             => clm3%g%l%c%cwbal%begwb
    soilpsi           => clm3%g%l%c%cps%soilpsi
    smp_l             => clm3%g%l%c%cws%smp_l
    hk_l              => clm3%g%l%c%cws%hk_l
    qflx_rsub_sat     => clm3%g%l%c%cwf%qflx_rsub_sat
    qflx_runoff       => clm3%g%l%c%cwf%qflx_runoff
    qflx_runoff_u     => clm3%g%l%c%cwf%qflx_runoff_u
    qflx_runoff_r     => clm3%g%l%c%cwf%qflx_runoff_r
    t_grnd_u          => clm3%g%l%c%ces%t_grnd_u
    t_grnd_r          => clm3%g%l%c%ces%t_grnd_r
    snot_top          => clm3%g%l%c%cps%snot_top
    dTdz_top          => clm3%g%l%c%cps%dTdz_top
    snw_rds           => clm3%g%l%c%cps%snw_rds    
    snw_rds_top       => clm3%g%l%c%cps%snw_rds_top
    sno_liq_top       => clm3%g%l%c%cps%sno_liq_top
    frac_sno          => clm3%g%l%c%cps%frac_sno
    h2osno_top        => clm3%g%l%c%cps%h2osno_top
    mss_bcpho         => clm3%g%l%c%cps%mss_bcpho
    mss_bcphi         => clm3%g%l%c%cps%mss_bcphi
    mss_bctot         => clm3%g%l%c%cps%mss_bctot
    mss_bc_col        => clm3%g%l%c%cps%mss_bc_col
    mss_bc_top        => clm3%g%l%c%cps%mss_bc_top
    mss_cnc_bcphi     => clm3%g%l%c%cps%mss_cnc_bcphi
    mss_cnc_bcpho     => clm3%g%l%c%cps%mss_cnc_bcpho
    mss_ocpho         => clm3%g%l%c%cps%mss_ocpho
    mss_ocphi         => clm3%g%l%c%cps%mss_ocphi
    mss_octot         => clm3%g%l%c%cps%mss_octot
    mss_oc_col        => clm3%g%l%c%cps%mss_oc_col
    mss_oc_top        => clm3%g%l%c%cps%mss_oc_top
    mss_cnc_ocphi     => clm3%g%l%c%cps%mss_cnc_ocphi
    mss_cnc_ocpho     => clm3%g%l%c%cps%mss_cnc_ocpho
    mss_dst1          => clm3%g%l%c%cps%mss_dst1
    mss_dst2          => clm3%g%l%c%cps%mss_dst2
    mss_dst3          => clm3%g%l%c%cps%mss_dst3
    mss_dst4          => clm3%g%l%c%cps%mss_dst4
    mss_dsttot        => clm3%g%l%c%cps%mss_dsttot
    mss_dst_col       => clm3%g%l%c%cps%mss_dst_col
    mss_dst_top       => clm3%g%l%c%cps%mss_dst_top
    mss_cnc_dst1      => clm3%g%l%c%cps%mss_cnc_dst1
    mss_cnc_dst2      => clm3%g%l%c%cps%mss_cnc_dst2
    mss_cnc_dst3      => clm3%g%l%c%cps%mss_cnc_dst3
    mss_cnc_dst4      => clm3%g%l%c%cps%mss_cnc_dst4
    do_capsnow        => clm3%g%l%c%cps%do_capsnow
    qflx_snwcp_ice    => clm3%g%l%c%cwf%pwf_a%qflx_snwcp_ice
    qflx_glcice       => clm3%g%l%c%cwf%qflx_glcice
    smpmin            => clm3%g%l%c%cps%smpmin
    qflx_glcice_frz   => clm3%g%l%c%cwf%qflx_glcice_frz

    ! Determine time step and step size

    nstep = get_nstep()
    dtime = get_step_size()

    ! Determine initial snow/no-snow filters (will be modified possibly by
    ! routines CombineSnowLayers and DivideSnowLayers below

    call BuildSnowFilter(lbc, ubc, num_nolakec, filter_nolakec, &
         num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

    ! Determine the change of snow mass and the snow water onto soil

    call SnowWater(lbc, ubc, num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

    ! Determine soil hydrology
#if (defined VICHYDRO)
    ! mapping soilmoist from CLM to VIC layers for runoff calculations
    call CLMVICMap(lbc, ubc, num_hydrologyc, filter_hydrologyc)
#endif

    ! moved vol_liq from SurfaceRunoff to Infiltration
    call SurfaceRunoff(lbc, ubc, lbp, ubp, num_hydrologyc, filter_hydrologyc, &
                       num_urbanc, filter_urbanc, icefrac )

    call Infiltration(lbc, ubc,  num_hydrologyc, filter_hydrologyc, &
                      num_urbanc, filter_urbanc, vol_liq)

    call SoilWater(lbc, ubc, num_hydrologyc, filter_hydrologyc, &
                   num_urbanc, filter_urbanc, dwat, hk, dhkdw)

#if (defined VICHYDRO)
    ! mapping soilmoist from CLM to VIC layers for runoff calculations
    call CLMVICMap(lbc, ubc, num_hydrologyc, filter_hydrologyc)
#endif

    call Drainage(lbc, ubc, num_hydrologyc, filter_hydrologyc, &
                  num_urbanc, filter_urbanc, vol_liq, icefrac)

    if (.not. is_perpetual()) then

       ! Natural compaction and metamorphosis.

       call SnowCompaction(lbc, ubc, num_snowc, filter_snowc)

       ! Combine thin snow elements

       call CombineSnowLayers(lbc, ubc, num_snowc, filter_snowc)

       ! Divide thick snow elements

       call DivideSnowLayers(lbc, ubc, num_snowc, filter_snowc)

    else

       do fc = 1, num_snowc
          c = filter_snowc(fc)
          h2osno(c) = 0._r8
       end do
       do j = -nlevsno+1,0
          do fc = 1, num_snowc
             c = filter_snowc(fc)
             if (j >= snl(c)+1) then
                h2osno(c) = h2osno(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
             end if
          end do
       end do

    end if

    ! Set empty snow layers to zero

    do j = -nlevsno+1,0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
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

    call BuildSnowFilter(lbc, ubc, num_nolakec, filter_nolakec, &
         num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)

    ! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice
    ! over all snow layers for history output

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       snowice(c) = 0._r8
       snowliq(c) = 0._r8
    end do

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             snowice(c) = snowice(c) + h2osoi_ice(c,j)
             snowliq(c) = snowliq(c) + h2osoi_liq(c,j)
          end if
       end do
    end do

! Calculate column average snow depth
    do c = lbc,ubc
       snowdp(c) = snow_depth(c) * frac_sno_eff(c)
    end do

   ! Determine ground temperature, ending water balance and volumetric soil water
    ! Calculate soil temperature and total water (liq+ice) in top 10cm of soil
    ! Calculate soil temperature and total water (liq+ice) in top 17cm of soil
    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       if (ityplun(l) /= isturb) then
          t_soi_10cm(c) = 0._r8
          tsoi17(c) = 0._r8
          h2osoi_liqice_10cm(c) = 0._r8
       end if
    end do
    do j = 1, nlevsoi
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (ityplun(l) /= isturb) then
    ! soil T at top 17 cm added by F. Li and S. Levis
            if (zi(c,j) <= 0.17_r8) then
              fracl = 1._r8
              tsoi17(c) = tsoi17(c) + t_soisno(c,j)*dz(c,j)*fracl
            else
              if (zi(c,j) > 0.17_r8 .and. zi(c,j-1) .lt. 0.17_r8) then 
                fracl = (0.17_r8 - zi(c,j-1))/dz(c,j)
                tsoi17(c) = tsoi17(c) + t_soisno(c,j)*dz(c,j)*fracl
               end if
            end if

            if (zi(c,j) <= 0.1_r8) then
              fracl = 1._r8
              t_soi_10cm(c) = t_soi_10cm(c) + t_soisno(c,j)*dz(c,j)*fracl
              h2osoi_liqice_10cm(c) = h2osoi_liqice_10cm(c) + &
              (h2osoi_liq(c,j)+h2osoi_ice(c,j))* &
                                       fracl
            else
              if (zi(c,j) > 0.1_r8 .and. zi(c,j-1) .lt. 0.1_r8) then
                 fracl = (0.1_r8 - zi(c,j-1))/dz(c,j)
                 t_soi_10cm(c) = t_soi_10cm(c) + t_soisno(c,j)*dz(c,j)*fracl
                 h2osoi_liqice_10cm(c) = h2osoi_liqice_10cm(c) + &
                 (h2osoi_liq(c,j)+h2osoi_ice(c,j))* &
                                          fracl
              end if
            end if
          end if
       end do
    end do

    do fc = 1, num_nolakec
       
       c = filter_nolakec(fc)
       l = clandunit(c)

       ! t_grnd is weighted average of exposed soil and snow
       if (snl(c) < 0) then
          t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
               + (1 - frac_sno_eff(c)- frac_h2osfc(c)) * t_soisno(c,1) &
               + frac_h2osfc(c) * t_h2osfc(c)
       else
          t_grnd(c) = (1 - frac_h2osfc(c)) * t_soisno(c,1) + frac_h2osfc(c) * t_h2osfc(c)
       endif

       if (ityplun(l)==isturb) then
          t_grnd_u(c) = t_soisno(c,snl(c)+1)
       else
          t_soi_10cm(c) = t_soi_10cm(c)/0.1_r8
           tsoi17(c) =  tsoi17(c)/0.17_r8         ! F. Li and S. Levis
       end if
       if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
         t_grnd_r(c) = t_soisno(c,snl(c)+1)
       end if
       if (ctype(c) == icol_roof .or. ctype(c) == icol_sunwall &
          .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_road_imperv) then
         endwb(c) = h2ocan(c) + h2osno(c)
       else
          ! add h2osfc to water balance
          endwb(c) = h2ocan(c) + h2osno(c) + h2osfc(c) + wa(c)

       end if
    end do

    do j = 1, nlevgrnd
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
               .or. ctype(c) == icol_roof) .and. j > nlevurb) then
          else
            endwb(c) = endwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
            h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
          end if
       end do
    end do

    ! Determine wetland and land ice hydrology (must be placed here
    ! since need snow updated from CombineSnowLayers)

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       g = cgridcell(c)
       if (ityplun(l)==istwet .or. ityplun(l)==istice      &
                              .or. ityplun(l)==istice_mec) then
          qflx_drain(c)         = 0._r8
          qflx_drain_perched(c) = 0._r8
          qflx_h2osfc_surf(c)   = 0._r8
          qflx_irrig(c)         = 0._r8
          qflx_surf(c)          = 0._r8
          qflx_infl(c)          = 0._r8
          ! add flood water flux to runoff for wetlands/glaciers
          qflx_qrgwl(c) = forc_rain(g) + forc_snow(g) + qflx_floodg(g) - qflx_evap_tot(c) - qflx_snwcp_ice(c) - &
                          (endwb(c)-begwb(c))/dtime
          ! For dynamic topography, add meltwater from glacier_mec ice to the runoff.
          ! (Negative qflx_glcice => positive contribution to runoff)
          ! Note: The meltwater contribution is computed in PhaseChanges (part of Biogeophysics2).
          !       This code will not work if Hydrology2 is called before Biogeophysics2, or if
          !        qflx_snwcp_ice has alread been included in qflx_glcice.
          !       (The snwcp flux is added to qflx_glcice later in this subroutine.)

          if (glc_dyntopo .and. ityplun(l)==istice_mec) then
             qflx_qrgwl(c) = qflx_qrgwl(c) - qflx_glcice(c)   ! meltwater from melted ice
          endif
          fcov(c)       = spval
          fsat(c)       = spval
          qcharge(c)    = spval
          qflx_rsub_sat(c) = spval
       else if (ityplun(l) == isturb .and. ctype(c) /= icol_road_perv) then
          fcov(c)               = spval
          fsat(c)               = spval
          qflx_drain_perched(c) = 0._r8
          qflx_h2osfc_surf(c)   = 0._r8
          qcharge(c)            = spval
          qflx_rsub_sat(c)      = spval
       end if
       ! If snow exceeds the thickness limit in glacier_mec columns, convert to an ice flux.
       ! For dynamic glacier topography, remove qflx_snwcp_ice from the runoff.
       ! Note that qflx_glcice can also have a negative component from melting of bare ice,
       !  as computed in SoilTemperatureMod.F90

       if (ityplun(l)==istice_mec) then

          qflx_glcice_frz(c) = qflx_snwcp_ice(c)
          qflx_glcice(c) = qflx_glcice(c) + qflx_glcice_frz(c)

          ! For dynamic topography, set qflx_snwcp_ice = 0 so that this ice mass does not run off.
          ! For static topography, qflx_glc_ice is passed to the ice sheet model, but the 
          !  CLM runoff terms are not changed.
 
          if (glc_dyntopo) qflx_snwcp_ice(c) = 0._r8

       endif   ! istice_mec


       qflx_runoff(c) = qflx_drain(c) + qflx_surf(c)  + qflx_h2osfc_surf(c) + qflx_qrgwl(c) + qflx_drain_perched(c)

       if ((ityplun(l)==istsoil .or. ityplun(l)==istcrop) &
           .and. cactive(c)) then
          qflx_runoff(c) = qflx_runoff(c) - qflx_irrig(c)
       end if
       if (ityplun(l)==isturb) then
         qflx_runoff_u(c) = qflx_runoff(c)
       else if (ityplun(l)==istsoil .or. ityplun(l)==istcrop) then
         qflx_runoff_r(c) = qflx_runoff(c)
       end if

    end do

#if (defined CN) 
    ! Update soilpsi.
    ! ZMS: Note this could be merged with the following loop updating smp_l in the future.
    do j = 1, nlevgrnd
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          
          if (h2osoi_liq(c,j) > 0._r8) then

             vwc = h2osoi_liq(c,j)/(dz(c,j)*denh2o)
             
             ! the following limit set to catch very small values of 
             ! fractional saturation that can crash the calculation of psi

             ! use the same contants used in the supercool so that psi for frozen soils is consistent
             fsattmp = max(vwc/watsat(c,j), 0.001_r8)
             psi = sucsat(c,j) * (-9.8e-6_r8) * (fsattmp)**(-bsw(c,j))  ! Mpa
             soilpsi(c,j) = min(max(psi,-15.0_r8),0._r8)
             
          else 
             soilpsi(c,j) = -15.0_r8
          end if
       end do
    end do
#endif

    ! Update smp_l for history and for ch4Mod.
    ! ZMS: Note, this form, which seems to be the same as used in SoilWater, DOES NOT distinguish between
    ! ice and water volume, in contrast to the soilpsi calculation above. It won't be used in ch4Mod if
    ! t_soisno <= tfrz, though.
    do j = 1, nlevgrnd
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
          s_node = min(1.0_r8, s_node)

          smp_l(c,j) = -sucsat(c,j)*s_node**(-bsw(c,j))
          smp_l(c,j) = max(smpmin(c), smp_l(c,j))
       end do
    end do

#if (defined CN)
    ! Available soil water up to a depth of 0.05 m.
    ! Potentially available soil water (=whc) up to a depth of 0.05 m.
    ! Water content as fraction of whc up to a depth of 0.05 m.

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       rwat(c) = 0._r8
       swat(c) = 0._r8
       rz(c)   = 0._r8
    end do

    do j = 1, nlevgrnd
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          !if (z(c,j)+0.5_r8*dz(c,j) <= 0.5_r8) then
          if (z(c,j)+0.5_r8*dz(c,j) <= 0.05_r8) then
             watdry = watsat(c,j) * (316230._r8/sucsat(c,j)) ** (-1._r8/bsw(c,j))
             rwat(c) = rwat(c) + (h2osoi_vol(c,j)-watdry) * dz(c,j)
             swat(c) = swat(c) + (watsat(c,j)    -watdry) * dz(c,j)
             rz(c) = rz(c) + dz(c,j)
          end if
       end do
    end do

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if (rz(c) /= 0._r8) then
          tsw  = rwat(c)/rz(c)
          stsw = swat(c)/rz(c)
       else
          watdry = watsat(c,1) * (316230._r8/sucsat(c,1)) ** (-1._r8/bsw(c,1))
          tsw = h2osoi_vol(c,1) - watdry
          stsw = watsat(c,1) - watdry
       end if
       wf(c) = tsw/stsw
    end do
    
      do j = 1, nlevgrnd
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if (z(c,j)+0.5_r8*dz(c,j) <= 0.17_r8) then
             watdry = watsat(c,j) * (316230._r8/sucsat(c,j)) ** (-1._r8/bsw(c,j))
             rwat(c) = rwat(c) + (h2osoi_vol(c,j)-watdry) * dz(c,j)
             swat(c) = swat(c) + (watsat(c,j)    -watdry) * dz(c,j)
             rz(c) = rz(c) + dz(c,j)
          end if
       end do
    end do
 
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if (rz(c) /= 0._r8) then
          tsw  = rwat(c)/rz(c)
          stsw = swat(c)/rz(c)
       else
          watdry = watsat(c,1) * (316230._r8/sucsat(c,1)) ** (-1._r8/bsw(c,1))
          tsw = h2osoi_vol(c,1) - watdry
          stsw = watsat(c,1) - watdry
       end if
       wf2(c) = tsw/stsw
    end do
#endif

    !  Calculate column-integrated aerosol masses, and
    !  mass concentrations for radiative calculations and output
    !  (based on new snow level state, after SnowFilter is rebuilt.
    !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there 
    !  can be zero snow layers but an active column in filter)

    do fc = 1, num_snowc
       c = filter_snowc(fc)

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
                snowcap_scl_fct = snowmass / (snowmass+(qflx_snwcp_ice(c)*dtime))

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
    do fc = 1, num_nosnowc
       c = filter_nosnowc(fc)
            
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

  end subroutine Hydrology2

end module Hydrology2Mod
