module SurfaceAlbedoMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Performs surface albedo calculations
  !
  ! !PUBLIC TYPES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use abortutils        , only : endrun
  use decompMod         , only : bounds_type
  use landunit_varcon   , only : istsoil, istcrop, istdlak
  use elm_varcon        , only : grlnd, namep, namet
  use elm_varpar        , only : numrad, nlevcan, nlevsno, nlevcan
  use elm_varctl        , only : fsurdat, iulog, subgridflag, use_snicar_frc, use_fates, use_snicar_ad, use_top_solar_rad
  use VegetationPropertiesType    , only : veg_vp
  use SnowSnicarMod     , only : sno_nbr_aer, SNICAR_RT, SNICAR_AD_RT, DO_SNO_AER, DO_SNO_OC
  use AerosolType       , only : aerosol_type
  use CanopyStateType   , only : canopystate_type
  use LakeStateType     , only : lakestate_type
  use SurfaceAlbedoType , only : surfalb_type
  use GridcellType      , only : grc_pp
  use LandunitType      , only : lun_pp
  use ColumnType        , only : col_pp
  use ColumnDataType    , only : col_es, col_ws
  use VegetationType    , only : veg_pp
  use VegetationDataType, only : veg_es, veg_ws

  use SurfaceAlbedoType , only : lake_melt_icealb, alblak, alblakwi
  use SurfaceAlbedoType , only : albice, albsat, albdry, isoicol
  use elm_instMod , only : alm_fates
  use topounit_varcon   , only : max_topounits  
  use TopounitType      , only : top_pp
  !
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceAlbedo  ! Surface albedo and two-stream fluxes
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: SoilAlbedo    ! Determine ground surface albedo
  private :: TwoStream     ! Two-stream fluxes for canopy radiative transfer
  private :: Albedo_TOP_Adjustment     ! TOP solar radiation parameterization for non-vegetation

  !
  !
  ! Coefficient for calculating ice "fraction" for lake surface albedo
  ! From D. Mironov (2010) Boreal Env. Research
  real(r8), parameter :: calb = 95.6_r8
  !$acc declare copyin(calb)

contains


  !-----------------------------------------------------------------------
  subroutine SurfaceAlbedo(bounds,     &
        num_nourbanc, filter_nourbanc, &
        num_nourbanp, filter_nourbanp, &
        num_urbanc   , filter_urbanc,  &
        num_urbanp   , filter_urbanp,  &
        nextsw_cday  , declinp1,       &
        aerosol_vars, canopystate_vars, &
        lakestate_vars, surfalb_vars &
        )
    ! !DESCRIPTION:
    ! Surface albedo and two-stream fluxes
    ! Surface albedos. Also fluxes (per unit incoming direct and diffuse
    ! radiation) reflected, transmitted, and absorbed by vegetation.
    ! Calculate sunlit and shaded fluxes as described by
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
    ! a multi-layer canopy to calculate APAR profile
     !
    ! The calling sequence is:
    ! -> SurfaceAlbedo:  albedos for next time step
    !    -> SoilAlbedo:  soil/lake/glacier/wetland albedos
    !    -> SNICAR_RT:   snow albedos: direct beam (SNICAR)
    !       or
    !       SNICAR_AD_RT: snow albedos: direct beam (SNICAR-AD)
    !    -> SNICAR_RT:   snow albedos: diffuse (SNICAR)
    !       or
    !       SNICAR_AD_RT:   snow albedos: diffuse (SNICAR-AD)
    !    -> TwoStream:   absorbed, reflected, transmitted solar fluxes (vis dir,vis dif, nir dir, nir dif)
    !
    ! Note that this is called with the "inactive_and_active" version of the filters, because
    ! the variables computed here are needed over inactive points that might later become
    ! active (due to landuse change). Thus, this routine cannot depend on variables that are
    ! only computed over active points.
    !
    ! !USES:
      !$acc routine seq
    use elm_varctl         , only : iulog, subgridflag, use_snicar_frc, use_fates, use_snicar_ad, use_top_solar_rad
    use shr_orb_mod

    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds    ! bounds
    integer , intent(in) :: num_nourbanc       ! number of columns in non-urban filter
    integer , intent(in) :: filter_nourbanc(:) ! column filter for non-urban points
    integer , intent(in) :: num_nourbanp       ! number of patches in non-urban filter
    integer , intent(in) :: filter_nourbanp(:) ! patch filter for non-urban points
    integer , intent(in) :: num_urbanc         ! number of columns in urban filter
    integer , intent(in) :: filter_urbanc(:)   ! column filter for urban points
    integer , intent(in) :: num_urbanp         ! number of patches in urban filter
    integer , intent(in) :: filter_urbanp(:)   ! patch filter for rban points
    real(r8), intent(in) :: nextsw_cday        ! calendar day at Greenwich (1.00, ..., days/year)
    real(r8), intent(in) :: declinp1           ! declination angle (radians) for next time step
    type(aerosol_type)     , intent(in)    :: aerosol_vars
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(lakestate_type)   , intent(in)    :: lakestate_vars
    type(surfalb_type)     , intent(inout) :: surfalb_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: i                                                                         ! index for layers [idx]
    integer  :: aer                                                                       ! index for sno_nbr_aer
    real(r8) :: extkn                      ! nitrogen allocation coefficient
    integer  :: fp,fc,g,c,p,iv             ! indices
    integer  :: ib                         ! band index
    integer  :: ic                         ! 0=unit incoming direct; 1=unit incoming diffuse
    real(r8) :: dinc                       ! lai+sai increment for canopy layer
    real(r8) :: dincmax                    ! maximum lai+sai increment for canopy layer
    real(r8) :: dincmax_sum                ! cumulative sum of maximum lai+sai increment for canopy layer
    real(r8) :: laisum                     ! sum of canopy layer lai for error check
    real(r8) :: saisum                     ! sum of canopy layer sai for error check
    integer  :: flg_slr                                                                   ! flag for SNICAR (=1 if direct, =2 if diffuse)
    integer  :: flg_snw_ice                                                               ! flag for SNICAR (=1 when called from ELM, =2 when called from sea-ice)
    integer  :: num_vegsol                                                                ! number of vegetated patches where coszen>0
    integer  :: num_novegsol                                                              ! number of vegetated patches where coszen>0
    integer  :: filter_vegsol   (bounds%endp-bounds%begp+1)                               ! patch filter where vegetated and coszen>0
    integer  :: filter_novegsol (bounds%endp-bounds%begp+1)                               ! patch filter where vegetated and coszen>0
    real(r8) :: wl              (bounds%begp:bounds%endp)                                 ! fraction of LAI+SAI that is LAI
    real(r8) :: ws              (bounds%begp:bounds%endp)                                 ! fraction of LAI+SAI that is SAI
    real(r8) :: blai(bounds%begp:bounds%endp)              ! lai buried by snow: tlai - elai
    real(r8) :: bsai(bounds%begp:bounds%endp)              ! sai buried by snow: tsai - esai
    real(r8) :: coszen_gcell    (bounds%begg:bounds%endg)                                 ! cosine solar zenith angle for next time step (grc)
    real(r8) :: coszen_patch    (bounds%begp:bounds%endp)                                 ! cosine solar zenith angle for next time step (pft)
    real(r8) :: rho(bounds%begp:bounds%endp,numrad)        ! leaf/stem refl weighted by fraction LAI and SAI
    real(r8) :: tau(bounds%begp:bounds%endp,numrad)        ! leaf/stem tran weighted by fraction LAI and SAI
    real(r8) :: albsfc          (bounds%begc:bounds%endc,numrad)                          ! albedo of surface underneath snow (col,bnd)
    real(r8) :: albsnd(bounds%begc:bounds%endc,numrad)     ! snow albedo (direct)
    real(r8) :: albsni(bounds%begc:bounds%endc,numrad)     ! snow albedo (diffuse)
    real(r8) :: albsnd_pur      (bounds%begc:bounds%endc,numrad)                          ! direct pure snow albedo (radiative forcing)
    real(r8) :: albsni_pur      (bounds%begc:bounds%endc,numrad)                          ! diffuse pure snow albedo (radiative forcing)
    real(r8) :: albsnd_bc       (bounds%begc:bounds%endc,numrad)                          ! direct snow albedo without BC (radiative forcing)
    real(r8) :: albsni_bc       (bounds%begc:bounds%endc,numrad)                          ! diffuse snow albedo without BC (radiative forcing)
    real(r8) :: albsnd_oc       (bounds%begc:bounds%endc,numrad)                          ! direct snow albedo without OC (radiative forcing)
    real(r8) :: albsni_oc       (bounds%begc:bounds%endc,numrad)                          ! diffuse snow albedo without OC (radiative forcing)
    real(r8) :: albsnd_dst      (bounds%begc:bounds%endc,numrad)                          ! direct snow albedo without dust (radiative forcing)
    real(r8) :: albsni_dst      (bounds%begc:bounds%endc,numrad)                          ! diffuse snow albedo without dust (radiative forcing)
    real(r8) :: flx_absd_snw    (bounds%begc:bounds%endc,-nlevsno+1:1,numrad)             ! flux absorption factor for just snow (direct) [frc]
    real(r8) :: flx_absi_snw    (bounds%begc:bounds%endc,-nlevsno+1:1,numrad)             ! flux absorption factor for just snow (diffuse) [frc]
    real(r8) :: foo_snw         (bounds%begc:bounds%endc,-nlevsno+1:1,numrad)             ! dummy array for forcing calls
    real(r8) :: h2osno_liq      (bounds%begc:bounds%endc,-nlevsno+1:0)                    ! liquid snow content (col,lyr) [kg m-2]
    real(r8) :: h2osno_ice      (bounds%begc:bounds%endc,-nlevsno+1:0)                    ! ice content in snow (col,lyr) [kg m-2]
    integer  :: snw_rds_in      (bounds%begc:bounds%endc,-nlevsno+1:0)                    ! snow grain size sent to SNICAR (col,lyr) [microns]
    real(r8) :: mss_cnc_aer_in_frc_pur (bounds%begc:bounds%endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for forcing calculation (zero) (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_bc  (bounds%begc:bounds%endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for BC forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_oc  (bounds%begc:bounds%endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for OC forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_dst (bounds%begc:bounds%endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for dust forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_fdb     (bounds%begc:bounds%endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of all aerosol species for feedback calculation (col,lyr,aer) [kg kg-1]
    real(r8), parameter :: mpe = 1.e-06_r8                                                ! prevents overflow for division by zero
    integer , parameter :: nband =numrad                                                  ! number of solar radiation waveband classes
  !-----------------------------------------------------------------------

   associate(&
          rhol          =>    veg_vp%rhol                     , & ! Input:  [real(r8)  (:,:) ]  leaf reflectance: 1=vis, 2=nir
          rhos          =>    veg_vp%rhos                     , & ! Input:  [real(r8)  (:,:) ]  stem reflectance: 1=vis, 2=nir
          taul          =>    veg_vp%taul                     , & ! Input:  [real(r8)  (:,:) ]  leaf transmittance: 1=vis, 2=nir
          taus          =>    veg_vp%taus                     , & ! Input:  [real(r8)  (:,:) ]  stem transmittance: 1=vis, 2=nir

          tlai          =>    canopystate_vars%tlai_patch     , & ! Input:  [real(r8)  (:)   ]  one-sided leaf area index, no burying by snow
          tsai          =>    canopystate_vars%tsai_patch     , & ! Input:  [real(r8)  (:)   ]  one-sided stem area index, no burying by snow
          elai          =>    canopystate_vars%elai_patch     , & ! Input:  [real(r8)  (:)   ]  one-sided leaf area index with burying by snow
          esai          =>    canopystate_vars%esai_patch     , & ! Input:  [real(r8)  (:)   ]  one-sided stem area index with burying by snow

          frac_sno      =>    col_ws%frac_sno        , & ! Input:  [real(r8)  (:)   ]  fraction of ground covered by snow (0 to 1)
          h2osno        =>    col_ws%h2osno          , & ! Input:  [real(r8)  (:)   ]  snow water (mm H2O)
          h2osoi_liq    =>    col_ws%h2osoi_liq      , & ! Input:  [real(r8)  (:,:) ]  liquid water content (col,lyr) [kg/m2]
          h2osoi_ice    =>    col_ws%h2osoi_ice      , & ! Input:  [real(r8)  (:,:) ]  ice lens content (col,lyr) [kg/m2]
          snw_rds       =>    col_ws%snw_rds         , & ! Input:  [real(r8)  (:,:) ]  snow grain radius (col,lyr) [microns]

          mss_cnc_bcphi =>    aerosol_vars%mss_cnc_bcphi_col      , & ! Input:  [real(r8)  (:,:) ]  mass concentration of hydrophilic BC (col,lyr) [kg/kg]
          mss_cnc_bcpho =>    aerosol_vars%mss_cnc_bcpho_col      , & ! Input:  [real(r8)  (:,:) ]  mass concentration of hydrophobic BC (col,lyr) [kg/kg]
          mss_cnc_ocphi =>    aerosol_vars%mss_cnc_ocphi_col      , & ! Input:  [real(r8)  (:,:) ]  mass concentration of hydrophilic OC (col,lyr) [kg/kg]
          mss_cnc_ocpho =>    aerosol_vars%mss_cnc_ocpho_col      , & ! Input:  [real(r8)  (:,:) ]  mass concentration of hydrophobic OC (col,lyr) [kg/kg]
          mss_cnc_dst1  =>    aerosol_vars%mss_cnc_dst1_col       , & ! Input:  [real(r8)  (:,:) ]  mass concentration of dust aerosol species 1 (col,lyr) [kg/kg]
          mss_cnc_dst2  =>    aerosol_vars%mss_cnc_dst2_col       , & ! Input:  [real(r8)  (:,:) ]  mass concentration of dust aerosol species 2 (col,lyr) [kg/kg]
          mss_cnc_dst3  =>    aerosol_vars%mss_cnc_dst3_col       , & ! Input:  [real(r8)  (:,:) ]  mass concentration of dust aerosol species 3 (col,lyr) [kg/kg]
          mss_cnc_dst4  =>    aerosol_vars%mss_cnc_dst4_col       , & ! Input:  [real(r8)  (:,:) ]  mass concentration of dust aerosol species 4 (col,lyr) [kg/kg]

          fsun_z        =>    surfalb_vars%fsun_z_patch           , & ! Output:  [real(r8) (:,:) ]  sunlit fraction of canopy layer
          tlai_z        =>    surfalb_vars%tlai_z_patch           , & ! Output:  [real(r8) (:,:) ]  tlai increment for canopy layer
          tsai_z        =>    surfalb_vars%tsai_z_patch           , & ! Output:  [real(r8) (:,:) ]  tsai increment for canopy layer
          vcmaxcintsun  =>    surfalb_vars%vcmaxcintsun_patch     , & ! Output:  [real(r8) (:)   ]  leaf to canopy scaling coefficient, sunlit leaf vcmax
          vcmaxcintsha  =>    surfalb_vars%vcmaxcintsha_patch     , & ! Output:  [real(r8) (:)   ]  leaf to canopy scaling coefficient, shaded leaf vcmax
          ncan          =>    surfalb_vars%ncan_patch             , & ! Output:  [integer  (:)   ]  number of canopy layers
          nrad          =>    surfalb_vars%nrad_patch             , & ! Output:  [integer  (:)   ]  number of canopy layers, above snow for radiative transfer
          coszen_col    =>    surfalb_vars%coszen_col             , & ! Output:  [real(r8) (:)   ]  cosine of solar zenith angle
          albgrd        =>    surfalb_vars%albgrd_col             , & ! Output:  [real(r8) (:,:) ]  ground albedo (direct)
          albgri        =>    surfalb_vars%albgri_col             , & ! Output:  [real(r8) (:,:) ]  ground albedo (diffuse)
          albsod        =>    surfalb_vars%albsod_col             , & ! Output:  [real(r8) (:,:) ]  direct-beam soil albedo (col,bnd) [frc]
          albsoi        =>    surfalb_vars%albsoi_col             , & ! Output:  [real(r8) (:,:) ]  diffuse soil albedo (col,bnd) [frc]
          albgrd_pur    =>    surfalb_vars%albgrd_pur_col         , & ! Output:  [real(r8) (:,:) ]  pure snow ground albedo (direct)
          albgri_pur    =>    surfalb_vars%albgri_pur_col         , & ! Output:  [real(r8) (:,:) ]  pure snow ground albedo (diffuse)
          albgrd_bc     =>    surfalb_vars%albgrd_bc_col          , & ! Output:  [real(r8) (:,:) ]  ground albedo without BC (direct)
          albgri_bc     =>    surfalb_vars%albgri_bc_col          , & ! Output:  [real(r8) (:,:) ]  ground albedo without BC (diffuse)
          albgrd_oc     =>    surfalb_vars%albgrd_oc_col          , & ! Output:  [real(r8) (:,:) ]  ground albedo without OC (direct)
          albgri_oc     =>    surfalb_vars%albgri_oc_col          , & ! Output:  [real(r8) (:,:) ]  ground albedo without OC (diffuse)
          albgrd_dst    =>    surfalb_vars%albgrd_dst_col         , & ! Output:  [real(r8) (:,:) ]  ground albedo without dust (direct)
          albgri_dst    =>    surfalb_vars%albgri_dst_col         , & ! Output:  [real(r8) (:,:) ]  ground albedo without dust (diffuse)
          albsnd_hst    =>    surfalb_vars%albsnd_hst_col         , & ! Output:  [real(r8) (:,:) ]  snow albedo, direct, for history files (col,bnd) [frc]
          albsni_hst    =>    surfalb_vars%albsni_hst_col         , & ! Output:  [real(r8) (:,:) ]  snow ground albedo, diffuse, for history files (col,bnd) [frc]
          albd          =>    surfalb_vars%albd_patch             , & ! Output:  [real(r8) (:,:) ]  surface albedo (direct)
          albi          =>    surfalb_vars%albi_patch             , & ! Output:  [real(r8) (:,:) ]  surface albedo (diffuse)
          fabd          =>    surfalb_vars%fabd_patch             , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux
          fabd_sun      =>    surfalb_vars%fabd_sun_patch         , & ! Output:  [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit direct flux
          fabd_sha      =>    surfalb_vars%fabd_sha_patch         , & ! Output:  [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit direct flux
          fabi          =>    surfalb_vars%fabi_patch             , & ! Output:  [real(r8) (:,:) ]  flux absorbed by canopy per unit diffuse flux
          fabi_sun      =>    surfalb_vars%fabi_sun_patch         , & ! Output:  [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit diffuse flux
          fabi_sha      =>    surfalb_vars%fabi_sha_patch         , & ! Output:  [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit diffuse flux
          ftdd          =>    surfalb_vars%ftdd_patch             , & ! Output:  [real(r8) (:,:) ]  down direct flux below canopy per unit direct flux
          ftid          =>    surfalb_vars%ftid_patch             , & ! Output:  [real(r8) (:,:) ]  down diffuse flux below canopy per unit direct flux
          ftii          =>    surfalb_vars%ftii_patch             , & ! Output:  [real(r8) (:,:) ]  down diffuse flux below canopy per unit diffuse flux
          flx_absdv     =>    surfalb_vars%flx_absdv_col          , & ! Output:  [real(r8) (:,:) ]  direct flux absorption factor (col,lyr): VIS [frc]
          flx_absdn     =>    surfalb_vars%flx_absdn_col          , & ! Output:  [real(r8) (:,:) ]  direct flux absorption factor (col,lyr): NIR [frc]
          flx_absiv     =>    surfalb_vars%flx_absiv_col          , & ! Output:  [real(r8) (:,:) ]  diffuse flux absorption factor (col,lyr): VIS [frc]
          flx_absin     =>    surfalb_vars%flx_absin_col          , & ! Output:  [real(r8) (:,:) ]  diffuse flux absorption factor (col,lyr): NIR [frc]
          fabd_sun_z    =>    surfalb_vars%fabd_sun_z_patch       , & ! Output:  [real(r8) (:,:) ]  absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabd_sha_z    =>    surfalb_vars%fabd_sha_z_patch       , & ! Output:  [real(r8) (:,:) ]  absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabi_sun_z    =>    surfalb_vars%fabi_sun_z_patch       , & ! Output:  [real(r8) (:,:) ]  absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
          fabi_sha_z    =>    surfalb_vars%fabi_sha_z_patch         & ! Output:  [real(r8) (:,:) ]  absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
          )

    ! Cosine solar zenith angle for next time step

    do g = bounds%begg,bounds%endg
       coszen_gcell(g) = shr_orb_cosz (nextsw_cday, grc_pp%lat(g), grc_pp%lon(g), declinp1)
    end do
    do c = bounds%begc,bounds%endc
       g = col_pp%gridcell(c)
       coszen_col(c) = coszen_gcell(g)
    end do
    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)
       g = veg_pp%gridcell(p)
          coszen_patch(p) = coszen_gcell(g)
    end do

    ! Initialize output because solar radiation only done if coszen > 0

    do ib = 1, numrad
       do fc = 1,num_nourbanc
          c = filter_nourbanc(fc)
          albsod(c,ib)     = 0._r8
          albsoi(c,ib)     = 0._r8
          albgrd(c,ib)     = 0._r8
          albgri(c,ib)     = 0._r8
          albgrd_pur(c,ib) = 0._r8
          albgri_pur(c,ib) = 0._r8
          albgrd_bc(c,ib)  = 0._r8
          albgri_bc(c,ib)  = 0._r8
          albgrd_oc(c,ib)  = 0._r8
          albgri_oc(c,ib)  = 0._r8
          albgrd_dst(c,ib) = 0._r8
          albgri_dst(c,ib) = 0._r8
          do i=-nlevsno+1,1,1
             flx_absdv(c,i) = 0._r8
             flx_absdn(c,i) = 0._r8
             flx_absiv(c,i) = 0._r8
             flx_absin(c,i) = 0._r8
          enddo
       end do

       do fp = 1,num_nourbanp
          p = filter_nourbanp(fp)
          albd(p,ib) = 1._r8
          albi(p,ib) = 1._r8
          fabd(p,ib) = 0._r8
          fabd_sun(p,ib) = 0._r8
          fabd_sha(p,ib) = 0._r8
          fabi(p,ib) = 0._r8
          fabi_sun(p,ib) = 0._r8
          fabi_sha(p,ib) = 0._r8
          ftdd(p,ib) = 0._r8
          ftid(p,ib) = 0._r8
          ftii(p,ib) = 0._r8

       end do

    end do  ! end of numrad loop

    ! SoilAlbedo called before SNICAR_RT/SNICAR_AD_RT
    ! so that reflectance of soil beneath snow column is known
    ! ahead of time for snow RT calculation.

    ! Snow albedos
    ! Note that snow albedo routine will only compute nonzero snow albedos
    ! where h2osno> 0 and coszen > 0

    ! Ground surface albedos
    ! Note that ground albedo routine will only compute nonzero snow albedos
    ! where coszen > 0

       call SoilAlbedo(bounds, &
            num_nourbanc, filter_nourbanc, &
         coszen_col(bounds%begc:bounds%endc), &
         albsnd(bounds%begc:bounds%endc, :), &
            albsni(bounds%begc:bounds%endc, :), &
            lakestate_vars, surfalb_vars)

    ! set variables to pass to SNICAR.

    flg_snw_ice = 1   ! calling from ELM, not CSIM
    do c=bounds%begc,bounds%endc
       albsfc(c,:)     = albsoi(c,:)
       h2osno_liq(c,:) = h2osoi_liq(c,-nlevsno+1:0)
       h2osno_ice(c,:) = h2osoi_ice(c,-nlevsno+1:0)
       snw_rds_in(c,:) = nint(snw_rds(c,:))
    end do

    ! zero aerosol input arrays
    do aer = 1, sno_nbr_aer
       do i = -nlevsno+1, 0
          do c = bounds%begc, bounds%endc
             mss_cnc_aer_in_frc_pur(c,i,aer) = 0._r8
             mss_cnc_aer_in_frc_bc(c,i,aer)  = 0._r8
             mss_cnc_aer_in_frc_oc(c,i,aer)  = 0._r8
             mss_cnc_aer_in_frc_dst(c,i,aer) = 0._r8
             mss_cnc_aer_in_fdb(c,i,aer)     = 0._r8
          end do
       end do
    end do

    ! Set aerosol input arrays
    ! feedback input arrays have been zeroed
    ! set soot and dust aerosol concentrations:
    if (DO_SNO_AER) then
       mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,1) = mss_cnc_bcphi(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,2) = mss_cnc_bcpho(bounds%begc:bounds%endc,:)

       ! DO_SNO_OC is set in SNICAR_varpar. Default case is to ignore OC concentrations because:
       !  1) Knowledge of their optical properties is primitive
       !  2) When 'water-soluble' OPAC optical properties are applied to OC in snow,
       !     it has a negligible darkening effect.
       if (DO_SNO_OC) then
          mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,3) = mss_cnc_ocphi(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,4) = mss_cnc_ocpho(bounds%begc:bounds%endc,:)
       endif

       mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,5) = mss_cnc_dst1(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,6) = mss_cnc_dst2(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,7) = mss_cnc_dst3(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_fdb(bounds%begc:bounds%endc,:,8) = mss_cnc_dst4(bounds%begc:bounds%endc,:)
    endif

    ! If radiative forcing is being calculated, first estimate clean-snow albedo

    if (use_snicar_frc) then
       ! 1. BC input array:
       !  set dust and (optionally) OC concentrations, so BC_FRC=[(BC+OC+dust)-(OC+dust)]
       mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,5) = mss_cnc_dst1(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,6) = mss_cnc_dst2(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,7) = mss_cnc_dst3(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,8) = mss_cnc_dst4(bounds%begc:bounds%endc,:)
       if (DO_SNO_OC) then
          mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,3) = mss_cnc_ocphi(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc,:,4) = mss_cnc_ocpho(bounds%begc:bounds%endc,:)
       endif

       ! BC FORCING CALCULATIONS
          flg_slr = 1; ! direct-beam
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                                coszen_col(bounds%begc:bounds%endc), &
                                flg_slr, &
                                h2osno_liq(bounds%begc:bounds%endc, :), &
                                h2osno_ice(bounds%begc:bounds%endc, :), &
                                snw_rds_in(bounds%begc:bounds%endc, :), &
                                mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc, :, :), &
                                albsfc(bounds%begc:bounds%endc, :), &
                                albsnd_bc(bounds%begc:bounds%endc, :), &
                                foo_snw(bounds%begc:bounds%endc, :, :) )
          else
                call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                             coszen_col(bounds%begc:bounds%endc), &
                             flg_slr, &
                             h2osno_liq(bounds%begc:bounds%endc, :), &
                             h2osno_ice(bounds%begc:bounds%endc, :), &
                             snw_rds_in(bounds%begc:bounds%endc, :), &
                             mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc, :, :), &
                             albsfc(bounds%begc:bounds%endc, :), &
                             albsnd_bc(bounds%begc:bounds%endc, :), &
                             foo_snw(bounds%begc:bounds%endc, :, :) )

          endif ! end if use_snicar_ad

          flg_slr = 2; ! diffuse
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                                coszen_col(bounds%begc:bounds%endc), &
                                flg_slr, &
                                h2osno_liq(bounds%begc:bounds%endc, :), &
                                h2osno_ice(bounds%begc:bounds%endc, :), &
                                snw_rds_in(bounds%begc:bounds%endc, :), &
                                mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc, :, :), &
                                albsfc(bounds%begc:bounds%endc, :), &
                                albsni_bc(bounds%begc:bounds%endc, :), &
                                foo_snw(bounds%begc:bounds%endc, :, :) )
          else
              call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                             coszen_col(bounds%begc:bounds%endc), &
                             flg_slr, &
                             h2osno_liq(bounds%begc:bounds%endc, :), &
                             h2osno_ice(bounds%begc:bounds%endc, :), &
                             snw_rds_in(bounds%begc:bounds%endc, :), &
                             mss_cnc_aer_in_frc_bc(bounds%begc:bounds%endc, :, :), &
                             albsfc(bounds%begc:bounds%endc, :), &
                             albsni_bc(bounds%begc:bounds%endc, :), &
                             foo_snw(bounds%begc:bounds%endc, :, :) )
          endif ! end if use_snicar_ad

       ! 2. OC input array:
       !  set BC and dust concentrations, so OC_FRC=[(BC+OC+dust)-(BC+dust)]
       if (DO_SNO_OC) then
          mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,1) = mss_cnc_bcphi(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,2) = mss_cnc_bcpho(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,5) = mss_cnc_dst1(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,6) = mss_cnc_dst2(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,7) = mss_cnc_dst3(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc,:,8) = mss_cnc_dst4(bounds%begc:bounds%endc,:)

       ! OC FORCING CALCULATIONS
          flg_slr = 1; ! direct-beam
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                                coszen_col(bounds%begc:bounds%endc), &
                                flg_slr, &
                                h2osno_liq(bounds%begc:bounds%endc, :), &
                                h2osno_ice(bounds%begc:bounds%endc, :), &
                                snw_rds_in(bounds%begc:bounds%endc, :), &
                                mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc, :, :), &
                                albsfc(bounds%begc:bounds%endc, :), &
                                albsnd_oc(bounds%begc:bounds%endc, :), &
                                foo_snw(bounds%begc:bounds%endc, :, :) )
          else
              call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                             coszen_col(bounds%begc:bounds%endc), &
                             flg_slr, &
                             h2osno_liq(bounds%begc:bounds%endc, :), &
                             h2osno_ice(bounds%begc:bounds%endc, :), &
                             snw_rds_in(bounds%begc:bounds%endc, :), &
                             mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc, :, :), &
                             albsfc(bounds%begc:bounds%endc, :), &
                             albsnd_oc(bounds%begc:bounds%endc, :), &
                             foo_snw(bounds%begc:bounds%endc, :, :) )
          endif ! end if use_snicar_ad

          flg_slr = 2; ! diffuse
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                                coszen_col(bounds%begc:bounds%endc), &
                                flg_slr, &
                                h2osno_liq(bounds%begc:bounds%endc, :), &
                                h2osno_ice(bounds%begc:bounds%endc, :), &
                                snw_rds_in(bounds%begc:bounds%endc, :), &
                                mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc, :, :), &
                                albsfc(bounds%begc:bounds%endc, :), &
                                albsni_oc(bounds%begc:bounds%endc, :), &
                                foo_snw(bounds%begc:bounds%endc, :, :) )
          else
              call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                             coszen_col(bounds%begc:bounds%endc), &
                             flg_slr, &
                             h2osno_liq(bounds%begc:bounds%endc, :), &
                             h2osno_ice(bounds%begc:bounds%endc, :), &
                             snw_rds_in(bounds%begc:bounds%endc, :), &
                             mss_cnc_aer_in_frc_oc(bounds%begc:bounds%endc, :, :), &
                             albsfc(bounds%begc:bounds%endc, :), &
                             albsni_oc(bounds%begc:bounds%endc, :), &
                             foo_snw(bounds%begc:bounds%endc, :, :) )
          endif ! end if use_snicar_ad
       endif

       ! 3. DUST input array:
       ! set BC and OC concentrations, so DST_FRC=[(BC+OC+dust)-(BC+OC)]
       mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc,:,1) = mss_cnc_bcphi(bounds%begc:bounds%endc,:)
       mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc,:,2) = mss_cnc_bcpho(bounds%begc:bounds%endc,:)
       if (DO_SNO_OC) then
          mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc,:,3) = mss_cnc_ocphi(bounds%begc:bounds%endc,:)
          mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc,:,4) = mss_cnc_ocpho(bounds%begc:bounds%endc,:)
       endif

       ! DUST FORCING CALCULATIONS
          flg_slr = 1; ! direct-beam
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                                coszen_col(bounds%begc:bounds%endc), &
                                flg_slr, &
                                h2osno_liq(bounds%begc:bounds%endc, :), &
                                h2osno_ice(bounds%begc:bounds%endc, :), &
                                snw_rds_in(bounds%begc:bounds%endc, :), &
                                mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc, :, :), &
                                albsfc(bounds%begc:bounds%endc, :), &
                                albsnd_dst(bounds%begc:bounds%endc, :), &
                                foo_snw(bounds%begc:bounds%endc, :, :) )
          else
              call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                             coszen_col(bounds%begc:bounds%endc), &
                             flg_slr, &
                             h2osno_liq(bounds%begc:bounds%endc, :), &
                             h2osno_ice(bounds%begc:bounds%endc, :), &
                             snw_rds_in(bounds%begc:bounds%endc, :), &
                             mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc, :, :), &
                             albsfc(bounds%begc:bounds%endc, :), &
                             albsnd_dst(bounds%begc:bounds%endc, :), &
                             foo_snw(bounds%begc:bounds%endc, :, :) )
          endif ! end if use_snicar_ad

          flg_slr = 2; ! diffuse
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                                coszen_col(bounds%begc:bounds%endc), &
                                flg_slr, &
                                h2osno_liq(bounds%begc:bounds%endc, :), &
                                h2osno_ice(bounds%begc:bounds%endc, :), &
                                snw_rds_in(bounds%begc:bounds%endc, :), &
                                mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc, :, :), &
                                albsfc(bounds%begc:bounds%endc, :), &
                                albsni_dst(bounds%begc:bounds%endc, :), &
                                foo_snw(bounds%begc:bounds%endc, :, :) )
          else
              call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                             coszen_col(bounds%begc:bounds%endc), &
                             flg_slr, &
                             h2osno_liq(bounds%begc:bounds%endc, :), &
                             h2osno_ice(bounds%begc:bounds%endc, :), &
                             snw_rds_in(bounds%begc:bounds%endc, :), &
                             mss_cnc_aer_in_frc_dst(bounds%begc:bounds%endc, :, :), &
                             albsfc(bounds%begc:bounds%endc, :), &
                             albsni_dst(bounds%begc:bounds%endc, :), &
                             foo_snw(bounds%begc:bounds%endc, :, :)  )
          endif ! end if use_snicar_ad

       ! 4. ALL AEROSOL FORCING CALCULATION
       ! (pure snow albedo)
          flg_slr = 1; ! direct-beam
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                                coszen_col(bounds%begc:bounds%endc), &
                                flg_slr, &
                                h2osno_liq(bounds%begc:bounds%endc, :), &
                                h2osno_ice(bounds%begc:bounds%endc, :), &
                                snw_rds_in(bounds%begc:bounds%endc, :), &
                                mss_cnc_aer_in_frc_pur(bounds%begc:bounds%endc, :, :), &
                                albsfc(bounds%begc:bounds%endc, :), &
                                albsnd_pur(bounds%begc:bounds%endc, :), &
                                foo_snw(bounds%begc:bounds%endc, :, :) )
          else
              call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                             coszen_col(bounds%begc:bounds%endc), &
                             flg_slr, &
                             h2osno_liq(bounds%begc:bounds%endc, :), &
                             h2osno_ice(bounds%begc:bounds%endc, :), &
                             snw_rds_in(bounds%begc:bounds%endc, :), &
                             mss_cnc_aer_in_frc_pur(bounds%begc:bounds%endc, :, :), &
                             albsfc(bounds%begc:bounds%endc, :), &
                             albsnd_pur(bounds%begc:bounds%endc, :), &
                             foo_snw(bounds%begc:bounds%endc, :, :) )
          endif ! end if use_snicar_ad

          flg_slr = 2; ! diffuse
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                                coszen_col(bounds%begc:bounds%endc), &
                                flg_slr, &
                                h2osno_liq(bounds%begc:bounds%endc, :), &
                                h2osno_ice(bounds%begc:bounds%endc, :), &
                                snw_rds_in(bounds%begc:bounds%endc, :), &
                                mss_cnc_aer_in_frc_pur(bounds%begc:bounds%endc, :, :), &
                                albsfc(bounds%begc:bounds%endc, :), &
                                albsni_pur(bounds%begc:bounds%endc, :), &
                                foo_snw(bounds%begc:bounds%endc, :, :) )
          else
              call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                             coszen_col(bounds%begc:bounds%endc), &
                             flg_slr, &
                             h2osno_liq(bounds%begc:bounds%endc, :), &
                             h2osno_ice(bounds%begc:bounds%endc, :), &
                             snw_rds_in(bounds%begc:bounds%endc, :), &
                             mss_cnc_aer_in_frc_pur(bounds%begc:bounds%endc, :, :), &
                             albsfc(bounds%begc:bounds%endc, :), &
                             albsni_pur(bounds%begc:bounds%endc, :), &
                             foo_snw(bounds%begc:bounds%endc, :, :) )
          endif ! end if use_snicar_ad
    end if !end if use_snicar_frc


    ! CLIMATE FEEDBACK CALCULATIONS, ALL AEROSOLS:
       flg_slr = 1; ! direct-beam
       if (use_snicar_ad) then
           call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                             coszen_col(bounds%begc:bounds%endc), &
                             flg_slr, &
                             h2osno_liq(bounds%begc:bounds%endc, :), &
                             h2osno_ice(bounds%begc:bounds%endc, :), &
                             snw_rds_in(bounds%begc:bounds%endc, :), &
                             mss_cnc_aer_in_fdb(bounds%begc:bounds%endc, :, :), &
                             albsfc(bounds%begc:bounds%endc, :), &
                             albsnd(bounds%begc:bounds%endc, :), &
                             flx_absd_snw(bounds%begc:bounds%endc, :, :) )
       else
            call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                           coszen_col(bounds%begc:bounds%endc), &
                           flg_slr, &
                           h2osno_liq(bounds%begc:bounds%endc, :), &
                           h2osno_ice(bounds%begc:bounds%endc, :), &
                           snw_rds_in(bounds%begc:bounds%endc, :), &
                           mss_cnc_aer_in_fdb(bounds%begc:bounds%endc, :, :), &
                           albsfc(bounds%begc:bounds%endc, :), &
                           albsnd(bounds%begc:bounds%endc, :), &
                           flx_absd_snw(bounds%begc:bounds%endc, :, :) )
       endif ! end if use_snicar_ad

       flg_slr = 2; ! diffuse
       if (use_snicar_ad) then
           call SNICAR_AD_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                             coszen_col(bounds%begc:bounds%endc), &
                             flg_slr, &
                             h2osno_liq(bounds%begc:bounds%endc, :), &
                             h2osno_ice(bounds%begc:bounds%endc, :), &
                             snw_rds_in(bounds%begc:bounds%endc, :), &
                             mss_cnc_aer_in_fdb(bounds%begc:bounds%endc, :, :), &
                             albsfc(bounds%begc:bounds%endc, :), &
                             albsni(bounds%begc:bounds%endc, :), &
                             flx_absi_snw(bounds%begc:bounds%endc, :, :) )
       else
            call SNICAR_RT(flg_snw_ice, bounds, num_nourbanc, filter_nourbanc,    &
                           coszen_col(bounds%begc:bounds%endc), &
                           flg_slr, &
                           h2osno_liq(bounds%begc:bounds%endc, :), &
                           h2osno_ice(bounds%begc:bounds%endc, :), &
                           snw_rds_in(bounds%begc:bounds%endc, :), &
                           mss_cnc_aer_in_fdb(bounds%begc:bounds%endc, :, :), &
                           albsfc(bounds%begc:bounds%endc, :), &
                           albsni(bounds%begc:bounds%endc, :), &
                           flx_absi_snw(bounds%begc:bounds%endc, :, :) )
       endif ! end if use_snicar_ad

    ! ground albedos and snow-fraction weighting of snow absorption factors
    do ib = 1, nband
       do fc = 1,num_nourbanc
          c = filter_nourbanc(fc)
             if (coszen_col(c) > 0._r8) then
             ! ground albedo was originally computed in SoilAlbedo, but is now computed here
             ! because the order of SoilAlbedo and SNICAR_RT/SNICAR_AD_RT was switched for SNICAR/SNICAR_AD_RT.
             albgrd(c,ib) = albsod(c,ib)*(1._r8-frac_sno(c)) + albsnd(c,ib)*frac_sno(c)
             albgri(c,ib) = albsoi(c,ib)*(1._r8-frac_sno(c)) + albsni(c,ib)*frac_sno(c)

             ! albedos for radiative forcing calculations:
             if (use_snicar_frc) then
                ! BC forcing albedo
                albgrd_bc(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_bc(c,ib)*frac_sno(c)
                albgri_bc(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_bc(c,ib)*frac_sno(c)

                if (DO_SNO_OC) then
                   ! OC forcing albedo
                   albgrd_oc(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_oc(c,ib)*frac_sno(c)
                   albgri_oc(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_oc(c,ib)*frac_sno(c)
                endif

                ! dust forcing albedo
                albgrd_dst(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_dst(c,ib)*frac_sno(c)
                albgri_dst(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_dst(c,ib)*frac_sno(c)

                ! pure snow albedo for all-aerosol radiative forcing
                albgrd_pur(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_pur(c,ib)*frac_sno(c)
                albgri_pur(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_pur(c,ib)*frac_sno(c)
             end if

             ! also in this loop (but optionally in a different loop for vectorized code)
             !  weight snow layer radiative absorption factors based on snow fraction and soil albedo
             !  (NEEDED FOR ENERGY CONSERVATION)
             do i = -nlevsno+1,1,1
              if (subgridflag == 0 .or. lun_pp%itype(col_pp%landunit(c)) == istdlak) then
                if (ib == 1) then
                   flx_absdv(c,i) = flx_absd_snw(c,i,ib)*frac_sno(c) + &
                        ((1.-frac_sno(c))*(1-albsod(c,ib))*(flx_absd_snw(c,i,ib)/(1.-albsnd(c,ib))))
                   flx_absiv(c,i) = flx_absi_snw(c,i,ib)*frac_sno(c) + &
                        ((1.-frac_sno(c))*(1-albsoi(c,ib))*(flx_absi_snw(c,i,ib)/(1.-albsni(c,ib))))
                elseif (ib == 2) then
                   flx_absdn(c,i) = flx_absd_snw(c,i,ib)*frac_sno(c) + &
                        ((1.-frac_sno(c))*(1-albsod(c,ib))*(flx_absd_snw(c,i,ib)/(1.-albsnd(c,ib))))
                   flx_absin(c,i) = flx_absi_snw(c,i,ib)*frac_sno(c) + &
                        ((1.-frac_sno(c))*(1-albsoi(c,ib))*(flx_absi_snw(c,i,ib)/(1.-albsni(c,ib))))
                endif
             else
                if (ib == 1) then
                   flx_absdv(c,i) = flx_absd_snw(c,i,ib)*(1.-albsnd(c,ib))
                   flx_absiv(c,i) = flx_absi_snw(c,i,ib)*(1.-albsni(c,ib))
                elseif (ib == 2) then
                   flx_absdn(c,i) = flx_absd_snw(c,i,ib)*(1.-albsnd(c,ib))
                   flx_absin(c,i) = flx_absi_snw(c,i,ib)*(1.-albsni(c,ib))
                endif
             endif
             enddo
          endif
       enddo
    enddo

       ! For diagnostics, set snow albedo to spval over non-snow non-urban points
       ! so that it is not averaged in history buffer (OPTIONAL)
       ! TODO - this is set to 0 not spval - seems wrong since it will be averaged in

    do ib = 1, nband
       do fc = 1,num_nourbanc
          c = filter_nourbanc(fc)
             if ((coszen_col(c) > 0._r8) .and. (h2osno(c) > 0._r8)) then
             albsnd_hst(c,ib) = albsnd(c,ib)
             albsni_hst(c,ib) = albsni(c,ib)
          else
             albsnd_hst(c,ib) = 0._r8
             albsni_hst(c,ib) = 0._r8
          endif
       enddo
    enddo

    ! Create solar-vegetated filter for the following calculations

    num_vegsol = 0
    num_novegsol = 0
    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)
          if (coszen_patch(p) > 0._r8) then
             if ((lun_pp%itype(veg_pp%landunit(p)) == istsoil .or.  &
                  lun_pp%itype(veg_pp%landunit(p)) == istcrop     ) &
                 .and. (elai(p) + esai(p)) > 0._r8) then
                    num_vegsol = num_vegsol + 1
                    filter_vegsol(num_vegsol) = p
             else
                num_novegsol = num_novegsol + 1
                filter_novegsol(num_novegsol) = p
             end if
          end if
    end do

    ! Weight reflectance/transmittance by lai and sai
       ! Only perform on vegetated patches where coszen > 0

    do fp = 1,num_vegsol
       p = filter_vegsol(fp)
       wl(p) = elai(p) / max( elai(p)+esai(p), mpe )
       ws(p) = esai(p) / max( elai(p)+esai(p), mpe )
    end do

    do ib = 1, numrad
       do fp = 1,num_vegsol
          p = filter_vegsol(fp)
          rho(p,ib) = max( rhol(veg_pp%itype(p),ib)*wl(p) + rhos(veg_pp%itype(p),ib)*ws(p), mpe )
          tau(p,ib) = max( taul(veg_pp%itype(p),ib)*wl(p) + taus(veg_pp%itype(p),ib)*ws(p), mpe )
       end do
    end do

    ! Diagnose number of canopy layers for radiative transfer, in increments of dincmax.
    ! Add to number of layers so long as cumulative leaf+stem area does not exceed total
    ! leaf+stem area. Then add any remaining leaf+stem area to next layer and exit the loop.
    ! Do this first for elai and esai (not buried by snow) and then for the part of the
    ! canopy that is buried by snow.
    ! ------------------
    ! tlai_z = leaf area increment for a layer
    ! tsai_z = stem area increment for a layer
    ! nrad   = number of canopy layers above snow
    ! ncan   = total number of canopy layers
    !
    ! tlai_z summed from 1 to nrad = elai
    ! tlai_z summed from 1 to ncan = tlai

    ! tsai_z summed from 1 to nrad = esai
    ! tsai_z summed from 1 to ncan = tsai
    ! ------------------
    !
    ! Canopy layering needs to be done for all "num_nourbanp" not "num_vegsol"
    ! because layering is needed for all time steps regardless of radiation
    !
    ! Sun/shade big leaf code uses only one layer (nrad = ncan = 1), triggered by
    ! nlevcan = 1

    dincmax = 0.25_r8
    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)

       if (nlevcan == 1) then
          nrad(p) = 1
          ncan(p) = 1
          tlai_z(p,1) = elai(p)
          tsai_z(p,1) = esai(p)
       else if (nlevcan > 1) then
          if (elai(p)+esai(p) == 0._r8) then
             nrad(p) = 0
          else
             dincmax_sum = 0._r8
             do iv = 1, nlevcan
                dincmax_sum = dincmax_sum + dincmax
                if (((elai(p)+esai(p))-dincmax_sum) > 1.e-06_r8) then
                   nrad(p) = iv
                   dinc = dincmax
                   tlai_z(p,iv) = dinc * elai(p) / max(elai(p)+esai(p), mpe)
                   tsai_z(p,iv) = dinc * esai(p) / max(elai(p)+esai(p), mpe)
                else
                   nrad(p) = iv
                   dinc = dincmax - (dincmax_sum - (elai(p)+esai(p)))
                   tlai_z(p,iv) = dinc * elai(p) / max(elai(p)+esai(p), mpe)
                   tsai_z(p,iv) = dinc * esai(p) / max(elai(p)+esai(p), mpe)
                   exit
                end if
             end do

             ! Mimumum of 4 canopy layers

             if (nrad(p) < 4) then
                nrad(p) = 4
                do iv = 1, nrad(p)
                   tlai_z(p,iv) = elai(p) / nrad(p)
                   tsai_z(p,iv) = esai(p) / nrad(p)
                end do
             end if
          end if
       end if

       ! Error check: make sure cumulative of increments does not exceed total

       laisum = 0._r8
       saisum = 0._r8
       do iv = 1, nrad(p)
          laisum = laisum + tlai_z(p,iv)
          saisum = saisum + tsai_z(p,iv)
       end do
       if (abs(laisum-elai(p)) > 1.e-06_r8 .or. abs(saisum-esai(p)) > 1.e-06_r8) then
#ifndef _OPENACC
          write (iulog,*) 'multi-layer canopy error 01 in SurfaceAlbedo: ',&
               nrad(p),elai(p),laisum,esai(p),saisum
          call endrun(decomp_index=p, elmlevel=namep, msg=errmsg(__FILE__, __LINE__))
#endif
       end if

       ! Repeat to find canopy layers buried by snow

       if (nlevcan > 1) then
          blai(p) = tlai(p) - elai(p)
          bsai(p) = tsai(p) - esai(p)
          if (blai(p)+bsai(p) == 0._r8) then
             ncan(p) = nrad(p)
          else
             dincmax_sum = 0._r8
             do iv = nrad(p)+1, nlevcan
                dincmax_sum = dincmax_sum + dincmax
                if (((blai(p)+bsai(p))-dincmax_sum) > 1.e-06_r8) then
                   ncan(p) = iv
                   dinc = dincmax
                   tlai_z(p,iv) = dinc * blai(p) / max(blai(p)+bsai(p), mpe)
                   tsai_z(p,iv) = dinc * bsai(p) / max(blai(p)+bsai(p), mpe)
                else
                   ncan(p) = iv
                   dinc = dincmax - (dincmax_sum - (blai(p)+bsai(p)))
                   tlai_z(p,iv) = dinc * blai(p) / max(blai(p)+bsai(p), mpe)
                   tsai_z(p,iv) = dinc * bsai(p) / max(blai(p)+bsai(p), mpe)
                   exit
                end if
             end do
          end if

          ! Error check: make sure cumulative of increments does not exceed total

          laisum = 0._r8
          saisum = 0._r8
          do iv = 1, ncan(p)
             laisum = laisum + tlai_z(p,iv)
             saisum = saisum + tsai_z(p,iv)
          end do
          if (abs(laisum-tlai(p)) > 1.e-06_r8 .or. abs(saisum-tsai(p)) > 1.e-06_r8) then
#ifndef _OPENACC
             write (iulog,*) 'multi-layer canopy error 02 in SurfaceAlbedo: ',nrad(p),ncan(p)
             write (iulog,*) tlai(p),elai(p),blai(p),laisum,tsai(p),esai(p),bsai(p),saisum
             call endrun(decomp_index=p, elmlevel=namep, msg=errmsg(__FILE__, __LINE__))
#endif
          end if
       end if

    end do

    ! Zero fluxes for active canopy layers

    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)
       do iv = 1, nrad(p)
          fabd_sun_z(p,iv) = 0._r8
          fabd_sha_z(p,iv) = 0._r8
          fabi_sun_z(p,iv) = 0._r8
          fabi_sha_z(p,iv) = 0._r8
          fsun_z(p,iv) = 0._r8
       end do
    end do

    ! Default leaf to canopy scaling coefficients, used when coszen <= 0.
    ! This is the leaf nitrogen profile integrated over the full canopy.
    ! Integrate exp(-kn*x) over x=0 to x=elai and assign to shaded canopy,
    ! because sunlit fraction is 0. Canopy scaling coefficients are set in
    ! TwoStream for coszen > 0. So kn must be set here and in TwoStream.

    extkn = 0.30_r8
    do fp = 1,num_nourbanp
       p = filter_nourbanp(fp)
       if (nlevcan == 1) then
          vcmaxcintsun(p) = 0._r8
          vcmaxcintsha(p) = (1._r8 - exp(-extkn*elai(p))) / extkn
          if (elai(p) > 0._r8) then
             vcmaxcintsha(p) = vcmaxcintsha(p) / elai(p)
          else
             vcmaxcintsha(p) = 0._r8
          end if
       else if (nlevcan > 1) then
          vcmaxcintsun(p) = 0._r8
          vcmaxcintsha(p) = 0._r8
       end if
    end do

    ! Calculate surface albedos and fluxes
    ! Only perform on vegetated patches where coszen > 0

    ! Calculate surface albedos and fluxes
    ! Only perform on vegetated pfts where coszen > 0
    if(use_fates)then
#ifndef _OPENACC

       call alm_fates%wrap_canopy_radiation(bounds, &
            num_vegsol, filter_vegsol, &
            coszen_patch(bounds%begp:bounds%endp), &
            surfalb_vars)
#endif
    else
    
      if (use_top_solar_rad) then
        call Albedo_TOP_Adjustment(bounds, num_vegsol, filter_vegsol, nextsw_cday, &
                                  coszen_patch(bounds%begp:bounds%endp), declinp1, surfalb_vars, .true.)
      endif

      call TwoStream (bounds, filter_vegsol, num_vegsol, &
               coszen_patch(bounds%begp:bounds%endp), &
               rho(bounds%begp:bounds%endp, :), &
               tau(bounds%begp:bounds%endp, :), &
               canopystate_vars, surfalb_vars, &
               nextsw_cday, declinp1)

    endif

       ! Determine values for non-vegetated patches where coszen > 0

    do ib = 1,numrad
       do fp = 1,num_novegsol
          p = filter_novegsol(fp)
          c = veg_pp%column(p)
          fabd(p,ib) = 0._r8
          fabd_sun(p,ib) = 0._r8
          fabd_sha(p,ib) = 0._r8
          fabi(p,ib) = 0._r8
          fabi_sun(p,ib) = 0._r8
          fabi_sha(p,ib) = 0._r8
          ftdd(p,ib) = 1._r8
          ftid(p,ib) = 0._r8
          ftii(p,ib) = 1._r8
          albd(p,ib) = albgrd(c,ib)
          albi(p,ib) = albgri(c,ib)
       end do
    end do

    if (use_top_solar_rad) then
       call Albedo_TOP_Adjustment(bounds, num_novegsol, filter_novegsol, nextsw_cday, &
                                  coszen_patch(bounds%begp:bounds%endp), declinp1, surfalb_vars, .false.)
    endif

     end associate

   end subroutine SurfaceAlbedo

   !-----------------------------------------------------------------------
   subroutine SoilAlbedo (bounds, &
        num_nourbanc, filter_nourbanc, &
        coszen, albsnd, albsni, &
        lakestate_vars, surfalb_vars)
     !
     ! !DESCRIPTION:
     ! Determine ground surface albedo, accounting for snow
     !
     ! !USES:
      !$acc routine seq
    use elm_varpar, only : numrad
    use elm_varcon      , only : tfrz
    use landunit_varcon, only : istice, istice_mec, istdlak
    use LakeCon         , only : lakepuddling
    use domainMod       , only : ldomain
    !
    ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds
     integer , intent(in) :: num_nourbanc               ! number of columns in non-urban points in column filter
     integer , intent(in) :: filter_nourbanc(:)          ! column filter for non-urban points
     real(r8), intent(in) :: coszen( bounds%begc: )      ! cos solar zenith angle next time step [col]
     real(r8), intent(in) :: albsnd( bounds%begc: , 1: ) ! snow albedo (direct) [col, numrad]
     real(r8), intent(in) :: albsni( bounds%begc: , 1: ) ! snow albedo (diffuse) [col, numrad]
     type(lakestate_type)   , intent(in)    :: lakestate_vars
     type(surfalb_type)     , intent(inout) :: surfalb_vars
     !
     ! !LOCAL VARIABLES:
     !
    integer, parameter :: nband =numrad ! number of solar radiation waveband classes
    integer  :: fc            ! non-urban filter column index
    integer  :: c,l,t,g,gi,pi,mi           ! indices
    integer  :: ib            ! waveband number (1=vis, 2=nir)
    real(r8) :: inc           ! soil water correction factor for soil albedo
    integer  :: soilcol       ! soilcolor
    real(r8) :: sicefr        ! Lake surface ice fraction (based on D. Mironov 2010)
    !-----------------------------------------------------------------------

     ! Enforce expected array sizes

   associate(&
          snl          => col_pp%snl                         , & ! Input:  [integer  (:)   ]  number of snow layers

          t_grnd       => col_es%t_grnd     , & ! Input:  [real(r8) (:)   ]  ground temperature (Kelvin)

          h2osoi_vol   => col_ws%h2osoi_vol  , & ! Input:  [real(r8) (:,:) ]  volumetric soil water [m3/m3]

          lake_icefrac => lakestate_vars%lake_icefrac_col , & ! Input:  [real(r8) (:,:) ]  mass fraction of lake layer that is frozen

          albgrd       => surfalb_vars%albgrd_col         , & ! Output: [real(r8) (:,:) ]  ground albedo (direct)
          albgri       => surfalb_vars%albgri_col         , & ! Output: [real(r8) (:,:) ]  ground albedo (diffuse)
          albsod       => surfalb_vars%albsod_col         , & ! Output: [real(r8) (:,:) ]  soil albedo (direct)
          albsoi       => surfalb_vars%albsoi_col           & ! Output: [real(r8) (:,:) ]  soil albedo (diffuse)
   )

    ! Compute soil albedos

    do ib = 1, nband
       do fc = 1,num_nourbanc
          c = filter_nourbanc(fc)
          g = col_pp%gridcell(c)
          t = col_pp%topounit(c)
          gi = grc_pp%gindex(g)
          mi = ldomain%mask(g)
          pi = ldomain%pftm(g)
          if (coszen(c) > 0._r8) then
             l = col_pp%landunit(c)
             if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop .and. top_pp%active(t))  then ! soil
                inc    = max(0.11_r8-0.40_r8*h2osoi_vol(c,1), 0._r8)
                soilcol = isoicol(c)
                ! changed from local variable to elm_type:
                !albsod = min(albsat(soilcol,ib)+inc, albdry(soilcol,ib))
                !albsoi = albsod
                albsod(c,ib) = min(albsat(soilcol,ib)+inc, albdry(soilcol,ib))
                albsoi(c,ib) = albsod(c,ib)
             else if (lun_pp%itype(l) == istice .or. lun_pp%itype(l) == istice_mec .and. top_pp%active(t))  then  ! land ice
                ! changed from local variable to clm_type:
                !albsod = albice(ib)
                !albsoi = albsod
                albsod(c,ib) = albice(ib)
                albsoi(c,ib) = albsod(c,ib)
             ! unfrozen lake, wetland
             else if (t_grnd(c) > tfrz .or. (lakepuddling .and. lun_pp%itype(l) == istdlak .and. t_grnd(c) == tfrz .and. &
                      lake_icefrac(c,1) < 1._r8 .and. lake_icefrac(c,2) > 0._r8)  .and. top_pp%active(t)) then

                albsod(c,ib) = 0.05_r8/(max(0.001_r8,coszen(c)) + 0.15_r8)
                ! This expression is apparently from BATS according to Yongjiu Dai.

                ! The diffuse albedo should be an average over the whole sky of an angular-dependent direct expression.
                ! The expression above may have been derived to encompass both (e.g. Henderson-Sellers 1986),
                ! but I'll assume it applies more appropriately to the direct form for now.

                ! ZMS: Attn EK, currently restoring this for wetlands even though it is wrong in order to try to get
                ! bfb baseline comparison when no lakes are present. I'm assuming wetlands will be phased out anyway.
                if (lun_pp%itype(l) == istdlak) then
                   albsoi(c,ib) = 0.10_r8
                else
                   albsoi(c,ib) = albsod(c,ib)
                end if

             else                                     ! frozen lake, wetland
                ! Introduce crude surface frozen fraction according to D. Mironov (2010)
                ! Attn EK: This formulation is probably just as good for "wetlands" if they are not phased out.
                ! Tenatively I'm restricting this to lakes because I haven't tested it for wetlands. But if anything
                ! the albedo should be lower when melting over frozen ground than a solid frozen lake.
                !
                if (lun_pp%itype(l) == istdlak .and. .not. lakepuddling .and. snl(c) == 0  .and. top_pp%active(t)) then
                    ! Need to reference snow layers here because t_grnd could be over snow or ice
                                      ! but we really want the ice surface temperature with no snow
                   sicefr = 1._r8 - exp(-calb * (tfrz - t_grnd(c))/tfrz)
                   albsod(c,ib) = sicefr*alblak(ib) + (1._r8-sicefr)*max(alblakwi(ib), &
                                  0.05_r8/(max(0.001_r8,coszen(c)) + 0.15_r8))
                   albsoi(c,ib) = sicefr*alblak(ib) + (1._r8-sicefr)*max(alblakwi(ib), 0.10_r8)
                   ! Make sure this is no less than the open water albedo above.
                   ! Setting lake_melt_icealb(:) = alblak(:) in namelist reverts the melting albedo to the cold
                   ! snow-free value.
                else
                   albsod(c,ib) = alblak(ib)
                   albsoi(c,ib) = albsod(c,ib)
                end if
             end if

             ! Weighting is done in SurfaceAlbedo, after the call to SNICAR_RT/SNICAR_AD_RT
             ! This had to be done, because SoilAlbedo is called before SNICAR_RT/SNICAR_AD_RT, so at
             ! this point, snow albedo is not yet known.
          end if
       end do
    end do

    end associate
   end subroutine SoilAlbedo

   !-----------------------------------------------------------------------
   subroutine TwoStream(bounds, &
        filter_vegsol, num_vegsol, &
        coszen, rho, tau, &
        canopystate_vars, surfalb_vars, nextsw_cday, decl)
     !
     ! !DESCRIPTION:
     ! Two-stream fluxes for canopy radiative transfer
     ! Use two-stream approximation of Dickinson (1983) Adv Geophysics
     ! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
     ! to calculate fluxes absorbed by vegetation, reflected by vegetation,
     ! and transmitted through vegetation for unit incoming direct or diffuse
     ! flux given an underlying surface with known albedo.
     ! Calculate sunlit and shaded fluxes as described by
     ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
     ! a multi-layer canopy to calculate APAR profile
     !
     ! !USES:
      !$acc routine seq
     use elm_varpar, only : numrad, nlevcan
     use elm_varcon, only : omegas, tfrz, betads, betais
     use elm_varctl, only : iulog, use_top_solar_rad
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds
     integer                , intent(in)    :: filter_vegsol (:)        ! filter for vegetated patches with coszen>0
     integer                , intent(in)    :: num_vegsol               ! number of vegetated patches where coszen>0
     real(r8)               , intent(in)    :: coszen( bounds%begp: )   ! cosine solar zenith angle for next time step [pft]
     real(r8)               , intent(in)    :: rho( bounds%begp: , 1: ) ! leaf/stem refl weighted by fraction LAI and SAI [pft, numrad]
     real(r8)               , intent(in)    :: tau( bounds%begp: , 1: ) ! leaf/stem tran weighted by fraction LAI and SAI [pft, numrad]
     type(canopystate_type) , intent(in)    :: canopystate_vars
     type(surfalb_type)     , intent(inout) :: surfalb_vars
     real(r8)               , intent(in)    :: nextsw_cday              ! calendar day at Greenwich (1.00, ..., days/year)
     real(r8)               , intent(in)    :: decl                     ! declination angle (radians) for next time step
     
     !
     ! !LOCAL VARIABLES:
     integer  :: fp,p,c,iv        ! array indices
     integer  :: ib               ! waveband number
     real(r8) :: cosz             ! 0.001 <= coszen <= 1.000
     real(r8) :: asu              ! single scattering albedo
     real(r8) :: chil(bounds%begp:bounds%endp)    ! -0.4 <= xl <= 0.6
     real(r8) :: gdir(bounds%begp:bounds%endp)    ! leaf projection in solar direction (0 to 1)
     real(r8) :: twostext(bounds%begp:bounds%endp)! optical depth of direct beam per unit leaf area
     real(r8) :: avmu(bounds%begp:bounds%endp)    ! average diffuse optical depth
     real(r8) :: omega(bounds%begp:bounds%endp,numrad)   ! fraction of intercepted radiation that is scattered (0 to 1)
     real(r8) :: omegal           ! omega for leaves
     real(r8) :: betai            ! upscatter parameter for diffuse radiation
     real(r8) :: betail           ! betai for leaves
     real(r8) :: betad            ! upscatter parameter for direct beam radiation
     real(r8) :: betadl           ! betad for leaves
     real(r8) :: tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9 ! temporary
     real(r8) :: p1,p2,p3,p4,s1,s2,u1,u2,u3                        ! temporary
     real(r8) :: b,c1,d,d1,d2,f,h,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10   ! temporary
     real(r8) :: phi1,phi2,sigma                                   ! temporary
     real(r8) :: temp1                                             ! temporary
     real(r8) :: temp0    (bounds%begp:bounds%endp)                ! temporary
     real(r8) :: temp2(bounds%begp:bounds%endp)                    ! temporary
     real(r8) :: t1                                                ! temporary
     real(r8) :: albd_adjust,albi_adjust,fd_prime,fi_prime         ! temporary
     real(r8) :: a1,a2                                             ! parameter for sunlit/shaded leaf radiation absorption
     real(r8) :: v,dv,u,du                                         ! temporary for flux derivatives
     real(r8) :: dh2,dh3,dh5,dh6,dh7,dh8,dh9,dh10                  ! temporary for flux derivatives
     real(r8) :: da1,da2                                           ! temporary for flux derivatives
     real(r8) :: d_ftid,d_ftii                                     ! ftid, ftii derivative with respect to lai+sai
     real(r8) :: d_fabd,d_fabi                                     ! fabd, fabi derivative with respect to lai+sai
     real(r8) :: d_fabd_sun,d_fabd_sha                             ! fabd_sun, fabd_sha derivative with respect to lai+sai
     real(r8) :: d_fabi_sun,d_fabi_sha                             ! fabi_sun, fabi_sha derivative with respect to lai+sai
     real(r8) :: laisum                                            ! cumulative lai+sai for canopy layer (at middle of layer)
     real(r8) :: extkb                                             ! direct beam extinction coefficient
     real(r8) :: extkn                                             ! nitrogen allocation coefficient
     !-----------------------------------------------------------------------

     ! Enforce expected array sizes

     associate(&
          xl            =>    veg_vp%xl                           , & ! Input:  [real(r8) (:)   ]  ecophys const - leaf/stem orientation index

          t_veg         =>    veg_es%t_veg                        , & ! Input:  [real(r8) (:)   ]  vegetation temperature (Kelvin)

          fwet          =>    veg_ws%fwet                         , & ! Input:  [real(r8) (:)   ]  fraction of canopy that is wet (0 to 1)

          elai          =>    canopystate_vars%elai_patch         , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow
          esai          =>    canopystate_vars%esai_patch         , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow

          tlai_z        =>    surfalb_vars%tlai_z_patch           , & ! Input:  [real(r8) (:,:) ]  tlai increment for canopy layer
          tsai_z        =>    surfalb_vars%tsai_z_patch           , & ! Input:  [real(r8) (:,:) ]  tsai increment for canopy layer
          nrad          =>    surfalb_vars%nrad_patch             , & ! Input:  [integer  (:)   ]  number of canopy layers, above snow for radiative transfer
          albgrd        =>    surfalb_vars%albgrd_col             , & ! Input:  [real(r8) (:,:) ]  ground albedo (direct) (column-level)
          albgri        =>    surfalb_vars%albgri_col             , & ! Input:  [real(r8) (:,:) ]  ground albedo (diffuse)(column-level)
          
          fd_top_adjust =>    surfalb_vars%fd_top_adjust          , & ! Input: TOP adjusted factor for direct radiation 
	  fi_top_adjust =>    surfalb_vars%fi_top_adjust          , & ! Input: TOP adjusted factor for diffuse radiation
          f_dir         =>    surfalb_vars%f_dir                  , & ! Input: TOP adjusted factor for direct radiation 
	  f_rdir        =>    surfalb_vars%f_rdir                 , & ! Input: TOP adjusted factor for reflected-direct radiation 
          f_dif         =>    surfalb_vars%f_dif                  , & ! Input: TOP adjusted factor for diffuse radiation 
	  f_rdif        =>    surfalb_vars%f_rdif                 , & ! Input: TOP adjusted factor for reflected-diffuse radiation 
          
          fsun_z        =>    surfalb_vars%fsun_z_patch           , & ! Output: [real(r8) (:,:) ]  sunlit fraction of canopy layer
          vcmaxcintsun  =>    surfalb_vars%vcmaxcintsun_patch     , & ! Output: [real(r8) (:)   ]  leaf to canopy scaling coefficient, sunlit leaf vcmax
          vcmaxcintsha  =>    surfalb_vars%vcmaxcintsha_patch     , & ! Output: [real(r8) (:)   ]  leaf to canopy scaling coefficient, shaded leaf vcmax
          fabd_sun_z    =>    surfalb_vars%fabd_sun_z_patch       , & ! Output: [real(r8) (:,:) ]  absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabd_sha_z    =>    surfalb_vars%fabd_sha_z_patch       , & ! Output: [real(r8) (:,:) ]  absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabi_sun_z    =>    surfalb_vars%fabi_sun_z_patch       , & ! Output: [real(r8) (:,:) ]  absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
          fabi_sha_z    =>    surfalb_vars%fabi_sha_z_patch       , & ! Output: [real(r8) (:,:) ]  absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
          albd          =>    surfalb_vars%albd_patch             , & ! Output: [real(r8) (:,:) ]  surface albedo (direct)
          albi          =>    surfalb_vars%albi_patch             , & ! Output: [real(r8) (:,:) ]  surface albedo (diffuse)
          fabd          =>    surfalb_vars%fabd_patch             , & ! Output: [real(r8) (:,:) ]  flux absorbed by canopy per unit direct flux
          fabd_sun      =>    surfalb_vars%fabd_sun_patch         , & ! Output: [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit direct flux
          fabd_sha      =>    surfalb_vars%fabd_sha_patch         , & ! Output: [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit direct flux
          fabi          =>    surfalb_vars%fabi_patch             , & ! Output: [real(r8) (:,:) ]  flux absorbed by canopy per unit diffuse flux
          fabi_sun      =>    surfalb_vars%fabi_sun_patch         , & ! Output: [real(r8) (:,:) ]  flux absorbed by sunlit canopy per unit diffuse flux
          fabi_sha      =>    surfalb_vars%fabi_sha_patch         , & ! Output: [real(r8) (:,:) ]  flux absorbed by shaded canopy per unit diffuse flux
          ftdd          =>    surfalb_vars%ftdd_patch             , & ! Output: [real(r8) (:,:) ]  down direct flux below canopy per unit direct flx
          ftid          =>    surfalb_vars%ftid_patch             , & ! Output: [real(r8) (:,:) ]  down diffuse flux below canopy per unit direct flx
          ftii          =>    surfalb_vars%ftii_patch               & ! Output: [real(r8) (:,:) ]  down diffuse flux below canopy per unit diffuse flx
          )

    ! Calculate two-stream parameters that are independent of waveband:
    ! chil, gdir, twostext, avmu, and temp0 and temp2 (used for asu)

    do fp = 1,num_vegsol
       p = filter_vegsol(fp)

       ! note that the following limit only acts on cosz values > 0 and less than
       ! 0.001, not on values cosz = 0, since these zero have already been filtered
       ! out in filter_vegsol
       cosz = max(0.001_r8, coszen(p))

       chil(p) = min( max(xl(veg_pp%itype(p)), -0.4_r8), 0.6_r8 )
       if (abs(chil(p)) <= 0.01_r8) chil(p) = 0.01_r8
       phi1 = 0.5_r8 - 0.633_r8*chil(p) - 0.330_r8*chil(p)*chil(p)
       phi2 = 0.877_r8 * (1._r8-2._r8*phi1)
       gdir(p) = phi1 + phi2*cosz
       twostext(p) = gdir(p)/cosz
       avmu(p) = ( 1._r8 - phi1/phi2 * log((phi1+phi2)/phi1) ) / phi2
       temp0(p) = gdir(p) + phi2*cosz
       temp1 = phi1*cosz
       temp2(p) = ( 1._r8 - temp1/temp0(p) * log((temp1+temp0(p))/temp1) )
    end do

   ! Loop over all wavebands to calculate for the full canopy the scattered fluxes
   ! reflected upward and transmitted downward by the canopy and the flux absorbed by the
   ! canopy for a unit incoming direct beam and diffuse flux at the top of the canopy given
   ! an underlying surface of known albedo.
   !
   ! Output:
   ! ------------------
   ! Direct beam fluxes
   ! ------------------
   ! albd       - Upward scattered flux above canopy (per unit direct beam flux)
   ! ftid       - Downward scattered flux below canopy (per unit direct beam flux)
   ! ftdd       - Transmitted direct beam flux below canopy (per unit direct beam flux)
   ! fabd       - Flux absorbed by canopy (per unit direct beam flux)
   ! fabd_sun   - Sunlit portion of fabd
   ! fabd_sha   - Shaded portion of fabd
   ! fabd_sun_z - absorbed sunlit leaf direct PAR (per unit sunlit lai+sai) for each canopy layer
   ! fabd_sha_z - absorbed shaded leaf direct PAR (per unit shaded lai+sai) for each canopy layer
   ! ------------------
   ! Diffuse fluxes
   ! ------------------
   ! albi       - Upward scattered flux above canopy (per unit diffuse flux)
   ! ftii       - Downward scattered flux below canopy (per unit diffuse flux)
   ! fabi       - Flux absorbed by canopy (per unit diffuse flux)
   ! fabi_sun   - Sunlit portion of fabi
   ! fabi_sha   - Shaded portion of fabi
   ! fabi_sun_z - absorbed sunlit leaf diffuse PAR (per unit sunlit lai+sai) for each canopy layer
   ! fabi_sha_z - absorbed shaded leaf diffuse PAR (per unit shaded lai+sai) for each canopy layer

    do ib = 1, numrad
       do fp = 1,num_vegsol
          p = filter_vegsol(fp)
          c = veg_pp%column(p)

          ! Calculate two-stream parameters omega, betad, and betai.
          ! Omega, betad, betai are adjusted for snow. Values for omega*betad
          ! and omega*betai are calculated and then divided by the new omega
          ! because the product omega*betai, omega*betad is used in solution.
          ! Also, the transmittances and reflectances (tau, rho) are linear
          ! weights of leaf and stem values.
          cosz = max(0.001_r8, coszen(p))
          
          omegal = rho(p,ib) + tau(p,ib)
          asu = 0.5_r8*omegal*gdir(p)/temp0(p) *temp2(p)
          betadl = (1._r8+avmu(p)*twostext(p))/(omegal*avmu(p)*twostext(p))*asu
          betail = 0.5_r8 * ((rho(p,ib)+tau(p,ib)) + (rho(p,ib)-tau(p,ib)) &
                 * ((1._r8+chil(p))/2._r8)**2) / omegal

          ! Adjust omega, betad, and betai for intercepted snow

          ! NOTE: If improvements to the canopy snow module emerge,
          ! Please notify the FATES team (RGK,GL,RAF,CDK,etc) to help us adapt
          ! changes into the FATES albedo and canopy radiation calculations.
          ! At present (01-2022), we are using a similar calculation here,
          ! where fwet = fsnow if t_veg <= tfrz, and 0 otherwise.

          if (t_veg(p) > tfrz) then                             !no snow
             tmp0 = omegal
             tmp1 = betadl
             tmp2 = betail
          else
             tmp0 =   (1._r8-fwet(p))*omegal        + fwet(p)*omegas(ib)
             tmp1 = ( (1._r8-fwet(p))*omegal*betadl + fwet(p)*omegas(ib)*betads ) / tmp0
             tmp2 = ( (1._r8-fwet(p))*omegal*betail + fwet(p)*omegas(ib)*betais ) / tmp0
          end if
          omega(p,ib) = tmp0
          betad = tmp1
          betai = tmp2

          ! Common terms

          b = 1._r8 - omega(p,ib) + omega(p,ib)*betai
          c1 = omega(p,ib)*betai
          tmp0 = avmu(p)*twostext(p)
          d = tmp0 * omega(p,ib)*betad
          f = tmp0 * omega(p,ib)*(1._r8-betad)
          tmp1 = b*b - c1*c1
          h = sqrt(tmp1) / avmu(p)
          sigma = tmp0*tmp0 - tmp1
          p1 = b + avmu(p)*h
          p2 = b - avmu(p)*h
          p3 = b + tmp0
          p4 = b - tmp0

          ! Absorbed, reflected, transmitted fluxes per unit incoming radiation
          ! for full canopy

          t1 = min(h*(elai(p)+esai(p)), 40._r8)
          s1 = exp(-t1)
          t1 = min(twostext(p)*(elai(p)+esai(p)), 40._r8)
          s2 = exp(-t1)

          ! Direct beam

          u1 = b - c1/albgrd(c,ib)
          u2 = b - c1*albgrd(c,ib)
          u3 = f + c1*albgrd(c,ib)
          tmp2 = u1 - avmu(p)*h
          tmp3 = u1 + avmu(p)*h
          d1 = p1*tmp2/s1 - p2*tmp3*s1
          tmp4 = u2 + avmu(p)*h
          tmp5 = u2 - avmu(p)*h
          d2 = tmp4/s1 - tmp5*s1
          h1 = -d*p4 - c1*f
          tmp6 = d - h1*p3/sigma
          tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
          h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
          h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
          h4 = -f*p3 - c1*d
          tmp8 = h4/sigma
          tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
          h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
          h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2

          albd(p,ib) = h1/sigma + h2 + h3
          ftid(p,ib) = h4*s2/sigma + h5*s1 + h6/s1
          ftdd(p,ib) = s2
          fabd(p,ib) = 1._r8 - albd(p,ib) - (1._r8-albgrd(c,ib))*ftdd(p,ib) - (1._r8-albgri(c,ib))*ftid(p,ib)
  
          if (use_top_solar_rad) then            
            f_rdir(p, ib) = f_rdir(p, ib) * (albd(p,ib) / 0.1_r8)     
            fd_prime = 1._r8 + f_dir(p, ib) + f_rdir(p, ib)
            albd_adjust = fd_prime * albd(p,ib) - (fd_prime-1._r8)

            if (albd_adjust <= 0._r8) then 
               albd_adjust = 0._r8
               fd_prime = 1._r8 / (1._r8 - albd(p,ib))
            endif

            albd(p,ib) = albd_adjust
            fabd(p,ib) = fabd(p,ib) * fd_prime
            ftdd(p,ib) = ftdd(p,ib) * fd_prime
            ftid(p,ib) = ftid(p,ib) * fd_prime
            fd_top_adjust(p,ib) = fd_prime
          endif

          a1 = h1 / sigma * (1._r8 - s2*s2) / (2._r8 * twostext(p)) &
             + h2         * (1._r8 - s2*s1) / (twostext(p) + h) &
             + h3         * (1._r8 - s2/s1) / (twostext(p) - h)

          a2 = h4 / sigma * (1._r8 - s2*s2) / (2._r8 * twostext(p)) &
             + h5         * (1._r8 - s2*s1) / (twostext(p) + h) &
             + h6         * (1._r8 - s2/s1) / (twostext(p) - h)

          if (use_top_solar_rad) then
             fabd_sun(p,ib) = (1._r8 - omega(p,ib)) * ( 1._r8 - s2 + 1._r8 / avmu(p) * (a1 + a2) ) * fd_top_adjust(p,ib)
          else
             fabd_sun(p,ib) = (1._r8 - omega(p,ib)) * ( 1._r8 - s2 + 1._r8 / avmu(p) * (a1 + a2) )
          endif
          fabd_sha(p,ib) = fabd(p,ib) - fabd_sun(p,ib)
          ! Diffuse

          u1 = b - c1/albgri(c,ib)
          u2 = b - c1*albgri(c,ib)
          tmp2 = u1 - avmu(p)*h
          tmp3 = u1 + avmu(p)*h
          d1 = p1*tmp2/s1 - p2*tmp3*s1
          tmp4 = u2 + avmu(p)*h
          tmp5 = u2 - avmu(p)*h
          d2 = tmp4/s1 - tmp5*s1
          h7 = (c1*tmp2) / (d1*s1)
          h8 = (-c1*tmp3*s1) / d1
          h9 = tmp4 / (d2*s1)
          h10 = (-tmp5*s1) / d2

          albi(p,ib) = h7 + h8
          ftii(p,ib) = h9*s1 + h10/s1
          fabi(p,ib) = 1._r8 - albi(p,ib) - (1._r8-albgri(c,ib))*ftii(p,ib)

          if (use_top_solar_rad) then
            f_rdif(p,ib) = f_rdif(p,ib) * (albi(p,ib) / 0.1_r8)
            fi_prime = 1._r8 + f_dif(p,ib) + f_rdif(p,ib)
            albi_adjust = fi_prime * albi(p,ib) - (fi_prime-1._r8)  

            if (albi_adjust <= 0._r8) then
               albi_adjust = 0._r8
               fi_prime = 1._r8 / (1._r8 - albi(p,ib))
            endif

            albi(p,ib) = albi_adjust
            fabi(p,ib) = fabi(p,ib) * fi_prime
            ftii(p,ib) = ftii(p,ib) * fi_prime
            fi_top_adjust(p,ib) = fi_prime
          endif
          
          a1 = h7 * (1._r8 - s2*s1) / (twostext(p) + h) +  h8 * (1._r8 - s2/s1) / (twostext(p) - h)
          a2 = h9 * (1._r8 - s2*s1) / (twostext(p) + h) + h10 * (1._r8 - s2/s1) / (twostext(p) - h)

          if (use_top_solar_rad) then
            fabi_sun(p,ib) = (1._r8 - omega(p,ib)) / avmu(p) * (a1 + a2) * fi_top_adjust(p,ib)
          else
            fabi_sun(p,ib) = (1._r8 - omega(p,ib)) / avmu(p) * (a1 + a2)
          endif
          fabi_sha(p,ib) = fabi(p,ib) - fabi_sun(p,ib)

          ! Repeat two-stream calculations for each canopy layer to calculate derivatives.
          ! tlai_z and tsai_z are the leaf+stem area increment for a layer. Derivatives are
          ! calculated at the center of the layer. Derivatives are needed only for the
          ! visible waveband to calculate absorbed PAR (per unit lai+sai) for each canopy layer.
          ! Derivatives are calculated first per unit lai+sai and then normalized for sunlit
          ! or shaded fraction of canopy layer.

          ! Sun/shade big leaf code uses only one layer, with canopy integrated values from above
          ! and also canopy-integrated scaling coefficients

          if (ib == 1) then
             if (nlevcan == 1) then

                ! sunlit fraction of canopy
                fsun_z(p,1) = (1._r8 - s2) / t1

                ! absorbed PAR (per unit sun/shade lai+sai)
                laisum = elai(p)+esai(p)
                fabd_sun_z(p,1) = fabd_sun(p,ib) / (fsun_z(p,1)*laisum)
                fabi_sun_z(p,1) = fabi_sun(p,ib) / (fsun_z(p,1)*laisum)
                fabd_sha_z(p,1) = fabd_sha(p,ib) / ((1._r8 - fsun_z(p,1))*laisum)
                fabi_sha_z(p,1) = fabi_sha(p,ib) / ((1._r8 - fsun_z(p,1))*laisum)

                ! leaf to canopy scaling coefficients
                extkn = 0.30_r8
                extkb = twostext(p)
                vcmaxcintsun(p) = (1._r8 - exp(-(extkn+extkb)*elai(p))) / (extkn + extkb)
                vcmaxcintsha(p) = (1._r8 - exp(-extkn*elai(p))) / extkn - vcmaxcintsun(p)
                if (elai(p)  >  0._r8) then
                  vcmaxcintsun(p) = vcmaxcintsun(p) / (fsun_z(p,1)*elai(p))
                  vcmaxcintsha(p) = vcmaxcintsha(p) / ((1._r8 - fsun_z(p,1))*elai(p))
                else
                  vcmaxcintsun(p) = 0._r8
                  vcmaxcintsha(p) = 0._r8
                end if

             else if (nlevcan > 1) then
                do iv = 1, nrad(p)

                ! Cumulative lai+sai at center of layer

                if (iv == 1) then
                   laisum = 0.5_r8 * (tlai_z(p,iv)+tsai_z(p,iv))
                else
                   laisum = laisum + 0.5_r8 * ((tlai_z(p,iv-1)+tsai_z(p,iv-1))+(tlai_z(p,iv)+tsai_z(p,iv)))
                end if

                ! Coefficients s1 and s2 depend on cumulative lai+sai. s2 is the sunlit fraction

                t1 = min(h*laisum, 40._r8)
                s1 = exp(-t1)
                t1 = min(twostext(p)*laisum, 40._r8)
                s2 = exp(-t1)
                fsun_z(p,iv) = s2

                ! ===============
                ! Direct beam
                ! ===============

                ! Coefficients h1-h6 and a1,a2 depend of cumulative lai+sai

                u1 = b - c1/albgrd(c,ib)
                u2 = b - c1*albgrd(c,ib)
                u3 = f + c1*albgrd(c,ib)
                tmp2 = u1 - avmu(p)*h
                tmp3 = u1 + avmu(p)*h
                d1 = p1*tmp2/s1 - p2*tmp3*s1
                tmp4 = u2 + avmu(p)*h
                tmp5 = u2 - avmu(p)*h
                d2 = tmp4/s1 - tmp5*s1
                h1 = -d*p4 - c1*f
                tmp6 = d - h1*p3/sigma
                tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
                h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
                h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
                h4 = -f*p3 - c1*d
                tmp8 = h4/sigma
                tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
                h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
                h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2

                a1 = h1 / sigma * (1._r8 - s2*s2) / (2._r8 * twostext(p)) &
                   + h2         * (1._r8 - s2*s1) / (twostext(p) + h) &
                   + h3         * (1._r8 - s2/s1) / (twostext(p) - h)

                a2 = h4 / sigma * (1._r8 - s2*s2) / (2._r8 * twostext(p)) &
                   + h5         * (1._r8 - s2*s1) / (twostext(p) + h) &
                   + h6         * (1._r8 - s2/s1) / (twostext(p) - h)

                ! Derivatives for h2, h3, h5, h6 and a1, a2

                v = d1
                dv = h * p1 * tmp2 / s1 + h * p2 * tmp3 * s1

                u = tmp6 * tmp2 / s1 - p2 * tmp7
                du = h * tmp6 * tmp2 / s1 + twostext(p) * p2 * tmp7
                dh2 = (v * du - u * dv) / (v * v)

                u = -tmp6 * tmp3 * s1 + p1 * tmp7
                du = h * tmp6 * tmp3 * s1 - twostext(p) * p1 * tmp7
                dh3 = (v * du - u * dv) / (v * v)

                v = d2
                dv = h * tmp4 / s1 + h * tmp5 * s1

                u = -h4/sigma * tmp4 / s1 - tmp9
                du = -h * h4/sigma * tmp4 / s1 + twostext(p) * tmp9
                dh5 = (v * du - u * dv) / (v * v)

                u = h4/sigma * tmp5 * s1 + tmp9
                du = -h * h4/sigma * tmp5 * s1 - twostext(p) * tmp9
                dh6 = (v * du - u * dv) / (v * v)

                da1 = h1/sigma * s2*s2 + h2 * s2*s1 + h3 * s2/s1 &
                    + (1._r8 - s2*s1) / (twostext(p) + h) * dh2 &
                    + (1._r8 - s2/s1) / (twostext(p) - h) * dh3
                da2 = h4/sigma * s2*s2 + h5 * s2*s1 + h6 * s2/s1 &
                    + (1._r8 - s2*s1) / (twostext(p) + h) * dh5 &
                    + (1._r8 - s2/s1) / (twostext(p) - h) * dh6

                ! Flux derivatives

                d_ftid = -twostext(p)*h4/sigma*s2 - h*h5*s1 + h*h6/s1 + dh5*s1 + dh6/s1
                d_fabd = -(dh2+dh3) + (1._r8-albgrd(c,ib))*twostext(p)*s2 - (1._r8-albgri(c,ib))*d_ftid
                if (use_top_solar_rad) then
                   d_fabd_sun = (1._r8 - omega(p,ib)) * (twostext(p)*s2 + 1._r8 / avmu(p) * (da1 + da2)) * fd_top_adjust(p,ib)
                else
                   d_fabd_sun = (1._r8 - omega(p,ib)) * (twostext(p)*s2 + 1._r8 / avmu(p) * (da1 + da2))
                endif
                d_fabd_sha = d_fabd - d_fabd_sun

                fabd_sun_z(p,iv) = max(d_fabd_sun, 0._r8)
                fabd_sha_z(p,iv) = max(d_fabd_sha, 0._r8)

                ! Flux derivatives are APARsun and APARsha per unit (LAI+SAI). Need
                ! to normalize derivatives by sunlit or shaded fraction to get
                ! APARsun per unit (LAI+SAI)sun and APARsha per unit (LAI+SAI)sha

                fabd_sun_z(p,iv) = fabd_sun_z(p,iv) / fsun_z(p,iv)
                fabd_sha_z(p,iv) = fabd_sha_z(p,iv) / (1._r8 - fsun_z(p,iv))

                ! ===============
                ! Diffuse
                ! ===============

                ! Coefficients h7-h10 and a1,a2 depend of cumulative lai+sai

                u1 = b - c1/albgri(c,ib)
                u2 = b - c1*albgri(c,ib)
                tmp2 = u1 - avmu(p)*h
                tmp3 = u1 + avmu(p)*h
                d1 = p1*tmp2/s1 - p2*tmp3*s1
                tmp4 = u2 + avmu(p)*h
                tmp5 = u2 - avmu(p)*h
                d2 = tmp4/s1 - tmp5*s1
                h7 = (c1*tmp2) / (d1*s1)
                h8 = (-c1*tmp3*s1) / d1
                h9 = tmp4 / (d2*s1)
                h10 = (-tmp5*s1) / d2

                a1 = h7 * (1._r8 - s2*s1) / (twostext(p) + h) +  h8 * (1._r8 - s2/s1) / (twostext(p) - h)
                a2 = h9 * (1._r8 - s2*s1) / (twostext(p) + h) + h10 * (1._r8 - s2/s1) / (twostext(p) - h)

                ! Derivatives for h7, h8, h9, h10 and a1, a2

                v = d1
                dv = h * p1 * tmp2 / s1 + h * p2 * tmp3 * s1

                u = c1 * tmp2 / s1
                du = h * c1 * tmp2 / s1
                dh7 = (v * du - u * dv) / (v * v)

                u = -c1 * tmp3 * s1
                du = h * c1 * tmp3 * s1
                dh8 = (v * du - u * dv) / (v * v)

                v = d2
                dv = h * tmp4 / s1 + h * tmp5 * s1

                u = tmp4 / s1
                du = h * tmp4 / s1
                dh9 = (v * du - u * dv) / (v * v)

                u = -tmp5 * s1
                du = h * tmp5 * s1
                dh10 = (v * du - u * dv) / (v * v)

                da1 = h7*s2*s1 +  h8*s2/s1 + (1._r8-s2*s1)/(twostext(p)+h)*dh7 + (1._r8-s2/s1)/(twostext(p)-h)*dh8
                da2 = h9*s2*s1 + h10*s2/s1 + (1._r8-s2*s1)/(twostext(p)+h)*dh9 + (1._r8-s2/s1)/(twostext(p)-h)*dh10

                ! Flux derivatives

                d_ftii = -h * h9 * s1 + h * h10 / s1 + dh9 * s1 + dh10 / s1
                d_fabi = -(dh7+dh8) - (1._r8-albgri(c,ib))*d_ftii
                if (use_top_solar_rad) then
                   d_fabi_sun = (1._r8 - omega(p,ib)) / avmu(p) * (da1 + da2) * fi_top_adjust(p,ib)
                else
                   d_fabi_sun = (1._r8 - omega(p,ib)) / avmu(p) * (da1 + da2)
                endif
                d_fabi_sha = d_fabi - d_fabi_sun

                fabi_sun_z(p,iv) = max(d_fabi_sun, 0._r8)
                fabi_sha_z(p,iv) = max(d_fabi_sha, 0._r8)

                ! Flux derivatives are APARsun and APARsha per unit (LAI+SAI). Need
                ! to normalize derivatives by sunlit or shaded fraction to get
                ! APARsun per unit (LAI+SAI)sun and APARsha per unit (LAI+SAI)sha

                fabi_sun_z(p,iv) = fabi_sun_z(p,iv) / fsun_z(p,iv)
                fabi_sha_z(p,iv) = fabi_sha_z(p,iv) / (1._r8 - fsun_z(p,iv))

                end do   ! end of canopy layer loop
             end if
          end if

       end do   ! end of pft loop
     end do   ! end of radiation band loop

     end associate

    end subroutine TwoStream
    
!----------------------------------------------------------------------
  subroutine Albedo_TOP_Adjustment(bounds, num_pft, filter_pft, &
                                  nextsw_cday, coszen, decl, surfalb_vars, is_veg)
! !DESCRIPTION:
! Adjust surface albedo by accounting for sub-grid topographic effect on surface solar radiation
!

! !USES:
    use shr_orb_mod
    use elm_varctl  , only: iulog

!
! !ARGUMENTS:
    implicit none
    type(bounds_type)  , intent(in)   :: bounds                     ! bounds
    integer            , intent(in)   :: num_pft                    ! number of pfts
    integer            , intent(in)   :: filter_pft(:)              ! bounds%endg-bounds%endg+1 pft filter
    real(r8)           , intent(in)   :: nextsw_cday                ! calendar day at Greenwich (1.00, ..., days/year)
    real(r8)           , intent(in)   :: coszen( bounds%begp: )     ! cos solar zenith angle next time step [gridcell]
    real(r8)           , intent(in)   :: decl                       ! declination angle (radians) for next time step
    logical            , intent(in)   :: is_veg                     ! True: vegetation; False: non-vegetation
    type(surfalb_type) , intent(inout):: surfalb_vars

 
!
! !CALLED FROM:
! subroutine SurfaceAlbedo

! !OTHER LOCAL VARIABLES:
!
    real(r8), parameter :: mpe = 1.e-06_r8                ! prevents overflow for division by zero
    real(r8), parameter :: pi = 3.14159265358979323846_r8 ! pi
    integer  :: fp,fc,g,c,p                               ! indices
    integer  :: ib                                        ! band index
    integer  :: ic                                        ! 0=unit incoming direct; 1=unit incoming diffuse
    integer  :: izen
    real(r8) :: lon_180                                   ! lon starting from -180
    real(r8) :: cosz                                      ! cosine solar zenith angle for next time step
    real(r8) :: sinz                                      ! sine of solar zenith angle
    real(r8) :: azi_angle                                 ! solar azimuth angle
    real(r8) :: next_tod                                  ! time of day for nextsw_cday in second
    real(r8) :: solar_inc                                 ! (solar incident angle) / cos(slope) / cosz
    real(r8) :: fd_prime                                  ! temp adjustment factor for direct flux
    real(r8) :: fi_prime                                  ! temp adjustment factor for diffuse flux
    real(r8) :: albd_adjust                               ! adjusted albedo for direct flux
    real(r8) :: albi_adjust                               ! adjusted albedo for diffuse flux
    real(r8) :: ftemp1, ftemp2, dzen1, dzen2
    real(r8) :: local_timeofday                           ! local time of day (second)
    real(r8) :: coeff_dir(3,0:7)                          ! regression coefficients for f_dir
    real(r8) :: coeff_rdir(3,0:7)                         ! regression coefficients for f_rdir
    real(r8) :: coeff_dif(4,0:7)                          ! regression coefficients for f_dif
    real(r8) :: coeff_rdif(3,0:7)                         ! regression coefficients for f_rdif

    data coeff_dir(:,1) /2.045E+1_r8, 6.792E-1_r8, -2.103E+1_r8/
    data coeff_dir(:,2) /1.993E+0_r8, 9.284E-1_r8, -2.911E+0_r8/
    data coeff_dir(:,3) /5.900E-2_r8, 9.863E-1_r8, -1.045E+0_r8/
    data coeff_dir(:,4) /5.270E-3_r8, 9.942E-1_r8, -9.995E-1_r8/
    data coeff_dir(:,5) /2.977E-3_r8, 9.959E-1_r8, -9.990E-1_r8/
    data coeff_dir(:,6) /2.977E-3_r8, 9.959E-1_r8, -9.990E-1_r8/
    data coeff_dir(:,7) /8.347E-3_r8,       0._r8, -8.393E-3_r8/

    data coeff_rdir(:,1) / 2.351E-1_r8, 1.590E-1_r8, -2.332E-1_r8/
    data coeff_rdir(:,2) / 1.368E-1_r8, 1.642E-1_r8, -1.358E-1_r8/
    data coeff_rdir(:,3) / 1.254E-1_r8, 1.653E-1_r8, -1.247E-1_r8/
    data coeff_rdir(:,4) / 1.274E-1_r8, 1.635E-1_r8, -1.267E-1_r8/
    data coeff_rdir(:,5) / 1.314E-1_r8, 1.623E-1_r8, -1.307E-1_r8/
    data coeff_rdir(:,6) / 1.359E-1_r8, 1.620E-1_r8, -1.352E-1_r8/
    data coeff_rdir(:,7) /-4.463E-6_r8, 1.556E-1_r8,  1.287E-3_r8/

    data coeff_dif(:,1) /3.146E-7_r8, 4.385E+0_r8, 6.723E-3_r8, -4.382E+0_r8/
    data coeff_dif(:,2) /6.001E-7_r8, 4.068E+0_r8, 2.456E-2_r8, -4.085E+0_r8/
    data coeff_dif(:,3) /7.436E-7_r8, 3.911E+0_r8, 5.606E-2_r8, -3.960E+0_r8/
    data coeff_dif(:,4) /7.806E-7_r8, 3.763E+0_r8, 1.049E-1_r8, -3.863E+0_r8/
    data coeff_dif(:,5) /7.581E-7_r8, 3.559E+0_r8, 1.734E-1_r8, -3.727E+0_r8/
    data coeff_dif(:,6) /7.015E-7_r8, 3.298E+0_r8, 2.543E-1_r8, -3.547E+0_r8/
    data coeff_dif(:,7) /6.359E-7_r8, 2.984E+0_r8,       0._r8, -2.984E+0_r8/

    data coeff_rdif(:,1) / 1.493E-1_r8, 1.621E-1_r8, -1.483E-1_r8/
    data coeff_rdif(:,2) / 1.462E-1_r8, 1.654E-1_r8, -1.454E-1_r8/
    data coeff_rdif(:,3) / 1.454E-1_r8, 1.673E-1_r8, -1.446E-1_r8/
    data coeff_rdif(:,4) / 1.465E-1_r8, 1.683E-1_r8, -1.457E-1_r8/
    data coeff_rdif(:,5) / 1.443E-1_r8, 1.682E-1_r8, -1.435E-1_r8/
    data coeff_rdif(:,6) / 1.446E-1_r8, 1.686E-1_r8, -1.439E-1_r8/
    data coeff_rdif(:,7) /-3.427E-6_r8, 1.576E-1_r8,  1.199E-3_r8/

     ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(coszen) == (/bounds%endp/)),         errMsg(__FILE__, __LINE__))

    ! Assign local pointers to derived subtypes components (gridcell-level)

    associate(&
          lat            =>    grc_pp%lat                         , & ! Input:   latitude             
          lon            =>    grc_pp%lon                         , & ! Input:   longitude               
          pgridcell      =>    veg_pp%gridcell                    , & ! Input:   gridcell 
          pcolumn        =>    veg_pp%column                      , & ! Input:   column 
          stdev_elev     =>    grc_pp%stdev_elev                  , & ! Input:   standard deviation of elevation 
          sky_view       =>    grc_pp%sky_view                    , & ! Input:   sky view factor
          terrain_config =>    grc_pp%terrain_config              , & ! Input:   terrain configuration factor
	  sinsl_cosas    =>    grc_pp%sinsl_cosas                 , & ! Input:   sin(slope) * cos(aspect)
          sinsl_sinas    =>    grc_pp%sinsl_sinas                 , & ! Input:   sin(slope) * sin(aspect)
          albd           =>    surfalb_vars%albd_patch            , & ! Output:  surface albedo (direct)               
          albi           =>    surfalb_vars%albi_patch            , & ! Output:  surface albedo (diffuse)              
          fabd           =>    surfalb_vars%fabd_patch            , & ! Output:  flux absorbed by canopy per unit direct flux
          fabi           =>    surfalb_vars%fabi_patch            , & ! Output:  flux absorbed by canopy per unit diffuse flux
          ftdd           =>    surfalb_vars%ftdd_patch            , & ! Output:  down direct flux below canopy per unit direct flux
          ftid           =>    surfalb_vars%ftid_patch            , & ! Output:  down diffuse flux below canopy per unit direct flux
          ftii           =>    surfalb_vars%ftii_patch            , & ! Output:  down diffuse flux below canopy per unit diffuse flux
          fd_top_adjust  =>    surfalb_vars%fd_top_adjust         , & ! Output:  TOP adjusted factor for direct radiation
          fi_top_adjust  =>    surfalb_vars%fi_top_adjust         , & ! Output:  TOP adjusted factor for diffuse radiation
          f_dir          =>    surfalb_vars%f_dir                 , & ! Output:  TOP adjusted factor for direct radiation 
	  f_rdir         =>    surfalb_vars%f_rdir                , & ! Output:  TOP adjusted factor for reflected_direct radiation 
          f_dif          =>    surfalb_vars%f_dif                 , & ! Output:  TOP adjusted factor for diffuse radiation 
	  f_rdif         =>    surfalb_vars%f_rdif                  & ! Output:  TOP adjusted factor for reflected_diffuse radiation 
          )



     coeff_dir(:,0) = coeff_dir(:,1)
     coeff_rdir(:,0) = coeff_rdir(:,1)
     coeff_dif(:,0) = coeff_dif(:,1)
     coeff_rdif(:,0) = coeff_rdif(:,1)

     next_tod = 86400._r8 * (nextsw_cday - int(nextsw_cday))

     do fp = 1,num_pft
        p = filter_pft(fp)
        g = pgridcell(p)
        
        cosz = coszen(p)
        fd_top_adjust(p,1:numrad) = 1._r8
        fi_top_adjust(p,1:numrad) = 1._r8
	   
       ! make sure the lon is between 0-180 
        lon_180 = lon(g)
         if (lon_180 > pi) lon_180 = lon_180-2._r8*pi    
    
         if (cosz > 0._r8 .and. abs(lat(g)) < 1.047_r8 .and. stdev_elev(g) > 0._r8) then
            local_timeofday = next_tod + lon_180 / pi * 180._r8 * 240._r8

            if (local_timeofday >= 86400._r8) then
               local_timeofday = local_timeofday - 86400._r8
            endif
    
            if (local_timeofday < 0._r8) then
               local_timeofday = local_timeofday + 86400._r8
            endif
          
            if (cosz == 1._r8) then
               azi_angle = 0._r8
               sinz = 0._r8
               solar_inc = 1._r8
            else
               sinz = sqrt(1._r8-cosz*cosz)
               azi_angle = (sin(lat(g))*cosz-sin(decl)) / (cos(lat(g))*sinz) !decl
               azi_angle = max(-1._r8,min(1._r8,azi_angle))
               azi_angle = acos(-azi_angle)
               if (local_timeofday >=43200._r8) then
                  azi_angle = 2._r8*pi - azi_angle
               endif
               azi_angle = pi / 2._r8 - azi_angle
               !write(iulog,*)  'lon180, ',azi_angle !test          
               solar_inc = 1._r8 + (sinz/cosz)*(cos(azi_angle)*sinsl_cosas(g)+sin(azi_angle)*sinsl_sinas(g))
            endif

            izen = int((cosz + 0.05_r8) / 0.15_r8)
            dzen1 = (cosz - (izen * 0.15_r8 - 0.05_r8)) / 0.15_r8
            dzen2 = 1._r8 - dzen1

            ftemp1 = coeff_dir(1,izen) * sky_view(g) + &
                     coeff_dir(2,izen) * solar_inc + coeff_dir(3,izen)
            ftemp2 = coeff_dir(1,izen+1) * sky_view(g) + &
                     coeff_dir(2,izen+1) * solar_inc + coeff_dir(3,izen+1)
            f_dir(p,1:numrad) = ftemp2 * dzen1 + ftemp1 * dzen2
            f_dir(p,1:numrad) = max(-1._r8,f_dir(p,1:numrad))

            ftemp1 = coeff_rdir(1,izen) * sky_view(g) + &
                     coeff_rdir(2,izen) * terrain_config(g) + coeff_rdir(3,izen)
            ftemp2 = coeff_rdir(1,izen+1) * sky_view(g) + &
                     coeff_rdir(2,izen+1) * terrain_config(g) + coeff_rdir(3,izen+1)
            f_rdir(p,1:numrad) = ftemp2 * dzen1 + ftemp1 * dzen2

            ftemp1 = coeff_dif(1,izen) * stdev_elev(g) + &
                     coeff_dif(2,izen) * sky_view(g) + &
                     coeff_dif(3,izen) * solar_inc + coeff_dif(4,izen)
            ftemp2 = coeff_dif(1,izen+1) * stdev_elev(g) + &
                     coeff_dif(2,izen+1) * sky_view(g) + &
                     coeff_dif(3,izen+1) * solar_inc + coeff_dif(4,izen+1)
            f_dif(p,1:numrad) = ftemp2 * dzen1 + ftemp1 * dzen2
            f_dif(p,1:numrad) = max(-1._r8,f_dif(p,1:numrad)) 

            ftemp1 = coeff_rdif(1,izen) * sky_view(g) + &
                     coeff_rdif(2,izen) * terrain_config(g) + coeff_rdif(3,izen)
            ftemp2 = coeff_rdif(1,izen+1) * sky_view(g) + &
                     coeff_rdif(2,izen+1) * terrain_config(g) + coeff_rdif(3,izen+1)
            f_rdif(p,1:numrad) = ftemp2 * dzen1 + ftemp1 * dzen2
            
            if(.not. is_veg) then
               do ib = 1, numrad
                 f_rdir(p,ib) = f_rdir(p,ib) * (albd(p,ib) / 0.1_r8)
                 f_rdif(p,ib) = f_rdif(p,ib) * (albi(p,ib) / 0.1_r8)

                 fd_prime = 1._r8 + f_dir(p,ib) + f_rdir(p,ib)
                 fi_prime = 1._r8 + f_dif(p,ib) + f_rdif(p,ib)

                 albd_adjust = fd_prime * albd(p,ib) - (fd_prime-1._r8)
                 albi_adjust = fi_prime * albi(p,ib) - (fi_prime-1._r8)

                 if (albd_adjust <= 0._r8) then 
                    albd_adjust = 0._r8
                    fd_prime = 1._r8 / (1._r8 - albd(p,ib))
                 endif

                 if (albi_adjust <= 0._r8) then
                    albi_adjust = 0._r8
                    fi_prime = 1._r8 / (1._r8 - albi(p,ib))
                 endif

                 albd(p,ib) = albd_adjust
                 fabd(p,ib) = fabd(p,ib) * fd_prime
                 ftdd(p,ib) = ftdd(p,ib) * fd_prime
                 ftid(p,ib) = ftid(p,ib) * fd_prime
                 fd_top_adjust(p,ib) = fd_prime

                 albi(p,ib) = albi_adjust
                 fabi(p,ib) = fabi(p,ib) * fi_prime
                 ftii(p,ib) = ftii(p,ib) * fi_prime
                 fi_top_adjust(p,ib) = fi_prime
               enddo
            endif
            
         endif
      enddo
	
     end associate

  end subroutine Albedo_TOP_Adjustment


end module SurfaceAlbedoMod
