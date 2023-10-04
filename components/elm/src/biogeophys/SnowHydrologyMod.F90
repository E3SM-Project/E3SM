module SnowHydrologyMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate snow hydrology.
  ! - Using as input aerosol deposition from atmosphere model calculate
  !   aerosol fluxes and masses in each layer - need for surface albedo calculation
  ! - Change of snow mass and the snow water onto soil
  ! - Change in snow layer thickness due to compaction
  ! - Combine snow layers less than a min thickness
  ! - Subdivide snow layers if they exceed maximum thickness
  ! - Construct snow/no-snow filters

  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use elm_varpar      , only : nlevsno
  use elm_varctl      , only : iulog, use_extrasnowlayers
  use elm_varcon      , only : namec, h2osno_max
  use atm2lndType     , only : atm2lnd_type
  use AerosolType     , only : aerosol_type
  use TopounitDataType, only : topounit_atmospheric_state
  use LandunitType    , only : lun_pp
  use ColumnType      , only : col_pp
  use ColumnDataType  , only : col_es, col_ef, col_ws, col_wf

  use timeinfoMod
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SnowWater                  ! Change of snow mass and the snow water onto soil
  public :: SnowCompaction             ! Change in snow layer thickness due to compaction
  public :: CombineSnowLayers          ! Combine snow layers less than a min thickness
  public :: DivideSnowLayers           ! Subdivide snow layers if they exceed maximum thickness
  public :: DivideExtraSnowLayers      ! Subdivide up to 16 snow layers if they exceed maximum thickness
  public :: BuildSnowFilter            ! Construct snow/no-snow filters
  public :: InitSnowLayers             ! Initialize cold-start snow layer thickness  
  public :: NewSnowBulkDensity         ! Compute bulk density of any newly-fallen snow
  public :: SnowCapping                ! Remove snow mass for capped columns
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: Combo            ! Returns the combined variables: dz, t, wliq, wice.
  !
  ! !PUBLIC DATA MEMBERS:
  !  Aerosol species indices:
  !  1= hydrophillic (bulk model) or within-ice (modal model) black carbon
  !  2= hydrophobic (bulk model) or external (modal model) black carbon
!  !  1= hydrophillic black carbon
!  !  2= hydrophobic black carbon
  !  3= hydrophilic organic carbon
  !  4= hydrophobic organic carbon
  !  5= dust species 1
  !  6= dust species 2
  !  7= dust species 3
  !  8= dust species 4
  !
  real(r8), public, parameter :: scvng_fct_mlt_bcphi = 0.20_r8 ! scavenging factor for hydrophillic BC inclusion in meltwater [frc]
  real(r8), public, parameter :: scvng_fct_mlt_bcpho = 0.03_r8 ! scavenging factor for hydrophobic BC inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_ocphi = 0.20_r8 ! scavenging factor for hydrophillic OC inclusion in meltwater [frc]
  real(r8), public, parameter :: scvng_fct_mlt_ocpho = 0.03_r8 ! scavenging factor for hydrophobic OC inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst1  = 0.02_r8 ! scavenging factor for dust species 1 inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst2  = 0.02_r8 ! scavenging factor for dust species 2 inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst3  = 0.01_r8 ! scavenging factor for dust species 3 inclusion in meltwater  [frc]
  real(r8), public, parameter :: scvng_fct_mlt_dst4  = 0.01_r8 ! scavenging factor for dust species 4 inclusion in meltwater  [frc]

  !$acc declare copyin(scvng_fct_mlt_bcphi)
  !$acc declare copyin(scvng_fct_mlt_bcpho)
  !$acc declare copyin(scvng_fct_mlt_ocphi)
  !$acc declare copyin(scvng_fct_mlt_ocpho)
  !$acc declare copyin(scvng_fct_mlt_dst1 )
  !$acc declare copyin(scvng_fct_mlt_dst2 )
  !$acc declare copyin(scvng_fct_mlt_dst3 )
  !$acc declare copyin(scvng_fct_mlt_dst4 )

  !-----------------------------------------------------------------------
  !H. Wang ++
  !  "Rfast" parameters used by Flanner et al (2012, ACP)
  !real(r8), public, parameter :: scvng_fct_mlt_bcphi = 1.00_r8   ! scavenging factor for hydrophillic BC inclusion in meltwater [frc]
  !real(r8), public, parameter :: scvng_fct_mlt_bcpho = 0.03_r8   ! scavenging factor for hydrophobic BC inclusion in meltwater [frc]
  !real(r8), public, parameter :: scvng_fct_mlt_ocphi = 1.00_r8   ! scavenging factor for hydrophillic OC inclusion in meltwater [frc]
  !real(r8), public, parameter :: scvng_fct_mlt_ocpho = 0.03_r8   ! scavenging factor for hydrophobic OC inclusion in meltwater [frc]
  !real(r8), public, parameter :: scvng_fct_mlt_dst1  = 0.03_r8   ! scavenging factor for dust species 1 inclusion in meltwater [frc]
  !real(r8), public, parameter :: scvng_fct_mlt_dst2  = 0.03_r8   ! scavenging factor for dust species 2 inclusion in meltwater [frc]
  !real(r8), public, parameter :: scvng_fct_mlt_dst3  = 0.03_r8   ! scavenging factor for dust species 3 inclusion in meltwater [frc]
  !real(r8), public, parameter :: scvng_fct_mlt_dst4  = 0.03_r8   ! scavenging factor for dust species 4 inclusion in meltwater [frc]
  ! ++
  
  ! Definition of snow pack vertical structure
  ! (adapted from CLMv5)
  ! Hardcoded maximum of 16 snowlayers
  ! The bottom layer has no limit on thickness, hence the last element of the dzmax_*
  ! arrays is 'huge'.
  real(r8), parameter :: dzmin16(16) = &      ! minimum of top snow layer
               (/ 0.010_r8, 0.015_r8, 0.025_r8, 0.055_r8, 0.115_r8, 0.235_r8, &
                  0.475_r8, 0.955_r8, 1.915_r8, 1.915_r8, 1.915_r8, 1.915_r8, &
                  1.915_r8, 1.915_r8, 1.915_r8, 1.915_r8 /)
  real(r8), parameter :: dzmax_l(16) = &    ! maximum thickness of layer when no layers beneath
               (/ 0.03_r8, 0.07_r8, 0.18_r8, 0.41_r8, 0.88_r8, 1.83_r8, &
                  3.74_r8, 7.57_r8, 15.24_r8, 15.24_r8, 15.24_r8, 15.24_r8, &
                  15.24_r8, 15.24_r8, 15.24_r8, huge(1._r8)  /)
  real(r8), parameter :: dzmax_u(16) = &    ! maximum thickness of layer when layers beneath
               (/ 0.02_r8, 0.05_r8, 0.11_r8, 0.23_r8, 0.47_r8, 0.95_r8, &
                  1.91_r8, 3.83_r8, 7.67_r8, 7.67_r8, 7.67_r8, 7.67_r8, & 
                  7.67_r8, 7.67_r8, 7.67_r8, huge(1._r8)  /)

contains

  !-----------------------------------------------------------------------
  subroutine SnowWater(bounds, &
       num_snowc, filter_snowc, num_nosnowc, filter_nosnowc, &
       atm2lnd_vars, aerosol_vars)
    !
    ! !DESCRIPTION:
    ! Evaluate the change of snow mass and the snow water onto soil.
    ! Water flow within snow is computed by an explicit and non-physical
    ! based scheme, which permits a part of liquid water over the holding
    ! capacity (a tentative value is used, i.e. equal to 0.033*porosity) to
    ! percolate into the underlying layer.  Except for cases where the
    ! porosity of one of the two neighboring layers is less than 0.05, zero
    ! flow is assumed. The water flow out of the bottom of the snow pack will
    ! participate as the input of the soil water and runoff.  This subroutine
    ! uses a filter for columns containing snow which must be constructed prior
    ! to being called.
    !
    ! !USES:
    use elm_varcon        , only : denh2o, denice, wimp, ssi
    use landunit_varcon   , only : istsoil
    use AerosolMod        , only : AerosolFluxes
    use elm_varctl        , only : use_vsfm
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds
    integer               , intent(in)    :: num_snowc         ! number of snow points in column filter
    integer               , intent(in)    :: filter_snowc(:)   ! column filter for snow points
    integer               , intent(in)    :: num_nosnowc       ! number of non-snow points in column filter
    integer               , intent(in)    :: filter_nosnowc(:) ! column filter for non-snow points
    type(atm2lnd_type)    , intent(in)    :: atm2lnd_vars
    type(aerosol_type)    , intent(inout) :: aerosol_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: g                                                  ! gridcell loop index
    integer  :: c, j, fc, l                                        ! do loop/array indices
    real(r8) :: qin         (1:num_snowc)                       ! water flow into the elmement (mm/s)
    real(r8) :: qout        (1:num_snowc)                      ! water flow out of the elmement (mm/s)
    real(r8) :: qin_bc_phi  (1:num_snowc)              ! flux of hydrophilic BC into   layer [kg]
    real(r8) :: qout_bc_phi (1:num_snowc)              ! flux of hydrophilic BC out of layer [kg]
    real(r8) :: qin_bc_pho  (1:num_snowc)              ! flux of hydrophobic BC into   layer [kg]
    real(r8) :: qout_bc_pho (1:num_snowc)              ! flux of hydrophobic BC out of layer [kg]
    real(r8) :: qin_oc_phi  (1:num_snowc)              ! flux of hydrophilic OC into   layer [kg]
    real(r8) :: qout_oc_phi (1:num_snowc)              ! flux of hydrophilic OC out of layer [kg]
    real(r8) :: qin_oc_pho  (1:num_snowc)              ! flux of hydrophobic OC into   layer [kg]
    real(r8) :: qout_oc_pho (1:num_snowc)              ! flux of hydrophobic OC out of layer [kg]
    real(r8) :: qin_dst1    (1:num_snowc)              ! flux of dust species 1 into   layer [kg]
    real(r8) :: qout_dst1   (1:num_snowc)              ! flux of dust species 1 out of layer [kg]
    real(r8) :: qin_dst2    (1:num_snowc)              ! flux of dust species 2 into   layer [kg]
    real(r8) :: qout_dst2   (1:num_snowc)              ! flux of dust species 2 out of layer [kg]
    real(r8) :: qin_dst3    (1:num_snowc)              ! flux of dust species 3 into   layer [kg]
    real(r8) :: qout_dst3   (1:num_snowc)              ! flux of dust species 3 out of layer [kg]
    real(r8) :: qin_dst4    (1:num_snowc)              ! flux of dust species 4 into   layer [kg]
    real(r8) :: qout_dst4   (1:num_snowc)              ! flux of dust species 4 out of layer [kg]
    real(r8) :: wgdif                                              ! ice mass after minus sublimation
    real(r8) :: vol_liq(1:num_snowc,-nlevsno+1:0)      ! partial volume of liquid water in layer
    real(r8) :: vol_ice(1:num_snowc,-nlevsno+1:0)      ! partial volume of ice lens in layer
    real(r8) :: eff_porosity(1:num_snowc,-nlevsno+1:0) ! effective porosity = porosity - vol_ice
    real(r8) :: mss_liqice(1:num_snowc,-nlevsno+1:0)   ! mass of liquid+ice in a layer
    !-----------------------------------------------------------------------
    !mgf++
    real(r8) :: refrzsnow                   ! re-frozen snow [kg m-2]
    real(r8) :: subsnow                     ! sublimated snow [kg m-2]
    real(r8) :: frc_refrz                   ! fraction of layer mass that is re-frozen snow [frc]
    real(r8) :: frc_sub                     ! fraction of layer mass that has sublimated [frc]
    real(r8) :: frc_transfer                ! frc_refrz + frc_sub
    real(r8) :: dm_int                      ! mass transfer [kg]
    !mgf--


    associate(                                                 &
         dz             => col_pp%dz                            , & ! Input:  [real(r8) (:,:) ] layer depth (m)
         snl            => col_pp%snl                           , & ! Input:  [integer  (:)   ] number of snow layers

         do_capsnow     => col_ws%do_capsnow    , & ! Input:  [logical  (:)   ] true => do snow capping
         frac_sno_eff   => col_ws%frac_sno_eff  , & ! Input:  [real(r8) (:)   ] eff. fraction of ground covered by snow (0 to 1)
         frac_sno       => col_ws%frac_sno      , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
         h2osno         => col_ws%h2osno        , & ! Input:  [real(r8) (:)   ] snow water (mm H2O)
         int_snow       => col_ws%int_snow      , & ! Output: [real(r8) (:)   ] integrated snowfall [mm]
         h2osoi_ice     => col_ws%h2osoi_ice    , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)
         h2osoi_liq     => col_ws%h2osoi_liq    , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)

         qflx_snomelt   => col_wf%qflx_snomelt   , & ! Input:  [real(r8) (:)   ] snow melt (mm H2O /s)
         qflx_rain_grnd => col_wf%qflx_rain_grnd , & ! Input:  [real(r8) (:)   ] rain on ground after interception (mm H2O/s) [+]
         qflx_sub_snow  => col_wf%qflx_sub_snow  , & ! Input:  [real(r8) (:)   ] sublimation rate from snow pack (mm H2O /s) [+]
         qflx_dew_snow  => col_wf%qflx_dew_snow  , & ! Input:  [real(r8) (:)   ] surface dew added to snow pack (mm H2O /s) [+]
         qflx_evap_grnd => col_wf%qflx_evap_grnd , & ! Input:  [real(r8) (:)   ] ground surface evaporation rate (mm H2O/s) [+]
         qflx_dew_grnd  => col_wf%qflx_dew_grnd  , & ! Input:  [real(r8) (:)   ] ground surface dew formation (mm H2O /s) [+]
         qflx_snow_melt => col_wf%qflx_snow_melt , & ! Output: [real(r8) (:)   ] net snow melt
         qflx_top_soil  => col_wf%qflx_top_soil  , & ! Output: [real(r8) (:)   ] net water input into soil from top (mm/s)
         ! qflx_snofrz_lyr => cwf%qflx_snofrz_lyr     , & ! HW+++ snow freezing rate (col,lyr) [kg m-2 s-1]

         mflx_neg_snow_col_1d =>  col_wf%mflx_neg_snow_1d , & ! Output:  [real(r8) (:)   ]  mass flux from top soil layer due to negative water content in snow layers (kg H2O /s)

         mss_bcphi      => aerosol_vars%mss_bcphi_col        , & ! Output: [real(r8) (:,:) ] hydrophillic BC mass in snow (col,lyr) [kg]
         mss_bcpho      => aerosol_vars%mss_bcpho_col        , & ! Output: [real(r8) (:,:) ] hydrophobic  BC mass in snow (col,lyr) [kg]
         mss_ocphi      => aerosol_vars%mss_ocphi_col        , & ! Output: [real(r8) (:,:) ] hydrophillic OC mass in snow (col,lyr) [kg]
         mss_ocpho      => aerosol_vars%mss_ocpho_col        , & ! Output: [real(r8) (:,:) ] hydrophobic  OC mass in snow (col,lyr) [kg]
         mss_dst1       => aerosol_vars%mss_dst1_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 1 in snow (col,lyr) [kg]
         mss_dst2       => aerosol_vars%mss_dst2_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 2 in snow (col,lyr) [kg]
         mss_dst3       => aerosol_vars%mss_dst3_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 3 in snow (col,lyr) [kg]
         mss_dst4       => aerosol_vars%mss_dst4_col         , & ! Output: [real(r8) (:,:) ] mass of dust species 4 in snow (col,lyr) [kg]

         begc           => bounds%begc                       , &
         endc           => bounds%endc                         &
         )
    !$acc enter data create(&
    !$acc qin(:), &
    !$acc qout(:), &
    !$acc qin_bc_phi(:), &
    !$acc qout_bc_phi(:), &
    !$acc qin_bc_pho(:), &
    !$acc qout_bc_pho(:), &
    !$acc qin_oc_phi(:), &
    !$acc qout_oc_phi(:), &
    !$acc qin_oc_pho(:), &
    !$acc qout_oc_pho(:), &
    !$acc qin_dst1(:), &
    !$acc qout_dst1(:), &
    !$acc qin_dst2(:), &
    !$acc qout_dst2(:), &
    !$acc qin_dst3(:), &
    !$acc qout_dst3(:), &
    !$acc qin_dst4(:), &
    !$acc qout_dst4(:), &
    !$acc vol_liq(:,:), &
    !$acc vol_ice(:,:), &
    !$acc eff_porosity(:,:), &
    !$acc mss_liqice(:,:))



      ! Renew the mass of ice lens (h2osoi_ice) and liquid (h2osoi_liq) in the
      ! surface snow layer resulting from sublimation (frost) / evaporation (condense)

      !TODO: check where else this is needed 
      !mflx_neg_snow_col_1d(:) = 0._r8

      !$acc parallel loop independent gang vector default(present)
      do fc = 1,num_snowc
         c = filter_snowc(fc)
         l=col_pp%landunit(c)

         if (do_capsnow(c) .and. .not. use_extrasnowlayers) then
            wgdif = h2osoi_ice(c,snl(c)+1) - frac_sno_eff(c)*qflx_sub_snow(c)*dtime_mod
            h2osoi_ice(c,snl(c)+1) = wgdif
            if (wgdif < 0._r8) then
               h2osoi_ice(c,snl(c)+1) = 0._r8
               h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) + wgdif
            end if
            h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) &
                 - frac_sno_eff(c)*qflx_evap_grnd(c) * dtime_mod
         else
            wgdif = h2osoi_ice(c,snl(c)+1) &
                 + frac_sno_eff(c) * (qflx_dew_snow(c) - qflx_sub_snow(c)) * dtime_mod
            h2osoi_ice(c,snl(c)+1) = wgdif
            if (wgdif < 0._r8) then
               h2osoi_ice(c,snl(c)+1) = 0._r8
               h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) + wgdif
            end if
            h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) +  &
                 frac_sno_eff(c) * (qflx_rain_grnd(c) + qflx_dew_grnd(c) &
                 - qflx_evap_grnd(c)) * dtime_mod
         end if
         ! if negative, reduce deeper layer's liquid water content sequentially
         if(h2osoi_liq(c,snl(c)+1) < 0._r8) then
            !$acc loop seq
            do j = snl(c)+1, 1
               wgdif=h2osoi_liq(c,j)
               if (wgdif >= 0._r8) exit
               h2osoi_liq(c,j) = 0._r8
               if (.not.(j+1 > 0 .and. use_vsfm)) then
                  h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + wgdif
               else
                  mflx_neg_snow_col_1d(c-bounds%begc+1) = wgdif/dtime_mod
               endif
            enddo
         end if
      end do

      ! Porosity and partial volume

      !$acc parallel loop independent gang vector default(present) collapse(2) 
      do j = -nlevsno+1, 0
         do fc = 1, num_snowc
            c = filter_snowc(fc)
            if (j >= snl(c)+1) then
               ! need to scale dz by frac_sno to convert to grid cell average depth
               vol_ice(fc,j)      = min(1._r8, h2osoi_ice(c,j)/(dz(c,j)*frac_sno_eff(c)*denice))
               eff_porosity(fc,j) = 1._r8 - vol_ice(fc,j)
               vol_liq(fc,j)      = min(eff_porosity(fc,j),h2osoi_liq(c,j)/(dz(c,j)*frac_sno_eff(c)*denh2o))
            end if
         end do
      end do

      ! Capillary forces within snow are usually two or more orders of magnitude
      ! less than those of gravity. Only gravity terms are considered.
      ! the genernal expression for water flow is "K * ss**3", however,
      ! no effective parameterization for "K".  Thus, a very simple consideration
      ! (not physically based) is introduced:
      ! when the liquid water of layer exceeds the layer's holding
      ! capacity, the excess meltwater adds to the underlying neighbor layer.

      ! Also compute aerosol fluxes through snowpack in this loop:
      ! 1) compute aerosol mass in each layer
      ! 2) add aerosol mass flux from above layer to mass of this layer
      ! 3) qout_xxx is mass flux of aerosol species xxx out bottom of
      !    layer in water flow, proportional to (current) concentration
      !    of aerosol in layer multiplied by a scavenging ratio.
      ! 4) update mass of aerosol in top layer, accordingly
      ! 5) update mass concentration of aerosol accordingly

      !$acc parallel loop independent gang vector default(present)
      do fc = 1,num_snowc
         c = filter_snowc(fc)
         qin(fc)         = 0._r8
         qin_bc_phi(fc) = 0._r8
         qin_bc_pho(fc) = 0._r8
         qin_oc_phi(fc) = 0._r8
         qin_oc_pho(fc) = 0._r8
         qin_dst1(fc) = 0._r8
         qin_dst2(fc) = 0._r8
         qin_dst3(fc) = 0._r8
         qin_dst4(fc) = 0._r8
      end do

      !NOTE : look into qout race condition?

      !$acc parallel loop independent gang vector default(present)
      do fc = 1, num_snowc
         c = filter_snowc(fc)
         !$acc loop seq 
         do j = -nlevsno+1, 0
            if (j >= snl(c)+1) then

               h2osoi_liq(c,j) = h2osoi_liq(c,j) + qin(fc)  

               mss_bcphi(c,j) = mss_bcphi(c,j) + qin_bc_phi(fc)
               mss_bcpho(c,j) = mss_bcpho(c,j) + qin_bc_pho(fc)
               mss_ocphi(c,j) = mss_ocphi(c,j) + qin_oc_phi(fc)
               mss_ocpho(c,j) = mss_ocpho(c,j) + qin_oc_pho(fc)

               mss_dst1(c,j)  = mss_dst1(c,j) + qin_dst1(fc)
               mss_dst2(c,j)  = mss_dst2(c,j) + qin_dst2(fc)
               mss_dst3(c,j)  = mss_dst3(c,j) + qin_dst3(fc)
               mss_dst4(c,j)  = mss_dst4(c,j) + qin_dst4(fc)

               if (j <= -1) then
                  ! No runoff over snow surface, just ponding on surface
                  if (eff_porosity(fc,j) < wimp .OR. eff_porosity(fc,j+1) < wimp) then
                     qout(fc) = 0._r8
                  else
                     ! dz must be scaled by frac_sno to obtain gridcell average value
                     qout(fc) = max(0._r8,(vol_liq(fc,j) &
                          - ssi*eff_porosity(fc,j))*dz(c,j)*frac_sno_eff(c))
                     qout(fc) = min(qout(fc),(1._r8-vol_ice(fc,j+1) &
                          - vol_liq(fc,j+1))*dz(c,j+1)*frac_sno_eff(c))
                  end if
               else
                  qout(fc) = max(0._r8,(vol_liq(fc,j) &
                       - ssi*eff_porosity(fc,j))*dz(c,j)*frac_sno_eff(c))
               end if
               qout(fc) = qout(fc)*1000._r8
               h2osoi_liq(c,j) = h2osoi_liq(c,j) - qout(fc)
               qin(fc) = qout(fc)

               ! mass of ice+water: in extremely rare circumstances, this can
               ! be zero, even though there is a snow layer defined. In
               ! this case, set the mass to a very small value to
               ! prevent division by zero.

               mss_liqice(fc,j) = h2osoi_liq(c,j)+h2osoi_ice(c,j)
               if (mss_liqice(fc,j) < 1E-30_r8) then
                  mss_liqice(fc,j) = 1E-30_r8
               endif

               ! BCPHI:
               ! 1. flux with meltwater:
               qout_bc_phi(fc) = qout(fc)*scvng_fct_mlt_bcphi*(mss_bcphi(c,j)/mss_liqice(fc,j))
               if (qout_bc_phi(fc) > mss_bcphi(c,j)) then
                  qout_bc_phi(fc) = mss_bcphi(c,j)
               endif
               mss_bcphi(c,j) = mss_bcphi(c,j) - qout_bc_phi(fc)
               qin_bc_phi(fc) = qout_bc_phi(fc)

               ! BCPHO:
               ! 1. flux with meltwater:
               qout_bc_pho(fc) = qout(fc)*scvng_fct_mlt_bcpho*(mss_bcpho(c,j)/mss_liqice(fc,j))
               if (qout_bc_pho(fc) > mss_bcpho(c,j)) then
                  qout_bc_pho(fc) = mss_bcpho(c,j)
               endif
               mss_bcpho(c,j) = mss_bcpho(c,j) - qout_bc_pho(fc)
               qin_bc_pho(fc) = qout_bc_pho(fc)

               ! OCPHI:
               ! 1. flux with meltwater:
               qout_oc_phi(fc) = qout(fc)*scvng_fct_mlt_ocphi*(mss_ocphi(c,j)/mss_liqice(fc,j))
               if (qout_oc_phi(fc) > mss_ocphi(c,j)) then
                  qout_oc_phi(fc) = mss_ocphi(c,j)
               endif
               mss_ocphi(c,j) = mss_ocphi(c,j) - qout_oc_phi(fc)
               qin_oc_phi(fc) = qout_oc_phi(fc)

               ! OCPHO:
               ! 1. flux with meltwater:
               qout_oc_pho(fc) = qout(fc)*scvng_fct_mlt_ocpho*(mss_ocpho(c,j)/mss_liqice(fc,j))
               if (qout_oc_pho(fc) > mss_ocpho(c,j)) then
                  qout_oc_pho(fc) = mss_ocpho(c,j)
               endif
               mss_ocpho(c,j) = mss_ocpho(c,j) - qout_oc_pho(fc)
               qin_oc_pho(fc) = qout_oc_pho(fc)

               ! DUST 1:
               ! 1. flux with meltwater:
               qout_dst1(fc) = qout(fc)*scvng_fct_mlt_dst1*(mss_dst1(c,j)/mss_liqice(fc,j))
               if (qout_dst1(fc) > mss_dst1(c,j)) then
                  qout_dst1(fc) = mss_dst1(c,j)
               endif
               mss_dst1(c,j) = mss_dst1(c,j) - qout_dst1(fc)
               qin_dst1(fc) = qout_dst1(fc)

               ! DUST 2:
               ! 1. flux with meltwater:
               qout_dst2(fc) = qout(fc)*scvng_fct_mlt_dst2*(mss_dst2(c,j)/mss_liqice(fc,j))
               if (qout_dst2(fc) > mss_dst2(c,j)) then
                  qout_dst2(fc) = mss_dst2(c,j)
               endif
               mss_dst2(c,j) = mss_dst2(c,j) - qout_dst2(fc)
               qin_dst2(fc) = qout_dst2(fc)

               ! DUST 3:
               ! 1. flux with meltwater:
               qout_dst3(fc) = qout(fc)*scvng_fct_mlt_dst3*(mss_dst3(c,j)/mss_liqice(fc,j))
               if (qout_dst3(fc) > mss_dst3(c,j)) then
                  qout_dst3(fc) = mss_dst3(c,j)
               endif
               mss_dst3(c,j) = mss_dst3(c,j) - qout_dst3(fc)
               qin_dst3(fc) = qout_dst3(fc)

               ! DUST 4:
               ! 1. flux with meltwater:
               qout_dst4(fc) = qout(fc)*scvng_fct_mlt_dst4*(mss_dst4(c,j)/mss_liqice(fc,j))
               if (qout_dst4(fc) > mss_dst4(c,j)) then
                  qout_dst4(fc) = mss_dst4(c,j)
               endif
               mss_dst4(c,j) = mss_dst4(c,j) - qout_dst4(fc)
               qin_dst4(fc) = qout_dst4(fc)

            end if
         end do
      end do

      ! Compute aerosol fluxes through snowpack and aerosol deposition fluxes into top layer
      call AerosolFluxes(bounds, num_snowc, filter_snowc, &
           atm2lnd_vars, aerosol_vars)

      ! Adjust layer thickness for any water+ice content changes in excess of previous
      ! layer thickness. Strictly speaking, only necessary for top snow layer, but doing
      ! it for all snow layers will catch problems with older initial files.
      ! Layer interfaces (zi) and node depths (z) do not need adjustment here because they
      ! are adjusted in CombineSnowLayers and are not used up to that point.

      !$acc parallel loop independent gang vector default(present) collapse(2) 
      do j = -nlevsno+1, 0
         do fc = 1, num_snowc
            c = filter_snowc(fc)
            if (j >= snl(c)+1) then
               dz(c,j) = max(dz(c,j),h2osoi_liq(c,j)/denh2o + h2osoi_ice(c,j)/denice)
            end if
         end do
      end do

      !$acc parallel loop independent gang vector default(present)
      do fc = 1, num_snowc
         c = filter_snowc(fc)
         ! Qout from snow bottom
         qflx_snow_melt(c) = qflx_snow_melt(c) + (qout(fc) / dtime_mod)

         qflx_top_soil(c) = (qout(fc) / dtime_mod) &
              + (1.0_r8 - frac_sno_eff(c)) * qflx_rain_grnd(c)
         int_snow(c) = int_snow(c) + frac_sno_eff(c) &
                       * (qflx_dew_snow(c) + qflx_dew_grnd(c) + qflx_rain_grnd(c)) * dtime_mod
      end do

      !$acc parallel loop independent gang vector default(present)
      do fc = 1, num_nosnowc
         c = filter_nosnowc(fc)
         qflx_snow_melt(c) = qflx_snomelt(c)

         qflx_top_soil(c) = qflx_rain_grnd(c) + qflx_snomelt(c)
         ! reset accumulated snow when no snow present
         if (h2osno(c) <= 0) int_snow(c) = 0.
         if (h2osno(c) <= 0) frac_sno(c) = 0.
      end do

#ifdef MODAL_AER
    !mgf++
    !
    ! Transfer BC and OC from the within-ice state to the external
    ! state based on snow sublimation and re-freezing of liquid water.
    ! Re-freezing effect is inactived by default because of
    ! uncertainty in how this process operates.
    !$acc parallel loop independent gang vector default(present) collapse(2) 
    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             !! snow that has re-frozen [kg/m2]
             !refrzsnow = max(0._r8, (qflx_snofrz_lyr(c,j)*dtime_mod))
             !
             !! fraction of layer mass that is re-frozen
             !if ((h2osoi_liq(c,j) + h2osoi_ice(c,j)) > 0._r8) then
             !   frc_refrz = refrzsnow / (h2osoi_liq(c,j) + h2osoi_ice(c,j))
             !else
             !   frc_refrz = 0._r8
             !endif

             if (j == snl(c)+1) then
                ! snow that has sublimated [kg/m2] (top layer only)
                subsnow = max(0._r8, (qflx_sub_snow(c)*dtime_mod))

                ! fraction of layer mass that has sublimated:
                if ((h2osoi_liq(c,j) + h2osoi_ice(c,j)) > 0._r8) then
                   frc_sub = subsnow / (h2osoi_liq(c,j) + h2osoi_ice(c,j))
                else
                   frc_sub = 0._r8
                endif
             else
                ! prohibit sublimation effect to operate on
                ! sub-surface layers:
                frc_sub = 0._r8
             endif

             ! fraction of layer mass transformed (sublimation only)
             !frc_transfer = frc_refrz + frc_sub
             frc_transfer = frc_sub

             ! cap the fraction at 1
             if (frc_transfer > 1._r8) then
                frc_transfer = 1._r8
             endif

             ! transfer proportionate mass of BC and OC:
             dm_int         = mss_bcphi(c,j)*frc_transfer
             mss_bcphi(c,j) = mss_bcphi(c,j) - dm_int
             mss_bcpho(c,j) = mss_bcpho(c,j) + dm_int

             dm_int         = mss_ocphi(c,j)*frc_transfer
             mss_ocphi(c,j) = mss_ocphi(c,j) - dm_int
             mss_ocpho(c,j) = mss_ocpho(c,j) + dm_int

          end if
       end do
    end do
    !mgf--
#endif




    !$acc exit data delete(&
    !$acc qin(:), &
    !$acc qout(:), &
    !$acc qin_bc_phi(:), &
    !$acc qout_bc_phi(:), &
    !$acc qin_bc_pho(:), &
    !$acc qout_bc_pho(:), &
    !$acc qin_oc_phi(:), &
    !$acc qout_oc_phi(:), &
    !$acc qin_oc_pho(:), &
    !$acc qout_oc_pho(:), &
    !$acc qin_dst1(:), &
    !$acc qout_dst1(:), &
    !$acc qin_dst2(:), &
    !$acc qout_dst2(:), &
    !$acc qin_dst3(:), &
    !$acc qout_dst3(:), &
    !$acc qin_dst4(:), &
    !$acc qout_dst4(:), &
    !$acc vol_liq(:,:), &
    !$acc vol_ice(:,:), &
    !$acc eff_porosity(:,:), &
    !$acc mss_liqice(:,:))

    end associate

   end subroutine SnowWater

   !-----------------------------------------------------------------------
   subroutine SnowCompaction(bounds, num_snowc, filter_snowc, &
         top_as, dtime)
     !
     ! !DESCRIPTION:
     ! Determine the change in snow layer thickness due to compaction and
     ! settling.
     ! Three metamorphisms of changing snow characteristics are implemented,
     ! i.e., destructive, overburden, and melt. The treatments of the former
     ! two are from SNTHERM.89 and SNTHERM.99 (1991, 1999). The contribution
     ! due to melt metamorphism is simply taken as a ratio of snow ice
     ! fraction after the melting versus before the melting.
     !
     ! !USES:
     use elm_varcon      , only : denice, denh2o, tfrz, rpi, grav, rgas
     use landunit_varcon , only : istice_mec, istdlak, istsoil, istcrop
     use elm_varctl      , only : subgridflag
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in) :: bounds
     integer                , intent(in) :: num_snowc       ! number of column snow points in column filter
     integer                , intent(in) :: filter_snowc(:) ! column filter for snow points
     type(topounit_atmospheric_state), intent(in) :: top_as
     real(r8), intent(in)  :: dtime
     !
     ! !LOCAL VARIABLES:
     integer :: j, l, c, fc, t                   ! indices
     ! parameters
     real(r8), parameter :: c2 = 23.e-3_r8       ! [m3/kg]
     real(r8), parameter :: c3 = 2.777e-6_r8     ! [1/s]
     real(r8), parameter :: c3_ams = 5.8e-7_r8   ! Schneider et al., 2020 [1/s]
     real(r8), parameter :: c4 = 0.04_r8         ! [1/K]
     real(r8), parameter :: c5 = 2.0_r8          !
     real(r8), parameter :: dm = 100.0_r8        ! Upper Limit on Destructive Metamorphism Compaction [kg/m3]
     real(r8), parameter :: rho_dm = 150.0_r8    ! Upper limit on destructive metamorphism compaction [kg/m3] (Anderson, 1976; Schneider et al., 2020)
     real(r8), parameter :: eta0 = 9.e+5_r8      ! The Viscosity Coefficient Eta0 [kg-s/m2]
     real(r8), parameter :: k_creep_snow = 1.4e-9_r8 ! Creep coefficient for snow (bi < 550 kg / m3) [m3-s/kg]
     real(r8), parameter :: k_creep_firn = 1.2e-9_r8 ! Creep coefficient for firn (bi > 550 kg / m3)
     !
     real(r8) :: p_gls                           ! grain load stress [kg / m-s2]
     real(r8) :: burden (1:num_snowc) ! pressure of overlying snow [kg/m2]
     real(r8) :: zpseudo(1:num_snowc)! wind drift compaction / pseudo depth
     logical  :: mobile (1:num_snowc) ! current snow layer is mobile, i.e. susceptible to wind drift
     real(r8) :: snw_ssa                         ! Equivalent snow specific surface area [m2/kg] 
     real(r8) :: ddz1_fresh                      ! Rate of settling of dendritic snowpack (Lehning et al., 2002) [1/s]
     real(r8) :: ddz1                            ! Rate of settling of snowpack due to destructive metamorphism.
     real(r8) :: ddz2                            ! Rate of compaction of snowpack due to overburden. [1/s]
     real(r8) :: ddz3                            ! Rate of compaction of snowpack due to melt [1/s]
     real(r8) :: ddz4                            ! Rate of compaction of snowpack due to wind drift [1/s]
     real(r8) :: dexpf                           ! expf=exp(-c4*(273.15-t_soisno)).
     real(r8) :: fi                              ! Fraction of ice relative to the total water content at current time step
     real(r8) :: td                              ! tfrz - t_soisno [K]
     real(r8) :: pdzdtc                          ! Nodal rate of change in fractional-thickness due to compaction [fraction/s]
     real(r8) :: void                            ! void (1 - vol_ice - vol_liq)
     real(r8) :: wx                              ! water mass (ice+liquid) [kg/m2]
     real(r8) :: bi                              ! partial density of ice [kg/m3]
     real(r8) :: wsum                            ! snowpack total water mass (ice+liquid) [kg/m2]
     real(r8) :: fsno_melt
     real(r8) :: burden_noextra 
     !-----------------------------------------------------------------------

     associate(                                              &
          snl          => col_pp%snl                          , & ! Input:  [integer (:)    ] number of snow layers
          n_melt       => col_pp%n_melt                       , & ! Input:  [real(r8) (:)   ] SCA shape parameter
          ltype        => lun_pp%itype                        , & ! Input:  [integer (:)    ] landunit type

          forc_wind    => top_as%windbot                 , & ! Input:  [real(r8) (:) ]  atmospheric wind speed (m/s)
          t_soisno     => col_es%t_soisno    , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)
          imelt        => col_ef%imelt       , & ! Input:  [integer (:,:)  ] flag for melting (=1), freezing (=2), Not=0

          snow_depth   => col_ws%snow_depth   , & ! Input:  [real(r8) (:)   ] snow height (m)
          frac_sno     => col_ws%frac_sno_eff , & ! Input:  [real(r8) (:)   ] snow covered fraction
          swe_old      => col_ws%swe_old      , & ! Input:  [real(r8) (:,:) ] initial swe values
          int_snow     => col_ws%int_snow     , & ! Input:  [real(r8) (:)   ] integrated snowfall [mm]
          frac_iceold  => col_ws%frac_iceold  , & ! Input:  [real(r8) (:,:) ] fraction of ice relative to the tot water
          h2osoi_ice   => col_ws%h2osoi_ice   , & ! Input:  [real(r8) (:,:) ] ice lens (kg/m2)                       
          h2osoi_liq   => col_ws%h2osoi_liq   , & ! Input:  [real(r8) (:,:) ] liquid water (kg/m2)                   
          snw_rds      => col_ws%snw_rds      , & ! Output: [real(r8) (:,:) ] effective snow grain radius (col,lyr) [microns, m^-6]
          dz           => col_pp%dz                             & ! Output: [real(r8) (: ,:) ] layer depth (m)                        
          )
     !$acc enter data create(&
     !$acc burden(:), &
     !$acc zpseudo(:), &
     !$acc mobile(:), &
     !$acc ddz1_fresh, &
     !$acc ddz1, &
     !$acc ddz3, &
     !$acc burden_noextra)




       ! Begin calculation - note that the following column loops are only invoked if snl(c) < 0
          if (use_extrasnowlayers) then
          do fc = 1, num_snowc
             burden(fc) = 0._r8
             zpseudo(fc) = 0._r8
             mobile(fc) = .true.
          end do
      !  else
      !    !$acc parallel loop independent gang vector default(present)
      !    do fc = 1, num_snowc
      !       burden(fc) = 0._r8
      !    end do 
       end if

       if( use_extrasnowlayers) then 
         do j = -nlevsno+1, 0
            do fc = 1, num_snowc
               c = filter_snowc(fc)
               t = col_pp%topounit(c)
               l = col_pp%landunit(c)
               if (j >= snl(c)+1) then

                  wx = (h2osoi_ice(c,j) + h2osoi_liq(c,j))
                  void = 1._r8 - (h2osoi_ice(c,j)/denice + h2osoi_liq(c,j)/denh2o)&
                       /(frac_sno(c) * dz(c,j))
                  
                  ! Allow compaction only for non-saturated node and higher ice lens node.
                  if (void > 0.001_r8 .and. h2osoi_ice(c,j) > .1_r8) then

                     bi = h2osoi_ice(c,j) / (frac_sno(c) * dz(c,j))
                     fi = h2osoi_ice(c,j) / wx
                     td = tfrz-t_soisno(c,j)
                     dexpf = exp(-c4*td)

                     ! Settling as a result of destructive metamorphism
                     ddz1_fresh = (-grav * (burden(fc) + wx/2._r8)) / &
                                   (0.007_r8 * bi**(4.75_r8 + td/40._r8))
                      snw_ssa = 3.e6_r8 / (denice * snw_rds(c,j))
                      if (snw_ssa < 50._r8) then
                          ddz1_fresh = ddz1_fresh * exp(-46.e-2_r8 * (50._r8 - snw_ssa))
                      endif
                      ddz1 = -c3_ams*dexpf
                      if (bi > rho_dm) ddz1 = ddz1*exp(-46.0e-3_r8*(bi-rho_dm))
                      ddz1 = ddz1 + ddz1_fresh
                      ! Liquid water term

                     if (h2osoi_liq(c,j) > 0.01_r8*dz(c,j)*frac_sno(c)) ddz1=ddz1*c5

                     ! Compaction due to overburden
                     p_gls = max(denice / bi, 1._r8) * grav * (burden(fc) + wx/2._r8)
                      if (bi <= 550._r8) then ! Low density, i.e. snow
                         ddz2 = (-k_creep_snow * (max(denice / bi, 1._r8) - 1._r8) * &
                                 exp(-60.e6_r8 / (rgas * t_soisno(c,j))) * p_gls) / &
                                 (snw_rds(c,j) * 1.e-6_r8 * snw_rds(c,j) * 1.e-6_r8) - &
                                 2.02e-10_r8
                      else ! High density, i.e. firn
                         ddz2 = (-k_creep_firn * (max(denice / bi, 1._r8) - 1._r8) * &
                                 exp(-60.e6_r8 / (rgas * t_soisno(c,j))) * p_gls) / &
                                 (snw_rds(c,j) * 1.e-6_r8 * snw_rds(c,j) * 1.e-6_r8) - &
                                 2.7e-11_r8
                      endif

                     ! Compaction occurring during melt
                     if (imelt(c,j) == 1) then
                        if(subgridflag==1 .and. (ltype(col_pp%landunit(c)) == istsoil .or. ltype(col_pp%landunit(c)) == istcrop)) then
                           ! first term is delta mass over mass
                           ddz3 = max(0._r8,min(1._r8,(swe_old(c,j) - wx)/wx))

                           ! 2nd term is delta fsno over fsno, allowing for negative values for ddz3
                           if ((swe_old(c,j) - wx) > 0._r8) then
                              wsum = sum(h2osoi_liq(c,snl(c)+1:0)+h2osoi_ice(c,snl(c)+1:0))
                              fsno_melt = 1. - (acos(2.*min(1._r8,wsum/int_snow(c)) - 1._r8)/rpi)**(n_melt(c))

                              ddz3 = ddz3 - max(0._r8,(fsno_melt - frac_sno(c))/frac_sno(c))
                           endif
                           ddz3 = -1._r8/dtime * ddz3
                        else
                           ddz3 = - 1._r8/dtime * max(0._r8,(frac_iceold(c,j) - fi)/frac_iceold(c,j))
                        endif
                     else
                        ddz3 = 0._r8
                     end if

                     ! Compaction occurring due to wind drift
                     call WindDriftCompaction( &
                     bi = bi, &
                     forc_wind = forc_wind(t), &
                     dz = dz(c,j), &
                     zpseudo = zpseudo(fc), &
                     mobile = mobile(fc), &
                     compaction_rate = ddz4)

                      ! Time rate of fractional change in dz (units of s-1)

                   pdzdtc = ddz1 + ddz2 + ddz3 + ddz4

                   ! The change in dz due to compaction
                   ! Limit compaction to be no greater than fully saturated layer thickness

                   dz(c,j) = max(dz(c,j) * (1._r8+pdzdtc*dtime),(h2osoi_ice(c,j)/denice+ h2osoi_liq(c,j)/denh2o)/frac_sno(c))
                
               else ! from CLMv5
                   ! saturated node is immobile
                   !
                   ! This is only needed if wind_dependent_snow_density is true, but it's
                   ! simplest just to update mobile always
                   mobile(fc) = .false.
               end if

                  ! Pressure of overlying snow

                  burden(fc) = burden(fc) + wx
               end if 
            end do
         end do  
      else 

         !NOTE: revisit these loops 
         !$acc parallel loop independent gang vector default(present) private(burden_noextra)
         do fc = 1, num_snowc
            c = filter_snowc(fc)
            burden_noextra = 0._r8 
            !$acc loop seq 
            do j = snl(c)+1, 0
               !if (j >= ) then
                  wx = (h2osoi_ice(c,j) + h2osoi_liq(c,j))
                  void = 1._r8 - (h2osoi_ice(c,j)/denice + h2osoi_liq(c,j)/denh2o)&
                       /(frac_sno(c) * dz(c,j))
                  ! If void is negative, then increase dz such that void = 0.
                  ! This should be done for any landunit, but for now is done only for glacier_mec 1andunits.
                
                  ! I don't think the next 5 lines are necessary (removed in CLMv5)
                  l = col_pp%landunit(c)
                  if (ltype(l)==istice_mec .and. void < 0._r8) then
                     dz(c,j) = h2osoi_ice(c,j)/denice + h2osoi_liq(c,j)/denh2o
                     void = 0._r8
                  endif

                  ! Allow compaction only for non-saturated node and higher ice lens node.
                  if (void > 0.001_r8 .and. h2osoi_ice(c,j) > .1_r8) then

                      bi = h2osoi_ice(c,j) / (frac_sno(c) * dz(c,j))
                      fi = h2osoi_ice(c,j) / wx
                      td = tfrz-t_soisno(c,j)
                      dexpf = exp(-c4*td)

                      ! Settling as a result of destructive metamorphism
                      ddz1 = -c3*dexpf 
                      if (bi > dm) ddz1 = ddz1*exp(-46.0e-3_r8*(bi-dm))

                      ! Liquid water term

                      if (h2osoi_liq(c,j) > 0.01_r8*dz(c,j)*frac_sno(c)) ddz1=ddz1*c5

                      ! Compaction due to overburden
                     ddz2 = -(burden_noextra+wx/2._r8)*exp(-0.08_r8*td - c2*bi)/eta0 

                      ! Compaction occurring during melt

                      if (imelt(c,j) == 1) then
                         if(subgridflag==1 .and. (ltype(col_pp%landunit(c)) == istsoil .or. ltype(col_pp%landunit(c)) == istcrop)) then
                            ! first term is delta mass over mass
                            ddz3 = max(0._r8,min(1._r8,(swe_old(c,j) - wx)/wx))

                            ! 2nd term is delta fsno over fsno, allowing for negative values for ddz3
                            if ((swe_old(c,j) - wx) > 0._r8) then
                               wsum = sum(h2osoi_liq(c,snl(c)+1:0)+h2osoi_ice(c,snl(c)+1:0))
                               fsno_melt = 1. - (acos(2.*min(1._r8,wsum/int_snow(c)) - 1._r8)/rpi)**(n_melt(c))

                               ddz3 = ddz3 - max(0._r8,(fsno_melt - frac_sno(c))/frac_sno(c))
                            endif
                            ddz3 = -1._r8/dtime * ddz3
                         else
                            ddz3 = - 1._r8/dtime * max(0._r8,(frac_iceold(c,j) - fi)/frac_iceold(c,j))
                         endif
                      else
                         ddz3 = 0._r8
                      end if

                     ddz4 = 0.0_r8

                     ! Time rate of fractional change in dz (units of s-1)
                      pdzdtc = ddz1 + ddz2 + ddz3 + ddz4

                      ! The change in dz due to compaction
                      ! Limit compaction to be no greater than fully saturated layer thickness

                      dz(c,j) = max(dz(c,j) * (1._r8+pdzdtc*dtime),(h2osoi_ice(c,j)/denice+ h2osoi_liq(c,j)/denh2o)/frac_sno(c))
                  
                  ! NOTE: This seems irrelevant unless extra snow layers is true?
                  ! else ! from CLMv5
                  !    ! saturated node is immobile
                  !    !
                  !    ! This is only needed if wind_dependent_snow_density is true, but it's
                  !    ! simplest just to update mobile always
                  !    mobile(fc) = .false.

                  end if

                  ! Pressure of overlying snow
                  burden_noextra = burden_noextra + wx
               !end if 
            end do !j loop now 
         end do  ! filter loop 
      end if



     !$acc exit data delete(&
     !$acc burden(:), &
     !$acc zpseudo(:), &
     !$acc mobile(:), &
     !$acc ddz1_fresh, &
     !$acc ddz1, &
     !$acc ddz3, &
     !$acc burden_noextra)

     end associate
   end subroutine SnowCompaction

   !-----------------------------------------------------------------------
   subroutine CombineSnowLayers(bounds, num_snowc, filter_snowc, &
        aerosol_vars, dtime)
     !
     ! !DESCRIPTION:
     ! Combine snow layers that are less than a minimum thickness or mass
     ! If the snow element thickness or mass is less than a prescribed minimum,
     ! then it is combined with a neighboring element.  The subroutine
     ! clm\_combo.f90 then executes the combination of mass and energy.
     !
     ! !USES:
     use landunit_varcon  , only : istsoil, istdlak, istsoil, istwet, istice, istice_mec, istcrop
     use LakeCon          , only : lsadz
     use elm_varcon       , only : denh2o
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds
     integer                , intent(inout) :: num_snowc       ! number of column snow points in column filter
     integer                , intent(inout) :: filter_snowc(:) ! column filter for snow points
     type(aerosol_type)     , intent(inout) :: aerosol_vars
     real(r8), intent(in)  :: dtime
     !
     ! !LOCAL VARIABLES:
     integer :: c, fc                            ! column indices
     integer :: i,k                              ! loop indices
     integer :: j,l                              ! node indices
     integer :: msn_old(1:num_snowc) ! number of top snow layer
     integer :: mssi   (1:num_snowc)    ! node index
     integer :: neibor                           ! adjacent node selected for combination
     real(r8):: dzminloc_mssi_c                  ! dzminloc evaluated at mssi(c)
     real(r8):: zwice (1:num_snowc)   ! total ice mass in snow
     real(r8):: zwliq (1:num_snowc)  ! total liquid water in snow
     real(r8):: dzmin(5)                         ! minimum of top snow layer
     real(r8):: dzminloc(5)                ! minimum of top snow layer (local)
     real(r8):: dzminloc16(16)             ! minimum of top snow layer (local)
     real(r8) :: sum1, sum2, sum3, sum4 

     data dzmin /0.010_r8, 0.015_r8, 0.025_r8, 0.055_r8, 0.115_r8/

     !-----------------------------------------------------------------------

     associate(                                                     &
          ltype            => lun_pp%itype                           , & ! Input:  [integer  (:)   ] landunit type
          urbpoi           => lun_pp%urbpoi                          , & ! Input:  [logical  (:)   ] true => landunit is an urban point

          t_soisno         => col_es%t_soisno       , & ! Output: [real(r8) (:,:) ] soil temperature (Kelvin)

          mss_bcphi        => aerosol_vars%mss_bcphi_col          , & ! Output: [real(r8) (:,:) ] hydrophilic BC mass in snow (col,lyr) [kg]
          mss_bcpho        => aerosol_vars%mss_bcpho_col          , & ! Output: [real(r8) (:,:) ] hydrophobic BC mass in snow (col,lyr) [kg]
          mss_ocphi        => aerosol_vars%mss_ocphi_col          , & ! Output: [real(r8) (:,:) ] hydrophilic OC mass in snow (col,lyr) [kg]
          mss_ocpho        => aerosol_vars%mss_ocpho_col          , & ! Output: [real(r8) (:,:) ] hydrophobic OC mass in snow (col,lyr) [kg]
          mss_dst1         => aerosol_vars%mss_dst1_col           , & ! Output: [real(r8) (:,:) ] dust species 1 mass in snow (col,lyr) [kg]
          mss_dst2         => aerosol_vars%mss_dst2_col           , & ! Output: [real(r8) (:,:) ] dust species 2 mass in snow (col,lyr) [kg]
          mss_dst3         => aerosol_vars%mss_dst3_col           , & ! Output: [real(r8) (:,:) ] dust species 3 mass in snow (col,lyr) [kg]
          mss_dst4         => aerosol_vars%mss_dst4_col           , & ! Output: [real(r8) (:,:) ] dust species 4 mass in snow (col,lyr) [kg]

          frac_sno         => col_ws%frac_sno        , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
          frac_sno_eff     => col_ws%frac_sno_eff    , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
          snow_depth       => col_ws%snow_depth      , & ! Input:  [real(r8) (:)   ] snow height (m)
          int_snow         => col_ws%int_snow        , & ! Input:  [real(r8) (:)   ] integrated snowfall [mm]
          h2osno           => col_ws%h2osno          , & ! Output: [real(r8) (:)   ] snow water (mm H2O)
          h2osoi_ice       => col_ws%h2osoi_ice      , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)
          h2osoi_liq       => col_ws%h2osoi_liq      , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
          snw_rds          => col_ws%snw_rds         , & ! Output: [real(r8) (:,:) ] effective snow grain radius (col,lyr) [microns, m^-6]
          mflx_snowlyr_col => col_wf%mflx_snowlyr    , & ! Output: [real(r8) (:)   ]  mass flux to top soil layer due to disappearance of snow (kg H2O /s)

          qflx_sl_top_soil => col_wf%qflx_sl_top_soil , & ! Output: [real(r8) (:)   ] liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
          qflx_snow2topsoi => col_wf%qflx_snow2topsoi , & ! Output: [real(r8) (:)   ] liquid water merged into top soil from snow

          snl              => col_pp%snl                             , & ! Output: [integer  (:)   ] number of snow layers
          dz               => col_pp%dz                              , & ! Output: [real(r8) (:,:) ] layer depth (m)
          zi               => col_pp%zi                              , & ! Output: [real(r8) (:,:) ] interface level below a "z" level (m)
          z                => col_pp%z                                 & ! Output: [real(r8) (:,:) ] layer thickness (m)
          )
     !$acc enter data create(&
     !$acc msn_old(:) , &
     !$acc mssi(:)    , &
     !$acc zwice(:)   , &
     !$acc zwliq(:)   , &
     !$acc dzmin(:)   , &
     !$acc dzminloc(:), &
     !$acc dzminloc16(:), &
     !$acc sum1, sum2, sum3, sum4)

       ! Check the mass of ice lens of snow, when the total is less than a small value,
       ! combine it with the underlying neighbor.
       

       ! dzmin will stay constant between timesteps
       if (.not. use_extrasnowlayers) then
         dzminloc(:) = dzmin(:)
       else
         dzminloc16(:) = dzmin16(:)
       endif

       ! Add lsadz to dzmin for lakes
       ! Determine whether called from LakeHydrology
       ! Note: this assumes that this function is called separately with the lake-snow and non-lake-snow filters.
       if (num_snowc > 0) then
          c = filter_snowc(1)
          l = col_pp%landunit(c)
          if (ltype(l) == istdlak) then ! Called from LakeHydrology
            if (.not. use_extrasnowlayers) then
               dzminloc(:) = dzmin(:) + lsadz
            else
               dzminloc16(:) = dzmin16(:) + lsadz
            end if
          end if
       end if

       !$acc parallel loop independent gang vector default(present)
       do fc = 1, num_snowc
          c = filter_snowc(fc)

          msn_old(fc) = snl(c)
          qflx_sl_top_soil(c) = 0._r8
          qflx_snow2topsoi(c) = 0._r8
          mflx_snowlyr_col(c) = 0._r8
       end do

       ! The following loop is NOT VECTORIZED
       !$acc parallel loop independent gang vector default(present)
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          l = col_pp%landunit(c)
          !$acc loop seq
          do j = msn_old(fc)+1,0
             ! use 0.01 to avoid runaway ice buildup
             if (h2osoi_ice(c,j) <= .01_r8) then
                if (ltype(l) == istsoil .or. urbpoi(l) .or. ltype(l) == istcrop) then
                   h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + h2osoi_liq(c,j)
                   h2osoi_ice(c,j+1) = h2osoi_ice(c,j+1) + h2osoi_ice(c,j)

                   if (j == 0) then
                      qflx_sl_top_soil(c) = (h2osoi_liq(c,j) + h2osoi_ice(c,j))/dtime
                      mflx_snowlyr_col(c) = mflx_snowlyr_col(c) + qflx_sl_top_soil(c)
                   end if

                   if (j /= 0) dz(c,j+1) = dz(c,j+1) + dz(c,j)

                   ! NOTE: Temperature, and similarly snw_rds, of the
                   ! underlying snow layer are NOT adjusted in this case.
                   ! Because the layer being eliminated has a small mass,
                   ! this should not make a large difference, but it
                   ! would be more thorough to do so.
                   if (j /= 0) then
                      mss_bcphi(c,j+1) = mss_bcphi(c,j+1)  + mss_bcphi(c,j)
                      mss_bcpho(c,j+1) = mss_bcpho(c,j+1)  + mss_bcpho(c,j)
                      mss_ocphi(c,j+1) = mss_ocphi(c,j+1)  + mss_ocphi(c,j)
                      mss_ocpho(c,j+1) = mss_ocpho(c,j+1)  + mss_ocpho(c,j)
                      mss_dst1(c,j+1)  = mss_dst1(c,j+1)   + mss_dst1(c,j)
                      mss_dst2(c,j+1)  = mss_dst2(c,j+1)   + mss_dst2(c,j)
                      mss_dst3(c,j+1)  = mss_dst3(c,j+1)   + mss_dst3(c,j)
                      mss_dst4(c,j+1)  = mss_dst4(c,j+1)   + mss_dst4(c,j)
                   end if

                else if (ltype(l) /= istsoil .and. .not. urbpoi(l) .and. ltype(l) /= istcrop .and. j /= 0) then

                   h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + h2osoi_liq(c,j)
                   h2osoi_ice(c,j+1) = h2osoi_ice(c,j+1) + h2osoi_ice(c,j)
                   dz(c,j+1) = dz(c,j+1) + dz(c,j)

                   mss_bcphi(c,j+1) = mss_bcphi(c,j+1)  + mss_bcphi(c,j)
                   mss_bcpho(c,j+1) = mss_bcpho(c,j+1)  + mss_bcpho(c,j)
                   mss_ocphi(c,j+1) = mss_ocphi(c,j+1)  + mss_ocphi(c,j)
                   mss_ocpho(c,j+1) = mss_ocpho(c,j+1)  + mss_ocpho(c,j)
                   mss_dst1(c,j+1)  = mss_dst1(c,j+1)   + mss_dst1(c,j)
                   mss_dst2(c,j+1)  = mss_dst2(c,j+1)   + mss_dst2(c,j)
                   mss_dst3(c,j+1)  = mss_dst3(c,j+1)   + mss_dst3(c,j)
                   mss_dst4(c,j+1)  = mss_dst4(c,j+1)   + mss_dst4(c,j)

                end if

                ! shift all elements above this down one.
                if (j > snl(c)+1 .and. snl(c) < -1) then
                  !$acc loop seq 
                   do i = j, snl(c)+2, -1
                      ! If the layer closest to the surface is less than 0.1 mm and the ltype is not
                      ! urban, soil or crop, the h2osoi_liq and h2osoi_ice associated with this layer is sent
                      ! to qflx_qrgwl later on in the code.  To keep track of this for the snow balance
                      ! error check, we add this to qflx_sl_top_soil here
                      if (ltype(l) /= istsoil .and. ltype(l) /= istcrop .and. .not. urbpoi(l) .and. i == 0) then
                         qflx_sl_top_soil(c) = (h2osoi_liq(c,i) + h2osoi_ice(c,i))/dtime
                      end if

                      t_soisno(c,i)   = t_soisno(c,i-1)
                      h2osoi_liq(c,i) = h2osoi_liq(c,i-1)
                      h2osoi_ice(c,i) = h2osoi_ice(c,i-1)

                      mss_bcphi(c,i)   = mss_bcphi(c,i-1)
                      mss_bcpho(c,i)   = mss_bcpho(c,i-1)
                      mss_ocphi(c,i)   = mss_ocphi(c,i-1)
                      mss_ocpho(c,i)   = mss_ocpho(c,i-1)
                      mss_dst1(c,i)    = mss_dst1(c,i-1)
                      mss_dst2(c,i)    = mss_dst2(c,i-1)
                      mss_dst3(c,i)    = mss_dst3(c,i-1)
                      mss_dst4(c,i)    = mss_dst4(c,i-1)
                      snw_rds(c,i)     = snw_rds(c,i-1)
                      dz(c,i)         = dz(c,i-1)
                   end do
                end if
                snl(c) = snl(c) + 1
             end if
          end do
       end do

       !$acc parallel loop independent gang worker default(present) private(sum1,sum2,sum3,sum4)
       do fc = 1, num_snowc
         c = filter_snowc(fc)
         sum1 = 0._r8; sum2 = 0._r8; 
         sum3 = 0._r8; sum4 = 0._r8;
         !$acc loop vector reduction(+:sum1,sum2,sum3,sum4)
         do j = -nlevsno+1,0
             if (j >= snl(c)+1) then
                sum1 = sum1 + h2osoi_ice(c,j) + h2osoi_liq(c,j)
                sum2 = sum2 + dz(c,j)
                sum3 = sum3 + h2osoi_ice(c,j)
                sum4 = sum4 + h2osoi_liq(c,j)
             end if
         end do
         h2osno(c) = sum1 
         snow_depth(c) = sum2
         zwice(fc)  = sum3 
         zwliq(fc) = sum4  
       end do

       ! Check the snow depth - all snow gone
       ! The liquid water assumes ponding on soil surface.

       !$acc parallel loop independent gang vector default(present)
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          l = col_pp%landunit(c)
          if (snow_depth(c) > 0._r8) then
             if ((ltype(l) == istdlak .and. snow_depth(c) < 0.01_r8 + lsadz ) .or. &
                  ((ltype(l) /= istdlak) .and. ((frac_sno_eff(c)*snow_depth(c) < 0.01_r8)  &
                  .or. (h2osno(c)/(frac_sno_eff(c)*snow_depth(c)) < 50._r8)))) then

                snl(c) = 0
                h2osno(c) = zwice(fc)

                mss_bcphi(c,:) = 0._r8
                mss_bcpho(c,:) = 0._r8
                mss_ocphi(c,:) = 0._r8
                mss_ocpho(c,:) = 0._r8
                mss_dst1(c,:)  = 0._r8
                mss_dst2(c,:)  = 0._r8
                mss_dst3(c,:)  = 0._r8
                mss_dst4(c,:)  = 0._r8

                if (h2osno(c) <= 0._r8) snow_depth(c) = 0._r8
                ! this is where water is transfered from layer 0 (snow) to layer 1 (soil)
                if (ltype(l) == istsoil .or. urbpoi(l) .or. ltype(l) == istcrop) then
                   h2osoi_liq(c,0) = 0.0_r8
                   h2osoi_liq(c,1) = h2osoi_liq(c,1) + zwliq(fc)
                   qflx_snow2topsoi(c) = zwliq(fc)/dtime
                   mflx_snowlyr_col(c) = mflx_snowlyr_col(c) + zwliq(fc)/dtime
                end if
                if (ltype(l) == istwet) then
                   h2osoi_liq(c,0) = 0.0_r8
                endif
                if (ltype(l) == istice .or. ltype(l)==istice_mec) then
                   h2osoi_liq(c,0) = 0.0_r8
                endif
             endif
          end if
          if (h2osno(c) <= 0._r8) then
             snow_depth(c) = 0._r8
             frac_sno(c) = 0._r8
             frac_sno_eff(c) = 0._r8
             int_snow(c) = 0._r8
          endif
       end do

       ! Check the snow depth - snow layers combined
       ! The following loop IS NOT VECTORIZED

       do fc = 1, num_snowc
          c = filter_snowc(fc)

          ! Two or more layers

          if (snl(c) < -1) then

             msn_old(fc) = snl(c)
             mssi(fc) = 1

             do i = msn_old(fc)+1,0
                if (.not. use_extrasnowlayers) then
                    dzminloc_mssi_c = dzminloc(mssi(fc))
                else
                    dzminloc_mssi_c = dzminloc16(mssi(fc))
                end if
                if ((frac_sno_eff(c)*dz(c,i) < dzminloc_mssi_c) .or. &
                     ((h2osoi_ice(c,i) + h2osoi_liq(c,i))/(frac_sno_eff(c)*dz(c,i)) < 50._r8)) then
                   if (i == snl(c)+1) then
                      ! If top node is removed, combine with bottom neighbor.
                      neibor = i + 1
                   else if (i == 0) then
                      ! If the bottom neighbor is not snow, combine with the top neighbor.
                      neibor = i - 1
                   else
                      ! If none of the above special cases apply, combine with the thinnest neighbor
                      neibor = i + 1
                      if ((dz(c,i-1)+dz(c,i)) < (dz(c,i+1)+dz(c,i))) neibor = i-1
                   end if

                   ! Node l and j are combined and stored as node j.
                   if (neibor > i) then
                      j = neibor
                      l = i
                   else
                      j = i
                      l = neibor
                   end if

                   ! this should be included in 'Combo' for consistency,
                   ! but functionally it is the same to do it here
                   mss_bcphi(c,j)=mss_bcphi(c,j)+mss_bcphi(c,l)
                   mss_bcpho(c,j)=mss_bcpho(c,j)+mss_bcpho(c,l)
                   mss_ocphi(c,j)=mss_ocphi(c,j)+mss_ocphi(c,l)
                   mss_ocpho(c,j)=mss_ocpho(c,j)+mss_ocpho(c,l)
                   mss_dst1(c,j)=mss_dst1(c,j)+mss_dst1(c,l)
                   mss_dst2(c,j)=mss_dst2(c,j)+mss_dst2(c,l)
                   mss_dst3(c,j)=mss_dst3(c,j)+mss_dst3(c,l)
                   mss_dst4(c,j)=mss_dst4(c,j)+mss_dst4(c,l)

                   ! mass-weighted combination of effective grain size:

                   snw_rds(c,j) = (snw_rds(c,j)*(h2osoi_liq(c,j)+h2osoi_ice(c,j)) + &
                        snw_rds(c,l)*(h2osoi_liq(c,l)+h2osoi_ice(c,l))) / &
                        (h2osoi_liq(c,j)+h2osoi_ice(c,j)+h2osoi_liq(c,l)+h2osoi_ice(c,l))

                   call Combo (dz(c,j), h2osoi_liq(c,j), h2osoi_ice(c,j), &
                        t_soisno(c,j), dz(c,l), h2osoi_liq(c,l), h2osoi_ice(c,l), t_soisno(c,l) )

                   ! Now shift all elements above this down one.
                   if (j-1 > snl(c)+1) then

                      do k = j-1, snl(c)+2, -1
                         t_soisno(c,k) = t_soisno(c,k-1)
                         h2osoi_ice(c,k) = h2osoi_ice(c,k-1)
                         h2osoi_liq(c,k) = h2osoi_liq(c,k-1)

                         mss_bcphi(c,k) = mss_bcphi(c,k-1)
                         mss_bcpho(c,k) = mss_bcpho(c,k-1)
                         mss_ocphi(c,k) = mss_ocphi(c,k-1)
                         mss_ocpho(c,k) = mss_ocpho(c,k-1)
                         mss_dst1(c,k)  = mss_dst1(c,k-1)
                         mss_dst2(c,k)  = mss_dst2(c,k-1)
                         mss_dst3(c,k)  = mss_dst3(c,k-1)
                         mss_dst4(c,k)  = mss_dst4(c,k-1)
                         snw_rds(c,k)   = snw_rds(c,k-1)

                         dz(c,k) = dz(c,k-1)
                      end do
                   end if

                   ! Decrease the number of snow layers
                   snl(c) = snl(c) + 1
                   if (snl(c) >= -1) EXIT

                else

                   ! The layer thickness is greater than the prescribed minimum value
                   mssi(fc) = mssi(fc) + 1

                end if
             end do

          end if

       end do

       ! Reset the node depth and the depth of layer interface
       !NOTE: racecondition? Testing out two loop structure:
       !$acc parallel default(present) 
      !$acc loop independent gang vector 
      do fc = 1, num_snowc
         c = filter_snowc(fc)
         !NOTE: can replace the loop bounds with [0,snl(c)+1] here now? 
         !$acc loop seq 
         do j = 0, -nlevsno+1, -1 
            if(j >= snl(c) + 1) then 
               zi(c,j-1) = zi(c,j) - dz(c,j)
            end if 
         end do 
      end do 
             
      !$acc loop independent gang vector collapse(2) 
      do j = 0, -nlevsno+1, -1
         do fc = 1, num_snowc
            c = filter_snowc(fc)
            if (j >= snl(c) + 1) then
               z(c,j) = zi(c,j) - 0.5_r8*dz(c,j)
            end if
         end do
      end do
       !$acc end parallel 
    


     !$acc exit data delete(&
     !$acc msn_old(:), &
     !$acc mssi(:), &
     !$acc zwice(:), &
     !$acc zwliq(:), &
     !$acc dzmin(:), &
     !$acc dzminloc(:), &
     !$acc dzminloc16(:), &
     !$acc sum1, &
     !$acc sum2, &
     !$acc sum3, &
     !$acc sum4)

    end associate

   end subroutine CombineSnowLayers

   !-----------------------------------------------------------------------
   subroutine DivideSnowLayers(bounds, num_snowc, filter_snowc, &
        aerosol_vars, is_lake)
     !
     ! !DESCRIPTION:
     ! Subdivides snow layers if they exceed their prescribed maximum thickness.
     !
     ! !USES:
     use elm_varcon,  only : tfrz
     use LakeCon   ,  only : lsadz
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds
     integer                , intent(in)    :: num_snowc       ! number of column snow points in column filter
     integer                , intent(in)    :: filter_snowc(:) ! column filter for snow points
     type(aerosol_type)     , intent(inout) :: aerosol_vars
     logical                , intent(in)    :: is_lake  !TODO - this should be examined and removed in the future
     !
     ! !LOCAL VARIABLES:
     integer  :: j, c, fc, k                  ! indices
     real(r8) :: drr                          ! thickness of the combined [m]
     integer  :: msno                         ! number of snow layer 1 (top) to msno (bottom)
     real(r8) :: dzsno(1:num_snowc,nlevsno)   ! Snow layer thickness [m]
     real(r8) :: swice(1:num_snowc,nlevsno)   ! Partial volume of ice [m3/m3]
     real(r8) :: swliq(1:num_snowc,nlevsno)   ! Partial volume of liquid water [m3/m3]
     real(r8) :: tsno(1:num_snowc ,nlevsno)   ! Nodel temperature [K]
     real(r8) :: zwice                        ! temporary
     real(r8) :: zwliq                        ! temporary
     real(r8) :: propor                       ! temporary
     real(r8) :: dtdz                         ! temporary
     ! temporary variables mimicking the structure of other layer division variables
     real(r8) :: mbc_phi(1:num_snowc,nlevsno) ! mass of BC in each snow layer
     real(r8) :: zmbc_phi                     ! temporary
     real(r8) :: mbc_pho(1:num_snowc,nlevsno) ! mass of BC in each snow layer
     real(r8) :: zmbc_pho                     ! temporary
     real(r8) :: moc_phi(1:num_snowc,nlevsno) ! mass of OC in each snow layer
     real(r8) :: zmoc_phi                     ! temporary
     real(r8) :: moc_pho(1:num_snowc,nlevsno) ! mass of OC in each snow layer
     real(r8) :: zmoc_pho                     ! temporary
     real(r8) :: mdst1(1:num_snowc,nlevsno)   ! mass of dust 1 in each snow layer
     real(r8) :: zmdst1                       ! temporary
     real(r8) :: mdst2(1:num_snowc,nlevsno)   ! mass of dust 2 in each snow layer
     real(r8) :: zmdst2                       ! temporary
     real(r8) :: mdst3(1:num_snowc,nlevsno)   ! mass of dust 3 in each snow layer
     real(r8) :: zmdst3                       ! temporary
     real(r8) :: mdst4(1:num_snowc,nlevsno)   ! mass of dust 4 in each snow layer
     real(r8) :: zmdst4                       ! temporary
     real(r8) :: rds(1:num_snowc,nlevsno)
     ! Variables for consistency check
     real(r8) :: dztot(1:num_snowc)
     real(r8) :: snwicetot(1:num_snowc)
     real(r8) :: snwliqtot(1:num_snowc)
     real(r8) :: offset ! temporary
     real(r8) :: sum1, sum2, sum3 
     integer :: snl_idx
     !-----------------------------------------------------------------------

     associate(                                            &
          t_soisno   => col_es%t_soisno    , & ! Output: [real(r8) (:,:) ] soil temperature (Kelvin)

          h2osoi_ice => col_ws%h2osoi_ice   , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)
          h2osoi_liq => col_ws%h2osoi_liq   , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
          frac_sno   => col_ws%frac_sno_eff , & ! Output: [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
          snw_rds    => col_ws%snw_rds      , & ! Output: [real(r8) (:,:) ] effective snow grain radius (col,lyr) [microns, m^-6]

          mss_bcphi  => aerosol_vars%mss_bcphi_col  , & ! Output: [real(r8) (:,:) ] hydrophilic BC mass in snow (col,lyr) [kg]
          mss_bcpho  => aerosol_vars%mss_bcpho_col  , & ! Output: [real(r8) (:,:) ] hydrophobic BC mass in snow (col,lyr) [kg]
          mss_ocphi  => aerosol_vars%mss_ocphi_col  , & ! Output: [real(r8) (:,:) ] hydrophilic OC mass in snow (col,lyr) [kg]
          mss_ocpho  => aerosol_vars%mss_ocpho_col  , & ! Output: [real(r8) (:,:) ] hydrophobic OC mass in snow (col,lyr) [kg]
          mss_dst1   => aerosol_vars%mss_dst1_col   , & ! Output: [real(r8) (:,:) ] dust species 1 mass in snow (col,lyr) [kg]
          mss_dst2   => aerosol_vars%mss_dst2_col   , & ! Output: [real(r8) (:,:) ] dust species 2 mass in snow (col,lyr) [kg]
          mss_dst3   => aerosol_vars%mss_dst3_col   , & ! Output: [real(r8) (:,:) ] dust species 3 mass in snow (col,lyr) [kg]
          mss_dst4   => aerosol_vars%mss_dst4_col   , & ! Output: [real(r8) (:,:) ] dust species 4 mass in snow (col,lyr) [kg]

          snl        => col_pp%snl                  , & ! Output: [integer  (:)   ] number of snow layers
          dz         => col_pp%dz                   , & ! Output: [real(r8) (:,:) ] layer depth (m)
          zi         => col_pp%zi                   , & ! Output: [real(r8) (:,:) ] interface level below a "z" level (m)
          z          => col_pp%z                      & ! Output: [real(r8) (:,:) ] layer thickness (m)
          )
     !$acc enter data create(&
     !$acc dzsno(:,:), &
     !$acc swice(:,:), &
     !$acc swliq(:,:), &
     !$acc tsno(:,:), &
     !$acc mbc_phi(:,:), &
     !$acc mbc_pho(:,:), &
     !$acc moc_phi(:,:), &
     !$acc moc_pho(:,:), &
     !$acc mdst1(:,:), &
     !$acc mdst2(:,:), &
     !$acc mdst3(:,:), &
     !$acc mdst4(:,:), &
     !$acc rds(:,:), &
     !$acc dztot(:), &
     !$acc snwicetot(:), &
     !$acc snwliqtot(:), &
     !$acc sum1, &
     !$acc sum2, &
     !$acc sum3)

       if ( is_lake ) then
          ! Initialize for consistency check
         !$acc parallel loop independent gang vector default(present) private(sum1,sum2,sum3)
         do fc = 1, num_snowc
            c = filter_snowc(fc)
            sum1 = 0._r8; sum2 = 0._r8; sum3 =0._r8  
            !$acc loop vector reduction(+:sum1,sum2,sum3)
            do j = -nlevsno+1,0
                if (j >= snl(c)+1) then
                   sum1 = sum1 + dz(c,j)
                   sum2 = sum2 + h2osoi_ice(c,j)
                   sum3 = sum3 + h2osoi_liq(c,j)
                end if
             end do
             dztot(fc)    = sum1  
            snwicetot(fc) = sum2 
            snwliqtot(fc) = sum3
          end do
       end if

       ! Begin calculation - note that the following column loops are only invoked
       ! for snow-covered columns

       !$acc parallel loop independent gang vector default(present) collapse(2) 
       !!!$acc present(mss_bcphi(:,:),mss_bcpho(:,:),mss_ocphi(:,:),mss_ocpho(:,:) &
       !!!$acc ,mss_dst1(:,:),mss_dst2(:,:),mss_dst3(:,:),mss_dst4(:,:) ,snw_rds(:,:), &
       !!!$acc  h2osoi_ice(:,:),h2osoi_liq(:,:), t_soisno(:,:),dz(:,:),snl(:))
       do j = 1,nlevsno
          do fc = 1, num_snowc
             c = filter_snowc(fc)
             snl_idx = j+snl(c)
             if (j <= abs(snl(c))) then
                if (is_lake) then
                   dzsno(fc,j) = dz(c,snl_idx)
                else
                   dzsno(fc,j) = frac_sno(c)*dz(c,snl_idx)
                end if
                swice(fc,j) = h2osoi_ice(c,snl_idx)
                swliq(fc,j) = h2osoi_liq(c,snl_idx)
                tsno (fc,j) = t_soisno(c,snl_idx)

                mbc_phi(fc,j) = mss_bcphi(c,snl_idx)
                mbc_pho(fc,j) = mss_bcpho(c,snl_idx)
                moc_phi(fc,j) = mss_ocphi(c,snl_idx)
                moc_pho(fc,j) = mss_ocpho(c,snl_idx)
                mdst1(fc,j)   = mss_dst1(c,snl_idx)
                mdst2(fc,j)   = mss_dst2(c,snl_idx)
                mdst3(fc,j)   = mss_dst3(c,snl_idx)
                mdst4(fc,j)   = mss_dst4(c,snl_idx)
                rds(fc,j)     = snw_rds(c,snl_idx)
             end if
          end do
       end do

       !$acc parallel loop independent gang vector default(present)
       do fc = 1, num_snowc
           c = filter_snowc(fc)

           msno = abs(snl(c))

           if (msno == 1) then
              ! Specify a new snow layer
              if (is_lake) then
                 offset = 2._r8 * lsadz
              else
                 offset = 0._r8
              end if
              if (dzsno(fc,1) > 0.03_r8 + offset) then
                 msno = 2
                 dzsno(fc,1) = dzsno(fc,1)/2._r8
                 swice(fc,1) = swice(fc,1)/2._r8
                 swliq(fc,1) = swliq(fc,1)/2._r8
                 dzsno(fc,2) = dzsno(fc,1)
                 swice(fc,2) = swice(fc,1)
                 swliq(fc,2) = swliq(fc,1)
                 tsno(fc,2)  = tsno(fc,1)

                 mbc_phi(fc,1) = mbc_phi(fc,1)/2._r8
                 mbc_phi(fc,2) = mbc_phi(fc,1)
                 mbc_pho(fc,1) = mbc_pho(fc,1)/2._r8
                 mbc_pho(fc,2) = mbc_pho(fc,1)
                 moc_phi(fc,1) = moc_phi(fc,1)/2._r8
                 moc_phi(fc,2) = moc_phi(fc,1)
                 moc_pho(fc,1) = moc_pho(fc,1)/2._r8
                 moc_pho(fc,2) = moc_pho(fc,1)
                 mdst1(fc,1) = mdst1(fc,1)/2._r8
                 mdst1(fc,2) = mdst1(fc,1)
                 mdst2(fc,1) = mdst2(fc,1)/2._r8
                 mdst2(fc,2) = mdst2(fc,1)
                 mdst3(fc,1) = mdst3(fc,1)/2._r8
                 mdst3(fc,2) = mdst3(fc,1)
                 mdst4(fc,1) = mdst4(fc,1)/2._r8
                 mdst4(fc,2) = mdst4(fc,1)
                 rds(fc,2) = rds(fc,1)

              end if
           end if

           if (msno > 1) then
              if (is_lake) then
                 offset = lsadz
              else
                 offset = 0._r8
              end if
              if (dzsno(fc,1) > 0.02_r8 + offset) then
                 if (is_lake) then
                    drr = dzsno(fc,1) - 0.02_r8 - lsadz
                 else
                    drr = dzsno(fc,1) - 0.02_r8
                 end if
                 propor = drr/dzsno(fc,1)
                 zwice = propor*swice(fc,1)
                 zwliq = propor*swliq(fc,1)

                 zmbc_phi = propor*mbc_phi(fc,1)
                 zmbc_pho = propor*mbc_pho(fc,1)
                 zmoc_phi = propor*moc_phi(fc,1)
                 zmoc_pho = propor*moc_pho(fc,1)
                 zmdst1 = propor*mdst1(fc,1)
                 zmdst2 = propor*mdst2(fc,1)
                 zmdst3 = propor*mdst3(fc,1)
                 zmdst4 = propor*mdst4(fc,1)

                 if (is_lake) then
                    propor = (0.02_r8+lsadz)/dzsno(fc,1)
                 else
                    propor = 0.02_r8/dzsno(fc,1)
                 endif 

                 swice(fc,1) = propor*swice(fc,1)
                 swliq(fc,1) = propor*swliq(fc,1)

                 mbc_phi(fc,1) = propor*mbc_phi(fc,1)
                 mbc_pho(fc,1) = propor*mbc_pho(fc,1)
                 moc_phi(fc,1) = propor*moc_phi(fc,1)
                 moc_pho(fc,1) = propor*moc_pho(fc,1)
                 mdst1(fc,1) = propor*mdst1(fc,1)
                 mdst2(fc,1) = propor*mdst2(fc,1)
                 mdst3(fc,1) = propor*mdst3(fc,1)
                 mdst4(fc,1) = propor*mdst4(fc,1)

                 if (is_lake) then
                    dzsno(fc,1) = 0.02_r8 + lsadz
                 else
                    dzsno(fc,1) = 0.02_r8
                 end if

                 mbc_phi(fc,2) = mbc_phi(fc,2)+zmbc_phi  ! (combo)
                 mbc_pho(fc,2) = mbc_pho(fc,2)+zmbc_pho  ! (combo)
                 moc_phi(fc,2) = moc_phi(fc,2)+zmoc_phi  ! (combo)
                 moc_pho(fc,2) = moc_pho(fc,2)+zmoc_pho  ! (combo)
                 mdst1(fc,2) = mdst1(fc,2)+zmdst1  ! (combo)
                 mdst2(fc,2) = mdst2(fc,2)+zmdst2  ! (combo)
                 mdst3(fc,2) = mdst3(fc,2)+zmdst3  ! (combo)
                 mdst4(fc,2) = mdst4(fc,2)+zmdst4  ! (combo)
#ifdef MODAL_AER
              !mgf++ bugfix
              rds(fc,2) = (rds(fc,2)*(swliq(fc,2)+swice(fc,2)) + rds(fc,1)*(zwliq+zwice))/(swliq(fc,2)+swice(fc,2)+zwliq+zwice)
                if ((rds(fc,2) < 30.) .or. (rds(fc,2) > 1500.)) then
                   write (iulog,*) "2. SNICAR ERROR: snow grain radius of",rds(fc,2),rds(fc,1)
                   write (iulog,*) "swliq, swice, zwliq, zwice", swliq(fc,2), swice(fc,2),zwliq, zwice
                   write (iulog,*) "propor: ", propor
                   write (iulog,*) "dzsno : ", dzsno(fc,1)
                   write (iulog,*) "layers ", msno
                endif
              !mgf--
#else
              rds(fc,2) = rds(fc,1) ! (combo)
#endif

                 call Combo (dzsno(fc,2), swliq(fc,2), swice(fc,2), tsno(fc,2), drr, &
                      zwliq, zwice, tsno(fc,1))

                 ! Subdivide a new layer
                 if (is_lake) then
                    offset = 2._r8*lsadz
                 else
                    offset = 0._r8
                 end if
                 if (msno <= 2 .and. dzsno(fc,2) > 0.07_r8 + offset) then
                    msno = 3
                    dtdz = (tsno(fc,1) - tsno(fc,2))/((dzsno(fc,1)+dzsno(fc,2))/2._r8) 
                    dzsno(fc,2) = dzsno(fc,2)/2._r8
                    swice(fc,2) = swice(fc,2)/2._r8
                    swliq(fc,2) = swliq(fc,2)/2._r8
                    dzsno(fc,3) = dzsno(fc,2)
                    swice(fc,3) = swice(fc,2)
                    swliq(fc,3) = swliq(fc,2)
                    tsno(fc,3) = tsno(fc,2) - dtdz*dzsno(fc,2)/2._r8
                    if (tsno(fc,3) >= tfrz) then 
                       tsno(fc,3)  = tsno(fc,2)
                    else
                       tsno(fc,2) = tsno(fc,2) + dtdz*dzsno(fc,2)/2._r8 
                    endif

                    mbc_phi(fc,2) = mbc_phi(fc,2)/2._r8
                    mbc_phi(fc,3) = mbc_phi(fc,2)
                    mbc_pho(fc,2) = mbc_pho(fc,2)/2._r8
                    mbc_pho(fc,3) = mbc_pho(fc,2)
                    moc_phi(fc,2) = moc_phi(fc,2)/2._r8
                    moc_phi(fc,3) = moc_phi(fc,2)
                    moc_pho(fc,2) = moc_pho(fc,2)/2._r8
                    moc_pho(fc,3) = moc_pho(fc,2)
                    mdst1(fc,2) = mdst1(fc,2)/2._r8
                    mdst1(fc,3) = mdst1(fc,2)
                    mdst2(fc,2) = mdst2(fc,2)/2._r8
                    mdst2(fc,3) = mdst2(fc,2)
                    mdst3(fc,2) = mdst3(fc,2)/2._r8
                    mdst3(fc,3) = mdst3(fc,2)
                    mdst4(fc,2) = mdst4(fc,2)/2._r8
                    mdst4(fc,3) = mdst4(fc,2)
                    rds(fc,3) = rds(fc,2)

                 end if
              end if
           end if

           if (msno > 2) then
              if (is_lake) then
                 offset = lsadz
              else
                 offset = 0._r8
              end if
              if (dzsno(fc,2) > 0.05_r8+offset) then
                 if (is_lake) then
                    drr = dzsno(fc,2) - 0.05_r8 - lsadz
                 else
                    drr = dzsno(fc,2) - 0.05_r8
                 end if
                 propor = drr/dzsno(fc,2)
                 zwice = propor*swice(fc,2)
                 zwliq = propor*swliq(fc,2)

                 zmbc_phi = propor*mbc_phi(fc,2)
                 zmbc_pho = propor*mbc_pho(fc,2)
                 zmoc_phi = propor*moc_phi(fc,2)
                 zmoc_pho = propor*moc_pho(fc,2)
                 zmdst1 = propor*mdst1(fc,2)
                 zmdst2 = propor*mdst2(fc,2)
                 zmdst3 = propor*mdst3(fc,2)
                 zmdst4 = propor*mdst4(fc,2)

                 if (is_lake) then
                    propor = (0.05_r8+lsadz)/dzsno(fc,2)
                 else
                    propor = 0.05_r8/dzsno(fc,2)
                 end if
                 swice(fc,2) = propor*swice(fc,2)
                 swliq(fc,2) = propor*swliq(fc,2)

                 mbc_phi(fc,2) = propor*mbc_phi(fc,2)
                 mbc_pho(fc,2) = propor*mbc_pho(fc,2)
                 moc_phi(fc,2) = propor*moc_phi(fc,2)
                 moc_pho(fc,2) = propor*moc_pho(fc,2)
                 mdst1(fc,2) = propor*mdst1(fc,2)
                 mdst2(fc,2) = propor*mdst2(fc,2)
                 mdst3(fc,2) = propor*mdst3(fc,2)
                 mdst4(fc,2) = propor*mdst4(fc,2)

                 if (is_lake) then
                    dzsno(fc,2) = 0.05_r8+lsadz
                 else
                    dzsno(fc,2) = 0.05_r8
                 end if

                 mbc_phi(fc,3) = mbc_phi(fc,3)+zmbc_phi  ! (combo)
                 mbc_pho(fc,3) = mbc_pho(fc,3)+zmbc_pho  ! (combo)
                 moc_phi(fc,3) = moc_phi(fc,3)+zmoc_phi  ! (combo)
                 moc_pho(fc,3) = moc_pho(fc,3)+zmoc_pho  ! (combo)
                 mdst1(fc,3) = mdst1(fc,3)+zmdst1  ! (combo)
                 mdst2(fc,3) = mdst2(fc,3)+zmdst2  ! (combo)
                 mdst3(fc,3) = mdst3(fc,3)+zmdst3  ! (combo)
                 mdst4(fc,3) = mdst4(fc,3)+zmdst4  ! (combo)
#ifdef MODAL_AER
              !mgf++ bugfix
              rds(fc,3) = (rds(fc,3)*(swliq(fc,3)+swice(fc,3)) + rds(fc,2)*(zwliq+zwice))/(swliq(fc,3)+swice(fc,3)+zwliq+zwice)
                if ((rds(fc,3) < 30.) .or. (rds(fc,3) > 1500.)) then
#ifndef _OPENACC
                   write (iulog,*) "3. SNICAR ERROR: snow grain radius of",rds(fc,3),rds(fc,2)
                   write (iulog,*) "swliq, swice, zwliq, zwice", swliq(fc,3), swice(fc,3),zwliq, zwice
                   write (iulog,*) "layers ", msno
#endif
                endif
              !mgf--
#else
              rds(fc,3) = rds(fc,2) ! (combo)
#endif

                 call Combo (dzsno(fc,3), swliq(fc,3), swice(fc,3), tsno(fc,3), drr, &
                      zwliq, zwice, tsno(fc,2))

                 ! Subdivided a new layer
                 if (is_lake) then
                    offset = 2._r8*lsadz
                 else
                    offset = 0._r8
                 end if
                 if (msno <= 3 .and. dzsno(fc,3) > 0.18_r8+offset) then
                    msno =  4
                    dtdz = (tsno(fc,2) - tsno(fc,3))/((dzsno(fc,2)+dzsno(fc,3))/2._r8) 
                    dzsno(fc,3) = dzsno(fc,3)/2._r8
                    swice(fc,3) = swice(fc,3)/2._r8
                    swliq(fc,3) = swliq(fc,3)/2._r8
                    dzsno(fc,4) = dzsno(fc,3)
                    swice(fc,4) = swice(fc,3)
                    swliq(fc,4) = swliq(fc,3)
                    tsno(fc,4) = tsno(fc,3) - dtdz*dzsno(fc,3)/2._r8
                    if (tsno(fc,4) >= tfrz) then 
                       tsno(fc,4)  = tsno(fc,3)
                    else
                       tsno(fc,3) = tsno(fc,3) + dtdz*dzsno(fc,3)/2._r8 
                    endif

                    mbc_phi(fc,3) = mbc_phi(fc,3)/2._r8
                    mbc_phi(fc,4) = mbc_phi(fc,3)
                    mbc_pho(fc,3) = mbc_pho(fc,3)/2._r8
                    mbc_pho(fc,4) = mbc_pho(fc,3)
                    moc_phi(fc,3) = moc_phi(fc,3)/2._r8
                    moc_phi(fc,4) = moc_phi(fc,3)
                    moc_pho(fc,3) = moc_pho(fc,3)/2._r8
                    moc_pho(fc,4) = moc_pho(fc,3)
                    mdst1(fc,3) = mdst1(fc,3)/2._r8
                    mdst1(fc,4) = mdst1(fc,3)
                    mdst2(fc,3) = mdst2(fc,3)/2._r8
                    mdst2(fc,4) = mdst2(fc,3)
                    mdst3(fc,3) = mdst3(fc,3)/2._r8
                    mdst3(fc,4) = mdst3(fc,3)
                    mdst4(fc,3) = mdst4(fc,3)/2._r8
                    mdst4(fc,4) = mdst4(fc,3)
                    rds(fc,4) = rds(fc,3)

                 end if
              end if
           end if

           if (msno > 3) then
              if (is_lake) then
                 offset = lsadz
              else
                 offset = 0._r8
              end if
              if (dzsno(fc,3) > 0.11_r8 + offset) then
                 if (is_lake) then
                    drr = dzsno(fc,3) - 0.11_r8 - lsadz
                 else
                    drr = dzsno(fc,3) - 0.11_r8
                 end if
                 propor = drr/dzsno(fc,3)
                 zwice = propor*swice(fc,3)
                 zwliq = propor*swliq(fc,3)

                 zmbc_phi = propor*mbc_phi(fc,3)
                 zmbc_pho = propor*mbc_pho(fc,3)
                 zmoc_phi = propor*moc_phi(fc,3)
                 zmoc_pho = propor*moc_pho(fc,3)
                 zmdst1 = propor*mdst1(fc,3)
                 zmdst2 = propor*mdst2(fc,3)
                 zmdst3 = propor*mdst3(fc,3)
                 zmdst4 = propor*mdst4(fc,3)

                 if (is_lake) then
                    propor = (0.11_r8+lsadz)/dzsno(fc,3)
                 else
                    propor = 0.11_r8/dzsno(fc,3)
                 end if
                 swice(fc,3) = propor*swice(fc,3)
                 swliq(fc,3) = propor*swliq(fc,3)

                 mbc_phi(fc,3) = propor*mbc_phi(fc,3)
                 mbc_pho(fc,3) = propor*mbc_pho(fc,3)
                 moc_phi(fc,3) = propor*moc_phi(fc,3)
                 moc_pho(fc,3) = propor*moc_pho(fc,3)
                 mdst1(fc,3) = propor*mdst1(fc,3)
                 mdst2(fc,3) = propor*mdst2(fc,3)
                 mdst3(fc,3) = propor*mdst3(fc,3)
                 mdst4(fc,3) = propor*mdst4(fc,3)

                 if (is_lake) then
                    dzsno(fc,3) = 0.11_r8 + lsadz
                 else
                    dzsno(fc,3) = 0.11_r8
                 end if

                 mbc_phi(fc,4) = mbc_phi(fc,4)+zmbc_phi  ! (combo)
                 mbc_pho(fc,4) = mbc_pho(fc,4)+zmbc_pho  ! (combo)
                 moc_phi(fc,4) = moc_phi(fc,4)+zmoc_phi  ! (combo)
                 moc_pho(fc,4) = moc_pho(fc,4)+zmoc_pho  ! (combo)
                 mdst1(fc,4) = mdst1(fc,4)+zmdst1  ! (combo)
                 mdst2(fc,4) = mdst2(fc,4)+zmdst2  ! (combo)
                 mdst3(fc,4) = mdst3(fc,4)+zmdst3  ! (combo)
                 mdst4(fc,4) = mdst4(fc,4)+zmdst4  ! (combo)
#ifdef MODAL_AER
              !mgf++ bugfix
              rds(fc,4) = (rds(fc,4)*(swliq(fc,4)+swice(fc,4)) + rds(fc,3)*(zwliq+zwice))/(swliq(fc,4)+swice(fc,4)+zwliq+zwice)
                if ((rds(fc,4) < 30.) .or. (rds(fc,4) > 1500.)) then
                   write (iulog,*) "4. SNICAR ERROR: snow grain radius of",rds(fc,4),rds(fc,3)
                   write (iulog,*) "swliq, swice, zwliq, zwice", swliq(fc,4), swice(fc,4),zwliq, zwice
                   write (iulog,*) "layers ", msno
                endif
              !mgf--
#else
              rds(fc,4) = rds(fc,3) ! (combo)
#endif

                 call Combo (dzsno(fc,4), swliq(fc,4), swice(fc,4), tsno(fc,4), drr, &
                      zwliq, zwice, tsno(fc,3))

                 ! Subdivided a new layer
                 if (is_lake) then
                    offset = 2._r8*lsadz
                 else
                    offset = 0._r8
                 end if
                 if (msno <= 4 .and. dzsno(fc,4) > 0.41_r8 + offset) then
                    msno = 5
                    dtdz = (tsno(fc,3) - tsno(fc,4))/((dzsno(fc,3)+dzsno(fc,4))/2._r8) 
                    dzsno(fc,4) = dzsno(fc,4)/2._r8
                    swice(fc,4) = swice(fc,4)/2._r8
                    swliq(fc,4) = swliq(fc,4)/2._r8
                    dzsno(fc,5) = dzsno(fc,4)
                    swice(fc,5) = swice(fc,4)
                    swliq(fc,5) = swliq(fc,4)
                    tsno(fc,5) = tsno(fc,4) - dtdz*dzsno(fc,4)/2._r8 
                    if (tsno(fc,5) >= tfrz) then 
                       tsno(fc,5)  = tsno(fc,4)
                    else
                       tsno(fc,4) = tsno(fc,4) + dtdz*dzsno(fc,4)/2._r8 
                    endif

                    mbc_phi(fc,4) = mbc_phi(fc,4)/2._r8
                    mbc_phi(fc,5) = mbc_phi(fc,4)
                    mbc_pho(fc,4) = mbc_pho(fc,4)/2._r8
                    mbc_pho(fc,5) = mbc_pho(fc,4)              
                    moc_phi(fc,4) = moc_phi(fc,4)/2._r8
                    moc_phi(fc,5) = moc_phi(fc,4)
                    moc_pho(fc,4) = moc_pho(fc,4)/2._r8
                    moc_pho(fc,5) = moc_pho(fc,4)
                    mdst1(fc,4) = mdst1(fc,4)/2._r8
                    mdst1(fc,5) = mdst1(fc,4)
                    mdst2(fc,4) = mdst2(fc,4)/2._r8
                    mdst2(fc,5) = mdst2(fc,4)
                    mdst3(fc,4) = mdst3(fc,4)/2._r8
                    mdst3(fc,5) = mdst3(fc,4)
                    mdst4(fc,4) = mdst4(fc,4)/2._r8
                    mdst4(fc,5) = mdst4(fc,4)
                    rds(fc,5) = rds(fc,4)

                 end if
              end if
           end if

           if (msno > 4) then
              if (is_lake) then
                 offset = lsadz
              else
                 offset = 0._r8
              end if
              if (dzsno(fc,4) > 0.23_r8+offset) then
                 if (is_lake) then
                    drr = dzsno(fc,4) - 0.23_r8 - lsadz
                 else
                    drr = dzsno(fc,4) - 0.23_r8
                 end if
                 propor = drr/dzsno(fc,4)
                 zwice = propor*swice(fc,4)
                 zwliq = propor*swliq(fc,4)

                 zmbc_phi = propor*mbc_phi(fc,4)
                 zmbc_pho = propor*mbc_pho(fc,4)
                 zmoc_phi = propor*moc_phi(fc,4)
                 zmoc_pho = propor*moc_pho(fc,4)
                 zmdst1 = propor*mdst1(fc,4)
                 zmdst2 = propor*mdst2(fc,4)
                 zmdst3 = propor*mdst3(fc,4)
                 zmdst4 = propor*mdst4(fc,4)

                 if (is_lake) then
                    propor = (0.23_r8+lsadz)/dzsno(fc,4)
                 else
                    propor = 0.23_r8/dzsno(fc,4)
                 end if
                 swice(fc,4) = propor*swice(fc,4)
                 swliq(fc,4) = propor*swliq(fc,4)

                 mbc_phi(fc,4) = propor*mbc_phi(fc,4)
                 mbc_pho(fc,4) = propor*mbc_pho(fc,4)
                 moc_phi(fc,4) = propor*moc_phi(fc,4)
                 moc_pho(fc,4) = propor*moc_pho(fc,4)
                 mdst1(fc,4) = propor*mdst1(fc,4)
                 mdst2(fc,4) = propor*mdst2(fc,4)
                 mdst3(fc,4) = propor*mdst3(fc,4)
                 mdst4(fc,4) = propor*mdst4(fc,4)

                 if (is_lake) then
                    dzsno(fc,4) = 0.23_r8 + lsadz
                 else
                    dzsno(fc,4) = 0.23_r8
                 end if

                 mbc_phi(fc,5) = mbc_phi(fc,5)+zmbc_phi  ! (combo)
                 mbc_pho(fc,5) = mbc_pho(fc,5)+zmbc_pho  ! (combo)
                 moc_phi(fc,5) = moc_phi(fc,5)+zmoc_phi  ! (combo)
                 moc_pho(fc,5) = moc_pho(fc,5)+zmoc_pho  ! (combo)
                 mdst1(fc,5) = mdst1(fc,5)+zmdst1  ! (combo)
                 mdst2(fc,5) = mdst2(fc,5)+zmdst2  ! (combo)
                 mdst3(fc,5) = mdst3(fc,5)+zmdst3  ! (combo)
                 mdst4(fc,5) = mdst4(fc,5)+zmdst4  ! (combo)
#ifdef MODAL_AER
              !mgf++ bugfix
              rds(fc,5) = (rds(fc,5)*(swliq(fc,5)+swice(fc,5)) + rds(fc,4)*(zwliq+zwice))/(swliq(fc,5)+swice(fc,5)+zwliq+zwice)
                if ((rds(fc,5) < 30.) .or. (rds(fc,5) > 1500.)) then
                   write (iulog,*) "5. SNICAR ERROR: snow grain radius of",rds(c,5),rds(c,4)
                   write (iulog,*) "swliq, swice, zwliq, zwice", swliq(c,5), swice(c,5),zwliq, zwice
                   write (iulog,*) "layers ", msno
                endif
              !mgf--
#else
              rds(fc,5) = rds(fc,4) ! (combo)
#endif

                 call Combo (dzsno(fc,5), swliq(fc,5), swice(fc,5), tsno(fc,5), drr, &
                      zwliq, zwice, tsno(fc,4))
              end if
           end if

           snl(c) = -msno

       end do

       !$acc parallel loop independent gang vector default(present) collapse(2) 
       do j = -nlevsno+1,0
          do fc = 1, num_snowc
             c = filter_snowc(fc)
             if (j >= snl(c)+1) then
                if (is_lake) then
                   dz(c,j) = dzsno(fc,j-snl(c))
                else
                   dz(c,j) = dzsno(fc,j-snl(c))/frac_sno(c)
                end if
                h2osoi_ice(c,j) = swice(fc,j-snl(c))
                h2osoi_liq(c,j) = swliq(fc,j-snl(c))
                t_soisno(c,j)   = tsno(fc,j-snl(c))
                mss_bcphi(c,j)   = mbc_phi(fc,j-snl(c))
                mss_bcpho(c,j)   = mbc_pho(fc,j-snl(c))
                mss_ocphi(c,j)   = moc_phi(fc,j-snl(c))
                mss_ocpho(c,j)   = moc_pho(fc,j-snl(c))
                mss_dst1(c,j)    = mdst1(fc,j-snl(c))
                mss_dst2(c,j)    = mdst2(fc,j-snl(c))
                mss_dst3(c,j)    = mdst3(fc,j-snl(c))
                mss_dst4(c,j)    = mdst4(fc,j-snl(c))
                snw_rds(c,j)     = rds(fc,j-snl(c))

             end if
          end do
       end do

       ! Consistency check
       if (is_lake) then

         do fc = 1, num_snowc
            c = filter_snowc(fc)
            
            do j = -nlevsno + 1, 0

                if (j >= snl(c)+1) then
                   dztot(fc) = dztot(fc) - dz(c,j)
                   snwicetot(fc) = snwicetot(fc) - h2osoi_ice(c,j)
                   snwliqtot(fc) = snwliqtot(fc) - h2osoi_liq(c,j)
                end if

                if (j == 0) then
                   if ( abs(dztot(fc)) > 1.e-10_r8 .or. abs(snwicetot(fc)) > 1.e-7_r8 .or. &
                        abs(snwliqtot(fc)) > 1.e-7_r8 ) then
                      write(iulog,*)'Inconsistency in SnowDivision_Lake! c, remainders', &
                           'dztot, snwicetot, snwliqtot = ',c,dztot(c),snwicetot(c),snwliqtot(c)
                      call endrun(decomp_index=c, elmlevel=namec, msg=errmsg(__FILE__, __LINE__))
                   end if
                end if
            end do
         end do
       end if

       !$acc parallel loop independent gang vector default(present) collapse(2) 
       do j = 0, -nlevsno+1, -1
          do fc = 1, num_snowc
             c = filter_snowc(fc)
             if (j >= snl(c)+1) then
                z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
                zi(c,j-1) = zi(c,j) - dz(c,j)
             end if
          end do
       end do

     !$acc exit data delete(&
     !$acc dzsno(:,:), &
     !$acc swice(:,:), &
     !$acc swliq(:,:), &
     !$acc tsno(:,:), &
     !$acc mbc_phi(:,:), &
     !$acc mbc_pho(:,:), &
     !$acc moc_phi(:,:), &
     !$acc moc_pho(:,:), &
     !$acc mdst1(:,:), &
     !$acc mdst2(:,:), &
     !$acc mdst3(:,:), &
     !$acc mdst4(:,:), &
     !$acc rds(:,:), &
     !$acc dztot(:), &
     !$acc snwicetot(:), &
     !$acc snwliqtot(:), &
     !$acc sum1, &
     !$acc sum2, &
     !$acc sum3)

     end associate

   end subroutine DivideSnowLayers
   
   !-----------------------------------------------------------------------
   subroutine DivideExtraSnowLayers(bounds, num_snowc, filter_snowc, &
        aerosol_vars,  is_lake)
     !
     ! !DESCRIPTION:
     ! Subdivides up to 16 snow layers if they exceed their prescribed maximum thickness.
     ! This alternate subroutine is needed if "use_extrasnowlayers" is true.
     !
     ! !USES:
     use elm_varcon,  only : tfrz 
     use LakeCon   ,  only : lsadz
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds          
     integer                , intent(in)    :: num_snowc       ! number of column snow points in column filter
     integer                , intent(in)    :: filter_snowc(:) ! column filter for snow points
     type(aerosol_type)     , intent(inout) :: aerosol_vars
     logical                , intent(in)    :: is_lake  !TODO - this should be examined and removed in the future
     !
     ! !LOCAL VARIABLES:
     integer  :: j, c, fc, k                              ! indices
     real(r8) :: drr                                      ! thickness of the combined [m]
     integer  :: msno                                     ! number of snow layer 1 (top) to msno (bottom)
     real(r8) :: dzsno(1:num_snowc,nlevsno)   ! Snow layer thickness [m]
     real(r8) :: swice(1:num_snowc,nlevsno)   ! Partial volume of ice [m3/m3]
     real(r8) :: swliq(1:num_snowc,nlevsno)   ! Partial volume of liquid water [m3/m3]
     real(r8) :: tsno(1:num_snowc ,nlevsno)   ! Nodel temperature [K]
     real(r8) :: zwice                                    ! temporary
     real(r8) :: zwliq                                    ! temporary
     real(r8) :: propor                                   ! temporary
     real(r8) :: dtdz                                     ! temporary
     ! temporary variables mimicking the structure of other layer division variables
     real(r8) :: mbc_phi(1:num_snowc,nlevsno) ! mass of BC in each snow layer
     real(r8) :: zmbc_phi                                 ! temporary
     real(r8) :: mbc_pho(1:num_snowc,nlevsno) ! mass of BC in each snow layer
     real(r8) :: zmbc_pho                                 ! temporary
     real(r8) :: moc_phi(1:num_snowc,nlevsno) ! mass of OC in each snow layer
     real(r8) :: zmoc_phi                                 ! temporary
     real(r8) :: moc_pho(1:num_snowc,nlevsno) ! mass of OC in each snow layer
     real(r8) :: zmoc_pho                                 ! temporary
     real(r8) :: mdst1(1:num_snowc,nlevsno)   ! mass of dust 1 in each snow layer
     real(r8) :: zmdst1                                   ! temporary
     real(r8) :: mdst2(1:num_snowc,nlevsno)   ! mass of dust 2 in each snow layer
     real(r8) :: zmdst2                                   ! temporary
     real(r8) :: mdst3(1:num_snowc,nlevsno)   ! mass of dust 3 in each snow layer
     real(r8) :: zmdst3                                   ! temporary
     real(r8) :: mdst4(1:num_snowc,nlevsno)   ! mass of dust 4 in each snow layer
     real(r8) :: zmdst4                                   ! temporary
     real(r8) :: rds(1:num_snowc,nlevsno)
     ! Variables for consistency check
     real(r8) :: dztot(1:num_snowc)
     real(r8) :: snwicetot(1:num_snowc)
     real(r8) :: snwliqtot(1:num_snowc)
     real(r8) :: offset ! temporary
     !-----------------------------------------------------------------------
     
     associate(                                            & 
          t_soisno   => col_es%t_soisno    , & ! Output: [real(r8) (:,:) ] soil temperature (Kelvin)              

          h2osoi_ice => col_ws%h2osoi_ice   , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)                       
          h2osoi_liq => col_ws%h2osoi_liq   , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)                   
          frac_sno   => col_ws%frac_sno_eff , & ! Output: [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
          snw_rds    => col_ws%snw_rds      , & ! Output: [real(r8) (:,:) ] effective snow grain radius (col,lyr) [microns, m^-6]

          mss_bcphi  => aerosol_vars%mss_bcphi_col       , & ! Output: [real(r8) (:,:) ] hydrophilic BC mass in snow (col,lyr) [kg]
          mss_bcpho  => aerosol_vars%mss_bcpho_col       , & ! Output: [real(r8) (:,:) ] hydrophobic BC mass in snow (col,lyr) [kg]
          mss_ocphi  => aerosol_vars%mss_ocphi_col       , & ! Output: [real(r8) (:,:) ] hydrophilic OC mass in snow (col,lyr) [kg]
          mss_ocpho  => aerosol_vars%mss_ocpho_col       , & ! Output: [real(r8) (:,:) ] hydrophobic OC mass in snow (col,lyr) [kg]
          mss_dst1   => aerosol_vars%mss_dst1_col        , & ! Output: [real(r8) (:,:) ] dust species 1 mass in snow (col,lyr) [kg]
          mss_dst2   => aerosol_vars%mss_dst2_col        , & ! Output: [real(r8) (:,:) ] dust species 2 mass in snow (col,lyr) [kg]
          mss_dst3   => aerosol_vars%mss_dst3_col        , & ! Output: [real(r8) (:,:) ] dust species 3 mass in snow (col,lyr) [kg]
          mss_dst4   => aerosol_vars%mss_dst4_col        , & ! Output: [real(r8) (:,:) ] dust species 4 mass in snow (col,lyr) [kg]

          snl        => col_pp%snl                          , & ! Output: [integer  (:)   ] number of snow layers                     
          dz         => col_pp%dz                           , & ! Output: [real(r8) (:,:) ] layer depth (m)                        
          zi         => col_pp%zi                           , & ! Output: [real(r8) (:,:) ] interface level below a "z" level (m)  
          z          => col_pp%z                              & ! Output: [real(r8) (:,:) ] layer thickness (m)                   
          )
     !$acc enter data create(&
     !$acc dzsno(:,:), &
     !$acc swice(:,:), &
     !$acc swliq(:,:), &
     !$acc tsno(:,:), &
     !$acc mbc_phi(:,:), &
     !$acc mbc_pho(:,:), &
     !$acc moc_phi(:,:), &
     !$acc moc_pho(:,:), &
     !$acc mdst1(:,:), &
     !$acc mdst2(:,:), &
     !$acc mdst3(:,:), &
     !$acc mdst4(:,:), &
     !$acc rds(:,:), &
     !$acc dztot(:), &
     !$acc snwicetot(:), &
     !$acc snwliqtot(:), &
     !$acc msno)

      
         

       if ( is_lake ) then
          ! Initialize for consistency check
          do j = -nlevsno+1,0
             do fc = 1, num_snowc
                c = filter_snowc(fc)
                
                if (j == -nlevsno+1) then
                   dztot(fc) = 0._r8
                   snwicetot(fc) = 0._r8
                   snwliqtot(fc) = 0._r8
                end if
                
                if (j >= snl(c)+1) then
                   dztot(fc) = dztot(fc) + dz(c,j)
                   snwicetot(fc) = snwicetot(fc) + h2osoi_ice(c,j)
                   snwliqtot(fc) = snwliqtot(fc) + h2osoi_liq(c,j)
                end if
             end do
          end do
       end if

       ! Begin calculation - note that the following column loops are only invoked
       ! for snow-covered columns

       !$acc parallel loop independent gang vector default(present) collapse(2) 
       do j = 1,nlevsno
          do fc = 1, num_snowc
             c = filter_snowc(fc)
             if (j <= abs(snl(c))) then
                if (is_lake) then
                   dzsno(fc,j) = dz(c,j+snl(c))
                else
                   dzsno(fc,j) = frac_sno(c)*dz(c,j+snl(c))
                end if
                swice(fc,j) = h2osoi_ice(c,j+snl(c))
                swliq(fc,j) = h2osoi_liq(c,j+snl(c))
                tsno(fc,j)  = t_soisno(c,j+snl(c))

                mbc_phi(fc,j) = mss_bcphi(c,j+snl(c))
                mbc_pho(fc,j) = mss_bcpho(c,j+snl(c))
                moc_phi(fc,j) = mss_ocphi(c,j+snl(c))
                moc_pho(fc,j) = mss_ocpho(c,j+snl(c))
                mdst1(fc,j)   = mss_dst1(c,j+snl(c))
                mdst2(fc,j)   = mss_dst2(c,j+snl(c))
                mdst3(fc,j)   = mss_dst3(c,j+snl(c))
                mdst4(fc,j)   = mss_dst4(c,j+snl(c))
                rds(fc,j)     = snw_rds(c,j+snl(c))
             end if
          end do
       end do

       !$acc parallel loop independent gang vector default(present)
       loop_snowcolumns: do fc = 1, num_snowc
           c = filter_snowc(fc)

           msno = abs(snl(c))

           ! Now traverse layers from top to bottom in a dynamic way, as the total
           ! number of layers (msno) may increase during the loop.
           ! Impose k < nlevsno; the special case 'k == nlevsno' is not relevant,
           ! as it is neither allowed to subdivide nor does it have layers below.
           k = 1
           loop_layers: do while( k <= msno .and. k < nlevsno )

              ! Current layer is bottom layer
              if (k == msno) then

                 if (is_lake) then
                    offset = 2._r8 * lsadz
                 else
                    offset = 0._r8
                 end if

                 if (dzsno(fc,k) > dzmax_l(k) + offset) then
                    ! Subdivide layer into two layers with equal thickness, water
                    ! content, ice content and temperature
                    msno = msno + 1
                    dzsno(fc,k)     = dzsno(fc,k) / 2.0_r8
                    dzsno(fc,k+1)   = dzsno(fc,k)
                    swice(fc,k)     = swice(fc,k) / 2.0_r8
                    swice(fc,k+1)   = swice(fc,k)
                    swliq(fc,k)     = swliq(fc,k) / 2.0_r8
                    swliq(fc,k+1)   = swliq(fc,k)

                    if (k == 1) then
                       ! special case
                       tsno(fc,k+1)    = tsno(fc,k)
                    else
                       ! use temperature gradient
                       dtdz           = (tsno(fc,k-1) - tsno(fc,k))/((dzsno(fc,k-1)+2*dzsno(fc,k))/2.0_r8)
                       tsno(fc,k+1) = tsno(fc,k) - dtdz*dzsno(fc,k)/2.0_r8
                       if (tsno(fc,k+1) >= tfrz) then
                          tsno(fc,k+1)  = tsno(fc,k)
                       else
                          tsno(fc,k) = tsno(fc,k) + dtdz*dzsno(fc,k)/2.0_r8
                       endif
                    end if

                    mbc_phi(fc,k)   = mbc_phi(fc,k) / 2.0_r8
                    mbc_phi(fc,k+1) = mbc_phi(fc,k)
                    mbc_pho(fc,k)   = mbc_pho(fc,k) / 2.0_r8
                    mbc_pho(fc,k+1) = mbc_pho(fc,k)
                    moc_phi(fc,k)   = moc_phi(fc,k) / 2.0_r8
                    moc_phi(fc,k+1) = moc_phi(fc,k)
                    moc_pho(fc,k)   = moc_pho(fc,k) / 2.0_r8
                    moc_pho(fc,k+1) = moc_pho(fc,k)
                    mdst1(fc,k)     = mdst1(fc,k) / 2.0_r8
                    mdst1(fc,k+1)   = mdst1(fc,k)
                    mdst2(fc,k)     = mdst2(fc,k) / 2.0_r8
                    mdst2(fc,k+1)   = mdst2(fc,k)
                    mdst3(fc,k)     = mdst3(fc,k) / 2.0_r8
                    mdst3(fc,k+1)   = mdst3(fc,k)
                    mdst4(fc,k)     = mdst4(fc,k) / 2.0_r8
                    mdst4(fc,k+1)   = mdst4(fc,k)

                    rds(fc,k+1)     = rds(fc,k)
                 end if
              end if

              ! There are layers below (note this is not exclusive with previous
              ! if-statement, since msno may have increased in the previous if-statement)
              if (k < msno) then

                 if (is_lake) then
                    offset = lsadz
                 else
                    offset = 0._r8
                 end if

                 if (dzsno(fc,k) > dzmax_u(k) + offset ) then
                    ! Only dump excess snow to underlying layer in a conservative fashion.
                    ! Other quantities will depend on the height of the excess snow: a ratio is used for this.
                    drr      = dzsno(fc,k) - dzmax_u(k) - offset

                    propor   = drr/dzsno(fc,k)
                    zwice    = propor*swice(fc,k)
                    zwliq    = propor*swliq(fc,k)
                    zmbc_phi = propor*mbc_phi(fc,k)
                    zmbc_pho = propor*mbc_pho(fc,k)
                    zmoc_phi = propor*moc_phi(fc,k)
                    zmoc_pho = propor*moc_pho(fc,k)
                    zmdst1   = propor*mdst1(fc,k)
                    zmdst2   = propor*mdst2(fc,k)
                    zmdst3   = propor*mdst3(fc,k)
                    zmdst4   = propor*mdst4(fc,k)

                    propor         = (dzmax_u(k)+offset)/dzsno(fc,k)
                    swice(fc,k)     = propor*swice(fc,k)
                    swliq(fc,k)     = propor*swliq(fc,k)
                    mbc_phi(fc,k)   = propor*mbc_phi(fc,k)
                    mbc_pho(fc,k)   = propor*mbc_pho(fc,k)
                    moc_phi(fc,k)   = propor*moc_phi(fc,k)
                    moc_pho(fc,k)   = propor*moc_pho(fc,k)
                    mdst1(fc,k)     = propor*mdst1(fc,k)
                    mdst2(fc,k)     = propor*mdst2(fc,k)
                    mdst3(fc,k)     = propor*mdst3(fc,k)
                    mdst4(fc,k)     = propor*mdst4(fc,k)

                    ! Set depth layer k to maximum allowed value
                    dzsno(fc,k)  = dzmax_u(k)  + offset

                    mbc_phi(fc,k+1) = mbc_phi(fc,k+1)+zmbc_phi  ! (combo)
                    mbc_pho(fc,k+1) = mbc_pho(fc,k+1)+zmbc_pho  ! (combo)
                    moc_phi(fc,k+1) = moc_phi(fc,k+1)+zmoc_phi  ! (combo)
                    moc_pho(fc,k+1) = moc_pho(fc,k+1)+zmoc_pho  ! (combo)
                    mdst1(fc,k+1)   = mdst1(fc,k+1)+zmdst1  ! (combo)
                    mdst2(fc,k+1)   = mdst2(fc,k+1)+zmdst2  ! (combo)
                    mdst3(fc,k+1)   = mdst3(fc,k+1)+zmdst3  ! (combo)
                    mdst4(fc,k+1)   = mdst4(fc,k+1)+zmdst4  ! (combo)

                    ! Mass-weighted combination of radius
                    rds(fc,k+1) = MassWeightedSnowRadius( rds(fc,k), rds(fc,k+1), &
                         (swliq(fc,k+1)+swice(fc,k+1)), (zwliq+zwice) )

                    call Combo (dzsno(fc,k+1), swliq(fc,k+1), swice(fc,k+1), tsno(fc,k+1), drr, &
                         zwliq, zwice, tsno(fc,k))
                 end if
              end if
              k = k+1
           end do loop_layers

           snl(c) = -msno

       end do loop_snowcolumns

       !$acc parallel loop independent gang vector default(present) collapse(2) 
       do j = -nlevsno+1,0
          do fc = 1, num_snowc
             c = filter_snowc(fc)
             if (j >= snl(c)+1) then
                if (is_lake) then
                   dz(c,j) = dzsno(fc,j-snl(c))
                else
                   dz(c,j) = dzsno(fc,j-snl(c))/frac_sno(c)
                end if
                h2osoi_ice(c,j) = swice(fc,j-snl(c))
                h2osoi_liq(c,j) = swliq(fc,j-snl(c))
                t_soisno(c,j)   = tsno(fc,j-snl(c))
                mss_bcphi(c,j)   = mbc_phi(fc,j-snl(c))
                mss_bcpho(c,j)   = mbc_pho(fc,j-snl(c))
                mss_ocphi(c,j)   = moc_phi(fc,j-snl(c))
                mss_ocpho(c,j)   = moc_pho(fc,j-snl(c))
                mss_dst1(c,j)    = mdst1(fc,j-snl(c))
                mss_dst2(c,j)    = mdst2(fc,j-snl(c))
                mss_dst3(c,j)    = mdst3(fc,j-snl(c))
                mss_dst4(c,j)    = mdst4(fc,j-snl(c))
                snw_rds(c,j)     = rds(fc,j-snl(c))
             end if
          end do
       end do

       ! Consistency check
       if (is_lake) then
          do j = -nlevsno + 1, 0
             do fc = 1, num_snowc
                c = filter_snowc(fc)

                if (j >= snl(c)+1) then
                   dztot(fc) = dztot(fc) - dz(c,j)
                   snwicetot(fc) = snwicetot(fc) - h2osoi_ice(c,j)
                   snwliqtot(fc) = snwliqtot(fc) - h2osoi_liq(c,j)
                end if

                if (j == 0) then
                   if ( abs(dztot(fc)) > 1.e-10_r8 .or. abs(snwicetot(fc)) > 1.e-7_r8 .or. &
                        abs(snwliqtot(fc)) > 1.e-7_r8 ) then
#ifndef _OPENACC
                      write(iulog,*)'Inconsistency in SnowDivision_Lake! c, remainders', &
                           'dztot, snwicetot, snwliqtot = ',c,dztot(fc),snwicetot(fc),snwliqtot(fc)
                      call endrun(decomp_index=c, elmlevel=namec, msg=errmsg(__FILE__, __LINE__))
#endif
                   end if
                end if
             end do
          end do
       end if

       !$acc parallel loop independent gang vector default(present) collapse(2) 
       do j = 0, -nlevsno+1, -1
          do fc = 1, num_snowc
             c = filter_snowc(fc)
             if (j >= snl(c)+1) then
                z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
                zi(c,j-1) = zi(c,j) - dz(c,j)
             end if
          end do
       end do

     !$acc exit data delete(&
     !$acc dzsno(:,:), &
     !$acc swice(:,:), &
     !$acc swliq(:,:), &
     !$acc tsno(:,:), &
     !$acc mbc_phi(:,:), &
     !$acc mbc_pho(:,:), &
     !$acc moc_phi(:,:), &
     !$acc moc_pho(:,:), &
     !$acc mdst1(:,:), &
     !$acc mdst2(:,:), &
     !$acc mdst3(:,:), &
     !$acc mdst4(:,:), &
     !$acc rds(:,:), &
     !$acc dztot(:), &
     !$acc snwicetot(:), &
     !$acc snwliqtot(:), &
     !$acc msno)

     end associate

   end subroutine DivideExtraSnowLayers   
   
   !-----------------------------------------------------------------------
   subroutine InitSnowLayers (bounds, snow_depth)
    ! (from CLMv5)
    ! !DESCRIPTION:
    ! Initialize snow layer depth from specified total depth.
    !
    ! !USES:
    use elm_varcon         , only : spval
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    real(r8)               , intent(in)    :: snow_depth(bounds%begc:)
    !
    ! LOCAL VARAIBLES:
    integer               :: c,l,j              ! indices
    real(r8)              :: minbound, maxbound ! helper variables
   
    !SHR_ASSERT_ALL((ubound(snow_depth)  == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate( &
         snl => col_pp%snl,   & ! Output: [integer (:)    ]  number of snow layers
         dz  => col_pp%dz,    & ! Output: [real(r8) (:,:) ]  layer thickness (m)  (-nlevsno+1:nlevgrnd)
         z   => col_pp%z,     & ! Output: [real(r8) (:,:) ]  layer depth (m) (-nlevsno+1:nlevgrnd)
         zi  => col_pp%zi     & ! Output: [real(r8) (:,:) ]  interface level below a "z" level (m) (-nlevsno+0:nlevgrnd)
    )
    loop_columns: do c = bounds%begc,bounds%endc
       l = col_pp%landunit(c)

       dz(c,-nlevsno+1: 0) = spval
       z (c,-nlevsno+1: 0) = spval
       zi(c,-nlevsno  :-1) = spval
    
       ! Special case: lake
       if (lun_pp%lakpoi(l)) then
          snl(c)              = 0
          dz(c,-nlevsno+1:0)  = 0._r8
          z(c,-nlevsno+1:0)   = 0._r8
          zi(c,-nlevsno+0:0)  = 0._r8
          cycle
       end if
       
       ! LvK 9-JUN-2015: in CanopyHydrologyMod , snow_depth is scaled with frac_sno
       ! Here we do not apply scaling to snow_depth, so inconsistent? TODO
       
       ! Special case: too little snow for snowpack existence
       if (snow_depth(c) < dzmin16(1)) then
          snl(c)              = 0
          dz(c,-nlevsno+1:0)  = 0._r8
          z(c,-nlevsno+1:0)   = 0._r8
          zi(c,-nlevsno+0:0)  = 0._r8
          cycle
       end if
       
       ! There has to be at least one snow layer
       snl(c)   = -1
       minbound = dzmin16(1)
       maxbound = dzmax_l(1)

       if (snow_depth(c) >= minbound .and. snow_depth(c) <= maxbound) then
          ! Special case: single layer
          dz(c,0) = snow_depth(c)

       else
          ! Search for appropriate number of layers (snl) by increasing the
          ! number of layers and check for matching bounds.
          snl(c) = snl(c) - 1
          minbound = maxbound
          maxbound = sum(dzmax_u(1:-snl(c)))

          do while(snow_depth(c) > maxbound .and. -snl(c) < nlevsno )
             snl(c) = snl(c) - 1
             minbound = maxbound
             maxbound = sum(dzmax_u(1:-snl(c)))
          end do

          ! Set thickness of all layers except bottom two
          do j = 1, -snl(c)-2
             dz(c,j+snl(c))  = dzmax_u(j)
          end do

          ! Determine whether the two bottom layers should be equal in size,
          ! or not. The rule here is: always create equal size when possible.
          if (snow_depth(c) <= sum(dzmax_u(1:-snl(c)-2)) + 2 * dzmax_u(-snl(c)-1)) then
             dz(c,-1) = (snow_depth(c) - sum(dzmax_u(1:-snl(c)-2))) / 2._r8
             dz(c,0)  = dz(c,-1)
          else
             dz(c,-1) = dzmax_u(-snl(c)-1)
             dz(c,0)  = snow_depth(c) - sum(dzmax_u(1:-snl(c)-1))
          endif
       endif

       ! Initialize the node depth and the depth of layer interface
       do j = 0, snl(c)+1, -1
          z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
          zi(c,j-1) = zi(c,j) - dz(c,j)
       end do

    end do loop_columns

    end associate
   end subroutine InitSnowLayers
   
   !-----------------------------------------------------------------------
   subroutine SnowCapping(bounds, num_nolakec, filter_initc, num_snowc, filter_snowc, &
        aerosol_vars  )
     ! from CLM (clm4_5_12_r190)
     ! !DESCRIPTION:
     ! Removes mass from bottom snow layer for columns that exceed the maximum snow depth.
     ! This routine is called twice: once for non-lake columns and once for lake columns. 
     ! The initialization of the snow capping fluxes should only be done ONCE for each group,
     ! therefore they are a passed as an extra argument (filter_initc). 
     ! Density and temperature of the layer are conserved (density needs some work, temperature is a state
     ! variable)
     !
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds          
     integer                , intent(in)    :: num_nolakec       ! number of column points that need to be initialized
     integer                , intent(in)    :: filter_initc(:) ! column filter for points that need to be initialized
     integer                , intent(in)    :: num_snowc       ! number of column snow points in column filter
     integer                , intent(in)    :: filter_snowc(:) ! column filter for snow points
     type(aerosol_type)     , intent(inout) :: aerosol_vars
     !
     ! !LOCAL VARIABLES:
     real(r8)   :: dtime                            ! land model time step (sec)
     real(r8)   :: mss_snwcp_tot                    ! total snow capping mass [kg/m2] 
     real(r8)   :: mss_snow_bottom_lyr              ! total snow mass (ice+liquid) in bottom layer [kg/m2]
     real(r8)   :: icefrac                          ! fraction of ice mass w.r.t. total mass [unitless]
     real(r8)   :: frac_adjust                      ! fraction of mass remaining after capping
     real(r8)   :: rho                              ! partial density of ice (not scaled with frac_sno) [kg/m3]
     integer    :: fc, c                            ! counters
     ! Always keep at least this fraction of the bottom snow layer when doing snow capping
     ! This needs to be slightly greater than 0 to avoid roundoff problems
     real(r8), parameter :: min_snow_to_keep = 1.e-9  ! fraction of bottom snow layer to keep with capping

     !-----------------------------------------------------------------------
     associate( &
         qflx_snwcp_ice     => col_wf%qflx_snwcp_ice               , & ! Output: [real(r8) (:)   ]  excess solid h2o due to snow capping (outgoing) (mm H2O /s) [+]
         qflx_snwcp_liq     => col_wf%qflx_snwcp_liq               , & ! Output: [real(r8) (:)   ]  excess liquid h2o due to snow capping (outgoing) (mm H2O /s) [+]
         h2osoi_ice         => col_ws%h2osoi_ice                   , & ! In/Out: [real(r8) (:,:) ] ice lens (kg/m2)                       
         h2osoi_liq         => col_ws%h2osoi_liq                   , & ! In/Out: [real(r8) (:,:) ] liquid water (kg/m2)                   
         h2osno             => col_ws%h2osno                       , & ! Input:  [real(r8) (:)   ] snow water (mm H2O)
         mss_bcphi          => aerosol_vars%mss_bcphi_col          , & ! In/Out: [real(r8) (:,:) ] hydrophilic BC mass in snow (col,lyr) [kg]
         mss_bcpho          => aerosol_vars%mss_bcpho_col          , & ! In/Out: [real(r8) (:,:) ] hydrophobic BC mass in snow (col,lyr) [kg]
         mss_ocphi          => aerosol_vars%mss_ocphi_col          , & ! In/Out: [real(r8) (:,:) ] hydrophilic OC mass in snow (col,lyr) [kg]
         mss_ocpho          => aerosol_vars%mss_ocpho_col          , & ! In/Out: [real(r8) (:,:) ] hydrophobic OC mass in snow (col,lyr) [kg]
         mss_dst1           => aerosol_vars%mss_dst1_col           , & ! In/Out: [real(r8) (:,:) ] dust species 1 mass in snow (col,lyr) [kg]
         mss_dst2           => aerosol_vars%mss_dst2_col           , & ! In/Out: [real(r8) (:,:) ] dust species 2 mass in snow (col,lyr) [kg]
         mss_dst3           => aerosol_vars%mss_dst3_col           , & ! In/Out: [real(r8) (:,:) ] dust species 3 mass in snow (col,lyr) [kg]
         mss_dst4           => aerosol_vars%mss_dst4_col           , & ! In/Out: [real(r8) (:,:) ] dust species 4 mass in snow (col,lyr) [kg]
         dz                 => col_pp%dz                             & ! In/Out: [real(r8) (:,:) ] layer thickness (m)
     )

     ! Determine model time step
     dtime = dtime_mod 
     ! Initialize capping fluxes for all columns in domain (lake or non-lake)
     do fc = 1, num_nolakec
        c = filter_initc(fc)
        qflx_snwcp_ice(c) = 0.0_r8
        qflx_snwcp_liq(c) = 0.0_r8
     end do

     loop_columns: do fc = 1, num_snowc
        c = filter_snowc(fc)

        if (h2osno(c) > h2osno_max) then
           mss_snow_bottom_lyr = h2osoi_ice(c,0) + h2osoi_liq(c,0) 
           mss_snwcp_tot = min(h2osno(c)-h2osno_max, mss_snow_bottom_lyr * (1._r8 - min_snow_to_keep)) ! Can't remove more mass than available

           ! Ratio of snow/liquid in bottom layer determines partitioning of runoff fluxes
           icefrac = h2osoi_ice(c,0) / mss_snow_bottom_lyr
           qflx_snwcp_ice(c) = mss_snwcp_tot/dtime * icefrac
           qflx_snwcp_liq(c) = mss_snwcp_tot/dtime * (1._r8 - icefrac)

           rho = h2osoi_ice(c,0) / dz(c,0) ! ice only

           ! Adjust water content
           h2osoi_ice(c,0) = h2osoi_ice(c,0) - qflx_snwcp_ice(c)*dtime
           h2osoi_liq(c,0) = h2osoi_liq(c,0) - qflx_snwcp_liq(c)*dtime

           ! Scale dz such that ice density (or: pore space) is conserved
           !
           ! Avoid scaling dz for very low ice densities. This can occur, in principle, if
           ! the layer is mostly liquid water. Furthermore, this check is critical in the
           ! unlikely event that rho is 0, which can happen if the layer is entirely liquid
           ! water.
           if (rho > 1.0_r8) then
             dz(c,0) = h2osoi_ice(c,0) / rho 
           end if

           ! Check that water capacity is still positive
           if (h2osoi_ice(c,0) < 0._r8 .or. h2osoi_liq(c,0) < 0._r8 ) then
              write(iulog,*)'ERROR: capping procedure failed (negative mass remaining) c = ',c
              write(iulog,*)'h2osoi_ice = ', h2osoi_ice(c,0), ' h2osoi_liq = ', h2osoi_liq(c,0)
              call endrun(decomp_index=c, elmlevel=namec, msg=errmsg(__FILE__, __LINE__))
           end if

           ! Correct the top layer aerosol mass to account for snow capping.
           ! This approach conserves the aerosol mass concentration but not aerosol mass. 
           frac_adjust = (mss_snow_bottom_lyr - mss_snwcp_tot) / mss_snow_bottom_lyr
           mss_bcphi(c,0)   = mss_bcphi(c,0) * frac_adjust 
           mss_bcpho(c,0)   = mss_bcpho(c,0) * frac_adjust
           mss_ocphi(c,0)   = mss_ocphi(c,0) * frac_adjust
           mss_ocpho(c,0)   = mss_ocpho(c,0) * frac_adjust
           mss_dst1(c,0)    = mss_dst1(c,0) * frac_adjust
           mss_dst2(c,0)    = mss_dst2(c,0) * frac_adjust
           mss_dst3(c,0)    = mss_dst3(c,0) * frac_adjust
           mss_dst4(c,0)    = mss_dst4(c,0) * frac_adjust
        end if

     end do loop_columns

     end associate
   end subroutine SnowCapping
   
   !-----------------------------------------------------------------------
   subroutine NewSnowBulkDensity(bounds, num_c, filter_c, top_as, bifall)
      ! (subroutine from CLMv5)
      ! !DESCRIPTION:
      ! Compute the bulk density of any newly-fallen snow.
      !
      ! The return value is placed in bifall. Only columns within the given filter are set:
      ! all other columns remain at their original values.
      !
      ! !USES:
      use elm_varcon,  only : tfrz
      !
      ! !ARGUMENTS:
      type(bounds_type)  , intent(in)    :: bounds
      integer            , intent(in)    :: num_c                ! number of columns in filterc
      integer            , intent(in)    :: filter_c(:)          ! column-level filter to operate on
      type(topounit_atmospheric_state) , intent(in)   :: top_as  
      real(r8)           , intent(inout) :: bifall(bounds%begc:) ! bulk density of newly fallen dry snow [kg/m3]
      !
      ! !LOCAL VARIABLES:
      integer :: fc, c, t!, g
      real(r8) :: t_for_bifall_degC  ! temperature to use in bifall equation (deg C)

      character(len=*), parameter :: subname = 'NewSnowBulkDensity'
      !-----------------------------------------------------------------------

      associate(forc_t    => top_as%tbot   , & ! Input:  [real(r8) (:) ]  atmospheric temperature (Kelvin)
                forc_wind => top_as%windbot  & ! Input:  [real(r8) (:) ]  atmospheric wind speed (m/s)
                )

      do fc = 1, num_c
         c = filter_c(fc)
         t = col_pp%topounit(c)

         if (forc_t(t) > tfrz + 2._r8) then
            bifall(c) = 50._r8 + 1.7_r8*(17.0_r8)**1.5_r8
         else if (forc_t(t) > tfrz - 15._r8) then
            bifall(c) = 50._r8 + 1.7_r8*(forc_t(t) - tfrz + 15._r8)**1.5_r8
         else ! below comment from CLMv5
            ! Andrew Slater: A temp of about -15C gives the nicest
            ! "blower" powder, but as you get colder the flake size decreases so
            ! density goes up. e.g. the smaller snow crystals from the Arctic and Antarctic
            ! winters
            if (forc_t(t) > tfrz - 57.55_r8) then
               t_for_bifall_degC = (forc_t(t)-tfrz)
            else ! below comment from CLMv5
               ! Below -57.55 deg C, the following function starts to decrease with
               ! decreasing temperatures. Limit the function to avoid this turning over.
               t_for_bifall_degC = -57.55_r8
            end if
            bifall(c) = -(50._r8/15._r8 + 0.0333_r8*15_r8)*t_for_bifall_degC - 0.0333_r8*t_for_bifall_degC**2
         end if

         if (forc_wind(t) > 0.1_r8) then
            ! Density offset for wind-driven compaction, initial ideas based on Liston et. al (2007) J. Glaciology,
            ! 53(181), 241-255. Modified for a continuous wind impact and slightly more sensitive
            ! to wind - Andrew Slater, 2016
            bifall(c) = bifall(c) + (266.861_r8 * ((1._r8 + TANH(forc_wind(t)/5.0_r8))/2._r8)**8.8_r8)
         end if

      end do

      end associate

   end subroutine NewSnowBulkDensity
   
   !-----------------------------------------------------------------------
   subroutine WindDriftCompaction(bi, forc_wind, dz, &
        zpseudo, mobile, compaction_rate)
     ! from CLMv5
     ! !DESCRIPTION:
     !
     ! Compute wind drift compaction for a single column and level.
     !
     ! Also updates zpseudo and mobile for this column. However, zpseudo remains unchanged
     ! if mobile is already false or becomes false within this subroutine.
     !
     ! The structure of the updates done here for zpseudo and mobile requires that this
     ! subroutine be called first for the top layer of snow, then for the 2nd layer down,
     ! etc. - and finally for the bottom layer. Before beginning the loops over layers,
     ! mobile should be initialized to .true. and zpseudo should be initialized to 0.
     !
     ! !USES:
     !
     ! !ARGUMENTS:
      !$acc routine seq 
     real(r8) , intent(in)    :: bi              ! partial density of ice [kg/m3]
     real(r8) , intent(in)    :: forc_wind       ! atmospheric wind speed [m/s]
     real(r8) , intent(in)    :: dz              ! layer depth for this column and level [m]
     real(r8) , intent(inout) :: zpseudo         ! wind drift compaction / pseudo depth for this column at this layer
     logical  , intent(inout) :: mobile          ! whether this snow column is still mobile at this layer (i.e., susceptible to wind drift)
     real(r8) , intent(out)   :: compaction_rate ! rate of compaction of snowpack due to wind drift, for the current column and layer
     !
     ! !LOCAL VARIABLES:
     real(r8) :: Frho        ! Mobility density factor [-]
     real(r8) :: MO          ! Mobility index [-]
     real(r8) :: SI          ! Driftability index [-]
     real(r8) :: gamma_drift ! Scaling factor for wind drift time scale [-]
     real(r8) :: tau_inverse ! Inverse of the effective time scale [1/s]

     real(r8), parameter :: rho_min = 50._r8      ! wind drift compaction / minimum density [kg/m3]
     real(r8), parameter :: rho_max = 350._r8     ! wind drift compaction / maximum density [kg/m3]
     real(r8), parameter :: drift_gs = 0.35e-3_r8 ! wind drift compaction / grain size (fixed value for now)
     real(r8), parameter :: drift_sph = 1.0_r8    ! wind drift compaction / sphericity
     real(r8), parameter :: tau_ref = 48._r8 * 3600._r8  ! wind drift compaction / reference time [s]

     character(len=*), parameter :: subname = 'WindDriftCompaction'
     !-----------------------------------------------------------------------

     if (mobile) then
        Frho = 1.25_r8 - 0.0042_r8*(max(rho_min, bi)-rho_min)
        ! assuming dendricity = 0, sphericity = 1, grain size = 0.35 mm Non-dendritic snow
        MO = 0.34_r8 * (-0.583_r8*drift_gs - 0.833_r8*drift_sph + 0.833_r8) + 0.66_r8*Frho
        SI = -2.868_r8 * exp(-0.085_r8*forc_wind) + 1._r8 + MO

        if (SI > 0.0_r8) then
           SI = min(SI, 3.25_r8)
           ! Increase zpseudo (wind drift / pseudo depth) to the middle of
           ! the pseudo-node for the sake of the following calculation
           zpseudo = zpseudo + 0.5_r8 * dz * (3.25_r8 - SI)
           gamma_drift = SI*exp(-zpseudo/0.1_r8)
           tau_inverse = gamma_drift / tau_ref
           compaction_rate = -max(0.0_r8, rho_max-bi) * tau_inverse
           ! Further increase zpseudo to the bottom of the pseudo-node for
           ! the sake of calculations done on the underlying layer (i.e.,
           ! the next time through the j loop).
           zpseudo = zpseudo + 0.5_r8 * dz * (3.25_r8 - SI)
        else  ! SI <= 0
           mobile = .false.
           compaction_rate = 0._r8
        end if
     else  ! .not. mobile
        compaction_rate = 0._r8
     end if

   end subroutine WindDriftCompaction 
   
   !-----------------------------------------------------------------------
   subroutine Combo(dz,  wliq,  wice, t, dz2, wliq2, wice2, t2)
     !
     ! !DESCRIPTION:
     ! Combines two elements and returns the following combined
     ! variables: dz, t, wliq, wice.
     ! The combined temperature is based on the equation:
     ! the sum of the enthalpies of the two elements =
     ! that of the combined element.
     !
     ! !USES:
      !$acc routine seq
     use elm_varcon,  only : cpice, cpliq, tfrz, hfus
     !
     ! !ARGUMENTS:
     implicit none
     real(r8), intent(in)    :: dz2   ! nodal thickness of 2 elements being combined [m]
     real(r8), intent(in)    :: wliq2 ! liquid water of element 2 [kg/m2]
     real(r8), intent(in)    :: wice2 ! ice of element 2 [kg/m2]
     real(r8), intent(in)    :: t2    ! nodal temperature of element 2 [K]
     real(r8), intent(inout) :: dz    ! nodal thickness of 1 elements being combined [m]
     real(r8), intent(inout) :: wliq  ! liquid water of element 1
     real(r8), intent(inout) :: wice  ! ice of element 1 [kg/m2]
     real(r8), intent(inout) :: t     ! nodel temperature of elment 1 [K]
     !
     ! !LOCAL VARIABLES:
     real(r8) :: dzc   ! Total thickness of nodes 1 and 2 (dzc=dz+dz2).
     real(r8) :: wliqc ! Combined liquid water [kg/m2]
     real(r8) :: wicec ! Combined ice [kg/m2]
     real(r8) :: tc    ! Combined node temperature [K]
     real(r8) :: h     ! enthalpy of element 1 [J/m2]
     real(r8) :: h2    ! enthalpy of element 2 [J/m2]
     real(r8) :: hc    ! temporary
     !-----------------------------------------------------------------------

     dzc = dz+dz2
     wicec = (wice+wice2)
     wliqc = (wliq+wliq2)
     h = (cpice*wice+cpliq*wliq) * (t-tfrz)+hfus*wliq
     h2= (cpice*wice2+cpliq*wliq2) * (t2-tfrz)+hfus*wliq2

     hc = h + h2
     tc = tfrz + (hc - hfus*wliqc) / (cpice*wicec + cpliq*wliqc)

     dz = dzc
     wice = wicec
     wliq = wliqc
     t = tc

   end subroutine Combo
   !-----------------------------------------------------------------------
   
   function MassWeightedSnowRadius( rds1, rds2, swtot, zwtot ) result(mass_weighted_snowradius)
     ! (from CLMv5)
     ! !DESCRIPTION:
     ! Calculate the mass weighted snow radius when two layers are combined
     !$acc routine seq 
     ! !USES:
     use AerosolMod   , only : snw_rds_min
     use SnowSnicarMod, only : snw_rds_max
     implicit none
     ! !ARGUMENTS:
     real(r8), intent(IN) :: rds1         ! Layer 1 radius
     real(r8), intent(IN) :: rds2         ! Layer 2 radius
     real(r8), intent(IN) :: swtot        ! snow water total layer 2
     real(r8), intent(IN) :: zwtot        ! snow water total layer 1
     real(r8) :: mass_weighted_snowradius ! resulting bounded mass weighted snow radius

     mass_weighted_snowradius = (rds2*swtot + rds1*zwtot)/(swtot+zwtot)

     if (      mass_weighted_snowradius > snw_rds_max ) then
        mass_weighted_snowradius = snw_rds_max
     else if ( mass_weighted_snowradius < snw_rds_min ) then
        mass_weighted_snowradius = snw_rds_min
     end if
   end function MassWeightedSnowRadius
   !-----------------------------------------------------------------------
   
   subroutine BuildSnowFilter(num_nolakec, filter_nolakec, &
        num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)
     !
     ! !DESCRIPTION:
     ! Constructs snow filter for use in vectorized loops for snow hydrology.
     !
     ! !USES:
     !
     ! !ARGUMENTS:
     integer           , intent(in)  :: num_nolakec       ! number of column non-lake points in column filter
     integer           , intent(in)  :: filter_nolakec(:) ! column filter for non-lake points
     integer           , intent(out) :: num_snowc         ! number of column snow points in column filter
     integer           , intent(out) :: filter_snowc(:)   ! column filter for snow points
     integer           , intent(out) :: num_nosnowc       ! number of column non-snow points in column filter
     integer           , intent(out) :: filter_nosnowc(:) ! column filter for non-snow points
     !
     ! !LOCAL VARIABLES:
     integer  :: fc, c, fsnow, fnosnow 
     integer :: snow_tot, nosnow_tot
     !-----------------------------------------------------------------------

     ! Build snow/no-snow filters for other subroutines
     !$acc enter data create(fsnow, fnosnow)

     snow_tot = 0
     nosnow_tot = 0
     !$acc parallel loop independent gang vector present(filter_snowc(:),filter_nosnowc(:),filter_nolakec(:)) &
     !$acc private(fsnow,fnosnow) copy(snow_tot, nosnow_tot)
     do fc = 1, num_nolakec
      
        c = filter_nolakec(fc)
        if (col_pp%snl(c) < 0) then
            !$acc atomic capture 
            snow_tot = snow_tot + 1
            fsnow = snow_tot 
            !$acc end atomic 
            filter_snowc(fsnow) = c
        else
            !$acc atomic capture 
            nosnow_tot = nosnow_tot + 1
            fnosnow = nosnow_tot 
            !$acc end atomic 
            filter_nosnowc(fnosnow) = c
        end if
     end do

     num_snowc = snow_tot 
     num_nosnowc = nosnow_tot

     !$acc exit data delete(fsnow, fnosnow)

   end subroutine BuildSnowFilter

end module SnowHydrologyMod
