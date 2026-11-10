! SNICAR with the bulk (non-modal) aerosol treatment, for compsets without BC (no -DMODAL_AER).
! Host only: there are no OpenACC directives here and SNICAR_AD_RT_nonmodal stops if the GPU is enabled.
! Callers select it in their use statements; removing the non-modal path means deleting this file,
! those #ifdef blocks and the use_modal_aer branches.
module SnowSnicarNonModalMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate albedo of snow containing impurities
  ! and the evolution of snow effective radius
  !
  ! !USES:
  use elm_varctl, only : elm_device_ready  ! compute constructs run on the host before the device copy
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_flush
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use elm_varctl      , only : iulog, use_firn_percolation_and_compaction
  use elm_varcon      , only : namec 
  use shr_const_mod   , only : SHR_CONST_RHOICE
  use abortutils      , only : endrun
  use decompMod       , only : bounds_type
  use AerosolMod      , only : snw_rds_min
  use AerosolType     , only : aerosol_type
  use GridcellType    , only : grc_pp
  use LandunitType    , only : lun_pp
  use ColumnType      , only : col_pp
  use ColumnDataType  , only : col_es, col_ws, col_wf
  !
  use timeinfoMod

  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SnowAge_grain_nonmodal    ! Snow effective grain size evolution
  public :: SnowAge_init_nonmodal     ! Initial read in of snow-aging file
  public :: SnowOptics_init_nonmodal  ! Initial read in of snow-optics file
  public :: SNICAR_AD_RT_nonmodal ! Snow albedo and vertically-resolved solar absorption
  !
  ! !PUBLIC DATA MEMBERS:
  integer, parameter :: sno_nbr_aer =   8        ! number of aerosol species in snowpack
                                                          ! (indices described above) [nbr]
  logical, parameter :: DO_SNO_OC =    .false.   ! parameter to include organic carbon (OC)
                                                          ! in snowpack radiative calculations
  logical, parameter :: DO_SNO_AER =   .true.    ! parameter to include aerosols in snowpack radiative calculations
  ! !PRIVATE DATA MEMBERS:
  integer,  parameter :: numrad_snw  =   5               ! number of spectral bands used in snow model [nbr]
  integer,  parameter :: nir_bnd_bgn =   2               ! first band index in near-IR spectrum [idx]
  integer,  parameter :: nir_bnd_end =   5               ! ending near-IR band index [idx]
  integer,  parameter :: idx_Mie_snw_mx = 1471           ! number of effective radius indices used in Mie lookup table [idx]
  integer,  parameter :: idx_T_max      = 11             ! maxiumum temperature index used in aging lookup table [idx]
  integer,  parameter :: idx_T_min      = 1              ! minimum temperature index used in aging lookup table [idx]
  integer,  parameter :: idx_Tgrd_max   = 31             ! maxiumum temperature gradient index used in aging lookup table [idx]
  integer,  parameter :: idx_Tgrd_min   = 1              ! minimum temperature gradient index used in aging lookup table [idx]
  integer,  parameter :: idx_rhos_max   = 8              ! maxiumum snow density index used in aging lookup table [idx]
  integer,  parameter :: idx_rhos_min   = 1              ! minimum snow density index used in aging lookup table [idx]

  !mgf++
  !mgf--

  integer,  parameter :: snw_rds_max_tbl = 1500          ! maximum effective radius defined in Mie lookup table [microns]
  integer,  parameter :: snw_rds_min_tbl = 30            ! minimium effective radius defined in Mie lookup table [microns]
  real(r8), parameter :: snw_rds_max     = 1500._r8      ! maximum allowed snow effective radius [microns]
  real(r8)            :: snw_rds_refrz   = 1000._r8      ! effective radius of re-frozen snow [microns]
  real(r8), parameter :: min_snw = 1.0E-30_r8            ! minimum snow mass required for SNICAR RT calculation [kg m-2]
  !real(r8), parameter :: C1_liq_Brun89 = 1.28E-17_r8    ! constant for liquid water grain growth [m3 s-1],
                                                         ! from Brun89
  real(r8), parameter :: C1_liq_Brun89 = 0._r8           ! constant for liquid water grain growth [m3 s-1],
                                                         ! from Brun89: zeroed to accomodate dry snow aging
  real(r8), parameter :: C2_liq_Brun89 = 4.22E-13_r8     ! constant for liquid water grain growth [m3 s-1],
                                                         ! from Brun89: corrected for LWC in units of percent

  real(r8), parameter :: tim_cns_bc_rmv  = 2.2E-8_r8     ! time constant for removal of BC in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)
  real(r8), parameter :: tim_cns_oc_rmv  = 2.2E-8_r8     ! time constant for removal of OC in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)
  real(r8), parameter :: tim_cns_dst_rmv = 2.2E-8_r8     ! time constant for removal of dust in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)

  ! scaling of the snow aging rate (tuning option):
  logical, parameter :: flg_snoage_scl    = .false.                 ! flag for scaling the snow aging rate by some arbitrary factor
  real(r8), parameter :: xdrdt = 1.0_r8                  ! arbitrary factor applied to snow aging rate
  ! snow and aerosol Mie parameters:
  ! (arrays declared here, but are set in iniTimeConst)
  ! (idx_Mie_snw_mx is number of snow radii with defined parameters (i.e. from 30um to 1500um))

  ! direct-beam weighted ice optical properties
  real(r8) :: ss_alb_snw_drc     (idx_Mie_snw_mx,numrad_snw);
  real(r8) :: asm_prm_snw_drc    (idx_Mie_snw_mx,numrad_snw);
  real(r8) :: ext_cff_mss_snw_drc(idx_Mie_snw_mx,numrad_snw);

  ! diffuse radiation weighted ice optical properties
  real(r8) :: ss_alb_snw_dfs     (idx_Mie_snw_mx,numrad_snw);
  real(r8) :: asm_prm_snw_dfs    (idx_Mie_snw_mx,numrad_snw);
  real(r8) :: ext_cff_mss_snw_dfs(idx_Mie_snw_mx,numrad_snw);

  ! direct & diffuse flux
  real(r8) :: flx_wgt_dir    (6, 90,numrad_snw) ! direct flux, six atmospheric types, 0-89 SZA
  real(r8) :: flx_wgt_dif    (6, numrad_snw)    ! diffuse flux, six atmospheric types
  
  ! snow grain shape
  integer, parameter :: snow_shape_sphere          = 1
  integer, parameter :: snow_shape_spheroid        = 2
  integer, parameter :: snow_shape_hexagonal_plate = 3
  integer, parameter :: snow_shape_koch_snowflake  = 4
  
  ! atmospheric condition for SNICAR-AD
  integer, parameter :: atm_type_default             = 0
  integer, parameter :: atm_type_mid_latitude_winter = 1
  integer, parameter :: atm_type_mid_latitude_summer = 2
  integer, parameter :: atm_type_sub_Arctic_winter   = 3
  integer, parameter :: atm_type_sub_Arctic_summer   = 4
  integer, parameter :: atm_type_summit_Greenland    = 5
  integer, parameter :: atm_type_high_mountain       = 6
  

  
  !mgf++
  ! Size-dependent BC optical properties. Currently a fixed BC size is
  ! assumed, but this framework enables optical properties to be
  ! assigned based on the BC effective radius, should this be
  ! implemented in the future.
  !
  ! within-ice BC (i.e., BC that was deposited within hydrometeors)
  real(r8) :: ss_alb_bc1     (numrad_snw);
  real(r8) :: asm_prm_bc1    (numrad_snw);
  real(r8) :: ext_cff_mss_bc1(numrad_snw);

  ! external BC
  real(r8) :: ss_alb_bc2     (numrad_snw);
  real(r8) :: asm_prm_bc2    (numrad_snw);
  real(r8) :: ext_cff_mss_bc2(numrad_snw);
  !mgf--

!  ! hydrophiliic BC
!  real(r8) :: ss_alb_bc1(numrad_snw)
!  real(r8) :: asm_prm_bc1(numrad_snw)
!  real(r8) :: ext_cff_mss_bc1(numrad_snw)
!
!  ! hydrophobic BC
!  real(r8) :: ss_alb_bc2(numrad_snw)
!  real(r8) :: asm_prm_bc2(numrad_snw)
!  real(r8) :: ext_cff_mss_bc2(numrad_snw)

  ! hydrophobic OC
  real(r8) :: ss_alb_oc1     (numrad_snw);
  real(r8) :: asm_prm_oc1    (numrad_snw);
  real(r8) :: ext_cff_mss_oc1(numrad_snw);

  ! hydrophilic OC
  real(r8) :: ss_alb_oc2     (numrad_snw);
  real(r8) :: asm_prm_oc2    (numrad_snw);
  real(r8) :: ext_cff_mss_oc2(numrad_snw);

  ! dust species 1:
  real(r8) :: ss_alb_dst1     (numrad_snw);
  real(r8) :: asm_prm_dst1    (numrad_snw);
  real(r8) :: ext_cff_mss_dst1(numrad_snw);

  ! dust species 2:
  real(r8) :: ss_alb_dst2     (numrad_snw);
  real(r8) :: asm_prm_dst2    (numrad_snw);
  real(r8) :: ext_cff_mss_dst2(numrad_snw);

  ! dust species 3:
  real(r8) :: ss_alb_dst3     (numrad_snw);
  real(r8) :: asm_prm_dst3    (numrad_snw);
  real(r8) :: ext_cff_mss_dst3(numrad_snw);

  ! dust species 4:
  real(r8) :: ss_alb_dst4     (numrad_snw);
  real(r8) :: asm_prm_dst4    (numrad_snw);
  real(r8) :: ext_cff_mss_dst4(numrad_snw);


  !mgf++
  ! Absorption enhancement factors for within-ice BC
  !mgf--

  ! best-fit parameters for snow aging defined over:
  !  11 temperatures from 225 to 273 K
  !  31 temperature gradients from 0 to 300 K/m
  !   8 snow densities from 0 to 350 kg/m3
  ! (arrays declared here, but are set in iniTimeConst)
  real(r8), pointer :: snowage_tau(:,:,:)   => null()! (idx_rhos_max,idx_Tgrd_max,idx_T_max)
  real(r8), pointer :: snowage_kappa(:,:,:) => null()! (idx_rhos_max,idx_Tgrd_max,idx_T_max)
  real(r8), pointer :: snowage_drdt0(:,:,:) => null()! idx_rhos_max,idx_Tgrd_max,idx_T_max)
  !
  ! !REVISION HISTORY:
  ! Created by Mark Flanner
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine SnowAge_grain_nonmodal(bounds, &
       num_snowc, filter_snowc, num_nosnowc, filter_nosnowc)
    !
    ! !DESCRIPTION:
    ! Updates the snow effective grain size (radius).
    ! Contributions to grain size evolution are from:
    !   1. vapor redistribution (dry snow)
    !   2. liquid water redistribution (wet snow)
    !   3. re-freezing of liquid water
    !
    ! Vapor redistribution: Method is to retrieve 3 best-bit parameters that
    ! depend on snow temperature, temperature gradient, and density,
    ! that are derived from the microphysical model described in:
    ! Flanner and Zender (2006), Linking snowpack microphysics and albedo
    ! evolution, J. Geophys. Res., 111, D12208, doi:10.1029/2005JD006834.
    ! The parametric equation has the form:
    ! dr/dt = drdt_0*(tau/(dr_fresh+tau))^(1/kappa), where:
    !   r is the effective radius,
    !   tau and kappa are best-fit parameters,
    !   drdt_0 is the initial rate of change of effective radius, and
    !   dr_fresh is the difference between the current and fresh snow states
    !  (r_current - r_fresh).
    !
    ! Liquid water redistribution: Apply the grain growth function from:
    !   Brun, E. (1989), Investigation of wet-snow metamorphism in respect of
    !   liquid-water content, Annals of Glaciology, 13, 22-26.
    !   There are two parameters that describe the grain growth rate as
    !   a function of snow liquid water content (LWC). The "LWC=0" parameter
    !   is zeroed here because we are accounting for dry snowing with a
    !   different representation
    !
    ! Re-freezing of liquid water: Assume that re-frozen liquid water clumps
    !   into an arbitrarily large effective grain size (snw_rds_refrz).
    !   The phenomenon is observed (Grenfell), but so far unquantified, as far as
    !   I am aware.
    !
    ! !USES:
    use elm_varpar       , only : nlevsno, nlevsno_max
    use elm_varcon       , only : spval
    use shr_const_mod    , only : SHR_CONST_RHOICE, SHR_CONST_PI
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_snowc         ! number of column snow points in column filter
    integer                , intent(in)    :: filter_snowc(:)   ! column filter for snow points
    integer                , intent(in)    :: num_nosnowc       ! number of column non-snow points in column filter
    integer                , intent(in)    :: filter_nosnowc(:) ! column filter for non-snow points

    !
    ! !LOCAL VARIABLES:
    integer :: snl_top                      ! top snow layer index [idx]
    integer :: snl_btm                      ! bottom snow layer index [idx]
    integer :: i                            ! layer index [idx]
    integer :: c_idx                        ! column index [idx]
    integer :: fc                           ! snow column filter index [idx]
    integer :: T_idx                        ! snow aging lookup table temperature index [idx]
    integer :: Tgrd_idx                     ! snow aging lookup table temperature gradient index [idx]
    integer :: rhos_idx                     ! snow aging lookup table snow density index [idx]
    real(r8) :: t_snotop                    ! temperature at upper layer boundary [K]
    real(r8) :: t_snobtm                    ! temperature at lower layer boundary [K]
    real(r8) :: dTdz(bounds%begc:bounds%endc,-nlevsno:0)    ! snow temperature gradient (col,lyr) [K m-1]
    real(r8) :: bst_tau                     ! snow aging parameter retrieved from lookup table [hour]
    real(r8) :: bst_kappa                   ! snow aging parameter retrieved from lookup table [unitless]
    real(r8) :: bst_drdt0                   ! snow aging parameter retrieved from lookup table [um hr-1]
    real(r8) :: dr                          ! incremental change in snow effective radius [um]
    real(r8) :: dr_wet                      ! incremental change in snow effective radius from wet growth [um]
    real(r8) :: dr_fresh                    ! difference between fresh snow r_e and current r_e [um]
    real(r8) :: newsnow                     ! fresh snowfall [kg m-2]
    real(r8) :: refrzsnow                   ! re-frozen snow [kg m-2]
    real(r8) :: frc_newsnow                 ! fraction of layer mass that is new snow [frc]
    real(r8) :: frc_oldsnow                 ! fraction of layer mass that is old snow [frc]
    real(r8) :: frc_refrz                   ! fraction of layer mass that is re-frozen snow [frc]
    real(r8) :: frc_liq                     ! fraction of layer mass that is liquid water[frc]
    real(r8)  :: dtime                       ! land model time step [sec]
    real(r8) :: rhos                        ! snow density [kg m-3]
    real(r8) :: h2osno_lyr                  ! liquid + solid H2O in snow layer [kg m-2]
    real(r8) :: cdz(-nlevsno_max+1:0)       ! column average layer thickness [m] (private per thread: fixed size)
    !--------------------------------------------------------------------------!

    associate(                                                      &
         snl                => col_pp%snl                            , & ! Input:  [integer  (:)   ]  negative number of snow layers (col) [nbr]
         dz                 => col_pp%dz                             , & ! Input:  [real(r8) (:,:) ]  layer thickness (col,lyr) [m]

         qflx_snow_grnd_col => col_wf%qflx_snow_grnd  , & ! Input:  [real(r8) (:)   ]  snow on ground after interception (col) [kg m-2 s-1]
         qflx_snwcp_ice     => col_wf%qflx_snwcp_ice  , & ! Input:  [real(r8) (:)   ]  excess precipitation due to snow capping [kg m-2 s-1]
         qflx_snofrz_lyr    => col_wf%qflx_snofrz_lyr , & ! Input:  [real(r8) (:,:) ]  snow freezing rate (col,lyr) [kg m-2 s-1]

         do_capsnow         => col_ws%do_capsnow     , & ! Input:  [logical  (:)   ]  true => do snow capping
         frac_sno           => col_ws%frac_sno_eff   , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
         h2osno             => col_ws%h2osno         , & ! Input:  [real(r8) (:)   ]  snow water (col) [mm H2O]
         h2osoi_liq         => col_ws%h2osoi_liq     , & ! Input:  [real(r8) (:,:) ]  liquid water content (col,lyr) [kg m-2]
         h2osoi_ice         => col_ws%h2osoi_ice     , & ! Input:  [real(r8) (:,:) ]  ice content (col,lyr) [kg m-2]
         snw_rds            => col_ws%snw_rds        , & ! Output: [real(r8) (:,:) ]  effective grain radius (col,lyr) [microns, m-6]
         snw_rds_top        => col_ws%snw_rds_top    , & ! Output: [real(r8) (:)   ]  effective grain radius, top layer (col) [microns, m-6]
         sno_liq_top        => col_ws%sno_liq_top    , & ! Output: [real(r8) (:)   ]  liquid water fraction (mass) in top snow layer (col) [frc]

         t_soisno           => col_es%t_soisno      , & ! Input:  [real(r8) (:,:) ]  soil and snow temperature (col,lyr) [K]
         t_grnd             => col_es%t_grnd        , & ! Input:  [real(r8) (:)   ]  ground temperature (col) [K]
         snot_top           => col_es%snot_top      , & ! Output: [real(r8) (:)   ]  temperature in top snow layer (col) [K]
         dTdz_top           => col_es%dTdz_top        & ! Output: [real(r8) (:)   ]  temperature gradient in top layer (col) [K m-1]
         )


      ! set timestep and step interval
      dtime = dtime_mod

      if (use_firn_percolation_and_compaction) then
          snw_rds_refrz   = 1500._r8
      else
          snw_rds_refrz   = 1000._r8
      endif

      ! loop over columns that have at least one snow layer (each writes its own column)
      do fc = 1, num_snowc
         c_idx = filter_snowc(fc)

         snl_btm = 0
         snl_top = snl(c_idx) + 1

         cdz(snl_top:snl_btm)=frac_sno(c_idx)*dz(c_idx,snl_top:snl_btm)

         ! loop over snow layers
         do i=snl_top,snl_btm,1
            !
            !**********  1. DRY SNOW AGING  ***********
            !
            h2osno_lyr = h2osoi_liq(c_idx,i) + h2osoi_ice(c_idx,i)

            ! temperature gradient
            if (i == snl_top) then
               ! top layer
               t_snotop = t_soisno(c_idx,snl_top)
               t_snobtm = (t_soisno(c_idx,i+1)*dz(c_idx,i) &
                    + t_soisno(c_idx,i)*dz(c_idx,i+1)) &
                    / (dz(c_idx,i)+dz(c_idx,i+1))
            else
               t_snotop = (t_soisno(c_idx,i-1)*dz(c_idx,i) &
                    + t_soisno(c_idx,i)*dz(c_idx,i-1)) &
                    / (dz(c_idx,i)+dz(c_idx,i-1))
               t_snobtm = (t_soisno(c_idx,i+1)*dz(c_idx,i) &
                    + t_soisno(c_idx,i)*dz(c_idx,i+1)) &
                    / (dz(c_idx,i)+dz(c_idx,i+1))
            endif

            dTdz(c_idx,i) = abs((t_snotop - t_snobtm) / cdz(i))

            ! snow density
            rhos = (h2osoi_liq(c_idx,i)+h2osoi_ice(c_idx,i)) / cdz(i)

            ! make sure rhos doesn't drop below 50 (see rhos_idx below)
            rhos=max(50._r8,rhos)

            ! best-fit table indecies
            T_idx    = nint((t_soisno(c_idx,i)-223) / 5) + 1
            Tgrd_idx = nint(dTdz(c_idx,i) / 10) + 1
            rhos_idx = nint((rhos-50) / 50) + 1

            ! boundary check:
            if (T_idx < idx_T_min) then
               T_idx = idx_T_min
            endif
            if (T_idx > idx_T_max) then
               T_idx = idx_T_max
            endif
            if (Tgrd_idx < idx_Tgrd_min) then
               Tgrd_idx = idx_Tgrd_min
            endif
            if (Tgrd_idx > idx_Tgrd_max) then
               Tgrd_idx = idx_Tgrd_max
            endif
            if (rhos_idx < idx_rhos_min) then
               rhos_idx = idx_rhos_min
            endif
            if (rhos_idx > idx_rhos_max) then
               rhos_idx = idx_rhos_max
            endif

            ! best-fit parameters
            bst_tau   = snowage_tau(rhos_idx,Tgrd_idx,T_idx)
            bst_kappa = snowage_kappa(rhos_idx,Tgrd_idx,T_idx)
            bst_drdt0 = snowage_drdt0(rhos_idx,Tgrd_idx,T_idx)


            ! change in snow effective radius, using best-fit parameters
            ! added checks suggested by mgf. --HW 10/15/2015
            dr_fresh = snw_rds(c_idx,i)-snw_rds_min

            dr = (bst_drdt0*(bst_tau/(dr_fresh+bst_tau))**(1/bst_kappa)) * (dtime/3600)
            !
            !**********  2. WET SNOW AGING  ***********
            !
            ! We are assuming wet and dry evolution occur simultaneously, and
            ! the contributions from both can be summed.
            ! This is justified by setting the linear offset constant C1_liq_Brun89 to zero [Brun, 1989]

            ! liquid water faction
            frc_liq = min(0.1_r8, (h2osoi_liq(c_idx,i) / (h2osoi_liq(c_idx,i)+h2osoi_ice(c_idx,i))))

            !dr_wet = 1E6_r8*(dtime*(C1_liq_Brun89 + C2_liq_Brun89*(frc_liq**(3))) / (4*SHR_CONST_PI*(snw_rds(c_idx,i)/1E6)**(2)))
            !simplified, units of microns:
            dr_wet = 1E18_r8*(dtime*(C2_liq_Brun89*(frc_liq**(3))) / (4*SHR_CONST_PI*snw_rds(c_idx,i)**(2)))

            dr = dr + dr_wet

            !
            !**********  3. SNOWAGE SCALING (TURNED OFF BY DEFAULT)  *************
            !
            ! Multiply rate of change of effective radius by some constant, xdrdt
            if (flg_snoage_scl) then
               dr = dr*xdrdt
            endif


            !
            !**********  4. INCREMENT EFFECTIVE RADIUS, ACCOUNTING FOR:  ***********
            !               DRY AGING
            !               WET AGING
            !               FRESH SNOW
            !               RE-FREEZING
            !
            ! new snowfall [kg/m2]
            if (do_capsnow(c_idx) .and. .not. use_firn_percolation_and_compaction) then
               newsnow = max(0._r8, (qflx_snwcp_ice(c_idx)*dtime))
            else
               newsnow = max(0._r8, (qflx_snow_grnd_col(c_idx)*dtime))
            endif

            ! snow that has re-frozen [kg/m2]
            refrzsnow = max(0._r8, (qflx_snofrz_lyr(c_idx,i)*dtime))

            ! fraction of layer mass that is re-frozen
            frc_refrz = refrzsnow / h2osno_lyr

            ! fraction of layer mass that is new snow
            if (i == snl_top) then
               frc_newsnow = newsnow / h2osno_lyr
            else
               frc_newsnow = 0._r8
            endif

            if ((frc_refrz + frc_newsnow) > 1._r8) then
               frc_refrz = frc_refrz / (frc_refrz + frc_newsnow)
               frc_newsnow = 1._r8 - frc_refrz
               frc_oldsnow = 0._r8
            else
               frc_oldsnow = 1._r8 - frc_refrz - frc_newsnow
            endif

            ! mass-weighted mean of fresh snow, old snow, and re-frozen snow effective radius
            snw_rds(c_idx,i) = (snw_rds(c_idx,i)+dr)*frc_oldsnow + snw_rds_min*frc_newsnow + snw_rds_refrz*frc_refrz
            !
            !**********  5. CHECK BOUNDARIES   ***********
            !
            ! boundary check
            if (snw_rds(c_idx,i) < snw_rds_min) then
               snw_rds(c_idx,i) = snw_rds_min
            endif

            if (snw_rds(c_idx,i) > snw_rds_max) then
               snw_rds(c_idx,i) = snw_rds_max
            end if

            ! set top layer variables for history files
            if (i == snl_top) then
               snot_top(c_idx)    = t_soisno(c_idx,i)
               dTdz_top(c_idx)    = dTdz(c_idx,i)
               snw_rds_top(c_idx) = snw_rds(c_idx,i)
               sno_liq_top(c_idx) = h2osoi_liq(c_idx,i) / (h2osoi_liq(c_idx,i)+h2osoi_ice(c_idx,i))
            endif

         enddo
      enddo

      ! Special case: snow on ground, but not enough to have defined a snow layer:
      !   set snw_rds to fresh snow grain size:
      do fc = 1, num_nosnowc
         c_idx = filter_nosnowc(fc)
         if (h2osno(c_idx) > 0._r8) then
            snw_rds(c_idx,0) = snw_rds_min
         endif
      enddo

    end associate

  end subroutine SnowAge_grain_nonmodal

  !-----------------------------------------------------------------------
     subroutine SnowOptics_init_nonmodal( )

      use fileutils  , only : getfil
      use elm_varctl , only : fsnowoptics, snicar_atm_type
      use spmdMod    , only : masterproc
      use ncdio_pio  , only : file_desc_t, ncd_io, ncd_pio_openfile, ncd_pio_closefile
      use ncdio_pio  , only : ncd_pio_openfile, ncd_inqfdims, ncd_pio_closefile, ncd_inqdid, ncd_inqdlen

      type(file_desc_t)  :: ncid                        ! netCDF file id
      character(len=256) :: locfn                       ! local filename
      character(len= 32) :: subname = 'SnowOptics_init_nonmodal' ! subroutine name
      integer            :: ier                         ! error status
      integer            :: atm_type_index              ! index for atmospheric type
	  
     !mgf++
     logical :: readvar      ! determine if variable was read from NetCDF file
     !mgf--

      atm_type_index = atm_type_default
      ! Define atmospheric type
      if (trim(snicar_atm_type) == 'default') then
        atm_type_index = atm_type_default
      elseif (trim(snicar_atm_type) == 'mid-latitude_winter') then
        atm_type_index = atm_type_mid_latitude_winter
      elseif (trim(snicar_atm_type) == 'mid-latitude_summer') then
        atm_type_index = atm_type_mid_latitude_summer
      elseif (trim(snicar_atm_type) == 'sub-Arctic_winter') then
        atm_type_index = atm_type_sub_Arctic_winter
      elseif (trim(snicar_atm_type) == 'sub-Arctic_summer') then
        atm_type_index = atm_type_sub_Arctic_summer
      elseif (trim(snicar_atm_type) == 'summit_Greenland') then
        atm_type_index = atm_type_summit_Greenland
      elseif (trim(snicar_atm_type) == 'high_mountain') then
        atm_type_index = atm_type_high_mountain
      else
	write(iulog,*) "snicar_atm_type = ", snicar_atm_type
        call endrun( "snicar_atm_type is unknown" )
      endif
	  
      !
      ! Open optics file:
      if(masterproc) write(iulog,*) 'Attempting to read snow optical properties .....'
      call getfil (fsnowoptics, locfn, 0)
      call ncd_pio_openfile(ncid, locfn, 0)
      if(masterproc) write(iulog,*) subname,trim(fsnowoptics)

      ! direct-beam snow Mie parameters:
      call ncd_io('ss_alb_ice_drc', ss_alb_snw_drc,            'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'asm_prm_ice_drc',asm_prm_snw_drc,          'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'ext_cff_mss_ice_drc', ext_cff_mss_snw_drc, 'read', ncid, posNOTonfile=.true.)
      !
      ! diffuse snow Mie parameters
      call ncd_io( 'ss_alb_ice_dfs', ss_alb_snw_dfs,           'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'asm_prm_ice_dfs', asm_prm_snw_dfs,         'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'ext_cff_mss_ice_dfs', ext_cff_mss_snw_dfs, 'read', ncid, posNOTonfile=.true.)
	  
      !!! Direct and diffuse flux under different atmospheric conditions
      ! Direct-beam incident spectral flux: 
      call ncd_io( 'flx_wgt_dir', flx_wgt_dir,           'read', ncid, readvar=readvar, posNOTonfile=.true.)
      if ((atm_type_index > 0) .and. (.not. readvar)) call endrun( "ERROR: error in reading in flx_wgt_dir data" )
      ! Diffuse incident spectral flux:
      call ncd_io( 'flx_wgt_dif', flx_wgt_dif,           'read', ncid, readvar=readvar, posNOTonfile=.true.)
      if ((atm_type_index > 0) .and. (.not. readvar)) call endrun( "ERROR: error in reading in flx_wgt_dif data" )

      !
     ! bulk aerosol treatment
      ! BC species 1 Mie parameters
      call ncd_io( 'ss_alb_bcphil', ss_alb_bc1,           'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'asm_prm_bcphil', asm_prm_bc1,         'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'ext_cff_mss_bcphil', ext_cff_mss_bc1, 'read', ncid, posNOTonfile=.true.)
      !
      ! BC species 2 Mie parameters
      call ncd_io( 'ss_alb_bcphob', ss_alb_bc2,           'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'asm_prm_bcphob', asm_prm_bc2,         'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'ext_cff_mss_bcphob', ext_cff_mss_bc2, 'read', ncid, posNOTonfile=.true.)
      !
      ! OC species 1 Mie parameters
      call ncd_io( 'ss_alb_ocphil',      ss_alb_oc1,      'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'asm_prm_ocphil',     asm_prm_oc1,     'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'ext_cff_mss_ocphil', ext_cff_mss_oc1, 'read', ncid, posNOTonfile=.true.)
      !
      ! OC species 2 Mie parameters
      call ncd_io( 'ss_alb_ocphob', ss_alb_oc2,           'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'asm_prm_ocphob', asm_prm_oc2,         'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'ext_cff_mss_ocphob', ext_cff_mss_oc2, 'read', ncid, posNOTonfile=.true.)
      !
      ! dust species 1 Mie parameters
      call ncd_io( 'ss_alb_dust01', ss_alb_dst1,           'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'asm_prm_dust01', asm_prm_dst1,         'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'ext_cff_mss_dust01', ext_cff_mss_dst1, 'read', ncid, posNOTonfile=.true.)
      !
      ! dust species 2 Mie parameters
      call ncd_io( 'ss_alb_dust02', ss_alb_dst2,           'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'asm_prm_dust02', asm_prm_dst2,         'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'ext_cff_mss_dust02', ext_cff_mss_dst2, 'read', ncid, posNOTonfile=.true.)
      !
      ! dust species 3 Mie parameters
      call ncd_io( 'ss_alb_dust03', ss_alb_dst3,           'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'asm_prm_dust03', asm_prm_dst3,         'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'ext_cff_mss_dust03', ext_cff_mss_dst3, 'read', ncid, posNOTonfile=.true.)
      !
      ! dust species 4 Mie parameters
      call ncd_io( 'ss_alb_dust04', ss_alb_dst4,           'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'asm_prm_dust04', asm_prm_dst4,         'read', ncid, posNOTonfile=.true.)
      call ncd_io( 'ext_cff_mss_dust04', ext_cff_mss_dst4, 'read', ncid, posNOTonfile=.true.)
      !
      !

      call ncd_pio_closefile(ncid)
      if (masterproc) then
        !
         write(iulog,*) 'Successfully read snow optical properties'
         ! print some diagnostics:
         write (iulog,*) 'SNICAR: Mie single scatter albedos for direct-beam ice, rds=100um: ', &
              ss_alb_snw_drc(71,1), ss_alb_snw_drc(71,2), ss_alb_snw_drc(71,3),     &
              ss_alb_snw_drc(71,4), ss_alb_snw_drc(71,5)
         write (iulog,*) 'SNICAR: Mie single scatter albedos for diffuse ice, rds=100um: ',     &
              ss_alb_snw_dfs(71,1), ss_alb_snw_dfs(71,2), ss_alb_snw_dfs(71,3),     &
              ss_alb_snw_dfs(71,4), ss_alb_snw_dfs(71,5)
         if (DO_SNO_OC) then
            write (iulog,*) 'SNICAR: Including OC aerosols from snow radiative transfer calculations'
         else
            write (iulog,*) 'SNICAR: Excluding OC aerosols from snow radiative transfer calculations'
         endif
         !
         write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophillic BC: ', &
              ss_alb_bc1(1), ss_alb_bc1(2), ss_alb_bc1(3), ss_alb_bc1(4), ss_alb_bc1(5)
         write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophobic BC: ', &
              ss_alb_bc2(1), ss_alb_bc2(2), ss_alb_bc2(3), ss_alb_bc2(4), ss_alb_bc2(5)
         if (DO_SNO_OC) then
            write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophillic OC: ', &
                 ss_alb_oc1(1), ss_alb_oc1(2), ss_alb_oc1(3), ss_alb_oc1(4), ss_alb_oc1(5)
            write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophobic OC: ', &
                 ss_alb_oc2(1), ss_alb_oc2(2), ss_alb_oc2(3), ss_alb_oc2(4), ss_alb_oc2(5)
         endif
         write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 1: ', &
              ss_alb_dst1(1), ss_alb_dst1(2), ss_alb_dst1(3), ss_alb_dst1(4), ss_alb_dst1(5)
         write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 2: ', &
              ss_alb_dst2(1), ss_alb_dst2(2), ss_alb_dst2(3), ss_alb_dst2(4), ss_alb_dst2(5)
         write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 3: ', &
              ss_alb_dst3(1), ss_alb_dst3(2), ss_alb_dst3(3), ss_alb_dst3(4), ss_alb_dst3(5)
         write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 4: ', &
              ss_alb_dst4(1), ss_alb_dst4(2), ss_alb_dst4(3), ss_alb_dst4(4), ss_alb_dst4(5)
         write(iulog,*)
      end if
      !
    end subroutine SnowOptics_init_nonmodal

   !-----------------------------------------------------------------------
   subroutine SnowAge_init_nonmodal( )
     use ELM_varctl      , only : fsnowaging
     use fileutils       , only : getfil
     use spmdMod         , only : masterproc
     use ncdio_pio       , only : file_desc_t, ncd_io, ncd_pio_openfile, ncd_pio_closefile

     type(file_desc_t)  :: ncid                        ! netCDF file id
     character(len=256) :: locfn                       ! local filename
     character(len= 32) :: subname = 'SnowOptics_init_nonmodal' ! subroutine name
     integer            :: varid                       ! netCDF id's
     integer            :: ier                         ! error status
     !
     ! Open snow aging (effective radius evolution) file:
     allocate(snowage_tau(idx_rhos_max,idx_Tgrd_max,idx_T_max))
     allocate(snowage_kappa(idx_rhos_max,idx_Tgrd_max,idx_T_max))
     allocate(snowage_drdt0(idx_rhos_max,idx_Tgrd_max,idx_T_max))
     !
     if(masterproc)  write(iulog,*) 'Attempting to read snow aging parameters .....'
     call getfil (fsnowaging, locfn, 0)
     call ncd_pio_openfile(ncid, locfn, 0)
     if(masterproc) write(iulog,*) subname,trim(fsnowaging)
     !
     ! snow aging parameters
     !
     call ncd_io('tau', snowage_tau,       'read', ncid, posNOTonfile=.true.)
     call ncd_io('kappa', snowage_kappa,   'read', ncid, posNOTonfile=.true.)
     call ncd_io('drdsdt0', snowage_drdt0, 'read', ncid, posNOTonfile=.true.)
     !
     call ncd_pio_closefile(ncid)
     if (masterproc) then
       !
        write(iulog,*) 'Successfully read snow aging properties'
        !
        ! print some diagnostics:
        write (iulog,*) 'SNICAR: snowage tau for T=263K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_tau(3,11,9)
        write (iulog,*) 'SNICAR: snowage kappa for T=263K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_kappa(3,11,9)
        write (iulog,*) 'SNICAR: snowage dr/dt_0 for T=263K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_drdt0(3,11,9)
     endif

    end subroutine SnowAge_init_nonmodal

   !-----------------------------------------------------------------------
   subroutine SNICAR_AD_RT_nonmodal (bounds, num_nourbanc, filter_nourbanc, coszen, ncall, call_slr, call_aer, &
                         nabs, aerosol_vars, albsfc, albout, flx_abs)
     !
     ! !DESCRIPTION:
     ! Determine reflectance of, and vertically-resolved solar absorption in,
     ! snow with impurities, with updated shortwave scheme
     !
     ! The multi-layer solution for multiple-scattering used here is from:
     ! Briegleb, P. and Light, B.: A Delta-Eddington mutiple scattering
     ! parameterization for solar radiation in the sea ice component of the
     ! community climate system model, 2007.
     !
     ! The implementation of the SNICAR-AD model in ELM is described in:
     ! Dang et al., Inter-comparison and improvement of 2-stream shortwave
     ! radiative transfer models for unified treatment of cryospheric surfaces
     ! in ESMs, in review, 2019
     !
     ! To use this subtroutine, set use_snicar_ad = true in ELM
     !
     ! if config_use_snicar_ad = true in MPAS-seaice
     ! Snow on land and snow on sea ice will be treated
     ! with the same model for their solar radiative properties.
     !
     ! The inputs and outputs are the same to subroutine SNICAR_RT
     !
     ! !USES:
     use elm_varpar       , only : nlevsno, nlevsno_max, numrad
     use elm_time_manager , only : get_nstep
     use shr_const_mod    , only : SHR_CONST_PI
     use elm_varctl       , only : snow_shape, snicar_atm_type, use_dust_snow_internal_mixing, use_finetop_rad
     !
     ! !ARGUMENTS:
     ! One call computes `ncall` cases at once (the radiative-forcing variants and the feedback
     ! case of SurfaceAlbedo): case k is direct (call_slr(k) = 1) or diffuse (2) incident flux on
     ! snow holding the aerosol species call_aer(:,k); the first nabs cases also return their
     ! absorbed-flux profiles (flx_abs(:,:,:,k)). Snow layers and aerosol concentrations are
     ! read from the column state; results are indexed by position in the filter.
     type (bounds_type), intent(in)  :: bounds
     integer           , intent(in)  :: num_nourbanc                                       ! number of columns in non-urban filter
     integer           , intent(in)  :: filter_nourbanc(:)                                 ! column filter for non-urban points
     real(r8)          , intent(in)  :: coszen         ( bounds%begc: )                    ! cosine of solar zenith angle for next time step (col) [unitless]
     integer           , intent(in)  :: ncall                                              ! number of cases
     integer           , intent(in)  :: call_slr(:)                                        ! (case) 1 = direct-beam, 2 = diffuse incident flux
     logical           , intent(in)  :: call_aer(:,:)                                      ! (aerosol species, case) species included
     integer           , intent(in)  :: nabs                                               ! cases 1..nabs return flx_abs
     type(aerosol_type), intent(in)  :: aerosol_vars
     real(r8)          , intent(in)  :: albsfc         ( bounds%begc: , 1: )               ! albedo of surface underlying snow (col,bnd) [frc]
     real(r8)          , intent(out) :: albout         ( 1: , 1: , 1: )                    ! snow albedo, 2 bands (=0 if no sun or no snow) (fc,bnd,case) [frc]
     real(r8)          , intent(out) :: flx_abs        ( 1: , -nlevsno+1: , 1: , 1: )      ! absorbed flux per unit incident flux (fc,lyr,bnd,case)
     !
     ! !LOCAL VARIABLES:
     !
     ! variables for snow radiative transfer calculations

     ! Local variables representing single-column values of arrays:
     integer :: snl_lcl                            ! negative number of snow layers [nbr]
     integer :: snw_rds_lcl(-nlevsno_max+1:0)          ! snow effective radius [m^-6]
     real(r8):: flx_slrd_b                         ! direct beam incident irradiance of the band [W/m2] (set to 1)
     real(r8):: flx_slri_b                         ! diffuse incident irradiance of the band [W/m2] (set to 1)
     real(r8):: mss_cnc_aer_lcl(-nlevsno_max+1:0,1:sno_nbr_aer) ! aerosol mass concentration (lyr,aer_nbr) [kg/kg]
     real(r8):: h2osno_lcl                         ! total column snow mass [kg/m2]
     real(r8):: h2osno_liq_lcl(-nlevsno_max+1:0)       ! liquid water mass [kg/m2]
     real(r8):: h2osno_ice_lcl(-nlevsno_max+1:0)       ! ice mass [kg/m2]
     real(r8):: ss_alb_snw_lcl(-nlevsno_max+1:0)       ! single-scatter albedo of ice grains (lyr) [frc]
     real(r8):: asm_prm_snw_lcl(-nlevsno_max+1:0)      ! asymmetry parameter of ice grains (lyr) [frc]
     real(r8):: ext_cff_mss_snw_lcl(-nlevsno_max+1:0)  ! mass extinction coefficient of ice grains (lyr) [m2/kg]
     real(r8):: ss_alb_aer_lcl(sno_nbr_aer)        ! single-scatter albedo of aerosol species (aer_nbr) [frc]
     real(r8):: asm_prm_aer_lcl(sno_nbr_aer)       ! asymmetry parameter of aerosol species (aer_nbr) [frc]
     real(r8):: ext_cff_mss_aer_lcl(sno_nbr_aer)   ! mass extinction coefficient of aerosol species (aer_nbr) [m2/kg]

     !mgf++
     !mgf--


     ! Other local variables
     integer :: DELTA                              ! flag to use Delta approximation (Joseph, 1976)
                                                   ! (1= use, 0= don't use)
     real(r8):: flx_wgt(1:numrad_snw)              ! weights applied to spectral bands,
                                                   ! specific to direct and diffuse cases (bnd) [frc]
     integer :: flg_nosnl                          ! flag: =1 if there is snow, but zero snow layers,
                                                   ! =0 if at least 1 snow layer [flg]
     !integer :: trip                               ! flag: =1 to redo RT calculation if result is unrealistic
     !integer :: flg_dover                          ! defines conditions for RT redo (explained below)

     real(r8):: slope_rad, deg2rad
     real(r8):: albedo                             ! temporary snow albedo [frc]
     real(r8):: flx_sum                            ! temporary summation variable for NIR weighting
     real(r8):: albout_b                           ! snow albedo of the band [frc]
     real(r8):: flx_abs_b(-nlevsno_max+1:1)            ! absorbed flux per unit incident flux at top of snowpack (lyr) [frc]

     real(r8):: L_snw(-nlevsno_max+1:0)                ! h2o mass (liquid+solid) in snow layer (lyr) [kg/m2]
     real(r8):: tau_snw(-nlevsno_max+1:0)              ! snow optical depth (lyr) [unitless]
     real(r8):: L_aer(-nlevsno_max+1:0,sno_nbr_aer)    ! aerosol mass in snow layer (lyr,nbr_aer) [kg/m2]
     real(r8):: tau_aer(-nlevsno_max+1:0,sno_nbr_aer)  ! aerosol optical depth (lyr,nbr_aer) [unitless]
     real(r8):: tau_sum                            ! cumulative (snow+aerosol) optical depth [unitless]
     real(r8):: tau_elm(-nlevsno_max+1:0)              ! column optical depth from layer bottom to snowpack top (lyr) [unitless]
     real(r8):: omega_sum                          ! temporary summation of single-scatter albedo of all aerosols [frc]
     real(r8):: g_sum                              ! temporary summation of asymmetry parameter of all aerosols [frc]

     real(r8):: tau(-nlevsno_max+1:0)                  ! weighted optical depth of snow+aerosol layer (lyr) [unitless]
     real(r8):: omega(-nlevsno_max+1:0)                ! weighted single-scatter albedo of snow+aerosol layer (lyr) [frc]
     real(r8):: g(-nlevsno_max+1:0)                    ! weighted asymmetry parameter of snow+aerosol layer (lyr) [frc]
     real(r8):: tau_star(-nlevsno_max+1:0)             ! transformed (i.e. Delta-Eddington) optical depth of snow+aerosol layer
                                                   ! (lyr) [unitless]
     real(r8):: omega_star(-nlevsno_max+1:0)           ! transformed (i.e. Delta-Eddington) SSA of snow+aerosol layer (lyr) [frc]
     real(r8):: g_star(-nlevsno_max+1:0)               ! transformed (i.e. Delta-Eddington) asymmetry paramater of snow+aerosol layer
                                                   ! (lyr) [frc]

     integer :: nstep                              ! current timestep [nbr] (debugging only)
     integer :: g_idx, c_idx, l_idx                ! gridcell, column, and landunit indices [idx]
     integer :: bnd_idx                            ! spectral band index (1 <= bnd_idx <= numrad_snw) [idx]
     integer :: rds_idx                            ! snow effective radius index for retrieving
                                                   ! Mie parameters from lookup table [idx]
     integer :: snl_btm                            ! index of bottom snow layer (0) [idx]
     integer :: snl_top                            ! index of top snow layer (-4 to 0) [idx]
     integer :: fc                                 ! column filter index
     integer :: k                                  ! case index
     integer :: flg_slr_in                         ! =1 for direct-beam incident flux, =2 for diffuse
     integer :: rds_top                            ! snow grain radius of the top layer [microns]
     integer :: nbad                               ! elements with an error
     integer :: err_w(num_nourbanc,numrad_snw,ncall)  ! error code (fc,bnd,case): 0 none, 1 grain size out of the table,
                                                   ! 2 negative absorption, 3 energy conservation, 4 albedo > 1
     real(r8):: errv_w(num_nourbanc,numrad_snw,ncall) ! the offending value
     real(r8):: alb_w(num_nourbanc,numrad_snw,ncall)                  ! band albedo (fc,bnd,case)
     real(r8):: abs_w(num_nourbanc,-nlevsno+1:1,numrad_snw,nabs)          ! band absorbed flux (fc,lyr,bnd,case)
     integer :: i                                  ! layer index [idx]
     integer :: j                                  ! aerosol number index [idx]
     integer :: m                                  ! secondary layer index [idx]
     integer :: nint_snw_rds_min                   ! nearest integer value of snw_rds_min

     real(r8):: F_abs(-nlevsno_max+1:0)                ! net absorbed radiative energy (lyr) [W/m^2]
     real(r8):: F_abs_sum                          ! total absorbed energy in column [W/m^2]
     real(r8):: F_sfc_pls                          ! upward radiative flux at snowpack top [W/m^2]
     real(r8):: F_btm_net                          ! net flux at bottom of snowpack [W/m^2]
     real(r8):: energy_sum                         ! sum of all energy terms; should be 0.0 [W/m^2]
     real(r8):: mu_not                             ! cosine of solar zenith angle (used locally) [frc]

     integer :: err_idx                            ! counter for number of times through error loop [nbr]
     real(r8):: lat_coord                          ! gridcell latitude (debugging only)
     real(r8):: lon_coord                          ! gridcell longitude (debugging only)
     integer :: sfctype                            ! underlying surface type (debugging only)
     real(r8):: pi                                 ! 3.1415...

     integer :: snw_shp                            ! Snow grain shape option (every layer):
                                                   ! 1=sphere; 2=spheroid; 3=hexagonal plate; 4=koch snowflake
     real(r8), parameter :: snw_fs = 0._r8         ! Shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
                                                   ! 0=use recommended default value
                                                   ! others(0<fs<1)= use user-specified value
                                                   ! only activated when sno_shp > 1 (i.e. nonspherical)
     real(r8), parameter :: snw_ar = 0._r8         ! % Aspect ratio: ratio of grain width to length
                                                   ! 0=use recommended default value
                                                   ! others(0.1<fs<20)= use user-specified value
                                                   ! only activated when sno_shp > 1 (i.e. nonspherical)
     real(r8):: &
         diam_ice           , & ! effective snow grain diameter
         fs_sphd            , & ! shape factor for spheroid
         fs_hex0            , & ! shape factor for hexagonal plate
         fs_hex             , & ! shape factor for hexagonal plate (reference)
         fs_koch            , & ! shape factor for koch snowflake
         AR_tmp             , & ! aspect ratio for spheroid
         g_ice_Cg_tmp(7)    , & ! temporary for calculation of asymetry factor
         gg_ice_F07_tmp(7)  , & ! temporary for calculation of asymetry factor
         g_ice_F07          , & ! temporary for calculation of asymetry factor
         g_ice              , & ! asymmetry factor
         gg_F07_intp        , & ! temporary for calculation of asymetry factor (interpolated)
         g_Cg_intp          , & ! temporary for calculation of asymetry factor  (interpolated)
	 R_1_omega_tmp      , & ! temporary for dust-snow mixing calculation 
         C_dust_total           ! dust concentration
	
     integer :: atm_type_index  ! index for atmospheric type
     integer :: slr_zen         ! integer value of solar zenith angle

     ! SNICAR_AD new variables, follow sea-ice shortwave conventions
     real(r8):: &
        trndir(-nlevsno_max+1:1)  , & ! solar beam down transmission from top
        trntdr(-nlevsno_max+1:1)  , & ! total transmission to direct beam for layers above
        trndif(-nlevsno_max+1:1)  , & ! diffuse transmission to diffuse beam for layers above
        rupdir(-nlevsno_max+1:1)  , & ! reflectivity to direct radiation for layers below
        rupdif(-nlevsno_max+1:1)  , & ! reflectivity to diffuse radiation for layers below
        rdndif(-nlevsno_max+1:1)  , & ! reflectivity to diffuse radiation for layers above
        dfdir(-nlevsno_max+1:1)   , & ! down-up flux at interface due to direct beam at top surface
        dfdif(-nlevsno_max+1:1)   , & ! down-up flux at interface due to diffuse beam at top surface
        dftmp(-nlevsno_max+1:1)       ! temporary variable for down-up flux at interface

     real(r8):: &
        rdir(-nlevsno_max+1:0)       , & ! layer reflectivity to direct radiation
        rdif_a(-nlevsno_max+1:0)     , & ! layer reflectivity to diffuse radiation from above
        rdif_b(-nlevsno_max+1:0)     , & ! layer reflectivity to diffuse radiation from below
        tdir(-nlevsno_max+1:0)       , & ! layer transmission to direct radiation (solar beam + diffuse)
        tdif_a(-nlevsno_max+1:0)     , & ! layer transmission to diffuse radiation from above
        tdif_b(-nlevsno_max+1:0)     , & ! layer transmission to diffuse radiation from below
        trnlay(-nlevsno_max+1:0)         ! solar beam transm for layer (direct beam only)

     real(r8):: &
         ts       , & ! layer delta-scaled extinction optical depth
         ws       , & ! layer delta-scaled single scattering albedo
         gs       , & ! layer delta-scaled asymmetry parameter
         extins   , & ! extinction
         alp      , & ! temporary for alpha
         gam      , & ! temporary for agamm
         amg      , & ! alp - gam
         apg      , & ! alp + gam
         ue       , & ! temporary for u
         refk     , & ! interface multiple scattering
         refkp1   , & ! interface multiple scattering for k+1
         refkm1   , & ! interface multiple scattering for k-1
         tdrrdir  , & ! direct tran times layer direct ref
         tdndif       ! total down diffuse = tot tran - direct tran

     real(r8) :: &
         alpha    , & ! term in direct reflectivity and transmissivity
         agamm    , & ! term in direct reflectivity and transmissivity
         el       , & ! term in alpha,agamm,n,u
         taus     , & ! scaled extinction optical depth
         omgs     , & ! scaled single particle scattering albedo
         asys     , & ! scaled asymmetry parameter
         u        , & ! term in diffuse reflectivity and transmissivity
         n        , & ! term in diffuse reflectivity and transmissivity
         lm       , & ! temporary for el
         mu       , & ! cosine solar zenith for either snow or water
         ne           ! temporary for n

     ! perpendicular and parallel relative to plane of incidence and scattering
     real(r8) :: &
         R1       , & ! perpendicular polarization reflection amplitude
         R2       , & ! parallel polarization reflection amplitude
         T1       , & ! perpendicular polarization transmission amplitude
         T2       , & ! parallel polarization transmission amplitude
         Rf_dir_a , & ! fresnel reflection to direct radiation
         Tf_dir_a , & ! fresnel transmission to direct radiation
         Rf_dif_a , & ! fresnel reflection to diff radiation from above
         Rf_dif_b , & ! fresnel reflection to diff radiation from below
         Tf_dif_a , & ! fresnel transmission to diff radiation from above
         Tf_dif_b     ! fresnel transmission to diff radiation from below

     real(r8) :: &
         gwt      , & ! gaussian weight
         swt      , & ! sum of weights
         trn      , & ! layer transmission
         rdr      , & ! rdir for gaussian integration
         tdr      , & ! tdir for gaussian integration
         smr      , & ! accumulator for rdif gaussian integration
         smt      , & ! accumulator for tdif gaussian integration
         exp_min      ! minimum exponential value

     integer :: &
         ng             , & ! gaussian integration index
         snl_btm_itf        ! index of bottom snow layer interfaces (1) [idx]
     integer, parameter :: ngmax = 8 ! gaussian integration index

     ! Gaussian integration angle and coefficients
     real(r8), parameter :: difgauspt(1:8) = &   ! gaussian angles (radians)
         (/ 0.9894009_r8,  0.9445750_r8, 0.8656312_r8,  0.7554044_r8, &
            0.6178762_r8,  0.4580168_r8, 0.2816036_r8,  0.0950125_r8/)
     real(r8), parameter :: difgauswt(1:8) = &   ! gaussian weights
         (/ 0.0271525_r8,  0.0622535_r8, 0.0951585_r8,  0.1246290_r8, &
            0.1495960_r8,  0.1691565_r8, 0.1826034_r8,  0.1894506_r8/)
     ! real(r8),  dimension (1:8) :: &
     !     dif_gauspt     & ! gaussian angles (radians)
     !       = (/ 0.9894009_r8,  0.9445750_r8, &
     !            0.8656312_r8,  0.7554044_r8, &
     !            0.6178762_r8,  0.4580168_r8, &
     !            0.2816036_r8,  0.0950125_r8/) , &
     !     dif_gauswt     & ! gaussian weights
     !       = (/ 0.0271525_r8,  0.0622535_r8, &
     !            0.0951585_r8,  0.1246290_r8, &
     !            0.1495960_r8,  0.1691565_r8, &
     !            0.1826034_r8,  0.1894506_r8/)

     ! constants used in algorithm
     real(r8), parameter :: &
         c0      = 0.0_r8     , &
         c1      = 1.0_r8     , &
         c3      = 3.0_r8     , &
         c4      = 4.0_r8     , &
         c6      = 6.0_r8     , &
         cp01    = 0.01_r8    , &
         cp5     = 0.5_r8     , &
         cp75    = 0.75_r8    , &
         c1p5    = 1.5_r8     , &
         trmin   = 0.001_r8   , &
         argmax  = 10.0_r8       ! maximum argument of exponential

     ! cconstant coefficients used for SZA parameterization
     real(r8), parameter :: &
         sza_a0 =  0.085730_r8 , &
         sza_a1 = -0.630883_r8 , &
         sza_a2 =  1.303723_r8 , &
         sza_b0 =  1.467291_r8 , &
         sza_b1 = -3.338043_r8 , &
         sza_b2 =  6.807489_r8 , &
         puny   =  1.0e-11_r8  , &
         mu_75  =  0.2588_r8       ! cosine of 75 degree

     ! coefficients used for SZA parameterization
     real(r8) :: &
         sza_c1          , & ! coefficient, SZA parameteirzation
         sza_c0          , & ! coefficient, SZA parameterization
         sza_factor      , & ! factor used to adjust NIR direct albedo
         flx_sza_adjust  , & ! direct NIR flux adjustment from sza_factor
         mu0                 ! incident solar zenith angle

     ! Delta-Eddington solution expressions
     ! alpha(w,uu,gg,e) = p75*w*uu*((c1 + gg*(c1-w))/(c1 - e*e*uu*uu))
     ! agamm(w,uu,gg,e) = p5*w*((c1 + c3*gg*(c1-w)*uu*uu)/(c1-e*e*uu*uu))
     ! n(uu,et)         = ((uu+c1)*(uu+c1)/et ) - ((uu-c1)*(uu-c1)*et)
     ! u(w,gg,e)        = c1p5*(c1 - w*gg)/e
     ! el(w,gg)         = sqrt(c3*(c1-w)*(c1 - w*gg))

     !-----------------------------------------------------------------------
         !mgf++
         !mgf--

      ! Constants for non-spherical ice particles and dust-snow internal mixing
      real(r8), parameter :: g_b2(7) = (/-2.66792E-04_r8,1.14088E-03_r8,2.37011E-04_r8,-2.35905E-04_r8,8.40449E-04_r8,-4.71484E-04_r8,-4.71484E-04_r8/)
      real(r8), parameter :: g_b1(7) = (/5.21042E-01_r8,4.96181E-01_r8,1.83711E-01_r8,1.37082E-01_r8,5.50598E-02_r8,8.48743E-02_r8,8.48743E-02_r8/)
      real(r8), parameter :: g_b0(7) = (/9.76029E-01_r8,9.67798E-01_r8,1.00111E+00_r8,1.00224E+00_r8,9.64295E-01_r8,9.97475E-01_r8,9.97475E-01_r8/)
      real(r8), parameter :: g_F07_c2(7) = (/1.349959E-1_r8,1.115697E-1_r8,9.853958E-2_r8,5.557793E-2_r8,-1.233493E-1_r8,0.0_r8,0.0_r8/)
      real(r8), parameter :: g_F07_c1(7) = (/-3.987320E-1_r8,-3.723287E-1_r8,-3.924784E-1_r8,-3.259404E-1_r8,4.429054E-2_r8,-1.726586E-1_r8,-1.726586E-1_r8/)
      real(r8), parameter :: g_F07_c0(7) = (/7.938904E-1_r8,8.030084E-1_r8,8.513932E-1_r8,8.692241E-1_r8,7.085850E-1_r8,6.412701E-1_r8,6.412701E-1_r8/)
      real(r8), parameter :: g_F07_p2(7) = (/3.165543E-3_r8,2.014810E-3_r8,1.780838E-3_r8,6.987734E-4_r8,-1.882932E-2_r8,-2.277872E-2_r8,-2.277872E-2_r8/)
      real(r8), parameter :: g_F07_p1(7) = (/1.140557E-1_r8,1.143152E-1_r8,1.143814E-1_r8,1.071238E-1_r8,1.353873E-1_r8,1.914431E-1_r8,1.914431E-1_r8/)
      real(r8), parameter :: g_F07_p0(7) = (/5.292852E-1_r8,5.425909E-1_r8,5.601598E-1_r8,6.023407E-1_r8,6.473899E-1_r8,4.634944E-1_r8,4.634944E-1_r8/)
      real(r8), parameter :: dust_clear_d0(3) = (/1.0413E+00_r8,1.0168E+00_r8,1.0189E+00_r8/)
      real(r8), parameter :: dust_clear_d1(3) = (/1.0016E+00_r8,1.0070E+00_r8,1.0840E+00_r8/)
      real(r8), parameter :: dust_clear_d2(3) = (/2.4208E-01_r8,1.5300E-03_r8,1.1230E-04_r8/)
      real(r8), parameter :: dust_cloudy_d0(3) = (/1.0388E+00_r8,1.0167E+00_r8,1.0189E+00_r8/)
      real(r8), parameter :: dust_cloudy_d1(3) = (/1.0015E+00_r8,1.0061E+00_r8,1.0823E+00_r8/)
      real(r8), parameter :: dust_cloudy_d2(3) = (/2.5973E-01_r8,1.6200E-03_r8,1.1721E-04_r8/)
      !!! factors for considering snow grain shape
      !!! factors for considring dust-snow internal mixing
     ! Enforce expected array sizes

     associate(&
          snl         =>   col_pp%snl           , & ! Input:  [integer (:)]  negative number of snow layers (col) [nbr]
          h2osno      =>   col_ws%h2osno        , & ! Input:  [real(r8) (:)]  snow liquid water equivalent (col) [kg/m2]
          h2osoi_liq  =>   col_ws%h2osoi_liq    , & ! Input:  [real(r8) (:,:)]  liquid water content (col,lyr) [kg/m2]
          h2osoi_ice  =>   col_ws%h2osoi_ice    , & ! Input:  [real(r8) (:,:)]  ice content (col,lyr) [kg/m2]
          snw_rds     =>   col_ws%snw_rds       , & ! Input:  [real(r8) (:,:)]  snow grain radius (col,lyr) [microns]
          bcphi       =>   aerosol_vars%mss_cnc_bcphi_col , & ! Input: [real(r8) (:,:)] aerosol mass concentrations (col,lyr) [kg/kg]
          bcpho       =>   aerosol_vars%mss_cnc_bcpho_col , &
          ocphi       =>   aerosol_vars%mss_cnc_ocphi_col , &
          ocpho       =>   aerosol_vars%mss_cnc_ocpho_col , &
          dst1        =>   aerosol_vars%mss_cnc_dst1_col  , &
          dst2        =>   aerosol_vars%mss_cnc_dst2_col  , &
          dst3        =>   aerosol_vars%mss_cnc_dst3_col  , &
          dst4        =>   aerosol_vars%mss_cnc_dst4_col    &
          )

       ! Define constants
       pi = SHR_CONST_PI
       nint_snw_rds_min = nint(snw_rds_min)
       deg2rad = SHR_CONST_PI/180._r8

       ! always use Delta approximation for snow
       DELTA = 1

       snw_shp = snow_shape_sphere
       atm_type_index = atm_type_default
       
       ! Define snow grain shape
       if (trim(snow_shape) == 'sphere') then
         snw_shp = snow_shape_sphere
       elseif (trim(snow_shape) == 'spheroid') then
         snw_shp = snow_shape_spheroid
       elseif (trim(snow_shape) == 'hexagonal_plate') then
	 snw_shp = snow_shape_hexagonal_plate
       elseif (trim(snow_shape) == 'koch_snowflake') then
         snw_shp = snow_shape_koch_snowflake
       else
	 write(iulog,*) "snow_shape = ", snow_shape
         call endrun( "snow_shape is unknown" )
       endif
	  	    
       ! Define atmospheric type
       if (trim(snicar_atm_type) == 'default') then
         atm_type_index = atm_type_default
       elseif (trim(snicar_atm_type) == 'mid-latitude_winter') then
         atm_type_index = atm_type_mid_latitude_winter
       elseif (trim(snicar_atm_type) == 'mid-latitude_summer') then
         atm_type_index = atm_type_mid_latitude_summer
       elseif (trim(snicar_atm_type) == 'sub-Arctic_winter') then
         atm_type_index = atm_type_sub_Arctic_winter
       elseif (trim(snicar_atm_type) == 'sub-Arctic_summer') then
         atm_type_index = atm_type_sub_Arctic_summer
       elseif (trim(snicar_atm_type) == 'summit_Greenland') then
         atm_type_index = atm_type_summit_Greenland
       elseif (trim(snicar_atm_type) == 'high_mountain') then
         atm_type_index = atm_type_high_mountain
       else
	 write(iulog,*) "snicar_atm_type = ", snicar_atm_type
         call endrun( "snicar_atm_type is unknown" )
       endif
	  
       exp_min = exp(-argmax)
      if (elm_device_ready) call endrun(msg="the non-modal aerosol path (compsets without BC) is not supported with GPU enabled")

       ! Each (case, column, band) is an independent radiative transfer problem: the column
       ! set-up below is repeated per band, and the bands are combined in a second loop.
       ! one thread per (case, column, band); the private arrays are the thread's layer profiles
       do k = 1, ncall
       do fc = 1,num_nourbanc
       do bnd_idx = 1,numrad_snw
          flg_slr_in = call_slr(k)
          c_idx = filter_nourbanc(fc)
          g_idx = col_pp%gridcell(c_idx)
          flx_abs_b(:) = 0._r8
          err_w(fc,bnd_idx,k) = 0
          h2osno_lcl = h2osno(c_idx)

          ! Qualifier for computing snow RT: sunlight from the atmosphere model and a minimum
          ! amount of snow on the ground (otherwise the second loop sets the snow albedo)
          if ((coszen(c_idx) > 0._r8) .and. (h2osno_lcl > min_snw) ) then

             ! If there is snow, but zero snow layers, we must create a layer locally.
             ! This layer is presumed to have the fresh snow effective radius.
             if (snl(c_idx) > -1) then
                flg_nosnl         =  1
                snl_lcl           =  -1
                h2osno_ice_lcl(0) =  h2osno_lcl
                h2osno_liq_lcl(0) =  0._r8
                snw_rds_lcl(0)    =  nint_snw_rds_min
             else
                flg_nosnl         =  0
                snl_lcl           =  snl(c_idx)
                h2osno_liq_lcl(-nlevsno+1:0) =  h2osoi_liq(c_idx,-nlevsno+1:0)
                h2osno_ice_lcl(-nlevsno+1:0) =  h2osoi_ice(c_idx,-nlevsno+1:0)
                snw_rds_lcl(-nlevsno+1:0)    =  nint(snw_rds(c_idx,-nlevsno+1:0))
             endif

             snl_btm   = 0
             snl_top   = snl_lcl+1

             ! for debugging only
             l_idx     = col_pp%landunit(c_idx)
             sfctype   = lun_pp%itype(l_idx)
             lat_coord = grc_pp%latdeg(g_idx)
             lon_coord = grc_pp%londeg(g_idx)

             if (use_finetop_rad) then
             slope_rad = grc_pp%slope_deg(g_idx) * deg2rad
                
             if (snl(c_idx) > -1) then
             h2osno_liq_lcl(0) = h2osno_liq_lcl(0) * cos(slope_rad)
             h2osno_ice_lcl(0) = h2osno_ice_lcl(0) * cos(slope_rad)
             else
             h2osno_liq_lcl(-nlevsno+1:0) = h2osno_liq_lcl(-nlevsno+1:0) * cos(slope_rad)
             h2osno_ice_lcl(-nlevsno+1:0) = h2osno_ice_lcl(-nlevsno+1:0) * cos(slope_rad)
             endif

             h2osno_lcl = h2osno_lcl * cos(slope_rad)
             endif


             ! the case's aerosol species
             mss_cnc_aer_lcl(:,:) = 0._r8
             if (call_aer(1,k)) mss_cnc_aer_lcl(-nlevsno+1:0,1) = bcphi(c_idx,:)
             if (call_aer(2,k)) mss_cnc_aer_lcl(-nlevsno+1:0,2) = bcpho(c_idx,:)
             if (call_aer(3,k)) mss_cnc_aer_lcl(-nlevsno+1:0,3) = ocphi(c_idx,:)
             if (call_aer(4,k)) mss_cnc_aer_lcl(-nlevsno+1:0,4) = ocpho(c_idx,:)
             if (call_aer(5,k)) mss_cnc_aer_lcl(-nlevsno+1:0,5) = dst1(c_idx,:)
             if (call_aer(6,k)) mss_cnc_aer_lcl(-nlevsno+1:0,6) = dst2(c_idx,:)
             if (call_aer(7,k)) mss_cnc_aer_lcl(-nlevsno+1:0,7) = dst3(c_idx,:)
             if (call_aer(8,k)) mss_cnc_aer_lcl(-nlevsno+1:0,8) = dst4(c_idx,:)

             ! Error check for snow grain size (reported after the loops)
             do i=snl_top,snl_btm,1
                if ((snw_rds_lcl(i) < snw_rds_min_tbl) .or. (snw_rds_lcl(i) > snw_rds_max_tbl)) then
                   if (err_w(fc,bnd_idx,k) == 0) then
                      err_w(fc,bnd_idx,k) = 1
                      errv_w(fc,bnd_idx,k) = snw_rds_lcl(i)
                   endif
                endif
             enddo
             if (err_w(fc,bnd_idx,k) == 0) then


               ! note that we can remove flg_dover since this algorithm is
               ! stable for mu_not > 0.01

               ! mu_not is cosine solar zenith angle above the fresnel level; make
               ! sure mu_not is large enough for stable and meaningful radiation
               ! solution: .01 is like sun just touching horizon with its lower edge
               ! equivalent to mu0 in sea-ice shortwave model ice_shortwave.F90
                mu_not = max(coszen(c_idx), cp01)


                   ! Set direct or diffuse incident irradiance to 1
                   ! (This has to be within the bnd loop because mu_not is adjusted in rare cases)
                   if (flg_slr_in == 1) then
                      flx_slrd_b = 1._r8/(mu_not*pi) ! this corresponds to incident irradiance of 1.0
                      flx_slri_b = 0._r8
                   else
                      flx_slrd_b = 0._r8
                      flx_slri_b = 1._r8
                   endif

                   ! Pre-emptive error handling: aerosols can reap havoc on these absorptive bands.
                   ! Since extremely high soot concentrations have a negligible effect on these bands, zero them.
                   if ( (numrad_snw == 5).and.((bnd_idx == 5).or.(bnd_idx == 4)) ) then
                      mss_cnc_aer_lcl(:,:) = 0._r8
                   endif

                   if ( (numrad_snw == 3).and.(bnd_idx == 3) ) then
                      mss_cnc_aer_lcl(:,:) = 0._r8
                   endif

                   ! Define local Mie parameters based on snow grain size and aerosol species,
                   !  retrieved from a lookup table.
                   if (flg_slr_in == 1) then
                      do i=snl_top,snl_btm,1
                         rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
                         ! snow optical properties (direct radiation)
                         ss_alb_snw_lcl(i)      = ss_alb_snw_drc(rds_idx,bnd_idx)
                         asm_prm_snw_lcl(i)     = asm_prm_snw_drc(rds_idx,bnd_idx)
                         ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_drc(rds_idx,bnd_idx)
                      enddo
                   elseif (flg_slr_in == 2) then
                      do i=snl_top,snl_btm,1
                         rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
                         ! snow optical properties (diffuse radiation)
                         ss_alb_snw_lcl(i)      = ss_alb_snw_dfs(rds_idx,bnd_idx)
                         asm_prm_snw_lcl(i)     = asm_prm_snw_dfs(rds_idx,bnd_idx)
                         ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_dfs(rds_idx,bnd_idx)
                      enddo
                   endif
				   
                  ! Calculate the asymetry factors under different snow grain shapes
                   do i=snl_top,snl_btm,1
                      if(snw_shp == snow_shape_spheroid) then ! spheroid
                         diam_ice = 2._r8*snw_rds_lcl(i)
                         if(snw_fs == 0._r8) then
                            fs_sphd = 0.929_r8
                         else
                            fs_sphd = snw_fs               
                         endif
                         fs_hex = 0.788_r8 
                         if(snw_ar == 0._r8) then
                            AR_tmp = 0.5_r8
                         else
                            AR_tmp = snw_ar              
                         endif
                         g_ice_Cg_tmp = g_b0 * ((fs_sphd/fs_hex)**g_b1) * (diam_ice**g_b2)
                         gg_ice_F07_tmp = g_F07_c0 + g_F07_c1 * AR_tmp + g_F07_c2 * (AR_tmp**2)			 
                      elseif(snw_shp == snow_shape_hexagonal_plate) then ! hexagonal plate
                         diam_ice = 2._r8*snw_rds_lcl(i)
                         if(snw_fs == 0._r8) then
                            fs_hex0 = 0.788_r8
                         else
                            fs_hex0 = snw_fs               
                         endif
                         fs_hex = 0.788_r8 
                         if(snw_ar == 0._r8) then
                           AR_tmp = 2.5_r8
                         else
                           AR_tmp = snw_ar              
                         endif
                         g_ice_Cg_tmp = g_b0 * ((fs_hex0/fs_hex)**g_b1) * (diam_ice**g_b2)
                         gg_ice_F07_tmp = g_F07_p0 + g_F07_p1 * log(AR_tmp) + g_F07_p2 * ((log(AR_tmp))**2)
                      elseif(snw_shp == snow_shape_koch_snowflake) then ! Koch snowflake
                         diam_ice = 2._r8 * snw_rds_lcl(i) /0.544_r8
                         if(snw_fs == 0._r8) then
                            fs_koch = 0.712_r8
                         else
                            fs_koch = snw_fs               
                         endif
                         fs_hex = 0.788_r8 
                         if(snw_ar == 0._r8) then
                            AR_tmp = 2.5_r8
                         else
                            AR_tmp = snw_ar              
                         endif
                         g_ice_Cg_tmp = g_b0 * ((fs_koch/fs_hex)**g_b1) * (diam_ice**g_b2)
                         gg_ice_F07_tmp = g_F07_p0 + g_F07_p1 * log(AR_tmp) + g_F07_p2 * ((log(AR_tmp))**2)	 
                     endif

                     ! Linear interpolation for calculating the asymetry factor at band_idx.
                     if(snw_shp > 1) then
                       if(bnd_idx == 1) then
                         g_Cg_intp = (g_ice_Cg_tmp(2)-g_ice_Cg_tmp(1))/(1.055_r8-0.475_r8)*(0.5_r8-0.475_r8)+g_ice_Cg_tmp(1)
                         gg_F07_intp = (gg_ice_F07_tmp(2)-gg_ice_F07_tmp(1))/(1.055_r8-0.475_r8)*(0.5_r8-0.475_r8)+gg_ice_F07_tmp(1)
                       elseif(bnd_idx == 2) then 
                         g_Cg_intp = (g_ice_Cg_tmp(2)-g_ice_Cg_tmp(1))/(1.055_r8-0.475_r8)*(0.85_r8-0.475_r8)+g_ice_Cg_tmp(1)
                         gg_F07_intp = (gg_ice_F07_tmp(2)-gg_ice_F07_tmp(1))/(1.055_r8-0.475_r8)*(0.85_r8-0.475_r8)+gg_ice_F07_tmp(1)
                       elseif(bnd_idx == 3) then 
                         g_Cg_intp = (g_ice_Cg_tmp(3)-g_ice_Cg_tmp(2))/(1.655_r8-1.055_r8)*(1.1_r8-1.055_r8)+g_ice_Cg_tmp(2)
                         gg_F07_intp = (gg_ice_F07_tmp(3)-gg_ice_F07_tmp(2))/(1.655_r8-1.055_r8)*(1.1_r8-1.055_r8)+gg_ice_F07_tmp(2)
                       elseif(bnd_idx == 4) then 
                         g_Cg_intp = (g_ice_Cg_tmp(3)-g_ice_Cg_tmp(2))/(1.655_r8-1.055_r8)*(1.35_r8-1.055_r8)+g_ice_Cg_tmp(2)
                         gg_F07_intp = (gg_ice_F07_tmp(3)-gg_ice_F07_tmp(2))/(1.655_r8-1.055_r8)*(1.35_r8-1.055_r8)+gg_ice_F07_tmp(2)
                       elseif(bnd_idx == 5) then
                         g_Cg_intp = (g_ice_Cg_tmp(6)-g_ice_Cg_tmp(5))/(3.75_r8-3.0_r8)*(3.25_r8-3.0_r8)+g_ice_Cg_tmp(5)
                         gg_F07_intp = (gg_ice_F07_tmp(6)-gg_ice_F07_tmp(5))/(3.75_r8-3.0_r8)*(3.25_r8-3.0_r8)+gg_ice_F07_tmp(5)
                       endif
                       g_ice_F07 = gg_F07_intp + (1._r8 - gg_F07_intp) / ss_alb_snw_lcl(i) / 2._r8
                       g_ice = g_ice_F07 * g_Cg_intp
                       asm_prm_snw_lcl(i) = g_ice
                     endif

                     if(asm_prm_snw_lcl(i) > 0.99_r8) then 
                       asm_prm_snw_lcl(i) = 0.99_r8
                     endif                        

                  enddo
                  !!!-end


   !H. Wang
                   ! aerosol species 1 optical properties
                  ! ss_alb_aer_lcl(1)        = ss_alb_bc1(bnd_idx)
                  ! asm_prm_aer_lcl(1)       = asm_prm_bc1(bnd_idx)
                  ! ext_cff_mss_aer_lcl(1)   = ext_cff_mss_bc1(bnd_idx)

                   ! aerosol species 2 optical properties
                  ! ss_alb_aer_lcl(2)        = ss_alb_bc2(bnd_idx)
                  ! asm_prm_aer_lcl(2)       = asm_prm_bc2(bnd_idx)
                  ! ext_cff_mss_aer_lcl(2)   = ext_cff_mss_bc2(bnd_idx)
   !H. Wang
                   ! aerosol species 3 optical properties
                   ss_alb_aer_lcl(3)        = ss_alb_oc1(bnd_idx)
                   asm_prm_aer_lcl(3)       = asm_prm_oc1(bnd_idx)
                   ext_cff_mss_aer_lcl(3)   = ext_cff_mss_oc1(bnd_idx)

                   ! aerosol species 4 optical properties
                   ss_alb_aer_lcl(4)        = ss_alb_oc2(bnd_idx)
                   asm_prm_aer_lcl(4)       = asm_prm_oc2(bnd_idx)
                   ext_cff_mss_aer_lcl(4)   = ext_cff_mss_oc2(bnd_idx)

                   ! aerosol species 5 optical properties
                   ss_alb_aer_lcl(5)        = ss_alb_dst1(bnd_idx)
                   asm_prm_aer_lcl(5)       = asm_prm_dst1(bnd_idx)
                   ext_cff_mss_aer_lcl(5)   = ext_cff_mss_dst1(bnd_idx)

                   ! aerosol species 6 optical properties
                   ss_alb_aer_lcl(6)        = ss_alb_dst2(bnd_idx)
                   asm_prm_aer_lcl(6)       = asm_prm_dst2(bnd_idx)
                   ext_cff_mss_aer_lcl(6)   = ext_cff_mss_dst2(bnd_idx)

                   ! aerosol species 7 optical properties
                   ss_alb_aer_lcl(7)        = ss_alb_dst3(bnd_idx)
                   asm_prm_aer_lcl(7)       = asm_prm_dst3(bnd_idx)
                   ext_cff_mss_aer_lcl(7)   = ext_cff_mss_dst3(bnd_idx)

                   ! aerosol species 8 optical properties
                   ss_alb_aer_lcl(8)        = ss_alb_dst4(bnd_idx)
                   asm_prm_aer_lcl(8)       = asm_prm_dst4(bnd_idx)
                   ext_cff_mss_aer_lcl(8)   = ext_cff_mss_dst4(bnd_idx)


                   ! 1. snow and aerosol layer column mass (L_snw, L_aer [kg/m^2])
                   ! 2. optical Depths (tau_snw, tau_aer)
                   ! 3. weighted Mie properties (tau, omega, g)

                   ! Weighted Mie parameters of each layer
                   do i=snl_top,snl_btm,1
                    ! bulk aerosol treatment (BC optical properties independent
                    ! of BC and ice grain size)
                    ! aerosol species 1 optical properties (within-ice BC)
                    ss_alb_aer_lcl(1)        = ss_alb_bc1(bnd_idx)
                    asm_prm_aer_lcl(1)       = asm_prm_bc1(bnd_idx)
                    ext_cff_mss_aer_lcl(1)   = ext_cff_mss_bc1(bnd_idx)

                    ! aerosol species 2 optical properties
                    ss_alb_aer_lcl(2)        = ss_alb_bc2(bnd_idx)
                    asm_prm_aer_lcl(2)       = asm_prm_bc2(bnd_idx)
                    ext_cff_mss_aer_lcl(2)   = ext_cff_mss_bc2(bnd_idx)

                    ! Calculate single-scattering albedo for internal mixing of dust-snow
                    if (use_dust_snow_internal_mixing) then
                        if (bnd_idx < 4) then
                           C_dust_total = mss_cnc_aer_lcl(i,5) + mss_cnc_aer_lcl(i,6) + mss_cnc_aer_lcl(i,7) + mss_cnc_aer_lcl(i,8)
                           C_dust_total = C_dust_total * 1.0E+06_r8 
                           if(C_dust_total > 0._r8) then
                              if (flg_slr_in == 1) then
                                 R_1_omega_tmp = dust_clear_d0(bnd_idx) + dust_clear_d2(bnd_idx)*(C_dust_total**dust_clear_d1(bnd_idx))                   
                              else
                                 R_1_omega_tmp = dust_cloudy_d0(bnd_idx) + dust_cloudy_d2(bnd_idx)*(C_dust_total**dust_cloudy_d1(bnd_idx))   
                              endif					  
                              ss_alb_snw_lcl(i) = 1.0_r8 - (1.0_r8 - ss_alb_snw_lcl(i)) *R_1_omega_tmp
                           endif
                        endif
                        do j = 5,8,1
                           ss_alb_aer_lcl(j)        = 0._r8
                           asm_prm_aer_lcl(j)       = 0._r8
                           ext_cff_mss_aer_lcl(j)   = 0._r8
                        enddo
                    endif
					
                    !mgf--

                      L_snw(i)   = h2osno_ice_lcl(i)+h2osno_liq_lcl(i)
                      tau_snw(i) = L_snw(i)*ext_cff_mss_snw_lcl(i)

                      do j=1,sno_nbr_aer
                         if (use_dust_snow_internal_mixing .and. (j >= 5)) then
                           L_aer(i,j)  = 0._r8
                         else
                           L_aer(i,j)   = L_snw(i)*mss_cnc_aer_lcl(i,j)
                         endif
                         tau_aer(i,j) = L_aer(i,j)*ext_cff_mss_aer_lcl(j)
                      enddo

                      tau_sum   = 0._r8
                      omega_sum = 0._r8
                      g_sum     = 0._r8

                      do j=1,sno_nbr_aer
                         tau_sum    = tau_sum + tau_aer(i,j)
                         omega_sum  = omega_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j))
                         g_sum      = g_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j)*asm_prm_aer_lcl(j))
                      enddo

                      tau(i)    = tau_sum + tau_snw(i)
                      omega(i)  = (1/tau(i))*(omega_sum+(ss_alb_snw_lcl(i)*tau_snw(i)))
                      g(i)      = (1/(tau(i)*omega(i)))*(g_sum+ (asm_prm_snw_lcl(i)*ss_alb_snw_lcl(i)*tau_snw(i)))
                   enddo ! endWeighted Mie parameters of each layer

                   ! DELTA transformations, if requested
                   if (DELTA == 1) then
                      do i=snl_top,snl_btm,1
                         g_star(i)     = g(i)/(1+g(i))
                         omega_star(i) = ((1-(g(i)**2))*omega(i)) / (1-(omega(i)*(g(i)**2)))
                         tau_star(i)   = (1-(omega(i)*(g(i)**2)))*tau(i)
                      enddo
                   else
                      do i=snl_top,snl_btm,1
                         g_star(i)     = g(i)
                         omega_star(i) = omega(i)
                         tau_star(i)   = tau(i)
                      enddo
                   endif

                   ! Begin radiative transfer solver
                   ! Given input vertical profiles of optical properties, evaluate the
                   ! monochromatic Delta-Eddington adding-doubling solution

                   ! note that trndir, trntdr, trndif, rupdir, rupdif, rdndif
                   ! are variables at the layer interface,
                   ! for snow with layers rangeing from snl_top to snl_btm
                   ! there are snl_top to snl_btm+1 layer interface
                   snl_btm_itf = snl_btm + 1

                   do i = snl_top,snl_btm_itf,1
                      trndir(i) = c0
                      trntdr(i) = c0
                      trndif(i) = c0
                      rupdir(i) = c0
                      rupdif(i) = c0
                      rdndif(i) = c0
                   enddo

                   ! initialize top interface of top layer
                   trndir(snl_top) = c1
                   trntdr(snl_top) = c1
                   trndif(snl_top) = c1
                   rdndif(snl_top) = c0

                  ! begin main level loop
                  ! for layer interfaces except for the very bottom
                  do i = snl_top,snl_btm,1

                     ! initialize all layer apparent optical properties to 0
                     rdir  (i) = c0
                     rdif_a(i) = c0
                     rdif_b(i) = c0
                     tdir  (i) = c0
                     tdif_a(i) = c0
                     tdif_b(i) = c0
                     trnlay(i) = c0

                     ! compute next layer Delta-eddington solution only if total transmission
                     ! of radiation to the interface just above the layer exceeds trmin.

                     if (trntdr(i) > trmin ) then

                        ! calculation over layers with penetrating radiation

                        ! delta-transformed single-scattering properties
                        ! of this layer
                        ts = tau_star(i)
                        ws = omega_star(i)
                        gs = g_star(i)

                       ! Delta-Eddington solution expressions
                        ! n(uu,et)         = ((uu+c1)*(uu+c1)/et ) - ((uu-c1)*(uu-c1)*et)
                        ! u(w,gg,e)        = c1p5*(c1 - w*gg)/e
                        ! el(w,gg)         = sqrt(c3*(c1-w)*(c1 - w*gg))
                        lm = sqrt(c3*(c1-ws)*(c1 - ws*gs))  !lm = el(ws,gs)
                        ue = c1p5*(c1 - ws*gs)/lm           !ue = u(ws,gs,lm)
                        extins = max(exp_min, exp(-lm*ts))
                        ne = ((ue+c1)*(ue+c1)/extins) - ((ue-c1)*(ue-c1)*extins) !ne = n(ue,extins)

                        ! first calculation of rdif, tdif using Delta-Eddington formulas
                        ! rdif_a(k) = (ue+c1)*(ue-c1)*(c1/extins - extins)/ne
                        rdif_a(i) = (ue**2-c1)*(c1/extins - extins)/ne
                        tdif_a(i) = c4*ue/ne

                        ! evaluate rdir,tdir for direct beam
                        trnlay(i) = max(exp_min, exp(-ts/mu_not))

                        ! Delta-Eddington solution expressions
                        ! alpha(w,uu,gg,e) = p75*w*uu*((c1 + gg*(c1-w))/(c1 - e*e*uu*uu))
                        ! agamm(w,uu,gg,e) = p5*w*((c1 + c3*gg*(c1-w)*uu*uu)/(c1-e*e*uu*uu))
                        ! alp = alpha(ws,mu_not,gs,lm)
                        ! gam = agamm(ws,mu_not,gs,lm)
                        alp = cp75*ws*mu_not*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu_not*mu_not))
                        gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu_not*mu_not)/(c1-lm*lm*mu_not*mu_not))
                        apg = alp + gam
                        amg = alp - gam

                        rdir(i) = apg*rdif_a(i) +  amg*(tdif_a(i)*trnlay(i) - c1)
                        tdir(i) = apg*tdif_a(i) + (amg* rdif_a(i)-apg+c1)*trnlay(i)

                        ! recalculate rdif,tdif using direct angular integration over rdir,tdir,
                        ! since Delta-Eddington rdif formula is not well-behaved (it is usually
                        ! biased low and can even be negative); use ngmax angles and gaussian
                        ! integration for most accuracy:
                        R1 = rdif_a(i) ! use R1 as temporary
                        T1 = tdif_a(i) ! use T1 as temporary
                        swt = c0
                        smr = c0
                        smt = c0
                        do ng=1,ngmax
                           mu  = difgauspt(ng)
                           gwt = difgauswt(ng)
                           swt = swt + mu*gwt
                           trn = max(exp_min, exp(-ts/mu))
                           ! alp = alpha(ws,mu,gs,lm)
                           ! gam = agamm(ws,mu,gs,lm)
                           alp = cp75*ws*mu*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu*mu))
                           gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu*mu)/(c1-lm*lm*mu*mu))
                           apg = alp + gam
                           amg = alp - gam
                           rdr = apg*R1 + amg*T1*trn - amg
                           tdr = apg*T1 + amg*R1*trn - apg*trn + trn
                           smr = smr + mu*rdr*gwt
                           smt = smt + mu*tdr*gwt
                        enddo      ! ng
                        rdif_a(i) = smr/swt
                        tdif_a(i) = smt/swt

                        ! homogeneous layer
                        rdif_b(i) = rdif_a(i)
                        tdif_b(i) = tdif_a(i)

                      endif ! trntdr(k) > trmin

                      ! Calculate the solar beam transmission, total transmission, and
                      ! reflectivity for diffuse radiation from below at interface i,
                      ! the top of the current layer k:
                      !
                      !              layers       interface
                      !
                      !       ---------------------  i-1
                      !                i-1
                      !       ---------------------  i
                      !                 i
                      !       ---------------------

                      trndir(i+1) = trndir(i)*trnlay(i)
                      refkm1      = c1/(c1 - rdndif(i)*rdif_a(i))
                      tdrrdir     = trndir(i)*rdir(i)
                      tdndif      = trntdr(i) - trndir(i)
                      trntdr(i+1) = trndir(i)*tdir(i) + &
                           (tdndif + tdrrdir*rdndif(i))*refkm1*tdif_a(i)
                      rdndif(i+1) = rdif_b(i) + &
                           (tdif_b(i)*rdndif(i)*refkm1*tdif_a(i))
                      trndif(i+1) = trndif(i)*refkm1*tdif_a(i)

                  enddo       ! i    end main level loop


                  ! compute reflectivity to direct and diffuse radiation for layers
                  ! below by adding succesive layers starting from the underlying
                  ! ground and working upwards:
                  !
                  !              layers       interface
                  !
                  !       ---------------------  i
                  !                 i
                  !       ---------------------  i+1
                  !                i+1
                  !       ---------------------

                  ! set the underlying ground albedo == albedo of near-IR
                  ! unless bnd_idx == 1, for visible
                  rupdir(snl_btm_itf) = albsfc(c_idx,2)
                  rupdif(snl_btm_itf) = albsfc(c_idx,2)
                  if (bnd_idx == 1) then
                      rupdir(snl_btm_itf) = albsfc(c_idx,1)
                      rupdif(snl_btm_itf) = albsfc(c_idx,1)
                  endif

                  do i=snl_btm,snl_top,-1
                     ! interface scattering
                     refkp1        = c1/( c1 - rdif_b(i)*rupdif(i+1))
                     ! dir from top layer plus exp tran ref from lower layer, interface
                     ! scattered and tran thru top layer from below, plus diff tran ref
                     ! from lower layer with interface scattering tran thru top from below
                     rupdir(i) = rdir(i) &
                          + (        trnlay(i)  *rupdir(i+1) &
                          +  (tdir(i)-trnlay(i))*rupdif(i+1))*refkp1*tdif_b(i)
                     ! dif from top layer from above, plus dif tran upwards reflected and
                     ! interface scattered which tran top from below
                     rupdif(i) = rdif_a(i) + tdif_a(i)*rupdif(i+1)*refkp1*tdif_b(i)
                  enddo       ! i

                  ! net flux (down-up) at each layer interface from the
                  ! snow top (i = snl_top) to bottom interface above land (i = snl_btm_itf)
                  ! the interface reflectivities and transmissivities required
                  ! to evaluate interface fluxes are returned from solution_dEdd;
                  ! now compute up and down fluxes for each interface, using the
                  ! combined layer properties at each interface:
                  !
                  !              layers       interface
                  !
                  !       ---------------------  i
                  !                 i
                  !       ---------------------

                  do i = snl_top, snl_btm_itf
                     ! interface scattering
                     refk          = c1/(c1 - rdndif(i)*rupdif(i))
                     ! dir tran ref from below times interface scattering, plus diff
                     ! tran and ref from below times interface scattering
                     ! fdirup(i) = (trndir(i)*rupdir(i) + &
                     !                 (trntdr(i)-trndir(i))  &
                     !                 *rupdif(i))*refk
                     ! dir tran plus total diff trans times interface scattering plus
                     ! dir tran with up dir ref and down dif ref times interface scattering
                     ! fdirdn(i) = trndir(i) + (trntdr(i) &
                     !               - trndir(i) + trndir(i)  &
                     !               *rupdir(i)*rdndif(i))*refk
                     ! diffuse tran ref from below times interface scattering
                     ! fdifup(i) = trndif(i)*rupdif(i)*refk
                     ! diffuse tran times interface scattering
                     ! fdifdn(i) = trndif(i)*refk

                     ! netflux, down - up
                     ! dfdir = fdirdn - fdirup
                     dfdir(i) = trndir(i) &
                                 + (trntdr(i)-trndir(i)) * (c1 - rupdif(i)) * refk &
                                 -  trndir(i)*rupdir(i)  * (c1 - rdndif(i)) * refk
                     if (dfdir(i) < puny) dfdir(i) = c0
                     ! dfdif = fdifdn - fdifup
                     dfdif(i) = trndif(i) * (c1 - rupdif(i)) * refk
                     if (dfdif(i) < puny) dfdif(i) = c0
                  enddo       ! k

                  ! SNICAR_AD_RT_nonmodal is called twice for direct and diffuse incident fluxes
                  ! direct incident
                  if (flg_slr_in == 1) then
                    albedo = rupdir(snl_top)
                    dftmp  = dfdir
                    refk   = c1/(c1 - rdndif(snl_top)*rupdif(snl_top))
                    F_sfc_pls = (trndir(snl_top)*rupdir(snl_top) + &
                                (trntdr(snl_top)-trndir(snl_top))  &
                                 *rupdif(snl_top))*refk
                  !diffuse incident
                  else
                    albedo = rupdif(snl_top)
                    dftmp  = dfdif
                    refk   = c1/(c1 - rdndif(snl_top)*rupdif(snl_top))
                    F_sfc_pls = trndif(snl_top)*rupdif(snl_top)*refk
                  endif

                  ! Absorbed flux in each layer
                  do i=snl_top,snl_btm,1
                    F_abs(i) = dftmp(i)-dftmp(i+1)
                    flx_abs_b(i) = F_abs(i)
                    ! ERROR check: negative absorption
                    if (flx_abs_b(i) < -0.00001) then
                       if (err_w(fc,bnd_idx,k) == 0) then
                          err_w(fc,bnd_idx,k) = 2
                          errv_w(fc,bnd_idx,k) = flx_abs_b(i)
                       endif
                    endif
                  enddo

                  ! absobed flux by the underlying ground
                  F_btm_net = dftmp(snl_btm_itf)

                  ! note here, snl_btm_itf = 1 by snow column set up in CLM
                  flx_abs_b(1) = F_btm_net

                 if (flg_nosnl == 1) then
                    ! If there are no snow layers (but still snow), all absorbed energy must be in top soil layer
                    !flx_abs_b(:) = 0._r8
                    !flx_abs_b(1) = F_abs(0) + F_btm_net

                    ! changed on 20070408:
                    ! OK to put absorbed energy in the fictitous snow layer because routine SurfaceRadiation
                    ! handles the case of no snow layers. Then, if a snow layer is addded between now and
                    ! SurfaceRadiation (called in CanopyHydrology), absorbed energy will be properly distributed.
                    flx_abs_b(0) = F_abs(0)
                    flx_abs_b(1) = F_btm_net
                 endif

                 !Underflow check (we've already tripped the error condition above)
                 do i=snl_top,1,1
                    if (flx_abs_b(i) < 0._r8) then
                       flx_abs_b(i) = 0._r8
                    endif
                 enddo

                 F_abs_sum = 0._r8
                 do i=snl_top,snl_btm,1
                    F_abs_sum = F_abs_sum + F_abs(i)
                 enddo

                !enddo !enddo while (flg_dover > 0)

                ! Energy conservation check:
                ! Incident direct+diffuse radiation equals (absorbed+bulk_transmitted+bulk_reflected)
                energy_sum = (mu_not*pi*flx_slrd_b) + flx_slri_b - (F_abs_sum + F_btm_net + F_sfc_pls)
                if (abs(energy_sum) > 0.00001_r8) then
                   if (err_w(fc,bnd_idx,k) == 0) then
                      err_w(fc,bnd_idx,k) = 3
                      errv_w(fc,bnd_idx,k) = energy_sum
                   endif
                endif
                albout_b = albedo
                ! Check that albedo is less than 1
                if (albout_b > 1.0) then
                   if (err_w(fc,bnd_idx,k) == 0) then
                      err_w(fc,bnd_idx,k) = 4
                      errv_w(fc,bnd_idx,k) = albout_b
                   endif
                endif


                alb_w(fc,bnd_idx,k) = albout_b
                if (k <= nabs) abs_w(fc,:,bnd_idx,k) = flx_abs_b(-nlevsno+1:1)
             endif  ! no error
          endif
       enddo
       enddo
       enddo

       ! Errors found above: the first in (case, column, band) order, then abort
       nbad = 0
       do k = 1, ncall
       do fc = 1,num_nourbanc
       do bnd_idx = 1,numrad_snw
          if (err_w(fc,bnd_idx,k) /= 0) nbad = nbad + 1
       enddo
       enddo
       enddo
       if (nbad > 0) then
          ! spel acc-sync: device data this host code uses
          call snicar_ad_report(num_nourbanc, filter_nourbanc, ncall, call_slr, err_w, errv_w, nbad)
          call endrun(msg="SNICAR_AD_RT_nonmodal: "//errmsg(__FILE__, __LINE__))
       endif

       ! Combine the bands into VIS and NIR, per case and column
       do k = 1, ncall
       do fc = 1,num_nourbanc
          flg_slr_in = call_slr(k)
          c_idx = filter_nourbanc(fc)
          if (k <= nabs) flx_abs(fc,:,:,k) = 0._r8
          h2osno_lcl = h2osno(c_idx)
          if ((coszen(c_idx) > 0._r8) .and. (h2osno_lcl > min_snw) ) then
             if (snl(c_idx) > -1) then
                snl_top = 0
                rds_top = nint_snw_rds_min
             else
                snl_top = snl(c_idx)+1
                rds_top = nint(snw_rds(c_idx,snl_top))
             endif
             mu_not = max(coszen(c_idx), cp01)

             ! Incident flux weighting parameters
             !  - sum of all VIS bands must equal 1
             !  - sum of all NIR bands must equal 1
             !
             ! Spectral bands (5-band case)
             !  Band 1: 0.3-0.7um (VIS)
             !  Band 2: 0.7-1.0um (NIR)
             !  Band 3: 1.0-1.2um (NIR)
             !  Band 4: 1.2-1.5um (NIR)
             !  Band 5: 1.5-5.0um (NIR)
             !
             ! The following weights are appropriate for surface-incident flux in a mid-latitude winter atmosphere
             !
             ! 3-band weights
             if (numrad_snw==3) then
                ! Direct:
                if (flg_slr_in == 1) then
                   flx_wgt(1) = 1._r8
                   flx_wgt(2) = 0.66628670195247_r8
                   flx_wgt(3) = 0.33371329804753_r8
                   ! Diffuse:
                elseif (flg_slr_in == 2) then
                   flx_wgt(1) = 1._r8
                   flx_wgt(2) = 0.77887652162877_r8
                   flx_wgt(3) = 0.22112347837123_r8
                endif

                ! 5-band weights
             elseif(numrad_snw==5) then
                ! Direct:
                if (flg_slr_in == 1) then
                  if (atm_type_index == atm_type_default) then
                     flx_wgt(1) = 1._r8
                     flx_wgt(2) = 0.49352158521175_r8
                     flx_wgt(3) = 0.18099494230665_r8
                     flx_wgt(4) = 0.12094898498813_r8
                     flx_wgt(5) = 0.20453448749347_r8
                  else                 
                     slr_zen = nint(acos(coszen(c_idx)) * 180._r8 / pi)
                     if (slr_zen>89) then
                        slr_zen = 89
                     endif
                     flx_wgt(1) = 1._r8
                     flx_wgt(2) = flx_wgt_dir(atm_type_index, slr_zen+1, 2)
                     flx_wgt(3) = flx_wgt_dir(atm_type_index, slr_zen+1, 3)
                     flx_wgt(4) = flx_wgt_dir(atm_type_index, slr_zen+1, 4)
                     flx_wgt(5) = flx_wgt_dir(atm_type_index, slr_zen+1, 5)  
                  endif
				  
                   ! Diffuse:
                elseif (flg_slr_in == 2) then
                   if  (atm_type_index == atm_type_default) then
                     flx_wgt(1) = 1._r8
                     flx_wgt(2) = 0.58581507618433_r8
                     flx_wgt(3) = 0.20156903770812_r8
                     flx_wgt(4) = 0.10917889346386_r8
                     flx_wgt(5) = 0.10343699264369_r8
                  else
                     flx_wgt(1) = 1._r8
                     flx_wgt(2) = flx_wgt_dif(atm_type_index, 2)
                     flx_wgt(3) = flx_wgt_dif(atm_type_index, 3)
                     flx_wgt(4) = flx_wgt_dif(atm_type_index, 4)
                     flx_wgt(5) = flx_wgt_dif(atm_type_index, 5)
                  endif
                endif
             endif ! end if numrad_snw

             ! Weight output NIR albedo appropriately
             albout(fc,1,k) = alb_w(fc,1,k)
             flx_sum        = 0._r8
             do bnd_idx= nir_bnd_bgn,nir_bnd_end
                flx_sum = flx_sum + flx_wgt(bnd_idx)*alb_w(fc,bnd_idx,k)
             enddo
             albout(fc,2,k) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))

             ! Weight output NIR absorbed layer fluxes (flx_abs) appropriately
             if (k <= nabs) then
                flx_abs(fc,:,1,k) = abs_w(fc,:,1,k)
                do i=snl_top,1,1
                   flx_sum = 0._r8
                   do bnd_idx= nir_bnd_bgn,nir_bnd_end
                      flx_sum = flx_sum + flx_wgt(bnd_idx)*abs_w(fc,i,bnd_idx,k)
                   enddo
                   flx_abs(fc,i,2,k) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
                enddo
             endif

             ! near-IR direct albedo/absorption adjustment for high solar zenith angles
             ! solar zenith angle parameterization
             ! calculate the scaling factor for NIR direct albedo if SZA>75 degree
             if ((mu_not < mu_75) .and. (flg_slr_in == 1)) then
                sza_c1 = sza_a0 + sza_a1 * mu_not + sza_a2 * mu_not**2
                sza_c0 = sza_b0 + sza_b1 * mu_not + sza_b2 * mu_not**2
                sza_factor = sza_c1 * (log10(rds_top * c1) - c6) + sza_c0
                flx_sza_adjust  = albout(fc,2,k) * (sza_factor-c1) * sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
                albout(fc,2,k) = albout(fc,2,k) * sza_factor
                if (k <= nabs) flx_abs(fc,snl_top,2,k) = flx_abs(fc,snl_top,2,k) - flx_sza_adjust
             endif

             ! If snow < minimum_snow, but > 0, and there is sun, set albedo to underlying surface albedo
          elseif ( (coszen(c_idx) > 0._r8) .and. (h2osno_lcl < min_snw) .and. (h2osno_lcl > 0._r8) ) then
             albout(fc,1,k) = albsfc(c_idx,1)
             albout(fc,2,k) = albsfc(c_idx,2)

             ! There is either zero snow, or no sun
          else
             albout(fc,1,k) = 0._r8
             albout(fc,2,k) = 0._r8
          endif    ! if column has snow and coszen > 0
       enddo
       enddo


     end associate

   end subroutine SNICAR_AD_RT_nonmodal

   !-----------------------------------------------------------------------
   subroutine snicar_ad_report(num_nourbanc, filter_nourbanc, ncall, call_slr, err_w, errv_w, nbad)
     ! Host report of SNICAR_AD_RT_nonmodal's first error in (case, column, band) order.
     use elm_varpar, only : nlevsno
     integer , intent(in) :: num_nourbanc, filter_nourbanc(:), ncall, call_slr(:), nbad
     integer , intent(in) :: err_w(:,:,:)
     real(r8), intent(in) :: errv_w(:,:,:)
     character(len=40), parameter :: what(4) = (/ "snow grain radius out of bounds        ", &
          "negative absorption                    ", "energy conservation error              ", &
          "albedo > 1                             " /)
     integer :: k, fc, b, c
     do k = 1, ncall
        do fc = 1, num_nourbanc
           do b = 1, numrad_snw
              if (err_w(fc,b,k) /= 0) then
                 c = filter_nourbanc(fc)
                 write (iulog,*) "SNICAR ERROR: ", trim(what(err_w(fc,b,k))), ": ", errv_w(fc,b,k)
                 write (iulog,*) "SNICAR STATS: column ", c, " band ", b, " case ", k, " flg_slr ", call_slr(k), &
                      " (", nbad, " elements with errors)"
                 write (iulog,*) "SNICAR STATS: snl= ", col_pp%snl(c), " h2osno= ", col_ws%h2osno(c)
                 write (iulog,*) "SNICAR STATS: snw_rds= ", col_ws%snw_rds(c,-nlevsno+1:0)
                 return
              endif
           enddo
        enddo
     enddo
   end subroutine snicar_ad_report



 end module SnowSnicarNonModalMod
