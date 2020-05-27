module FanUpdateMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !
  ! This module interfaces the FAN (Flow of Agricultural Nitrogen) process model with
  ! ELM. This includes the main driver routine (fan_eval), initialization, and
  ! use-modifiable parameters, some of which are read from namelist. The remaining,
  ! internal parameters are set in FanMod.
  !
  ! FAN is implemented on column level. Synthetic fertilizer application is taken from the
  ! ELM crop model and remains in crop columns. Manure N is read from a stream (by
  ! fanStreamMod) which distinguishes pastoral and mixed/landless livestock
  ! systems. Manure in pastures is allocated to the native soil column. The mixed/landless
  ! systems are associated with crop columns, however, some N may be transferred to the
  ! native column due to manure spreading or seasonal grazing.
  !
  ! Within FAN, the nitrogen is distributed to several pools which represent different
  ! types of input (manures, fertilizers) and different "age" (time since
  ! fertilizer/manure application). The age determines properties like pH. The pools of
  ! same type but different age are called age classes in the FAN description paper. The
  ! model includes 4 slurry (manure) age classes, 3 grazing manure age classes, 2 urea age
  ! classes, 3 age classes for ammonium produced from urea, and 1 age class for non-urea
  ! NH4 fertilizer N.

  use FanMod
  use shr_kind_mod, only : r8 => shr_kind_r8, CL => shr_kind_cl
  use decompMod                       , only : bounds_type
  use atm2lndType                     , only : atm2lnd_type
  use TemperatureType                 , only : temperature_type
  use FrictionVelocityType            , only : frictionvel_type
  use shr_infnan_mod                  , only : isnan => shr_infnan_isnan
  use CNNitrogenStateType             , only : nitrogenstate_type
  use CNNitrogenFluxType              , only : nitrogenflux_type
!  use CNVegNitrogenFluxType           , only : cnveg_nitrogenflux_type
  use WaterStateType                  , only : waterstate_type
  use WaterFluxType                   , only : waterflux_type
  use SoilStateType                   , only : soilstate_type
  use ColumnType                      , only : col_pp                
  use VegetationType                  , only : veg_pp                
  use clm_varctl                      , only : iulog
  use subgridAveMod                   , only : p2c
  use VegetationDataType              , only : veg_nf, veg_ns, veg_es
  use ColumnDataType                  , only : col_nf, col_ns, col_ws, col_wf, col_es


  implicit none

  private

  ! !PUBLIC MEMBER FUNCTIONS:
!  public fan_readnml
  public fan_eval
!  public fan_to_sminn

  ! Structure of FAN TAN pools: number of age classes for N each type:
  integer, parameter :: num_cls_slr = 4 ! slurry (S0,S1,S2,S3)
  integer, parameter :: num_cls_grz = 3 ! grazing (G1, G2, G3)
  integer, parameter :: num_cls_urea = 2 ! urea before hydrolysis (U1, U2)
  integer, parameter :: num_cls_fert = 3 ! tan formed from urea (F1, F2, F3)
  integer, parameter :: num_cls_otherfert = 1 ! F4
  integer, parameter :: num_cls_max = 4 ! max of above

  ! Hydrogen ion concentration in TAN pools, mol/l == 10**-pH
  !
  ! Pastures and slurry. The last age class gets soil pH (from the FAN stream),
  ! so sizes
  ! are one less than the number of classes.
  real(r8), parameter :: Hconc_grz_def(num_cls_grz-1) = 10**(/-8.5_r8, -8.0_r8/)
  real(r8), parameter :: Hconc_slr_def(num_cls_slr-1) = 10**(/-8.0_r8, -8.0_r8, -8.0_r8/)
  ! Urea fertilizer. The other fertilizer (F4) pool gets soil pH.
  real(r8), parameter :: Hconc_fert(num_cls_fert) = 10**(/-7.0_r8, -8.5_r8, -8.0_r8/)

  ! Active layer thickness used by FAN. This is assumed to match the topmost CLM
  ! layer. If
  ! this is not the case, handling of the soil moisture becomes inconsistent. 
  real(r8), parameter :: dz_layer_fert = 0.02_r8 ! m, fertilizer
  real(r8), parameter :: dz_layer_grz = 0.02_r8  ! m, grazing
  real(r8), parameter :: dz_layer_slr = 0.02_r8  ! m, slurry

  ! Manure N composition
  real(r8) :: fract_tan = 0.6_r8 ! fraction of total ammoniacal nitrogen
  ! The following are fractions of the remaining non-TAN N:
  real(r8), parameter :: &
       fract_resist = 0.45_r8, & ! resistant organic N
       fract_unavail = 0.05_r8, & ! unvavailable organic N 
       fract_avail = 0.5_r8 ! available organic N

  ! application rate in meters water:
  real(r8), parameter :: water_init_grz = 0.006_r8 ! urine patch depth (m)
  real(r8), parameter :: depth_slurry = 0.005_r8   ! slurry application rate (m)
  real(r8), parameter :: water_init_fert = 1e-9_r8   ! water in fertilizer (assumed very little).
  real(r8), parameter :: cnc_nh3_air = 0.0_r8
  ! Slurry infiltration time
  real(r8), parameter :: slurry_infiltr_time = 6.0_r8*3600_r8 ! seconds
  ! Reduction factor for fertilizer due to mechanical incorporation.
  ! N available for volatilization becomes multiplied by (1-fert_incorp_reduct).
  real(r8) :: fert_incorp_reduct = 0.25_r8 !?BAD

  ! TAN pool age ranges (sec). 
  real(r8), parameter ::          &
       poolranges_grz(num_cls_grz) = (/24*3600.0_r8, 10*24*3600.0_r8, 360*24*3600.0_r8/), &
       poolranges_fert(num_cls_fert) = (/2.36_r8*24*3600.0_r8, 24*3600.0_r8, 360*24*3600.0_r8/), &
       poolranges_slr(num_cls_slr) = (/slurry_infiltr_time, 24*3600.0_r8, 10*24*3600.0_r8, 360*24*3600.0_r8/), &
       poolrange_otherfert(num_cls_otherfert) = (/360*24*3600.0_r8/) !?BAD

  ! soil pH for crops restricted between these limits: 
  real(r8), parameter :: pH_crop_min = 5.5_r8
  real(r8), parameter :: pH_crop_max = 7.5_r8

  ! Parameters for grazing in mixed/landless systems:
  real(r8), parameter :: tempr_min_grazing = 283.0_r8 ! Lowest 10-day daily-min temperature for grazing, K
  ! Fraction of ruminants grazing when permitted by temperature
  real(r8), parameter :: max_grazing_fract = 0.65_r8
  ! Normalization constants for barn and storage emissions.  The defaults are calibrated
  ! to roughly reproduce the EMEP emission factors under European climate.
  real(r8), parameter :: volat_coef_barns_open = 0.03_r8, volat_coef_barns_closed = 0.025_r8, volat_coef_stores = 0.025_r8

!  ! Fraction of manure N moved from crop to native columns (manure spreading)
!  real(r8) :: fract_spread_grass = 1.0_r8 ! BAD?
!
!  ! Fan coupling to soil BGC. Can be set on separately for crop and other
!  ! columns.
!  logical, public :: fan_to_bgc_crop = .false. !BAD remove?
!  logical, public :: fan_to_bgc_veg = .false.  !BAD remove?
!
!  ! Whether manure N in mixed/landless systems (manure_sgrz and manure_ngrz
!  ! streams) is
!  ! defined per crop or land area:
!  logical :: crop_man_is4crop_area = .true. !BAD?

contains

!  subroutine fan_readnml(NLFilename) !BAD could delete, I don't think we have these in elm
!    !
!    ! Read FAN namelist and set the module variables according to it. 
!    ! 
!    use spmdMod        , only : masterproc, mpicom
!    use fileutils      , only : getavu, relavu, opnfil
!    use clm_varctl     , only : use_fan
!    use shr_log_mod    , only : errMsg => shr_log_errMsg
!    use shr_nl_mod     , only : shr_nl_find_group_name
!    use abortutils     , only : endrun
!    use shr_mpi_mod    , only : shr_mpi_bcast
!    use FanStreamMod   , only : set_bcast_fanstream_pars
!    use FanMod         , only : nh4_ads_coef
!    character(len=*), intent(in) :: NLFilename ! Namelist filename
!
!    integer :: ierr                 ! error code
!    integer :: unitn                ! unit for namelist file
!    character(len=*), parameter :: subname = 'fan_readnml'
!    character(len=*), parameter :: nmlname = 'fan_nml'
!    integer :: stream_year_first_fan      ! first year in stream to use
!    integer :: stream_year_last_fan       ! last year in stream to use
!    integer :: model_year_align_fan       ! align stream_year_firstndep2 with 
!    character(len=CL)  :: stream_fldFileName_fan
!    character(len=CL)  :: fan_mapalgo
!
!    namelist /fan_nml/ fan_to_bgc_crop, fan_to_bgc_veg, stream_year_first_fan, &
!         stream_year_last_fan, model_year_align_fan, fan_mapalgo, stream_fldFileName_fan, &
!         fract_spread_grass, nh4_ads_coef
!
!    if (.not. use_fan) return
!
!    if (masterproc) then
!       unitn = getavu()
!       write(iulog, *) 'Read in ' // nmlname // '  namelist'
!       call opnfil(NLFilename, unitn, 'F')
!       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
!       if (ierr == 0) then
!          read(unitn, nml=fan_nml, iostat=ierr)
!          if (ierr /= 0) then
!             call endrun(msg="ERROR reading " // nmlname // "namelist" // errmsg(__FILE__, __LINE__))
!          end if
!       else
!          call endrun(msg="ERROR could NOT find " // nmlname // "namelist" // errmsg(__FILE__, __LINE__))
!       end if
!       call relavu(unitn)
!    end if
!
!    call set_bcast_fanstream_pars(stream_year_first_fan, stream_year_last_fan, &
!         model_year_align_fan, fan_mapalgo, stream_fldFileName_fan, crop_man_is4crop_area)
!
!    call shr_mpi_bcast(fan_to_bgc_crop, mpicom)
!    call shr_mpi_bcast(fan_to_bgc_veg, mpicom)
!    call shr_mpi_bcast(fract_spread_grass, mpicom)
!    call shr_mpi_bcast(nh4_ads_coef, mpicom)
!
!    if (nh4_ads_coef < 0) then
!       call endrun(msg="ERROR invalid nh4_ads_coef")
!    end if
!    if (fract_spread_grass > 1 .or. fract_spread_grass < 0) then
!       call endrun(msg="ERROR invalid fract_spread_grass")
!    end if
!
!  end subroutine fan_readnml

  !************************************************************************************

  subroutine fan_eval(bounds, num_soilc, filter_soilc, &
       atm2lnd_vars, &
!       cnveg_nitrogenflux_inst, &
       nitrogenflux_vars, &
       nitrogenstate_vars, &
       waterstate_vars, soilstate_vars, temperature_vars, &
       waterflux_vars, frictionvel_vars)
    use clm_time_manager, only: get_step_size, get_curr_date, get_curr_calday, get_nstep
    use clm_varpar, only: max_patch_per_col
    use LandunitType, only: lun_pp
    use shr_sys_mod, only : shr_sys_flush
    use GridcellType, only: grc_pp
    use abortutils, only : endrun
    use pftvarcon, only : nc4_grass, nc3_nonarctic_grass, nc3_arctic_grass
    use landunit_varcon, only:  istsoil, istcrop
    use clm_varcon, only : spval, ispval
    use decompMod, only : bounds_type
!    use subgridAveMod, only: p2c
!    use VegetationDataType, only: veg_nf, veg_ns
!    use ColumnDataType, only: col_nf, col_ns
    !
    ! Evaluate the N fluxes and update the state of the FAN pools. Uses the fertilization
    ! flux determined in CropPhenology; the manure inputs come from the FAN stream. The
    ! CLM soil N pools are not changed here but in fan_to_sminn.
    ! 
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
!    type(wateratm2lndbulk_type), intent(in)  :: wateratm2lndbulk_inst
!    type(cnveg_nitrogenflux_type)          , intent(in) :: cnveg_nitrogenflux_inst
    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type), intent(inout) :: nitrogenstate_vars
    type(waterstate_type)              , intent(in)    :: waterstate_vars
    type(soilstate_type)                   , intent(in)    :: soilstate_vars
    type(temperature_type)                 , intent(in) :: temperature_vars
    type(waterflux_type)               , intent(in)    :: waterflux_vars
    type(frictionvel_type)                 , intent(in) :: frictionvel_vars

    integer, parameter :: &
         ! Use this many sub-steps. This improves numerical accuracy but is perhaps not
         ! essential, because FAN includes an ad-hoc fixer for negative fluxes.
         num_substeps = 4, &
         ! FAN includes a separate nitrogen conservation check, which can be done every
         ! nth time step. This is mostly redundant, because the FAN pools and fluxes are
         ! now included in the main CLM soil N balance check. The FAN check has more
         ! detail.
         balance_check_freq = 1000
    ! Number of organic N types (available, unavailable, resistant)
    integer, parameter :: num_org_n_types = 3

    ! Organic (non-urea) manure N: production, pools, flux to TAN.  Named indices to
    ! denote each N fraction.
    real(r8) :: ndep_org(num_org_n_types) ! gN/m2/sec
    real(r8) :: orgpools(num_org_n_types) ! gN/m2
    real(r8) :: tanprod(num_org_n_types)  ! gN/m2/sec
    ! Temporary arrays for N fluxes and pools. Named indices to denote different fluxes.
    ! Numerical indices to denote the age class.
    real(r8) :: tanpools(num_cls_max)  ! gN/m2
    real(r8) :: ureapools(num_cls_urea) ! gN/m2
    real(r8) :: fluxes_tmp(num_fluxes) ! gN/m2/sec
    real(r8) :: fluxes(num_fluxes, num_cls_max) ! gN/m2/sec
    ! TAN production flux from urea. Include one extra flux for the residual flux out of U2.
    real(r8) :: tanprod_from_urea(num_cls_urea + 1) ! gN/m2/sec
    ! timestep-cumulative (gN/m2) aging flux out of the oldest class
    real(r8) :: n_residual, n_residual_total, urea_resid
    ! H ion concentrations for grz and slr pools (== prescribed + soil pH for oldest class)
    real(r8) :: Hconc_grz(num_cls_grz)
    real(r8) :: Hconc_slr(num_cls_slr)

    real(r8) :: dt ! timestep, sec
    real(r8) :: watertend ! time derivative of soil water content, m/s
    real(r8) :: ratm      ! total resistance Ra + Rb s/m
    real(r8) :: pH_soil, pH_crop ! background soil pH in native veg. and crops
    real(r8) :: tg ! soil temprature, K
    real(r8) :: theta ! volumetric soil moisture, m3/m3
    real(r8) :: thetasat ! volumetric soil moisture at saturation (porosity), m3/m3
    real(r8) :: infiltr_m_s, evap_m_s, runoff_m_s ! infiltration, evaporation, runoff, m/s
    real(r8) :: soilpsi ! soil matric potential, MPa
    real(r8) :: org_n_tot ! organic N in slurry application, gN/m2/sec
    real(r8) :: fert_total, fert_urea, fert_tan, fert_no3, soilflux_org
    real(r8) :: ngrz ! manure N flux in grazing, gN/m2/sec 
    real(r8) :: fert_inc_tan ! urea and NH4 fertilizer incorporated directly to soil, gN/m2/sec
    real(r8) :: fert_generic ! non-urea, non-no3 fertilizer N applied, gN/m2/sec

    ! index and auxiliary variables
    real(r8) :: soilph_min, soilph_max, bsw
    real(r8) :: nsoilman_old, nsoilfert_old, nstored_old
    integer :: def_ph_count
    integer :: c, g, p, l, fc, ind_substep, patchcounter, status
    logical :: do_balance_checks

    ! Set the constant pHs. The column-dependent pHs will be set below.
    Hconc_grz(1:num_cls_grz-1) = Hconc_grz_def
    Hconc_slr(1:num_cls_slr-1) = Hconc_slr_def

    soilph_min = 999
    soilph_max = -999
    def_ph_count = 0
    dt = real(get_step_size(), r8)
    do_balance_checks = balance_check_freq > 0 .and. mod(get_nstep(), balance_check_freq) == 0

    associate(&
         ns => nitrogenstate_vars, &
         nf => nitrogenflux_vars, &
         forc_ndep_urea => atm2lnd_vars%forc_ndep_urea_grc, & ! Fraction of urea in fertilizer N
         forc_ndep_no3 => atm2lnd_vars%forc_ndep_nitr_grc, &   ! Fraction of NO3 in fertilizer N
         ram1 => frictionvel_vars%ram1_patch, & ! Aerodynamic resistance, s/m
         rb1 => frictionvel_vars%rb1_patch)     ! Quasi-laminar layer resistance, s/m

    ! Convert the patch fertilizer application (from crop model) to column flux:
    call p2c(bounds, num_soilc, filter_soilc, &
         veg_nf%fert(bounds%begp:bounds%endp), &
         col_nf%fert_n_appl(bounds%begc:bounds%endc))

    if (do_balance_checks) then
       nstored_old = get_total_n(ns, nf, 'pools_storage')
       nsoilman_old = get_total_n(ns, nf, 'pools_manure')
       nsoilfert_old = get_total_n(ns, nf, 'pools_fertilizer')
    end if

    ! Assign the "pastoral" manure entirely to the natural vegetation column
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       l = col_pp%landunit(c)
       if (.not. col_pp%active(c) .or. col_pp%wtgcell(c) < 1e-6) cycle

       g = col_pp%gridcell(c)
       if (lun_pp%itype(l) == istsoil) then
          col_nf%manure_n_grz(c) &
               = atm2lnd_vars%forc_ndep3_grc(g) / col_pp%wtgcell(c) * 1e3 ! kg to g 
          if (col_nf%manure_n_appl(c) > 0) then
             write(iulog, *) 'manure_n_appl:', col_nf%manure_n_appl(c)
             call endrun(msg='Found fertilizer in soil column')
          end if
       else
          col_nf%manure_n_grz(c) = 0.0
       end if
    end do

    call handle_storage(bounds, temperature_vars, frictionvel_vars, dt, &
!         atm2lnd_vars%forc_ndep_sgrz_grc, atm2lnd_vars%forc_ndep_ngrz_grc, & !BAD do we need both of these?
         atm2lnd_vars%forc_ndep2_grc, & ! BAD substituted this instead
         col_ns%manure_n_stored, col_ns%manure_tan_stored, &
         col_nf%manure_n_appl, col_nf%manure_tan_appl, &
         col_nf%manure_n_grz, col_nf%manure_n_mix, &
         col_nf%nh3_stores, col_nf%nh3_barns, &
         col_nf%manure_n_transf, col_ns%fan_grz_fract, &
         col_nf%manure_n_barns, &
         fract_tan, &
         filter_soilc, num_soilc)

    if (debug_fan) then
       if (any(isnan(col_nf%nh3_stores))) then
          call endrun('nan nh3 stores')
       end if
       if (any(isnan(col_nf%nh3_barns))) then
          call endrun('nan nh3 barns')
       end if
       if (any(isnan(col_nf%manure_n_appl))) then
          call endrun('nan nh3 appl')
       end if
       if (any(isnan(col_nf%manure_n_mix))) then
          call endrun('nan nh3 appl')
       end if
    end if

    do fc = 1, num_soilc
       c = filter_soilc(fc)
       l = col_pp%landunit(c)
       g = col_pp%gridcell(c)
       if (.not. (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop)) cycle
       if (.not. col_pp%active(c) .or. col_pp%wtgcell(c) < 1e-15) cycle

       ! Find and average the atmospheric resistances Rb and Ra.
       ! 
       if (lun_pp%itype(col_pp%landunit(c)) == istcrop) then
          ! Crop column, only one patch
          p = col_pp%pfti(c)
          if (p /= col_pp%pftf(c)) call endrun(msg='Strange patch for crop')
          ratm = ram1(p) + rb1(p)
       else
          ! if natural, find average over grasses
          ratm = 0.0
          patchcounter = 0
          do p = col_pp%pfti(c), col_pp%pftf(c)
             if (      veg_pp%itype(p) == nc4_grass &
                  .or. veg_pp%itype(p) == nc3_nonarctic_grass &
                  .or. veg_pp%itype(p) == nc3_arctic_grass) then
                if (.not. veg_pp%active(p) .or. ram1(p) == spval .or. rb1(p) == spval) cycle
                ratm = ratm + ram1(p) + rb1(p)
                patchcounter = patchcounter + 1
             end if
          end do
          if (patchcounter > 0) then
             ratm = ratm / patchcounter
          else
             ! grass not found, take average over everything
             do p = col_pp%pfti(c), col_pp%pftf(c)
                if (.not. veg_pp%active(p) .or. ram1(p) == spval .or. rb1(p) == spval) cycle
                ratm = ratm + ram1(p) + rb1(p)
                patchcounter = patchcounter + 1
             end do
             if (patchcounter == 0) then
                call endrun(msg='Could not find any useful pft for ram1')
             end if
             ratm = ratm / patchcounter
          end if
          col_ns%fan_grz_fract(c) = 1.0_r8 ! for crops handled by handle_storage
       end if

!       watertend = waterstate_vars%h2osoi_tend_tsl(c) * 1e-3 ! to m/s
       watertend = 0.0_r8

       tg = col_es%t_grnd(c)
       theta = col_ws%h2osoi_vol(c,1)
       thetasat = soilstate_vars%watsat_col(c,1)
       bsw = soilstate_vars%bsw_col(c,1) !BAD?
       theta = min(theta, 0.98_r8*thetasat)
       infiltr_m_s = max(col_wf%qflx_infl(c), 0.0) * 1e-3 
       evap_m_s = col_wf%qflx_evap_grnd(c) * 1e-3
       runoff_m_s = max(col_wf%qflx_surf(c), 0.0) * 1e-3
       if (runoff_m_s > 1.0_r8) then
          runoff_m_s = 0.0_r8
       end if
       soilpsi = soilstate_vars%soilpsi_col(c,1) !BAD?

       ! grazing
       !
       ngrz = col_nf%manure_n_grz(c)
       ndep_org(ind_avail) = ngrz * (1.0_r8-fract_tan) * fract_avail
       ndep_org(ind_resist) = ngrz * (1.0_r8-fract_tan) * fract_resist
       ndep_org(ind_unavail) = ngrz * (1.0_r8-fract_tan) * fract_unavail

       orgpools(ind_avail) = col_ns%manure_a_grz(c)
       orgpools(ind_resist) = col_ns%manure_r_grz(c)
       orgpools(ind_unavail) = col_ns%manure_u_grz(c)
       call update_org_n(ndep_org, tg, orgpools, dt, tanprod, soilflux_org)
!       call update_org_n(ndep_org, tg, soilpsi, orgpools, dt, dz_layer_grz, &
!            tanprod, soilflux_org, status)
       col_ns%manure_a_grz(c) = orgpools(ind_avail)
       col_ns%manure_r_grz(c) = orgpools(ind_resist) 
       col_ns%manure_u_grz(c) = orgpools(ind_unavail)

       tanpools(1) = col_ns%tan_g1(c)
       tanpools(2) = col_ns%tan_g2(c)
       tanpools(3) = col_ns%tan_g3(c)

       ph_soil = atm2lnd_vars%forc_soilph_grc(g)
       if (ph_soil < 0.1) then
          ! Missing values in the pH, eg. Antarctica.
          ph_soil = 6.5_r8
          def_ph_count = def_ph_count + 1
       end if
       Hconc_grz(3) = 10**(-ph_soil)
       soilph_max = max(soilph_max, ph_soil)
       soilph_min = min(soilph_min, ph_soil)

       fluxes_tmp = 0.0
       n_residual_total = 0.0
       fluxes = 0.0
       n_residual = 0
       do ind_substep = 1, num_substeps
          call update_npool(tg, ratm, &
               theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_vars%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, ngrz * fract_tan, &
               (/0.0_r8, 0.0_r8, sum(tanprod)/), water_init_grz, &
               cnc_nh3_air, poolranges_grz, Hconc_grz, dz_layer_grz, tanpools(1:num_cls_grz), &
               fluxes(1:num_fluxes,1:num_cls_grz), &
               n_residual, dt/num_substeps, status, num_cls_grz)
!          call update_npool(tg, ratm, &
!               theta, thetasat, infiltr_m_s, evap_m_s, &
!               atm2lnd_vars%forc_q_downscaled(c), watertend, &
!               runoff_m_s, ngrz * fract_tan, &
!               (/0.0_r8, 0.0_r8, sum(tanprod)/), & ! all TAN procuced from org N goes to G3
!               water_init_grz, &
!               bsw, poolranges_grz, Hconc_grz, dz_layer_grz, tanpools(1:num_cls_grz), &
!               fluxes(1:num_fluxes,1:num_cls_grz), &
!               n_residual, dt/num_substeps, num_cls_grz, num_fluxes, status)
          if (status /= 0) then
             write(iulog, *) 'status = ', status, tanpools(1:num_cls_grz), &
                  ratm, theta, thetasat, ngrz * fract_tan, tanprod
             call endrun(msg='update_npool status /= 0')
          end if
          fluxes_tmp = fluxes_tmp + sum(fluxes(:,1:num_cls_grz), dim=2)
          n_residual_total = n_residual_total + n_residual
       end do
       fluxes_tmp = fluxes_tmp / num_substeps

       col_ns%tan_g1(c) = tanpools(1)
       col_ns%tan_g2(c) = tanpools(2)
       col_ns%tan_g3(c) = tanpools(3)

       col_nf%nh3_grz(c) = fluxes_tmp(iflx_air)
       col_nf%manure_nh4_runoff(c) = fluxes_tmp(iflx_roff)
       col_nf%manure_no3_to_soil(c) = fluxes_tmp(iflx_no3)
       col_nf%manure_nh4_to_soil(c) &
            = fluxes_tmp(iflx_soild) + fluxes_tmp(iflx_soilq) &
              + n_residual_total / dt + soilflux_org

       ! Manure application

       org_n_tot = col_nf%manure_n_appl(c) - col_nf%manure_tan_appl(c)
       ! Use the the same fractionation of organic N as for grazing, after removing the
       ! "explicitly" calculated TAN.
       ndep_org(ind_avail) = org_n_tot * fract_avail
       ndep_org(ind_resist) = org_n_tot * fract_resist 
       ndep_org(ind_unavail) = org_n_tot * fract_unavail 

       orgpools(ind_avail) = col_ns%manure_a_app(c)
       orgpools(ind_resist) = col_ns%manure_r_app(c)
       orgpools(ind_unavail) = col_ns%manure_u_app(c)
       call update_org_n(ndep_org, tg, orgpools, dt, tanprod, soilflux_org)
!       call update_org_n(ndep_org, tg, soilpsi, orgpools, dt, dz_layer_slr, &
!            tanprod, soilflux_org, status)
       col_ns%manure_a_app(c) = orgpools(ind_avail)
       col_ns%manure_r_app(c) = orgpools(ind_resist)
       col_ns%manure_u_app(c) = orgpools(ind_unavail)

       tanpools(1) = col_ns%tan_s0(c)
       tanpools(2) = col_ns%tan_s1(c)
       tanpools(3) = col_ns%tan_s2(c)
       tanpools(4) = col_ns%tan_s3(c)

       ph_crop = min(max(ph_soil, ph_crop_min), ph_crop_max)
       Hconc_slr(4) = 10**(-ph_crop)

       fluxes_tmp = 0.0
       n_residual_total = 0.0
       fluxes = 0.0
       do ind_substep = 1, num_substeps
          call update_4pool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_vars%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, col_nf%manure_tan_appl(c), sum(tanprod), cnc_nh3_air, depth_slurry, &
               poolranges_slr, tanpools(1:num_cls_slr), Hconc_slr, &
               fluxes(1:num_fluxes, 1:num_cls_slr), &
               n_residual, dt / num_substeps, status)
!          call update_4pool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
!               atm2lnd_vars%forc_q_downscaled(c), watertend, &
!               runoff_m_s, col_nf%manure_tan_appl(c), sum(tanprod), bsw, depth_slurry, &
!               poolranges_slr, tanpools(1:num_cls_slr), Hconc_slr, &
!               fluxes(1:num_fluxes, 1:num_cls_slr), &
!               n_residual, dt / num_substeps, dz_layer_slr, status)
          if (status /= 0) then
             write(iulog, *) 'status and tanpools: ', status, tanpools(1:num_cls_slr)
             write(iulog, *) 'tg, ratm, theta, thetasat:', tg, ratm, theta, thetasat
             write(iulog, *) 'tanfluxes:', col_nf%manure_tan_appl(c), tanprod, fluxes(1:num_fluxes,1:num_cls_slr)
             call endrun(msg='update_4pool status /= 0')
          end if
          fluxes_tmp = fluxes_tmp + sum(fluxes(:,1:num_cls_slr), dim=2)
          n_residual_total = n_residual_total + n_residual
       end do
       fluxes_tmp = fluxes_tmp / num_substeps

       col_ns%tan_s0(c) = tanpools(1)
       col_ns%tan_s1(c) = tanpools(2)
       col_ns%tan_s2(c) = tanpools(3)
       col_ns%tan_s3(c) = tanpools(4)

       col_nf%nh3_manure_app(c) = fluxes_tmp(iflx_air)
       col_nf%manure_nh4_runoff(c) = col_nf%manure_nh4_runoff(c) + fluxes_tmp(iflx_roff)
       col_nf%manure_no3_to_soil(c) = col_nf%manure_no3_to_soil(c) + fluxes_tmp(iflx_no3)
       col_nf%manure_nh4_to_soil(c) &
            = col_nf%manure_nh4_to_soil(c) + fluxes_tmp(iflx_soild) + fluxes_tmp(iflx_soilq) &
            + n_residual_total / dt + soilflux_org

       ! Fertilizer
       ! split fertilizer N berween urea, no3 and the remaining other ammonium-N.
       fert_total = col_nf%fert_n_appl(c)

       ! A fraction of fertilizer N is made unavailable by mechanical incorporation, this
       ! will be added directly to the to-soil flux (tan) or no3 production (no3) below.
       fert_inc_tan = fert_total * fert_incorp_reduct * (1.0 - forc_ndep_no3(g))

       if (forc_ndep_urea(g) < 0 .or. forc_ndep_no3(g) < 0 .or. forc_ndep_urea(g) + forc_ndep_no3(g) > 1) then
          call endrun('bad fertilizer fractions')
       end if
       fert_urea = fert_total * forc_ndep_urea(g) * (1.0_r8 - fert_incorp_reduct)

       ! Fertilizer nitrate goes straight to the no3_prod, incorporated or not.
       fert_no3 = fert_total * forc_ndep_no3(g)
       fert_generic = fert_total * (1.0_r8 - forc_ndep_urea(g) - forc_ndep_no3(g)) &
                                 * (1.0_r8 - fert_incorp_reduct)
       ! Below also includes the incorporated N:
       col_nf%otherfert_n_appl(c) = fert_total * (1.0_r8 - forc_ndep_urea(g)) 

       ! Urea decomposition 
       ! 
       ureapools(1) = col_ns%fert_u1(c)
       ureapools(2) = col_ns%fert_u2(c)
       fluxes = 0.0
       call update_urea(tg, theta, thetasat, infiltr_m_s, evap_m_s, watertend, &
            runoff_m_s, fert_urea, ureapools, fluxes(1:num_fluxes,1:num_cls_urea), &
            urea_resid, poolranges_fert(1:num_cls_urea), &
            dt, status, 2)
!       call update_urea(tg, theta, thetasat, infiltr_m_s, evap_m_s, watertend, &
!            runoff_m_s, fert_urea, bsw, ureapools, fluxes(1:num_fluxes,1:num_cls_urea), &
!            urea_resid, poolranges_fert(1:num_cls_urea), &
!            dt, dz_layer_fert, status)
       if (status /= 0) then
          call endrun(msg='Bad status after update_urea for fertilizer')
       end if
       ! Nitrogen fluxes from the urea pool. Be sure to not zero below!
       fluxes_tmp = sum(fluxes(:,1:num_cls_urea), dim=2)

       col_ns%fert_u1(c) = ureapools(1)
       col_ns%fert_u2(c) = ureapools(2)
       ! Collect the formed ammonia for updating the TAN pools
       tanprod_from_urea(1:num_cls_urea) = fluxes(iflx_to_tan, 1:num_cls_urea)
       ! There is no urea pool corresponding to tan_f2, because most of the urea will
       ! have decomposed. Here whatever remains gets sent to tan_f2. 
       tanprod_from_urea(num_cls_urea+1) = urea_resid / dt 

       tanpools(1) = col_ns%tan_f1(c)
       tanpools(2) = col_ns%tan_f2(c)
       tanpools(3) = col_ns%tan_f3(c)         
       n_residual_total = 0.0
       fluxes = 0.0
       col_nf%nh3_otherfert(c) = 0.0
       do ind_substep = 1, num_substeps
          ! Fertilizer pools f1...f3
          call update_npool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_vars%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, 0.0_r8, tanprod_from_urea, water_init_fert, cnc_nh3_air, &
               poolranges_fert, Hconc_fert, dz_layer_fert, &
               tanpools(1:num_cls_fert), fluxes(1:num_fluxes,1:num_cls_fert), &
               n_residual, dt/num_substeps, status, num_cls_fert)
!          call update_npool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
!               atm2lnd_vars%forc_q_downscaled(c), watertend, &
!               runoff_m_s, 0.0_r8, tanprod_from_urea, water_init_fert, bsw, &
!               poolranges_fert, Hconc_fert, dz_layer_fert, &
!               tanpools(1:num_cls_fert), fluxes(1:num_fluxes,1:num_cls_fert), &
!               n_residual, dt/num_substeps, num_cls_fert, num_fluxes, status)
          if (status /= 0) then
             write(iulog, *) 'status:', status, tanpools(1:num_cls_fert), col_nf%fert_n_appl(c)
             call endrun(msg='Bad status after npool for fertilizer')
          end if
          fluxes_tmp = fluxes_tmp + sum(fluxes(:,1:num_cls_fert), dim=2) /num_substeps
          n_residual_total = n_residual_total + n_residual

          ! Fertilizer pool f4
          call update_npool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_vars%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, fert_generic, (/0.0_r8/), water_init_fert, cnc_nh3_air, &
               poolrange_otherfert, (/10**(-ph_crop)/), dz_layer_fert, &
               col_ns%tan_f4(c:c), fluxes(1:num_fluxes,1:1), &
               n_residual, dt/num_substeps, status, 1)
!          call update_npool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
!               atm2lnd_vars%forc_q_downscaled(c), watertend, &
!               runoff_m_s, fert_generic, (/0.0_r8/), water_init_fert, bsw, &
!               poolrange_otherfert, (/10**(-ph_crop)/), dz_layer_fert, &
!               col_ns%tan_f4(c:c), fluxes(1:num_fluxes,1:1), &
!               n_residual, dt/num_substeps, 1, num_fluxes, status)
          if (status /= 0) then
             write(iulog, *) 'status:', status, col_ns%tan_f4(c:c), col_nf%fert_n_appl(c)
             call endrun(msg='Bad status after npool for generic')
          end if
          fluxes_tmp = fluxes_tmp + fluxes(:, 1) / num_substeps
          n_residual_total = n_residual_total + n_residual
          col_nf%nh3_otherfert(c) = col_nf%nh3_otherfert(c) + fluxes(iflx_air, 1) / num_substeps
       end do

       col_ns%tan_f1(c) = tanpools(1)
       col_ns%tan_f2(c) = tanpools(2)
       col_ns%tan_f3(c) = tanpools(3)
       ! !!tan_f4 already updated above by update_npool!!

       col_nf%nh3_fert(c) = fluxes_tmp(iflx_air)
       col_nf%fert_nh4_runoff(c) = fluxes_tmp(iflx_roff)
       col_nf%fert_no3_to_soil(c) = fluxes_tmp(iflx_no3) + fert_no3
       col_nf%fert_nh4_to_soil(c) = fluxes_tmp(iflx_soild) + fluxes_tmp(iflx_soilq) &
            + n_residual_total/dt + fert_inc_tan

       ! Total flux
       ! 
       col_nf%nh3_total(c) = col_nf%nh3_fert(c) + col_nf%nh3_manure_app(c) &
            + col_nf%nh3_grz(c) + col_nf%nh3_stores(c) +  col_nf%nh3_barns(c)
       if (col_nf%nh3_total(c) < -1e15) then
          call endrun(msg='ERROR: FAN, negative total emission')
       end if
    end do

    if (do_balance_checks) then
       call balance_check('Storage', nstored_old, &
            get_total_n(ns, nf, 'pools_storage'), get_total_n(ns, nf, 'fluxes_storage'))
       call balance_check('Manure', nsoilman_old, &
            get_total_n(ns, nf, 'pools_manure'), get_total_n(ns, nf, 'fluxes_manure'))
       call balance_check('Fertilizer', nsoilfert_old, &
            get_total_n(ns, nf, 'pools_fertilizer'), get_total_n(ns, nf, 'fluxes_fertilizer'))
    end if

!    call update_summary(ns, nf, filter_soilc, num_soilc)

    end associate

  contains

    real(r8) function get_total_n(ns, nf, which) result(total)
      type(nitrogenstate_type), intent(in) :: ns
      type(nitrogenflux_type), intent(in) :: nf
      character(len=*), intent(in) :: which

      total = 0

      associate(soilc => filter_soilc(1:num_soilc))

      select case(which)
      case('pools_storage')
    !     total = 0.0_r8 ! no explicit storage in this version
         total = sum(col_ns%manure_n_stored(soilc))

      case('fluxes_storage')
         total = sum(col_nf%manure_n_mix(soilc))
         total = total - sum(col_nf%nh3_stores(soilc))
         total = total - sum(col_nf%nh3_barns(soilc)) - sum(col_nf%manure_n_transf(soilc))

      case('pools_manure')
         total = total + sum(col_ns%tan_g1(soilc)) + sum(col_ns%tan_g2(soilc)) + sum(col_ns%tan_g3(soilc)) 
         total = total + sum(col_ns%manure_u_grz(soilc)) &
              + sum(col_ns%manure_a_grz(soilc)) + sum(col_ns%manure_r_grz(soilc))
         total = total + sum(col_ns%tan_s0(soilc)) &
              + sum(col_ns%tan_s1(soilc)) + sum(col_ns%tan_s2(soilc)) + sum(col_ns%tan_s3(soilc))
         total = total + sum(col_ns%manure_u_app(soilc)) &
              + sum(col_ns%manure_a_app(soilc)) + sum(col_ns%manure_r_app(soilc))

      case('fluxes_manure')
         total = sum(col_nf%manure_n_grz(soilc)) + sum(col_nf%manure_n_barns(soilc)) 
         total = total - sum(col_nf%nh3_manure_app(soilc)) &
              - sum(col_nf%nh3_grz(soilc)) - sum(col_nf%manure_nh4_runoff(soilc))
         total = total - sum(col_nf%manure_no3_to_soil(soilc)) - sum(col_nf%manure_nh4_to_soil(soilc))
         total = total - sum(col_nf%manure_n_transf(soilc)) - sum(col_nf%nh3_stores(soilc)) - sum(col_nf%nh3_barns(soilc))

      case('pools_fertilizer')
         total = sum(col_ns%tan_f1((soilc))) + sum(col_ns%tan_f2((soilc))) + sum(col_ns%tan_f3(soilc)) &
              + sum(col_ns%tan_f4(soilc))
         total = total + sum(col_ns%fert_u1(soilc)) + sum(col_ns%fert_u2(soilc))

      case('fluxes_fertilizer')
         total = sum(col_nf%fert_n_appl(soilc))
         total = total - sum(col_nf%nh3_fert(soilc)) - sum(col_nf%fert_nh4_runoff(soilc))
         total = total - sum(col_nf%fert_no3_to_soil(soilc)) - sum(col_nf%fert_nh4_to_soil(soilc))

      case default
         call endrun(msg='Bad argument to get_total_n')

      end select

      end associate

    end function get_total_n

    subroutine balance_check(label, total_old, total_new, flux)
      ! Check and report that the net flux equals the accumulated mass in pools. The
      ! total pools and fluxes can be evaluated by the function get_total_n.
      character(len=*), intent(in) :: label
      real(r8), intent(in) :: total_old, total_new, flux

      real(r8) :: diff, accflux
      real(r8) :: tol = 1e-6_r8

      diff = total_new - total_old
      accflux = flux*dt

    end subroutine balance_check

  end subroutine fan_eval

  !************************************************************************************

  subroutine handle_storage(bounds, temperature_vars, frictionvel_vars, dt,  &
       ndep_mixed_grc, n_stored, tan_stored, &
!       ndep_sgrz_grc, ndep_ngrz_grc, &
       n_manure_spread, tan_manure_spread, &
       n_manure_graze, n_manure_mixed, &
       nh3_flux_stores, nh3_flux_barns, manure_n_transf, &
       grz_fract, manure_n_barns, tan_fract_excr, &
       filter_soilc, num_soilc)
    !
    ! Evaluate storage losse for manure in the mixed/landless production systems
    ! associated with crop columns in CLM. The N remaining after storage losses is applied
    ! on soil either in the same column, or a fraction may be moved to the native
    ! vegetation column within the gridcell (grasslands). This subroutine also evaluates
    ! seasonal grazing of livestock, and the manure N on pastures is also moved into the
    ! native vegatation.
    ! 
    use landunit_varcon, only : max_lunit
    use pftvarcon, only : nc4_grass, nc3_nonarctic_grass, nc3_arctic_grass
    use clm_varcon, only : ispval
    use landunit_varcon,      only:  istsoil, istcrop
    use abortutils     , only : endrun
    use LandunitType   , only: lun_pp
    use GridcellType   , only: grc_pp
    use clm_varctl     , only : iulog
    use VegetationType           , only : veg_pp

    implicit none
    type(bounds_type), intent(in)    :: bounds
    type(temperature_type) , intent(in) :: temperature_vars
    type(frictionvel_type) , intent(in) :: frictionvel_vars
    real(r8), intent(in) :: dt ! timestep, sec

    ! N excreted in manure, gN/m2:
!    real(r8), intent(in) :: ndep_sgrz_grc(bounds%begg:) ! seasonally grazing animals
!    real(r8), intent(in) :: ndep_ngrz_grc(bounds%begg:) ! non-grazing animals
    real(r8), intent(in) :: ndep_mixed_grc(bounds%begg:bounds%endg)
    real(r8), intent(inout) :: n_stored(bounds%begc:bounds%endc), tan_stored(bounds%begc:bounds%endc) ! N, TAN currently stored, gN/m2
    ! N, TAN spread on grasslands, gN/m2/s:
    real(r8), intent(inout) :: n_manure_spread(bounds%begc:) 
    real(r8), intent(inout) :: tan_manure_spread(bounds%begc:) ! output, calculated from the above and stored manure
    ! N excreted by animals allocated to mixed production systems temporarily
    ! grazing on grasslands:
    real(r8), intent(inout) :: n_manure_graze(bounds%begc:)
    ! N excreted by animals in mixed systems, total
    real(r8), intent(inout) :: n_manure_mixed(bounds%begc:)
    ! NH3 emission fluxes from manure storage and housings, gN/m2/s
    real(r8), intent(inout) :: nh3_flux_stores(bounds%begc:), nh3_flux_barns(bounds%begc:)
    ! total nitrogen flux transferred out of a crop column (manure spreading + temporary grazing)
    real(r8), intent(inout) :: manure_n_transf(bounds%begc:)
    ! Total nitrogen excreted in barns
    real(r8), intent(inout) :: manure_n_barns(bounds%begc:)
    ! fraction of manure excreted when grazing
    real(r8), intent(inout) :: grz_fract(bounds%begc:)
    ! TAN fraction in excreted N
    real(r8), intent(in) :: tan_fract_excr
    integer, intent(in)    :: num_soilc       ! number of soil columns in filter
    integer, intent(in)    :: filter_soilc(:) ! filter for soil columns

    real(r8), parameter :: kg_to_g = 1e3_r8
    ! FAN allows a fraction of manure N diverted before storage; this option is currently
    ! not used.
    real(r8), parameter :: fract_direct = 0.0_r8

    ! N fluxes, gN/m2/sec:
    !
!    real(r8) :: flux_avail_rum ! Total ruminant manure in barns (not grazing)
!    real(r8) :: flux_avail_mg  ! Non-ruminant manure (only barns)
    real(r8) :: flux_avail
    real(r8) :: flux_grazing   ! Ruminant manure, grazing in mixed/landless systems Manure
    ! Manure N collected from crop columns and moved to nat. veg. column within the gridcell;
    ! grazing:
    real(r8) :: flux_grass_graze
    ! Manure N collected from crop columns and moved to nat. veg. column within the gridcell;
    ! manure application:
    real(r8) :: flux_grass_spread ! Organic + TAN
    real(r8) :: flux_grass_spread_tan ! TAN only
    real(r8) :: total_to_store ! total N remaining after losses in barns
    real(r8) :: total_to_store_tan ! TAN remaining after losses in barns

    ! scaling factor for converting from per-land-area to per-crop area if needed:
!    real(r8) :: invscale

    ! N fluxes evaluated by eval_fluxes_storage. Indices are in FanMode.F9.  Dimensions
    ! are (type of flux).
!    real(r8) :: fluxes_nitr(num_fluxes,2), fluxes_tan(num_fluxes,2) 
    real(r8) :: fluxes_nitr(num_fluxes), fluxes_tan(num_fluxes)
    ! Auxiliary and index variables:
    logical :: is_grass
    integer :: begg, endg, g, l, c, il, counter, col_grass, status, p, fc
    real(r8) :: cumflux, totalinput 

    begg = bounds%begg; endg = bounds%endg

    associate(&
      t_ref2m => veg_es%t_ref2m, & ! 2m temperature, K
      u10 => frictionvel_vars%u10_patch, & ! 10m wind speed, m/s
      t_a10min => veg_es%t_a10min) ! 10-day running mean 2m temperature (K)

    totalinput = 0.0
    cumflux = 0.0

    do fc = 1, num_soilc
       ! Zero the manure N spread arrays because the nat veg columns receive N additively (see below)
       c = filter_soilc(fc)
       tan_manure_spread(c) = 0.0_r8
       n_manure_spread(c) = 0.0_r8
    end do

    do g = begg, endg
       ! First find out if there are grasslands in this cell. If yes, a fraction of manure
       ! can be diverted to them before storage. At this time, we assume that the
       ! grasslands are always found in the natural vegetation column.
       col_grass = ispval
       l = grc_pp%landunit_indices(istsoil, g)
       do c = lun_pp%coli(l), lun_pp%colf(l)
          if (col_pp%itype(c) == istsoil) then
             col_grass = c
             exit
          end if
       end do
       if (col_grass == ispval) then
          call endrun(msg='Failed to find net veg column')
       end if

       ! Transfer of manure from all crop columns to the natural vegetation column:
       flux_grass_graze = 0.0_r8
       flux_grass_spread = 0.0_r8
       flux_grass_spread_tan = 0.0_r8

       l = grc_pp%landunit_indices(istcrop, g)
       if (l /= ispval) then
          ! flux_avail = manure excreted per m2 of crops (ndep_mixed_grc = per m2 / all land units)
          do c = lun_pp%coli(l), lun_pp%colf(l)
             if (.not. col_pp%active(c)) cycle

!             if (crop_man_is4crop_area) then
!                invscale = 1.0_r8
!             else
!                invscale = 1.0_r8 / lun_pp%wtgcell(l)
!             end if

!             n_manure_mixed(c) = (ndep_ngrz_grc(g) + ndep_sgrz_grc(g)) * kg_to_g * invscale
              flux_avail = ndep_mixed_grc(g) * kg_to_g / lun_pp%wtgcell(l)

             n_manure_mixed(c) = flux_avail
             totalinput = totalinput + flux_avail

             if (t_a10min(col_pp%pfti(c)) > tempr_min_grazing) then
                ! fraction of animals grazing -> allocate some manure to grasslands before barns
                flux_grazing = max_grazing_fract * flux_avail !ndep_sgrz_grc(g) * kg_to_g * invscale
                flux_avail = flux_avail - flux_grazing
!                flux_avail_rum = (ndep_sgrz_grc(g)*(1.0_r8 - max_grazing_fract)) * kg_to_g * invscale
                grz_fract(c) = max_grazing_fract
             else
                flux_grazing = 0.0_r8
!                flux_avail_rum = ndep_sgrz_grc(g) * kg_to_g * invscale
                grz_fract(c) = 0.0_r8
             end if
!             flux_avail_mg = ndep_ngrz_grc(g) * kg_to_g * invscale
             flux_grass_graze = flux_grass_graze + flux_grazing*col_pp%wtgcell(c)
!
!             totalinput = totalinput + flux_avail_rum + flux_avail_mg

             counter = 0
             if (col_grass == c) call endrun('Something wrong with the indices')
             if (col_pp%pfti(c) /= col_pp%pftf(c)) then
                call endrun(msg="ERROR crop column has multiple patches")
             end if

!             man_n_barns(c) = flux_avail_rum + flux_avail_mg
             manure_n_barns(c) = flux_avail

             ! Evaluate the NH3 losses, separate for ruminants (open barns) and others
             ! (poultry and pigs, closed barns). Note the slicing of fluxes(:,:) and
             ! fluxes_tan(:,:).
             nh3_flux_stores(c) = 0.0

!             if (flux_avail_rum < 0) then 
!                write(iulog, *) 'flux:', flux_avail_rum, ndep_sgrz_grc(g), (1.0_r8 - max_grazing_fract) * kg_to_g * invscale
!                call endrun(msg='negat flux_avail for ruminants')
!             end if

!             ! Ruminants
!             call eval_fluxes_storage(flux_avail_rum, 'open', &
!                  t_ref2m(col_pp%pfti(c)), u10(col_pp%pfti(c)), &
!                  fract_direct, volat_coef_barns_open, volat_coef_stores, &
!                  tan_fract_excr, fluxes_nitr(:,1), fluxes_tan(:,1), status)
!             if (status /=0) then 
!                write(iulog, *) 'status = ', status
!                call endrun(msg='eval_fluxes_storage failed for ruminants')
!             end if
!
!             ! Others
!             call eval_fluxes_storage(flux_avail_mg, 'closed', &
!                  t_ref2m(col_pp%pfti(c)), u10(col_pp%pfti(c)), &
!                  fract_direct, volat_coef_barns_closed, volat_coef_stores, &
!                  tan_fract_excr, fluxes_nitr(:,2), fluxes_tan(:,2), status)
!             if (status /=0) then 
!                write(iulog, *) 'status = ', status
!                call endrun(msg='eval_fluxes_storage failed for other livestock')
!             end if
             call eval_fluxes_storage(flux_avail, &
                  t_ref2m(col_pp%pfti(c)), u10(col_pp%pfti(c)), &
                  fract_direct, volat_coef_barns_closed, volat_coef_stores, &
                  tan_fract_excr, fluxes_nitr, fluxes_tan, status)
             if (status /=0) then 
                write(iulog, *) 'status = ', status
                call endrun(msg='eval_fluxes_storage failed')
             end if

             cumflux = cumflux + sum(fluxes_nitr)

             if (fluxes_tan(iflx_to_store) < 0) then
                call endrun(msg="ERROR too much manure lost")
             end if
             ! Manure storage is not evaluated explicitly, instead, the flux to storage
             ! will be spread "immediately".
!             total_to_store = sum(fluxes_nitr(iflx_to_store,:))
!             total_to_store_tan = sum(fluxes_tan(iflx_to_store,:))

!             n_manure_spread(c) = (1.0_r8 - fract_spread_grass) * total_to_store
!             tan_manure_spread(c) = (1.0_r8 - fract_spread_grass) * total_to_store_tan

!             flux_grass_spread = flux_grass_spread &
!                  + fract_spread_grass * total_to_store*col_pp%wtgcell(c)
!             flux_grass_spread_tan = flux_grass_spread_tan &
!                  + fract_spread_grass * total_to_store_tan*col_pp%wtgcell(c)
!
!             man_n_transf(c) = flux_grazing + fract_spread_grass*total_to_store
!
!             nh3_flux_stores(c) = sum(fluxes_nitr(iflx_air_stores,:))
!             nh3_flux_barns(c) = sum(fluxes_nitr(iflx_air_barns,:))

              flux_grass_spread = flux_grass_spread + fluxes_nitr(iflx_to_store)*col_pp%wtgcell(c)
              flux_grass_spread_tan = flux_grass_spread_tan + fluxes_tan(iflx_to_store)*col_pp%wtgcell(c)

              manure_n_transf(c) = flux_grazing + fluxes_nitr(iflx_to_store)

              nh3_flux_stores(c) = fluxes_nitr(iflx_air_stores)
              nh3_flux_barns(c) = fluxes_nitr(iflx_air_barns)

          end do ! column
       end if ! land unit not ispval

       if (col_grass /= ispval) then
          n_manure_spread(col_grass) = n_manure_spread(col_grass) &
               + flux_grass_spread / col_pp%wtgcell(col_grass)
          tan_manure_spread(col_grass) = tan_manure_spread(col_grass) &
               + flux_grass_spread_tan / col_pp%wtgcell(col_grass)
          n_manure_graze(col_grass) = n_manure_graze(col_grass) + flux_grass_graze / col_pp%wtgcell(col_grass)
          if (tan_manure_spread(col_grass) > 1) then
             ! In principle this could happen if col_pp%wtgcell(col_grass) is very small.
             write(iulog, *) 'Warning (FAN): suspicious manure N spread flux to natural vegetation column; flux, icol:', &
                  flux_grass_spread_tan, col_grass, tan_manure_spread(col_grass)
          end if
       else if (flux_grass_spread > 0) then
          ! There was no column that had a grass pft:
          write(iulog, *) 'Warning (FAN): fract_spread_grass > 0 not possible in this cell:', g
       end if

    end do ! grid
  end associate
  end subroutine handle_storage

  !************************************************************************************
!BAD commenting out for now, but this can go back in later when we get the 2-way coupling in

!  subroutine update_summary(ns, nf, filter_soilc, num_soilc)
!    ! Collect FAN fluxes and pools into aggregates used by the
!    ! NitrogenBalanceCheck.
!    use ColumnType, only : col
!    use LandunitType   , only: lun
!    use landunit_varcon, only : istcrop
!
!    type(soilbiogeochem_nitrogenstate_type), intent(inout) :: ns
!    type(soilbiogeochem_nitrogenflux_type), intent(inout) :: nf
!    integer, intent(in)    :: num_soilc       ! number of soil columns in filter
!    integer, intent(in)    :: filter_soilc(:) ! filter for soil columns
!
!    integer :: c, fc
!    real(r8) :: total, fluxout, fluxin, flux_loss    
!
!    do fc = 1, num_soilc
!       c = filter_soilc(fc)
!       total = col_ns%tan_g1(c) + col_ns%tan_g2(c) + col_ns%tan_g3(c)
!       total = total + col_ns%manure_u_grz(c) + col_ns%manure_a_grz(c) + col_ns%manure_r_grz(c)
!       total = total + col_ns%tan_s0(c) + col_ns%tan_s1(c) + col_ns%tan_s2(c) + col_ns%tan_s3(c)
!       total = total + col_ns%manure_u_app(c) + col_ns%manure_a_app(c) + col_ns%manure_r_app(c)
!       total = total + col_ns%tan_f1(c) + col_ns%tan_f2(c) + col_ns%tan_f3(c) + col_ns%tan_f4(c)
!       total = total + col_ns%fert_u1(c) + col_ns%fert_u2(c)
!       col_ns%fan_totn(c) = total
!
!       if (lun_pp%itype(col_pp%landunit(c)) == istcrop) then
!          ! no grazing, manure_n_appl is from the same column and not counted as
!          ! input
!          fluxin = col_nf%manure_n_mix(c) + col_nf%fert_n_appl(c)
!       else
!          ! no barns or fertilization. manure_n_appl is transferred from crop
!          ! columns and not
!          ! included in the other inputs.
!          fluxin = col_nf%manure_n_grz(c) + col_nf%manure_n_appl(c)
!       end if
!
!       flux_loss = col_nf%nh3_manure_app(c) + col_nf%nh3_grz(c) + col_nf%manure_nh4_runoff(c) &
!            + col_nf%nh3_stores(c) + col_nf%nh3_barns(c) &
!            + col_nf%nh3_fert(c) + col_nf%fert_nh4_runoff(c)
!       fluxout = col_nf%fert_no3_to_soil(c) + col_nf%fert_nh4_to_soil(c) &
!            + col_nf%manure_no3_to_soil(c) + col_nf%manure_nh4_to_soil(c) &
!            + col_nf%manure_n_transf(c) + flux_loss
!
!       col_nf%fan_totnin(c) = fluxin
!       col_nf%fan_totnout(c) = fluxout
!       col_nf%manure_n_total(c) = col_nf%manure_n_grz(c) + col_nf%manure_n_barns(c)
!    end do
!
!  end subroutine update_summary

  !************************************************************************************

!  subroutine fan_to_sminn(bounds, filter_soilc, num_soilc, sbgc_nf, nfertilization_patch)
!    !
!    ! Collect the FAN fluxes into totals which are either passed to the CLM N cycle
!    ! (depending on the fan_to_bgc_ switches) or used diagnostically.
!    !
!    use ColumnType, only : col
!    use LandunitType   , only: lun
!    use landunit_varcon, only : istcrop, istsoil
!
!    type(bounds_type), intent(in) :: bounds
!    integer, intent(in) :: filter_soilc(:)
!    integer, intent(in) :: num_soilc
!    type(soilbiogeochem_nitrogenflux_type), intent(inout) :: sbgc_nf
!    ! patch level fertilizer application + manure production 
!    real(r8), intent(inout) :: nfertilization_patch(bounds%begp:) 
!
!    integer :: c, fc, p
!    real(r8) :: flux_manure, flux_fert, manure_prod
!    logical :: included
!
!    if (.not. (fan_to_bgc_veg .or. fan_to_bgc_crop)) return
!
!    do fc = 1, num_soilc
!       c = filter_soilc(fc)
!       flux_manure = sbgc_nf%manure_no3_to_soil(c) + sbgc_nf%manure_nh4_to_soil(c)
!       flux_fert = sbgc_nf%fert_no3_to_soil(c) + sbgc_nf%fert_nh4_to_soil(c)
!       manure_prod = sbgc_nf%manure_n_barns(c) + sbgc_nf%manure_n_grz(c)
!
!       included = (lun_pp%itype(col_pp%landunit(c)) == istcrop .and. fan_to_bgc_crop) &
!             .or. (lun_pp%itype(col_pp%landunit(c)) == istsoil .and. fan_to_bgc_veg)
!
!       if (included) then
!          sbgc_nf%fert_to_sminn(c) = flux_fert + flux_manure
!          sbgc_nf%manure_n_to_sminn(c) = flux_manure
!          sbgc_nf%synthfert_n_to_sminn(c) = flux_fert
!          do p = col_pp%pfti(c), col_pp%pftf(c)
!             ! NFERTILIZATION gets the fertilizer applied + all manure N produced in the
!             ! column. Note that if fract_spread_grass > 0 then some of this N might be
!             ! still moved to the native veg. column.
!             nfertilization_patch(p) = manure_prod + sbgc_nf%fert_n_appl(c)
!          end do
!       end if
!    end do
!
!  end subroutine fan_to_sminn

end module FanUpdateMod
