module elm_driver

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides the main ELM driver physics calling sequence.  Most
  ! computations occurs over ``clumps'' of gridcells (and associated subgrid
  ! scale entities) assigned to each MPI process. Computation is further
  ! parallelized by looping over clumps on each process using shared memory OpenMP.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_sys_mod            , only : shr_sys_flush
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use elm_varpar             , only : nlevtrc_soil, nlevsoi
  use elm_varctl             , only : wrtdia, iulog, create_glacier_mec_landunit, use_fates, use_betr, use_firn_percolation_and_compaction
  use elm_varctl             , only : use_cn, use_lch4, use_voc, use_noio, use_c13, use_c14
  use elm_varctl             , only : use_erosion, use_fates_sp, use_fan
  use elm_varctl             , only : mpi_sync_nstep_freq
  use elm_time_manager       , only : get_step_size, get_curr_date, get_ref_date, get_nstep, is_beg_curr_day, get_curr_time_string
  use elm_time_manager       , only : get_curr_calday, get_days_per_year
  use elm_varpar             , only : nlevsno, nlevgrnd, crop_prog
  use spmdMod                , only : masterproc, mpicom
  use decompMod              , only : get_proc_clumps, get_clump_bounds, get_proc_bounds, bounds_type
  use filterMod              , only : filter, filter_inactive_and_active
  use histFileMod            , only : hist_update_hbuf, hist_htapes_wrapup
  use restFileMod            , only : restFile_write, restFile_filename
  use abortutils             , only : endrun
  !
  use dynSubgridDriverMod    , only : dynSubgrid_driver
  use BalanceCheckMod        , only : BeginColWaterBalance, ColWaterBalanceCheck
  use BalanceCheckMod        , only : BeginGridWaterBalance, GridBalanceCheck
  !
  use CanopyTemperatureMod   , only : CanopyTemperature ! (formerly Biogeophysics1Mod)
  use SoilTemperatureMod     , only : SoilTemperature
  use LakeTemperatureMod     , only : LakeTemperature
  !
  use BareGroundFluxesMod    , only : BareGroundFluxes
  use CanopyFluxesMod        , only : CanopyFluxes
  use SedYieldMod            , only : SoilErosion
  use SoilFluxesMod          , only : SoilFluxes ! (formerly Biogeophysics2Mod)
  use UrbanFluxesMod         , only : UrbanFluxes
  use LakeFluxesMod          , only : LakeFluxes
  !
  use HydrologyNoDrainageMod , only : HydrologyNoDrainage ! (formerly Hydrology2Mod)
  use HydrologyDrainageMod   , only : HydrologyDrainage   ! (formerly Hydrology2Mod)
  use CanopyHydrologyMod     , only : CanopyHydrology     ! (formerly Hydrology1Mod)
  use LakeHydrologyMod       , only : LakeHydrology
  !
  use AerosolMod             , only : AerosolMasses
  use SnowSnicarMod          , only : SnowAge_grain
  use SurfaceAlbedoMod       , only : SurfaceAlbedo
  use UrbanAlbedoMod         , only : UrbanAlbedo
  !
  use SurfaceRadiationMod    , only : SurfaceRadiation, CanopySunShadeFractions
  use UrbanRadiationMod      , only : UrbanRadiation
  !
  use SedFluxType            , only : sedflux_type
  !clm_interface
  use EcosystemDynMod      , only : EcosystemDynNoLeaching1, EcosystemDynNoLeaching2

  use EcosystemDynMod      , only : EcosystemDynLeaching
  use VegStructUpdateMod   , only : VegStructUpdate
  use AnnualUpdateMod      , only : AnnualUpdate
  use EcosystemBalanceCheckMod      , only : BeginColCBalance, BeginColNBalance, ColCBalanceCheck, ColNBalanceCheck
  use EcosystemBalanceCheckMod      , only : BeginColPBalance, ColPBalanceCheck
  use EcosystemBalanceCheckMod      , only : BeginGridCBalance, GridCBalanceCheck
  use EcosystemBalanceCheckMod      , only : BeginGridNBalance
  use EcosystemBalanceCheckMod      , only : BeginGridPBalance
  use EcosystemBalanceCheckMod      , only : EndGridCBalanceAfterDynSubgridDriver
  use EcosystemBalanceCheckMod      , only : EndGridNBalanceAfterDynSubgridDriver
  use EcosystemBalanceCheckMod      , only : EndGridPBalanceAfterDynSubgridDriver
  use VerticalProfileMod   , only : decomp_vertprofiles
  use FireMod              , only : FireInterp
  use SatellitePhenologyMod  , only : SatellitePhenology, interpMonthlyVeg
  use ndepStreamMod          , only : ndep_interp
  use pdepStreamMod          , only : pdep_interp
  use ActiveLayerMod         , only : alt_calc
  use CH4Mod                 , only : CH4
  use DUSTMod                , only : DustDryDep, DustEmission
  use VOCEmissionMod         , only : VOCEmission

  !
  use filterMod              , only : setFilters
  !
  use atm2lndMod             , only : downscale_forcings
  use lnd2atmMod             , only : lnd2atm
  use lnd2glcMod             , only : lnd2glc_type
  !
  use seq_drydep_mod         , only : n_drydep, drydep_method, DD_XLND
  use DryDepVelocity         , only : depvel_compute
  !
  use DaylengthMod           , only : UpdateDaylength
  use perf_mod
  !
  use elm_instMod            , only : ch4_vars, ep_betr
  use elm_instMod            , only : carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars
  use elm_instMod            , only : carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars
  use elm_instMod            , only : nitrogenstate_vars
  use elm_instMod            , only : nitrogenflux_vars
  use elm_instMod            , only : phosphorusstate_vars
  use elm_instMod            , only : phosphorusflux_vars
  use elm_instMod            , only : crop_vars
  use elm_instMod            , only : cnstate_vars
  use elm_instMod            , only : dust_vars
  use elm_instMod            , only : vocemis_vars
  use elm_instMod            , only : drydepvel_vars
  use elm_instMod            , only : aerosol_vars
  use elm_instMod            , only : canopystate_vars
  use elm_instMod            , only : energyflux_vars
  use elm_instMod            , only : frictionvel_vars
  use elm_instMod            , only : lakestate_vars
  use elm_instMod            , only : photosyns_vars
  use elm_instMod            , only : sedflux_vars
  use elm_instMod            , only : soilstate_vars
  use elm_instMod            , only : soilhydrology_vars
  use elm_instMod            , only : solarabs_vars
  use elm_instMod            , only : soilhydrology_vars
  use elm_instMod            , only : surfalb_vars
  use elm_instMod            , only : surfrad_vars
  use elm_instMod            , only : temperature_vars
  use elm_instMod            , only : col_es
  use elm_instMod            , only : waterflux_vars
  use elm_instMod            , only : waterstate_vars
  use elm_instMod            , only : atm2lnd_vars
  use elm_instMod            , only : lnd2atm_vars
  use elm_instMod            , only : glc2lnd_vars
  use elm_instMod            , only : lnd2glc_vars
  use elm_instMod            , only : soil_water_retention_curve
  use elm_instMod            , only : chemstate_vars
  use elm_instMod            , only : alm_fates
  use elm_instMod            , only : PlantMicKinetics_vars
  use elm_instMod            , only : sedflux_vars
  use tracer_varcon          , only : is_active_betr_bgc
  use CNEcosystemDynBetrMod  , only : CNEcosystemDynBetr, CNFluxStateBetrSummary
  use UrbanParamsType        , only : urbanparams_vars

  use GridcellType           , only : grc_pp
  use GridcellDataType       , only : grc_cs, c13_grc_cs, c14_grc_cs
  use GridcellDataType       , only : grc_cf, c13_grc_cf, c14_grc_cf
  use GridcellDataType       , only : grc_ns, grc_nf
  use GridcellDataType       , only : grc_ps, grc_pf
  use TopounitDataType       , only : top_as, top_af
  use LandunitType           , only : lun_pp
  use ColumnType             , only : col_pp
  use ColumnDataType         , only : col_es, col_ef, col_ws, col_wf
  use ColumnDataType         , only : col_cs, c13_col_cs, c14_col_cs
  use ColumnDataType         , only : col_cf, c13_col_cf, c14_col_cf
  use ColumnDataType         , only : col_ns, col_nf
  use ColumnDataType         , only : col_ps, col_pf
  use VegetationType         , only : veg_pp
  use VegetationDataType     , only : veg_es, veg_ws, veg_wf, veg_cf
  use VegetationDataType     , only : veg_cs, c13_veg_cs, c14_veg_cs
  use VegetationDataType     , only : veg_ns, veg_nf
  use VegetationDataType     , only : veg_ps, veg_pf
  use FanStreamMod           , only : fanstream_interp

  !----------------------------------------------------------------------------
  ! bgc interface & pflotran:
  use elm_varctl             , only : use_elm_interface
  use elm_instMod            , only : elm_interface_data
  use elm_interface_funcsMod , only : get_elm_data
  ! (1) clm_bgc through interface
  use elm_varctl             , only : use_elm_bgc
  use elm_interface_funcsMod , only : elm_bgc_run, update_bgc_data_elm2elm
  ! (2) pflotran
  use elm_time_manager            , only : nsstep, nestep
  use elm_varctl                  , only : use_pflotran, pf_cmode, pf_hmode, pf_tmode
  use elm_interface_funcsMod      , only : update_bgc_data_pf2elm, update_th_data_pf2elm
  use elm_interface_pflotranMod   , only : elm_pf_run, elm_pf_write_restart
  use elm_interface_pflotranMod   , only : elm_pf_finalize
  !----------------------------------------------------------------------------
  use WaterBudgetMod              , only : WaterBudget_Reset, WaterBudget_Run, WaterBudget_Accum, WaterBudget_Print
  use WaterBudgetMod              , only : WaterBudget_SetBeginningMonthlyStates
  use WaterBudgetMod              , only : WaterBudget_SetEndingMonthlyStates
  use CNPBudgetMod                , only : CNPBudget_Run, CNPBudget_Accum, CNPBudget_Print, CNPBudget_Reset
  use CNPBudgetMod                , only : CNPBudget_SetBeginningMonthlyStates, CNPBudget_SetEndingMonthlyStates
  use elm_varctl                  , only : do_budgets, budget_inst, budget_daily, budget_month
  use elm_varctl                  , only : budget_ann, budget_ltann, budget_ltend

  use timeinfoMod
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: elm_drv            ! Main elm driver
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: elm_drv_patch2col
  private :: elm_drv_init      ! Initialization of variables needed from previous timestep
  private :: write_diagnostic  ! Write diagnostic information to log file
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine elm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate)
    !
    ! !DESCRIPTION:
    !
    ! First phase of the clm driver calling the clm physics. An outline of
    ! the calling tree is given in the description of this module.
    !
    ! !USES:
     use elm_varctl            , only : fates_spitfire_mode
     use elm_varctl            , only : fates_seeddisp_cadence
     use FATESFireFactoryMod   , only : scalar_lightning
     use FatesInterfaceTypesMod, only : fates_dispersal_cadence_none
     
    ! !ARGUMENTS:
    implicit none
    logical ,        intent(in) :: doalb       ! true if time for surface albedo calc
    real(r8),        intent(in) :: nextsw_cday ! calendar day for nstep+1
    real(r8),        intent(in) :: declinp1    ! declination angle for next time step
    real(r8),        intent(in) :: declin      ! declination angle for current time step
    logical,         intent(in) :: rstwr       ! true => write restart file this step
    logical,         intent(in) :: nlend       ! true => end of run on this step
    character(len=*),intent(in) :: rdate       ! restart file time stamp for name
    !
    ! !LOCAL VARIABLES:
    integer              :: nstep                   ! time step number
    real(r8)             :: dtime                   ! land model time step (sec)
    integer              :: nc, c, p, l, g          ! indices
    integer              :: nclumps                 ! number of clumps on this processor
    integer              :: yrp1                    ! year (0, ...) for nstep+1
    integer              :: monp1                   ! month (1, ..., 12) for nstep+1
    integer              :: dayp1                   ! day of month (1, ..., 31) for nstep+1
    integer              :: secp1                   ! seconds into current date for nstep+1
    integer              :: yr                      ! year (0, ...)
    integer              :: mon                     ! month (1, ..., 12)
    integer              :: day                     ! day of month (1, ..., 31)
    integer              :: sec                     ! seconds of the day
    integer              :: ncdate                  ! current date
    integer              :: nbdate                  ! base date (reference date)
    integer              :: kyr                     ! thousand years, equals 2 at end of first year
    character(len=256)   :: filer                   ! restart file name
    integer              :: ier                     ! error code
    character(len=256)   :: dateTimeString
    type(bounds_type)    :: bounds_clump
    type(bounds_type)    :: bounds_proc
    !-----------------------------------------------------------------------

    call get_curr_time_string(dateTimeString)
    if (masterproc) then
       write(iulog,*)'Beginning timestep   : ',trim(dateTimeString)
       call shr_sys_flush(iulog)
    endif
    ! Determine processor bounds and clumps for this processor

    call get_proc_bounds(bounds_proc)
    nclumps = get_proc_clumps()
    nstep_mod = get_nstep()
    dtime_mod = real(get_step_size(),r8)
    call get_curr_date(year_curr,mon_curr, day_curr,secs_curr)
    dayspyr_mod = get_days_per_year()
    jday_mod = get_curr_calday()

    if (mpi_sync_nstep_freq > 0) then
       if (mod(nstep_mod,mpi_sync_nstep_freq) == 0) then
          call MPI_Barrier(mpicom, ier)
          if (masterproc) then
             write(iulog,*)'                       A MPI_Barrier is added in this timestep.'
          end if
       end if
    end if

    if (do_budgets) then
       call WaterBudget_Reset()

       if (use_cn) then
          call CNPBudget_Reset()
       end if
    end if


    ! ============================================================================
    ! Specified phenology
    ! ============================================================================

    if (use_cn) then
       ! For dry-deposition need to call CLMSP so that mlaidiff is obtained
       if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
          call t_startf('interpMonthlyVeg')
          call interpMonthlyVeg(bounds_proc, canopystate_vars)
          call t_stopf('interpMonthlyVeg')
       endif

    elseif(use_fates) then
       if(use_fates_sp) then
       
          ! For FATES satellite phenology mode interpolate the weights for
          ! time-interpolation of monthly vegetation data (as in SP mode below)
          ! Also for FATES with dry-deposition as above need to call CLMSP so that mlaidiff is obtained
          !if ( use_fates_sp .or. (n_drydep > 0 .and. drydep_method == DD_XLND ) ) then
          ! Replace with this when we have dry-deposition working
          ! For now don't allow for dry-deposition because of issues in #1044 EBK Jun/17/2022
          call t_startf('interpMonthlyVeg')
          call interpMonthlyVeg(bounds_proc, canopystate_vars)
          call t_stopf('interpMonthlyVeg')
       end if
    else
       ! Determine weights for time interpolation of monthly vegetation data.
       ! This also determines whether it is time to read new monthly vegetation and
       ! obtain updated leaf area index [mlai1,mlai2], stem area index [msai1,msai2],
       ! vegetation top [mhvt1,mhvt2] and vegetation bottom [mhvb1,mhvb2]. The
       ! weights obtained here are used in subroutine SatellitePhenology to obtain time
       ! interpolated values.
       if (doalb .or. ( n_drydep > 0 .and. drydep_method == DD_XLND )) then
          call t_startf('interpMonthlyVeg')
          call interpMonthlyVeg(bounds_proc, canopystate_vars)
          call t_stopf('interpMonthlyVeg')
       end if
    end if

    ! ==================================================================================
    ! Determine decomp vertical profiles
    !
    ! These routines (alt_calc & decomp_vertprofiles) need to be called before
    ! pftdyn_cnbal, and it appears that they need to be called before pftdyn_interp and
    ! the associated filter updates, too (otherwise we get a carbon balance error)
    ! ==================================================================================

    !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       call t_startf("decomp_vert")
       call alt_calc(filter(nc)%num_soilc, filter(nc)%soilc, &
            temperature_vars, canopystate_vars)

       !  Note (WJS, 6-12-13): Because of this routine's placement in the driver sequence
       !  (it is called very early in each timestep, before weights are adjusted and
       !  filters are updated), it may be necessary for this routine to compute values over
       !  inactive as well as active points (since some inactive points may soon become
       !  active) - so that's what is done now. Currently, it seems to be okay to do this,
       !  because the variables computed here seem to only depend on quantities that are
       !  valid over inactive as well as active points.

       if(use_fates .or. use_cn) then

          call decomp_vertprofiles(bounds_clump, &
               filter_inactive_and_active(nc)%num_soilc, &
               filter_inactive_and_active(nc)%soilc, &
               filter_inactive_and_active(nc)%num_soilp, &
               filter_inactive_and_active(nc)%soilp, &
               soilstate_vars, canopystate_vars, cnstate_vars)
       end if
       call t_stopf("decomp_vert")
    end do
    !$OMP END PARALLEL DO

    ! ============================================================================
    ! Zero fluxes for transient land cover
    ! ============================================================================

    !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       call t_startf('beggridwbal')
       call BeginGridWaterBalance(bounds_clump,               &
            filter(nc)%num_nolakec, filter(nc)%nolakec,       &
            filter(nc)%num_lakec, filter(nc)%lakec,           &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
            soilhydrology_vars )
       call t_stopf('beggridwbal')

       if (use_betr) then
         dtime=get_step_size(); nstep=get_nstep()
         call ep_betr%SetClock(dtime= dtime, nelapstep=nstep)
         call ep_betr%BeginMassBalanceCheck(bounds_clump)
       endif

       call t_startf('cnpinit')

       if (use_cn) then
          call t_startf('cnpvegzero')

          call veg_cs%ZeroDwt(bounds_clump)
          if (use_c13) then
             call c13_grc_cf%ZeroDWT(bounds_clump)
             call c13_col_cf%ZeroDWT(bounds_clump)
          end if
          if (use_c14) then
             call c14_grc_cf%ZeroDWT(bounds_clump)
             call c14_col_cf%ZeroDWT(bounds_clump)
          end if
          call veg_ns%ZeroDWT(bounds_clump)
          call veg_ps%ZeroDWT(bounds_clump)
          call t_stopf('cnpvegzero')
       end if

       if (use_cn .or. use_fates) then
          call t_startf('cnpzero')
          call grc_cf%ZeroDWT(bounds_clump)
          call col_cf%ZeroDWT(bounds_clump)
          call grc_nf%ZeroDWT(bounds_clump)
          call col_nf%ZeroDWT(bounds_clump)
          call grc_pf%ZeroDWT(bounds_clump)
          call col_pf%ZeroDWT(bounds_clump)
          call t_stopf('cnpzero')
       end if

       call t_startf('cnpvegsumm')
       if(use_cn) then
          call veg_cs%Summary(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, col_cs)
          call veg_ns%Summary(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, col_ns)
          call veg_ps%Summary(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, col_ps)

       elseif(use_fates)then
          ! In this scenario, we simply zero all of the
          ! column level variables that would had been upscaled
          ! in the veg summary with p2c
          call col_cs%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
          call col_ns%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
          call col_ps%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
       end if
       call t_stopf('cnpvegsumm')

       if(use_cn .or. use_fates)then

          call col_cs%Summary(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc)

          call col_ns%Summary(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc)

          call col_ps%Summary(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc)

          call BeginGridCBalance(bounds_clump, col_cs, grc_cs)
          call BeginGridNBalance(bounds_clump, col_ns, grc_ns)
          call BeginGridPBalance(bounds_clump, col_ps, grc_ps)

       end if

       call t_stopf('cnpinit')


    end do
    !$OMP END PARALLEL DO

    ! ============================================================================
    ! Update subgrid weights with dynamic landcover (prescribed transient patches,
    ! and or dynamic landunits), and do related adjustments. Note that this
    ! call needs to happen outside loops over nclumps.
    ! ============================================================================

    call t_startf('dyn_subgrid')
    call dynSubgrid_driver(bounds_proc,                                      &
       urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
       energyflux_vars, canopystate_vars, photosyns_vars, cnstate_vars,                       &
       veg_cs, c13_veg_cs, c14_veg_cs,         &
       col_cs, c13_col_cs, c14_col_cs, col_cf,  &
       grc_cs, grc_cf , glc2lnd_vars,  crop_vars)
    call t_stopf('dyn_subgrid')

    if (use_cn  .or. use_fates) then
       nstep = get_nstep()

       if (nstep < 2 )then
          if (masterproc) then
             write(iulog,*) '--WARNING-- skipping CN balance check for first timestep'
          end if
       else
          call t_startf('cnbalchk_at_grid')

          !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
          do nc = 1,nclumps
             call get_clump_bounds(nc, bounds_clump)

             if(use_cn) then
                call veg_cs%Summary(bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     filter(nc)%num_soilp, filter(nc)%soilp, col_cs)

                call veg_ns%Summary(bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     filter(nc)%num_soilp, filter(nc)%soilp, col_ns)

                call veg_ps%Summary(bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     filter(nc)%num_soilp, filter(nc)%soilp, col_ps)

             elseif(use_fates)then
                ! In this scenario, we simply zero all of the
                ! column level variables that would had been upscaled
                ! in the veg summary with p2c
                call col_cs%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
                call col_ns%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
                call col_ps%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
             end if

             call col_cs%Summary(bounds_clump, &
                  filter(nc)%num_soilc, filter(nc)%soilc)

             call col_ns%Summary(bounds_clump, &
                  filter(nc)%num_soilc, filter(nc)%soilc)

             call col_ps%Summary(bounds_clump, &
                  filter(nc)%num_soilc, filter(nc)%soilc)

             call EndGridCBalanceAfterDynSubgridDriver(bounds_clump, &
                  filter(nc)%num_soilc, filter(nc)%soilc, &
                  col_cs, grc_cs, grc_cf)

             call EndGridNBalanceAfterDynSubgridDriver(bounds_clump, &
                  filter(nc)%num_soilc, filter(nc)%soilc, &
                  col_ns, grc_ns, grc_nf)

             call EndGridPBalanceAfterDynSubgridDriver(bounds_clump, &
                  filter(nc)%num_soilc, filter(nc)%soilc, &
                  col_ps, grc_ps, grc_pf)

          end do
          !$OMP END PARALLEL DO
          call t_stopf('cnbalchk_at_grid')

       end if
    end if


    ! ============================================================================
    ! Initialize the mass balance checks for water.
    !
    ! Currently, I believe this needs to be done after weights are updated for
    ! prescribed transient patches, because column-level water is not
    ! generally conserved when weights change (instead the difference is put in
    ! the grid cell-level terms, qflx_liq_dynbal, etc.). In the future, we may
    ! want to change the balance checks to ensure that the grid cell-level water
    ! is conserved, considering qflx_liq_dynbal; in this case, the call to
    ! BeginWaterBalance should be moved to before the weight updates.
    !
    ! For CNP: This needs to be done after dynSubgrid_driver, because the
    ! changes due to dynamic area adjustments can break column-level conservation
    ! ============================================================================

    !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       call t_startf('begwbal')
       call BeginColWaterBalance(bounds_clump,                &
            filter(nc)%num_nolakec, filter(nc)%nolakec,       &
            filter(nc)%num_lakec, filter(nc)%lakec,           &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
            soilhydrology_vars )
       call t_stopf('begwbal')


       call t_startf('begcnpbal')
       ! call veg summary before col summary, for p2c
       if (use_cn) then
          call veg_cs%Summary(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, col_cs)
          call veg_ns%Summary(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, col_ns)
          call veg_ps%Summary(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, col_ps)
       elseif(use_fates)then
          ! In this scenario, we simply zero all of the
          ! column level variables that would had been upscaled
          ! in the veg summary with p2c
          call col_cs%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
          call col_ns%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
          call col_ps%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
       end if
       call t_stopf('begcnpbal')


       if (use_cn  .or. use_fates) then
          call t_startf('begcnpbalwf')
          call col_cs%Summary(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc)
          call col_ns%Summary(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc)
          call col_ps%Summary(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc)
          call BeginColCBalance(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               col_cs)
          call BeginColNBalance(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               col_ns)
          call BeginColPBalance(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               col_ps)

          call t_stopf('begcnpbalwf')
       end if


       if (do_budgets) then
          call WaterBudget_SetBeginningMonthlyStates(bounds_clump )
          if (use_cn) then
             call CNPBudget_SetBeginningMonthlyStates(bounds_clump, col_cs, grc_cs)
          endif
       endif

    end do
    !$OMP END PARALLEL DO



#ifndef CPL_BYPASS

    if (use_cn .or. use_fates) then

       ! ============================================================================
       ! Update dynamic N deposition field, on albedo timestep
       ! currently being done outside clumps loop, but no reason why it couldn't be
       ! re-written to go inside.
       ! ============================================================================

       call t_startf('ndep_interp')
       ! PET: switching CN timestep
       call ndep_interp(bounds_proc, atm2lnd_vars)
       call t_stopf('ndep_interp')
    end if


    if (use_cn) then
       call t_startf('fireinterp')
       call FireInterp(bounds_proc)
       call t_stopf('fireinterp')
    elseif (use_fates) then
       ! fates_spitfire_mode is assigned an integer value in the namelist
       ! see bld/namelist_files/namelist_definition.xml for details
       if (fates_spitfire_mode > scalar_lightning) then
          call alm_fates%InterpFileInputs(bounds_proc)
       end if
    end if

    if (use_cn .or. use_fates) then
       ! ============================================================================
       ! Update dynamic P deposition field, on albedo timestep
       ! currently being done outside clumps loop, but no reason why it couldn't be
       ! re-written to go inside.
       ! ============================================================================

       call t_startf('pdep_interp')
       ! PET: switching CN timestep
       call pdep_interp(bounds_proc, atm2lnd_vars)
       call t_stopf('pdep_interp')
    end if

#endif

    if (use_fan) then
       call fanstream_interp(bounds_proc, atm2lnd_vars)
    end if

    ! ============================================================================
    ! Initialize variables from previous time step, downscale atm forcings, and
    ! Determine canopy interception and precipitation onto ground surface.
    ! Determine the fraction of foliage covered by water and the fraction
    ! of foliage that is dry and transpiring. Initialize snow layer if the
    ! snow accumulation exceeds 10 mm.
    ! ============================================================================

    !$OMP PARALLEL DO PRIVATE (nc,l,c, bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       call t_startf('drvinit')

       call UpdateDaylength(bounds_clump, declin)

       ! Initialze variables needed for new driver time step
       call elm_drv_init(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_nolakep, filter(nc)%nolakep, &
            filter(nc)%num_soilp  , filter(nc)%soilp,   &
            canopystate_vars, energyflux_vars)

       call downscale_forcings(bounds_clump, &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c, &
            atm2lnd_vars)

       call t_stopf('drvinit')

       ! ============================================================================
       ! Canopy Hydrology
       ! (1) water storage of intercepted precipitation
       ! (2) direct throughfall and canopy drainage of precipitation
       ! (3) fraction of foliage covered by water and the fraction is dry and transpiring
       ! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
       ! ============================================================================

       call t_startf('canhydro')
       call CanopyHydrology(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_nolakep, filter(nc)%nolakep, &
            atm2lnd_vars, canopystate_vars, &
            aerosol_vars )
       call t_stopf('canhydro')

       ! ============================================================================
       ! Surface Radiation
       ! ============================================================================

       call t_startf('surfrad')

       ! Surface Radiation primarily for non-urban columns

       ! Most of the surface radiation calculations are agnostic to the forest-model
       ! but the calculations of the fractions of sunlit and shaded canopies
       ! are specific, calculate them first.
       ! The nourbanp filter is set in dySubgrid_driver (earlier in this call)
       ! over the patch index range defined by bounds_clump%begp:bounds_proc%endp

       if(use_fates) then
          call alm_fates%wrap_sunfrac(bounds_clump, top_af, canopystate_vars)
       else
          call CanopySunShadeFractions(filter(nc)%num_nourbanp, filter(nc)%nourbanp,    &
                                       atm2lnd_vars, surfalb_vars, canopystate_vars,    &
                                       solarabs_vars)
       end if

       call SurfaceRadiation(bounds_clump,                                 &
            filter(nc)%num_nourbanp, filter(nc)%nourbanp,                  &
            filter(nc)%num_urbanp, filter(nc)%urbanp    ,                  &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                      &
            atm2lnd_vars, canopystate_vars, surfalb_vars, &
            solarabs_vars, surfrad_vars)

       ! Surface Radiation for only urban columns

       call UrbanRadiation(bounds_clump,                                       &
            filter(nc)%num_nourbanl, filter(nc)%nourbanl,                      &
            filter(nc)%num_urbanl, filter(nc)%urbanl,                          &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                          &
            filter(nc)%num_urbanp, filter(nc)%urbanp,                          &
            atm2lnd_vars, urbanparams_vars, &
            solarabs_vars, surfalb_vars, energyflux_vars)

       call t_stopf('surfrad')

       ! ============================================================================
       ! Determine leaf temperature and surface fluxes based on ground
       ! temperature from previous time step.
       ! ============================================================================

       call t_startf('bgp1')
       call CanopyTemperature(bounds_clump,                                   &
            filter(nc)%num_nolakec, filter(nc)%nolakec,                       &
            filter(nc)%num_nolakep, filter(nc)%nolakep,                       &
            atm2lnd_vars, canopystate_vars, soilstate_vars, frictionvel_vars, &
            energyflux_vars)
       call t_stopf('bgp1')

       ! ============================================================================
       ! Determine fluxes
       ! ============================================================================

       call t_startf('bgflux')

       call col_wf%Reset(bounds_clump, filter(nc)%num_nolakec , filter(nc)%nolakec)

       ! Bareground fluxes for all patches except lakes and urban landunits

       call BareGroundFluxes(bounds_clump,                                 &
            filter(nc)%num_nolakeurbanp, filter(nc)%nolakeurbanp,          &
            atm2lnd_vars, canopystate_vars, soilstate_vars,                &
            frictionvel_vars, ch4_vars  )
       call t_stopf('bgflux')

       ! non-bareground fluxes for all patches except lakes and urban landunits
       ! Calculate canopy temperature, latent and sensible fluxes from the canopy,
       ! and leaf water change by evapotranspiration

       call t_startf('canflux')
       call CanopyFluxes(bounds_clump,                                                   &
            filter(nc)%num_nolakeurbanp, filter(nc)%nolakeurbanp,                        &
            atm2lnd_vars, canopystate_vars, cnstate_vars, energyflux_vars,               &
            frictionvel_vars, soilstate_vars, solarabs_vars, surfalb_vars,               &
            ch4_vars, photosyns_vars )
       call t_stopf('canflux')

       ! Fluxes for all urban landunits

       call t_startf('uflux')
       call UrbanFluxes(bounds_clump,                        &
            filter(nc)%num_nourbanl, filter(nc)%nourbanl,    &
            filter(nc)%num_urbanl, filter(nc)%urbanl,        &
            filter(nc)%num_urbanc, filter(nc)%urbanc,        &
            filter(nc)%num_urbanp, filter(nc)%urbanp,        &
            atm2lnd_vars, urbanparams_vars, soilstate_vars,  &
            frictionvel_vars, energyflux_vars)
       call t_stopf('uflux')

       ! Fluxes for all lake landunits

       call t_startf('bgplake')
       call LakeFluxes(bounds_clump,                                         &
            filter(nc)%num_lakec, filter(nc)%lakec,                          &
            filter(nc)%num_lakep, filter(nc)%lakep,                          &
            atm2lnd_vars, solarabs_vars, frictionvel_vars,  &
            energyflux_vars, lakestate_vars)

       ! ============================================================================
       ! DUST and VOC emissions
       ! ============================================================================

       call t_startf('bgc')

       ! Dust mobilization (C. Zender's modified codes)
       call DustEmission(bounds_clump,                      &
            filter(nc)%num_nolakep, filter(nc)%nolakep,     &
            atm2lnd_vars, soilstate_vars, canopystate_vars, &
            frictionvel_vars, dust_vars)

       ! Dust dry deposition (C. Zender's modified codes)
       call DustDryDep(bounds_clump, &
            atm2lnd_vars, frictionvel_vars, dust_vars)

       ! VOC emission (A. Guenther's MEGAN (2006) model)
       if (use_voc) then
          call VOCEmission(bounds_clump,                                         &
               filter(nc)%num_soilp, filter(nc)%soilp,                           &
               atm2lnd_vars, canopystate_vars, photosyns_vars, temperature_vars, &
               vocemis_vars)
       end if

       call t_stopf('bgc')

       ! ============================================================================
       ! Determine temperatures
       ! ============================================================================
       if(use_betr)then
         call ep_betr%BeTRSetBiophysForcing(bounds_clump, col_pp, veg_pp, 1, nlevsoi, waterstate_vars=col_ws)
         call ep_betr%PreDiagSoilColWaterFlux(filter(nc)%num_nolakec , filter(nc)%nolakec)
       endif
       ! Set lake temperature

       call LakeTemperature(bounds_clump,             &
            filter(nc)%num_lakec, filter(nc)%lakec,   &
            filter(nc)%num_lakep, filter(nc)%lakep,   &
            solarabs_vars, soilstate_vars,  ch4_vars, &
            energyflux_vars, lakestate_vars)
       call t_stopf('bgplake')

       ! Set soil/snow temperatures including ground temperature

       call t_startf('soiltemperature')
       call SoilTemperature(bounds_clump,                     &
            filter(nc)%num_urbanl  , filter(nc)%urbanl,       &
            filter(nc)%num_nolakec , filter(nc)%nolakec,      &
            atm2lnd_vars, urbanparams_vars, canopystate_vars, &
            solarabs_vars, soilstate_vars, energyflux_vars )
       call t_stopf('soiltemperature')


       if(use_betr)then
         call ep_betr%BeTRSetBiophysForcing(bounds_clump, col_pp, veg_pp, 1, nlevsoi, waterstate_vars=col_ws)
         call ep_betr%DiagnoseDtracerFreezeThaw(bounds_clump, filter(nc)%num_nolakec , filter(nc)%nolakec, col_pp, lun_pp)
       endif
       ! ============================================================================
       ! update surface fluxes for new ground temperature.
       ! ============================================================================

       call t_startf('bgp2')
       call SoilFluxes(bounds_clump,                       &
            filter(nc)%num_urbanl,  filter(nc)%urbanl,     &
            filter(nc)%num_nolakec, filter(nc)%nolakec,    &
            filter(nc)%num_nolakep, filter(nc)%nolakep,    &
            atm2lnd_vars, solarabs_vars, canopystate_vars, &
            energyflux_vars )
       call t_stopf('bgp2')

       ! ============================================================================
       ! Perform averaging from patch level to column level
       ! ============================================================================

       call t_startf('patch2col')
       call elm_drv_patch2col(bounds_clump, filter(nc)%num_nolakec, filter(nc)%nolakec, energyflux_vars)
       call t_stopf('patch2col')

       ! ============================================================================
       ! Vertical (column) soil and surface hydrology
       ! ============================================================================

       ! Note that filter_snowc and filter_nosnowc are returned by
       ! LakeHydrology after the new snow filter is built

       call t_startf('hydro without drainage')

       call HydrologyNoDrainage(bounds_clump,                                &
            filter(nc)%num_nolakec, filter(nc)%nolakec,                      &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc,                &
            filter(nc)%num_hydrononsoic, filter(nc)%hydrononsoic,            &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                        &
            filter(nc)%num_snowc, filter(nc)%snowc,                          &
            filter(nc)%num_nosnowc, filter(nc)%nosnowc,canopystate_vars,     &
            atm2lnd_vars, lnd2atm_vars, soilstate_vars, energyflux_vars,     &
            soilhydrology_vars, aerosol_vars )

       !  Calculate column-integrated aerosol masses, and
       !  mass concentrations for radiative calculations and output
       !  (based on new snow level state, after SnowFilter is rebuilt.
       !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there
       !  can be zero snow layers but an active column in filter)

       call AerosolMasses( bounds_clump,                                   &
            num_on=filter(nc)%num_snowc, filter_on=filter(nc)%snowc,       &
            num_off=filter(nc)%num_nosnowc, filter_off=filter(nc)%nosnowc, &
            aerosol_vars=aerosol_vars)

       call t_stopf('hydro without drainage')

       ! ============================================================================
       ! Lake hydrology
       ! ============================================================================

       ! Note that filter_lakesnowc and filter_lakenosnowc are returned by
       ! LakeHydrology after the new snow filter is built

       call t_startf('hylake')
       call LakeHydrology(bounds_clump,                                                      &
            filter(nc)%num_lakec, filter(nc)%lakec,                                          &
            filter(nc)%num_lakep, filter(nc)%lakep,                                          &
            filter(nc)%num_lakesnowc, filter(nc)%lakesnowc,                                  &
            filter(nc)%num_lakenosnowc, filter(nc)%lakenosnowc,                              &
            atm2lnd_vars, soilstate_vars,  &
            energyflux_vars, aerosol_vars, lakestate_vars)

       !  Calculate column-integrated aerosol masses, and
       !  mass concentrations for radiative calculations and output
       !  (based on new snow level state, after SnowFilter is rebuilt.
       !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there
       !  can be zero snow layers but an active column in filter)

       call AerosolMasses(bounds_clump,                                            &
            num_on=filter(nc)%num_lakesnowc, filter_on=filter(nc)%lakesnowc,       &
            num_off=filter(nc)%num_lakenosnowc, filter_off=filter(nc)%lakenosnowc, &
            aerosol_vars=aerosol_vars)

       ! Must be done here because must use a snow filter for lake columns

       call SnowAge_grain(bounds_clump,                         &
            filter(nc)%num_lakesnowc, filter(nc)%lakesnowc,     &
            filter(nc)%num_lakenosnowc, filter(nc)%lakenosnowc )


       call t_stopf('hylake')

       ! ============================================================================
       ! ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)
       ! ============================================================================

       do c = bounds_clump%begc,bounds_clump%endc
          l = col_pp%landunit(c)
          if (lun_pp%urbpoi(l)) then
             ! Urban landunit use Bonan 1996 (LSM Technical Note)
             col_ws%frac_sno(c) = min( col_ws%snow_depth(c)/0.05_r8, 1._r8)
          end if
       end do

       ! ============================================================================
       ! Snow aging routine based on Flanner and Zender (2006), Linking snowpack
       ! microphysics and albedo evolution, JGR, and Brun (1989), Investigation of
       ! wet-snow metamorphism in respect of liquid-water content, Ann. Glaciol.
       ! ============================================================================
       ! Note the snow filters here do not include lakes
       ! TODO: move this up

       call t_startf('snow_init')
       call SnowAge_grain(bounds_clump,                 &
            filter(nc)%num_snowc, filter(nc)%snowc,     &
            filter(nc)%num_nosnowc, filter(nc)%nosnowc )
       call t_stopf('snow_init')

       ! ============================================================================
       ! Update sediment fluxes from land unit
       ! ============================================================================

       if (use_cn .and. use_erosion) then
          call t_startf('erosion')
          call SoilErosion(bounds_clump, filter(nc)%num_soilc, filter(nc)%soilc, &
               canopystate_vars, cnstate_vars, soilstate_vars, sedflux_vars)
          call t_stopf('erosion')
       end if

       ! ============================================================================
       ! Ecosystem dynamics: Uses CN, or static parameterizations
       ! ============================================================================

       call t_startf('ecosysdyn')
       if (use_cn)then
          call crop_vars%CropIncrementYear(filter(nc)%num_pcropp, filter(nc)%pcropp)
       endif

       if(use_betr)then
         !right now betr bgc is intended only for non-ed mode

         if(is_active_betr_bgc)then
           !this returns the plant nutrient demand to soil bgc
           call CNEcosystemDynBetr(bounds_clump,                                &
                 filter(nc)%num_soilc, filter(nc)%soilc,                        &
                 filter(nc)%num_soilp, filter(nc)%soilp,                        &
                 filter(nc)%num_pcropp, filter(nc)%pcropp, doalb,               &
                 cnstate_vars, carbonflux_vars, carbonstate_vars,               &
                 c13_carbonflux_vars, c13_carbonstate_vars,                     &
                 c14_carbonflux_vars, c14_carbonstate_vars,                     &
                 nitrogenflux_vars, nitrogenstate_vars,                         &
                 atm2lnd_vars, waterstate_vars, waterflux_vars,                 &
                 canopystate_vars, soilstate_vars, temperature_vars, crop_vars, &
                 photosyns_vars, soilhydrology_vars, energyflux_vars,&
                 PlantMicKinetics_vars,                                         &
                 phosphorusflux_vars, phosphorusstate_vars, frictionvel_vars)

           call AnnualUpdate(bounds_clump,            &
                  filter(nc)%num_soilc, filter(nc)%soilc, &
                  filter(nc)%num_soilp, filter(nc)%soilp, &
                  cnstate_vars)
         endif
       endif

       ! FIX(SPM,032414)  push these checks into the routines below and/or make this consistent.

       if( .not. is_active_betr_bgc) then

          if (use_cn .or. use_fates) then

             ! fully prognostic canopy structure and C-N biogeochemistry
             ! - crop model:  crop algorithms called from within CNEcosystemDyn

             !===========================================================================================
             ! clm_interface: 'EcosystemDynNoLeaching' is divided into 2 subroutines (1 & 2): BEGIN
             ! EcosystemDynNoLeaching1 is called before clm_interface
             ! EcosystemDynNoLeaching2 is called after clm_interface
             !===========================================================================================
             call EcosystemDynNoLeaching1(bounds_clump,         &
                       filter(nc)%num_soilc, filter(nc)%soilc,  &
                       filter(nc)%num_soilp, filter(nc)%soilp,  &
                       filter(nc)%num_pcropp, filter(nc)%pcropp, &
                       cnstate_vars,            &
                       atm2lnd_vars,            &
                       canopystate_vars, soilstate_vars, crop_vars,   &
                       ch4_vars, photosyns_vars, frictionvel_vars )

             !--------------------------------------------------------------------------------
             if (use_elm_interface) then
                 ! STEP-1: pass data from CLM to elm_interface_data (INTERFACE DATA TYPE)
                 call get_elm_data(elm_interface_data,bounds_clump,                     &
                           filter(nc)%num_soilc, filter(nc)%soilc,                      &
                           filter(nc)%num_soilp, filter(nc)%soilp,                      &
                           atm2lnd_vars, soilstate_vars,                                &
                           waterstate_vars, waterflux_vars,                             &
                           temperature_vars, energyflux_vars,                           &
                           cnstate_vars, carbonflux_vars, carbonstate_vars,             &
                           nitrogenflux_vars, nitrogenstate_vars,                       &
                           phosphorusflux_vars, phosphorusstate_vars,                   &
                           ch4_vars)


                 if (use_pflotran .and. pf_cmode) then
                    call t_startf('pflotran')
                    ! -------------------------------------------------------------------------
                    ! PFLOTRAN calling for solving below-ground and ground-surface processes,
                    ! including thermal, hydrological and biogeochemical processes
                    ! STEP-2: (1) pass data from elm_interface_data to pflotran
                    ! STEP-2: (2) run pflotran
                    ! STEP-2: (3) update elm_interface_data from pflotran
                    ! -------------------------------------------------------------------------
                    call elm_pf_run(elm_interface_data, bounds_clump, filter, nc)

                    ! STEP-3: update CLM from elm_interface_data
                    call update_bgc_data_pf2elm(elm_interface_data%bgc,         &
                           bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc, &
                           filter(nc)%num_soilp, filter(nc)%soilp,              &
                           cnstate_vars, carbonflux_vars, carbonstate_vars,     &
                           nitrogenflux_vars, nitrogenstate_vars,               &
                           phosphorusflux_vars, phosphorusstate_vars,           &
                           ch4_vars)

                    call t_stopf('pflotran')

                 elseif (use_elm_bgc) then
                    call t_startf('elm-bgc via interface')
                    ! -------------------------------------------------------------------------
                    ! run elm-bgc (SoilLittDecompAlloc) through interface
                    ! STEP-2: (1) pass data from elm_interface_data to SoilLittDecompAlloc
                    ! STEP-2: (2) run SoilLittDecompAlloc
                    ! STEP-2: (3) update elm_interface_data from SoilLittDecompAlloc
                    ! -------------------------------------------------------------------------
                    call elm_bgc_run(elm_interface_data, bounds_clump,          &
                           filter(nc)%num_soilc, filter(nc)%soilc,              &
                           filter(nc)%num_soilp, filter(nc)%soilp,              &
                           canopystate_vars, soilstate_vars,                    &
                           temperature_vars, waterstate_vars,                   &
                           cnstate_vars, ch4_vars,                              &
                           carbonstate_vars, carbonflux_vars,                   &
                           nitrogenstate_vars, nitrogenflux_vars,               &
                           phosphorusstate_vars,phosphorusflux_vars)

                    ! STEP-3: update CLM from elm_interface_data
                    call update_bgc_data_elm2elm(elm_interface_data%bgc,        &
                           bounds_clump, filter(nc)%num_soilc, filter(nc)%soilc,&
                           filter(nc)%num_soilp, filter(nc)%soilp,              &
                           cnstate_vars, carbonflux_vars, carbonstate_vars,     &
                           nitrogenflux_vars, nitrogenstate_vars,               &
                           phosphorusflux_vars, phosphorusstate_vars,           &
                           ch4_vars)
                    call t_stopf('elm-bgc via interface')
                 end if !if (use_pflotran .and. pf_cmode)
             end if !if (use_elm_interface)
             !--------------------------------------------------------------------------------

             call EcosystemDynNoLeaching2(bounds_clump,                &
                   filter(nc)%num_soilc, filter(nc)%soilc,             &
                   filter(nc)%num_soilp, filter(nc)%soilp,             &
                   filter(nc)%num_pcropp, filter(nc)%pcropp, doalb,    &
                   filter(nc)%num_ppercropp, filter(nc)%ppercropp,     &
                   cnstate_vars,  atm2lnd_vars,          &
                   canopystate_vars, soilstate_vars, crop_vars, ch4_vars, &
                   photosyns_vars, soilhydrology_vars, energyflux_vars,   &
                   sedflux_vars, solarabs_vars)

             !===========================================================================================
             ! clm_interface: 'EcosystemDynNoLeaching' is divided into 2 subroutines (1 & 2): END
             !===========================================================================================
             if(.not.use_fates)then
               call AnnualUpdate(bounds_clump,            &
                    filter(nc)%num_soilc, filter(nc)%soilc, &
                    filter(nc)%num_soilp, filter(nc)%soilp, &
                    cnstate_vars)
             end if
             
             if (use_fates_sp) then
               call SatellitePhenology(bounds_clump,               &
               filter_inactive_and_active(nc)%num_soilp, filter_inactive_and_active(nc)%soilp,    &
               waterstate_vars, canopystate_vars)
             endif
             
          else ! not ( if-use_cn   or if-use_fates)
             if (doalb) then
                ! Prescribed biogeography - prescribed canopy structure, some prognostic carbon fluxes
                call SatellitePhenology(bounds_clump,               &
                     filter(nc)%num_nolakep, filter(nc)%nolakep,    &
                     waterstate_vars, canopystate_vars)
             end if
          end if  ! end of if-use_cn   or if-use_fates
       end if ! end of is_active_betr_bgc

       call t_stopf('ecosysdyn')

       ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)
       call t_startf('depvel')
       if(.not.use_fates)then
         call depvel_compute(bounds_clump, &
              atm2lnd_vars, canopystate_vars, frictionvel_vars, &
              photosyns_vars, drydepvel_vars)
       end if
       call t_stopf('depvel')

       if (use_betr)then
          call ep_betr%CalcSmpL(bounds_clump, 1, nlevsoi, filter(nc)%num_soilc, filter(nc)%soilc, &
               col_es%t_soisno(bounds_clump%begc:bounds_clump%endc,1:nlevsoi), &
               soilstate_vars, col_ws, soil_water_retention_curve)

          call ep_betr%SetBiophysForcing(bounds_clump, col_pp, veg_pp,                         &
             carbonflux_vars=col_cf,     pf_carbonflux_vars=veg_cf,                          &
             waterstate_vars=col_ws,         waterflux_vars=col_wf, pf_waterflux_vars=veg_wf,        &
             temperature_vars=col_es, pf_temperature_vars=veg_es,  soilhydrology_vars=soilhydrology_vars, &
             atm2lnd_vars=atm2lnd_vars,               canopystate_vars=canopystate_vars,     &
             chemstate_vars=chemstate_vars,           soilstate_vars=soilstate_vars, &
             cnstate_vars = cnstate_vars, carbonstate_vars=col_cs)

          if(is_active_betr_bgc)then
             call ep_betr%PlantSoilBGCSend(bounds_clump, col_pp, veg_pp, &
                  filter(nc)%num_soilc,  filter(nc)%soilc, cnstate_vars, &
               col_cs, col_cf, c13_col_cs, c13_col_cf, c14_col_cs, c14_col_cf, &
               col_ns, col_nf, col_ps, col_pf,&
               PlantMicKinetics_vars)                  
          endif
          call ep_betr%StepWithoutDrainage(bounds_clump, col_pp, veg_pp)
       endif  !end use_betr

       if (use_lch4 .and. .not. is_active_betr_bgc) then
          !warning: do not call ch4 before AnnualUpdate, which will fail the ch4 model
          call t_startf('ch4')
          call CH4 (bounds_clump,                                                                  &
               filter(nc)%num_soilc, filter(nc)%soilc,                                             &
               filter(nc)%num_lakec, filter(nc)%lakec,                                             &
               filter(nc)%num_soilp, filter(nc)%soilp,                                             &
               atm2lnd_vars, lakestate_vars, canopystate_vars, soilstate_vars, soilhydrology_vars, &
               energyflux_vars, ch4_vars, lnd2atm_vars, alm_fates)
          call t_stopf('ch4')
       end if

       ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)
       call t_startf('depvel')
       if(.not.use_fates)then
         call depvel_compute(bounds_clump, &
              atm2lnd_vars, canopystate_vars, frictionvel_vars, &
              photosyns_vars, drydepvel_vars)
       end if
       call t_stopf('depvel')
       ! ============================================================================
       ! Calculate soil/snow hydrology with drainage (subsurface runoff)
       ! ============================================================================

       call t_startf('hydro2 drainage')

       if (use_elm_interface .and. (use_pflotran .and. pf_hmode)) then
         ! pflotran only works on 'soilc' (already done above).
         ! here for non-soil hydrology columns
         call HydrologyDrainage(bounds_clump,                     &
            filter(nc)%num_nolakec, filter(nc)%nolakec,           &
            filter(nc)%num_hydrononsoic, filter(nc)%hydrononsoic, &
            filter(nc)%num_urbanc, filter(nc)%urbanc,             &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,         &
            atm2lnd_vars, glc2lnd_vars,        &
            soilhydrology_vars, soilstate_vars)

       else

         call HydrologyDrainage(bounds_clump,                 &
            filter(nc)%num_nolakec, filter(nc)%nolakec,       &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
            filter(nc)%num_urbanc, filter(nc)%urbanc,         &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,     &
            atm2lnd_vars, glc2lnd_vars,      &
            soilhydrology_vars, soilstate_vars)

       end if

       call t_stopf('hydro2 drainage')

       if (use_betr) then
          call t_startf('betr drainage')
          call ep_betr%StepWithDrainage(bounds_clump, col_pp)
          call t_stopf('betr drainage')

          call t_startf('betr balchk')
          call ep_betr%MassBalanceCheck(bounds_clump)
          call t_stopf('betr balchk')
          call ep_betr%HistRetrieval(filter(nc)%num_nolakec, filter(nc)%nolakec)

          if(is_active_betr_bgc)then

            !extract nitrogen pool and flux from betr
            call ep_betr%PlantSoilBGCRecv(bounds_clump, col_pp, veg_pp, filter(nc)%num_soilc, filter(nc)%soilc,&
               col_cs, col_cf, veg_cf, c13_col_cs, c13_col_cf, &
               c14_col_cs, c14_col_cf, &
               col_ns, veg_ns, col_nf, veg_nf, col_ps, col_pf, veg_pf)
            !summarize total column nitrogen and carbon
            call CNFluxStateBetrSummary(bounds_clump, col_pp, veg_pp, &
                 filter(nc)%num_soilc, filter(nc)%soilc,                       &
                 filter(nc)%num_soilp, filter(nc)%soilp,                       &
                 carbonflux_vars, carbonstate_vars,                            &
                 c13_carbonflux_vars, c13_carbonstate_vars,                    &
                 c14_carbonflux_vars, c14_carbonstate_vars,                    &
                 nitrogenflux_vars, nitrogenstate_vars,                        &
                 phosphorusflux_vars, phosphorusstate_vars)
          endif
       endif  !end use_betr


       if (use_cn .or. use_fates) then
          if (.not. is_active_betr_bgc)then
            call EcosystemDynLeaching(bounds_clump,                &
               filter(nc)%num_soilc, filter(nc)%soilc,             &
               filter(nc)%num_soilp, filter(nc)%soilp,             &
               filter(nc)%num_pcropp, filter(nc)%pcropp, doalb,    &
               cnstate_vars,  frictionvel_vars, canopystate_vars )
         end if
       end if

       ! ============================================================================
       ! Update Vegetation
       ! ============================================================================

       ! Execute FATES dynamics
       if ( use_fates ) then
          
          ! FATES has its own running mean functions, such as 24hr
          ! vegetation temperature and exponential moving averages
          ! for leaf photosynthetic acclimation temperature. These
          ! moving averages are updated here
          call alm_fates%WrapUpdateFatesRmean(nc)
          
           ! Update high-frequency history diagnostics for FATES
           call alm_fates%wrap_update_hifrq_hist(bounds_clump)
           if ( is_beg_curr_day() ) then ! run ED at the start of each day
               call alm_fates%dynamics_driv( bounds_clump, top_as,          &
                    top_af, atm2lnd_vars, soilstate_vars, &
                    canopystate_vars, frictionvel_vars, soil_water_retention_curve)
           end if
       end if

       if (use_cn .and. doalb) then
           call VegStructUpdate(filter(nc)%num_soilp, filter(nc)%soilp,   &
                frictionvel_vars, cnstate_vars, &
                canopystate_vars, crop_vars)
       end if


       ! ============================================================================
       ! Check the energy and water balance, also carbon and nitrogen balance
       ! ============================================================================

       call t_startf('balchk')
       call ColWaterBalanceCheck(bounds_clump, &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c, &
            atm2lnd_vars, glc2lnd_vars, solarabs_vars,  &
            energyflux_vars, canopystate_vars)
       call t_stopf('balchk')

       call t_startf('gridbalchk')
       call GridBalanceCheck(bounds_clump                  , &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c   , &
            atm2lnd_vars, glc2lnd_vars, solarabs_vars,       &
            energyflux_vars, canopystate_vars              , &
            soilhydrology_vars)
       call t_stopf('gridbalchk')

       if (do_budgets) then
          call WaterBudget_SetEndingMonthlyStates(bounds_clump)
          if (use_cn) then
             call CNPBudget_SetEndingMonthlyStates(bounds_clump, col_cs, grc_cs)
          endif
       endif

       if (use_cn .or. use_fates) then

          call t_startf('cnbalchk')

          call ColCBalanceCheck(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               col_cs, col_cf)

          call ColNBalanceCheck(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               col_ns, col_nf)

          call ColPBalanceCheck(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               col_ps, col_pf)

          call GridCBalanceCheck(bounds_clump, col_cs, col_cf, grc_cs, grc_cf)

          call t_stopf('cnbalchk')
       end if


       ! ============================================================================
       ! Determine albedos for next time step
       ! ============================================================================

       if (doalb) then

          ! Albedos for non-urban columns
          call t_startf('surfalb')
          call SurfaceAlbedo(bounds_clump,                      &
               filter_inactive_and_active(nc)%num_nourbanc,     &
               filter_inactive_and_active(nc)%nourbanc,         &
               filter_inactive_and_active(nc)%num_nourbanp,     &
               filter_inactive_and_active(nc)%nourbanp,         &
               filter_inactive_and_active(nc)%num_urbanc,       &
               filter_inactive_and_active(nc)%urbanc,           &
               filter_inactive_and_active(nc)%num_urbanp,       &
               filter_inactive_and_active(nc)%urbanp,           &
               nextsw_cday, declinp1,                           &
               aerosol_vars, canopystate_vars, &
               lakestate_vars, surfalb_vars )
          call t_stopf('surfalb')

          ! Albedos for urban columns
          if (filter_inactive_and_active(nc)%num_urbanl > 0) then
             call t_startf('urbsurfalb')
             call UrbanAlbedo(bounds_clump,                  &
                  filter_inactive_and_active(nc)%num_urbanl, &
                  filter_inactive_and_active(nc)%urbanl,     &
                  filter_inactive_and_active(nc)%num_urbanc, &
                  filter_inactive_and_active(nc)%urbanc,     &
                  filter_inactive_and_active(nc)%num_urbanp, &
                  filter_inactive_and_active(nc)%urbanp,     &
                  urbanparams_vars, solarabs_vars, surfalb_vars)
             call t_stopf('urbsurfalb')
          end if

       end if

    end do
    !$OMP END PARALLEL DO

    ! Pass fates seed dispersal information to all nodes
    if (use_fates) then
       if (fates_seeddisp_cadence /= fates_dispersal_cadence_none) call alm_fates%WrapGlobalSeedDispersal()
    end if
    
    ! ============================================================================
    ! Determine gridcell averaged properties to send to atm
    ! ============================================================================

    if(use_betr)then
      call ep_betr%DiagnoseLnd2atm(bounds_proc, col_pp, lnd2atm_vars)
    endif

    call t_startf('lnd2atm')
    call lnd2atm(bounds_proc,                                   &
         atm2lnd_vars, surfalb_vars, frictionvel_vars,          &
         energyflux_vars, solarabs_vars, drydepvel_vars,        &
         vocemis_vars, dust_vars, ch4_vars, soilhydrology_vars, &
         sedflux_vars, lnd2atm_vars)
    call t_stopf('lnd2atm')

    ! ============================================================================
    ! Determine gridcell averaged properties to send to glc
    ! ============================================================================

    if (create_glacier_mec_landunit) then
       call t_startf('lnd2glc')
       !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
       do nc = 1,nclumps
          call get_clump_bounds(nc, bounds_clump)
          call lnd2glc_vars%update_lnd2glc(bounds_clump,       &
               filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,   &
               init=.false.)
       end do
       !$OMP END PARALLEL DO
       call t_stopf('lnd2glc')
    end if

    ! ============================================================================
    ! Write global average diagnostics to standard output
    ! ============================================================================

    nstep = get_nstep()
    if (wrtdia) call mpi_barrier(mpicom,ier)
    call t_startf('wrtdiag')
    call write_diagnostic(bounds_proc, wrtdia, nstep, lnd2atm_vars)
    call t_stopf('wrtdiag')

    ! ============================================================================
    ! Update accumulators
    ! ============================================================================

    ! FIX(SPM,032414) double check why this isn't called for ED

    if (nstep > 0) then

       call t_startf('accum')

       call atm2lnd_vars%UpdateAccVars(bounds_proc)

       call top_as%UpdateAccVars(bounds_proc)

       call top_af%UpdateAccVars(bounds_proc)

       call veg_es%UpdateAccVars(bounds_proc)

       call canopystate_vars%UpdateAccVars(bounds_proc)

       if (crop_prog) then
          call crop_vars%UpdateAccVars(bounds_proc, temperature_vars)
       end if

       call cnstate_vars%UpdateAccVars(bounds_proc)
       
       if(use_fates) then
          call alm_fates%UpdateAccVars(bounds_proc)
       end if

       call t_stopf('accum')

    end if

    ! ============================================================================
    ! Update history buffer
    ! ============================================================================

    call t_startf('hbuf')
    call hist_update_hbuf(bounds_proc)
    call t_stopf('hbuf')

    ! ============================================================================
    ! Compute water budget
    ! ============================================================================
    if (get_nstep()>0 .and. do_budgets) then
       call WaterBudget_Run(bounds_proc, atm2lnd_vars, lnd2atm_vars, &
            soilhydrology_vars)
       call WaterBudget_Accum()
       call WaterBudget_Print(budget_inst,  budget_daily,  budget_month,  &
            budget_ann,  budget_ltann,  budget_ltend)

       if (use_cn .and. do_budgets) then
          call CNPBudget_Run(bounds_proc, atm2lnd_vars, lnd2atm_vars, grc_cs, grc_cf)
          call CNPBudget_Accum()
          call CNPBudget_Print(budget_inst,  budget_daily,  budget_month,  &
               budget_ann,  budget_ltann,  budget_ltend)
      end if
    endif

    ! ============================================================================
    ! History/Restart output
    ! ============================================================================

    if (.not. use_noio) then

       call t_startf('elm_drv_io')

       ! Create history and write history tapes if appropriate
       call t_startf('elm_drv_io_htapes')

       call hist_htapes_wrapup( rstwr, nlend, bounds_proc,                    &
            soilstate_vars%watsat_col(bounds_proc%begc:bounds_proc%endc, 1:), &
            soilstate_vars%sucsat_col(bounds_proc%begc:bounds_proc%endc, 1:), &
            soilstate_vars%bsw_col(bounds_proc%begc:bounds_proc%endc, 1:),    &
            soilstate_vars%hksat_col(bounds_proc%begc:bounds_proc%endc, 1:))

       call t_stopf('elm_drv_io_htapes')
       ! Write restart/initial files if appropriate
       if (rstwr) then
          call t_startf('elm_drv_io_wrest')
          filer = restFile_filename(rdate=rdate)

          call restFile_write( bounds_proc, filer,                          &
               atm2lnd_vars, aerosol_vars, canopystate_vars, cnstate_vars,  &
               ch4_vars, energyflux_vars, frictionvel_vars, lakestate_vars, &
               photosyns_vars, soilhydrology_vars,     &
               soilstate_vars, solarabs_vars, surfalb_vars,  &
               sedflux_vars, ep_betr, alm_fates, crop_vars, rdate=rdate )

         !----------------------------------------------
         ! pflotran (off now)
         ! if (use_pflotran) then
         !     call elm_pf_write_restart(rdate)
         ! end if
         !----------------------------------------------


          call t_stopf('elm_drv_io_wrest')
       end if
       call t_stopf('elm_drv_io')

    end if

    if (use_pflotran .and. nstep>=nestep) then
       call elm_pf_finalize()
    end if

  end subroutine elm_drv

  !-----------------------------------------------------------------------
  subroutine elm_drv_init(bounds, &
       num_nolakec, filter_nolakec, &
       num_nolakep, filter_nolakep, &
       num_soilp  , filter_soilp, &
       canopystate_vars,  energyflux_vars)
    !
    ! !DESCRIPTION:
    ! Initialization of clm driver variables needed from previous timestep
    !
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use elm_varpar         , only : nlevsno
    use elm_varcon         , only : h2osno_max
    use landunit_varcon    , only : istice_mec
    use CanopyStateType    , only : canopystate_type
    use EnergyFluxType     , only : energyflux_type
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds
    integer               , intent(in)    :: num_nolakec       ! number of non-lake points in column filter
    integer               , intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
    integer               , intent(in)    :: num_nolakep       ! number of non-lake points in patch filter
    integer               , intent(in)    :: filter_nolakep(:) ! patch filter for non-lake points
    integer               , intent(in)    :: num_soilp         ! number of soil points in patch filter
    integer               , intent(in)    :: filter_soilp(:)   ! patch filter for soil points
    type(canopystate_type), intent(inout) :: canopystate_vars
    type(energyflux_type) , intent(inout) :: energyflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: l, c, p, f, j         ! indices
    integer :: fp, fc                  ! filter indices
    !-----------------------------------------------------------------------

    associate(                                                             &
         snl                => col_pp%snl                                   , & ! Input:  [integer  (:)   ]  number of snow layers

         h2osno             => col_ws%h2osno                , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)
         h2osoi_ice         => col_ws%h2osoi_ice            , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq         => col_ws%h2osoi_liq            , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         do_capsnow         => col_ws%do_capsnow            , & ! Output: [logical  (:)   ]  true => do snow capping
         h2osno_old         => col_ws%h2osno_old            , & ! Output: [real(r8) (:)   ]  snow water (mm H2O) at previous time step
         frac_iceold        => col_ws%frac_iceold           , & ! Output: [real(r8) (:,:) ]  fraction of ice relative to the tot water

         elai               => canopystate_vars%elai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow
         esai               => canopystate_vars%esai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow
         frac_veg_nosno     => canopystate_vars%frac_veg_nosno_patch     , & ! Output: [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         frac_veg_nosno_alb => canopystate_vars%frac_veg_nosno_alb_patch , & ! Output: [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]

         qflx_glcice        => col_wf%qflx_glcice            , & ! Output: [real(r8) (:)   ]  flux of new glacier ice (mm H2O/s) [+ = ice grows]

         qflx_glcice_diag   => col_wf%qflx_glcice_diag       , & ! Output: [real(r8) (:)   ]  flux of new glacier ice (mm H2O/s) [+ = ice grows]

         eflx_bot           => col_ef%eflx_bot              , & ! Output: [real(r8) (:)   ]  heat flux from beneath soil/ice column (W/m**2)

         cisun_z            => photosyns_vars%cisun_z_patch              , & ! Output: [real(r8) (:)   ]  intracellular sunlit leaf CO2 (Pa)
         cisha_z            => photosyns_vars%cisha_z_patch                & ! Output: [real(r8) (:)   ]  intracellular shaded leaf CO2 (Pa)
         )

      ! Initialize intracellular CO2 (Pa) parameters each timestep for use in VOCEmission
      do p = bounds%begp,bounds%endp
         cisun_z(p,:) = -999._r8
         cisha_z(p,:) = -999._r8
      end do

      do c = bounds%begc,bounds%endc
         l = col_pp%landunit(c)

         ! Save snow mass at previous time step
         h2osno_old(c) = h2osno(c)

         if (.not. use_firn_percolation_and_compaction) then
            ! Decide whether to cap snow
            if (h2osno(c) > h2osno_max) then
               do_capsnow(c) = .true.
            else
               do_capsnow(c) = .false.
            end if
         ! else, snow capping subroutine in SnowHydrologyMod
         end if

         ! Reset flux from beneath soil/ice column
         eflx_bot(c)  = 0._r8

         ! Initialize qflx_glcice everywhere, to zero.
         qflx_glcice(c) = 0._r8
         qflx_glcice_diag(c) = 0._r8

      end do

      ! Initialize fraction of vegetation not covered by snow

      do p = bounds%begp,bounds%endp
         if (veg_pp%active(p)) then
            frac_veg_nosno(p) = frac_veg_nosno_alb(p)
         else
            frac_veg_nosno(p) = 0._r8
         end if
      end do

      ! Initialize set of previous time-step variables
      ! Ice fraction of snow at previous time step

      do j = -nlevsno+1,0
         do f = 1, num_nolakec
            c = filter_nolakec(f)
            if (j >= snl(c) + 1) then
               frac_iceold(c,j) = h2osoi_ice(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))
            end if
         end do
      end do

    end associate

  end subroutine elm_drv_init

  !-----------------------------------------------------------------------
  subroutine elm_drv_patch2col (bounds, num_nolakec, filter_nolakec, &
         energyflux_vars )
    !
    ! !DESCRIPTION:
    ! Averages over all patchs for variables defined over both soil and lake
    ! to provide the column-level averages of state and flux variables
    ! defined at the patch level.
    !
    ! !USES:
    use EnergyFluxType , only : energyflux_type
    use subgridAveMod  , only : p2c
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds
    integer               , intent(in)    :: num_nolakec       ! number of column non-lake points in column filter
    integer               , intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
    type(energyflux_type) , intent(inout) :: energyflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,fc              ! indices
    integer :: num_allc          ! number of active column points
    integer :: filter_allc(bounds%endp-bounds%begp+1)    ! filter for all active column points
    ! -----------------------------------------------------------------
    associate( &
      h2ocan_patch       => veg_ws%h2ocan , &
      h2ocan_col         => col_ws%h2ocan , &
      qflx_ev_snow_patch => veg_wf%qflx_ev_snow , &
      qflx_ev_snow_col   => col_wf%qflx_ev_snow, &
      qflx_ev_soil_patch =>veg_wf%qflx_ev_soil ,&
      qflx_ev_soil_col   =>col_wf%qflx_ev_soil ,&
      qflx_ev_h2osfc_patch => veg_wf%qflx_ev_h2osfc , &
      qflx_ev_h2osfc_col   => col_wf%qflx_ev_h2osfc , &
      qflx_evap_soi_patch => veg_wf%qflx_evap_soi   , &
      qflx_evap_soi_col   => col_wf%qflx_evap_soi   , &
      qflx_evap_tot_patch => veg_wf%qflx_evap_tot   , &
      qflx_evap_tot_col   => col_wf%qflx_evap_tot   , &
      qflx_rain_grnd_patch => veg_wf%qflx_rain_grnd , &
      qflx_rain_grnd_col   => col_wf%qflx_rain_grnd , &
      qflx_snow_grnd_patch => veg_wf%qflx_snow_grnd , &
      qflx_snow_grnd_col   => col_wf%qflx_snow_grnd , &
      qflx_snwcp_liq_patch => veg_wf%qflx_snwcp_liq , &
      qflx_snwcp_liq_col   => col_wf%qflx_snwcp_liq , &
      qflx_snwcp_ice_patch => veg_wf%qflx_snwcp_ice , &
      qflx_snwcp_ice_col   => col_wf%qflx_snwcp_ice , &
      qflx_tran_veg_patch  => veg_wf%qflx_tran_veg  , &
      qflx_tran_veg_col    => col_wf%qflx_tran_veg  , &
      qflx_evap_grnd_patch => veg_wf%qflx_evap_grnd , &
      qflx_evap_grnd_col   => col_wf%qflx_evap_grnd , &
      qflx_prec_grnd_patch => veg_wf%qflx_prec_grnd , &
      qflx_prec_grnd_col   => col_wf%qflx_prec_grnd , &
      qflx_dew_grnd_patch  => veg_wf%qflx_dew_grnd  , &
      qflx_dew_grnd_col    => col_wf%qflx_dew_grnd  , &
      qflx_dirct_rain_patch => veg_wf%qflx_dirct_rain , &
      qflx_dirct_rain_col   => col_wf%qflx_dirct_rain , &
      qflx_leafdrip_patch   => veg_wf%qflx_leafdrip   , &
      qflx_leafdrip_col     => col_wf%qflx_leafdrip   , &
      qflx_sub_snow_patch   => veg_wf%qflx_sub_snow   , &
      qflx_sub_snow_col     => col_wf%qflx_sub_snow   , &
      qflx_dew_snow_patch   => veg_wf%qflx_dew_snow   , &
      qflx_dew_snow_col     => col_wf%qflx_dew_snow   , &
      qflx_irrig_patch     => veg_wf%qflx_irrig_patch , &
      qflx_irrig_col       => col_wf%qflx_irrig       , &
      qflx_evap_veg_patch  => veg_wf%qflx_evap_veg    , &
      qflx_evap_veg_col    => col_wf%qflx_evap_veg     &
      )
    ! Set up a filter for all active column points

    fc = 0
    do c = bounds%begc,bounds%endc
       if (col_pp%active(c)) then
          fc = fc + 1
          filter_allc(fc) = c
       end if
    end do
    num_allc = fc

    ! Note: lake points are excluded from many of the following
    ! averages. For some fields, this is because the field doesn't
    ! apply over lakes. However, for many others, this is because the
    ! field is computed in LakeHydrologyMod, which is called after
    ! this routine; thus, for lakes, the column-level values of these
    ! fields are explicitly set in LakeHydrologyMod. (The fields that
    ! are included here for lakes are computed elsewhere, e.g., in
    ! LakeFluxesMod.)

    ! Averaging for patch water state variables

    call p2c (bounds, num_nolakec, filter_nolakec, &
         h2ocan_patch(bounds%begp:bounds%endp), &
         h2ocan_col(bounds%begc:bounds%endc))

    ! Averaging for patch evaporative flux variables

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_ev_snow_patch(bounds%begp:bounds%endp), &
         qflx_ev_snow_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_ev_soil_patch(bounds%begp:bounds%endp), &
         qflx_ev_soil_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_ev_h2osfc_patch(bounds%begp:bounds%endp), &
         qflx_ev_h2osfc_col  (bounds%begc:bounds%endc))

    ! Averaging for patch water flux variables

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_evap_soi_patch(bounds%begp:bounds%endp), &
         qflx_evap_soi_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_evap_tot_patch(bounds%begp:bounds%endp), &
         qflx_evap_tot_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_rain_grnd_patch(bounds%begp:bounds%endp), &
         qflx_rain_grnd_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_snow_grnd_patch(bounds%begp:bounds%endp), &
         qflx_snow_grnd_col  (bounds%begc:bounds%endc))

    if (.not. use_firn_percolation_and_compaction) then
       call p2c (bounds, num_allc, filter_allc, &
            veg_wf%qflx_snwcp_liq(bounds%begp:bounds%endp), &
            col_wf%qflx_snwcp_liq(bounds%begc:bounds%endc))
            !TODO - WJS has suggested that at this point qflx_snwcp_liq_patch should
            ! now be set to nan in order to ensure that this variable is not used
            ! for the remainder of the timestep - other variables where this should
            ! occur in this routine should be examined as well

            ! For lakes, this field is initially set in LakeFluxesMod (which
            ! is called before this routine; hence it is appropriate to
            ! include lake columns in this p2c call.  However, it is later
            ! overwritten in LakeHydrologyMod, both on the patch and the column
            ! level.

       call p2c (bounds, num_allc, filter_allc, &
            veg_wf%qflx_snwcp_ice(bounds%begp:bounds%endp), &
            col_wf%qflx_snwcp_ice(bounds%begc:bounds%endc))
    end if

    call p2c (bounds, num_allc, filter_allc, &
         qflx_snwcp_liq_patch(bounds%begp:bounds%endp), &
         qflx_snwcp_liq_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_tran_veg_patch(bounds%begp:bounds%endp), &
         qflx_tran_veg_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_evap_grnd_patch(bounds%begp:bounds%endp), &
         qflx_evap_grnd_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_allc, filter_allc, &
         qflx_evap_soi_patch(bounds%begp:bounds%endp), &
         qflx_evap_soi_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_prec_grnd_patch(bounds%begp:bounds%endp), &
         qflx_prec_grnd_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_dew_grnd_patch(bounds%begp:bounds%endp), &
         qflx_dew_grnd_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_dirct_rain_patch(bounds%begp:bounds%endp), &
         qflx_dirct_rain_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_leafdrip_patch(bounds%begp:bounds%endp), &
         qflx_leafdrip_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_sub_snow_patch(bounds%begp:bounds%endp), &
         qflx_sub_snow_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_dew_snow_patch(bounds%begp:bounds%endp), &
         qflx_dew_snow_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_irrig_patch(bounds%begp:bounds%endp), &
         qflx_irrig_col  (bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_tran_veg_patch(bounds%begp:bounds%endp), &
         qflx_tran_veg_col  (bounds%begc:bounds%endc) )

    call p2c (bounds, num_nolakec, filter_nolakec, &
         qflx_evap_veg_patch(bounds%begp:bounds%endp), &
         qflx_evap_veg_col  (bounds%begc:bounds%endc))

    end associate

  end subroutine elm_drv_patch2col

  !------------------------------------------------------------------------
  subroutine write_diagnostic (bounds, wrtdia, nstep, lnd2atm_vars)
    !
    ! !DESCRIPTION:
    ! Write diagnostic surface temperature output each timestep.  Written to
    ! be fast but not bit-for-bit because order of summations can change each
    ! timestep.
    !
    ! !USES:
    use decompMod  , only : get_proc_global
    use spmdMod    , only : masterproc, npes, MPI_REAL8
    use spmdMod    , only : MPI_STATUS_SIZE, mpicom, MPI_SUM
    use shr_sys_mod, only : shr_sys_flush
    use abortutils , only : endrun
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use lnd2atmType, only : lnd2atm_type
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in) :: bounds
    logical            , intent(in) :: wrtdia     !true => write diagnostic
    integer            , intent(in) :: nstep      !model time step
    type(lnd2atm_type) , intent(in) :: lnd2atm_vars
    !
    ! !REVISION HISTORY:
    ! Created by Mariana Vertenstein
    !
    ! !LOCAL VARIABLES:
    integer :: p                       ! loop index
    integer :: numg                    ! total number of gridcells across all processors
    integer :: ier                     ! error status
    real(r8):: psum                    ! partial sum of ts
    real(r8):: tsum                    ! sum of ts
    real(r8):: tsxyav                  ! average ts for diagnostic output
    integer :: status(MPI_STATUS_SIZE) ! mpi status
    character(len=256)   :: dateTimeString
    !------------------------------------------------------------------------

    call get_proc_global(ng=numg)

    if (wrtdia) then

       call t_barrierf('sync_write_diag', mpicom)
       psum = sum(lnd2atm_vars%t_rad_grc(bounds%begg:bounds%endg))
       call mpi_reduce(psum, tsum, 1, MPI_REAL8, MPI_SUM, 0, mpicom, ier)
       if (ier/=0) then
          write(iulog,*) 'write_diagnostic: Error in mpi_reduce()'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       if (masterproc) then
          tsxyav = tsum / numg
          write(iulog,1000) nstep, tsxyav
          call shr_sys_flush(iulog)
       end if

    else

#ifndef CPL_BYPASS
       call get_curr_time_string(dateTimeString)
       if (masterproc) then
          write(iulog,*)'   Completed timestep: ',trim(dateTimeString)
          call shr_sys_flush(iulog)
       end if
#endif

    endif

1000 format (1x,'nstep = ',i10,'   TS = ',f21.15)

  end subroutine write_diagnostic

end module elm_driver
