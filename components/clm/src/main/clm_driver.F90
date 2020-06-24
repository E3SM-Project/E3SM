module clm_driver
  !
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides the main CLM driver physics calling sequence.  Most
  ! computations occurs over ``clumps'' of gridcells (and associated subgrid
  ! scale entities) assigned to each MPI process. Computation is further
  ! parallelized by looping over clumps on each process using shared memory OpenMP.
  !
  ! !USES:
  use clm_varcon , only : zsoi
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_sys_mod            , only : shr_sys_flush
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varctl             , only : wrtdia, iulog, create_glacier_mec_landunit, use_fates
  use clm_varpar             , only : nlevtrc_soil, nlevsoi, nlevdecomp_full
  use clm_varctl             , only : wrtdia, iulog, create_glacier_mec_landunit, use_fates, use_betr
  use clm_varctl             , only : use_cn, use_lch4, use_voc, use_noio, use_c13, use_c14
  use clm_varctl             , only : use_erosion
  use clm_time_manager       , only : get_step_size, get_curr_date, get_ref_date, get_nstep, is_beg_curr_day, get_curr_time_string
  use clm_varpar             , only : nlevsno, nlevgrnd, crop_prog
  use spmdMod                , only : masterproc, mpicom
  use decompMod              , only : get_proc_clumps, get_clump_bounds, get_proc_bounds, bounds_type
  use decompMod              , only : get_clump_bounds_gpu
  use filterMod              , only : filter, filter_inactive_and_active
  use histFileMod            , only : hist_update_hbuf, hist_htapes_wrapup
  use restFileMod            , only : restFile_write, restFile_filename
  use abortutils             , only : endrun
  !
  use dynSubgridDriverMod    , only : dynSubgrid_driver,dynSubgrid_wrapup_weight_changes
  use dynSubgridDriverMod    , only : prior_weights, column_state_updater, patch_state_updater
  use BalanceCheckMod        , only : BeginColWaterBalance, ColWaterBalanceCheck
  use BalanceCheckMod        , only : BeginGridWaterBalance, GridBalanceCheck
  !
  use CanopyTemperatureMod   , only : CanopyTemperature ! (formerly Biogeophysics1Mod)
  use SoilTemperatureMod     , only : SoilTemperature
  use LakeTemperatureMod     , only : LakeTemperature
  !
  use BareGroundFluxesMod    , only : BareGroundFluxes
  use CanopyFluxesMod        , only : CanopyFluxes
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
  !clm_interface
  use EcosystemDynMod      , only : EcosystemDynNoLeaching1, EcosystemDynNoLeaching2

  use EcosystemDynMod      , only : EcosystemDynLeaching
  use VegStructUpdateMod   , only : VegStructUpdate
  use AnnualUpdateMod      , only : AnnualUpdate
  use EcosystemBalanceCheckMod      , only : BeginColCBalance, BeginColNBalance, ColCBalanceCheck, ColNBalanceCheck
  use EcosystemBalanceCheckMod      , only : BeginColPBalance, ColPBalanceCheck
  use EcosystemBalanceCheckMod      , only : BeginGridCBalanceBeforeDynSubgridDriver
  use EcosystemBalanceCheckMod      , only : BeginGridNBalanceBeforeDynSubgridDriver
  use EcosystemBalanceCheckMod      , only : BeginGridPBalanceBeforeDynSubgridDriver
  use EcosystemBalanceCheckMod      , only : EndGridCBalanceAfterDynSubgridDriver
  use EcosystemBalanceCheckMod      , only : EndGridNBalanceAfterDynSubgridDriver
  use EcosystemBalanceCheckMod      , only : EndGridPBalanceAfterDynSubgridDriver
  use VerticalProfileMod   , only : decomp_vertprofiles,decomp_vertprofiles2
  use FireMod              , only : FireInterp
  use SatellitePhenologyMod  , only : SatellitePhenology, interpMonthlyVeg
  use ndepStreamMod          , only : ndep_interp
  use pdepStreamMod          , only : pdep_interp
  use ActiveLayerMod         , only : alt_calc
  use CH4Mod                 , only : CH4
  use DUSTMod                , only : DustDryDep, DustEmission
  use VOCEmissionMod         , only : VOCEmission
  use FatesBGCDynMod         , only : FatesBGCDyn
  !
  use filterMod              , only : setFilters
  !
  use atm2lndMod             , only : downscale_forcings
  use lnd2atmMod             , only : lnd2atm
  use lnd2glcMod             , only : lnd2glc_type, update_lnd2glc_GPU
  !
  use seq_drydep_mod_elm     , only : n_drydep, drydep_method, DD_XLND
  use DryDepVelocity         , only : depvel_compute
  !
  use DaylengthMod           , only : UpdateDaylength, first_step
  use perf_mod
  use SedYieldMod            , only : SoilErosion
  !
  use clm_instMod            , only : ch4_vars, ep_betr
  use clm_instMod            , only : carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars
  use clm_instMod            , only : carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars
  use clm_instMod            , only : nitrogenstate_vars
  use clm_instMod            , only : nitrogenflux_vars
  use clm_instMod            , only : phosphorusstate_vars
  use clm_instMod            , only : phosphorusflux_vars
  use clm_instMod            , only : sedflux_vars
  use clm_instMod            , only : crop_vars
  use clm_instMod            , only : cnstate_vars
  use clm_instMod            , only : dust_vars
  use clm_instMod            , only : vocemis_vars
  use clm_instMod            , only : drydepvel_vars
  use clm_instMod            , only : aerosol_vars
  use clm_instMod            , only : canopystate_vars
  use clm_instMod            , only : energyflux_vars
  use clm_instMod            , only : frictionvel_vars
  use clm_instMod            , only : lakestate_vars
  use clm_instMod            , only : photosyns_vars
  use clm_instMod            , only : soilstate_vars
  use clm_instMod            , only : soilhydrology_vars
  use clm_instMod            , only : solarabs_vars
  use clm_instMod            , only : soilhydrology_vars
  use clm_instMod            , only : surfalb_vars
  use clm_instMod            , only : surfrad_vars
  use clm_instMod            , only : temperature_vars
  use clm_instMod            , only : waterflux_vars
  use clm_instMod            , only : waterstate_vars
  use clm_instMod            , only : atm2lnd_vars
  use clm_instMod            , only : lnd2atm_vars
  use clm_instMod            , only : glc2lnd_vars
  use clm_instMod            , only : lnd2glc_vars
  use clm_instMod            , only : soil_water_retention_curve
  use clm_instMod            , only : chemstate_vars
  use clm_instMod            , only : alm_fates
  use clm_instMod            , only : PlantMicKinetics_vars
  use tracer_varcon          , only : is_active_betr_bgc
  use CNEcosystemDynBetrMod  , only : CNEcosystemDynBetr, CNFluxStateBetrSummary
  use UrbanParamsType        , only : urbanparams_vars

  use GridcellType             , only : grc_pp
  use GridcellDataType         , only : grc_cs, c13_grc_cs, c14_grc_cs
  use GridcellDataType         , only : grc_cf, c13_grc_cf, c14_grc_cf
  use GridcellDataType         , only : grc_nf, grc_pf, grc_ef, grc_wf, grc_ws
  use GridcellDataType         , only : grc_es, grc_ns, grc_ps
  use TopounitDataType         , only : top_as, top_af
  use TopounitType             , only : top_pp
  use LandunitType             , only : lun_pp
  use ColumnType               , only : col_pp
  use ColumnDataType           , only : col_es, col_ef, col_ws, col_wf
  use ColumnDataType           , only : col_cs, c13_col_cs, c14_col_cs
  use ColumnDataType           , only : col_cf, c13_col_cf, c14_col_cf
  use ColumnDataType           , only : col_ns, col_nf
  use ColumnDataType           , only : col_ps, col_pf
  use VegetationType           , only : veg_pp
  use VegetationDataType       , only : veg_es, veg_ws, veg_wf
  use VegetationDataType       , only : veg_cs, c13_veg_cs, c14_veg_cs
  use VegetationDataType       , only : c13_veg_cf, c14_veg_cf
  use VegetationDataType       , only : veg_ns, veg_nf
  use VegetationDataType       , only : veg_ps, veg_pf
  use VegetationDataType       , only : veg_cf, veg_ef
  use VegetationPropertiesType , only : veg_vp
  use PhotosynthesisMod        , only : params_inst
  use SharedParamsMod          , only : ParamsShareInst
  use CH4Mod                   , only : CH4ParamsInst
  use CNDecompCascadeConType   , only : decomp_cascade_con
  use DecompCascadeBGCMod      , only : DecompBGCParamsInst
  use GapMortalityMod          , only : cngapmortparamsinst
  use DecompCascadeCNMod       , only : DecompCNParamsInst
  use NitrifDenitrifMod        , only : NitrifDenitrifParamsInst
  use SoilLittDecompMod        , only : cndecompparamsinst
  use AllocationMod            , only : AllocParamsInst
  use LandunitDataType         , only : lun_ef, lun_es, lun_ws, lun_wf
  use glc2lndMod                , only : glc2lnd_vars_update_glc2lnd_acc


  !----------------------------------------------------------------------------
  ! bgc interface & pflotran:
  use clm_varctl             , only : use_clm_interface
  use clm_instMod            , only : clm_interface_data
  use clm_interface_funcsMod , only : get_clm_data
  ! (1) clm_bgc through interface
  use clm_varctl             , only : use_clm_bgc
  use clm_interface_funcsMod , only : clm_bgc_run, update_bgc_data_clm2clm
  ! (2) pflotran
  use clm_time_manager            , only : nsstep, nestep
  use clm_varctl                  , only : use_pflotran, pf_cmode, pf_hmode, pf_tmode
  use clm_interface_funcsMod      , only : update_bgc_data_pf2clm, update_th_data_pf2clm
  use clm_interface_pflotranMod   , only : clm_pf_run, clm_pf_write_restart
  use clm_interface_pflotranMod   , only : clm_pf_finalize
  !----------------------------------------------------------------------------
  use WaterBudgetMod              , only : WaterBudget_Reset, WaterBudget_Run, WaterBudget_Accum, WaterBudget_Print
  use WaterBudgetMod              , only : WaterBudget_SetBeginningMonthlyStates
  use WaterBudgetMod              , only : WaterBudget_SetEndingMonthlyStates
  use clm_varctl                  , only : do_budgets, budget_inst, budget_daily, budget_month
  use clm_varctl                  , only : budget_ann, budget_ltann, budget_ltend
  use clm_varctl                  , only : spinup_state, nyears_ad_carbon_only, spinup_mortality_factor
  use clm_varctl                  , only : carbon_only , carbonphosphorus_only, carbonnitrogen_only
  use decompMod                   , only : clumps, procinfo
  use domainMod                   , only : ldomain
  use verificationMod
  use update_accMod
  use cudafor
  use Method_procs_acc
  use histGPUMod
  use NitrogenDynamicsMod, only : CNNDynamicsParamsInst
  use dynSubgridControlMod , only : dyn_subgrid_control_inst
  use dynInitColumnsMod    , only : initialize_new_columns
  use dynConsBiogeophysMod , only : dyn_hwcontent_init, dyn_hwcontent_final
  use dynConsBiogeochemMod , only : dyn_cnbal_patch, dyn_cnbal_column
  use reweightMod          , only : reweight_wrapup
  use subgridWeightsMod    , only : compute_higher_order_weights, set_subgrid_diagnostic_fields
  use subgridWeightsMod    , only : subgrid_weights_diagnostics
  use CarbonStateUpdate1Mod   , only : CarbonStateUpdateDynPatch
  use NitrogenStateUpdate1Mod   , only : NitrogenStateUpdateDynPatch
  use PhosphorusStateUpdate1Mod     , only : PhosphorusStateUpdateDynPatch
  use dynUpdateModAcc , only : update_column_state_acc, update_patch_stateAcc
  use dynUpdateModAcc , only : patch_set_new_weightsAcc, column_set_new_weightsAcc
  use dynPriorWeightsMod , only : set_prior_weights_acc
  use dynPatchStateUpdaterMod  , only : patch_state_set_old_weights_acc
  use dynColumnStateUpdaterMod , only : column_state_set_old_weights_acc
  use ForcingUpdateMod         , only : update_forcings_CPLBYPASS
  use clm_varctl

  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_drv            ! Main clm driver
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: clm_drv_patch2col
  private :: clm_drv_init      ! Initialization of variables needed from previous timestep
  private :: write_diagnostic  ! Write diagnostic information to log file
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine clm_drv(step_count, rstwr, nlend, rdate)
    !
    ! !DESCRIPTION:
    !
    ! First phase of the clm driver calling the clm physics. An outline of
    ! the calling tree is given in the description of this module.
    !
    ! !USES:
    !
    use shr_const_mod   , only : SHR_CONST_TKFRZ
    use clm_time_manager
    use timeinfoMod
    use histfileMod   , only : clmptr_rs, tape, clmptr_ra
    use accumulGPUMod
    use decompMod     , only : init_proc_clump_info, gpu_clumps, gpu_procinfo
    use clm_varorb
    use shr_orb_mod_elm
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(inout)  :: step_count
    logical           :: rstwr       ! true => write restart file this step
    logical,         intent(in) :: nlend       ! true => end of run on this step
    character(len=*),intent(in) :: rdate       ! restart file time stamp for name
    !
    ! !LOCAL VARIABLES:
    integer              :: nc, c, p, l, g, fc, j   ! indices
    integer              :: nclumps                 ! number of clumps on this processor
    character(len=256)   :: filer                   ! restart file name
    integer              :: ier                     ! error code
    character(len=256)   :: dateTimeString
    type(bounds_type)    :: bounds_clump
    type(bounds_type)    :: bounds_proc
    integer   ::  mygpu, ngpus, cid, fp, idle
    logical :: found_thawlayer
    integer :: k_frz
    real*8    :: sto
    integer , parameter :: gpu = 1, numdays = 1
    #if _CUDA
    integer(kind=cuda_count_kind) :: heapsize,free1,free2,total
    integer  :: istat, val
    #endif
    !-----------------------------------------------------------------------
    call get_curr_time_string(dateTimeString)
    if (masterproc) then
       write(iulog,*)'Beginning timestep   : ',trim(dateTimeString)
       call shr_sys_flush(iulog)
    endif

    idle = -1
    if(step_count == -24 ) idle = 0
    do while (idle == 0)
        idle = idle + 0
        call sleep(1)
    end do
    #if _CUDA
    istat = cudaDeviceGetLimit(heapsize, cudaLimitMallocHeapSize)
    print *, "SETTING Heap Limit from", heapsize
    heapsize = 10_8*1024_8*1024_8
    print *, "TO:",heapsize
    istat = cudaDeviceSetLimit(cudaLimitMallocHeapSize,heapsize)
    istat = cudaMemGetInfo(free1, total)
    print *, "Free1:",free1
    #endif
    ! Determine processor bounds and clumps for this processor
    call get_proc_bounds(bounds_proc)
    nclumps = get_proc_clumps()
    print *, "step:", step_count
    if(step_count == 0 ) then
      print *, "transferring data to GPU"
       call init_proc_clump_info()
       !$acc update device( &
       !$acc        spinup_state            &
       !$acc       , nyears_ad_carbon_only   &
       !$acc       , spinup_mortality_factor &
       !$acc       , carbon_only &
       !$acc       , carbonphosphorus_only &
       !$acc       , carbonnitrogen_only &
       !$acc       ,use_crop            &
       !$acc       ,use_snicar_frc      &
       !$acc       ,use_snicar_ad       &
       !$acc       ,use_vancouver       &
       !$acc       ,use_mexicocity      &
       !$acc       ,use_noio            &
       !$acc       ,use_var_soil_thick  &
       !$acc       ,NFIX_PTASE_plant &
       !$acc       ,tw_irr &
       !$acc       ,use_erosion &
       !$acc       ,ero_ccycle  &
       !$acc       ,anoxia &
       !$acc       , glc_do_dynglacier &
       !$acc       , all_active &
       !$acc       , co2_ppmv &
       !$acc       , const_climate_hist &
       !$acc     )
       !$acc update device(first_step, nlevgrnd, eccen, obliqr, lambm0, mvelpp )
       call update_acc_variables()
       !
       !$acc enter data copyin(filter(:), &
       !$acc filter_inactive_and_active(:), gpu_clumps(:), gpu_procinfo )
       !
        !$acc enter data copyin(&
        !$acc aerosol_vars     , &
        !$acc AllocParamsInst  , &
        !$acc atm2lnd_vars     , &
        !$acc c13_col_cf     , &
        !$acc c13_col_cs     , &
        !$acc c13_grc_cf     , &
        !$acc c13_veg_cf     , &
        !$acc c13_veg_cs     , &
        !$acc c14_col_cf     , &
        !$acc c14_col_cs     , &
        !$acc c14_grc_cf     , &
        !$acc c14_veg_cf     , &
        !$acc c14_veg_cs     , &
        !$acc canopystate_vars, &
        !$acc CH4ParamsInst     , &
        !$acc ch4_vars          , &
        !$acc CNDecompParamsInst     , &
        !$acc CNGapMortParamsInst     , &
        !$acc CNNDynamicsParamsInst     , &
        !$acc cnstate_vars     , &
        !$acc column_state_updater     , &
        !$acc col_cf     , &
        !$acc col_cs     , &
        !$acc col_ef     , &
        !$acc col_es     , &
        !$acc col_nf     , &
        !$acc col_ns     , &
        !$acc col_pf     , &
        !$acc col_pp     , &
        !$acc col_ps     , &
        !$acc col_wf     , &
        !$acc col_ws     , &
        !$acc crop_vars     , &
        !$acc dyn_subgrid_control_inst , &
        !$acc subgrid_weights_diagnostics, &
        !$acc DecompBGCParamsInst     , &
        !$acc DecompCNParamsInst     , &
        !$acc decomp_cascade_con     , &
        !$acc drydepvel_vars     , &
        !$acc dust_vars     , &
        !$acc energyflux_vars     , &
        !$acc frictionvel_vars     , &
        !$acc glc2lnd_vars     , &
        !$acc grc_cf     , &
        !$acc grc_cs     , &
        !$acc grc_ef     , &
        !$acc grc_es     , &
        !$acc grc_nf     , &
        !$acc grc_ns     , &
        !$acc grc_pf     , &
        !$acc grc_pp     , &
        !$acc grc_ps     , &
        !$acc grc_wf     , &
        !$acc grc_ws     , &
        !$acc lakestate_vars , &
        !$acc ldomain  ,&
        !$acc lnd2glc_vars   , &
        !$acc lnd2atm_vars   , &
        !$acc lun_ef     , &
        !$acc lun_es     , &
        !$acc lun_pp     , &
        !$acc lun_ws     , &
        !$acc NitrifDenitrifParamsInst     , &
        !$acc ParamsShareInst     , &
        !$acc params_inst     , &
        !$acc patch_state_updater     , &
        !$acc photosyns_vars     , &
        !$acc prior_weights     , &
        !$acc sedflux_vars     , &
        !$acc soilhydrology_vars     , &
        !$acc soilstate_vars     , &
        !$acc solarabs_vars     , &
        !$acc surfalb_vars     , &
        !$acc surfrad_vars     , &
        !$acc top_af     , &
        !$acc top_as     , &
        !$acc top_pp     , &
        !$acc urbanparams_vars     , &
        !$acc veg_cf     , &
        !$acc veg_cs     , &
        !$acc veg_ef     , &
        !$acc veg_es     , &
        !$acc veg_nf     , &
        !$acc veg_ns     , &
        !$acc veg_pf     , &
        !$acc veg_pp     , &
        !$acc veg_ps     , &
        !$acc veg_vp     , &
        !$acc veg_wf     , &
        !$acc veg_ws     , &
        !$acc vocemis_vars &
        !$acc   )
        !$acc enter data copyin(tape_gpu,clmptr_ra,clmptr_rs)
        !$acc enter data copyin( doalb, declinp1, declin )
        
        #if _CUDA
              istat = cudaMemGetInfo(free2, total)
              print *, "Transferred:", free1-free2
              print *, "Total:",total
              print *, "Free:", free2
        #endif
      end if
     
    if (do_budgets) call WaterBudget_Reset()

    ! ============================================================================
    ! Specified phenology
    ! ============================================================================

    if (.not.use_fates) then
       if (use_cn) then
          ! For dry-deposition need to call CLMSP so that mlaidiff is obtained
          if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
             call t_startf('interpMonthlyVeg')
             print *, "TURNED OFF INTERPMONTHLYVEG"
             !call interpMonthlyVeg(bounds_proc, canopystate_vars)
             call t_stopf('interpMonthlyVeg')
          endif

       else
          ! Determine weights for time interpolation of monthly vegetation data.
          ! This also determines whether it is time to read new monthly vegetation and
          ! obtain updated leaf area index [mlai1,mlai2], stem area index [msai1,msai2],
          ! vegetation top [mhvt1,mhvt2] and vegetation bottom [mhvb1,mhvb2]. The
          ! weights obtained here are used in subroutine SatellitePhenology to obtain time
          ! interpolated values.
          if (doalb .or. ( n_drydep > 0 .and. drydep_method == DD_XLND )) then
             print *, "TURNED OFF INTERPMONTHLYVEG"
              call t_startf('interpMonthlyVeg')
             !call interpMonthlyVeg(bounds_proc, canopystate_vars)
             call t_stopf('interpMonthlyVeg')
          end if
       end if
    end if

    !$acc serial default(present)
    call increment_time_vars()
    call shr_orb_decl(thiscalday_mod , eccen, mvelpp, lambm0, obliqr, declin  , eccf )
    call shr_orb_decl(nextsw_cday_mod, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )
    !$acc end serial
    
   print *, "first loop!"
  !$acc parallel default(present)
   !$acc loop independent gang private(nc, bounds_clump)
    do nc = 1, nclumps
      call get_clump_bounds_gpu(nc, bounds_clump)
      !
      if (nstep_mod > 1) then
        call update_forcings_CPLBYPASS(bounds_clump, atm2lnd_vars, int(dtime_mod),&
                  thiscalday_mod, secs_curr, year_curr, mon_curr, nstep_mod)
      end if
       ! ==================================================================================
       ! Determine decomp vertical profiles
       !
       ! These routines (alt_calc & decomp_vertprofiles) need to be called before
       ! pftdyn_cnbal, and it appears that they need to be called before pftdyn_interp and
       ! the associated filter updates, too (otherwise we get a carbon balance error)
       ! ==================================================================================
      call alt_calc(filter(nc)%num_soilc, filter(nc)%soilc, canopystate_vars, &
                      year_curr, mon_curr,day_curr,secs_curr, dtime_mod)
      if (use_cn) then
          !  Note (WJS, 6-12-13): Because of this routine's placement in the driver sequence
          !  (it is called very early in each timestep, before weights are adjusted and
          !  filters are updated), it may be necessary for this routine to compute values over
          !  inactive as well as active points (since some inactive points may soon become
          !  active) - so that's what is done now. Currently, it seems to be okay to do this,
          !  because the variables computed here seem to only depend on quantities that are
          !  valid over inactive as well as active points.
          call decomp_vertprofiles2(bounds_clump, &
               filter_inactive_and_active(nc), &
               soilstate_vars, canopystate_vars, cnstate_vars)
       end if
       ! ============================================================================
       ! Zero fluxes for transient land cover
       ! ============================================================================

       call BeginGridWaterBalance(bounds_clump,               &
            filter(nc)%num_nolakec, filter(nc)%nolakec,       &
            filter(nc)%num_lakec, filter(nc)%lakec,           &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
            soilhydrology_vars)

      call elm_zero_fluxes(bounds_clump)

      if (use_cn) then
        call vegcs_Summary_acc(veg_cs,bounds_clump, &
             filter(nc)%num_soilc, filter(nc)%soilc, &
             filter(nc)%num_soilp, filter(nc)%soilp, col_cs)
        call colcs_Summary_acc(col_cs,bounds_clump, &
             filter(nc)%num_soilc, filter(nc)%soilc)

        call vegns_Summary_acc(veg_ns,bounds_clump, &
             filter(nc)%num_soilc, filter(nc)%soilc, &
             filter(nc)%num_soilp, filter(nc)%soilp, col_ns)

        call colns_Summary_acc(col_ns, bounds_clump, &
             filter(nc)%num_soilc, filter(nc)%soilc)

        call vegps_Summary_acc(veg_ps, bounds_clump, &
             filter(nc)%num_soilc, filter(nc)%soilc, &
             filter(nc)%num_soilp, filter(nc)%soilp, col_ps)
        call colps_Summary_acc(col_ps, bounds_clump, &
             filter(nc)%num_soilc, filter(nc)%soilc)

        call BeginGridCBalanceBeforeDynSubgridDriver(bounds_clump, col_cs, grc_cs)
        call BeginGridNBalanceBeforeDynSubgridDriver(bounds_clump)
        call BeginGridPBalanceBeforeDynSubgridDriver(bounds_clump)

     end if

   end do
  !$acc end parallel
     
  !$acc parallel default(present)

      !$acc loop independent gang private(nc, bounds_clump)
      do nc = 1, nclumps
          call get_clump_bounds_gpu(nc,bounds_clump)

        ! ============================================================================
        ! Update subgrid weights with dynamic landcover (prescribed transient patches,
        ! and or dynamic landunits), and do related adjustments. Note that this
        ! call needs to happen outside loops over nclumps.
        ! ============================================================================

       call dyn_hwcontent_init(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_lakec, filter(nc)%lakec, &
            urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
            energyflux_vars)

       call set_prior_weights_acc(prior_weights, bounds_clump)
       call patch_state_set_old_weights_acc  (patch_state_updater,bounds_clump)
       call column_state_set_old_weights_acc(column_state_updater,bounds_clump)

       if (create_glacier_mec_landunit) then
          call glc2lnd_vars_update_glc2lnd_acc(glc2lnd_vars ,bounds_clump)
       end if

       ! Everything following this point in this loop only needs to be called if we have
       ! actually changed some weights in this time step. This is also required in the
       ! first time step of the run to update filters to reflect state of CISM
       ! (particularly mask that is past through coupler).
       call dynSubgrid_wrapup_weight_changes(bounds_clump, glc2lnd_vars)

       call patch_set_new_weightsAcc (patch_state_updater ,bounds_clump)
       call column_set_new_weightsAcc(column_state_updater,bounds_clump, nc)

       call set_subgrid_diagnostic_fields(bounds_clump)

       call initialize_new_columns(bounds_clump, &
                prior_weights%cactive(bounds_clump%begc:bounds_clump%endc) , &
                soilhydrology_vars)

       call dyn_hwcontent_final(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_lakec, filter(nc)%lakec, &
            urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
            energyflux_vars, dtime_mod)

       if (use_cn) then
          call dyn_cnbal_patch(bounds_clump, &
               filter_inactive_and_active(nc)%num_soilp, filter_inactive_and_active(nc)%soilp, &
               filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc, &
               prior_weights, &
               patch_state_updater, &
               canopystate_vars, photosyns_vars, cnstate_vars, &
               veg_cs, c13_veg_cs, c14_veg_cs, &
               veg_ns, veg_ps, dtime_mod)

          call CarbonStateUpdateDynPatch(bounds_clump, &
               filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc, &
               grc_cs, grc_cf, col_cs, col_cf, dtime_mod)

          call NitrogenStateUpdateDynPatch(bounds_clump, &
               filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc,dtime_mod)

          call PhosphorusStateUpdateDynPatch(bounds_clump, &
               filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc,dtime_mod)

          call dyn_cnbal_column(bounds_clump, nc, column_state_updater, &
               col_cs, c13_col_cs, c14_col_cs, &
               col_ns, col_ps )
       end if


    end do
    !$acc end parallel
    
    !call dynSubgrid_driver(bounds_proc,                                      &
  !   urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
  !   energyflux_vars,  canopystate_vars, photosyns_vars, cnstate_vars,     &
  !   veg_cs, c13_veg_cs, c14_veg_cs,         &
  !   col_cs, c13_col_cs, c14_col_cs, col_cf,  &
  !   grc_cs, grc_cf, glc2lnd_vars, crop_vars, dtime_mod )

    if (.not. use_fates)then
       if (use_cn) then
          if (nstep_mod < 2 )then
             if (masterproc) then
                write(iulog,*) '--WARNING-- skipping CN balance check for first timestep'
             end if
          else
            !$acc parallel default(present)

              !$acc loop independent gang private(nc, bounds_clump)

             do nc = 1,nclumps
                call get_clump_bounds_gpu(nc, bounds_clump)
                call vegcs_summary_acc(veg_cs, bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     filter(nc)%num_soilp, filter(nc)%soilp, col_cs)

                call colcs_Summary_acc(col_cs, bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc)

                call vegns_Summary_acc(veg_ns,bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     filter(nc)%num_soilp, filter(nc)%soilp, col_ns)

                call colns_Summary_acc(col_ns,bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc)

                call vegps_Summary_acc(veg_ps,bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     filter(nc)%num_soilp, filter(nc)%soilp, col_ps)

                call colps_Summary_acc(col_ps,bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc)

                call EndGridCBalanceAfterDynSubgridDriver(bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     col_cs, grc_cs, dtime_mod)

                call EndGridNBalanceAfterDynSubgridDriver(bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     dtime_mod)

                call EndGridPBalanceAfterDynSubgridDriver(bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     dtime_mod)

             end do
            !$acc end parallel
          end if
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
    !$acc parallel default(present)

      !$acc loop independent gang private(nc, bounds_clump)
     do nc = 1,nclumps
       call get_clump_bounds_gpu(nc, bounds_clump)

       call BeginColWaterBalance(bounds_clump,                &
            filter(nc)%num_nolakec, filter(nc)%nolakec,       &
            filter(nc)%num_lakec, filter(nc)%lakec,           &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
            soilhydrology_vars)

       if (use_cn) then
          ! call veg summary before col summary, for p2c
          call vegcs_Summary_acc(veg_cs,bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, col_cs)
          call colcs_Summary_acc(col_cs,bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc)

          call vegns_Summary_acc(veg_ns,bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, col_ns)
          call colns_Summary_acc(col_ns,bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc)

          call vegps_Summary_acc(veg_ps,bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, col_ps)

          call colps_Summary_acc(col_ps, bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc)

          call BeginColCBalance(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               col_cs)

          call BeginColNBalance(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc)

          call BeginColPBalance(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc)
       end if

       if (do_budgets) then
          call WaterBudget_SetBeginningMonthlyStates(bounds_clump , nstep_mod, day_curr, secs_curr &
                                        , day_prev, secs_prev)
       endif

    end do
   !$acc end parallel

    ! ============================================================================
    ! Update dynamic N deposition field, on albedo timestep
    ! currently being done outside clumps loop, but no reason why it couldn't be
    ! re-written to go inside.
    ! ============================================================================

#ifndef CPL_BYPASS
    if (use_cn) then
       call t_startf('ndep_interp')
       ! PET: switching CN timestep
       call ndep_interp(bounds_proc, atm2lnd_vars)
       call FireInterp(bounds_proc)
       call t_stopf('ndep_interp')
    end if

    ! ============================================================================
    ! Update dynamic P deposition field, on albedo timestep
    ! currently being done outside clumps loop, but no reason why it couldn't be
    ! re-written to go inside.
    ! ============================================================================

    if (use_cn) then
       ! PET: switching CN timestep
       call pdep_interp(bounds_proc, atm2lnd_vars)
    end if

#endif

    ! ============================================================================
    ! Initialize variables from previous time step, downscale atm forcings, and
    ! Determine canopy interception and precipitation onto ground surface.
    ! Determine the fraction of foliage covered by water and the fraction
    ! of foliage that is dry and transpiring. Initialize snow layer if the
    ! snow accumulation exceeds 10 mm.
    ! ============================================================================


    print *, "main loop"
  !$acc parallel vector_length(32) default(present)

    !$acc loop independent gang private(nc, bounds_clump)
    do nc = 1,nclumps

       call get_clump_bounds_gpu(nc, bounds_clump)
       
       call UpdateDaylength(bounds_clump, declin)
       ! Initialze variables needed for new driver time step
       call clm_drv_init(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_nolakep, filter(nc)%nolakep, &
            filter(nc)%num_soilp  , filter(nc)%soilp,   &
            canopystate_vars, photosyns_vars, &
            col_wf, col_ef)

       call downscale_forcings(bounds_clump, &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c, &
            atm2lnd_vars)
       ! ============================================================================
       ! Canopy Hydrology
       ! (1) water storage of intercepted precipitation
       ! (2) direct throughfall and canopy drainage of precipitation
       ! (3) fraction of foliage covered by water and the fraction is dry and transpiring
       ! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
       ! ============================================================================
       call CanopyHydrology(bounds_clump, &
                  filter(nc)%num_nolakec, filter(nc)%nolakec, &
                  filter(nc)%num_nolakep, filter(nc)%nolakep, &
                  atm2lnd_vars, canopystate_vars, &
                  aerosol_vars, dtime_mod)

       ! ============================================================================
       ! Surface Radiation
       ! ============================================================================

       ! Surface Radiation primarily for non-urban columns

       ! Most of the surface radiation calculations are agnostic to the forest-model
       ! but the calculations of the fractions of sunlit and shaded canopies
       ! are specific, calculate them first.
       ! The nourbanp filter is set in dySubgrid_driver (earlier in this call)
       ! over the patch index range defined by bounds_clump%begp:bounds_proc%endp

          call CanopySunShadeFractions(filter(nc)%num_nourbanp, filter(nc)%nourbanp,    &
                                       atm2lnd_vars, surfalb_vars, canopystate_vars,    &
                                       solarabs_vars)

          call SurfaceRadiation(bounds_clump,                        &
                   filter(nc)%num_nourbanp, filter(nc)%nourbanp,  &
                   filter(nc)%num_urbanp, filter(nc)%urbanp    ,  &
                   filter(nc)%num_urbanc, filter(nc)%urbanc,      &
                   atm2lnd_vars, canopystate_vars, surfalb_vars, &
                   solarabs_vars, surfrad_vars, dtime_mod, secs_curr)

       ! Surface Radiation for only urban columns
       call UrbanRadiation(bounds_clump,                                       &
                 filter(nc)%num_nourbanl, filter(nc)%nourbanl,                      &
                 filter(nc)%num_urbanl, filter(nc)%urbanl,                          &
                 filter(nc)%num_urbanc, filter(nc)%urbanc,                          &
                 filter(nc)%num_urbanp, filter(nc)%urbanp,                          &
                 atm2lnd_vars, urbanparams_vars, &
                 solarabs_vars, surfalb_vars, energyflux_vars)


       ! ============================================================================
       ! Determine leaf temperature and surface fluxes based on ground
       ! temperature from previous time step.
       ! ============================================================================
       call CanopyTemperature(bounds_clump,                                   &
                   filter(nc)%num_nolakec, filter(nc)%nolakec,                       &
                   filter(nc)%num_nolakep, filter(nc)%nolakep,                       &
                   atm2lnd_vars, canopystate_vars, soilstate_vars, frictionvel_vars, &
                   energyflux_vars)
       ! Determine fluxes
       ! =======================================================================
       call col_wf_reset(col_wf,filter(nc)%num_nolakec, filter(nc)%nolakec)
       !========================================================================
       call BareGroundFluxes(bounds_clump,                               &
                filter(nc)%num_nolakeurbanp, filter(nc)%nolakeurbanp,    &
                atm2lnd_vars, canopystate_vars, soilstate_vars,          &
                frictionvel_vars, ch4_vars )

       ! non-bareground fluxes for all patches except lakes and urban landunits
       ! Calculate canopy temperature, latent and sensible fluxes from the canopy,
       ! and leaf water change by evapotranspiration
       !call t_startf('canflux')
       call CanopyFluxes(bounds_clump,                                                   &
                filter(nc)%num_nolakeurbanp, filter(nc)%nolakeurbanp,                        &
                atm2lnd_vars, canopystate_vars, cnstate_vars, energyflux_vars,               &
                frictionvel_vars, soilstate_vars, solarabs_vars, surfalb_vars,               &
                ch4_vars, photosyns_vars, dtime_mod, year_curr, mon_curr, day_curr, secs_curr)

       ! Fluxes for all urban landunits
       call UrbanFluxes(bounds_clump,                            &
                filter(nc)%num_nourbanl, filter(nc)%nourbanl,    &
                filter(nc)%num_urbanl, filter(nc)%urbanl,        &
                filter(nc)%num_urbanc, filter(nc)%urbanc,        &
                filter(nc)%num_urbanp, filter(nc)%urbanp,        &
                atm2lnd_vars, urbanparams_vars, soilstate_vars,  &
                frictionvel_vars, energyflux_vars, &
                nstep_mod, dtime_mod, year_curr, mon_curr, day_curr, secs_curr)
       ! Fluxes for all lake landunits
       call LakeFluxes(bounds_clump,                                         &
                   filter(nc)%num_lakec, filter(nc)%lakec,                          &
                   filter(nc)%num_lakep, filter(nc)%lakep,                          &
                   atm2lnd_vars, solarabs_vars, frictionvel_vars, &
                   energyflux_vars,  lakestate_vars)

       ! ============================================================================
       ! DUST and VOC emissions
       ! ============================================================================

       ! Dust mobilization (C. Zender's modified codes)
       call DustEmission(bounds_clump,                        &
              filter(nc)%num_nolakep, filter(nc)%nolakep,     &
              atm2lnd_vars, soilstate_vars, canopystate_vars, &
              frictionvel_vars, dust_vars)
       ! Dust dry deposition (C. Zender's modified codes)
       call DustDryDep(bounds_clump,atm2lnd_vars, frictionvel_vars, dust_vars)

       ! ============================================================================
       ! Determine temperatures
       ! ============================================================================

       ! Set lake temperature
       call LakeTemperature(bounds_clump,                  &
                  filter(nc)%num_lakec, filter(nc)%lakec,  &
                  filter(nc)%num_lakep, filter(nc)%lakep,  &
                  solarabs_vars, soilstate_vars, ch4_vars, &
                  energyflux_vars, lakestate_vars, dtime_mod)
       !call t_stopf('bgplake')
       ! Set soil/snow temperatures including ground temperature
       call SoilTemperature(bounds_clump,                            &
                   filter(nc)%num_urbanl  , filter(nc)%urbanl,       &
                   filter(nc)%num_nolakec , filter(nc)%nolakec,      &
                   atm2lnd_vars, urbanparams_vars, canopystate_vars, &
                   solarabs_vars, soilstate_vars, energyflux_vars, dtime_mod)

       ! ============================================================================
       ! update surface fluxes for new ground temperature.
       ! ============================================================================
       call SoilFluxes(bounds_clump,                                &
                   filter(nc)%num_urbanl,  filter(nc)%urbanl,       &
                   filter(nc)%num_nolakec, filter(nc)%nolakec,      &
                   filter(nc)%num_nolakep, filter(nc)%nolakep,      &
                   atm2lnd_vars, solarabs_vars,  canopystate_vars,  &
                   energyflux_vars, dtime_mod)
       ! ============================================================================
       ! Perform averaging from patch level to column level
       ! ============================================================================

       call clm_drv_patch2col(bounds_clump, filter(nc)%num_nolakec, filter(nc)%nolakec, &
            energyflux_vars)
       ! ============================================================================
       ! Vertical (column) soil and surface hydrology
       ! ============================================================================

       ! Note that filter_snowc and filter_nosnowc are returned by
       ! LakeHydrology after the new snow filter is built

       call HydrologyNoDrainage(bounds_clump,                     &
            filter(nc)%num_nolakec, filter(nc)%nolakec,           &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc,     &
            filter(nc)%num_hydrononsoic, filter(nc)%hydrononsoic, &
            filter(nc)%num_urbanc, filter(nc)%urbanc,             &
            filter(nc)%num_snowc, filter(nc)%snowc,               &
            filter(nc)%num_nosnowc, filter(nc)%nosnowc,           &
            canopystate_vars,atm2lnd_vars, soilstate_vars, energyflux_vars,        &
            soilhydrology_vars, aerosol_vars, dtime_mod)

       !  Calculate column-integrated aerosol masses, and
       !  mass concentrations for radiative calculations and output
       !  (based on new snow level state, after SnowFilter is rebuilt.
       !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there
       !  can be zero snow layers but an active column in filter)

       call AerosolMasses( bounds=bounds_clump,                                   &
                   num_on=filter(nc)%num_snowc, filter_on=filter(nc)%snowc,       &
                   num_off=filter(nc)%num_nosnowc, filter_off=filter(nc)%nosnowc, &
                   aerosol_vars=aerosol_vars, dtime=dtime_mod)

       ! ============================================================================
       ! Lake hydrology
       ! ============================================================================
       ! Note that filter_lakesnowc and filter_lakenosnowc are returned by
       ! LakeHydrology after the new snow filter is built
      call LakeHydrology(bounds_clump,                                  &
                  filter(nc)%num_lakec, filter(nc)%lakec,                      &
                  filter(nc)%num_lakep, filter(nc)%lakep,                      &
                  filter(nc)%num_lakesnowc, filter(nc)%lakesnowc,              &
                  filter(nc)%num_lakenosnowc, filter(nc)%lakenosnowc,          &
                  atm2lnd_vars,  soilstate_vars,  &
                  energyflux_vars, aerosol_vars, lakestate_vars, dtime_mod)


       !  Calculate column-integrated aerosol masses, and
       !  mass concentrations for radiative calculations and output
       !  (based on new snow level state, after SnowFilter is rebuilt.
       !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there
       !  can be zero snow layers but an active column in filter)

       call AerosolMasses(bounds=bounds_clump,                                   &
                   num_on=filter(nc)%num_lakesnowc,    filter_on=filter(nc)%lakesnowc,   &
                   num_off=filter(nc)%num_lakenosnowc, filter_off=filter(nc)%lakenosnowc, &
                   aerosol_vars=aerosol_vars,  dtime=dtime_mod)

           call SnowAge_grain(bounds_clump,                         &
                   filter(nc)%num_lakesnowc, filter(nc)%lakesnowc,     &
                   filter(nc)%num_lakenosnowc, filter(nc)%lakenosnowc, dtime_mod)

       call set_fracsno(bounds_clump)

       ! ============================================================================
       ! Snow aging routine based on Flanner and Zender (2006), Linking snowpack
       ! microphysics and albedo evolution, JGR, and Brun (1989), Investigation of
       ! wet-snow metamorphism in respect of liquid-water content, Ann. Glaciol.
       ! ============================================================================
       ! Note the snow filters here do not include lakes
       ! TODO: move this up
       call SnowAge_grain(bounds_clump,                         &
                   filter(nc)%num_snowc, filter(nc)%snowc,     &
                   filter(nc)%num_nosnowc, filter(nc)%nosnowc, dtime_mod)
       ! ============================================================================
       ! Update sediment fluxes from land unit
       ! ============================================================================

       if (use_erosion) then
          call SoilErosion(bounds_clump, filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
               atm2lnd_vars, canopystate_vars, soilstate_vars, sedflux_vars, dtime_mod)
       end if
       ! ============================================================================
       ! Ecosystem dynamics: Uses CN, or static parameterizations
       ! ============================================================================

       if ((mon_curr == 1 .and. day_curr == 1 .and. secs_curr == 0) .and. nstep_mod > 1) then
          call crop_vars_CropIncrementYear(filter(nc)%num_pcropp,filter(nc)%pcropp,crop_vars)
       end if

         ! FIX(SPM,032414)  push these checks into the routines below and/or make this consistent.
       if (.not. use_fates) then
         if( .not. is_active_betr_bgc) then
           if (use_cn) then

             ! fully prognostic canopy structure and C-N biogeochemistry
             ! - crop model:  crop algorithms called from within CNEcosystemDyn

             !===========================================================================================
             ! clm_interface: 'EcosystemDynNoLeaching' is divided into 2 subroutines (1 & 2): BEGIN
             ! EcosystemDynNoLeaching1 is called before clm_interface
             ! EcosystemDynNoLeaching2 is called after clm_interface
             !===========================================================================================
             call EcosystemDynNoLeaching1(bounds_clump,                               &
                                    filter(nc)%num_soilc, filter(nc)%soilc,                          &
                                    filter(nc)%num_soilp, filter(nc)%soilp,                          &
                                    cnstate_vars, atm2lnd_vars,                    &
                                    canopystate_vars, soilstate_vars, crop_vars,   &
                                    ch4_vars, photosyns_vars,                                        &
                                    dtime_mod, dayspyr_mod,year_curr, mon_curr, day_curr, secs_curr)
             !--------------------------------------------------------------------------------

             !--------------------------------------------------------------------------------

             call EcosystemDynNoLeaching2(bounds_clump,                                   &
                                 filter(nc)%num_soilc, filter(nc)%soilc,                                  &
                                 filter(nc)%num_soilp, filter(nc)%soilp,                                  &
                                 filter(nc)%num_pcropp, filter(nc)%pcropp, doalb,                         &
                                 cnstate_vars, atm2lnd_vars,canopystate_vars,&
                                 soilstate_vars, crop_vars, ch4_vars, &
                                 photosyns_vars, soilhydrology_vars, energyflux_vars, sedflux_vars, &
                                 year_curr, mon_curr, day_curr, secs_curr, secs_curr, 0, dayspyr_mod, dtime_mod, nstep_mod)

             ! clm_interface: 'EcosystemDynNoLeaching' is divided into 2 subroutines (1 & 2): END
             !===========================================================================================

             call AnnualUpdate(bounds_clump,            &
                  filter(nc)%num_soilc, filter(nc)%soilc, &
                  filter(nc)%num_soilp, filter(nc)%soilp, &
                  cnstate_vars, dtime_mod)
           else ! not use_cn

             if (doalb) then
                ! Prescribed biogeography - prescribed canopy structure, some prognostic carbon fluxes
                !!NOTE: what is use_lai_streams?? requires mpicom for shr_strdata_advance !!
                call SatellitePhenology(bounds_clump,               &
                     filter(nc)%num_nolakep, filter(nc)%nolakep,    &
                     canopystate_vars)
             end if

          end if  ! end of if-use_cn
          end if  ! end of is_active_betr_bgc
        end if    ! end of if-use_fates



         ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)
         call depvel_compute(bounds_clump, &
              atm2lnd_vars, canopystate_vars, frictionvel_vars, &
              photosyns_vars, drydepvel_vars)


         if (use_lch4 .and. .not. is_active_betr_bgc) then
           !warning: do not call ch4 before AnnualUpdate, which will fail the ch4 model
           call CH4 (bounds_clump,                                                                  &
               filter(nc)%num_soilc, filter(nc)%soilc,                                             &
               filter(nc)%num_lakec, filter(nc)%lakec,                                             &
               filter(nc)%num_soilp, filter(nc)%soilp,                                             &
               atm2lnd_vars, lakestate_vars, canopystate_vars, soilstate_vars, soilhydrology_vars, &
               energyflux_vars, ch4_vars, lnd2atm_vars, dtime_mod )
         end if

       ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)
       call depvel_compute(bounds_clump, &
            atm2lnd_vars, canopystate_vars, frictionvel_vars, &
            photosyns_vars, drydepvel_vars)
       ! ============================================================================
       ! Calculate soil/snow hydrology with drainage (subsurface runoff)
       ! ============================================================================

         call HydrologyDrainage(bounds_clump,                 &
            filter(nc)%num_nolakec, filter(nc)%nolakec,       &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
            filter(nc)%num_urbanc, filter(nc)%urbanc,         &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,     &
            atm2lnd_vars, glc2lnd_vars,     &
            soilhydrology_vars, soilstate_vars, dtime_mod)

       if (.not. use_fates) then

          if (use_cn) then

            if (.not. is_active_betr_bgc)then
             ! FIX(SPM,032414) there are use_fates checks in this routine...be consistent
             ! (see comment above re: no leaching
             call EcosystemDynLeaching( bounds_clump,                &
                        filter(nc)%num_soilc, filter(nc)%soilc,      &
                        filter(nc)%num_soilp, filter(nc)%soilp,      &
                        filter(nc)%num_pcropp,filter(nc)%pcropp, doalb,  &
                        cnstate_vars, frictionvel_vars, canopystate_vars,&
                        dtime_mod )
             end if

             if (doalb) then
                call VegStructUpdate(filter(nc)%num_soilp, filter(nc)%soilp,   &
                     frictionvel_vars, cnstate_vars, &
                     canopystate_vars, crop_vars, dtime_mod)
             end if

          end if
       end if

       call ColWaterBalanceCheck(bounds_clump, &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c, &
            atm2lnd_vars, glc2lnd_vars, solarabs_vars, &
            energyflux_vars, canopystate_vars,surfalb_vars, dtime_mod, nstep_mod)

       call GridBalanceCheck(bounds_clump                , &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c , &
            atm2lnd_vars, glc2lnd_vars, solarabs_vars, &
            energyflux_vars, canopystate_vars        , &
            soilhydrology_vars,surfalb_vars, dtime_mod)

       call WaterBudget_SetEndingMonthlyStates(bounds_clump, nstep_mod, year_curr, mon_curr, day_curr, secs_curr)

       if (.not. use_fates)then
          if (use_cn) then
             if (nstep_mod < 2 )then
                if (masterproc) then
                   write(iulog,*) '--WARNING-- skipping CN balance check for first timestep'
                end if
             else

                call ColCBalanceCheck(bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     col_cs, dtime_mod)

                call ColNBalanceCheck(bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     dtime_mod, year_curr, mon_curr, day_curr, secs_curr)

                call ColPBalanceCheck(bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     dtime_mod, year_curr, mon_curr, day_curr, secs_curr)

             end if
          end if
       end if

       ! ============================================================================
       ! Determine albedos for next time step
       ! ============================================================================

       if (doalb) then
          ! Albedos for non-urban columns
          call SurfaceAlbedo(bounds_clump,                      &
               filter_inactive_and_active(nc)%num_nourbanc,     &
               filter_inactive_and_active(nc)%nourbanc,         &
               filter_inactive_and_active(nc)%num_nourbanp,     &
               filter_inactive_and_active(nc)%nourbanp,         &
               filter_inactive_and_active(nc)%num_urbanc,       &
               filter_inactive_and_active(nc)%urbanc,           &
               filter_inactive_and_active(nc)%num_urbanp,       &
               filter_inactive_and_active(nc)%urbanp,           &
               nextsw_cday_mod, declinp1,                           &
               aerosol_vars, canopystate_vars, &
               lakestate_vars,  surfalb_vars )
          ! Albedos for urban columns
          if (filter_inactive_and_active(nc)%num_urbanl > 0) then
             call UrbanAlbedo(bounds_clump,                  &
                  filter_inactive_and_active(nc)%num_urbanl, &
                  filter_inactive_and_active(nc)%urbanl,     &
                  filter_inactive_and_active(nc)%num_urbanc, &
                  filter_inactive_and_active(nc)%urbanc,     &
                  filter_inactive_and_active(nc)%num_urbanp, &
                  filter_inactive_and_active(nc)%urbanp,     &
                  urbanparams_vars, solarabs_vars, surfalb_vars)
          end if

       end if
      ! ============================================================================
      ! Determine gridcell averaged properties to send to atm
      ! ===========================================================================
       call lnd2atm(bounds_clump,                          &
                atm2lnd_vars, surfalb_vars,  frictionvel_vars, &
               energyflux_vars, solarabs_vars, drydepvel_vars,&
               vocemis_vars, dust_vars, ch4_vars, soilhydrology_vars, lnd2atm_vars)

        ! ============================================================================
        ! Determine gridcell averaged properties to send to glc
        ! ============================================================================

        if (create_glacier_mec_landunit) then
             call update_lnd2glc_GPU(lnd2glc_vars, bounds_clump,       &
                  filter(nc)%num_do_smb_c, filter(nc)%do_smb_c, .false.)
        end if
        ! ============================================================================
        ! Update accumulators
        ! ============================================================================

        ! FIX(SPM,032414) double check why this isn't called for ED

        if (nstep_mod > 0) then

              call atm2lnd_UpdateAccVars(atm2lnd_vars, bounds_clump, nstep_mod)
              !call top_as%UpdateAccVars(bounds_proc)  !!only needed if use_fates
              call update_acc_vars_top_afGPU (top_af, bounds_clump, nstep_mod)

              call update_acc_vars_veg_es_GPU(veg_es, bounds_clump, int(dtime_mod), nstep_mod, &
                        end_cd_mod, mon_curr,day_curr, secs_curr)

              call CanopyState_UpdateAccVars (canopystate_vars, bounds_clump, int(dtime_mod),nstep_mod)

              if (crop_prog) then
                call crop_vars_UpdateAccVars(crop_vars, bounds_clump, int(dtime_mod), nstep_mod)
              end if

              call CNState_UpdateAccVars (cnstate_vars, bounds_clump, nstep_mod)
        end if

    end do

    !$acc end parallel
    ! ============================================================================
    ! Write global average diagnostics to standard output
    ! ============================================================================
    

    !!if (wrtdia) call mpi_barrier( mpicom,ier)
    !!call t_startf('wrtdiag')
    !!call write_diagnostic(bounds_proc, wrtdia, nstep_mod, lnd2atm_vars)
    !!call t_stopf('wrtdiag')

    ! ============================================================================
    ! Update history buffer
    ! ============================================================================
    
    !if(step_count == 24) then
        call hist_update_hbuf_gpu(step_count,24*numdays, nclumps)
    !    stop
    !end if
    ! ============================================================================
    ! Compute water budget
    ! ============================================================================
    if (get_nstep()>0 .and. do_budgets) then
       call WaterBudget_Run(bounds_clump, atm2lnd_vars, lnd2atm_vars,  &
            soilhydrology_vars)
       call WaterBudget_Accum()
       call WaterBudget_Print(budget_inst,  budget_daily,  budget_month,  &
            budget_ann,  budget_ltann,  budget_ltend,year_curr, mon_curr, day_curr, secs_curr)
    endif

    ! ============================================================================
    ! History/Restart output
    ! ============================================================================
    !! MOVE outside of clm_driver ?
    if (.not. use_noio) then
        
       call t_startf('clm_drv_io')
       ! Create history and write history tapes if appropriate
       call t_startf('clm_drv_io_htapes')

       if(step_count == 2400)  rstwr = .true. 
         
       call hist_htapes_wrapup(step_count,24*numdays, rstwr, nlend, bounds_proc,                    &
            soilstate_vars%watsat_col(bounds_proc%begc:bounds_proc%endc, 1:), &
            soilstate_vars%sucsat_col(bounds_proc%begc:bounds_proc%endc, 1:), &
            soilstate_vars%bsw_col(bounds_proc%begc:bounds_proc%endc, 1:),    &
            soilstate_vars%hksat_col(bounds_proc%begc:bounds_proc%endc, 1:))
        call set_gpu_tape()
        call t_stopf('clm_drv_io_htapes')
        
       ! Write restart/initial files if appropriate
       if (rstwr) then
          call t_startf('clm_drv_io_wrest')
          filer = restFile_filename(rdate=rdate)
          print *, "Copying out data for restart"
         !$acc exit data copyout(atm2lnd_vars, canopystate_vars, energyflux_vars,&
         !$acc col_ef, veg_ef, frictionvel_vars, lakestate_vars, photosyns_vars, &
         !$acc  soilhydrology_vars, soilstate_vars, solarabs_vars, grc_wf,col_wf,&
         !$acc  veg_wf, lun_es, col_es,ch4_vars, veg_es,clmptr_rs,clmptr_ra)
         !$acc exit data copyout(lun_ws, col_ws, veg_ws, aerosol_vars,surfalb_vars,&
         !$acc cnstate_vars, col_cs, veg_cs, c13_col_cs, c13_veg_cs,c14_col_cs,&
         !$acc c14_veg_cs, c14_veg_cs, col_cf, veg_cf, col_ns, veg_ns, &
         !$acc col_nf,veg_nf,col_ps,veg_ps, col_pf,veg_pf, crop_vars ) 
         
          print *, "calling restFile_write:" 
         call restFile_write( bounds_proc, filer,                                            &
               atm2lnd_vars, aerosol_vars, canopystate_vars, cnstate_vars,                    &
               carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, carbonflux_vars, &
               ch4_vars, energyflux_vars, frictionvel_vars, lakestate_vars,        &
               nitrogenstate_vars, nitrogenflux_vars, photosyns_vars, soilhydrology_vars,     &
               soilstate_vars, solarabs_vars, surfalb_vars, temperature_vars,                 &
               waterflux_vars, waterstate_vars, sedflux_vars,                                 &
               phosphorusstate_vars,phosphorusflux_vars,                                      &
               ep_betr, alm_fates, crop_vars, rdate=rdate )

         !----------------------------------------------
         ! pflotran (off now)
         ! if (use_pflotran) then
         !     call clm_pf_write_restart(rdate)
         ! end if
         !----------------------------------------------


          call t_stopf('clm_drv_io_wrest')
       end if
       call t_stopf('clm_drv_io')

    end if

    if (use_pflotran .and. nstep_mod>=nestep) then
       call clm_pf_finalize()
    end if
    step_count = step_count + 1
     
  end subroutine clm_drv

  !-----------------------------------------------------------------------
    subroutine clm_drv_init(bounds, &
       num_nolakec, filter_nolakec, &
       num_nolakep, filter_nolakep, &
       num_soilp  , filter_soilp, &
       canopystate_vars, photosyns_vars,&
       col_wf, col_ef)
    !
    !$acc routine seq
    ! !DESCRIPTION:
    ! Initialization of clm driver variables needed from previous timestep
    !
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
  !  use shr_infnan_mod     , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar         , only : nlevsno
    use clm_varcon         , only : h2osno_max
    use landunit_varcon    , only : istice_mec
    use CanopyStateType    , only : canopystate_type
    use ColumnDataType     , only : column_energy_flux, column_water_flux
    use PhotosynthesisType , only : photosyns_type
  !  use VegetationType,      only : veg_pp
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds
    integer               , intent(in)    :: num_nolakec       ! number of non-lake points in column filter
    integer               , intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
    integer               , intent(in)    :: num_nolakep       ! number of non-lake points in patch filter
    integer               , intent(in)    :: filter_nolakep(:) ! patch filter for non-lake points
    integer               , intent(in)    :: num_soilp         ! number of soil points in patch filter
    integer               , intent(in)    :: filter_soilp(:)   ! patch filter for soil points

    type(canopystate_type)  , intent(inout) :: canopystate_vars
    type(photosyns_type)    , target, intent(inout) :: photosyns_vars
    type(column_water_flux) , target, intent(inout) :: col_wf
    type(column_energy_flux), target, intent(inout) :: col_ef
    !
    ! !LOCAL VARIABLES:
    integer :: l, c, p, f, j         ! indices
    integer :: fp, fc                  ! filter indices
    !-----------------------------------------------------------------------

      ! Initialize intracellular CO2 (Pa) parameters each timestep for use in VOCEmission
      do p = bounds%begp,bounds%endp
         photosyns_vars%cisun_z_patch(p,:) = -999._r8
         photosyns_vars%cisha_z_patch(p,:) = -999._r8

      end do


      do c = bounds%begc,bounds%endc
         l = col_pp%landunit(c)

         ! Save snow mass at previous time step
         col_ws%h2osno_old(c) = col_ws%h2osno(c)

         ! Decide whether to cap snow
         if (col_ws%h2osno(c) > h2osno_max) then
            col_ws%do_capsnow(c) = .true.
         else
            col_ws%do_capsnow(c) = .false.
         end if

         ! Reset flux from beneath soil/ice column
         col_ef%eflx_bot(c)  = 0._r8

         ! Initialize col_wf%qflx_glcice everywhere, to zero.
         col_wf%qflx_glcice(c) = 0._r8

      end do

     ! ! Initialize fraction of vegetation not covered by snow
      do p = bounds%begp, bounds%endp
         if (veg_pp%active(p)) then
           canopystate_vars%frac_veg_nosno_patch(p) = canopystate_vars%frac_veg_nosno_alb_patch(p)
         else
           canopystate_vars%frac_veg_nosno_patch(p) = 0
         end if
      end do

     ! ! Initialize set of previous time-step variables
     ! ! Ice fraction of snow at previous time step
      do j = -nlevsno+1,0
         do f = 1, num_nolakec
            c = filter_nolakec(f)
            if (j >= col_pp%snl(c) + 1) then
               col_ws%frac_iceold(c,j) = col_ws%h2osoi_ice(c,j)/(col_ws%h2osoi_liq(c,j)+col_ws%h2osoi_ice(c,j))
            end if
         end do
      end do

  end subroutine clm_drv_init

  !-----------------------------------------------------------------------
  subroutine clm_drv_patch2col (bounds, num_nolakec, filter_nolakec, &
       energyflux_vars)
    !
    ! !DESCRIPTION:
    ! Averages over all patchs for variables defined over both soil and lake
    ! to provide the column-level averages of state and flux variables
    ! defined at the patch level.
    !
    ! !USES:
      !$acc routine seq
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

    call p2c (bounds, num_allc, filter_allc, &
         qflx_snwcp_liq_patch(bounds%begp:bounds%endp), &
         qflx_snwcp_liq_col  (bounds%begc:bounds%endc))

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
         qflx_snwcp_ice_patch(bounds%begp:bounds%endp), &
         qflx_snwcp_ice_col  (bounds%begc:bounds%endc))

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

  end subroutine clm_drv_patch2col


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

end module clm_driver
