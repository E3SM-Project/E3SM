module elm_driver
   
   !-----------------------------------------------------------------------
   ! !DESCRIPTION:
   ! This module provides the main ELM driver physics calling sequence.  Most
   ! computations occurs over ``clumps'' of gridcells (and associated subgrid
   ! scale entities) assigned to each MPI process. Computation is further
   ! parallelized by looping over clumps on each process using shared memory OpenMP.
   !
   ! !USES:
   use write_filterMod, only: write_filter 
   use shr_kind_mod           , only : r8 => shr_kind_r8
   use shr_sys_mod            , only : shr_sys_flush
   use shr_log_mod            , only : errMsg => shr_log_errMsg
   use elm_varpar             , only : nlevtrc_soil, nlevsoi
   use elm_varpar             , only : nlevsno, nlevgrnd, crop_prog
   use elm_varctl             , only : wrtdia, iulog, create_glacier_mec_landunit, use_fates, use_betr, use_extrasnowlayers
   use elm_varctl             , only : use_cn, use_lch4, use_voc, use_noio, use_c13, use_c14
   use elm_varctl             , only : use_erosion
   use clm_time_manager       , only : get_step_size, get_curr_date, get_ref_date, get_nstep, is_beg_curr_day, get_curr_time_string
   use clm_time_manager       , only : get_curr_calday, get_days_per_year, get_prev_date
   use spmdMod                , only : masterproc, mpicom, iam
   use decompMod              , only : get_proc_clumps, get_clump_bounds, get_proc_bounds, bounds_type
   use filterMod              , only : filter, filter_inactive_and_active, proc_filter
   use filterMod              , only : setProcFilters, createProcessorFilter, proc_filter_inactive_and_active
   use filterMod              , only : updateFracNoSnoFilters 
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
   !use EcosystemDynMod      , only : EcosystemDynNoLeaching1, EcosystemDynNoLeaching2
   
   use EcosystemDynMod      , only : EcosystemDynLeaching
   use VegStructUpdateMod   , only : VegStructUpdate
   use AnnualUpdateMod      , only : AnnualUpdate
   use EcosystemBalanceCheckMod      , only : BeginColCBalance, BeginColNBalance, ColCBalanceCheck, ColNBalanceCheck
   use EcosystemBalanceCheckMod      , only : BeginColPBalance, ColPBalanceCheck
   use EcosystemBalanceCheckMod      , only : BeginGridCBalance,GridCBalanceCheck
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
   use DaylengthMod           , only : UpdateDaylength, first_step
   use perf_mod
   !
   use elm_instMod            , only : ch4_vars, ep_betr
   use elm_instMod            , only : carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars
   use elm_instMod            , only : carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars
   use elm_instMod            , only : nitrogenstate_vars,   nitrogenflux_vars
   use elm_instMod            , only : phosphorusstate_vars, phosphorusflux_vars
   use elm_instMod            , only : crop_vars, cnstate_vars
   use elm_instMod            , only : dust_vars, vocemis_vars
   use elm_instMod            , only : drydepvel_vars,aerosol_vars
   use elm_instMod            , only : canopystate_vars,energyflux_vars
   use elm_instMod            , only : frictionvel_vars,lakestate_vars
   use elm_instMod            , only : photosyns_vars,sedflux_vars
   use elm_instMod            , only : soilstate_vars,soilhydrology_vars
   use elm_instMod            , only : solarabs_vars,soilhydrology_vars
   use elm_instMod            , only : surfalb_vars,surfrad_vars
   use elm_instMod            , only : temperature_vars
   use elm_instMod            , only : waterflux_vars
   use elm_instMod            , only : waterstate_vars
   use elm_instMod            , only : atm2lnd_vars,lnd2atm_vars
   use elm_instMod            , only : glc2lnd_vars,lnd2glc_vars
   use elm_instMod            , only : soil_water_retention_curve
   use elm_instMod            , only : chemstate_vars
   use elm_instMod            , only : alm_fates
   use elm_instMod            , only : PlantMicKinetics_vars
   use tracer_varcon          , only : is_active_betr_bgc
   use CNEcosystemDynBetrMod  , only : CNEcosystemDynBetr, CNFluxStateBetrSummary
   use UrbanParamsType        , only : urbanparams_vars
   
   use elm_instMod            , only : col_es ! Why is this being used via elm_instMod??
   
   use GridcellType           , only : grc_pp
   use GridcellDataType       , only : grc_cs, c13_grc_cs, c14_grc_cs
   use GridcellDataType       , only : grc_cf, c13_grc_cf, c14_grc_cf
   use GridcellDataType       , only : grc_nf, grc_pf, grc_ef, grc_wf, grc_ws
   use GridcellDataType       , only : grc_es, grc_ns, grc_ps
   use TopounitDataType         , only : top_as, top_af
   use TopounitType             , only : top_pp
   use LandunitType           , only : lun_pp
   use ColumnType             , only : col_pp
   use ColumnDataType         , only : col_es, col_ef, col_ws, col_wf
   use ColumnDataType         , only : col_cs, c13_col_cs, c14_col_cs
   use ColumnDataType         , only : col_cf, c13_col_cf, c14_col_cf
   use ColumnDataType         , only : col_ns, col_nf
   use ColumnDataType         , only : col_ps, col_pf
   use VegetationType         , only : veg_pp
   use VegetationDataType       , only : veg_es, veg_ws, veg_wf
   use VegetationDataType       , only : veg_cs, c13_veg_cs, c14_veg_cs
   use VegetationDataType       , only : c13_veg_cf, c14_veg_cf
   use VegetationDataType       , only : veg_ns, veg_nf
   use VegetationDataType       , only : veg_ps, veg_pf
   use VegetationDataType       , only : veg_cf, veg_ef
   use VegetationPropertiesType , only : veg_vp
   
   !----------------------------------------------------------------------------
   ! bgc interface & pflotran:
   use elm_varctl             , only : use_elm_interface
   use elm_instMod            , only : elm_interface_data
   use elm_interface_funcsMod , only : get_elm_data
   ! (1) clm_bgc through interface
   use elm_varctl             , only : use_elm_bgc
   use elm_interface_funcsMod , only : elm_bgc_run, update_bgc_data_elm2elm
   ! (2) pflotran
   use clm_time_manager            , only : nsstep, nestep
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
   use elm_varctl
   use PhotosynthesisMod        , only : params_inst
   use SharedParamsMod          , only : ParamsShareInst
   use CH4Mod                   , only : CH4ParamsInst
   use CNDecompCascadeConType   , only : decomp_cascade_con
   use DecompCascadeBGCMod      , only : DecompBGCParamsInst
   use GapMortalityMod          , only : cngapmortparamsinst
   use DecompCascadeCNMod       , only : DecompCNParamsInst
   use NitrifDenitrifMod        , only : NitrifDenitrifParamsInst
   use SoilLittDecompMod        , only : cndecompparamsinst
   use AllocationMod            , only : AllocParamsInst, nu_com_phosphatase,nu_com_nfix
   use LandunitDataType         , only : lun_ef, lun_es, lun_ws
   
   use NitrogenDynamicsMod, only : CNNDynamicsParamsInst
   use dynSubgridControlMod , only : dyn_subgrid_control_inst
   use dynInitColumnsMod    , only : initialize_new_columns
   use dynConsBiogeophysMod , only : dyn_hwcontent_init, dyn_hwcontent_final
   use dynConsBiogeochemMod , only : dyn_cnbal_patch, dyn_cnbal_column
   use reweightMod          , only : reweight_wrapup
   use subgridWeightsMod    , only : compute_higher_order_weights, set_subgrid_diagnostic_fields
   use subgridWeightsMod    , only : subgrid_weights_diagnostics
   use CarbonStateUpdate1Mod   , only : CarbonStateUpdateDynPatch
   use NitrogenStateUpdate1Mod , only : NitrogenStateUpdateDynPatch
   use PhosphorusStateUpdate1Mod     , only : PhosphorusStateUpdateDynPatch
   
   use dynSubgridDriverMod    , only : dynSubgrid_driver,dynSubgrid_wrapup_weight_changes
   use dynSubgridDriverMod    , only : prior_weights  
   !use dynColumnStateUpdaterMod , only: column_state_updater
   use elm_instMod , only : patch_state_updater
   use elm_instMod , only : column_state_updater
   use domainMod , only : ldomain_gpu
   use histGPUMod , only : tape_gpu, htape_gpu_init, hist_update_hbuf_gpu
   use histFileMod , only : elmptr_ra, elmptr_rs, tape
   use update_accMod
   use timeinfoMod
   use ColumnDataType, only : col_cs_summary, col_ns_summary,col_ps_summary
   use VegetationDataType, only : veg_cs_summary, veg_ns_summary, veg_ps_summary
   use pftvarcon
   use decompMod , only : clumps, procinfo
   use ColumnWorkRoutinesMod 
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
      ! First phase of the elm driver calling the elm physics. An outline of
      ! the calling tree is given in the description of this module.
      !
      ! !USES:
      #if _CUDA
      use cudafor
      #endif
      use shr_orb_mod_elm
      use decompMod , only : init_proc_clump_info, gpu_clumps, gpu_procinfo
      use decompMod , only : get_clump_bounds_gpu
      use ColumnDataType, only : col_cf_setvalues, col_nf_SetValues, col_pf_SetValues
      use VegetationDataType, only : veg_cf_summary_rr
      use VegetationSummaryRoutinesMod
      use NitrogenDynamicsMod, only : NitrogenDeposition, NitrogenFixation, NitrogenFixation_balance, NitrogenFert
      use NitrogenDynamicsMod, only : CNSoyFix
      use MaintenanceRespMod , only : MaintenanceResp
      use PhosphorusDynamicsMod,only: PhosphorusWeathering,PhosphorusBiochemMin,PhosphorusBiochemMin_balance,PhosphorusDeposition
      use DecompCascadeCNMod   , only: decomp_rate_constants_cn
      use AllocationMod, only : Allocation1_PlantNPDemand
      use SoilLittDecompMod, only: SoilLittDecompAlloc, SoilLittDecompAlloc2
      use PhenologyMod, only : Phenology, CNLitterToColumn
      use GrowthRespMod, only : GrowthResp
      use CarbonStateUpdate1Mod, only:CarbonStateUpdate0,CarbonStateUpdate_Phase1_col, CarbonStateUpdate_Phase1_pft
      use NitrogenStateUpdate1Mod,only:NitrogenStateUpdate_Phase1_col, NitrogenStateUpdate_Phase1_pft
      use PhosphorusStateUpdate1Mod,only:PhosphorusStateUpdate_Phase1_col, PhosphorusStateUpdate_Phase1_pft
      use CarbonStateUpdate1Mod , only : CarbonStateDynGridUpdate
      use NitrogenStateUpdate1Mod , only : NitrogenStateDynGridUpdate
      use PhosphorusStateUpdate1Mod , only : PhosphorusStateDynGridUpdate
      use SoilLittVertTranspMod,only:SoilLittVertTransp, createLitterTransportList, transport_ptr_list 
      use GapMortalityMod,only:GapMortality
      use NitrogenStateUpdate3Mod   , only: NitrogenStateUpdate3
      use PhosphorusStateUpdate3Mod , only: PhosphorusStateUpdate3
      use CarbonStateUpdate2Mod     , only: CarbonStateUpdate2, CarbonStateUpdate2h
      use NitrogenStateUpdate2Mod   , only: NitrogenStateUpdate2, NitrogenStateUpdate2h
      use PhosphorusStateUpdate2Mod , only: PhosphorusStateUpdate2, PhosphorusStateUpdate2h
      use FireMod              , only: FireArea, FireFluxes
      use C14DecayMod          , only: C14Decay, C14BombSpike
      use WoodProductsMod      , only: WoodProducts
      use CropHarvestPoolsMod  , only: CropHarvestPools
      !
      use dynSubgridControlMod, only : get_flanduse_timeseries
      use dynSubgridControlMod, only : get_do_transient_pfts, get_do_transient_crops
      use dynSubgridControlMod, only : get_do_harvest
      use dynPriorWeightsMod  , only : set_prior_weights
      use dynUpdateModAcc 
      use glc2lndMod  , only : glc2lnd_vars_update_glc2lnd_acc
      use ColumnWorkRoutinesMod
      use dynpftFileMod, only : dynpft_interp
      use dynHarvestMod, only : dynHarvest_interp
      use dyncropFileMod      , only : dyncrop_init, dyncrop_interp
      use writeMod , only : write_vars 
      use histFileMod, only : ntapes  
      use dynSubgridAdjustmentsMod, only : dyn_col_cs_Adjustments,dyn_col_ns_Adjustments,dyn_col_ps_Adjustments
      use verificationMod
      use ForcingUpdateMod , only : update_forcings_cplbypass
      use elm_instMod , only : cpl_bypass_input   
      use elm_varpar , only : ndecomp_pools 
      use UrbanParamsType, only : urban_hac_int 
      #ifdef _OPENACC 
      #define gpuflag 1
      #else
      #define gpuflag 0
      #endif

      !
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
      integer            :: nstep                   ! time step number
      real(r8)           :: dtime                   ! land model time step (sec)
      integer            :: nc, c, p, l, g, fp, j,t ! indices
      integer            :: nclumps                 ! number of clumps on this processor
      integer            :: yrp1                    ! year (0, ...) for nstep+1
      integer            :: monp1                   ! month (1, ..., 12) for nstep+1
      integer            :: dayp1                   ! day of month (1, ..., 31) for nstep+1
      integer            :: secp1                   ! seconds into current date for nstep+1
      integer            :: yr                      ! year (0, ...)
      integer            :: mon                     ! month (1, ..., 12)
      integer            :: day                     ! day of month (1, ..., 31)
      integer            :: sec                     ! seconds of the day
      integer            :: ncdate                  ! current date
      integer            :: nbdate                  ! base date (reference date)
      integer            :: kyr                     ! thousand years, equals 2 at end of first year
      character(len=256) :: filer                   ! restart file name
      integer            :: ier                     ! error code
      character(len=256) :: dateTimeString
      type(bounds_type)  :: bounds_clump
      type(bounds_type)  :: bounds_proc
      integer :: fc,fl
      integer :: num_soilc,num_nourbanl 
      real :: startt, stopt, eco_startt, eco_stopt,outer_start, outer_stop 
      #if _CUDA
      integer(kind=cuda_count_kind) :: heapsize,free1,free2,total
      integer  :: istat, val
      #endif
      integer :: maxcols, begg,endg, i,maxpfts,begc,begp, endc, endp 
      integer, allocatable :: ncols(:) 
      integer, allocatable :: npfts(:) 
      character(len=8) :: desc
      logical :: transfer_tapes  
      !-----------------------------------------------------------------------
      call get_curr_time_string(dateTimeString)
      if (masterproc) then
         write(iulog,*)'Beginning timestep   : ',trim(dateTimeString)
         write(iulog,*) 'doalb :', doalb 
         call shr_sys_flush(iulog)
      endif
      
      call get_proc_bounds(bounds_proc)
      nclumps = get_proc_clumps()
      !Update time variables
      nstep_mod = get_nstep()
      dtime_mod = real(get_step_size(),r8)
      call get_curr_date(year_curr,mon_curr,day_curr,secs_curr)
      call get_prev_date(year_prev,mon_prev,day_prev,secs_prev)
      write(iulog,*) iam, "STEP :", nstep_mod 
      print *, year_curr,mon_curr, day_curr, secs_curr 
      print *, year_prev, mon_prev, day_prev, secs_prev 
      dayspyr_mod = get_days_per_year()
      jday_mod = get_curr_calday()
      
      if (do_budgets) then
         call WaterBudget_Reset()
         
         if (use_cn) then
            call CNPBudget_Reset()
         end if
      end if
      ! Determine processor bounds and clumps for this processor
      write(iulog,*) "TOTAL Gridcells/clumps:",bounds_proc%endg,nclumps
      if(nstep_mod == 0 ) then
        #if _CUDA
        istat = cudaDeviceGetLimit(heapsize, cudaLimitMallocHeapSize)
        write(iulog,*) "SETTING Heap Limit from", heapsize
        heapsize = 2000_8*1024_8*1024_8
        write(iulog,*) "TO:",heapsize
        istat = cudaDeviceSetLimit(cudaLimitMallocHeapSize,heapsize)
        istat = cudaMemGetInfo(free1, total)
        write(iulog,*) iam,"Free1:",free1
        #endif 
         write(iulog,*) "transferring data to GPU"
         call init_proc_clump_info()
         call createLitterTransportList() 

         begg = bounds_proc%begg; endg = bounds_proc%endg 
         begc = bounds_proc%begc; endc = bounds_proc%endc 
         begp = bounds_proc%begp; endp = bounds_proc%endp 
         
         allocate(ncols(begg:endg)); ncols(:) = 0 
          do c = bounds_proc%begc, bounds_proc%endc
              g = col_pp%gridcell(c) 
              ncols(g) = ncols(g) + 1  
          end do 
          maxcols = maxval(ncols)
          ncols(:) = 1 

          allocate(grc_pp%cols(begg:endg,maxcols))  
          do c = begc, endc 
              g = col_pp%gridcell(c) 
              i = ncols(g) 
              grc_pp%cols(g,i) = c 
              grc_pp%ncolumns(g) = i 
              ncols(g) = ncols(g) + 1
          end do 
           
          deallocate(ncols(:))  

          allocate(npfts(begg:endg)); 
          npfts(:) = 0
          allocate(grc_pp%npfts(begg:endg)) 

          do p = begp, endp
              c = veg_pp%column(p) 
              g = col_pp%gridcell(c) 
              npfts(g) = npfts(g) + 1 
              grc_pp%npfts(g) = npfts(g)  
          end do 

          maxpfts = maxval(npfts)
          write(iulog,*) "maxpfts",maxpfts 
          
          npfts(:) = 1 
          
          allocate(grc_pp%pfts(maxpfts,begg:endg)) 

          do p = begp, endp
              c = veg_pp%column(p)  
              g = col_pp%gridcell(c) 
              i = npfts(g) 
              grc_pp%pfts(i,g) = p
              npfts(g) = npfts(g) + 1
          end do

         deallocate(npfts(:)) 
         !$acc update device( &
         !$acc        spinup_state            &
         !$acc       ,nyears_ad_carbon_only   &
         !$acc       ,spinup_mortality_factor &
         !$acc       ,carbon_only &
         !$acc       ,carbonphosphorus_only &
         !$acc       ,carbonnitrogen_only &
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
         !$acc       , urban_hac_int &
         !$acc       , const_climate_hist &
         !$acc     )
         !$acc update device(first_step, nlevgrnd) ! TODO!!!!!!, eccen, obliqr, lambm0, mvelpp )
         call update_acc_variables()
         
         !$acc update device ( &
         !$acc       noveg                &
         !$acc       ,ndllf_evr_tmp_tree   &
         !$acc       ,ndllf_evr_brl_tree   &
         !$acc       ,ndllf_dcd_brl_tree   &
         !$acc       ,nbrdlf_evr_trp_tree  &
         !$acc       ,nbrdlf_evr_tmp_tree  &
         !$acc       ,nbrdlf_dcd_trp_tree  &
         !$acc       ,nbrdlf_dcd_tmp_tree  &
         !$acc       ,nbrdlf_dcd_brl_tree  &
         !$acc       ,ntree                &
         !$acc       ,nbrdlf_evr_shrub     &
         !$acc       ,nbrdlf_dcd_tmp_shrub &
         !$acc       ,nbrdlf_dcd_brl_shrub &
         !$acc       ,nc3_arctic_grass     &
         !$acc       ,nc3_nonarctic_grass  &
         !$acc       ,nc4_grass            &
         !$acc       ,npcropmin            &
         !$acc       ,ncorn                &
         !$acc       ,ncornirrig           &
         !$acc       ,nscereal             &
         !$acc       ,nscerealirrig        &
         !$acc       ,nwcereal             &
         !$acc       ,nwcerealirrig        &
         !$acc       ,nsoybean             &
         !$acc       ,nsoybeanirrig        &
         !$acc       ,npcropmax            &
         !$acc       ,nc3crop              &
         !$acc       ,nc3irrig             &
         !$acc       ,nmiscanthus          &
         !$acc       ,nmiscanthusirrig     &
         !$acc       ,nswitchgrass         &
         !$acc       ,nswitchgrassirrig    &
         !$acc       ,num_cfts_known_to_model )
         #if _CUDA 
         istat = cudaMemGetInfo(free2, total)
         write(iulog,*) "Transferred:", free1-free2
         write(iulog,*) "Total:",total
         write(iulog,*) iam,"Free:", free2/1.E+9
         #endif
        call createProcessorFilter(nclumps, bounds_proc, proc_filter, glc2lnd_vars%icemask_grc)
        call createProcessorFilter(nclumps, bounds_proc, proc_filter_inactive_and_active, glc2lnd_vars%icemask_grc)
        
        !$acc enter data copyin(&
      !$acc aerosol_vars     , &
      !$acc AllocParamsInst  , &
      !$acc atm2lnd_vars     , &
      !$acc canopystate_vars, &
      !$acc CH4ParamsInst     , &
      !$acc ch4_vars          , &
      !$acc CNDecompParamsInst     , &
      !$acc CNGapMortParamsInst     , &
      !$acc CNNDynamicsParamsInst     , &
      !$acc cnstate_vars      )
     
      !$acc enter data copyin(&
      !$acc photosyns_vars     , &
      !$acc sedflux_vars     , &
      !$acc soilhydrology_vars     , &
      !$acc soilstate_vars     , &
      !$acc solarabs_vars     , &
      !$acc surfalb_vars     , &
      !$acc surfrad_vars )   
      
      #if _CUDA 
      istat = cudaMemGetInfo(free2, total)
      write(iulog,*) iam,"Free after  moremore_vars:", free2/1.E9
      #endif
      !$acc enter data copyin(&
      !$acc patch_state_updater     , &
      !$acc column_state_updater , &
      !$acc prior_weights ) 
            
      #if _CUDA 
      istat = cudaMemGetInfo(free2, total)
      write(iulog,*) iam,"Free 1st:", free2/1.E9
      #endif
      !$acc enter data copyin(&
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
      !$acc col_ws     ) 
      
      #if _CUDA 
      istat = cudaMemGetInfo(free2, total)
      write(iulog,*) iam,"Free after Col:", free2/1.E9
      #endif
      !$acc enter data copyin( &
      !$acc crop_vars  , &
      !$acc DecompBGCParamsInst     , &
      !$acc DecompCNParamsInst     , &
      !$acc decomp_cascade_con     , &
      !$acc drydepvel_vars     , &
      !$acc dust_vars     , &
      !$acc energyflux_vars     , &
      !$acc frictionvel_vars     , &
      !$acc glc2lnd_vars  ) 
      
      #if _CUDA 
      istat = cudaMemGetInfo(free2, total)
      write(iulog,*) iam,"Free after more vars:", free2/1.E9
      #endif
      !$acc enter data copyin(&
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
      !$acc grc_ws     ) 
      
      #if _CUDA 
      istat = cudaMemGetInfo(free2, total)
      write(iulog,*) iam,"Free after grc_vars:", free2/1.E9
      #endif
      !$acc enter data copyin(&
      !$acc lakestate_vars , &
      !$acc ldomain_gpu  ,&
      !$acc lnd2glc_vars   , &
      !$acc lnd2atm_vars   , &
      !$acc lun_ef     , &
      !$acc lun_es     , &
      !$acc lun_pp     , &
      !$acc lun_ws     , &
      !$acc NitrifDenitrifParamsInst     , &
      !$acc ParamsShareInst     , &
      !$acc params_inst     ) 
      
      #if _CUDA 
      istat = cudaMemGetInfo(free2, total)
      write(iulog,*) iam,"Free after lun/lnd2atm/domain/params:", free2/1.E9
      #endif
      !$acc enter data copyin( subgrid_weights_diagnostics, &
      !$acc top_af     , &
      !$acc top_as     , &
      !$acc top_pp     , &
      !$acc urbanparams_vars , &
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
      !$acc veg_ws      &
      !$acc   )
      #if _CUDA 
      istat = cudaMemGetInfo(free2, total)
      write(iulog,*) iam,"Free after top/veg vars:", free2/1.E9
      #endif
      call htape_gpu_init() 
      !$acc enter data copyin(tape_gpu(:),elmptr_ra(:),elmptr_rs(:))
      !$acc enter data copyin( doalb, declinp1, declin )
      !$acc enter data copyin(filter(:), gpu_clumps(:), gpu_procinfo, proc_filter,proc_filter_inactive_and_active  )
      !$acc enter data copyin(filter_inactive_and_active(:),bounds_proc )
      !$acc enter data copyin(transport_ptr_list(:)) 
      #if _CUDA 
      istat = cudaMemGetInfo(free2, total)
      write(iulog,*) iam,"Free after final copyin:", free2/1.E9
      #endif
      call cpu_time(startt)
      print *, "active only:" 
      call setProcFilters(bounds_proc, proc_filter, .false., glc2lnd_vars%icemask_grc)
      write(iulog,*) "inactive and active:" 
      call setProcFilters(bounds_proc, proc_filter_inactive_and_active, .true., glc2lnd_vars%icemask_grc)
      call cpu_time(stopt) 
      write(iulog,*) iam,"TIMING SetProcFilters :: ",(stopt-startt)*1.E+3, "ms"
      !$acc enter data copyin(cpl_bypass_input%atm_input(:,:,:,1:5))
      end if
      
      !$acc enter data copyin(nstep_mod, dtime_mod, &
      !$acc   year_curr,mon_curr,day_curr,secs_curr,&
      !$acc   year_prev,mon_prev,day_prev,secs_prev, dayspyr_mod,jday_mod)
      
     write(iulog,*) "update_forcings cplbypass : " 
     !call update_forcings_cplbypass(bounds_proc, atm2lnd_vars, cpl_bypass_input, &
     !       dtime_mod, thiscalday_mod,secs_curr, year_curr, mon_curr, nstep_mod) 
      if (do_budgets) call WaterBudget_Reset()
      
      ! ============================================================================
      ! Specified phenology
      ! ============================================================================
      
      if (.not.use_fates) then
         if (use_cn) then
            ! For dry-deposition need to call CLMSP so that mlaidiff is obtained
            if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
               call t_startf('interpMonthlyVeg')
               call interpMonthlyVeg(bounds_proc, canopystate_vars)
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
               call t_startf('interpMonthlyVeg')
               call interpMonthlyVeg(bounds_proc, canopystate_vars)
               call t_stopf('interpMonthlyVeg')
            end if
         end if
      end if
      
      ! ==================================================================================
      ! Determine decomp vertical profiles
      !
      ! These routines (alt_calc & decomp_vertprofiles) need to be called before
      ! pftdyn_cnbal, and it appears that they need to be called before pftdyn_interp and
      ! the associated filter updates, too (otherwise we get a carbon balance error)
      ! ==================================================================================
     !$acc parallel loop independent gang vector default(present) 
     do nc = 1,nclumps
        call alt_calc(filter(nc)%num_soilc, filter(nc)%soilc,canopystate_vars)
     end do
      
      !  Note (WJS, 6-12-13): Because of this routine's placement in the driver sequence
      !  (it is called very early in each timestep, before weights are adjusted and
      !  filters are updated), it may be necessary for this routine to compute values over
      !  inactive as well as active points (since some inactive points may soon become
      !  active) - so that's what is done now. Currently, it seems to be okay to do this,
      !  because the variables computed here seem to only depend on quantities that are
      !  valid over inactive as well as active points.
      if(use_fates .or. use_cn) then 
         call cpu_time(startt)
         call decomp_vertprofiles(bounds_proc,      &
            proc_filter_inactive_and_active%num_soilc, proc_filter_inactive_and_active%soilc, &
            proc_filter_inactive_and_active%num_soilp, proc_filter_inactive_and_active%soilp, &
            soilstate_vars, canopystate_vars, cnstate_vars)
         call cpu_time(stopt)
         write(iulog,*) iam, "TIMING vert_profiles",(stopt-startt)*1.E+3,"ms"
      end if
      ! ============================================================================
      ! Zero fluxes for transient land cover
      ! ============================================================================
      call t_startf('cnpinit')
      call cpu_time(startt) 
         
      call t_startf('beggridwbal')
      call BeginGridWaterBalance(bounds_proc,               &
         proc_filter%num_nolakec, proc_filter%nolakec,       &
         proc_filter%num_lakec, proc_filter%lakec,           &
         proc_filter%num_hydrologyc, proc_filter%hydrologyc, &
         soilhydrology_vars )
      call t_stopf('beggridwbal')
         
      do nc = 1,nclumps
         call get_clump_bounds_gpu(nc, bounds_clump)
         if (use_betr) then
            dtime=get_step_size(); nstep=get_nstep()
            call ep_betr%SetClock(dtime= dtime, nelapstep=nstep)
            call ep_betr%BeginMassBalanceCheck(bounds_clump)
         endif
      end do
      
     call cpu_time(stopt) 
     write(iulog,*) "TIMING BeginGridWaterBalance :: ",(stopt-startt)*1.E+3,"ms"

     call cpu_time(startt)

      if(use_cn) then
         call t_startf('cnpvegzero')
         call zero_elm_weights(bounds_proc)
         call t_stopf('cnpvegzero')
      end if 
      
      call t_startf('cnpvegsumm')
      if(use_cn) then
        call veg_cs_Summary_acc(veg_cs,proc_filter%num_soilp, proc_filter%soilp)
        call veg_ns_Summary_acc(veg_ns,proc_filter%num_soilp, proc_filter%soilp)
        call veg_ps_Summary_acc(veg_ps,proc_filter%num_soilp, proc_filter%soilp)
         
         call summary_veg_state_p2c ( proc_filter%num_soilc, proc_filter%soilc, &
         veg_cs, col_cs, veg_ps,col_ps, veg_ns, col_ns)
         
         ! elseif(use_fates)then
         !      ! In this scenario, we simply zero all of the
         !      ! column level variables that would had been upscaled
         !      ! in the veg summary with p2c
         !      call col_cs%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
         !      call col_ns%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
         !      call col_ps%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
      end if
      call cpu_time(stopt) 
      write(iulog,*) iam,"TIMING cnpvegsum ",(stopt-startt)*1.E+3,"ms"
      
      call t_stopf('cnpvegsumm')
      
      if(use_cn .or. use_fates)then
         !NOTE: potential optimization is to combine these into one subroutine to minimize
         !      number of kernels launched
         call cpu_time(startt) 

         call col_cs_summary_acc(col_cs,proc_filter%num_soilc, proc_filter%soilc)
         call col_ns_Summary_acc(col_ns,proc_filter%num_soilc, proc_filter%soilc)
         call col_ps_Summary_acc(col_ps,proc_filter%num_soilc, proc_filter%soilc)
         call cpu_time(stopt) 
         write(iulog,*) iam, "TIMING col_summary :: ",(stopt-startt)*1.E+3,"ms" 
         
         call cpu_time(startt) 
         call BeginGridCBalance(bounds_proc, col_cs, grc_cs)
         call BeginGridNBalance(bounds_proc, col_ns, grc_ns)
         call BeginGridPBalance(bounds_proc, col_ps, grc_ps)
         call cpu_time(stopt) 
         write(iulog,*) iam,"TIMING cnp_grid_balance :: ",(stopt-startt)*1.E+3,"ms"
      end if
      call t_stopf('cnpinit')
      ! ============================================================================
      ! Update subgrid weights with dynamic landcover (prescribed transient patches,
      ! and or dynamic landunits), and do related adjustments. Note that this
      ! call needs to happen outside loops over nclumps.
      ! ============================================================================
      
     call cpu_time(startt)
     call dyn_hwcontent_init(bounds_proc, &
      proc_filter%num_nolakec, proc_filter%nolakec, &
      proc_filter%num_lakec, proc_filter%lakec, &
      urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars )

     !$acc parallel loop independent gang vector default(present) private(nc,bounds_clump) 
     do nc = 1, nclumps
        call get_clump_bounds_gpu(nc, bounds_clump)
        call set_prior_weights(prior_weights, bounds_clump)
        call set_old_patch_weightsAcc  (patch_state_updater,bounds_clump)
        call set_old_column_weightsAcc (column_state_updater,bounds_clump)
     end do
     call cpu_time(stopt) 
     write(iulog,*) iam, "TIMING dynhwcontent_init :: ",(stopt-startt)*1.E+3,"ms"
     
      ! ==========================================================================
      ! Do land cover change that requires I/O, and thus must be outside a threaded region
      ! ==========================================================================
      
      if (get_do_transient_pfts()) then
         call dynpft_interp(bounds_proc)
      end if
      if (get_do_transient_crops()) then
         call dyncrop_interp(bounds_proc,crop_vars)
      end if
      if (get_do_harvest()) then
         call dynHarvest_interp(bounds_proc)
      end if
      
      ! ==========================================================================
      ! Do everything else related to land cover change
      ! ==========================================================================
      !$acc parallel loop independent gang vector default(present) private(nc, bounds_clump)
      do nc = 1, nclumps
         call get_clump_bounds_gpu(nc, bounds_clump)
         #ifndef _OPENACC 
         if (use_fates) then
         end if
         #endif
         if (create_glacier_mec_landunit) then
            call glc2lnd_vars_update_glc2lnd_acc(glc2lnd_vars ,bounds_clump)
         end if
         
         ! Everything following this point in this loop only needs to be called if we have
         ! actually changed some weights in this time step. This is also required in the
         ! first time step of the run to update filters to reflect state of CISM
         ! (particularly mask that is past through coupler).
         call dynSubgrid_wrapup_weight_changes(bounds_clump, glc2lnd_vars)
         call column_set_new_weightsAcc(column_state_updater,bounds_clump, nc)
      end do
      
      call setProcFilters(bounds_proc, proc_filter, .false.,glc2lnd_vars%icemask_grc)
      
      call patch_set_new_weightsAcc(patch_state_updater ,bounds_proc)
       
      call cpu_time(startt)  
      !$acc parallel loop independent gang vector default(present) private(bounds_clump)
      do nc = 1, nclumps 
         call get_clump_bounds_gpu(nc, bounds_clump) 
         call set_subgrid_diagnostic_fields(bounds_clump)
         
         call initialize_new_columns(bounds_clump, &
            prior_weights%cactive, soilhydrology_vars )
      end do 

      call dyn_hwcontent_final(bounds_proc, &
            proc_filter%num_nolakec, proc_filter%nolakec, &
            proc_filter%num_lakec, proc_filter%lakec, &
            urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
            dtime_mod)

      call cpu_time(stopt) 
      write(iulog,*) iam,"TIMING dyn_hwcontent_final :: ",(stopt-startt)*1.E+3,"ms" 
      call shr_sys_flush(iulog)  
      if (use_cn) then
         call cpu_time(startt) 
         call dyn_cnbal_patch(bounds_proc, &
            proc_filter_inactive_and_active%num_soilp, proc_filter_inactive_and_active%soilp, &
            proc_filter_inactive_and_active%num_soilc, proc_filter_inactive_and_active%soilc, &
            prior_weights, &
            patch_state_updater, &
            canopystate_vars, photosyns_vars, cnstate_vars, &
            veg_cs, c13_veg_cs, c14_veg_cs, &
            veg_ns, veg_ps, dtime_mod) 
          
         call cpu_time(stopt)
         write(iulog,*) iam,"TIMING dyn_cnbal_patch :: ",(stopt-startt)*1.E+3,"ms" 
         if(.not. use_fates ) then 
            call cpu_time(startt) 
            call CarbonStateDynGridUpdate(bounds_proc ,dtime_mod)
            call PhosphorusStateDynGridUpdate(bounds_proc, dtime_mod)
            call NitrogenStateDynGridUpdate(bounds_proc, dtime_mod)
            
            ! Transfer root/seed litter C/N/P to decomposer pools
            call CarbonStateUpdateDynPatch(proc_filter_inactive_and_active%num_soilc,&
            proc_filter_inactive_and_active%soilc,dtime_mod)
            
            call NitrogenStateUpdateDynPatch(proc_filter_inactive_and_active%num_soilc, &
            proc_filter_inactive_and_active%soilc,dtime_mod)
            
            call PhosphorusStateUpdateDynPatch(proc_filter_inactive_and_active%num_soilc, &
            proc_filter_inactive_and_active%soilc,dtime_mod)

            call cpu_time(stopt) 
            write(iulog,*) iam,"TIMING DynGridPatch :: ",(stopt-startt)*1.E+3, "ms"
         end if 
         
      end if
      
      if(use_cn .or. use_fates)then
         call cpu_time(startt) 
        call dyn_col_cs_Adjustments(bounds_proc%begc,bounds_proc%endc, nclumps, column_state_updater, col_cs)
        call dyn_col_ns_Adjustments(bounds_proc%begc,bounds_proc%endc, nclumps, column_state_updater, col_ns)
        call dyn_col_ps_Adjustments(bounds_proc%begc,bounds_proc%endc, nclumps, column_state_updater, col_ps)
         call cpu_time(stopt) 
         write(iulog,*) iam,"TIMING dyn_cnbal_column :: ",(stopt-startt)*1.E+3,"ms"
      end if
      
      if (use_cn  .or. use_fates) then
         nstep = get_nstep()
         
         if (nstep < 2 )then
            if (masterproc) then
               write(iulog,*) '--WARNING-- skipping CN balance check for first timestep'
            end if
         else
            call t_startf('cnbalchk_at_grid')

            call cpu_time(startt) 
            
            if(use_cn) then
              call veg_cs_Summary_acc(veg_cs,proc_filter%num_soilp, proc_filter%soilp)
              call veg_ns_Summary_acc(veg_ns,proc_filter%num_soilp, proc_filter%soilp)
              call veg_ps_Summary_acc(veg_ps,proc_filter%num_soilp, proc_filter%soilp)
               
              call summary_veg_state_p2c ( proc_filter%num_soilc, proc_filter%soilc, &
                     veg_cs, col_cs, veg_ps,col_ps, veg_ns, col_ns)
               
               
            elseif(use_fates)then
               do nc = 1,nclumps
                  call get_clump_bounds(nc, bounds_clump)
                  
                  ! In this scenario, we simply zero all of the
                  ! column level variables that would had been upscaled
                  ! in the veg summary with p2c
                  call col_cs%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
                  call col_ns%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
                  call col_ps%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
               end do 
            end if
            !
            call col_cs_summary_acc(col_cs, proc_filter%num_soilc, proc_filter%soilc)
            call col_ns_Summary_acc(col_ns, proc_filter%num_soilc, proc_filter%soilc)
            call col_ps_Summary_acc(col_ps, proc_filter%num_soilc, proc_filter%soilc)
            
           !$acc parallel loop independent gang vector default(present) private(bounds_clump)
           do nc =1, nclumps 
              call get_clump_bounds_gpu(nc,bounds_clump) 
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
            call t_stopf('cnbalchk_at_grid')
            call cpu_time(stopt) 
            write(iulog,*) iam,"TIMING cnbalchk_at_grid :: ",(stopt-startt)*1.E+3, "ms"
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
         call cpu_time(startt) 
         call t_startf('begwbal')
         !$acc parallel loop independent gang vector default(present) private(nc, bounds_clump)
         do nc = 1,nclumps
            call get_clump_bounds_gpu(nc, bounds_clump)
            
            call BeginColWaterBalance(bounds_clump,                &
               filter(nc)%num_nolakec, filter(nc)%nolakec,       &
               filter(nc)%num_lakec, filter(nc)%lakec,           &
               filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
               soilhydrology_vars )
         end do 
         call t_stopf('begwbal')
            
         call t_startf('begcnpbal')
         ! call veg summary before col summary, for p2c
         if (use_cn) then
              call veg_cs_Summary_acc(veg_cs,proc_filter%num_soilp, proc_filter%soilp)
              call veg_ns_Summary_acc(veg_ns,proc_filter%num_soilp, proc_filter%soilp)
              call veg_ps_Summary_acc(veg_ps,proc_filter%num_soilp, proc_filter%soilp)
            
            call summary_veg_state_p2c(proc_filter%num_soilc,proc_filter%soilc,veg_cs,&
                     col_cs, veg_ps, col_ps, veg_ns, col_ns)

         else if(use_fates)then
               ! In this scenario, we simply zero all of the
               ! column level variables that would had been upscaled
               ! in the veg summary with p2c
               do nc = 1, nclumps 
                  call get_clump_bounds(nc,bounds_clump) 
                  call col_cs%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
                  call col_ns%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
                  call col_ps%ZeroForFates(bounds_clump,filter(nc)%num_soilc, filter(nc)%soilc)
               end do 
         endif
         call t_stopf('begcnpbal')
               
         if (use_cn  .or. use_fates) then
            call t_startf('begcnpbalwf')
            call col_cs_summary_acc(col_cs,proc_filter%num_soilc, proc_filter%soilc)
            call col_ns_Summary_acc(col_ns,proc_filter%num_soilc, proc_filter%soilc)
            call col_ps_Summary_acc(col_ps,proc_filter%num_soilc, proc_filter%soilc)
            
            num_soilc = proc_filter%num_soilc
            !$acc parallel loop independent gang vector default(present) 
            do fc = 1, num_soilc 
               c = proc_filter%soilc(fc) 
               col_cs%begcb(c) = col_cs%totcolc(c) 
               col_ns%begnb(c) = col_ns%totcoln(c) 
               col_ps%begpb(c) = col_ps%totcolp(c) 
            end do 
            
            call t_stopf('begcnpbalwf')
            call cpu_time(stopt) 
            write(iulog,*) iam,"TIMING ColBalanceCheck :: ",(stopt-startt)*1.E+3,"ms"
         end if
         if (do_budgets) then
            call WaterBudget_SetBeginningMonthlyStates(bounds_proc )
            if (use_cn) then
               call CNPBudget_SetBeginningMonthlyStates(bounds_proc, col_cs, grc_cs)
            endif
         endif
               
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

            ! ============================================================================
            ! Initialize variables from previous time step, downscale atm forcings, and
            ! Determine canopy interception and precipitation onto ground surface.
            ! Determine the fraction of foliage covered by water and the fraction
            ! of foliage that is dry and transpiring. Initialize snow layer if the
            ! snow accumulation exceeds 10 mm.
            ! ============================================================================
            
            call cpu_time(startt)
           !$acc parallel loop independent gang private(nc,bounds_clump)
            do nc = 1,nclumps
               call get_clump_bounds_gpu(nc, bounds_clump)
               
               call UpdateDaylength(bounds_clump, declin)
               
               ! Initialze variables needed for new driver time step
               ! elm_drv_init uses bounds_clump in the routine.
               call elm_drv_init(bounds_clump, &
                  filter(nc)%num_nolakec, filter(nc)%nolakec, &
                  filter(nc)%num_nolakep, filter(nc)%nolakep, &
                  filter(nc)%num_soilp  , filter(nc)%soilp,   &
                  canopystate_vars, energyflux_vars)
               
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
                  aerosol_vars )
            end do
            call cpu_time(stopt) 
            write(iulog,*) iam,"TIMING Biogeophys_setup :: ", (stopt-startt)*1.E+3, "ms"
            call cpu_time(startt) 
            !NOTE: canopystate_vars%frac_veg_nosno_alb_patch may be sufficient here?
            call updateFracNoSnoFilters(bounds_proc,proc_filter, canopystate_vars%frac_veg_nosno_patch)
            call cpu_time(stopt)
            write(iulog,*) iam,"TIMING updateFracNoSnoFilters :: ",(stopt-startt)*1.E+3,"ms"

            call cpu_time(startt)
            ! ============================================================================
            ! Surface Radiation
            ! ============================================================================
            
            !call t_startf('surfrad')
            
            ! Surface Radiation primarily for non-urban columns
            
            ! Most of the surface radiation calculations are agnostic to the forest-model
            ! but the calculations of the fractions of sunlit and shaded canopies
            ! are specific, calculate them first.
            ! The nourbanp filter is set in dySubgrid_driver (earlier in this call)
            ! over the patch index range defined by bounds_clump%begp:bounds_proc%endp
            
            ! if(use_fates) then
            !    call alm_fates%wrap_sunfrac(bounds_clump, top_af, canopystate_vars)
            ! else
            call CanopySunShadeFractions(proc_filter%num_nourbanp, proc_filter%nourbanp,    &
                surfalb_vars, canopystate_vars, solarabs_vars)
    
            call SurfaceRadiation(bounds_proc,                &
                proc_filter%num_nourbanp, proc_filter%nourbanp, &
                proc_filter%num_urbanp, proc_filter%urbanp    , &
                proc_filter%num_urbanc, proc_filter%urbanc,     &
                canopystate_vars, surfalb_vars, &
                solarabs_vars, surfrad_vars)
             
             ! Surface Radiation for only urban columns
            call UrbanRadiation(                   &
                proc_filter%num_nourbanl, proc_filter%nourbanl,    &
                proc_filter%num_urbanl, proc_filter%urbanl,        &
                urbanparams_vars, solarabs_vars, surfalb_vars)
            
            call cpu_time(stopt) 
            write(iulog,*) iam,"TIMING Radiation :: ",(stopt-startt)*1.E+3, "ms" 
            call shr_sys_flush(iulog) 
            call cpu_time(startt)
            ! ============================================================================ !
            ! Determine leaf temperature and surface fluxes based on ground                !
            ! temperature from previous time step.                                         !
            ! ============================================================================ !
            call t_startf('bgp1')
            call CanopyTemperature(bounds_proc,                     &
               proc_filter%num_nolakec, proc_filter%nolakec,        &
               proc_filter%num_nolakep, proc_filter%nolakep,        &
               atm2lnd_vars, canopystate_vars, soilstate_vars, frictionvel_vars, &
               energyflux_vars)
            call t_stopf('bgp1')
            
            !$acc parallel loop independent gang vector default(present) 
            do fc = 1, proc_filter%num_nolakec
               c = proc_filter%nolakec(fc)
               col_wf%qflx_snow2topsoi     (c)   = 0._r8
               col_wf%qflx_h2osfc2topsoi   (c)   = 0._r8
            enddo
            
            call cpu_time(stopt)
            write(iulog,*) iam,"TIMING CanopyTemp :: ",(stopt-startt)*1.E+3,"ms"
            ! Bareground fluxes for all patches except lakes and urban landunits
            ! ============================================================================
            ! Determine fluxes
            ! ============================================================================
            if(proc_filter%num_nolu_barep > 0) then 
               call BareGroundFluxes(proc_filter%num_nolu_barep,proc_filter%nolu_barep,&
                  canopystate_vars, soilstate_vars, &
                  frictionvel_vars, ch4_vars)
            endif 
            
            call cpu_time(startt)
            call CanopyFluxes(bounds_proc, proc_filter%num_nolu_barep, proc_filter%nolu_barep, &
                   proc_filter%num_nolu_vegp, proc_filter%nolu_vegp , &
                   canopystate_vars, cnstate_vars  , energyflux_vars, &
                   frictionvel_vars, soilstate_vars, solarabs_vars, surfalb_vars, &
                   ch4_vars, photosyns_vars)
            call cpu_time(stopt)

            write(iulog,*) iam,"TIMING CanopyFluxes :: ",(stopt-startt)*1.E+3,"ms"
            call cpu_time(startt)
            ! Define fields that appear on the restart file for non-urban landunits
            num_nourbanl = proc_filter%num_nourbanl
            !$acc parallel loop independent gang vector default(present) private(l)
            do fl = 1, num_nourbanl
               l = proc_filter%nourbanl(fl)
               lun_es%taf(l) = spval
               lun_ws%qaf(l) = spval
            end do
            
            if(proc_filter%num_urbanl > 0) then
               call UrbanFluxes(bounds_proc,       &
                  proc_filter%num_nourbanl, proc_filter%nourbanl,    &
                  proc_filter%num_urbanl, proc_filter%urbanl,        &
                  proc_filter%num_urbanc, proc_filter%urbanc,        &
                  proc_filter%num_urbanp, proc_filter%urbanp,        &
                  urbanparams_vars, soilstate_vars,  &
                  frictionvel_vars )
            endif
            
            call cpu_time(stopt)
            write(iulog,*) iam, "TIMING UrbanFluxes :: ",(stopt-startt)*1.E+3,"ms"
            call shr_sys_flush(iulog) 
            call cpu_time(startt)
            call LakeFluxes(proc_filter%num_lakep, proc_filter%lakep, &
                       solarabs_vars, frictionvel_vars, &
                       energyflux_vars, lakestate_vars)
            call cpu_time(stopt)
            write(iulog,*)iam, "TIMING LakeFluxes :",(stopt-startt)*1.E+3,"ms"
            ! ============================================================================
            ! DUST and VOC emissions
            ! ============================================================================
            call cpu_time(startt)

            !call t_startf('bgc')
            ! Dust mobilization (C. Zender's modified codes)
            call DustEmission(bounds_proc,      &
            proc_filter%num_nolakep, proc_filter%nolakep, &
            soilstate_vars, canopystate_vars, &
            frictionvel_vars, dust_vars)
            
            ! Dust dry deposition (C. Zender's modified codes)
            call DustDryDep(bounds_proc, frictionvel_vars, dust_vars)
            
            ! VOC emission (A. Guenther's MEGAN (2006) model)
            call shr_sys_flush(iulog) 
            !if (use_voc) then
            !   call VOCEmission(bounds_clump,                                         &
            !        filter(nc)%num_soilp, filter(nc)%soilp,                           &
            !        atm2lnd_vars, canopystate_vars, photosyns_vars, temperature_vars, &
            !        vocemis_vars)
            !end if
            
            call cpu_time(stopt)
            write(iulog,*)iam, "TIMING Dust :: ",(stopt-startt)*1.E+3, "ms"
            
            ! ============================================================================
            ! Determine temperatures
            ! ============================================================================
            call cpu_time(startt)
            call LakeTemperature(bounds_proc, proc_filter%num_lakec, proc_filter%lakec,&
               proc_filter%num_lakep, proc_filter%lakep, &
               solarabs_vars, soilstate_vars, ch4_vars, &
               lakestate_vars)
            call cpu_time(stopt) 
            write(iulog,*) iam, "TIMING LakeTemps :: ",(stopt-startt)*1.E+3,"ms"
            
            ! Set soil/snow temperatures including ground temperature
            call t_startf('soiltemperature')
            call cpu_time(startt)
            call SoilTemperature(bounds_proc,  &
               proc_filter%num_urbanl  , proc_filter%urbanl,   &
               proc_filter%num_nolakec , proc_filter%nolakec,  &
               urbanparams_vars, canopystate_vars, &
               solarabs_vars, soilstate_vars )
            call cpu_time(stopt) 
            write(iulog,*) iam,"TIMING SoilTemp ",(stopt-startt)*1.E+3,"ms"

            call t_stopf('soiltemperature')
            
            call cpu_time(startt) 
            !$acc parallel loop independent gang vector private(nc,bounds_clump)
            do nc = 1,nclumps
               call get_clump_bounds_gpu(nc, bounds_clump)
               ! ============================================================================
               ! update surface fluxes for new ground temperature.
               ! ============================================================================
               !call t_startf('bgp2')
               call SoilFluxes(bounds_clump,                 &
               filter(nc)%num_urbanl,  filter(nc)%urbanl,     &
               filter(nc)%num_nolakec, filter(nc)%nolakec,    &
               filter(nc)%num_nolakep, filter(nc)%nolakep,    &
               atm2lnd_vars, solarabs_vars, canopystate_vars, &
               energyflux_vars )
               !call t_stopf('bgp2')
               
            end do

            !$acc parallel loop independent gang vector private(nc,bounds_clump)
            do nc = 1,nclumps
               call get_clump_bounds_gpu(nc, bounds_clump)
               ! ============================================================================
               ! Perform averaging from patch level to column level
               ! ============================================================================
               !call t_startf('patch2col')
               call elm_drv_patch2col(bounds_clump, filter(nc)%num_nolakec, filter(nc)%nolakec, energyflux_vars)
               !call t_stopf('patch2col')
            end do 
            call cpu_time(stopt) 
            write(iulog,*) iam, "TIMING SoilFlux/p2c ",(stopt-startt)*1.E+3,"ms"

            call cpu_time(outer_start)
            call cpu_time(startt) 
            ! ============================================================================
            ! Vertical (column) soil and surface hydrology
            ! ============================================================================
            ! Note that filter_snowc and filter_nosnowc are returned by
            ! LakeHydrology after the new snow filter is built
            !call t_startf('hydro without drainage')

            call HydrologyNoDrainage(bounds_proc,                   &
                proc_filter%num_nolakec, proc_filter%nolakec,           &
                proc_filter%num_hydrologyc, proc_filter%hydrologyc,     &
                proc_filter%num_hydrononsoic, proc_filter%hydrononsoic, &
                proc_filter%num_urbanc, proc_filter%urbanc,             &
                proc_filter%num_snowc, proc_filter%snowc,               &
                proc_filter%num_nosnowc, proc_filter%nosnowc,           &
                canopystate_vars, atm2lnd_vars, soilstate_vars,  &
                soilhydrology_vars, aerosol_vars )
               
            call cpu_time(stopt) 
            write(iulog,*) iam, "TIMING HydroNoDrainage :: ",(stopt-startt)*1.E+3,"ms" 

            ! Calculate column-integrated aerosol masses, and
            ! mass concentrations for radiative calculations and output
            ! (based on new snow level state, after SnowFilter is rebuilt.
            ! NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there
            ! can be zero snow layers but an active column in filter)
            
            call cpu_time(startt)
            call AerosolMasses( bounds_proc,                                   &
            num_on=proc_filter%num_snowc, filter_on=proc_filter%snowc,       &
            num_off=proc_filter%num_nosnowc, filter_off=proc_filter%nosnowc, &
            aerosol_vars=aerosol_vars)
            call cpu_time(stopt)
            write(iulog,*) iam, "TIMING AerosolMasses :: ",(stopt-startt)*1.E+3,"ms" 

            !call t_stopf('hydro without drainage')
               
            ! ============================================================================
            ! Lake hydrology
            ! ============================================================================
            ! Note that filter_lakesnowc and filter_lakenosnowc are returned by
            ! LakeHydrology after the new snow filter is built
            !call t_startf('hylake')
            
            call cpu_time(startt)
            call LakeHydrology(bounds_proc,                       &
            proc_filter%num_lakec, proc_filter%lakec,                &
            proc_filter%num_lakep, proc_filter%lakep,                &
            proc_filter%num_lakesnowc, proc_filter%lakesnowc,        &
            proc_filter%num_lakenosnowc, proc_filter%lakenosnowc,    &
            atm2lnd_vars, soilstate_vars,  &
            aerosol_vars, lakestate_vars)
            call cpu_time(stopt) 

            !  Calculate column-integrated aerosol masses, and
            !  mass concentrations for radiative calculations and output
            !  (based on new snow level state, after SnowFilter is rebuilt.
            !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there
            !  can be zero snow layers but an active column in filter)
            call cpu_time(startt)
            call AerosolMasses(bounds_proc,                                       &
            num_on=proc_filter%num_lakesnowc, filter_on=proc_filter%lakesnowc,       &
            num_off=proc_filter%num_lakenosnowc, filter_off=proc_filter%lakenosnowc, &
            aerosol_vars=aerosol_vars)
            call cpu_time(stopt) 
            write(iulog,*) iam, "TIMING AerosolMasses :: ",(stopt-startt)*1.E+3,"ms" 

            call cpu_time(startt) 
            ! Must be done here because must use a snow filter for lake columns
            call SnowAge_grain(bounds_proc,                &
            proc_filter%num_lakesnowc, proc_filter%lakesnowc, &
            proc_filter%num_lakenosnowc, proc_filter%lakenosnowc )
               
            call cpu_time(stopt)
            call cpu_time(outer_stop)
            write(iulog,*) iam, "TIMING SnowAgeGrain :: ",(stopt-startt)*1.E+3,"ms" 
            write(iulog,*) iam, "TIMING Hydro-Aerosol :: ",(outer_stop-outer_start)*1.E+3,"ms" 
            ! ============================================================================
            ! ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)
            ! ============================================================================
            call cpu_time(startt) 
            !$acc parallel loop independent gang vector private(c,l)
            do c = bounds_proc%begc,bounds_proc%endc
               l = col_pp%landunit(c)
               if (lun_pp%urbpoi(l)) then
                  ! Urban landunit use Bonan 1996 (LSM Technical Note)
                  col_ws%frac_sno(c) = min(col_ws%snow_depth(c)/0.05_r8, 1._r8)
               end  if
            end do
            call cpu_time(stopt) 
            write(iulog,*) iam,"TIMING frac_sno :: ",(stopt-startt)*1.E+3,"ms" 

            call cpu_time(startt) 
            ! ============================================================================
            ! Snow aging routine based on Flanner and Zender (2006), Linking snowpack
            ! microphysics and albedo evolution, JGR, and Brun (1989), Investigation of
            ! wet-snow metamorphism in respect of liquid-water content, Ann. Glaciol.
            ! ============================================================================
            ! Note the snow filters here do not include lakes
            ! TODO: move this up
            call t_startf('snow_init')
            call SnowAge_grain(bounds_proc,                 &
            proc_filter%num_snowc, proc_filter%snowc,     &
            proc_filter%num_nosnowc, proc_filter%nosnowc )
            call t_stopf('snow_init')
            ! ============================================================================
            ! Update sediment fluxes from land unit
            ! ============================================================================
            if (use_erosion) then
               call t_startf('erosion')
               call SoilErosion(bounds_clump, filter(nc)%num_soilc, filter(nc)%soilc, &
                    atm2lnd_vars, canopystate_vars, soilstate_vars,  sedflux_vars)
               call t_stopf('erosion')
            end if
            ! ============================================================================
            ! Ecosystem dynamics: Uses CN, or static parameterizations
            ! ============================================================================
            call t_startf('ecosysdyn')
            !if (use_cn)then
            !    call crop_vars%CropIncrementYear(filter(nc)%num_pcropp, filter(nc)%pcropp)
            !endif
            call cpu_time(stopt) 
            write(iulog,*) iam,"TIMING SnowAge_grain :: ",(stopt-startt)*1.E+3,"ms" 
            
            ! fully prognostic canopy structure and C-N biogeochemistry
            ! - crop model:  crop algorithms called from within CNEcosystemDyn
            
            !===========================================================================================
            ! elm_interface: 'EcosystemDynNoLeaching' is divided into 2 subroutines (1 & 2): BEGIN
            ! EcosystemDynNoLeaching1 is called before clm_interface
            ! EcosystemDynNoLeaching2 is called after clm_interface
            !===========================================================================================
            call cpu_time(eco_startt) 

            call cpu_time(startt)
            call col_cf_SetValues_acc(col_cf,proc_filter%num_soilc, proc_filter%soilc)
            call col_nf_SetValues_acc(col_nf,proc_filter%num_soilc, proc_filter%soilc)
            call col_pf_SetValues_acc(col_pf,proc_filter%num_soilc, proc_filter%soilc)
            
            call veg_cf_SetValues_acc(veg_cf,proc_filter%num_soilp, proc_filter%soilp, 0._r8)
            call veg_nf_setvalues_acc(veg_nf,proc_filter%num_soilp, proc_filter%soilp, 0._r8)
            call veg_pf_SetValues_acc(veg_pf,proc_filter%num_soilp, proc_filter%soilp, 0._r8)
            call cpu_time(stopt) 
            write(iulog,*) iam, "TIMING SetValuesCols ", (stopt-startt)*1.E+3,"ms"

            call cpu_time(startt)
            call NitrogenDeposition(bounds_proc, atm2lnd_vars)
            
            if ( (.not. nu_com_nfix) .or. use_fates) then
               #ifndef OPENACC
               !soilc loop
               call NitrogenFixation( proc_filter%num_soilc, proc_filter%soilc, dayspyr_mod)
               #endif
            else
               ! nu_com_nfix is true
               !soilc loop
               call NitrogenFixation_balance( proc_filter%num_soilc, proc_filter%soilc, cnstate_vars )
            end if
            call cpu_time(stopt)
            write(iulog,*) iam, "TIMING NDeposition: ",(stopt-startt)*1.E+3,"ms"
            
            call get_proc_bounds(bounds_proc)
            call cpu_time(startt)
            if (crop_prog) then
               call NitrogenFert(proc_filter%num_soilc,proc_filter%soilc )
               
               call CNSoyfix( proc_filter%num_soilc, proc_filter%soilc,&
               proc_filter%num_soilp, proc_filter%soilp, &
               crop_vars, cnstate_vars )
            end if
            ! This is auto-trophic respiration, thus don't call this for FATES
            call MaintenanceResp( proc_filter%num_soilc, proc_filter%soilc, &
            proc_filter%num_soilp, proc_filter%soilp, &
            canopystate_vars, soilstate_vars,  photosyns_vars )
            
            if ( nu_com .ne. 'RD') then
               ! for P competition purpose, calculate P fluxes that will potentially increase solution P pool
               ! then competitors take up solution P
               #ifdef OPENACC 
               call endrun("only RD is supported with OpenACC")
               #endif 
               call PhosphorusWeathering(proc_filter%num_soilc, proc_filter%soilc, cnstate_vars, dtime_mod)
               !NOTE:  nu_com_phosphatase is FALSE for RD
               if (.not. nu_com_phosphatase) then
                  call PhosphorusBiochemMin(proc_filter%num_soilc, proc_filter%soilc, &
                  cnstate_vars, dtime_mod)
               else
                  call PhosphorusBiochemMin_balance(bounds_proc,proc_filter%num_soilc, proc_filter%soilc, &
                  cnstate_vars, dtime_mod)
               end if
            end if
            
            ! --------------------------------------------------
            ! Phosphorus Deposition ! X.SHI
            ! --------------------------------------------------
            call PhosphorusDeposition(bounds_proc,  atm2lnd_vars )
            call cpu_time(stopt)
            write(iulog,*) iam, "TIMING MR/PDeposition: ",(stopt-startt)*1.E+3, "ms"
            !-------------------------------------------------------------------------------------------------
            ! plfotran: 'decomp_rate_constants' must be calculated before entering "clm_interface"
            call cpu_time(startt)
            
            call decomp_rate_constants_cn( proc_filter%num_soilc, proc_filter%soilc, &
            canopystate_vars, soilstate_vars,  ch4_vars, cnstate_vars)
            call cpu_time(stopt)
            write(iulog,*) iam, "TIMING DecompRate:",(stopt-startt)*1.E+3,"ms"
            !-------------------------------------------------------------------------------------------------
            ! 'decomp_vertprofiles' (calc nfixation_prof) is moved from SoilLittDecompAlloc:
            ! ------------------------------------------------------------------------------------------------
            call cpu_time(startt)
            call decomp_vertprofiles(bounds_proc,     &
            proc_filter%num_soilc, proc_filter%soilc, &
            proc_filter%num_soilp, proc_filter%soilp, &
            soilstate_vars, canopystate_vars, cnstate_vars)
            call cpu_time(stopt)
            write(iulog,*) iam, "TIMING EcoVert_profiles",(stopt-startt)*1.E+3,"ms"
            
            ! Allocation1 is always called (w/ or w/o use_elm_interface)
            ! pflotran: call 'Allocation1' to obtain potential N demand for support initial GPP
            if(.not.use_fates)then
               call cpu_time(startt)
               call Allocation1_PlantNPDemand ( &
               proc_filter%num_soilc, proc_filter%soilc, &
               proc_filter%num_soilp, proc_filter%soilp, &
               photosyns_vars, crop_vars, canopystate_vars, cnstate_vars, &
               dtime_mod, year_curr )
               call cpu_time(stopt)
               write(iulog,*) iam, "TIMING AllocationPhase1 : ",(stopt-startt)*1.E+3,"ms"
            end if
            call cpu_time(startt)
            ! directly run elm-bgc
            ! if (use_elm_interface & use_elm_bgc), then CNDecomAlloc is called in elm_driver
           call SoilLittDecompAlloc (bounds_proc, proc_filter%num_soilc, proc_filter%soilc,    &
                 proc_filter%num_soilp, proc_filter%soilp,    &
                 canopystate_vars, soilstate_vars,            &
                 cnstate_vars, ch4_vars, dtime_mod )
            call cpu_time(stopt)
            write(iulog,*) iam, "TIMING SoilLittDecompAlloc:",(stopt-startt)*1.E+3,"ms"
            
            !----------------------------------------------------------------
            ! SoilLittDecompAlloc2 is called by both elm-bgc & pflotran
            ! pflotran: call 'SoilLittDecompAlloc2' to calculate some diagnostic variables and 'fpg' for plant N uptake
            ! pflotran & elm-bgc : 'Allocation3_AG' and vertically integrate net and gross mineralization fluxes
            call cpu_time(startt)
            call SoilLittDecompAlloc2 (proc_filter%num_soilc, proc_filter%soilc,&
            proc_filter%num_soilp, proc_filter%soilp,  &
            canopystate_vars, soilstate_vars,          &
            cnstate_vars, crop_vars, dtime_mod )
            call cpu_time(stopt)
            write(iulog,*) iam, "TIMING SoilLittDecompAlloc2:",(stopt-startt)*1.E+3,"ms"
            
            !--------------------------------------------
            ! Phenology
            !--------------------------------------------
            ! Phenology needs to be called after SoilLittDecompAlloc, because it
            ! depends on current time-step fluxes to new growth on the last
            ! litterfall timestep in deciduous systems
            ! event = 'Phenology'
            ! call t_start_lnd(event)
            call cpu_time(startt)
            call Phenology(proc_filter%num_soilc, proc_filter%soilc, &
            proc_filter%num_soilp, proc_filter%soilp, &
            proc_filter%num_pcropp, proc_filter%pcropp,filter(1)%num_ppercropp,filter(1)%ppercropp, doalb, &
            crop_vars, canopystate_vars, soilstate_vars, cnstate_vars )
            call cpu_time(stopt)
            write(iulog,*) iam, "TIMING Phenology :",(stopt-startt)*1.E+3,"ms"
           
            call cpu_time(startt)  
            call GrowthResp(proc_filter%num_soilp, proc_filter%soilp)
            call veg_cf_summary_rr(veg_cf, proc_filter%num_soilp, proc_filter%soilp, &
                                   proc_filter%num_soilc, proc_filter%soilc, col_cf)
            
            call CarbonStateUpdate0(proc_filter%num_soilp,proc_filter%soilp,veg_cs,veg_cf, dtime_mod)
            !if ( use_c13 ) then
            !    call CarbonStateUpdate0(p,c13_veg_cs,c13_veg_cf, dtime_mod)
            !end if
            !if ( use_c14 ) then
            !    call CarbonStateUpdate0(p,c14_veg_cs,c14_veg_cf, dtime_mod)
            !end if
            call CNLitterToColumn(proc_filter%num_soilc, proc_filter%soilc, cnstate_vars )
            call cpu_time(stopt) 
            write(iulog,*) iam, "TIMING Summary ::",(stopt-startt)*1.E+3,"ms" 
            call cpu_time(startt) 
            call CarbonStateUpdate_Phase1_col(proc_filter%num_soilc, proc_filter%soilc, col_cs, col_cf, dtime_mod)
            call NitrogenStateUpdate_Phase1_col(proc_filter%num_soilc, proc_filter%soilc, cnstate_vars, dtime_mod)
            call PhosphorusStateUpdate_Phase1_col(proc_filter%num_soilc, proc_filter%soilc, cnstate_vars, dtime_mod)
            call cpu_time(stopt) 
            write(iulog,*) "Phase1_col", (stopt-startt)*1.E+3, "ms" 
            if(.not. use_fates) then 
               call CarbonStateUpdate_Phase1_PFT(proc_filter%num_soilp,proc_filter%soilp,crop_vars, veg_cs, veg_cf,dtime_mod)
               call NitrogenStateUpdate_Phase1_pft(proc_filter%num_soilp,proc_filter%soilp, dtime_mod)
               call PhosphorusStateUpdate_Phase1_pft(proc_filter%num_soilp,proc_filter%soilp, dtime_mod) 
            end if
            
            call cpu_time(startt)
            call SoilLittVertTransp( proc_filter%num_soilc, proc_filter%soilc, &
               canopystate_vars, cnstate_vars )
            call cpu_time(stopt)
            write(iulog,*) iam, "TIMING SoilLittVertTransp: ",(stopt-startt)*1.E+3,"ms"

            call cpu_time(startt) 
            call GapMortality( proc_filter%num_soilc, proc_filter%soilc, &
                   proc_filter%num_soilp, proc_filter%soilp,&
                   cnstate_vars )
            !--------------------------------------------
            ! Update2
            !--------------------------------------------
            call CarbonStateUpdate2( proc_filter%num_soilc, proc_filter%soilc, &
                     proc_filter%num_soilp, proc_filter%soilp, &
                     col_cs, veg_cs, col_cf, veg_cf)

            call NitrogenStateUpdate2(proc_filter%num_soilc, proc_filter%soilc, &
            proc_filter%num_soilp, proc_filter%soilp )
            
            call PhosphorusStateUpdate2(proc_filter%num_soilc, proc_filter%soilc, &
            proc_filter%num_soilp, proc_filter%soilp )
            
            call CarbonStateUpdate2h( proc_filter%num_soilc, proc_filter%soilc,  &
            proc_filter%num_soilp, proc_filter%soilp, &
            col_cs, veg_cs, col_cf, veg_cf)
            
            call NitrogenStateUpdate2h(proc_filter%num_soilc, proc_filter%soilc, &
            proc_filter%num_soilp, proc_filter%soilp)
            
            call PhosphorusStateUpdate2h(proc_filter%num_soilc, proc_filter%soilc, proc_filter%num_soilp, proc_filter%soilp)
            call WoodProducts(proc_filter%num_soilc, proc_filter%soilc )
            call CropHarvestPools(proc_filter%num_soilc, proc_filter%soilc, dtime_mod)
            call cpu_time(stopt) 
            write(iulog,*) "TIMING StateUpdate :: ",(stopt-startt)*1.E+3,"ms" 
            
            call cpu_time(startt) 
            call FireArea( proc_filter%num_soilc, proc_filter%soilc, &
                    proc_filter%num_soilp, proc_filter%soilp, &
                    atm2lnd_vars, energyflux_vars, soilhydrology_vars, &
                    cnstate_vars )
             
            call FireFluxes(proc_filter%num_soilc, proc_filter%soilc, &
                    proc_filter%num_soilp, proc_filter%soilp, cnstate_vars)

            call cpu_time(stopt) 
            write(iulog,*) iam, "TIMING FireMod :: ",(stopt-startt)*1.E+3,"ms" 
            !===========================================================================================
            ! elm_interface: 'EcosystemDynNoLeaching' is divided into 2 subroutines (1 & 2): END
            !===========================================================================================
            call cpu_time(startt)
           !$acc parallel loop independent gang vector private(nc,bounds_clump)
           do nc = 1,nclumps
              call get_clump_bounds_gpu(nc, bounds_clump)
              call AnnualUpdate(bounds_clump,            &
              filter(nc)%num_soilc, filter(nc)%soilc, &
              filter(nc)%num_soilp, filter(nc)%soilp, &
              cnstate_vars)
              !call t_stopf('ecosysdyn')
              
              ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)
              !call t_startf('depvel')
              if(.not.use_fates)then
                call depvel_compute(bounds_clump, &
                     atm2lnd_vars, canopystate_vars, frictionvel_vars, &
                     photosyns_vars, drydepvel_vars)
              end if
              !call t_stopf('depvel')
           end do
            call cpu_time(stopt) 
            write(iulog,*) iam, "TIMING AnnualUpdate: ", (stopt-startt)*1.E+3,"ms"
            
            if (use_lch4 .and. .not. is_active_betr_bgc) then
               !warning: do not call ch4 before AnnualUpdate, which will fail the ch4 model
               call t_startf('ch4')
               call CH4 (bounds_proc,                        &
                  proc_filter%num_soilc, proc_filter%soilc,  &
                  proc_filter%num_lakec, proc_filter%lakec,  &
                  proc_filter%num_soilp, proc_filter%soilp,  &
                  lakestate_vars, canopystate_vars, &
                  soilstate_vars, soilhydrology_vars, &
                  energyflux_vars, ch4_vars, lnd2atm_vars)
               call t_stopf('ch4')
            end if

          call t_startf('depvel')
            call cpu_time(startt) 
          !$acc parallel default(present)
          !$acc loop independent gang vector private(nc,bounds_clump)
          do nc = 1,nclumps
             call get_clump_bounds_gpu(nc, bounds_clump)
             ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)
             !if(.not.use_fates)then
             call depvel_compute(bounds_clump, &
             atm2lnd_vars, canopystate_vars, frictionvel_vars, &
             photosyns_vars, drydepvel_vars)
             !end if
          end do
          !$acc end parallel  
           call t_stopf('depvel')
           ! ============================================================================
           ! Calculate soil/snow hydrology with drainage (subsurface runoff)
           ! ============================================================================
           
           call t_startf('hydro2 drainage')
           
           call HydrologyDrainage(bounds_proc,                 &
              proc_filter%num_nolakec, proc_filter%nolakec,       &
              proc_filter%num_hydrologyc, proc_filter%hydrologyc, &
              proc_filter%num_urbanc, proc_filter%urbanc,         &
              proc_filter%num_do_smb_c, proc_filter%do_smb_c,     &
              atm2lnd_vars, glc2lnd_vars,      &
              soilhydrology_vars, soilstate_vars)
            
           call t_stopf('hydro2 drainage')
            call cpu_time(stopt)
            write(iulog,*) iam,"TIMING Depvel/HydroDrainage :: ",(stopt-startt)*1.E+3,"ms"
 
            if (use_cn .or. use_fates) then
               call cpu_time(startt)
               call EcosystemDynLeaching(               &
                  proc_filter%num_soilc, proc_filter%soilc,             &
                  proc_filter%num_soilp, proc_filter%soilp,             &
                  cnstate_vars )
               call cpu_time(stopt)
               write(iulog,*) iam,"TIMING EcosystemDynLeaching :: ",(stopt-startt)*1.E+3,"ms"
            end if
            ! ============================================================================
            ! Update Vegetation
            ! ============================================================================
            ! Execute FATES dynamics
            if ( use_fates ) then
                ! Update high-frequency history diagnostics for FATES
                call alm_fates%wrap_update_hifrq_hist(bounds_clump)
                if ( is_beg_curr_day() ) then ! run ED at the start of each day
                    call alm_fates%dynamics_driv( bounds_clump, top_as,          &
                         top_af, atm2lnd_vars, soilstate_vars, temperature_vars, &
                         canopystate_vars, frictionvel_vars)
                end if
            end if
            call cpu_time(startt)
            if (use_cn .and. doalb) then
                call VegStructUpdate(proc_filter%num_soilp, proc_filter%soilp,   &
                     frictionvel_vars, cnstate_vars, &
                     canopystate_vars, crop_vars, dtime_mod)
            end if
            ! ============================================================================
            ! Check the energy and water balance, also carbon and nitrogen balance
            ! ============================================================================
            !$acc parallel  default(present) 
            !$acc loop independent gang private(nc,bounds_clump)
            do nc = 1,nclumps
               call get_clump_bounds_gpu(nc, bounds_clump)
               !call t_startf('balchk')
               call ColWaterBalanceCheck(bounds_clump, &
               filter(nc)%num_do_smb_c, filter(nc)%do_smb_c, &
               atm2lnd_vars, glc2lnd_vars, solarabs_vars,  &
               energyflux_vars, canopystate_vars)
               !call t_stopf('balchk')
               
               !call t_startf('gridbalchk')
               call GridBalanceCheck(bounds_clump                  , &
               filter(nc)%num_do_smb_c, filter(nc)%do_smb_c   , &
               atm2lnd_vars, glc2lnd_vars, solarabs_vars,       &
               energyflux_vars, canopystate_vars              , &
               soilhydrology_vars)
               !call t_stopf('gridbalchk')
               
               ! if (do_budgets) then
               !    call WaterBudget_SetEndingMonthlyStates(bounds_clump)
               !    if (use_cn) then
               !       call CNPBudget_SetEndingMonthlyStates(bounds_clump, col_cs, grc_cs)
               !    endif
               ! endif
               
               !if (use_cn .or. use_fates) then
               
               !call t_startf('cnbalchk')
               
               call ColCBalanceCheck( &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               col_cs, col_cf)
               
               call ColNBalanceCheck(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               col_ns, col_nf)
               
               call ColPBalanceCheck(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               col_ps, col_pf)
               
               call GridCBalanceCheck(bounds_clump, col_cs, col_cf, grc_cs, grc_cf)
               
               !call t_stopf('cnbalchk')
               !end if
            end do
            !$acc end parallel
            call cpu_time(stopt)
            write(iulog,*) iam,"TIMING VegStruct/BalanceCheck :: ",(stopt-startt)*1.E+3,"ms" 
            ! ============================================================================
            ! Determine albedos for next time step
            ! ============================================================================
            
            if (doalb) then
               !if(nstep_mod >1 ) call write_vars() 
               !$acc parallel loop independent gang private(nc,bounds_clump)
               do nc = 1,nclumps
                  call get_clump_bounds_gpu(nc, bounds_clump)
                  ! Albedos for  non-urban columns
                  !call t_startf('surfalb')
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
                  !call t_stopf('surfalb')
                  
                  ! Albedos for urban columns
                  if (filter_inactive_and_active(nc)%num_urbanl > 0) then
                     !call t_startf('urbsurfalb')
                     call UrbanAlbedo(         &
                     filter_inactive_and_active(nc)%num_urbanl, &
                     filter_inactive_and_active(nc)%urbanl,     &
                     filter_inactive_and_active(nc)%num_urbanc, &
                     filter_inactive_and_active(nc)%urbanc,     &
                     filter_inactive_and_active(nc)%num_urbanp, &
                     filter_inactive_and_active(nc)%urbanp,     &
                     urbanparams_vars, solarabs_vars, surfalb_vars)
                     !call t_stopf('urbsurfalb')
                  end if
               end do
            end if
            
            #ifdef _CUDA
            istat = cudaMemGetInfo(free2, total)
            write(iulog,*) iam,"Free:", free2/1.E9
            #endif
            ! ============================================================================
            ! Determine gridcell averaged properties to send to atm
            ! ============================================================================
            if(use_betr)then
               call ep_betr%DiagnoseLnd2atm(bounds_proc, col_pp, lnd2atm_vars)
            endif
            call t_startf('lnd2atm')
            call cpu_time(startt) 
            call lnd2atm(bounds_proc,       &
               atm2lnd_vars, surfalb_vars, frictionvel_vars,    &
               energyflux_vars, solarabs_vars, drydepvel_vars,  &
               dust_vars, ch4_vars, soilhydrology_vars, lnd2atm_vars)
            call cpu_time(stopt) 
            write(iulog,*) iam, "TIMING lnd2atm :: ",(stopt-startt)*1.E+3,"ms"
            call t_stopf('lnd2atm')
            ! ============================================================================
            ! Determine gridcell averaged properties to send to glc
            ! ============================================================================
            
            if (create_glacier_mec_landunit) then
               call t_startf('lnd2glc')
               do nc = 1,nclumps
                  call get_clump_bounds(nc, bounds_clump)
                  call lnd2glc_vars%update_lnd2glc(bounds_clump,       &
                  filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,   &
                  init=.false.)
               end do
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
            
            if (nstep_mod > 0) then
               
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
               
               call t_stopf('accum')
               
            end if
            
            ! ============================================================================
            ! Update history buffer
            ! ============================================================================
            
            ! Determine if end of history interval
            transfer_tapes = .false.
            do t = 1, ntapes 
               if (tape(t)%nhtfrq==0) then   !monthly average
                  if (mon_curr /= mon_prev) then 
                     tape(t)%is_endhist = .true.
                     transfer_tapes = .true. 
                  end if 
               else
                  write(iulog,*) "nhtfrq: ", tape(t)%nhtfrq 
                  if (mod(nstep_mod,tape(t)%nhtfrq) == 0) then
                     tape(t)%is_endhist = .true.
                     transfer_tapes = .true.
                  end if 
               end if
            end do 
            call t_startf('hbuf')
            !call hist_update_hbuf_gpu(nstep_mod,transfer_tapes, nclumps)
            call hist_update_hbuf(bounds_proc)  
            call t_stopf('hbuf')
            write(iulog,*) iam, "TIMING hist_update_hbuf :: ",(stopt-startt)*1.E+3,"ms"
            
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
            ! Initialization of elm driver variables needed from previous timestep
            !$acc routine seq
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
            
            associate(                                       &
               snl                => col_pp%snl          , & ! Input:  [integer  (:)   ]  number of snow layers
               h2osno             => col_ws%h2osno       , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)
               h2osoi_ice         => col_ws%h2osoi_ice   , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
               h2osoi_liq         => col_ws%h2osoi_liq   , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
               do_capsnow         => col_ws%do_capsnow   , & ! Output: [logical  (:)   ]  true => do snow capping
               h2osno_old         => col_ws%h2osno_old   , & ! Output: [real(r8) (:)   ]  snow water (mm H2O) at previous time step
               frac_iceold        => col_ws%frac_iceold  , & ! Output: [real(r8) (:,:) ]  fraction of ice relative to the tot water
               elai               => canopystate_vars%elai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow
               esai               => canopystate_vars%esai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow
               frac_veg_nosno     => canopystate_vars%frac_veg_nosno_patch     , & ! Output: [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
               frac_veg_nosno_alb => canopystate_vars%frac_veg_nosno_alb_patch , & ! Output: [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
               qflx_glcice        => col_wf%qflx_glcice   , & ! Output: [real(r8) (:)   ]  flux of new glacier ice (mm H2O/s) [+ = ice grows]
               eflx_bot           => col_ef%eflx_bot      , & ! Output: [real(r8) (:)   ]  heat flux from beneath soil/ice column (W/m**2)
               cisun_z            => photosyns_vars%cisun_z_patch , & ! Output: [real(r8) (:)   ]  intracellular sunlit leaf CO2 (Pa)
               cisha_z            => photosyns_vars%cisha_z_patch   & ! Output: [real(r8) (:)   ]  intracellular shaded leaf CO2 (Pa)
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
                  
                  if (.not. use_extrasnowlayers) then
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
               ! NOTE:  Test creation of new p2c_1d_filter that is called within a parallel
               !        loop
               !$acc routine seq
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
                  #ifndef _OPENACC
                  if (.not. use_extrasnowlayers) then
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
                  #ENDIF
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
               
               subroutine zero_elm_weights(bounds)
                  ! Zeroes the weights of variables 
                  ! at start of elm timestep 
                  implicit none 
                  type(bounds_type), intent(in) :: bounds 
                  integer :: p, g ,c, j
                  ! Vegetation
                  !$acc parallel loop independent gang vector default(present) 
                  do p = bounds%begp,bounds%endp
                     !C
                     veg_cs%dispvegc(p) = 0._r8
                     veg_cs%storvegc(p) = 0._r8
                     veg_cs%totpftc(p)  = 0._r8
                     !N
                     veg_ns%dispvegn(p) = 0._r8
                     veg_ns%storvegn(p) = 0._r8
                     veg_ns%totvegn(p)  = 0._r8
                     veg_ns%totpftn(p)  = 0._r8
                     !P
                     veg_ps%dispvegp(p) = 0._r8
                     veg_ps%storvegp(p) = 0._r8
                     veg_ps%totvegp(p)  = 0._r8
                     veg_ps%totpftp(p)  = 0._r8
                  end do
                  
                  ! if (use_c13) then
                  !      call c13_grc_cf%ZeroDWT(bounds_clump)
                  !      call c13_col_cf%ZeroDWT(bounds_clump)
                  ! end if
                  ! if (use_c14) then
                  !      call c14_grc_cf%ZeroDWT(bounds_clump)
                  !      call c14_col_cf%ZeroDWT(bounds_clump)
                  ! end if
                  !$acc parallel loop independent gang vector default(present)
                  do g = bounds%begg, bounds%endg
                     grc_cf%dwt_seedc_to_leaf(g)         = 0._r8
                     grc_cf%dwt_seedc_to_deadstem(g)     = 0._r8
                     grc_cf%dwt_conv_cflux(g)            = 0._r8
                     grc_cf%dwt_prod10c_gain(g)          = 0._r8
                     grc_cf%dwt_prod100c_gain(g)         = 0._r8
                     grc_cf%hrv_deadstemc_to_prod10c(g)  = 0._r8
                     grc_cf%hrv_deadstemc_to_prod100c(g) = 0._r8
                     !
                     grc_nf%dwt_seedn_to_leaf(g)     = 0._r8
                     grc_nf%dwt_seedn_to_deadstem(g) = 0._r8
                     grc_nf%dwt_conv_nflux(g)        = 0._r8
                     grc_nf%dwt_seedn_to_npool(g)    = 0._r8
                     grc_nf%dwt_prod10n_gain(g)      = 0._r8
                     grc_nf%dwt_prod100n_gain(g)     = 0._r8
                     !
                     grc_pf%dwt_seedp_to_leaf(g)     = 0._r8
                     grc_pf%dwt_seedp_to_deadstem(g) = 0._r8
                     grc_pf%dwt_conv_pflux(g)        = 0._r8
                     grc_pf%dwt_seedp_to_ppool(g)    = 0._r8
                     grc_pf%dwt_prod10p_gain(g)      = 0._r8
                     grc_pf%dwt_prod100p_gain(g)     = 0._r8
                     !
                  end do
                  
                  !COLUMN VARIABLES
                  !$acc parallel loop independent gang vector default(present) 
                  do c = bounds%begc,bounds%endc
                     col_cf%dwt_conv_cflux(c)           = 0._r8
                     col_cf%dwt_prod10c_gain(c)         = 0._r8
                     col_cf%dwt_prod100c_gain(c)        = 0._r8
                     col_cf%dwt_crop_productc_gain(c)   = 0._r8
                     col_cf%dwt_slash_cflux(c)          = 0._r8
                     
                     col_nf%dwt_conv_nflux(c)        = 0._r8
                     col_nf%dwt_prod10n_gain(c)      = 0._r8
                     col_nf%dwt_prod100n_gain(c)     = 0._r8
                     col_nf%dwt_crop_productn_gain(c)= 0._r8
                     col_nf%dwt_slash_nflux(c)       = 0._r8
                     !
                     col_pf%dwt_conv_pflux(c)        = 0._r8
                     col_pf%dwt_prod10p_gain(c)      = 0._r8
                     col_pf%dwt_prod100p_gain(c)     = 0._r8
                     col_pf%dwt_crop_productp_gain(c) = 0._r8
                     col_pf%dwt_slash_pflux(c)       = 0._r8
                  end do
                  
                  !$acc parallel loop independent gang worker default(present)
                  do j = 1, nlevdecomp_full
                     !$acc loop vector independent 
                     do c = bounds%begc,bounds%endc
                        col_cf%dwt_frootc_to_litr_met_c(c,j)    = 0._r8
                        col_cf%dwt_frootc_to_litr_cel_c(c,j)    = 0._r8
                        col_cf%dwt_frootc_to_litr_lig_c(c,j)    = 0._r8
                        col_cf%dwt_livecrootc_to_cwdc(c,j)      = 0._r8
                        col_cf%dwt_deadcrootc_to_cwdc(c,j)      = 0._r8
                        !
                        col_nf%dwt_frootn_to_litr_met_n(c,j) = 0._r8
                        col_nf%dwt_frootn_to_litr_cel_n(c,j) = 0._r8
                        col_nf%dwt_frootn_to_litr_lig_n(c,j) = 0._r8
                        col_nf%dwt_livecrootn_to_cwdn(c,j)   = 0._r8
                        col_nf%dwt_deadcrootn_to_cwdn(c,j)   = 0._r8
                        !
                        col_pf%dwt_frootp_to_litr_met_p(c,j) = 0._r8
                        col_pf%dwt_frootp_to_litr_cel_p(c,j) = 0._r8
                        col_pf%dwt_frootp_to_litr_lig_p(c,j) = 0._r8
                        col_pf%dwt_livecrootp_to_cwdp(c,j)   = 0._r8
                        col_pf%dwt_deadcrootp_to_cwdp(c,j)   = 0._r8
                     end do
                  end do
                  
                  
               end subroutine zero_elm_weights
               
            end module elm_driver
            
            
