module clm_driver

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides the main CLM driver physics calling sequence.  Most
  ! computations occurs over ``clumps'' of gridcells (and associated subgrid
  ! scale entities) assigned to each MPI process. Computation is further
  ! parallelized by looping over clumps on each process using shared memory OpenMP.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use clm_varctl             , only : wrtdia, iulog, create_glacier_mec_landunit, use_ed
  use clm_varpar             , only : nlevtrc_soil,nlevsoi
  use clm_varctl             , only : wrtdia, iulog, create_glacier_mec_landunit, use_ed, use_betr  
  use clm_varctl             , only : use_cn, use_cndv, use_lch4, use_voc, use_noio, use_c13, use_c14
  use clm_time_manager       , only : get_step_size, get_curr_date, get_ref_date, get_nstep, is_beg_curr_day, get_curr_time_string
  use clm_varpar             , only : nlevsno, nlevgrnd, crop_prog
  use spmdMod                , only : masterproc, mpicom
  use decompMod              , only : get_proc_clumps, get_clump_bounds, get_proc_bounds, bounds_type
  use filterMod              , only : filter, filter_inactive_and_active
  use histFileMod            , only : hist_update_hbuf, hist_htapes_wrapup
  use restFileMod            , only : restFile_write, restFile_filename
  use abortutils             , only : endrun
  !
  use dynSubgridDriverMod    , only : dynSubgrid_driver
  use BalanceCheckMod        , only : BeginWaterBalance, BalanceCheck
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
  use SurfaceRadiationMod    , only : SurfaceRadiation
  use UrbanRadiationMod      , only : UrbanRadiation
  !clm_bgc_interface
  use CNEcosystemDynMod      , only : CNEcosystemDynNoLeaching1, CNEcosystemDynNoLeaching2

  use CNEcosystemDynMod      , only : CNEcosystemDynLeaching !!CNEcosystemDynNoLeaching,
  use CNVegStructUpdateMod   , only : CNVegStructUpdate 
  use CNAnnualUpdateMod      , only : CNAnnualUpdate
  use CNBalanceCheckMod      , only : BeginCBalance, BeginNBalance, CBalanceCheck, NBalanceCheck
  use CNBalanceCheckMod      , only : BeginPBalance, PBalanceCheck
  use CNVerticalProfileMod   , only : decomp_vertprofiles
  use CNFireMod              , only : CNFireInterp
  use CNDVDriverMod          , only : CNDVDriver, CNDVHIST
  use SatellitePhenologyMod  , only : SatellitePhenology, interpMonthlyVeg
  use ndepStreamMod          , only : ndep_interp
  use pdepStreamMod          , only : pdep_interp
  use ActiveLayerMod         , only : alt_calc
  use ch4Mod                 , only : ch4
  use DUSTMod                , only : DustDryDep, DustEmission
  use VOCEmissionMod         , only : VOCEmission
  !
  use filterMod              , only : setFilters
  use EDmainMod              , only : edmodel
  use EDCLMLinkMod           , only : clm_ed_link
  use EDtypesMod             , only : gridCellEdState
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
  use clm_instMod            , only : ch4_vars, ep_betr
  use clm_instMod            , only : carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars
  use clm_instMod            , only : carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars
  use clm_instMod            , only : nitrogenstate_vars
  use clm_instMod            , only : nitrogenflux_vars
  use clm_instMod            , only : phosphorusstate_vars
  use clm_instMod            , only : phosphorusflux_vars
  use clm_instMod            , only : dgvs_vars
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
  use clm_instMod            , only : urbanparams_vars
  use clm_instMod            , only : waterflux_vars
  use clm_instMod            , only : waterstate_vars
  use clm_instMod            , only : atm2lnd_vars
  use clm_instMod            , only : lnd2atm_vars
  use clm_instMod            , only : glc2lnd_vars
  use clm_instMod            , only : lnd2glc_vars
  use clm_instMod            , only : EDbio_vars
  use clm_instMod            , only : soil_water_retention_curve
  use clm_instMod            , only : chemstate_vars
!  use betr_initializeMod     , only : betrtracer_vars
!  use betr_initializeMod     , only : tracercoeff_vars
!  use betr_initializeMod     , only : tracerflux_vars
!  use betr_initializeMod     , only : tracerState_vars
!  use betr_initializeMod     , only : tracerboundarycond_vars
!  use betr_initializeMod     , only : bgc_reaction
!  use betr_initializeMod     , only : plantsoilnutrientflux_vars
!  use BetrBGCMod             , only : run_betr_one_step_without_drainage
!  use BetrBGCMod             , only : run_betr_one_step_with_drainage
!  use TracerBalanceMod       , only : betr_tracer_massbalance_check
!  use TracerBalanceMod       , only : begin_betr_tracer_massbalance
  use tracer_varcon          , only : is_active_betr_bgc, do_betr_leaching
!  use CNEcosystemDynBetrMod  , only : CNEcosystemDynBetrVeg, CNEcosystemDynBetrSummary, CNFluxStateBetrSummary
  use GridcellType           , only : grc                
  use LandunitType           , only : lun                
  use ColumnType             , only : col                
  use PatchType              , only : pft
  use shr_sys_mod            , only : shr_sys_flush

  !!----------------------------------------------------------------------------
  !! bgc interface & pflotran:
  use clm_varctl             , only : use_bgc_interface
  use clm_instMod            , only : clm_bgc_data
  use clm_bgc_interfaceMod   , only : get_clm_bgc_data
  !! (1) clm_bgc through interface
  use clm_varctl             , only : use_clm_bgc
  use clm_bgc_interfaceMod   , only : clm_bgc_run, update_bgc_data_clm2clm
  !! (2) pflotran
  use clm_varctl             , only : use_pflotran, pf_cmode, pf_hmode, pf_tmode
  use clm_bgc_interfaceMod   , only : update_bgc_data_pf2clm
  use clm_pflotran_interfaceMod   , only : clm_pf_run, clm_pf_write_restart
!  use clm_pflotran_interfaceMod   , only : clm_pf_finalize
  !!----------------------------------------------------------------------------

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
  subroutine clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate)
    !
    ! !DESCRIPTION:
    !
    ! First phase of the clm driver calling the clm physics. An outline of
    ! the calling tree is given in the description of this module.
    !
    ! !USES:
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

    ! Update time-related info

    call cnstate_vars%CropRestIncYear()

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

       if (use_cn) then
          !  Note (WJS, 6-12-13): Because of this routine's placement in the driver sequence
          !  (it is called very early in each timestep, before weights are adjusted and
          !  filters are updated), it may be necessary for this routine to compute values over
          !  inactive as well as active points (since some inactive points may soon become
          !  active) - so that's what is done now. Currently, it seems to be okay to do this,
          !  because the variables computed here seem to only depend on quantities that are
          !  valid over inactive as well as active points.

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
    ! Initialize the mass balance checks for carbon and nitrogen, and zero fluxes for
    ! transient land cover
    ! ============================================================================

    !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       if (use_betr) then
         dtime=get_step_size(); nstep=get_nstep()
         call ep_betr%SetClock(dtime= dtime, nelapstep=nstep)
         call ep_betr%BeginMassBalanceCheck(bounds_clump)
       endif
       
       if (use_cn) then
          call t_startf('begcnbal')

          call BeginCBalance(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               carbonstate_vars)

          call BeginNBalance(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               nitrogenstate_vars)

          

          call BeginPBalance(bounds_clump, &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               phosphorusstate_vars)


          call carbonflux_vars%ZeroDWT(bounds_clump)
          if (use_c13) then
             call c13_carbonflux_vars%ZeroDWT(bounds_clump)
          end if
          if (use_c14) then
             call c14_carbonflux_vars%ZeroDWT(bounds_clump)
          end if
          call nitrogenflux_vars%ZeroDWT(bounds_clump)
          call phosphorusflux_vars%ZeroDWT(bounds_clump)

          call carbonstate_vars%ZeroDWT(bounds_clump)
          call nitrogenstate_vars%ZeroDWT(bounds_clump)
          call phosphorusstate_vars%ZeroDWT(bounds_clump)

          call t_stopf('begcnbal')
       end if
    end do
    !$OMP END PARALLEL DO

    ! ============================================================================
    ! Update subgrid weights with dynamic landcover (prescribed transient patches,
    ! CNDV, and or dynamic landunits), and do related adjustments. Note that this
    ! call needs to happen outside loops over nclumps.
    ! ============================================================================

    call t_startf('dyn_subgrid')
    call dynSubgrid_driver(bounds_proc,                                      &
       urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
       waterstate_vars, waterflux_vars, temperature_vars, energyflux_vars,   &
       canopystate_vars, photosyns_vars, cnstate_vars, dgvs_vars,            &
       carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars,         &
       carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars,            &
       nitrogenstate_vars, nitrogenflux_vars, glc2lnd_vars,                  &
       phosphorusstate_vars,phosphorusflux_vars)
    call t_stopf('dyn_subgrid')

    ! ============================================================================
    ! Initialize the mass balance checks for water.
    !
    ! Currently, I believe this needs to be done after weights are updated for
    ! prescribed transient patches or CNDV, because column-level water is not
    ! generally conserved when weights change (instead the difference is put in
    ! the grid cell-level terms, qflx_liq_dynbal, etc.). In the future, we may
    ! want to change the balance checks to ensure that the grid cell-level water
    ! is conserved, considering qflx_liq_dynbal; in this case, the call to
    ! BeginWaterBalance should be moved to before the weight updates.
    ! ============================================================================

    !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       call t_startf('begwbal')
       call BeginWaterBalance(bounds_clump,                   &
            filter(nc)%num_nolakec, filter(nc)%nolakec,       &
            filter(nc)%num_lakec, filter(nc)%lakec,           &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
            soilhydrology_vars, waterstate_vars)
       call t_stopf('begwbal')
    end do
    !$OMP END PARALLEL DO

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
       call CNFireInterp(bounds_proc)
       call t_stopf('ndep_interp')
    end if

    ! ============================================================================
    ! Update dynamic P deposition field, on albedo timestep
    ! currently being done outside clumps loop, but no reason why it couldn't be
    ! re-written to go inside.
    ! ============================================================================

    if (use_cn) then
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

    !$OMP PARALLEL DO PRIVATE (nc,l,c, bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       call t_startf('drvinit')

       call UpdateDaylength(bounds_clump, declin)

       ! Initialze variables needed for new driver time step 
       call clm_drv_init(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_nolakep, filter(nc)%nolakep, &
            filter(nc)%num_soilp  , filter(nc)%soilp,   &
            canopystate_vars, waterstate_vars, waterflux_vars, energyflux_vars)

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
            atm2lnd_vars, canopystate_vars, temperature_vars, &
            aerosol_vars, waterstate_vars, waterflux_vars)
       call t_stopf('canhydro')

       ! ============================================================================
       ! Surface Radiation
       ! ============================================================================

       call t_startf('surfrad')

       ! Surface Radiation primarily for non-urban columns 

       call SurfaceRadiation(bounds_clump,                                 &
            filter(nc)%num_nourbanp, filter(nc)%nourbanp,                  &
            filter(nc)%num_urbanp, filter(nc)%urbanp    ,                  &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                      &
            atm2lnd_vars, waterstate_vars, canopystate_vars, surfalb_vars, &
            solarabs_vars, surfrad_vars)

       ! Surface Radiation for only urban columns

       call UrbanRadiation(bounds_clump,                                       &
            filter(nc)%num_nourbanl, filter(nc)%nourbanl,                      &
            filter(nc)%num_urbanl, filter(nc)%urbanl,                          &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                          &
            filter(nc)%num_urbanp, filter(nc)%urbanp,                          &
            atm2lnd_vars, waterstate_vars, temperature_vars, urbanparams_vars, &
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
            waterstate_vars, waterflux_vars, energyflux_vars, temperature_vars)
       call t_stopf('bgp1')

       ! ============================================================================
       ! Determine fluxes
       ! ============================================================================

       call t_startf('bgflux')

       call waterflux_vars%Reset(bounds_clump, filter(nc)%num_nolakec , filter(nc)%nolakec)       
       ! Bareground fluxes for all patches except lakes and urban landunits

       call BareGroundFluxes(bounds_clump,                                 &
            filter(nc)%num_nolakeurbanp, filter(nc)%nolakeurbanp,          &
            atm2lnd_vars, canopystate_vars, soilstate_vars,                &
            frictionvel_vars, ch4_vars, energyflux_vars, temperature_vars, &
            waterflux_vars, waterstate_vars)
       call t_stopf('bgflux')

       ! non-bareground fluxes for all patches except lakes and urban landunits
       ! Calculate canopy temperature, latent and sensible fluxes from the canopy,
       ! and leaf water change by evapotranspiration

       call t_startf('canflux')
       call CanopyFluxes(bounds_clump,                                                   &
            filter(nc)%num_nolakeurbanp, filter(nc)%nolakeurbanp,                        &
            atm2lnd_vars, canopystate_vars, cnstate_vars, energyflux_vars,               &
            frictionvel_vars, soilstate_vars, solarabs_vars, surfalb_vars,               &
            temperature_vars, waterflux_vars, waterstate_vars, ch4_vars, photosyns_vars, &
            EDbio_vars, soil_water_retention_curve, nitrogenstate_vars,phosphorusstate_vars) 
       call t_stopf('canflux')

       ! Fluxes for all urban landunits

       call t_startf('uflux')
       call UrbanFluxes(bounds_clump,                                         &
            filter(nc)%num_nourbanl, filter(nc)%nourbanl,                     &
            filter(nc)%num_urbanl, filter(nc)%urbanl,                         &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                         &
            filter(nc)%num_urbanp, filter(nc)%urbanp,                         &
            atm2lnd_vars, urbanparams_vars, soilstate_vars, temperature_vars, &
            waterstate_vars, frictionvel_vars, energyflux_vars, waterflux_vars) 
       call t_stopf('uflux')

       ! Fluxes for all lake landunits

       call t_startf('bgplake')
       call LakeFluxes(bounds_clump,                                         &
            filter(nc)%num_lakec, filter(nc)%lakec,                          &
            filter(nc)%num_lakep, filter(nc)%lakep,                          &
            atm2lnd_vars, solarabs_vars, frictionvel_vars, temperature_vars, &
            energyflux_vars, waterstate_vars, waterflux_vars, lakestate_vars) 

       ! ============================================================================
       ! DUST and VOC emissions
       ! ============================================================================

       call t_startf('bgc')

       ! Dust mobilization (C. Zender's modified codes)
       call DustEmission(bounds_clump,                                       &
            filter(nc)%num_nolakep, filter(nc)%nolakep,                      &
            atm2lnd_vars, soilstate_vars, canopystate_vars, waterstate_vars, &
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

       ! Set lake temperature 

       call LakeTemperature(bounds_clump,                                             &
            filter(nc)%num_lakec, filter(nc)%lakec,                                   &
            filter(nc)%num_lakep, filter(nc)%lakep,                                   & 
            solarabs_vars, soilstate_vars, waterstate_vars, waterflux_vars, ch4_vars, &
            energyflux_vars, temperature_vars, lakestate_vars)
       call t_stopf('bgplake')

       ! Set soil/snow temperatures including ground temperature

       call t_startf('soiltemperature')
       call SoilTemperature(bounds_clump,                                                      &
            filter(nc)%num_urbanl  , filter(nc)%urbanl,                                        &
            filter(nc)%num_nolakec , filter(nc)%nolakec,                                       &
            atm2lnd_vars, urbanparams_vars, canopystate_vars, waterstate_vars, waterflux_vars, &
            solarabs_vars, soilstate_vars, energyflux_vars,  temperature_vars)
       call t_stopf('soiltemperature')

       ! ============================================================================
       ! update surface fluxes for new ground temperature.
       ! ============================================================================

       call t_startf('bgp2')
       call SoilFluxes(bounds_clump,                                                          &
            filter(nc)%num_urbanl,  filter(nc)%urbanl,                                        &
            filter(nc)%num_nolakec, filter(nc)%nolakec,                                       &
            filter(nc)%num_nolakep, filter(nc)%nolakep,                                       &
            atm2lnd_vars, solarabs_vars, temperature_vars, canopystate_vars, waterstate_vars, &
            energyflux_vars, waterflux_vars)            
       call t_stopf('bgp2')

       ! ============================================================================
       ! Perform averaging from patch level to column level
       ! ============================================================================

       call t_startf('patch2col')
       call clm_drv_patch2col(bounds_clump, filter(nc)%num_nolakec, filter(nc)%nolakec, &
            waterstate_vars, energyflux_vars, waterflux_vars)
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
            filter(nc)%num_urbanc, filter(nc)%urbanc,                        &
            filter(nc)%num_snowc, filter(nc)%snowc,                          &
            filter(nc)%num_nosnowc, filter(nc)%nosnowc,                      &
            atm2lnd_vars, soilstate_vars, energyflux_vars, temperature_vars, &
            waterflux_vars, waterstate_vars, soilhydrology_vars, aerosol_vars, &
            soil_water_retention_curve,ep_betr) 

       !  Calculate column-integrated aerosol masses, and
       !  mass concentrations for radiative calculations and output
       !  (based on new snow level state, after SnowFilter is rebuilt.
       !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there 
       !  can be zero snow layers but an active column in filter)
      
       call AerosolMasses( bounds_clump,                                   &
            num_on=filter(nc)%num_snowc, filter_on=filter(nc)%snowc,       &
            num_off=filter(nc)%num_nosnowc, filter_off=filter(nc)%nosnowc, &
            waterflux_vars=waterflux_vars,                                 &
            waterstate_vars=waterstate_vars,                               &
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
            atm2lnd_vars, temperature_vars, soilstate_vars, waterstate_vars, waterflux_vars, &
            energyflux_vars, aerosol_vars, lakestate_vars)
       
       !  Calculate column-integrated aerosol masses, and
       !  mass concentrations for radiative calculations and output
       !  (based on new snow level state, after SnowFilter is rebuilt.
       !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there 
       !  can be zero snow layers but an active column in filter)

       call AerosolMasses(bounds_clump,                                            &
            num_on=filter(nc)%num_lakesnowc, filter_on=filter(nc)%lakesnowc,       &
            num_off=filter(nc)%num_lakenosnowc, filter_off=filter(nc)%lakenosnowc, &
            waterflux_vars=waterflux_vars,                                         &
            waterstate_vars=waterstate_vars,                                       &
            aerosol_vars=aerosol_vars)                      

       ! Must be done here because must use a snow filter for lake columns

       call SnowAge_grain(bounds_clump,                         &
            filter(nc)%num_lakesnowc, filter(nc)%lakesnowc,     &
            filter(nc)%num_lakenosnowc, filter(nc)%lakenosnowc, &
            waterflux_vars, waterstate_vars, temperature_vars)

       call t_stopf('hylake')

       ! ============================================================================
       ! ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)
       ! ============================================================================

       do c = bounds_clump%begc,bounds_clump%endc
          l = col%landunit(c)
          if (lun%urbpoi(l)) then
             ! Urban landunit use Bonan 1996 (LSM Technical Note)
             waterstate_vars%frac_sno_col(c) = min( waterstate_vars%snow_depth_col(c)/0.05_r8, 1._r8)
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
            filter(nc)%num_nosnowc, filter(nc)%nosnowc, &
            waterflux_vars, waterstate_vars, temperature_vars)
       call t_stopf('snow_init')

       ! ============================================================================
       ! Ecosystem dynamics: Uses CN, CNDV, or static parameterizations
       ! ============================================================================

       call t_startf('ecosysdyn')

       if(is_active_betr_bgc)then
         !right now betr bgc is intended only for non-ed mode
         
         !this returns the plant nutrient demand to soil bgc
!         call CNEcosystemDynBetrVeg(bounds_clump,                                    &
!                  filter(nc)%num_soilc, filter(nc)%soilc,                        &
!                  filter(nc)%num_soilp, filter(nc)%soilp,                        &
!                  filter(nc)%num_pcropp, filter(nc)%pcropp, doalb,               &
!                  cnstate_vars, carbonflux_vars, carbonstate_vars,               &
!                  c13_carbonflux_vars, c13_carbonstate_vars,                     &
!                  c14_carbonflux_vars, c14_carbonstate_vars,                     &
!                  nitrogenflux_vars, nitrogenstate_vars,                         &
!!                  atm2lnd_vars, waterstate_vars, waterflux_vars,                 &
!                  canopystate_vars, soilstate_vars, temperature_vars, crop_vars, &
!                  dgvs_vars, photosyns_vars, soilhydrology_vars, energyflux_vars,&
!                  plantsoilnutrientflux_vars,                                    &
!                  phosphorusflux_vars, phosphorusstate_vars)

         !do belowground bgc and transport         
         call t_startf('betr_nodrain')

!         call run_betr_one_step_without_drainage(bounds_clump, 1, nlevtrc_soil,      &
!              filter(nc)%num_soilc, filter(nc)%soilc,                                &
!              filter(nc)%num_soilp, filter(nc)%soilp,                                &
!              col, atm2lnd_vars,                                                     &
!              soilhydrology_vars, soilstate_vars, waterstate_vars, temperature_vars, &
!              waterflux_vars, chemstate_vars, cnstate_vars, canopystate_vars,        &
!              carbonstate_vars, carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars,&
!              betrtracer_vars, bgc_reaction,     &
!              tracerboundarycond_vars, tracercoeff_vars, tracerstate_vars,           &
!              tracerflux_vars, plantsoilnutrientflux_vars)

         call t_stopf('betr_nodrain')
         
         !do ecosystem variable summary
!         call CNEcosystemDynBetrSummary(bounds_clump,                                &
!                  filter(nc)%num_soilc, filter(nc)%soilc,                        &
!                  filter(nc)%num_soilp, filter(nc)%soilp,                        &
!                  filter(nc)%num_pcropp, filter(nc)%pcropp, doalb,               &
!                  cnstate_vars, carbonflux_vars, carbonstate_vars,               &
!                  c13_carbonflux_vars, c13_carbonstate_vars,                     &
!                  c14_carbonflux_vars, c14_carbonstate_vars,                     &
!                  nitrogenflux_vars, nitrogenstate_vars,                         &
!                  atm2lnd_vars, waterstate_vars, waterflux_vars,                 &
!                  canopystate_vars, soilstate_vars, temperature_vars, crop_vars, &
!                  dgvs_vars, photosyns_vars, soilhydrology_vars, energyflux_vars,&
!                  plantsoilnutrientflux_vars, phosphorusstate_vars)  
       else       
       
         ! FIX(SPM,032414)  push these checks into the routines below and/or make this consistent.
         if (.not. use_ed) then 
           if (use_cn) then 

             ! fully prognostic canopy structure and C-N biogeochemistry
             ! - CNDV defined: prognostic biogeography; else prescribed
             ! - crop model:  crop algorithms called from within CNEcosystemDyn
             
             !!===========================================================================================
             !! clm_bgc_interface: 'CNEcosystemDynNoLeaching' is divided into 2 subroutines (1 & 2): BEGIN
             !! CNEcosystemDynNoLeaching1 is called before clm_bgc_interface
             !! CNEcosystemDynNoLeaching2 is called after clm_bgc_interface
             !!===========================================================================================
             call CNEcosystemDynNoLeaching1(bounds_clump,                               &
                       filter(nc)%num_soilc, filter(nc)%soilc,                          &
                       filter(nc)%num_soilp, filter(nc)%soilp,                          &
                       cnstate_vars, carbonflux_vars, carbonstate_vars,                 &
                       c13_carbonflux_vars,                                             &
                       c14_carbonflux_vars,                                             &
                       nitrogenflux_vars, nitrogenstate_vars,                           &
                       atm2lnd_vars, waterstate_vars, waterflux_vars,                   &
                       canopystate_vars, soilstate_vars, temperature_vars, crop_vars,   &
                       ch4_vars, photosyns_vars,                                        &
                       phosphorusflux_vars,phosphorusstate_vars)

             !!--------------------------------------------------------------------------------
             if (use_bgc_interface) then
                 !! STEP-1: pass data from CLM to clm_bgc_data (INTERFACE DATA TYPE)
                 call get_clm_bgc_data(clm_bgc_data,bounds_clump,                       &
                           filter(nc)%num_soilc, filter(nc)%soilc,                      &
                           filter(nc)%num_soilp, filter(nc)%soilp,                      &
                           atm2lnd_vars, waterstate_vars, waterflux_vars,               &
                           soilstate_vars,  temperature_vars, energyflux_vars,          &
                           soilhydrology_vars, soil_water_retention_curve,              &
                           cnstate_vars, carbonflux_vars, carbonstate_vars,             &
                           nitrogenflux_vars, nitrogenstate_vars,                       &
                           phosphorusflux_vars, phosphorusstate_vars,                   &
                           canopystate_vars, ch4_vars)

                 if (use_pflotran .and. pf_cmode) then
                    call t_startf('pflotran')
                    ! -------------------------------------------------------------------------
                    ! PFLOTRAN calling for solving below-ground and ground-surface processes,
                    ! including thermal, hydrological and biogeochemical processes
                    !! STEP-2: (1) pass data from clm_bgc_data to pflotran
                    !! STEP-2: (2) run pflotran
                    !! STEP-2: (3) update clm_bgc_data from pflotran
                    ! -------------------------------------------------------------------------
                    call clm_pf_run(clm_bgc_data,bounds_clump,                  &
                           filter(nc)%num_soilc, filter(nc)%soilc)

                    !! STEP-3: update CLM from clm_bgc_data
                    call update_bgc_data_pf2clm(clm_bgc_data,bounds_clump,      &
                           filter(nc)%num_soilc, filter(nc)%soilc,              &
                           filter(nc)%num_soilp, filter(nc)%soilp,              &
                           atm2lnd_vars,                                        &
                           waterstate_vars, waterflux_vars,                     &
                           soilstate_vars,  temperature_vars, energyflux_vars,  &
                           soilhydrology_vars, soil_water_retention_curve,      &
                           cnstate_vars, carbonflux_vars, carbonstate_vars,     &
                           nitrogenflux_vars, nitrogenstate_vars,               &
                           phosphorusflux_vars, phosphorusstate_vars,           &
                           ch4_vars)

                    call t_stopf('pflotran')

                 elseif (use_clm_bgc) then
                    call t_startf('clm-bgc via interface')
                    ! -------------------------------------------------------------------------
                    !! run clm-bgc (CNDecompAlloc) through interface
                    !! STEP-2: (1) pass data from clm_bgc_data to CNDecompAlloc
                    !! STEP-2: (2) run CNDecompAlloc
                    !! STEP-2: (3) update clm_bgc_data from CNDecompAlloc
                    ! -------------------------------------------------------------------------
                    call clm_bgc_run(clm_bgc_data, bounds_clump,                &
                           filter(nc)%num_soilc, filter(nc)%soilc,              &
                           filter(nc)%num_soilp, filter(nc)%soilp,              &
                           canopystate_vars, soilstate_vars,                    &
                           temperature_vars, waterstate_vars,                   &
                           cnstate_vars, ch4_vars,                              &
                           carbonstate_vars, carbonflux_vars,                   &
                           nitrogenstate_vars, nitrogenflux_vars,               &
                           phosphorusstate_vars,phosphorusflux_vars)

                    !! STEP-3: update CLM from clm_bgc_data
                    call update_bgc_data_clm2clm(clm_bgc_data, bounds_clump,    &
                           filter(nc)%num_soilc, filter(nc)%soilc,              &
                           filter(nc)%num_soilp, filter(nc)%soilp,              &
                           atm2lnd_vars,                                        &
                           waterstate_vars, waterflux_vars,                     &
                           soilstate_vars,  temperature_vars, energyflux_vars,  &
                           soilhydrology_vars, soil_water_retention_curve,      &
                           cnstate_vars, carbonflux_vars, carbonstate_vars,     &
                           nitrogenflux_vars, nitrogenstate_vars,               &
                           phosphorusflux_vars, phosphorusstate_vars,           &
                           ch4_vars)
                    call t_stopf('clm-bgc via interface')
                 end if !!if (use_pflotran .and. pf_cmode)
             end if !!if (use_bgc_interface)
             !!--------------------------------------------------------------------------------

             call CNEcosystemDynNoLeaching2(bounds_clump,                                   &
                   filter(nc)%num_soilc, filter(nc)%soilc,                                  &
                   filter(nc)%num_soilp, filter(nc)%soilp,                                  &
                   filter(nc)%num_pcropp, filter(nc)%pcropp, doalb,                         &
                   cnstate_vars, carbonflux_vars, carbonstate_vars,                         &
                   c13_carbonflux_vars, c13_carbonstate_vars,                               &
                   c14_carbonflux_vars, c14_carbonstate_vars,                               &
                   nitrogenflux_vars, nitrogenstate_vars,                                   &
                   atm2lnd_vars, waterstate_vars, waterflux_vars,                           &
                   canopystate_vars, soilstate_vars, temperature_vars, crop_vars, ch4_vars, &
                   dgvs_vars, photosyns_vars, soilhydrology_vars, energyflux_vars,          &
                   phosphorusflux_vars,phosphorusstate_vars)

             !!===========================================================================================
             !! clm_bgc_interface: 'CNEcosystemDynNoLeaching' is divided into 2 subroutines (1 & 2): END
             !!===========================================================================================

             call CNAnnualUpdate(bounds_clump,            &
                  filter(nc)%num_soilc, filter(nc)%soilc, &
                  filter(nc)%num_soilp, filter(nc)%soilp, &
                  cnstate_vars, carbonflux_vars)
          else ! not use_cn

             if (doalb) then
                ! Prescribed biogeography - prescribed canopy structure, some prognostic carbon fluxes

                call SatellitePhenology(bounds_clump,               &
                     filter(nc)%num_nolakep, filter(nc)%nolakep,    &
                     waterstate_vars, canopystate_vars)
             end if

           end if  ! end of if-use_cn

         else ! use_ed

           call carbonflux_vars%SetValues(&
               filter(nc)%num_soilp, filter(nc)%soilp, 0._r8, filter(nc)%num_soilc, filter(nc)%soilc, 0._r8)
           if ( use_c13 ) then
             call c13_carbonflux_vars%SetValues(&
                  filter(nc)%num_soilp, filter(nc)%soilp, 0._r8, filter(nc)%num_soilc, filter(nc)%soilc, 0._r8)
           end if
           if ( use_c14 ) then
             call c14_carbonflux_vars%SetValues(&
                  filter(nc)%num_soilp, filter(nc)%soilp, 0._r8, filter(nc)%num_soilc, filter(nc)%soilc, 0._r8)
           end if
           call nitrogenflux_vars%SetValues(&
                  filter(nc)%num_soilp, filter(nc)%soilp, 0._r8, filter(nc)%num_soilc, filter(nc)%soilc, 0._r8)

           call EDbio_vars%SetValues( 0._r8 )
          
         end if  ! end of if-use_ed

         call t_stopf('ecosysdyn')

         ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)
         call t_startf('depvel')
         call depvel_compute(bounds_clump, &
              atm2lnd_vars, canopystate_vars, waterstate_vars, frictionvel_vars, &
              photosyns_vars, drydepvel_vars)
         call t_stopf('depvel')

         if (use_betr)then

           call ep_betr%CalcSmpL(bounds_clump, 1, nlevsoi, filter(nc)%num_soilc, filter(nc)%soilc, &
              temperature_vars%t_soisno_col, soilstate_vars, waterstate_vars, soil_water_retention_curve)

           call ep_betr%SetBiophysForcing(bounds_clump, col, pft,                               &
             carbonflux_vars=carbonflux_vars,                                                &
             waterstate_vars=waterstate_vars,         waterflux_vars=waterflux_vars,         &
             temperature_vars=temperature_vars,       soilhydrology_vars=soilhydrology_vars, &
             atm2lnd_vars=atm2lnd_vars,               canopystate_vars=canopystate_vars,     &
             chemstate_vars=chemstate_vars,           soilstate_vars=soilstate_vars, &
             cnstate_vars = cnstate_vars)

           call ep_betr%StepWithoutDrainage(bounds_clump, col, pft)
         endif
         
         if (use_lch4) then
           !warning: do not call ch4 before CNAnnualUpdate, which will fail the ch4 model
           call t_startf('ch4')
           call ch4 (bounds_clump,                                                                  &
               filter(nc)%num_soilc, filter(nc)%soilc,                                             &
               filter(nc)%num_lakec, filter(nc)%lakec,                                             &
               filter(nc)%num_soilp, filter(nc)%soilp,                                             &
               atm2lnd_vars, lakestate_vars, canopystate_vars, soilstate_vars, soilhydrology_vars, &
               temperature_vars, energyflux_vars, waterstate_vars, waterflux_vars,                 &
               carbonstate_vars, carbonflux_vars, nitrogenflux_vars, ch4_vars, lnd2atm_vars)
           call t_stopf('ch4')
         end if
       endif !end of if is_active_betr_bgc

       ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)
       call t_startf('depvel')
       call depvel_compute(bounds_clump, &
            atm2lnd_vars, canopystate_vars, waterstate_vars, frictionvel_vars, &
            photosyns_vars, drydepvel_vars)
       call t_stopf('depvel')     
       ! ============================================================================
       ! Calculate soil/snow hydrology with drainage (subsurface runoff)
       ! ============================================================================

       call t_startf('hydro2 drainage')

       call HydrologyDrainage(bounds_clump,                   &
            filter(nc)%num_nolakec, filter(nc)%nolakec,       &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
            filter(nc)%num_urbanc, filter(nc)%urbanc,         &                 
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,     &                
            atm2lnd_vars, glc2lnd_vars, temperature_vars,     &
            soilhydrology_vars, soilstate_vars, waterstate_vars, waterflux_vars)

       call t_stopf('hydro2 drainage')     

       if (use_betr) then

          call t_startf('betr drainage')       
          call ep_betr%BeTRSetBiophysForcing(bounds_clump, col, pft, 1, nlevsoi,&
                waterflux_vars=waterflux_vars )
          call ep_betr%StepWithDrainage(bounds_clump, col)
          call t_stopf('betr drainage')

          call t_startf('betr balchk')
          call ep_betr%MassBalanceCheck(bounds_clump)
          call t_stopf('betr balchk')  

!          call bgc_reaction%betr_alm_flux_statevar_feedback(bounds_clump, &
!               filter(nc)%num_soilc, filter(nc)%soilc,                    &
!               carbonstate_vars, nitrogenstate_vars, nitrogenflux_vars,   &
!               tracerstate_vars, tracerflux_vars,  betrtracer_vars)          
       endif          

       ! ============================================================================
       ! Check the energy and water balance, also carbon and nitrogen balance
       ! ============================================================================

       if (.not. use_ed) then

          if (use_cn) then
            
            if (is_active_betr_bgc)then
               !extract nitrogen pool and flux from betr
               !summarize total column nitrogen and carbon
!              call CNFluxStateBetrSummary(bounds_clump, filter(nc)%num_soilc, filter(nc)%soilc,   &
!                filter(nc)%num_soilp, filter(nc)%soilp,                                           &
!                carbonflux_vars, carbonstate_vars,                                                &
!                c13_carbonflux_vars, c13_carbonstate_vars,                                        &
!                c14_carbonflux_vars, c14_carbonstate_vars, nitrogenflux_vars, nitrogenstate_vars)! &
!                betrtracer_vars, tracerflux_vars, tracerstate_vars)                  
            else
             ! FIX(SPM,032414) there are use_ed checks in this routine...be consistent 
             ! (see comment above re: no leaching
               call CNEcosystemDynLeaching(bounds_clump,                &
                  filter(nc)%num_soilc, filter(nc)%soilc,               &
                  filter(nc)%num_soilp, filter(nc)%soilp,               &
                  filter(nc)%num_pcropp, filter(nc)%pcropp, doalb,      &
                  cnstate_vars, carbonflux_vars, carbonstate_vars,      &
                  c13_carbonflux_vars, c13_carbonstate_vars,            &
                  c14_carbonflux_vars, c14_carbonstate_vars, dgvs_vars, &
                  nitrogenflux_vars, nitrogenstate_vars,                &
                  waterstate_vars, waterflux_vars, frictionvel_vars,    &
                  canopystate_vars,                                     &
                  phosphorusflux_vars,phosphorusstate_vars)
             end if

             if (doalb) then   
                call CNVegStructUpdate(filter(nc)%num_soilp, filter(nc)%soilp,   &
                     waterstate_vars, frictionvel_vars, dgvs_vars, cnstate_vars, &
                     carbonstate_vars, canopystate_vars)
            end if
               
          end if
       end if

       call t_startf('balchk')
       call BalanceCheck(bounds_clump, &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c, &
            atm2lnd_vars, glc2lnd_vars, solarabs_vars, waterflux_vars, &
            waterstate_vars, energyflux_vars, canopystate_vars)
       call t_stopf('balchk')

       if (.not. use_ed)then
          if (use_cn) then
             nstep = get_nstep()

             if (nstep < 2 )then
                if (masterproc) then
                   write(iulog,*) '--WARNING-- skipping CN balance check for first timestep'
                end if
             else
                call t_startf('cnbalchk')

                call CBalanceCheck(bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     carbonstate_vars, carbonflux_vars)

                call NBalanceCheck(bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     nitrogenstate_vars, nitrogenflux_vars)


                call PBalanceCheck(bounds_clump, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     phosphorusstate_vars, phosphorusflux_vars)

                call t_stopf('cnbalchk')
             end if
          end if
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
               aerosol_vars, canopystate_vars, waterstate_vars, &
               lakestate_vars, temperature_vars, surfalb_vars)
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
                  waterstate_vars, urbanparams_vars, solarabs_vars, surfalb_vars) 
             call t_stopf('urbsurfalb')
          end if

       end if

    end do
    !$OMP END PARALLEL DO

    ! ============================================================================
    ! Determine gridcell averaged properties to send to atm
    ! ============================================================================

    call t_startf('lnd2atm')
    call lnd2atm(bounds_proc,                                            &
         atm2lnd_vars, surfalb_vars, temperature_vars, frictionvel_vars, &
         waterstate_vars, waterflux_vars, energyflux_vars,               &
         solarabs_vars, carbonflux_vars, drydepvel_vars,                 &
         vocemis_vars, dust_vars, ch4_vars, lnd2atm_vars) 
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
               temperature_vars, waterflux_vars,               &
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
       !
       ! FIX(SPM, 082814) - in the ED branch RF and I comment out the if(.not.
       ! use_ed) then statement ... double check if this is required and why
       !
       !if (.not. use_ed) then
          call t_startf('accum')

          call atm2lnd_vars%UpdateAccVars(bounds_proc)

          call temperature_vars%UpdateAccVars(EDBio_vars, bounds_proc)

          call canopystate_vars%UpdateAccVars(bounds_proc)

          if (use_cndv) then
             call dgvs_vars%UpdateCNDVAccVars(bounds_proc, temperature_vars)
          end if

          if (crop_prog) then
             call crop_vars%UpdateAccVars(bounds_proc, temperature_vars, cnstate_vars)
          end if

          call t_stopf('accum')
       !end if
    end if

    ! ============================================================================
    ! Update history buffer
    ! ============================================================================

    call t_startf('hbuf')
    call hist_update_hbuf(bounds_proc)
    call t_stopf('hbuf')

    ! ============================================================================
    ! Call dv (dynamic vegetation) at last time step of year
    ! NOTE: monp1, dayp1, and secp1 correspond to nstep+1
    ! ============================================================================

    if (use_cndv) then
       call t_startf('d2dgvm')
       dtime = get_step_size()
       call get_curr_date(yrp1, monp1, dayp1, secp1, offset=int(dtime))
       if (monp1==1 .and. dayp1==1 .and. secp1==dtime .and. nstep>0)  then

          ! Get date info.  kyr is used in lpj().  At end of first year, kyr = 2.
          call get_curr_date(yr, mon, day, sec)
          ncdate = yr*10000 + mon*100 + day
          call get_ref_date(yr, mon, day, sec)
          nbdate = yr*10000 + mon*100 + day
          kyr = ncdate/10000 - nbdate/10000 + 1

          if (masterproc) write(iulog,*) 'End of year. CNDV called now: ncdate=', &
               ncdate,' nbdate=',nbdate,' kyr=',kyr,' nstep=', nstep

          nclumps = get_proc_clumps()

          !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
          do nc = 1,nclumps
             call get_clump_bounds(nc, bounds_clump)
             call CNDVDriver(bounds_clump,                                           &
                  filter(nc)%num_natvegp, filter(nc)%natvegp, kyr,                   &
                  atm2lnd_vars, carbonflux_vars, carbonstate_vars, dgvs_vars)
          end do
          !$OMP END PARALLEL DO
       end if
       call t_stopf('d2dgvm')
    end if

    ! ============================================================================
    ! Call ED model on daily timestep
    ! ============================================================================
    
    if ( use_ed ) then
       if ( is_beg_curr_day() ) then ! run ED at the start of each day

          if ( masterproc ) write(iulog,*)  'edtime ed call edmodel ',get_nstep()
          nclumps = get_proc_clumps()

          !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
          do nc = 1,nclumps

             call get_clump_bounds(nc, bounds_clump)

             call edmodel( bounds_clump, & 
                  atm2lnd_vars, soilstate_vars, temperature_vars, waterstate_vars )

             ! link to CLM structures
             call clm_ed_link( bounds_clump, gridCellEdState, &
                  waterstate_vars, canopystate_vars, EDbio_vars, &
                  carbonstate_vars, nitrogenstate_vars, carbonflux_vars) 

             call setFilters( bounds_clump, glc2lnd_vars%icemask_grc )

             !reset surface albedo fluxes in case there is a mismatch between elai and canopy absorbtion. 
             call SurfaceAlbedo(bounds_clump,                      &
                filter_inactive_and_active(nc)%num_nourbanc,       &
                filter_inactive_and_active(nc)%nourbanc,           &
                filter_inactive_and_active(nc)%num_nourbanp,       &
                filter_inactive_and_active(nc)%nourbanp,           &
                filter_inactive_and_active(nc)%num_urbanc,         &
                filter_inactive_and_active(nc)%urbanc,             &
                filter_inactive_and_active(nc)%num_urbanp,         & 
                filter_inactive_and_active(nc)%urbanp,             &
                nextsw_cday, declinp1,                             &
                aerosol_vars, canopystate_vars, waterstate_vars,   &
                lakestate_vars, temperature_vars, surfalb_vars)

          end do
          !$OMP END PARALLEL DO
       end if
    end if ! use_ed branch

    ! ============================================================================
    ! History/Restart output
    ! ============================================================================

    if (.not. use_noio) then

       call t_startf('clm_drv_io')

       ! Create history and write history tapes if appropriate
       call t_startf('clm_drv_io_htapes')

       call hist_htapes_wrapup( rstwr, nlend, bounds_proc,                    &
            soilstate_vars%watsat_col(bounds_proc%begc:bounds_proc%endc, 1:), &
            soilstate_vars%sucsat_col(bounds_proc%begc:bounds_proc%endc, 1:), &
            soilstate_vars%bsw_col(bounds_proc%begc:bounds_proc%endc, 1:),    &
            soilstate_vars%hksat_col(bounds_proc%begc:bounds_proc%endc, 1:))

       !----------------------------------------------
       ! pflotran
       if (use_pflotran) then
            call clm_pf_write_restart(rdate)
       end if
       !----------------------------------------------

       call t_stopf('clm_drv_io_htapes')

       ! Write to CNDV history buffer if appropriate
       if (use_cndv) then
          if (monp1==1 .and. dayp1==1 .and. secp1==dtime .and. nstep>0)  then
             call t_startf('clm_drv_io_hdgvm')
             call CNDVHist( bounds_proc, dgvs_vars )
             if (masterproc) write(iulog,*) 'Annual CNDV calculations are complete'
             call t_stopf('clm_drv_io_hdgvm')
          end if
       end if

       ! Write restart/initial files if appropriate
       if (rstwr) then
          call t_startf('clm_drv_io_wrest')
          filer = restFile_filename(rdate=rdate)

          call restFile_write( bounds_proc, filer,                                            &
               atm2lnd_vars, aerosol_vars, canopystate_vars, cnstate_vars,                    &
               carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, carbonflux_vars, &
               ch4_vars, dgvs_vars, energyflux_vars, frictionvel_vars, lakestate_vars,        &
               nitrogenstate_vars, nitrogenflux_vars, photosyns_vars, soilhydrology_vars,     &
               soilstate_vars, solarabs_vars, surfalb_vars, temperature_vars,                 &
               waterflux_vars, waterstate_vars, EDbio_vars,                                   &
               phosphorusstate_vars,phosphorusflux_vars,                                      &
               ep_betr,                                                                       &
               rdate=rdate)

          call t_stopf('clm_drv_io_wrest')
       end if
       call t_stopf('clm_drv_io')

    end if

  end subroutine clm_drv

  !-----------------------------------------------------------------------
  subroutine clm_drv_init(bounds, &
       num_nolakec, filter_nolakec, &
       num_nolakep, filter_nolakep, &
       num_soilp  , filter_soilp, &
       canopystate_vars, waterstate_vars, waterflux_vars, energyflux_vars)
    !
    ! !DESCRIPTION:
    ! Initialization of clm driver variables needed from previous timestep
    !
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use shr_infnan_mod     , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar         , only : nlevsno
    use clm_varcon         , only : h2osno_max
    use landunit_varcon    , only : istice_mec
    use CanopyStateType    , only : canopystate_type
    use WaterStateType     , only : waterstate_type
    use WaterFluxType      , only : waterflux_type
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
    type(waterstate_type) , intent(inout) :: waterstate_vars
    type(waterflux_type)  , intent(inout) :: waterflux_vars
    type(energyflux_type) , intent(inout) :: energyflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: l, c, p, f, j         ! indices
    integer :: fp, fc                  ! filter indices
    !-----------------------------------------------------------------------

    associate(                                                             & 
         snl                => col%snl                                   , & ! Input:  [integer  (:)   ]  number of snow layers                    
        
         h2osno             => waterstate_vars%h2osno_col                , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                     
         h2osoi_ice         => waterstate_vars%h2osoi_ice_col            , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                      
         h2osoi_liq         => waterstate_vars%h2osoi_liq_col            , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                  
         do_capsnow         => waterstate_vars%do_capsnow_col            , & ! Output: [logical  (:)   ]  true => do snow capping                  
         h2osno_old         => waterstate_vars%h2osno_old_col            , & ! Output: [real(r8) (:)   ]  snow water (mm H2O) at previous time step
         frac_iceold        => waterstate_vars%frac_iceold_col           , & ! Output: [real(r8) (:,:) ]  fraction of ice relative to the tot water

         elai               => canopystate_vars%elai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow    
         esai               => canopystate_vars%esai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow    
         frac_veg_nosno     => canopystate_vars%frac_veg_nosno_patch     , & ! Output: [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         frac_veg_nosno_alb => canopystate_vars%frac_veg_nosno_alb_patch , & ! Output: [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]

         qflx_glcice        => waterflux_vars%qflx_glcice_col            , & ! Output: [real(r8) (:)   ]  flux of new glacier ice (mm H2O/s) [+ = ice grows]

         eflx_bot           => energyflux_vars%eflx_bot_col              , & ! Output: [real(r8) (:)   ]  heat flux from beneath soil/ice column (W/m**2)

         cisun_z            => photosyns_vars%cisun_z_patch              , & ! Output: [real(r8) (:)   ]  intracellular sunlit leaf CO2 (Pa)
         cisha_z            => photosyns_vars%cisha_z_patch                & ! Output: [real(r8) (:)   ]  intracellular shaded leaf CO2 (Pa)
         )

      ! Initialize intracellular CO2 (Pa) parameters each timestep for use in VOCEmission
      do p = bounds%begp,bounds%endp
         cisun_z(p,:) = -999._r8
         cisha_z(p,:) = -999._r8
      end do

      do c = bounds%begc,bounds%endc
         l = col%landunit(c)

         ! Save snow mass at previous time step
         h2osno_old(c) = h2osno(c)

         ! Decide whether to cap snow
         if (h2osno(c) > h2osno_max) then
            do_capsnow(c) = .true.
         else
            do_capsnow(c) = .false.
         end if

         ! Reset flux from beneath soil/ice column 
         eflx_bot(c)  = 0._r8

         ! Initialize qflx_glcice everywhere, to zero.
         qflx_glcice(c) = 0._r8     

      end do

      ! Initialize fraction of vegetation not covered by snow 

      do p = bounds%begp,bounds%endp
         if (pft%active(p)) then
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

  end subroutine clm_drv_init
  
  !-----------------------------------------------------------------------
  subroutine clm_drv_patch2col (bounds, num_nolakec, filter_nolakec, &
       waterstate_vars, energyflux_vars, waterflux_vars)
    !
    ! !DESCRIPTION:
    ! Averages over all patchs for variables defined over both soil and lake
    ! to provide the column-level averages of state and flux variables
    ! defined at the patch level.
    !
    ! !USES:
    use WaterStateType , only : waterstate_type
    use WaterFluxType  , only : waterflux_type
    use EnergyFluxType , only : energyflux_type
    use subgridAveMod  , only : p2c
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds  
    integer               , intent(in)    :: num_nolakec       ! number of column non-lake points in column filter
    integer               , intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
    type(waterstate_type) , intent(inout) :: waterstate_vars
    type(waterflux_type)  , intent(inout) :: waterflux_vars
    type(energyflux_type) , intent(inout) :: energyflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,fc              ! indices
    integer :: num_allc          ! number of active column points
    integer :: filter_allc(bounds%endp-bounds%begp+1)    ! filter for all active column points
    ! -----------------------------------------------------------------

    ! Set up a filter for all active column points

    fc = 0
    do c = bounds%begc,bounds%endc
       if (col%active(c)) then
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
         waterstate_vars%h2ocan_patch(bounds%begp:bounds%endp), &
         waterstate_vars%h2ocan_col(bounds%begc:bounds%endc))

    ! Averaging for patch evaporative flux variables

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_ev_snow_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_ev_snow_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_ev_soil_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_ev_soil_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_ev_h2osfc_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_ev_h2osfc_col(bounds%begc:bounds%endc))

    ! Averaging for patch water flux variables

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_evap_soi_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_evap_soi_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_evap_tot_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_evap_tot_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_rain_grnd_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_rain_grnd_col(bounds%begc:bounds%endc))
    
    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_snow_grnd_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_snow_grnd_col(bounds%begc:bounds%endc))
    
    call p2c (bounds, num_allc, filter_allc, &
         waterflux_vars%qflx_snwcp_liq_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_snwcp_liq_col(bounds%begc:bounds%endc))
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
         waterflux_vars%qflx_snwcp_ice_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_snwcp_ice_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_tran_veg_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_tran_veg_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_evap_grnd_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_evap_grnd_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_allc, filter_allc, &
         waterflux_vars%qflx_evap_soi_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_evap_soi_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_prec_grnd_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_prec_grnd_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_dew_grnd_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_dew_grnd_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_sub_snow_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_sub_snow_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_dew_snow_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_dew_snow_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_irrig_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_irrig_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_tran_veg_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_tran_veg_col(bounds%begc:bounds%endc) )

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_vars%qflx_evap_veg_patch(bounds%begp:bounds%endp), &
         waterflux_vars%qflx_evap_veg_col (bounds%begc:bounds%endc))
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
