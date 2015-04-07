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
  use clm_varctl             , only : use_cn, use_cndv, use_lch4, use_voc, use_noio, use_c13, use_c14
  use clm_time_manager       , only : get_step_size, get_curr_date, get_ref_date, get_nstep, is_beg_curr_day
  use clm_time_manager       , only : get_prev_date
  use clm_varpar             , only : nlevsno, nlevgrnd, crop_prog
  use spmdMod                , only : masterproc, mpicom
  use decompMod              , only : get_proc_clumps, get_clump_bounds, get_proc_bounds, bounds_type
  use filterMod              , only : filter, filter_inactive_and_active
  use filterMod              , only : setExposedvegpFilter
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
  !
  use CNDriverMod            , only : CNDriverNoLeaching, CNDriverLeaching, CNDriverSummary
  use CNVegStructUpdateMod   , only : CNVegStructUpdate 
  use CNAnnualUpdateMod      , only : CNAnnualUpdate
  use CNBalanceCheckMod      , only : BeginCNBalance, CBalanceCheck, NBalanceCheck
  use SoilBiogeochemVerticalProfileMod   , only : SoilBiogeochemVerticalProfile
  use CNFireMod              , only : CNFireInterp
  use CNDVDriverMod          , only : CNDVDriver, CNDVHIST
  use SatellitePhenologyMod  , only : SatellitePhenology, interpMonthlyVeg
  use ndepStreamMod          , only : ndep_interp
  use ActiveLayerMod         , only : alt_calc
  use ch4Mod                 , only : ch4
  use DUSTMod                , only : DustDryDep, DustEmission
  use VOCEmissionMod         , only : VOCEmission
  use EDMainMod              , only : ed_driver
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
  use clm_initializeMod      , only : nutrient_competition_method
  use GridcellType           , only : grc                
  use LandunitType           , only : lun                
  use ColumnType             , only : col                
  use PatchType              , only : patch                
  use clm_instMod
  use clm_initializeMod      , only : soil_water_retention_curve
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
    integer              :: yr_prev                 ! year (0, ...) at start of timestep
    integer              :: mon_prev                ! month (1, ..., 12) at start of timestep
    integer              :: day_prev                ! day of month (1, ..., 31) at start of timestep
    integer              :: sec_prev                ! seconds of the day at start of timestep
    integer              :: ncdate                  ! current date
    integer              :: nbdate                  ! base date (reference date)
    integer              :: kyr                     ! thousand years, equals 2 at end of first year
    character(len=256)   :: filer                   ! restart file name
    integer              :: ier                     ! error code
    type(bounds_type)    :: bounds_clump    
    type(bounds_type)    :: bounds_proc     
    !-----------------------------------------------------------------------

    ! Determine processor bounds and clumps for this processor

    call get_proc_bounds(bounds_proc)
    nclumps = get_proc_clumps()

    ! Update time-related info

    call cnveg_state_inst%CropRestIncYear()

    ! ============================================================================
    ! Specified phenology
    ! ============================================================================

    if (use_cn) then 
       ! For dry-deposition need to call CLMSP so that mlaidiff is obtained
       if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
          call t_startf('interpMonthlyVeg')
          call interpMonthlyVeg(bounds_proc, canopystate_inst)
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
          call interpMonthlyVeg(bounds_proc, canopystate_inst)
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
            temperature_inst, canopystate_inst) 

       if (use_cn) then
          !  Note (WJS, 6-12-13): Because of this routine's placement in the driver sequence
          !  (it is called very early in each timestep, before weights are adjusted and
          !  filters are updated), it may be necessary for this routine to compute values over
          !  inactive as well as active points (since some inactive points may soon become
          !  active) - so that's what is done now. Currently, it seems to be okay to do this,
          !  because the variables computed here seem to only depend on quantities that are
          !  valid over inactive as well as active points.

          call SoilBiogeochemVerticalProfile(bounds_clump                                       , &
               filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc   , &
               filter_inactive_and_active(nc)%num_soilp, filter_inactive_and_active(nc)%soilp   , &
	       canopystate_inst, soilstate_inst, soilbiogeochem_state_inst)
       end if

       call t_stopf("decomp_vert")
    end do
    !$OMP END PARALLEL DO

    ! ============================================================================
    ! Initialize the mass balance checks for carbon and nitrogen, and zero fluxes for
    ! transient land cover
    ! ============================================================================

    if (use_cn) then
       !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
       do nc = 1,nclumps
          call get_clump_bounds(nc, bounds_clump)

          call t_startf('begcnbal')

          call BeginCNBalance(bounds_clump, filter(nc)%num_soilc, filter(nc)%soilc, &
               cnveg_carbonstate_inst, cnveg_nitrogenstate_inst)

          call cnveg_carbonflux_inst%ZeroDWT(bounds_clump)
          if (use_c13) then
             call c13_cnveg_carbonflux_inst%ZeroDWT(bounds_clump)
          end if
          if (use_c14) then
             call c14_cnveg_carbonflux_inst%ZeroDWT(bounds_clump)
          end if
          call cnveg_nitrogenflux_inst%ZeroDWT(bounds_clump)
          call cnveg_carbonstate_inst%ZeroDWT(bounds_clump)
          call cnveg_nitrogenstate_inst%ZeroDWT(bounds_clump)

          call soilbiogeochem_carbonflux_inst%ZeroDWT(bounds_clump)
          if (use_c13) then
             call c13_soilbiogeochem_carbonflux_inst%ZeroDWT(bounds_clump)
          end if
          if (use_c14) then
             call c14_soilbiogeochem_carbonflux_inst%ZeroDWT(bounds_clump)
          end if

          call t_stopf('begcnbal')
       end do
       !$OMP END PARALLEL DO
    end if

    ! ============================================================================
    ! Update subgrid weights with dynamic landcover (prescribed transient patches,
    ! CNDV, and or dynamic landunits), and do related adjustments. Note that this
    ! call needs to happen outside loops over nclumps.
    ! ============================================================================

    call t_startf('dyn_subgrid')
    call dynSubgrid_driver(bounds_proc,                                                  &
         urbanparams_inst, soilstate_inst, soilhydrology_inst, lakestate_inst,           &
         waterstate_inst, waterflux_inst, temperature_inst, energyflux_inst,             &
         canopystate_inst, photosyns_inst, dgvs_inst, glc2lnd_inst, cnveg_state_inst,    &
         cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
         cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst,    &
         cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,                              &
         soilbiogeochem_state_inst, soilbiogeochem_carbonflux_inst)
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
            soilhydrology_inst, waterstate_inst)
       call t_stopf('begwbal')
    end do
    !$OMP END PARALLEL DO

    ! ============================================================================
    ! Update dynamic N deposition field, on albedo timestep
    ! currently being done outside clumps loop, but no reason why it couldn't be
    ! re-written to go inside.
    ! ============================================================================

    if (use_cn) then
       call t_startf('ndep_interp')
       call ndep_interp(bounds_proc, atm2lnd_inst)
       call CNFireInterp(bounds_proc)
       call t_stopf('ndep_interp')
    end if

    ! ============================================================================
    ! Initialize variables from previous time step, downscale atm forcings, and
    ! Determine canopy interception and precipitation onto ground surface.
    ! Determine the fraction of foliage covered by water and the fraction
    ! of foliage that is dry and transpiring. Initialize snow layer if the
    ! snow accumulation exceeds 10 mm.
    ! ============================================================================

    ! Get time as of beginning of time step
    call get_prev_date(yr_prev, mon_prev, day_prev, sec_prev)

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
            canopystate_inst, waterstate_inst, waterflux_inst, energyflux_inst)

       call downscale_forcings(bounds_clump, &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c, &
            atm2lnd_inst)

       call t_stopf('drvinit')

       ! Update filters that depend on variables set in clm_drv_init
       
       call setExposedvegpFilter(bounds_clump, &
            canopystate_inst%frac_veg_nosno_patch(bounds_clump%begp:bounds_clump%endp))

       ! Irrigation flux

       call irrigation_inst%ApplyIrrigation(bounds_clump)

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
            atm2lnd_inst, canopystate_inst, temperature_inst, &
            aerosol_inst, waterstate_inst, waterflux_inst, &
            irrigation_inst)
       call t_stopf('canhydro')

       ! ============================================================================
       ! Surface Radiation
       ! ============================================================================

       call t_startf('surfrad')

       ! Surface Radiation primarily for non-urban columns 

       call SurfaceRadiation(bounds_clump,                           &
            filter(nc)%num_nourbanp, filter(nc)%nourbanp,            &
            filter(nc)%num_urbanp, filter(nc)%urbanp,                &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                &
            ed_allsites_inst(bounds_clump%begg:bounds_clump%endg), atm2lnd_inst, &
            waterstate_inst, canopystate_inst, surfalb_inst, solarabs_inst, surfrad_inst)

       ! Surface Radiation for only urban columns

       call UrbanRadiation(bounds_clump,                                       &
            filter(nc)%num_nourbanl, filter(nc)%nourbanl,                      &
            filter(nc)%num_urbanl, filter(nc)%urbanl,                          &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                          &
            filter(nc)%num_urbanp, filter(nc)%urbanp,                          &
            atm2lnd_inst, waterstate_inst, temperature_inst, urbanparams_inst, &
            solarabs_inst, surfalb_inst, energyflux_inst)

       call t_stopf('surfrad')

       ! ============================================================================
       ! Determine leaf temperature and surface fluxes based on ground
       ! temperature from previous time step.
       ! ============================================================================

       call t_startf('bgp1')
       call CanopyTemperature(bounds_clump,                                   &
            filter(nc)%num_nolakec, filter(nc)%nolakec,                       &
            filter(nc)%num_nolakep, filter(nc)%nolakep,                       &
            atm2lnd_inst, canopystate_inst, soilstate_inst, frictionvel_inst, &
            waterstate_inst, waterflux_inst, energyflux_inst, temperature_inst)
       call t_stopf('bgp1')

       ! ============================================================================
       ! Determine fluxes
       ! ============================================================================

       call t_startf('bgflux')

       ! Bareground fluxes for all patches except lakes and urban landunits

       call BareGroundFluxes(bounds_clump,                                 &
            filter(nc)%num_noexposedvegp, filter(nc)%noexposedvegp,          &
            atm2lnd_inst, soilstate_inst,                &
            frictionvel_inst, ch4_inst, energyflux_inst, temperature_inst, &
            waterflux_inst, waterstate_inst, photosyns_inst, humanindex_inst)
       call t_stopf('bgflux')

       ! non-bareground fluxes for all patches except lakes and urban landunits
       ! Calculate canopy temperature, latent and sensible fluxes from the canopy,
       ! and leaf water change by evapotranspiration

       call t_startf('canflux')
       call CanopyFluxes(bounds_clump,                                                      &
            filter(nc)%num_exposedvegp, filter(nc)%exposedvegp,                             &
            ed_allsites_inst(bounds_clump%begg:bounds_clump%endg),                          &
            atm2lnd_inst, canopystate_inst, cnveg_state_inst,                               &
            energyflux_inst, frictionvel_inst, soilstate_inst, solarabs_inst, surfalb_inst, &
            temperature_inst, waterflux_inst, waterstate_inst, ch4_inst, photosyns_inst,    &
            humanindex_inst, soil_water_retention_curve) 
       call t_stopf('canflux')

       ! Fluxes for all urban landunits

       call t_startf('uflux')
       call UrbanFluxes(bounds_clump,                                           &
            filter(nc)%num_nourbanl, filter(nc)%nourbanl,                       &
            filter(nc)%num_urbanl, filter(nc)%urbanl,                           &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                           &
            filter(nc)%num_urbanp, filter(nc)%urbanp,                           &
            atm2lnd_inst, urbanparams_inst, soilstate_inst, temperature_inst,   &
            waterstate_inst, frictionvel_inst, energyflux_inst, waterflux_inst, &
            humanindex_inst)
       call t_stopf('uflux')

       ! Fluxes for all lake landunits

       call t_startf('bgplake')
       call LakeFluxes(bounds_clump,                                         &
            filter(nc)%num_lakec, filter(nc)%lakec,                          &
            filter(nc)%num_lakep, filter(nc)%lakep,                          &
            atm2lnd_inst, solarabs_inst, frictionvel_inst, temperature_inst, &
            energyflux_inst, waterstate_inst, waterflux_inst, lakestate_inst,&  
            humanindex_inst) 

       ! ============================================================================
       ! Determine irrigation needed for future time steps
       ! ============================================================================

       ! This needs to be called after btran is computed

       call irrigation_inst%CalcIrrigationNeeded( &
            bounds             = bounds_clump, &
            num_exposedvegp    = filter(nc)%num_exposedvegp, &
            filter_exposedvegp = filter(nc)%exposedvegp, &
            time_prev          = sec_prev, &
            elai               = canopystate_inst%elai_patch(bounds_clump%begp:bounds_clump%endp), &
            btran              = energyflux_inst%btran_patch(bounds_clump%begp:bounds_clump%endp), &
            rootfr             = soilstate_inst%rootfr_patch(bounds_clump%begp:bounds_clump%endp    , 1:nlevgrnd), &
            t_soisno           = temperature_inst%t_soisno_col(bounds_clump%begc:bounds_clump%endc  , 1:nlevgrnd), &
            eff_porosity       = soilstate_inst%eff_porosity_col(bounds_clump%begc:bounds_clump%endc, 1:nlevgrnd), &
            h2osoi_liq         = waterstate_inst%h2osoi_liq_col(bounds_clump%begc:bounds_clump%endc , 1:nlevgrnd))

       ! ============================================================================
       ! DUST and VOC emissions
       ! ============================================================================

       call t_startf('bgc')

       ! Dust mobilization (C. Zender's modified codes)
       call DustEmission(bounds_clump,                                       &
            filter(nc)%num_nolakep, filter(nc)%nolakep,                      &
            atm2lnd_inst, soilstate_inst, canopystate_inst, waterstate_inst, &
            frictionvel_inst, dust_inst)

       ! Dust dry deposition (C. Zender's modified codes)
       call DustDryDep(bounds_clump, &
            atm2lnd_inst, frictionvel_inst, dust_inst)

       ! VOC emission (A. Guenther's MEGAN (2006) model)
       if (use_voc) then
          call VOCEmission(bounds_clump,                                         &
               filter(nc)%num_soilp, filter(nc)%soilp,                           &
               atm2lnd_inst, canopystate_inst, photosyns_inst, temperature_inst, &
               vocemis_inst)
       end if

       call t_stopf('bgc')

       ! ============================================================================
       ! Determine temperatures
       ! ============================================================================

       ! Set lake temperature 

       call LakeTemperature(bounds_clump,                                             &
            filter(nc)%num_lakec, filter(nc)%lakec,                                   &
            filter(nc)%num_lakep, filter(nc)%lakep,                                   & 
            solarabs_inst, soilstate_inst, waterstate_inst, waterflux_inst, ch4_inst, &
            energyflux_inst, temperature_inst, lakestate_inst)
       call t_stopf('bgplake')

       ! Set soil/snow temperatures including ground temperature

       call t_startf('soiltemperature')
       call SoilTemperature(bounds_clump,                                                      &
            filter(nc)%num_urbanl  , filter(nc)%urbanl,                                        &
            filter(nc)%num_nolakec , filter(nc)%nolakec,                                       &
            atm2lnd_inst, urbanparams_inst, canopystate_inst, waterstate_inst, waterflux_inst, &
            solarabs_inst, soilstate_inst, energyflux_inst,  temperature_inst)
       call t_stopf('soiltemperature')

       ! ============================================================================
       ! update surface fluxes for new ground temperature.
       ! ============================================================================

       call t_startf('bgp2')
       call SoilFluxes(bounds_clump,                                                          &
            filter(nc)%num_urbanl,  filter(nc)%urbanl,                                        &
            filter(nc)%num_nolakec, filter(nc)%nolakec,                                       &
            filter(nc)%num_nolakep, filter(nc)%nolakep,                                       &
            atm2lnd_inst, solarabs_inst, temperature_inst, canopystate_inst, waterstate_inst, &
            energyflux_inst, waterflux_inst)            
       call t_stopf('bgp2')

       ! ============================================================================
       ! Perform averaging from patch level to column level
       ! ============================================================================

       call t_startf('patch2col')
       call clm_drv_patch2col(bounds_clump, filter(nc)%num_nolakec, filter(nc)%nolakec, &
            waterstate_inst, energyflux_inst, waterflux_inst)
       call t_stopf('patch2col')

       ! ============================================================================
       ! Vertical (column) soil and surface hydrology
       ! ============================================================================

       ! Note that filter_snowc and filter_nosnowc are returned by
       ! LakeHydrology after the new snow filter is built

       call t_startf('hydro without drainage')

       call HydrologyNoDrainage(bounds_clump,                                  &
            filter(nc)%num_nolakec, filter(nc)%nolakec,                        &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc,                  &
            filter(nc)%num_urbanc, filter(nc)%urbanc,                          &
            filter(nc)%num_snowc, filter(nc)%snowc,                            &
            filter(nc)%num_nosnowc, filter(nc)%nosnowc,                        &
            atm2lnd_inst, soilstate_inst, energyflux_inst, temperature_inst,   &
            waterflux_inst, waterstate_inst, soilhydrology_inst, aerosol_inst, &
            soil_water_retention_curve)

       !  Calculate column-integrated aerosol masses, and
       !  mass concentrations for radiative calculations and output
       !  (based on new snow level state, after SnowFilter is rebuilt.
       !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there 
       !  can be zero snow layers but an active column in filter)
      
       call AerosolMasses( bounds_clump,                                   &
            num_on=filter(nc)%num_snowc, filter_on=filter(nc)%snowc,       &
            num_off=filter(nc)%num_nosnowc, filter_off=filter(nc)%nosnowc, &
            waterflux_inst=waterflux_inst,                                 &
            waterstate_inst=waterstate_inst,                               &
            aerosol_inst=aerosol_inst)                      

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
            atm2lnd_inst, temperature_inst, soilstate_inst, waterstate_inst, waterflux_inst, &
            energyflux_inst, aerosol_inst, lakestate_inst)
       
       !  Calculate column-integrated aerosol masses, and
       !  mass concentrations for radiative calculations and output
       !  (based on new snow level state, after SnowFilter is rebuilt.
       !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there 
       !  can be zero snow layers but an active column in filter)

       call AerosolMasses(bounds_clump,                                            &
            num_on=filter(nc)%num_lakesnowc, filter_on=filter(nc)%lakesnowc,       &
            num_off=filter(nc)%num_lakenosnowc, filter_off=filter(nc)%lakenosnowc, &
            waterflux_inst=waterflux_inst,                                         &
            waterstate_inst=waterstate_inst,                                       &
            aerosol_inst=aerosol_inst)                      

       ! Must be done here because must use a snow filter for lake columns

       call SnowAge_grain(bounds_clump,                         &
            filter(nc)%num_lakesnowc, filter(nc)%lakesnowc,     &
            filter(nc)%num_lakenosnowc, filter(nc)%lakenosnowc, &
            waterflux_inst, waterstate_inst, temperature_inst)

       call t_stopf('hylake')

       ! ============================================================================
       ! ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)
       ! ============================================================================

       do c = bounds_clump%begc,bounds_clump%endc
          l = col%landunit(c)
          if (lun%urbpoi(l)) then
             ! Urban landunit use Bonan 1996 (LSM Technical Note)
             waterstate_inst%frac_sno_col(c) = min( waterstate_inst%snow_depth_col(c)/0.05_r8, 1._r8)
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
            waterflux_inst, waterstate_inst, temperature_inst)
       call t_stopf('snow_init')

       ! ============================================================================
       ! Ecosystem dynamics: Uses CN, CNDV, or static parameterizations
       ! ============================================================================

       ! FIX(SPM,032414)  push these checks into the routines below and/or make this consistent.

       ! fully prognostic canopy structure and C-N biogeochemistry
       ! - CNDV defined: prognostic biogeography; else prescribed
       ! - crop model:  crop algorithms called from within CNDriver

       if (use_cn) then 
          call t_startf('ecosysdyn')
          call CNDriverNoLeaching(bounds_clump,                                         &
               filter(nc)%num_soilc, filter(nc)%soilc,                                  &
               filter(nc)%num_soilp, filter(nc)%soilp,                                  &
               filter(nc)%num_pcropp, filter(nc)%pcropp, doalb,                         &
               cnveg_state_inst,                                                        &
               cnveg_carbonflux_inst, cnveg_carbonstate_inst,                           &
               c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst,                   &
               c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst,                   &
               cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst,                       &
               soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst,         &
               c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
               c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
               soilbiogeochem_state_inst,                                               &
               soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst,     &
               atm2lnd_inst, waterstate_inst, waterflux_inst,                           &
               canopystate_inst, soilstate_inst, temperature_inst, crop_inst, ch4_inst, &
               dgvs_inst, photosyns_inst, soilhydrology_inst, energyflux_inst,          &
               nutrient_competition_method)

          call CNAnnualUpdate(bounds_clump,            &
               filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_soilp, filter(nc)%soilp, &
               cnveg_state_inst, cnveg_carbonflux_inst)
          call t_stopf('ecosysdyn')

       end if

       ! Prescribed biogeography - prescribed canopy structure, some prognostic carbon fluxes

       if ((.not. use_cn) .and. (.not. use_ed) .and. (doalb)) then 
          call SatellitePhenology(bounds_clump, filter(nc)%num_nolakep, filter(nc)%nolakep, &
               waterstate_inst, canopystate_inst)
       end if  
          
       ! Ecosystem demography

       if (use_ed) then
          call ed_clm_inst%SetValues( bounds_clump, 0._r8 )
       end if  

       ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)

       call t_startf('depvel')
       call depvel_compute(bounds_clump, &
            atm2lnd_inst, canopystate_inst, waterstate_inst, frictionvel_inst, &
            photosyns_inst, drydepvel_inst)
       call t_stopf('depvel')

       ! Calculation of methane fluxes

       if (use_lch4) then
          call t_startf('ch4')
          call ch4 (bounds_clump,                                                                  &
               filter(nc)%num_soilc, filter(nc)%soilc,                                             &
               filter(nc)%num_lakec, filter(nc)%lakec,                                             &
               filter(nc)%num_soilp, filter(nc)%soilp,                                             &
               atm2lnd_inst, lakestate_inst, canopystate_inst, soilstate_inst, soilhydrology_inst, &
               temperature_inst, energyflux_inst, waterstate_inst, waterflux_inst,                 &
               cnveg_carbonflux_inst, soilbiogeochem_carbonflux_inst,                          &
               soilbiogeochem_nitrogenflux_inst, ch4_inst, lnd2atm_inst)
          call t_stopf('ch4')
       end if

       ! ============================================================================
       ! Calculate soil/snow hydrology with drainage (subsurface runoff)
       ! ============================================================================

       call t_startf('hydro2 drainage')

       call HydrologyDrainage(bounds_clump,                   &
            filter(nc)%num_nolakec, filter(nc)%nolakec,       &
            filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
            filter(nc)%num_urbanc, filter(nc)%urbanc,         &                 
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,     &                
            atm2lnd_inst, glc2lnd_inst, temperature_inst,     &
            soilhydrology_inst, soilstate_inst, waterstate_inst, waterflux_inst, &
            irrigation_inst)

       call t_stopf('hydro2 drainage')     

       ! ============================================================================
       ! - Update the nitrogen leaching rate as a function of soluble mineral N 
       !   and total soil water outflow.
       ! - Call to all CN summary routines
       ! - On the radiation time step, use C state variables to diagnose
       !   vegetation structure (LAI, SAI, height)
       ! ============================================================================

       if (use_cn) then

          ! Update the nitrogen leaching rate as a function of soluble mineral N 
          ! and total soil water outflow.

          call CNDriverLeaching(bounds_clump,                     &
               filter(nc)%num_soilc, filter(nc)%soilc,            &
               filter(nc)%num_soilp, filter(nc)%soilp,            &
               waterstate_inst, waterflux_inst,                   &
               cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, &
               SoilBiogeochem_nitrogenflux_inst, SoilBiogeochem_nitrogenstate_inst)

          ! Call to all CN summary routines
          
          call  CNDriverSummary(bounds_clump,                                           &
               filter(nc)%num_soilc, filter(nc)%soilc,                                  &
               filter(nc)%num_soilp, filter(nc)%soilp,                                  &
               cnveg_state_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
               c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, &
               c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, &
               cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, &
               soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
               c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
               c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
               soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)

          ! On the radiation time step, use C state variables to calculate
          ! vegetation structure (LAI, SAI, height)

          if (doalb) then   
             call CNVegStructUpdate(filter(nc)%num_soilp, filter(nc)%soilp, &
                  waterstate_inst, frictionvel_inst, dgvs_inst, cnveg_state_inst, &
                  cnveg_carbonstate_inst, canopystate_inst)
          end if

       end if

       ! ============================================================================
       ! Check the energy and water balance and also carbon and nitrogen balance
       ! ============================================================================

       call t_startf('balchk')
       call BalanceCheck(bounds_clump, &
            filter(nc)%num_do_smb_c, filter(nc)%do_smb_c, &
            atm2lnd_inst, glc2lnd_inst, solarabs_inst, waterflux_inst, &
            waterstate_inst, irrigation_inst, energyflux_inst, canopystate_inst)
       call t_stopf('balchk')

       ! ============================================================================
       ! Check the carbon and nitrogen balance
       ! ============================================================================

       if (use_cn) then
          nstep = get_nstep()
          if (nstep < 2 )then
             if (masterproc) then
                write(iulog,*) '--WARNING-- skipping CN balance check for first timestep'
             end if
          else
             call t_startf('cnbalchk')

             call CBalanceCheck(bounds_clump, filter(nc)%num_soilc, filter(nc)%soilc, &
                  soilbiogeochem_carbonflux_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst)

             call NBalanceCheck(bounds_clump, filter(nc)%num_soilc, filter(nc)%soilc, &
                  soilbiogeochem_nitrogenflux_inst, cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst)

             call t_stopf('cnbalchk')
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
               ed_allsites_inst(bounds_clump%begg:bounds_clump%endg), &
               aerosol_inst, canopystate_inst, waterstate_inst, &
               lakestate_inst, temperature_inst, surfalb_inst)
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
                  waterstate_inst, urbanparams_inst,         &
                  solarabs_inst, surfalb_inst)
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
         atm2lnd_inst, surfalb_inst, temperature_inst, frictionvel_inst, &
         waterstate_inst, waterflux_inst, energyflux_inst,               &
         solarabs_inst, cnveg_carbonflux_inst, drydepvel_inst,       &
         vocemis_inst, dust_inst, ch4_inst, lnd2atm_inst) 
    call t_stopf('lnd2atm')

    ! ============================================================================
    ! Determine gridcell averaged properties to send to glc
    ! ============================================================================

    if (create_glacier_mec_landunit) then
       call t_startf('lnd2glc')
       !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
       do nc = 1,nclumps
          call get_clump_bounds(nc, bounds_clump)
          call lnd2glc_inst%update_lnd2glc(bounds_clump,       &
               filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,   &
               temperature_inst, waterflux_inst,               &
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
    call write_diagnostic(bounds_proc, wrtdia, nstep, lnd2atm_inst)
    call t_stopf('wrtdiag')

    ! ============================================================================
    ! Update accumulators
    ! ============================================================================

    ! FIX(SPM,032414) double check why this isn't called for ED
    ! FIX(SPM, 082814) - in the ED branch RF and I commented out the if(.not.
    ! use_ed) then statement ... double check if this is required and why

    if (nstep > 0) then
       call t_startf('accum')

       call atm2lnd_inst%UpdateAccVars(bounds_proc)

       call temperature_inst%UpdateAccVars(bounds_proc)

       call canopystate_inst%UpdateAccVars(bounds_proc)

       if (use_ed) then
          call ed_phenology_inst%accumulateAndExtract(bounds_proc, &
               temperature_inst%t_ref2m_patch(bounds_proc%begp:bounds_proc%endp), &
               patch%gridcell(bounds_proc%begp:bounds_proc%endp), &
               grc%latdeg(bounds_proc%begg:bounds_proc%endg), &
               mon, day, sec)
       endif

       if (use_cndv) then
          call dgvs_inst%UpdateAccVars(bounds_proc, &
               t_a10_patch=temperature_inst%t_a10_patch(bounds_proc%begp:bounds_proc%endp), &
               t_ref2m_patch=temperature_inst%t_ref2m_patch(bounds_proc%begp:bounds_proc%endp))
       end if

       if (crop_prog) then
          call crop_inst%CropUpdateAccVars(bounds_proc, &
               temperature_inst%t_ref2m_patch, temperature_inst%t_soisno_col, cnveg_state_inst)
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
                  atm2lnd_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, dgvs_inst)
          end do
          !$OMP END PARALLEL DO
       end if
       call t_stopf('d2dgvm')
    end if

    ! ============================================================================
    ! Call ED model on daily timestep
    ! ============================================================================
    
    if ( use_ed  .and. is_beg_curr_day() ) then ! run ED at the start of each day

       if ( masterproc ) then
          write(iulog,*)  'edtime ed call edmodel ',get_nstep()
       end if

       !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
       do nc = 1, nclumps

          call get_clump_bounds(nc, bounds_clump)

          call ed_driver( bounds_clump,                               &
               ed_allsites_inst(bounds_clump%begg:bounds_clump%endg), &
               ed_clm_inst, ed_phenology_inst,                        &
               atm2lnd_inst, soilstate_inst, temperature_inst,        &
               waterstate_inst, canopystate_inst)

          call setFilters( bounds_clump, glc2lnd_inst%icemask_grc )

          !reset surface albedo fluxes in case there is a mismatch between elai and canopy absorbtion. 
          call SurfaceAlbedo(bounds_clump,                        &
               filter_inactive_and_active(nc)%num_nourbanc,       &
               filter_inactive_and_active(nc)%nourbanc,           &
               filter_inactive_and_active(nc)%num_nourbanp,       &
               filter_inactive_and_active(nc)%nourbanp,           &
               filter_inactive_and_active(nc)%num_urbanc,         &
               filter_inactive_and_active(nc)%urbanc,             &
               filter_inactive_and_active(nc)%num_urbanp,         & 
               filter_inactive_and_active(nc)%urbanp,             &
               nextsw_cday, declinp1,                             &
               ed_allsites_inst(bounds_clump%begg:bounds_clump%endg), &
               aerosol_inst, canopystate_inst, waterstate_inst, &
               lakestate_inst, temperature_inst, surfalb_inst)

       end do
       !$OMP END PARALLEL DO

    end if ! use_ed branch

    ! ============================================================================
    ! History/Restart output
    ! ============================================================================

    if (.not. use_noio) then

       call t_startf('clm_drv_io')

       ! Create history and write history tapes if appropriate
       call t_startf('clm_drv_io_htapes')

       call hist_htapes_wrapup( rstwr, nlend, bounds_proc,                    &
            soilstate_inst%watsat_col(bounds_proc%begc:bounds_proc%endc, 1:), &
            soilstate_inst%sucsat_col(bounds_proc%begc:bounds_proc%endc, 1:), &
            soilstate_inst%bsw_col(bounds_proc%begc:bounds_proc%endc, 1:),    &
            soilstate_inst%hksat_col(bounds_proc%begc:bounds_proc%endc, 1:))

       call t_stopf('clm_drv_io_htapes')

       ! Write to CNDV history buffer if appropriate
       if (use_cndv) then
          if (monp1==1 .and. dayp1==1 .and. secp1==dtime .and. nstep>0)  then
             call t_startf('clm_drv_io_hdgvm')
             call CNDVHist( bounds_proc, dgvs_inst )
             if (masterproc) write(iulog,*) 'Annual CNDV calculations are complete'
             call t_stopf('clm_drv_io_hdgvm')
          end if
       end if

       ! Write restart/initial files if appropriate
       if (rstwr) then
          call t_startf('clm_drv_io_wrest')
          filer = restFile_filename(rdate=rdate)

          call restFile_write( bounds_proc, filer, rdate=rdate )

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
       canopystate_inst, waterstate_inst, waterflux_inst, energyflux_inst)
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
    type(canopystate_type), intent(inout) :: canopystate_inst
    type(waterstate_type) , intent(inout) :: waterstate_inst
    type(waterflux_type)  , intent(inout) :: waterflux_inst
    type(energyflux_type) , intent(inout) :: energyflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: l, c, p, f, j         ! indices
    integer :: fp, fc                  ! filter indices
    !-----------------------------------------------------------------------

    associate(                                                             & 
         snl                => col%snl                                   , & ! Input:  [integer  (:)   ]  number of snow layers                    
        
         h2osno             => waterstate_inst%h2osno_col                , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                     
         h2osoi_ice         => waterstate_inst%h2osoi_ice_col            , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                      
         h2osoi_liq         => waterstate_inst%h2osoi_liq_col            , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                  
         do_capsnow         => waterstate_inst%do_capsnow_col            , & ! Output: [logical  (:)   ]  true => do snow capping                  
         h2osno_old         => waterstate_inst%h2osno_old_col            , & ! Output: [real(r8) (:)   ]  snow water (mm H2O) at previous time step
         frac_iceold        => waterstate_inst%frac_iceold_col           , & ! Output: [real(r8) (:,:) ]  fraction of ice relative to the tot water

         elai               => canopystate_inst%elai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow    
         esai               => canopystate_inst%esai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow    
         frac_veg_nosno     => canopystate_inst%frac_veg_nosno_patch     , & ! Output: [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         frac_veg_nosno_alb => canopystate_inst%frac_veg_nosno_alb_patch , & ! Output: [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]

         qflx_glcice        => waterflux_inst%qflx_glcice_col            , & ! Output: [real(r8) (:)   ]  flux of new glacier ice (mm H2O/s) [+ = ice grows]

         eflx_bot           => energyflux_inst%eflx_bot_col              , & ! Output: [real(r8) (:)   ]  heat flux from beneath soil/ice column (W/m**2)

         cisun_z            => photosyns_inst%cisun_z_patch              , & ! Output: [real(r8) (:)   ]  intracellular sunlit leaf CO2 (Pa)
         cisha_z            => photosyns_inst%cisha_z_patch                & ! Output: [real(r8) (:)   ]  intracellular shaded leaf CO2 (Pa)
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
         if (patch%active(p)) then
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
       waterstate_inst, energyflux_inst, waterflux_inst)
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
    type(waterstate_type) , intent(inout) :: waterstate_inst
    type(waterflux_type)  , intent(inout) :: waterflux_inst
    type(energyflux_type) , intent(inout) :: energyflux_inst
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
         waterstate_inst%h2ocan_patch(bounds%begp:bounds%endp), &
         waterstate_inst%h2ocan_col(bounds%begc:bounds%endc))

    ! Averaging for patch evaporative flux variables

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_inst%qflx_ev_snow_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_ev_snow_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_inst%qflx_ev_soil_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_ev_soil_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_inst%qflx_ev_h2osfc_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_ev_h2osfc_col(bounds%begc:bounds%endc))

    ! Averaging for patch water flux variables

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_inst%qflx_evap_soi_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_evap_soi_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_inst%qflx_evap_tot_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_evap_tot_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_inst%qflx_rain_grnd_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_rain_grnd_col(bounds%begc:bounds%endc))
    
    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_inst%qflx_snow_grnd_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_snow_grnd_col(bounds%begc:bounds%endc))
    
    call p2c (bounds, num_allc, filter_allc, &
         waterflux_inst%qflx_snwcp_liq_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_snwcp_liq_col(bounds%begc:bounds%endc))
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
         waterflux_inst%qflx_snwcp_ice_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_snwcp_ice_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_inst%qflx_tran_veg_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_tran_veg_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_inst%qflx_evap_grnd_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_evap_grnd_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_allc, filter_allc, &
         waterflux_inst%qflx_evap_soi_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_evap_soi_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_inst%qflx_prec_grnd_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_prec_grnd_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_inst%qflx_dew_grnd_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_dew_grnd_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_inst%qflx_sub_snow_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_sub_snow_col(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         waterflux_inst%qflx_dew_snow_patch(bounds%begp:bounds%endp), &
         waterflux_inst%qflx_dew_snow_col(bounds%begc:bounds%endc))

  end subroutine clm_drv_patch2col

  !------------------------------------------------------------------------
  subroutine write_diagnostic (bounds, wrtdia, nstep, lnd2atm_inst)
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
    type(lnd2atm_type) , intent(in) :: lnd2atm_inst
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
    !------------------------------------------------------------------------

    call get_proc_global(ng=numg)

    if (wrtdia) then

       call t_barrierf('sync_write_diag', mpicom)
       psum = sum(lnd2atm_inst%t_rad_grc(bounds%begg:bounds%endg))
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

       if (masterproc) then
          write(iulog,*)'clm: completed timestep ',nstep
          call shr_sys_flush(iulog)
       end if

    endif

1000 format (1x,'nstep = ',i10,'   TS = ',f21.15)

  end subroutine write_diagnostic

end module clm_driver
