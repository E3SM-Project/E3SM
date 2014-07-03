module clm_driver

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides the main CLM driver physics calling sequence.  Most
  ! computations occurs over ``clumps'' of gridcells (and associated subgrid
  ! scale entities) assigned to each MPI process. Computation is further
  ! parallelized by looping over clumps on each process using shared memory OpenMP.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clmtype
  use clm_varctl          , only : wrtdia, iulog, create_glacier_mec_landunit
  use clm_varctl          , only : use_cn, use_cndv, use_lch4, use_voc, use_noio
  use spmdMod             , only : masterproc,mpicom
  use decompMod           , only : get_proc_clumps, get_clump_bounds, get_proc_bounds, bounds_type
  use filterMod           , only : filter, filter_inactive_and_active
  use dynSubgridDriverMod , only : dynSubgrid_driver
  use CNDVMod             , only : dv, histCNDV
  use clm_varcon          , only : zlnd
  use clm_time_manager    , only : get_step_size,get_curr_date,get_ref_date,get_nstep
  use CropRestMod         , only : CropRestIncYear
  use histFileMod         , only : hist_update_hbuf, hist_htapes_wrapup
  use restFileMod         , only : restFile_write, restFile_filename
  use accFldsMod          , only : updateAccFlds
  use clm_driverInitMod   , only : clm_driverInit
  use BalanceCheckMod     , only : BeginWaterBalance, BalanceCheck
  use SurfaceRadiationMod , only : SurfaceRadiation
  use Hydrology1Mod       , only : Hydrology1
  use Hydrology2Mod       , only : Hydrology2NoDrainage, Hydrology2Drainage
  use SLakeFluxesMod      , only : SLakeFluxes
  use SLakeTemperatureMod , only : SLakeTemperature
  use SLakeHydrologyMod   , only : SLakeHydrology
  use Biogeophysics1Mod   , only : Biogeophysics1
  use BareGroundFluxesMod , only : BareGroundFluxes
  use CanopyFluxesMod     , only : CanopyFluxes
  use Biogeophysics2Mod   , only : Biogeophysics2
  use SurfaceAlbedoMod    , only : SurfaceAlbedo
  use pft2colMod          , only : pft2col
  use CNSetValueMod       , only : CNZeroFluxes_dwt
  use CNEcosystemDynMod   , only : CNEcosystemDynNoLeaching,CNEcosystemDynLeaching
  use CNAnnualUpdateMod   , only : CNAnnualUpdate
  use CNBalanceCheckMod   , only : BeginCBalance, BeginNBalance, CBalanceCheck, NBalanceCheck
  use ndepStreamMod       , only : ndep_interp
  use CNVerticalProfileMod, only : decomp_vertprofiles
  use CNFireMod           , only : CNFireInterp
  use ActiveLayerMod      , only : alt_calc
  use DUSTMod             , only : DustDryDep, DustEmission
  use VOCEmissionMod      , only : VOCEmission
  use seq_drydep_mod      , only : n_drydep, drydep_method, DD_XLND
  use DryDepVelocity      , only : depvel_compute
  use ch4Mod              , only : ch4
  use abortutils          , only : endrun
  use UrbanMod            , only : UrbanAlbedo, UrbanRadiation, UrbanFluxes 
  use SNICARMod           , only : SnowAge_grain
  use DaylengthMod        , only : UpdateDaylength
  use clm_atmlnd          , only : clm_map2gcell, downscale_forcings
  use clm_glclnd          , only : update_clm_s2x
  use SatellitePhenologyMod, only : SatellitePhenology, interpMonthlyVeg
  use perf_mod
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_drv                 ! clm physics,history, restart writes
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: write_diagnostic        ! Write diagnostic information to log file
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
  integer  :: nstep                    ! time step number
  real(r8) :: dtime                    ! land model time step (sec)
  real(r8) :: t1, t2, t3               ! temporary for mass balance checks
  integer  :: nc, fc, c, fp, p, l, g   ! indices
  integer  :: nclumps                  ! number of clumps on this processor
  integer  :: yrp1                     ! year (0, ...) for nstep+1
  integer  :: monp1                    ! month (1, ..., 12) for nstep+1
  integer  :: dayp1                    ! day of month (1, ..., 31) for nstep+1
  integer  :: secp1                    ! seconds into current date for nstep+1
  integer  :: yr                       ! year (0, ...)
  integer  :: mon                      ! month (1, ..., 12)
  integer  :: day                      ! day of month (1, ..., 31)
  integer  :: sec                      ! seconds of the day
  integer  :: ncdate                   ! current date
  integer  :: nbdate                   ! base date (reference date)
  integer  :: kyr                      ! thousand years, equals 2 at end of first year
  character(len=256) :: filer          ! restart file name
  integer :: ier                       ! error code
  type(bounds_type) :: bounds_clump    ! bounds
  type(bounds_type) :: bounds_proc     ! bounds
  !-----------------------------------------------------------------------

  ! Determine processor bounds and clumps for this processor

  call get_proc_bounds(bounds_proc)
  nclumps = get_proc_clumps()

  ! Update time-related info

  call CropRestIncYear()

  ! ============================================================================
  ! Specified phenology
  ! ============================================================================
  
  if (use_cn) then 
     ! For dry-deposition need to call CLMSP so that mlaidiff is obtained
     if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
        call t_startf('interpMonthlyVeg')
        call interpMonthlyVeg(bounds_proc)
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
        call interpMonthlyVeg(bounds_proc)
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
      call alt_calc(filter(nc)%num_soilc, filter(nc)%soilc)

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
              filter_inactive_and_active(nc)%soilp)
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
      
      if (use_cn) then
         call t_startf('begcnbal')
         call BeginCBalance(bounds_clump, filter(nc)%num_soilc, filter(nc)%soilc)
         call BeginNBalance(bounds_clump, filter(nc)%num_soilc, filter(nc)%soilc)
         call CNZeroFluxes_dwt(bounds_clump)
         call t_stopf('begcnbal')
      end if
  end do
  !$OMP END PARALLEL DO

  ! ============================================================================
  ! Update subgrid weights with dynamic landcover (prescribed transient PFTs,
  ! CNDV, and or dynamic landunits), and do related adjustments. Note that this
  ! call needs to happen outside loops over nclumps.
  ! ============================================================================

  call t_startf('dyn_subgrid')
  call dynSubgrid_driver(bounds_proc)
  call t_stopf('dyn_subgrid')

  ! ============================================================================
  ! Initialize the mass balance checks for water.
  !
  ! Currently, I believe this needs to be done after weights are updated for
  ! prescribed transient PFTs or CNDV, because column-level water is not
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
      call BeginWaterBalance(bounds_clump, &
           filter(nc)%num_nolakec, filter(nc)%nolakec, &
           filter(nc)%num_lakec, filter(nc)%lakec, &
           filter(nc)%num_hydrologyc, filter(nc)%hydrologyc)
      call t_stopf('begwbal')
   end do
  

  ! ============================================================================
  ! Update dynamic N deposition field, on albedo timestep
  ! currently being done outside clumps loop, but no reason why it couldn't be
  ! re-written to go inside.
  ! ============================================================================
  
  call t_startf('ndep_interp')
  if (use_cn) then
     ! PET: switching CN timestep
     call ndep_interp(bounds_proc)
     call CNFireInterp(bounds_proc)
  end if
  call t_stopf('ndep_interp')

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
     ! initialize intracellular CO2 (Pa) parameters each timestep for use in VOCEmission
     pcf%cisun_z(bounds_clump%begp:bounds_clump%endp, :) = -999._r8
     pcf%cisha_z(bounds_clump%begp:bounds_clump%endp, :) = -999._r8

     call UpdateDaylength(bounds_clump, declin)
     
     call clm_driverInit(bounds_clump, &
          filter(nc)%num_nolakec, filter(nc)%nolakec)

     call downscale_forcings(bounds_clump, &
          filter(nc)%num_do_smb_c, filter(nc)%do_smb_c)

     call t_stopf('drvinit')

     ! ============================================================================
     ! Hydrology1
     ! ============================================================================

     call t_startf('hydro1')
     call Hydrology1(bounds_clump, &
          filter(nc)%num_nolakec, filter(nc)%nolakec, &
          filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('hydro1')

     ! ============================================================================
     ! Surface Radiation
     ! ============================================================================

     call t_startf('surfrad')

     ! Surface Radiation for non-urban columns

     call SurfaceRadiation(bounds_clump, &
          filter(nc)%num_nourbanp, filter(nc)%nourbanp)

     ! Surface Radiation for urban columns

     call UrbanRadiation(bounds_clump, &
          filter(nc)%num_nourbanl, filter(nc)%nourbanl, &
          filter(nc)%num_urbanl, filter(nc)%urbanl, &
          filter(nc)%num_urbanc, filter(nc)%urbanc, &
          filter(nc)%num_urbanp, filter(nc)%urbanp)

     call t_stopf('surfrad')

     ! ============================================================================
     ! Determine leaf temperature and surface fluxes based on ground
     ! temperature from previous time step.
     ! ============================================================================

     call t_startf('bgp1')
     call Biogeophysics1(bounds_clump, &
          filter(nc)%num_nolakec, filter(nc)%nolakec, &
          filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('bgp1')

     ! ============================================================================
     ! Determine bare soil or snow-covered vegetation surface temperature and fluxes
     ! Calculate Ground fluxes (frac_veg_nosno is either 1 or 0)
     ! ============================================================================

     call t_startf('bgflux')

     ! BareGroundFluxes for all pfts except lakes and urban landunits

     call BareGroundFluxes(bounds_clump, &
          filter(nc)%num_nolakeurbanp, filter(nc)%nolakeurbanp)
     call t_stopf('bgflux')

     ! Fluxes for all Urban landunits

     call t_startf('uflux')
     call UrbanFluxes(bounds_clump, &
          filter(nc)%num_nourbanl, filter(nc)%nourbanl, &
          filter(nc)%num_urbanl, filter(nc)%urbanl, &
          filter(nc)%num_urbanc, filter(nc)%urbanc, &
          filter(nc)%num_urbanp, filter(nc)%urbanp)
     call t_stopf('uflux')

     ! ============================================================================
     ! Determine non snow-covered vegetation surface temperature and fluxes
     ! Calculate canopy temperature, latent and sensible fluxes from the canopy,
     ! and leaf water change by evapotranspiration
     ! ============================================================================

     call t_startf('canflux')
     call CanopyFluxes(bounds_clump, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('canflux')

     ! ============================================================================
     ! Determine lake temperature and surface fluxes
     ! ============================================================================

     call t_startf('bgplake')
     call SLakeFluxes(bounds_clump, &
          filter(nc)%num_lakec, filter(nc)%lakec, &
          filter(nc)%num_lakep, filter(nc)%lakep)
     call SLakeTemperature(bounds_clump, &
          filter(nc)%num_lakec, filter(nc)%lakec, &
          filter(nc)%num_lakep, filter(nc)%lakep)
     call t_stopf('bgplake')

     ! ============================================================================
     ! DUST and VOC emissions
     ! ============================================================================

     call t_startf('bgc')

     ! Dust mobilization (C. Zender's modified codes)
     call DustEmission(bounds_clump, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep)

     ! Dust dry deposition (C. Zender's modified codes)
     call DustDryDep(bounds_clump)

     ! VOC emission (A. Guenther's MEGAN (2006) model)
     if (use_voc) then
        call VOCEmission(bounds_clump, &
             filter(nc)%num_soilp, filter(nc)%soilp)
     end if

     call t_stopf('bgc')


     ! ============================================================================
     ! Determine soil/snow temperatures including ground temperature and
     ! update surface fluxes for new ground temperature.
     ! ============================================================================

     call t_startf('bgp2')
     call Biogeophysics2(bounds_clump, &
          filter(nc)%num_urbanl,  filter(nc)%urbanl, &
          filter(nc)%num_nolakec, filter(nc)%nolakec, &
          filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('bgp2')

     ! ============================================================================
     ! Perform averaging from PFT level to column level
     ! ============================================================================

     call t_startf('pft2col')
     call pft2col(bounds_clump, filter(nc)%num_nolakec, filter(nc)%nolakec)
     call t_stopf('pft2col')


     ! ============================================================================
     ! Vertical (column) soil and surface hydrology
     ! ============================================================================

     call t_startf('hydro2 without drainage')
     call Hydrology2NoDrainage(bounds_clump, &
          filter(nc)%num_nolakec, filter(nc)%nolakec, &
          filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
          filter(nc)%num_urbanc, filter(nc)%urbanc, &
          filter(nc)%num_snowc, filter(nc)%snowc, &
          filter(nc)%num_nosnowc, filter(nc)%nosnowc)
                     
     call t_stopf('hydro2 without drainage')
     
     ! ============================================================================
     ! Lake hydrology
     ! ============================================================================

     call t_startf('hylake')
     call SLakeHydrology(bounds_clump, &
          filter(nc)%num_lakec, filter(nc)%lakec, &
          filter(nc)%num_lakep, filter(nc)%lakep)
     call t_stopf('hylake')

     ! ============================================================================
     ! ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)
     ! ============================================================================

     call t_startf('snow_init')
     do c = bounds_clump%begc,bounds_clump%endc
        l = col%landunit(c)
        if (lun%urbpoi(l)) then
           ! Urban landunit use Bonan 1996 (LSM Technical Note)
           cps%frac_sno(c) = min( cps%snow_depth(c)/0.05_r8, 1._r8)
        end if
     end do

     ! ============================================================================
     ! Snow aging routine based on Flanner and Zender (2006), Linking snowpack 
     ! microphysics and albedo evolution, JGR, and Brun (1989), Investigation of 
     ! wet-snow metamorphism in respect of liquid-water content, Ann. Glaciol.
     ! ============================================================================
     ! Note the snow filters here do not include lakes; SnowAge_grain is called
     ! for lakes from SLakeHydrology.

     call SnowAge_grain(bounds_clump, &
          filter(nc)%num_snowc, filter(nc)%snowc, &
          filter(nc)%num_nosnowc, filter(nc)%nosnowc)
     call t_stopf('snow_init')

     ! ============================================================================
     ! Ecosystem dynamics: Uses CN, CNDV, or static parameterizations
     ! ============================================================================
     call t_startf('ecosysdyn')

     if (use_cn) then

        ! fully prognostic canopy structure and C-N biogeochemistry
        ! - CNDV defined: prognostic biogeography; else prescribed
        ! - crop model:   crop algorithms called from within CNEcosystemDyn
        call CNEcosystemDynNoLeaching(bounds_clump, &
             filter(nc)%num_soilc,&
             filter(nc)%soilc, filter(nc)%num_soilp, &
             filter(nc)%soilp, filter(nc)%num_pcropp, &
             filter(nc)%pcropp, doalb)
        
        call CNAnnualUpdate(bounds_clump, &
             filter(nc)%num_soilc,&
             filter(nc)%soilc, filter(nc)%num_soilp, &
             filter(nc)%soilp)

     else

        ! Prescribed biogeography,
        ! prescribed canopy structure, some prognostic carbon fluxes
        call SatellitePhenology(bounds_clump, &
             filter(nc)%num_nolakep, filter(nc)%nolakep, &
             doalb)

     end if
     call t_stopf('ecosysdyn')

     ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)
     call t_startf('depvel')
     call depvel_compute(bounds_clump)
     call t_stopf('depvel')

     if (use_lch4) then
        call t_startf('ch4')
        call ch4 (bounds_clump, &
             filter(nc)%num_soilc, filter(nc)%soilc, &
             filter(nc)%num_lakec, filter(nc)%lakec, &
             filter(nc)%num_soilp, filter(nc)%soilp)
        call t_stopf('ch4')
     end if
     
     
     call t_startf('hydro2 drainage')
     call Hydrology2Drainage(bounds_clump, &
                     filter(nc)%num_nolakec, filter(nc)%nolakec, &
                     filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
                     filter(nc)%num_urbanc, filter(nc)%urbanc, &
                     filter(nc)%num_do_smb_c, filter(nc)%do_smb_c)                
     
     call t_stopf('hydro2 drainage')     
     ! ============================================================================
     ! Check the energy and water balance, also carbon and nitrogen balance
     ! ============================================================================
     if (use_cn) then
        call CNEcosystemDynLeaching(bounds_clump, &
             filter(nc)%num_soilc,&
             filter(nc)%soilc, filter(nc)%num_soilp, &
             filter(nc)%soilp, filter(nc)%num_pcropp, &
             filter(nc)%pcropp, doalb)
     end if
     
     call t_startf('balchk')
     call BalanceCheck(bounds_clump,filter(nc)%num_do_smb_c,filter(nc)%do_smb_c)
     call t_stopf('balchk')

     if (use_cn) then
        nstep = get_nstep()
        if (nstep .lt. 2 )then
           if (masterproc) write(iulog,*) '--WARNING-- skipping CN balance check for first timestep'
        else
           call t_startf('cnbalchk')
           call CBalanceCheck(bounds_clump, filter(nc)%num_soilc, filter(nc)%soilc)
           call NBalanceCheck(bounds_clump, filter(nc)%num_soilc, filter(nc)%soilc)
           call t_stopf('cnbalchk')
        end if
     end if

     ! ============================================================================
     ! Determine albedos for next time step
     ! ============================================================================

     if (doalb) then

        ! Albedos for non-urban columns
        call t_startf('surfalb')
        call SurfaceAlbedo(bounds_clump, &
             filter_inactive_and_active(nc)%num_nourbanc, &
             filter_inactive_and_active(nc)%nourbanc, &
             filter_inactive_and_active(nc)%num_nourbanp, &
             filter_inactive_and_active(nc)%nourbanp, &
             nextsw_cday, declinp1)
        call t_stopf('surfalb')

        ! Albedos for urban columns
        if (filter_inactive_and_active(nc)%num_urbanl > 0) then
           call t_startf('urbsurfalb')
           call UrbanAlbedo(bounds_clump, &
                filter_inactive_and_active(nc)%num_urbanl, &
                filter_inactive_and_active(nc)%urbanl, &
                filter_inactive_and_active(nc)%num_urbanc, &
                filter_inactive_and_active(nc)%urbanc, &
                filter_inactive_and_active(nc)%num_urbanp, &
                filter_inactive_and_active(nc)%urbanp)
           call t_stopf('urbsurfalb')
        end if

     end if
     
  end do
  !$OMP END PARALLEL DO

  ! ============================================================================
  ! Determine gridcell averaged properties to send to atm to send to glc
  ! ============================================================================

  call t_startf('clm_map2gcell')
  call clm_map2gcell(bounds_proc)
  call t_stopf('clm_map2gcell')
  
  !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
  do nc = 1,nclumps
     call get_clump_bounds(nc, bounds_clump)
  
     if (create_glacier_mec_landunit) then  
        call t_startf('create_s2x')
        call update_clm_s2x(bounds_clump, &
             filter(nc)%num_do_smb_c, filter(nc)%do_smb_c, &
             init=.false.)
        call t_stopf('create_s2x')
     end if
  end do
  !$OMP END PARALLEL DO

  ! ============================================================================
  ! Write global average diagnostics to standard output
  ! ============================================================================

  nstep = get_nstep()
  if (wrtdia) call mpi_barrier(mpicom,ier)
  call t_startf('wrtdiag')
  call write_diagnostic(wrtdia, nstep)
  call t_stopf('wrtdiag')

  ! ============================================================================
  ! Update accumulators
  ! ============================================================================

  call t_startf('accum')
  call updateAccFlds(bounds_proc)
  call t_stopf('accum')

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
           call dv(bounds_clump, &
                filter(nc)%num_natvegp, filter(nc)%natvegp, kyr)
        end do
        !$OMP END PARALLEL DO
     end if
     call t_stopf('d2dgvm')
  end if

  
  ! ============================================================================
  ! History/Restart output
  ! ============================================================================

  if (.not. use_noio) then

     call t_startf('clm_drv_io')

     ! Create history and write history tapes if appropriate
     call t_startf('clm_drv_io_htapes')
     call hist_htapes_wrapup( rstwr, nlend )
     call t_stopf('clm_drv_io_htapes')

     ! Write to CNDV history buffer if appropriate
     if (use_cndv) then
        if (monp1==1 .and. dayp1==1 .and. secp1==dtime .and. nstep>0)  then
           call t_startf('clm_drv_io_hdgvm')
           call histCNDV( bounds_proc )
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

!------------------------------------------------------------------------
subroutine write_diagnostic (wrtdia, nstep)
  !
  ! !DESCRIPTION:
  ! Write diagnostic surface temperature output each timestep.  Written to
  ! be fast but not bit-for-bit because order of summations can change each
  ! timestep.
  !
  ! !USES:
  use clm_atmlnd , only : clm_l2a
  use decompMod  , only : get_proc_bounds, get_proc_global
  use spmdMod    , only : masterproc, npes, MPI_REAL8, MPI_ANY_SOURCE
  use spmdMod    , only : MPI_STATUS_SIZE, mpicom, MPI_SUM
  use shr_sys_mod, only : shr_sys_flush
  use abortutils , only : endrun
  use shr_log_mod, only : errMsg => shr_log_errMsg
  !
  ! !ARGUMENTS:
  implicit none
  logical, intent(in) :: wrtdia     !true => write diagnostic
  integer, intent(in) :: nstep      !model time step
  !
  ! !REVISION HISTORY:
  ! Created by Mariana Vertenstein
  !
  ! !LOCAL VARIABLES:
  integer :: p                       ! loop index
  integer :: numg                    ! total number of gridcells across all processors
  integer :: numl                    ! total number of landunits across all processors
  integer :: numc                    ! total number of columns across all processors
  integer :: nump                    ! total number of pfts across all processors
  integer :: ier                     ! error status
  real(r8):: psum                    ! partial sum of ts
  real(r8):: tsum                    ! sum of ts
  real(r8):: tsxyav                  ! average ts for diagnostic output
  integer :: status(MPI_STATUS_SIZE) ! mpi status
  type(bounds_type) :: bounds        ! bounds
  logical,parameter :: old_sendrecv = .false.  ! Flag if should use old send/receive method rather than MPI reduce
!------------------------------------------------------------------------

  call get_proc_bounds(bounds)
  call get_proc_global(numg, numl, numc, nump)

  if (wrtdia) then

     call t_barrierf('sync_write_diag', mpicom)
     psum = sum(clm_l2a%t_rad(bounds%begg:bounds%endg))
     if (old_sendrecv) then
        if (masterproc) then
           tsum = psum
           do p = 1, npes-1
              call mpi_recv(psum, 1, MPI_REAL8, p, 999, mpicom, status, ier)
              if (ier/=0) then
                 write(iulog,*) 'write_diagnostic: Error in mpi_recv()'
                 call endrun(msg=errMsg(__FILE__, __LINE__))
              end if
              tsum = tsum + psum
           end do
        else
           call mpi_send(psum, 1, MPI_REAL8, 0, 999, mpicom, ier)
           if (ier/=0) then
              write(iulog,*) 'write_diagnostic: Error in mpi_send()'
              call endrun(msg=errMsg(__FILE__, __LINE__))
           end if
        end if
     else
        call mpi_reduce(psum, tsum, 1, MPI_REAL8, MPI_SUM, 0, mpicom, ier)
        if (ier/=0) then
           write(iulog,*) 'write_diagnostic: Error in mpi_reduce()'
           call endrun(msg=errMsg(__FILE__, __LINE__))
        end if
     endif
     if (masterproc) then
        tsxyav = tsum / numg
        write(iulog,1000) nstep, tsxyav
        call shr_sys_flush(iulog)
     end if

  else

     if (masterproc) then
        write(iulog,*)'clm2: completed timestep ',nstep
        call shr_sys_flush(iulog)
     end if

  endif

1000 format (1x,'nstep = ',i10,'   TS = ',f21.15)

end subroutine write_diagnostic

end module clm_driver
