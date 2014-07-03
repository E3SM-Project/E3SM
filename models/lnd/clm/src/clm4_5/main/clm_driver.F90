module clm_driver

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_driver
!
! !DESCRIPTION:
! This module provides the main CLM driver physics calling sequence.  Most
! computations occurs over ``clumps'' of gridcells (and associated subgrid
! scale entities) assigned to each MPI process. Computation is further
! parallelized by looping over clumps on each process using shared memory OpenMP.
!
! The main CLM driver physics calling sequence for clm_driver1 is as follows:
! \begin{verbatim}
!
!     + interpMonthlyVeg      interpolate monthly vegetation data        [! CN or ! CNDV]
!       + readMonthlyVegetation read vegetation data for two months      [! CN or ! CNDV]
!
! ==== Begin Loop over clumps ====
!  -> dynland_hwcontent   Get initial heat, water content
!     + pftdyn_interp                                                    [pftdyn]
!     + dynland_hwcontent   Get new heat, water content                  [pftdyn]
! ==== End Loop over clumps  ====
!
! ==== Begin Loop over clumps ====
!  -> clm_driverInit      save of variables from previous time step
!  -> Hydrology1          canopy interception and precip on ground
!     -> FracWet          fraction of wet vegetated surface and dry elai
!  -> SurfaceRadiation    surface solar radiation
!  -> UrbanRadiation      surface solar and longwave radiation for Urban landunits
!  -> Biogeophysics1      leaf temperature and surface fluxes
!  -> BareGroundFluxes    surface fluxes for bare soil or snow-covered
!                         vegetation patches
!  -> UrbanFluxes         surface fluxes for urban landunits
!     -> MoninObukIni     first-guess Monin-Obukhov length and wind speed
!     -> FrictionVelocity friction velocity and potential temperature and
!                         humidity profiles
!  -> CanopyFluxes        leaf temperature and surface fluxes for vegetated
!                         patches
!     -> QSat             saturated vapor pressure, specific humidity, &
!                         derivatives at leaf surface
!     -> MoninObukIni     first-guess Monin-Obukhov length and wind speed
!     -> FrictionVelocity friction velocity and potential temperature and
!                         humidity profiles
!     -> Stomata          stomatal resistance and photosynthesis for
!                         sunlit leaves
!     -> Stomata          stomatal resistance and photosynthesis for
!                         shaded leaves
!     -> QSat             recalculation of saturated vapor pressure,
!                         specific humidity, & derivatives at leaf surface
!   + DustEmission        Dust mobilization
!   + DustDryDep          Dust dry deposition
!  -> SLakeFluxes         lake surface fluxes
!  -> SLakeTemperature    lake temperature
!   + VOCEmission         compute VOC emission                          [VOC]
!  -> Biogeophysics2      soil/snow & ground temp and update surface fluxes
!  -> pft2col             Average from PFT level to column level
!  -> Hydrology2          surface and soil hydrology
!  -> SLakeHydrology   lake hydrology
!  -> SnowAge_grain       update snow effective grain size for snow radiative transfer
!   + CNEcosystemDyn      Carbon Nitrogen model ecosystem dynamics:     [CN]
!                         vegetation phenology and soil carbon  
!   + EcosystemDyn        "static" ecosystem dynamics:                  [! CN ]
!                         vegetation phenology and soil carbon  
!  -> BalanceCheck        check for errors in energy and water balances
!  -> SurfaceAlbedo       albedos for next time step
!  -> UrbanAlbedo         Urban landunit albedos for next time step
!  ====  End Loop over clumps  ====
!
! Second phase of the clm main driver, for handling history and restart file output.
!
!  -> write_diagnostic    output diagnostic if appropriate
!   + inicPerp            initial snow and soil moisture               [is_perpetual]
!  -> updateAccFlds       update accumulated fields
!  -> hist_update_hbuf    accumulate history fields for time interval
!  -> htapes_wrapup       write history tapes if appropriate
!  -> restFile_write      write restart file if appropriate
! \end{verbatim}
!
! Optional subroutines are denoted by an plus (+) with the associated
! CPP token or variable in brackets at the end of the line.
!
! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clmtype
  use clm_varctl          , only : wrtdia, fpftdyn, iulog, create_glacier_mec_landunit
  use spmdMod             , only : masterproc,mpicom
  use decompMod           , only : get_proc_clumps, get_clump_bounds, get_proc_bounds
  use filterMod           , only : filter
  use reweightMod         , only : reweightWrapup
#if (defined CNDV)
  use CNDVMod             , only : dv, histCNDV
  use pftdynMod           , only : pftwt_interp
#endif
  use pftdynMod           , only : pftdyn_interp, pftdyn_wbal_init, pftdyn_wbal
#ifdef CN
  use pftdynMod           , only : pftdyn_cnbal
#endif
  use dynlandMod          , only : dynland_hwcontent
  use clm_varcon          , only : zlnd, isturb
  use clm_time_manager    , only : get_step_size,get_curr_date,get_ref_date,get_nstep,is_perpetual
  use histFileMod         , only : hist_update_hbuf, hist_htapes_wrapup
  use restFileMod         , only : restFile_write, restFile_filename
  use inicPerpMod         , only : inicPerp  
  use accFldsMod          , only : updateAccFlds
  use clm_driverInitMod   , only : clm_driverInit
  use BalanceCheckMod     , only : BeginWaterBalance, BalanceCheck
  use SurfaceRadiationMod , only : SurfaceRadiation
  use Hydrology1Mod       , only : Hydrology1
  use Hydrology2Mod       , only : Hydrology2
  use SLakeFluxesMod   , only : SLakeFluxes
  use SLakeTemperatureMod, only : SLakeTemperature
  use SLakeHydrologyMod, only : SLakeHydrology
  use Biogeophysics1Mod   , only : Biogeophysics1
  use BareGroundFluxesMod , only : BareGroundFluxes
  use CanopyFluxesMod     , only : CanopyFluxes
  use Biogeophysics2Mod   , only : Biogeophysics2
  use SurfaceAlbedoMod    , only : SurfaceAlbedo
  use pft2colMod          , only : pft2col
#if (defined CN)
  use CNSetValueMod       , only : CNZeroFluxes_dwt
  use CNEcosystemDynMod   , only : CNEcosystemDyn
  use CNAnnualUpdateMod   , only : CNAnnualUpdate
  use CNBalanceCheckMod   , only : BeginCBalance, BeginNBalance, &
                                   CBalanceCheck, NBalanceCheck
  use ndepStreamMod       , only : ndep_interp
  use CNVerticalProfileMod, only: decomp_vertprofiles
  use CNFireMod           , only : CNFireInterp
#else
  use STATICEcosysDynMod  , only : EcosystemDyn
#endif
  use ActiveLayerMod      , only : alt_calc
  use DUSTMod             , only : DustDryDep, DustEmission
  use VOCEmissionMod      , only : VOCEmission
  use seq_drydep_mod      , only : n_drydep, drydep_method, DD_XLND
  use STATICEcosysDynMod  , only : interpMonthlyVeg
  use DryDepVelocity      , only : depvel_compute
#if (defined LCH4)
  use ch4Mod              , only : ch4
#endif
  use abortutils          , only : endrun
  use UrbanMod            , only : UrbanAlbedo, UrbanRadiation, UrbanFluxes 
  use SNICARMod           , only : SnowAge_grain
  use clm_atmlnd          , only : clm_map2gcell
  use clm_glclnd          , only : create_clm_s2x
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
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: clm_drv
!
! !INTERFACE:
subroutine clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate)
!
! !DESCRIPTION:
!
! First phase of the clm driver calling the clm physics. An outline of
! the calling tree is given in the description of this module.
!
! !USES:

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
! !REVISION HISTORY:
! 2002.10.01  Mariana Vertenstein latest update to new data structures
! 11/26/03, Peter Thornton: Added new call for SurfaceRadiationSunShade when
!  cpp directive SUNSHA is set, for sunlit/shaded canopy radiation.
! 4/25/05, Peter Thornton: Made the sun/shade routine the default, no longer
!  need to have SUNSHA defined.  
! Oct/05 & Jul/07 Sam Levis: Starting dates of CNDV and crop model work
! 2/29/08, Dave Lawrence: Revised snow cover fraction according to Niu and Yang, 2007
! 3/6/09, Peter Thornton: Added declin as new argument, for daylength control on Vcmax
! 2008.11.12  B. Kauffman: morph routine casa() in casa_ecosytemDyn(), so casa
!    is more similar to CN & DGVM
! 2/25/2012 M. Vertenstein: Removed CASA references 
!
!EOP
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
  integer , pointer :: clandunit(:) ! landunit index associated with each column
  integer , pointer :: itypelun(:)  ! landunit type
!
! !OTHER LOCAL VARIABLES:
  integer  :: nstep                    ! time step number
  real(r8) :: dtime                    ! land model time step (sec)
  real(r8) :: t1, t2, t3               ! temporary for mass balance checks
  integer  :: nc, fc, c, fp, p, l, g   ! indices
  integer  :: nclumps                  ! number of clumps on this processor
  integer  :: begg, endg               ! clump beginning and ending gridcell indices
  integer  :: begl, endl               ! clump beginning and ending landunit indices
  integer  :: begc, endc               ! clump beginning and ending column indices
  integer  :: begp, endp               ! clump beginning and ending pft indices
  integer  :: begg_proc, endg_proc     ! proc beginning and ending gridcell indices
  integer  :: begl_proc, endl_proc     ! proc beginning and ending landunit indices
  integer  :: begc_proc, endc_proc     ! proc beginning and ending column indices
  integer  :: begp_proc, endp_proc     ! proc beginning and ending pft indices
  type(column_type), pointer :: cptr   ! pointer to column derived subtype
#if (defined CNDV)
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
#endif
  character(len=256) :: filer          ! restart file name
  integer :: ier                       ! error code
!-----------------------------------------------------------------------

  ! Assign local pointers to derived subtypes components (landunit-level)

  itypelun            => clm3%g%l%itype

  ! Assign local pointers to derived subtypes components (column-level)

  clandunit           => clm3%g%l%c%landunit

  ! Set pointers into derived type

  cptr => clm3%g%l%c

#ifdef CN
  ! For dry-deposition need to call CLMSP so that mlaidiff is obtained
  if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
     call t_startf('interpMonthlyVeg')
     call interpMonthlyVeg()
     call t_stopf('interpMonthlyVeg')
  endif
#else
  ! ============================================================================
  ! Determine weights for time interpolation of monthly vegetation data.
  ! This also determines whether it is time to read new monthly vegetation and
  ! obtain updated leaf area index [mlai1,mlai2], stem area index [msai1,msai2],
  ! vegetation top [mhvt1,mhvt2] and vegetation bottom [mhvb1,mhvb2]. The
  ! weights obtained here are used in subroutine ecosystemdyn to obtain time
  ! interpolated values.
  ! ============================================================================
  if (doalb .or. ( n_drydep > 0 .and. drydep_method == DD_XLND )) then
     call t_startf('interpMonthlyVeg')
     call interpMonthlyVeg()
     call t_stopf('interpMonthlyVeg')
  end if
#endif

  ! ============================================================================
  ! Loop over clumps
  ! ============================================================================

  nclumps = get_proc_clumps()

#if (defined _OPENMP)
  call endrun( 'ERROR: OpenMP does NOT work on the CLM4.5 science branch. Set the number of threads for LND to 1 and rerun.' )
#endif

!$OMP PARALLEL DO PRIVATE (nc,g,begg,endg,begl,endl,begc,endc,begp,endp)
  do nc = 1,nclumps

     ! ============================================================================
     ! Determine clump boundaries
     ! ============================================================================

     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

     ! ============================================================================
     ! change pft weights and compute associated heat & water fluxes
     ! ============================================================================

     call t_startf("init_hwcontent")
     ! initialize heat and water content and dynamic balance fields to zero
     do g = begg,endg
        clm3%g%gwf%qflx_liq_dynbal(g) = 0._r8
        clm3%g%gws%gc_liq2(g)         = 0._r8
        clm3%g%gws%gc_liq1(g)         = 0._r8
        clm3%g%gwf%qflx_ice_dynbal(g) = 0._r8
        clm3%g%gws%gc_ice2(g)         = 0._r8 
        clm3%g%gws%gc_ice1(g)         = 0._r8
        clm3%g%gef%eflx_dynbal(g)     = 0._r8
        clm3%g%ges%gc_heat2(g)        = 0._r8
        clm3%g%ges%gc_heat1(g)        = 0._r8
     enddo

     !--- get initial heat,water content ---
      call dynland_hwcontent( begg, endg, clm3%g%gws%gc_liq1(begg:endg), &
                              clm3%g%gws%gc_ice1(begg:endg), clm3%g%ges%gc_heat1(begg:endg) )
     call t_stopf("init_hwcontent")
   end do
!$OMP END PARALLEL DO


   ! ==================================================================================
   ! Determine decomp vertical profiles
   !
   ! These routines (alt_calc & decomp_vertprofiles) need to be called before
   ! pftdyn_cnbal, and it appears that they need to be called before pftdyn_interp and
   ! the associated filter updates, too (otherwise we get a carbon balance error)
   ! ==================================================================================

!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
   do nc = 1,nclumps
      call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
      
      call t_startf("decomp_vert")
      call alt_calc(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)

#if (defined CN)
      call decomp_vertprofiles(begp, endp, begc, endc, filter(nc)%num_soilc, filter(nc)%soilc, filter(nc)%num_soilp, filter(nc)%soilp)
#endif

      call t_stopf("decomp_vert")
   end do
!$OMP END PARALLEL DO


#if (!defined CNDV)
   if (fpftdyn /= ' ') then
      call t_startf("pftdyn_interp")
      call pftdyn_interp  ! change the pft weights
      
!$OMP PARALLEL DO PRIVATE (nc,g,begg,endg,begl,endl,begc,endc,begp,endp)
      do nc = 1,nclumps
         call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

         ! do stuff that needs to be done after changing weights
         ! This call should be made as soon as possible after pftdyn_interp
         call reweightWrapup(nc)

         !--- get new heat,water content: (new-old)/dt = flux into lnd model ---
         call dynland_hwcontent( begg, endg, clm3%g%gws%gc_liq2(begg:endg), &
                                 clm3%g%gws%gc_ice2(begg:endg), clm3%g%ges%gc_heat2(begg:endg) )
         dtime = get_step_size()
         do g = begg,endg
            clm3%g%gwf%qflx_liq_dynbal(g) = (clm3%g%gws%gc_liq2 (g) - clm3%g%gws%gc_liq1 (g))/dtime
            clm3%g%gwf%qflx_ice_dynbal(g) = (clm3%g%gws%gc_ice2 (g) - clm3%g%gws%gc_ice1 (g))/dtime
            clm3%g%gef%eflx_dynbal    (g) = (clm3%g%ges%gc_heat2(g) - clm3%g%ges%gc_heat1(g))/dtime
         enddo
      end do
!$OMP END PARALLEL DO
      call t_stopf("pftdyn_interp")
   end if
#endif
      
!$OMP PARALLEL DO PRIVATE (nc,g,begg,endg,begl,endl,begc,endc,begp,endp)
  do nc = 1,nclumps
     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

     ! ============================================================================
     ! Initialize the mass balance checks: water, carbon, and nitrogen
     ! ============================================================================

     call t_startf('begwbal')
     call BeginWaterBalance(begc, endc, begp, endp, &
          filter(nc)%num_nolakec, filter(nc)%nolakec, filter(nc)%num_lakec, filter(nc)%lakec, &
          filter(nc)%num_hydrologyc, filter(nc)%hydrologyc)
     call t_stopf('begwbal')

#if (defined CN)
     call t_startf('begcnbal')
     call BeginCBalance(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)
     call BeginNBalance(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)
     call t_stopf('begcnbal')
#endif

  end do
!$OMP END PARALLEL DO

  ! ============================================================================
  ! Initialize h2ocan_loss to zero
  ! ============================================================================
  call t_startf('pftdynwts')
!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
  do nc = 1,nclumps
     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
     call pftdyn_wbal_init( begc, endc )

#if (defined CNDV)
     ! NOTE: Currently CNDV and fpftdyn /= ' ' are incompatible
     call CNZeroFluxes_dwt( begc, endc, begp, endp )
     call pftwt_interp( begp, endp )
     call reweightWrapup(nc)
     call pftdyn_wbal( begg, endg, begc, endc, begp, endp )
     call pftdyn_cnbal( begc, endc, begp, endp )
#else

#if (defined CN)
     call CNZeroFluxes_dwt( begc, endc, begp, endp )
     if (fpftdyn /= ' ') then
        call pftdyn_cnbal( begc, endc, begp, endp )
     end if
#endif
#endif
  end do
!$OMP END PARALLEL DO


#if (defined CN)
  ! ============================================================================
  ! Update dynamic N deposition field, on albedo timestep
  ! currently being done outside clumps loop, but no reason why it couldn't be
  ! re-written to go inside.
  ! ============================================================================
  ! PET: switching CN timestep
  call ndep_interp()
  call CNFireInterp()
#endif
  call t_stopf('pftdynwts')

!$OMP PARALLEL DO PRIVATE (nc,l,c,begg,endg,begl,endl,begc,endc,begp,endp)
  do nc = 1,nclumps

     ! ============================================================================
     ! Determine clump boundaries
     ! ============================================================================

     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

     ! ============================================================================
     ! Initialize variables from previous time step and
     ! Determine canopy interception and precipitation onto ground surface.
     ! Determine the fraction of foliage covered by water and the fraction
     ! of foliage that is dry and transpiring. Initialize snow layer if the
     ! snow accumulation exceeds 10 mm.
     ! ============================================================================
     
     call t_startf('drvinit')
     ! initialize intracellular CO2 (Pa) parameters each timestep for use in VOCEmission
     clm3%g%l%c%p%pcf%cisun_z(begp:endp,:) = -999._r8
     clm3%g%l%c%p%pcf%cisha_z(begp:endp,:) = -999._r8

     ! initialize declination for current timestep
     do c = begc,endc
        clm3%g%l%c%cps%decl(c) = declin
     end do
     
     call clm_driverInit(begc, endc, begp, endp, &
          filter(nc)%num_nolakec, filter(nc)%nolakec, filter(nc)%num_lakec, filter(nc)%lakec)
     call t_stopf('drvinit')

     ! ============================================================================
     ! Hydrology1
     ! ============================================================================

     call t_startf('hydro1')
     call Hydrology1(begc, endc, begp, endp, &
                     filter(nc)%num_nolakec, filter(nc)%nolakec, &
                     filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('hydro1')

     ! ============================================================================
     ! Surface Radiation
     ! ============================================================================

     call t_startf('surfrad')

     ! Surface Radiation for non-urban columns

     call SurfaceRadiation(begp, endp, &
                           filter(nc)%num_nourbanp, filter(nc)%nourbanp)

     ! Surface Radiation for urban columns

     call UrbanRadiation(nc, begl, endl, begc, endc, begp, endp, &
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
     call Biogeophysics1(begg, endg, begc, endc, begp, endp, &
                         filter(nc)%num_nolakec, filter(nc)%nolakec, &
                         filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('bgp1')

     ! ============================================================================
     ! Determine bare soil or snow-covered vegetation surface temperature and fluxes
     ! Calculate Ground fluxes (frac_veg_nosno is either 1 or 0)
     ! ============================================================================

     call t_startf('bgflux')

     ! BareGroundFluxes for all pfts except lakes and urban landunits

     call BareGroundFluxes(begp, endp, &
                           filter(nc)%num_nolakeurbanp, filter(nc)%nolakeurbanp)
     call t_stopf('bgflux')

     ! Fluxes for all Urban landunits

     call t_startf('uflux')
     call UrbanFluxes(nc, begp, endp, begl, endl, begc, endc, &
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
     call CanopyFluxes(begg, endg, begc, endc, begp, endp, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('canflux')

     ! ============================================================================
     ! Determine lake temperature and surface fluxes
     ! ============================================================================

     call t_startf('bgplake')
     call SLakeFluxes(begc, endc, begp, endp, &
                         filter(nc)%num_lakec, filter(nc)%lakec, &
                         filter(nc)%num_lakep, filter(nc)%lakep)
     call SLakeTemperature(begc, endc, begp, endp, &
                              filter(nc)%num_lakec, filter(nc)%lakec, &
                              filter(nc)%num_lakep, filter(nc)%lakep)
     call t_stopf('bgplake')

     ! ============================================================================
     ! DUST and VOC emissions
     ! ============================================================================

     call t_startf('bgc')

     ! Dust mobilization (C. Zender's modified codes)
     call DustEmission(begp, endp, begc, endc, begl, endl, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep)

     ! Dust dry deposition (C. Zender's modified codes)
     call DustDryDep(begp, endp)

     ! VOC emission (A. Guenther's MEGAN (2006) model)
     call VOCEmission(begp, endp, &
                      filter(nc)%num_soilp, filter(nc)%soilp)

     call t_stopf('bgc')

     ! ============================================================================
     ! Determine soil/snow temperatures including ground temperature and
     ! update surface fluxes for new ground temperature.
     ! ============================================================================

     call t_startf('bgp2')
     call Biogeophysics2(begl, endl, begc, endc, begp, endp, &
                         filter(nc)%num_urbanl,  filter(nc)%urbanl, &
                         filter(nc)%num_nolakec, filter(nc)%nolakec, &
                         filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('bgp2')

     ! ============================================================================
     ! Perform averaging from PFT level to column level
     ! ============================================================================

     call t_startf('pft2col')
     call pft2col(begc, endc, filter(nc)%num_nolakec, filter(nc)%nolakec)
     call t_stopf('pft2col')

     ! ============================================================================
     ! Vertical (column) soil and surface hydrology
     ! ============================================================================

     call t_startf('hydro2')
     call Hydrology2(begc, endc, begp, endp, &
                     filter(nc)%num_nolakec, filter(nc)%nolakec, &
                     filter(nc)%num_hydrologyc, filter(nc)%hydrologyc, &
                     filter(nc)%num_urbanc, filter(nc)%urbanc, &
                     filter(nc)%num_snowc, filter(nc)%snowc, &
                     filter(nc)%num_nosnowc, filter(nc)%nosnowc)
     call t_stopf('hydro2')

     ! ============================================================================
     ! Lake hydrology
     ! ============================================================================

     call t_startf('hylake')
     call SLakeHydrology(begc, endc, begp, endp, filter(nc)%num_lakec, filter(nc)%lakec, &
                            filter(nc)%num_lakep, filter(nc)%lakep)
     call t_stopf('hylake')

     ! ============================================================================
     ! ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)
     ! ============================================================================

     call t_startf('snow_init')
     do c = begc,endc
        l = clandunit(c)
        if (itypelun(l) == isturb) then
           ! Urban landunit use Bonan 1996 (LSM Technical Note)
           cptr%cps%frac_sno(c) = min( cptr%cps%snow_depth(c)/0.05_r8, 1._r8)
        end if
     end do

     ! ============================================================================
     ! Snow aging routine based on Flanner and Zender (2006), Linking snowpack 
     ! microphysics and albedo evolution, JGR, and Brun (1989), Investigation of 
     ! wet-snow metamorphism in respect of liquid-water content, Ann. Glaciol.
     ! ============================================================================
     ! Note the snow filters here do not include lakes; SnowAge_grain is called
     ! for lakes from SLakeHydrology.

     call SnowAge_grain(begc, endc, &
          filter(nc)%num_snowc, filter(nc)%snowc, &
          filter(nc)%num_nosnowc, filter(nc)%nosnowc)
     call t_stopf('snow_init')

     ! ============================================================================
     ! Ecosystem dynamics: Uses CN, CNDV, or static parameterizations
     ! ============================================================================
     call t_startf('ecosysdyn')

#if (defined CN)
     ! fully prognostic canopy structure and C-N biogeochemistry
     ! - CNDV defined: prognostic biogeography; else prescribed
     ! - crop model:   crop algorithms called from within CNEcosystemDyn
     call CNEcosystemDyn(begc,endc,begp,endp,filter(nc)%num_soilc,&
                  filter(nc)%soilc, filter(nc)%num_soilp, &
                  filter(nc)%soilp, filter(nc)%num_pcropp, &
                  filter(nc)%pcropp, doalb)
     call CNAnnualUpdate(begc,endc,begp,endp,filter(nc)%num_soilc,&
                  filter(nc)%soilc, filter(nc)%num_soilp, &
                  filter(nc)%soilp)
#else
     ! Prescribed biogeography,
     ! prescribed canopy structure, some prognostic carbon fluxes
     call EcosystemDyn(begp, endp, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep, &
                       doalb)
#endif
     call t_stopf('ecosysdyn')

     ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)
     call t_startf('depvel')
     call depvel_compute(begp,endp)
     call t_stopf('depvel')

     ! ============================================================================
     ! Check the energy and water balance, also carbon and nitrogen balance
     ! ============================================================================

     call t_startf('balchk')
     call BalanceCheck(begp, endp, begc, endc, begl, endl, begg, endg)
     call t_stopf('balchk')

#if (defined CN)
     nstep = get_nstep()
     if (nstep .lt. 2 )then
        if (masterproc) write(iulog,*) '--WARNING-- skipping CN balance check for first timestep'
     else
        call t_startf('cnbalchk')
        call CBalanceCheck(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)
        call NBalanceCheck(begc, endc, filter(nc)%num_soilc, filter(nc)%soilc)
        call t_stopf('cnbalchk')
     end if
#endif

   ! CH4
#if (defined LCH4)
     call t_startf('ch4')
     call ch4 (begg, endg, begl, endl, begc, endc, begp, endp, filter(nc)%num_soilc, filter(nc)%soilc, &
               filter(nc)%num_lakec, filter(nc)%lakec, filter(nc)%num_soilp, filter(nc)%soilp)
     call t_stopf('ch4')
#endif

     ! ============================================================================
     ! Determine albedos for next time step
     ! ============================================================================

     if (doalb) then
        call t_startf('surfalb')

        ! Albedos for non-urban columns

        call SurfaceAlbedo(begg, endg, begc, endc, begp, endp, &
                           filter(nc)%num_nourbanc, filter(nc)%nourbanc, &
                           filter(nc)%num_nourbanp, filter(nc)%nourbanp, &
                           nextsw_cday, declinp1)

        call t_stopf('surfalb')

        ! Albedos for urban columns

        call t_startf('urbsurfalb')

        if (filter(nc)%num_urbanl > 0) then
           call UrbanAlbedo(nc, begl, endl, begc, endc, begp, endp,   &
                            filter(nc)%num_urbanl, filter(nc)%urbanl, &
                            filter(nc)%num_urbanc, filter(nc)%urbanc, &
                            filter(nc)%num_urbanp, filter(nc)%urbanp)
        end if

        call t_stopf('urbsurfalb')

     end if

  end do
!$OMP END PARALLEL DO

  ! ============================================================================
  ! Determine gridcell averaged properties to send to atm (l2as and l2af derived types)
  ! ============================================================================

  call t_startf('clm_map2gcell')
  call clm_map2gcell( )
  call t_stopf('clm_map2gcell')

  ! ============================================================================
  ! Determine fields to send to glc
  ! ============================================================================
  
  if (create_glacier_mec_landunit) then
     call t_startf('create_s2x')
     call create_clm_s2x(init=.false.)
     call t_stopf('create_s2x')
  end if
  

  ! ============================================================================
  ! Write global average diagnostics to standard output
  ! ============================================================================

  nstep = get_nstep()
  if (wrtdia) call mpi_barrier(mpicom,ier)
  call t_startf('wrtdiag')
  call write_diagnostic(wrtdia, nstep)
  call t_stopf('wrtdiag')

  ! ============================================================================
  ! Read initial snow and soil moisture data at each time step
  ! ============================================================================

  call t_startf('inicperp')
  if (is_perpetual()) then
     call inicPerp()
  end if
  call t_stopf('inicperp')

  ! ============================================================================
  ! Update accumulators
  ! ============================================================================

  call t_startf('accum')
  call updateAccFlds()
  call t_stopf('accum')

  ! ============================================================================
  ! Update history buffer
  ! ============================================================================

  call t_startf('hbuf')
  call hist_update_hbuf()
  call t_stopf('hbuf')

  ! ============================================================================
  ! Call dv (dynamic vegetation) at last time step of year
  ! NOTE: monp1, dayp1, and secp1 correspond to nstep+1
  ! ============================================================================

#if (defined CNDV)
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

!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
     do nc = 1,nclumps
        call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
        call dv(begg, endg, begp, endp,                  &
           filter(nc)%num_natvegp, filter(nc)%natvegp, kyr)
     end do
!$OMP END PARALLEL DO
  end if
  call t_stopf('d2dgvm')
#endif

  ! ============================================================================
  ! Create history and write history tapes if appropriate
  ! ============================================================================

  call t_startf('clm_drv_io')
#ifndef _NOIO

  call t_startf('clm_drv_io_htapes')
  call hist_htapes_wrapup( rstwr, nlend )
  call t_stopf('clm_drv_io_htapes')

  ! ============================================================================
  ! Write to CNDV history buffer if appropriate
  ! ============================================================================

#if (defined CNDV)
  if (monp1==1 .and. dayp1==1 .and. secp1==dtime .and. nstep>0)  then
     call t_startf('clm_drv_io_hdgvm')
     call histCNDV()
     if (masterproc) write(iulog,*) 'Annual CNDV calculations are complete'
     call t_stopf('clm_drv_io_hdgvm')
  end if
#endif

  ! ============================================================================
  ! Write restart/initial files if appropriate
  ! ============================================================================

  if (rstwr) then
     call t_startf('clm_drv_io_wrest')
     filer = restFile_filename(rdate=rdate)
     call restFile_write( filer, nlend, rdate=rdate )
     call t_stopf('clm_drv_io_wrest')
  end if

#endif
  call t_stopf('clm_drv_io')

end subroutine clm_drv

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diagnostic
!
! !INTERFACE:
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
  use spmdMod    , only : masterproc, npes, MPI_REAL8, MPI_ANY_SOURCE, &
                          MPI_STATUS_SIZE, mpicom, MPI_SUM
  use shr_sys_mod, only : shr_sys_flush
  use abortutils , only : endrun
!
! !ARGUMENTS:
  implicit none
  logical, intent(in) :: wrtdia     !true => write diagnostic
  integer, intent(in) :: nstep      !model time step
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: p                       ! loop index
  integer :: begp, endp              ! per-proc beginning and ending pft indices
  integer :: begc, endc              ! per-proc beginning and ending column indices
  integer :: begl, endl              ! per-proc beginning and ending landunit indices
  integer :: begg, endg              ! per-proc gridcell ending gridcell indices
  integer :: numg                    ! total number of gridcells across all processors
  integer :: numl                    ! total number of landunits across all processors
  integer :: numc                    ! total number of columns across all processors
  integer :: nump                    ! total number of pfts across all processors
  integer :: ier                     ! error status
  real(r8):: psum                    ! partial sum of ts
  real(r8):: tsum                    ! sum of ts
  real(r8):: tsxyav                  ! average ts for diagnostic output
  integer :: status(MPI_STATUS_SIZE) ! mpi status
  logical,parameter :: old_sendrecv = .false.  ! Flag if should use old send/receive method rather than MPI reduce
!------------------------------------------------------------------------

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
  call get_proc_global(numg, numl, numc, nump)

  if (wrtdia) then

     call t_barrierf('sync_write_diag', mpicom)
     psum = sum(clm_l2a%t_rad(begg:endg))
     if (old_sendrecv) then
        if (masterproc) then
           tsum = psum
           do p = 1, npes-1
              call mpi_recv(psum, 1, MPI_REAL8, p, 999, mpicom, status, ier)
              if (ier/=0) then
                 write(iulog,*) 'write_diagnostic: Error in mpi_recv()'
                 call endrun
              end if
              tsum = tsum + psum
           end do
        else
           call mpi_send(psum, 1, MPI_REAL8, 0, 999, mpicom, ier)
           if (ier/=0) then
              write(iulog,*) 'write_diagnostic: Error in mpi_send()'
              call endrun
           end if
        end if
     else
        call mpi_reduce(psum, tsum, 1, MPI_REAL8, MPI_SUM, 0, mpicom, ier)
        if (ier/=0) then
           write(iulog,*) 'write_diagnostic: Error in mpi_reduce()'
           call endrun
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
