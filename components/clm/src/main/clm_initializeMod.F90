module clm_initializeMod

  !-----------------------------------------------------------------------
  ! Performs land model initialization
  !
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use spmdMod          , only : masterproc, iam
  use shr_sys_mod      , only : shr_sys_flush
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use decompMod        , only : bounds_type, get_proc_bounds, get_proc_clumps, get_clump_bounds
  use abortutils       , only : endrun
  use clm_varctl       , only : nsrest, nsrStartup, nsrContinue, nsrBranch
  use clm_varctl       , only : create_glacier_mec_landunit, iulog
  use clm_varctl       , only : use_lch4, use_cn, use_cndv, use_voc, use_c13, use_c14, use_fates, use_betr  
  use clm_varsur       , only : wt_lunit, urban_valid, wt_nat_patch, wt_cft, wt_glc_mec, topo_glc_mec
  use clm_varsur       , only : fert_cft
  use perf_mod         , only : t_startf, t_stopf
  !use readParamsMod    , only : readParameters
  use readParamsMod    , only : readSharedParameters, readPrivateParameters
  use ncdio_pio        , only : file_desc_t
  use FatesInterfaceMod, only : set_fates_global_elements
  use BeTRSimulationALM, only : create_betr_simulation_alm
  ! 
  !-----------------------------------------
  ! Definition of component types
  !-----------------------------------------
  use GridcellType           , only : grc_pp
  use TopounitType           , only : top_pp, top_as, top_af, top_es
  use LandunitType           , only : lun_pp                
  use ColumnType             , only : col_pp                
  use VegetationType         , only : veg_pp                
  use clm_instMod
  use WaterBudgetMod         , only : WaterBudget_Reset
  use clm_varctl             , only : do_budgets
  !
  implicit none
  save
  !
  public :: initialize1  ! Phase one initialization
  public :: initialize2  ! Phase two initialization
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initialize1( )
    !
    ! !DESCRIPTION:
    ! CLM initialization first phase 
    !
    ! !USES:
    use clm_varpar                , only: clm_varpar_init, natpft_lb, natpft_ub, cft_lb, cft_ub, maxpatch_glcmec
    use clm_varcon                , only: clm_varcon_init
    use landunit_varcon           , only: landunit_varcon_init, max_lunit, istice_mec
    use column_varcon             , only: col_itype_to_icemec_class
    use clm_varctl                , only: fsurdat, fatmlndfrc, flndtopo, fglcmask, noland, version  
    use pftvarcon                 , only: pftconrd
    use soilorder_varcon          , only: soilorder_conrd
    use decompInitMod             , only: decompInit_lnd, decompInit_clumps, decompInit_gtlcp
    use domainMod                 , only: domain_check, ldomain, domain_init
    use surfrdMod                 , only: surfrd_get_globmask, surfrd_get_grid, surfrd_get_topo, surfrd_get_data
    use controlMod                , only: control_init, control_print, NLFilename
    use ncdio_pio                 , only: ncd_pio_init
    use initGridCellsMod          , only: initGridCells, initGhostGridCells
    use CH4varcon                 , only: CH4conrd
    use UrbanParamsType           , only: UrbanInput
    use CLMFatesParamInterfaceMod , only: FatesReadPFTs
    use surfrdMod                 , only: surfrd_get_grid_conn
    use clm_varctl                , only: lateral_connectivity, domain_decomp_type
    use decompInitMod             , only: decompInit_lnd_using_gp, decompInit_ghosts
    use domainLateralMod          , only: ldomain_lateral, domainlateral_init
    use SoilTemperatureMod        , only: init_soil_temperature
    use ExternalModelInterfaceMod , only: EMI_Determine_Active_EMs
    use dynSubgridControlMod      , only: dynSubgridControl_init
    use filterMod                 , only: allocFilters
    use reweightMod               , only: reweight_wrapup
    !
    ! !LOCAL VARIABLES:
    integer           :: ier                     ! error status
    integer           :: i,j,n,k,c,l,g           ! indices
    integer           :: nl                      ! gdc and glo lnd indices
    integer           :: ns, ni, nj              ! global grid sizes
    integer           :: begg, endg              ! processor bounds
    integer           :: icemec_class            ! current icemec class (1..maxpatch_glcmec)
    type(bounds_type) :: bounds_proc             
    type(bounds_type) :: bounds_clump            ! clump bounds
    integer ,pointer  :: amask(:)                ! global land mask
    integer ,pointer  :: cellsOnCell(:,:)        ! grid cell level connectivity
    integer ,pointer  :: edgesOnCell(:,:)        ! index to determine distance between neighbors from dcEdge
    integer ,pointer  :: nEdgesOnCell(:)         ! number of edges
    real(r8), pointer :: dcEdge(:)               ! distance between centroids of grid cells
    real(r8), pointer :: dvEdge(:)               ! distance between vertices
    real(r8), pointer :: areaCell(:)             ! area of grid cells [m^2]
    integer           :: nCells_loc              ! number of grid cell level connectivity saved locally
    integer           :: nEdges_loc              ! number of edge length saved locally
    integer           :: maxEdges                ! max number of edges/neighbors
    integer           :: nclumps                 ! number of clumps on this processor
    integer           :: nc                      ! clump index
    character(len=32) :: subname = 'initialize1' ! subroutine name
    !-----------------------------------------------------------------------

    call t_startf('clm_init1')

    ! ------------------------------------------------------------------------
    ! Initialize run control variables, timestep
    ! ------------------------------------------------------------------------

    if ( masterproc )then
       write(iulog,*) trim(version)
       write(iulog,*)
       write(iulog,*) 'Attempting to initialize the land model .....'
       write(iulog,*)
       call shr_sys_flush(iulog)
    endif

    call control_init()
    call clm_varpar_init()
    call clm_varcon_init()
    call landunit_varcon_init()
    call ncd_pio_init()
    call clm_petsc_init()
    call init_soil_temperature()

    if (masterproc) call control_print()

    call dynSubgridControl_init(NLFilename)

    ! ------------------------------------------------------------------------
    ! Read in global land grid and land mask (amask)- needed to set decomposition
    ! ------------------------------------------------------------------------

    ! global memory for amask is allocate in surfrd_get_glomask - must be
    ! deallocated below
    if (masterproc) then
       write(iulog,*) 'Attempting to read global land mask from ',trim(fatmlndfrc)
       call shr_sys_flush(iulog)
    endif
    call surfrd_get_globmask(filename=fatmlndfrc, mask=amask, ni=ni, nj=nj)

    ! Exit early if no valid land points
    if ( all(amask == 0) )then
       if (masterproc) write(iulog,*) trim(subname)//': no valid land points do NOT run clm'
       noland = .true.
       return
    end if


    ! ------------------------------------------------------------------------
    ! If specified, read the grid level connectivity
    ! ------------------------------------------------------------------------

    if (lateral_connectivity) then
       call surfrd_get_grid_conn(fatmlndfrc, cellsOnCell, edgesOnCell, &
            nEdgesOnCell, areaCell, dcEdge, dvEdge, &
            nCells_loc, nEdges_loc, maxEdges)
    else
       nullify(cellsOnCell)
       nCells_loc = 0
       maxEdges   = 0
    endif

    ! ------------------------------------------------------------------------
    ! Determine clm gridcell decomposition and processor bounds for gridcells
    ! ------------------------------------------------------------------------

    select case (trim(domain_decomp_type))
    case ("round_robin")
       call decompInit_lnd(ni, nj, amask)
       deallocate(amask)
    case ("graph_partitioning")
       call decompInit_lnd_using_gp(ni, nj, cellsOnCell, nCells_loc, maxEdges, amask)
    case default
       call endrun(msg='ERROR clm_initializeMod: '//&
            'Unsupported domain_decomp_type = ' // trim(domain_decomp_type))
    end select

    if (lateral_connectivity) then
       call domainlateral_init(ldomain_lateral, cellsOnCell, edgesOnCell, &
            nEdgesOnCell, areaCell, dcEdge, dvEdge, &
            nCells_loc, nEdges_loc, maxEdges)
    endif

    ! *** Get JUST gridcell processor bounds ***
    ! Remaining bounds (landunits, columns, patches) will be determined 
    ! after the call to decompInit_glcp - so get_proc_bounds is called
    ! twice and the gridcell information is just filled in twice

    call get_proc_bounds(begg, endg)

    ! ------------------------------------------------------------------------
    ! Get grid and land fraction (set ldomain)
    ! ------------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read ldomain from ',trim(fatmlndfrc)
       call shr_sys_flush(iulog)
    endif
    if (create_glacier_mec_landunit) then
       call surfrd_get_grid(begg, endg, ldomain, fatmlndfrc, fglcmask)
    else
       call surfrd_get_grid(begg, endg, ldomain, fatmlndfrc)
    endif
    if (masterproc) then
       call domain_check(ldomain)
    endif
    ldomain%mask = 1  !!! TODO - is this needed?

    ! Get topo if appropriate (set ldomain%topo)

    if (flndtopo /= " ") then
       if (masterproc) then
          write(iulog,*) 'Attempting to read atm topo from ',trim(flndtopo)
          call shr_sys_flush(iulog)
       endif
       call surfrd_get_topo(ldomain, flndtopo)  
    endif

    ! Initialize urban model input (initialize urbinp data structure)
    ! This needs to be called BEFORE the call to surfrd_get_data since
    ! that will call surfrd_get_special which in turn calls check_urban

    call UrbanInput(begg, endg, mode='initialize')

    ! Allocate surface grid dynamic memory (just gridcell bounds dependent)

    allocate (wt_lunit     (begg:endg, max_lunit           ))
    allocate (urban_valid  (begg:endg                      ))
    allocate (wt_nat_patch (begg:endg, natpft_lb:natpft_ub ))
    allocate (wt_cft       (begg:endg, cft_lb:cft_ub       ))
    allocate (fert_cft     (begg:endg, cft_lb:cft_ub       ))
    if (create_glacier_mec_landunit) then
       allocate (wt_glc_mec  (begg:endg, maxpatch_glcmec))
       allocate (topo_glc_mec(begg:endg, maxpatch_glcmec))
    else
       allocate (wt_glc_mec  (1,1))
       allocate (topo_glc_mec(1,1))
    endif

    ! Read list of Patches and their corresponding parameter values
    ! Independent of model resolution, Needs to stay before surfrd_get_data

    call pftconrd()
    call soilorder_conrd()

    ! Read in FATES parameter values early in the call sequence as well
    ! The PFT file, specifically, will dictate how many pfts are used
    ! in fates, and this will influence the amount of memory we
    ! request from the model, which is relevant in set_fates_global_elements()
    if (use_fates) then
       call FatesReadPFTs()
    end if

    ! Read surface dataset and set up subgrid weight arrays
    call surfrd_get_data(begg, endg, ldomain, fsurdat)

    ! ------------------------------------------------------------------------
    ! Ask Fates to evaluate its own dimensioning needs.
    ! This determines the total amount of space it requires in its largest
    ! dimension.  We are currently calling that the "cohort" dimension, but
    ! it is really a utility dimension that captures the models largest
    ! size need.
    ! Sets:
    ! fates_maxElementsPerPatch
    ! fates_maxElementsPerSite (where a site is roughly equivalent to a column)
    ! 
    ! (Note: fates_maxELementsPerSite is the critical variable used by CLM
    ! to allocate space)
    ! ------------------------------------------------------------------------
    call set_fates_global_elements(use_fates)
    

    ! ------------------------------------------------------------------------
    ! Determine decomposition of subgrid scale landunits, topounits, columns, patches
    ! ------------------------------------------------------------------------

    if (create_glacier_mec_landunit) then
       call decompInit_clumps(ldomain%glcmask)
       call decompInit_ghosts(ldomain%glcmask)
    else
       call decompInit_clumps()
       call decompInit_ghosts()
    endif

    ! *** Get ALL processor bounds - for gridcells, landunit, columns and patches ***

    call get_proc_bounds(bounds_proc)
    
    ! Allocate memory for subgrid data structures
    ! This is needed here BEFORE the following call to initGridcells
    ! Note that the assumption is made that none of the subgrid initialization
    ! can depend on other elements of the subgrid in the calls below

    ! Initialize the gridcell data types
    call grc_pp%Init (bounds_proc%begg_all, bounds_proc%endg_all)
    
    ! Initialize the topographic unit data types
    call top_pp%Init (bounds_proc%begt_all, bounds_proc%endt_all) ! topology and physical properties
    call top_as%Init (bounds_proc%begt_all, bounds_proc%endt_all) ! atmospheric state variables (forcings)
    call top_af%Init (bounds_proc%begt_all, bounds_proc%endt_all) ! atmospheric flux variables (forcings)
    call top_es%Init (bounds_proc%begt_all, bounds_proc%endt_all) ! energy state
    
    ! Initialize the landunit data types
    call lun_pp%Init (bounds_proc%begl_all, bounds_proc%endl_all)
    
    ! Initialize the column data types
    call col_pp%Init (bounds_proc%begc_all, bounds_proc%endc_all)
    
    ! Initialize the vegetation (PFT) data types
    call veg_pp%Init (bounds_proc%begp_all, bounds_proc%endp_all)
    
    ! Initialize the cohort data types (nothing here yet)
    ! ...to be added later...

    ! Determine the number of active external models.
    call EMI_Determine_Active_EMs()

    ! Build hierarchy and topological info for derived types
    ! This is needed here for the following call to decompInit_glcp

    call initGridCells()

    ! Set global seg maps for gridcells, landlunits, columns and patches

    if (create_glacier_mec_landunit) then
       call decompInit_gtlcp(ns, ni, nj, ldomain%glcmask)
    else
       call decompInit_gtlcp(ns, ni, nj)
    endif

    ! Set filters

    call t_startf('init_filters')
    call allocFilters()
    call t_stopf('init_filters')

    nclumps = get_proc_clumps()
    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)
       call reweight_wrapup(bounds_clump, &
            ldomain%glcmask(bounds_clump%begg:bounds_clump%endg)*1._r8)
    end do
    !$OMP END PARALLEL DO

    ! ------------------------------------------------------------------------
    ! Remainder of initialization1
    ! ------------------------------------------------------------------------

    ! Set CH4 Model Parameters from namelist.
    ! Need to do before initTimeConst so that it knows whether to 
    ! look for several optional parameters on surfdata file.

    if (use_lch4) then
       call CH4conrd()
    end if

    ! Deallocate surface grid dynamic memory for variables that aren't needed elsewhere.
    ! Some things are kept until the end of initialize2; urban_valid is kept through the
    ! end of the run for error checking.

    deallocate (wt_lunit, wt_cft, wt_glc_mec)

    call t_stopf('clm_init1')

    ! initialize glc_topo
    ! TODO - does this belong here?
    do c = bounds_proc%begc, bounds_proc%endc
       l = col_pp%landunit(c)
       g = col_pp%gridcell(c)

       if (lun_pp%itype(l) == istice_mec) then
          ! For ice_mec landunits, initialize glc_topo based on surface dataset; this
          ! will get overwritten in the run loop by values sent from CISM
          icemec_class = col_itype_to_icemec_class(col_pp%itype(c))
          col_pp%glc_topo(c) = topo_glc_mec(g, icemec_class)
       else
          ! For other landunits, arbitrarily initialize glc_topo to 0 m; for landunits
          ! where this matters, this will get overwritten in the run loop by values sent
          ! from CISM
          col_pp%glc_topo(c) = 0._r8
       end if
    end do

  end subroutine initialize1


  !-----------------------------------------------------------------------
  subroutine initialize2( )
    !
    ! !DESCRIPTION:
    ! CLM initialization - second phase
    !
    ! !USES:
    use shr_orb_mod           , only : shr_orb_decl
    use shr_scam_mod          , only : shr_scam_getCloseLatLon
    use seq_drydep_mod        , only : n_drydep, drydep_method, DD_XLND
    use clm_varpar            , only : nlevsno, numpft, crop_prog, nlevsoi,max_patch_per_col
    use clm_varcon            , only : h2osno_max, bdsno, spval
    use landunit_varcon       , only : istice, istice_mec, istsoil
    use clm_varctl            , only : finidat, finidat_interp_source, finidat_interp_dest, fsurdat
    use clm_varctl            , only : use_century_decomp, single_column, scmlat, scmlon, use_cn, use_fates
    use clm_varorb            , only : eccen, mvelpp, lambm0, obliqr
    use clm_time_manager      , only : get_step_size, get_curr_calday
    use clm_time_manager      , only : get_curr_date, get_nstep, advance_timestep 
    use clm_time_manager      , only : timemgr_init, timemgr_restart_io, timemgr_restart
    use controlMod            , only : nlfilename
    use decompMod             , only : get_proc_clumps, get_proc_bounds, get_clump_bounds, bounds_type
    use domainMod             , only : ldomain
    use initInterpMod         , only : initInterp
    use DaylengthMod          , only : InitDaylength, daylength
    use fileutils             , only : getfil
    use filterMod             , only : filter
    use dynSubgridDriverMod   , only : dynSubgrid_init
    use reweightMod           , only : reweight_wrapup
    use subgridWeightsMod     , only : init_subgrid_weights_mod
    use histFileMod           , only : hist_htapes_build, htapes_fieldlist
    use histFileMod           , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    use restFileMod           , only : restFile_getfile, restFile_open, restFile_close
    use restFileMod           , only : restFile_read, restFile_write 
    use accumulMod            , only : print_accum_fields 
    use ndepStreamMod         , only : ndep_init, ndep_interp
    use EcosystemDynMod     , only : EcosystemDynInit
    use pdepStreamMod         , only : pdep_init, pdep_interp
    use DecompCascadeBGCMod , only : init_decompcascade_bgc
    use DecompCascadeCNMod  , only : init_decompcascade_cn
    use CNDecompCascadeContype, only : init_decomp_cascade_constants
    use VegetationPropertiesType        , only : veg_vp 
    use SoilorderConType      , only : soilorderconInit 
    use LakeCon               , only : LakeConInit 
    use SatellitePhenologyMod , only : SatellitePhenologyInit, readAnnualVegetation, interpMonthlyVeg
    use SnowSnicarMod         , only : SnowAge_init, SnowOptics_init
    use initVerticalMod       , only : initVertical
    use lnd2atmMod            , only : lnd2atm_minimal
    use glc2lndMod            , only : glc2lnd_type
    use lnd2glcMod            , only : lnd2glc_type 
    use SoilWaterRetentionCurveFactoryMod   , only : create_soil_water_retention_curve
    use clm_varctl                          , only : use_clm_interface, use_pflotran
    use clm_interface_pflotranMod           , only : clm_pf_interface_init !, clm_pf_set_restart_stamp
    use tracer_varcon         , only : is_active_betr_bgc    
    use clm_time_manager      , only : is_restart
    use ALMbetrNLMod          , only : betr_namelist_buffer
    !
    ! !ARGUMENTS    
    implicit none
    !
    ! !LOCAL VARIABLES:
    integer               :: c,i,g,j,k,l,p! indices
    integer               :: yr           ! current year (0, ...)
    integer               :: mon          ! current month (1 -> 12)
    integer               :: day          ! current day (1 -> 31)
    integer               :: ncsec        ! current time of day [seconds]
    integer               :: nc           ! clump index
    integer               :: nclumps      ! number of clumps on this processor
    character(len=256)    :: fnamer       ! name of netcdf restart file 
    character(len=256)    :: pnamer       ! full pathname of netcdf restart file
    character(len=256)    :: locfn        ! local file name
    type(file_desc_t)     :: ncid         ! netcdf id
    real(r8)              :: dtime        ! time step increment (sec)
    integer               :: nstep        ! model time step
    real(r8)              :: calday       ! calendar day for nstep
    real(r8)              :: caldaym1     ! calendar day for nstep-1
    real(r8)              :: declin       ! solar declination angle in radians for nstep
    real(r8)              :: declinm1     ! solar declination angle in radians for nstep-1
    real(r8)              :: eccf         ! earth orbit eccentricity factor
    type(bounds_type)     :: bounds_proc  ! processor bounds
    type(bounds_type)     :: bounds_clump ! clump bounds
    logical               :: lexist
    integer               :: closelatidx,closelonidx
    real(r8)              :: closelat,closelon
    real(r8), allocatable :: h2osno_col(:)
    real(r8), allocatable :: snow_depth_col(:)
    real(r8)              :: max_decl      ! temporary, for calculation of max_dayl
    integer               :: begp, endp
    integer               :: begc, endc
    integer               :: begl, endl
    real(r8), pointer     :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    character(len=32)     :: subname = 'initialize2' 
    !----------------------------------------------------------------------

    call t_startf('clm_init2')

    if (do_budgets) call WaterBudget_Reset('all')

    ! ------------------------------------------------------------------------
    ! Determine processor bounds and clumps for this processor
    ! ------------------------------------------------------------------------

    call get_proc_bounds(bounds_proc)
    nclumps = get_proc_clumps()

    ! ------------------------------------------------------------------------
    ! Read in shared parameters files
    ! ------------------------------------------------------------------------

    call readSharedParameters()

    ! ------------------------------------------------------------------------
    ! Initialize time manager
    ! ------------------------------------------------------------------------

    if (nsrest == nsrStartup) then  
       call timemgr_init()
    else
       call restFile_getfile(file=fnamer, path=pnamer)
       call restFile_open( flag='read', file=fnamer, ncid=ncid )
       call timemgr_restart_io( ncid=ncid, flag='read' )
       call restFile_close( ncid=ncid )
       call timemgr_restart()
    end if

    ! ------------------------------------------------------------------------
    ! Initialize daylength from the previous time step (needed so prev_dayl can be set correctly)
    ! ------------------------------------------------------------------------

    call t_startf('init_orbd')

    calday = get_curr_calday()
    call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr, declin, eccf )

    dtime = get_step_size()
    caldaym1 = get_curr_calday(offset=-int(dtime))
    call shr_orb_decl( caldaym1, eccen, mvelpp, lambm0, obliqr, declinm1, eccf )

    call t_stopf('init_orbd')
    
    call InitDaylength(bounds_proc, declin=declin, declinm1=declinm1)
             
    ! Initialize maximum daylength, based on latitude and maximum declination
    ! maximum declination hardwired for present-day orbital parameters, 
    ! +/- 23.4667 degrees = +/- 0.409571 radians, use negative value for S. Hem

    do g = bounds_proc%begg,bounds_proc%endg
       max_decl = 0.409571
       if (grc_pp%lat(g) < 0._r8) max_decl = -max_decl
       grc_pp%max_dayl(g) = daylength(grc_pp%lat(g), max_decl)
    end do

    ! History file variables

    if (use_cn) then
       call hist_addfld1d (fname='DAYL',  units='s', &
            avgflag='A', long_name='daylength', &
            ptr_gcell=grc_pp%dayl, default='inactive')

       call hist_addfld1d (fname='PREV_DAYL', units='s', &
            avgflag='A', long_name='daylength from previous timestep', &
            ptr_gcell=grc_pp%prev_dayl, default='inactive')
    end if

    ! ------------------------------------------------------------------------
    ! Initialize component data structures 
    ! ------------------------------------------------------------------------

    ! Note: new logic is in place that sets all the history fields to spval so
    ! there is no guesswork in the initialization to nans of the allocated variables

    ! First put in history calls for subgrid data structures - these cannot appear in the
    ! module for the subgrid data definition due to circular dependencies that are introduced
    
    data2dptr => col_pp%dz(:,-nlevsno+1:0)
    col_pp%dz(bounds_proc%begc:bounds_proc%endc,:) = spval
    call hist_addfld2d (fname='SNO_Z', units='m', type2d='levsno',  &
         avgflag='A', long_name='Snow layer thicknesses', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    col_pp%zii(bounds_proc%begc:bounds_proc%endc) = spval
    call hist_addfld1d (fname='ZII', units='m', &
         avgflag='A', long_name='convective boundary height', &
         ptr_col=col_pp%zii, default='inactive')

    call clm_inst_biogeophys(bounds_proc)

    if(use_betr)then
      !allocate memory for betr simulator
      allocate(ep_betr, source=create_betr_simulation_alm())
      !set internal filters for betr
      call ep_betr%BeTRSetFilter(maxpft_per_col=max_patch_per_col, boffline=.false.)
      call ep_betr%InitOnline(bounds_proc, lun_pp, col_pp, veg_pp, waterstate_vars, betr_namelist_buffer, masterproc)
      is_active_betr_bgc = ep_betr%do_soibgc()
    else
      allocate(ep_betr, source=create_betr_simulation_alm())
    endif
    
    call SnowOptics_init( ) ! SNICAR optical parameters:

    call SnowAge_init( )    ! SNICAR aging   parameters:

    ! ------------------------------------------------------------------------
    ! Read in private parameters files, this should be preferred for mulitphysics
    ! implementation, jinyun Tang, Feb. 11, 2015
    ! ------------------------------------------------------------------------
    if(use_cn .or. use_fates) then
       call init_decomp_cascade_constants()
    endif
    !read bgc implementation specific parameters when needed
    call readPrivateParameters()

    if (use_cn .or. use_fates) then
       if (.not. is_active_betr_bgc)then
          if (use_century_decomp) then
           ! Note that init_decompcascade_bgc needs cnstate_vars to be initialized
             call init_decompcascade_bgc(bounds_proc, cnstate_vars, soilstate_vars)
          else 
           ! Note that init_decompcascade_cn needs cnstate_vars to be initialized
             call init_decompcascade_cn(bounds_proc, cnstate_vars)
          end if
       endif
    endif

    ! FATES is instantiated in the following call.  The global is in clm_inst
    call clm_inst_biogeochem(bounds_proc)

    ! ------------------------------------------------------------------------
    ! Initialize accumulated fields
    ! ------------------------------------------------------------------------

    ! The time manager needs to be initialized before thes called is made, since
    ! the step size is needed. 

    call t_startf('init_accflds')

    call atm2lnd_vars%initAccBuffer(bounds_proc)
    
    call top_as%InitAccBuffer(bounds_proc)
    
    call top_af%InitAccBuffer(bounds_proc)

    call temperature_vars%initAccBuffer(bounds_proc)

    call canopystate_vars%initAccBuffer(bounds_proc)

    if (use_cndv) then
       call dgvs_vars%initAccBuffer(bounds_proc)
    end if

    if (crop_prog) then
       call crop_vars%initAccBuffer(bounds_proc)
    end if

    call cnstate_vars%initAccBuffer(bounds_proc)

    call print_accum_fields()

    call t_stopf('init_accflds')

    ! ------------------------------------------------------------------------
    ! Initializate dynamic subgrid weights (for prescribed transient Patches, CNDV
    ! and/or dynamic landunits); note that these will be overwritten in a
    ! restart run
    ! ------------------------------------------------------------------------

    call t_startf('init_dyn_subgrid')
    call init_subgrid_weights_mod(bounds_proc)
    call dynSubgrid_init(bounds_proc, dgvs_vars, glc2lnd_vars, crop_vars)
    call t_stopf('init_dyn_subgrid')

    ! ------------------------------------------------------------------------
    ! Initialize modules (after time-manager initialization in most cases)
    ! ------------------------------------------------------------------------

    if (use_cn) then
       call EcosystemDynInit(bounds_proc)
    else
       call SatellitePhenologyInit(bounds_proc)
    end if

    if (use_cn .and. n_drydep > 0 .and. drydep_method == DD_XLND) then
       ! Must do this also when drydeposition is used so that estimates of monthly 
       ! differences in LAI can be computed
       call SatellitePhenologyInit(bounds_proc)
    end if

    ! ------------------------------------------------------------------------
    ! On restart only - process the history namelist. 
    ! ------------------------------------------------------------------------

    ! Later the namelist from the restart file will be used.  This allows basic
    ! checking to make sure you didn't try to change the history namelist on restart.

    if (nsrest == nsrContinue ) then
       call htapes_fieldlist()
    end if

    ! ------------------------------------------------------------------------
    ! Read restart/initial info 
    ! ------------------------------------------------------------------------

    if (nsrest == nsrStartup) then

       if (finidat == ' ') then
          if (finidat_interp_source == ' ') then
             if (masterproc) then
                write(iulog,*)'Using cold start initial conditions '
             end if
          else 
             if (masterproc) then
                write(iulog,*)'Interpolating initial conditions from ',trim(finidat_interp_source),&
                     ' and creating new initial conditions ', trim(finidat_interp_dest)
             end if
          end if
       else 
          if (masterproc) then
             write(iulog,*)'Reading initial conditions from ',trim(finidat)
          end if
          call getfil( finidat, fnamer, 0 )
          call restFile_read(bounds_proc, fnamer,                                             &
               atm2lnd_vars, aerosol_vars, canopystate_vars, cnstate_vars,                    &
               carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, carbonflux_vars, &
               ch4_vars, dgvs_vars, energyflux_vars, frictionvel_vars, lakestate_vars,        &
               nitrogenstate_vars, nitrogenflux_vars, photosyns_vars, soilhydrology_vars,     &
               soilstate_vars, solarabs_vars, surfalb_vars, temperature_vars,                 &
               waterflux_vars, waterstate_vars,                                               &
               phosphorusstate_vars,phosphorusflux_vars,                                      &
               ep_betr,                                                                       &
               alm_fates, glc2lnd_vars, crop_vars)
       end if

    else if ((nsrest == nsrContinue) .or. (nsrest == nsrBranch)) then
       if (masterproc) then
          write(iulog,*)'Reading restart file ',trim(fnamer)
       end if
       call restFile_read(bounds_proc, fnamer,                                             &
            atm2lnd_vars, aerosol_vars, canopystate_vars, cnstate_vars,                    &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, carbonflux_vars, &
            ch4_vars, dgvs_vars, energyflux_vars, frictionvel_vars, lakestate_vars,        &
            nitrogenstate_vars, nitrogenflux_vars, photosyns_vars, soilhydrology_vars,     &
            soilstate_vars, solarabs_vars, surfalb_vars, temperature_vars,                 &
            waterflux_vars, waterstate_vars,                                               &
            phosphorusstate_vars,phosphorusflux_vars,                                      &
            ep_betr,                                                                       &
            alm_fates, glc2lnd_vars, crop_vars)

    end if
       
    ! ------------------------------------------------------------------------
    ! If appropriate, create interpolated initial conditions
    ! ------------------------------------------------------------------------

    if (nsrest == nsrStartup .and. finidat_interp_source /= ' ') then

       ! Check that finidat is not cold start - abort if it is
       if (finidat /= ' ') then
          call endrun(msg='ERROR clm_initializeMod: '//&
               'finidat and finidat_interp_source cannot both be non-blank')
       end if

       !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
       do nc = 1, nclumps
          call get_clump_bounds(nc, bounds_clump)
          call reweight_wrapup(bounds_clump, &
               glc2lnd_vars%icemask_grc(bounds_clump%begg:bounds_clump%endg))
       end do
       !$OMP END PARALLEL DO
       if(use_betr)then
         call ep_betr%set_active(bounds_proc, col_pp)
       endif
       ! Create new template file using cold start
       call restFile_write(bounds_proc, finidat_interp_dest,                               &
            atm2lnd_vars, aerosol_vars, canopystate_vars, cnstate_vars,                    &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, carbonflux_vars, &
            ch4_vars, dgvs_vars, energyflux_vars, frictionvel_vars, lakestate_vars,        &
            nitrogenstate_vars, nitrogenflux_vars, photosyns_vars, soilhydrology_vars,     &
            soilstate_vars, solarabs_vars, surfalb_vars, temperature_vars,                 &
            waterflux_vars, waterstate_vars,                                               &
            phosphorusstate_vars,phosphorusflux_vars,                                      &
            ep_betr,                                                                       &
            alm_fates, crop_vars)

       ! Interpolate finidat onto new template file
       call getfil( finidat_interp_source, fnamer,  0 )
       call initInterp(filei=fnamer, fileo=finidat_interp_dest, bounds=bounds_proc)

       ! Read new interpolated conditions file back in
       call restFile_read(bounds_proc, finidat_interp_dest,                                &
            atm2lnd_vars, aerosol_vars, canopystate_vars, cnstate_vars,                    &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, carbonflux_vars, &
            ch4_vars, dgvs_vars, energyflux_vars, frictionvel_vars, lakestate_vars,        &
            nitrogenstate_vars, nitrogenflux_vars, photosyns_vars, soilhydrology_vars,     &
            soilstate_vars, solarabs_vars, surfalb_vars, temperature_vars,                 &
            waterflux_vars, waterstate_vars,                                               &
            phosphorusstate_vars,phosphorusflux_vars,                                      &
            ep_betr,                                                                       &
            alm_fates, glc2lnd_vars, crop_vars)

       ! Reset finidat to now be finidat_interp_dest 
       ! (to be compatible with routines still using finidat)
       finidat = trim(finidat_interp_dest)

    end if

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)
       call reweight_wrapup(bounds_clump, &
            glc2lnd_vars%icemask_grc(bounds_clump%begg:bounds_clump%endg))
    end do
    !$OMP END PARALLEL DO

    if(use_betr)then
      call ep_betr%set_active(bounds_proc, col_pp)
    endif 
    ! ------------------------------------------------------------------------
    ! Initialize nitrogen deposition
    ! ------------------------------------------------------------------------

    if (use_cn) then
       call t_startf('init_ndep')
       call ndep_init(bounds_proc)
       call ndep_interp(bounds_proc, atm2lnd_vars)
       call t_stopf('init_ndep')
    end if
    
    ! ------------------------------------------------------------------------
    ! Initialize phosphorus deposition
    ! ------------------------------------------------------------------------

    if (use_cn) then
       call t_startf('init_pdep')
       call pdep_init(bounds_proc)
       call pdep_interp(bounds_proc, atm2lnd_vars)
       call t_stopf('init_pdep')
    end if
 

    ! ------------------------------------------------------------------------
    ! Initialize active history fields. 
    ! ------------------------------------------------------------------------

    ! This is only done if not a restart run. If a restart run, then this 
    ! information has already been obtained from the restart data read above. 
    ! Note that routine hist_htapes_build needs time manager information,
    ! so this call must be made after the restart information has been read.

    if (nsrest /= nsrContinue) then
       call hist_htapes_build()
    end if

    ! ------------------------------------------------------------------------
    ! Initialize variables that are associated with accumulated fields.
    ! ------------------------------------------------------------------------

    ! The following is called for both initial and restart runs and must
    ! must be called after the restart file is read 

    call atm2lnd_vars%initAccVars(bounds_proc)
    call top_as%InitAccVars(bounds_proc)
    call top_af%InitAccVars(bounds_proc)
    call temperature_vars%initAccVars(bounds_proc)
    call canopystate_vars%initAccVars(bounds_proc)
    if (use_cndv) then
       call dgvs_vars%initAccVars(bounds_proc)
    end if
    if (crop_prog) then
       call crop_vars%initAccVars(bounds_proc)
    end if
    call cnstate_vars%initAccVars(bounds_proc)

    !------------------------------------------------------------       
    ! Read monthly vegetation
    !------------------------------------------------------------       

    ! Even if CN is on, and dry-deposition is active, read CLMSP annual vegetation 
    ! to get estimates of monthly LAI

    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       call readAnnualVegetation(bounds_proc, canopystate_vars)
       if (nsrest == nsrStartup .and. finidat /= ' ') then
          ! Call interpMonthlyVeg for dry-deposition so that mlaidiff will be calculated
          ! This needs to be done even if CN or CNDV is on!
          call interpMonthlyVeg(bounds_proc, canopystate_vars)
       end if
    end if

    !------------------------------------------------------------       
    ! Determine gridcell averaged properties to send to atm
    !------------------------------------------------------------       

    if (nsrest == nsrStartup) then
       call t_startf('init_map2gc')
       call lnd2atm_minimal(bounds_proc, &
            waterstate_vars, surfalb_vars, energyflux_vars, lnd2atm_vars)
       call t_stopf('init_map2gc')
    end if

    !------------------------------------------------------------       
    ! Initialize sno export state to send to glc
    !------------------------------------------------------------       

    if (create_glacier_mec_landunit) then  
       !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
       do nc = 1,nclumps
          call get_clump_bounds(nc, bounds_clump)

          call t_startf('init_lnd2glc')
          call lnd2glc_vars%update_lnd2glc(bounds_clump,       &
               filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,   &
               temperature_vars, waterflux_vars,               &
               init=.true.)
          call t_stopf('init_lnd2glc')
       end do
       !$OMP END PARALLEL DO
    end if

    !------------------------------------------------------------       
    ! Deallocate wt_nat_patch
    !------------------------------------------------------------       

    ! wt_nat_patch was allocated in initialize1, but needed to be kept around through
    ! initialize2 for some consistency checking; now it can be deallocated

    deallocate(wt_nat_patch)

    ! --------------------------------------------------------------
    ! Initialise the FATES model state structure cold-start
    ! --------------------------------------------------------------
   
    if ( use_fates .and. .not.is_restart() .and. finidat == ' ') then
       call alm_fates%init_coldstart(waterstate_vars,canopystate_vars, &
                                     soilstate_vars, frictionvel_vars)
    end if

    ! topo_glc_mec was allocated in initialize1, but needed to be kept around through
    ! initialize2 because it is used to initialize other variables; now it can be
    ! deallocated
    


    ! topo_glc_mec was allocated in initialize1, but needed to be kept around through
    ! initialize2 because it is used to initialize other variables; now it can be
    ! deallocated

    deallocate(topo_glc_mec)

    !------------------------------------------------------------
    ! initialize clm_bgc_interface_data_type
    call t_startf('init_clm_interface_data & pflotran')
    if (use_clm_interface) then
        call clm_interface_data%Init(bounds_proc)
        ! PFLOTRAN initialization
        if (use_pflotran) then
            call clm_pf_interface_init(bounds_proc)
        end if
    end if
    call t_stopf('init_clm_interface_data & pflotran')
    !------------------------------------------------------------

    !------------------------------------------------------------       
    ! Write log output for end of initialization
    !------------------------------------------------------------       

    call t_startf('init_wlog')
    if (masterproc) then
       write(iulog,*) 'Successfully initialized the land model'
       if (nsrest == nsrStartup) then
          write(iulog,*) 'begin initial run at: '
       else
          write(iulog,*) 'begin continuation run at:'
       end if
       call get_curr_date(yr, mon, day, ncsec)
       write(iulog,*) '   nstep= ',get_nstep(), ' year= ',yr,' month= ',mon,&
            ' day= ',day,' seconds= ',ncsec
       write(iulog,*)
       write(iulog,'(72a1)') ("*",i=1,60)
       write(iulog,*)
    endif
    call t_stopf('init_wlog')

    call t_stopf('clm_init2')

  end subroutine initialize2

  !-----------------------------------------------------------------------
  subroutine initialize3( )
    !
    ! !DESCRIPTION:
    ! CLM initialization - third phase
    !
    ! !USES:
    use clm_varpar               , only : nlevsoi, nlevgrnd, nlevsno, max_patch_per_col
    use landunit_varcon          , only : istsoil, istcrop, istice_mec, istice_mec
    use landunit_varcon          , only : istice, istdlak, istwet, max_lunit
    use column_varcon            , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv
    use clm_varctl               , only : use_vsfm, vsfm_use_dynamic_linesearch
    use clm_varctl               , only : vsfm_include_seepage_bc, vsfm_satfunc_type
    use clm_varctl               , only : vsfm_lateral_model_type
    use clm_varctl               , only : use_petsc_thermal_model
    use clm_varctl               , only : lateral_connectivity
    use clm_varctl               , only : finidat
    use decompMod                , only : get_proc_clumps
    use mpp_varpar               , only : mpp_varpar_init
    use mpp_varcon               , only : mpp_varcon_init_landunit
    use mpp_varcon               , only : mpp_varcon_init_column
    use mpp_varctl               , only : mpp_varctl_init_vsfm
    use mpp_varctl               , only : mpp_varctl_init_petsc_thermal
    use mpp_bounds               , only : mpp_bounds_init_proc_bounds
    use mpp_bounds               , only : mpp_bounds_init_clump
    use ExternalModelInterfaceMod, only : EMI_Init_EM
    use ExternalModelConstants   , only : EM_ID_VSFM
    use ExternalModelConstants   , only : EM_ID_PTM

    implicit none

    type(bounds_type) :: bounds_proc
    logical           :: restart_vsfm          ! does VSFM need to be restarted

    call t_startf('clm_init3')

    ! Is this a restart run?
    restart_vsfm = .false.
    if (nsrest == nsrStartup) then
       if (finidat == ' ') then
          restart_vsfm = .false.
       else
          restart_vsfm = .true.
       end if
    else if ((nsrest == nsrContinue) .or. (nsrest == nsrBranch)) then
       restart_vsfm = .true.
    end if

    call mpp_varpar_init (nlevsoi, nlevgrnd, nlevsno, max_patch_per_col)

    call mpp_varcon_init_landunit   (istsoil, istcrop, istice, istice_mec, &
           istdlak, istwet, max_lunit)

    call mpp_varcon_init_column(icol_roof, icol_sunwall, icol_shadewall, &
      icol_road_imperv, icol_road_perv)

    call mpp_varctl_init_vsfm(use_vsfm, vsfm_use_dynamic_linesearch, &
      vsfm_include_seepage_bc, lateral_connectivity, restart_vsfm, &
      vsfm_satfunc_type, vsfm_lateral_model_type)

    call mpp_varctl_init_petsc_thermal(use_petsc_thermal_model)

    call get_proc_bounds(bounds_proc)
    call mpp_bounds_init_proc_bounds(bounds_proc%begg    , bounds_proc%endg,     &
                                     bounds_proc%begg_all, bounds_proc%endg_all, &
                                     bounds_proc%begc    , bounds_proc%endc,     &
                                     bounds_proc%begc_all, bounds_proc%endc_all)

    call mpp_bounds_init_clump(get_proc_clumps())

    if (use_vsfm) then
       call EMI_Init_EM(EM_ID_VSFM)
    endif

    if (use_petsc_thermal_model) then
       call EMI_Init_EM(EM_ID_PTM)
    endif

    call t_stopf('clm_init3')


  end subroutine initialize3

  !-----------------------------------------------------------------------
  subroutine clm_petsc_init()
    !
    ! !DESCRIPTION:
    ! Initialize PETSc
    !
#ifdef USE_PETSC_LIB
#include <petsc/finclude/petsc.h>
#endif
    ! !USES:
    use spmdMod    , only : mpicom
    use clm_varctl , only : use_vsfm
    use clm_varctl , only : lateral_connectivity
    use clm_varctl , only : use_petsc_thermal_model
#ifdef USE_PETSC_LIB
    use petscsys
#endif
    !
    implicit none
    !
    ! !LOCAL VARIABLES:
#ifdef USE_PETSC_LIB
    PetscErrorCode        :: ierr                  ! get error code from PETSc
#endif

    if ( (.not. use_vsfm)               .and. &
         (.not. lateral_connectivity)   .and. &
         (.not. use_petsc_thermal_model) ) return

#ifdef USE_PETSC_LIB
    ! Initialize PETSc
    PETSC_COMM_WORLD = mpicom
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRQ(ierr)

    PETSC_COMM_SELF  = MPI_COMM_SELF
    PETSC_COMM_WORLD = mpicom
#else
    call endrun(msg='ERROR clm_petsc_init: '//&
         'PETSc required but the code was not compiled using -DUSE_PETSC_LIB')
#endif

  end subroutine clm_petsc_init


end module clm_initializeMod
