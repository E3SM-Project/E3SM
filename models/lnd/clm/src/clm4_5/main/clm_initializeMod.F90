module clm_initializeMod

  !-----------------------------------------------------------------------
  ! Performs land model initialization
  !
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use spmdMod         , only : masterproc
  use shr_sys_mod     , only : shr_sys_flush
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use abortutils      , only : endrun
  use clm_varctl      , only : nsrest, nsrStartup, nsrContinue, nsrBranch
  use clm_varctl      , only : create_glacier_mec_landunit, iulog, use_lch4, use_cn, use_cndv, use_voc
  use clm_varsur      , only : wt_lunit, urban_valid, wt_nat_pft, wt_cft, wt_glc_mec, topo_glc_mec
  use perf_mod        , only : t_startf, t_stopf
  use readParamsMod   , only : readParameters
  use ncdio_pio
  use mct_mod

  implicit none
  save
  private    ! By default everything is private

  public :: initialize1  ! Phase one initialization
  public :: initialize2  ! Phase two initialization
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initialize1( )
    !
    ! !DESCRIPTION:
    ! CLM initialization1
    !
    ! !USES:
    use clmtypeInitMod   ,only: initClmtype
    use clm_varpar       ,only: clm_varpar_init, natpft_lb, natpft_ub, cft_lb, cft_ub, maxpatch_glcmec
    use clm_varcon       ,only: clm_varcon_init, max_lunit
    use clm_varctl       ,only: fsurdat, fatmlndfrc, flndtopo, fglcmask, noland, version  
    use pftvarcon        ,only: pftconrd
    use decompInitMod    ,only: decompInit_lnd, decompInit_clumps, decompInit_glcp
    use decompMod        ,only: bounds_type, get_proc_bounds 
    use domainMod        ,only: domain_check, ldomain, domain_init
    use surfrdMod        ,only: surfrd_get_globmask, surfrd_get_grid, surfrd_get_topo, surfrd_get_data 
    use controlMod       ,only: control_init, control_print, nlfilename
    use UrbanInputMod    ,only: UrbanInput
    use ncdio_pio        ,only: ncd_pio_init
    use clm_atmlnd       ,only: init_atm2lnd, init_lnd2atm
    use clm_glclnd       ,only: init_glc2lnd_type, init_lnd2glc_type, clm_x2s, clm_s2x
    use initGridCellsMod ,only: initGridCells
    use ch4varcon        ,only: ch4conrd
    !
    ! !LOCAL VARIABLES:
    integer  :: ier                              ! error status
    integer  :: i,j,n,k                          ! loop indices
    integer  :: nl                               ! gdc and glo lnd indices
    integer  :: ns, ni, nj                       ! global grid sizes
    integer  :: begg, endg                       ! processor bounds
    type(bounds_type) :: bounds_proc             ! bounds
    integer ,pointer  :: amask(:)                ! global land mask
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
    call ncd_pio_init()

    if (masterproc) call control_print()

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
    ! Determine clm gridcell decomposition and processor bounds for gridcells
    ! ------------------------------------------------------------------------

    call decompInit_lnd(ni, nj, amask)
    deallocate(amask)

    ! *** Get JUST gridcell processor bounds ***
    ! Remaining bounds (landunits, columns, pfts) will be determined 
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

    call UrbanInput(begg, endg, mode='initialize')

    ! Allocate surface grid dynamic memory (just gridcell bounds dependent)

    allocate (wt_lunit   (begg:endg, max_lunit))
    allocate (urban_valid(begg:endg))
    allocate (wt_nat_pft (begg:endg, natpft_lb:natpft_ub))
    allocate (wt_cft     (begg:endg, cft_lb:cft_ub))
    if (create_glacier_mec_landunit) then
       allocate (wt_glc_mec  (begg:endg, maxpatch_glcmec))
       allocate (topo_glc_mec(begg:endg, maxpatch_glcmec))
    else
       allocate (wt_glc_mec  (1,1))
       allocate (topo_glc_mec(1,1))
    endif

    ! Read list of PFTs and their corresponding parameter values
    ! Independent of model resolution, Needs to stay before surfrd_get_data

    call pftconrd()

    ! Read surface dataset and set up subgrid weight arrays
    
    call surfrd_get_data(begg, endg, ldomain, fsurdat)

    ! ------------------------------------------------------------------------
    ! Determine decomposition of subgrid scale landunits, columns, pfts
    ! ------------------------------------------------------------------------

    if (create_glacier_mec_landunit) then
       call decompInit_clumps(ns, ni, nj, ldomain%glcmask)
    else
       call decompInit_clumps(ns, ni, nj)
    endif

    ! *** Get ALL processor bounds - for gridcells, landunit, columns and pfts ***

    call get_proc_bounds(bounds_proc)
    
    ! Allocate memory and initialize values of clmtype data structures
    ! This is needed here for the following call to initGridcells

    call initClmtype(bounds_proc)

    ! Build hierarchy and topological info for derived types
    ! This is needed here for the following call to decompInit_glcp

    call initGridCells()

    ! Set global seg maps for gridcells, landlunits, columns and pfts

    if (create_glacier_mec_landunit) then
       call decompInit_glcp(ns, ni, nj, ldomain%glcmask)
    else
       call decompInit_glcp(ns, ni, nj)
    endif

    ! Initialize atm->lnd, lnd->atm, glc->lnd and lnd->glc data structures

    call init_atm2lnd(bounds_proc)
    call init_lnd2atm(bounds_proc)
    if (create_glacier_mec_landunit) then
       call init_glc2lnd_type(bounds_proc, clm_x2s)
       call init_lnd2glc_type(bounds_proc, clm_s2x)
    endif

    ! ------------------------------------------------------------------------
    ! Remainder of initialization1
    ! ------------------------------------------------------------------------

    ! Set CH4 Model Parameters from namelist.
    ! Need to do before iniTimeConst so that it knows whether to 
    ! look for several optional parameters on surfdata file.

    if (use_lch4) then
       call ch4conrd()
    end if

    ! Deallocate surface grid dynamic memory for variables that aren't needed elsewhere.
    ! Some things are kept until the end of initialize2; urban_valid is kept through the
    ! end of the run for error checking.

    deallocate (wt_lunit, wt_cft, wt_glc_mec)

    call t_stopf('clm_init1')

  end subroutine initialize1


  !-----------------------------------------------------------------------
  subroutine initialize2( )
    !
    ! !DESCRIPTION:
    ! CLM initialization2
    !
    ! !USES:
    use shr_orb_mod           , only : shr_orb_decl
    use seq_drydep_mod        , only : n_drydep, drydep_method, DD_XLND
    use clm_atmlnd            , only : clm_map2gcell_minimal
    use clm_glclnd            , only : update_clm_s2x
    use clm_glclnd            , only : init_glc2lnd_type, init_lnd2glc_type, clm_x2s, clm_s2x
    use clm_varctl            , only : finidat, finidat_interp_source, finidat_interp_dest
    use clm_varorb            , only : eccen, mvelpp, lambm0, obliqr
    use clm_time_manager      , only : get_step_size, get_curr_calday
    use clm_time_manager      , only : get_curr_date, get_nstep, advance_timestep 
    use clm_time_manager      , only : timemgr_init, timemgr_restart_io, timemgr_restart
    use decompMod             , only : get_proc_clumps, get_proc_bounds, get_clump_bounds, bounds_type
    use filterMod             , only : allocFilters, filter
    use reweightMod           , only : reweight_wrapup
    use subgridWeightsMod     , only : init_subgrid_weights_mod
    use histFldsMod           , only : hist_initFlds
    use histFileMod           , only : hist_htapes_build, htapes_fieldlist
    use restFileMod           , only : restFile_getfile, restFile_open, restFile_close
    use restFileMod           , only : restFile_read, restFile_write 
    use accFldsMod            , only : initAccFlds, initAccClmtype
    use ndepStreamMod         , only : ndep_init, ndep_interp
    use dynSubgridDriverMod   , only : dynSubgrid_init
    use CNEcosystemDynMod     , only : CNEcosystemDynInit
    use SatellitePhenologyMod , only : SatellitePhenologyInit, readAnnualVegetation, interpMonthlyVeg
    use DustMod               , only : Dustini
    use VOCEmissionMod        , only : VOCEmission_init
    use initTimeConstMod      , only : initTimeConst
    use UrbanInitMod          , only : initTimeConstUrban
    use SLakeInitMod          , only : initTimeConstSLake
    use initColdMod           , only : initCold
    use initInterpMod         , only : initInterp
    use DaylengthMod          , only : InitDaylength
    use BiogeophysRestMod     , only : bound_h2osoi
    use fileutils             , only : getfil
    !
    ! !ARGUMENTS    
    implicit none
    !
    ! !LOCAL VARIABLES:
    integer            :: i,j,k        ! indices
    integer            :: yr           ! current year (0, ...)
    integer            :: mon          ! current month (1 -> 12)
    integer            :: day          ! current day (1 -> 31)
    integer            :: ncsec        ! current time of day [seconds]
    integer            :: nc           ! clump index
    integer            :: nclumps      ! number of clumps on this processor
    character(len=256) :: fnamer       ! name of netcdf restart file 
    character(len=256) :: pnamer       ! full pathname of netcdf restart file
    type(file_desc_t)  :: ncid         ! netcdf id
    real(r8)           :: dtime        ! time step increment (sec)
    integer            :: nstep        ! model time step
    real(r8)           :: calday       ! calendar day for nstep
    real(r8)           :: caldaym1     ! calendar day for nstep-1
    real(r8)           :: declin       ! solar declination angle in radians for nstep
    real(r8)           :: declinm1     ! solar declination angle in radians for nstep-1
    real(r8)           :: eccf         ! earth orbit eccentricity factor
    type(bounds_type)  :: bounds_proc  ! processor bounds
    type(bounds_type)  :: bounds_clump ! clump bounds
    logical            :: lexist
    character(len=32)  :: subname = 'initialize2' 
    !----------------------------------------------------------------------

    call t_startf('clm_init2')

    ! ------------------------------------------------------------------------
    ! Determine processor bounds and clumps for this processor
    ! ------------------------------------------------------------------------

    call get_proc_bounds(bounds_proc)
    nclumps = get_proc_clumps()

    ! ------------------------------------------------------------------------
    ! Initialize CLMSP ecosystem dynamics 
    ! ------------------------------------------------------------------------

    ! Must do this also when drydeposition is used so that estimates of monthly 
    ! differences in LAI can be computed

    call t_startf('init_sphen')
    if ((.not. use_cn) .or. ((use_cn) .and. (n_drydep > 0 .and. drydep_method == DD_XLND ))) then
       call SatellitePhenologyInit(bounds_proc)
    end if
    call t_stopf('init_sphen')

    ! ------------------------------------------------------------------------
    ! Read in parameters files
    ! ------------------------------------------------------------------------

    call readParameters()

    ! ------------------------------------------------------------------------
    ! Initialize dust emissions model 
    ! ------------------------------------------------------------------------

    call Dustini(bounds_proc)
   
    ! ------------------------------------------------------------------------
    ! Initialize MEGAN emissions model 
    ! ------------------------------------------------------------------------

    if (use_voc) then
       call VOCEmission_init()
    end if

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
             
    ! ------------------------------------------------------------------------
    ! Initialize CN Ecosystem Dynamics (must be after time-manager initialization)
    ! ------------------------------------------------------------------------

    if (use_cn) then 
       call CNEcosystemDynInit(bounds_proc)
    end if

    ! ------------------------------------------------------------------------
    ! Initialize accumulated fields
    ! ------------------------------------------------------------------------

    ! The time manager needs to be initialized before this called is made, since
    ! the step size is needed.

    call t_startf('init_accflds')
    call initAccFlds(bounds_proc)
    call t_stopf('init_accflds')

    ! ------------------------------------------------------------------------
    ! Initialization of dynamic subgrid weights (for prescribed transient PFTs,
    ! CNDV, and/or dynamic landunits); note that these will be overwritten in a
    ! restart run
    ! ------------------------------------------------------------------------

    call t_startf('init_dyn_subgrid')
    call init_subgrid_weights_mod(bounds_proc)
    call dynSubgrid_init(bounds_proc)
    call t_stopf('init_dyn_subgrid')

    ! ------------------------------------------------------------------------
    ! Initialize time constant variables
    ! Must be done BEFORE calling initCold()
    ! ------------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*)'Initializing time constant variables '
    end if

    ! Note that initTimeConst depends on variables that are calculated in 
    ! initTimeConstUrban - so initTimeConstUrban must be called first
    call initTimeConstUrban(bounds_proc)
    call initTimeConst(bounds_proc) 
    call initTimeConstSlake(bounds_proc)

    ! ------------------------------------------------------------------------
    ! Create cold start initial conditions (will be overwritten by 
    ! initial/restart file read if approprate
    ! **** NOTE - this is always called****
    ! ------------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*)'Creating cold start initial conditions '
    end if

    call initCold(bounds_proc)

    ! ------------------------------------------------------------------------
    ! Initialize master history list. 
    ! ------------------------------------------------------------------------

    call t_startf('hist_initFlds')
    call hist_initFlds()

    ! ------------------------------------------------------------------------
    ! On restart only - process the history namelist. 
    ! ------------------------------------------------------------------------

    ! Later the namelist from the restart file will be used.  This allows basic
    ! checking to make sure you didn't try to change the history namelist on restart.

    if (nsrest == nsrContinue ) then
       call htapes_fieldlist()
    end if
    call t_stopf('hist_initFlds')

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
          call restFile_read(bounds_proc, fnamer)
       end if
    else if ((nsrest == nsrContinue) .or. (nsrest == nsrBranch)) then
       if (masterproc) then
          write(iulog,*)'Reading restart file ',trim(fnamer)
       end if
       call restFile_read(bounds_proc, fnamer)
    end if

    ! ------------------------------------------------------------------------
    ! Initialize filters and weights
    ! ------------------------------------------------------------------------
    
    call t_startf('init_filters')
    call allocFilters()
    call t_stopf('init_filters')

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
          call reweight_wrapup(bounds_clump)
       end do
       !$OMP END PARALLEL DO

       ! Create new template file using cold start
       call restFile_write(bounds_proc, finidat_interp_dest)

       ! Interpolate finidat onto new template file
       call getfil( finidat_interp_source, fnamer,  0 )
       call initInterp(filei=fnamer, fileo=finidat_interp_dest, bounds=bounds_proc)

       ! Read new interpolated conditions file back in
       call restFile_read(bounds_proc, finidat_interp_dest)

       ! Reset finidat to now be finidat_interp_dest 
       ! (to be compatible with routines still using finidat)
       finidat = trim(finidat_interp_dest)

    end if

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)
       call reweight_wrapup(bounds_clump)
    end do
    !$OMP END PARALLEL DO

    ! ------------------------------------------------------------------------
    ! Initialize nitrogen deposition
    ! ------------------------------------------------------------------------

    if (use_cn) then
       call t_startf('init_ndep')
       call ndep_init(bounds_proc)
       call ndep_interp(bounds_proc)
       call t_stopf('init_ndep')
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
    ! Initialize clmtype variables that are obtained from accumulated fields.
    ! ------------------------------------------------------------------------

    ! This routine is called for both initial and restart runs and must
    ! must be called after the restart file is read 

    call initAccClmtype(bounds_proc)

    !------------------------------------------------------------       
    ! Read monthly vegetation
    !------------------------------------------------------------       

    ! Even if CN is on, and dry-deposition is active, read CLMSP annual vegetation 
    ! to get estimates of monthly LAI

    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       call readAnnualVegetation(bounds_proc)
       if (nsrest == nsrStartup .and. finidat /= ' ') then
          ! Call interpMonthlyVeg for dry-deposition so that mlaidiff will be calculated
          ! This needs to be done even if CN or CNDV is on!
          call interpMonthlyVeg(bounds_proc)
       end if
    end if

    !------------------------------------------------------------       
    ! Determine gridcell averaged properties to send to atm
    !------------------------------------------------------------       

    if (nsrest == nsrStartup) then
       call t_startf('init_map2gc')
       call clm_map2gcell_minimal(bounds_proc)
       call t_stopf('init_map2gc')
    end if

    !------------------------------------------------------------       
    ! Initialize sno export state to send to glc
    !------------------------------------------------------------       

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       if (create_glacier_mec_landunit) then  
          call t_startf('create_s2x')
          call update_clm_s2x(bounds_clump, &
               filter(nc)%num_do_smb_c, filter(nc)%do_smb_c, &
               init=.true.)
          call t_stopf('create_s2x')
       end if
    end do
    !$OMP END PARALLEL DO

    !------------------------------------------------------------       
    ! Deallocate wt_nat_pft
    !------------------------------------------------------------       

    ! wt_nat_pft was allocated in initialize1, but needed to be kept around through
    ! initialize2 for some consistency checking; now it can be deallocated

    deallocate(wt_nat_pft)

    ! topo_glc_mec was allocated in initialize1, but needed to be kept around through
    ! initialize2 because it is used to initialize other variables; now it can be
    ! deallocated

    deallocate(topo_glc_mec)

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

end module clm_initializeMod
