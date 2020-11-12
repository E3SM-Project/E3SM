module controlMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module which initializes run control variables. The following possible
  ! namelist variables are set default values and possibly read in on startup
  !
  ! Note: For definitions of namelist variables see
  !       ../../bld/namelist_files/namelist_definition.xml
  !       Display the file in a browser to see it neatly formatted in html.
  !
  ! !USES:
  use elm_varctl   
  use shr_kind_mod            , only: r8 => shr_kind_r8, SHR_KIND_CL
  use shr_nl_mod              , only: shr_nl_find_group_name
  use shr_const_mod           , only: SHR_CONST_CDAY
  use shr_log_mod             , only: errMsg => shr_log_errMsg
  use abortutils              , only: endrun
  use spmdMod                 , only: masterproc
  use decompMod               , only: clump_pproc
  use elm_varpar              , only: maxpatch_pft, maxpatch_glcmec, more_vertlayers
  use histFileMod             , only: max_tapes, max_namlen 
  use histFileMod             , only: hist_empty_htapes, hist_dov2xy, hist_avgflag_pertape, hist_type1d_pertape 
  use histFileMod             , only: hist_nhtfrq, hist_ndens, hist_mfilt, hist_fincl1, hist_fincl2, hist_fincl3
  use histFileMod             , only: hist_fincl4, hist_fincl5, hist_fincl6, hist_fexcl1, hist_fexcl2, hist_fexcl3
  use histFileMod             , only: hist_fexcl4, hist_fexcl5, hist_fexcl6
  use LakeCon                 , only: deepmixing_depthcrit, deepmixing_mixfact 
  use AllocationMod         , only: suplnitro
  use AllocationMod         , only: suplphos
  use ColumnDataType          , only: nfix_timeconst
  use CNNitrifDenitrifMod     , only: no_frozen_nitrif_denitrif
  use C14DecayMod           , only: use_c14_bombspike, atm_c14_filename
  use SoilLittVertTranspMod , only: som_adv_flux, max_depth_cryoturb
  use VerticalProfileMod    , only: exponential_rooting_profile, rootprof_exp 
  use VerticalProfileMod    , only: surfprof_exp, pftspecific_rootingprofile  
  use SharedParamsMod       , only: anoxia_wtsat
  use CanopyfluxesMod         , only: perchroot, perchroot_alt
  use CanopyHydrologyMod      , only: CanopyHydrology_readnl
  use SurfaceAlbedoMod        , only: albice, lake_melt_icealb
  use UrbanParamsType         , only: urban_hac, urban_traffic
  use elm_varcon              , only: h2osno_max
  use elm_varctl              , only: use_dynroot
  use AllocationMod         , only: nu_com_phosphatase,nu_com_nfix 
  use elm_varctl              , only: nu_com, use_var_soil_thick
  use seq_drydep_mod          , only: drydep_method, DD_XLND, n_drydep
  use elm_varctl              , only: forest_fert_exp
  use elm_varctl              , only: ECA_Pconst_RGspin
  use elm_varctl              , only: NFIX_PTASE_plant
  use elm_varctl              , only : use_pheno_flux_limiter
  use elm_varctl              , only: startdate_add_temperature, startdate_add_co2
  use elm_varctl              , only: add_temperature, add_co2
  use elm_varctl              , only: const_climate_hist
 !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: control_setNL ! Set namelist filename
  public :: control_init  ! initial run control information
  public :: control_print ! print run control information

  !
  !
  ! !PRIVATE TYPES:
  character(len=  7) :: runtyp(4)                        ! run type
  character(len=SHR_KIND_CL) :: NLFilename = 'lnd.stdin' ! Namelist filename

#if (defined _OPENMP)
  integer, external :: omp_get_max_threads  ! max number of threads that can execute concurrently in a single parallel region
#endif
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine control_setNL( NLfile )
    !
    ! !DESCRIPTION:
    ! Set the namelist filename to use
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(IN) :: NLFile ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'control_setNL'  ! subroutine name
    logical :: lexist                               ! File exists
    !------------------------------------------------------------------------

    ! Error checking...
    if ( len_trim(NLFile) == 0 )then
       call endrun(msg=' error: nlfilename entered is not set'//errMsg(__FILE__, __LINE__))
    end if
    inquire (file = trim(NLFile), exist = lexist)
    if ( .not. lexist )then
       call endrun(msg=' error: NLfilename entered does NOT exist:'//&
            trim(NLFile)//errMsg(__FILE__, __LINE__))
    end if
    if ( len_trim(NLFile) > len(NLFilename) )then
       call endrun(msg=' error: entered NLFile is too long'//errMsg(__FILE__, __LINE__))
    end if
    ! Set the filename
    NLFilename = NLFile
    NLFilename_in = NLFilename   ! For use in external namelists and to avoid creating dependencies on controlMod
  end subroutine control_setNL

  !------------------------------------------------------------------------
  subroutine control_init( )
    !
    ! !DESCRIPTION:
    ! Initialize CLM run control information
    !
    ! !USES:
    use clm_time_manager          , only : set_timemgr_init, get_timemgr_defaults
    use fileutils                 , only : getavu, relavu
    use shr_string_mod            , only : shr_string_getParentDir
    use elm_interface_pflotranMod , only : elm_pf_readnl
    use ALMBeTRNLMod              , only : betr_readNL
    !
    implicit none
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: starttype ! infodata start type
    integer :: i,j,n                ! loop indices
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    integer :: dtime                ! Integer time-step
    integer :: override_nsrest      ! If want to override the startup type sent from driver
    character(len=32) :: subname = 'control_init'  ! subroutine name
    !------------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Namelist Variables
    ! ----------------------------------------------------------------------

    ! Time step
    namelist / elm_inparm/ &
         dtime

    ! CLM namelist settings

    namelist /elm_inparm / &
         fatmlndfrc, finidat, nrevsn, &
         finidat_interp_source, finidat_interp_dest

    ! Input datasets

    namelist /elm_inparm/  &
         fsurdat, fatmtopo, flndtopo, &
         paramfile, fsnowoptics, fsnowaging,fsoilordercon


    ! History, restart options

    namelist /elm_inparm/  &
         hist_empty_htapes, hist_dov2xy, &
         hist_avgflag_pertape, hist_type1d_pertape, &
         hist_nhtfrq,  hist_ndens, hist_mfilt, &
         hist_fincl1,  hist_fincl2, hist_fincl3, &
         hist_fincl4,  hist_fincl5, hist_fincl6, &
         hist_fexcl1,  hist_fexcl2, hist_fexcl3, &
         hist_fexcl4,  hist_fexcl5, hist_fexcl6
    namelist /elm_inparm/ hist_wrtch4diag

    ! BGC info
    namelist /elm_inparm/  &
         nu_com
    namelist /elm_inparm/  &
         nu_com_phosphatase
    namelist /elm_inparm/  &
         nu_com_nfix
    namelist /elm_inparm/ &
         forest_fert_exp
    namelist /elm_inparm/ &
         ECA_Pconst_RGspin
    namelist /elm_inparm/ &
         NFIX_PTASE_plant

    ! For experimental manipulations
    namelist /elm_inparm/ &
         startdate_add_temperature
    namelist /elm_inparm/ &
         startdate_add_co2
    namelist /elm_inparm/ &
         add_temperature
    namelist /elm_inparm/ &
         add_co2

    namelist /elm_inparm/ &
         use_pheno_flux_limiter
         
    namelist /elm_inparm/  &
         suplnitro,suplphos
    namelist /elm_inparm/ &
         nfix_timeconst
    namelist /elm_inparm/ &
         spinup_state, override_bgc_restart_mismatch_dump
    namelist /elm_inparm/ &
         nyears_ad_carbon_only, spinup_mortality_factor

    namelist /elm_inparm / &
         co2_type

    namelist /elm_inparm / &
         perchroot, perchroot_alt

    namelist /elm_inparm / &
         anoxia, anoxia_wtsat

    namelist /elm_inparm / &
         deepmixing_depthcrit, deepmixing_mixfact, lake_melt_icealb
    ! lake_melt_icealb is of dimension numrad

    ! Glacier_mec info
    namelist /elm_inparm/ &    
         maxpatch_glcmec, glc_smb, glc_do_dynglacier, glcmec_downscale_rain_snow_convert, &
         glcmec_downscale_longwave, glc_snow_persistence_max_days, glc_grid, fglcmask 

    ! Other options

    namelist /elm_inparm/  &
         clump_pproc, wrtdia, &
         create_crop_landunit, nsegspc, co2_ppmv, override_nsrest, &
         albice, more_vertlayers, subgridflag, irrigate, tw_irr, extra_gw_irr, firrig_data, all_active
    ! Urban options

    namelist /elm_inparm/  &
         urban_hac, urban_traffic

    ! vertical soil mixing variables
    namelist /elm_inparm/  &
         som_adv_flux, max_depth_cryoturb

    ! C and N input vertical profiles
    namelist /elm_inparm/  & 
          exponential_rooting_profile, rootprof_exp, surfprof_exp, pftspecific_rootingprofile

    namelist /elm_inparm / no_frozen_nitrif_denitrif

    namelist /elm_inparm / use_c13, use_c14

    namelist /elm_inparm/ fates_paramfile, use_fates,      &
          fates_spitfire_mode, use_fates_logging,        &
          use_fates_planthydro, use_fates_ed_st3,       &
          use_fates_cohort_age_tracking,                &
          use_fates_ed_prescribed_phys,                 &
          use_fates_inventory_init,                     &
          fates_inventory_ctrl_filename,                &
          use_fates_fixed_biogeog, &
          fates_parteh_mode

    namelist /elm_inparm / use_betr
        
    namelist /elm_inparm / use_lai_streams

    namelist /elm_inparm/  &
         use_c14_bombspike, atm_c14_filename

    ! All old cpp-ifdefs are below and have been converted to namelist variables 

    ! max number of plant functional types in naturally vegetated landunit
    namelist /elm_inparm/ maxpatch_pft

    namelist /elm_inparm/ &
         use_nofire, use_lch4, use_nitrif_denitrif, use_vertsoilc, use_extralakelayers, &
         use_vichydro, use_century_decomp, use_cn, use_crop, use_snicar_frc, &
         use_snicar_ad, use_vancouver, use_mexicocity, use_noio

    ! cpl_bypass variables
    namelist /elm_inparm/ metdata_type, metdata_bypass, metdata_biases, &
         co2_file, aero_file,const_climate_hist

    ! bgc & pflotran interface
    namelist /elm_inparm/ use_elm_interface, use_elm_bgc, use_pflotran

    namelist /elm_inparm/ use_dynroot

    namelist /elm_inparm/ use_var_soil_thick

    namelist /elm_inparm / &
         use_vsfm, vsfm_satfunc_type, vsfm_use_dynamic_linesearch, &
         vsfm_lateral_model_type, vsfm_include_seepage_bc

    namelist /elm_inparm/ use_hydrstress

    namelist /elm_inparm/ &
       lateral_connectivity, domain_decomp_type

    namelist /elm_inparm/ &
         use_petsc_thermal_model

    namelist /elm_inparm/ &
         do_budgets, budget_inst, budget_daily, budget_month, &
         budget_ann, budget_ltann, budget_ltend

    namelist /elm_inparm/ &
         use_erosion, ero_ccycle

    ! ----------------------------------------------------------------------
    ! Default values
    ! ----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to initialize run control settings .....'
    endif

    runtyp(:)               = 'missing'
    runtyp(nsrStartup  + 1) = 'initial'
    runtyp(nsrContinue + 1) = 'restart'
    runtyp(nsrBranch   + 1) = 'branch '

    ! Set clumps per procoessor

#if (defined _OPENMP)
    clump_pproc = omp_get_max_threads()
#else
    clump_pproc = 1
#endif

    override_nsrest = nsrest

    if (masterproc) then

       ! ----------------------------------------------------------------------
       ! Read namelist from standard input. 
       ! ----------------------------------------------------------------------

       if ( len_trim(NLFilename) == 0  )then
          call endrun(msg=' error: nlfilename not set'//errMsg(__FILE__, __LINE__))
       end if
       unitn = getavu()
       write(iulog,*) 'Read in elm_inparm namelist from: ', trim(NLFilename)
       open( unitn, file=trim(NLFilename), status='old' )
       print*,trim(NLFilename),"X.YANG debug"
       call shr_nl_find_group_name(unitn, 'elm_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, elm_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg='ERROR reading elm_inparm namelist'//errMsg(__FILE__, __LINE__))
          end if
       end if
       
       print*,"X.YANG debug SUPL NITROGEN and PHOSPHORUS ",suplnitro,suplphos
       call relavu( unitn )

       ! ----------------------------------------------------------------------
       ! Consistency checks on input namelist.
       ! ----------------------------------------------------------------------

       call set_timemgr_init( dtime_in=dtime )

       if (urban_traffic) then
          write(iulog,*)'Urban traffic fluxes are not implemented currently'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       ! History and restart files

       do i = 1, max_tapes
          if (hist_nhtfrq(i) < 0) then
             hist_nhtfrq(i) = nint(-hist_nhtfrq(i)*SHR_CONST_CDAY/(24._r8*dtime))
          endif
       end do

       ! Override start-type (can only override to branch (3)  and only 
       ! if the driver is a startup type
       if ( override_nsrest /= nsrest )then
           if ( override_nsrest /= nsrBranch .and. nsrest /= nsrStartup )then
              call endrun(msg= ' ERROR: can ONLY override elm start-type ' // &
                   'to branch type and ONLY if driver is a startup type'// &
                   errMsg(__FILE__, __LINE__))
           end if
           call elm_varctl_set( nsrest_in=override_nsrest )
       end if
       
       if (maxpatch_glcmec > 0) then
          create_glacier_mec_landunit = .true.
       else
          create_glacier_mec_landunit = .false.
       end if

       ! Check compatibility with the FATES model 
       if ( use_fates ) then

          use_voc = .false.

          if ( use_cn) then
             call endrun(msg=' ERROR: use_cn and use_fates cannot both be set to true.'//&
                   errMsg(__FILE__, __LINE__))
          end if
          
          if ( use_crop ) then
             call endrun(msg=' ERROR: use_crop and use_fates cannot both be set to true.'//&
                   errMsg(__FILE__, __LINE__))
          end if
          
          if( use_lch4 ) then
             call endrun(msg=' ERROR: use_lch4 (methane) and use_fates cannot both be set to true.'//&
                   errMsg(__FILE__, __LINE__))
          end if

          if ( n_drydep > 0 .and. drydep_method /= DD_XLND ) then
             call endrun(msg=' ERROR: dry deposition via ML Welsey is not compatible with FATES.'//&
                   errMsg(__FILE__, __LINE__))
          end if

       end if


       if (use_crop .and. (use_c13 .or. use_c14)) then
          call endrun(msg=' ERROR:: CROP and C13/C14 can NOT be on at the same time'//&
            errMsg(__FILE__, __LINE__))
       end if
       
       if (use_crop .and. .not. create_crop_landunit) then
          call endrun(msg=' ERROR: prognostic crop Patches require create_crop_landunit=.true.'//&
            errMsg(__FILE__, __LINE__))
       end if
       
       if (.not. use_erosion .and. ero_ccycle) then
          call endrun(msg=' ERROR: ero_ccycle = .true. requires erosion model active.'//&
            errMsg(__FILE__, __LINE__))
       end if
       
       if (use_lch4 .and. use_vertsoilc) then 
          anoxia = .true.
       else
          anoxia = .false.
       end if

       ! ----------------------------------------------------------------------
       !TURN OFF MEGAN VOCs if crop prognostic is on
       ! This is a temporary place holder and should be removed once MEGAN VOCs and
       ! crop ar compatible
       if (use_crop) then
          use_voc = .false.
       end if

       ! ----------------------------------------------------------------------
       ! bgc & pflotran interface
       if(.not.use_elm_interface) then
            use_elm_bgc     = .false.
            use_pflotran    = .false.
       else
       ! use_elm_interface
            if (use_elm_bgc) then
                use_pflotran = .false.
            end if

            if (use_pflotran) then
                use_elm_bgc = .false.
                ! enable 'use_nitrif_denitrif' to initilize Nh4 & NO3 pools,
                ! but NOT to implement 'nitrif_denitrif'
                use_nitrif_denitrif = .true.
            end if
       end if

    endif   ! end of if-masterproc if-block

    ! ----------------------------------------------------------------------
    ! Read in other namelists for other modules
    ! ----------------------------------------------------------------------
    !I call init_hydrology to set up default hydrology sub-module methods.
    !For future version, I suggest to  put the following two calls inside their
    !own modules, which are called from their own initializing methods
    call init_hydrology( NLFilename )
    
    call CanopyHydrology_readnl( NLFilename )

    ! ----------------------------------------------------------------------
    ! Broadcast all control information if appropriate
    ! ----------------------------------------------------------------------

    call control_spmd()
    
    if (use_pflotran) then
       call elm_pf_readnl(NLFilename)
    end if

    if (use_betr) then
       call betr_readNL( NLFilename, use_c13, use_c14)
    endif    

    ! ----------------------------------------------------------------------
    ! consistency checks
    ! ----------------------------------------------------------------------

    ! Consistency settings for co2 type
    if (co2_type /= 'constant' .and. co2_type /= 'prognostic' .and. co2_type /= 'diagnostic') then
       write(iulog,*)'co2_type = ',co2_type,' is not supported'
       call endrun(msg=' ERROR:: choices are constant, prognostic or diagnostic'//&
            errMsg(__FILE__, __LINE__))
    end if

    ! Check on run type
    if (nsrest == iundef) then
       call endrun(msg=' ERROR:: must set nsrest'//& 
            errMsg(__FILE__, __LINE__))
    end if
    if (nsrest == nsrBranch .and. nrevsn == ' ') then
       call endrun(msg=' ERROR: need to set restart data file name'//&
            errMsg(__FILE__, __LINE__))
    end if

    ! Consistency settings for co2_ppvm
    if ( (co2_ppmv <= 0.0_r8) .or. (co2_ppmv > 3000.0_r8) ) then
       call endrun(msg=' ERROR: co2_ppmv is out of a reasonable range'//& 
            errMsg(__FILE__, __LINE__))
    end if

    ! Consistency settings for nrevsn

    if (nsrest == nsrStartup ) nrevsn = ' '
    if (nsrest == nsrContinue) nrevsn = 'set by restart pointer file file'
    if (nsrest /= nsrStartup .and. nsrest /= nsrContinue .and. nsrest /= nsrBranch ) then
       call endrun(msg=' ERROR: nsrest NOT set to a valid value'//&
            errMsg(__FILE__, __LINE__))
    end if

    ! Single Column
    if ( single_column .and. (scmlat == rundef  .or. scmlon == rundef ) ) then
       call endrun(msg=' ERROR:: single column mode on -- but scmlat and scmlon are NOT set'//&
            errMsg(__FILE__, __LINE__))
       if (.not. use_lch4 .and. anoxia) then
          call endrun(msg='ERROR:: anoxia is turned on, but this currently requires turning on the CH4 submodel'//&
            errMsg(__FILE__, __LINE__))
       end if
    end if

    ! Consistency settings for vsfm settings
    if (vsfm_satfunc_type /= 'brooks_corey'             .and. &
        vsfm_satfunc_type /= 'smooth_brooks_corey_bz2'  .and. &
        vsfm_satfunc_type /= 'smooth_brooks_corey_bz3'  .and. &
        vsfm_satfunc_type /= 'van_genuchten') then
       write(iulog,*)'vsfm_satfunc_type = ',vsfm_satfunc_type,' is not supported'
       call endrun(msg=' ERROR:: choices are brooks_corey, smooth_brooks_corey_bz2, '//&
            'smooth_brooks_corey_bz3 or van_genuchten'//&
            errMsg(__FILE__, __LINE__))
    end if

    if (vsfm_lateral_model_type /= 'none'        .and. &
        vsfm_lateral_model_type /= 'source_sink' .and. &
        vsfm_lateral_model_type /= 'three_dimensional' ) then
       write(iulog,*)'vsfm_lateral_model_type = ',trim(vsfm_lateral_model_type), ' is not supported'
       call endrun(msg=' ERROR:: choices are source_sink or three_dimensional ' // &
            errMsg(__FILE__, __LINE__))
    endif

    ! Lateral connectivity
    if (.not.lateral_connectivity) then

       if (vsfm_lateral_model_type /= 'none') then
          call endrun(msg=' ERROR:: Lateral flow in VSFM requires lateral_connectivity to be true '// &
               errMsg(__FILE__, __LINE__))
       endif

       if (trim(domain_decomp_type) == 'graph_partitioning') then
          call endrun(msg=' ERROR: domain_decomp_type = graph_partitioning requires ' // &
               'lateral_connectivity to be true.'                                     // &
               errMsg(__FILE__, __LINE__))
       endif
    endif

    if (masterproc) then
       write(iulog,*) 'Successfully initialized run control settings'
       write(iulog,*)
    endif

  end subroutine control_init

  !------------------------------------------------------------------------
  subroutine control_spmd()
    !
    ! !DESCRIPTION:
    ! Distribute namelist data all processors. All program i/o is 
    ! funnelled through the master processor. Processor 0 either 
    ! reads restart/history data from the disk and distributes 
    ! it to all processors, or collects data from
    ! all processors and writes it to disk.
    !
    ! !USES:
    !
    use spmdMod,    only : mpicom, MPI_CHARACTER, MPI_INTEGER, MPI_LOGICAL, MPI_REAL8
    use elm_varpar, only : numrad
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !LOCAL VARIABLES:
    integer ier       !error code
    !-----------------------------------------------------------------------

    ! run control variables
    call mpi_bcast (caseid, len(caseid), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (ctitle, len(ctitle), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (version, len(version), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hostname, len(hostname), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (username, len(username), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (nsrest, 1, MPI_INTEGER, 0, mpicom, ier)

    call mpi_bcast (use_nofire, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_lch4, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_nitrif_denitrif, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_vertsoilc, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_extralakelayers, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_vichydro, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_century_decomp, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_cn, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_crop, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_voc, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_snicar_frc, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_snicar_ad, 1, MPI_LOGICAL, 0, mpicom, ier)   
    call mpi_bcast (use_vancouver, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_mexicocity, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_noio, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! initial file variables
    call mpi_bcast (nrevsn, len(nrevsn), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (finidat, len(finidat), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (finidat_interp_source, len(finidat_interp_source), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (finidat_interp_dest, len(finidat_interp_dest), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsurdat, len(fsurdat), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fatmlndfrc,len(fatmlndfrc),MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fatmtopo, len(fatmtopo) ,MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (flndtopo, len(flndtopo) ,MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (paramfile, len(paramfile) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsoilordercon, len(fsoilordercon) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsnowoptics, len(fsnowoptics),  MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsnowaging,  len(fsnowaging),   MPI_CHARACTER, 0, mpicom, ier)

    ! Irrigation
    call mpi_bcast(irrigate, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast(tw_irr, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast(extra_gw_irr, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast(firrig_data, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! Landunit generation
    call mpi_bcast(create_crop_landunit, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! Other subgrid logic
    call mpi_bcast(all_active, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! max number of plant functional types in naturally vegetated landunit
    call mpi_bcast(maxpatch_pft, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! BGC
    call mpi_bcast (co2_type, len(co2_type), MPI_CHARACTER, 0, mpicom, ier)
    if (use_cn) then
       call mpi_bcast (suplnitro, len(suplnitro), MPI_CHARACTER, 0, mpicom, ier)
       call mpi_bcast (nfix_timeconst, 1, MPI_REAL8, 0, mpicom, ier)
       call mpi_bcast (spinup_state, 1, MPI_INTEGER, 0, mpicom, ier)
       call mpi_bcast (nyears_ad_carbon_only, 1, MPI_INTEGER, 0, mpicom, ier)
       call mpi_bcast (spinup_mortality_factor, 1, MPI_REAL8, 0, mpicom, ier)
       call mpi_bcast (override_bgc_restart_mismatch_dump, 1, MPI_LOGICAL, 0, mpicom, ier)
    end if
    call mpi_bcast (suplphos, len(suplphos), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (nu_com, len(nu_com), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (nu_com_phosphatase, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (nu_com_nfix, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (forest_fert_exp, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (ECA_Pconst_RGspin, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (NFIX_PTASE_plant, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_pheno_flux_limiter, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (startdate_add_temperature, 1, MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (startdate_add_co2, 1, MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (add_co2, 1, MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (add_temperature, 1, MPI_REAL8, 0, mpicom, ier)

    ! isotopes
    call mpi_bcast (use_c13, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_c14, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (use_fates, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (fates_spitfire_mode, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (fates_paramfile, len(fates_paramfile) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (use_fates_logging, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_planthydro, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_cohort_age_tracking, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_ed_st3, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_fixed_biogeog, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_ed_prescribed_phys,  1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_fates_inventory_init, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (fates_inventory_ctrl_filename, len(fates_inventory_ctrl_filename), &
          MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fates_parteh_mode, 1, MPI_INTEGER, 0, mpicom, ier)


    call mpi_bcast (use_betr, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (use_lai_streams, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (use_dynroot, 1, MPI_LOGICAL, 0, mpicom, ier)

    if (use_cn .and. use_vertsoilc) then
       ! vertical soil mixing variables
       call mpi_bcast (som_adv_flux, 1, MPI_REAL8,  0, mpicom, ier)
       call mpi_bcast (max_depth_cryoturb, 1, MPI_REAL8,  0, mpicom, ier)

       ! C and N input vertical profiles
       call mpi_bcast (exponential_rooting_profile,       1, MPI_LOGICAL,  0, mpicom, ier)
       call mpi_bcast (rootprof_exp,            1, MPI_REAL8,  0, mpicom, ier)
       call mpi_bcast (surfprof_exp,            1, MPI_REAL8,  0, mpicom, ier)
       call mpi_bcast (pftspecific_rootingprofile,        1, MPI_LOGICAL,  0, mpicom, ier)
    end if

    if (use_cn .and. use_nitrif_denitrif) then 
       call mpi_bcast (no_frozen_nitrif_denitrif,  1, MPI_LOGICAL, 0, mpicom, ier)
    end if

    if (use_cn) then
       call mpi_bcast (use_c14_bombspike,  1, MPI_LOGICAL, 0, mpicom, ier)
       call mpi_bcast (atm_c14_filename,  len(atm_c14_filename), MPI_CHARACTER, 0, mpicom, ier)
    end if

    call mpi_bcast (perchroot, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (perchroot_alt, 1, MPI_LOGICAL, 0, mpicom, ier)
    if (use_lch4) then
       call mpi_bcast (anoxia, 1, MPI_LOGICAL, 0, mpicom, ier)
       call mpi_bcast (anoxia_wtsat, 1, MPI_LOGICAL, 0, mpicom, ier)
    end if

    ! lakes
    call mpi_bcast (deepmixing_depthcrit,  1, MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (deepmixing_mixfact,    1, MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (lake_melt_icealb, numrad, MPI_REAL8, 0, mpicom, ier)

    ! physics variables
    call mpi_bcast (urban_hac, len(urban_hac), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (urban_traffic , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (nsegspc, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (subgridflag , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (wrtdia, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (single_column,1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (scmlat, 1, MPI_REAL8,0, mpicom, ier)
    call mpi_bcast (scmlon, 1, MPI_REAL8,0, mpicom, ier)
    call mpi_bcast (co2_ppmv, 1, MPI_REAL8,0, mpicom, ier)
    call mpi_bcast (albice, 2, MPI_REAL8,0, mpicom, ier)
    call mpi_bcast (more_vertlayers,1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (const_climate_hist, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! glacier_mec variables
    call mpi_bcast (create_glacier_mec_landunit, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (maxpatch_glcmec, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (glc_smb, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (glc_do_dynglacier, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (glcmec_downscale_rain_snow_convert, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (glcmec_downscale_longwave, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (glc_snow_persistence_max_days, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (glc_grid, len(glc_grid), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fglcmask, len(fglcmask), MPI_CHARACTER, 0, mpicom, ier)

    ! history file variables
    call mpi_bcast (hist_empty_htapes, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_dov2xy, size(hist_dov2xy), MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_nhtfrq, size(hist_nhtfrq), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_mfilt, size(hist_mfilt), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_ndens, size(hist_ndens), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_avgflag_pertape, size(hist_avgflag_pertape), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_type1d_pertape, max_namlen*size(hist_type1d_pertape), MPI_CHARACTER, 0, mpicom, ier)
    if (use_lch4) then
       call mpi_bcast (hist_wrtch4diag, 1, MPI_LOGICAL, 0, mpicom, ier)
    end if
    call mpi_bcast (hist_fexcl1, max_namlen*size(hist_fexcl1), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl2, max_namlen*size(hist_fexcl2), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl3, max_namlen*size(hist_fexcl3), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl4, max_namlen*size(hist_fexcl4), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl5, max_namlen*size(hist_fexcl5), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl6, max_namlen*size(hist_fexcl6), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl1, (max_namlen+2)*size(hist_fincl1), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl2, (max_namlen+2)*size(hist_fincl2), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl3, (max_namlen+2)*size(hist_fincl3), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl4, (max_namlen+2)*size(hist_fincl4), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl5, (max_namlen+2)*size(hist_fincl5), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl6, (max_namlen+2)*size(hist_fincl6), MPI_CHARACTER, 0, mpicom, ier)

    ! restart file variables

    call mpi_bcast (rpntfil, len(rpntfil), MPI_CHARACTER, 0, mpicom, ier)

    ! clump decomposition variables

    call mpi_bcast (clump_pproc, 1, MPI_INTEGER, 0, mpicom, ier)

    ! lateral connectivity
    call mpi_bcast (lateral_connectivity, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (domain_decomp_type, len(domain_decomp_type), MPI_CHARACTER, 0, mpicom, ier)

    ! bgc & pflotran interface
    call mpi_bcast (use_elm_interface, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_elm_bgc, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_pflotran, 1, MPI_LOGICAL, 0, mpicom, ier)
    
    !cpl_bypass
     call mpi_bcast (metdata_type,   len(metdata_type),   MPI_CHARACTER, 0, mpicom, ier)
     call mpi_bcast (metdata_bypass, len(metdata_bypass), MPI_CHARACTER, 0, mpicom, ier)
     call mpi_bcast (metdata_biases, len(metdata_biases), MPI_CHARACTER, 0, mpicom, ier)
     call mpi_bcast (co2_file,       len(co2_file),       MPI_CHARACTER, 0, mpicom, ier)
     call mpi_bcast (aero_file,      len(aero_file),      MPI_CHARACTER, 0, mpicom, ier)

    ! plant hydraulics
    call mpi_bcast (use_hydrstress, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! VSFM variable

    call mpi_bcast (use_vsfm                   , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (vsfm_use_dynamic_linesearch, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (vsfm_include_seepage_bc    , 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (vsfm_satfunc_type      , len(vsfm_satfunc_type)      , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (vsfm_lateral_model_type, len(vsfm_lateral_model_type), MPI_CHARACTER, 0, mpicom, ier)

    ! PETSc-based thermal model
    call mpi_bcast (use_petsc_thermal_model, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! soil erosion
    call mpi_bcast (use_erosion, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (ero_ccycle , 1, MPI_LOGICAL, 0, mpicom, ier)

    ! Budget
    call mpi_bcast (do_budgets   , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (budget_inst  , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (budget_daily , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (budget_month , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (budget_ann   , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (budget_ltann , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (budget_ltend , 1, MPI_INTEGER, 0, mpicom, ier)

  end subroutine control_spmd

  !------------------------------------------------------------------------
  subroutine control_print ()
    !
    ! !DESCRIPTION:
    ! Write out the clm namelist run control variables
    !
    ! !USES:
    !
    use AllocationMod, only : suplnitro, suplnNon
    use AllocationMod, only : suplphos, suplpNon
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !LOCAL VARIABLES:
    integer i  !loop index
    character(len=32) :: subname = 'control_print'  ! subroutine name
    !------------------------------------------------------------------------

    write(iulog,*) 'define run:'
    write(iulog,*) '   source                = ',trim(source)
    write(iulog,*) '   model_version         = ',trim(version)
    write(iulog,*) '   run type              = ',runtyp(nsrest+1)
    write(iulog,*) '   case title            = ',trim(ctitle)
    write(iulog,*) '   username              = ',trim(username)
    write(iulog,*) '   hostname              = ',trim(hostname)
    write(iulog,*) 'process control parameters:'
    write(iulog,*) '    use_nofire = ', use_nofire
    write(iulog,*) '    use_lch4 = ', use_lch4
    write(iulog,*) '    use_nitrif_denitrif = ', use_nitrif_denitrif
    write(iulog,*) '    use_vertsoilc = ', use_vertsoilc
    write(iulog,*) '    use_var_soil_thick = ', use_var_soil_thick
    write(iulog,*) '    use_extralakelayers = ', use_extralakelayers
    write(iulog,*) '    use_vichydro = ', use_vichydro
    write(iulog,*) '    use_century_decomp = ', use_century_decomp
    write(iulog,*) '    use_cn = ', use_cn
    write(iulog,*) '    use_crop = ', use_crop
    write(iulog,*) '    irrigate = ', irrigate
    write(iulog,*) '    two-way irrigation = ', tw_irr
    write(iulog,*) '    use_snicar_frc = ', use_snicar_frc
    write(iulog,*) '    use_snicar_ad = ', use_snicar_ad
    write(iulog,*) '    use_vancouver = ', use_vancouver
    write(iulog,*) '    use_mexicocity = ', use_mexicocity
    write(iulog,*) '    use_noio = ', use_noio
    write(iulog,*) '    use_betr = ', use_betr
    write(iulog,*) 'input data files:'
    write(iulog,*) '   PFT physiology and parameters file = ',trim(paramfile)
    write(iulog,*) '   Soil order dependent parameters file = ',trim(fsoilordercon)
    if (fsurdat == ' ') then
       write(iulog,*) '   fsurdat, surface dataset not set'
    else
       write(iulog,*) '   surface data   = ',trim(fsurdat)
    end if
    if (fatmlndfrc == ' ') then
       write(iulog,*) '   fatmlndfrc not set, setting frac/mask to 1'
    else
       write(iulog,*) '   land frac data = ',trim(fatmlndfrc)
    end if
    if (flndtopo == ' ') then
       write(iulog,*) '   flndtopo not set'
    else
       write(iulog,*) '   land topographic data = ',trim(flndtopo)
    end if
    if (fatmtopo == ' ') then
       write(iulog,*) '   fatmtopo not set'
    else
       write(iulog,*) '   atm topographic data = ',trim(fatmtopo)
    end if
    if (use_cn) then
       if (suplnitro /= suplnNon)then
          write(iulog,*) '   Supplemental Nitrogen mode is set to run over Patches: ', &
               trim(suplnitro)
       end if
       
       if (nfix_timeconst /= 0._r8) then
          write(iulog,*) '   nfix_timeconst, timescale for smoothing npp in N fixation term: ', nfix_timeconst
       else
          write(iulog,*) '   nfix_timeconst == zero, use standard N fixation scheme. '
       end if
       
       write(iulog,*) '   spinup_state, (0 = normal mode; 1 = AD spinup)         : ', spinup_state
       if ( spinup_state .eq. 0 ) then
          write(iulog,*) '   model is currently NOT in AD spinup mode.'
       else if ( spinup_state .eq. 1 ) then
          write(iulog,*) '   model is currently in AD spinup mode.'
          write(iulog,*) '   nyears in carbon only mode = ', nyears_ad_carbon_only
          write(iulog,*) '   dead wood mortality acceleration = ',spinup_mortality_factor
       else
          call endrun(msg=' error: spinup_state can only have integer value of 0 or 1'//&
               errMsg(__FILE__, __LINE__))
       end if
       
       write(iulog,*) '   override_bgc_restart_mismatch_dump                     : ', override_bgc_restart_mismatch_dump
    end if
       if (suplphos /= suplpNon)then
          write(iulog,*) '   Supplemental Phosphorus mode is set to run over Patches: ', &
               trim(suplphos)
       end if

    if (use_cn .and. use_vertsoilc) then
       write(iulog, *) '   som_adv_flux, the advection term in soil mixing (m/s) : ', som_adv_flux
       write(iulog, *) '   max_depth_cryoturb (m)                                : ', max_depth_cryoturb
       
       write(iulog, *) '   exponential_rooting_profile                           : ', exponential_rooting_profile
       write(iulog, *) '   rootprof_exp                                          : ', rootprof_exp
       write(iulog, *) '   surfprof_exp                                          : ', surfprof_exp
       write(iulog, *) '   pftspecific_rootingprofile                            : ', pftspecific_rootingprofile
       write(iulog, *) '   dynamic roots                                         : ', use_dynroot
    end if
       
    if (use_cn .and. .not. use_nitrif_denitrif) then
       write(iulog, *) '   no_frozen_nitrif_denitrif                             : ', no_frozen_nitrif_denitrif
    end if

    if (use_cn) then
       write(iulog, *) '  use_c13                                                : ', use_c13
       write(iulog, *) '  use_c14                                                : ', use_c14
       write(iulog, *) '  use_c14_bombspike                                      : ', use_c14_bombspike
       write(iulog, *) '  atm_c14_filename                                       : ', atm_c14_filename
    end if

    if (fsnowoptics == ' ') then
       write(iulog,*) '   snow optical properties file NOT set'
    else
       write(iulog,*) '   snow optical properties file = ',trim(fsnowoptics)
    endif
    if (fsnowaging == ' ') then
       write(iulog,*) '   snow aging parameters file NOT set'
    else
       write(iulog,*) '   snow aging parameters file = ',trim(fsnowaging)
    endif

    if (create_glacier_mec_landunit) then
       write(iulog,*) '   glc number of elevation classes =', maxpatch_glcmec
       write(iulog,*) '   glc grid for glacier mask file = ',trim(glc_grid)
       write(iulog,*) '   glc glacier mask file = ',trim(fglcmask)
       
       write(iulog,*) '   Max snow depth (mm) =', h2osno_max
       if (glcmec_downscale_rain_snow_convert) then
          write(iulog,*) '   Rain and snow will be converted based on surface temperature'
       else
          write(iulog,*) '   Rain and snow will NOT be converted based on surface temperature'
       endif
       if (glcmec_downscale_longwave) then
          write(iulog,*) '   Longwave radiation will be downscaled'
       else
          write(iulog,*) '   Longwave radiation will NOT be downscaled'
       endif
       if (glc_do_dynglacier) then
          write(iulog,*) '   glc ELM glacier areas and topography WILL evolve dynamically'
       else
          write(iulog,*) '   glc ELM glacier areas and topography will NOT evolve dynamically'
       end if
       if (glc_smb) then
          write(iulog,*) '   glc surface mass balance will be passed to ice sheet model'
       else
          write(iulog,*) '   glc positive-degree-day info will be passed to ice sheet model'
       endif
       write(iulog,*) '   glc snow persistence max days = ', glc_snow_persistence_max_days
    endif

    if (nsrest == nsrStartup .and. finidat == ' ') write(iulog,*) '   initial data created by model'
    if (nsrest == nsrStartup .and. finidat /= ' ') write(iulog,*) '   initial data   = ',trim(finidat)
    if (nsrest /= nsrStartup) write(iulog,*) '   restart data   = ',trim(nrevsn)
    write(iulog,*) '   atmospheric forcing data is from cesm atm model'
    write(iulog,*) 'Restart parameters:'
    write(iulog,*)'   restart pointer file directory     = ',trim(rpntdir)
    write(iulog,*)'   restart pointer file name          = ',trim(rpntfil)
    write(iulog,*) 'model physics parameters:'

    if ( trim(co2_type) == 'constant' )then
       write(iulog,*) '   CO2 volume mixing ratio   (umol/mol)   = ', co2_ppmv
    else
       write(iulog,*) '   CO2 volume mixing ratio                = ', co2_type
    end if

    write(iulog,*) '   constant historical climate during transient simulation = ', const_climate_hist

    write(iulog,*) '   land-ice albedos      (unitless 0-1)   = ', albice
    write(iulog,*) '   urban air conditioning/heating and wasteheat   = ', urban_hac
    write(iulog,*) '   urban traffic flux   = ', urban_traffic
    write(iulog,*) '   more vertical layers = ', more_vertlayers
    if (nsrest == nsrContinue) then
       write(iulog,*) 'restart warning:'
       write(iulog,*) '   Namelist not checked for agreement with initial run.'
       write(iulog,*) '   Namelist should not differ except for ending time step and run type'
    end if
    if (nsrest == nsrBranch) then
       write(iulog,*) 'branch warning:'
       write(iulog,*) '   Namelist not checked for agreement with initial run.'
       write(iulog,*) '   Surface data set and reference date should not differ from initial run'
    end if
    write(iulog,*) '   maxpatch_pft         = ',maxpatch_pft
    write(iulog,*) '   nsegspc              = ',nsegspc
    ! New fields
    write(iulog,*) ' perchroot (plant water stress based on unfrozen layers only) = ',perchroot
    write(iulog,*) ' perchroot (plant water stress based on time-integrated active layer only) = ',perchroot
    if (use_lch4) then
       write(iulog,*) ' anoxia (applied to soil decomposition)             = ',anoxia
       write(iulog,*) ' anoxia_wtsat (weight anoxia by inundated fraction) = ',anoxia_wtsat
    end if
    ! Lakes
    write(iulog,*)
    write(iulog,*) 'Lake Model Namelists:'
    write(iulog,*) 'Increased mixing relative to Hostetler wind-driven eddy expression ',&
                   'will be used for deep lakes exceeding depth ', deepmixing_depthcrit,&
                      ' by a factor of ', deepmixing_mixfact, '.'
    write(iulog,*) 'Albedo over melting lakes will approach values (visible, NIR):', lake_melt_icealb, &
                   'as compared with 0.60, 0.40 for cold frozen lakes with no snow.'

    ! FATES
    write(iulog, *) '    use_fates = ', use_fates
    if (use_fates) then
       write(iulog, *) '    fates_spitfire_mode = ', fates_spitfire_mode
       write(iulog, *) '    use_fates_logging = ', use_fates_logging
       write(iulog, *) '    fates_paramfile = ', fates_paramfile
       write(iulog, *) '    use_fates_planthydro = ', use_fates_planthydro
       write(iulog, *) '    use_fates_cohort_age_tracking = ',use_fates_cohort_age_tracking
       write(iulog, *) '    fates_parteh_mode = ', fates_parteh_mode
       write(iulog, *) '    use_fates_ed_st3 = ',use_fates_ed_st3
       write(iulog, *) '    use_fates_ed_prescribed_phys = ',use_fates_ed_prescribed_phys
       write(iulog, *) '    use_fates_inventory_init = ',use_fates_inventory_init
       write(iulog, *) '    use_fates_fixed_biogeog = ', use_fates_fixed_biogeog
       write(iulog, *) '    fates_inventory_ctrl_filename = ',fates_inventory_ctrl_filename
    end if

    ! VSFM
    if (use_vsfm) then
       write(iulog,*)
       write(iulog,*) 'VSFM Namelists:'
       write(iulog, *) '  vsfm_satfunc_type                                      : ', vsfm_satfunc_type
       write(iulog, *) '  vsfm_use_dynamic_linesearch                            : ', vsfm_use_dynamic_linesearch
       write(iulog,*) '  vsfm_lateral_model_type                                 : ', vsfm_lateral_model_type
    endif

  end subroutine control_print

end module controlMod
