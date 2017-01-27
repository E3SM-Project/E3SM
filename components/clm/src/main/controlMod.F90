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
  use clm_varctl   
  use shr_kind_mod            , only: r8 => shr_kind_r8, SHR_KIND_CL
  use shr_nl_mod              , only: shr_nl_find_group_name
  use shr_const_mod           , only: SHR_CONST_CDAY
  use shr_log_mod             , only: errMsg => shr_log_errMsg
  use abortutils              , only: endrun
  use spmdMod                 , only: masterproc
  use decompMod               , only: clump_pproc
  use clm_varpar              , only: maxpatch_pft, maxpatch_glcmec, more_vertlayers
  use histFileMod             , only: max_tapes, max_namlen 
  use histFileMod             , only: hist_empty_htapes, hist_dov2xy, hist_avgflag_pertape, hist_type1d_pertape 
  use histFileMod             , only: hist_nhtfrq, hist_ndens, hist_mfilt, hist_fincl1, hist_fincl2, hist_fincl3
  use histFileMod             , only: hist_fincl4, hist_fincl5, hist_fincl6, hist_fexcl1, hist_fexcl2, hist_fexcl3
  use histFileMod             , only: hist_fexcl4, hist_fexcl5, hist_fexcl6
  use LakeCon                 , only: deepmixing_depthcrit, deepmixing_mixfact 
  use CNAllocationMod         , only: suplnitro
  use CNAllocationMod         , only: suplphos
  use CNCarbonFluxType        , only: nfix_timeconst
  use CNNitrifDenitrifMod     , only: no_frozen_nitrif_denitrif
  use CNC14DecayMod           , only: use_c14_bombspike, atm_c14_filename
  use CNSoilLittVertTranspMod , only: som_adv_flux, max_depth_cryoturb
  use CNVerticalProfileMod    , only: exponential_rooting_profile, rootprof_exp 
  use CNVerticalProfileMod    , only: surfprof_exp, pftspecific_rootingprofile  
  use CNSharedParamsMod       , only: anoxia_wtsat
  use CanopyfluxesMod         , only: perchroot, perchroot_alt
  use CanopyHydrologyMod      , only: CanopyHydrology_readnl
  use SurfaceAlbedoMod        , only: albice, lake_melt_icealb
  use UrbanParamsType         , only: urban_hac, urban_traffic
  use clm_varcon              , only: h2osno_max
  use clm_varctl              , only: use_dynroot
  use CNAllocationMod         , only: nu_com_phosphatase,nu_com_nfix 
  use clm_varctl              , only: nu_com, do_varsoil
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
    use clm_pflotran_interfaceMod , only : clm_pf_readnl
    use betr_initializeMod        , only : betr_readNL
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
    namelist / clm_inparm/ &
    dtime

    ! CLM namelist settings

    namelist /clm_inparm / &
         fatmlndfrc, finidat, nrevsn, &
         finidat_interp_source, finidat_interp_dest

    ! Input datasets

    namelist /clm_inparm/  &
         fsurdat, fatmtopo, flndtopo, &
         paramfile, flanduse_timeseries,  fsnowoptics, fsnowaging,fsoilordercon


    ! History, restart options

    namelist /clm_inparm/  &
         hist_empty_htapes, hist_dov2xy, &
         hist_avgflag_pertape, hist_type1d_pertape, &
         hist_nhtfrq,  hist_ndens, hist_mfilt, &
         hist_fincl1,  hist_fincl2, hist_fincl3, &
         hist_fincl4,  hist_fincl5, hist_fincl6, &
         hist_fexcl1,  hist_fexcl2, hist_fexcl3, &
         hist_fexcl4,  hist_fexcl5, hist_fexcl6
    namelist /clm_inparm/ hist_wrtch4diag

    ! BGC info
    namelist /clm_inparm/  &
         nu_com
    namelist /clm_inparm/  &
         nu_com_phosphatase
    namelist /clm_inparm/  &
         nu_com_nfix
         
    namelist /clm_inparm/  &
         suplnitro,suplphos
    namelist /clm_inparm/ &
         nfix_timeconst
    namelist /clm_inparm/ &
         spinup_state, override_bgc_restart_mismatch_dump
    namelist /clm_inparm/ &
         nyears_ad_carbon_only, spinup_mortality_factor

    namelist /clm_inparm / &
         co2_type

    namelist /clm_inparm / &
         perchroot, perchroot_alt

    namelist /clm_inparm / &
         anoxia, anoxia_wtsat

    namelist /clm_inparm / &
         deepmixing_depthcrit, deepmixing_mixfact, lake_melt_icealb
    ! lake_melt_icealb is of dimension numrad

    ! Glacier_mec info
    namelist /clm_inparm/ &    
         maxpatch_glcmec, glc_smb, glc_do_dynglacier, glcmec_downscale_rain_snow_convert, &
         glcmec_downscale_longwave, glc_snow_persistence_max_days, glc_grid, fglcmask 

    ! Other options

    namelist /clm_inparm/  &
         clump_pproc, wrtdia, &
         create_crop_landunit, nsegspc, co2_ppmv, override_nsrest, &
         albice, more_vertlayers, subgridflag, irrigate, all_active
    ! Urban options

    namelist /clm_inparm/  &
         urban_hac, urban_traffic

    ! vertical soil mixing variables
    namelist /clm_inparm/  &
         som_adv_flux, max_depth_cryoturb

    ! C and N input vertical profiles
    namelist /clm_inparm/  & 
          exponential_rooting_profile, rootprof_exp, surfprof_exp, pftspecific_rootingprofile

    namelist /clm_inparm / no_frozen_nitrif_denitrif

    namelist /clm_inparm / use_c13, use_c14

    namelist /clm_inparm / use_ed, use_ed_spit_fire
    
    namelist /clm_inparm / use_betr
        
    namelist /clm_inparm / use_lai_streams

    namelist /clm_inparm/  &
         use_c14_bombspike, atm_c14_filename

    ! All old cpp-ifdefs are below and have been converted to namelist variables 

    ! max number of plant functional types in naturally vegetated landunit
    namelist /clm_inparm/ maxpatch_pft

    namelist /clm_inparm/ &
         use_nofire, use_lch4, use_nitrif_denitrif, use_vertsoilc, use_extralakelayers, &
         use_vichydro, use_century_decomp, use_cn, use_cndv, use_crop, use_snicar_frc, &
         use_vancouver, use_mexicocity, use_noio

    ! cpl_bypass variables
    namelist /clm_inparm/ metdata_type, metdata_bypass, metdata_biases, &
         co2_file, aero_file

    ! bgc & pflotran interface
    namelist /clm_inparm/ use_bgc_interface, use_clm_bgc, use_pflotran

    namelist /clm_inparm/ use_dynroot

    namelist /clm_inparm/ do_varsoil

    namelist /clm_inparm / &
         use_vsfm, vsfm_satfunc_type, vsfm_use_dynamic_linesearch

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
       write(iulog,*) 'Read in clm_inparm namelist from: ', trim(NLFilename)
       open( unitn, file=trim(NLFilename), status='old' )
       print*,trim(NLFilename),"X.YANG debug"
       call shr_nl_find_group_name(unitn, 'clm_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg='ERROR reading clm_inparm namelist'//errMsg(__FILE__, __LINE__))
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
          if (hist_nhtfrq(i) == 0) then
             hist_mfilt(i) = 1
          else if (hist_nhtfrq(i) < 0) then
             hist_nhtfrq(i) = nint(-hist_nhtfrq(i)*SHR_CONST_CDAY/(24._r8*dtime))
          endif
       end do

       ! Override start-type (can only override to branch (3)  and only 
       ! if the driver is a startup type
       if ( override_nsrest /= nsrest )then
           if ( override_nsrest /= nsrBranch .and. nsrest /= nsrStartup )then
              call endrun(msg= ' ERROR: can ONLY override clm start-type ' // &
                   'to branch type and ONLY if driver is a startup type'// &
                   errMsg(__FILE__, __LINE__))
           end if
           call clm_varctl_set( nsrest_in=override_nsrest )
       end if
       
       if (maxpatch_glcmec > 0) then
          create_glacier_mec_landunit = .true.
       else
          create_glacier_mec_landunit = .false.
       end if

       if (use_crop .and. (use_c13 .or. use_c14)) then
          call endrun(msg=' ERROR:: CROP and C13/C14 can NOT be on at the same time'//&
            errMsg(__FILE__, __LINE__))
       end if
       
       if (use_crop .and. .not. create_crop_landunit) then
          call endrun(msg=' ERROR: prognostic crop Patches require create_crop_landunit=.true.'//&
            errMsg(__FILE__, __LINE__))
       end if
       
       if (use_crop .and. flanduse_timeseries /= ' ') then
          call endrun(msg=' ERROR: prognostic crop is incompatible with transient landuse'//&
               errMsg(__FILE__, __LINE__))
       end if
       
       if (.not. use_crop .and. irrigate) then
          call endrun(msg=' ERROR: irrigate = .true. requires CROP model active.'//&
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
       !! bgc & pflotran interface
       if(.not.use_bgc_interface) then
            use_clm_bgc     = .true.
            use_pflotran    = .false.
       end if

       if (use_clm_bgc) then
            use_pflotran = .false.
       end if

       if (use_pflotran) then
            use_clm_bgc = .false.
            !! enable 'use_nitrif_denitrif' to initilize Nh4 & NO3 pools, NOT to implement 'nitrif_denitrif'
            use_nitrif_denitrif = .true.
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
       call clm_pf_readnl(NLFilename)
    end if

    if (use_betr) then
       call betr_readNL( NLFilename )
    endif    

    ! ----------------------------------------------------------------------
    ! consistency checks
    ! ----------------------------------------------------------------------

    if (flanduse_timeseries /= ' ' .and. create_crop_landunit) then
       call endrun(msg=' ERROR:: dynamic landuse is currently not supported with create_crop_landunit option'//&
            errMsg(__FILE__, __LINE__))
    end if
    if (flanduse_timeseries /= ' ' .and. use_cndv) then
       call endrun(msg=' ERROR:: dynamic landuse is currently not supported with CNDV option'//&
            errMsg(__FILE__, __LINE__))
    end if

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

    ! Consistency settings for co2 type
    if (vsfm_satfunc_type /= 'brooks_corey'             .and. &
        vsfm_satfunc_type /= 'smooth_brooks_corey_bz2'  .and. &
        vsfm_satfunc_type /= 'smooth_brooks_corey_bz3'  .and. &
        vsfm_satfunc_type /= 'van_genuchten') then
       write(iulog,*)'vsfm_satfunc_type = ',vsfm_satfunc_type,' is not supported'
       call endrun(msg=' ERROR:: choices are brooks_corey, smooth_brooks_corey_bz2, '//&
            'smooth_brooks_corey_bz3 or van_genuchten'//&
            errMsg(__FILE__, __LINE__))
    end if

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
    use clm_varpar, only : numrad
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
    call mpi_bcast (use_cndv, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_crop, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_voc, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_snicar_frc, 1, MPI_LOGICAL, 0, mpicom, ier)
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
    call mpi_bcast (flanduse_timeseries , len(flanduse_timeseries) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsnowoptics, len(fsnowoptics),  MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsnowaging,  len(fsnowaging),   MPI_CHARACTER, 0, mpicom, ier)

    ! Irrigation
    call mpi_bcast(irrigate, 1, MPI_LOGICAL, 0, mpicom, ier)

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

    ! isotopes
    call mpi_bcast (use_c13, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_c14, 1, MPI_LOGICAL, 0, mpicom, ier)

    call mpi_bcast (use_ed, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_ed_spit_fire, 1, MPI_LOGICAL, 0, mpicom, ier)

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

    ! bgc & pflotran interface
    call mpi_bcast (use_bgc_interface, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_clm_bgc, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (use_pflotran, 1, MPI_LOGICAL, 0, mpicom, ier)
    
    !cpl_bypass
     call mpi_bcast (metdata_type,   len(metdata_type),   MPI_CHARACTER, 0, mpicom, ier)
     call mpi_bcast (metdata_bypass, len(metdata_bypass), MPI_CHARACTER, 0, mpicom, ier)
     call mpi_bcast (metdata_biases, len(metdata_biases), MPI_CHARACTER, 0, mpicom, ier)
     call mpi_bcast (co2_file,       len(co2_file),       MPI_CHARACTER, 0, mpicom, ier)
     call mpi_bcast (aero_file,      len(aero_file),      MPI_CHARACTER, 0, mpicom, ier)


    ! VSFM variable

    call mpi_bcast (use_vsfm, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (vsfm_satfunc_type, len(vsfm_satfunc_type), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (vsfm_use_dynamic_linesearch, 1, MPI_LOGICAL, 0, mpicom, ier)

  end subroutine control_spmd

  !------------------------------------------------------------------------
  subroutine control_print ()
    !
    ! !DESCRIPTION:
    ! Write out the clm namelist run control variables
    !
    ! !USES:
    !
    use CNAllocationMod, only : suplnitro, suplnNon
    use CNAllocationMod, only : suplphos, suplpNon
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
    write(iulog,*) '    do_varsoil = ', do_varsoil
    write(iulog,*) '    use_extralakelayers = ', use_extralakelayers
    write(iulog,*) '    use_vichydro = ', use_vichydro
    write(iulog,*) '    use_century_decomp = ', use_century_decomp
    write(iulog,*) '    use_cn = ', use_cn
    write(iulog,*) '    use_cndv = ', use_cndv
    write(iulog,*) '    use_crop = ', use_crop
    write(iulog,*) '    use_snicar_frc = ', use_snicar_frc
    write(iulog,*) '    use_vancouver = ', use_vancouver
    write(iulog,*) '    use_mexicocity = ', use_mexicocity
    write(iulog,*) '    use_noio = ', use_noio

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
          write(iulog,*) '   glc CLM glacier areas and topography WILL evolve dynamically'
       else
          write(iulog,*) '   glc CLM glacier areas and topography will NOT evolve dynamically'
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

    ! VSFM
    if (use_vsfm) then
       write(iulog,*)
       write(iulog,*) 'VSFM Namelists:'
       write(iulog, *) '  vsfm_satfunc_type                                      : ', vsfm_satfunc_type
       write(iulog, *) '  vsfm_use_dynamic_linesearch                            : ', vsfm_use_dynamic_linesearch
    endif

  end subroutine control_print

end module controlMod
