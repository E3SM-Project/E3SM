module controlMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: controlMod
!
! !DESCRIPTION:
! Module which initializes run control variables. The following possible
! namelist variables are set default values and possibly read in on startup
!
! Note: For definitions of namelist variables see
!       ../../bld/namelist_files/namelist_definition.xml
!       Display the file in a browser to see it neatly formatted in html.
!

! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8, SHR_KIND_CL
  use clm_varpar   , only : maxpatch_pft, maxpatch_glcmec, more_vertlayers
  use clm_varctl   , only : caseid, ctitle, nsrest, brnch_retain_casename, hostname, &
                            model_version=>version,    &
                            iulog, outnc_large_files, finidat, fsurdat, fatmlndfrc,  &
                            fatmtopo, flndtopo, fpftdyn, fpftcon, nrevsn, &
                            create_crop_landunit, allocate_all_vegpfts,   &
                            co2_type, wrtdia, co2_ppmv, nsegspc, pertlim,       &
                            username, fsnowaging, fsnowoptics, fglcmask, &
                            create_glacier_mec_landunit, glc_dyntopo, glc_smb, &
                            glc_topomax, glc_grid, subgridflag, &
                            use_c13, use_c14, irrigate, &
                            spinup_state, override_bgc_restart_mismatch_dump
  use CanopyFluxesMod , only : perchroot, perchroot_alt
#if (defined LCH4) && (defined CN)
  use clm_varctl   , only : anoxia
#ifndef CENTURY_DECOMP
  use CNDecompCascadeMod_BGC , only : anoxia_wtsat
#else
  use CNDecompCascadeMod_CENTURY , only : anoxia_wtsat
#endif
#endif
! Lakes
  use SLakeCon, only : deepmixing_depthcrit, deepmixing_mixfact, lake_melt_icealb
  ! lake_use_old_fcrit_minz0, lakepuddling, lake_puddle_thick, and lake_no_ed are currently hardwired.
!
  use SurfaceAlbedoMod, only : albice
#ifdef CN
  use CNAllocationMod, only : suplnitro
  use CNNDynamicsMod, only : nfix_timeconst
#endif
  use spmdMod      , only : masterproc
  use decompMod    , only : clump_pproc
  use histFileMod  , only : max_tapes, max_namlen, &
                            hist_empty_htapes, hist_dov2xy, &
                            hist_avgflag_pertape, hist_type1d_pertape, &
                            hist_nhtfrq, hist_ndens, hist_mfilt, &
                            hist_fincl1, hist_fincl2, hist_fincl3, &
                            hist_fincl4, hist_fincl5, hist_fincl6, &
                            hist_fexcl1, hist_fexcl2, hist_fexcl3, &
                            hist_fexcl4, hist_fexcl5, hist_fexcl6
#ifdef LCH4
  use histFileMod  , only : hist_wrtch4diag
#endif
  use shr_const_mod, only : SHR_CONST_CDAY
  use abortutils   , only : endrun
  use UrbanMod     , only : urban_hac, urban_traffic

#if (defined CN) && (defined VERTSOILC)
  use CNSoilLittVertTranspMod, only: som_diffus, som_adv_flux, cryoturb_diffusion_k, max_altdepth_cryoturbation, max_depth_cryoturb

  use CNVerticalProfileMod, only: exponential_rooting_profile, rootprof_exp, surfprof_exp, pftspecific_rootingprofile

#ifndef CENTURY_DECOMP
  use CNDecompCascadeMod_BGC, only: decomp_depth_efolding, froz_q10
#else
  use CNDecompCascadeMod_CENTURY, only: decomp_depth_efolding, froz_q10
#endif

#endif         

#if (defined CN) && (defined NITRIF_DENITRIF)
  use CNNitrifDenitrifMod, only: no_frozen_nitrif_denitrif
#endif

#if (defined CN)
  !!! C14
  use CNC14DecayMod, only: use_c14_bombspike, atm_c14_filename
#endif


  use SurfaceAlbedoMod, only : albice
#ifdef CN
  use CNAllocationMod , only : suplnitro
#endif
  use clm_nlUtilsMod  , only : find_nlgroup_name
  use Hydrology1Mod   , only : Hydrology1_readnl
  use SoilHydrologyMod, only : SoilHydrology_readnl
 
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
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !PRIVATE TYPES:
! Namelist variables only used locally
  character(len=  7) :: runtyp(4)                        ! run type
  character(len=SHR_KIND_CL) :: NLFilename = 'lnd.stdin' ! Namelist filename
#if (defined _OPENMP)
   integer, external :: omp_get_max_threads  ! max number of threads that can execute
                                             ! concurrently in a single parallel region
#endif
!EOP
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_setNL
!
! !INTERFACE:
  subroutine control_setNL( NLfile )
!
! !USES:
  use clm_varctl , only : NLFileName_in
!
!
    implicit none
!
! !DESCRIPTION:
! Set the namelist filename to use
!
! !ARGUMENTS:
  character(len=*), intent(IN) :: NLFile ! Namelist filename
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=32) :: subname = 'control_setNL'  ! subroutine name
    logical :: lexist                               ! File exists

    ! Error checking...
    if ( len_trim(NLFile) == 0 )then
       call endrun( subname//' error: nlfilename entered is not set' )
    end if
    inquire (file = trim(NLFile), exist = lexist)
    if ( .not. lexist )then
       call endrun( subname//' error: NLfilename entered does NOT exist:'//trim(NLFile) )
    end if
    if ( len_trim(NLFile) > len(NLFilename) )then
       call endrun( subname//' error: entered NLFile is too long' )
    end if
    ! Set the filename
    NLFilename = NLFile
    NLFilename_in = NLFilename   ! For use in external namelists and to avoid creating dependencies on controlMod
  end subroutine control_setNL

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_init
!
! !INTERFACE:
  subroutine control_init( )
!
! !DESCRIPTION:
! Initialize CLM run control information
!
! !USES:
    use clm_time_manager , only : set_timemgr_init, is_perpetual, get_timemgr_defaults
    use fileutils        , only : getavu, relavu
    use shr_string_mod   , only : shr_string_getParentDir
    use clm_varctl       , only : clmvarctl_init, set_clmvarctl, nsrBranch, nsrStartup, &
                                  nsrContinue
    use clm_cpl_indices  , only : glc_nec

    implicit none
!
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
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
         fatmlndfrc, finidat, nrevsn

    ! Input datasets

    namelist /clm_inparm/  &
         fsurdat, fatmtopo, flndtopo, &
         fpftcon, fpftdyn,  fsnowoptics, fsnowaging

    ! History, restart options

    namelist /clm_inparm/  &
         hist_empty_htapes, hist_dov2xy, &
         hist_avgflag_pertape, hist_type1d_pertape, &
         hist_nhtfrq,  hist_ndens, hist_mfilt, &
         hist_fincl1,  hist_fincl2, hist_fincl3, &
         hist_fincl4,  hist_fincl5, hist_fincl6, &
         hist_fexcl1,  hist_fexcl2, hist_fexcl3, &
         hist_fexcl4,  hist_fexcl5, hist_fexcl6, &
         outnc_large_files
#ifdef LCH4
    namelist /clm_inparm/ hist_wrtch4diag
#endif

    ! BGC info

#if (defined CN)
    namelist /clm_inparm/  &
         suplnitro
    namelist /clm_inparm/ &
         nfix_timeconst
    namelist /clm_inparm/ &
         spinup_state, override_bgc_restart_mismatch_dump
#endif

    namelist /clm_inparm / &
         co2_type

    namelist /clm_inparm / perchroot, perchroot_alt
#ifdef LCH4
    namelist /clm_inparm / anoxia, anoxia_wtsat
#endif
    namelist /clm_inparm / deepmixing_depthcrit, deepmixing_mixfact,  lake_melt_icealb
                                                                      ! lake_melt_icealb is of dimension numrad

    ! Glacier_mec info
    namelist /clm_inparm / &    
         maxpatch_glcmec, glc_smb, glc_dyntopo, glc_grid, fglcmask 

    ! Other options

    namelist /clm_inparm/  &
         clump_pproc, wrtdia, pertlim, &
         create_crop_landunit, nsegspc, co2_ppmv, override_nsrest, &
         albice, more_vertlayers, subgridflag, irrigate
    ! Urban options

    namelist /clm_inparm/  &
         urban_hac, urban_traffic

#if (defined CN) && (defined VERTSOILC)
    ! vertical soil mixing variables
    namelist /clm_inparm/  &
         som_diffus, som_adv_flux, cryoturb_diffusion_k, max_altdepth_cryoturbation, max_depth_cryoturb

    ! depth inhibition of decomposition paramters
    namelist /clm_inparm/  &
         decomp_depth_efolding, froz_q10

    ! C and N input vertical profiles
    namelist /clm_inparm/  & 
          exponential_rooting_profile, rootprof_exp, surfprof_exp, pftspecific_rootingprofile
#endif

#if (defined CN) && (defined NITRIF_DENITRIF)
  namelist /clm_inparm / no_frozen_nitrif_denitrif
#endif

    namelist /clm_inparm / use_c13, use_c14

#if (defined CN)
    !!! C14
    namelist /clm_inparm/  &
         use_c14_bombspike, atm_c14_filename
#endif


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
          call endrun( subname//' error: nlfilename not set' )
       end if
       unitn = getavu()
       write(iulog,*) 'Read in clm_inparm namelist from: ', trim(NLFilename)
       open( unitn, file=trim(NLFilename), status='old' )
       call find_nlgroup_name(unitn, 'clm_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading clm_inparm namelist')
          end if
       end if
       call relavu( unitn )

       ! ----------------------------------------------------------------------
       ! Consistency checks on input namelist.
       ! ----------------------------------------------------------------------

       call set_timemgr_init( dtime_in=dtime )

#if (defined CNDV)
       if (is_perpetual()) then
          write(iulog,*)'RTM or CNDV cannot be defined in perpetual mode'
          call endrun()
       end if
#endif
       if (is_perpetual()) then
          if (finidat == ' ') then
             write(iulog,*)'must specify initial dataset for perpetual mode'
             call endrun()
          end if
       end if

       if (urban_traffic) then
          write(iulog,*)'Urban traffic fluxes are not implemented currently'
          call endrun()
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
              call endrun( subname//' ERROR: can ONLY override clm start-type ' // &
                           'to branch type and ONLY if driver is a startup type' )
           end if
           call set_clmvarctl( nsrest_in=override_nsrest )
       end if
       
       if (maxpatch_glcmec > 0) then
          create_glacier_mec_landunit = .true.
       else
          create_glacier_mec_landunit = .false.
       end if
       
    endif   ! end of if-masterproc if-block

    call clmvarctl_init( masterproc, dtime )

    ! ----------------------------------------------------------------------
    ! Read in other namelists for other modules
    ! ----------------------------------------------------------------------

    call Hydrology1_readnl(    NLFilename )
    call SoilHydrology_readnl( NLFilename )


    ! ----------------------------------------------------------------------
    ! Broadcast all control information if appropriate
    ! ----------------------------------------------------------------------

    call control_spmd()
    
    if (masterproc) then
       write(iulog,*) 'Successfully initialized run control settings'
       write(iulog,*)
    endif

  end subroutine control_init


!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_spmd
!
! !INTERFACE:
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
    use clm_varctl, only : single_column, scmlat, scmlon, rpntfil
    use clm_varpar, only : numrad
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer ier       !error code
!-----------------------------------------------------------------------

    ! run control variables

    call mpi_bcast (caseid,         len(caseid),        MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (ctitle,         len(ctitle),        MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (model_version,  len(model_version), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hostname,       len(hostname),      MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (username,       len(username),      MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (nsrest,                     1,      MPI_INTEGER  , 0, mpicom, ier)

    ! initial file variables

    call mpi_bcast (nrevsn  , len(nrevsn)  , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (finidat , len(finidat) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsurdat , len(fsurdat) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fatmlndfrc,len(fatmlndfrc),MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fatmtopo, len(fatmtopo) ,MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (flndtopo, len(flndtopo) ,MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fpftcon , len(fpftcon) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fpftdyn , len(fpftdyn) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsnowoptics,  len(fsnowoptics),  MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsnowaging,   len(fsnowaging),   MPI_CHARACTER, 0, mpicom, ier)

    ! Irrigation

    call mpi_bcast(irrigate,             1, MPI_LOGICAL, 0, mpicom, ier)

    ! Landunit generation

    call mpi_bcast(create_crop_landunit, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast(allocate_all_vegpfts, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! BGC

    call mpi_bcast (co2_type, len(co2_type), MPI_CHARACTER, 0, mpicom, ier)
#ifdef CN
    call mpi_bcast (suplnitro, len(suplnitro), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (nfix_timeconst,             1, MPI_REAL8,     0, mpicom, ier)
    call mpi_bcast (spinup_state,               1, MPI_INTEGER,   0, mpicom, ier)
    call mpi_bcast (override_bgc_restart_mismatch_dump, 1, MPI_LOGICAL,   0, mpicom, ier)
#endif

    ! isotopes
    call mpi_bcast (use_c13,          1, MPI_LOGICAL,     0, mpicom, ier)
    call mpi_bcast (use_c14,          1, MPI_LOGICAL,     0, mpicom, ier)

#if (defined CN) && (defined VERTSOILC)
    ! vertical soil mixing variables
    call mpi_bcast (som_diffus,                 1, MPI_REAL8,     0, mpicom, ier)
    call mpi_bcast (som_adv_flux,               1, MPI_REAL8,     0, mpicom, ier)
    call mpi_bcast (cryoturb_diffusion_k,       1, MPI_REAL8,     0, mpicom, ier)
    call mpi_bcast (max_altdepth_cryoturbation, 1, MPI_REAL8,     0, mpicom, ier)
    call mpi_bcast (max_depth_cryoturb,         1, MPI_REAL8,     0, mpicom, ier)

    ! depth inhibition of decomposition paramters
    call mpi_bcast (decomp_depth_efolding,          1, MPI_REAL8,     0, mpicom, ier)
    call mpi_bcast (froz_q10,                       1, MPI_REAL8,     0, mpicom, ier)

    ! C and N input vertical profiles
    call mpi_bcast (exponential_rooting_profile,          1, MPI_LOGICAL,     0, mpicom, ier)
    call mpi_bcast (rootprof_exp,                         1, MPI_REAL8,     0, mpicom, ier)
    call mpi_bcast (surfprof_exp,                         1, MPI_REAL8,     0, mpicom, ier)
    call mpi_bcast (pftspecific_rootingprofile,           1, MPI_LOGICAL,     0, mpicom, ier)
#endif

#if (defined CN) && (defined NITRIF_DENITRIF)
    call mpi_bcast (no_frozen_nitrif_denitrif,  1, MPI_LOGICAL, 0, mpicom, ier)
#endif

#if (defined CN)
    !!! C14
    call mpi_bcast (use_c14_bombspike,  1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (atm_c14_filename,  len(atm_c14_filename), MPI_CHARACTER, 0, mpicom, ier)
#endif

    call mpi_bcast (perchroot, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (perchroot_alt, 1, MPI_LOGICAL, 0, mpicom, ier)
#ifdef LCH4
    call mpi_bcast (anoxia, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (anoxia_wtsat, 1, MPI_LOGICAL, 0, mpicom, ier)
#endif
! Lakes
    call mpi_bcast (deepmixing_depthcrit,     1, MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (deepmixing_mixfact,       1, MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (lake_melt_icealb,    numrad, MPI_REAL8, 0, mpicom, ier)

    ! physics variables

    call mpi_bcast (urban_hac     , len(urban_hac), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (urban_traffic , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (nsegspc     , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (subgridflag , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (wrtdia      , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (single_column,1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (scmlat,       1, MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast (scmlon,       1, MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast (co2_ppmv    , 1, MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast (albice      , 2, MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast (more_vertlayers,1, MPI_LOGICAL, 0, mpicom, ier)

    ! glacier_mec variables

    call mpi_bcast (create_glacier_mec_landunit, 1, MPI_LOGICAL  , 0, mpicom, ier)
    call mpi_bcast (maxpatch_glcmec             ,1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (glc_smb,                     1, MPI_LOGICAL  , 0, mpicom, ier)
    call mpi_bcast (glc_dyntopo,                 1, MPI_LOGICAL  , 0, mpicom, ier)
    call mpi_bcast (glc_grid,        len(glc_grid), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fglcmask,        len(fglcmask), MPI_CHARACTER, 0, mpicom, ier)

    ! history file variables

    call mpi_bcast (outnc_large_files, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_empty_htapes, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_dov2xy, size(hist_dov2xy), MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_nhtfrq, size(hist_nhtfrq), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_mfilt, size(hist_mfilt), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_ndens, size(hist_ndens), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_avgflag_pertape, size(hist_avgflag_pertape), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_type1d_pertape, max_namlen*size(hist_type1d_pertape), MPI_CHARACTER, 0, mpicom, ier)
#ifdef LCH4
    call mpi_bcast (hist_wrtch4diag, 1, MPI_LOGICAL, 0, mpicom, ier)
#endif
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

    ! error growth perturbation limit
    call mpi_bcast (pertlim, 1, MPI_REAL8, 0, mpicom, ier)

  end subroutine control_spmd

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_print
!
! !INTERFACE:
  subroutine control_print ()
!
! !DESCRIPTION:
! Write out the clm namelist run control variables
!
! !USES:
!
    use clm_varctl,      only : source, rpntdir, rpntfil, nsrStartup, nsrBranch, &
                                nsrContinue
#ifdef CN
    use CNAllocationMod, only : suplnitro, suplnNon
#endif
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer i  !loop index
    character(len=32) :: subname = 'control_print'  ! subroutine name
!------------------------------------------------------------------------

    write(iulog,*) 'define run:'
    write(iulog,*) '   source                = ',trim(source)
    write(iulog,*) '   model_version         = ',trim(model_version)
    write(iulog,*) '   run type              = ',runtyp(nsrest+1)
    write(iulog,*) '   case title            = ',trim(ctitle)
    write(iulog,*) '   username              = ',trim(username)
    write(iulog,*) '   hostname              = ',trim(hostname)
    write(iulog,*) 'input data files:'
    write(iulog,*) '   PFT physiology = ',trim(fpftcon)
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
#ifdef CN
    if (suplnitro /= suplnNon)then
        write(iulog,*) '   Supplemental Nitrogen mode is set to run over PFTs: ', &
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
    else
       call endrun( subname//' error: spinup_state can only have integer value of 0 or 1' )
    end if

    write(iulog,*) '   override_bgc_restart_mismatch_dump                     : ', override_bgc_restart_mismatch_dump
#endif

#if (defined CN) && (defined VERTSOILC)
    write(iulog, *) '   som_adv_flux, the advection term in soil mixing (m/s) : ', som_adv_flux
    write(iulog, *) '   som_diffus, the diffusion term in soil mixing (m/s^2) : ', som_diffus
    write(iulog, *) '   cryoturb_diffusion_k  (m/s^2)                         : ', cryoturb_diffusion_k
    write(iulog, *) '   max_altdepth_cryoturbation (m)                        : ', max_altdepth_cryoturbation
    write(iulog, *) '   max_depth_cryoturb (m)                                : ', max_depth_cryoturb

    write(iulog, *) '   decomp_depth_efolding                                 : ', decomp_depth_efolding
    write(iulog, *) '   froz_q10                                              : ', froz_q10

    write(iulog, *) '   exponential_rooting_profile                           : ', exponential_rooting_profile
    write(iulog, *) '   rootprof_exp                                          : ', rootprof_exp
    write(iulog, *) '   surfprof_exp                                          : ', surfprof_exp
    write(iulog, *) '   pftspecific_rootingprofile                            : ', pftspecific_rootingprofile
#endif

#if (defined CN) && (defined NITRIF_DENITRIF)
    write(iulog, *) '   no_frozen_nitrif_denitrif                             : ', no_frozen_nitrif_denitrif
#endif

#if (defined CN)
    write(iulog, *) '  use_c13                                                : ', use_c13
    write(iulog, *) '  use_c14                                                : ', use_c14
    !!! C14
    write(iulog, *) '  use_c14_bombspike                                      : ', use_c14_bombspike
    write(iulog, *) '  atm_c14_filename                                       : ', atm_c14_filename
#endif

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
       if (glc_dyntopo) then
          write(iulog,*) '   glc CLM glacier topography will evolve dynamically'
       else
          write(iulog,*) '   glc CLM glacier topography will NOT evolve dynamically'
       endif
       if (glc_smb) then
          write(iulog,*) '   glc surface mass balance will be passed to ice sheet model'
       else
          write(iulog,*) '   glc positive-degree-day info will be passed to ice sheet model'
       endif
    endif

    if (nsrest == nsrStartup .and. finidat == ' ') write(iulog,*) '   initial data created by model'
    if (nsrest == nsrStartup .and. finidat /= ' ') write(iulog,*) '   initial data   = ',trim(finidat)
    if (nsrest /= nsrStartup) write(iulog,*) '   restart data   = ',trim(nrevsn)
    write(iulog,*) '   atmospheric forcing data is from cesm atm model'
    write(iulog,*) 'Restart parameters:'
    write(iulog,*)'   restart pointer file directory     = ',trim(rpntdir)
    write(iulog,*)'   restart pointer file name          = ',trim(rpntfil)
    if ( outnc_large_files ) then
       write(iulog,*)'Large file support for output files is ON'
    end if
    write(iulog,*) 'model physics parameters:'
#if (defined PERGRO)
    write(iulog,*) '   flag for random perturbation test is set'
#else
    write(iulog,*) '   flag for random perturbation test is not set'
#endif
    write(iulog,*) '   CO2 volume mixing ratio   (umol/mol)   = ', co2_ppmv
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
    if ( pertlim /= 0.0_r8 ) &
    write(iulog,*) '   perturbation limit   = ',pertlim
    write(iulog,*) '   maxpatch_pft         = ',maxpatch_pft
    write(iulog,*) '   allocate_all_vegpfts = ',allocate_all_vegpfts
    write(iulog,*) '   nsegspc              = ',nsegspc
! New fields
    write(iulog,*) ' perchroot (plant water stress based on unfrozen layers only) = ',perchroot
    write(iulog,*) ' perchroot (plant water stress based on time-integrated active layer only) = ',perchroot
#ifdef LCH4
    write(iulog,*) ' anoxia (applied to soil decomposition)             = ',anoxia
    write(iulog,*) ' anoxia_wtsat (weight anoxia by inundated fraction) = ',anoxia_wtsat
#endif
! Lakes
    write(iulog,*)
    write(iulog,*) 'Lake Model Namelists:'
    write(iulog,*) 'Increased mixing relative to Hostetler wind-driven eddy expression ',&
                   'will be used for deep lakes exceeding depth ', deepmixing_depthcrit,&
                      ' by a factor of ', deepmixing_mixfact, '.'
    write(iulog,*) 'Albedo over melting lakes will approach values (visible, NIR):', lake_melt_icealb, &
                   'as compared with 0.60, 0.40 for cold frozen lakes with no snow.'

  end subroutine control_print

end module controlMod
