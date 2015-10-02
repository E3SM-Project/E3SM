module clm_varctl

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varctl
!
! !DESCRIPTION:
! Module containing run control variables
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: clm_varctl_set     ! Set variables
  public :: clmvarctl_init    ! Initialize and check values after namelist input

  private
  save
!
! !PUBLIC TYPES:
!
  integer, parameter, private ::  iundef = -9999999
  integer, parameter, private ::  rundef = -9999999._r8
!
! Run control variables
!
  character(len=256), public :: caseid  = ' '                            ! case id
  character(len=256), public :: ctitle  = ' '                            ! case title
  integer, public :: nsrest             = iundef                         ! Type of run
  integer, public, parameter :: nsrStartup  = 0                          ! Startup from initial conditions
  integer, public, parameter :: nsrContinue = 1                          ! Continue from restart files
  integer, public, parameter :: nsrBranch   = 2                          ! Branch from restart files
  logical, public :: brnch_retain_casename = .false.                     ! true => allow case name to remain the same for branch run
                                                                         ! by default this is not allowed
  logical, public :: noland = .false.                                    ! true => no valid land points -- do NOT run
  character(len=256), public :: hostname = ' '                           ! Hostname of machine running on
  character(len=256), public :: username = ' '                           ! username of user running program
  character(len=256), public :: source   = "Community Land Model CLM4.0" ! description of this source
  character(len=256), public :: version  = " "                           ! version of program
  character(len=256), public :: conventions = "CF-1.0"                   ! dataset conventions
!
! Unit Numbers
!
  integer, public :: iulog = 6        ! "stdout" log file unit number, default is 6
!
! Output NetCDF files
!
  logical, public :: outnc_large_files = .true.         ! large file support for output NetCDF files
!
! Run input files
!
  character(len=256), public :: finidat    = ' '        ! initial conditions file name
  character(len=256), public :: fsurdat    = ' '        ! surface data file name
  character(len=256), public :: fatmgrid   = ' '        ! atm grid file name
  character(len=256), public :: fatmlndfrc = ' '        ! lnd frac file on atm grid
  character(len=256), public :: fatmtopo   = ' '        ! topography on atm grid
  character(len=256), public :: flndtopo   = ' '        ! topography on lnd grid
  character(len=256), public :: flanduse_timeseries    = ' '        ! dynamic landuse dataset
  character(len=256), public :: fpftcon    = ' '        ! ASCII data file with PFT physiological constants
  character(len=256), public :: nrevsn     = ' '        ! restart data file name for branch run
  character(len=256), public :: fsnowoptics  = ' '      ! snow optical properties file name
  character(len=256), public :: fsnowaging   = ' '      ! snow aging parameters file name

!
! Landunit logic
!
  logical, public :: create_crop_landunit = .false.     ! true => separate crop landunit is not created by default
  logical, public :: allocate_all_vegpfts = .false.     ! true => allocate memory for all possible vegetated pfts on
                                                        ! vegetated landunit if at least one pft has nonzero weight
!
! BGC logic and datasets
!
  character(len=16), public :: co2_type = 'constant'    ! values of 'prognostic','diagnostic','constant'
!
! Physics
!
  logical,  public :: wrtdia       = .false.            ! true => write global average diagnostics to std out
  real(r8), public :: co2_ppmv     = 355._r8            ! atmospheric CO2 molar ratio (by volume) (umol/mol)

 ! C isotopes
  logical, public :: use_c13 = .false.                  ! true => use C-13 model
  logical, public :: use_c14 = .false.                  ! true => use C-14 model

! glacier_mec control variables: default values (may be overwritten by namelist)
! NOTE: glc_smb must have the same values for CLM and GLC

  logical , public :: create_glacier_mec_landunit = .false. ! glacier_mec landunit is not created (set in controlMod)
  logical , public :: glc_smb = .true.                      ! if true, pass surface mass balance info to GLC
                                                            ! if false, pass positive-degree-day info to GLC
  logical , public :: glc_dyntopo = .false.                 ! true => CLM glacier topography changes dynamically
  real(r8), public, allocatable :: glc_topomax(:)           ! upper limit of each class (m)  (set in surfrd)
  character(len=256), public :: glc_grid = ' '              ! glc_grid used to determine fglcmask  
  character(len=256), public :: fglcmask = ' '              ! glacier mask file name (based on glc_grid)
!
! single column control variables
!
  logical,  public :: single_column = .false.           ! true => single column mode
  real(r8), public :: scmlat        = rundef            ! single column lat
  real(r8), public :: scmlon        = rundef            ! single column lon
!
! instance control
!
  integer, public :: inst_index
  character(len=16), public :: inst_name
  character(len=16), public :: inst_suffix
!
! Decomp control variables
!
  integer, public :: nsegspc = 20                       ! number of segments per clump for decomp
!
! Derived variables (run, history and restart file)
!
  character(len=256), public :: rpntdir = '.'            ! directory name for local restart pointer file
  character(len=256), public :: rpntfil = 'rpointer.lnd' ! file name for local restart pointer file
!
! Migration of CPP variables
! 
#if (defined CN)
  logical, public :: use_cn = .true.
#else
  logical, public :: use_cn = .false.
#endif
#if (defined CNDV)
  logical, public :: use_cndv = .true.
#else
  logical, public :: use_cndv = .false.
#endif
#if (defined CROP)
  logical, public :: use_crop = .true.
#else
  logical, public :: use_crop = .false.
#endif
#if (defined SNICAR_FRC)
  logical, public :: use_snicar_frc = .true.
#else
  logical, public :: use_snicar_frc = .false.
#endif
#if (defined NOFIRE)
  logical, public :: use_nofire = .true.
#else
  logical, public :: use_nofire = .false.
#endif
#if (defined VANCOUVER)
  logical, public :: use_vancouver = .true.
#else
  logical, public :: use_vancouver = .false.
#endif
#if (defined MEXICOCITY)
  logical, public :: use_mexicocity = .true.
#else
  logical, public :: use_mexicocity = .false.
#endif
#if (defined AD_SPINUP) 
  logical, public :: use_ad_spinup = .true.
#else
  logical, public :: use_ad_spinup = .false.
#endif
#if (defined EXIT_SPINUP) 
  logical, public :: use_exit_spinup = .true.
#else
  logical, public :: use_exit_spinup = .false.
#endif
  !needed for compatibility with changes in clm4_5 and reference dy lnd_comp_mct
  !however use_voc is not used anywhere inside the clm4_0 code
  logical, public :: use_voc = .true. 
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein and Gordon Bonan
! 1 June 2004, Peter Thornton: added fnedpdat for nitrogen deposition data
!
!EOP
!-----------------------------------------------------------------------
  logical, private :: clmvarctl_isset = .false.

!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_varctl_set
!
! !INTERFACE:
  subroutine clm_varctl_set( caseid_in, ctitle_in, brnch_retain_casename_in,    &
                            single_column_in, scmlat_in, scmlon_in, nsrest_in, &
                            version_in, hostname_in, username_in)
!
! !DESCRIPTION:
!      Set input control variables.
!
! !USES:
  use shr_sys_mod, only : shr_sys_abort
!
! !ARGUMENTS:
  character(len=256), optional, intent(IN) :: caseid_in                ! case id
  character(len=256), optional, intent(IN) :: ctitle_in                ! case title
  logical,            optional, intent(IN) :: brnch_retain_casename_in ! true => allow case name to remain the same for branch run
  logical,            optional, intent(IN) :: single_column_in         ! true => single column mode
  real(r8),           optional, intent(IN) :: scmlat_in                ! single column lat
  real(r8),           optional, intent(IN) :: scmlon_in                ! single column lon
  integer,            optional, intent(IN) :: nsrest_in                ! 0: initial run. 1: restart: 3: branch
  character(len=256), optional, intent(IN) :: version_in               ! model version
  character(len=256), optional, intent(IN) :: hostname_in              ! hostname running on
  character(len=256), optional, intent(IN) :: username_in              ! username running job

!
! !LOCAL VARIABLES:
   character(len=32) :: subname = 'clm_varctl_set'  ! subroutine name
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!EOP
!-----------------------------------------------------------------------
    if ( clmvarctl_isset )then
       call shr_sys_abort( subname//' ERROR:: control variables already set -- can not call this subroutine' )
    end if
    if ( present(caseid_in       ) ) caseid        = caseid_in
    if ( present(ctitle_in       ) ) ctitle        = ctitle_in
    if ( present(single_column_in) ) single_column = single_column_in
    if ( present(scmlat_in       ) ) scmlat        = scmlat_in
    if ( present(scmlon_in       ) ) scmlon        = scmlon_in
    if ( present(nsrest_in       ) ) nsrest        = nsrest_in
    if ( present(brnch_retain_casename_in) ) brnch_retain_casename = brnch_retain_casename_in
    if ( present(version_in      ) ) version       = version_in
    if ( present(username_in     ) ) username      = username_in
    if ( present(hostname_in     ) ) hostname      = hostname_in

  end subroutine clm_varctl_set

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clmvarctl_init
!
! !INTERFACE:
  subroutine clmvarctl_init( masterproc, dtime )
!
! !DESCRIPTION:
!      Check that values are correct, and finish setting variables based on other variables.
!
! !USES:
  use shr_sys_mod  , only : shr_sys_abort
  use clm_varpar   , only : maxpatch_pft, numpft
!
! !ARGUMENTS:
  logical, intent(IN) :: masterproc  ! proc 0 logical for printing msgs
  integer, intent(IN) :: dtime       ! timestep in seconds
!
! !LOCAL VARIABLES:
   character(len=32) :: subname = 'clmvarctl_init'  ! subroutine name
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!EOP
!-----------------------------------------------------------------------

    ! landunit generation

    if (maxpatch_pft == numpft+1) then
       allocate_all_vegpfts = .true.
    else
       allocate_all_vegpfts = .false.
       if (use_crop) then
          write(iulog,*)'maxpatch_pft = ',maxpatch_pft,&
               ' does NOT equal numpft+1 = ',numpft+1
          call shr_sys_abort( subname//' ERROR:: Can NOT turn CROP on without all PFTs' )
       end if
    end if

    if (masterproc) then

       ! Consistency settings for co2 type

       if (co2_type /= 'constant' .and. co2_type /= 'prognostic' .and. co2_type /= 'diagnostic') then
          write(iulog,*)'co2_type = ',co2_type,' is not supported'
          call shr_sys_abort( subname//' ERROR:: choices are constant, prognostic or diagnostic' )
       end if

       ! Consistency settings for dynamic land use, etc.

       if (flanduse_timeseries /= ' ' .and. create_crop_landunit) &
          call shr_sys_abort( subname//' ERROR:: dynamic landuse is currently not supported with create_crop_landunit option' )
       if (create_crop_landunit .and. .not.allocate_all_vegpfts) &
          call shr_sys_abort( subname//' ERROR:: maxpft<numpft+1 is currently not supported with create_crop_landunit option' )
       if (flanduse_timeseries /= ' ') then
          if (use_cndv) then
             call shr_sys_abort( subname//' ERROR:: dynamic landuse is currently not supported with CNDV option' )
          end if
       end if

       ! Check on run type
       if (nsrest == iundef) call shr_sys_abort( subname//' ERROR:: must set nsrest' )
       if (nsrest == nsrBranch .and. nrevsn == ' ') &
          call shr_sys_abort( subname//' ERROR: need to set restart data file name' )

       ! Model physics

       if ( (co2_ppmv <= 0.0_r8) .or. (co2_ppmv > 3000.0_r8) ) &
          call shr_sys_abort( subname//' ERROR: co2_ppmv is out of a reasonable range' )

       if (nsrest == nsrStartup ) nrevsn = ' '
       if (nsrest == nsrContinue) nrevsn = 'set by restart pointer file file'
       if (nsrest /= nsrStartup .and. nsrest /= nsrContinue .and. nsrest /= nsrBranch ) &
          call shr_sys_abort( subname//' ERROR: nsrest NOT set to a valid value' )

       if ( single_column .and. (scmlat == rundef  .or. scmlon == rundef ) ) &
          call shr_sys_abort( subname//' ERROR:: single column mode on -- but scmlat and scmlon are NOT set' )

    endif   ! end of if-masterproc if-block

    clmvarctl_isset = .true.

  end subroutine clmvarctl_init

end module clm_varctl
