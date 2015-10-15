!===============================================================================
! SVN $Id$
! SVN $URL$
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_inputInfo_mod --- Module for input data shared between CCSM components
!
! !DESCRIPTION:
!
!     A module to create a derived type and a namelist read for data that will
!  be shared between CCSM components. This data is also written out to a NetCDF restart
!  file and read in when restarting.
!
! Typical usage:
!
! ! Set the default values of the object...
! shr_inputInfo_initSetDefaults( CCSMInit )
! ! Change some of the values in the object...
! shr_inputInfo_initPutData( CCSMInit, start_type="continue", aqua_planet=.true. )
! ! Set data for the object based on input namelist
! shr_inputInfo_initRead( NLFileName, LogPrint=.true., mpicom=mpicom, &
!                         MasterTask=mastertask, CCSMInitOut=CCSMInit )
! ! Print the object out
! shr_inputInfo_initPrint( CCSMInit )
! ! Get data from the object
! shr_inputInfo_initGetData( CCSMInit, case_name=case_name, case_desc=case_desc )
! ! Check some of the logical values...
! if ( shr_inputInfo_initIsRestart( CCSMInit ) ) ....
! if ( shr_inputInfo_initIsContinue( CCSMInit ) ) ....
! if ( shr_inputInfo_initIsBranch( CCSMInit ) ) ....
! ! Test if a given model should run based on the contents of the object
! if ( shr_inputInfo_initIsModelToRun( CCSMInit, model="atm" ) ) ....
!
! ! Read in restart pointer file and then restart file for inputinfo...
! shr_inputInfo_initRPointerRead( restart_file, MPICom, MasterTask, &
!                                 CCSMInitOut )  ! Also gets file from archive
! shr_inputInfo_initRestRead( restart_file, MPICom, MasterTask, &
!                             CCSMInitOut )
! .... read any other object data from this same driver restart file
!
! ! Now write out a driver restart file
! ! First write out the driver restart pointer file and get file-name with time-stamp
! shr_inputInfo_initRPointerWrite( ymd, tod, MPICom, MasterTask, &
!                                  CCSMInit, restart_file )
! shr_inputInfo_initRestWrite( restart_file, MPICom, MasterTask, &
!                              CCSMInit )
! .... write any other object data to this same driver restart file
! ! Archive the file to MSS
! shr_inputInfo_initRestArchive( MasterTask, CCSMInit, async=.true., remove=.true., &
!                                restart_file=restart_file )
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2005-Nov-11 - E. Kluzek - creation
!
! !INTERFACE: ------------------------------------------------------------------

MODULE shr_inputInfo_mod

! !USES:

   USE shr_kind_mod,      ONLY : SHR_KIND_CS, SHR_KIND_CL, SHR_KIND_IN,      &
                                 SHR_KIND_R8, SHR_KIND_I8
   USE shr_sys_mod,       ONLY : shr_sys_flush, shr_sys_abort, shr_sys_getenv
   USE shr_ncio_mod,      ONLY : shr_ncio_descripType, shr_ncio_open,        &
                                 shr_ncio_close

   implicit none

   private  ! default private

! !PUBLIC TYPES:

   public :: shr_inputInfo_initType        ! Derived type for shared data

! !PUBLIC MEMBER FUNCTIONS

   public :: shr_inputInfo_initSetDefaults     ! Set default values
   public :: shr_inputInfo_initPrint           ! Print out values
   public :: shr_inputInfo_initGetData         ! Get values from object
   public :: shr_inputInfo_initPutData         ! Change values
   public :: shr_inputInfo_initRead            ! Input namelist read method
   public :: shr_inputInfo_initIsStartup       ! Return if startup case or not
   public :: shr_inputInfo_initIsRestart       ! Return if restart (branch OR continue)
   public :: shr_inputInfo_initIsContinue      ! Return if continue case or not
   public :: shr_inputInfo_initIsBranch        ! Return if branch case or not
   public :: shr_inputInfo_initIsSame          ! Check that the two input objects are same
   public :: shr_inputInfo_initRestRead        ! Read in restart information
   public :: shr_inputInfo_initRestWrite       ! Write out restart information
   public :: shr_inputInfo_initIsModelToRun    ! Run model component
   public :: shr_inputInfo_initRPointerRead    ! Read restart pointer file get filename
   public :: shr_inputInfo_initRPointerWrite   ! Write restart pointer file get filename
   public :: shr_inputInfo_initRestArchive     ! Archive the restart file
   public :: shr_inputInfo_initCheck           ! Check object for consistency/validity

! !PUBLIC DATA MEMBERS:

   ! Shared namelist name
   character(len=*), public, parameter :: &
                     shr_inputInfo_NLFilename = "namelist"
   ! Strings of valid start_type options
   character(len=*), public, parameter :: &
                     shr_inputInfo_start_type_start = "startup"
   character(len=*), public, parameter :: &
                     shr_inputInfo_start_type_cont  = "continue"
   character(len=*), public, parameter :: &
                     shr_inputInfo_start_type_brnch = "branch"

!EOP

   integer, parameter :: NVars = 10       ! Number of variables saved on the restart file
   ! Names of variables from InitType to save on restart file
   character(len=*), parameter :: shr_inputInfo_saveRestart(NVars) = (/ &
      "case_desc      ", &
      "case_name      ", &
      "atm_adiabatic  ", &
      "atm_ideal_phys ", &
      "aqua_planet    ", &
      "single_column  ", &
      "scmlat         ", &
      "scmlon         ", &
      "mss_irt        ", &
      "mss_wpass      " /)
   !
   ! InputInfo derived type
   !
   type shr_inputInfo_initType
      private     ! This type is opaque
      !-----------------------------------------------------------------------------
      ! Information shared by all model components
      !-----------------------------------------------------------------------------
      ! ----------- Information not saved on restart --------------------
      ! Start type: startup, branch or continue
      character(SHR_KIND_CL) :: start_type      ! Type of startup
      character(SHR_KIND_CS) :: logFilePostFix  ! postfix for output log files
      character(SHR_KIND_CS) :: outPathRoot     ! root for output log files

      ! ----------- Information saved on restart ------------------------
      ! Case description and name (save these on restart)
      character(SHR_KIND_CL)  :: case_desc      ! Long description of this case
      character(SHR_KIND_CS)  :: case_name      ! Short case identification
      ! Special configurations (save these on restart)
      logical                 :: atm_adiabatic  ! No surface models and atm adiabatic mode
      logical                 :: atm_ideal_phys ! No surface models and atm ideal-physics
      logical                 :: aqua_planet    ! No ice/lnd, analytic ocn, perpetual time
      logical                 :: single_column  ! single column mode
      real (SHR_KIND_R8)      :: scmlat         ! single column lat
      real (SHR_KIND_R8)      :: scmlon         ! single column lon
      ! MSS information needed for CAM now, but not in the future? (save these on restart)
      integer                 :: mss_irt        ! MSS retention period
      character(SHR_KIND_CL)  :: mss_wpass      ! MSS write password

      !-----------------------------------------------------------------------------
      ! Information ONLY used by the top level application driver.
      !-----------------------------------------------------------------------------
      ! ----------- Information not saved on restart --------------------
      ! Restart file info
      character(SHR_KIND_CL)  :: restart_pfile       ! Restart pointer file
      character(SHR_KIND_CL)  :: restart_file        ! Full archive path to restart file
      character(SHR_KIND_CL)  :: archive_dir         ! Directory to archive to...
      character(SHR_KIND_CL)  :: restart_file_override ! List to override on restart file
      logical :: brnch_retain_casename               ! If branch and can use same casename

      !-----------------------------------------------------------------------------
      ! Variables that are written out to restart file from above
      !-----------------------------------------------------------------------------
      type(shr_ncio_descripType) :: var(NVars)
   end type shr_inputInfo_InitType
   !
   ! List of values read in from restart_file that can be over-ridden on namelist
   !
   character(len=*), private, parameter :: restartFileOverRideList = &
              "mss_irt"    &    ! Archival retention period
           //":mss_wpass"  &    ! Archival write password
           //":case_name"  &    ! Case name (only for branch)
           //":case_desc"       ! Case description
   !
   ! Prefix to put in front of variable names on restart file
   !
   character(len=*), parameter :: prefix = "shr_inputInfo_"
   character(len=*), parameter :: LogPrefix = "(shr_inputInfo) "
   character(len=*), parameter :: FAA = "('"//LogPrefix//"',A)"

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
!===============================================================================
! !IROUTINE: shr_inputInfo_initSetDefaults -- Set defaults for shr_CCSMInit data
!   
! !DESCRIPTION:
!   
!     Set the defaults for the input CCSMInit object.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_inputInfo_initSetDefaults( CCSMInit )

! !USES:

   USE shr_ncio_mod,   ONLY: shr_ncio_descripInit, shr_ncio_descripSetDefault
   USE shr_string_mod, ONLY: shr_string_listGetNum, shr_string_listIntersect, &
                             shr_string_listGetName

   implicit none
    
! !INPUT/OUTPUT PARAMETERS:

   type(shr_inputInfo_InitType), intent(OUT) :: CCSMInit   ! Output CCSM structure

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(shr_inputInfo_initSetDefaults) '
    character(len=SHR_KIND_CL) :: restartList ! Items in restart-list and override
    character(len=SHR_KIND_CL) :: name        ! name of list item
    integer :: n                              ! Index
    integer :: rc                             ! Return code
    logical :: init                           ! If MPI initialized or not

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Set default values for each item
    !---------------------------------------------------------------------------
    CCSMInit%case_desc             = ' '
    CCSMInit%case_name             = ' '
    CCSMInit%start_type            = ' '
    CCSMInit%archive_dir           = ' '
    CCSMInit%atm_ideal_phys        = .false.
    CCSMInit%atm_adiabatic         = .false.
    CCSMInit%aqua_planet           = .false.
    CCSMInit%single_column         = .false.
    CCSMInit%scmlat                = -999.
    CCSMInit%scmlon                = -999.
    CCSMInit%restart_pfile         = ' '
    CCSMInit%restart_file          = ' '
    CCSMInit%restart_file_override = ' '
    CCSMInit%mss_irt               = 365
    CCSMInit%mss_wpass             = ' '
    CCSMInit%brnch_retain_casename = .false.
    CCSMInit%logFilePostFix        = '.log'
    CCSMInit%outPathRoot           = './'
    !---------------------------------------------------------------------------
    ! Prepare names of data to write to restart file
    !---------------------------------------------------------------------------
    call shr_ncio_descripSetDefault( NVars, CCSMInit%var )
    do n = 1, NVars
       name = shr_inputInfo_saveRestart(n)
       if (      trim(name) == "atm_adiabatic" )then
          call shr_ncio_descripInit( CCSMInit%var(n), name,                    &
                                 LongName="Atmosphere runs in adiabatic mode", &
                                 LogicalData=.true. )
       else if ( trim(name) == "atm_ideal_phys" )then
          call shr_ncio_descripInit( CCSMInit%var(n), name,                    &
                       LongName="Atmosphere runs in idealized physics mode",   &
                       LogicalData=.true. )
       else if ( trim(name) == "aqua_planet" )then
          call shr_ncio_descripInit( CCSMInit%var(n), name,                    &
                           LongName="Aqua_planet mode (no-land/ice, "//        &
                           "perpetual-date, and analytic SST)",                &
                           LogicalData=.true. )
       else if ( trim(name) == "single_column" )then
          call shr_ncio_descripInit( CCSMInit%var(n), name,                    &
                           LongName="CCSM run in Single Column Mode",           &
                           LogicalData=.true. )
       else if ( trim(name) == "scmlat" )then
          call shr_ncio_descripInit( CCSMInit%var(n), name,                    &
                           LongName="Single Column Mode Latitude",           &
                           units="degrees_north", RealR8Data=.true. )
       else if ( trim(name) == "scmlon" )then
          call shr_ncio_descripInit( CCSMInit%var(n), name,                    &
                           LongName="Single Column Mode Longitude",           &
                           units="degrees_east", RealR8Data=.true. )
       else if ( trim(name) == "mss_irt" )then
          call shr_ncio_descripInit( CCSMInit%var(n), name,                    &
                                     LongName="Mass Store Retention period",   &
                                     units="days", IntegerData=.true. )
       else if ( trim(name) == "case_desc" )then
          call shr_ncio_descripInit( CCSMInit%var(n), name,                    &
                                     LongName="Case description",              &
                                     StringData=.true. )
       else if ( trim(name) == "case_name" )then
          call shr_ncio_descripInit( CCSMInit%var(n), name,                    &
                                     LongName="Case identification name",      &
                                     StringData=.true. )
       else if ( trim(name) == "mss_wpass" )then
          call shr_ncio_descripInit( CCSMInit%var(n), name,                    &
                                     LongName="Mass Store write password",     &
                                     StringData=.true. )
       else
          call shr_sys_abort( subname//': unrecognized variable to store ' // &
                              'on restart file: '//name )
       end if
    end do

END SUBROUTINE shr_inputInfo_initSetDefaults

!===============================================================================
!===============================================================================
! !IROUTINE: shr_inputInfo_initPrint -- Print out values to log file
!   
! !DESCRIPTION:
!   
!     Print derivied type out to screen.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_inputInfo_initPrint( CCSMInit )
    
   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_inputInfo_initType), intent(IN) :: CCSMInit  ! Output CCSM structure

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(shr_inputInfo_initPrint) '
    character(len=*), parameter :: F0A = "('"//LogPrefix//"',A,T45,' = ', A)"
    character(len=*), parameter :: F0L = "('"//LogPrefix//"',A,T45,' = ', L1)"
    character(len=*), parameter :: F0I = "('"//LogPrefix//"',A,T45,' = ', I4)"
    character(len=*), parameter :: F01 = "('"//LogPrefix//"',A,T45,' = ', f18.13)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    write(6,F0A)  'Case description:           ', trim(CCSMInit%case_desc)
    write(6,F0A)  'Case name:                  ', trim(CCSMInit%case_name)
    write(6,F0A)  'Start type:                 ', trim(CCSMInit%start_type)
    write(6,F0A)  'Archive directory:          ', trim(CCSMInit%archive_dir)
    write(6,F0A)  'Log output ending name:     ', trim(CCSMInit%logFilePostFix)
    write(6,F0A)  'Output path root directory: ', trim(CCSMInit%outPathRoot)
    write(6,F0L)  'atm idealized-physics mode: ', CCSMInit%atm_ideal_phys
    write(6,F0L)  'atm adiabatic mode:         ', CCSMInit%atm_adiabatic
    write(6,F0L)  'aqua_planet mode:           ', CCSMInit%aqua_planet
    write(6,F0L)  'single_column:              ', CCSMInit%single_column
    write(6,F01)  'scmlat:                     ', CCSMInit%scmlat
    write(6,F01)  'scmlon:                     ', CCSMInit%scmlon
    write(6,F0A)  'Restart pointer file:       ', trim(CCSMInit%restart_pfile)
    write(6,F0I)  'Archive retention period:   ', CCSMInit%mss_irt
    write(6,F0A)  'Archive write password:     ', trim(CCSMInit%mss_wpass)
    ! --- Some values only printed out on restart ------
    if ( shr_inputInfo_initIsRestart( CCSMInit ) )then
       if ( shr_inputInfo_initIsBranch( CCSMInit ) ) &
          write(6,F0A)  'Restart file (full path):   ', trim(CCSMInit%restart_file)
       if ( len_trim(CCSMInit%restart_file_override) == 0 )then
          write(6,FAA)  'Do NOT override any values from restart_file by namelist'
       else
          write(6,F0A)  'Restart file override list: ', &
                    trim(CCSMInit%restart_file_override)
       end if
    end if

END SUBROUTINE shr_inputInfo_initPrint

!===============================================================================
!===============================================================================
! !IROUTINE: shr_inputInfo_initGetData -- Get values from CCSMInit object
!   
! !DESCRIPTION:
!   
!     Get values out of the CCSMInit object.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_inputInfo_initGetData( CCSMInit, case_name, case_desc,       &
                                      archive_dir, mss_irt, mss_wpass,      &
                                      start_type, restart_pfile,            &
                                      restart_file, perpetual_run,          &
                                      perpetual_ymd, aqua_planet,           &
                                      atm_ideal_phys, atm_adiabatic,        &
                                      brnch_retain_casename, single_column, &
                                      scmlat,scmlon, logFilePostFix,        &
                                      outPathRoot      )
    
   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_inputInfo_InitType),  intent(IN) :: CCSMInit    ! Input CCSM structure
   character(len=*), optional, intent(OUT) :: case_name     ! Short case identification
   character(len=*), optional, intent(OUT) :: case_desc     ! Long case description
   character(len=*), optional, intent(OUT) :: archive_dir   ! Archival directory
   integer,          optional, intent(OUT) :: mss_irt       ! MSS retention period
   character(len=*), optional, intent(OUT) :: mss_wpass     ! Mass Store write pass
   character(len=*), optional, intent(OUT) :: start_type    ! Start type
   character(len=*), optional, intent(OUT) :: restart_pfile ! Restart pointer file
   character(len=*), optional, intent(OUT) :: restart_file  ! Restart file pathname
   logical,          optional, intent(OUT) :: perpetual_run ! If this is perpetual
   integer,          optional, intent(OUT) :: perpetual_ymd ! If perpetual, date
   logical,          optional, intent(OUT) :: aqua_planet   ! aqua_planet mode
   logical,          optional, intent(OUT) :: atm_adiabatic ! atm adiabatic mode
   logical,          optional, intent(OUT) :: atm_ideal_phys! atm idealized-physics mode
   character(len=*), optional, intent(OUT) :: logFilePostFix! output log file postfix
   character(len=*), optional, intent(OUT) :: outPathRoot   ! output file root
   ! --- If branch and can retain the restart file casename ------
   logical,          optional, intent(OUT) :: brnch_retain_casename
   logical,          optional, intent(OUT) :: single_column
   real (SHR_KIND_R8),  optional, intent(OUT) :: scmlat
   real (SHR_KIND_R8),  optional, intent(OUT) :: scmlon

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(shr_inputInfo_initGetData) '
    integer, parameter :: aqua_perpetual_ymd = 321       ! Default perpetual date

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    if ( present(aqua_planet)   )   aqua_planet    = CCSMInit%aqua_planet
    if ( present(single_column)   ) single_column  = CCSMInit%single_column
    if ( present(scmlat)   )        scmlat         = CCSMInit%scmlat
    if ( present(scmlon)   )        scmlon         = CCSMInit%scmlon
    if ( present(atm_adiabatic)  )  atm_adiabatic  = CCSMInit%atm_adiabatic
    if ( present(atm_ideal_phys) )  atm_ideal_phys = CCSMInit%atm_ideal_phys
    if ( present(brnch_retain_casename) ) brnch_retain_casename = &
                                          CCSMInit%brnch_retain_casename
    if ( present(start_type)   )    start_type     = CCSMInit%start_type
    if ( present(restart_pfile) )   restart_pfile  = CCSMInit%restart_pfile
    if ( present(restart_file)  )   restart_file   = CCSMInit%restart_file
    if ( present(case_name)  )      case_name      = CCSMInit%case_name
    if ( present(case_desc)  )      case_desc      = CCSMInit%case_desc
    if ( present(archive_dir)  )    archive_dir    = CCSMInit%archive_dir
    if ( present(mss_irt)  )        mss_irt        = CCSMInit%mss_irt
    if ( present(mss_wpass)  )      mss_wpass      = CCSMInit%mss_wpass
    if ( present(logFilePostFix) ) logFilePostFix= CCSMInit%logFilePostFix
    if ( present(outPathRoot) )    outPathRoot   = CCSMInit%outPathRoot
    
    ! --- Perpetual run mode and date (if aqua_planet mode set else false) ------
    if ( present(perpetual_run)  )then
       if ( .not. present(perpetual_ymd)  )then
          call shr_sys_abort( subname//':: asking for perpetual_run but not '//&
                              'perpetual_ymd' )
       end if
       if ( CCSMInit%aqua_planet )then
          perpetual_run  = .true.
          perpetual_ymd  = aqua_perpetual_ymd
       else
          perpetual_run  = .false.
          perpetual_ymd  = 0
       end if
     end if

END SUBROUTINE shr_inputInfo_initGetData

!===============================================================================
!===============================================================================
! !IROUTINE: shr_inputInfo_initPutData -- Change values of input CCSMInit object
!   
! !DESCRIPTION:
!   
!     Change the values of the input CCSMInit object.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_inputInfo_initPutData( CCSMInit, case_name, case_desc,     &
                                      start_type, mss_irt, atm_adiabatic, &
                                      atm_ideal_phys, aqua_planet,        &
                                      archive_dir, restart_file,          &
                                      restart_pfile, mss_wpass,           &
                                      restart_file_override,              &
                                      brnch_retain_casename,single_column,&
                                      scmlat,scmlon, logFilePostFix,      &
                                      outPathRoot         )

   implicit none
    
! !INPUT/OUTPUT PARAMETERS:

    type(shr_inputInfo_InitType), intent(INOUT)  :: CCSMInit  ! Output CCSM structure
    character(len=*), intent(IN), optional :: case_name       ! Short case identification
    character(len=*), intent(IN), optional :: case_desc       ! Long case description
    character(len=*), intent(IN), optional :: start_type      ! Type of startup
    character(len=*), intent(IN), optional :: restart_file    ! Restart file pathname
    character(len=*), intent(IN), optional :: restart_pfile   ! Restart pointer file
    character(len=*), intent(IN), optional :: archive_dir     ! Archive directory name
    character(len=*), intent(IN), optional :: logFilePostFix  ! Output log file postfix
    character(len=*), intent(IN), optional :: outPathRoot     ! Output file root
    ! -- Items to override on restart_file --- 
    character(len=*), intent(IN), optional :: restart_file_override 
    character(len=*), intent(IN), optional :: mss_wpass       ! MSS write password
    integer,          intent(IN), optional :: mss_irt         ! MSS retention period
    logical,          intent(IN), optional :: atm_adiabatic   ! no surf & atm adiabatic
    logical,          intent(IN), optional :: atm_ideal_phys  ! no surf & atm ideal-phyps
    logical,          intent(IN), optional :: aqua_planet     ! Pure ocean surface
    logical,          intent(IN), optional :: single_column   ! single column mode
    real (SHR_KIND_R8),intent(IN), optional :: scmlat          ! single column mode latitude
    real (SHR_KIND_R8),intent(IN), optional :: scmlon          ! single column mode longitude
    ! --- If retain previous casename or not ------
    logical,          intent(IN), optional :: brnch_retain_casename

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(shr_inputInfo_initPutData) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    if ( present(case_name)        ) CCSMInit%case_name      = case_name
    if ( present(case_desc)        ) CCSMInit%case_desc      = case_desc
    if ( present(start_type)       ) CCSMInit%start_type     = start_type
    if ( present(mss_irt)          ) CCSMInit%mss_irt        = mss_irt
    if ( present(mss_wpass)        ) CCSMInit%mss_wpass      = mss_wpass
    if ( present(atm_adiabatic)    ) CCSMInit%atm_adiabatic  = atm_adiabatic
    if ( present(atm_ideal_phys)   ) CCSMInit%atm_ideal_phys = atm_ideal_phys
    if ( present(aqua_planet)      ) CCSMInit%aqua_planet    = aqua_planet
    if ( present(single_column)    ) CCSMInit%single_column  = single_column
    if ( present(scmlat)           ) CCSMInit%scmlat         = scmlat
    if ( present(scmlon)           ) CCSMInit%scmlon         = scmlon
    if ( present(restart_file)     ) CCSMInit%restart_file   = restart_file
    if ( present(restart_pfile)    ) CCSMInit%restart_pfile  = restart_pfile
    if ( present(archive_dir)      ) CCSMInit%archive_dir    = archive_dir
    if ( present(logFilePostFix)   ) CCSMInit%logFilePostFix = logFilePostFix 
    if ( present(outPathRoot)      ) CCSMInit%outPathRoot    = outPathRoot 
    if ( present(restart_file_override) ) &
          CCSMInit%restart_file_override = restart_file_override
    if ( present(brnch_retain_casename)) CCSMInit%brnch_retain_casename = &
                                         brnch_retain_casename

END SUBROUTINE shr_inputInfo_initPutData

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_inputInfo_initRead -- read in CCSM shared namelist
!
! !DESCRIPTION:
!
!     Read in input from ccsm_inparm namelist, output ccsm derived type for
!     miscillaneous info.
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_inputInfo_initRead( NLFilename, LogPrint, mpicom, masterTask, &
                                   CCSMInitOut )

! !USES:

   USE shr_file_mod,   ONLY : shr_file_getUnit, shr_file_freeUnit
   USE shr_string_mod, ONLY : shr_string_toUpper, shr_string_listAppend
   USE shr_mpi_mod,    ONLY : shr_mpi_bcast

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(len=*),            intent(IN) :: NLFilename      ! Name-list filename
   logical,                     intent(IN) :: LogPrint        ! If print out to log file
   integer, optional,           intent(IN) :: mpicom          ! MPI communicator
   logical, optional,           intent(IN) :: MasterTask      ! If MPI master task
   type(shr_inputInfo_InitType), intent(INOUT) :: CCSMInitOut ! Output CCSM structure

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(shr_inputInfo_initRead) '
    logical                    :: MasterTask2        ! If MPI master task
    character(len=SHR_KIND_CL) :: case_desc          ! Case long description
    character(SHR_KIND_CS)     :: case_name          ! Case short name
    character(SHR_KIND_CL)     :: start_type         ! Startup-type: startup, continue, branch
    character(SHR_KIND_CL)     :: archive_dir        ! Archive directory
    character(SHR_KIND_CL)     :: restart_pfile      ! Restart pointer filename
    character(SHR_KIND_CL)     :: restart_file       ! Restart filename
    ! ----- Namelist items to override from the restart-file ------
    character(SHR_KIND_CL)     :: restart_file_override     
    logical                    :: aqua_planet        ! Aqua-planet mode (surface is all ocean)
    logical                    :: single_column      ! single column mode
    real (SHR_KIND_R8)         :: scmlat             ! single column mode latitude
    real (SHR_KIND_R8)         :: scmlon             ! single column mode longitude
    logical                    :: atm_adiabatic      ! Atmosphere adiabatic physics mode
    logical                    :: atm_ideal_phys     ! Atmosphere idealized physics mode
    logical                    :: brnch_retain_casename ! If retain casename for branch
    character(len=SHR_KIND_CS) :: logFilePostFix     ! postfix for output log files
    character(len=SHR_KIND_CS) :: outPathRoot        ! root output files
    integer(SHR_KIND_IN)       :: mss_irt            ! Archive retention period
    character(len=SHR_KIND_CL) :: mss_wpass          ! Archive write password
    character(len=SHR_KIND_CL) :: home_dir           ! Home directory
    character(len=SHR_KIND_CL) :: logname            ! user name
    character(len=SHR_KIND_CL) :: upCaseLogName      ! user name in upper-case
    integer                    :: rcode              ! Return code
    integer                    :: ierr               ! I/O error code
    integer                    :: unitn              ! Namelist unit number to read

    namelist /ccsm_inparm/  &
         case_desc, case_name, start_type, archive_dir, &
         atm_adiabatic, atm_ideal_phys, aqua_planet,    &
         mss_irt, mss_wpass, brnch_retain_casename,     &
         restart_pfile, restart_file,                   &
         restart_file_override,single_column, scmlat,   &
         scmlon, logFilePostFix, outPathRoot 

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    if ( present(MasterTask) )then
       MasterTask2 = MasterTask
    else
       MasterTask2 = .true.
    end if
    if (       present(MasterTask) .and. .not. present(mpicom) )then
       call shr_sys_abort( subname//':: if MasterTask input, so must mpicom' )
    endif
    if ( .not. present(MasterTask) .and.      present(mpicom) )then
       call shr_sys_abort( subname//':: if mpicom input, so must MasterTask' )
    endif
    !---------------------------------------------------------------------------
    ! Set defaults based on current value of input object
    !---------------------------------------------------------------------------
    case_desc             = CCSMInitOut%case_desc
    case_name             = CCSMInitOut%case_name
    start_type            = CCSMInitOut%start_type
    archive_dir           = CCSMInitOut%archive_dir
    atm_ideal_phys        = CCSMInitOut%atm_ideal_phys
    atm_adiabatic         = CCSMInitOut%atm_adiabatic
    aqua_planet           = CCSMInitOut%aqua_planet
    single_column         = CCSMInitOut%single_column
    scmlat                = CCSMInitOut%scmlat
    scmlon                = CCSMInitOut%scmlon
    mss_irt               = CCSMInitOut%mss_irt
    mss_wpass             = CCSMInitOut%mss_wpass
    restart_pfile         = CCSMInitOut%restart_pfile
    restart_file          = CCSMInitOut%restart_file
    restart_file_override = CCSMInitOut%restart_file_override
    brnch_retain_casename = CCSMInitOut%brnch_retain_casename
    logFilePostFix        = CCSMInitOut%logFilePostFix
    outPathRoot           = CCSMInitOut%outPathRoot
    !---------------------------------------------------------------------------
    ! Read in namelist
    !---------------------------------------------------------------------------
    if ( MasterTask2 )then
       unitn = shr_file_getUnit()
       write(6,FAA) 'Read in ccsm_inparm namelist from: '//trim(NLFilename)
       open( unitn, file=trim(NLFilename), status='old' )
       ierr = 1
       do while( ierr /= 0 )
          read(unitn,nml=ccsm_inparm,iostat=ierr)
          if (ierr < 0) then
             call shr_sys_abort( subname//':: namelist read returns an'// &
                                 ' end of file or end of record condition' )
          end if
       end do
       close(unitn)
       call shr_file_freeUnit( unitn )
    end if
    !
    ! If mss_irt too high truncate to reasonable value
    !
    if ( mss_irt > 4096 ) mss_irt = 4096
    !---------------------------------------------------------------------------
    ! If MPI broadcast values to all processors
    !---------------------------------------------------------------------------
    if ( present(mpicom) )then
       call shr_mpi_bcast( case_desc,             MPICom )
       call shr_mpi_bcast( case_name,             MPICom )
       call shr_mpi_bcast( start_type,            MPICom )
       call shr_mpi_bcast( archive_dir,           MPICom )
       call shr_mpi_bcast( atm_ideal_phys,        MPICom )
       call shr_mpi_bcast( atm_adiabatic,         MPICom )
       call shr_mpi_bcast( aqua_planet,           MPICom )
       call shr_mpi_bcast( single_column,         MPICom )
       call shr_mpi_bcast( scmlat,                MPICom )
       call shr_mpi_bcast( scmlon,                MPICom )
       call shr_mpi_bcast( mss_irt,               MPICom )
       call shr_mpi_bcast( mss_wpass,             MPICom )
       call shr_mpi_bcast( restart_pfile,         MPICom )
       call shr_mpi_bcast( restart_file,          MPICom )
       call shr_mpi_bcast( restart_file_override, MPICom )
       call shr_mpi_bcast( brnch_retain_casename, MPICom )
       call shr_mpi_bcast( logFilePostFix,        MPICom )
       call shr_mpi_bcast( outPathRoot,           MPICom )
    end if
    !---------------------------------------------------------------------------
    ! Put into ccsm init derived type
    !---------------------------------------------------------------------------
    CCSMInitOut%case_desc             = case_desc
    CCSMInitOut%case_name             = case_name
    CCSMInitOut%start_type            = start_type
    CCSMInitOut%archive_dir           = archive_dir
    CCSMInitOut%atm_ideal_phys        = atm_ideal_phys
    CCSMInitOut%atm_adiabatic         = atm_adiabatic
    CCSMInitOut%aqua_planet           = aqua_planet
    CCSMInitOut%single_column         = single_column
    CCSMInitOut%scmlat                = scmlat
    CCSMInitOut%scmlon                = scmlon
    CCSMInitOut%mss_irt               = mss_irt
    CCSMInitOut%mss_wpass             = mss_wpass
    CCSMInitOut%restart_pfile         = restart_pfile
    CCSMInitOut%restart_file          = restart_file
    CCSMInitOut%restart_file_override = restart_file_override
    CCSMInitOut%brnch_retain_casename = brnch_retain_casename
    CCSMInitOut%logFilePostFix        = logFilePostFix
    CCSMInitOut%outPathRoot           = outPathRoot
    !---------------------------------------------------------------------------
    ! Set restart pointer file name if wasn't set on namelist
    !---------------------------------------------------------------------------
    call shr_sys_getenv ('HOME',home_dir,rcode)
    if (rcode /= 0) then
       call shr_sys_abort( subname//': Cannot find HOME environment variable')
    end if
    if ( len_trim(CCSMInitOut%restart_pfile) == 0 )then
       CCSMInitOut%restart_pfile = trim(home_dir) //'/drv.'//  &
                                   trim(CCSMInitOut%case_name) //'.RPointer'
    end if
    !---------------------------------------------------------------------------
    ! Set archive directory name
    !---------------------------------------------------------------------------
    if ( len_trim(CCSMInitOut%archive_dir) == 0 )then
      logname = ' '
      call shr_sys_getenv ('LOGNAME',logname,rcode)
      if (rcode /= 0) then
         call shr_sys_abort( subname//': Cannot find LOGNAME environment '//&
                             'variable')
      end if
      upCaseLogName = shr_string_toUpper(logname)
      CCSMInitOut%archive_dir   = '/'//trim(upCaseLogName)//'/csm/'// &
                                       trim(CCSMInitOut%case_name)//'/drv/'
      ! This behavior should go away at Stage-2, and all archiving yanked out
      if ( CCSMInitOut%mss_irt > 0 )then
         CCSMInitOut%archive_dir   = "mss:"//CCSMInitOut%archive_dir
      else
         CCSMInitOut%archive_dir   = "null:"//CCSMInitOut%archive_dir
      end if
    end if
    !---------------------------------------------------------------------------
    ! Append values to override list depending on type
    ! (For branch, override, case_name, for restart don't)
    !---------------------------------------------------------------------------
    if ( shr_inputInfo_initIsBranch( CCSMInitOut ) .and.  &
         (.not. CCSMInitOut%brnch_retain_casename) )then
       call shr_string_listAppend( CCSMInitOut%restart_file_override, &
                                   "case_name" )
    end if
    ! --- Error check the input values and print them out ------
    call shr_inputInfo_initCheck( CCSMInitOut )
    if ( LogPrint) call shr_inputInfo_initPrint( CCSMInitOut )

END SUBROUTINE shr_inputInfo_initRead

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_inputInfo_initIsStartup -- Return true if this is a startup case
!
! !DESCRIPTION:
!
! Check if this is a startup case.
!
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION shr_inputInfo_initIsStartup( CCSMInit )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_inputInfo_InitType), intent(IN) :: CCSMInit    ! Input CCSM structure

!EOP

    !----- local -----

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  ! --- Return if startup or not -----
   if ( trim(CCSMInit%start_type) == shr_inputInfo_start_type_start )then
      shr_inputInfo_initIsStartup = .true.
   else
      shr_inputInfo_initIsStartup = .false.
   end if

END FUNCTION shr_inputInfo_initIsStartup

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_inputInfo_initIsRestart -- Return true if this is continue or branch
!
! !DESCRIPTION:
!
! Check if this is a continue or "branch" case"
!
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION shr_inputInfo_initIsRestart( CCSMInit )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_inputInfo_InitType), intent(IN) :: CCSMInit    ! Input CCSM structure

!EOP

    !----- local -----

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  ! --- Return if restart or not (negation of startup case)
  shr_inputInfo_initIsRestart = .not. shr_inputInfo_initIsStartup( CCSMInit )

END FUNCTION shr_inputInfo_initIsRestart

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_inputInfo_initIsContinue -- Return true if this is a continue case
!
! !DESCRIPTION:
!
! Check if this is a continue case.
!
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION shr_inputInfo_initIsContinue( CCSMInit )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_inputInfo_InitType), intent(IN) :: CCSMInit    ! Input CCSM structure

!EOP

    !----- local -----

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  ! --- Return if continue or not ------
   if ( trim(CCSMInit%start_type) == shr_inputInfo_start_type_cont )then
      shr_inputInfo_initIsContinue  = .true.
   else
      shr_inputInfo_initIsContinue  = .false.
   end if

END FUNCTION shr_inputInfo_initIsContinue

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_inputInfo_initIsBranch -- Return true if this is a branch case
!
! !DESCRIPTION:
!
! Check if this is a branch case.
!
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION shr_inputInfo_initIsBranch( CCSMInit )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_inputInfo_InitType), intent(IN) :: CCSMInit    ! Input CCSM structure

!EOP

    !----- local -----

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  ! --- Return if branch or not -----
   if ( trim(CCSMInit%start_type) == shr_inputInfo_start_type_brnch )then
      shr_inputInfo_initIsBranch = .true.
   else
      shr_inputInfo_initIsBranch = .false.
   end if

END FUNCTION shr_inputInfo_initIsBranch

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_inputInfo_initIsSame - Return true if input objects same
!
! !DESCRIPTION:
!
! Check if the two input init objects have the same values.
!
! !INTERFACE: ------------------------------------------------------------------
logical FUNCTION shr_inputInfo_initIsSame( CCSMInit, CCSMInit2 )

  implicit none

! !INPUT/OUTPUT PARAMETERS:
  type(shr_inputInfo_InitType), intent(IN) :: CCSMInit    ! Input CCSM structure
  type(shr_inputInfo_InitType), intent(IN) :: CCSMInit2   ! Second Input CCSM structure
!EOP
    !----- local -----

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    shr_inputInfo_initIsSame = .false.
    if ( CCSMInit%start_type            /= CCSMInit2%start_type            ) return
    if ( CCSMInit%case_desc             /= CCSMInit2%case_desc             ) return
    if ( CCSMInit%case_name             /= CCSMInit2%case_name             ) return
    if ( CCSMInit%atm_adiabatic     .neqv. CCSMInit2%atm_adiabatic         ) return
    if ( CCSMInit%atm_ideal_phys    .neqv. CCSMInit2%atm_ideal_phys        ) return
    if ( CCSMInit%aqua_planet       .neqv. CCSMInit2%aqua_planet           ) return
    if ( CCSMInit%mss_irt               /= CCSMInit2%mss_irt               ) return
    if ( CCSMInit%mss_wpass             /= CCSMInit2%mss_wpass             ) return
    if ( CCSMInit%restart_pfile         /= CCSMInit2%restart_pfile         ) return
    if ( CCSMInit%restart_file          /= CCSMInit2%restart_file          ) return
    if ( CCSMInit%archive_dir           /= CCSMInit2%archive_dir           ) return
    if ( CCSMInit%restart_file_override /= CCSMInit2%restart_file_override ) return
    if ( CCSMInit%logFilePostFix        /= CCSMInit2%logFilePostFix        ) return
    if ( CCSMInit%outPathRoot           /= CCSMInit2%outPathRoot           ) return
    if ( CCSMInit%brnch_retain_casename .neqv. CCSMInit2%brnch_retain_casename ) return
    shr_inputInfo_initIsSame = .true.
    return
END FUNCTION shr_inputInfo_initIsSame
!===============================================================================

!BOP ===========================================================================
!
! !IROUTINE: shr_inputInfo_initRestRead -- Read in NetCDF restart file
!
! !DESCRIPTION:
!
! Read in restart file information from netCDF input restart file
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_inputInfo_initRestRead( restart_file, MPICom, MasterTask, &
                                          CCSMInitOut )

!USES:

  use netcdf
  use shr_ncio_mod,   only: shr_ncio_descripRead, shr_ncio_descripName, &
                            shr_ncio_descripGetString,                  &
                            shr_ncio_descripGetRealR8,                   &
                            shr_ncio_descripGetInteger,                 &
                            shr_ncio_descripGetLogical
  use shr_string_mod, only: shr_string_listGetIndexF

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(len=SHR_KIND_CL),   intent(IN) :: restart_file   ! Restart local filename
  integer, optional,            intent(IN) :: MPICom         ! MPI Communicator
  logical, optional,            intent(IN) :: MasterTask     ! If MPI master task
  type(shr_inputInfo_InitType), intent(INOUT) :: CCSMInitOut ! Output CCSM structure

!EOP

   !----- local -----
  character(len=*), parameter :: subname = "(shr_inputInfo_initRestRead) "
  logical :: MasterTask2                    ! If MPI master task
  integer :: i                              ! Index
  integer :: list                           ! Index in override list
  integer :: ncId                           ! NetCDF file handle
  logical :: exists                         ! If file exists or not

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------
  if ( present(MasterTask) )then
     MasterTask2 = MasterTask
  else
     MasterTask2 = .true.
  end if
  if (       present(MasterTask) .and. .not. present(mpicom) )then
     call shr_sys_abort( subname//':: if MasterTask input, so must mpicom' )
  endif
  if ( .not. present(MasterTask) .and.      present(mpicom) )then
     call shr_sys_abort( subname//':: if mpicom input, so must MasterTask' )
  endif

  if ( .not. shr_inputInfo_initIsRestart( CCSMInitOut ) ) return
  if ( shr_inputInfo_initIsBranch( CCSMInitOut ) .and. len_trim(restart_file) == 0 &
  )then
     call shr_sys_abort( subname//': on a branch run, restart_file MUST be '// &
                         'set on namelist' )
  end if
  ! --- Read the restart file -----
  call shr_ncio_open( restart_file, MasterTask2, FileType=prefix//"restart_file", &
                      ncId=ncId, exists=exists )
  if ( present(MPICom) )then
     call shr_ncio_descripRead( ncId, NVars, prefix, MPICom, MasterTask2, &
                         var=CCSMInitOut%var )
  else
     call shr_ncio_descripRead( ncId, NVars, prefix, var=CCSMInitOut%var )
  end if
  call shr_ncio_close( ncId, MasterTask2, type=prefix//"restart_file", &
                       NCFileName=restart_file )
  !-----------------------------------------------------------------------------
  ! Put input read into derived type
  ! Check if item is on over-ride list (if list_index != 0)
  ! If it is not over-ride list -- use the value from the restart file
  !-----------------------------------------------------------------------------
  do i = 1, NVars
     if (      trim(shr_ncio_descripName(CCSMInitOut%var(i))) == "case_desc" )then
        list = shr_string_listGetIndexF( CCSMInitOut%restart_file_override, &
                         trim(shr_ncio_descripName(CCSMInitOut%var(i)) ) )
        if ( list == 0 ) CCSMInitOut%case_desc = &
                         shr_ncio_descripGetString(CCSMInitOut%var(i))
     else if ( trim(shr_ncio_descripName(CCSMInitOut%var(i))) == "case_name" )then
        list = shr_string_listGetIndexF( CCSMInitOut%restart_file_override, &
                                trim(shr_ncio_descripName(CCSMInitOut%var(i)) ) )
        if ( shr_inputInfo_initIsBranch( CCSMInitOut ) .and. &
             (.not. CCSMInitOut%brnch_retain_casename) )then
           if ( trim(CCSMInitOut%case_name) == &
                 trim(shr_ncio_descripGetString(CCSMInitOut%var(i))) )then
              call shr_sys_abort( subname//': Branch and case_name same as '// &
                          'on restart_file  and brnch_retain_casename not set' )
           end if
        end if
        if ( list == 0 ) CCSMInitOut%case_name = &
                      shr_ncio_descripGetString(CCSMInitOut%var(i))
     else if ( trim(shr_ncio_descripName(CCSMInitOut%var(i))) == "mss_wpass" )then
        list = shr_string_listGetIndexF( CCSMInitOut%restart_file_override, &
                            trim(shr_ncio_descripName(CCSMInitOut%var(i)) ) )
        if ( list == 0 ) CCSMInitOut%mss_wpass = &
                         shr_ncio_descripGetString(CCSMInitOut%var(i))
     else if ( trim(shr_ncio_descripName(CCSMInitOut%var(i))) == "atm_adiabatic" )then
        list = shr_string_listGetIndexF( CCSMInitOut%restart_file_override, &
                           trim(shr_ncio_descripName(CCSMInitOut%var(i)) ) )
        if ( list == 0 ) CCSMInitOut%atm_adiabatic = &
                                shr_ncio_descripGetLogical(CCSMInitOut%var(i))
     else if ( trim(shr_ncio_descripName(CCSMInitOut%var(i))) == "atm_ideal_phys" )then
        list = shr_string_listGetIndexF( CCSMInitOut%restart_file_override, &
                               trim(shr_ncio_descripName(CCSMInitOut%var(i)) ) )
        if ( list == 0 ) CCSMInitOut%atm_ideal_phys = &
                         shr_ncio_descripGetLogical(CCSMInitOut%var(i))
     else if ( trim(shr_ncio_descripName(CCSMInitOut%var(i))) == "aqua_planet" )then
        list = shr_string_listGetIndexF( CCSMInitOut%restart_file_override, &
                            trim(shr_ncio_descripName(CCSMInitOut%var(i)) ) )
        if ( list == 0 ) CCSMInitOut%aqua_planet  = &
                         shr_ncio_descripGetLogical(CCSMInitOut%var(i))
     else if ( trim(shr_ncio_descripName(CCSMInitOut%var(i))) == "single_column" )then
        list = shr_string_listGetIndexF( CCSMInitOut%restart_file_override, &
                            trim(shr_ncio_descripName(CCSMInitOut%var(i)) ) )
        if ( list == 0 ) CCSMInitOut%single_column  = &
                         shr_ncio_descripGetLogical(CCSMInitOut%var(i))
     else if ( trim(shr_ncio_descripName(CCSMInitOut%var(i))) == "scmlat" )then
        list = shr_string_listGetIndexF( CCSMInitOut%restart_file_override, &
                            trim(shr_ncio_descripName(CCSMInitOut%var(i)) ) )
        if ( list == 0 ) CCSMInitOut%scmlat = &
                         shr_ncio_descripGetRealR8(CCSMInitOut%var(i))
     else if ( trim(shr_ncio_descripName(CCSMInitOut%var(i))) == "scmlon" )then
        list = shr_string_listGetIndexF( CCSMInitOut%restart_file_override, &
                           trim(shr_ncio_descripName(CCSMInitOut%var(i)) ) )
        if ( list == 0 ) CCSMInitOut%scmlon = &
                         shr_ncio_descripGetRealR8(CCSMInitOut%var(i))
     else if ( trim(shr_ncio_descripName(CCSMInitOut%var(i))) == "mss_irt" )then
        list = shr_string_listGetIndexF( CCSMInitOut%restart_file_override, &
                             trim(shr_ncio_descripName(CCSMInitOut%var(i)) ) )
        if ( list == 0 ) CCSMInitOut%mss_irt = &
                         shr_ncio_descripGetInteger(CCSMInitOut%var(i))
     else
        call shr_sys_abort( subname//': unrecognized variable: '// &
                            trim(shr_ncio_descripName(CCSMInitOut%var(i))) )
     end if
  end do
  ! --- Check to make sure the structure is valid ------
  call shr_inputInfo_initCheck( CCSMInitOut )

END SUBROUTINE shr_inputInfo_initRestRead  

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_inputInfo_initRestWrite -- Write out the NetCDF restart file
!
! !DESCRIPTION:
!
! Write out restart file information to netCDF restart file
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_inputInfo_initRestWrite( restart_file, MPICom, MasterTask, &
                                           CCSMInit )
! !USES:
  use netcdf
  use shr_ncio_mod, only: shr_ncio_descripWrite, shr_ncio_descripPutData, &
                          shr_ncio_descripName

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(len=SHR_KIND_CL),   intent(IN) :: restart_file ! Restart local filename
  integer, optional,            intent(IN) :: MPICom       ! MPI Communicator
  logical, optional,            intent(IN) :: MasterTask   ! If MPI master task
  type(shr_inputInfo_InitType), intent(INOUT) :: CCSMInit  ! Output CCSM structure

!EOP

   !----- local -----
  character(len=*), parameter :: subname = "(shr_inputInfo_initRestWrite) "
  logical :: MasterTask2   ! If MPI master task
  logical :: exists        ! If file exists already or not
  integer :: ncId          ! NetCDF file handle
  integer :: i             ! Index

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  if ( present(MasterTask) )then
     MasterTask2 = MasterTask
  else
     MasterTask2 = .true.
  end if
  if (       present(MasterTask) .and. .not. present(mpicom) )then
     call shr_sys_abort( subname//':: if MasterTask input, so must mpicom' )
  endif
  if ( .not. present(MasterTask) .and.      present(mpicom) )then
     call shr_sys_abort( subname//':: if mpicom input, so must MasterTask' )
  endif
  !-----------------------------------------------------------------------------
  ! Loop over variables, set the data for each structure element
  !-----------------------------------------------------------------------------
  do i = 1, NVars
      if (      trim(shr_ncio_descripName( CCSMInit%var(i) ) ) == "atm_adiabatic" )then
         call shr_ncio_descripPutData( CCSMInit%var(i), "atm_adiabatic", &
                                   LogicalData=CCSMInit%atm_adiabatic )
      else if ( trim(shr_ncio_descripName( CCSMInit%var(i) ) ) == "atm_ideal_phys" &
      )then
         call shr_ncio_descripPutData( CCSMInit%var(i), "atm_ideal_phys", &
                                   LogicalData=CCSMInit%atm_ideal_phys )
      else if ( trim(shr_ncio_descripName( CCSMInit%var(i) ) ) == "aqua_planet" )then
         call shr_ncio_descripPutData( CCSMInit%var(i), "aqua_planet", &
                                   LogicalData=CCSMInit%aqua_planet )
      else if ( trim(shr_ncio_descripName( CCSMInit%var(i) ) ) == "single_column" )then
         call shr_ncio_descripPutData( CCSMInit%var(i), "single_column", &
                                   LogicalData=CCSMInit%single_column )
      else if ( trim(shr_ncio_descripName( CCSMInit%var(i) ) ) == "scmlat" )then
         call shr_ncio_descripPutData( CCSMInit%var(i), "scmlat", &
                                   RealR8Data=CCSMInit%scmlat )
      else if ( trim(shr_ncio_descripName( CCSMInit%var(i) ) ) == "scmlon" )then
         call shr_ncio_descripPutData( CCSMInit%var(i), "scmlon", &
                                   RealR8Data=CCSMInit%scmlon )
      else if ( trim(shr_ncio_descripName( CCSMInit%var(i) ) ) == "mss_irt" )then
         call shr_ncio_descripPutData( CCSMInit%var(i), "mss_irt",  &
                                   IntegerData=CCSMInit%mss_irt )
      else if ( trim(shr_ncio_descripName( CCSMInit%var(i) ) ) == "case_desc" )then
         call shr_ncio_descripPutData( CCSMInit%var(i), "case_desc", &
                                   StringData=CCSMInit%case_desc )
      else if ( trim(shr_ncio_descripName( CCSMInit%var(i) ) ) == "case_name" )then
         call shr_ncio_descripPutData( CCSMInit%var(i), "case_name", &
                                   StringData=CCSMInit%case_name )
      else if ( trim(shr_ncio_descripName( CCSMInit%var(i) ) ) == "mss_wpass" )then
         call shr_ncio_descripPutData( CCSMInit%var(i), "mss_wpass", &
                                   StringData=CCSMInit%mss_wpass )
      else
         call shr_sys_abort( subname//': unknown variable name: '// &
                             trim(shr_ncio_descripName( CCSMInit%var(i) ) ) )
      end if
  end do
  ! --- Actually write the data to the file -------
  call shr_ncio_open( NCFileName=restart_file, MasterTask=masterTask2,  &
                      FileType=trim(prefix)//"restart_file",            &
                      ncId=ncId, exists=exists, writing=.true. )
  if ( present(MPICom) )then
     call shr_ncio_descripWrite( ncId, NVars, prefix, mpicom, MasterTask2, &
                                 exists=exists, var=CCSMInit%var )
  else
     call shr_ncio_descripWrite( ncId, NVars, prefix, exists=exists, &
                                 var=CCSMInit%var )
  end if
  call shr_ncio_close( ncId, MasterTask2, type=prefix//"restart_file",  &
                       NCFileName=restart_file )

END SUBROUTINE shr_inputInfo_initRestWrite  

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_inputInfo_initIsModelToRun -- Return true if given model should run.
!
! !DESCRIPTION:
!
! Return true if should run model component. Query the InputInfo data structure for
! options that determine if the given type of model component should run or not.
!
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION shr_inputInfo_initIsModelToRun( CCSMInit, model )

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(shr_inputInfo_InitType), intent(IN) :: CCSMInit ! Output CCSM structure
  character(len=*),             intent(IN) :: model    ! Which model component

!EOP

  !----- local -----
  character(len=*), parameter :: subname = "(shr_inputInfo_initIsModelToRun) "
  integer :: iModel                      ! Model index
  integer, parameter :: nModels = 4      ! Number of model components
  integer, parameter :: LND = 1          ! Land model type
  integer, parameter :: OCN = 2          ! Ocean model type
  integer, parameter :: ICE = 3          ! Sea-ice model type
  integer, parameter :: ATM = 4          ! Atmosphere model type
  logical :: aqua_planet(nModels)        ! Run table for aqua_planet
  logical :: single_column(nModels)      ! Run table for single column mode
  logical :: atm_adiabatic(nModels)      ! Run table for atm_adiabatic
  logical :: atm_ideal_phys(nModels)     ! Run table for atm_ideal_phys

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  if (      trim(model) == "atm" )then
     iModel = ATM
  else if ( trim(model) == "ice" )then
     iModel = ICE
  else if ( trim(model) == "lnd" )then
     iModel = LND
  else if ( trim(model) == "ocn" )then
     iModel = OCN
  else
     call shr_sys_abort( subname//": error: bad model type" )
  end if
  aqua_planet(LND)    = .false.
  aqua_planet(OCN)    = .true.
  aqua_planet(ICE)    = .false.
  aqua_planet(ATM)    = .true.
  atm_adiabatic(LND)  = .false.
  atm_adiabatic(OCN)  = .false.
  atm_adiabatic(ICE)  = .false.
  atm_adiabatic(ATM)  = .true.
  atm_ideal_phys(LND) = .false.
  atm_ideal_phys(OCN) = .false.
  atm_ideal_phys(ICE) = .false.
  atm_ideal_phys(ATM) = .true.
  single_column(LND)    = .true.
  single_column(OCN)    = .true.
  single_column(ICE)    = .true.
  single_column(ATM)    = .true.

  if (      CCSMInit%aqua_planet    )then
     shr_inputInfo_initIsModelToRun = aqua_planet(iModel)
  else if ( CCSMInit%single_column  )then
     shr_inputInfo_initIsModelToRun = single_column(iModel)
  else if ( CCSMInit%atm_adiabatic  )then
     shr_inputInfo_initIsModelToRun = atm_adiabatic(iModel)
  else if ( CCSMInit%atm_ideal_phys )then
     shr_inputInfo_initIsModelToRun = atm_ideal_phys(iModel)
  else
     shr_inputInfo_initIsModelToRun = .true.
  end if

END FUNCTION shr_inputInfo_initIsModelToRun

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_inputInfo_initRPointerRead  -- Read in the restart pointer file
!
! !DESCRIPTION:
!
! Read in restart pointer file, also get the local restart filename
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_inputInfo_initRPointerRead( MPICom, MasterTask, CCSMInit, &
                                           restart_file )
! !USES:

  USE shr_file_mod,  ONLY : shr_file_getUnit, shr_file_freeUnit, shr_file_get
  USE shr_mpi_mod,    ONLY : shr_mpi_bcast

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer, optional,            intent(IN)    :: MPICom       ! MPI communicator
  logical, optional,            intent(IN)    :: MasterTask   ! If MPI master task
  type(shr_inputInfo_InitType), intent(INOUT) :: CCSMInit     ! Output CCSM structure
  character(len=SHR_KIND_CL),   intent(OUT)   :: restart_file ! Local restart filename

!EOP

  !----- local -----
  character(len=*), parameter :: subname = "(shr_inputInfo_initRPointerRead) "
  integer :: rUnit                                ! Fortran unit to open
  logical :: MasterTask2                          ! If MPI master task
  character(len=SHR_KIND_CL) :: rest_file_archive ! Restart filename full archive path
  integer :: ierr                                 ! Return code
  integer :: i                                    ! Index

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  if ( present(MasterTask) )then
     MasterTask2 = MasterTask
  else
     MasterTask2 = .true.
  end if
  if (       present(MasterTask) .and. .not. present(mpicom) )then
     call shr_sys_abort( subname//':: if MasterTask input, so must mpicom' )
  endif
  if ( .not. present(MasterTask) .and.      present(mpicom) )then
     call shr_sys_abort( subname//':: if mpicom input, so must MasterTask' )
  endif
  ! --- Only continue if this is a continue or restart case
  if ( .not. shr_inputInfo_initIsRestart( CCSMInit ) )then
     call shr_sys_abort( subname//":: can only read RPointer file if "// &
                         "restart or continue" )
  end if

  ! --- Only read restart pointer file -- if this is a continue case -----
  if ( shr_inputInfo_initIsContinue( CCSMInit ) )then
     ! --- Read restart pointer file only on MasterTask -----
     if ( MasterTask2 )then
        rUnit = shr_file_getUnit()
        write(6,FAA) 'Read ccsm driver restart pointer file: '// &
                     trim(CCSMInit%restart_pfile)
        open(rUnit, file=CCSMInit%restart_pfile, form='FORMATTED', status='old')
        read(rUnit,'(a)') rest_file_archive
        close(rUnit)
        call shr_file_freeUnit( rUnit )
     end if
     ! --- If MPI broadcast archive filename to all tasks -------
     if ( present(MPICom) ) call shr_mpi_bcast( rest_file_archive, MPICom )
  else
  ! --- For branch case get restart_file archive name from input object ---
     rest_file_archive = CCSMInit%restart_file
  end if
  CCSMInit%restart_file = rest_file_archive
  ! --- Get restart file from archival system, and get local filename ----
  i = index( rest_file_archive, "/", back=.true. )
  if ( i > 0 )then
     restart_file = rest_file_archive(i+1:len_trim(rest_file_archive))
  else
     restart_file = rest_file_archive
  end if
  if ( MasterTask2 ) call shr_file_get( ierr, restart_file, rest_file_archive )

END SUBROUTINE shr_inputInfo_initRPointerRead

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_inputInfo_initRPointerWrite  -- Write out the restart pointer file
!
! !DESCRIPTION:
!
! Write out restart pointer file also get local restart filename
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_inputInfo_initRPointerWrite( ymd, tod, MPICom, MasterTask, &
                                            FileType, CCSMInit, restart_file )
! !USES:

  use shr_file_mod,   only: shr_file_getUnit, shr_file_freeUnit, shr_file_get

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,                      intent(IN)    :: ymd          ! Simulation date (YYYYMMDD)
  integer,                      intent(IN)    :: tod          ! Simulation time of day (s)
  integer,          optional,   intent(IN)    :: MPICom       ! MPI communicator
  logical,          optional,   intent(IN)    :: MasterTask   ! If MPI master task
  character(len=*), optional,   intent(IN)    :: FileType     ! File-type descriptor
  type(shr_inputInfo_InitType), intent(INOUT) :: CCSMInit     ! Output CCSM structure
  character(len=SHR_KIND_CL),   intent(OUT)   :: restart_file ! Restart local filename

!EOP

  !----- local -----
  character(len=*), parameter :: subname = "(shr_inputInfo_initRPointerWrite) "
  integer :: rUnit                                   ! Restart pointer file unit number`
  character(len=SHR_KIND_CL) :: restart_file_archive ! Restart filename full archive path
  character(len=*), parameter :: filename_format = &
                   "(a,'.',a,'.', i4.4,'-',i2.2,'-',i2.2,'-',i5.5,'.nc')"
  character(len=SHR_KIND_CS) :: type          ! String to describe type of restart
  logical :: MasterTask2                      ! If MPI master task
  integer :: ierr                             ! Error code
  integer :: yr                               ! Simulation year
  integer :: mon                              ! Simulation month
  integer :: day                              ! Simulation day of month

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  if ( present(MasterTask) )then
     MasterTask2 = MasterTask
  else
     MasterTask2 = .true.
  end if
  if ( present(FileType) )then
     type = trim(FileType)
  else
     type = "drv.r"
  end if
  
  ! --- Figure out restart filename -------
  yr   = ymd / 10000
  mon  = (ymd - yr*10000) / 100
  day  =  ymd - yr*10000 - mon*100
  write(restart_file,filename_format) trim(CCSMInit%case_name), trim(type), &
                                      yr, mon, day, tod
  restart_file_archive  = trim(CCSMInit%archive_dir)//trim(restart_file)
  CCSMInit%restart_file = restart_file_archive
  ! --- Only write file out if master-task -------
  if ( MasterTask2 )then
     write(6,FAA) 'Write ccsm driver restart pointer file: '// &
                   trim(CCSMInit%restart_pfile)
     rUnit = shr_file_getUnit()
     open(rUnit, file=CCSMInit%restart_pfile, form='FORMATTED', &
          status='UNKNOWN')
     write(rUnit,'(a)') trim(restart_file)
     close(rUnit)
     call shr_file_freeUnit( rUnit )
     call shr_file_get( ierr, restart_file, restart_file_archive )
  end if

END SUBROUTINE shr_inputInfo_initRPointerWrite

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_inputInfo_initRestArchive  -- Archive the restart file
!
! !DESCRIPTION:
!
! Archive restart file to the archival system.
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_inputInfo_initRestArchive( MasterTask, CCSMInit, async, remove, &
                                          restart_file )

! !USES:

  use shr_file_mod, only : shr_file_put

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  logical, intent(IN) :: MasterTask                         ! If MPI master task
  type(shr_inputInfo_InitType), intent(IN) :: CCSMInit      ! Input CCSM structure
  logical, intent(IN) :: async                              ! If archive asynchronously
  logical, intent(IN) :: remove                             ! If remove file afterwards
  character(len=SHR_KIND_CL), intent(INOUT) :: restart_file ! Restart local filename

!EOP

  !----- local -----
  character(len=*), parameter :: subname = '(shr_inputInfo_initRestArchive) '
  integer :: ierr        ! Error code

!-------------------------------------------------------------------------------
! Notes:
! This should go completely away at Stage-2. The archiving should go away by then,
! and this is hard-coded to assume MSS type archiving (hence triggered by mss_irt).
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Only archive from master-task, and only if retention period if greater than zero
  !-----------------------------------------------------------------------------
  if ( CCSMInit%mss_irt > 0 .and. MasterTask )then
     ! --- If have a write password -------
     if ( len_trim(CCSMInit%mss_wpass) > 0 )then
        call shr_file_put( ierr, restart_file, CCSMInit%restart_file,        &
                           passwd=CCSMInit%mss_wpass, rtpd=CCSMInit%mss_irt, &
                           async=async, remove=remove )
     ! --- If don't have a write password -------
     else
        call shr_file_put( ierr, restart_file, CCSMInit%restart_file,         &
                           rtpd=CCSMInit%mss_irt, async=async, remove=remove )
     end if
     ! --- Die if an error occurred in the archival process ------
     if ( ierr /= 0 ) &
        call shr_sys_abort( subname//": error on archival of restart file: "// &
                            trim(CCSMInit%restart_file) )
  end if
  ! --- Set restart_filename to null, so can't add more data to it afterwards... ------
  restart_file = ""
  
END SUBROUTINE shr_inputInfo_initRestArchive

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_inputInfo_initCheck  -- Check that input InputInfo derived type is valid
!
! !DESCRIPTION:
!
! Check that input CCSMInit object has reasonable values
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_inputInfo_initCheck( CCSMInit )

! !USES:

  use shr_string_mod,   only: shr_string_listIntersect

  implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(shr_inputInfo_InitType), intent(INOUT) :: CCSMInit    ! Output CCSM structure

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(shr_inputInfo_initCheck) '
    ! --- Intersection of over-ride list to master override list ------
    character(len=SHR_KIND_CL) :: OverrideList 
    integer :: lastchar                        ! Last character index
    integer :: rc                              ! Return code
    character(len=*), parameter :: F00 = "('"//LogPrefix//"',A,I3,A)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------
    
    if (CCSMInit%case_name==' ') then
       call shr_sys_abort( subname//': variable case_name must be set')
    end if
    ! --- Override values ------------------------------------------------------
    if ( len_trim(CCSMInit%restart_file_override) > 0 )then
       call shr_string_listIntersect(    &
                    CCSMInit%restart_file_override,restartFileOverRideList, &
                    OverRideList, rc)
       if ( trim(OverRideList) /= trim(CCSMInit%restart_file_override) )then
          write(6,FAA)'ERROR: the only values that can be overridden are: '// &
                      trim(restartFileOverRideList)
          call shr_sys_abort( subname//': list of values to override ' // &
                              'includes values that can NOT' &
                              //" be overridden" )
       end if
    end if
    ! --- Case name ------
    lastchar = len(CCSMInit%case_name)
    if ( len_trim(CCSMInit%case_name) == 0) then
       call shr_sys_abort( subname//': variable case_name must be set, aborting')
    end if
    if (CCSMInit%case_name(lastchar:lastchar) /= ' ') then
       write(6,F00)'ERROR: case_name must not exceed ', len(CCSMInit%case_name)-1, &
                 ' characters'
       call shr_sys_abort( subname//': variable case_name must be set, aborting')
    end if
    ! --- Special configurations ----- 
    if ( CCSMInit%atm_adiabatic .and. (CCSMInit%atm_ideal_phys .or. &
         CCSMInit%aqua_planet) )then
       call shr_sys_abort( subname//': only one of atm_adiabatic, ' // &
                           'atm_ideal_phys or aqua_planet can be set' )
    end if
    ! --- Restart pointer file -----
    if ( len_trim(CCSMInit%restart_pfile) == 0 ) then
       call shr_sys_abort( subname//': restart_pfile must be set' )
    end if
    ! --- Archive directory ----- 
    if ( len_trim(CCSMInit%archive_dir) == 0 ) then
       call shr_sys_abort( subname//': archive_dir must be set' )
    end if
    if ( index(CCSMInit%archive_dir,"/",back=.true.) /= &
         len_trim(CCSMInit%archive_dir) ) then
       call shr_sys_abort( subname//': archive_dir must end with a slash' )
    end if
    ! --- LogFile ending name ----- 
    if ( len_trim(CCSMInit%logFilePostFix) == 0 ) then
       call shr_sys_abort( subname//': logFilePostFix  must be set to something not blank' )
    end if
    ! --- Output path root directory ----- 
    if ( len_trim(CCSMInit%outPathRoot) == 0 ) then
       call shr_sys_abort( subname//': outPathRoot  must be set' )
    end if
    if ( index(CCSMInit%outPathRoot,"/",back=.true.) /= &
         len_trim(CCSMInit%outPathRoot) ) then
       call shr_sys_abort( subname//': outPathRoot must end with a slash' )
    end if
    ! --- Start-type ------
    if ( len_trim(CCSMInit%start_type) == 0 ) then
       call shr_sys_abort( subname//': variable start_type must be set')
    end if
    if ( (trim(CCSMInit%start_type) /= shr_inputInfo_start_type_start  )  &
    .and. (trim(CCSMInit%start_type) /= shr_inputInfo_start_type_cont ) &
    .and. (trim(CCSMInit%start_type) /= shr_inputInfo_start_type_brnch)    &
    ) then
       call shr_sys_abort( subname//': variable start_type can only be ' // &
                           'set to:' // &
                           shr_inputInfo_start_type_start//', '//    &
                           shr_inputInfo_start_type_cont//', or'// &
                           shr_inputInfo_start_type_brnch )
    end if

END SUBROUTINE shr_inputInfo_initCheck

!===============================================================================
!===============================================================================

END MODULE shr_inputInfo_mod
