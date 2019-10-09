module cam_history
  !-------------------------------------------------------------------------------------------
  ! 
  ! The cam_history module provides the user interface for CAM's history output capabilities.
  ! It maintains the lists of fields that are written to each history file, and the associated
  ! metadata for those fields such as descriptive names, physical units, time axis properties,
  ! etc.  It also contains the programmer interface which provides routines that are called from
  ! the physics and dynamics initialization routines to define the fields that are produced by
  ! the model and are available for output, and the routine that is called from the corresponding
  ! run method to add the field values into a history buffer so that they may be output to disk.
  !
  ! There are two special history files.  The initial file and the satellite track file.
  !
  ! Public functions/subroutines:
  !   addfld, add_default
  !   intht
  !   write_restart_history
  !   read_restart_history
  !   outfld
  !   wshist
  !   initialize_iop_history
  !-----------------------------------------------------------------------

   use shr_kind_mod,    only: r8 => shr_kind_r8, r4 => shr_kind_r4
   use shr_sys_mod,     only: shr_sys_flush
   use spmd_utils,      only: masterproc
   use ppgrid,          only: pcols, psubcols
   use cam_instance,    only: inst_suffix
   use filenames,       only: interpret_filename_spec
   use filenames,       only: ncdata, bnd_topo
   use cam_abortutils,  only: endrun

   use pio,          only: file_desc_t, var_desc_t, pio_setframe, pio_write,  &
                           pio_noerr, pio_bcast_error, pio_internal_error,    &
                           pio_seterrorhandling, pio_get_var, pio_clobber,    &
                           pio_int, pio_real, pio_double, pio_char,           &
                           pio_offset_kind, pio_unlimited, pio_global,        &
                           pio_inq_dimlen, pio_def_var, pio_enddef,           &
                           pio_put_att, pio_put_var, pio_get_att
                           

   use perf_mod,            only: t_startf, t_stopf
   use cam_logfile,         only: iulog	
   use cam_history_support, only: max_fieldname_len, fieldname_suffix_len,    &
                                  max_chars, ptapes, fieldname_len,           &
                                  max_string_len, date2yyyymmdd, pflds,       &
                                  fieldname_lenp2, sec2hms,                   &
                                  field_info, active_entry, hentry,           &
                                  horiz_only, write_hist_coord_attrs,         &
                                  write_hist_coord_vars, interp_info_t,       &
                                  lookup_hist_coord_indices, get_hist_coord_index
   use sat_hist,            only: is_satfile
   use mo_solar_parms,      only: solar_parms_get, solar_parms_on

  implicit none
  private
  save

  ! Forward common parameters to present unified interface to cam_history
  public :: fieldname_len, horiz_only

  !
  ! master_entry: elements of an entry in the master field list
  !
  type master_entry
    type (field_info)        :: field            ! field information
    character(len=max_fieldname_len) :: meridional_field = '' ! for vector fields
    character(len=max_fieldname_len) :: zonal_field = '' ! for vector fields
    character(len=1)         :: avgflag(ptapes)  ! averaging flag
    character(len=max_chars) :: time_op(ptapes)  ! time operator (e.g. max, min, avg)
    logical               :: act_sometape     ! Field is active on some tape
    logical               :: actflag(ptapes)  ! Per tape active/inactive flag
    integer               :: htapeindx(ptapes)! This field's index on particular history tape
    type(master_entry), pointer :: next_entry => null() ! The next master entry
  end type master_entry

  type (master_entry), pointer :: masterlinkedlist => null()   ! master field linkedlist top

  type master_list
    type(master_entry), pointer :: thisentry => null()
  end type master_list

  type (master_list), pointer :: masterlist(:) => null() ! master field array for hash lookup of field

  ! history tape info
  type (active_entry), pointer :: tape(:) => null()          ! history tapes
  type (active_entry), target,allocatable :: history_tape(:) ! history tapes
  type (active_entry), target, allocatable :: restarthistory_tape(:) ! restart history tapes

  type rvar_id
    type(var_desc_t), pointer :: vdesc => null()
    integer :: type
    integer :: ndims
    integer :: dims(4)
    character(len=fieldname_lenp2) :: name
  end type rvar_id
  type rdim_id
    integer :: len
    integer :: dimid
    character(len=fieldname_lenp2) :: name      
  end type rdim_id
  !
  !   The size of these parameters should match the assignments in restart_vars_setnames and restart_dims_setnames below
  !   
  integer, parameter :: restartvarcnt              = 37
  integer, parameter :: restartdimcnt              = 10
  type(rvar_id)      :: restartvars(restartvarcnt)
  type(rdim_id)      :: restartdims(restartdimcnt)
  integer, parameter :: ptapes_dim_ind             =  1
  integer, parameter :: max_string_len_dim_ind     =  2
  integer, parameter :: fieldname_lenp2_dim_ind    =  3
  integer, parameter :: pflds_dim_ind              =  4
  integer, parameter :: max_chars_dim_ind          =  5
  integer, parameter :: max_fieldname_len_dim_ind  =  6
  integer, parameter :: maxnflds_dim_ind           =  7
  integer, parameter :: maxvarmdims_dim_ind        =  8
  integer, parameter :: registeredmdims_dim_ind    =  9
  integer, parameter :: max_hcoordname_len_dim_ind = 10

  integer :: nfmaster = 0             ! number of fields in master field list
  integer :: nflds(ptapes)            ! number of fields per tape

  ! per tape sampling frequency (0=monthly avg)

  integer :: i                        ! index for nhtfrq initialization
  integer :: nhtfrq(ptapes) = (/0, (-24, i=2,ptapes)/)  ! history write frequency (0 = monthly)
  integer :: mfilt(ptapes) = 30       ! number of time samples per tape
  integer :: nfils(ptapes)            ! Array of no. of files on current h-file
  integer :: ndens(ptapes) = 2        ! packing density (double (1) or real (2))
  integer :: ncprec(ptapes) = -999    ! netcdf packing parameter based on ndens
  real(r8) :: beg_time(ptapes)        ! time at beginning of an averaging interval

  logical :: rgnht(ptapes) = .false.  ! flag array indicating regeneration volumes
  logical :: hstwr(ptapes) = .false.  ! Flag for history writes
  logical :: empty_htapes  = .false.  ! Namelist flag indicates no default history fields
  logical :: htapes_defined = .false. ! flag indicates history contents have been defined

  ! NB: This name must match the group name in namelist_definition.xml
  character(len=*), parameter   :: history_namelist = 'cam_history_nl'
  character(len=max_string_len) :: hrestpath(ptapes) = (/(' ',i=1,ptapes)/) ! Full history restart pathnames
  character(len=max_string_len) :: nfpath(ptapes) = (/(' ',i=1,ptapes)/) ! Array of first pathnames, for header
  character(len=max_string_len) :: cpath(ptapes)                   ! Array of current pathnames
  character(len=max_string_len) :: nhfil(ptapes)                   ! Array of current file names
  character(len=1)  :: avgflag_pertape(ptapes) = (/(' ',i=1,ptapes)/) ! per tape averaging flag
  character(len=16)  :: logname             ! user name
  character(len=16)  :: host                ! host name
  character(len=max_string_len) :: ctitle = ' '      ! Case title
  character(len=8)   :: inithist = 'YEARLY' ! If set to '6-HOURLY, 'DAILY', 'MONTHLY' or
  ! 'YEARLY' then write IC file 
  logical            :: inithist_all = .false. ! Flag to indicate set of fields to be 
                                          ! included on IC file
                                          !  .false.  include only required fields
                                          !  .true.   include required *and* optional fields
  character(len=fieldname_lenp2) :: fincl(pflds,ptapes) ! List of fields to add to primary h-file
  character(len=max_chars)       :: fincllonlat(pflds,ptapes) ! List of fields to add to primary h-file
  character(len=fieldname_lenp2) :: fexcl(pflds,ptapes) ! List of fields to rm from primary h-file
  character(len=fieldname_lenp2) :: fwrtpr(pflds,ptapes) ! List of fields to change default history output prec
  character(len=fieldname_suffix_len ) :: fieldname_suffix = '&IC' ! Suffix appended to field names for IC file

  ! Parameters for interpolated output tapes
  logical, public       :: interpolate_output(ptapes) = .false.
  ! The last two history files are not supported for interpolation
  type(interp_info_t)   :: interpolate_info(ptapes - 2)

  ! Allowed history averaging flags
  ! This should match namelist_definition.xml => avgflag_pertape (+ ' ')
  ! The presence of 'ABI' and 'XML' in this string is a coincidence
  character(len=7), parameter    :: HIST_AVG_FLAGS = ' ABIXML'
  character(len=22) ,parameter   :: LT_DESC = 'mean (over local time)' ! local time description
  logical :: collect_column_output(ptapes)

  integer :: maxvarmdims=1
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  Hashing.
  !
  !  Accelerate outfld processing by using a hash function of the field name
  !  to index masterlist and determine whehter the particular field is to
  !  be written to any history tape.
  !
  !
  !  Note: the outfld hashing logic will fail if any of the following are true:
  !
  !         1) The lower bound on the dimension of 'masterlist' is less than 1.
  !
  !         2) 'outfld' is called with field names that are not defined on
  !            masterlist.  This applies to both initial/branch and restart
  !            runs.
  !
  !         3) An inconsistency between a field's tape active flag
  !            'masterlist(ff)%actflag(t)' and active fields read from
  !            restart files.
  !
  !         4) Invoking function 'gen_hash_key' before the primary and secondary
  !            hash tables have been created (routine bld_outfld_hash_tbls).
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !
  !  User definable constants for hash and overflow tables.
  !  Define size of primary hash table (specified as 2**size).
  !
  integer, parameter :: tbl_hash_pri_sz_lg2 = 16
  !
  !  Define size of overflow hash table % of primary hash table.
  !
  integer, parameter :: tbl_hash_oflow_percent = 20
  !
  !  Do *not* modify the parameters below.
  !
  integer, parameter :: tbl_hash_pri_sz = 2**tbl_hash_pri_sz_lg2
  integer, parameter :: tbl_hash_oflow_sz = tbl_hash_pri_sz * (tbl_hash_oflow_percent/100.0_r8) 
  !
  !  The primary and overflow tables are organized to mimimize space (read:
  !  try to maximimze cache line usage).
  !
  !  gen_hash_key(fieldname) will return an index on the interval
  !  [0 ... tbl_hash_pri_sz-1].
  !
  !
  !  Primary:
  !  gen_hash_key(fieldname)-------+     +----------+
  !                                |     |   -ii    | 1 ------>tbl_hash_oflow(ii)
  !                                |     +----------+
  !                                +-->  |    ff    | 2 ------>masterlist(ff)
  !                                      +----------+
  !                                      |          | ...
  !                                      +----------+
  !                                      |          | tbl_hash_pri_sz
  !                                      +----------+
  !
  !  Overflow (if tbl_hash_pri() < 0):
  !  tbl_hash_pri(gen_hash_key(fieldname))
  !                         |
  !                         |            +----------+
  !                         |            |     1    | 1  (one entry on O.F. chain)
  !                         |            +----------+
  !                         |            |    ff_m  | 2
  !                         |            +----------+
  !                         +--------->  |     3    | 3  (three entries on chain)
  !                                      +----------+
  !                                      |    ff_x  | 4
  !                                      +----------+
  !                                      |    ff_y  | 5
  !                                      +----------+
  !                                      |    ff_z  | 6
  !                                      +----------+
  !                                      |          | ...
  !                                      +----------+
  !                                      |          | tbl_hash_oflow_sz
  !                                      +----------+
  !
  !
  integer, dimension(0:tbl_hash_pri_sz-1) :: tbl_hash_pri ! Primary hash table
  integer, dimension(tbl_hash_oflow_sz) :: tbl_hash_oflow ! Overflow hash table
  !
  !  Constants used in hashing function gen_hash_key.
  !  Note: if the constants in table 'tbl_gen_hash_key' below are modified,
  !        changes are required to routine 'gen_hash_key' because of specific
  !        logic in the routine that optimizes character strings of length 8.
  !

  integer, parameter :: gen_hash_key_offset = z'000053db'

  integer, parameter :: tbl_max_idx = 15  ! 2**N - 1
  integer, dimension(0:tbl_max_idx) :: tbl_gen_hash_key = &
       (/61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1/)

  !
  ! Filename specifiers for history, initial files and restart history files
  ! (%c = caseid, $y = year, $m = month, $d = day, $s = seconds in day, %t = tape number)
  !
  character(len=max_string_len) :: rhfilename_spec = '%c.cam.rh%t.%y-%m-%d-%s.nc' ! history restart
  character(len=max_string_len) :: hfilename_spec(ptapes) = (/ (' ', i=1, ptapes) /) ! filename specifyer


  interface addfld
    module procedure addfld_1d
    module procedure addfld_nd
  end interface

  ! needed by history_default
  public :: init_masterlinkedlist

  ! Needed by runtime_opts
  public :: hfilename_spec, mfilt, fincl, nhtfrq, avgflag_pertape, ctitle
  ! Needed by cam_diagnostics
  public :: inithist_all

  integer :: lcltod_start(ptapes) ! start time of day for local time averaging (sec)
  integer :: lcltod_stop(ptapes)  ! stop time of day for local time averaging, stop > start is wrap around (sec)

  ! Needed by stepon
  public :: hstwr
  public :: nfils

  ! Functions
  public :: history_readnl            ! Namelist reader for CAM history
  public :: init_restart_history      ! Write restart history data
  public :: write_restart_history     ! Write restart history data
  public :: read_restart_history      ! Read restart history data
  public :: wshist                    ! Write files out
  public :: outfld                    ! Output a field
  public :: intht                     ! Initialization
  public :: wrapup                    ! process history files at end of run
  public :: write_inithist            ! logical flag to allow dump of IC history buffer to IC file
  public :: addfld                    ! Add a field to history file
  public :: add_default               ! Add the default fields
  public :: register_vector_field     ! Register vector field set for interpolated output
  public :: get_hfilepath             ! Return history filename
  public :: get_ptapes                ! Return the number of tapes being used
  public :: get_hist_restart_filepath ! Return the full filepath to the history restart file
  public :: hist_fld_active           ! Determine if a field is active on any history file
  public :: hist_fld_col_active       ! Determine if a field is active on any history file at
  ! each column in a chunk

CONTAINS

  subroutine init_masterlinkedlist()

    nullify(masterlinkedlist)
    nullify(tape)

  end subroutine init_masterlinkedlist


  subroutine intht ()
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Initialize history file handler for initial or continuation run.
    !          For example, on an initial run, this routine initializes "ptapes"
    !          history files.  On a restart or regeneration  run, this routine 
    !          only initializes history files declared beyond what existed on the 
    !          previous run.  Files which already existed on the previous run have 
    !          already been initialized (i.e. named and opened) in routine RESTRT.
    ! 
    ! Method: Loop over tapes and fields per tape setting appropriate variables and
    !         calling appropriate routines
    ! 
    ! Author: CCM Core Group
    ! 
    !-----------------------------------------------------------------------
    use ioFileMod
    use shr_sys_mod,     only: shr_sys_getenv
    use time_manager,    only: get_prev_time, get_curr_time
    use cam_control_mod, only: nsrest
    use sat_hist,        only: sat_hist_init
#if (defined SPMD)
    use spmd_utils,     only: mpichar, mpicom, masterprocid, mpicom, mpi_character
#endif
    !
    !-----------------------------------------------------------------------
    !
    ! Local workspace
    !
    integer :: t, f              ! tape, field indices
    integer :: begdim1           ! on-node dim1 start index
    integer :: enddim1           ! on-node dim1 end index
    integer :: begdim2           ! on-node dim2 start index
    integer :: enddim2           ! on-node dim2 end index
    integer :: begdim3           ! on-node chunk or lat start index
    integer :: enddim3           ! on-node chunk or lat end index
    integer :: day, sec          ! day and seconds from base date
    integer :: rcode             ! shr_sys_getenv return code
    type(master_entry), pointer :: listentry
    character(len=32) :: fldname ! temp variable used to produce a left justified field name
    ! in the formatted logfile output

    !
    ! Print master field list
    !

    if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*)' ******* MASTER FIELD LIST *******'
    end if
    listentry=>masterlinkedlist
    f=0
    do while(associated(listentry))
      f=f+1
      if(masterproc) then
        fldname = listentry%field%name
        write(iulog,9000) f, fldname, listentry%field%units, listentry%field%numlev, &
             listentry%avgflag(1), trim(listentry%field%long_name)
9000    format(i5, 1x, a32, 1x, a16, 1x, i4, 1x, a1, 2x, a)
      end if
      listentry=>listentry%next_entry
    end do
    nfmaster = f
    if(masterproc) write(iulog,*)'intht:nfmaster=',nfmaster

    !
    !  Now that masterlinkedlist is defined and we are performing a restart run
    !  (after all addfld calls), construct primary and secondary hashing tables.
    !
    if (nsrest == 1) then
       call print_active_fldlst()
       call bld_outfld_hash_tbls()
       call bld_htapefld_indices()
       return
    end if
    !
    ! Get users logname and machine hostname
    !
    if ( masterproc )then
      logname = ' '
      call shr_sys_getenv ('LOGNAME',logname,rcode)
      host = ' '
      call shr_sys_getenv ('HOST',host,rcode)
    end if
#ifdef SPMD
    ! PIO requires netcdf attributes have consistant values on all tasks
    call mpi_bcast(logname, len(logname), mpi_character, masterprocid, mpicom, rcode)
    call mpi_bcast(host,    len(host),    mpi_character, masterprocid, mpicom, rcode)
#endif
    !
    ! Override averaging flag for all fields on a particular tape if namelist input so specifies
    !
    do t=1,ptapes
      if (avgflag_pertape(t) /= ' ') then
        call h_override (t)
      end if
    end do
    !
    ! Define field list information for all history files.  
    !
    call fldlst ()
    !
    ! Loop over max. no. of history files permitted  
    !
    if (nsrest == 3) then
      call get_prev_time(day, sec)  ! elapased time since reference date
    else
      call get_curr_time(day, sec)  ! elapased time since reference date
    end if
    do t=1,ptapes
      nfils(t) = 0            ! no. of time samples in hist. file no. t

      ! Time at beginning of current averaging interval.

      beg_time(t) = day + sec/86400._r8
    end do

    !
    ! Initialize history variables
    !
    do t=1,ptapes
      do f=1,nflds(t)
        begdim1  = tape(t)%hlist(f)%field%begdim1
        enddim1  = tape(t)%hlist(f)%field%enddim1
        begdim2  = tape(t)%hlist(f)%field%begdim2
        enddim2  = tape(t)%hlist(f)%field%enddim2
        begdim3  = tape(t)%hlist(f)%field%begdim3
        enddim3  = tape(t)%hlist(f)%field%enddim3
        allocate(tape(t)%hlist(f)%hbuf(begdim1:enddim1,begdim2:enddim2,begdim3:enddim3))
        tape(t)%hlist(f)%hbuf = 0._r8
        if(tape(t)%hlist(f)%field%flag_xyfill .or. (avgflag_pertape(t) .eq. 'L')) then
          allocate (tape(t)%hlist(f)%nacs(begdim1:enddim1,begdim3:enddim3))
        else
          allocate (tape(t)%hlist(f)%nacs(1,begdim3:enddim3))
        end if
        tape(t)%hlist(f)%nacs(:,:) = 0
        tape(t)%hlist(f)%field%meridional_complement = -1
        tape(t)%hlist(f)%field%zonal_complement = -1
      end do
    end do
    ! Setup vector pairs for unstructured grid interpolation
    call setup_interpolation_and_define_vector_complements()
    !  Initialize the sat following history subsystem
    call sat_hist_init()

    return
  end subroutine intht

  subroutine history_readnl(nlfile, dtime)

    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use spmd_utils,     only: masterproc, masterprocid, mpicom
    use spmd_utils,     only: mpi_integer, mpi_logical, mpi_character
    use shr_string_mod, only: shr_string_toUpper

    ! Dummy argument
    character(len=*), intent(in)   :: nlfile  ! filepath of namelist input file
    integer,          intent(in)   :: dtime   ! Step time in seconds

    !
    ! Local variables
    integer                        :: unitn, ierr, f, t
    character(len=8)               :: ctemp      ! Temporary character string

    character(len=fieldname_lenp2) :: fincl1(pflds)
    character(len=fieldname_lenp2) :: fincl2(pflds)
    character(len=fieldname_lenp2) :: fincl3(pflds)
    character(len=fieldname_lenp2) :: fincl4(pflds)
    character(len=fieldname_lenp2) :: fincl5(pflds)
    character(len=fieldname_lenp2) :: fincl6(pflds)
    character(len=fieldname_lenp2) :: fincl7(pflds)
    character(len=fieldname_lenp2) :: fincl8(pflds)
    character(len=fieldname_lenp2) :: fincl9(pflds)
    character(len=fieldname_lenp2) :: fincl10(pflds)

    character(len=max_chars)       :: fincl1lonlat(pflds)
    character(len=max_chars)       :: fincl2lonlat(pflds)
    character(len=max_chars)       :: fincl3lonlat(pflds)
    character(len=max_chars)       :: fincl4lonlat(pflds)
    character(len=max_chars)       :: fincl5lonlat(pflds)
    character(len=max_chars)       :: fincl6lonlat(pflds)
    character(len=max_chars)       :: fincl7lonlat(pflds)
    character(len=max_chars)       :: fincl8lonlat(pflds)
    character(len=max_chars)       :: fincl9lonlat(pflds)
    character(len=max_chars)       :: fincl10lonlat(pflds)

    character(len=fieldname_len)   :: fexcl1(pflds)
    character(len=fieldname_len)   :: fexcl2(pflds)
    character(len=fieldname_len)   :: fexcl3(pflds)
    character(len=fieldname_len)   :: fexcl4(pflds)
    character(len=fieldname_len)   :: fexcl5(pflds)
    character(len=fieldname_len)   :: fexcl6(pflds)
    character(len=fieldname_len)   :: fexcl7(pflds)
    character(len=fieldname_len)   :: fexcl8(pflds)
    character(len=fieldname_len)   :: fexcl9(pflds)
    character(len=fieldname_len)   :: fexcl10(pflds)

    character(len=fieldname_lenp2) :: fwrtpr1(pflds)
    character(len=fieldname_lenp2) :: fwrtpr2(pflds)
    character(len=fieldname_lenp2) :: fwrtpr3(pflds)
    character(len=fieldname_lenp2) :: fwrtpr4(pflds)
    character(len=fieldname_lenp2) :: fwrtpr5(pflds)
    character(len=fieldname_lenp2) :: fwrtpr6(pflds)
    character(len=fieldname_lenp2) :: fwrtpr7(pflds)
    character(len=fieldname_lenp2) :: fwrtpr8(pflds)
    character(len=fieldname_lenp2) :: fwrtpr9(pflds)
    character(len=fieldname_lenp2) :: fwrtpr10(pflds)

    integer                        :: interpolate_nlat(size(interpolate_info))
    integer                        :: interpolate_nlon(size(interpolate_info))
    integer                        :: interpolate_gridtype(size(interpolate_info))
    integer                        :: interpolate_type(size(interpolate_info))

    ! History namelist items
    namelist /cam_history_nl/ ndens, nhtfrq, mfilt, inithist, inithist_all,    &
         avgflag_pertape, empty_htapes, lcltod_start, lcltod_stop,             &
         fincl1lonlat, fincl2lonlat, fincl3lonlat, fincl4lonlat, fincl5lonlat, &
         fincl6lonlat, fincl7lonlat, fincl8lonlat, fincl9lonlat,               &
         fincl10lonlat, collect_column_output, hfilename_spec,                 &
         fincl1,  fincl2,  fincl3,  fincl4,  fincl5,                           &
         fincl6,  fincl7,  fincl8,  fincl9,  fincl10,                          &
         fexcl1,  fexcl2,  fexcl3,  fexcl4,  fexcl5,                           &
         fexcl6,  fexcl7,  fexcl8,  fexcl9,  fexcl10,                          &
         fwrtpr1, fwrtpr2, fwrtpr3, fwrtpr4, fwrtpr5,                          &
         fwrtpr6, fwrtpr7, fwrtpr8, fwrtpr9, fwrtpr10,                         &
         interpolate_nlat, interpolate_nlon,                                   &
         interpolate_gridtype, interpolate_type, interpolate_output

    ! Set namelist defaults (these should match initial values if given)
    fincl(:,:)               = ' '
    fincllonlat(:,:)         = ' '
    fexcl(:,:)               = ' '
    fwrtpr(:,:)              = ' '
    collect_column_output(:) = .false.
    avgflag_pertape(:)       = ' '
    ndens                    = 2
    nhtfrq(1)                = 0
    nhtfrq(2:)               = -24
    mfilt                    = 30
    inithist                 = 'YEARLY'
    inithist_all             = .false.
    empty_htapes             = .false.
    lcltod_start(:)          = 0
    lcltod_stop(:)           = 0
    hfilename_spec(:)        = ' '
    interpolate_nlat(:)      = 0
    interpolate_nlon(:)      = 0
    interpolate_gridtype(:)  = 1
    interpolate_type(:)      = 1
    interpolate_output(:)    = .false.

    ! Initialize namelist 'temporary variables'
    do f = 1, pflds
      fincl1(f)        = ' '         
      fincl2(f)        = ' '         
      fincl3(f)        = ' '         
      fincl4(f)        = ' '         
      fincl5(f)        = ' '         
      fincl6(f)        = ' '         
      fincl7(f)        = ' '         
      fincl8(f)        = ' '         
      fincl9(f)        = ' '         
      fincl10(f)       = ' '         
      fincl1lonlat(f)  = ' '
      fincl2lonlat(f)  = ' '
      fincl3lonlat(f)  = ' '
      fincl4lonlat(f)  = ' '
      fincl5lonlat(f)  = ' '
      fincl6lonlat(f)  = ' '
      fincl7lonlat(f)  = ' '
      fincl8lonlat(f)  = ' '
      fincl9lonlat(f)  = ' '
      fincl10lonlat(f) = ' '
      fexcl1(f)        = ' '
      fexcl2(f)        = ' '
      fexcl3(f)        = ' '
      fexcl4(f)        = ' '
      fexcl5(f)        = ' '
      fexcl6(f)        = ' '
      fexcl7(f)        = ' '
      fexcl8(f)        = ' '
      fexcl9(f)        = ' '
      fexcl10(f)       = ' '
      fwrtpr1(f)       = ' '
      fwrtpr2(f)       = ' '
      fwrtpr3(f)       = ' '
      fwrtpr4(f)       = ' '
      fwrtpr5(f)       = ' '
      fwrtpr6(f)       = ' '
      fwrtpr7(f)       = ' '
      fwrtpr8(f)       = ' '
      fwrtpr9(f)       = ' '
      fwrtpr10(f)      = ' '
    end do

    if (trim(history_namelist) /= 'cam_history_nl') then
      call endrun('HISTORY_READNL: CAM history namelist mismatch')
    end if
    if (masterproc) then
      write(iulog, *) 'Read in ',history_namelist,' namelist from: ',trim(nlfile)
      unitn = getunit()
      open(unitn, file=trim(nlfile), status='old')
      call find_group_name(unitn, history_namelist, status=ierr)
      if (ierr == 0) then
        read(unitn, cam_history_nl, iostat=ierr)
        if (ierr /= 0) then
          call endrun('history_readnl: ERROR reading namelist, '//trim(history_namelist))
        end if
      end if
      close(unitn)
      call freeunit(unitn)

      do f = 1, pflds
        fincl(f, 1) = fincl1(f)
        fincl(f, 2) = fincl2(f)
        fincl(f, 3) = fincl3(f)
        fincl(f, 4) = fincl4(f)
        fincl(f, 5) = fincl5(f)
        fincl(f, 6) = fincl6(f)
        fincl(f, 7) = fincl7(f)
        fincl(f, 8) = fincl8(f)
        fincl(f, 9) = fincl9(f)
        fincl(f,10) = fincl10(f)

        fincllonlat(f, 1) = fincl1lonlat(f)
        fincllonlat(f, 2) = fincl2lonlat(f)
        fincllonlat(f, 3) = fincl3lonlat(f)
        fincllonlat(f, 4) = fincl4lonlat(f)
        fincllonlat(f, 5) = fincl5lonlat(f)
        fincllonlat(f, 6) = fincl6lonlat(f)
        fincllonlat(f, 7) = fincl7lonlat(f)
        fincllonlat(f, 8) = fincl8lonlat(f)
        fincllonlat(f, 9) = fincl9lonlat(f)
        fincllonlat(f,10) = fincl10lonlat(f)

        fexcl(f, 1) = fexcl1(f)
        fexcl(f, 2) = fexcl2(f)
        fexcl(f, 3) = fexcl3(f)
        fexcl(f, 4) = fexcl4(f)
        fexcl(f, 5) = fexcl5(f)
        fexcl(f, 6) = fexcl6(f)
        fexcl(f, 7) = fexcl7(f)
        fexcl(f, 8) = fexcl8(f)
        fexcl(f, 9) = fexcl9(f)
        fexcl(f,10) = fexcl10(f)

        fwrtpr(f, 1) = fwrtpr1(f)
        fwrtpr(f, 2) = fwrtpr2(f)
        fwrtpr(f, 3) = fwrtpr3(f)
        fwrtpr(f, 4) = fwrtpr4(f)
        fwrtpr(f, 5) = fwrtpr5(f)
        fwrtpr(f, 6) = fwrtpr6(f)
        fwrtpr(f, 7) = fwrtpr7(f)
        fwrtpr(f, 8) = fwrtpr8(f)
        fwrtpr(f, 9) = fwrtpr9(f)
        fwrtpr(f,10) = fwrtpr10(f)
      end do

      !
      ! If generate an initial conditions history file as an auxillary tape:
      !
      ctemp = shr_string_toUpper(inithist) 
      inithist = trim(ctemp)
      if ( (inithist /= '6-HOURLY') .and. (inithist /= 'DAILY')  .and.        &
           (inithist /= 'MONTHLY')  .and. (inithist /= 'YEARLY') .and.        &
           (inithist /= 'CAMIOP')   .and. (inithist /= 'ENDOFRUN')) then
        inithist = 'NONE'
      end if
      ! 
      ! History file write times
      ! Convert write freq. of hist files from hours to timesteps if necessary.
      ! 
      do t = 1, ptapes
        if (nhtfrq(t) < 0) then
          nhtfrq(t) = nint((-nhtfrq(t) * 3600._r8) / dtime)
        end if
      end do
      !
      ! Initialize the filename specifier if not already set
      ! This is the format for the history filenames:
      ! %c= caseid, %t=tape no., %y=year, %m=month, %d=day, %s=second, %%=%
      ! See the filenames module for more information
      !
      do t = 1, ptapes
        if ( len_trim(hfilename_spec(t)) == 0 )then
          if ( nhtfrq(t) == 0 )then
            ! Monthly files
            hfilename_spec(t) = '%c.cam' // trim(inst_suffix) // '.h%t.%y-%m.nc'
          else
            hfilename_spec(t) = '%c.cam' // trim(inst_suffix) // '.h%t.%y-%m-%d-%s.nc'
          end if
        end if
        !
        ! Only one time sample allowed per monthly average file
        !
        if (nhtfrq(t) == 0) then
          mfilt(t) = 1
        end if
      end do
    end if ! masterproc

    ! Print per-tape averaging flags
    if (masterproc) then
      do t = 1, ptapes
        if (avgflag_pertape(t) /= ' ') then
          write(iulog,*)'Unless overridden by namelist input on a per-field basis (FINCL),'
          write(iulog,*)'All fields on history file ',t,' will have averaging flag ',avgflag_pertape(t)
        end if
        ! Enforce no interpolation for satellite files
        if (is_satfile(t) .and. interpolate_output(t)) then
          write(iulog, *) 'WARNING: Interpolated output not supported for a satellite history file, ignored'
          interpolate_output(t) = .false.
        end if
        ! Enforce no interpolation for IC files
        if (is_initfile(t) .and. interpolate_output(t)) then
          write(iulog, *) 'WARNING: Interpolated output not supported for a satellite history file, ignored'
          interpolate_output(t) = .false.
        end if
      end do
    end if

    ! Write out inithist info
    if (masterproc) then
      if (inithist == '6-HOURLY' ) then
        write(iulog,*)'Initial conditions history files will be written 6-hourly.'
      else if (inithist == 'DAILY' ) then
        write(iulog,*)'Initial conditions history files will be written daily.'
      else if (inithist == 'MONTHLY' ) then
        write(iulog,*)'Initial conditions history files will be written monthly.'
      else if (inithist == 'YEARLY' ) then
        write(iulog,*)'Initial conditions history files will be written yearly.'
      else if (inithist == 'CAMIOP' ) then
        write(iulog,*)'Initial conditions history files will be written for IOP.'
      else if (inithist == 'ENDOFRUN' ) then
        write(iulog,*)'Initial conditions history files will be written at end of run.'
      else
        write(iulog,*)'Initial conditions history files will not be created'
      end if
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpi_bcast(ndens, ptapes, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(nhtfrq, ptapes, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(mfilt, ptapes, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(inithist,len(inithist), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(inithist_all,1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(lcltod_start, ptapes, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(lcltod_stop,  ptapes, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(collect_column_output, ptapes, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(empty_htapes,1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(avgflag_pertape, ptapes, mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(hfilename_spec, len(hfilename_spec(1))*ptapes, mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(fincl, len(fincl (1,1))*pflds*ptapes, mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(fexcl, len(fexcl (1,1))*pflds*ptapes, mpi_character, masterprocid, mpicom, ierr)

    call mpi_bcast(fincllonlat, len(fincllonlat (1,1))*pflds*ptapes, mpi_character, masterprocid, mpicom, ierr)

    call mpi_bcast(fwrtpr, len(fwrtpr(1,1))*pflds*ptapes, mpi_character, masterprocid, mpicom, ierr)
    t = size(interpolate_nlat, 1)
    call mpi_bcast(interpolate_nlat, t, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(interpolate_nlon, t, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(interpolate_gridtype, t, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(interpolate_type, t, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(interpolate_output, ptapes, mpi_logical, masterprocid, mpicom, ierr)
#endif

    ! Setup the interpolate_info structures
    do t = 1, size(interpolate_info)
      interpolate_info(t)%interp_type = interpolate_type(t)
      interpolate_info(t)%interp_gridtype = interpolate_gridtype(t)
      interpolate_info(t)%interp_nlat = interpolate_nlat(t)
      interpolate_info(t)%interp_nlon = interpolate_nlon(t)
    end do

  end subroutine history_readnl

  subroutine set_field_dimensions(field)
    use cam_history_support, only: hist_coord_size
    use cam_grid_support,    only: cam_grid_get_array_bounds, cam_grid_is_block_indexed
    ! Dummy arguments
    type(field_info), intent(inout) :: field

    ! Local variables
    integer                         :: i
    integer                         :: msize
    integer                         :: dimbounds(2,2)

    call cam_grid_get_array_bounds(field%decomp_type, dimbounds)
    field%begdim1  = dimbounds(1,1)
    field%enddim1  = dimbounds(1,2)
    field%begdim2  = 1
    if (associated(field%mdims)) then
      if (size(field%mdims) > 0) then
        field%enddim2  = 1
        do i = 1, size(field%mdims)
          msize = hist_coord_size(field%mdims(i))
          if (msize <= 0) then
            call endrun('set_field_dimensions: mdim size must be > 0')
          end if
          field%enddim2 = field%enddim2 * msize
        end do
      else
        if (field%numlev < 1) then
          if (masterproc) then
            write(iulog, *) 'SET_FIELD_DIMENSIONS WARNING: illegal numlev for ', trim(field%name)
          end if
          field%numlev = 1
        end if
        field%enddim2 = field%numlev
      end if
    else
      if (field%numlev < 1) then
        if (masterproc) then
          write(iulog, *) 'SET_FIELD_DIMENSIONS WARNING: illegal numlev for ', trim(field%name)
        end if
        field%numlev = 1
      end if
      field%enddim2 = field%numlev
    end if
    field%begdim3  = dimbounds(2,1)
    field%enddim3  = dimbounds(2,2)
    field%colperchunk = cam_grid_is_block_indexed(field%decomp_type)

  end subroutine set_field_dimensions

  subroutine setup_interpolation_and_define_vector_complements()
    use interp_mod, only: setup_history_interpolation

    ! Local variables
    integer :: hf, f, ii, ff
    logical :: interp_ok
    character(len=max_fieldname_len) :: mname
    character(len=max_fieldname_len) :: zname
    character(len=*), parameter      :: subname='setup_interpolation_and_define_vector_complements'

    ! Do not interpolate IC history and sat hist files
    if (any(interpolate_output)) then
      call setup_history_interpolation(interp_ok, ptapes-2,                   &
           interpolate_output, interpolate_info)
      do hf = 1, ptapes - 2
        if((.not. is_satfile(hf)) .and. (.not. is_initfile(hf))) then
          do f = 1, nflds(hf)
            if (field_part_of_vector(trim(tape(hf)%hlist(f)%field%name),      &
                 mname, zname)) then
              if (len_trim(mname) > 0) then
                ! This field is a zonal part of a set, find the meridional partner
                do ff = 1, nflds(hf)
                  if (trim(mname) == trim(tape(hf)%hlist(ff)%field%name)) then
                    tape(hf)%hlist(f)%field%meridional_complement = ff
                    tape(hf)%hlist(ff)%field%zonal_complement     = f
                    exit
                  end if
                  if (ff == nflds(hf)) then
                    call endrun(trim(subname)//': No meridional match for '//trim(tape(hf)%hlist(f)%field%name))
                  end if
                end do
              else if (len_trim(zname) > 0) then
                ! This field is a meridional part of a set, find the zonal partner
                do ff = 1, nflds(hf)
                  if (trim(zname) == trim(tape(hf)%hlist(ff)%field%name)) then
                    tape(hf)%hlist(f)%field%zonal_complement       = ff
                    tape(hf)%hlist(ff)%field%meridional_complement = f
                    exit
                  end if
                  if (ff == nflds(hf)) then
                    call endrun(trim(subname)//': No zonal match for '//trim(tape(hf)%hlist(f)%field%name))
                  end if
                end do
              else
                call endrun(trim(subname)//': INTERNAL ERROR, bad vector field')
              end if
            end if
          end do
        end if
      end do
    end if
  end subroutine setup_interpolation_and_define_vector_complements

  subroutine restart_vars_setnames()

    ! Local variable
    integer :: rvindex

    rvindex = 1
    restartvars(rvindex)%name = 'rgnht'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1     
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'nhtfrq'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'nflds'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'nfils'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'mfilt'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'nfpath'
    restartvars(rvindex)%type = pio_char
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = max_string_len_dim_ind
    restartvars(rvindex)%dims(2) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'cpath'
    restartvars(rvindex)%type = pio_char
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = max_string_len_dim_ind
    restartvars(rvindex)%dims(2) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'nhfil'
    restartvars(rvindex)%type = pio_char
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = max_string_len_dim_ind
    restartvars(rvindex)%dims(2) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'ndens'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'fincllonlat'
    restartvars(rvindex)%type = pio_char
    restartvars(rvindex)%ndims = 3
    restartvars(rvindex)%dims(1) = max_chars_dim_ind
    restartvars(rvindex)%dims(2) = pflds_dim_ind
    restartvars(rvindex)%dims(3) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'ncprec'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'beg_time'
    restartvars(rvindex)%type = pio_double
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'fincl'
    restartvars(rvindex)%type = pio_char
    restartvars(rvindex)%ndims = 3
    restartvars(rvindex)%dims(1) = fieldname_lenp2_dim_ind
    restartvars(rvindex)%dims(2) = pflds_dim_ind
    restartvars(rvindex)%dims(3) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'fexcl'
    restartvars(rvindex)%type = pio_char
    restartvars(rvindex)%ndims = 3
    restartvars(rvindex)%dims(1) = fieldname_lenp2_dim_ind
    restartvars(rvindex)%dims(2) = pflds_dim_ind
    restartvars(rvindex)%dims(3) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'field_name'
    restartvars(rvindex)%type = pio_char
    restartvars(rvindex)%ndims = 3
    restartvars(rvindex)%dims(1) = max_fieldname_len_dim_ind
    restartvars(rvindex)%dims(2) = maxnflds_dim_ind
    restartvars(rvindex)%dims(3) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'decomp_type'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = maxnflds_dim_ind
    restartvars(rvindex)%dims(2) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'numlev'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = maxnflds_dim_ind
    restartvars(rvindex)%dims(2) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'hrestpath'
    restartvars(rvindex)%type = pio_char
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = max_string_len_dim_ind
    restartvars(rvindex)%dims(2) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'hwrt_prec'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = maxnflds_dim_ind
    restartvars(rvindex)%dims(2) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'avgflag'
    restartvars(rvindex)%type = pio_char
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = maxnflds_dim_ind
    restartvars(rvindex)%dims(2) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'sampling_seq'
    restartvars(rvindex)%type = pio_char
    restartvars(rvindex)%ndims = 3
    restartvars(rvindex)%dims(1) = max_chars_dim_ind
    restartvars(rvindex)%dims(2) = maxnflds_dim_ind
    restartvars(rvindex)%dims(3) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'long_name'
    restartvars(rvindex)%type = pio_char
    restartvars(rvindex)%ndims = 3
    restartvars(rvindex)%dims(1) = max_chars_dim_ind
    restartvars(rvindex)%dims(2) = maxnflds_dim_ind
    restartvars(rvindex)%dims(3) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'units'
    restartvars(rvindex)%type = pio_char
    restartvars(rvindex)%ndims = 3
    restartvars(rvindex)%dims(1) = max_chars_dim_ind
    restartvars(rvindex)%dims(2) = maxnflds_dim_ind
    restartvars(rvindex)%dims(3) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'xyfill'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = maxnflds_dim_ind
    restartvars(rvindex)%dims(2) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'lcltod_start'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'lcltod_stop'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'fillvalue'
    restartvars(rvindex)%type = pio_double
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = maxnflds_dim_ind
    restartvars(rvindex)%dims(2) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'mdims'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 3
    restartvars(rvindex)%dims(1) = maxvarmdims_dim_ind
    restartvars(rvindex)%dims(2) = maxnflds_dim_ind
    restartvars(rvindex)%dims(3) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'mdimnames'
    restartvars(rvindex)%type = pio_char
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = max_hcoordname_len_dim_ind
    restartvars(rvindex)%dims(2) = registeredmdims_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'is_subcol'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = maxnflds_dim_ind
    restartvars(rvindex)%dims(2) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'interpolate_output'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'interpolate_type'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'interpolate_gridtype'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'interpolate_nlat'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'interpolate_nlon'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 1
    restartvars(rvindex)%dims(1) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'meridional_complement'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = maxnflds_dim_ind
    restartvars(rvindex)%dims(2) = ptapes_dim_ind

    rvindex = rvindex + 1
    restartvars(rvindex)%name = 'zonal_complement'
    restartvars(rvindex)%type = pio_int
    restartvars(rvindex)%ndims = 2
    restartvars(rvindex)%dims(1) = maxnflds_dim_ind
    restartvars(rvindex)%dims(2) = ptapes_dim_ind

  end subroutine restart_vars_setnames

  subroutine restart_dims_setnames()
    use cam_grid_support,    only: max_hcoordname_len
    use cam_history_support, only: registeredmdims

    restartdims(ptapes_dim_ind)%name = 'ptapes'
    restartdims(ptapes_dim_ind)%len  = ptapes

    restartdims(max_string_len_dim_ind)%name = 'max_string_len'
    restartdims(max_string_len_dim_ind)%len  = max_string_len

    restartdims(fieldname_lenp2_dim_ind)%name = 'fieldname_lenp2'
    restartdims(fieldname_lenp2_dim_ind)%len  = fieldname_lenp2

    restartdims(pflds_dim_ind)%name = 'pflds'
    restartdims(pflds_dim_ind)%len  = pflds

    restartdims(max_chars_dim_ind)%name = 'max_chars'
    restartdims(max_chars_dim_ind)%len  = max_chars

    restartdims(max_fieldname_len_dim_ind)%name = 'max_fieldname_len'
    restartdims(max_fieldname_len_dim_ind)%len  = max_fieldname_len

    restartdims(maxnflds_dim_ind)%name = 'maxnflds'
    restartdims(maxnflds_dim_ind)%len  = maxval(nflds)

    restartdims(maxvarmdims_dim_ind)%name = 'maxvarmdims'
    restartdims(maxvarmdims_dim_ind)%len  = maxvarmdims

    restartdims(registeredmdims_dim_ind)%name = 'registeredmdims'
    restartdims(registeredmdims_dim_ind)%len  = registeredmdims

    restartdims(max_hcoordname_len_dim_ind)%name = 'max_hcoordname_len'
    restartdims(max_hcoordname_len_dim_ind)%len  = max_hcoordname_len

  end subroutine restart_dims_setnames


  subroutine init_restart_history (File)
    use cam_pio_utils,  only: cam_pio_def_dim
    use cam_pio_utils,  only: cam_pio_handle_error

    !---------------------------------------------------------------------------
    !
    ! Arguments
    !
    type(file_desc_t), intent(inout) :: File                 ! Pio file Handle
    !
    ! Local 
    !
    integer :: dimids(4), ndims
    integer :: ierr, i, k

    ! Don't need to write restart data if we have written the file this step
    where (hstwr(:))
      rgnht(:) = .false.
    elsewhere
      rgnht(:) = .true.
    end where

    if(maxval(nflds)>0) then
      call restart_vars_setnames()
      call restart_dims_setnames()

      do i=1,restartdimcnt
        ! it's possible that one or more of these have been defined elsewhere
        call cam_pio_def_dim(File, restartdims(i)%name, restartdims(i)%len,   &
             restartdims(i)%dimid, existOK=.true.)
      end do

      do i=1,restartvarcnt
        ndims= restartvars(i)%ndims
        do k=1,ndims
          dimids(k)=restartdims(restartvars(i)%dims(k))%dimid
        end do
        allocate(restartvars(i)%vdesc)
        ierr = pio_def_var(File, restartvars(i)%name, restartvars(i)%type, dimids(1:ndims), restartvars(i)%vdesc)
        call cam_pio_handle_error(ierr, 'INIT_RESTART_HISTORY: Error defining '//trim(restartvars(i)%name))

      end do
    end if
  end subroutine init_restart_history

  function restartvar_getdesc(name) result(vdesc)
    character(len=*), intent(in) :: name
    type(var_desc_t), pointer :: vdesc
    character(len=max_chars) :: errmsg
    integer :: i

    nullify(vdesc)
    do i=1,restartvarcnt
      if(name .eq. restartvars(i)%name) then
        vdesc=>restartvars(i)%vdesc
        exit
      end if
    end do
    if(.not.associated(vdesc)) then
      errmsg = 'Could not find restart variable '//name
      call endrun(errmsg)
    end if
  end function restartvar_getdesc


  !#######################################################################

  subroutine write_restart_history ( File, & 
       yr_spec, mon_spec, day_spec, sec_spec )
    use cam_history_support, only: hist_coord_name, registeredmdims

    implicit none
    !--------------------------------------------------------------------------------------------------
    !
    ! Arguments
    !
    type(file_desc_t), intent(inout) :: file         ! PIO restart file pointer
    integer, intent(in), optional :: yr_spec         ! Simulation year
    integer, intent(in), optional :: mon_spec        ! Simulation month
    integer, intent(in), optional :: day_spec        ! Simulation day
    integer, intent(in), optional :: sec_spec        ! Seconds into current simulation day
    !
    ! Local workspace
    !
    integer :: ierr, t, f
    integer :: rgnht_int(ptapes), start(2), startc(3)
    type(var_desc_t), pointer :: vdesc

    ! PIO variable descriptors
    type(var_desc_t), pointer ::  field_name_desc   ! Restart field names
    type(var_desc_t), pointer ::  decomp_type_desc
    type(var_desc_t), pointer ::  numlev_desc
    type(var_desc_t), pointer ::  avgflag_desc
    type(var_desc_t), pointer ::  sseq_desc
    type(var_desc_t), pointer ::  longname_desc
    type(var_desc_t), pointer ::  units_desc
    type(var_desc_t), pointer ::  hwrt_prec_desc
    type(var_desc_t), pointer ::  xyfill_desc
    type(var_desc_t), pointer ::  mdims_desc        ! mdim name indices
    type(var_desc_t), pointer ::  mdimname_desc     ! mdim names
    type(var_desc_t), pointer ::  issubcol_desc
    type(var_desc_t), pointer ::  fillval_desc
    type(var_desc_t), pointer ::  interpolate_output_desc
    type(var_desc_t), pointer ::  interpolate_type_desc
    type(var_desc_t), pointer ::  interpolate_gridtype_desc
    type(var_desc_t), pointer ::  interpolate_nlat_desc
    type(var_desc_t), pointer ::  interpolate_nlon_desc
    type(var_desc_t), pointer ::  meridional_complement_desc
    type(var_desc_t), pointer ::  zonal_complement_desc

    integer, allocatable      ::  allmdims(:,:,:)
    integer, allocatable      ::  xyfill(:,:)
    integer, allocatable      ::  is_subcol(:,:)
    integer, allocatable      ::  interp_output(:)

    integer                   ::  maxnflds


    maxnflds = maxval(nflds)
    allocate(xyfill(maxnflds, ptapes))
    xyfill = 0
    allocate(is_subcol(maxnflds, ptapes))
    is_subcol = 0
    allocate(interp_output(ptapes))
    interp_output = 0

    !
    !-----------------------------------------------------------------------
    ! Write the history restart data if necessary
    !-----------------------------------------------------------------------

    rgnht_int(:) = 0

    if(.not.allocated(restarthistory_tape)) allocate(restarthistory_tape(ptapes))

    do t=1,ptapes
      ! No need to write history IC restart because it is always instantaneous
      if (is_initfile(file_index=t)) rgnht(t) = .false.
      ! No need to write restart data for empty files
      if (nflds(t) == 0) rgnht(t) = .false.
      if(rgnht(t)) then
        rgnht_int(t) = 1
        restarthistory_tape(t)%hlist => history_tape(t)%hlist

        if(associated(history_tape(t)%grid_ids)) then
          restarthistory_tape(t)%grid_ids => history_tape(t)%grid_ids
        end if
        if(associated(history_tape(t)%patches)) then
          restarthistory_tape(t)%patches => history_tape(t)%patches
        end if
      end if
    end do

    if(maxval(nflds)<=0) return

    call wshist(rgnht)

    vdesc => restartvar_getdesc('fincl')
    ierr= pio_put_var(File, vdesc, fincl(:,1:ptapes))

    vdesc => restartvar_getdesc('fincllonlat')
    ierr= pio_put_var(File, vdesc, fincllonlat(:,1:ptapes))

    vdesc => restartvar_getdesc('fexcl')
    ierr= pio_put_var(File, vdesc, fexcl(:,1:ptapes))

    vdesc => restartvar_getdesc('rgnht')
    ierr= pio_put_var(File, vdesc, rgnht_int(1:ptapes))

    vdesc => restartvar_getdesc('nhtfrq')
    ierr= pio_put_var(File, vdesc, nhtfrq(1:ptapes))

    vdesc => restartvar_getdesc('nflds')
    ierr= pio_put_var(File, vdesc, nflds(1:ptapes))

    vdesc => restartvar_getdesc('nfils')
    ierr= pio_put_var(File, vdesc, nfils(1:ptapes))

    vdesc => restartvar_getdesc('mfilt')
    ierr= pio_put_var(File, vdesc, mfilt(1:ptapes))

    vdesc => restartvar_getdesc('nfpath')
    ierr= pio_put_var(File, vdesc, nfpath(1:ptapes))

    vdesc => restartvar_getdesc('cpath')
    ierr= pio_put_var(File, vdesc,  cpath(1:ptapes))

    vdesc => restartvar_getdesc('nhfil')
    ierr= pio_put_var(File, vdesc, nhfil(1:ptapes))

    vdesc => restartvar_getdesc('ndens')
    ierr= pio_put_var(File, vdesc, ndens(1:ptapes))
    vdesc => restartvar_getdesc('ncprec')
    ierr= pio_put_var(File, vdesc, ncprec(1:ptapes))
    vdesc => restartvar_getdesc('beg_time')
    ierr= pio_put_var(File, vdesc, beg_time(1:ptapes))

    vdesc => restartvar_getdesc('hrestpath')
    ierr = pio_put_var(File, vdesc, hrestpath(1:ptapes))

    vdesc => restartvar_getdesc('lcltod_start')
    ierr = pio_put_var(File, vdesc, lcltod_start(1:ptapes))

    vdesc => restartvar_getdesc('lcltod_stop')
    ierr = pio_put_var(File, vdesc, lcltod_stop(1:ptapes))

    field_name_desc => restartvar_getdesc('field_name')
    decomp_type_desc => restartvar_getdesc('decomp_type')
    numlev_desc => restartvar_getdesc('numlev')
    hwrt_prec_desc => restartvar_getdesc('hwrt_prec')

    sseq_desc => restartvar_getdesc('sampling_seq')
    longname_desc => restartvar_getdesc('long_name')
    units_desc => restartvar_getdesc('units')
    avgflag_desc => restartvar_getdesc('avgflag')
    xyfill_desc => restartvar_getdesc('xyfill')
    issubcol_desc => restartvar_getdesc('is_subcol')

    interpolate_output_desc => restartvar_getdesc('interpolate_output')
    interpolate_type_desc => restartvar_getdesc('interpolate_type')
    interpolate_gridtype_desc => restartvar_getdesc('interpolate_gridtype')
    interpolate_nlat_desc => restartvar_getdesc('interpolate_nlat')
    interpolate_nlon_desc => restartvar_getdesc('interpolate_nlon')

    meridional_complement_desc => restartvar_getdesc('meridional_complement')
    zonal_complement_desc => restartvar_getdesc('zonal_complement')

    mdims_desc => restartvar_getdesc('mdims')
    mdimname_desc => restartvar_getdesc('mdimnames')
    fillval_desc => restartvar_getdesc('fillvalue')

    tape=>history_tape

    ! allmdims specifies the mdim indices for each field
    allocate(allmdims(maxvarmdims,maxval(nflds),ptapes))
    allmdims=-1

    startc(1)=1
    do t = 1,ptapes
      start(2)=t
      startc(3)=t
      do f=1,nflds(t)
        start(1)=f
        startc(2)=f
        ierr = pio_put_var(File, field_name_desc,startc,tape(t)%hlist(f)%field%name)
        ierr = pio_put_var(File, decomp_type_desc,start,tape(t)%hlist(f)%field%decomp_type)
        ierr = pio_put_var(File, numlev_desc,start,tape(t)%hlist(f)%field%numlev)

        ierr = pio_put_var(File, hwrt_prec_desc,start,tape(t)%hlist(f)%hwrt_prec)
        ierr = pio_put_var(File, sseq_desc,startc,tape(t)%hlist(f)%field%sampling_seq)
        ierr = pio_put_var(File, longname_desc,startc,tape(t)%hlist(f)%field%long_name)
        ierr = pio_put_var(File, units_desc,startc,tape(t)%hlist(f)%field%units)
        ierr = pio_put_var(File, avgflag_desc,start, tape(t)%hlist(f)%avgflag)

        ierr = pio_put_var(File, fillval_desc,start, tape(t)%hlist(f)%field%fillvalue)
        ierr = pio_put_var(File, meridional_complement_desc,start, tape(t)%hlist(f)%field%meridional_complement)
        ierr = pio_put_var(File, zonal_complement_desc,start, tape(t)%hlist(f)%field%zonal_complement)
        if(associated(tape(t)%hlist(f)%field%mdims)) then
          allmdims(1:size(tape(t)%hlist(f)%field%mdims),f,t) = tape(t)%hlist(f)%field%mdims
        else
        end if
        if(tape(t)%hlist(f)%field%flag_xyfill) then
           xyfill(f,t) = 1
        end if
        if(tape(t)%hlist(f)%field%is_subcol) then
           is_subcol(f,t) = 1
        end if
      end do
      if (interpolate_output(t)) then
        interp_output(t) = 1
      end if
    end do
    ierr = pio_put_var(File, xyfill_desc, xyfill)
    ierr = pio_put_var(File, mdims_desc, allmdims)
    ierr = pio_put_var(File, issubcol_desc, is_subcol)
    !! Interpolated output variables
    ierr = pio_put_var(File, interpolate_output_desc, interp_output)
    interp_output = 1
    do t = 1, size(interpolate_info)
      interp_output(t) = interpolate_info(t)%interp_type
    end do
    ierr = pio_put_var(File, interpolate_type_desc, interp_output)
    interp_output = 1
    do t = 1, size(interpolate_info)
      interp_output(t) = interpolate_info(t)%interp_gridtype
    end do
    ierr = pio_put_var(File, interpolate_gridtype_desc, interp_output)
    interp_output = 0
    do t = 1, size(interpolate_info)
      interp_output(t) = interpolate_info(t)%interp_nlat
    end do
    ierr = pio_put_var(File, interpolate_nlat_desc, interp_output)
    interp_output = 0
    do t = 1, size(interpolate_info)
      interp_output(t) = interpolate_info(t)%interp_nlon
    end do
    ierr = pio_put_var(File, interpolate_nlon_desc, interp_output)
    ! Registered history coordinates
    start(1) = 1
    do f = 1, registeredmdims
      start(2) = f
      ierr = pio_put_var(File, mdimname_desc, start, hist_coord_name(f))
    end do

    deallocate(xyfill, allmdims)
    return

  end subroutine write_restart_history


  !#######################################################################

  subroutine read_restart_history (File)
    use pio,                 only: pio_inq_dimid
    use pio,                 only: pio_inq_varid, pio_inq_dimname
    use cam_pio_utils,       only: cam_pio_openfile, cam_pio_closefile
    use cam_pio_utils,       only: cam_pio_var_info
    use ioFileMod,           only: getfil
    use sat_hist,            only: sat_hist_define, sat_hist_init
    use cam_grid_support,    only: cam_grid_read_dist_array, cam_grid_num_grids
    use cam_history_support, only: get_hist_coord_index, add_hist_coord

    use shr_sys_mod,         only: shr_sys_getenv
#if (defined SPMD)
    use spmd_utils,          only: mpicom, mpichar
#endif
    !
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(file_desc_t), intent(inout) :: File            ! unit number
    !
    ! Local workspace
    !
    integer t, f, ff                 ! tape, field indices
    integer begdim2                  ! on-node vert start index
    integer enddim2                  ! on-node vert end index
    integer begdim1                  ! on-node dim1 start index
    integer enddim1                  ! on-node dim1 end index
    integer begdim3                  ! on-node chunk or lat start index
    integer enddim3                  ! on-node chunk or lat end index


    integer rgnht_int(ptapes)
    integer :: ierr

    character(len=max_string_len)  :: locfn       ! Local filename
    character(len=max_fieldname_len), allocatable :: tmpname(:,:)
    integer, allocatable :: decomp(:,:), tmpnumlev(:,:)
    integer, pointer :: nacs(:,:)    ! accumulation counter
    character(len=max_fieldname_len) :: fname_tmp ! local copy of field name
    character(len=max_fieldname_len) :: dname_tmp ! local copy of dim name

    integer :: i, ptapes_dimid

    type(var_desc_t)                 :: vdesc
    type(var_desc_t)                 :: longname_desc
    type(var_desc_t)                 :: units_desc
    type(var_desc_t)                 :: avgflag_desc
    type(var_desc_t)                 :: sseq_desc
    type(var_desc_t)                 :: fillval_desc
    type(var_desc_t)                 :: meridional_complement_desc
    type(var_desc_t)                 :: zonal_complement_desc
    integer,            allocatable  :: tmpprec(:,:)
    integer,            allocatable  :: xyfill(:,:)
    integer,            allocatable  :: allmdims(:,:,:)
    integer,            allocatable  :: is_subcol(:,:)
    integer,            allocatable  :: interp_output(:)
    integer                          :: nacsdimcnt, nacsval
    integer                          :: maxnflds, dimid

    ! List of active grids (first dim) for each tape (second dim)
    ! An active grid is one for which there is a least one field being output
    !    on that grid.
    integer, allocatable        :: gridsontape(:,:)

    character(len=16),  allocatable  :: mdimnames(:) ! Names of all hist coords (inc. vertical)
    integer                          :: ndims, dimids(8)
    integer                          :: tmpdims(8), dimcnt
    integer                          :: dimlens(7)
    integer                          :: mtapes, mdimcnt
    integer                          :: fdims(3)         ! Field dims
    integer                          :: nfdims           ! 2 or 3 (for 2D,3D)
    integer                          :: fdecomp          ! Grid ID for field

    !
    ! Get users logname and machine hostname
    !
    if ( masterproc )then
      logname = ' '
      call shr_sys_getenv ('LOGNAME',logname,ierr)
      host = ' '
      call shr_sys_getenv ('HOST',host,ierr)
    end if
#ifdef SPMD
    ! PIO requires netcdf attributes have consistant values on all tasks
    call mpibcast(logname, len(logname), mpichar, 0, mpicom)
    call mpibcast(host, len(host), mpichar, 0, mpicom)
#endif

    call pio_seterrorhandling(File, PIO_BCAST_ERROR)

    ierr = pio_inq_dimid(File, 'ptapes', ptapes_dimid)
    if(ierr/= PIO_NOERR) then
      if(masterproc) write(iulog,*) 'Not reading history info from restart file', ierr
      return   ! no history info in restart file
    end if
    call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)

    ierr = pio_inq_dimlen(File, ptapes_dimid, mtapes)

    ierr = pio_inq_dimid(File, 'maxnflds', dimid)
    ierr = pio_inq_dimlen(File, dimid, maxnflds)

    ierr = pio_inq_dimid(File, 'maxvarmdims', dimid)
    ierr = pio_inq_dimlen(File, dimid, maxvarmdims)

    ierr = pio_inq_varid(File, 'rgnht', vdesc)
    ierr = pio_get_var(File, vdesc, rgnht_int(1:mtapes))      

    ierr = pio_inq_varid(File, 'nhtfrq', vdesc)
    ierr = pio_get_var(File, vdesc, nhtfrq(1:mtapes))

    ierr = pio_inq_varid(File, 'nflds', vdesc)
    ierr = pio_get_var(File, vdesc, nflds(1:mtapes))
    ierr = pio_inq_varid(File, 'nfils', vdesc)
    ierr = pio_get_var(File, vdesc, nfils(1:mtapes))
    ierr = pio_inq_varid(File, 'mfilt', vdesc)
    ierr = pio_get_var(File, vdesc, mfilt(1:mtapes))

    ierr = pio_inq_varid(File, 'nfpath', vdesc)
    ierr = pio_get_var(File, vdesc, nfpath(1:mtapes))
    ierr = pio_inq_varid(File, 'cpath', vdesc)
    ierr = pio_get_var(File, vdesc, cpath(1:mtapes))
    ierr = pio_inq_varid(File, 'nhfil', vdesc)
    ierr = pio_get_var(File, vdesc, nhfil(1:mtapes))
    ierr = pio_inq_varid(File, 'hrestpath', vdesc)
    ierr = pio_get_var(File, vdesc, hrestpath(1:mtapes))


    ierr = pio_inq_varid(File, 'ndens', vdesc)
    ierr = pio_get_var(File, vdesc, ndens(1:mtapes))
    ierr = pio_inq_varid(File, 'ncprec', vdesc)
    ierr = pio_get_var(File, vdesc, ncprec(1:mtapes))
    ierr = pio_inq_varid(File, 'beg_time', vdesc)
    ierr = pio_get_var(File, vdesc, beg_time(1:mtapes))


    ierr = pio_inq_varid(File, 'fincl', vdesc)
    ierr = pio_get_var(File, vdesc, fincl(:,1:mtapes))

    ierr = pio_inq_varid(File, 'fincllonlat', vdesc)
    ierr = pio_get_var(File, vdesc, fincllonlat(:,1:mtapes))

    ierr = pio_inq_varid(File, 'fexcl', vdesc)
    ierr = pio_get_var(File, vdesc, fexcl(:,1:mtapes))

    ierr = pio_inq_varid(File, 'lcltod_start', vdesc)
    ierr = pio_get_var(File, vdesc, lcltod_start(1:mtapes))

    ierr = pio_inq_varid(File, 'lcltod_stop', vdesc)
    ierr = pio_get_var(File, vdesc, lcltod_stop(1:mtapes))




    allocate(tmpname(maxnflds, mtapes), decomp(maxnflds, mtapes), tmpnumlev(maxnflds,mtapes))
    ierr = pio_inq_varid(File, 'field_name', vdesc)
    ierr = pio_get_var(File, vdesc, tmpname)

    ierr = pio_inq_varid(File, 'decomp_type', vdesc)
    ierr = pio_get_var(File, vdesc, decomp)
    ierr = pio_inq_varid(File, 'numlev', vdesc)
    ierr = pio_get_var(File, vdesc, tmpnumlev)

    allocate(tmpprec(maxnflds,mtapes))
    ierr = pio_inq_varid(File, 'hwrt_prec',vdesc)
    ierr = pio_get_var(File, vdesc, tmpprec(:,:))

    allocate(xyfill(maxnflds,mtapes))
    ierr = pio_inq_varid(File, 'xyfill', vdesc)
    ierr = pio_get_var(File, vdesc, xyfill)

    allocate(is_subcol(maxnflds,mtapes))
    ierr = pio_inq_varid(File, 'is_subcol', vdesc)
    ierr = pio_get_var(File, vdesc, is_subcol)

    !! interpolated output
    allocate(interp_output(mtapes))
    ierr = pio_inq_varid(File, 'interpolate_output', vdesc)
    ierr = pio_get_var(File, vdesc, interp_output)
    interpolate_output(1:mtapes) = interp_output(1:mtapes) > 0
    if (ptapes > mtapes) then
      interpolate_output(mtapes+1:ptapes) = .false.
    end if
    ierr = pio_inq_varid(File, 'interpolate_type', vdesc)
    ierr = pio_get_var(File, vdesc, interp_output)
    do t = 1, mtapes
      if (interpolate_output(t)) then
        interpolate_info(t)%interp_type = interp_output(t)
      end if
    end do
    ierr = pio_inq_varid(File, 'interpolate_gridtype', vdesc)
    ierr = pio_get_var(File, vdesc, interp_output)
    do t = 1, mtapes
      if (interpolate_output(t)) then
        interpolate_info(t)%interp_gridtype = interp_output(t)
      end if
    end do
    ierr = pio_inq_varid(File, 'interpolate_nlat', vdesc)
    ierr = pio_get_var(File, vdesc, interp_output)
    do t = 1, mtapes
      if (interpolate_output(t)) then
        interpolate_info(t)%interp_nlat = interp_output(t)
      end if
    end do
    ierr = pio_inq_varid(File, 'interpolate_nlon', vdesc)
    ierr = pio_get_var(File, vdesc, interp_output)
    do t = 1, mtapes
      if (interpolate_output(t)) then
        interpolate_info(t)%interp_nlon = interp_output(t)
      end if
    end do

    !! mdim indices
    allocate(allmdims(maxvarmdims,maxnflds,mtapes))
    ierr = pio_inq_varid(File, 'mdims', vdesc)
    ierr = pio_get_var(File, vdesc, allmdims)

    !! mdim names
    ! Read the hist coord names to make sure they are all registered
    ierr = pio_inq_varid(File, 'mdimnames', vdesc)
    call cam_pio_var_info(File, vdesc, ndims, dimids, dimlens)
    mdimcnt = dimlens(2)
    allocate(mdimnames(mdimcnt))
    ierr = pio_get_var(File, vdesc, mdimnames)
    do f = 1, mdimcnt
      ! Check to see if the mdim is registered
      if (get_hist_coord_index(trim(mdimnames(f))) <= 0) then
        ! We need to register this mdim (hist_coord)
        call add_hist_coord(trim(mdimnames(f)))
      end if
    end do

    ierr = pio_inq_varid(File, 'avgflag', avgflag_desc)

    ierr = pio_inq_varid(File, 'long_name', longname_desc)
    ierr = pio_inq_varid(File, 'units', units_desc)
    ierr = pio_inq_varid(File, 'sampling_seq', sseq_desc)

    ierr = pio_inq_varid(File, 'fillvalue', fillval_desc)
    ierr = pio_inq_varid(File, 'meridional_complement', meridional_complement_desc)
    ierr = pio_inq_varid(File, 'zonal_complement', zonal_complement_desc)

    rgnht(:)=.false.

    allocate(history_tape(mtapes))

    tape => history_tape

    do t=1,mtapes

      if(rgnht_int(t)==1) rgnht(t)=.true.


      call strip_null(nfpath(t))
      call strip_null(cpath(t))
      call strip_null(hrestpath(t))
      allocate(tape(t)%hlist(nflds(t)))

      do f=1,nflds(t)
        if (associated(tape(t)%hlist(f)%field%mdims)) then
          deallocate(tape(t)%hlist(f)%field%mdims)
        end if
        nullify(tape(t)%hlist(f)%field%mdims)
        ierr = pio_get_var(File,fillval_desc, (/f,t/), tape(t)%hlist(f)%field%fillvalue)
        ierr = pio_get_var(File,meridional_complement_desc, (/f,t/), tape(t)%hlist(f)%field%meridional_complement)
        ierr = pio_get_var(File,zonal_complement_desc, (/f,t/), tape(t)%hlist(f)%field%zonal_complement)
        ierr = pio_get_var(File,avgflag_desc, (/f,t/), tape(t)%hlist(f)%avgflag)
        ierr = pio_get_var(File,longname_desc, (/1,f,t/), tape(t)%hlist(f)%field%long_name)
        ierr = pio_get_var(File,units_desc, (/1,f,t/), tape(t)%hlist(f)%field%units)
        tape(t)%hlist(f)%field%sampling_seq(1:max_chars) = ' '
        ierr = pio_get_var(File,sseq_desc, (/1,f,t/), tape(t)%hlist(f)%field%sampling_seq)
        call strip_null(tape(t)%hlist(f)%field%sampling_seq)
        if(xyfill(f,t) ==1) then
          tape(t)%hlist(f)%field%flag_xyfill=.true.
        else
          tape(t)%hlist(f)%field%flag_xyfill=.false.
        end if
        if(is_subcol(f,t) ==1) then
           tape(t)%hlist(f)%field%is_subcol=.true.
        else
           tape(t)%hlist(f)%field%is_subcol=.false.
        end if
        call strip_null(tmpname(f,t))
        tape(t)%hlist(f)%field%name = tmpname(f,t)
        tape(t)%hlist(f)%field%decomp_type = decomp(f,t)
        tape(t)%hlist(f)%field%numlev = tmpnumlev(f,t)
        tape(t)%hlist(f)%hwrt_prec = tmpprec(f,t)

        mdimcnt = count(allmdims(:,f,t) > 0)
        if(mdimcnt > 0) then
          allocate(tape(t)%hlist(f)%field%mdims(mdimcnt))
          do i = 1, mdimcnt
            tape(t)%hlist(f)%field%mdims(i) = get_hist_coord_index(mdimnames(allmdims(i,f,t)))
          end do
        end if

      end do
    end do
    deallocate(tmpname, tmpnumlev, tmpprec, decomp, xyfill, is_subcol)
    deallocate(mdimnames)

    allocate(gridsontape(cam_grid_num_grids() + 1, ptapes))
    gridsontape = -1
    do t = 1, ptapes
      do f = 1, nflds(t)
        call set_field_dimensions(tape(t)%hlist(f)%field)

        begdim1 = tape(t)%hlist(f)%field%begdim1
        enddim1 = tape(t)%hlist(f)%field%enddim1
        begdim2 = tape(t)%hlist(f)%field%begdim2
        enddim2 = tape(t)%hlist(f)%field%enddim2
        begdim3 = tape(t)%hlist(f)%field%begdim3
        enddim3 = tape(t)%hlist(f)%field%enddim3

        allocate(tape(t)%hlist(f)%hbuf(begdim1:enddim1,begdim2:enddim2,begdim3:enddim3))

        if (associated(tape(t)%hlist(f)%varid)) then
          deallocate(tape(t)%hlist(f)%varid)
        end if
        nullify(tape(t)%hlist(f)%varid)
        if (associated(tape(t)%hlist(f)%nacs)) then
          deallocate(tape(t)%hlist(f)%nacs)
        end if
        nullify(tape(t)%hlist(f)%nacs)
        if(tape(t)%hlist(f)%field%flag_xyfill .or. (avgflag_pertape(t)=='L')) then
          allocate (tape(t)%hlist(f)%nacs(begdim1:enddim1,begdim3:enddim3))
        else
          allocate(tape(t)%hlist(f)%nacs(1,begdim3:enddim3))
        end if
        ! initialize all buffers to zero - this will be overwritten later by the
        ! data in the history restart file if it exists.
        call h_zero(f,t)

        ! Make sure this field's decomp is listed on the tape
        fdecomp = tape(t)%hlist(f)%field%decomp_type
        do ff = 1, size(gridsontape, 1)
          if (fdecomp == gridsontape(ff, t)) then
            exit
          else if (gridsontape(ff, t) < 0) then
            gridsontape(ff, t) = fdecomp
            exit
          end if
        end do

      end do
    end do
    !
    !-----------------------------------------------------------------------
    ! Read history restart files
    !-----------------------------------------------------------------------
    !
    ! Loop over the total number of history files declared and
    ! read the pathname for any history restart files
    ! that are present (if any). Test to see if the run is a restart run
    ! AND if any history buffer regen files exist (rgnht=.T.). Note, rgnht 
    ! is preset to false, reset to true in routine WSDS if hbuf restart files
    ! are written and saved in the master restart file. Each history buffer
    ! restart file is then obtained.
    ! Note: some f90 compilers (e.g. SGI) complain about I/O of 
    ! derived types which have pointer components, so explicitly read each one.
    ! 
    do t=1,mtapes
      if (rgnht(t)) then
        !
        ! Open history restart file
        !
        call getfil (hrestpath(t), locfn)
        call cam_pio_openfile(tape(t)%File, locfn, 0)
        !
        ! Read history restart file
        !
        do f = 1, nflds(t)  

          fname_tmp = strip_suffix(tape(t)%hlist(f)%field%name)
          if(masterproc) write(iulog, *) 'Reading history variable ',fname_tmp
          ierr = pio_inq_varid(tape(t)%File, fname_tmp, vdesc)

          call cam_pio_var_info(tape(t)%File, vdesc, ndims, dimids, dimlens)
          if(.not. associated(tape(t)%hlist(f)%field%mdims)) then
            dimcnt = 0
            do i=1,ndims
              ierr = pio_inq_dimname(tape(t)%File, dimids(i), dname_tmp)
              dimid = get_hist_coord_index(dname_tmp)
              if(dimid >= 1) then
                dimcnt = dimcnt + 1
                tmpdims(dimcnt) = dimid
              ! No else, just looking for mdims (grid dims won't be hist coords)
              end if
            end do
            if(dimcnt > 0) then
              allocate(tape(t)%hlist(f)%field%mdims(dimcnt))
              tape(t)%hlist(f)%field%mdims(:) = tmpdims(1:dimcnt)
              if(dimcnt > maxvarmdims) maxvarmdims=dimcnt
            end if
          end if
          call set_field_dimensions(tape(t)%hlist(f)%field)
          begdim1    =  tape(t)%hlist(f)%field%begdim1
          enddim1    =  tape(t)%hlist(f)%field%enddim1
          fdims(1)   =  enddim1 - begdim1 + 1
          begdim2    =  tape(t)%hlist(f)%field%begdim2
          enddim2    =  tape(t)%hlist(f)%field%enddim2
          fdims(2)   =  enddim2 - begdim2 + 1
          begdim3    =  tape(t)%hlist(f)%field%begdim3
          enddim3    =  tape(t)%hlist(f)%field%enddim3
          fdims(3)   =  enddim3 - begdim3 + 1
          if (fdims(2) > 1) then
            nfdims = 3
          else
            nfdims = 2
            fdims(2) = fdims(3)
          end if
          fdecomp = tape(t)%hlist(f)%field%decomp_type
          if (nfdims > 2) then
            call cam_grid_read_dist_array(tape(t)%File, fdecomp,              &
                 fdims(1:nfdims), dimlens(1:ndims), tape(t)%hlist(f)%hbuf, vdesc)
          else
            call cam_grid_read_dist_array(tape(t)%File, fdecomp,              &
                 fdims(1:nfdims), dimlens(1:ndims), tape(t)%hlist(f)%hbuf(:,1,:), vdesc)
          end if

          ierr = pio_inq_varid(tape(t)%File, trim(fname_tmp)//'_nacs', vdesc)
          call cam_pio_var_info(tape(t)%File, vdesc, nacsdimcnt, dimids, dimlens)

          if(nacsdimcnt > 0) then
            if (nfdims > 2) then
              ! nacs only has 2 dims (no levels)
              fdims(2) = fdims(3)
            end if
            allocate(tape(t)%hlist(f)%nacs(begdim1:enddim1,begdim3:enddim3))
            nacs       => tape(t)%hlist(f)%nacs(:,:)
            call cam_grid_read_dist_array(tape(t)%File, fdecomp, fdims(1:2),  &
                 dimlens(1:nacsdimcnt), nacs, vdesc)
          else
            allocate(tape(t)%hlist(f)%nacs(1,begdim3:enddim3))
            ierr = pio_get_var(tape(t)%File, vdesc, nacsval)
            tape(t)%hlist(f)%nacs(1,:)= nacsval
          end if

        end do
        !          
        ! Done reading this history restart file
        !
        call cam_pio_closefile(tape(t)%File)

      end if  ! rgnht(t)

      ! (re)create the master list of grid IDs
      ff = 0
      do f = 1, size(gridsontape, 1)
        if (gridsontape(f, t) > 0) then
          ff = ff + 1
        end if
      end do
      allocate(tape(t)%grid_ids(ff))
      ff = 1
      do f = 1, size(gridsontape, 1)
        if (gridsontape(f, t) > 0) then
          tape(t)%grid_ids(ff) = gridsontape(f, t)
          ff = ff + 1
        end if
      end do
      call patch_init(t)
    end do     ! end of do mtapes loop

    !
    ! If the history files are partially complete (contain less than
    ! mfilt(t) time samples, then get the files and open them.)
    !
    ! NOTE:  No need to perform this operation for IC history files or empty files
    !

    do t=1,mtapes
      if (is_initfile(file_index=t)) then
        ! Initialize filename specifier for IC file
        hfilename_spec(t) = '%c.cam' // trim(inst_suffix) // '.i.%y-%m-%d-%s.nc'
        nfils(t) = 0
      else if (nflds(t) == 0) then
        nfils(t) = 0
      else
        if (nfils(t) > 0) then
          call getfil (cpath(t), locfn)
          call cam_pio_openfile(tape(t)%File, locfn, PIO_WRITE)
          call h_inquire (t)
          if(is_satfile(t)) then
            !  Initialize the sat following history subsystem
            call sat_hist_init()
            call sat_hist_define(tape(t)%File)
          end if
        end if
        !
        ! If the history file is full, close the current unit
        !
        if (nfils(t) >= mfilt(t)) then
          if (masterproc) then
            write(iulog,*)'READ_RESTART_HISTORY: nf_close(',t,')=',nhfil(t), mfilt(t)
          end if
          do f=1,nflds(t)
            deallocate(tape(t)%hlist(f)%varid)
            nullify(tape(t)%hlist(f)%varid)
          end do
          call cam_pio_closefile(tape(t)%File)
          nfils(t) = 0
        end if
      end if
    end do

    ! Setup vector pairs for unstructured grid interpolation
    call setup_interpolation_and_define_vector_complements()

    if(mtapes/=ptapes .and. masterproc) then
      write(iulog,*) ' WARNING: Restart file ptapes setting ',mtapes,' not equal to model setting ',ptapes
    end if

    return
  end subroutine read_restart_history

  !#######################################################################

  character(len=max_string_len) function get_hfilepath( tape )
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Return full filepath of history file for given tape number
    ! This allows public read access to the filenames without making
    ! the filenames public data.
    !
    !----------------------------------------------------------------------- 
    !
    integer, intent(in) :: tape  ! Tape number

    get_hfilepath = cpath( tape )
  end function get_hfilepath

  !#######################################################################

  character(len=max_string_len) function get_hist_restart_filepath( tape )
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Return full filepath of restart file for given tape number
    ! This allows public read access to the filenames without making
    ! the filenames public data.
    !
    !----------------------------------------------------------------------- 
    !
    integer, intent(in) :: tape  ! Tape number

    get_hist_restart_filepath = hrestpath( tape )
  end function get_hist_restart_filepath

  !#######################################################################

  integer function get_ptapes( )
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Return the number of tapes being used.
    ! This allows public read access to the number of tapes without making
    ! ptapes public data.
    !
    !----------------------------------------------------------------------- 
    !
    get_ptapes = ptapes
  end function get_ptapes

  !#######################################################################

  recursive function get_entry_by_name(listentry, name) result(entry)
    type(master_entry),  pointer :: listentry
    character(len=*), intent(in) :: name ! variable name
    type(master_entry), pointer :: entry

    if(associated(listentry)) then
      if(listentry%field%name .eq. name) then
        entry => listentry
      else
        entry=>get_entry_by_name(listentry%next_entry, name)
      end if
    else
      nullify(entry)
    end if
  end function get_entry_by_name

  !#######################################################################

  subroutine AvgflagToString(avgflag, time_op)
    ! Dummy arguments
    character(len=1),           intent(in)  :: avgflag ! averaging flag
    character(len=max_chars),   intent(out) :: time_op ! time op (e.g. max)

    ! Local variable
    character(len=*), parameter             :: subname = 'AvgflagToString'

    select case (avgflag)
    case ('A')
      time_op(:) = 'mean'
    case ('B')
      time_op(:) = 'mean00z'
    case ('I')
      time_op(:) = ' '
    case ('X')
      time_op(:) = 'maximum'
    case ('M')
      time_op(:) = 'minimum'
    case('L')
      time_op(:) = LT_DESC
    case default
      call endrun(subname//': unknown avgflag = '//avgflag)
    end select
  end subroutine AvgflagToString

  !#######################################################################

  subroutine fldlst ()

    use cam_grid_support, only: cam_grid_num_grids
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Define the contents of each history file based on namelist input for initial or branch
    ! run, and restart data if a restart run.
    !          
    ! Method: Use arrays fincl and fexcl to modify default history tape contents.
    !         Then sort the result alphanumerically for later use by OUTFLD to
    !         allow an n log n search time.
    !
    !---------------------------Local variables-----------------------------
    !
    integer t, f                   ! tape, field indices
    integer ff                     ! index into include, exclude and fprec list
    integer num_patches            ! number of column areas
    character(len=fieldname_len) :: name ! field name portion of fincl (i.e. no avgflag separator)
    character(len=max_fieldname_len) :: mastername ! name from masterlist field
    character(len=max_chars) :: errormsg ! error output field
    character(len=1) :: avgflag    ! averaging flag
    character(len=1) :: prec_wrt   ! history buffer write precision flag

    type (hentry) :: tmp           ! temporary used for swapping

    type(master_entry), pointer :: listentry
    logical                     :: fieldontape      ! .true. iff field on tape

    ! List of active grids (first dim) for each tape (second dim)
    ! An active grid is one for which there is a least one field being output
    !    on that grid.
    integer, allocatable        :: gridsontape(:,:)

    !
    ! First ensure contents of fincl, fexcl, and fwrtpr are all valid names
    !
    do t=1,ptapes
      f = 1
      do while (f < pflds .and. fincl(f,t) /= ' ')
        name = getname (fincl(f,t))
        mastername=''
        listentry => get_entry_by_name(masterlinkedlist, name)
        if(associated(listentry)) mastername = listentry%field%name
        if (name /= mastername) then
          write(iulog,*)'FLDLST: ', trim(name), ' in fincl(', f, ') not found'
          call endrun
        end if
        f = f + 1
      end do

      f = 1
      do while (f < pflds .and. fexcl(f,t) /= ' ')
        mastername=''
        listentry => get_entry_by_name(masterlinkedlist, fexcl(f,t))
        if(associated(listentry)) mastername = listentry%field%name

        if (fexcl(f,t) /= mastername) then
          write(iulog,*)'FLDLST: ', fexcl(f,t), ' in fexcl(', f, ') not found'
          call endrun
        end if
        f = f + 1
      end do

      f = 1
      do while (f < pflds .and. fwrtpr(f,t) /= ' ')
        name = getname (fwrtpr(f,t))
        mastername=''
        listentry => get_entry_by_name(masterlinkedlist, name)
        if(associated(listentry)) mastername = listentry%field%name
        if (name /= mastername) then
          write(iulog,*)'FLDLST: ', trim(name), ' in fwrtpr(', f, ') not found'
          call endrun
        end if
        do ff=1,f-1                 ! If duplicate entry is found, stop
          if (trim(name) == trim(getname(fwrtpr(ff,t)))) then
            write(iulog,*)'FLDLST: Duplicate field ', name, ' in fwrtpr'
            call endrun
          end if
        end do
        f = f + 1
      end do
    end do

    nflds(:) = 0
    ! IC history file is to be created, set properties
    if(is_initfile()) then
      hfilename_spec(ptapes) = '%c.cam' // trim(inst_suffix) // '.i.%y-%m-%d-%s.nc'

      ncprec(ptapes) = pio_double
      ndens (ptapes) = 1
      mfilt (ptapes) = 1
    end if



    allocate(gridsontape(cam_grid_num_grids() + 1, ptapes))
    gridsontape = -1
    do t=1,ptapes
      !
      ! Add the field to the tape if specified via namelist (FINCL[1-ptapes]), or if
      ! it is on by default and was not excluded via namelist (FEXCL[1-ptapes]).
      ! Also set history buffer accumulation and output precision values according
      ! to the values specified via namelist (FWRTPR[1-ptapes])
      ! or, if not on the list, to the default values given by ndens(t).
      !
      listentry => masterlinkedlist
      do while(associated(listentry))
        mastername = listentry%field%name
        call list_index (fincl(1,t), mastername, ff)

        fieldontape = .false.
        if (ff > 0) then
          fieldontape = .true.
        else if ((.not. empty_htapes) .or. (is_initfile(file_index=t))) then
          call list_index (fexcl(1,t), mastername, ff)
          if (ff == 0 .and. listentry%actflag(t)) then
            fieldontape = .true.
          end if
        end if
        if (fieldontape) then
          ! The field is active so increment the number fo fields and add
          ! its decomp type to the list of decomp types on this tape
          nflds(t) = nflds(t) + 1
          do ff = 1, size(gridsontape, 1)
            if (listentry%field%decomp_type == gridsontape(ff, t)) then
              exit
            else if (gridsontape(ff, t) < 0) then
              gridsontape(ff, t) = listentry%field%decomp_type
              exit
            end if
          end do
        end if
        listentry=>listentry%next_entry
      end do
    end do
    !
    ! Determine total number of active history tapes
    !
    if (masterproc) then
      do t=1,ptapes
        if (nflds(t)  ==  0) then
          write(iulog,*)'FLDLST: Tape ',t,' is empty'
        end if
      end do
    endif
    allocate(history_tape(ptapes))
    tape=>history_tape


    do t=1,ptapes
      nullify(tape(t)%hlist)
      ! Now we have a field count and can allocate
      if(nflds(t) > 0) then
        ! Allocate the correct number of hentry slots
        allocate(tape(t)%hlist(nflds(t)))
        ! Count up the number of grids output on this tape
        ff = 0
        do f = 1, size(gridsontape, 1)
          if (gridsontape(f, t) > 0) then
            ff = ff + 1
          end if
        end do
        allocate(tape(t)%grid_ids(ff))
        ff = 1
        do f = 1, size(gridsontape, 1)
          if (gridsontape(f, t) > 0) then
            tape(t)%grid_ids(ff) = gridsontape(f, t)
            ff = ff + 1
          end if
        end do
      end if
      do ff=1,nflds(t)
        nullify(tape(t)%hlist(ff)%hbuf)
        nullify(tape(t)%hlist(ff)%nacs)
        nullify(tape(t)%hlist(ff)%varid)
      end do


      nflds(t) = 0 ! recount to support array based method
      listentry => masterlinkedlist
      do while(associated(listentry))
        mastername = listentry%field%name

        call list_index (fwrtpr(1,t), mastername, ff)
        if (ff > 0) then
          prec_wrt = getflag(fwrtpr(ff,t))
        else
          prec_wrt = ' '
        end if

        call list_index (fincl(1,t), mastername, ff)

        if (ff > 0) then
          avgflag = getflag (fincl(ff,t))
          call inifld (t, listentry, avgflag,  prec_wrt)
        else if ((.not. empty_htapes) .or. (is_initfile(file_index=t))) then
          call list_index (fexcl(1,t), mastername, ff)
          if (ff == 0 .and. listentry%actflag(t)) then
            call inifld (t, listentry, ' ', prec_wrt)
          else
            listentry%actflag(t) = .false.
          end if
        else
          listentry%actflag(t) = .false.
        end if
        listentry=>listentry%next_entry

      end do
      !
      ! If column output is specified make sure there are some fields defined
      ! for that tape
      !
      if (nflds(t) .eq. 0 .and. fincllonlat(1,t) .ne. ' ') then
        write(errormsg,'(a,i2,a)') 'FLDLST: Column output is specified for tape ',t,' but no fields defined for that tape.'
        call endrun(errormsg)
      else
        call patch_init(t)
      end if
      !
      ! Specification of tape contents now complete.  Sort each list of active 
      ! entries for efficiency in OUTFLD.  Simple bubble sort.
      !
!!XXgoldyXX: v In the future, we will sort according to decomp to speed I/O
      do f=nflds(t)-1,1,-1
        do ff=1,f

          if (tape(t)%hlist(ff)%field%name > tape(t)%hlist(ff+1)%field%name) then

            tmp = tape(t)%hlist(ff)
            tape(t)%hlist(ff  ) = tape(t)%hlist(ff+1)
            tape(t)%hlist(ff+1) = tmp

          else if (tape(t)%hlist(ff  )%field%name == tape(t)%hlist(ff+1)%field%name) then

            write(errormsg,'(2a,2(a,i3))') 'FLDLST: Duplicate field: ', &
                 trim(tape(t)%hlist(ff)%field%name),', tape = ', t, ', ff = ', ff
            call endrun(errormsg)

          end if

        end do
      end do

    end do    ! do t=1,ptapes
    deallocate(gridsontape)

    call print_active_fldlst()

    !
    ! Packing density, ndens: With netcdf, only 1 (nf_double) and 2 (pio_real)
    ! are allowed
    !
    do t=1,ptapes
      if (ndens(t) == 1) then
        ncprec(t) = pio_double
      else if (ndens(t) == 2) then
        ncprec(t) = pio_real
      else
        call endrun ('FLDLST: ndens must be 1 or 2')
      end if

    end do
    !
    !  Now that masterlinkedlist is defined, construct primary and secondary hashing
    !  tables.
    !
    call bld_outfld_hash_tbls()
    call bld_htapefld_indices()

    return
  end subroutine fldlst

!#########################################################################################

subroutine print_active_fldlst()

   integer :: f, ff, i, t
   integer :: num_patches

   character(len=6) :: prec_str
   character(len=max_chars) :: fldname, fname_tmp

   type(active_entry), pointer :: hfile(:) => null()  ! history files

   if (masterproc) then

      hfile=>history_tape

      do t=1,ptapes

         if (nflds(t) > 0) then
            write(iulog,*) ' '
            write(iulog,*)'FLDLST: History file ', t, ' contains ', nflds(t), ' fields'

            if (is_initfile(file_index=t)) then
               write(iulog,*) ' Write frequency:                 ',inithist,' (INITIAL CONDITIONS)'
            else
               if (nhtfrq(t) == 0) then
                  write(iulog,*) ' Write frequency:                  MONTHLY'
               else
                  write(iulog,*) ' Write frequency:                 ',nhtfrq(t)
               end if
            end if

            write(iulog,*) ' Filename specifier:              ', trim(hfilename_spec(t))

            prec_str = 'double'
            if (ndens(t) == 2) prec_str = 'single'
            write(iulog,*) ' Output precision:                ', prec_str
            write(iulog,*) ' Number of time samples per file: ', mfilt(t)

            ! grid info
            if (associated(hfile(t)%patches)) then
               write(iulog,*) ' Fields are represented on columns (FIELD_LON_LAT)'
            else if (associated(hfile(t)%grid_ids)) then
               write(iulog,*) ' Fields are represented on global grids:'
               do i = 1, size(hfile(t)%grid_ids)
                  write(iulog,*) ' ', hfile(t)%grid_ids(i)
               end do
            else
               call endrun('print_active_fldlst: error in active_entry object')
            end if

            write(iulog,*)' Included fields are:'

         end if

         do f = 1, nflds(t)
            if (associated(hfile(t)%patches)) then
               num_patches = size(hfile(t)%patches)
               fldname = strip_suffix(hfile(t)%hlist(f)%field%name)
               do i = 1, num_patches
                  ff = (f-1)*num_patches + i
                  fname_tmp = trim(fldname)
                  call hfile(t)%patches(i)%field_name(fname_tmp)
                  write(iulog,9000) ff, fname_tmp, hfile(t)%hlist(f)%field%units, &
                     hfile(t)%hlist(f)%field%numlev, hfile(t)%hlist(f)%avgflag,   &
                     trim(hfile(t)%hlist(f)%field%long_name)
               end do
            else
               fldname = hfile(t)%hlist(f)%field%name
               write(iulog,9000) f, fldname, hfile(t)%hlist(f)%field%units,  &
                  hfile(t)%hlist(f)%field%numlev, hfile(t)%hlist(f)%avgflag, &
                  trim(hfile(t)%hlist(f)%field%long_name)
            end if

         end do

      end do

   end if

9000 format(i5, 1x, a32, 1x, a16, 1x, i4, 1x, a1, 2x, 256a)

end subroutine print_active_fldlst

!#########################################################################################

  subroutine inifld (t, listentry, avgflag, prec_wrt)
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Add a field to the active list for a history tape
    ! 
    ! Method: Copy the data from the master field list to the active list for the tape
    !         Also: define mapping arrays from (col,chunk) -> (lon,lat)
    ! 
    ! Author: CCM Core Group
    ! 
    !-----------------------------------------------------------------------


    !
    ! Arguments
    !
    integer, intent(in) :: t            ! history tape index

    type(master_entry), pointer :: listentry

    character*1, intent(in) :: avgflag  ! averaging flag
    character*1, intent(in) :: prec_wrt ! history output precision flag
    !
    ! Local workspace
    !
    integer :: n                  ! field index on defined tape


    !
    ! Ensure that it is not to late to add a field to the history tape
    !
    if (htapes_defined) then
      call endrun ('INIFLD: Attempt to add field '//listentry%field%name//' after history files set')
    end if


    nflds(t) = nflds(t) + 1
    n = nflds(t)
    !
    ! Copy field info.
    !
    if(n > size(tape(t)%hlist)) then
      write(iulog,*) 'tape field miscount error ', n, size(tape(t)%hlist)
      call endrun()
    end if

    tape(t)%hlist(n)%field = listentry%field

    select case (prec_wrt)
    case (' ')
      if (ndens(t) == 1) then
        tape(t)%hlist(n)%hwrt_prec = 8
      else
        tape(t)%hlist(n)%hwrt_prec = 4
      end if
    case ('4')
      tape(t)%hlist(n)%hwrt_prec = 4
      if (masterproc) then
        write(iulog,*) 'INIFLD: Output data type for ', tape(t)%hlist(n)%field%name, &
             ' is real*4'
      end if
    case ('8')
      tape(t)%hlist(n)%hwrt_prec = 8
      if (masterproc) then
        write(iulog,*) 'INIFLD: Output data type for ', tape(t)%hlist(n)%field%name, &
             ' is real*8'
      end if
    case default
      call endrun ('INIFLD: unknown prec_wrt='//prec_wrt)
    end select
    !
    ! Override the default averaging (masterlist) averaging flag if non-blank
    !
    if (avgflag == ' ') then
      tape(t)%hlist(n)%avgflag = listentry%avgflag(t)
      tape(t)%hlist(n)%time_op = listentry%time_op(t)
    else
      tape(t)%hlist(n)%avgflag = avgflag
      call AvgflagToString(avgflag, tape(t)%hlist(n)%time_op)
    end if

#ifdef HDEBUG
    write(iulog,*)'HDEBUG: ',__LINE__,' field ', tape(t)%hlist(n)%field%name, ' added as ', 'field number ', n,' on tape ', t
    write(iulog,*)'units=',tape(t)%hlist(n)%field%units
    write(iulog,*)'numlev=',tape(t)%hlist(n)%field%numlev
    write(iulog,*)'avgflag=',tape(t)%hlist(n)%avgflag
    write(iulog,*)'time_op=',tape(t)%hlist(n)%time_op
    write(iulog,*)'hwrt_prec=',tape(t)%hlist(n)%hwrt_prec
#endif

    return
  end subroutine inifld


  subroutine patch_init(t)
    use cam_history_support, only: history_patch_t
    use cam_grid_support,    only: cam_grid_compute_patch

    ! Dummy arguments
    integer, intent(in)               :: t     ! Current tape

    ! Local variables
    integer                           :: ff    ! Loop over fincllonlat entries
    integer                           :: i     ! General loop index
    integer                           :: npatches
    type(history_patch_t), pointer    :: patchptr
    
    character(len=max_chars)          :: errormsg
    character(len=max_chars)          :: lonlatname(pflds)
    real(r8)                          :: beglon, beglat, endlon, endlat

    !
    ! Setup column information if this field will be written as group
    ! First verify the column information in the namelist
    ! Duplicates are an error, but we can just ignore them
    !

    ! I know, this shouldn't happen . . . yet: (better safe than sorry)
    if (associated(tape(t)%patches)) then
      do i = 1, size(tape(t)%patches)
        call tape(t)%patches(i)%deallocate()
      end do
      deallocate(tape(t)%patches)
      nullify(tape(t)%patches)
    end if

    ! First, count the number of patches and check for duplicates
    ff = 1  ! Index of fincllonlat entry
    npatches = 0   ! Number of unique patches in namelist entry
    do while (len_trim(fincllonlat(ff, t)) > 0)
      npatches = npatches + 1
      lonlatname(npatches) = trim(fincllonlat(ff, t))
      ! Check for duplicates
      do i = 1, npatches - 1
        if (trim(lonlatname(i)) == trim(lonlatname(npatches))) then
          write(errormsg, '(a,i0,3a)') 'Duplicate fincl', t, 'lonlat entry.', &
               'Duplicate entry is ', trim(lonlatname(i))
          write(iulog, *) 'patch_init: WARNING: '//errormsg
          ! Remove the new entry
          lonlatname(npatches) = ''
          npatches = npatches - 1
          exit
        end if
      end do
      ff = ff + 1
    end do

    ! Now we know how many patches, allocate space
    if (npatches > 0) then
      if (collect_column_output(t)) then
        allocate(tape(t)%patches(1))
      else
        allocate(tape(t)%patches(npatches))
      end if

      ! For each lat/lon specification, parse and create a patch for each grid
      do ff = 1, npatches
        if (collect_column_output(t)) then
          ! For colleccted column output, we only have one patch
          patchptr => tape(t)%patches(1)
        else
          patchptr => tape(t)%patches(ff)
          patchptr%namelist_entry = trim(lonlatname(ff))
        end if
        ! We need to set up one patch per (active) grid
        patchptr%collected_output = collect_column_output(t)
        call parseLonLat(lonlatname(ff),                                      &
             beglon, endlon, patchptr%lon_axis_name,                          &
             beglat, endlat, patchptr%lat_axis_name)
        if (associated(patchptr%patches)) then
          ! One last sanity check
          if (.not. collect_column_output(t)) then
            write(errormsg, '(a,i0,2a)') 'Attempt to overwrite fincl', t,     &
                 'lonlat entry, ', trim(patchptr%namelist_entry)
            call endrun('patch_init: '//errormsg)
          end if
        else
          allocate(patchptr%patches(size(tape(t)%grid_ids)))
        end if
        do i = 1, size(tape(t)%grid_ids)
          call cam_grid_compute_patch(tape(t)%grid_ids(i), patchptr%patches(i),&
               beglon, endlon, beglat, endlat)
        end do
        nullify(patchptr)
      end do
    end if
    ! We are done processing this tape's fincl#lonlat entries. Now, 
    ! compact each patch so that the output variables have no holes
    ! We wait until now for when collect_column_output(t) is .true. since
    !    all the fincl#lonlat entries are concatenated
    if (associated(tape(t)%patches)) then
      do ff = 1, size(tape(t)%patches)
        call tape(t)%patches(ff)%compact()
      end do
    end if

  end subroutine patch_init

  !#######################################################################

  subroutine strip_null(str)
    character(len=*), intent(inout) :: str
    do i=1,len(str)
      if(ichar(str(i:i))==0) str(i:i)=' '
    end do
  end subroutine strip_null

  character(len=max_fieldname_len) function strip_suffix (name)
    !
    !---------------------------------------------------------- 
    ! 
    ! Purpose:  Strip "&IC" suffix from fieldnames if it exists
    !          
    !----------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*), intent(in) :: name
    !
    ! Local workspace
    !
    integer :: n
    !
    !-----------------------------------------------------------------------
    !
    strip_suffix = ' '

    do n = 1,fieldname_len
      strip_suffix(n:n) = name(n:n)
      if(name(n+1:n+1         ) == ' '                       ) return
      if(name(n+1:n+fieldname_suffix_len) == fieldname_suffix) return
    end do

    strip_suffix(fieldname_len+1:max_fieldname_len) = name(fieldname_len+1:max_fieldname_len)

    return

  end function strip_suffix

  !#######################################################################

  character(len=fieldname_len) function getname (inname)
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: retrieve name portion of inname
    !          
    ! Method:  If an averaging flag separater character is present (":") in inname, 
    !          lop it off
    ! 
    !-------------------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*), intent(in) :: inname
    !
    ! Local workspace
    !
    integer :: length
    integer :: i

    length = len (inname)

    if (length < fieldname_len .or. length > fieldname_lenp2) then
      write(iulog,*) 'GETNAME: bad length=',length
      call endrun
    end if

    getname = ' '
    do i=1,fieldname_len
      if (inname(i:i) == ':') exit
      getname(i:i) = inname(i:i)
    end do

    return
  end function getname

  !#######################################################################

  ! parseRangeString: Parse either a coordinate descriptor (e.g., 10S) or a
  !                   coordinate range (e.g., 10e:20e)
  !                   chars represents the allowed coordinate character.
  !                   NB: Does not validate numerical values (e.g., lat <= 90)
  subroutine parseRangeString(rangestr, chars, begval, begchar, begname, endval, endchar, endname)

    ! Dummy arguments
    character(len=*),       intent(in)    :: rangestr
    character(len=*),       intent(in)    :: chars
    real(r8),               intent(out)   :: begval
    character,              intent(out)   :: begchar
    character(len=*),       intent(out)   :: begname
    real(r8),               intent(out)   :: endval
    character,              intent(out)   :: endchar
    character(len=*),       intent(out)   :: endname

    ! Local variables
    character(len=128)                    :: errormsg
    integer                               :: colonpos
    integer                               :: beglen, endlen

    ! First, see if we have a position or a range
    colonpos = scan(rangestr, ':')
    if (colonpos == 0) then
      begname = trim(rangestr)
      beglen = len_trim(begname)
      endname = trim(begname)
    else
      beglen = colonpos - 1
      begname = rangestr(1:beglen)
      endname = trim(rangestr(colonpos+1:))
      endlen = len_trim(endname)
    end if
    ! begname should be a number (integer or real) followed by a character
    if (verify(begname, '0123456789.') /= beglen) then
      write(errormsg, *) 'Coordinate range must begin with number, ', begname
      call endrun('parseRangeString: '//errormsg)
    end if
    if (verify(begname(beglen:beglen), chars) /= 0) then
      write(errormsg, *) 'Coordinate range must end with character in the ',  &
           'set [', trim(chars), '] ', begname
      call endrun('parseRangeString: '//errormsg)
    end if
    ! begname parses so collect the values
    read(begname(1:beglen-1), *) begval
    begchar = begname(beglen:beglen)
    if (colonpos /= 0) then
      ! endname should be a number (integer or real) followed by a character
      if (verify(endname, '0123456789.') /= endlen) then
        write(errormsg, *) 'Coordinate range must begin with number, ', endname
        call endrun('parseRangeString: '//errormsg)
      end if
      if (verify(endname(endlen:endlen), chars) /= 0) then
        write(errormsg, *) 'Coordinate range must end with character in the ',&
             'set [', trim(chars), '] ', endname
        call endrun('parseRangeString: '//errormsg)
      end if
      ! endname parses so collect the values
      read(endname(1:endlen-1), *) endval
      endchar = endname(endlen:endlen)
    else
      endval = begval
      endchar = begchar
    end if

  end subroutine parseRangeString

  ! parseLonLat: Parse a lon_lat description allowed by the fincllonlat(n)
  !              namelist entries. Returns the starting and ending values of
  !              the point or range specified.
  !              NB: Does not validate the range against any particular grid
  subroutine parseLonLat(lonlatname, beglon, endlon, lonname, beglat, endlat, latname)

    ! Dummy arguments
    character(len=*),       intent(in)    :: lonlatname
    real(r8),               intent(out)   :: beglon
    real(r8),               intent(out)   :: endlon
    character(len=*),       intent(out)   :: lonname
    real(r8),               intent(out)   :: beglat
    real(r8),               intent(out)   :: endlat
    character(len=*),       intent(out)   :: latname

    ! Local variables
    character(len=128)                    :: errormsg
    character(len=MAX_CHARS)              :: lonstr, latstr
    character(len=MAX_CHARS)              :: begname, endname
    character                             :: begchar, endchar
    integer                               :: underpos

    !
    ! make sure _ separator is present
    !      
    underpos = scan(lonlatname, '_')
    if (underpos == 0) then
      write(errormsg,*) 'Improperly formatted fincllonlat string. ',          &
           'Missing underscore character (xxxE_yyyS) ', lonlatname
      call endrun('parseLonLat: '//errormsg)
    end if

    ! Break out the longitude and latitude sections
    lonstr = lonlatname(:underpos-1)
    latstr = trim(lonlatname(underpos+1:))

    ! Parse the longitude section
    call parseRangeString(lonstr, 'eEwW', beglon, begchar, begname, endlon, endchar, endname)
    ! Convert longitude to degrees East
    if ((begchar == 'w') .or. (begchar == 'W')) then
      beglon = 360._r8 - beglon
    end if
    if ((beglon < 0._r8) .or. (beglon > 360._r8)) then
      write(errormsg, *) 'Longitude specification out of range, ', trim(begname)
      call endrun('parseLonLat: '//errormsg)
    end if
    if ((endchar == 'w') .or. (endchar == 'W')) then
      endlon = 360._r8 - endlon
    end if
    if ((endlon < 0._r8) .or. (endlon > 360._r8)) then
      write(errormsg, *) 'Longitude specification out of range, ', trim(endname)
      call endrun('parseLonLat: '//errormsg)
    end if
    if (beglon == endlon) then
      lonname = trim(begname)
    else
      lonname = trim(begname)//'_to_'//trim(endname)
    end if

    ! Parse the latitude section
    call parseRangeString(latstr, 'nNsS', beglat, begchar, begname, endlat, endchar, endname)
    ! Convert longitude to degrees East
    if ((begchar == 's') .or. (begchar == 'S')) then
      beglat = (-1._r8) * beglat
    end if
    if ((beglat < -90._r8) .or. (beglat > 90._r8)) then
      write(errormsg, *) 'Latitude specification out of range, ', trim(begname)
      call endrun('parseLonLat: '//errormsg)
    end if
    if ((endchar == 's') .or. (endchar == 'S')) then
      endlat = (-1._r8) * endlat
    end if
    if ((endlat < -90._r8) .or. (endlat > 90._r8)) then
      write(errormsg, *) 'Latitude specification out of range, ', trim(endname)
      call endrun('parseLonLat: '//errormsg)
    end if
    if (beglat == endlat) then
      latname = trim(begname)
    else
      latname = trim(begname)//'_to_'//trim(endname)
    end if

  end subroutine parseLonLat


  !#######################################################################

  character(len=1) function getflag (inname)
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: retrieve flag portion of inname
    !          
    ! Method:  If an averaging flag separater character is present (":") in inname, 
    !          return the character after it as the flag
    ! 
    !-------------------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*), intent(in) :: inname   ! character string
    !
    ! Local workspace
    !
    integer :: length         ! length of inname
    integer :: i              ! loop index

    length = len (inname)

    if (length /= fieldname_lenp2) then
      write(iulog,*) 'GETFLAG: bad length=',length
      call endrun
    end if

    getflag = ' '
    do i=1,fieldname_lenp2-1
      if (inname(i:i) == ':') then
        getflag = inname(i+1:i+1)
        exit
      end if
    end do

    return
  end function getflag

  !#######################################################################

  subroutine list_index (list, name, index)
    !
    ! Input arguments
    !
    character(len=*), intent(in) :: list(pflds) ! input list of names, possibly ":" delimited
    character(len=max_fieldname_len), intent(in) :: name ! name to be searched for
    !
    ! Output arguments
    !
    integer, intent(out) :: index               ! index of "name" in "list"
    !
    ! Local workspace
    !
    character(len=fieldname_len) :: listname    ! input name with ":" stripped off.
    integer f                       ! field index

    index = 0
    do f=1,pflds
      !
      ! Only list items
      !
      listname = getname (list(f))
      if (listname == ' ') exit
      if (listname == name) then
        index = f
        exit
      end if
    end do

    return
  end subroutine list_index

  !#######################################################################

  recursive subroutine outfld (fname, field, idim, c, avg_subcol_field)
    use cam_history_buffers, only: hbuf_accum_inst, hbuf_accum_add,  &
         hbuf_accum_add00z, hbuf_accum_max, hbuf_accum_min,          &
         hbuf_accum_addlcltime
    use cam_history_support, only: dim_index_2d
    use subcol_utils,        only: subcol_unpack
    use cam_grid_support,    only: cam_grid_id

    interface
      subroutine subcol_field_avg_handler(idim, field_in, c, field_out)
        use shr_kind_mod, only: r8 => shr_kind_r8
        integer,  intent(in)  :: idim
        real(r8), intent(in)  :: field_in(idim, *)
        integer,  intent(in)  :: c
        real(r8), intent(out) :: field_out(:,:)
      end subroutine subcol_field_avg_handler
    end interface

    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Accumulate (or take min, max, etc. as appropriate) input field
    !          into its history buffer for appropriate tapes
    ! 
    ! Method: Check 'masterlist' whether the requested field 'fname' is active
    !         on one or more history tapes, and if so do the accumulation.
    !         If not found, return silently.
    !         subcol_field_avg_handler:
    !            An interface into subcol_field_avg without creating a dependency as
    !            this would cause a dependency loop. See subcol.F90
    ! Note: We cannot know a priori if field is a grid average field or a subcolumn
    !       field because many fields passed to outfld are defined on ncol rather
    !       than pcols or psetcols. Therefore, we use the avg_subcol_field input
    !       to determine whether to average the field input before accumulation.
    !       NB: If output is on a subcolumn grid (requested in addfle), it is
    !           an error to use avg_subcol_field. A subcolumn field is assumed and
    !           subcol_unpack is called before accumulation.
    ! 
    ! Author: CCM Core Group
    ! 
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*), intent(in) :: fname ! Field name--should be 8 chars long

    ! For structured grids, idim is the local longitude dimension.
    ! For unstructured grids, idim is the local column dimension
    ! For phys_decomp, it should be pcols or pcols*psubcols
    integer, intent(in)           :: idim
    real(r8), intent(in)          :: field(idim,*) ! Array containing field values
    integer, intent(in)           :: c             ! chunk (physics) or latitude (dynamics) index
    logical, optional, intent(in) :: avg_subcol_field
    !
    ! Local variables
    !
    integer               :: t, f          ! tape, field indices

    character*1           :: avgflag       ! averaging flag

    type (active_entry), pointer :: otape(:) ! Local history_tape pointer
    real(r8),pointer      :: hbuf(:,:)     ! history buffer
    integer, pointer      :: nacs(:)       ! accumulation counter
    integer               :: begdim2, enddim2, endi
    integer               :: phys_decomp
    type (dim_index_2d)   :: dimind        ! 2-D dimension index
    logical               :: flag_xyfill   ! non-applicable xy points flagged with fillvalue
    real(r8)              :: fillvalue
    real(r8), allocatable :: afield(:,:)   ! Averaged field values
    real(r8), allocatable :: ufield(:,:,:) ! Unpacked field values
    integer               :: ff            ! masterlist index pointer
    integer               :: i, j
    logical               :: found
    logical               :: avg_subcols   ! average subcols before accum
    !-----------------------------------------------------------------------

    call get_field_properties(fname, found, tape_out=otape, ff_out=ff)
    phys_decomp = cam_grid_id('physgrid')

    ! If this field is not active, return now
    if (.not. found) then
      return
    end if

    !
    ! Note, the field may be on any or all of the history files (primary
    ! and auxiliary).
    !
    !      write(iulog,*)'fname_loc=',fname_loc
    do t = 1, ptapes
      if ( .not. masterlist(ff)%thisentry%actflag(t)) cycle
      f = masterlist(ff)%thisentry%htapeindx(t)
      !
      ! Update history buffer
      !
      flag_xyfill = otape(t)%hlist(f)%field%flag_xyfill
      fillvalue = otape(t)%hlist(f)%field%fillvalue
      avgflag = otape(t)%hlist(f)%avgflag
      nacs   => otape(t)%hlist(f)%nacs(:,c)
      hbuf => otape(t)%hlist(f)%hbuf(:,:,c)

      dimind = otape(t)%hlist(f)%field%get_dims(c)

      ! See notes above about validity of avg_subcol_field
      if (otape(t)%hlist(f)%field%is_subcol) then
        if (present(avg_subcol_field)) then
          call endrun('OUTFLD: Cannot average '//trim(fname)//', subcolumn output was requested in addfld')
        end if
        avg_subcols = .false.
      else if (otape(t)%hlist(f)%field%decomp_type == phys_decomp) then
        if (present(avg_subcol_field)) then
          avg_subcols = avg_subcol_field
        else
          avg_subcols = .false.
        end if
      else ! Any dynamics decomposition
        if (present(avg_subcol_field)) then
          call endrun('OUTFLD: avg_subcol_field only valid for physgrid')
        else
          avg_subcols = .false.
        end if
      end if

      begdim2 = otape(t)%hlist(f)%field%begdim2
      enddim2 = otape(t)%hlist(f)%field%enddim2
      if (avg_subcols) then
        allocate(afield(pcols, begdim2:enddim2))
        call subcol_field_avg_handler(idim, field, c, afield)
        ! Hack! Avoid duplicating select statement below
        call outfld(fname, afield, pcols, c)
        deallocate(afield)
      else if (otape(t)%hlist(f)%field%is_subcol) then
        ! We have to assume that using mdimnames (e.g., psubcols) is
        ! incompatible with the begdimx, enddimx usage (checked in addfld)
        ! Since psubcols is included in levels, take that out
        endi = (enddim2 - begdim2 + 1) / psubcols
        allocate(ufield(pcols, psubcols, endi))
        allocate(afield(pcols*psubcols, endi))
        do j = 1, endi
          do i = 1, idim
            afield(i, j) = field(i, j)
          end do
        end do
        ! Initialize unused aray locations.
        if (idim < pcols*psubcols) then
          if (flag_xyfill) then
            afield(idim+1:pcols*psubcols, :) = fillvalue
          else
            afield(idim+1:pcols*psubcols, :) = 0.0_r8
          end if
        end if
        if (flag_xyfill) then
          call subcol_unpack(c, afield, ufield, fillvalue)
        else
          call subcol_unpack(c, afield, ufield)
        end if
        deallocate(afield)
        select case (avgflag)

        case ('I') ! Instantaneous
          call hbuf_accum_inst(hbuf, ufield, nacs, dimind, pcols,        &
               flag_xyfill, fillvalue)
          
        case ('A') ! Time average
          call hbuf_accum_add(hbuf, ufield, nacs, dimind, pcols,         &
               flag_xyfill, fillvalue)

        case ('B') ! Time average only 00z values
          call hbuf_accum_add00z(hbuf, ufield, nacs, dimind, pcols,      &
               flag_xyfill, fillvalue)

        case ('X') ! Maximum over time
          call hbuf_accum_max (hbuf, ufield, nacs, dimind, pcols,        &
               flag_xyfill, fillvalue)

        case ('M') ! Minimum over time
          call hbuf_accum_min(hbuf, ufield, nacs, dimind, pcols,         &
               flag_xyfill, fillvalue)

        case ('L')
          call hbuf_accum_addlcltime(hbuf, ufield, nacs, dimind, pcols,   &
               flag_xyfill, fillvalue, c,                                &
               otape(t)%hlist(f)%field%decomp_type,                      &
               lcltod_start(t), lcltod_stop(t))

        case default
          call endrun ('OUTFLD: invalid avgflag='//avgflag)

        end select
        deallocate(ufield)
      else
        select case (avgflag)

        case ('I') ! Instantaneous
          call hbuf_accum_inst(hbuf, field, nacs, dimind, idim,          &
               flag_xyfill, fillvalue)

        case ('A') ! Time average
          call hbuf_accum_add(hbuf, field, nacs, dimind, idim,           &
               flag_xyfill, fillvalue)

        case ('B') ! Time average only 00z values
          call hbuf_accum_add00z(hbuf, field, nacs, dimind, idim,        &
               flag_xyfill, fillvalue)

        case ('X') ! Maximum over time
          call hbuf_accum_max (hbuf, field, nacs, dimind, idim,          &
               flag_xyfill, fillvalue)

        case ('M') ! Minimum over time
          call hbuf_accum_min(hbuf, field, nacs, dimind, idim,           &
               flag_xyfill, fillvalue)

        case ('L')
          call hbuf_accum_addlcltime(hbuf, field, nacs, dimind, idim,    &
               flag_xyfill, fillvalue, c,                                &
               otape(t)%hlist(f)%field%decomp_type,                      &
               lcltod_start(t), lcltod_stop(t))

        case default
          call endrun ('OUTFLD: invalid avgflag='//avgflag)

        end select
      end if

    end do

    return
  end subroutine outfld

  !#######################################################################

  subroutine get_field_properties(fname, found, tape_out, ff_out)

    implicit none
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: If fname is active, lookup and return field information
    ! 
    ! Method: Check 'masterlist' whether the requested field 'fname' is active
    !         on one or more history tapes, and if so, return the requested
    !         field information
    ! 
    ! Author: goldy
    ! 
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*),   intent(in)  :: fname ! Field name--should be 8 chars long
    logical,            intent(out) :: found ! Set to true if fname is active
    type(active_entry), pointer, optional :: tape_out(:)
    integer,            intent(out), optional :: ff_out

    !
    ! Local variables
    !
    character*(max_fieldname_len) :: fname_loc  ! max-char equivalent of fname
    integer :: t, ff          ! tape, masterindex indices
     !-----------------------------------------------------------------------

    ! Need to re-cast the field name so that the hashing works #hackalert
    fname_loc = fname
    ff = get_masterlist_indx(fname_loc)

    ! Set found to .false. so we can return early if fname is not active
    found = .false.
    if (present(tape_out)) then
      nullify(tape_out)
    end if
    if (present(ff_out)) then
      ff_out = -1
    end if

    !
    !  If ( ff < 0 ), the field is not defined on the masterlist. This check
    !  is necessary because of coding errors calling outfld without first defining
    !  the field on masterlist.
    !
    if ( ff < 0 ) then
      return
    end if
    !
    !  Next, check to see whether this field is active on one or more history
    !  tapes.
    !
    if ( .not. masterlist(ff)%thisentry%act_sometape )  then
      return
    end if
    !
    ! Note, the field may be on any or all of the history files (primary
    ! and auxiliary).
    !

    do t=1, ptapes
      if (masterlist(ff)%thisentry%actflag(t)) then
        found    =  .true.
        if (present(tape_out)) then
          tape_out => history_tape
        end if
        if (present(ff_out)) then
          ff_out   =  ff
        end if
        ! We found the info so we are done with the loop
        exit
      end if
    end do

  end subroutine get_field_properties

  !#######################################################################

  logical function is_initfile (file_index)
    !
    !------------------------------------------------------------------------ 
    ! 
    ! Purpose: to determine:
    !
    !   a) if an IC file is active in this model run at all
    !       OR,
    !   b) if it is active, is the current file index referencing the IC file
    !      (IC file is always at ptapes)
    ! 
    !------------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(in), optional :: file_index ! index of file in question

    is_initfile = .false.

    if (present(file_index)) then
      if (inithist /= 'NONE' .and. file_index == ptapes) is_initfile = .true.
    else
      if (inithist /= 'NONE'                           ) is_initfile = .true.
    end if

    return

  end function is_initfile

  !#######################################################################

  integer function strcmpf (name1, name2)
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Return the lexical difference between two strings
    ! 
    ! Method: Use ichar() intrinsic as we loop through the names
    ! 
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=max_fieldname_len), intent(in) :: name1, name2 ! strings to compare
    integer n                                     ! loop index

    do n=1,max_fieldname_len
      strcmpf = ichar(name1(n:n)) - ichar(name2(n:n))
      if (strcmpf /= 0) exit
    end do

    return
  end function strcmpf

  !#######################################################################

  subroutine h_inquire (t)
    use pio,           only: pio_inq_varid, pio_inq_attlen
    use cam_pio_utils, only: cam_pio_handle_error
   !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Ensure that the proper variables are on a history file
    ! 
    ! Method: Issue the appropriate netcdf wrapper calls
    ! 
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(in) :: t   ! tape index
    !
    ! Local workspace
    !
    integer                  :: f            ! field index
    integer                  :: ierr
    integer                  :: i
    integer                  :: num_patches
    integer(pio_offset_kind) :: mdimsize
    character(len=max_chars) :: fldname, fname_tmp, basename

    !
    !
    ! Dimension id's
    !
    tape => history_tape



    !
    ! Create variables for model timing and header information 
    !
    if(.not. is_satfile(t)) then
      ierr=pio_inq_varid (tape(t)%File,'ndcur   ',    tape(t)%ndcurid)
      ierr=pio_inq_varid (tape(t)%File,'nscur   ',    tape(t)%nscurid)
      ierr=pio_inq_varid (tape(t)%File,'nsteph  ',    tape(t)%nstephid)

      ierr=pio_inq_varid (tape(t)%File,'time_bnds',   tape(t)%tbndid)
      ierr=pio_inq_varid (tape(t)%File,'date_written',tape(t)%date_writtenid)
      ierr=pio_inq_varid (tape(t)%File,'time_written',tape(t)%time_writtenid)
#if ( defined BFB_CAM_SCAM_IOP )
      ierr=pio_inq_varid (tape(t)%File,'tsec    ',tape(t)%tsecid)
      ierr=pio_inq_varid (tape(t)%File,'bdate   ',tape(t)%bdateid)
#endif
      if (.not. is_initfile(file_index=t) ) then
        ! Don't write the GHG/Solar forcing data to the IC file.  It is never
        ! read from that file so it's confusing to have it there.
        ierr=pio_inq_varid (tape(t)%File,'co2vmr  ',    tape(t)%co2vmrid)
        ierr=pio_inq_varid (tape(t)%File,'ch4vmr  ',    tape(t)%ch4vmrid)
        ierr=pio_inq_varid (tape(t)%File,'n2ovmr  ',    tape(t)%n2ovmrid)
        ierr=pio_inq_varid (tape(t)%File,'f11vmr  ',    tape(t)%f11vmrid)
        ierr=pio_inq_varid (tape(t)%File,'f12vmr  ',    tape(t)%f12vmrid)
        ierr=pio_inq_varid (tape(t)%File,'sol_tsi ',    tape(t)%sol_tsiid)
        if (solar_parms_on) then
          ierr=pio_inq_varid (tape(t)%File,'f107    ',    tape(t)%f107id)
          ierr=pio_inq_varid (tape(t)%File,'f107a   ',    tape(t)%f107aid)
          ierr=pio_inq_varid (tape(t)%File,'kp      ',    tape(t)%kpid)
          ierr=pio_inq_varid (tape(t)%File,'ap      ',    tape(t)%apid)
        endif
      end if
    end if
    ierr=pio_inq_varid (tape(t)%File,'date    ',    tape(t)%dateid)
    ierr=pio_inq_varid (tape(t)%File,'datesec ',    tape(t)%datesecid)
    ierr=pio_inq_varid (tape(t)%File,'time    ',    tape(t)%timeid)


    !
    ! Obtain variable name from ID which was read from restart file
    !
    do f=1,nflds(t)
      if(.not. associated(tape(t)%hlist(f)%varid)) then
        if (associated(tape(t)%patches)) then
          allocate(tape(t)%hlist(f)%varid(size(tape(t)%patches)))
        else
          allocate(tape(t)%hlist(f)%varid(1))
        end if
      end if
      !
      ! If this field will be put out as columns then get column names for field
      !
      if (associated(tape(t)%patches)) then
        num_patches = size(tape(t)%patches)
        fldname = strip_suffix(tape(t)%hlist(f)%field%name)
        do i = 1, num_patches
          fname_tmp = trim(fldname)
          call tape(t)%patches(i)%field_name(fname_tmp)
          ierr = pio_inq_varid(tape(t)%File, trim(fname_tmp), tape(t)%hlist(f)%varid(i))
          call cam_pio_handle_error(ierr, 'H_INQUIRE: Error getting ID for '//trim(fname_tmp))
          ierr = pio_get_att(tape(t)%File, tape(t)%hlist(f)%varid(i), 'basename', basename)
          call cam_pio_handle_error(ierr, 'H_INQUIRE: Error getting basename for '//trim(fname_tmp))
          if (trim(fldname) /= trim(basename)) then
            call endrun('H_INQUIRE: basename ('//trim(basename)//') does not match fldname ('//trim(fldname)//')')
          end if
        end do
      else
        fldname = tape(t)%hlist(f)%field%name
        ierr = pio_inq_varid(tape(t)%File, trim(fldname), tape(t)%hlist(f)%varid(1))
        call cam_pio_handle_error(ierr, 'H_INQUIRE: Error getting ID for '//trim(fldname))
      end if
      if(tape(t)%hlist(f)%field%numlev>1) then
        ierr = pio_inq_attlen(tape(t)%File,tape(t)%hlist(f)%varid(1),'mdims', mdimsize)
        if(.not. associated(tape(t)%hlist(f)%field%mdims)) then
          allocate(tape(t)%hlist(f)%field%mdims(mdimsize))
        end if
        ierr=pio_get_att(tape(t)%File,tape(t)%hlist(f)%varid(1),'mdims', &
             tape(t)%hlist(f)%field%mdims(1:mdimsize))
        if(mdimsize>maxvarmdims) maxvarmdims=mdimsize
      end if

    end do

    if(masterproc) then
      write(iulog,*)'H_INQUIRE: Successfully opened netcdf file '
    end if

    return
  end subroutine h_inquire

  !#######################################################################

  subroutine add_default (name, tindex, flag)
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Add a field to the default "on" list for a given history file
    ! 
    ! Method: 
    ! 
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*), intent(in) :: name  ! field name
    character(len=1), intent(in) :: flag  ! averaging flag

    integer, intent(in) :: tindex         ! history tape index
    !
    ! Local workspace
    !
    integer :: t            ! file index
    type(master_entry), pointer :: listentry

    if (htapes_defined) then
      call endrun ('ADD_DEFAULT: Attempt to add hist default '//trim(name)//' after history files set')
    end if
    !
    ! Check validity of input arguments
    !
    if (tindex > ptapes) then
      write(iulog,*)'ADD_DEFAULT: tape index=', tindex, ' is too big'
      call endrun
    end if

    ! Add to IC file if tindex = 0, reset to ptapes
    if (tindex == 0) then
      t = ptapes
      if ( .not. is_initfile(file_index=t) ) return
    else
      t = tindex
    end if

    if (verify(flag, HIST_AVG_FLAGS) /= 0) then
      call endrun ('ADD_DEFAULT: unknown averaging flag='//flag)
    end if
    !
    ! Look through master list for input field name.  When found, set active
    ! flag for that tape to true.  Also set averaging flag if told to use other
    ! than default.
    !
    listentry => get_entry_by_name(masterlinkedlist, trim(name))
    if(.not.associated(listentry)) then
      call endrun ('ADD_DEFAULT: field='//name//' not found')
    end if
    listentry%actflag(t) = .true.
    if (flag /= ' ') then
      listentry%avgflag(t) = flag
      call AvgflagToString(flag, listentry%time_op(t))
    end if

    return
  end subroutine add_default

  !#######################################################################

  subroutine h_override (t)
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Override default history tape contents for a specific tape
    !
    ! Method: Copy the flag into the master field list
    ! 
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(in) :: t         ! history tape index
    !
    ! Local workspace
    !
    character(len=1) :: avgflg       ! lcl equiv of avgflag_pertape(t) (to address xlf90 compiler bug)

    type(master_entry), pointer :: listentry


    avgflg = avgflag_pertape(t)


    listentry=>masterlinkedlist
    do while(associated(listentry))
      call AvgflagToString(avgflg, listentry%time_op(t))
      listentry%avgflag(t) = avgflag_pertape(t)
      listentry=>listentry%next_entry
    end do

  end subroutine h_override

  !#######################################################################

  subroutine h_define (t, restart)
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Define contents of history file t
    ! 
    ! Method: Issue the required netcdf wrapper calls to define the history file contents
    ! 
    !-----------------------------------------------------------------------
    use cam_grid_support, only: cam_grid_header_info_t
    use cam_grid_support, only: cam_grid_write_attr, cam_grid_write_var
    use time_manager,     only: get_step_size, get_ref_date, timemgr_get_calendar_cf
    use filenames,        only: caseid
    use cam_abortutils,   only: endrun
    use cam_pio_utils,    only: vdesc_ptr, cam_pio_handle_error, cam_pio_def_dim
    use cam_pio_utils,    only: cam_pio_createfile, cam_pio_def_var
    use sat_hist,         only: sat_hist_define

    !-----------------------------------------------------------------------

    !
    ! Input arguments
    !
    integer, intent(in) :: t   ! tape index
    logical, intent(in) :: restart
    !
    ! Local workspace
    !
    integer :: i, j            ! longitude, latitude indices
    integer :: grd             ! indices for looping through grids
    integer :: f               ! field index
    integer :: ncreal          ! real data type for output
    integer :: dtime           ! timestep size
    integer :: ndbase = 0      ! days component of base time
    integer :: nsbase = 0      ! seconds component of base time
    integer :: nbdate          ! base date in yyyymmdd format
    integer :: nbsec           ! time of day component of base date [seconds]
    integer :: yr, mon, day    ! year, month, day components of a date

    character(len=max_chars) :: str       ! character temporary 
    character(len=max_chars) :: fname_tmp ! local copy of field name
    character(len=max_chars) :: calendar  ! Calendar type
    character(len=16)        :: time_per_freq
    character(len=128)       :: errormsg

    integer :: ret                        ! function return value

    !
    ! netcdf dimensions
    !
    integer :: chardim            ! character dimension id
    integer :: dimenchar(2)       ! character dimension ids
    integer :: nacsdims(2)        ! dimension ids for nacs (used in restart file)
    integer :: bnddim             ! bounds dimension id
    integer :: timdim             ! unlimited dimension id

    integer :: dimindex(8)        ! dimension ids for variable declaration
    integer :: dimids_tmp(8)      ! dimension ids for variable declaration

    !
    ! netcdf variables
    !
    ! A structure to hold the horizontal dimension and coordinate info
    type(cam_grid_header_info_t), allocatable :: header_info(:)
    ! For satellite files and column output
    type(vdesc_ptr), allocatable :: latvar(:)    ! latitude variable ids
    type(vdesc_ptr), allocatable :: lonvar(:)    ! longitude variable ids

    type(var_desc_t), pointer        :: varid => NULL() ! temporary variable descriptor
    integer                          :: num_hdims, fdims
    integer                          :: num_patches ! How many entries for a field on this tape?
    integer,          pointer        :: mdims(:) => NULL()
    integer                          :: mdimsize
    integer                          :: ierr
    integer,          allocatable    :: mdimids(:)
    integer                          :: amode
    logical                          :: interpolate
    logical                          :: patch_output

    if(restart) then
      tape => restarthistory_tape
      if(masterproc) write(iulog,*)'Opening netcdf history restart file ', trim(hrestpath(t))
    else
      tape => history_tape
      if(masterproc) write(iulog,*)'Opening netcdf history file ', trim(nhfil(t))
    end if

    amode = PIO_CLOBBER

    if(restart) then
      call cam_pio_createfile (tape(t)%File, hrestpath(t), amode)
    else
      call cam_pio_createfile (tape(t)%File, nhfil(t), amode)
    end if
    if(is_satfile(t)) then
      interpolate = .false. ! !!XXgoldyXX: Do we ever want to support this?
      patch_output = .false.
      call cam_pio_def_dim(tape(t)%File, 'ncol', pio_unlimited, timdim)
      call cam_pio_def_dim(tape(t)%File, 'nbnd', 2, bnddim)

      allocate(latvar(1), lonvar(1))
      allocate(latvar(1)%vd, lonvar(1)%vd)
      call cam_pio_def_var(tape(t)%File, 'lat', pio_double, (/timdim/),       &
           latvar(1)%vd)
      ierr=pio_put_att (tape(t)%File, latvar(1)%vd, 'long_name', 'latitude')
      ierr=pio_put_att (tape(t)%File, latvar(1)%vd, 'units', 'degrees_north')

      call cam_pio_def_var(tape(t)%File, 'lon', pio_double, (/timdim/),       &
           lonvar(1)%vd)
      ierr=pio_put_att (tape(t)%File, lonvar(1)%vd,'long_name','longitude')
      ierr=pio_put_att (tape(t)%File, lonvar(1)%vd,'units','degrees_east')

    else
      !
      ! Setup netcdf file - create the dimensions of lat,lon,time,level
      !     
      ! interpolate is only supported for unstructured dycores
      interpolate = (interpolate_output(t) .and. (.not. restart))
      patch_output = (associated(tape(t)%patches) .and. (.not. restart))

      ! First define the horizontal grid dims
      ! Interpolation is special in that we ignore the native grids
      if(interpolate) then
        allocate(header_info(1))
        call cam_grid_write_attr(tape(t)%File, interpolate_info(t)%grid_id, header_info(1))
      else if (patch_output) then
        ! We are doing patch (column) output
        if (allocated(header_info)) then
          ! We shouldn't have any header_info yet
          call endrun('H_DEFINE: header_info should not be allocated for patch output')
        end if
        do i = 1, size(tape(t)%patches)
          call tape(t)%patches(i)%write_attrs(tape(t)%File)
        end do
      else
        allocate(header_info(size(tape(t)%grid_ids)))
        do i = 1, size(tape(t)%grid_ids)
          call cam_grid_write_attr(tape(t)%File, tape(t)%grid_ids(i), header_info(i))
        end do
      end if   ! interpolate

      ! Define the unlimited time dim
      call cam_pio_def_dim(tape(t)%File, 'time', pio_unlimited, timdim)
      call cam_pio_def_dim(tape(t)%File, 'nbnd', 2, bnddim)
      call cam_pio_def_dim(tape(t)%File, 'chars', 8, chardim)
    end if   ! is satfile

    ! Populate the history coordinate (well, mdims anyway) attributes
    ! This routine also allocates the mdimids array
    call write_hist_coord_attrs(tape(t)%File, bnddim, mdimids, restart)

    call get_ref_date(yr, mon, day, nbsec)
    nbdate = yr*10000 + mon*100 + day
    ierr=pio_def_var (tape(t)%File,'time',pio_double,(/timdim/),tape(t)%timeid)
    ierr=pio_put_att (tape(t)%File, tape(t)%timeid, 'long_name', 'time')
    str = 'days since ' // date2yyyymmdd(nbdate) // ' ' // sec2hms(nbsec)
    ierr=pio_put_att (tape(t)%File, tape(t)%timeid, 'units', trim(str))

    calendar = timemgr_get_calendar_cf()
    ierr=pio_put_att (tape(t)%File, tape(t)%timeid, 'calendar', trim(calendar))


    ierr=pio_def_var (tape(t)%File,'date    ',pio_int,(/timdim/),tape(t)%dateid)
    str = 'current date (YYYYMMDD)'
    ierr=pio_put_att (tape(t)%File, tape(t)%dateid, 'long_name', trim(str))


    ierr=pio_def_var (tape(t)%File,'datesec ',pio_int,(/timdim/), tape(t)%datesecid)
    str = 'current seconds of current date'
    ierr=pio_put_att (tape(t)%File, tape(t)%datesecid, 'long_name', trim(str))

    !     
    ! Character header information 
    !
    str = 'CF-1.0'
    ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'Conventions', trim(str))
    ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'source', 'CAM')
#if ( defined BFB_CAM_SCAM_IOP )
    ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'CAM_GENERATED_FORCING','create SCAM IOP dataset')
#endif
    ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'case',caseid)
    ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'title',ctitle)
    ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'logname',logname)
    ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'host', host)
    ierr= pio_put_att (tape(t)%File, PIO_GLOBAL, 'Version', &
         '$Name$')
    ierr= pio_put_att (tape(t)%File, PIO_GLOBAL, 'revision_Id', &
         '$Id$')
    ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'initial_file', ncdata)
    ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'topography_file', bnd_topo)

    ! Determine what time period frequency is being output for each file
    ! Note that nhtfrq is now in timesteps
    dtime = get_step_size()
    if (nhtfrq(t) == 0) then                                !month 
      time_per_freq = 'month_1'
    else if (mod(nhtfrq(t)*dtime,86400) == 0) then          ! day
      write(time_per_freq,999) 'day_',nhtfrq(t)*dtime/86400
    else if (mod(nhtfrq(t)*dtime,3600) == 0) then           ! hour
      write(time_per_freq,999) 'hour_',(nhtfrq(t)*dtime)/3600
    else if (mod(nhtfrq(t)*dtime,60) == 0) then           ! hour
      write(time_per_freq,999) 'minute_',(nhtfrq(t)*dtime)/60
    else                                                    ! second
      write(time_per_freq,999) 'second_',nhtfrq(t)*dtime
    end if
999 format(a,i0)

    ierr=pio_put_att (tape(t)%File, PIO_GLOBAL, 'time_period_freq', trim(time_per_freq))

    if(.not. is_satfile(t)) then

      ierr=pio_put_att (tape(t)%File, tape(t)%timeid, 'bounds', 'time_bnds')

      ierr=pio_def_var (tape(t)%File,'time_bnds',pio_double,(/bnddim,timdim/),tape(t)%tbndid)
      ierr=pio_put_att (tape(t)%File, tape(t)%tbndid, 'long_name', 'time interval endpoints')
      !
      ! Character
      !
      dimenchar(1) = chardim
      dimenchar(2) = timdim
      ierr=pio_def_var (tape(t)%File,'date_written',PIO_CHAR,dimenchar, tape(t)%date_writtenid)
      ierr=pio_def_var (tape(t)%File,'time_written',PIO_CHAR,dimenchar, tape(t)%time_writtenid)
      !
      ! Integer Header
      !

      ierr=pio_def_var (tape(t)%File,'ndbase',PIO_INT,tape(t)%ndbaseid)
      str = 'base day'
      ierr=pio_put_att (tape(t)%File, tape(t)%ndbaseid, 'long_name', trim(str))

      ierr=pio_def_var (tape(t)%File,'nsbase',PIO_INT,tape(t)%nsbaseid)
      str = 'seconds of base day'
      ierr=pio_put_att (tape(t)%File, tape(t)%nsbaseid, 'long_name', trim(str))

      ierr=pio_def_var (tape(t)%File,'nbdate',PIO_INT,tape(t)%nbdateid)
      str = 'base date (YYYYMMDD)'
      ierr=pio_put_att (tape(t)%File, tape(t)%nbdateid, 'long_name', trim(str))

#if ( defined BFB_CAM_SCAM_IOP )
      ierr=pio_def_var (tape(t)%File,'bdate',PIO_INT,tape(t)%bdateid)
      str = 'base date (YYYYMMDD)'
      ierr=pio_put_att (tape(t)%File, tape(t)%bdateid, 'long_name', trim(str))
#endif
      ierr=pio_def_var (tape(t)%File,'nbsec',PIO_INT,tape(t)%nbsecid)
      str = 'seconds of base date'
      ierr=pio_put_att (tape(t)%File, tape(t)%nbsecid, 'long_name', trim(str))

      ierr=pio_def_var (tape(t)%File,'mdt',PIO_INT,tape(t)%mdtid)
      ierr=pio_put_att (tape(t)%File, tape(t)%mdtid, 'long_name', 'timestep')
      ierr=pio_put_att (tape(t)%File, tape(t)%mdtid, 'units', 's')

      !
      ! Create variables for model timing and header information 
      !

      ierr=pio_def_var (tape(t)%File,'ndcur   ',pio_int,(/timdim/),tape(t)%ndcurid)
      str = 'current day (from base day)'
      ierr=pio_put_att (tape(t)%File, tape(t)%ndcurid, 'long_name', trim(str))

      ierr=pio_def_var (tape(t)%File,'nscur   ',pio_int,(/timdim/),tape(t)%nscurid)
      str = 'current seconds of current day'
      ierr=pio_put_att (tape(t)%File, tape(t)%nscurid, 'long_name', trim(str))


      if (.not. is_initfile(file_index=t)) then
        ! Don't write the GHG/Solar forcing data to the IC file.
        ierr=pio_def_var (tape(t)%File,'co2vmr  ',pio_double,(/timdim/),tape(t)%co2vmrid)
        str = 'co2 volume mixing ratio'
        ierr=pio_put_att (tape(t)%File, tape(t)%co2vmrid, 'long_name', trim(str))

        ierr=pio_def_var (tape(t)%File,'ch4vmr  ',pio_double,(/timdim/),tape(t)%ch4vmrid)
        str = 'ch4 volume mixing ratio'
        ierr=pio_put_att (tape(t)%File, tape(t)%ch4vmrid, 'long_name', trim(str))

        ierr=pio_def_var (tape(t)%File,'n2ovmr  ',pio_double,(/timdim/),tape(t)%n2ovmrid)
        str = 'n2o volume mixing ratio'
        ierr=pio_put_att (tape(t)%File, tape(t)%n2ovmrid, 'long_name', trim(str))

        ierr=pio_def_var (tape(t)%File,'f11vmr  ',pio_double,(/timdim/),tape(t)%f11vmrid)
        str = 'f11 volume mixing ratio'
        ierr=pio_put_att (tape(t)%File, tape(t)%f11vmrid, 'long_name', trim(str))

        ierr=pio_def_var (tape(t)%File,'f12vmr  ',pio_double,(/timdim/),tape(t)%f12vmrid)
        str = 'f12 volume mixing ratio'
        ierr=pio_put_att (tape(t)%File, tape(t)%f12vmrid, 'long_name', trim(str))

        ierr=pio_def_var (tape(t)%File,'sol_tsi ',pio_double,(/timdim/),tape(t)%sol_tsiid)
        str = 'total solar irradiance'
        ierr=pio_put_att (tape(t)%File, tape(t)%sol_tsiid, 'long_name', trim(str))
        str = 'W/m2'
        ierr=pio_put_att (tape(t)%File, tape(t)%sol_tsiid, 'units', trim(str))

        if (solar_parms_on) then
          ! solar / geomagetic activity indices...
          ierr=pio_def_var (tape(t)%File,'f107',pio_double,(/timdim/),tape(t)%f107id)
          str = '10.7 cm solar radio flux (F10.7)'
          ierr=pio_put_att (tape(t)%File, tape(t)%f107id, 'long_name', trim(str))
          str = '10^-22 W m^-2 Hz^-1'
          ierr=pio_put_att (tape(t)%File, tape(t)%f107id, 'units', trim(str))

          ierr=pio_def_var (tape(t)%File,'f107a',pio_double,(/timdim/),tape(t)%f107aid)
          str = '81-day centered mean of 10.7 cm solar radio flux (F10.7)'
          ierr=pio_put_att (tape(t)%File, tape(t)%f107aid, 'long_name', trim(str))

          ierr=pio_def_var (tape(t)%File,'kp',pio_double,(/timdim/),tape(t)%kpid)
          str = 'Daily planetary K geomagnetic index'
          ierr=pio_put_att (tape(t)%File, tape(t)%kpid, 'long_name', trim(str))

          ierr=pio_def_var (tape(t)%File,'ap',pio_double,(/timdim/),tape(t)%apid)
          str = 'Daily planetary A geomagnetic index'
          ierr=pio_put_att (tape(t)%File, tape(t)%apid, 'long_name', trim(str))
        endif
      end if


#if ( defined BFB_CAM_SCAM_IOP )
      ierr=pio_def_var (tape(t)%File,'tsec ',pio_int,(/timdim/), tape(t)%tsecid)
      str = 'current seconds of current date needed for scam'
      ierr=pio_put_att (tape(t)%File, tape(t)%tsecid, 'long_name', trim(str))
#endif
      ierr=pio_def_var (tape(t)%File,'nsteph  ',pio_int,(/timdim/),tape(t)%nstephid)
      str = 'current timestep'
      ierr=pio_put_att (tape(t)%File, tape(t)%nstephid, 'long_name', trim(str))
    end if ! .not. is_satfile

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Create variables and attributes for field list
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do f = 1, nflds(t)

      !! Collect some field properties
      call AvgflagToString(tape(t)%hlist(f)%avgflag, tape(t)%hlist(f)%time_op)

      if ((tape(t)%hlist(f)%hwrt_prec == 8) .or. restart) then
        ncreal = pio_double
      else
        ncreal = pio_real 
      end if

      if(associated(tape(t)%hlist(f)%field%mdims)) then
        mdims => tape(t)%hlist(f)%field%mdims
        mdimsize = size(mdims)
      else if(tape(t)%hlist(f)%field%numlev > 1) then
        call endrun('mdims not defined for variable '//trim(tape(t)%hlist(f)%field%name))
      else
        mdimsize=0   
      end if

      ! num_patches will loop through the number of patches (or just one
      !             for the whole grid) for this field for this tape
      if (patch_output) then
        num_patches = size(tape(t)%patches)
      else
        num_patches = 1
      end if
      if(.not.associated(tape(t)%hlist(f)%varid)) then
        allocate(tape(t)%hlist(f)%varid(num_patches))
      end if
      fname_tmp = strip_suffix(tape(t)%hlist(f)%field%name)

      if(is_satfile(t)) then
        num_hdims=0
        nfils(t)=1
        call sat_hist_define(tape(t)%File)
      else if (interpolate) then
        ! Interpolate can't use normal grid code since we are forcing fields
        ! to use interpolate decomp
        if (.not. allocated(header_info)) then
          ! Safety check
          call endrun('h_define: header_info not allocated')
        end if
        num_hdims = 2
        do i = 1, num_hdims
          dimindex(i) = header_info(1)%get_hdimid(i)
          nacsdims(i) = header_info(1)%get_hdimid(i)
        end do
      else if (patch_output) then
        ! All patches for this variable should be on the same grid
        num_hdims = tape(t)%patches(1)%num_hdims(tape(t)%hlist(f)%field%decomp_type)
      else
        ! Normal grid output
        ! Find appropriate grid in header_info
        if (.not. allocated(header_info)) then
          ! Safety check
          call endrun('h_define: header_info not allocated')
        end if
        grd = -1
        do i = 1, size(header_info)
          if (header_info(i)%get_gridid() == tape(t)%hlist(f)%field%decomp_type) then
            grd = i
            exit
          end if
        end do
        if (grd < 0) then
          write(errormsg, '(a,i0,2a)') 'grid, ',tape(t)%hlist(f)%field%decomp_type,', not found for ',trim(fname_tmp)
          call endrun('H_DEFINE: '//errormsg)
        end if
        num_hdims = header_info(grd)%num_hdims()
        do i = 1, num_hdims
          dimindex(i) = header_info(grd)%get_hdimid(i)
          nacsdims(i) = header_info(grd)%get_hdimid(i)
        end do
      end if     ! is_satfile

      !
      !  Create variables and atributes for fields written out as columns
      !         

      do i = 1, num_patches
        fname_tmp = strip_suffix(tape(t)%hlist(f)%field%name)
        varid => tape(t)%hlist(f)%varid(i)
        dimids_tmp = dimindex
        ! Figure the dimension ID array for this field
        ! We have defined the horizontal grid dimensions in dimindex
        fdims = num_hdims
        do j = 1, mdimsize
          fdims = fdims + 1
          dimids_tmp(fdims) = mdimids(mdims(j))
        end do
        if(.not. restart) then
          ! Only add time dimension if this is not a restart history tape
          fdims = fdims + 1
          dimids_tmp(fdims) = timdim
        end if
        if (patch_output) then
          ! For patch output, we need new dimension IDs and a different name
          call tape(t)%patches(i)%get_var_data(fname_tmp,                     &
               dimids_tmp(1:fdims), tape(t)%hlist(f)%field%decomp_type)
        end if
        ! Define the variable
        call cam_pio_def_var(tape(t)%File, trim(fname_tmp), ncreal,           &
             dimids_tmp(1:fdims), varid)
        if (mdimsize > 0) then
          ierr = pio_put_att(tape(t)%File, varid, 'mdims', mdims(1:mdimsize))
          call cam_pio_handle_error(ierr, 'h_define: cannot define mdims for '//trim(fname_tmp))
        end if
        str = tape(t)%hlist(f)%field%sampling_seq
        if (len_trim(str) > 0) then
          ierr = pio_put_att(tape(t)%File, varid, 'Sampling_Sequence', trim(str))
          call cam_pio_handle_error(ierr, 'h_define: cannot define Sampling_Sequence for '//trim(fname_tmp))
        end if

        if (tape(t)%hlist(f)%field%flag_xyfill) then
          ! Add both _FillValue and missing_value to cover expectations
          !     of various applications.
          ! The attribute type must match the data type.
          if ((tape(t)%hlist(f)%hwrt_prec == 8) .or. restart) then
            ierr = pio_put_att(tape(t)%File, varid, '_FillValue',             &
                 tape(t)%hlist(f)%field%fillvalue)
            call cam_pio_handle_error(ierr,                                   &
                 'h_define: cannot define _FillValue for '//trim(fname_tmp))
            ierr = pio_put_att(tape(t)%File, varid, 'missing_value',          &
                 tape(t)%hlist(f)%field%fillvalue)
            call cam_pio_handle_error(ierr,                                   &
                 'h_define: cannot define missing_value for '//trim(fname_tmp))
          else
            ierr = pio_put_att(tape(t)%File, varid, '_FillValue',             &
                 REAL(tape(t)%hlist(f)%field%fillvalue,r4))
            call cam_pio_handle_error(ierr,                                   &
                 'h_define: cannot define _FillValue for '//trim(fname_tmp))
            ierr = pio_put_att(tape(t)%File, varid, 'missing_value',          &
                 REAL(tape(t)%hlist(f)%field%fillvalue,r4))
            call cam_pio_handle_error(ierr,                                   &
                 'h_define: cannot define missing_value for '//trim(fname_tmp))
          end if
        end if

        str = tape(t)%hlist(f)%field%units
        if (len_trim(str) > 0) then
          ierr=pio_put_att (tape(t)%File, varid, 'units', trim(str))
          call cam_pio_handle_error(ierr,                                     &
               'h_define: cannot define units for '//trim(fname_tmp))
        end if

        str = tape(t)%hlist(f)%field%long_name
        ierr=pio_put_att (tape(t)%File, varid, 'long_name', trim(str))
        call cam_pio_handle_error(ierr,                                       &
             'h_define: cannot define long_name for '//trim(fname_tmp))
        !
        ! Assign field attributes defining valid levels and averaging info
        !
        str = tape(t)%hlist(f)%time_op
        select case (str)
        case ('mean', 'maximum', 'minimum' )
          ierr = pio_put_att(tape(t)%File, varid, 'cell_methods', 'time: '//str)
          call cam_pio_handle_error(ierr,                                     &
               'h_define: cannot define cell_methods for '//trim(fname_tmp))
        end select
        if (patch_output) then
          ierr = pio_put_att(tape(t)%File, varid, 'basename',                 &
               tape(t)%hlist(f)%field%name)
          call cam_pio_handle_error(ierr,                                     &
               'h_define: cannot define basename for '//trim(fname_tmp))
        end if

        if (restart) then
          ! For restart history files, we need to save accumulation counts
          fname_tmp = trim(fname_tmp)//'_nacs'
          if (.not. associated(tape(t)%hlist(f)%nacs_varid)) then
            allocate(tape(t)%hlist(f)%nacs_varid)
          end if
          if (size(tape(t)%hlist(f)%nacs, 1) > 1) then
            call cam_pio_def_var(tape(t)%File, trim(fname_tmp), pio_int,      &
                 nacsdims(1:num_hdims), tape(t)%hlist(f)%nacs_varid)
          else
            ! Save just one value representing all chunks
            call cam_pio_def_var(tape(t)%File, trim(fname_tmp), pio_int,      &
                 tape(t)%hlist(f)%nacs_varid)
          end if
        end if
      end do ! Loop over output patches
    end do   ! Loop over fields
    !
    deallocate(mdimids)
    ret = pio_enddef(tape(t)%File)

    if(masterproc) then
      write(iulog,*)'H_DEFINE: Successfully opened netcdf file '
    endif
    !
    ! Write time-invariant portion of history header
    !
    if(.not. is_satfile(t)) then
      if(interpolate) then
        call cam_grid_write_var(tape(t)%File, interpolate_info(t)%grid_id)
      else if((.not. patch_output) .or. restart) then
        do i = 1, size(tape(t)%grid_ids)
          call cam_grid_write_var(tape(t)%File, tape(t)%grid_ids(i))
        end do
      else
        ! Patch output
        do i = 1, size(tape(t)%patches)
          call tape(t)%patches(i)%write_vals(tape(t)%File)
        end do
      end if ! interpolate
      if (allocated(lonvar)) then
        deallocate(lonvar)
      end if
      if (allocated(latvar)) then
        deallocate(latvar)
      end if

      dtime = get_step_size()
      ierr = pio_put_var(tape(t)%File, tape(t)%mdtid, (/dtime/))
      call cam_pio_handle_error(ierr, 'h_define: cannot put mdt')
      !
      ! Model date info
      !
      ierr = pio_put_var(tape(t)%File, tape(t)%ndbaseid, (/ndbase/))
      call cam_pio_handle_error(ierr, 'h_define: cannot put ndbase')
      ierr = pio_put_var(tape(t)%File, tape(t)%nsbaseid, (/nsbase/))
      call cam_pio_handle_error(ierr, 'h_define: cannot put nsbase')

      ierr = pio_put_var(tape(t)%File, tape(t)%nbdateid, (/nbdate/))
      call cam_pio_handle_error(ierr, 'h_define: cannot put nbdate')
#if ( defined BFB_CAM_SCAM_IOP )
      ierr = pio_put_var(tape(t)%File, tape(t)%bdateid, (/nbdate/))
      call cam_pio_handle_error(ierr, 'h_define: cannot put bdate')
#endif
      ierr = pio_put_var(tape(t)%File, tape(t)%nbsecid, (/nbsec/))
      call cam_pio_handle_error(ierr, 'h_define: cannot put nbsec')
      !
      ! Reduced grid info
      !

    end if ! .not. is_satfile

    if (allocated(header_info)) then
      do i = 1, size(header_info)
        call header_info(i)%deallocate()
      end do
      deallocate(header_info)
    end if

    ! Write the mdim variable data
    call write_hist_coord_vars(tape(t)%File, restart)

  end subroutine h_define

  !#######################################################################

  subroutine h_normalize (f, t)

    use cam_history_support, only: dim_index_2d

    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Normalize fields on a history file by the number of accumulations
    ! 
    ! Method: Loop over fields on the tape.  Need averaging flag and number of
    !         accumulations to perform normalization.
    ! 
    !-----------------------------------------------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: f       ! field index
    integer, intent(in) :: t       ! tape index
    !
    ! Local workspace
    !
    type (dim_index_2d) :: dimind  ! 2-D dimension index
    integer     :: c               ! chunk (or lat) index
    integer     :: ib, ie    ! beginning and ending indices of first dimension
    integer     :: jb, je    ! beginning and ending indices of second dimension
    integer     :: begdim3, enddim3 ! Chunk or block bounds
    integer     :: k         ! level

    logical     :: flag_xyfill ! non-applicable xy points flagged with fillvalue
    character*1 :: avgflag     ! averaging flag

    call t_startf ('h_normalize')

    call tape(t)%hlist(f)%field%get_bounds(3, begdim3, enddim3)

    !
    ! normalize by number of accumulations for averaged case
    !
    flag_xyfill = tape(t)%hlist(f)%field%flag_xyfill
    avgflag = tape(t)%hlist(f)%avgflag

    do c = begdim3, enddim3
      dimind = tape(t)%hlist(f)%field%get_dims(c)

      ib = dimind%beg1
      ie = dimind%end1
      jb = dimind%beg2
      je = dimind%end2

      if (flag_xyfill) then
        do k = jb, je
          where (tape(t)%hlist(f)%nacs(ib:ie, c) == 0)
            tape(t)%hlist(f)%hbuf(ib:ie,k,c) = tape(t)%hlist(f)%field%fillvalue
          endwhere
        end do
      end if

      if (avgflag == 'A' .or. avgflag == 'B' .or. avgflag == 'L') then
        if (size(tape(t)%hlist(f)%nacs, 1) > 1) then
          do k = jb, je
            where (tape(t)%hlist(f)%nacs(ib:ie,c) /= 0)
              tape(t)%hlist(f)%hbuf(ib:ie,k,c) = &
                   tape(t)%hlist(f)%hbuf(ib:ie,k,c) &
                   / tape(t)%hlist(f)%nacs(ib:ie,c)
            endwhere
          end do
        else if(tape(t)%hlist(f)%nacs(1,c) > 0) then
          do k=jb,je
            tape(t)%hlist(f)%hbuf(ib:ie,k,c) = &
                 tape(t)%hlist(f)%hbuf(ib:ie,k,c) &
                 / tape(t)%hlist(f)%nacs(1,c)
          end do
        end if
      end if
    end do

    call t_stopf ('h_normalize')

    return
  end subroutine h_normalize

  !#######################################################################

  subroutine h_zero (f, t)
    use cam_history_support, only: dim_index_2d
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Zero out accumulation buffers for a tape
    ! 
    ! Method: Loop through fields on the tape
    ! 
    !-----------------------------------------------------------------------
    !
    integer, intent(in) :: f     ! field index
    integer, intent(in) :: t     ! tape index
    !
    ! Local workspace
    !
    type (dim_index_2d) :: dimind   ! 2-D dimension index
    integer :: c                    ! chunk index
    integer :: begdim3              ! on-node chunk or lat start index
    integer :: enddim3              ! on-node chunk or lat end index  

    call t_startf ('h_zero')

    call tape(t)%hlist(f)%field%get_bounds(3, begdim3, enddim3)

    do c = begdim3, enddim3
      dimind = tape(t)%hlist(f)%field%get_dims(c)
      tape(t)%hlist(f)%hbuf(dimind%beg1:dimind%end1,dimind%beg2:dimind%end2,c)=0._r8
    end do
    tape(t)%hlist(f)%nacs(:,:) = 0

    call t_stopf ('h_zero')

    return
  end subroutine h_zero

  !#######################################################################

  subroutine dump_field (f, t, restart)
    use cam_history_support, only: history_patch_t, dim_index_3d
    use cam_grid_support,    only: cam_grid_write_dist_array, cam_grid_dimensions
    use interp_mod,       only : write_interpolated

    ! Dummy arguments
    integer,     intent(in)    :: f
    integer,     intent(in)    :: t
    logical,     intent(in)    :: restart
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Write a variable to a history tape using PIO
    !          For restart tapes, also write the accumulation buffer (nacs)
    !
    !-----------------------------------------------------------------------
    ! Local variables
    integer                          :: ierr
    type(var_desc_t),      pointer   :: varid      ! PIO ID for var
    type(var_desc_t),      pointer   :: compid     ! PIO ID for vector comp.
    integer                          :: compind    ! index of vector comp.
    integer                          :: fdims(8)   ! Field file dim sizes
    integer                          :: frank      ! Field file rank
    integer                          :: nacsrank   ! Field file rank for nacs
    type(dim_index_3d)               :: dimind     ! 3-D dimension index
    integer                          :: adims(3)   ! Field array dim sizes
    integer                          :: nadims     ! # of used adims
    integer                          :: fdecomp
    integer                          :: num_patches
    integer                          :: mdimsize   ! Total # on-node elements
    logical                          :: interpolate
    logical                          :: patch_output
    type(history_patch_t), pointer   :: patchptr

    interpolate = (interpolate_output(t) .and. (.not. restart))
    patch_output = (associated(tape(t)%patches) .and. (.not. restart))

    !!! Get the field's shape and decomposition

    ! Shape on disk
    call tape(t)%hlist(f)%field%get_shape(fdims, frank)

    ! Shape of array
    dimind = tape(t)%hlist(f)%field%get_dims()
    call dimind%dim_sizes(adims)
    if (adims(2) <= 1) then
      adims(2) = adims(3)
      nadims = 2
    else
      nadims = 3
    end if
    fdecomp = tape(t)%hlist(f)%field%decomp_type

    ! num_patches will loop through the number of patches (or just one
    !             for the whole grid) for this field for this tape
    if (patch_output) then
      num_patches = size(tape(t)%patches)
    else
      num_patches = 1
    end if

    do i = 1, num_patches
      varid => tape(t)%hlist(f)%varid(i)

      if (restart) then
        call pio_setframe(tape(t)%File, varid, int(-1,kind=PIO_OFFSET_KIND))
      else
        call pio_setframe(tape(t)%File, varid, int(max(1,nfils(t)),kind=PIO_OFFSET_KIND))
      end if
      if (patch_output) then
        ! We are outputting patches
        patchptr => tape(t)%patches(i)
        if (interpolate) then
          call endrun('dump_field: interpolate incompatible with regional output')
        end if
        call patchptr%write_var(tape(t)%File, fdecomp, adims(1:nadims),       &
             pio_double, tape(t)%hlist(f)%hbuf, varid)
      else
        ! We are doing output via the field's grid
        if (interpolate) then
          mdimsize = tape(t)%hlist(f)%field%enddim2 - tape(t)%hlist(f)%field%begdim2 + 1
          if (mdimsize == 0) then
            mdimsize = tape(t)%hlist(f)%field%numlev
          end if
          if (tape(t)%hlist(f)%field%meridional_complement > 0) then
            compind = tape(t)%hlist(f)%field%meridional_complement
            compid => tape(t)%hlist(compind)%varid(i)
            ! We didn't call set frame on the meridional complement field
            call pio_setframe(tape(t)%File, compid, int(max(1,nfils(t)),kind=PIO_OFFSET_KIND))
            call write_interpolated(tape(t)%File, varid, compid,              &
                 tape(t)%hlist(f)%hbuf, tape(t)%hlist(compind)%hbuf,          &
                 mdimsize, PIO_DOUBLE, fdecomp)
          else if (tape(t)%hlist(f)%field%zonal_complement > 0) then
            ! We don't want to double write so do nothing here
!            compind = tape(t)%hlist(f)%field%zonal_complement
!            compid => tape(t)%hlist(compind)%varid(i)
!            call write_interpolated(tape(t)%File, compid, varid,              &
!                 tape(t)%hlist(compind)%hbuf, tape(t)%hlist(f)%hbuf,          &
!                 mdimsize, PIO_DOUBLE, fdecomp)
          else
            ! Scalar field
            call write_interpolated(tape(t)%File, varid,                      &
                 tape(t)%hlist(f)%hbuf, mdimsize, PIO_DOUBLE, fdecomp)
          end if
        else if (nadims == 2) then
          ! Special case for 2D field (no levels) due to hbuf structure
          call cam_grid_write_dist_array(tape(t)%File, fdecomp,               &
               adims(1:nadims), fdims(1:frank), tape(t)%hlist(f)%hbuf(:,1,:), varid)
        else
          call cam_grid_write_dist_array(tape(t)%File, fdecomp, adims,        &
               fdims(1:frank), tape(t)%hlist(f)%hbuf, varid)
        end if
      end if
    end do
    !! NACS
    if(restart) then
      if (size(tape(t)%hlist(f)%nacs, 1) > 1) then
        if (nadims > 2) then
          adims(2) = adims(3)
          nadims = 2
        end if
        call cam_grid_dimensions(fdecomp, fdims(1:2), nacsrank)
        call cam_grid_write_dist_array(tape(t)%File, fdecomp, adims(1:nadims),&
             fdims(1:nacsrank), tape(t)%hlist(f)%nacs, tape(t)%hlist(f)%nacs_varid)

        else
          ierr = pio_put_var(tape(t)%File, tape(t)%hlist(f)%nacs_varid,     &
               tape(t)%hlist(f)%nacs(:, tape(t)%hlist(f)%field%begdim3:tape(t)%hlist(f)%field%enddim3))
        end if
      end if

    return
  end subroutine dump_field

  !#######################################################################

  logical function write_inithist ()
    !
    !-----------------------------------------------------------------------
    ! 
    ! Purpose: Set flags that will initiate dump to IC file when OUTFLD and
    ! WSHIST are called
    ! 
    !-----------------------------------------------------------------------
    !
    use time_manager, only: get_nstep, get_curr_date, get_step_size, is_last_step
    !
    ! Local workspace
    !
    integer :: yr, mon, day      ! year, month, and day components of
    ! a date
    integer :: nstep             ! current timestep number
    integer :: ncsec             ! current time of day [seconds]
    integer :: dtime             ! timestep size

    !-----------------------------------------------------------------------

    write_inithist  = .false.

    if(is_initfile()) then

      nstep = get_nstep()
      call get_curr_date(yr, mon, day, ncsec)

      if    (inithist == '6-HOURLY') then
        dtime  = get_step_size()
        write_inithist = nstep /= 0 .and. mod( nstep, nint((6._r8*3600._r8)/dtime) ) == 0
      elseif(inithist == 'DAILY'   ) then
        write_inithist = nstep /= 0 .and. ncsec == 0
      elseif(inithist == 'MONTHLY' ) then
        write_inithist = nstep /= 0 .and. ncsec == 0 .and. day == 1
      elseif(inithist == 'YEARLY'  ) then
        write_inithist = nstep /= 0 .and. ncsec == 0 .and. day == 1 .and. mon == 1
      elseif(inithist == 'CAMIOP'  ) then
        write_inithist = nstep == 0 
      elseif(inithist == 'ENDOFRUN'  ) then
        write_inithist = nstep /= 0 .and. is_last_step()
      end if
    end if

    return
  end function write_inithist

  !#######################################################################

  subroutine wshist (rgnht_in)
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Driver routine to write fields on history tape t
    ! 
    ! Method: For variables which do not need to be gathered (SPMD) just issue the netcdf call
    !         For those that do need to be gathered, call "dump_field" to do the operation.
    !         Finally, zero the history buffers for each field written.
    ! 
    ! Author: CCM Core Group
    ! 
    !-----------------------------------------------------------------------
    use time_manager,  only: get_nstep, get_curr_date, get_curr_time, get_step_size
    use chem_surfvals, only: chem_surfvals_get, chem_surfvals_co2_rad
    use solar_data,    only: sol_tsi
    use sat_hist,      only: sat_hist_write
    use interp_mod,    only: set_interp_hfile
    use cam_pio_utils, only: cam_pio_closefile

    logical, intent(in), optional :: rgnht_in(ptapes)
    !
    ! Local workspace
    !
    character(len=8) :: cdate  ! system date
    character(len=8) :: ctime  ! system time

    logical  :: rgnht(ptapes), restart
    integer t, f               ! tape, field indices
    integer start              ! starting index required by nf_put_vara
    integer count1             ! count values required by nf_put_vara
    integer startc(2)          ! start values required by nf_put_vara (character)
    integer countc(2)          ! count values required by nf_put_vara (character)
#ifdef HDEBUG
    !      integer begdim3
    !      integer enddim3
#endif

    integer :: yr, mon, day      ! year, month, and day components of a date
    integer :: nstep             ! current timestep number
    integer :: ncdate            ! current date in integer format [yyyymmdd]
    integer :: ncsec             ! current time of day [seconds]
    integer :: ndcur             ! day component of current time
    integer :: nscur             ! seconds component of current time
    real(r8) :: time             ! current time
    real(r8) :: tdata(2)         ! time interval boundaries
    character(len=max_string_len) :: fname ! Filename
    logical :: prev              ! Label file with previous date rather than current
    integer :: ierr
#if ( defined BFB_CAM_SCAM_IOP )
    integer :: tsec             ! day component of current time
    integer :: dtime            ! seconds component of current time
#endif
    real(r8) :: kp, ap, f107, f107a

    !-----------------------------------------------------------------------
    ! get solar/geomagnetic activity data...
    call solar_parms_get( f107_s=f107, f107a_s=f107a, ap_s=ap, kp_s=kp )

    if(present(rgnht_in)) then
      rgnht=rgnht_in
      restart=.true.
      tape => restarthistory_tape
    else
      rgnht=.false.
      restart=.false.
      tape => history_tape
    end if

    nstep = get_nstep()
    call get_curr_date(yr, mon, day, ncsec)
    ncdate = yr*10000 + mon*100 + day
    call get_curr_time(ndcur, nscur)
    !
    ! Write time-varying portion of history file header
    !
    do t=1,ptapes
      if (nflds(t) == 0 .or. (restart .and.(.not.rgnht(t)))) cycle
      !
      ! Check if this is the IC file and if it's time to write.
      ! Else, use "nhtfrq" to determine if it's time to write
      ! the other history files.
      !
      if((.not. restart) .or. rgnht(t)) then
        if( is_initfile(file_index=t) ) then
          hstwr(t) =  write_inithist()
          prev     = .false.
        else 
          if (nhtfrq(t) == 0) then
            hstwr(t) = nstep /= 0 .and. day == 1 .and. ncsec == 0
            prev     = .true.
          else
            hstwr(t) = mod(nstep,nhtfrq(t)) == 0
            prev     = .false.
          end if
        end if
      end if
      if (hstwr(t) .or. (restart .and. rgnht(t))) then
        if(masterproc) then
          if(is_initfile(file_index=t)) then
            write(iulog,100) yr,mon,day,ncsec
100         format('WSHIST: writing time sample to Initial Conditions h-file', &
                 ' DATE=',i4.4,'/',i2.2,'/',i2.2,' NCSEC=',i6)
          else if(is_satfile(t)) then
            write(iulog,150) nfils(t),t,yr,mon,day,ncsec
150         format('WSHIST: writing sat columns ',i6,' to h-file ', &
                 i1,' DATE=',i4.4,'/',i2.2,'/',i2.2,' NCSEC=',i6)
          else if(hstwr(t)) then
            write(iulog,200) nfils(t),t,yr,mon,day,ncsec
200         format('WSHIST: writing time sample ',i3,' to h-file ', &
                 i1,' DATE=',i4.4,'/',i2.2,'/',i2.2,' NCSEC=',i6)
          else if(restart .and. rgnht(t)) then
            write(iulog,300) nfils(t),t,yr,mon,day,ncsec
300         format('WSHIST: writing history restart ',i3,' to hr-file ', &
                 i1,' DATE=',i4.4,'/',i2.2,'/',i2.2,' NCSEC=',i6)
          end if
          write(iulog,*)
        end if
        !
        ! Starting a new volume => define the metadata
        !
        if (nfils(t)==0 .or. (restart.and.rgnht(t))) then
          if(restart) then
            rhfilename_spec = '%c.cam' // trim(inst_suffix) // '.rh%t.%y-%m-%d-%s.nc' 
            fname = interpret_filename_spec( rhfilename_spec, number=(t-1))
            hrestpath(t)=fname
          else if(is_initfile(file_index=t)) then
            fname = interpret_filename_spec( hfilename_spec(t) )
          else
            fname = interpret_filename_spec( hfilename_spec(t), number=(t-1), &
                 prev=prev )
          end if
          !
          ! Check that this new filename isn't the same as a previous or current filename
          !
          do f = 1, ptapes
            if (masterproc.and. trim(fname) == trim(nhfil(f)) )then
              write(iulog,*)'WSHIST: New filename same as old file = ', trim(fname)
              write(iulog,*)'Is there an error in your filename specifiers?'
              write(iulog,*)'hfilename_spec(', t, ') = ', hfilename_spec(t)
              if ( t /= f )then
                write(iulog,*)'hfilename_spec(', f, ') = ', hfilename_spec(f)
              end if
              call endrun
            end if
          end do
          if(.not. restart) then
            nhfil(t) = fname
            if(masterproc) write(iulog,*)'WSHIST: nhfil(',t,')=',trim(nhfil(t))
            cpath(t) = nhfil(t)
            if ( len_trim(nfpath(t)) == 0 ) nfpath(t) = cpath(t)
          end if
          call h_define (t, restart)
        end if

        if(is_satfile(t)) then
          call sat_hist_write( tape(t), nflds(t), nfils(t))
        else
          if(restart) then
            start=1
          else
            nfils(t) = nfils(t) + 1
            start = nfils(t)
          end if
          count1 = 1
          ! Setup interpolation data if history file is interpolated
          if (interpolate_output(t) .and. (.not. restart)) then
            call set_interp_hfile(t, interpolate_info)
          end if

          ierr = pio_put_var (tape(t)%File, tape(t)%ndcurid,(/start/), (/count1/),(/ndcur/))
          ierr = pio_put_var (tape(t)%File, tape(t)%nscurid,(/start/), (/count1/),(/nscur/))
          ierr = pio_put_var (tape(t)%File, tape(t)%dateid,(/start/), (/count1/),(/ncdate/))

          if (.not. is_initfile(file_index=t)) then
            ! Don't write the GHG/Solar forcing data to the IC file.
            ierr=pio_put_var (tape(t)%File, tape(t)%co2vmrid,(/start/), (/count1/),(/chem_surfvals_co2_rad(vmr_in=.true.)/))
            ierr=pio_put_var (tape(t)%File, tape(t)%ch4vmrid,(/start/), (/count1/),(/chem_surfvals_get('CH4VMR')/))
            ierr=pio_put_var (tape(t)%File, tape(t)%n2ovmrid,(/start/), (/count1/),(/chem_surfvals_get('N2OVMR')/))
            ierr=pio_put_var (tape(t)%File, tape(t)%f11vmrid,(/start/), (/count1/),(/chem_surfvals_get('F11VMR')/))
            ierr=pio_put_var (tape(t)%File, tape(t)%f12vmrid,(/start/), (/count1/),(/chem_surfvals_get('F12VMR')/))
            ierr=pio_put_var (tape(t)%File, tape(t)%sol_tsiid,(/start/), (/count1/),(/sol_tsi/))

            if (solar_parms_on) then
              ierr=pio_put_var (tape(t)%File, tape(t)%f107id, (/start/), (/count1/),(/ f107 /) )
              ierr=pio_put_var (tape(t)%File, tape(t)%f107aid,(/start/), (/count1/),(/ f107a /) )
              ierr=pio_put_var (tape(t)%File, tape(t)%kpid,   (/start/), (/count1/),(/ kp /) )
              ierr=pio_put_var (tape(t)%File, tape(t)%apid,   (/start/), (/count1/),(/ ap /) )
            endif

          end if

          ierr = pio_put_var (tape(t)%File, tape(t)%datesecid,(/start/),(/count1/),(/ncsec/))
#if ( defined BFB_CAM_SCAM_IOP )
          dtime = get_step_size()
          tsec=dtime*nstep
          ierr = pio_put_var (tape(t)%File, tape(t)%tsecid,(/start/),(/count1/),(/tsec/))
#endif
          ierr = pio_put_var (tape(t)%File, tape(t)%nstephid,(/start/),(/count1/),(/nstep/))
          time = ndcur + nscur/86400._r8
          ierr=pio_put_var (tape(t)%File, tape(t)%timeid, (/start/),(/count1/),(/time/))

          startc(1) = 1
          startc(2) = start
          countc(1) = 2
          countc(2) = 1
          if (is_initfile(file_index=t)) then
            tdata = time   ! Inithist file is always instantanious data
          else
            tdata(1) = beg_time(t)
            tdata(2) = time
          end if
          ierr=pio_put_var (tape(t)%File, tape(t)%tbndid, startc, countc, tdata)
          if(.not.restart) beg_time(t) = time  ! update beginning time of next interval
          startc(1) = 1
          startc(2) = start
          countc(1) = 8
          countc(2) = 1
          call datetime (cdate, ctime)
          ierr = pio_put_var (tape(t)%File, tape(t)%date_writtenid, startc, countc, (/cdate/))
          ierr = pio_put_var (tape(t)%File, tape(t)%time_writtenid, startc, countc, (/ctime/))

          if(.not. restart) then
            !$OMP PARALLEL DO PRIVATE (F)
            do f=1,nflds(t)
              ! Normalized averaged fields
              if (tape(t)%hlist(f)%avgflag /= 'I') then
                call h_normalize (f, t)
              end if
            end do
          end if
          !
          ! Write field to history tape.  Note that this is NOT threaded due to netcdf limitations
          !
          call t_startf ('dump_field')
          do f=1,nflds(t)
            call dump_field(f, t, restart)
          end do
          call t_stopf ('dump_field')
          !
          ! Zero history buffers and accumulators now that the fields have been written.
          !



          if(restart) then
            do f=1,nflds(t)
              if(associated(tape(t)%hlist(f)%varid)) then
                deallocate(tape(t)%hlist(f)%varid)
                nullify(tape(t)%hlist(f)%varid)
              end if
            end do
            call cam_pio_closefile(tape(t)%File)
          else
            !$OMP PARALLEL DO PRIVATE (F)
            do f=1,nflds(t)
              call h_zero (f, t)
            end do
          end if
        end if
      end if

    end do

    return
  end subroutine wshist

  !#######################################################################

  subroutine addfld_1d(fname, vdim_name, avgflag, units, long_name,           &
       gridname, flag_xyfill, sampling_seq, standard_name, fill_value)

    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Add a field to the master field list
    ! 
    ! Method: Put input arguments of field name, units, number of levels,
    !         averaging flag, and long name into a type entry in the global
    !         master field list (masterlist).
    ! 
    !-----------------------------------------------------------------------

    !
    ! Arguments
    !
    character(len=*), intent(in)  :: fname      ! field name (max_fieldname_len)
    character(len=*), intent(in)  :: vdim_name  ! NetCDF dimension name (or scalar coordinate)
    character(len=1), intent(in)  :: avgflag    ! averaging flag
    character(len=*), intent(in)  :: units      ! units of fname (max_chars)
    character(len=*), intent(in)  :: long_name  ! long name of field (max_chars)

    character(len=*), intent(in), optional :: gridname    ! decomposition type
    logical, intent(in), optional :: flag_xyfill ! non-applicable xy points flagged with fillvalue
    character(len=*), intent(in), optional :: sampling_seq ! sampling sequence - if not every timestep, 
    ! how often field is sampled:  
    ! every other; only during LW/SW radiation calcs, etc.
    character(len=*), intent(in), optional :: standard_name  ! CF standard name (max_chars)
    real(r8),         intent(in), optional :: fill_value

    !
    ! Local workspace
    !
    character(len=max_chars), allocatable :: dimnames(:)
    integer                               :: index

    if (trim(vdim_name) == trim(horiz_only)) then
      allocate(dimnames(0))
    else
      index = get_hist_coord_index(trim(vdim_name))
      if (index < 1) then
        call endrun('ADDFLD: Invalid coordinate, '//trim(vdim_name))
      end if
      allocate(dimnames(1))
      dimnames(1) = trim(vdim_name)
    end if
    call addfld(fname, dimnames, avgflag, units, long_name, gridname,         &
         flag_xyfill, sampling_seq, standard_name, fill_value)

  end subroutine addfld_1d

  subroutine addfld_nd(fname, dimnames, avgflag, units, long_name,            &
       gridname, flag_xyfill, sampling_seq, standard_name, fill_value)

    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Add a field to the master field list
    ! 
    ! Method: Put input arguments of field name, units, number of levels,
    !         averaging flag, and long name into a type entry in the global
    !         master field list (masterlist).
    ! 
    !-----------------------------------------------------------------------
    use cam_history_support, only: fillvalue, hist_coord_find_levels
    use cam_grid_support,    only: cam_grid_id

    !
    ! Arguments
    !
    character(len=*), intent(in)  :: fname      ! field name (max_fieldname_len)
    character(len=*), intent(in)  :: dimnames(:) ! NetCDF dimension names (except grid dims)
    character(len=1), intent(in)  :: avgflag    ! averaging flag
    character(len=*), intent(in)  :: units      ! units of fname (max_chars)
    character(len=*), intent(in)  :: long_name  ! long name of field (max_chars)

    character(len=*), intent(in), optional :: gridname    ! decomposition type
    logical, intent(in), optional :: flag_xyfill ! non-applicable xy points flagged with fillvalue
    character(len=*), intent(in), optional :: sampling_seq ! sampling sequence - if not every timestep, 
    ! how often field is sampled:  
    ! every other; only during LW/SW radiation calcs, etc.
    character(len=*), intent(in), optional :: standard_name  ! CF standard name (max_chars)
    real(r8),         intent(in), optional :: fill_value

    !
    ! Local workspace
    !
    character(len=max_fieldname_len) :: fname_tmp ! local copy of fname
    character(len=128)               :: errormsg
    type(master_entry), pointer      :: listentry

    integer :: dimcnt

    if (htapes_defined) then
      call endrun ('ADDFLD: Attempt to add field '//trim(fname)//' after history files set')
    end if

    !
    ! Ensure that new field name is not all blanks
    !
    if (len_trim(fname)==0) then
      call endrun('ADDFLD: blank field name not allowed')
    end if
    !
    ! Ensure that new field name is not longer than allowed
    ! (strip "&IC" suffix if it exists)
    !
    fname_tmp  = fname
    fname_tmp  = strip_suffix(fname_tmp)

    if (len_trim(fname_tmp) > fieldname_len) then
      write(iulog,*)'ADDFLD: field name cannot be longer than ',fieldname_len,' characters long'
      write(iulog,*)'Field name:  ',fname
      write(errormsg, *) 'Field name, "', trim(fname), '" is too long'
      call endrun('ADDFLD: '//trim(errormsg))
    end if
    !
    ! Ensure that new field doesn't already exist
    !
    listentry => get_entry_by_name(masterlinkedlist, fname)
    if(associated(listentry)) then
      call endrun ('ADDFLD:  '//fname//' already on list')
    end if

    !
    ! Add field to Master Field List arrays fieldn and iflds
    !
    allocate(listentry)
    listentry%field%name        = fname
    listentry%field%long_name   = long_name
    listentry%field%numlev      = 1        ! Will change if lev or ilev in shape
    listentry%field%units       = units
    listentry%field%meridional_complement = -1
    listentry%field%zonal_complement      = -1
    listentry%htapeindx(:) = -1
    listentry%act_sometape = .false.
    listentry%actflag(:) = .false.

    ! Make sure we have a valid gridname
    if (present(gridname)) then
      listentry%field%decomp_type = cam_grid_id(trim(gridname))
    else
      listentry%field%decomp_type = cam_grid_id('physgrid')
    end if
    if (listentry%field%decomp_type < 0) then
      write(errormsg, *) 'Invalid grid name, "', trim(gridname), '" for ',    &
           trim(fname)
      call endrun('ADDFLD: '//trim(errormsg))
    end if

    !
    ! Indicate sampling sequence of field (i.e., how often "outfld" is called)
    ! If not every timestep (default), then give a descriptor indicating the
    ! sampling pattern.  Currently, the only valid value is "rad_lwsw" for sampling
    ! during LW/SW radiation timesteps only
    !
    if (present(sampling_seq)) then
      listentry%field%sampling_seq = sampling_seq
    else
      listentry%field%sampling_seq = ' '
    end if
    !
    ! Whether to apply xy fillvalue: default is false
    !
    if (present(flag_xyfill)) then
      listentry%field%flag_xyfill = flag_xyfill
    else
      listentry%field%flag_xyfill = .false.
    end if

    !
    !    Allow external packages to have fillvalues different than default
    !

    if(present(fill_value)) then
      listentry%field%fillvalue = fill_value
    else
      listentry%field%fillvalue = fillvalue
    endif

    !
    ! Process shape
    !

    if (associated(listentry%field%mdims)) then
      deallocate(listentry%field%mdims)
    end if
    nullify(listentry%field%mdims)
    dimcnt = size(dimnames)
    allocate(listentry%field%mdims(dimcnt))
    call lookup_hist_coord_indices(dimnames, listentry%field%mdims)
    if(dimcnt > maxvarmdims) then
      maxvarmdims = dimcnt
    end if
    ! Check for subcols (currently limited to first dimension)
    listentry%field%is_subcol = .false.
    if (size(dimnames) > 0) then
      if (trim(dimnames(1)) == 'psubcols') then
        if (listentry%field%decomp_type /= cam_grid_id('physgrid')) then
          write(errormsg, *) "Cannot add ", trim(fname),                      &
               "Subcolumn history output only allowed on physgrid"
          call endrun("ADDFLD: "//errormsg)
          listentry%field%is_subcol = .true.
        end if
      end if
    end if
    ! Levels
    listentry%field%numlev = hist_coord_find_levels(dimnames)
    if (listentry%field%numlev <= 0) then
      listentry%field%numlev = 1
    end if

    !
    ! Dimension history info based on decomposition type (grid)
    !
    call set_field_dimensions(listentry%field)

    !
    ! These 2 fields are used only in master field list, not runtime field list
    !
    listentry%avgflag(:) = avgflag
    listentry%actflag(:) = .false.

    do dimcnt = 1, ptapes
      call AvgflagToString(avgflag, listentry%time_op(dimcnt))
    end do

    nullify(listentry%next_entry)

    call add_entry_to_master(listentry)
    return
  end subroutine addfld_nd

  !#######################################################################

  ! field_part_of_vector: Determinie if fname is part of a vector set
  !       Optionally fill in the names of the vector set fields
  logical function field_part_of_vector(fname, meridional_name, zonal_name)

    ! Dummy arguments
    character(len=*),           intent(in)  :: fname
    character(len=*), optional, intent(out) :: meridional_name
    character(len=*), optional, intent(out) :: zonal_name

    ! Local variables
    type(master_entry), pointer             :: listentry

    listentry => get_entry_by_name(masterlinkedlist, fname)
    if (associated(listentry)) then
      if ( (len_trim(listentry%meridional_field) > 0) .or.                     &
           (len_trim(listentry%zonal_field) > 0)) then
        field_part_of_vector = .true.
        if (present(meridional_name)) then
          meridional_name = listentry%meridional_field
        end if
        if (present(zonal_name)) then
          zonal_name = listentry%zonal_field
        end if
      else
        field_part_of_vector = .false.
      end if
    else
      field_part_of_vector = .false.
    end if
    if (.not. field_part_of_vector) then
      if (present(meridional_name)) then
        meridional_name = ''
      end if
      if (present(zonal_name)) then
        zonal_name = ''
      end if
    end if

  end function field_part_of_vector


  ! register_vector_field: Register a pair of history field names as
  !           being a vector complement set.
  !           This information is used to set up interpolated history output.
  ! NB: register_vector_field must be called after both fields are defined
  !     with addfld
  subroutine register_vector_field(zonal_field_name, meridional_field_name)

    ! Dummy arguments
    character(len=*),             intent(in) :: zonal_field_name
    character(len=*),             intent(in) :: meridional_field_name

    ! Local variables
    type(master_entry), pointer      :: mlistentry
    type(master_entry), pointer      :: zlistentry
    character(len=*),   parameter    :: subname = 'REGISTER_VECTOR_FIELD'
    character(len=max_chars)         :: errormsg
    
    if (htapes_defined) then
      write(errormsg, '(5a)') ': Attempt to register vector field (',         &
           trim(zonal_field_name), ', ', trim(meridional_field_name),         &
           ') after history files set'
      call endrun (trim(subname)//errormsg)
    end if

    ! Look for the field IDs
    zlistentry => get_entry_by_name(masterlinkedlist, zonal_field_name)
    mlistentry => get_entry_by_name(masterlinkedlist, meridional_field_name)
    ! Has either of these fields been previously registered?
    if (associated(mlistentry)) then
      if (len_trim(mlistentry%meridional_field) > 0) then
        write(errormsg, '(9a)') ': ERROR attempting to register vector ',     &
             'field (', trim(zonal_field_name), ', ',                         &
             trim(meridional_field_name), '), ', trim(meridional_field_name), &
             ' has been registered as part of a vector field with ',          &
             trim(mlistentry%meridional_field)
        call endrun (trim(subname)//errormsg)
      else if (len_trim(mlistentry%zonal_field) > 0) then
        write(errormsg, '(9a)') ': ERROR attempting to register vector ',     &
             'field (', trim(zonal_field_name), ', ',                         &
             trim(meridional_field_name), '), ', trim(meridional_field_name), &
             ' has been registered as part of a vector field with ',          &
             trim(mlistentry%zonal_field)
        call endrun (trim(subname)//errormsg)
      end if
    end if
    if (associated(zlistentry)) then
      if (len_trim(zlistentry%meridional_field) > 0) then
        write(errormsg, '(9a)') ': ERROR attempting to register vector ',     &
             'field (', trim(zonal_field_name), ', ',                         &
             trim(meridional_field_name), '), ', trim(zonal_field_name),      &
             ' has been registered as part of a vector field with ',          &
             trim(zlistentry%meridional_field)
        call endrun (trim(subname)//errormsg)
      else if (len_trim(zlistentry%zonal_field) > 0) then
        write(errormsg, '(9a)') ': ERROR attempting to register vector ',     &
             'field (', trim(zonal_field_name), ', ',                         &
             trim(meridional_field_name), '), ', trim(zonal_field_name),      &
             ' has been registered as part of a vector field with ',          &
             trim(zlistentry%meridional_field)
        call endrun (trim(subname)//errormsg)
      end if
    end if
    if(associated(mlistentry) .and. associated(zlistentry)) then
      zlistentry%meridional_field = mlistentry%field%name
      zlistentry%zonal_field      = ''
      mlistentry%meridional_field = ''
      mlistentry%zonal_field      = zlistentry%field%name
    else if (associated(mlistentry)) then
      write(errormsg, '(7a)') ': ERROR attempting to register vector field (',&
           trim(zonal_field_name), ', ', trim(meridional_field_name),         &
           '), ', trim(zonal_field_name), ' is not defined'
      call endrun (trim(subname)//errormsg)
    else if (associated(zlistentry)) then
      write(errormsg, '(7a)') ': ERROR attempting to register vector field (',&
           trim(zonal_field_name), ', ', trim(meridional_field_name),         &
           '), ', trim(meridional_field_name), ' is not defined'
      call endrun (trim(subname)//errormsg)
    else
      write(errormsg, '(5a)') ': ERROR attempting to register vector field (',&
           trim(zonal_field_name), ', ', trim(meridional_field_name),         &
           '), neither field is defined'
      call endrun (trim(subname)//errormsg)
    end if
  end subroutine register_vector_field

  subroutine add_entry_to_master( newentry)
    type(master_entry), target, intent(in) :: newentry
    type(master_entry), pointer :: listentry

    if(associated(masterlinkedlist)) then
      listentry => masterlinkedlist
      do while(associated(listentry%next_entry))
        listentry=>listentry%next_entry
      end do
      listentry%next_entry=>newentry
    else
      masterlinkedlist=>newentry
    end if

  end subroutine add_entry_to_master

  !#######################################################################

  subroutine wrapup (rstwr, nlend)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose: 
    ! Close history files.
    ! 
    ! Method: 
    ! This routine will close any full hist. files
    ! or any hist. file that has data on it when restart files are being 
    ! written.
    ! If a partially full history file was disposed (for restart 
    ! purposes), then wrapup will open that unit back up and position 
    ! it for appending new data. 
    !
    ! Original version: CCM2
    !
    !-----------------------------------------------------------------------
    !
    use pio, only : pio_file_is_open
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use pspect
    use ioFileMod
    use time_manager,  only: get_nstep, get_curr_date, get_curr_time
    use cam_pio_utils, only: cam_pio_openfile, cam_pio_closefile

    !
    ! Input arguments
    !
    logical, intent(in) :: rstwr   ! true => restart files are written this timestep
    logical, intent(in) :: nlend   ! Flag if time to end

    !
    ! Local workspace
    !
    integer  :: nstep            ! current timestep number
    integer  :: ncdate           ! current date in integer format [yyyymmdd]
    integer  :: ncsec            ! time of day relative to current date [secs]
    integer  :: ndcur            ! days component of current time
    integer  :: nscur            ! seconds component of current time
    integer  :: yr, mon, day     ! year, month, day components of a date

    logical  :: lfill   (ptapes) ! Is history file ready to dispose?
    logical  :: lhdisp           ! true => history file is disposed
    logical  :: lhfill           ! true => history file is full

    integer  :: t                ! History file number
    integer  :: f
    real(r8) :: tday             ! Model day number for printout
    !-----------------------------------------------------------------------

    tape => history_tape

    nstep = get_nstep()
    call get_curr_date(yr, mon, day, ncsec)
    ncdate = yr*10000 + mon*100 + day
    call get_curr_time(ndcur, nscur)
    !
    !-----------------------------------------------------------------------
    ! Dispose history files.
    !-----------------------------------------------------------------------
    !
    ! Begin loop over ptapes (the no. of declared history files - primary
    ! and auxiliary).  This loop disposes a history file to Mass Store
    ! when appropriate.
    !
    do t=1,ptapes
      if (nflds(t) == 0) cycle

      lfill(t) = .false.
      !
      ! Find out if file is full
      !
      if (hstwr(t) .and. nfils(t) >= mfilt(t)) then
        lfill(t) = .true.
      endif
      !
      ! Dispose history file if 
      !    1) file is filled or 
      !    2) this is the end of run and file has data on it or
      !    3) restarts are being put out and history file has data on it
      !
      if (lfill(t) .or. (nlend .and. nfils(t) >= 1) .or. (rstwr .and. nfils(t) >= 1)) then
        !
        ! Dispose history file
        !
        !
        ! Is this the 0 timestep data of a monthly run?
        ! If so, just close primary unit do not dispose. 
        !
        if (masterproc) write(iulog,*)'WRAPUP: nf_close(',t,')=',trim(nhfil(t))
        if(pio_file_is_open(tape(t)%File)) then
          if (nlend .or. lfill(t)) then
            do f=1,nflds(t)
              if (associated(tape(t)%hlist(f)%varid)) then
                deallocate(tape(t)%hlist(f)%varid)
                nullify(tape(t)%hlist(f)%varid)
              end if
            end do
          end if
          call cam_pio_closefile(tape(t)%File)
        end if
        if (nhtfrq(t) /= 0 .or. nstep > 0) then

          ! 
          ! Print information concerning model output.
          ! Model day number = iteration number of history file data * delta-t / (seconds per day)
          ! 
          tday = ndcur + nscur/86400._r8
          if(masterproc) then
            if (t==1) then
              write(iulog,*)'   Primary history file'
            else
              write(iulog,*)'   Auxiliary history file number ', t-1
            end if
            write(iulog,9003)nstep,nfils(t),tday
            write(iulog,9004)
          end if
          !                      
          ! Auxilary files may have been closed and saved off without being full. 
          ! We must reopen the files and position them for more data.
          ! Must position auxiliary files if not full
          !              
          if (.not.nlend .and. .not.lfill(t)) then
            call cam_PIO_openfile (tape(t)%File, nhfil(t), PIO_WRITE)
            call h_inquire(t)
          end if
        endif                 ! if 0 timestep of montly run****
      end if                      ! if time dispose history fiels***   
    end do                         ! do ptapes
    !
    ! Reset number of files on each history tape
    !
    do t=1,ptapes
      if (nflds(t) == 0) cycle
      lhfill = hstwr(t) .and. nfils(t) >= mfilt(t)
      lhdisp = lhfill .or. (nlend .and. nfils(t) >= 1) .or. &
           (rstwr .and. nfils(t) >= 1)
      if (lhfill.and.lhdisp) then
        nfils(t) = 0
      endif
    end do
    return
9003 format('    Output at NSTEP     = ',i10,/, &
         '    Number of time samples on this file = ',i10,/, &
         '    Model Day           = ',f10.2)
9004 format('---------------------------------------')
  end subroutine wrapup


  integer function gen_hash_key(string)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose: Generate a hash key on the interval [0 .. tbl_hash_pri_sz-1]
    !          given a character string.
    !
    ! Algorithm is a variant of perl's internal hashing function.
    !
    !-----------------------------------------------------------------------
    !
    implicit none
    !
    !  Arguments:
    !
    character(len=*), intent(in) :: string
    !
    !  Local.
    !
    integer :: hash
    integer :: i

    hash = gen_hash_key_offset

    if ( len(string) /= 19 ) then
      !
      !     Process arbitrary string length.
      !
      do i = 1, len(string)
        hash = ieor(hash, (ichar(string(i:i)) * tbl_gen_hash_key(iand(i-1,tbl_max_idx))))
      end do
    else
      !
      !     Special case string length = 19
      !
      hash = ieor(hash , ichar(string(1:1))   * 61)
      hash = ieor(hash , ichar(string(2:2))   * 59)
      hash = ieor(hash , ichar(string(3:3))   * 53)
      hash = ieor(hash , ichar(string(4:4))   * 47)
      hash = ieor(hash , ichar(string(5:5))   * 43)
      hash = ieor(hash , ichar(string(6:6))   * 41)
      hash = ieor(hash , ichar(string(7:7))   * 37)
      hash = ieor(hash , ichar(string(8:8))   * 31)
      hash = ieor(hash , ichar(string(9:9))   * 29)
      hash = ieor(hash , ichar(string(10:10)) * 23)
      hash = ieor(hash , ichar(string(11:11)) * 17)
      hash = ieor(hash , ichar(string(12:12)) * 13)
      hash = ieor(hash , ichar(string(13:13)) * 11)
      hash = ieor(hash , ichar(string(14:14)) * 7)
      hash = ieor(hash , ichar(string(15:15)) * 3)
      hash = ieor(hash , ichar(string(16:16)) * 1)
      hash = ieor(hash , ichar(string(17:17)) * 61)
      hash = ieor(hash , ichar(string(18:18)) * 59)
      hash = ieor(hash , ichar(string(19:19)) * 53)
    end if

    gen_hash_key = iand(hash, tbl_hash_pri_sz-1)

    return

  end function gen_hash_key

  !#######################################################################

  integer function get_masterlist_indx(fldname)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose: Return the the index of the field's name on the master file list.
    !
    !          If the field is not found on the masterlist, return -1.
    !
    !-----------------------------------------------------------------------
    !
    !  Arguments:
    !
    character(len=*), intent(in) :: fldname
    !
    !  Local.
    !
    integer :: hash_key
    integer :: ff
    integer :: ii
    integer :: io   ! Index of overflow chain in overflow table
    integer :: in   ! Number of entries on overflow chain

    hash_key = gen_hash_key(fldname)
    ff = tbl_hash_pri(hash_key)
    if ( ff < 0 ) then
      io = abs(ff)
      in = tbl_hash_oflow(io)
      do ii = 1, in
        ff = tbl_hash_oflow(io+ii)
        if ( masterlist(ff)%thisentry%field%name == fldname ) exit
      end do
    end if

    if (ff == 0) then
      ! fldname generated a hash key that doesn't have an entry in tbl_hash_pri.
      ! This means that fldname isn't in the masterlist
      call endrun ('GET_MASTERLIST_INDX: attemping to output field '//fldname//' not on master list')
    end if

    if (associated(masterlist(ff)%thisentry) .and. masterlist(ff)%thisentry%field%name /= fldname ) then
      call endrun ('GET_MASTERLIST_INDX: error finding field '//fldname//' on master list')
    end if

    get_masterlist_indx = ff
    return
  end function get_masterlist_indx
  !#######################################################################

  subroutine bld_outfld_hash_tbls()
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose: Build primary and overflow hash tables for outfld processing.
    !
    ! Steps:
    !  1) Foreach field on masterlist, find all collisions.
    !  2) Given the number of collisions, verify overflow table has sufficient
    !     space.
    !  3) Build primary and overflow indices.
    !
    !-----------------------------------------------------------------------
    !
    !  Local.
    !
    integer :: ff
    integer :: ii
    integer :: itemp
    integer :: ncollisions
    integer :: hash_key
    type(master_entry), pointer :: listentry
    !
    !  1) Find all collisions.
    !
    tbl_hash_pri = 0

    ff=0
    allocate(masterlist(nfmaster))
    listentry=>masterlinkedlist
    do while(associated(listentry))
      ff=ff+1
      masterlist(ff)%thisentry=>listentry
      listentry=>listentry%next_entry
    end do
    if(ff /= nfmaster) then
      write(iulog,*) 'nfmaster = ',nfmaster, ' ff=',ff
      call endrun('mismatch in expected size of nfmaster')
    end if


    do ff = 1, nfmaster
      hash_key = gen_hash_key(masterlist(ff)%thisentry%field%name)
      tbl_hash_pri(hash_key) = tbl_hash_pri(hash_key) + 1
    end do

    !
    !  2) Count number of collisions and define start of a individual
    !     collision's chain in overflow table. A collision is defined to be any
    !     location in tbl_hash_pri that has a value > 1.
    !
    ncollisions = 0
    do ii = 0, tbl_hash_pri_sz-1
      if ( tbl_hash_pri(ii) > 1 ) then  ! Define start of chain in O.F. table
        itemp = tbl_hash_pri(ii)
        tbl_hash_pri(ii) = -(ncollisions + 1)
        ncollisions = ncollisions + itemp + 1
      end if
    end do

    if ( ncollisions > tbl_hash_oflow_sz ) then
      write(iulog,*) 'BLD_OUTFLD_HASH_TBLS: ncollisions > tbl_hash_oflow_sz', &
           ncollisions, tbl_hash_oflow_sz
      call endrun()
    end if

    !
    !  3) Build primary and overflow tables.
    !     i - set collisions in tbl_hash_pri to point to their respective
    !         chain in the overflow table.
    !
    tbl_hash_oflow = 0

    do ff = 1, nfmaster
      hash_key = gen_hash_key(masterlist(ff)%thisentry%field%name)
      if ( tbl_hash_pri(hash_key) < 0 ) then
        ii = abs(tbl_hash_pri(hash_key))
        tbl_hash_oflow(ii) = tbl_hash_oflow(ii) + 1
        tbl_hash_oflow(ii+tbl_hash_oflow(ii)) = ff
      else
        tbl_hash_pri(hash_key) = ff
      end if
    end do

    !
    !  Dump out primary and overflow hashing tables.
    !
    !   if ( masterproc ) then
    !      do ii = 0, tbl_hash_pri_sz-1
    !         if ( tbl_hash_pri(ii) /= 0 ) write(iulog,666) 'tbl_hash_pri', ii, tbl_hash_pri(ii)
    !      end do
    !
    !      do ii = 1, tbl_hash_oflow_sz
    !         if ( tbl_hash_oflow(ii) /= 0 ) write(iulog,666) 'tbl_hash_oflow', ii, tbl_hash_oflow(ii)
    !      end do
    !
    !      itemp = 0
    !      ii = 1
    !      do 
    !         if ( tbl_hash_oflow(ii) == 0 ) exit
    !         itemp = itemp + 1
    !         write(iulog,*) 'Overflow chain ', itemp, ' has ', tbl_hash_oflow(ii), ' entries:'
    !         do ff = 1, tbl_hash_oflow(ii)  ! dump out colliding names on this chain
    !            write(iulog,*) '     ', ff, ' = ', tbl_hash_oflow(ii+ff), &
    !                       ' ', masterlist(tbl_hash_oflow(ii+ff))%thisentry%field%name
    !         end do
    !         ii = ii + tbl_hash_oflow(ii) +1 !advance pointer to start of next chain
    !      end do
    !   end if

    return
666 format(1x, a, '(', i4, ')', 1x, i6)

  end subroutine bld_outfld_hash_tbls

  !#######################################################################

  subroutine bld_htapefld_indices
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose: Set history tape field indicies in masterlist for each
    !          field defined on every tape.
    !
    ! Note: because of restart processing, the actflag field is cleared and
    !       then set only for active output fields on the different history
    !       tapes.
    !
    !-----------------------------------------------------------------------
    !
    !  Arguments:
    !

    !
    !  Local.
    !
    integer :: f
    integer :: t

    !
    !  Initialize htapeindx to an invalid value.
    !
    type(master_entry), pointer :: listentry

    ! reset all the active flags to false 
    ! this is needed so that restarts work properly -- fvitt
    listentry=>masterlinkedlist
    do while(associated(listentry))
      listentry%actflag(:) = .false.
      listentry%act_sometape = .false.
      listentry=>listentry%next_entry
    end do

    do t = 1, ptapes
      do f = 1, nflds(t)
        listentry => get_entry_by_name(masterlinkedlist, tape(t)%hlist(f)%field%name)
        if(.not.associated(listentry)) then
          write(iulog,*) 'BLD_HTAPEFLD_INDICES: something wrong, field not found on masterlist'
          write(iulog,*) 'BLD_HTAPEFLD_INDICES: t, f, ff = ', t, f
          write(iulog,*) 'BLD_HTAPEFLD_INDICES: tape%name = ', tape(t)%hlist(f)%field%name
          call endrun
        end if
        listentry%act_sometape = .true.
        listentry%actflag(t) = .true.
        listentry%htapeindx(t) = f
      end do
    end do

    !
    ! set flag indicating h-tape contents are now defined (needed by addfld)
    !
    htapes_defined = .true.      

    return
  end subroutine bld_htapefld_indices

  !#######################################################################

  logical function hist_fld_active(fname)
    !
    !------------------------------------------------------------------------
    !
    ! Purpose: determine if a field is active on any history file
    !
    !------------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*), intent(in) :: fname ! Field name
    !
    ! Local variables
    !
    character*(max_fieldname_len) :: fname_loc  ! max-char equivalent of fname
    integer :: ff                  ! masterlist index pointer
    !-----------------------------------------------------------------------

    fname_loc = fname
    ff = get_masterlist_indx(fname_loc)
    if ( ff < 0 ) then
      hist_fld_active = .false.
    else 
      hist_fld_active = masterlist(ff)%thisentry%act_sometape
    end if

  end function hist_fld_active

  !#######################################################################

  function hist_fld_col_active(fname, lchnk, numcols)
    use cam_history_support, only: history_patch_t

    ! Determine whether each column in a field is active on any history file.
    ! The purpose of this routine is to provide information which would allow
    ! a diagnostic physics parameterization to only be run on a subset of
    ! columns in the case when only column or regional output is requested.
    !
    ! **N.B.** The field is assumed to be using the physics decomposition.

    ! Arguments
    character(len=*), intent(in) :: fname   ! Field name
    integer,          intent(in) :: lchnk   ! chunk ID
    integer,          intent(in) :: numcols ! Size of return array

    ! Return value
    logical :: hist_fld_col_active(numcols)

    ! Local variables
    integer                         :: ff          ! masterlist index pointer
    integer                         :: i
    integer                         :: t           ! history file (tape) index
    integer                         :: f           ! field index
    integer                         :: decomp
    logical                         :: activeloc(numcols)
    integer                         :: num_patches
    logical                         :: patch_output
    logical                         :: found
    type(history_patch_t), pointer  :: patchptr

    type (active_entry),   pointer  :: tape(:)

    !-----------------------------------------------------------------------

    ! Initialize to false.  Then look to see if and where active.
    hist_fld_col_active = .false.

    ! Check for name in the master list.
    call get_field_properties(fname, found, tape_out=tape, ff_out=ff)

    ! If not in master list then return.
    if (.not. found) return

    ! If in master list, but not active on any file then return
    if (.not. masterlist(ff)%thisentry%act_sometape) return

    ! Loop over history files and check for the field/column in each one
    do t = 1, ptapes

      ! Is the field active in this file?  If not the cycle to next file.
      if (.not. masterlist(ff)%thisentry%actflag(t)) cycle

      f = masterlist(ff)%thisentry%htapeindx(t)
      decomp = tape(t)%hlist(f)%field%decomp_type
      patch_output = associated(tape(t)%patches)

      ! Check whether this file has patch (column) output.
      if (patch_output) then
        num_patches = size(tape(t)%patches)

        do i = 1, num_patches
          patchptr => tape(t)%patches(i)
          activeloc = .false.
          call patchptr%active_cols(decomp, lchnk, activeloc)
          hist_fld_col_active = hist_fld_col_active .or. activeloc
        end do
      else

        ! No column output has been requested.  In that case the field has
        ! global output which implies all columns are active.  No need to 
        ! check any other history files.
        hist_fld_col_active = .true.
        exit

      end if

    end do ! history files

  end function hist_fld_col_active

end module cam_history
