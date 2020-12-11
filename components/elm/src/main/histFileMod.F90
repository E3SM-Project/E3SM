module histFileMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing methods to for CLM history file handling.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_sys_mod    , only : shr_sys_flush
  use spmdMod        , only : masterproc
  use abortutils     , only : endrun
  use elm_varctl     , only : iulog, use_vertsoilc, use_fates
  use elm_varcon     , only : spval, ispval, dzsoi_decomp 
  use elm_varcon     , only : grlnd, nameg, namet, namel, namec, namep
  use decompMod      , only : get_proc_bounds, get_proc_global, bounds_type
  use GridcellType   , only : grc_pp                
  use LandunitType   , only : lun_pp                
  use ColumnType     , only : col_pp                
  use VegetationType , only : veg_pp                
  use ncdio_pio 
  use EDTypesMod        , only : nclmax_fates     => nclmax
  use EDTypesMod        , only : nlevleaf_fates   => nlevleaf
  use FatesInterfaceTypesMod , only : nlevsclass_fates => nlevsclass
  use FatesInterfaceTypesMod , only : nlevage_fates    => nlevage
  use FatesInterfaceTypesMod , only : nlevheight_fates => nlevheight
  use FatesInterfaceTypesMod , only : nlevcoage
  use EDTypesMod        , only : nfsc_fates       => nfsc
  use FatesLitterMod    , only : ncwd_fates       => ncwd
  use FatesInterfaceTypesMod , only : numpft_fates     => numpft
  use EDTypesMod        , only : nelements_fates  => num_elements

  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  !
  ! Constants
  !
  integer , public, parameter :: max_tapes = 6          ! max number of history tapes
  integer , public, parameter :: max_flds = 2500        ! max number of history fields
  integer , public, parameter :: max_namlen = 64        ! maximum number of characters for field name
  integer , private, parameter :: hist_dim_name_length = 16 ! lenngth of character strings in dimension names

  ! Possible ways to treat multi-layer snow fields at times when no snow is present in a
  ! given layer. Note that the public parameters are the only ones that can be used by
  ! calls to hist_addfld2d; the private parameters are just used internally by the
  ! histFile implementation.
  integer , private, parameter :: no_snow_MIN = 1                 ! minimum valid value for this flag
  integer , public , parameter :: no_snow_normal = 1              ! normal treatment, which should be used for most fields (use spval when snow layer not present)
  integer , public , parameter :: no_snow_zero = 2                ! average in a 0 value for times when the snow layer isn't present
  integer , private, parameter :: no_snow_MAX = 2                 ! maximum valid value for this flag
  integer , private, parameter :: no_snow_unset = no_snow_MIN - 1 ! flag specifying that field is NOT a multi-layer snow field
  !
  ! Counters
  !
  integer , public :: ntapes = 0         ! index of max history file requested
  !
  ! Namelist
  !
  integer :: ni                          ! implicit index below
  logical, public :: &
       hist_empty_htapes  = .false.      ! namelist: flag indicates no default history fields
  integer, public :: &
       hist_ndens(max_tapes) = 2         ! namelist: output density of netcdf history files
  integer, public :: &
       hist_mfilt(max_tapes) = (/ 1, (30, ni=2, max_tapes)/)        ! namelist: number of time samples per tape
  logical, public :: &
       hist_dov2xy(max_tapes) = (/.true.,(.true.,ni=2,max_tapes)/) ! namelist: true=> do grid averaging
  integer, public :: &
       hist_nhtfrq(max_tapes) = (/0, (-24, ni=2,max_tapes)/)        ! namelist: history write freq(0=monthly)
  character(len=1), public :: &
       hist_avgflag_pertape(max_tapes) = (/(' ',ni=1,max_tapes)/)   ! namelist: per tape averaging flag
  character(len=max_namlen), public :: &
       hist_type1d_pertape(max_tapes)  = (/(' ',ni=1,max_tapes)/)   ! namelist: per tape type1d

  character(len=max_namlen+2), public :: &
       fincl(max_flds,max_tapes)         ! namelist-equivalence list of fields to add

  character(len=max_namlen+2), public :: &
       hist_fincl1(max_flds) = ' '       ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl2(max_flds) = ' '       ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl3(max_flds) = ' '       ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl4(max_flds) = ' '       ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl5(max_flds) = ' '       ! namelist: list of fields to add
  character(len=max_namlen+2), public :: &
       hist_fincl6(max_flds) = ' '       ! namelist: list of fields to add

  character(len=max_namlen+2), public :: &
       fexcl(max_flds,max_tapes)         ! namelist-equivalence list of fields to remove

  character(len=max_namlen+2), public :: &
       hist_fexcl1(max_flds) = ' ' ! namelist: list of fields to remove
  character(len=max_namlen+2), public :: &
       hist_fexcl2(max_flds) = ' ' ! namelist: list of fields to remove
  character(len=max_namlen+2), public :: &
       hist_fexcl3(max_flds) = ' ' ! namelist: list of fields to remove
  character(len=max_namlen+2), public :: &
       hist_fexcl4(max_flds) = ' ' ! namelist: list of fields to remove
  character(len=max_namlen+2), public :: &
       hist_fexcl5(max_flds) = ' ' ! namelist: list of fields to remove
  character(len=max_namlen+2), public :: &
       hist_fexcl6(max_flds) = ' ' ! namelist: list of fields to remove

  logical, private :: if_disphist(max_tapes)   ! restart, true => save history file
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: hist_addfld1d        ! Add a 1d single-level field to the master field list
  public :: hist_addfld2d        ! Add a 2d multi-level field to the master field list
  public :: hist_addfld_decomp   ! Add a 2d multi-level field to the master field list
  public :: hist_add_subscript   ! Add a 2d subscript dimension
  public :: hist_printflds       ! Print summary of master field list
  public :: hist_htapes_build    ! Initialize history file handler for initial or continue run
  public :: hist_update_hbuf     ! Updates history buffer for all fields and tapes
  public :: hist_htapes_wrapup   ! Write history tape(s)
  public :: hist_restart_ncd     ! Read/write history file restart data
  public :: htapes_fieldlist     ! Define the contents of each history file based on namelist
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: masterlist_make_active    ! Add a field to a history file default "on" list
  private :: masterlist_addfld         ! Add a field to the master field list
  private :: masterlist_change_timeavg ! Override default history tape contents for specific tape
  private :: htape_addfld              ! Add a field to the active list for a history tape
  private :: htape_create              ! Define contents of history file t
  private :: htape_add_ltype_metadata  ! Add global metadata defining landunit types
  private :: htape_add_natpft_metadata ! Add global metadata defining natpft types
  private :: htape_add_cft_metadata    ! Add global metadata defining cft types
  private :: htape_timeconst           ! Write time constant values to history tape
  private :: htape_timeconst3D         ! Write time constant 3D values to primary history tape
  private :: hfields_normalize         ! Normalize history file fields by number of accumulations
  private :: hfields_zero              ! Zero out accumulation and hsitory buffers for a tape
  private :: hfields_write             ! Write a variable to a history tape
  private :: hfields_1dinfo            ! Define/output 1d subgrid info if appropriate
  private :: hist_update_hbuf_field_1d ! Updates history buffer for specific field and tape
  private :: hist_update_hbuf_field_2d ! Updates history buffer for specific field and tape 
  private :: hist_set_snow_field_2d    ! Set values in history field dimensioned by levsno
  private :: list_index                ! Find index of field in exclude list
  private :: set_hist_filename         ! Determine history dataset filenames
  private :: getname                   ! Retrieve name portion of input "inname"
  private :: getflag                   ! Retrieve flag
  private :: pointer_index             ! Track data pointer indices
  private :: max_nFields               ! The max number of fields on any tape
  !
  ! !PRIVATE TYPES:
  ! Constants
  !
  integer, parameter :: max_chars = 256        ! max chars for char variables
  integer, parameter :: max_subs = 100         ! max number of subscripts
  integer            :: num_subs = 0           ! actual number of subscripts
  character(len=32)  :: subs_name(max_subs)    ! name of subscript
  integer            :: subs_dim(max_subs)     ! dimension of subscript
  !
  type field_info
     character(len=max_namlen) :: name         ! field name
     character(len=max_chars)  :: long_name    ! long name
     character(len=max_chars)  :: units        ! units
     character(len=hist_dim_name_length) :: type1d                ! pointer to first dimension type from data type (nameg, etc)
     character(len=hist_dim_name_length) :: type1d_out            ! hbuf first dimension type from data type (nameg, etc)
     character(len=hist_dim_name_length) :: type2d                ! hbuf second dimension type ["levgrnd","levlak","numrad","ltype","natpft","cft","glc_nec","elevclas","subname(n)","month"]
     integer :: beg1d                          ! on-node 1d clm pointer start index
     integer :: end1d                          ! on-node 1d clm pointer end index
     integer :: num1d                          ! size of clm pointer first dimension (all nodes)
     integer :: numdims                        ! the actual number of dimensions, this allows
                                               ! for 2D arrays, where the second dimension is allowed
                                               ! to be 1
     integer :: beg1d_out                      ! on-node 1d hbuf pointer start index
     integer :: end1d_out                      ! on-node 1d hbuf pointer end index
     integer :: num1d_out                      ! size of hbuf first dimension (all nodes)
     integer :: num2d                          ! size of hbuf second dimension (e.g. number of vertical levels)
     integer :: hpindex                        ! history pointer index 
     character(len=8) :: p2c_scale_type        ! scale factor when averaging pft to column
     character(len=8) :: c2l_scale_type        ! scale factor when averaging column to landunit
     character(len=8) :: l2g_scale_type        ! scale factor when averaging landunit to gridcell
     character(len=8) :: t2g_scale_type        ! scale factor when averaging topounit to gridcell
     integer :: no_snow_behavior               ! for multi-layer snow fields, flag saying how to treat times when a given snow layer is absent
  end type field_info

  type master_entry
     type (field_info)  :: field               ! field information
     logical            :: actflag(max_tapes)  ! active/inactive flag
     character(len=1)   :: avgflag(max_tapes)  ! time averaging flag ("X","A","M" or "I",)
  end type master_entry

  type history_entry
     type (field_info) :: field                ! field information
     character(len=1)  :: avgflag              ! time averaging flag
     real(r8), pointer :: hbuf(:,:)            ! history buffer (dimensions: dim1d x num2d)
     integer , pointer :: nacs(:,:)            ! accumulation counter (dimensions: dim1d x num2d)
  end type history_entry

  type history_tape
     integer  :: nflds                         ! number of active fields on tape
     integer  :: ntimes                        ! current number of time samples on tape
     integer  :: mfilt                         ! maximum number of time samples per tape
     integer  :: nhtfrq                        ! number of time samples per tape
     integer  :: ncprec                        ! netcdf output precision
     logical  :: dov2xy                        ! true => do xy average for all fields
     logical  :: is_endhist                    ! true => current time step is end of history interval
     real(r8) :: begtime                       ! time at beginning of history averaging interval
     type (history_entry) :: hlist(max_flds)   ! array of active history tape entries
  end type history_tape

  type elmpoint_rs                             ! Pointer to real scalar data (1D)
     real(r8), pointer :: ptr(:)
  end type elmpoint_rs
  type elmpoint_ra                             ! Pointer to real array data (2D)
     real(r8), pointer :: ptr(:,:)
  end type elmpoint_ra

  ! Pointers into datatype  arrays
  integer, parameter :: max_mapflds = 2500     ! Maximum number of fields to track
  type (elmpoint_rs) :: elmptr_rs(max_mapflds) ! Real scalar data (1D)
  type (elmpoint_ra) :: elmptr_ra(max_mapflds) ! Real array data (2D)
  !
  ! Master list: an array of master_entry entities
  !
  type (master_entry) :: masterlist(max_flds)  ! master field list
  !
  ! History tape: an array of history_tape entities (only active fields)
  !
  type (history_tape) :: tape(max_tapes)       ! array history tapes
  !
  ! Namelist input
  !
  ! Counters
  !
  integer :: nfmaster = 0                        ! number of fields in master field list
  !
  ! Other variables
  !
  character(len=max_chars) :: locfnh(max_tapes)  ! local history file names
  character(len=max_chars) :: locfnhr(max_tapes) ! local history restart file names
  logical :: htapes_defined = .false.            ! flag indicates history contents have been defined
  !
  ! NetCDF  Id's
  !
  type(file_desc_t) :: nfid(max_tapes)       ! file ids
  type(file_desc_t) :: ncid_hist(max_tapes)  ! file ids for history restart files
  integer :: time_dimid                      ! time dimension id
  integer :: hist_interval_dimid             ! time bounds dimension id
  integer :: strlen_dimid                    ! string dimension id
  !
  ! Time Constant variable names and filename
  !
  character(len=max_chars) :: TimeConst3DVars_Filename = ' '
  character(len=max_chars) :: TimeConst3DVars          = ' '
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine hist_printflds()
    !
    ! !DESCRIPTION:
    ! Print summary of master field list.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer nf
    character(len=*),parameter :: subname = 'ELM_hist_printflds'
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) trim(subname),' : number of master fields = ',nfmaster
       write(iulog,*)' ******* MASTER FIELD LIST *******'
       do nf = 1,nfmaster
          write(iulog,9000)nf, masterlist(nf)%field%name, masterlist(nf)%field%units
9000      format (i5,1x,a32,1x,a16)
       end do
       call shr_sys_flush(iulog)
    end if

  end subroutine hist_printflds

  !-----------------------------------------------------------------------
  subroutine masterlist_addfld (fname, type1d, type1d_out, &
        type2d, numdims, num2d, units, avgflag, long_name, hpindex, &
        p2c_scale_type, c2l_scale_type, l2g_scale_type, t2g_scale_type, &
        no_snow_behavior)
    !
    ! !DESCRIPTION:
    ! Add a field to the master field list. Put input arguments of
    ! field name, units, number of levels, averaging flag, and long name
    ! into a type entry in the global master field list (masterlist).
    !
    ! The optional argument no_snow_behavior should be given when this is a multi-layer
    ! snow field, and should be absent otherwise. It should take on one of the no_snow_*
    ! parameters defined above
    !
    ! !ARGUMENTS:
    character(len=*), intent(in)  :: fname            ! field name
    character(len=*), intent(in)  :: type1d           ! 1d data type
    character(len=*), intent(in)  :: type1d_out       ! 1d output type
    character(len=*), intent(in)  :: type2d           ! 2d output type
    integer         , intent(in)  :: numdims          ! number of dimensions
    integer         , intent(in)  :: num2d            ! size of second dimension (e.g. number of vertical levels)
    character(len=*), intent(in)  :: units            ! units of field
    character(len=1), intent(in)  :: avgflag          ! time averaging flag
    character(len=*), intent(in)  :: long_name        ! long name of field
    integer         , intent(in)  :: hpindex          ! data type index for history buffer output
    character(len=*), intent(in)  :: p2c_scale_type   ! scale type for subgrid averaging of pfts to column
    character(len=*), intent(in)  :: c2l_scale_type   ! scale type for subgrid averaging of columns to landunits
    character(len=*), intent(in)  :: l2g_scale_type   ! scale type for subgrid averaging of landunits to gridcells
    character(len=*), intent(in)  :: t2g_scale_type   ! scale type for subgrid averaging of topounits to gridcells
    integer, intent(in), optional :: no_snow_behavior ! if a multi-layer snow field, behavior to use for absent snow layers
    !
    ! !LOCAL VARIABLES:
    integer :: n            ! loop index
    integer :: f            ! masterlist index
    integer :: numa         ! total number of atm cells across all processors
    integer :: numg         ! total number of gridcells across all processors
    integer :: numt         ! total number of topounits across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
    type(bounds_type) :: bounds                  
    character(len=*),parameter :: subname = 'masterlist_addfld'
    !------------------------------------------------------------------------

    ! Determine bounds

    call get_proc_bounds(bounds)
    call get_proc_global(ng=numg, nt=numt, nl=numl, nc=numc, np=nump)

    ! Ensure that new field is not all blanks

    if (fname == ' ') then
       write(iulog,*) trim(subname),' ERROR: blank field name not allowed'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Ensure that new field doesn't already exist

    do n = 1,nfmaster
       if (masterlist(n)%field%name == fname) then
          write(iulog,*) trim(subname),' ERROR:', fname, ' already on list'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    end do

    ! Increase number of fields on master field list

    nfmaster = nfmaster + 1
    f = nfmaster

    ! Check number of fields in master list against maximum number for master list

    if (nfmaster > max_flds) then
       write(iulog,*) trim(subname),' ERROR: too many fields for primary history file ', &
            '-- max_flds,nfmaster=', max_flds, nfmaster
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Add field to master list

    masterlist(f)%field%name           = fname
    masterlist(f)%field%long_name      = long_name
    masterlist(f)%field%units          = units
    masterlist(f)%field%type1d         = type1d
    masterlist(f)%field%type1d_out     = type1d_out
    masterlist(f)%field%type2d         = type2d
    masterlist(f)%field%numdims        = numdims
    masterlist(f)%field%num2d          = num2d
    masterlist(f)%field%hpindex        = hpindex
    masterlist(f)%field%p2c_scale_type = p2c_scale_type
    masterlist(f)%field%c2l_scale_type = c2l_scale_type
    masterlist(f)%field%l2g_scale_type = l2g_scale_type
    masterlist(f)%field%t2g_scale_type = t2g_scale_type

    select case (type1d)
    case (grlnd)
       masterlist(f)%field%beg1d = bounds%begg
       masterlist(f)%field%end1d = bounds%endg
       masterlist(f)%field%num1d = numg
    case (nameg)
       masterlist(f)%field%beg1d = bounds%begg
       masterlist(f)%field%end1d = bounds%endg
       masterlist(f)%field%num1d = numg
    case (namet)
       masterlist(f)%field%beg1d = bounds%begt
       masterlist(f)%field%end1d = bounds%endt
       masterlist(f)%field%num1d = numt
    case (namel)
       masterlist(f)%field%beg1d = bounds%begl
       masterlist(f)%field%end1d = bounds%endl
       masterlist(f)%field%num1d = numl
    case (namec)
       masterlist(f)%field%beg1d = bounds%begc
       masterlist(f)%field%end1d = bounds%endc
       masterlist(f)%field%num1d = numc
    case (namep)
       masterlist(f)%field%beg1d = bounds%begp
       masterlist(f)%field%end1d = bounds%endp
       masterlist(f)%field%num1d = nump
    case default
       write(iulog,*) trim(subname),' ERROR: unknown 1d output type= ',type1d
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    if (present(no_snow_behavior)) then
       masterlist(f)%field%no_snow_behavior = no_snow_behavior
    else
       masterlist(f)%field%no_snow_behavior = no_snow_unset
    end if

    ! The following two fields are used only in master field list,
    ! NOT in the runtime active field list
    ! ALL FIELDS IN THE MASTER LIST ARE INITIALIZED WITH THE ACTIVE
    ! FLAG SET TO FALSE

    masterlist(f)%avgflag(:) = avgflag
    masterlist(f)%actflag(:) = .false.

  end subroutine masterlist_addfld

  !-----------------------------------------------------------------------
  subroutine hist_htapes_build ()
    !
    ! !DESCRIPTION:
    ! Initialize history file for initial or continuation run.  For example,
    ! on an initial run, this routine initializes ``ntapes'' history files.
    ! On a restart run, this routine only initializes history files declared
    ! beyond what existed on the previous run.  Files which already existed on
    ! the previous run have already been initialized (i.e. named and opened)
    ! in routine restart\_history.  Loop over tapes and fields per tape setting
    ! appropriate variables and calling appropriate routines
    !
    ! !USES:
    use clm_time_manager, only: get_prev_time
    use elm_varcon      , only: secspday
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: i                   ! index
    integer :: ier                 ! error code
    integer :: t, f                ! tape, field indices
    integer :: day, sec            ! day and seconds from base date
    character(len=*),parameter :: subname = 'hist_htapes_build'
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*)  trim(subname),' Initializing elm history files'
       write(iulog,'(72a1)') ("-",i=1,60)
       call shr_sys_flush(iulog)
    endif

    ! Define field list information for all history files.
    ! Update ntapes to reflect number of active history files
    ! Note - branch runs can have additional auxiliary history files
    ! declared).

    call htapes_fieldlist()

    ! Determine if gridcell (xy) averaging is done for all fields on tape

    do t=1,ntapes
       tape(t)%dov2xy = hist_dov2xy(t)
       if (masterproc) then
          write(iulog,*)trim(subname),' hist tape = ',t,&
               ' written with dov2xy= ',tape(t)%dov2xy
       end if
    end do

    ! Set number of time samples in each history file and
    ! Note - the following entries will be overwritten by history restart
    ! Note - with netcdf, only 1 (ncd_double) and 2 (ncd_float) are allowed

    do t=1,ntapes
       tape(t)%ntimes = 0
       tape(t)%dov2xy = hist_dov2xy(t)
       tape(t)%nhtfrq = hist_nhtfrq(t)
       tape(t)%mfilt = hist_mfilt(t)
       if (hist_ndens(t) == 1) then
          tape(t)%ncprec = ncd_double
       else
          tape(t)%ncprec = ncd_float
       endif
    end do

    ! Set time of beginning of current averaging interval
    ! First etermine elapased time since reference date

    call get_prev_time(day, sec)
    do t=1,ntapes
       tape(t)%begtime = day + sec/secspday
    end do

    if (masterproc) then
       write(iulog,*)  trim(subname),' Successfully initialized elm history files'
       write(iulog,'(72a1)') ("-",i=1,60)
       call shr_sys_flush(iulog)
    endif

  end subroutine hist_htapes_build

  !-----------------------------------------------------------------------
  subroutine masterlist_make_active (name, tape_index, avgflag)
    !
    ! !DESCRIPTION:
    ! Add a field to the default ``on'' list for a given history file.
    ! Also change the default time averaging flag if requested.
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: name          ! field name
    integer, intent(in) :: tape_index             ! history tape index
    character(len=1), intent(in), optional :: avgflag  ! time averaging flag
    !
    ! !LOCAL VARIABLES:
    integer :: f            ! field index
    logical :: found        ! flag indicates field found in masterlist
    character(len=*),parameter :: subname = 'masterlist_make_active'
    !-----------------------------------------------------------------------

    ! Check validity of input arguments

    if (tape_index > max_tapes) then
       write(iulog,*) trim(subname),' ERROR: tape index=', tape_index, ' is too big'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (present(avgflag)) then
       if ( avgflag /= ' ' .and. &
            avgflag /= 'A' .and. avgflag /= 'I' .and. &
            avgflag /= 'X' .and. avgflag /= 'M') then
          write(iulog,*) trim(subname),' ERROR: unknown averaging flag=', avgflag
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
    end if

    ! Look through master list for input field name.
    ! When found, set active flag for that tape to true.
    ! Also reset averaging flag if told to use other than default.

    found = .false.
    do f = 1,nfmaster
       if (trim(name) == trim(masterlist(f)%field%name)) then
          masterlist(f)%actflag(tape_index) = .true.
          if (present(avgflag)) then
             if (avgflag/= ' ') masterlist(f)%avgflag(tape_index) = avgflag
          end if
          found = .true.
          exit
       end if
    end do
    if (.not. found) then
       write(iulog,*) trim(subname),' ERROR: field=', name, ' not found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

  end subroutine masterlist_make_active

  !-----------------------------------------------------------------------
  subroutine masterlist_change_timeavg (t)
    !
    ! !DESCRIPTION:
    ! Override default history tape contents for a specific tape.
    ! Copy the flag into the master field list.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t         ! history tape index
    !
    ! !LOCAL VARIABLES:
    integer :: f                     ! field index
    character(len=1) :: avgflag      ! lcl equiv of hist_avgflag_pertape(t)
    character(len=*),parameter :: subname = 'masterlist_change_timeavg'
    !-----------------------------------------------------------------------

    avgflag = hist_avgflag_pertape(t)

    do f = 1,nfmaster
       select case (avgflag)
       case ('A')
          masterlist(f)%avgflag(t) = avgflag
       case ('I')
          masterlist(f)%avgflag(t) = avgflag
       case ('X')
          masterlist(f)%avgflag(t) = avgflag
       case ('M')
          masterlist(f)%avgflag(t) = avgflag
       case default
          write(iulog,*) trim(subname),' ERROR: unknown avgflag=',avgflag
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    end do

  end subroutine masterlist_change_timeavg

  !-----------------------------------------------------------------------
  subroutine htapes_fieldlist()
    !
    ! !DESCRIPTION:
    ! Define the contents of each history file based on namelist
    ! input for initial or branch run, and restart data if a restart run.
    ! Use arrays fincl and fexcl to modify default history tape contents.
    ! Then sort the result alphanumerically.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: t, f                         ! tape, field indices
    integer :: ff                           ! index into include, exclude and fprec list
    character(len=max_namlen) :: name       ! field name portion of fincl (i.e. no avgflag separator)
    character(len=max_namlen) :: mastername ! name from masterlist field
    character(len=1)  :: avgflag            ! averaging flag
    character(len=1)  :: prec_acc           ! history buffer precision flag
    character(len=1)  :: prec_wrt           ! history buffer write precision flag
    type (history_entry) :: tmp             ! temporary used for swapping
    character(len=*),parameter :: subname = 'htapes_fieldlist'
    !-----------------------------------------------------------------------

    ! Override averaging flag for all fields on a particular tape
    ! if namelist input so specifies

    do t=1,max_tapes
       if (hist_avgflag_pertape(t) /= ' ') then
          call masterlist_change_timeavg (t)
       end if
    end do

    fincl(:,1) = hist_fincl1(:)
    fincl(:,2) = hist_fincl2(:)
    fincl(:,3) = hist_fincl3(:)
    fincl(:,4) = hist_fincl4(:)
    fincl(:,5) = hist_fincl5(:)
    fincl(:,6) = hist_fincl6(:)

    fexcl(:,1) = hist_fexcl1(:)
    fexcl(:,2) = hist_fexcl2(:)
    fexcl(:,3) = hist_fexcl3(:)
    fexcl(:,4) = hist_fexcl4(:)
    fexcl(:,5) = hist_fexcl5(:)
    fexcl(:,6) = hist_fexcl6(:)


    ! First ensure contents of fincl and fexcl are valid names

    do t = 1,max_tapes
       f = 1
       do while (f < max_flds .and. fincl(f,t) /= ' ')
          name = getname (fincl(f,t))
          do ff = 1,nfmaster
             mastername = masterlist(ff)%field%name
             if (name == mastername) exit
          end do
          if (name /= mastername) then
             write(iulog,*) trim(subname),' ERROR: ', trim(name), ' in fincl(', f, ') ',&
                  'for history tape ',t,' not found'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if
          f = f + 1
       end do

       f = 1
       do while (f < max_flds .and. fexcl(f,t) /= ' ')
          do ff = 1,nfmaster
             mastername = masterlist(ff)%field%name
             if (fexcl(f,t) == mastername) exit
          end do
          if (fexcl(f,t) /= mastername) then
             write(iulog,*) trim(subname),' ERROR: ', fexcl(f,t), ' in fexcl(', f, ') ', &
                  'for history tape ',t,' not found'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if
          f = f + 1
       end do
    end do

    tape(:)%nflds = 0
    do t = 1,max_tapes

       ! Loop through the masterlist set of field names and determine if any of those
       ! are in the FINCL or FEXCL arrays
       ! The call to list_index determines the index in the FINCL or FEXCL arrays
       ! that the masterlist field corresponds to
       ! Add the field to the tape if specified via namelist (FINCL[1-max_tapes]),
       ! or if it is on by default and was not excluded via namelist (FEXCL[1-max_tapes]).

       do f = 1,nfmaster
          mastername = masterlist(f)%field%name
          call list_index (fincl(1,t), mastername, ff)

          if (ff > 0) then

             ! if field is in include list, ff > 0 and htape_addfld
             ! will not be called for field

             avgflag = getflag (fincl(ff,t))
             call htape_addfld (t, f, avgflag)

          else if (.not. hist_empty_htapes) then

             ! find index of field in exclude list

             call list_index (fexcl(1,t), mastername, ff)

             ! if field is in exclude list, ff > 0 and htape_addfld
             ! will not be called for field
             ! if field is not in exclude list, ff =0 and htape_addfld
             ! will be called for field (note that htape_addfld will be
             ! called below only if field is not in exclude list OR in
             ! include list

             if (ff == 0 .and. masterlist(f)%actflag(t)) then
                call htape_addfld (t, f, ' ')
             end if

          end if
       end do

       ! Specification of tape contents now complete.
       ! Sort each list of active entries

       do f = tape(t)%nflds-1,1,-1
          do ff = 1,f
             if (tape(t)%hlist(ff)%field%name > tape(t)%hlist(ff+1)%field%name) then

                tmp = tape(t)%hlist(ff)
                tape(t)%hlist(ff  ) = tape(t)%hlist(ff+1)
                tape(t)%hlist(ff+1) = tmp

             else if (tape(t)%hlist(ff)%field%name == tape(t)%hlist(ff+1)%field%name) then

                write(iulog,*) trim(subname),' ERROR: Duplicate field ', &
                   tape(t)%hlist(ff)%field%name, &
                   't,ff,name=',t,ff,tape(t)%hlist(ff+1)%field%name
                call endrun(msg=errMsg(__FILE__, __LINE__))

             end if
          end do
       end do

       if (masterproc) then
          if (tape(t)%nflds > 0) then
             write(iulog,*) trim(subname),' : Included fields tape ',t,'=',tape(t)%nflds
          end if
          do f = 1,tape(t)%nflds
             write(iulog,*) f,' ',tape(t)%hlist(f)%field%name, &
                  tape(t)%hlist(f)%field%num2d,' ',tape(t)%hlist(f)%avgflag
          end do
          call shr_sys_flush(iulog)
       end if
    end do

    ! Determine total number of active history tapes

    ntapes = 0
    do t = max_tapes,1,-1
       if (tape(t)%nflds > 0) then
          ntapes = t
          exit
       end if
    end do

    ! Ensure there are no "holes" in tape specification, i.e. empty tapes.
    ! Enabling holes should not be difficult if necessary.

    do t = 1,ntapes
       if (tape(t)%nflds  ==  0) then
          write(iulog,*) trim(subname),' ERROR: Tape ',t,' is empty'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    end do

    ! Check that the number of history files declared does not exceed
    ! the maximum allowed.

    if (ntapes > max_tapes) then
       write(iulog,*) trim(subname),' ERROR: Too many history files declared, max_tapes=',max_tapes
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Change 1d output per tape output flag if requested - only for history
    ! tapes where 2d xy averaging is not enabled

    do t = 1,ntapes
       if (hist_type1d_pertape(t) /= ' ' .and. (.not. hist_dov2xy(t))) then
          select case (trim(hist_type1d_pertape(t)))
          case ('PFTS','COLS', 'LAND', 'TOPO', 'GRID')
             if ( masterproc ) &
             write(iulog,*)'history tape ',t,' will have 1d output type of ',hist_type1d_pertape(t)
          case default
             write(iulog,*) trim(subname),' ERROR: unknown namelist type1d per tape=',hist_type1d_pertape(t)
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select
       end if
    end do

    if (masterproc) then
       write(iulog,*) 'There will be a total of ',ntapes,' history tapes'
       do t=1,ntapes
          write(iulog,*)
          if (hist_nhtfrq(t) == 0) then
             write(iulog,*)'History tape ',t,' write frequency is MONTHLY'
          else
             write(iulog,*)'History tape ',t,' write frequency = ',hist_nhtfrq(t)
          endif
          if (hist_dov2xy(t)) then
             write(iulog,*)'All fields on history tape ',t,' are grid averaged'
          else
             write(iulog,*)'All fields on history tape ',t,' are not grid averaged'
          end if
          write(iulog,*)'Number of time samples on history tape ',t,' is ',hist_mfilt(t)
          write(iulog,*)'Output precision on history tape ',t,'=',hist_ndens(t)
          write(iulog,*)
       end do
       call shr_sys_flush(iulog)
    end if

    ! Set flag indicating h-tape contents are now defined (needed by masterlist_addfld)

    htapes_defined = .true.


  end subroutine htapes_fieldlist

  !-----------------------------------------------------------------------
  subroutine htape_addfld (t, f, avgflag)
    !
    ! !DESCRIPTION:
    ! Add a field to the active list for a history tape. Copy the data from
    ! the master field list to the active list for the tape.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t                 ! history tape index
    integer, intent(in) :: f                 ! field index from master field list
    character(len=1), intent(in) :: avgflag  ! time averaging flag
    !
    ! !LOCAL VARIABLES:
    integer :: n                    ! field index on defined tape
    character(len=hist_dim_name_length) :: type1d      ! clm pointer 1d type
    character(len=hist_dim_name_length) :: type1d_out  ! history buffer 1d type
    integer :: numa                 ! total number of atm cells across all processors
    integer :: numg                 ! total number of gridcells across all processors
    integer :: numt                 ! total number of topounits across all processors
    integer :: numl                 ! total number of landunits across all processors
    integer :: numc                 ! total number of columns across all processors
    integer :: nump                 ! total number of pfts across all processors
    integer :: num2d                ! size of second dimension (e.g. .number of vertical levels)
    integer :: beg1d_out,end1d_out  ! history output per-proc 1d beginning and ending indices
    integer :: num1d_out            ! history output 1d size
    type(bounds_type) :: bounds     
    character(len=*),parameter :: subname = 'htape_addfld'
    !-----------------------------------------------------------------------

    ! Ensure that it is not to late to add a field to the history tape

    if (htapes_defined) then
       write(iulog,*) trim(subname),' ERROR: attempt to add field ', &
            masterlist(f)%field%name, ' after history files are set'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    tape(t)%nflds = tape(t)%nflds + 1
    n = tape(t)%nflds

    ! Copy field information

    tape(t)%hlist(n)%field = masterlist(f)%field

    ! Determine bounds

    call get_proc_bounds(bounds)
    call get_proc_global(ng=numg, nt=numt, nl=numl, nc=numc, np=nump)

    ! Modify type1d_out if necessary

    if (hist_dov2xy(t)) then

       ! If xy output averaging is requested, set output 1d type to grlnd
       ! ***NOTE- the following logic is what permits non lat/lon grids to
       ! be written to clm history file

       type1d = tape(t)%hlist(n)%field%type1d

       if (type1d == nameg .or. &
           type1d == namet .or. &
           type1d == namel .or. &
           type1d == namec .or. &
           type1d == namep) then
          tape(t)%hlist(n)%field%type1d_out = grlnd
       end if
       if (type1d == grlnd) then
          tape(t)%hlist(n)%field%type1d_out = grlnd
       end if

    else if (hist_type1d_pertape(t) /= ' ') then

       ! Set output 1d type  based on namelist setting of  hist_type1d_pertape
       ! Only applies to tapes when xy output is not required

       type1d = tape(t)%hlist(n)%field%type1d

       select case (trim(hist_type1d_pertape(t)))
       case('GRID')
          tape(t)%hlist(n)%field%type1d_out = nameg
       case('TOPO')
          tape(t)%hlist(n)%field%type1d_out = namet
       case('LAND')
          tape(t)%hlist(n)%field%type1d_out = namel
       case('COLS')
          tape(t)%hlist(n)%field%type1d_out = namec
       case ('PFTS')
          tape(t)%hlist(n)%field%type1d_out = namep
       case default
          write(iulog,*) trim(subname),' ERROR: unknown input hist_type1d_pertape= ', hist_type1d_pertape(t)
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select

    endif

    ! Determine output 1d dimensions

    type1d_out = tape(t)%hlist(n)%field%type1d_out
    if (type1d_out == grlnd) then
       beg1d_out = bounds%begg
       end1d_out = bounds%endg
       num1d_out = numg
    else if (type1d_out == nameg) then
       beg1d_out = bounds%begg
       end1d_out = bounds%endg
       num1d_out = numg
    else if (type1d_out == namet) then
       beg1d_out = bounds%begt
       end1d_out = bounds%endt
       num1d_out = numt
    else if (type1d_out == namel) then
       beg1d_out = bounds%begl
       end1d_out = bounds%endl
       num1d_out = numl
    else if (type1d_out == namec) then
       beg1d_out = bounds%begc
       end1d_out = bounds%endc
       num1d_out = numc
    else if (type1d_out == namep) then
       beg1d_out = bounds%begp
       end1d_out = bounds%endp
       num1d_out = nump
    else
       write(iulog,*) trim(subname),' ERROR: incorrect value of type1d_out= ',type1d_out
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    tape(t)%hlist(n)%field%beg1d_out = beg1d_out
    tape(t)%hlist(n)%field%end1d_out = end1d_out
    tape(t)%hlist(n)%field%num1d_out = num1d_out
    
    ! Allocate and initialize history buffer and related info

    num2d = tape(t)%hlist(n)%field%num2d
    allocate (tape(t)%hlist(n)%hbuf(beg1d_out:end1d_out,num2d))
    allocate (tape(t)%hlist(n)%nacs(beg1d_out:end1d_out,num2d))
    tape(t)%hlist(n)%hbuf(:,:) = 0._r8
    tape(t)%hlist(n)%nacs(:,:) = 0

    ! Set time averaging flag based on masterlist setting or
    ! override the default averaging flag with namelist setting

    select case (avgflag)
    case (' ')
       tape(t)%hlist(n)%avgflag = masterlist(f)%avgflag(t)
    case ('A','I','X','M')
       tape(t)%hlist(n)%avgflag = avgflag
    case default
       write(iulog,*) trim(subname),' ERROR: unknown avgflag=', avgflag
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine htape_addfld

  !-----------------------------------------------------------------------
  subroutine hist_update_hbuf(bounds)
    !
    ! !DESCRIPTION:
    ! Accumulate (or take min, max, etc. as appropriate) input field
    ! into its history buffer for appropriate tapes.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds   
    !
    ! !LOCAL VARIABLES:
    integer :: t                   ! tape index
    integer :: f                   ! field index
    integer :: numdims             ! number of dimensions
    integer :: num2d               ! size of second dimension (e.g. number of vertical levels)
    character(len=*),parameter :: subname = 'hist_update_hbuf'
    !-----------------------------------------------------------------------

    do t = 1,ntapes
!$OMP PARALLEL DO PRIVATE (f, numdims, num2d)
       do f = 1,tape(t)%nflds
          numdims = tape(t)%hlist(f)%field%numdims
          if ( numdims == 1) then   
             call hist_update_hbuf_field_1d (t, f, bounds)
          else
             num2d = tape(t)%hlist(f)%field%num2d
             call hist_update_hbuf_field_2d (t, f, bounds, num2d)
          end if
       end do
!$OMP END PARALLEL DO
    end do

  end subroutine hist_update_hbuf

  !-----------------------------------------------------------------------
  subroutine hist_update_hbuf_field_1d (t, f, bounds)
    !
    ! !DESCRIPTION:
    ! Accumulate (or take min, max, etc. as appropriate) input field
    ! into its history buffer for appropriate tapes.
    !
    ! This canNOT be called from within a threaded region (see comment below regarding the
    ! call to p2g, and the lack of explicit bounds on its arguments; see also bug 1786)
    !
    ! !USES:
    use subgridAveMod   , only : p2g, c2g, l2g, t2g
    use landunit_varcon , only : istice_mec
    use decompMod       , only : BOUNDS_LEVEL_PROC
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t            ! tape index
    integer, intent(in) :: f            ! field index
    type(bounds_type), intent(in) :: bounds         
    !
    ! !LOCAL VARIABLES:
    integer  :: hpindex                 ! history pointer index
    integer  :: k                       ! gridcell, landunit, column or pft index
    integer  :: beg1d,end1d             ! beginning and ending indices
    logical  :: check_active            ! true => check 'active' flag of each point (this refers to a point being active, NOT a history field being active)
    logical  :: valid                   ! true => history operation is valid
    logical  :: map2gcell               ! true => map clm pointer field to gridcell
    character(len=hist_dim_name_length)  :: type1d         ! 1d clm pointerr type   ["gridcell","landunit","column","pft"]
    character(len=hist_dim_name_length)  :: type1d_out     ! 1d history buffer type ["gridcell","landunit","column","pft"]
    character(len=1)  :: avgflag        ! time averaging flag
    character(len=8)  :: p2c_scale_type ! scale type for subgrid averaging of pfts to column
    character(len=8)  :: c2l_scale_type ! scale type for subgrid averaging of columns to landunits
    character(len=8)  :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
    character(len=8)  :: t2g_scale_type ! scale type for subgrid averaging of topounits to gridcells
    real(r8), pointer :: hbuf(:,:)      ! history buffer
    integer , pointer :: nacs(:,:)      ! accumulation counter
    real(r8), pointer :: field(:)       ! clm 1d pointer field
    logical , pointer :: active(:)      ! flag saying whether each point is active (used for type1d = landunit/column/pft) (this refers to a point being active, NOT a history field being active)
    real(r8) :: field_gcell(bounds%begg:bounds%endg)  ! gricell level field (used if mapping to gridcell is done)
    integer j
    character(len=*),parameter :: subname = 'hist_update_hbuf_field_1d'
    integer k_offset                    ! offset for mapping sliced subarray pointers when outputting variables in PFT/col vector form
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, errMsg(__FILE__, __LINE__))

    avgflag        =  tape(t)%hlist(f)%avgflag
    nacs           => tape(t)%hlist(f)%nacs
    hbuf           => tape(t)%hlist(f)%hbuf
    beg1d          =  tape(t)%hlist(f)%field%beg1d
    end1d          =  tape(t)%hlist(f)%field%end1d
    type1d         =  tape(t)%hlist(f)%field%type1d
    type1d_out     =  tape(t)%hlist(f)%field%type1d_out
    p2c_scale_type =  tape(t)%hlist(f)%field%p2c_scale_type
    c2l_scale_type =  tape(t)%hlist(f)%field%c2l_scale_type
    l2g_scale_type =  tape(t)%hlist(f)%field%l2g_scale_type
    t2g_scale_type =  tape(t)%hlist(f)%field%t2g_scale_type
    hpindex        =  tape(t)%hlist(f)%field%hpindex
    field          => elmptr_rs(hpindex)%ptr

    ! set variables to check weights when allocate all pfts

    map2gcell = .false.
    if (type1d_out == nameg .or. type1d_out == grlnd) then
       if (type1d == namep) then
          ! In this and the following calls, we do NOT explicitly subset field using
	  ! bounds (e.g., we do NOT do field(bounds%begp:bounds%endp). This is because,
          ! for some fields, the lower bound has been reset to 1 due to taking a pointer
          ! to an array slice. Thus, this code will NOT work properly if done within a
          ! threaded region! (See also bug 1786)
          call p2g(bounds, &
               field, &
               field_gcell(bounds%begg:bounds%endg), &
               p2c_scale_type, c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namec) then
          call c2g(bounds, &
               field, &
               field_gcell(bounds%begg:bounds%endg), &
               c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namel) then
          call l2g(bounds, &
               field, &
               field_gcell(bounds%begg:bounds%endg), &
               l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namet) then
          call t2g(bounds, &
               field, &
               field_gcell(bounds%begg:bounds%endg), &
               t2g_scale_type)
          map2gcell = .true.
       end if
    end if

    if (map2gcell) then  ! Map to gridcell

       ! note that in this case beg1d = begg and end1d=endg
       select case (avgflag)
       case ('I') ! Instantaneous
          do k = bounds%begg,bounds%endg
             if (field_gcell(k) /= spval) then
                hbuf(k,1) = field_gcell(k)
             else
                hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case ('A') ! Time average
          do k = bounds%begg,bounds%endg
             if (field_gcell(k) /= spval) then
                if (nacs(k,1) == 0) hbuf(k,1) = 0._r8
                hbuf(k,1) = hbuf(k,1) + field_gcell(k)
                nacs(k,1) = nacs(k,1) + 1
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
          end do
       case ('X') ! Maximum over time
          do k = bounds%begg,bounds%endg
             if (field_gcell(k) /= spval) then
                if (nacs(k,1) == 0) hbuf(k,1) = -1.e50_r8
                hbuf(k,1) = max( hbuf(k,1), field_gcell(k) )
             else
                hbuf(k,1) = spval
             endif
             nacs(k,1) = 1
          end do
       case ('M') ! Minimum over time
          do k = bounds%begg,bounds%endg
             if (field_gcell(k) /= spval) then
                if (nacs(k,1) == 0) hbuf(k,1) = +1.e50_r8
                hbuf(k,1) = min( hbuf(k,1), field_gcell(k) )
             else
                hbuf(k,1) = spval
             endif
             nacs(k,1) = 1
          end do
       case default
          write(iulog,*) trim(subname),' ERROR: invalid time averaging flag ', avgflag
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select

    else  ! Do not map to gridcell

       ! For data defined on the pft, col, and landunit we need to check if a point is active
       ! to determine whether that point should be assigned spval
       if (type1d == namep) then
          check_active = .true.
          active => veg_pp%active
       else if (type1d == namec) then
          check_active = .true.
          active => col_pp%active
       else if (type1d == namel) then
          check_active = .true.
          active =>lun_pp%active
       else
          check_active = .false.
       end if

       select case (avgflag)
       case ('I') ! Instantaneous
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k) /= spval) then
                   hbuf(k,1) = field(k)
                else
                   hbuf(k,1) = spval
                end if
             else
                hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case ('A') ! Time average
          ! create mappings for array slice pointers (which go from 1 to size(field) rather than beg1d to end1d)
          if ( end1d .eq. ubound(field,1) ) then
             k_offset = 0
          else
             k_offset = 1 - beg1d 
          endif
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k+k_offset) /= spval) then   ! add k_offset
                   if (nacs(k,1) == 0) hbuf(k,1) = 0._r8
                   hbuf(k,1) = hbuf(k,1) + field(k+k_offset)   ! add k_offset
                   nacs(k,1) = nacs(k,1) + 1
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
          end do
       case ('X') ! Maximum over time
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k) /= spval) then
                   if (nacs(k,1) == 0) hbuf(k,1) = -1.e50_r8
                   hbuf(k,1) = max( hbuf(k,1), field(k) )
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case ('M') ! Minimum over time
          do k = beg1d,end1d
             valid = .true.
             if (check_active) then
                if (.not. active(k)) valid = .false.
             end if
             if (valid) then
                if (field(k) /= spval) then
                   if (nacs(k,1) == 0) hbuf(k,1) = +1.e50_r8
                   hbuf(k,1) = min( hbuf(k,1), field(k) )
                else
                   if (nacs(k,1) == 0) hbuf(k,1) = spval
                end if
             else
                if (nacs(k,1) == 0) hbuf(k,1) = spval
             end if
             nacs(k,1) = 1
          end do
       case default
          write(iulog,*) trim(subname),' ERROR: invalid time averaging flag ', avgflag
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    end if

  end subroutine hist_update_hbuf_field_1d

  !-----------------------------------------------------------------------
  subroutine hist_update_hbuf_field_2d (t, f, bounds, num2d)
    !
    ! !DESCRIPTION:
    ! Accumulate (or take min, max, etc. as appropriate) input field
    ! into its history buffer for appropriate tapes.
    !
    ! This canNOT be called from within a threaded region (see comment below regarding the
    ! call to p2g, and the lack of explicit bounds on its arguments; see also bug 1786)
    !
    ! !USES:
    use subgridAveMod   , only : p2g, c2g, l2g, t2g
    use landunit_varcon , only : istice_mec
    use decompMod       , only : BOUNDS_LEVEL_PROC
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t            ! tape index
    integer, intent(in) :: f            ! field index
    type(bounds_type), intent(in) :: bounds         
    integer, intent(in) :: num2d        ! size of second dimension
    !
    ! !LOCAL VARIABLES:
    integer  :: hpindex                 ! history pointer index
    integer  :: k                       ! gridcell, landunit, column or pft index
    integer  :: j                       ! level index
    integer  :: beg1d,end1d             ! beginning and ending indices
    logical  :: check_active            ! true => check 'active' flag of each point (this refers to a point being active, NOT a history field being active)
    logical  :: valid                   ! true => history operation is valid
    logical  :: map2gcell               ! true => map clm pointer field to gridcell
    character(len=hist_dim_name_length)  :: type1d         ! 1d clm pointerr type   ["gridcell","landunit","column","pft"]
    character(len=hist_dim_name_length)  :: type1d_out     ! 1d history buffer type ["gridcell","landunit","column","pft"]
    character(len=1)  :: avgflag        ! time averaging flag
    character(len=8)  :: p2c_scale_type ! scale type for subgrid averaging of pfts to column
    character(len=8)  :: c2l_scale_type ! scale type for subgrid averaging of columns to landunits
    character(len=8)  :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
    character(len=8)  :: t2g_scale_type ! scale type for subgrid averaging of topounits to gridcells
    integer  :: no_snow_behavior        ! for multi-layer snow fields, behavior to use when a given layer is absent
    real(r8), pointer :: hbuf(:,:)      ! history buffer
    integer , pointer :: nacs(:,:)      ! accumulation counter
    real(r8), pointer :: field(:,:)     ! clm 2d pointer field
    logical           :: field_allocated! whether 'field' was allocated here
    logical , pointer :: active(:)      ! flag saying whether each point is active (used for type1d = landunit/column/pft) 
                                        !(this refers to a point being active, NOT a history field being active)
    real(r8) :: field_gcell(bounds%begg:bounds%endg,num2d) ! gricell level field (used if mapping to gridcell is done)
    character(len=*),parameter :: subname = 'hist_update_hbuf_field_2d'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, errMsg(__FILE__, __LINE__))

    avgflag             =  tape(t)%hlist(f)%avgflag
    nacs                => tape(t)%hlist(f)%nacs
    hbuf                => tape(t)%hlist(f)%hbuf
    beg1d               =  tape(t)%hlist(f)%field%beg1d
    end1d               =  tape(t)%hlist(f)%field%end1d
    type1d              =  tape(t)%hlist(f)%field%type1d
    type1d_out          =  tape(t)%hlist(f)%field%type1d_out
    p2c_scale_type      =  tape(t)%hlist(f)%field%p2c_scale_type
    c2l_scale_type      =  tape(t)%hlist(f)%field%c2l_scale_type
    l2g_scale_type      =  tape(t)%hlist(f)%field%l2g_scale_type
    t2g_scale_type      =  tape(t)%hlist(f)%field%t2g_scale_type
    no_snow_behavior    =  tape(t)%hlist(f)%field%no_snow_behavior
    hpindex             =  tape(t)%hlist(f)%field%hpindex

    if (no_snow_behavior /= no_snow_unset) then
       ! For multi-layer snow fields, build a special output variable that handles
       ! missing snow layers appropriately

       ! Note, regarding bug 1786: The following allocation is not what we would want if
       ! this routine were operating in a threaded region (or, more generally, within a
       ! loop over nclumps) - in that case we would want to use the bounds information for
       ! this clump. But currently that's not possible because the bounds of some fields
       ! have been reset to 1 - see also bug 1786. Similarly, if we wanted to allow
       ! operation within a loop over clumps, we would need to pass 'bounds' to
       ! hist_set_snow_field_2d rather than relying on beg1d & end1d (which give the proc,
       ! bounds not the clump bounds)

       allocate(field(lbound(elmptr_ra(hpindex)%ptr, 1) : ubound(elmptr_ra(hpindex)%ptr, 1), 1:num2d))
       field_allocated = .true.

       call hist_set_snow_field_2d(field, elmptr_ra(hpindex)%ptr, no_snow_behavior, type1d, &
            beg1d, end1d)
    else
       field => elmptr_ra(hpindex)%ptr(:,1:num2d)
       field_allocated = .false.
    end if

    ! set variables to check weights when allocate all pfts

    map2gcell = .false.
    if (type1d_out == nameg .or. type1d_out == grlnd) then
       if (type1d == namep) then
          ! In this and the following calls, we do NOT explicitly subset field using
	  ! (e.g., we do NOT do field(bounds%begp:bounds%endp). This is because,
          ! for some fields, the lower bound has been reset to 1 due to taking a pointer
          ! to an array slice. Thus, this code will NOT work properly if done within a
          ! threaded region! (See also bug 1786)
          call p2g(bounds, num2d, &
               field, &
               field_gcell(bounds%begg:bounds%endg, :), &
               p2c_scale_type, c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namec) then
          call c2g(bounds, num2d, &
               field, &
               field_gcell(bounds%begg:bounds%endg, :), &
               c2l_scale_type, l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namel) then
          call l2g(bounds, num2d, &
               field, &
               field_gcell(bounds%begg:bounds%endg, :), &
               l2g_scale_type)
          map2gcell = .true.
       else if (type1d == namet) then
          call t2g(bounds, num2d, &
               field, &
               field_gcell(bounds%begg:bounds%endg, :), &
               t2g_scale_type)
          map2gcell = .true.
       end if
    end if

    if (map2gcell) then  ! Map to gridcell

       ! note that in this case beg1d = begg and end1d=endg
       select case (avgflag)
       case ('I') ! Instantaneous
          do j = 1,num2d
             do k = bounds%begg,bounds%endg
                if (field_gcell(k,j) /= spval) then
                   hbuf(k,j) = field_gcell(k,j)
                else
                   hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case ('A') ! Time average
          do j = 1,num2d
             do k = bounds%begg,bounds%endg
                if (field_gcell(k,j) /= spval) then
                   if (nacs(k,j) == 0) hbuf(k,j) = 0._r8
                   hbuf(k,j) = hbuf(k,j) + field_gcell(k,j)
                   nacs(k,j) = nacs(k,j) + 1
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                endif
             end do
          end do
       case ('X') ! Maximum over time
          do j = 1,num2d
             do k = bounds%begg,bounds%endg
                if (field_gcell(k,j) /= spval) then
                   if (nacs(k,j) == 0) hbuf(k,j) = -1.e50_r8
                   hbuf(k,j) = max( hbuf(k,j), field_gcell(k,j) )
                else
                   hbuf(k,j) = spval
                endif
                nacs(k,j) = 1
             end do
          end do
       case ('M') ! Minimum over time
          do j = 1,num2d
             do k = bounds%begg,bounds%endg
                if (field_gcell(k,j) /= spval) then
                   if (nacs(k,j) == 0) hbuf(k,j) = +1.e50_r8
                   hbuf(k,j) = min( hbuf(k,j), field_gcell(k,j) )
                else
                   hbuf(k,j) = spval
                endif
                nacs(k,j) = 1
             end do
          end do
       case default
          write(iulog,*) trim(subname),' ERROR: invalid time averaging flag ', avgflag
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select

    else  ! Do not map to gridcell

       ! For data defined on the pft, col or landunit, we need to check if a point is active
       ! to determine whether that point should be assigned spval
       if (type1d == namep) then
          check_active = .true.
          active => veg_pp%active
       else if (type1d == namec) then
          check_active = .true.
          active => col_pp%active
       else if (type1d == namel) then
          check_active = .true.
          active =>lun_pp%active
       else
          check_active = .false.
       end if

       ! Note that since field points to an array section the
       ! bounds are field(1:end1d-beg1d+1, num2d) - therefore
       ! need to do the shifting below

       select case (avgflag)
       case ('I') ! Instantaneous
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      hbuf(k,j) = field(k-beg1d+1,j)
                   else
                      hbuf(k,j) = spval
                   end if
                else
                   hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case ('A') ! Time average
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      if (nacs(k,j) == 0) hbuf(k,j) = 0._r8
                      hbuf(k,j) = hbuf(k,j) + field(k-beg1d+1,j)
                      nacs(k,j) = nacs(k,j) + 1
                   else
                      if (nacs(k,j) == 0) hbuf(k,j) = spval
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
             end do
          end do
       case ('X') ! Maximum over time
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      if (nacs(k,j) == 0) hbuf(k,j) = -1.e50_r8
                      hbuf(k,j) = max( hbuf(k,j), field(k-beg1d+1,j) )
                   else
                      if (nacs(k,j) == 0) hbuf(k,j) = spval
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case ('M') ! Minimum over time
          do j = 1,num2d
             do k = beg1d,end1d
                valid = .true.
                if (check_active) then
                   if (.not. active(k)) valid = .false.
                end if
                if (valid) then
                   if (field(k-beg1d+1,j) /= spval) then
                      if (nacs(k,j) == 0) hbuf(k,j) = +1.e50_r8
                      hbuf(k,j) = min( hbuf(k,j), field(k-beg1d+1,j))
                   else
                      if (nacs(k,j) == 0) hbuf(k,j) = spval
                   end if
                else
                   if (nacs(k,j) == 0) hbuf(k,j) = spval
                end if
                nacs(k,j) = 1
             end do
          end do
       case default
          write(iulog,*) trim(subname),' ERROR: invalid time averaging flag ', avgflag
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    end if

    if (field_allocated) then
       deallocate(field)
    end if

  end subroutine hist_update_hbuf_field_2d

  !-----------------------------------------------------------------------
  subroutine hist_set_snow_field_2d (field_out, field_in, no_snow_behavior, type1d, beg1d, end1d)
    !
    ! !DESCRIPTION:
    ! Set values in history field dimensioned by levsno. 
    !
    ! This routine handles what to do when a given snow layer doesn't exist for a given
    ! point, based on the no_snow_behavior argument. Options are:
    !
    ! - no_snow_normal: This is the normal behavior, which applies to most snow fields:
    !   Use spval (missing value flag). This means that temporal averages will just
    !   consider times when a particular snow layer actually existed
    !
    ! - no_snow_zero: Average in a 0 value for times when the snow layer isn't present
    !
    ! Input and output fields can be defined at the pft or column level
    !
    ! !ARGUMENTS:
    integer         , intent(in)  :: beg1d                    ! beginning spatial index
    integer         , intent(in)  :: end1d                    ! ending spatial index
    real(r8)        , intent(out) :: field_out( beg1d: , 1: ) ! output field [point, lev]
    real(r8)        , intent(in)  :: field_in ( beg1d: , 1: ) ! input field [point, lev]
    integer         , intent(in)  :: no_snow_behavior         ! behavior to use when a snow layer is absent
    character(len=*), intent(in)  :: type1d                   ! 1d clm pointer type ("column" or "pft")
    !
    ! !LOCAL VARIABLES:
    integer :: num_levels             ! total number of possible snow layers
    integer :: point
    integer :: level
    integer :: num_snow_layers        ! number of snow layers that exist at a point
    integer :: num_nonexistent_layers
    integer :: c                      ! column index
    real(r8):: no_snow_val            ! value to use when a snow layer is missing
    character(len=*), parameter :: subname = 'hist_set_snow_field_2d'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(field_out, 1) == end1d), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(field_in , 1) == end1d), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(field_out, 2) == ubound(field_in, 2)), errMsg(__FILE__, __LINE__))

    associate(&
    snl            => col_pp%snl  &   ! Input: [integer (:)] number of snow layers (negative)
    )

    num_levels = ubound(field_in, 2)

    ! Determine no_snow_val
    select case (no_snow_behavior)
    case (no_snow_normal)
       no_snow_val = spval
    case (no_snow_zero)
       no_snow_val = 0._r8
    case default
       write(iulog,*) trim(subname), ' ERROR: unrecognized no_snow_behavior: ', &
            no_snow_behavior
       call endrun()
    end select

    do point = beg1d, end1d

       ! Get number of snow layers at this point

       if (type1d == namec) then
          c = point
       else if (type1d == namep) then
          c = veg_pp%column(point)
       else
          write(iulog,*) trim(subname), ' ERROR: Only implemented for pft and col-level fields'
          write(iulog,*) 'type1d = ', trim(type1d)
          call endrun()
       end if

       num_snow_layers = abs(snl(c))
       num_nonexistent_layers = num_levels - num_snow_layers
       
       ! Fill output field appropriately for each layer
       ! When only a subset of snow layers exist, it is the LAST num_snow_layers that exist

       do level = 1, num_nonexistent_layers
          field_out(point, level) = no_snow_val
       end do
       do level = (num_nonexistent_layers + 1), num_levels
          field_out(point, level) = field_in(point, level)
       end do
          
    end do

    end associate

  end subroutine hist_set_snow_field_2d


  !-----------------------------------------------------------------------
  subroutine hfields_normalize (t)
    !
    ! !DESCRIPTION:
    ! Normalize fields on a history file by the number of accumulations.
    ! Loop over fields on the tape.  Need averaging flag and number of
    ! accumulations to perform normalization.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t       ! tape index
    !
    ! !LOCAL VARIABLES:
    integer :: f                   ! field index
    integer :: k                   ! 1d index
    integer :: j                   ! 2d index
    logical :: aflag               ! averaging flag
    integer :: beg1d_out,end1d_out ! hbuf 1d beginning and ending indices
    integer :: num2d               ! hbuf size of second dimension (e.g. number of vertical levels)
    character(len=1)  :: avgflag   ! averaging flag
    real(r8), pointer :: hbuf(:,:) ! history buffer
    integer , pointer :: nacs(:,:) ! accumulation counter
    character(len=*),parameter :: subname = 'hfields_normalize'
    !-----------------------------------------------------------------------

    ! Normalize by number of accumulations for time averaged case

    do f = 1,tape(t)%nflds
       avgflag   =  tape(t)%hlist(f)%avgflag
       beg1d_out =  tape(t)%hlist(f)%field%beg1d_out
       end1d_out =  tape(t)%hlist(f)%field%end1d_out
       num2d     =  tape(t)%hlist(f)%field%num2d
       nacs      => tape(t)%hlist(f)%nacs
       hbuf      => tape(t)%hlist(f)%hbuf

       if (avgflag == 'A') then
          aflag = .true.
       else
          aflag = .false.
       end if

       do j = 1, num2d
          do k = beg1d_out, end1d_out
             if (aflag .and. nacs(k,j) /= 0) then
                hbuf(k,j) = hbuf(k,j) / float(nacs(k,j))
             end if
          end do
       end do
    end do

  end subroutine hfields_normalize

  !-----------------------------------------------------------------------
  subroutine hfields_zero (t)
    !
    ! !DESCRIPTION:
    ! Zero out accumulation and history buffers for a given history tape.
    ! Loop through fields on the tape.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t     ! tape index
    !
    ! !LOCAL VARIABLES:
    integer :: f                 ! field index
    character(len=*),parameter :: subname = 'hfields_zero'
    !-----------------------------------------------------------------------

    do f = 1,tape(t)%nflds
       tape(t)%hlist(f)%hbuf(:,:) = 0._r8
       tape(t)%hlist(f)%nacs(:,:) = 0
    end do

  end subroutine hfields_zero

  !-----------------------------------------------------------------------
  subroutine htape_create (t, histrest)
    !
    ! !DESCRIPTION:
    ! Define contents of history file t. Issue the required netcdf
    ! wrapper calls to define the history file contents.
    !
    ! !USES:
    use elm_varpar      , only : nlevgrnd, nlevsno, nlevlak, nlevurb, numrad, nmonth
    use elm_varpar      , only : natpft_size, cft_size, maxpatch_glcmec, nlevdecomp_full, nlevtrc_full, nvegwcs
    use landunit_varcon , only : max_lunit
    use elm_varctl      , only : caseid, ctitle, fsurdat, finidat, paramfile
    use elm_varctl      , only : version, hostname, username, conventions, source
    use domainMod       , only : ldomain
    use fileutils       , only : get_filename
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t                   ! tape index
    logical, intent(in), optional :: histrest  ! if creating the history restart file
    !
    ! !LOCAL VARIABLES:
    integer :: f                   ! field index
    integer :: p,c,l,n             ! indices
    integer :: ier                 ! error code
    integer :: num2d               ! size of second dimension (e.g. number of vertical levels)
    integer :: dimid               ! dimension id temporary
    integer :: dim1id(1)           ! netCDF dimension id
    integer :: dim2id(2)           ! netCDF dimension id
    integer :: ndims               ! dimension counter
    integer :: omode               ! returned mode from netCDF call
    integer :: ncprec              ! output netCDF write precision
    integer :: ret                 ! netCDF error status
    integer :: nump                ! total number of pfts across all processors
    integer :: numc                ! total number of columns across all processors
    integer :: numl                ! total number of landunits across all processors
    integer :: numt                ! total number of topounits across all processors
    integer :: numg                ! total number of gridcells across all processors
    integer :: numa                ! total number of atm cells across all processors
    logical :: avoid_pnetcdf       ! whether we should avoid using pnetcdf
    logical :: lhistrest           ! local history restart flag
    type(file_desc_t) :: lnfid     ! local file id
    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time
    character(len=256) :: name     ! name of attribute
    character(len=256) :: units    ! units of attribute
    character(len=256) :: str      ! global attribute string
    character(len=  1) :: avgflag  ! time averaging flag
    character(len=*),parameter :: subname = 'htape_create'
    !-----------------------------------------------------------------------

    if ( present(histrest) )then
       lhistrest = histrest
    else
       lhistrest = .false.
    end if

    ! Determine necessary indices

    call get_proc_global(ng=numg, nt=numt, nl=numl, nc=numc, np=nump)

    ! define output write precsion for tape

    ncprec = tape(t)%ncprec
    
    ! BUG(wjs, 2014-10-20, bugz 1730) Workaround for
    ! http://bugs.cgd.ucar.edu/show_bug.cgi?id=1730
    ! - 1-d hist files have problems with pnetcdf. A better workaround in terms of
    ! performance is to keep pnetcdf, but set PIO_BUFFER_SIZE_LIMIT=0, but that can't be
    ! done on a per-file basis.
    if (.not. tape(t)%dov2xy) then
       avoid_pnetcdf = .true.
    else
       avoid_pnetcdf = .false.
    end if

    ! Create new netCDF file. It will be in define mode

    if ( .not. lhistrest )then
       if (masterproc) then
          write(iulog,*) trim(subname),' : Opening netcdf htape ', &
                                      trim(locfnh(t))
          call shr_sys_flush(iulog)
       end if
       call ncd_pio_createfile(lnfid, trim(locfnh(t)), avoid_pnetcdf=avoid_pnetcdf)
       call ncd_putatt(lnfid, ncd_global, 'title', 'ELM History file information' )
    else
       if (masterproc) then
          write(iulog,*) trim(subname),' : Opening netcdf rhtape ', &
                                      trim(locfnhr(t))
          call shr_sys_flush(iulog)
       end if
       call ncd_pio_createfile(lnfid, trim(locfnhr(t)), avoid_pnetcdf=avoid_pnetcdf)
       call ncd_putatt(lnfid, ncd_global, 'title', &
          'ELM Restart History information, required to continue a simulation' )
       call ncd_putatt(lnfid, ncd_global, 'comment', &
                       "This entire file NOT needed for startup or branch simulations")
    end if

    ! Create global attributes. Attributes are used to store information
    ! about the data set. Global attributes are information about the
    ! data set as a whole, as opposed to a single variable

    call ncd_putatt(lnfid, ncd_global, 'source'  , trim(source))
    call ncd_putatt(lnfid, ncd_global, 'case', trim(caseid))
    call ncd_putatt(lnfid, ncd_global, 'username', trim(username))
    call ncd_putatt(lnfid, ncd_global, 'hostname', trim(hostname))
    call ncd_putatt(lnfid, ncd_global, 'git_version' , trim(version))
    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call ncd_putatt(lnfid, ncd_global, 'history' , trim(str))
    call ncd_putatt(lnfid, ncd_global, 'institution_id', 'E3SM-Project')
    call ncd_putatt(lnfid, ncd_global, 'contact', &
          'e3sm-data-support@listserv.llnl.gov')
    call ncd_putatt(lnfid, ncd_global, 'Conventions', trim(conventions))
    call ncd_putatt(lnfid, ncd_global, 'comment', &
          "NOTE: None of the variables are weighted by land fraction!" )
    str = get_filename(fsurdat)
    call ncd_putatt(lnfid, ncd_global, 'Surface_dataset', trim(str))
    if (finidat == ' ') then
       str = 'arbitrary initialization'
    else
       str = get_filename(finidat)
    endif
    call ncd_putatt(lnfid, ncd_global, 'Initial_conditions_dataset', trim(str))
    str = get_filename(paramfile)
    call ncd_putatt(lnfid, ncd_global, 'PFT_physiological_constants_dataset', trim(str))

    ! Define dimensions.
    ! Time is an unlimited dimension. Character string is treated as an array of characters.

    ! Global uncompressed dimensions (including non-land points)
    if (ldomain%isgrid2d) then
       call ncd_defdim(lnfid, 'lon'   , ldomain%ni, dimid)
       call ncd_defdim(lnfid, 'lat'   , ldomain%nj, dimid)
    else
       call ncd_defdim(lnfid, trim(grlnd), ldomain%ns, dimid)
    end if

    ! Global compressed dimensions (not including non-land points)
    call ncd_defdim(lnfid, trim(nameg), numg, dimid)
    call ncd_defdim(lnfid, trim(namet), numt, dimid)
    call ncd_defdim(lnfid, trim(namel), numl, dimid)
    call ncd_defdim(lnfid, trim(namec), numc, dimid)
    call ncd_defdim(lnfid, trim(namep), nump, dimid)

    ! "level" dimensions
    call ncd_defdim(lnfid, 'levgrnd', nlevgrnd, dimid)
    if (nlevurb > 0) then
       call ncd_defdim(lnfid, 'levurb' , nlevurb, dimid)
    end if
    call ncd_defdim(lnfid, 'levlak' , nlevlak, dimid)
    call ncd_defdim(lnfid, 'numrad' , numrad , dimid)
    call ncd_defdim(lnfid, 'month'  , nmonth,  dimid)
    call ncd_defdim(lnfid, 'levsno' , nlevsno , dimid)
    call ncd_defdim(lnfid, 'ltype', max_lunit, dimid)
    call ncd_defdim(lnfid, 'nvegwcs',nvegwcs, dimid)
    call htape_add_ltype_metadata(lnfid)
    call ncd_defdim(lnfid, 'natpft', natpft_size, dimid)
    if (cft_size > 0) then
       call ncd_defdim(lnfid, 'cft', cft_size, dimid)
       call htape_add_cft_metadata(lnfid)
    end if
    if (maxpatch_glcmec > 0) then
       call ncd_defdim(lnfid, 'glc_nec' , maxpatch_glcmec , dimid)
       ! elevclas (in contrast to glc_nec) includes elevation class 0 (bare land)
       ! (although on the history file it will go 1:(nec+1) rather than 0:nec)
       call ncd_defdim(lnfid, 'elevclas' , maxpatch_glcmec + 1, dimid)
    end if

    do n = 1,num_subs
       call ncd_defdim(lnfid, subs_name(n), subs_dim(n), dimid)
    end do
    call ncd_defdim(lnfid, 'string_length', hist_dim_name_length, strlen_dimid)
    call ncd_defdim( lnfid, 'levdcmp', nlevdecomp_full, dimid)
    call ncd_defdim( lnfid, 'levtrc', nlevtrc_full, dimid)    
    
    if(use_fates)then
       call ncd_defdim(lnfid, 'fates_levscag', nlevsclass_fates * nlevage_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levscls', nlevsclass_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levcacls',nlevcoage, dimid)
       call ncd_defdim(lnfid, 'fates_levpft', numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levage', nlevage_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levfuel', nfsc_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levcwdsc', ncwd_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levscpf', nlevsclass_fates*numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levcapf',  nlevcoage*numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levcan', nclmax_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levcnlf', nlevleaf_fates * nclmax_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levcnlfpf', nlevleaf_fates * nclmax_fates * numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levscagpf', nlevsclass_fates * nlevage_fates * numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levagepft', nlevage_fates * numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levheight', nlevheight_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levelem', nelements_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levelpft', nelements_fates * numpft_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levelcwd', nelements_fates * ncwd_fates, dimid)
       call ncd_defdim(lnfid, 'fates_levelage', nelements_fates * nlevage_fates, dimid)
    end if

    if ( .not. lhistrest )then
       call ncd_defdim(lnfid, 'hist_interval', 2, hist_interval_dimid)
       call ncd_defdim(lnfid, 'time', ncd_unlimited, time_dimid)
       nfid(t) = lnfid
       if (masterproc)then
          write(iulog,*) trim(subname), &
                          ' : Successfully defined netcdf history file ',t
          call shr_sys_flush(iulog)
       end if
    else
       ncid_hist(t) = lnfid
       if (masterproc)then
          write(iulog,*) trim(subname), &
                          ' : Successfully defined netcdf restart history file ',t
          call shr_sys_flush(iulog)
       end if
    end if

  end subroutine htape_create

  !-----------------------------------------------------------------------
  subroutine htape_add_ltype_metadata(lnfid)
    !
    ! !DESCRIPTION:
    ! Add global metadata defining landunit types
    !
    ! !USES:
    use landunit_varcon, only : max_lunit, landunit_names, landunit_name_length
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: lnfid ! local file id
    !
    ! !LOCAL VARIABLES:
    integer :: ltype  ! landunit type
    character(len=*), parameter :: att_prefix = 'ltype_'  ! prefix for attributes
    character(len=len(att_prefix)+landunit_name_length) :: attname ! attribute name

    character(len=*), parameter :: subname = 'htape_add_ltype_metadata'
    !-----------------------------------------------------------------------
    
    do ltype = 1, max_lunit
       attname = att_prefix // landunit_names(ltype)
       call ncd_putatt(lnfid, ncd_global, attname, ltype)
    end do

  end subroutine htape_add_ltype_metadata

  !-----------------------------------------------------------------------
  subroutine htape_add_natpft_metadata(lnfid)
    !
    ! !DESCRIPTION:
    ! Add global metadata defining natpft types
    !
    ! !USES:
    use elm_varpar, only : natpft_lb, natpft_ub
    use pftvarcon , only : pftname_len, pftname
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: lnfid ! local file id
    !
    ! !LOCAL VARIABLES:
    integer :: ptype  ! pft type
    integer :: ptype_1_indexing ! pft type, translated to 1 indexing
    character(len=*), parameter :: att_prefix = 'natpft_' ! prefix for attributes
    character(len=len(att_prefix)+pftname_len) :: attname ! attribute name

    character(len=*), parameter :: subname = 'htape_add_natpft_metadata'
    !-----------------------------------------------------------------------
    
    do ptype = natpft_lb, natpft_ub
       ptype_1_indexing = ptype + (1 - natpft_lb)
       attname = att_prefix // pftname(ptype)
       call ncd_putatt(lnfid, ncd_global, attname, ptype_1_indexing)
    end do

  end subroutine htape_add_natpft_metadata

  !-----------------------------------------------------------------------
  subroutine htape_add_cft_metadata(lnfid)
    !
    ! !DESCRIPTION:
    ! Add global metadata defining natpft types
    !
    ! !USES:
    use elm_varpar, only : cft_lb, cft_ub
    use pftvarcon , only : pftname_len, pftname
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: lnfid ! local file id
    !
    ! !LOCAL VARIABLES:
    integer :: ptype  ! pft type
    integer :: ptype_1_indexing ! pft type, translated to 1 indexing
    character(len=*), parameter :: att_prefix = 'cft_'    ! prefix for attributes
    character(len=len(att_prefix)+pftname_len) :: attname ! attribute name

    character(len=*), parameter :: subname = 'htape_add_cft_metadata'
    !-----------------------------------------------------------------------
    
    do ptype = cft_lb, cft_ub
       ptype_1_indexing = ptype + (1 - cft_lb)
       attname = att_prefix // pftname(ptype)
       call ncd_putatt(lnfid, ncd_global, attname, ptype_1_indexing)
    end do

  end subroutine htape_add_cft_metadata

  !-----------------------------------------------------------------------
  subroutine htape_timeconst3D(t, &
       bounds, watsat_col, sucsat_col, bsw_col, hksat_col, mode)
    !
    ! !DESCRIPTION:
    ! Write time constant 3D variables to history tapes.
    ! Only write out when this subroutine is called (normally only for
    ! primary history files at very first time-step, nstep=0).
    ! Issue the required netcdf wrapper calls to define the history file
    ! contents.
    !
    ! !USES:
    use subgridAveMod  , only : c2g
    use elm_varpar     , only : nlevgrnd ,nlevlak
    use shr_string_mod , only : shr_string_listAppend
    use domainMod      , only : ldomain
    !
    ! !ARGUMENTS:
    integer           , intent(in) :: t    ! tape index
    type(bounds_type) , intent(in) :: bounds           
    real(r8)          , intent(in) :: watsat_col( bounds%begc:,1: ) 
    real(r8)          , intent(in) :: sucsat_col( bounds%begc:,1: ) 
    real(r8)          , intent(in) :: bsw_col( bounds%begc:,1: ) 
    real(r8)          , intent(in) :: hksat_col( bounds%begc:,1: ) 
    character(len=*)  , intent(in) :: mode ! 'define' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: c,l,lev,ifld               ! indices
    integer :: ier                        ! error status
    character(len=max_chars) :: long_name ! variable long name
    character(len=max_namlen):: varname   ! variable name
    character(len=max_namlen):: units     ! variable units
    character(len=8) :: l2g_scale_type    ! scale type for subgrid averaging of landunits to grid cells
    !
    real(r8), pointer :: histi(:,:)       ! temporary
    real(r8), pointer :: histo(:,:)       ! temporary
    integer, parameter :: nflds = 6       ! Number of 3D time-constant fields
    character(len=*),parameter :: subname = 'htape_timeconst3D'
    character(len=*),parameter :: varnames(nflds) = (/ &
                                                        'ZSOI  ', &
                                                        'DZSOI ', &
                                                        'WATSAT', &
                                                        'SUCSAT', &
                                                        'BSW   ', &
                                                        'HKSAT '  &
                                                    /)
    real(r8), pointer :: histil(:,:)      ! temporary
    real(r8), pointer :: histol(:,:)
    integer, parameter :: nfldsl = 2
    character(len=*),parameter :: varnamesl(nfldsl) = (/ &
                                                          'ZLAKE ', &
                                                          'DZLAKE' &
                                                      /)
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(watsat_col) == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sucsat_col) == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bsw_col)    == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(hksat_col)  == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

    !-------------------------------------------------------------------------------
    !***      Non-time varying 3D fields                    ***
    !***      Only write out when this subroutine is called ***
    !***       Normally only called once for primary tapes  ***
    !-------------------------------------------------------------------------------

    if (mode == 'define') then

       do ifld = 1,nflds
          ! Field indices MUST match varnames array order above!
          if (ifld == 1) then
             long_name='soil depth'; units = 'm'
          else if (ifld == 2) then
             long_name='soil thickness'; units = 'm'
          else if (ifld == 3) then
             long_name='saturated soil water content (porosity)';  units = 'mm3/mm3'
          else if (ifld == 4) then
             long_name='saturated soil matric potential'; units = 'mm'
          else if (ifld == 5) then
             long_name='slope of soil water retention curve'; units = 'unitless'
          else if (ifld == 6) then
             long_name='saturated hydraulic conductivity'; units = 'unitless'
          else
             call endrun(msg=' ERROR: bad 3D time-constant field index'//errMsg(__FILE__, __LINE__))
          end if
          if (tape(t)%dov2xy) then
             if (ldomain%isgrid2d) then
                call ncd_defvar(ncid=nfid(t), varname=trim(varnames(ifld)), xtype=tape(t)%ncprec,&
                     dim1name='lon', dim2name='lat', dim3name='levgrnd', &
                     long_name=long_name, units=units, missing_value=spval, fill_value=spval)
             else
                call ncd_defvar(ncid=nfid(t), varname=trim(varnames(ifld)), xtype=tape(t)%ncprec, &
                        dim1name=grlnd, dim2name='levgrnd', &
                     long_name=long_name, units=units, missing_value=spval, fill_value=spval)
             end if
          else
             call ncd_defvar(ncid=nfid(t), varname=trim(varnames(ifld)), xtype=tape(t)%ncprec, &
                  dim1name=namec, dim2name='levgrnd', &
                  long_name=long_name, units=units, missing_value=spval, fill_value=spval)
          end if
          call shr_string_listAppend(TimeConst3DVars,varnames(ifld))
       end do

    else if (mode == 'write') then

       allocate(histi(bounds%begc:bounds%endc,nlevgrnd), stat=ier)
       if (ier /= 0) then
          write(iulog,*) trim(subname),' ERROR: allocation error for histi'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       ! Write time constant fields

       if (tape(t)%dov2xy) then
          allocate(histo(bounds%begg:bounds%endg,nlevgrnd), stat=ier)
          if (ier /= 0) then
             write(iulog,*)  trim(subname),' ERROR: allocation error for histo'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if
       end if

       do ifld = 1,nflds

          ! WJS (10-25-11): Note about l2g_scale_type in the following: ZSOI & DZSOI are
          ! currently constant in space, except for urban points, so their scale type
          ! doesn't matter at the moment as long as it excludes urban points. I am using
          ! 'nonurb' so that the values are output everywhere where the fields are
          ! constant (i.e., everywhere except urban points). For the other fields, I am
          ! using 'veg' to be consistent with the l2g_scale_type that is now used for many
          ! of the 3-d time-variant fields; in theory, though, one might want versions of
          ! these variables output for different landunits.

          ! Field indices MUST match varnames array order above!
          if      (ifld == 1) then  ! ZSOI
             l2g_scale_type = 'nonurb'
          else if (ifld == 2) then  ! DZSOI
             l2g_scale_type = 'nonurb'
          else if (ifld == 3) then  ! WATSAT
             l2g_scale_type = 'veg'
          else if (ifld == 4) then  ! SUCSAT
             l2g_scale_type = 'veg'
          else if (ifld == 5) then  ! BSW
             l2g_scale_type = 'veg'
          else if (ifld == 6) then  ! HKSAT
             l2g_scale_type = 'veg'
          end if

          histi(:,:) = spval
          do lev = 1,nlevgrnd
             do c = bounds%begc,bounds%endc
                l = col_pp%landunit(c)
                   ! Field indices MUST match varnames array order above!
                   if (ifld ==1) histi(c,lev) = col_pp%z(c,lev)
                   if (ifld ==2) histi(c,lev) = col_pp%dz(c,lev)
                   if (ifld ==3) histi(c,lev) = watsat_col(c,lev)
                   if (ifld ==4) histi(c,lev) = sucsat_col(c,lev)
                   if (ifld ==5) histi(c,lev) = bsw_col(c,lev)
                   if (ifld ==6) histi(c,lev) = hksat_col(c,lev)
             end do
          end do
          if (tape(t)%dov2xy) then
             histo(:,:) = spval

             call c2g(bounds, nlevgrnd, &
                  histi(bounds%begc:bounds%endc, :), &
                  histo(bounds%begg:bounds%endg, :), &
                  c2l_scale_type='unity', l2g_scale_type=l2g_scale_type)

             if (ldomain%isgrid2d) then
                call ncd_io(varname=trim(varnames(ifld)), dim1name=grlnd, &
                     data=histo, ncid=nfid(t), flag='write')
             else
                call ncd_io(varname=trim(varnames(ifld)), dim1name=grlnd, &
                     data=histo, ncid=nfid(t), flag='write')
             end if
          else
             call ncd_io(varname=trim(varnames(ifld)), dim1name=namec, &
                  data=histi, ncid=nfid(t), flag='write')
          end if
       end do

       if (tape(t)%dov2xy) deallocate(histo)
       deallocate(histi)

    end if  ! (define/write mode

    if (mode == 'define') then
       do ifld = 1,nfldsl
          ! Field indices MUST match varnamesl array order above!
          if (ifld == 1) then
             long_name='lake layer node depth'; units = 'm'
          else if (ifld == 2) then
             long_name='lake layer thickness'; units = 'm'
          else
             call endrun(msg=' ERROR: bad 3D time-constant field index'//errMsg(__FILE__, __LINE__))
          end if
          if (tape(t)%dov2xy) then
             if (ldomain%isgrid2d) then
                call ncd_defvar(ncid=nfid(t), varname=trim(varnamesl(ifld)), xtype=tape(t)%ncprec,&
                     dim1name='lon', dim2name='lat', dim3name='levlak', &
                     long_name=long_name, units=units, missing_value=spval, fill_value=spval)
             else
                call ncd_defvar(ncid=nfid(t), varname=trim(varnamesl(ifld)), xtype=tape(t)%ncprec, &
                        dim1name=grlnd, dim2name='levlak', &
                     long_name=long_name, units=units, missing_value=spval, fill_value=spval)
             end if
          else
             call ncd_defvar(ncid=nfid(t), varname=trim(varnamesl(ifld)), xtype=tape(t)%ncprec, &
                  dim1name=namec, dim2name='levlak', &
                  long_name=long_name, units=units, missing_value=spval, fill_value=spval)
          end if
          call shr_string_listAppend(TimeConst3DVars,varnamesl(ifld))
       end do

    else if (mode == 'write') then

       allocate(histil(bounds%begc:bounds%endc,nlevlak), stat=ier)
       if (ier /= 0) then
          write(iulog,*) trim(subname),' ERROR: allocation error for histil'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       ! Write time constant fields

       if (tape(t)%dov2xy) then
          allocate(histol(bounds%begg:bounds%endg,nlevlak), stat=ier)
          if (ier /= 0) then
             write(iulog,*)  trim(subname),' ERROR: allocation error for histol'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if
       end if

       do ifld = 1,nfldsl
          histil(:,:) = spval
          do lev = 1,nlevlak
             do c = bounds%begc,bounds%endc
                l = col_pp%landunit(c)
                if (lun_pp%lakpoi(l)) then
                   ! Field indices MUST match varnamesl array order above!
                   if (ifld ==1) histil(c,lev) = col_pp%z_lake(c,lev) 
                   if (ifld ==2) histil(c,lev) = col_pp%dz_lake(c,lev)
                end if
             end do
          end do
          if (tape(t)%dov2xy) then
             histol(:,:) = spval
             call c2g(bounds, nlevlak, &
                  histil(bounds%begc:bounds%endc, :), &
                  histol(bounds%begg:bounds%endg, :), &
                  c2l_scale_type='unity', l2g_scale_type='lake')
             if (ldomain%isgrid2d) then
                call ncd_io(varname=trim(varnamesl(ifld)), dim1name=grlnd, &
                     data=histol, ncid=nfid(t), flag='write')
             else
                call ncd_io(varname=trim(varnamesl(ifld)), dim1name=grlnd, &
                     data=histol, ncid=nfid(t), flag='write')
             end if
          else
             call ncd_io(varname=trim(varnamesl(ifld)), dim1name=namec,  &
                  data=histil, ncid=nfid(t), flag='write')
          end if
       end do

       if (tape(t)%dov2xy) deallocate(histol)
       deallocate(histil)

    end if  ! (define/write mode

  end subroutine htape_timeconst3D

  !-----------------------------------------------------------------------
  subroutine htape_timeconst(t, mode)
    !
    ! !DESCRIPTION:
    ! Write time constant values to primary history tape.
    ! Issue the required netcdf wrapper calls to define the history file
    ! contents.
    !
    ! !USES:
    use elm_varcon      , only : zsoi, zlak, secspday
    use domainMod       , only : ldomain, lon1d, lat1d
    use clm_time_manager, only : get_nstep, get_curr_date, get_curr_time
    use clm_time_manager, only : get_ref_date, get_calendar, NO_LEAP_C, GREGORIAN_C
    use FatesInterfaceTypesMod, only : fates_hdim_levsclass
    use FatesInterfaceTypesMod, only : fates_hdim_pfmap_levscpf
    use FatesInterfaceTypesMod, only : fates_hdim_scmap_levscpf
    use FatesInterfaceTypesMod, only : fates_hdim_levcoage
    use FatesInterfaceTypesMod, only : fates_hdim_pfmap_levcapf
    use FatesInterfaceTypesMod, only : fates_hdim_camap_levcapf
    use FatesInterfaceTypesMod, only : fates_hdim_levage
    use FatesInterfaceTypesMod, only : fates_hdim_levpft
    use FatesInterfaceTypesMod, only : fates_hdim_scmap_levscag
    use FatesInterfaceTypesMod, only : fates_hdim_agmap_levscag
    use FatesInterfaceTypesMod, only : fates_hdim_levfuel
    use FatesInterfaceTypesMod, only : fates_hdim_levcwdsc
    use FatesInterfaceTypesMod, only : fates_hdim_levcan
    use FatesInterfaceTypesMod, only : fates_hdim_canmap_levcnlf
    use FatesInterfaceTypesMod, only : fates_hdim_lfmap_levcnlf
    use FatesInterfaceTypesMod, only : fates_hdim_canmap_levcnlfpf
    use FatesInterfaceTypesMod, only : fates_hdim_lfmap_levcnlfpf
    use FatesInterfaceTypesMod, only : fates_hdim_pftmap_levcnlfpf
    use FatesInterfaceTypesMod, only : fates_hdim_levheight
    use FatesInterfaceTypesMod, only : fates_hdim_scmap_levscagpft
    use FatesInterfaceTypesMod, only : fates_hdim_agmap_levscagpft
    use FatesInterfaceTypesMod, only : fates_hdim_pftmap_levscagpft
    use FatesInterfaceTypesMod, only : fates_hdim_agmap_levagepft
    use FatesInterfaceTypesMod, only : fates_hdim_pftmap_levagepft
    use FatesInterfaceTypesMod, only : fates_hdim_levelem
    use FatesInterfaceTypesMod, only : fates_hdim_elmap_levelpft
    use FatesInterfaceTypesMod, only : fates_hdim_pftmap_levelpft
    use FatesInterfaceTypesMod, only : fates_hdim_elmap_levelcwd
    use FatesInterfaceTypesMod, only : fates_hdim_cwdmap_levelcwd
    use FatesInterfaceTypesMod, only : fates_hdim_elmap_levelage
    use FatesInterfaceTypesMod, only : fates_hdim_agemap_levelage



    !
    ! !ARGUMENTS:
    integer, intent(in) :: t              ! tape index
    character(len=*), intent(in) :: mode  ! 'define' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: vid,n,i,j,m                ! indices
    integer :: nstep                      ! current step
    integer :: mcsec                      ! seconds of current date
    integer :: mdcur                      ! current day
    integer :: mscur                      ! seconds of current day
    integer :: mcdate                     ! current date
    integer :: yr,mon,day,nbsec           ! year,month,day,seconds components of a date
    integer :: hours,minutes,secs         ! hours,minutes,seconds of hh:mm:ss
    character(len= 10) :: basedate        ! base date (yyyymmdd)
    character(len=  8) :: basesec         ! base seconds
    character(len=  8) :: cdate           ! system date
    character(len=  8) :: ctime           ! system time
    real(r8):: time                       ! current time
    real(r8):: timedata(2)                ! time interval boundaries
    integer :: dim1id(1)                  ! netCDF dimension id
    integer :: dim2id(2)                  ! netCDF dimension id
    integer :: varid                      ! netCDF variable id
    character(len=max_chars) :: long_name ! variable long name
    character(len=max_namlen):: varname   ! variable name
    character(len=max_namlen):: units     ! variable units
    character(len=max_namlen):: cal       ! calendar from the time-manager
    character(len=max_namlen):: caldesc   ! calendar description to put on file
    character(len=256):: str              ! global attribute string
    real(r8), pointer :: histo(:,:)       ! temporary
    integer :: status
    real(r8) :: zsoi_1d(1)
    character(len=*),parameter :: subname = 'htape_timeconst'
    !-----------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    !***     Time constant grid variables only on first time-sample of file ***
    !-------------------------------------------------------------------------------

    if (tape(t)%ntimes == 1) then
       if (mode == 'define') then
          call ncd_defvar(varname='levgrnd', xtype=tape(t)%ncprec, &
                dim1name='levgrnd', &
                long_name='coordinate soil levels', units='m', ncid=nfid(t))
          call ncd_defvar(varname='levlak', xtype=tape(t)%ncprec, &
                dim1name='levlak', &
                long_name='coordinate lake levels', units='m', ncid=nfid(t))
          call ncd_defvar(varname='levdcmp', xtype=tape(t)%ncprec, dim1name='levdcmp', &
                long_name='coordinate soil levels', units='m', ncid=nfid(t))

          if(use_fates)then

             call ncd_defvar(varname='fates_levscls', xtype=tape(t)%ncprec, dim1name='fates_levscls', &
                   long_name='FATES diameter size class lower bound', units='cm', ncid=nfid(t))
             call ncd_defvar(varname='fates_scmap_levscag', xtype=ncd_int, dim1name='fates_levscag', &
                   long_name='FATES size-class map into size x patch age', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_agmap_levscag', xtype=ncd_int, dim1name='fates_levscag', &
                   long_name='FATES age-class map into size x patch age', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_pftmap_levscpf',xtype=ncd_int, dim1name='fates_levscpf', &
                   long_name='FATES pft index of the combined pft-size class dimension', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_scmap_levscpf',xtype=ncd_int, dim1name='fates_levscpf', &
                  long_name='FATES size index of the combined pft-size class dimension', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_levcacls', xtype=tape(t)%ncprec, dim1name='fates_levcacls', &
                  long_name='FATES cohort age class lower bound', units='years', ncid=nfid(t))
             call ncd_defvar(varname='fates_pftmap_levcapf',xtype=ncd_int, dim1name='fates_levcapf', &
                  long_name='FATES pft index of the combined pft-cohort age class dimension', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_camap_levcapf',xtype=ncd_int, dim1name='fates_levcapf', &
                  long_name='FATES cohort age index of the combined pft-cohort age dimension', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_levage',xtype=tape(t)%ncprec, dim1name='fates_levage', &
                   long_name='FATES patch age (yr)', ncid=nfid(t))
             call ncd_defvar(varname='fates_levpft',xtype=ncd_int, dim1name='fates_levpft', &
                   long_name='FATES pft number', ncid=nfid(t))
             call ncd_defvar(varname='fates_levfuel',xtype=ncd_int, dim1name='fates_levfuel', &
                   long_name='FATES fuel index', ncid=nfid(t))
             call ncd_defvar(varname='fates_levcwdsc',xtype=ncd_int, dim1name='fates_levcwdsc', &
                   long_name='FATES cwd size class', ncid=nfid(t))
             call ncd_defvar(varname='fates_levcan',xtype=ncd_int, dim1name='fates_levcan', &
                   long_name='FATES canopy level', ncid=nfid(t))
             call ncd_defvar(varname='fates_canmap_levcnlf',xtype=ncd_int, dim1name='fates_levcnlf', &
                   long_name='FATES canopy level of combined canopy-leaf dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_lfmap_levcnlf',xtype=ncd_int, dim1name='fates_levcnlf', &
                   long_name='FATES leaf level of combined canopy-leaf dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_canmap_levcnlfpf',xtype=ncd_int, dim1name='fates_levcnlfpf', &
                   long_name='FATES canopy level of combined canopy x leaf x pft dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_lfmap_levcnlfpf',xtype=ncd_int, dim1name='fates_levcnlfpf', &
                   long_name='FATES leaf level of combined canopy x leaf x pft dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_pftmap_levcnlfpf',xtype=ncd_int, dim1name='fates_levcnlfpf', &
                   long_name='FATES PFT level of combined canopy x leaf x pft dimension', ncid=nfid(t))
             call ncd_defvar(varname='fates_levheight',xtype=tape(t)%ncprec, dim1name='fates_levheight', &
                  long_name='FATES height (m)', ncid=nfid(t))
             call ncd_defvar(varname='fates_scmap_levscagpft', xtype=ncd_int, dim1name='fates_levscagpf', &
                  long_name='FATES size-class map into size x patch age x pft', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_agmap_levscagpft', xtype=ncd_int, dim1name='fates_levscagpf', &
                  long_name='FATES age-class map into size x patch age x pft', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_pftmap_levscagpft', xtype=ncd_int, dim1name='fates_levscagpf', &
                  long_name='FATES pft map into size x patch age x pft', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_pftmap_levagepft', xtype=ncd_int, dim1name='fates_levagepft', &
                  long_name='FATES pft map into patch age x pft', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_agmap_levagepft', xtype=ncd_int, dim1name='fates_levagepft', &
                   long_name='FATES age-class map into patch age x pft', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_levelem',xtype=ncd_int, dim1name='fates_levelem', &
                   long_name='FATES element (C,N,P,...) identifier', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_elmap_levelpft', xtype=ncd_int, dim1name='fates_levelpft', &
                   long_name='FATES element map into element x pft   ', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_pftmap_levelpft', xtype=ncd_int, dim1name='fates_levelpft', &
                   long_name='FATES pft map into element x pft', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_elmap_levelcwd', xtype=ncd_int, dim1name='fates_levelcwd', &
                   long_name='FATES element map into element x cwd', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_cwdmap_levelcwd', xtype=ncd_int, dim1name='fates_levelcwd', &
                   long_name='FATES cwd map into element x cwd', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_elmap_levelage', xtype=ncd_int, dim1name='fates_levelage', &
                   long_name='FATES element map into age x pft', units='-', ncid=nfid(t))
             call ncd_defvar(varname='fates_agemap_levelage', xtype=ncd_int, dim1name='fates_levelage', &
                   long_name='FATES element map into age x pft', units='-', ncid=nfid(t))

          end if

       elseif (mode == 'write') then
          if ( masterproc ) write(iulog, *) ' zsoi:',zsoi
          call ncd_io(varname='levgrnd', data=zsoi, ncid=nfid(t), flag='write')
          call ncd_io(varname='levlak' , data=zlak, ncid=nfid(t), flag='write')
          if (use_vertsoilc) then
             call ncd_io(varname='levdcmp', data=zsoi, ncid=nfid(t), flag='write')
          else
             zsoi_1d(1) = 1._r8
             call ncd_io(varname='levdcmp', data=zsoi_1d, ncid=nfid(t), flag='write')
          end if
          if(use_fates)then
             call ncd_io(varname='fates_scmap_levscag',data=fates_hdim_scmap_levscag, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_agmap_levscag',data=fates_hdim_agmap_levscag, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levscls',data=fates_hdim_levsclass, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levcacls',data=fates_hdim_levcoage, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_pftmap_levcapf',data=fates_hdim_pfmap_levcapf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_camap_levcapf',data=fates_hdim_camap_levcapf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_pftmap_levscpf',data=fates_hdim_pfmap_levscpf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_scmap_levscpf',data=fates_hdim_scmap_levscpf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levage',data=fates_hdim_levage, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levpft',data=fates_hdim_levpft, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levfuel',data=fates_hdim_levfuel, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levcwdsc',data=fates_hdim_levcwdsc, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levcan',data=fates_hdim_levcan, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_canmap_levcnlf',data=fates_hdim_canmap_levcnlf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_lfmap_levcnlf',data=fates_hdim_lfmap_levcnlf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_canmap_levcnlfpf',data=fates_hdim_canmap_levcnlfpf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_lfmap_levcnlfpf',data=fates_hdim_lfmap_levcnlfpf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_pftmap_levcnlfpf',data=fates_hdim_pftmap_levcnlfpf, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levheight',data=fates_hdim_levheight, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_scmap_levscagpft',data=fates_hdim_scmap_levscagpft, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_agmap_levscagpft',data=fates_hdim_agmap_levscagpft, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_pftmap_levscagpft',data=fates_hdim_pftmap_levscagpft, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_pftmap_levagepft',data=fates_hdim_pftmap_levagepft, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_agmap_levagepft',data=fates_hdim_agmap_levagepft, ncid=nfid(t), flag='write')
             call ncd_io(varname='fates_levelem',data=fates_hdim_levelem, ncid=nfid(t),flag='write')
             call ncd_io(varname='fates_elmap_levelpft',data=fates_hdim_elmap_levelpft, ncid=nfid(t),flag='write')
             call ncd_io(varname='fates_pftmap_levelpft',data=fates_hdim_pftmap_levelpft, ncid=nfid(t),flag='write')
             call ncd_io(varname='fates_elmap_levelcwd',data=fates_hdim_elmap_levelcwd, ncid=nfid(t),flag='write')
             call ncd_io(varname='fates_cwdmap_levelcwd',data=fates_hdim_cwdmap_levelcwd, ncid=nfid(t),flag='write')
             call ncd_io(varname='fates_elmap_levelage',data=fates_hdim_elmap_levelage, ncid=nfid(t),flag='write')
             call ncd_io(varname='fates_agemap_levelage',data=fates_hdim_agemap_levelage, ncid=nfid(t),flag='write')

          end if

       endif
    endif

    !-------------------------------------------------------------------------------
    !***     Time definition variables ***
    !-------------------------------------------------------------------------------

    ! For define mode -- only do this for first time-sample
    if (mode == 'define' .and. tape(t)%ntimes == 1) then
       call get_ref_date(yr, mon, day, nbsec)
       nstep = get_nstep()
       hours   = nbsec / 3600
       minutes = (nbsec - hours*3600) / 60
       secs    = (nbsec - hours*3600 - minutes*60)
       write(basedate,80) yr,mon,day
80     format(i4.4,'-',i2.2,'-',i2.2)
       write(basesec ,90) hours, minutes, secs
90     format(i2.2,':',i2.2,':',i2.2)

       dim1id(1) = time_dimid
       str = 'days since ' // basedate // " " // basesec
       call ncd_defvar(nfid(t), 'time', tape(t)%ncprec, 1, dim1id, varid, &
            long_name='time',units=str) 
       cal = get_calendar()
       if (      trim(cal) == NO_LEAP_C   )then
          caldesc = "noleap"
       else if ( trim(cal) == GREGORIAN_C )then
          caldesc = "gregorian"
       end if
       call ncd_putatt(nfid(t), varid, 'calendar', caldesc)
       call ncd_putatt(nfid(t), varid, 'bounds', 'time_bounds')

       dim1id(1) = time_dimid
       call ncd_defvar(nfid(t) , 'mcdate', ncd_int, 1, dim1id , varid, &
          long_name = 'current date (YYYYMMDD)')
       call ncd_defvar(nfid(t) , 'mcsec' , ncd_int, 1, dim1id , varid, &
          long_name = 'current seconds of current date', units='s')
       call ncd_defvar(nfid(t) , 'mdcur' , ncd_int, 1, dim1id , varid, &
          long_name = 'current day (from base day)')
       call ncd_defvar(nfid(t) , 'mscur' , ncd_int, 1, dim1id , varid, &
          long_name = 'current seconds of current day')
       call ncd_defvar(nfid(t) , 'nstep' , ncd_int, 1, dim1id , varid, &
          long_name = 'time step')

       dim2id(1) = hist_interval_dimid;  dim2id(2) = time_dimid
       call ncd_defvar(nfid(t), 'time_bounds', ncd_double, 2, dim2id, varid, &
          long_name = 'history time interval endpoints')

       dim2id(1) = strlen_dimid;  dim2id(2) = time_dimid
       call ncd_defvar(nfid(t), 'date_written', ncd_char, 2, dim2id, varid)
       call ncd_defvar(nfid(t), 'time_written', ncd_char, 2, dim2id, varid)

       if ( len_trim(TimeConst3DVars_Filename) > 0 )then
          call ncd_putatt(nfid(t), ncd_global, 'Time_constant_3Dvars_filename', &
                          trim(TimeConst3DVars_Filename))
       end if
       if ( len_trim(TimeConst3DVars)          > 0 )then
          call ncd_putatt(nfid(t), ncd_global, 'Time_constant_3Dvars',          &
                          trim(TimeConst3DVars))
       end if

    elseif (mode == 'write') then

       call get_curr_time (mdcur, mscur)
       call get_curr_date (yr, mon, day, mcsec)
       mcdate = yr*10000 + mon*100 + day
       nstep = get_nstep()

       call ncd_io('mcdate', mcdate, 'write', nfid(t), nt=tape(t)%ntimes)
       call ncd_io('mcsec' , mcsec , 'write', nfid(t), nt=tape(t)%ntimes)
       call ncd_io('mdcur' , mdcur , 'write', nfid(t), nt=tape(t)%ntimes)
       call ncd_io('mscur' , mscur , 'write', nfid(t), nt=tape(t)%ntimes)
       call ncd_io('nstep' , nstep , 'write', nfid(t), nt=tape(t)%ntimes)

       time = mdcur + mscur/secspday
       call ncd_io('time'  , time  , 'write', nfid(t), nt=tape(t)%ntimes)

       timedata(1) = tape(t)%begtime
       timedata(2) = time
       call ncd_io('time_bounds', timedata, 'write', nfid(t), nt=tape(t)%ntimes)

       call getdatetime (cdate, ctime)
       call ncd_io('date_written', cdate, 'write', nfid(t), nt=tape(t)%ntimes)

       call ncd_io('time_written', ctime, 'write', nfid(t), nt=tape(t)%ntimes)

    endif

    !-------------------------------------------------------------------------------
    !***     Grid definition variables ***
    !-------------------------------------------------------------------------------
    ! For define mode -- only do this for first time-sample
    if (mode == 'define' .and. tape(t)%ntimes == 1) then

       if (ldomain%isgrid2d) then
          call ncd_defvar(varname='lon', xtype=tape(t)%ncprec, dim1name='lon', &
              long_name='coordinate longitude', units='degrees_east', &
              ncid=nfid(t), missing_value=spval, fill_value=spval)
       else
          call ncd_defvar(varname='lon', xtype=tape(t)%ncprec, &
              dim1name=grlnd, &
              long_name='coordinate longitude', units='degrees_east', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       end if
       if (ldomain%isgrid2d) then
          call ncd_defvar(varname='lat', xtype=tape(t)%ncprec, dim1name='lat', &
              long_name='coordinate latitude', units='degrees_north', &
              ncid=nfid(t), missing_value=spval, fill_value=spval)
       else
          call ncd_defvar(varname='lat', xtype=tape(t)%ncprec, &
              dim1name=grlnd, &
              long_name='coordinate latitude', units='degrees_north', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       end if
       if (ldomain%isgrid2d) then
          call ncd_defvar(varname='area', xtype=tape(t)%ncprec, &
              dim1name='lon', dim2name='lat',&
              long_name='grid cell areas', units='km^2', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       else
          call ncd_defvar(varname='area', xtype=tape(t)%ncprec, &
              dim1name=grlnd, &
              long_name='grid cell areas', units='km^2', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       end if
       if (ldomain%isgrid2d) then
          call ncd_defvar(varname='topo', xtype=tape(t)%ncprec, &
              dim1name='lon', dim2name='lat',&
              long_name='grid cell topography', units='m', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       else
          call ncd_defvar(varname='topo', xtype=tape(t)%ncprec, &
              dim1name=grlnd, &
              long_name='grid cell topography', units='m', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       end if
       if (ldomain%isgrid2d) then
          call ncd_defvar(varname='landfrac', xtype=tape(t)%ncprec, &
              dim1name='lon', dim2name='lat', &
              long_name='land fraction', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       else
          call ncd_defvar(varname='landfrac', xtype=tape(t)%ncprec, &
              dim1name=grlnd, &
              long_name='land fraction', ncid=nfid(t), &
              missing_value=spval, fill_value=spval)
       end if
       if (ldomain%isgrid2d) then
          call ncd_defvar(varname='landmask', xtype=ncd_int, &
              dim1name='lon', dim2name='lat', &
              long_name='land/ocean mask (0.=ocean and 1.=land)', ncid=nfid(t), &
              imissing_value=ispval, ifill_value=ispval)
       else
          call ncd_defvar(varname='landmask', xtype=ncd_int, &
              dim1name=grlnd, &
              long_name='land/ocean mask (0.=ocean and 1.=land)', ncid=nfid(t), &
              imissing_value=ispval, ifill_value=ispval)
       end if
       if (ldomain%isgrid2d) then
          call ncd_defvar(varname='pftmask' , xtype=ncd_int, &
              dim1name='lon', dim2name='lat', &
              long_name='pft real/fake mask (0.=fake and 1.=real)', ncid=nfid(t), &
              imissing_value=ispval, ifill_value=ispval)
       else
          call ncd_defvar(varname='pftmask' , xtype=ncd_int, &
              dim1name=grlnd, &
              long_name='pft real/fake mask (0.=fake and 1.=real)', ncid=nfid(t), &
              imissing_value=ispval, ifill_value=ispval)
       end if

    else if (mode == 'write') then

       ! Most of this is constant and only needs to be done on tape(t)%ntimes=1
       ! But, some may change for dynamic PFT mode for example

       if (ldomain%isgrid2d) then
          call ncd_io(varname='lon', data=lon1d, ncid=nfid(t), flag='write')
          call ncd_io(varname='lat', data=lat1d, ncid=nfid(t), flag='write')
       else
          call ncd_io(varname='lon', data=ldomain%lonc, dim1name=grlnd, ncid=nfid(t), flag='write')
          call ncd_io(varname='lat', data=ldomain%latc, dim1name=grlnd, ncid=nfid(t), flag='write')
       end if
       call ncd_io(varname='topo'    , data=ldomain%topo, dim1name=grlnd, ncid=nfid(t), flag='write')
       call ncd_io(varname='area'    , data=ldomain%area, dim1name=grlnd, ncid=nfid(t), flag='write')
       call ncd_io(varname='landfrac', data=ldomain%frac, dim1name=grlnd, ncid=nfid(t), flag='write')
       call ncd_io(varname='landmask', data=ldomain%mask, dim1name=grlnd, ncid=nfid(t), flag='write')
       call ncd_io(varname='pftmask' , data=ldomain%pftm, dim1name=grlnd, ncid=nfid(t), flag='write')

    end if  ! (define/write mode

  end subroutine htape_timeconst

  !-----------------------------------------------------------------------
  subroutine hfields_write(t, mode)
    !
    ! !DESCRIPTION:
    ! Write history tape.  Issue the call to write the variable.
    !
    ! !USES:
    use domainMod , only : ldomain
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t                ! tape index
    character(len=*), intent(in) :: mode    ! 'define' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: f                         ! field index
    integer :: k                         ! 1d index
    integer :: topo,c,l,p                ! indices
    integer :: beg1d_out                 ! on-node 1d hbuf pointer start index
    integer :: end1d_out                 ! on-node 1d hbuf pointer end index
    integer :: num1d_out                 ! size of hbuf first dimension (overall all nodes)
    integer :: num2d                     ! hbuf second dimension size
    integer :: nt                        ! time index
    integer :: ier                       ! error status
    integer :: numdims                   ! number of dimensions
    character(len=1)         :: avgflag  ! time averaging flag
    character(len=max_chars) :: long_name! long name
    character(len=max_chars) :: units    ! units
    character(len=max_namlen):: varname  ! variable name
    character(len=32) :: avgstr          ! time averaging type
    character(len=hist_dim_name_length)  :: type1d_out      ! history output 1d type
    character(len=hist_dim_name_length)  :: type2d          ! history output 2d type
    character(len=32) :: dim1name        ! temporary
    character(len=32) :: dim2name        ! temporary
    real(r8), pointer :: histo(:,:)      ! temporary
    real(r8), pointer :: hist1do(:)      ! temporary
    character(len=*),parameter :: subname = 'hfields_write'
!-----------------------------------------------------------------------

    ! Write/define 1d topological info

    if (.not. tape(t)%dov2xy) then
       if (mode == 'define') then
          call hfields_1dinfo(t, mode='define')
       else if (mode == 'write') then
          call hfields_1dinfo(t, mode='write')
       end if
    end if

    ! Define time-dependent variables create variables and attributes for field list

    do f = 1,tape(t)%nflds

       ! Set history field variables

       varname    = tape(t)%hlist(f)%field%name
       long_name  = tape(t)%hlist(f)%field%long_name
       units      = tape(t)%hlist(f)%field%units
       avgflag    = tape(t)%hlist(f)%avgflag
       type1d_out = tape(t)%hlist(f)%field%type1d_out
       beg1d_out  = tape(t)%hlist(f)%field%beg1d_out
       end1d_out  = tape(t)%hlist(f)%field%end1d_out
       num1d_out  = tape(t)%hlist(f)%field%num1d_out
       type2d     = tape(t)%hlist(f)%field%type2d 
       numdims    = tape(t)%hlist(f)%field%numdims
       num2d      = tape(t)%hlist(f)%field%num2d
       nt         = tape(t)%ntimes

       if (mode == 'define') then

          select case (avgflag)
          case ('A')
             avgstr = 'mean'
          case ('I')
             avgstr = 'instantaneous'
          case ('X')
             avgstr = 'maximum'
          case ('M')
             avgstr = 'minimum'
          case default
             write(iulog,*) trim(subname),' ERROR: unknown time averaging flag (avgflag)=',avgflag
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select

          if (type1d_out == grlnd) then
             if (ldomain%isgrid2d) then
                dim1name = 'lon'      ; dim2name = 'lat'
             else
                dim1name = trim(grlnd); dim2name = 'undefined'
             end if
          else
             dim1name = type1d_out ; dim2name = 'undefined'
          endif

          if (dim2name == 'undefined') then
             if (numdims == 1) then
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name='time', &
                     long_name=long_name, units=units, cell_method=avgstr, &
                     missing_value=spval, fill_value=spval)
             else
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name=type2d, dim3name='time', &
                     long_name=long_name, units=units, cell_method=avgstr, &
                     missing_value=spval, fill_value=spval)
             end if
          else
             if (numdims == 1) then
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name=dim2name, dim3name='time', &
                     long_name=long_name, units=units, cell_method=avgstr, &
                     missing_value=spval, fill_value=spval)
             else
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name=dim1name, dim2name=dim2name, dim3name=type2d, dim4name='time', &
                     long_name=long_name, units=units, cell_method=avgstr, &
                     missing_value=spval, fill_value=spval)
             end if
          endif

       else if (mode == 'write') then

          ! Determine output buffer

          histo => tape(t)%hlist(f)%hbuf

          ! Allocate dynamic memory

          if (numdims == 1) then
             allocate(hist1do(beg1d_out:end1d_out), stat=ier)
             if (ier /= 0) then
                write(iulog,*) trim(subname),' ERROR: allocation'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end if
             hist1do(beg1d_out:end1d_out) = histo(beg1d_out:end1d_out,1)
          end if

          ! Write history output.  Always output land and ocean runoff on xy grid.

          if (numdims == 1) then
             call ncd_io(flag='write', varname=varname, &
                  dim1name=type1d_out, data=hist1do, ncid=nfid(t), nt=nt)
          else
             call ncd_io(flag='write', varname=varname, &
                  dim1name=type1d_out, data=histo, ncid=nfid(t), nt=nt)
          end if


          ! Deallocate dynamic memory

          if (numdims == 1) then
             deallocate(hist1do)
          end if

       end if

   end do

  end subroutine hfields_write

  !-----------------------------------------------------------------------
  subroutine hfields_1dinfo(t, mode)
    !
    ! !DESCRIPTION:
    ! Write/define 1d info for history tape.
    !
    ! !USES:
    use decompMod   , only : ldecomp
    use domainMod   , only : ldomain, ldomain
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t                ! tape index
    character(len=*), intent(in) :: mode    ! 'define' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: f                         ! field index
    integer :: k                         ! 1d index
    integer :: g,c,l,topo,p              ! indices
    integer :: ier                       ! errir status
    real(r8), pointer :: rgarr(:)        ! temporary
    real(r8), pointer :: rcarr(:)        ! temporary
    real(r8), pointer :: rlarr(:)        ! temporary
    real(r8), pointer :: rparr(:)        ! temporary
    integer , pointer :: igarr(:)        ! temporary
    integer , pointer :: icarr(:)        ! temporary
    integer , pointer :: ilarr(:)        ! temporary
    integer , pointer :: iparr(:)        ! temporary
    type(file_desc_t) :: ncid            ! netcdf file
    type(bounds_type) :: bounds          
    character(len=*),parameter :: subname = 'hfields_1dinfo'
!-----------------------------------------------------------------------

    call get_proc_bounds(bounds)

    ncid = nfid(t)

    if (mode == 'define') then

          ! Define gridcell info

          call ncd_defvar(varname='grid1d_lon', xtype=ncd_double, dim1name=nameg, &
               long_name='gridcell longitude', units='degrees_east', ncid=ncid)

          call ncd_defvar(varname='grid1d_lat', xtype=ncd_double,  dim1name=nameg, &
               long_name='gridcell latitude', units='degrees_north', ncid=ncid)

          call ncd_defvar(varname='grid1d_ixy', xtype=ncd_int, dim1name=nameg, &
               long_name='2d longitude index of corresponding gridcell', ncid=ncid)

          call ncd_defvar(varname='grid1d_jxy', xtype=ncd_int, dim1name=nameg, &
               long_name='2d latitude index of corresponding gridcell', ncid=ncid)

          ! Define landunit info

          call ncd_defvar(varname='land1d_lon', xtype=ncd_double, dim1name=namel, &
               long_name='landunit longitude', units='degrees_east', ncid=ncid)

          call ncd_defvar(varname='land1d_lat', xtype=ncd_double, dim1name=namel, &
               long_name='landunit latitude', units='degrees_north', ncid=ncid)

          call ncd_defvar(varname='land1d_ixy', xtype=ncd_int, dim1name=namel, &
               long_name='2d longitude index of corresponding landunit', ncid=ncid)

          call ncd_defvar(varname='land1d_jxy', xtype=ncd_int, dim1name=namel, &
               long_name='2d latitude index of corresponding landunit', ncid=ncid)

          ! --- EBK Do NOT write out indices that are incorrect 4/1/2011 --- Bug 1310
          !call ncd_defvar(varname='land1d_gi', xtype=ncd_int, dim1name='landunit', &
          !     long_name='1d grid index of corresponding landunit', ncid=ncid)
          ! ----------------------------------------------------------------

          call ncd_defvar(varname='land1d_wtgcell', xtype=ncd_double, dim1name=namel, &
               long_name='landunit weight relative to corresponding gridcell', ncid=ncid)

          call ncd_defvar(varname='land1d_ityplunit', xtype=ncd_int, dim1name=namel, &
               long_name='landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)', &
                  ncid=ncid)

          call ncd_defvar(varname='land1d_active', xtype=ncd_log, dim1name=namel, &
               long_name='true => do computations on this landunit', ncid=ncid)

          ! Define column info

          call ncd_defvar(varname='cols1d_lon', xtype=ncd_double, dim1name=namec, &
               long_name='column longitude', units='degrees_east', ncid=ncid)

          call ncd_defvar(varname='cols1d_lat', xtype=ncd_double, dim1name=namec, &
               long_name='column latitude', units='degrees_north', ncid=ncid)

          call ncd_defvar(varname='cols1d_ixy', xtype=ncd_int, dim1name=namec, &
               long_name='2d longitude index of corresponding column', ncid=ncid)

          call ncd_defvar(varname='cols1d_jxy', xtype=ncd_int, dim1name=namec, &
               long_name='2d latitude index of corresponding column', ncid=ncid)

          ! --- EBK Do NOT write out indices that are incorrect 4/1/2011 --- Bug 1310
          !call ncd_defvar(varname='cols1d_gi', xtype=ncd_int, dim1name='column', &
          !     long_name='1d grid index of corresponding column', ncid=ncid)

          !call ncd_defvar(varname='cols1d_li', xtype=ncd_int, dim1name='column', &
          !     long_name='1d landunit index of corresponding column', ncid=ncid)
          ! ----------------------------------------------------------------

          call ncd_defvar(varname='cols1d_wtgcell', xtype=ncd_double, dim1name=namec, &
               long_name='column weight relative to corresponding gridcell', ncid=ncid)

          call ncd_defvar(varname='cols1d_wtlunit', xtype=ncd_double, dim1name=namec, &
               long_name='column weight relative to corresponding landunit', ncid=ncid)

          call ncd_defvar(varname='cols1d_itype_lunit', xtype=ncd_int, dim1name=namec, &
               long_name='column landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)', &
                  ncid=ncid)

          call ncd_defvar(varname='cols1d_active', xtype=ncd_log, dim1name=namec, &
               long_name='true => do computations on this column', ncid=ncid)

          ! Define pft info

          call ncd_defvar(varname='pfts1d_lon', xtype=ncd_double, dim1name=namep, &
               long_name='pft longitude', units='degrees_east', ncid=ncid)

          call ncd_defvar(varname='pfts1d_lat', xtype=ncd_double, dim1name=namep, &
               long_name='pft latitude', units='degrees_north', ncid=ncid)

          call ncd_defvar(varname='pfts1d_ixy', xtype=ncd_int, dim1name=namep, &
               long_name='2d longitude index of corresponding pft', ncid=ncid)

          call ncd_defvar(varname='pfts1d_jxy', xtype=ncd_int, dim1name=namep, &
               long_name='2d latitude index of corresponding pft', ncid=ncid)

          ! --- EBK Do NOT write out indices that are incorrect 4/1/2011 --- Bug 1310
          !call ncd_defvar(varname='pfts1d_gi', xtype=ncd_int, dim1name='pft', &
          !     long_name='1d grid index of corresponding pft', ncid=ncid)

          !call ncd_defvar(varname='pfts1d_li', xtype=ncd_int, dim1name='pft', &
          !     long_name='1d landunit index of corresponding pft', ncid=ncid)

          !call ncd_defvar(varname='pfts1d_ci', xtype=ncd_int, dim1name='pft', &
          !     long_name='1d column index of corresponding pft', ncid=ncid)
          ! ----------------------------------------------------------------

          call ncd_defvar(varname='pfts1d_wtgcell', xtype=ncd_double, dim1name=namep, &
               long_name='pft weight relative to corresponding gridcell', ncid=ncid)

          call ncd_defvar(varname='pfts1d_wtlunit', xtype=ncd_double, dim1name=namep, &
               long_name='pft weight relative to corresponding landunit', ncid=ncid)

          call ncd_defvar(varname='pfts1d_wtcol', xtype=ncd_double, dim1name=namep, &
               long_name='pft weight relative to corresponding column', ncid=ncid)

          call ncd_defvar(varname='pfts1d_itype_veg', xtype=ncd_int, dim1name=namep, &
               long_name='pft vegetation type', ncid=ncid)

          call ncd_defvar(varname='pfts1d_itype_lunit', xtype=ncd_int, dim1name=namep, &
               long_name='pft landunit type (vegetated,urban,lake,wetland,glacier or glacier_mec)',  &
                  ncid=ncid)

          call ncd_defvar(varname='pfts1d_active', xtype=ncd_log, dim1name=namep, &
               long_name='true => do computations on this pft', ncid=ncid)

    else if (mode == 'write') then

       ! Determine bounds

       allocate(&
            rgarr(bounds%begg:bounds%endg),&
            rlarr(bounds%begl:bounds%endl),&
            rcarr(bounds%begc:bounds%endc),&
            rparr(bounds%begp:bounds%endp),&
            stat=ier)
       if (ier /= 0) then
          call endrun(msg=' hfields_1dinfo allocation error of rarrs'//errMsg(__FILE__, __LINE__))
       end if

       allocate(&
            igarr(bounds%begg:bounds%endg),&
            ilarr(bounds%begl:bounds%endl),&
            icarr(bounds%begc:bounds%endc),&
            iparr(bounds%begp:bounds%endp),stat=ier)
       if (ier /= 0) then
          call endrun(msg=' hfields_1dinfo allocation error of iarrs'//errMsg(__FILE__, __LINE__))
       end if

       ! Write gridcell info

       call ncd_io(varname='grid1d_lon', data=grc_pp%londeg, dim1name=nameg, ncid=ncid, flag='write')
       call ncd_io(varname='grid1d_lat', data=grc_pp%latdeg, dim1name=nameg, ncid=ncid, flag='write')
       do g = bounds%begg,bounds%endg
         igarr(g)= mod(ldecomp%gdc2glo(g)-1,ldomain%ni) + 1
       enddo
       call ncd_io(varname='grid1d_ixy', data=igarr      , dim1name=nameg, ncid=ncid, flag='write')
       do g = bounds%begg,bounds%endg
         igarr(g)= (ldecomp%gdc2glo(g) - 1)/ldomain%ni + 1
       enddo
       call ncd_io(varname='grid1d_jxy', data=igarr      , dim1name=nameg, ncid=ncid, flag='write')

       ! Write landunit info

       do l = bounds%begl,bounds%endl
         rlarr(l) = grc_pp%londeg(lun_pp%gridcell(l))
       enddo
       call ncd_io(varname='land1d_lon', data=rlarr, dim1name=namel, ncid=ncid, flag='write')
       do l = bounds%begl,bounds%endl
         rlarr(l) = grc_pp%latdeg(lun_pp%gridcell(l))
       enddo
       call ncd_io(varname='land1d_lat', data=rlarr, dim1name=namel, ncid=ncid, flag='write')
       do l= bounds%begl,bounds%endl
         ilarr(l) = mod(ldecomp%gdc2glo(lun_pp%gridcell(l))-1,ldomain%ni) + 1
       enddo
       call ncd_io(varname='land1d_ixy', data=ilarr, dim1name=namel, ncid=ncid, flag='write')
       do l = bounds%begl,bounds%endl
         ilarr(l) = (ldecomp%gdc2glo(lun_pp%gridcell(l))-1)/ldomain%ni + 1
       enddo
       call ncd_io(varname='land1d_jxy'      , data=ilarr        , dim1name=namel, ncid=ncid, flag='write')
       ! --- EBK Do NOT write out indices that are incorrect 4/1/2011 Bug 1310
       !call ncd_io(varname='land1d_gi'       , data=lun_pp%gridcell, dim1name=namel, ncid=ncid, flag='write')
       ! ----------------------------------------------------------------
       call ncd_io(varname='land1d_wtgcell'  , data=lun_pp%wtgcell , dim1name=namel, ncid=ncid, flag='write')
       call ncd_io(varname='land1d_ityplunit', data=lun_pp%itype   , dim1name=namel, ncid=ncid, flag='write')
       call ncd_io(varname='land1d_active'   , data=lun_pp%active  , dim1name=namel, ncid=ncid, flag='write')

       ! Write column info

       do c = bounds%begc,bounds%endc
         rcarr(c) = grc_pp%londeg(col_pp%gridcell(c))
       enddo
       call ncd_io(varname='cols1d_lon', data=rcarr, dim1name=namec, ncid=ncid, flag='write')
       do c = bounds%begc,bounds%endc
         rcarr(c) = grc_pp%latdeg(col_pp%gridcell(c))
       enddo
       call ncd_io(varname='cols1d_lat', data=rcarr, dim1name=namec, ncid=ncid, flag='write')
       do c = bounds%begc,bounds%endc
         icarr(c) = mod(ldecomp%gdc2glo(col_pp%gridcell(c))-1,ldomain%ni) + 1
       enddo
       call ncd_io(varname='cols1d_ixy', data=icarr, dim1name=namec, ncid=ncid, flag='write')
       do c = bounds%begc,bounds%endc
         icarr(c) = (ldecomp%gdc2glo(col_pp%gridcell(c))-1)/ldomain%ni + 1
       enddo
       call ncd_io(varname='cols1d_jxy'    , data=icarr         ,dim1name=namec, ncid=ncid, flag='write')
       ! --- EBK Do NOT write out indices that are incorrect 4/1/2011 Bug 1310
       !call ncd_io(varname='cols1d_gi'     , data=col_pp%gridcell, dim1name=namec, ncid=ncid, flag='write')
       !call ncd_io(varname='cols1d_li'     , data=col_pp%landunit, dim1name=namec, ncid=ncid, flag='write')
       ! ----------------------------------------------------------------
       call ncd_io(varname='cols1d_wtgcell', data=col_pp%wtgcell , dim1name=namec, ncid=ncid, flag='write')
       call ncd_io(varname='cols1d_wtlunit', data=col_pp%wtlunit , dim1name=namec, ncid=ncid, flag='write')
       do c = bounds%begc,bounds%endc
         icarr(c) = lun_pp%itype(col_pp%landunit(c))
       enddo
       call ncd_io(varname='cols1d_itype_lunit', data=icarr    , dim1name=namec, ncid=ncid, flag='write')
       call ncd_io(varname='cols1d_active' , data=col_pp%active  , dim1name=namec, ncid=ncid, flag='write')

       ! Write pft info

       do p = bounds%begp,bounds%endp
         rparr(p) = grc_pp%londeg(veg_pp%gridcell(p))
       enddo
       call ncd_io(varname='pfts1d_lon', data=rparr, dim1name=namep, ncid=ncid, flag='write')
       do p = bounds%begp,bounds%endp
         rparr(p) = grc_pp%latdeg(veg_pp%gridcell(p))
       enddo
       call ncd_io(varname='pfts1d_lat', data=rparr, dim1name=namep, ncid=ncid, flag='write')
       do p = bounds%begp,bounds%endp
         iparr(p) = mod(ldecomp%gdc2glo(veg_pp%gridcell(p))-1,ldomain%ni) + 1
       enddo
       call ncd_io(varname='pfts1d_ixy', data=iparr, dim1name=namep, ncid=ncid, flag='write')
       do p = bounds%begp,bounds%endp
         iparr(p) = (ldecomp%gdc2glo(veg_pp%gridcell(p))-1)/ldomain%ni + 1
       enddo
       call ncd_io(varname='pfts1d_jxy'      , data=iparr        , dim1name=namep, ncid=ncid, flag='write')
       ! --- EBK Do NOT write out indices that are incorrect 4/1/2011 --- Bug 1310
       !call ncd_io(varname='pfts1d_gi'       , data=veg_pp%gridcell, dim1name=namep, ncid=ncid, flag='write')
       !call ncd_io(varname='pfts1d_li'       , data=veg_pp%landunit, dim1name=namep, ncid=ncid, flag='write')
       !call ncd_io(varname='pfts1d_ci'       , data=veg_pp%column  , dim1name=namep, ncid=ncid, flag='write')
       ! ----------------------------------------------------------------
       call ncd_io(varname='pfts1d_wtgcell'  , data=veg_pp%wtgcell , dim1name=namep, ncid=ncid, flag='write')
       call ncd_io(varname='pfts1d_wtlunit'  , data=veg_pp%wtlunit , dim1name=namep, ncid=ncid, flag='write')
       call ncd_io(varname='pfts1d_wtcol'    , data=veg_pp%wtcol   , dim1name=namep, ncid=ncid, flag='write')
       call ncd_io(varname='pfts1d_itype_veg', data=veg_pp%itype   , dim1name=namep, ncid=ncid, flag='write')

       do p = bounds%begp,bounds%endp
          iparr(p) = lun_pp%itype(veg_pp%landunit(p))
       enddo
       call ncd_io(varname='pfts1d_itype_lunit', data=iparr      , dim1name=namep, ncid=ncid, flag='write')
       call ncd_io(varname='pfts1d_active'   , data=veg_pp%active  , dim1name=namep, ncid=ncid, flag='write')

       deallocate(rgarr,rlarr,rcarr,rparr)
       deallocate(igarr,ilarr,icarr,iparr)

    end if

  end subroutine hfields_1dinfo

  !-----------------------------------------------------------------------
  subroutine hist_htapes_wrapup( rstwr, nlend, bounds, &
       watsat_col, sucsat_col, bsw_col, hksat_col)
    !
    ! !DESCRIPTION:
    ! Write history tape(s)
    ! Determine if next time step is beginning of history interval and if so:
    !   increment the current time sample counter, open a new history file
    !   and if needed (i.e., when ntim = 1), write history data to current
    !   history file, reset field accumulation counters to zero.
    ! If primary history file is full or at the last time step of the simulation,
    !   write restart dataset and close all history fiels.
    ! If history file is full or at the last time step of the simulation:
    !   close history file
    !   and reset time sample counter to zero if file is full.
    ! Daily-averaged data for the first day in September are written on
    !   date = 00/09/02 with mscur = 0.
    ! Daily-averaged data for the first day in month mm are written on
    !   date = yyyy/mm/02 with mscur = 0.
    ! Daily-averaged data for the 30th day (last day in September) are written
    !   on date = 0000/10/01 mscur = 0.
    ! Daily-averaged data for the last day in month mm are written on
    !   date = yyyy/mm+1/01 with mscur = 0.
    !
    ! !USES:
    use clm_time_manager, only : get_nstep, get_curr_date, get_curr_time, get_prev_date
    use elm_varcon      , only : secspday
    use perf_mod        , only : t_startf, t_stopf
    use elm_varpar      , only : nlevgrnd
    !
    ! !ARGUMENTS:
    logical, intent(in) :: rstwr    ! true => write restart file this step
    logical, intent(in) :: nlend    ! true => end of run on this step
    type(bounds_type) , intent(in) :: bounds           
    real(r8)          , intent(in) :: watsat_col( bounds%begc:,1: ) 
    real(r8)          , intent(in) :: sucsat_col( bounds%begc:,1: ) 
    real(r8)          , intent(in) :: bsw_col( bounds%begc:,1: ) 
    real(r8)          , intent(in) :: hksat_col( bounds%begc:,1: ) 
    !
    ! !LOCAL VARIABLES:
    integer :: t                          ! tape index
    integer :: f                          ! field index
    integer :: ier                        ! error code
    integer :: nstep                      ! current step
    integer :: day                        ! current day (1 -> 31)
    integer :: mon                        ! current month (1 -> 12)
    integer :: yr                         ! current year (0 -> ...)
    integer :: mdcur                      ! current day
    integer :: mscur                      ! seconds of current day
    integer :: mcsec                      ! current time of day [seconds]
    integer :: daym1                      ! nstep-1 day (1 -> 31)
    integer :: monm1                      ! nstep-1 month (1 -> 12)
    integer :: yrm1                       ! nstep-1 year (0 -> ...)
    integer :: mcsecm1                    ! nstep-1 time of day [seconds]
    real(r8):: time                       ! current time
    character(len=256) :: str             ! global attribute string
    logical :: if_stop                    ! true => last time step of run
    logical, save :: do_3Dtconst = .true. ! true => write out 3D time-constant data
    character(len=*),parameter :: subname = 'hist_htapes_wrapup'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(watsat_col) == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sucsat_col) == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bsw_col)    == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(hksat_col)  == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

    ! get current step

    nstep = get_nstep()

    ! Set calendar for current time step

    call get_curr_date (yr, mon, day, mcsec)
    call get_curr_time (mdcur, mscur)
    time = mdcur + mscur/secspday

    ! Set calendar for current for previous time step

    call get_prev_date (yrm1, monm1, daym1, mcsecm1)

    ! Loop over active history tapes, create new history files if necessary
    ! and write data to history files if end of history interval.
    do t = 1, ntapes

       ! Skip nstep=0 if monthly average

       if (nstep==0 .and. tape(t)%nhtfrq==0) cycle

       ! Determine if end of history interval
       tape(t)%is_endhist = .false.
       if (tape(t)%nhtfrq==0) then   !monthly average
          if (mon /= monm1) tape(t)%is_endhist = .true.
       else
          if (mod(nstep,tape(t)%nhtfrq) == 0) tape(t)%is_endhist = .true.
       end if

       ! If end of history interval

       if (tape(t)%is_endhist) then

          ! Normalize history buffer if time averaged

          call hfields_normalize(t)

          ! Increment current time sample counter.

          tape(t)%ntimes = tape(t)%ntimes + 1

          ! Create history file if appropriate and build time comment

          ! If first time sample, generate unique history file name, open file,
          ! define dims, vars, etc.


          if (tape(t)%ntimes == 1) then
             call t_startf('hist_htapes_wrapup_define')
             locfnh(t) = set_hist_filename (hist_freq=tape(t)%nhtfrq, &
                                            hist_mfilt=tape(t)%mfilt, hist_file=t)
             if (masterproc) then
                write(iulog,*) trim(subname),' : Creating history file ', trim(locfnh(t)), &
                     ' at nstep = ',get_nstep()
                write(iulog,*)'calling htape_create for file t = ',t
             endif
             call htape_create (t)

             ! Define time-constant field variables
             call htape_timeconst(t, mode='define')

             ! Define 3D time-constant field variables only to first primary tape
             if ( do_3Dtconst .and. t == 1 ) then
                call htape_timeconst3D(t, &
                     bounds, watsat_col, sucsat_col, bsw_col, hksat_col, mode='define')
                TimeConst3DVars_Filename = trim(locfnh(t))
             end if

             ! Define model field variables
             call hfields_write(t, mode='define')

             ! Exit define model
             call ncd_enddef(nfid(t))
             call t_stopf('hist_htapes_wrapup_define')
          endif

          call t_startf('hist_htapes_wrapup_tconst')
          ! Write time constant history variables
          call htape_timeconst(t, mode='write')

          ! Write 3D time constant history variables only to first primary tape
          if ( do_3Dtconst .and. t == 1 .and. tape(t)%ntimes == 1 )then
             call htape_timeconst3D(t, &
                  bounds, watsat_col, sucsat_col, bsw_col, hksat_col, mode='write')
             do_3Dtconst = .false.
          end if

          if (masterproc) then
             write(iulog,*)
             write(iulog,*) trim(subname),' : Writing current time sample to local history file ', &
                  trim(locfnh(t)),' at nstep = ',get_nstep(), &
                  ' for history time interval beginning at ', tape(t)%begtime, &
                  ' and ending at ',time
             write(iulog,*)
             call shr_sys_flush(iulog)
          endif

          ! Update beginning time of next interval
          tape(t)%begtime = time
          call t_stopf('hist_htapes_wrapup_tconst')

          ! Write history time samples
          call t_startf('hist_htapes_wrapup_write')
          call hfields_write(t, mode='write')
          call t_stopf('hist_htapes_wrapup_write')

          ! Zero necessary history buffers
          call hfields_zero(t)

       end if

    end do  ! end loop over history tapes

    ! Determine if file needs to be closed

    call hist_do_disp (ntapes, tape(:)%ntimes, tape(:)%mfilt, if_stop, if_disphist, rstwr, nlend)

    ! Close open history file
    ! Auxilary files may have been closed and saved off without being full,
    ! must reopen the files

    do t = 1, ntapes
       if (if_disphist(t)) then
          if (tape(t)%ntimes /= 0) then
             if (masterproc) then
                write(iulog,*)
                write(iulog,*)  trim(subname),' : Closing local history file ',&
                     trim(locfnh(t)),' at nstep = ', get_nstep()
                write(iulog,*)
             endif
	     call ncd_pio_closefile(nfid(t))
             if (.not.if_stop .and. (tape(t)%ntimes/=tape(t)%mfilt)) then
                call ncd_pio_openfile (nfid(t), trim(locfnh(t)), ncd_write)
             end if
          else
             if (masterproc) then
                write(iulog,*) trim(subname),' : history tape ',t,': no open file to close'
             end if
          endif
       endif
    end do

    ! Reset number of time samples to zero if file is full 
    
    do t = 1, ntapes
       if (if_disphist(t) .and. tape(t)%ntimes==tape(t)%mfilt) then
          tape(t)%ntimes = 0
       end if
    end do
    
  end subroutine hist_htapes_wrapup

  !-----------------------------------------------------------------------
  subroutine hist_restart_ncd (bounds, ncid, flag, rdate)
    !
    ! !DESCRIPTION:
    ! Read/write history file restart data.
    ! If the current history file(s) are not full, file(s) are opened
    ! so that subsequent time samples are added until the file is full.
    ! A new history file is used on a branch run.
    !
    ! !USES:
    use elm_varctl      , only : nsrest, caseid, inst_suffix, nsrStartup, nsrBranch
    use fileutils       , only : getfil
    use domainMod       , only : ldomain
    use elm_varpar      , only : nlevgrnd, nlevlak, numrad, nlevdecomp_full, nmonth
    use clm_time_manager, only : is_restart
    use restUtilMod     , only : iflag_skip
    use pio
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds  
    type(file_desc_t), intent(inout) :: ncid     ! netcdf file
    character(len=*) , intent(in)    :: flag     !'read' or 'write'
    character(len=*) , intent(in), optional :: rdate    ! restart file time stamp for name
    !
    ! !LOCAL VARIABLES:
    integer :: max_nflds                     ! Max number of fields
    integer :: num1d,beg1d,end1d             ! 1d size, beginning and ending indices
    integer :: num1d_out,beg1d_out,end1d_out ! 1d size, beginning and ending indices
    integer :: num2d                         ! 2d size (e.g. number of vertical levels)
    integer :: numa                 ! total number of atm cells across all processors
    integer :: numg                 ! total number of gridcells across all processors
    integer :: numt                 ! total number of topounits across all processors
    integer :: numl                 ! total number of landunits across all processors
    integer :: numc                 ! total number of columns across all processors
    integer :: nump                 ! total number of pfts across all processors
    character(len=max_namlen) :: name            ! variable name
    character(len=max_namlen) :: name_acc        ! accumulator variable name
    character(len=max_namlen) :: long_name       ! long name of variable
    character(len=max_chars)  :: long_name_acc   ! long name for accumulator
    character(len=max_chars)  :: units           ! units of variable
    character(len=max_chars)  :: units_acc       ! accumulator units
    character(len=max_chars)  :: fname           ! full name of history file
    character(len=max_chars)  :: locrest(max_tapes) ! local history restart file names

    character(len=max_namlen),allocatable :: tname(:)
    character(len=max_chars), allocatable :: tunits(:),tlongname(:)
    character(len=hist_dim_name_length), allocatable :: tmpstr(:,:)
    character(len=1), allocatable :: tavgflag(:)
    integer :: start(2)

    character(len=1)   :: hnum                   ! history file index
    character(len=hist_dim_name_length)   :: type1d                 ! clm pointer 1d type
    character(len=hist_dim_name_length)   :: type1d_out             ! history buffer 1d type
    character(len=hist_dim_name_length)   :: type2d                 ! history buffer 2d type
    character(len=32)  :: dim1name               ! temporary
    character(len=32)  :: dim2name               ! temporary
    type(var_desc_t)   :: name_desc              ! variable descriptor for name
    type(var_desc_t)   :: longname_desc          ! variable descriptor for long_name
    type(var_desc_t)   :: units_desc             ! variable descriptor for units
    type(var_desc_t)   :: type1d_desc            ! variable descriptor for type1d
    type(var_desc_t)   :: type1d_out_desc        ! variable descriptor for type1d_out
    type(var_desc_t)   :: type2d_desc            ! variable descriptor for type2d
    type(var_desc_t)   :: avgflag_desc           ! variable descriptor for avgflag
    type(var_desc_t)   :: p2c_scale_type_desc    ! variable descriptor for p2c_scale_type
    type(var_desc_t)   :: c2l_scale_type_desc    ! variable descriptor for c2l_scale_type
    type(var_desc_t)   :: l2g_scale_type_desc    ! variable descriptor for l2g_scale_type
    type(var_desc_t)   :: t2g_scale_type_desc    ! variable descriptor for t2g_scale_type
    integer :: status                            ! error status
    integer :: dimid                             ! dimension ID
    integer :: k                                 ! 1d index
    integer :: ntapes_onfile                     ! number of history tapes on the restart file
    integer :: nflds_onfile                      ! number of history fields on the restart file
    integer :: t                                 ! tape index
    integer :: f                                 ! field index
    integer :: varid                             ! variable id
    integer, allocatable :: itemp2d(:,:)         ! 2D temporary
    real(r8), pointer :: hbuf(:,:)               ! history buffer
    real(r8), pointer :: hbuf1d(:)               ! 1d history buffer
    integer , pointer :: nacs(:,:)               ! accumulation counter
    integer , pointer :: nacs1d(:)               ! 1d accumulation counter
    integer           :: ier                     ! error code
    type(Var_desc_t)  :: vardesc                 ! netCDF variable description
    character(len=*),parameter :: subname = 'hist_restart_ncd'
!------------------------------------------------------------------------

    call get_proc_global(ng=numg, nt=numt, nl=numl, nc=numc, np=nump)

    ! If branch run, initialize file times and return

    if (flag == 'read') then
       if (nsrest == nsrBranch) then
          do t = 1,ntapes
             tape(t)%ntimes = 0
          end do
          return
       end if
       ! If startup run just return
       if (nsrest == nsrStartup) then
          RETURN
       end if
    endif

    ! Read history file data only for restart run (not for branch run)

    !
    ! First when writing out and in define mode, create files and define all variables
    !
    !================================================
    if (flag == 'define') then
    !================================================

       if (.not. present(rdate)) then
          call endrun(msg=' variable rdate must be present for writing restart files'//&
               errMsg(__FILE__, __LINE__))
       end if

       !
       ! On master restart file add ntapes/max_chars dimension
       ! and then add the history and history restart filenames
       !
       call ncd_defdim( ncid, 'ntapes'       , ntapes      , dimid)
       call ncd_defdim( ncid, 'max_chars'    , max_chars   , dimid)

       call ncd_defvar(ncid=ncid, varname='locfnh', xtype=ncd_char, &
            long_name="History filename",     &
            comment="This variable NOT needed for startup or branch simulations", &
            dim1name='max_chars', dim2name="ntapes" )
       ier = PIO_inq_varid(ncid, 'locfnh', vardesc)
       ier = PIO_put_att(ncid, vardesc%varid, 'interpinic_flag', iflag_skip)

       call ncd_defvar(ncid=ncid, varname='locfnhr', xtype=ncd_char, &
            long_name="Restart history filename",     &
            comment="This variable NOT needed for startup or branch simulations", &
            dim1name='max_chars', dim2name="ntapes" )
       ier = PIO_inq_varid(ncid, 'locfnhr', vardesc)
       ier = PIO_put_att(ncid, vardesc%varid, 'interpinic_flag', iflag_skip)

       ! max_nflds is the maximum number of fields on any tape
       ! max_flds is the maximum number possible number of fields 

       max_nflds = max_nFields()

       ! Loop over tapes - write out namelist information to each restart-history tape
       ! only read/write accumulators and counters if needed

       do t = 1,ntapes

          ! Create the restart history filename and open it
          write(hnum,'(i1.1)') t-1
          locfnhr(t) = "./" // trim(caseid) //".elm"// trim(inst_suffix) &
                        // ".rh" // hnum //"."// trim(rdate) //".nc"

          call htape_create( t, histrest=.true. )

          ! Add read/write accumultators and counters if needed
          if (.not. tape(t)%is_endhist) then
             do f = 1,tape(t)%nflds
                name           =  tape(t)%hlist(f)%field%name
                long_name      =  tape(t)%hlist(f)%field%long_name
                units          =  tape(t)%hlist(f)%field%units
                name_acc       =  trim(name) // "_acc"
                units_acc      =  "unitless positive integer"
                long_name_acc  =  trim(long_name) // " accumulator number of samples"
                type1d_out     =  tape(t)%hlist(f)%field%type1d_out
                type2d         =  tape(t)%hlist(f)%field%type2d
                num2d          =  tape(t)%hlist(f)%field%num2d
                nacs           => tape(t)%hlist(f)%nacs
                hbuf           => tape(t)%hlist(f)%hbuf
               
                if (type1d_out == grlnd) then
                   if (ldomain%isgrid2d) then
                      dim1name = 'lon'      ; dim2name = 'lat'
                   else
                      dim1name = trim(grlnd); dim2name = 'undefined'
                   end if
                else
                   dim1name = type1d_out ; dim2name = 'undefined'
                endif
                   
                if (dim2name == 'undefined') then
                   if (num2d == 1) then
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name), xtype=ncd_double, & 
                           dim1name=dim1name, &
                           long_name=trim(long_name), units=trim(units))
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name_acc), xtype=ncd_int,  &
                           dim1name=dim1name, &
                           long_name=trim(long_name_acc), units=trim(units_acc))
                   else
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name), xtype=ncd_double, &
                           dim1name=dim1name, dim2name=type2d, &
                           long_name=trim(long_name), units=trim(units))
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name_acc), xtype=ncd_int,  &
                           dim1name=dim1name, dim2name=type2d, &
                           long_name=trim(long_name_acc), units=trim(units_acc))
                   end if
                else
                   if (num2d == 1) then
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name), xtype=ncd_double, &
                           dim1name=dim1name, dim2name=dim2name, &
                           long_name=trim(long_name), units=trim(units))
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name_acc), xtype=ncd_int,  &
                           dim1name=dim1name, dim2name=dim2name, &
                           long_name=trim(long_name_acc), units=trim(units_acc))
                   else
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name), xtype=ncd_double, &
                           dim1name=dim1name, dim2name=dim2name, dim3name=type2d, &
                           long_name=trim(long_name), units=trim(units))
                      call ncd_defvar(ncid=ncid_hist(t), varname=trim(name_acc), xtype=ncd_int,  &
                           dim1name=dim1name, dim2name=dim2name, dim3name=type2d, &
                           long_name=trim(long_name_acc), units=trim(units_acc))
                   end if
                endif
             end do
          endif

          !
          ! Add namelist information to each restart history tape
          !
          call ncd_defdim( ncid_hist(t), 'fname_lenp2'  , max_namlen+2, dimid)
          call ncd_defdim( ncid_hist(t), 'fname_len'    , max_namlen  , dimid)
          call ncd_defdim( ncid_hist(t), 'len1'         , 1           , dimid)
          call ncd_defdim( ncid_hist(t), 'scalar'       , 1           , dimid)
          call ncd_defdim( ncid_hist(t), 'max_chars'    , max_chars   , dimid)
          call ncd_defdim( ncid_hist(t), 'max_nflds'    , max_nflds   ,  dimid)   
          call ncd_defdim( ncid_hist(t), 'max_flds'     , max_flds    , dimid)   
       
          call ncd_defvar(ncid=ncid_hist(t), varname='nhtfrq', xtype=ncd_int, &
               long_name="Frequency of history writes",               &
               comment="Namelist item", &
               units="absolute value of negative is in hours, 0=monthly, positive is time-steps",     &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='mfilt', xtype=ncd_int, &
               long_name="Number of history time samples on a file", units="unitless",     &
               comment="Namelist item", &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='ncprec', xtype=ncd_int, &
               long_name="Flag for data precision", flag_values=(/1,2/), &
               comment="Namelist item", &
               nvalid_range=(/1,2/), &
               flag_meanings=(/"single-precision", "double-precision"/), &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='dov2xy', xtype=ncd_log, &
               long_name="Output on 2D grid format (TRUE) or vector format (FALSE)", &
               comment="Namelist item", &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='fincl', xtype=ncd_char, &
               comment="Namelist item", &
               long_name="Fieldnames to include", &
               dim1name='fname_lenp2', dim2name='max_flds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='fexcl', xtype=ncd_char, &
               comment="Namelist item", &
               long_name="Fieldnames to exclude",  &
               dim1name='fname_lenp2', dim2name='max_flds' )

          call ncd_defvar(ncid=ncid_hist(t), varname='nflds', xtype=ncd_int, &
               long_name="Number of fields on file", units="unitless",        &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='ntimes', xtype=ncd_int, &
               long_name="Number of time steps on file", units="time-step",     &
               dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='is_endhist', xtype=ncd_log, &
               long_name="End of history file", dim1name='scalar')
          call ncd_defvar(ncid=ncid_hist(t), varname='begtime', xtype=ncd_double, &
               long_name="Beginning time", units="time units",     &
               dim1name='scalar')
   
          call ncd_defvar(ncid=ncid_hist(t), varname='num2d', xtype=ncd_int, &
               long_name="Size of second dimension", units="unitless",     &
               dim1name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='hpindex', xtype=ncd_int, &
               long_name="History pointer index", units="unitless",     &
               dim1name='max_nflds' )

          call ncd_defvar(ncid=ncid_hist(t), varname='avgflag', xtype=ncd_char, &
               long_name="Averaging flag", &
               units="A=Average, X=Maximum, M=Minimum, I=Instantaneous", &
               dim1name='len1', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='name', xtype=ncd_char, &
               long_name="Fieldnames",  &
               dim1name='fname_len', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='long_name', xtype=ncd_char, &
               long_name="Long descriptive names for fields", &
               dim1name='max_chars', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='units', xtype=ncd_char, &
               long_name="Units for each history field output", &
               dim1name='max_chars', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='type1d', xtype=ncd_char, &
               long_name="1st dimension type", &
               dim1name='string_length', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='type1d_out', xtype=ncd_char, &
               long_name="1st output dimension type", &
               dim1name='string_length', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='type2d', xtype=ncd_char, &
               long_name="2nd dimension type", &
               dim1name='string_length', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='p2c_scale_type', xtype=ncd_char, &
               long_name="PFT to column scale type", &
               dim1name='string_length', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='c2l_scale_type', xtype=ncd_char, &
               long_name="column to landunit scale type", &
               dim1name='string_length', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='l2g_scale_type', xtype=ncd_char, &
               long_name="landunit to gridpoint scale type", &
               dim1name='string_length', dim2name='max_nflds' )
          call ncd_defvar(ncid=ncid_hist(t), varname='t2g_scale_type', xtype=ncd_char, &
               long_name="topounit to gridpoint scale type", &
               dim1name='string_length', dim2name='max_nflds' )

          call ncd_enddef(ncid_hist(t))

       end do   ! end of ntapes loop   

       RETURN

    !
    ! First write out namelist information to each restart history file
    !
    !================================================
    else if (flag == 'write') then
    !================================================

       ! Add history filenames to master restart file
       do t = 1,ntapes
          call ncd_io('locfnh',  locfnh(t),  'write', ncid, nt=t)
          call ncd_io('locfnhr', locfnhr(t), 'write', ncid, nt=t)
       end do
       
       fincl(:,1) = hist_fincl1(:)
       fincl(:,2) = hist_fincl2(:)
       fincl(:,3) = hist_fincl3(:)
       fincl(:,4) = hist_fincl4(:)
       fincl(:,5) = hist_fincl5(:)
       fincl(:,6) = hist_fincl6(:)

       fexcl(:,1) = hist_fexcl1(:)
       fexcl(:,2) = hist_fexcl2(:)
       fexcl(:,3) = hist_fexcl3(:)
       fexcl(:,4) = hist_fexcl4(:)
       fexcl(:,5) = hist_fexcl5(:)
       fexcl(:,6) = hist_fexcl6(:)

       max_nflds = max_nFields()

       start(1)=1

       allocate(itemp2d(max_nflds,ntapes))

       !
       ! Add history namelist data to each history restart tape
       !
       do t = 1,ntapes
          call ncd_io(varname='fincl', data=fincl(:,t), ncid=ncid_hist(t), flag='write')

          call ncd_io(varname='fexcl', data=fexcl(:,t), ncid=ncid_hist(t), flag='write')

          call ncd_io(varname='is_endhist', data=tape(t)%is_endhist, ncid=ncid_hist(t), flag='write')

          call ncd_io(varname='dov2xy', data=tape(t)%dov2xy, ncid=ncid_hist(t), flag='write')

          itemp2d(:,:) = 0
          do f=1,tape(t)%nflds
             itemp2d(f,t) = tape(t)%hlist(f)%field%num2d
          end do
          call ncd_io(varname='num2d', data=itemp2d(:,t), ncid=ncid_hist(t), flag='write')

          itemp2d(:,:) = 0
          do f=1,tape(t)%nflds
             itemp2d(f,t) = tape(t)%hlist(f)%field%hpindex
          end do
          call ncd_io(varname='hpindex', data=itemp2d(:,t), ncid=ncid_hist(t), flag='write')

          call ncd_io('nflds',        tape(t)%nflds,   'write', ncid_hist(t) )
          call ncd_io('ntimes',       tape(t)%ntimes,  'write', ncid_hist(t) )
          call ncd_io('nhtfrq',  tape(t)%nhtfrq,  'write', ncid_hist(t) )
          call ncd_io('mfilt',   tape(t)%mfilt,   'write', ncid_hist(t) )
          call ncd_io('ncprec',  tape(t)%ncprec,  'write', ncid_hist(t) )
          call ncd_io('begtime',      tape(t)%begtime, 'write', ncid_hist(t) )
          allocate(tmpstr(tape(t)%nflds,7 ),tname(tape(t)%nflds), &
                   tavgflag(tape(t)%nflds),tunits(tape(t)%nflds),tlongname(tape(t)%nflds))
          do f=1,tape(t)%nflds
             tname(f)  = tape(t)%hlist(f)%field%name
             tunits(f) = tape(t)%hlist(f)%field%units
             tlongname(f) = tape(t)%hlist(f)%field%long_name
             tmpstr(f,1) = tape(t)%hlist(f)%field%type1d
             tmpstr(f,2) = tape(t)%hlist(f)%field%type1d_out
             tmpstr(f,3) = tape(t)%hlist(f)%field%type2d
             tavgflag(f) = tape(t)%hlist(f)%avgflag
             tmpstr(f,4) = tape(t)%hlist(f)%field%p2c_scale_type
             tmpstr(f,5) = tape(t)%hlist(f)%field%c2l_scale_type
             tmpstr(f,6) = tape(t)%hlist(f)%field%l2g_scale_type
             tmpstr(f,7) = tape(t)%hlist(f)%field%t2g_scale_type
          end do
          call ncd_io( 'name', tname, 'write',ncid_hist(t))
          call ncd_io('long_name', tlongname, 'write', ncid_hist(t))
          call ncd_io('units', tunits, 'write',ncid_hist(t))
          call ncd_io('type1d', tmpstr(:,1), 'write', ncid_hist(t))
          call ncd_io('type1d_out', tmpstr(:,2), 'write', ncid_hist(t))
          call ncd_io('type2d', tmpstr(:,3), 'write', ncid_hist(t))
          call ncd_io('avgflag',tavgflag , 'write', ncid_hist(t))
          call ncd_io('p2c_scale_type', tmpstr(:,4), 'write', ncid_hist(t))
          call ncd_io('c2l_scale_type', tmpstr(:,5), 'write', ncid_hist(t))
          call ncd_io('l2g_scale_type', tmpstr(:,6), 'write', ncid_hist(t))
          call ncd_io('t2g_scale_type', tmpstr(:,7), 'write', ncid_hist(t))
          deallocate(tname,tlongname,tunits,tmpstr,tavgflag)
       enddo       
       deallocate(itemp2d)

    !
    ! Read in namelist information
    !
    !================================================
    else if (flag == 'read') then
    !================================================

       call ncd_inqdlen(ncid,dimid,ntapes_onfile, name='ntapes')
       if ( is_restart() .and. ntapes_onfile /= ntapes )then
          write(iulog,*) 'ntapes = ', ntapes, ' ntapes_onfile = ', ntapes_onfile
          call endrun(msg=' ERROR: number of ntapes different than on restart file!,'// &
               ' you can NOT change history options on restart!' //&
               errMsg(__FILE__, __LINE__))
       end if
       if ( is_restart() .and. ntapes > 0 )then
          call ncd_io('locfnh',  locfnh(1:ntapes),  'read', ncid )
          call ncd_io('locfnhr', locrest(1:ntapes), 'read', ncid )
          do t = 1,ntapes
             call strip_null(locrest(t))
             call strip_null(locfnh(t))
          end do
       end if

       ! Determine necessary indices - the following is needed if model decomposition is different on restart
  
       start(1)=1

       if ( is_restart() )then
          do t = 1,ntapes

             call getfil( locrest(t), locfnhr(t), 0 )
             call ncd_pio_openfile (ncid_hist(t), trim(locfnhr(t)), ncd_nowrite)

             if ( t == 1 )then

                call ncd_inqdlen(ncid_hist(1),dimid,max_nflds,name='max_nflds')
   
                allocate(itemp2d(max_nflds,ntapes))
             end if

             call ncd_inqvid(ncid_hist(t), 'name',           varid, name_desc)
             call ncd_inqvid(ncid_hist(t), 'long_name',      varid, longname_desc)
             call ncd_inqvid(ncid_hist(t), 'units',          varid, units_desc)
             call ncd_inqvid(ncid_hist(t), 'type1d',         varid, type1d_desc)
             call ncd_inqvid(ncid_hist(t), 'type1d_out',     varid, type1d_out_desc)
             call ncd_inqvid(ncid_hist(t), 'type2d',         varid, type2d_desc)
             call ncd_inqvid(ncid_hist(t), 'avgflag',        varid, avgflag_desc)
             call ncd_inqvid(ncid_hist(t), 'p2c_scale_type', varid, p2c_scale_type_desc)
             call ncd_inqvid(ncid_hist(t), 'c2l_scale_type', varid, c2l_scale_type_desc)
             call ncd_inqvid(ncid_hist(t), 'l2g_scale_type', varid, l2g_scale_type_desc)
             call ncd_inqvid(ncid_hist(t), 't2g_scale_type', varid, t2g_scale_type_desc)

             call ncd_io(varname='fincl', data=fincl(:,t), ncid=ncid_hist(t), flag='read')

             call ncd_io(varname='fexcl', data=fexcl(:,t), ncid=ncid_hist(t), flag='read')

             call ncd_io('nflds',   nflds_onfile, 'read', ncid_hist(t) )
             if ( nflds_onfile /= tape(t)%nflds )then
                write(iulog,*) 'nflds = ', tape(t)%nflds, ' nflds_onfile = ', nflds_onfile
                call endrun(msg=' ERROR: number of fields different than on restart file!,'// &
                     ' you can NOT change history options on restart!' //&
                     errMsg(__FILE__, __LINE__))
             end if
             call ncd_io('ntimes',  tape(t)%ntimes, 'read', ncid_hist(t) )
             call ncd_io('nhtfrq',  tape(t)%nhtfrq, 'read', ncid_hist(t) )
             call ncd_io('mfilt',   tape(t)%mfilt, 'read', ncid_hist(t) )
             call ncd_io('ncprec',  tape(t)%ncprec, 'read', ncid_hist(t) )
             call ncd_io('begtime', tape(t)%begtime, 'read', ncid_hist(t) )

             call ncd_io(varname='is_endhist', data=tape(t)%is_endhist, ncid=ncid_hist(t), flag='read')
             call ncd_io(varname='dov2xy', data=tape(t)%dov2xy, ncid=ncid_hist(t), flag='read')
             call ncd_io(varname='num2d', data=itemp2d(:,t), ncid=ncid_hist(t), flag='read')
             do f=1,tape(t)%nflds
                tape(t)%hlist(f)%field%num2d = itemp2d(f,t)
             end do

             call ncd_io(varname='hpindex', data=itemp2d(:,t), ncid=ncid_hist(t), flag='read')
             do f=1,tape(t)%nflds
                tape(t)%hlist(f)%field%hpindex = itemp2d(f,t)
             end do

             do f=1,tape(t)%nflds
                start(2) = f
                call ncd_io( name_desc,           tape(t)%hlist(f)%field%name,       &
                             'read', ncid_hist(t), start )
                call ncd_io( longname_desc,       tape(t)%hlist(f)%field%long_name,  &
                             'read', ncid_hist(t), start )
                call ncd_io( units_desc,          tape(t)%hlist(f)%field%units,      &
                             'read', ncid_hist(t), start )
                call ncd_io( type1d_desc,         tape(t)%hlist(f)%field%type1d,     &
                             'read', ncid_hist(t), start )
                call ncd_io( type1d_out_desc,     tape(t)%hlist(f)%field%type1d_out, &
                             'read', ncid_hist(t), start )
                call ncd_io( type2d_desc,         tape(t)%hlist(f)%field%type2d,     &
                             'read', ncid_hist(t), start )
                call ncd_io( avgflag_desc,        tape(t)%hlist(f)%avgflag,          &
                             'read', ncid_hist(t), start )
                call ncd_io( p2c_scale_type_desc, tape(t)%hlist(f)%field%p2c_scale_type,   &
                             'read', ncid_hist(t), start )
                call ncd_io( c2l_scale_type_desc, tape(t)%hlist(f)%field%c2l_scale_type,   &
                             'read', ncid_hist(t), start )
                call ncd_io( l2g_scale_type_desc, tape(t)%hlist(f)%field%l2g_scale_type,   &
                             'read', ncid_hist(t), start )
                call ncd_io( t2g_scale_type_desc, tape(t)%hlist(f)%field%t2g_scale_type,   &
                             'read', ncid_hist(t), start )
                call strip_null(tape(t)%hlist(f)%field%name)
                call strip_null(tape(t)%hlist(f)%field%long_name)
                call strip_null(tape(t)%hlist(f)%field%units)
                call strip_null(tape(t)%hlist(f)%field%type1d)
                call strip_null(tape(t)%hlist(f)%field%type1d_out)
                call strip_null(tape(t)%hlist(f)%field%type2d)
                call strip_null(tape(t)%hlist(f)%field%p2c_scale_type)
                call strip_null(tape(t)%hlist(f)%field%c2l_scale_type)
                call strip_null(tape(t)%hlist(f)%field%l2g_scale_type)
                call strip_null(tape(t)%hlist(f)%field%t2g_scale_type)
                call strip_null(tape(t)%hlist(f)%avgflag)

                type1d_out = trim(tape(t)%hlist(f)%field%type1d_out)
                select case (trim(type1d_out))
                case (grlnd)
                   num1d_out = numg
                   beg1d_out = bounds%begg
                   end1d_out = bounds%endg
                case (nameg)
                   num1d_out = numg
                   beg1d_out = bounds%begg
                   end1d_out = bounds%endg
                case (namet)
                   num1d_out = numt
                   beg1d_out = bounds%begt
                   end1d_out = bounds%endt
                case (namel)
                   num1d_out = numl
                   beg1d_out = bounds%begl
                   end1d_out = bounds%endl
                case (namec)
                   num1d_out = numc
                   beg1d_out = bounds%begc
                   end1d_out = bounds%endc
                case (namep)
                   num1d_out = nump
                   beg1d_out = bounds%begp
                   end1d_out = bounds%endp
                case default
                   write(iulog,*) trim(subname),' ERROR: read unknown 1d output type=',trim(type1d_out)
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                end select

                tape(t)%hlist(f)%field%num1d_out = num1d_out
                tape(t)%hlist(f)%field%beg1d_out = beg1d_out
                tape(t)%hlist(f)%field%end1d_out = end1d_out

                num2d  = tape(t)%hlist(f)%field%num2d
                allocate (tape(t)%hlist(f)%hbuf(beg1d_out:end1d_out,num2d), &
                          tape(t)%hlist(f)%nacs(beg1d_out:end1d_out,num2d), &
                          stat=status)
                if (status /= 0) then
                   write(iulog,*) trim(subname),' ERROR: allocation error for hbuf,nacs at t,f=',t,f
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                endif
                tape(t)%hlist(f)%hbuf(:,:) = 0._r8
                tape(t)%hlist(f)%nacs(:,:) = 0

                type1d = tape(t)%hlist(f)%field%type1d
                select case (type1d)
                case (grlnd)
                   num1d = numg
                   beg1d = bounds%begg
                   end1d = bounds%endg
                case (nameg)
                   num1d = numg
                   beg1d = bounds%begg
                   end1d = bounds%endg
                case (namet)
                   num1d = numt
                   beg1d = bounds%begt
                   end1d = bounds%endt
                case (namel)
                   num1d = numl
                   beg1d = bounds%begl
                   end1d = bounds%endl
                case (namec)
                   num1d = numc
                   beg1d = bounds%begc
                   end1d = bounds%endc
                case (namep)
                   num1d = nump
                   beg1d = bounds%begp
                   end1d = bounds%endp
                case default
                   write(iulog,*) trim(subname),' ERROR: read unknown 1d type=',type1d
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                end select

                tape(t)%hlist(f)%field%num1d = num1d
                tape(t)%hlist(f)%field%beg1d = beg1d
                tape(t)%hlist(f)%field%end1d = end1d

             end do   ! end of flds loop

             ! If history file is not full, open it

             if (tape(t)%ntimes /= 0) then
                call ncd_pio_openfile (nfid(t), trim(locfnh(t)), ncd_write)
             end if

          end do  ! end of tapes loop

          hist_fincl1(:) = fincl(:,1)
          hist_fincl2(:) = fincl(:,2)
          hist_fincl3(:) = fincl(:,3)
          hist_fincl4(:) = fincl(:,4)
          hist_fincl5(:) = fincl(:,5)
          hist_fincl6(:) = fincl(:,6)

          hist_fexcl1(:) = fexcl(:,1)
          hist_fexcl2(:) = fexcl(:,2)
          hist_fexcl3(:) = fexcl(:,3)
          hist_fexcl4(:) = fexcl(:,4)
          hist_fexcl5(:) = fexcl(:,5)
          hist_fexcl6(:) = fexcl(:,6)

       end if
       
       if ( allocated(itemp2d) ) deallocate(itemp2d)

    end if

    !======================================================================
    ! Read/write history file restart data.
    ! If the current history file(s) are not full, file(s) are opened
    ! so that subsequent time samples are added until the file is full.
    ! A new history file is used on a branch run.
    !======================================================================
    
    if (flag == 'write') then     

       do t = 1,ntapes
          if (.not. tape(t)%is_endhist) then

             do f = 1,tape(t)%nflds
                name       =  tape(t)%hlist(f)%field%name
                name_acc   =  trim(name) // "_acc"
                type1d_out =  tape(t)%hlist(f)%field%type1d_out
                type2d     =  tape(t)%hlist(f)%field%type2d
                num2d      =  tape(t)%hlist(f)%field%num2d
                beg1d_out  =  tape(t)%hlist(f)%field%beg1d_out
                end1d_out  =  tape(t)%hlist(f)%field%end1d_out
                nacs       => tape(t)%hlist(f)%nacs
                hbuf       => tape(t)%hlist(f)%hbuf

                if (num2d == 1) then
                   allocate(hbuf1d(beg1d_out:end1d_out), &
                            nacs1d(beg1d_out:end1d_out), stat=status)
                   if (status /= 0) then
                      write(iulog,*) trim(subname),' ERROR: allocation'
                      call endrun(msg=errMsg(__FILE__, __LINE__))
                   end if
                
                   hbuf1d(beg1d_out:end1d_out) = hbuf(beg1d_out:end1d_out,1)
                   nacs1d(beg1d_out:end1d_out) = nacs(beg1d_out:end1d_out,1)

                   call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name), &
                        dim1name=type1d_out, data=hbuf1d)
                   call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name_acc), &
                        dim1name=type1d_out, data=nacs1d)

                   deallocate(hbuf1d)
                   deallocate(nacs1d)
                else
                   call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name), &
                        dim1name=type1d_out, data=hbuf)
                   call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name_acc), &
                        dim1name=type1d_out, data=nacs)
                end if

             end do

          end if  ! end of is_endhist block

          call ncd_pio_closefile(ncid_hist(t))

       end do   ! end of ntapes loop   

    else if (flag == 'read') then 

       ! Read history restart information if history files are not full

       do t = 1,ntapes

          if (.not. tape(t)%is_endhist) then

             do f = 1,tape(t)%nflds
                name       =  tape(t)%hlist(f)%field%name
                name_acc   =  trim(name) // "_acc"
                type1d_out =  tape(t)%hlist(f)%field%type1d_out
                type2d     =  tape(t)%hlist(f)%field%type2d
                num2d      =  tape(t)%hlist(f)%field%num2d
                beg1d_out  =  tape(t)%hlist(f)%field%beg1d_out
                end1d_out  =  tape(t)%hlist(f)%field%end1d_out
                nacs       => tape(t)%hlist(f)%nacs
                hbuf       => tape(t)%hlist(f)%hbuf
                
                if (num2d == 1) then
                   allocate(hbuf1d(beg1d_out:end1d_out), &
                        nacs1d(beg1d_out:end1d_out), stat=status)
                   if (status /= 0) then
                      write(iulog,*) trim(subname),' ERROR: allocation'
                      call endrun(msg=errMsg(__FILE__, __LINE__))
                   end if
                   
                   call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name), &
                        dim1name=type1d_out, data=hbuf1d)
                   call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name_acc), &
                        dim1name=type1d_out, data=nacs1d)
                   
                   hbuf(beg1d_out:end1d_out,1) = hbuf1d(beg1d_out:end1d_out)
                   nacs(beg1d_out:end1d_out,1) = nacs1d(beg1d_out:end1d_out)
                   
                   deallocate(hbuf1d)
                   deallocate(nacs1d)
                else
                   call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name), &
                        dim1name=type1d_out, data=hbuf)
                   call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name_acc), &
                        dim1name=type1d_out, data=nacs)
                end if
             end do

          end if
             
          call ncd_pio_closefile(ncid_hist(t))
             
       end do
       
    end if
    
  end subroutine hist_restart_ncd

  !-----------------------------------------------------------------------
  integer function max_nFields()
    !
    ! !DESCRIPTION:
    ! Get the maximum number of fields on all tapes.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: t  ! index
    character(len=*),parameter :: subname = 'max_nFields'
    !-----------------------------------------------------------------------

    max_nFields = 0
    do t = 1,ntapes
       max_nFields = max(max_nFields, tape(t)%nflds)
    end do
    return
  end function max_nFields
  
  !-----------------------------------------------------------------------
  character(len=max_namlen) function getname (inname)
    !
    ! !DESCRIPTION:
    ! Retrieve name portion of inname. If an averaging flag separater character
    ! is present (:) in inname, lop it off.
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: inname
    !
    ! !LOCAL VARIABLES:
    integer :: length
    integer :: i
    character(len=*),parameter :: subname = 'getname'
    !-----------------------------------------------------------------------

     length = len (inname)

     if (length < max_namlen .or. length > max_namlen+2) then
        write(iulog,*) trim(subname),' ERROR: bad length=',length
        call endrun(msg=errMsg(__FILE__, __LINE__))
     end if

     getname = ' '
     do i = 1,max_namlen
        if (inname(i:i) == ':') exit
        getname(i:i) = inname(i:i)
     end do

   end function getname

   !-----------------------------------------------------------------------
   character(len=1) function getflag (inname)
     !
     ! !DESCRIPTION:
     ! Retrieve flag portion of inname. If an averaging flag separater character
     ! is present (:) in inname, return the character after it as the flag
     !
     ! !ARGUMENTS:
     character(len=*) inname   ! character string
     !
     ! !LOCAL VARIABLES:
     integer :: length         ! length of inname
     integer :: i              ! loop index
     character(len=*),parameter :: subname = 'getflag'
     !-----------------------------------------------------------------------

     length = len (inname)

     if (length < max_namlen .or. length > max_namlen+2) then
        write(iulog,*) trim(subname),' ERROR: bad length=',length
        call endrun(msg=errMsg(__FILE__, __LINE__))
     end if

     getflag = ' '
     do i = 1,length
        if (inname(i:i) == ':') then
           getflag = inname(i+1:i+1)
           exit
        end if
     end do

   end function getflag

   !-----------------------------------------------------------------------
   subroutine list_index (list, name, index)
     !
     ! !ARGUMENTS:
     character(len=*), intent(in) :: list(max_flds)  ! input list of names, possibly ":" delimited
     character(len=max_namlen), intent(in) :: name   ! name to be searched for
     integer, intent(out) :: index                   ! index of "name" in "list"
     !
     ! !LOCAL VARIABLES:
     !EOP
     character(len=max_namlen) :: listname           ! input name with ":" stripped off.
     integer f                                       ! field index
     character(len=*),parameter :: subname = 'list_index'
     !-----------------------------------------------------------------------

     ! Only list items

     index = 0
     do f=1,max_flds
        listname = getname (list(f))
        if (listname == ' ') exit
        if (listname == name) then
           index = f
           exit
        end if
     end do

   end subroutine list_index

   !-----------------------------------------------------------------------
   character(len=256) function set_hist_filename (hist_freq, hist_mfilt, hist_file)
     !
     ! !DESCRIPTION:
     ! Determine history dataset filenames.
     !
     ! !USES:
     use elm_varctl, only : caseid, inst_suffix
     use clm_time_manager, only : get_curr_date, get_prev_date
     !
     ! !ARGUMENTS:
     integer, intent(in)  :: hist_freq   !history file frequency
     integer, intent(in)  :: hist_mfilt  !history file number of time-samples
     integer, intent(in)  :: hist_file   !history file index
     !
     ! !LOCAL VARIABLES:
     !EOP
     character(len=256) :: cdate       !date char string
     character(len=  1) :: hist_index  !p,1 or 2 (currently)
     integer :: day                    !day (1 -> 31)
     integer :: mon                    !month (1 -> 12)
     integer :: yr                     !year (0 -> ...)
     integer :: sec                    !seconds into current day
     character(len=*),parameter :: subname = 'set_hist_filename'
     !-----------------------------------------------------------------------

   if (hist_freq == 0 .and. hist_mfilt == 1) then   !monthly
      call get_prev_date (yr, mon, day, sec)
      write(cdate,'(i4.4,"-",i2.2)') yr,mon
   else                        !other
      call get_curr_date (yr, mon, day, sec)
      write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
   endif
   write(hist_index,'(i1.1)') hist_file - 1
   set_hist_filename = "./"//trim(caseid)//".elm"//trim(inst_suffix)//&
                       ".h"//hist_index//"."//trim(cdate)//".nc"

  end function set_hist_filename

  !-----------------------------------------------------------------------
  subroutine hist_addfld1d (fname, units, avgflag, long_name, type1d_out, &
                        ptr_gcell, ptr_topo, ptr_lunit, ptr_col, ptr_patch, ptr_lnd, &
                        ptr_atm, p2c_scale_type, c2l_scale_type, &
                        l2g_scale_type, t2g_scale_type, set_lake, set_nolake, set_urb, set_nourb, &
                        set_noglcmec, set_spec, default)
    !
    ! !DESCRIPTION:
    ! Initialize a single level history field. The pointer, ptrhist,
    ! is a pointer to the data type array that the history buffer will use.
    ! The value of type1d passed to masterlist\_add\_fld determines which of the
    ! 1d type of the output and the beginning and ending indices the history
    ! buffer field). Default history contents for given field on all tapes
    ! are set by calling [masterlist\_make\_active] for the appropriate tape.
    ! After the masterlist is built, routine [htapes\_build] is called for an
    ! initial or branch run to initialize the actual history tapes.
    !
    ! !ARGUMENTS:
    character(len=*), intent(in)           :: fname          ! field name
    character(len=*), intent(in)           :: units          ! units of field
    character(len=1), intent(in)           :: avgflag        ! time averaging flag
    character(len=*), intent(in)           :: long_name      ! long name of field
    character(len=*), optional, intent(in) :: type1d_out     ! output type (from data type)
    real(r8)        , optional, pointer    :: ptr_gcell(:)   ! pointer to gridcell array
    real(r8)        , optional, pointer    :: ptr_topo(:)    ! pointer to topounit array
    real(r8)        , optional, pointer    :: ptr_lunit(:)   ! pointer to landunit array
    real(r8)        , optional, pointer    :: ptr_col(:)     ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_patch(:)     ! pointer to pft array
    real(r8)        , optional, pointer    :: ptr_lnd(:)     ! pointer to lnd array
    real(r8)        , optional, pointer    :: ptr_atm(:)     ! pointer to atm array
    real(r8)        , optional, intent(in) :: set_lake       ! value to set lakes to
    real(r8)        , optional, intent(in) :: set_nolake     ! value to set non-lakes to
    real(r8)        , optional, intent(in) :: set_urb        ! value to set urban to
    real(r8)        , optional, intent(in) :: set_nourb      ! value to set non-urban to
    real(r8)        , optional, intent(in) :: set_noglcmec   ! value to set non-glacier_mec to
    real(r8)        , optional, intent(in) :: set_spec       ! value to set special to
    character(len=*), optional, intent(in) :: p2c_scale_type ! scale type for subgrid averaging of pfts to column
    character(len=*), optional, intent(in) :: c2l_scale_type ! scale type for subgrid averaging of columns to landunits
    character(len=*), optional, intent(in) :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
    character(len=*), optional, intent(in) :: t2g_scale_type ! scale type for subgrid averaging of topounits to gridcells
    character(len=*), optional, intent(in) :: default        ! if set to 'inactive, field will not appear on primary tape
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l,t,g               ! indices
    integer :: hpindex                 ! history buffer pointer index
    character(len=hist_dim_name_length) :: l_type1d       ! 1d data type
    character(len=hist_dim_name_length) :: l_type1d_out   ! 1d output type
    character(len=8) :: scale_type_p2c ! scale type for subgrid averaging of pfts to column
    character(len=8) :: scale_type_c2l ! scale type for subgrid averaging of columns to landunits
    character(len=8) :: scale_type_l2g ! scale type for subgrid averaging of landunits to gridcells
    character(len=8) :: scale_type_t2g ! scale type for subgrid averaging of topounits to gridcells
    type(bounds_type):: bounds         ! boudns 
    character(len=16):: l_default      ! local version of 'default'
    character(len=*),parameter :: subname = 'hist_addfld1d'
!------------------------------------------------------------------------

    ! Determine processor bounds

    call get_proc_bounds(bounds)

    ! History buffer pointer

    hpindex = pointer_index()

    if (present(ptr_lnd)) then
       l_type1d = grlnd
       l_type1d_out = grlnd
       elmptr_rs(hpindex)%ptr => ptr_lnd

    else if (present(ptr_gcell)) then
       l_type1d = nameg
       l_type1d_out = nameg
       elmptr_rs(hpindex)%ptr => ptr_gcell

    else if (present(ptr_topo)) then
       l_type1d = namet
       l_type1d_out = namet
       elmptr_rs(hpindex)%ptr => ptr_topo

    else if (present(ptr_lunit)) then
       l_type1d = namel
       l_type1d_out = namel
       elmptr_rs(hpindex)%ptr => ptr_lunit
       if (present(set_lake)) then
          do l = bounds%begl,bounds%endl
             if (lun_pp%lakpoi(l)) ptr_lunit(l) = set_lake
          end do
       end if
       if (present(set_nolake)) then
          do l = bounds%begl,bounds%endl
             if (.not.(lun_pp%lakpoi(l))) ptr_lunit(l) = set_nolake
          end do
       end if
       if (present(set_urb)) then
          do l = bounds%begl,bounds%endl
             if (lun_pp%urbpoi(l)) ptr_lunit(l) = set_urb
          end do
       end if
       if (present(set_nourb)) then
          do l = bounds%begl,bounds%endl
             if (.not.(lun_pp%urbpoi(l))) ptr_lunit(l) = set_nourb
          end do
       end if
       if (present(set_spec)) then
          do l = bounds%begl,bounds%endl
             if (lun_pp%ifspecial(l)) ptr_lunit(l) = set_spec
          end do
       end if

    else if (present(ptr_col)) then
       l_type1d = namec
       l_type1d_out = namec
       elmptr_rs(hpindex)%ptr => ptr_col
       if (present(set_lake)) then
          do c = bounds%begc,bounds%endc
             l =col_pp%landunit(c)
             if (lun_pp%lakpoi(l)) ptr_col(c) = set_lake
          end do
       end if
       if (present(set_nolake)) then
          do c = bounds%begc,bounds%endc
             l =col_pp%landunit(c)
             if (.not.(lun_pp%lakpoi(l))) ptr_col(c) = set_nolake
          end do
       end if
       if (present(set_urb)) then
          do c = bounds%begc,bounds%endc
             l =col_pp%landunit(c)
             if (lun_pp%urbpoi(l)) ptr_col(c) = set_urb
          end do
       end if
       if (present(set_nourb)) then
          do c = bounds%begc,bounds%endc
             l =col_pp%landunit(c)
             if (.not.(lun_pp%urbpoi(l))) ptr_col(c) = set_nourb
          end do
       end if
       if (present(set_spec)) then
          do c = bounds%begc,bounds%endc
             l =col_pp%landunit(c)
             if (lun_pp%ifspecial(l)) ptr_col(c) = set_spec
          end do
       end if
       if (present(set_noglcmec)) then
          do c = bounds%begc,bounds%endc
             l =col_pp%landunit(c)
             if (.not.(lun_pp%glcmecpoi(l))) ptr_col(c) = set_noglcmec
          end do
       endif

    else if (present(ptr_patch)) then
       l_type1d = namep
       l_type1d_out = namep
       elmptr_rs(hpindex)%ptr => ptr_patch
       if (present(set_lake)) then
          do p = bounds%begp,bounds%endp
             l =veg_pp%landunit(p)
             if (lun_pp%lakpoi(l)) ptr_patch(p) = set_lake
          end do
       end if
       if (present(set_nolake)) then
          do p = bounds%begp,bounds%endp
             l =veg_pp%landunit(p)
             if (.not.(lun_pp%lakpoi(l))) ptr_patch(p) = set_nolake
          end do
       end if
       if (present(set_urb)) then
          do p = bounds%begp,bounds%endp
             l =veg_pp%landunit(p)
             if (lun_pp%urbpoi(l)) ptr_patch(p) = set_urb
          end do
       end if
       if (present(set_nourb)) then
          do p = bounds%begp,bounds%endp
             l =veg_pp%landunit(p)
             if (.not.(lun_pp%urbpoi(l))) ptr_patch(p) = set_nourb
          end do
       end if
       if (present(set_spec)) then
          do p = bounds%begp,bounds%endp
             l =veg_pp%landunit(p)
             if (lun_pp%ifspecial(l)) ptr_patch(p) = set_spec
          end do
       end if
       if (present(set_noglcmec)) then
          do p = bounds%begp,bounds%endp
             l =veg_pp%landunit(p)
             if (.not.(lun_pp%glcmecpoi(l))) ptr_patch(p) = set_noglcmec
          end do
       end if
    else
       write(iulog,*) trim(subname),' ERROR: must specify a valid pointer index,', &
          ' choices are [ptr_atm, ptr_lnd, ptr_gcell, ptr_topo, ptr_lunit, ptr_col, ptr_patch] '
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end if

    ! Set scaling factor

    scale_type_p2c = 'unity'
    scale_type_c2l = 'unity'
    scale_type_l2g = 'unity'
    scale_type_t2g = 'unity'

    if (present(p2c_scale_type)) scale_type_p2c = p2c_scale_type
    if (present(c2l_scale_type)) scale_type_c2l = c2l_scale_type
    if (present(l2g_scale_type)) scale_type_l2g = l2g_scale_type
    if (present(t2g_scale_type)) scale_type_t2g = t2g_scale_type
    if (present(type1d_out)) l_type1d_out = type1d_out

    ! Add field to masterlist

    call masterlist_addfld (fname=trim(fname), type1d=l_type1d, type1d_out=l_type1d_out, &
         type2d='unset', numdims=1, num2d=1, &
         units=units, avgflag=avgflag, long_name=long_name, hpindex=hpindex, &
         p2c_scale_type=scale_type_p2c, c2l_scale_type=scale_type_c2l, l2g_scale_type=scale_type_l2g, &
         t2g_scale_type=scale_type_t2g)

    l_default = 'active'
    if (present(default)) then
       l_default = default
    end if
    if (trim(l_default) == 'inactive') then
       return
    else
       call masterlist_make_active (name=trim(fname), tape_index=1)
    end if

  end subroutine hist_addfld1d

  !-----------------------------------------------------------------------
  subroutine hist_addfld2d (fname, type2d, units, avgflag, long_name, type1d_out, &
                        ptr_gcell, ptr_topo, ptr_lunit, ptr_col, ptr_patch, ptr_lnd, ptr_atm, &
                        p2c_scale_type, c2l_scale_type, l2g_scale_type, t2g_scale_type, &
                        set_lake, set_nolake, set_urb, set_nourb, set_spec, &
                        no_snow_behavior, default)
    !
    ! !DESCRIPTION:
    ! Initialize a single level history field. The pointer, ptrhist,
    ! is a pointer to the data type array that the history buffer will use.
    ! The value of type1d passed to masterlist\_add\_fld determines which of the
    ! 1d type of the output and the beginning and ending indices the history
    ! buffer field). Default history contents for given field on all tapes
    ! are set by calling [masterlist\_make\_active] for the appropriatae tape.
    ! After the masterlist is built, routine [htapes\_build] is called for an
    ! initial or branch run to initialize the actual history tapes.
    !
    ! !USES:
    use elm_varpar      , only : nlevgrnd, nlevsno, nlevlak, numrad, nlevdecomp_full, nlevtrc_soil, nmonth, nvegwcs
    use elm_varpar      , only : natpft_size, cft_size, maxpatch_glcmec
    use landunit_varcon , only : max_lunit
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: fname                      ! field name
    character(len=*), intent(in) :: type2d                     ! 2d output type
    character(len=*), intent(in) :: units                      ! units of field
    character(len=1), intent(in) :: avgflag                    ! time averaging flag
    character(len=*), intent(in) :: long_name                  ! long name of field
    character(len=*), optional, intent(in) :: type1d_out       ! output type (from data type)
    real(r8)        , optional, pointer    :: ptr_atm(:,:)     ! pointer to atm array
    real(r8)        , optional, pointer    :: ptr_lnd(:,:)     ! pointer to lnd array
    real(r8)        , optional, pointer    :: ptr_gcell(:,:)   ! pointer to gridcell array
    real(r8)        , optional, pointer    :: ptr_topo(:,:)    ! pointer to topounit array
    real(r8)        , optional, pointer    :: ptr_lunit(:,:)   ! pointer to landunit array
    real(r8)        , optional, pointer    :: ptr_col(:,:)     ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_patch(:,:)     ! pointer to pft array
    real(r8)        , optional, intent(in) :: set_lake         ! value to set lakes to
    real(r8)        , optional, intent(in) :: set_nolake       ! value to set non-lakes to
    real(r8)        , optional, intent(in) :: set_urb          ! value to set urban to
    real(r8)        , optional, intent(in) :: set_nourb        ! value to set non-urban to
    real(r8)        , optional, intent(in) :: set_spec         ! value to set special to
    integer         , optional, intent(in) :: no_snow_behavior ! if a multi-layer snow field, behavior to use for absent snow layers (should be one of the public no_snow_* parameters defined above)
    character(len=*), optional, intent(in) :: p2c_scale_type   ! scale type for subgrid averaging of pfts to column
    character(len=*), optional, intent(in) :: c2l_scale_type   ! scale type for subgrid averaging of columns to landunits
    character(len=*), optional, intent(in) :: l2g_scale_type   ! scale type for subgrid averaging of landunits to gridcells
    character(len=*), optional, intent(in) :: t2g_scale_type   ! scale type for subgrid averaging of topounits to gridcells
    character(len=*), optional, intent(in) :: default          ! if set to 'inactive, field will not appear on primary tape
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l,t,g                 ! indices
    integer :: num2d                   ! size of second dimension (e.g. number of vertical levels)
    integer :: hpindex                 ! history buffer index
    character(len=hist_dim_name_length) :: l_type1d         ! 1d data type
    character(len=hist_dim_name_length) :: l_type1d_out     ! 1d output type
    character(len=8) :: scale_type_p2c ! scale type for subgrid averaging of pfts to column
    character(len=8) :: scale_type_c2l ! scale type for subgrid averaging of columns to landunits
    character(len=8) :: scale_type_l2g ! scale type for subgrid averaging of landunits to gridcells
    character(len=8) :: scale_type_t2g ! scale type for subgrid averaging of topounits to gridcells
    type(bounds_type):: bounds          
    character(len=16):: l_default      ! local version of 'default'
    character(len=*),parameter :: subname = 'hist_addfld2d'
!------------------------------------------------------------------------

    call get_proc_bounds(bounds)
    
    ! Error-check no_snow_behavior optional argument: It should be present if and only if
    ! type2d is 'levsno', and its value should be one of the public no_snow_* parameters
    ! defined above.
    if (present(no_snow_behavior)) then
       if (type2d /= 'levsno') then
          write(iulog,*) trim(subname), &
               ' ERROR: Only specify no_snow_behavior for fields with dimension levsno'
          call endrun()
       end if

       if (no_snow_behavior < no_snow_MIN .or. no_snow_behavior > no_snow_MAX) then
          write(iulog,*) trim(subname), &
               ' ERROR: Invalid value for no_snow_behavior: ', no_snow_behavior
          call endrun()
       end if

    else  ! no_snow_behavior is absent
       if (type2d == 'levsno') then
          write(iulog,*) trim(subname), &
               ' ERROR: must specify no_snow_behavior for fields with dimension levsno'
          call endrun()
       end if
    end if

    ! Determine second dimension size

    select case (type2d)
    case ('levgrnd')
       num2d = nlevgrnd
    case ('levlak')
       num2d = nlevlak
    case ('numrad')
       num2d = numrad
    case ('month')
       num2d = nmonth
    case ('levdcmp')
       num2d = nlevdecomp_full
    case ('levtrc')
       num2d = nlevtrc_soil       
    case('ltype')
       num2d = max_lunit
    case('natpft')
       num2d = natpft_size
    case ('fates_levelem')
       num2d = nelements_fates
    case ('fates_levelpft')
       num2d = nelements_fates*numpft_fates
    case ('fates_levelcwd')
       num2d = nelements_fates*ncwd_fates
    case ('fates_levelage')
       num2d = nelements_fates*nlevage_fates

    case('cft')
       if (cft_size > 0) then
          num2d = cft_size
       else
          write(iulog,*) trim(subname),' ERROR: 2d type =', trim(type2d), &
               ' only valid for cft_size > 0'
          call endrun()
       end if
    case ('glc_nec')
       if (maxpatch_glcmec > 0) then
          num2d = maxpatch_glcmec
       else
          write(iulog,*) trim(subname),' ERROR: 2d type =', trim(type2d), &
               ' only valid for maxpatch_glcmec > 0'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    case ('elevclas')
       if (maxpatch_glcmec > 0) then
          ! add one because indexing starts at 0 (elevclas, unlike glc_nec, includes the
          ! bare ground "elevation class")
          num2d = maxpatch_glcmec + 1
       else
          write(iulog,*) trim(subname),' ERROR: 2d type =', trim(type2d), &
               ' only valid for maxpatch_glcmec > 0'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    case ('levsno')
       num2d = nlevsno
    case ('nvegwcs')
        num2d = nvegwcs
    case ('fates_levscls')
       num2d = nlevsclass_fates
    case('fates_levcacls')
       num2d = nlevcoage
    case ('fates_levcapf')
       num2d = nlevcoage*numpft_fates
    case ('fates_levpft')
       num2d = numpft_fates
    case ('fates_levage')
       num2d = nlevage_fates
    case ('fates_levfuel')
       num2d = nfsc_fates
    case ('fates_levcwdsc')
       num2d = ncwd_fates
    case ('fates_levscpf')
       num2d = nlevsclass_fates*numpft_fates
    case ('fates_levscag')
       num2d = nlevsclass_fates*nlevage_fates
    case ('fates_levcan')
       num2d = nclmax_fates
    case ('fates_levcnlf')
       num2d = nlevleaf_fates * nclmax_fates
    case ('fates_levcnlfpf')
       num2d = nlevleaf_fates * nclmax_fates * numpft_fates
    case ('fates_levheight')
       num2d = nlevheight_fates
    case ('fates_levscagpf')
       num2d = nlevsclass_fates*nlevage_fates*numpft_fates
    case ('fates_levagepft')
       num2d = nlevage_fates*numpft_fates
    case default
       write(iulog,*) trim(subname),' ERROR: unsupported 2d type ',type2d, &
          ' currently supported types for multi level fields are: ', &
          '[levgrnd,levlak,numrad,nmonthlevdcmp,levtrc,ltype,natpft,cft,glc_nec,elevclas,levsno]'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    ! History buffer pointer

    hpindex = pointer_index()

    if (present(ptr_lnd)) then
       l_type1d = grlnd
       l_type1d_out = grlnd
       elmptr_ra(hpindex)%ptr => ptr_lnd

    else if (present(ptr_gcell)) then
       l_type1d = nameg
       l_type1d_out = nameg
       elmptr_ra(hpindex)%ptr => ptr_gcell

    else if (present(ptr_topo)) then
       l_type1d = namet
       l_type1d_out = namet
       elmptr_ra(hpindex)%ptr => ptr_topo

    else if (present(ptr_lunit)) then
       l_type1d = namel
       l_type1d_out = namel
       elmptr_ra(hpindex)%ptr => ptr_lunit
       if (present(set_lake)) then
          do l = bounds%begl,bounds%endl
             if (lun_pp%lakpoi(l)) ptr_lunit(l,:) = set_lake
          end do
       end if
       if (present(set_nolake)) then
          do l = bounds%begl,bounds%endl
             if (.not.(lun_pp%lakpoi(l))) ptr_lunit(l,:) = set_nolake
          end do
       end if
       if (present(set_urb)) then
          do l = bounds%begl,bounds%endl
             if (lun_pp%urbpoi(l)) ptr_lunit(l,:) = set_urb
          end do
       end if
       if (present(set_nourb)) then
          do l = bounds%begl,bounds%endl
             if (.not.(lun_pp%urbpoi(l))) ptr_lunit(l,:) = set_nourb
          end do
       end if
       if (present(set_spec)) then
          do l = bounds%begl,bounds%endl
             if (lun_pp%ifspecial(l)) ptr_lunit(l,:) = set_spec
          end do
       end if

    else if (present(ptr_col)) then
       l_type1d = namec
       l_type1d_out = namec
       elmptr_ra(hpindex)%ptr => ptr_col
       if (present(set_lake)) then
          do c = bounds%begc,bounds%endc
             l =col_pp%landunit(c)
             if (lun_pp%lakpoi(l)) ptr_col(c,:) = set_lake
          end do
       end if
       if (present(set_nolake)) then
          do c = bounds%begc,bounds%endc
             l =col_pp%landunit(c)
             if (.not.(lun_pp%lakpoi(l))) ptr_col(c,:) = set_nolake
          end do
       end if
       if (present(set_urb)) then
          do c = bounds%begc,bounds%endc
             l =col_pp%landunit(c)
             if (lun_pp%urbpoi(l)) ptr_col(c,:) = set_urb
          end do
       end if
       if (present(set_nourb)) then
          do c = bounds%begc,bounds%endc
             l =col_pp%landunit(c)
             if (.not.(lun_pp%urbpoi(l))) ptr_col(c,:) = set_nourb
          end do
       end if
       if (present(set_spec)) then
          do c = bounds%begc,bounds%endc
             l =col_pp%landunit(c)
             if (lun_pp%ifspecial(l)) ptr_col(c,:) = set_spec
          end do
       end if

    else if (present(ptr_patch)) then
       l_type1d = namep
       l_type1d_out = namep
       elmptr_ra(hpindex)%ptr => ptr_patch
       if (present(set_lake)) then
          do p = bounds%begp,bounds%endp
             l =veg_pp%landunit(p)
             if (lun_pp%lakpoi(l)) ptr_patch(p,:) = set_lake
          end do
       end if
       if (present(set_nolake)) then
          do p = bounds%begp,bounds%endp
             l =veg_pp%landunit(p)
             if (.not.(lun_pp%lakpoi(l))) ptr_patch(p,:) = set_nolake
          end do
       end if
       if (present(set_urb)) then
          do p = bounds%begp,bounds%endp
             l =veg_pp%landunit(p)
             if (lun_pp%urbpoi(l)) ptr_patch(p,:) = set_urb
          end do
       end if
       if (present(set_nourb)) then
          do p = bounds%begp,bounds%endp
             l =veg_pp%landunit(p)
             if (.not.(lun_pp%urbpoi(l))) ptr_patch(p,:) = set_nourb
          end do
       end if
       if (present(set_spec)) then
          do p = bounds%begp,bounds%endp
             l =veg_pp%landunit(p)
             if (lun_pp%ifspecial(l)) ptr_patch(p,:) = set_spec
          end do
       end if

    else
       write(iulog,*) trim(subname),' ERROR: must specify a valid pointer index,', &
          ' choices are ptr_atm, ptr_lnd, ptr_gcell, ptr_lunit, ptr_col, ptr_patch'
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end if

    ! Set scaling factor

    scale_type_p2c = 'unity'
    scale_type_c2l = 'unity'
    scale_type_l2g = 'unity'
    scale_type_t2g = 'unity'

    if (present(p2c_scale_type)) scale_type_p2c = p2c_scale_type
    if (present(c2l_scale_type)) scale_type_c2l = c2l_scale_type
    if (present(l2g_scale_type)) scale_type_l2g = l2g_scale_type
    if (present(t2g_scale_type)) scale_type_t2g = t2g_scale_type
    if (present(type1d_out)) l_type1d_out = type1d_out

    ! Add field to masterlist

    call masterlist_addfld (fname=trim(fname), type1d=l_type1d, type1d_out=l_type1d_out, &
         type2d=type2d, numdims=2, num2d=num2d, &
         units=units, avgflag=avgflag, long_name=long_name, hpindex=hpindex, &
         p2c_scale_type=scale_type_p2c, c2l_scale_type=scale_type_c2l, l2g_scale_type=scale_type_l2g, &
         t2g_scale_type=scale_type_t2g, no_snow_behavior=no_snow_behavior)

    l_default = 'active'
    if (present(default)) then
       l_default = default
    end if
    if (trim(l_default) == 'inactive') then
       return
    else
       call masterlist_make_active (name=trim(fname), tape_index=1)
    end if

  end subroutine hist_addfld2d

  !-----------------------------------------------------------------------
  subroutine hist_addfld_decomp (fname, type2d, units, avgflag, long_name, ptr_col, ptr_patch, default)

    !
    ! !USES:
    use elm_varpar  , only : nlevdecomp_full, crop_prog
    use elm_varctl  , only : iulog
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: fname                    ! field name
    character(len=*), intent(in) :: type2d                   ! 2d output type
    character(len=*), intent(in) :: units                    ! units of field
    character(len=1), intent(in) :: avgflag                  ! time averaging flag
    character(len=*), intent(in) :: long_name                ! long name of field
    real(r8)        , optional, pointer    :: ptr_col(:,:)   ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_patch(:,:)   ! pointer to pft array
    character(len=*), optional, intent(in) :: default        ! if set to 'inactive, field will not appear on primary tape
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer  :: ptr_1d(:)
    !-----------------------------------------------------------------------

    if (present(ptr_col)) then

       ! column-level data
       if (present(default)) then
          if ( nlevdecomp_full > 1 ) then
             call hist_addfld2d (fname=trim(fname), units=units, type2d=type2d, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_col=ptr_col, default=default)
          else
             ptr_1d => ptr_col(:,1)
             call hist_addfld1d (fname=trim(fname), units=units, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_col=ptr_1d, default=default)
          endif
       else
          if ( nlevdecomp_full > 1 ) then
             call hist_addfld2d (fname=trim(fname), units=units, type2d=type2d, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_col=ptr_col)
          else
             ptr_1d => ptr_col(:,1)
             call hist_addfld1d (fname=trim(fname), units=units, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_col=ptr_1d)
          endif
       endif

    else if (present(ptr_patch)) then

       ! pft-level data
       if (present(default)) then
          if ( nlevdecomp_full > 1 ) then
             call hist_addfld2d (fname=trim(fname), units=units, type2d=type2d, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_patch=ptr_patch, default=default)
          else
             ptr_1d => ptr_patch(:,1)
             call hist_addfld1d (fname=trim(fname), units=units, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_patch=ptr_1d, default=default)
          endif
       else
          if ( nlevdecomp_full > 1 ) then
             call hist_addfld2d (fname=trim(fname), units=units, type2d=type2d, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_patch=ptr_patch)
          else
             ptr_1d => ptr_patch(:,1)
             call hist_addfld1d (fname=trim(fname), units=units, &
                  avgflag=avgflag, long_name=long_name, &
                  ptr_patch=ptr_1d)
          endif
       endif

    else
       write(iulog, *) ' error: hist_addfld_decomp needs either pft or column level pointer'
       write(iulog, *) fname
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end subroutine hist_addfld_decomp

  !-----------------------------------------------------------------------
  integer function pointer_index ()
    !
    ! !DESCRIPTION:
    ! Set the current pointer index and increment the value of the index.
    !
    ! !ARGUMENTS:
    !
    integer, save :: lastindex = 1
    character(len=*),parameter :: subname = 'pointer_index'
    !-----------------------------------------------------------------------

    pointer_index = lastindex
    lastindex = lastindex + 1
    if (lastindex > max_mapflds) then
       write(iulog,*) trim(subname),' ERROR: ',&
            ' lastindex = ',lastindex,' greater than max_mapflds= ',max_mapflds
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end function pointer_index

  !-----------------------------------------------------------------------
  subroutine hist_add_subscript(name, dim)
    !
    ! !DESCRIPTION:
    ! Add a history variable to the output history tape.
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: name ! name of subscript
    integer         , intent(in) :: dim  ! dimension of subscript
    !
    ! !LOCAL VARIABLES:
    character(len=*),parameter :: subname = 'hist_add_subscript'
    !-----------------------------------------------------------------------

    num_subs = num_subs + 1
    if (num_subs > max_subs) then
       write(iulog,*) trim(subname),' ERROR: ',&
            ' num_subs = ',num_subs,' greater than max_subs= ',max_subs
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    subs_name(num_subs) = name
    subs_dim(num_subs) =  dim

  end subroutine hist_add_subscript

  !-----------------------------------------------------------------------

  subroutine strip_null(str)
    character(len=*), intent(inout) :: str
    integer :: i	
    do i=1,len(str)
       if(ichar(str(i:i))==0) str(i:i)=' '
    end do
  end subroutine strip_null
  
  !------------------------------------------------------------------------
  subroutine hist_do_disp (ntapes, hist_ntimes, hist_mfilt, if_stop, if_disphist, rstwr, nlend)
    !
    ! !DESCRIPTION:
    ! Determine logic for closeing and/or disposing history file
    ! Sets values for if_disphist, if_stop (arguments)
    ! Remove history files unless this is end of run or
    ! history file is not full.
    !
    ! !USES:
    use clm_time_manager, only : is_last_step
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: ntapes              !actual number of history tapes
    integer, intent(in)  :: hist_ntimes(ntapes) !current numbers of time samples on history tape
    integer, intent(in)  :: hist_mfilt(ntapes)  !maximum number of time samples per tape
    logical, intent(out) :: if_stop             !true => last time step of run
    logical, intent(out) :: if_disphist(ntapes) !true => save and dispose history file
    logical, intent(in)  :: rstwr
    logical, intent(in)  :: nlend	
    !
    ! !LOCAL VARIABLES:
    integer :: t                   ! history tape index
    logical :: rest_now            ! temporary
    logical :: stop_now            ! temporary
    !------------------------------------------------------------------------

    rest_now = .false.
    stop_now = .false.
    
    if (nlend) stop_now = .true.
    if (rstwr) rest_now = .true.
    
    if_stop = stop_now
    
    if (stop_now) then
       ! End of run -  dispose all history files
       
       if_disphist(1:ntapes) = .true.
       
    else if (rest_now) then
       ! Restart - dispose all history files
       
       do t = 1,ntapes
          if_disphist(t) = .true.
       end do
    else
       ! Dispose
       
       if_disphist(1:ntapes) = .false.
       do t = 1,ntapes
          if (hist_ntimes(t) ==  hist_mfilt(t)) then
             if_disphist(t) = .true.
          endif
       end do
    endif
    
  end subroutine hist_do_disp

end module histFileMod

