module RtmHistFile
!-----------------------------------------------------------------------
! !MODULE: RtmHistFileMod
!
! !DESCRIPTION:
! Module containing methods to for RTM history file handling.
!
! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_sys_mod   , only : shr_sys_flush, shr_sys_abort
  use RunoffMod     , only : runoff
  use RtmVar        , only : rtmlon, rtmlat, spval, ispval, secspday, frivinp_rtm, &   
                             iulog, nsrest, caseid, inst_suffix, nsrStartup, nsrBranch, & 
                             ctitle, version, hostname, username, conventions, source
  use RtmFileUtils  , only : get_filename, getfil
  use RtmTimeManager, only : get_nstep, get_curr_date, get_curr_time, get_ref_date, &
                             get_prev_time, get_prev_date, is_last_step
  use RtmSpmd       , only : masterproc
  use RtmIO
  use RtmDateTime

  implicit none
  save
  private

!
! !PUBLIC TYPES:
!
! Constants
!
  integer , public, parameter :: max_tapes  = 3     ! max number of history tapes
  integer , public, parameter :: max_flds   = 1500  ! max number of history fields
  integer , public, parameter :: max_namlen = 32    ! maximum number of characters for field name
!
! Counters
!
  integer , public :: ntapes = 0         ! index of max history file requested
!
! Namelist
!
  integer :: ni
  integer, public :: &
       rtmhist_ndens(max_tapes) = 2         ! namelist: output density of netcdf history files
  integer, public :: &
       rtmhist_mfilt(max_tapes) = 30        ! namelist: number of time samples per tape
  integer, public :: &
       rtmhist_nhtfrq(max_tapes) = (/0, -24, -24/)  ! namelist: history write freq(0=monthly)
  character(len=1), public :: &
       rtmhist_avgflag_pertape(max_tapes) = (/(' ',ni=1,max_tapes)/)   ! namelist: per tape averaging flag
  
  ! list of fields to add
  character(len=max_namlen+2), public :: rtmhist_fincl1(max_flds) = ' '       
  character(len=max_namlen+2), public :: rtmhist_fincl2(max_flds) = ' '
  character(len=max_namlen+2), public :: rtmhist_fincl3(max_flds) = ' '

  ! list of fields to remove
  character(len=max_namlen+2), public :: rtmhist_fexcl1(max_flds) = ' ' 
  character(len=max_namlen+2), public :: rtmhist_fexcl2(max_flds) = ' ' 
  character(len=max_namlen+2), public :: rtmhist_fexcl3(max_flds) = ' ' 

  ! equivalence list of fields to add/remove
  character(len=max_namlen+2), public :: fexcl(max_flds,max_tapes)         
  character(len=max_namlen+2), public :: fincl(max_flds,max_tapes)         

!! Restart
!
  logical, private :: if_close(max_tapes)   ! true => save history file
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: RtmHistAddfld        ! Add a field to the master field list
  public :: RtmHistPrintflds     ! Print summary of master field list
  public :: RtmHistHtapesBuild   ! Initialize history file handler for initial or continue run
  public :: RtmHistUpdateHbuf    ! Updates history buffer for all fields and tapes
  public :: RtmHistHtapesWrapup  ! Write history tape(s)
  public :: RtmHistRestart       ! Read/write history file restart data
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: htapes_fieldlist          ! Define the contents of each history file based on namelist
  private :: htape_addfld              ! Add a field to the active list for a history tape
  private :: htape_create              ! Define contents of history file t
  private :: htape_timeconst           ! Write time constant values to history tape
  private :: set_hist_filename         ! Determine history dataset filenames
  private :: list_index                ! Find index of field in exclude list
  private :: getname                   ! Retrieve name portion of input "inname"
  private :: getflag                   ! Retrieve flag
  private :: max_nFields               ! The max number of fields on any tape

! !PRIVATE TYPES:
! Constants
!
  integer, parameter :: max_chars = 128        ! max chars for char variables
!
! Subscript dimensions
!
  integer, parameter :: max_subs = 100         ! max number of subscripts
  character(len=32)  :: subs_name(max_subs)    ! name of subscript
  integer            :: subs_dim(max_subs)     ! dimension of subscript
!
! Derived types
!
  type field_info
     character(len=max_namlen) :: name         ! field name
     character(len=max_chars)  :: long_name    ! long name
     character(len=max_chars)  :: units        ! units
     integer :: hpindex                        ! history pointer index 
  end type field_info

  type master_entry
     type (field_info)  :: field               ! field information
     logical            :: actflag(max_tapes)  ! active/inactive flag
     character(len=1)   :: avgflag(max_tapes)  ! time averaging flag ("X","A","M" or "I",)
  end type master_entry

  type history_entry
     type (field_info) :: field                ! field information
     character(len=1)  :: avgflag              ! time averaging flag
     real(r8), pointer :: hbuf(:)              ! history buffer (dimensions: dim1d x 1)
     integer , pointer :: nacs(:)              ! accumulation counter (dimensions: dim1d x 1)
  end type history_entry

  type history_tape
     integer  :: nflds                         ! number of active fields on tape
     integer  :: ntimes                        ! current number of time samples on tape
     integer  :: mfilt                         ! maximum number of time samples per tape
     integer  :: nhtfrq                        ! number of time samples per tape
     integer  :: ncprec                        ! netcdf output precision
     logical  :: is_endhist                    ! true => current time step is end of history interval
     real(r8) :: begtime                       ! time at beginning of history averaging interval
     type (history_entry) :: hlist(max_flds)   ! array of active history tape entries
  end type history_tape

  type rtmpoint                             ! Pointer to real scalar data (1D)
     real(r8), pointer :: ptr(:)
  end type rtmpoint
!EOP
!
! Pointers
!
  integer, parameter :: max_mapflds = 1500     ! Maximum number of fields to track
  type (rtmpoint)    :: rtmptr(max_mapflds) ! Real scalar data (1D)
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
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

  subroutine RtmHistPrintflds()

    ! DESCRIPTION:
    ! Print summary of master field list.

    ! !ARGUMENTS:
    implicit none

    ! !LOCAL VARIABLES:
    integer nf
    character(len=*),parameter :: subname = 'RTM_hist_printflds'

    if (masterproc) then
       write(iulog,*) trim(subname),' : number of master fields = ',nfmaster
       write(iulog,*)' ******* MASTER FIELD LIST *******'
       do nf = 1,nfmaster
          write(iulog,9000)nf, masterlist(nf)%field%name, masterlist(nf)%field%units
9000      format (i5,1x,a32,1x,a16)
       end do
       call shr_sys_flush(iulog)
    end if

  end subroutine RtmHistPrintflds

!-----------------------------------------------------------------------

  subroutine RtmHistHtapesBuild ()

    ! !DESCRIPTION:
    ! Initialize ntapes history file for initial or branch run.  

    ! !ARGUMENTS:
    implicit none

    ! !LOCAL VARIABLES:
    integer :: i                   ! index
    integer :: ier                 ! error code
    integer :: t, f                ! tape, field indices
    integer :: day, sec            ! day and seconds from base date
    character(len=1) :: avgflag    ! lcl equiv of rtmhist_avgflag_pertape(t)
    character(len=*),parameter :: subname = 'hist_htapes_build'
    !----------------------------------------------------------

    if (masterproc) then
       write(iulog,*)  trim(subname),' Initializing RTM history files'
       write(iulog,'(72a1)') ("-",i=1,60)
       call shr_sys_flush(iulog)
    endif

    ! Override averaging flag for all fields on a particular tape
    ! if namelist input so specifies

    do t=1,max_tapes
       if (rtmhist_avgflag_pertape(t) /= ' ') then
          avgflag = rtmhist_avgflag_pertape(t)
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
                call shr_sys_abort ()
             end select
          end do
       end if
    end do

    fincl(:,1) = rtmhist_fincl1(:)
    fincl(:,2) = rtmhist_fincl2(:)
    fincl(:,3) = rtmhist_fincl3(:)

    fexcl(:,1) = rtmhist_fexcl1(:)
    fexcl(:,2) = rtmhist_fexcl2(:)
    fexcl(:,3) = rtmhist_fexcl3(:)

    ! Define field list information for all history files.
    ! Update ntapes to reflect number of active history files
    ! Note - branch runs can have additional auxiliary history files declared

    call htapes_fieldlist()

    ! Set number of time samples in each history file and
    ! Note - the following entries will be overwritten by history restart
    ! Note - with netcdf, only 1 (ncd_double) and 2 (ncd_float) are allowed

    do t=1,ntapes
       tape(t)%ntimes = 0
       tape(t)%nhtfrq = rtmhist_nhtfrq(t)
       tape(t)%mfilt = rtmhist_mfilt(t)
       if (rtmhist_ndens(t) == 1) then
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
       write(iulog,*)  trim(subname),' Successfully initialized RTM history files'
       write(iulog,'(72a1)') ("-",i=1,60)
       call shr_sys_flush(iulog)
    endif

  end subroutine RtmHistHtapesBuild

!-----------------------------------------------------------------------

  subroutine htapes_fieldlist()

    ! !DESCRIPTION:
    ! Define the contents of each history file based on namelist
    ! input for initial or branch run, and restart data if a restart run.
    ! Use arrays fincl and fexcl to modify default history tape contents.
    ! Then sort the result alphanumerically.

    ! !ARGUMENTS:
    implicit none
    !
    ! !LOCAL VARIABLES:
    integer :: t, f                ! tape, field indices
    integer :: ff                  ! index into include, exclude and fprec list
    character(len=max_namlen) :: name ! field name portion of fincl (i.e. no avgflag separator)
    character(len=max_namlen) :: mastername ! name from masterlist field
    character(len=1)  :: avgflag    ! averaging flag
    character(len=1)  :: prec_acc   ! history buffer precision flag
    character(len=1)  :: prec_wrt   ! history buffer write precision flag
    type (history_entry) :: tmp     ! temporary used for swapping
    character(len=*),parameter :: subname = 'htapes_fieldlist'
    !---------------------------------------------------------

    ! First ensure contents of fincl and fexcl are valid names
    do t = 1,max_tapes
       f = 1
       do while (f < max_flds .and. fincl(f,t) /= ' ')
          name = getname (fincl(f,t)) !namelist
          do ff = 1,nfmaster
             mastername = masterlist(ff)%field%name
             if (name == mastername) exit
          end do
          if (name /= mastername) then
             write(iulog,*) trim(subname),' ERROR: ', trim(name), ' in fincl(', f, ') ',&
                  'for history tape ',t,' not found'
             call shr_sys_abort()
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
             call shr_sys_abort()
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
          else 
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
                call shr_sys_abort()
             end if
          end do
       end do

       if (masterproc) then
          if (tape(t)%nflds > 0) then
             write(iulog,*) trim(subname),' : Included fields tape ',t,'=',tape(t)%nflds
          end if
          do f = 1,tape(t)%nflds
             write(iulog,*) f,' ',tape(t)%hlist(f)%field%name,' ',tape(t)%hlist(f)%avgflag
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
          call shr_sys_abort()
       end if
    end do

    ! Check that the number of history files declared does not exceed
    ! the maximum allowed.

    if (ntapes > max_tapes) then
       write(iulog,*) trim(subname),' ERROR: Too many history files declared, max_tapes=',max_tapes
       call shr_sys_abort()
    end if

    if (masterproc) then
       write(iulog,*) 'There will be a total of ',ntapes,'RTM  history tapes'
       do t=1,ntapes
          write(iulog,*)
          if (rtmhist_nhtfrq(t) == 0) then
             write(iulog,*)'RTM History tape ',t,' write frequency is MONTHLY'
          else
             write(iulog,*)'RTM History tape ',t,' write frequency = ',rtmhist_nhtfrq(t)
          endif
          write(iulog,*)'Number of time samples on RTM history tape ',t,' is ',rtmhist_mfilt(t)
          write(iulog,*)'Output precision on RTM history tape ',t,'=',rtmhist_ndens(t)
          write(iulog,*)
       end do
       call shr_sys_flush(iulog)
    end if

    ! Set flag indicating h-tape contents are now defined 

    htapes_defined = .true.

  end subroutine htapes_fieldlist

!-----------------------------------------------------------------------

  subroutine htape_addfld (t, f, avgflag)

    ! !DESCRIPTION:
    ! Add a field to the active list for a history tape. Copy the data from
    ! the master field list to the active list for the tape.

    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t                 ! history tape index
    integer, intent(in) :: f                 ! field index from master field list
    character(len=1), intent(in) :: avgflag  ! time averaging flag

    ! !LOCAL VARIABLES:
    integer :: n               ! field index on defined tape
    integer :: begrof          ! per-proc beginning land runoff index
    integer :: endrof          ! per-proc ending land runoff index
    integer :: numrtm          ! total number of rtm cells on all procs
    character(len=*),parameter :: subname = 'htape_addfld'
    !-------------------------------------------------------

    ! Ensure that it is not to late to add a field to the history tape
    if (htapes_defined) then
       write(iulog,*) trim(subname),' ERROR: attempt to add field ', &
            masterlist(f)%field%name, ' after history files are set'
       call shr_sys_abort()
    end if

    ! Determine bounds
    begrof = runoff%begr
    endrof = runoff%endr
    numrtm = runoff%numr

    tape(t)%nflds = tape(t)%nflds + 1
    n = tape(t)%nflds
    tape(t)%hlist(n)%field = masterlist(f)%field
    allocate (tape(t)%hlist(n)%hbuf(begrof:endrof))
    allocate (tape(t)%hlist(n)%nacs(begrof:endrof))
    tape(t)%hlist(n)%hbuf(:) = 0._r8
    tape(t)%hlist(n)%nacs(:) = 0

    ! Set time averaging flag based on masterlist setting or
    ! override the default averaging flag with namelist setting
    select case (avgflag)
    case (' ')
       tape(t)%hlist(n)%avgflag = masterlist(f)%avgflag(t)
    case ('A','I','X','M')
       tape(t)%hlist(n)%avgflag = avgflag
    case default
       write(iulog,*) trim(subname),' ERROR: unknown avgflag=', avgflag
       call shr_sys_abort()
    end select

  end subroutine htape_addfld

!-----------------------------------------------------------------------

  subroutine RtmHistUpdateHbuf()

    ! !DESCRIPTION:
    ! Accumulate (or take min, max, etc. as appropriate) input field
    ! into its history buffer for appropriate tapes.

    ! !ARGUMENTS:
    implicit none

    ! !LOCAL VARIABLES:
    integer :: t                   ! tape index
    integer :: f                   ! field index
    integer :: k                   ! index
    integer :: hpindex             ! history pointer index
    integer :: begrof,endrof       ! beginning and ending indices
    character(len=1)  :: avgflag   ! time averaging flag
    real(r8), pointer :: hbuf(:)   ! history buffer
    integer , pointer :: nacs(:)   ! accumulation counter
    real(r8), pointer :: field(:)  ! 1d pointer field
    integer j
    character(len=*),parameter :: subname = 'RtmHistUpdateHbuf'
    !----------------------------------------------------------

    begrof = runoff%begr
    endrof = runoff%endr

    do t = 1,ntapes
       do f = 1,tape(t)%nflds
          avgflag  =  tape(t)%hlist(f)%avgflag
          nacs     => tape(t)%hlist(f)%nacs
          hbuf     => tape(t)%hlist(f)%hbuf
          hpindex  =  tape(t)%hlist(f)%field%hpindex
          field    => rtmptr(hpindex)%ptr

          select case (avgflag)
          case ('I') ! Instantaneous
             do k = begrof,endrof
                if (field(k) /= spval) then
                   hbuf(k) = field(k)
                else
                   hbuf(k) = spval
                end if
                nacs(k) = 1
             end do
          case ('A') ! Time average
             do k = begrof,endrof
                if (field(k) /= spval) then
                   if (nacs(k) == 0) hbuf(k) = 0._r8
                   hbuf(k) = hbuf(k) + field(k)
                   nacs(k) = nacs(k) + 1
                else
                   if (nacs(k) == 0) hbuf(k) = spval
                end if
             end do
          case ('X') ! Maximum over time
             do k = begrof,endrof
                if (field(k) /= spval) then
                   if (nacs(k) == 0) hbuf(k) = -1.e50_r8
                   hbuf(k) = max( hbuf(k), field(k) )
                else
                   if (nacs(k) == 0) hbuf(k) = spval
                end if
                nacs(k) = 1
             end do
          case ('M') ! Minimum over time
             do k = begrof,endrof
                if (field(k) /= spval) then
                   if (nacs(k) == 0) hbuf(k) = +1.e50_r8
                   hbuf(k) = min( hbuf(k), field(k) )
                else
                   if (nacs(k) == 0) hbuf(k) = spval
                end if
                nacs(k) = 1
             end do
          case default
             write(iulog,*) trim(subname),' ERROR: invalid time averaging flag ', avgflag
             call shr_sys_abort()
          end select
       end do
    end do

  end subroutine RtmHistUpdateHbuf

!-----------------------------------------------------------------------

  subroutine htape_create (t, histrest)

    ! !DESCRIPTION:
    ! Define contents of history file t. Issue the required netcdf
    ! wrapper calls to define the history file contents.
    !

    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t                   ! tape index
    logical, intent(in), optional :: histrest  ! if creating the history restart file

    ! !LOCAL VARIABLES:
    integer :: f                   ! field index
    integer :: p,c,l,n             ! indices
    integer :: ier                 ! error code
    integer :: dimid               ! dimension id temporary
    integer :: dim1id(1)           ! netCDF dimension id
    integer :: dim2id(2)           ! netCDF dimension id
    integer :: ndims               ! dimension counter
    integer :: omode               ! returned mode from netCDF call
    integer :: ncprec              ! output netCDF write precision
    integer :: ret                 ! netCDF error status
    integer :: num_lndrof          ! total number of land runoff across all procs
    integer :: num_ocnrof          ! total number of ocean runoff across all procs
    integer :: numrtm              ! total number of rtm cells on all procs
    logical :: lhistrest           ! local history restart flag
    type(file_desc_t) :: lnfid     ! local file id
    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time
    character(len=256) :: name     ! name of attribute
    character(len=256) :: units    ! units of attribute
    character(len=256) :: str      ! global attribute string
    character(len=  1) :: avgflag  ! time averaging flag
    character(len=*),parameter :: subname = 'htape_create'
    !-----------------------------------------------------

    if ( present(histrest) )then
       lhistrest = histrest
    else
       lhistrest = .false.
    end if

    ! Define output write precsion for tape
    ncprec = tape(t)%ncprec

    ! Create new netCDF file. It will be in define mode
    if ( .not. lhistrest )then
       if (masterproc) then
          write(iulog,*) trim(subname),' : Opening netcdf htape ', &
                                      trim(locfnh(t))
          call shr_sys_flush(iulog)
       end if
       call ncd_pio_createfile(lnfid, trim(locfnh(t)))
       call ncd_putatt(lnfid, ncd_global, 'title', 'RTM History file information' )
       call ncd_putatt(lnfid, ncd_global, 'comment', &
          "NOTE: None of the variables are weighted by land fraction!" )
    else
       if (masterproc) then
          write(iulog,*) trim(subname),' : Opening netcdf rhtape ', &
                                      trim(locfnhr(t))
          call shr_sys_flush(iulog)
       end if
       call ncd_pio_createfile(lnfid, trim(locfnhr(t)))
       call ncd_putatt(lnfid, ncd_global, 'title', &
            'RTM Restart History information, required to continue a simulation' )
       call ncd_putatt(lnfid, ncd_global, 'comment', &
            "This entire file NOT needed for startup or branch simulations")
    end if

    ! Create global attributes. Attributes are used to store information
    ! about the data set. Global attributes are information about the
    ! data set as a whole, as opposed to a single variable

    call ncd_putatt(lnfid, ncd_global, 'Conventions', trim(conventions))
    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call ncd_putatt(lnfid, ncd_global, 'history' , trim(str))
    call ncd_putatt(lnfid, ncd_global, 'source'  , trim(source))
    call ncd_putatt(lnfid, ncd_global, 'hostname', trim(hostname))
    call ncd_putatt(lnfid, ncd_global, 'username', trim(username))
    call ncd_putatt(lnfid, ncd_global, 'version' , trim(version))

    str = &
    '$Id: histFileMod.F90 36692 2012-04-27 18:39:55Z tcraig $'
    call ncd_putatt(lnfid, ncd_global, 'revision_id', trim(str))
    call ncd_putatt(lnfid, ncd_global, 'case_title', trim(ctitle))
    call ncd_putatt(lnfid, ncd_global, 'case_id', trim(caseid))

    str = get_filename(frivinp_rtm)
    call ncd_putatt(lnfid, ncd_global, 'RTM_input_dataset', trim(str))

    ! Define dimensions.
    ! Time is an unlimited dimension. Character string is treated as an array of characters.

    ! Global uncompressed dimensions (including non-land points)
    numrtm     = runoff%numr
    num_lndrof = runoff%numrl
    num_ocnrof = runoff%numro
    call ncd_defdim( lnfid, 'lon', rtmlon    , dimid)
    call ncd_defdim( lnfid, 'lat', rtmlat    , dimid)
    call ncd_defdim( lnfid, 'ocnrof', num_ocnrof, dimid)
    call ncd_defdim( lnfid, 'lndrof', num_lndrof, dimid)
    call ncd_defdim( lnfid, 'allrof', numrtm    , dimid)

    call ncd_defdim(lnfid, 'string_length', 8, strlen_dimid)

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

  subroutine htape_timeconst(t, mode)

    ! !DESCRIPTION:
    ! Write time constant values to primary history tape.

    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t              ! tape index
    character(len=*), intent(in) :: mode  ! 'define' or 'write'

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
    type(Var_desc_t) :: vardesc           ! netCDF variable description
    character(len=max_chars) :: long_name ! variable long name
    character(len=max_namlen):: varname   ! variable name
    character(len=max_namlen):: units     ! variable units
    character(len=256):: str              ! global attribute string
    integer :: status
    logical, save :: writeTC = .true. ! true => write out time-constant data
    character(len=*),parameter :: subname = 'htape_timeconst'
    !--------------------------------------------------------

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
       call ncd_putatt(nfid(t), varid, 'calendar', 'noleap')
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

       call ncd_defvar(varname='lon', xtype=tape(t)%ncprec, dim1name='lon', &
            long_name='runoff coordinate longitude', units='degrees_east', ncid=nfid(t))
       call ncd_defvar(varname='lat', xtype=tape(t)%ncprec, dim1name='lat', &
            long_name='runoff coordinate latitude', units='degrees_north', ncid=nfid(t))

       if(writeTC .and. t == 1) then 
          call ncd_defvar(varname='fthresh', xtype=tape(t)%ncprec, dim1name='lon', &
               dim2name='lat', long_name='flooding threshold', &
               missing_value=spval, fill_value=spval, units='m3', ncid=nfid(t))
       endif
    else if (mode == 'write') then

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

       call ncd_io(varname='lon', data=runoff%rlon, ncid=nfid(t), flag='write')
       call ncd_io(varname='lat', data=runoff%rlat, ncid=nfid(t), flag='write')

       if(writeTC .and. t == 1 .and. tape(t)%ntimes == 1 ) then 
          call ncd_io(varname='fthresh', data=runoff%fthresh, ncid=nfid(t), flag='write' &
               ,dim1name='allrof')
          writeTC = .false.
       endif
    endif

  end subroutine htape_timeconst

!-----------------------------------------------------------------------

  subroutine RtmHistHtapesWrapup( rstwr, nlend )

    ! DESCRIPTION:
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


    ! !ARGUMENTS:
    implicit none
    logical, intent(in) :: rstwr   ! true => write restart file this step
    logical, intent(in) :: nlend   ! true => end of run on this step

    ! !LOCAL VARIABLES:
    integer :: begrof, endrof      ! beg and end rof indices
    integer :: t,f,k,nt            ! indices
    integer :: nstep               ! current step
    integer :: day                 ! current day (1 -> 31)
    integer :: mon                 ! current month (1 -> 12)
    integer :: yr                  ! current year (0 -> ...)
    integer :: mdcur               ! current day
    integer :: mscur               ! seconds of current day
    integer :: mcsec               ! current time of day [seconds]
    integer :: daym1               ! nstep-1 day (1 -> 31)
    integer :: monm1               ! nstep-1 month (1 -> 12)
    integer :: yrm1                ! nstep-1 year (0 -> ...)
    integer :: mcsecm1             ! nstep-1 time of day [seconds]
    real(r8):: time                ! current time
    character(len=256):: str       ! global attribute string
    character(len=1)  :: avgflag   ! averaging flag
    real(r8), pointer :: histo(:)  ! temporary
    real(r8), pointer :: hbuf(:)   ! history buffer
    integer , pointer :: nacs(:)   ! accumulation counter
    character(len=32) :: avgstr    ! time averaging type
    character(len=max_chars) :: long_name ! long name
    character(len=max_chars) :: units     ! units
    character(len=max_namlen):: varname   ! variable name
    character(len=*),parameter :: subname = 'hist_htapes_wrapup'
    !-----------------------------------------------------------

    begrof = runoff%begr
    endrof = runoff%endr

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
          if (mon /= monm1) then
             tape(t)%is_endhist = .true.
          end if
       else
          if (mod(nstep,tape(t)%nhtfrq) == 0) then
             tape(t)%is_endhist = .true.
          end if
       end if

       ! If end of history interval
       if (tape(t)%is_endhist) then

          ! Normalize by number of accumulations for time averaged case
          do f = 1,tape(t)%nflds
             avgflag =  tape(t)%hlist(f)%avgflag
             nacs    => tape(t)%hlist(f)%nacs
             hbuf    => tape(t)%hlist(f)%hbuf
             do k = begrof, endrof
                if ((avgflag == 'A') .and. nacs(k) /= 0) then
                   hbuf(k) = hbuf(k) / float(nacs(k))
                end if
             end do
          end do
          
          ! Increment current time sample counter.
          tape(t)%ntimes = tape(t)%ntimes + 1

          ! Create history file if appropriate and build time comment

          ! If first time sample, generate unique history file name, open file,
          ! define dims, vars, etc.

          if (tape(t)%ntimes == 1) then
             locfnh(t) = set_hist_filename (hist_freq=tape(t)%nhtfrq, &
                                            rtmhist_mfilt=tape(t)%mfilt, hist_file=t)
             if (masterproc) then
                write(iulog,*) trim(subname),' : Creating history file ', trim(locfnh(t)), &
                     ' at nstep = ',get_nstep()
                write(iulog,*)'calling htape_create for file t = ',t
             endif
             call htape_create (t)

             ! Define time-constant field variables
             call htape_timeconst(t, mode='define')

             ! Define model field variables

             do f = 1,tape(t)%nflds
                varname    = tape(t)%hlist(f)%field%name
                long_name  = tape(t)%hlist(f)%field%long_name
                units      = tape(t)%hlist(f)%field%units
                avgflag    = tape(t)%hlist(f)%avgflag
                
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
                   write(iulog,*) trim(subname),&
                        ' ERROR: unknown time averaging flag (avgflag)=',avgflag
                   call shr_sys_abort()
                end select
                
                call ncd_defvar(ncid=nfid(t), varname=varname, xtype=tape(t)%ncprec, &
                     dim1name='lon', dim2name='lat', dim3name='time', &
                     long_name=long_name, units=units, cell_method=avgstr, &
                     missing_value=spval, fill_value=spval)
             end do
                
             ! Exit define model
             call ncd_enddef(nfid(t))

          endif

          ! Write time constant history variables
          call htape_timeconst(t, mode='write')

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

          ! Write history time slice
          do f = 1,tape(t)%nflds
             varname =  tape(t)%hlist(f)%field%name
             nt      =  tape(t)%ntimes
             histo   => tape(t)%hlist(f)%hbuf
             call ncd_io(flag='write', varname=varname, dim1name='allrof', &
                  data=histo, ncid=nfid(t), nt=nt)
          end do

          ! Zero necessary history buffers
          do f = 1,tape(t)%nflds
             tape(t)%hlist(f)%hbuf(:) = 0._r8
             tape(t)%hlist(f)%nacs(:) = 0
          end do

       end if

    end do  ! end loop over history tapes

    ! Close open history files
    ! Auxilary files may have been closed and saved off without being full,
    ! must reopen the files

    do t = 1, ntapes
       if (nlend) then      
          if_close(t) = .true.
       else if (rstwr) then 
          if_close(t) = .true.
       else                  
          if (tape(t)%ntimes == tape(t)%mfilt) then
             if_close(t) = .true.
          else
             if_close(t) = .false.
          end if
       endif
       if (if_close(t)) then
          if (tape(t)%ntimes /= 0) then
             if (masterproc) then
                write(iulog,*)
                write(iulog,*)  trim(subname),' : Closing local history file ',&
                     trim(locfnh(t)),' at nstep = ', get_nstep()
                write(iulog,*)
             endif
	     call ncd_pio_closefile(nfid(t))
             if ((.not.nlend) .and. (tape(t)%ntimes/=tape(t)%mfilt)) then
                call ncd_pio_openfile (nfid(t), trim(locfnh(t)), ncd_write)
             end if
          else
             if (masterproc) then
                write(iulog,*) trim(subname),' : history tape ',t,': no open file to close'
             end if
          endif
          if (tape(t)%ntimes==tape(t)%mfilt) then
             tape(t)%ntimes = 0
          end if
       endif
    end do

  end subroutine RtmHistHtapesWrapup

!-----------------------------------------------------------------------

  subroutine RtmHistRestart (ncid, flag, rdate)

    ! !DESCRIPTION:
    ! Read/write history file restart data.
    ! If the current history file(s) are not full, file(s) are opened
    ! so that subsequent time samples are added until the file is full.
    ! A new history file is used on a branch run.

    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid     ! netcdf file
    character(len=*) , intent(in)    :: flag     !'read' or 'write'
    character(len=*) , intent(in), optional :: rdate    ! restart file time stamp for name

    ! !LOCAL VARIABLES:
    integer :: max_nflds                         ! max number of fields
    integer :: begrof                            ! per-proc beginning ocean runoff index
    integer :: endrof                            ! per-proc ending ocean runoff index
    character(len=max_namlen) :: name            ! variable name
    character(len=max_namlen) :: name_acc        ! accumulator variable name
    character(len=max_namlen) :: long_name       ! long name of variable
    character(len=max_chars)  :: long_name_acc   ! long name for accumulator
    character(len=max_chars)  :: units           ! units of variable
    character(len=max_chars)  :: units_acc       ! accumulator units
    character(len=max_chars)  :: fname           ! full name of history file
    character(len=max_chars)  :: locrest(max_tapes) ! local history restart file names
    character(len=1)   :: hnum                   ! history file index
    type(var_desc_t)   :: name_desc              ! variable descriptor for name
    type(var_desc_t)   :: longname_desc          ! variable descriptor for long_name
    type(var_desc_t)   :: units_desc             ! variable descriptor for units
    type(var_desc_t)   :: avgflag_desc           ! variable descriptor for avgflag
    integer :: status                            ! error status
    integer :: dimid                             ! dimension ID
    integer :: start(2)                          ! Start array index
    integer :: k                                 ! 1d index
    integer :: t                                 ! tape index
    integer :: f                                 ! field index
    integer :: varid                             ! variable id
    integer, allocatable :: itemp2d(:,:)         ! 2D temporary
    real(r8), pointer :: hbuf(:)                 ! history buffer
    integer , pointer :: nacs(:)                 ! accumulation counter
    character(len=*),parameter :: subname = 'hist_restart_ncd'
    !---------------------------------------------------------

    ! If branch run, initialize file times and return

    if (flag == 'read') then
       if (nsrest == nsrBranch) then
          do t = 1,ntapes
             tape(t)%ntimes = 0
          end do
          RETURN
       end if
       ! If startup run just return
       if (nsrest == nsrStartup) then
          RETURN
       end if
    endif

    ! Read history file data only for restart run (not for branch run)

    ! First when writing out and in define mode, create files and define all variables
    !================================================
    if (flag == 'define') then
    !================================================

       if (.not. present(rdate)) then
          call shr_sys_abort('variable rdate must be present for writing restart files')
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
       call ncd_defvar(ncid=ncid, varname='locfnhr', xtype=ncd_char, &
            long_name="Restart history filename",     &
            comment="This variable NOT needed for startup or branch simulations", &
            dim1name='max_chars', dim2name="ntapes" )

       ! max_nflds is the maximum number of fields on any tape
       ! max_flds is the maximum number possible number of fields 
       max_nflds = max_nFields()

       ! Loop over tapes - write out namelist information to each restart-history tape
       ! only read/write accumulators and counters if needed

       do t = 1,ntapes
          !
          ! Create the restart history filename and open it
          !
          write(hnum,'(i1.1)') t-1
          locfnhr(t) = "./" // trim(caseid) //".rtm"// trim(inst_suffix) &
                        // ".rh" // hnum //"."// trim(rdate) //".nc"
          call htape_create( t, histrest=.true. )
          !
          ! Add read/write accumultators and counters if needed
          !
          if (.not. tape(t)%is_endhist) then
             do f = 1,tape(t)%nflds
                name           =  tape(t)%hlist(f)%field%name
                long_name      =  tape(t)%hlist(f)%field%long_name
                units          =  tape(t)%hlist(f)%field%units
                name_acc       =  trim(name) // "_acc"
                units_acc      =  "unitless positive integer"
                long_name_acc  =  trim(long_name) // " accumulator number of samples"
                nacs           => tape(t)%hlist(f)%nacs
                hbuf           => tape(t)%hlist(f)%hbuf
               
                call ncd_defvar(ncid=ncid_hist(t), varname=trim(name), xtype=ncd_double, &
                     dim1name='lon', dim2name='lat', &
                     long_name=trim(long_name), units=trim(units))
                call ncd_defvar(ncid=ncid_hist(t), varname=trim(name_acc), xtype=ncd_int,  &
                     dim1name='lon', dim2name='lat', &
                     long_name=trim(long_name_acc), units=trim(units_acc))
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
          call ncd_enddef(ncid_hist(t))

       end do   ! end of ntapes loop   

       RETURN

    !================================================
    else if (flag == 'write') then
    !================================================

       ! Add history filenames to master restart file
       do t = 1,ntapes
          call ncd_io('locfnh',  locfnh(t),  'write', ncid, nt=t)
          call ncd_io('locfnhr', locfnhr(t), 'write', ncid, nt=t)
       end do
       
       fincl(:,1) = rtmhist_fincl1(:)
       fincl(:,2) = rtmhist_fincl2(:)
       fincl(:,3) = rtmhist_fincl3(:)

       fexcl(:,1) = rtmhist_fexcl1(:)
       fexcl(:,2) = rtmhist_fexcl2(:)
       fexcl(:,3) = rtmhist_fexcl3(:)

       max_nflds = max_nFields()

       start(1)=1


       ! Add history namelist data to each history restart tape
       allocate(itemp2d(max_nflds,ntapes))
       do t = 1,ntapes
          call ncd_inqvid(ncid_hist(t), 'name',           varid, name_desc)
          call ncd_inqvid(ncid_hist(t), 'long_name',      varid, longname_desc)
          call ncd_inqvid(ncid_hist(t), 'units',          varid, units_desc)
          call ncd_inqvid(ncid_hist(t), 'avgflag',        varid, avgflag_desc)

          call ncd_io(varname='fincl'     , data=fincl(:,t)        , ncid=ncid_hist(t), flag='write')
          call ncd_io(varname='fexcl'     , data=fexcl(:,t)        , ncid=ncid_hist(t), flag='write')
          call ncd_io(varname='is_endhist', data=tape(t)%is_endhist, ncid=ncid_hist(t), flag='write')

          itemp2d(:,:) = 0
          do f=1,tape(t)%nflds
             itemp2d(f,t) = tape(t)%hlist(f)%field%hpindex
          end do
          call ncd_io(varname='hpindex', data=itemp2d(:,t), ncid=ncid_hist(t), flag='write')

          call ncd_io('nflds' ,  tape(t)%nflds,   'write', ncid_hist(t))
          call ncd_io('ntimes',  tape(t)%ntimes,  'write', ncid_hist(t))
          call ncd_io('nhtfrq',  tape(t)%nhtfrq,  'write', ncid_hist(t))
          call ncd_io('mfilt' ,  tape(t)%mfilt,   'write', ncid_hist(t))
          call ncd_io('ncprec',  tape(t)%ncprec,  'write', ncid_hist(t))
          call ncd_io('begtime', tape(t)%begtime, 'write', ncid_hist(t))
          do f=1,tape(t)%nflds
             start(2) = f
             call ncd_io( name_desc,           tape(t)%hlist(f)%field%name,       &
                          'write', ncid_hist(t), start )
             call ncd_io( longname_desc,       tape(t)%hlist(f)%field%long_name,  &
                          'write', ncid_hist(t), start )
             call ncd_io( units_desc,          tape(t)%hlist(f)%field%units,      &
                          'write', ncid_hist(t), start )
             call ncd_io( avgflag_desc,        tape(t)%hlist(f)%avgflag,          &
                          'write', ncid_hist(t), start )
          end do
       end do
       deallocate(itemp2d)

    !================================================
    else if (flag == 'read') then
    !================================================

       call ncd_inqdlen(ncid,dimid,ntapes,   name='ntapes')
       call ncd_io('locfnh',  locfnh(1:ntapes),  'read', ncid )
       call ncd_io('locfnhr', locrest(1:ntapes), 'read', ncid )
       do t = 1,ntapes
          call strip_null(locrest(t))
          call strip_null(locfnh(t))
       end do

       ! Determine necessary indices - the following is needed if model decomposition 
       ! is different on restart
       begrof = runoff%begr
       endrof = runoff%endr
       
       start(1)=1
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
          call ncd_inqvid(ncid_hist(t), 'avgflag',        varid, avgflag_desc)

          call ncd_io(varname='fincl', data=fincl(:,t), ncid=ncid_hist(t), flag='read')
          call ncd_io(varname='fexcl', data=fexcl(:,t), ncid=ncid_hist(t), flag='read')

          call ncd_io('nflds',   tape(t)%nflds,  'read', ncid_hist(t) )
          call ncd_io('ntimes',  tape(t)%ntimes, 'read', ncid_hist(t) )
          call ncd_io('nhtfrq',  tape(t)%nhtfrq, 'read', ncid_hist(t) )
          call ncd_io('mfilt',   tape(t)%mfilt,  'read', ncid_hist(t) )
          call ncd_io('ncprec',  tape(t)%ncprec, 'read', ncid_hist(t) )
          call ncd_io('begtime', tape(t)%begtime,'read', ncid_hist(t) )

          call ncd_io(varname='is_endhist', data=tape(t)%is_endhist, ncid=ncid_hist(t), flag='read')
          call ncd_io(varname='hpindex'   , data=itemp2d(:,t)      , ncid=ncid_hist(t), flag='read')
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
             call ncd_io( avgflag_desc,        tape(t)%hlist(f)%avgflag,          &
                          'read', ncid_hist(t), start )
             call strip_null(tape(t)%hlist(f)%field%name)
             call strip_null(tape(t)%hlist(f)%field%long_name)
             call strip_null(tape(t)%hlist(f)%field%units)
             call strip_null(tape(t)%hlist(f)%avgflag)

             allocate (tape(t)%hlist(f)%hbuf(begrof:endrof), &
                       tape(t)%hlist(f)%nacs(begrof:endrof), stat=status)
             if (status /= 0) then
                write(iulog,*) trim(subname),' ERROR: allocation error for hbuf,nacs at t,f=',t,f
                call shr_sys_abort()
             endif
             tape(t)%hlist(f)%hbuf(:) = 0._r8
             tape(t)%hlist(f)%nacs(:) = 0
          end do   ! end of flds loop

          ! If history file is not full, open it

          if (tape(t)%ntimes /= 0) then
             call ncd_pio_openfile (nfid(t), trim(locfnh(t)), ncd_write)
          end if

       end do  ! end of tapes loop

       rtmhist_fincl1(:) = fincl(:,1)
       rtmhist_fincl2(:) = fincl(:,2)
       rtmhist_fincl3(:) = fincl(:,3)

       rtmhist_fexcl1(:) = fexcl(:,1)
       rtmhist_fexcl2(:) = fexcl(:,2)
       rtmhist_fexcl3(:) = fexcl(:,3)
       
       if ( allocated(itemp2d) ) deallocate(itemp2d)

    end if

    ! Read/write history file restart data.
    ! If the current history file(s) are not full, file(s) are opened
    ! so that subsequent time samples are added until the file is full.
    ! A new history file is used on a branch run.
    
    if (flag == 'write') then     

       do t = 1,ntapes
          if (.not. tape(t)%is_endhist) then
             do f = 1,tape(t)%nflds
                name       =  tape(t)%hlist(f)%field%name
                name_acc   =  trim(name) // "_acc"
                nacs       => tape(t)%hlist(f)%nacs
                hbuf       => tape(t)%hlist(f)%hbuf

                call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name), &
                     dim1name='allrof', data=hbuf)
                call ncd_io(ncid=ncid_hist(t), flag='write', varname=trim(name_acc), &
                     dim1name='allrof', data=nacs)
             end do
          end if  ! end of is_endhist block
          call ncd_pio_closefile(ncid_hist(t))
       end do  ! end of ntapes loop   

    else if (flag == 'read') then 

       ! Read history restart information if history files are not full
       do t = 1,ntapes
          if (.not. tape(t)%is_endhist) then
             do f = 1,tape(t)%nflds
                name       =  tape(t)%hlist(f)%field%name
                name_acc   =  trim(name) // "_acc"
                nacs       => tape(t)%hlist(f)%nacs
                hbuf       => tape(t)%hlist(f)%hbuf
                
                call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name), &
                     dim1name='allrof', data=hbuf)
                call ncd_io(ncid=ncid_hist(t), flag='read', varname=trim(name_acc), &
                     dim1name='allrof', data=nacs)
             end do
          end if
          call ncd_pio_closefile(ncid_hist(t))
       end do

    end if
    
  end subroutine RtmHistRestart

!-----------------------------------------------------------------------

  integer function max_nFields()

    ! DESCRIPTION:
    ! Get the maximum number of fields on all tapes.

    ! ARGUMENTS:
    implicit none

    ! LOCAL VARIABLES:
    integer :: t  ! index
    character(len=*),parameter :: subname = 'max_nFields'

    max_nFields = 0
    do t = 1,ntapes
       max_nFields = max(max_nFields, tape(t)%nflds)
    end do
    
  end function max_nFields

!-----------------------------------------------------------------------

  character(len=max_namlen) function getname (inname)

    ! DESCRIPTION:
    ! Retrieve name portion of inname. If an averaging flag separater character
    ! is present (:) in inname, lop it off.
    
    ! ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: inname
    
    integer :: length
    integer :: i
    character(len=*),parameter :: subname = 'getname'
    
    length = len (inname)
    if (length < max_namlen .or. length > max_namlen+2) then
       write(iulog,*) trim(subname),' ERROR: bad length=',length
       call shr_sys_abort()
    end if
    
    getname = ' '
    do i = 1,max_namlen
       if (inname(i:i) == ':') exit
       getname(i:i) = inname(i:i)
    end do
    
   end function getname
   
!-----------------------------------------------------------------------

   character(len=1) function getflag (inname)

     ! DESCRIPTION:
     ! Retrieve flag portion of inname. If an averaging flag separater character
     ! is present (:) in inname, return the character after it as the flag

     ! ARGUMENTS:
     implicit none
     character(len=*) inname   ! character string

     ! LOCAL VARIABLES:
     integer :: length         ! length of inname
     integer :: i              ! loop index
     character(len=*),parameter :: subname = 'getflag'

     length = len (inname)
     if (length < max_namlen .or. length > max_namlen+2) then
        write(iulog,*) trim(subname),' ERROR: bad length=',length
        call shr_sys_abort()
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

     ! ARGUMENTS:
     implicit none
     character(len=*), intent(in) :: list(max_flds)  ! input list of names, possibly ":" delimited
     character(len=max_namlen), intent(in) :: name   ! name to be searched for
     integer, intent(out) :: index                   ! index of "name" in "list"

     ! !LOCAL VARIABLES:
     character(len=max_namlen) :: listname           ! input name with ":" stripped off.
     integer f                                       ! field index
     character(len=*),parameter :: subname = 'list_index'

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

  character(len=256) function set_hist_filename (hist_freq, rtmhist_mfilt, hist_file)

    ! Determine history dataset filenames.
    
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: hist_freq   !history file frequency
    integer, intent(in)  :: rtmhist_mfilt  !history file number of time-samples
    integer, intent(in)  :: hist_file   !history file index
    
    ! !LOCAL VARIABLES:
    character(len=256) :: cdate       !date char string
    character(len=  1) :: hist_index  !p,1 or 2 (currently)
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
    character(len=*),parameter :: subname = 'set_hist_filename'
    
    if (hist_freq == 0 .and. rtmhist_mfilt == 1) then   !monthly
       call get_prev_date (yr, mon, day, sec)
       write(cdate,'(i4.4,"-",i2.2)') yr,mon
    else                        !other
       call get_curr_date (yr, mon, day, sec)
       write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    endif
    write(hist_index,'(i1.1)') hist_file - 1
    set_hist_filename = "./"//trim(caseid)//".rtm"//trim(inst_suffix)//&
                        ".h"//hist_index//"."//trim(cdate)//".nc"

  end function set_hist_filename

!------------------------------------------------------------------------

  subroutine RtmHistAddfld (fname, units, avgflag, long_name, ptr_rof, default)

    ! Initialize a single level history field. 

    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)           :: fname          ! field name
    character(len=*), intent(in)           :: units          ! units of field
    character(len=1), intent(in)           :: avgflag        ! time averaging flag
    character(len=*), intent(in)           :: long_name      ! long name of field
    real(r8)        , pointer              :: ptr_rof(:)     ! pointer to channel runoff
    character(len=*), optional, intent(in) :: default        ! if set to 'inactive, 
                                                             ! field will not appear on primary tape

    ! !LOCAL VARIABLES:
    integer :: n              ! loop index
    integer :: f              ! masterlist index
    integer :: hpindex        ! history buffer pointer index
    logical :: found          ! flag indicates field found in masterlist
    integer, save :: lastindex = 1
    character(len=*),parameter :: subname = 'RtmHistAddfld'
    !------------------------------------------------------

    ! History buffer pointer

    hpindex = lastindex
    rtmptr(hpindex)%ptr => ptr_rof
    lastindex = lastindex + 1
    if (lastindex > max_mapflds) then
       write(iulog,*) trim(subname),' ERROR: ',&
            ' lastindex = ',lastindex,' greater than max_mapflds= ',max_mapflds
       call shr_sys_abort()
    endif

    ! Add field to masterlist

    if (fname == ' ') then
       write(iulog,*) trim(subname),' ERROR: blank field name not allowed'
       call shr_sys_abort()
    end if
    do n = 1,nfmaster
       if (masterlist(n)%field%name == fname) then
          write(iulog,*) trim(subname),' ERROR:', fname, ' already on list'
          call shr_sys_abort()
       end if
    end do
    nfmaster = nfmaster + 1
    f = nfmaster
    if (nfmaster > max_flds) then
       write(iulog,*) trim(subname),' ERROR: too many fields for primary history file ', &
            '-- max_flds,nfmaster=', max_flds, nfmaster
       call shr_sys_abort()
    end if
    masterlist(f)%field%name      = fname
    masterlist(f)%field%long_name = long_name
    masterlist(f)%field%units     = units
    masterlist(f)%field%hpindex   = hpindex

    ! The next two fields are only in master field list, NOT in runtime active field list
    ! ALL FIELDS IN THE MASTER LIST ARE INITIALIZED WITH THE ACTIVE FLAG SET TO FALSE
    masterlist(f)%avgflag(:) = avgflag
    masterlist(f)%actflag(:) = .false.

    if (present(default)) then
       if (trim(default) == 'inactive') return
    else
       ! Look through master list for input field name.
       ! When found, set active flag for that tape to true.
       found = .false.
       do f = 1,nfmaster
          if (trim(fname) == trim(masterlist(f)%field%name)) then
             masterlist(f)%actflag(1) = .true.
             found = .true.
             exit
          end if
       end do
       if (.not. found) then
          write(iulog,*) trim(subname),' ERROR: field=', fname, ' not found'
          call shr_sys_abort()
       end if
    end if

  end subroutine RtmHistAddfld

!-----------------------------------------------------------------------

  subroutine strip_null(str)
    character(len=*), intent(inout) :: str
    integer :: i	
    do i=1,len(str)
       if(ichar(str(i:i))==0) str(i:i)=' '
    end do
  end subroutine strip_null
  
end module RtmHistFile

