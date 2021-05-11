!===============================================================================
! SVN $Id: shr_file_mod.F90 2476 2006-11-15 23:16:39Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070423/shr/shr_file_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_file_mod.F90 --- Module to handle various file utilily functions.
!
! !DESCRIPTION:
!
! Miscilaneous methods to handle file and directory utilities as well as FORTRAN
! unit control. Also put/get local files into/from archival location
!
! File utilites used with CCSM Message passing:
!
! shr_file_stdio is the main example here, it changes the working directory, 
!                changes stdin and stdout to a given filename.
!
! This is needed because some implementations of MPI with MPMD so that
! each executable can run in a different working directory and redirect
! output to different files.
!
! File name archival convention, eg.
!    call shr_file_put(rcode,"foo","mss:/USER/foo",rtpd=3650)
! is extensible -- the existence of the option file name prefix, eg. "mss:", 
! and optional arguments, eg. rtpd-3650 can be used to access site-specific 
! storage devices.  Based on CCM (atmosphere) getfile & putfile routines, but
! intended to be a more extensible, shared code.
! 
! !REVISION HISTORY:
!   2006-05-08 E. Kluzek, Add in shr_file_mod and getUnit, freeUnif methods.
!   2000-??-?? B. Kauffman, original version circa 2000
!
! !INTERFACE: ------------------------------------------------------------------

MODULE shr_file_mod

! !USES:

   use shr_kind_mod  ! defines kinds
   use shr_sys_mod   ! system calls
   
   IMPLICIT none

   PRIVATE           ! By default everything is private to this module
   
! !PUBLIC TYPES:               
   
   ! no public types

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_file_put          ! Put a file to an archive location
   public :: shr_file_get          ! Get a file from an archive location
   public :: shr_file_getUnit      ! Get a logical unit for reading or writing
   public :: shr_file_freeUnit     ! Free a logical unit
   public :: shr_file_stdio        ! change dir and stdin and stdout
   public :: shr_file_chDir        ! change current working directory
   public :: shr_file_dirio        ! change stdin and stdout
   public :: shr_file_chStdIn      ! change stdin  (attach to a file)
   public :: shr_file_chStdOut     ! change stdout (attach to a file)

! !PUBLIC DATA MEMBERS:

   ! no public data members

!EOP
   integer(SHR_KIND_IN),parameter :: shr_file_nUnits = 99       ! Total number of I/O units
   logical, save        :: UnitTag(0:shr_file_nUnits) = .false. ! Logical units in use

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_put -- Put a file to an archival location.
!
! !DESCRIPTION:
!    a generic, extensible put-local-file-into-archive routine
! USAGE:
!    call shr_file_put(rcode,"foo","/home/user/foo")
!    if ( rcode /= 0 ) call shr_sys_abort( "error copying foo" )
!    call shr_file_put(rcode,"foo","cp:/home/user/foo",remove=.true.)
!    if ( rcode /= 0 ) call shr_sys_abort( "error copying foo" )
!    call shr_file_put(rcode,"foo","mss:/USER/foo",rtpd=3650)
!    if ( rcode /= 0 ) call shr_sys_abort( "error archiving foo to MSS" )
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_put(rcode,loc_fn,rem_fn,passwd,rtpd,async,remove)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(out)         :: rcode  ! return code (non-zero -- error)
   character(*),        intent(in)          :: loc_fn ! local filename
   character(*),        intent(in)          :: rem_fn ! remote filename
   character(*),        intent(in),optional :: passwd ! password
   integer(SHR_KIND_IN),intent(in),optional :: rtpd   ! MSS retention period
   logical,             intent(in),optional :: async  ! true <=> asynchronous put
   logical,             intent(in),optional :: remove ! true <=> rm after put

!EOP

   !----- local -----  
   integer(SHR_KIND_IN)   :: rtpd2    ! MSS retention period
   logical                :: remove2  ! true <=> rm after put
   logical                :: async2   ! true <=> asynchronous put
   character(SHR_KIND_CL) :: passwd2  ! password
   character(SHR_KIND_CL) :: rfn      ! rem_fn without the destination prefix
   character(SHR_KIND_CL) :: cmd      ! command sent to system call

   !----- formats -----
   character(*),parameter :: subName = '(shr_file_put) '
   character(*),parameter :: F00 =   "('(shr_file_put) ',4a)"
   character(*),parameter :: F01 =   "('(shr_file_put) ',a,i3,2a)"
   character(*),parameter :: F02 =   "(a,i4)"

!-------------------------------------------------------------------------------
! Notes: 
! - On some machines the system call will not return a valid error code
! - when things are sent asynchronously, there probably won't be a error code
!   returned.
!-------------------------------------------------------------------------------

   remove2 =.false. ; if ( PRESENT(remove )) remove2 = remove
   async2  =.true.  ; if ( PRESENT(async  )) async2  = async
   passwd2 = " "    ; if ( PRESENT(passwd )) passwd2 = passwd
   rtpd2   = 365    ; if ( PRESENT(rtpd   )) rtpd2   = rtpd
   rcode   = 0

   if ( trim(rem_fn) == trim(loc_fn) ) then
      !------------------------------------------------------
      ! (remote file name) == (local file name) => do nothing
      !------------------------------------------------------
      cmd = 'do nothing: remote file = local file = '//trim(loc_fn)
      rcode = 0
   else if ( index(rem_fn,":") == 0 .or. rem_fn(1:3) == "cp:" ) then
      !------------------------------------------------------
      ! put via unix cp
      !------------------------------------------------------
      rfn = rem_fn
      if ( rem_fn(1:3) == "cp:") rfn = rem_fn(4:len_trim(rem_fn))
#if defined(CATAMOUNT)
      call shr_jlcp(trim(loc_fn),len_trim(loc_fn),trim(rfn),len_trim(rfn),rcode)
      if (remove2) call unlink(trim(loc_fn))
      if (async2) write(6,F00) 'Error: asynchronous copy not supported.'
      cmd = 'shr_jlcp -f '//trim(loc_fn)//' '//trim(rfn)
      rcode = 0
#else
      cmd = '/bin/cp -f '//trim(loc_fn)//' '//trim(rfn)
      if (remove2) cmd = trim(cmd)//' && /bin/rm -f '//trim(loc_fn)
      if (async2 ) cmd = trim(cmd)//' & '
      call shr_sys_system(trim(cmd),rcode)
#endif
   else if ( rem_fn(1:4) == "mss:" ) then
      !------------------------------------------------------
      ! put onto NCAR's MSS
      !------------------------------------------------------
      if (rtpd2 > 9999) rtpd2 = 9999
      write(cmd,F02) '/usr/local/bin/msrcp -period ',rtpd2
      if (async2 .and. (.not. remove2) ) cmd = trim(cmd)//' -async '
      if (len_trim(passwd2) > 0 ) cmd = trim(cmd)//' -wpwd '//trim(passwd)
      cmd = trim(cmd)//' '//trim(loc_fn)//' '//trim(rem_fn)
      if (remove2) cmd = trim(cmd)//' && /bin/rm -f '//trim(loc_fn)
      if (async2 .and. remove2 ) cmd = trim(cmd)//' & '
      call shr_sys_system(trim(cmd),rcode)
   else if ( rem_fn(1:5) == "hpss:") then
      !------------------------------------------------------
      ! put onto LANL's hpss
      !------------------------------------------------------
      rcode = -1
      cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
      write(6,F00) 'ERROR: hpss option not yet implemented'
   else if ( rem_fn(1:5) == "rcp:" ) then
      !------------------------------------------------------
      ! put via rcp
      !------------------------------------------------------
      rcode = -1
      cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
      write(6,F00) 'ERROR: rcp option not yet implemented'
   else if ( rem_fn(1:5) == "null:" ) then
      ! do nothing
      cmd = "null prefix => no file archival, do nothing"
      rcode = 0
   else 
      !------------------------------------------------------
      ! unrecognized remote file location
      !------------------------------------------------------
      rcode = -1
      cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
      write(6,F00) 'ERROR: unrecognized archive device = ',trim(rem_fn)
   end if

   write(6,F01) 'rcode =',rcode,' cmd = ', trim(cmd)

END SUBROUTINE shr_file_put

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_get -- Get a file from archival location.
!
! !DESCRIPTION:
!    a generic, extensible get-local-file-from-archive routine
!
! USAGE:
!    call shr_file_get(rcode,"foo","/home/user/foo")
!    if ( rcode /= 0 ) call shr_sys_abort( "error getting file foo" )
!    call shr_file_get(rcode,"foo","cp:/home/user/foo",remove=.true.)
!    if ( rcode /= 0 ) call shr_sys_abort( "error getting file foo" )
!    call shr_file_get(rcode,"foo","mss:/USER/foo",clobber=.true.)
!    if ( rcode /= 0 ) call shr_sys_abort( "error getting file foo from MSS" )
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_get(rcode,loc_fn,rem_fn,passwd,async,clobber)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(out)         :: rcode   ! return code (non-zero means error)
   character(*)        ,intent(in)          :: loc_fn  ! local filename
   character(*)        ,intent(in)          :: rem_fn  ! remote filename
   character(*)        ,intent(in),optional :: passwd  ! password
   logical             ,intent(in),optional :: async   ! true <=> asynchronous get
   logical             ,intent(in),optional :: clobber ! true <=> clobber existing file

!EOP

   !----- local -----
   logical                :: async2   ! true <=> asynchronous get
   logical                :: clobber2 ! true <=> clobber existing file
   logical                :: exists   ! true <=> local file a ready exists
   character(SHR_KIND_CL) :: passwd2  ! password
   character(SHR_KIND_CL) :: rfn      ! rem_fn without the destination prefix
   character(SHR_KIND_CL) :: cmd      ! command sent to system call

   !----- formats -----
   character(*),parameter :: subName = '(shr_file_get) '
   character(*),parameter :: F00 =   "('(shr_file_get) ',4a)"
   character(*),parameter :: F01 =   "('(shr_file_get) ',a,i3,2a)"

!-------------------------------------------------------------------------------
! Notes: 
! - On some machines the system call will not return a valid error code
! - When things are sent asynchronously, there probably won't be a error code
!   returned.
!-------------------------------------------------------------------------------

   passwd2  = " "     ; if (PRESENT(passwd )) passwd2  = passwd
   async2   = .false. ; if (PRESENT(async  )) async2   = async
   clobber2 = .false. ; if (PRESENT(clobber)) clobber2 = clobber
   rcode    = 0

   inquire(file=trim(loc_fn),exist=exists)

   if ( exists .and. .not. clobber2 ) then
      !------------------------------------------------------
      ! (file exists) and (don't clobber) => do nothing
      !------------------------------------------------------
      cmd = 'do nothing: file exists & no-clobber for '//trim(loc_fn)
      rcode = 0
   else if ( trim(rem_fn) == trim(loc_fn) ) then
      !------------------------------------------------------
      ! (remote file name) == (local file name) => do nothing
      !------------------------------------------------------
      cmd = 'do nothing: remote file = local file for '//trim(loc_fn)
      rcode = 0
   else if (index(rem_fn,":") == 0 .or. rem_fn(1:3) == "cp:") then
      !------------------------------------------------------
      ! get via unix cp
      !------------------------------------------------------
      rfn = rem_fn  ! remove prefix from this temp file name
      if (rem_fn(1:3) == "cp:") rfn = rem_fn(4:len_trim(rem_fn))
#if defined(CATAMOUNT)
      call shr_jlcp(trim(rfn),len(trim(rfn)),trim(loc_fn),len(trim(loc_fn)),rcode)
      if (async2) write(6,F00) 'Error: asynchronous copy not supported.'
      cmd = 'shr_jlcp -f '//trim(rfn)//' '//trim(loc_fn)
      rcode = 0
#else
      cmd = '/bin/cp -f '//trim(rfn)//' '//trim(loc_fn)
      if (async2) cmd = trim(cmd)//' & '
      call shr_sys_system(trim(cmd),rcode)
#endif
   else if ( rem_fn(1:4) == "mss:" ) then
      !------------------------------------------------------
      ! get from NCAR's MSS
      !------------------------------------------------------
      cmd = '/usr/local/bin/msrcp '
      if (async2)   cmd = trim(cmd)//' -async '
      cmd = trim(cmd)//' '//trim(rem_fn)//' '//trim(loc_fn)
      call shr_sys_system(trim(cmd),rcode)
   else if ( rem_fn(1:5) == "hpss:") then
      !------------------------------------------------------
      ! get from LANL's hpss
      !------------------------------------------------------
      rcode = -1
      cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
      write(6,F00) 'ERROR: hpss option not yet implemented'
   else if ( rem_fn(1:5) == "rcp:" ) then
      !------------------------------------------------------
      ! get via rcp
      !------------------------------------------------------
      rcode = -1
      cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
      write(6,F00) 'ERROR: rcp option not yet implemented'
   else if ( rem_fn(1:5) == "null:" ) then
      ! do nothing
      cmd = "null prefix => no file retrieval, do nothing"
      rcode = 0
   else 
      !------------------------------------------------------
      ! unrecognized remote file location
      !------------------------------------------------------
      rcode = -1
      cmd = 'rem_fn='//trim(rem_fn)//' loc_fn='//trim(loc_fn)
      write(6,F00) 'ERROR: unrecognized archive device = ',trim(rem_fn)
   end if

   write(6,F01) 'rcode =',rcode,' cmd = ', trim(cmd)

END SUBROUTINE shr_file_get

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_getUnit -- Get a free FORTRAN unit number
!
! !DESCRIPTION: Get the next free FORTRAN unit number.
!
! !REVISION HISTORY:
!     2005-Dec-14 - E. Kluzek - creation
!
! !INTERFACE: ------------------------------------------------------------------  

INTEGER FUNCTION shr_file_getUnit ( unit )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in),optional :: unit ! desired unit number

!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: n      ! loop index
   logical                :: opened ! If unit opened or not

   !----- formats -----
   character(*),parameter :: subName = '(shr_file_getUnit) '
   character(*),parameter :: F00   = "('(shr_file_getUnit) ',A,I4,A)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if (present (unit)) then
      inquire( unit, opened=opened )
      if (unit < 0 .or. unit > shr_file_nUnits) then
         write(6,F00) 'invalid unit number request:', unit
         call shr_sys_abort( 'ERROR: bad input unit number' )
      else if (opened .or. UnitTag(unit) .or. unit == 0 .or. unit == 5 &
      .or.     unit == 6) then
         write(6,F00) 'unit number ', unit, ' is already in use'
         call shr_sys_abort( 'ERROR: Input unit number already in use' )
      else
         shr_file_getUnit = unit
         UnitTag (unit)   = .true.
         return
      end if

   else
      ! --- Choose first available unit other than 0, 5, or 6 ------
      do n=1,shr_file_nUnits
         inquire( n, opened=opened )
         if (n == 5 .or. n == 6 .or. opened) then
            cycle
         end if
         if ( .not. UnitTag(n) ) then
            shr_file_getUnit = n
            UnitTag(n)       = .true.
            return
         end if
      end do
   end if

   call shr_sys_abort( subName//': Error: no available units found' )

END FUNCTION shr_file_getUnit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_freeUnit -- Free up a FORTRAN unit number
!
! !DESCRIPTION: Free up the given unit number
!
! !REVISION HISTORY:
!     2005-Dec-14 - E. Kluzek - creation
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_freeUnit ( unit)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: unit  ! unit number to be freed

!EOP

   !----- local -----

   !----- formats -----
   character(*), parameter :: subName = '(shr_file_freeUnit) '
   character(*), parameter :: F00 =   "('(shr_file_freeUnit) ',A,I4,A)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if (unit < 0 .or. unit > shr_file_nUnits) then
      write(6,F00) 'invalid unit number request:', unit
   else if (unit == 0 .or. unit == 5 .or. unit == 6) then
      call shr_sys_abort( subName//': Error: units 0, 5, and 6 must not be freed' )
   else if (UnitTag(unit)) then
      UnitTag (unit) = .false.
   else
      write(6,F00) 'unit ', unit, ' was not in use'
   end if

   return

END SUBROUTINE shr_file_freeUnit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_stdio -- Change working directory, and redirect stdin/stdout
!
! !DESCRIPTION:
!   1) change the cwd (current working directory) and 
!   2) redirect stdin & stdout (units 5 & 6) to named files,
!   where the desired cwd & files are specified by namelist file.
!
!   Normally this is done to work around limitations in the execution syntax
!   of common MPI implementations.  For example, SGI's mpirun syntax is not 
!   flexible enough to allow MPMD models to select different execution
!   directories or to redirect stdin & stdout on the command line.  
!   Such functionality is highly desireable for CCSM purposes.  
!   ie. mpirun can't handle this:
!   unix> cd /usr/tmp/jdoe/csm/case01/atm ; atm < atm.parm > atm.log &
!   unix> cd /usr/tmp/jdoe/csm/case01/cpl ; cpl < cpl.parm > cpl.log &
!   etc.
!
! ASSUMPTIONS:
! o if the cwd, stdin, or stdout are to be changed, there must be a namelist
!   file in the cwd named <model>_stdio.nml  where <model> is provided via
!   subroutine dummy argument. 
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_stdio(model)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in) :: model ! used to construct env varible name

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(shr_file_stdio) '
   character(*),parameter :: F00   = "('(shr_file_stdio) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   call shr_file_chdir   (model) ! changes cwd
   call shr_file_chStdOut(model) ! open units 5 & 6 to named files
   call shr_file_chStdIn (model) ! open units 5 & 6 to named files
 
END SUBROUTINE shr_file_stdio

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_chdir -- Change working directory.
!
! !DESCRIPTION:
!   change the cwd (current working directory), see shr_file_stdio for notes
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_chdir(model, rcodeOut)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: model    ! used to construct env varible name
   integer(SHR_KIND_IN),intent(out),optional :: rcodeOut ! Return error code

!EOP

   !--- local ---
   character(SHR_KIND_CL) :: dir      ! directory to cd to
   integer  (SHR_KIND_IN) :: rcode    ! Return error code
   character(SHR_KIND_CL) :: filename ! namelist file to read

   !--- formats ---
   character(*),parameter :: subName = '(shr_file_chdir) '
   character(*),parameter :: F00 =   "('(shr_file_chdir) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   call shr_file_stdioReadNL( model, filename, dirOut=dir, rcodeOut=rcode )
   if (dir /= "nochange") then
      call shr_sys_chdir(dir ,rcode)
      write(*,F00) "read ",trim(filename),", changed cwd to ",trim(dir)
   else
      write(*,F00) "read ",trim(filename),", cwd has *not* been changed"
      rcode = 1
   endif
   if ( present(rcodeOut) ) rcodeOut = rcode
 
END SUBROUTINE shr_file_chdir

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_dirio --- Change stdin and stdout.
!
! !DESCRIPTION:
!   change the stdin & stdout (units 5 & 6), see shr_file_stdio for notes
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_dirio(model)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in) :: model ! used to construct env varible name

!EOP

   !--- local ---

   !--- formats ---
   character(*),parameter :: subName = '(shr_file_dirio) '

!-------------------------------------------------------------------------------
! Notes:
!
!-------------------------------------------------------------------------------

   call shr_file_chStdIn (model)
   call shr_file_chStdOut(model)
 
END SUBROUTINE shr_file_dirio

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_chStdIn -- Change stdin
!
! !DESCRIPTION:
!   change the stdin (unit 5), see shr_file_stdio for notes
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_chStdIn( model, NLFilename, rcodeOut )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)          ,intent(in)           :: model      ! used to construct env var name
   character(SHR_KIND_CL),intent(out),optional :: NLFilename ! open unit 5 to this
   integer  (SHR_KIND_IN),intent(out),optional :: rcodeOut   ! return code

!EOP

   !--- local ---
   character(SHR_KIND_CL) :: stdin    ! open unit 5 to this file
   character(SHR_KIND_CL) :: nlfile   ! Namelist filename for model to read from
   character(SHR_KIND_CL) :: filename ! namelist file to read
   integer  (SHR_KIND_IN) :: rcode    ! return code

   !--- formats ---
   character(*),parameter :: subName = '(shr_file_chStdIn) '
   character(*),parameter :: F00   = "('(shr_file_chStdIn) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   call shr_file_stdioReadNL( model, filename, stdinOut=stdin, &
                              nlfileOut=nlfile, rcodeOut=rcode )
   if (stdin  /= "nochange") then
      open(unit=5,file=stdin ,status='UNKNOWN',iostat=rcode)
      if ( rcode /= 0 )then
         write(*,F00) "read ",trim(filename),', error opening file as unit 5:', &
                      trim(nlfile)
      else
         write(*,F00) "read ",trim(filename),', unit 5 connected to ', &
                      trim(stdin)
      end if
   else
      write(*,F00) "read ",trim(filename),', unit 5 has *not* been redirected'
   endif
   if ( len_trim(nlfile) > 0) then
      write(*,F00) "read ",trim(filename),', read namelist from file:',trim(nlfile)
      if ( .not. present(NLFilename) )then
         write(*,F00) "error: namelist filename NOT present"
         rcode = 7
      end if
   else
      write(*,F00) "read ",trim(filename),", "
      if ( present(NLFilename) )then
         write(*,F00) "error: namelist filename present, but returned value is empty"
         rcode = 8
      end if
   endif
   if ( present(NLFilename) ) NLFilename = nlfile
   if ( present(rcodeOut)   ) rcodeOut   = rcode
 
END SUBROUTINE shr_file_chStdIn

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_stdout -- Change stdout
!
! !DESCRIPTION:
!   change the stdout (unit 6), see shr_file_stdio for notes
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_chStdOut(model,rcodeOut)

   implicit none

! !INPUT/OUTPUT PARAMETERS:
   !--- arguments ---
   character(*)        ,intent(in)           :: model    ! used to construct env varible name
   integer(SHR_KIND_IN),intent(out),optional :: rcodeOut ! Return error code
!EOP

   !--- local ---
   character(SHR_KIND_CL) :: filename ! namelist file to read
   character(SHR_KIND_CL) :: stdout   ! open unit 6 to this file
   integer  (SHR_KIND_IN) :: rcode    ! return code

   !--- formats ---
   character(*),parameter :: subName = '(shr_file_chStdOut) '
   character(*),parameter :: F00   = "('(shr_file_chStdOut) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   call shr_file_stdioReadNL( model, filename, stdoutOut=stdout, &
                              rcodeOut=rcode )
   if (stdout /= "nochange") then
      close(6)
      open(unit=6,file=stdout,position='APPEND')
      write(6,F00) "read ",trim(filename),', unit 6 connected to ',trim(stdout)
      call shr_sys_flush(6)
   else
      write(*,F00) "read ",trim(filename),', unit 6 has *not* been redirected'
      rcode = 1
   endif
 
   if ( present(rcodeOut)   ) rcodeOut = rcode
 
END SUBROUTINE shr_file_chStdOut

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_file_stdioReadNL -- read in stdio namelist
!
! !DESCRIPTION:
!   Read in the stdio namelist for any given model type. Return any of the
! needed input namelist variables as optional arguments. Return "nochange" in
! dir, stdin, or stdout if shouldn't change.
!
! !INTERFACE: ------------------------------------------------------------------  

SUBROUTINE shr_file_stdioReadNL( model, filename, dirOut, stdinOut, stdoutOut, &
                                 NLFileOut, rcodeOut )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)          ,intent(in)           :: model     ! used to construct env varible name
   character(SHR_KIND_CL),intent(out)          :: filename  ! nml file to read from unit 5
   character(SHR_KIND_CL),intent(out),optional :: NLFileOut ! open unit 6 to this file
   character(SHR_KIND_CL),intent(out),optional :: dirOut    ! directory to cd to
   character(SHR_KIND_CL),intent(out),optional :: stdinOut  ! open unit 5 to this file
   character(SHR_KIND_CL),intent(out),optional :: stdoutOut ! open unit 6 to this file
   integer  (SHR_KIND_IN),intent(out),optional :: rcodeOut  ! return code

!EOP

   !--- local ---
   logical                :: exists   ! true iff file exists
   character(SHR_KIND_CL) :: dir      ! directory to cd to
   character(SHR_KIND_CL) :: stdin    ! open unit 5 to this file
   character(SHR_KIND_CL) :: stdout   ! open unit 6 to this file
   character(SHR_KIND_CL) :: NLFile   ! namelist file to read seperately
   integer  (SHR_KIND_IN) :: rcode    ! return code
   integer  (SHR_KIND_IN) :: unit     ! Unit to read from

   namelist / stdio / dir,stdin,stdout,NLFile

   !--- formats ---
   character(*),parameter :: subName = '(shr_file_stdioReadNL) '
   character(*),parameter :: F00   = "('(shr_file_stdioReadNL) ',4a)"
   character(*),parameter :: F01   = "('(shr_file_stdioReadNL) ',2a,i6)"

!-------------------------------------------------------------------------------
! Notes:
!
!-------------------------------------------------------------------------------

   rcode  = 0
   dir    = "nochange"
   stdin  = "nochange"
   stdout = "nochange"
   NLFile = " "

   filename = trim(model)//"_stdio.nml"  ! eg. file="cpl_stdio.nml"
   inquire(file=filename,exist=exists)

   if (.not. exists) then
      write(*,F00) "file ",trim(filename),& 
      & " doesn't exist, can not read stdio namelist from it"
      rcode = 9
   else
      unit = shr_file_getUnit()
      open (unit,file=filename,action="READ")
      read (unit,nml=stdio,iostat=rcode)
      close(unit)
      call shr_file_freeUnit( unit )
      if (rcode /= 0) then
         write(*,F01) 'ERROR: reading ',trim(filename),', iostat=',rcode
         call shr_sys_abort(subName//" ERROR reading "//trim(filename) )
      end if
   endif
   if ( len_trim(NLFile) > 0 .and. trim(stdin) /= "nochange" )then
      write(*,F00) "Error: input namelist:"
      write(*,nml=stdio)
      call shr_sys_abort(subName//" ERROR trying to both redirect AND "// &
                         "open namelist filename" )
   end if
   if ( present(NLFileOut) ) NLFileOut = NLFile
   if ( present(dirOut)    ) dirOut    = dir
   if ( present(stdinOut)  ) stdinOut  = stdin
   if ( present(stdoutOut) ) stdoutOut = stdout
   if ( present(rcodeOut)  ) rcodeOut  = rcode
 
END SUBROUTINE shr_file_stdioReadNL

!===============================================================================
!===============================================================================

END MODULE shr_file_mod
