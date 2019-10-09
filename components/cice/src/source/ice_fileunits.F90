!=======================================================================
!BOP
!
! !MODULE: ice_fileunits
!
! !DESCRIPTION:

!  This module contains an I/O unit manager for tracking, assigning
!  and reserving I/O unit numbers.
!
!  There are three reserved I/O units set as parameters in this
!  module.  The default units for standard input (stdin), standard
!  output (stdout) and standard error (stderr).  These are currently
!  set as units 5,6,6, respectively as that is the most commonly
!  used among vendors. However, the user may change these if those
!  default units are conflicting with other models or if the
!  vendor is using different values.
!
!  The maximum number of I/O units per node is currently set by
!  the parameter ice\_IOMaxUnits.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_fileunits.F90 58 2007-03-29 15:56:53Z eclare $
!
! author: Elizabeth C. Hunke, LANL
! 2006: ECH converted to free source form (F90)
! 2007: ECH added dynamic file units, modified from POP_IOUnitsMod.F90
!
! !INTERFACE:
!
      module ice_fileunits
!
! !USES:
      use ice_kinds_mod
#ifdef CCSMCOUPLED
      use shr_file_mod
#endif
!
!EOP
!=======================================================================

      implicit none
      save

      character (len=char_len) :: &
         diag_type               ! 'stdout' or 'file'

      integer (kind=int_kind) :: &
         nu_grid       , &  ! grid file
         nu_kmt        , &  ! land mask file
         nu_nml        , &  ! namelist input file
         nu_forcing    , &  ! forcing data file
         nu_dump       , &  ! dump file for restarting
         nu_restart    , &  ! restart input file
         nu_dump_aero  , &  ! dump file for restarting ice aerosol tracer MH
         nu_restart_aero,&  ! restart input file for ice aerosol tracer MH
         nu_dump_age   , &  ! dump file for restarting ice age tracer
         nu_restart_age, &  ! restart input file for ice age tracer
         nu_dump_FY    , &  ! dump file for restarting FY ice tracer
         nu_restart_FY , &  ! restart input file for FY ice tracer
         nu_dump_lvl   , &  ! dump file for restarting level ice tracer
         nu_restart_lvl, &  ! restart input file for level ice tracer
         nu_dump_pond  , &  ! dump file for restarting melt pond tracer
         nu_restart_pond,&  ! restart input file for melt pond tracer
         nu_rst_pointer, &  ! pointer to latest restart file
         nu_history    , &  ! binary history output file
         nu_hdr        , &  ! header file for binary history output
         nu_diag       , &  ! diagnostics output file
         nu_timing          ! timing output file

      integer (kind=int_kind) :: &
         diag_level         ! per-processor diagnostics level

      character (11) :: &
         nml_filename = 'ice_in' ! namelist input file name

      integer (kind=int_kind), parameter :: &
         ice_stdin  =  5, & ! reserved unit for standard input
         ice_stdout =  6, & ! reserved unit for standard output
         ice_stderr =  6    ! reserved unit for standard error


!EOP
!BOC
      integer (kind=int_kind), parameter :: &
         ice_IOUnitsMinUnits = 11,   & ! do not use unit numbers below this
         ice_IOUnitsMaxUnits = 99      ! maximum number of open units

      logical (kind=log_kind), dimension(ice_IOUnitsMaxUnits) :: &
         ice_IOUnitsInUse   ! flag=.true. if unit currently open

#ifdef CCSMCOUPLED
      ! instance control
      integer (kind=int_kind), public :: inst_index
      character(len=16)      , public :: inst_name
      character(len=16)      , public :: inst_suffix
#endif

!EOC
!=======================================================================

contains

!=======================================================================
!BOP
! !IROUTINE: init_fileunits
! !INTERFACE:

 subroutine init_fileunits

! !DESCRIPTION:
!  This routine grabs needed unit numbers. 
!  nu_diag is set to 6 (stdout) but may be reset later by the namelist. 
!  nu_nml is obtained separately.

         nu_diag = ice_stdout  ! default

#ifdef CCSMCOUPLED
         ice_IOUnitsInUse = .false.
         ice_IOUnitsInUse(ice_stdin)  = .true. ! reserve unit 5
         ice_IOUnitsInUse(ice_stdout) = .true. ! reserve unit 6
         ice_IOUnitsInUse(ice_stderr) = .true.

         nu_grid        = shr_file_getUnit()
         nu_kmt         = shr_file_getUnit()
         nu_forcing     = shr_file_getUnit()
         nu_dump        = shr_file_getUnit()
         nu_restart     = shr_file_getUnit()
         nu_dump_aero   = shr_file_getUnit() 
         nu_restart_aero= shr_file_getUnit()
         nu_dump_age    = shr_file_getUnit()
         nu_restart_age = shr_file_getUnit()
         nu_dump_FY     = shr_file_getUnit()
         nu_restart_FY  = shr_file_getUnit()
         nu_dump_lvl    = shr_file_getUnit()
         nu_restart_lvl = shr_file_getUnit()
         nu_dump_pond   = shr_file_getUnit()
         nu_restart_pond = shr_file_getUnit()
         nu_rst_pointer = shr_file_getUnit()
         nu_history     = shr_file_getUnit()
         nu_hdr         = shr_file_getUnit()
         nu_timing         = shr_file_getUnit()
#else
         call get_fileunit(nu_grid)
         call get_fileunit(nu_kmt)
         call get_fileunit(nu_forcing)
         call get_fileunit(nu_dump)
         call get_fileunit(nu_restart)
         call get_fileunit(nu_dump_aero)
         call get_fileunit(nu_restart_aero)
         call get_fileunit(nu_dump_age)
         call get_fileunit(nu_restart_age)
         call get_fileunit(nu_dump_FY)
         call get_fileunit(nu_restart_FY)
         call get_fileunit(nu_dump_lvl)
         call get_fileunit(nu_restart_lvl)
         call get_fileunit(nu_dump_pond)
         call get_fileunit(nu_restart_pond)
         call get_fileunit(nu_rst_pointer)
         call get_fileunit(nu_history)
         call get_fileunit(nu_hdr)
         call get_fileunit(nu_timing)
#endif

 end subroutine init_fileunits

!=======================================================================
!BOP
! !IROUTINE: get_fileunit
! !INTERFACE:

 subroutine get_fileunit(iunit)

! !DESCRIPTION:
!  This routine returns the next available I/O unit and marks it as
!  in use to prevent any later use.
!  Note that {\em all} processors must call this routine even if only
!  the master task is doing the I/O.  This is necessary insure that
!  the units remain synchronized for other parallel I/O functions.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (kind=int_kind), intent(out) :: &
      iunit                     ! next free I/O unit
!EOP
!BOC

   integer (kind=int_kind) :: n  ! dummy loop index

   logical (kind=log_kind) :: alreadyInUse

#ifdef CCSMCOUPLED
   iunit = shr_file_getUnit()
#else
   srch_units: do n=ice_IOUnitsMinUnits, ice_IOUnitsMaxUnits
      if (.not. ice_IOUnitsInUse(n)) then   ! I found one, I found one

         !*** make sure not in use by library or calling routines
         INQUIRE (unit=n,OPENED=alreadyInUse)

         if (.not. alreadyInUse) then
            iunit = n        ! return the free unit number
            ice_IOUnitsInUse(n) = .true.  ! mark iunit as being in use
            exit srch_units
         else
            !*** if inquire shows this unit in use, mark it as
            !***    in use to prevent further queries
            ice_IOUnitsInUse(iunit) = .true.
         endif
      endif
   end do srch_units

   if (iunit > ice_IOUnitsMaxUnits) stop 'ice_IOUnitsGet: No free units'
#endif

!EOC
 end subroutine get_fileunit

!=======================================================================
!BOP
! !IROUTINE: release_all_fileunits
! !INTERFACE:

 subroutine release_all_fileunits

! !DESCRIPTION:
!  This routine releases unit numbers at the end of a run. 

      call release_fileunit(nu_grid)
      call release_fileunit(nu_kmt)
      call release_fileunit(nu_forcing)
      call release_fileunit(nu_dump)
      call release_fileunit(nu_restart)
      call release_fileunit(nu_dump_aero)
      call release_fileunit(nu_restart_aero)
      call release_fileunit(nu_dump_age)
      call release_fileunit(nu_restart_age)
      call release_fileunit(nu_dump_FY)
      call release_fileunit(nu_restart_FY)
      call release_fileunit(nu_dump_lvl)
      call release_fileunit(nu_restart_lvl)
      call release_fileunit(nu_dump_pond)
      call release_fileunit(nu_restart_pond)
      call release_fileunit(nu_rst_pointer)
      call release_fileunit(nu_history)
      call release_fileunit(nu_hdr)
      call release_fileunit(nu_timing)
      if (nu_diag /= ice_stdout) call release_fileunit(nu_diag)

 end subroutine release_all_fileunits

!=======================================================================
!BOP
! !IROUTINE: release_fileunit
! !INTERFACE:

 subroutine release_fileunit(iunit)

! !DESCRIPTION:
!  This routine releases an I/O unit (marks it as available).
!  Note that {\em all} processors must call this routine even if only
!  the master task is doing the I/O.  This is necessary insure that
!  the units remain synchronized for other parallel I/O functions.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETER:

   integer (kind=int_kind), intent(in) :: &
      iunit                    ! I/O unit to be released
!EOP
!BOC

#ifdef CCSMCOUPLED
   call shr_file_freeUnit(iunit)
#else
!  check for proper unit number
   if (iunit < 1 .or. iunit > ice_IOUnitsMaxUnits) then
      stop 'release_fileunit: bad unit'
   endif

!  mark the unit as not in use
   ice_IOUnitsInUse(iunit) = .false.  !  that was easy...
#endif

!EOC
 end subroutine release_fileunit

!=======================================================================
!BOP
! !IROUTINE: flush_fileunit
! !INTERFACE:

 subroutine flush_fileunit(iunit)

! !DESCRIPTION:
!  This routine enables a user to flush the output from an IO unit
!  (typically stdout) to force output when the system is buffering
!  such output.  Because this system function is system dependent,
!  we only support this wrapper and users are welcome to insert the
!  code relevant to their local machine.  In the case where the CCSM
!  libraries are available, the shared routine for sys flush can be
!  used (and is provided here under a preprocessor option).
!
! !REVISION HISTORY:
!  same as module
!
! !USES:

#ifdef CCSMCOUPLED
      use shr_sys_mod, only : shr_sys_flush
#endif

! !INPUT PARAMETER:

   integer (kind=int_kind), intent(in) :: &
      iunit                    ! I/O unit to be flushed
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  insert your system code here
!
!-----------------------------------------------------------------------

#ifdef CCSMCOUPLED
   call shr_sys_flush(iunit)
#else
#if (defined IRIX64 || defined CRAY || defined OSF1 || defined SUNOS || defined LINUX || defined NEC_SX | defined UNICOSMP)
   call flush(iunit)
#endif
#if (defined AIX)
   call flush_(iunit)
#endif
#endif

!EOC
 end subroutine flush_fileunit

!=======================================================================

      end module ice_fileunits

!=======================================================================
