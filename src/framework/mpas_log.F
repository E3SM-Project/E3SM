! Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_log
!
!> \MPAS message logging manager
!> \author Matt Hoffman
!> \date   14 February 2017
!> \details
!>  This module contains the routines for managing the writing of messages
!>  to log files.
!>  The log module operates around a variable named mpas_log_info that contains
!>  all of the external information needed to perform logging operations
!>  for the current model configuration.
!>  mpas_log_info is a module level pointer to the active core's mpas_log_type
!>  instance.  This allows the core's domain instance to "own" the log manager,
!>  while the log module has internal access to it.  This approach has two
!>  benefits:
!>  1. Calls to mpas_write_log do not require any configuration metadata
!>  (e.g., a unit number) to be passed in, simplifying the call API.
!>  2. Because the log module uses a pointer to a mpas_log_type instance,
!>  different instances can be swapped in and out from higher levels of
!>  the modeling framework during execution.  This is necessary to allow
!>  multiple different MPAS cores to operate in the same climate model,
!>  but would also apply to multiple MPAS cores running in a single MPAS
!>  executable (if this is ever supported), or multiple domain instances
!>  of the same MPAS core running together.  (It is the responsibility of
!>  the higher-level driver code(s) to manage any swapping of mpas_log_type
!>  instances.)
!
!-----------------------------------------------------------------------


module mpas_log

   use mpas_derived_types
   use mpas_abort, only : mpas_dmpar_global_abort

   implicit none
   private

   !--------------------------------------------------------------------
   ! Public module variables
   !--------------------------------------------------------------------
   ! (Public parameters found in mpas_log_types.inc)

   type (mpas_log_type), pointer :: mpas_log_info => null() !< derived type holding the metadata for a logging manager instance
      !< Each core needs to set this up.  Note that this pointer can be externally
      !< reassigned if the logging instance needs to be swapped
      !< (e.g., when multiple logging instances are running in an ESM).

   !--------------------------------------------------------------------
   ! Public member functions
   !--------------------------------------------------------------------
   public :: mpas_log_init, &
             mpas_log_open, &
             mpas_log_finalize, &
             mpas_log_write

   !--------------------------------------------------------------------
   ! Private module variables
   !--------------------------------------------------------------------
   character(len=16), parameter, private :: messagePrefixWarning = "WARNING: "
   character(len=16), parameter, private :: messagePrefixError = "ERROR: "
   character(len=16), parameter, private :: messagePrefixCritical = "CRITICAL ERROR: "


!***********************************************************************


   contains


!***********************************************************************
!
!  routine mpas_log_init
!
!> \brief   Initializes log manager
!> \author  Matt Hoffman
!> \date    14 February 2017
!> \details
!>  This routine initializes the log manager for the active core on each task
!>  by assigning appropriate values to the members of mpas_log_info, based
!>  on the model configuration.  mpas_log_info is a module level pointer
!>  of the active core's mpas_log_type instance.  This allows the core
!>  domain instance to "own" the log manager, while the log module has
!>  internal access to it
!>  It also opens the log file(s) and sends a "Hello World"-like message.
!
!-----------------------------------------------------------------------

   subroutine mpas_log_init(coreLogInfo, domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (domain_type), pointer, intent(in) :: domain !< Input: domain information

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (mpas_log_type), pointer, intent(inout) :: coreLogInfo !< Input/Output: log manager info (already set up by the core)


      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------


      err = 0

      ! Create log instance
      allocate(coreLogInfo)

      ! ------
      ! Initialize all log info members
      ! ------
      coreLogInfo % logUnit = -1  ! initialize to invalid unit number
      coreLogInfo % logActive = .false.
      coreLogInfo % openedByLogModule = .false.
      coreLogInfo % errUnit = -1  ! initialize to invalid unit number
      coreLogInfo % errUnitOpen = .false.

      coreLogInfo % taskID = -1 ! initialize to invalid task ID
      coreLogInfo % nTasks = -1 ! initialize to invalid number of tasks

      ! initialize to generic name that won't conflict with old logging method
      coreLogInfo % logFileName = 'mpaslog.out'
      coreLogInfo % errFileName = 'mpaslog.err'

      ! Initialize counters
      coreLogInfo % outputMessageCount = 0
      coreLogInfo % warningMessageCount = 0
      coreLogInfo % errorMessageCount = 0
      coreLogInfo % criticalErrorMessageCount = 0


      ! ------
      ! Point the module-level log instance to the core's instance
      ! (From this point onward, only use the module level instance)
      ! ------
      mpas_log_info => coreLogInfo

      ! ------
      ! Update values for some members based on model configuration
      !   Note: these choices can be replaced by the core/driver after this routine returns
      ! ------
      ! Store the name of this core (makes it available for inserting into messages or filenames)
      mpas_log_info % coreName = domain % core % coreName

      ! Store the task ID (eliminates the need for dminfo later)
      !   This will be used for 1) checking if on master task, 2) inserting task ID into the name of the abort file
      mpas_log_info % taskID = domain % dminfo % my_proc_id
      mpas_log_info % nTasks = domain % dminfo % nprocs

      ! Set log file to be active or not based on master/nonmaster task and optimized/debug build
      !   * Optimized build: Only master task log is active
      !   * Debug build: All tasks active
#ifdef MPAS_DEBUG
      ! in debug mode, all tasks have active logs
      mpas_log_info % logActive = .true.
#else
      ! in non-debug mode, only task 0 has active log
      if (mpas_log_info % taskID == 0) then
         mpas_log_info % logActive = .true.
      else
         mpas_log_info % logActive = .false.
      endif
#endif

      ! Generate the log file name
      !    Note: If log is inactive, this won't actually be used, but it is always set.
      ! TODO: add routine to expand procid with various numbers of digits.  Unify with usage in global_abort routine.
      write(mpas_log_info % logFileName, fmt='(a, i4.4, a)') "stick.", mpas_log_info % taskID, ".out"

      ! Note: unit number is not set until file is opened

      ! Set up err file name now
      write(mpas_log_info % errFileName, fmt='(a, i4.4, a)') "stick.", mpas_log_info % taskID, ".err"

   !--------------------------------------------------------------------
   end subroutine mpas_log_init



!***********************************************************************
!
!  routine mpas_log_open
!
!> \brief   Connects to a log file for the log module to use
!> \author  Matt Hoffman
!> \date    16 February 2017
!> \details
!>  This routine connects to a log file for the log module to use.
!>  By default, a new file with a new file unit is opened.
!>  Optionally, a fileUnit can be passed in that must already be attached
!>  to an appropriate file.
!>  By default, an output log file is opened, but optionally this can
!>  be an error file.
!>  This routine assumes that the log manager has been initialized and all
!>  log settings are set correctly, including filename, active status, etc.
!>  It uses the module-level pointer instance of the log manager.
!
!-----------------------------------------------------------------------

   subroutine mpas_log_open(logUnit, openErrorFile, err)

      use mpas_io_units

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      integer, intent(in), optional :: logUnit !< Input - Optional: Unit (already open) to set log unit to
      logical, intent(in), optional :: openErrorFile !< Input - Optional: flag to indicate this task's error file should be opened
         !< Default = false, meaning open output log file

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out), optional :: err !< Output - Optional: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer :: err_tmp
      logical :: thisIsOutputLog
      logical :: activeFile
      character(len=StrKind) :: fileName
      integer :: unitNumber
      character (len=6) :: logTypeString
      ! Variables used to inquire about info on a pre-existing file unit (if supplied)
      logical :: isOpen
      character(len=16) :: writeSetting, sequentialSetting, formattedSetting ! length just needs to store 'YES'
      character(len=8) :: date
      character(len=10) :: time

      err = 0

      thisIsOutputLog = .true.
      if (present(openErrorFile)) then
         if (openErrorFile) then
            thisIsOutputLog = .false.
            activeFile = .true.  ! A request to open the error file is always honored.
         endif
      endif

      if (thisIsOutputLog) then
         if (mpas_log_info % logActive) then
            activeFile = .true.
         else
            activeFile = .false.
         endif
      endif

      if (thisIsOutputLog) then
         logTypeString = "Output"
      else
         logTypeString = "Error"
      endif
      ! Open file if it is set to be active
      if (activeFile) then

         ! If a pre-opened unit is provided, just use that
         if(present(logUnit)) then

            unitNumber = logUnit  ! Set local variable used below

            ! Save unit number provided to appropriate variable
            if (thisIsOutputLog) then
               mpas_log_info % logUnit = unitNumber
            else
               mpas_log_info % errUnit = unitNumber
            endif

            ! Retrieve metadata for this unit to get the file name, and to check it is an approriate unit for a log
            inquire(unitNumber, name=fileName, opened=isOpen, write=writeSetting, &
               sequential=sequentialSetting, formatted=formattedSetting, iostat=err_tmp)

            if (err_tmp > 0) then
               call mpas_dmpar_global_abort('ERROR: Log initialization failed because supplied log unit could not be inquired!')
            endif
            if (.not. isOpen) then
               call mpas_dmpar_global_abort('ERROR: Log initialization failed because supplied log unit is not attached'//&
                  ' to an open file!')
            endif
            if (trim(writeSetting) /= 'YES') then
               call mpas_dmpar_global_abort('ERROR: Log initialization failed because supplied log unit is not attached'//&
                 ' to a file with write settings!')
            endif
            if (trim(sequentialSetting) /= 'YES') then
               call mpas_dmpar_global_abort('ERROR: Log initialization failed because supplied log unit is not attached'//&
                 ' to a file with sequential access!')
            endif
            if (trim(formattedSetting) /= 'YES') then
               call mpas_dmpar_global_abort('ERROR: Log initialization failed because supplied log unit is not attached'//&
                 ' to a file connected for formatted output!')
            endif

            ! Save file name we retrieved to appropriate field
            if (thisIsOutputLog) then
               mpas_log_info % logFileName = fileName
            else
               mpas_log_info % errFileName = fileName
            endif


         else  ! If not, we need to set the unit number and open it ourselves

            ! Determine the unit number
            call mpas_new_unit(unitNumber)
            if (thisIsOutputLog) then
               mpas_log_info % logUnit = unitNumber ! save unit
               fileName = mpas_log_info % logFileName ! get desired filename
            else
               mpas_log_info % errUnit = unitNumber ! save unit
               fileName = mpas_log_info % errFileName ! get desired filenamee
            endif

            ! Open the file
            open (unit = unitNumber, file = fileName, action="WRITE", status="REPLACE", IOSTAT=err_tmp)
            if ( err_tmp /= 0 ) then
               call mpas_dmpar_global_abort('ERROR: Opening of log file failed for filename: ' // mpas_log_info % logFileName)
            endif
            err = ior(err, err_tmp)  ! This is unnecessary because we just die with an error
            ! TODO: revisit if we can relax the decision to die immediately here.

            mpas_log_info % openedByLogModule = .true.  ! indicate that we opened it

         endif ! if logUnit was provided

         ! -- Write a header message --
         !write(unitNumber,*) "Hello world!"
         write(unitNumber, '(a)') '----------------------------------------------------------------------'
         write(unitNumber, '(a,a,a,a,a,i7.1,a,i7.1)') 'Beginning MPAS-', trim(mpas_log_info % coreName), ' ', &
            trim(logTypeString), ' Log File for task ', mpas_log_info % taskID, ' of ', mpas_log_info % nTasks
         call date_and_time(date,time)
         write(unitNumber, '(a)') '    Opened at ' // date(1:4)//'/'//date(5:6)//'/'//date(7:8) // &
            ' ' // time(1:2)//':'//time(3:4)//':'//time(5:6)
         write(unitNumber, '(a)') '----------------------------------------------------------------------'
         write(unitNumber, '(a)') ''
         flush(unitNumber)  ! flush header message immediately to eliminate potential confusion

      endif  ! if activeFile

   !--------------------------------------------------------------------
   end subroutine mpas_log_open



!***********************************************************************
!
!  routine mpas_log_write
!
!> \brief   Writes a message to the log file
!> \author  Matt Hoffman
!> \date    14 February 2017
!> \details
!>  This routine writes a message to the log file.  Details...
!
!-----------------------------------------------------------------------

   subroutine mpas_log_write(message, messageType, masterOnly, flushNow, &
                 intArgs, realArgs, expArgs, logicArgs, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      character (len=*), intent(in) :: message  !< Input: message to be printed
      integer, intent(in), optional :: messageType !< Input: message "type"
           !<  output (default), warning, error, critical error
           !< (integer values defined by public module parameters above)
      logical, intent(in), optional :: masterOnly  !< Input: flag to only print message on master task
      logical, intent(in), optional :: flushNow  !< Input: flag to force a flush of the message buffer
      integer, dimension(:), intent(in), optional :: intArgs  !< Input: integer variable values to insert into message
      real(kind=RKIND), dimension(:), intent(in), optional :: realArgs  !< Input: real variable values to insert into message
      real(kind=RKIND), dimension(:), intent(in), optional :: expArgs
         !< Input: exponential notation variable values to insert into message
      logical, dimension(:), intent(in), optional :: logicArgs  !< Input: logical variable values to insert into message

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out), optional :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer :: ierr, ierr_tmp
      logical :: masterOnlyWrite  !< local version of masterOnly
      integer :: messageTypeHere !< local version of messageType
      character(len=16) :: messagePrefix  !< the "keyword" prepended before some message types

      ierr = 0

      if (present(messageType)) then
         messageTypeHere = messageType
      else
         messageTypeHere = MPAS_LOG_OUT  ! Default is 'output' message
      endif

      if (present(masterOnly)) then
         masterOnlyWrite = masterOnly
      else
         masterOnlyWrite = .false.
      endif
      ! Critical errors do not respect the masterOnly flag!
      if ((masterOnlyWrite) .and. (messageTypeHere == MPAS_LOG_CRIT)) then
         masterOnlyWrite = .false.
      endif


      ! Construct message by expanding variable values as needed and inserting message type prefix
      ! TODO

      ! Determine message prefix
      select case (messageTypeHere)
      case (MPAS_LOG_OUT)
         messagePrefix = ""  ! (not used)
      case (MPAS_LOG_WARN)
         messagePrefix = messagePrefixWarning
      case (MPAS_LOG_ERR)
         messagePrefix = messagePrefixError
      case (MPAS_LOG_CRIT)
         messagePrefix = messagePrefixCritical
      case default
         ! TODO handle this error?
      end select

      ! Write message to standard log file
      !   Accounting for debug/optimized settings (logActive)
      !   Accounting for optional head-node-only argument
      if (mpas_log_info % logActive) then
         if ( ((masterOnlyWrite) .and. (mpas_log_info % taskID == 0)) .or. (.not. masterOnlyWrite)) then

            ! --- Actually write the message here! ---
            !    We have to treat MPAS_LOG_OUT messages separately from others because that type
            !    does not have a prefix and we want to add a space after the trimmed prefix for other types.
            if (messageTypeHere == MPAS_LOG_OUT) then
               write(mpas_log_info % logUnit, '(a)') trim(message)
            else
               write(mpas_log_info % logUnit, '(a,a,a)') trim(messagePrefix), ' ', trim(message)
            endif

            ! Optionally flush the buffer (only attempt if the file is active)
            if (present(flushNow)) then
               if (flushNow) then
                  flush(mpas_log_info % logUnit)
               endif
            endif

            ! Increment appropriate message counter
            !   (counters only apply to the output log)
            select case (messageTypeHere)
            case (MPAS_LOG_OUT)
               mpas_log_info % outputMessageCount = mpas_log_info % outputMessageCount + 1
            case (MPAS_LOG_WARN)
               mpas_log_info % warningMessageCount = mpas_log_info % warningMessageCount + 1
            case (MPAS_LOG_ERR)
               mpas_log_info % errorMessageCount = mpas_log_info % errorMessageCount + 1
            case (MPAS_LOG_CRIT)
               mpas_log_info % criticalErrorMessageCount = mpas_log_info % criticalErrorMessageCount + 1
            case default
               ! TODO handle this error? (would have encountered it above...)
            end select


         endif
      endif


      ! Handle error/critical error logic if necessary
      if ((messageTypeHere == MPAS_LOG_ERR) .or. (messageTypeHere == MPAS_LOG_CRIT)) then
         ! Error messages also need to respect the masterOnly option (but critical errors do not - handled above)
         if ( ((masterOnlyWrite) .and. (mpas_log_info % taskID == 0)) .or. (.not. masterOnlyWrite)) then

            ! if err file is not already open, we need to open it
            if (.not. mpas_log_info % errUnitOpen) then
               call mpas_log_open(openErrorFile = .true., err = ierr_tmp)
               ierr = ior(ierr, ierr_tmp)
               mpas_log_info % errUnitOpen = .true.
            endif

            ! Now we need to repeat the message to the err file
            write(mpas_log_info % errUnit, '(a,a,a)') trim(messagePrefix), ' ', trim(message)

            ! If this is a critical error, we need to kill the model
            if (messageTypeHere == MPAS_LOG_CRIT) then
               ! Now abort
               call log_abort()
            endif

         endif
      endif


      ! Pass back err code if it was requested
      if (present(err)) then
         err = ierr
      endif

   !--------------------------------------------------------------------
   end subroutine mpas_log_write



!***********************************************************************
!
!  routine mpas_log_finalize
!
!> \brief   Finalizes the log manager
!> \author  Matt Hoffman
!> \date    14 February 2017
!> \details
!>  This routine finalizes the log manager.
!>  It prints out logging statistics and closes the log file.
!
!-----------------------------------------------------------------------

   subroutine mpas_log_finalize(err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      character(len=strKind) :: msg !< variable to build custom messages
      integer :: outputCount
      character(len=8) :: date
      character(len=10) :: time

      err = 0

      outputCount = mpas_log_info % outputMessageCount  ! save the count before we write the header lines

      ! Print out statistics of how many of each message type were printed
      !   (use mpas_log_write, rather than writing directly to file, so that logActive is applied correctly)
      call mpas_log_write('')
      call mpas_log_write('-----------------------------------------')
      call mpas_log_write('Total log messages printed:')
      write(msg, fmt='(a, i12.1)') '   Output messages =         ', outputCount
      call mpas_log_write(msg)
      write(msg, fmt='(a, i12.1)') '   Warning messages =        ', mpas_log_info % warningMessageCount
      call mpas_log_write(msg)
      write(msg, fmt='(a, i12.1)') '   Error messages =          ', mpas_log_info % errorMessageCount
      call mpas_log_write(msg)
      write(msg, fmt='(a, i12.1)') '   Critical error messages = ', mpas_log_info % criticalErrorMessageCount
      call mpas_log_write(msg)
      call mpas_log_write('-----------------------------------------')
      call date_and_time(date,time)
      write(msg, '(a)') 'Logging complete.  Closing file at ' // date(1:4)//'/'//date(5:6)//'/'//date(7:8) // &
         ' ' // time(1:2)//':'//time(3:4)//':'//time(5:6)
      call mpas_log_write(msg)

      ! Close the log file
      !   Only do this if:
      !   1) the log is active AND
      !   2) the log mgr opened the file (otherwise the driver that opened it should close it)
      if (mpas_log_info % logActive .and. mpas_log_info % openedByLogModule) then
         close(mpas_log_info % logUnit, iostat = err)
      endif

   !--------------------------------------------------------------------
   end subroutine mpas_log_finalize


!***********************************************************************
! PRIVATE ROUTINES
!***********************************************************************


!***********************************************************************
!
!  routine log_abort
!
!> \brief   Kills the model
!> \author  Matt Hoffman
!> \date    14 February 2017
!> \details
!>  This routine kills the model.
!
!-----------------------------------------------------------------------

   subroutine log_abort()

#ifdef _MPI
#ifndef NOMPIMOD
      use mpi
#endif
#endif

      implicit none

#ifdef _MPI
#ifdef NOMPIMOD
      include 'mpif.h'
#endif
#endif

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
#ifdef _MPI
      integer :: mpi_ierr, mpi_errcode, err
#endif
      character(len=8) :: date
      character(len=10) :: time


      ! Close the output log since we can (and provides stats)
      call mpas_log_finalize(err)

      ! Write final message to err file
      call date_and_time(date,time)
      write(mpas_log_info % errUnit, '(a)') 'Logging complete.  Closing file at ' // &
         date(1:4)//'/'//date(5:6)//'/'//date(7:8) // ' ' // time(1:2)//':'//time(3:4)//':'//time(5:6)

      ! Close the err log to be clean
      close(mpas_log_info % errUnit)

#ifdef _MPI
      call MPI_Abort(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
#else
      stop
#endif

   !--------------------------------------------------------------------
   end subroutine log_abort


end module mpas_log
