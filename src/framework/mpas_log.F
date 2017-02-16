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
             mpas_log_instance_create, &
             mpas_log_finalize, &
             mpas_log_write

   !--------------------------------------------------------------------
   ! Private module variables
   !--------------------------------------------------------------------


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

   subroutine mpas_log_init(coreLogInfo, domain, err, logUnit)

      use mpas_io_units

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (domain_type), pointer, intent(in) :: domain !< Input: domain information
      integer, intent(in), optional :: logUnit !< Input - Optional: Unit to set log unit to

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
      integer :: err_tmp

      err = 0

      ! Create log instance
      call mpas_log_instance_create(coreLogInfo, err_tmp)
      err = ior(err, err_tmp)

      ! Point the module-level log instance to the core's instance
      ! (From this point onward, only use the module level instance)
      mpas_log_info => coreLogInfo

      ! Store the name of this core (makes it available for inserting into messages or filenames)
      mpas_log_info % coreName = domain % core % coreName

      ! Store the task ID (eliminates the need for dminfo later)
      !   This will be used for 1) checking if on master task, 2) inserting task ID into the name of the abort file
      mpas_log_info % taskID = domain % dminfo % my_proc_id

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

      ! Open log file if it is active
      if (mpas_log_info % logActive) then

         ! Generate the log file name
         ! TODO: routine to expand procid with various numbers of digits.  Unify with usage in global_abort routine.
         write(mpas_log_info % logFileName, fmt='(a, i4.4, a)') "stick.", mpas_log_info % taskID, ".out"

         ! Determine the unit number
         if(present(logUnit)) then  ! TODO If unit number provided, that should assume the file has already been opened!
            mpas_log_info % logUnit = logUnit
         else
            call mpas_new_unit(mpas_log_info % logUnit)
         end if

         ! Open the file
         open (unit = mpas_log_info % logUnit, file = mpas_log_info % logFileName, &
            action="WRITE", status="REPLACE", IOSTAT=err_tmp)
         if ( err_tmp /= 0 ) then
            call mpas_dmpar_global_abort('ERROR: Opening of log file failed for filename: ' // mpas_log_info % logFileName)
         end if
         err = ior(err, err_tmp)  ! This is unnecessary because we just die with an error

         ! -- Test message --
         write(mpas_log_info % logUnit, '(a)') '------------------------------------------------------'
         write(mpas_log_info % logUnit, '(a,a,a,i5,a,i5)') 'Beginning MPAS-', trim(mpas_log_info % coreName), ' log file for task ', mpas_log_info % taskID, ' of ', domain % dminfo % nprocs
         write(mpas_log_info % logUnit, '(a)') '------------------------------------------------------'
         write(mpas_log_info % logUnit, '(a)') ''
         !write(mpas_log_info % logUnit,*) "Hello world!"
         flush(mpas_log_info % logUnit)

      endif


   !--------------------------------------------------------------------
   end subroutine mpas_log_init



!***********************************************************************
!
!  routine mpas_log_instance_create
!
!> \brief   Creates a log manager instance
!> \author  Matt Hoffman
!> \date    14 February 2017
!> \details
!>  This routine creates an instance of the log manager.
!>  It allocates an mpas_log_type instance and assigns initial values to
!>  all of its members.  The initial values are intended to be invalid
!>  for proper usage.  They should be subsequently populated by mpas_log_init.
!
!-----------------------------------------------------------------------

   subroutine mpas_log_instance_create(logInfo, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (MPAS_log_type), pointer, intent(inout) :: logInfo  !< log manager instance to create

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------


      err = 0

      allocate(logInfo)

      ! Initialize log info
      logInfo % logUnit = -1  ! initialize to invalid unit number
      logInfo % logActive = .false.
      logInfo % errUnit = -1  ! initialize to invalid unit number
      logInfo % errUnitOpen = .false.

      logInfo % taskID = -1 ! initialize to invalid task ID

      ! initialize to generic name that won't conflict with old logging method
      logInfo % logFileName = 'mpaslog.out'
      logInfo % errFileName = 'mpaslog.err'

      ! Initialize counters
      logInfo % standardMessageCount = 0
      logInfo % warningMessageCount = 0
      logInfo % errorMessageCount = 0
      logInfo % criticalErrorMessageCount = 0

   !--------------------------------------------------------------------
   end subroutine mpas_log_instance_create



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
      real(kind=RKIND), dimension(:), intent(in), optional :: expArgs  !< Input: exponential notation variable values to insert into message
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


      err = 0


      ! Construct message by expanding variable values as needed and inserting message type prefix

      ! Write message to standard log file
      !   Accounting for debug/optimized settings
      !   Accounting for optional head-node-only argument
      write(mpas_log_info % logUnit, *) message

      ! Optionally flush the buffer

      ! Increment appropriate message counter

      ! Handle error/critical error logic if necessary


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
!>  This routine finalizes the log manager. Details...
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


      err = 0

      !Calculate and print out statistics of how many of each message type were printed (which will be stored as module level private variables)

      !Close the log file

   !--------------------------------------------------------------------
   end subroutine mpas_log_finalize


end module mpas_log
