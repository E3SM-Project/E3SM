!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_ConfigMod

!BOP
!
! !MODULE:  POP_ConfigMod
!
! !DESCRIPTION:
!  This module contains routines for reading input configuration
!  data from a configuration file.  Variables for a specified module 
!  are read from an input file and broadcast to all processors.
!  A default value for a variable can be specified and will be
!  used if the variable is not found in the input file.
!
! !REVISION HISTORY:
!  SVN:$Id: POP_ConfigMod.F90 15 2006-08-21 20:04:13Z  $
!  2007-12-17: Phil Jones
!     initial version
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_IOUnitsMod
   use POP_CommMod
   use POP_BroadcastMod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_ConfigOpen,               &
             POP_ConfigClose,              &
             POP_ConfigRead

! !PUBLIC DATA MEMBERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  private config module variables
!
!-----------------------------------------------------------------------

   character (POP_charLength) :: &
      configFileDefault = 'pop2_in'

   interface POP_ConfigRead
      module procedure POP_ConfigReadI4,      &
                       POP_ConfigReadR4,      &
                       POP_ConfigReadR8,      &
                       POP_ConfigReadLogical, &
                       POP_ConfigReadCharacter
   end interface

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: POP_ConfigOpen
! !INTERFACE:

 subroutine POP_ConfigOpen(iunit, errorCode, configFileName)

! !DESCRIPTION:
!  This routine opens a configuration input file for reading.
!  If no filename is supplied, the default filename is used.
!  The unit number for the input file is returned for use by
!  the configuration read routines.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in), optional :: &
      configFileName

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      iunit,                   &! I/O unit for config file
      errorCode                 ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: ierr  ! status flag from open

!-----------------------------------------------------------------------
!
!  get unit number for input file to be read
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_IOUnitsGet(iunit)

!-----------------------------------------------------------------------
!
!  if this is master task, open the file for reading and check for 
!  errors
!
!-----------------------------------------------------------------------

   configFileDefault = 'pop2_in' // trim(inst_suffix)

   if (POP_myTask == POP_masterTask) then
      if (present(configFileName)) then
         open(unit=iunit, file=configFileName, form='formatted', &
              status='old', action='read', position='rewind',    &
              iostat=ierr)
      else
         open(unit=iunit, file=configFileDefault, form='formatted', &
              status='old', action='read', position='rewind',       &
              iostat=ierr)
      endif
   endif

   call POP_Broadcast(ierr, POP_masterTask, errorCode)

   if (ierr > 0 .or. errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_ConfigOpen: error opening config file')
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_ConfigOpen

!***********************************************************************
!BOP
! !IROUTINE: POP_ConfigClose
! !INTERFACE:

 subroutine POP_ConfigClose(iunit, errorCode)

! !DESCRIPTION:
!  This routine closes an open configuration file and releases
!  the assigned unit.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      iunit                     ! I/O unit for config file

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode                 ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  if this is master task, open the file for reading and check for 
!  errors
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (POP_myTask == POP_masterTask) close(iunit)

!-----------------------------------------------------------------------
!
!  release unit number
!
!-----------------------------------------------------------------------

   call POP_IOUnitsRelease(iunit)

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_ConfigClose

!***********************************************************************
!BOP
! !IROUTINE: POP_ConfigRead
! !INTERFACE:

 subroutine POP_ConfigReadI4(iunit, moduleName, variableName,   &
                             variable, defaultValue, errorCode, &
                             outStringBefore, outStringAfter)

! !DESCRIPTION:
!  This routine reads a variable from a configuration input file
!  that has already been opened with a ConfigOpen call.  Each variable
!  in the input file is associated with a module, so the module name
!  must also be supplied.  If the variable is not present in the
!  input file, the defaultValue is assigned.  After a successful read, 
!  the value for the variable is broadcasted to other processors.
!  Finally, the value is printing to stdout using either a generic
!  output string or user-specified output defined by outStringBefore
!  and outStringAfter.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETER:

   integer (POP_i4), intent(in) :: &
      iunit                    ! i/o unit of config file

   character (*), intent(in) :: &
      moduleName,               &! name of module where this var resides
      variableName               ! name of variable to be read

   integer (POP_i4), intent(in) :: &
      defaultValue               ! default value to assign to variable

   character (*), intent(in), optional :: &
      outStringBefore,   &! optional output string to precede variable value
      outStringAfter      ! optional output string to follow  variable value

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      variable                  ! variable to assing input value

   integer (POP_i4), intent(out) :: &
      errorCode                 ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (POP_logical) :: &
      moduleFound,          &! logical flag for module   search
      variableFound,        &! logical flag for variable search
      isOpen                 ! logical flag for file inquiry

   integer (POP_i4) :: &
      istat,                &! I/O status flag
      indx,                 &! index for manipulating string
      errVal                 ! internal error flag

   character (POP_charLength) :: &
      inputString,          &! temp for reading each record
      tmpString              ! temp for manipulating input string

!-----------------------------------------------------------------------
!
!  check to see if unit is open and rewind unit
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   errVal = 0

   if (POP_myTask == POP_masterTask) then
      inquire(unit=iunit, opened=isOpen)
      if (isOpen) then
         rewind(iunit)
      else
         errVal = -1
      endif
   endif

!-----------------------------------------------------------------------
!
!  look for module name
!
!-----------------------------------------------------------------------

   moduleFound = .false.
   if (POP_myTask == POP_masterTask .and. isOpen) then
      moduleSearch: do 

         ! read line from input file
         read(iunit, '(a100)', iostat=istat) inputString

         ! check for read errors
         if (istat < 0) then ! end of file
            exit moduleSearch
         else if (istat > 0) then ! error reading from file
            errVal = -2
            exit moduleSearch
         endif

         ! look for ampersand, signifying a module name
         ! then check module name for a match
         tmpString = adjustl(inputString)
         if (tmpString(1:1) == '&') then
            if (trim(tmpString(2:)) == trim(moduleName)) then
               moduleFound = .true.
               exit moduleSearch
            endif
         else
            cycle moduleSearch
         endif
         
      end do moduleSearch

      if (.not. moduleFound) then
         errVal = -3
      endif
   endif

!-----------------------------------------------------------------------
!
!  look for variable name
!
!-----------------------------------------------------------------------

   variableFound = .false.
   if (POP_myTask == POP_masterTask .and. moduleFound) then
      varSearch: do 

         ! read line from input file: should be name = value
         read(iunit, '(a100)', iostat=istat) inputString

         ! check for read errors
         if (istat < 0) then ! end of file
            exit varSearch
         else if (istat > 0) then ! error reading from file
            errVal = -2
            exit varSearch
         endif

         ! check for end of module block
         tmpString = adjustl(inputString)
         if (tmpString(1:1) == '/') exit varSearch

         ! then check for a variable name match
         indx = index(tmpString,'=')
         if (trim(adjustl(tmpString(1:indx-1))) == &
             trim(variableName)) then
            variableFound = .true.
            exit varSearch
         endif
      end do varSearch
   endif

!-----------------------------------------------------------------------
!
!  check for errors
!
!-----------------------------------------------------------------------

   call POP_Broadcast(errVal, POP_masterTask, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error broadcasting error value')
      return
   endif

   select case(errVal)
   case (-1)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: config file not opened for reading')
      return
   case (-2)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error reading record from config file')
      return
   case (-3)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: module name not found in config file')
      return
   case default
   end select

!-----------------------------------------------------------------------
!
!  extract value from input string or set value to default
!
!-----------------------------------------------------------------------

   if (POP_myTask == POP_masterTask) then

      if (variableFound) then
         read(tmpString(indx+1:),*) variable
      else
         variable = defaultValue
         write(POP_stdout, '(a37,a)') &
            '   Using default value for variable: ', variableName
      endif
   endif

!-----------------------------------------------------------------------
!
!  broadcast value to all processors
!
!-----------------------------------------------------------------------

   call POP_Broadcast(variable, POP_masterTask, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error broadcasting variable')
      return
   endif

!-----------------------------------------------------------------------
!
!  output value to stdout
!
!-----------------------------------------------------------------------

   if (POP_myTask == POP_masterTask) then
      if (present(outStringBefore)) then
         tmpString = outStringBefore
      else
         tmpString(:) = ' '
         write(tmpString,'(a,a3)') variableName,' = '
      endif

      if (present(outStringAfter)) then
         write(POP_stdout,'(a,a1,i10,a)') trim(tmpString), ' ', &
                                          variable, outStringAfter
      else
         write(POP_stdout,'(a,a1,i10)'  ) trim(tmpString), ' ', &
                                          variable
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_ConfigReadI4

!***********************************************************************
!BOP
! !IROUTINE: POP_ConfigRead
! !INTERFACE:

 subroutine POP_ConfigReadR4(iunit, moduleName, variableName,   &
                             variable, defaultValue, errorCode, &
                             outStringBefore, outStringAfter)

! !DESCRIPTION:
!  This routine reads a variable from a configuration input file
!  that has already been opened with a ConfigOpen call.  Each variable
!  in the input file is associated with a module, so the module name
!  must also be supplied.  If the variable is not present in the
!  input file, the defaultValue is assigned.  After a successful read, 
!  the value for the variable is broadcasted to other processors.
!  Finally, the value is printing to stdout using either a generic
!  output string or user-specified output defined by outStringBefore
!  and outStringAfter.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETER:

   integer (POP_i4), intent(in) :: &
      iunit                    ! i/o unit of config file

   character (*), intent(in) :: &
      moduleName,               &! name of module where this var resides
      variableName               ! name of variable to be read

   real (POP_r4), intent(in) :: &
      defaultValue               ! default value to assign to variable

   character (*), intent(in), optional :: &
      outStringBefore,   &! optional output string to precede variable value
      outStringAfter      ! optional output string to follow  variable value

! !OUTPUT PARAMETERS:

   real (POP_r4), intent(out) :: &
      variable                  ! variable to assing input value

   integer (POP_i4), intent(out) :: &
      errorCode                 ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (POP_logical) :: &
      moduleFound,          &! logical flag for module   search
      variableFound,        &! logical flag for variable search
      isOpen                 ! logical flag for file inquiry

   integer (POP_i4) :: &
      istat,                &! I/O status flag
      indx,                 &! index for manipulating string
      errVal                 ! internal error flag

   character (POP_charLength) :: &
      inputString,          &! temp for reading each record
      tmpString              ! temp for manipulating input string

!-----------------------------------------------------------------------
!
!  check to see if unit is open and rewind unit
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   errVal = 0

   if (POP_myTask == POP_masterTask) then
      inquire(unit=iunit, opened=isOpen)
      if (isOpen) then
         rewind(iunit)
      else
         errVal = -1
      endif
   endif

!-----------------------------------------------------------------------
!
!  look for module name
!
!-----------------------------------------------------------------------

   moduleFound = .false.
   if (POP_myTask == POP_masterTask .and. isOpen) then
      moduleSearch: do 

         ! read line from input file
         read(iunit, '(a100)', iostat=istat) inputString

         ! check for read errors
         if (istat < 0) then ! end of file
            exit moduleSearch
         else if (istat > 0) then ! error reading from file
            errVal = -2
            exit moduleSearch
         endif

         ! look for ampersand, signifying a module name
         ! then check module name for a match
         tmpString = adjustl(inputString)
         if (tmpString(1:1) == '&') then
            if (trim(tmpString(2:)) == trim(moduleName)) then
               moduleFound = .true.
               exit moduleSearch
            endif
         else
            cycle moduleSearch
         endif
         
      end do moduleSearch

      if (.not. moduleFound) then
         errVal = -3
      endif
   endif

!-----------------------------------------------------------------------
!
!  look for variable name
!
!-----------------------------------------------------------------------

   variableFound = .false.
   if (POP_myTask == POP_masterTask .and. moduleFound) then
      varSearch: do 

         ! read line from input file: should be name = value
         read(iunit, '(a100)', iostat=istat) inputString

         ! check for read errors
         if (istat < 0) then ! end of file
            exit varSearch
         else if (istat > 0) then ! error reading from file
            errVal = -2
            exit varSearch
         endif

         ! check for end of module block
         tmpString = adjustl(inputString)
         if (tmpString(1:1) == '/') exit varSearch

         ! then check for a variable name match
         indx = index(tmpString,'=')
         if (trim(adjustl(tmpString(1:indx-1))) == &
             trim(variableName)) then
            variableFound = .true.
            exit varSearch
         endif
      end do varSearch
   endif

!-----------------------------------------------------------------------
!
!  check for errors
!
!-----------------------------------------------------------------------

   call POP_Broadcast(errVal, POP_masterTask, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error broadcasting error value')
      return
   endif

   select case(errVal)
   case (-1)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: config file not opened for reading')
      return
   case (-2)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error reading record from config file')
      return
   case (-3)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: module name not found in config file')
      return
   case default
   end select

!-----------------------------------------------------------------------
!
!  extract value from input string or set value to default
!
!-----------------------------------------------------------------------

   if (POP_myTask == POP_masterTask) then

      if (variableFound) then
         read(tmpString(indx+1:),*) variable
      else
         variable = defaultValue
         write(POP_stdout, '(a37,a)') &
            '   Using default value for variable: ', variableName
      endif
   endif

!-----------------------------------------------------------------------
!
!  broadcast value to all processors
!
!-----------------------------------------------------------------------

   call POP_Broadcast(variable, POP_masterTask, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error broadcasting variable')
      return
   endif

!-----------------------------------------------------------------------
!
!  output value to stdout
!
!-----------------------------------------------------------------------

   if (POP_myTask == POP_masterTask) then
      if (present(outStringBefore)) then
         tmpString = outStringBefore
      else
         write(tmpString,'(a,a3)') variableName,' = '
      endif

      if (present(outStringAfter)) then
         write(POP_stdout,'(a,a1,1pe12.5,a)') trim(tmpString), ' ', &
                                              variable, outStringAfter
      else
         write(POP_stdout,'(a,a1,1pe12.5)'  ) trim(tmpString), ' ', &
                                              variable
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_ConfigReadR4

!***********************************************************************
!BOP
! !IROUTINE: POP_ConfigRead
! !INTERFACE:

 subroutine POP_ConfigReadR8(iunit, moduleName, variableName,   &
                             variable, defaultValue, errorCode, &
                             outStringBefore, outStringAfter)

! !DESCRIPTION:
!  This routine reads a variable from a configuration input file
!  that has already been opened with a ConfigOpen call.  Each variable
!  in the input file is associated with a module, so the module name
!  must also be supplied.  If the variable is not present in the
!  input file, the defaultValue is assigned.  After a successful read, 
!  the value for the variable is broadcasted to other processors.
!  Finally, the value is printing to stdout using either a generic
!  output string or user-specified output defined by outStringBefore
!  and outStringAfter.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETER:

   integer (POP_i4), intent(in) :: &
      iunit                    ! i/o unit of config file

   character (*), intent(in) :: &
      moduleName,               &! name of module where this var resides
      variableName               ! name of variable to be read

   real (POP_r8), intent(in) :: &
      defaultValue               ! default value to assign to variable

   character (*), intent(in), optional :: &
      outStringBefore,   &! optional output string to precede variable value
      outStringAfter      ! optional output string to follow  variable value

! !OUTPUT PARAMETERS:

   real (POP_r8), intent(out) :: &
      variable                  ! variable to assign input value

   integer (POP_i4), intent(out) :: &
      errorCode                 ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (POP_logical) :: &
      moduleFound,          &! logical flag for module   search
      variableFound,        &! logical flag for variable search
      isOpen                 ! logical flag for file inquiry

   integer (POP_i4) :: &
      istat,                &! I/O status flag
      indx,                 &! index for manipulating string
      errVal                 ! internal error flag

   character (POP_charLength) :: &
      inputString,          &! temp for reading each record
      tmpString              ! temp for manipulating input string

!-----------------------------------------------------------------------
!
!  check to see if unit is open and rewind unit
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   errVal = 0

   if (POP_myTask == POP_masterTask) then
      inquire(unit=iunit, opened=isOpen)
      if (isOpen) then
         rewind(iunit)
      else
         errVal = -1
      endif
   endif

!-----------------------------------------------------------------------
!
!  look for module name
!
!-----------------------------------------------------------------------

   moduleFound = .false.
   if (POP_myTask == POP_masterTask .and. isOpen) then
      moduleSearch: do 

         ! read line from input file
         read(iunit, '(a100)', iostat=istat) inputString

         ! check for read errors
         if (istat < 0) then ! end of file
            exit moduleSearch
         else if (istat > 0) then ! error reading from file
            errVal = -2
            exit moduleSearch
         endif

         ! look for ampersand, signifying a module name
         ! then check module name for a match
         tmpString = adjustl(inputString)
         if (tmpString(1:1) == '&') then
            if (trim(tmpString(2:)) == trim(moduleName)) then
               moduleFound = .true.
               exit moduleSearch
            endif
         else
            cycle moduleSearch
         endif
         
      end do moduleSearch

      if (.not. moduleFound) then
         errVal = -3
      endif
   endif

!-----------------------------------------------------------------------
!
!  look for variable name
!
!-----------------------------------------------------------------------

   variableFound = .false.
   if (POP_myTask == POP_masterTask .and. moduleFound) then
      varSearch: do 

         ! read line from input file: should be name = value
         read(iunit, '(a100)', iostat=istat) inputString

         ! check for read errors
         if (istat < 0) then ! end of file
            exit varSearch
         else if (istat > 0) then ! error reading from file
            errVal = -2
            exit varSearch
         endif

         ! check for end of module block
         tmpString = adjustl(inputString)
         if (tmpString(1:1) == '/') exit varSearch

         ! then check for a variable name match
         indx = index(tmpString,'=')
         if (trim(adjustl(tmpString(1:indx-1))) == &
             trim(variableName)) then
            variableFound = .true.
            exit varSearch
         endif
      end do varSearch
   endif

!-----------------------------------------------------------------------
!
!  check for errors
!
!-----------------------------------------------------------------------

   call POP_Broadcast(errVal, POP_masterTask, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error broadcasting error value')
      return
   endif

   select case(errVal)
   case (-1)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: config file not opened for reading')
      return
   case (-2)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error reading record from config file')
      return
   case (-3)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: module name not found in config file')
      return
   case default
   end select

!-----------------------------------------------------------------------
!
!  extract value from input string or set value to default
!
!-----------------------------------------------------------------------

   if (POP_myTask == POP_masterTask) then

      if (variableFound) then
         read(tmpString(indx+1:),*) variable
      else
         variable = defaultValue
         write(POP_stdout, '(a37,a)') &
            '   Using default value for variable: ', variableName
      endif
   endif

!-----------------------------------------------------------------------
!
!  broadcast value to all processors
!
!-----------------------------------------------------------------------

   call POP_Broadcast(variable, POP_masterTask, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error broadcasting variable')
      return
   endif

!-----------------------------------------------------------------------
!
!  output value to stdout
!
!-----------------------------------------------------------------------

   if (POP_myTask == POP_masterTask) then
      if (present(outStringBefore)) then
         tmpString = outStringBefore
      else
         write(tmpString,'(a,a3)') variableName,' = '
      endif

      if (present(outStringAfter)) then
         write(POP_stdout,'(a,a1,1pe22.15,a)') trim(tmpString), ' ', &
                                               variable, outStringAfter
      else
         write(POP_stdout,'(a,a1,1pe22.15)'  ) trim(tmpString), ' ', &
                                               variable
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_ConfigReadR8

!***********************************************************************
!BOP
! !IROUTINE: POP_ConfigRead
! !INTERFACE:

 subroutine POP_ConfigReadLogical(iunit, moduleName, variableName,  &
                             variable, defaultValue, errorCode,     &
                             outStringBefore, outStringAfter)

! !DESCRIPTION:
!  This routine reads a variable from a configuration input file
!  that has already been opened with a ConfigOpen call.  Each variable
!  in the input file is associated with a module, so the module name
!  must also be supplied.  If the variable is not present in the
!  input file, the defaultValue is assigned.  After a successful read, 
!  the value for the variable is broadcasted to other processors.
!  Finally, the value is printing to stdout using either a generic
!  output string or user-specified output defined by outStringBefore
!  and outStringAfter.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETER:

   integer (POP_i4), intent(in) :: &
      iunit                    ! i/o unit of config file

   character (*), intent(in) :: &
      moduleName,               &! name of module where this var resides
      variableName               ! name of variable to be read

   logical (POP_logical), intent(in) :: &
      defaultValue               ! default value to assign to variable

   character (*), intent(in), optional :: &
      outStringBefore,   &! optional output string to precede variable value
      outStringAfter      ! optional output string to follow  variable value

! !OUTPUT PARAMETERS:

   logical (POP_logical), intent(out) :: &
      variable                  ! variable to assign input value

   integer (POP_i4), intent(out) :: &
      errorCode                 ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (POP_logical) :: &
      moduleFound,          &! logical flag for module   search
      variableFound,        &! logical flag for variable search
      isOpen                 ! logical flag for file inquiry

   integer (POP_i4) :: &
      istat,                &! I/O status flag
      indx,                 &! index for manipulating string
      errVal                 ! internal error flag

   character (POP_charLength) :: &
      inputString,          &! temp for reading each record
      tmpString              ! temp for manipulating input string

!-----------------------------------------------------------------------
!
!  check to see if unit is open and rewind unit
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   errVal = 0

   if (POP_myTask == POP_masterTask) then
      inquire(unit=iunit, opened=isOpen)
      if (isOpen) then
         rewind(iunit)
      else
         errVal = -1
      endif
   endif

!-----------------------------------------------------------------------
!
!  look for module name
!
!-----------------------------------------------------------------------

   moduleFound = .false.
   if (POP_myTask == POP_masterTask .and. isOpen) then
      moduleSearch: do 

         ! read line from input file
         read(iunit, '(a100)', iostat=istat) inputString

         ! check for read errors
         if (istat < 0) then ! end of file
            exit moduleSearch
         else if (istat > 0) then ! error reading from file
            errVal = -2
            exit moduleSearch
         endif

         ! look for ampersand, signifying a module name
         ! then check module name for a match
         tmpString = adjustl(inputString)
         if (tmpString(1:1) == '&') then
            if (trim(tmpString(2:)) == trim(moduleName)) then
               moduleFound = .true.
               exit moduleSearch
            endif
         else
            cycle moduleSearch
         endif
         
      end do moduleSearch

      if (.not. moduleFound) then
         errVal = -3
      endif
   endif

!-----------------------------------------------------------------------
!
!  look for variable name
!
!-----------------------------------------------------------------------

   variableFound = .false.
   if (POP_myTask == POP_masterTask .and. moduleFound) then
      varSearch: do 

         ! read line from input file: should be name = value
         read(iunit, '(a100)', iostat=istat) inputString

         ! check for read errors
         if (istat < 0) then ! end of file
            exit varSearch
         else if (istat > 0) then ! error reading from file
            errVal = -2
            exit varSearch
         endif

         ! check for end of module block
         tmpString = adjustl(inputString)
         if (tmpString(1:1) == '/') exit varSearch

         ! then check for a variable name match
         indx = index(tmpString,'=')
         if (trim(adjustl(tmpString(1:indx-1))) == &
             trim(variableName)) then
            variableFound = .true.
            exit varSearch
         endif
      end do varSearch
   endif

!-----------------------------------------------------------------------
!
!  check for errors
!
!-----------------------------------------------------------------------

   call POP_Broadcast(errVal, POP_masterTask, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error broadcasting error value')
      return
   endif

   select case(errVal)
   case (-1)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: config file not opened for reading')
      return
   case (-2)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error reading record from config file')
      return
   case (-3)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: module name not found in config file')
      return
   case default
   end select

!-----------------------------------------------------------------------
!
!  extract value from input string or set value to default
!
!-----------------------------------------------------------------------

   if (POP_myTask == POP_masterTask) then

      if (variableFound) then
         read(tmpString(indx+1:),*) variable
      else
         variable = defaultValue
         write(POP_stdout, '(a37,a)') &
            '   Using default value for variable: ', variableName
      endif
   endif

!-----------------------------------------------------------------------
!
!  broadcast value to all processors
!
!-----------------------------------------------------------------------

   call POP_Broadcast(variable, POP_masterTask, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error broadcasting variable')
      return
   endif

!-----------------------------------------------------------------------
!
!  output value to stdout
!
!-----------------------------------------------------------------------

   if (POP_myTask == POP_masterTask) then
      if (present(outStringBefore)) then
         tmpString = outStringBefore
      else
         write(tmpString,'(a,a3)') variableName,' = '
      endif

      if (present(outStringAfter)) then
         if (variable) then
            write(POP_stdout,'(a,a5,a)') trim(tmpString), ' true', &
                                         outStringAfter
         else
            write(POP_stdout,'(a,a6,a)') trim(tmpString), ' false', &
                                         outStringAfter
         endif
      else
         if (variable) then
            write(POP_stdout,'(a,a5)') trim(tmpString), ' true'
         else
            write(POP_stdout,'(a,a6)') trim(tmpString), ' false'
         endif
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_ConfigReadLogical

!***********************************************************************
!BOP
! !IROUTINE: POP_ConfigRead
! !INTERFACE:

 subroutine POP_ConfigReadCharacter(iunit, moduleName, variableName, &
                             variable, defaultValue, errorCode,      &
                             outStringBefore, outStringAfter)

! !DESCRIPTION:
!  This routine reads a variable from a configuration input file
!  that has already been opened with a ConfigOpen call.  Each variable
!  in the input file is associated with a module, so the module name
!  must also be supplied.  If the variable is not present in the
!  input file, the defaultValue is assigned.  After a successful read, 
!  the value for the variable is broadcasted to other processors.
!  Finally, the value is printing to stdout using either a generic
!  output string or user-specified output defined by outStringBefore
!  and outStringAfter.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETER:

   integer (POP_i4), intent(in) :: &
      iunit                    ! i/o unit of config file

   character (*), intent(in) :: &
      moduleName,               &! name of module where this var resides
      variableName               ! name of variable to be read

   character (*), intent(in) :: &
      defaultValue               ! default value to assign to variable

   character (*), intent(in), optional :: &
      outStringBefore,   &! optional output string to precede variable value
      outStringAfter      ! optional output string to follow  variable value

! !OUTPUT PARAMETERS:

   character (POP_charLength), intent(out) :: &
      variable                  ! variable to assign input value

   integer (POP_i4), intent(out) :: &
      errorCode                 ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (POP_logical) :: &
      moduleFound,          &! logical flag for module   search
      variableFound,        &! logical flag for variable search
      isOpen                 ! logical flag for file inquiry

   integer (POP_i4) :: &
      istat,                &! I/O status flag
      indx,                 &! index for manipulating string
      errVal                 ! internal error flag

   character (POP_charLength) :: &
      inputString,          &! temp for reading each record
      tmpString              ! temp for manipulating input string

!-----------------------------------------------------------------------
!
!  check to see if unit is open and rewind unit
!
!-----------------------------------------------------------------------

   errorCode = POP_Success
   errVal = 0

   if (POP_myTask == POP_masterTask) then
      inquire(unit=iunit, opened=isOpen)
      if (isOpen) then
         rewind(iunit)
      else
         errVal = -1
      endif
   endif

!-----------------------------------------------------------------------
!
!  look for module name
!
!-----------------------------------------------------------------------

   moduleFound = .false.
   if (POP_myTask == POP_masterTask .and. isOpen) then
      moduleSearch: do 

         ! read line from input file
         read(iunit, '(a100)', iostat=istat) inputString

         ! check for read errors
         if (istat < 0) then ! end of file
            exit moduleSearch
         else if (istat > 0) then ! error reading from file
            errVal = -2
            exit moduleSearch
         endif

         ! look for ampersand, signifying a module name
         ! then check module name for a match
         tmpString = adjustl(inputString)
         if (tmpString(1:1) == '&') then
            if (trim(tmpString(2:)) == trim(moduleName)) then
               moduleFound = .true.
               exit moduleSearch
            endif
         else
            cycle moduleSearch
         endif
         
      end do moduleSearch

      if (.not. moduleFound) then
         errVal = -3
      endif
   endif

!-----------------------------------------------------------------------
!
!  look for variable name
!
!-----------------------------------------------------------------------

   variableFound = .false.
   if (POP_myTask == POP_masterTask .and. moduleFound) then
      varSearch: do 

         ! read line from input file: should be name = value
         read(iunit, '(a100)', iostat=istat) inputString

         ! check for read errors
         if (istat < 0) then ! end of file
            exit varSearch
         else if (istat > 0) then ! error reading from file
            errVal = -2
            exit varSearch
         endif

         ! check for end of module block
         tmpString = adjustl(inputString)
         if (tmpString(1:1) == '/') exit varSearch

         ! then check for a variable name match
         indx = index(tmpString,'=')
         if (trim(adjustl(tmpString(1:indx-1))) == &
             trim(variableName)) then
            variableFound = .true.
            exit varSearch
         endif
      end do varSearch
   endif

!-----------------------------------------------------------------------
!
!  check for errors
!
!-----------------------------------------------------------------------

   call POP_Broadcast(errVal, POP_masterTask, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error broadcasting error value')
      return
   endif

   select case(errVal)
   case (-1)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: config file not opened for reading')
      return
   case (-2)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error reading record from config file')
      return
   case (-3)
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: module name not found in config file')
      return
   case default
   end select

!-----------------------------------------------------------------------
!
!  extract value from input string or set value to default
!
!-----------------------------------------------------------------------

   if (POP_myTask == POP_masterTask) then

      if (variableFound) then
         read(tmpString(indx+1:),*) variable
      else
         variable = defaultValue
         write(POP_stdout, '(a37,a)') &
            '   Using default value for variable: ', variableName
      endif
   endif

!-----------------------------------------------------------------------
!
!  broadcast value to all processors
!
!-----------------------------------------------------------------------

   call POP_Broadcast(variable, POP_masterTask, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_ConfigRead: error broadcasting variable')
      return
   endif

!-----------------------------------------------------------------------
!
!  output value to stdout
!
!-----------------------------------------------------------------------

   if (POP_myTask == POP_masterTask) then
      if (present(outStringBefore)) then
         tmpString = outStringBefore
      else
         write(tmpString,'(a,a3)') variableName,' = '
      endif

      if (present(outStringAfter)) then
         write(POP_stdout,'(a,a1,a,a1,a)') trim(tmpString), ' ', &
                                           trim(variable),  ' ', &
                                           outStringAfter
      else
         write(POP_stdout,'(a,a1,a)'  ) trim(tmpString), ' ', &
                                        trim(variable)
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_ConfigReadCharacter

!***********************************************************************

 end module POP_ConfigMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
