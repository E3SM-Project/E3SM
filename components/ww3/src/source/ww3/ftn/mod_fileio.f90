!-----------------------------------------------------------------------------!
module m_fileio
!-----------------------------------------------------------------------------!
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 8 Feb. 2003
!   +---+ |   |
!         +---+
!
! Module for storing file i/o related variables
!
! The values for the parameter i_log, i_prt and iw_tst must be set
! in one of the routines of the host program or in subroutine sys_init
!
! Version 1.1   29 May  2000   Initial version
!         1.2   21 Sep. 2001   Form=binary added (B)
!         1.3    5 Oct. 2001   Form=direct access, unformatted, fixed record (R)
!         1.4   24 Aug. 2002   Bug fixed and restructure of test output
!         1.5    8 Feb. 2003   Error check included when incorrect path (Z_FILEIO)
!
!-----------------------------------------------------------------------------!
! The following two parameters must be set by the user
! They define the overall test level and the output channel
!
integer,parameter :: i_print=0  ! (0/1/2)  Test output printing off/on
!                               ! Output channel defined by i_out
!
integer,parameter :: i_out=6    ! Output channel to screen
!                               ! ==1 screen output for Unix/Linux systems
!                               ! ==6 screen output for Windows
!------------------------------------------------------------------------------
!
! Standard switches to activate Logging, Test and Print ouput
!
integer i_log        ! (0/1)      Logging off/on
integer i_prt        ! (0/1)      Printing off/on
integer i_tst        ! (0,1,2...) Test level off/on
!
!
! Standard unit numbers of input & output files
!
integer lu_err  ! standard error file
integer lu_inp  ! standard input file
integer lu_log  ! standard logging
integer lu_prt  ! standard print output
integer lu_tst  ! standard test output
!
contains
!-----------------------------------------------------------------------------!
subroutine z_fileio(filename,qual,iufind,iunit,iostat)                        !
!-----------------------------------------------------------------------------!
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+
!   +---+ |   |
!         +---+
!
implicit none
!
!  0. Update history
!
!     24/07/1999  First version
!     28/09/1999  Module name changed from FILEOPEN -> Z_FILEIO
!     27/10/1999  Option to delete an existing file added
!     18/11/1999  Argument IUNIT used to control use of Z_FLUNIT
!     22/11/1999  Parameter iunit not changed unless by z_flunit
!     28/12/1999  Interface with Z_FLUNIT updated and
!                 input parameter iufind added
!     14/04/2000  Module m_fileio included in this routine
!     25/05/2000  Module m_fileio excluded, if an already opened file is
!                 found, the corresponding unit number is assigned to output
!     21/09/2001  Form=binary added, extension to Fortran 95 standard
!      5/10/2001  Form=fixed Record length, as specified in input argument
!     17/06/2002  Initialisation of IUNIT=-1 included
!     24/08/2002  Bug fixed when routine called with IUFIND=0
!     08/02/2003  Bug fixed when file could not be created due to invalid path 
!
!  1. Purpose
!
!     Open file with name FILENAME and determine unit number IUNIT
!     With file type determined in QUAL
!
!     Depending on the value of IUFIND a search is performed for a
!     free unit number
!
!  2. Method
!
!     If file exists then
!       if QUAL = 'D'
!         delete file
!       Else
!         inquire if file opened
!         If opened
!           determine unit number
!         Else
!           If iunit >= 10 Find free unit number
!           Open file with unit number and file qualifier
!         End if
!       End if
!     Else
!       If QUAL='SNU'
!         If iunit >= 10 find free unit number
!         Open new file with unit number and qualifier
!       Else
!         Iunit = -1   File does not exist
!       End if
!     End if
!
!
!  3. Parameter list
!
!Type       I/O           Name           Description
!----------------------------------------------------
character(len=*), intent(in)  :: filename  ! File name
character(len=2), intent(in)  :: qual      ! File qualifyer
integer,   intent(in)         :: iufind    ! Indicator for search of unit number
integer,   intent(inout)      :: iunit     ! Unit number
integer,   intent(out)        :: iostat    ! Error indicator
!
!  4. Subroutines used
!
!     Z_FLUNIT
!
!  5. Error messages
!
!     IUNIT > 0    File exists, is (already) connected to unit number IUNIT, or is
!                  created and connected to unit number
!     IUNIT == 0   File has been deleted or does not exist
!           < 0    An error occurred, no file or unit number found
!
!     IOSTAT =  0   No errors detected
!              -1   Incorrect file qualifier
!              -2   Unit number does not exist
!              -3   Attempt to open non-existing file with status=OLD
!              -4   Attempt to open existing file with wrong FORMATTING
!              -5   Incorrect value for IUFIND: not in range [0,1]
!              -6   File could not be created due to,e.g. incorrect path
!
!  6. Remarks
!
!     1) Use of file qualifier:
!
!        1st char: O(ld),R(eplace),S(cratch),
!                  U(nknown),(D)elete
!        2nd char: F(ormatted),U(nformatted),B(inary)
!
!     2) Use of IUFIND
!
!        if IUFIND==0, No search is performed for a free unit number
!                 ==1, A search is performed in routine Z_FLUNIT
!
!     3) This routine is based on routine FOR from
!        SWAN version 40.00 of Delft University of Technology
!
!------------------------------------------------------------------------------
! Local variables
!
character(len=7)  :: cstat  ! string with status of file I/O
character(len=11) :: cform  ! string with format of file I/O
integer junit               ! temporary unit number
logical lexist              ! indicator if a file exists
logical lopen               ! indicator if a file is opened
integer iuerr               ! error indicator from Z_FLUNIT
!-------------------------------------------------------------------------------------
! initialisations
!-------------------------------------------------------------------------------------
iostat = 0
if(iufind==1) iunit  = -1
!
!
!  Check value of IUFIND
!
if(iufind/=0 .and. iufind/=1) then
  if(i_print >0) write(i_out,*) 'Z_FILEIO: Incorrect value for IUFIND:',iufind
  iostat = -5
  goto 9999
end if
!
!
!  check input argument QUAL
!
if(i_print>=1) write(i_out,*) 'Z_FILEIO/A:',trim(filename),' ',qual,iunit,iostat
!
if (index('ORSUD',qual(1:1)) ==0 .or. index('FUB',qual(2:2)) ==0) then
  if(i_print > 0) write(i_out,*) 'Incorrect file qualifier'
  iostat = -1
else
  if(qual(1:1) == 'O') cstat = 'old'
  if(qual(1:1) == 'R') cstat = 'replace'
  if(qual(1:1) == 'S') cstat = 'scratch'
  if(qual(1:1) == 'U') cstat = 'unknown'
  if(qual(1:1) == 'D') cstat = 'delete'
!
  if(qual(2:2) == 'F') cform = 'formatted'
  if(qual(2:2) == 'U') cform = 'unformatted'
  if(qual(2:2) == 'B') cform = 'binary'          ! extension to FORTRAN 95 standard
  if(qual(2:2) == 'R') cform = 'unformatted'
!
!  Check if file exists
!
  inquire(file=filename,exist=lexist)
  if(i_print >=2) write(i_out,*) 'Z_FILEIO  file exists?:',trim(filename),':',lexist
!
!  delete file if it exists and qual == 'D'
!
  if(lexist .and. qual(1:1)=='D') then
    inquire(file=filename,opened=lopen)
    if(lopen) then
      inquire(file=filename,number=junit)
    else
      if(iufind == 1) call z_flunit(iunit,iuerr)
      junit = iunit
      if(junit > 0) then
        open(file=filename,unit=junit,form=cform,iostat=iostat)
        if(iostat/=0) then
          iostat = -4
          goto 9999
        end if
      end if
    end if
    close(junit,status=cstat)
    goto 9999
  end if
!
!  if the file exists, check if it is opened
!
  if(lexist) then
    if(i_print >=2) write(i_out,*) 'Z_FILEIO: File exists:',trim(filename)
    inquire(file=filename,opened=lopen)
    if(lopen) then
      if(i_print >=2) write(i_out,*) 'Z_FILEIO: File is opened:',trim(filename)
!
!  determine unit number to which this file is connected
!  and assign it to the output number
!
      inquire(file=filename,number=junit)
      if(i_print >=2) write(i_out,*) 'Z_FILEIO: File is connected to unit:', junit
      iunit = junit
    else
!
!  if the file exists and not connected to a unit number, search a free unit number
!
      if(i_print >=2) write(i_out,*) 'Z_FILEIO: File is not connected to a unit number'
      if(iufind==0) then
        if(i_print >=2) write(i_out,*) 'Z_FILEIO: Assign user defined unit number:',iunit
      elseif(iufind==1) then
        call z_flunit(iunit,iuerr)
        if(i_print >=2) write(i_out,*) 'Z_FILEIO: New unit number IUNIT:',iunit
      end if
      junit = iunit
!
      if(junit > 0) then
        open(file=filename,unit=junit,form=cform,status=cstat)
      else
        iostat = -2
      end if
   end if
!
!  the file does not exist, so open it and find a free unit number
!
  else
!
    if(i_print>=2) then
       write(i_out,*) 'Z_FILEIO: File does not exist !'
       write(i_out,*) 'Z_FILEIO: Qual:',qual(1:1)
    end if
!
    if(index('SRU',qual(1:1)) > 0) then
      if(iufind==1) then
        call z_flunit(iunit,iuerr)
        if(i_print >=1) write(i_out,*) 'Z_FILEIO: New unit number IUNIT:',iunit
      end if
      junit = iunit
!
!  open file to IUNIT, if possible
!
      if(junit > 0) then
        open(file=filename,unit=junit,form=cform,iostat=iuerr)
!
! check added 8/2/2003
!
        if(iuerr/=0) then
          iunit = -1
          iostat = -6
        end if
      else
        iostat = -2
      end if
!
!  file cannot be opened because it does not exist
!
    elseif('O'==qual(1:1)) then   ! File should exist
      if(i_print>=2) write(i_out,*) 'Z_FILEIO: File cannot be opened because it does not exist'
      iostat = -3
    end if
  end if
end if
!
9999 continue
!
if(i_print>=1) write(i_out,*) 'Z_FILEIO/Z:',trim(filename),' ',qual,iunit,iostat
!
return
end  subroutine
!
!-----------------------------------------------------------------------------!
subroutine z_fclose(iunit)                                                    !
!-----------------------------------------------------------------------------!
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+
!   +---+ |   |
!         +---+
!
implicit none
!
!  0. Update history
!
!     0.01 24/08/2000  First version
!
!  1. Purpose
!
!     Close file with unit number IUNIT, and set IUNIT=-1
!
!  2. Method
!
!
!  3. Parameter list
!
!Type       I/O           Name           Description
!-----------------------------------------------------------------------------
integer, intent(inout)  :: iunit          ! Unit number
!-----------------------------------------------------------------------------
close(iunit)
iunit = -1
!
return
end subroutine
!
!-----------------------------------------------------------------------------!
subroutine z_flunit(iunit,ierr)                                               !
!-----------------------------------------------------------------------------!
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+
!   +---+ |   |
!         +---+
!
implicit none
!
!  0. Update history
!
!     Version   Date    Modification
!
!     0.01  24/07/1999  Initial version
!     0.02  01/10/1999  Extra check added to ensure maximum unit number
!     0.03  07/10/1999  Check of existence of uni number deleted,
!                       since this test produces different answer
!                       on Lahey compiler
!     0.04  25/11/1999  Intent added
!     0.05  24/12/1999  Module M_GENVAR added for information about range of unit numbers
!     0.06  27/12/1999  Module M_GENVAR replaced by M_FILEIO
!                       Check added for forbidden unit numbers
!     0.07  28/12/1999  Internal checks added and IERR added to parameter list
!     0.08  08/02/2000  User of lu_min & lu_max deleted
!     0.09  14/04/2000  Module m_fileio included in this routine
!
!  1. Purpose
!
!     Find a free unit number
!
!  2. Method
!
!     Starting at LU_MIN till LU_MAX are investigated until
!     a free (i.e. not connected to a file) is found.
!     Use is made of the standard fortran INQUIRE function.
!     The values of LU_MIN and LU_MAX should be specified
!     in an initialisation routine
!
!  3. Parameter list
!
!Type     I/O           Name         Description
!----------------------------------------------------------
integer, intent(out) :: iunit       ! resulting unit number
integer, intent(out) :: ierr        ! error level
!
!  4. Subroutines used
!
!     None
!
!  5. Error messages
!
!     ierr=0   No errors encountered
!          1   Invalud combination lu_low >= lu_high
!          2   Invalid value for lu_low
!          3   Invalid value for lu_high
!          4   No free unit number could be found
!
!  6. Remarks
!
!     If no free unit number if found in the range
!     lu_min - lu_high, then the function returns IUNIT = -1
!
!     The switch i_print can be used to generate test output
!
!----------------------------------------------------------------------------------
! local parameters
!
integer junit                       ! counter for unit numbers
logical lopen                       ! indicator if a unit number is connected to a file
logical lnot                        ! indicates if a forbidden unit number is checked
integer i_not                       ! counter to check forbidded unit numbers
!
!---------------------------------------------------------------------------------
!  range of unit numbers to search
!
integer, parameter :: lu_min=60     ! minimum unit number
integer, parameter :: lu_max=200    ! maximum unit number
!
! specification of forbidden unit numbers
!
integer, parameter :: lu_nr=3   ! number of forbidden unit numbers
integer lu_not(lu_nr)           ! list of forbidden unit numbers
!----------------------------------------------------------------------------------
lu_not(1) = 100
lu_not(2) = 101
lu_not(3) = 102
!-----------------------------------------------------------------------------------
!
ierr = 0
!
if(i_print >= 2) then
  write(i_out,*) 'Z_FLUNIT: forbidden     :',lu_not
  write(i_out,*) 'Z_FLUNIT: lu_min lu_max :',lu_min,lu_max
end if
!
!  check data specified in Module Z_FILEIO
!
if(lu_min >= lu_max) then
  ierr = 1
  write(i_out,*) 'Z_FLUNIT: Incorrect boundaries for LU_MIN & LU_MAX:',&
&            lu_min,lu_max
end if
!
junit = lu_min
!
iunit = -1
!
do while (iunit ==-1)
!
! Check if unit number is free, i.e. not in use by an opened file
!
   inquire(unit=junit,opened=lopen)
!
!  check if unit number is not a forbidden unit number
!
   lnot = .false.
   do i_not=1,lu_nr
     if(lu_not(i_not)==junit) then
       lnot = .true.
       if(i_print >= 1) write(i_out,*) 'Z_FLUNIT: a forbidden unit number was encountered:',junit
     end if
   end do
!
   if(lopen.or.lnot) then
      junit = junit + 1
   else
      iunit = junit
   end if
   if(junit > lu_max) exit
end do
!
if(iunit < 0) then
  write(i_out,*) 'ERROR in Z_FLUNIT: No free unit number could be found'
end if
!
return
end subroutine
!
end module
