module Input_Aux_module

  use Option_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public :: input_type 
    PetscInt :: fid
    PetscErrorCode :: ierr
    character(len=MAXWORDLENGTH) :: filename
    character(len=MAXSTRINGLENGTH) :: buf
    character(len=MAXSTRINGLENGTH) :: err_buf
    character(len=MAXSTRINGLENGTH) :: err_buf2
    PetscBool :: broadcast_read
    PetscBool :: force_units ! force user to declare units on datasets
    type(input_type), pointer :: parent
  end type input_type

  type :: input_dbase_type
    character(len=MAXWORDLENGTH), pointer :: icard(:)
    character(len=MAXWORDLENGTH), pointer :: rcard(:)
    character(len=MAXWORDLENGTH), pointer :: ccard(:)
    PetscInt, pointer :: ivalue(:)
    PetscReal, pointer :: rvalue(:)
    character(len=MAXWORDLENGTH), pointer :: cvalue(:)
  end type input_dbase_type

  type(input_dbase_type), pointer, public :: dbase => null()

  interface InputReadWord
    module procedure InputReadWord1
    module procedure InputReadWord2
  end interface
  
  interface InputReadNChars
    module procedure InputReadNChars1
    module procedure InputReadNChars2
  end interface
  
  interface InputReadInt
    module procedure InputReadInt1
    module procedure InputReadInt2
#if defined(PETSC_USE_64BIT_INDICES) && (PETSC_SIZEOF_MPI_FINT * PETSC_BITS_PER_BYTE != 64)
    ! If PetscInt and PetscMPIInt have different sizes (occurs for some builds
    ! with 64 bit indices), then we need to have additional routines for the 
    ! InputReadInt() generic subroutine.  (We use the above check instead of 
    ! directly checking to see if PetscInt and PetscMPIInt have the same size
    ! because the size of PetscInt is not included in the 
    ! $PETSC_DIR/$PETSC_ARCH/include/petscconf.h file.) If the two types have
    ! the same size, then these additional routines for type PetscMPIInt must
    ! *not* be defined, because then the interface becomes ambiguous, since 
    ! Fortran doesn't know the difference between PetscInt and PetscMPIInt if
    ! they are identically sized integers.  --RTM
    module procedure InputReadInt3
    module procedure InputReadInt4
#endif
  end interface
  
  interface InputReadDouble
    module procedure InputReadDouble1
    module procedure InputReadDouble2
  end interface
  
  interface InputReadNDoubles
    module procedure InputReadNDoubles1
    module procedure InputReadNDoubles2
  end interface
  
  interface InputError
    module procedure InputError1
    module procedure InputError2
  end interface
  
  interface InputErrorMsg
    module procedure InputErrorMsg1
    module procedure InputErrorMsg2
  end interface
  
  interface InputDefaultMsg
    module procedure InputDefaultMsg1
    module procedure InputDefaultMsg2
  end interface
  
  interface InputReadStringErrorMsg
    module procedure InputReadStringErrorMsg1
    module procedure InputReadStringErrorMsg2
  end interface
  
  interface InputFindStringInFile
    module procedure InputFindStringInFile1
    module procedure InputFindStringInFile2
  end interface

  interface InputKeywordUnrecognized
    module procedure InputKeywordUnrecognized1
    module procedure InputKeywordUnrecognized2
  end interface
  
  public :: InputCreate, InputDestroy, InputReadPflotranString, &
            InputReadWord, InputReadDouble, InputReadInt, InputCheckExit, &
            InputReadNDoubles, &
            InputSkipToEND, InputFindStringInFile, InputErrorMsg, &
            InputDefaultMsg, InputReadStringErrorMsg, &
            InputFindStringErrorMsg, InputError, &
            InputReadNChars, InputReadQuotedWord, &
            InputReadPath, &
            InputGetCommandLineInt, &
            InputGetCommandLineReal, &
            InputGetCommandLineTruth, &
            InputGetCommandLineString, &
            InputReadFilenames, &
            InputGetLineCount, &
            InputReadToBuffer, &
            InputReadASCIIDbase, &
            InputKeywordUnrecognized, &
            InputCheckMandatoryUnits, &
            InputDbaseDestroy, &
            InputPushExternalFile, &
            InputReadWordDbaseCompatible, &
            InputReadAndConvertUnits

contains

! ************************************************************************** !

function InputCreate(fid,filename,option)
  ! 
  ! Allocates and initializes a new Input object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  use Option_module

  implicit none
  
  PetscInt :: fid
  character(len=*) :: filename
  type(option_type) :: option
  
  type(input_type), pointer :: InputCreate
  PetscInt :: status  
  type(input_type), pointer :: input
  
  allocate(input)

  input%fid = fid
  input%filename = filename
  input%ierr = 0
  input%buf = ''
  input%err_buf = ''
  input%err_buf2 = ''
  input%broadcast_read = PETSC_FALSE
  input%force_units = PETSC_FALSE
  nullify(input%parent)
  
  if (fid == MAX_IN_UNIT) then
    option%io_buffer = 'MAX_IN_UNIT in pflotran_constants.h must be increased to' // &
      ' accommodate a larger number of embedded files.'
    call printErrMsg(option)
  endif

  open(unit=input%fid,file=filename,status="old",iostat=status)
  if (status /= 0) then
    if (len_trim(filename) == 0) filename = '<blank>'
    option%io_buffer = 'File: "' // trim(filename) // '" not found.'
    call printErrMsg(option)
  endif
  
  InputCreate => input
  
end function InputCreate

! ************************************************************************** !

subroutine InputDefaultMsg1(input,option,buffer)
  ! 
  ! If ierr /= 0, informs user that default value will be used.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=*) :: buffer

  if (InputError(input)) then
    input%err_buf = buffer
    call InputDefaultMsg(input,option)
  endif

end subroutine InputDefaultMsg1

! ************************************************************************** !

subroutine InputDefaultMsg2(input,option)
  ! 
  ! If ierr /= 0, informs user that default value will be used.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option

  if (InputError(input)) then
    option%io_buffer =  '"' // trim(input%err_buf) // &
                        '" set to default value.'
    call printMsg(option)
    input%ierr = 0
  endif

end subroutine InputDefaultMsg2

! ************************************************************************** !

subroutine InputErrorMsg1(input,option,buffer1,buffer2)
  ! 
  ! If ierr /= 0, If ierr /= 0, informs user of error and stops.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=*) :: buffer1, buffer2

  if (InputError(input)) then
    input%err_buf = buffer1
    input%err_buf2 = buffer2
    call InputErrorMsg(input,option)
  endif

end subroutine InputErrorMsg1

! ************************************************************************** !

subroutine InputErrorMsg2(input,option)
  ! 
  ! InputErrorMsg: If ierr /= 0, If ierr /= 0, informs user of error and stops.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option

  if (InputError(input)) then
    option%io_buffer = 'While reading "' // trim(input%err_buf) // &
                       '" under keyword: ' // trim(input%err_buf2) // '.'
    call printErrMsg(option)
  endif

end subroutine InputErrorMsg2

! ************************************************************************** !

subroutine InputReadStringErrorMsg1(input, option, buffer)
  ! 
  ! If ierr /= 0, informs user of error and stops.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=*) :: buffer

  if (InputError(input)) then
    input%err_buf = buffer
    call InputReadStringErrorMsg(input, option)
  endif

end subroutine InputReadStringErrorMsg1

! ************************************************************************** !

subroutine InputReadStringErrorMsg2(input, option)
  ! 
  ! If ierr /= 0, informs user of error and stops.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option

  if (InputError(input)) then
    option%io_buffer = 'While reading in string in "' // &
                       trim(input%err_buf) // '".'
    call printErrMsg(option)
  endif

end subroutine InputReadStringErrorMsg2

! ************************************************************************** !

subroutine InputFindStringErrorMsg(input, option, string)
  ! 
  ! If ierr /= 0, informs user of error and stops.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=*) :: string

  if (InputError(input)) then
    option%io_buffer = 'Card (' // trim(string) // ') not &
                       &found in file.'
    call printErrMsg(option)    
  endif

end subroutine InputFindStringErrorMsg

! ************************************************************************** !

subroutine InputReadInt1(input, option, int)
  ! 
  ! reads and removes an integer value from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscInt :: int

  character(len=MAXWORDLENGTH) :: word
  PetscBool :: found

  found = PETSC_FALSE
  if (associated(dbase)) then
    call InputParseDbaseForInt(input%buf,int,found,input%ierr)
  endif
  
  if (.not.found) then
    call InputReadWord(input%buf,word,PETSC_TRUE,input%ierr)
  
    if (.not.InputError(input)) then
      read(word,*,iostat=input%ierr) int
    endif
  endif

end subroutine InputReadInt1

! ************************************************************************** !

subroutine InputReadInt2(string, option, int, ierr)
  ! 
  ! reads and removes an integer value from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscInt :: int
  PetscErrorCode :: ierr

  character(len=MAXWORDLENGTH) :: word
  PetscBool :: found

  ierr = 0

  found = PETSC_FALSE
  if (associated(dbase)) then
    call InputParseDbaseForInt(string,int,found,ierr)
  endif
  
  if (.not.found) then
    call InputReadWord(string,word,PETSC_TRUE,ierr)
  
    if (.not.InputError(ierr)) then
      read(word,*,iostat=ierr) int
    endif
  endif

end subroutine InputReadInt2

#if defined(PETSC_USE_64BIT_INDICES) && (PETSC_SIZEOF_MPI_FINT * PETSC_BITS_PER_BYTE != 64)

! ************************************************************************** !

subroutine InputReadInt3(input, option, int)
  ! 
  ! InputReadInt3() and InputReadInt4() must only be defined if PetscInt and
  ! PetscMPIInt differ in size.  See notes above in the interface definition.
  ! --RTM
  ! reads and removes an integer value from a string
  ! authors: Glenn Hammond, Richard Mills
  ! 
  ! Date: 2/3/2012
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscMPIInt :: int

  character(len=MAXWORDLENGTH) :: word

  call InputReadWord(input%buf,word,PETSC_TRUE,input%ierr)
  
  if (.not.InputError(input)) then
    read(word,*,iostat=input%ierr) int
  endif

end subroutine InputReadInt3

! ************************************************************************** !

subroutine InputReadInt4(string, option, int, ierr)
  ! 
  ! reads and removes an integer value from a string
  ! authors: Glenn Hammond, Richard Mills
  ! 
  ! Date: 2/3/2012
  ! 

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscMPIInt :: int
  PetscErrorCode :: ierr

  character(len=MAXWORDLENGTH) :: word

  ierr = 0
  call InputReadWord(string,word,PETSC_TRUE,ierr)
  
  if (.not.InputError(ierr)) then
    read(word,*,iostat=ierr) int
  endif

end subroutine InputReadInt4

#endif
! End of defined(PETSC_USE_64BIT_INDICES) &&
! (PETSC_SIZEOF_MPI_FINT * PETSC_BITS_PER_BYTE != 64) conditional

! ************************************************************************** !

subroutine InputReadDouble1(input, option, double)
  ! 
  ! reads and removes a real value from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscReal :: double

  character(len=MAXWORDLENGTH) :: word
  PetscBool :: found

  found = PETSC_FALSE
  if (associated(dbase)) then
    call InputParseDbaseForDouble(input%buf,double,found,input%ierr)
  endif
  
  if (.not.found) then
    call InputReadWord(input%buf,word,PETSC_TRUE,input%ierr)
  
    if (.not.InputError(input)) then
      read(word,*,iostat=input%ierr) double
    endif
  endif

end subroutine InputReadDouble1

! ************************************************************************** !

subroutine InputReadDouble2(string, option, double, ierr)
  ! 
  ! reads and removes a real value from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscReal :: double
  PetscErrorCode :: ierr

  character(len=MAXWORDLENGTH) :: word
  PetscBool :: found

  ierr = 0
  
  found = PETSC_FALSE
  if (associated(dbase)) then
    call InputParseDbaseForDouble(string,double,found,ierr)
  endif
  
  if (.not.found) then
    call InputReadWord(string,word,PETSC_TRUE,ierr)
  
    if (.not.InputError(ierr)) then
      read(word,*,iostat=ierr) double
    endif
  endif

end subroutine InputReadDouble2

! ************************************************************************** !

subroutine InputReadNDoubles1(input, option, double, n)
  ! 
  ! reads and removes "n" real value from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/29/11
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscInt :: n
  PetscReal :: double(n)

  PetscInt :: i

  do i = 1, n
    call InputReadDouble(input,option,double(i))
    if (InputError(input)) return
  enddo

end subroutine InputReadNDoubles1

! ************************************************************************** !

subroutine InputReadNDoubles2(string, option, double, n, ierr)
  ! 
  ! reads and removes "n" real values from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/29/11
  ! 

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscInt :: n
  PetscReal :: double(n)
  PetscErrorCode :: ierr

  PetscInt :: i

  do i = 1, n
    call InputReadDouble(string,option,double(i),ierr)
    if (InputError(ierr)) return
  enddo

end subroutine InputReadNDoubles2

! ************************************************************************** !

subroutine InputReadPflotranString(input, option)
  ! 
  ! Reads a string (strlen characters long) from a
  ! file while avoiding commented or skipped lines.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  PetscInt :: flag

  if (input%broadcast_read) then
    if (option%myrank == option%io_rank) then
      call InputReadPflotranStringSlave(input, option)
    endif
    flag = input%ierr
    call MPI_Bcast(flag,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                   option%mycomm,ierr)
    input%ierr = flag
    if (.not.InputError(input)) then  
      call MPI_Bcast(input%buf,MAXSTRINGLENGTH,MPI_CHARACTER, &
                     option%io_rank,option%mycomm,ierr)      
    endif
  else
    call InputReadPflotranStringSlave(input, option)
  endif

end subroutine InputReadPflotranString

! ************************************************************************** !

subroutine InputReadPflotranStringSlave(input, option)
  ! 
  ! Reads a string (strlen characters long) from a
  ! file while avoiding commented or skipped lines.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  use String_module
  
  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) ::  tempstring
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: i
  PetscInt :: skip_count

  input%ierr = 0

! we initialize the word to blanks to avoid error reported by valgrind
!  do i=1,MAXWORDLENGTH
!     word(i:i) = ' '
!  enddo
  word = ''
  
  do
    read(input%fid,'(a512)',iostat=input%ierr) input%buf
    call StringAdjustl(input%buf)

    if (InputError(input)) then
      ! check to see if another file is on the stack
      if (InputPopExternalFile(input)) then
        cycle
      else
        exit
      endif
    endif

    if (input%buf(1:1) == '#' .or. input%buf(1:1) == '!') cycle

    tempstring = input%buf
    call InputReadWord(tempstring,word,PETSC_TRUE,input%ierr)
    call StringToUpper(word)
    
    if (word(1:13) == 'EXTERNAL_FILE') then
      ! have to stip the card 'EXTERNAL_FILE' from the buffer
      call InputReadWord(input,option,word,PETSC_TRUE)
      ! push a new input file to stack
      call InputPushExternalFile(input,option)
      cycle
    else if (word(1:4) == 'SKIP') then
      skip_count = 1
      do 
        read(input%fid,'(a512)',iostat=input%ierr) tempstring
        if (InputError(input)) then
          option%io_buffer = 'End of file reached in ' // &
              'InputReadPflotranStringSlave.  SKIP encountered ' // &
              'without a matching NOSKIP.'
          call printErrMsg(option)              
        endif
        call InputReadWord(tempstring,word,PETSC_FALSE,input%ierr)
        call StringToUpper(word)
        if (word(1:4) == 'SKIP') skip_count = skip_count + 1
        if (word(1:4) == 'NOSK') then
          skip_count = skip_count - 1
          if (skip_count == 0) exit
        endif
      enddo
      if (InputError(input)) exit
    else if (word(1:1) /= ' ' .and. word(1:4) /= 'NOSK') then
      exit
    endif
  enddo
  
  ! Check for comment midway along a string
  if (.not.InputError(input)) then
    tempstring = input%buf
    input%buf = repeat(' ',MAXSTRINGLENGTH)
    do i=1,len_trim(tempstring)
      if (tempstring(i:i) /= '#' .and. tempstring(i:i) /= '!') then
        input%buf(i:i) = tempstring(i:i)
      else
        exit
      endif
    enddo
  endif

end subroutine InputReadPflotranStringSlave

! ************************************************************************** !

subroutine InputReadWord1(input, option, word, return_blank_error)
  ! 
  ! reads and removes a word (consecutive characters) from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  PetscBool :: return_blank_error
  
  if (InputError(input)) return
  
  call InputReadWord2(input%buf, word, return_blank_error, input%ierr)

end subroutine InputReadWord1

! ************************************************************************** !

subroutine InputReadWord2(string, word, return_blank_error, ierr)
  ! 
  ! reads and removes a word (consecutive characters) from a
  ! string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  character(len=*) :: string
  character(len=*) :: word
  PetscBool :: return_blank_error
  PetscErrorCode :: ierr
  
  PetscInt :: i, begins, ends, length
  character(len=1), parameter :: tab = achar(9), backslash = achar(92)

  if (ierr /= 0) return

  ! Initialize character string to blank.
  ! Initialize character string to blank.  len_trim(word) is not
  ! defined if word is allocated but not initialized.  This works on
  ! most compilers, but may not work on some?  Holler if it
  ! errors... - etc
  word = ''
  ! do i=1,len_trim(word)
  !   word(i:i) = ' '
  ! enddo

  length = len_trim(string)
  
  if (length == 0) then
    if (return_blank_error) then
      ierr = 1
    else
      ierr = 0
    endif
    return
  else
    ierr = 0

    ! Remove leading blanks and tabs
    i=1
    do while((string(i:i) == ' ' .or. string(i:i) == ',' .or. &
             string(i:i) == tab) .and. i <= length) 
      i=i+1
    enddo

    if (i > length) then
      if (return_blank_error) then
        ierr = 1
      else
        ierr = 0
      endif
      return
    endif
    
    begins=i

    ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (string(i:i) /= ' ' .and. string(i:i) /= ',' .and. &
              string(i:i) /= tab .and. &
              (i == begins .or. string(i:i) /= backslash))
      i=i+1
    enddo

    ends=i-1

    ! Avoid copying beyond the end of the word (32 characters).
    if (ends-begins > (MAXWORDLENGTH-1)) ends = begins + (MAXWORDLENGTH-1)

    ! Copy (ends-begins) characters to 'word'
    word = string(begins:ends)
    ! Remove chars from string
    string = string(ends+1:)

  endif

end subroutine InputReadWord2

! ************************************************************************** !

subroutine InputReadWordDbaseCompatible(input, option, word, &
                                        return_blank_error)
  ! 
  ! reads a word and checks whether there is an entry in the Dbase with which
  ! to swap
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/22/16
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  PetscBool :: return_blank_error
  
  PetscBool :: found

  if (InputError(input)) return

  found = PETSC_FALSE
  if (associated(dbase)) then
    call InputParseDbaseForWord(input%buf,word,found,input%ierr)
  endif
  
  if (.not.found) then
    call InputReadWord(input%buf,word,PETSC_TRUE,input%ierr)
  endif

end subroutine InputReadWordDbaseCompatible

! ************************************************************************** !

subroutine InputReadNChars1(input, option, chars, n, return_blank_error)
  ! 
  ! reads and removes a specified number of characters from a
  ! string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/00
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscBool :: return_blank_error ! Return an error for a blank line
                                   ! Therefore, a blank line is not acceptable.
  
  PetscInt :: n, begins, ends
  character(len=n) :: chars

  if (InputError(input)) return

  call InputReadNChars2(input%buf, chars, n, return_blank_error, input%ierr)
  
end subroutine InputReadNChars1

! ************************************************************************** !

subroutine InputReadNChars2(string, chars, n, return_blank_error, ierr)
  ! 
  ! reads and removes a specified number of characters from a
  ! string
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/00
  ! 

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: return_blank_error ! Return an error for a blank line
                                   ! Therefore, a blank line is not acceptable.
  
  PetscInt :: i, n, begins, ends
  character(len=n) :: chars
  PetscErrorCode :: ierr
  character(len=1), parameter :: tab = achar(9), backslash = achar(92)

  if (InputError(ierr)) return

  ! Initialize character string to blank.
  chars(1:n) = repeat(' ',n)

  ierr = len_trim(string)
  if (.not.InputError(ierr)) then
    if (return_blank_error) then
      ierr = 1
    else
      ierr = 0
    endif
    return
  else
    ierr = 0

    ! Remove leading blanks and tabs
    i=1
    do while(string(i:i) == ' ' .or. string(i:i) == tab) 
      i=i+1
    enddo

    begins=i

    ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (string(i:i) /= ' ' .and. string(i:i) /= ',' .and. &
              string(i:i) /= tab  .and. &
              (i == begins .or. string(i:i) /= backslash))
      i=i+1
    enddo

    ends=i-1

    if (ends-begins+1 > n) then ! string read is too large for 'chars'
      ierr = 1
      return
    endif

    ! Copy (ends-begins) characters to 'chars'
    chars = string(begins:ends)
    ! Remove chars from string
    string = string(ends+1:)

  endif

end subroutine InputReadNChars2

! ************************************************************************** !

subroutine InputReadQuotedWord(input, option, word, return_blank_error)
  ! 
  ! reads and removes a word from a string, that is
  ! delimited by "'".
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/07/00
  ! 

  implicit none

  type(input_type) :: input
  type(option_type) :: option
  PetscInt :: i, begins, ends, realends, len_trim_word
  PetscBool :: return_blank_error ! Return an error for a blank line
                                ! Therefore, a blank line is not acceptable.
  character(len=*) :: word
  PetscBool :: openquotefound
  character(len=1), parameter :: tab = achar(9), backslash = achar(92)

  if (InputError(input)) return

  openquotefound = PETSC_FALSE
  ! Initialize character string to blank.
  len_trim_word = len_trim(word)
  word(1:len_trim_word) = repeat(' ',len_trim_word)
  
  if (len_trim(input%buf) == 0) then
    if (return_blank_error) then
      input%ierr = 1
    else
      input%ierr = 0
    endif
    return
  else
    input%ierr = 0  
    
    ! Remove leading blanks and tabs
    i=1
    do while(input%buf(i:i) == ' ' .or. input%buf(i:i) == tab) 
      i=i+1
    enddo

    if (input%buf(i:i) == "'") then
      openquotefound = PETSC_TRUE
      i=i+1
    endif

    begins=i

    if (openquotefound) then
      do while (input%buf(i:i) /= "'")
        if (i > (MAXWORDLENGTH-1)) exit
        i=i+1
      enddo
    else
    ! Count # of continuous characters (no blanks, commas, etc. in between)
      do while (input%buf(i:i) /= ' ' .and. input%buf(i:i) /= ',' .and. &
                input%buf(i:i) /= tab .and. &
                (i == begins .or. input%buf(i:i) /= backslash))
        i=i+1
      enddo
    endif

    realends = i
    ends=i-1

    ! Avoid copying beyond the end of the word (32 characters).
    if (ends-begins > (MAXWORDLENGTH-1)) ends = begins + (MAXWORDLENGTH-1)

    ! Copy (ends-begins) characters to 'chars'
    word = input%buf(begins:ends)
    ! Remove chars from string
    input%buf = input%buf(realends+1:)
  endif

end subroutine InputReadQuotedWord

! ************************************************************************** !

subroutine InputReadPath(string, word, return_blank_error, ierr)
  ! 
  ! reads and removes a words from a path
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/14/10
  ! 

  implicit none

  character(len=*) :: string
  character(len=*) :: word
  PetscBool :: return_blank_error
  PetscErrorCode :: ierr
  
  PetscInt :: i, begins, ends, len_trim_word
  character(len=1), parameter :: slash = achar(47), backslash = achar(92)

  if (ierr /= 0) return

  ! Initialize character string to blank.
  len_trim_word = len_trim(word)
  word(1:len_trim_word) = repeat(' ',len_trim_word)

  ierr = len_trim(string)
  
  if (ierr == 0) then
    if (return_blank_error) then
      ierr = 1
    else
      ierr = 0
    endif
    return
  else
    ierr = 0

    ! Remove leading blanks and tabs
    i=1
    do while(string(i:i) == ' ' .and. string(i:i) == slash) 
      i=i+1
    enddo

    begins=i

    ! Count # of characters (no slashes in between)
    do while (string(i:i) /= slash .and. &
              (i == begins .or. string(i:i) /= backslash))
      i=i+1
    enddo

    ends=i-1

    ! Avoid copying beyond the end of the word (32 characters).
    if (ends-begins > (MAXWORDLENGTH-1)) ends = begins + (MAXWORDLENGTH-1)

    ! Copy (ends-begins) characters to 'word'
    word = string(begins:ends)
    ! Remove chars from string
    string = string(ends+1:)

  endif
  
end subroutine InputReadPath

! ************************************************************************** !

subroutine InputFindStringInFile1(input, option, string)
  ! 
  ! Rewinds file and finds the first occurrence of
  ! 'string'.  Note that the line must start with 'string'
  ! in order to match and that line is NOT returned
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/07
  ! 

  use String_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  
  call InputFindStringInFile2(input, option, string, PETSC_TRUE)
  
end subroutine InputFindStringInFile1

! ************************************************************************** !

subroutine InputFindStringInFile2(input, option, string, print_warning)
  ! 
  ! Rewinds file and finds the first occurrence of
  ! 'string'.  Note that the line must start with 'string'
  ! in order to match and that line is NOT returned
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/07/07
  ! 

  use String_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: print_warning
  
  character(len=MAXWORDLENGTH) :: word
  PetscBool :: found = PETSC_FALSE
  PetscInt :: length1, length2

  input%ierr = 0

  length1 = len_trim(string)

  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    call InputReadWord(input,option,word,PETSC_TRUE)
    if (InputError(input)) exit
    length2 = len_trim(word)
    if (length1 == length2 .and. StringCompare(string,word,length1)) then
      found = PETSC_TRUE
      exit
    endif
  enddo
  
  ! if not found, rewind once and try again.  this approach avoids excessive 
  ! reading if successive searches for strings are in descending order in 
  ! the file.
  if (InputError(input)) then
    input%ierr = 0
    rewind(input%fid)
    do 
      call InputReadPflotranString(input,option)
      if (InputError(input)) exit
      call InputReadWord(input,option,word,PETSC_TRUE)
      if (InputError(input)) exit
      length2 = len_trim(word)
      if (length1 == length2 .and. StringCompare(string,word,length1)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
  endif    
  
  if (.not.found .and. print_warning) then
    option%io_buffer = 'Card (' // trim(string) // ') not found in input file.'
    call printWrnMsg(option)
    input%ierr = 1
  endif
  
end subroutine InputFindStringInFile2

! ************************************************************************** !

subroutine InputSkipToEND(input,option,string)
  ! 
  ! Skips to keyword END
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  ! 

  implicit none
  
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=*) :: string

  do
    call InputReadPflotranString(input,option)
    input%err_buf = 'End of file found before end of card ' // trim(string)
    call InputReadStringErrorMsg(input,option)
    if (InputCheckExit(input,option)) exit
  enddo

end subroutine InputSkipToEND

! ************************************************************************** !

function InputCheckExit(input,option)
  ! 
  ! Checks whether an end character (.,/,'END') has been found
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 

  use String_module
  
  implicit none

  type(input_type) :: input
  type(option_type) :: option  
  PetscInt :: i
  character(len=1) :: tab
  
  PetscBool :: InputCheckExit

  ! We must remove leading blanks and tabs. --RTM
  tab = achar(9)
  i=1
  do while(input%buf(i:i) == ' ' .or. input%buf(i:i) == tab) 
    i=i+1
  enddo

  if (input%buf(i:i) == '/' .or. &
!geh: this fails when the keyword starts with END
!geh      StringCompare(input%buf(i:),'END',THREE_INTEGER)) then
      StringCompare(input%buf(i:),'END') .or. &
      ! to end a block, e.g. END_SUBSURFACE
      StringStartsWith(input%buf(i:),'END_')) then
    InputCheckExit = PETSC_TRUE
  else
    InputCheckExit = PETSC_FALSE
  endif

end function InputCheckExit

! ************************************************************************** !

function InputError1(input)
  ! 
  ! Returns true if an error has occurred
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/08
  ! 

  implicit none

  type(input_type) :: input
  
  PetscBool :: InputError1

  if (input%ierr == 0) then
    InputError1 = PETSC_FALSE
  else
    InputError1 = PETSC_TRUE
  endif

end function InputError1

! ************************************************************************** !

function InputError2(ierr)
  ! 
  ! Returns true if an error has occurred
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/08
  ! 

  implicit none

  PetscErrorCode :: ierr
  
  PetscBool :: InputError2

  if (ierr == 0) then
    InputError2 = PETSC_FALSE
  else
    InputError2 = PETSC_TRUE
  endif

end function InputError2

! ************************************************************************** !

subroutine InputGetCommandLineInt(string,int_value,found,option)
  ! 
  ! Returns integer value associated with a command
  ! line argument
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/09
  ! 

  use String_module
  use Option_module

  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscBool :: found
  PetscInt :: int_value

  PetscInt :: iarg, narg
  character(len=MAXSTRINGLENGTH) :: string2
  PetscErrorCode :: ierr
  
  ierr = 0
  ! do not initialize int_value, as it may already have a value
  found = PETSC_FALSE
  narg = getCommandLineArgumentCount()
  string = adjustl(string)
  do iarg = 1, narg
    call getCommandLineArgument(iarg,string2)
    if (StringCompare(string,string2)) then
      found = PETSC_TRUE
      if (iarg+1 <= narg) then
        call getCommandLineArgument(iarg+1,string2)
        call InputReadInt(string2,option,int_value,ierr)
      else
        ierr = 1
      endif
      if (InputError(ierr)) then
        option%io_buffer = 'Integer argument for command line argument "' // &
                           trim(adjustl(string)) // '" not found.'
        call printErrMsg(option)
      endif
      exit
    endif
  enddo
  
end subroutine InputGetCommandLineInt

! ************************************************************************** !

subroutine InputGetCommandLineReal(string,double_value,found,option)
  ! 
  ! Returns real*8 value associated with a command
  ! line argument
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/09
  ! 

  use String_module
  use Option_module
  
  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscBool :: found
  PetscReal :: double_value

  PetscInt :: iarg, narg
  character(len=MAXSTRINGLENGTH) :: string2
  PetscErrorCode :: ierr
  
  ierr = 0
  ! do not initialize int_value, as it may already have a value
  found = PETSC_FALSE
  narg = getCommandLineArgumentCount()
  string = adjustl(string)
  do iarg = 1, narg
    call getCommandLineArgument(iarg,string2)
    if (StringCompare(string,string2)) then
      found = PETSC_TRUE
      if (iarg+1 <= narg) then
        call getCommandLineArgument(iarg+1,string2)
        call InputReadDouble(string2,option,double_value,ierr)
      else
        ierr = 1
      endif
      if (InputError(ierr)) then
        option%io_buffer = 'Real argument for command line argument "' // &
                           trim(adjustl(string)) // '" not found.'
        call printErrMsg(option)
      endif
      exit
    endif
  enddo
  
end subroutine InputGetCommandLineReal

! ************************************************************************** !

subroutine InputGetCommandLineString(string,string_value,found,option)
  ! 
  ! Returns a string associated with a command
  ! line argument
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/09
  ! 

  use String_module
  use Option_module
  
  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: string_value

  PetscInt :: iarg, narg
  character(len=MAXSTRINGLENGTH) :: string2
  PetscErrorCode :: ierr
  
  ierr = 0
  ! do not initialize int_value, as it may already have a value
  found = PETSC_FALSE
  narg = getCommandLineArgumentCount()
  string = adjustl(string)
  do iarg = 1, narg
    call getCommandLineArgument(iarg,string2)
    if (StringCompare(string,string2)) then
      found = PETSC_TRUE
      if (iarg+1 <= narg) then
        call getCommandLineArgument(iarg+1,string2)
        call InputReadNChars(string2,string_value,MAXSTRINGLENGTH, &
                             PETSC_TRUE,ierr)
        if (string_value(1:1) == '-') then
          ! no argument exists
          option%io_buffer = 'String argument (' // &
                             trim(adjustl(string_value)) // & 
                             ') for command line argument "' // &
                             trim(adjustl(string)) // '" not recognized.'
          call printErrMsg(option)
        endif
      else
        ierr = 1
      endif
      if (InputError(ierr)) then
        option%io_buffer = 'String argument for command line argument "' // &
                           trim(adjustl(string)) // '" not found.'
        call printErrMsg(option)
      endif
      exit
    endif
  enddo
  
end subroutine InputGetCommandLineString

! ************************************************************************** !

subroutine InputGetCommandLineTruth(string,truth_value,found,option)
  ! 
  ! Returns logical associated with a command
  ! line argument
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/09
  ! 

  use String_module
  use Option_module
  
  implicit none

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type) :: option
  PetscBool :: found
  PetscBool :: truth_value

  PetscInt :: iarg, narg
  character(len=MAXSTRINGLENGTH) :: string2
  character(len=MAXWORDLENGTH) :: word
  PetscErrorCode :: ierr
  
  ierr = 0
  ! do not initialize int_value, as it may already have a value
  found = PETSC_FALSE
  narg = getCommandLineArgumentCount()
  string = adjustl(string)
  do iarg = 1, narg
    call getCommandLineArgument(iarg,string2)
    if (StringCompare(string,string2)) then
      found = PETSC_TRUE
      if (iarg+1 <= narg) then
        call getCommandLineArgument(iarg+1,string2)
        call InputReadWord(string2,word,PETSC_TRUE,ierr)
      else
        ! check if no argument exists, which is valid and means 'true'
        truth_value = PETSC_TRUE
        exit
      endif    
      if (word(1:1) == '-') then
        ! no argument exists, which is valid and means 'true'
        truth_value = PETSC_TRUE
        exit
      endif
      call StringToLower(word)
      select case(trim(word))
        case('yes','true','1','on')
          truth_value = PETSC_TRUE
        case('no','false','0','off')
          truth_value = PETSC_FALSE
        case default
          option%io_buffer = 'Truth argument for command line argument "' // &
                             trim(adjustl(string)) // '" not recognized.'
          call printErrMsg(option)
      end select
    endif
  enddo
  
end subroutine InputGetCommandLineTruth

! ************************************************************************** !

function getCommandLineArgumentCount()
  ! 
  ! Returns the number of command line arguments
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/10
  ! 

  implicit none
  
  integer :: iargc
  
  PetscInt :: getCommandLineArgumentCount
  
  ! initialize to zero
  getCommandLineArgumentCount = 0
  
#if defined(PETSC_HAVE_FORTRAN_GET_COMMAND_ARGUMENT)
  getCommandLineArgumentCount = command_argument_count()
#elif defined(PETSC_HAVE_GETARG)
  getCommandLineArgumentCount = iargc()
#endif

end function getCommandLineArgumentCount

! ************************************************************************** !

subroutine getCommandLineArgument(i,arg)
  ! 
  ! Returns the ith command line argument
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/10
  ! 

  implicit none
  
  PetscInt :: i
  character(len=*) :: arg

  integer*4 :: fortran_int

  fortran_int = i
#if defined(PETSC_HAVE_FORTRAN_GET_COMMAND_ARGUMENT)
  call get_command_argument(fortran_int,arg)
#elif defined(PETSC_HAVE_GETARG)
  call getarg(fortran_int,arg)
#endif

end subroutine getCommandLineArgument

! ************************************************************************** !

subroutine InputReadFilenames(option,filenames)
  ! 
  ! Reads filenames for multi-simulation runs
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/11/09
  ! 

  use Option_module

  type(option_type) :: option
  character(len=MAXSTRINGLENGTH), pointer :: filenames(:)

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: filename_count
  type(input_type), pointer :: input
  PetscBool :: card_found

  input => InputCreate(IN_UNIT,option%input_filename,option)

  string = "FILENAMES"
  call InputFindStringInFile(input,option,string) 

  card_found = PETSC_FALSE
  if (InputError(input)) then
    ! if the FILENAMES card is not included, we will assume that only
    ! filenames exist in the file.
    rewind(input%fid)
  else
    card_found = PETSC_TRUE
  endif
    
  filename_count = 0     
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit  
    call InputReadNChars(input,option,filename,MAXSTRINGLENGTH,PETSC_FALSE)
    filename_count = filename_count + 1
  enddo
  
  allocate(filenames(filename_count))
  filenames = ''
  rewind(input%fid) 

  if (card_found) then
    string = "FILENAMES"
    call InputFindStringInFile(input,option,string) 
  endif
  
  filename_count = 0     
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit  
    call InputReadNChars(input,option,filename,MAXSTRINGLENGTH,PETSC_FALSE)
    filename_count = filename_count + 1
    filenames(filename_count) = filename
  enddo

  call InputDestroy(input)

end subroutine InputReadFilenames

! ************************************************************************** !

function InputGetLineCount(input)

  implicit none
  
  type(input_type), pointer :: input
  PetscInt :: line_count
  PetscInt :: InputGetLineCount

  rewind(input%fid)

  line_count = 0
  do
    read(input%fid, '(a512)', iostat=input%ierr)
    if (InputError(input)) exit
    line_count = line_count + 1
  enddo
  
  InputGetLineCount = line_count
  
end function InputGetLineCount

! ************************************************************************** !

subroutine InputReadToBuffer(input, buffer)

  implicit none
  
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: buffer(:)
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: line

  rewind(input%fid)
  line = 0
  do
    read(input%fid, '(a512)', iostat=input%ierr) string
    if (InputError(input)) exit
    line = line + 1
    buffer(line) = string
  end do
  
end subroutine InputReadToBuffer

! ************************************************************************** !

subroutine InputReadASCIIDbase(filename,option)
  ! 
  ! Read in an ASCII database
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/19/14
  ! 
  use Option_module
  use String_module
  
  implicit none
  
  character(len=*) :: filename
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH), allocatable :: words(:)
  character(len=MAXWORDLENGTH) :: object_name
  type(input_type), pointer :: input
  PetscInt :: icount
  PetscInt :: value_count
  PetscInt :: value_index
  PetscInt :: value_type
  PetscInt :: num_values_in_dataset
  PetscInt :: num_words, num_ints, num_reals
  
  input => InputCreate(IUNIT_TEMP,filename,option)
  
  icount = 0
  num_values_in_dataset = 0
  num_ints = 0
  num_reals = 0
  num_words = 0
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    call InputReadNChars(input,option,string,MAXSTRINGLENGTH,PETSC_FALSE)
    if (len_trim(string) > MAXWORDLENGTH) then
      option%io_buffer = 'ASCII DBASE object names must be shorter than &
        &32 characters: ' // trim(string)
      call printErrMsg(option)
    endif
    word = trim(string)
    if (StringStartsWithAlpha(word)) then
      icount = icount + 1
      if (icount == 1) then
        string = input%buf
        do
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (input%ierr /= 0) exit
          num_values_in_dataset = num_values_in_dataset + 1
        enddo
        input%buf = string
      endif
      input%ierr = 0
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'value','ASCII Dbase')
      select case(StringIntegerDoubleOrWord(word))
        case(STRING_IS_INTEGER)
          num_ints = num_ints + 1
        case(STRING_IS_DOUBLE)
          num_reals = num_reals + 1
        case(STRING_IS_WORD)
          num_words = num_words + 1
      end select
    endif
  enddo

  value_index = 1
  if (option%id > 0) then
    if (option%id > num_values_in_dataset) then
      write(word,*) num_values_in_dataset
        option%io_buffer = 'Data in DBASE_FILENAME "' // &
        trim(filename) // &
        '" is too small (' // trim(adjustl(word)) // &
        ') for number of realizations.'
      call printErrMsg(option)
    endif
    value_index = option%id
  endif
  allocate(words(num_values_in_dataset))
  words = ''
  
  rewind(input%fid)
  allocate(dbase)
  if (num_ints > 0) then
    allocate(dbase%icard(num_ints))
    dbase%icard = ''
    allocate(dbase%ivalue(num_ints))
    dbase%ivalue = UNINITIALIZED_INTEGER
  endif
  if (num_reals > 0) then
    allocate(dbase%rcard(num_reals))
    dbase%rcard = ''
    allocate(dbase%rvalue(num_reals))
    dbase%rvalue = UNINITIALIZED_DOUBLE
  endif
  if (num_words > 0) then
    allocate(dbase%ccard(num_words))
    dbase%ccard = ''
    allocate(dbase%cvalue(num_words))
    dbase%cvalue = '-999'
  endif
  num_ints = 0
  num_reals = 0
  num_words = 0
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    call InputReadWord(input,option,word,PETSC_FALSE)
    if (StringStartsWithAlpha(word)) then
      object_name = word
      words = ''
      value_count = 0
      do
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr /= 0) exit
        value_count = value_count + 1
        if (value_count <= num_values_in_dataset) &
          words(value_count) = word
      enddo
      if (value_count /= num_values_in_dataset) then
        write(word,*) value_count
        option%io_buffer = 'Data in DBASE_FILENAME "' // &
          trim(object_name) // &
          '" has an inconsistent number of values (' // &
          trim(adjustl(word)) // &
          ') for number of realizations ('
        write(word,*) num_values_in_dataset
        option%io_buffer = trim(option%io_buffer) // &
          trim(adjustl(word)) // ').'
        call printErrMsg(option)
      endif
      call StringToUpper(object_name)
      string = words(value_index)
      value_type = StringIntegerDoubleOrWord(string)
      string = words(value_index)
      select case(value_type)
        case(STRING_IS_INTEGER)
          num_ints = num_ints + 1
          dbase%icard(num_ints) = adjustl(object_name)
          call InputReadInt(string,option,dbase%ivalue(num_ints),input%ierr)
          call InputErrorMsg(input,option,'ivalue','ASCII Dbase '//object_name)
        case(STRING_IS_DOUBLE)
          num_reals = num_reals + 1
          dbase%rcard(num_reals) = adjustl(object_name)
          call InputReadDouble(string,option,dbase%rvalue(num_reals),input%ierr)
          call InputErrorMsg(input,option,'rvalue','ASCII Dbase '//object_name)
        case(STRING_IS_WORD)
          num_words = num_words + 1
          dbase%ccard(num_words) = adjustl(object_name)
          dbase%cvalue(num_words) = words(value_index)
      end select
    endif
  enddo
  deallocate(words)
  
  call InputDestroy(input)
  
end subroutine InputReadASCIIDbase

! ************************************************************************** !

subroutine InputParseDbaseForInt(buffer,value,found,ierr)
  ! 
  ! Parses database for an integer value
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/19/14
  ! 
  use String_module
  
  implicit none
  
  character(len=MAXSTRINGLENGTH) :: buffer
  PetscInt :: value
  PetscBool :: found
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: buffer_save
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: dbase_keyword = 'DBASE_VALUE'
  
  buffer_save = buffer
  found = PETSC_FALSE
  call InputReadWord(buffer,word,PETSC_TRUE,ierr)
  if (StringCompareIgnoreCase(word,dbase_keyword)) then
    call InputReadWord(buffer,word,PETSC_TRUE,ierr)
    call DbaseLookupInt(word,value,ierr)
    if (ierr == 0) then
      found = PETSC_TRUE
    endif
  else
    buffer = buffer_save
  endif
  
end subroutine InputParseDbaseForInt

! ************************************************************************** !

subroutine InputParseDbaseForDouble(buffer,value,found,ierr)
  ! 
  ! Parses database for an double precision value
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/19/14
  ! 
  use String_module
  
  implicit none
  
  character(len=MAXSTRINGLENGTH) :: buffer
  PetscReal :: value
  PetscBool :: found
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: buffer_save
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: dbase_keyword = 'DBASE_VALUE'
  
  buffer_save = buffer
  found = PETSC_FALSE
  call InputReadWord(buffer,word,PETSC_TRUE,ierr)
  if (StringCompareIgnoreCase(word,dbase_keyword)) then
    call InputReadWord(buffer,word,PETSC_TRUE,ierr)
    call DbaseLookupDouble(word,value,ierr)
    if (ierr == 0) then
      found = PETSC_TRUE
    endif
  else
    buffer = buffer_save
  endif
  
end subroutine InputParseDbaseForDouble

! ************************************************************************** !

subroutine InputParseDbaseForWord(buffer,value,found,ierr)
  ! 
  ! Parses database for a word
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/22/16
  ! 
  use String_module
  
  implicit none
  
  character(len=MAXSTRINGLENGTH) :: buffer
  character(len=MAXWORDLENGTH) :: value
  PetscBool :: found
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: buffer_save
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: dbase_keyword = 'DBASE_VALUE'
  
  buffer_save = buffer
  found = PETSC_FALSE
  call InputReadWord(buffer,word,PETSC_TRUE,ierr)
  if (StringCompareIgnoreCase(word,dbase_keyword)) then
    call InputReadWord(buffer,word,PETSC_TRUE,ierr)
    call DbaseLookupWord(word,value,ierr)
    if (ierr == 0) then
      found = PETSC_TRUE
    endif
  else
    buffer = buffer_save
  endif
  
end subroutine InputParseDbaseForWord

! ************************************************************************** !

subroutine DbaseLookupInt(keyword,value,ierr)
  ! 
  ! Looks up double precision value in database
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/19/14
  ! 
  use String_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: keyword
  PetscInt :: value
  PetscErrorCode :: ierr
  
  PetscInt :: i
  PetscBool :: found

  ierr = 0
  
  call StringToUpper(keyword)
  
  found = PETSC_FALSE
  if (associated(dbase%icard)) then
    do i = 1, size(dbase%icard)
      if (StringCompare(keyword,dbase%icard(i))) then
        found = PETSC_TRUE
        value = dbase%ivalue(i)
        exit
      endif
    enddo
  endif
  
  if (.not.found) then
    ierr = 1
  endif
  
end subroutine DbaseLookupInt

! ************************************************************************** !

subroutine DbaseLookupDouble(keyword,value,ierr)
  ! 
  ! Looks up double precision value in database
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/19/14
  ! 
  use String_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: keyword
  PetscReal :: value
  PetscErrorCode :: ierr
  
  PetscInt :: i
  PetscBool :: found

  ierr = 0
  
  call StringToUpper(keyword)
  
  found = PETSC_FALSE
  if (associated(dbase%rcard)) then
    do i = 1, size(dbase%rcard)
      if (StringCompare(keyword,dbase%rcard(i))) then
        found = PETSC_TRUE
        value = dbase%rvalue(i)
        exit
      endif
    enddo
  endif
  
  if (.not.found) then
    ierr = 1
  endif
  
end subroutine DbaseLookupDouble

! ************************************************************************** !

subroutine DbaseLookupWord(keyword,value,ierr)
  ! 
  ! Looks up double precision value in database
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/19/14
  ! 
  use String_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: value
  PetscErrorCode :: ierr
  
  PetscInt :: i
  PetscBool :: found

  ierr = 0
  
  call StringToUpper(keyword)
  
  found = PETSC_FALSE
  if (associated(dbase%ccard)) then
    do i = 1, size(dbase%ccard)
      if (StringCompare(keyword,dbase%ccard(i))) then
        found = PETSC_TRUE
        value = dbase%cvalue(i)
        exit
      endif
    enddo
  endif
  
  if (.not.found) then
    ierr = 1
  endif
  
end subroutine DbaseLookupWord

! ************************************************************************** !

subroutine InputKeywordUnrecognized1(keyword,string,option)
  ! 
  ! Looks up double precision value in database
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/19/14
  ! 
  use Option_module
  
  implicit none
  
  character(len=*) :: keyword
  character(len=*) :: string
  type(option_type) :: option

  character(len=1) :: null_string

  null_string = '' 
  call InputKeywordUnrecognized2(keyword,string,null_string,option)
  
end subroutine InputKeywordUnrecognized1

! ************************************************************************** !

subroutine InputKeywordUnrecognized2(keyword,string,string2,option)
  ! 
  ! Looks up double precision value in database
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/19/14
  ! 
  use Option_module
  
  implicit none
  
  character(len=*) :: keyword
  character(len=*) :: string
  character(len=*) :: string2
  type(option_type) :: option
  
  option%io_buffer = 'Keyword "' // &
                     trim(keyword) // &
                     '" not recognized in ' // &
                     trim(string) // '.'
  if (len_trim(string2) > 0) then
    option%io_buffer = trim(option%io_buffer) // ' ' // &
                     trim(string2) // '.'
  endif
  call printErrMsg(option)
  
end subroutine InputKeywordUnrecognized2

! ************************************************************************** !

subroutine InputCheckMandatoryUnits(input,option)
  ! 
  ! Looks up double precision value in database
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/19/14
  ! 
  use Option_module
  
  implicit none
  
  type(input_type) :: input
  type(option_type) :: option
  
  if (input%force_units) then
    option%io_buffer = 'Missing units'
    if (len_trim(input%err_buf) > 1) then
      option%io_buffer = trim(option%io_buffer) // ' in ' // &
                         trim(input%err_buf) // '.'
    endif
    call printErrMsg(option)
  endif
  
end subroutine InputCheckMandatoryUnits

! ************************************************************************** !

subroutine InputReadAndConvertUnits(input,double_value,internal_units, &
                                    keyword_string,option)
  ! 
  ! Reads units if they exist and returns the units conversion factor.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/16
  ! 
  use Option_module
  use Units_module
  
  implicit none
  
  type(input_type) :: input
  PetscReal :: double_value
  character(len=*) :: internal_units
  character(len=*) :: keyword_string
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: units
  character(len=MAXWORDLENGTH) :: internal_units_word
  character(len=MAXSTRINGLENGTH) :: string

  call InputReadWord(input,option,units,PETSC_TRUE)
  if (input%ierr == 0) then
    if (len_trim(internal_units) < 1) then
      option%io_buffer = 'No internal units provided in &
                         &InputReadAndConvertUnits()'
      call printErrMsg(option)
    endif
    internal_units_word = trim(internal_units)
    double_value = double_value * &
                   UnitsConvertToInternal(units,internal_units_word,option)
  else
    string = trim(keyword_string) // ' units'
    call InputDefaultMsg(input,option,string)
  endif
  
end subroutine InputReadAndConvertUnits

! ************************************************************************** !

subroutine InputPushExternalFile(input,option)
  ! 
  ! Looks up double precision value in database
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/19/14
  ! 
  use Option_module
  
  implicit none
  
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  type(input_type), pointer :: input_child
  
  call InputReadNChars(input,option,string,MAXSTRINGLENGTH,PETSC_TRUE)
  call InputErrorMsg(input,option,'filename','EXTERNAL_FILE')
  input_child => InputCreate(input%fid+1,string,option) 
  input_child%parent => input
  input => input_child

end subroutine InputPushExternalFile

! ************************************************************************** !

function InputPopExternalFile(input)
  ! 
  ! Looks up double precision value in database
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/19/14
  ! 
  
  implicit none
  
  type(input_type), pointer :: input

  PetscBool :: InputPopExternalFile
  type(input_type), pointer :: input_parent
  
  InputPopExternalFile = PETSC_FALSE
  if (associated(input%parent)) then
    input_parent => input%parent
    call InputDestroy(input)
    input => input_parent
    nullify(input_parent)
    InputPopExternalFile = PETSC_TRUE
  endif

end function InputPopExternalFile

! ************************************************************************** !

subroutine InputDbaseDestroy()
  ! 
  ! Destroys the input dbase and members
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/20/14
  ! 

  implicit none
  
  if (associated(dbase)) then
    ! due to circular dependencies, cannot use Utilty_module::DeallocateArray 
    if (associated(dbase%icard)) deallocate(dbase%icard)
    nullify(dbase%icard)
    if (associated(dbase%rcard)) deallocate(dbase%rcard)
    nullify(dbase%rcard)
    if (associated(dbase%ccard)) deallocate(dbase%ccard)
    nullify(dbase%ccard)
    if (associated(dbase%ivalue)) deallocate(dbase%ivalue)
    nullify(dbase%ivalue)
    if (associated(dbase%rvalue)) deallocate(dbase%rvalue)
    nullify(dbase%rvalue)
    if (associated(dbase%cvalue)) deallocate(dbase%cvalue)
    nullify(dbase%cvalue)
    deallocate(dbase)
    nullify(dbase)
  endif
  
end subroutine InputDbaseDestroy

! ************************************************************************** !

subroutine InputDestroy(input)
  ! 
  ! Deallocates an input object
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none
  
  type(input_type), pointer :: input
  
  if (input%fid /= 0) close(input%fid)
  input%fid = 0
  deallocate(input)
  nullify(input)
  
end subroutine InputDestroy

end module Input_Aux_module
