module String_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter, public :: STRING_IS_AN_INTEGER = 1
  PetscInt, parameter, public :: STRING_IS_A_DOUBLE = 2
  PetscInt, parameter, public :: STRING_IS_A_WORD = 3

  PetscInt, parameter, public :: STRING_YES = 1
  PetscInt, parameter, public :: STRING_NO = 0
  PetscInt, parameter, public :: STRING_OTHER = UNINITIALIZED_INTEGER

  public :: StringCompare, &
            StringCompareIgnoreCase, &
            StringToUpper, &
            StringToLower, &
            StringReadQuotedWord, &
            StringStartsWithAlpha, &
            StringStartsWith, &
            StringEndsWith, &
            StringAdjustl, &
            StringNull, &
            StringFindEntryInList, &
            StringSplit, &
            StringSwapChar, &
            StringFormatInt, &
            StringFormatDouble, &
            StringIntegerDoubleOrWord, &
            StringYesNoOther
  
  interface StringCompare
    module procedure StringCompare1
    module procedure StringCompare2
  end interface

  interface StringCompareIgnoreCase
    module procedure StringCompareIgnoreCase1
    module procedure StringCompareIgnoreCase2
  end interface

contains

! ************************************************************************** !

PetscBool function StringCompare1(string1,string2,n)
  ! 
  ! compares two strings
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  PetscInt :: i, n
  character(len=n) :: string1, string2
  
  do i=1,n
    if (string1(i:i) /= string2(i:i)) then
      StringCompare1 = PETSC_FALSE
      return
    endif
  enddo

  StringCompare1 = PETSC_TRUE
  return

end function StringCompare1

! ************************************************************************** !

PetscBool function StringCompare2(string1,string2)
  ! 
  ! compares two strings
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/11
  ! 

  implicit none

  PetscInt :: i, length1, length2
  character(len=*) :: string1, string2
  
  length1 = len_trim(string1)
  length2 = len_trim(string2)
  if (length1 /= length2) then
    StringCompare2 = PETSC_FALSE
    return
  endif

  do i=1,length1
    if (string1(i:i) /= string2(i:i)) then
      StringCompare2 = PETSC_FALSE
      return
    endif
  enddo

  StringCompare2 = PETSC_TRUE
  return

end function StringCompare2

! ************************************************************************** !

function StringCompareIgnoreCase1(string1,string2,n)
  ! 
  ! compares two strings
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  PetscInt :: i, n
  character(len=n) :: string1, string2
  
  character(len=n) :: upper1, upper2
  PetscBool :: StringCompareIgnoreCase1
  
  upper1 = string1
  upper2 = string2
  
  call StringToUpper(upper1)
  call StringToUpper(upper2)
  
  do i=1,n
    if (upper1(i:i) /= upper2(i:i)) then
      StringCompareIgnoreCase1 = PETSC_FALSE
      return
    endif
  enddo

  StringCompareIgnoreCase1 = PETSC_TRUE
  return

end function StringCompareIgnoreCase1

! ************************************************************************** !

function StringCompareIgnoreCase2(string1,string2)
  ! 
  ! StringCompare: compares two strings
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  PetscInt :: i, length1, length2
  character(len=*) :: string1, string2
  
  character(len=MAXSTRINGLENGTH) :: upper1, upper2
  PetscBool :: StringCompareIgnoreCase2
  
  length1 = len_trim(string1)
  length2 = len_trim(string2)
  if (length1 /= length2) then
    StringCompareIgnoreCase2 = PETSC_FALSE
    return
  endif

  upper1 = string1
  upper2 = string2
  
  call StringToUpper(upper1)
  call StringToUpper(upper2)
  
  do i=1,length1
    if (upper1(i:i) /= upper2(i:i)) then
      StringCompareIgnoreCase2 = PETSC_FALSE
      return
    endif
  enddo

  StringCompareIgnoreCase2 = PETSC_TRUE
  return

end function StringCompareIgnoreCase2

! ************************************************************************** !

subroutine StringToUpper(string)
  ! 
  ! converts lowercase characters in a card to uppercase
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 
      
  implicit none

  PetscInt :: i
  character(len=*) :: string

  do i=1,len_trim(string)
    if (string(i:i) >= 'a' .and. string(i:i) <= 'z') then
      string(i:i) = achar(iachar(string(i:i)) - 32)
    endif
  enddo

end subroutine StringToUpper

! ************************************************************************** !

subroutine StringToLower(string)
  ! 
  ! converts uppercase characters in a card to lowercase
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 
      
  implicit none

  PetscInt :: i
  character(len=*) :: string

  do i=1,len_trim(string)
    if (string(i:i) >= 'A' .and. string(i:i) <= 'Z') then
      string(i:i) = achar(iachar(string(i:i)) + 32)
    endif
  enddo

end subroutine StringToLower

! ************************************************************************** !

subroutine StringReadQuotedWord(string, name, return_blank_error, ierr)
  ! 
  ! reads and removes a name from a string read from the
  ! database.  "'" are used as delimiters.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/10/08
  ! 

  implicit none

  PetscInt :: i, begins, ends, realends, length
  PetscBool :: return_blank_error ! Return an error for a blank line
                                ! Therefore, a blank line is not acceptable.
  character(len=*) :: string
  character(len=*) :: name
  character(len=1), parameter :: tab = achar(9)
  PetscBool :: openquotefound
  PetscErrorCode :: ierr

  if (ierr /= 0) return

  openquotefound = PETSC_FALSE
  ! Initialize character string to blank.
  length = len_trim(name)
  name(1:length) = repeat(' ',length)

  ierr = 0
  length = len_trim(string)

  ! Remove leading blanks and tabs
  i=1
  do while(string(i:i) == ' ' .or. string(i:i) == tab) 
    i=i+1
  enddo

  if (string(i:i) == "'") then
    openquotefound = PETSC_TRUE
    i=i+1
  endif

  begins=i

  if (openquotefound) then
    do while (string(i:i) /= "'")
      if (i > length) exit
      i=i+1
    enddo
  else
  ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (string(i:i) /= ' ' .and. string(i:i) /= ',' .and. &
              string(i:i) /= tab) 
      i=i+1
    enddo
  endif

  realends = i
  ends=i-1

  ! Avoid copying beyond the end of the word (32 characters).
  if (ends-begins > MAXWORDLENGTH - 1) ends = begins + MAXWORDLENGTH - 1

  ! Copy (ends-begins) characters to 'chars'
  name = string(begins:ends)
  ! Remove chars from string
  string = string(realends+1:)

end subroutine StringReadQuotedWord

! ************************************************************************** !

function StringStartsWithAlpha(string)
  ! 
  ! Determines whether a string starts with an alpha char
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/07/10
  ! 
      
  implicit none

  character(len=*) :: string

  PetscBool :: StringStartsWithAlpha

  string = adjustl(string)

  if ((string(1:1) >= 'a' .and. string(1:1) <= 'z') .or. &
      (string(1:1) >= 'A' .and. string(1:1) <= 'Z')) then
    StringStartsWithAlpha = PETSC_TRUE
  else
    StringStartsWithAlpha = PETSC_FALSE
  endif

end function StringStartsWithAlpha

! ************************************************************************** !

function StringStartsWith(string,string2)
  ! 
  ! Determines whether a string starts with characters
  ! identical to another string
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/16/12
  ! 
      
  implicit none

  character(len=*) :: string
  character(len=*) :: string2

  PetscBool :: StringStartsWith
  
  
  PetscInt :: length, i

  length = min(len_trim(string),len_trim(string2))
  
  do i = 1, length
    if (string(i:i) /= string2(i:i)) then
      StringStartsWith = PETSC_FALSE
      return
    endif
  enddo
  
  StringStartsWith = PETSC_TRUE

end function StringStartsWith

! ************************************************************************** !

function StringEndsWith(string,string2)
  ! 
  ! Determines whether a string ends with characters identical to another 
  ! string
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/16/12
  ! 
      
  implicit none

  character(len=*) :: string
  character(len=*) :: string2

  PetscBool :: StringEndsWith
  
  PetscInt :: search_length, i, i1, i2, len1, len2

  len1 = len_trim(string)
  len2 = len_trim(string2)
  search_length = min(len1,len2)
  
  do i = 1, search_length
    ! search backward
    i1 = len1+1-i
    i2 = len2+1-i
    if (string(i1:i1) /= string2(i2:i2)) then
      StringEndsWith = PETSC_FALSE
      return
    endif
  enddo
  
  StringEndsWith = PETSC_TRUE

end function StringEndsWith

! ************************************************************************** !

subroutine StringAdjustl(string)
  ! 
  ! Left adjusts a string by removing leading spaces and tabs.
  ! This subroutine is needed because the adjustl() Fortran 90
  ! intrinsic will not remove leading tabs.
  ! 
  ! Author: Richard Tran Mills
  ! Date: 9/21/2010
  ! 

  implicit none

  character(len=*) :: string
  
  PetscInt :: i
  PetscInt :: string_length
  character(len=1), parameter :: tab = achar(9)

  ! We have to manually convert any leading tabs into spaces, as the 
  ! adjustl() intrinsic does not eliminate leading tabs.
  i=1
  string_length = len_trim(string)
  do while((string(i:i) == ' ' .or. string(i:i) == tab) .and. &
           i <= string_length)
    if (string(i:i) == tab) string(i:i) = ' '
    i=i+1
  enddo

  ! adjustl() will do what we want, now that tabs are removed.
  string = adjustl(string) 

end subroutine StringAdjustl

! ************************************************************************** !

function StringNull(string)
  ! 
  ! Returns PETSC_TRUE if a string is blank
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/07/10
  ! 
      
  implicit none

  character(len=*) :: string

  PetscBool :: StringNull
  PetscInt :: length

  length = len_trim(adjustl(string))
  if (length > 0) then
    StringNull = PETSC_FALSE
  else
    StringNull = PETSC_TRUE
  endif

end function StringNull

! ************************************************************************** !

function StringFindEntryInList(string,string_array)
  ! 
  ! Returns the index of a string if found in a list
  ! of strings
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  ! 
      
  implicit none

  character(len=*) :: string
  character(len=*) :: string_array(:)

  PetscInt :: StringFindEntryInList
  PetscInt :: i

  StringFindEntryInList = 0
  
  do i = 1, size(string_array)
    if (StringCompare(string,string_array(i))) then
      StringFindEntryInList = i
      exit
    endif
  enddo
  
end function StringFindEntryInList

! ************************************************************************** !

subroutine StringSwapChar(string,char_in,char_out)
  ! 
  ! Swaps a character from a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 

  implicit none
 
  character(len=*) :: string
  character(len=1) :: char_in
  character(len=1) :: char_out
 
  PetscInt :: i
 
  do i=1, len_trim(string)
   if (string(i:i) == char_in(1:1)) string(i:i) = char_out(1:1)
  enddo
 
end subroutine StringSwapChar

! ************************************************************************** !

function StringSplit(string,chars)
  ! 
  ! Splits a string based on a set of chars
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/28/13
  ! 
      
  implicit none

  character(len=*) :: string
  character(len=*) :: chars

  character(len=MAXSTRINGLENGTH), pointer :: strings(:), StringSplit(:)
  
  character(len=MAXSTRINGLENGTH) :: string1
  PetscInt :: i, icount, istart, iend, length, length_chars
  PetscInt :: last_index
  
  nullify(StringSplit)
  
  ! determine how many delimiting block in string
  length = len_trim(string)
  ! do not use len_trim for chars.  All characters including the blank and 
  ! trailing blanks (spaces) should be accounted for.
  length_chars = len(chars)
  icount = 0
  last_index = 1
  iend = length-length_chars+1
  do i = 1, iend
    string1 = string(i:i+length_chars-1)
    if (StringCompare(string1,chars,length_chars)) then
      last_index = i+1
      icount = icount + 1
    endif
  enddo
  
  ! check for characters after last delimiter; add a string if they exist
  if (last_index <= length) then
    if (.not.StringNull(string(last_index:))) then
      icount = icount + 1
    endif
  endif
  
  if (icount == 0) return
  
  ! allocate strings
  allocate(strings(icount))
  strings = ''

  ! split string into strings
  istart = 1
  icount = 0
  iend = length-length_chars+1
  i = 1
  do 
    if (i > iend) exit
    string1 = string(i:i+length_chars-1)
    if (StringCompare(string1,chars,length_chars)) then
      icount = icount + 1
      strings(icount) = adjustl(string(istart:i-1))
      i = i + length_chars
      istart = i
    else
      i = i + 1
    endif
  enddo 
  
  ! add remaining string
  if (icount < size(strings)) then
    icount = icount + 1
    strings(icount) = adjustl(string(istart:))
  endif  
  
  StringSplit => strings
  
end function StringSplit

! ************************************************************************** !

function StringFormatInt(int_value)
  ! 
  ! Writes a integer to a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/13/12
  ! 

  implicit none
  
  PetscInt :: int_value
  
  character(len=MAXWORDLENGTH) :: StringFormatInt

  write(StringFormatInt,'(1i12)') int_value
  
  StringFormatInt = adjustl(StringFormatInt)
  
end function StringFormatInt

! ************************************************************************** !

function StringFormatDouble(real_value)
  ! 
  ! Writes a double or real to a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/13/12
  ! 

  implicit none
  
  PetscReal :: real_value
  
  character(len=MAXWORDLENGTH) :: StringFormatDouble

  write(StringFormatDouble,'(1es13.5)') real_value
  
  StringFormatDouble = adjustl(StringFormatDouble)
  
end function StringFormatDouble

! ************************************************************************** !

function StringIntegerDoubleOrWord(string_in)
  ! 
  ! Returns whether a value read from a string is a double, integer, or word
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/13/12
  ! 

  implicit none
  
  character(len=*) :: string_in

  PetscInt :: StringIntegerDoubleOrWord

  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: d
  PetscInt :: i
  PetscBool :: double_syntax_found
  PetscErrorCode :: ierr

  string = trim(string_in)

  StringIntegerDoubleOrWord = -999
  ierr = 0
  double_syntax_found = (index(string,'.') > 0 .or. &
      index(string,'d') > 0 .or. index(string,'D') > 0 .or. &
      index(string,'e') > 0 .or. index(string,'E') > 0) 
  read(string,*,iostat=ierr) i
  if (ierr == 0) then
    ! the Intel compiler does not alway catch the misread of a double to an 
    ! integer
    if (double_syntax_found) then
      StringIntegerDoubleOrWord = STRING_IS_A_DOUBLE
      return
    endif
    StringIntegerDoubleOrWord = STRING_IS_AN_INTEGER
    return
  endif
  ierr = 0
  read(string,*,iostat=ierr) d
  if (ierr == 0) then
    StringIntegerDoubleOrWord = STRING_IS_A_DOUBLE
    return
  endif
  if (len_trim(string) > 0) StringIntegerDoubleOrWord = STRING_IS_A_WORD
  
end function StringIntegerDoubleOrWord

! ************************************************************************** !

function StringYesNoOther(string)
  ! 
  ! Returns PETSC_TRUE if a string is blank
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/07/10
  ! 
      
  implicit none

  character(len=*) :: string

  character(len=MAXSTRINGLENGTH) :: string2

  PetscInt :: StringYesNoOther

  string2 = adjustl(string)
  call StringToUpper(string2)

  StringYesNoOther = STRING_OTHER
  if (len_trim(string2) == 3 .and. StringStartsWith(string2,'YES')) then
    StringYesNoOther = STRING_YES
  elseif (len_trim(string2) == 2 .and. StringStartsWith(string2,'NO')) then
    StringYesNoOther = STRING_NO
  endif

end function StringYesNoOther

end module String_module
