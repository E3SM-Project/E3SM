! FILE: c_interface_module.f
! PURPOSE: Supplement ISO-C-Binding to provide type aliases and interfaces
!          to common ISO-C string functions to aid working with strings.
! AUTHOR: Joseph M. Krahn
! STATUS: Still in development. Reasonably complete, but somewhat limited testing.
!
! The idea is to provide type aliases for all ISO-C types, so that the
! Fortran interface code more explicitly defines the actual C interface.
! This should be updated to support F2008 variable-length allocatable
! strings.
!
! Entity names all have the "C_" prefix, as with ISO-C-Binding, with a
! few exceptions.
! 
module C_interface_mod
  use, intrinsic :: ISO_C_Binding, &
  ! C type aliases for pointer derived types:
      C_ptr => C_ptr , &
      C_char_ptr => C_ptr, &
      C_const_char_ptr => C_ptr, &
      C_void_ptr => C_ptr, &
      C_const_void_ptr => C_ptr

  implicit none

!----------------------------------------------------------------------------
! C type aliases for intrinsic type KIND parameters:

! NOTE: a C enum may not always be a standard C int
  integer, parameter :: C_enum = C_int

! Defining off_t is difficult, because it may depend on "LARGEFILE" selection.
!  integer, parameter :: C_off_t = ??

! C string terminator alais using the 3-letter ASCII name.
! The C_ prefix is not used because it is just an ASCII character.
  character(len=1,kind=C_char), parameter :: NUL = C_NULL_char

! NOTE: In C, "char" is distinct from "signed char", unlike integers.
! The plain "char" type is specific for text/string values, whereas
! "signed char" should indicate 1-byte integer data.
!
! Most ISO-C systems have wide chars "wchar_t", but Fortran compilers
! have limited support for different character kinds. UTF encoding
! adds more complexity. This should be updated as Fortran compilers
! include support for more character types.
! 

! Fortran does not (yet) support unsigned types.
  integer, parameter :: &
      C_unsigned = C_int, &
      C_unsigned_short = C_short, &
      C_unsigned_long = C_long, &
      C_unsigned_long_long = C_long_long, &
      C_unsigned_char = C_signed_char, &
      C_ssize_t = C_size_t, &
      C_uint8_t = C_int8_t, &
      C_uint16_t = C_int16_t, &
      C_uint32_t = C_int32_t, &
      C_uint64_t = C_int64_t, &
      C_uint_least8_t = C_int_least8_t, &
      C_uint_least16_t = C_int_least16_t, &
      C_uint_least32_t = C_int_least32_t, &
      C_uint_least64_t = C_int_least64_t, &
      C_uint_fast8_t = C_int_fast8_t, &
      C_uint_fast16_t = C_int_fast16_t, &
      C_uint_fast32_t = C_int_fast32_t, &
      C_uint_fast64_t = C_int_fast64_t, &
      C_uintmax_t = C_intmax_t
! Note: ptrdiff_t cannot be reliably defined from other types.
! When practical, it is larger than a pointer because it benefits
! from the full unsigned range in both positive and negative directions.

! Integer versions including 'int', where the 'int' is optional:
  integer, parameter :: &
      C_short_int = C_short, &
      C_long_int = C_long, &
      C_long_long_int = C_long_long, &
      C_unsigned_int = C_unsigned, &
      C_unsigned_short_int = C_short, &
      C_unsigned_long_int = C_long, &
      C_unsigned_long_long_int = C_long_long

  interface C_F_string
    module procedure C_F_string_ptr
    module procedure C_F_string_chars
  end interface C_F_string

  interface F_C_string
    module procedure F_C_string_ptr
    module procedure F_C_string_chars
  end interface F_C_string

!=======================================================================
! Some useful ISO C library string functions from <string.h>
! These are based on GCC header sections marked as NAMESPACE_STD
  interface

! Copy N bytes of SRC to DEST, no aliasing or overlapping allowed.
!extern void *memcpy (void *dest, const void *src, size_t n);
    function C_memcpy(dest, src, n) result(result) bind(C,name="memcpy")
      import C_void_ptr, C_size_t
      type(C_void_ptr) :: result
      type(C_void_ptr), value, intent(in) :: dest ! target=intent(out)
      type(C_void_ptr), value, intent(in) :: src  ! target=intent(in)
      integer(C_size_t), value, intent(in) :: n
    end function C_memcpy

! Copy N bytes of SRC to DEST, guaranteeing correct behavior for overlapping strings.
!extern void *memmove (void *dest, const void *src, size_t n)
    function C_memmove(dest, src, n) result(result) bind(C,name="memmove")
      import C_void_ptr, C_size_t
      type(C_void_ptr) :: result
      type(C_void_ptr), value, intent(in) :: dest ! target=intent(out)
      type(C_void_ptr), value, intent(in) :: src
      integer(C_size_t), value, intent(in) :: n
    end function C_memmove

! Set N bytes of S to C.
!extern void *memset (void *s, int c, size_t n)
    function C_memset(s, c, n) result(result) bind(C,name="memset")
      import C_void_ptr, C_int, C_size_t
      type(C_void_ptr) :: result
      type(C_void_ptr), value, intent(in) :: s ! target=intent(out)
      integer(C_int), value, intent(in) :: c
      integer(C_size_t), value, intent(in) :: n
    end function C_memset

! Compare N bytes of S1 and S2.
!extern int memcmp (const void *s1, const void *s2, size_t n)
    pure &
    function C_memcmp(s1, s2, n) result(result) bind(C,name="memcmp")
      import C_int, C_void_ptr, C_size_t
      integer(C_int) :: result
      type(C_void_ptr), value, intent(in) :: s1
      type(C_void_ptr), value, intent(in) :: s2
      integer(C_size_t), value, intent(in) :: n
    end function C_memcmp

! Search N bytes of S for C.
!extern void *memchr (const void *s, int c, size_t n)
    pure &
    function C_memchr(s, c, n) result(result) bind(C,name="memchr")
      import C_void_ptr, C_int, C_size_t
      type(C_void_ptr) :: result
      type(C_void_ptr), value, intent(in) :: s
      integer(C_int), value, intent(in) :: c
      integer(C_size_t), value, intent(in) :: n
    end function C_memchr

! Copy SRC to DEST.
!extern char *strcpy (char *dest, const char *src)
    function C_strcpy(dest, src) result(result) bind(C,name="strcpy")
      import C_char_ptr, C_size_t
      type(C_char_ptr) :: result
      type(C_char_ptr), value, intent(in) :: dest ! target=intent(out)
      type(C_char_ptr), value, intent(in) :: src
    end function C_strcpy

! Copy no more than N characters of SRC to DEST.
!extern char *strncpy (char *dest, const char *src, size_t n)
    function C_strncpy(dest, src, n) result(result) bind(C,name="strncpy")
      import C_char_ptr, C_size_t
      type(C_char_ptr) :: result
      type(C_char_ptr), value, intent(in) :: dest ! target=intent(out)
      type(C_char_ptr), value, intent(in) :: src
      integer(C_size_t), value, intent(in) :: n
    end function C_strncpy

! Append SRC onto DEST.
!extern char *strcat (char *dest, const char *src)
    function C_strcat(dest, src) result(result) bind(C,name="strcat")
      import C_char_ptr, C_size_t
      type(C_char_ptr) :: result
      type(C_char_ptr), value, intent(in) :: dest ! target=intent(out)
      type(C_char_ptr), value, intent(in) :: src
    end function C_strcat

! Append no more than N characters from SRC onto DEST.
!extern char *strncat (char *dest, const char *src, size_t n)
    function C_strncat(dest, src, n) result(result) bind(C,name="strncat")
      import C_char_ptr, C_size_t
      type(C_char_ptr) :: result
      type(C_char_ptr), value, intent(in) :: dest ! target=intent(out)
      type(C_char_ptr), value, intent(in) :: src
      integer(C_size_t), value, intent(in) :: n
    end function C_strncat

! Compare S1 and S2.
!extern int strcmp (const char *s1, const char *s2)
    pure &
    function C_strcmp(s1, s2) result(result) bind(C,name="strcmp")
      import C_int, C_char_ptr, C_size_t
      integer(C_int) :: result
      type(C_char_ptr), value, intent(in) :: s1
      type(C_char_ptr), value, intent(in) :: s2
    end function C_strcmp

! Compare N characters of S1 and S2.
!extern int strncmp (const char *s1, const char *s2, size_t n)
    pure &
    function C_strncmp(s1, s2, n) result(result) bind(C,name="strncmp")
      import C_int, C_char_ptr, C_size_t
      integer(C_int) :: result
      type(C_char_ptr), value, intent(in) :: s1
      type(C_char_ptr), value, intent(in) :: s2
      integer(C_size_t), value, intent(in) :: n
    end function C_strncmp

! Return the length of S.
!extern size_t strlen (const char *s)
    pure &
    function C_strlen(s) result(result) bind(C,name="strlen")
      import C_char_ptr, C_size_t
      integer(C_size_t) :: result
      type(C_char_ptr), value, intent(in) :: s  !character(len=*), intent(in)
    end function C_strlen

  end interface

! End of <string.h>
!=========================================================================
! Standard ISO-C malloc routines:
  interface

    ! void *calloc(size_t nmemb, size_t size);
    type(C_void_ptr) &
    function C_calloc(nmemb, size) bind(C,name="calloc")
      import C_void_ptr, C_size_t
      integer(C_size_t), value, intent(in) :: nmemb, size
    end function C_calloc

    ! void *malloc(size_t size);
    type(C_void_ptr) &
    function C_malloc(size) bind(C,name="malloc")
      import C_void_ptr, C_size_t
      integer(C_size_t), value, intent(in) :: size
    end function C_malloc

    ! void free(void *ptr);
    subroutine C_free(ptr) bind(C,name="free")
      import C_void_ptr
      type(C_void_ptr), value, intent(in) :: ptr
    end subroutine C_free

    ! void *realloc(void *ptr, size_t size);
    type(C_void_ptr) &
    function C_realloc(ptr,size) bind(C,name="realloc")
      import C_void_ptr, C_size_t
      type(C_void_ptr), value, intent(in) :: ptr
      integer(C_size_t), value, intent(in) :: size
    end function C_realloc

  end interface

  interface assignment(=)
    module procedure F_string_assign_C_string
  end interface assignment(=)

!==========================================================================

contains

  ! HACK: For some reason, C_associated was not defined as pure. 
  pure logical &
  function C_associated_pure(ptr) result(associated)
    type(C_ptr), intent(in) :: ptr
    integer(C_intptr_t) :: iptr
    iptr = transfer(ptr,iptr)
    associated = (iptr /= 0)
  end function C_associated_pure

! Set a fixed-length Fortran string to the value of a C string.
  subroutine F_string_assign_C_string(F_string, C_string)
    character(len=*), intent(out) :: F_string
    type(C_ptr), intent(in) :: C_string
    character(len=1,kind=C_char), pointer :: p_chars(:)
    integer :: i
    if (.not. C_associated(C_string) ) then
      F_string = ' '
    else
      call C_F_pointer(C_string,p_chars,[huge(0)])
      i=1
      do while(p_chars(i)/=NUL .and. i<=len(F_string))
        F_string(i:i) = p_chars(i)
        i=i+1
      end do
      if (i<len(F_string)) F_string(i:) = ' '
    end if
  end subroutine F_string_assign_C_string

! Copy a C string, passed by pointer, to a Fortran string.
! If the C pointer is NULL, the Fortran string is blanked.
! C_string must be NUL terminated, or at least as long as F_string.
! If C_string is longer, it is truncated. Otherwise, F_string is
! blank-padded at the end.
  subroutine C_F_string_ptr(C_string, F_string)
    type(C_ptr), intent(in) :: C_string
    character(len=*), intent(out) :: F_string
    character(len=1,kind=C_char), dimension(:), pointer :: p_chars
    integer :: i
    if (.not. C_associated(C_string)) then
      F_string = ' '
    else
      call C_F_pointer(C_string,p_chars,[huge(0)])
      i=1
      do while(p_chars(i)/=NUL .and. i<=len(F_string))
        F_string(i:i) = p_chars(i)
        i=i+1
      end do
      if (i<len(F_string)) F_string(i:) = ' '
    end if
  end subroutine C_F_string_ptr

! Copy a C string, passed as a char-array reference, to a Fortran string.
  subroutine C_F_string_chars(C_string, F_string)
    character(len=1,kind=C_char), intent(in) :: C_string(*)
    character(len=*), intent(out) :: F_string
    integer :: i
    i=1
    do while(C_string(i)/=NUL .and. i<=len(F_string))
      F_string(i:i) = C_string(i)
      i=i+1
    end do
    if (i<len(F_string)) F_string(i:) = ' '
  end subroutine C_F_string_chars

! Copy a Fortran string to an allocated C string pointer.
! If the C pointer is NULL, no action is taken. (Maybe auto allocate via libc call?)
! If the length is not passed, the C string must be at least: len(F_string)+1
! If the length is passed and F_string is too long, it is truncated.
  subroutine F_C_string_ptr(F_string, C_string, C_string_len)
    character(len=*), intent(in) :: F_string
    type(C_ptr), intent(in) :: C_string ! target = intent(out)
    integer, intent(in), optional :: C_string_len  ! Max string length,
                                                   ! INCLUDING THE TERMINAL NUL
    character(len=1,kind=C_char), dimension(:), pointer :: p_chars
    integer :: i, strlen
    strlen = len(F_string)
    if (present(C_string_len)) then
      if (C_string_len <= 0) return
      strlen = min(strlen,C_string_len-1)
    end if
    if (.not. C_associated(C_string)) then
      return
    end if
    call C_F_pointer(C_string,p_chars,[strlen+1])
    forall (i=1:strlen)
      p_chars(i) = F_string(i:i)
    end forall
    p_chars(strlen+1) = NUL
  end subroutine F_C_string_ptr

  pure &
  function C_strlen_safe(s) result(length)
    integer(C_size_t) :: length
    type(C_char_ptr), value, intent(in) :: s
    if (.not. C_associated_pure(s)) then
      length = 0
    else
      length = C_strlen(s)
    end if
  end function C_strlen_safe

  function C_string_value(C_string) result(F_string)
    type(C_ptr), intent(in) :: C_string
    character(len=C_strlen_safe(C_string)) :: F_string
    character(len=1,kind=C_char), dimension(:), pointer :: p_chars
    integer :: i, length
    length = len(F_string)
    if (length/=0) then
      call C_F_pointer(C_string,p_chars,[length])
      forall (i=1:length)
        F_string(i:i) = p_chars(i)
      end forall
    end if
  end function C_string_value

! Copy a Fortran string to a C string passed by char-array reference.
! If the length is not passed, the C string must be at least: len(F_string)+1
! If the length is passed and F_string is too long, it is truncated.
  subroutine F_C_string_chars(F_string, C_string, C_string_len)
    character(len=*), intent(in) :: F_string
    character(len=1,kind=C_char), dimension(*), intent(out) :: C_string
    integer, intent(in), optional :: C_string_len  ! Max string length,
                                                   ! INCLUDING THE TERMINAL NUL
    integer :: i, strlen
    strlen = len(F_string)
    if (present(C_string_len)) then
      if (C_string_len <= 0) return
      strlen = min(strlen,C_string_len-1)
    end if
    forall (i=1:strlen)
      C_string(i) = F_string(i:i)
    end forall
    C_string(strlen+1) = NUL
  end subroutine F_C_string_chars

! NOTE: Strings allocated here must be freed by the
! C library, such as via C_free() or C_string_free(),
  type(C_ptr) &
  function F_C_string_dup(F_string,length) result(C_string)
    character(len=*), intent(in) :: F_string
    integer, intent(in), optional :: length
    character(len=1,kind=C_char), pointer :: C_string_ptr(:)
    integer :: i
    integer(C_size_t) :: strlen
    if (present(length)) then
      strlen = length
    else
      strlen = len(F_string)
    end if
    if (strlen <= 0) then
      C_string = C_NULL_ptr
    else
      C_string = C_malloc(strlen+1)
      if (C_associated(C_string)) then
        call C_F_pointer(C_string,C_string_ptr,[strlen+1])
        forall (i=1:strlen)
          C_string_ptr(i) = F_string(i:i)
        end forall
        C_string_ptr(strlen+1) = NUL
      end if
    end if
  end function F_C_string_dup

! NOTE: Strings allocated here must be freed by the
! C library, such as via C_free() or C_string_free(),
  type(C_ptr) &
  function C_string_alloc(length) result(C_string)
    integer(C_size_t), intent(in) :: length
    character(len=1,kind=C_char), pointer :: C_charptr
    C_string = C_malloc(length+1)
    if (C_associated(C_string)) then
      call C_F_pointer(C_string,C_charptr)
      C_charptr = NUL
    end if
  end function C_string_alloc

  subroutine C_string_free(string)
    type(C_ptr), intent(inout) :: string
    if (C_associated(string)) then
      call C_free(string)
      string = C_NULL_ptr
    end if
  end subroutine C_string_free

end module C_interface_mod
