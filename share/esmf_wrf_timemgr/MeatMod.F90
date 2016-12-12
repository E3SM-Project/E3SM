
module MeatMod

#include <ESMF_TimeMgr.inc>

  use ESMF_BaseMod

  implicit none

  private

  public fraction_to_stringi8
  public fraction_to_string

!==============================================================================
contains
!==============================================================================

!==============================================================================

!==============================================================================
! Convert fraction to string with leading sign.
! If fraction simplifies to a whole number or if
! denominator is zero, return empty string.
! INTEGER*8 interface.
SUBROUTINE fraction_to_stringi8( numerator, denominator, frac_str )
  IMPLICIT NONE
  INTEGER(ESMF_KIND_I8), INTENT(IN) :: numerator
  INTEGER(ESMF_KIND_I8), INTENT(IN) :: denominator
  CHARACTER (LEN=*), INTENT(OUT) :: frac_str
  IF ( denominator > 0 ) THEN
    IF ( mod( numerator, denominator ) /= 0 ) THEN
      IF ( numerator > 0 ) THEN
        WRITE(frac_str,FMT="('+',I2.2,'/',I2.2)") abs(numerator), denominator
      ELSE   ! numerator < 0
        WRITE(frac_str,FMT="('-',I2.2,'/',I2.2)") abs(numerator), denominator
      ENDIF
    ELSE   ! includes numerator == 0 case
      frac_str = ''
    ENDIF
  ELSE   ! no-fraction case
    frac_str = ''
  ENDIF
END SUBROUTINE fraction_to_stringi8

!==============================================================================

! Convert fraction to string with leading sign.
! If fraction simplifies to a whole number or if
! denominator is zero, return empty string.
! INTEGER interface.
SUBROUTINE fraction_to_string( numerator, denominator, frac_str )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: numerator
  INTEGER, INTENT(IN) :: denominator
  CHARACTER (LEN=*), INTENT(OUT) :: frac_str
  ! locals
  INTEGER(ESMF_KIND_I8) :: numerator_i8, denominator_i8
  numerator_i8 = INT( numerator, ESMF_KIND_I8 )
  denominator_i8 = INT( denominator, ESMF_KIND_I8 )
  CALL fraction_to_stringi8( numerator_i8, denominator_i8, frac_str )
END SUBROUTINE fraction_to_string

!==============================================================================

end module MeatMod
