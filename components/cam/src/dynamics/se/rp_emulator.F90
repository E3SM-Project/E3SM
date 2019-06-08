! Copyright 2015-2016 Andrew Dawson, Peter Dueben
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

MODULE rp_emulator
! A reduced-precision emulator.
!
! The `rpe_var` type is a simple container for a double precision
! floating point value.
!
    IMPLICIT NONE

    ! All definitions are private by default.
    PRIVATE

!-----------------------------------------------------------------------
! Module parameters and variables:
!-----------------------------------------------------------------------

    !: The Fortran kind of a single-precision and double-precision
    !: floating-point number.
    INTEGER, PARAMETER, PUBLIC :: RPE_SINGLE_KIND = selected_real_kind(6), &
                                  RPE_DOUBLE_KIND = selected_real_kind(12)

    !: The Fortran kind of the real data type used by the emulator
    !: (usually 64-bit double-precision).
    INTEGER, PARAMETER, PUBLIC :: RPE_REAL_KIND = RPE_DOUBLE_KIND
    !: The Fortran kind of an alternate real data type (usually 32-bit
    !: single-precision), should be single-precision if RPE_REAL_KIND
    !: is double-precision or double-precision if RPE_REAL_KIND is
    !: single-precision.
    INTEGER, PARAMETER, PUBLIC :: RPE_ALTERNATE_KIND = RPE_SINGLE_KIND

    !: Logical flag for turning the emulator on/off.
    LOGICAL, PUBLIC :: RPE_ACTIVE = .TRUE.

    !: The default number of bits to use in the reduced-precision significand.
    INTEGER, PUBLIC :: RPE_DEFAULT_SBITS = 52

    !: Logical flag for determining if IEEE half-precision rules should
    !: be used when operating on values with 10 bits in the significand.
    LOGICAL, PUBLIC :: RPE_IEEE_HALF = .FALSE.

    !: An internal value used to represent the case where a reduced-precision
    !: number has no specified precision yet.
    INTEGER, PARAMETER, PRIVATE :: RPE_SBITS_UNSPECIFIED = -1

!-----------------------------------------------------------------------
! Module derived-type definitions:
!-----------------------------------------------------------------------

    PUBLIC :: rpe_var
    TYPE :: rpe_var
    ! A reduced-precision floating-point number.
    !
    ! This type is a container for a floating-point number which is
    ! operated on in reduced precision.
    !
        INTEGER :: sbits = RPE_SBITS_UNSPECIFIED
        REAL(KIND=RPE_REAL_KIND) :: val
    END TYPE

    ! Create a public interface for constructing literal reduced
    ! precision values (rpe_var instances).
    PUBLIC rpe_literal
    INTERFACE rpe_literal
        MODULE PROCEDURE rpe_literal_real
        MODULE PROCEDURE rpe_literal_alternate
        MODULE PROCEDURE rpe_literal_integer
        MODULE PROCEDURE rpe_literal_long
    END INTERFACE

    ! Make the core emulator routines importable.
    PUBLIC :: apply_truncation, significand_bits

    PUBLIC ASSIGNMENT(=)
    INTERFACE ASSIGNMENT(=)
    ! Interfaces for the assignment operator with `rpe_type` instances.
    !
        MODULE PROCEDURE assign_rpe_rpe
        MODULE PROCEDURE assign_rpe_real
        MODULE PROCEDURE assign_rpe_alternate
        MODULE PROCEDURE assign_rpe_integer
        MODULE PROCEDURE assign_rpe_long
        MODULE PROCEDURE assign_real_rpe
        MODULE PROCEDURE assign_alternate_rpe
        MODULE PROCEDURE assign_integer_rpe
        MODULE PROCEDURE assign_long_rpe
    END INTERFACE

!-----------------------------------------------------------------------
! Interfaces for overloaded operators (external):
!-----------------------------------------------------------------------

#include "./include/interface_operators.i"

!-----------------------------------------------------------------------
! Interfaces for overloaded intrinsic functions (external):
!-----------------------------------------------------------------------

#include "./include/interface_intrinsics.i"

!-----------------------------------------------------------------------
! Interfaces for other extensions (external):
!-----------------------------------------------------------------------

#include "./include/interface_extras.i"

CONTAINS

!-----------------------------------------------------------------------
! Core emulator procedures:
!-----------------------------------------------------------------------

    ELEMENTAL SUBROUTINE apply_truncation (x)
    ! Reduce the precision of a `rpe_type` instance.
    !
    ! Truncates the given floating-point number significand to the
    ! number of bits defined by the `sbits` member of the number. If the
    ! `sbits` attribute is not set it will truncate to the number of
    ! bits specified by the current value of `RPE_DEFAULT_SBITS`.
    !
    ! If the module variable RPE_ACTIVE is false this subroutine returns
    ! the unaltered input value, it only performs the bit truncation if
    ! RPE_ACTIVE is true.
    !
    ! Argument:
    !
    ! * x: type(rpe_type) [input/output]
    !     The `rpe_type` instance to truncate.
    !
        TYPE(rpe_var), INTENT(INOUT) :: x
        REAL(KIND=RPE_DOUBLE_KIND)   :: y
        INTEGER :: n
        IF (RPE_ACTIVE) THEN
            ! Cast the input to a double-precision value.
            y = REAL(x%val, RPE_DOUBLE_KIND)
            IF (x%sbits == RPE_SBITS_UNSPECIFIED) THEN
                ! If the input does not have a specified precision then assume
                ! the default precision. This is does not fix the precision of
                ! the input variable, it will still use whatever is specified
                ! as the default, even if that changes later.
                n = RPE_DEFAULT_SBITS
            ELSE
                n = x%sbits
            END IF
            x%val = truncate_significand(y, n)
        END IF
    END SUBROUTINE apply_truncation

    ELEMENTAL FUNCTION truncate_significand (x, n) RESULT (t)
    ! Truncate the significand of a double precision floating point
    ! number to a specified number of bits.
    !
    ! Arguments:
    !
    ! * x: real(kind=RPE_DOUBLE_KIND) [input]
    !     The double precision number to truncate.
    !
    ! * n: integer [input]
    !     The number of bits to truncate the significand to.
    !
    ! Returns:
    !
    ! * t: real(kind=RPE_DOUBLE_KIND)
    !     A double precision number representing `x` truncated to `n`
    !     bits in the significand.
    !
        REAL(KIND=RPE_DOUBLE_KIND), INTENT(IN) :: x
        INTEGER,                    INTENT(IN) :: n
        REAL(KIND=RPE_DOUBLE_KIND) :: t
        INTEGER                    :: lmtb
        INTEGER(KIND=8), PARAMETER :: two = 2
        INTEGER(KIND=8), PARAMETER :: zero_bits = 0
        INTEGER(KIND=8)            :: bits
        ! The left-most truncated bit is the last bit that will be truncated
        ! (counting from 0 at the right-most bit). Double precision values
        ! have 52 bits in their significand.
        lmtb = 52 - n - 1
        IF (lmtb >= 0) THEN
            ! Copy the double-precision bit representation of the input
            ! into an integer so it can be manipulated:
            bits = TRANSFER(x, bits)
            ! Round the number up first if required according to IEEE 754
            ! specifications.
            IF (BTEST(bits, lmtb)) THEN
                IF (IAND(bits, two ** (lmtb + 1) - 1) == two ** lmtb) THEN
                    ! We are truncating a number half-way between two
                    ! representations so we must round to the nearest even
                    ! representation.
                    IF (BTEST(bits, lmtb + 1)) THEN
                        bits = bits + two ** (lmtb + 1)
                    END IF
                ELSE
                    ! The left-most truncated bit is set and we are not
                    ! half-way between two representations so we need to
                    ! round to the nearest representation.
                    bits = bits + two ** (lmtb + 1)
                END IF
            END IF
            ! Move rounding_bit + 1 bits from the number zero (all bits
            ! set to zero) into the target to truncate at the given
            ! number of bits.
            CALL MVBITS (zero_bits, 0, lmtb + 1, bits, 0)
            t = TRANSFER(bits, t)
            ! Special case for IEEE half-precision representation.
            IF (n == 10 .AND. RPE_IEEE_HALF) THEN
                t = adjust_ieee_half(t)
            END IF
        ELSE
            t = x
        END IF
    END FUNCTION truncate_significand

    ELEMENTAL FUNCTION adjust_ieee_half (x) RESULT (y)
    ! Adjust a floating-point number according to the IEEE half-precision
    ! specification by ensuring that numbers outside the range of
    ! half-precision are either overflowed or denormalized.
    !
    ! Arguments:
    !
    ! * x: real(kind=RPE_DOUBLE_KIND) [input]
    !     The floating-point number to adjust.
    !
    ! Returns:
    !
    ! * y: real(kind=RPE_DOUBLE_KIND) [output]
    !     The adjusted floating-point value.
    !
        REAL(KIND=RPE_DOUBLE_KIND), INTENT(IN)  :: x
        REAL(KIND=RPE_DOUBLE_KIND)            :: y
        REAL(KIND=RPE_DOUBLE_KIND), PARAMETER :: two = 2.0_RPE_DOUBLE_KIND
        REAL(KIND=RPE_DOUBLE_KIND)            :: sx, d1, d2
        IF (ABS(x) > 65504.0_RPE_DOUBLE_KIND) THEN
        ! Handle half-precision overflows.
            sx = SIGN(1.0_RPE_DOUBLE_KIND, x)
            d1 = HUGE(d1) * sx
            d2 = HUGE(d2) * sx
            ! This will deliberately cause a floating-point overflow exception
            ! and yield an infinity value with the same sign as the input if
            ! the exception is not handled.
            y = d1 + d2
        ELSE IF (ABS(x) < two ** (-14)) THEN
        ! Handle half-precision subnormal values.
            d1 = two ** (-24)    ! step size for subnormal values
            d2 = MOD(ABS(x), d1)
            ! The new value is the old value rounded to the nearest subnormal
            ! step interval in the direction of zero.
            y = (ABS(x) - d2) * SIGN(1.0_RPE_DOUBLE_KIND, x)
            ! If the rounding should have gone away from zero then correct for
            ! it afterwards.
            IF (ABS(d2) > two ** (-25)) THEN
                y = y + d1 * SIGN(1.0_RPE_DOUBLE_KIND, y)
            END IF
        ELSE
        ! If the value is in range then just return it.
            y = x
        END IF
    END FUNCTION adjust_ieee_half

    ELEMENTAL FUNCTION significand_bits (x) RESULT (z)
    ! Retrieve the number of bits in a floating point significand.
    !
    ! This returns actual values for inputs of type `rpe_type` or `real`
    ! and 0 for anything else. This function is usually used to find the
    ! highest precision level involved in a floating-point calculation.
    !
    ! Arguments:
    !
    ! * x: class(*) [input]
    !     A scalar input of any type.
    !
    ! Returns:
    !
    ! * z: integer [output]
    !     The number of bits in the significand of the input floating-point
    !     value, or 0 if the input was not a floating-point value.
    !
        CLASS(*), INTENT(IN) :: x
        INTEGER :: z
        SELECT TYPE (x)
        TYPE IS (REAL(KIND=RPE_DOUBLE_KIND))
            z = 52
        TYPE IS (rpe_var)
            IF (x%sbits == RPE_SBITS_UNSPECIFIED) THEN
                z = RPE_DEFAULT_SBITS
            ELSE
                z = x%sbits
            END IF
        TYPE IS (REAL(KIND=RPE_SINGLE_KIND))
            z = 23
        CLASS DEFAULT
            z = 0
        END SELECT
    END FUNCTION significand_bits

    FUNCTION rpe_literal_real (x, n) RESULT (z)
    ! Create an `rpe_var` instance from a real literal.
    !
    ! Arguments:
    !
    ! * x: real(kind=RPE_REAL_KIND) [input]
    !     The literal to transform to a reduced precision `rpe_var` instance.
    !
    ! * n: integer [input, optional]
    !     The number of bits in the significand of the resulting reduced
    !     precision number. If not specified then the result will have the
    !     default precision.
    !
    ! Returns:
    !
    ! * z: rpe_var
    !     An `rpe_var` instance representing the input literal at the given
    !     precision.
    !
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: x
        INTEGER, OPTIONAL,        INTENT(IN) :: n
        TYPE(rpe_var) :: z
        IF (PRESENT(n)) THEN
            z%sbits = n
        END IF
        z = x
    END FUNCTION rpe_literal_real

    FUNCTION rpe_literal_alternate (x, n) RESULT (z)
    ! Create an `rpe_var` instance from a real literal.
    !
    ! Arguments:
    !
    ! * x: real(kind=RPE_ALTERNATE_KIND) [input]
    !     The literal to transform to a reduced precision `rpe_var` instance.
    !
    ! * n: integer [input, optional]
    !     The number of bits in the significand of the resulting reduced
    !     precision number. If not specified then the result will have the
    !     default precision.
    !
    ! Returns:
    !
    ! * z: rpe_var
    !     An `rpe_var` instance representing the input literal at the given
    !     precision.
    !
        REAL(KIND=RPE_ALTERNATE_KIND),           INTENT(IN) :: x
        INTEGER,                       OPTIONAL, INTENT(IN) :: n
        TYPE(rpe_var) :: z
        IF (PRESENT(n)) THEN
            z%sbits = n
        END IF
        z = x
    END FUNCTION rpe_literal_alternate

    FUNCTION rpe_literal_integer (x, n) RESULT (z)
    ! Create an `rpe_var` instance from an integer literal.
    !
    ! Arguments:
    !
    ! * x: integer [input]
    !     The literal to transform to a reduced precision `rpe_var` instance.
    !
    ! * n: integer [input, optional]
    !     The number of bits in the significand of the resulting reduced
    !     precision number. If not specified then the result will have the
    !     default precision.
    !
    ! Returns:
    !
    ! * z: rpe_var
    !     An `rpe_var` instance representing the input literal at the given
    !     precision.
    !
        INTEGER,           INTENT(IN) :: x
        INTEGER, OPTIONAL, INTENT(IN) :: n
        TYPE(rpe_var) :: z
        IF (PRESENT(n)) THEN
            z%sbits = n
        END IF
        z = x
    END FUNCTION rpe_literal_integer

    FUNCTION rpe_literal_long (x, n) RESULT (z)
    ! Create an `rpe_var` instance from a long integer literal.
    !
    ! Arguments:
    !
    ! * x: integer(KIND=8) [input]
    !     The literal to transform to a reduced precision `rpe_var` instance.
    !
    ! * n: integer [input, optional]
    !     The number of bits in the significand of the resulting reduced
    !     precision number. If not specified then the result will have the
    !     default precision.
    !
    ! Returns:
    !
    ! * z: rpe_var
    !     An `rpe_var` instance representing the input literal at the given
    !     precision.
    !
        INTEGER(KIND=8),           INTENT(IN) :: x
        INTEGER,         OPTIONAL, INTENT(IN) :: n
        TYPE(rpe_var) :: z
        IF (PRESENT(n)) THEN
            z%sbits = n
        END IF
        z = x
    END FUNCTION rpe_literal_long


!-----------------------------------------------------------------------
! Overloaded assignment definitions:
!-----------------------------------------------------------------------

    ELEMENTAL SUBROUTINE assign_rpe_rpe (r1, r2)
    ! Assign an `rpe_type` instance to another `rpe_type` instance.
    !
    ! Arguments:
    !
    ! * r1: class(rpe_type) [input/output]
    !       An `rpe_type` instance to assign to.
    !
    ! * r2: class(rpe_type) [input]
    !       An `rpe_type` instance whose value will be assigned to `r1`.
    !
        TYPE(rpe_var), INTENT(INOUT) :: r1
        TYPE(rpe_var), INTENT(IN)    :: r2
        r1%val = r2%val
        CALL apply_truncation (r1)
    END SUBROUTINE assign_rpe_rpe

    ELEMENTAL SUBROUTINE assign_rpe_real (rpe, x)
    ! Assign a real variable to an `rpe_type` instance.
    !
    ! Arguments:
    !
    ! * rpe: class(rpe_type) [input/output]
    !       An `rpe_type` instance to assign to.
    !
    ! * x: real(kind=RPE_REAL_KIND) [input]
    !       A real variable whose value will be assigned to `rpe`.
    !
        TYPE(rpe_var),            INTENT(INOUT) :: rpe
        REAL(KIND=RPE_REAL_KIND), INTENT(IN)    :: x
        rpe%val = x
        CALL apply_truncation (rpe)
    END SUBROUTINE assign_rpe_real

    ELEMENTAL SUBROUTINE assign_rpe_alternate (rpe, x)
    ! Assign a real variable to an `rpe_type` instance.
    !
    ! Arguments:
    !
    ! * rpe: class(rpe_type) [input/output]
    !       An `rpe_type` instance to assign to.
    !
    ! * x: real(kind=RPE_ALTERNATE_KIND) [input]
    !       A real variable whose value will be assigned to `rpe`.
    !
        TYPE(rpe_var),                 INTENT(INOUT) :: rpe
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN)    :: x
        rpe%val = x
        CALL apply_truncation (rpe)
    END SUBROUTINE assign_rpe_alternate

    ELEMENTAL SUBROUTINE assign_rpe_integer (rpe, x)
    ! Assign an integer variable to an `rpe_type` instance.
    !
    ! Arguments:
    !
    ! * rpe: class(rpe_type) [input/output]
    !       An `rpe_type` instance to assign to.
    !
    ! * x: integer(kind=4) [input]
    !       An integer variable whose value will be assigned to `rpe`.
    !
        TYPE(rpe_var),   INTENT(INOUT) :: rpe
        INTEGER(KIND=4), INTENT(IN)    :: x
        rpe%val = x
        CALL apply_truncation (rpe)
    END SUBROUTINE assign_rpe_integer

    ELEMENTAL SUBROUTINE assign_rpe_long (rpe, x)
    ! Assign a long integer variable to an `rpe_type` instance.
    !
    ! Arguments:
    !
    ! * rpe: class(rpe_type) [input/output]
    !       An `rpe_type` instance to assign to.
    !
    ! * x: integer(kind=8) [input]
    !       A long integer variable whose value will be assigned to `rpe`.
    !
        TYPE(rpe_var),   INTENT(INOUT) :: rpe
        INTEGER(KIND=8), INTENT(IN)    :: x
        rpe%val = x
        CALL apply_truncation (rpe)
    END SUBROUTINE assign_rpe_long

    ELEMENTAL SUBROUTINE assign_real_rpe (x, rpe)
    ! Assign an `rpe_type` instance to a real variable.
    !
    ! Arguments:
    !
    ! * x: real(kind=RPE_REAL_KIND) [input/output]
    !       A real variable assign to.
    !
    ! * rpe: class(rpe_type) [input]
    !       An `rpe_type` instance whose value will be assigned to `x`.
    !
        REAL(KIND=RPE_REAL_KIND), INTENT(INOUT) :: x
        TYPE(rpe_var),            INTENT(IN)    :: rpe
        x = rpe%val
    END SUBROUTINE assign_real_rpe

    ELEMENTAL SUBROUTINE assign_alternate_rpe (x, rpe)
    ! Assign an `rpe_type` instance to a real variable.
    !
    ! Arguments:
    !
    ! * x: real(kind=RPE_ALTERNATE_KIND) [input/output]
    !       A real variable assign to.
    !
    ! * rpe: class(rpe_type) [input]
    !       An `rpe_type` instance whose value will be assigned to `x`.
    !
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(INOUT) :: x
        TYPE(rpe_var),                 INTENT(IN)    :: rpe
        x = rpe%val
    END SUBROUTINE assign_alternate_rpe

    ELEMENTAL SUBROUTINE assign_integer_rpe (x, rpe)
    ! Assign an `rpe_type` instance to an integer variable.
    !
    ! Arguments:
    !
    ! * x: integer(kind=4) [input/output]
    !       An integer variable assign to.
    !
    ! * rpe: class(rpe_type) [input]
    !       An `rpe_type` instance whose value will be assigned to `x`.
    !
        INTEGER(KIND=4), INTENT(INOUT) :: x
        TYPE(rpe_var),   INTENT(IN)    :: rpe
        x = rpe%val
    END SUBROUTINE assign_integer_rpe

    ELEMENTAL SUBROUTINE assign_long_rpe (x, rpe)
    ! Assign an `rpe_type` instance to a long integer variable.
    !
    ! Arguments:
    !
    ! * x: integer(kind=8) [input/output]
    !       A long integer variable assign to.
    !
    ! * rpe: class(rpe_type) [input]
    !       An `rpe_type` instance whose value will be assigned to `x`.
    !
        INTEGER(KIND=8), INTENT(INOUT) :: x
        TYPE(rpe_var),   INTENT(IN)    :: rpe
        x = rpe%val
    END SUBROUTINE assign_long_rpe

!-----------------------------------------------------------------------
! Overloaded operator definitions (external):
!-----------------------------------------------------------------------

#include "./include/implementation_operators.f90"

!-----------------------------------------------------------------------
! Overloaded intrinsic function definitions (external):
!-----------------------------------------------------------------------

#include "./include/implementation_intrinsics.f90"

!-----------------------------------------------------------------------
! Other extensions (external):
!-----------------------------------------------------------------------

#include "./include/implementation_extras.f90"

END MODULE rp_emulator
