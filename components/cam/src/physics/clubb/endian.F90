!----------------------------------------------------------------------
! $Id: endian.F90 7920 2015-11-07 04:01:20Z raut@uwm.edu $
!----------------------------------------------------------------------
module endian

! Description:
!   big_endian and little_endian are parameters set at compile time
!   based on whether the architecture is big or little endian.

!   native_4byte_real is a portable byte re-ordering subroutine
!   native_8byte_real is a knock off of the other routine for 8 bytes
! References:
!   big_endian, little_endian from:
!   <http://www.star.le.ac.uk/~cgp/f90course/f90.html>
!----------------------------------------------------------------------

  implicit none

  interface byte_order_swap
    module procedure native_4byte_real, native_8byte_real 
  end interface
  
  public  :: big_endian, little_endian, byte_order_swap
  private :: native_4byte_real, native_8byte_real

  private ! Default scope
  ! External
  intrinsic :: selected_int_kind, ichar, transfer

  ! Parameters
  integer, parameter :: &
    i4  = 4, & ! 4 byte long integer
    ich = ichar( transfer( 1_i4, "a" ) )

  logical, parameter :: &
    big_endian    = ich == 0, &
    little_endian = .not. big_endian

  contains

!-------------------------------------------------------------------------------
!     SUBPROGRAM: native_4byte_real
!
!         AUTHOR: David Stepaniak, NCAR/CGD/CAS
! DATE INITIATED: 29 April 2003 
!  LAST MODIFIED: 19 April 2005
!
!       SYNOPSIS: Converts a 32 bit, 4 byte, REAL from big Endian to
!                 little Endian, or conversely from little Endian to big
!                 Endian.
!
!    DESCRIPTION: This subprogram allows one to convert a 32 bit, 4 byte,
!                 REAL data element that was generated with, say, a big
!                 Endian processor (e.g. Sun/sparc, SGI/R10000, etc.) to its
!                 equivalent little Endian representation for use on little
!                 Endian processors (e.g. PC/Pentium running Linux). The
!                 converse, little Endian to big Endian, also holds.
!                 This conversion is accomplished by writing the 32 bits of
!                 the REAL data element into a generic 32 bit INTEGER space
!                 with the TRANSFER intrinsic, reordering the 4 bytes with
!                 the MVBITS intrinsic, and writing the reordered bytes into
!                 a new 32 bit REAL data element, again with the TRANSFER
!                 intrinsic. The following schematic illustrates the
!                 reordering process
!
!
!                  --------    --------    --------    --------
!                 |    D   |  |    C   |  |    B   |  |    A   |  4 Bytes
!                  --------    --------    --------    --------
!                                                             |
!                                                              -> 1 bit
!                                       ||
!                                     MVBITS
!                                       ||
!                                       \/
!
!                  --------    --------    --------    --------
!                 |    A   |  |    B   |  |    C   |  |    D   |  4 Bytes
!                  --------    --------    --------    --------
!                         |           |           |           |
!                         24          16          8           0   <- bit
!                                                                 position
!
!          INPUT: realIn,  a single 32 bit, 4 byte REAL data element.
!         OUTPUT: realOut, a single 32 bit, 4 byte REAL data element, with
!                 reverse byte order to that of realIn.
!    RESTRICTION: It is assumed that the default REAL data element is
!                 32 bits / 4 bytes.
!
!-----------------------------------------------------------------------
  SUBROUTINE native_4byte_real( realInOut )

    IMPLICIT NONE

    ! Added by Eric Raut, Nov 2015
    integer, parameter :: &
      sp = selected_real_kind( 6 ), & ! 32-bit floating point kind
      int32 = selected_int_kind( 9 )  ! 32-bit integer kind

    REAL(KIND=sp), INTENT(INOUT):: realInOut      ! a single 32 bit, 4 byte
                                                   ! REAL data element
!   Modified 8/1/05 
!   I found transfer does not work on pgf90 when -r8 is used and the mold
!   is a literal constant real; Changed the mold "0.0" to "readInOut"
!   -dschanen
!
!   REAL, INTENT(IN):: realInOut
!   REAL, INTENT(OUT) :: realOut
!                                                   ! a single 32 bit, 4 byte
!                                                   ! REAL data element, with
!                                                   ! reverse byte order to
!                                                   ! that of realIn
!----------------------------------------------------------------------
! Local variables (generic 32 bit INTEGER spaces):

    INTEGER(KIND=int32)                               :: i_element
    INTEGER(KIND=int32)                               :: i_element_br
!----------------------------------------------------------------------
! Transfer 32 bits of realIn to generic 32 bit INTEGER space:
    i_element = TRANSFER( realInOut, i_element )
!----------------------------------------------------------------------
! Reverse order of 4 bytes in 32 bit INTEGER space:
    CALL MVBITS( i_element, 24, 8, i_element_br, 0  )
    CALL MVBITS( i_element, 16, 8, i_element_br, 8  )
    CALL MVBITS( i_element,  8, 8, i_element_br, 16 )
    CALL MVBITS( i_element,  0, 8, i_element_br, 24 )
!----------------------------------------------------------------------
! Transfer reversed order bytes to 32 bit REAL space (realOut):
    realInOut = TRANSFER( i_element_br, realInOut )

    RETURN
  END SUBROUTINE native_4byte_real

!-------------------------------------------------------------------------------
  subroutine native_8byte_real( realInOut )

! Description:
!   This is just a modification of the above routine for 64 bit data
!-------------------------------------------------------------------------------

    ! Added by Eric Raut, Nov 2015
    use clubb_precision, only: &
      dp      ! Constant (64-bit floating point kind)

    implicit none

    ! Added by Eric Raut, Nov 2015
    integer, parameter :: &
      int64 = selected_int_kind( 18 )  ! 64-bit integer kind

    ! External
    intrinsic :: mvbits, transfer

    real(kind=dp), intent(inout) :: realInOut   ! a single 64 bit, 8 byte
                                                   ! REAL data element
    ! Local variables (generic 64 bit INTEGER spaces):

    integer(kind=int64) :: i_element
    integer(kind=int64) :: i_element_br

!-------------------------------------------------------------------------------

    ! Transfer 64 bits of realIn to generic 64 bit INTEGER space:
    i_element = transfer( realInOut, i_element )

    ! Reverse order of 8 bytes in 64 bit INTEGER space:
    call mvbits( i_element, 56, 8, i_element_br, 0  )
    call mvbits( i_element, 48, 8, i_element_br, 8  )
    call mvbits( i_element, 40, 8, i_element_br, 16 )
    call mvbits( i_element, 32, 8, i_element_br, 24 )
    call mvbits( i_element, 24, 8, i_element_br, 32 )
    call mvbits( i_element, 16, 8, i_element_br, 40 )
    call mvbits( i_element,  8, 8, i_element_br, 48 )
    call mvbits( i_element,  0, 8, i_element_br, 56 )

    ! Transfer reversed order bytes to 64 bit REAL space (realOut):
    realInOut = transfer( i_element_br, realInOut )

    return
  end subroutine native_8byte_real
!-------------------------------------------------------------------------------

end module endian

!-------------------------------------------------------------------------------
