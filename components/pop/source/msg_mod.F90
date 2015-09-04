!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module msg_mod

!BOP
! !MODULE: msg_mod
!
! !DESCRIPTION:
!  This module prints messages. 
!  All msg_write calls have a routine name as the first argument
!  positional notation added to arguments to get module to compile
!-----------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  SVN:$Id: msg_mod.F90 808 2006-04-28 17:06:38Z njn01 $


! !USES:
  use kinds_mod

  implicit none
  private
  save


! !PUBLIC MEMBER FUNCTIONS:
  public :: &
       msg_set_state, &
       msg_get_state, &
       msg_set_iunit, &
       msg_get_iunit, &
       msg_write

!EOP
!BOC

  logical(kind=log_kind) :: msg_state = .TRUE.

  integer(kind=int_kind) :: msg_iunit = 6

  !-----------------------------------------------------------------------------
  !   generic interfaces
  !-----------------------------------------------------------------------------

  INTERFACE msg_write
     MODULE PROCEDURE &
          msg_write_A, &
          msg_write_AA, &
          msg_write_AI, &
          msg_write_AAA, &
          msg_write_AAI, &
          msg_write_AIA, &
          msg_write_AAAA, &
          msg_write_AIAI, &
          msg_write_AAAI, &
          msg_write_AAIA, &
          msg_write_AAAIAI
  END INTERFACE

!EOC
!*****************************************************************************

contains

!*****************************************************************************

  subroutine msg_set_state(state)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    logical(kind=log_kind), intent(in) :: state

    msg_state = state

  end subroutine msg_set_state

!*****************************************************************************

  subroutine msg_get_state(state)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    logical(kind=log_kind), intent(out) :: state

    state = msg_state

  end subroutine msg_get_state

!*****************************************************************************

  subroutine msg_set_iunit(msg_iunit_in)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    integer(kind=int_kind), intent(in) :: msg_iunit_in

    msg_iunit = msg_iunit_in

  end subroutine msg_set_iunit

!*****************************************************************************

  subroutine msg_get_iunit(msg_iunit_out)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    integer(kind=int_kind), intent(out) :: msg_iunit_out

    msg_iunit_out = msg_iunit

  end subroutine msg_get_iunit

!*****************************************************************************

  subroutine msg_write_A(sub_name, A1_1)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    character(len=*), intent(in) :: sub_name, A1_1

    if (msg_state) then
       write(unit=msg_iunit, fmt="('(', A, ') ', A)") &
            sub_name, A1_1
    end if

  end subroutine msg_write_A

!*****************************************************************************

  subroutine msg_write_AA(sub_name, A1_1, A2_2)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    character(len=*), intent(in) :: sub_name, A1_1, A2_2

    if (msg_state) then
       write(unit=msg_iunit, fmt="('(', A, ') ', A, A)") &
            sub_name, A1_1, A2_2
    end if

  end subroutine msg_write_AA

!*****************************************************************************

  subroutine msg_write_AI(sub_name, A1_1, I1_2)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    character(len=*), intent(in) :: sub_name, A1_1
    integer(kind=int_kind), intent(in) :: I1_2

    if (msg_state) then
       write(unit=msg_iunit, fmt="('(', A, ') ', A, I6)") &
            sub_name, A1_1, I1_2
    end if

  end subroutine msg_write_AI

!*****************************************************************************

  subroutine msg_write_AAA(sub_name, A1_1, A2_2, A3_3)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    character(len=*), intent(in) :: sub_name, A1_1, A2_2, A3_3

    if (msg_state) then
       write(unit=msg_iunit, fmt="('(', A, ') ', A, A, A)") &
            sub_name, A1_1, A2_2, A3_3
    end if

  end subroutine msg_write_AAA

!*****************************************************************************

  subroutine msg_write_AAI(sub_name, A1_1, A2_2, I1_3)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    character(len=*), intent(in) :: sub_name, A1_1, A2_2
    integer(kind=int_kind), intent(in) :: I1_3

    if (msg_state) then
       write(unit=msg_iunit, fmt="('(', A, ') ', A, A, I6)") &
            sub_name, A1_1, A2_2, I1_3
    end if

  end subroutine msg_write_AAI

!*****************************************************************************

  subroutine msg_write_AIA(sub_name, A1_1, I1_2, A2_3)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    character(len=*), intent(in) :: sub_name, A1_1, A2_3
    integer(kind=int_kind), intent(in) :: I1_2

    if (msg_state) then
       write(unit=msg_iunit, fmt="('(', A, ') ', A, I6, A)") &
            sub_name, A1_1, I1_2, A2_3
    end if

  end subroutine msg_write_AIA

!*****************************************************************************

  subroutine msg_write_AAAA(sub_name, A1_1, A2_2, A3_3, A4_4)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    character(len=*), intent(in) :: sub_name, A1_1, A2_2, A3_3, A4_4

    if (msg_state) then
       write(unit=msg_iunit, fmt="('(', A, ') ', A, A, A, A)") &
            sub_name, A1_1, A2_2, A3_3, A4_4
    end if

  end subroutine msg_write_AAAA

!*****************************************************************************

  subroutine msg_write_AIAI(sub_name, A1_1, I1_2, A2_3, I2_4)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    character(len=*), intent(in) :: sub_name, A1_1, A2_3
    integer(kind=int_kind), intent(in) :: I1_2, I2_4

    if (msg_state) then
       write(unit=msg_iunit, fmt="('(', A, ') ', A, I6, A, I6)") &
            sub_name, A1_1, I1_2, A2_3, I2_4
    end if

  end subroutine msg_write_AIAI

!*****************************************************************************

  subroutine msg_write_AAAI(sub_name, A1_1, A2_2, A3_3, I1_4)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    character(len=*), intent(in) :: sub_name, A1_1, A2_2, A3_3
    integer(kind=int_kind), intent(in) :: I1_4

    if (msg_state) then
       write(unit=msg_iunit, fmt="('(', A, ') ', A, A, A, I6)") &
            sub_name, A1_1, A2_2, A3_3, I1_4
    end if

  end subroutine msg_write_AAAI

!*****************************************************************************

  subroutine msg_write_AAIA(sub_name, A1_1, A2_2, I1_3, A3_4)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    character(len=*), intent(in) :: sub_name, A1_1, A2_2, A3_4
    integer(kind=int_kind), intent(in) :: I1_3

    if (msg_state) then
       write(unit=msg_iunit, fmt="('(', A, ') ', A, A, I6, A)") &
            sub_name, A1_1, A2_2, I1_3, A3_4
    end if

  end subroutine msg_write_AAIA

!*****************************************************************************

  subroutine msg_write_AAAIAI(sub_name, A1_1, A2_2, A3_3, I1_4, A4_5, I2_6)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    character(len=*), intent(in) :: sub_name, A1_1, A2_2, A3_3, A4_5
    integer(kind=int_kind), intent(in) :: I1_4, I2_6

    if (msg_state) then
       write(unit=msg_iunit, &
            fmt="('(', A, ') ', A, A, A, I6, A, I6)") &
            sub_name, A1_1, A2_2, A3_3, I1_4, A4_5, I2_6
    end if

  end subroutine msg_write_AAAIAI

!*****************************************************************************

end module msg_mod
