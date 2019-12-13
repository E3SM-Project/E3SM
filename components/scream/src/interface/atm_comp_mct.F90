module atm_comp_mct

  ! These are used across the module, so we issue a single use clause here
  use mct_mod,         only: mct_aVect
  use ESMF_ClockMod,   only: ESMF_Clock
  use seq_cdata_mod,   only: seq_cdata
  use seq_timemgr_mod, only: seq_timemgr_EClockGetData
  use shr_kind_mod,    only: SHR_KIND_IN

  implicit none

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: atm_init_mct
  public :: atm_run_mct
  public :: atm_final_mct

contains

!================================================================================

  subroutine atm_init_mct( EClock, cdata_a, x2a_a, a2x_a, NLFilename )
    !
    ! F90<->CXX interfaces
    !
    interface
      subroutine scream_init(f_comm,start_ymd,start_tod) bind(c)
        use iso_c_binding, only: c_int
        !
        ! Arguments
        !
        integer (kind=c_int), intent(in) :: start_tod, start_ymd, f_comm
      end subroutine scream_init
    end interface
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(inout)              :: EClock
    type(seq_cdata), intent(inout)              :: cdata_a
    type(mct_aVect), intent(inout)              :: x2a_a
    type(mct_aVect), intent(inout)              :: a2x_a
    character(len=*), optional,   intent(IN)    :: NLFilename ! Namelist filename
    !
    ! Locals
    !
    integer :: f_mpicomm
    integer (kind=SHR_KIND_IN) :: start_tod, start_ymd

    print *, "Hello world, from scream initialization function"

    ! Unpack data stored in cdata_a
    call seq_cdata_setptrs(cdata_a,mpicom=f_mpicomm)

    ! Get the time clock data
    ! TODO: what does scream need to init its time stamps?
    call seq_timemgr_EClockGetData(EClock, start_ymd=start_ymd, start_tod=start_tod)

    ! Init scream
    call scream_init (f_mpiocomm, INT(start_ycm, KIND=c_int), INT(start_tod, KIND=c_int))
 end subroutine atm_init_mct

!================================================================================

 subroutine atm_run_mct( EClock, cdata_a, x2a_a, a2x_a)
    !
    ! F90<->CXX interfaces
    !
    interface
      subroutine scream_run (dt) bind(c)
        use iso_c_binding, only: c_double
        !
        ! Arguments
        !
        integer(kind=c_double), intent(in) :: dt
      end subroutine scream_run
    end interface
    !
    ! Arguments
    !
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata_a
    type(mct_aVect)             ,intent(inout) :: x2a_a
    type(mct_aVect)             ,intent(inout) :: a2x_a

    print *, "Hello world, from scream run function"

    call scream_run(dt)

  end subroutine atm_run_mct

!================================================================================

  subroutine atm_final_mct( EClock, cdata_a, x2a_a, a2x_a)
    !
    ! F90<->CXX interfaces
    !
    interface
      subroutine scream_finalize () bind(c)
      end subroutine scream_finalize
    end interface
    !
    ! Arguments
    !
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata_a
    type(mct_aVect)             ,intent(inout) :: x2a_a
    type(mct_aVect)             ,intent(inout) :: a2x_a

    print *, "Hello world, from scream finalization function"

    call scream_finalize()

  end subroutine atm_final_mct

end module atm_comp_mct
