module atm_comp_mct

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
    use ESM_ClockMod,  only: EClock
    use seq_cdata_mod, only: seq_cdata
    use mct_mod,       only: mct_aVect
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(inout)              :: EClock
    type(seq_cdata), intent(inout)              :: cdata_a
    type(mct_aVect), intent(inout)              :: x2a_a
    type(mct_aVect), intent(inout)              :: a2x_a   
    character(len=*), optional,   intent(IN)    :: NLFilename ! Namelist filename

    print *, "Hello world, from scream initialization function"

 end subroutine atm_init_mct

!================================================================================

 subroutine atm_run_mct( EClock, cdata_a, x2a_a, a2x_a)
    use ESM_ClockMod,  only: EClock
    use seq_cdata_mod, only: seq_cdata
    use mct_mod,       only: mct_aVect
    !
    ! Arguments
    !
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata_a
    type(mct_aVect)             ,intent(inout) :: x2a_a
    type(mct_aVect)             ,intent(inout) :: a2x_a

    print *, "Hello world, from scream run function"
    
  end subroutine atm_run_mct

!================================================================================

  subroutine atm_final_mct( EClock, cdata_a, x2a_a, a2x_a)
    use ESM_ClockMod,  only: EClock
    use seq_cdata_mod, only: seq_cdata
    use mct_mod,       only: mct_aVect
    !
    ! Arguments
    !
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata_a
    type(mct_aVect)             ,intent(inout) :: x2a_a
    type(mct_aVect)             ,intent(inout) :: a2x_a

    print *, "Hello world, from scream finalization function"

  end subroutine atm_final_mct

end module atm_comp_mct
