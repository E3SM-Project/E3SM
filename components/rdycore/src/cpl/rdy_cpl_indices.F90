module rdy_cpl_indices

  use shr_sys_mod, only : shr_sys_abort
  implicit none

  save
  private

  public :: rdy_cpl_indices_set

  integer, public :: index_x2rdy_Flrl_qsur = 0  ! lnd->rdy liquid surface runoff forcing from land
  integer, public :: index_x2rdy_Flrl_qsub = 0  ! lnd->rdy liquid subsurface runoff from land

contains

  subroutine rdy_cpl_indices_set()
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Set the coupler indices needed by the rdycore model coupler interface.
    !
    ! !USES:
    use seq_flds_mod , only : seq_flds_x2r_fields
    use mct_mod      , only : mct_aVect, mct_aVect_init, mct_avect_indexra, &
                              mct_aVect_clean
    !
    implicit none
    !
    ! !LOCAL VARIABLES:
    type(mct_aVect)   :: avtmp                            ! temporary
    character(len=32) :: subname = 'rdy_cpl_indicies_set' ! subroutine name

    call mct_aVect_init(avtmp, rList=seq_flds_x2r_fields, lsize=1)

    index_x2rdy_Flrl_qsur = mct_avect_indexra(avtmp,'Flrl_rofsur')
    index_x2rdy_Flrl_qsub = mct_avect_indexra(avtmp,'Flrl_rofsub')

    call mct_aVect_clean(avtmp)

  end subroutine rdy_cpl_indices_set

end module rdy_cpl_indices
