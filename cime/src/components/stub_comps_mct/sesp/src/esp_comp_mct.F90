module esp_comp_mct

  ! !USES:

  use mct_mod,          only: mct_aVect
  use esmf,             only: ESMF_Clock
  use seq_cdata_mod,    only: seq_cdata
  use seq_infodata_mod, only: seq_infodata_PutData

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: esp_init_mct
  public :: esp_run_mct
  public :: esp_final_mct
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: esp_init_mct
  !
  ! !DESCRIPTION:
  !     stub esp model init
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine esp_init_mct( EClock, cdata, x2d, d2x, NLFilename )

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2d, d2x
    character(len=*), optional  , intent(in)    :: NLFilename

    !EOP
    !-------------------------------------------------------------------------------

    call seq_infodata_PutData(cdata%infodata, esp_present=.false., &
         esp_prognostic=.false., esp_phase=1)

  end subroutine esp_init_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: esp_run_mct
  !
  ! !DESCRIPTION:
  !     stub esp model run
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine esp_run_mct( EClock, cdata, x2d, d2x)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d
    type(mct_aVect)             ,intent(inout) :: d2x

    !EOP
    !-------------------------------------------------------------------------------

  end subroutine esp_run_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: esp_final_mct
  !
  ! !DESCRIPTION:
  !     stub esp model finalize
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------
  !
  subroutine esp_final_mct( EClock, cdata, x2d, d2x)

    implicit none

    !----- arguments -----
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d
    type(mct_aVect)             ,intent(inout) :: d2x

    !EOP
    !-------------------------------------------------------------------------------

  end subroutine esp_final_mct

  !===============================================================================

end module esp_comp_mct
