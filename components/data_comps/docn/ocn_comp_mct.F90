module ocn_comp_mct

  ! !USES:

  use seq_cdata_mod
  use esmf
  use mct_mod

  use docn_comp_mod

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: ocn_init_mct
  !
  ! !DESCRIPTION:
  !     initialize data ocn model
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine ocn_init_mct( EClock, cdata, x2o, o2x, NLFilename )

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2o, o2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

    !EOP

    character(*), parameter :: subName = "(ocn_init_mct) "
    !-------------------------------------------------------------------------------


    if (present(NLFilename)) then
       call docn_comp_init(EClock, cdata, x2o, o2x, NLFilename)
    else
       call docn_comp_init(EClock, cdata, x2o, o2x)
    endif

  end subroutine ocn_init_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: ocn_run_mct
  !
  ! !DESCRIPTION:
  !     run method for dead ocn model
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine ocn_run_mct( EClock, cdata,  x2o, o2x)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2o        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: o2x        ! dead   -> driver

    !EOP

    character(*), parameter :: subName = "(ocn_run_mct) "
    !-------------------------------------------------------------------------------

    call docn_comp_run(EClock, cdata, x2o, o2x)

  end subroutine ocn_run_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: ocn_final_mct
  !
  ! !DESCRIPTION:
  !     finalize method for dead ocn model
  !
  ! !REVISION HISTORY:
  ! 
  ! !INTERFACE: ------------------------------------------------------------------
  !
  subroutine ocn_final_mct(EClock, cdata, x2d, d2x)

    implicit none

    !----- arguments -----

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

    !EOP

    !--- formats ---
    character(*), parameter :: subName = "(ocn_final_mct) "
    !-------------------------------------------------------------------------------

    call docn_comp_final()

  end subroutine ocn_final_mct
  !===============================================================================
  !===============================================================================


end module ocn_comp_mct
