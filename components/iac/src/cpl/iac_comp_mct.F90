module iac_comp_mct
  
!========================================================================
! DESCRIPTION:
! Interface of integrated assessment component with the main E3SM driver. 
! This is a thin interface taking E3SM driver information
! in MCT (Model Coupling Toolkit) format and converting it to use by GCAM
! Shamelessly copied and stolen from other components (mostly
! RTM/ROF) in the ACME development tree, so if something looks wrong
! it's probably because I didn't understand correctly what the
! original code did.

! Note: Iac single letter descriptor is 'z', because all the other logical ones were
! already taken.

  use seq_flds_mod
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                shr_file_getLogUnit, shr_file_getLogLevel, &
                                shr_file_getUnit, shr_file_setIO
  use shr_const_mod    , only : SHR_CONST_REARTH
  use seq_cdata_mod    , only : seq_cdata, seq_cdata_setptrs
  use seq_timemgr_mod  , only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn, &
                                seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
  use seq_infodata_mod , only : seq_infodata_type, seq_infodata_GetData, seq_infodata_PutData, &
                                seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                                seq_infodata_start_type_brnch
  use seq_comm_mct     , only : seq_comm_suffix, seq_comm_inst, seq_comm_name

! Stub in other moduals for running iac
! use IacMod
! use IacVar

!
! PUBLIC MEMBER FUNCTIONS:
  implicit none
  SAVE
  private                              ! By default make data private
!
! PUBLIC MEMBER FUNCTIONS:

  public :: iac_init_mct               ! iac initialization
  public :: iac_run_mct                ! iac run phase
  public :: iac_final_mct              ! iac finalization/cleanup
!
! PRIVATE MEMBER FUNCTIONS:

! REVISION HISTORY:
! Author: Tim Shippert
!===============================================================
contains
!===============================================================

!========================================================================

  subroutine iac_init_mct( EClock, cdata_z, x2z_z, z2x_z, NLFilename)

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Initialize iac model and hook in arrays from lnd module
    !
    ! !ARGUMENTS:
    type(ESMF_Clock),           intent(inout) :: EClock     ! Input synchronization clock
    type(seq_cdata),            intent(inout) :: cdata_z    ! Input iac driver data
    type(mct_aVect) ,           intent(inout) :: x2z_z      ! Iac import state
    type(mct_aVect),            intent(inout) :: z2x_z      ! Iac export state
    character(len=*), optional, intent(in)    :: NLFilename ! Namelist filename to read
    !
    ! !LOCAL VARIABLES:
    ! 

  end subroutine iac_init_mct

!---------------------------------------------------------------------------

  subroutine iac_run_mct( EClock, cdata_z, x2z_z, z2x_z)

    !-------------------------------------------------------
    ! DESCRIPTION:
    ! Run IAC model

    ! ARGUMENTS:
    implicit none
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_z   ! Input driver data for runoff model
    type(mct_aVect)  , intent(inout) :: x2z_z     ! Import state from runoff model
    type(mct_aVect)  , intent(inout) :: z2x_z     ! Export state from runoff model

    ! LOCAL VARIABLES:

  end subroutine iac_run_mct

!===============================================================================

  subroutine iac_final_mct( EClock, cdata_z, x2z_z, z2x_z)

    !-----------------------------------------------------
    ! DESCRIPTION:
    ! Finalize iac surface model
    !
    ! ARGUMENTS:
    implicit none
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_z   ! Input driver data for runoff model
    type(mct_aVect)  , intent(inout) :: x2z_z     ! Import state from runoff model
    type(mct_aVect)  , intent(inout) :: z2x_z     ! Export state from runoff model
    !-----------------------------------------------------

   ! fill this in
  end subroutine iac_final_mct

!===============================================================================

!==============================================================
! Local functions
!==============================================================

end module iac_comp_mct
