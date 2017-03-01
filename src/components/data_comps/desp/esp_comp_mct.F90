module esp_comp_mct

  ! !USES:

  use esmf,             only: ESMF_Clock
  use mct_mod,          only: mct_aVect
  use seq_cdata_mod,    only: seq_cdata
  use seq_infodata_mod, only: seq_infodata_type

  ! !PUBLIC TYPES:
  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: esp_init_mct
  public :: esp_run_mct
  public :: esp_final_mct

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !============================================================================
  !BOP ========================================================================
  !
  ! !IROUTINE: esp_init_mct
  !
  ! !DESCRIPTION:
  !     initialize data esp model
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ---------------------------------------------------------------

  subroutine esp_init_mct(EClock, cdata, x2a, a2x, NLFilename)
    use desp_comp_mod,    only: desp_comp_init
    use seq_cdata_mod,    only: seq_cdata_setptrs
    use seq_infodata_mod, only: seq_infodata_putData, seq_infodata_GetData

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock),                  intent(inout) :: EClock
    type(seq_cdata),                   intent(inout) :: cdata
    type(mct_aVect),                   intent(inout) :: x2a, a2x   ! Not used
    character(len=*),        optional, intent(in)    :: NLFilename ! Not used

    !EOP

    integer                                          :: ESPID
    integer                                          :: mpicom_esp
    integer                                          :: esp_phase
    logical                                          :: esp_present
    logical                                          :: esp_prognostic
    logical                                          :: read_restart
    type(seq_infodata_type), pointer                 :: infodata
    character(len=*),        parameter               :: subName = "(esp_init_mct) "
    !--------------------------------------------------------------------------


    ! Retrieve  info for init method
    call seq_cdata_setptrs(cdata, ID=ESPID, mpicom=mpicom_esp,                &
         infodata=infodata)
    call seq_infodata_getData(infodata, esp_phase=esp_phase,                  &
         read_restart=read_restart)

    call desp_comp_init(EClock, ESPID, mpicom_esp, esp_phase, read_restart,   &
         esp_present, esp_prognostic)

    ! Set the ESP model state
    call seq_infodata_PutData(infodata,                                       &
         esp_present=esp_present, esp_prognostic=esp_prognostic)

  end subroutine esp_init_mct

  !============================================================================
  !BOP ========================================================================
  !
  ! !IROUTINE: esp_run_mct
  !
  ! !DESCRIPTION:
  !     run method for data esp model
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine esp_run_mct( EClock, cdata,  x2a, a2x)
    use seq_comm_mct,     only: num_inst_atm, num_inst_lnd, num_inst_rof
    use seq_comm_mct,     only: num_inst_ocn, num_inst_ice, num_inst_glc
    use seq_comm_mct,     only: num_inst_wav
    use desp_comp_mod,    only: desp_comp_run, atm_ind, lnd_ind, ocn_ind
    use desp_comp_mod,    only: ice_ind, glc_ind, rof_ind, wav_ind, cpl_ind, max_ind
    use seq_cdata_mod,    only: seq_cdata_setptrs
    use seq_infodata_mod, only: seq_infodata_putData, seq_infodata_GetData
    use shr_kind_mod,     only: CL=>SHR_KIND_CL

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(inout)    :: EClock
    type(seq_cdata),  intent(inout)    :: cdata
    type(mct_aVect),  intent(inout)    :: x2a ! Not used
    type(mct_aVect),  intent(inout)    :: a2x ! Not used

    !EOP

    type(seq_infodata_type), pointer   :: infodata
    logical                            :: pause_sig(max_ind)
    character(len=CL)                  :: atm_resume(num_inst_atm)
    character(len=CL)                  :: lnd_resume(num_inst_lnd)
    character(len=CL)                  :: rof_resume(num_inst_rof)
    character(len=CL)                  :: ocn_resume(num_inst_ocn)
    character(len=CL)                  :: ice_resume(num_inst_ice)
    character(len=CL)                  :: glc_resume(num_inst_glc)
    character(len=CL)                  :: wav_resume(num_inst_wav)
    character(len=CL)                  :: cpl_resume
    character(len=CL)                  :: case_name       
    character(len=*),        parameter :: subName = "(esp_run_mct) "
    !--------------------------------------------------------------------------

    ! Grab infodata
    call seq_cdata_setptrs(cdata, infodata=infodata)
    ! Find out if we should be running
    ! The data ESP component only runs during a pause alarm
    call seq_infodata_GetData(infodata, atm_pause=pause_sig(atm_ind),         &
         lnd_pause=pause_sig(lnd_ind), ocn_pause=pause_sig(ocn_ind),          &
         ice_pause=pause_sig(ice_ind), glc_pause=pause_sig(glc_ind),          &
         rof_pause=pause_sig(rof_ind), wav_pause=pause_sig(wav_ind),          &
         cpl_pause=pause_sig(cpl_ind), case_name=case_name)

    call desp_comp_run(EClock, case_name, pause_sig, atm_resume, lnd_resume,  &
         rof_resume, ocn_resume, ice_resume, glc_resume, wav_resume,          &
         cpl_resume)
    ! Set any resume signals resulting from data ESP run
    call seq_infodata_PutData(infodata, atm_resume=atm_resume,                &
         lnd_resume=lnd_resume, ocn_resume=ocn_resume, ice_resume=ice_resume, &
         glc_resume=glc_resume, rof_resume=rof_resume, wav_resume=wav_resume, &
         cpl_resume=cpl_resume)

  end subroutine esp_run_mct

  !============================================================================
  !BOP ========================================================================
  !
  ! !IROUTINE: esp_final_mct
  !
  ! !DESCRIPTION:
  !     finalize method for data esp model
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ---------------------------------------------------------------
  !
  subroutine esp_final_mct(EClock, cdata, x2d, d2x)
    use desp_comp_mod, only: desp_comp_final

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> data
    type(mct_aVect)             ,intent(inout) :: d2x        ! data   -> driver

    !EOP

    !--- formats ---
    character(*), parameter :: subName = "(esp_final_mct) "
    !--------------------------------------------------------------------------

    call desp_comp_final()

  end subroutine esp_final_mct
  !============================================================================
  !============================================================================


end module esp_comp_mct
