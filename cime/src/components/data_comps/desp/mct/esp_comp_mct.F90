module esp_comp_mct

  ! !USES:

  use esmf,              only: ESMF_Clock
  use mct_mod,           only: mct_aVect
  use seq_cdata_mod,     only: seq_cdata
  use seq_infodata_mod,  only: seq_infodata_type
  use desp_comp_mod,     only: desp_num_comps

  ! !PUBLIC TYPES:
  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: esp_init_mct
  public :: esp_run_mct
  public :: esp_final_mct

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------
  integer :: comp_index(desp_num_comps)

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
    use desp_comp_mod,    only: desp_comp_init, comp_names
    use seq_cdata_mod,    only: seq_cdata_setptrs
    use seq_infodata_mod, only: seq_infodata_putData, seq_infodata_GetData
    use seq_timemgr_mod,  only: seq_timemgr_pause_component_index


    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock),                  intent(inout) :: EClock
    type(seq_cdata),                   intent(inout) :: cdata
    type(mct_aVect),                   intent(inout) :: x2a, a2x   ! Not used
    character(len=*),        optional, intent(in)    :: NLFilename ! Not used

    !EOP

    integer                                          :: ind
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

    ! Retrieve component indices from the time manager
    do ind = 1, desp_num_comps
       comp_index(ind) = seq_timemgr_pause_component_index(comp_names(ind))
    end do

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
    use shr_kind_mod,     only: CL=>SHR_KIND_CL
    use seq_comm_mct,     only: num_inst_atm, num_inst_lnd, num_inst_rof
    use seq_comm_mct,     only: num_inst_ocn, num_inst_ice, num_inst_glc
    use seq_comm_mct,     only: num_inst_wav
    use desp_comp_mod,    only: desp_comp_run
    use seq_cdata_mod,    only: seq_cdata_setptrs
    use seq_infodata_mod, only: seq_infodata_putData, seq_infodata_GetData
    use seq_timemgr_mod,  only: seq_timemgr_pause_component_active

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(inout)    :: EClock
    type(seq_cdata),  intent(inout)    :: cdata
    type(mct_aVect),  intent(inout)    :: x2a ! Not used
    type(mct_aVect),  intent(inout)    :: a2x ! Not used

    !EOP

    integer                            :: ind
    type(seq_infodata_type), pointer   :: infodata
    logical                            :: pause_sig(desp_num_comps)
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

    ! Grab infodata and case name
    call seq_cdata_setptrs(cdata, infodata=infodata)
    call seq_infodata_GetData(infodata, case_name=case_name)
    ! Find out if we should be running
    do ind = 1, desp_num_comps
       pause_sig(ind) = seq_timemgr_pause_component_active(comp_index(ind))
    end do
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
