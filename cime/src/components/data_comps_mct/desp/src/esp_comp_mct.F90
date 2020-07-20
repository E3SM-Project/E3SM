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
    use seq_comm_mct,     only: seq_comm_inst, seq_comm_name, seq_comm_suffix
    use seq_infodata_mod, only: seq_infodata_putData, seq_infodata_GetData
    use seq_timemgr_mod,  only: seq_timemgr_pause_component_index


    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock),           intent(inout) :: EClock
    type(seq_cdata),            intent(inout) :: cdata
    type(mct_aVect),            intent(inout) :: x2a, a2x   ! Not used
    character(len=*), optional, intent(in)    :: NLFilename ! Not used

                                                                   !EOP

    integer                            :: inst_index  ! (e.g., 1)
    character(len=16)                  :: inst_name   ! (e.g. "exp_0001")
    character(len=16)                  :: inst_suffix ! (e.g. "_0001" or "")
    integer                            :: ind
    integer                            :: ESPID
    integer                            :: mpicom_esp
    integer                            :: esp_phase
    logical                            :: esp_present
    logical                            :: esp_prognostic
    logical                            :: read_restart
    type(seq_infodata_type), pointer   :: infodata
    character(len=*),        parameter :: subName = "(esp_init_mct) "
    !--------------------------------------------------------------------------

    ! Retrieve  info for init method
    call seq_cdata_setptrs(cdata, ID=ESPID, mpicom=mpicom_esp,                &
         infodata=infodata)
    call seq_infodata_getData(infodata, esp_phase=esp_phase,                  &
         read_restart=read_restart)

    inst_name   = seq_comm_name(ESPID)
    inst_index  = seq_comm_inst(ESPID)
    inst_suffix = seq_comm_suffix(ESPID)

    call desp_comp_init(EClock, ESPID, mpicom_esp, esp_phase, read_restart,   &
         inst_name, inst_index, inst_suffix, esp_present, esp_prognostic)

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
    use shr_kind_mod,        only: CL=>SHR_KIND_CL
    use seq_comm_mct,        only: num_inst_atm, num_inst_lnd, num_inst_rof
    use seq_comm_mct,        only: num_inst_ocn, num_inst_ice, num_inst_glc
    use seq_comm_mct,        only: num_inst_wav, num_inst_max, num_inst_driver
    use desp_comp_mod,       only: desp_comp_run, desp_bcast_res_files
    use seq_cdata_mod,       only: seq_cdata_setptrs
    use seq_infodata_mod,    only: seq_infodata_putData, seq_infodata_GetData
    use seq_pauseresume_mod, only: seq_resume_get_files, seq_resume_store_comp
    use seq_timemgr_mod,     only: seq_timemgr_pause_component_active

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(inout)    :: EClock
    type(seq_cdata),  intent(inout)    :: cdata
    type(mct_aVect),  intent(inout)    :: x2a ! Not used
    type(mct_aVect),  intent(inout)    :: a2x ! Not used

    !EOP

    integer                            :: ind
    type(seq_infodata_type), pointer   :: infodata
    logical                            :: pause_sig(desp_num_comps)
    character(len=CL),       pointer   :: atm_resume(:)
    character(len=CL),       pointer   :: lnd_resume(:)
    character(len=CL),       pointer   :: rof_resume(:)
    character(len=CL),       pointer   :: ocn_resume(:)
    character(len=CL),       pointer   :: ice_resume(:)
    character(len=CL),       pointer   :: glc_resume(:)
    character(len=CL),       pointer   :: wav_resume(:)
    character(len=CL),       pointer   :: cpl_resume(:)
    character(len=CL)                  :: case_name
    character(len=*),        parameter :: subName = "(esp_run_mct) "
    !--------------------------------------------------------------------------

    ! Grab infodata and case name
    call seq_cdata_setptrs(cdata, infodata=infodata)
    call seq_infodata_GetData(infodata, case_name=case_name)
    ! Grab any active resume filenames
    call seq_resume_get_files('a', atm_resume, bcast=desp_bcast_res_files('a'))
    call seq_resume_get_files('l', lnd_resume, bcast=desp_bcast_res_files('l'))
    call seq_resume_get_files('o', ocn_resume, bcast=desp_bcast_res_files('o'))
    call seq_resume_get_files('i', ice_resume, bcast=desp_bcast_res_files('i'))
    call seq_resume_get_files('r', rof_resume, bcast=desp_bcast_res_files('r'))
    call seq_resume_get_files('g', glc_resume, bcast=desp_bcast_res_files('g'))
    call seq_resume_get_files('w', wav_resume, bcast=desp_bcast_res_files('w'))
    call seq_resume_get_files('x', cpl_resume, bcast=desp_bcast_res_files('x'))
    ! Find out if we should be running
    do ind = 1, desp_num_comps
       pause_sig(ind) = seq_timemgr_pause_component_active(comp_index(ind))
    end do
    call desp_comp_run(EClock, case_name, pause_sig, atm_resume, lnd_resume,  &
         rof_resume, ocn_resume, ice_resume, glc_resume, wav_resume,          &
         cpl_resume)

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
