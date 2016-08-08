module esmfshr_infodata_mod
#ifdef USE_ESMF_LIB

use ESMF
use seq_infodata_mod
use shr_kind_mod, only : SHR_KIND_CS, SHR_KIND_CL, SHR_KIND_IN,      &
                         SHR_KIND_R8, SHR_KIND_I8
use shr_sys_mod, only : shr_sys_abort

implicit none

public esmfshr_infodata_state2infodata
public esmfshr_infodata_infodata2state
private

contains

!--------------------------------------------------------------------
subroutine esmfshr_infodata_infodata2state(infodata, state, ID, rc)

    implicit none

    !inout parameters
    type(seq_infodata_type), intent(inout) :: infodata
    type(ESMF_State), intent(inout)     :: state
    integer, intent(in),  optional      :: ID
    integer, intent(out), optional      :: rc

    integer                    :: localrc, lID
    character(len=*),parameter :: subname = 'esmfshr_infodata_infodata2state'

    if(present(rc)) rc = ESMF_SUCCESS

    lID = -1
    if (present(ID)) then
       lID = ID
    endif

    call esmfshr_infodata_convert(infodata,lID,state,'i2s',localrc)

end subroutine esmfshr_infodata_infodata2state
!--------------------------------------------------------------------
subroutine esmfshr_infodata_state2infodata(state, infodata, ID, rc)

    implicit none

    !inout parameters
    type(seq_infodata_type), intent(inout) :: infodata
    type(ESMF_State), intent(inout)     :: state
    integer, intent(out), optional      :: ID
    integer, intent(out), optional      :: rc

    integer                    :: localrc, lID
    character(len=*),parameter :: subname = 'esmfshr_infodata_state2infodata'

    if(present(rc)) rc = ESMF_SUCCESS

    call esmfshr_infodata_convert(infodata,lID,state,'s2i',localrc)

    if (present(ID)) then
       ID = lID
    endif

end subroutine esmfshr_infodata_state2infodata
!--------------------------------------------------------------------

subroutine esmfshr_infodata_convert(infodata,ID,state,direction,rc)
    
    type(seq_infodata_type), intent(inout) :: infodata
    integer, intent(inout)                 :: ID
    type(ESMF_State), intent(inout)        :: state
    character(len=*), intent(in)           :: direction
    integer, intent(out), optional         :: rc

    ! local variables
    integer                 :: int_buf
    logical                 :: log_buf
    real(SHR_KIND_R8)       :: real_buf
    character(SHR_KIND_CL)  :: char_buf
    logical                 :: i2s

    integer                 :: localrc
    character(len=*),parameter :: subname = 'esmfshr_infodata_convert'

    if(present(rc)) rc = ESMF_SUCCESS

    if (trim(direction) == 'i2s') then
       i2s = .true.
    elseif (trim(direction) == 's2i') then
       i2s = .false.
    else
       write(6,*) subname,' ERROR: unknown direction ',trim(direction)
       call shr_sys_abort()
    endif

    if (i2s) then
       call ESMF_AttributeSet(state, name="ID", value=ID, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ID", value=ID, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, start_type=char_buf)
       call ESMF_AttributeSet(state, name="start_type", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="start_type", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, start_type=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, case_name=char_buf)
       call ESMF_AttributeSet(state, name="case_name", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="case_name", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, case_name=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, case_desc=char_buf)
       call ESMF_AttributeSet(state, name="case_desc", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="case_desc", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, case_desc=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, model_version=char_buf)
       call ESMF_AttributeSet(state, name="model_version", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="model_version", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, model_version=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, username=char_buf)
       call ESMF_AttributeSet(state, name="username", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="username", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, username=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, hostname=char_buf)
       call ESMF_AttributeSet(state, name="hostname", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="hostname", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, hostname=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, timing_dir=char_buf)
       call ESMF_AttributeSet(state, name="timing_dir", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="timing_dir", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, timing_dir=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, tchkpt_dir=char_buf)
       call ESMF_AttributeSet(state, name="tchkpt_dir", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="tchkpt_dir", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, tchkpt_dir=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, atm_adiabatic=log_buf)
       call ESMF_AttributeSet(state, name="atm_adiabatic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="atm_adiabatic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, atm_adiabatic=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, atm_ideal_phys=log_buf)
       call ESMF_AttributeSet(state, name="atm_ideal_phys", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="atm_ideal_phys", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, atm_ideal_phys=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, aqua_planet=log_buf)
       call ESMF_AttributeSet(state, name="aqua_planet", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="aqua_planet", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, aqua_planet=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, run_barriers=log_buf)
       call ESMF_AttributeSet(state, name="run_barriers", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="run_barriers", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, run_barriers=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, brnch_retain_casename=log_buf)
       call ESMF_AttributeSet(state, name="brnch_retain_casename", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="brnch_retain_casename", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, brnch_retain_casename=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, read_restart=log_buf)
       call ESMF_AttributeSet(state, name="read_restart", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="read_restart", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, read_restart=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, restart_pfile=char_buf)
       call ESMF_AttributeSet(state, name="restart_pfile", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="restart_pfile", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, restart_pfile=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, restart_file=char_buf)
       call ESMF_AttributeSet(state, name="restart_file", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="restart_file", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, restart_file=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, single_column=log_buf)
       call ESMF_AttributeSet(state, name="single_column", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="single_column", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, single_column=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, scmlat=real_buf)
       call ESMF_AttributeSet(state, name="scmlat", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="scmlat", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, scmlat=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, scmlon=real_buf)
       call ESMF_AttributeSet(state, name="scmlon", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="scmlon", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, scmlon=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, logFilePostFix=char_buf)
       call ESMF_AttributeSet(state, name="logFilePostFix", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="logFilePostFix", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, logFilePostFix=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, outPathRoot=char_buf)
       call ESMF_AttributeSet(state, name="outPathRoot", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="outPathRoot", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, outPathRoot=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, perpetual=log_buf)
       call ESMF_AttributeSet(state, name="perpetual", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="perpetual", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, perpetual=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, perpetual_ymd=int_buf)
       call ESMF_AttributeSet(state, name="perpetual_ymd", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="perpetual_ymd", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, perpetual_ymd=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, orb_mode=char_buf)
       call ESMF_AttributeSet(state, name="orb_mode", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="orb_mode", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, orb_mode=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, orb_iyear=int_buf)
       call ESMF_AttributeSet(state, name="orb_iyear", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="orb_iyear", value=int_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, orb_iyear=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, orb_iyear_align=int_buf)
       call ESMF_AttributeSet(state, name="orb_iyear_align", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="orb_iyear_align", value=int_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, orb_iyear_align=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, orb_eccen=real_buf)
       call ESMF_AttributeSet(state, name="orb_eccen", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="orb_eccen", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, orb_eccen=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, orb_obliqr=real_buf)
       call ESMF_AttributeSet(state, name="orb_obliqr", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="orb_obliqr", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, orb_obliqr=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, orb_lambm0=real_buf)
       call ESMF_AttributeSet(state, name="orb_lambm0", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="orb_lambm0", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, orb_lambm0=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, orb_mvelpp=real_buf)
       call ESMF_AttributeSet(state, name="orb_mvelpp", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="orb_mvelpp", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, orb_mvelpp=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, orb_obliq=real_buf)
       call ESMF_AttributeSet(state, name="orb_obliq", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="orb_obliq", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, orb_obliq=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, orb_mvelp=real_buf)
       call ESMF_AttributeSet(state, name="orb_mvelp", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="orb_mvelp", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, orb_mvelp=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, flux_epbal=char_buf)
       call ESMF_AttributeSet(state, name="flux_epbal", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="flux_epbal", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, flux_epbal=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, flux_albav=log_buf)
       call ESMF_AttributeSet(state, name="flux_albav", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="flux_albav", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, flux_albav=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, wall_time_limit=real_buf)
       call ESMF_AttributeSet(state, name="wall_time_limit", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="wall_time_limit", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, wall_time_limit=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, force_stop_at=char_buf)
       call ESMF_AttributeSet(state, name="force_stop_at", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="force_stop_at", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, force_stop_at=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, atm_gnam=char_buf)
       call ESMF_AttributeSet(state, name="atm_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="atm_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, atm_gnam=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, lnd_gnam=char_buf)
       call ESMF_AttributeSet(state, name="lnd_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="lnd_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, lnd_gnam=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, ocn_gnam=char_buf)
       call ESMF_AttributeSet(state, name="ocn_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ocn_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, ocn_gnam=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, ice_gnam=char_buf)
       call ESMF_AttributeSet(state, name="ice_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ice_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, ice_gnam=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, rof_gnam=char_buf)
       call ESMF_AttributeSet(state, name="rof_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="rof_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, rof_gnam=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, glc_gnam=char_buf)
       call ESMF_AttributeSet(state, name="glc_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="glc_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, glc_gnam=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, wav_gnam=char_buf)
       call ESMF_AttributeSet(state, name="wav_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="wav_gnam", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, wav_gnam=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, shr_map_dopole=log_buf)
       call ESMF_AttributeSet(state, name="shr_map_dopole", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="shr_map_dopole", value=log_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, shr_map_dopole=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, vect_map=char_buf)
       call ESMF_AttributeSet(state, name="vect_map", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="vect_map", value=char_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, vect_map=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, aoflux_grid=char_buf)
       call ESMF_AttributeSet(state, name="aoflux_grid", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="aoflux_grid", value=char_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, aoflux_grid=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, cpl_seq_option=char_buf)
       call ESMF_AttributeSet(state, name="cpl_seq_option", value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="cpl_seq_option", value=char_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, cpl_seq_option=char_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, cpl_cdf64=log_buf)
       call ESMF_AttributeSet(state, name="cpl_cdf64", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="cpl_cdf64", value=log_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, cpl_cdf64=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, do_budgets=log_buf)
       call ESMF_AttributeSet(state, name="do_budgets", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="do_budgets", value=log_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, do_budgets=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, budget_inst=int_buf)
       call ESMF_AttributeSet(state, name="budget_inst", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="budget_inst", value=int_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, budget_inst=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, budget_daily=int_buf)
       call ESMF_AttributeSet(state, name="budget_daily", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="budget_daily", value=int_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, budget_daily=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, budget_month=int_buf)
       call ESMF_AttributeSet(state, name="budget_month", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="budget_month", value=int_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, budget_month=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, budget_ann=int_buf)
       call ESMF_AttributeSet(state, name="budget_ann", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="budget_ann", value=int_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, budget_ann=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, budget_ltann=int_buf)
       call ESMF_AttributeSet(state, name="budget_ltann", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="budget_ltann", value=int_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, budget_ltann=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, budget_ltend=int_buf)
       call ESMF_AttributeSet(state, name="budget_ltend", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="budget_ltend", value=int_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, budget_ltend=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, drv_threading=log_buf)
       call ESMF_AttributeSet(state, name="drv_threading", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="drv_threading", value=log_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, drv_threading=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, eps_frac=real_buf)
       call ESMF_AttributeSet(state, name="eps_frac", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="eps_frac", value=real_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, eps_frac=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, eps_amask=real_buf)
       call ESMF_AttributeSet(state, name="eps_amask", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="eps_amask", value=real_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, eps_amask=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, eps_agrid=real_buf)
       call ESMF_AttributeSet(state, name="eps_agrid", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="eps_agrid", value=real_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, eps_agrid=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, eps_aarea=real_buf)
       call ESMF_AttributeSet(state, name="eps_aarea", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="eps_aarea", value=real_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, eps_aarea=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, eps_omask=real_buf)
       call ESMF_AttributeSet(state, name="eps_omask", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="eps_omask", value=real_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, eps_omask=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, eps_ogrid=real_buf)
       call ESMF_AttributeSet(state, name="eps_ogrid", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="eps_ogrid", value=real_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, eps_ogrid=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, eps_oarea=real_buf)
       call ESMF_AttributeSet(state, name="eps_oarea", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="eps_oarea", value=real_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, eps_oarea=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, info_debug=int_buf)
       call ESMF_AttributeSet(state, name="info_debug", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="info_debug", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, info_debug=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, bfbflag=log_buf)
       call ESMF_AttributeSet(state, name="bfbflag", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="bfbflag", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, bfbflag=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, atm_present=log_buf)
       call ESMF_AttributeSet(state, name="atm_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="atm_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, atm_present=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, atm_prognostic=log_buf)
       call ESMF_AttributeSet(state, name="atm_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="atm_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, atm_prognostic=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, lnd_present=log_buf)
       call ESMF_AttributeSet(state, name="lnd_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="lnd_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, lnd_present=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, lnd_prognostic=log_buf)
       call ESMF_AttributeSet(state, name="lnd_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="lnd_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, lnd_prognostic=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, rof_present=log_buf)
       call ESMF_AttributeSet(state, name="rof_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="rof_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, rof_present=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, rofice_present=log_buf)
       call ESMF_AttributeSet(state, name="rofice_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="rofice_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, rofice_present=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, rof_prognostic=log_buf)
       call ESMF_AttributeSet(state, name="rof_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="rof_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, rof_prognostic=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, flood_present=log_buf)
       call ESMF_AttributeSet(state, name="flood_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="flood_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, flood_present=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, ocn_present=log_buf)
       call ESMF_AttributeSet(state, name="ocn_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ocn_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, ocn_present=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, ocn_prognostic=log_buf)
       call ESMF_AttributeSet(state, name="ocn_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ocn_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, ocn_prognostic=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, ocnrof_prognostic=log_buf)
       call ESMF_AttributeSet(state, name="ocnrof_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ocnrof_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, ocnrof_prognostic=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, ice_present=log_buf)
       call ESMF_AttributeSet(state, name="ice_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ice_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, ice_present=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, ice_prognostic=log_buf)
       call ESMF_AttributeSet(state, name="ice_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ice_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, ice_prognostic=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, iceberg_prognostic=log_buf)
       call ESMF_AttributeSet(state, name="iceberg_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="iceberg_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, iceberg_prognostic=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, glc_present=log_buf)
       call ESMF_AttributeSet(state, name="glc_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="glc_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, glc_present=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, glclnd_present=log_buf)
       call ESMF_AttributeSet(state, name="glclnd_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="glclnd_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, glclnd_present=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, glcocn_present=log_buf)
       call ESMF_AttributeSet(state, name="glcocn_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="glcocn_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, glcocn_present=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, glcice_present=log_buf)
       call ESMF_AttributeSet(state, name="glcice_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="glcice_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, glcice_present=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, glc_prognostic=log_buf)
       call ESMF_AttributeSet(state, name="glc_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="glc_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, glc_prognostic=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, wav_present=log_buf)
       call ESMF_AttributeSet(state, name="wav_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="wav_present", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, wav_present=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, wav_prognostic=log_buf)
       call ESMF_AttributeSet(state, name="wav_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="wav_prognostic", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, wav_prognostic=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, dead_comps=log_buf)
       call ESMF_AttributeSet(state, name="dead_comps", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="dead_comps", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, dead_comps=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, atm_nx=int_buf)
       call ESMF_AttributeSet(state, name="atm_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="atm_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, atm_nx=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, atm_ny=int_buf)
       call ESMF_AttributeSet(state, name="atm_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="atm_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, atm_ny=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, lnd_nx=int_buf)
       call ESMF_AttributeSet(state, name="lnd_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="lnd_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, lnd_nx=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, lnd_ny=int_buf)
       call ESMF_AttributeSet(state, name="lnd_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="lnd_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, lnd_ny=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, ice_nx=int_buf)
       call ESMF_AttributeSet(state, name="ice_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ice_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, ice_nx=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, ice_ny=int_buf)
       call ESMF_AttributeSet(state, name="ice_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ice_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, ice_ny=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, ocn_nx=int_buf)
       call ESMF_AttributeSet(state, name="ocn_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ocn_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, ocn_nx=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, ocn_ny=int_buf)
       call ESMF_AttributeSet(state, name="ocn_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ocn_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, ocn_ny=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, rof_nx=int_buf)
       call ESMF_AttributeSet(state, name="rof_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="rof_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, rof_nx=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, rof_ny=int_buf)
       call ESMF_AttributeSet(state, name="rof_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="rof_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, rof_ny=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, glc_nx=int_buf)
       call ESMF_AttributeSet(state, name="glc_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="glc_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, glc_nx=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, glc_ny=int_buf)
       call ESMF_AttributeSet(state, name="glc_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="glc_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, glc_ny=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, wav_nx=int_buf)
       call ESMF_AttributeSet(state, name="wav_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="wav_nx", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, wav_nx=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, wav_ny=int_buf)
       call ESMF_AttributeSet(state, name="wav_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="wav_ny", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, wav_ny=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, nextsw_cday=real_buf)
       call ESMF_AttributeSet(state, name="nextsw_cday", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="nextsw_cday", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, nextsw_cday=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, precip_fact=real_buf)
       call ESMF_AttributeSet(state, name="precip_fact", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="precip_fact", value=real_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, precip_fact=real_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, atm_phase=int_buf)
       call ESMF_AttributeSet(state, name="atm_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="atm_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, atm_phase=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, lnd_phase=int_buf)
       call ESMF_AttributeSet(state, name="lnd_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="lnd_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, lnd_phase=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, rof_phase=int_buf)
       call ESMF_AttributeSet(state, name="rof_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="rof_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, rof_phase=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, ice_phase=int_buf)
       call ESMF_AttributeSet(state, name="ice_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ice_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, ice_phase=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, ocn_phase=int_buf)
       call ESMF_AttributeSet(state, name="ocn_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="ocn_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, ocn_phase=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, glc_phase=int_buf)
       call ESMF_AttributeSet(state, name="glc_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="glc_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, glc_phase=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, wav_phase=int_buf)
       call ESMF_AttributeSet(state, name="wav_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="wav_phase", value=int_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, wav_phase=int_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, glcrun_alarm=log_buf)
       call ESMF_AttributeSet(state, name="glcrun_alarm", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="glcrun_alarm", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, glcrun_alarm=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, atm_aero=log_buf)
       call ESMF_AttributeSet(state, name="atm_aero", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="atm_aero", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, atm_aero=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, glc_g2lupdate=log_buf)
       call ESMF_AttributeSet(state, name="glc_g2lupdate", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeGet(state, name="glc_g2lupdate", value=log_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
       call seq_infodata_Putdata(infodata, glc_g2lupdate=log_buf)
    endif

    if (i2s) then
       call seq_infodata_Getdata(infodata, rest_case_name=char_buf)
       call ESMF_AttributeSet(state, name="rest_case_name",  value=char_buf, rc=localrc)
       if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
    else
        call ESMF_AttributeGet(state, name="rest_case_name",  value=char_buf, rc=localrc)
        if(localrc /= ESMF_SUCCESS) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
        call seq_infodata_Putdata(infodata, rest_case_name=char_buf)
    endif

end subroutine esmfshr_infodata_convert

!--------------------------------------------------------------------

#endif
end module esmfshr_infodata_mod
