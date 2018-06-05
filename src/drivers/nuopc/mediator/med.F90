module MED

  !-----------------------------------------------------------------------------
  ! Mediator Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Mediator, only: &
    mediator_routine_SS             => SetServices, &
    mediator_routine_Run            => routine_Run, &
    mediator_label_DataInitialize   => label_DataInitialize, &
    mediator_label_Advance          => label_Advance, &
    mediator_label_CheckImport      => label_CheckImport, &
    mediator_label_TimestampExport  => label_TimestampExport, &
    mediator_label_SetRunClock      => label_SetRunClock, &
    NUOPC_MediatorGet

  use shr_kind_mod              , only: SHR_KIND_CX, SHR_KIND_CL, SHR_KIND_CS
  use shr_sys_mod               , only: shr_sys_flush, shr_sys_abort
  use esmFlds                   , only: flds_scalar_name
  use esmFlds                   , only: flds_scalar_num
  use esmFlds                   , only: fldListFr, fldListTo
  use esmFlds                   , only: ncomps, compmed, compatm, compocn
  use esmFlds                   , only: compice, complnd, comprof, compwav, compglc, compname
  use esmFlds                   , only: fldListMed_ocnalb_o, fldListMed_aoflux_a, fldListMed_aoflux_o
  use shr_nuopc_fldList_mod     , only: shr_nuopc_fldList_Realize
  use shr_nuopc_fldList_mod     , only: shr_nuopc_fldList_GetFldNames
  use shr_nuopc_fldList_mod     , only: shr_nuopc_fldList_GetNumFlds
  use shr_nuopc_fldList_mod     , only: shr_nuopc_fldList_GetFldInfo
  use shr_nuopc_methods_mod     , only: shr_nuopc_methods_FB_Init
  use shr_nuopc_methods_mod     , only: shr_nuopc_methods_FB_Reset
  use shr_nuopc_methods_mod     , only: shr_nuopc_methods_FB_Clean
  use shr_nuopc_methods_mod     , only: shr_nuopc_methods_FB_Copy
  use shr_nuopc_methods_mod     , only: shr_nuopc_methods_FB_GetFldPtr
  use shr_nuopc_methods_mod     , only: shr_nuopc_methods_Field_GeomPrint
  use shr_nuopc_methods_mod     , only: shr_nuopc_methods_State_GeomPrint
  use shr_nuopc_methods_mod     , only: shr_nuopc_methods_State_GeomWrite
  use shr_nuopc_methods_mod     , only: shr_nuopc_methods_State_reset
  use shr_nuopc_methods_mod     , only: shr_nuopc_methods_State_getNumFields
  use shr_nuopc_methods_mod     , only: shr_nuopc_methods_State_Diagnose
  use shr_nuopc_methods_mod     , only: shr_nuopc_methods_clock_timeprint
  use shr_nuopc_methods_mod     , only: shr_nuopc_methods_ChkErr
  use med_infodata_mod          , only: med_infodata_CopyStateToInfodata
  use med_infodata_mod          , only: med_infodata
  use med_internalstate_mod     , only: InternalState, llogunit=>logunit
  use med_internalstate_mod     , only: med_coupling_allowed
  use med_connectors_mod        , only: med_connectors_prep_med2atm
  use med_connectors_mod        , only: med_connectors_prep_med2ocn
  use med_connectors_mod        , only: med_connectors_prep_med2ice
  use med_connectors_mod        , only: med_connectors_prep_med2lnd
  use med_connectors_mod        , only: med_connectors_prep_med2rof
  use med_connectors_mod        , only: med_connectors_prep_med2wav
  use med_connectors_mod        , only: med_connectors_prep_med2glc
  use med_connectors_mod        , only: med_connectors_post_atm2med
  use med_connectors_mod        , only: med_connectors_post_ocn2med
  use med_connectors_mod        , only: med_connectors_post_ice2med
  use med_connectors_mod        , only: med_connectors_post_lnd2med
  use med_connectors_mod        , only: med_connectors_post_rof2med
  use med_connectors_mod        , only: med_connectors_post_wav2med
  use med_connectors_mod        , only: med_connectors_post_glc2med
  use med_phases_mod            , only: med_phases_init 
  use med_phases_prep_ocn_mod   , only: med_phases_prep_ocn_map
  use med_phases_prep_ocn_mod   , only: med_phases_prep_ocn_merge
  use med_phases_prep_ocn_mod   , only: med_phases_prep_ocn_accum_fast
  use med_phases_prep_ocn_mod   , only: med_phases_prep_ocn_accum_avg
  use med_phases_prep_atm_mod   , only: med_phases_prep_atm
  use med_phases_prep_ice_mod   , only: med_phases_prep_ice
  use med_phases_prep_lnd_mod   , only: med_phases_prep_lnd
  use med_phases_prep_rof_mod   , only: med_phases_prep_rof
  use med_phases_prep_wav_mod   , only: med_phases_prep_wav
  use med_phases_prep_glc_mod   , only: med_phases_prep_glc
  use med_phases_ocnalb_mod     , only: med_phases_ocnalb_init 
  use med_phases_ocnalb_mod     , only: med_phases_ocnalb_run
  use med_phases_aofluxes_mod   , only: med_phases_aofluxes_init 
  use med_phases_aofluxes_mod   , only: med_phases_aofluxes_run
  use med_phases_history_mod    , only: med_phases_history
  use med_fraction_mod          , only: med_fraction_init, med_fraction_set
  use med_constants_mod         , only: med_constants_dbug_flag
  use med_constants_mod         , only: med_constants_spval_init
  use med_constants_mod         , only: med_constants_spval
  use med_constants_mod         , only: med_constants_czero
  use med_constants_mod         , only: med_constants_ispval_mask
  use med_constants_mod         , only: med_constants_spval_rhfile
  use med_map_mod               , only: med_map_RouteHandles_init
  use med_map_mod               , only: med_map_MapNorm_init
  use med_io_mod                , only: med_io_cpl_init

  implicit none
  private

  integer            :: dbrc
  integer            :: stat
  character(len=1024):: msgString
  type(ESMF_VM)      :: vm
  integer            :: localPet
  logical            :: mastertask
  integer            :: dbug_flag = med_constants_dbug_flag

  character(len=*)  , parameter :: grid_arbopt = "grid_reg"   ! grid_reg or grid_arb
  real(ESMF_KIND_R8), parameter :: spval_init  = med_constants_spval_init
  real(ESMF_KIND_R8), parameter :: spval       = med_constants_spval
  real(ESMF_KIND_R8), parameter :: czero       = med_constants_czero
  integer           , parameter :: ispval_mask = med_constants_ispval_mask
  character(*)      , parameter :: u_FILE_u    = __FILE__

  public  SetServices

  private InitializeP0
  private InitializeIPDv03p1 ! advertise fields
  private InitializeIPDv03p3 ! realize connected Fields with transfer action "provide"
  private InitializeIPDv03p4 ! optionally modify the decomp/distr of transferred Grid/Mesh
  private InitializeIPDv03p5 ! realize all Fields with transfer action "accept"
  private DataInitialize     ! finish initialization and resolve data dependencies
  private SetRunClock

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    !------------------
    ! the NUOPC model component mediator_routine_SS will register the generic methods
    !------------------

    call NUOPC_CompDerive(gcomp, mediator_routine_SS, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! set entry point for methods that require specific implementation
    ! Provide InitializeP0 to switch from default IPDv00 to IPDv03
    !------------------

    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         InitializeP0, phase=0, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! IPDv03p1: advertise Fields
    !------------------

    ! Mediator advertises its import and export Fields and sets the TransferOfferGeomObject Attribute.
    ! The TransferOfferGeomObject is a String value indicating a component's
    ! intention to transfer the underlying Grid or Mesh on which an advertised Field object is defined.
    ! The valid values are: [will provide, can provide, cannot provide]

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeIPDv03p1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! IPDv03p3: realize connected Fields with transfer action "provide"
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeIPDv03p3, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! IPDv03p4: optionally modify the decomp/distr of transferred Grid/Mesh
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p4"/), userRoutine=InitializeIPDv03p4, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! IPDv03p5: realize all Fields with transfer action "accept"
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p5"/), userRoutine=InitializeIPDv03p5, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! attach specializing method for DataInitialize
    !------------------

    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_DataInitialize, &
         specRoutine=DataInitialize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! setup mediator history phase
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_history"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_history", specRoutine=med_phases_history, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep and post phases for connectors
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_connectors_prep_med2atm"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_connectors_prep_med2atm", specRoutine=med_connectors_prep_med2atm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_connectors_post_atm2med"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_connectors_post_atm2med", specRoutine=med_connectors_post_atm2med, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_connectors_prep_med2ocn"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_connectors_prep_med2ocn", specRoutine=med_connectors_prep_med2ocn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_connectors_post_ocn2med"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_connectors_post_ocn2med", specRoutine=med_connectors_post_ocn2med, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
      phaseLabelList=(/"med_connectors_prep_med2ice"/), &
      userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
      specPhaseLabel="med_connectors_prep_med2ice", specRoutine=med_connectors_prep_med2ice, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_connectors_post_ice2med"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_connectors_post_ice2med", specRoutine=med_connectors_post_ice2med, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_connectors_prep_med2lnd"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_connectors_prep_med2lnd", specRoutine=med_connectors_prep_med2lnd, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_connectors_post_lnd2med"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_connectors_post_lnd2med", specRoutine=med_connectors_post_lnd2med, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_connectors_prep_med2rof"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_connectors_prep_med2rof", specRoutine=med_connectors_prep_med2rof, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
      phaseLabelList=(/"med_connectors_post_rof2med"/), &
      userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
      specPhaseLabel="med_connectors_post_rof2med", specRoutine=med_connectors_post_rof2med, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_connectors_prep_med2wav"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_connectors_prep_med2wav", specRoutine=med_connectors_prep_med2wav, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_connectors_post_wav2med"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
      specPhaseLabel="med_connectors_post_wav2med", specRoutine=med_connectors_post_wav2med, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_connectors_prep_med2glc"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_connectors_prep_med2glc", specRoutine=med_connectors_prep_med2glc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_connectors_post_glc2med"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_connectors_post_glc2med", specRoutine=med_connectors_post_glc2med, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep routines for atm
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_atm"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_atm", specRoutine=med_phases_prep_atm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep routines for ocn
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_ocn_map"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_ocn_map", specRoutine=med_phases_prep_ocn_map, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_ocn_merge"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_ocn_merge", specRoutine=med_phases_prep_ocn_merge, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_ocn_accum_fast"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_ocn_accum_fast", specRoutine=med_phases_prep_ocn_accum_fast, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_ocn_accum_avg"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_ocn_accum_avg", specRoutine=med_phases_prep_ocn_accum_avg, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep routines for ice
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_ice"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_ice", specRoutine=med_phases_prep_ice, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep routines for lnd
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_lnd"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_lnd", specRoutine=med_phases_prep_lnd, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep routines for rof
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_rof"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_rof", specRoutine=med_phases_prep_rof, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep routines for wav
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_wav"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_wav", specRoutine=med_phases_prep_wav, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! prep routines for glc
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_prep_glc"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_prep_glc", specRoutine=med_phases_prep_glc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! phase routine for ocean albedo computation
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_ocnalb_run"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_ocnalb_run", specRoutine=med_phases_ocnalb_run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! phase routine for ocn/atm flux computation 
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_phases_aofluxes_run"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_phases_aofluxes_run", specRoutine=med_phases_aofluxes_run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! phase routine for updating fractions
    !------------------

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         phaseLabelList=(/"med_fraction_set"/), userRoutine=mediator_routine_Run, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_Advance, &
         specPhaseLabel="med_fraction_set", specRoutine=med_fraction_set, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! attach specializing method(s)
    ! -> NUOPC specializes by default --->>> first need to remove the default
    !------------------

    call ESMF_MethodRemove(gcomp, mediator_label_CheckImport, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_CheckImport, specRoutine=NUOPC_NoOp, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! attach specializing method(s)
    ! -> NUOPC specializes by default --->>> first need to remove the default
    !------------------

    call ESMF_MethodRemove(gcomp, mediator_label_SetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=mediator_label_SetRunClock, specRoutine=SetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    ! local variables
    character(len=*),parameter :: subname='(module_MED:InitializeP0)'
    character(len=128)         :: value
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    mastertask = .false.
    if (localPet == 0) mastertask=.true.

    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, defaultValue="max", &
         convention="NUOPC", purpose="Instance", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": Mediator verbosity is "//trim(value), ESMF_LOGMSG_INFO, rc=dbrc)

    dbug_flag = ESMF_UtilString2Int(value, &
         specialStringList=(/"min","max","high"/), specialValueList=(/0,255,255/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    write(msgString,'(A,i6)') trim(subname)//' dbug_flag = ',dbug_flag
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    ! Switch to IPDv03 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv03p"/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine InitializeP0

  !-----------------------------------------------------------------------

  subroutine InitializeIPDv03p1(gcomp, importState, exportState, clock, rc)

    ! Mediator advertises its import and export Fields and sets the
    ! TransferOfferGeomObject Attribute.

    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=SHR_KIND_CS) :: stdname, shortname
    logical                    :: activefld
    integer                    :: n, n1, n2, ncomp, nflds
    character(len=SHR_KIND_CS) :: transferOffer
    type(InternalState)        :: is_local
    character(len=*),parameter :: subname='(module_MED:InitializeIPDv03p1)'
    !-----------------------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    !------------------
    ! Allocate memory for the internal state and set it in the Component.
    !------------------

    allocate(is_local%wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
         msg="Allocation of the internal state memory failed.", line=__LINE__, file=u_FILE_u)) then
       return  ! bail out
    end if

    call ESMF_GridCompSetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! add a namespace (i.e. nested state)  for each import and export component state in the mediator's InternalState
    !------------------

    ! Namespaces are implemented via nested states. This creates a nested state inside of
    ! state. The nested state is returned as nestedState. nestedStateName will be used to name the
    ! newly created nested state.

    call NUOPC_AddNamespace(importState, namespace="ATM", nestedStateName="NestedState-AtmImp", &
         nestedState=is_local%wrap%NStateImp(compatm), rc=rc)
    call NUOPC_AddNamespace(importState, namespace="OCN", nestedStateName="NestedState-OcnImp", &
         nestedState=is_local%wrap%NStateImp(compocn), rc=rc)
    call NUOPC_AddNamespace(importState, namespace="ICE", nestedStateName="NestedState-IceImp", &
         nestedState=is_local%wrap%NStateImp(compice), rc=rc)
    call NUOPC_AddNamespace(importState, namespace="LND", nestedStateName="NestedState-LndImp", &
         nestedState=is_local%wrap%NStateImp(complnd), rc=rc)
    call NUOPC_AddNamespace(importState, namespace="ROF", nestedStateName="NestedState-RofImp", &
         nestedState=is_local%wrap%NStateImp(comprof), rc=rc)
    call NUOPC_AddNamespace(importState, namespace="WAV", nestedStateName="NestedState-WavImp", &
         nestedState=is_local%wrap%NStateImp(compwav), rc=rc)
    call NUOPC_AddNamespace(importState, namespace="GLC", nestedStateName="NestedState-GlcImp", &
         nestedState=is_local%wrap%NStateImp(compglc), rc=rc)
    call NUOPC_AddNamespace(exportState, namespace="ATM", nestedStateName="NestedState-AtmExp", &
         nestedState=is_local%wrap%NStateExp(compatm), rc=rc)
    call NUOPC_AddNamespace(exportState, namespace="OCN", nestedStateName="NestedState-OcnExp", &
         nestedState=is_local%wrap%NStateExp(compocn), rc=rc)
    call NUOPC_AddNamespace(exportState, namespace="ICE", nestedStateName="NestedState-IceExp", &
         nestedState=is_local%wrap%NStateExp(compice), rc=rc)
    call NUOPC_AddNamespace(exportState, namespace="LND", nestedStateName="NestedState-LndExp", &
         nestedState=is_local%wrap%NStateExp(complnd), rc=rc)
    call NUOPC_AddNamespace(exportState, namespace="ROF", nestedStateName="NestedState-RofExp", &
         nestedState=is_local%wrap%NStateExp(comprof), rc=rc)
    call NUOPC_AddNamespace(exportState, namespace="WAV", nestedStateName="NestedState-WavExp", &
         nestedState=is_local%wrap%NStateExp(compwav), rc=rc)
    call NUOPC_AddNamespace(exportState, namespace="GLC", nestedStateName="NestedState-GlcExp", &
         nestedState=is_local%wrap%NStateExp(compglc), rc=rc)

    !------------------
    ! Advertise import/export mediator field names
    !------------------

    do ncomp = 1,ncomps
       if (ncomp /= compmed) then
          nflds = shr_nuopc_fldList_GetNumFlds(fldListFr(ncomp))
          do n = 1,nflds
             call shr_nuopc_fldList_GetFldInfo(fldListFr(ncomp), n, activefld, stdname, shortname)
             ! Skip if field is not active
             if (activefld) then
                if (trim(shortname) == flds_scalar_name) then
                   transferOffer = 'will provide'
                else
                   transferOffer = 'cannot provide'
                end if
                call NUOPC_Advertise(is_local%wrap%NStateImp(ncomp), standardName=stdname, shortname=shortname, name=shortname, &
                     TransferOfferGeomObject=transferOffer)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                call ESMF_LogWrite(subname//':Fr_'//trim(compname(ncomp))//': '//trim(shortname), ESMF_LOGMSG_INFO)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end do

          nflds = shr_nuopc_fldList_GetNumFlds(fldListTo(ncomp))
          do n = 1,nflds
             call shr_nuopc_fldList_GetFldInfo(fldListTo(ncomp), n, activefld, stdname, shortname)
             ! Skip if field is not active
             if (activefld) then
                if (trim(shortname) == flds_scalar_name) then
                   transferOffer = 'will provide'
                else
                   transferOffer = 'cannot provide'
                end if
                call NUOPC_Advertise(is_local%wrap%NStateExp(ncomp), standardName=stdname, shortname=shortname, name=shortname, &
                     TransferOfferGeomObject=transferOffer)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                call ESMF_LogWrite(subname//':To_'//trim(compname(ncomp))//': '//trim(shortname), ESMF_LOGMSG_INFO)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end do
       end if
    end do ! end of ncomps loop

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine InitializeIPDv03p1

  !-----------------------------------------------------------------------------

  subroutine InitializeIPDv03p3(gcomp, importState, exportState, clock, rc)

    ! Realize connected Fields with transfer action "provide"

    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    integer                         :: i, j
    real(kind=ESMF_KIND_R8),pointer :: lonPtr(:), latPtr(:)
    type(InternalState)             :: is_local
    integer                         :: lmpicom
    real(ESMF_KIND_R8)              :: intervalSec
    type(ESMF_TimeInterval)         :: timeStep
    ! tcx XGrid
    ! type(ESMF_Field)              :: fieldX, fieldA, fieldO
    ! type(ESMF_XGrid)              :: xgrid
    integer                         :: n, n1, n2
    character(SHR_KIND_CL)          :: cvalue
    logical                         :: connected
    character(len=*),parameter      :: subname='(module_MED:InitializeIPDv03p3)'
    !-----------------------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize the internal state members
    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call MPI_Comm_Dup(lmpicom, is_local%wrap%mpicom, stat)

    ! Realize States
    do n = 1,ncomps
      if (ESMF_StateIsCreated(is_local%wrap%NStateImp(n), rc=rc)) then
         call shr_nuopc_fldList_Realize(is_local%wrap%NStateImp(n), fldListFr(n), flds_scalar_name, flds_scalar_num, &
              tag=subname//':Fr_'//trim(compname(n)), rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
      if (ESMF_StateIsCreated(is_local%wrap%NStateExp(n), rc=rc)) then
         call shr_nuopc_fldList_Realize(is_local%wrap%NStateExp(n), fldListTo(n), flds_scalar_name, flds_scalar_num, &
              tag=subname//':To_'//trim(compname(n)), rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
    enddo

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine InitializeIPDv03p3

  !-----------------------------------------------------------------------------

  subroutine InitializeIPDv03p4(gcomp, importState, exportState, clock, rc)

    ! Optionally modify the decomp/distr of transferred Grid/Mesh

    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer :: n1,n2
    !    type(ESMF_Field)              :: field
    !    type(ESMF_Grid)               :: grid
    !    integer                       :: localDeCount
    !    type(ESMF_DistGrid)           :: distgrid
    !    integer                       :: dimCount, tileCount, petCount
    !    integer                       :: deCountPTile, extraDEs
    !    integer, allocatable          :: minIndexPTile(:,:), maxIndexPTile(:,:)
    !    integer, allocatable          :: regDecompPTile(:,:)
    !    integer                       :: i, j, n, n1
    character(len=*),parameter :: subname='(module_MED:realizeConnectedGrid)'
    !-----------------------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! Get the internal state from the mediator gridded component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! Recieve Grids
    !------------------

    do n1 = 1,ncomps
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": calling for component "//trim(compname(n1)), ESMF_LOGMSG_INFO, rc=dbrc)
       end if
       if (ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then
          call realizeConnectedGrid(is_local%wrap%NStateImp(n1), trim(compname(n1))//'Imp', rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
       if (ESMF_StateIsCreated(is_local%wrap%NStateExp(n1),rc=rc)) then
          call realizeConnectedGrid(is_local%wrap%NStateExp(n1), trim(compname(n1))//'Exp', rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": finished for component "//trim(compname(n1)), ESMF_LOGMSG_INFO, rc=dbrc)
       end if
    enddo

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine realizeConnectedGrid(State,string,rc)

      type(ESMF_State)   , intent(inout) :: State
      character(len=*)   , intent(in)    :: string
      integer            , intent(out)   :: rc

      ! local variables
      type(ESMF_Field)              :: field
      type(ESMF_Grid)               :: grid
      integer                       :: localDeCount

      type(ESMF_DistGrid)           :: distgrid
      type(ESMF_DistGridConnection), allocatable :: connectionList(:)
      integer                       :: arbDimCount
      integer                       :: dimCount, tileCount, petCount
      integer                       :: connectionCount
      integer                       :: deCountPTile, extraDEs
      integer, allocatable          :: minIndexPTile(:,:), maxIndexPTile(:,:)
      integer, allocatable          :: regDecompPTile(:,:)
      integer                       :: i, j, n, n1, fieldCount, nxg, i1, i2
      type(ESMF_GeomType_Flag)      :: geomtype
      character(ESMF_MAXSTR),allocatable :: fieldNameList(:)
      type(ESMF_FieldStatus_Flag)   :: fieldStatus
      character(len=*),parameter :: subname='(module_MEDIATOR:realizeConnectedGrid)'

      !NOTE: All of the Fields that set their TransferOfferGeomObject Attribute
      !NOTE: to "cannot provide" should now have the accepted Grid available.
      !NOTE: Go and pull out this Grid for one of a representative Field and
      !NOTE: modify the decomposition and distribution of the Grid to match the
      !NOTE: Mediator PETs.

      !TODO: quick implementation, do it for each field one by one
      !TODO: commented out below are application to other fields

      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
      endif
      rc = ESMF_Success

      call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      allocate(fieldNameList(fieldCount))
      call ESMF_StateGet(State, itemNameList=fieldNameList, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !     do n=1, fieldCount
      do n=1, min(fieldCount,1)

         call ESMF_StateGet(State, field=field, itemName=fieldNameList(n), rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

         call ESMF_FieldGet(field, status=fieldStatus, rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

         if (fieldStatus==ESMF_FIELDSTATUS_GRIDSET) then

            ! while this is still an empty field, it does now hold a Grid with DistGrid
            call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            if (geomtype == ESMF_GEOMTYPE_GRID) then

               call shr_nuopc_methods_Field_GeomPrint(field,trim(fieldNameList(n))//'_orig',rc)
               if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

               call ESMF_AttributeGet(field, name="ArbDimCount", value=arbDimCount, &
                    convention="NUOPC", purpose="Instance", rc=rc)
               if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

               if (dbug_flag > 1) then
                  call ESMF_LogWrite(trim(subname)//": geomtype is ESMF_GEOMTYPE_GRID for "//trim(fieldnameList(n)), &
                       ESMF_LOGMSG_INFO, rc=dbrc)
                  write(msgString,'(A,i8)') trim(subname)//':arbdimcount =',arbdimcount
                  call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=dbrc)
               endif

               ! make decision on whether the incoming Grid is arbDistr or not
               if (arbDimCount>0) then
                  ! The provider defined an arbDistr grid
                  !
                  ! Need to make a choice here to either represent the grid as a
                  ! regDecomp grid on the acceptor side, or to stay with arbDistr grid:
                  !
                  ! Setting the PRECIP_REGDECOMP macro will set up a regDecomp grid on the
                  ! acceptor side.
                  !
                  ! Not setting the PRECIP_REGDECOMP macro will default into keeping the
                  ! original arbDistr Grid.

                  if (grid_arbopt == "grid_reg") then

                     if (dbug_flag > 1) then
                        call ESMF_LogWrite(trim(subname)//trim(string)//": accept arb2reg grid for "//trim(fieldNameList(n)), &
                             ESMF_LOGMSG_INFO, rc=dbrc)
                     endif

                     ! Use a regDecomp representation for the grid
                     ! first get tile min/max, only single tile supported for arbDistr Grid
                     allocate(minIndexPTile(arbDimCount,1),maxIndexPTile(arbDimCount,1))
                     call ESMF_AttributeGet(field, name="MinIndex", &
                          valueList=minIndexPTile(:,1), &
                          convention="NUOPC", purpose="Instance", rc=rc)
                     if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                     call ESMF_AttributeGet(field, name="MaxIndex", &
                          valueList=maxIndexPTile(:,1), &
                          convention="NUOPC", purpose="Instance", rc=rc)
                     if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                     ! create default regDecomp DistGrid
                     distgrid = ESMF_DistGridCreate(minIndexPTile=minIndexPTile, &
                          maxIndexPTile=maxIndexPTile, connectionList=connectionList, rc=rc)
                     if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                     ! Create default regDecomp Grid
                     grid = ESMF_GridCreate(distgrid, rc=rc)
                     if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                     ! swap out the transferred grid for the newly created one
                     call ESMF_FieldEmptySet(field, grid=grid, rc=rc)
                     if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                     do i1 = 1,arbDimCount
                        write(msgString,'(A,3i8)') trim(subname)//':PTile =',i1,minIndexPTile(i1,1),maxIndexPTile(i1,1)
                        call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
                     enddo
                     deallocate(minIndexPTile,maxIndexPTile)

                  elseif (grid_arbopt == "grid_arb") then

                     ! Stick with the arbDistr representation of the grid:
                     ! There is nothing to do here if the same number of DEs is kept on the
                     ! acceptor side. Alternatively, the acceptor side could set up a more
                     ! natural number of DEs (maybe same number as acceptor PETs), and then
                     ! redistribute the arbSeqIndexList. Here simply keep the DEs of the
                     ! provider Grid.
                     if (dbug_flag > 1) then
                        call ESMF_LogWrite(trim(subname)//trim(string)//": accept arb2arb grid for "//trim(fieldNameList(n)), &
                             ESMF_LOGMSG_INFO, rc=dbrc)
                     endif

                  else   ! grid_arbopt

                     call ESMF_LogWrite(trim(subname)//trim(string)//": ERROR grid_arbopt setting = "//trim(grid_arbopt), &
                          ESMF_LOGMSG_INFO, rc=rc)
                     rc = ESMF_FAILURE
                     if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                  endif  ! grid_arbopt


               else   ! arbdimcount <= 0

                  ! The provider defined as non arb grid

                  ! access localDeCount to show this is a real Grid
                  if (dbug_flag > 1) then
                     call ESMF_LogWrite(trim(subname)//trim(string)//": accept reg2reg grid for "//&
                          trim(fieldNameList(n)), ESMF_LOGMSG_INFO, rc=dbrc)
                  endif

                  call ESMF_FieldGet(field, grid=grid, rc=rc)
                  if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                  call ESMF_GridGet(grid, localDeCount=localDeCount, distgrid=distgrid, rc=rc)
                  if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                  ! Create a custom DistGrid, based on the minIndex, maxIndex of the
                  ! accepted DistGrid, but with a default regDecomp for the current VM
                  ! that leads to 1DE/PET.

                  ! get dimCount and tileCount
                  call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, &
                       connectionCount=connectionCount, rc=rc)
                  if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                  ! allocate minIndexPTile and maxIndexPTile accord. to dimCount and tileCount
                  allocate(minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount))
                  allocate(connectionList(connectionCount))

                  ! get minIndex and maxIndex arrays, and connectionList
                  call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
                       maxIndexPTile=maxIndexPTile, connectionList=connectionList, rc=rc)
                  if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                  ! construct a default regDecompPTile -> TODO: move this into ESMF as default
                  call ESMF_GridCompGet(gcomp, petCount=petCount, rc=rc)
                  if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                  allocate(regDecompPTile(dimCount, tileCount))
                  deCountPTile = petCount/tileCount
                  extraDEs = max(0, petCount-deCountPTile)
                  do i=1, tileCount
                     if (i<=extraDEs) then
                        regDecompPTile(1, i) = deCountPTile + 1
                     else
                        regDecompPTile(1, i) = deCountPTile
                     endif
                     do j=2, dimCount
                        regDecompPTile(j, i) = 1
                     enddo
                  enddo

                  do i2 = 1,tileCount
                     do i1 = 1,dimCount
                        write(msgString,'(A,5i8)') trim(subname)//':PTile =',i2,i1,minIndexPTile(i1,i2),&
                             maxIndexPTile(i1,i2),regDecompPTile(i1,i2)
                        call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
                     enddo
                  enddo

                  !--- tcraig, hardwire i direction wraparound, temporary
                  !--- tcraig, now getting info from model distgrid, see above
                  !              allocate(connectionList(1))
                  !              nxg = maxIndexPTile(1,1) - minIndexPTile(1,1) + 1
                  !              write(msgstring,*) trim(subname)//trim(string),': connlist nxg = ',nxg
                  !              call ESMF_LogWrite(trim(msgstring), ESMF_LOGMSG_INFO, rc=rc)
                  !              if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                  !              call ESMF_DistGridConnectionSet(connectionList(1), tileIndexA=1, &
                  !                tileIndexB=1, positionVector=(/nxg, 0/), rc=rc)
                  !              if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                  ! create the new DistGrid with the same minIndexPTile and maxIndexPTile,
                  ! but with a default regDecompPTile
                  ! tcraig, force connectionlist and gridEdge arguments to fix wraparound
                  ! need ESMF fixes to implement properly.
                  if (dimcount == 2) then
                     distgrid = ESMF_DistGridCreate(minIndexPTile=minIndexPTile, &
                          maxIndexPTile=maxIndexPTile, regDecompPTile=regDecompPTile, &
                          connectionList=connectionList, rc=rc)
                     if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                     if (dbug_flag > 1) then
                        call ESMF_LogWrite(trim(subname)//trim(string)//': distgrid with dimcount=2', ESMF_LOGMSG_INFO, rc=rc)
                        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                     endif

                     ! Create a new Grid on the new DistGrid and swap it in the Field
                     grid = ESMF_GridCreate(distgrid, gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), rc=rc)
                     if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                  else
                     distgrid = ESMF_DistGridCreate(minIndexPTile=minIndexPTile, &
                          maxIndexPTile=maxIndexPTile, regDecompPTile=regDecompPTile, rc=rc)
                     if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                     if (dbug_flag > 1) then
                        call ESMF_LogWrite(trim(subname)//trim(string)//': distgrid with dimcount=1', ESMF_LOGMSG_INFO, rc=rc)
                        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                     endif

                     ! Create a new Grid on the new DistGrid and swap it in the Field
                     grid = ESMF_GridCreate(distgrid, gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/), rc=rc)
                     if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                  endif

                  ! local clean-up
                  deallocate(connectionList)
                  deallocate(minIndexPTile, maxIndexPTile, regDecompPTile)

               endif  ! arbdimCount

               ! Swap all the Grids in the State

               ! do n1=n,n
               do n1=1, fieldCount
                  ! access a field in the State and set the Grid
                  call ESMF_StateGet(State, field=field, itemName=fieldNameList(n1), rc=rc)
                  if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                  call ESMF_FieldGet(field, status=fieldStatus, rc=rc)
                  if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                  if (fieldStatus==ESMF_FIELDSTATUS_EMPTY) then
                     call ESMF_FieldEmptySet(field, grid=grid, rc=rc)
                     if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                  endif

                  if (dbug_flag > 1) then
                     call ESMF_LogWrite(trim(subname)//trim(string)//": attach grid for "//trim(fieldNameList(n1)), &
                          ESMF_LOGMSG_INFO, rc=rc)
                     if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                  endif

                  call shr_nuopc_methods_Field_GeomPrint(field,trim(fieldNameList(n1))//'_new',rc)
                  if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
               enddo

            elseif (geomtype == ESMF_GEOMTYPE_MESH) then

               if (dbug_flag > 1) then
                  call ESMF_LogWrite(trim(subname)//": geomtype is ESMF_GEOMTYPE_MESH for "//trim(fieldnameList(n)), &
                       ESMF_LOGMSG_INFO, rc=dbrc)
               end if

               call shr_nuopc_methods_Field_GeomPrint(field,trim(fieldNameList(n))//'_orig',rc)
               if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            else  ! geomtype

               call ESMF_LogWrite(trim(subname)//": ERROR geomtype not supported ", ESMF_LOGMSG_INFO, rc=rc)
               rc=ESMF_FAILURE
               return

            endif ! geomtype

         elseif (fieldStatus==ESMF_FIELDSTATUS_EMPTY) then

            call ESMF_LogWrite(trim(subname)//trim(string)//": provide grid for "//trim(fieldNameList(n)), &
                 ESMF_LOGMSG_INFO, rc=dbrc)

         elseif (fieldStatus==ESMF_FIELDSTATUS_COMPLETE) then

            call ESMF_LogWrite(trim(subname)//trim(string)//": no grid provided for "//trim(fieldNameList(n)), &
                 ESMF_LOGMSG_INFO, rc=dbrc)

         else

            call ESMF_LogWrite(trim(subname)//": ERROR fieldStatus not supported ", ESMF_LOGMSG_INFO, rc=rc)
            rc=ESMF_FAILURE
            return

         endif   ! fieldStatus

      enddo   ! nflds

      deallocate(fieldNameList)

      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
      endif

    end subroutine realizeConnectedGrid

  end subroutine InitializeIPDv03p4

  !-----------------------------------------------------------------------------

  subroutine InitializeIPDv03p5(gcomp, importState, exportState, clock, rc)

    !----------------------------------------------------------
    ! realize all Fields with transfer action "accept"
    !----------------------------------------------------------

    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n1,n2
    character(len=*),parameter  :: subname='(module_MED:InitializeIPDv03p5)'
    !-----------------------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    rc = ESMF_SUCCESS

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--- Finish initializing the State Fields
    !--- Write out grid information

    do n1 = 1,ncomps

      if (ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then
        if (dbug_flag > 5) then
           call ESMF_LogWrite(trim(subname)//": calling completeFieldInitialize import states from "//trim(compname(n1)), &
                ESMF_LOGMSG_INFO, rc=dbrc)
        end if
        call completeFieldInitialization(is_local%wrap%NStateImp(n1), rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        call shr_nuopc_methods_State_reset(is_local%wrap%NStateImp(n1), value=spval_init, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif

      if (ESMF_StateIsCreated(is_local%wrap%NStateExp(n1),rc=rc)) then
        if (dbug_flag > 5) then
           call ESMF_LogWrite(trim(subname)//": calling completeFieldInitialize export states to "//trim(compname(n1)), &
                ESMF_LOGMSG_INFO, rc=dbrc)
        end if
        call completeFieldInitialization(is_local%wrap%NStateExp(n1), rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        call shr_nuopc_methods_State_reset(is_local%wrap%NStateExp(n1), value=spval_init, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        call shr_nuopc_methods_State_GeomPrint(is_local%wrap%NStateExp(n1),'gridExp'//trim(compname(n1)),rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        call shr_nuopc_methods_State_GeomWrite(is_local%wrap%NStateExp(n1), 'grid_med_'//trim(compname(n1)), rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
    enddo

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine completeFieldInitialization(State,rc)

      type(ESMF_State)   , intent(inout) :: State
      integer            , intent(out)   :: rc

      integer                     :: n, fieldCount
      character(ESMF_MAXSTR)      :: fieldName
      type(ESMF_Grid)             :: grid
      type(ESMF_Mesh)             :: mesh
      type(ESMF_Field)            :: meshField
      type(ESMF_Field),pointer    :: fieldList(:)
      type(ESMF_FieldStatus_Flag) :: fieldStatus
      type(ESMF_GeomType_Flag)    :: geomtype
      character(len=*),parameter  :: subname='(module_MED:completeFieldInitialization)'

      if (dbug_flag > 5) then
        call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
      endif

      rc = ESMF_Success

      call shr_nuopc_methods_State_GetNumFields(State, fieldCount, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      if (fieldCount > 0) then
        nullify(fieldList)
        call NUOPC_getStateMemberLists(State, fieldList=fieldList, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        do n=1, fieldCount

          call ESMF_FieldGet(fieldList(n), status=fieldStatus, name=fieldName, &
            geomtype=geomtype, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          if (geomtype == ESMF_GEOMTYPE_GRID .and. fieldName /= flds_scalar_name) then
            ! Grab grid
            call shr_nuopc_methods_Field_GeomPrint(fieldList(n),trim(fieldName)//'_premesh',rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
            call ESMF_FieldGet(fieldList(n), grid=grid, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            ! Convert grid to mesh
            mesh = ESMF_GridToMeshCell(grid,rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            meshField = ESMF_FieldCreate(mesh, typekind=ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, name=fieldName, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            ! Swap grid for mesh, at this point, only connected fields are in the state
            call NUOPC_Realize(State, field=meshField, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif

          if (fieldStatus==ESMF_FIELDSTATUS_GRIDSET) then
            if (dbug_flag > 1) then
              call ESMF_LogWrite(subname//" is allocating field memory for field "//trim(fieldName), &
                   ESMF_LOGMSG_INFO, rc=rc)
              if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
            endif
            call ESMF_FieldEmptyComplete(fieldList(n), typekind=ESMF_TYPEKIND_R8, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif   ! fieldStatus

          call shr_nuopc_methods_Field_GeomPrint(fieldList(n), trim(subname)//':'//trim(fieldName), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        enddo
        deallocate(fieldList)
      endif

      if (dbug_flag > 5) then
        call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
      endif

    end subroutine completeFieldInitialization

  end subroutine InitializeIPDv03p5

  !-----------------------------------------------------------------------------

  subroutine DataInitialize(gcomp, rc)

    !----------------------------------------------------------
    ! Finish initialization and resolve data dependencies
    ! There will be multiple passes
    ! For first time through:
    !   Do not assume any import fields are connected, just allocate space and such
    !   -- Check present flags
    !   -- Check for active coupling interactions
    !   -- Initialize connector count arrays in med_internal_state
    !   -- Create FBs: FBImp, FBExp, FBImpAccum, FBExpAccum
    !   -- Create mediator specific field bundles (not part of import/export states)
    !   -- Initialize med_infodata, Accums (to zero), and FBImp (from NStateImp)
    !   -- Read mediator restarts
    !   -- Initialize route handles 
    !   -- Initialize field bundles for normalization
    !   -- return!
    ! For second loop:
    !   -- Copy import fields to local FBs
    !   -- Create FBfrac and initialize fractions
    ! Once the ocean is ready:
    !   -- Copy import fields to local FBs
    !   -- Re-initialize fractions
    !   -- Carry out ocnalb_init 
    !   -- Carry out aoffluxes_init
    ! Once the atm is ready:
    !   -- Copy import fields to local FBs
    !----------------------------------------------------------

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)                :: is_local
    type(ESMF_Clock)                   :: clock
    type(ESMF_State)                   :: importState, exportState
    type(ESMF_Time)                    :: time
    type(ESMF_Field)                   :: field
    type(ESMF_StateItem_Flag)          :: itemType
    logical                            :: atCorrectTime, connected
    integer                            :: n1,n2,n
    integer                            :: cntn1, cntn2
    integer                            :: fieldCount
    character(SHR_KIND_CL), pointer    :: fldnames(:)
    character(ESMF_MAXSTR),allocatable :: fieldNameList(:)
    character(len=128)                 :: value
    character(SHR_KIND_CL)             :: cvalue
    logical                            :: LocalDone
    logical,save                       :: atmDone = .false.
    logical,save                       :: ocnDone = .false.
    logical,save                       :: allDone = .false.
    logical,save                       :: first_call = .true.
    character(len=*), parameter        :: subname='(module_MED:DataInitialize)'
    !-----------------------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    rc = ESMF_SUCCESS

    call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="false", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=time, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Beginning  of first_call block
    !---------------------------------------

    if (first_call) then

      ! initialize the present flags in the mediator
      if (dbug_flag > 1) then
        call ESMF_LogWrite("Starting to initialize present flags", ESMF_LOGMSG_INFO)
        call ESMF_LogFlush()
      endif

      !----------------------------------------------------------
      !--- Check present flags
      !----------------------------------------------------------

      do n1 = 1,ncomps
        call ESMF_AttributeGet(gcomp, name=trim(compname(n1))//"_present", value=value, defaultValue="false", &
             convention="NUOPC", purpose="Instance", rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        is_local%wrap%comp_present(n1) = (value == "true")
        write(msgString,'(A,L4)') trim(subname)//' comp_present(comp'//trim(compname(n1))//') = ',&
             is_local%wrap%comp_present(n1)
        call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
      enddo

      !----------------------------------------------------------
      !--- Check for active coupling interactions
      !    must be allowed, bundles created, and both sides have some fields
      !----------------------------------------------------------

      if (dbug_flag > 1) then
        call ESMF_LogWrite("Starting to initialize active flags", ESMF_LOGMSG_INFO)
        call ESMF_LogFlush()
      endif

      ! initialize med_coupling_active
      is_local%wrap%med_coupling_active(:,:) = .false.

      do n1 = 1,ncomps
        if (is_local%wrap%comp_present(n1) .and. ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then
          call shr_nuopc_methods_State_GetNumFields(is_local%wrap%NStateImp(n1), cntn1, rc=rc) ! Import Field Count
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          if (cntn1 > 0) then
             do n2 = 1,ncomps
                if (is_local%wrap%comp_present(n2) .and. ESMF_StateIsCreated(is_local%wrap%NStateExp(n2),rc=rc) &
                     .and. med_coupling_allowed(n1,n2)) then
                   call shr_nuopc_methods_State_GetNumFields(is_local%wrap%NStateExp(n2), cntn2, rc=rc) ! Import Field Count
                   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                   if (cntn2 > 0) then
                      is_local%wrap%med_coupling_active(n1,n2) = .true.
                   endif
                endif
             enddo
          end if
        endif
      enddo

      ! create tables of output
      if (mastertask) then
         if (dbug_flag > 5) then
            write(llogunit,*) ' '
            write(llogunit,'(A)') subname//' Allowed coupling flags'
            write(llogunit,'(2x,A10,20(A5))') '|from to->',(compname(n2),n2=1,ncomps)
            do n1 = 1,ncomps
               write(msgString,'(2x,a1,A,5x,20(L5))') '|',trim(compname(n1)),(med_coupling_allowed(n1,n2),n2=1,ncomps)
               do n2 = 1,len_trim(msgString)
                  if (msgString(n2:n2) == 'F') msgString(n2:n2)='-'
               enddo
               write(llogunit,'(A)') trim(msgString)
            enddo
            write(llogunit,*) ' '
            call shr_sys_flush(llogunit)
         endif

         if (dbug_flag >= 0) then
            write(llogunit,*) ' '
            write(llogunit,'(A)') subname//' Active coupling flags'
            write(llogunit,'(2x,A10,20(A5))') '|from to->',(compname(n2),n2=1,ncomps)
            do n1 = 1,ncomps
               write(msgString,'(2x,a1,A,5x,20(L5))') '|',trim(compname(n1)),&
                    (is_local%wrap%med_coupling_active(n1,n2),n2=1,ncomps)
               do n2 = 1,len_trim(msgString)
                  if (msgString(n2:n2) == 'F') msgString(n2:n2)='-'
               enddo
               write(llogunit,'(A)') trim(msgString)
            enddo
            write(llogunit,*) ' '
            call shr_sys_flush(llogunit)
         endif
      endif

      !----------------------------------------------------------
      ! Initialize connector count
      !----------------------------------------------------------

      if (dbug_flag > 1) then
        call ESMF_LogWrite("Starting to Create FBs", ESMF_LOGMSG_INFO)
        call ESMF_LogFlush()
      endif

      is_local%wrap%conn_prep_cnt(:) = 0
      is_local%wrap%conn_post_cnt(:) = 0

      !----------------------------------------------------------
      ! Create various FBs, FBImp, FBExp, FBImpAccum, FBExpAccum
      !----------------------------------------------------------

      do n1 = 1,ncomps
         if (is_local%wrap%comp_present(n1) .and. &
              ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc) .and. &
              ESMF_StateIsCreated(is_local%wrap%NStateExp(n1),rc=rc)) then

            if (mastertask) write(llogunit,*) subname,' initializing FBs for '//trim(compname(n1))

            call shr_nuopc_methods_FB_init(is_local%wrap%FBImp(n1,n1), flds_scalar_name, &
                 STgeom=is_local%wrap%NStateImp(n1), &
                 STflds=is_local%wrap%NStateImp(n1), &
                 name='FBImp'//trim(compname(n1)), rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            call shr_nuopc_methods_FB_init(is_local%wrap%FBExp(n1), flds_scalar_name, &
                 STgeom=is_local%wrap%NStateExp(n1), &
                 STflds=is_local%wrap%NStateExp(n1), &
                 name='FBExp'//trim(compname(n1)), rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            call shr_nuopc_methods_FB_init(is_local%wrap%FBImpAccum(n1), flds_scalar_name, &
                 STgeom=is_local%wrap%NStateImp(n1), &
                 STflds=is_local%wrap%NStateImp(n1), &
                 name='FBImpAccum'//trim(compname(n1)), rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            call shr_nuopc_methods_FB_init(is_local%wrap%FBExpAccum(n1), flds_scalar_name, &
                 STgeom=is_local%wrap%NStateExp(n1), &
                 STflds=is_local%wrap%NStateExp(n1), &
                 name='FBExpAccum'//trim(compname(n1)), rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            is_local%wrap%FBImpAccumCnt(n1) = 0
            is_local%wrap%FBExpAccumCnt(n1) = 0

            call shr_nuopc_methods_FB_reset(is_local%wrap%FBImpAccum(n1), value=czero, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            call shr_nuopc_methods_FB_reset(is_local%wrap%FBExpAccum(n1), value=czero, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         endif
         if (mastertask) call shr_sys_flush(llogunit)

         ! These are the FBImp mapped to different grids, FBImp(n1,n1) is handled above
         do n2 = 1,ncomps
            if (n1 /= n2 .and. &
                 is_local%wrap%med_coupling_active(n1,n2) .and. &
                 ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc) .and. &
                 ESMF_StateIsCreated(is_local%wrap%NStateExp(n2),rc=rc)) then
               if (mastertask) write(llogunit,*) subname,' initializing FBs for '//trim(compname(n1))//'_'//trim(compname(n2))

               ! TODO:
               ! The NStateImp(n2) should be used here rather than NStateExp(n2), since
               ! the export state might only contain control data and no grid information if
               ! if the target component (n2) is not prognostic only receives control data back
               ! But if STgeom=is_local%wrap%NStateImp(n2) is substituted for STgeom=is_local%wrap%NStateExp(n2) 
               ! then an error occurs as follows

               call shr_nuopc_methods_FB_init(is_local%wrap%FBImp(n1,n2), flds_scalar_name, &
                    STgeom=is_local%wrap%NStateExp(n2), &
                    STflds=is_local%wrap%NStateImp(n1), &
                    name='FBImp'//trim(compname(n1))//'_'//trim(compname(n2)), rc=rc)
               if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
            endif
         enddo
      enddo
      if (mastertask) call shr_sys_flush(llogunit)

#if (1 == 0)
      !---------------------------------------
      ! read mediator restarts
      !---------------------------------------
      !---tcraig, turn if on to force no mediator restarts for testing
      !if (.not.coldstart) then
        call Mediator_restart(gcomp,'read','mediator',rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      !endif

      ! default initialize s_surf to work around limitations of current initialization sequence
      call ESMF_StateGet(is_local%wrap%NStateExp(compice), itemName='s_surf', itemType=itemType, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (itemType /= ESMF_STATEITEM_NOTFOUND) then
        if (NUOPC_IsConnected(is_local%wrap%NStateExp(compice),'s_surf',rc=rc)) then
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call State_SetFldPtr(is_local%wrap%NStateExp(compice), 's_surf', 34.0_ESMF_KIND_R8, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
      endif
#endif

      !---------------------------------------
      !--- Initialize route handles and required normalization field bunds
      !---------------------------------------
      
      call med_map_RouteHandles_init(gcomp, llogunit, rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      
      call med_map_MapNorm_init(gcomp, llogunit, rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      
      !---------------------------------------
      ! Initialize field bundles needed for ocn albedo and ocn/atm flux calculations
      !---------------------------------------

      if (is_local%wrap%med_coupling_active(compocn,compatm) .and. &
          is_local%wrap%med_coupling_active(compatm,compocn)) then

         ! NOTE: the NStateImp(compocn) or NStateImp(compatm) used below
         ! rather than NStateExp(n2), since the export state might only
         ! contain control data and no grid information if if the target
         ! component (n2) is not prognostic only receives control data back

         ! Create field bundles for ocean albedo computation

         fieldCount = shr_nuopc_fldList_GetNumFlds(fldListMed_ocnalb_o)
         allocate(fldnames(fieldCount))
         call shr_nuopc_fldList_getfldnames(fldListMed_ocnalb_o%flds, fldnames)

         call shr_nuopc_methods_FB_init(is_local%wrap%FBMed_ocnalb_a, flds_scalar_name, &
            STgeom=is_local%wrap%NStateImp(compatm), fieldnamelist=fldnames, name='FBMed_ocnalb_a', rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

         call shr_nuopc_methods_FB_init(is_local%wrap%FBMed_ocnalb_o, flds_scalar_name, &
            STgeom=is_local%wrap%NStateImp(compocn), fieldnamelist=fldnames, name='FBMed_ocnalb_o', rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         deallocate(fldnames)

         ! Create field bundles for ocean/atmosphere flux computation

         fieldCount = shr_nuopc_fldList_GetNumFlds(fldListMed_aoflux_o)
         allocate(fldnames(fieldCount))
         call shr_nuopc_fldList_getfldnames(fldListMed_aoflux_a%flds, fldnames)

         call shr_nuopc_methods_FB_init(is_local%wrap%FBMed_aoflux_a, flds_scalar_name, &
            STgeom=is_local%wrap%NStateImp(compatm), fieldnamelist=fldnames, name='FBMed_aoflux_a', rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

         call shr_nuopc_methods_FB_init(is_local%wrap%FBMed_aoflux_o, flds_scalar_name, &
            STgeom=is_local%wrap%NStateImp(compocn), fieldnamelist=fldnames, name='FBMed_aoflux_o', rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         deallocate(fldnames)
      end if

      !----------------------------------------------------------
      ! Create mediator specific field bundles needed in phases routines
      ! TODO: this needs to be filled in
      !----------------------------------------------------------

      ! FBs for lnd <-> glc accumulation and elevation class downscaling
      if (is_local%wrap%comp_present(complnd) .and. is_local%wrap%comp_present(compglc)) then
         ! call shr_nuopc_methods_FB_init(is_local%wrap%FBMed_l2x_to_glc_accum, &
         !      STgeom=is_local%wrap%NStateImp(complnd), fieldnamelist=flds_l2x_to_glc, name='FBMed_l2g_l_accum', rc=rc)
         ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

         ! call shr_nuopc_methods_FB_init(is_local%wrap%FBMed_g2x_to_lnd, &
         !      STgeom=is_local%wrap%NStateImp(complnd), fieldnamelist=flds_g2x_to_lnd, name='FBMed_g2x_to_lnd', rc=rc)
         ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      end if

      first_call = .false.
      call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="false", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      return
    endif  ! end first_call if-block

    !---------------------------------------
    ! Initialize mediator fields and infodata
    ! This is called every loop around DataInitialize
    !---------------------------------------

    do n1 = 1,ncomps
       LocalDone = .true.
       if (is_local%wrap%comp_present(n1) .and. ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then

          call ESMF_StateGet(is_local%wrap%NStateImp(n1), itemCount=fieldCount, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
 
          allocate(fieldNameList(fieldCount))
          call ESMF_StateGet(is_local%wrap%NStateImp(n1), itemNameList=fieldNameList, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          do n=1, fieldCount
             call ESMF_StateGet(is_local%wrap%NStateImp(n1), itemName=fieldNameList(n), field=field, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             atCorrectTime = NUOPC_IsAtTime(field, time, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             if (atCorrectTime) then
                if (fieldNameList(n) == flds_scalar_name) then
                   call med_infodata_CopyStateToInfodata(is_local%wrap%NStateImp(n1), med_infodata, &
                        trim(compname(n1))//'2cpli', is_local%wrap%mpicom, rc)
                   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                   call ESMF_LogWrite(trim(subname)//" MED - Initialize-Data-Dependency CSTI "//trim(compname(n1)), &
                        ESMF_LOGMSG_INFO, rc=rc)
                   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                endif
             else
                LocalDone=.false.
             endif
          enddo
          deallocate(fieldNameList)

          if (LocalDone) then
             call shr_nuopc_methods_FB_copy(is_local%wrap%FBImp(n1,n1), is_local%wrap%NStateImp(n1), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             call ESMF_LogWrite(trim(subname)//" MED - Initialize-Data-Dependency Copy Import "//trim(compname(n1)), ESMF_LOGMSG_INFO, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       endif
    enddo

    !----------------------------------------------------------
    ! Create FBfrac field bundles and initialize fractions
    ! This has some complex dependencies on fractions from import States
    !  and appropriate checks are not implemented.  These fractions are needed
    !  also in the ocean ocnalb_init and ocnaoflux_init.  We might need to split 
    !  out the fraction FB allocation and the fraction initialization
    !----------------------------------------------------------
    
    call med_fraction_init(gcomp,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Carry out data dependency for ocn initialization if needed
    !---------------------------------------

    if (.not. is_local%wrap%comp_present(compocn)) then
       ocnDone = .true.
    endif

    if (.not. is_local%wrap%comp_present(compocn)) then
       atmDone = .true.
    endif

    if (.not. ocnDone .and. is_local%wrap%comp_present(compocn)) then

      ocnDone = .true.  ! reset if an item is found that is not done

      if (is_local%wrap%med_coupling_active(compocn,compatm) .and. &
          is_local%wrap%med_coupling_active(compatm,compocn)) then

         call ESMF_StateGet(is_local%wrap%NStateImp(compocn), itemCount=fieldCount, rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

         allocate(fieldNameList(fieldCount))
         call ESMF_StateGet(is_local%wrap%NStateImp(compocn), itemNameList=fieldNameList, rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

         do n=1, fieldCount
            call ESMF_StateGet(is_local%wrap%NStateImp(compocn), itemName=fieldNameList(n), field=field, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            atCorrectTime = NUOPC_IsAtTime(field, time, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            if (.not. atCorrectTime) then
               ! If any ocn import fields are not time stamped correctly, then dependency is not satisified - must return to ocn
               call ESMF_LogWrite("MED - Initialize-Data-Dependency from OCN NOT YET SATISFIED!!!", ESMF_LOGMSG_INFO, rc=rc)
               if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
               ocnDone = .false.
               exit  ! break out of the loop when first not satisfied found
            endif
         enddo
         deallocate(fieldNameList)

         if (ocnDone) then
            !---------------------------------------
            ! Initialize the atm/ocean fluxes and compute the ocean albedos
            !---------------------------------------
            call ESMF_LogWrite("MED - initialize atm/ocn fluxes and compute ocean albedo", ESMF_LOGMSG_INFO, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            ! Copy the NstateImp(compocn) to FBImp(compocn)
            call med_connectors_post_ocn2med(gcomp, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            ! Initialize the atm/ocean fluxes and compute the ocean albedos
            call ESMF_LogWrite("MED - initialize atm/ocn fluxes and compute ocean albedo", ESMF_LOGMSG_INFO, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            ! Update fractions again in case any import fields have changed
            call med_fraction_init(gcomp,rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            ! Initialize ocean albedo module and compute ocean albedos
            ! This will update the relevant module arrays in med_phases_ocnalb_mod
            call med_phases_ocnalb_init(gcomp, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            ! Initialize atm/ocn fluxes module
            ! This will update the relevant module arrays in med_phases_aoflux_mod
            call med_phases_aofluxes_init(gcomp, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         endif
      end if

    endif

    !---------------------------------------
    ! Carry out data dependency for atm initialization if needed
    !---------------------------------------

    if (.not. atmDone .and. ocnDone .and. is_local%wrap%comp_present(compatm)) then

       atmDone = .true.  ! reset if an item is found that is not done

       call ESMF_StateGet(is_local%wrap%NStateImp(compatm), itemCount=fieldCount, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       allocate(fieldNameList(fieldCount))
       call ESMF_StateGet(is_local%wrap%NStateImp(compatm), itemNameList=fieldNameList, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       do n=1, fieldCount 
          call ESMF_StateGet(is_local%wrap%NStateImp(compatm), itemName=fieldNameList(n), field=field, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          atCorrectTime = NUOPC_IsAtTime(field, time, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          if (.not. atCorrectTime) then
             ! If any atm import fields are not time stamped correctly, then dependency is not satisified - must return to atm
             call ESMF_LogWrite("MED - Initialize-Data-Dependency from ATM NOT YET SATISFIED!!!", ESMF_LOGMSG_INFO, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             atmdone = .false.
             exit  ! break out of the loop when first not satisfied found
          endif
       enddo
       deallocate(fieldNameList)

       if (.not. atmdone) then  ! atmdone is not true

          ! Update fractions again in case any import fields have changed
          call med_fraction_init(gcomp,rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          ! initialize fractions
          call med_fraction_set(gcomp, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          ! do the merge to the atmospheric component
          call med_phases_prep_atm(gcomp, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          ! copy the FBExp(compatm) to NstatExp(compatm)
          call med_connectors_prep_med2atm(gcomp, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          ! change 'Updated' attribute to true for ALL exportState fields
          call ESMF_StateGet(is_local%wrap%NStateExp(compatm), itemCount=fieldCount, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(fieldNameList(fieldCount))
          call ESMF_StateGet(is_local%wrap%NStateExp(compatm), itemNameList=fieldNameList, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          do n=1, fieldCount
             call ESMF_StateGet(is_local%wrap%NStateExp(compatm), itemName=fieldNameList(n), field=field, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end do
          deallocate(fieldNameList)

          ! Connectors will be automatically called between the mediator and atm until allDone is true
          call ESMF_LogWrite("MED - Initialize-Data-Dependency Sending Data to ATM", ESMF_LOGMSG_INFO, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
    end if

    if (atmDone .and. ocnDone) then
       if (is_local%wrap%comp_present(compatm)) then
          ! Copy the NstateImp(compatm) to FBImp(compatm)
          call med_connectors_post_atm2med(gcomp, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    endif

    allDone = .true.
    do n1 = 1,ncomps
       if (is_local%wrap%comp_present(n1) .and. ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then

          call ESMF_StateGet(is_local%wrap%NStateImp(n1), itemCount=fieldCount, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
 
          allocate(fieldNameList(fieldCount))
          call ESMF_StateGet(is_local%wrap%NStateImp(n1), itemNameList=fieldNameList, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          do n=1, fieldCount
             call ESMF_StateGet(is_local%wrap%NStateImp(n1), itemName=fieldNameList(n), field=field, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             atCorrectTime = NUOPC_IsAtTime(field, time, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             if (.not. atCorrectTime) then
                allDone=.false.
             endif
          enddo
          deallocate(fieldNameList)
       endif
    enddo

    ! set InitializeDataComplete Component Attribute to "true", indicating
    ! to the driver that this Component has fully initialized its data

    if (allDone) then
       call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_LogWrite("MED - Initialize-Data-Dependency allDone check Passed", ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call med_io_cpl_init()
    else
       call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="false", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_LogWrite("MED - Initialize-Data-Dependency allDone check Failed, another loop is required", ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_LogWrite("MED - Initialize-Data-Dependency from OCN is SATISFIED!!!", ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine DataInitialize

  !-----------------------------------------------------------------------------

  subroutine SetRunClock(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)           :: mediatorClock, driverClock
    type(ESMF_Time)            :: currTime
    type(ESMF_TimeInterval)    :: timeStep
    character(len=*),parameter :: subname='(module_MED:SetRunClock)'
    !-----------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    ! query the Mediator for clocks
    call NUOPC_MediatorGet(gcomp, mediatorClock=mediatorClock, &
      driverClock=driverClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call shr_nuopc_methods_Clock_TimePrint(driverClock  ,trim(subname)//'driver clock1',rc)
       call shr_nuopc_methods_Clock_TimePrint(mediatorClock,trim(subname)//'mediat clock1',rc)
    endif

    ! set the mediatorClock to have the current start time as the driverClock
    call ESMF_ClockGet(driverClock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockSet(mediatorClock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call shr_nuopc_methods_Clock_TimePrint(driverClock  ,trim(subname)//'driver clock2',rc)
       call shr_nuopc_methods_Clock_TimePrint(mediatorClock,trim(subname)//'mediat clock2',rc)
    endif

    ! check and set the component clock against the driver clock
    call NUOPC_CompCheckSetClock(gcomp, driverClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine SetRunClock

  !-----------------------------------------------------------------------------
#if (1 == 0)

  subroutine Mediator_restart(gcomp,mode,bfname,rc)
    !
    ! read/write mediator restart file
    !
    type(ESMF_GridComp)  :: gcomp
    character(len=*), intent(in)    :: mode
    character(len=*), intent(in)    :: bfname
    integer         , intent(inout) :: rc

    type(InternalState)  :: is_local
    character(len=1280)  :: fname
    integer              :: funit
    logical              :: fexists
    character(len=*),parameter :: subname='(module_MED:Mediator_restart)'
    !-----------------------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    if (mode /= 'write' .and. mode /= 'read') then
       call ESMF_LogWrite(trim(subname)//": ERROR mode not allowed "//trim(mode), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    fname = trim(bfname)//'_FBaccum(compatm)_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBaccum(compatm),read_rest_FBaccum(compatm),rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    fname = trim(bfname)//'_FBaccum(compocn)_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBaccum(compocn),read_rest_FBaccum(compocn),rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    fname = trim(bfname)//'_FBaccum(compice)_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBaccum(compice),read_rest_FBaccum(compice),rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    fname = trim(bfname)//'_FBaccum(complnd)_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBaccum(complnd),read_rest_FBaccum(complnd),rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    fname = trim(bfname)//'_FBaccum(comprof)_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBaccum(comprof),read_rest_FBaccum(comprof),rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    fname = trim(bfname)//'_FBaccum(compwav)_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBaccum(compwav),read_rest_FBaccum(compwav),rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    fname = trim(bfname)//'_FBaccum(compglc)_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBaccum(compglc),read_rest_FBaccum(compglc),rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    fname = trim(bfname)//'_FBaccumAOflux_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBaccumAOflux,read_rest_FBaccumAOflux,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    fname = trim(bfname)//'_FBAtm_a_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBImp(compatm,compatm),read_rest_FBAtm_a,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mode == 'read') then
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateImp(compatm), is_local%wrap%FBImp(compatm,compatm), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    fname = trim(bfname)//'_FBImp(compice,compice)_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBImp(compice,compice),read_rest_FBImp(compice,compice),rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mode == 'read') then
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateImp(compice), is_local%wrap%FBImp(compice,compice), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    fname = trim(bfname)//'_FBImp(compocn,compocn)_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBImp(compocn,compocn),read_rest_FBImp(compocn,compocn),rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mode == 'read') then
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateImp(compocn), is_local%wrap%FBImp(compocn,compocn), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    fname = trim(bfname)//'_FBImp(complnd,complnd)_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBImp(complnd,complnd),read_rest_FBImp(complnd,complnd),rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mode == 'read') then
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateImp(complnd), is_local%wrap%FBImp(complnd,complnd), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    fname = trim(bfname)//'_FBImp(comprof,comprof)_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBImp(comprof,comprof),read_rest_FBImp(comprof,comprof),rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mode == 'read') then
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateImp(comprof), is_local%wrap%FBImp(comprof,comprof), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    fname = trim(bfname)//'_FBImp(compwav,comprof)_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBImp(compwav,comprof),read_rest_FBImp(compwav,comprof),rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mode == 'read') then
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateImp(compwav), is_local%wrap%FBImp(compwav,comprof), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    fname = trim(bfname)//'_FBImp(compglc,comprof)_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBImp(compglc,comprof),read_rest_FBImp(compglc,comprof),rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mode == 'read') then
      call shr_nuopc_methods_FB_copy(is_local%wrap%NStateImp(compglc), is_local%wrap%FBImp(compglc,comprof), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    fname = trim(bfname)//'_FBAOFlux_o_restart.nc'
    call FieldBundle_RWFields(mode,fname,is_local%wrap%FBaoflux_o,read_rest_FBaoflux_o,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    funit = 1101
    fname = trim(bfname)//'_scalars_restart.txt'
    if (mode == 'write') then
      call ESMF_LogWrite(trim(subname)//": write "//trim(fname), ESMF_LOGMSG_INFO, rc=dbrc)
      open(funit,file=fname,form='formatted')
      write(funit,*) is_local%wrap%FBaccumcnt(compatm)
      write(funit,*) is_local%wrap%FBaccumcnt(compocn)
      write(funit,*) is_local%wrap%FBaccumcnt(compice)
      write(funit,*) is_local%wrap%FBaccumcntAOflux
      write(funit,*) is_local%wrap%FBaccumcnt(complnd)
      write(funit,*) is_local%wrap%FBaccumcnt(comprof)
      write(funit,*) is_local%wrap%FBaccumcnt(compwav)
      write(funit,*) is_local%wrap%FBaccumcnt(compglc)
      close(funit)
    elseif (mode == 'read') then
      inquire(file=fname,exist=fexists)
      if (fexists) then
        call ESMF_LogWrite(trim(subname)//": read "//trim(fname), ESMF_LOGMSG_INFO, rc=dbrc)
        open(funit,file=fname,form='formatted')
        ! DCR - temporary skip reading Lnd and Rof until components are added to test case
        !       restart files
        is_local%wrap%FBaccumcnt(compatm)=0
        is_local%wrap%FBaccumcnt(compocn)=0
        is_local%wrap%FBaccumcnt(compice)=0
        is_local%wrap%FBaccumcntAOflux=0
        is_local%wrap%FBaccumcnt(complnd)=0
        is_local%wrap%FBaccumcnt(comprof)=0
        is_local%wrap%FBaccumcnt(compwav)=0
        is_local%wrap%FBaccumcnt(compglc)=0
        read (funit,*) is_local%wrap%FBaccumcnt(compatm)
        read (funit,*) is_local%wrap%FBaccumcnt(compocn)
        read (funit,*) is_local%wrap%FBaccumcnt(compice)
        read (funit,*) is_local%wrap%FBaccumcntAOflux
        read (funit,*) is_local%wrap%FBaccumcnt(complnd)
        read (funit,*) is_local%wrap%FBaccumcnt(comprof)
        read (funit,*) is_local%wrap%FBaccumcnt(compwav)
        read (funit,*) is_local%wrap%FBaccumcnt(compglc)
        close(funit)
      else
        read_rest_FBaccum(compatm) = .false.
        read_rest_FBaccum(compocn) = .false.
        read_rest_FBaccum(compice) = .false.
        read_rest_FBaccum(complnd) = .false.
        read_rest_FBaccum(comprof) = .false.
        read_rest_FBaccum(compwav) = .false.
        read_rest_FBaccum(compglc) = .false.
        read_rest_FBaccumAOflux    = .false.
      endif
    endif

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine Mediator_restart
#endif

  !-----------------------------------------------------------------------------

end module MED
