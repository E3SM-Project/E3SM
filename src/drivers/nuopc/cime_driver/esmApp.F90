program esmApp

  !-----------------------------------------------------------------------------
  ! Generic ESM application driver
  !-----------------------------------------------------------------------------

  use ESMF, only             : ESMF_Initialize, ESMF_CALKIND_GREGORIAN, ESMF_LOGKIND_MULTI
  use ESMF, only             : ESMF_END_ABORT, ESMF_LogFoundError, ESMF_Finalize, ESMF_LOGERR_PASSTHRU
  use ESMF, only             : ESMF_GridCompSetServices, ESMF_GridCompFinalize, ESMF_LogSet, ESMF_LogWrite
  use ESMF, only             : ESMF_GridCompDestroy, ESMF_LOGMSG_INFO, ESMF_GridComp, ESMF_GridCompRun
  use ESMF, only             : ESMF_GridCompFinalize, ESMF_GridCompCreate, ESMF_GridCompInitialize
  use ESMF, only             : ESMF_LOGKIND_MULTI_ON_ERROR
  use NUOPC, only            : NUOPC_FieldDictionarySetup
  use ensemble_driver,  only : SetServices

  implicit none

  integer                    :: rc, urc
  type(ESMF_GridComp)        :: ensemble_driver_comp

  !-----------------------------------------------------------------------------
  ! Initialize ESMF
  !-----------------------------------------------------------------------------
#ifdef DEBUG
  call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, logappendflag=.false., &
    defaultCalkind=ESMF_CALKIND_GREGORIAN, ioUnitLBound=5001, ioUnitUBound=5101, rc=rc)
#else
  call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI_ON_ERROR, logappendflag=.false., &
    defaultCalkind=ESMF_CALKIND_GREGORIAN, ioUnitLBound=5001, ioUnitUBound=5101, rc=rc)
#endif
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_LogSet(flush=.true.)

  call ESMF_LogWrite("esmApp STARTING", ESMF_LOGMSG_INFO, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  ! Operate on the NUOPC Field dictionary
  !-----------------------------------------------------------------------------

  call NUOPC_FieldDictionarySetup("fd.yaml", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  ! Create the earth system ensemble driver Component
  !-----------------------------------------------------------------------------

  ensemble_driver_comp = ESMF_GridCompCreate(name="ensemble", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  ! SetServices for the ensemble driver Component
  !-----------------------------------------------------------------------------
  call ESMF_GridCompSetServices(ensemble_driver_comp, SetServices, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  ! Call Initialize for the earth system ensemble Component
  !-----------------------------------------------------------------------------

  call ESMF_GridCompInitialize(ensemble_driver_comp, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  ! Call Run  for the ensemble driver
  !-----------------------------------------------------------------------------
  call ESMF_GridCompRun(ensemble_driver_comp, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  ! Call Finalize for the ensemble driver
  ! Destroy the ensemble driver
  !-----------------------------------------------------------------------------
  call ESMF_GridCompFinalize(ensemble_driver_comp, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_LogWrite("ESMF_GridCompDestroy called", ESMF_LOGMSG_INFO, rc=rc)
!  call ESMF_LogSet(flush=.true., trace=.true., rc=rc)
  call ESMF_GridCompDestroy(ensemble_driver_comp, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_LogWrite("ESMF_GridCompDestroy finished", ESMF_LOGMSG_INFO, rc=rc)

  !-----------------------------------------------------------------------------
  ! Finalize ESMF
  !-----------------------------------------------------------------------------

  call ESMF_LogWrite("esmApp FINISHED", ESMF_LOGMSG_INFO, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
  ! Finalize ESMF
  call ESMF_Finalize()

end program
