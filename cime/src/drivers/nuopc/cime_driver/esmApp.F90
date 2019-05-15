program esmApp

  !-----------------------------------------------------------------------------
  ! Generic ESM application driver
  !-----------------------------------------------------------------------------

  use ESMF,            only : ESMF_Initialize, ESMF_CALKIND_GREGORIAN, ESMF_LOGKIND_MULTI
  use ESMF,            only : ESMF_END_ABORT, ESMF_LogFoundError, ESMF_Finalize, ESMF_LOGERR_PASSTHRU
  use ESMF,            only : ESMF_GridCompSetServices, ESMF_GridCompFinalize, ESMF_LogSet, ESMF_LogWrite
  use ESMF,            only : ESMF_GridCompDestroy, ESMF_LOGMSG_INFO, ESMF_GridComp, ESMF_GridCompRun
  use ESMF,            only : ESMF_GridCompFinalize, ESMF_GridCompCreate, ESMF_GridCompInitialize
  use ESMF,            only : ESMF_LOGKIND_MULTI_ON_ERROR
  use mpi,             only : MPI_COMM_WORLD, MPI_COMM_NULL, MPI_Init, MPI_FINALIZE
  use NUOPC,           only : NUOPC_FieldDictionarySetup
  use ensemble_driver, only : SetServices
  use shr_pio_mod,     only : shr_pio_init1, shr_pio_init2

  implicit none

  ! local variables
  integer                    :: COMP_COMM
  integer                    :: rc, urc
  type(ESMF_GridComp)        :: ensemble_driver_comp

  !-----------------------------------------------------------------------------
  ! Initiallize MPI
  !-----------------------------------------------------------------------------

  call MPI_init(rc)
  COMP_COMM = MPI_COMM_WORLD

  !-----------------------------------------------------------------------------
  ! Initialize PIO
  !-----------------------------------------------------------------------------

  ! For planned future use of async io using pio2.  The IO tasks are seperated from the compute tasks here
  ! and COMP_COMM will be MPI_COMM_NULL on the IO tasks which then call shr_pio_init2 and do not return until
  ! the model completes.  All other tasks call ESMF_Initialize.  8 is the maximum number of component models
  ! supported

  call shr_pio_init1(8, "drv_in", COMP_COMM)

  if (COMP_COMM .eq. MPI_COMM_NULL) then
     ! call shr_pio_init2(
     call mpi_finalize(ierror=rc)
     stop
  endif

  !-----------------------------------------------------------------------------
  ! Initialize ESMF
  !-----------------------------------------------------------------------------

#ifdef DEBUG
  call ESMF_Initialize(mpiCommunicator=COMP_COMM, logkindflag=ESMF_LOGKIND_MULTI, logappendflag=.false., &
       defaultCalkind=ESMF_CALKIND_GREGORIAN, ioUnitLBound=5001, ioUnitUBound=5101, rc=rc)
#else
  call ESMF_Initialize(mpiCommunicator=COMP_COMM, logkindflag=ESMF_LOGKIND_MULTI_ON_ERROR, logappendflag=.false., &
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
  ! call ESMF_LogSet(flush=.true., trace=.true., rc=rc)
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
