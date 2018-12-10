module ESM

  !-----------------------------------------------------------------------------
  ! Code that specializes generic ESM Component code.
  !-----------------------------------------------------------------------------

  use ESMF                  , only : ESMF_Clock
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use shr_kind_mod          , only : SHR_KIND_R8, SHR_KIND_CS, SHR_KIND_CL
  use shr_log_mod           , only : shr_log_Unit, shr_log_Level
  use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
  use med_internalstate_mod , only : logunit, loglevel
  use seq_comm_mct, only : num_inst_total
  use shr_nuopc_utils_mod, only : shr_nuopc_memcheck

  implicit none
  private

  character(len=512)             :: msgstr
  logical                        :: mastertask ! master processor for driver gcomp
  integer                        :: componentCount
  character(len=8)               :: atm_present, lnd_present, ocn_present
  character(len=8)               :: ice_present, rof_present, wav_present
  character(len=8)               :: glc_present, med_present
  character(*), parameter        :: nlfilename = "drv_in" ! input namelist filename
  character(*), parameter        :: u_FILE_u = &
       __FILE__

  public  :: SetServices
  ! used in ensemble_driver
  public :: ReadAttributes

  private :: AddAttributes
  private :: SetModelServices
  private :: SetRunSequence
  private :: ModifyCplLists
  private :: InitAttributes
  private :: CheckAttributes
  private :: InitAdvertize

!================================================================================
  contains
!================================================================================

  subroutine SetServices(driver, rc)

    use NUOPC        , only : NUOPC_CompDerive, NUOPC_CompSpecialize, NUOPC_CompSetInternalEntryPoint
    use NUOPC_Driver , only : driver_routine_SS             => SetServices
    use NUOPC_Driver , only : driver_label_SetModelServices => label_SetModelServices
    use NUOPC_Driver , only : driver_label_SetRunSequence   => label_SetRunSequence
    use NUOPC_Driver , only : driver_label_Finalize         => label_Finalize
    use ESMF         , only : ESMF_GridComp, ESMF_Config, ESMF_GridCompSet, ESMF_ConfigLoadFile
    use ESMF         , only : ESMF_METHOD_INITIALIZE
    use ESMF         , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO

    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    integer           :: dbrc
    type(ESMF_Config)    :: runSeq
    character(len=*), parameter :: subname = "(esm.F90:SetServices)"
    !---------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    ! NUOPC_Driver registers the generic methods
    call NUOPC_CompDerive(driver, driver_routine_SS, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(driver, specLabel=driver_label_SetModelServices, &
         specRoutine=SetModelServices, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(driver, specLabel=driver_label_SetRunSequence, &
         specRoutine=SetRunSequence, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! register an internal initialization method
    call NUOPC_CompSetInternalEntryPoint(driver, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p2"/), userRoutine=ModifyCplLists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !
    ! This prevents the driver trying to "auto" connect to the ensemble_driver
    ! by default the FieldTransferPolicy is "transferall" and we need "transfernone"
    !
    call NUOPC_CompSetInternalEntryPoint(driver, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv05p1"/), userRoutine=InitAdvertize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set a finalize method
    call NUOPC_CompSpecialize(driver, specLabel=driver_label_Finalize, &
         specRoutine=esm_finalize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create, open and set the config

    call ESMF_GridCompSet(driver, configFile="nuopc.runconfig", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine SetServices

  !================================================================================

  subroutine SetModelServices(driver, rc)

    use ESMF                  , only : ESMF_GridComp, ESMF_VM, ESMF_Config, ESMF_VMBarrier
    use ESMF                  , only : ESMF_GridCompGet, ESMF_VMGet, ESMF_ConfigGetAttribute
    use ESMF                  , only : ESMF_ConfigGetLen, ESMF_RC_NOT_VALID, ESMF_LogFoundAllocError
    use ESMF                  , only : ESMF_LogSetError, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use ESMF                  , only : ESMF_GridCompSet, ESMF_SUCCESS, ESMF_METHOD_INITIALIZE
    use ESMF                  , only : ESMF_VMisCreated, ESMF_GridCompIsPetLocal
    use ESMF                  , only : ESMF_RC_FILE_OPEN, ESMF_RC_FILE_READ
    use ESMF                  , only : ESMF_AttributeUpdate, ESMF_VMBroadcast
    use ESMF                  , only : ESMF_MethodAdd
    use NUOPC                 , only : NUOPC_CompSetInternalEntryPoint, NUOPC_CompAttributeGet
    use NUOPC                 , only : NUOPC_CompAttributeAdd, NUOPC_CompAttributeSet
    use NUOPC_Driver          , only : NUOPC_DriverAddComp, NUOPC_DriverGetComp

    use shr_nuopc_methods_mod , only : shr_nuopc_methods_Clock_TimePrint
    use shr_file_mod          , only : shr_file_setLogunit, shr_file_getunit
    use med                   , only : med_SS         => SetServices
    use atm_comp_nuopc        , only : ATMSetServices => SetServices
    use ice_comp_nuopc        , only : ICESetServices => SetServices
    use lnd_comp_nuopc        , only : LNDSetServices => SetServices
    use ocn_comp_nuopc        , only : OCNSetServices => SetServices
    use wav_comp_nuopc        , only : WAVSetServices => SetServices
    use rof_comp_nuopc        , only : ROFSetServices => SetServices
    use glc_comp_nuopc        , only : GLCSetServices => SetServices
    use shr_pio_mod           , only : shr_pio_init1
    use pio                   , only : pio_file_is_open, pio_closefile, file_desc_t
    use perf_mod              , only : t_initf
    use shr_nuopc_time_mod    , only : shr_nuopc_time_clockInit
    use shr_log_mod           , only : shrlogunit=> shr_log_unit

    type(ESMF_GridComp)    :: driver
    integer, intent(out)   :: rc


    ! local variables
    type(ESMF_GridComp)    :: mediator
    type(ESMF_Clock)       :: clock
    type(ESMF_VM)          :: vm
    type(ESMF_GridComp)    :: child
    type(ESMF_Config)      :: config
    integer                :: compid
    type(file_desc_t)      :: pioid
    integer                :: n, n1, stat
    integer, pointer       :: petList(:)
    character(len=20)      :: model, prefix
    integer                :: petCount, i
    integer                :: localPet, medpet
    logical                :: is_set
    character(SHR_KIND_CS) :: cvalue
    character(len=512)     :: diro
    character(len=512)     :: logfile
    integer                :: global_comm
    logical                :: isPresent
    integer                :: maxthreads
    integer                :: dbrc
    integer                :: readunit, iostat
    logical                :: is_restart
    character(len=512)     :: restartfile
    character(len=*), parameter    :: subname = "(esm.F90:SetModelServices)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    !-------------------------------------------
    ! Set the io logunit to the value defined in ensemble_driver
    ! it may be corrected below if the med mastertask is not the driver mastertask
    !-------------------------------------------
    call shr_file_setLogunit(logunit)

    !-------------------------------------------
    ! Get the config and vm objects from the driver
    !-------------------------------------------

    call ESMF_GridCompGet(driver, vm=vm, config=config, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=localPet, mpiCommunicator=global_comm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (localPet == 0) then
       mastertask = .true.
    else
       mastertask = .false.
    end if

    !-------------------------------------------
    ! determine the generic component labels
    !-------------------------------------------

    componentCount = ESMF_ConfigGetLen(config,label="CESM_component_list:", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (componentCount == 0) then
      write (msgstr, *) "No models were specified in CESM_component_list "
      call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
      return  ! bail out
    endif

    !-------------------------------------------
    ! Obtain driver attributes
    !-------------------------------------------

    call ReadAttributes(driver, config, "DRIVER_attributes::", formatprint=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(driver, config, "FLDS_attributes::", formatprint=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(driver, config, "CLOCK_attributes::", formatprint=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(driver, config, "ALLCOMP_attributes::", formatprint=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(driver, config, "PELAYOUT_attributes::", formatprint=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call CheckAttributes(driver, rc)

    !-------------------------------------------
    ! Initialize communicators and PIO
    !-------------------------------------------

    ! Call first phase of pio initialization The call to pio_init1
    ! should be the first routine called after mpi_init.  It reads the
    ! pio default settings (pio_default_inparm) from the file
    ! 'nlfilename' and if the namelist variable 'pio_async_interface'
    ! is true, it splits the IO tasks away from the Compute tasks.  It
    ! then returns the new compute comm in Global_Comm and sets module
    ! variable io_comm.
    ! TODO: this must be reconciled with having the asynchronous io
    ! processors just be a separate gridded component in NUOPC
    ! TODO: global_comm should be the same as mpicom for the driver vm - however
    ! cannot set this to mpicom since the call to shr_pio_init1 can possibly change this

    call shr_pio_init1(num_inst_total, nlfilename, global_comm)

    !-------------------------------------------
    ! Initialize other attributes (after initializing driver clock)
    !-------------------------------------------

    call InitAttributes(driver, rc)

    !-------------------------------------------
    ! Initialize component pe layouts
    !-------------------------------------------

    call esm_init_pelayout(driver, maxthreads)

    !-------------------------------------------
    ! Reset log unit for mediator
    !-------------------------------------------
    call NUOPC_DriverGetComp(driver, 'MED', comp=mediator, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    medPet = -1
    if (ESMF_GridCompIsPetLocal(mediator, rc=rc)) then
       call ESMF_GridCompGet(mediator, vm=vm, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMGet(vm, localPet=medPet, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       ! ensemble_driver set the med pet to task 0 of the member driver, correct it here
       if(medPet == 0 .and. localPet /= 0) then
          call NUOPC_CompAttributeGet(driver, name="diro", value=diro, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeGet(driver, name="logfile", value=logfile, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          logunit = shr_file_getUnit()
          open(logunit,file=trim(diro)//"/"//trim(logfile), position='append')
          mastertask = .true.
       else
          logUnit = shrlogunit
          mastertask = .false.
       endif
    else
       logUnit = shrlogunit
       mastertask = .false.
    endif
          ! Print out present flags to mediator log file
    if (medPet==0) then
       write(logunit,*) trim(subname)//":atm_present="//trim(atm_present)
       write(logunit,*) trim(subname)//":lnd_present="//trim(lnd_present)
       write(logunit,*) trim(subname)//":ocn_present="//trim(ocn_present)
       write(logunit,*) trim(subname)//":ice_present="//trim(ice_present)
       write(logunit,*) trim(subname)//":rof_present="//trim(rof_present)
       write(logunit,*) trim(subname)//":wav_present="//trim(wav_present)
       write(logunit,*) trim(subname)//":glc_present="//trim(glc_present)
       write(logunit,*) trim(subname)//":med_present="//trim(med_present)
    end if

    !-------------------------------------------
    ! Timer initialization (has to be after pelayouts are determined)
    !-------------------------------------------

    call t_initf(nlfilename, LogPrint=.true., mpicom=global_comm, &
         mastertask=mastertask, MaxThreads=maxthreads)

    ! finish the pio initialization (this calls shr_pio_init2)
    call InitPIO(driver, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------
    ! Perform restarts if appropriate
    !-------------------------------------------

    call InitRestart(driver, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (pio_file_is_open(pioid)) then
       call pio_closefile(pioid)
    end if

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine SetModelServices

  !================================================================================

  subroutine SetRunSequence(driver, rc)
    use ESMF                  , only : ESMF_GridComp, ESMF_LogWrite, ESMF_SUCCESS, ESMF_LOGMSG_INFO
    use ESMF                  , only : ESMF_Time, ESMF_TimeInterval, ESMF_Clock, ESMF_Config
    use ESMF                  , only : ESMF_GridCompGet, ESMF_ConfigLoadFile, ESMF_ConfigCreate
    use NUOPC                 , only : NUOPC_FreeFormat, NUOPC_FreeFormatPrint, NUOPC_FreeFormatDestroy
    use NUOPC                 , only : NUOPC_FreeFormatCreate
    use NUOPC_Driver          , only : NUOPC_DriverIngestRunSequence, NUOPC_DriverSetRunSequence
    use NUOPC_Driver          , only : NUOPC_DriverPrint

    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    integer                 :: localrc
    type(ESMF_Config)       :: runSeq
    type(NUOPC_FreeFormat)  :: runSeqFF
    integer                 :: dbrc
    character(len=*), parameter :: subname = "(esm.F90:SetRunSequence)"

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    !--------
    ! Run Sequence and Connectors
    !--------

    ! read free format run sequence

    runSeq = ESMF_ConfigCreate(rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ConfigLoadFile(runSeq, "nuopc.runseq", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    runSeqFF = NUOPC_FreeFormatCreate(runSeq, label="runSeq::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_DriverIngestRunSequence(driver, runSeqFF, autoAddConnectors=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Uncomment these to add debugging information for driver
    ! call NUOPC_DriverPrint(driver, orderflag=.true.)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !   line=__LINE__, &
    !   file=__FILE__)) &
    !   return  ! bail out

    if (mastertask) then
       call NUOPC_FreeFormatPrint(runSeqFF, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call NUOPC_FreeFormatDestroy(runSeqFF, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine SetRunSequence

  !================================================================================

  recursive subroutine ModifyCplLists(driver, importState, exportState, clock, rc)

    use ESMF         , only : ESMF_GridComp, ESMF_State, ESMF_Clock, ESMF_LogWrite
    use ESMF         , only : ESMF_LOGMSG_INFO, ESMF_CplComp, ESMF_SUCCESS
    use NUOPC        , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
    use NUOPC_Driver , only : NUOPC_DriverGetComp

    type(ESMF_GridComp)  :: driver
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(ESMF_CplComp), pointer     :: connectorList(:)
    integer                         :: i, j, cplListSize
    character(len=160), allocatable :: cplList(:)
    character(len=160)              :: tempString
    integer                         :: dbrc
    character(len=*), parameter :: subname = "(esm.F90:ModifyCplLists)"

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_LogWrite("Driver is in ModifyCplLists()", ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    nullify(connectorList)
    call NUOPC_DriverGetComp(driver, compList=connectorList, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    write (msgstr,*) "Found ", size(connectorList), " Connectors."// " Modifying CplList Attribute...."
    call ESMF_LogWrite(trim(msgstr), ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    do i=1, size(connectorList)

      ! query the cplList for connector i
      call NUOPC_CompAttributeGet(connectorList(i), name="CplList", itemCount=cplListSize, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      if (cplListSize>0) then
        allocate(cplList(cplListSize))

        call NUOPC_CompAttributeGet(connectorList(i), name="CplList", valueList=cplList, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        ! go through all of the entries in the cplList and set the mapping method to "redist"
        do j=1, cplListSize
           !tempString = trim(cplList(j))//":REMAPMETHOD=bilinear"//&
           !":SrcTermProcessing=1:DUMPWEIGHTS=true:TermOrder=SrcSeq"

           tempString = trim(cplList(j))//":remapmethod=redist"
           cplList(j) = trim(tempString)
        enddo

        ! store the modified cplList in CplList attribute of connector i
        call NUOPC_CompAttributeSet(connectorList(i), name="CplList", valueList=cplList, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        deallocate(cplList)
      endif
    enddo

    deallocate(connectorList)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine ModifyCplLists

  !================================================================================

  subroutine InitPIO(driver, rc)

    !----------------------------------------------------------
    ! Initialize MPI communicators,  IO and threading
    !----------------------------------------------------------
    use ESMF         , only : ESMF_GridComp, ESMF_LOGMSG_INFO, ESMF_LogWrite, ESMF_SUCCESS
    use NUOPC        , only : NUOPC_CompAttributeGet
    use shr_pio_mod  , only : shr_pio_init2
    use mpi          , only : MPI_COMM_NULL
    use seq_comm_mct , only : CPLID, GLOID, ATMID, LNDID, OCNID, ICEID, GLCID, ROFID, WAVID, ESPID
    use seq_comm_mct , only : num_inst_atm, num_inst_lnd, num_inst_rof
    use seq_comm_mct , only : num_inst_ocn, num_inst_ice, num_inst_glc
    use seq_comm_mct , only : num_inst_wav, num_inst_esp
    use shr_mem_mod  , only : shr_mem_init
    use seq_comm_mct , only : seq_comm_inst, seq_comm_name, seq_comm_suffix
    use seq_comm_mct , only : seq_comm_setnthreads, seq_comm_getnthreads
    use seq_comm_mct , only : seq_comm_iamin, seq_comm_name, seq_comm_namelen, seq_comm_iamroot
    use seq_comm_mct , only : seq_comm_getinfo => seq_comm_setptrs

    ! input/output variables
    type(ESMF_GridComp), intent(inout) :: driver
    integer            , intent(out)   :: rc

    ! local variables
    integer                         :: it,n
    integer                         :: mpicom_GLOID          ! MPI global communicator
    integer                         :: mpicom_OCNID          ! MPI ocn communicator for ensemble member 1
    integer                         :: nthreads_GLOID        ! OMP global number of threads
    integer                         :: nthreads_CPLID        ! OMP cpl number of threads
    integer                         :: nthreads_ATMID        ! OMP atm number of threads
    integer                         :: nthreads_LNDID        ! OMP lnd number of threads
    integer                         :: nthreads_ICEID        ! OMP ice number of threads
    integer                         :: nthreads_OCNID        ! OMP ocn number of threads
    integer                         :: nthreads_GLCID        ! OMP glc number of threads
    integer                         :: nthreads_ROFID        ! OMP glc number of threads
    integer                         :: nthreads_WAVID        ! OMP wav number of threads
    integer                         :: nthreads_ESPID        ! OMP esp number of threads
    integer                         :: pethreads_GLOID       ! OMP number of threads per task
    logical                         :: drv_threading         ! driver threading control
    integer                         :: maxthreads
    integer                         :: comp_id(num_inst_total)
    integer                         :: comp_comm(num_inst_total)
    integer                         :: comp_comm_iam(num_inst_total)
    logical                         :: comp_iamin(num_inst_total)
    logical                         :: iamroot_med
    character(SHR_KIND_CL)          :: cvalue
    character(len=seq_comm_namelen) :: comp_name(num_inst_total)
    integer :: dbrc
    character(len=*) , parameter    :: subname = "(esm.F90:InitPIO)"
    !----------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    ! set task based threading counts
    call seq_comm_getinfo(GLOID, pethreads=pethreads_GLOID)
    call seq_comm_setnthreads(pethreads_GLOID)
    call seq_comm_getinfo(GLOID, mpicom=mpicom_GLOID, nthreads=nthreads_GLOID)

    ! determine arrays comp_id, comp_comm, comp_iam and comp_name - these are needed for call
    ! to shr_pio_init2

    comp_comm(:) = MPI_COMM_NULL

    it=1
    comp_id(it)    = CPLID
    comp_iamin(it) = seq_comm_iamin(comp_id(it))
    comp_name(it)  = seq_comm_name(comp_id(it))
    call seq_comm_getinfo(CPLID, mpicom=comp_comm(it), nthreads=nthreads_CPLID, iam=comp_comm_iam(it))

    do n = 1,num_inst_atm
       it=it+1
       comp_id(it)    = ATMID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(ATMID(n), mpicom=comp_comm(it), nthreads=nthreads_ATMID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_lnd
       it=it+1
       comp_id(it)    = LNDID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(LNDID(n), mpicom=comp_comm(it), nthreads=nthreads_LNDID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_ocn
       it=it+1
       comp_id(it)    = OCNID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(OCNID(n), mpicom=comp_comm(it), nthreads=nthreads_OCNID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_ice
       it=it+1
       comp_id(it)    = ICEID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(ICEID(n), mpicom=comp_comm(it), nthreads=nthreads_ICEID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_glc
       it=it+1
       comp_id(it)    = GLCID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(GLCID(n), mpicom=comp_comm(it), nthreads=nthreads_GLCID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_rof
       it=it+1
       comp_id(it)    = ROFID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(ROFID(n), mpicom=comp_comm(it), nthreads=nthreads_ROFID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_wav
       it=it+1
       comp_id(it)    = WAVID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(WAVID(n), mpicom=comp_comm(it), nthreads=nthreads_WAVID, iam=comp_comm_iam(it))
    enddo
    do n = 1,num_inst_esp
       it=it+1
       comp_id(it)    = ESPID(n)
       comp_iamin(it) = seq_comm_iamin(comp_id(it))
       comp_name(it)  = seq_comm_name(comp_id(it))
       call seq_comm_getinfo(ESPID(n), mpicom=comp_comm(it), nthreads=nthreads_ESPID, iam=comp_comm_iam(it))
    enddo
    ! ESP components do not use the coupler (they are 'external')

    !----------------------------------------------------------
    ! Log info about the environment settings
    !----------------------------------------------------------

    !  When using io servers (pio_async_interface=.true.) the server tasks do not return from
    !  shr_pio_init2
    call shr_pio_init2(comp_id, comp_name, comp_iamin, comp_comm, comp_comm_iam)

    ! Timer initialization (has to be after mpi init)

    maxthreads = max(nthreads_GLOID, nthreads_CPLID, nthreads_ATMID,                 &
                     nthreads_LNDID, nthreads_ICEID, nthreads_OCNID, nthreads_GLCID, &
                     nthreads_ROFID, nthreads_WAVID, nthreads_ESPID, pethreads_GLOID )

    !----------------------------------------------------------
    ! Initialize timing library
    !----------------------------------------------------------

    !----------------------------------------------------------
    ! Test Threading Setup in driver happens to be valid on all pes for all IDs
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="drv_threading", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) drv_threading

    if (drv_threading) then
       if (mastertask) write(logunit,*) ' '
       if (mastertask) write(logunit,'(2A)    ') subname,' Test Threading in driver'
       call seq_comm_setnthreads(nthreads_GLOID)
       if (mastertask) write(logunit,'(2A,2I4)') subname,'    nthreads_GLOID = ',nthreads_GLOID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_CPLID)
       if (mastertask) write(logunit,'(2A,2I4)') subname,'    nthreads_CPLID = ',nthreads_CPLID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_ATMID)
       if (mastertask) write(logunit,'(2A,2I4)') subname,'    nthreads_ATMID = ',nthreads_ATMID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_LNDID)
       if (mastertask) write(logunit,'(2A,2I4)') subname,'    nthreads_LNDID = ',nthreads_LNDID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_OCNID)
       if (mastertask) write(logunit,'(2A,2I4)') subname,'    nthreads_OCNID = ',nthreads_OCNID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_ICEID)
       if (mastertask) write(logunit,'(2A,2I4)') subname,'    nthreads_ICEID = ',nthreads_ICEID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_GLCID)
       if (mastertask) write(logunit,'(2A,2I4)') subname,'    nthreads_GLCID = ',nthreads_GLCID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_ROFID)
       if (mastertask) write(logunit,'(2A,2I4)') subname,'    nthreads_ROFID = ',nthreads_ROFID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_WAVID)
       if (mastertask) write(logunit,'(2A,2I4)') subname,'    nthreads_WAVID = ',nthreads_WAVID,seq_comm_getnthreads()
       call seq_comm_setnthreads(nthreads_ESPID)
       if (mastertask) write(logunit,'(2A,2I4)') subname,'    nthreads_ESPID = ',nthreads_ESPID,seq_comm_getnthreads()
       if (mastertask) write(logunit,*) ' '
       call seq_comm_setnthreads(nthreads_GLOID)
    endif

    !----------------------------------------------------------
    ! Memory test
    !----------------------------------------------------------

    call seq_comm_getinfo(CPLID, iamroot=iamroot_med)
    call shr_mem_init(prt=iamroot_med)

  end subroutine InitPIO

  !================================================================================

  function IsRestart(gcomp, rc)

    use ESMF         , only : ESMF_GridComp, ESMF_SUCCESS
    use ESMF         , only : ESMF_LogSetError, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_RC_NOT_VALID
    use NUOPC        , only : NUOPC_CompAttributeGet

    ! input/output variables
    logical                                :: IsRestart
    type(ESMF_GridComp)    , intent(inout) :: gcomp
    integer                , intent(out)   :: rc

    ! locals
    character(len=*) , parameter :: start_type_start = "startup"
    character(len=*) , parameter :: start_type_cont  = "continue"
    character(len=*) , parameter :: start_type_brnch = "branch"
    character(SHR_KIND_CL)       :: start_type     ! Type of startup
    character(len=*), parameter  :: subname = "(esm.F90:IsRestart)"

    rc = ESMF_SUCCESS

    ! First Determine if restart is read
    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=start_type, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if ((trim(start_type) /= start_type_start) .and.  &
        (trim(start_type) /= start_type_cont ) .and.  &
        (trim(start_type) /= start_type_brnch)) then
       write (msgstr, *) subname//': start_type invalid = '//trim(start_type)
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
    end if

    !TODO: this is hard-wired to CIME start/continue types in terms of gcomp
    IsRestart = .false.
    if (trim(start_type) == trim(start_type_cont) .or. trim(start_type) == trim(start_type_brnch)) then
       IsRestart = .true.
    end if

  end function IsRestart

  !================================================================================

  subroutine InitRestart(driver, rc)

    !-----------------------------------------------------
    ! Determine if will restart and read pointer file
    ! if appropriate
    !-----------------------------------------------------
    use ESMF         , only : ESMF_GridComp, ESMF_VM, ESMF_GridCompGet, ESMF_VMGet, ESMF_SUCCESS
    use ESMF         , only : ESMF_LogSetError, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_RC_NOT_VALID
    use NUOPC        , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
    use shr_sys_mod  , only : shr_sys_abort
    use shr_file_mod , only : shr_file_getUnit, shr_file_freeUnit
    use shr_mpi_mod  , only : shr_mpi_bcast

    ! input/output variables
    type(ESMF_GridComp)    , intent(inout) :: driver
    integer                , intent(out)   :: rc

    ! local variables
    character(SHR_KIND_CL)       :: cvalue         ! temporary
    logical                      :: read_restart   ! read the restart file, based on start_type
    character(SHR_KIND_CL)       :: rest_case_name ! Short case identification
    character(len=*) , parameter :: subname = "(esm.F90:InitRestart)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    endif

    !-----------------------------------------------------
    ! Carry out restart if appropriate
    !-----------------------------------------------------

    read_restart = IsRestart(driver, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Add rest_case_name and read_restart to driver attributes
    call NUOPC_CompAttributeAdd(driver, attrList=(/'rest_case_name','read_restart'/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    rest_case_name = ' '
    call NUOPC_CompAttributeSet(driver, name='rest_case_name', value=rest_case_name, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    write(cvalue,*) read_restart
    call NUOPC_CompAttributeSet(driver, name='read_restart', value=trim(cvalue), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitRestart

  !================================================================================

  subroutine InitAttributes(driver, rc)

    use shr_sys_mod      , only : shr_sys_abort
    use ESMF             , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF             , only : ESMF_Clock, ESMF_ClockGet, ESMF_Time, ESMF_TimeGet
    use ESMF             , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LogSetError, ESMF_LOGMSG_INFO
    use ESMF             , only : ESMF_RC_NOT_VALID
    use ESMF             , only : ESMF_GridCompIsPetLocal, ESMF_VMBroadcast
    use NUOPC            , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
    use shr_orb_mod      , only : shr_orb_params, SHR_ORB_UNDEF_INT, SHR_ORB_UNDEF_REAL
    use seq_comm_mct     , only : CPLID, OCNID
    use seq_comm_mct     , only : seq_comm_getinfo => seq_comm_setptrs
    use shr_assert_mod   , only : shr_assert_in_domain
    use shr_cal_mod      , only : shr_cal_date2ymd
    use shr_const_mod    , only : shr_const_tkfrz, shr_const_tktrip
    use shr_const_mod    , only : shr_const_mwwv, shr_const_mwdair
    use shr_frz_mod      , only : shr_frz_freezetemp_init
    use shr_reprosum_mod , only : shr_reprosum_setopts
    use shr_wv_sat_mod   , only : shr_wv_sat_set_default, shr_wv_sat_init
    use shr_wv_sat_mod   , only : shr_wv_sat_make_tables, ShrWVSatTableSpec
    use shr_wv_sat_mod   , only : shr_wv_sat_get_scheme_idx, shr_wv_sat_valid_idx
   !use shr_scam_mod     , only : shr_scam_checkSurface

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: driver
    integer             , intent(out)   :: rc                    ! return code

    ! local variables
    type(ESMF_Clock)                :: clock
    type(ESMF_Time)                 :: currTime
    character(SHR_KIND_CL)          :: errstring
    character(SHR_KIND_CL)          :: cvalue
    integer                         :: mpicom_OCNID          ! MPI ocn communicator for ensemble member 1
    logical                         :: drv_threading         ! driver threading control
    logical                         :: reprosum_use_ddpdd    ! setup reprosum, use ddpdd
    real(SHR_KIND_R8)               :: reprosum_diffmax      ! setup reprosum, set rel_diff_max
    logical                         :: reprosum_recompute    ! setup reprosum, recompute if tolerance exceeded
    integer                         :: year                  ! Current date (YYYY)
    character(SHR_KIND_CS)          :: tfreeze_option        ! Freezing point calculation
    character(SHR_KIND_CL)          :: orb_mode              ! orbital mode
    integer                         :: orb_iyear             ! orbital year
    integer                         :: orb_iyear_align       ! associated with model year
    integer                         :: orb_cyear             ! orbital year for current orbital computation
    integer                         :: orb_nyear             ! orbital year associated with currrent model year
    integer                         :: orbitmp(4)            ! array for integer parameter broadcast
    real(SHR_KIND_R8)               :: orbrtmp(6)            ! array for real parameter broadcast
    real(SHR_KIND_R8)               :: orb_eccen             ! orbital eccentricity
    real(SHR_KIND_R8)               :: orb_obliq             ! obliquity in degrees
    real(SHR_KIND_R8)               :: orb_mvelp             ! moving vernal equinox long
    real(SHR_KIND_R8)               :: orb_obliqr            ! Earths obliquity in rad
    real(SHR_KIND_R8)               :: orb_lambm0            ! Mean long of perihelion at vernal equinox (radians)
    real(SHR_KIND_R8)               :: orb_mvelpp            ! moving vernal equinox long
    real(SHR_KIND_R8)               :: wall_time_limit       ! wall time limit in hours
    logical                         :: single_column         ! scm mode logical
    real(SHR_KIND_R8)               :: scmlon                ! single column lon
    real(SHR_KIND_R8)               :: scmlat                ! single column lat
    integer                         :: maxthreads
    character(SHR_KIND_CS)          :: wv_sat_scheme
    real(SHR_KIND_R8)               :: wv_sat_transition_start
    logical                         :: wv_sat_use_tables
    real(SHR_KIND_R8)               :: wv_sat_table_spacing
    type(ShrWVSatTableSpec)         :: liquid_spec
    type(ShrWVSatTableSpec)         :: ice_spec
    type(ShrWVSatTableSpec)         :: mixed_spec
    logical                         :: flag
    integer                         :: i, it, n
    integer                         :: unitn                 ! Namelist unit number to read
    integer                         :: dbrc
    integer                         :: localPet, rootpe_med
    integer          , parameter    :: ens1=1                ! use first instance of ensemble only
    integer          , parameter    :: fix1=1                ! temporary hard-coding to first ensemble, needs to be fixed
    real(SHR_KIND_R8), parameter    :: epsilo = shr_const_mwwv/shr_const_mwdair
    character(len=*) , parameter    :: orb_fixed_year       = 'fixed_year'
    character(len=*) , parameter    :: orb_variable_year    = 'variable_year'
    character(len=*) , parameter    :: orb_fixed_parameters = 'fixed_parameters'
    character(len=*) , parameter    :: subname = '(InitAttributes)'

    !----------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call shr_nuopc_memcheck(subname, 0, mastertask)

    !----------------------------------------------------------
    ! Initialize options for reproducible sums
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="reprosum_use_ddpdd", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) reprosum_use_ddpdd

    call NUOPC_CompAttributeGet(driver, name="reprosum_diffmax", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) reprosum_diffmax

    call NUOPC_CompAttributeGet(driver, name="reprosum_recompute", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) reprosum_recompute

    call shr_reprosum_setopts(repro_sum_use_ddpdd_in=reprosum_use_ddpdd, &
         repro_sum_rel_diff_max_in=reprosum_diffmax, repro_sum_recompute_in=reprosum_recompute)

    !----------------------------------------------------------
    ! Initialize freezing point calculation for all components
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="tfreeze_option", value=tfreeze_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_frz_freezetemp_init(tfreeze_option)

    !----------------------------------------------------------
    ! Initialize orbital related values
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="orb_mode", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mode

    call NUOPC_CompAttributeGet(driver, name="orb_iyear", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear

    call NUOPC_CompAttributeGet(driver, name="orb_iyear_align", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear_align

    call NUOPC_CompAttributeGet(driver, name="orb_obliq", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_obliq

    call NUOPC_CompAttributeGet(driver, name="orb_eccen", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_eccen

    call NUOPC_CompAttributeGet(driver, name="orb_mvelp", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mvelp

    if (trim(orb_mode) == trim(orb_fixed_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT) then
          write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
          write(logunit,*) trim(subname),' ERROR: fixed_year settings = ',orb_iyear
          write (msgstr, *) ' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    elseif (trim(orb_mode) == trim(orb_variable_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT .or. orb_iyear_align == SHR_ORB_UNDEF_INT) then
          write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
          write(logunit,*) trim(subname),' ERROR: variable_year settings = ',orb_iyear, orb_iyear_align
          write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    elseif (trim(orb_mode) == trim(orb_fixed_parameters)) then
       !-- force orb_iyear to undef to make sure shr_orb_params works properly
       orb_iyear = SHR_ORB_UNDEF_INT
       orb_iyear_align = SHR_ORB_UNDEF_INT
       if (orb_eccen == SHR_ORB_UNDEF_REAL .or. &
           orb_obliq == SHR_ORB_UNDEF_REAL .or. &
           orb_mvelp == SHR_ORB_UNDEF_REAL) then
          write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
          write(logunit,*) trim(subname),' ERROR: orb_eccen = ',orb_eccen
          write(logunit,*) trim(subname),' ERROR: orb_obliq = ',orb_obliq
          write(logunit,*) trim(subname),' ERROR: orb_mvelp = ',orb_mvelp
          write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    else
       write (msgstr, *) subname//' ERROR: invalid orb_mode '//trim(orb_mode)
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    endif

    call NUOPC_CompAttributeGet(driver, name='cpl_rootpe', value=cvalue, rc=rc)
    read(cvalue, *) rootpe_med
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_GridCompGet(driver, localPet=localPet, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! Determine orbital params
    if (trim(orb_mode) == trim(orb_variable_year)) then
       call ESMF_GridCompGet(driver, clock=clock, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_TimeGet(CurrTime, yy=year, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       orb_cyear = orb_iyear + (year - orb_iyear_align)
       call shr_orb_params(orb_cyear, orb_eccen, orb_obliq, orb_mvelp, &
            orb_obliqr, orb_lambm0, orb_mvelpp, localPet==rootpe_med )
    else
       call shr_orb_params(orb_iyear, orb_eccen, orb_obliq, orb_mvelp, &
            orb_obliqr, orb_lambm0, orb_mvelpp, localPet==rootpe_med )
    end if

    if (orb_eccen  == SHR_ORB_UNDEF_REAL .or. &
         orb_obliqr == SHR_ORB_UNDEF_REAL .or. &
         orb_mvelpp == SHR_ORB_UNDEF_REAL .or. &
         orb_lambm0 == SHR_ORB_UNDEF_REAL) then
       write (msgstr, *) subname//' ERROR: orb params incorrect'
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    endif

    ! Add updated orbital params to driver attributes

    call NUOPC_CompAttributeAdd(driver, attrList=(/'orb_obliqr', 'orb_lambm0', 'orb_mvelpp'/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    write(cvalue,*) orb_eccen
    call NUOPC_CompAttributeSet(driver, name="orb_eccen", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    write(cvalue,*) orb_obliqr
    call NUOPC_CompAttributeSet(driver, name="orb_obliqr", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    write(cvalue,*) orb_lambm0
    call NUOPC_CompAttributeSet(driver, name="orb_lambm0", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    write(cvalue,*) orb_mvelpp
    call NUOPC_CompAttributeSet(driver, name="orb_mvelpp", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO: need to update orbital parameters during run time - actually - each component needs to update its orbital
    ! parameters to be consistent

    !----------------------------------------------------------
    ! Initialize water vapor info
    !----------------------------------------------------------

    ! TODO: this does not seem to belong here - where should it go?

    call NUOPC_CompAttributeGet(driver, name="wv_sat_scheme", value=wv_sat_scheme, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (.not. shr_wv_sat_valid_idx(shr_wv_sat_get_scheme_idx(trim(wv_sat_scheme)))) then
       call shr_sys_abort(subname//': "'//trim(wv_sat_scheme)//'" is not a recognized saturation vapor pressure scheme name')
    end if
    if (.not. shr_wv_sat_set_default(wv_sat_scheme)) then
       call shr_sys_abort('Invalid wv_sat_scheme.')
    end if

    call NUOPC_CompAttributeGet(driver, name="wv_sat_transition_start", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) wv_sat_transition_start

    call shr_assert_in_domain(wv_sat_transition_start, &
         ge=0._SHR_KIND_R8, le=40._SHR_KIND_R8, &
         varname="wv_sat_transition_start", msg="Invalid transition temperature range.")

    call NUOPC_CompAttributeGet(driver, name="wv_sat_use_tables", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) wv_sat_use_tables

    call NUOPC_CompAttributeGet(driver, name="wv_sat_table_spacing", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) wv_sat_table_spacing

    ! A transition range averaging method in CAM is only valid for:
    ! -40 deg C <= T <= 0 deg C
    ! shr_wv_sat_mod itself checks for values with the wrong sign, but we
    ! have to check that the range is no more than 40 deg C here. Even
    ! though this is a CAM-specific restriction, it's not really likely
    ! that any other parameterization will be dealing with mixed-phase
    ! water below 40 deg C anyway.

    call shr_wv_sat_init(shr_const_tkfrz, shr_const_tktrip, wv_sat_transition_start, epsilo, errstring)
    if (errstring /= "") then
       call shr_sys_abort('shr_wv_sat_init: '//trim(errstring))
    end if

    ! The below produces internal lookup tables in the range 175-374K for
    ! liquid water, and 125-274K for ice, with a resolution set by the
    ! option wv_sat_table_spacing.
    ! In theory these ranges could be specified in the namelist, but in
    ! practice users will want to change them *very* rarely if ever, which
    ! is why only the spacing is in the namelist.

    if (wv_sat_use_tables) then
       liquid_spec = ShrWVSatTableSpec(ceiling(200._SHR_KIND_R8/wv_sat_table_spacing), 175._SHR_KIND_R8, wv_sat_table_spacing)
       ice_spec    = ShrWVSatTableSpec(ceiling(150._SHR_KIND_R8/wv_sat_table_spacing), 125._SHR_KIND_R8, wv_sat_table_spacing)
       mixed_spec  = ShrWVSatTableSpec(ceiling(250._SHR_KIND_R8/wv_sat_table_spacing), 125._SHR_KIND_R8, wv_sat_table_spacing)
       call shr_wv_sat_make_tables(liquid_spec, ice_spec, mixed_spec)
    end if

    !----------------------------------------------------------
    ! Set single_column flags
    ! If in single column mode, overwrite flags according to focndomain file
    ! in ocn_in namelist. SCAM can reset the "present" flags for lnd,
    ! ocn, ice, rof, and flood.
    !----------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="single_column", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) single_column

    ! NOTE: cam stand-alone aqua-planet model will no longer be supported here - only the data model aqua-planet
    ! will be supported
    if (single_column) then

       call NUOPC_CompAttributeGet(driver, name="scmlon", value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) scmlon

       call NUOPC_CompAttributeGet(driver, name="scmlat", value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) scmlat

       call seq_comm_getinfo(OCNID(ens1), mpicom=mpicom_OCNID)

       ! TODO: Single column mode needs to be re-implemented - previously all of the xxx_present flags were set
       ! in med_infodata calls, reset here and the copied back into med_infodata - this is no longer the case
       ! call shr_scam_checkSurface(scmlon, scmlat, &
       !      OCNID(ens1), mpicom_OCNID,            &
       !      lnd_present=lnd_present,              &
       !      ocn_present=ocn_present,              &
       !      ice_present=ice_present,              &
       !      rof_present=rof_present,              &
       !      flood_present=flood_present,          &
       !      rofice_present=rofice_present)

    endif

  end subroutine InitAttributes

  !================================================================================

  subroutine CheckAttributes( driver, rc )

    ! !DESCRIPTION: Check that input driver config values have reasonable values
    use shr_sys_mod           , only : shr_sys_abort
    use ESMF, only : ESMF_GridComp, ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use NUOPC, only : NUOPC_CompAttributeGet

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_GridComp) , intent(inout) :: driver
    integer             , intent(out)   :: rc

    !----- local -----
    character(SHR_KIND_CL) :: cvalue         ! temporary
    character(SHR_KIND_CL) :: start_type     ! Type of startup
    character(SHR_KIND_CL) :: rest_case_name ! Short case identification
    character(SHR_KIND_CS) :: logFilePostFix ! postfix for output log files
    character(SHR_KIND_CL) :: outPathRoot    ! root for output log files
    character(SHR_KIND_CS) :: cime_model
    integer :: dbrc
    character(len=*), parameter :: subname = '(driver_attributes_check) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call NUOPC_CompAttributeGet(driver, name="cime_model", value=cime_model, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if ( trim(cime_model) /= 'cesm') then
       call shr_sys_abort( subname//': cime_model must be set to cesm, aborting')
    end if

    ! --- LogFile ending name -----
    call NUOPC_CompAttributeGet(driver, name="logFilePostFix", value=logFilePostFix, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( len_trim(logFilePostFix) == 0 ) then
       call shr_sys_abort( subname//': logFilePostFix  must be set to something not blank' )
    end if

    ! --- Output path root directory -----
    call NUOPC_CompAttributeGet(driver, name="outPathRoot", value=outPathRoot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( len_trim(outPathRoot) == 0 ) then
       call shr_sys_abort( subname//': outPathRoot  must be set' )
    end if
    if ( index(outPathRoot, "/", back=.true.) /= len_trim(outPathRoot) ) then
       call shr_sys_abort( subname//': outPathRoot must end with a slash' )
    end if

    ! --- Case name and restart case name ------
    ! call NUOPC_CompAttributeGet(driver, name="rest_case_name", value=rest_case_name, rc=rc)
    ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! if ((trim(start_type) == start_type_cont ) .and. (trim(case_name)  /= trim(rest_case_name))) then
    !    write(logunit,'(10a)') subname,' case_name =',trim(case_name),':',' rest_case_name =',trim(rest_case_name),':'
    !    call shr_sys_abort(subname//': invalid continue restart case name = '//trim(rest_case_name))
    ! endif

  end subroutine CheckAttributes

  !===============================================================================

  subroutine AddAttributes(gcomp, driver, config, compid, compname, inst_suffix, rc)

    ! Add specific set of attributes to gcomp from driver attributes
    use ESMF, only : ESMF_GridComp, ESMF_Config, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF, only : ESMF_LogFoundAllocError, ESMF_ConfigGetLen, ESMF_ConfigGetAttribute
    use NUOPC, only : NUOPC_CompAttributeAdd, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
    use seq_comm_mct, only : seq_comm_inst, seq_comm_name, seq_comm_suffix

    ! input/output parameters
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_GridComp) , intent(in)    :: driver
    type(ESMF_Config)   , intent(inout)    :: config
    integer             , intent(in)    :: compid
    character(len=*)    , intent(in)    :: compname
    character(len=*)    , intent(in)    :: inst_suffix
    integer             , intent(inout) :: rc

    ! locals
    integer                     :: n
    integer                     :: stat
    character(len=SHR_KIND_CL)  :: cvalue
    character(len=32), allocatable :: compLabels(:)
    integer         , parameter :: nattrlist = 5
    character(len=*), parameter :: attrList(nattrlist) = &
         (/"read_restart", "orb_eccen", "orb_obliqr", "orb_lambm0", "orb_mvelpp"/)
    integer :: dbrc
    character(len=*), parameter :: subname = "(esm.F90:AddAttributes)"
    !-------------------------------------------

    rc = ESMF_Success
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    ! First add compid to gcomp attributes

    write(cvalue,*) compid
    call NUOPC_CompAttributeAdd(gcomp, attrList=(/'MCTID'/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, name='MCTID', value=trim(cvalue), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Now add all the other attributes in AttrList (which have already been added to driver attributes)
    call NUOPC_CompAttributeAdd(gcomp, attrList=attrList, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1,nattrlist
       call NUOPC_CompAttributeGet(driver, name=trim(attrList(n)), value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name=trim(attrList(n)), value=trim(cvalue), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    ! Now add component specific attributes
    call ReadAttributes(gcomp, config, trim(compname)//"_attributes::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ReadAttributes(gcomp, config, "ALLCOMP_attributes::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(gcomp, config, trim(compname)//"_modelio"//trim(inst_suffix)//"::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(gcomp, config, "CLOCK_attributes::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (compname == 'MED') then

       call ReadAttributes(gcomp, config, "MED_history_attributes::", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ReadAttributes(gcomp, config, "FLDS_attributes::", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeAdd(gcomp, &
            attrList=(/'atm_present','lnd_present','ocn_present','ice_present',&
                       'rof_present','wav_present','glc_present','med_present'/), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       componentCount = ESMF_ConfigGetLen(config,label="CESM_component_list:", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       allocate(compLabels(componentCount), stat=stat)
       if (ESMF_LogFoundAllocError(statusToCheck=stat, msg="Allocation of compLabels failed.", &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)) return

       call ESMF_ConfigGetAttribute(config, valueList=compLabels, label="CESM_component_list:", count=componentCount, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       med_present = "false"
       atm_present = "false"
       lnd_present = "false"
       ocn_present = "false"
       ice_present = "false"
       rof_present = "false"
       wav_present = "false"
       glc_present = "false"
       do n=1, componentCount
          if (trim(compLabels(n)) == "MED") med_present = "true"
          if (trim(compLabels(n)) == "ATM") atm_present = "true"
          if (trim(compLabels(n)) == "LND") lnd_present = "true"
          if (trim(compLabels(n)) == "OCN") ocn_present = "true"
          if (trim(compLabels(n)) == "ICE") ice_present = "true"
          if (trim(compLabels(n)) == "ROF") rof_present = "true"
          if (trim(compLabels(n)) == "WAV") wav_present = "true"
          if (trim(compLabels(n)) == "GLC") glc_present = "true"
       enddo
       deallocate(compLabels)
       call NUOPC_CompAttributeSet(gcomp, name="atm_present", value=atm_present, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name="lnd_present", value=lnd_present, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name="ocn_present", value=ocn_present, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name="ice_present", value=ice_present, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name="rof_present", value=rof_present, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name="wav_present", value=wav_present, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name="glc_present", value=glc_present, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name="med_present", value=med_present, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    ! inst_name and inst_index are no longer required for cime internal components
    call NUOPC_CompAttributeAdd(gcomp, attrList=(/'inst_name','inst_index'/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeSet(gcomp, name='inst_name', value=trim(seq_comm_name(compid)), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    write(cvalue,*) seq_comm_inst(compid)
    call NUOPC_CompAttributeSet(gcomp, name='inst_index', value=trim(cvalue), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if(len_trim(inst_suffix) > 0) then
       call NUOPC_CompAttributeAdd(gcomp, attrList=(/'inst_suffix'/), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name='inst_suffix', value=inst_suffix, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
  end subroutine AddAttributes

  !================================================================================

  subroutine ReadAttributes(gcomp, config, label, relaxedflag, formatprint, rc)
    use ESMF, only : ESMF_GridComp, ESMF_Config, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use NUOPC, only : NUOPC_FreeFormatCreate, NUOPC_FreeFormatPrint, NUOPC_CompAttributeIngest
    use NUOPC, only : NUOPC_FreeFormatDestroy, NUOPC_FreeFormat

    type(ESMF_GridComp) , intent(inout)        :: gcomp
    type(ESMF_Config)   , intent(in)           :: config
    character(len=*)    , intent(in)           :: label
    logical             , intent(in), optional :: relaxedflag
    logical             , intent(in), optional :: formatprint
    integer             , intent(inout)        :: rc

    ! local variables
    type(NUOPC_FreeFormat)      :: attrFF
    integer :: dbrc
    character(len=*), parameter :: subname = "(esm.F90:ReadAttributes)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    if (present(relaxedflag)) then
       attrFF = NUOPC_FreeFormatCreate(config, label=trim(label), relaxedflag=.true., rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       attrFF = NUOPC_FreeFormatCreate(config, label=trim(label), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call NUOPC_CompAttributeIngest(gcomp, attrFF, addFlag=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (present (formatprint)) then
       if (mastertask) then
          call NUOPC_FreeFormatPrint(attrFF, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine ReadAttributes

! This empty InitAdvertise is needed because it overrides the behavior
! of the default InitAdvertise inside the generic NUOPC_Driver.F90. The
! default behavior tries to mirror the fields up the hierarchy (i.e., up
! to the ensemble driver). This would be used if we needed to
! communicate between the ensemble members. Since we do not need that
! right now, we turn it off with this empty subroutine.

  subroutine InitAdvertize(driver, importState, exportState, clock, rc)
    use ESMF, only : ESMF_GridComp, ESMF_State, ESMF_Clock
    use ESMF, only : ESMF_LogWrite, ESMF_SUCCESS, ESMF_LOGMSG_INFO
    type(ESMF_GridComp)  :: driver
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    integer :: dbrc
    character(len=*), parameter :: subname = "(esm.F90:InitAdvertize)"

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine InitAdvertize

  subroutine esm_init_pelayout(driver, maxthreads)
    use ESMF, only : ESMF_GridComp, ESMF_GridCompGet, ESMF_VM, ESMF_VMGet
    use ESMF, only : ESMF_LogWrite, ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_Config
    use ESMF, only : ESMF_ConfigGetLen, ESMF_LogFoundAllocError, ESMF_ConfigGetAttribute
    use ESMF, only : ESMF_RC_NOT_VALID, ESMF_LogSetError
    use ESMF, only : ESMF_GridCompIsPetLocal, ESMF_MethodAdd
    use NUOPC, only : NUOPC_CompAttributeGet
    use NUOPC_Driver, only: NUOPC_DriverAddComp
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_utils_mod, only : shr_nuopc_abort
    use shr_string_mod, only : toLower => shr_string_toLower
    use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag, CS, CL
    use atm_comp_nuopc   , only : ATMSetServices => SetServices
    use ice_comp_nuopc   , only : ICESetServices => SetServices
    use lnd_comp_nuopc   , only : LNDSetServices => SetServices
    use ocn_comp_nuopc   , only : OCNSetServices => SetServices
    use wav_comp_nuopc   , only : WAVSetServices => SetServices
    use rof_comp_nuopc   , only : ROFSetServices => SetServices
    use glc_comp_nuopc   , only : GLCSetServices => SetServices
    use MED              , only : MEDSetServices => SetServices
    use mpi, only  : MPI_COMM_NULL
    ! These should be removed
    use seq_comm_mct, only: GLOID, CPLID, ATMID, LNDID, OCNID, ICEID, GLCID, ROFID, WAVID, ESPID
    use seq_comm_mct, only: seq_comm_setcomm
    use mct_mod, only : mct_world_init

    type(ESMF_GridComp) :: driver
    integer, intent(out) :: maxthreads ! maximum number of threads any component
    type(ESMF_GridComp) :: child
    type(ESMF_VM) :: vm
    type(ESMF_Config) :: config
    integer :: componentcount
    integer :: PetCount
    integer :: LocalPet
    integer :: ntasks, rootpe, nthrds, stride
    integer :: ntask, cnt
    integer :: rc
    integer :: i
    integer :: stat
    character(len=32), allocatable :: compLabels(:)
    character(CS) :: namestr
    character(CL) :: msgstr
    integer :: pelist(3,1)       ! start, stop, stride for group (mct initialization to be removed)
    integer, allocatable :: petlist(:)
    integer, pointer :: comms(:), comps(:)
    integer :: Global_Comm
    logical :: isPresent
    character(len=5) inst_suffix
    character(CL)          :: cvalue

    character(len=*), parameter :: subname = "(esm_pelayout.F90:esm_init_pelayout)"

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    maxthreads = 1
    call ESMF_GridCompGet(driver, vm=vm, config=config, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(driver, config, "PELAYOUT_attributes::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, petCount=petCount, localPet=localPet, mpiCommunicator=Global_Comm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    componentCount = ESMF_ConfigGetLen(config,label="CESM_component_list:", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(compLabels(componentCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg="Allocation of compLabels failed.", &
         line=__LINE__, file=u_FILE_u, rcToReturn=rc)) return

    call ESMF_ConfigGetAttribute(config, valueList=compLabels, label="CESM_component_list:", count=componentCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(driver, name="inst_suffix", value=inst_suffix, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       inst_suffix = ""
    endif

! Initialize mct ID's - these should be removed
    GLOID = 1
    pelist(1,1) = 0
    pelist(2,1) = PetCount - 1
    pelist(3,1) = 1

    call seq_comm_setcomm(GLOID, pelist, iname='GLOBAL', comm_in=Global_Comm)
    allocate(comms(componentCount+1), comps(componentCount+1))
    comps(1) = GLOID
    comms(1) = Global_Comm

    do i=1,componentCount
       namestr = toLower(compLabels(i))
       if (namestr == 'med') namestr = 'cpl'
       call NUOPC_CompAttributeGet(driver, name=trim(namestr)//'_ntasks', value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) ntasks

       if (ntasks < 0 .or. ntasks > PetCount) then
          write (msgstr, *) "Invalid NTASKS value specified for component: ",namestr, ' ntasks: ',ntasks
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
       endif

       call NUOPC_CompAttributeGet(driver, name=trim(namestr)//'_nthreads', value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) nthrds

       if(nthrds > maxthreads) maxthreads = nthrds

       call NUOPC_CompAttributeGet(driver, name=trim(namestr)//'_rootpe', value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) rootpe
       if (rootpe < 0 .or. rootpe > PetCount) then
          write (msgstr, *) "Invalid Rootpe value specified for component: ",namestr, ' rootpe: ',rootpe
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
       endif
       if(rootpe+ntasks > PetCount) then
          write (msgstr, *) "Invalid pelayout value specified for component: ",namestr, ' rootpe+ntasks: ',rootpe+ntasks
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
       endif

       call NUOPC_CompAttributeGet(driver, name=trim(namestr)//'_pestride', value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stride
       if (stride < 1 .or. rootpe+ntasks*stride > PetCount) then
          write (msgstr, *) "Invalid pestride value specified for component: ",namestr, ' rootpe: ',rootpe, ' pestride: ', stride
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
       endif

       if (allocated(petlist) .and. size(petlist) .ne. ntasks) then
          deallocate(petlist)
       endif
       if(.not. allocated(petlist)) then
          allocate(petlist(ntasks))
       endif

       cnt=1
       do ntask = rootpe, (rootpe+ntasks*stride)-1, stride
          petlist(cnt) = ntask
          cnt=cnt+1
       enddo
! Initialize mct comm stuff - to be removed
       comps(i+1) = i+1
       if (trim(compLabels(i)) .eq. 'MED') then
          CPLID = i+1
          pelist(1,1) = rootpe
          pelist(2,1) = rootpe+ntasks-1
          pelist(3,1) = stride
          call seq_comm_setcomm(CPLID, pelist, nthreads=nthrds, iname='CPL')
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), MEDSetServices, petList=petlist, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       elseif(trim(compLabels(i)) .eq. 'ATM') then
          ATMID(1) = i+1
          pelist(1,1) = rootpe
          pelist(2,1) = rootpe+ntasks-1
          pelist(3,1) = stride
          call seq_comm_setcomm(ATMID(1), pelist, nthreads=nthrds, iname=trim(compLabels(i)))
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), ATMSetServices, petList=petlist, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) call shr_nuopc_abort()
       elseif(trim(compLabels(i)) .eq. 'LND') then
          LNDID(1) = i+1
          pelist(1,1) = rootpe
          pelist(2,1) = rootpe+ntasks-1
          pelist(3,1) = stride
          call seq_comm_setcomm(LNDID(1), pelist, nthreads=nthrds, iname=trim(compLabels(i)))
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), LNDSetServices, PetList=petlist, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       elseif(trim(compLabels(i)) .eq. 'OCN') then
          OCNID(1) = i+1
          pelist(1,1) = rootpe
          pelist(2,1) = rootpe+ntasks-1
          pelist(3,1) = stride
          call seq_comm_setcomm(OCNID(1), pelist, nthreads=nthrds, iname=trim(compLabels(i)))
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), OCNSetServices, PetList=petlist, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        elseif(trim(compLabels(i)) .eq. 'ICE') then
          ICEID(1) = i+1
          pelist(1,1) = rootpe
          pelist(2,1) = rootpe+ntasks-1
          pelist(3,1) = stride
          call seq_comm_setcomm(ICEID(1), pelist, nthreads=nthrds, iname=trim(compLabels(i)))
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), ICESetServices, PetList=petlist, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       elseif(trim(compLabels(i)) .eq. 'GLC') then
          GLCID(1) = i+1
          pelist(1,1) = rootpe
          pelist(2,1) = rootpe+ntasks-1
          pelist(3,1) = stride
          call seq_comm_setcomm(GLCID(1), pelist, nthreads=nthrds, iname=trim(compLabels(i)))
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), GLCSetServices, PetList=petlist, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       elseif(trim(compLabels(i)) .eq. 'ROF') then
          ROFID(1) = i+1
          pelist(1,1) = rootpe
          pelist(2,1) = rootpe+ntasks-1
          pelist(3,1) = stride
          call seq_comm_setcomm(ROFID(1), pelist, nthreads=nthrds, iname=trim(compLabels(i)))
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), ROFSetServices, PetList=petlist, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       elseif(trim(compLabels(i)) .eq. 'WAV') then
          WAVID(1) = i+1
          pelist(1,1) = rootpe
          pelist(2,1) = rootpe+ntasks-1
          pelist(3,1) = stride
          call seq_comm_setcomm(WAVID(1), pelist, nthreads=nthrds, iname=trim(compLabels(i)))
          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), WAVSetServices, PetList=petlist, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
!       elseif(trim(compLabels(i)) .eq. 'ESP') then
!          ESPID(1) = i+1
!          pelist(1,1) = rootpe
!          pelist(2,1) = rootpe+ntasks-1
!          pelist(3,1) = stride
!          call seq_comm_setcomm(ESPID(1), pelist, nthreads=nthrds, iname=trim(compLabels(i)))
!          call NUOPC_DriverAddComp(driver, trim(compLabels(i)), ESPSetServices, PetList=petlist, comp=child, rc=rc)
       endif
       call AddAttributes(child, driver, config, i+1, trim(compLabels(i)), inst_suffix, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) call shr_nuopc_abort()
       if (ESMF_GridCompIsPetLocal(child, rc=rc)) then
          call ESMF_GridCompGet(child, vm=vm, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_VMGet(vm, mpiCommunicator=comms(i+1), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Attach methods for handling reading/writing of restart pointer file
          call ESMF_MethodAdd(child, label="GetRestartFileToWrite", &
               userRoutine=GetRestartFileToWrite, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_MethodAdd(child, label="GetRestartFileToRead", &
               userRoutine=GetRestartFileToRead, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       else
          comms(i+1) = MPI_COMM_NULL
       endif
    enddo
    call mct_world_init(componentCount, GLOBAL_COMM, comms, comps)
    deallocate(petlist, comms, comps)

  end subroutine esm_init_pelayout

  subroutine esm_finalize(driver, rc)
    use ESMF,  only : ESMF_GridComp, ESMF_GridCompGet, ESMF_VM, ESMF_VMGet
    use ESMF,  only : ESMF_SUCCESS
    use NUOPC, only : NUOPC_CompAttributeGet
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_ChkErr
    use perf_mod, only : t_prf, t_finalizef
    use med_constants_mod, only : CL

    type(ESMF_GridComp) :: driver
    integer, intent(out) :: rc
    character(CL) :: timing_dir        ! timing directory
    character(len=5) :: inst_suffix
    logical :: isPresent

    rc = ESMF_SUCCESS

    call NUOPC_CompAttributeGet(driver, name="timing_dir",value=timing_dir, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(driver, name="inst_suffix", value=inst_suffix, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       inst_suffix = ""
    endif

    call t_finalizef()

  end subroutine esm_finalize

  !================================================================================

  ! Method to be attached to components to handle
  ! CESM specific ways of writing restart files
  !
  ! This is used with MOM6 now and may need to be
  ! extended or generalized to other components
  subroutine GetRestartFileToWrite(gcomp, rc)
    use ESMF,         only: ESMF_GridComp, ESMF_GridCompGet
    use ESMF,         only: ESMF_LogSetError, ESMF_SUCCESS, ESMF_RC_FILE_OPEN
    use ESMF,         only: ESMF_RC_ATTR_NOTSET
    use ESMF,         only: ESMF_Time, ESMF_TimeGet
    use ESMF,         only: ESMF_Clock, ESMF_ClockGetNextTime
    use ESMF,         only: ESMF_VM, ESMF_VMGet
    use ESMF,         only: ESMF_MAXSTR, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use NUOPC,        only: NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
    use shr_file_mod, only: shr_file_getUnit, shr_file_freeUnit

    type(ESMF_GridComp)                 :: gcomp
    integer            , intent(out)    :: rc

    type(ESMF_VM)           :: vm
    integer                 :: localPet, nu, iostat
    type(ESMF_Clock)        :: clock
    type(ESMF_Time)         :: nextTime
    character(ESMF_MAXSTR)  :: casename, restartname
    logical                 :: isPresent, isSet
    integer                 :: year, month, day, seconds
    character(len=*), parameter :: subname='GetRestartFileToWrite'

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    rc = ESMF_SUCCESS

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=casename, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (.not. isPresent .or. .not. isSet) then
       call ESMF_LogSetError(ESMF_RC_ATTR_NOTSET, &
            msg=subname//": case_name attribute must be set to generate restart filename",  &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
    endif

    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Need to use next time step since clock is
    ! not advanced until the end of the time interval
    call ESMF_ClockGetNextTime(clock, nextTime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(nextTime, yy=year, mm=month, dd=day, s=seconds, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    write(restartname,'(A,".mom6.r.",I4.4,"-",I2.2,"-",I2.2,"-",I5.5)') &
         trim(casename), year, month, day, seconds

    call NUOPC_CompAttributeSet(gcomp, name="RestartFileToWrite", &
         value=trim(restartname), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (localPet == 0) then
       ! Write name of restart file in the rpointer file
       ! This is currently hard-coded for the ocean
       nu = shr_file_getUnit()
       open(nu, file='rpointer.ocn', form='formatted', &
            status='unknown', iostat=iostat)
       if (iostat /= 0) then
          call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
               msg=subname//' ERROR opening rpointer.ocn', &
               line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return
       endif
       write(nu,'(a)') trim(restartname)//'.nc'
       close(nu)
       call shr_file_freeUnit(nu)
    endif
    call ESMF_LogWrite(trim(subname)//": returning", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine GetRestartFileToWrite

  !================================================================================

  subroutine GetRestartFileToRead(gcomp, rc)

    use ESMF,         only: ESMF_GridComp, ESMF_GridCompGet
    use ESMF,         only: ESMF_LogSetError, ESMF_SUCCESS, ESMF_RC_FILE_OPEN
    use ESMF,         only: ESMF_RC_FILE_READ
    use ESMF,         only: ESMF_VM, ESMF_VMGet, ESMF_VMBroadcast
    use ESMF,         only: ESMF_MAXSTR, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use NUOPC,        only: NUOPC_CompAttributeSet
    use shr_file_mod, only: shr_file_getUnit, shr_file_freeUnit

    type(ESMF_GridComp)                 :: gcomp
    integer            , intent(out)    :: rc

    type(ESMF_VM)           :: vm
    integer                 :: localPet, readunit, iostat
    logical                 :: is_restart
    character(ESMF_MAXSTR)  :: restartname
    character(len=*), parameter :: subname='GetRestartFileToRead'

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=rc)
    endif

    is_restart = IsRestart(gcomp, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (is_restart) then
       restartname = ""

       call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_VMGet(vm, localPet=localPet, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (localPet == 0) then
          readunit = shr_file_getUnit()
          ! this hard coded for rpointer.ocn right now
          open(readunit, file='rpointer.ocn', form='formatted', status='old', iostat=iostat)
          if (iostat /= 0) then
             call ESMF_LogSetError(ESMF_RC_FILE_OPEN, msg=subname//' ERROR opening rpointer.ocn', &
                  line=__LINE__, file=u_FILE_u, rcToReturn=rc)
             return
          endif
          read(readunit,'(a)', iostat=iostat) restartname
          if (iostat /= 0) then
             call ESMF_LogSetError(ESMF_RC_FILE_READ, msg=subname//' ERROR reading rpointer.ocn', &
                  line=__LINE__, file=u_FILE_u, rcToReturn=rc)
             return
          endif
          close(readunit)
       endif

       ! broadcast attribute set on master task to all tasks
       call ESMF_VMBroadcast(vm, restartname, count=ESMF_MAXSTR-1, rootPet=0, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       !write(logunit,*) trim(subname)//":restartfile after broadcast = "//trim(restartfile)

       call NUOPC_CompAttributeSet(gcomp, name='RestartFileToRead', &
            value=trim(restartname), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    call ESMF_LogWrite(trim(subname)//": returning", ESMF_LOGMSG_INFO, rc=rc)

  end subroutine GetRestartFileToRead

end module ESM
