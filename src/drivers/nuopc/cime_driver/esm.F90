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
  character(len=32), allocatable :: compLabels(:)
  character(len=8)               :: atm_present, lnd_present, ocn_present
  character(len=8)               :: ice_present, rof_present, wav_present
  character(len=8)               :: glc_present, med_present
  character(*), parameter        :: nlfilename = "drv_in" ! input namelist filename
  character(*), parameter        :: u_FILE_u = &
       __FILE__

  public  :: SetServices
  public :: ReadAttributes

  private :: SetModelServices
  private :: SetRunSequence
  private :: ModifyCplLists
  private :: InitAttributes
  private :: CheckAttributes
  private :: AddAttributes
  private :: InitAdvertize

!================================================================================
  contains
!================================================================================

  subroutine SetServices(driver, rc)

    use NUOPC        , only : NUOPC_CompDerive, NUOPC_CompSpecialize, NUOPC_CompSetInternalEntryPoint
    use NUOPC_Driver , only : driver_routine_SS             => SetServices
    use NUOPC_Driver , only : driver_label_SetModelServices => label_SetModelServices
    use NUOPC_Driver , only : driver_label_SetRunSequence   => label_SetRunSequence
    use ESMF         , only : ESMF_GridComp, ESMF_Config, ESMF_GridCompSet, ESMF_ConfigLoadFile
    use ESMF         , only : ESMF_ConfigCreate, ESMF_METHOD_INITIALIZE
    use ESMF         , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO

    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Config) :: config
    integer           :: dbrc
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

    ! Create, open and set the config
    config = ESMF_ConfigCreate(rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ConfigLoadFile(config, "cesm.runconfig", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompSet(driver, config=config, rc=rc)
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
    use NUOPC                 , only : NUOPC_CompSetInternalEntryPoint, NUOPC_CompAttributeGet
    use NUOPC_Driver          , only : NUOPC_DriverAddComp
    use seq_comm_mct          , only : CPLID, GLOID, ATMID, LNDID, OCNID, ICEID, GLCID, ROFID, WAVID, ESPID
    use seq_comm_mct          , only : seq_comm_init, seq_comm_printcomms, seq_comm_petlist
    use seq_comm_mct          , only : seq_comm_getinfo => seq_comm_setptrs
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_Clock_TimePrint
    use shr_file_mod          , only : shr_file_getlogunit, shr_file_setLogunit
    use shr_file_mod          , only : shr_file_getlogLevel, shr_file_setLogLevel
    use shr_file_mod          , only : shr_file_getUnit, shr_file_freeUnit
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
    use shr_nuopc_time_mod    , only : shr_nuopc_time_clockInit
    use mpi                   , only : MPI_COMM_WORLD

    type(ESMF_GridComp)    :: driver
    integer, intent(out)   :: rc


    ! local variables
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
    integer                :: localPet
    logical                :: is_set
    character(SHR_KIND_CS) :: cvalue
    character(len=5) inst_suffix
    character(len=512)     :: diro
    character(len=512)     :: logfile
    integer                :: shrlogunit  ! original log unit
    integer                :: shrloglev   ! original log level
    integer                :: global_comm
    logical                :: iamroot_med ! mediator masterproc
    logical                :: isPresent
    integer                :: dbrc
    character(len=*), parameter    :: subname = "(esm.F90:SetModelServices)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    !-------------------------------------------
    ! Get the config and vm objects from the driver
    !-------------------------------------------

    call ESMF_GridCompGet(driver, vm=vm, config=config, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=localPet, mpiCommunicator=global_comm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (localPet == 0) then
       mastertask=.true.
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

    allocate(compLabels(componentCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg="Allocation of compLabels failed.", &
          line=__LINE__, file=u_FILE_u, rcToReturn=rc)) return

    call ESMF_ConfigGetAttribute(config, valueList=compLabels, label="CESM_component_list:", count=componentCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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

    ! NOTE: if pio_async_interface is true global_comm is MPI_COMM_NULL on the servernodes
    ! and server nodes do not return from shr_pio_init2
    ! NOTE: if (global_comm /= MPI_COMM_NULL) then the following call also initializes
    ! MCT which is still needed for some models
    call seq_comm_init(global_comm, nlfilename)

    call NUOPC_CompAttributeGet(driver, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(driver, name="inst_suffix", value=inst_suffix, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       inst_suffix = ""
    endif

    !-------------------------------------------
    ! Perform restarts if appropriate
    !-------------------------------------------

    ! finish the pio initialization (this calls shr_pio_init2)
    call InitPIO(driver, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call InitRestart(driver, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (pio_file_is_open(pioid)) then
       call pio_closefile(pioid)
    end if

    !-------------------------------------------
    ! Initialize driver clock
    !-------------------------------------------

    call shr_nuopc_time_clockInit(driver, logunit, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------
    ! Initialize other attributes (after initializing driver clock)
    !-------------------------------------------

    call InitAttributes(driver, rc)

    !-------------------------------------------
    ! Determine information for each component and add to the driver
    !-------------------------------------------

    do n=1, componentCount

      !--- construct component prefix
      prefix=trim(compLabels(n))

      !--- read in model instance name
      call ESMF_ConfigGetAttribute(config, model, label=trim(prefix)//"_model:", default="none", rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--- check that there was a model instance specified
      if (trim(model) == "none") then
         write (msgstr, *) "No model was specified for component: ",trim(prefix)
         call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
         return
      endif

#if (1 == 0)
      ! read in petList bounds
      call ESMF_ConfigGetAttribute(config, petListBounds, label=trim(prefix)//"_petlist_bounds:", default=-1, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=trim(name)//":"//__FILE__)) return
      ! handle the default situation
      if (petListBounds(1)==-1 .or. petListBounds(2)==-1) then
         petListBounds(1) = 0
         petListBounds(2) = petCount - 1
      endif
      ! set petList for this component
      allocate(petList(petListBounds(2)-petListBounds(1)+1))
      do j=petListBounds(1), petListBounds(2)
         petList(j-petListBounds(1)+1) = j ! PETs are 0 based
      enddo
#endif

      !--------
      ! ATM
      !--------

      if (trim(prefix) == "ATM") then

        compid = ATMID(1)
        call seq_comm_petlist(compid, petList)

        is_set = .false.
        call NUOPC_DriverAddComp(driver, "ATM", ATMSetServices, petList=petList, comp=child, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        is_set = .true.

        if (.not. is_set) then
           call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=subname//' model unavailable = ATM:'//trim(model), &
                line=__LINE__, file=u_FILE_u, rcToReturn=rc)
           return
        end if

        call AddAttributes(child, driver, config, compid, 'ATM', inst_suffix, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! OCN
      !--------

      elseif (trim(prefix) == "OCN") then

         compid = OCNID(1)
         call seq_comm_petlist(compid,petList)

         is_set = .false.
         call NUOPC_DriverAddComp(driver, "OCN", OCNSetServices, petList=petList, comp=child, rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         is_set = .true.
         if (.not. is_set) then
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=subname//' model unavailable = OCN:'//trim(model), &
                 line=__LINE__, file=u_FILE_u, rcToReturn=rc)
            return
         end if

        call AddAttributes(child, driver, config, compid, 'OCN', inst_suffix, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! ICE
      !--------

      elseif (trim(prefix) == "ICE") then

        compid = ICEID(1)
        call seq_comm_petlist(compid, petList)

        is_set = .false.
        call NUOPC_DriverAddComp(driver, "ICE", ICESetServices, petList=petList, comp=child, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        is_set = .true.
        if (.not. is_set) then
           call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=subname//' model unavailable = ICE:'//trim(model), &
                line=__LINE__, file=u_FILE_u, rcToReturn=rc)
           return
        end if

        call AddAttributes(child, driver, config, compid, 'ICE', inst_suffix, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! LND
      !--------

      elseif (trim(prefix) == "LND") then

        compid = LNDID(1)
        call seq_comm_petlist(compid, petList)

        is_set = .false.
        call NUOPC_DriverAddComp(driver, "LND", LNDSetServices, petList=petList, comp=child, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        is_set = .true.
        if (.not. is_set) then
           call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=subname//' model unavailable = LND:'//trim(model), &
                line=__LINE__, file=u_FILE_u, rcToReturn=rc)
           return
        end if

        call AddAttributes(child, driver, config, compid, 'LND', inst_suffix, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! WAV
      !--------

      elseif (trim(prefix) == "WAV") then

        compid = WAVID(1)
        call seq_comm_petlist(compid, petList)

        is_set = .false.
        call NUOPC_DriverAddComp(driver, "WAV", WAVSetServices, petList=petList, comp=child, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        is_set = .true.
        if (.not. is_set) then
           call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=subname//' model unavailable = WAV:'//trim(model), &
                line=__LINE__, file=u_FILE_u, rcToReturn=rc)
           return
        end if

        call AddAttributes(child, driver, config, compid, 'WAV', inst_suffix, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! GLC
      !--------

      elseif (trim(prefix) == "GLC") then

        compid = GLCID(1)
        call seq_comm_petlist(compid, petList)

        call NUOPC_DriverAddComp(driver, "GLC", GLCSetServices, petList=petList, comp=child, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        call AddAttributes(child, driver, config, compid, 'GLC', inst_suffix, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! ROF
      !--------

      elseif (trim(prefix) == "ROF") then

        compid = ROFID(1)
        call seq_comm_petlist(compid, petList)

        call NUOPC_DriverAddComp(driver, "ROF", ROFSetServices, petList=petList, comp=child, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        call AddAttributes(child, driver, config, compid, 'ROF', inst_suffix, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      !--------
      ! MED
      !--------

      elseif (trim(prefix) == "MED") then

        compid = CPLID
        call seq_comm_petlist(compid, petList)

        if (trim(model) == "cesm") then
          call NUOPC_DriverAddComp(driver, "MED", med_SS, petList=petList, comp=child, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        else
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=subname//' invalid model = MED:'//trim(model), &
               line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif

        call AddAttributes(child, driver, config, compid, 'MED', inst_suffix, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        iamroot_med = localPet == petList(1)
        call shr_file_getLogUnit (shrlogunit)
        if (iamroot_med) then
           logunit = shr_file_getUnit()
           open(logunit,file=trim(diro)//"/"//trim(logfile))
        else
           logUnit = shrlogunit
        endif

        call shr_file_getLogLevel(shrloglev)
        call shr_file_setLogLevel(max(shrloglev,1))
        call shr_file_setLogUnit (logunit)
        if(iamroot_med) print *,__FILE__,__LINE__,logunit, shrlogunit

        ! Print out present flags to mediator log file
        if (iamroot_med) then
           write(logunit,*) trim(subname)//":atm_present="//trim(atm_present)
           write(logunit,*) trim(subname)//":lnd_present="//trim(lnd_present)
           write(logunit,*) trim(subname)//":ocn_present="//trim(ocn_present)
           write(logunit,*) trim(subname)//":ice_present="//trim(ice_present)
           write(logunit,*) trim(subname)//":rof_present="//trim(rof_present)
           write(logunit,*) trim(subname)//":wav_present="//trim(wav_present)
           write(logunit,*) trim(subname)//":glc_present="//trim(glc_present)
           write(logunit,*) trim(subname)//":med_present="//trim(med_present)
        end if

      else

        call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
             msg=subname//' invalid component = '//trim(prefix), line=__LINE__, file=u_FILE_u, rcToReturn=rc)
        return  ! bail out

      endif

    enddo

    !--------
    ! clean-up
    !--------
    deallocate(compLabels)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

  end subroutine SetModelServices

  !================================================================================

  subroutine SetRunSequence(driver, rc)
    use ESMF                  , only : ESMF_GridComp, ESMF_LogWrite, ESMF_SUCCESS, ESMF_LOGMSG_INFO
    use ESMF                  , only : ESMF_Time, ESMF_TimeInterval, ESMF_Clock, ESMF_Config
    use ESMF                  , only : ESMF_GridCompGet
    use NUOPC                 , only : NUOPC_FreeFormat, NUOPC_FreeFormatPrint, NUOPC_FreeFormatDestroy
    use NUOPC                 , only : NUOPC_FreeFormatCreate
    use NUOPC_Driver          , only : NUOPC_DriverIngestRunSequence, NUOPC_DriverSetRunSequence
    use NUOPC_Driver          , only : NUOPC_DriverPrint

    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    integer                 :: localrc
    type(ESMF_Config)       :: config
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

    ! read free format run sequence from config
    call ESMF_GridCompGet(driver, config=config, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    runSeqFF = NUOPC_FreeFormatCreate(config, label="runSeq::", rc=rc)
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
    use perf_mod     , only : t_initf

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

    call t_initf(nlfilename, LogPrint=.true., mpicom=mpicom_GLOID, mastertask=mastertask, MaxThreads=maxthreads)

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
    type(ESMF_VM)                :: vm
    character(SHR_KIND_CL)       :: cvalue         ! temporary
    integer                      :: ierr           ! error return
    integer                      :: lmpicom        ! driver mpi communicator
    integer                      :: unitn          ! Namelist unit number to read
    character(SHR_KIND_CL)       :: restart_file   ! Full archive path to restart file
    character(SHR_KIND_CL)       :: restart_pfile  ! Restart pointer file
    character(SHR_KIND_CL)       :: rest_case_name ! Short case identification
    character(len=*) , parameter :: sp_str = 'str_undefined'
    integer :: dbrc
    character(len=*) , parameter :: subname = "(esm.F90:InitRestart)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_GridCompGet(driver, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------
    ! Carry out restart if appropriate
    !-----------------------------------------------------


    ! Add rest_case_name to driver attributes
    call NUOPC_CompAttributeAdd(driver, attrList=(/'rest_case_name'/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    rest_case_name = ' '
    call NUOPC_CompAttributeSet(driver, name='rest_case_name', value=rest_case_name, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitRestart

  !================================================================================

  subroutine InitAttributes(driver, rc)

    use shr_sys_mod      , only : shr_sys_abort
    use ESMF             , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF             , only : ESMF_Clock, ESMF_ClockGet, ESMF_Time, ESMF_TimeGet
    use ESMF             , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LogSetError, ESMF_LOGMSG_INFO
    use ESMF             , only : ESMF_RC_NOT_VALID
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
    logical                         :: iamroot_med           ! mediator masterproc
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

    ! Determine orbital params

    call seq_comm_getinfo(CPLID, iamroot=iamroot_med)
    if (trim(orb_mode) == trim(orb_variable_year)) then
       call ESMF_GridCompGet(driver, clock=clock, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_TimeGet(CurrTime, yy=year, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       orb_cyear = orb_iyear + (year - orb_iyear_align)
       call shr_orb_params(orb_cyear, orb_eccen, orb_obliq, orb_mvelp, &
                           orb_obliqr, orb_lambm0, orb_mvelpp, iamroot_med)
    else
       call shr_orb_params(orb_iyear, orb_eccen, orb_obliq, orb_mvelp, &
                           orb_obliqr, orb_lambm0, orb_mvelpp, iamroot_med)
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
    use NUOPC, only : NUOPC_CompAttributeAdd, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
    use seq_comm_mct, only : seq_comm_inst, seq_comm_name, seq_comm_suffix

    ! input/output parameters
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_GridComp) , intent(in)    :: driver
    type(ESMF_Config)   , intent(in)    :: config
    integer             , intent(in)    :: compid
    character(len=*)    , intent(in)    :: compname
    character(len=*)    , intent(in)    :: inst_suffix
    integer             , intent(inout) :: rc

    ! locals
    integer                     :: n
    character(len=SHR_KIND_CL)  :: cvalue
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
    end if
    if(len_trim(inst_suffix) > 0) then
       call NUOPC_CompAttributeAdd(gcomp, attrList=(/'inst_suffix'/), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name='inst_suffix', value=inst_suffix, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeSet(gcomp, name='inst_string', value=inst_string, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (len_trim(seq_comm_suffix(compid)) > 0) then
          call NUOPC_CompAttributeAdd(gcomp, attrList=(/'inst_suffix'/), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeSet(gcomp, name='inst_suffix', value=trim(seq_comm_suffix(compid)), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
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
    use ESMF, only: ESMF_LogWrite, ESMF_SUCCESS, ESMF_LOGMSG_INFO
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

end module ESM
