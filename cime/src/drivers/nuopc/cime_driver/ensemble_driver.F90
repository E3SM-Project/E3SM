module Ensemble_driver

  !-----------------------------------------------------------------------------
  ! Code that creates the ensemble driver layer above the esm driver.
  ! The ensmeble driver is configured to run a single clock cycle in nuopc with time step
  ! length of stop_time - start_time.  It's purpose is to instantiate NINST copies of the
  ! esm driver and its components layed out concurently across mpi tasks.
  !-----------------------------------------------------------------------------
  use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag, CL
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr

  implicit none
  private

  public  :: SetServices
  private :: SetModelServices

  character(*),parameter :: u_FILE_u = __FILE__

!================================================================================
  contains
!================================================================================

  subroutine SetServices(ensemble_driver, rc)
    use NUOPC        , only : NUOPC_CompDerive, NUOPC_CompSpecialize
    use NUOPC_Driver , only : driver_routine_SS             => SetServices
    use NUOPC_Driver , only : ensemble_label_SetModelServices => label_SetModelServices
    use ESMF         , only : ESMF_GridComp, ESMF_Config, ESMF_GridCompSet, ESMF_ConfigLoadFile
    use ESMF         , only : ESMF_ConfigCreate
    use ESMF         , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO
    type(ESMF_GridComp)  :: ensemble_driver
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Config)    :: config
    integer              :: dbrc
    character(len=*), parameter :: subname = "(ensemble_driver.F90:SetServices)"
    !---------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    ! NUOPC_Driver registers the generic methods
    call NUOPC_CompDerive(ensemble_driver, driver_routine_SS, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(ensemble_driver, specLabel=ensemble_label_SetModelServices, &
         specRoutine=SetModelServices, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create, open and set the config
    config = ESMF_ConfigCreate(rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ConfigLoadFile(config, "nuopc.runconfig", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompSet(ensemble_driver, config=config, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine SetServices

  !================================================================================

  subroutine SetModelServices(ensemble_driver, rc)
    use ESMF                  , only : ESMF_GridComp, ESMF_VM, ESMF_Config, ESMF_Clock, ESMF_VMGet
    use ESMF                  , only : ESMF_GridCompGet, ESMF_VMGet, ESMF_ConfigGetAttribute
    use ESMF                  , only : ESMF_ConfigGetLen, ESMF_RC_NOT_VALID, ESMF_LogFoundAllocError
    use ESMF                  , only : ESMF_LogSetError, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use ESMF                  , only : ESMF_GridCompSet, ESMF_SUCCESS, ESMF_METHOD_INITIALIZE, ESMF_RC_ARG_BAD
    use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
    use NUOPC_Driver          , only : NUOPC_DriverAddComp
    use esm, only : ESMSetServices => SetServices, ReadAttributes
    use shr_nuopc_time_mod    , only : shr_nuopc_time_clockInit
    use med_internalstate_mod , only : logunit  ! initialized here
    use shr_log_mod           , only : shrloglev=>shr_log_level, shrlogunit=> shr_log_unit
    use shr_file_mod          , only : shr_file_getUnit, shr_file_getLoglevel
    use shr_file_mod          , only : shr_file_setloglevel, shr_file_setlogunit

    type(ESMF_GridComp)    :: ensemble_driver
    integer, intent(out)   :: rc

    ! local variables
    type(ESMF_VM)          :: vm
    type(ESMF_GridComp)    :: driver, gridcomptmp
    type(ESMF_Config)      :: config
    integer                :: n, n1, stat
    integer, pointer       :: petList(:)
    character(len=20)      :: model, prefix
    integer                :: petCount, i
    integer                :: localPet
    integer                :: rootpe_med
    logical                :: is_set
    character(len=512)     :: diro
    character(len=512)     :: logfile
    integer                :: global_comm
    logical                :: iamroot_med ! mediator masterproc
    logical                :: read_restart
    integer                :: dbrc
    integer                :: inst
    integer                :: number_of_members
    integer                :: ntasks_per_member
    character(CL)          :: start_type     ! Type of startup
    character(len=7)       :: drvrinst
    character(len=5)       :: inst_suffix
    character(len=CL)      :: msgstr
    character(len=CL)      :: cvalue
    character(len=*) , parameter :: start_type_start = "startup"
    character(len=*) , parameter :: start_type_cont  = "continue"
    character(len=*) , parameter :: start_type_brnch = "branch"
    character(len=*) , parameter :: subname = "(ensemble_driver.F90:SetModelServices)"

    !-------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_GridCompGet(ensemble_driver, config=config, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=localPet, PetCount=PetCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------
    ! Initialize clocks
    !-------------------------------------------
    call ReadAttributes(ensemble_driver, config, "ALLCOMP_attributes::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(ensemble_driver, config, "CLOCK_attributes::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(ensemble_driver, config, "PELAYOUT_attributes::", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Check valid values of start type
    call NUOPC_CompAttributeGet(ensemble_driver, name="start_type", value=start_type, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if ((trim(start_type) /= start_type_start) .and.  &
        (trim(start_type) /= start_type_cont ) .and.  &
        (trim(start_type) /= start_type_brnch)) then
       write (msgstr, *) subname//': start_type invalid = '//trim(start_type)
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
    end if

    call InitRestart(ensemble_driver, read_restart, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Need to get number of ensemble members here
    call NUOPC_CompAttributeGet(ensemble_driver, name="ninst", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue, *) number_of_members
    !-------------------------------------------
    ! Extract the config object from the ensemble_driver
    !-------------------------------------------
    ntasks_per_member = PetCount/number_of_members
    if(ntasks_per_member*number_of_members /= PetCount) then
       write (msgstr, *) "PetCount must be evenly divisable by number of members "
       call ESMF_LogSetError(ESMF_RC_ARG_BAD, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
    endif

    allocate(petList(ntasks_per_member))

    call NUOPC_CompAttributeGet(ensemble_driver, name='cpl_rootpe', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue, *) rootpe_med

    do inst=1,number_of_members

       petList(1) = (inst-1) * ntasks_per_member
       do n=2,ntasks_per_member
          petList(n) = petList(n-1) + 1
       enddo
       write(drvrinst,'(a,i4.4)') "ESM",inst
       call NUOPC_DriverAddComp(ensemble_driver, drvrinst, ESMSetServices, petList=petList, comp=gridcomptmp, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (localpet >= petlist(1) .and. localpet <= petlist(ntasks_per_member)) then
          driver = gridcomptmp
          if(number_of_members > 1) then
             call NUOPC_CompAttributeAdd(driver, attrList=(/'inst_suffix'/), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             write(inst_suffix,'(a,i4.4)') '_',inst
             call NUOPC_CompAttributeSet(driver, name='inst_suffix', value=inst_suffix, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             inst_suffix = ''
          endif
          write(cvalue,*) read_restart
          call NUOPC_CompAttributeAdd(driver, attrList=(/'read_restart'/), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeSet(driver, name='read_restart', value=trim(cvalue), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call ReadAttributes(driver, config, "MED_attributes::", rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call ReadAttributes(driver, config, "CLOCK_attributes::", rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call ReadAttributes(driver, config, "MED_modelio"//trim(inst_suffix)//"::", rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          if (mod(localPet, ntasks_per_member) == rootpe_med) then
             call NUOPC_CompAttributeGet(driver, name="diro", value=diro, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             call NUOPC_CompAttributeGet(driver, name="logfile", value=logfile, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             logunit = shr_file_getUnit()
             open(logunit,file=trim(diro)//"/"//trim(logfile))
          else
             logUnit = shrlogunit
          endif
          call shr_file_getLogLevel(shrloglev)
          call shr_file_setLogLevel(max(shrloglev,1))
          call shr_file_setLogUnit (logunit)
       endif
    enddo
    call shr_nuopc_time_clockInit(ensemble_driver, driver, logunit, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(petList)

  end subroutine SetModelServices

  subroutine InitRestart(ensemble_driver, read_restart, rc)

    !-----------------------------------------------------
    ! Determine if will restart and read pointer file
    ! if appropriate
    !-----------------------------------------------------
    use ESMF         , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_SUCCESS
    use ESMF         , only : ESMF_LogSetError, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_RC_NOT_VALID
    use NUOPC        , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
    ! input/output variables
    type(ESMF_GridComp)    , intent(inout) :: ensemble_driver
    logical                , intent(out)   :: read_restart   ! read the restart file, based on start_type
    integer                , intent(out)   :: rc

    ! local variables
    character(len=CL)       :: cvalue         ! temporary
    integer                 :: ierr           ! error return
    integer                 :: unitn          ! Namelist unit number to read

    character(len=CL)       :: restart_file   ! Full archive path to restart file
    character(len=CL)       :: restart_pfile  ! Restart pointer file
    character(len=CL)       :: rest_case_name ! Short case identification
    character(len=CL)       :: start_type     ! Type of startup
    character(len=CL)       :: msgstr
    character(len=*) , parameter :: start_type_start = "startup"
    character(len=*) , parameter :: start_type_cont  = "continue"
    character(len=*) , parameter :: start_type_brnch = "branch"
    character(len=*) , parameter :: sp_str = 'str_undefined'
    integer :: dbrc
    character(len=*) , parameter :: subname = "(esm.F90:InitRestart)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    !-----------------------------------------------------
    ! Carry out restart if appropriate
    !-----------------------------------------------------

    ! First Determine if restart is read
    call NUOPC_CompAttributeGet(ensemble_driver, name='start_type', value=start_type, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Check valid values of start type

    if ((trim(start_type) /= start_type_start) .and.  &
        (trim(start_type) /= start_type_cont ) .and.  &
        (trim(start_type) /= start_type_brnch)) then
       write (msgstr, *) subname//': start_type invalid = '//trim(start_type)
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
    end if

    !TODO: this is hard-wired to CIME start/continue types in terms of gcomp
    read_restart = .false.
    if (trim(start_type) == trim(start_type_cont) .or. trim(start_type) == trim(start_type_brnch)) then
       read_restart = .true.
    endif

    ! Add rest_case_name and read_restart to ensemble_driver attributes
    call NUOPC_CompAttributeAdd(ensemble_driver, attrList=(/'rest_case_name','read_restart'/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    rest_case_name = ' '
    call NUOPC_CompAttributeSet(ensemble_driver, name='rest_case_name', value=rest_case_name, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    write(cvalue,*) read_restart
    call NUOPC_CompAttributeSet(ensemble_driver, name='read_restart', value=trim(cvalue), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitRestart

end module ENSEMBLE_DRIVER
