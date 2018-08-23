module Ensemble_driver

  !-----------------------------------------------------------------------------
  ! Code that creates the ensemble driver layer above the esm driver.
  ! The ensmeble driver is configured to run a single clock cycle in nuopc with time step
  ! length of stop_time - start_time.  It's purpose is to instantiate NINST copies of the
  ! esm driver and its components layed out concurently across mpi tasks.  
  !-----------------------------------------------------------------------------
  use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr

  implicit none
  private

  public  :: SetServices
  private :: SetModelServices
  logical :: mastertask

  character(*),parameter :: u_FILE_u = __FILE__

!================================================================================
  contains
!================================================================================

  subroutine SetServices(ensemble_driver, rc)
    use NUOPC        , only : NUOPC_CompDerive, NUOPC_CompSpecialize
    use NUOPC_Driver , only : driver_routine_SS             => SetServices
    use NUOPC_Driver , only : driver_label_SetModelServices => label_SetModelServices
    use ESMF         , only : ESMF_GridComp, ESMF_Config, ESMF_GridCompSet, ESMF_ConfigLoadFile
    use ESMF         , only : ESMF_ConfigCreate
    use ESMF         , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO

    type(ESMF_GridComp)  :: ensemble_driver
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Config) :: config
    integer           :: dbrc
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
    call NUOPC_CompSpecialize(ensemble_driver, specLabel=driver_label_SetModelServices, &
         specRoutine=SetModelServices, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create, open and set the config
    config = ESMF_ConfigCreate(rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ConfigLoadFile(config, "cesm.runconfig", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompSet(ensemble_driver, config=config, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine SetServices

  !================================================================================

  subroutine SetModelServices(ensemble_driver, rc)
    use ESMF                  , only : ESMF_GridComp, ESMF_VM, ESMF_Config, ESMF_Clock
    use ESMF                  , only : ESMF_GridCompGet, ESMF_VMGet, ESMF_ConfigGetAttribute
    use ESMF                  , only : ESMF_ConfigGetLen, ESMF_RC_NOT_VALID, ESMF_LogFoundAllocError
    use ESMF                  , only : ESMF_LogSetError, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use ESMF                  , only : ESMF_GridCompSet, ESMF_SUCCESS, ESMF_METHOD_INITIALIZE, ESMF_RC_ARG_BAD
    use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
    use NUOPC_Driver          , only : NUOPC_DriverAddComp

    use med_constants_mod, only : CL
    use esm, only : ESMSetServices => SetServices, ReadAttributes
    use esm, only : EClock_d, EClock_a, EClock_l, EClock_o, EClock_i, EClock_g, EClock_r, Eclock_w, EClock_e

    type(ESMF_GridComp)    :: ensemble_driver
    type(ESMF_Clock) :: Eclock_ensemble
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
    logical                :: is_set
    character(len=512)     :: diro
    character(len=512)     :: logfile
    integer                :: shrlogunit  ! original log unit
    integer                :: shrloglev   ! original log level
    integer                :: global_comm
    logical                :: iamroot_med ! mediator masterproc
    integer                :: dbrc
    integer :: inst
    integer :: number_of_members
    integer :: ntasks_per_member
    logical :: read_restart
    character(CL)       :: start_type     ! Type of startup
    character(len=*) , parameter :: start_type_start = "startup"
    character(len=*) , parameter :: start_type_cont  = "continue"
    character(len=*) , parameter :: start_type_brnch = "branch"


    character(CL) :: msgstr, cvalue
    character(len=7) :: drvrinst
    character(len=4) :: inst_string

    character(len=*), parameter    :: subname = "(ensemble_driver.F90:SetModelServices)"

    !-------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_GridCompGet(ensemble_driver, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=localPet, PetCount=PetCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (localPet == 0) then
       mastertask=.true.
    else
       mastertask = .false.
    end if
    call ESMF_GridCompGet(ensemble_driver, config=config, vm=vm, rc=rc)
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

    !TODO: this is hard-wired to CIME start/continue types in terms of gcomp
    read_restart = .false.
    if (trim(start_type) == trim(start_type_cont) .or. trim(start_type) == trim(start_type_brnch)) then
       read_restart = .true.
    endif
    ! Add read_restart to driver attributes
    call NUOPC_CompAttributeAdd(ensemble_driver, attrList=(/'read_restart'/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    write(cvalue,*) read_restart
    call NUOPC_CompAttributeSet(ensemble_driver, name='read_restart', value=trim(cvalue), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

! Need to get number of ensemble members here

    ! Check valid values of start type
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

    write(cvalue,*) read_restart

    do inst=1,number_of_members

       petList(1) = (inst-1) * ntasks_per_member
       do n=2,ntasks_per_member
          petList(n) = petList(n-1) + 1
       enddo
       write(drvrinst,'(a,i4.4)') "ESM",inst
       call NUOPC_DriverAddComp(ensemble_driver, drvrinst, ESMSetServices, petList=petList, comp=gridcomptmp, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if(number_of_members == 1) then
          inst_string = "_"
       else
          write(inst_string,'(i4.4)') inst
       endif
       if (localpet >= petlist(1) .and. localpet <= petlist(ntasks_per_member)) then
          driver = gridcomptmp
          call NUOPC_CompAttributeAdd(driver, attrList=(/'INST'/), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeSet(driver, name='INST', value=inst_string, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          
          call NUOPC_CompAttributeAdd(driver, attrList=(/'read_restart'/), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeSet(driver, name='read_restart', value=trim(cvalue), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
    enddo
    deallocate(petList)



    call InitClocks(ensemble_driver, EClock_ensemble, Eclock_d, Eclock_a, Eclock_l, Eclock_o, &
         Eclock_i, Eclock_g, Eclock_r, Eclock_w, Eclock_e, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------
    ! Set baseline clock
    !--------
    call ESMF_GridCompSet(ensemble_driver, clock=Eclock_ensemble, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine SetModelServices


  subroutine InitClocks(driver, EClock_ensemble, Eclock_d, Eclock_a, Eclock_l, Eclock_o, &
       Eclock_i, Eclock_g, Eclock_r, Eclock_w, Eclock_e, rc)

    !----------------------------------------------------------
    ! Initialize time manager
    !----------------------------------------------------------
    use ESMF            , only : ESMF_GridComp, ESMF_Clock, ESMF_VM, ESMF_GridCompGet, ESMF_VMGet
    use ESMF            , only : ESMF_LogWrite, ESMF_SUCCESS, ESMF_LOGMSG_INFO
    use pio             , only : pio_file_is_open, pio_closefile, file_desc_t
    use seq_timemgr_mod , only : seq_timemgr_clockInit

    ! INPUT/OUTPUT PARAMETERS:
    type(ESMF_GridComp)    , intent(inout) :: driver
    type(ESMF_Clock)       , intent(inout) :: EClock_ensemble
    type(ESMF_Clock)       , intent(inout) :: EClock_d
    type(ESMF_Clock)       , intent(inout) :: EClock_a
    type(ESMF_Clock)       , intent(inout) :: EClock_l
    type(ESMF_Clock)       , intent(inout) :: EClock_o
    type(ESMF_Clock)       , intent(inout) :: EClock_i
    type(ESMF_Clock)       , intent(inout) :: EClock_g
    type(ESMF_Clock)       , intent(inout) :: EClock_r
    type(ESMF_Clock)       , intent(inout) :: EClock_w
    type(ESMF_Clock)       , intent(inout) :: EClock_e
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_VM)     :: vm
    type(file_desc_t) :: pioid
    integer           :: lmpicom
    integer :: dbrc
    character(len=*) , parameter    :: subname = '(InitClocks)'
    !----------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    !----------------------------------------------------------
    ! Initialize time manager
    !----------------------------------------------------------

    call ESMF_GridCompGet(driver, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call seq_timemgr_clockInit(driver, mastertask, pioid, lmpicom, &
         Eclock_ensemble, EClock_d, EClock_a, EClock_l, EClock_o, &
         EClock_i, Eclock_g, Eclock_r, Eclock_w, Eclock_e)

    if (pio_file_is_open(pioid)) then
       call pio_closefile(pioid)
    endif

  end subroutine InitClocks

end module ENSEMBLE_DRIVER
