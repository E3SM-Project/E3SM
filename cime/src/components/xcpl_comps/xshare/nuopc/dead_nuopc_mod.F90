module dead_nuopc_mod

  use ESMF              , only : ESMF_Gridcomp, ESMF_State, ESMF_StateGet
  use ESMF              , only : ESMF_Clock, ESMF_Time, ESMF_TimeInterval, ESMF_Alarm
  use ESMF              , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_ClockSet, ESMF_ClockAdvance, ESMF_AlarmSet
  use ESMF              , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_METHOD_INITIALIZE
  use ESMF              , only : ESMF_FAILURE, ESMF_LOGMSG_ERROR
  use ESMF              , only : ESMF_VMGetCurrent, ESMF_VM, ESMF_VMBroadcast, ESMF_VMGet
  use ESMF              , only : ESMF_VM, ESMF_VMGetCurrent, ESMF_VmGet
  use ESMF              , only : operator(/=), operator(==), operator(+)
  use shr_kind_mod      , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod       , only : shr_sys_abort
  use dead_methods_mod  , only : chkerr, alarmInit

  implicit none
  private

  public :: dead_init_nuopc
  public :: dead_final_nuopc
  public :: dead_meshinit
  public :: ModelInitPhase
  public :: ModelSetRunClock
  public :: fld_list_add
  public :: fld_list_realize

  ! !PUBLIC DATA MEMBERS:
  integer, public :: dead_grid_lat    = 1   ! lat from component
  integer, public :: dead_grid_lon    = 2   ! lon from component
  integer, public :: dead_grid_area   = 3   ! area from component
  integer, public :: dead_grid_mask   = 4   ! mask, 0 = inactive cell
  integer, public :: dead_grid_frac   = 5   ! fractional area coverage
  integer, public :: dead_grid_index  = 6   ! global index
  integer, public :: dead_grid_total  = 6

  type fld_list_type
    character(len=128) :: stdname
     integer :: ungridded_lbound = 0
     integer :: ungridded_ubound = 0
  end type fld_list_type
  public :: fld_list_type

  integer, parameter, public :: fldsMax = 100
  integer                    :: dbug_flag = 0
  character(*), parameter    :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dead_read_inparms(model,  inst_suffix, logunit, &
       nxg, nyg, decomp_type, nproc_x, seg_len)

    ! input/output variables
    character(len=*) , intent(in)    :: model
    character(len=*) , intent(in)    :: inst_suffix ! char string associated with instance
    integer          , intent(in)    :: logunit     ! logging unit number
    integer    	     , intent(out)   :: nproc_x
    integer          , intent(out)   :: seg_len
    integer          , intent(out)   :: nxg         ! global dim i-direction
    integer          , intent(out)   :: nyg         ! global dim j-direction
    integer          , intent(out)   :: decomp_type ! decomposition type

    ! local variables
    type(ESMF_VM)           :: vm
    character(CL)           :: fileName ! generic file name
    integer                 :: nunit    ! unit number
    integer                 :: unitn    ! Unit for namelist file
    integer                 :: tmp(5)   ! array for broadcast
    integer                 :: localPet ! mpi id of current task in current context
    integer                 :: rc       ! return code
    character(*), parameter :: F00   = "('(dead_read_inparms) ',8a)"
    character(*), parameter :: F01   = "('(dead_read_inparms) ',a,a,4i8)"
    character(*), parameter :: F03   = "('(dead_read_inparms) ',a,a,i8,a)"
    character(*), parameter :: subName = "(dead_read_inpamrs) "
    !-------------------------------------------------------------------------------

    ! read the input parms (used to configure model)
    nxg            =  -9999
    nyg            =  -9999
    nproc_x        =  -9999
    seg_len        =  -9999
    decomp_type    =  -9999

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (localPet==0) then
       open(newunit=unitn, file='x'//model//'_in'//trim(inst_suffix), status='old' )
       read(unitn,*) nxg
       read(unitn,*) nyg
       read(unitn,*) decomp_type
       read(unitn,*) nproc_x
       read(unitn,*) seg_len
       close (unitn)
    endif

    tmp(1) = nxg
    tmp(2) = nyg
    tmp(3) = decomp_type
    tmp(4) = nproc_x
    tmp(5) = seg_len

    call ESMF_VMBroadcast(vm, tmp, 6, 0, rc=rc)

    nxg         = tmp(1)
    nyg         = tmp(2)
    decomp_type = tmp(3)
    nproc_x     = tmp(4)
    seg_len     = tmp(5)

    if (localPet==0) then
       write(logunit,*)' Read in X'//model//' input from file= x'//model//'_in'
       write(logunit,F00) model
       write(logunit,F00) model,'         Model  :  ',model
       write(logunit,F01) model,'           NGX  :  ',nxg
       write(logunit,F01) model,'           NGY  :  ',nyg
       write(logunit,F01) model,' Decomposition  :  ',decomp_type
       write(logunit,F03) model,' Num pes in X   :  ',nproc_x,'  (type 3 only)'
       write(logunit,F03) model,' Segment Length :  ',seg_len,'  (type 11 only)'
       write(logunit,F00) model,'    inst_suffix :  ',trim(inst_suffix)
       write(logunit,F00) model
    end if

  end subroutine dead_read_inparms

  !===============================================================================

  subroutine dead_setNewGrid(decomp_type, nxg, nyg, logunit, lsize, gbuf, seg_len, nproc_x)

    ! This sets up some defaults.  The user may want to overwrite some
    ! of these fields in the main program after initialization in complete.

    use shr_const_mod , only : shr_const_pi, shr_const_rearth

    ! input/output parameters:
    integer , intent(in)          :: decomp_type !
    integer , intent(in)          :: nxg,nyg     ! global grid sizes
    integer , intent(in)          :: logunit     ! output logunit
    integer , intent(out)         :: lsize       ! local grid sizes
    real(R8), pointer             :: gbuf(:,:)   ! output data
    integer , intent(in),optional :: seg_len     ! seg len decomp setting
    integer , intent(in),optional :: nproc_x     ! 2d decomp setting

    ! local
    type(ESMF_VM)           :: vm
    integer                 :: rc
    integer                 :: mype
    integer                 :: totpe ! total number of pes
    logical                 :: found
    integer                 :: i,j,ig,jg
    integer                 :: n,ng,is,ie,js,je,nx,ny
    integer                 :: npesx,npesy,mypex,mypey,nxp,nyp
    real(R8)                :: hscore,bscore
    real(R8)                :: dx,dy,deg2rad,ys,yc,yn,area,re
    integer, allocatable    :: gindex(:)
    character(*), parameter :: F00   = "('(dead_setNewGrid) ',8a)"
    character(*), parameter :: F01   = "('(dead_setNewGrid) ',a,4i8)"
    character(*), parameter :: subName = "(dead_setNewGrid) "
    !-------------------------------------------------------------------------------

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=mype, peCount=totpe, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if ( decomp_type == 1 .or. &
         decomp_type == 2 .or. &
         decomp_type == 3 .or. &
         decomp_type == 4 .or. &
         decomp_type == 11) then
    else
       ! invalid decomposition type
       if (mype == 0) then
          write(logunit,F01) 'ERROR: invalid decomp_type = ',decomp_type
       end if
       call shr_sys_abort(subName//'invalid decomp_type')
    endif

    if (nxg*nyg == 0) then
       lsize = 0
       allocate(gbuf(lsize,dead_grid_total))
       if (mype == 0) then
          write(logunit,*) subname,' grid size is zero, lsize = ',lsize
       end if
       return
    endif

    found = .false.

    if (decomp_type == 1) then  ! 1d decomp by lat
       npesx = 1
       npesy = totpe
       found = .true.
    elseif (decomp_type == 2) then  ! 1d decomp by lon
       npesx = totpe
       npesy = 1
       found = .true.
    elseif (decomp_type == 3) then  ! 2d decomp
       if (present(nproc_x)) then
          if ( nproc_x > 0 ) then
             npesx=nproc_x
             npesy=totpe/npesx
             if ( npesx*npesy /= totpe) then
                write(logunit,F00) 'ERROR: uneven decomposition'
                call shr_sys_abort(subName//'uneven decomp')
             end if
             found = .true.
          endif
       endif
       if (.not.found) then  ! narrow blocks
          do nx = 1,totpe
             ny = totpe/nx
             if (nx*ny == totpe) then
                npesx = nx
                npesy = ny
                found = .true.
             endif
          enddo
       endif
    elseif (decomp_type == 4) then  ! 2d evenly divisible square block decomp
       hscore = nxg*nyg
       do nx = 1,totpe
          ny = totpe/nx
          if (nx*ny == totpe .and. mod(nxg,nx) == 0 .and. mod(nyg,ny) == 0) then
             bscore = ((nxg*ny*1.0_r8) / (nyg*nx*1.0_r8)) - 1.0_r8
             bscore = bscore * bscore
             if (bscore < hscore .or. .not.found) then
                hscore = bscore
                npesx = nx
                npesy = ny
                found = .true.
             endif
          endif
       enddo
    endif

    if (found) then
       nx = nxg/npesx
       mypex = mod(mype,npesx)
       mypey = mype/npesx
       is = (mypex    ) * (nx) + 1
       ie = (mypex + 1) * (nx)

       ny = nyg/npesy
       js = (mypey    ) * (ny) + 1
       je = (mypey + 1) * (ny)

       nxp = nxg - (nx*npesx)       ! extra lons not accounted for yet
       nyp = nyg - (ny*npesy)       ! extra lats not accounted for yet

       is = is + min(mypex,nxp)      ! add them to first few pes and shift everything
       ie = ie + min(mypex+1,nxp)
       js = js + min(mypey,nyp)      ! add them to first few pes and shift everything
       je = je + min(mypey+1,nyp)

       lsize = (ie - is + 1) * (je - js + 1)

       allocate(gindex(lsize))
       n = 0
       do j = js,je
          do i = is,ie
             n = n + 1
             gindex(n) = (j-1)*nxg + i
          enddo
       enddo
    endif

    if (.not.found) then
       !-------------------------------------------------------------------------
       ! type 11 general segment decomp
       !-------------------------------------------------------------------------
       nx = nxg*nyg / (totpe*13) + 1   ! 13 segments per pe (arbitrary)
       ! nx override with seg_len
       if (present(seg_len)) then
          if (seg_len > 0) nx = seg_len
       endif

       n = 0
       i = 0
       lsize = 0
       do while (n < nxg*nyg)
          ny = min(nx,nxg*nyg-n)
          do j = 1,ny
             n = n + 1
             if (mype == mod(i,totpe)) then
                lsize = lsize + 1
             endif
          enddo
          i = i + 1
       enddo

       allocate(gindex(lsize))

       n = 0
       i = 0
       lsize = 0
       do while (n < nxg*nyg)
          ny = min(nx,nxg*nyg-n)
          do j = 1,ny
             n = n + 1
             if (mype == mod(i,totpe)) then
                lsize = lsize + 1
                gindex(lsize) = n
             endif
          enddo
          i = i + 1
       enddo

       if (mype == 0) then
          write(logunit,*) 'dead_setNewGrid decomp seg ',mype,lsize,nx
       end if

       found = .true.

    endif

    if ( .not.found ) then
       write(logunit,F01) 'ERROR: with decomp nxg,nyg,totpe=',nxg,nyg,totpe
       call shr_sys_abort(subName//'decomp')
    endif

    deg2rad = shr_const_pi / 180.0_R8
    re = shr_const_rearth

    allocate(gbuf(lsize,dead_grid_total))
    gbuf = -888.0_R8
    if (mype == 0) then
       write(logunit,*) subname,' Decomp is ',decomp_type,' lsize = ',lsize
    end if

    n=0
    dx = 360.0_R8/nxg * deg2rad
    do n = 1,lsize
       ig = mod((gindex(n)-1),nxg) + 1
       jg =     (gindex(n)-1)/nxg  + 1

       ys = -90.0_R8 + (jg-1.0_R8)*180.0_R8/(nyg)
       yc = -90.0_R8 + (jg-0.5_R8)*180.0_R8/(nyg)
       yn = -90.0_R8 + (jg-0.0_R8)*180.0_R8/(nyg)
       dy = sin(yn*deg2rad) - sin(ys*deg2rad)
       area = dx*dy*re*re

       gbuf(n,dead_grid_lon  ) = (ig-1.0_R8)*360.0_R8/(nxg)
       gbuf(n,dead_grid_lat  ) = yc
       gbuf(n,dead_grid_index) = gindex(n)
       gbuf(n,dead_grid_area ) = area
       gbuf(n,dead_grid_mask ) = 0
       gbuf(n,dead_grid_frac ) = 1.0_R8
    enddo

    deallocate(gindex)

  end subroutine dead_setNewGrid

  !===============================================================================

  subroutine dead_init_nuopc(model, inst_suffix, logunit, lsize, gbuf, nxg, nyg)

    ! input/output parameters:
    character(len=*) , intent(in)    :: model
    character(len=*) , intent(in)    :: inst_suffix ! char string associated with instance
    integer          , intent(in)    :: logunit     ! logging unit number
    integer          , intent(out)   :: lsize       ! logging unit number
    real(r8)         , pointer       :: gbuf(:,:)   ! model grid
    integer          , pointer       :: gindex(:)   ! global index space
    integer          , intent(out)   :: nxg         ! global dim i-direction
    integer          , intent(out)   :: nyg         ! global dim j-direction

    !--- local variables ---
    integer                  :: local_comm    ! local communicator
    integer                  :: nproc_x
    integer                  :: seg_len
    integer                  :: decomp_type
    character(*), parameter  :: subName = "(dead_init_nuopc) "
    !-------------------------------------------------------------------------------

    ! Read input parms
    call dead_read_inparms(model, inst_suffix, logunit, nxg, nyg, decomp_type, nproc_x, seg_len)

    ! Initialize grid
    call dead_setNewGrid(decomp_type, nxg, nyg, logunit, lsize, gbuf, seg_len, nproc_x)

  end subroutine dead_init_nuopc

  !===============================================================================

  subroutine dead_final_nuopc(model, logunit)

    ! finalize method for xcpl component

    ! input/output parameters:
    character(len=*) , intent(in) :: model
    integer          , intent(in) :: logunit     ! logging unit number

    ! local variables
    type(ESMF_VM) :: vm
    integer       :: rc
    integer       :: localPet
    character(*), parameter :: F00   = "('(dead_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(dead_comp_final) ',73('-'))"
    character(*), parameter :: subName = "(dead_comp_final) "
    !-------------------------------------------------------------------------------

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (localPet==0) then
       write(logunit,F91)
       write(logunit,F00) trim(model),': end of main integration loop'
       write(logunit,F91)
    end if

  end subroutine dead_final_nuopc

  !===============================================================================

  subroutine fld_list_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer                    , intent(inout) :: num
    type(fld_list_type)        , intent(inout) :: fldlist(:)
    character(len=*)           , intent(in)    :: stdname
    integer,          optional , intent(in)    :: ungridded_lbound
    integer,          optional , intent(in)    :: ungridded_ubound

    ! local variables
    character(len=*), parameter :: subname='(dead_nuopc_mod:fld_list_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information
    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
       return
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fld_list_add

  !===============================================================================

  subroutine fld_list_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer           :: n
    type(ESMF_Field)  :: field
    character(len=80) :: stdname
    integer           :: gridtoFieldMap=2
    character(len=*),parameter  :: subname='(dead_nuopc_mod:fld_list_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO)
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                     ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                     ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                     gridToFieldMap=(/gridToFieldMap/), rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
             end if
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------

      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(dead_nuopc_mod:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fld_list_realize

  !===============================================================================

  subroutine ModelInitPhase(gcomp, importState, exportState, clock, rc)

    use NUOPC, only : NUOPC_CompFilterPhaseMap

    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine ModelInitPhase

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    use ESMF               , only : ESMF_ClockGetAlarmList, ESMF_ALARMLIST_ALL
    use NUOPC_Model        , only : NUOPC_ModelGet
    use NUOPC              , only : NUOPC_CompAttributeGet

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option       ! Restart option units
    integer                  :: restart_n            ! Number until restart interval
    integer                  :: restart_ymd          ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname='dead_nuopc_mod:(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart alarm
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================
  
  subroutine dead_meshinit(gcomp, nx_global, ny_global, gindex, lon, lat, Emesh, rc)

    !-----------------------------------------
    ! create an Emesh object for Fields
    !-----------------------------------------

    use ESMF , only : ESMF_GridComp, ESMF_VM, ESMF_Mesh
    use ESMF , only : ESMF_VMGet, ESMF_GridCompGet, ESMF_VMBroadCast, ESMF_VMAllGatherV
    use ESMF , only : ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_LogWrite
    use ESMF , only : ESMF_VMGather, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
    use ESMF , only : ESMF_MeshCreate, ESMF_COORDSYS_SPH_DEG, ESMF_REDUCE_SUM
    use ESMF , only : ESMF_VMAllReduce, ESMF_MESHELEMTYPE_QUAD

    ! input/output arguments
    type(ESMF_GridComp)               :: gcomp
    integer           , intent(in)    :: nx_global
    integer           , intent(in)    :: ny_global
    integer           , intent(in)    :: gindex(:)
    real(r8), pointer , intent(in)    :: lon(:)
    real(r8), pointer , intent(in)    :: lat(:)
    type(ESMF_Mesh)   , intent(inout) :: Emesh
    integer           , intent(inout) :: rc

    ! local variables
    integer          :: n,n1,n2,de
    integer          :: iam
    integer          :: lsize
    integer          :: numTotElems, numNodes, numConn, nodeindx
    integer          :: iur,iul,ill,ilr
    integer          :: xid, yid, xid0, yid0
    real(r8)         :: lonur, lonul, lonll, lonlr
    integer, pointer :: iurpts(:)
    integer, pointer :: elemIds(:)
    integer, pointer :: elemTypes(:)
    integer, pointer :: elemConn(:)
    real(r8),pointer :: elemCoords(:)
    integer, pointer :: nodeIds(:)
    integer, pointer :: nodeOwners(:)
    real(r8),pointer :: nodeCoords(:)
    real(r8),pointer :: latG(:)
    real(r8),pointer :: lonG(:)
    integer ,pointer :: pes_local(:)
    integer ,pointer :: pes_global(:)
    integer, pointer :: recvOffsets(:)
    integer, pointer :: recvCounts(:)
    integer          :: sendData(1)
    type(ESMF_VM)    :: vm
    integer          :: petCount
    character(len=*),parameter :: subname='(dead_MeshInit)'
    !--------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname, ESMF_LOGMSG_INFO, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lsize = size(gindex)

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, petCount=petCount, localpet=iam, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    allocate(latG(nx_global*ny_global))
    allocate(lonG(nx_global*ny_global))

    allocate(recvoffsets(petCount))
    allocate(recvCounts(petCount))

    sendData(1) = lsize
    call ESMF_VMGather(vm, sendData=sendData, recvData=recvCounts, count=1, rootPet=0, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMBroadCast(vm, bcstData=recvCounts, count=petCount, rootPet=0, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    recvoffsets(1) = 0
    do n = 2,petCount
       recvoffsets(n) = recvoffsets(n-1) + recvCounts(n-1)
    end do

    call ESMF_VMAllGatherV(vm, lat, lsize, latG, recvCounts, recvOffsets, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMAllGatherV(vm, lon, lsize, lonG, recvCounts, recvOffsets, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    deallocate(recvoffsets)
    deallocate(recvCounts)

    ! assumes quadrilaterals for each gridcell (element)
    ! element index matches gsmap index value
    ! nodeid at lower left of each gridcell matches gsmap index value
    ! assumes wrap around in x direction but no wrap in y direction
    ! node ids need to be described in counter clockwise direction
    ! node id associated with lower left cell is assigned to local PET
    ! node ids at top of y boundary assigned to the element to the right

    numTotElems = lsize

    allocate(elemIds(numTotElems))
    allocate(elemTypes(numTotElems))
    elemTypes=(/ESMF_MESHELEMTYPE_QUAD/)
    allocate(elemConn(4*numTotElems))
    allocate(elemCoords(2*numTotElems))

    allocate(nodeIds(numTotElems*4))
    nodeIds = -99

    elemIds(:) = gindex(:)
    numNodes = 0
    numConn = 0

    do n = 1,numTotElems
       elemTypes(n) = ESMF_MESHELEMTYPE_QUAD
       elemCoords(2*n-1) = lon(n)
       elemCoords(2*n)   = lat(n)

       do n1 = 1,4

          numNodes = numNodes + 1
          nodeindx = numNodes
          if (n1 == 1 .or. n1 == 3) xid = mod(elemIds(n)-1,nx_global) + 1
          if (n1 == 2 .or. n1 == 4) xid = mod(elemIds(n)  ,nx_global) + 1
          if (n1 == 1 .or. n1 == 2) yid = (elemIds(n)-1)/nx_global + 1
          if (n1 == 3 .or. n1 == 4) yid = (elemIds(n)-1)/nx_global + 2
          nodeIds(numNodes) = (yid-1) * nx_global + xid
          n2 = 0
          do while (n2 < numNodes - 1 .and. nodeindx == numNodes)
             n2 = n2 + 1
             if (nodeIds(numNodes) == nodeIds(n2)) nodeindx = n2
          enddo
          if (nodeindx /= numNodes) then
             numNodes = numNodes - 1
          endif

          numConn = numConn + 1
          elemConn(numConn) = nodeindx
       enddo
    enddo


    allocate(nodeCoords(2*numNodes))
    allocate(nodeOwners(numNodes))
    allocate(iurpts(numNodes))

    do n = 1,numNodes

       xid0 = mod(nodeIds(n)-1, nx_global) + 1
       yid0 = (nodeIds(n)-1) / nx_global + 1

       xid = xid0
       yid = max(min(yid0,ny_global),1)
       iur = (yid-1) * nx_global + xid
       iurpts(n) = iur

       xid = mod(xid0 - 2 + nx_global, nx_global) + 1
       yid = max(min(yid0,ny_global),1)
       iul = (yid-1) * nx_global + xid

       xid = mod(xid0 - 2 + nx_global, nx_global) + 1
       yid = max(min(yid0-1,ny_global),1)
       ill = (yid-1) * nx_global + xid

       xid = xid0
       yid = max(min(yid0-1,ny_global),1)
       ilr = (yid-1) * nx_global + xid

       ! write(tmpstr,'(2a,8i6)') subname,' nodecoord = ',n,nodeIds(n),xid0,yid0,iur,iul,ill,ilr
       ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

       ! need to normalize lon values to same 360 degree setting, use lonur as reference value
       lonur = lonG(iur)
       lonul = lonG(iul)
       lonll = lonG(ill)
       lonlr = lonG(ilr)

       if (abs(lonul + 360._r8 - lonur) < abs(lonul - lonur)) lonul = lonul + 360._r8
       if (abs(lonul - 360._r8 - lonur) < abs(lonul - lonur)) lonul = lonul - 360._r8
       if (abs(lonll + 360._r8 - lonur) < abs(lonll - lonur)) lonll = lonll + 360._r8
       if (abs(lonll - 360._r8 - lonur) < abs(lonll - lonur)) lonll = lonll - 360._r8
       if (abs(lonlr + 360._r8 - lonur) < abs(lonlr - lonur)) lonlr = lonlr + 360._r8
       if (abs(lonlr - 360._r8 - lonur) < abs(lonlr - lonur)) lonlr = lonlr - 360._r8

       nodeCoords(2*n-1) = 0.25_r8 * (lonur + lonul + lonll + lonlr)
       nodeCoords(2*n)   = 0.25_r8 * (latG(iur) + latG(iul) + latG(ill) + latG(ilr))
    enddo

    deallocate(lonG)
    deallocate(latG)

    ! Determine the pes that own each index of iurpts (nodeOwners)

    allocate(pes_local(nx_global*ny_global))
    allocate(pes_global(nx_global*ny_global))
    pes_local(:) = 0
    do n = 1,lsize
       pes_local(gindex(n)) = iam
    end do

    call ESMF_VMAllReduce(vm, sendData=pes_local, recvData=pes_global, count=nx_global*ny_global, &
         reduceflag=ESMF_REDUCE_SUM, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1,numNodes
       nodeOwners(n) = pes_global(iurpts(n))
    end do
    deallocate(pes_local)
    deallocate(pes_global)

    Emesh = ESMF_MeshCreate(parametricDim=2, &
         spatialDim=2, &
         coordSys=ESMF_COORDSYS_SPH_DEG, &
         nodeIds=nodeIds(1:numNodes), &
         nodeCoords=nodeCoords, &
         nodeOwners=nodeOwners, &
         elementIds=elemIds,&
         elementTypes=elemTypes, &
         elementConn=elemConn, &
         elementCoords=elemCoords, &
         rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    deallocate(iurpts)
    deallocate(nodeIds, nodeCoords, nodeOwners)
    deallocate(elemIds, elemTypes, elemConn, elemCoords)

  end subroutine dead_meshinit

end module dead_nuopc_mod
