
module dead_nuopc_mod

  use ESMF                  , only : ESMF_Gridcomp, ESMF_State, ESMF_StateGet
  use ESMF                  , only : ESMF_Clock, ESMF_Time, ESMF_TimeInterval, ESMF_Alarm
  use ESMF                  , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_ClockSet, ESMF_ClockAdvance, ESMF_AlarmSet
  use ESMF                  , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_METHOD_INITIALIZE
  use ESMF                  , only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF                  , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE
  use ESMF                  , only : operator(/=), operator(==), operator(+)
  use med_constants_mod, only : IN, R8, CS, CL
  use shr_file_mod          , only : shr_file_getunit, shr_file_freeunit
  use shr_sys_mod           , only : shr_sys_abort
  use shr_nuopc_utils_mod , only : shr_nuopc_utils_ChkErr

  implicit none
  private

!  public :: dead_setNewGrid
!  public :: dead_read_inparms
  public :: dead_init_nuopc
  public :: dead_run_nuopc
  public :: dead_final_nuopc
  public :: ModelInitPhase
  public :: ModelSetRunClock
  public :: fld_list_add
  public :: fld_list_realize
  public :: state_getimport
  public :: state_setexport
  public :: Print_FieldExchInfo

  private :: state_getfldptr

  ! !PUBLIC DATA MEMBERS:
  integer, public :: dead_grid_lat    = 1   ! lat from component
  integer, public :: dead_grid_lon    = 2   ! lon from component
  integer, public :: dead_grid_area   = 3   ! area from component
  integer, public :: dead_grid_mask   = 4   ! mask, 0 = inactive cell
  integer, public :: dead_grid_frac   = 5   ! fractional area coverage
  integer, public :: dead_grid_aream  = 6   ! area from mapping file
  integer, public :: dead_grid_index  = 7   ! global index
  integer, public :: dead_grid_pid    = 8   ! proc id number
  integer, public :: dead_grid_total  = 8

  type fld_list_type
    character(len=128) :: stdname
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
       nxg, nyg, decomp_type, nproc_x, seg_len, flood)
    use ESMF, only : ESMF_VMGetCurrent, ESMF_VM, ESMF_VMBroadcast, ESMF_VMGet
    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*)       , intent(in)    :: model
    character(len=*)      , intent(in)    :: inst_suffix ! char string associated with instance
    integer(IN)            , intent(in)    :: logunit     ! logging unit number
    integer(IN)		   , intent(out)   :: nproc_x
    integer(IN)            , intent(out)   :: seg_len
    integer(IN)            , intent(out)   :: nxg         ! global dim i-direction
    integer(IN)            , intent(out)   :: nyg         ! global dim j-direction
    integer(IN)            , intent(out)   :: decomp_type ! decomposition type
    logical                , intent(out)   :: flood       ! rof flood flag

    !--- local variables ---
    type(ESMF_VM) :: vm
    character(CL) :: fileName    ! generic file name
    integer(IN)   :: nunit       ! unit number
    integer(IN)   :: ierr        ! error code
    integer(IN)   :: unitn       ! Unit for namelist file
    integer(IN) :: tmp(6)     ! array for broadcast
    integer(IN) :: localPet   ! mpi id of current task in current context
    integer :: rc                  ! EMSF return code
    !--- formats ---
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
    flood = .false.

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return


    if (localPet==0) then
       unitn = shr_file_getUnit()
       open(unitn, file='x'//model//'_in'//trim(inst_suffix), status='old' )
       read(unitn,*) nxg
       read(unitn,*) nyg
       read(unitn,*) decomp_type
       read(unitn,*) nproc_x
       read(unitn,*) seg_len
       if (model.eq.'rof') then
          read(unitn,*) flood
       end if
       close (unitn)
       call shr_file_freeunit(unitn)
    endif

    tmp(1) = nxg
    tmp(2) = nyg
    tmp(3) = decomp_type
    tmp(4) = nproc_x
    tmp(5) = seg_len
    if (model.eq.'rof' .and. flood) then
       tmp(6) = 1
    else
       tmp(6) = 0
    endif
    call ESMF_VMBroadcast(vm, tmp, 6, 0, rc=rc)
    nxg = tmp(1)
    nyg = tmp(2)
    decomp_type = tmp(3)
    nproc_x = tmp(4)
    seg_len = tmp(5)
    if(tmp(6) == 1) then
       flood = .true.
    endif

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
       if (model.eq.'rof') then
          write(logunit,F01) ' Flood mode     :  ',flood
       endif
       write(logunit,F00) model
    end if
  end subroutine dead_read_inparms

  !===============================================================================
  subroutine dead_setNewGrid(decomp_type, nxg, nyg, logunit, lsize,  &
                             gbuf, seg_len, nproc_x)
    use ESMF, only : ESMF_VM, ESMF_VMGetCurrent, ESMF_VmGet
    use shr_const_mod, only : shr_const_pi, shr_const_rearth

    ! DESCRIPTION:
    ! This sets up some defaults.  The user may want to overwrite some
    ! of these fields in the main program after initialization in complete.

    ! input/output parameters:
    integer(IN) ,intent(in)          :: decomp_type !
    integer(IN) ,intent(in)          :: nxg,nyg     ! global grid sizes
    integer(IN) ,intent(in)          :: logunit     ! output logunit
    integer(IN) ,intent(out)         :: lsize       ! local grid sizes
    real(R8)    ,pointer             :: gbuf(:,:)   ! output data
    integer(IN) ,intent(in),optional :: seg_len     ! seg len decomp setting
    integer(IN) ,intent(in),optional :: nproc_x     ! 2d decomp setting

    !--- local ---
    type(ESMF_VM) :: vm
    integer(IN)   :: rc
    integer(IN) ::  mype
    integer(IN) :: totpe       ! total number of pes
    integer(IN)             :: ierr            ! error code
    logical                 :: found
    integer(IN)             :: i,j,ig,jg
    integer(IN)             :: n,ng,is,ie,js,je,nx,ny      ! indices
    integer(IN)             :: npesx,npesy,mypex,mypey,nxp,nyp
    real   (R8)             :: hscore,bscore
    real   (R8)             :: dx,dy,deg2rad,ys,yc,yn,area,re
    integer(IN),allocatable :: gindex(:)

    !--- formats ---
    character(*), parameter :: F00   = "('(dead_setNewGrid) ',8a)"
    character(*), parameter :: F01   = "('(dead_setNewGrid) ',a,4i8)"
    character(*), parameter :: subName = "(dead_setNewGrid) "
    !-------------------------------------------------------------------------------

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=mype, peCount=totpe, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (decomp_type == 1 .or. &
         decomp_type == 2 .or. &
         decomp_type == 3 .or. &
         decomp_type == 4 .or. &
         decomp_type == 11) then
       ! valid
    else
       !-------------------------------------------------------------------------
       ! invalid decomposition type
       !-------------------------------------------------------------------------
       if (mype == 0) then
          write(logunit,F01) 'ERROR: invalid decomp_type = ',decomp_type
       end if
       call shr_sys_abort(subName//'invalid decomp_type')
    endif

    if (nxg*nyg == 0) then
       lsize = 0
       allocate(gbuf(lsize,dead_grid_total))
       !      gbuf = -888.0_R8
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
       gbuf(n,dead_grid_mask ) = 1
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
    integer                  :: ierr          ! error code
    integer                  :: local_comm    ! local communicator
    integer                  :: mype          ! pe info
    integer                  :: totpe         ! total number of pes
    integer                  :: nproc_x
    integer                  :: seg_len
    integer                  :: decomp_type
    logical                  :: flood=.false. ! rof flood flag
    character(*), parameter  :: subName = "(dead_init_nuopc) "
    !-------------------------------------------------------------------------------

    ! Read input parms

    call dead_read_inparms(model, inst_suffix, logunit, &
         nxg, nyg, decomp_type, nproc_x, seg_len, flood)

    ! Initialize grid

    call dead_setNewGrid(decomp_type, nxg, nyg, logunit, &
         lsize, gbuf, seg_len, nproc_x)

  end subroutine dead_init_nuopc

  !===============================================================================
  subroutine dead_run_nuopc(model, d2x, gbuf, flds_d2x)

    use shr_const_mod  , only : shr_const_pi
    use shr_string_mod , only : shr_string_listGetIndexF

    ! DESCRIPTION: run method for dead model

    ! input/output parameters:
    character(len=*) , intent(in)    :: model
    real(r8)         , intent(inout) :: d2x(:,:)    ! dead   -> driver
    real(r8)         , pointer       :: gbuf(:,:)   ! model grid
    character(len=*) , intent(in)    :: flds_d2x    ! list of fields to dead -> driver

    !--- local ---
    integer                 :: n                 ! index
    integer                 :: nf                ! fields loop index
    integer                 :: ki                ! index
    integer                 :: lsize             ! size of AttrVect
    real(R8)                :: lat               ! latitude
    real(R8)                :: lon               ! longitude
    integer                 :: nflds_d2x
    integer                 :: ncomp
    character(*), parameter :: F04   = "('(',a,'_run_nuopc) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(dead_run_nuopc) "
    !-------------------------------------------------------------------------------

    nflds_d2x = size(d2x, dim=1)

    ! PACK (currently no unpacking)

    selectcase(model)
    case('atm')
       ncomp = 1
    case('lnd')
       ncomp = 2
    case('ice')
       ncomp = 3
    case('ocn')
       ncomp = 4
    case('glc')
       ncomp = 5
    case('rof')
       ncomp = 6
    case('wav')
       ncomp = 7
    end select

    nflds_d2x = size(d2x, dim=1)
    lsize = size(d2x, dim=2)

    if (model.eq.'rof') then

       do nf=1,nflds_d2x
          do n=1,lsize
             d2x(nf,n) = (nf+1) * 1.0_r8
          enddo
       enddo

    else if (model.eq.'glc') then

       do nf=1,nflds_d2x
          do n=1,lsize
             lon = gbuf(n,dead_grid_lon)
             lat = gbuf(n,dead_grid_lat)
             d2x(nf,n) = (nf*100)                          &
                  *  cos (SHR_CONST_PI*lat/180.0_R8)       &
                  *  cos (SHR_CONST_PI*lat/180.0_R8)       &
                  *  sin (SHR_CONST_PI*lon/180.0_R8)       &
                  *  sin (SHR_CONST_PI*lon/180.0_R8)       &
                  + (ncomp*10.0_R8)
          enddo
       enddo

    else

       do nf=1,nflds_d2x
          do n=1,lsize
             lon = gbuf(n,dead_grid_lon)
             lat = gbuf(n,dead_grid_lat)
             d2x(nf,n) = (nf*100)                          &
                  *  cos (SHR_CONST_PI*lat/180.0_R8)       &
                  *  sin((SHR_CONST_PI*lon/180.0_R8)       &
                  -      (ncomp-1)*(SHR_CONST_PI/3.0_R8) ) &
                  + (ncomp*10.0_R8)
          enddo
       enddo

    endif

    selectcase(model)
    case('ice')

       ki = shr_string_listGetIndexF(flds_d2x, "Si_ifrac")
       d2x(ki,:) = min(1.0_R8,max(0.0_R8,d2x(ki,:)))

       ki = shr_string_listGetIndexF(flds_d2x, "Si_imask")
       d2x(ki,:) = float(nint(min(1.0_R8,max(0.0_R8,d2x(ki,:)))))

    case('ocn')

       ki = shr_string_listGetIndexF(flds_d2x, "So_omask")
       d2x(ki,:) = float(nint(min(1.0_R8,max(0.0_R8,d2x(ki,:)))))

    case('lnd')

       ki = shr_string_listGetIndexF(flds_d2x, "Sl_lfrin")
       d2x(ki,:) = 1.0_R8

    case('glc')

       ki = shr_string_listGetIndexF(flds_d2x, "Sg_icemask")
       d2x(ki,:) = 1.0_R8

       ki = shr_string_listGetIndexF(flds_d2x, "Sg_icemask_coupled_fluxes")
       d2x(ki,:) = 1.0_R8

       ki = shr_string_listGetIndexF(flds_d2x, "Sg_ice_covered")
       d2x(ki,:) = 1.0_R8

    end select

  end subroutine dead_run_nuopc

  !===============================================================================
  subroutine dead_final_nuopc(model, logunit)
    use ESMF, only : ESMF_VM, ESMF_VMGetCurrent, ESMF_VMGet
    ! DESCRIPTION: finalize method for datm model

    ! input/output parameters:
    character(len=*) , intent(in) :: model
    integer          , intent(in) :: logunit     ! logging unit number

    !-- local --
    type(ESMF_VM) :: vm
    integer :: rc
    integer :: localPet

    !--- formats ---
    character(*), parameter :: F00   = "('(dead_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(dead_comp_final) ',73('-'))"
    character(*), parameter :: subName = "(dead_comp_final) "
    !-------------------------------------------------------------------------------

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (localPet==0) then
       write(logunit,F91)
       write(logunit,F00) trim(model),': end of main integration loop'
       write(logunit,F91)
    end if

  end subroutine dead_final_nuopc

  !===============================================================================
  subroutine fld_list_add(num, fldlist, stdname, flds_concat)

    use ESMF, only : ESMF_LogWrite, ESMF_LOGMSG_ERROR

    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname
    character(len=*), optional, intent(inout) :: flds_concat

    ! local variables
    integer :: rc
    integer :: dbrc
    character(len=*), parameter :: subname='(dshr_nuopc_mod:fld_list_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
       return
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(flds_concat)) then
       if (len_trim(flds_concat) + len_trim(stdname) + 1 >= len(flds_concat)) then
          call ESMF_LogWrite(subname//': ERROR: max len of flds_concat has been exceeded', &
               ESMF_LOGMSG_ERROR, line=__LINE__, file= u_FILE_u, rc=dbrc)
       end if
       if (trim(flds_concat) == '') then
          flds_concat = trim(stdname)
       else
          flds_concat = trim(flds_concat)//':'//trim(stdname)
       end if
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
    integer                :: dbrc
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(dshr_nuopc_mod:fld_list_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the field
             field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
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
      character(len=*), parameter :: subname='(dshr_nuopc_mod:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), rc=rc)
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
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine ModelInitPhase

  !===============================================================================
  subroutine ModelSetRunClock(gcomp, rc)

    use shr_nuopc_time_mod , only : shr_nuopc_time_alarmInit
    use ESMF               , only : ESMF_ClockGetAlarmList, ESMF_ALARMLIST_ALL
    use NUOPC_Model        , only : NUOPC_ModelGet
    use NUOPC              , only : NUOPC_CompAttributeGet

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
    integer                  :: dbrc
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname='dshr_nuopc_mod:(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart alarm
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO, rc=dbrc)

       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call shr_nuopc_time_alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelSetRunClock

  !===============================================================================
  subroutine state_getimport(state, fldname, output, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    character(len=*)    , intent(in)    :: fldname
    real(r8)            , intent(out)   :: output(:)
    integer             , intent(out)   :: rc

    ! local variables
    integer                     :: g, i
    real(R8), pointer           :: fldptr(:)
    type(ESMF_StateItem_Flag)   :: itemFlag
    integer                     :: dbrc
    character(len=*), parameter :: subname='(lnd_import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine if field with name fldname exists in state
    call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! if field exists then create output array - else do nothing
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then

       ! get field pointer
       call state_getfldptr(state, trim(fldname), fldptr, rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! determine output array
       do g = 1,size(fldptr)
          output(g) = fldptr(g)
       end do
    end if

  end subroutine state_getimport

  !===============================================================================

  subroutine state_setexport(state, fldname, input,  rc)
    ! ----------------------------------------------
    ! Map input array to export state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    character(len=*)    , intent(in)    :: fldname
    real(r8)            , intent(in)    :: input(:)
    integer             , intent(out)   :: rc

    ! local variables
    integer                     :: g, i
    real(R8), pointer           :: fldptr(:)
    type(ESMF_StateItem_Flag)   :: itemFlag
    integer                     :: dbrc
    character(len=*), parameter :: subname='(lnd_import_export:state_setexport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine if field with name fldname exists in state
    call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! if field exists then create output array - else do nothing
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then

       ! get field pointer
       call state_getfldptr(state, trim(fldname), fldptr,  rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! set fldptr values to input array
       do g = 1,size(fldptr)
          fldptr(g) = input(g)
       end do
    end if

  end subroutine state_setexport

  !===============================================================================

  subroutine state_getfldptr(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------
    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    type(ESMF_State),  intent(in)    :: State
    character(len=*),  intent(in)    :: fldname
    real(R8), pointer, intent(out)   :: fldptr(:)
    integer,           intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    type(ESMF_Mesh)             :: lmesh
    integer                     :: dbrc
    integer                     :: nnodes, nelements
    character(len=*), parameter :: subname='(lnd_import_export:state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, status=status, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    else
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (nnodes == 0 .and. nelements == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", ESMF_LOGMSG_INFO, rc=dbrc)
          rc = ESMF_FAILURE
          return
       end if

       call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif  ! status

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine state_getfldptr

  !===============================================================================

  subroutine Print_FieldExchInfo(values, logunit, fldlist, nflds, istr)

    use med_constants_mod   , only : R8
    use ESMF                , only : ESMF_MAXSTR

    ! !DESCRIPTION:
    ! Print out information about values to stdount
    ! - flag sets the level of information:
    ! - print out names of fields in values 2d array
    ! - also print out local max and min of data in values 2d array
    ! If optional argument istr is present, it will be output before any of the information.


    ! input/output parameters:
    real(R8)            , intent(in)          :: values(:,:) ! arrays sent to/recieved from mediator
    integer             , intent(in)          :: logunit
    type(fld_list_type) , intent(in)          :: fldlist(:)
    integer             , intent(in)          :: nflds
    character(*)        , intent(in),optional :: istr  ! string for print

    !--- local ---
    integer                    :: n           ! generic indicies
    integer                    :: nsize       ! grid point in values array
    real(R8)                   :: minl(nflds) ! local min
    real(R8)                   :: maxl(nflds) ! local max
    character(len=ESMF_MAXSTR) :: name

    !--- formats ---
    character(*),parameter :: subName = '(print_FieldExchInfo) '
    character(*),parameter :: F00 = "('(print_FieldExchInfo) ',8a)"
    character(*),parameter :: F01 = "('(print_FieldExchInfo) ',a,i9)"
    character(*),parameter :: F02 = "('(print_FieldExchInfo) ',a,2es11.3,i4,2x,a)"
    !-------------------------------------------------------------------------------

    if (present(istr)) write(logunit,*) trim(istr)
    nsize = size(values, dim=2)
    write(logunit,F01) "local size =",nsize
    do n = 1, nflds
       minl(n) = minval(values(n,:))
       maxl(n) = maxval(values(n,:))
       write(logunit,F02) 'l min/max ',minl(n),maxl(n),n,fldlist(n)%stdname
    enddo

  end subroutine Print_FieldExchInfo

end module dead_nuopc_mod
