module dshr_nuopc_mod

  use NUOPC
  use NUOPC_Model
  use ESMF
  use mct_mod          , only : mct_gsmap_init, mct_avect_lsize, mct_avect_indexra
  use dshr_methods_mod , only : alarmInit, chkerr, state_getfldptr
  use dshr_methods_mod , only : get_component_instance
  use shr_strdata_mod  , only : shr_strdata_type
  use shr_kind_mod     , only : r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx
  use shr_string_mod   , only : shr_string_listGetIndex
  use shr_sys_mod      , only : shr_sys_abort

  implicit none
  public

  ! remove the following when no longer eaded
  public :: ModelInitPhase   
  public :: ModelSetRunClock 
  public :: dshr_avect_add
  public :: dshr_translate_add
  public :: dshr_import
  public :: dshr_export
  ! end remove

  public :: dshr_set_runclock
  public :: dshr_model_initphase
  public :: dshr_advertise
  public :: dshr_fld_add
  public :: dshr_dfield_add
  public :: dshr_realize
  public :: dshr_sdat_init
  public :: dshr_streams_copy
  public :: dshr_restart_read
  public :: dshr_restart_write
  public :: dshr_check_mesh

  interface dshr_dfield_add
     module procedure dshr_dfield_add_1d
     module procedure dshr_dfield_add_2d
  end interface dshr_dfield_add

  type fld_list_type
    character(len=128) :: stdname
    integer :: ungridded_lbound = 0
    integer :: ungridded_ubound = 0
  end type fld_list_type

  type dfield_type
     character(CS)              :: state_fldname
     ! state data
     real(r8), pointer          :: state_data1d(:) => null()
     real(r8), pointer          :: state_data2d(:,:) => null()
     ! stream data input (always assumed to be 1d for now)
     real(r8), pointer          :: stream_data1d(:) => null()
     ! stream data pointers for 1d state data
     integer                    :: sdat_stream_index = 0
     integer                    :: sdat_avect_index = 0
     character(CS)              :: stream_fldname = 'unset'
     ! stream data pointers for 2d state data
     integer, pointer           :: sdat_stream_indices(:) => null()
     integer, pointer           :: sdat_avect_indices(:) => null()
     character(CS), pointer     :: stream_fldnames(:)  => null()
  end type dfield_type

  ! Note that gridTofieldMap = 2, therefore the ungridded dimension is innermost

  integer                 :: iunset = -999
  integer     , parameter :: fldsMax = 100
  integer     , parameter :: dbug = 10
  character(*), parameter :: modName =  "(dhsr_nuopc_mod)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dshr_advertise(gcomp, mpicom, my_task,  inst_index, inst_suffix, &
       flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    integer          , intent(inout) :: mpicom
    integer          , intent(out) :: my_task
    integer          , intent(out) :: inst_index
    character(len=*) , intent(out) :: inst_suffix
    character(len=*) , intent(out) :: flds_scalar_name
    integer          , intent(out) :: flds_scalar_num
    integer          , intent(out) :: flds_scalar_index_nx
    integer          , intent(out) :: flds_scalar_index_ny
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_VM)      :: vm
    character(len=CL)  :: cvalue
    character(len=CL)  :: logmsg
    logical            :: isPresent, isSet
    character(len=*),parameter  :: subname='(dshr_advertise)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! generate local mpi comm
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, mpiCommunicator=mpicom, localPet=my_task, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine instance information
    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get scalar attributes
    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

  end subroutine dshr_advertise

!===============================================================================

  subroutine dshr_sdat_init(mpicom, compid, my_task, master_task, logunit, &
       scmmode, scmlon, scmlat, clock, mesh,  model_name, sdat, &
       dmodel_domain_fracname_from_stream, reset_domain_mask, rc)

    ! ----------------------------------------------
    ! Initialize sdat
    ! ----------------------------------------------

    use shr_strdata_mod , only : shr_strdata_pioinit
    use shr_strdata_mod , only : shr_strdata_init_model_domain
    use shr_strdata_mod , only : shr_strdata_init_streams
    use shr_strdata_mod , only : shr_strdata_init_mapping
    use shr_strdata_mod , only : shr_strdata_print
    use shr_cal_mod     , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_calendarname

    ! input/output variables
    integer                    , intent(in)    :: mpicom   ! mpi communicator
    integer                    , intent(in)    :: compid
    integer                    , intent(in)    :: my_task
    integer                    , intent(in)    :: master_task
    integer                    , intent(in)    :: logunit
    logical                    , intent(in)    :: scmMode
    real(r8)                   , intent(in)    :: scmlat
    real(r8)                   , intent(in)    :: scmlon 
    type(ESMF_Clock)           , intent(in)    :: clock
    type(ESMF_Mesh)            , intent(in)    :: mesh
    character(len=*)           , intent(in)    :: model_name
    type(shr_strdata_type)     , intent(inout) :: sdat
    character(len=*), optional , intent(in)    :: dmodel_domain_fracname_from_stream
    logical         , optional , intent(in)    :: reset_domain_mask
    integer                    , intent(out)   :: rc
    
    ! local varaibles
    integer                      :: n,k          ! generic counters
    integer                      :: lsize        ! local size
    integer                      :: gsize
    type(ESMF_DistGrid)          :: distGrid
    integer, allocatable, target :: gindex(:)
    integer                      :: dimCount
    integer                      :: tileCount
    integer                      :: deCount
    integer, allocatable         :: elementCountPTile(:)
    integer, allocatable         :: indexCountPDE(:,:)
    type(ESMF_CalKind_Flag)      :: esmf_caltype ! esmf calendar type
    character(len=CS)            :: calendar     ! calendar name
    character(len=*), parameter  :: subname='(dshr_nuopc_mod:dshr_sdat_init)'
    character(*)    , parameter  :: F01="('(dshr_init_strdata) ',a,2f10.4)"
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! obtain the distgrid from the mesh that was read in
    call ESMF_MeshGet(mesh, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine local size on my processor
    call ESMF_DistGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine global index space for my processor
    allocate(gindex(lsize))
    call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine global size of distgrid
    call ESMF_DistGridGet(distGrid, dimCount=dimCount, deCount=deCount, tileCount=tileCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(elementCountPTile(tileCount))
    call ESMF_distGridGet(distGrid, elementCountPTile=elementCountPTile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    gsize = 0
    do n = 1,size(elementCountPTile)
       gsize = gsize + elementCountPTile(n)
    end do
    deallocate(elementCountPTile)

    ! initialize sdat%gsmap (the data mdel gsmap)
    call mct_gsMap_init(sdat%gsmap, gindex, mpicom, compid, lsize, gsize)
    deallocate(gindex)

    ! initialize shr_strdata pio
    call shr_strdata_pioinit(sdat, compid)

    ! initialize sdat calendar
    call ESMF_ClockGet(clock, calkindflag=esmf_caltype, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call shr_sys_abort(subname//" ERROR bad ESMF calendar name "//trim(calendar))
    end if
    sdat%calendar = trim(shr_cal_calendarName(trim(calendar)))

    ! initialize sdat domain (sdat%grid)
    if (my_task == master_task) then
       write(logunit,*) ' scm mode, lon lat = ',scmmode, scmlon,scmlat
       if (present(reset_domain_mask)) then
          write(logunit,*) ' resetting domain mask'
       end if
       if (present(dmodel_domain_fracname_from_stream)) then
          write(logunit,*)' reading fracname ',trim(dmodel_domain_fracname_from_stream),&
               ' from the domain of the first stream'
       end if
    end if
    call shr_strdata_init_model_domain(sdat, mpicom, compid, my_task, &
         scmmode=scmmode, scmlon=scmlon, scmlat=scmlat, gsmap=sdat%gsmap, &
         dmodel_domain_fracname_from_stream=dmodel_domain_fracname_from_stream, &
         reset_domain_mask=reset_domain_mask)

    ! initialize sdat attributes for streams and mapping of streams to model domain
    call shr_strdata_init_streams(sdat, compid, mpicom, my_task)
    call shr_strdata_init_mapping(sdat, compid, mpicom, my_task)

    if (my_task == master_task) then
       call shr_strdata_print(sdat,'SDAT data from '//trim(model_name))
    endif

  end subroutine dshr_sdat_init

!===============================================================================

  subroutine dshr_check_mesh (mesh, sdat, model_name, rc)

    type(ESMF_MESH)        , intent(in)  :: mesh
    type(shr_strdata_type) , intent(in)  :: sdat
    character(len=*)       , intent(in)  :: model_name
    integer                , intent(out) :: rc

    ! local variables
    integer           :: n,lsize
    integer           :: klat, klon, kfrac  ! AV indices
    real(R8)          :: domlon,domlat      ! domain lats and lots
    integer           :: spatialDim         ! number of dimension in mesh
    integer           :: numOwnedElements   ! size of mesh
    real(R8), pointer :: ownedElemCoords(:) ! mesh lat and lons
    real(r8), pointer :: xc(:), yc(:)       ! mesh lats and lons
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! get local size from first stream attribute vector of sdat
    lsize = mct_avect_lsize(sdat%avs(1))

    ! obtain mesh lats and lons
    call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    allocate(xc(numOwnedElements), yc(numOwnedElements))
    call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (numOwnedElements /= lsize) then
       call shr_sys_abort('ERROR: numOwnedElements is not equal to lsize')
    end if
    do n = 1,lsize
       xc(n) = ownedElemCoords(2*n-1)
       yc(n) = ownedElemCoords(2*n)
    end do

    ! obtain sdat lat and lons
    klon = mct_aVect_indexRA(sdat%grid%data,'lon')
    klat = mct_aVect_indexRA(sdat%grid%data,'lat')
    do n = 1, lsize
       domlon = sdat%grid%data%rattr(klon,n)
       domlat = sdat%grid%data%rattr(klat,n)
       if (abs( domlon - xc(n)) > 1.e-10 .and. domlon /= 0.0_r8) then
          write(6,100) 'ERROR: '//trim(model_name)//' n, dom_lon, mesh_lon, diff_lon = ',n, domlon, xc(n), abs(xc(n)-domlon)
          call shr_sys_abort()
       end if
       if (abs( domlat - yc(n)) > 1.e-10 .and. domlat /= 0.0_r8) then
          write(6,100) 'ERROR: '//trim(model_name)//' n, dom_lat, mesh_lat, diff_lat = ',n, domlat, yc(n), abs(yc(n)-domlat)
          call shr_sys_abort()
       end if
100    format(a,i6,2(f21.13,3x),d21.5)
       !SDAT%grid%data%rattr(klon,n) = xc(n)
       !SDAT%grid%data%rattr(klat,n) = yc(n)
    end do
    deallocate(xc, yc)

  end subroutine dshr_check_mesh

!===============================================================================

  subroutine dshr_fld_add(fldname, fldlist_num, fldlist, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    character(len=*)               , intent(in)    :: fldname
    integer                        , intent(inout) :: fldlist_num
    type(fld_list_type)            , intent(inout) :: fldlist(:)
    integer             , optional , intent(in)    :: ungridded_lbound
    integer             , optional , intent(in)    :: ungridded_ubound

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(dshr_nuopc_mod:dshr_fld_add)'
    ! ----------------------------------------------

    ! Set up a list of field information

    fldlist_num = fldlist_num + 1
    if (fldlist_num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(fldname), ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE
      return
    endif
    fldlist(fldlist_num)%stdname = trim(fldname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(fldlist_num)%ungridded_lbound = ungridded_lbound
       fldlist(fldlist_num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine dshr_fld_add

!===============================================================================

  subroutine dshr_dfield_add_1d(sdat, state, dfields, dfields_num, state_fld, strm_fld, strm_ptr, state_ptr, rc)

    ! Set 1d dfield values

    type(shr_strdata_type) , intent(in)    :: sdat
    type(ESMF_State)       , intent(inout) :: state
    type(dfield_type)      , intent(inout) :: dfields(:)
    integer                , intent(inout) :: dfields_num
    character(len=*)       , intent(in)    :: state_fld
    character(len=*)       , intent(in)    :: strm_fld
    real(r8), optional     , pointer       :: strm_ptr(:)
    real(r8), optional     , pointer       :: state_ptr(:)
    integer                , intent(out)   :: rc

    ! local variables
    integer :: ns, kf
    integer :: lsize, num
    logical :: found
    character(len=*), parameter :: subname='(dfield_add_1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! error checks
    dfields_num = dfields_num + 1
    if (dfields_num > fldsMax) then
       call shr_sys_abort(trim(subname)//": ERROR num > fldsMax for "//trim(state_fld))
    endif
    num = dfields_num

    ! set the import or export state field name
    dfields(num)%state_fldname = trim(state_fld)


    ! Initialize strm_ptr and state_ptr if it is present
    ! These will be set to valid values if the relevant fields are found
    if (present(strm_ptr)) strm_ptr => null()
    if (present(state_ptr)) state_ptr => null()

    ! determine local size
    lsize = mct_avect_lsize(sdat%avs(1))

    ! always allocate memory for stream_data and initialize it to 0
    ! note that if the attribute vector is not found in the streams
    allocate(dfields(num)%stream_data1d(lsize))
    dfields(num)%stream_data1d(:) = 0._r8

    ! loop over all input streams
    found = .false.
    do ns = 1, sdat%nstreams
       ! determine if the strm_fld is in the attribute vector of stream ns
       kf  = mct_aVect_indexRA(sdat%avs(ns), trim(strm_fld), perrWith='quiet')
       if (kf > 0) then
          ! the field is in the attribute vector - so set the following values
          dfields(num)%sdat_stream_index = ns
          dfields(num)%sdat_avect_index = kf
          ! set strm_ptr if argument is present
          if (present(strm_ptr)) strm_ptr => dfields(num)%stream_data1d
          found = .true.
          exit
       end if
    end do
    if (.not. found) then
       dfields(num)%sdat_stream_index = iunset
       dfields(num)%sdat_avect_index = iunset
    end if

    ! set export state pointer array values
    if (trim(state_fld) /= 'unset') then

       ! Set export state array pointer
       call state_getfldptr(State, fldname=trim(state_fld), fldptr1=dfields(num)%state_data1d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       dfields(num)%state_data1d = 0.0_r8

       ! Return array pointer if argument is present
       if (present(state_ptr)) then
          state_ptr => dfields(num)%state_data1d
       end if
    end if

  end subroutine dshr_dfield_add_1d

  !===============================================================================

  subroutine dshr_dfield_add_2d(sdat, state, dfields, dfields_num, state_fld, strm_flds, state_ptr, rc)

    ! input/output variables
    type(shr_strdata_type) , intent(in)    :: sdat
    type(ESMF_State)       , intent(inout) :: state
    type(dfield_type)      , intent(inout) :: dfields(:)
    integer                , intent(inout) :: dfields_num
    character(len=*)       , intent(in)    :: state_fld
    character(len=*)       , intent(in)    :: strm_flds(:)
    real(r8), optional     , pointer       :: state_ptr(:,:)
    integer                , intent(out)   :: rc

    ! local variables
    integer :: n, i, kf, ns, nf
    integer :: nflds, lsize, num
    logical :: found
    character(len=*), parameter :: subname='(dfield_add_2d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! error checks
    dfields_num = dfields_num + 1
    if (dfields_num > fldsMax) then
       call shr_sys_abort(trim(subname)//": ERROR num > fldsMax for "//trim(state_fld))
    endif
    num = dfields_num

    ! assume only one only two-dimensional state field
    dfields(num)%state_fldname = state_fld

    ! determine stream fldnames array
    nflds = size(strm_flds)
    allocate(dfields(num)%sdat_stream_indices(nflds))
    allocate(dfields(num)%sdat_avect_indices(nflds))
    allocate(dfields(num)%stream_fldnames(nflds))
    do n = 1, nflds
       dfields(num)%stream_fldnames(n) = trim(strm_flds(n))
    end do

    ! determine local size
    lsize = mct_avect_lsize(sdat%avs(1))

    ! loop through input array of stream field names
    do nf = 1, nflds
       ! loop through input streams
       do ns = 1, sdat%nstreams
          ! determine if the strm_flds(nf) is in the attribute vector of stream ns
          kf = mct_aVect_indexRA(sdat%avs(ns), trim(strm_flds(nf)), perrWith='quiet')
          if (kf > 0) then
             dfields(num)%sdat_stream_indices(nf) = ns
             dfields(num)%sdat_avect_indices(nf) = kf
          end if
       end do
    end do

    ! set export state pointer array values
    if (trim(state_fld) /= 'unset') then

       ! Set export state array pointer
       !allocate(dfields(num)%state_data2d(nflds,lsize))
       dfields(num)%state_data2d(nflds,lsize) = 0._r8       
       call state_getfldptr(State, fldname=trim(state_fld), fldptr2=dfields(num)%state_data2d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Return array pointer if argument is present
       if (present(state_ptr)) then
          state_ptr => dfields(num)%state_data2d
       end if
    end if

  end subroutine dshr_dfield_add_2d

  !===============================================================================

  subroutine dshr_streams_copy(dfields, dfields_num, sdat, rc)

    ! Copy stream data into dfield data type for each element of dfields

    ! input/output variables
    type(dfield_type)      , intent(inout) :: dfields(:)
    integer                , intent(in)    :: dfields_num
    type(shr_strdata_type) , intent(in)    :: sdat
    integer                , intent(out)   :: rc

    ! local variables
    integer :: n, i, k
    integer :: avect_index
    integer :: stream_index
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Loop over all dfield entries and fill in stream_data and state_data1d or state_data2d arrays
    do n = 1,dfields_num

       ! Map the stream data to the state data
       if (associated(dfields(n)%state_data1d)) then

          stream_index = dfields(n)%sdat_stream_index
          avect_index =  dfields(n)%sdat_avect_index
          if (stream_index /= iunset .and. avect_index /= iunset) then
             dfields(n)%stream_data1d(:) = sdat%avs(stream_index)%rattr(avect_index,:)
             dfields(n)%state_data1d(:) = dfields(n)%stream_data1d(:)
          end if

       else if (associated(dfields(n)%state_data2d)) then

          do i = 1,size(dfields(n)%stream_fldnames)
             stream_index = dfields(n)%sdat_stream_indices(i)
             avect_index = dfields(n)%sdat_avect_indices(i)
             dfields(n)%state_data2d(i,:) = sdat%avs(stream_index)%rattr(avect_index,:)
          end do

       end if
    end do

  end subroutine dshr_streams_copy

!===============================================================================

  subroutine dshr_avect_add(model_fld, model_fld_concat, av_index)

    ! input/output variables
    character(len=*)               , intent(in)    :: model_fld
    character(len=*)               , intent(inout) :: model_fld_concat
    integer             , optional , intent(out)   :: av_index

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(dshr_nuopc_mod:dshr_dmodel_add1)'
    ! ----------------------------------------------

    if (len_trim(model_fld_concat) + len_trim(model_fld) + 1 >= len(model_fld_concat)) then
       call ESMF_LogWrite(subname//': ERROR: max len of model_fld_concat has been exceeded', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    if (trim(model_fld_concat) == '') then
       model_fld_concat = trim(model_fld)
    else
       model_fld_concat = trim(model_fld_concat)//':'//trim(model_fld)
    end if

    if (present(av_index)) then
       call shr_string_listGetIndex(trim(model_fld_concat), trim(model_fld),  av_index)
    end if

  end subroutine dshr_avect_add

  !===============================================================================

  subroutine dshr_translate_add( fldname1, fldname2, fldnames1, fldnames2)

    !----------------------------------
    ! Create new character arrays fldnames1 and fldnames2
    !----------------------------------

    ! input/output variables
    character(len=*) , intent(in)    :: fldname1
    character(len=*) , intent(in)    :: fldname2
    character(len=*) , pointer       :: fldnames1(:)
    character(len=*) , pointer       :: fldnames2(:)

    ! local variables
    integer                     :: rc
    integer                     :: n, oldsize, id
    character(len=CS), pointer  :: new_fldnames1(:)
    character(len=CS), pointer  :: new_fldnames2(:)
    character(len=*), parameter :: subname='(dshr_nuopc_mod:dshr_translate_add) '
    ! ----------------------------------------------

    ! 1) determine new index
    if (associated(fldnames1)) then
       oldsize = size(fldnames1)
    else
       oldsize = 0
    end if
    id = oldsize + 1

    ! 2) allocate new_fldnames1 and new_fldnames2 to one element larger than input
    allocate(new_fldnames1(id))
    allocate(new_fldnames2(id))

    ! 3) copy fldnames1 and fldnames2 into first N-1 elements offldnames1 and fldnames2
    do n = 1,oldsize
       new_fldnames1(n) = fldnames1(n)
       new_fldnames2(n) = fldnames2(n)
    end do

    ! 4) deallocate / nullify stream_fldnames and model_fldnames
    if (oldsize >  0) then
       deallocate(fldnames1)
       deallocate(fldnames2)
       nullify(fldnames1)
       nullify(fldnames2)
    end if

    ! 5) point fldnames1 => new_fldnames1 and fldnames2 => new_fldnames2 and update info for new entry
    fldnames1  => new_fldnames1
    fldnames2  => new_fldnames2
    fldnames1(id) = trim(fldname1)
    fldnames2(id) = trim(fldname2)

  end subroutine dshr_translate_add

  !===============================================================================

  subroutine dshr_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
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
                  ESMF_LOGMSG_INFO)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          else
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                     ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                     ungriddedUbound=(/fldlist(n)%ungridded_ubound/), gridToFieldMap=(/2/), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
             end if
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO)
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
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine dshr_realize

  !===============================================================================

  subroutine dshr_model_initphase(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine dshr_model_initphase

  !===============================================================================

  subroutine ModelInitPhase(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine ModelInitPhase

  !===============================================================================

  subroutine dshr_set_runclock(gcomp, rc)
  
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
    character(len=*),parameter :: subname='dshr_nuopc_mod:(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart alarm
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine dshr_set_runclock

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

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
    character(len=*),parameter :: subname='dshr_nuopc_mod:(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart alarm
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine dshr_export(array, state, fldname, ungridded_index, rc)

    ! ----------------------------------
    ! copy array data to state fields
    ! ----------------------------------

    ! input/otuput variables
    real(r8)         , intent(inout) :: array(:)
    type(ESMF_State) , intent(inout) :: state
    character(len=*) , intent(in)    :: fldname
    integer, optional, intent(in)    :: ungridded_index
    integer          , intent(out)   :: rc

    ! local variables
    integer           :: lsize, n
    type(ESMF_Field)  :: lfield
    real(R8), pointer :: farray1d(:)
    real(R8), pointer :: farray2d(:,:)
    character(*),parameter :: subName = "(dshr_nuopc_mod: dshr_export)"
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemName=trim(fldname), field=lfield, rc=rc)
    if (.not. ChkErr(rc,__LINE__,u_FILE_u)) then
       call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" copy", ESMF_LOGMSG_INFO)

       lsize = size(array)
       if (present(ungridded_index)) then
          call ESMF_FieldGet(lfield, farrayPtr=farray2d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n = 1,lsize
             farray2d(ungridded_index,n) = array(n)
          end do
       else
          call ESMF_FieldGet(lfield, farrayPtr=farray1d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n = 1,lsize
             farray1d(n) = array(n)
          enddo
       end if
    end if

  end subroutine dshr_export

  !===============================================================================

  subroutine dshr_import(state, fldname, array, ungridded_index, rc)

    ! ----------------------------------
    ! copy state field to array data
    ! ----------------------------------

    ! input/output variables
    type(ESMF_State)  , intent(in)    :: state
    character(len=*)  , intent(in)    :: fldname
    real(r8)          , intent(inout) :: array(:)
    integer, optional , intent(in)    :: ungridded_index
    integer           , intent(out)   :: rc

    ! local variables
    integer           :: lsize, n
    type(ESMF_Field)  :: lfield
    real(R8), pointer :: farray1d(:)
    real(R8), pointer :: farray2d(:,:)
    character(*),parameter :: subName = "(dshr_nuopc_mod: dshr_import)"
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemName=trim(fldname), field=lfield, rc=rc)
    if (.not. ChkErr(rc,__LINE__,u_FILE_u)) then
       call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" copy", ESMF_LOGMSG_INFO)

       lsize = size(array)
       if (present(ungridded_index)) then
          call ESMF_FieldGet(lfield, farrayPtr=farray2d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n = 1,lsize
             array(n) = farray2d(ungridded_index,n)
          enddo
       else
          call ESMF_FieldGet(lfield, farrayPtr=farray1d,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n = 1,lsize
             array(n) = farray1d(n)
          enddo
       end if
    end if

   end subroutine dshr_import

   !===============================================================================

   subroutine dshr_restart_read(rest_file, rest_file_strm, rpfile, inst_suffix, nullstr, &
        logunit, my_task, master_task, mpicom, sdat, fld, fldname)

     use shr_pcdf_mod    , only : shr_pcdf_readwrite
     use shr_strdata_mod , only : shr_strdata_restRead, shr_strdata_type
     use shr_mpi_mod     , only : shr_mpi_bcast

     ! input/output arguments
     character(len=*)            , intent(inout) :: rest_file
     character(len=*)            , intent(inout) :: rest_file_strm
     character(len=*)            , intent(in)    :: rpfile
     character(len=*)            , intent(in)    :: inst_suffix
     character(len=*)            , intent(in)    :: nullstr
     integer                     , intent(in)    :: logunit
     integer                     , intent(in)    :: my_task
     integer                     , intent(in)    :: master_task
     integer                     , intent(in)    :: mpicom
     type(shr_strdata_type)      , intent(inout) :: sdat
     real(r8)         , optional , pointer       :: fld(:)
     character(len=*) , optional , intent(in)    :: fldname

     ! local variables
     integer :: nu
     logical :: exists  ! file existance
     character(*), parameter :: F00   = "('(dshr_restart_read) ',8a)"
     character(*), parameter :: subName = "(dshr_restart_read) "
     !-------------------------------------------------------------------------------

     if (trim(rest_file) == trim(nullstr) .and. trim(rest_file_strm) == trim(nullstr)) then
        if (my_task == master_task) then
           write(logunit,F00) ' restart filenames from rpointer'
           inquire(file=trim(rpfile)//trim(inst_suffix), exist=exists)
           if (.not.exists) then
              write(logunit, F00) ' ERROR: rpointer file does not exist'
              call shr_sys_abort(trim(subname)//' ERROR: rpointer file missing')
           endif
           open(newunit=nu, file=trim(rpfile)//trim(inst_suffix), form='formatted')
           read(nu, '(a)') rest_file
           read(nu, '(a)') rest_file_strm
           close(nu)
           inquire(file=trim(rest_file_strm), exist=exists)
        endif
        call shr_mpi_bcast(rest_file, mpicom, 'rest_file')
        call shr_mpi_bcast(rest_file_strm, mpicom, 'rest_file_strm')
     else
        ! use namelist already read
        if (my_task == master_task) then
           write(logunit, F00) ' restart filenames from namelist '
           inquire(file=trim(rest_file_strm), exist=exists)
        endif
     endif
     call shr_mpi_bcast(exists, mpicom, 'exists')
     if (exists) then
        if (my_task == master_task) write(logunit, F00) ' reading ', trim(rest_file_strm)
        if (present(fld) .and. present(fldname)) then
           call shr_pcdf_readwrite('read', sdat%pio_subsystem, sdat%io_type, trim(rest_file), &
                mpicom, sdat%gsmap, clobber=.true., rf1=fld, rf1n=trim(fldname), io_format=sdat%io_format)
        end if
        call shr_strdata_restRead(trim(rest_file_strm), sdat, mpicom)
     else
        if (my_task == master_task) write(logunit, F00) ' file not found, skipping ',trim(rest_file_strm)
     endif
   end subroutine dshr_restart_read

   !===============================================================================

   subroutine dshr_restart_write(rpfile, case_name, model_name, inst_suffix, ymd, tod, &
        logunit, mpicom, my_task, master_task, sdat, fld, fldname) 

     use shr_pcdf_mod    , only : shr_pcdf_readwrite
     use shr_cal_mod     , only : shr_cal_datetod2string
     use shr_strdata_mod , only : shr_strdata_restWrite, shr_strdata_type

     ! input/output variables
     character(len=*)            , intent(in)    :: rpfile
     character(len=*)            , intent(in)    :: case_name
     character(len=*)            , intent(in)    :: model_name
     character(len=*)            , intent(in)    :: inst_suffix
     integer                     , intent(in)    :: ymd       ! model date
     integer                     , intent(in)    :: tod       ! model sec into model date
     integer                     , intent(in)    :: logunit
     integer                     , intent(in)    :: my_task
     integer                     , intent(in)    :: master_task
     integer                     , intent(in)    :: mpicom
     type(shr_strdata_type)      , intent(inout) :: sdat
     real(r8)         , optional , pointer       :: fld(:)
     character(len=*) , optional , intent(in)    :: fldname

     ! local variables
     character(len=CL) :: rest_file
     character(len=CL) :: rest_file_strm
     character(len=CS) :: date_str
     integer           :: nu
     character(*), parameter  :: F00   = "('(dshr_restart_write) ',a,2f10.4)"
     !-------------------------------------------------------------------------------

     call shr_cal_datetod2string(date_str, ymd, tod)
     write(rest_file,"(7a)")      trim(case_name),'.', trim(model_name),trim(inst_suffix),'.r.'  , trim(date_str),'.nc'
     write(rest_file_strm,"(7a)") trim(case_name),'.', trim(model_name),trim(inst_suffix),'.rs1.', trim(date_str),'.bin'
     if (my_task == master_task) then
        open(newunit=nu, file=trim(rpfile)//trim(inst_suffix), form='formatted')
        write(nu,'(a)') rest_file
        write(nu,'(a)') rest_file_strm
        close(nu)
     endif
     if (my_task == master_task) then
        write(logunit,F00) ' writing ',trim(rest_file_strm), ymd, tod
     end if
     if (present(fld) .and. present(fldname)) then
        call shr_pcdf_readwrite('write', sdat%pio_subsystem, sdat%io_type,&
             trim(rest_file), mpicom, sdat%gsmap, clobber=.true., rf1=fld, rf1n=trim(fldname))
     end if
     call shr_strdata_restWrite(trim(rest_file_strm), sdat, mpicom, trim(case_name), 'SDAT strdata from '//trim(model_name))

   end subroutine dshr_restart_write

end module dshr_nuopc_mod
