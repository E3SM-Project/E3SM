module shr_nuopc_methods_mod

  !-----------------------------------------------------------------------------
  ! Generic operation methods used by the Mediator Component.
  !-----------------------------------------------------------------------------

  use ESMF               , only : operator(<), operator(/=), operator(+), operator(-), operator(*) , operator(>=)
  use ESMF               , only : operator(<=), operator(>), operator(==)
  use ESMF               , only : ESMF_GeomType_Flag, ESMF_FieldStatus_Flag, ESMF_PoleMethod_Flag
  use ESMF               , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF               , only : ESMF_LOGERR_PASSTHRU, ESMF_LogFoundError, ESMF_LOGMSG_ERROR
  use ESMF               , only : ESMF_MAXSTR, ESMF_LOGMSG_WARNING, ESMF_POLEMETHOD_ALLAVG
  use med_constants_mod  , only : dbug_flag => med_constants_dbug_flag

  use shr_nuopc_utils_mod, only : shr_nuopc_methods_ChkErr => shr_nuopc_utils_ChkErr, shr_nuopc_utils_ChkErr
  implicit none
  private

  interface shr_nuopc_methods_FB_accum ; module procedure &
    shr_nuopc_methods_FB_accumFB2FB, &
    shr_nuopc_methods_FB_accumFB2ST, &
    shr_nuopc_methods_FB_accumST2FB
  end interface

  interface shr_nuopc_methods_FB_copy ; module procedure &
    shr_nuopc_methods_FB_copyFB2FB, &
    shr_nuopc_methods_FB_copyFB2ST, &
    shr_nuopc_methods_FB_copyST2FB
  end interface

  interface shr_nuopc_methods_FieldPtr_compare ; module procedure &
    shr_nuopc_methods_FieldPtr_compare1, &
    shr_nuopc_methods_FieldPtr_compare2
  end interface

  ! tcraig, interfaces cannot be differentiated, revisit later
  !  interface shr_nuopc_methods_FB_FieldMerge ; module procedure &
  !    shr_nuopc_methods_FB_FieldMerge_1D, &
  !    shr_nuopc_methods_FB_FieldMerge_2D
  !  end interface

  interface shr_nuopc_methods_UpdateTimestamp; module procedure &
    shr_nuopc_methods_State_UpdateTimestamp, &
    shr_nuopc_methods_Field_UpdateTimestamp
  end interface

  ! used/reused in module

  logical                               :: isPresent
  logical                               :: rhprint_flag = .false.
  integer                               :: srcTermProcessing_Value = 0
  character(len=1024)                   :: msgString
  type(ESMF_GeomType_Flag)              :: geomtype
  type(ESMF_FieldStatus_Flag)           :: status
  type(ESMF_PoleMethod_Flag), parameter :: polemethod=ESMF_POLEMETHOD_ALLAVG
  character(*)      , parameter         :: u_FILE_u = &
       __FILE__

  public shr_nuopc_methods_FB_copy
  public shr_nuopc_methods_FB_accum
  public shr_nuopc_methods_FB_average
  public shr_nuopc_methods_FB_init
  public shr_nuopc_methods_FB_reset
  public shr_nuopc_methods_FB_clean
  public shr_nuopc_methods_FB_diagnose
  public shr_nuopc_methods_FB_Regrid
  public shr_nuopc_methods_FB_FldChk
  public shr_nuopc_methods_FB_FieldRegrid
  public shr_nuopc_methods_FB_FieldMerge_1D
  public shr_nuopc_methods_FB_FieldMerge_2D
  public shr_nuopc_methods_FB_GetFldPtr
  public shr_nuopc_methods_FB_getNameN
  public shr_nuopc_methods_FB_getFieldN
  public shr_nuopc_methods_FB_Field_diagnose
  public shr_nuopc_methods_State_reset
  public shr_nuopc_methods_State_diagnose
  public shr_nuopc_methods_State_GeomPrint
  public shr_nuopc_methods_State_GeomWrite
  public shr_nuopc_methods_State_GetFldPtr
  public shr_nuopc_methods_State_SetScalar
  public shr_nuopc_methods_State_GetScalar
  public shr_nuopc_methods_State_GetNumFields
  public shr_nuopc_methods_State_getFieldN
  public shr_nuopc_methods_Field_GeomPrint
  public shr_nuopc_methods_Clock_TimePrint
  public shr_nuopc_methods_UpdateTimestamp
  public shr_nuopc_methods_ChkErr
  public shr_nuopc_methods_Distgrid_Match
  public shr_nuopc_methods_Print_FieldExchInfo

  private shr_nuopc_methods_Grid_Write
  private shr_nuopc_methods_Grid_Print
  private shr_nuopc_methods_Mesh_Print
  private shr_nuopc_methods_Mesh_Write
  private shr_nuopc_methods_Field_GetFldPtr
  private shr_nuopc_methods_Field_GeomWrite
  private shr_nuopc_methods_Field_UpdateTimestamp
  private shr_nuopc_methods_FB_GeomPrint
  private shr_nuopc_methods_FB_GeomWrite
  private shr_nuopc_methods_FB_RWFields
  private shr_nuopc_methods_FB_getFieldByName
  private shr_nuopc_methods_FB_FieldCopy
  private shr_nuopc_methods_FB_SetFldPtr
  private shr_nuopc_methods_FB_copyFB2FB
  private shr_nuopc_methods_FB_copyFB2ST
  private shr_nuopc_methods_FB_copyST2FB
  private shr_nuopc_methods_FB_accumFB2FB
  private shr_nuopc_methods_FB_accumST2FB
  private shr_nuopc_methods_FB_accumFB2ST
  private shr_nuopc_methods_FieldPtr_compare1
  private shr_nuopc_methods_FieldPtr_compare2
  private shr_nuopc_methods_State_UpdateTimestamp
  private shr_nuopc_methods_State_getNameN
  private shr_nuopc_methods_State_getFieldByName
  private shr_nuopc_methods_State_SetFldPtr
  private shr_nuopc_methods_Array_diagnose

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_RWFields(mode,fname,FB,flag,rc)
    ! ----------------------------------------------
    ! Read or Write Field Bundles
    ! ----------------------------------------------
    use ESMF, only : ESMF_Field, ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_FieldBundleWrite
    use ESMF, only : ESMF_FieldRead, ESMF_IOFMT_NETCDF, ESMF_FILESTATUS_REPLACE

    character(len=*)       :: mode
    character(len=*)       :: fname
    type(ESMF_FieldBundle) :: FB
    logical,optional       :: flag
    integer,optional       :: rc

    ! local variables
    type(ESMF_Field)           :: field
    character(len=ESMF_MAXSTR) :: name
    integer                    :: fieldcount, n
    logical                    :: fexists
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_FB_RWFields)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//trim(fname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    if (mode == 'write') then
      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//": write "//trim(fname), ESMF_LOGMSG_INFO, rc=dbrc)
      end if
      call ESMF_FieldBundleWrite(FB, fname, &
        singleFile=.true., status=ESMF_FILESTATUS_REPLACE, iofmt=ESMF_IOFMT_NETCDF, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_diagnose(FB, 'write '//trim(fname), rc)

    elseif (mode == 'read') then
      inquire(file=fname,exist=fexists)
      if (fexists) then
        if (dbug_flag > 5) then
           call ESMF_LogWrite(trim(subname)//": read "//trim(fname), ESMF_LOGMSG_INFO, rc=dbrc)
        end if
        !-----------------------------------------------------------------------------------------------------
        ! tcraig, ESMF_FieldBundleRead fails if a field is not on the field bundle, but we really want to just
        ! ignore that field and read the rest, so instead read each field one at a time through ESMF_FieldRead
        !        call ESMF_FieldBundleRead (FB, fname, &
        !          singleFile=.true., iofmt=ESMF_IOFMT_NETCDF, rc=rc)
        !        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        !-----------------------------------------------------------------------------------------------------
        call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        do n = 1,fieldCount
          call shr_nuopc_methods_FB_getFieldByName(FB, name, field, rc)
          if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldRead (field, fname, iofmt=ESMF_IOFMT_NETCDF, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=u_FILE_u)) call ESMF_LogWrite(trim(subname)//' WARNING missing field '//trim(name),rc=dbrc)
        enddo

        call shr_nuopc_methods_FB_diagnose(FB, 'read '//trim(fname), rc)
	if (present(flag)) flag = .true.
      endif

    else
      call ESMF_LogWrite(trim(subname)//": mode WARNING "//trim(fname)//" mode="//trim(mode), ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//trim(fname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_RWFields

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_init(FBout, flds_scalar_name, fieldNameList, FBgeom, STgeom, FBflds, STflds, name, rc)
    use ESMF              , only : ESMF_Field, ESMF_FieldBundle, ESMF_FieldBundleCreate, ESMF_FieldBundleGet
    use ESMF              , only : ESMF_State, ESMF_Grid, ESMF_Mesh, ESMF_StaggerLoc, ESMF_MeshLoc
    use ESMF              , only : ESMF_StateGet, ESMF_FieldGet, ESMF_FieldBundleAdd, ESMF_FieldCreate
    use ESMF              , only : ESMF_TYPEKIND_R8, ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID
    use ESMF              , only : ESMF_FIELDSTATUS_EMPTY
    use med_constants_mod , only : spval_init => med_constants_spval_init

    ! ----------------------------------------------
    ! Create FBout from fieldNameList, FBflds, STflds, FBgeom or STgeom in that order or priority
    ! Pass in FBgeom OR STgeom, get grid/mesh from that object
    ! ----------------------------------------------
    type(ESMF_FieldBundle), intent(inout) :: FBout
    character(len=*)      , intent(in)    :: flds_scalar_name
    character(len=*)      , intent(in), optional :: fieldNameList(:)
    type(ESMF_FieldBundle), intent(in), optional :: FBgeom
    type(ESMF_State)      , intent(in), optional :: STgeom
    type(ESMF_FieldBundle), intent(in), optional :: FBflds
    type(ESMF_State)      , intent(in), optional :: STflds
    character(len=*)      , intent(in), optional :: name
    integer               , intent(out) :: rc

    ! local variables
    integer                    :: i,j,n,n1
    integer                    :: fieldCount,fieldCountgeom
    logical                    :: found
    character(ESMF_MAXSTR)     :: lname
    character(ESMF_MAXSTR),allocatable :: lfieldNameList(:)
    type(ESMF_Field)           :: field,lfield
    type(ESMF_Grid)            :: lgrid
    type(ESMF_Mesh)            :: lmesh
    type(ESMF_StaggerLoc)      :: staggerloc
    type(ESMF_MeshLoc)         :: meshloc
    integer :: dbrc
    character(len=*),parameter :: subname='(shr_nuopc_methods_FB_init)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    lname = 'undefined'
    if (present(name)) then
      lname = trim(name)
    endif
    lname = 'FB '//trim(lname)

    !---------------------------------
    ! check argument consistency and
    ! verify that geom argument has a field
    !---------------------------------

    if (present(fieldNameList) .and. present(FBflds) .and. present(STflds)) then
      call ESMF_LogWrite(trim(subname)//": ERROR only fieldNameList, FBflds, or STflds can be an argument", &
            ESMF_LOGMSG_INFO, rc=rc)
      rc = ESMF_FAILURE
      return
    endif

    if (present(FBgeom) .and. present(STgeom)) then
       call ESMF_LogWrite(trim(subname)//": ERROR FBgeom and STgeom cannot both be arguments", &
            ESMF_LOGMSG_INFO, rc=rc)
      rc = ESMF_FAILURE
      return
    endif

    if (.not.present(FBgeom) .and. .not.present(STgeom)) then
       call ESMF_LogWrite(trim(subname)//": ERROR FBgeom or STgeom must be an argument", &
            ESMF_LOGMSG_INFO, rc=rc)
      rc = ESMF_FAILURE
      return
    endif

    if (present(FBgeom)) then
      call ESMF_FieldBundleGet(FBgeom, fieldCount=fieldCountGeom, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    elseif (present(STgeom)) then
      call ESMF_StateGet(STgeom, itemCount=fieldCountGeom, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//": ERROR FBgeom or STgeom must be passed", ESMF_LOGMSG_INFO, rc=rc)
      rc = ESMF_FAILURE
      return
    endif

    !---------------------------------
    ! determine the names of fields that will be in FBout
    !---------------------------------

    if (present(fieldNameList)) then
      fieldcount = size(fieldNameList)
      allocate(lfieldNameList(fieldcount))
      lfieldNameList = fieldNameList
      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" fieldNameList from argument", ESMF_LOGMSG_INFO, rc=rc)
      end if
    elseif (present(FBflds)) then
      call ESMF_FieldBundleGet(FBflds, fieldCount=fieldCount, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      allocate(lfieldNameList(fieldCount))
      call ESMF_FieldBundleGet(FBflds, fieldNameList=lfieldNameList, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" fieldNameList from FBflds", ESMF_LOGMSG_INFO, rc=rc)
      end if
    elseif (present(STflds)) then
      call ESMF_StateGet(STflds, itemCount=fieldCount, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      allocate(lfieldNameList(fieldCount))
      call ESMF_StateGet(STflds, itemNameList=lfieldNameList, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" fieldNameList from STflds", ESMF_LOGMSG_INFO, rc=rc)
      end if
    elseif (present(FBgeom)) then
      call ESMF_FieldBundleGet(FBgeom, fieldCount=fieldCount, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      allocate(lfieldNameList(fieldCount))
      call ESMF_FieldBundleGet(FBgeom, fieldNameList=lfieldNameList, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" fieldNameList from FBgeom", ESMF_LOGMSG_INFO, rc=rc)
      end if
    elseif (present(STgeom)) then
      call ESMF_StateGet(STgeom, itemCount=fieldCount, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      allocate(lfieldNameList(fieldCount))
      call ESMF_StateGet(STgeom, itemNameList=lfieldNameList, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (dbug_flag > 5) then
         call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" fieldNameList from STflds", ESMF_LOGMSG_INFO, rc=rc)
      end if
    else
       call ESMF_LogWrite(trim(subname)//": ERROR fieldNameList, FBflds, STflds, FBgeom, or STgeom must be passed", &
            ESMF_LOGMSG_INFO, rc=rc)
      rc = ESMF_FAILURE
      return
    endif

    !---------------------------------
    ! remove scalar field and blank fields from field bundle
    !---------------------------------

    do n = 1, fieldCount
      if (trim(lfieldnamelist(n)) == trim(flds_scalar_name) .or. &
          trim(lfieldnamelist(n)) == '') then
        do n1 = n, fieldCount-1
          lfieldnamelist(n1) = lfieldnamelist(n1+1)
        enddo
        fieldCount = fieldCount - 1
      endif
    enddo  ! n

    !---------------------------------
    ! create the grid (lgrid) or mesh(lmesh)
    ! that will be used for FBout fields
    !---------------------------------

    if (fieldcount > 0 .and. fieldcountgeom > 0) then

      ! Look at only the first field in either the FBgeom and STgeom to get the grid
      if (present(FBgeom)) then
        call shr_nuopc_methods_FB_getFieldN(FBgeom, 1, lfield, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 5) then
           call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" grid/mesh from FBgeom", ESMF_LOGMSG_INFO, rc=rc)
        end if
      elseif (present(STgeom)) then
        call shr_nuopc_methods_State_getFieldN(STgeom, 1, lfield, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 5) then
           call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" grid/mesh from STgeom", ESMF_LOGMSG_INFO, rc=rc)
        end if
      else
        call ESMF_LogWrite(trim(subname)//": ERROR FBgeom or STgeom must be passed", ESMF_LOGMSG_INFO, rc=rc)
        rc = ESMF_FAILURE
        return
      endif

      ! Make sure the field is not empty - if it is return with an error
      call ESMF_FieldGet(lfield, status=status, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (status == ESMF_FIELDSTATUS_EMPTY) then
         call ESMF_LogWrite(trim(subname)//":"//trim(lname)//": ERROR field does not have a geom yet ", &
              ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
        rc = ESMF_FAILURE
        return
      endif

      ! Determine if first field in either FBgeom or STgeom is on a grid or a mesh
      call ESMF_FieldGet(lfield, geomtype=geomtype, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      if (geomtype == ESMF_GEOMTYPE_GRID) then
        call ESMF_FieldGet(lfield, grid=lgrid, staggerloc=staggerloc, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 5) then
           call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" use grid", ESMF_LOGMSG_INFO, rc=rc)
        end if
      elseif (geomtype == ESMF_GEOMTYPE_MESH) then
        call ESMF_FieldGet(lfield, mesh=lmesh, meshloc=meshloc, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 5) then
           call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" use mesh", ESMF_LOGMSG_INFO, rc=rc)
        end if
      else  ! geomtype
        call ESMF_LogWrite(trim(subname)//": ERROR geomtype not supported ", ESMF_LOGMSG_INFO, rc=rc)
        rc = ESMF_FAILURE
        return
      endif ! geomtype

    endif  ! fieldcount > 0

    !---------------------------------
    ! create FBout
    !---------------------------------

    FBout = ESMF_FieldBundleCreate(name=trim(lname), rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (fieldcountgeom > 0) then

      ! Now loop over all the fields in either FBgeom or STgeom
      do n = 1, fieldCount

         ! Create the field on either lgrid or lmesh
        if (geomtype == ESMF_GEOMTYPE_GRID) then
          field = ESMF_FieldCreate(lgrid, ESMF_TYPEKIND_R8, staggerloc=staggerloc, name=lfieldNameList(n), rc=rc)
          if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        elseif (geomtype == ESMF_GEOMTYPE_MESH) then
          field = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_R8, meshloc=meshloc, name=lfieldNameList(n), rc=rc)
          if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        else  ! geomtype
          call ESMF_LogWrite(trim(subname)//": ERROR no grid/mesh for field ", ESMF_LOGMSG_INFO, rc=rc)
          rc = ESMF_FAILURE
          return
        endif

        ! Add the created field bundle FBout
        call ESMF_FieldBundleAdd(FBout, (/field/), rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (dbug_flag > 1) then
          call ESMF_LogWrite(trim(subname)//":"//trim(lname)//" adding field "//trim(lfieldNameList(n)), &
               ESMF_LOGMSG_INFO, rc=dbrc)
        endif

      enddo  ! fieldCount

    endif  ! fieldcountgeom

    deallocate(lfieldNameList)

    call shr_nuopc_methods_FB_reset(FBout, value=spval_init, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_init

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_getNameN(FB, fieldnum, fieldname, rc)
    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet
    ! ----------------------------------------------
    ! Get name of field number fieldnum in FB
    ! ----------------------------------------------
    type(ESMF_FieldBundle), intent(in)    :: FB
    integer               , intent(in)    :: fieldnum
    character(len=*)      , intent(out)   :: fieldname
    integer               , intent(out)   :: rc

    ! local variables
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    integer                         :: dbrc
    character(len=*),parameter      :: subname='(shr_nuopc_methods_FB_getNameN)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    fieldname = ' '

    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (fieldnum > fieldCount) then
      call ESMF_LogWrite(trim(subname)//": ERROR fieldnum > fieldCount ", ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    fieldname = lfieldnamelist(fieldnum)

    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_getNameN

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_getFieldN(FB, fieldnum, field, rc)
    use ESMF, only : ESMF_Field, ESMF_FieldBundle, ESMF_FieldBundleGet

    ! ----------------------------------------------
    ! Get field number fieldnum out of FB
    ! ----------------------------------------------
    type(ESMF_FieldBundle), intent(in)    :: FB
    integer               , intent(in)    :: fieldnum
    type(ESMF_Field)      , intent(inout) :: field
    integer               , intent(out)   :: rc

    ! local variables
    character(len=ESMF_MAXSTR) :: name
    integer                    :: dbrc
    character(len=*),parameter :: subname='(shr_nuopc_methods_FB_getFieldN)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call shr_nuopc_methods_FB_getNameN(FB, fieldnum, name, rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(FB, fieldName=name, field=field, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_getFieldN

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_getFieldByName(FB, fieldname, field, rc)
    use ESMF, only : ESMF_Field, ESMF_FieldBundle, ESMF_FieldBundleGet
    ! ----------------------------------------------
    ! Get field associated with fieldname out of FB
    ! ----------------------------------------------
    type(ESMF_FieldBundle), intent(in)    :: FB
    character(len=*)      , intent(in)    :: fieldname
    type(ESMF_Field)      , intent(inout) :: field
    integer               , intent(out)   :: rc

    ! local variables
    integer :: dbrc
    character(len=*),parameter :: subname='(shr_nuopc_methods_FB_getFieldByName)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(FB, fieldName=fieldname, field=field, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_getFieldByName

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_State_getNameN(State, fieldnum, fieldname, rc)
    use ESMF, only : ESMF_State, ESMF_StateGet
    ! ----------------------------------------------
    ! Get field number fieldnum name out of State
    ! ----------------------------------------------
    type(ESMF_State), intent(in)    :: State
    integer         , intent(in)    :: fieldnum
    character(len=*), intent(out)   :: fieldname
    integer         , intent(out)   :: rc

    ! local variables
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    integer                         :: dbrc
    character(len=*),parameter      :: subname='(shr_nuopc_methods_State_getNameN)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    fieldname = ' '

    call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (fieldnum > fieldCount) then
      call ESMF_LogWrite(trim(subname)//": ERROR fieldnum > fieldCount ", ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    allocate(lfieldnamelist(fieldCount))
    call ESMF_StateGet(State, itemNameList=lfieldnamelist, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    fieldname = lfieldnamelist(fieldnum)

    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_State_getNameN

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_State_getNumFields(State, fieldnum, rc)
    use NUOPC , only : NUOPC_GetStateMemberLists
    use ESMF  , only : ESMF_State, ESMF_Field, ESMF_StateGet, ESMF_STATEITEM_FIELD
    use ESMF  , only : ESMF_StateItem_Flag
    ! ----------------------------------------------
    ! Get field number fieldnum name out of State
    ! ----------------------------------------------
    type(ESMF_State), intent(in)    :: State
    integer         , intent(inout) :: fieldnum
    integer         , intent(out)   :: rc

    ! local variables
    integer                            :: n,itemCount
    type(ESMF_Field), pointer          :: fieldList(:)
    type(ESMF_StateItem_Flag), pointer :: itemTypeList(:)
    logical, parameter                 :: use_NUOPC_method = .true.
    integer                            :: dbrc
    character(len=*),parameter         :: subname='(shr_nuopc_methods_State_getNumFields)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    if (use_NUOPC_method) then

      nullify(fieldList)
      call NUOPC_GetStateMemberLists(state, fieldList=fieldList, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      fieldnum = 0
      if (associated(fieldList)) then
        fieldnum = size(fieldList)
        deallocate(fieldList)
      endif

    else

      fieldnum = 0
      call ESMF_StateGet(State, itemCount=itemCount, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      if (itemCount > 0) then
        allocate(itemTypeList(itemCount))
        call ESMF_StateGet(State, itemTypeList=itemTypeList, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

        do n = 1,itemCount
          if (itemTypeList(n) == ESMF_STATEITEM_FIELD) fieldnum=fieldnum+1
        enddo
        deallocate(itemTypeList)
      endif

    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_State_getNumFields

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_State_getFieldN(State, fieldnum, field, rc)
    use ESMF, only : ESMF_State, ESMF_Field, ESMF_StateGet
    ! ----------------------------------------------
    ! Get field number fieldnum in State
    ! ----------------------------------------------
    type(ESMF_State), intent(in)    :: State
    integer         , intent(in)    :: fieldnum
    type(ESMF_Field), intent(inout) :: field
    integer         , intent(out)   :: rc

    ! local variables
    character(len=ESMF_MAXSTR) :: name
    integer :: dbrc
    character(len=*),parameter :: subname='(shr_nuopc_methods_State_getFieldN)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call shr_nuopc_methods_State_getNameN(State, fieldnum, name, rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(State, itemName=name, field=field, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_State_getFieldN

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_State_getFieldByName(State, fieldname, field, rc)
    ! ----------------------------------------------
    ! Get field associated with fieldname from State
    ! ----------------------------------------------
    use ESMF, only : ESMF_State, ESMF_Field, ESMF_StateGet

    type(ESMF_State), intent(in)    :: State
    character(len=*), intent(in)    :: fieldname
    type(ESMF_Field), intent(inout) :: field
    integer         , intent(out)   :: rc

    ! local variables
    integer :: dbrc
    character(len=*),parameter :: subname='(shr_nuopc_methods_State_getFieldByName)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=fieldname, field=field, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_State_getFieldByName

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_clean(FB, rc)
    ! ----------------------------------------------
    ! Destroy fields in FB and FB
    ! ----------------------------------------------
    use med_constants_mod, only : R8
    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_FieldDestroy
    use ESMF, only : ESMF_FieldBundleDestroy, ESMF_Field

    type(ESMF_FieldBundle), intent(inout) :: FB
    integer               , intent(out)   :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    type(ESMF_Field)                :: field
    integer :: dbrc
    character(len=*),parameter      :: subname='(shr_nuopc_methods_FB_clean)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1, fieldCount
      call ESMF_FieldBundleGet(FB, fieldName=lfieldnamelist(n), field=field, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldDestroy(field, rc=rc, noGarbage=.true.)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo
    call ESMF_FieldBundleDestroy(FB, rc=rc, noGarbage=.true.)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_clean

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_reset(FB, value, rc)
    ! ----------------------------------------------
    ! Set all fields to value in FB
    ! If value is not provided, reset to 0.0
    ! ----------------------------------------------
    use med_constants_mod , only : czero => med_constants_czero
    use med_constants_mod , only : R8
    use ESMF              , only : ESMF_FieldBundle, ESMF_FieldBundleGet

    type(ESMF_FieldBundle), intent(inout)        :: FB
    real(R8)    , intent(in), optional :: value
    integer               , intent(out)          :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    real(R8)              :: lvalue
    integer :: dbrc
    character(len=*),parameter      :: subname='(shr_nuopc_methods_FB_reset)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    lvalue = czero
    if (present(value)) then
      lvalue = value
    endif

    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount
      call shr_nuopc_methods_FB_SetFldPtr(FB, lfieldnamelist(n), lvalue, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_reset

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_FieldCopy(FBin,fldin,FBout,fldout,rc)
    use med_constants_mod, only : R8
    use ESMF, only : ESMF_FieldBundle, ESMF_FAILURE, ESMF_LOGMSG_ERROR
    ! ----------------------------------------------
    ! Copy a field in a field bundle to another field in a field bundle
    ! ----------------------------------------------
    type(ESMF_FieldBundle), intent(inout) :: FBin
    character(len=*)      , intent(in)    :: fldin
    type(ESMF_FieldBundle), intent(inout) :: FBout
    character(len=*)      , intent(in)    :: fldout
    integer               , intent(out)   :: rc

    ! local
    real(R8), pointer :: dataPtrIn1(:)
    real(R8), pointer :: dataPtrOut1(:)
    real(R8), pointer :: dataPtrIn2(:,:)
    real(R8), pointer :: dataPtrOut2(:,:)
    integer                     :: lrankIn, lrankOut
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_FB_FieldCopy)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    if (shr_nuopc_methods_FB_FldChk(FBin , trim(fldin) , rc=rc) .and. &
        shr_nuopc_methods_FB_FldChk(FBout, trim(fldout), rc=rc)) then

      call shr_nuopc_methods_FB_GetFldPtr(FBin, trim(fldin), dataPtrIn1, dataPtrIn2, lrankIn, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_GetFldPtr(FBout, trim(fldout), dataPtrOut1, dataPtrOut2, lrankOut, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      if (lrankIn /= lrankOut) then
        call ESMF_LogWrite(trim(subname)//": ERROR FBin and FBout different rank", ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
        rc = ESMF_FAILURE
        return
      endif

      if (lrankIn == 1 .and. lrankOut == 1) then
        if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtrIn1, dataPtrOut1, subname, rc)) then
          call ESMF_LogWrite(trim(subname)//": ERROR data1d different sizes", ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return
        endif

        dataPtrOut1 = dataPtrIn1

      elseif (lrankIn == 2 .and. lrankOut == 2) then
        if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtrIn2, dataPtrOut2, subname, rc)) then
          call ESMF_LogWrite(trim(subname)//": ERROR data2d different sizes", ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return
        endif

        dataPtrOut2 = dataPtrIn2

      elseif (lrankIn == 2 .and. lrankOut == 2) then
        call ESMF_LogWrite(trim(subname)//": ERROR raskIn and rankOut invalid", ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
        rc = ESMF_FAILURE
        return

      endif

    else

       if (dbug_flag > 1) then
         call ESMF_LogWrite(trim(subname)//" field not found: "//trim(fldin)//","//trim(fldout), ESMF_LOGMSG_INFO, rc=dbrc)
       endif

    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_FieldCopy

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_Regrid(shortnames, flds_scalar_name, mappings, &
       FBin, FBout, &
       consfmap, consdmap, bilnrmap, patchmap, &
       fcopymap, string, rc)

    use ESMF, only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_FAILURE
    use ESMF, only : ESMF_FieldBundle, ESMF_RouteHandle, ESMF_RouteHandleIsCreated

    character(len=*)       , intent(in)               :: shortnames(:)
    character(len=*)       , intent(in)               :: mappings(:)
    character(len=*)       , intent(in)               :: flds_scalar_name
    type(ESMF_FieldBundle) , intent(inout)            :: FBin
    type(ESMF_FieldBundle) , intent(inout)            :: FBout
    type(ESMF_Routehandle) , intent(inout) , optional :: consfmap
    type(ESMF_Routehandle) , intent(inout) , optional :: consdmap
    type(ESMF_Routehandle) , intent(inout) , optional :: bilnrmap
    type(ESMF_Routehandle) , intent(inout) , optional :: patchmap
    type(ESMF_Routehandle) , intent(inout) , optional :: fcopymap
    character(len=*)       , intent(in)    , optional :: string
    integer                , intent(out)   , optional :: rc

    ! local variables
    integer           :: n
    character(len=64) :: lstring
    logical           :: okconsf, okconsd, okbilnr, okpatch, okfcopy
    integer :: dbrc
    character(len=*),parameter :: subname='(shr_nuopc_methods_FB_Regrid)'
    ! ----------------------------------------------

    if (present(string)) then
      lstring = trim(string)
    else
      lstring = " "
    endif

    if (.not.present(rc)) then
      call ESMF_LogWrite(trim(subname)//trim(lstring)//": ERROR rc expected", ESMF_LOGMSG_INFO, rc=rc)
      rc = ESMF_FAILURE
      return
    endif

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//trim(lstring)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    okconsf = .false.
    if (present(consfmap)) then
      if (ESMF_RouteHandleIsCreated(consfmap, rc=rc)) okconsf = .true.
    endif

    okconsd = .false.
    if (present(consdmap)) then
      if (ESMF_RouteHandleIsCreated(consdmap, rc=rc)) okconsd = .true.
    endif

    okbilnr = .false.
    if (present(bilnrmap)) then
      if (ESMF_RouteHandleIsCreated(bilnrmap, rc=rc)) okbilnr = .true.
    endif

    okpatch = .false.
    if (present(patchmap)) then
      if (ESMF_RouteHandleIsCreated(patchmap, rc=rc)) okpatch = .true.
    endif

    okfcopy = .false.
    if (present(fcopymap)) then
      if (ESMF_RouteHandleIsCreated(fcopymap, rc=rc)) okfcopy = .true.
    endif

    ! Loop over all fields in the field bundle - determine the mapping for the target field from
    ! the mappings setting and apply it
    do n = 1,size(shortnames)

      if (shortnames(n) == trim(flds_scalar_name)) then
        if (dbug_flag > 1) then
           call ESMF_LogWrite(trim(subname)//trim(lstring)//": skip : fld="//trim(shortnames(n)), &
                ESMF_LOGMSG_INFO, rc=dbrc)
        endif

      elseif (shr_nuopc_methods_FB_FldChk(FBin , shortnames(n), rc=rc) .and. &
              shr_nuopc_methods_FB_FldChk(FBout, shortnames(n), rc=rc)) then

        if (dbug_flag > 1) then
          call ESMF_LogWrite(trim(subname)//trim(lstring)//": map="//trim(mappings(n))// &
            ": fld="//trim(shortnames(n)), ESMF_LOGMSG_INFO, rc=dbrc)
        endif

        if (mappings(n) == 'bilinear') then
          if (.not. okbilnr) then
            call ESMF_LogWrite(trim(subname)//trim(lstring)//": ERROR RH not available for "//trim(mappings(n))// &
              ": fld="//trim(shortnames(n)), ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
            rc = ESMF_FAILURE
            return
          endif
          call shr_nuopc_methods_FB_FieldRegrid(FBin, shortnames(n), FBout, shortnames(n), bilnrmap,rc)
          if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

        elseif (mappings(n) == "conservefrac") then
          if (.not. okconsf) then
            call ESMF_LogWrite(trim(subname)//trim(lstring)//": ERROR RH not available for "//trim(mappings(n))// &
              ": fld="//trim(shortnames(n)), ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
            rc = ESMF_FAILURE
            return
          endif
          call shr_nuopc_methods_FB_FieldRegrid(FBin, shortnames(n), FBout,shortnames(n), consfmap, rc)
          if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

        elseif (mappings(n) == "conservedst") then
          if (.not. okconsd) then
            call ESMF_LogWrite(trim(subname)//trim(lstring)//": ERROR RH not available for "//trim(mappings(n))// &
              ": fld="//trim(shortnames(n)), ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
            rc = ESMF_FAILURE
            return
          endif
          call shr_nuopc_methods_FB_FieldRegrid(FBin, shortnames(n), FBout,shortnames(n), consdmap, rc)
          if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

        elseif (mappings(n) == 'patch') then
          if (.not. okpatch) then
            call ESMF_LogWrite(trim(subname)//trim(lstring)//": ERROR RH not available for "//trim(mappings(n))// &
              ": fld="//trim(shortnames(n)), ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
            rc = ESMF_FAILURE
            return
          endif
          call shr_nuopc_methods_FB_FieldRegrid(FBin, shortnames(n), FBout,shortnames(n), patchmap,rc)
          if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

        elseif (mappings(n) == 'copy') then
          !-------------------------------------------
          ! copy will not exist for some grid combinations
          ! so fall back to conservative frac as a secondary option
          !-------------------------------------------
          if (.not. okfcopy) then
            if (.not. okconsf) then
              call ESMF_LogWrite(trim(subname)//trim(lstring)//": ERROR RH not available for "//trim(mappings(n))// &
                ": fld="//trim(shortnames(n)), ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
              rc = ESMF_FAILURE
              return
            else
              if (dbug_flag > 1) then
                 call ESMF_LogWrite(trim(subname)//trim(lstring)//": NOTE using conservative instead of copy for"// &
                      " fld="//trim(shortnames(n)), ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
              end if
              call shr_nuopc_methods_FB_FieldRegrid(FBin ,shortnames(n), FBout, shortnames(n), consfmap,rc)
              if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
            endif
          else
            call shr_nuopc_methods_FB_FieldRegrid(FBin ,shortnames(n), FBout,shortnames(n), fcopymap,rc)
            if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif

        else
          call ESMF_LogWrite(trim(subname)//trim(lstring)//": ERROR unrecognized mapping "//trim(mappings(n))// &
            ": fld="//trim(shortnames(n)), ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return
        endif

      else
        if (dbug_flag > 1) then
          call ESMF_LogWrite(trim(subname)//" field not found in FB: "//trim(shortnames(n)), ESMF_LOGMSG_INFO, rc=dbrc)
        endif
      endif
    enddo

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//trim(lstring)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_Regrid

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_FieldRegrid(FBin,fldin,FBout,fldout,RH,rc,debug)

    ! ----------------------------------------------
    ! Regrid a field in a field bundle to another field in a field bundle
    ! ----------------------------------------------
    use ESMF, only : ESMF_FieldBundle, ESMF_RouteHandle, ESMF_FieldRegrid, ESMF_Field
    use ESMF, only : ESMF_TERMORDER_SRCSEQ, ESMF_FieldRegridStore, ESMF_SparseMatrixWrite
    use med_constants_mod, only : R8
    type(ESMF_FieldBundle), intent(inout) :: FBin
    character(len=*)      , intent(in)    :: fldin
    type(ESMF_FieldBundle), intent(inout) :: FBout
    character(len=*)      , intent(in)    :: fldout
    type(ESMF_RouteHandle), intent(inout) :: RH
    integer               , intent(out)   :: rc
    logical, intent(in), optional :: debug
    ! local
    real(R8),             pointer :: factorList(:)
    integer,          pointer :: factorIndexList(:,:)
    type(ESMF_Field) :: field1, field2
    integer :: dbrc
    character(len=*),parameter :: subname='(shr_nuopc_methods_FB_FieldRegrid)'
    integer :: rank
    character(len=8) :: filename
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    if (shr_nuopc_methods_FB_FldChk(FBin , trim(fldin) , rc=rc) .and. &
        shr_nuopc_methods_FB_FldChk(FBout, trim(fldout), rc=rc)) then

      call shr_nuopc_methods_FB_getFieldByName(FBin, trim(fldin), field1, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_FB_getFieldByName(FBout, trim(fldout), field2, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_FieldRegrid(field1, field2, routehandle=RH, &
        termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=.true., rc=rc)

      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
   else

      if (dbug_flag > 1) then
        call ESMF_LogWrite(trim(subname)//" field not found: "//trim(fldin)//","//trim(fldout), ESMF_LOGMSG_INFO, rc=dbrc)
      endif

    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_FieldRegrid

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_FieldMerge_2D(FBout, fnameout, &
                                    FBinA, fnameA, wgtA, &
                                    FBinB, fnameB, wgtB, &
                                    FBinC, fnameC, wgtC, &
                                    FBinD, fnameD, wgtD, &
                                    FBinE, fnameE, wgtE, rc)

    ! ----------------------------------------------
    ! Supports up to a five way merge
    ! ----------------------------------------------
    use med_constants_mod , only : R8
    use ESMF              , only : ESMF_FieldBundle
    use med_constants_mod , only : czero => med_constants_czero

    type(ESMF_FieldBundle), intent(inout)        :: FBout
    character(len=*)      , intent(in)           :: fnameout
    type(ESMF_FieldBundle), intent(in), optional :: FBinA
    character(len=*)      , intent(in), optional :: fnameA
    real(R8)    , intent(in), optional, pointer  :: wgtA(:,:)
    type(ESMF_FieldBundle), intent(in), optional :: FBinB
    character(len=*)      , intent(in), optional :: fnameB
    real(R8)    , intent(in), optional, pointer  :: wgtB(:,:)
    type(ESMF_FieldBundle), intent(in), optional :: FBinC
    character(len=*)      , intent(in), optional :: fnameC
    real(R8)    , intent(in), optional, pointer  :: wgtC(:,:)
    type(ESMF_FieldBundle), intent(in), optional :: FBinD
    character(len=*)      , intent(in), optional :: fnameD
    real(R8)    , intent(in), optional, pointer  :: wgtD(:,:)
    type(ESMF_FieldBundle), intent(in), optional :: FBinE
    character(len=*)      , intent(in), optional :: fnameE
    real(R8)    , intent(in), optional, pointer  :: wgtE(:,:)
    integer               , intent(out)          :: rc

    ! local variables
    real(R8), pointer          :: dataOut(:,:)
    real(R8), pointer          :: dataPtr(:,:)
    real(R8), pointer          :: wgt(:,:)
    integer                    :: lb1,ub1,lb2,ub2,i,j,n
    logical                    :: wgtfound, FBinfound
    integer                    :: dbrc
    character(len=*),parameter :: subname='(shr_nuopc_methods_FB_FieldMerge_2D)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc=ESMF_SUCCESS

    if (.not. shr_nuopc_methods_FB_FldChk(FBout, trim(fnameout), rc=rc)) then
      call ESMF_LogWrite(trim(subname)//": WARNING field not in FBout, skipping merge "//&
           trim(fnameout), ESMF_LOGMSG_WARNING, line=__LINE__, file=u_FILE_u, rc=dbrc)
      return
    endif

    call shr_nuopc_methods_FB_GetFldPtr(FBout, trim(fnameout), fldptr2=dataOut, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    lb1 = lbound(dataOut,1)
    ub1 = ubound(dataOut,1)
    lb2 = lbound(dataOut,2)
    ub2 = ubound(dataOut,2)
    allocate(wgt(lb1:ub1,lb2:ub2))

    dataOut = czero

    ! check each field has a fieldname passed in
    if ((present(FBinA) .and. .not.present(fnameA)) .or. &
        (present(FBinB) .and. .not.present(fnameB)) .or. &
        (present(FBinC) .and. .not.present(fnameC)) .or. &
        (present(FBinD) .and. .not.present(fnameD)) .or. &
        (present(FBinE) .and. .not.present(fnameE))) then
      call ESMF_LogWrite(trim(subname)//": ERROR fname not present with FBin", &
           ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    ! check that each field passed in actually exists, if not DO NOT do any merge
    FBinfound = .true.
    if (present(FBinA)) then
      if (.not. shr_nuopc_methods_FB_FldChk(FBinA, trim(fnameA), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinB)) then
      if (.not. shr_nuopc_methods_FB_FldChk(FBinB, trim(fnameB), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinC)) then
      if (.not. shr_nuopc_methods_FB_FldChk(FBinC, trim(fnameC), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinD)) then
      if (.not. shr_nuopc_methods_FB_FldChk(FBinD, trim(fnameD), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinE)) then
      if (.not. shr_nuopc_methods_FB_FldChk(FBinE, trim(fnameE), rc=rc)) FBinfound = .false.
    endif
    if (.not. FBinfound) then
      call ESMF_LogWrite(trim(subname)//": WARNING field not found in FBin, skipping merge "//trim(fnameout), &
           ESMF_LOGMSG_WARNING, line=__LINE__, file=u_FILE_u, rc=dbrc)
      return
    endif

    ! n=1,5 represents adding A to E inputs if they exist
    do n = 1,5
      FBinfound = .false.
      wgtfound = .false.

      if (n == 1 .and. present(FBinA)) then
        FBinfound = .true.
        call shr_nuopc_methods_FB_GetFldPtr(FBinA, trim(fnameA), fldptr2=dataPtr, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (present(wgtA)) then
          wgtfound = .true.
          wgt => wgtA
        endif

      elseif (n == 2 .and. present(FBinB)) then
        FBinfound = .true.
        call shr_nuopc_methods_FB_GetFldPtr(FBinB, trim(fnameB), fldptr2=dataPtr, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (present(wgtB)) then
          wgtfound = .true.
          wgt => wgtB
        endif

      elseif (n == 3 .and. present(FBinC)) then
        FBinfound = .true.
        call shr_nuopc_methods_FB_GetFldPtr(FBinC, trim(fnameC), fldptr2=dataPtr, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (present(wgtC)) then
          wgtfound = .true.
          wgt => wgtC
        endif

      elseif (n == 4 .and. present(FBinD)) then
        FBinfound = .true.
        call shr_nuopc_methods_FB_GetFldPtr(FBinD, trim(fnameD), fldptr2=dataPtr, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (present(wgtD)) then
          wgtfound = .true.
          wgt => wgtD
        endif

      elseif (n == 5 .and. present(FBinE)) then
        FBinfound = .true.
        call shr_nuopc_methods_FB_GetFldPtr(FBinE, trim(fnameE), fldptr2=dataPtr, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (present(wgtE)) then
          wgtfound = .true.
          wgt => wgtE
        endif

      endif

      if (FBinfound) then
        if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtr, dataOut, subname, rc)) then
          call ESMF_LogWrite(trim(subname)//": ERROR FBin wrong size", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return
        endif

        if (wgtfound) then
          if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtr, wgt, subname, rc)) then
            call ESMF_LogWrite(trim(subname)//": ERROR wgt wrong size", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
            rc = ESMF_FAILURE
            return
          endif
          do j = lb2,ub2
          do i = lb1,ub1
            dataOut(i,j) = dataOut(i,j) + dataPtr(i,j) * wgt(i,j)
          enddo
          enddo
        else
          do j = lb2,ub2
          do i = lb1,ub1
            dataOut(i,j) = dataOut(i,j) + dataPtr(i,j)
          enddo
          enddo
        endif  ! wgtfound

      endif  ! FBin found
    enddo  ! n

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_FieldMerge_2D

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_FieldMerge_1D(FBout, fnameout, &
                                    FBinA, fnameA, wgtA, &
                                    FBinB, fnameB, wgtB, &
                                    FBinC, fnameC, wgtC, &
                                    FBinD, fnameD, wgtD, &
                                    FBinE, fnameE, wgtE, rc)
    ! ----------------------------------------------
    ! Supports up to a five way merge
    ! ----------------------------------------------
    use med_constants_mod , only : R8
    use med_constants_mod , only : czero => med_constants_czero
    use ESMF              , only : ESMF_FieldBundle

    type(ESMF_FieldBundle), intent(inout)                 :: FBout
    character(len=*)      , intent(in)                    :: fnameout
    type(ESMF_FieldBundle), intent(in), optional          :: FBinA
    character(len=*)      , intent(in), optional          :: fnameA
    real(R8)    , intent(in), optional, pointer :: wgtA(:)
    type(ESMF_FieldBundle), intent(in), optional          :: FBinB
    character(len=*)      , intent(in), optional          :: fnameB
    real(R8)    , intent(in), optional, pointer :: wgtB(:)
    type(ESMF_FieldBundle), intent(in), optional          :: FBinC
    character(len=*)      , intent(in), optional          :: fnameC
    real(R8)    , intent(in), optional, pointer :: wgtC(:)
    type(ESMF_FieldBundle), intent(in), optional          :: FBinD
    character(len=*)      , intent(in), optional          :: fnameD
    real(R8)    , intent(in), optional, pointer :: wgtD(:)
    type(ESMF_FieldBundle), intent(in), optional          :: FBinE
    character(len=*)      , intent(in), optional          :: fnameE
    real(R8)    , intent(in), optional, pointer :: wgtE(:)
    integer               , intent(out)                   :: rc

    ! local variables
    real(R8), pointer :: dataOut(:)
    real(R8), pointer :: dataPtr(:)
    real(R8), pointer :: wgt(:)
    integer                     :: lb1,ub1,i,j,n
    logical                     :: wgtfound, FBinfound
    integer :: dbrc
    character(len=*),parameter  :: subname='(shr_nuopc_methods_FB_FieldMerge_1D)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc=ESMF_SUCCESS

    if (.not. shr_nuopc_methods_FB_FldChk(FBout, trim(fnameout), rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": WARNING field not in FBout, skipping merge "//trim(fnameout), &
            ESMF_LOGMSG_WARNING, line=__LINE__, file=u_FILE_u, rc=dbrc)
      return
    endif

    call shr_nuopc_methods_FB_GetFldPtr(FBout, trim(fnameout), fldptr1=dataOut, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    lb1 = lbound(dataOut,1)
    ub1 = ubound(dataOut,1)
    allocate(wgt(lb1:ub1))

    dataOut = czero

    ! check each field has a fieldname passed in
    if ((present(FBinA) .and. .not.present(fnameA)) .or. &
        (present(FBinB) .and. .not.present(fnameB)) .or. &
        (present(FBinC) .and. .not.present(fnameC)) .or. &
        (present(FBinD) .and. .not.present(fnameD)) .or. &
        (present(FBinE) .and. .not.present(fnameE))) then
       call ESMF_LogWrite(trim(subname)//": ERROR fname not present with FBin", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    ! check that each field passed in actually exists, if not DO NOT do any merge
    FBinfound = .true.
    if (present(FBinA)) then
      if (.not. shr_nuopc_methods_FB_FldChk(FBinA, trim(fnameA), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinB)) then
      if (.not. shr_nuopc_methods_FB_FldChk(FBinB, trim(fnameB), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinC)) then
      if (.not. shr_nuopc_methods_FB_FldChk(FBinC, trim(fnameC), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinD)) then
      if (.not. shr_nuopc_methods_FB_FldChk(FBinD, trim(fnameD), rc=rc)) FBinfound = .false.
    endif
    if (present(FBinE)) then
      if (.not. shr_nuopc_methods_FB_FldChk(FBinE, trim(fnameE), rc=rc)) FBinfound = .false.
    endif
    if (.not. FBinfound) then
      call ESMF_LogWrite(trim(subname)//": WARNING field not found in FBin, skipping merge "//trim(fnameout), &
           ESMF_LOGMSG_WARNING, line=__LINE__, file=u_FILE_u, rc=dbrc)
      return
    endif

    ! n=1,5 represents adding A to E inputs if they exist
    do n = 1,5
      FBinfound = .false.
      wgtfound = .false.

      if (n == 1 .and. present(FBinA)) then
        FBinfound = .true.
        call shr_nuopc_methods_FB_GetFldPtr(FBinA, trim(fnameA), fldptr1=dataPtr, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (present(wgtA)) then
          wgtfound = .true.
          wgt => wgtA
        endif

      elseif (n == 2 .and. present(FBinB)) then
        FBinfound = .true.
        call shr_nuopc_methods_FB_GetFldPtr(FBinB, trim(fnameB), fldptr1=dataPtr, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (present(wgtB)) then
          wgtfound = .true.
          wgt => wgtB
        endif

      elseif (n == 3 .and. present(FBinC)) then
        FBinfound = .true.
        call shr_nuopc_methods_FB_GetFldPtr(FBinC, trim(fnameC), fldptr1=dataPtr, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (present(wgtC)) then
          wgtfound = .true.
          wgt => wgtC
        endif

      elseif (n == 4 .and. present(FBinD)) then
        FBinfound = .true.
        call shr_nuopc_methods_FB_GetFldPtr(FBinD, trim(fnameD), fldptr1=dataPtr, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (present(wgtD)) then
          wgtfound = .true.
          wgt => wgtD
        endif

      elseif (n == 5 .and. present(FBinE)) then
        FBinfound = .true.
        call shr_nuopc_methods_FB_GetFldPtr(FBinE, trim(fnameE), fldptr1=dataPtr, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (present(wgtE)) then
          wgtfound = .true.
          wgt => wgtE
        endif

      endif

      if (FBinfound) then
        if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtr, dataOut, subname, rc)) then
          call ESMF_LogWrite(trim(subname)//": ERROR FBin wrong size", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return
        endif

        if (wgtfound) then
          if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtr, wgt, subname, rc)) then
            call ESMF_LogWrite(trim(subname)//": ERROR wgt wrong size", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
            rc = ESMF_FAILURE
            return
          endif
          do i = lb1,ub1
            dataOut(i) = dataOut(i) + dataPtr(i) * wgt(i)
          enddo
        else
          do i = lb1,ub1
            dataOut(i) = dataOut(i) + dataPtr(i)
          enddo
        endif  ! wgtfound

      endif  ! FBin found
    enddo  ! n

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_FieldMerge_1D

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_State_reset(State, value, rc)

    ! ----------------------------------------------
    ! Set all fields to value in State
    ! If value is not provided, reset to 0.0
    ! ----------------------------------------------
    use med_constants_mod , only : czero => med_constants_czero
    use med_constants_mod , only : R8
    use ESMF              , only : ESMF_State, ESMF_StateGet

    type(ESMF_State) , intent(inout)        :: State
    real(R8)         , intent(in), optional :: value
    integer          , intent(out)          :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    real(R8)              :: lvalue
    integer :: dbrc
    character(len=*),parameter      :: subname='(shr_nuopc_methods_State_reset)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    lvalue = czero
    if (present(value)) then
      lvalue = value
    endif

    call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))
    call ESMF_StateGet(State, itemNameList=lfieldnamelist, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount
      call shr_nuopc_methods_State_SetFldPtr(State, lfieldnamelist(n), lvalue, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_State_reset

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_average(FB, count, rc)
    ! ----------------------------------------------
    ! Set all fields to zero in FB
    ! ----------------------------------------------
    use med_constants_mod , only : R8
    use ESMF              , only : ESMF_FieldBundle, ESMF_FieldBundleGet

    type(ESMF_FieldBundle), intent(inout) :: FB
    integer               , intent(in)    :: count
    integer               , intent(out)   :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    real(R8), pointer               :: dataPtr1(:)
    real(R8), pointer               :: dataPtr2(:,:)
    integer                         :: dbrc
    character(len=*),parameter      :: subname='(shr_nuopc_methods_FB_average)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    if (count == 0) then

       if (dbug_flag > 10) then
          call ESMF_LogWrite(trim(subname)//": WARNING count is 0", ESMF_LOGMSG_INFO, rc=dbrc)
       end if
       !call ESMF_LogWrite(trim(subname)//": WARNING count is 0 set avg to spval", ESMF_LOGMSG_INFO, rc=dbrc)
       !call shr_nuopc_methods_FB_reset(FB, value=spval, rc=rc)
       !if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    else

      call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      allocate(lfieldnamelist(fieldCount))
      call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      do n = 1, fieldCount
        call shr_nuopc_methods_FB_GetFldPtr(FB, lfieldnamelist(n), dataPtr1, dataPtr2, lrank, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

        if (lrank == 0) then
          ! no local data
        elseif (lrank == 1) then
          do i=lbound(dataptr1,1),ubound(dataptr1,1)
            dataptr1(i) = dataptr1(i) / real(count, R8)
          enddo
        elseif (lrank == 2) then
          do j=lbound(dataptr2,2),ubound(dataptr2,2)
          do i=lbound(dataptr2,1),ubound(dataptr2,1)
            dataptr2(i,j) = dataptr2(i,j) / real(count, R8)
          enddo
          enddo
        else
          call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return
        endif
      enddo
      deallocate(lfieldnamelist)

    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_average

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_diagnose(FB, string, rc)
    ! ----------------------------------------------
    ! Diagnose status of FB
    ! ----------------------------------------------
    use med_constants_mod , only : R8
    use ESMF              , only : ESMF_FieldBundle, ESMF_FieldBundleGet

    type(ESMF_FieldBundle) , intent(inout)        :: FB
    character(len=*)       , intent(in), optional :: string
    integer                , intent(out)          :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR), pointer :: lfieldnamelist(:)
    character(len=64)               :: lstring
    real(R8), pointer               :: dataPtr1d(:)
    real(R8), pointer               :: dataPtr2d(:,:)
    integer                         :: dbrc
    character(len=*), parameter     :: subname='(shr_nuopc_methods_FB_diagnose)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    lstring = ''
    if (present(string)) then
       lstring = trim(string) // ' '
    endif

    ! Determine number of fields in field bundle and allocate memory for lfieldnamelist
    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    ! Get the fields in the field bundle
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! For each field in the bundle, get its memory location and print out the field
    do n = 1, fieldCount
       call shr_nuopc_methods_FB_GetFldPtr(FB, lfieldnamelist(n), &
            fldptr1=dataPtr1d, fldptr2=dataPtr2d, rank=lrank, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data

       elseif (lrank == 1) then
          if (size(dataPtr1d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  " no data"
          endif

       elseif (lrank == 2) then
          if (size(dataPtr2d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  " no data"
          endif

       else
          call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return
       endif
       if (dbug_flag > 1) then
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
       end if
    enddo

    ! Deallocate memory
    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_diagnose

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_Array_diagnose(array, string, rc)
    ! ----------------------------------------------
    ! Diagnose status of Array
    ! ----------------------------------------------
    use med_constants_mod, only : R8
    use ESMF, only : ESMF_Array, ESMF_ArrayGet

    type(ESMF_Array), intent(inout)        :: array
    character(len=*), intent(in), optional :: string
    integer         , intent(out)          :: rc

    ! local variables
    character(len=64)           :: lstring
    real(R8), pointer :: dataPtr3d(:,:,:)
    integer :: dbrc
    character(len=*),parameter  :: subname='(shr_nuopc_methods_Array_diagnose)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! this is not working yet, not sure about dataPtr dim/type
    return

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call ESMF_ArrayGet(Array, farrayPtr=dataPtr3d, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    write(msgString,'(A,3g14.7)') trim(subname)//' '//trim(lstring), &
        minval(dataPtr3d), maxval(dataPtr3d), sum(dataPtr3d)

    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
    end if

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_Array_diagnose

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_State_diagnose(State, string, rc)
    ! ----------------------------------------------
    ! Diagnose status of State
    ! ----------------------------------------------
    use med_constants_mod, only : R8
    use ESMF, only : ESMF_State, ESMF_StateGet

    type(ESMF_State), intent(in)           :: State
    character(len=*), intent(in), optional :: string
    integer         , intent(out)          :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    character(len=64)               :: lstring
    real(R8), pointer     :: dataPtr1d(:)
    real(R8), pointer     :: dataPtr2d(:,:)
    integer :: dbrc
    character(len=*),parameter      :: subname='(shr_nuopc_methods_State_diagnose)'
    ! ----------------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    call ESMF_StateGet(State, itemNameList=lfieldnamelist, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount
       call shr_nuopc_methods_State_GetFldPtr(State, lfieldnamelist(n), &
            fldptr1=dataPtr1d, fldptr2=dataPtr2d, rank=lrank, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          if (lrank == 0) then
             ! no local data

          elseif (lrank == 1) then
             if (size(dataPtr1d) > 0) then
                write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                     minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
             else
                write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                     " no data"
             endif

          elseif (lrank == 2) then
             if (size(dataPtr2d) > 0) then
                write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                     minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
             else
                write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                     " no data"
             endif

          else
             call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR, line=__LINE__, &
                  file=u_FILE_u, rc=dbrc)
             rc = ESMF_FAILURE
             return
          endif

          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
       end if
    enddo

    deallocate(lfieldnamelist)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
   endif

  end subroutine shr_nuopc_methods_State_diagnose

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_Field_diagnose(FB, fieldname, string, rc)
    ! ----------------------------------------------
    ! Diagnose status of State
    ! ----------------------------------------------
    use med_constants_mod, only : R8
    use ESMF, only : ESMF_FieldBundle

    type(ESMF_FieldBundle), intent(inout)  :: FB
    character(len=*), intent(in)           :: fieldname
    character(len=*), intent(in), optional :: string
    integer         , intent(out)          :: rc

    ! local variables
    integer                         :: lrank
    character(len=64)               :: lstring
    real(R8), pointer     :: dataPtr1d(:)
    real(R8), pointer     :: dataPtr2d(:,:)
    integer :: dbrc
    character(len=*),parameter      :: subname='(shr_nuopc_methods_FB_FieldDiagnose)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call shr_nuopc_methods_FB_GetFldPtr(FB, fieldname, dataPtr1d, dataPtr2d, lrank, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (lrank == 0) then
       ! no local data
    elseif (lrank == 1) then
       if (size(dataPtr1d) > 0) then
          write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
       else
          write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               " no data"
       endif
    elseif (lrank == 2) then
       if (size(dataPtr2d) > 0) then
          write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
       else
          write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               " no data"
       endif
    else
       call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR, line=__LINE__, &
            file=u_FILE_u, rc=dbrc)
       rc = ESMF_FAILURE
       return
    endif
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
    end if

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_Field_diagnose

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_copyFB2FB(FBout, FBin, rc)
    ! ----------------------------------------------
    ! Copy common field names from FBin to FBout
    ! ----------------------------------------------
    use ESMF, only : ESMF_FieldBundle

    type(ESMF_FieldBundle), intent(inout) :: FBout
    type(ESMF_FieldBundle), intent(in)    :: FBin
    integer               , intent(out)   :: rc

    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_FB_copyFB2FB)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call shr_nuopc_methods_FB_accum(FBout, FBin, copy=.true., rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_copyFB2FB

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_copyFB2ST(STout, FBin, rc)
    ! ----------------------------------------------
    ! Copy common field names from FBin to STout
    ! ----------------------------------------------
    use ESMF, only : ESMF_State, ESMF_FieldBundle

    type(ESMF_State)      , intent(inout) :: STout
    type(ESMF_FieldBundle), intent(in)    :: FBin
    integer               , intent(out)   :: rc

    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_FB_copyFB2ST)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call shr_nuopc_methods_FB_accum(STout, FBin, copy=.true., rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_copyFB2ST

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_copyST2FB(FBout, STin, rc)
    ! ----------------------------------------------
    ! Copy common field names from STin to FBout
    ! ----------------------------------------------
    use ESMF, only : ESMF_State, ESMF_FieldBundle

    type(ESMF_FieldBundle), intent(inout) :: FBout
    type(ESMF_State)      , intent(in)    :: STin
    integer               , intent(out)   :: rc
    integer :: dbrc
    character(len=*),parameter :: subname='(shr_nuopc_methods_FB_copyST2FB)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call shr_nuopc_methods_FB_accum(FBout, STin, copy=.true., rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_copyST2FB

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_accumFB2FB(FBout, FBin, copy, rc)
    ! ----------------------------------------------
    ! Accumulate common field names from FBin to FBout
    ! If copy is passed in and true, the this is a copy
    ! ----------------------------------------------
    use med_constants_mod , only : R8
    use ESMF              , only : ESMF_FieldBundle
    use ESMF              , only : ESMF_FieldBundleGet

    type(ESMF_FieldBundle), intent(inout) :: FBout
    type(ESMF_FieldBundle), intent(in)    :: FBin
    logical, optional     , intent(in)    :: copy
    integer               , intent(out)   :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount, lranki, lranko
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    logical                         :: exists
    logical                         :: lcopy
    real(R8), pointer               :: dataPtri1(:)  , dataPtro1(:)
    real(R8), pointer               :: dataPtri2(:,:), dataPtro2(:,:)
    integer                         :: dbrc
    character(len=*), parameter     :: subname='(shr_nuopc_methods_FB_accumFB2FB)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    lcopy = .false.  ! accumulate by default
    if (present(copy)) then
      lcopy = copy
    endif

    call ESMF_FieldBundleGet(FBout, fieldCount=fieldCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FBout, fieldNameList=lfieldnamelist, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount
      call ESMF_FieldBundleGet(FBin, fieldName=lfieldnamelist(n), isPresent=exists, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (exists) then
        call shr_nuopc_methods_FB_GetFldPtr(FBin,  lfieldnamelist(n), dataPtri1, dataPtri2, lranki, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_GetFldPtr(FBout, lfieldnamelist(n), dataPtro1, dataPtro2, lranko, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

        if (lranki == 1 .and. lranko == 1) then

          if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtro1, dataPtri1, subname, rc)) then
            call ESMF_LogWrite(trim(subname)//": ERROR in dataPtr1 size ", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=rc)
            rc = ESMF_FAILURE
            return
          endif

          if (lcopy) then
            do i=lbound(dataPtri1,1),ubound(dataPtri1,1)
              dataPtro1(i) = dataPtri1(i)
            enddo
          else
            do i=lbound(dataPtri1,1),ubound(dataPtri1,1)
              dataPtro1(i) = dataPtro1(i) + dataPtri1(i)
            enddo
          endif

        elseif (lranki == 2 .and. lranko == 2) then

          if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtro2, dataPtri2, subname, rc)) then
            call ESMF_LogWrite(trim(subname)//": ERROR in dataPtr2 size ", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=rc)
            rc = ESMF_FAILURE
            return
          endif

          if (lcopy) then
            do j=lbound(dataPtri2,2),ubound(dataPtri2,2)
            do i=lbound(dataPtri2,1),ubound(dataPtri2,1)
              dataPtro2(i,j) = dataPtri2(i,j)
            enddo
            enddo
          else
            do j=lbound(dataPtri2,2),ubound(dataPtri2,2)
            do i=lbound(dataPtri2,1),ubound(dataPtri2,1)
              dataPtro2(i,j) = dataPtro2(i,j) + dataPtri2(i,j)
            enddo
            enddo
          endif

        else

          write(msgString,'(a,2i8)') trim(subname)//": ranki, ranko = ",lranki,lranko
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
          call ESMF_LogWrite(trim(subname)//": ERROR ranki ranko not supported "//trim(lfieldnamelist(n)), &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return

        endif

      endif
    enddo

    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_accumFB2FB
  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_accumST2FB(FBout, STin, copy, rc)
    ! ----------------------------------------------
    ! Accumulate common field names from State to FB
    ! If copy is passed in and true, the this is a copy
    ! ----------------------------------------------
    use med_constants_mod, only : R8
    use ESMF, only : ESMF_State, ESMF_FieldBundle
    use ESMF, only : ESMF_StateGet, ESMF_FieldBundleGet
    use ESMF, only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag

    type(ESMF_FieldBundle), intent(inout) :: FBout
    type(ESMF_State)      , intent(in)    :: STin
    logical, optional     , intent(in)    :: copy
    integer               , intent(out)   :: rc

    ! local variables
    integer                     :: i,j,n
    integer                     :: fieldCount, lrankS, lrankB
    logical                     :: lcopy
    character(ESMF_MAXSTR) ,pointer  :: lfieldnamelist(:)
    type(ESMF_StateItem_Flag)   :: itemType
    real(R8), pointer :: dataPtrS1(:)  , dataPtrB1(:)
    real(R8), pointer :: dataPtrS2(:,:), dataPtrB2(:,:)
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_FB_accumST2FB)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    lcopy = .false.
    if (present(copy)) then
      lcopy = copy
    endif

    call ESMF_FieldBundleGet(FBout, fieldCount=fieldCount, rc=rc)
    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FBout, fieldNameList=lfieldnamelist, rc=rc)
    do n = 1, fieldCount
      call ESMF_StateGet(STin, itemName=lfieldnamelist(n), itemType=itemType, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (itemType /= ESMF_STATEITEM_NOTFOUND) then

        call shr_nuopc_methods_State_GetFldPtr(STin, lfieldnamelist(n), dataPtrS1, dataPtrS2, lrankS, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_FB_GetFldPtr(FBout, lfieldnamelist(n), dataPtrB1, dataPtrB2, lrankB, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

        if (lrankB == 0 .and. lrankS == 0) then

          ! no local data

        elseif (lrankS == 1 .and. lrankB == 1) then

          if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtrS1, dataPtrB1, subname, rc)) then
            call ESMF_LogWrite(trim(subname)//": ERROR in dataPtr1 size ", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=rc)
            rc = ESMF_FAILURE
            return
          endif

          if (lcopy) then
            do i=lbound(dataPtrB1,1),ubound(dataPtrB1,1)
              dataPtrB1(i) = dataPtrS1(i)
            enddo
          else
            do i=lbound(dataPtrB1,1),ubound(dataPtrB1,1)
              dataPtrB1(i) = dataPtrB1(i) + dataPtrS1(i)
            enddo
          endif

        elseif (lrankS == 2 .and. lrankB == 2) then

          if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtrS2, dataPtrB2, subname, rc)) then
            call ESMF_LogWrite(trim(subname)//": ERROR in dataPtr2 size ", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=rc)
            rc = ESMF_FAILURE
            return
          endif

          if (lcopy) then
            do j=lbound(dataPtrB2,2),ubound(dataPtrB2,2)
            do i=lbound(dataPtrB2,1),ubound(dataPtrB2,1)
              dataPtrB2(i,j) = dataPtrS2(i,j)
            enddo
            enddo
          else
            do j=lbound(dataPtrB2,2),ubound(dataPtrB2,2)
            do i=lbound(dataPtrB2,1),ubound(dataPtrB2,1)
              dataPtrB2(i,j) = dataPtrB2(i,j) + dataPtrS2(i,j)
            enddo
            enddo
          endif

        else

          write(msgString,'(a,2i8)') trim(subname)//": rankB, ranks = ",lrankB,lrankS
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
          call ESMF_LogWrite(trim(subname)//": ERROR rankB rankS not supported "//trim(lfieldnamelist(n)), &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return

        endif

      endif  ! statefound
    enddo  ! fieldCount

    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_accumST2FB

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_accumFB2ST(STout, FBin, copy, rc)
    ! ----------------------------------------------
    ! Accumulate common field names from FB to State
    ! If copy is passed in and true, the this is a copy
    ! ----------------------------------------------
    use med_constants_mod , only : R8
    use ESMF              , only : ESMF_State, ESMF_FieldBundle
    use ESMF              , only : ESMF_StateGet, ESMF_FieldBundleGet
    use ESMF              , only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag

    type(ESMF_State)      , intent(inout) :: STout
    type(ESMF_FieldBundle), intent(in)    :: FBin
    logical, optional     , intent(in)    :: copy
    integer               , intent(out)   :: rc

    ! local variables
    integer                     :: i,j,n
    integer                     :: fieldCount, lrankS, lrankB
    logical                     :: lcopy
    character(ESMF_MAXSTR) ,pointer  :: lfieldnamelist(:)
    type(ESMF_StateItem_Flag)   :: itemType
    real(R8), pointer :: dataPtrS1(:), dataPtrB1(:)
    real(R8), pointer :: dataPtrS2(:,:), dataPtrB2(:,:)
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_FB_accumFB2ST)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    lcopy = .false.
    if (present(copy)) then
      lcopy = copy
    endif

    call ESMF_FieldBundleGet(FBin, fieldCount=fieldCount, rc=rc)
    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FBin, fieldNameList=lfieldnamelist, rc=rc)
    do n = 1, fieldCount
      call ESMF_StateGet(STout, itemName=lfieldnamelist(n), itemType=itemType, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (itemType /= ESMF_STATEITEM_NOTFOUND) then

        call shr_nuopc_methods_FB_GetFldPtr(FBin, lfieldnamelist(n), dataPtrB1, dataPtrB2, lrankB, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        call shr_nuopc_methods_State_GetFldPtr(STout, lfieldnamelist(n), dataPtrS1, dataPtrS2, lrankS, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

        if (lrankB == 0 .and. lrankS == 0) then

          ! no local data

        elseif (lrankB == 1 .and. lrankS == 1) then

          if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtrS1, dataPtrB1, subname, rc)) then
            call ESMF_LogWrite(trim(subname)//": ERROR in dataPtr1 size ", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=rc)
            rc = ESMF_FAILURE
            return
          endif

          if (lcopy) then
            do i=lbound(dataPtrB1,1),ubound(dataPtrB1,1)
              dataPtrS1(i) = dataPtrB1(i)
            enddo
          else
            do i=lbound(dataPtrB1,1),ubound(dataPtrB1,1)
              dataPtrS1(i) = dataPtrS1(i) + dataPtrB1(i)
            enddo
          endif

        elseif (lrankB == 2 .and. lrankS == 2) then

          if (.not.shr_nuopc_methods_FieldPtr_Compare(dataPtrS2, dataPtrB2, subname, rc)) then
            call ESMF_LogWrite(trim(subname)//": ERROR in dataPtr2 size ", &
                 ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=rc)
            rc = ESMF_FAILURE
            return
          endif

          if (lcopy) then
            do j=lbound(dataPtrB2,2),ubound(dataPtrB2,2)
            do i=lbound(dataPtrB2,1),ubound(dataPtrB2,1)
              dataPtrS2(i,j) = dataPtrB2(i,j)
            enddo
            enddo
          else
            do j=lbound(dataPtrB2,2),ubound(dataPtrB2,2)
            do i=lbound(dataPtrB2,1),ubound(dataPtrB2,1)
              dataPtrS2(i,j) = dataPtrS2(i,j) + dataPtrB2(i,j)
            enddo
            enddo
          endif

        else

          write(msgString,'(a,2i8)') trim(subname)//": rankB, ranks = ",lrankB,lrankS
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
          call ESMF_LogWrite(trim(subname)//": ERROR rankB rankS not supported "//trim(lfieldnamelist(n)), &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return

        endif

      endif  ! statefound
    enddo  ! fieldCount

    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_accumFB2ST

  !-----------------------------------------------------------------------------

  logical function shr_nuopc_methods_FB_FldChk(FB, fldname, rc)
    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet

    type(ESMF_FieldBundle), intent(in)  :: FB
    character(len=*)      , intent(in)  :: fldname
    integer               , intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_FB_FldChk)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    shr_nuopc_methods_FB_FldChk = .false.

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), isPresent=isPresent, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       shr_nuopc_methods_FB_FldChk = .true.
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end function shr_nuopc_methods_FB_FldChk

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_Field_GetFldPtr(field, fldptr1, fldptr2, rank, abort, rc)
    ! ----------------------------------------------
    ! for a field, determine rank and return fldptr1 or fldptr2
    ! abort is true by default and will abort if fldptr is not yet allocated in field
    ! rank returns 0, 1, or 2.  0 means fldptr not allocated and abort=false
    ! ----------------------------------------------
    use med_constants_mod , only : R8
    use ESMF              , only : ESMF_Field,ESMF_Mesh, ESMF_FieldGet, ESMF_MeshGet
    use ESMF              , only : ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID, ESMF_FIELDSTATUS_COMPLETE

    type(ESMF_Field)  , intent(in)              :: field
    real(R8), pointer , intent(inout), optional :: fldptr1(:)
    real(R8), pointer , intent(inout), optional :: fldptr2(:,:)
    integer           , intent(out)  , optional :: rank
    logical           , intent(in)   , optional :: abort
    integer           , intent(out)  , optional :: rc

    ! local variables
    type(ESMF_Mesh) :: lmesh
    integer         :: lrank, nnodes, nelements
    logical         :: labort
    integer         :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_Field_GetFldPtr)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    if (.not.present(rc)) then
       call ESMF_LogWrite(trim(subname)//": ERROR rc not present ", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    rc = ESMF_SUCCESS

    labort = .true.
    if (present(abort)) then
      labort = abort
    endif
    lrank = -99

    call ESMF_FieldGet(field, status=status, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
      lrank = 0
      if (labort) then
        call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
        rc = ESMF_FAILURE
        return
      else
        call ESMF_LogWrite(trim(subname)//": WARNING data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
      endif
    else

      call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      if (geomtype == ESMF_GEOMTYPE_GRID) then
        call ESMF_FieldGet(field, rank=lrank, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      elseif (geomtype == ESMF_GEOMTYPE_MESH) then
        lrank = 1
        call ESMF_FieldGet(field, mesh=lmesh, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (nnodes == 0 .and. nelements == 0) lrank = 0
      else  ! geomtype
         call ESMF_LogWrite(trim(subname)//": ERROR geomtype not supported ", &
              ESMF_LOGMSG_INFO, rc=rc)
        rc = ESMF_FAILURE
        return
      endif ! geomtype

      if (lrank == 0) then
         call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", &
              ESMF_LOGMSG_INFO, rc=dbrc)
      elseif (lrank == 1) then
        if (.not.present(fldptr1)) then
           call ESMF_LogWrite(trim(subname)//": ERROR missing rank=1 array ", &
                ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return
        endif
        call ESMF_FieldGet(field, farrayPtr=fldptr1, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      elseif (lrank == 2) then
        if (.not.present(fldptr2)) then
           call ESMF_LogWrite(trim(subname)//": ERROR missing rank=2 array ", &
                ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return
        endif
        call ESMF_FieldGet(field, farrayPtr=fldptr2, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      else
         call ESMF_LogWrite(trim(subname)//": ERROR in rank ", &
              ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
        rc = ESMF_FAILURE
        return
      endif

    endif  ! status

    if (present(rank)) then
      rank = lrank
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_Field_GetFldPtr

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_GetFldPtr(FB, fldname, fldptr1, fldptr2, rank, rc)
    use med_constants_mod , only : R8
    use ESMF              , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_Field

    ! ----------------------------------------------
    ! Get pointer to a field bundle field
    ! ----------------------------------------------
    type(ESMF_FieldBundle) , intent(in)              :: FB
    character(len=*)       , intent(in)              :: fldname
    real(R8), pointer      , intent(inout), optional :: fldptr1(:)
    real(R8), pointer      , intent(inout), optional :: fldptr2(:,:)
    integer                , intent(out),   optional :: rank
    integer                , intent(out),   optional :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer          :: lrank
    integer          :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_FB_GetFldPtr)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    if (.not.present(rc)) then
       call ESMF_LogWrite(trim(subname)//": ERROR rc not present "//trim(fldname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    rc = ESMF_SUCCESS

    if (.not. shr_nuopc_methods_FB_FldChk(FB, trim(fldname), rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": ERROR field "//trim(fldname)//" not in FB ", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), field=lfield, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_Field_GetFldPtr(lfield, &
         fldptr1=fldptr1, fldptr2=fldptr2, rank=lrank, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (present(rank)) then
      rank = lrank
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_GetFldPtr

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_SetFldPtr(FB, fldname, val, rc)
    use med_constants_mod, only : R8
    use ESMF, only : ESMF_FieldBundle, ESMF_Field

    type(ESMF_FieldBundle), intent(in)  :: FB
    character(len=*)      , intent(in)  :: fldname
    real(R8)    , intent(in)  :: val
    integer               , intent(out) :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer          :: lrank
    real(R8), pointer :: fldptr1(:)
    real(R8), pointer :: fldptr2(:,:)
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_FB_SetFldPtr)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call shr_nuopc_methods_FB_GetFldPtr(FB, fldname, fldptr1, fldptr2, lrank, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (lrank == 0) then
      ! no local data
    elseif (lrank == 1) then
      fldptr1 = val
    elseif (lrank == 2) then
      fldptr2 = val
    else
       call ESMF_LogWrite(trim(subname)//": ERROR in rank "//trim(fldname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_SetFldPtr

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_State_GetFldPtr(ST, fldname, fldptr1, fldptr2, rank, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------
    use med_constants_mod, only : R8
    use ESMF, only : ESMF_State, ESMF_Field, ESMF_StateGet

    type(ESMF_State),            intent(in)              :: ST
    character(len=*),            intent(in)              :: fldname
    real(R8), pointer, intent(inout), optional :: fldptr1(:)
    real(R8), pointer, intent(inout), optional :: fldptr2(:,:)
    integer         ,            intent(out),   optional :: rank
    integer         ,            intent(out),   optional :: rc

    ! local variables
    type(ESMF_Field)           :: lfield
    integer                    :: lrank
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_State_GetFldPtr)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    if (.not.present(rc)) then
       call ESMF_LogWrite(trim(subname)//": ERROR rc not present "//trim(fldname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    rc = ESMF_SUCCESS

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_Field_GetFldPtr(lfield, &
         fldptr1=fldptr1, fldptr2=fldptr2, rank=lrank, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (present(rank)) then
      rank = lrank
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_State_GetFldPtr

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_State_SetFldPtr(ST, fldname, val, rc)
    use med_constants_mod, only : R8
    use ESMF, only : ESMF_State, ESMF_Field

    type(ESMF_State)  , intent(in)  :: ST
    character(len=*)  , intent(in)  :: fldname
    real(R8), intent(in)  :: val
    integer           , intent(out) :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer          :: lrank
    real(R8), pointer :: fldptr1(:)
    real(R8), pointer :: fldptr2(:,:)
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_State_SetFldPtr)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call shr_nuopc_methods_State_GetFldPtr(ST, fldname, fldptr1, fldptr2, lrank, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (lrank == 0) then
      ! no local data
    elseif (lrank == 1) then
      fldptr1 = val
    elseif (lrank == 2) then
      fldptr2 = val
    else
       call ESMF_LogWrite(trim(subname)//": ERROR in rank "//trim(fldname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_State_SetFldPtr

  !-----------------------------------------------------------------------------

  logical function shr_nuopc_methods_FieldPtr_Compare1(fldptr1, fldptr2, cstring, rc)
    use med_constants_mod, only : R8
    real(R8), pointer, intent(in)  :: fldptr1(:)
    real(R8), pointer, intent(in)  :: fldptr2(:)
    character(len=*)           , intent(in)  :: cstring
    integer                    , intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_FieldPtr_Compare1)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    shr_nuopc_methods_FieldPtr_Compare1 = .false.
    if (lbound(fldptr2,1) /= lbound(fldptr1,1) .or. &
        ubound(fldptr2,1) /= ubound(fldptr1,1)) then
      call ESMF_LogWrite(trim(subname)//": ERROR in data size "//trim(cstring), ESMF_LOGMSG_ERROR, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      write(msgString,*) trim(subname)//': fldptr1 ',lbound(fldptr1),ubound(fldptr1)
      call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
      write(msgString,*) trim(subname)//': fldptr2 ',lbound(fldptr2),ubound(fldptr2)
      call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
    else
      shr_nuopc_methods_FieldPtr_Compare1 = .true.
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end function shr_nuopc_methods_FieldPtr_Compare1

  !-----------------------------------------------------------------------------

  logical function shr_nuopc_methods_FieldPtr_Compare2(fldptr1, fldptr2, cstring, rc)
    use med_constants_mod, only : R8
    real(R8), pointer, intent(in)  :: fldptr1(:,:)
    real(R8), pointer, intent(in)  :: fldptr2(:,:)
    character(len=*)           , intent(in)  :: cstring
    integer                    , intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_FieldPtr_Compare2)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    shr_nuopc_methods_FieldPtr_Compare2 = .false.
    if (lbound(fldptr2,2) /= lbound(fldptr1,2) .or. &
        lbound(fldptr2,1) /= lbound(fldptr1,1) .or. &
        ubound(fldptr2,2) /= ubound(fldptr1,2) .or. &
        ubound(fldptr2,1) /= ubound(fldptr1,1)) then
      call ESMF_LogWrite(trim(subname)//": ERROR in data size "//trim(cstring), ESMF_LOGMSG_ERROR, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      write(msgString,*) trim(subname)//': fldptr2 ',lbound(fldptr2),ubound(fldptr2)
      call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
      write(msgString,*) trim(subname)//': fldptr1 ',lbound(fldptr1),ubound(fldptr1)
      call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
    else
      shr_nuopc_methods_FieldPtr_Compare2 = .true.
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end function shr_nuopc_methods_FieldPtr_Compare2

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_State_GeomPrint(state, string, rc)
    use ESMF, only : ESMF_State, ESMF_Field, ESMF_StateGet
    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    type(ESMF_Field)  :: lfield
    integer           :: fieldcount
    integer :: dbrc
    character(len=*),parameter  :: subname='(shr_nuopc_methods_State_GeomPrint)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemCount=fieldCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (fieldCount > 0) then
      call shr_nuopc_methods_State_GetFieldN(state, 1, lfield, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_Field_GeomPrint(lfield, string, rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//":"//trim(string)//": no fields", ESMF_LOGMSG_INFO, rc=dbrc)
    endif  ! fieldCount > 0

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_State_GeomPrint

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_GeomPrint(FB, string, rc)
    use ESMF, only : ESMF_FieldBundle, ESMF_Field, ESMF_FieldBundleGet

    type(ESMF_FieldBundle), intent(in)  :: FB
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    type(ESMF_Field)  :: lfield
    integer           :: fieldcount
    integer :: dbrc
    character(len=*),parameter  :: subname='(shr_nuopc_methods_FB_GeomPrint)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (fieldCount > 0) then

      call shr_nuopc_methods_Field_GeomPrint(lfield, string, rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//":"//trim(string)//": no fields", ESMF_LOGMSG_INFO, rc=dbrc)
    endif  ! fieldCount > 0

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_GeomPrint

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_Field_GeomPrint(field, string, rc)
    use med_constants_mod, only : R8
    use ESMF, only : ESMF_Field, ESMF_Grid, ESMF_Mesh
    use ESMF, only : ESMF_FieldGet, ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID, ESMF_FIELDSTATUS_EMPTY

    type(ESMF_Field), intent(in)  :: field
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    type(ESMF_Grid)     :: lgrid
    type(ESMF_Mesh)     :: lmesh
    integer             :: lrank
    real(R8), pointer :: dataPtr1(:)
    real(R8), pointer :: dataPtr2(:,:)
    integer :: dbrc
    character(len=*),parameter  :: subname='(shr_nuopc_methods_Field_GeomPrint)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, status=status, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (status == ESMF_FIELDSTATUS_EMPTY) then
       call ESMF_LogWrite(trim(subname)//":"//trim(string)//": ERROR field does not have a geom yet ", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (geomtype == ESMF_GEOMTYPE_GRID) then
      call ESMF_FieldGet(field, grid=lgrid, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_Grid_Print(lgrid, string, rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    elseif (geomtype == ESMF_GEOMTYPE_MESH) then
      call ESMF_FieldGet(field, mesh=lmesh, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_Mesh_Print(lmesh, string, rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call shr_nuopc_methods_Field_GetFldPtr(field, &
         fldptr1=dataPtr1, fldptr2=dataPtr2, rank=lrank, abort=.false., rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (lrank == 0) then
      ! no local data
    elseif (lrank == 1) then
      write (msgString,*) trim(subname)//":"//trim(string)//": dataptr bounds dim=1 ",lbound(dataptr1,1),ubound(dataptr1,1)
      call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    elseif (lrank == 2) then
      write (msgString,*) trim(subname)//":"//trim(string)//": dataptr bounds dim=1 ",lbound(dataptr2,1),ubound(dataptr2,1)
      call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      write (msgString,*) trim(subname)//":"//trim(string)//": dataptr bounds dim=2 ",lbound(dataptr2,2),ubound(dataptr2,2)
      call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    elseif (lrank == 0) then
      ! means data allocation does not exist yet
      continue
    else
       call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_Field_GeomPrint

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_Mesh_Print(mesh, string, rc)
    use ESMF, only: ESMF_Mesh, ESMF_DistGrid, ESMF_MeshGet, ESMF_DistGridGet
    use ESMF, only: ESMF_DELayoutGet, ESMF_DELayout
    use ESMF, only: ESMF_MeshStatus_Flag, ESMF_MeshStatus_Complete
    type(ESMF_Mesh) , intent(in)  :: mesh
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    type(ESMF_Distgrid)         :: distgrid
    type(ESMF_DELayout)         :: delayout
    integer                     :: pdim, sdim, nnodes, nelements
    integer                     :: localDeCount
    integer                     :: DeCount
    integer                     :: dimCount, tileCount
    integer, allocatable        :: minIndexPTile(:,:), maxIndexPTile(:,:)
    type(ESMF_MeshStatus_Flag)  :: meshStatus
    logical                     :: elemDGPresent, nodeDGPresent
    integer :: dbrc
    character(len=*),parameter  :: subname='(shr_nuopc_methods_Mesh_Print)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_MeshGet(mesh, elementDistGridIsPresent=elemDGPresent, &
         nodalDistgridIsPresent=nodeDGPresent, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MeshGet(mesh, status=meshStatus, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! first get the distgrid, which should be available
    if (elemDGPresent) then
       call ESMF_MeshGet(mesh, elementDistgrid=distgrid, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//": distGrid=element"
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_DistGridGet(distgrid, deLayout=deLayout, dimCount=dimCount, &
            tileCount=tileCount, deCount=deCount, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    dimCount=", dimCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    tileCount=", tileCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    deCount=", deCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_DELayoutGet(deLayout, localDeCount=localDeCount, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    localDeCount=", localDeCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! allocate minIndexPTile and maxIndexPTile accord. to dimCount and tileCount
       allocate(minIndexPTile(dimCount, tileCount), &
            maxIndexPTile(dimCount, tileCount))

       ! get minIndex and maxIndex arrays
       call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
            maxIndexPTile=maxIndexPTile, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    minIndexPTile=", minIndexPTile
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    maxIndexPTile=", maxIndexPTile
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       deallocate(minIndexPTile, maxIndexPTile)

    endif

    if (nodeDGPresent) then
       call ESMF_MeshGet(mesh, nodalDistgrid=distgrid, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//": distGrid=nodal"
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_DistGridGet(distgrid, deLayout=deLayout, dimCount=dimCount, &
            tileCount=tileCount, deCount=deCount, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    dimCount=", dimCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    tileCount=", tileCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    deCount=", deCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_DELayoutGet(deLayout, localDeCount=localDeCount, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    localDeCount=", localDeCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! allocate minIndexPTile and maxIndexPTile accord. to dimCount and tileCount
       allocate(minIndexPTile(dimCount, tileCount), &
            maxIndexPTile(dimCount, tileCount))

       ! get minIndex and maxIndex arrays
       call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
            maxIndexPTile=maxIndexPTile, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    minIndexPTile=", minIndexPTile
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    maxIndexPTile=", maxIndexPTile
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       deallocate(minIndexPTile, maxIndexPTile)

    endif

    if (.not. elemDGPresent .and. .not. nodeDGPresent) then
       call ESMF_LogWrite(trim(subname)//": cannot print distgrid from mesh", &
            ESMF_LOGMSG_WARNING, rc=rc)
       return
    endif

    ! if mesh is complete, also get additional parameters
    if (meshStatus==ESMF_MESHSTATUS_COMPLETE) then
       ! access localDeCount to show this is a real Grid
       call ESMF_MeshGet(mesh, parametricDim=pdim, spatialDim=sdim, &
            numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//": parametricDim=", pdim
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       write (msgString,*) trim(subname)//":"//trim(string)//": spatialDim=", sdim
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       write (msgString,*) trim(subname)//":"//trim(string)//": numOwnedNodes=", nnodes
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       write (msgString,*) trim(subname)//":"//trim(string)//": numOwnedElements=", nelements
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_Mesh_Print

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_Grid_Print(grid, string, rc)
    use med_constants_mod, only : R8
    use ESMF, only : ESMF_Grid, ESMF_DistGrid, ESMF_StaggerLoc
    use ESMF, only : ESMF_GridGet, ESMF_DistGridGet, ESMF_GridGetCoord
    use ESMF, only : ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER
    type(ESMF_Grid) , intent(in)  :: grid
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    type(ESMF_Distgrid) :: distgrid
    integer                     :: localDeCount
    integer                     :: DeCount
    integer                     :: dimCount, tileCount
    integer                     :: staggerlocCount, arbdimCount, rank
    type(ESMF_StaggerLoc)       :: staggerloc
    character(len=32)           :: staggerstr
    integer, allocatable        :: minIndexPTile(:,:), maxIndexPTile(:,:)
    real(R8), pointer :: fldptr1(:)
    real(R8), pointer :: fldptr2(:,:)
    integer                     :: n1,n2,n3
    integer :: dbrc
    character(len=*),parameter  :: subname='(shr_nuopc_methods_Grid_Print)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! access localDeCount to show this is a real Grid
    call ESMF_GridGet(grid, localDeCount=localDeCount, distgrid=distgrid, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    write (msgString,*) trim(subname)//":"//trim(string)//": localDeCount=", localDeCount
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get dimCount and tileCount
    call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, deCount=deCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    write (msgString,*) trim(subname)//":"//trim(string)//": dimCount=", dimCount
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    write (msgString,*) trim(subname)//":"//trim(string)//": tileCount=", tileCount
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    write (msgString,*) trim(subname)//":"//trim(string)//": deCount=", deCount
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! allocate minIndexPTile and maxIndexPTile accord. to dimCount and tileCount
    allocate(minIndexPTile(dimCount, tileCount), &
             maxIndexPTile(dimCount, tileCount))

    ! get minIndex and maxIndex arrays
    call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
       maxIndexPTile=maxIndexPTile, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    write (msgString,*) trim(subname)//":"//trim(string)//": minIndexPTile=", minIndexPTile
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    write (msgString,*) trim(subname)//":"//trim(string)//": maxIndexPTile=", maxIndexPTile
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(minIndexPTile, maxIndexPTile)

    ! get staggerlocCount, arbDimCount
!    call ESMF_GridGet(grid, staggerlocCount=staggerlocCount, rc=rc)
!    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

!    write (msgString,*) trim(subname)//":"//trim(string)//": staggerlocCount=", staggerlocCount
!    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
!    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

!    call ESMF_GridGet(grid, arbDimCount=arbDimCount, rc=rc)
!    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

!    write (msgString,*) trim(subname)//":"//trim(string)//": arbDimCount=", arbDimCount
!    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
!    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get rank
    call ESMF_GridGet(grid, rank=rank, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    write (msgString,*) trim(subname)//":"//trim(string)//": rank=", rank
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n1 = 1,2
      if (n1 == 1) then
        staggerloc = ESMF_STAGGERLOC_CENTER
        staggerstr = 'ESMF_STAGGERLOC_CENTER'
      elseif (n1 == 2) then
        staggerloc = ESMF_STAGGERLOC_CORNER
        staggerstr = 'ESMF_STAGGERLOC_CORNER'
      else
        rc = ESMF_FAILURE
        call ESMF_LogWrite(trim(subname)//":staggerloc failure", ESMF_LOGMSG_INFO, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
      call ESMF_GridGetCoord(grid, staggerloc=staggerloc, isPresent=isPresent, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      write (msgString,*) trim(subname)//":"//trim(staggerstr)//" present=",isPresent
      call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (isPresent) then
        do n3 = 0,localDECount-1
        do n2 = 1,dimCount
          if (rank == 1) then
            call ESMF_GridGetCoord(grid,coordDim=n2,localDE=n3,staggerloc=staggerloc,farrayPtr=fldptr1,rc=rc)
            if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
            write (msgString,'(a,2i4,2f16.8)') trim(subname)//":"//trim(staggerstr)//" coord=",n2,n3,minval(fldptr1),maxval(fldptr1)
          endif
          if (rank == 2) then
            call ESMF_GridGetCoord(grid,coordDim=n2,localDE=n3,staggerloc=staggerloc,farrayPtr=fldptr2,rc=rc)
            if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
            write (msgString,'(a,2i4,2f16.8)') trim(subname)//":"//trim(staggerstr)//" coord=",n2,n3,minval(fldptr2),maxval(fldptr2)
          endif
          call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
          if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
        enddo
        enddo
      endif
    enddo

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_Grid_Print

!-----------------------------------------------------------------------------
  subroutine shr_nuopc_methods_Clock_TimePrint(clock,string,rc)
    use ESMF, only : ESMF_Clock, ESMF_Time, ESMF_TimeInterval
    use ESMF, only : ESMF_ClockGet, ESMF_TimeGet, ESMF_TimeIntervalGet

    type(ESMF_Clock),intent(in) :: clock
    character(len=*),intent(in),optional :: string
    integer, intent(out) :: rc

    type(ESMF_Time)      :: time
    type(ESMF_TimeInterval) :: timeStep
    character(len=64)    :: timestr
    character(len=512)   :: lstring
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_Clock_TimePrint)'

    rc = ESMF_SUCCESS

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    if (present(string)) then
      lstring = trim(subname)//":"//trim(string)
    else
      lstring = trim(subname)
    endif

    call ESMF_ClockGet(clock,currtime=time,rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(time,timestring=timestr,rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(lstring)//": currtime = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_ClockGet(clock,starttime=time,rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(time,timestring=timestr,rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(lstring)//": startime = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_ClockGet(clock,stoptime=time,rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(time,timestring=timestr,rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(lstring)//": stoptime = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_ClockGet(clock,timestep=timestep,rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet(timestep,timestring=timestr,rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(lstring)//": timestep = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_Clock_TimePrint

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_Mesh_Write(mesh, string, rc)
    use med_constants_mod, only : R8
    use ESMF, only : ESMF_Mesh, ESMF_MeshGet, ESMF_Array, ESMF_ArrayWrite, ESMF_DistGrid

    type(ESMF_Mesh) ,intent(in)  :: mesh
    character(len=*),intent(in)  :: string
    integer         ,intent(out) :: rc

    ! local
    integer  :: n,l,i,lsize,ndims
    character(len=64) :: name
    type(ESMF_DISTGRID)         :: distgrid
    type(ESMF_Array)            :: array
    real(R8), pointer :: rawdata(:)
    real(R8), pointer :: coord(:)
    integer :: dbrc
    character(len=*),parameter  :: subname='(shr_nuopc_methods_Mesh_Write)'

    rc = ESMF_SUCCESS
    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

#if (1 == 0)
    !--- elements ---

    call ESMF_MeshGet(mesh, spatialDim=ndims, numownedElements=lsize, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(rawdata(ndims*lsize))
    allocate(coord(lsize))

    call ESMF_MeshGet(mesh, elementDistgrid=distgrid, ownedElemCoords=rawdata, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1,ndims
      name = "unknown"
      if (n == 1) name = "lon_element"
      if (n == 2) name = "lat_element"
    do l = 1,lsize
      i = 2*(l-1) + n
      coord(l) = rawdata(i)
      array = ESMF_ArrayCreate(distgrid, farrayPtr=coord, name=name, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call shr_nuopc_methods_Array_diagnose(array, string=trim(string)//"_"//trim(name), rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_ArrayWrite(array, trim(string)//"_"//trim(name)//".nc", overwrite=.true., rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo
    enddo

    deallocate(rawdata,coord)

    !--- nodes ---

    call ESMF_MeshGet(mesh, spatialDim=ndims, numownedNodes=lsize, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(rawdata(ndims*lsize))
    allocate(coord(lsize))

    call ESMF_MeshGet(mesh, nodalDistgrid=distgrid, ownedNodeCoords=rawdata, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1,ndims
      name = "unknown"
      if (n == 1) name = "lon_nodes"
      if (n == 2) name = "lat_nodes"
    do l = 1,lsize
      i = 2*(l-1) + n
      coord(l) = rawdata(i)
      array = ESMF_ArrayCreate(distgrid, farrayPtr=coord, name=name, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call shr_nuopc_methods_Array_diagnose(array, string=trim(string)//"_"//trim(name), rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_ArrayWrite(array, trim(string)//"_"//trim(name)//".nc", overwrite=.true., rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo
    enddo

    deallocate(rawdata,coord)
#else
      call ESMF_LogWrite(trim(subname)//": turned off right now", ESMF_LOGMSG_INFO, rc=dbrc)
#endif

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_Mesh_Write

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_State_GeomWrite(state, string, rc)
    use ESMF, only : ESMF_State, ESMF_Field, ESMF_StateGet
    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    type(ESMF_Field)  :: lfield
    integer           :: fieldcount
    integer :: dbrc
    character(len=*),parameter  :: subname='(shr_nuopc_methods_State_GeomWrite)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemCount=fieldCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (fieldCount > 0) then
      call shr_nuopc_methods_State_getFieldN(state, 1, lfield, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_Field_GeomWrite(lfield, string, rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//":"//trim(string)//": no fields", ESMF_LOGMSG_INFO, rc=dbrc)
    endif  ! fieldCount > 0

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_State_GeomWrite

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_FB_GeomWrite(FB, string, rc)
    use ESMF, only : ESMF_Field, ESMF_FieldBundle, ESMF_FieldBundleGet

    type(ESMF_FieldBundle), intent(in)  :: FB
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    type(ESMF_Field)  :: lfield
    integer           :: fieldcount
    integer :: dbrc
    character(len=*),parameter  :: subname='(shr_nuopc_methods_FB_GeomWrite)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (fieldCount > 0) then
      call shr_nuopc_methods_FB_getFieldN(FB, 1, lfield, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_Field_GeomWrite(lfield, string, rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//":"//trim(string)//": no fields", ESMF_LOGMSG_INFO, rc=dbrc)
    endif  ! fieldCount > 0

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_FB_GeomWrite

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_Field_GeomWrite(field, string, rc)
    use ESMF, only : ESMF_Field, ESMF_Grid, ESMF_Mesh, ESMF_FIeldGet, ESMF_FIELDSTATUS_EMPTY
    use ESMF, only : ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID

    type(ESMF_Field), intent(in)  :: field
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    type(ESMF_Grid)     :: lgrid
    type(ESMF_Mesh)     :: lmesh
    integer :: dbrc
    character(len=*),parameter  :: subname='(shr_nuopc_methods_Field_GeomWrite)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, status=status, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (status == ESMF_FIELDSTATUS_EMPTY) then
      call ESMF_LogWrite(trim(subname)//":"//trim(string)//": ERROR field does not have a geom yet ", ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (geomtype == ESMF_GEOMTYPE_GRID) then
      call ESMF_FieldGet(field, grid=lgrid, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_Grid_Write(lgrid, string, rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    elseif (geomtype == ESMF_GEOMTYPE_MESH) then
      call ESMF_FieldGet(field, mesh=lmesh, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      call shr_nuopc_methods_Mesh_Write(lmesh, string, rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_Field_GeomWrite

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_Grid_Write(grid, string, rc)
    use ESMF, only : ESMF_Grid, ESMF_Array, ESMF_GridGetCoord, ESMF_ArraySet
    use ESMF, only : ESMF_ArrayWrite, ESMF_GridGetItem, ESMF_GridGetCoord

    use ESMF, only : ESMF_GRIDITEM_AREA, ESMF_GRIDITEM_MASK
    use ESMF, only : ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER

    type(ESMF_Grid) ,intent(in)  :: grid
    character(len=*),intent(in)  :: string
    integer         ,intent(out) :: rc

    ! local
    type(ESMF_Array)            :: array
    character(len=64)           :: name
    integer :: dbrc
    character(len=*),parameter  :: subname='(shr_nuopc_methods_Grid_Write)'

    rc = ESMF_SUCCESS
    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    ! -- centers --

    call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
      name = "lon_center"
      call ESMF_GridGetCoord(grid, coordDim=1, staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_ArraySet(array, name=name, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call shr_nuopc_methods_Array_diagnose(array, string=trim(string)//"_"//trim(name), rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_ArrayWrite(array, trim(string)//"_"//trim(name)//".nc", overwrite=.true., rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      name = "lat_center"
      call ESMF_GridGetCoord(grid, coordDim=2, staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_ArraySet(array, name=name, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call shr_nuopc_methods_Array_diagnose(array, string=trim(string)//"_"//trim(name), rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_ArrayWrite(array, trim(string)//"_"//trim(name)//".nc", overwrite=.true., rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    ! -- corners --

    call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, isPresent=isPresent, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
      name = "lon_corner"
      call ESMF_GridGetCoord(grid, coordDim=1, staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
      if (.not. ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) then
        call ESMF_ArraySet(array, name=name, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

        call shr_nuopc_methods_Array_diagnose(array, string=trim(string)//"_"//trim(name), rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

        call ESMF_ArrayWrite(array, trim(string)//"_"//trim(name)//".nc", overwrite=.true., rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif

      name = "lat_corner"
      call ESMF_GridGetCoord(grid, coordDim=2, staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
      if (.not. ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) then
        call ESMF_ArraySet(array, name=name, rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

        call shr_nuopc_methods_Array_diagnose(array, string=trim(string)//"_"//trim(name), rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

        call ESMF_ArrayWrite(array, trim(string)//"_"//trim(name)//".nc", overwrite=.true., rc=rc)
        if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
    endif

    ! -- mask --

    name = "mask"
    call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
      call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, itemflag=ESMF_GRIDITEM_MASK, array=array, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_ArraySet(array, name=name, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call shr_nuopc_methods_Array_diagnose(array, string=trim(string)//"_"//trim(name), rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_ArrayWrite(array, trim(string)//"_"//trim(name)//".nc", overwrite=.true., rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    ! -- area --

    name = "area"
    call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
      call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, itemflag=ESMF_GRIDITEM_AREA, array=array, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_ArraySet(array, name=name, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call shr_nuopc_methods_Array_diagnose(array, string=trim(string)//trim(name), rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_ArrayWrite(array, trim(string)//"_"//trim(name)//".nc", overwrite=.true., rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_methods_Grid_Write

  !-----------------------------------------------------------------------------

  logical function shr_nuopc_methods_Distgrid_Match(distGrid1, distGrid2, rc)
    use ESMF, only : ESMF_DistGrid, ESMF_DistGridGet
    ! Arguments
    type(ESMF_DistGrid), intent(in)     :: distGrid1
    type(ESMF_DistGrid), intent(in)     :: distGrid2
    integer, intent(out), optional  :: rc

    ! Local Variables
    integer                         :: dimCount1, dimCount2
    integer                         :: tileCount1, tileCount2
    integer, allocatable            :: minIndexPTile1(:,:), minIndexPTile2(:,:)
    integer, allocatable            :: maxIndexPTile1(:,:), maxIndexPTile2(:,:)
    integer, allocatable            :: elementCountPTile1(:), elementCountPTile2(:)
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_Distgrid_Match)'

    if (dbug_flag > 10) then
      call ESMF_LogWrite(subname//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    if(present(rc)) rc = ESMF_SUCCESS
    shr_nuopc_methods_Distgrid_Match = .true.

    call ESMF_DistGridGet(distGrid1, &
      dimCount=dimCount1, tileCount=tileCount1, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_DistGridGet(distGrid2, &
      dimCount=dimCount2, tileCount=tileCount2, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( dimCount1 /= dimCount2) then
      shr_nuopc_methods_Distgrid_Match = .false.
      if (dbug_flag > 1) then
        call ESMF_LogWrite(trim(subname)//": Grid dimCount MISMATCH ", &
          ESMF_LOGMSG_INFO, rc=dbrc)
      endif
    endif

    if ( tileCount1 /= tileCount2) then
      shr_nuopc_methods_Distgrid_Match = .false.
      if (dbug_flag > 1) then
        call ESMF_LogWrite(trim(subname)//": Grid tileCount MISMATCH ", &
          ESMF_LOGMSG_INFO, rc=dbrc)
      endif
    endif

    allocate(elementCountPTile1(tileCount1))
    allocate(elementCountPTile2(tileCount2))
    allocate(minIndexPTile1(dimCount1,tileCount1))
    allocate(minIndexPTile2(dimCount2,tileCount2))
    allocate(maxIndexPTile1(dimCount1,tileCount1))
    allocate(maxIndexPTile2(dimCount2,tileCount2))

    call ESMF_DistGridGet(distGrid1, &
      elementCountPTile=elementCountPTile1, &
      minIndexPTile=minIndexPTile1, &
      maxIndexPTile=maxIndexPTile1, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_DistGridGet(distGrid2, &
      elementCountPTile=elementCountPTile2, &
      minIndexPTile=minIndexPTile2, &
      maxIndexPTile=maxIndexPTile2, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( ANY((elementCountPTile1 - elementCountPTile2) .NE. 0) ) then
      shr_nuopc_methods_Distgrid_Match = .false.
      if (dbug_flag > 1) then
        call ESMF_LogWrite(trim(subname)//": Grid elementCountPTile MISMATCH ", &
          ESMF_LOGMSG_INFO, rc=dbrc)
      endif
    endif

    if ( ANY((minIndexPTile1 - minIndexPTile2) .NE. 0) ) then
      shr_nuopc_methods_Distgrid_Match = .false.
      if (dbug_flag > 1) then
        call ESMF_LogWrite(trim(subname)//": Grid minIndexPTile MISMATCH ", &
          ESMF_LOGMSG_INFO, rc=dbrc)
      endif
    endif

    if ( ANY((maxIndexPTile1 - maxIndexPTile2) .NE. 0) ) then
      shr_nuopc_methods_Distgrid_Match = .false.
      if (dbug_flag > 1) then
        call ESMF_LogWrite(trim(subname)//": Grid maxIndexPTile MISMATCH ", &
          ESMF_LOGMSG_INFO, rc=dbrc)
      endif
    endif

    deallocate(elementCountPTile1)
    deallocate(elementCountPTile2)
    deallocate(minIndexPTile1)
    deallocate(minIndexPTile2)
    deallocate(maxIndexPTile1)
    deallocate(maxIndexPTile2)

    ! TODO: Optionally Check Coordinates


    if (dbug_flag > 10) then
      call ESMF_LogWrite(subname//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end function shr_nuopc_methods_Distgrid_Match

!================================================================================

  subroutine shr_nuopc_methods_State_GetScalar(State, scalar_id, value, flds_scalar_name, flds_scalar_num, rc)
    use med_constants_mod , only : R8
    use ESMF              , only : ESMF_SUCCESS, ESMF_State, ESMF_StateGet, ESMF_Field, ESMF_FieldGet
    use ESMF              , only : ESMF_FAILURE, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_LogWrite
    use ESMF              , only : ESMF_LOGMSG_INFO, ESMF_VM, ESMF_VMBroadCast, ESMF_VMGetCurrent

    ! ----------------------------------------------
    ! Get scalar data from State for a particular name and broadcast it to all other pets
    ! ----------------------------------------------

    type(ESMF_State), intent(in)     :: State
    integer,          intent(in)     :: scalar_id
    real(R8),         intent(out)    :: value
    character(len=*), intent(in)     :: flds_scalar_name
    integer,          intent(in)     :: flds_scalar_num
    integer,          intent(inout)  :: rc

    ! local variables
    integer           :: mytask, ierr, len
    type(ESMF_VM)     :: vm
    type(ESMF_Field)  :: field
    real(R8), pointer :: farrayptr(:,:)
    integer           :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_State_GetScalar)'

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (mytask == 0) then
      call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
        call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
        rc = ESMF_FAILURE
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
      endif
    endif
    call ESMF_VMBroadCast(vm, farrayptr(:,scalar_id), 1, 0, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    value = farrayptr(1,scalar_id)


  end subroutine shr_nuopc_methods_State_GetScalar

!================================================================================

  subroutine shr_nuopc_methods_State_SetScalar(value, scalar_id, State, flds_scalar_name, flds_scalar_num,  rc)
    ! ----------------------------------------------
    ! Set scalar data from State for a particular name
    ! ----------------------------------------------
    use med_constants_mod , only : R8
    use ESMF              , only : ESMF_Field, ESMF_State, ESMF_StateGet, ESMF_FieldGet
    use ESMF              , only : ESMF_VM, ESMF_VMGetCurrent, ESMF_VMGet
    real(R8),         intent(in)     :: value
    integer,          intent(in)     :: scalar_id
    type(ESMF_State), intent(inout)  :: State
    character(len=*), intent(in)     :: flds_scalar_name
    integer,          intent(in)     :: flds_scalar_num
    integer,          intent(inout)  :: rc

    ! local variables
    integer           :: mytask
    type(ESMF_Field)  :: field
    type(ESMF_VM)     :: vm
    real(R8), pointer :: farrayptr(:,:)
    integer           :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_State_SetScalar)'

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (mytask == 0) then
      call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
      if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
        call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
        rc = ESMF_FAILURE
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
      endif
      farrayptr(1,scalar_id) = value
    endif

  end subroutine shr_nuopc_methods_State_SetScalar

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_State_UpdateTimestamp(state, time, rc)
    use NUOPC , only : NUOPC_GetStateMemberLists
    use ESMF  , only : ESMF_State, ESMF_Time, ESMF_Field, ESMF_SUCCESS

    type(ESMF_State) , intent(inout) :: state
    type(ESMF_Time)  , intent(in)    :: time
    integer          , intent(out)   :: rc

    ! local variables
    integer                  :: i
    type(ESMF_Field),pointer :: fieldList(:)
    integer                  :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_State_UpdateTimestamp)'

    rc = ESMF_SUCCESS

    call NUOPC_GetStateMemberLists(state, fieldList=fieldList, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    do i=1, size(fieldList)
      call shr_nuopc_methods_Field_UpdateTimestamp(fieldList(i), time, rc=rc)
      if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

  end subroutine shr_nuopc_methods_State_UpdateTimestamp

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_Field_UpdateTimestamp(field, time, rc)
    use ESMF, only : ESMF_Field, ESMF_Time, ESMF_TimeGet, ESMF_AttributeSet, ESMF_ATTNEST_ON, ESMF_SUCCESS

    type(ESMF_Field) , intent(inout) :: field
    type(ESMF_Time)  , intent(in)    :: time
    integer          , intent(out)   :: rc

    ! local variables
    integer :: yy, mm, dd, h, m, s, ms, us, ns
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_Field_UpdateTimestamp)'

    rc = ESMF_SUCCESS

    call ESMF_TimeGet(time, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, ms=ms, us=us, &
      ns=ns, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AttributeSet(field, &
      name="TimeStamp", valueList=(/yy,mm,dd,h,m,s,ms,us,ns/), &
      convention="NUOPC", purpose="Instance", &
      attnestflag=ESMF_ATTNEST_ON, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine shr_nuopc_methods_Field_UpdateTimestamp

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_methods_Print_FieldExchInfo(flag, values, logunit, fldlist, nflds, istr)
    use shr_nuopc_utils_mod , only : shr_nuopc_string_listGetName
    use med_constants_mod   , only : R8
    use ESMF                , only : ESMF_MAXSTR

    ! !DESCRIPTION:
    ! Print out information about values to stdount
    ! - flag sets the level of information:
    ! - print out names of fields in values 2d array
    ! - also print out local max and min of data in values 2d array
    ! If optional argument istr is present, it will be output before any of the information.


    ! !INPUT/OUTPUT PARAMETERS:
    integer          , intent(in)          :: flag  ! info level flag
    real(R8)         , intent(in)          :: values(:,:) ! arrays sent to/recieved from mediator
    integer          , intent(in)          :: logunit
    character(len=*) , intent(in)          :: fldlist
    integer          , intent(in)          :: nflds
    character(*)     , intent(in),optional :: istr  ! string for print

    !--- local ---
    integer                    :: n           ! generic indicies
    integer                    :: nsize       ! grid point in values array
    real(R8)                   :: minl(nflds) ! local min
    real(R8)                   :: maxl(nflds) ! local max
    character(len=ESMF_MAXSTR) :: name

    !--- formats ---
    character(*),parameter :: subName = '(shr_nuopc_methods_Print_FieldExchInfo) '
    character(*),parameter :: F00 = "('(shr_nuopc_methods_Print_FieldExchInfo) ',8a)"
    character(*),parameter :: F01 = "('(shr_nuopc_methods_Print_FieldExchInfo) ',a,i9)"
    character(*),parameter :: F02 = "('(shr_nuopc_methods_Print_FieldExchInfo) ',240a)"
    character(*),parameter :: F03 = "('(shr_nuopc_methods_Print_FieldExchInfo) ',a,2es11.3,i4,2x,a)"
    !-------------------------------------------------------------------------------

    if (flag >= 1) then
       if (present(istr)) then
          write(logunit,*) trim(istr)
       endif
       nsize = size(values, dim=2)
       write(logunit,F01) "local size =",nsize
       write(logunit,F02) "Fldlist = ",trim(fldlist)
    endif

    if (flag >= 2) then
       do n = 1, nflds
          minl(n) = minval(values(n,:))
          maxl(n) = maxval(values(n,:))
          call shr_nuopc_string_listGetName(fldlist, n, name)
          write(logunit,F03) 'l min/max ',minl(n),maxl(n),n,trim(name)
       enddo
    endif

  end subroutine shr_nuopc_methods_Print_FieldExchInfo

end module shr_nuopc_methods_mod
