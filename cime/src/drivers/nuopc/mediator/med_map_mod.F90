module med_map_mod

  use med_constants_mod , only : CX, CS, CL, R8
  use med_constants_mod , only : ispval_mask => med_constants_ispval_mask
  use med_constants_mod , only : czero => med_constants_czero
  use med_constants_mod , only : dbug_flag => med_constants_dbug_flag

  implicit none
  private

  ! public routines
  public :: med_map_RouteHandles_init
  public :: med_map_Fractions_init
  public :: med_map_MapNorm_init
  public :: med_map_FB_Regrid_Norm

  ! private module variables

  character(*)      , parameter :: u_FILE_u    = __FILE__
  ! should this be a module variable?
  integer                       :: srcTermProcessing_Value = 0
  logical                       :: mastertask

!================================================================================
contains
!================================================================================

  subroutine med_map_RouteHandles_init(gcomp, llogunit, rc)

    !---------------------------------------------
    ! Initialize route handles in the mediator
    ! Assumptions:
    !   -  Route handles are created per target field bundles NOT
    !      per individual fields in the bundle
    !   -  ALL fields in the bundle are on identical grids
    !   -  MULTIPLE route handles are going to be generated for
    !      given field bundle source and destination grids
    !    - Route handles will ONLY be created if coupling is active
    !      between n1 and n2
    ! Algorithm
    !     n1=source component index
    !     n2=destination component index
    !     nf=field index
    !     fldListFr(n)%flds(nf) is queried to determine the mapindex and mapfile
    !     and the appropriate route handle is created
    !
    ! Regridding is done on a per-field basis AND only for those fields that have a
    ! valid mapping index for the destination component
    !     n = source field index field index
    !     destcomp = destination component index
    !     The fldsSrc input argument is queried for the mapping type of the field
    !     for the desination component
    !        mapindex = fldsSrc(n)%mapindex(destcomp)
    !     If the mapindex is 0 (there is no valid mapping) then NO mapping is done
    !        for the field
    !---------------------------------------------

    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_LogFlush
    use ESMF                  , only : ESMF_GridComp, ESMF_VM, ESMF_Field, ESMF_PoleMethod_Flag, ESMF_POLEMETHOD_ALLAVG
    use ESMF                  , only : ESMF_GridCompGet, ESMF_VMGet, ESMF_FieldSMMStore, ESMF_RouteHandleIsCreated
    use ESMF                  , only : ESMF_FieldRedistStore, ESMF_FieldRegridStore, ESMF_REGRIDMETHOD_BILINEAR
    use ESMF                  , only : ESMF_UNMAPPEDACTION_IGNORE, ESMF_REGRIDMETHOD_CONSERVE, ESMF_NORMTYPE_FRACAREA
    use ESMF                  , only : ESMF_NORMTYPE_DSTAREA, ESMF_REGRIDMETHOD_PATCH, ESMF_RouteHandlePrint
    use NUOPC                 , only : NUOPC_Write
    use esmFlds               , only : ncomps, compice, compocn, compname
    use esmFlds               , only : fldListFr, fldListTo
    use shr_nuopc_fldList_mod , only : mapnames
    use shr_nuopc_fldList_mod , only : mapbilnr, mapconsf, mapconsd, mappatch, mapfcopy, mapunset, mapfiler
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFieldN
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use med_internalstate_mod , only : InternalState
    use perf_mod              , only : t_startf, t_stopf
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(in)  :: llogunit
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)              :: is_local
    type(ESMF_VM)                    :: vm
    type(ESMF_Field)                 :: fldsrc
    type(ESMF_Field)                 :: flddst
    integer                          :: localPet
    integer                          :: n,n1,n2,m,nf,nflds,ncomp
    integer                          :: SrcMaskValue
    integer                          :: DstMaskValue
    character(len=128)               :: value
    character(len=128)               :: rhname
    character(len=128)               :: rhname_file
    character(len=CS)       :: mapname
    character(len=CX)       :: mapfile
    character(len=CS)       :: string
    integer                          :: mapindex
    logical                          :: rhprint_flag = .false.
    real(R8)     , pointer :: factorList(:)
    character(CL) , pointer :: fldnames(:)
    type(ESMF_PoleMethod_Flag), parameter :: polemethod=ESMF_POLEMETHOD_ALLAVG
    character(len=*), parameter :: subname=' (module_med_map: RouteHandles_init) '
    integer                       :: dbrc
    !-----------------------------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 1) then
       call ESMF_LogWrite("Starting to initialize RHs", ESMF_LOGMSG_INFO)
       call ESMF_LogFlush()
    endif

    rc = ESMF_SUCCESS

    ! Determine mastertask
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    mastertask = .false.
    if (localPet == 0) mastertask=.true.
    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create the necessary route handles
    if (mastertask) write(llogunit,*) ' '
    do n1 = 1, ncomps
       do n2 = 1, ncomps

          dstMaskValue = ispval_mask
          srcMaskValue = ispval_mask
          if (n1 == compocn .or. n1 == compice) srcMaskValue = 0
          if (n2 == compocn .or. n2 == compice) dstMaskValue = 0

          !--- get single fields from bundles
          !--- 1) ASSUMES all fields in the bundle are on identical grids
          !--- 2) MULTIPLE route handles are going to be generated for
          !---    given field bundle source and destination grids

          if (n1 /= n2) then

             ! Determine route handle names
             rhname = trim(compname(n1))//"2"//trim(compname(n2))

             if (is_local%wrap%med_coupling_active(n1,n2)) then ! If coupling is active between n1 and n2

                call shr_nuopc_methods_FB_getFieldN(is_local%wrap%FBImp(n1,n1), 1, fldsrc, rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                call shr_nuopc_methods_FB_getFieldN(is_local%wrap%FBImp(n1,n2), 1, flddst, rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                ! Loop over fields
                do nf = 1,size(fldListFr(n1)%flds)

                   ! Determine the mapping type for mapping field nf from n1 to n2
                   mapindex = fldListFr(n1)%flds(nf)%mapindex(n2)

                   ! separate check first since Fortran does not have short-circuit evaluation
                   if (mapindex == mapunset) cycle

                   ! Create route handle for target mapindex if route handle is required
                   ! (i.e. mapindex /= mapunset) and route handle has not already been created
                   if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(n1,n2,mapindex), rc=rc)) then

                      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                      mapname  = trim(mapnames(mapindex))
                      mapfile  = trim(fldListFr(n1)%flds(nf)%mapfile(n2))
                      string   = trim(rhname)//'_weights'

                      if (mapindex == mapfiler .and. mapfile /= 'unset') then
                         ! TODO: actually error out if mapfile is unset in this case
                         if (mastertask) then
                            write(llogunit,'(4A)') subname,trim(string),' RH '//trim(mapname)//' via input file ',&
                                 trim(mapfile)
                         end if
                         call ESMF_LogWrite(subname // trim(string) //&
                              ' RH '//trim(mapname)//' via input file '//trim(mapfile), ESMF_LOGMSG_INFO, rc=dbrc)
                         call ESMF_FieldSMMStore(fldsrc, flddst, mapfile, &
                              routehandle=is_local%wrap%RH(n1,n2,mapindex), &
                              ignoreUnmatchedIndices=.true., &
                              srcTermProcessing=srcTermProcessing_Value, rc=rc)
                         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                      else if (mapindex == mapfcopy) then
                         if (mastertask) then
                            write(llogunit,'(3A)') subname,trim(string),' RH redist '
                         end if
                         call ESMF_LogWrite(trim(subname) // trim(string) // ' RH redist ', ESMF_LOGMSG_INFO, rc=dbrc)
                         call ESMF_FieldRedistStore(fldsrc, flddst, &
                              routehandle=is_local%wrap%RH(n1,n2,mapindex), &
                              ignoreUnmatchedIndices = .true., rc=rc)
                         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                      else if (mapfile /= 'unset') then
                         if (mastertask) then
                            write(llogunit,'(4A)') subname,trim(string),' RH '//trim(mapname)//' via input file ',&
                                 trim(mapfile)
                         end if
                         call ESMF_LogWrite(subname // trim(string) //&
                              ' RH '//trim(mapname)//' via input file '//trim(mapfile), ESMF_LOGMSG_INFO, rc=dbrc)
                         call ESMF_FieldSMMStore(fldsrc, flddst, mapfile, &
                              routehandle=is_local%wrap%RH(n1,n2,mapindex), &
                              ignoreUnmatchedIndices=.true., &
                              srcTermProcessing=srcTermProcessing_Value, rc=rc)
                         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                      else
                         if (mastertask) write(llogunit,'(3A)') subname,trim(string),&
                              ' RH regrid for '//trim(mapname)//' computed on the fly'
                         call ESMF_LogWrite(subname // trim(string) //&
                              ' RH regrid for '//trim(mapname)//' computed on the fly', ESMF_LOGMSG_INFO, rc=dbrc)
                         if (mapindex == mapbilnr) then
                            srcTermProcessing_Value = 0
                            call ESMF_FieldRegridStore(fldsrc, flddst, &
                                 routehandle=is_local%wrap%RH(n1,n2,mapindex), &
                                 srcMaskValues=(/srcMaskValue/), &
                                 dstMaskValues=(/dstMaskValue/), &
                                 regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
                                 polemethod=polemethod, &
                                 srcTermProcessing=srcTermProcessing_Value, &
                                 factorList=factorList, &
                                 ignoreDegenerate=.true., &
                                 unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
                         else if (mapindex == mapconsf) then
                            call ESMF_FieldRegridStore(fldsrc, flddst, &
                                 routehandle=is_local%wrap%RH(n1,n2,mapindex), &
                                 srcMaskValues=(/srcMaskValue/), &
                                 dstMaskValues=(/dstMaskValue/), &
                                 regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
                                 normType=ESMF_NORMTYPE_FRACAREA, &
                                 srcTermProcessing=srcTermProcessing_Value, &
                                 factorList=factorList, &
                                 ignoreDegenerate=.true., &
                                 unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
                         else if (mapindex == mapconsd) then
                            call ESMF_FieldRegridStore(fldsrc, flddst, &
                                 routehandle=is_local%wrap%RH(n1,n2,mapindex), &
                                 srcMaskValues=(/srcMaskValue/), &
                                 dstMaskValues=(/dstMaskValue/), &
                                 regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
                                 normType=ESMF_NORMTYPE_DSTAREA, &
                                 srcTermProcessing=srcTermProcessing_Value, &
                                 factorList=factorList, &
                                 ignoreDegenerate=.true., &
                                 unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
                         else if (mapindex == mappatch) then
                            call ESMF_FieldRegridStore(fldsrc, flddst, &
                                 routehandle=is_local%wrap%RH(n1,n2,mapindex), &
                                 srcMaskValues=(/srcMaskValue/), &
                                 dstMaskValues=(/dstMaskValue/), &
                                 regridmethod=ESMF_REGRIDMETHOD_PATCH, &
                                 polemethod=polemethod, &
                                 srcTermProcessing=srcTermProcessing_Value, &
                                 factorList=factorList, &
                                 ignoreDegenerate=.true., &
                                 unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
                         end if
                         if (rhprint_flag) then
                            call NUOPC_Write(factorList, "array_med_"//trim(string)//"_consf.nc", rc)
                            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                         endif
                      end if
                      if (rhprint_flag) then
                         call ESMF_LogWrite(trim(subname)//trim(string)//": printing  RH for "//trim(mapname), &
                              ESMF_LOGMSG_INFO, rc=dbrc)

                         call ESMF_RouteHandlePrint(is_local%wrap%RH(n1,n2,mapindex), rc=rc)
                         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                      endif
                      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                      ! Check that a valid route handle has been created
                      if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(n1,n2,mapindex), rc=rc)) then
                         call ESMF_LogWrite(trim(subname)//trim(string)//": failed   RH "//trim(mapname), &
                              ESMF_LOGMSG_INFO, rc=dbrc)
                      endif
                   end if
                end do ! loop over fields
             end if ! if coupling is active between n1 and n2
          end if ! if n1 not equal to n2
       end do ! loop over n2
    end do ! loop over n1

    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_map_RouteHandles_init

  !================================================================================

  subroutine med_map_Fractions_init(gcomp, n1, n2, FBSrc, FBDst, RouteHandle, rc)

    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_LogFlush
    use ESMF                  , only : ESMF_GridComp, ESMF_FieldBundle, ESMF_RouteHandle, ESMF_Field
    use ESMF                  , only : ESMF_FieldRedistStore, ESMF_FieldSMMStore, ESMF_FieldRegridStore
    use ESMF                  , only : ESMF_UNMAPPEDACTION_IGNORE, ESMF_REGRIDMETHOD_CONSERVE, ESMF_NORMTYPE_FRACAREA
    use NUOPC                 , only : NUOPC_CompAttributeGet
    use esmFlds               , only : ncomps, compice, compocn, compname
    use shr_nuopc_fldList_mod , only : mapnames, mapconsf
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFieldN
    use perf_mod              , only : t_startf, t_stopf
    !---------------------------------------------
    ! Initialize initialize additional route handles
    ! for mapping fractions
    !---------------------------------------------

    type(ESMF_GridComp)                    :: gcomp
    integer                , intent(in)    :: n1
    integer                , intent(in)    :: n2
    type(ESMF_FieldBundle) , intent(in)    :: FBSrc
    type(ESMF_FieldBundle) , intent(in)    :: FBDst
    type(ESMF_RouteHandle) , intent(inout) :: RouteHandle
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)   :: fldsrc
    type(ESMF_Field)   :: flddst
    character(len=128) :: rhname
    character(len=CS)  :: mapname
    character(len=CX)  :: mapfile
    character(len=CS)  :: string
    integer            :: SrcMaskValue
    integer            :: DstMaskValue
    real(R8), pointer  :: factorList(:)
    integer            :: dbrc
    character(len=*), parameter :: subname=' (med_map_fractions_init: ) '
    !---------------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 1) then
       call ESMF_LogWrite("Initializing RHs not yet created and needed for mapping fractions", &
            ESMF_LOGMSG_INFO)
       call ESMF_LogFlush()
    endif

    rc = ESMF_SUCCESS

    call shr_nuopc_methods_FB_getFieldN(FBsrc, 1, fldsrc, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_FB_getFieldN(FBDst, 1, flddst, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    dstMaskValue = ispval_mask
    srcMaskValue = ispval_mask
    if (n1 == compocn .or. n1 == compice) srcMaskValue = 0
    if (n2 == compocn .or. n2 == compice) dstMaskValue = 0

    rhname = trim(compname(n1))//"2"//trim(compname(n2))
    string   = trim(rhname)//'_weights'
    if ( (n1 == compocn .and. n2 == compice) .or. (n1 == compice .and. n2 == compocn)) then
       mapfile = 'idmap'
    else
       call ESMF_LogWrite("Querying for attribute "//trim(rhname)//"_fmapname = ", ESMF_LOGMSG_INFO)
       call NUOPC_CompAttributeGet(gcomp, name=trim(rhname)//"_fmapname", value=mapfile, rc=rc)
       mapname = trim(mapnames(mapconsf))
    end if

    if (mapfile == 'idmap') then
       call ESMF_LogWrite(trim(subname) // trim(string) //&
            ' RH '//trim(mapname)// ' is redist', ESMF_LOGMSG_INFO, rc=dbrc)
       call ESMF_FieldRedistStore(fldsrc, flddst, &
            routehandle=RouteHandle, &
            ignoreUnmatchedIndices = .true., rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (mapfile /= 'unset') then
       call ESMF_LogWrite(subname // trim(string) //&
            ' RH '//trim(mapname)//' via input file '//trim(mapfile), ESMF_LOGMSG_INFO, rc=dbrc)
       call ESMF_FieldSMMStore(fldsrc, flddst, mapfile, &
            routehandle=RouteHandle, &
            ignoreUnmatchedIndices=.true., &
            srcTermProcessing=srcTermProcessing_Value, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_LogWrite(subname // trim(string) //&
            ' RH '//trim(mapname)//' computed on the fly '//trim(mapfile), ESMF_LOGMSG_INFO, rc=dbrc)
       call ESMF_FieldRegridStore(fldsrc, flddst, &
            routehandle=RouteHandle, &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
            normType=ESMF_NORMTYPE_FRACAREA, &
            srcTermProcessing=srcTermProcessing_Value, &
            factorList=factorList, &
            ignoreDegenerate=.true., &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    end if

    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_map_Fractions_init

!================================================================================

  subroutine med_map_MapNorm_init(gcomp, llogunit, rc)

    !---------------------------------------
    ! Initialize unity normalization field bundle
    ! and do the mapping for unity normalization up front
    !---------------------------------------

    use ESMF                  , only: ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_LogFlush
    use ESMF                  , only: ESMF_GridComp, ESMF_FieldBundle, ESMF_RouteHandleIsCreated
    use esmFlds               , only: ncomps, compice, compocn, compname
    use med_internalstate_mod , only: InternalState
    use shr_nuopc_scalars_mod , only: flds_scalar_name, flds_scalar_num
    use shr_nuopc_fldList_mod , only: mapnames, nmappers
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_Init
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_Reset
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_Clean
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_FieldRegrid
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_ChkErr
    use perf_mod              , only: t_startf, t_stopf
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(in)  :: llogunit
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    type(ESMF_FieldBundle)     :: FBTmp
    integer                    :: n1, n2, m
    character(len=CS)          :: normname
    character(len=1)           :: cn1,cn2,cm
    real(R8), pointer          :: dataptr(:)
    integer                    :: dbrc
    character(len=*),parameter :: subname='(module_MED_MAP:MapNorm_init)'
    !-----------------------------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 1) then
       call ESMF_LogWrite("Starting to initialize unity map normalizations", ESMF_LOGMSG_INFO)
       call ESMF_LogFlush()
    endif

    rc = ESMF_SUCCESS

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create the normalization field bundles
    normname = 'one'
    do n1 = 1,ncomps
       do n2 = 1,ncomps
          if (n1 /= n2) then
             if (is_local%wrap%med_coupling_active(n1,n2)) then
                do m = 1,nmappers
                   if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(n1,n2,m), rc=rc)) then
                      if (dbug_flag > 1) then
                         write(cn1,'(i1)') n1; write(cn2,'(i1)') n2; write(cm ,'(i1)') m
                         call ESMF_LogWrite(trim(subname)//":"//'creating FBMapNormOne for '&
                              //compname(n1)//'->'//compname(n2)//'with mapping '//mapnames(m), &
                              ESMF_LOGMSG_INFO, rc=dbrc)
                      endif
                      call shr_nuopc_methods_FB_init(FBout=is_local%wrap%FBNormOne(n1,n2,m), &
                           flds_scalar_name=flds_scalar_name, &
                           FBgeom=is_local%wrap%FBImp(n1,n2), &
                           fieldNameList=(/trim(normname)/), name='FBNormOne', rc=rc)
                      if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

                      call shr_nuopc_methods_FB_reset(is_local%wrap%FBNormOne(n1,n2,m), value=czero, rc=rc)
                      if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

                      call shr_nuopc_methods_FB_init(FBout=FBTmp, &
                           flds_scalar_name=flds_scalar_name, &
                           STgeom=is_local%wrap%NStateImp(n1), &
                           fieldNameList=(/trim(normname)/), name='FBTmp', rc=rc)
                      if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

                      call shr_nuopc_methods_FB_GetFldPtr(FBTmp, trim(normname), fldptr1=dataPtr, rc=rc)
                      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                      dataptr(:) = 1.0_R8

                      call shr_nuopc_methods_FB_FieldRegrid(&
                           FBTmp                           , normname, &
                           is_local%wrap%FBNormOne(n1,n2,m), normname, &
                           is_local%wrap%RH(n1,n2,m), rc)
                      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                      call shr_nuopc_methods_FB_clean(FBTmp, rc=rc)
                      if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return
                   end if
                end do
             end if
          end if
       end do
    end do

    if (dbug_flag > 1) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_map_MapNorm_init

  !================================================================================

  subroutine med_map_FB_Regrid_Norm(fldsSrc, srccomp, destcomp, &
       FBSrc, FBDst, FBFrac, FBNormOne, RouteHandles, string, rc)

    ! ----------------------------------------------
    ! Map field bundles with appropriate fraction weighting
    ! ----------------------------------------------

    use NUOPC                 , only: NUOPC_IsConnected
    use ESMF                  , only: ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only: ESMF_LOGMSG_ERROR, ESMF_FAILURE
    use ESMF                  , only: ESMF_FieldBundle, ESMF_FieldBundleIsCreated, ESMF_FieldBundleGet
    use ESMF                  , only: ESMF_RouteHandle, ESMF_RouteHandleIsCreated, ESMF_Field
    use esmFlds               , only: compname
    use shr_nuopc_scalars_mod , only: flds_scalar_name
    use shr_nuopc_fldList_mod , only: mapnames, mapfcopy
    use shr_nuopc_fldList_mod , only: shr_nuopc_fldList_entry_type
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_Init
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_Reset
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_Clean
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_FieldRegrid
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_FldChk
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_FB_Field_diagnose
    use shr_nuopc_methods_mod , only: shr_nuopc_methods_ChkErr
    use shr_nuopc_utils_mod   , only: shr_nuopc_memcheck
    use perf_mod              , only: t_startf, t_stopf
    ! ----------------------------------------------
    ! Map field bundles with appropriate fraction weighting
    ! ----------------------------------------------

    ! input/output variables
    type(shr_nuopc_fldList_entry_type) , pointer       :: fldsSrc(:)
    integer                            , intent(in)    :: srccomp
    integer                            , intent(in)    :: destcomp
    type(ESMF_FieldBundle)             , intent(inout) :: FBSrc
    type(ESMF_FieldBundle)             , intent(inout) :: FBDst
    type(ESMF_FieldBundle)             , intent(in)    :: FBFrac
    type(ESMF_FieldBundle)             , intent(in)    :: FBNormOne(:)
    type(ESMF_RouteHandle)             , intent(inout) :: RouteHandles(:)
    character(len=*), optional         , intent(in)    :: string
    integer                            , intent(out)   :: rc

    ! local variables
    integer                     :: i, n
    type(ESMF_FieldBundle)      :: FBSrcTmp        ! temporary
    type(ESMF_FieldBundle)      :: FBNormSrc       ! temporary
    type(ESMF_FieldBundle)      :: FBNormDst       ! temporary
    integer                     :: mapindex
    character(len=CS)  :: lstring
    character(len=CS)  :: mapnorm
    character(len=CS)  :: fldname
    character(len=CS)  :: csize1, csize2
    real(R8), pointer :: data_srctmp(:)  ! temporary
    real(R8), pointer :: data_src(:)     ! temporary
    real(R8), pointer :: data_dst(:)     ! temporary
    real(R8), pointer :: data_srcnorm(:) ! temporary
    real(R8), pointer :: data_dstnorm(:) ! temporary
    real(R8), pointer :: data_frac(:)    ! temporary
    real(R8), pointer :: data_norm(:)    ! temporary
    character(len=*), parameter :: subname='(module_MED_Map:med_map_Regrid_Norm)'
    integer :: dbrc

    !-------------------------------------------------------------------------------
    call t_startf('MED:'//subname)
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)
    call shr_nuopc_memcheck(subname, 1, mastertask)

    !---------------------------------------

    if (present(string)) then
      lstring = trim(string)
    else
      lstring = " "
    endif

    rc = ESMF_SUCCESS

    !---------------------------------------
    ! First - reset the field bundle on the destination grid to zero
    !---------------------------------------

    call shr_nuopc_methods_FB_reset(FBDst, value=czero, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Loop over all fields in the source field bundle and map them to
    ! the destination field bundle accordingly
    !---------------------------------------

    call ESMF_LogWrite(trim(subname)//" *** mapping from "//trim(compname(srccomp))//" to "//&
         trim(compname(destcomp))//" ***", ESMF_LOGMSG_INFO, rc=dbrc)

    do n = 1,size(fldsSrc)
       ! Determine if field is a scalar - and if so go to next iternation
       fldname  = fldsSrc(n)%shortname
       if (fldname == flds_scalar_name) CYCLE

       ! Determine if there is a map index and if its zero go to next iteration
       mapindex = fldsSrc(n)%mapindex(destcomp)
       if (mapindex == 0) CYCLE
       mapnorm  = fldsSrc(n)%mapnorm(destcomp)

       ! Error checks
       if (.not. shr_nuopc_methods_FB_FldChk(FBSrc, fldname, rc=rc)) then
          call ESMF_LogWrite(trim(subname)//" field not found in FBSrc: "//trim(fldname), ESMF_LOGMSG_INFO, rc=dbrc)
       else if (.not. shr_nuopc_methods_FB_FldChk(FBDst, fldname, rc=rc)) then
          call ESMF_LogWrite(trim(subname)//" field not found in FBDst: "//trim(fldname), ESMF_LOGMSG_INFO, rc=dbrc)
       else if (.not. ESMF_RouteHandleIsCreated(RouteHandles(mapindex), rc=rc)) then
          call ESMF_LogWrite(trim(subname)//trim(lstring)//&
               ": ERROR RH not available for "//mapnames(mapindex)//": fld="//trim(fldname), &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
          rc = ESMF_FAILURE
          return
       end if

       ! Determine if field is FBSrc or FBDst or connected - and if not go to next iteration
       if (.not. shr_nuopc_methods_FB_FldChk(FBSrc, trim(fldname), rc=rc)) then
          if (dbug_flag > 5) then
             call ESMF_LogWrite(trim(subname)//" field not found in FBSrc: "//trim(fldname), ESMF_LOGMSG_INFO, rc=dbrc)
          end if
          CYCLE
       else if (.not. shr_nuopc_methods_FB_FldChk(FBDst, trim(fldname), rc=rc)) then
          if (dbug_flag > 5) then
             call ESMF_LogWrite(trim(subname)//" field not found in FBDst: "//trim(fldname), ESMF_LOGMSG_INFO, rc=dbrc)
          end if
          CYCLE
       end if

       ! Determine the normalization for the map
       mapnorm  = fldsSrc(n)%mapnorm(destcomp)

       ! Do the mapping
       if (mapindex == mapfcopy) then

          call shr_nuopc_methods_FB_FieldRegrid(FBSrc, fldname, FBDst, fldname, RouteHandles(mapindex), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       else

          if (.not. ESMF_FieldBundleIsCreated(FBSrcTmp)) then
             ! Create a new temporary field bundle, FBSrcTmp if needed
             call shr_nuopc_methods_FB_init(FBSrcTmp, flds_scalar_name, FBgeom=FBSrc, &
                  fieldNameList=(/'data_srctmp'/), name='data_srctmp', rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             ! Get pointer to source field data in FBSrcTmp
             call shr_nuopc_methods_FB_GetFldPtr(FBSrcTmp, 'data_srctmp', data_srctmp, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          ! Get pointer to source field data in FBSrc
          call shr_nuopc_methods_FB_GetFldPtr(FBSrc, fldname, data_src, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          if ( trim(mapnorm) /= 'unset' .and. trim(mapnorm) /= 'one' .and. trim(mapnorm) /= 'none') then

             !-------------------------------------------------
             ! fractional normalization
             !-------------------------------------------------

             ! create a temporary field bundle that will contain normalization on the source grid
             if (.not. ESMF_FieldBundleIsCreated(FBNormSrc)) then
                call shr_nuopc_methods_FB_init(FBout=FBNormSrc, flds_scalar_name=flds_scalar_name, &
                     FBgeom=FBSrc, fieldNameList=(/trim(mapnorm)/), name='normsrc', rc=rc)
                if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return
             endif
             call shr_nuopc_methods_FB_reset(FBNormSrc, value=czero, rc=rc)
             if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

             call shr_nuopc_methods_FB_GetFldPtr(FBNormSrc, trim(mapnorm), data_srcnorm, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             ! create a temporary field bundle that will contain normalization on the destination grid
             if (.not. ESMF_FieldBundleIsCreated(FBNormDst)) then
                call shr_nuopc_methods_FB_init(FBout=FBNormDst, flds_scalar_name=flds_scalar_name, &
                     FBgeom=FBDst, fieldNameList=(/trim(mapnorm)/), name='normdst', rc=rc)
                if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return
             endif
             call shr_nuopc_methods_FB_reset(FBNormDst, value=czero, rc=rc)
             if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return

             ! get a pointer to the array of the normalization on the source grid - this must
             ! be the same size is as fraction on the source grid
             call shr_nuopc_methods_FB_GetFldPtr(FBFrac, trim(mapnorm), data_frac, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             ! error checks
             if (size(data_srcnorm) /= size(data_frac)) then
                call ESMF_LogWrite(trim(subname)//" fldname= "//trim(fldname)//" mapnorm= "//trim(mapnorm), &
                     ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
                write(csize1,'(i8)') size(data_srcnorm)
                write(csize2,'(i8)') size(data_frac)
                call ESMF_LogWrite(trim(subname)//": ERROR data_normsrc size "//trim(csize1)//&
                     " and data_frac size "//trim(csize2)//" are inconsistent", &
                     ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
                rc = ESMF_FAILURE
                return
             else if (size(data_srcnorm) /= size(data_srctmp)) then
                write(csize1,'(i8)') size(data_srcnorm)
                write(csize2,'(i8)') size(data_srctmp)
                call ESMF_LogWrite(trim(subname)//": ERROR data_srcnorm size "//trim(csize1)//&
                     " and data_srctmp size "//trim(csize2)//" are inconsistent", &
                     ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
                rc = ESMF_FAILURE
                return
             end if

             ! now fill in the values for data_srcnorm and data_srctmp - these are the two arrays needed for normalization
             ! Note that FBsrcTmp will now have the data_srctmp value
             do i = 1,size(data_frac)
                data_srcnorm(i) = data_frac(i)
                data_srctmp(i)  = data_src(i) * data_frac(i)  ! Multiply initial field by data_frac
             end do

             ! regrid FBSrcTmp to FBDst
             if (trim(fldname) == trim(flds_scalar_name)) then
                call ESMF_LogWrite(trim(subname)//trim(lstring)//": skip : fld="//trim(fldname), &
                     ESMF_LOGMSG_INFO, rc=dbrc)
             else
                call shr_nuopc_methods_FB_FieldRegrid( FBSrcTmp, 'data_srctmp', FBDst, fldname, RouteHandles(mapindex), rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             end if

             call shr_nuopc_methods_FB_FieldRegrid(FBNormSrc, mapnorm, FBNormDst, mapnorm, RouteHandles(mapindex), rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             ! multiply interpolated field (FBDst) by reciprocal of fraction on destination grid (FBNormDst)
             call shr_nuopc_methods_FB_GetFldPtr(FBNormDst, trim(mapnorm), data_dstnorm, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             call shr_nuopc_methods_FB_GetFldPtr(FBDst, trim(fldname), data_dst, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             do i= 1,size(data_dst)
                if (data_dstnorm(i) == 0.0_R8) then
                   data_dst(i) = 0.0_R8
                else
                   data_dst(i) = data_dst(i)/data_dstnorm(i)
                endif
             end do

          else if (trim(mapnorm) == 'one' .or. trim(mapnorm) == 'none') then

             !-------------------------------------------------
             ! unity or no normalization
             !-------------------------------------------------

             ! map source field to destination grid
             mapindex = fldsSrc(n)%mapindex(destcomp)
             call shr_nuopc_methods_FB_FieldRegrid(FBSrc, fldname, FBDst, fldname, RouteHandles(mapindex), rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             ! obtain unity normalization factor and multiply interpolated field by reciprocal of normalization factor
             if (trim(mapnorm) == 'one') then
                call shr_nuopc_methods_FB_GetFldPtr(FBNormOne(mapindex), 'one', data_norm, rc=rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

                call shr_nuopc_methods_FB_GetFldPtr(FBDst, trim(fldname), data_dst, rc=rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                do i= 1,size(data_dst)
                   if (data_norm(i) == 0.0_R8) then
                      data_dst(i) = 0.0_R8
                   else
                      data_dst(i) = data_dst(i)/data_norm(i)
                   endif
                enddo
             end if ! mapnorm is 'one'

          end if ! mapnorm is 'one' or 'nne'
       end if ! mapindex is not mapfcopy and field exists

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_Field_diagnose(FBDst, fldname, &
               string=trim(subname) //' FBImp('//trim(compname(srccomp))//','//trim(compname(destcomp))//') ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

    end do  ! loop over fields

    ! Clean up temporary field bundles
    if (ESMF_FieldBundleIsCreated(FBSrcTmp)) then
       call shr_nuopc_methods_FB_clean(FBSrcTmp, rc=rc)
       if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return
    end if
    if (ESMF_FieldBundleIsCreated(FBNormSrc)) then
       call shr_nuopc_methods_FB_clean(FBNormSrc, rc=rc)
       if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return
    end if
    if (ESMF_FieldBundleIsCreated(FBNormDst)) then
       call shr_nuopc_methods_FB_clean(FBNormDst, rc=rc)
       if (shr_nuopc_methods_chkerr(rc,__line__,u_file_u)) return
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_map_FB_Regrid_Norm

!================================================================================


end module med_map_mod
