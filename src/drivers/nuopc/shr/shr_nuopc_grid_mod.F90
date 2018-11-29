!================================================================================
module shr_nuopc_grid_mod
    use shr_nuopc_utils_mod, only : shr_nuopc_utils_ChkErr
  implicit none
  private

  public :: shr_nuopc_grid_MeshInit
  public :: shr_nuopc_grid_ArrayToState
  public :: shr_nuopc_grid_StateToArray

  character(len=*), parameter :: u_FILE_u = &
    __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
  subroutine shr_nuopc_grid_MeshInit(gcomp, nx_global, ny_global, gindex, lon, lat, Emesh, rc)

    !-----------------------------------------
    ! create an Emesh object for Fields
    !-----------------------------------------
    use shr_kind_mod, only : R8=>shr_kind_r8
    use ESMF, only : ESMF_GridComp, ESMF_VM, ESMF_Mesh
    use ESMF, only : ESMF_VMGet, ESMF_GridCompGet, ESMF_VMBroadCast, ESMF_VMAllGatherV
    use ESMF, only : ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_LogWrite
    use ESMF, only : ESMF_VMGather, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
    use ESMF, only : ESMF_MeshCreate, ESMF_COORDSYS_SPH_DEG, ESMF_REDUCE_SUM
    use ESMF, only : ESMF_VMAllReduce, ESMF_MESHELEMTYPE_QUAD

    type(ESMF_GridComp)               :: gcomp
    integer           , intent(in)    :: nx_global
    integer           , intent(in)    :: ny_global
    integer           , intent(in)    :: gindex(:)
    real(r8), pointer , intent(in)    :: lon(:)
    real(r8), pointer , intent(in)    :: lat(:)
    type(ESMF_Mesh)   , intent(inout) :: Emesh
    integer           , intent(inout) :: rc

    !--- local ---
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
    character(len=*),parameter :: subname='(shr_nuopc_grid_MeshInit)'
    !--------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname, ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    lsize = size(gindex)

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, petCount=petCount, localpet=iam, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(latG(nx_global*ny_global))
    allocate(lonG(nx_global*ny_global))

    allocate(recvoffsets(petCount))
    allocate(recvCounts(petCount))

    sendData(1) = lsize
    call ESMF_VMGather(vm, sendData=sendData, recvData=recvCounts, count=1, rootPet=0, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMBroadCast(vm, bcstData=recvCounts, count=petCount, rootPet=0, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    recvoffsets(1) = 0
    do n = 2,petCount
       recvoffsets(n) = recvoffsets(n-1) + recvCounts(n-1)
    end do

    call ESMF_VMAllGatherV(vm, lat, lsize, latG, recvCounts, recvOffsets, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMAllGatherV(vm, lon, lsize, lonG, recvCounts, recvOffsets, rc=rc)
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

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
       ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

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
    if (shr_nuopc_utils_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1,numNodes
       nodeOwners(n) = pes_global(iurpts(n))
    end do
    deallocate(pes_local)
    deallocate(pes_global)

    ! do n = 1,numtotelems
    !   write(tmpstr,'(2a,2i8,2g13.6)') subname,' elemA = ',n,elemIds(n),elemCoords(2*n-1:2*n)
    !   call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    !   write(tmpstr,'(2a,6i8)') subname,' elemB = ',n,elemIds(n),nodeIds(elemConn(4*n-3)),&
    !      nodeIds(elemConn(4*n-2)),nodeIds(elemConn(4*n-1)),nodeIds(elemConn(4*n))
    !   call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    ! enddo
    ! do n = 1,numNodes
    !   write(tmpstr,'(2a,3i8,2g13.6)') subname,' nodesA = ',n,nodeIds(n),nodeOwners(n),nodeCoords(2*n-1:2*n)
    !   call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    ! enddo

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
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    deallocate(iurpts)
    deallocate(nodeIds, nodeCoords, nodeOwners)
    deallocate(elemIds, elemTypes, elemConn, elemCoords)

  end subroutine shr_nuopc_grid_MeshInit

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_grid_ArrayToState(array, rList, state, grid_option, rc)

    ! copy array data to state fields
    use ESMF                  , only : ESMF_State, ESMF_Field, ESMF_SUCCESS
    use ESMF                  , only : ESMF_LogWrite, ESMF_FieldGet, ESMF_StateGet
    use ESMF                  , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_LOGMSG_INFO
    use shr_kind_mod          , only : R8=>shr_kind_r8, CS=>shr_kind_cs, IN=>shr_kind_in
    use shr_string_mod        , only : shr_string_listGetName, shr_string_listGetNum
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_reset
    use med_constants_mod, only : CL
    !----- arguments -----
    real(r8)         , intent(inout) :: array(:,:)
    character(len=*) , intent(in)    :: rList
    type(ESMF_State) , intent(inout) :: state
    character(len=*) , intent(in)    :: grid_option
    integer          , intent(out)   :: rc

    !----- local -----
    integer(IN)       :: nflds, lsize, n, nf
    character(len=CS) :: fldname
    type(ESMF_Field)  :: lfield
    real(R8), pointer :: farray1(:)
    integer           :: dbrc
    character(len=CL) :: tmpstr
    character(*),parameter :: subName = "(shr_nuopc_grid_ArrayToState)"
    !----------------------------------------------------------

    rc = ESMF_SUCCESS
    call shr_nuopc_methods_State_reset(state, value = 0.0_r8, rc=rc)

    lsize = size(array, dim=2)
    nflds = shr_string_listGetNum(rList)
    do nf = 1,nflds
      call shr_string_listGetName(rList, nf, fldname, dbrc)

      call ESMF_StateGet(state, itemName=trim(fldname), field=lfield, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) then
         call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" not found on state", &
              ESMF_LOGMSG_INFO, rc=dbrc)
      else
         call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" copy", ESMF_LOGMSG_INFO, rc=dbrc)
         call ESMF_FieldGet(lfield, farrayPtr=farray1, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
         do n = 1,lsize
            farray1(n) = array(nf,n)
         enddo
         write(tmpstr,'(a,3g13.6)') trim(subname)//":"//trim(fldname)//"=",minval(farray1),maxval(farray1),sum(farray1)
         call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
      end if
    enddo

  end subroutine shr_nuopc_grid_ArrayToState

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_grid_StateToArray(state, array, rList, grid_option, rc)

    ! copy state fields to array data
    use ESMF           , only : ESMF_State, ESMF_Field
    use ESMF           , only : ESMF_StateGet, ESMF_FieldGet, ESMF_LogFoundError, ESMF_LogWrite
    use ESMF           , only : ESMF_LOGERR_PASSTHRU, ESMF_SUCCESS, ESMF_LOGMSG_INFO
    use shr_kind_mod   , only : R8=>shr_kind_r8, CS=>shr_kind_CS, IN=>shr_kind_in
    use shr_string_mod , only : shr_string_listGetName, shr_string_listGetNum
    use med_constants_mod, only : CL

    !----- arguments -----
    type(ESMF_State) , intent(in)    :: state
    real(r8)         , intent(inout) :: array(:,:)
    character(len=*) , intent(in)    :: rList
    character(len=*) , intent(in)    :: grid_option
    integer          , intent(out)   :: rc

    !----- local -----
    integer(IN)       :: nflds, lsize, n, nf
    character(len=CS) :: fldname
    type(ESMF_Field)  :: lfield
    real(R8), pointer :: farray1(:)
    integer           :: dbrc
    character(len=CL) :: tmpstr
    character(*),parameter :: subName = "(shr_nuopc_grid_StateToArray)"
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    nflds = shr_string_listGetNum(rList)
    lsize = size(array, dim=2)

    do nf = 1,nflds
      call shr_string_listGetName(rList, nf, fldname, dbrc)
      call ESMF_StateGet(state, itemName=trim(fldname), field=lfield, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) then
         call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" not found on state", &
              ESMF_LOGMSG_INFO, rc=dbrc)
      else
        call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" copy", ESMF_LOGMSG_INFO, rc=dbrc)
        call ESMF_FieldGet(lfield, farrayPtr=farray1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
        do n = 1,lsize
          array(nf,n) = farray1(n)
        enddo
        write(tmpstr,'(a,3g13.6)') trim(subname)//":"//trim(fldname)//"=",&
             minval(farray1),maxval(farray1),sum(farray1)
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

      endif

    enddo

  end subroutine shr_nuopc_grid_StateToArray

end module shr_nuopc_grid_mod
