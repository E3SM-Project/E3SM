!================================================================================
module shr_nuopc_grid_mod

  implicit none
  private

  public :: shr_nuopc_grid_ArbInit
  public :: shr_nuopc_grid_DEInit
  public :: shr_nuopc_grid_RegInit
  public :: shr_nuopc_grid_MeshInit
  public :: shr_nuopc_grid_ArrayToState
  public :: shr_nuopc_grid_StateToArray
  public :: shr_nuopc_grid_CreateCoords
  public :: shr_nuopc_grid_CopyCoord
  public :: shr_nuopc_grid_CopyItem

 !integer             :: dbug_flag = 6
  integer             :: dbug_flag = 0
  character(len=1024) :: tmpstr
  character(len=1024) :: msgString
  character(len=*), parameter :: u_FILE_u = &
    __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine shr_nuopc_grid_ArbInit(nx_global, ny_global, mpicom, gindex, EGrid, rc)

    !-----------------------------------------
    ! create a Egrid object for Fields
    !-----------------------------------------
    use ESMF, only : ESMF_Grid, ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU
    use ESMF, only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_GridAddItem
    use ESMF, only : ESMF_GridCreate1PeriDim, ESMF_GridAddCoord, ESMF_LogFoundError
    use ESMF, only : ESMF_GRIDITEM_AREA, ESMF_GRIDITEM_MASK, ESMF_COORDSYS_SPH_DEG
    use ESMF, only : ESMF_STAGGERLOC_CENTER, ESMF_GridCreate1PeriDim
    use ESMF, only : ESMF_TYPEKIND_R8, ESMF_TYPEKIND_I4
    use mpi, only : mpi_comm_rank
    integer         , intent(in)    :: nx_global
    integer         , intent(in)    :: ny_global
    integer         , intent(in)    :: mpicom
    integer         , intent(in)    :: gindex(:)
    type(ESMF_Grid) , intent(inout) :: Egrid
    integer         , intent(inout) :: rc

    !--- local ---
    integer          :: n
    integer          :: iam,ierr
    integer          :: lsize
    integer, pointer :: localArbIndex(:,:)
    integer :: dbrc
    character(len=*),parameter :: subname='(shr_nuopc_grid_ArbInit)'
    !--------------------------------------------------------------

    rc = ESMF_SUCCESS

    call MPI_COMM_RANK(mpicom, iam, ierr)
    call ESMF_LogWrite(subname, ESMF_LOGMSG_INFO, rc=dbrc)

    lsize = size(gindex)
    allocate(localArbIndex(lsize,2))
    do n = 1,lsize
       localArbIndex(n,1) = mod(gindex(n)-1,nx_global) + 1
       localArbIndex(n,2) = (gindex(n)-1)/nx_global + 1
    enddo

    Egrid=ESMF_GridCreate1PeriDim(&
         minIndex = (/1,1/), &
         maxIndex = (/nx_global,ny_global/), &
         arbIndexCount = lsize, &
         arbIndexList = localArbIndex, &
         periodicDim = 1, &
         coordSys = ESMF_COORDSYS_SPH_DEG, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    deallocate(localArbIndex)

    call ESMF_GridAddCoord(Egrid, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    ! call ESMF_GridAddCoord(Egrid, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    call ESMF_GridAddItem(Egrid, itemFlag=ESMF_GRIDITEM_MASK, itemTypeKind=ESMF_TYPEKIND_I4, &
         staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    call ESMF_GridAddItem(Egrid, itemFlag=ESMF_GRIDITEM_AREA, itemTypeKind=ESMF_TYPEKIND_R8, &
         staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

  end subroutine shr_nuopc_grid_ArbInit

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_grid_DEInit(gcomp, nx_global, ny_global, mpicom, gindex, Egrid, rc)

    !-----------------------------------------
    ! create a Egrid object for Fields
    !-----------------------------------------
    use shr_kind_mod, only : R8=>shr_kind_r8
    use ESMF, only : ESMF_GridComp, ESMF_Grid, ESMF_GridCompGet, ESMF_VM, ESMF_VMGet
    use ESMF, only : ESMF_GridGetCoord, ESMF_GRIDITEM_MASK, ESMF_GridGet, ESMF_GridGetItem
    use ESMF, only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_LogWrite
    use ESMF, only : ESMF_VMAllReduce, ESMF_DistGridConnectionSet, ESMF_DistGrid, ESMF_DELayout
    use ESMF, only : ESMF_DistGridGet, ESMF_DistGridPrint, ESMF_GridAddCoord, ESMF_GridAddItem
    use ESMF, only : ESMF_GRIDITEM_AREA, ESMF_TYPEKIND_R8, ESMF_TYPEKIND_I4
    use ESMF, only : ESMF_STAGGERLOC_CENTER, ESMF_DELayoutCreate, ESMF_VMAllReduce
    use ESMF, only : ESMF_REDUCE_SUM, ESMF_LOGMSG_INFO, ESMF_DistGridConnection, ESMF_Grid
    use ESMF, only : ESMF_GridCreate, ESMF_COORDSYS_SPH_DEG
    use ESMF, only : ESMF_DistGridCreate
    use mpi, only : mpi_comm_rank
    use shr_sys_mod, only : shr_sys_abort

    type(ESMF_GridComp)             :: gcomp
    integer         , intent(in)    :: nx_global
    integer         , intent(in)    :: ny_global
    integer         , intent(in)    :: mpicom
    integer         , intent(in)    :: gindex(:)
    type(ESMF_Grid) , intent(inout) :: Egrid
    integer         , intent(inout) :: rc

    !--- local ---
    integer             :: n,n1,n2,ig,jg,cnt
    integer             :: de,decount,dimcount
    integer             :: iam,ierr
    integer             :: lsize,gsize,nblocks_tot,ngseg
    integer             :: lbnd(2),ubnd(2)
    integer             :: global_index
    integer, pointer    :: indexList(:)
    integer, pointer    :: deBlockList(:,:,:)
    integer, pointer    :: petMap(:)
    real(r8),pointer    :: falon(:),falat(:)
    real(r8),pointer    :: famask(:),faarea(:)
    integer, pointer    :: iarray2(:,:)
    real(r8),pointer    :: farray2(:,:)
    type(ESMF_DELayout) :: delayout
    type(ESMF_DistGrid) :: distgrid
    integer ,pointer    :: pes_local(:)
    integer ,pointer    :: pes_global(:)
    type(ESMF_VM)       :: vm
    integer             :: petCount
    type(ESMF_DistGridConnection), allocatable :: connectionList(:)
    integer :: dbrc
    character(len=*),parameter :: subname='(shr_nuopc_grid_DEInit)'
    !--------------------------------------------------------------

    call MPI_COMM_RANK(mpicom, iam, ierr)

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u))  return  ! bail out

    call ESMF_VMGet(vm, petCount=petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return  ! bail out

    ! 1 gridcell per DE

    lsize = size(gindex)
    gsize = nx_global * ny_global

    write(tmpstr,'(a,4i8)') subname//' nx,ny,lsize = ',nx_global,ny_global,lsize,gsize
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    nblocks_tot = gsize

    allocate(deBlockList(2,2,nblocks_tot))
    allocate(petMap(nblocks_tot))

    write(tmpstr,'(a,1i8)') subname//' nblocks = ',nblocks_tot
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    allocate(pes_local(gsize))
    allocate(pes_global(gsize))

    pes_local(:) = 0
    do n = 1,lsize
       pes_local(gindex(n)) = iam
    end do

    call ESMF_VMAllReduce(vm, sendData=pes_local, recvData=pes_global, count=nx_global*ny_global, &
         reduceflag=ESMF_REDUCE_SUM, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u))  return  ! bail out

    ! Note that below this is a global search overall all the points - not just the ones
    ! passed in my gindex(:)
    n = 0
    do global_index = 1,gsize
       ig = mod(global_index-1,nx_global) + 1
       jg = (global_index-1)/nx_global + 1
       deBlockList(1,1,n) = ig
       deBlockList(1,2,n) = ig
       deBlockList(2,1,n) = jg
       deBlockList(2,2,n) = jg
       petMap(global_index) = pes_global(global_index)
       ! write(tmpstr,'(a,2i8)') subname//' IDs  = ',n,petMap(n)
       ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       ! write(tmpstr,'(a,3i8)') subname//' iglo = ',n,deBlockList(1,1,n),deBlockList(1,2,n)
       ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       ! write(tmpstr,'(a,3i8)') subname//' jglo = ',n,deBlockList(2,1,n),deBlockList(2,2,n)
       ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    enddo

    deallocate(pes_local)
    deallocate(pes_global)

    delayout = ESMF_DELayoutCreate(petMap, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    allocate(connectionList(1))
    call ESMF_DistGridConnectionSet(connectionList(1), tileIndexA=1, &
         tileIndexB=1, positionVector=(/nx_global, 0/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=u_FILE_u)) &
         return  ! bail out

    distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/nx_global,ny_global/), &
         !          indexflag = ESMF_INDEX_DELOCAL, &
         deBlockList=deBlockList, &
         delayout=delayout, &
         connectionList=connectionList, &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    deallocate(deBlockList)
    deallocate(petMap)
    deallocate(connectionList)

    call ESMF_DistGridPrint(distgrid=distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=u_FILE_u)) &
         return  ! bail out

    call ESMF_DistGridGet(distgrid=distgrid, localDE=0, elementCount=cnt, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    allocate(indexList(cnt))
    !      write(tmpstr,'(a,i8)') subname//' distgrid cnt= ',cnt
    !      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call ESMF_DistGridGet(distgrid=distgrid, localDE=0, seqIndexList=indexList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    !      write(tmpstr,'(a,4i8)') subname//' distgrid list= ',indexList(1),indexList(cnt),minval(indexList), maxval(indexList)
    !      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    deallocate(IndexList)

    Egrid = ESMF_GridCreate(distgrid=distgrid, &
         coordSys = ESMF_COORDSYS_SPH_DEG, &
         gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), &
         rc = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    call ESMF_GridGet(Egrid, localDEcount=DEcount, dimCount=dimCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    !      write(tmpstr,'(a,2i8)') subname//' localDEcount = ',DEcount,lsize
    !      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    !      write(tmpstr,'(a,2i8)') subname//' dimCount = ',dimCount
    !      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    if (DEcount /= lsize) then
       call shr_sys_abort(subname//' DEcount /= lsize')
    endif

    call ESMF_GridAddCoord(Egrid, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    !      call ESMF_GridAddCoord(Egrid, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc)
    !      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    call ESMF_GridAddItem(Egrid, itemFlag=ESMF_GRIDITEM_MASK, itemTypeKind=ESMF_TYPEKIND_I4, &
         staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    call ESMF_GridAddItem(Egrid, itemFlag=ESMF_GRIDITEM_AREA, itemTypeKind=ESMF_TYPEKIND_R8, &
         staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    allocate(falon(lsize),falat(lsize),famask(lsize),faarea(lsize))

    do n = 1,lsize
       DE = n-1

       !         write(tmpstr,'(a,3i8)') subname//' n,DE,lsize ',n,DE,lsize
       !         call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       !         write(tmpstr,'(a,i8,4g13.6)') subname//' grid values ',DE,falon(n),falat(n),famask(n),faarea(n)
       !         call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

       call ESMF_GridGetCoord(Egrid, coordDim=1, localDE=DE, staggerLoc=ESMF_STAGGERLOC_CENTER, &
            computationalLBound=lbnd, computationalUBound=ubnd, &
            farrayPtr=farray2, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       farray2(1,1) = falon(n)

       !         write(tmpstr,'(a,5i8)') subname//' lbnd ubnd ',DE,lbnd,ubnd
       !         call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

       call ESMF_GridGetCoord(Egrid, coordDim=2, localDE=DE, staggerLoc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=farray2, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       farray2(1,1) = falat(n)

       call ESMF_GridGetItem(Egrid, itemflag=ESMF_GRIDITEM_MASK, localDE=DE, staggerloc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=iarray2, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       iarray2(1,1) = nint(famask(n))

       call ESMF_GridGetItem(Egrid, itemflag=ESMF_GRIDITEM_AREA, localDE=DE, staggerloc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=farray2, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       farray2(1,1) = faarea(n)

    enddo
    deallocate(falon,falat,famask,faarea)

  end subroutine shr_nuopc_grid_DEInit

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_grid_RegInit(nx_global, ny_global, mpicom, EGrid, rc)

    !-----------------------------------------
    ! create a Grid object for Fields
    !-----------------------------------------
    use shr_kind_mod, only : R8=>shr_kind_r8
    use mpi, only : mpi_comm_rank
    use ESMF, only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
    use ESMF, only : ESMF_GridCreateNoPeriDimUfrm, ESMF_COORDSYS_SPH_DEG, ESMF_STAGGERLOC_CENTER
    use ESMF, only : ESMF_Grid, ESMF_SUCCESS
    integer         , intent(in)    :: nx_global
    integer         , intent(in)    :: ny_global
    integer         , intent(in)    :: mpicom
    type(ESMF_Grid) ,intent(inout)  :: Egrid
    integer         , intent(inout) :: rc

    !--- local ---
    integer :: iam,ierr
    integer :: dbrc
    character(len=*),parameter :: subname='(shr_nuopc_grid_RegInit)'
    !--------------------------------------------------------------

    rc = ESMF_SUCCESS

    call MPI_COMM_RANK(mpicom, iam, ierr)
    call ESMF_LogWrite(subname, ESMF_LOGMSG_INFO, rc=dbrc)

    Egrid = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/nx_global, ny_global/), &
         minCornerCoord=(/0._R8, -180._R8/), &
         maxCornerCoord=(/360._R8, 180._R8/), &
         coordSys=ESMF_COORDSYS_SPH_DEG, staggerLocList=(/ESMF_STAGGERLOC_CENTER/), &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=u_FILE_u)) &
         return  ! bail out

  end subroutine shr_nuopc_grid_RegInit

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_grid_MeshInit(gcomp, nx_global, ny_global, mpicom, gindex, lon, lat, Emesh, rc)

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
    use mpi, only : mpi_comm_rank

    type(ESMF_GridComp)               :: gcomp
    integer           , intent(in)    :: nx_global
    integer           , intent(in)    :: ny_global
    integer           , intent(in)    :: mpicom
    integer           , intent(in)    :: gindex(:)
    real(r8), pointer , intent(in)    :: lon(:)
    real(r8), pointer , intent(in)    :: lat(:)
    type(ESMF_Mesh)   , intent(inout) :: Emesh
    integer           , intent(inout) :: rc

    !--- local ---
    integer          :: n,n1,n2,de
    integer          :: iam,ierr
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
    integer :: dbrc
    character(len=*),parameter :: subname='(shr_nuopc_grid_MeshInit)'
    !--------------------------------------------------------------

    rc = ESMF_SUCCESS

    call MPI_COMM_RANK(mpicom, iam, ierr)
    call ESMF_LogWrite(subname, ESMF_LOGMSG_INFO, rc=dbrc)

    lsize = size(gindex)

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u))  return  ! bail out

    call ESMF_VMGet(vm, petCount=petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return  ! bail out

    allocate(latG(nx_global*ny_global))
    allocate(lonG(nx_global*ny_global))

    allocate(recvoffsets(petCount))
    allocate(recvCounts(petCount))

    sendData(1) = lsize
    call ESMF_VMGather(vm, sendData=sendData, recvData=recvCounts, count=1, rootPet=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return  ! bail out

    call ESMF_VMBroadCast(vm, bcstData=recvCounts, count=petCount, rootPet=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return  ! bail out

    recvoffsets(1) = 0
    do n = 2,petCount
       recvoffsets(n) = recvoffsets(n-1) + recvCounts(n-1)
    end do

    call ESMF_VMAllGatherV(vm, lat, lsize, latG, recvCounts, recvOffsets, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u))  return  ! bail out

    call ESMF_VMAllGatherV(vm, lon, lsize, lonG, recvCounts, recvOffsets, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u))  return  ! bail out

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
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u))  return  ! bail out

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

    !----- arguments -----
    real(r8)         , intent(inout) :: array(:,:)
    character(len=*) , intent(in)    :: rList
    type(ESMF_State) , intent(inout) :: state
    character(len=*) , intent(in)    :: grid_option
    integer          , intent(out)   :: rc

    !----- local -----
    integer(IN)       :: nflds, lsize, n, nf, DE
    character(len=CS) :: fldname
    type(ESMF_Field)  :: lfield
    real(R8), pointer :: farray2(:,:)
    real(R8), pointer :: farray1(:)
    integer           :: dbrc
    character(*),parameter :: subName = "(shr_nuopc_grid_ArrayToState)"
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    call shr_nuopc_methods_State_reset(state, value = -9999._R8, rc=rc)

    nflds = shr_string_listGetNum(rList)
    lsize = size(array, dim=2)

    do nf = 1,nflds

      rc = ESMF_SUCCESS
      write(6,*)'DEBUG: rlist = ',trim(rList)
      write(6,*)'DBUG: fldname = ',trim(fldname)
      call shr_string_listGetName(rList, nf, fldname, dbrc)
      call ESMF_StateGet(state, itemName=trim(fldname), field=lfield, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) then

        if (dbug_flag > 2) then
           call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" not found on state", &
                ESMF_LOGMSG_INFO, rc=dbrc)
        end if

      elseif (grid_option == "grid_de") then

        if (dbug_flag > 2) then
           call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" copy", ESMF_LOGMSG_INFO, rc=dbrc)
        end if
        do n = 1,lsize
          DE = n-1
          call ESMF_FieldGet(lfield, localDE=DE, farrayPtr=farray2, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          farray2(1,1) = array(nf,n)
        enddo

      elseif (grid_option == 'mesh' .or. grid_option == 'arb') then

        if (dbug_flag > 2) then
           call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" copy", ESMF_LOGMSG_INFO, rc=dbrc)
        end if
        call ESMF_FieldGet(lfield, farrayPtr=farray1, rc=rc)
        do n = 1,lsize
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
           farray1(n) = array(nf,n)
        enddo
        if (dbug_flag > 2) then
           write(tmpstr,'(a,3g13.6)') trim(subname)//":"//trim(fldname)//"=",&
                minval(farray1),maxval(farray1),sum(farray1)
           call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
        end if

      else

         if (dbug_flag > 2) then
            call ESMF_LogWrite(trim(subname)//": fldname = "//&
                 trim(fldname)//" copy skipped due to grid_option", ESMF_LOGMSG_INFO, rc=dbrc)
         end if
     endif

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

    !----- arguments -----
    type(ESMF_State) , intent(in)    :: state
    real(r8)         , intent(inout) :: array(:,:)
    character(len=*) , intent(in)    :: rList
    character(len=*) , intent(in)    :: grid_option
    integer          , intent(out)   :: rc

    !----- local -----
    integer(IN)       :: nflds, lsize, n, nf, DE
    character(len=CS) :: fldname
    type(ESMF_Field)  :: lfield
    real(R8), pointer :: farray2(:,:)
    real(R8), pointer :: farray1(:)
    integer           :: dbrc
    character(*),parameter :: subName = "(shr_nuopc_grid_StateToArray)"
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    nflds = shr_string_listGetNum(rList)
    lsize = size(array, dim=2)

    do nf = 1,nflds

      rc = ESMF_SUCCESS
      call shr_string_listGetName(rList, nf, fldname, dbrc)
      call ESMF_StateGet(state, itemName=trim(fldname), field=lfield, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU)) then

        if (dbug_flag > 2) then
           call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" not found on state", ESMF_LOGMSG_INFO, rc=dbrc)
        end if

      elseif (grid_option == "grid_de") then

        if (dbug_flag > 2) then
           call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" copy", ESMF_LOGMSG_INFO, rc=dbrc)
        end if
        do n = 1,lsize
          DE = n-1
          call ESMF_FieldGet(lfield, localDE=DE, farrayPtr=farray2, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          array(nf,n) = farray2(1,1)
        enddo

      elseif (grid_option == 'mesh' .or. grid_option == 'arb') then

        if (dbug_flag > 2) then
           call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" copy", ESMF_LOGMSG_INFO, rc=dbrc)
        end if
        call ESMF_FieldGet(lfield, farrayPtr=farray1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
        do n = 1,lsize
          array(nf,n) = farray1(n)
        enddo
        if (dbug_flag > 2) then
           write(tmpstr,'(a,3g13.6)') trim(subname)//":"//trim(fldname)//"=",&
                minval(farray1),maxval(farray1),sum(farray1)
           call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
        end if

      else

        if (dbug_flag > 2) then
           call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//&
                " copy skipped due to grid_option", ESMF_LOGMSG_INFO, rc=dbrc)
        end if

      endif

    enddo

  end subroutine shr_nuopc_grid_StateToArray

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_grid_CreateCoords(gridNew, gridOld, rc)
    use ESMF                  , only : ESMF_Grid, ESMF_DistGrid, ESMF_CoordSys_Flag, ESMF_Index_Flag
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_GridGet, ESMF_GridAddCoord
    use ESMF                  , only : ESMF_GridGetCoord, ESMF_GridAddCoord, ESMF_STAGGERLOC_CENTER
    use ESMF                  , only : ESMF_GridCreate, ESMF_SUCCESS, ESMF_STAGGERLOC_CORNER
    use shr_kind_mod          , only : R8=>shr_kind_r8
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr

    type(ESMF_Grid), intent(inout) :: gridNew
    type(ESMF_Grid), intent(inout) :: gridOld
    integer        , intent(out)   :: rc

    ! local variables
    integer                    :: localDE, localDECount
    type(ESMF_DistGrid)        :: distgrid
    type(ESMF_CoordSys_Flag)   :: coordSys
    type(ESMF_Index_Flag)      :: indexflag
    real(R8),pointer           :: dataPtr1(:,:), dataPtr2(:,:)
    integer                    :: dimCount
    integer, pointer           :: gridEdgeLWidth(:), gridEdgeUWidth(:)
    character(len=*),parameter :: subname='(shr_nuopc_methods_grid_createcoords)'
    integer :: dbrc
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(trim(subname)//": tcxA", ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_GridGet(gridold, dimCount=dimCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(gridEdgeLWidth(dimCount),gridEdgeUWidth(dimCount))
    call ESMF_GridGet(gridold,distgrid=distgrid, coordSys=coordSys, indexflag=indexflag, dimCount=dimCount, &
       gridEdgeLWidth=gridEdgeLWidth, gridEdgeUWidth=gridEdgeUWidth, localDECount=localDECount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//": tcxB", ESMF_LOGMSG_INFO, rc=dbrc)

    write(msgString,*) trim(subname)//' localDECount = ',localDECount
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
    write(msgString,*) trim(subname)//' dimCount = ',dimCount
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
    write(msgString,*) trim(subname)//' size(gELW) = ',size(gridEdgeLWidth)
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
    write(msgString,*) trim(subname)//' gridEdgeLWidth = ',gridEdgeLWidth
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)
    write(msgString,*) trim(subname)//' gridEdgeUWidth = ',gridEdgeUWidth
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_LogWrite(trim(subname)//": tcxC", ESMF_LOGMSG_INFO, rc=dbrc)

    gridnew = ESMF_GridCreate(distgrid=distgrid, coordSys=coordSys, indexflag=indexflag, &
       gridEdgeLWidth=gridEdgeLWidth, gridEdgeUWidth=gridEdgeUWidth, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    deallocate(gridEdgeLWidth, gridEdgeUWidth)

    call ESMF_GridAddCoord(gridnew, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_GridAddCoord(gridnew, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    do localDE = 0,localDeCount-1

      call ESMF_GridGetCoord(gridold, coordDim=1, localDE=localDE,  &
        staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=dataPtr1, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_GridGetCoord(gridnew, coordDim=1, localDE=localDE,  &
        staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=dataPtr2, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      dataPtr2 = dataPtr1

      call ESMF_GridGetCoord(gridold, coordDim=2, localDE=localDE,  &
        staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=dataPtr1, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_GridGetCoord(gridnew, coordDim=2, localDE=localDE,  &
        staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=dataPtr2, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      dataPtr2 = dataPtr1

      call ESMF_GridGetCoord(gridold, coordDim=1, localDE=localDE,  &
        staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=dataPtr1, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_GridGetCoord(gridnew, coordDim=1, localDE=localDE,  &
        staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=dataPtr2, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      dataPtr2 = dataPtr1

      call ESMF_GridGetCoord(gridold, coordDim=2, localDE=localDE,  &
        staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=dataPtr1, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_GridGetCoord(gridnew, coordDim=2, localDE=localDE,  &
        staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=dataPtr2, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      dataPtr2 = dataPtr1

    enddo

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_grid_CreateCoords

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_grid_CopyCoord(gridcomp, gridSrc, gridDst, staggerloc, &
       tolerance, compare, invert, rc)
    use shr_kind_mod, only : I8=>shr_kind_i8, I4=> shr_kind_in
    use ESMF, only : ESMF_GridComp, ESMF_Grid, ESMF_StaggerLoc, ESMF_VM, ESMF_DistGrid, ESMF_Array
    use ESMF, only : ESMF_TypeKind_Flag, ESMF_CoordSys_Flag, ESMF_RouteHandle
    use ESMF, only : ESMF_GridGet, ESMF_LogSetError, ESMF_VMGet, ESMF_GridAddCoord, ESMF_GridGetCoord
    use ESMF, only : ESMF_RC_ARG_BAD, ESMF_ArrayGet, ESMF_DistGridGet
    use ESMF, only : ESMF_LogWrite, ESMF_SUCCESS, ESMF_LOGMSG_INFO
    use ESMF, only : ESMF_ArraySMMStore, ESMF_GridCompGet
    use ESMF, only : ESMF_ArraySMM, ESMF_ArraySMMRelease
    use ESMF, only : ESMF_ArrayRedist, ESMF_ArrayRedistStore, ESMF_ArrayRedistRelease
    use ESMF, only : ESMF_RC_NOT_IMPL, ESMF_LogSetError
    use ESMF, only : operator(/=)
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_Distgrid_Match, shr_nuopc_methods_ChkErr

    ! Arguments
    type(ESMF_GridComp) ,  intent(in)           :: gridcomp
    type(ESMF_Grid),       intent(in)           :: gridSrc
    type(ESMF_Grid),       intent(in)           :: gridDst
    type(ESMF_StaggerLoc), intent(in)           :: staggerloc(:)
    real,                  intent(in), optional :: tolerance
    logical,               intent(in), optional :: compare
    integer,               intent(in), optional :: invert(:)
    integer,               intent(out),optional :: rc

    ! Local Variables
    real                              :: l_tolerance
    logical                           :: l_compare
    integer, allocatable              :: l_invert(:)
    integer                           :: i
    type(ESMF_VM)                     :: vm
    type(ESMF_DistGrid)               :: distGridSrc, distGridDst
    type(ESMF_Array)                  :: coordArraySrc, coordArrayDst
    integer(I4),allocatable :: factorList(:)
    integer, allocatable              :: factorIndexList(:,:)
    type(ESMF_RouteHandle)            :: routehandle
    integer                           :: dimCountSrc, dimCountDst
    integer                           :: deCountDst
    integer, allocatable              :: elementCountPDeDst(:)
    integer(I8)             :: sumElementCountPDeDst
    type(ESMF_TypeKind_Flag)          :: coordTypeKindSrc, coordTypeKindDst
    type(ESMF_CoordSys_Flag)          :: coordSysSrc, coordSysDst
    logical                           :: isPresentSrc, isPresentDst
    integer                           :: dimIndex, staggerlocIndex
    integer                           :: localPet
    character(len=10)                 :: numString
    integer :: dbrc
    character(len=*), parameter       :: subname='(shr_nuopc_methods_Grid_CopyCoord)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(subname//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    l_tolerance = 0.0
    if (present(tolerance)) l_tolerance = tolerance
    l_compare = .FALSE.
    if (present(compare)) l_compare = compare
    if (present(invert)) then
      allocate(l_invert(size(invert)))
      l_invert = invert
    else
      allocate(l_invert(1))
      l_invert = -1
    endif

    call ESMF_GridGet(gridSrc, distGrid=distGridSrc, &
      dimCount=dimCountSrc, coordTypeKind=coordTypeKindSrc, &
      coordSys=coordSysSrc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridGet(gridDst, distGrid=distGridDst, &
      dimCount=dimCountDst, coordTypeKind=coordTypeKindDst, &
      coordSys=coordSysDst, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (.NOT. shr_nuopc_methods_Distgrid_Match(distGrid1=distGridSrc, distGrid2=distGridDst, rc=rc)) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg=subname//": Unable to redistribute coordinates. DistGrids do not match.", &
        line=__LINE__, file=u_FILE_u, rcToReturn=rc)
      return  ! bail out
    endif

    if ( dimCountSrc /= dimCountDst) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg=subname//": DIMCOUNT MISMATCH", &
        line=__LINE__, file=u_FILE_u, rcToReturn=rc)
      return  ! bail out
    endif

    if ( coordTypeKindSrc /= coordTypeKindDst) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg=subname//": COORDTYPEKIND MISMATCH", &
        line=__LINE__, file=u_FILE_u, rcToReturn=rc)
      return  ! bail out
    endif

    if ( coordSysSrc /= coordSysDst) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg=subname//": COORDSYS MISMATCH", &
        line=__LINE__, file=u_FILE_u, rcToReturn=rc)
      return  ! bail out
    endif

    do dimIndex=1, dimCountDst
    do staggerlocIndex=1, size(staggerloc)
      call ESMF_GridGetCoord(gridSrc, staggerloc=staggerloc(staggerlocIndex), &
        isPresent=isPresentSrc, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      if(isPresentSrc) then
        call ESMF_GridGetCoord(gridSrc, coordDim=dimIndex, &
          staggerloc=staggerloc(staggerlocIndex), &
          array=coordArraySrc, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call ESMF_GridGetCoord(gridDst, &
          staggerloc=staggerloc(staggerlocIndex), &
          isPresent=isPresentDst, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if(.NOT.isPresentDst) then
          call ESMF_GridAddCoord(gridDst, &
            staggerLoc=staggerloc(staggerlocIndex), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        else
          if(l_compare .EQV. .TRUE.) then
            ! TODO: Compare existing coordinates
            call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
              msg=subname//": Cannot compare existing coordinates.", &
              line=__LINE__, file=u_FILE_u, rcToReturn=rc)
              return  ! bail out
          end if
        endif
        call ESMF_GridGetCoord(gridDst, coordDim=dimIndex, &
          staggerloc=staggerloc(staggerlocIndex), &
          array=coordArrayDst, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArrayGet(coordArraySrc, distGrid=distGridSrc, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArrayGet(coordArrayDst, distGrid=distGridDst, &
          dimCount=dimCountDst, deCount=deCountDst, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (.NOT. shr_nuopc_methods_Distgrid_Match(distGrid1=distGridSrc, distGrid2=distGridDst, rc=rc)) then
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
          msg=subname//": Unable to redistribute coordinates. DistGrids do not match.", &
          line=__LINE__, file=u_FILE_u, rcToReturn=rc)
          return  ! bail out
        endif

        if ( ANY( l_invert == dimIndex )) then
          call ESMF_GridCompGet(gridcomp, vm=vm, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_VMGet(vm, localPet=localPet, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          if (localPet == 0) then
            call ESMF_DistGridGet(distGridDst, deCount=deCountDst, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            allocate(elementCountPDeDst(deCountDst))
            call ESMF_DistGridGet(distGridDst, &
              elementCountPDe=elementCountPDeDst, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

            sumElementCountPDeDst = SUM(elementCountPDeDst)
            if (dbug_flag >= 2) then
              write (numString, "(I10)") sumElementCountPDeDst
              call ESMF_LogWrite(subname//": sumElementCountPDeDst: "//trim(adjustl(numString)), ESMF_LOGMSG_INFO, rc=dbrc)
            endif

            allocate(factorList(sumElementCountPDeDst))
            allocate(factorIndexList(2,sumElementCountPDeDst))

            factorList(:) = 1
            factorIndexList(1,:) = (/(i, i=1, sumElementCountPDeDst, 1)/)
            factorIndexList(2,:) = (/(i, i=sumElementCountPDeDst, 1, -1)/)

            if (dbug_flag >= 2) then
              write (numString, "(I10)") factorIndexList(1,1)
              write (msgString, "(A)") "Src=>Dst: "//trim(adjustl(numString))//"=>"
              write (numString, "(I10)") factorIndexList(2,1)
              write (msgString, "(A)") trim(msgString)//trim(adjustl(numString))
              write (numString, "(I10)") factorIndexList(1,sumElementCountPDeDst)
              write (msgString, "(A)") trim(msgString)//" "//trim(adjustl(numString))//"=>"
              write (numString, "(I10)") factorIndexList(2,sumElementCountPDEDst)
     	      write (msgString, "(A)") trim(msgString)//trim(adjustl(numString))
              call ESMF_LogWrite(subname//": Invert Mapping: "//msgString, ESMF_LOGMSG_INFO, rc=dbrc)
            endif

            call ESMF_ArraySMMStore(srcArray=coordArraySrc, dstArray=coordArrayDst, &
              routehandle=routehandle, factorList=factorList, &
              factorIndexList=factorIndexList, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
            deallocate(elementCountPDeDst)
            deallocate(factorList)
            deallocate(factorIndexList)
          else
            call ESMF_ArraySMMStore(srcArray=coordArraySrc, dstArray=coordArrayDst, &
              routehandle=routehandle, rc=rc)
            if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif

          call ESMF_ArraySMM(srcArray=coordArraySrc, dstArray=coordArrayDst, &
            routehandle=routehandle, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ArraySMMRelease(routehandle=routehandle, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

        else
          call ESMF_ArrayRedistStore(coordArraySrc, coordArrayDst, routehandle=routehandle, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ArrayRedist(coordArraySrc, coordArrayDst, routehandle=routehandle, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ArrayRedistRelease(routehandle=routehandle, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
      else
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
          msg=subname//": SOURCE GRID MISSING STAGGER LOCATION", &
          line=__LINE__, file=u_FILE_u, rcToReturn=rc)
        return  ! bail out
      endif
    enddo
    enddo

    deallocate(l_invert)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(subname//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_grid_CopyCoord

  !-----------------------------------------------------------------------------

  subroutine shr_nuopc_grid_CopyItem(gridcomp, gridSrc, gridDst, item, &
       tolerance, compare, invert, rc)

    use ESMF, only : ESMF_GridComp, ESMF_Grid, ESMF_GridItem_Flag, ESMF_StaggerLoc
    use ESMF, only : ESMF_GridGetItem, ESMF_GridAddItem
    use ESMF, only : ESMF_DistGrid, ESMF_Array, ESMF_RouteHandle, ESMF_TypeKind_Flag
    use ESMF, only : ESMF_CoordSys_Flag, ESMF_LogWrite, ESMF_LOGMSG_INFO
    use ESMF, only : ESMF_LogSetError, ESMF_GridGet, ESMF_ArrayRedist, ESMF_ArrayRedistStore
    use ESMF, only : ESMF_ArrayRedistRelease, ESMF_GridGetItem, ESMF_STAGGERLOC_CENTER
    use ESMF, only : ESMF_RC_ARG_BAD, ESMF_RC_NOT_IMPL, ESMF_ArrayGet
    use ESMF, only : operator(/=)
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_Distgrid_Match, shr_nuopc_methods_ChkErr

    ! ----------------------------------------------
    type(ESMF_GridComp),      intent(in)           :: gridcomp
    type(ESMF_Grid),          intent(in)           :: gridSrc
    type(ESMF_Grid),          intent(in)           :: gridDst
    type(ESMF_GridItem_Flag), intent(in)           :: item(:)
    real,                     intent(in), optional :: tolerance
    logical,                  intent(in), optional :: compare
    integer,                  intent(in), optional :: invert(:)
    integer,                  intent(out),optional :: rc

    ! Local Variables
    real                        :: l_tolerance
    logical                     :: l_compare
    integer, allocatable        :: l_invert(:)
    type(ESMF_StaggerLoc)       :: l_staggerloc
    type(ESMF_DistGrid)         :: distGridSrc, distGridDst
    type(ESMF_Array)            :: itemArraySrc, itemArrayDst
    type(ESMF_RouteHandle)      :: routehandle
    integer                     :: dimCountSrc, dimCountDst
    type(ESMF_TypeKind_Flag)    :: coordTypeKindSrc, coordTypeKindDst
    type(ESMF_CoordSys_Flag)    :: coordSysSrc, coordSysDst
    logical                     :: isPresentSrc, isPresentDst
    integer                     :: itemIndex
    integer                     :: localPet
    character(len=10)           :: numString
    integer :: dbrc
    character(len=*), parameter :: subname='(shr_nuopc_methods_Grid_CopyItem)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(subname//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    l_tolerance = 0.0
    if (present(tolerance)) l_tolerance = tolerance
    l_compare = .FALSE.
    if (present(compare)) l_compare = compare
    if (present(invert)) then
      allocate(l_invert(size(invert)))
      l_invert = invert
    else
      allocate(l_invert(1))
      l_invert = -1
    endif
    l_staggerloc = ESMF_STAGGERLOC_CENTER

    call ESMF_GridGet(gridSrc, distGrid=distGridSrc, &
      dimCount=dimCountSrc, coordTypeKind=coordTypeKindSrc, &
      coordSys=coordSysSrc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridGet(gridDst, distGrid=distGridDst, &
      dimCount=dimCountDst, coordTypeKind=coordTypeKindDst, &
      coordSys=coordSysDst, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (.NOT. shr_nuopc_methods_Distgrid_Match(distGrid1=distGridSrc, distGrid2=distGridDst, rc=rc)) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg=subname//": Unable to redistribute coordinates. DistGrids do not match.", &
        line=__LINE__, file=u_FILE_u, rcToReturn=rc)
      return  ! bail out
    endif

    if ( dimCountSrc /= dimCountDst) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg=subname//": DIMCOUNT MISMATCH", &
        line=__LINE__, file=u_FILE_u, rcToReturn=rc)
      return  ! bail out
    endif

    if ( coordTypeKindSrc /= coordTypeKindDst) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg=subname//": COORDTYPEKIND MISMATCH", &
        line=__LINE__, file=u_FILE_u, rcToReturn=rc)
      return  ! bail out
    endif

    if ( coordSysSrc /= coordSysDst) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg=subname//": COORDSYS MISMATCH", &
        line=__LINE__, file=u_FILE_u, rcToReturn=rc)
      return  ! bail out
    endif

    do itemIndex=1, size(item)
      call ESMF_GridGetItem(gridSrc, itemflag=item(itemIndex), &
        staggerloc=l_staggerloc, isPresent=isPresentSrc, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      if(isPresentSrc) then
        call ESMF_GridGetItem(gridSrc, itemflag=item(itemIndex), &
          staggerloc=l_staggerloc, array=itemArraySrc, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call ESMF_GridGetItem(gridDst, itemflag=item(itemIndex), &
          staggerloc=l_staggerloc, isPresent=isPresentDst, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if(.NOT.isPresentDst) then
          call ESMF_GridAddItem(gridDst, itemflag=item(itemIndex), &
            staggerLoc=l_staggerloc, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        else
          if(l_compare .EQV. .TRUE.) then
            ! TODO: Compare existing coordinates
            call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
              msg=subname//": Cannot compare existing coordinates.", &
              line=__LINE__, file=u_FILE_u, rcToReturn=rc)
              return  ! bail out
          end if
        endif
        call ESMF_GridGetItem(gridDst, itemflag=item(itemIndex), &
          staggerloc=l_staggerloc, array=itemArrayDst, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArrayGet(itemArraySrc, distGrid=distGridSrc, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        call ESMF_ArrayGet(itemArrayDst, distGrid=distGridDst, rc=rc)
        if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        if (.NOT. shr_nuopc_methods_Distgrid_Match(distGrid1=distGridSrc, distGrid2=distGridDst, rc=rc)) then
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg=subname//": Unable to redistribute coordinates. DistGrids do not match.", &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
            return  ! bail out
        endif

        if ( ANY( l_invert > 0 )) then
          ! TODO: Invert Item
            call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
              msg=subname//": Cannot invert item.", &
              line=__LINE__, file=u_FILE_u, rcToReturn=rc)
              return  ! bail out
        else
          call ESMF_ArrayRedistStore(itemArraySrc, itemArrayDst, routehandle=routehandle, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ArrayRedist(itemArraySrc, itemArrayDst, routehandle=routehandle, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_ArrayRedistRelease(routehandle=routehandle, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
        endif
      else
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
          msg=subname//": SOURCE GRID MISSING ITEM", &
          line=__LINE__, file=u_FILE_u, rcToReturn=rc)
        return  ! bail out
      endif
    enddo

    deallocate(l_invert)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(subname//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine shr_nuopc_grid_CopyItem

end module shr_nuopc_grid_mod
