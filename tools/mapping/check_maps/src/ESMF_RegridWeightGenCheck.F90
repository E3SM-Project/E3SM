! $Id: ESMF_RegridWeightGenCheck.F90 59441 2014-04-22 22:51:36Z mlevy@ucar.edu $
! $URL: https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/trunk_tags/mapping_140422a/check_maps/src/ESMF_RegridWeightGenCheck.F90 $
!
! Earth System Modeling Framework
! Copyright 2002-2012, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!
!===============================================================================
!                            ESMF_RegridWeightGenCheck.F90
! 
! This is the source file for the RegridWeightGen application in ESMF
!===============================================================================

! Offline Regrid Test Program
program OfflineTester

!------------------------------------------------------------------------------
! SPECIFICATION
!------------------------------------------------------------------------------
! !USES:
  use ESMF

  use netcdf

  implicit none

  integer,parameter :: R8 = selected_real_kind(12) ! 8 byte real
  integer,parameter :: R4 = selected_real_kind( 6) ! 4 byte real
  integer,parameter :: RN = kind(1.0)              ! native real
  integer,parameter :: I8 = selected_int_kind (13) ! 8 byte integer
  integer,parameter :: I4 = selected_int_kind ( 6) ! 4 byte integer
  integer,parameter :: IN = kind(1)                ! native integer
  integer,parameter :: CS = 80                     ! short char
  integer,parameter :: CL = 256                    ! long char
  integer,parameter :: CX = 512                    ! extra-long char

!------------------------------------------------------------------------------
! EXECUTION
!------------------------------------------------------------------------------
      integer :: localPet, nPet
      integer :: failCnt, totCnt, numarg
      integer :: status, rc
#ifdef VERBOSE
      logical, parameter :: Verbose=.true.
      logical, parameter :: Debug=.true.
#else
      logical, parameter :: Verbose=.false.
      logical, parameter :: Debug=.false.
#endif
      integer, parameter :: Num_Tests=5
      logical, parameter :: esmf_smm = .false.
      logical :: success, successful_map

      type(ESMF_VM) :: vm

      character(ESMF_MAXSTR) :: wgtfile, title

      real(ESMF_KIND_R8), pointer :: factorList(:)
      integer, pointer            :: factorIndexList(:,:)

      real(ESMF_KIND_R8), pointer :: src_lat(:), src_lon(:), &
                                     dst_lat(:), dst_lon(:), &
                                     src_area(:), dst_area(:), & 
                                     src_mask(:), dst_mask(:), &
                                     src_frac(:), dst_frac(:)

      !--- for local mapping ---
      integer  :: ns, col, row
      real(r8) :: wgt

      integer :: src_dim, dst_dim, nxs, nys, nxd, nyd
      integer :: n, ni, nj, i, j, src, dst
      integer :: ncid, dimid, dimids(2), dimidd(2), vid
      character(len=128) :: fname,fname1,sname,dname
      real(ESMF_KIND_R8),allocatable :: srcfld(:,:),dstfld(:,:)

      real(ESMF_KIND_R8), parameter :: one = 1.0
      real(ESMF_KIND_R8), parameter :: two = 2.0
      real(ESMF_KIND_R8), parameter :: d2r = 3.141592653589793238/180
      real(ESMF_KIND_R8), parameter :: UNINITVAL = 422397696

      real(ESMF_KIND_R8), allocatable :: FsrcArray(:)
      real(ESMF_KIND_R8), allocatable :: FdstArray(:), FdstArrayX(:)

      real(ESMF_KIND_R8), allocatable :: FsrcMatrix(:,:)
      real(ESMF_KIND_R8), allocatable :: FdstMatrix(:,:), FdstMatrixX(:,:)

      type(ESMF_DistGrid) :: src_distgrid, dst_distgrid
      type(ESMF_ArraySpec):: src_arrayspec, dst_arrayspec
      type(ESMF_Array) :: srcArray, dstArray
      type(ESMF_RouteHandle) :: routehandle
      type(ESMF_Grid) :: srcGrid, dstGrid
      type(ESMF_Field) :: srcField, dstField

      real(ESMF_KIND_R8) :: reltotError, reltwoError, avgError
      real(ESMF_KIND_R8) :: reltotErrorBound, reltwoErrorBound
      real(ESMF_KIND_R8) :: totArea, totAreaBound
      real(ESMF_KIND_R8) :: totErrDif, twoErrDif, twoErrX
      real(ESMF_KIND_R8) :: err, maxneg, maxpos
      real(ESMF_KIND_R8) :: maxerr, minerr, maxerr2, minerr2
      real(ESMF_KIND_R8) :: grid1min, grid1max, grid2min, grid2max
      real(ESMF_KIND_R8) :: srcfrac_min, srcfrac_max, dstfrac_min, dstfrac_max 
      real(ESMF_KIND_R8), pointer :: grid1area(:), grid2area(:)
      real(ESMF_KIND_R8), pointer :: grid1areaXX(:), grid2areaXX(:)

      ! Initialize ESMF
      call ESMF_Initialize (defaultCalKind=ESMF_CALKIND_GREGORIAN, &
        defaultlogfilename="RegridWeightGenCheck.Log", &
        logkindflag=ESMF_LOGKIND_MULTI, rc=status)
      if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      if (Debug) print*, "ESMF_Initialized!"

      ! set log to flush after every message
      call ESMF_LogSet(flush=.true., rc=status)
      if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      if (Debug) print*, "ESMF_LogSet!"

      ! get all vm information
      call ESMF_VMGetGlobal(vm, rc=status)
      if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)

      ! set up local pet info
      call ESMF_VMGet(vm, localPet=localPet, petCount=nPet, rc=status)
      if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)

      ! Usage:  ESMF_FieldRegridOffline weight_file
      call ESMF_UtilGetArgC (numarg)
      if (numarg < 1) then
        if (nPet == 0) then
          print *, 'ERROR: insufficient arguments'
          print *, 'USAGE: ESMF_FieldRegridOfflineUTest weight_file'
          print *, 'The weight_file is the output weight file in SCRIP format'
        endif
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif
      call ESMF_UtilGetArg(1, argvalue=wgtfile)
      if (Debug) print*, "ESMF_UtilGetArg!"
      if (verbose) write(6,*) 'opening ',trim(wgtfile)

      !Set finalrc to success
      rc = ESMF_SUCCESS
      failCnt = 0
      totCnt = 0

      ! read in the grid dimensions
      call NCFileInquire(wgtfile, title, src_dim, nxs, nys, &
                                         dst_dim, nxd, nyd, localrc=status)
      if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      if (Debug) print*, "NCFileInquire!"
      if (verbose) write(6,*) 'NCFileInquire done'
      if (verbose) write(6,*) '  src grid size = ',nxs,nys
      if (verbose) write(6,*) '  dst grid size = ',nxd,nyd

    !  only read the data on PET 0 until we get ArrayRead going...
    if (localPet == 0) then

      allocate(src_lat(src_dim))
      allocate(src_lon(src_dim))
      allocate(src_area(src_dim))
      allocate(src_mask(src_dim))
      allocate(src_frac(src_dim))
      allocate(dst_lat(dst_dim))
      allocate(dst_lon(dst_dim))
      allocate(dst_area(dst_dim))
      allocate(dst_mask(dst_dim))
      allocate(dst_frac(dst_dim))

      call GridReadCoords(wgtfile, src_lat, src_lon, src_area, src_mask, src_frac, &
        dst_lat, dst_lon, dst_area, dst_mask, dst_frac, localrc=status)
      if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      if (Debug) print*, "GridReadCoords!"
      if (verbose) write(6,*) 'GridReadCoords done'

      ! create Fortran arrays
      allocate(FsrcArray(src_dim))
      allocate(FdstArray(dst_dim))
      allocate(FdstArrayX(dst_dim))

      ! create Fortran matrices
      allocate(FsrcMatrix(Num_Tests,src_dim))
      allocate(FdstMatrix(Num_Tests,dst_dim))
      allocate(FdstMatrixX(Num_Tests,dst_dim))

      do j=1,Num_Tests
        ! Initialize data
        ! Avoid memory issues by only populating values in rows up to Num_Tests
        if (j.eq.1) then
          ! Test 1: Initial Test (k=2)
          FsrcMatrix(1,:) = two + cos(src_lat)**2*cos(two*src_lon)
          FdstMatrixX(1,:) = two + cos(dst_lat)**2*cos(two*dst_lon)
          if (Verbose) &
            print *, "Test ", j, ": Initial test -- 2 + cos^2(lat) * cos(2 * lon)"
        else if (j.eq.2) then
          ! Test 2: Constant
          FsrcMatrix(2,:) = two
          FdstMatrixX(2,:) = two
          if (Verbose) print *, "Test ", j, ": Constant test -- 2"
        else
          ! Test 3+: Higher Wave Number (k=2^j)
          FsrcMatrix(j,:) = two + cos(src_lat)**2*cos((two**j)*src_lon)
          FdstMatrixX(j,:) = two + cos(dst_lat)**2*cos((two**j)*dst_lon)
          if (Verbose) &
            print *, "Test ", j, &
            ": Higher Wave Number  test -- 2 + cos^2(lat) * cos(",2**j,"* lon)"
        end if
      end do

      FdstMatrix = UNINITVAL
      if (verbose) write(6,*) 'Test Fields Setup done'

      ! deallocate arrays
      deallocate(src_lat)
      deallocate(src_lon)
      deallocate(dst_lat)
      deallocate(dst_lon)
    endif

    if (esmf_smm) then
      if (verbose) write(6,*) 'using esmf_smm'
      do j=1,Num_Tests
        if (verbose) write(6,'(a,i4,a,i4)') 'mapping field ',j,' of ',Num_Tests
        FsrcArray = FsrcMatrix(j,:)
        FdstArray = FdstMatrix(j,:)

        ! create DistGrids for the ESMF Arrays
        src_distgrid = ESMF_DistGridCreate(minIndex=(/1/), &
                                           maxIndex=(/src_dim/), rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        dst_distgrid = ESMF_DistGridCreate(minIndex=(/1/), &
                                           maxIndex=(/dst_dim/), rc=status)      
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        if (Debug) print*, j, ": ESMF_DistGridCreate!"

        ! create dummy grids for fields
        srcGrid = ESMF_GridCreate(src_distgrid, rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT) 
        dstGrid = ESMF_GridCreate(dst_distgrid, rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT) 
        if (Debug) print*, j, ": ESMF_GridCreate!"

        ! create ArraySpecs for the ESMF Arrays
        call ESMF_ArraySpecSet(src_arrayspec, typekind=ESMF_TYPEKIND_R8, rank=1, rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        call ESMF_ArraySpecSet(dst_arrayspec, typekind=ESMF_TYPEKIND_R8, rank=1, rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        if (Debug) print*, j, ": ESMF_ArraySpecSet!"

        ! create the ESMF Arrays
        srcArray = ESMF_ArrayCreate(arrayspec=src_arrayspec, distgrid=src_distgrid, rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        dstArray = ESMF_ArrayCreate(farray=FdstArray, distgrid=dst_distgrid, &
                                    indexflag=ESMF_INDEX_DELOCAL, rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        if (Debug) print*, j, ": ESMF_ArrayCreate!"

        ! Scatter the ESMF Arrays
        call ESMF_ArrayScatter(srcArray, farray=FsrcArray, rootPet=0, rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        if (Debug) print*, j, ": ESMF_ArrayScatter!"

        ! create fields on the empty grid and arrays
        srcField = ESMF_FieldCreate(srcGrid, srcArray, rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)

        ! create fields on the empty grid and arrays
        dstField = ESMF_FieldCreate(dstGrid, dstArray, rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        if (Debug) print*, j, ": ESMF_FieldCreate!"

        if (localPet == 0) then
          ! read in the regriding weights from specified file -> factorList and factorIndex list
          call ESMF_FieldRegridReadSCRIPFileP(wgtfile, factorList, factorIndexList, rc=status)
          if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) &
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          if (Debug) print*, j, ": ESMF_FieldRegridReadSCRIPFileP!"
        endif

        ! Field and Grid way of doing things
        ! store the factorList and factorIndex list into a routehandle for SMM
        call ESMF_FieldSMMStore(srcField=srcField, dstField=dstField, routehandle=routehandle, &
                                factorList=factorList, factorIndexList=factorIndexList, rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        if (Debug) print*, j, ": ESMF_FieldSMMStore!"

        ! compute a Regrid from srcField to dstField
        call ESMF_FieldRegrid(srcField, dstField, routehandle, &
                              zeroregion=ESMF_REGION_SELECT, rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        if (Debug) print*, j, ": ESMF_FieldRegrid!"

        ! ArrayGather the dst array
        call ESMF_ArrayGather(dstArray, farray=FdstArray, rootPet=0, rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        if (Debug) print*, j, ": ESMF_ArrayGather!"
        FdstMatrix(j,:) = FdstArray
        if (verbose) write(6,*) 'mapping field done ',j
      end do
    else    ! esmf_smm
      if (verbose) write(6,*) 'using local smm'

!      status = nf90_open(trim(wgtfile),NF90_NOWRITE,ncid)
!      if(status /= nf90_NoErr) write(6,*) 'nf90_open:',trim(nf90_strerror(status))
!      status = nf90_inq_dimid(ncid,'n_s',dimid)
!      if(status /= nf90_NoErr) write(6,*) 'nf90_inq_dimid n_s:',trim(nf90_strerror(status))
!      status = nf90_inquire_dimension(ncid,dimid,len=ns)
!      if(status /= nf90_NoErr) write(6,*) 'nf90_inq_dim_len n_s:',trim(nf90_strerror(status))
!      allocate(row(ns),col(ns),wgt(ns))
!      status = nf90_inq_varid(ncid,'row',vid)
!      if(status /= nf90_NoErr) write(6,*) 'nf90_varid row:',trim(nf90_strerror(status))
!      status = nf90_get_var(ncid,vid,row)
!      if(status /= nf90_NoErr) write(6,*) 'nf90_get_var row:',trim(nf90_strerror(status))
!      status = nf90_inq_varid(ncid,'col',vid)
!      if(status /= nf90_NoErr) write(6,*) 'nf90_varid col:',trim(nf90_strerror(status))
!      status = nf90_get_var(ncid,vid,col)
!      if(status /= nf90_NoErr) write(6,*) 'nf90_get_var col:',trim(nf90_strerror(status))
!      status = nf90_inq_varid(ncid,'S',vid)
!      if(status /= nf90_NoErr) write(6,*) 'nf90_varid wgt:',trim(nf90_strerror(status))
!      status = nf90_get_var(ncid,vid,wgt)
!      if(status /= nf90_NoErr) write(6,*) 'nf90_get_var wgt:',trim(nf90_strerror(status))
!      status = nf90_close(ncid)
!      if(status /= nf90_NoErr) write(6,*) 'nf90_close:',trim(nf90_strerror(status))

      call ESMF_FieldRegridReadSCRIPFileP(wgtfile, factorList, factorIndexList, rc=status)
      if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      if (Debug) print*, j, ": ESMF_FieldRegridReadSCRIPFileP!"

      ns = size(factorList)
      if (verbose) write(6,*) 'number of weights = ',ns

      ! set zeros in FdstMatrix dst indices as starting point
      do n = 1,ns
         row = factorIndexList(2,n)
         FdstMatrix(:,row) = 0.0_r8
      enddo

      do n = 1,ns
      do j = 1,Num_Tests
         col = factorIndexList(1,n)
         row = factorIndexList(2,n)
         wgt = factorList(n)
         FdstMatrix(j,row) = FdstMatrix(j,row) + FsrcMatrix(j,col)*wgt
      enddo
      enddo
    endif

      ! -----------------------------------------------------------------------
      ! netcdf output
      ! -----------------------------------------------------------------------

    if (localPet == 0) then
       allocate(srcfld(nxs,nys),dstfld(nxd,nyd))
       fname1 = trim(wgtfile)
       n = index(fname1,'/',back=.true.)
       fname = './test_'//fname1(n+1:len_trim(fname1))
       if (verbose) write(6,*) 'writing ',trim(fname)
       status = nf90_create(trim(fname),NF90_CLOBBER,ncid)
       if(status /= nf90_NoErr) write(6,*) 'nf90_create:',trim(nf90_strerror(status))
       status = nf90_def_dim(ncid,'nxs',nxs,dimids(1))
       if(status /= nf90_NoErr) write(6,*) 'nf90_def_dim nxs:',trim(nf90_strerror(status))
       status = nf90_def_dim(ncid,'nys',nys,dimids(2))
       if(status /= nf90_NoErr) write(6,*) 'nf90_def_dim nys:',trim(nf90_strerror(status))
       status = nf90_def_dim(ncid,'nxd',nxd,dimidd(1))
       if(status /= nf90_NoErr) write(6,*) 'nf90_def_dim nxd:',trim(nf90_strerror(status))
       status = nf90_def_dim(ncid,'nyd',nyd,dimidd(2))
       if(status /= nf90_NoErr) write(6,*) 'nf90_def_dim nyd:',trim(nf90_strerror(status))
       do j = 1,Num_Tests
          write(sname,'(a,i2.2)') 'src',j
          write(dname,'(a,i2.2)') 'dst',j
          status = nf90_def_var(ncid,trim(sname),nf90_double,dimids,vid)
          if(status /= nf90_NoErr) write(6,*) 'nf90_def_var sname:',trim(nf90_strerror(status))
          status = nf90_put_att(ncid,vid,'_FillValue',UNINITVAL)
          if(status /= nf90_NoErr) write(6,*) 'nf90_put_att fill sname:',trim(nf90_strerror(status))
          status = nf90_def_var(ncid,trim(dname),nf90_double,dimidd,vid)
          if(status /= nf90_NoErr) write(6,*) 'nf90_def_var dname:',trim(nf90_strerror(status))
          status = nf90_put_att(ncid,vid,'_FillValue',UNINITVAL)
          if(status /= nf90_NoErr) write(6,*) 'nf90_put_att dname:',trim(nf90_strerror(status))
       enddo
       status = nf90_enddef(ncid)
       if(status /= nf90_NoErr) write(6,*) 'nf90_enddef:',trim(nf90_strerror(status))
       do j = 1,Num_Tests
          write(sname,'(a,i2.2)') 'src',j
          write(dname,'(a,i2.2)') 'dst',j
          n = 0
          do nj = 1,nys
          do ni = 1,nxs
             n = n + 1
             srcfld(ni,nj) = FsrcMatrix(j,n)
          enddo
          enddo
          n = 0
          do nj = 1,nyd
          do ni = 1,nxd
             n = n + 1
             dstfld(ni,nj) = FdstMatrix(j,n)
          enddo
          enddo
          status = nf90_inq_varid(ncid,trim(sname),vid)
          if(status /= nf90_NoErr) write(6,*) 'nf90_inq_varid sname:',trim(nf90_strerror(status))
          status = nf90_put_var(ncid,vid,srcfld)
          if(status /= nf90_NoErr) write(6,*) 'nf90_put_var srcfld:',trim(nf90_strerror(status))
          status = nf90_inq_varid(ncid,trim(dname),vid)
          if(status /= nf90_NoErr) write(6,*) 'nf90_inq_varid dname:',trim(nf90_strerror(status))
          status = nf90_put_var(ncid,vid,dstfld)
          if(status /= nf90_NoErr) write(6,*) 'nf90_put_var dstfld:',trim(nf90_strerror(status))
       enddo
       status = nf90_close(ncid)
       if(status /= nf90_NoErr) write(6,*) 'nf90_close:',trim(nf90_strerror(status))
       deallocate(srcfld,dstfld)
    endif     

      ! -----------------------------------------------------------------------
      ! ERROR ANALYSIS - serial
      ! -----------------------------------------------------------------------
    if (localPet == 0) then
      success = .true.
      do j=1,Num_Tests
        if (verbose) write(6,'(a,i4,a,i4)') 'analyzing field ',j,' of ',Num_Tests
        FsrcArray = FsrcMatrix(j,:)
        FdstArray = FdstMatrix(j,:)
        FdstArrayX = FdstMatrixX(j,:)
        totErrDif = 0
        twoErrDif = 0
        twoErrX = 0
        maxerr = 0
        minerr = 1
        maxerr2 = 0
        minerr2 = 1

        ! source error
        grid1min = UNINITVAL
        grid1max = 0
        srcfrac_min = UNINITVAL
        srcfrac_max = 0
        do i=1,src_dim
          if (src_mask(i) /=0) then
            if(FsrcArray(i) < grid1min) grid1min = FsrcArray(i) ! test 2 is original
            if(FsrcArray(i) > grid1max) grid1max = FsrcArray(i) ! test 2 is original
            ! Find min / max frac_a
            if(src_frac(i) < srcfrac_min) srcfrac_min = src_frac(i)
            if(src_frac(i) > srcfrac_max) srcfrac_max = src_frac(i)
          endif
        enddo

        ! destination error
        grid2min = UNINITVAL
        grid2max = 0
        dstfrac_min = UNINITVAL
        dstfrac_max = 0
        successful_map = .false.
        do i=1,dst_dim
          ! Find min / max frac_b
          if(dst_frac(i) < dstfrac_min) dstfrac_min = dst_frac(i)
          if(dst_frac(i) > dstfrac_max) dstfrac_max = dst_frac(i)

          ! don't look in masked cells
          ! if frac is below .999, then a significant portion of this cell is 
          ! missing from the weight calculation and error is misleading here
          ! also don't look in unitialized cells, for the regional to global cases
          if (dst_mask(i) /= 0 .and. dst_frac(i) > .999 .and. FdstArray(i) /= UNINITVAL) then
            successful_map = .true.
            err = FdstArray(i) - FdstArrayX(i)
            totErrDif = totErrDif + abs(err)
            twoErrDif = twoErrDif + err**2
            twoErrX = twoErrX + FdstArrayX(i)**2
            if(err < minerr) minerr = err
            if(err > maxerr) maxerr = err
            if(abs(err) < minerr2) minerr2 = abs(err)
            if(abs(err) > maxerr2) maxerr2 = abs(err)

            ! masking will screw this one up
            if (FdstArray(i) < grid2min) grid2min = FdstArray(i) ! test 2 is original
            if (FdstArray(i) > grid2max) grid2max = FdstArray(i) ! test 2 is original
          endif
        enddo

        if (successful_map) then
          ! maximum negative weight
          maxneg = 0
          maxneg = minval(factorList)
          if (maxneg > 0) maxneg = 0

          ! maximum positive weight
          maxpos = 0
          maxpos = maxval(factorList)

          ! relative error
          reltotError = totErrDif/sum(abs(FdstArrayX))
          reltwoError = sqrt(twoErrDif)/sqrt(twoErrX)
          avgError = (minerr + maxerr)/2

          ! area calculations
          ! use one of src_ or dst_frac, but NOT both!
          allocate(grid1area(src_dim))
          allocate(grid2area(dst_dim))
          grid2area=0.0
          allocate(grid1areaXX(src_dim))
          allocate(grid2areaXX(dst_dim))
          grid1area = FsrcArray*src_area*src_frac

          ! Only calculate dst area over region that is unmasked and initialized
          do i=1,dst_dim
           if ((dst_mask(i) /= 0) .and. (FdstArray(i) /=UNINITVAL)) then
              grid2area(i) = FdstArray(i)*dst_area(i)
           endif
          enddo

          grid1areaXX = FsrcArray*src_area
          grid2areaXX = FdstArray*dst_area*dst_frac

          if ((j.eq.1).and.(sum(grid1area).gt.0)) then
              totCnt = totCnt + 2
            if (Verbose) then
              print *, ""
              print *, "Need smallest weight >= -0.8 ..."
            endif
            if (maxneg.lt.-8d-4) then
              success = .false.
              failCnt = failCnt + 1
              print *, "FAILED: max neg weight = ", maxneg
              print *, "Note: only a concern for conservative map"
            else
              if (Verbose) print*, "PASSED: min weight = ", maxneg
            endif

            if (Verbose) print *, "Need max positive weight < 1.00046 ..."
            if (maxpos-1.gt.4.6d-4) then
              success = .false.
              failCnt = failCnt + 1
              print *, "FAILED: max weight = ", maxpos
              print*, "Note: only a concern for conservative map"
            else
              if (Verbose) print*, "PASSED: max weight = ", maxpos
            endif
            if (Verbose) print *, ""
          endif

          if (j.eq.1) then
              totCnt = totCnt + 4
            if (Verbose) print *, "Need min(frac_a) > -1e-6 ..."
            if (srcfrac_min.lt.-1d-6) then
              success = .false.
              failCnt = failCnt + 1
              print*, "FAILED: min(frac_a) = ", srcfrac_min
            else
              if (Verbose) print*, "PASSED: min(frac_a) = ", srcfrac_min
            endif

            if (Verbose) print *, "Need max(frac_a) < 1+1e-6 ..."
            if (srcfrac_max.gt.1+1d-6) then
              success = .false.
              failCnt = failCnt + 1
              print*, "FAILED: max(frac_a) = ", srcfrac_max
            else
              if (Verbose) print*, "PASSED: max(frac_a) = ", srcfrac_max
            endif

            if (Verbose) print *, "Need min(frac_b) > -1e-6 ..."
            if (dstfrac_min.lt.-1d-6) then
              success = .false.
              failCnt = failCnt + 1
              print*, "FAILED: min(frac_b) = ", dstfrac_min
            else
              if (Verbose) print*, "PASSED: min(frac_b) = ", dstfrac_min
            endif

            if (Verbose) print *, "Need max(frac_b) < 1+1e-6 ..."
            if (dstfrac_max.gt.1+1d-6) then
              success = .false.
              failCnt = failCnt + 1
              print*, "FAILED: max(frac_b) = ", dstfrac_max
            else
              if (Verbose) print*, "PASSED: max(frac_b) = ", dstfrac_max
            endif

            if (Verbose) print *, ""
          endif

          if (Verbose) then
            print *, "  Test ", j
            print *, "+--------"
          endif

          select case (j)
            ! Note: totAreaBound only used for aave
            ! Using max values for bounds rather than statistics
            ! But mean and std-dev are provided anyway
            case (1)
              !   L1 error - avg = 5.411e-4, sigma = 8.211e-4
              !   L2 error - avg = 1.062e-3, sigma = 1.831e-3
              !   Area error - avg = 4.074e-13, sigma = 7.382e-13
              reltotErrorBound = 4.3e-3
              reltwoErrorBound = 1.5d-2
              totAreaBound = 5.3d-12
            case (2)
              !   L1 error - avg = 1.004e-7, sigma = 2.312e-7
              !   L2 error - avg = 4.934e-6, sigma = 8.707e-6
              !   Area error - avg = 5.472e-12, sigma = 1.383e-11
              reltotErrorBound = 2.0d-6
              reltwoErrorBound = 4.0d-5
              totAreaBound = 6.8d-11
            case (3)
              !   L1 error - avg = 2.203e-3, sigma = 2.857e-3
              !   L2 error - avg = 4.475e-3, sigma = 6.544e-3
              !   Area error - avg = 3.853e-13, sigma = 5.844e-13
              reltotErrorBound = 2.0d-2
              reltwoErrorBound = 5.1d-2
              totAreaBound = 2.6d-12
            case (4)
              !   L1 error - avg = 5.764e-3, sigma = 6.382e-3
              !   L2 error - avg = 1.177e-2, sigma = 1.482e-2
              !   Area error - avg = 4.878e-13, sigma = 9.986e-13
              reltotErrorBound = 3.4d-2
              reltwoErrorBound = 1.2d-1
              totAreaBound = 6.5d-12
            case (5)
              !   L1 error - avg = 1.686e-2, sigma = 1.751e-2
              !   L2 error - avg = 3.344e-2, sigma = 3.452e-2
              !   Area error - avg = 5.322e-13, sigma = 1.279e-12
              reltotErrorBound = 8.9d-2
              reltwoErrorBound = 1.4d-1
              totAreaBound = 1.1d-11
            case default
              print *, "Note: statistics only gathered for 5 tests... ", &
                       "no bounds for test ", j
              reltotErrorBound = UNINITVAL 
              reltwoErrorBound = UNINITVAL
              totAreaBound = UNINITVAL 
          end select

          totCnt = totCnt + 2
          if (Verbose) print *, "| Need L1 Error < ", reltotErrorBound, " ..."
          if (reltotError.gt.reltotErrorBound) then
            success = .false.
            failCnt = failCnt + 1
            if (Verbose) then
              print *, "| FAILED: L1 error = ", reltotError
            else
              print *, "FAILED: L1 error = ", reltotError, " in test ", j
            endif
          else
            if (Verbose) print *, "| PASSED: L1 error = ", reltotError
          endif

          if (Verbose) print *, "| Need L2 Error < ", reltwoErrorBound, " ..."
          if (reltwoError.gt.reltwoErrorBound) then
            success = .false.
            failCnt = failCnt + 1
            if (Verbose) then
              print *, "| FAILED: L2 error = ", reltwoError
            else
              print *, "FAILED: L2 error = ", reltwoError, " in test ", j
            endif
          else
            if (Verbose) print *, "| PASSED: L2 error = ", reltwoError
          endif
          if (sum(grid1area).gt.0) then
            totArea = abs(sum(grid2area)-sum(grid1area))
            totCnt = totCnt + 1

            if (Verbose) print *, "| Need area error < ", totAreaBound, " ..."
            if (totArea.gt.totAreaBound) then
              success = .false.
              failCnt = failCnt + 1
              if (Verbose) then
                print *, "FAILED: conservation error = ", totArea
              else
                print *, "FAILED: conservation error = ", totArea, &
                         " in test ", j
              endif
            else
              if (Verbose) print *, "| PASSED: Area error = ", totArea
            endif
      !      print *, " "
      !      print *, "reverse fracs  - Grid 1 area = ", sum(grid1areaXX(j,:))
      !      print *, "reverse fracs  - Grid 2 area = ", sum(grid2areaXX(j,:))
      !      print *, "reverse - Conservation error = ", abs(sum(grid2areaXX(j,:))-sum(grid1areaXX(j,:)))
          endif
          if (Verbose) then
            print *, "+--------"
            print *, " "
          end if
          deallocate(grid1area)
          deallocate(grid2area)
          deallocate(grid1areaXX)
          deallocate(grid2areaXX)
        else
          success = .false.
          print *, "ERROR: the test did not successfully map any values ", &
                   "from the source grid to the destination grid"
          exit
        endif
      enddo
      if (success) then
        print *, "All ", totCnt, " tests passed!"
      else
        print *, failCnt, " of ", totCnt, " tests failed. See above for details."
      endif
      deallocate(src_area)
      deallocate(src_frac)
      deallocate(dst_area)
      deallocate(FsrcArray)
      deallocate(FdstArray)
      deallocate(FdstArrayX)
    endif

      ! destroy and deallocate
        call ESMF_ArrayDestroy(srcArray, rc=status)
        call ESMF_ArrayDestroy(dstArray, rc=status)
        if (ESMF_LogFoundError(rcToCheck=status, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_Finalize()

contains

!***********************************************************************************
! Read in the grid dimensions info from the weights file.
! The weights file should have the source and destination grid information
! provided.
!***********************************************************************************
  subroutine NCFileInquire (wgtfile, title, src_dim, nxs, nys, dst_dim, nxd, nyd, localrc)
 
    character(ESMF_MAXSTR), intent(in)  :: wgtfile
    character(ESMF_MAXSTR), intent(out)  :: title
    integer, intent(out)                :: src_dim
    integer, intent(out)                :: nxs, nys
    integer, intent(out)                :: dst_dim
    integer, intent(out)                :: nxd, nyd
    integer, intent(out), optional      :: localrc
      
    integer :: ncstat,  nc_file_id,  nc_srcdim_id, nc_dstdim_id, srcdim, dstdim
    integer :: gdims(2), dim_ids(1)
    integer :: titleLen

    character(ESMF_MAXSTR) :: msg

      !-----------------------------------------------------------------
      ! open netcdf file
      !-----------------------------------------------------------------

      ncstat = nf90_open(wgtfile, NF90_NOWRITE, nc_file_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_open error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif

      !-----------------------------------------------------------------
      ! source grid dimensions
      !-----------------------------------------------------------------

      ncstat = nf90_inquire_attribute(nc_file_id, nf90_global, 'title', len=titleLen)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inquire_attribute error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      if(len(title) < titleLen) then
        print *, "Not enough space to put title."
        return
      end if
      ncstat = nf90_get_att(nc_file_id, nf90_global, "title", title)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_att error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      !-----------------------------------------------------------------
      ! source grid dimensions
      !-----------------------------------------------------------------

      ncstat = nf90_inq_dimid(nc_file_id, 'n_a', nc_srcdim_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inq_dimid error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      ncstat = nf90_inquire_dimension(nc_file_id, nc_srcdim_id, len=src_dim)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inquire_variable error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      nxs = src_dim
      nys = 1

      !-----------------------------------------------------------------
      ! destination grid dimensions
      !-----------------------------------------------------------------

      ncstat = nf90_inq_dimid(nc_file_id, 'n_b', nc_dstdim_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inq_dimid error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      ncstat = nf90_inquire_dimension(nc_file_id, nc_dstdim_id, len=dst_dim)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inquire_variable error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      nxd = dst_dim
      nyd = 1

      !------------------------------------------------------------------------
      !     get 2d grid sizes if we can, overwrite src_dim and dst_dim above
      !------------------------------------------------------------------------

      ncstat = nf90_inq_varid(nc_file_id, 'src_grid_dims', nc_srcdim_id)
      ncstat = nf90_inquire_variable(nc_file_id, nc_srcdim_id, dimids=dim_ids)
      if (ncstat.eq.0) then
        ncstat = nf90_inquire_dimension(nc_file_id, dim_ids(1), len=srcdim)
      else
        print*, "ERROR READING DIMENSION OF src_grid_dims"
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      if (srcdim.eq.2) then
        ncstat = nf90_get_var(ncid=nc_file_id, varid=nc_srcdim_id, values=gdims)
        nxs = gdims(1)
        nys = gdims(2)
      endif

      ncstat = nf90_inq_varid(nc_file_id, 'dst_grid_dims', nc_dstdim_id)
      ncstat = nf90_inquire_variable(nc_file_id, nc_dstdim_id, dimids=dim_ids)
      if (ncstat.eq.0) then
        ncstat = nf90_inquire_dimension(nc_file_id, dim_ids(1), len=dstdim)
      else
        print*, "ERROR READING DIMENSION OF dst_grid_dims"
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      if (dstdim.eq.2) then
        ncstat = nf90_get_var(ncid=nc_file_id, varid=nc_dstdim_id, values=gdims)
        nxd = gdims(1)
        nyd = gdims(2)
      end if

      !------------------------------------------------------------------------
      !     close input file
      !------------------------------------------------------------------------

      ncstat = nf90_close(nc_file_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_close error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      if(present(localrc)) localrc = ESMF_SUCCESS

  end subroutine NCFileInquire

!***********************************************************************************
! Read in the grid info from the weights file.
! The weights file should have the source and destination grid information
! provided.
!***********************************************************************************
  subroutine GridReadCoords (wgtfile, src_lat, src_lon, src_area, &
                             src_mask, src_frac, &
                             dst_lat, dst_lon, dst_area, dst_mask, &
                             dst_frac, localrc)
 
    character(ESMF_MAXSTR), intent(in)  :: wgtfile
    real(ESMF_KIND_R8), pointer         :: src_lat(:)
    real(ESMF_KIND_R8), pointer         :: src_lon(:)
    real(ESMF_KIND_R8), pointer         :: src_area(:)
    real(ESMF_KIND_R8), pointer         :: src_mask(:)
    real(ESMF_KIND_R8), pointer         :: src_frac(:)
    real(ESMF_KIND_R8), pointer         :: dst_lat(:)
    real(ESMF_KIND_R8), pointer         :: dst_lon(:)
    real(ESMF_KIND_R8), pointer         :: dst_area(:)
    real(ESMF_KIND_R8), pointer         :: dst_mask(:)
    real(ESMF_KIND_R8), pointer         :: dst_frac(:)
    integer, intent(out), optional      :: localrc
      
    integer :: ncstat,  nc_file_id
    integer :: nc_srcgridlat_id, nc_srcgridlon_id, &
               nc_dstgridlat_id, nc_dstgridlon_id, &
               nc_srcarea_id, nc_dstarea_id, &
               nc_srcmask_id, nc_dstmask_id, &
               nc_srcfrac_id, nc_dstfrac_id 
    integer :: unitsLen

    character(ESMF_MAXSTR) :: units, buffer
    character(ESMF_MAXSTR) :: msg

    real(ESMF_KIND_R8), parameter :: d2r = 3.141592653589793238/180

      !-----------------------------------------------------------------
      ! open netcdf file
      !-----------------------------------------------------------------

      ncstat = nf90_open(wgtfile, NF90_NOWRITE, nc_file_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_open error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
      endif

      !-----------------------------------------------------------------
      ! get the grid coordinates
      !-----------------------------------------------------------------

      ncstat = nf90_inq_varid(nc_file_id, 'yc_a', nc_srcgridlat_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inq_varid error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      ncstat = nf90_get_var(ncid=nc_file_id, varid=nc_srcgridlat_id, &
        values=src_lat)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_var error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      ! get the units of the grid coordinates
      ncstat = nf90_inquire_attribute(nc_file_id, nc_srcgridlat_id, 'units', &
        len=unitsLen)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inquire_attribute error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      if(len(units) < unitsLen) then
        print *, "Not enough space to get units."
        return
      endif
      ncstat = nf90_get_att(nc_file_id, nc_srcgridlat_id, "units", buffer)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_att error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      units = buffer(1:unitsLen)
      ! convert to radians if coordinates are in degrees
      if ((trim(units)==trim("degrees")).or.(trim(units)==trim("degrees north"))) then
        src_lat = src_lat*d2r
      else if (trim(units)/=trim("radians")) then
        write (msg, '(a,i4)') "- units are not 'degrees' or 'radians'"
        call ESMF_LogSetError(ESMF_RC_OBJ_BAD, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      ncstat = nf90_inq_varid(nc_file_id, 'xc_a', nc_srcgridlon_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inq_varid error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      ncstat = nf90_get_var(ncid=nc_file_id, varid=nc_srcgridlon_id, &
        values=src_lon)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_var error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      ! get the units of the grid coordinates
      ncstat = nf90_inquire_attribute(nc_file_id, nc_srcgridlon_id, 'units', &
        len=unitsLen)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inquire_attribute error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      if(len(units) < unitsLen) then
        print *, "Not enough space to get units."
        return
      endif
      ncstat = nf90_get_att(nc_file_id, nc_srcgridlon_id, "units", buffer)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_att error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      units = buffer(1:unitsLen)
      ! convert to radians if coordinates are in degrees
      if ((trim(units)==trim("degrees")).or.(trim(units)==trim("degrees east"))) then
        src_lon = src_lon*d2r
      else if (trim(units)/=trim("radians")) then
        write (msg, '(a,i4)') "- units are not 'degrees' or 'radians'"
        call ESMF_LogSetError(ESMF_RC_OBJ_BAD, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      !-----------------------------------------------------------------
      ! get the grid coordinates
      !-----------------------------------------------------------------
      ncstat = nf90_inq_varid(nc_file_id, 'yc_b', nc_dstgridlat_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inq_varid error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      ncstat = nf90_get_var(ncid=nc_file_id, varid=nc_dstgridlat_id, &
        values=dst_lat)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_var error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      ! get the units of the grid coordinates
      ncstat = nf90_inquire_attribute(nc_file_id, nc_dstgridlat_id, 'units', &
        len=unitsLen)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inquire_attribute error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      if(len(units) < unitsLen) then
        print *, "Not enough space to get units."
        return
      endif
      ncstat = nf90_get_att(nc_file_id, nc_dstgridlat_id, "units", buffer)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_att error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      units = buffer(1:unitsLen)
      ! convert to radians if coordinates are in degrees
      if ((trim(units)==trim("degrees")).or.(trim(units)==trim("degrees north"))) then
        dst_lat = dst_lat*d2r
      else if (trim(units)/=trim("radians")) then
        write (msg, '(a,i4)') "- units are not 'degrees' or 'radians'"
        call ESMF_LogSetError(ESMF_RC_OBJ_BAD, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      ncstat = nf90_inq_varid(nc_file_id, 'xc_b', nc_dstgridlon_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inq_varid error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      ncstat = nf90_get_var(ncid=nc_file_id, varid=nc_dstgridlon_id, &
        values=dst_lon)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_var error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      ! get the units of the grid coordinates
      ncstat = nf90_inquire_attribute(nc_file_id, nc_dstgridlon_id, 'units', &
        len=unitsLen)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inquire_attribute error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      if(len(units) < unitsLen) then
        print *, "Not enough space to get units."
        return
      endif
      ncstat = nf90_get_att(nc_file_id, nc_dstgridlon_id, "units", buffer)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_att error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      units = buffer(1:unitsLen)
      ! convert to radians if coordinates are in degrees
      if ((trim(units)==trim("degrees")).or.(trim(units)==trim("degrees east"))) then
        dst_lon = dst_lon*d2r
      else if (trim(units)/=trim("radians")) then
        write (msg, '(a,i4)') "- units are not 'degrees' or 'radians'"
        call ESMF_LogSetError(ESMF_RC_OBJ_BAD, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      !-----------------------------------------------------------------
      ! get the grid areas
      !-----------------------------------------------------------------
      ncstat = nf90_inq_varid(nc_file_id, 'area_a', nc_srcarea_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inq_varid error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      ncstat = nf90_get_var(ncid=nc_file_id, varid=nc_srcarea_id, &
        values=src_area)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_var error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      ncstat = nf90_inq_varid(nc_file_id, 'area_b', nc_dstarea_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inq_varid error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      ncstat = nf90_get_var(ncid=nc_file_id, varid=nc_dstarea_id, &
        values=dst_area)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_var error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      !-----------------------------------------------------------------
      ! get the grid masks
      !-----------------------------------------------------------------
      ncstat = nf90_inq_varid(nc_file_id, 'mask_a', nc_srcmask_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inq_varid error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      ncstat = nf90_get_var(ncid=nc_file_id, varid=nc_srcmask_id, &
        values=src_mask)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_var error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      ncstat = nf90_inq_varid(nc_file_id, 'mask_b', nc_dstmask_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inq_varid error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      ncstat = nf90_get_var(ncid=nc_file_id, varid=nc_dstmask_id, &
        values=dst_mask)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_var error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      !-----------------------------------------------------------------
      ! get the grid fracs
      !-----------------------------------------------------------------
      ncstat = nf90_inq_varid(nc_file_id, 'frac_a', nc_srcfrac_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inq_varid error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      ncstat = nf90_get_var(ncid=nc_file_id, varid=nc_srcfrac_id, &
        values=src_frac)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_var error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      ncstat = nf90_inq_varid(nc_file_id, 'frac_b', nc_dstfrac_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_inq_varid error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif
      ncstat = nf90_get_var(ncid=nc_file_id, varid=nc_dstfrac_id, &
        values=dst_frac)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_get_var error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      !------------------------------------------------------------------------
      !     close input file
      !------------------------------------------------------------------------

      ncstat = nf90_close(nc_file_id)
      if(ncstat /= 0) then
        write (msg, '(a,i4)') "- nf90_close error:", ncstat
        call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
          line=__LINE__, file=__FILE__ , rcToReturn=rc)
        return
      endif

      if(present(localrc)) localrc = ESMF_SUCCESS

  end subroutine GridReadCoords

!------------------------------------------------------------------------------

   subroutine ESMF_FieldRegridReadSCRIPFileP(remapFile, factorList, factorIndexList, rc)
!------------------------------------------------------------------------
!     call arguments
!------------------------------------------------------------------------

     character (ESMF_MAXSTR), intent(in)         :: remapFile
     real(ESMF_KIND_R8), pointer                 :: factorList(:)
     integer, pointer                            :: factorIndexList(:,:)
     integer, intent(out), optional              :: rc

!------------------------------------------------------------------------
!     local variables
!------------------------------------------------------------------------

     integer :: ncstat,  nc_file_id,  nc_numlinks_id, nc_numwgts_id, &
     nc_dstgrdadd_id, nc_srcgrdadd_id, nc_rmpmatrix_id

     integer :: num_links, num_wts

     character (ESMF_MAXSTR) :: nm, msg

     integer, allocatable  :: address(:), localSize(:), localOffset(:)
     type(ESMF_VM)         :: vm
     integer               :: i, localpet, npet, nlinksPPet, FlocalPet

     ! get lpe number
     call ESMF_VMGetCurrent(vm, rc=rc)
     if(rc /= ESMF_SUCCESS) then
       write (msg, '(a,i4)') "- failed to get current vm", ncstat
       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
     endif

     call ESMF_VMGet(vm, localPet=localPet, petCount=npet, rc=rc)
     if(rc /= ESMF_SUCCESS) then
       write (msg, '(a,i4)') "- failed to get current vm", ncstat
       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
     endif

    ! set the localPet to localPet +1 for fortran array indices
    FlocalPet = localPet+1

!-----------------------------------------------------------------------
!     open remap file and read meta data
!-----------------------------------------------------------------------
     !-----------------------------------------------------------------
     ! open netcdf file
     !-----------------------------------------------------------------

     ncstat = nf90_open(remapFile, NF90_NOWRITE, nc_file_id)
     if(ncstat /= 0) then
       write (msg, '(a,i4)') "- nf90_open error:", ncstat
       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
     endif

!-----------------------------------------------------------------------
! read source grid meta data for consistency check
!-----------------------------------------------------------------------
     !-----------------------------------------------------------------
     ! number of address pairs in the remappings
     !-----------------------------------------------------------------

     ncstat = nf90_inq_dimid(nc_file_id, 'n_s', nc_numlinks_id)
     if(ncstat /= 0) then
       write (msg, '(a,i4)') "- nf90_inq_dimid error:", ncstat
       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
     endif
     ncstat = nf90_inquire_dimension(nc_file_id, nc_numlinks_id, len=num_links)
     if(ncstat /= 0) then
       write (msg, '(a,i4)') "- nf90_inquire_dimension error:", ncstat
       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
     endif

     !-----------------------------------------------------------------
     ! number of weights per point/order of interpolation method
!     ncstat = nf90_inq_dimid(nc_file_id, 'num_wgts', nc_numwgts_id)
!     if(ncstat /= 0) then
!       write (msg, '(a,i4)') "- nf90_inq_dimid error:", ncstat
!       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
!         line=__LINE__, file=__FILE__, rcToReturn=rc)
!       return
!     endif
!     ncstat = nf90_inquire_dimension(nc_file_id, nc_numwgts_id, len=num_wts)
!     if(ncstat /= 0) then
!       write (msg, '(a,i4)') "- nf90_inquire_dimension error:", ncstat
!       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
!         line=__LINE__, file=__FILE__, rcToReturn=rc)
!       return
!     endif

     ! split the input data between PETs
     ! allocate factorList and factorIndexList
#if 0
     allocate( localSize(npet), localOffset(npet) )
     nlinksPPet = num_links/npet
     localSize(:) = nlinksPPet

     do i = 1, npet
         localOffset(i) = 1 + (i-1)*nlinksPPet
     enddo
     localSize(npet) = nlinksPPet+MOD(num_links, npet)

     allocate( factorIndexList(2,localSize(FlocalPet)) )
     allocate( factorList(localSize(FlocalPet)) )
#endif
     allocate( factorIndexList(2,num_links) )
     allocate( factorList(num_links) )
     !-----------------------------------------------------------------
     ! source addresses for weights
     !-----------------------------------------------------------------

     allocate( address(num_links) )
     ncstat = nf90_inq_varid(nc_file_id, 'col', nc_srcgrdadd_id)
     if(ncstat /= 0) then
       write (msg, '(a,i4)') "- nf90_inq_varid error:", ncstat
       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
     endif
     ncstat = nf90_get_var(ncid=nc_file_id, varid=nc_srcgrdadd_id, &
       values=address)
     if(ncstat /= 0) then
       write (msg, '(a,i4)') "- nf90_get_var error:", ncstat
       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
     endif
     factorIndexList(1,:) = address

     !-----------------------------------------------------------------
     ! destination addresss for weights
     !-----------------------------------------------------------------

     ncstat = nf90_inq_varid(nc_file_id, 'row', nc_dstgrdadd_id)
     if(ncstat /= 0) then
       write (msg, '(a,i4)') "- nf90_inq_varid error:", ncstat
       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
     endif
     ncstat = nf90_get_var(nc_file_id, nc_dstgrdadd_id, address)
     if(ncstat /= 0) then
       write (msg, '(a,i4)') "- nf90_get_var error:", ncstat
       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
     endif
     factorIndexList(2,:) = address
     deallocate( address )

     !-----------------------------------------------------------------
     !     read all variables
     !-----------------------------------------------------------------

     ncstat = nf90_inq_varid(nc_file_id, 'S', nc_rmpmatrix_id)
     if(ncstat /= 0) then
       write (msg, '(a,i4)') "- nf90_inq_varid error:", ncstat
       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
     endif
       write (msg, '(a,i4)') "- nf90_get_var error:", ncstat
     ncstat = nf90_get_var(nc_file_id, nc_rmpmatrix_id, factorList)
!       localOffset(FlocalPet), localSize(FlocalPet))
     if(ncstat /= 0) then
       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
     endif
!------------------------------------------------------------------------
!     close input file
!------------------------------------------------------------------------

     ncstat = nf90_close(nc_file_id)
     if(ncstat /= 0) then
       write (msg, '(a,i4)') "- nf90_close error:", ncstat
       call ESMF_LogSetError(ESMF_RC_SYS, msg=msg, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
     endif

     if(present(rc)) rc = ESMF_SUCCESS

   end subroutine ESMF_FieldRegridReadSCRIPFileP

!------------------------------------------------------------------------------

end program OfflineTester
