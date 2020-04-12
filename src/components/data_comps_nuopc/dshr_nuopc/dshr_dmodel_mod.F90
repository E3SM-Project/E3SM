module dshr_dmodel_mod

  ! Obtain the model domain and the stream domain for each stream
  ! For the model domain 
  !  - if single column - read in the data model domain file - and find the nearest neighbor
  !  - if not single column - will obtain it directly from the mesh input but will still need 
  !    to read in the domain file to obtain the global nx and ny that need to be passed as scalar
  !    data back to the mediator

  use ESMF
  use shr_sys_mod
  use shr_kind_mod  , only : R8=>SHR_KIND_R8
  use shr_kind_mod  , only : CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_kind_mod  , only : CX=>SHR_KIND_CX, CXX=>SHR_KIND_CXX
  use shr_log_mod   , only : logunit => shr_log_Unit
  use shr_mpi_mod   , only : shr_mpi_bcast
  use shr_const_mod , only : shr_const_cDay
  use shr_ncread_mod, only : shr_ncread_open, shr_ncread_close, shr_ncread_varDimSizes, shr_ncread_tField
  use shr_ncread_mod, only : shr_ncread_domain, shr_ncread_vardimsizes
  use shr_ncread_mod, only : shr_ncread_varexists, shr_ncread_vardimnum,  shr_ncread_field4dG
  use shr_ncread_mod, only : shr_ncread_open, shr_ncread_close, shr_ncread_varDimSizes
  use shr_ncread_mod, only : shr_ncread_tField
  use dshr_methods_mod , only : chkerr
  use shr_map_mod
  use dshr_stream_mod
  use mct_mod
  use perf_mod
  use pio

  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public  :: shr_dmodel_gGridCompare
  public  :: shr_dmodel_mapSet
  public  :: shr_dmodel_readLBUB

  integer,parameter,public :: CompareXYabs      = 1 ! X,Y  relative error
  integer,parameter,public :: CompareXYrel      = 2 ! X,Y  absolute error
  integer,parameter,public :: CompareAreaAbs    = 3 ! area relative error
  integer,parameter,public :: CompareAreaRel    = 4 ! area absolute error
  integer,parameter,public :: CompareMaskIdent  = 5 ! masks are identical
  integer,parameter,public :: CompareMaskZeros  = 6 ! masks have same zeros
  integer,parameter,public :: CompareMaskSubset = 7 ! mask is subset of other

  integer,parameter,public :: CompareXYabsMask  = 101 ! X,Y  relative error
  integer,parameter,public :: iotype_std_netcdf = -99 ! non pio option

  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine shr_dmodel_readLBUB(stream,pio_subsystem,pio_iotype,pio_iodesc,&
       mDate,mSec,mpicom,gsMap, &
       avLB,mDateLB,mSecLB,avUB,mDateUB,mSecUB,avFile,readMode, &
       newData,rmOldFile,istr)

    !-------------------------------------------------------------------------
    ! Read LB and UB of stream data
    !-------------------------------------------------------------------------

    !----- arguments -----
    type(shr_stream_streamType) ,intent(inout) :: stream
    type(iosystem_desc_t)       ,intent(inout), target :: pio_subsystem
    integer                     ,intent(in)    :: pio_iotype
    type(io_desc_t)             ,intent(inout) :: pio_iodesc
    integer                     ,intent(in)    :: mDate  ,mSec
    integer                     ,intent(in)    :: mpicom
    type(mct_gsMap)             ,intent(in)    :: gsMap
    type(mct_aVect)             ,intent(inout) :: avLB
    integer                     ,intent(inout) :: mDateLB,mSecLB
    type(mct_aVect)             ,intent(inout) :: avUB
    integer                     ,intent(inout) :: mDateUB,mSecUB
    type(mct_aVect)             ,intent(inout) :: avFile
    character(len=*)            ,intent(in)    :: readMode
    logical                     ,intent(out)   :: newData
    logical          ,optional  ,intent(in)    :: rmOldFile
    character(len=*) ,optional  ,intent(in)    :: istr

    !----- local -----
    integer           :: my_task, master_task
    integer           :: ierr       ! error code
    integer           :: rCode      ! return code
    logical           :: localCopy,fileexists
    integer           :: ivals(6)   ! bcast buffer
    integer           :: oDateLB,oSecLB,dDateLB,oDateUB,oSecUB,dDateUB
    real(R8)          :: rDateM,rDateLB,rDateUB  ! model,LB,UB dates with fractional days
    integer           :: n_lb, n_ub
    character(CL)     :: fn_lb,fn_ub,fn_next,fn_prev
    character(CL)     :: path
    character(len=32) :: lstr
    real(R8)          :: spd
    character(*), parameter :: subname = '(shr_dmodel_readLBUB) '
    character(*), parameter :: F00   = "('(shr_dmodel_readLBUB) ',8a)"
    character(*), parameter :: F01   = "('(shr_dmodel_readLBUB) ',a,5i8)"
    !-------------------------------------------------------------------------

    lstr = 'shr_dmodel_readLBUB'
    if (present(istr)) then
       lstr = trim(istr)
    endif

    call t_startf(trim(lstr)//'_setup')
    call MPI_COMM_RANK(mpicom,my_task,ierr)
    master_task = 0
    spd = shr_const_cday

    newData = .false.
    n_lb = -1
    n_ub = -1
    fn_lb = 'undefinedlb'
    fn_ub = 'undefinedub'

    oDateLB = mDateLB
    oSecLB  = mSecLB
    oDateUB = mDateUB
    oSecUB  = mSecUB

    rDateM  = real(mDate  ,R8) + real(mSec  ,R8)/spd
    rDateLB = real(mDateLB,R8) + real(mSecLB,R8)/spd
    rDateUB = real(mDateUB,R8) + real(mSecUB,R8)/spd
    call t_stopf(trim(lstr)//'_setup')

    if (rDateM < rDateLB .or. rDateM > rDateUB) then
       call t_startf(trim(lstr)//'_fbound')
       if (my_task == master_task) then
          call shr_stream_findBounds(stream,mDate,mSec,                 &
               ivals(1),dDateLB,ivals(2),ivals(5),fn_lb, &
               ivals(3),dDateUB,ivals(4),ivals(6),fn_ub  )
          call shr_stream_getFilePath(stream,path)
       endif
       call t_stopf(trim(lstr)//'_fbound')
       call t_startf(trim(lstr)//'_bcast')

       !    --- change 4 bcasts to a single bcast and copy for performance ---
       !     call shr_mpi_bcast(mDateLB,mpicom)
       !     call shr_mpi_bcast(mSecLB,mpicom)
       !     call shr_mpi_bcast(mDateUB,mpicom)
       !     call shr_mpi_bcast(mSecUB,mpicom)
       call shr_mpi_bcast(stream%calendar,mpicom)
       call shr_mpi_bcast(ivals,mpicom)
       mDateLB = ivals(1)
       mSecLB  = ivals(2)
       mDateUB = ivals(3)
       mSecUB  = ivals(4)
       n_lb    = ivals(5)
       n_ub    = ivals(6)
       call t_stopf(trim(lstr)//'_bcast')
    endif

    if (mDateLB /= oDateLB .or. mSecLB /= oSecLB) then
       newdata = .true.
       if (mDateLB == oDateUB .and. mSecLB == oSecUB) then
          call t_startf(trim(lstr)//'_LB_copy')
          avLB%rAttr(:,:) = avUB%rAttr(:,:)
          call t_stopf(trim(lstr)//'_LB_copy')
       else

          select case(readMode)
          case ('single')
             call shr_dmodel_readstrm(stream, pio_subsystem, pio_iotype, pio_iodesc, gsMap, avLB, mpicom, &
                  path, fn_lb, n_lb,istr=trim(lstr)//'_LB', boundstr = 'lb')
          case ('full_file')
             call shr_dmodel_readstrm_fullfile(stream, pio_subsystem, pio_iotype, &
                  gsMap, avLB, avFile, mpicom, &
                  path, fn_lb, n_lb,istr=trim(lstr)//'_LB', boundstr = 'lb')
          case default
             write(logunit,F00) "ERROR: Unsupported readmode : ", trim(readMode)
             call shr_sys_abort(subName//"ERROR: Unsupported readmode: "//trim(readMode))
          end select
       endif
    endif

    if (mDateUB /= oDateUB .or. mSecUB /= oSecUB) then
       newdata = .true.

       select case(readMode)
       case ('single')
          call shr_dmodel_readstrm(stream, pio_subsystem, pio_iotype, pio_iodesc, gsMap, avUB, mpicom, &
               path, fn_ub, n_ub,istr=trim(lstr)//'_UB', boundstr = 'ub')
       case ('full_file')
          call shr_dmodel_readstrm_fullfile(stream, pio_subsystem, pio_iotype, &
               gsMap, avUB, avFile, mpicom, &
               path, fn_ub, n_ub,istr=trim(lstr)//'_UB', boundstr = 'ub')
       case default
          write(logunit,F00) "ERROR: Unsupported readmode : ", trim(readMode)
          call shr_sys_abort(subName//"ERROR: Unsupported readmode: "//trim(readMode))
       end select

    endif

    call t_startf(trim(lstr)//'_filemgt')
    !--- determine previous & next data files in list of files ---
    if (my_task == master_task .and. newdata) then
       call shr_stream_getFilePath(stream,path)
       call shr_stream_getPrevFileName(stream,fn_lb,fn_prev,path)
       call shr_stream_getNextFileName(stream,fn_ub,fn_next,path)
       inquire(file=trim(fn_next),exist=fileExists)
       if ( trim(fn_next) == "unknown" .or. fileExists) then
          ! do nothing
       end if
    endif
    call t_stopf(trim(lstr)//'_filemgt')

  end subroutine shr_dmodel_readLBUB

  !===============================================================================

  subroutine shr_dmodel_readstrm(stream, pio_subsystem, pio_iotype, pio_iodesc, gsMap, av, mpicom, &
       path, fn, nt, istr, boundstr)

    !----- arguments -----
    type(shr_stream_streamType) ,intent(inout)         :: stream
    type(iosystem_desc_t)       ,intent(inout), target :: pio_subsystem
    integer                     ,intent(in)            :: pio_iotype
    type(io_desc_t)             ,intent(inout)         :: pio_iodesc
    type(mct_gsMap)             ,intent(in)            :: gsMap
    type(mct_aVect)             ,intent(inout)         :: av
    integer                     ,intent(in)            :: mpicom
    character(len=*)            ,intent(in)            :: path
    character(len=*)            ,intent(in)            :: fn
    integer                     ,intent(in)            :: nt
    character(len=*),optional   ,intent(in)            :: istr
    character(len=*),optional   ,intent(in)            :: boundstr

    !----- local -----
    integer              :: my_task
    integer              :: master_task
    integer              :: ierr
    logical              :: localCopy,fileexists
    type(mct_avect)      :: avG
    integer              :: gsize,nx,ny,nz
    integer              :: k
    integer              :: fid
    integer              :: rCode      ! return code
    real(R8),allocatable :: data2d(:,:)
    real(R8),allocatable :: data3d(:,:,:)
    logical              :: d3dflag
    character(CL)        :: fileName
    character(CL)        :: sfldName
    type(mct_avect)      :: avtmp
    character(len=32)    :: lstr
    character(len=32)    :: bstr
    logical              :: fileopen
    character(CL)        :: currfile
    integer              :: ndims
    integer,pointer      :: dimid(:)
    type(file_desc_t)    :: pioid
    type(var_desc_t)     :: varid
    integer(kind=pio_offset_kind) :: frame
    character(*), parameter :: subname = '(shr_dmodel_readstrm) '
    character(*), parameter :: F00   = "('(shr_dmodel_readstrm) ',8a)"
    character(*), parameter :: F01   = "('(shr_dmodel_readstrm) ',a,5i8)"
    character(*), parameter :: F02   = "('(shr_dmodel_readstrm) ',2a,i8)"
    !-------------------------------------------------------------------------------

    lstr = 'shr_dmodel_readstrm'
    if (present(istr)) then
       lstr = trim(istr)
    endif

    bstr = ''
    if (present(boundstr)) then
       bstr = trim(boundstr)
    endif

    call t_barrierf(trim(lstr)//'_BARRIER',mpicom)
    call t_startf(trim(lstr)//'_setup')
    call MPI_COMM_RANK(mpicom,my_task,ierr)
    master_task = 0

    gsize = mct_gsmap_gsize(gsMap)

    if (my_task == master_task) then
       fileName = trim(path)//trim(fn)
       inquire(file=trim(fileName),exist=fileExists)
       if (.not. fileExists) then
          write(logunit,F00) "ERROR: file does not exist: ", trim(fileName)
          call shr_sys_abort(subName//"ERROR: file does not exist: "//trim(fileName))
       end if
    endif

    if (my_task == master_task) then
       call shr_stream_getFileFieldName(stream,1,sfldName)
    endif

    call t_stopf(trim(lstr)//'_setup')

    if (pio_iotype == iotype_std_netcdf) then

       call t_startf(trim(lstr)//'_readcdf')
       if (my_task == master_task) then
          call shr_ncread_varDimSizes(trim(fileName),trim(sfldName),nx,ny,nz)
          if (gsize == nx*ny) then
             d3dflag = .false.
             allocate(data2d(nx,ny))
          elseif (gsize == nx*ny*nz) then
             d3dflag = .true.
             allocate(data3d(nx,ny,nz))
          else
             write(logunit,F01) "ERROR in data sizes ",nx,ny,nz,gsize
             call shr_sys_abort(subname//"ERROR in data sizes")
          endif
          call mct_aVect_init(avG,av,gsize)
          call shr_ncread_open(trim(fileName),fid,rCode)
          do k = 1,mct_aVect_nRAttr(av)
             call shr_stream_getFileFieldName(stream,k,sfldName)
             if (d3dflag) then
                call shr_ncread_tField(fileName,nt,sfldName,data3d,fidi=fid,rc=rCode)
                avG%rAttr(k,:) = reshape(data3d, (/gsize/))
             else
                call shr_ncread_tField(fileName,nt,sfldName,data2d,fidi=fid,rc=rCode)
                avG%rAttr(k,:) = reshape(data2d, (/gsize/))
             endif
          enddo
          call shr_ncread_close(fid,rCode)
          if (d3dflag) then
             deallocate(data3d)
          else
             deallocate(data2d)
          endif
       endif
       call t_stopf(trim(lstr)//'_readcdf')
       call t_barrierf(trim(lstr)//'_scatter'//'_BARRIER',mpicom)
       call t_startf(trim(lstr)//'_scatter')
       call mct_aVect_scatter(avG,avtmp,gsMap,master_task,mpicom)
       call mct_aVect_copy(avtmp,av)
       if (my_task == master_task) call mct_aVect_clean(avG)
       call mct_aVect_clean(avtmp)
       call t_stopf(trim(lstr)//'_scatter')

    else

       call t_startf(trim(lstr)//'_readpio')
       call shr_mpi_bcast(sfldName,mpicom,'sfldName')
       call shr_mpi_bcast(filename,mpicom,'filename')

       call shr_stream_getCurrFile(stream,fileopen=fileopen,currfile=currfile,currpioid=pioid)

       if (fileopen .and. currfile==filename) then
          ! don't reopen file, all good
       else
          ! otherwise close the old file if open and open new file
          if (fileopen) then
             if (my_task == master_task) then
                write(logunit,F00) 'close  : ',trim(currfile)
                call shr_sys_flush(logunit)
             endif
             call pio_closefile(pioid)
          endif
          if (my_task == master_task) then
             write(logunit,F00) 'open   : ',trim(filename)
             call shr_sys_flush(logunit)
          endif

          rcode = pio_openfile(pio_subsystem, pioid, pio_iotype, trim(filename), pio_nowrite)
          call shr_stream_setCurrFile(stream,fileopen=.true.,currfile=trim(filename),currpioid=pioid)
       endif

       if (my_task == master_task) then
          write(logunit,F02) 'file ' // trim(bstr) //': ',trim(filename),nt
          call shr_sys_flush(logunit)
       endif

       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)

       rcode = pio_inq_varid(pioid,trim(sfldName),varid)
       rcode = pio_inq_varndims(pioid, varid, ndims)
       allocate(dimid(ndims))
       rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
       if (ndims >= 1) rcode = pio_inq_dimlen(pioid, dimid(1), nx)
       if (ndims >= 2) rcode = pio_inq_dimlen(pioid, dimid(2), ny)
       if (ndims >= 3) rcode = pio_inq_dimlen(pioid, dimid(3), nz)
       deallocate(dimid)
       if (gsize /= nx*ny .and. gsize /= nx*ny*nz) then
          write(logunit,F01) "ERROR in data sizes ",nx,ny,nz,gsize
          call shr_sys_abort(subname//"ERROR in data sizes")
       endif

       do k = 1,mct_aVect_nRAttr(av)
          if (my_task == master_task) then
             call shr_stream_getFileFieldName(stream,k,sfldName)
          endif
          call shr_mpi_bcast(sfldName,mpicom,'sfldName')
          rcode = pio_inq_varid(pioid,trim(sfldName),varid)
          frame = nt
          call pio_setframe(pioid,varid,frame)
          call pio_read_darray(pioid,varid,pio_iodesc,av%rattr(k,:),rcode)
       enddo

       call t_stopf(trim(lstr)//'_readpio')

    endif

  end subroutine shr_dmodel_readstrm

  !===============================================================================

  subroutine shr_dmodel_readstrm_fullfile(stream, pio_subsystem, pio_iotype, &
       gsMap, av, avFile, mpicom, &
       path, fn, nt, istr, boundstr)

    !----- arguments -----
    type      (shr_stream_streamType) ,intent(inout)         :: stream
    type      (iosystem_desc_t)       ,intent(inout), target :: pio_subsystem
    integer                       ,intent(in)            :: pio_iotype
    type      (mct_gsMap)             ,intent(in)            :: gsMap
    type      (mct_aVect)             ,intent(inout)         :: av
    type      (mct_aVect)             ,intent(inout)         :: avFile
    integer                       ,intent(in)            :: mpicom
    character (len=*)                 ,intent(in)            :: path
    character (len=*)                 ,intent(in)            :: fn
    integer                       ,intent(in)            :: nt
    character (len=*)                 ,intent(in) ,optional  :: istr
    character(len=*)                  ,intent(in) ,optional  :: boundstr

    !----- local -----
    integer                       :: my_task
    integer                       :: master_task
    integer                       :: ierr
    logical                       :: localCopy,fileexists
    integer                       :: gsize,nx,ny,nz
    integer                       :: k
    integer                       :: rCode   ! return code
    character(CL)                 :: fileName
    character(CL)                 :: sfldName
    character(len=32)             :: lstr
    character(len=32)             :: bstr
    logical                       :: fileopen
    character(CL)                 :: currfile
    character(CXX)                :: fldList ! list of fields
    integer                       :: ndims
    integer,pointer               :: dimid(:)
    type(file_desc_t)             :: pioid
    type(var_desc_t)              :: varid
    integer(kind=pio_offset_kind) :: frame
    type(io_desc_t)               :: pio_iodesc_local
    integer                       :: avFile_beg, avFile_end
    integer                       :: lsize, cnt,m,n
    integer, allocatable          :: count(:), compDOF(:)
    integer, pointer,dimension(:) :: gsmOP   ! gsmap ordered points
    character(*), parameter :: subname = ' (shr_dmodel_readstrm_fullfile) '
    character(*), parameter :: F00   = "(' (shr_dmodel_readstrm_fullfile) ',8a)"
    character(*), parameter :: F01   = "(' (shr_dmodel_readstrm_fullfile) ',a,5i8)"
    character(*), parameter :: F02   = "(' (shr_dmodel_readstrm_fullfile) ',2a,2i8)"
    !-------------------------------------------------------------------------------

    lstr = 'shr_dmodel_readstrm_fullfile'
    if (present(istr)) then
       lstr = trim(istr)
    endif

    bstr = ''
    if (present(boundstr)) then
       bstr = trim(boundstr)
    endif

    call t_barrierf(trim(lstr)//'_BARRIER',mpicom)
    call t_startf(trim(lstr)//'_setup')
    call MPI_COMM_RANK(mpicom,my_task,ierr)
    master_task = 0

    gsize = mct_gsmap_gsize(gsMap)
    lsize = mct_gsmap_lsize(gsMap,mpicom)

    if (my_task == master_task) then
       fileName = trim(path) // trim(fn)
       inquire(file=trim(fileName),exist=fileExists)
       if (.not. fileExists) then
          write(logunit,F00) "ERROR: file does not exist: ", trim(fileName)
          call shr_sys_abort(subName//"ERROR: file does not exist: "//trim(fileName))
       end if
    endif

    if (my_task == master_task) then
       call shr_stream_getFileFieldName(stream,1,sfldName)
    endif

    call t_stopf(trim(lstr)//'_setup')

    if (pio_iotype == iotype_std_netcdf) then

       write(logunit,F01) "shr_dmodel_readstrm_fullfile: not supported for iotype_std_netcdf"
       call shr_sys_abort(subname//"ERROR extend shr_dmodel_readstrm_fullfile")

    else

       call t_startf(trim(lstr)//'_readpio')
       call shr_mpi_bcast(sfldName,mpicom,'sfldName')
       call shr_mpi_bcast(filename,mpicom,'filename')

       call shr_stream_getCurrFile(stream,fileopen=fileopen,currfile=currfile,currpioid=pioid)

       if (fileopen .and. currfile==filename) then
          ! don't reopen file, all good
       else
          ! otherwise close the old file if open, open the new file,
          ! and read all time slices of a temporal dataset within the new file.
          if (fileopen) then
             if (my_task == master_task) then
                write(logunit,F00) 'close  : ',trim(currfile)
                call shr_sys_flush(logunit)
             endif
             call pio_closefile(pioid)
          endif
          if (my_task == master_task) then
             write(logunit,F00) 'open   : ',trim(filename)
             call shr_sys_flush(logunit)
          endif
          rcode = pio_openfile(pio_subsystem, pioid, pio_iotype, trim(filename), pio_nowrite)
          call shr_stream_setCurrFile(stream,fileopen=.true.,currfile=trim(filename),currpioid=pioid)

          call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)

          rcode = pio_inq_varid(pioid,trim(sfldName),varid)
          rcode = pio_inq_varndims(pioid, varid, ndims)
          allocate(dimid(ndims))

          rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))

          nx = 1
          ny = 1
          nz = 1
          if (ndims >= 1) rcode = pio_inq_dimlen(pioid, dimid(1), nx)
          if (ndims >= 2) rcode = pio_inq_dimlen(pioid, dimid(2), ny)
          if (ndims >= 3) rcode = pio_inq_dimlen(pioid, dimid(3), nz)
          deallocate(dimid)

          if (gsize /= nx*ny) then
             write(logunit,F01) "ERROR in data sizes ",nx,ny,nz,gsize
             call shr_sys_abort(subname//"ERROR in data sizes")
          endif

          if (my_task == master_task) then
             call shr_stream_getModelFieldList(stream,fldList)
          endif
          call shr_mpi_bcast(fldList,mpicom)

          call mct_avect_clean(avFile)
          call mct_aVect_init(avFile,rlist=fldList,lsize=nx*ny*nz)

          call mct_gsmap_orderedPoints(gsMap,my_task,gsmOP)

          allocate(count(3))
          allocate(compDOF(lsize*nz),stat=rcode)
          if (rcode /= 0) call shr_sys_abort(subname//"ERROR insufficient memory")

          count(1) = nx
          count(2) = ny
          count(3) = nz

          if (my_task == master_task) then
             write(logunit,F02) 'file ' // trim(bstr) //': ',trim(filename),1,nz
             call shr_sys_flush(logunit)
          endif

          ! Create a 3D MCT component DOF corresponding to "2D(=gsmOP) x nz"
          cnt = 0
          do n = 1,nz
             do m = 1,lsize
                cnt = cnt + 1
                compDOF(cnt) = (n-1)*gsize + gsmOP(m)
             enddo
          enddo

          ! Initialize the decomposition
          call pio_initdecomp(pio_subsystem, pio_double, count, compDOF, pio_iodesc_local)

          ! For each attribute, read all frames in one go
          frame = 1
          do k = 1, mct_aVect_nRAttr(avFile)
             if (my_task == master_task) then
                call shr_stream_getFileFieldName(stream,k,sfldName)
             endif
             call shr_mpi_bcast(sfldName,mpicom,'sfldName')
             rcode = pio_inq_varid(pioid,trim(sfldName),varid)

             call pio_setframe(pioid,varid,frame)

             call pio_read_darray(pioid, varid, pio_iodesc_local, avFile%rattr(k,:), rcode)
          enddo

          call pio_freedecomp(pio_subsystem, pio_iodesc_local)

          deallocate(count)
          deallocate(compDOF)

       endif

       ! Copy the `nt` time slice data from avFile into av
       avFile_beg = lsize*(nt-1) + 1
       avFile_end = lsize*nt
       do k = 1, mct_aVect_nRAttr(avFile)
          av%rattr(k,1:lsize) = avFile%rattr(k,avFile_beg:avFile_end)
       enddo

       call t_stopf(trim(lstr)//'_readpio')

    endif

  end subroutine shr_dmodel_readstrm_fullfile

  !===============================================================================

  logical function shr_dmodel_gGridCompare(ggrid1,gsmap1,ggrid2,gsmap2,method,mpicom,eps)

    ! Returns TRUE if two domains are the the same (within tolerance).

    ! input/output variables
    type(mct_gGrid)   ,intent(in)  :: ggrid1   ! 1st ggrid
    type(mct_gsmap)   ,intent(in)  :: gsmap1   ! 1st gsmap
    type(mct_gGrid)   ,intent(in)  :: ggrid2   ! 2nd ggrid
    type(mct_gsmap)   ,intent(in)  :: gsmap2   ! 2nd gsmap
    integer           ,intent(in)  :: method   ! selects what to compare
    integer           ,intent(in)  :: mpicom   ! mpicom
    real(R8),optional ,intent(in)  :: eps      ! epsilon compare value

    ! local variables
    real(R8)        :: leps         ! local epsilon
    integer         :: n            ! counters
    integer         :: my_task,master_task
    integer         :: gsize
    integer         :: ierr
    integer         :: nlon1, nlon2, nlat1, nlat2, nmask1, nmask2  ! av field indices
    logical         :: compare               ! local compare logical
    real(R8)        :: lon1,lon2             ! longitudes to compare
    real(R8)        :: lat1,lat2             ! latitudes to compare
    real(R8)        :: msk1,msk2             ! masks to compare
    integer         :: nx,ni1,ni2            ! i grid size, i offset for 1 vs 2 and 2 vs 1
    integer         :: n1,n2,i,j             ! local indices
    type(mct_aVect) :: avG1                  ! global av
    type(mct_aVect) :: avG2                  ! global av
    integer         :: lmethod               ! local method
    logical         :: maskmethod, maskpoint ! masking on method
    character(*),parameter :: subName = '(shr_dmodel_gGridCompare) '
    character(*),parameter :: F01     = "('(shr_dmodel_gGridCompare) ',4a)"
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ierr)
    master_task = 0

    leps = 1.0e-6_R8
    if (present(eps)) leps = eps

    lmethod = mod(method,100)
    if (method > 100) then
       maskmethod=.true.
    else
       maskmethod=.false.
    endif

    call mct_aVect_gather(gGrid1%data,avG1,gsmap1,master_task,mpicom)
    call mct_aVect_gather(gGrid2%data,avG2,gsmap2,master_task,mpicom)

    if (my_task == master_task) then

       compare = .true.
       gsize = mct_aVect_lsize(avG1)
       if (gsize /= mct_aVect_lsize(avG2)) then
          compare = .false.
       endif

       if (.not. compare ) then
          !--- already failed the comparison test, check no futher ---
       else
          nlon1 = mct_aVect_indexRA(avG1,'lon')
          nlat1 = mct_aVect_indexRA(avG1,'lat')
          nlon2 = mct_aVect_indexRA(avG2,'lon')
          nlat2 = mct_aVect_indexRA(avG2,'lat')
          nmask1 = mct_aVect_indexRA(avG1,'mask')
          nmask2 = mct_aVect_indexRA(avG2,'mask')

          ! To compare, want to be able to treat longitude wraparound generally.
          ! So we need to compute i index offset and we need to compute the size of the nx dimension
          ! First adjust the lon so it's in the range [0,360), add 1440 to lon to take into
          ! accounts lons less than 1440.  if any lon is less than -1440, abort.  1440 is arbitrary
          ! Next, comute ni1 and ni2.  These are the offsets of grid1 relative to grid2 and
          ! grid2 relative to grid1.  The sum of those offsets is nx.  Use ni1 to offset grid2
          ! in comparison and compute new grid2 index from ni1 and nx.  If ni1 is zero, then
          ! there is no offset, don't need to compute ni2, and nx can be anything > 0.

          !--- compute offset of grid2 compared to pt 1 of grid 1
          lon1 = minval(avG1%rAttr(nlon1,:))
          lon2 = minval(avG2%rAttr(nlon2,:))
          if ((lon1 < -1440.0_R8) .or. (lon2 < -1440.0_R8)) then
             write(logunit,*) subname,' ERROR: lon1 lon2 lt -1440 ',lon1,lon2
             call shr_sys_abort(subname//' ERROR: lon1 lon2 lt -1440')
          endif

          lon1 = mod(avG1%rAttr(nlon1,1)+1440.0_R8,360.0_R8)
          lat1 = avG1%rAttr(nlat1,1)
          ni1 = -1
          do n = 1,gsize
             lon2 = mod(avG2%rAttr(nlon2,n)+1440.0_R8,360.0_R8)
             lat2 = avG2%rAttr(nlat2,n)
             if ((ni1 < 0) .and. abs(lon1-lon2) <= leps .and. abs(lat1-lat2) <= leps) then
                ni1 = n - 1  ! offset, compare to first gridcell in grid 1
             endif
          enddo

          if (ni1 < 0) then        ! no match for grid point 1, so fails without going further
             compare = .false.
          elseif (ni1 == 0) then   ! no offset, set nx to anything > 0
             nx = 1
          else                     ! now compute ni2
             ! compute offset of grid1 compared to pt 1 of grid 2
             lon2 = mod(avG2%rAttr(nlon2,1)+1440.0_R8,360.0_R8)
             lat2 = avG2%rAttr(nlat2,1)
             ni2 = -1
             do n = 1,gsize
                lon1 = mod(avG1%rAttr(nlon1,n)+1440.0_R8,360.0_R8)
                lat1 = avG1%rAttr(nlat1,n)
                if ((ni2 < 0) .and. abs(lon1-lon2) <= leps .and. abs(lat1-lat2) <= leps) then
                   ni2 = n - 1  ! offset, compare to first gridcell in grid 1
                endif
             enddo
             if (ni2 < 0) then
                write(logunit,*) subname,' ERROR in ni2 ',ni1,ni2
                call shr_sys_abort(subname//' ERROR in ni2')
             endif
             nx = ni1 + ni2
          endif

          if (compare) then
             do n = 1,gsize
                j = ((n-1)/nx) + 1
                i = mod(n-1,nx) + 1
                n1 = (j-1)*nx + mod(n-1,nx) + 1
                n2 = (j-1)*nx + mod(n-1+ni1,nx) + 1
                if (n1 /= n) then    ! sanity check, could be commented out
                   write(logunit,*) subname,' ERROR in n1 n2 ',n,i,j,n1,n2
                   call shr_sys_abort(subname//' ERROR in n1 n2')
                endif
                lon1 = mod(avG1%rAttr(nlon1,n1)+1440.0_R8,360.0_R8)
                lat1 = avG1%rAttr(nlat1,n1)
                lon2 = mod(avG2%rAttr(nlon2,n2)+1440.0_R8,360.0_R8)
                lat2 = avG2%rAttr(nlat2,n2)
                msk1 = avG1%rAttr(nmask1,n1)
                msk2 = avG2%rAttr(nmask2,n2)

                maskpoint = .true.
                if (maskmethod .and. (msk1 == 0 .or. msk2 == 0)) then
                   maskpoint = .false.
                endif

                if (maskpoint) then
                   if (lmethod == CompareXYabs      ) then
                      if (abs(lon1 - lon2) > leps) compare = .false.
                      if (abs(lat1 - lat2) > leps) compare = .false.
                   else if (lmethod == CompareXYrel      ) then
                      if (rdiff(lon1,lon2) > leps) compare = .false.
                      if (rdiff(lat1,lat2) > leps) compare = .false.
                   else if (lmethod == CompareMaskIdent  ) then
                      if (msk1 /= msk2)compare = .false.
                   else if (lmethod == CompareMaskZeros  ) then
                      if (msk1 == 0 .and. msk2 /= 0) compare = .false.
                      if (msk1 /= 0 .and. msk2 == 0) compare = .false.
                   else if (lmethod == CompareMaskSubset ) then
                      if (msk1 /= 0 .and. msk2 == 0) compare = .false.
                   else
                      write(logunit,F01) "ERROR: compare method not recognized, method = ",method
                      call shr_sys_abort(subName//"ERROR: compare method not recognized")
                   endif  ! lmethod
                endif  ! maskpoint
             enddo ! gsize
          endif  ! compare
       endif   ! compare
    endif   ! master_task

    if (my_task == master_task) call mct_avect_clean(avG1)
    if (my_task == master_task) call mct_avect_clean(avG2)

    call shr_mpi_bcast(compare,mpicom)
    shr_dmodel_gGridCompare = compare

  contains   ! internal subprogram

    real(R8) function rdiff(v1,v2) ! internal function
      !------------------------------------------
      real(R8),intent(in) :: v1,v2                 ! two values to compare
      real(R8),parameter  :: c0           = 0.0_R8 ! zero
      real(R8),parameter  :: large_number = 1.0e20_R8 ! infinity
      !------------------------------------------
      if (v1 == v2) then
         rdiff = c0
      elseif (v1 == c0 .and. v2 /= c0) then
         rdiff = large_number
      elseif (v2 == c0 .and. v1 /= c0) then
         rdiff = large_number
      else
         !        rdiff = abs((v2-v1)/v1)   ! old version, but rdiff(v1,v2) /= vdiff(v2,v1)
         rdiff = abs(2.0_R8*(v2-v1)/(v1+v2))
      endif
      !------------------------------------------
    end function rdiff

  end function shr_dmodel_gGridCompare

  !===============================================================================

  subroutine shr_dmodel_mapSet(smatp,&
       ggridS,gsmapS,nxgS,nygS,&
       ggridD,gsmapD,nxgD,nygD, &
       name,type,algo,mask,vect,&
       compid,mpicom,strategy)

    !----- arguments -----
    type(mct_sMatP)  , intent(inout)       :: smatp
    type(mct_gGrid)  , intent(in)          :: ggridS
    type(mct_gsmap)  , intent(in)          :: gsmapS
    integer          , intent(in)          :: nxgS
    integer          , intent(in)          :: nygS
    type(mct_gGrid)  , intent(in)          :: ggridD
    type(mct_gsmap)  , intent(in)          :: gsmapD
    integer          , intent(in)          :: nxgD
    integer          , intent(in)          :: nygD
    character(len=*) , intent(in)          :: name
    character(len=*) , intent(in)          :: type
    character(len=*) , intent(in)          :: algo
    character(len=*) , intent(in)          :: mask
    character(len=*) , intent(in)          :: vect
    integer          , intent(in)          :: compid
    integer          , intent(in)          :: mpicom
    character(len=*) , intent(in),optional :: strategy

    !----- local -----

    integer :: n,i,j
    integer :: lsizeS,gsizeS,lsizeD,gsizeD
    integer :: nlon,nlat,nmsk
    integer :: my_task,master_task,ierr

    real(R8)   , pointer :: Xsrc(:,:)
    real(R8)   , pointer :: Ysrc(:,:)
    integer, pointer :: Msrc(:,:)
    real(R8)   , pointer :: Xdst(:,:)
    real(R8)   , pointer :: Ydst(:,:)
    integer, pointer :: Mdst(:,:)
    type(shr_map_mapType) :: shrmap
    type(mct_aVect) :: AVl
    type(mct_aVect) :: AVg

    character(len=32) :: lstrategy
    integer :: nsrc,ndst,nwts
    integer, pointer :: isrc(:)
    integer, pointer :: idst(:)
    real(R8)   , pointer :: wgts(:)
    type(mct_sMat) :: sMat0

    character(*), parameter :: subname = '(shr_dmodel_mapSet) '
    character(*), parameter :: F00   = "('(shr_dmodel_mapSet) ',8a)"
    character(*), parameter :: F01   = "('(shr_dmodel_mapSet) ',a,5i8)"

    !-------------------------------------------------------------------------------
    ! PURPOSE:  Initialize sMatP from mct gGrid
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ierr)
    master_task = 0

    !--- get sizes and allocate for SRC ---

    lsizeS = mct_aVect_lsize(ggridS%data)
    call mct_avect_init(AVl,rList='lon:lat:mask',lsize=lsizeS)
    call mct_avect_copy(ggridS%data,AVl,rList='lon:lat:mask')
    call mct_avect_gather(AVl,AVg,gsmapS,master_task,mpicom)

    if (my_task == master_task) then
       gsizeS = mct_aVect_lsize(AVg)
       if (gsizeS /= nxgS*nygS) then
          write(logunit,F01) ' ERROR: gsizeS ',gsizeS,nxgS,nygS
          call shr_sys_abort(subname//' ERROR gsizeS')
       endif
       allocate(Xsrc(nxgS,nygS),Ysrc(nxgS,nygS),Msrc(nxgS,nygS))

       nlon = mct_avect_indexRA(AVg,'lon')
       nlat = mct_avect_indexRA(AVg,'lat')
       nmsk = mct_avect_indexRA(AVg,'mask')

       n = 0
       Msrc = 1
       do j = 1,nygS
          do i = 1,nxgS
             n = n + 1
             Xsrc(i,j) = AVg%rAttr(nlon,n)
             Ysrc(i,j) = AVg%rAttr(nlat,n)
             if (abs(AVg%rAttr(nmsk,n)) < 0.5_R8) Msrc(i,j) = 0
          enddo
       enddo
    endif

    if (my_task == master_task) call mct_aVect_clean(AVg)
    call mct_aVect_clean(AVl)

    !--- get sizes and allocate for DST ---

    lsizeD = mct_aVect_lsize(ggridD%data)
    call mct_avect_init(AVl,rList='lon:lat:mask',lsize=lsizeD)
    call mct_avect_copy(ggridD%data,AVl,rList='lon:lat:mask')
    call mct_avect_gather(AVl,AVg,gsmapD,master_task,mpicom)

    if (my_task == master_task) then
       gsizeD = mct_aVect_lsize(AVg)
       if (gsizeD /= nxgD*nygD) then
          write(logunit,F01) ' ERROR: gsizeD ',gsizeD,nxgD,nygD
          call shr_sys_abort(subname//' ERROR gsizeD')
       endif
       allocate(Xdst(nxgD,nygD),Ydst(nxgD,nygD),Mdst(nxgD,nygD))

       nlon = mct_avect_indexRA(AVg,'lon')
       nlat = mct_avect_indexRA(AVg,'lat')
       nmsk = mct_avect_indexRA(AVg,'mask')

       n = 0
       Mdst = 1
       do j = 1,nygD
          do i = 1,nxgD
             n = n + 1
             Xdst(i,j) = AVg%rAttr(nlon,n)
             Ydst(i,j) = AVg%rAttr(nlat,n)
             if (abs(AVg%rAttr(nmsk,n)) < 0.5_R8) Mdst(i,j) = 0
          enddo
       enddo
    endif

    if (my_task == master_task) call mct_aVect_clean(AVg)
    call mct_aVect_clean(AVl)

    !--- set map ---

    if (my_task == master_task) then
       call shr_map_mapSet(shrmap,Xsrc,Ysrc,Msrc,Xdst,Ydst,Mdst, &
            trim(name),trim(type),trim(algo),trim(mask),trim(vect))
       deallocate(Xsrc,Ysrc,Msrc)
       deallocate(Xdst,Ydst,Mdst)
    endif

    !--- convert map to sMatP ---

    lstrategy = 'Xonly'
    if (present(strategy)) then
       lstrategy = trim(strategy)
    endif

    if (my_task == master_task) then
       call shr_map_get(shrmap,shr_map_fs_nsrc,nsrc)
       call shr_map_get(shrmap,shr_map_fs_ndst,ndst)
       call shr_map_get(shrmap,shr_map_fs_nwts,nwts)
       allocate(isrc(nwts),idst(nwts),wgts(nwts))
       call shr_map_get(shrmap,isrc,idst,wgts)
       call shr_map_clean(shrmap)

       call mct_sMat_init(sMat0,ndst,nsrc,nwts)

       call mct_sMat_ImpGColI (sMat0,isrc,nwts)
       call mct_sMat_ImpGRowI (sMat0,idst,nwts)
       call mct_sMat_ImpMatrix(sMat0,wgts,nwts)
       deallocate(isrc,idst,wgts)
    endif

    call mct_sMatP_Init(sMatP,sMat0,gsmapS,gsmapD,lstrategy,master_task,mpicom,compid)

    if (my_task == master_task) then
       call mct_sMat_clean(sMat0)
    endif

  end subroutine shr_dmodel_mapSet

  !===============================================================================

end module dshr_dmodel_mod
