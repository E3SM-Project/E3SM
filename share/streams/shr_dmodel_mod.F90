module shr_dmodel_mod

  ! !USES:

  use shr_sys_mod
  use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, &
       R4=>SHR_KIND_R4, &
       CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, &
       CX=>SHR_KIND_CX, CXX=>SHR_KIND_CXX
  use shr_log_mod, only: loglev  => shr_log_Level
  use shr_log_mod, only: logunit => shr_log_Unit
  use shr_mpi_mod, only: shr_mpi_bcast
  use mct_mod
  use perf_mod
  use pio
  use m_SparseMatrixPlus, only : exportStrategyToChar

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: shr_dmodel_gsmapCreate
  public :: shr_dmodel_readLBUB
  public :: shr_dmodel_readgrid
  public :: shr_dmodel_gGridCompare
  public :: shr_dmodel_mapSet
  public :: shr_dmodel_translateAV
  public :: shr_dmodel_translateAV_list
  public :: shr_dmodel_translate_list
  public :: shr_dmodel_rearrGGrid

  interface shr_dmodel_gsmapCreate; module procedure &
       shr_dmodel_gsmapCreate_gsize, &
       shr_dmodel_gsmapCreate_nxnynz
  end interface shr_dmodel_gsmapCreate

  interface shr_dmodel_mapSet; module procedure &
       !shr_dmodel_mapSet_dest, &
       shr_dmodel_mapSet_global
  end interface shr_dmodel_mapSet

  integer(IN),parameter,public :: shr_dmodel_gGridCompareXYabs      = 1 ! X,Y  relative error
  integer(IN),parameter,public :: shr_dmodel_gGridCompareXYrel      = 2 ! X,Y  absolute error
  integer(IN),parameter,public :: shr_dmodel_gGridCompareAreaAbs    = 3 ! area relative error
  integer(IN),parameter,public :: shr_dmodel_gGridCompareAreaRel    = 4 ! area absolute error
  integer(IN),parameter,public :: shr_dmodel_gGridCompareMaskIdent  = 5 ! masks are identical
  integer(IN),parameter,public :: shr_dmodel_gGridCompareMaskZeros  = 6 ! masks have same zeros
  integer(IN),parameter,public :: shr_dmodel_gGridCompareMaskSubset = 7 ! mask is subset of other

  ! masked methods

  integer(IN),parameter,public :: shr_dmodel_gGridCompareXYabsMask      = 101 ! X,Y  relative error
  integer(IN),parameter,public :: shr_dmodel_gGridCompareXYrelMask      = 102 ! X,Y  absolute error
  integer(IN),parameter,public :: shr_dmodel_gGridCompareAreaAbsMask    = 103 ! area relative error
  integer(IN),parameter,public :: shr_dmodel_gGridCompareAreaRelMask    = 104 ! area absolute error
  integer(IN),parameter,public :: shr_dmodel_gGridCompareMaskIdentMask  = 105 ! masks are identical
  integer(IN),parameter,public :: shr_dmodel_gGridCompareMaskZerosMask  = 106 ! masks have same zeros
  integer(IN),parameter,public :: shr_dmodel_gGridCompareMaskSubsetMask = 107 ! mask is subset of other

  integer(IN),parameter,public :: iotype_std_netcdf = -99 ! non pio option

!===============================================================================
CONTAINS
!===============================================================================

  subroutine shr_dmodel_gsmapCreate_gsize(gsmap,gsize,compid,mpicom,decomp)

    type(mct_gsMap), intent(inout) :: gsmap
    integer(IN)    , intent(in)    :: gsize
    integer(IN)    , intent(in)    :: compid
    integer(IN)    , intent(in)    :: mpicom
    character(len=*),intent(in)    :: decomp

    ! local
    integer(IN) :: n,npes,ierr
    integer(IN), pointer :: start(:)     ! for gsmap initialization
    integer(IN), pointer :: length(:)    ! for gsmap initialization
    integer(IN), pointer :: pe_loc(:)    ! for gsmap initialization
    character(*), parameter :: subname = '(shr_dmodel_gsmapCreate_gsize) '
    character(*), parameter :: F00   = "('(shr_dmodel_gsmapCreate_gsize) ',8a)"
    character(*), parameter :: F01   = "('(shr_dmodel_gsmapCreate_gsize) ',a,5i8)"
    ! ---------------------------------------------

    if (gsize > 0) then
       call mpi_comm_size(mpicom,npes,ierr)
       allocate(start(npes),length(npes),pe_loc(npes))

       start = 0
       length = 0
       do n = 1,npes
          if (trim(decomp) == '1d') then
             length(n)  = gsize/npes
             if (n <= mod(gsize,npes)) length(n) = length(n) + 1
          elseif (trim(decomp) == 'root') then
             length = 0
             length(1) = gsize
          else
             write(logunit,F00) ' ERROR: decomp not allowed, ',trim(decomp)
             call shr_sys_abort(subname//' ERROR decomp')
          endif
          if (n == 1) then
             start(n) = 1
          else
             start(n) = start(n-1) + length(n-1)
          endif
          pe_loc(n) = n-1
       enddo
       call mct_gsMap_init( gsMap, COMPID, npes, gsize, start, length, pe_loc)
       deallocate(start,length,pe_loc)
    endif

  end subroutine shr_dmodel_gsmapCreate_gsize

  !===============================================================================

  subroutine shr_dmodel_gsmapCreate_nxnynz(gsmap,nxg,nyg,nzg,compid,mpicom,decomp)

    type(mct_gsMap), intent(inout) :: gsmap
    integer(IN)    , intent(in)    :: nxg,nyg,nzg
    integer(IN)    , intent(in)    :: compid
    integer(IN)    , intent(in)    :: mpicom
    character(len=*),intent(in)    :: decomp

    ! local

    integer(IN) :: n,nz,nb,npes,ierr,gsize,dsize,ngseg,lnzg
    integer(IN), pointer :: start(:)     ! for gsmap initialization
    integer(IN), pointer :: length(:)    ! for gsmap initialization
    integer(IN), pointer :: pe_loc(:)    ! for gsmap initialization
    character(*), parameter :: subname = '(shr_dmodel_gsmapCreate_nxnynz) '
    character(*), parameter :: F00   = "('(shr_dmodel_gsmapCreate_nxnynz) ',8a)"
    character(*), parameter :: F01   = "('(shr_dmodel_gsmapCreate_nxnynz) ',a,5i8)"

    ! ---------------------------------------------

    gsize = abs(nxg*nyg*nzg)
    dsize = nxg*nyg
    lnzg = 1
    if (nzg > 1) lnzg = nzg  ! check for 3d

    if (gsize > 0) then
       call mpi_comm_size(mpicom,npes,ierr)

       !--- 1d decomp of 2d grid plus 3rd dim if exists ---
       if (trim(decomp) == '2d1d') then
          ngseg = npes*lnzg
          allocate(start(ngseg),length(ngseg),pe_loc(ngseg))
          start = 0
          length = 0
          pe_loc = 0
          do n = 1,npes
             length(n)  = dsize/npes
             if (n <= mod(dsize,npes)) length(n) = length(n) + 1
             if (n == 1) then
                start(n) = 1
             else
                start(n) = start(n-1) + length(n-1)
             endif
             pe_loc(n) = n-1
             do nz = 2,lnzg
                nb = (nz-1)*npes + n
                start(nb)  = start(n) + (nz-1)*dsize
                length(nb) = length(n)
                pe_loc(nb) = pe_loc(n)
             enddo
          enddo

          !--- all data on root ---
       elseif (trim(decomp) == 'root') then
          ngseg = 1
          allocate(start(ngseg),length(ngseg),pe_loc(ngseg))
          start(1) = 1
          length(1) = gsize
          pe_loc(1) = 0

       else
          write(logunit,F00) ' ERROR: decomp not allowed, ',trim(decomp)
          call shr_sys_abort(subname//' ERROR decomp')
       endif

       call mct_gsMap_init( gsMap, COMPID, ngseg, gsize, start, length, pe_loc)
       deallocate(start,length,pe_loc)
    endif

  end subroutine shr_dmodel_gsmapCreate_nxnynz

  !===============================================================================

  subroutine shr_dmodel_readgrid( gGrid, gsMap, nxgo, nygo, nzgo, filename, compid, mpicom, &
       decomp, lonname, latname, hgtname, maskname, areaname, fracname, readfrac, &
       scmmode, scm_multcols, scmlon, scmlat, scm_nx, scm_ny)

    use shr_file_mod  , only : shr_file_noprefix, shr_file_queryprefix, shr_file_get
    use shr_string_mod, only : shr_string_lastindex
    use shr_ncread_mod, only : shr_ncread_domain, shr_ncread_vardimsizes
    use shr_ncread_mod, only : shr_ncread_varexists, shr_ncread_vardimnum,  shr_ncread_field4dG

    !----- arguments -----
    type(mct_gGrid)  ,          intent(inout) :: gGrid
    type(mct_gsMap)  ,          intent(inout) :: gsMap
    integer(IN)      ,          intent(out)   :: nxgo
    integer(IN)      ,          intent(out)   :: nygo
    integer(IN)      ,          intent(out)   :: nzgo
    character(len=*) ,          intent(in)    :: filename
    integer(IN)      ,          intent(in)    :: compid
    integer(IN)      ,          intent(in)    :: mpicom
    character(len=*) ,optional, intent(in)    :: decomp   ! decomp strategy for gsmap
    character(len=*) ,optional, intent(in)    :: lonname  ! name of  lon variable in file
    character(len=*) ,optional, intent(in)    :: latname  ! name of  lat variable in file
    character(len=*) ,optional, intent(in)    :: hgtname  ! name of  hgt variable in file
    character(len=*) ,optional, intent(in)    :: maskname ! name of mask variable in file
    character(len=*) ,optional, intent(in)    :: areaname ! name of area variable in file
    character(len=*) ,optional, intent(in)    :: fracname ! name of frac variable in file
    logical          ,optional, intent(in)    :: readfrac ! T <=> also read frac  in file
    logical          ,optional, intent(in)    :: scmmode  ! single column mode
    logical          ,optional, intent(in)    :: scm_multcols ! SCM mode for multiple columns
    real(R8)         ,optional, intent(in)    :: scmlon   ! single column lon
    real(R8)         ,optional, intent(in)    :: scmlat   ! single column lat
    integer(IN)      ,optional, intent(in)    :: scm_nx   ! points in x direction for SCM functionality
    integer(IN)      ,optional, intent(in)    :: scm_ny   ! points in y direction for SCM functionality

    !----- local -----
    integer(IN)   :: n,k,j,i     ! indices
    integer(IN)   :: lsize       ! lsize
    integer(IN)   :: gsize       ! gsize
    integer(IN)   :: my_task, master_task
    integer(IN)   :: ierr        ! error code
    logical       :: fileexists  !
    integer(IN)   :: rCode       ! return code
    character(CL) :: remoteFn    ! input file name (possibly at an archival location)
    character(CL) :: localFn     ! file name to be opened (possibly a local copy)
    character(CS) :: prefix      ! file prefix
    character(CS) :: ldecomp     ! decomp strategy
    character(CS) :: llonname    ! name of  lon variable
    character(CS) :: llatname    ! name of  lat variable
    character(CS) :: lhgtname    ! name of  hgt variable
    character(CS) :: lmaskname   ! name of mask variable
    character(CS) :: lareaname   ! name of area variable
    character(CS) :: lfracname   ! name of area variable
    logical       :: lreadfrac   ! read fraction
    logical       :: maskexists  ! is mask on dataset
    integer(IN)   :: nxg,nyg,nzg ! size of input fields
    integer(IN)   :: ndims       ! number of dims
    integer(IN)   :: nlon,nlat,narea,nmask,nfrac,nhgt
    logical       :: lscmmode    ! local scm mode
    logical       :: lscm_multcols! local scm multiple column mode
    logical       :: lscmgrid    ! have scm grid information
    real(R8)      :: dist,mind   ! scmmode point search
    integer(IN)   :: ni,nj       ! scmmode point search
    real(R8)      :: lscmlon     ! local copy of scmlon
    integer(IN)   :: i_scm, j_scm

    real   (R8),allocatable ::  lon(:,:) ! temp array for domain lon  info
    real   (R8),allocatable ::  lat(:,:) ! temp array for domain lat  info
    integer(IN),allocatable :: mask(:,:) ! temp array for domain mask info
    real   (R8),allocatable :: area(:,:) ! temp array for domain area info
    real   (R8),allocatable :: frac(:,:) ! temp array for domain frac info
    real   (R8),allocatable ::  hgt(:)   ! temp array for domain height info
    real   (R8),allocatable ::  a4d(:,:,:,:) ! temp array for reading generic stuff

    integer(IN), pointer :: idata(:)   ! temporary
    type(mct_ggrid)      :: gGridRoot       ! global mct ggrid

    character(*), parameter :: subname = '(shr_dmodel_readgrid) '
    character(*), parameter :: F00   = "('(shr_dmodel_readgrid) ',8a)"
    character(*), parameter :: F01   = "('(shr_dmodel_readgrid) ',a,5i8)"

    !-------------------------------------------------------------------------------
    ! PURPOSE:  Read MCT ggrid and set gsmap if not input
    !----------------------------------------------------------------------------
    ! Notes:
    ! o as per shr_file_get(), the file name format is expected to be
    !   remoteFn = [location:][directory path]localFn
    !   eg. "foobar.nc"  "/home/user/foobar.nc"  "mss:/USER/fobar.nc"
    ! o assumes a very specific netCDF domain file format wrt var names, etc.
    !
    ! TO DO: have the calling routine select/input the domain's file name
    !----------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ierr)
    master_task = 0

    lscmmode = .false.
    lscm_multcols = .false.
    lscmgrid = .false.
    if (present(scmmode)) then
       lscmmode = scmmode
       lscm_multcols = scm_multcols
       if (lscm_multcols .and. .not. lscmmode) then
         write(logunit,*) subname, ' ERROR: SCM mode for multiple columns must be run with SCM functionality'
         call shr_sys_abort(subname//' ERROR: SCM multiple column mode must be in SCM mode')
       endif
       if (lscmmode) then
          if (.not.present(scmlon) .or. .not.present(scmlat)) then
             write(logunit,*) subname,' ERROR: scmmode must supply scmlon and scmlat'
             call shr_sys_abort(subname//' ERROR: scmmode1 lon lat')
          endif
          if (my_task > 0 .and. .not. lscm_multcols) then
             write(logunit,*) subname,' ERROR: scmmode must be run on one pe'
             call shr_sys_abort(subname//' ERROR: scmmode2 tasks')
          endif
       endif
       if (lscm_multcols) then
         if (present(scm_nx) .and. present(scm_ny)) then
           lscmgrid = .true.
         endif
       endif
    endif

    if (my_task == master_task) then
       if ( shr_file_queryPrefix(fileName,prefix=prefix) /= shr_file_noPrefix ) then
          n        = max(len_trim(prefix),shr_string_lastIndex(fileName,"/"))
          remoteFn = fileName
          localFn  = fileName(n+1: len_trim(fileName) )
          call shr_file_get(rCode,localFn,remoteFn)
       else
          remoteFn = "undefined" ! this isn't needed
          localFn  = fileName    ! file to open
       end if
       inquire(file=trim(localFn),exist=fileExists)
       if (.not. fileExists) then
          write(logunit,F00) "ERROR: file does not exist: ", trim(localFn)
          call shr_sys_abort(subName//"ERROR: file does not exist: "//trim(localFn))
       end if
    endif

    lreadfrac = .false.
    ldecomp   = "2d1d"
    llonname  = "xc"  ! default values / standard data model domain file format
    llatname  = "yc"
    lhgtname  = "hgt"
    lmaskname = "mask"
    lareaname = "area"
    lfracname = "frac"
    if (present(  decomp))   ldecomp =   decomp
    if (present(readfrac)) lreadfrac = readfrac
    if (present( lonname))  llonname =  lonname
    if (present( latname))  llatname =  latname
    if (present( hgtname))  lhgtname =  hgtname
    if (present(maskname)) lmaskname = maskname
    if (present(areaname)) lareaname = areaname
    if (present(fracname)) lfracname = fracname

    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (ocean), 0 (non-ocean)

    if (my_task == master_task) then
       if (shr_ncread_varexists(localFn,lmaskname)) then
          maskexists = .true.
          call shr_ncread_varDimSizes(localFn,lmaskname,nxg,nyg)
       else
          maskexists = .false.
          call shr_ncread_varDimNum(localFn,llonName,ndims)
          if (ndims == 1) then
             call shr_ncread_varDimSizes(localFn,llonName,nxg)
             call shr_ncread_varDimSizes(localFn,llatName,nyg)
          else
             call shr_ncread_varDimSizes(localFn,llonName,nxg,nyg)
          endif
       endif
       if (shr_ncread_varexists(localFn,lhgtName)) then
          call shr_ncread_varDimSizes(localFn,lhgtname,nzg)
       else
          nzg = -1
       endif
    endif
    call shr_mpi_bcast(nxg,mpicom)
    call shr_mpi_bcast(nyg,mpicom)
    call shr_mpi_bcast(nzg,mpicom)
    if (lscmmode .and. .not. lscm_multcols) then
       ! Standard SCM mode
       nxgo = 1
       nygo = 1
       nzgo = -1
       gsize = 1
    else
       if (lscmgrid) then
         ! Run on user defined grid with SCM functionality
         nxgo = scm_nx
         nygo = scm_ny
         nzgo = -1
         gsize = abs(scm_nx*scm_nx*nzgo)
       else
         nxgo = nxg
         nygo = nyg
         nzgo = nzg
         gsize = abs(nxg*nyg*nzg)
       endif
       if (gsize < 1) return
    endif

    ! Create gsmap if input gsmap is not given

    if (trim(ldecomp) == 'use_input_gsmap') then
       if (my_task == master_task) then
          write(logunit,*)trim(subname) // ': Using input gsmap'
       end if
    else
       call shr_dmodel_gsmapCreate(gsMap, nxgo, nygo, nzgo, compid, mpicom, trim(ldecomp))
    end if

    ! Create gGrid using either created or input gsmap

    lsize = mct_gsMap_lsize(gsMap, mpicom)
    call mct_gGrid_init(GGrid=Ggrid, CoordChars='lat:lon:hgt', OtherChars='area:aream:mask:frac', lsize=lsize )

    ! Determine global gridpoint number attribute, GlobGridNum, automatically in ggrid

    call mct_gsMap_orderedPoints(gsMap, my_task, idata)
    call mct_gGrid_importIAttr(gGrid,'GlobGridNum',idata,lsize)
    deallocate(idata)

    ! Initialize attribute vector with special value

    gGrid%data%rAttr = -9999.0_R8

    ! Load file data into domG then scatter to ggrid

    if (my_task == master_task) then

       allocate(lon(nxg,nyg))
       allocate(lat(nxg,nyg))
       allocate(area(nxg,nyg))
       allocate(mask(nxg,nyg))
       allocate(frac(nxg,nyg))
       allocate(hgt(abs(nzg)))

       if (.not.maskexists) then
          call shr_ncread_domain(localFn,llonName,lon,llatName,lat)
          mask = 1
          frac = 1.0_R8
          area = 1.0e36_R8
       else
          if (lreadfrac) then
             call shr_ncread_domain(localFn,llonName,lon,llatName,lat, &
                  lmaskName,mask,lareaName,area,lfracName,frac)
          else ! assume frac = 1.0
             call shr_ncread_domain(localFn,llonName,lon,llatName,lat, &
                  lmaskName,mask,lareaName,area)
             where (mask == 0)
                frac = 0.0_R8
             elsewhere
                frac = 1.0_R8
             end where
          endif
       endif

       if (nzg > 1) then
          allocate(a4d(nzg,1,1,1))
          call shr_ncread_field4dG(localFn,hgtName,rfld=a4d)
          hgt(:) = a4d(:,1,1,1)
          deallocate(a4d)
       else
          hgt = 1
       endif

       call mct_gGrid_init(gGridRoot,gGrid,gsize)

       ! initialize gGridRoot to avoid errors when using strict compiler checks
       gGridRoot%data%rAttr = -9999.0_R8

       nlon  = mct_aVect_indexRA(gGridRoot%data,'lon')
       nlat  = mct_aVect_indexRA(gGridRoot%data,'lat')
       narea = mct_aVect_indexRA(gGridRoot%data,'area')
       nmask = mct_aVect_indexRA(gGridRoot%data,'mask')
       nfrac = mct_aVect_indexRA(gGridRoot%data,'frac')
       nhgt  = mct_aVect_indexRA(gGridRoot%data,'hgt')

       if (lscmmode) then
          !--- assumes regular 2d grid for compatability with shr_scam_getCloseLatLon ---
          !--- want lon values between 0 and 360, assume 1440 is enough ---
          lscmlon = mod(scmlon+1440.0_r8,360.0_r8)
          lon     = mod(lon   +1440.0_r8,360.0_r8)

          !--- start with wraparound ---

          !--- determine whether dealing with 2D input files (typical of Eulerian
          !--- dynamical core) or 1D files (typical of Spectral Element)

          if (nyg .ne. 1) then
            ni = 1
            mind = abs(lscmlon - (lon(1,1)+360.0_r8))
            do i=1,nxg
              dist = abs(lscmlon - lon(i,1))
              if (dist < mind) then
                mind = dist
                ni = i
              endif
            enddo

            nj = -1
            mind = 1.0e20
            do j=1,nyg
              dist = abs(scmlat - lat(1,j))
              if (dist < mind) then
                mind = dist
                nj = j
              endif
            enddo

            j = nj
            j_scm = nj

          else ! lat and lon are on 1D arrays

            !--- to deal with spectral element grids
            mind = 1.0e20
            do i=1,nxg
              dist=abs(lscmlon - lon(i,1)) + abs(scmlat - lat(i,1))
              if (dist < mind) then
                mind = dist
                ni = i
              endif
            enddo

            j = 1
            j_scm = 1

          endif

          i_scm = ni

          if (lscm_multcols) then

            ! If using SCM functionality across multiple columns,
            !   then we want the surface to be covered homogeneously, with
            !   the same lat and lon as close to the lat/lon in IOP forcing
            !   file as possible
            i_scm = ni

            n=0
            do k=1,abs(nzg)
              do j=1,nygo
                do i=1,nxgo
                  n=n+1
                  gGridRoot%data%rAttr(nlat ,n) = lat(i_scm,j_scm)
                  gGridRoot%data%rAttr(nlon ,n) = lon(i_scm,j_scm)
                  gGridRoot%data%rAttr(narea,n) = area(i_scm,j_scm)
                  gGridRoot%data%rAttr(nmask,n) = real(mask(i_scm,j_scm),R8)
                  gGridRoot%data%rAttr(nfrac,n) = frac(i_scm,j_scm)
                  gGridRoot%data%rAttr(nhgt ,n) = hgt(k)
                enddo
              enddo
            enddo

            else

            i = ni
            n = 1

            gGridRoot%data%rAttr(nlat ,n) = lat(i,j)
            gGridRoot%data%rAttr(nlon ,n) = lon(i,j)
            gGridRoot%data%rAttr(narea,n) = area(i,j)
            gGridRoot%data%rAttr(nmask,n) = real(mask(i,j),R8)
            gGridRoot%data%rAttr(nfrac,n) = frac(i,j)
            gGridRoot%data%rAttr(nhgt, n) = 1

          endif

       else
          n=0
          do k=1,abs(nzg)
             do j=1,nyg
                do i=1,nxg
                   n=n+1
                   gGridRoot%data%rAttr(nlat ,n) = lat(i,j)
                   gGridRoot%data%rAttr(nlon ,n) = lon(i,j)
                   gGridRoot%data%rAttr(narea,n) = area(i,j)
                   gGridRoot%data%rAttr(nmask,n) = real(mask(i,j),R8)
                   gGridRoot%data%rAttr(nfrac,n) = frac(i,j)
                   gGridRoot%data%rAttr(nhgt ,n) = hgt(k)
                enddo
             enddo
          enddo
       endif
    endif

    if (my_task == master_task) then
       deallocate(lon)
       deallocate(lat)
       deallocate(area)
       deallocate(mask)
       deallocate(frac)
       deallocate(hgt)
    endif

    call mct_gGrid_scatter(gGridRoot, gGrid, gsMap, master_task, mpicom)
    if (my_task == master_task) call mct_gGrid_clean(gGridRoot)

  end subroutine shr_dmodel_readgrid

  !===============================================================================

  subroutine shr_dmodel_readLBUB(stream,pio_subsystem,pio_iotype,pio_iodesc_r8, pio_iodesc_r4, pio_iodesc_int, &
       mDate,mSec,mpicom,gsMap, &
       avLB,mDateLB,mSecLB,avUB,mDateUB,mSecUB,avFile,readMode, &
       newData,rmOldFile,istr)

    use shr_file_mod, only : shr_file_noprefix, shr_file_queryprefix, shr_file_get
    use shr_const_mod, only : shr_const_cDay
    use shr_stream_mod

    !----- arguments -----
    type(shr_stream_streamType),intent(inout) :: stream
    type(iosystem_desc_t)      ,intent(inout), target :: pio_subsystem
    integer(IN)                ,intent(in)    :: pio_iotype
    type(io_desc_t)            ,intent(inout) :: pio_iodesc_r8, pio_iodesc_r4, pio_iodesc_int
    integer(IN)                ,intent(in)    :: mDate  ,mSec
    integer(IN)                ,intent(in)    :: mpicom
    type(mct_gsMap)            ,intent(in)    :: gsMap
    type(mct_aVect)            ,intent(inout) :: avLB
    integer(IN)                ,intent(inout) :: mDateLB,mSecLB
    type(mct_aVect)            ,intent(inout) :: avUB
    integer(IN)                ,intent(inout) :: mDateUB,mSecUB
    type(mct_aVect)            ,intent(inout) :: avFile
    character(len=*)           ,intent(in)    :: readMode
    logical                    ,intent(out)   :: newData
    logical,optional           ,intent(in)    :: rmOldFile
    character(len=*),optional  ,intent(in)    :: istr

    !----- local -----
    integer(IN)   :: my_task, master_task
    integer(IN)   :: ierr       ! error code
    integer(IN)   :: rCode      ! return code
    logical       :: localCopy,fileexists
    integer(IN)   :: ivals(6)   ! bcast buffer

    integer(IN)   :: oDateLB,oSecLB,dDateLB,oDateUB,oSecUB,dDateUB
    real(R8)      :: rDateM,rDateLB,rDateUB  ! model,LB,UB dates with fractional days
    integer(IN)   ::  n_lb, n_ub
    character(CL) :: fn_lb,fn_ub,fn_next,fn_prev
    character(CL) :: path
    character(len=32) :: lstr

    real(R8)      :: spd

    character(*), parameter :: subname = '(shr_dmodel_readLBUB) '
    character(*), parameter :: F00   = "('(shr_dmodel_readLBUB) ',8a)"
    character(*), parameter :: F01   = "('(shr_dmodel_readLBUB) ',a,5i8)"

    !-------------------------------------------------------------------------------
    ! PURPOSE:  Read LB and UB stream data
    !----------------------------------------------------------------------------

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

    if (rDateM < rDateLB .or. rDateM >= rDateUB) then
       call t_startf(trim(lstr)//'_fbound')
       if (my_task == master_task) then
          !       call shr_stream_findBounds(stream,mDate,mSec,                 &
          !                                  mDateLB,dDateLB,mSecLB,n_lb,fn_lb, &
          !                                  mDateUB,dDateUB,mSecUB,n_ub,fn_ub  )
          call shr_stream_findBounds(stream,mDate,mSec,                 &
               ivals(1),dDateLB,ivals(2),ivals(5),fn_lb, &
               ivals(3),dDateUB,ivals(4),ivals(6),fn_ub  )
          call shr_stream_getFilePath(stream,path)
          localCopy = (shr_file_queryPrefix(path) /= shr_file_noPrefix )
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
             call shr_dmodel_readstrm(stream, pio_subsystem, pio_iotype, pio_iodesc_r8, pio_iodesc_r4, &
                  pio_iodesc_int, gsMap, avLB, mpicom, &
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
          call shr_dmodel_readstrm(stream, pio_subsystem, pio_iotype, pio_iodesc_r8, pio_iodesc_r4, &
               pio_iodesc_int, gsMap, avUB, mpicom, &
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
       localCopy = (shr_file_queryPrefix(path) /= shr_file_noPrefix )
       if (localCopy) then
          call shr_stream_getPrevFileName(stream,fn_lb,fn_prev,path)
          call shr_stream_getNextFileName(stream,fn_ub,fn_next,path)
          inquire(file=trim(fn_next),exist=fileExists)
          if ( trim(fn_next) == "unknown" .or. fileExists) then
             ! do nothing
          else
             call shr_file_get(rCode,fn_next,trim(path)//fn_next,async=.true.)
             write(logunit,F00) "get next file: ",trim(fn_next)
             call shr_sys_flush(logunit)
          end if

          !---  remove the old file? (only if acquiring local copies) ---
          if ( rmOldFile .and. &
               fn_prev/=fn_lb .and. fn_prev/=fn_ub .and. fn_prev/=fn_next ) then
             !--- previous file is not in use and is not next in list ---
             inquire(file=trim(fn_prev),exist=fileExists)
             if ( fileExists ) then
                call shr_sys_system(" rm "//trim(fn_prev),rCode)
                write(logunit,F00) "rm  prev file: ",trim(fn_prev)
                call shr_sys_flush(logunit)
             end if
          end if
       endif
    endif
    call t_stopf(trim(lstr)//'_filemgt')

  end subroutine shr_dmodel_readLBUB

  !===============================================================================

  subroutine shr_dmodel_readstrm(stream, pio_subsystem, pio_iotype, pio_iodesc_r8, pio_iodesc_r4, pio_iodesc_int, &
       gsMap, av, mpicom, path, fn, nt, istr, boundstr)

    use shr_file_mod, only : shr_file_noprefix, shr_file_queryprefix, shr_file_get
    use shr_stream_mod
    use shr_ncread_mod, only: shr_ncread_open, shr_ncread_close, shr_ncread_varDimSizes, shr_ncread_tField

    !----- arguments -----
    type(shr_stream_streamType),intent(inout) :: stream
    type(iosystem_desc_t),intent(inout), target  :: pio_subsystem
    integer(IN)     ,intent(in)    :: pio_iotype
    type(io_desc_t)      ,intent(inout)  :: pio_iodesc_r8, pio_iodesc_r4, pio_iodesc_int
    type(mct_gsMap) ,intent(in)    :: gsMap
    type(mct_aVect) ,intent(inout) :: av
    integer(IN)     ,intent(in)    :: mpicom
    character(len=*),intent(in)    :: path
    character(len=*),intent(in)    :: fn
    integer(IN)     ,intent(in)    :: nt
    character(len=*),optional  ,intent(in)    :: istr
    character(len=*),optional  ,intent(in)    :: boundstr

    !----- local -----
    integer(IN) :: my_task
    integer(IN) :: master_task
    integer(IN) :: ierr
    logical     :: localCopy,fileexists
    type(mct_avect) :: avG
    integer(IN) :: gsize,nx,ny,nz
    integer(IN) :: k
    integer(IN) :: fid
    integer(IN) :: rCode      ! return code
    real(R8),allocatable :: data2d(:,:)
    real(R8),allocatable :: data3d(:,:,:)
    real(r4), allocatable :: rdata(:)
    integer, allocatable :: idata(:)
    logical     :: d3dflag
    character(CL) :: fileName
    character(CL) :: sfldName
    type(mct_avect) :: avtmp
    character(len=32) :: lstr
    character(len=32) :: bstr
    logical :: fileopen
    character(CL) :: currfile
    integer :: vtype
    integer(in) :: ndims
    integer(in),pointer :: dimid(:)
    type(file_desc_t) :: pioid
    type(var_desc_t) :: varid
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
       localCopy = (shr_file_queryPrefix(path) /= shr_file_noPrefix )
       if (localCopy) then
          call shr_file_get(rCode,fn,trim(path)//fn)
          fileName = fn
       else                 ! DON'T acquire a local copy of the data file
          fileName = trim(path)//fn
       end if
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
          rcode = pio_inq_vartype(pioid, varid, vtype)
          frame = nt
          call pio_setframe(pioid,varid,frame)
          if(vtype == PIO_DOUBLE) then
             call pio_read_darray(pioid,varid,pio_iodesc_r8,av%rattr(k,:),rcode)
          else if(vtype == PIO_REAL) then
             allocate(rdata(size(av%rattr(k,:))))
             call pio_read_darray(pioid,varid,pio_iodesc_r4,rdata,rcode)
             av%rattr(k,:) = real(rdata, kind=r8)
             deallocate(rdata)
          else
             allocate(idata(size(av%rattr(k,:))))
             call pio_read_darray(pioid,varid,pio_iodesc_int,idata,rcode)
             av%rattr(k,:) = real(idata, kind=r8)
             deallocate(idata)
          endif
       enddo

       call t_stopf(trim(lstr)//'_readpio')

    endif

  end subroutine shr_dmodel_readstrm

  !===============================================================================

  subroutine shr_dmodel_readstrm_fullfile(stream, pio_subsystem, pio_iotype, &
       gsMap, av, avFile, mpicom, &
       path, fn, nt, istr, boundstr)

    use shr_file_mod, only : shr_file_noprefix, shr_file_queryprefix, shr_file_get
    use shr_stream_mod
    use shr_ncread_mod, only: shr_ncread_open, shr_ncread_close, shr_ncread_varDimSizes
    use shr_ncread_mod, only: shr_ncread_tField

    !----- arguments -----
    type      (shr_stream_streamType) ,intent(inout)         :: stream
    type      (iosystem_desc_t)       ,intent(inout), target :: pio_subsystem
    integer   (IN)                    ,intent(in)            :: pio_iotype
    type      (mct_gsMap)             ,intent(in)            :: gsMap
    type      (mct_aVect)             ,intent(inout)         :: av
    type      (mct_aVect)             ,intent(inout)         :: avFile
    integer   (IN)                    ,intent(in)            :: mpicom
    character (len=*)                 ,intent(in)            :: path
    character (len=*)                 ,intent(in)            :: fn
    integer   (IN)                    ,intent(in)            :: nt
    character (len=*)                 ,intent(in) ,optional  :: istr
    character(len=*)                  ,intent(in) ,optional  :: boundstr

    !----- local -----
    integer(IN)                   :: my_task
    integer(IN)                   :: master_task
    integer(IN)                   :: ierr
    logical                       :: localCopy,fileexists
    integer(IN)                   :: gsize,nx,ny,nz
    integer(IN)                   :: k
    integer(IN)                   :: rCode   ! return code
    character(CL)                 :: fileName
    character(CL)                 :: sfldName
    character(len=32)             :: lstr
    character(len=32)             :: bstr
    logical                       :: fileopen
    character(CL)                 :: currfile
    character(CXX)                :: fldList ! list of fields

    integer(in)                   :: ndims
    integer(in),pointer           :: dimid(:)
    type(file_desc_t)             :: pioid
    type(var_desc_t)              :: varid
    integer(kind=pio_offset_kind) :: frame
    type(io_desc_t)               :: pio_iodesc_r8_local
    type(io_desc_t)               :: pio_iodesc_r4_local
    type(io_desc_t)               :: pio_iodesc_int_local
    integer(IN)                   :: avFile_beg, avFile_end

    real(r4), allocatable         :: rdata(:)
    integer, allocatable          :: idata(:)
    integer                       :: lsize, cnt,m,n
    integer, allocatable          :: count(:), compDOF(:)
    integer, pointer,dimension(:) :: gsmOP   ! gsmap ordered points
    integer                       :: vtype ! variable type on stream file
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
       localCopy = (shr_file_queryPrefix(path) /= shr_file_noPrefix )
       if (localCopy) then
          call shr_file_get(rCode,fn,trim(path)//fn)
          fileName = fn
       else                 ! DON'T acquire a local copy of the data file
          fileName = trim(path)//fn
       end if
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
          call pio_initdecomp(pio_subsystem, pio_double, count, compDOF, pio_iodesc_r8_local)
          call pio_initdecomp(pio_subsystem, pio_real, count, compDOF, pio_iodesc_r4_local)
          call pio_initdecomp(pio_subsystem, pio_int, count, compDOF, pio_iodesc_int_local)

          ! For each attribute, read all frames in one go
          frame = 1
          do k = 1, mct_aVect_nRAttr(avFile)
             if (my_task == master_task) then
                call shr_stream_getFileFieldName(stream,k,sfldName)
             endif
             call shr_mpi_bcast(sfldName,mpicom,'sfldName')
             rcode = pio_inq_varid(pioid,trim(sfldName),varid)
             rcode = pio_inq_vartype(pioid, varid, vtype)
             call pio_setframe(pioid,varid,frame)
             if(vtype == PIO_DOUBLE) then
                call pio_read_darray(pioid, varid, pio_iodesc_r8_local, avFile%rattr(k,:), rcode)
             else if(vtype == PIO_REAL) then
                allocate(rdata(size(avFile%rattr(k,:))))
                call pio_read_darray(pioid, varid, pio_iodesc_r4_local, rdata, rcode)
                avFile%rattr(k,:) = real(rdata, kind=r8)
                deallocate(rdata)
             else if(vtype == PIO_INT) then
                allocate(idata(size(avFile%rattr(k,:))))
                call pio_read_darray(pioid, varid, pio_iodesc_int_local, idata, rcode)
                avFile%rattr(k,:) = real(idata, kind=r8)
                deallocate(idata)
             endif
          enddo

          call pio_freedecomp(pio_subsystem, pio_iodesc_r8_local)
          call pio_freedecomp(pio_subsystem, pio_iodesc_r4_local)
          call pio_freedecomp(pio_subsystem, pio_iodesc_int_local)

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

    ! !DESCRIPTION:
    !    Returns TRUE if two domains are the the same (within tolerance).

    ! !INPUT/OUTPUT PARAMETERS:
    type(mct_gGrid)     ,intent(in)  :: ggrid1   ! 1st ggrid
    type(mct_gsmap)     ,intent(in)  :: gsmap1   ! 1st gsmap
    type(mct_gGrid)     ,intent(in)  :: ggrid2   ! 2nd ggrid
    type(mct_gsmap)     ,intent(in)  :: gsmap2   ! 2nd gsmap
    integer(IN)         ,intent(in)  :: method   ! selects what to compare
    integer(IN)         ,intent(in)  :: mpicom   ! mpicom
    real(R8)   ,optional,intent(in)  :: eps      ! epsilon compare value

    !--- local ---
    real(R8)    :: leps         ! local epsilon
    integer(IN) :: n            ! counters
    integer(IN) :: my_task,master_task
    integer(IN) :: gsize
    integer(IN) :: ierr
    integer(IN) :: nlon1, nlon2, nlat1, nlat2, nmask1, nmask2  ! av field indices
    logical     :: compare      ! local compare logical
    real(R8)    :: lon1,lon2    ! longitudes to compare
    real(R8)    :: lat1,lat2    ! latitudes to compare
    real(R8)    :: msk1,msk2    ! masks to compare
    integer(IN) :: nx,ni1,ni2   ! i grid size, i offset for 1 vs 2 and 2 vs 1
    integer(IN) :: n1,n2,i,j    ! local indices
    type(mct_aVect) :: avG1     ! global av
    type(mct_aVect) :: avG2     ! global av
    integer(IN) :: lmethod      ! local method
    logical     :: maskmethod, maskpoint ! masking on method
    character(*),parameter :: subName = '(shr_dmodel_gGridCompare) '
    character(*),parameter :: F01     = "('(shr_dmodel_gGridCompare) ',4a)"
    !-------------------------------------------------------------------------------
    !
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
                   if (lmethod == shr_dmodel_gGridCompareXYabs      ) then
                      if (abs(lon1 - lon2) > leps) compare = .false.
                      if (abs(lat1 - lat2) > leps) compare = .false.
                   else if (lmethod == shr_dmodel_gGridCompareXYrel      ) then
                      if (rdiff(lon1,lon2) > leps) compare = .false.
                      if (rdiff(lat1,lat2) > leps) compare = .false.
                   else if (lmethod == shr_dmodel_gGridCompareMaskIdent  ) then
                      if (msk1 /= msk2)compare = .false.
                   else if (lmethod == shr_dmodel_gGridCompareMaskZeros  ) then
                      if (msk1 == 0 .and. msk2 /= 0) compare = .false.
                      if (msk1 /= 0 .and. msk2 == 0) compare = .false.
                   else if (lmethod == shr_dmodel_gGridCompareMaskSubset ) then
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

    return

    !-------------------------------------------------------------------------------
  contains   ! internal subprogram
    !-------------------------------------------------------------------------------

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

  subroutine shr_dmodel_mapSet_global(smatp,&
       ggridS,gsmapS,nxgS,nygS,&
       ggridD,gsmapD,nxgD,nygD, &
       name,type,algo,mask,vect,&
       compid,mpicom,strategy)

    use shr_map_mod

    !----- arguments -----
    type(mct_sMatP), intent(inout) :: smatp
    type(mct_gGrid), intent(in)    :: ggridS
    type(mct_gsmap), intent(in)    :: gsmapS
    integer(IN)    , intent(in)    :: nxgS
    integer(IN)    , intent(in)    :: nygS
    type(mct_gGrid), intent(in)    :: ggridD
    type(mct_gsmap), intent(in)    :: gsmapD
    integer(IN)    , intent(in)    :: nxgD
    integer(IN)    , intent(in)    :: nygD
    character(len=*),intent(in)    :: name
    character(len=*),intent(in)    :: type
    character(len=*),intent(in)    :: algo
    character(len=*),intent(in)    :: mask
    character(len=*),intent(in)    :: vect
    integer(IN)    , intent(in)    :: compid
    integer(IN)    , intent(in)    :: mpicom
    character(len=*),intent(in),optional :: strategy

    !----- local -----

    integer(IN) :: n,i,j
    integer(IN) :: lsizeS,gsizeS,lsizeD,gsizeD
    integer(IN) :: nlon,nlat,nmsk
    integer(IN) :: my_task,master_task,ierr

    real(R8)   , pointer :: Xsrc(:,:)
    real(R8)   , pointer :: Ysrc(:,:)
    integer(IN), pointer :: Msrc(:,:)
    real(R8)   , pointer :: Xdst(:,:)
    real(R8)   , pointer :: Ydst(:,:)
    integer(IN), pointer :: Mdst(:,:)
    type(shr_map_mapType) :: shrmap
    type(mct_aVect) :: AVl
    type(mct_aVect) :: AVg

    character(len=32) :: lstrategy
    integer(IN) :: nsrc,ndst,nwts
    integer(IN), pointer :: isrc(:)
    integer(IN), pointer :: idst(:)
    real(R8)   , pointer :: wgts(:)
    type(mct_sMat) :: sMat0
    character(*), parameter :: subname = '(shr_dmodel_mapSet_global) '
    character(*), parameter :: F00   = "('(shr_dmodel_mapSet_global) ',8a)"
    character(*), parameter :: F01   = "('(shr_dmodel_mapSet_global) ',a,5i8)"

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

       call shr_map_getARptr(shrmap,isrc,idst,wgts)
       call mct_sMat_init(sMat0,ndst,nsrc,nwts)

       call mct_sMat_ImpGColI (sMat0,isrc,nwts)
       call mct_sMat_ImpGRowI (sMat0,idst,nwts)
       call mct_sMat_ImpMatrix(sMat0,wgts,nwts)
       call shr_map_clean(shrmap)
    endif

    call mct_sMatP_Init(sMatP,sMat0,gsmapS,gsmapD,lstrategy,master_task,mpicom,compid)

    if (my_task == master_task) then
       call mct_sMat_clean(sMat0)
    endif

  end subroutine shr_dmodel_mapSet_global

  !===============================================================================

  subroutine shr_dmodel_mapSet_dest(smatp,ggridS,gsmapS,nxgS,nygS,ggridD,gsmapD,nxgD,nygD, &
       name,type,algo,mask,vect,compid,mpicom,strategy)

    use shr_map_mod

    !----- arguments -----
    type(mct_sMatP), intent(inout) :: smatp
    type(mct_gGrid), intent(in)    :: ggridS
    type(mct_gsmap), intent(in)    :: gsmapS
    integer(IN)    , intent(in)    :: nxgS
    integer(IN)    , intent(in)    :: nygS
    type(mct_gGrid), intent(in)    :: ggridD
    type(mct_gsmap), intent(in)    :: gsmapD
    integer(IN)    , intent(in)    :: nxgD
    integer(IN)    , intent(in)    :: nygD
    character(len=*),intent(in)    :: name
    character(len=*),intent(in)    :: type
    character(len=*),intent(in)    :: algo
    character(len=*),intent(in)    :: mask
    character(len=*),intent(in)    :: vect
    integer(IN)    , intent(in)    :: compid
    integer(IN)    , intent(in)    :: mpicom
    character(len=*),intent(in),optional :: strategy

    !----- local -----

    integer(IN) :: n,i,j
    integer(IN) :: lsizeS,gsizeS,lsizeD
    integer(IN) :: nlon,nlat,nmsk
    integer(IN) :: my_task,master_task,ierr

    real(R8)   , pointer :: Xsrc(:,:)
    real(R8)   , pointer :: Ysrc(:,:)
    integer(IN), pointer :: Msrc(:,:)
    real(R8)   , pointer :: Xdst(:)
    real(R8)   , pointer :: Ydst(:)
    integer(IN), pointer :: Mdst(:)
    type(shr_map_mapType) :: shrmap
    type(mct_aVect) :: AVl
    type(mct_aVect) :: AVg

    character(len=32) :: lstrategy
    integer(IN) :: nsrc,ndst,nwts
    integer(IN), pointer :: points(:)
    integer(IN), pointer :: isrc(:)
    integer(IN), pointer :: idst(:)
    real(R8)   , pointer :: wgts(:)
    type(mct_sMat) :: sMat0

    character(*), parameter :: subname = '(shr_dmodel_mapSet_dest) '
    character(*), parameter :: F00   = "('(shr_dmodel_mapSet_dest) ',8a)"
    character(*), parameter :: F01   = "('(shr_dmodel_mapSet_dest) ',a,5i8)"

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

    allocate(Xsrc(nxgS,nygS),Ysrc(nxgS,nygS),Msrc(nxgS,nygS))
    if (my_task == master_task) then
       gsizeS = mct_aVect_lsize(AVg)
       if (gsizeS /= nxgS*nygS) then
          write(logunit,F01) ' ERROR: gsizeS ',gsizeS,nxgS,nygS
          call shr_sys_abort(subname//' ERROR gsizeS')
       endif

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
    call shr_mpi_bcast(Xsrc,mpicom)
    call shr_mpi_bcast(Ysrc,mpicom)
    call shr_mpi_bcast(Msrc,mpicom)

    if (my_task == master_task) call mct_aVect_clean(AVg)
    call mct_aVect_clean(AVl)

    !--- get sizes and allocate for DST ---

    lsizeD = mct_aVect_lsize(ggridD%data)
    call mct_avect_init(AVl,rList='lon:lat:mask',lsize=lsizeD)
    call mct_avect_copy(ggridD%data,AVl,rList='lon:lat:mask')

    allocate(Xdst(lsizeD),Ydst(lsizeD),Mdst(lsizeD))

    nlon = mct_avect_indexRA(AVl,'lon')
    nlat = mct_avect_indexRA(AVl,'lat')
    nmsk = mct_avect_indexRA(AVl,'mask')

    Mdst = 1
    do n = 1,lsizeD
       Xdst(n) = AVl%rAttr(nlon,n)
       Ydst(n) = AVl%rAttr(nlat,n)
       if (abs(AVl%rAttr(nmsk,n)) < 0.5_R8) Mdst(n) = 0
    enddo

    call mct_aVect_clean(AVl)

    !--- set map ---

    nsrc = nxgS*nygS
    ndst = nxgD*nygD
    call mct_gsmap_orderedPoints(gsmapD,my_task,points)
    if (size(points) /= size(Xdst)) then
       write(logunit,F01) ' ERROR: gsizeD ',size(points),size(Xdst)
       call shr_sys_abort(subname//' ERROR points size')
    endif
    call shr_map_mapSet(shrmap,Xsrc,Ysrc,Msrc,Xdst,Ydst,Mdst,ndst,points, &
         trim(name),trim(type),trim(algo),trim(mask),trim(vect))
    deallocate(points)
    deallocate(Xsrc,Ysrc,Msrc)
    deallocate(Xdst,Ydst,Mdst)

    !--- convert map to sMatP ---

    lstrategy = 'Xonly'
    if (present(strategy)) then
       lstrategy = trim(strategy)
    endif

    call shr_map_get(shrmap,shr_map_fs_nwts,nwts)
    allocate(isrc(nwts),idst(nwts),wgts(nwts))
    call shr_map_get(shrmap,isrc,idst,wgts)
    call shr_map_clean(shrmap)

    call mct_sMat_init(sMat0,ndst,nsrc,nwts)

    call mct_sMat_ImpLColI (sMat0,isrc,nwts)
    call mct_sMat_ImpLRowI (sMat0,idst,nwts)
    call mct_sMat_ImpMatrix(sMat0,wgts,nwts)
    deallocate(isrc,idst,wgts)

    call mct_sMatP_Init(sMatP,sMat0,gsmapS,gsmapD,master_task,mpicom,compid)

    call mct_sMat_clean(sMat0)

  end subroutine shr_dmodel_mapSet_dest

  !===============================================================================

  subroutine shr_dmodel_rearrGGrid( ggridi, ggrido, gsmap, rearr, mpicom )

    !----- arguments -----
    type(mct_ggrid), intent(in)    :: ggridi
    type(mct_ggrid), intent(inout) :: ggrido
    type(mct_gsmap), intent(in)    :: gsmap
    type(mct_rearr), intent(in)    :: rearr
    integer(IN)    , intent(in)    :: mpicom

    !----- local -----
    integer(IN)          :: lsize      ! lsize
    real(R8)   , pointer :: data(:)    ! temporary
    integer(IN), pointer :: idata(:)   ! temporary
    integer(IN)          :: my_task    ! local pe number
    integer(IN)          :: ier        ! error code
    character(*), parameter :: subname = '(shr_dmodel_rearrGGrid) '

    !-------------------------------------------------------------------------------
    ! PURPOSE:  Determine MCT ggrid
    !-------------------------------------------------------------------------------

    ! Initialize mct ggrid type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (ocean), 0 (non-ocean)

    call mpi_comm_rank(mpicom,my_task,ier)

    lsize = mct_gsMap_lsize(gsMap, mpicom)
    call mct_gGrid_init( ggrido, ggridi, lsize=lsize )

    ! Determine global gridpoint number attribute, GlobGridNum, automatically in ggrid

    call mct_gsMap_orderedPoints(gsMap, my_task, idata)
    call mct_gGrid_importIAttr(ggrido,'GlobGridNum',idata,lsize)
    deallocate(idata)

    ! Initialize attribute vector with special value

    allocate(data(lsize))

    data(:) = -9999.0_R8
    call mct_gGrid_importRAttr(ggrido,"lat"  ,data,lsize)
    call mct_gGrid_importRAttr(ggrido,"lon"  ,data,lsize)
    call mct_gGrid_importRAttr(ggrido,"area" ,data,lsize)
    call mct_gGrid_importRAttr(ggrido,"aream",data,lsize)
    data(:) = 0.0_R8
    call mct_gGrid_importRAttr(ggrido,"mask",data,lsize)
    call mct_gGrid_importRAttr(ggrido,"frac",data,lsize)

    deallocate(data)

    call mct_rearr_rearrange(ggridi%data, ggrido%data, rearr)

  end subroutine shr_dmodel_rearrGGrid

  !===============================================================================

  subroutine shr_dmodel_translateAV( avi, avo, avifld, avofld, rearr )

    !----- arguments -----
    type(mct_aVect), intent(in)    :: avi       ! input av
    type(mct_aVect), intent(inout) :: avo       ! output av
    character(len=*),intent(in)    :: avifld(:) ! input field names for translation
    character(len=*),intent(in)    :: avofld(:) ! output field names for translation
    type(mct_rearr), intent(in),optional :: rearr     ! rearranger for diff decomp

    !----- local -----
    integer(IN)      :: k,ka,kb,kc,cnt  ! indices
    integer(IN)      :: lsize      ! lsize
    integer(IN)      :: nflds      ! number of fields in avi

    type(mct_aVect)  :: avtri,avtro ! translated av on input/output grid
    character(CXX)   :: ilist      ! input list for translation
    character(CXX)   :: olist      ! output list for translation
    character(CX)    :: cfld       ! character field name
    type(mct_string) :: sfld       ! string field
    integer(IN)      :: ktrans
    character(*), parameter :: subname = '(shr_dmodel_translateAV) '

    !-------------------------------------------------------------------------------
    ! PURPOSE:  Fill avo from avi
    !-------------------------------------------------------------------------------

    if (size(avifld) /= size(avofld)) then
       write(logunit,*) subname,' ERROR": avi and avo fld list ',size(avifld),size(avofld)
       call shr_sys_flush(logunit)
    endif
    ktrans = size(avifld)

    ! generate fld lists
    nflds = mct_aVect_nRattr(avi)
    cnt = 0
    do ka = 1,nflds
       call mct_aVect_getRList(sfld,ka,avi)
       cfld = mct_string_toChar(sfld)
       call mct_string_clean(sfld)

       k = 0
       kb = 0
       kc = 0
       do while (kb == 0 .and. k < ktrans)
          k = k + 1
          if (trim(avifld(k)) == trim(cfld)) then
             kb = k
             kc = mct_aVect_indexRA(avo,trim(avofld(kb)),perrWith='quiet')
             if (ka > 0 .and. kc > 0) then
                cnt = cnt + 1
                if (cnt == 1) then
                   ilist = trim(avifld(kb))
                   olist = trim(avofld(kb))
                else
                   ilist = trim(ilist)//':'//trim(avifld(kb))
                   olist = trim(olist)//':'//trim(avofld(kb))
                endif
             endif
          endif
       enddo
    enddo

    if (cnt > 0) then
       lsize = mct_avect_lsize(avi)
       call mct_avect_init(avtri,rlist=olist,lsize=lsize)
       call mct_avect_Copy(avi,avtri,rList=ilist,TrList=olist)

       if (present(rearr)) then
          lsize = mct_avect_lsize(avo)
          call mct_avect_init(avtro,rlist=olist,lsize=lsize)
          call mct_avect_zero(avtro)
          call mct_rearr_rearrange(avtri, avtro, rearr)
          call mct_avect_Copy(avtro,avo)
          call mct_aVect_clean(avtro)
       else
          call mct_avect_Copy(avtri,avo)
       endif

       call mct_aVect_clean(avtri)
    endif

  end subroutine shr_dmodel_translateAV

  !===============================================================================

  subroutine shr_dmodel_translate_list( avi, avo, avifld, avofld, ilist, olist, cnt)

    !----- arguments -----
    type(mct_aVect), intent(in)    :: avi       ! input av
    type(mct_aVect), intent(inout) :: avo       ! output av
    character(len=*),intent(in)    :: avifld(:) ! input field names for translation
    character(len=*),intent(in)    :: avofld(:) ! output field names for translation
    character(CL)   ,intent(out)   :: ilist     ! input list for translation
    character(CL)   ,intent(out)   :: olist     ! output list for translation
    integer(IN)     ,intent(out)   :: cnt       ! indices


    !----- local -----
    integer(IN)      :: k,ka,kb,kc ! indices
    integer(IN)      :: nflds        ! number of fields in avi
    character(CL)    :: cfld         ! character field name
    type(mct_string) :: sfld         ! string field
    integer(IN)      :: ktrans
    character(*), parameter :: subname = '(shr_dmodel_translateAV) '

    !-------------------------------------------------------------------------------
    ! PURPOSE:  Fill avo from avi
    !-------------------------------------------------------------------------------

    if (size(avifld) /= size(avofld)) then
       write(logunit,*) subname,' ERROR": avi and avo fld list ',size(avifld),size(avofld)
       call shr_sys_flush(logunit)
    endif
    ktrans = size(avifld)

    ! generate fld lists
    nflds = mct_aVect_nRattr(avi)
    cnt = 0
    do ka = 1,nflds
       call mct_aVect_getRList(sfld,ka,avi)
       cfld = mct_string_toChar(sfld)
       call mct_string_clean(sfld)

       k = 0
       kb = 0
       kc = 0
       do while (kb == 0 .and. k < ktrans)
          k = k + 1
          if (trim(avifld(k)) == trim(cfld)) then
             kb = k
             kc = mct_aVect_indexRA(avo,trim(avofld(kb)),perrWith='quiet')
             if (ka > 0 .and. kc > 0) then
                cnt = cnt + 1
                if (cnt == 1) then
                   ilist = trim(avifld(kb))
                   olist = trim(avofld(kb))
                else
                   ilist = trim(ilist)//':'//trim(avifld(kb))
                   olist = trim(olist)//':'//trim(avofld(kb))
                endif
             endif
          endif
       enddo
    enddo

  end subroutine shr_dmodel_translate_list

  !===============================================================================

  subroutine shr_dmodel_translateAV_list( avi, avo, ilist, olist, rearr )

    !----- arguments -----
    type(mct_aVect), intent(in)    :: avi       ! input av
    type(mct_aVect), intent(inout) :: avo       ! output av
    character(CL)   ,intent(in)    :: ilist     ! input list for translation
    character(CL)   ,intent(in)    :: olist     ! output list for translation
    type(mct_rearr), intent(in),optional :: rearr  ! rearranger for diff decomp

    !----- local -----
    integer(IN)      :: lsize       ! lsize
    type(mct_aVect)  :: avtri,avtro ! translated av on input/output grid
    character(*), parameter :: subname = '(shr_dmodel_translateAV) '

    !-------------------------------------------------------------------------------
    ! PURPOSE:  Fill avo from avi
    !-------------------------------------------------------------------------------

    lsize = mct_avect_lsize(avi)
    call mct_avect_init(avtri,rlist=olist,lsize=lsize)
    call mct_avect_Copy(avi,avtri,rList=ilist,TrList=olist)

    if (present(rearr)) then
       lsize = mct_avect_lsize(avo)
       call mct_avect_init(avtro,rlist=olist,lsize=lsize)
       call mct_avect_zero(avtro)
       call mct_rearr_rearrange(avtri, avtro, rearr)
       call mct_avect_Copy(avtro,avo)
       call mct_aVect_clean(avtro)
    else
       call mct_avect_Copy(avtri,avo)
    endif

    call mct_aVect_clean(avtri)

  end subroutine shr_dmodel_translateAV_list

  !===============================================================================

end module shr_dmodel_mod
