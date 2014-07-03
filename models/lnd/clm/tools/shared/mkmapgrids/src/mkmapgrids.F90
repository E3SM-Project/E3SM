program mkmapgrids

!-----------------------------------------------------------------------
!BOP
!
! !PROGRAM: mkmapgrids
!
! !DESCRIPTION:
! Routines to create land model grid
!
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use domainMod   

  implicit none

  ! !LOCAL VARIABLES:
  type(domain_type) :: tdomain            ! local domain
  character(len=256):: locfn              ! local dataset file name
  character(len=256):: fname_in, fname_out 
  integer :: ier

  ! read input and write domain in scrip input format

  namelist /mkmapgrids_in/ fname_in, fname_out          	
  read(5, mkmapgrids_in, iostat=ier)
  if (ier /= 0) then
     write(6,*)'error: namelist input resulted in error code ',ier
     call abort()
  endif
  
  call read_domain(tdomain,trim(fname_in))
  call write_domain(tdomain, trim(fname_in), grid_file_out=trim(fname_out))

contains

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_domain
!
! !INTERFACE:
  subroutine read_domain(domain,fname)
!
! !DESCRIPTION:
! Read a grid file
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain
    character(len=*) ,intent(in)    :: fname
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    include 'netcdf.inc'
    real(r8), parameter :: flandmin = 0.001 ! minimum land frac for land cell
    integer :: nlon,nlat                    ! size
    character(len=256) :: locfn             ! local file name
    integer :: i,j                          ! indexes
    integer :: ncid                         ! netCDF file id
    integer :: dimid                        ! netCDF dimension id
    integer :: varid                        ! netCDF variable id
    logical :: dimset                       ! local ni,nj
    logical :: lonlatset                    ! local lon(:,:), lat(:,:)
    logical :: edgeNESWset                  ! local EDGE[NESW]
    logical :: llneswset                    ! local lat[ns],lon[we]
    logical :: landfracset                  ! local landfrac
    logical :: maskset                      ! local mask
    integer :: ndims                        ! number of dims for variable
    integer :: ier                          ! error status
    real(r8):: edgen                        ! edge north
    real(r8):: edgee                        ! edge east
    real(r8):: edges                        ! edge south
    real(r8):: edgew                        ! edge west
    character(len= 32) :: subname = 'read_domain'
!-----------------------------------------------------------------

    write(6,*) ' ' 

    dimset      = .false. 
    lonlatset   = .false. 
    edgeNESWset = .false. 
    llneswset   = .false. 
    landfracset = .false. 
    maskset     = .false. 

    ! Read domain file and compute stuff as needed

    call check_ret(nf_open(fname, 0, ncid), subname)

    ier = nf_inq_dimid (ncid, 'lon', dimid)
    if (ier == NF_NOERR) then
       if (dimset) write(6,*) trim(subname),' WARNING, overwriting dims'
       dimset = .true.
       write(6,*) trim(subname),' read lon and lat dims'
       call check_ret(nf_inq_dimid  (ncid, 'lon', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlon), subname)
       call check_ret(nf_inq_dimid  (ncid, 'lat', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlat), subname)
    endif

    ier = nf_inq_dimid (ncid, 'ni', dimid)
    if (ier == NF_NOERR) then
       if (dimset) write(6,*) trim(subname),' WARNING, overwriting dims'
       dimset = .true.
       write(6,*) trim(subname),' read ni and nj dims'
       call check_ret(nf_inq_dimid  (ncid, 'ni', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlon), subname)
       call check_ret(nf_inq_dimid  (ncid, 'nj', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlat), subname)
    endif

    ier = nf_inq_dimid (ncid, 'lsmlon', dimid)
    if (ier == NF_NOERR) then
       if (dimset) write(6,*) trim(subname),' WARNING, overwriting dims'
       dimset = .true.
       write(6,*) trim(subname),' read lsmlon and lsmlat dims'
       call check_ret(nf_inq_dimid  (ncid, 'lsmlon', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlon), subname)
       call check_ret(nf_inq_dimid  (ncid, 'lsmlat', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlat), subname)
    endif

    if (dimset) then
       write(6,*) trim(subname),' initialized domain'
       call domain_init(domain,nlon*nlat,nlon,nlat)
    else
       write(6,*) trim(subname),' ERROR: ns not set for domain_init'
       stop
    endif

    ier = nf_inq_varid (ncid, 'xc', varid)
    if (ier == NF_NOERR) then
       if (lonlatset) write(6,*) trim(subname),' WARNING, overwriting lat,lon'
       lonlatset = .true.
       write(6,*) trim(subname),' read xc and yc fields'
       call check_ret(nf_inq_varid (ncid, 'xc', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%lonc), subname)
       call check_ret(nf_inq_varid (ncid, 'yc', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%latc), subname)
    endif

    ier = nf_inq_varid (ncid, 'lon', varid)
    if (ier == NF_NOERR) then
       if (lonlatset) write(6,*) trim(subname),' WARNING, overwriting lat,lon'
       lonlatset = .true.
       ier = nf_inq_varndims(ncid,varid,ndims)
       write(6,*) trim(subname),' read lon and lat fields'
       call check_ret(nf_inq_varid (ncid, 'lat', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%latc), subname)
       call check_ret(nf_inq_varid (ncid, 'lon', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%lonc), subname)
    endif

    ier = nf_inq_varid (ncid, 'LATIXY', varid)
    if (ier == NF_NOERR) then
       if (lonlatset) write(6,*) trim(subname),' WARNING, overwriting lat,lon'
       lonlatset = .true.
       write(6,*) trim(subname),' read LONGXY and LATIXY fields'
       call check_ret(nf_inq_varid (ncid, 'LONGXY', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%lonc), subname)
       call check_ret(nf_inq_varid (ncid, 'LATIXY', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%latc), subname)
    endif

    ier = nf_inq_varid (ncid, 'LATN', varid)
    if (ier == NF_NOERR) then
       if (llneswset) write(6,*) trim(subname),' WARNING, overwriting lat[ns],lon[we]'
       llneswset = .true.
       write(6,*) trim(subname),' read LAT[NS],LON[WE]'
       call check_ret(nf_inq_varid (ncid, 'LATN', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%latn), subname)
       call check_ret(nf_inq_varid (ncid, 'LONE', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%lone), subname)
       call check_ret(nf_inq_varid (ncid, 'LATS', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%lats), subname)
       call check_ret(nf_inq_varid (ncid, 'LONW', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%lonw), subname)
    endif

    domain%frac = 1.
    ier = nf_inq_varid (ncid, 'frac', varid)
    if (ier == NF_NOERR) then
       if (landfracset) write(6,*) trim(subname),' WARNING, overwriting frac'
       landfracset = .true.
       write(6,*) trim(subname),' read frac'
       call check_ret(nf_inq_varid (ncid, 'frac', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%frac), subname)
    endif
    ier = nf_inq_varid (ncid, 'LANDFRAC', varid)
    if (ier == NF_NOERR) then
       if (landfracset) write(6,*) trim(subname),' WARNING, overwriting frac'
       landfracset = .true.
       write(6,*) trim(subname),' read LANDFRAC'
       call check_ret(nf_inq_varid (ncid, 'LANDFRAC', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%frac), subname)
    endif

    domain%mask = 1.
    ier = nf_inq_varid (ncid, 'mask', varid)
    if (ier == NF_NOERR) then
       if (maskset) write(6,*) trim(subname),' WARNING, overwriting mask'
       maskset = .true.
       write(6,*) trim(subname),' read mask'
       call check_ret(nf_inq_varid (ncid, 'mask', varid), subname)
       call check_ret(nf_get_var_int (ncid, varid, domain%mask), subname)
    endif
    ier = nf_inq_varid (ncid, 'LANDMASK', varid)
    if (ier == NF_NOERR) then
       if (maskset) write(6,*) trim(subname),' WARNING, overwriting mask'
       maskset = .true.
       write(6,*) trim(subname),' read LANDMASK'
       call check_ret(nf_inq_varid (ncid, 'LANDMASK', varid), subname)
       call check_ret(nf_get_var_int (ncid, varid, domain%mask), subname)
    endif

    call check_ret(nf_close(ncid), subname)

    if (.not.lonlatset) then
       write(6,*) trim(subname),' ERROR: long and lati not set '
       call abort()
    endif

    if (.not.maskset.and.landfracset) then
       maskset = .true.
       where (domain%frac < flandmin)
          domain%mask = 0     !ocean
       elsewhere
          domain%mask = 1     !land
       endwhere
    endif

    if (.not.llneswset) then
       if (edgeneswset) then
          write(6,*) trim(subname),' compute lat[ns],lon[we] from edge[nesw]'
          call domain_celledge_regional(domain, edgen, edgee, edges, edgew)
       else
          write(6,*) trim(subname),' compute lat[ns],lon[we]'
          call domain_celledge_global(domain)
       endif
    endif

  end subroutine read_domain

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_domain
!
! !INTERFACE:
  subroutine write_domain(domain,fname,grid_file_out)
!
! !USES:
    use shr_sys_mod, only : shr_sys_getenv
!
! !DESCRIPTION:
! Write a scrip netcdf grid file

! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(domain_type),intent(inout) :: domain
    character(len=*) ,intent(in)    :: fname
    character(len=*) ,intent(in)    :: grid_file_out
!
! !REVISION HISTORY:
! Author: T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer ::  grid_rank=2, grid_corners=4

    integer , allocatable :: grid_imask(:)

    real(r8), allocatable :: grid_center_lat(:), &  ! lat/lon coordinates for
                             grid_center_lon(:)     ! each grid center in degrees
    
    real(r8), allocatable :: grid_corner_lat(:,:),& ! lat/lon coordinates for
                             grid_corner_lon(:,:)   ! each grid corner in degrees
    

    integer :: ncstat,          &  ! general netCDF status variable
               ncid,      &  ! netCDF grid dataset id
               gridsize_id,  &  ! netCDF grid size dim id
               gridcorn_id,  &  ! netCDF grid corner dim id
               gridrank_id,  &  ! netCDF grid rank dim id
               griddims_id,  &  ! netCDF grid dimension size id
               grdcntrlat_id,&  ! netCDF grid center lat id
               grdcntrlon_id,&  ! netCDF grid center lon id
               grdimask_id,  &  ! netCDF grid mask id
               grdcrnrlat_id,&  ! netCDF grid corner lat id
               grdcrnrlon_id    ! netCDF grid corner lon id
    
    integer :: dims2_id(2)  ! netCDF dim id array for 2-d arrays
    integer :: dims1_id(1) 
    integer :: nx,ny
    integer :: omode           ! netCDF output mode
    integer :: i,j,n,ns
    integer :: grid_size
    integer :: values(8)
    character(len= 8) :: date
    character(len=10) :: time
    character(len= 5) :: zone
    character(len=18) :: datetime
    character(len=16) :: logname
    character(len=16) :: hostname
    character(len=32) :: subname = 'write_domain'
    character(len=32) :: str
    character(len=256) :: version = &
    "$HeadURL: https://svn-ccsm-models.cgd.ucar.edu/clm2/trunk_tags/clm4_5_71/models/lnd/clm/tools/shared/mkmapgrids/src/mkmapgrids.F90 $"
    character(len=256) :: revision_id = &
    "$Id: mkmapgrids.F90 46863 2013-05-08 03:25:21Z sacks $"
    character(len=1500) :: hist
    integer :: ier, varid
    integer,  allocatable :: grid_dims(:)
!-----------------------------------------------------------------

    ! longitudes and latitudes of cell centers and corners.

    ns = domain%ns
    nx = domain%ni
    ny = domain%nj

    grid_size = ns
    allocate(grid_dims(grid_rank))
    allocate(grid_center_lat(grid_size))
    allocate(grid_center_lon(grid_size))
    allocate(grid_corner_lat(grid_corners,grid_size))
    allocate(grid_corner_lon(grid_corners,grid_size))
    allocate(grid_imask(grid_size))

    grid_dims(1) = nx
    grid_dims(2) = ny

    do n=1,ns
       grid_center_lat(n  ) = domain%latc(n)
       grid_corner_lat(1,n) = domain%lats(n)
       grid_corner_lat(2,n) = domain%lats(n)
       grid_corner_lat(3,n) = domain%latn(n) 
       grid_corner_lat(4,n) = domain%latn(n) 

       grid_center_lon(n  ) = domain%lonc(n)
       grid_corner_lon(1,n) = domain%lonw(n)
       grid_corner_lon(2,n) = domain%lone(n)
       grid_corner_lon(3,n) = domain%lone(n)
       grid_corner_lon(4,n) = domain%lonw(n)
       
       grid_imask(n) = domain%mask(n)       
    end do

    ! create netCDF dataset for this grid

    call check_ret(nf_create(trim(grid_file_out), NF_CLOBBER, ncid), subname)
    !BUG uncommenting the following results in 0 for grid_dims
    ! call check_ret(nf_set_fill(ncid, NF_NOFILL, omode), subname)

    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, 'title', &
         len_trim(grid_file_out), grid_file_out), subname)

    call date_and_time (date, time, zone, values)

    datetime(1:8) =        date(5:6) // '/' // date(7:8) // '/' // date(3:4)
    datetime(9:)  = ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6) // ' '

    call shr_sys_getenv ('LOGNAME', logname,  ier)
    call shr_sys_getenv ('HOST',    hostname, ier)

    str = 'NCAR-CESM:CF-1.0'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
       'Conventions', len_trim(str), trim(str)), subname)


    hist = datetime // trim (logname) // ':' // trim (hostname)

    call check_ret(nf_put_att_text (ncid, nf_global, 'history', len_trim(hist), trim(hist) ), subname )


    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, 'input_file', &
         len_trim(fname), trim(fname)), subname)

    call check_ret(nf_put_att_text (ncid, nf_global, 'version', &
         len_trim(version), version), subname)
    call check_ret(nf_put_att_text (ncid, nf_global, 'revision_Id', &
         len_trim(revision_id), revision_id), subname)

#ifdef OPT
    str = 'TRUE'
#else
    str = 'FALSE'
#endif

    call check_ret(nf_put_att_text (ncid, nf_global, 'Compiler_Optimized', &
         len_trim(str), str), subname)

    ! define grid size, grid corner and grid rank dimensions

    call check_ret(nf_def_dim(ncid, 'grid_size'   , grid_size   , gridsize_id), subname)
    call check_ret(nf_def_dim(ncid, 'grid_corners', grid_corners, gridcorn_id), subname)
    call check_ret(nf_def_dim(ncid, 'grid_rank   ', grid_rank   , gridrank_id), subname)

    ! define grid dimension size array
    call check_ret(nf_def_var (ncid, 'grid_dims', NF_INT, &
         1, gridrank_id, griddims_id), subname)

    ! define grid center latitude array
    call check_ret(nf_def_var (ncid, 'grid_center_lat', NF_DOUBLE, &
         1, gridsize_id, grdcntrlat_id), subname)

    call check_ret(nf_put_att_text (ncid, grdcntrlat_id, 'units', &
         7, 'degrees'), subname)
    
    ! define grid center longitude array
    call check_ret(nf_def_var (ncid, 'grid_center_lon', NF_DOUBLE, &
         1, gridsize_id, grdcntrlon_id), subname)

    call check_ret(nf_put_att_text (ncid, grdcntrlon_id, 'units', &
         7, 'degrees'), subname)

    ! define grid mask array
    call check_ret(nf_def_var (ncid, 'grid_imask', NF_INT, &
         1, gridsize_id, grdimask_id), subname)

    call check_ret(nf_put_att_text (ncid, grdimask_id, 'units', &
         8, 'unitless'), subname)

    ! define grid corner latitude array
    dims2_id(1) = gridcorn_id
    dims2_id(2) = gridsize_id
    
    call check_ret(nf_def_var (ncid, 'grid_corner_lat', NF_DOUBLE, &
         2, dims2_id, grdcrnrlat_id), subname)
    
    call check_ret(nf_put_att_text (ncid, grdcrnrlat_id, 'units', &
         7, 'degrees'), subname)
    
    ! define grid corner longitude array
    call check_ret(nf_def_var (ncid, 'grid_corner_lon', NF_DOUBLE, &
         2, dims2_id, grdcrnrlon_id), subname)
    
    call check_ret(nf_put_att_text (ncid, grdcrnrlon_id, 'units', &
         7, 'degrees'), subname)
    
    ! end definition stage
    call check_ret(nf_enddef(ncid), subname)

    ! write grid data

    call check_ret(nf_put_var_int   (ncid, griddims_id  , grid_dims      ), subname)
    call check_ret(nf_put_var_int   (ncid, grdimask_id  , grid_imask     ), subname)
    call check_ret(nf_put_var_double(ncid, grdcntrlat_id, grid_center_lat), subname)
    call check_ret(nf_put_var_double(ncid, grdcntrlon_id, grid_center_lon), subname)
    call check_ret(nf_put_var_double(ncid, grdcrnrlat_id, grid_corner_lat), subname)
    call check_ret(nf_put_var_double(ncid, grdcrnrlon_id, grid_corner_lon), subname)

    ! Close grid dataset
    call check_ret(nf_close(ncid), subname)
    
    write(6,*)'Successfully created SCRIP grid file'

  end subroutine write_domain

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_ret
!
! !INTERFACE:
  subroutine check_ret(ret, calling)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer, intent(in) :: ret
    character(len=*) :: calling
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------

    if (ret /= NF_NOERR) then
       write(6,*)'netcdf error from ',trim(calling), ' rcode = ', ret, &
                 ' error = ', NF_STRERROR(ret)
       call abort()
    end if

  end subroutine check_ret

end program mkmapgrids
