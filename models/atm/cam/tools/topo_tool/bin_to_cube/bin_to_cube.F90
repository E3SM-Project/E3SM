!
!  DATE CODED:   Nov 7, 2011
!
!  DESCRIPTION:  This program reads USGS 30-sec terrain dataset from NetCDF file and
!                bins it to an approximately 3km cubed-sphere grid and outputs the
!                data in netCDF format.
!
!                The LANDM_COSLAT field is read in from a separate netCDF file and linearly
!                interpolated to the 3km cubed-sphere grid.
!
!  Author: Peter Hjort Lauritzen (pel@ucar.edu) 
!
!  ROUTINES CALLED:
!       netcdf routines
!
!  COMPILING:
!
program convterr
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
#     include         <netcdf.inc>
  !
  integer :: im, jm
  
  integer,  parameter :: ncube = 3000 !dimension of cubed-sphere grid
!  integer,  parameter :: ncube = 540 !dimension of cubed-sphere grid
  !        integer,  parameter :: ncube = 361 ! for debugging
  
  integer*2,  allocatable, dimension(:,:) :: terr     ! global 30-sec terrain data
  integer*1,  allocatable, dimension(:,:) :: landfrac ! global 30-sec land fraction
  
  integer :: alloc_error,dealloc_error  
  integer :: i,j,n,k,index                                ! index
  integer*2, allocatable, dimension(:,:)  :: iterr        ! terrain data for 30-sec tile
  integer ncid,status, dimlatid,dimlonid, landid, topoid  ! for netCDF USGS data file
  integer :: srcid,dstid                                  ! for netCDF weight file
  
  real(r8), allocatable, dimension(:)   :: lon  , lat
  real(r8), allocatable, dimension(:)   :: lon_landm  , lat_landm
  real(r8), allocatable, dimension(:,:) :: landm_coslat
  integer :: im_landm, jm_landm
  integer :: lonid, latid
  integer :: lon_vid, lat_vid
  
  REAL    (r8), PARAMETER :: tiny  = 1.0E-10
  REAL    (r8), PARAMETER :: pi    = 3.14159265358979323846264338327
  REAL    (r8), PARAMETER :: piq   = 0.25*pi
  REAL    (r8), PARAMETER :: rad2deg   = 180.0/pi
  REAL    (r8), PARAMETER :: deg2rad   = pi/180.0
  
  real(r8) :: alpha, beta,da,wt,dlat
  integer  :: ipanel,icube,jcube
  real(r8), allocatable, dimension(:,:,:)   :: weight,terr_cube,landfrac_cube,sgh30_cube
  real(r8), allocatable, dimension(:,:,:)   :: landm_coslat_cube
  integer , allocatable, dimension(:,:)     :: idx,idy,idp
  !
  real(r8) :: dx,dy
  !
  ! for "bi-linear" interpolation
  !
  real(r8) :: lambda,theta,wx,wy
  integer :: ilon,ilat,ip1,jp1
  !
  ! variable for regridding
  !
  integer ::  src_grid_dim  ! for netCDF weight file
  !
  ! this is only used if target grid is a lat-lon grid
  !
  integer , parameter :: im_target = 360 , jm_target = 180
  logical , parameter :: ltarget_rll = .TRUE.
  !
  ! this is only used if target grid is not a lat-lon grid
  !
  real(r8), allocatable, dimension(:) :: lon_target, lat_target
  !
  ! compute volume of surface topography
  !
  real(r8) :: vol,dx_rad,vol_cube,area_latlon,darea_latlon       ! latitude array
  real(r8), allocatable, dimension(:,:) :: darea_cube

  !
  ! read in USGS data from netCDF file
  !
  !        status = nf_open('topo-lowres.nc', 0, ncid) !for debugging
  status = nf_open('usgs-rawdata.nc', 0, ncid)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
  
  status = NF_INQ_DIMID(ncid, 'lat', dimlatid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_INQ_DIMLEN(ncid, dimlatid, jm)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  status = NF_INQ_DIMID(ncid, 'lon', dimlonid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_INQ_DIMLEN(ncid, dimlonid, im)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  WRITE(*,*) "lon-lat dimensions: ",im,jm
  
  allocate ( landfrac(im,jm),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for landfrac'
    stop
  end if
  
  allocate ( terr(im,jm),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for terr'
    stop
  end if
  
  allocate ( lon(im),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for landfrac'
    stop
  end if
  
  allocate ( lat(jm),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for landfrac'
    stop
  end if
  
  terr = -999999
  landfrac = -99.0
  
  status = NF_INQ_VARID(ncid, 'landfract', landid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  status = NF_GET_VAR_INT1(ncid, landid,landfrac)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  WRITE(*,*) "min/max of 30sec land fraction",MINVAL(landfrac),MAXVAL(landfrac)
  
  
  status = NF_INQ_VARID(ncid, 'htopo', topoid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  WRITE(*,*) "read terrain data"
  status = NF_GET_VAR_INT2(ncid, topoid,terr)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  status = NF_INQ_VARID(ncid, 'lon', lonid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  WRITE(*,*) "read lon"
  status = NF_GET_VAR_DOUBLE(ncid, lonid,lon)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  status = NF_INQ_VARID(ncid, 'lat', latid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  WRITE(*,*) "read lat"
  status = NF_GET_VAR_DOUBLE(ncid, latid,lat)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  print *,"close file topo.nc"
  status = nf_close (ncid)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  WRITE(*,*) 'done reading in USGS data from netCDF file'
  
  WRITE(*,*) "Adjustments to land fraction: Extend land fraction for Ross Ice shelf by"
  WRITE(*,*) "setting all landfractions south of 79S to 1"
  DO j=1,jm
    IF (lat(j)<-79.0) THEN
      DO i=1,im
        landfrac(i,j) = 1
      END DO
    END IF
  END DO
  
  WRITE(*,*) "compute volume for USGS raw data"
  vol = 0.0
  dx = (lon(2)-lon(1))
  dx_rad = dx*deg2rad
  do j=1,jm
    do i=1,im
      darea_latlon = dx_rad*(SIN(deg2rad*(-90.0+dx*j))-SIN(deg2rad*(-90.0+dx*(j-1))))
      vol = vol+DBLE(terr(i,j))*darea_latlon
      area_latlon = area_latlon + darea_latlon
    end do
  end do
  vol = vol/area_latlon
  WRITE(*,*) "consistency of lat-lon area",area_latlon-4.0*pi
  WRITE(*,*) "volume of topography about sea-level (raw usgs data)",vol

  
  !
  !****************************************************
  !
  ! read LANDM_COSLAT
  !
  !****************************************************
  !
  WRITE(*,*) "read LANDM_COSLAT from file"
  status = nf_open('landm_coslat.nc', 0, ncid)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
  
  status = NF_INQ_DIMID(ncid, 'lat', dimlatid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_INQ_DIMLEN(ncid, dimlatid, jm_landm)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  status = NF_INQ_DIMID(ncid, 'lon', dimlonid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_INQ_DIMLEN(ncid, dimlonid, im_landm)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  WRITE(*,*) "lon-lat dimensions: ",im_landm,jm_landm
  
  allocate ( landm_coslat(im_landm,jm_landm),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for landfrac'
    stop
  end if
  
  allocate ( lon_landm(im_landm),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for landfrac'
    stop
  end if
  
  allocate ( lat_landm(jm_landm),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for landfrac'
    stop
  end if
  
  do j = 1, jm_landm
    do i = 1, im_landm
      landm_coslat(i,j) = -999999.99
    end do
  end do
  
  status = NF_INQ_VARID(ncid, 'LANDM_COSLAT', landid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  status = NF_GET_VAR_DOUBLE(ncid, landid,landm_coslat)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  WRITE(*,*) "min/max of landm_coslat",MINVAL(landm_coslat),MAXVAL(landm_coslat)
  
  status = NF_INQ_VARID(ncid, 'lon', lonid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  WRITE(*,*) "read lon"
  status = NF_GET_VAR_DOUBLE(ncid, lonid,lon_landm)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  status = NF_INQ_VARID(ncid, 'lat', latid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  WRITE(*,*) "read lat"
  status = NF_GET_VAR_DOUBLE(ncid, latid,lat_landm)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  print *,"close file"
  status = nf_close (ncid)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  WRITE(*,*) 'done reading in LANDM_COSLAT data from netCDF file'
  
  !
  ! bin data on cubed-sphere grid
  !
  da   = pi / DBLE(2*ncube)!equal-angle cubed-sphere grid spacing
  lon  = deg2rad*lon
  lat  = deg2rad*lat
  dlat = pi/DBLE(jm)
  allocate ( weight(ncube,ncube,6),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for weight'
    stop
  end if
  weight    = 0.0
  allocate ( terr_cube(ncube,ncube,6),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for terr_cube'
    stop
  end if
  terr_cube = 0.0
  allocate ( landfrac_cube(ncube,ncube,6),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for terr_cube'
    stop
  end if
  landfrac_cube = 0.0
  allocate ( landm_coslat_cube(ncube,ncube,6),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for terr_cube'
    stop
  end if
  landm_coslat_cube = 0.0
  
  
  allocate ( idx(im,jm),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for idx'
    stop
  end if
  allocate ( idy(im,jm),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for idy'
    stop
  end if
  allocate ( idp(im,jm),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for idp'
    stop
  end if
  
  WRITE(*,*) "bin lat-lon data on cubed-sphere"

  !
  ! for debugging ONLY
  !
!  DO j=1,jm
!    DO i=1,im
!!      terr(i,j) = 10000.0*(2.0+cos(lat(j))*cos(lat(j))*cos(2.0*lon(i)))!Y22
!!      terr(i,j) = 10000.0*(2.0+(sin(2.0*lat(j))**16)*cos(16.0*lon(i))) !Y16_32
!      terr(i,j) = 10000.0*(2.0+cos(16.0*lon(i))) !Y16_32
!    END DO
!  END DO

  DO j=1,jm
    DO i=1,im
!      WRITE(*,*) "bin to cube ",100.0*FLOAT(i+(j-1)*im)/FLOAT(im*jm),"% done"
      call CubedSphereABPFromRLL(lon(i), lat(j), alpha, beta, ipanel)            
      icube = CEILING((alpha + piq) / da)
      jcube = CEILING((beta  + piq) / da)
      IF (icube<1.OR.icube>ncube.OR.jcube<1.OR.jcube>ncube) THEN
        WRITE(*,*) "fatal error in search algorithm"
        WRITE(*,*) "icube or jcube out of range: ",icube,jcube
        STOP
      END IF
      wt    = SIN( lat(j)+0.5*dlat ) - SIN( lat(j)-0.5*dlat )
      weight(icube,jcube,ipanel) = weight(icube,jcube,ipanel)+wt
      !
      terr_cube    (icube,jcube,ipanel)     = terr_cube    (icube,jcube,ipanel)+wt*DBLE(terr(i,j))
      landfrac_cube(icube,jcube,ipanel)     = landfrac_cube(icube,jcube,ipanel)+wt*DBLE(landfrac(i,j))
      !
      ! save "index-association" for variance computation
      !
      idx(i,j) = icube
      idy(i,j) = jcube
      idp(i,j) = ipanel
    END DO
  END DO
  
  dx = deg2rad*(lon_landm(2)-lon_landm(1))
  !
  ! lat_landm is not exactly equally spaced so a search is needed in the loop below
  !
  dy = deg2rad*(lat_landm(2)-lat_landm(1))
  DO k=1,6
    DO j=1,ncube          
      DO i=1,ncube
        IF (ABS(weight(i,j,k))<1.0E-9) THEN
          WRITE(*,*) "there is no lat-lon grid point in cubed sphere cell ",i,j,k
          WRITE(*,*) "fatal error"
          STOP
        ELSE
          terr_cube        (i,j,k) = terr_cube        (i,j,k)/weight(i,j,k)                
          landfrac_cube    (i,j,k) = landfrac_cube    (i,j,k)/weight(i,j,k)                
        END IF
        !
        ! linearly interpolate landm_coslat
        !
        alpha = -piq+(i-0.5)*da
        beta  = -piq+(j-0.5)*da
        CALL CubedSphereRLLFromABP(alpha, beta, k, lambda, theta)   
        IF (theta>lat_landm(jm_landm)*deg2rad-tiny) THEN
          landm_coslat_cube(i,j,k) = 0.0
        ELSE IF (theta<lat_landm(1)*deg2rad+tiny) THEN
          landm_coslat_cube(i,j,k) = 1.0
        ELSE
          !
          ! this code assumes data is equally spaced in longitude
          !
          ilon = MAX(MIN(INT((lambda-lon_landm(1)*deg2rad)/dx)+1,im_landm),1)
          ip1  = MOD(ilon,im_landm)+1
          wx = (lambda-lon_landm(ilon)*deg2rad)/dx
          !
          ilat = MAX(MIN(INT((theta -lat_landm(1)*deg2rad)/dy)+1,jm_landm-1),1)
          jp1  = ilat+1
          wy = (theta -lat_landm(ilat)*deg2rad)/(lat_landm(jp1)-lat_landm(ilat))
          !
          ! since LANDM_COSLAT is not equally spaced in latitude a search is needed
          ! 
          DO WHILE (wy>1.0.OR.wy<0.0)
            jp1  = ilat+1
            wy = (theta -lat_landm(ilat)*deg2rad)/((lat_landm(jp1)-lat_landm(ilat))*deg2rad)
            IF (wy>1.0) THEN
              ilat=ilat+1
            ELSE IF (wy<0.0) THEN
              ilat=ilat-1
            END IF
          END DO
          
          IF (wx>1.0+tiny.OR.wx<0.0-tiny) THEN
            WRITE(*,*) "wx out of range",wx
            stop
          END IF
          IF (wy>1.0+tiny.OR.wy<0.0-tiny) THEN
            WRITE(*,*) "wy out of range",wy
            stop
          END IF
          !
          ! "crude" bi-linear interpolation
          !
          landm_coslat_cube(i,j,k) =&
               (1.0-wx)*(1.0-wy)*landm_coslat(ilon,ilat)+   wx *(1-wy)*landm_coslat(ip1,ilat)+&
               (1.0-wx)*     wy *landm_coslat(ilon,jp1 )+   wx * wy   *landm_coslat(ip1,jp1)
        END IF
      END DO
    END DO
  END DO
  WRITE(*,*) "min/max value of terr_cube:", MINVAL(terr_cube), MAXVAL(terr_cube)
  WRITE(*,*) "min/max value of landm_coslat_cube:", MINVAL(landm_coslat_cube), MAXVAL(landm_coslat_cube)
  !
  ! compute volume of topography on cubed-sphere
  !
  WRITE(*,*) "compute volume for cubed-sphere binned data"
  allocate (darea_cube(ncube,ncube),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for idp'
    stop
  end if
  CALL EquiangularAllAreas(ncube, darea_cube)
  vol_cube = 0.0
  do ipanel=1,6
    do j=1,ncube
      do i=1,ncube
        vol_cube = vol_cube+terr_cube(i,j,ipanel)*darea_cube(i,j)
      end do
    end do
  end do 
  vol_cube=vol_cube/(4.0*pi)
  deallocate(darea_cube)
  WRITE(*,*) "mean height (globally) of topography about sea-level (3km cube data)",vol_cube,(vol_cube-vol)/vol
  !*********************************************************
  !
  ! compute variance
  !
  !*********************************************************
  !
  allocate ( sgh30_cube(ncube,ncube,6),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for sgh30_cube'
    stop
  end if
  sgh30_cube = 0.0
  DO j=1,jm
    DO i=1,im
      icube  = idx(i,j)
      jcube  = idy(i,j)
      ipanel = idp(i,j)
      wt    = SIN( lat(j)+0.5*dlat ) - SIN( lat(j)-0.5*dlat )
      sgh30_cube(icube,jcube,ipanel) = sgh30_cube(icube,jcube,ipanel) + &
           (wt*(terr_cube(icube,jcube,ipanel)-terr(i,j))**2)/weight(icube,jcube,ipanel)
    END DO
  END DO
  !        sgh30_cube=sgh30_cube/weight
  WRITE(*,*) "min/max value of sgh30_cube:", MINVAL(sgh30_cube), MAXVAL(sgh30_cube)
  !
  ! write data to NetCDF file
  !
  CALL wrt_cube(ncube,terr_cube,landfrac_cube,landm_coslat_cube,sgh30_cube)
  DEALLOCATE(weight,terr,landfrac,idx,idy,idp,lat,lon)
  WRITE(*,*) "done writing cubed sphere data"
end program convterr


!************************************************************************
!!handle_err
!************************************************************************
!
!!ROUTINE:      handle_err
!!DESCRIPTION:  error handler
!--------------------------------------------------------------------------

subroutine handle_err(status)
  
  implicit         none
  
#     include          <netcdf.inc>
  
  integer          status
  
  if (status .ne. nf_noerr) then
    print *, nf_strerror(status)
    stop 'Stopped'
  endif
  
end subroutine handle_err


!------------------------------------------------------------------------------
! SUBROUTINE CubedSphereABPFromRLL
!
! Description:
!   Determine the (alpha,beta,panel) coordinate of a point on the sphere from
!   a given regular lat lon coordinate.
!
! Parameters:
!   lon - Coordinate longitude
!   lat - Coordinate latitude
!   alpha (OUT) - Alpha coordinate
!   beta (OUT) - Beta coordinate
!   ipanel (OUT) - Face panel
!------------------------------------------------------------------------------
SUBROUTINE CubedSphereABPFromRLL(lon, lat, alpha, beta, ipanel)
  use shr_kind_mod, only: r8 => shr_kind_r8
  IMPLICIT NONE
  
  REAL    (R8), INTENT(IN)  :: lon, lat
  REAL    (R8), INTENT(OUT) :: alpha, beta
  INTEGER, INTENT(OUT) :: ipanel
  REAL    (r8), PARAMETER :: pi   = 3.14159265358979323846264338327
  REAL    (r8), PARAMETER :: piq   = 0.25*pi
  REAL    (r8), PARAMETER :: rotate_cube = 0.0
  
  ! Local variables
  REAL    (R8) :: xx, yy, zz, pm
  REAL    (R8) :: sx, sy, sz
  INTEGER  :: ix, iy, iz
  
  ! Translate to (x,y,z) space
  xx = COS(lon-rotate_cube) * COS(lat)
  yy = SIN(lon-rotate_cube) * COS(lat)
  zz = SIN(lat)
  
  pm = MAX(ABS(xx), ABS(yy), ABS(zz))
  
  ! Check maximality of the x coordinate
  IF (pm == ABS(xx)) THEN
    IF (xx > 0) THEN; ix = 1; ELSE; ix = -1; ENDIF
  ELSE
    ix = 0
  ENDIF
  
  ! Check maximality of the y coordinate
  IF (pm == ABS(yy)) THEN
    IF (yy > 0) THEN; iy = 1; ELSE; iy = -1; ENDIF
  ELSE
    iy = 0
  ENDIF
      
  ! Check maximality of the z coordinate
  IF (pm == ABS(zz)) THEN
    IF (zz > 0) THEN; iz = 1; ELSE; iz = -1; ENDIF
  ELSE
    iz = 0
  ENDIF
  
  ! Panel assignments
  IF (iz  ==  1) THEN
    ipanel = 6; sx = yy; sy = -xx; sz = zz
    
  ELSEIF (iz  == -1) THEN
    ipanel = 5; sx = yy; sy = xx; sz = -zz
    
  ELSEIF ((ix == 1) .AND. (iy /= 1)) THEN
    ipanel = 1; sx = yy; sy = zz; sz = xx
    
  ELSEIF ((ix == -1) .AND. (iy /= -1)) THEN
    ipanel = 3; sx = -yy; sy = zz; sz = -xx
    
  ELSEIF ((iy == 1) .AND. (ix /= -1)) THEN
    ipanel = 2; sx = -xx; sy = zz; sz = yy
    
  ELSEIF ((iy == -1) .AND. (ix /=  1)) THEN
    ipanel = 4; sx = xx; sy = zz; sz = -yy
    
  ELSE
    WRITE(*,*) 'Fatal Error: CubedSphereABPFromRLL failed'
    WRITE(*,*) '(xx, yy, zz) = (', xx, ',', yy, ',', zz, ')'
    WRITE(*,*) 'pm =', pm, ' (ix, iy, iz) = (', ix, ',', iy, ',', iz, ')'
    STOP
  ENDIF
  
  ! Use panel information to calculate (alpha, beta) coords
  alpha = ATAN(sx / sz)
  beta = ATAN(sy / sz)
  
END SUBROUTINE CubedSphereABPFromRLL



!
! write netCDF file
! 
subroutine wrt_cube(ncube,terr_cube,landfrac_cube,landm_coslat_cube,sgh30_cube)
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
#     include         <netcdf.inc>
  
  !
  ! Dummy arguments
  !
  integer, intent(in) :: ncube
  real (r8), dimension(6*ncube*ncube), intent(in) :: terr_cube,landfrac_cube,sgh30_cube,landm_coslat_cube
  !
  ! Local variables
  !
  !-----------------------------------------------------------------------
  !
  !     grid coordinates and masks
  !
  !-----------------------------------------------------------------------
  
  real (r8), dimension(6*ncube*ncube) :: grid_center_lat  ! lat/lon coordinates for
  real (r8), dimension(6*ncube*ncube) :: grid_center_lon  ! each grid center in degrees
  
  integer  :: ncstat             ! general netCDF status variable
  integer  :: nc_grid_id         ! netCDF grid dataset id
  integer  :: nc_gridsize_id     ! netCDF grid size dim id
  integer  :: nc_gridrank_id     ! netCDF grid rank dim id
  integer  :: nc_griddims_id     ! netCDF grid dimension size id
  integer  :: nc_grdcntrlat_id   ! netCDF grid center lat id
  integer  :: nc_grdcntrlon_id   ! netCDF grid center lon id
  integer  :: nc_terr_id
  integer  :: nc_landfrac_id
  integer  :: nc_landm_coslat_id
  integer  :: nc_var_id
  
  
  integer, dimension(2) :: nc_dims2_id ! netCDF dim id array for 2-d arrays
  integer :: grid_dims
  
  character(18), parameter :: grid_file_out = 'USGS-topo-cube.nc'
  character(90), parameter :: grid_name = 'equi-angular gnomonic cubed sphere grid'
  
  character (len=32) :: fout       ! NetCDF output file
  integer            :: foutid     ! Output file id
  integer            :: lonid, lonvid
  integer            :: latid, latvid
  integer            :: status    ! return value for error control of netcdf routin
  integer            :: i,j,k
  character (len=8)  :: datestring
  
  integer  :: atm_add,n
  real(r8) :: xgno_ce,lon,ygno_ce,lat
  
  REAL    (r8), PARAMETER :: pi    = 3.14159265358979323846264338327
  REAL    (r8), PARAMETER :: piq   = 0.25*pi
  REAL    (r8), PARAMETER :: rad2deg   = 180.0/pi
  
  real(r8) :: da, a1,a2,a3,a4,dbg_area,max_size
  real(r8), dimension(2,2) :: ang
  real(r8) :: tmp_lon,min_lon,max_lon!,sum,lflag_value
  logical :: lflag
  
  grid_dims = 6*ncube*ncube
  
  dbg_area = 0.0
  
  da = pi / DBLE(2*ncube)        
  atm_add = 1
  do k=1,6
    do j=1,ncube
      ygno_ce = -piq + da * (DBLE(j-1)+0.5) !center of cell
      do i=1,ncube
        xgno_ce = -piq + da * (DBLE(i-1)+0.5)
        call CubedSphereRLLFromABP(xgno_ce, ygno_ce, k, lon, lat)
        grid_center_lon(atm_add  ) = lon*rad2deg                          
        grid_center_lat(atm_add  ) = lat*rad2deg                                    
        atm_add = atm_add+1
      end do
    end do
  end do
  
  WRITE(*,*) "Create NetCDF file for output"
  ncstat = nf_create (grid_file_out, NF_64BIT_OFFSET,nc_grid_id)
  call handle_err(ncstat)
  
  ncstat = nf_put_att_text (nc_grid_id, NF_GLOBAL, 'title',len_trim(grid_name), grid_name)
  call handle_err(ncstat)
  
  WRITE(*,*) "define grid size dimension"
  ncstat = nf_def_dim (nc_grid_id, 'grid_size', 6*ncube*ncube, nc_gridsize_id)
  call handle_err(ncstat)
  
  WRITE(*,*) "define grid rank dimension"
  ncstat = nf_def_dim (nc_grid_id, 'grid_rank', 1, nc_gridrank_id)
  call handle_err(ncstat)
  
  WRITE(*,*) "define grid dimension size array"
  ncstat = nf_def_var (nc_grid_id, 'grid_dims', NF_INT,1, nc_gridrank_id, nc_griddims_id)
  call handle_err(ncstat)
  
  WRITE(*,*) "define grid center latitude array"
  ncstat = nf_def_var (nc_grid_id, 'lat', NF_DOUBLE,1, nc_gridsize_id, nc_grdcntrlat_id)
  call handle_err(ncstat)        
  ncstat = nf_put_att_text (nc_grid_id, nc_grdcntrlat_id, 'units',13, 'degrees_north')
  call handle_err(ncstat)
  
  WRITE(*,*) "define grid center longitude array"
  ncstat = nf_def_var (nc_grid_id, 'lon', NF_DOUBLE,1, nc_gridsize_id, nc_grdcntrlon_id)
  call handle_err(ncstat)
  ncstat = nf_put_att_text (nc_grid_id, nc_grdcntrlon_id, 'units',12, 'degrees_east')
  call handle_err(ncstat)
  
  WRITE(*,*) "define terr_cube array"
  ncstat = nf_def_var (nc_grid_id, 'terr', NF_DOUBLE,1, nc_gridsize_id, nc_terr_id)
  call handle_err(ncstat)
  ncstat = nf_put_att_text (nc_grid_id, nc_terr_id, 'units',1, 'm')
  call handle_err(ncstat)
  
  WRITE(*,*) "define landfrac_cube array"
  ncstat = nf_def_var (nc_grid_id, 'LANDFRAC', NF_DOUBLE,1, nc_gridsize_id, nc_landfrac_id)
  call handle_err(ncstat)
  ncstat = nf_put_att_text (nc_grid_id, nc_landfrac_id, 'long_name',70,&
       'land ocean transition mask: ocean (0), continent (1), transition (0-1)')
  call handle_err(ncstat)
  
  WRITE(*,*) "define landm_coslat_cube array"
  ncstat = nf_def_var (nc_grid_id, 'LANDM_COSLAT', NF_DOUBLE,1, nc_gridsize_id, nc_landm_coslat_id)
  call handle_err(ncstat)
  ncstat = nf_put_att_text (nc_grid_id, nc_landm_coslat_id, 'long_name',35,'smoothed land ocean transition mask')
  call handle_err(ncstat)
  
  WRITE(*,*) "define sgh30_cube array"
  ncstat = nf_def_var (nc_grid_id, 'SGH30', NF_DOUBLE,1, nc_gridsize_id, nc_var_id)
  call handle_err(ncstat)
  ncstat = nf_put_att_text (nc_grid_id, nc_var_id, 'units',12, 'm')
  call handle_err(ncstat)
  ncstat = nf_put_att_text (nc_grid_id, nc_var_id, 'long_name',58,&
       'variance of elevation from 30s lat-lon to 3km cubed-sphere')
  
  WRITE(*,*) "end definition stage"
  ncstat = nf_enddef(nc_grid_id)
  call handle_err(ncstat)
  
  !-----------------------------------------------------------------------
  !
  !     write grid data
  !
  !-----------------------------------------------------------------------
  
  
  WRITE(*,*) "write grid data"        
  ncstat = nf_put_var_int(nc_grid_id, nc_griddims_id, grid_dims)
  call handle_err(ncstat)
  
  ncstat = nf_put_var_double(nc_grid_id, nc_grdcntrlat_id, grid_center_lat)
  call handle_err(ncstat)
  
  ncstat = nf_put_var_double(nc_grid_id, nc_grdcntrlon_id, grid_center_lon)
  call handle_err(ncstat)
  
  ncstat = nf_put_var_double(nc_grid_id, nc_terr_id, terr_cube)
  call handle_err(ncstat)
  
  ncstat = nf_put_var_double(nc_grid_id, nc_landfrac_id, landfrac_cube)
  call handle_err(ncstat)
  
  ncstat = nf_put_var_double(nc_grid_id, nc_landm_coslat_id, landm_coslat_cube)
  call handle_err(ncstat)
  
  ncstat = nf_put_var_double(nc_grid_id, nc_var_id, sgh30_cube)
  call handle_err(ncstat)
  
  WRITE(*,*) "Close output file"
  ncstat = nf_close(nc_grid_id)
  call handle_err(ncstat)
end subroutine wrt_cube


!------------------------------------------------------------------------------
! SUBROUTINE EquiangularAllAreas
!
! Description:
!   Compute the area of all cubed sphere grid cells, storing the results in
!   a two dimensional array.
!
! Parameters: 
!   icube - Resolution of the cubed sphere
!   dA (OUT) - Output array containing the area of all cubed sphere grid cells
!------------------------------------------------------------------------------
SUBROUTINE EquiangularAllAreas(icube, dA)
  use shr_kind_mod, only: r8 => shr_kind_r8        
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)                           :: icube
  REAL (r8), DIMENSION(icube,icube), INTENT(OUT) :: dA
  
  ! Local variables
  INTEGER                       :: k, k1, k2
  REAL (r8)                          :: a1, a2, a3, a4
  REAL (r8), DIMENSION(icube+1,icube+1)  :: ang
  REAL (r8), DIMENSION(icube+1)      :: gp
  
  REAL    (r8), PARAMETER :: pi   = 3.14159265358979323846264338327
  REAL    (r8), PARAMETER :: piq   = 0.25*pi
  
  
  !#ifdef DBG 
  REAL (r8)   :: dbg1 !DBG
  !#endif
  
  ! Recall that we are using equi-angular spherical gridding
  !   Compute the angle between equiangular cubed sphere projection grid lines.
  DO k = 1, icube+1
    gp(k) = -piq + (pi/DBLE(2*(icube))) * DBLE(k-1)
  ENDDO
  
  DO k2=1,icube+1
    DO k1=1,icube+1
      ang(k1,k2) =ACOS(-SIN(gp(k1)) * SIN(gp(k2)))
    ENDDO
  ENDDO
  
  DO k2=1,icube
    DO k1=1,icube
      a1 =      ang(k1  , k2  )
      a2 = pi - ang(k1+1, k2  )
      a3 = pi - ang(k1  , k2+1)
      a4 =      ang(k1+1, k2+1)
      
      ! area = r*r*(-2*pi+sum(interior angles))
      DA(k1,k2) = -2.0*pi+a1+a2+a3+a4
    ENDDO
  ENDDO
  
  !#ifdef DBG 
  ! Only for debugging - test consistency
  dbg1 = 0.0                           !DBG
  DO k2=1,icube
    DO k1=1,icube
      dbg1 = dbg1 + DA(k1,k2)         !DBG
    ENDDO
  ENDDO
  write(*,*) 'DAcube consistency: ',dbg1-4.0*pi/6.0 !DBG
  !#endif
END SUBROUTINE EquiangularAllAreas


!------------------------------------------------------------------------------
! SUBROUTINE CubedSphereRLLFromABP
!
! Description:
!   Determine the lat lon coordinate of a point on a sphere given its
!   (alpha,beta,panel) coordinate.
!
! Parameters:
!   alpha - Alpha coordinate
!   beta - Beta coordinate
!   panel - Cubed sphere panel id
!   lon (OUT) - Calculated longitude
!   lat (OUT) - Calculated latitude
!------------------------------------------------------------------------------
SUBROUTINE CubedSphereRLLFromABP(alpha, beta, ipanel, lon, lat)
  use shr_kind_mod, only: r8 => shr_kind_r8        
  IMPLICIT NONE        
  REAL    (r8), INTENT(IN)  :: alpha, beta
  INTEGER     , INTENT(IN)  :: ipanel
  REAL    (r8), INTENT(OUT) :: lon, lat        
  ! Local variables
  REAL    (r8) :: xx, yy, zz, rotate_cube
  REAL    (r8), PARAMETER :: pi   = 3.14159265358979323846264338327
  REAL    (r8), PARAMETER :: piq  = 0.25*pi
  
  rotate_cube = 0.0
  ! Convert to cartesian coordinates
  CALL CubedSphereXYZFromABP(alpha, beta, ipanel, xx, yy, zz)        
  ! Convert back to lat lon
  lat = ASIN(zz)
  if (xx==0.0.and.yy==0.0) THEN
    lon = 0.0
  else
    lon = ATAN2(yy, xx) +rotate_cube 
    IF (lon<0.0) lon=lon+2.0*pi
    IF (lon>2.0*pi) lon=lon-2.0*pi
  end if
END SUBROUTINE CubedSphereRLLFromABP

!------------------------------------------------------------------------------
! SUBROUTINE CubedSphereXYZFromABP
!
! Description:
!   Determine the Cartesian coordinate of a point on a sphere given its
!   (alpha,beta,panel) coordinate.
!
! Parameters:
!   alpha - Alpha coordinate
!   beta - Beta coordinate
!   panel - Cubed sphere panel id
!   xx (OUT) - Calculated x coordinate
!   yy (OUT) - Calculated y coordinate
!   zz (OUT) - Calculated z coordinate
!------------------------------------------------------------------------------
SUBROUTINE CubedSphereXYZFromABP(alpha, beta, ipanel, xx, yy, zz)
  use shr_kind_mod, only: r8 => shr_kind_r8        
  IMPLICIT NONE
  
  REAL    (r8), INTENT(IN)  :: alpha, beta
  INTEGER     , INTENT(IN)  :: ipanel
  REAL    (r8), INTENT(OUT) :: xx, yy, zz        
  ! Local variables
  REAL    (r8) :: a1, b1, pm
  REAL    (r8) :: sx, sy, sz       
  
  ! Convert to Cartesian coordinates
  a1 = TAN(alpha)
  b1 = TAN(beta)
  
  sz = (1.0 + a1 * a1 + b1 * b1)**(-0.5)
  sx = sz * a1
  sy = sz * b1        
  ! Panel assignments
  IF (ipanel == 6) THEN
    yy = sx; xx = -sy; zz = sz          
  ELSEIF (ipanel == 5) THEN
    yy = sx; xx = sy; zz = -sz          
  ELSEIF (ipanel == 1) THEN
    yy = sx; zz = sy; xx = sz          
  ELSEIF (ipanel == 3) THEN
    yy = -sx; zz = sy; xx = -sz          
  ELSEIF (ipanel == 2) THEN
    xx = -sx; zz = sy; yy = sz          
  ELSEIF (ipanel == 4) THEN
    xx = sx; zz = sy; yy = -sz          
  ELSE
    WRITE(*,*) 'Fatal Error: Panel out of range in CubedSphereXYZFromABP'
    WRITE(*,*) '(alpha, beta, panel) = (', alpha, ',', beta, ',', ipanel, ')'
    STOP
  ENDIF
END SUBROUTINE CubedSphereXYZFromABP


