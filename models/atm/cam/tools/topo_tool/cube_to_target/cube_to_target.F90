!
!  DATE CODED:   Nov 7, 2011 to Oct 15, 2012
!  DESCRIPTION:  Remap topo data from cubed-sphere grid to target grid using rigorous remapping
!                (Lauritzen, Nair and Ullrich, 2010, J. Comput. Phys.)
!
!  Author: Peter Hjort Lauritzen (pel@ucar.edu), AMP/CGD/NESL/NCAR 
!
program convterr
  use shr_kind_mod, only: r8 => shr_kind_r8
  use reconstruct
  implicit none
#     include         <netcdf.inc>

  !**************************************
  !
  ! USER SETTINGS BELOW
  !
  !**************************************
  !
  !
  ! if smoothed PHIS is available SGH needs to be recomputed  to account for the sub-grid-scale
  ! variability introduced by the smoothing
  !
  logical :: lsmooth_terr = .FALSE. 
  !
  ! PHIS is smoothed by other software/dynamical core
  !
  logical :: lexternal_smooth_terr = .FALSE. ! lexternal_smooth_terr = .FALSE. is NOT supported currently
  !
  ! set PHIS=0.0 if LANDFRAC<0.01
  !
  logical :: lzero_out_ocean_point_phis = .FALSE.
  !
  ! For internal smoothing (experimental at this point)
  ! ===================================================
  !
  ! if smoothing is internal (lexternal_smooth_terr=.FALSE.) choose coarsening factor
  !
  ! recommendation: 2*(target resolution)/(0.03 degree)
  !
  ! factor must be an even integer
  !
  integer, parameter :: factor = 60 !coarse grid = 2.25 degrees
  integer, parameter :: norder = 2
  integer, parameter :: nmono  = 0
  integer, parameter :: npd    = 1
  !
  !**********************************************************************
  !
  ! END OF USER SETTINS BELOW
  ! (do not edit beyond this point unless you know what you are doing!)
  !
  !**********************************************************************
  !
  integer :: im, jm, ncoarse
  integer :: ncube !dimension of cubed-sphere grid
  
  real(r8),  allocatable, dimension(:) :: landm_coslat, landfrac, terr, sgh30
  real(r8),  allocatable, dimension(:) :: terr_coarse !for internal smoothing
  
  integer :: alloc_error,dealloc_error
  integer :: i,j,n,k,index                               
  integer*2, allocatable, dimension(:,:)  :: iterr        ! terrain data for 30-sec tile
  integer ncid,status, dimlatid,dimlonid, landid, topoid  ! for netCDF USGS data file
  integer :: srcid,dstid,  jm_dbg  ! for netCDF weight file
  integer, dimension(2) ::  src_grid_dims  ! for netCDF weight file
  
  integer :: dimid
  
  logical :: ldbg
  real(r8), allocatable, dimension(:)   :: lon  , lat
  real(r8), allocatable, dimension(:)   :: lon_landm  , lat_landm
  real(r8), allocatable, dimension(:)   :: area
  integer :: im_landm, jm_landm
  integer :: lonid, latid, phisid
  !
  ! constants
  !  
  REAL    (r8), PARAMETER :: pi        = 3.14159265358979323846264338327
  REAL    (r8), PARAMETER :: piq       = 0.25*pi
  REAL    (r8), PARAMETER :: pih       = 0.50*pi
  REAL    (r8), PARAMETER :: deg2rad   = pi/180.0
  
  real(r8) :: wt,dlat
  integer  :: ipanel,icube,jcube
  real(r8), allocatable, dimension(:,:,:)   :: weight,terr_cube,landfrac_cube,sgh30_cube
  real(r8), allocatable, dimension(:,:,:)   :: landm_coslat_cube
  integer, allocatable, dimension(:,:) :: idx,idy,idp
  integer :: npatch, isub,jsub, itmp, iplm1,jmin,jmax
  real(r8) :: sum,dx,scale,dmax,arad,jof,term,s1,c1,clon,iof,dy,s2,c2,dist
  !
  ! for linear interpolation
  !
  real(r8) :: lambda,theta,wx,wy,offset
  integer :: ilon,ilat,ip1,jp1
  !
  ! variable for regridding
  !
  integer :: src_grid_dim  ! for netCDF weight file
  integer :: n_a,n_b,n_s,n_aid,n_bid,n_sid
  integer :: count
  real(r8), allocatable, dimension(:) :: landfrac_target, terr_target, sgh30_target, sgh_target
  real(r8), allocatable, dimension(:) :: landm_coslat_target, area_target
  !
  ! this is only used if target grid is a lat-lon grid
  !
  integer , parameter :: im_target = 360 , jm_target = 180
  !
  ! this is only used if target grid is not a lat-lon grid
  !
  real(r8), allocatable, dimension(:) :: lon_target, lat_target
  !
  ! new
  !
  integer :: ntarget, ntarget_id, ncorner, ncorner_id, nrank, nrank_id
  integer :: ntarget_smooth
  real(r8), allocatable, dimension(:,:):: target_corner_lon, target_corner_lat
  real(r8), allocatable, dimension(:)  :: target_center_lon, target_center_lat, target_area
  integer :: ii,ip,jx,jy,jp
  real(r8), dimension(:), allocatable  :: xcell, ycell, xgno, ygno
  real(r8), dimension(:), allocatable  :: gauss_weights,abscissae
  integer, parameter :: ngauss = 3
  integer :: jmax_segments,jall
  real(r8) :: tmp
  
  real(r8), allocatable, dimension(:,:) :: weights_all
  integer , allocatable, dimension(:,:) :: weights_eul_index_all
  integer , allocatable, dimension(:)   :: weights_lgr_index_all
  integer :: ix,iy
  !
  ! volume of topography
  !
  real(r8) :: vol_target, vol_target_un, area_target_total,vol_source,vol_tmp
  integer :: nlon,nlon_smooth,nlat,nlat_smooth
  logical :: ltarget_latlon,lpole
  real(r8), allocatable, dimension(:,:)   :: terr_smooth
  !
  ! for internal filtering
  !
  real(r8), allocatable, dimension(:,:) :: weights_all_coarse
  integer , allocatable, dimension(:,:) :: weights_eul_index_all_coarse
  integer , allocatable, dimension(:)   :: weights_lgr_index_all_coarse
  real(r8), allocatable, dimension(:)   :: area_target_coarse
  real(r8), allocatable, dimension(:,:) :: da_coarse,da
  real(r8), allocatable, dimension(:,:) :: recons,centroids
  integer :: nreconstruction
  
  integer :: jmax_segments_coarse,jall_coarse,ncube_coarse
  real(r8) :: all_weights

  !
  ! turn extra debugging on/off
  !
  ldbg = .FALSE.
  
  nreconstruction = 1
  !
  !*********************************************************
  !
  ! read in target grid
  !
  !*********************************************************
  !
  status = nf_open('target.nc', 0, ncid)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
  
  status = NF_INQ_DIMID(ncid, 'grid_size', ntarget_id)
  status = NF_INQ_DIMLEN(ncid, ntarget_id, ntarget)
  WRITE(*,*) "dimension of target grid: ntarget=",ntarget
  
  status = NF_INQ_DIMID(ncid, 'grid_corners', ncorner_id)
  status = NF_INQ_DIMLEN(ncid, ncorner_id, ncorner)
  WRITE(*,*) "maximum number of corners: ncorner=",ncorner
  
  status = NF_INQ_DIMID(ncid, 'grid_rank', nrank_id);status = NF_INQ_DIMLEN(ncid, nrank_id, nrank)
  WRITE(*,*) "grid rank: nrank=",nrank
  IF (nrank==2) THEN
    WRITE(*,*) "target grid is a lat-lon grid"
    ltarget_latlon = .TRUE.
    status = NF_INQ_DIMID(ncid, 'nlon', ntarget_id)
    status = NF_INQ_DIMLEN(ncid, ntarget_id, nlon)
    status = NF_INQ_DIMID(ncid, 'nlat', ntarget_id)
    status = NF_INQ_DIMLEN(ncid, ntarget_id, nlat)
    status = NF_INQ_DIMID(ncid, 'lpole', ntarget_id)
    status = NF_INQ_DIMLEN(ncid, ntarget_id, lpole)
    WRITE(*,*) "nlon=",nlon,"nlat=",nlat
    IF (lpole) THEN
      WRITE(*,*) "center of most Northern grid cell is lat=90; similarly for South pole"
    ELSE
      WRITE(*,*) "center of most Northern grid cell is NOT lat=90; similarly for South pole"
    END IF
  ELSE IF (nrank==1) THEN
    ltarget_latlon = .FALSE.
  ELSE
    WRITE(*,*) "nrank out of range",nrank
    STOP
  ENDIF
  
  allocate ( target_corner_lon(ncorner,ntarget),stat=alloc_error)
  allocate ( target_corner_lat(ncorner,ntarget),stat=alloc_error)
  
  status = NF_INQ_VARID(ncid, 'grid_corner_lon', lonid)
  status = NF_GET_VAR_DOUBLE(ncid, lonid,target_corner_lon)
  IF (maxval(target_corner_lon)>10.0) target_corner_lon = deg2rad*target_corner_lon
  
  status = NF_INQ_VARID(ncid, 'grid_corner_lat', latid)
  status = NF_GET_VAR_DOUBLE(ncid, latid,target_corner_lat)
  IF (maxval(target_corner_lat)>10.0) target_corner_lat = deg2rad*target_corner_lat
  !
  ! for writing remapped data on file at the end of the program
  !
  allocate ( target_center_lon(ntarget),stat=alloc_error)
  allocate ( target_center_lat(ntarget),stat=alloc_error)
  allocate ( target_area      (ntarget),stat=alloc_error)!dbg
  
  status = NF_INQ_VARID(ncid, 'grid_center_lon', lonid)
  status = NF_GET_VAR_DOUBLE(ncid, lonid,target_center_lon)
  
  status = NF_INQ_VARID(ncid, 'grid_center_lat', latid)
  status = NF_GET_VAR_DOUBLE(ncid, latid,target_center_lat)
  
  status = NF_INQ_VARID(ncid, 'grid_area', latid)
  status = NF_GET_VAR_DOUBLE(ncid, latid,target_area)
  
  status = nf_close (ncid)
  if (status .ne. NF_NOERR) call handle_err(status)          
  !
  !****************************************************
  !
  ! get dimension of cubed-sphere grid
  !
  !****************************************************
  !
  WRITE(*,*) "get dimension of cubed-sphere data from file"
  status = nf_open('USGS-topo-cube.nc', 0, ncid)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
  
  status = NF_INQ_DIMID(ncid, 'grid_size', dimid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_INQ_DIMLEN(ncid, dimid, n)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  ncube = INT(SQRT(DBLE(n/6)))
  WRITE(*,*) "cubed-sphere dimension: ncube = ",ncube
  WRITE(*,*) "average grid-spacing at the Equator (degrees):" ,90.0/ncube
  
  status = nf_close (ncid)
  if (status .ne. NF_NOERR) call handle_err(status)          
  !
  !****************************************************
  !
  ! compute weights for remapping
  !
  !****************************************************
  !
  jall = ncube*ncube*12*10 !anticipated number of weights (cab be tweaked)
  jmax_segments = 100000   !can be tweaked
  
  allocate (weights_all(jall,nreconstruction),stat=alloc_error )
  allocate (weights_eul_index_all(jall,3),stat=alloc_error )
  allocate (weights_lgr_index_all(jall),stat=alloc_error )
  
  CALL overlap_weights(weights_lgr_index_all,weights_eul_index_all,weights_all,&
       jall,ncube,ngauss,ntarget,ncorner,jmax_segments,target_corner_lon,target_corner_lat,nreconstruction)
  !
  !****************************************************
  !
  ! read cubed-sphere 3km data
  !
  !****************************************************
  !
  WRITE(*,*) "read cubed-sphere 3km data from file"
  status = nf_open('USGS-topo-cube.nc', 0, ncid)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
  
  status = NF_INQ_DIMID(ncid, 'grid_size', dimid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_INQ_DIMLEN(ncid, dimid, n)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  ncube = INT(SQRT(DBLE(n/6)))
  WRITE(*,*) "cubed-sphere dimension, ncube: ",ncube
  
  allocate ( landm_coslat(n),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for landfrac'
    stop
  end if
  
  status = NF_INQ_VARID(ncid, 'LANDM_COSLAT', landid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

  status = NF_GET_VAR_DOUBLE(ncid, landid,landm_coslat)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  WRITE(*,*) "min/max of landm_coslat",MINVAL(landm_coslat),MAXVAL(landm_coslat)
  !
  ! read LANDFRAC
  !
  allocate ( landfrac(n),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for landfrac'
    stop
  end if
  
  status = NF_INQ_VARID(ncid, 'LANDFRAC', landid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  status = NF_GET_VAR_DOUBLE(ncid, landid,landfrac)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  WRITE(*,*) "min/max of landfrac",MINVAL(landfrac),MAXVAL(landfrac)
  !
  ! read terr
  !
  allocate ( terr(n),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for landfrac'
    stop
  end if
  
  status = NF_INQ_VARID(ncid, 'terr', landid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  status = NF_GET_VAR_DOUBLE(ncid, landid,terr)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  WRITE(*,*) "min/max of terr",MINVAL(terr),MAXVAL(terr)
  !
  !
  !
  allocate ( sgh30(n),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for landfrac'
    stop
  end if
  
  status = NF_INQ_VARID(ncid, 'SGH30', landid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  
  status = NF_GET_VAR_DOUBLE(ncid, landid,sgh30)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  WRITE(*,*) "min/max of sgh30",MINVAL(sgh30),MAXVAL(sgh30)
  print *,"close file"
  status = nf_close (ncid)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  WRITE(*,*) 'done reading in LANDM_COSLAT data from netCDF file'
  !
  !*********************************************************
  !
  ! do actual remapping
  !
  !*********************************************************
  !
  allocate (terr_target(ntarget),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for terr_target'
    stop
  end if
  allocate (landfrac_target(ntarget),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for landfrac_target'
    stop
  end if
  allocate (landm_coslat_target(ntarget),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for landfrac_target'
    stop
  end if
  allocate (sgh30_target(ntarget),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for sgh30_target'
    stop
  end if
  allocate (area_target(ntarget),stat=alloc_error )
  terr_target     = 0.0
  landfrac_target = 0.0
  sgh30_target    = 0.0
  landm_coslat_target = 0.0
  area_target = 0.0
  
  tmp = 0.0
  do count=1,jall
    i    = weights_lgr_index_all(count)
    wt = weights_all(count,1)
    area_target        (i) = area_target(i) + wt
  end do
  
  do count=1,jall
    i    = weights_lgr_index_all(count)
    
    ix  = weights_eul_index_all(count,1)
    iy  = weights_eul_index_all(count,2)
    ip  = weights_eul_index_all(count,3)
    !
    ! convert to 1D indexing of cubed-sphere
    !
    ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix
    
    wt = weights_all(count,1)
    
    terr_target        (i) = terr_target        (i) + wt*terr        (ii)/area_target(i)
    landfrac_target    (i) = landfrac_target    (i) + wt*landfrac    (ii)/area_target(i)
    landm_coslat_target(i) = landm_coslat_target(i) + wt*landm_coslat(ii)/area_target(i)
    sgh30_target       (i) = sgh30_target       (i) + wt*sgh30       (ii)/area_target(i)
    
    tmp = tmp+wt*terr(ii)
  end do
  
  
  write(*,*) "tmp", tmp
  WRITE(*,*) "max difference between target grid area and remapping software area",&
              MAXVAL(target_area-area_target)
  
  do count=1,ntarget
    if (terr_target(count)>8848.0) then
      !
      ! max height is higher than Mount Everest
      !
      write(*,*) "FATAL error: max height is higher than Mount Everest!"
      write(*,*) "terr_target",count,terr_target(count)
      write(*,*) "(lon,lat) locations of vertices of cell with excessive max height::"
      do i=1,ncorner
        write(*,*) target_corner_lon(i,count),target_corner_lat(i,count)
      end do
      STOP
    else if (terr_target(count)<-423.0) then
      !
      ! min height is lower than Dead Sea
      !
      write(*,*) "FATAL error: min height is lower than Dead Sea!"
      write(*,*) "terr_target",count,terr_target(count)
      write(*,*) "(lon,lat) locations of vertices of cell with excessive min height::"
      do i=1,ncorner
        write(*,*) target_corner_lon(i,count),target_corner_lat(i,count)
      end do
      STOP
    else 
      
    end if
  end do
  WRITE(*,*) "Elevation data passed min/max consistency check!"
  WRITE(*,*) 
  
  WRITE(*,*) "min/max of unsmoothed terr_target        : ",MINVAL(terr_target    ),MAXVAL(terr_target    )
  WRITE(*,*) "min/max of landfrac_target               : ",MINVAL(landfrac_target),MAXVAL(landfrac_target)
  WRITE(*,*) "min/max of landm_coslat_target           : ",&
       MINVAL(landm_coslat_target),MAXVAL(landm_coslat_target)
  WRITE(*,*) "min/max of var30_target                  : ",MINVAL(sgh30_target   ),MAXVAL(sgh30_target   )
  !
  ! compute mean height (globally) of topography about sea-level for target grid unfiltered elevation
  !
  vol_target_un     = 0.0
  area_target_total = 0.0
  DO i=1,ntarget
    area_target_total = area_target_total+area_target(i)
    vol_target_un = vol_target_un+terr_target(i)*area_target(i)
  END DO
  WRITE(*,*) "mean height (globally) of topography about sea-level for target grid unfiltered elevation",&
       vol_target_un/area_target_total
  
  !
  ! diagnostics
  !
  vol_source     = 0.0
  allocate ( dA(ncube,ncube),stat=alloc_error )
  CALL EquiangularAllAreas(ncube, dA)
  DO jp=1,6
    DO jy=1,ncube
      DO jx=1,ncube
        ii = (jp-1)*ncube*ncube+(jy-1)*ncube+jx
        vol_source = vol_source+terr(ii)*dA(jx,jy)
      END DO
    END DO
  END DO
  WRITE(*,*) "volume of input cubed-sphere terrain           :",vol_source
  WRITE(*,*) "average elevation of input cubed-sphere terrain:",vol_source/(4.0*pi)
  
  DEALLOCATE(dA)
  !
  !
  !
  allocate (sgh_target(ntarget),stat=alloc_error )
  if( alloc_error /= 0 ) then
    print*,'Program could not allocate space for sgh_target'
    stop
  end if
  !
  ! compute variance with respect to cubed-sphere data
  !
  WRITE(*,*) "compute variance with respect to 3km cubed-sphere data: SGH"
  
  IF (lsmooth_terr) THEN
    WRITE(*,*) "smoothing PHIS"
    IF (lexternal_smooth_terr) THEN
      WRITE(*,*) "using externally generated smoothed topography"
      
      status = nf_open('phis-smooth.nc', 0, ncid)
      IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)           
      !
      IF (.NOT.ltarget_latlon) THEN
        !
        !*********************************************************
        !
        ! read in smoothed topography
        !
        !*********************************************************
        !
        status = NF_INQ_DIMID (ncid, 'ncol', ntarget_id    )
        status = NF_INQ_DIMLEN(ncid, ntarget_id , ntarget_smooth)
        IF (ntarget.NE.ntarget_smooth) THEN
          WRITE(*,*) "mismatch in smoothed data-set and target grid specification"
          WRITE(*,*) ntarget, ntarget_smooth
          STOP
        END IF
        status = NF_INQ_VARID(ncid, 'PHIS', phisid)
        !
        ! overwrite terr_target with smoothed version
        !
        status = NF_GET_VAR_DOUBLE(ncid, phisid,terr_target)
        terr_target = terr_target/9.80616
      ELSE
        !
        ! read in smoothed lat-lon topography
        !
        status = NF_INQ_DIMID(ncid, 'lon', ntarget_id)
        status = NF_INQ_DIMLEN(ncid, ntarget_id, nlon_smooth)
        status = NF_INQ_DIMID(ncid, 'lat', ntarget_id)
        status = NF_INQ_DIMLEN(ncid, ntarget_id, nlat_smooth)
        IF (nlon.NE.nlon_smooth.OR.nlat.NE.nlat_smooth) THEN
          WRITE(*,*) "smoothed topography dimensions do not match target grid dimensions"
          WRITE(*,*) "target grid  : nlon       ,nlat        =",nlon,nlat
          WRITE(*,*) "smoothed topo: nlon_smooth,nlat_smooth =",nlon_smooth,nlat_smooth
          STOP
        END IF
        ALLOCATE(terr_smooth(nlon_smooth,nlat_smooth),stat=alloc_error)
        status = NF_INQ_VARID(ncid, 'PHIS', phisid)
        status = NF_GET_VAR_DOUBLE(ncid, phisid,terr_smooth)
        !
        ! overwrite terr_target with smoothed version
        !
        ii=1
        DO j=1,nlat
          DO i=1,nlon
            terr_target(ii) = terr_smooth(i,j)/9.80616                  
            ii=ii+1
          END DO
        END DO
        DEALLOCATE(terr_smooth)
      END IF
    ELSE
      WRITE(*,*) "unstested software - uncomment this line of you know what you are doing!"
      STOP
      !
      !*****************************************************
      !
      ! smoothing topography internally
      !
      !*****************************************************
      !
      WRITE(*,*) "internally smoothing orography"
      !            CALL smooth(terr_target,ntarget,target_corner_lon,target_corner_lat)
      !
      ! smooth topography internally
      !            
      ncoarse = n/(factor*factor)
      !
      ! 
      !
      ncube_coarse = ncube/factor
      WRITE(*,*) "resolution of coarse grid", 90.0/ncube_coarse
      allocate ( terr_coarse(ncoarse),stat=alloc_error )
      if( alloc_error /= 0 ) then
        print*,'Program could not allocate space for landfrac'
        stop
      end if
      WRITE(*,*) "coarsening"
      allocate ( dA_coarse(ncube_coarse,ncube_coarse),stat=alloc_error )
      CALL coarsen(terr,terr_coarse,factor,n,dA_coarse)
      !
      !
      !
      vol_tmp     = 0.0
      DO jp=1,6
        DO jy=1,ncube_coarse
          DO jx=1,ncube_coarse
            ii = (jp-1)*ncube_coarse*ncube_coarse+(jy-1)*ncube_coarse+jx
            vol_tmp = vol_tmp+terr_coarse(ii)*dA_coarse(jx,jy)
          END DO
        END DO
      END DO
      WRITE(*,*) "volume of coarsened cubed-sphere terrain           :",vol_source
      WRITE(*,*) "difference between coarsened cubed-sphere data and input cubed-sphere data",&
           vol_tmp-vol_source
      
      
      
      WRITE(*,*) "done coarsening"
      
      nreconstruction = 1
      IF (norder>1) THEN
        IF (norder == 2) THEN
          nreconstruction = 3
        ELSEIF (norder == 3) THEN
          nreconstruction = 6
        END IF
        ALLOCATE(recons   (nreconstruction, ncoarse), STAT=status)
        ALLOCATE(centroids(nreconstruction, ncoarse), STAT=status)
        CALL get_reconstruction(terr_coarse,norder, nmono, recons, npd,da_coarse,&
             ncube_coarse+1,nreconstruction,centroids)
        SELECT CASE (nmono) 
        CASE (0)
          WRITE(*,*) "coarse grid reconstructions are not filtered with shape-preesrving filter"
        CASE (1)
          WRITE(*,*) "coarse grid reconstructions are filtered with shape-preserving filter"
        CASE DEFAULT
          WRITE(*,*) "nmono out of range: ",nmono
          STOP
        END SELECT
        SELECT CASE (0)
        CASE (0)
          WRITE(*,*) "coarse grid reconstructions are not filtered with positive definite filter"
        CASE (1)
          WRITE(*,*) "coarse grid reconstructions filtered with positive definite filter"
        CASE DEFAULT
          WRITE(*,*) "npd out of range: ",npd
          STOP
        END SELECT
      END IF

      jall_coarse = (ncube*ncube*12) !anticipated number of weights
      jmax_segments_coarse = jmax_segments!/factor !
      WRITE(*,*) "anticipated",jall_coarse
      allocate (weights_all_coarse(jall_coarse,nreconstruction),stat=alloc_error )
      allocate (weights_eul_index_all_coarse(jall_coarse,3),stat=alloc_error )
      allocate (weights_lgr_index_all_coarse(jall_coarse),stat=alloc_error )
      !
      !
      !
      CALL overlap_weights(weights_lgr_index_all_coarse,weights_eul_index_all_coarse,weights_all_coarse,&
           jall_coarse,ncube_coarse,ngauss,ntarget,ncorner,jmax_segments_coarse,target_corner_lon,&
           target_corner_lat,nreconstruction)            
      WRITE(*,*) "MIN/MAX of area-weight [0:1]: ",&
           MINVAL(weights_all_coarse(:,1)),MAXVAL(weights_all_coarse(:,1))
      !
      ! compute new weights
      !
      
      ! 
      ! do mapping
      !
      terr_target = 0.0
      tmp = 0.0
      allocate ( area_target_coarse(ntarget),stat=alloc_error)
      all_weights = 0.0
      area_target_coarse = 0.0
      do count=1,jall_coarse
        i    = weights_lgr_index_all_coarse(count)
        wt = weights_all_coarse(count,1)
        area_target_coarse        (i) = area_target_coarse(i) + wt
        all_weights = all_weights+wt
      end do
      WRITE(*,*) "sum of all weights (coarse to target) minus area of sphere : ",all_weights-4.0*pi
      WRITE(*,*) "MIN/MAX of area_target_coarse [0:1]:",&
           MINVAL(area_target_coarse),MAXVAL(area_target_coarse)
      IF (norder==1) THEN
        do count=1,jall_coarse
          i    = weights_lgr_index_all_coarse(count)
          
          ix  = weights_eul_index_all_coarse(count,1)
          iy  = weights_eul_index_all_coarse(count,2)
          ip  = weights_eul_index_all_coarse(count,3)
          !
          ! convert to 1D indexing of cubed-sphere
          !
          ii = (ip-1)*ncube_coarse*ncube_coarse+(iy-1)*ncube_coarse+ix
          
          wt = weights_all_coarse(count,1)
          
          terr_target(i) = terr_target(i) + wt*terr_coarse(ii)/area_target_coarse(i)
          tmp = tmp+wt*terr_coarse(ii)
        end do
      ELSE IF (norder==2) THEN
        do count=1,jall_coarse
          i    = weights_lgr_index_all_coarse(count)
          IF (i>jall_coarse.OR.i<1) THEN
            WRITE(*,*) i,jall_coarse
            STOP
          END IF
          ix  = weights_eul_index_all_coarse(count,1)
          iy  = weights_eul_index_all_coarse(count,2)
          ip  = weights_eul_index_all_coarse(count,3)
          !
          ! convert to 1D indexing of cubed-sphere
          !
          ii = (ip-1)*ncube_coarse*ncube_coarse+(iy-1)*ncube_coarse+ix
          
          terr_target(i) = terr_target(i) + (weights_all_coarse(count,1)*(&
               !
               ! all constant terms 
               !
               terr_coarse(ii) &
               - recons(1,ii)*centroids(1,ii) &
               - recons(2,ii)*centroids(2,ii) &
               !
               !                     + recons(3,ii)*(2.0*centroids(1,ii)**2-centroids(3,ii))&
               !                     + recons(4,ii)*(2.0*centroids(2,ii)**2-centroids(4,ii))&
               !
               !                     + recons(5,ii)*(2.0*centroids(1,ii)*centroids(2,ii)-centroids(5,ii))&
               )+&
               !
               ! linear terms
               !
               weights_all_coarse(count,2)*(&
               
               recons(1,ii)&
               
               !                     - recons(3,ii)*2.0*centroids(1,ii)&
               !                     - recons(5,ii)*    centroids(2,ii)&
               )+&
               !
               weights_all_coarse(count,3)*(&
               recons(2,ii)&
               !
               !                     - recons(4,ii)*2.0*centroids(2,ii)&
               !                     - recons(5,ii)*    centroids(1,ii)&
               )&
               !
               ! quadratic terms
               !
               !                     weights_all_coarse(count,4)*recons(3,ii)+&
               !                     weights_all_coarse(count,5)*recons(4,ii)+&
               !                     weights_all_coarse(count,6)*recons(5,ii)
               )/area_target_coarse(i)
        end do
        DEALLOCATE(centroids)
        DEALLOCATE(recons)
        DEALLOCATE(weights_all_coarse)
        
      ELSE IF (norder==3) THEN
        !              recons(4,:) = 0.0
        !              recons(5,:) = 0.0
        do count=1,jall_coarse
          i    = weights_lgr_index_all_coarse(count)
          IF (i>jall_coarse.OR.i<1) THEN
            WRITE(*,*) i,jall_coarse
            STOP
          END IF
          ix  = weights_eul_index_all_coarse(count,1)
          iy  = weights_eul_index_all_coarse(count,2)
          ip  = weights_eul_index_all_coarse(count,3)
          !
          ! convert to 1D indexing of cubed-sphere
          !
          ii = (ip-1)*ncube_coarse*ncube_coarse+(iy-1)*ncube_coarse+ix
          
          !                terr_target(i) = terr_target(i) + wt*terr_coarse(ii)/area_target_coarse(i)
          
          !                WRITE(*,*) count,area_target_coarse(i)
          !                terr_target(i) = terr_target(i) + area_target_coarse(i)
          !
          terr_target(i) = terr_target(i) + (weights_all_coarse(count,1)*(&
               
               
               !                     centroids(5,ii))/area_target_coarse(i))
               !                     centroids(1,ii)/area_target_coarse(i))
               !                     /area_target_coarse(i))
               
               
               
               
               !
               ! all constant terms 
               !
               terr_coarse(ii) &
               - recons(1,ii)*centroids(1,ii) &
               - recons(2,ii)*centroids(2,ii) &
               !
               + recons(3,ii)*(2.0*centroids(1,ii)**2-centroids(3,ii))&
               + recons(4,ii)*(2.0*centroids(2,ii)**2-centroids(4,ii))&
               !
               + recons(5,ii)*(2.0*centroids(1,ii)*centroids(2,ii)-centroids(5,ii))&
               )+&
               !
               ! linear terms
               !
               weights_all_coarse(count,2)*(&
               
               recons(1,ii)&
               
               - recons(3,ii)*2.0*centroids(1,ii)&
               - recons(5,ii)*    centroids(2,ii)&
               )+&
               !
               weights_all_coarse(count,3)*(&
               recons(2,ii)&
               !
               - recons(4,ii)*2.0*centroids(2,ii)&
               - recons(5,ii)*    centroids(1,ii)&
               )+&
               !
               ! quadratic terms
               !
               weights_all_coarse(count,4)*recons(3,ii)+&
               weights_all_coarse(count,5)*recons(4,ii)+&
               weights_all_coarse(count,6)*recons(5,ii))/area_target_coarse(i)
        end do
        DEALLOCATE(centroids)
        DEALLOCATE(recons)
        DEALLOCATE(weights_all_coarse)
      END IF
      DEALLOCATE(area_target_coarse)
      WRITE(*,*) "done smoothing"
    END IF
    !
    ! compute mean height (globally) of topography about sea-level for target grid filtered elevation
    !
    vol_target = 0.0
    DO i=1,ntarget
      vol_target = vol_target+terr_target(i)*area_target(i)
      !            if (ABS(area_target(i)-area_target_coarse(i))>0.000001) THEN
      !              WRITE(*,*) "xxx",area_target(i),area_target_coarse(i),area_target(i)-area_target_coarse(i)
      !              STOP
      !            END IF
    END DO
    WRITE(*,*) "mean height (globally) of topography about sea-level for target grid filtered elevation",&
         vol_target/area_target_total
    WRITE(*,*) "percentage change in mean height between filtered and unfiltered elevations",&
         100.0*(vol_target-vol_target_un)/vol_target_un
    WRITE(*,*) "percentage change in mean height between input cubed-sphere and unfiltered elevations",&
         100.0*(vol_source-vol_target_un)/vol_source
    
  END IF
  !
  ! Done internal smoothing
  !
  WRITE(*,*) "min/max of terr_target     : ",MINVAL(terr_target),MAXVAL(terr_target)
  
  if (lzero_out_ocean_point_phis) then
    WRITE(*,*) "if ocean mask PHIS=0.0"
  end if
  
  
  sgh_target=0.0
  do count=1,jall
    i    = weights_lgr_index_all(count)!!
    !
    ix  = weights_eul_index_all(count,1)
    iy  = weights_eul_index_all(count,2)
    ip  = weights_eul_index_all(count,3)
    !
    ! convert to 1D indexing of cubed-sphere
    !
    ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix!
    
    wt = weights_all(count,1)
    
    if (lzero_out_ocean_point_phis.AND.landfrac_target(i).lt.0.01_r8) then
      terr_target(i) = 0.0_r8   !5*terr_target(i)
    end if
    sgh_target(i) = sgh_target(i)+wt*((terr_target(i)-terr(ii))**2)/area_target(i)
  end do
  
  
  
  !
  ! zero out small values
  !
  DO i=1,ntarget
    IF (landfrac_target(i)<.001_r8) landfrac_target(i) = 0.0
    IF (sgh_target(i)<0.5) sgh_target(i) = 0.0
    IF (sgh30_target(i)<0.5) sgh30_target(i) = 0.0
  END DO
  sgh_target = SQRT(sgh_target)
  sgh30_target = SQRT(sgh30_target)
  WRITE(*,*) "min/max of sgh_target     : ",MINVAL(sgh_target),MAXVAL(sgh_target)
  WRITE(*,*) "min/max of sgh30_target   : ",MINVAL(sgh30_target),MAXVAL(sgh30_target)
  
  DEALLOCATE(terr,weights_all,weights_eul_index_all,landfrac,landm_coslat)
  
  
  IF (ltarget_latlon) THEN
    CALL wrtncdf_rll(nlon,nlat,lpole,ntarget,terr_target,landfrac_target,sgh_target,sgh30_target,&
         landm_coslat_target,target_center_lon,target_center_lat,.true.)
  ELSE
    CALL wrtncdf_unstructured(ntarget,terr_target,landfrac_target,sgh_target,sgh30_target,&
         landm_coslat_target,target_center_lon,target_center_lat)
  END IF
  DEALLOCATE(terr_target,landfrac_target,sgh30_target,sgh_target,landm_coslat_target)
  
end program convterr
  
!
!
!
subroutine wrtncdf_unstructured(n,terr,landfrac,sgh,sgh30,landm_coslat,lon,lat)
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  
#     include         <netcdf.inc>
  
  !
  ! Dummy arguments
  !
  integer, intent(in) :: n
  real(r8),dimension(n)  , intent(in) :: terr, landfrac,sgh,sgh30,lon, lat, landm_coslat
  !
  ! Local variables
  !
  character (len=64) :: fout       ! NetCDF output file
  integer            :: foutid     ! Output file id
  integer            :: lonid, lonvid
  integer            :: latid, latvid
  integer            :: terrid,nid
  integer            :: terrdim,landfracid,sghid,sgh30id,landm_coslatid
  integer            :: status    ! return value for error control of netcdf routin
  integer            :: i,j
  integer, dimension(2) :: nc_lat_vid,nc_lon_vid
  character (len=8)  :: datestring
  integer :: nc_gridcorn_id, lat_vid, lon_vid
  
  real(r8), parameter :: fillvalue = 1.d36
  
  fout='new-topo-file.nc'
  !
  !  Create NetCDF file for output
  !
  print *,"Create NetCDF file for output"
  status = nf_create (fout, NF_64BIT_OFFSET , foutid)
  if (status .ne. NF_NOERR) call handle_err(status)
  !
  ! Create dimensions for output
  !
  status = nf_def_dim (foutid, 'ncol', n, nid)
  if (status .ne. NF_NOERR) call handle_err(status)
  !
  ! Create variable for output
  !
  print *,"Create variable for output"
  status = nf_def_var (foutid,'PHIS', NF_DOUBLE, 1, nid, terrid)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  status = nf_def_var (foutid,'LANDFRAC', NF_DOUBLE, 1, nid, landfracid)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  status = nf_def_var (foutid,'SGH', NF_DOUBLE, 1, nid, sghid)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  status = nf_def_var (foutid,'SGH30', NF_DOUBLE, 1, nid, sgh30id)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  status = nf_def_var (foutid,'LANDM_COSLAT', NF_DOUBLE, 1, nid, landm_coslatid)
  if (status .ne. NF_NOERR) call handle_err(status)
  !
  status = nf_def_var (foutid,'lat', NF_DOUBLE, 1, nid, latvid)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  status = nf_def_var (foutid,'lon', NF_DOUBLE, 1, nid, lonvid)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  !
  ! Create attributes for output variables
  !
  status = nf_put_att_text (foutid,terrid,'long_name', 21, 'surface geopotential')
  status = nf_put_att_text (foutid,terrid,'units', 5, 'm2/s2')
  status = nf_put_att_double (foutid, terrid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, terrid, '_FillValue'   , nf_double, 1, fillvalue)
  !        status = nf_put_att_text (foutid,terrid,'filter', 35, 'area averaged from USGS 30-sec data')
  
  status = nf_put_att_double (foutid, sghid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, sghid, '_FillValue'   , nf_double, 1, fillvalue)
  status = nf_put_att_text   (foutid, sghid, 'long_name' , 48, &
       'standard deviation of 3km cubed-sphere elevation and target grid elevation')
  status = nf_put_att_text   (foutid, sghid, 'units'     , 1, 'm')
  !        status = nf_put_att_text   (foutid, sghid, 'filter'    , 4, 'none')
  
  status = nf_put_att_double (foutid, sgh30id, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, sgh30id, '_FillValue'   , nf_double, 1, fillvalue)
  status = nf_put_att_text   (foutid, sgh30id, 'long_name' , 49, &
       'standard deviation of 30s elevation from 3km cubed-sphere cell average height')
  status = nf_put_att_text   (foutid, sgh30id, 'units'     , 1, 'm')
  !        status = nf_put_att_text   (foutid, sgh30id, 'filter'    , 4, 'none')
  
  status = nf_put_att_double (foutid, landm_coslatid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, landm_coslatid, '_FillValue'   , nf_double, 1, fillvalue)
  status = nf_put_att_text   (foutid, landm_coslatid, 'long_name' , 23, 'smoothed land fraction')
  status = nf_put_att_text   (foutid, landm_coslatid, 'filter'    , 4, 'none')
  
  status = nf_put_att_double (foutid, landfracid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, landfracid, '_FillValue'   , nf_double, 1, fillvalue)
  status = nf_put_att_text   (foutid, landfracid, 'long_name', 21, 'gridbox land fraction')
  !        status = nf_put_att_text   (foutid, landfracid, 'filter', 40, 'area averaged from 30-sec USGS raw data')
  
  
  status = nf_put_att_text (foutid,latvid,'long_name', 8, 'latitude')
  if (status .ne. NF_NOERR) call handle_err(status)
  status = nf_put_att_text (foutid,latvid,'units', 13, 'degrees_north')
  if (status .ne. NF_NOERR) call handle_err(status)
  !        status = nf_put_att_text (foutid,latvid,'units', 21, 'cell center locations')
  !        if (status .ne. NF_NOERR) call handle_err(status)
  
  status = nf_put_att_text (foutid,lonvid,'long_name', 9, 'longitude')
  if (status .ne. NF_NOERR) call handle_err(status)
  status = nf_put_att_text (foutid,lonvid,'units', 12, 'degrees_east')
  if (status .ne. NF_NOERR) call handle_err(status)
  !        status = nf_put_att_text (foutid,lonvid,'units' , 21, 'cell center locations')
  !        if (status .ne. NF_NOERR) call handle_err(status)
  
  status = nf_put_att_text (foutid,NF_GLOBAL,'source', 50, 'USGS 30-sec dataset binned to ncube3000 (cube-sphere) grid')
  if (status .ne. NF_NOERR) call handle_err(status)
  status = nf_put_att_text (foutid,NF_GLOBAL,'title',  24, '30-second USGS topo data')
  if (status .ne. NF_NOERR) call handle_err(status)
  call DATE_AND_TIME(DATE=datestring)
  status = nf_put_att_text (foutid,NF_GLOBAL,'history',25, 'Written on date: ' // datestring )
  if (status .ne. NF_NOERR) call handle_err(status)
  
  !
  ! End define mode for output file
  !
  status = nf_enddef (foutid)
  if (status .ne. NF_NOERR) call handle_err(status)
  !
  ! Write variable for output
  !
  print*,"writing terrain data",MINVAL(terr),MAXVAL(terr)
  status = nf_put_var_double (foutid, terrid, terr*9.80616)
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing terrain data"
  
  print*,"writing landfrac data",MINVAL(landfrac),MAXVAL(landfrac)
  status = nf_put_var_double (foutid, landfracid, landfrac)
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing landfrac data"
  
  print*,"writing sgh data",MINVAL(sgh),MAXVAL(sgh)
  status = nf_put_var_double (foutid, sghid, sgh)
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing sgh data"
  
  print*,"writing sgh30 data",MINVAL(sgh30),MAXVAL(sgh30)
  status = nf_put_var_double (foutid, sgh30id, sgh30)
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing sgh30 data"
  
  print*,"writing landm_coslat data",MINVAL(landm_coslat),MAXVAL(landm_coslat)
  status = nf_put_var_double (foutid, landm_coslatid, landm_coslat)
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing sgh30 data"
  !
  print*,"writing lat data"
  status = nf_put_var_double (foutid, latvid, lat)
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing lat data"
  
  print*,"writing lon data"
  status = nf_put_var_double (foutid, lonvid, lon)
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing lon data"
  !
  ! Close output file
  !
  print *,"close file"
  status = nf_close (foutid)
  if (status .ne. NF_NOERR) call handle_err(status)
end subroutine wrtncdf_unstructured
!
!**************************************************************     
! 
! if target grid is lat-lon output structured
!
!**************************************************************     
!
subroutine wrtncdf_rll(nlon,nlat,lpole,n,terr_in,landfrac_in,sgh_in,sgh30_in,landm_coslat_in,lon,lat,lprepare_fv_smoothing_routine)
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  
#     include         <netcdf.inc>
  
  !
  ! Dummy arguments
  !
  integer, intent(in) :: n,nlon,nlat
  !
  ! lprepare_fv_smoothing_routine is to make a NetCDF file that can be used with the CAM-FV smoothing software
  !
  logical , intent(in) :: lpole,lprepare_fv_smoothing_routine
  real(r8),dimension(n)  , intent(in) :: terr_in, landfrac_in,sgh_in,sgh30_in,lon, lat, landm_coslat_in
  !
  ! Local variables
  !
  character (len=32) :: fout       ! NetCDF output file
  integer            :: foutid     ! Output file id
  integer            :: lonid, lonvid
  integer            :: latid, latvid
  integer            :: terrid,nid
  integer            :: terrdim,landfracid,sghid,sgh30id,landm_coslatid
  integer            :: status    ! return value for error control of netcdf routin
  integer            :: i,j
  integer, dimension(2) :: nc_lat_vid,nc_lon_vid
  character (len=8)  :: datestring
  integer :: nc_gridcorn_id, lat_vid, lon_vid
  real(r8), parameter :: fillvalue = 1.d36
  real(r8) :: ave
  
  real(r8),dimension(nlon) :: lonar       ! longitude array
  real(r8),dimension(nlat) :: latar       ! latitude array
  
  integer, dimension(2) :: htopodim,landfdim,sghdim,sgh30dim,landmcoslatdim
  real(r8),dimension(n) :: terr, landfrac,sgh,sgh30,landm_coslat
  
  IF (nlon*nlat.NE.n) THEN
    WRITE(*,*) "inconsistent input for wrtncdf_rll"
    STOP
  END IF
  !
  ! we assume that the unstructured layout of the lat-lon grid is ordered in latitude rows, that is,
  ! unstructured index n is given by
  !
  !   n = (j-1)*nlon+i
  !
  ! where j is latitude index and i longitude index
  !
  do i = 1,nlon
    lonar(i)=  lon(i)
  enddo
  do j = 1,nlat
    latar(j)= lat((j-1)*nlon+1)
  enddo
  
  terr = terr_in
  sgh=sgh_in
  sgh30 =sgh30_in
  landfrac = landfrac_in
  landm_coslat = landm_coslat_in
  
  if (lpole) then
    write(*,*) "average pole control volume"
    !
    ! North pole - terr
    !
    ave = 0.0
    do i=1,nlon
      ave = ave + terr_in(i)
    end do
    terr(1:nlon) = ave/DBLE(nlon)
    !
    ! South pole
    !
    ave = 0.0
    do i=n-(nlon+1),n
      ave = ave + terr_in(i)
    end do
    terr(n-(nlon+1):n) = ave/DBLE(nlon)
    
    !
    ! North pole - sgh
    !
    ave = 0.0
    do i=1,nlon
      ave = ave + sgh_in(i)
    end do
    sgh(1:nlon) = ave/DBLE(nlon)
    !
    ! South pole
    !
    ave = 0.0
    do i=n-(nlon+1),n
      ave = ave + sgh_in(i)
    end do
    sgh(n-(nlon+1):n) = ave/DBLE(nlon)
    
    !
    ! North pole - sgh30
    !
    ave = 0.0
    do i=1,nlon
      ave = ave + sgh30_in(i)
    end do
    sgh30(1:nlon) = ave/DBLE(nlon)
    !
    ! South pole
    !
    ave = 0.0
    do i=n-(nlon+1),n
      ave = ave + sgh30_in(i)
    end do
    sgh30(n-(nlon+1):n) = ave/DBLE(nlon)
    
    !
    ! North pole - landfrac
    !
    ave = 0.0
    do i=1,nlon
      ave = ave + landfrac_in(i)
    end do
    landfrac(1:nlon) = ave/DBLE(nlon)
    !
    ! South pole
    !
    ave = 0.0
    do i=n-(nlon+1),n
      ave = ave + landfrac_in(i)
    end do
    landfrac(n-(nlon+1):n) = ave/DBLE(nlon)
    
    !
    ! North pole - landm_coslat
    !
    ave = 0.0
    do i=1,nlon
      ave = ave + landm_coslat_in(i)
    end do
    landm_coslat(1:nlon) = ave/DBLE(nlon)
    !
    ! South pole
    !
    ave = 0.0
    do i=n-(nlon+1),n
      ave = ave + landm_coslat_in(i)
    end do
    landm_coslat(n-(nlon+1):n) = ave/DBLE(nlon)
  end if
  
  
  fout='final.nc'
  !
  !  Create NetCDF file for output
  !
  print *,"Create NetCDF file for output"
  status = nf_create (fout, NF_64BIT_OFFSET , foutid)
  if (status .ne. NF_NOERR) call handle_err(status)
  !
  ! Create dimensions for output
  !
  print *,"Create dimensions for output"
  status = nf_def_dim (foutid, 'lon', nlon, lonid)
  if (status .ne. NF_NOERR) call handle_err(status)
  status = nf_def_dim (foutid, 'lat', nlat, latid)
  if (status .ne. NF_NOERR) call handle_err(status)
  !
  ! Create variable for output
  !
  print *,"Create variable for output"
  
  htopodim(1)=lonid
  htopodim(2)=latid

  if (lprepare_fv_smoothing_routine) then
    status = nf_def_var (foutid,'htopo', NF_DOUBLE, 2, htopodim, terrid)
  else
    status = nf_def_var (foutid,'PHIS', NF_DOUBLE, 2, htopodim, terrid)
  end if
  if (status .ne. NF_NOERR) call handle_err(status)
  
  landfdim(1)=lonid
  landfdim(2)=latid
  
  if (lprepare_fv_smoothing_routine) then
    status = nf_def_var (foutid,'ftopo', NF_DOUBLE, 2, landfdim, landfracid)
  else
    status = nf_def_var (foutid,'LANDFRAC', NF_DOUBLE, 2, landfdim, landfracid)
  end if

  if (status .ne. NF_NOERR) call handle_err(status)
  
  sghdim(1)=lonid
  sghdim(2)=latid
  
  status = nf_def_var (foutid,'SGH', NF_DOUBLE, 2, sghdim, sghid)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  sgh30dim(1)=lonid
  sgh30dim(2)=latid
  
  status = nf_def_var (foutid,'SGH30', NF_DOUBLE, 2, sgh30dim, sgh30id)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  landmcoslatdim(1)=lonid
  landmcoslatdim(2)=latid
  
  status = nf_def_var (foutid,'LANDM_COSLAT', NF_DOUBLE, 2, landmcoslatdim, landm_coslatid)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  status = nf_def_var (foutid,'lat', NF_DOUBLE, 1, latid, latvid)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  status = nf_def_var (foutid,'lon', NF_DOUBLE, 1, lonid, lonvid)
  if (status .ne. NF_NOERR) call handle_err(status)
  
  !
  ! Create attributes for output variables
  !
  status = nf_put_att_text (foutid,terrid,'long_name', 21, 'surface geopotential')
  status = nf_put_att_text (foutid,terrid,'units', 5, 'm2/s2')
  status = nf_put_att_text (foutid,terrid,'filter', 35, 'area averaged from ncube3000 data')
  status = nf_put_att_double (foutid, terrid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, terrid, '_FillValue'   , nf_double, 1, fillvalue)
  
  
  status = nf_put_att_double (foutid, sghid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, sghid, '_FillValue'   , nf_double, 1, fillvalue)
  status = nf_put_att_text   (foutid, sghid, 'long_name' , 48, &
       'standard deviation of 3km cubed-sphere elevation and target grid elevation')
  status = nf_put_att_text   (foutid, sghid, 'units'     , 1, 'm')
  status = nf_put_att_text   (foutid, sghid, 'filter'    , 4, 'none')
  
  status = nf_put_att_double (foutid, sgh30id, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, sgh30id, '_FillValue'   , nf_double, 1, fillvalue)
  status = nf_put_att_text   (foutid, sgh30id, 'long_name' , 49, &
       'standard deviation of 30s elevation from 3km cubed-sphere cell average height')
  status = nf_put_att_text   (foutid, sgh30id, 'units'     , 1, 'm')
  status = nf_put_att_text   (foutid, sgh30id, 'filter'    , 4, 'none')
  
  status = nf_put_att_double (foutid, landm_coslatid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, landm_coslatid, '_FillValue'   , nf_double, 1, fillvalue)
  status = nf_put_att_text   (foutid, landm_coslatid, 'long_name' , 23, 'smoothed land fraction')
  status = nf_put_att_text   (foutid, landm_coslatid, 'filter'    , 4, 'none')
  
  status = nf_put_att_double (foutid, landfracid, 'missing_value', nf_double, 1, fillvalue)
  status = nf_put_att_double (foutid, landfracid, '_FillValue'   , nf_double, 1, fillvalue)
  status = nf_put_att_text   (foutid, landfracid, 'long_name', 21, 'gridbox land fraction')
  status = nf_put_att_text   (foutid, landfracid, 'filter', 40, 'area averaged from 30-sec USGS raw data')
  
  
  status = nf_put_att_text (foutid,latvid,'long_name', 8, 'latitude')
  if (status .ne. NF_NOERR) call handle_err(status)
  status = nf_put_att_text (foutid,latvid,'units', 13, 'degrees_north')
  if (status .ne. NF_NOERR) call handle_err(status)
  !        status = nf_put_att_text (foutid,latvid,'units', 21, 'cell center locations')
  !        if (status .ne. NF_NOERR) call handle_err(status)
  
  status = nf_put_att_text (foutid,lonvid,'long_name', 9, 'longitude')
  if (status .ne. NF_NOERR) call handle_err(status)
  status = nf_put_att_text (foutid,lonvid,'units', 12, 'degrees_east')
  if (status .ne. NF_NOERR) call handle_err(status)
  !        status = nf_put_att_text (foutid,lonvid,'units' , 21, 'cell center locations')
  !        if (status .ne. NF_NOERR) call handle_err(status)
  
  status = nf_put_att_text (foutid,NF_GLOBAL,'source', 27, 'USGS 30-sec dataset GTOPO30')
  if (status .ne. NF_NOERR) call handle_err(status)
  status = nf_put_att_text (foutid,NF_GLOBAL,'title',  24, '30-second USGS topo data')
  if (status .ne. NF_NOERR) call handle_err(status)
  call DATE_AND_TIME(DATE=datestring)
  status = nf_put_att_text (foutid,NF_GLOBAL,'history',25, 'Written on date: ' // datestring )
  if (status .ne. NF_NOERR) call handle_err(status)
  
  !
  ! End define mode for output file
  !
  status = nf_enddef (foutid)
  if (status .ne. NF_NOERR) call handle_err(status)
  !
  ! Write variable for output
  !
  print*,"writing terrain data",MINVAL(terr),MAXVAL(terr)
  if (lprepare_fv_smoothing_routine) then
    status = nf_put_var_double (foutid, terrid, terr)
  else
    status = nf_put_var_double (foutid, terrid, terr*9.80616)
  end if
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing terrain data"
  
  print*,"writing landfrac data",MINVAL(landfrac),MAXVAL(landfrac)
  status = nf_put_var_double (foutid, landfracid, landfrac)
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing landfrac data"
  
  print*,"writing sgh data",MINVAL(sgh),MAXVAL(sgh)
  status = nf_put_var_double (foutid, sghid, sgh)
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing sgh data"
  
  print*,"writing sgh30 data",MINVAL(sgh30),MAXVAL(sgh30)
  status = nf_put_var_double (foutid, sgh30id, sgh30)
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing sgh30 data"
  
  print*,"writing landm_coslat data",MINVAL(landm_coslat),MAXVAL(landm_coslat)
  status = nf_put_var_double (foutid, landm_coslatid, landm_coslat)
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing sgh30 data"
  !
  print*,"writing lat data"
  status = nf_put_var_double (foutid, latvid, latar)
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing lat data"
  
  print*,"writing lon data"
  status = nf_put_var_double (foutid, lonvid, lonar)
  if (status .ne. NF_NOERR) call handle_err(status)
  print*,"done writing lon data"
  !
  ! Close output file
  !
  print *,"close file"
  status = nf_close (foutid)
  if (status .ne. NF_NOERR) call handle_err(status)
end subroutine wrtncdf_rll
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


SUBROUTINE coarsen(f,fcoarse,nf,n,dA_coarse)
  use shr_kind_mod, only: r8 => shr_kind_r8
  IMPLICIT NONE
  REAL (R8), DIMENSION(n)       , INTENT(IN)  :: f
  REAL (R8), DIMENSION(n/nf), INTENT(OUT) :: fcoarse             
  INTEGER, INTENT(in) :: n,nf
  REAL(R8), DIMENSION(INT(SQRT(DBLE(n/6)))/nf,INT(SQRT(DBLE(n/6)))/nf),INTENT(OUT) :: dA_coarse
  !must be an even number
  !
  ! local workspace
  !
  ! ncube = INT(SQRT(DBLE(n/6)))
  
  REAL(R8), DIMENSION(INT(SQRT(DBLE(n/6))),INT(SQRT(DBLE(n/6)))):: dA        
  REAL (R8)    :: sum, sum_area,tmp
  INTEGER  :: jx,jy,jp,ii,ii_coarse,coarse_ncube,ncube
  INTEGER  :: jx_coarse,jy_coarse,jx_s,jy_s
  
  
  !        REAL(R8), DIMENSION(INT(SQRT(DBLE(n/6)))/nf,INT(SQRT(DBLE(n/6)))/nf) :: dAtmp
  
  ncube = INT(SQRT(DBLE(n/6)))
  coarse_ncube = ncube/nf
  
  IF (ABS(DBLE(ncube)/DBLE(nf)-coarse_ncube)>0.000001) THEN
    WRITE(*,*) "ncube/nf must be an integer"
    WRITE(*,*) "ncube and nf: ",ncube,nf
    STOP
  END IF
  
  da_coarse = 0.0
  
  WRITE(*,*) "compute all areas"
  CALL EquiangularAllAreas(ncube, dA)
  !        CALL EquiangularAllAreas(coarse_ncube, dAtmp)!dbg
  tmp = 0.0
  DO jp=1,6
    DO jy_coarse=1,coarse_ncube
      DO jx_coarse=1,coarse_ncube
        !
        ! inner loop
        !
        sum      = 0.0
        sum_area = 0.0
        DO jy_s=1,nf
          jy = (jy_coarse-1)*nf+jy_s
          DO jx_s=1,nf
            jx = (jx_coarse-1)*nf+jx_s
            ii = (jp-1)*ncube*ncube+(jy-1)*ncube+jx
            sum      = sum     +f(ii)*dA(jx,jy)
            sum_area = sum_area+dA(jx,jy)
            !                  WRITE(*,*) "jx,jy",jx,jy
          END DO
        END DO
        tmp = tmp+sum_area
        da_coarse(jx_coarse,jy_coarse) = sum_area
        !              WRITE(*,*) "jx_coarse,jy_coarse",jx_coarse,jy_coarse,&
        !                   da_coarse(jx_coarse,jy_coarse)-datmp(jx_coarse,jy_coarse)
        ii_coarse = (jp-1)*coarse_ncube*coarse_ncube+(jy_coarse-1)*coarse_ncube+jx_coarse
        fcoarse(ii_coarse) = sum/sum_area 
      END DO
    END DO
  END DO
  WRITE(*,*) "coarsened surface area",tmp-4.0*3.141592654
END SUBROUTINE COARSEN

SUBROUTINE overlap_weights(weights_lgr_index_all,weights_eul_index_all,weights_all,&
     jall,ncube,ngauss,ntarget,ncorner,jmax_segments,target_corner_lon,target_corner_lat,nreconstruction)
  use shr_kind_mod, only: r8 => shr_kind_r8
  use remap
  IMPLICIT NONE
  
  
  INTEGER, INTENT(INOUT) :: jall !anticipated number of weights
  INTEGER, INTENT(IN)    :: ncube, ngauss, ntarget, jmax_segments, ncorner, nreconstruction
  
  INTEGER, DIMENSION(jall,3), INTENT(OUT) :: weights_eul_index_all
  REAL(R8), DIMENSION(jall,nreconstruction)  , INTENT(OUT) :: weights_all
  INTEGER, DIMENSION(jall)  , INTENT(OUT) :: weights_lgr_index_all
  
  REAL(R8), DIMENSION(ncorner,ntarget), INTENT(IN) :: target_corner_lon, target_corner_lat
  
  INTEGER,  DIMENSION(ncorner+1) :: ipanel_array, ipanel_tmp
  REAL(R8), DIMENSION(ncorner)  :: lat, lon
  REAL(R8), DIMENSION(0:ncube+2):: xgno, ygno
  REAL(R8), DIMENSION(0:ncorner+1) :: xcell, ycell
  
  REAL(R8), DIMENSION(ngauss) :: gauss_weights, abscissae
  
  REAL(R8) :: da, tmp, alpha, beta
  REAL    (r8), PARAMETER :: pi    = 3.14159265358979323846264338327
  REAL    (r8), PARAMETER :: piq   = 0.25*pi
  REAL    (r8), PARAMETER :: pih   = 0.50*pi
  INTEGER :: i, j,ncorner_this_cell,k,ip,ipanel,ii,jx,jy,jcollect
  integer :: alloc_error
  
  REAL    (r8), PARAMETER :: rad2deg   = 180.0/pi
  
  real(r8), allocatable, dimension(:,:) :: weights
  integer , allocatable, dimension(:,:) :: weights_eul_index
  
  
  LOGICAL:: ldbg = .FAlSE.
  
  INTEGER :: jall_anticipated
  
  jall_anticipated = jall
  
  ipanel_array = -99
  !
  da = pih/DBLE(ncube)
  xgno(0) = -bignum
  DO i=1,ncube+1
    xgno(i) = TAN(-piq+(i-1)*da)
  END DO
  xgno(ncube+2) = bignum
  ygno = xgno
  
  CALL glwp(ngauss,gauss_weights,abscissae)
  
  
  allocate (weights(jmax_segments,nreconstruction),stat=alloc_error )
  allocate (weights_eul_index(jmax_segments,2),stat=alloc_error )
  
  tmp = 0.0
  jall = 1
  DO i=1,ntarget
    WRITE(*,*) "cell",i,"  ",100.0*DBLE(i)/DBLE(ntarget),"% done"
    !
    !---------------------------------------------------          
    !
    ! determine how many vertices the cell has
    !
    !---------------------------------------------------
    !
    CALL remove_duplicates_latlon(ncorner,target_corner_lon(:,i),target_corner_lat(:,i),&
         ncorner_this_cell,lon,lat,1.0E-10,ldbg)
    
    IF (ldbg) THEN
      WRITE(*,*) "number of vertices ",ncorner_this_cell
      WRITE(*,*) "vertices locations lon,",lon(1:ncorner_this_cell)*rad2deg
      WRITE(*,*) "vertices locations lat,",lat(1:ncorner_this_cell)*rad2deg
      DO j=1,ncorner_this_cell
        WRITE(*,*) lon(j)*rad2deg, lat(j)*rad2deg
      END DO
      WRITE(*,*) "  "
    END IF
    !
    !---------------------------------------------------
    !
    ! determine how many and which panels the cell spans
    !
    !---------------------------------------------------          
    !
    DO j=1,ncorner_this_cell
      CALL CubedSphereABPFromRLL(lon(j), lat(j), alpha, beta, ipanel_tmp(j), .TRUE.)
      IF (ldbg) WRITE(*,*) "ipanel for corner ",j," is ",ipanel_tmp(j)
    END DO
    ipanel_tmp(ncorner_this_cell+1) = ipanel_tmp(1)
    ! make sure to include possible overlap areas not on the face the vertices are located
    IF (MINVAL(lat(1:ncorner_this_cell))<-pi/6.0) THEN
      ! include South-pole panel in search
      ipanel_tmp(ncorner_this_cell+1) = 5
      IF (ldbg) WRITE(*,*)  "add panel 5 to search"
    END IF
    IF (MAXVAL(lat(1:ncorner_this_cell))>pi/6.0) THEN
      ! include North-pole panel in search
      ipanel_tmp(ncorner_this_cell+1) = 6
      IF (ldbg) WRITE(*,*)  "add panel 6 to search"
    END IF
    !
    ! remove duplicates in ipanel_tmp
    !
    CALL remove_duplicates_integer(ncorner_this_cell+1,ipanel_tmp(1:ncorner_this_cell+1),&
         k,ipanel_array(1:ncorner_this_cell+1))
    !
    !---------------------------------------------------
    !
    ! loop over panels with possible overlap areas
    !
    !---------------------------------------------------          
    !
    DO ip = 1,k
      ipanel = ipanel_array(ip)
      DO j=1,ncorner_this_cell
        ii = ipanel
        CALL CubedSphereABPFromRLL(lon(j), lat(j), alpha, beta, ii,.FALSE.)            
        IF (j==1) THEN
          jx = CEILING((alpha + piq) / da)
          jy = CEILING((beta  + piq) / da)
        END IF
        xcell(ncorner_this_cell+1-j) = TAN(alpha)
        ycell(ncorner_this_cell+1-j) = TAN(beta)
      END DO
      xcell(0) = xcell(ncorner_this_cell)
      ycell(0) = ycell(ncorner_this_cell)
      xcell(ncorner_this_cell+1) = xcell(1)
      ycell(ncorner_this_cell+1) = ycell(1)
      
      jx = MAX(MIN(jx,ncube+1),0)
      jy = MAX(MIN(jy,ncube+1),0)
      
      CALL compute_weights_cell(xcell(0:ncorner_this_cell+1),ycell(0:ncorner_this_cell+1),&
           jx,jy,nreconstruction,xgno,ygno,&
           1, ncube+1, 1,ncube+1, tmp,&
           ngauss,gauss_weights,abscissae,weights,weights_eul_index,jcollect,jmax_segments,&
           ncube,0,ncorner_this_cell,ldbg)
      
      weights_all(jall:jall+jcollect-1,1:nreconstruction)  = weights(1:jcollect,1:nreconstruction)
      
      weights_eul_index_all(jall:jall+jcollect-1,1:2) = weights_eul_index(1:jcollect,:)
      weights_eul_index_all(jall:jall+jcollect-1,  3) = ipanel
      weights_lgr_index_all(jall:jall+jcollect-1    ) = i
      
      jall = jall+jcollect
      IF (jall>jall_anticipated) THEN
        WRITE(*,*) "more weights than anticipated"
        WRITE(*,*) "increase jall"
        STOP
      END IF
      IF (ldbg) WRITE(*,*) "jcollect",jcollect
    END DO
  END DO
  jall = jall-1
  WRITE(*,*) "sum of all weights divided by surface area of sphere  =",tmp/(4.0*pi)
  WRITE(*,*) "actual number of weights",jall
  WRITE(*,*) "anticipated number of weights",jall_anticipated
  IF (jall>jall_anticipated) THEN
    WRITE(*,*) "anticipated number of weights < actual number of weights"
    WRITE(*,*) "increase jall!"
    STOP
  END IF
  WRITE(*,*) MINVAL(weights_all(1:jall,1)),MAXVAL(weights_all(1:jall,1))
  IF (ABS(tmp/(4.0*pi))-1.0>0.001) THEN
    WRITE(*,*) "sum of all weights does not match the surface area of the sphere"
    WRITE(*,*) "sum of all weights is : ",tmp
    WRITE(*,*) "surface area of sphere: ",4.0*pi
    STOP
  END IF
END SUBROUTINE overlap_weights


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
SUBROUTINE CubedSphereABPFromRLL(lon, lat, alpha, beta, ipanel, ldetermine_panel)
  use shr_kind_mod, only: r8 => shr_kind_r8
  IMPLICIT NONE
  
  REAL    (R8), INTENT(IN)  :: lon, lat
  REAL    (R8), INTENT(OUT) :: alpha, beta
  INTEGER :: ipanel
  LOGICAL, INTENT(IN) :: ldetermine_panel
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
  IF (ldetermine_panel) THEN
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
  ELSE
    IF (ipanel  ==  6) THEN
      sx = yy; sy = -xx; sz = zz
    ELSEIF (ipanel  == 5) THEN
      sx = yy; sy = xx; sz = -zz
    ELSEIF (ipanel == 1) THEN
      sx = yy; sy = zz; sz = xx        
    ELSEIF (ipanel == 3) THEN
      sx = -yy; sy = zz; sz = -xx
    ELSEIF (ipanel == 2) THEN
      sx = -xx; sy = zz; sz = yy
    ELSEIF (ipanel == 4) THEN
      sx = xx; sy = zz; sz = -yy
    ELSE
      WRITE(*,*) "ipanel out of range",ipanel
      STOP
    END IF
  END IF
  
  ! Use panel information to calculate (alpha, beta) coords
  alpha = ATAN(sx / sz)
  beta = ATAN(sy / sz)
  
END SUBROUTINE CubedSphereABPFromRLL

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


SUBROUTINE remove_duplicates_integer(n_in,f_in,n_out,f_out)
  use shr_kind_mod, only: r8 => shr_kind_r8
  integer, intent(in) :: n_in
  integer,dimension(n_in), intent(in) :: f_in
  integer, intent(out) :: n_out
  integer,dimension(n_in), intent(out) :: f_out
  !
  ! local work space
  !
  integer :: k,i,j
  !
  ! remove duplicates in ipanel_tmp
  !
  k = 1
  f_out(1) = f_in(1)
  outer: do i=2,n_in
    do j=1,k
      !            if (f_out(j) == f_in(i)) then
      if (ABS(f_out(j)-f_in(i))<1.0E-10) then
        ! Found a match so start looking again
        cycle outer
      end if
    end do
    ! No match found so add it to the output
    k = k + 1
    f_out(k) = f_in(i)
  end do outer
  n_out = k
END SUBROUTINE remove_duplicates_integer

SUBROUTINE remove_duplicates_latlon(n_in,lon_in,lat_in,n_out,lon_out,lat_out,tiny,ldbg)
  use shr_kind_mod, only: r8 => shr_kind_r8
  integer, intent(in) :: n_in
  real(r8),dimension(n_in), intent(inout) :: lon_in,lat_in
  real, intent(in) :: tiny
  integer, intent(out) :: n_out
  real(r8),dimension(n_in), intent(out) :: lon_out,lat_out
  logical :: ldbg
  !
  ! local work space
  !
  integer :: k,i,j
  REAL    (r8), PARAMETER :: pi        = 3.14159265358979323846264338327
  REAL    (r8), PARAMETER :: pih       = 0.50*pi
  !
  ! for pole points: make sure the longitudes are identical so that algorithm below works properly
  !
  do i=2,n_in
    if (abs(lat_in(i)-pih)<tiny.or.abs(lat_in(i)+pih)<tiny) then 
      lon_in(i) = lon_in(i-1)    
      write(*,*) "pole fix"
    end if
  end do

  lon_out = -9999999.9
  lat_out = -9999999.9
  !
  k = 1
  lon_out(1) = lon_in(1)
  lat_out(1) = lat_in(1)
  outer: do i=2,n_in
    do j=1,k
      if (ABS(lon_out(j)-lon_in(i))<tiny.AND.ABS(lat_out(j)-lat_in(i))<tiny) then
        ! Found a match so start looking again
        cycle outer
      end if
    end do
    ! No match found so add it to the output
    k = k + 1
    lon_out(k) = lon_in(i)
    lat_out(k) = lat_in(i)
  end do outer
  n_out = k
END SUBROUTINE remove_duplicates_latlon
