! Code that computes control volumes with the same area as the GLL weights
! (for SCRIP) uses analytic area formula.

module comp_gll_ctr_vol
  use shr_kind_mod,           only: r8=>shr_kind_r8, shr_kind_cl
  use cam_abortutils,         only: endrun
  use cam_logfile,            only: iulog
  use shr_sys_mod,            only: shr_sys_flush
  use global_norms_mod,       only: wrap_repro_sum
  use physconst,              only: pi
  use infnan,                 only: isnan

  use coordinate_systems_mod, only: cartesian3d_t, cartesian2d_t
  use coordinate_systems_mod, only: spherical_polar_t, change_coordinates
  use coordinate_systems_mod, only: cubedsphere2cart, cube_face_number_from_cart
  use coordinate_systems_mod, only: distance, sphere_tri_area
  use parallel_mod,           only: global_shared_buf, global_shared_sum
  use edgetype_mod,           only: EdgeBuffer_t, Ghostbuffer3D_t
  use dimensions_mod,         only: np, ne
  use control_mod,            only: fine_ne
  use reduction_mod,          only: red_sum, parallelmin, parallelmax

  implicit none
  private
  save

  character(len=16),  public :: se_write_gll_grid = "no"

  ! nv_max will be set to 2*max_elements_attached_to_node
  !        This works out to 6 for the regular case and 14 for refined meshes
  integer :: nv_max=-99

  type, public :: ctrlvol_t
    real(r8)                             :: vol(np,np)         ! area of the unit sphere covered (local)
    real(r8)                             :: totvol(np,np)      ! area of the unit sphere covered (local)
    real(r8)                             :: invvol(np,np)      ! inverse area (includes neigbors)
    type(cartesian2d_t)                  :: cartp_dual(0:np,0:np)
    type(cartesian3d_t)                  :: cart3d_dual(0:np,0:np)
    type(cartesian3D_t),     allocatable :: vert(:,:,:)        ! bounding box for the polygon
    type(spherical_polar_t), allocatable :: vert_latlon(:,:,:) ! bounding box for the polygon
    integer,                 allocatable :: face_no(:,:,:)     ! face_no of cv vertex.  0 if on cube edge
    integer                              :: nvert(np,np)       ! number of vertex per polygon
  end type ctrlvol_t

  ! tho options:
  ! (1) for NP<>4 or Refined Meshes  (this is less accurate)
  !     build control volumes out of lines which are
  !     always gnomonic coordinate lines.  results in hexagon control volumes
  !     at cube corners and edges.  control volumes at cube-sphere edges are
  !     non-convex, which breaks SCRIP.
  ! iterative option for NP=4 only:
  ! (2) (USE_PENTAGONS)
  !     iterate to minimize difference between spherical area and GLL weight
  !     introduce pentagons in the center of each element to make areas agree
  !     control volumes are triangles, squares or pentagons

  type(ctrlvol_t),    allocatable, target :: cvlist(:)
  type(EdgeBuffer_t)                      :: edge1
  type(GhostBuffer3D_t)                   :: ghost_buf

  ! User interface
  public :: gll_grid_write ! Write the grid in SCRIP format and exit

  ! Private interfaces
  private:: InitControlVolumesData   ! allocates internal data structure
  private:: InitControlVolumes       ! Inits all surfaces: vol,totvol, invvol

  private:: GetVolumeLocal
  private:: GetVolume

  logical, private :: initialized = .false.

CONTAINS

  subroutine gll_grid_write(elem, grid_format, filename_in)
    use netcdf,                 only: nf90_strerror
    use spmd_utils,             only: masterproc, mpicom
    use pio,                    only: var_desc_t, file_desc_t
    use pio,                    only: pio_int, pio_double, PIO_NOERR
    use pio,                    only: pio_put_att, pio_put_var, pio_enddef
    use cam_pio_utils,          only: cam_pio_createfile, cam_pio_closefile
    use cam_pio_utils,          only: cam_pio_def_dim, cam_pio_def_var
    use cam_grid_support,       only: cam_grid_id, cam_grid_dimensions
    use cam_grid_support,       only: cam_grid_write_dist_array
 !!XXgoldyXX: v debug only
#define USE_PIO3D
#ifdef USE_PIO3D
    use pio,                    only: io_desc_t, pio_write_darray, PIO_OFFSET_KIND
    use cam_pio_utils,          only: cam_pio_newdecomp
    use spmd_utils,             only: iam
#endif
 !!XXgoldyXX: ^ debug only

    use hybrid_mod,             only: hybrid_t, config_thread_region
    use parallel_mod,           only: par
    use dimensions_mod,         only: nelem, nelemd
    use control_mod,            only: refined_mesh, fine_ne
    use element_mod,            only: element_t
    use dof_mod,                only: UniquePoints
    use coordinate_systems_mod, only: cart2spherical
    
    ! Inputs
    type(element_t),   intent(in) :: elem(:)
    character(len=*),  intent(in) :: grid_format
    character(len=*),  intent(in) :: filename_in
    
    real(r8), parameter :: rad2deg = 180._r8/pi
    
    ! Local variables
!!XXgoldyXX: v debug only
#ifdef USE_PIO3D
    integer(PIO_OFFSET_KIND), allocatable :: ldof(:)
    integer :: ii, jj
    type(io_desc_t), pointer :: iodesc
#endif
!!XXgoldyXX: ^ debug only
    integer                       :: i, j, ie, ierror, status, ivtx, index
    integer                       :: grid_corners_id, grid_rank_id, grid_size_id
    type(var_desc_t)              :: grid_dims_id, grid_area_id, grid_center_lat_id
    type(var_desc_t)              :: grid_center_lon_id, grid_corner_lat_id
    type(var_desc_t)              :: grid_corner_lon_id, grid_imask_id

    type(file_desc_t)             :: file
    integer                       :: gll_grid       ! Grid ID
    integer                       :: gridsize       ! Total number of unique grid columns
    integer                       :: arr_dims2d(2)  ! (/ np*np, nelemed)
    integer                       :: file_dims2d(1) ! (/ gridsize /)
    integer                       :: arr_dims3d(3)  ! (/ np*np, nv_max, nelemed)
    integer                       :: file_dims3d(2) ! (/ nv_max, gridsize /)

    real(r8),         allocatable :: gwork(:,:,:)   ! np*np, nv_max, nelemd
    type(hybrid_t)                :: hybrid
    character(len=256)            :: errormsg
    character(len=shr_kind_cl)    :: filename
    type(spherical_polar_t)       :: sphere
    character(len=*), parameter   :: subname = 'gll_grid_write'

    !! Check to see if we are doing grid output
    if (trim(grid_format) == "no") then
      if (masterproc) then
        write(iulog, *) subname, ': Not writing phys_grid file.'
      end if
      return
    else if (trim(grid_format) /= 'SCRIP') then
      write(errormsg, *) subname, ': ERROR, bad value for se_write_grid, ', &
           trim(grid_format)
      call endrun(errormsg)
    end if

    ! Set up the control volumes
    if (refined_mesh) then
      nv_max = 14
    else
      nv_max = 5
    end if
    if (masterproc) then
      write(iulog, *) subname, ': computing GLL dual grid for control volumes:'
    end if
    call InitControlVolumesData(par,elem,nelemd)
    ! single thread
    hybrid = config_thread_region(par,'serial')
    call InitControlVolumes(elem,hybrid,1,nelemd)
    if (masterproc) then
      write(6, *) subname, ': done computing GLL dual grid for control volumes.'
    end if

    ! Create the NetCDF file
    if (len_trim(filename_in) == 0) then
      if (refined_mesh) then
        if (fine_ne <= 0) then
          call endrun('gll_grid_write: refined_mesh selected but fine_ne not set')
        else
          write(filename,'(a,i0,a,a,3a)') "ne0np", np, "_refined_", trim(grid_format), ".nc"
        end if
      else
        write(filename, '(a,i0,a,i0,a,a,3a)') "ne", ne, "np", np,                 &
             "_", trim(grid_format), ".nc"
      end if
    else
      filename = trim(filename_in)
    end if
    if (masterproc) then
      write(iulog, *) 'Writing gll SCRIP grid file: ', trim(filename)
      call shr_sys_flush(iulog)
    end if

    call cam_pio_createfile(file, trim(filename))
    gll_grid = cam_grid_id('GLL')
    call cam_grid_dimensions(gll_grid, file_dims3d)
    gridsize = file_dims3d(1)
    file_dims2d(1) = gridsize
    file_dims3d(1) = nv_max
    file_dims3d(2) = gridsize
    arr_dims2d(1) = np*np
    arr_dims2d(2) = nelemd
    arr_dims3d(1) = np*np
    arr_dims3d(2) = nv_max
    arr_dims3d(3) = nelemd
    call cam_pio_def_dim(file, "grid_corners", nv_max,   grid_corners_id)
    call cam_pio_def_dim(file, "grid_rank",    1,        grid_rank_id)
    call cam_pio_def_dim(file, "grid_size",    gridsize, grid_size_id)
    ! Define the coordinate variables
    call cam_pio_def_var(file, "grid_dims", pio_int, (/ grid_rank_id /),  &
           grid_dims_id)

    ! Define grid area
    call cam_pio_def_var(file, "grid_area", pio_double,                 &
           (/grid_size_id/), grid_area_id)
    status = pio_put_att(file, grid_area_id, "units", "radians^2")
    if (status /= pio_noerr) then
      write(iulog, *) subname,': Error defining units attribute for grid_area'
      call shr_sys_flush(iulog)
      call endrun(subname//": "//trim(nf90_strerror(status)))
    end if
    status = pio_put_att(file, grid_area_id, "long_name", "area weights")
    if (status /= pio_noerr) then
      write(iulog, *) subname,': Error defining long_name attribute for grid_area'
      call shr_sys_flush(iulog)
      call endrun(subname//": "//trim(nf90_strerror(status)))
    end if

    ! Define grid center latitudes
    call cam_pio_def_var(file, "grid_center_lat", pio_double,           &
           (/grid_size_id/), grid_center_lat_id)
    status = pio_put_att(file, grid_center_lat_id, "units", "degrees")
    if (status /= pio_noerr) then
      write(iulog, *) subname,': Error defining units attribute for grid_center_lat'
      call shr_sys_flush(iulog)
      call endrun(subname//": "//trim(nf90_strerror(status)))
    end if

    ! Define grid center longitudes
    call cam_pio_def_var(file, "grid_center_lon", pio_double,           &
           (/grid_size_id/), grid_center_lon_id)
    status = pio_put_att(file, grid_center_lon_id, "units", "degrees")
    if (status /= pio_noerr) then
      write(iulog, *) subname,': Error defining units attribute for grid_center_lon'
      call shr_sys_flush(iulog)
      call endrun(subname//": "//trim(nf90_strerror(status)))
    end if

    ! Define grid corner latitudes
    call cam_pio_def_var(file, "grid_corner_lat", pio_double,           &
           (/grid_corners_id, grid_size_id/), grid_corner_lat_id)
    status = pio_put_att(file, grid_corner_lat_id, "units", "degrees")
    if (status /= pio_noerr) then
      write(iulog, *) subname,': Error defining units attribute for grid_corner_lon'
      call shr_sys_flush(iulog)
      call endrun(subname//": "//trim(nf90_strerror(status)))
    end if

    ! Grid corner longitudes
    call cam_pio_def_var(file, "grid_corner_lon", pio_double,           &
           (/grid_corners_id, grid_size_id/), grid_corner_lon_id)
    status = pio_put_att(file, grid_corner_lon_id, "units", "degrees")
    if (status /= pio_noerr) then
      write(iulog, *) subname,': Error defining units attribute for grid_corner_lon'
      call shr_sys_flush(iulog)
      call endrun(subname//": "//trim(nf90_strerror(status)))
    end if

    ! Grid mask
    call cam_pio_def_var(file, "grid_imask", pio_double,                &
           (/grid_size_id/), grid_imask_id)

    ! End of NetCDF definitions
    status = PIO_enddef(file)
    if (status /= pio_noerr) then
      write(iulog, *) subname, ': Error calling enddef'
      call shr_sys_flush(iulog)
      call endrun(subname//": "//trim(nf90_strerror(status)))
    end if

    ! Work array to gather info before writing
    allocate(gwork(np*np, nv_max, nelemd))

    ! Write grid size
    status = pio_put_var(file, grid_dims_id, (/ gridsize /))
    if (status /= pio_noerr) then
      write(iulog, *) subname, ': Error writing variable grid_dims'
      call shr_sys_flush(iulog)
      call endrun(subname//": "//trim(nf90_strerror(status)))
    end if
    ! Write GLL grid areas
    do ie = 1, nelemd
      index = 1
      do j = 1, np
        do i = 1, np
          gwork(index, 1, ie) = cvlist(ie)%vol(i,j)
          index = index + 1
        end do
      end do
    end do
    call cam_grid_write_dist_array(file, gll_grid, arr_dims2d, file_dims2d,   &
         gwork(:,1,:), grid_area_id)
    ! Write GLL grid cell center latitude
    do ie = 1, nelemd
      index = 1
      do j = 1, np
        do i = 1, np
          gwork(index, 1, ie) = elem(ie)%spherep(i,j)%lat * rad2deg
          index = index + 1
        end do
      end do
    end do
    call cam_grid_write_dist_array(file, gll_grid, arr_dims2d, file_dims2d,   &
         gwork(:,1,:), grid_center_lat_id)
    ! Write GLL grid cell center longitude
    do ie = 1, nelemd
      index = 1
      do j = 1, np
        do i = 1, np
          gwork(index, 1, ie) = elem(ie)%spherep(i,j)%lon * rad2deg
          index = index + 1
        end do
      end do
    end do
    call cam_grid_write_dist_array(file, gll_grid, arr_dims2d, file_dims2d,   &
         gwork(:,1,:), grid_center_lon_id)

    ! GLL grid corners
    ! Collect all information for the grid corner latitude (counter-clockwise)
    do ie = 1, nelemd
      do ivtx = 1, nv_max
        index = 1
        do j = 1, np
          do i = 1, np
            gwork(index, ivtx, ie) = cvlist(ie)%vert_latlon(ivtx,i,j)%lat * rad2deg
            index = index + 1
          end do
        end do
      end do
    end do
!!XXgoldyXX: v debug only
#ifdef USE_PIO3D
allocate(ldof(np*np*nelemd*nv_max))
ldof = 0
do ie = 1, nelemd
  do index = 1, elem(ie)%idxP%NumUniquePts
    i = elem(ie)%idxP%ia(index)
    j = elem(ie)%idxP%ja(index)
    ii = (i - 1) + ((j - 1) * np) + ((ie - 1) * np * np * nv_max) + 1
    jj = (elem(ie)%idxP%UniquePtOffset + index - 2) * nv_max
    do ivtx = 1, nv_max
      ldof(ii) = jj + ivtx
      if ((jj+ivtx < 1) .or. (jj+ivtx > gridsize*nv_max)) then
        write(errormsg, '(4(a,i0))') ' ERROR (',iam,'): ldof(',ii,') = ',jj + ivtx,' > ',gridsize*nv_max
        call endrun(subname//trim(errormsg))
      end if
      ii = ii + np*np
    end do
  end do
end do
allocate(iodesc)
call cam_pio_newdecomp(iodesc, (/ nv_max, gridsize /), ldof, PIO_double)
call pio_write_darray(file, grid_corner_lat_id, iodesc, gwork, status)
#else
!!XXgoldyXX: ^ debug only
    call cam_grid_write_dist_array(file, gll_grid, arr_dims3d, file_dims3d,   &
         gwork, grid_corner_lat_id)
!!XXgoldyXX: v debug only
#endif
!!XXgoldyXX: ^ debug only
    ! Collect all information for the grid corner longitude (counter-clockwise)
    do ie = 1, nelemd
      do ivtx = 1, nv_max
        index = 1
        do j = 1, np
          do i = 1, np
            gwork(index, ivtx, ie) = cvlist(ie)%vert_latlon(ivtx,i,j)%lon * rad2deg
            index = index + 1
          end do
        end do
      end do
    end do
!!XXgoldyXX: v debug only
#ifdef USE_PIO3D
call pio_write_darray(file, grid_corner_lon_id, iodesc, gwork, status)
#else
!!XXgoldyXX: ^ debug only
    call cam_grid_write_dist_array(file, gll_grid, arr_dims3d, file_dims3d,   &
         gwork, grid_corner_lon_id)
!!XXgoldyXX: v debug only
#endif
!!XXgoldyXX: ^ debug only
      ! Grid imask
      gwork(:,1,:) = 1.0_r8
      call cam_grid_write_dist_array(file, gll_grid, arr_dims2d, file_dims2d, &
           gwork(:,1,:), grid_imask_id)

    call mpi_barrier(mpicom, ierror)
    ! Close the file
    call cam_pio_closefile(file)
    if(masterproc) then
      write(iulog, *) 'Finished writing physics grid file: ', trim(filename)
      call shr_sys_flush(iulog)
    end if

  end subroutine gll_grid_write

  ! elemid is the local element id (in nets:nete)
  function GetVolume(elemid) result(vol)

    integer,       intent(in) :: elemid
    real(kind=r8), pointer    :: vol(:,:)

    if(.not. initialized) then
      call endrun('Attempt to use volumes prior to initializing')
    end if
    vol => cvlist(elemid)%totvol

  end function GetVolume

  function GetVolumeLocal(elemid) result(vol)

    integer,          intent(in) :: elemid
    real(r8), pointer            :: vol(:,:)

    if(.not. initialized) then
      call endrun('Attempt to use volumes prior to initializing')
    end if
    vol => cvlist(elemid)%vol

  end function GetVolumeLocal

  subroutine InitControlVolumesData(par, elem, nelemd)
    use edge_mod,     only: initedgebuffer, initGhostBuffer3D
    use parallel_mod, only: parallel_t, HME_BNDRY_P2P
    use element_mod,  only: element_t
    use thread_mod,   only: horz_num_threads

    type(parallel_t), intent(in) :: par
    type(element_t),  intent(in) :: elem(:)
    integer,          intent(in) :: nelemd

    integer                      :: ie

    ! Cannot be done in a threaded region
    allocate(cvlist(nelemd))
    do ie = 1, nelemd
      allocate(cvlist(ie)%vert(nv_max, np,np))
      allocate(cvlist(ie)%vert_latlon(nv_max,np,np))
      allocate(cvlist(ie)%face_no(nv_max,np,np))
    end do

    call initedgebuffer(par,edge1,elem,3,bndry_type=HME_BNDRY_P2P, nthreads=horz_num_threads)
    call initGhostBuffer3D(ghost_buf,3,np+1,1)
  end subroutine InitControlVolumesData

  subroutine VerifyAreas(elem,hybrid,nets,nete)

    use element_mod,        only: element_t
    use hybrid_mod,         only: hybrid_t

    integer,         intent(in)         :: nets,nete
    type(element_t), intent(in), target :: elem(:)
    type(hybrid_t),  intent(in)         :: hybrid

    integer                             :: i, j, ie, k, kptr, kmax
    real(r8)                            :: rspheremp(np,np)
    real(r8)                            :: invvol(np,np)
    real(r8)                            :: error, max_error, max_invvol, maxrsphere

    error = 0
    max_error = 0
    do ie=nets,nete
      rspheremp = elem(ie)%rspheremp
      invvol    = cvlist(ie)%invvol
      do j=1,np
        do i=1,np
          error = 100*ABS(rspheremp(i,j)-invvol(i,j))/invvol(i,j)
          if (max_error.lt.error) then
            max_error  = error
            max_invvol = invvol(i,j)
            maxrsphere = rspheremp(i,j)
          end if
        end do
      end do
    end do
    print '(A,F16.4 )',"Control Volume Stats: Max error percent:", max_error
    print '(A,F16.12)',"                     Value From Element:",1/maxrsphere
    print '(A,F16.12)',"              Value From Control Volume:",1/max_invvol
    max_error = parallelmax(max_error,hybrid)
    if (hybrid%masterthread) then
      write(6, '(a,f16.4)') "Control volume area vs. gll area: max error (percent):", max_error
    end if

  end subroutine VerifyAreas


  subroutine InitControlVolumes(elem, hybrid,nets,nete)
    use element_mod,  only: element_t
    use hybrid_mod,   only: hybrid_t
    use control_mod,  only: refined_mesh

    integer,         intent(in)         :: nets,nete
    type(element_t), intent(in), target :: elem(:)
    type(hybrid_t),  intent(in)         :: hybrid

    if (refined_mesh .or. (np /= 4)) then
      call InitControlVolumes_duel(elem, hybrid,nets,nete)
    else
      call InitControlVolumes_gll(elem, hybrid,nets,nete)
      call VerifVolumes(elem, hybrid,nets,nete)
    end if
  end subroutine InitControlVolumes

  subroutine InitControlVolumes_duel(elem, hybrid,nets,nete)
    use bndry_mod,   only: bndry_exchange
    use edge_mod,    only: edgeVpack, edgeVunpack, freeedgebuffer, freeghostbuffer3D
    use element_mod, only: element_t, element_var_coordinates, element_var_coordinates3d
    use hybrid_mod,  only: hybrid_t

    use quadrature_mod,         only: quadrature_t, gausslobatto
    use coordinate_systems_mod, only: cube_face_number_from_sphere

    integer,         intent(in)         :: nets,nete
    type(element_t), intent(in), target :: elem(:)
    type(hybrid_t),  intent(in)         :: hybrid

    type(quadrature_t)  :: gll_pts
    type(cartesian3d_t) :: quad(4),corners3d(4)
    real(r8)            :: cv_pts(0:np) !was kind=longdouble_kind in HOMME
    real(r8)            :: test(np,np,1)

    integer             :: i, j, ie, k, kmax2, kk

    gll_pts = gausslobatto(np)
    ! gll points
    cv_pts(0)=-1
    do i=1,np
      cv_pts(i) = cv_pts(i-1) + gll_pts%weights(i)
    end do
    cv_pts(np)=1
    do i=1,np-1
      if (gll_pts%points(i) > cv_pts(i) .or. cv_pts(i) > gll_pts%points(i+1)) then
        call endrun("Error: CV and GLL points not interleaved")
      end if
    end do


    ! intialize local element areas
    test = 0
    do ie=nets,nete
      cvlist(ie)%cart3d_dual(0:np,0:np) = element_var_coordinates3D(elem(ie)%corners3D, cv_pts)

      ! compute true area of element and SEM area
      cvlist(ie)%vol=0
      do i=1,np
        do j=1,np
          ! (gnomonic coordinate lines only), more accurate
          quad(1) = cvlist(ie)%cart3d_dual(i-1,j-1)
          quad(2) = cvlist(ie)%cart3d_dual(i,j-1)
          quad(3) = cvlist(ie)%cart3d_dual(i,j)
          quad(4) = cvlist(ie)%cart3d_dual(i-1,j)
          cvlist(ie)%vol(i,j) = surfarea(quad,4)
        end do
      end do
      test(:,:,1) = cvlist(ie)%vol(:,:)
      call edgeVpack(edge1,test,1,0,ie)
    end do

    call bndry_exchange(hybrid, edge1)

    test = 0
    do ie=nets,nete
      test(:,:,1) = cvlist(ie)%vol(:,:)
      call edgeVunpack(edge1, test, 1, 0, ie)
      cvlist(ie)%totvol(:,:) = test(:,:,1)
      cvlist(ie)%invvol(:,:)=1.0_r8/cvlist(ie)%totvol(:,:)
    end do

    call VerifyAreas(elem, hybrid, nets, nete)

    ! construct the global CV grid and global CV areas from the
    ! local dual grid (cvlist()%cart_dual) and local areas (cvlist()%vol)
    call construct_cv_duel(elem, hybrid, nets, nete)
    ! compute output needed for SCRIP:  lat/lon coordinates, and for the
    ! control volume with only 3 corners, repeat the last point to make a
    ! degenerate quad.
    kmax2 = 0
    do ie = nets, nete
      kmax2 = MAX(kmax2, MAXVAL(cvlist(ie)%nvert))
    end do
    do ie = nets, nete
      do j = 1, np
        do i = 1, np
          cvlist(ie)%vert_latlon(:,i,j)%lat = 0.0_r8
          cvlist(ie)%vert_latlon(:,i,j)%lon = 0.0_r8
          k = cvlist(ie)%nvert(i,j)
          !
          ! follow SCRIP protocol - of kk>k then repeat last vertex
          !
          do kk = k+1, nv_max
            cvlist(ie)%vert(kk, i, j) = cvlist(ie)%vert(k,i,j)
          end do
          do kk = 1, nv_max
            cvlist(ie)%vert_latlon(kk, i, j) = change_coordinates(cvlist(ie)%vert(kk, i, j))
            cvlist(ie)%face_no(kk, i, j) = cube_face_number_from_sphere(cvlist(ie)%vert_latlon(kk, i, j))
          end do

        end do
      end do
    end do
    ! Release memory
    if(hybrid%masterthread) then
      call freeedgebuffer(edge1)
      call FreeGhostBuffer3D(ghost_buf)
    end if

    initialized=.true.
  end subroutine  InitControlVolumes_duel

  function average(t, n) result(a)

    integer,             intent(in) :: n
    type(cartesian3d_t), intent(in) :: t(n)
    type(cartesian3d_t)             :: a
    integer                         :: i

    a%x = 0._r8
    a%y = 0._r8
    a%z = 0._r8
    do i = 1, n
      a%x = a%x + t(i)%x
      a%y = a%y + t(i)%y
      a%z = a%z + t(i)%z
    end do
    a%x = a%x / n
    a%y = a%y / n
    a%z = a%z / n
    return
  end function  average

  function make_unique(a, n) result(m)

    integer,  intent(in)    :: n
    real(r8), intent(inout) :: a(n)
    integer                 :: m
    integer                 :: i,j
    real(r8)                :: delta

    do i=1,n-1
      do j=i+1,n
        !        if (ABS(a(j)-a(i)).lt. 1e-6)  a(j) = 9999
        delta = abs(a(j)-a(i))
        if (delta < 1.e-6_r8)  a(j) = 9999.0_r8
        if (abs((2.0_r8*pi) - delta) < 1.0e-6_r8)  a(j) = 9999.0_r8
      end do
    end do
    m = 0
    do i=1,n
      if (a(i) < 9000.0_r8) m = m + 1
    end do
    if (mod(m,2).ne.0) then
      do i=1,n
        print *,'angle with centroid: ',i,a(i),mod(a(i),2*pi)
      end do
      call endrun("Error: Found an odd number or nodes for cv element. Should be even.")
    end if
    return
  end function  make_unique

  function SortNodes(t3, n) result(m)
    use coordinate_systems_mod,  only: cube_face_number_from_cart, cart2cubedsphere, change_coordinates


    integer,             intent(in)    :: n
    type(cartesian3d_t), intent(inout) :: t3(n)

    type(cartesian3d_t)                :: c3, t(n)
    type(cartesian2d_t)                :: c2, t2
    real(r8)                           :: angle(n)
    integer                            :: i,j,k,m,f
    integer                            :: ip(n)

    c3 = average(t3, n)
    f  = cube_face_number_from_cart(c3)
    c2 = cart2cubedsphere(c3, f)

    do i=1,n
      t2 = cart2cubedsphere(t3(i), f)
      t2%x = t2%x - c2%x
      t2%y = t2%y - c2%y
      angle(i) = atan2(t2%y, t2%x)
    end do
    m = make_unique(angle,n)
    do i=1,m
      k = 1
      do j=2,n
        if (angle(j)<angle(k)) k=j
      end do
      angle(k) = 9999 ! greater than pi
      ip(i)=k
    end do
    t(1:m) = t3(ip(1:m))
    t3(1:m) = t(1:m)
    return
  end function SortNodes

  subroutine construct_cv_duel(elem,hybrid,nets,nete)
    !
    ! construct global dual grid from local element dual grid cvlist(ie)%cart3d_dual(:,:)
    ! all control volumes will be squares or triangles (at cube corners)
    !
    ! 10/2009: added option to make hexagon control volumes at cube edges and corners
    !
    use element_mod,    only: element_t
    use hybrid_mod,     only: hybrid_t
    use edge_mod,       only: ghostVpack3D,ghostVunpack3d
    use bndry_mod    ,  only: ghost_exchangeVfull
    use dimensions_mod, only: max_corner_elem
    use control_mod,    only: north, south, east, west, neast, nwest, seast, swest

    type(element_t), intent(in), target :: elem(:)
    type(hybrid_t),  intent(in)         :: hybrid
    integer                             :: nets,nete
    !   local
    integer             :: i,j,k,m,n,o,p,ie,m2
    real(r8)            :: vertpack  (    0:np,       0:np,    3)
    real(r8)            :: vertpack2 (    1:np+1,     1:np+1,  3)
    real(r8)            :: vertunpack(   -1:np+1,    -1:np+1,  3)
    real(r8)            :: vertunpack2(   0:np+2,     0:np+2,  3)
    real(r8)            :: sw(   -1:   0,    -1:   0,  3, max_corner_elem-1)
    real(r8)            :: se(   np:np+1,    -1:   0,  3, max_corner_elem-1)
    real(r8)            :: ne(   np:np+1,    np:np+1,  3, max_corner_elem-1)
    real(r8)            :: nw(   -1:   0,    np:np+1,  3, max_corner_elem-1)
    real(r8)            :: sw2(    0:   1,     0:   1,  3, max_corner_elem-1)
    real(r8)            :: se2( np+1:np+2,     0:   1,  3, max_corner_elem-1)
    real(r8)            :: ne2( np+1:np+2,  np+1:np+2,  3, max_corner_elem-1)
    real(r8)            :: nw2(    0:   1,  np+1:np+2,  3, max_corner_elem-1)

    type(cartesian3d_t) :: vert(2*nv_max)
    type(cartesian3d_t) :: pt_3d
    type(cartesian3d_t) :: cv   (-1:np+1,   -1:np+1)
    type(cartesian3d_t) :: cv_sw(-1:   0,   -1:   0,  max_corner_elem-1)
    type(cartesian3d_t) :: cv_se(np:np+1,   -1:   0,  max_corner_elem-1)
    type(cartesian3d_t) :: cv_ne(np:np+1,   np:np+1,  max_corner_elem-1)
    type(cartesian3d_t) :: cv_nw(-1:   0,   np:np+1,  max_corner_elem-1)
    integer             :: mlt(5:8)


    vertpack   = 0
    vertunpack = 0
    vertpack2  = 0
    vertunpack2= 0
    sw         = 0
    se         = 0
    ne         = 0
    nw         = 0
    do ie=nets,nete
      do j=0,np
        do i=0,np
          pt_3d = cvlist(ie)%cart3d_dual(i,j)
          vertpack(i,j,1) = pt_3d%x
          vertpack(i,j,2) = pt_3d%y
          vertpack(i,j,3) = pt_3d%z
        end do
      end do
      do j=0,np
        do i=0,np
          do k=1,3
            vertpack2(i+1,j+1,k) = vertpack(i,j,k)
          end do
        end do
      end do
      call ghostVpack3D(ghost_buf, vertpack2, 3, 0, elem(ie)%desc)
    end do

    call ghost_exchangeVfull(hybrid%par,hybrid%ithr,ghost_buf)
    do ie=nets,nete
      do j=0,np
        do i=0,np
          pt_3d = cvlist(ie)%cart3d_dual(i,j)
          vertunpack(i,j,1) = pt_3d%x
          vertunpack(i,j,2) = pt_3d%y
          vertunpack(i,j,3) = pt_3d%z
        end do
      end do
      do j=0,np
        do i=0,np
          do k=1,3
            vertunpack2(i+1,j+1,k) = vertunpack(i,j,k)
          end do
        end do
      end do
      sw2=0
      se2=0
      nw2=0
      ne2=0
      call ghostVunpack3d(ghost_buf, vertunpack2, 3, 0, elem(ie)%desc, sw2,se2,nw2,ne2,mlt)
      do j=0,np+2
        do i=0,np+2
          do k=1,3
            vertunpack(i-1,j-1,k) = vertunpack2(i,j,k)
          end do
        end do
      end do
      sw=0
      se=0
      nw=0
      ne=0
      do m=1,mlt(swest)-1
        do k=1,3
          do j=0,1
            do i=0,1
              sw(i-1,j-1,k,m) = sw2(i,j,k,m)
            end do
          end do
        end do
      end do
      do m=1,mlt(seast)-1
        do k=1,3
          do j=0,1
            do i=np+1,np+2
              se(i-1,j-1,k,m) = se2(i,j,k,m)
            end do
          end do
        end do
      end do
      do m=1,mlt(nwest)-1
        do k=1,3
          do j=np+1,np+2
            do i=0,1
              nw(i-1,j-1,k,m) = nw2(i,j,k,m)
            end do
          end do
        end do
      end do
      do m=1,mlt(neast)-1
        do k=1,3
          do j=np+1,np+2
            do i=np+1,np+2
              ne(i-1,j-1,k,m) = ne2(i,j,k,m)
            end do
          end do
        end do
      end do
      ! Count and orient vert array
      ! Positive: 1->2->3->4 is counter clockwise on the sphere
      ! Negative: clockwise orientation

      do j=1,np
        do i=1,np
          cvlist(ie)%vert(:,i,j)%x = 0.0_r8
          cvlist(ie)%vert(:,i,j)%y = 0.0_r8
          cvlist(ie)%vert(:,i,j)%z = 0.0_r8
        end do
      end do

      do j=-1,np+1
        do i=-1,np+1
          cv(i,j)%x = vertunpack(i,j,1)
          cv(i,j)%y = vertunpack(i,j,2)
          cv(i,j)%z = vertunpack(i,j,3)
        end do
      end do

      do j=-1,0
        do i=-1,0
          do k=1,mlt(swest)-1
            cv_sw(i,j,k)%x = sw(i,j,1,k)
            cv_sw(i,j,k)%y = sw(i,j,2,k)
            cv_sw(i,j,k)%z = sw(i,j,3,k)
          end do
        end do
      end do
      do j=-1,0
        do i=np,np+1
          do k=1,mlt(seast)-1
            cv_se(i,j,k)%x = se(i,j,1,k)
            cv_se(i,j,k)%y = se(i,j,2,k)
            cv_se(i,j,k)%z = se(i,j,3,k)
          end do
        end do
      end do
      do j=np,np+1
        do i=-1,0
          do k=1,mlt(nwest)-1
            cv_nw(i,j,k)%x = nw(i,j,1,k)
            cv_nw(i,j,k)%y = nw(i,j,2,k)
            cv_nw(i,j,k)%z = nw(i,j,3,k)
          end do
        end do
      end do
      do j=np,np+1
        do i=np,np+1
          do k=1,mlt(neast)-1
            cv_ne(i,j,k)%x = ne(i,j,1,k)
            cv_ne(i,j,k)%y = ne(i,j,2,k)
            cv_ne(i,j,k)%z = ne(i,j,3,k)
          end do
        end do
      end do

      do j=2,np-1
        do i=2,np-1
          ! internal vertex on Cubed sphere
          ! Here is the order:
          !
          ! 4NW <- 3NE
          !  |      ^
          !  v      |
          ! 1SW -> 2SE
          vert(1) = cv(i-1, j-1)
          vert(2) = cv(i  , j-1)
          vert(3) = cv(i  , j  )
          vert(4) = cv(i-1, j  )
          cvlist(ie)%vert(1:4,i,j) = vert(1:4)
          cvlist(ie)%nvert(i,j) = 4
          m=4
        end do
      end do

      do j=0,np,np
        do i=2,np-1
          vert(1) = cv(i-1, j-1)
          vert(2) = cv(i  , j-1)
          vert(3) = cv(i  , j  )
          vert(4) = cv(i  , j+1)
          vert(5) = cv(i-1, j+1)
          vert(6) = cv(i-1, j  )
          p = j
          if (p.eq.0) p=1
          cvlist(ie)%vert(1:6,i,p) = vert(1:6)
          cvlist(ie)%nvert(i,p) = 6
          m=6
        end do
      end do

      do j=2,np-1
        do i=0,np,np
          vert(1) = cv(i-1, j-1)
          vert(2) = cv(i  , j-1)
          vert(3) = cv(i+1, j-1)
          vert(4) = cv(i+1, j  )
          vert(5) = cv(i  , j  )
          vert(6) = cv(i-1, j  )
          o = i
          if (o.eq.0) o=1
          cvlist(ie)%vert(1:6,o,j) = vert(1:6)
          cvlist(ie)%nvert(o,j) = 6
          m=6
        end do
      end do
      do j=0,np,np
        do i=0,np,np
          m = 0
          vert(:)%x = 0
          vert(:)%y = 0
          vert(:)%z = 0
          if (i.eq.0.and.j.eq.0) then
            ! counterclockwise from lower right
            vert(m+1) = cv(i+1, j-1)  !     5       4
            vert(m+2) = cv(i+1, j  )  !  (-1,+1)  (0,+1)  (+1,+1)  3
            vert(m+3) = cv(i+1, j+1)  !
            vert(m+4) = cv(i  , j+1)  !  (-1, 0)  (i, j)  (+1, 0)  2
            vert(m+5) = cv(i-1, j+1)  !
            vert(m+6) = cv(i-1, j  )  !     X       X     (+1,-1)  1
            m = m + 6
            if (mlt(swest).ne.0) then
              vert(m+1) = cv(i-1, j-1)
              vert(m+2) = cv(i  , j-1)
              m = m+2
              do k=1,mlt(swest)-1   ! Bummer, toss in (-1,0) because transpose is undetectable
                vert(m+1) = cv_sw(i-1, j  , k)
                vert(m+2) = cv_sw(i-1, j-1, k)
                vert(m+3) = cv_sw(i  , j-1, k)
                m=m+3
              end do
            end if
          end if
          if (i.eq.np.and.j.eq.0) then
            if (mlt(seast).ne.0) then
              vert(m+1) = cv(i+1, j-1)
              vert(m+2) = cv(i+1, j  )
              m = m+2
              do k=1,mlt(seast)-1
                vert(m+1) = cv_se(i  , j-1, k)
                vert(m+2) = cv_se(i+1, j-1, k)
                vert(m+3) = cv_se(i+1, j  , k)
                m=m+3
              end do
            end if
            vert(m+1) = cv(i+1, j+1)
            vert(m+2) = cv(i  , j+1)
            vert(m+3) = cv(i-1, j+1)
            vert(m+4) = cv(i-1, j  )
            vert(m+5) = cv(i-1, j-1)
            vert(m+6) = cv(i  , j-1)
            m = m + 6
          end if
          if (i.eq.np.and.j.eq.np) then
            vert(1) = cv(i+1, j-1)
            vert(2) = cv(i+1, j  )
            m = m + 2
            if (mlt(neast).ne.0) then
              vert(m+1) = cv(i+1, j+1)
              vert(m+2) = cv(i  , j+1)
              m = m+2
              do k=1,mlt(neast)-1
                vert(m+1) = cv_ne(i+1, j  , k)
                vert(m+2) = cv_ne(i+1, j+1, k)
                vert(m+3) = cv_ne(i  , j+1, k)
                m=m+3
              end do
            end if
            vert(m+1) = cv(i-1, j+1)
            vert(m+2) = cv(i-1, j  )
            vert(m+3) = cv(i-1, j-1)
            vert(m+4) = cv(i  , j-1)
            m = m + 4
          end if
          if (i.eq.0.and.j.eq.np) then
            vert(m+1) = cv(i+1, j-1)
            vert(m+2) = cv(i+1, j  )
            vert(m+3) = cv(i+1, j+1)
            vert(m+4) = cv(i  , j+1)
            m = m + 4
            if (mlt(nwest).ne.0) then
              vert(m+1) = cv(i-1, j+1)
              vert(m+2) = cv(i-1, j  )
              m = m+2
              do k=1,mlt(nwest)-1
                vert(m+1) = cv_nw(i  , j+1, k)
                vert(m+2) = cv_nw(i-1, j+1, k)
                vert(m+3) = cv_nw(i-1, j  , k)
                m=m+3
              end do
            end if
            vert(m+1) = cv(i-1, j-1)
            vert(m+2) = cv(i  , j-1)
            m = m + 2
          end if
          o = i
          p = j
          if (o.eq.0) o=1
          if (p.eq.0) p=1
          m2=m
          if (8 < m) then
            m = SortNodes(vert, m2)
          end if
          if (m > nv_max) then
            call endrun("error: vert dimensioned too small")
          end if
          cvlist(ie)%vert(1:m,o,p) = vert(1:m)
          cvlist(ie)%nvert(o,p) = m
        end do
      end do
    end do
  end subroutine construct_cv_duel

  function SurfArea( cv, nvert ) result(area)

    type(cartesian3D_t), intent(in) :: cv(:)
    integer,             intent(in) :: nvert

    real(kind=r8)                   :: area, area1, area2, area3

    if (abs(nvert) == 3 ) then
      area2 = 0.0_r8
      area3 = 0.0_r8
      if (cv(1)%x == 0) then
        call sphere_tri_area(cv(2), cv(3), cv(4), area1)
      else if (cv(2)%x == 0) then
        call sphere_tri_area(cv(1), cv(3), cv(4), area1)
      else if (cv(3)%x == 0) then
        call sphere_tri_area(cv(1), cv(2), cv(4), area1)
      else if (cv(4)%x == 0) then
        call sphere_tri_area(cv(1), cv(2), cv(3), area1)
      else
        write(iulog, *) cv(1)%x, cv(1)%y
        write(iulog, *) cv(2)%x, cv(2)%y
        write(iulog, *) cv(3)%x, cv(3)%y
        write(iulog, *) cv(4)%x, cv(4)%y
        write(iulog, *) 'SurfArea error: should never happen'
        call shr_sys_flush(iulog)
        call endrun('SurfArea: invalid cv coordinates')
      end if
    else if (abs(nvert)==4) then
      call sphere_tri_area(cv(1), cv(2), cv(3), area1)
      call sphere_tri_area(cv(1), cv(3), cv(4), area2)
      area3 = 0.0_r8

    else if (abs(nvert)==5) then
      call sphere_tri_area(cv(1),cv(2),cv(3),area1)
      call sphere_tri_area(cv(1),cv(3),cv(4),area2)
      call sphere_tri_area(cv(1),cv(4),cv(5),area3)
    else
      call endrun('SurfArea: nvert > 5 not yet supported')
    end if
    area = area1 + area2 + area3
  end function SurfArea

  !   ^
  !   |dy  o
  !   |
  ! (x,y) ---->dx
  function SurfArea_dxdy(dx, dy, corner) result(integral)
    use quadrature_mod, only: quadrature_t

    real(r8),            intent(in) :: dx, dy
    type(cartesian2d_t), intent(in) :: corner
    real(r8)                        :: integral

    real(r8)                        :: alpha, beta, a1, a2, a3, a4

    ! cubed-sphere cell area, from Lauritzen & Nair MWR 2008
    ! central angles:
    ! cube face: -pi/4,-pi/4 -> pi/4,pi/4
    ! this formula gives 2   so normalize by 4pi/6 / 2 = pi/3
    alpha = corner%x
    beta  = corner%y
    a1 =  acos(-sin(alpha)*sin(beta))             !  2.094
    a2 = -acos(-sin(alpha+dx)*sin(beta) )         ! -1.047
    a3 =- acos(-sin(alpha)*sin(beta+dy) )         ! -1.047
    a4 =  acos(-sin(alpha+dx)*sin(beta+dy) )      !  2.094
    integral = (a1+a2+a3+a4)
    return
  end function SurfArea_dxdy

  function find_intersect(x1in, x2in, y1in, y2in) result(sect)

    type(cartesian2D_t), intent(in) :: x1in, x2in, y1in, y2in
    type(cartesian2D_t)             :: sect

    type(cartesian2D_t)             :: x, y, b, x1, x2, y1, y2
    real(kind=r8)                   :: s1, s2, detA

    !  x1 + (x2-x1)*s1  = y1 + (y2-y1)*s2
    ! b = y1-x1
    ! x=x2-x1
    ! y=y2-y1
    !  x s1 - y s2 = b
    !  x(1) s1 - y(1) s2 = b(1)
    !  x(2) s1 - y(2) s2 = b(2)
    !
    !  x(1) -y(1)   s1   =  b(1)        A s = b
    !  x(2) -y(2)   s2   =  b(2)
    !
    !  A2=  -y(2)  y(1)
    !       -x(2)  x(1)                s = A2 * b /detA

    ! convert to gnomonic
    x1%x = tan(x1in%x)
    x2%x = tan(x2in%x)
    y1%x = tan(y1in%x)
    y2%x = tan(y2in%x)
    x1%y = tan(x1in%y)
    x2%y = tan(x2in%y)
    y1%y = tan(y1in%y)
    y2%y = tan(y2in%y)

    x%x = x2%x-x1%x
    x%y = x2%y-x1%y
    y%x = y2%x-y1%x
    y%y = y2%y-y1%y
    b%x = y1%x-x1%x
    b%y = y1%y-x1%y

    detA = -x%x*y%y + x%y*y%x

    s1 =  (-y%y*b%x + y%x*b%y )/detA
    s2 =  (-x%y*b%x + x%x*b%y )/detA

    sect%x = x1%x + (x2%x-x1%x)*s1
    sect%y = x1%y + (x2%y-x1%y)*s1

    sect%x = (sect%x + y1%x + (y2%x-y1%x)*s2)/2
    sect%y = (sect%y + y1%y + (y2%y-y1%y)*s2)/2

    if (s1<0 .or. s1>1) then
      write(iulog, *) 'failed: intersection: ',s1,s2
      call shr_sys_flush(iulog)
      call endrun('find_intersect: intersection failure')
    end if

    ! convert back to equal angle:
    sect%x = atan(sect%x)
    sect%y = atan(sect%y)
  end function find_intersect

  subroutine pentagon_iteration(sq1,sq2,pent,asq1,asq2,apent,faceno,anorm)
    !               sq2
    !              4  3
    !              1  2
    !
    !    sq1       4  3
    !    2  1    5  pent
    !    3  4    1    2
    !
    !
    !   d/dt sq1(1)  =    (area(sq1)-asq1) * [ com(sq1)-sq1(1) ]
    !                    +(area(sq2)-asq2) * [ com(sq2)-sq1(1) ]
    !                    +(area(pent)-apent) * [ com(pent)-sq1(1) ]
    !
    !
    !
    type(cartesian2d_t), intent(inout) :: sq1(4), sq2(4), pent(5)
    real(r8),            intent(in)    :: asq1, asq2, apent, anorm
    integer,             intent(in)    :: faceno

    type(cartesian3D_t) :: sq1_3d(4), sq2_3d(4), pent_3d(5)
    real(r8)            :: isq1, isq2, ipent, diff1, diff2, diffp, err
    real(r8), parameter :: dt = .5_r8
    real(r8), parameter :: tol_pentagon_iteration = 1.0e-10_r8
    type(cartesian2d_t) :: sq1com, sq2com, pentcom, ds1, ds2
    integer             :: i, iter
    integer,  parameter :: iter_max = 10000

    ! compute center of mass:
    sq1com%x = sum(sq1(:)%x)/4
    sq1com%y = sum(sq1(:)%y)/4
    sq2com%x = sum(sq2(:)%x)/4
    sq2com%y = sum(sq2(:)%y)/4
    pentcom%x = sum(pent(:)%x)/5
    pentcom%y = sum(pent(:)%y)/5

    do i = 1, 4
      sq1_3d(i)=cubedsphere2cart(sq1(i),faceno  )
      sq2_3d(i)=cubedsphere2cart(sq2(i),faceno  )
      pent_3d(i)=cubedsphere2cart(pent(i),faceno  )
    end do
    pent_3d(5)=cubedsphere2cart(pent(5),faceno  )

    do iter = 1, iter_max
      isq1 = SurfArea(sq1_3d,4)
      isq2 = SurfArea(sq2_3d,4)
      ipent = SurfArea(pent_3d,5)

      !   d/dt sq1(1)  =    (area(sq1)-asq1) * [ com(sq1)-sq1(1) ]
      !                    +(area(sq2)-asq2) * [ com(sq2)-sq1(1) ]
      !                    +(area(pent)-apent) * [ com(pent)-sq1(1) ]
      !
      diff1 = (isq1-asq1)/anorm
      diff2 = (isq2-asq2)/anorm
      diffp = (ipent-apent)/anorm

      err = abs(diff1) + abs(diff2) + abs(diffp)
      if (err < tol_pentagon_iteration) exit
      if (mod(iter,1000) == 0) then
        write(iulog, '(i5,3e18.5)') iter, err
        call shr_sys_flush(iulog)
      end if

      ds1%x = diff1* ( sq1com%x - sq1(1)%x )
      ds1%y = diff1* ( sq1com%y - sq1(1)%y )
      ds1%x = ds1%x + diffp* ( pentcom%x - sq1(1)%x )
      ds1%y = ds1%y + diffp* ( pentcom%y - sq1(1)%y )

      ds2%x = diff2* ( sq2com%x - sq2(1)%x )
      ds2%y = diff2* ( sq2com%y - sq2(1)%y )
      ds2%x = ds2%x + diffp* ( pentcom%x - sq2(1)%x )
      ds2%y = ds2%y + diffp* ( pentcom%y - sq2(1)%y )

      sq1(1)%x = sq1(1)%x + dt*ds1%x
      sq1(1)%y = sq1(1)%y + dt*ds1%y
      sq2(1)%x = sq2(1)%x + dt*ds2%x
      sq2(1)%y = sq2(1)%y + dt*ds2%y
      pent(4)=sq2(1)
      pent(5)=sq1(1)
      sq1_3d(1)=cubedsphere2cart(sq1(1),faceno  )
      sq2_3d(1)=cubedsphere2cart(sq2(1),faceno  )
      pent_3d(4)=sq2_3d(1)
      pent_3d(5)=sq1_3d(1)
    end do
    if (iter >= iter_max) then
      write(iulog, *) 'pentagon iteration did not converge err=', err
      call shr_sys_flush(iulog)
    end if
  end subroutine pentagon_iteration

  subroutine InitControlVolumes_gll(elem, hybrid,nets,nete)
    use edge_mod,               only: freeedgebuffer
    use element_mod,            only: element_t,element_coordinates
    use hybrid_mod,             only: hybrid_t

    use quadrature_mod,         only: quadrature_t, gausslobatto
    use dimensions_mod,         only: nlev
    use cube_mod,               only: convert_gbl_index
    use coordinate_systems_mod, only: cart2cubedsphere_failsafe, cart2cubedsphere
    use coordinate_systems_mod, only: cube_face_number_from_sphere

    integer,         intent(in)         :: nets,nete
    type(element_t), intent(in), target :: elem(:)
    type(hybrid_t),  intent(in)         :: hybrid

    type(cartesian2d_t)     :: cartp_com(np,np)  ! center of mass
    type(cartesian2d_t)     :: cartp_nm1(0:np,0:np)
    real(r8)                :: delx_k,dely_k,sum_dbg,r
    integer                 :: i,j,ie,k,kptr,gllpts,nvert,k2,ie1,je1,face_no,kinsert
    integer                 :: iter,iter_max,i1,j1
    real(r8)                :: diff(np,np),diffy(np-1,np-1),diffx(np-1,np-1)
    real(r8)                :: dx,dy,a1(nets:nete),a2(nets:nete),d1(nets:nete),d1mid(nets:nete)
    real(r8)                :: d2,d1_global,d1_global_mid,sphere1,sphere2,diff2,diff3
    real(r8)                :: diff23,diff32,diff33,diff22
    real(r8)                :: gllnm1(0:np) !was longdouble_kind in HOMME
    type(cartesian2d_t)     :: corner,start,endd,cv_loc_2d(4,np,np),cvnew_loc_2d(4,np,np)
    type(cartesian3D_t)     :: cart,cv_loc_3d(nv_max,np,np)
    type(cartesian3D_t)     :: temp3d(nv_max)
    type(cartesian2d_t)     :: cartp2d(np,np)
    type(cartesian2d_t)     :: x1,x2,x3,x
    type(cartesian2d_t)     :: sq1(4),sq2(4),pent(5)
    type(cartesian3D_t)     :: x1_3d,x2_3d,x3_3d
    type(quadrature_t)      :: gll
    type(cartesian2d_t)     :: dir,dirsum
    type(spherical_polar_t) :: polar_tmp(0:np,0:np)
    real(r8)                :: rvert,area1,area2,ave,lat(4),lon(4)
    real(r8)                :: s,ds,triarea,triarea_target
    real(r8)                :: xp1,xm1,yp1,ym1,sumdiff
    real(r8)                :: tiny = 1e-11_r8,norm
    real(r8)                :: tol = 2.e-11_r8  ! convergece outer iteration
    real(r8)                :: tol_pentagons = 1.e-13_r8  ! convergece pentagon iteration

    ! area difference to trigger pentagons.
    ! if it is too small, we will have pentagons with 1 very short edges
    ! accuracy of surfarea() with very thin triangles seems poor (1e-11)
    ! ne=30  1e-3:  add 648 pentagons.  area ratio:  1.003
    ! ne=30  1e-4:  add 696 pentagons.  area ratio:  1.000004102
    ! ne=30  1e-5:  add 696 pentagons.  area ratio:  1.000004102
    ! ne=240 1e-4:  add 5688/ 345600 pentagons, area ratio: 1.0004
    ! ne=240 1e-5:  add 5736/ 345600 pentagons, area ratio: 1.000000078
    real(r8)                :: tol_use_pentagons=1.0e-5_r8
    logical                 :: Debug=.FALSE.,keep

    integer                 :: face1,face2,found,ie_max,movex,movey,moved,ii,kmax,kk
    integer                 :: nskip,npent
    integer                 :: nskipie(nets:nete), npentie(nets:nete)
    type(cartesian2d_t)     :: vert1_2d, vert_2d,vert2_2d
    type(cartesian3D_t)     :: vert1,vert2,vert_inserted(7)

    kmax=4

    gll = gausslobatto(np)
    ! mid point rule:
    do i=1,np-1
      gllnm1(i) = ( gll%points(i) + gll%points(i+1) ) /2
    end do
    ! check that gll(i) < gllnm1(i) < gll(i+1)
    do i=1,np-1
      if (gll%points(i) > gllnm1(i) .or. gllnm1(i) > gll%points(i+1)) then
        call endrun("InitControlVolumes_gll: CV and GLL points not interleaved")
      end if
    end do
    gllnm1(0)=-1
    gllnm1(np)=1

    ! MNL: dx and dy are no longer part of element_t
    !      but they are easily computed for the
    !      uniform case
    dx = pi/(2.0d0*dble(ne))
    dy = dx

    ! intialize local element dual grid, local element areas

    do ie=nets,nete

      call convert_gbl_index(elem(ie)%vertex%number,ie1,je1,face_no)
      start%x=-pi/4 + ie1*dx
      start%y=-pi/4 + je1*dy
      endd%x  =start%x + dx
      endd%y  =start%y + dy
      cartp_nm1(0:np,0:np) = element_coordinates(start,endd,gllnm1)
      cvlist(ie)%cartp_dual = cartp_nm1

      ! compute true area of element and SEM area
      a1(ie) = SurfArea_dxdy(dx,dy,elem(ie)%cartp(1,1))
      a2(ie) = sum(elem(ie)%spheremp(:,:))
      do i=1,np
        do j=1,np
          ! (gnomonic coordinate lines only), more accurate
          delx_k = cartp_nm1(i,j-1)%x - cartp_nm1(i-1,j-1)%x
          dely_k = cartp_nm1(i-1,j)%y - cartp_nm1(i-1,j-1)%y
          cvlist(ie)%vol(i,j) = SurfArea_dxdy(delx_k,dely_k,cartp_nm1(i-1,j-1))
        end do
      end do
      global_shared_buf(ie,1) = a1(ie)
      global_shared_buf(ie,2) = a2(ie)
    end do
    call wrap_repro_sum(nvars=2, comm=hybrid%par%comm)
    sphere1 = global_shared_sum(1)
    sphere2 = global_shared_sum(2)

    ! construct the global CV grid and global CV areas from the
    ! local dual grid (cvlist()%cart_dual) and local areas (cvlist()%vol)
    call construct_cv_gll(elem,hybrid,nets,nete)

    iter_max=2000
    if (iter_max>0) then
      ! areas computed from eleemnts on boundaries are from hexagons and pentagons
      ! compute new areas where all CVs are squares or triangles
      do ie=nets,nete
        do i=1,np
          do j=1,np
            ! ifort bug if we try this:
            ! area2 = surfarea(cvlist(ie)%vert(1:4,i,j),cvlist(ie)%nvert(i,j))
            cv_loc_3d(:,i,j)=cvlist(ie)%vert(:,i,j)
            area2 = surfarea(cv_loc_3d(:,i,j),cvlist(ie)%nvert(i,j))
            cvlist(ie)%totvol(i,j)=area2
          end do
        end do
      end do
    end if
    ! iteration over cvlist(ie)%totvol
    d1_global=0
    do iter=1,iter_max
      ie_max=-1
      do ie=nets,nete
        ! we want at each point, the gll_area = true_area
        ! but sum(gll_area) = a2 and sum(true_area)=a1
        ! so normalize so that: gll_area/a2 = true_area/a1, or gll_area = area*a2/a1

        ! requires more iterations, but the total volume within an
        ! element is always correct
        diff(:,:) = ( cvlist(ie)%vol(:,:) - elem(ie)%spheremp(:,:)*a1(ie)/a2(ie) )
        sumdiff=sum( cvlist(ie)%vol(:,:)) - a1(ie)
        diff(:,:) = diff(:,:)/(a1(ie)/(np*np))



        ! set boundary values (actually not used)
        cartp_nm1 = cvlist(ie)%cartp_dual(0:np,0:np)
        ! convert 9 cv corners in this element into cart_nm1 cubed-sphere coordiantes
        do i=1,np-1
          do j=1,np-1
            cartp_nm1(i,j) = cart2cubedsphere( cvlist(ie)%vert(3,i,j),elem(ie)%FaceNum  )
          end do
        end do
        ! compute center of mass of control volumes:
        ! todo: move points towards GLL points, not center of mass
        ! center of mass could send up a feedback with CV points!
        do i=1,np
          do j=1,np
            cart%x = sum( cvlist(ie)%vert(:,i,j)%x )/abs(cvlist(ie)%nvert(i,j))
            cart%y = sum( cvlist(ie)%vert(:,i,j)%y )/abs(cvlist(ie)%nvert(i,j))
            cart%z = sum( cvlist(ie)%vert(:,i,j)%z )/abs(cvlist(ie)%nvert(i,j))
            cartp_com(i,j) = cart2cubedsphere( cart,elem(ie)%FaceNum  )
          end do
        end do
        d2=0
        do i=1,np-1
          do j=1,np-1
            dirsum%x=0
            dirsum%y=0
            movex=1
            movey=1
            moved=0

            do i1=0,1
              do j1=0,1
                ! keep=.true. : .85/1.05
                ! corners only: .93/1.07
                ! corners and edges:  .89/1.11
                keep=.false.
                ! corner volumes
                if (i==1 .and. j==1) then
                  if (i1==0 .and. j1==0) keep=.true.
                  moved=1
                else if (i==np-1 .and. j==1) then
                  if (i1==1 .and. j1==0) keep=.true.
                  moved=-1
                else if (i==1 .and. j==np-1) then
                  if (i1==0 .and. j1==1) keep=.true.
                  moved=-1
                else if (i==np-1 .and. j==np-1) then
                  if (i1==1 .and. j1==1) keep=.true.
                  moved=1
                  ! edge volumes


                else if (i==1) then
                  if (i1==0) keep=.true.
                else if (i==np-1) then
                  if (i1==1) keep=.true.
                else if (j==1) then
                  if (j1==0) keep=.true.
                else if (j==np-1) then
                  if (j1==1) keep=.true.
                else
                  keep=.true.
                end if
                if (keep) then
                  ! error weighted direction towards center of mass of area
                  ! move towards grid point
                  dir%x =  (elem(ie)%cartp(i+i1,j+j1)%x - cartp_nm1(i,j)%x )*(abs(diff(i+i1,j+j1)))
                  dir%y =  (elem(ie)%cartp(i+i1,j+j1)%y - cartp_nm1(i,j)%y )*(abs(diff(i+i1,j+j1)))
                  if (moved==1) then
                    ! project onto (1,1)/sqrt(2)
                    dir%x = dir%x/sqrt(2d0) + dir%y/sqrt(2d0)
                    dir%y = dir%x
                  end if
                  if (moved==-1) then
                    ! project onto (-1,1)/sqrt(2)
                    dir%y = -dir%x/sqrt(2d0) + dir%y/sqrt(2d0)
                    dir%x = -dir%y
                  end if


                  if ( diff(i+i1,j+j1) > 0 ) then
                    ! this volume is too big, so move cv point towards grid center
                    ! weighted by length error
                    dirsum%x = dirsum%x + movex*dir%x
                    dirsum%y = dirsum%y + movey*dir%y
                  else
                    dirsum%x = dirsum%x - movex*dir%x
                    dirsum%y = dirsum%y - movey*dir%y
                  end if
                end if
              end do
            end do
            d2 = d2 + dirsum%x**2 + dirsum%y**2
            cartp_nm1(i,j)%x = cartp_nm1(i,j)%x + 0.25_r8*dirsum%x
            cartp_nm1(i,j)%y = cartp_nm1(i,j)%y + 0.25_r8*dirsum%y

          end do
        end do
        cvlist(ie)%cartp_dual(0:np,0:np) = cartp_nm1
        d2=sqrt(d2)

        d1(ie)=sqrt(sum(diff**2))

        d1mid(ie)=d1(ie)
        ! ignore center cv's:
        diff(2:3,2:3)=0
        d1mid(ie)=sqrt(sum(diff**2))

      end do  ! ie loop
      dx=maxval(d1)
      d1_global = ParallelMax(dx,hybrid)
      dx=maxval(d1mid)
      d1_global_mid = ParallelMax(dx,hybrid)
      if (mod(iter-1,250).eq.0) then
        if (hybrid%masterthread) write(iulog, *) iter,"max d1=",d1_global,d1_global_mid
      end if
      ! compute new global CV  (cvlist(ie)%vert from cvlist(ie)%cartp_dual).
      ! cvlist()%totarea incorrect since local volumes not computed above
      call construct_cv_gll(elem,hybrid,nets,nete)

      ! update totvol (area of multi-element cv)
      do ie=nets,nete
        do i=1,np
          do j=1,np
            ! ifort bug if we try this:
            ! area2 = surfarea(cvlist(ie)%vert(1:4,i,j),cvlist(ie)%nvert(i,j))
            cv_loc_3d(:,i,j)=cvlist(ie)%vert(:,i,j)
            area2 = surfarea(cv_loc_3d(:,i,j),cvlist(ie)%nvert(i,j))
            cvlist(ie)%totvol(i,j) = area2
            if (isnan(area2)) then
              write(iulog, *) 'ie,i,j',ie,i,j
              write(iulog, *) cvlist(ie)%nvert(i,j)
              write(iulog, *) cv_loc_3d(1,i,j)
              write(iulog, *) cv_loc_3d(2,i,j)
              write(iulog, *) cv_loc_3d(3,i,j)
              write(iulog, *) cv_loc_3d(4,i,j)
              call shr_sys_flush(iulog)
              call endrun('InitControlVolumes_gll: area = NaN')
            end if
          end do
        end do
      end do

      ! update %vol (local control volume within each element)
      do ie=nets,nete
        cartp2d = elem(ie)%cartp
        do i=1,np
          do j=1,np
            ! ifort bug if we try this:
            ! area2 = surfarea(cvlist(ie)%vert(1:4,i,j),cvlist(ie)%nvert(i,j))

            do ii=1,4
              !
              ! if we do not use _failsafe version of cart2cubedsphere code will fail with "-debug"
              !
              cv_loc_2d(ii,i,j) = cart2cubedsphere_failsafe( cvlist(ie)%vert(ii,i,j),elem(ie)%FaceNum  )
            end do
            if (i==1 .and. j==1) then
              cv_loc_2d(1,i,j)=cartp2d(i,j)
            end if
            if (i==np .and. j==1) then
              cv_loc_2d(2,i,j)=cartp2d(i,j)
            end if
            if (i==1 .and. j==np) then
              cv_loc_2d(4,i,j)=cartp2d(i,j)
            end if
            if (i==np .and. j==np) then
              cv_loc_2d(3,i,j)=cartp2d(i,j)
            end if


            cvnew_loc_2d(:,i,j)=cv_loc_2d(:,i,j)

            !
            ! 4NW <- 3NE
            !  |      ^
            !  v      |
            ! 1SW -> 2SE
            if (i==1) then
              ! replace points with x< elem(ie)%vert(i,j)%x
              if (cv_loc_2d(1,i,j)%x < cartp2d(i,j)%x) then
                cvnew_loc_2d(1,i,j) = find_intersect(&
                     cv_loc_2d(1,i,j), cv_loc_2d(2,i,j),&
                     elem(ie)%cartp(i,1),elem(ie)%cartp(i,np))
              end if
              if (cv_loc_2d(4,i,j)%x < cartp2d(i,j)%x) then
                cvnew_loc_2d(4,i,j) = find_intersect(&
                     cv_loc_2d(4,i,j), cv_loc_2d(3,i,j),&
                     elem(ie)%cartp(i,1),elem(ie)%cartp(i,np))
              end if
            end if

            if (i==np) then
              ! replace points with x> elem(ie)%vert(i,j)%x
              if (cv_loc_2d(2,i,j)%x > cartp2d(i,j)%x) then
                cvnew_loc_2d(2,i,j) = find_intersect(&
                     cv_loc_2d(1,i,j), cv_loc_2d(2,i,j),&
                     elem(ie)%cartp(i,1),elem(ie)%cartp(i,np))
              end if
              if (cv_loc_2d(3,i,j)%x > cartp2d(i,j)%x) then
                cvnew_loc_2d(3,i,j) = find_intersect(&
                     cv_loc_2d(4,i,j), cv_loc_2d(3,i,j),&
                     elem(ie)%cartp(i,1),elem(ie)%cartp(i,np))
              end if
            end if
            !
            ! 4NW <- 3NE
            !  |      ^
            !  v      |
            ! 1SW -> 2SE
            if (j==1) then
              ! replace points with y < elem(ie)%vert(i,j)%y
              if (cv_loc_2d(1,i,j)%y < cartp2d(i,j)%y) then
                cvnew_loc_2d(1,i,j) = find_intersect(&
                     cv_loc_2d(1,i,j), cv_loc_2d(4,i,j),&
                     elem(ie)%cartp(1,j),elem(ie)%cartp(np,j))
              end if
              if (cv_loc_2d(2,i,j)%y < cartp2d(i,j)%y) then
                cvnew_loc_2d(2,i,j) = find_intersect(&
                     cv_loc_2d(2,i,j), cv_loc_2d(3,i,j),&
                     elem(ie)%cartp(1,j),elem(ie)%cartp(np,j))
              end if
            end if
            if (j==np) then
              ! replace points with y > elem(ie)%vert(i,j)%y
              if (cv_loc_2d(4,i,j)%y > cartp2d(i,j)%y) then
                cvnew_loc_2d(4,i,j) = find_intersect(&
                     cv_loc_2d(1,i,j), cv_loc_2d(4,i,j),&
                     elem(ie)%cartp(1,j),elem(ie)%cartp(np,j))
              end if
              if (cv_loc_2d(3,i,j)%y > cartp2d(i,j)%y) then
                cvnew_loc_2d(3,i,j) = find_intersect(&
                     cv_loc_2d(2,i,j), cv_loc_2d(3,i,j),&
                     elem(ie)%cartp(1,j),elem(ie)%cartp(np,j))
              end if
            end if
            do ii=1,4
              cv_loc_3d(ii,i,j)=cubedsphere2cart(cvnew_loc_2d(ii,i,j),elem(ie)%FaceNum  )
            end do
            area2 = surfarea(cv_loc_3d(:,i,j),4)
            cvlist(ie)%vol(i,j) = area2
            if (isnan(area2)) then
              write(iulog, *) 'ie,i,j',ie,i,j
              write(iulog, *) cvlist(ie)%nvert(i,j)
              write(iulog, *) cv_loc_3d(1,i,j)
              write(iulog, *) cv_loc_3d(2,i,j)
              write(iulog, *) cv_loc_3d(3,i,j)
              write(iulog, *) cv_loc_3d(4,i,j)
              call shr_sys_flush(iulog)
              call endrun('InitControlVolumes_gll: area = NaN')
            end if
          end do
        end do
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( d1_global > 10.0_r8 .or. d1_global_mid < tol) then
        if (hybrid%masterthread) then
          write(iulog, *) 'first iteration stopping:'
          write(iulog, *) iter, "max error=", d1_global_mid
          call shr_sys_flush(iulog)
        end if
        exit
      end if
    end do  ! iteration loop

    kmax=5

    nskip=0
    npent=0
    nskipie(:) = 0
    npentie(:) = 0
    do ie=nets,nete
      diff = ( cvlist(ie)%vol(:,:) - elem(ie)%spheremp(:,:)*a1(ie)/a2(ie) )
      if ( maxval(abs(diff(2:3,2:3)))/a1(ie)  > tol_use_pentagons ) then
        npent=npent+1
        npentie(ie) = npentie(ie) + 1
        !
        ! 4NW <- 3NE
        !  |      ^
        !  v      |             23   33
        ! 1SW -> 2SE            22   32
        if (diff(2,2)>0 .and. diff(3,3)>0) then
          x1 = cart2cubedsphere( cvlist(ie)%vert(3,2,2),elem(ie)%FaceNum  )
          x2 = cart2cubedsphere( cvlist(ie)%vert(1,2,2),elem(ie)%FaceNum  )
          s = .99_r8
          x3%x = x2%x + (x1%x-x2%x)*s
          x3%y = x2%y + (x1%y-x2%y)*s

          sq1(1) = x3
          sq1(2) = cart2cubedsphere( cvlist(ie)%vert(4,2,2),elem(ie)%FaceNum  )
          sq1(3) = cart2cubedsphere( cvlist(ie)%vert(1,2,2),elem(ie)%FaceNum  )
          sq1(4) = cart2cubedsphere( cvlist(ie)%vert(2,2,2),elem(ie)%FaceNum  )

          x2 = cart2cubedsphere( cvlist(ie)%vert(3,3,3),elem(ie)%FaceNum  )
          s = .99_r8
          x3%x = x2%x + (x1%x-x2%x)*s
          x3%y = x2%y + (x1%y-x2%y)*s

          sq2(1) = x3
          sq2(2) = cart2cubedsphere( cvlist(ie)%vert(2,3,3),elem(ie)%FaceNum  )
          sq2(3) = cart2cubedsphere( cvlist(ie)%vert(3,3,3),elem(ie)%FaceNum  )
          sq2(4) = cart2cubedsphere( cvlist(ie)%vert(4,3,3),elem(ie)%FaceNum  )

          pent(1) = cart2cubedsphere( cvlist(ie)%vert(1,3,2),elem(ie)%FaceNum  )
          pent(2) = cart2cubedsphere( cvlist(ie)%vert(2,3,2),elem(ie)%FaceNum  )
          pent(3) = cart2cubedsphere( cvlist(ie)%vert(3,3,2),elem(ie)%FaceNum  )
          pent(4) = sq2(1)
          pent(5) = sq1(1)

          call pentagon_iteration(sq1,sq2,pent,&
               elem(ie)%spheremp(2,2)*a1(ie)/a2(ie), &
               elem(ie)%spheremp(3,3)*a1(ie)/a2(ie), &
               elem(ie)%spheremp(3,2)*a1(ie)/a2(ie),elem(ie)%FaceNum,a1(ie))

          x2_3d=cubedsphere2cart(sq1(1),elem(ie)%FaceNum  )
          x3_3d=cubedsphere2cart(sq2(1),elem(ie)%FaceNum  )

          cvlist(ie)%vert(3,2,2)=x2_3d
          cvlist(ie)%vert(1,3,3)=x3_3d

          cvlist(ie)%vert(5,2,3)=cvlist(ie)%vert(4,2,3)
          cvlist(ie)%vert(4,2,3)=cvlist(ie)%vert(3,2,3)
          cvlist(ie)%vert(2,2,3)=x2_3d
          cvlist(ie)%vert(3,2,3)=x3_3d

          cvlist(ie)%vert(5,3,2)=x2_3d
          cvlist(ie)%vert(4,3,2)=x3_3d

          cvlist(ie)%nvert(2,3)=sign(5,cvlist(ie)%nvert(2,3))
          cvlist(ie)%nvert(3,2)=sign(5,cvlist(ie)%nvert(3,2))
        else if (diff(2,3) >0 .and. diff(3,2)>0) then
          !
          ! 4NW <- 3NE
          !  |      ^
          !  v      |             23   33
          ! 1SW -> 2SE            22   32
          x1 = cart2cubedsphere( cvlist(ie)%vert(2,2,3),elem(ie)%FaceNum  )
          x2 = cart2cubedsphere( cvlist(ie)%vert(4,2,3),elem(ie)%FaceNum  )
          s = .99_r8
          x3%x = x2%x + (x1%x-x2%x)*s
          x3%y = x2%y + (x1%y-x2%y)*s

          sq1(1) = x3
          sq1(2) = cart2cubedsphere( cvlist(ie)%vert(3,2,3),elem(ie)%FaceNum  )
          sq1(3) = cart2cubedsphere( cvlist(ie)%vert(4,2,3),elem(ie)%FaceNum  )
          sq1(4) = cart2cubedsphere( cvlist(ie)%vert(1,2,3),elem(ie)%FaceNum  )

          x2 = cart2cubedsphere( cvlist(ie)%vert(2,3,2),elem(ie)%FaceNum  )
          s = .99_r8
          x3%x = x2%x + (x1%x-x2%x)*s
          x3%y = x2%y + (x1%y-x2%y)*s

          sq2(1) = x3
          sq2(2) = cart2cubedsphere( cvlist(ie)%vert(1,3,2),elem(ie)%FaceNum  )
          sq2(3) = cart2cubedsphere( cvlist(ie)%vert(2,3,2),elem(ie)%FaceNum  )
          sq2(4) = cart2cubedsphere( cvlist(ie)%vert(3,3,2),elem(ie)%FaceNum  )

          pent(1) = cart2cubedsphere( cvlist(ie)%vert(4,2,2),elem(ie)%FaceNum  )
          pent(2) = cart2cubedsphere( cvlist(ie)%vert(1,2,2),elem(ie)%FaceNum  )
          pent(3) = cart2cubedsphere( cvlist(ie)%vert(2,2,2),elem(ie)%FaceNum  )
          pent(4) = sq2(1)
          pent(5) = sq1(1)

          call pentagon_iteration(sq1,sq2,pent,&
               elem(ie)%spheremp(2,3)*a1(ie)/a2(ie), &
               elem(ie)%spheremp(3,2)*a1(ie)/a2(ie), &
               elem(ie)%spheremp(2,2)*a1(ie)/a2(ie),elem(ie)%FaceNum,a1(ie))

          x2_3d=cubedsphere2cart(sq1(1),elem(ie)%FaceNum  )
          x3_3d=cubedsphere2cart(sq2(1),elem(ie)%FaceNum  )

          cvlist(ie)%vert(2,2,3)=x2_3d

          cvlist(ie)%vert(4,3,2)=x3_3d

          cvlist(ie)%vert(5,2,2)=cvlist(ie)%vert(4,2,2)
          cvlist(ie)%vert(3,2,2)=x3_3d
          cvlist(ie)%vert(4,2,2)=x2_3d


          cvlist(ie)%vert(1,3,3)=x3_3d
          cvlist(ie)%vert(5,3,3)=x2_3d

          cvlist(ie)%nvert(3,3)=sign(5,cvlist(ie)%nvert(3,3))
          cvlist(ie)%nvert(2,2)=sign(5,cvlist(ie)%nvert(2,2))
        else
          if (hybrid%masterthread) then
            write(iulog, *) ie,'bad type'
            call shr_sys_flush(iulog)
          end if
          call endrun('InitControlVolumes_gll: bad type')
        end if
        ! recompute areas:
        do i=2,3
          do j=2,3
            nvert=abs(cvlist(ie)%nvert(i,j))
            temp3d(1:nvert)=cvlist(ie)%vert(1:nvert,i,j)
            cvlist(ie)%vol(i,j)=surfarea(temp3d,nvert)
            cvlist(ie)%totvol(i,j)=cvlist(ie)%vol(i,j)
          end do
        end do
      else
        !write(iulog, *) 'skipping pentagon procedure ie=',ie
        !write(iulog, *) 'maxval diff: ',maxval(abs(diff(2:3,2:3)))/a1(ie)
        nskip=nskip+1
        nskipie(ie) = nskipie(ie) + 1
      end if
      global_shared_buf(ie,1) = nskipie(ie)
      global_shared_buf(ie,2) = npentie(ie)
    end do
    call wrap_repro_sum(nvars=2, comm=hybrid%par%comm)
    nskip = global_shared_sum(1)
    npent = global_shared_sum(2)
    if (hybrid%masterthread) then
      write(*,'(a,i7,a,i7)') "no. elements where pentagons were added: ",npent,"/",npent+nskip
    end if

    ! compute output needed for SCRIP:  lat/lon coordinates, and for the
    ! control volume with only 3 corners, repeat the last point to make a
    ! degenerate quad.
    do ie=nets,nete
      do j=1,np
        do i=1,np
          cvlist(ie)%vert_latlon(:,i,j)%lat = 0._r8
          cvlist(ie)%vert_latlon(:,i,j)%lon = 0._r8
          do k = 1, kmax
            rvert = cvlist(ie)%vert(k,i,j)%x**2+cvlist(ie)%vert(k,i,j)%y**2+cvlist(ie)%vert(k,i,j)%z**2
            if(rvert > 0.9_r8) then
              cvlist(ie)%vert_latlon(k,i,j) = change_coordinates(cvlist(ie)%vert(k,i,j))
            else
              ! coordinates = 0, this corner was not set above because this point
              ! only has 3 cells (corner point) pick either neighbor to make a degenerate quad
              k2 = k - 1
              if (k2 == 0) then
                k2 = 2  ! can only happen for corner point with data in 2,3,4
              end if
              cvlist(ie)%vert_latlon(k,i,j) = change_coordinates(cvlist(ie)%vert(k2,i,j))
              cvlist(ie)%vert(k,i,j) = cvlist(ie)%vert(k2,i,j)
            end if
          end do
        end do
      end do
    end do
    ! Release memory
    if(hybrid%masterthread) then
      call freeedgebuffer(edge1)
    end if

    initialized=.true.
  end subroutine  InitControlVolumes_gll

  subroutine construct_cv_gll(elem,hybrid,nets,nete)
    !
    ! construct global dual grid from local element dual grid cvlist(ie)%cartp_dual(:,:)
    ! all control volumes will be squares or triangles (at cube corners)
    !
    ! 10/2009: added option to make hexagon control volumes at cube edges and corners
    !
    use bndry_mod,   only: bndry_exchange
    use element_mod, only: element_t
    use hybrid_mod,  only: hybrid_t
    use edge_mod,    only: edgeVpack, edgeVunpack, edgeVunpackVert

    type(element_t), intent(in), target :: elem(:)
    type(hybrid_t),  intent(in)         :: hybrid
    integer,         intent(in)         :: nets,nete
    !   local
    integer             :: i,j,k,ie,kptr,nvert,ie2
    logical             :: corner
    real(r8)            :: test(np,np,1),vertpack(np,np,3),rvert
    type(cartesian2d_t) :: vert(4)
    type(cartesian2d_t) :: cartp_nm1(0:np,0:np)

    test(:,:,:) = 0

    do ie=nets,nete
      ! now construct the dual grid

      cartp_nm1 = cvlist(ie)%cartp_dual

      do j=1,np
        do i=1,np
          cvlist(ie)%vert(:,i,j)%x = 0_r8
          cvlist(ie)%vert(:,i,j)%y = 0_r8
          cvlist(ie)%vert(:,i,j)%z = 0_r8
        end do
      end do

      ! interior

      do j=2,np-1
        do i=2,np-1

          ! internal vertex on Cubed sphere
          ! Here is the order:
          !
          ! 4NW <- 3NE
          !  |      ^
          !  v      |
          ! 1SW -> 2SE
          vert(1)%x = cartp_nm1(i-1,j-1)%x
          vert(1)%y = cartp_nm1(i-1,j-1)%y
          vert(2)%x = cartp_nm1(i  ,j-1)%x
          vert(2)%y = cartp_nm1(i  ,j-1)%y
          vert(3)%x = cartp_nm1(i  ,j  )%x
          vert(3)%y = cartp_nm1(i  ,j  )%y
          vert(4)%x = cartp_nm1(i-1,j  )%x
          vert(4)%y = cartp_nm1(i-1,j  )%y

          do k=1,4
            cvlist(ie)%vert(k,i,j) = cubedsphere2cart(vert(k),elem(ie)%FaceNum)
          end do
          cvlist(ie)%nvert(i,j) = 4

        end do
      end do

      ! Compute everything on the edges and then sum
      do i=2,np-1
        j=1
        !
        ! 4NW <- 3NE
        !  |      ^
        !  v      |
        ! 1SW -> 2SE
        !
        !
        ! only pack top two nodes.
        ! leave other data zero, filled in by edgeexchange
        cvlist(ie)%vert(4,i,j)%x = cvlist(ie)%vert(1,i,j+1)%x
        cvlist(ie)%vert(4,i,j)%y = cvlist(ie)%vert(1,i,j+1)%y
        cvlist(ie)%vert(4,i,j)%z = cvlist(ie)%vert(1,i,j+1)%z
        cvlist(ie)%vert(3,i,j)%x = cvlist(ie)%vert(2,i,j+1)%x
        cvlist(ie)%vert(3,i,j)%y = cvlist(ie)%vert(2,i,j+1)%y
        cvlist(ie)%vert(3,i,j)%z = cvlist(ie)%vert(2,i,j+1)%z


        j=np

        cvlist(ie)%vert(1,i,j)%x = cvlist(ie)%vert(4,i,j-1)%x
        cvlist(ie)%vert(1,i,j)%y = cvlist(ie)%vert(4,i,j-1)%y
        cvlist(ie)%vert(1,i,j)%z = cvlist(ie)%vert(4,i,j-1)%z
        cvlist(ie)%vert(2,i,j)%x = cvlist(ie)%vert(3,i,j-1)%x
        cvlist(ie)%vert(2,i,j)%y = cvlist(ie)%vert(3,i,j-1)%y
        cvlist(ie)%vert(2,i,j)%z = cvlist(ie)%vert(3,i,j-1)%z

      end do

      do j=2,np-1
        i=1

        cvlist(ie)%vert(2,i,j)%x = cvlist(ie)%vert(1,i+1,j)%x
        cvlist(ie)%vert(2,i,j)%y = cvlist(ie)%vert(1,i+1,j)%y
        cvlist(ie)%vert(2,i,j)%z = cvlist(ie)%vert(1,i+1,j)%z
        cvlist(ie)%vert(3,i,j)%x = cvlist(ie)%vert(4,i+1,j)%x
        cvlist(ie)%vert(3,i,j)%y = cvlist(ie)%vert(4,i+1,j)%y
        cvlist(ie)%vert(3,i,j)%z = cvlist(ie)%vert(4,i+1,j)%z

        i=np

        cvlist(ie)%vert(4,i,j)%x = cvlist(ie)%vert(3,i-1,j)%x
        cvlist(ie)%vert(4,i,j)%y = cvlist(ie)%vert(3,i-1,j)%y
        cvlist(ie)%vert(4,i,j)%z = cvlist(ie)%vert(3,i-1,j)%z
        cvlist(ie)%vert(1,i,j)%x = cvlist(ie)%vert(2,i-1,j)%x
        cvlist(ie)%vert(1,i,j)%y = cvlist(ie)%vert(2,i-1,j)%y
        cvlist(ie)%vert(1,i,j)%z = cvlist(ie)%vert(2,i-1,j)%z

      end do

      ! Corners
      ! SW
      cvlist(ie)%vert(3,1,1)%x = cvlist(ie)%vert(1,2,2)%x
      cvlist(ie)%vert(3,1,1)%y = cvlist(ie)%vert(1,2,2)%y
      cvlist(ie)%vert(3,1,1)%z = cvlist(ie)%vert(1,2,2)%z

      ! SE
      cvlist(ie)%vert(4,np,1)%x = cvlist(ie)%vert(2,np-1,2)%x
      cvlist(ie)%vert(4,np,1)%y = cvlist(ie)%vert(2,np-1,2)%y
      cvlist(ie)%vert(4,np,1)%z = cvlist(ie)%vert(2,np-1,2)%z

      ! NE
      cvlist(ie)%vert(1,np,np)%x = cvlist(ie)%vert(3,np-1,np-1)%x
      cvlist(ie)%vert(1,np,np)%y = cvlist(ie)%vert(3,np-1,np-1)%y
      cvlist(ie)%vert(1,np,np)%z = cvlist(ie)%vert(3,np-1,np-1)%z

      ! NW
      cvlist(ie)%vert(2,1,np)%x = cvlist(ie)%vert(4,2,np-1)%x
      cvlist(ie)%vert(2,1,np)%y = cvlist(ie)%vert(4,2,np-1)%y
      cvlist(ie)%vert(2,1,np)%z = cvlist(ie)%vert(4,2,np-1)%z

      kptr=0
      test(:,:,1) = cvlist(ie)%vol(:,:)
      call edgeVpack(edge1,test(1,1,1),1,kptr,ie)

      cvlist(ie)%invvol(:,:) = cvlist(ie)%vol(:,:)

    end do  ! loop over NE

    call bndry_exchange(hybrid,edge1)

    do ie=nets,nete
      kptr=0
      call edgeVunpack(edge1, cvlist(ie)%invvol(1,1),1, kptr, ie)
      cvlist(ie)%totvol(:,:)=cvlist(ie)%invvol(:,:)
      cvlist(ie)%invvol(:,:)=1.0_r8/cvlist(ie)%invvol(:,:)
    end do

    ! Create the polygon at the edges of the element


    if(.NOT.(MODULO(np,2)==0)) then
      call endrun("surfaces_mod: NV odd not implemented")
    end if
    vertpack = 0
    do ie=nets,nete
      ! Special messed up copy
      !
      !ASC should be replaced by a edgepack
      ! S+N
      do i=1,np/2
        j=1
        vertpack(i,j,1) = cvlist(ie)%vert(3,i,j)%x
        vertpack(i,j,2) = cvlist(ie)%vert(3,i,j)%y
        vertpack(i,j,3) = cvlist(ie)%vert(3,i,j)%z
        j=np
        vertpack(i,j,1) = cvlist(ie)%vert(2,i,j)%x
        vertpack(i,j,2) = cvlist(ie)%vert(2,i,j)%y
        vertpack(i,j,3) = cvlist(ie)%vert(2,i,j)%z
      end do

      do i=np/2+1,np
        j=1
        vertpack(i,j,1) = cvlist(ie)%vert(4,i,j)%x
        vertpack(i,j,2) = cvlist(ie)%vert(4,i,j)%y
        vertpack(i,j,3) = cvlist(ie)%vert(4,i,j)%z
        j=np
        vertpack(i,j,1) = cvlist(ie)%vert(1,i,j)%x
        vertpack(i,j,2) = cvlist(ie)%vert(1,i,j)%y
        vertpack(i,j,3) = cvlist(ie)%vert(1,i,j)%z
      end do

      ! E+W
      do j=2,np/2
        i=1
        vertpack(i,j,1) = cvlist(ie)%vert(3,i,j)%x
        vertpack(i,j,2) = cvlist(ie)%vert(3,i,j)%y
        vertpack(i,j,3) = cvlist(ie)%vert(3,i,j)%z
        i=np
        vertpack(i,j,1) = cvlist(ie)%vert(4,i,j)%x
        vertpack(i,j,2) = cvlist(ie)%vert(4,i,j)%y
        vertpack(i,j,3) = cvlist(ie)%vert(4,i,j)%z
      end do

      do j=np/2+1,np-1
        i=1
        vertpack(i,j,1) = cvlist(ie)%vert(2,i,j)%x
        vertpack(i,j,2) = cvlist(ie)%vert(2,i,j)%y
        vertpack(i,j,3) = cvlist(ie)%vert(2,i,j)%z
        i=np
        vertpack(i,j,1) = cvlist(ie)%vert(1,i,j)%x
        vertpack(i,j,2) = cvlist(ie)%vert(1,i,j)%y
        vertpack(i,j,3) = cvlist(ie)%vert(1,i,j)%z
      end do

      do j=2,np-1
        do i=2,np-1
          vertpack(i,j,1) =0_r8
          vertpack(i,j,2) =0_r8
          vertpack(i,j,3) =0_r8
        end do
      end do

      kptr=0
      call edgeVpack(edge1,vertpack,3,kptr,ie)
    end do

    call bndry_exchange(hybrid,edge1)

    do ie=nets,nete
      kptr=0
      call edgeVunpackVert(edge1, cvlist(ie)%vert,ie)
      ! Count and orient vert array
      ! nvert is an integer: -4,-3,3,4
      ! Positive: 1->2->3->4 is counter clockwise on the sphere
      ! Negative: clockwise orientation
      do j=1,np
        do i=1,np
          nvert=0
          do k=1,4
            rvert = cvlist(ie)%vert(k,i,j)%x**2+cvlist(ie)%vert(k,i,j)%y**2+cvlist(ie)%vert(k,i,j)%z**2
            if(rvert>0.9_r8)nvert=nvert+1
          end do
          if(.NOT.Orientation(cvlist(ie)%vert(:,i,j),elem(ie)%FaceNum))nvert=-nvert
          cvlist(ie)%nvert(i,j) = nvert
          corner = ( ((i==1) .and. (j==1)) .or. &
               ((i==1) .and. (j==np)) .or. &
               ((i==np) .and. (j==1)) .or. &
               ((i==np) .and. (j==np)) )
          if (abs(nvert)/=4) then
            if (abs(nvert)/=3) then
              write(iulog, *) 'i,j,nvert=',i,j,nvert
              call shr_sys_flush(iulog)
              call endrun('construct_cv_gll: bad value of nvert')
            end if
            if (.not. corner) then
              write(iulog, *) 'non-corner node with only 3 verticies'
              write(iulog, *) 'ie,i,j,nvert,corner=',ie,i,j,nvert,corner
              write(iulog, *) cvlist(ie)%vert(1,i,j)%x
              write(iulog, *) cvlist(ie)%vert(2,i,j)%x
              write(iulog, *) cvlist(ie)%vert(3,i,j)%x
              write(iulog, *) cvlist(ie)%vert(4,i,j)%x
              !write(iulog, *) 'dual:'
              !do ie2=nets,nete
              !   write(iulog, *) ie2,maxval(cvlist(ie2)%cartp_dual(:,:)%x)
              !   write(iulog, *) ie2,maxval(cvlist(ie2)%cartp_dual(:,:)%y)
              !end do
              call shr_sys_flush(iulog)
              call endrun('construct_cv_gll: corner point should have nvert=3')
            end if
            ! nvert=3.  we are at a cube corner.  One of the control volume
            ! nodes from the 'missing' corner element should be all zeros:
            if (cvlist(ie)%vert(1,i,j)%x==0) then
              ! ok
            else if (cvlist(ie)%vert(2,i,j)%x==0) then
              ! ok
            else if (cvlist(ie)%vert(3,i,j)%x==0) then
              ! ok
            else if (cvlist(ie)%vert(4,i,j)%x==0) then
              ! ok
            else
              write(iulog, *) 'cube corner node with 4 neighbors'
              write(iulog, *) 'ie,i,j,nvert,corner=',ie,i,j,nvert,corner
              write(iulog, *) cvlist(ie)%vert(1,i,j)%x
              write(iulog, *) cvlist(ie)%vert(2,i,j)%x
              write(iulog, *) cvlist(ie)%vert(3,i,j)%x
              write(iulog, *) cvlist(ie)%vert(4,i,j)%x
              call shr_sys_flush(iulog)
              call endrun('construct_cv_gll: control volume at cube corner should be a triangle')
            end if

          end if
        end do
      end do
    end do
  end subroutine  construct_cv_gll

  logical function Orientation(v, FaceNum) result(orient)

    type(cartesian3d_t), intent(in) :: v(3)
    integer,             intent(in) :: FaceNum

    type(cartesian3D_t)             :: v12, v23
    real(r8)                        :: test, cart(3,3)

    orient = .FALSE.

    if ((FaceNum == 5).OR.(FaceNum == 6)) then

      cart(1,1) = v(1)%x
      cart(2,1) = v(1)%y
      cart(3,1) = v(1)%z

      cart(1,2) = v(2)%x
      cart(2,2) = v(2)%y
      cart(3,2) = v(2)%z

      cart(1,3) = v(3)%x
      cart(2,3) = v(3)%y
      cart(3,3) = v(3)%z

      v12%x = cart(1,2) - cart(1,1)
      v12%y = cart(2,2) - cart(2,1)
      v12%z = cart(3,2) - cart(3,1)

      v23%x = cart(1,3) - cart(1,2)
      v23%y = cart(2,3) - cart(2,2)
      v23%z = cart(3,3) - cart(3,2)

      test = (v12%y*v23%z - v12%z*v23%y)*v12%x &
           - (v12%x*v23%z - v12%z*v23%x)*v12%y &
           + (v12%x*v23%y - v12%y*v23%x)*v12%z

      if (test > 0_r8)then
        orient=.TRUE.
      end if

    else
      orient=.TRUE.
    end if

  end function Orientation

  subroutine VerifVolumes(elem, hybrid,nets,nete)
    use hybrid_mod,  only: hybrid_t
    use element_mod, only: element_t

    type(element_t), intent(in) :: elem(:)
    integer,         intent(in) :: nets,nete
    type(hybrid_t),  intent(in) :: hybrid

    real(r8)          :: psum,ptot,Vol_tmp(1),corr,maxelem_variation
    real(r8)          :: vol(np,np,nets:nete),r,rmin,rmax,a1,a2,locmin,locmax,emin,emax,dx,dy
    integer           :: i,j,ie,kptr,face

    real(r8), pointer :: locvol(:,:)

    dx = pi/(2.0d0*dble(ne))
    dy = dx

    if(.not. initialized) then
      call endrun('VerifyVolumes: Attempt to use volumes prior to initializing')
    end if
    rmin=2
    rmax=0
    maxelem_variation=0
    do ie=nets,nete
      locvol => GetVolume(ie)
      locmin = minval(locvol(:,:)*elem(ie)%rspheremp(:,:))
      locmax = maxval(locvol(:,:)*elem(ie)%rspheremp(:,:))
      rmin = min(rmin,locmin)
      rmax = max(rmax,locmax)

      if (locmax > 1.01_r8) then
        write(iulog, *) 'locmin(:,i)=',ie,locvol(1,1),1/elem(ie)%rspheremp(1,1)
      end if


      if (locmax-locmin > maxelem_variation) then
        maxelem_variation = locmax-locmin
        emin=locmin
        emax=locmax
      end if
    end do
    rmin = ParallelMin(rmin,hybrid)
    rmax = ParallelMax(rmax,hybrid)
    if(hybrid%masterthread) then
      write(iulog,'(a,2e14.7)') "Min/max ratio between spherical and GLL area:",rmin,rmax
    end if
    if (maxelem_variation == ParallelMax(maxelem_variation,hybrid) ) then
      write(iulog,'(a,2e14.7)') "Min/max ratio element with largest variation:",emin,emax
    end if
    call shr_sys_flush(iulog)

    rmin=2
    rmax=0
    do ie=nets,nete
      a1 = SurfArea_dxdy(dx,dy,elem(ie)%cartp(1,1))
      a2 = sum(elem(ie)%spheremp(:,:))
      r=a1/a2
      rmin = min(r,rmin)
      rmax = max(r,rmax)
    end do
    rmin = ParallelMin(rmin,hybrid)
    rmax = ParallelMax(rmax,hybrid)
    if(hybrid%masterthread) then
      write(*,'(a,2f12.9)') "Min/max ratio spherical and GLL element area:",rmin,rmax
    end if

    do ie=nets,nete
      global_shared_buf(ie,1:6) = 0.d0
      face = elem(ie)%FaceNum
      locvol => GetVolumeLocal(ie)
      do j=1,np
        do i=1,np
          global_shared_buf(ie,face) = global_shared_buf(ie,face) + locvol(i,j)
        end do
      end do
    end do
    call wrap_repro_sum(nvars=6, comm=hybrid%par%comm)

    ptot=0_r8
    do face=1,6
      red_sum%buf(1) = global_shared_sum(face)
      psum = red_sum%buf(1)

      ptot = ptot + psum

      if(hybrid%masterthread) then
        write(*,'(a,i2,a,2e23.15)') "cube face:",face," : SURFACE FV =",&
             6_r8*psum/(4_r8 * pi), &
             6_r8*psum/(4_r8 * pi)-1
      end if
    end do

    if(hybrid%masterthread) then
      write(iulog, *) "SURFACE FV (total)= ", ptot/(4_r8 * pi)
    end if

  end subroutine VerifVolumes

end module comp_gll_ctr_vol
