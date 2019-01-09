! Separate dynamics and physics grids

module dp_mapping
  use dimensions_mod,         only: np, npsq, fv_nphys
  use shr_kind_mod,           only: r8=>shr_kind_r8, shr_kind_cl
  use coordinate_systems_mod, only: spherical_polar_t
  use shr_const_mod,          only: pi => shr_const_pi
  use fvm_control_volume_mod, only: fvm_struct

  implicit none
  private
  save

  public :: dp_init
  public :: dp_reoorder
  public :: dp_write

  ! Total number of physics points per spectral element
  ! no physgrid:   nphys_pts = npsq   (physics on GLL grid)
  !    physgrid:   nphys_pts = nphys2 (physics on CSLAM grid)
  ! Value is set when se_fv_nphys namelist variable is read
  integer,            public :: nphys_pts = npsq

  ! NOTE:  dp_gid() is in space filling curve rank order
  !        all other global arrays are in block id (global id) order
  !
  !        dp_gid() is used to re-order data collected on root via mpi_gatherv
  !          into block id ordering
  !
  ! j=dp_gid(i)   i = element space filling curve rank
  !               j = element global id = block id = history file ordering
  !
  integer,        allocatable,dimension(:)         :: dp_gid  ! NE=240, integer*4 = 1.3MB
  integer, public,allocatable,dimension(:)         :: dp_owner

  real (r8),public,allocatable,dimension(:,:,:)  :: weights_all_fvm2phys
  integer  ,public,allocatable,dimension(:,:,:)  :: weights_eul_index_all_fvm2phys,weights_lgr_index_all_fvm2phys
  real (r8),public,allocatable,dimension(:,:,:)  :: weights_all_phys2fvm
  integer  ,public,allocatable,dimension(:,:,:)  :: weights_eul_index_all_phys2fvm,weights_lgr_index_all_phys2fvm
  integer  ,public,allocatable,dimension(:)      :: jall_fvm2phys,jall_phys2fvm
  integer  ,public                               :: num_weights_fvm2phys,num_weights_phys2fvm


contains
  subroutine dp_init(elem,fvm)
    use cam_abortutils, only: endrun
    use dimensions_mod, only: nelemd,nc,irecons_tracer
    use element_mod, only: element_t
    use spmd_utils, only: masterproc
    use cam_logfile, only: iulog
    use thread_mod, only: horz_num_threads

    implicit none
    type(element_t)  , dimension(nelemd), intent(in) :: elem
    type (fvm_struct), dimension(nelemd), intent(in) :: fvm

    num_weights_phys2fvm = 0
    num_weights_fvm2phys = 0
    if (fv_nphys>0) then
      num_weights_phys2fvm = (nc+fv_nphys)**2
      num_weights_fvm2phys = (nc+fv_nphys)**2

       allocate(weights_all_fvm2phys(num_weights_fvm2phys,irecons_tracer,nelemd))
       allocate(weights_eul_index_all_fvm2phys(num_weights_fvm2phys,2,nelemd))
       allocate(weights_lgr_index_all_fvm2phys(num_weights_fvm2phys,2,nelemd))

       allocate(weights_all_phys2fvm(num_weights_phys2fvm,irecons_tracer,nelemd))
       allocate(weights_eul_index_all_phys2fvm(num_weights_phys2fvm,2,nelemd))
       allocate(weights_lgr_index_all_phys2fvm(num_weights_phys2fvm,2,nelemd))
       allocate(jall_fvm2phys(nelemd))
       allocate(jall_phys2fvm(nelemd))

       call fvm2phys_init(elem,fvm,nc,fv_nphys,irecons_tracer,&
            weights_all_fvm2phys,weights_eul_index_all_fvm2phys,weights_lgr_index_all_fvm2phys,&
            weights_all_phys2fvm,weights_eul_index_all_phys2fvm,weights_lgr_index_all_phys2fvm,&
            jall_fvm2phys,jall_phys2fvm)

      call dp_replicated_init(elem)

      if (masterproc) then
        write(iulog, *) 'dp_init: Initialized phys2fvm/fvm2phys mapping vars'
      end if

    end if
  end subroutine dp_init

  subroutine dp_reoorder(before,after)
    use cam_abortutils, only: endrun
    use dimensions_mod, only: nelem
    !XXgoldyXX
    use cam_logfile,    only: iulog
    use spmd_utils,     only: masterproc
    use shr_sys_mod,    only: shr_sys_flush
    !XXgoldyXX
    implicit none
    real(r8), dimension(fv_nphys*fv_nphys,*), intent(in)  :: before
    real(r8), dimension(fv_nphys*fv_nphys,*), intent(out) :: after
    integer :: ie

    ! begin
    do ie = 1,nelem
    !XXgoldyXX
       if (dp_gid(ie) < 0) then
          if (masterproc) then
             write(iulog,*) 'ie =',ie,', dp_gid(ie) =',dp_gid(ie)
             call shr_sys_flush(iulog)
          end if
          call endrun('Bad element remap in dp_reoorder')
       end if
    !XXgoldyXX
       after(:,dp_gid(ie)) = before(:,ie)
    end do
  end subroutine dp_reoorder

  !!!

  subroutine dp_replicated_init(elem)
    use dimensions_mod, only: nelem, nelemd
    use element_mod,    only: element_t
    use cam_abortutils, only: endrun
    use spmd_utils,     only: masterproc, masterprocid, npes
    use spmd_utils,     only: mpicom, mpi_integer

    implicit none
    type(element_t),dimension(nelemd),intent(in) :: elem

    integer                          :: i,j,ierror
    integer,dimension(nelemd)        :: lgid
    integer,dimension(:),allocatable :: displs,recvcount

    ! begin

    allocate(displs(npes))
    allocate(dp_gid(nelem))
    allocate(recvcount(npes))
    call mpi_gather(nelemd, 1, mpi_integer, recvcount, 1, mpi_integer,        &
         masterprocid, mpicom, ierror)
    lgid(:) = elem(:)%globalid
    if (masterproc) then
      displs(1) = 0
      do i = 2,npes
        displs(i) = displs(i-1)+recvcount(i-1)
      end do
    end if
    call mpi_gatherv(lgid, nelemd, mpi_integer, dp_gid, recvcount, displs,    &
         mpi_integer, masterprocid, mpicom, ierror)
    if (masterproc) then
      allocate(dp_owner(nelem))
      dp_owner(:) = -1
      do i = 1,npes
        do j = displs(i)+1,displs(i)+recvcount(i)
          dp_owner(dp_gid(j)) = i-1
        end do
      end do
    end if
    deallocate(displs)
    deallocate(recvcount)
    ! minimize global memory use
    call mpi_barrier(mpicom,ierror)
    if (.not.masterproc) then
      allocate(dp_owner(nelem))
    end if
    call mpi_bcast(dp_gid,nelem,mpi_integer,masterprocid,mpicom,ierror)
    call mpi_bcast(dp_owner,nelem,mpi_integer,masterprocid,mpicom,ierror)
  end subroutine dp_replicated_init

  !!!


  !!!

  subroutine dp_write(elem, fvm, grid_format, filename_in)
    use cam_abortutils,         only: endrun
    use dimensions_mod,         only: nelem, nelemd
    use element_mod,            only: element_t
    use netcdf,                 only: nf90_create, nf90_close, nf90_enddef
    use netcdf,                 only: nf90_def_dim, nf90_def_var, nf90_put_var
    use netcdf,                 only: nf90_double, nf90_int, nf90_put_att
    use netcdf,                 only: nf90_noerr, nf90_strerror, nf90_clobber
    use spmd_utils,             only: masterproc, masterprocid, mpicom, npes
    use spmd_utils,             only: mpi_integer, mpi_real8
    use cam_logfile,            only: iulog
    use shr_sys_mod,            only: shr_sys_flush
    use dimensions_mod,         only: ne
    use coordinate_systems_mod, only: cart2spherical
    
    ! Inputs
    type(element_t),   intent(in) :: elem(:)
    type (fvm_struct), intent(in) :: fvm(:)
    character(len=*),  intent(in) :: grid_format
    character(len=*),  intent(in) :: filename_in
    
    real(r8), parameter :: rad2deg = 180._r8/pi
    
    ! Local variables
    integer           :: i, ie, ierror, j, status, ivtx
    integer           :: grid_corners_id, grid_rank_id, grid_size_id
    character(len=256) :: errormsg
    character(len=shr_kind_cl) :: filename
    integer           :: ncid
    integer           :: grid_dims_id, grid_area_id, grid_center_lat_id
    integer           :: grid_center_lon_id, grid_corner_lat_id
    integer           :: grid_corner_lon_id, grid_imask_id
    integer           :: gridsize
    integer           :: IOrootID
    logical           :: IOroot
    integer,allocatable,dimension(:) :: displs,recvcount

    real(r8), dimension(fv_nphys, fv_nphys, nelemd, 4, 2)  :: corners
    real(r8), dimension(fv_nphys, fv_nphys, nelemd)        :: lwork
    real(r8), allocatable, dimension(:)                    :: recvbuf
    real(r8), allocatable, dimension(:,:)                  :: gwork
    real(r8)                                               :: x, y
    type (spherical_polar_t)                               :: sphere

    ! begin

    !! Check to see if we are doing grid output
    if (trim(grid_format) == "no") then
      if (masterproc) then
        write(iulog, *) 'dp_write: Not writing phys_grid file.'
      end if
      return
    else if (trim(grid_format) /= 'SCRIP') then
      if (masterproc) then
        write(errormsg, *) 'dp_write: ERROR, bad value for se_write_grid, ',&
             trim(grid_format)
        call endrun(errormsg)
      end if
    end if
    
    ! Create the NetCDF file
    if (len_trim(filename_in) == 0) then
      write(filename, '(3(a,i0),3a)') "ne", ne, "np", np, ".pg", fv_nphys,    &
           "_", trim(grid_format), ".nc"
    else
      filename = trim(filename_in)
    end if
    status = nf90_create(trim(filename), nf90_clobber, ncid)
    if (status /= nf90_noerr) then
      call endrun("dp_write: "//trim(nf90_strerror(status)))
    end if
    
    ! PIO_put_var puts from its root node, find that (so we do our work there)
    IOrootID = masterprocid
    IOroot = masterproc
    
    ! Allocate workspace and calculate PE displacement information
    if (IOroot) then
      allocate(displs(npes))
      allocate(recvcount(npes))
    else
      allocate(displs(0))
      allocate(recvcount(0))
    end if
    gridsize = nelem * fv_nphys*fv_nphys
    if(masterproc) then
      write(iulog, *) 'Writing physics SCRIP grid file: ', trim(filename)
      write(iulog, '(a,i7,a,i3)') 'nelem = ', nelem, ', fv_nphys = ', fv_nphys
      call shr_sys_flush(iulog)
    end if
    call mpi_gather(nelemd*fv_nphys*fv_nphys, 1, mpi_integer, recvcount, 1,           &
         mpi_integer, IOrootID, mpicom, ierror)

    if (IOroot) then
      displs(1) = 0
      do i = 2, npes
        displs(i) = displs(i-1)+recvcount(i-1)
      end do
      allocate(recvbuf(gridsize))
    else
      allocate(recvbuf(0))
    end if
    allocate(gwork(4, gridsize))

    if (IOroot) then
      ! Define the horizontal grid dimensions for SCRIP output
      status = nf90_def_dim(ncid, "grid_corners", 4,        grid_corners_id)
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining dimension, grid_corners'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if
      status = nf90_def_dim(ncid, "grid_rank",    1,        grid_rank_id)
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining dimension, grid_rank'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if
      status = nf90_def_dim(ncid, "grid_size",    gridsize, grid_size_id)
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining dimension, grid_size'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      ! Define the coordinate variables
      status = nf90_def_var(ncid, "grid_dims", nf90_int, (/grid_rank_id/),  &
           grid_dims_id)
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining variable grid_dims'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      status = nf90_def_var(ncid, "grid_area", nf90_double,                 &
           (/grid_size_id/), grid_area_id)
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining variable grid_area'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      status = nf90_put_att(ncid, grid_area_id, "units", "radians^2")
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining attributes for grid_area'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      status = nf90_put_att(ncid, grid_area_id, "long_name", "area weights")
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining attributes for grid_area'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      status = nf90_def_var(ncid, "grid_center_lat", nf90_double,           &
           (/grid_size_id/), grid_center_lat_id)
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining variable grid_center_lat'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      status = nf90_put_att(ncid, grid_center_lat_id, "units", "degrees")
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining attributes for grid_center_lat'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      status = nf90_def_var(ncid, "grid_center_lon", nf90_double,           &
           (/grid_size_id/), grid_center_lon_id)
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining variable grid_center_lon'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      status = nf90_put_att(ncid, grid_center_lon_id, "units", "degrees")
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining attributes for grid_center_lon'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      status = nf90_def_var(ncid, "grid_corner_lat", nf90_double,           &
           (/grid_corners_id, grid_size_id/), grid_corner_lat_id)
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining grid_corner_lat'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      status = nf90_put_att(ncid, grid_corner_lat_id, "units", "degrees")
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining attributes for grid_corner_lat'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      status = nf90_def_var(ncid, "grid_corner_lon", nf90_double,           &
           (/grid_corners_id, grid_size_id/), grid_corner_lon_id)
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining variable grid_corner_lon'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      status = nf90_put_att(ncid, grid_corner_lon_id, "units", "degrees")
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining attributes for grid_corner_lon'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      status = nf90_def_var(ncid, "grid_imask", nf90_double,                &
           (/grid_size_id/), grid_imask_id)
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error defining variable grid_imask'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if

      ! End of NetCDF definitions
      status = nf90_enddef(ncid)
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error calling enddef'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if
    end if ! IOroot

    if (IOroot) then
      status = nf90_put_var(ncid, grid_dims_id, (/ gridsize /))
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error writing variable grid_dims'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if
    end if

    do ie=1,nelemd
      lwork(:,:,ie) = fvm(ie)%area_sphere_physgrid(:,:)
    end do
    call mpi_gatherv(lwork, size(lwork), mpi_real8, recvbuf, recvcount, &
         displs, mpi_real8, IOrootID, mpicom, ierror)
    if (IOroot) then
      call dp_reoorder(recvbuf, gwork(1,:))
      status = nf90_put_var(ncid, grid_area_id, gwork(1,:))
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error writing variable grid_area'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if
    end if
    do ie=1,nelemd
      lwork(:,:,ie) = rad2deg*fvm(ie)%center_cart_physgrid(:,:)%lat
    end do
    call mpi_gatherv(lwork, size(lwork), mpi_real8, recvbuf, recvcount,       &
         displs, mpi_real8, IOrootID, mpicom, ierror)
    if (IOroot) then
      call dp_reoorder(recvbuf, gwork(1,:))
      status = nf90_put_var(ncid, grid_center_lat_id, gwork(1,:))
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error writing variable grid_center_lat'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if
    end if

    do ie=1,nelemd
      lwork(:,:,ie) = rad2deg*fvm(ie)%center_cart_physgrid(:,:)%lon
    end do
    call mpi_gatherv(lwork, size(lwork), mpi_real8, recvbuf, recvcount,       &
         displs, mpi_real8, IOrootID, mpicom, ierror)
    if (IOroot) then
      call dp_reoorder(recvbuf, gwork(1,:))
      status = nf90_put_var(ncid, grid_center_lon_id, gwork(1,:))
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error writing variable grid_center_lon'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if
    end if
    ! compute physgrid grid corners
    do ie=1,nelemd
      do j=1,fv_nphys
        do i=1,fv_nphys
          do ivtx=1,4
            x      = fvm(ie)%vtx_cart_physgrid(ivtx,1,i,j)
            y      = fvm(ie)%vtx_cart_physgrid(ivtx,2,i,j)
            sphere = cart2spherical(x,y,elem(ie)%FaceNum)
            corners(i,j,ie,ivtx,1) = rad2deg * sphere%lat
            corners(i,j,ie,ivtx,2) = rad2deg * sphere%lon
          end do
        end do
      end do
    end do
    ! Collect all information for the grid corner latitude (counter-clockwise)
    do ivtx=1,4
      call mpi_gatherv(corners(:,:,:,ivtx,1), size(corners(:,:,:,ivtx,1)), mpi_real8, recvbuf, recvcount,     &
           displs, mpi_real8, IOrootID, mpicom, ierror)
      if (IOroot) then
        call dp_reoorder(recvbuf, gwork(ivtx,:))
      end if
    end do
    if (IOroot) then
      status = nf90_put_var(ncid, grid_corner_lat_id, gwork)
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error writing variable grid_corner_lat'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if
    end if
    ! Collect all information for the grid corner longitudes (counter-clockwise)
    do ivtx=1,4
      call mpi_gatherv(corners(:,:,:,ivtx,2), size(corners(:,:,:,ivtx,2)), mpi_real8, recvbuf, recvcount,     &
           displs, mpi_real8, IOrootID, mpicom, ierror)
      if (IOroot) then
        call dp_reoorder(recvbuf, gwork(ivtx,:))
      end if
    end do
    if (IOroot) then
      status = nf90_put_var(ncid, grid_corner_lon_id, gwork)
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error writing variable grid_corner_lon'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if
    end if

    if (IOroot) then
      gwork(1,:) = 1._r8
      status = nf90_put_var(ncid, grid_imask_id, gwork(1,:))
      if (status /= nf90_noerr) then
        write(iulog, *) 'dp_write: Error writing variable grid_imask'
        call shr_sys_flush(iulog)
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if
    end if

!    call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
    ! Close the file
    call mpi_barrier(mpicom, ierror)
    if (IOroot) then
      status = nf90_close(ncid)
      if (status /= nf90_noerr) then
        call endrun("dp_write: "//trim(nf90_strerror(status)))
      end if
    end if

    call mpi_barrier(mpicom, ierror)
    if(masterproc) then
      write(iulog, *) 'Finished writing physics grid file: ', trim(filename)
      call shr_sys_flush(iulog)
    end if

  end subroutine dp_write
  !!!

  subroutine fvm2phys_init(elem,fvm,fvm_nc,phys_nc,irecons,&
            weights_all_fvm2phys,weights_eul_index_all_fvm2phys,weights_lgr_index_all_fvm2phys,&
            weights_all_phys2fvm,weights_eul_index_all_phys2fvm,weights_lgr_index_all_phys2fvm,&
            jall_fvm2phys,jall_phys2fvm)
    use dimensions_mod  , only: ngpc,nelemd
    use fvm_overlap_mod , only: compute_weights_cell
    use element_mod     , only: element_t
    type(element_t)  , dimension(nelemd), intent(in) :: elem
    type (fvm_struct), dimension(nelemd), intent(in) :: fvm
    integer        , intent(in) :: fvm_nc, phys_nc, irecons
    real (kind=r8)       :: dalpha,dbeta
    real (kind=r8), dimension(0:phys_nc+2):: xgno_phys,ygno_phys
    real (kind=r8), dimension(0:fvm_nc+2) :: xgno_fvm,ygno_fvm

    real (kind=r8), dimension(ngpc):: gauss_weights, abscissae !dimension(ngauss)

    integer :: i,j,h
    integer, parameter :: nvertex = 4
    real (kind=r8), dimension(nvertex)            :: xcell,ycell

    real (kind=r8)   , dimension(num_weights_fvm2phys,irecons,nelemd),intent(out) :: weights_all_fvm2phys
    integer,  dimension(num_weights_fvm2phys,2,nelemd),intent(out) :: weights_eul_index_all_fvm2phys
    integer,  dimension(num_weights_fvm2phys,2,nelemd),intent(out) :: weights_lgr_index_all_fvm2phys

    real (kind=r8)   , dimension(num_weights_phys2fvm,irecons,nelemd),intent(out) :: weights_all_phys2fvm
    integer,  dimension(num_weights_phys2fvm,2,nelemd),intent(out) :: weights_eul_index_all_phys2fvm
    integer,  dimension(num_weights_phys2fvm,2,nelemd),intent(out) :: weights_lgr_index_all_phys2fvm

    integer                , dimension(nelemd)                       ,intent(out) :: jall_fvm2phys,jall_phys2fvm

    integer, parameter :: jmax_segments_cell = 50
    real (kind=r8)   , dimension(jmax_segments_cell,irecons)  :: weights_cell
    integer , dimension(jmax_segments_cell,2)        :: weights_eul_index_cell
    integer :: jcollect_cell,ie
    real(kind=r8), dimension(phys_nc,phys_nc) :: phys_area, factor
    real(kind=r8), dimension(fvm_nc,fvm_nc)   :: fvm_area, facfvm

    xgno_phys(0) = -1D20; xgno_phys(phys_nc+2) = 1D20
    xgno_fvm(0)  = -1D20; xgno_fvm(fvm_nc+2)   = 1D20
    do ie=1,nelemd 
      dalpha         = abs(elem(ie)%corners(1)%x-elem(ie)%corners(2)%x)/phys_nc !in alpha
      dbeta          = abs(elem(ie)%corners(1)%y-elem(ie)%corners(4)%y)/phys_nc  !in beta
      do i=1,phys_nc+1
        xgno_phys(i) = tan(elem(ie)%corners(1)%x+(i-1)*dalpha)
        ygno_phys(i) = tan(elem(ie)%corners(1)%y+(i-1)*dbeta )
      end do
      
      dalpha         = abs(elem(ie)%corners(1)%x-elem(ie)%corners(2)%x)/fvm_nc !in alpha
      dbeta          = abs(elem(ie)%corners(1)%y-elem(ie)%corners(4)%y)/fvm_nc  !in beta
      do i=1,fvm_nc+1
        xgno_fvm(i) = tan(elem(ie)%corners(1)%x+(i-1)*dalpha)
        ygno_fvm(i) = tan(elem(ie)%corners(1)%y+(i-1)*dbeta )
      end do
      
      !
      ! compute area using line-integrals
      !
      !       do j=1,phys_nc
      !          do i=1,phys_nc
      !             da_phys(i,j) = (I_00(xgno_phys(i+1),ygno_phys(j+1)) - I_00(xgno_phys(i  ),ygno_phys(j+1)) + &
      !                  I_00(xgno_phys(i  ),ygno_phys(j  )) - I_00(xgno_phys(i+1),ygno_phys(j  )))
      !          end do
      !       end do
      !
      !       do j=1,fvm_nc
      !          do i=1,fvm_nc
      !             da_fvm(i,j) = (I_00(xgno_fvm(i+1),ygno_fvm(j+1)) - I_00(xgno_fvm(i  ),ygno_fvm(j+1)) + &
      !                  I_00(xgno_fvm(i  ),ygno_fvm(j  )) - I_00(xgno_fvm(i+1),ygno_fvm(j  )))
      !          end do
      !       end do
      
      gauss_weights = 0.0D0; abscissae=0.0D0!not used since line-segments are parallel to coordinate
      
      jall_fvm2phys(ie)=1
      do j=1,phys_nc
        do i=1,phys_nc
          xcell(1) = xgno_phys(i)  ; ycell(1) = ygno_phys(j)
          xcell(2) = xgno_phys(i)  ; ycell(2) = ygno_phys(j+1)
          xcell(3) = xgno_phys(i+1); ycell(3) = ygno_phys(j+1)
          xcell(4) = xgno_phys(i+1); ycell(4) = ygno_phys(j)
          
          call compute_weights_cell(nvertex,.true.,&
               xcell,ycell,i,j,irecons,xgno_fvm,ygno_fvm,0,fvm_nc+2,&
               1,fvm_nc+1,1,fvm_nc+1,&
               ngpc,gauss_weights,abscissae,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          
          if (jcollect_cell>0) then
            weights_all_fvm2phys(jall_fvm2phys(ie):jall_fvm2phys(ie)+jcollect_cell-1,:,ie) = &
                 weights_cell(1:jcollect_cell,:)!/fvm(ie)%area_sphere_physgrid(i,j)!da_phys(i,j) 
            
            weights_eul_index_all_fvm2phys(jall_fvm2phys(ie):jall_fvm2phys(ie)+jcollect_cell-1,:,ie) = &
                 weights_eul_index_cell(1:jcollect_cell,:)
            weights_lgr_index_all_fvm2phys(jall_fvm2phys(ie):jall_fvm2phys(ie)+jcollect_cell-1,1,ie) = i
            weights_lgr_index_all_fvm2phys(jall_fvm2phys(ie):jall_fvm2phys(ie)+jcollect_cell-1,2,ie) = j
            jall_fvm2phys(ie) = jall_fvm2phys(ie)+jcollect_cell
          endif
        end do
      enddo
      jall_fvm2phys(ie)=jall_fvm2phys(ie)-1
      !
      ! make sure sum of area overlap weights exactly match fvm%%area_sphere_physgrid
      !
      phys_area = 0.0_r8
      do h=1,jall_fvm2phys(ie)
        i  = weights_lgr_index_all_fvm2phys(h,1,ie); j  = weights_lgr_index_all_fvm2phys(h,2,ie)
        phys_area(i,j) = phys_area(i,j) +weights_all_fvm2phys(h,1,ie)
      end do
      factor(:,:) = fvm(ie)%area_sphere_physgrid(:,:)/phys_area(:,:)
      do h=1,jall_fvm2phys(ie)
        i  = weights_lgr_index_all_fvm2phys(h,1,ie); j  = weights_lgr_index_all_fvm2phys(h,2,ie)
        weights_all_fvm2phys(h,1,ie) = weights_all_fvm2phys(h,1,ie)*factor(i,j)
      end do
      
      jall_phys2fvm(ie)=1
      do j=1,fvm_nc
        do i=1,fvm_nc
          xcell(1) = xgno_fvm(i)  ; ycell(1) = ygno_fvm(j)
          xcell(2) = xgno_fvm(i)  ; ycell(2) = ygno_fvm(j+1)
          xcell(3) = xgno_fvm(i+1); ycell(3) = ygno_fvm(j+1)
          xcell(4) = xgno_fvm(i+1); ycell(4) = ygno_fvm(j)
          
          call compute_weights_cell(nvertex,.true.,&
               xcell,ycell,i,j,irecons,xgno_phys,ygno_phys,0,phys_nc+2,&
               1,phys_nc+1,1,phys_nc+1,&
               ngpc,gauss_weights,abscissae,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          
          if (jcollect_cell>0) then
            weights_all_phys2fvm(jall_phys2fvm(ie):jall_phys2fvm(ie)+jcollect_cell-1,:,ie) &
                 = weights_cell(1:jcollect_cell,:)!/fvm(ie)%area_sphere(i,j)!da_fvm(i,j)
            
            weights_eul_index_all_phys2fvm(jall_phys2fvm(ie):jall_phys2fvm(ie)+jcollect_cell-1,:,ie) = &
                 weights_eul_index_cell(1:jcollect_cell,:)
            weights_lgr_index_all_phys2fvm(jall_phys2fvm(ie):jall_phys2fvm(ie)+jcollect_cell-1,1,ie) = i
            weights_lgr_index_all_phys2fvm(jall_phys2fvm(ie):jall_phys2fvm(ie)+jcollect_cell-1,2,ie) = j
            jall_phys2fvm(ie) = jall_phys2fvm(ie)+jcollect_cell
          endif
        end do
      enddo
      jall_phys2fvm(ie)=jall_phys2fvm(ie)-1
      !
      ! make sure sum of area overlap weights exactly matches fvm%%area_sphere_physgrid
      !
      fvm_area = 0.0_r8
      do h=1,jall_phys2fvm(ie)
        i  = weights_lgr_index_all_phys2fvm(h,1,ie); j  = weights_lgr_index_all_phys2fvm(h,2,ie)
        fvm_area(i,j) = fvm_area(i,j) +weights_all_phys2fvm(h,1,ie)
      end do
      fvm_area(:,:) = fvm(ie)%area_sphere(:,:)/fvm_area(:,:)
      do h=1,jall_phys2fvm(ie)
        i  = weights_lgr_index_all_phys2fvm(h,1,ie); j  = weights_lgr_index_all_phys2fvm(h,2,ie)
        weights_all_phys2fvm(h,1,ie) = weights_all_phys2fvm(h,1,ie)*fvm_area(i,j)
      end do
    end do
    end subroutine fvm2phys_init
end module dp_mapping
