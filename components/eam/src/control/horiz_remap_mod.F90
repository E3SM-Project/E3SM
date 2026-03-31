module horiz_remap_mod
  !-------------------------------------------------------------------------------------------
  !
  ! EAM wrapper for shared horizontal remapping (shr_horiz_remap_mod).
  !
  ! Delegates map reading, CRS SpMV, and source-gather MPI to the shared
  ! module.  Only EAM-specific concerns remain here: chunk/icol indexing,
  ! cam_grid registration, PIO output decomposition, and packing source
  ! data from the chunk layout into the send buffer.
  !
  !-------------------------------------------------------------------------------------------

  use shr_kind_mod,        only: r8 => shr_kind_r8
  use shr_horiz_remap_mod, only: shr_horiz_remap_t, shr_horiz_remap_read_mapfile, &
                                 shr_horiz_remap_build_comm, shr_horiz_remap_apply
  use cam_logfile,         only: iulog
  use cam_abortutils,      only: endrun
  use spmd_utils,          only: masterproc, iam, npes, mpicom
  use cam_history_support, only: ptapes
  use ppgrid,              only: pcols, begchunk, endchunk

  implicit none
  private
  save

  public :: horiz_remap_init
  public :: horiz_remap_field
  public :: horiz_remap_write
  public :: horiz_remap_is_active
  public :: horiz_remap_get_grid_id

  type eam_horiz_remap_t
    type(shr_horiz_remap_t) :: shared
    integer :: grid_id = 0

    ! Send-side column mapping: for each entry in send_gcol_list,
    ! the (chunk, col-in-chunk) pair for packing the send buffer.
    integer, allocatable :: send_cols_chunk(:)
    integer, allocatable :: send_cols_icol(:)

    ! Persistent send buffer workspace
    real(r8), allocatable :: ws_send_buf(:)

    ! PIO decomposition cache (keyed on tape + numlev + data_type)
    logical :: iodesc_2d_valid = .false.
    integer :: iodesc_2d_dtype = 0
    logical :: iodesc_3d_valid = .false.
    integer :: iodesc_3d_nlev = 0
    integer :: iodesc_3d_dtype = 0
  end type eam_horiz_remap_t

  type(eam_horiz_remap_t), target :: remap_data(ptapes)

CONTAINS

  !-------------------------------------------------------------------------------------------
  logical function horiz_remap_is_active(t)
    integer, intent(in) :: t
    horiz_remap_is_active = remap_data(t)%shared%initialized
  end function horiz_remap_is_active

  !-------------------------------------------------------------------------------------------
  integer function horiz_remap_get_grid_id(t)
    integer, intent(in) :: t
    horiz_remap_get_grid_id = remap_data(t)%grid_id
  end function horiz_remap_get_grid_id

  !-------------------------------------------------------------------------------------------
  subroutine horiz_remap_init(t, mapfile)
    use pio,              only: iosystem_desc_t
    use shr_pio_mod,      only: shr_pio_getiosys
    use cam_instance,     only: atm_id
    use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap, &
                                cam_grid_register, cam_grid_attribute_register
    use phys_grid,        only: get_ncols_p, get_gcol_p, chunks, knuhcs

    integer,          intent(in) :: t
    character(len=*), intent(in) :: mapfile

    type(shr_horiz_remap_t), pointer :: rd
    type(eam_horiz_remap_t), pointer :: erd
    type(iosystem_desc_t), pointer :: pio_subsystem
    integer :: ierr, i, cnt, idx, lchnk, ncols, icol
    integer :: n_my_cols
    integer, allocatable :: gcol_to_rank(:), send_gcol_list(:)
    integer, allocatable :: gcol_to_myidx(:), my_gcol_chunk(:), my_gcol_icol(:)
    type(horiz_coord_t), pointer :: lat_coord, lon_coord
    integer(iMap), pointer :: grid_map(:,:)
    character(len=32) :: gridname
    character(len=*), parameter :: subname = 'horiz_remap_init'

    erd => remap_data(t)
    rd  => erd%shared

    if (masterproc) then
      write(iulog,*) trim(subname), ': Reading mapping file for tape ', t, ': ', trim(mapfile)
    end if

    ! Phase 1: read map file
    pio_subsystem => shr_pio_getiosys(atm_id)
    call shr_horiz_remap_read_mapfile(rd, mapfile, mpicom, iam, npes, pio_subsystem, ierr)
    if (ierr /= 0) then
      call endrun(trim(subname)//': Error reading mapping file')
    end if

    ! Validate n_a
    if (rd%n_a /= size(knuhcs)) then
      if (masterproc) then
        write(iulog,*) trim(subname), ': ERROR: map n_a=', rd%n_a, &
             ' but EAM has ', size(knuhcs), ' columns'
      end if
      call endrun(trim(subname)//': n_a mismatch between map file and EAM grid')
    end if

    if (masterproc) then
      write(iulog,*) trim(subname), ': n_a=', rd%n_a, ' n_b=', rd%n_b
      write(iulog,*) trim(subname), ': target grid nlat=', rd%nlat, ' nlon=', rd%nlon
    end if

    ! Fill gcol_to_rank(1:n_a) using EAM decomposition
    allocate(gcol_to_rank(rd%n_a))
    do i = 1, rd%n_a
      gcol_to_rank(i) = chunks(knuhcs(i)%chunkid)%owner
    end do

    ! Phase 2: build comm pattern + CRS matrix
    call shr_horiz_remap_build_comm(rd, gcol_to_rank, mpicom, iam, npes, &
         send_gcol_list, ierr)
    deallocate(gcol_to_rank)

    if (ierr /= 0) then
      call endrun(trim(subname)//': Error building comm pattern')
    end if

    ! Convert send_gcol_list to (chunk, icol) pairs
    n_my_cols = 0
    do lchnk = begchunk, endchunk
      n_my_cols = n_my_cols + get_ncols_p(lchnk)
    end do

    allocate(gcol_to_myidx(rd%n_a))
    allocate(my_gcol_chunk(n_my_cols), my_gcol_icol(n_my_cols))
    gcol_to_myidx(:) = 0
    cnt = 0
    do lchnk = begchunk, endchunk
      ncols = get_ncols_p(lchnk)
      do icol = 1, ncols
        cnt = cnt + 1
        gcol_to_myidx(get_gcol_p(lchnk, icol)) = cnt
        my_gcol_chunk(cnt) = lchnk
        my_gcol_icol(cnt) = icol
      end do
    end do

    allocate(erd%send_cols_chunk(rd%n_send_total))
    allocate(erd%send_cols_icol(rd%n_send_total))
    do i = 1, rd%n_send_total
      idx = gcol_to_myidx(send_gcol_list(i))
      if (idx == 0) then
        call endrun(trim(subname)//': Requested gcol not owned by this rank')
      end if
      erd%send_cols_chunk(i) = my_gcol_chunk(idx)
      erd%send_cols_icol(i) = my_gcol_icol(idx)
    end do

    deallocate(send_gcol_list, gcol_to_myidx, my_gcol_chunk, my_gcol_icol)

    ! Register the lat-lon output grid
    nullify(grid_map)
    erd%grid_id = 300 + t
    write(gridname, '(a,i0)') 'horiz_remap_', t

    lat_coord => horiz_coord_create('lat', '', rd%nlat, 'latitude', 'degrees_north', &
         1, rd%nlat, rd%lat)
    lon_coord => horiz_coord_create('lon', '', rd%nlon, 'longitude', 'degrees_east', &
         1, rd%nlon, rd%lon)
    call cam_grid_register(trim(gridname), erd%grid_id, &
         lat_coord, lon_coord, grid_map, unstruct=.false.)
    call cam_grid_attribute_register(trim(gridname), &
         'horiz_remap_file', trim(mapfile))

    if (masterproc) then
      write(iulog,*) trim(subname), ': Initialized for tape ', t, &
           ', n_b_local=', rd%n_b_local
    end if

  end subroutine horiz_remap_init

  !-------------------------------------------------------------------------------------------
  subroutine horiz_remap_field(t, hbuf, numlev, fld_out)
    !
    ! Remap a field from the physics grid to the target lat-lon grid.
    !
    ! Input: hbuf(pcols, numlev, begchunk:endchunk)
    ! Output: fld_out(n_b_local, numlev) - allocated here
    !
    integer,  intent(in)    :: t
    real(r8), intent(in)    :: hbuf(:,:,:)
    integer,  intent(in)    :: numlev
    real(r8), allocatable, intent(out) :: fld_out(:,:)

    type(eam_horiz_remap_t), pointer :: erd
    type(shr_horiz_remap_t), pointer :: rd
    integer :: i, k, lchnk, icol, ierr, needed
    character(len=*), parameter :: subname = 'horiz_remap_field'

    erd => remap_data(t)
    rd  => erd%shared

    if (.not. rd%initialized) then
      call endrun(trim(subname)//': Remapping not initialized for this tape')
    end if

    allocate(fld_out(rd%n_b_local, numlev))

    ! Grow persistent send workspace if needed
    needed = rd%n_send_total * numlev
    if (needed > 0) then
      if (.not. allocated(erd%ws_send_buf) .or. size(erd%ws_send_buf) < needed) then
        if (allocated(erd%ws_send_buf)) deallocate(erd%ws_send_buf)
        allocate(erd%ws_send_buf(needed))
      end if
    end if

    ! Pack send buffer from EAM chunk layout
    do i = 1, rd%n_send_total
      lchnk = erd%send_cols_chunk(i)
      icol  = erd%send_cols_icol(i)
      do k = 1, numlev
        erd%ws_send_buf((i-1)*numlev + k) = hbuf(icol, k, lchnk - begchunk + 1)
      end do
    end do

    ! Apply remap: alltoallv + CRS SpMV
    call shr_horiz_remap_apply(rd, erd%ws_send_buf, numlev, fld_out, mpicom, npes, ierr)

  end subroutine horiz_remap_field

  !-------------------------------------------------------------------------------------------
  subroutine horiz_remap_write(t, File, varid, fld_out, numlev, data_type)
    use pio,              only: file_desc_t, var_desc_t, io_desc_t, &
                                pio_initdecomp, pio_freedecomp, &
                                pio_write_darray, iosystem_desc_t, &
                                PIO_OFFSET_KIND
    use shr_pio_mod,      only: shr_pio_getiosys
    use cam_instance,     only: atm_id

    integer,           intent(in)    :: t
    type(file_desc_t), intent(inout) :: File
    type(var_desc_t),  intent(inout) :: varid
    real(r8),          intent(in)    :: fld_out(:,:)
    integer,           intent(in)    :: numlev
    integer,           intent(in)    :: data_type

    type(eam_horiz_remap_t), pointer :: erd
    type(shr_horiz_remap_t), pointer :: rd
    type(io_desc_t), save :: iodesc_2d_cache(ptapes)
    type(io_desc_t), save :: iodesc_3d_cache(ptapes)
    type(iosystem_desc_t), pointer :: pio_subsystem
    integer(PIO_OFFSET_KIND), allocatable :: idof(:)
    integer :: i, k, global_row, ilon, ilat, ierr
    integer :: nlat, nlon

    erd => remap_data(t)
    rd  => erd%shared
    nlat = rd%nlat
    nlon = rd%nlon

    pio_subsystem => shr_pio_getiosys(atm_id)

    if (numlev <= 1) then
      if (erd%iodesc_2d_valid .and. erd%iodesc_2d_dtype /= data_type) then
        call pio_freedecomp(File, iodesc_2d_cache(t))
        erd%iodesc_2d_valid = .false.
      end if
      if (.not. erd%iodesc_2d_valid) then
        allocate(idof(rd%n_b_local))
        do i = 1, rd%n_b_local
          global_row = rd%row_start + i - 1
          ilon = mod(global_row - 1, nlon) + 1
          ilat = (global_row - 1) / nlon + 1
          idof(i) = int(ilon + nlon * (ilat - 1), PIO_OFFSET_KIND)
        end do
        call pio_initdecomp(pio_subsystem, data_type, (/nlon, nlat/), idof, iodesc_2d_cache(t))
        deallocate(idof)
        erd%iodesc_2d_valid = .true.
        erd%iodesc_2d_dtype = data_type
      end if
      call pio_write_darray(File, varid, iodesc_2d_cache(t), fld_out(:,1), ierr)
    else
      if (erd%iodesc_3d_valid .and. &
          (erd%iodesc_3d_nlev /= numlev .or. erd%iodesc_3d_dtype /= data_type)) then
        call pio_freedecomp(File, iodesc_3d_cache(t))
        erd%iodesc_3d_valid = .false.
      end if
      if (.not. erd%iodesc_3d_valid) then
        allocate(idof(rd%n_b_local * numlev))
        do i = 1, rd%n_b_local
          global_row = rd%row_start + i - 1
          ilon = mod(global_row - 1, nlon) + 1
          ilat = (global_row - 1) / nlon + 1
          do k = 1, numlev
            idof((k-1)*rd%n_b_local + i) = int(ilon + nlon*(ilat-1) + nlon*nlat*(k-1), PIO_OFFSET_KIND)
          end do
        end do
        call pio_initdecomp(pio_subsystem, data_type, (/nlon, nlat, numlev/), idof, iodesc_3d_cache(t))
        deallocate(idof)
        erd%iodesc_3d_valid = .true.
        erd%iodesc_3d_nlev = numlev
        erd%iodesc_3d_dtype = data_type
      end if
      call pio_write_darray(File, varid, iodesc_3d_cache(t), fld_out, ierr)
    end if

  end subroutine horiz_remap_write

end module horiz_remap_mod
