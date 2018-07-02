module ncdio_atm

  !----------------------------------------------------------------------- 
  !BOP
  !
  ! !MODULE: ncdio_atm
  !
  ! !DESCRIPTION: 
  ! Generic interfaces to write fields to PIO files
  !
  ! !USES:

  use pio, only: pio_offset_kind, file_desc_t, var_desc_t, pio_double,        &
       pio_inq_dimlen, pio_inq_dimid, pio_max_var_dims,                       &
       pio_inq_vardimid, pio_inq_varndims, io_desc_t, pio_setframe
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_sys_mod,    only: shr_sys_flush      ! Standardized system subroutines
  use shr_scam_mod,   only: shr_scam_getCloseLatLon  ! Standardized system subroutines
  use spmd_utils,     only: masterproc
  use cam_abortutils, only: endrun
  use scamMod,        only: scmlat,scmlon,single_column
  use cam_logfile,    only: iulog
  !
  ! !PUBLIC TYPES:
  implicit none

  PRIVATE

  save

  logical :: debug = .false.

  !
  !EOP
  !
  interface infld
     module procedure infld_real_1d_2d
     module procedure infld_real_2d_2d
     module procedure infld_real_2d_3d
     module procedure infld_real_3d_3d
  end interface


  public :: infld

  integer STATUS
  real(r8) surfdat
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: infld_real_1d_2d
  !
  ! !INTERFACE:
  subroutine infld_real_1d_2d(varname, ncid, dimname1,                        &
       dim1b, dim1e, dim2b, dim2e, field, readvar, gridname, timelevel)
    !
    ! !DESCRIPTION: 
    ! Netcdf I/O of initial real field from netCDF file
    ! Read a 1-D field (or slice) into a 2-D variable
    !
    ! !USES
    !

    use pio,              only: pio_get_var, pio_read_darray, pio_setdebuglevel
    use pio,              only: PIO_MAX_NAME, pio_inquire, pio_inq_dimname
    use cam_grid_support, only: cam_grid_check, cam_grid_get_decomp, cam_grid_id
    use cam_pio_utils,    only: cam_pio_check_var

    !
    ! !ARGUMENTS:
    implicit none
    character(len=*),  intent(in)     :: varname  ! variable name
    type(file_desc_t), intent(inout)  :: ncid     ! input unit
    character(len=*),  intent(in)     :: dimname1 ! name of 1st array dimensions of field on file (array order)
    integer,           intent(in)     :: dim1b    ! start of first  dimension of array to be returned
    integer,           intent(in)     :: dim1e    ! end   of first  dimension of array to be returned
    integer,           intent(in)     :: dim2b    ! start of second dimension of array to be returned
    integer,           intent(in)     :: dim2e    ! end   of second dimension of array to be returned
    real(r8), target,  intent(out)    :: field(dim1b:dim1e,dim2b:dim2e) ! array to be returned (decomposed or global)
    logical,           intent(out)    :: readvar  ! true => variable is on initial dataset
    character(len=*), optional, intent(in) :: gridname ! Name of variable's grid
    integer, optional, intent(in)     :: timelevel
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    type(io_desc_t), pointer  :: iodesc
    integer                   :: grid_map  ! grid ID for data mapping
    integer                   :: i, j      ! indices
    integer                   :: ierr      ! error status
    type(var_desc_t)          :: varid     ! variable id

    integer                   :: arraydimsize(2) ! field dimension lengths
    integer                   :: arraydimid      ! Dimension ID

    integer                   :: ndims ! number of dimensions
    integer                   :: dimids(PIO_MAX_VAR_DIMS) ! file variable dims
    integer                   :: dimlens(PIO_MAX_VAR_DIMS) ! file variable shape

    ! Offsets for reading global variables
    integer                   :: strt(1) = 1 ! start ncol index for netcdf 1-d
    integer                   :: cnt (1) = 1 ! ncol count for netcdf 1-d
    character(len=PIO_MAX_NAME) :: tmpname
    character(len=128)        :: errormsg

    logical                   :: readvar_tmp ! if true, variable is on tape
    character(len=32)         :: subname='INFLD_REAL_1D_2D' ! subroutine name

    ! For SCAM
    real(r8)                  :: closelat, closelon
    integer                   :: lonidx, latidx

    nullify(iodesc)

    !
    !-----------------------------------------------------------------------
    !
    !    call pio_setdebuglevel(3)

    !
    ! Error conditions
    !
    if (present(gridname)) then
      grid_map = cam_grid_id(trim(gridname))
    else
      grid_map = cam_grid_id('physgrid')
    end if
    if (.not. cam_grid_check(grid_map)) then
      if(masterproc) then
        if (present(gridname)) then
          write(errormsg, *)': invalid gridname, "',trim(gridname),'", specified for field ',trim(varname)
        else
          write(errormsg, *)': Internal error, no "physgrid" gridname'
        end if
      end if
      call endrun(trim(subname)//errormsg)
    end if

    if (debug .and. masterproc) then
      if (present(gridname)) then
        write(iulog, '(5a)') trim(subname),': field = ',trim(varname),', grid = ',trim(gridname)
      else
        write(iulog, '(4a)') trim(subname),': field = ',trim(varname),', grid = physgrid'
      end if
      call shr_sys_flush(iulog)
    end if
    !
    ! Read netCDF file
    !
    !
    ! Check if field is on file; get netCDF variable id
    !
    call cam_pio_check_var(ncid, varname, varid, ndims, dimids, dimlens, readvar_tmp)
    !
    ! If field is on file:
    !
    if (readvar_tmp) then
        if (debug .and. masterproc) then
          write(iulog, '(2a,5(i0,a))') trim(subname),': field(',              &
               dim1b,':',dim1e,',',dim2b,':',dim2e, '), file(',dimlens(1),')'
          call shr_sys_flush(iulog)
        end if
      !
      ! Get array dimension id's and sizes
      !
      ierr = PIO_inq_dimid(ncid, dimname1, arraydimid)
      arraydimsize(1) = (dim1e - dim1b + 1)
      arraydimsize(2) = (dim2e - dim2b + 1)
      do j = 1, 2
        if (arraydimsize(j) /= size(field, j)) then
          write(errormsg, *) ': Mismatch between array bounds and field size for ', &
               trim(varname), ', dimension', j
          call endrun(trim(subname)//errormsg)
        end if
      end do

      if (ndims > 2) then
        call endrun(trim(subname)//': too many dimensions for '//trim(varname))
      else if (ndims < 1) then
        call endrun(trim(subname)//': too few dimensions for '//trim(varname))
      else
        ! Check to make sure that the second dimension is time
        if (ndims == 2) then
          ierr = pio_inq_dimname(ncid, dimids(2), tmpname)
          if (trim(tmpname) /= 'time') then
            call endrun(trim(subname)//': dimension mismatch for '//trim(varname))
          end if
        end if
      end if

      if(ndims == 2) then
        if(present(timelevel)) then
          call pio_setframe(ncid, varid, int(timelevel,kind=pio_offset_kind))
        end if
      end if

      if (single_column .and. dim1e == 1) then
        ! Specifically, this condition is for when the single column model 
        !  is run in the Spectral Element dycore
        cnt(1) = 1 
        call shr_scam_getCloseLatLon(ncid%fh,scmlat,scmlon,closelat,closelon,latidx,lonidx)
        strt(1) = lonidx
        ierr = pio_get_var(ncid, varid, strt, cnt, field)

      else

      ! NB: strt and cnt were initialized to 1
        ! All distributed array processing
        call cam_grid_get_decomp(grid_map, arraydimsize, dimlens(1:ndims),    &
             pio_double, iodesc)
        call pio_read_darray(ncid, varid, iodesc, field, ierr)
      end if
    end if  ! end of readvar_tmp

    readvar = readvar_tmp

    return

  end subroutine infld_real_1d_2d

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: infld_real_2d_2d
  !
  ! !INTERFACE:
  subroutine infld_real_2d_2d(varname, ncid, dimname1, dimname2,              &
       dim1b, dim1e, dim2b, dim2e, field, readvar, gridname, timelevel)
    !
    ! !DESCRIPTION: 
    ! Netcdf I/O of initial real field from netCDF file
    ! Read a 2-D field (or slice) into a 2-D variable
    !
    ! !USES
    !

    use pio,              only: pio_get_var, pio_read_darray, pio_setdebuglevel
    use pio,              only: PIO_MAX_NAME, pio_inquire, pio_inq_dimname
    use cam_grid_support, only: cam_grid_check, cam_grid_get_decomp, cam_grid_id
    use cam_pio_utils,    only: cam_permute_array, calc_permutation, cam_pio_check_var

    !
    ! !ARGUMENTS:
    implicit none
    character(len=*),  intent(in)     :: varname  ! variable name
    type(file_desc_t), intent(inout)  :: ncid     ! input unit
    character(len=*),  intent(in)     :: dimname1 ! name of 1st array dimensions of field on file (array order)
    character(len=*),  intent(in)     :: dimname2 ! name of 2nd array dimensions of field on file (array order)
    integer,           intent(in)     :: dim1b    ! start of first  dimension of array to be returned
    integer,           intent(in)     :: dim1e    ! end   of first  dimension of array to be returned
    integer,           intent(in)     :: dim2b    ! start of second dimension of array to be returned
    integer,           intent(in)     :: dim2e    ! end   of second dimension of array to be returned
    real(r8), target,  intent(out)    :: field(dim1b:dim1e,dim2b:dim2e) ! array to be returned (decomposed or global)
    logical,           intent(out)    :: readvar  ! true => variable is on initial dataset
    character(len=*), optional, intent(in) :: gridname ! Name of variable's grid
    integer, optional, intent(in)     :: timelevel
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    type(io_desc_t), pointer  :: iodesc
    integer                   :: grid_map  ! grid ID for data mapping
    integer                   :: i, j      ! indices
    integer                   :: ierr      ! error status
    type(var_desc_t)          :: varid     ! variable id

    integer                   :: arraydimsize(2) ! field dimension lengths
    integer                   :: arraydimids(2) ! Dimension IDs
    integer                   :: permutation(2)
    logical                   :: ispermuted

    integer                   :: ndims ! number of dimensions on file
    integer                   :: dimids(PIO_MAX_VAR_DIMS) ! file variable dims
    integer                   :: dimlens(PIO_MAX_VAR_DIMS) ! file variable shape

    ! Offsets for reading global variables
    integer                   :: strt(2) ! start lon, lat indices for netcdf 2-d
    integer                   :: cnt (2) ! lon, lat counts for netcdf 2-d
    character(len=PIO_MAX_NAME) :: tmpname
    character(len=128)        :: errormsg

    real(r8), pointer         :: tmp2d(:,:) ! input data for permutation

    logical                   :: readvar_tmp ! if true, variable is on tape
    character(len=32)         :: subname='INFLD_REAL_2D_2D' ! subroutine name
    character(len=PIO_MAX_NAME) :: field_dnames(2)

    ! For SCAM
    real(r8)                  :: closelat, closelon
    integer                   :: lonidx, latidx

    nullify(iodesc)

    !
    !-----------------------------------------------------------------------
    !
    !    call pio_setdebuglevel(3)

    ! Should we be using a different interface?
    if ((trim(dimname1) == trim(dimname2)) .or. (len_trim(dimname2) == 0)) then
      call infld(varname, ncid, dimname1, dim1b, dim1e, dim2b, dim2e,         &
           field, readvar, gridname, timelevel)
    else

      !
      ! Error conditions
      !
      if (present(gridname)) then
        grid_map = cam_grid_id(trim(gridname))
      else
        grid_map = cam_grid_id('physgrid')
      end if
      if (.not. cam_grid_check(grid_map)) then
        if(masterproc) then
          if (present(gridname)) then
            write(errormsg, *)': invalid gridname, "',trim(gridname),'", specified for field ',trim(varname)
          else
            write(errormsg, *)': Internal error, no "physgrid" gridname'
          end if
        end if
        call endrun(trim(subname)//errormsg)
      end if

    if (debug .and. masterproc) then
      if (present(gridname)) then
        write(iulog, '(5a)') trim(subname),': field = ',trim(varname),', grid = ',trim(gridname)
      else
        write(iulog, '(4a)') trim(subname),': field = ',trim(varname),', grid = physgrid'
      end if
      call shr_sys_flush(iulog)
    end if

      !
      ! Read netCDF file
      !
      !
      ! Check if field is on file; get netCDF variable id
      !
      call cam_pio_check_var(ncid, varname, varid, ndims, dimids, dimlens, readvar_tmp)
      !
      ! If field is on file:
      !
      if (readvar_tmp) then
        if (debug .and. masterproc) then
          write(iulog, '(2a,6(i0,a))') trim(subname),': field(',              &
               dim1b,':',dim1e,',',dim2b,':',dim2e,                           &
               '), file(',dimlens(1),',',dimlens(2),')'
          call shr_sys_flush(iulog)
        end if
        !
        ! Get array dimension id's and sizes
        !
        ierr = PIO_inq_dimid(ncid, dimname1, arraydimids(1))
        ierr = PIO_inq_dimid(ncid, dimname2, arraydimids(2))
        arraydimsize(1) = (dim1e - dim1b + 1)
        arraydimsize(2) = (dim2e - dim2b + 1)
        do j = 1, 2
          if (arraydimsize(j) /= size(field, j)) then
            write(errormsg, *) ': Mismatch between array bounds and field size for ', &
                 trim(varname), ', dimension', j
            call endrun(trim(subname)//errormsg)
          end if
        end do

        if (ndims > 3) then
          call endrun(trim(subname)//': too many dimensions for '//trim(varname))
        else if (ndims < 2) then
          call endrun(trim(subname)//': too few dimensions for '//trim(varname))
        else
          ! Check to make sure that the third dimension is time
          if (ndims == 3) then
            ierr = pio_inq_dimname(ncid, dimids(3), tmpname)
            if (trim(tmpname) /= 'time') then
              call endrun(trim(subname)//': dimension mismatch for '//trim(varname))
            end if
          end if
        end if

        if(ndims == 3) then
          if(present(timelevel)) then
            call pio_setframe(ncid, varid, int(timelevel,kind=pio_offset_kind))
          end if
        end if

        field_dnames(1) = dimname1
        field_dnames(2) = dimname2
        if (single_column) then	
          ! This could be generalized but for now only handles a single point
          strt(1) = dim1b
          strt(2) = dim2b
          cnt = arraydimsize
          call shr_scam_getCloseLatLon(ncid%fh,scmlat,scmlon,closelat,closelon,latidx,lonidx)
          if (trim(field_dnames(1)) == 'lon') then
            strt(1) = lonidx ! First dim always lon for Eulerian dycore
          else
            call endrun(trim(subname)//': lon should be first dimension for '//trim(varname))
          end if
          if (trim(field_dnames(2)) == 'lat') then
            strt(2) = latidx
          else
            call endrun(trim(subname)//': lat dimension not found for '//trim(varname))
          end if

          ! Check for permuted dimensions ('out of order' array)
          call calc_permutation(dimids, arraydimids, permutation, ispermuted)
          if (ispermuted) then
            call cam_permute_array(strt, permutation)
            call cam_permute_array(cnt, permutation)
            allocate(tmp2d(1:cnt(1), 1:cnt(2)))
            ierr = pio_get_var(ncid, varid, strt, cnt, tmp2d)
            do j = dim2b, dim2e
              do i = dim1b, dim1e
                ! We don't need strt anymore, reuse it
                strt(1) = i - dim1b + 1
                strt(2) = j - dim2b + 1
                call cam_permute_array(strt, permutation)
                field(i,j) = tmp2d(strt(1), strt(2))
              end do
            end do
            deallocate(tmp2d)
          else
            ierr = pio_get_var(ncid, varid, strt, cnt, field)
          end if
        else
          ! All distributed array processing
          call cam_grid_get_decomp(grid_map, arraydimsize, dimlens(1:2),      &
               pio_double, iodesc, field_dnames=field_dnames)
          call pio_read_darray(ncid, varid, iodesc, field, ierr)
        end if
      end if  ! end of readvar_tmp

      readvar = readvar_tmp

    end if ! end of call infld_real_1d_2d instead

  end subroutine infld_real_2d_2d


  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: infld_real_2d_3d
  !
  ! !INTERFACE:
  subroutine infld_real_2d_3d(varname, ncid, dimname1, dimname2,              &
       dim1b, dim1e, dim2b, dim2e, dim3b, dim3e,                              &
       field, readvar, gridname, timelevel)
    !
    ! !DESCRIPTION: 
    ! Netcdf I/O of initial real field from netCDF file
    ! Read a 2-D field (or slice) into a 3-D variable
    !
    ! !USES
    !

    use pio,              only: pio_get_var, pio_read_darray, pio_setdebuglevel
    use pio,              only: PIO_MAX_NAME, pio_inquire, pio_inq_dimname
    use cam_grid_support, only: cam_grid_check, cam_grid_get_decomp, cam_grid_id
    use cam_pio_utils,    only: cam_permute_array, calc_permutation, cam_pio_check_var

    !
    ! !ARGUMENTS:
    implicit none
    character(len=*),  intent(in)     :: varname  ! variable name
    type(file_desc_t), intent(inout)  :: ncid     ! input unit
    character(len=*),  intent(in)     :: dimname1 ! name of 1st array dimensions of field on file (array order)
    character(len=*),  intent(in)     :: dimname2 ! name of 2nd array dimensions of field on file (array order)
    integer,           intent(in)     :: dim1b    ! start of first  dimension of array to be returned
    integer,           intent(in)     :: dim1e    ! end   of first  dimension of array to be returned
    integer,           intent(in)     :: dim2b    ! start of second dimension of array to be returned
    integer,           intent(in)     :: dim2e    ! end   of second dimension of array to be returned
    integer,           intent(in)     :: dim3b    ! start of third  dimension of array to be returned
    integer,           intent(in)     :: dim3e    ! end   of third  dimension of array to be returned
    real(r8), target,  intent(out)    :: field(dim1b:dim1e,dim2b:dim2e,dim3b:dim3e) ! array to be returned (decomposed or global)
    logical,           intent(out)    :: readvar  ! true => variable is on initial dataset
    character(len=*), optional, intent(in) :: gridname ! Name of variable's grid
    integer, optional, intent(in)     :: timelevel
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    type(io_desc_t), pointer  :: iodesc
    integer                   :: grid_map  ! grid ID for data mapping
    integer                   :: i, j, k   ! indices
    integer                   :: ierr      ! error status
    type(var_desc_t)          :: varid     ! variable id

    integer                   :: arraydimsize(3) ! field dimension lengths
    integer                   :: arraydimids(2) ! Dimension IDs
    integer                   :: permutation(2)
    logical                   :: ispermuted

    integer                   :: ndims ! number of dimensions
    integer                   :: dimids(PIO_MAX_VAR_DIMS) ! file variable dims
    integer                   :: dimlens(PIO_MAX_VAR_DIMS) ! file variable shape

    ! Offsets for reading global variables
    integer                   :: strt(3) = 1 ! start ncol, lev indices for netcdf 2-d
    integer                   :: cnt (3) = 1 ! ncol, lev counts for netcdf 2-d
    character(len=PIO_MAX_NAME) :: tmpname

    real(r8), pointer         :: tmp3d(:,:,:) ! input data for permutation

    logical                   :: readvar_tmp ! if true, variable is on tape
    character(len=32)         :: subname='INFLD_REAL_2D_3D' ! subroutine name
    character(len=128)        :: errormsg
    character(len=PIO_MAX_NAME) :: field_dnames(2)

    ! For SCAM
    real(r8)                  :: closelat, closelon
    integer                   :: lonidx, latidx

    nullify(iodesc)

    !
    !-----------------------------------------------------------------------
    !
    !    call pio_setdebuglevel(3)

    !
    ! Error conditions
    !
    if (present(gridname)) then
      grid_map = cam_grid_id(trim(gridname))
    else
      grid_map = cam_grid_id('physgrid')
    end if
    if (.not. cam_grid_check(grid_map)) then
      if(masterproc) then
        if (present(gridname)) then
          write(errormsg, *)': invalid gridname, "',trim(gridname),'", specified for field ',trim(varname)
        else
          write(errormsg, *)': Internal error, no "physgrid" gridname'
        end if
      end if
      call endrun(trim(subname)//errormsg)
    end if

    if (debug .and. masterproc) then
      if (present(gridname)) then
        write(iulog, '(5a)') trim(subname),': field = ',trim(varname),', grid = ',trim(gridname)
      else
        write(iulog, '(4a)') trim(subname),': field = ',trim(varname),', grid = physgrid'
      end if
      call shr_sys_flush(iulog)
    end if

    !
    ! Read netCDF file
    !
    !
    ! Check if field is on file; get netCDF variable id
    !
    call cam_pio_check_var(ncid, varname, varid, ndims, dimids, dimlens, readvar_tmp)
    !
    ! If field is on file:
    !
    if (readvar_tmp) then
      if (debug .and. masterproc) then
        write(iulog, '(2a,8(i0,a))') trim(subname),': field(',                &
             dim1b,':',dim1e,',',dim2b,':',dim2e,',',dim3b,':',dim3e,         &
             '), file(',dimlens(1),',',dimlens(2),')'
        call shr_sys_flush(iulog)
      end if
      !
      ! Get array dimension id's and sizes
      !
      ierr = PIO_inq_dimid(ncid, dimname1, arraydimids(1))
      ierr = PIO_inq_dimid(ncid, dimname2, arraydimids(2))
      arraydimsize(1) = (dim1e - dim1b + 1)
      arraydimsize(2) = (dim2e - dim2b + 1)
      arraydimsize(3) = (dim3e - dim3b + 1)
      do j = 1, 3
        if (arraydimsize(j) /= size(field, j)) then
          write(errormsg, *) ': Mismatch between array bounds and field size for ', &
               trim(varname), ', dimension', j
          call endrun(trim(subname)//errormsg)
        end if
      end do

      if (ndims > 3) then
        call endrun(trim(subname)//': too many dimensions for '//trim(varname))
      else if (ndims < 2) then
        call endrun(trim(subname)//': too few dimensions for '//trim(varname))
      else
        ! Check to make sure that the fourth dimension is time
        if (ndims == 3) then
          ierr = pio_inq_dimname(ncid, dimids(3), tmpname)
          if (trim(tmpname) /= 'time') then
            call endrun(trim(subname)//': dimension mismatch for '//trim(varname))
          end if
        end if
      end if

      if(ndims == 3) then
        if(present(timelevel)) then
          call pio_setframe(ncid, varid, int(timelevel,kind=pio_offset_kind))
        end if
      end if

      field_dnames(1) = dimname1
      field_dnames(2) = dimname2
      ! NB: strt and cnt were initialized to 1
      call cam_grid_get_decomp(grid_map, arraydimsize, dimlens(1:2),        &
             pio_double, iodesc, field_dnames=field_dnames)
      call pio_read_darray(ncid, varid, iodesc, field, ierr)
    end if  ! end of readvar_tmp

    readvar = readvar_tmp

    return

  end subroutine infld_real_2d_3d

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: infld_real_3d_3d
  !
  ! !INTERFACE:
  subroutine infld_real_3d_3d(varname, ncid, dimname1, dimname2, dimname3,    &
       dim1b, dim1e, dim2b, dim2e, dim3b, dim3e,                              &
       field, readvar, gridname, timelevel)
    !
    ! !DESCRIPTION: 
    ! Netcdf I/O of initial real field from netCDF file
    ! Read a 3-D field (or slice) into a 3-D variable
    !
    ! !USES
    !

    use pio,              only: pio_get_var, pio_read_darray, pio_setdebuglevel
    use pio,              only: PIO_MAX_NAME, pio_inquire, pio_inq_dimname
    use cam_grid_support, only: cam_grid_check, cam_grid_get_decomp, cam_grid_id
    use cam_pio_utils,    only: cam_permute_array, calc_permutation, cam_pio_check_var

    !
    ! !ARGUMENTS:
    implicit none
    character(len=*),  intent(in)     :: varname  ! variable name
    type(file_desc_t), intent(inout)  :: ncid     ! input unit
    character(len=*),  intent(in)     :: dimname1 ! name of 1st array dimensions of field on file (array order)
    character(len=*),  intent(in)     :: dimname2 ! name of 2nd array dimensions of field on file (array order)
    character(len=*),  intent(in)     :: dimname3 ! name of 3rd array dimensions of field on file (array order)
    integer,           intent(in)     :: dim1b    ! start of first  dimension of array to be returned
    integer,           intent(in)     :: dim1e    ! end   of first  dimension of array to be returned
    integer,           intent(in)     :: dim2b    ! start of second dimension of array to be returned
    integer,           intent(in)     :: dim2e    ! end   of second dimension of array to be returned
    integer,           intent(in)     :: dim3b    ! start of third  dimension of array to be returned
    integer,           intent(in)     :: dim3e    ! end   of third  dimension of array to be returned
    real(r8), target,  intent(out)    :: field(dim1b:dim1e,dim2b:dim2e,dim3b:dim3e) ! array to be returned (decomposed or global)
    logical,           intent(out)    :: readvar  ! true => variable is on initial dataset
    character(len=*), optional, intent(in) :: gridname ! Name of variable's grid
    integer, optional, intent(in)     :: timelevel
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    type(io_desc_t), pointer  :: iodesc
    integer                   :: grid_map  ! grid ID for data mapping
    integer                   :: i, j, k   ! indices
    integer                   :: ierr      ! error status
    type(var_desc_t)          :: varid     ! variable id

    integer                   :: arraydimsize(3) ! field dimension lengths
    integer                   :: arraydimids(3) ! Dimension IDs
    integer                   :: permutation(3)
    logical                   :: ispermuted

    integer                   :: ndims ! number of dimensions
    integer                   :: pdims ! number of dimensions w/o timeslice
    integer                   :: dimids(PIO_MAX_VAR_DIMS) ! file variable dims
    integer                   :: dimlens(PIO_MAX_VAR_DIMS) ! file variable shape

    ! Offsets for reading global variables
    integer                   :: strt(3)     ! start lon, lev, lat indices for netcdf 3-d
    integer                   :: cnt (3)     ! lon, lat counts for netcdf 3-d
    character(len=PIO_MAX_NAME) :: tmpname

    real(r8), pointer         :: tmp3d(:,:,:) ! input data for permutation

    logical                   :: readvar_tmp ! if true, variable is on tape
    character(len=32)         :: subname='INFLD_REAL_3D_3D' ! subroutine name
    character(len=128)        :: errormsg
    character(len=PIO_MAX_NAME) :: field_dnames(3)
    character(len=PIO_MAX_NAME) :: file_dnames(4)

    ! For SCAM
    real(r8)                  :: closelat, closelon
    integer                   :: lonidx, latidx

    nullify(iodesc)

    !
    !-----------------------------------------------------------------------
    !
    !    call pio_setdebuglevel(3)

    ! Should we be using a different interface?
    if ((trim(dimname1) == trim(dimname2)) .or. (len_trim(dimname2) == 0)) then
      call infld(varname, ncid, dimname1, dimname3,                           &
           dim1b, dim1e, dim2b, dim2e, dim3b, dim3e,                          &
           field, readvar, gridname, timelevel)
    else if ((trim(dimname1) == trim(dimname3)) .or. (len_trim(dimname3) == 0)) then
      call infld(varname, ncid, dimname1, dimname2,                           &
           dim1b, dim1e, dim2b, dim2e, dim3b, dim3e,                          &
           field, readvar, gridname, timelevel)
    else

      !
      ! Error conditions
      !
      if (present(gridname)) then
        grid_map = cam_grid_id(trim(gridname))
      else
        grid_map = cam_grid_id('physgrid')
      end if
      if (.not. cam_grid_check(grid_map)) then
        if(masterproc) then
          if (present(gridname)) then
            write(errormsg, *)': invalid gridname, "',trim(gridname),'", specified for field ',trim(varname)
          else
            write(errormsg, *)': Internal error, no "physgrid" gridname'
          end if
        end if
        call endrun(trim(subname)//errormsg)
      end if

      if (debug .and. masterproc) then
        if (present(gridname)) then
          write(iulog, '(5a)') trim(subname),': field = ',trim(varname),', grid = ',trim(gridname)
        else
          write(iulog, '(4a)') trim(subname),': field = ',trim(varname),', grid = physgrid'
        end if
        call shr_sys_flush(iulog)
      end if

      !
      ! Read netCDF file
      !
      !
      ! Check if field is on file; get netCDF variable id
      !
      call cam_pio_check_var(ncid, varname, varid, ndims, dimids, dimlens,    &
           readvar_tmp, dimnames=file_dnames)
      !
      ! If field is on file:
      !
      if (readvar_tmp) then
        if (debug .and. masterproc) then
          write(iulog, '(2a,9(i0,a))') trim(subname),': field(',              &
               dim1b,':',dim1e,',',dim2b,':',dim2e,',',dim3b,':',dim3e,       &
               '), file(',dimlens(1),',',dimlens(2),',',dimlens(3),')'
          call shr_sys_flush(iulog)
        end if
        !
        ! Get array dimension id's and sizes
        !
        ierr = PIO_inq_dimid(ncid, dimname1, arraydimids(1))
        ierr = PIO_inq_dimid(ncid, dimname2, arraydimids(2))
        ierr = PIO_inq_dimid(ncid, dimname3, arraydimids(3))
        arraydimsize(1) = (dim1e - dim1b + 1)
        arraydimsize(2) = (dim2e - dim2b + 1)
        arraydimsize(3) = (dim3e - dim3b + 1)

        do j = 1, 3
          if (arraydimsize(j) /= size(field, j)) then
            write(errormsg, *) ': Mismatch between array bounds and field size for ', &
                 trim(varname), ', dimension', j
            call endrun(trim(subname)//errormsg)
          end if
        end do

        pdims = ndims
        if (ndims > 4) then
          call endrun(trim(subname)//': too many dimensions for '//trim(varname))
        else if (ndims < 3) then
          call endrun(trim(subname)//': too few dimensions for '//trim(varname))
        else
          ! Check to make sure that the fourth dimension is time
          if (ndims == 4) then
            ierr = pio_inq_dimname(ncid, dimids(4), tmpname)
            if (trim(tmpname) /= 'time') then
              call endrun(trim(subname)//': dimension mismatch for '//trim(varname))
            end if
            pdims = 3
          end if
        end if

        if(ndims == 4) then
          if(present(timelevel)) then
            call pio_setframe(ncid, varid, int(timelevel,kind=pio_offset_kind))
          end if
        end if

        field_dnames(1) = dimname1
        field_dnames(2) = dimname2
        field_dnames(3) = dimname3

        if (single_column) then	
          ! This could be generalized but for now only handles a single point
          strt(1) = dim1b
          strt(2) = dim2b
          strt(3) = dim3b
          cnt = arraydimsize
          call shr_scam_getCloseLatLon(ncid%fh,scmlat,scmlon,closelat,closelon,latidx,lonidx)
          if (trim(field_dnames(1)) == 'lon') then
            strt(1) = lonidx ! First dim always lon for Eulerian dycore
          else
            call endrun(trim(subname)//': lon should be first dimension for '//trim(varname))
          end if
          if (trim(field_dnames(2)) == 'lat') then
            strt(2) = latidx
          else if (trim(field_dnames(3)) == 'lat') then
            strt(3) = latidx
          else
            call endrun(trim(subname)//': lat dimension not found for '//trim(varname))
          end if

          ! Check for permuted dimensions ('out of order' array)
          call calc_permutation(dimids, arraydimids, permutation, ispermuted)
          if (ispermuted) then
            call cam_permute_array(strt, permutation)
            call cam_permute_array(cnt, permutation)
            allocate(tmp3d(1:cnt(1), 1:cnt(2), 1:cnt(3)))
            ierr = pio_get_var(ncid, varid, strt, cnt, tmp3d)
            do k = dim3b, dim3e
              do j = dim2b, dim2e
                do i = dim1b, dim1e
                  ! We don't need strt anymore, reuse it
                  strt(1) = i - dim1b + 1
                  strt(2) = j - dim2b + 1
                  strt(3) = k - dim3b + 1
                  call cam_permute_array(strt, permutation)
                  field(i,j,k) = tmp3d(strt(1), strt(2), strt(3))
                end do
              end do
            end do
            deallocate(tmp3d)
          else
            ierr = pio_get_var(ncid, varid, strt, cnt, field)
          end if
        else
          ! All distributed array processing
          call cam_grid_get_decomp(grid_map, arraydimsize, dimlens(1:pdims),  &
               pio_double, iodesc, field_dnames=field_dnames, file_dnames=file_dnames(1:3))
          call pio_read_darray(ncid, varid, iodesc, field, ierr)
        end if ! end of single column
      end if  ! end of readvar_tmp

      readvar = readvar_tmp

    end if ! end of call infld_real_2d_3d instead

  end subroutine infld_real_3d_3d


end module ncdio_atm



