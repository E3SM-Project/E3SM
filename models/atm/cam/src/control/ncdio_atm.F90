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
  use pio, only : pio_offset, file_desc_t, var_desc_t, pio_noerr, pio_double, &
       pio_inq_varid, pio_inq_dimlen, pio_max_name, pio_bcast_error, &
       pio_internal_error, pio_inq_dimid, pio_max_var_dims, pio_inq_vardimid, &
       pio_inq_varndims, io_desc_t, pio_setframe, pio_read_darray, pio_seterrorhandling
  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_sys_mod  , only: shr_sys_flush         ! Standardized system subroutines
  use shr_scam_mod  , only: shr_scam_getCloseLatLon  ! Standardized system subroutines
  use spmd_utils   , only: masterproc
  use abortutils   , only: endrun
  use dycore,        only: dycore_is
  use scamMod,      only: initTimeIdx,scmlat,scmlon,single_column
  use cam_logfile,  only: iulog
  !
  ! !PUBLIC TYPES:
  implicit none

  PRIVATE

  save
  public :: check_var   ! determine if variable is on PIO file
  public :: check_dim   ! validity check on dimension
  !
  !EOP
  !
  interface infld
     module procedure infld_real_2d
     module procedure infld_real_3d
     module procedure infld_real_2dncol
     module procedure infld_real_3dncol
  end interface


  public :: infld

  integer STATUS
  real(r8) surfdat
  !-----------------------------------------------------------------------

contains
  !BOP
  !
  ! !IROUTINE: check_dim
  !
  ! !INTERFACE:
  subroutine check_dim(ncid, dimname, value)
    use pio, only : pio_inq_dimid, pio_inq_dimlen
    !
    ! !DESCRIPTION: 
    ! Validity check on dimension
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer, intent(in) :: value
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    integer :: dimid, dimlen    ! temporaries
    integer :: ierr
    !-----------------------------------------------------------------------

    ierr = PIO_inq_dimid (ncid, trim(dimname), dimid)
    ierr = PIO_inq_dimlen (ncid, dimid, dimlen)
    if (dimlen /= value) then
       write(iulog,*) 'CHECK_DIM error: mismatch of input dimension ',dimlen, &
            ' with expected value ',value,' for variable ',trim(dimname)
       call endrun()
    end if

  end subroutine check_dim

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: check_var
  !
  ! !INTERFACE:
  subroutine check_var(ncid, varname, varid, readvar)
    !
    ! !DESCRIPTION: 
    ! Check if variable is on netcdf file
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout)          :: ncid
    character(len=*), intent(in) :: varname
    type(var_desc_t), intent(out)         :: varid
    logical, intent(out)         :: readvar 
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    integer :: ret     ! return value
    !-----------------------------------------------------------------------
    call pio_seterrorhandling( ncid, PIO_BCAST_ERROR)

    readvar = .true.
    ret = PIO_inq_varid (ncid, varname, varid)
    if (ret/=PIO_NOERR) then
       if (masterproc) then
          write(iulog,*)'CHECK_VAR INFO: variable ',trim(varname),' is not on file'
          call shr_sys_flush(iulog)
       end if
       readvar = .false.
    end if
    call pio_seterrorhandling( ncid, PIO_INTERNAL_ERROR)

  end subroutine check_var

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: infld_real_2d
  !
  ! !INTERFACE:
  subroutine infld_real_2d(varname ,ncid   , lonnam  ,latnam  , dim1b   , &
       dim1e   ,dim2b   ,dim2e   ,field   , readvar , &
       grid_map, timelevel)

    use dyn_grid, only: get_block_gcol_d, get_block_gcol_cnt_d, get_horiz_grid_dim_d
    use phys_grid, only: get_ncols_p, get_gcol_all_p
    use xpavg_mod, only: xpavg
    use pio, only : pio_get_var, pio_read_darray
    use cam_pio_utils, only : get_phys_decomp, get_dyn_decomp !, dyn_decomp
    !
    ! !DESCRIPTION: 
    ! Netcdf i/o of 2d initial real field from netCDF file
    !
    ! !USES
    !
    use string_utils, only: to_upper
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname  ! variable name
    type(file_desc_t)        , intent(inout)  :: ncid     ! input unit
    character(len=*), intent(in)  :: lonnam   ! name of longitude dimension of field on file
    character(len=*), intent(in)  :: latnam   ! name of latitude  dimension of field on file
    integer         , intent(in)  :: dim1b    ! start of first  dimension of array to be returned
    integer         , intent(in)  :: dim1e    ! end   of first  dimension of array to be returned
    integer         , intent(in)  :: dim2b    ! start of second dimension of array to be returned
    integer         , intent(in)  :: dim2e    ! end   of second dimension of array to be returned
    real(r8)        , intent(out) :: field(dim1b:dim1e,dim2b:dim2e) ! array to be returned (decomposed or global)
    logical         , intent(out) :: readvar  ! true => variable is on initial dataset
    character(len=*), intent(in)  :: grid_map ! flag indicating which grid to map data to
    integer, optional, intent(in) :: timelevel ! Timelevel to read
    !
    !EOP
    !
    ! !LOCAL VARIABLES:

    type(io_desc_t), pointer :: iodesc
    integer :: i                        ! index
    integer :: ier                      ! error status
    type(var_desc_t) :: varid                    ! variable id
    integer :: dimlon, dimlat           ! lon, lat, lev dimension lengths
    integer tmptype
    integer ndims                       ! number of dimensions
    integer dims(PIO_MAX_VAR_DIMS)       ! variable shape
    integer londimid, latdimid          ! Dimension ID's
    integer strt(3)                     ! start lon, lat, time indices for netcdf 2-d
    integer cnt (3)                     ! lon, lat, time counts for netcdf 2-d
    data strt/3*1/                      ! 
    data cnt /1,1,1/                    ! 2-d arrs
    real(r8), pointer :: tmp(:)       ! input data

    character(len=PIO_MAX_NAME) tmpname
    character(len= 8) :: grid_map_tmp   ! Local character string
    character(len=32) :: subname='INFLD_REAL_2D' ! subroutine name
    real(r8) closelat,closelon
    integer latidx,lonidx
    integer :: gcols(dim1b:dim1e), j
    integer :: hdim1_d, hdim2_d


    !
    !-----------------------------------------------------------------------
    !
    grid_map_tmp = trim( to_upper(trim(grid_map)) )

    !
    ! Error conditions
    !
    if ( grid_map_tmp /= 'GLOBAL' .and. grid_map_tmp /= 'PHYS' .and. grid_map_tmp /= 'DYN' ) then
       if(masterproc) then
          write(iulog,*) trim(subname), ' Error:  invalid grid-map flag specified for field ',trim(varname)
          write(iulog,*) '                            grid_map = ', grid_map_tmp
       end if
       call endrun()
    end if

    !
    ! Read netCDF file
    !
    !
    ! Check if field is on file; get netCDF variable id
    !
    call check_var(ncid, varname, varid, readvar)
    !
    ! If field is on file:
    !
    if (readvar) then

       if(present(timelevel)) then
          call pio_setframe(varid, int(timelevel,kind=PIO_OFFSET))
       else
          call pio_setframe(varid, 1_PIO_OFFSET)
       end if

       if(lonnam.eq.'ncol') then
          call get_phys_decomp(iodesc, 1,1,1,pio_double)
          call pio_read_darray(ncid, varid, iodesc, field, ier)
          return
       end if
       !
       ! Get dimension id's and sizes
       !
       ier = PIO_inq_dimid  (ncid, lonnam  , londimid)
       ier = PIO_inq_dimlen (ncid, londimid, dimlon)
       ier = PIO_inq_dimid  (ncid, latnam  , latdimid)
       ier = PIO_inq_dimlen (ncid, latdimid, dimlat)

       !
       ! Check order of dimensions in variable
       !
       ier = PIO_inq_varndims(ncid, varid, ndims)
       ier = PIO_inq_vardimid(ncid, varid, dims(1:ndims))
       if (dims(1) /= londimid .or. dims(2) /= latdimid .or. ndims > 3) then
          write(iulog,*) trim(subname), ' Error: Bad number of dims or ordering while reading field ', trim(varname)
          call endrun()
       end if
       if ( single_column ) then
          call shr_scam_getCloseLatLon(ncid%fh,scmlat,scmlon,closelat,closelon,latidx,lonidx)
          dimlon=1
          dimlat=1
          strt(1)=lonidx
          strt(2)=latidx
          strt(3)=1
          cnt(3) = 1
       endif

       if ( grid_map_tmp == 'GLOBAL' .or. single_column) then
          !
          ! Allocate memory and read variable
          !

          cnt(1) = dimlon
          cnt(2) = dimlat

          ier = PIO_GET_VAR(ncid, varid, strt, cnt, field)

       else if(grid_map_tmp == 'PHYS' ) then
          call get_phys_decomp(iodesc, 1,1,1,pio_double)
          if(.not.associated(iodesc)) then
             call endrun('error getting iodesc')
          end if

          call pio_read_darray(ncid, varid, iodesc, field, ier)

       else
          call get_horiz_grid_dim_d(hdim1_d, hdim2_d)

          if(latnam .eq. 'slat') then
             call get_dyn_decomp(iodesc, hdim1_d, hdim2_d-1,1,0,pio_double)
          else
             call get_dyn_decomp(iodesc, hdim1_d, hdim2_d,1,0,pio_double)
          end if

          if(.not.associated(iodesc)) then
             call endrun('error getting iodesc')
          end if


          call pio_read_darray(ncid, varid, iodesc, field, ier)

       end if  ! end of grid_map_tmp
    end if  ! end readvar

    return
  end subroutine infld_real_2d

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: infld_real_2dncol
  !
  ! !INTERFACE:
  subroutine infld_real_2dncol(varname ,ncid , iodesc, field, readvar, timelevel )
    !
    ! !DESCRIPTION: 
    ! Netcdf i/o of 2d initial real field from netCDF file
    !
    ! !USES
    !
    use string_utils, only: to_upper
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname  ! variable name
    type(file_desc_t)         , intent(inout)  :: ncid     ! input unit
    type(io_desc_t)                         :: iodesc
    real(r8)        , intent(out) :: field(:) ! array to be returned (decomposed or global)
    logical         , intent(out) :: readvar  ! true => variable is on initial dataset
    integer, optional, intent(in) :: timelevel
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    integer :: i                        ! index
    integer :: ier                      ! error status
    type(var_desc_T) :: varid                    ! variable id
    integer :: dimncol                  ! ncol dimension length
    integer :: maxsiz
    integer tmptype
    integer ndims                       ! number of dimensions
    integer dims(PIO_MAX_VAR_DIMS)       ! variable shape
    integer ncoldimid                     ! Dimension ID's
    integer strt(2)                     ! start ncol, time indices for netcdf 2-d
    integer cnt (2)                     ! ncol, time counts for netcdf 2-d
    logical :: readvar_tmp              ! if true, variable is on tape
    character*(PIO_MAX_NAME) tmpname
    character(len=32) :: subname='INFLD_REAL_2Dncol' ! subroutine name


    !
    ! Read netCDF file
    !
    !
    ! Check if field is on file; get netCDF variable id
    !
    call check_var(ncid, varname, varid, readvar_tmp)
    !
    ! If field is on file:
    !
    if (readvar_tmp) then

       !
       ! Get dimension id's and sizes
       !
       ier = PIO_inq_dimid  (ncid, 'ncol'  , ncoldimid)
       ier = PIO_inq_dimlen (ncid, ncoldimid, dimncol)
       !
       ! Check order of dimensions in variable
       !
       ier = PIO_inq_varndims(ncid, varid, ndims)
       ier = PIO_inq_vardimid(ncid, varid, dims(1:ndims))
       if (dims(1) /= ncoldimid .or. ndims > 2) then
          write(iulog,*) trim(subname), ' Error: dim mismatch while reading field ', trim(varname)
          call endrun(subname)
       end if

       if(present(timelevel)) then
          call pio_setframe(varid, int(timelevel,kind=PIO_OFFSET))
       else
          call pio_setframe(varid, 1_PIO_OFFSET)
       end if

       !
       ! Read variable
       !
       call pio_read_darray(ncid, varid, iodesc, field, ier)

    end if  ! end of readvar_tmp

    readvar = readvar_tmp

    return

  end subroutine infld_real_2dncol

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: infld_real_3d
  !
  ! !INTERFACE:
  subroutine infld_real_3d(varname ,ncid    ,lonnam  ,levnam  ,latnam  , &
       dim1b   ,dim1e   ,dim2b   ,dim2e   ,dim3b   , &
       dim3e   ,field   ,readvar ,grid_map, array_order_in, timelevel)
    !
    ! !DESCRIPTION: 
    ! Netcdf i/o of 3d initial real field from netCDF file
    !
    ! !USES
    !
    use string_utils, only: to_upper
    use dyn_grid, only: get_block_gcol_d, get_dyn_grid_parm, get_horiz_grid_dim_d
    use phys_grid, only: get_ncols_p, get_gcol_all_p
    use constituents,     only: qmin
    use xpavg_mod, only: xpavg
    use pio, only : pio_get_var, pio_read_darray, pio_setdebuglevel
    use cam_pio_utils, only : get_phys_decomp , get_dyn_decomp

    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname  ! variable name
    type(file_desc_t), intent(inout)  :: ncid     ! input unit
    character(len=*), intent(in)  :: lonnam   ! name of longitude dimension of field on file
    character(len=*), intent(in)  :: levnam   ! name of level     dimension of field on file
    character(len=*), intent(in)  :: latnam   ! name of latitude  dimension of field on file
    integer         , intent(in)  :: dim1b    ! start of first  dimension of array to be returned
    integer         , intent(in)  :: dim1e    ! end   of first  dimension of array to be returned
    integer         , intent(in)  :: dim2b    ! start of second dimension of array to be returned
    integer         , intent(in)  :: dim2e    ! end   of second dimension of array to be returned
    integer         , intent(in)  :: dim3b    ! start of third  dimension of array to be returned
    integer         , intent(in)  :: dim3e    ! end   of third  dimension of array to be returned
    real(r8),target , intent(out) :: field(dim1b:dim1e,dim2b:dim2e,dim3b:dim3e) ! array to be returned (decomposed or global)
    logical         , intent(out) :: readvar  ! true => variable is on initial dataset
    character(len=*), intent(in)  :: grid_map ! flag indicating which grid to map data to
    character(len=3), optional, intent(in) :: array_order_in
    integer, optional, intent(in) :: timelevel
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    type(io_desc_t), pointer :: iodesc
    integer :: i, k, lat                ! indices
    integer :: ier                      ! error status
    type(var_desc_t) :: varid                    ! variable id
    integer :: dimlon, dimlat, dimlev   ! lon, lat, lev dimension lengths
    integer tmptype
    integer ndims                       ! number of dimensions
    integer dims(PIO_MAX_VAR_DIMS)       ! variable shape
    integer londimid, latdimid, levdimid ! Dimension ID's
    integer :: strt(4)                     ! start lon, lat, time indices for netcdf 2-d
    integer :: cnt (4)                     ! lon, lat, time counts for netcdf 2-d
    data strt/4*1/                      ! 
    data cnt /1,1,1,1/                  ! 3-d arrs
    real(r8), pointer :: tmp(:), tmp3d(:,:,:) ! input data
    logical :: readvar_tmp              ! if true, variable is on tape
    character*(PIO_MAX_NAME) tmpname
    character(len= 8) :: grid_map_tmp   ! Local character string
    character(len=32) :: subname='INFLD_REAL_3D' ! subroutine name
    real(r8) closelat,closelon
    integer lonidx,latidx
    integer :: gcols(dim1b:dim1e), ncols_d, ncols, j, ii, ierr
    integer :: hdim1_d, hdim2_d
    character(len=3) :: array_order


    nullify(iodesc)

    if(present(array_order_in)) then
       array_order=array_order_in
    else	
       array_order='xzy'
    end if
    !
    !-----------------------------------------------------------------------
    !
    !    call pio_setdebuglevel(3)
    grid_map_tmp = trim( to_upper(grid_map) )
    !
    ! Error conditions
    !
    if ( grid_map_tmp /= 'GLOBAL' .and. grid_map_tmp /= 'PHYS' .and. grid_map_tmp /= 'DYN') then
       if(masterproc) then
          write(iulog,*) trim(subname), ' Error:  invalid grid-map flag specified for field ',trim(varname)
          write(iulog,*) '                            grid_map = ', grid_map_tmp
       end if
       call endrun()
    end if
    if(lonnam.eq.'ncol') then
       call infld_real_3dncolphys(varname, ncid, ncols_d, gcols, levnam, field, readvar )
       return 
    end if
    !
    ! Read netCDF file
    !
    !
    ! Check if field is on file; get netCDF variable id
    !
    call check_var(ncid, varname, varid, readvar_tmp)
    !
    ! If field is on file:
    !
    if (readvar_tmp) then
       !
       ! Get dimension id's and sizes
       !
       ier = PIO_inq_dimid  (ncid, lonnam  , londimid)
       ier = PIO_inq_dimid  (ncid, levnam  , levdimid)
       ier = PIO_inq_dimid  (ncid, latnam  , latdimid)
       ier = PIO_inq_dimlen (ncid, londimid, dimlon)
       ier = PIO_inq_dimlen (ncid, latdimid, dimlat)
       ier = PIO_inq_dimlen (ncid, levdimid, dimlev)

       if ( single_column ) then	
          call shr_scam_getCloseLatLon(ncid%fh,scmlat,scmlon,closelat,closelon,latidx,lonidx)
          dimlon=1
          dimlat=1
       endif
       !
       ! Check order of dimensions in variable, allocate memory, and read
       ! (reverse lat/lev dimensions if field is dimensioned backwards on file)
       !
       ier = PIO_inq_varndims(ncid, varid, ndims)
       ier = PIO_inq_vardimid(ncid, varid, dims(1:ndims))

       if(ndims==4) then
          if(present(timelevel)) then
             call pio_setframe(varid, int(timelevel,kind=PIO_OFFSET))
          else
             call pio_setframe(varid, 1_PIO_OFFSET)
          end if
       end if

       if     (dims(1) == londimid .and. dims(2) == latdimid .and. &
            dims(3) == levdimid .and. ndims   <= 4) then
          if(grid_map_tmp .eq. 'GLOBAL' .or. single_column) then
             if ( single_column ) then
                strt(1)=lonidx
                strt(2)=latidx
                strt(3)=1
                strt(4)=1
                cnt(4) = 1
             endif
             cnt(1) = dimlon
             cnt(2) = dimlat
             cnt(3) = dimlev
             if(array_order.eq.'xzy') then
                allocate(tmp3d(dimlon,dimlat,dimlev))
             else
                tmp3d => field
             end if
             ierr = pio_get_var(ncid, varid, strt, cnt, tmp3d)

             if(array_order.eq.'xzy') then

                do lat=dim3b,dim3e
                   do k=dim2b,dim2e
                      do i=dim1b,dim1e
                         field(i,k,lat) = tmp3d(i,lat-dim3b+1,k)
                      end do
                   end do
                end do
                deallocate(tmp3d)
             end if
          else if(grid_map_tmp .eq. 'PHYS') then
             call get_phys_decomp(iodesc, 1,dimlev,1,pio_double)

             call pio_read_darray(ncid, varid, iodesc, field, ierr)

          else 
             call get_horiz_grid_dim_d(hdim1_d, hdim2_d)
             if(latnam .eq. 'slat') then
                call get_dyn_decomp(iodesc, hdim1_d,hdim2_d-1,dimlev,0,pio_double)
             else
                call get_dyn_decomp(iodesc, hdim1_d,hdim2_d,dimlev,0, pio_double)
             end if

             call pio_read_darray(ncid, varid, iodesc, field, ierr)
          end if

       elseif (dims(1) == londimid .and. dims(2) == levdimid .and. &
            dims(3) == latdimid .and. ndims   <= 4) then

          if(grid_map_tmp .eq. 'GLOBAL' .or. single_column) then
             if ( single_column ) then
                strt(1)=lonidx
                strt(2)=1
                strt(3)=latidx
                strt(4)=1
                cnt(4) = 1
             endif
             cnt(1) = dimlon
             cnt(2) = dimlev
             cnt(3) = dimlat

             ierr = pio_get_var(ncid, varid, strt, cnt, field)

          else if(grid_map_tmp .eq. 'PHYS') then
             call get_phys_decomp(iodesc, 1, dimlev, 1, pio_double, 'xzy')
             call pio_read_darray(ncid, varid, iodesc, field, ierr)
          else

             call get_horiz_grid_dim_d(hdim1_d, hdim2_d)
             if(latnam .eq. 'slat') then
                call get_dyn_decomp(iodesc, hdim1_d,hdim2_d-1,dimlev,0,pio_double, fileorder_in='xzy')
             else
                call get_dyn_decomp(iodesc, hdim1_d,hdim2_d,dimlev,0, pio_double, fileorder_in='xzy')
             end if
             call pio_read_darray(ncid, varid, iodesc, field, ierr)
          end if
       else
          if(masterproc) then
             write(iulog,*) trim(subname), ' Error: Bad number of dims or ordering while reading field ', trim(varname)
             call endrun()
          end if
       end if
    end if  ! end of readvar_tmp

    readvar = readvar_tmp

    return

  end subroutine infld_real_3d
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: infld_real_3dncol
  !
  ! !INTERFACE:
  subroutine infld_real_3dncol(varname ,ncid , iodesc, ncols, levnam, field, readvar, timelevel )
    !
    ! !DESCRIPTION: 
    ! Netcdf i/o of 3d initial real field from netCDF file
    !
    ! !USES
    !
    use string_utils, only: to_upper

    !
    ! !ARGUMENTS:
    implicit none
    integer :: ncols
    character(len=*), intent(in)      :: varname  ! variable name
    type(file_desc_t), intent(inout)  :: ncid     ! input unit
    type(io_desc_t)                   :: iodesc
    character(len=*), intent(in)      :: levnam   ! name of level     dimension of field on file
    real(r8),  intent(out)     :: field(:,:) ! array to be returned (decomposed or global)
    logical         , intent(out)     :: readvar  ! true => variable is on initial dataset
    integer, optional, intent(in) :: timelevel
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    integer :: i                        ! index
    integer :: ier                      ! error status
    type(var_desc_t) :: varid                    ! variable id
    integer :: dimncol                  ! ncol dimension length
    integer :: dimlev                   ! level dimension length
    integer :: maxsiz
    integer tmptype
    integer ndims                       ! number of dimensions
    integer dims(PIO_MAX_VAR_DIMS)       ! variable shape
    integer ncoldimid                     ! Dimension ID's
    integer levdimid                     ! Dimension ID's
    real(r8), pointer :: tmp(:)       ! input data
    logical :: readvar_tmp              ! if true, variable is on tape
    character*(PIO_MAX_NAME) tmpname
    character(len=32) :: subname='INFLD_REAL_3Dncol' ! subroutine name
    integer :: lsize
    !
    ! Read netCDF file
    !
    !
    ! Check if field is on file; get netCDF variable id
    !
    call check_var(ncid, varname, varid, readvar_tmp)
    !
    ! If field is on file:
    !
    if (readvar_tmp) then

       !
       ! Get dimension id's and sizes
       !
       ier = PIO_inq_dimid  (ncid, 'ncol'  , ncoldimid)
       ier = PIO_inq_dimlen (ncid, ncoldimid, dimncol)
       ier = PIO_inq_dimid  (ncid, levnam  , levdimid)
       ier = PIO_inq_dimlen (ncid, levdimid, dimlev)

       !
       ! Check order of dimensions in variable
       !
       ier = PIO_inq_varndims(ncid, varid, ndims)
       ier = PIO_inq_vardimid(ncid, varid, dims(1:ndims))
       if (dims(1) /= ncoldimid .or. dims(2) /= levdimid .or. ndims > 3) then
          write(iulog,*) trim(subname), ' Error: dim mismatch while reading field ', trim(varname)
          call endrun(subname)
       end if
       !
       ! Read variable
       !

       if(present(timelevel)) then
          call pio_setframe(varid, int(timelevel,kind=PIO_OFFSET))
       else
          call pio_setframe(varid, 1_PIO_OFFSET)
       end if


       call pio_read_darray(ncid, varid, iodesc, field, ier)

    end if  ! end of readvar_tmp

    readvar = readvar_tmp

    return

  end subroutine infld_real_3dncol
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: infld_real_3dncol
  !
  ! !INTERFACE:
  subroutine infld_real_3dncolphys(varname ,ncid , ncols, gcol, levnam, field, readvar , timelevel)
    !
    ! !DESCRIPTION: 
    ! Netcdf i/o of 3d initial real field from netCDF file
    !
    ! !USES
    !
    use phys_grid, only : get_ncols_p
    use ppgrid, only : begchunk, endchunk
    use string_utils, only: to_upper
    use cam_pio_utils, only : get_phys_decomp
    use pio, only : pio_freedecomp
    !
    ! !ARGUMENTS:
    implicit none
    integer :: ncols
    integer :: gcol(:)
    character(len=*), intent(in)  :: varname  ! variable name
    type(file_desc_t), intent(inout)  :: ncid     ! input unit
    character(len=*), intent(in)  :: levnam   ! name of level     dimension of field on file
    real(r8), intent(out) :: field(:,:,:) ! array to be returned (decomposed or global)
    logical         , intent(out) :: readvar  ! true => variable is on initial dataset
    integer, optional, intent(in) :: timelevel
    !
    !EOP
    !
    ! !LOCAL VARIABLES:
    type(io_desc_t), pointer :: iodesc
    integer :: i                        ! index
    integer :: ier                      ! error status
    type(var_desc_t) :: varid                    ! variable id
    integer :: dimncol                  ! ncol dimension length
    integer :: dimlev                   ! level dimension length
    integer :: maxsiz
    integer tmptype
    integer ndims                       ! number of dimensions
    integer dims(PIO_MAX_VAR_DIMS)       ! variable shape
    integer ncoldimid                     ! Dimension ID's
    integer levdimid                     ! Dimension ID's
    real(r8), pointer :: tmp(:)       ! input data
    logical :: readvar_tmp              ! if true, variable is on tape
    character*(PIO_MAX_NAME) tmpname
    character(len=32) :: subname='INFLD_REAL_3Dncol' ! subroutine name
    integer :: ii, j, k, lchnk
    !
    ! Read netCDF file
    !
    !
    ! Check if field is on file; get netCDF variable id
    !
    call check_var(ncid, varname, varid, readvar_tmp)
    !
    ! If field is on file:
    !
    if (readvar_tmp) then
       !
       ! Get dimension id's and sizes
       !
       ier = PIO_inq_dimid  (ncid, 'ncol'  , ncoldimid)
       ier = PIO_inq_dimlen (ncid, ncoldimid, dimncol)
       ier = PIO_inq_dimid  (ncid, levnam  , levdimid)
       ier = PIO_inq_dimlen (ncid, levdimid, dimlev)
       call get_phys_decomp(iodesc, 1,dimlev,1,pio_double)

       !
       ! Check order of dimensions in variable
       !
       ier = PIO_inq_varndims(ncid, varid, ndims)
       ier = PIO_inq_vardimid(ncid, varid, dims(1:ndims))
       if (dims(1) /= ncoldimid .or. dims(2) /= levdimid .or. ndims > 3) then
          write(iulog,*) trim(subname), ' Error: dim mismatch while reading field ', trim(varname)
          call endrun(subname)
       end if
       !
       ! Read variable
       !


       if(present(timelevel)) then
          call pio_setframe(varid, int(timelevel,kind=PIO_OFFSET))
       else
          call pio_setframe(varid, 1_PIO_OFFSET)
       end if

       call pio_read_darray(ncid, varid, iodesc, field, ier)

    end if  ! end of readvar_tmp

    readvar = readvar_tmp

    return

  end subroutine infld_real_3dncolphys


end module ncdio_atm



