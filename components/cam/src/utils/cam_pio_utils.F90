! Utility functions in support of PIO io interface
module cam_pio_utils

  use pio,          only: io_desc_t, iosystem_desc_t, file_desc_t, var_desc_t
  use pio,          only: pio_freedecomp, pio_rearr_subset
  use shr_kind_mod, only: r8=>shr_kind_r8
  use cam_logfile,  only: iulog
  use perf_mod,     only: t_startf, t_stopf
  use spmd_utils,   only: masterproc

  implicit none
  private
  save

  public :: cam_pio_createfile   ! Create a new NetCDF file for PIO writing
  public :: cam_pio_openfile     ! Open an existing NetCDF file
  public :: cam_pio_closefile    ! Close an open PIO file handle
  public :: cam_pio_fileexists   ! Check if file exists
  public :: cam_pio_newdecomp    ! Create a new PIO decompsition (mapping)
  public :: init_pio_subsystem   ! called from cam_comp
  public :: cam_pio_get_decomp   ! Find an existing decomp or create a new one
  public :: cam_pio_handle_error ! If error, print a custom error message

  public :: cam_permute_array
  public :: calc_permutation

  ! Convenience interfaces
  public :: cam_pio_def_dim
  public :: cam_pio_def_var
  public :: cam_pio_get_var

  ! General utility
  public :: cam_pio_var_info
  public :: cam_pio_find_var
  public :: cam_pio_check_var

  public :: clean_iodesc_list

  integer            :: pio_iotype
  integer, parameter :: pio_rearranger = pio_rearr_subset

  ! This variable should be private ?
  type(iosystem_desc_t), pointer, public :: pio_subsystem => null()

  ! Some private string length parameters
  integer, parameter :: errormsg_str_len = 128

  ! The iodesc_list allows us to cache existing PIO decompositions
  ! The tag needs the dim lengths, the dtype and map id (+ optional permutation)
  integer, parameter      :: tag_len           = 48
  type iodesc_list
    character(tag_len)          :: tag
    type(io_desc_t),    pointer :: iodesc => NULL()
    type(iodesc_list),  pointer :: next => NULL()
  end type iodesc_list

  type(iodesc_list), target :: iodesc_list_top

  ! Create a special type to hold a var_desc_t pointer so we can have an
  ! array of them
  type, public :: vdesc_ptr
    type(var_desc_t), pointer :: vd => NULL()
  end type vdesc_ptr

  interface cam_pio_def_var
    module procedure cam_pio_def_var_0d
    module procedure cam_pio_def_var_md
  end interface

  interface cam_pio_get_var
    module procedure cam_pio_get_var_2d_r8
    module procedure cam_pio_get_var_2d_r8_perm
    module procedure cam_pio_get_var_3d_r8
    module procedure cam_pio_get_var_3d_r8_perm
  end interface

  interface calc_permutation
    module procedure calc_permutation_int
    module procedure calc_permutation_char
  end interface

  interface cam_permute_array
    module procedure permute_array_int
    module procedure permute_array_r8
  end interface

contains

  ! use_scam_limits is a private interface used to gather information about
  !    single-column usage and limits for use by the cam_pio_get_var interfaces
  ! This still only works for lat/lon dycores
  logical function use_scam_limits(File, start, kount, dimnames)
    use shr_scam_mod,   only: shr_scam_getCloseLatLon
    use scamMod,        only: scmlat, scmlon, single_column
    use cam_abortutils, only: endrun

    ! Dummy arguments
    type(file_desc_t),          intent(inout) :: File
    integer,                    intent(inout) :: start(:)
    integer,                    intent(inout) :: kount(:)
    character(len=*), optional, intent(in)    :: dimnames(:)

    ! Local variables
    character(len=*),   parameter             :: subname='USE_SCAM_LIMITS'
    real(r8)                                  :: closelat, closelon
    integer                                   :: latidx,   lonidx
    integer                                   :: i
    logical                                   :: latfound

    use_scam_limits = single_column
    if (use_scam_limits) then
      call shr_scam_getCloseLatLon(File, scmlat, scmlon, closelat, closelon,  &
           latidx, lonidx)
      if (present(dimnames)) then
        if (trim(dimnames(1)) == 'lon') then
          start(1) = lonidx ! First dim always lon for Eulerian dycore
          ! This could be generalized -- for now, stick with single column
          kount(1) = 1
        else
          call endrun(trim(subname)//': lon should be first dimension')
        end if
        latfound = .false.
        do i = 2, size(dimnames)
          if (size(start) < i) then
            call endrun(trim(subname)//': start too small')
          end if
          if (trim(dimnames(i)) == 'lat') then
            start(i) = latidx
            ! This could be generalized -- for now, stick with single column
            kount(i) = 1
            latfound = .true.
          end if
        end do
        if (.not. latfound) then
          call endrun(trim(subname)//': lat dimension not found')
        end if
      else
        ! No dimnames, assume standard positions (lon,lat)
        start(1) = lonidx
        start(2) = latidx
        ! This could be generalized -- for now, stick with single column
        kount(1:2) = 1
      end if
    end if

  end function use_scam_limits

  ! calc_permutation: Calculate a permutation array if filedims and arraydims
  !                   are in a different order
  !    E.g.: If filedims is  (lon, lat, lev, time) and
  !             arraydims is (lon, lev, lat), then
  !             perm is      (1,   3,   2) and isperm is set to .true.
  subroutine calc_permutation_int(filedims, arraydims, perm, isperm)
    use cam_abortutils,   only: endrun

    ! Dummy variables
    integer,         intent(in)  :: filedims(:)
    integer,         intent(in)  :: arraydims(:)
    integer,         intent(out) :: perm(:)
    logical,         intent(out) :: isperm

    ! Local variables
    character(len=*), parameter  :: subname='CALC_PERMUTATION_INT'
    integer                      :: i, j
    integer                      :: adims, fdims

    perm = 0
    isperm = .false.
    adims = size(arraydims)
    fdims = size(filedims)

    if (size(perm) < adims) then
      call endrun(trim(subname)//': perm smaller than arraydims')
    end if

    if (fdims < adims) then
      call endrun(trim(subname)//': filedims smaller than arraydims')
    end if

    do i = 1, adims
      if (arraydims(i) == filedims(i)) then
        perm(i) = i
      else
        isperm = .true.
        do j = 1, fdims
          if (arraydims(i) == filedims(j)) then
            perm(i) = j
            exit
          else if (j == fdims) then
            call endrun(trim(subname)//': No match for array dimension')
          ! No else, just try the next j index
          end if
        end do
      end if
    end do

  end subroutine calc_permutation_int

  subroutine calc_permutation_char(filedims, arraydims, perm, isperm)
    use cam_abortutils,   only: endrun

    ! Dummy variables
    character(len=*),   intent(in)  :: filedims(:)
    character(len=*),   intent(in)  :: arraydims(:)
    integer,            intent(out) :: perm(:)
    logical,            intent(out) :: isperm

    ! Local variables
    character(len=*),   parameter   :: subname='CALC_PERMUTATION_CHAR'
    integer                         :: i, j
    integer                         :: adims, fdims

    perm = 0
    isperm = .false.
    adims = size(arraydims)
    fdims = size(filedims)

    if (size(perm) < adims) then
      call endrun(trim(subname)//': perm smaller than arraydims')
    end if

    if (fdims < adims) then
      call endrun(trim(subname)//': filedims smaller than arraydims')
    end if

    ILOOP : do i = 1, adims
      if (trim(arraydims(i)) == trim(filedims(i))) then
        perm(i) = i
      else
        isperm = .true.
        do j = 1, fdims
          if (trim(arraydims(i)) == trim(filedims(j))) then
            perm(i) = j
            exit
          else if (j == fdims) then
            ! We have no match but for character strings, just say no perm
            isperm = .false.
            exit ILOOP
          ! No else, just try the next j index
          end if
        end do
      end if
    end do ILOOP

  end subroutine calc_permutation_char

  subroutine permute_array_int(array, perm)

    ! Dummy arguments
    integer, intent(inout) :: array(:)
    integer, intent(in)    :: perm(:)

    ! Local variables
    integer, allocatable   :: temp(:)
    integer                :: nelem, i

    nelem = size(array)
    allocate(temp(nelem))
    temp = array
    do i = 1, nelem
      array(i) = temp(perm(i))
    end do

    deallocate(temp)
  end subroutine permute_array_int

  subroutine permute_array_r8(array, perm)

    ! Dummy arguments
    real(r8), intent(inout) :: array(:)
    integer,  intent(in)    :: perm(:)

    ! Local variables
    real(r8), allocatable   :: temp(:)
    integer                 :: nelem, i

    nelem = size(array)
    allocate(temp(nelem))
    temp = array
    do i = 1, nelem
      array(i) = temp(perm(i))
    end do

    deallocate(temp)
  end subroutine permute_array_r8

  subroutine cam_pio_handle_error(ierr, errorstr)
    use cam_abortutils,   only: endrun
    use pio,          only: pio_noerr

    ! Dummy arguments
    integer,          intent(in)  :: ierr
    character(len=*), intent(in)  :: errorstr

    ! Local variables
    character(len=256) :: errormsg

    if (ierr /= PIO_NOERR) then
      write(errormsg, '(a,i6,2a)') '(PIO:', ierr, ') ', trim(errorstr)
      call endrun(errormsg)
    end if
    
  end subroutine cam_pio_handle_error

  !-----------------------------------------------------------------------
  !
  ! cam_pio_var_info: Retrieve variable properties
  !
  !-----------------------------------------------------------------------
  subroutine cam_pio_var_info(ncid, varid, ndims, dimids, dimlens, dimnames, varname)
    use pio,        only: PIO_inq_varndims, PIO_inq_vardimid, PIO_inq_dimlen, PIO_inq_dimname
    use pio,        only: PIO_seterrorhandling, PIO_BCAST_ERROR
    use cam_abortutils, only: endrun


    ! Dummy arguments
    type(file_desc_t),           intent(inout) :: ncid
    type(var_desc_t),            intent(in)    :: varid
    integer,                     intent(out)   :: ndims
    integer,                     intent(out)   :: dimids(:)
    integer,                     intent(out)   :: dimlens(:)
    character(len=*),  optional, intent(out)   :: dimnames(:)
    character(len=*),  optional, intent(in)    :: varname

    ! Local variables
    integer                                    :: ret     ! PIO return value
    integer                                    :: i
    integer                                    :: err_handling
    character(len=128)                         :: errsuff
    !-----------------------------------------------------------------------
    ! We will handle errors for this routine
    !!XXgoldyXX: This hack should be replaced with the PIO interface
    !err_handling = ncid%iosystem%error_handling !! Hack
    call PIO_seterrorhandling(ncid, PIO_BCAST_ERROR, err_handling)

    dimids = -1
    ndims = 0
    dimlens = 0

    if (present(varname)) then
      errsuff = ' for '//trim(varname)
    else
      errsuff = ''
    end if
    ! Check dimensions
    ret = PIO_inq_varndims(ncid, varid, ndims)
    call cam_pio_handle_error(ret, 'CAM_PIO_VAR_INFO: Error with num dimensions')
    if (size(dimids) < ndims) then
      call endrun('CAM_PIO_VAR_INFO: dimids too small'//trim(errsuff))
    end if
    ret = PIO_inq_vardimid(ncid, varid, dimids(1:ndims))
    call cam_pio_handle_error(ret, 'CAM_PIO_VAR_INFO: Error with inq dim ids'//trim(errsuff))
    if (size(dimlens) < ndims) then
      call endrun('CAM_PIO_VAR_INFO: dimlens too small'//trim(errsuff))
    end if
    do i = 1, ndims
      ret = PIO_inq_dimlen(ncid, dimids(i), dimlens(i))
      call cam_pio_handle_error(ret, 'CAM_PIO_VAR_INFO: Error with inq dimlens')
      if (present(dimnames)) then
        ret = PIO_inq_dimname(ncid, dimids(i), dimnames(i))
        call cam_pio_handle_error(ret, 'CAM_PIO_VAR_INFO: Error with inq dimnames')
      end if
    end do
    call PIO_seterrorhandling(ncid, err_handling)

  end subroutine cam_pio_var_info

  subroutine cam_pio_find_var(ncid, varname, varid, found)
    use pio,            only: pio_inq_varid, pio_noerr
    use pio,            only: PIO_seterrorhandling, PIO_BCAST_ERROR

    ! Dummy arguments
    type(file_desc_t),           intent(inout) :: ncid
    character(len=*),            intent(in)    :: varname
    type(var_desc_t),            intent(out)   :: varid
    logical,                     intent(out)   :: found

    ! Local variables
    integer                                    :: ret     ! PIO return value
    integer                                    :: err_handling

    !-----------------------------------------------------------------------
    ! We will handle errors for this routine
    !!XXgoldyXX: This hack should be replaced with the PIO interface
    !err_handling = ncid%iosystem%error_handling !! Hack
    call PIO_seterrorhandling(ncid, PIO_BCAST_ERROR, err_handling)
    ret = PIO_inq_varid(ncid, trim(varname), varid)
    found = (ret == PIO_NOERR)
    call PIO_seterrorhandling(ncid, err_handling)

  end subroutine cam_pio_find_var


  !-----------------------------------------------------------------------
  !
  ! cam_pio_check_var: Make sure var exists and retrieve properties
  !
  !-----------------------------------------------------------------------
  subroutine cam_pio_check_var(ncid, varname, varid, ndims, dimids, dimlens,  &
       readvar, dimnames)
    use pio,         only: PIO_inq_varid, PIO_NOERR
    use pio,         only: PIO_seterrorhandling, PIO_BCAST_ERROR
    use shr_sys_mod, only: shr_sys_flush ! Standardized system subroutines

    ! Dummy arguments
    type(file_desc_t),          intent(inout) :: ncid
    character(len=*),           intent(in)    :: varname
    type(var_desc_t),           intent(out)   :: varid
    integer,                    intent(out)   :: ndims
    integer,                    intent(out)   :: dimids(:)
    integer,                    intent(out)   :: dimlens(:)
    logical,                    intent(out)   :: readvar
    character(len=*), optional, intent(out)   :: dimnames(:)

    ! Local variables
    integer                                   :: ret     ! PIO return value
    integer                                   :: err_handling

    !-----------------------------------------------------------------------
    ! We will handle errors for this routine
    !!XXgoldyXX: This hack should be replaced with the PIO interface
    !err_handling = ncid%iosystem%error_handling !! Hack
    call pio_seterrorhandling(ncid, PIO_BCAST_ERROR, err_handling)

    dimids = -1
    ndims = 0
    dimlens = 0
    ret = PIO_inq_varid(ncid, trim(varname), varid)
    if (ret /= PIO_NOERR) then
      readvar = .false.
      if (masterproc) then
        write(iulog,*)'CAM_PIO_CHECK_VAR INFO: variable ',trim(varname),' is not on file'
        call shr_sys_flush(iulog)
      end if
    else
      readvar = .true.
      call cam_pio_var_info(ncid, varid, ndims, dimids, dimlens,              &
           dimnames=dimnames, varname=varname)
    end if
    call pio_seterrorhandling(ncid, err_handling)

  end subroutine cam_pio_check_var

  subroutine init_pio_subsystem(nlfilename)
    use shr_pio_mod,   only: shr_pio_getiosys, shr_pio_getiotype
    use cam_instance, only: atm_id

    ! Dummy argument
    character(len=*) nlfilename

    pio_subsystem => shr_pio_getiosys(atm_id)
    pio_iotype =  shr_pio_getiotype(atm_id)

  end subroutine init_pio_subsystem

  ! cam_pio_get_decomp: retrieve or create a PIO decomposition for the field
  !                     described by ldims and dtype where dims is the field's
  !                     local shape.
  !                     fdims is the shape of the field in a NetCDF file.
  !                     map describes the mapping of the distributed dimensions
  !                     field_dist_in is used if the dimensions of the
  !                        field array are not in map order
  !                     file_dist_in is used if the dimensions of the
  !                        field on file are not in map order
  !                     
  subroutine cam_pio_get_decomp(iodesc, ldims, fdims, dtype, map,             &
       field_dist_in, file_dist_in, permute)
    use pio,            only: pio_offset_kind
    use cam_abortutils, only: endrun
    use cam_map_utils,  only: cam_filemap_t

    ! Dummy arguments
    type(io_desc_t),           pointer             :: iodesc   ! intent(out)
    integer,                           intent(in)  :: ldims(:) ! Local array
    integer,                           intent(in)  :: fdims(:) ! File dims
    integer,                           intent(in)  :: dtype
    type(cam_filemap_t),       target, intent(in)  :: map
    integer,         optional,         intent(in)  :: field_dist_in(:)
    integer,         optional,         intent(in)  :: file_dist_in(:)
    integer,         optional,         intent(in)  :: permute(:)

    ! Local variables
    logical                                        :: found
    integer(PIO_OFFSET_KIND),  pointer             :: dof(:)
    type(iodesc_list),         pointer             :: iodesc_p

    call t_startf('get_decomp')

    nullify(iodesc_p)
    nullify(dof)
    call find_iodesc(ldims, fdims, dtype, map, iodesc_p, found, perm=permute)

    if (.not. found) then
      ! Create a new iodesc
      if(masterproc) then
        write(iulog,*) 'Creating new decomp: ', iodesc_p%tag
      end if

      call t_startf('get_filemap')
      call map%get_filemap(ldims, fdims, dof,                                 &
           src_in=field_dist_in, dest_in=file_dist_in, permutation_in=permute)
      call t_stopf('get_filemap')
      if (associated(iodesc_p%iodesc)) then
        ! Quick sanity check
        call endrun('cam_pio_get_decomp: iodesc already allocated')
      end if
      allocate(iodesc_p%iodesc)
      call t_startf('newdecomp')
      call cam_pio_newdecomp(iodesc_p%iodesc, fdims, dof, dtype)
      call t_stopf('newdecomp')

      deallocate(dof)
      nullify(dof)
    end if
    ! At this point, we should have a decomp, assign iodesc
    iodesc => iodesc_p%iodesc
    nullify(iodesc_p)

    call t_stopf('get_decomp')

  end subroutine cam_pio_get_decomp

  subroutine cam_pio_newdecomp(iodesc, dims, dof, dtype)
    use pio,          only: pio_initdecomp, pio_offset_kind

    type(io_desc_t),          pointer              :: iodesc
    integer,                           intent(in)  :: dims(:)
    integer(kind=PIO_OFFSET_KIND),     intent(in)  :: dof(:)
    integer,                           intent(in)  :: dtype

    call pio_initdecomp(pio_subsystem, dtype, dims, dof, iodesc,              &
         rearr=pio_rearranger)

  end subroutine cam_pio_newdecomp

  subroutine find_iodesc(ldimlens, fdimlens, dtype, map, iodesc_p, found, perm)
    use cam_abortutils,    only: endrun
    use cam_map_utils,     only: cam_filemap_t

    ! Dummy arguments
    integer,                    intent(in)    :: ldimlens(:)
    integer,                    intent(in)    :: fdimlens(:)
    integer,                    intent(in)    :: dtype
    type(cam_filemap_t),        intent(in)    :: map
    type(iodesc_list), pointer                :: iodesc_p
    logical,                    intent(out)   :: found
    integer, optional,          intent(in)    :: perm(:)

    ! Local variables
    type(iodesc_list),      pointer :: curr, prev
    integer                         :: i
    integer                         :: lcnt
    integer                         :: fcnt
    integer                         :: mapind
    integer                         :: nperm
    character(len=128)              :: form
    character(len=tag_len)          :: tag
    character(len=*), parameter     :: formc = 'i0,"(i0,""!""),""!"",",'
    character(len=*), parameter     :: forme = '"""d"",i0,""!i"",i0,""!"""'
    character(len=*), parameter     :: form2 = '("(",'//formc//formc//forme//',")")'
    character(len=*), parameter     :: form3 = '("(",'//formc//formc//formc//forme//',")")'

    found = .false.
    curr => iodesc_list_top

    ! Retrieve the (hopefully) unique tag for this iodesc
    ! If a decomp was created using an earlier version of the map (hey, that
    ! might happen), we won't find it using this search because the current
    ! index is part of the search tag
    mapind = map%get_index()
    lcnt = size(ldimlens)
    fcnt = size(fdimlens)
    if (present(perm)) then
      if (size(perm) /= lcnt) then
        write(form, '(i0,a,i0)')  size(perm), ', should be ', lcnt
        call endrun('FIND_IODESC: perm has wrong size, '//form)
      end if
      nperm = lcnt
    else
      nperm = 0
    end if
    if (present(perm)) then
      write(form, form3) lcnt, fcnt, nperm
      write(tag, form) (ldimlens(i),i=1,lcnt), (fdimlens(i),i=1,fcnt), (perm(i),i=1,lcnt), dtype, mapind
    else
      write(form, form2) lcnt, fcnt
      write(tag, form) (ldimlens(i),i=1,lcnt), (fdimlens(i),i=1,fcnt), dtype, mapind
    end if

    do while(associated(curr) .and. (.not. found))
      if(trim(tag) == trim(curr%tag)) then
        found  =  .true.
        iodesc_p => curr
      else
        prev => curr
        curr => curr%next
      end if
    end do
    if(.not. found) then
      ! We didn't find a match, make sure there is an unused iodesc_list
      !    object at the end of the list for the new decomp to be stored
      curr => prev
      if(associated(curr%iodesc)) then
        allocate(curr%next)
        curr => curr%next
        nullify(curr%iodesc) ! Should already be null but . . .
        nullify(curr%next)   ! Should already be null but . . .
      end if
      ! This should be an unused object at the end of the list
      curr%tag = tag
      iodesc_p => curr
    end if
!    if(masterproc) write(iulog,*) 'Using decomp: ',curr%tag
    
  end subroutine find_iodesc


  ! cam_pio_def_dim: Define a NetCDF dimension using the PIO interface
  subroutine cam_pio_def_dim(File, name, size, dimid, existOK)
    use cam_abortutils,   only: endrun
    use pio, only: pio_inq_dimid, pio_def_dim, pio_inq_dimlen, PIO_NOERR
    use pio, only: PIO_seterrorhandling, PIO_BCAST_ERROR

    ! Dummy arguments
    type(file_desc_t),      intent(inout)  :: File    ! PIO file Handle
    character(len=*),       intent(in)     :: name    ! Dimension name
    integer,                intent(in)     :: size    ! Dimension length
    integer,                intent(out)    :: dimid   ! NetCDF dimension ID
    logical, optional,      intent(in)     :: existOK ! OK if dim defined

    ! Local variables
    logical                                :: ok_if_dim_exists
    integer                                :: ierr
    integer                                :: err_handling
    integer                                :: dimlen
    character(len=errormsg_str_len)        :: errormsg
    character(len=*), parameter            :: subname = 'cam_pio_def_dim'

    if (present(existOK)) then
      ok_if_dim_exists = existOK
    else
      ok_if_dim_exists = .false.
    end if

    ! We will handle errors for this routine
    !!XXgoldyXX: This hack should be replaced with the PIO interface
    !err_handling = File%iosystem%error_handling !! Hack
    call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)

    ierr = pio_inq_dimid(File, trim(name), dimid)
    if (ierr == PIO_NOERR) then
      if (.not. ok_if_dim_exists) then
        write(errormsg, *) ': A dimension already exists for ', trim(name)
        call endrun(trim(subname)//errormsg)
      else
        ! It is OK for the dimension to exist but it better have the same size
        ierr = pio_inq_dimlen(File, dimid, dimlen)
        if (ierr /= PIO_NOERR) then
          write(errormsg, '(2a,i0,2a)') trim(subname), ': Error ', ierr,      &
               ' finding dimension length for ', trim(name)
          call endrun(errormsg)
        else if (dimlen /= size) then
          write(errormsg, '(3a,2(i0,a))') ': Size mismatch for dimension, ',  &
               trim(name), ': ', dimlen, ' (current), ', size, ' (desired)'
          call endrun(trim(subname)//errormsg)
          ! No else, existing dimension is OK
        end if
      end if
    else
      ! inq_dimid returned an error, define the dimension
      ierr = pio_def_dim(File, trim(name), size, dimid)
      call cam_pio_handle_error(ierr, trim(subname)//': Unable to define dimension '//trim(name))
    end if

    ! Back to whatever error handling was running before this routine
    call pio_seterrorhandling(File, err_handling)

  end subroutine cam_pio_def_dim

  ! cam_pio_def_var_0d: Define a NetCDF variable using the PIO interface
  subroutine cam_pio_def_var_0d(File, name, dtype, vardesc, existOK)

    ! Dummy arguments
    type(file_desc_t),         intent(inout) :: File      ! PIO file Handle
    character(len=*),          intent(in)    :: name      ! Variable name
    integer,                   intent(in)    :: dtype     ! e.g., pio_int
    type(var_desc_t),          intent(inout) :: vardesc   ! Variable descriptor
    logical,          optional, intent(in)   :: existOK   ! OK if var defined

    ! Local variables
    integer            :: dimids(0)

    call cam_pio_def_var(File, trim(name), dtype, dimids, vardesc, existOK)
  end subroutine cam_pio_def_var_0d

  ! cam_pio_def_var_md: Define a NetCDF variable using the PIO interface
  subroutine cam_pio_def_var_md(File, name, dtype, dimids, vardesc, existOK)
    use cam_abortutils,   only: endrun
    use pio, only: pio_inq_varid, pio_def_var, PIO_NOERR
    use pio, only: PIO_seterrorhandling, PIO_BCAST_ERROR

    ! Dummy arguments
    type(file_desc_t),          intent(inout) :: File      ! PIO file Handle
    character(len=*),           intent(in)    :: name      ! Variable name
    integer,                    intent(in)    :: dtype     ! e.g., pio_int
    integer,                    intent(in)    :: dimids(:) ! NetCDF dim IDs
    type(var_desc_t),           intent(inout) :: vardesc   ! Var descriptor
    logical, optional,          intent(in)    :: existOK   ! OK if var defined

    ! Local variables
    integer                                   :: ierr
    integer                                   :: err_handling
    logical                                   :: ok_if_var_exists
    character(len=errormsg_str_len)           :: errormsg
    character(len=*), parameter               :: subname = 'cam_pio_def_var'

    if (present(existOK)) then
      ok_if_var_exists = existOK
    else
      ok_if_var_exists = .false.
    end if

    ! We will handle errors for this routine
    !!XXgoldyXX: This hack should be replaced with the PIO interface
    !err_handling = File%iosystem%error_handling !! Hack
    call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)

    ! Check to see if the variable already exists in the file
    ierr = pio_inq_varid(File, name, vardesc)
    if (ierr == PIO_NOERR) then
      if (.not. ok_if_var_exists) then
        write(errormsg, *) ': A variable already exists for ', trim(name)
        call endrun(trim(subname)//errormsg)
      end if
    else
      ! OK to define the variable
      if (size(dimids) > 0) then
        ierr = pio_def_var(File, trim(name), dtype, dimids, vardesc)
      else
        ierr = pio_def_var(File, trim(name), dtype, vardesc)
      end if
      call cam_pio_handle_error(ierr, trim(subname)//': Unable to define variable '//trim(name))
    end if

    ! Back to whatever error handling was running before this routine
    call pio_seterrorhandling(File, err_handling)

  end subroutine cam_pio_def_var_md

  subroutine cam_pio_get_var_2d_r8(varname, File, field, start, kount, found)
    use cam_abortutils, only: endrun
    use pio,            only: file_desc_t, var_desc_t, pio_get_var, PIO_MAX_NAME
    use pio,            only: pio_inq_dimname

    ! Dummy arguments
    character(len=*),  intent(in)    :: varname
    type(file_desc_t), intent(inout) :: File      ! PIO file Handle
    real(r8),          intent(inout) :: field(:,:)
    integer, optional, intent(in)    :: start(2)
    integer, optional, intent(in)    :: kount(2)
    logical, optional, intent(out)   :: found
    
    ! Local variables
    character(len=*), parameter      :: subname = 'cam_pio_get_var_2d_r8'
    character(len=PIO_MAX_NAME)      :: tmpname
    type(var_desc_t)                 :: varid   ! Var descriptor
    integer                          :: ierr
    integer                          :: strt(3)
    integer                          :: cnt(3)
    integer                          :: ndims
    integer                          :: dimids(3)
    logical                          :: exists

    if ( (present(start) .and. (.not. present(kount))) .or.                   &
         (present(kount) .and. (.not. present(start)))) then
      call endrun(trim(subname)//': start and kount must both be present')
    end if
      
    call cam_pio_find_var(File, trim(varname), varid, exists)
    if (present(found)) then
      found = exists
    else if (.not. exists) then
      call endrun(trim(subname)//': '//trim(varname)//' not found')
    end if
    if (exists) then
      call cam_pio_var_info(File, varid, ndims, dimids, cnt, varname=varname)

      if (present(start)) then
        ! start and kount override other options and are not error checked
        strt(1:2) = start(1:2)
        strt(3) = 1
        cnt(1:2) = kount(1:2)
        cnt(3) = 1
      else
        strt = 1     ! cnt set by cam_pio_var_info
        exists = use_scam_limits(File, strt, cnt)
      end if
      if (ndims == 3) then
        ierr = pio_inq_dimname(File, dimids(3), tmpname)
        if (trim(tmpname) /= 'time') then
          call endrun(trim(subname)//': dimension mismatch for '//trim(varname))
        else
          ierr = pio_get_var(File, varid, strt, cnt, field)
        end if
      else if (ndims == 2) then
        ierr = pio_get_var(File, varid, strt, cnt, field)
      else if (ndims == 1) then
        ierr = pio_get_var(File, varid, strt(1:1), cnt(1:1), field(:,1))
      else
        call endrun(trim(subname)//': Incorrect variable rank')
      end if
    end if

  end subroutine cam_pio_get_var_2d_r8

  subroutine cam_pio_get_var_2d_r8_perm(varname, File, arraydims, field,      &
       start, kount, found)
    use cam_abortutils, only: endrun
    use pio,            only: file_desc_t, var_desc_t, pio_get_var, PIO_MAX_NAME

    ! Dummy arguments
    character(len=*),  intent(in)    :: varname
    type(file_desc_t), intent(inout) :: File      ! PIO file Handle
    character(len=*),  intent(in)    :: arraydims(2)
    real(r8),          intent(inout) :: field(:,:)
    integer, optional, intent(in)    :: start(2)
    integer, optional, intent(in)    :: kount(2)
    logical, optional, intent(out)   :: found
    
    ! Local variables
    character(len=*), parameter      :: subname = 'cam_pio_get_var_2d_r8_perm'
    type(var_desc_t)                 :: varid   ! Var descriptor
    integer                          :: ierr
    integer                          :: i, j, ind(2)
    integer                          :: strt(3)
    integer                          :: cnt(3)
    integer                          :: ndims
    integer                          :: dimids(3)
    integer                          :: perm(2)
    logical                          :: isperm
    logical                          :: exists
    real(r8), allocatable            :: tmp_fld(:,:)
    character(len=PIO_MAX_NAME)      :: filedims(3)

    if ( (present(start) .and. (.not. present(kount))) .or.                   &
         (present(kount) .and. (.not. present(start)))) then
      call endrun(trim(subname)//': start and kount must both be present')
    end if

    call cam_pio_find_var(File, trim(varname), varid, exists)

    if (present(found)) then
      found = exists
    else if (.not. exists) then
      call endrun(trim(subname)//': '//trim(varname)//' not found')
    end if
    if (exists) then
      call cam_pio_var_info(File, varid, ndims, dimids, cnt,                  &
           dimnames=filedims, varname=varname)

      if (present(start)) then
        ! start and kount override other options and are not error checked
        strt(1:2) = start
        strt(3) = 1
        cnt(1:2) = kount
      else
        strt = 1   ! cnt set by cam_pio_var_info
        exists = use_scam_limits(File, strt, cnt)
      end if
      if ( ((ndims == 2) .and. (trim(filedims(2)) /= 'time')) .or.            &
           ((ndims == 3) .and. (trim(filedims(3)) == 'time'))) then
        call calc_permutation(filedims(1:2), arraydims, perm, isperm)
        if (isperm) then
          allocate(tmp_fld(cnt(1), cnt(2)))
          ierr = pio_get_var(File, varid, strt(1:ndims), cnt(1:ndims), tmp_fld)
          do j = 1, cnt(2)
            ind(2) = j
            do i = 1, cnt(1)
              ind(1) = i
              field(ind(perm(1)), ind(perm(2))) = tmp_fld(i, j)
            end do
          end do
        else
          ierr = pio_get_var(File, varid, strt(1:ndims), cnt(1:ndims), field)
        end if
      else
        call endrun(trim(subname)//': Incorrect variable rank')
      end if
    end if

  end subroutine cam_pio_get_var_2d_r8_perm

  subroutine cam_pio_get_var_3d_r8(varname, File, field, start, kount, found)
    use cam_abortutils, only: endrun
    use pio,            only: file_desc_t, var_desc_t, pio_get_var, PIO_MAX_NAME
    use pio,            only: pio_inq_dimname

    ! Dummy arguments
    character(len=*),  intent(in)    :: varname
    type(file_desc_t), intent(inout) :: File          ! PIO file Handle
    real(r8),          intent(inout) :: field(:,:,:)
    integer, optional, intent(in)    :: start(3)
    integer, optional, intent(in)    :: kount(3)
    logical, optional, intent(out)   :: found
    
    ! Local variables
    character(len=*), parameter      :: subname = 'cam_pio_get_var_3d_r8'
    character(len=PIO_MAX_NAME)      :: tmpname
    type(var_desc_t)                 :: varid   ! Var descriptor
    integer                          :: ierr
    integer                          :: strt(4)
    integer                          :: cnt(4)
    integer                          :: ndims
    integer                          :: dimids(4)
    logical                          :: exists

    if ( (present(start) .and. (.not. present(kount))) .or.                   &
         (present(kount) .and. (.not. present(start)))) then
      call endrun(trim(subname)//': start and kount must both be present')
    end if

   call cam_pio_find_var(File, trim(varname), varid, exists)

    if (present(found)) then
      found = exists
    else if (.not. exists) then
      call endrun(trim(subname)//': '//trim(varname)//' not found')
    end if
    if (exists) then
      call cam_pio_var_info(File, varid, ndims, dimids, cnt, varname=varname)

      if (present(start)) then
        ! start and kount override other options and are not error checked
        strt(1:3) = start(1:3)
        strt(4) = 1
        cnt(1:3) = kount(1:3)
        cnt(4) = 1
      else
        strt = 1    ! cnt set by cam_pio_var_info
        exists = use_scam_limits(File, strt, cnt)
      end if

      if (ndims == 4) then
        ierr = pio_inq_dimname(File, dimids(4), tmpname)
        if (trim(tmpname) /= 'time') then
          call endrun(trim(subname)//': dimension mismatch for '//trim(varname))
        else
          ierr = pio_get_var(File, varid, strt, cnt, field)
        end if
      else if (ndims == 3) then
        ierr = pio_get_var(File, varid, strt, cnt, field)
      else if (ndims == 2) then
        ierr = pio_get_var(File, varid, strt(1:ndims), cnt(1:ndims), field(:,:,1))
      else
        call endrun(trim(subname)//': Incorrect variable rank')
      end if
    end if

  end subroutine cam_pio_get_var_3d_r8

  subroutine cam_pio_get_var_3d_r8_perm(varname, File, arraydims, field,      &
       start, kount, found)
    use cam_abortutils, only: endrun
    use pio,            only: file_desc_t, var_desc_t, pio_get_var, PIO_MAX_NAME

    ! Dummy arguments
    character(len=*),  intent(in)    :: varname
    type(file_desc_t), intent(inout) :: File      ! PIO file Handle
    character(len=*),  intent(in)    :: arraydims(3)
    real(r8),          intent(inout) :: field(:,:,:)
    integer, optional, intent(in)    :: start(3)
    integer, optional, intent(in)    :: kount(3)
    logical, optional, intent(out)   :: found
    
    ! Local variables
    character(len=*), parameter      :: subname = 'cam_pio_get_var_3d_r8_perm'
    type(var_desc_t)                 :: varid   ! Var descriptor
    integer                          :: ierr
    integer                          :: i, j, k, ind(3)
    integer                          :: strt(4)
    integer                          :: cnt(4)
    integer                          :: ndims
    integer                          :: dimids(4)
    integer                          :: perm(3)
    logical                          :: exists
    logical                          :: isperm
    real(r8), allocatable            :: tmp_fld(:,:,:)
    character(len=PIO_MAX_NAME)      :: filedims(4)

    if ( (present(start) .and. (.not. present(kount))) .or.                   &
         (present(kount) .and. (.not. present(start)))) then
      call endrun(trim(subname)//': start and kount must both be present')
    end if

    call cam_pio_find_var(File, trim(varname), varid, exists)

    if (present(found)) then
      found = exists
    else if (.not. exists) then
      call endrun(trim(subname)//': '//trim(varname)//' not found')
    end if
    if (exists) then
      call cam_pio_var_info(File, varid, ndims, dimids, cnt,                  &
           dimnames=filedims, varname=varname)

      if (present(start)) then
        ! start and kount override other options and are not error checked
        strt(1:3) = start
        strt(4) = 1
        cnt(1:3) = kount
      else
        strt = 1   ! cnt set by cam_pio_var_info
        exists = use_scam_limits(File, strt, cnt)
      end if

      if ( ((ndims == 3) .and. (trim(filedims(3)) /= 'time')) .or.            &
           ((ndims == 4) .and. (trim(filedims(4)) == 'time'))) then
        call calc_permutation(filedims(1:3), arraydims, perm, isperm)
        if (isperm) then
          allocate(tmp_fld(cnt(1), cnt(2), cnt(3)))
          ierr = pio_get_var(File, varid, strt(1:3), cnt(1:3), tmp_fld)
          do k = 1, cnt(3)
            ind(3) = k
            do j = 1, cnt(2)
              ind(2) = j
              do i = 1, cnt(1)
                ind(1) = i
                field(ind(perm(1)), ind(perm(2)), ind(perm(3))) = tmp_fld(i, j, k)
              end do
            end do
          end do
        else
          ierr = pio_get_var(File, varid, strt(1:3), cnt(1:3), field)
        end if
      else
        call endrun(trim(subname)//': Incorrect variable rank')
      end if
    end if

  end subroutine cam_pio_get_var_3d_r8_perm

  ! clean_iodesc_list: Deallocate all entries in the iodesc list
  subroutine clean_iodesc_list()
    type(iodesc_list), pointer :: this, prev

    if(associated(iodesc_list_top%iodesc)) then
      ! iodesc_list_top is not allocated so leave it (just empty)
      this => iodesc_list_top
      iodesc_list_top%tag = ''
      call pio_freedecomp(pio_subsystem, this%iodesc)
      deallocate(this%iodesc)
      nullify(this%iodesc)
      this => this%next
      nullify(iodesc_list_top%next)
       
      ! All the other list items were allocated, blow them away
      do while(associated(this))
        call pio_freedecomp(pio_subsystem, this%iodesc)
        deallocate(this%iodesc)
        prev => this
        this => this%next
        deallocate(prev)
      end do
    end if
  end subroutine clean_iodesc_list

  subroutine cam_pio_createfile(file, fname, mode_in)
    use pio, only : pio_createfile, file_desc_t, pio_noerr, pio_clobber,      &
         pio_64bit_offset, pio_iotask_rank
    use cam_abortutils, only : endrun

    ! Dummy arguments
    type(file_desc_t),          intent(inout) :: file
    character(len=*),           intent(in)    :: fname
    integer,          optional, intent(in)    :: mode_in

    ! Local variables
    integer                                   :: ierr
    integer                                   :: mode
    
    mode = ior(PIO_CLOBBER, PIO_64BIT_OFFSET)
    if (present(mode_in)) then
      mode = ior(mode, mode_in)
    end if

    ierr = pio_createfile(pio_subsystem, file, pio_iotype, fname, mode)

    if(ierr /= PIO_NOERR) then
       call endrun('Failed to open file,'//trim(fname)//', to write')
    else if(pio_iotask_rank(pio_subsystem) == 0) then
       write(iulog, *) 'Opened file ', trim(fname),  ' to write', file%fh
    end if

  end subroutine cam_pio_createfile

  subroutine cam_pio_openfile(file, fname, mode)
    use pio,            only: pio_openfile, file_desc_t, pio_noerr, pio_iotask_rank
    use cam_abortutils, only: endrun

    type(file_desc_t), intent(inout), target :: file
    character(len=*), intent(in) :: fname
    integer, intent(in) :: mode

    integer :: ierr

    ierr = pio_openfile(pio_subsystem, file, pio_iotype, fname, mode)

    if(ierr/= PIO_NOERR) then
       call endrun('Failed to open restart file to read')
    else if(pio_iotask_rank(pio_subsystem) == 0) then
       write(iulog,*) 'Opened existing file ', trim(fname), file%fh
    end if

  end subroutine cam_pio_openfile

  subroutine cam_pio_closefile(file)

    use pio, only : pio_closefile, file_desc_t

    type(file_desc_t), intent(inout), target :: file

    call pio_closefile(file)

  end subroutine cam_pio_closefile

  logical function cam_pio_fileexists(fname)
    use pio,            only: pio_openfile, file_desc_t, pio_noerr, PIO_NOWRITE
    use pio,            only: pio_seterrorhandling, PIO_BCAST_ERROR
    use pio,            only : pio_closefile

    character(len=*), intent(in) :: fname

    type(file_desc_t)            :: file
    integer                      :: ierr
    integer                      :: err_handling

    ! We will handle errors for this routine
    !!XXgoldyXX: This hack should be replaced with the PIO interface
    !err_handling = pio_subsystem%error_handling !! Hack
    call pio_seterrorhandling(pio_subsystem, PIO_BCAST_ERROR, err_handling)

    ierr = pio_openfile(pio_subsystem, file, pio_iotype, fname, PIO_NOWRITE)
    cam_pio_fileexists = (ierr == PIO_NOERR)
    if (cam_pio_fileexists) then
      call pio_closefile(file)
    end if

    ! Back to whatever error handling was running before this routine
    call pio_seterrorhandling(File, err_handling)

  end function cam_pio_fileexists

end module cam_pio_utils
