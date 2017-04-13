module esp_utils

! !USES:

  use shr_kind_mod, only: r8=>shr_kind_r8, CL=>SHR_KIND_CL
  use shr_sys_mod,  only: shr_sys_abort, shr_sys_flush

  implicit none
  private

  public :: esp_pio_modify_variable

CONTAINS

  subroutine esp_modify_array(array, pertlim, absolute)
    real(r8), intent(inout) :: array(:)
    real(r8), intent(in)    :: pertlim
    logical,  intent(in)    :: absolute

    integer                 :: ind
    integer                 :: rndm_seed_sz
    integer, allocatable    :: rndm_seed(:)
    real(r8)                :: pertval

    call random_seed(size=rndm_seed_sz)
    allocate(rndm_seed(rndm_seed_sz))
    rndm_seed = 8.53973422267357_r8
    call random_seed(put=rndm_seed)
    do ind = 1, size(array)
      call random_number(pertval)
      pertval = 2._r8 * pertlim * (0.5_r8 - pertval)
      if (absolute) then
        array(ind) = pertval
      else
        array(ind) = array(ind) * (1.0_r8 + pertval)
      end if
    end do

  end subroutine esp_modify_array

  subroutine esp_pio_handle_error(ierr, errorstr)
    use pio,          only: pio_noerr

    ! Dummy arguments
    integer,          intent(in)  :: ierr
    character(len=*), intent(in)  :: errorstr

    ! Local variables
    character(len=256) :: errormsg

    if (ierr /= PIO_NOERR) then
      write(errormsg, '(a,i6,2a)') '(PIO:', ierr, ') ', trim(errorstr)
      call shr_sys_abort(errormsg)
    end if
    
  end subroutine esp_pio_handle_error

  subroutine esp_pio_openfile(file, fname, piosys, iotype, mode, iulog)
    use pio, only: pio_openfile, file_desc_t, pio_iotask_rank, iosystem_desc_t
    use pio, only: PIO_NOWRITE, PIO_WRITE, PIO_NOERR

    type(file_desc_t),     target,   intent(inout) :: file
    character(len=*),                intent(in)    :: fname
    type(iosystem_desc_t), pointer                 :: piosys
    integer,                         intent(in)    :: iotype
    integer,                         intent(in)    :: mode
    integer,               optional, intent(in)    :: iulog

    integer                                        :: ierr
    character(len=CL)                              :: errmsg
    character(len=*),      parameter               :: subname = 'ESP_PIO_OPENFILE: '

    ierr = pio_openfile(piosys, file, iotype, fname, mode)

    if(ierr/= PIO_NOERR) then
      if (mode == PIO_nowrite) then
        write(errmsg, '(3a,i0)') 'Failed to open ',trim(fname),' to read, error = ',ierr
        call shr_sys_abort(subname//trim(errmsg))
      else if (mode == PIO_write) then
        write(errmsg, '(3a,i0)') 'Failed to open ',trim(fname),' to write, error = ',ierr
        call shr_sys_abort(subname//trim(errmsg))
      else
        write(errmsg, '(3a,i0,a,i0)') 'Failed to open ',trim(fname),' with mode = ',mode,', error = ',ierr
        call shr_sys_abort(subname//trim(errmsg))
      end if
    else if((pio_iotask_rank(piosys) == 0) .and. present(iulog)) then
      write(iulog,'(2a)') 'Opened existing file, ', trim(fname)
      call shr_sys_flush(iulog)
    end if

  end subroutine esp_pio_openfile

  subroutine esp_pio_closefile(file)

    use pio, only : pio_closefile, file_desc_t

    type(file_desc_t), intent(inout), target :: file

    call pio_closefile(file)

  end subroutine esp_pio_closefile

  logical function esp_pio_fileexists(fname, piosys, iotype)
    use pio, only: pio_openfile, file_desc_t, iosystem_desc_t
    use pio, only: pio_seterrorhandling, PIO_BCAST_ERROR
    use pio, only: pio_closefile, PIO_NOERR, PIO_NOWRITE

    character(len=*), intent(in)   :: fname
    type(iosystem_desc_t), pointer :: piosys
    integer,          intent(in)   :: iotype

    type(file_desc_t)              :: file
    integer                        :: ierr
    integer                        :: err_handling

    ! We will handle errors for this routine

    call pio_seterrorhandling(piosys, PIO_BCAST_ERROR, err_handling)

    ierr = pio_openfile(piosys, file, iotype, fname, PIO_NOWRITE)
    esp_pio_fileexists = (ierr == PIO_NOERR)
    if (esp_pio_fileexists) then
      call pio_closefile(file)
    end if

    ! Back to whatever error handling was running before this routine
    call pio_seterrorhandling(File, err_handling)

  end function esp_pio_fileexists

  !-----------------------------------------------------------------------
  !
  ! esp_pio_var_info: Retrieve variable properties
  !
  !-----------------------------------------------------------------------
  subroutine esp_pio_var_info(ncid, varid, ndims, dimids, dimlens, dimnames, varname, unlimDimID)
    use pio,        only: file_desc_t, var_desc_t
    use pio,        only: PIO_inq_varndims, PIO_inq_vardimid, PIO_inq_dimlen
    use pio,        only: PIO_inquire, PIO_inq_dimname
    use pio,        only: PIO_seterrorhandling, PIO_BCAST_ERROR


    ! Dummy arguments
    type(file_desc_t),           intent(inout) :: ncid
    type(var_desc_t),            intent(in)    :: varid
    integer,                     intent(out)   :: ndims
    integer,                     intent(out)   :: dimids(:)
    integer,                     intent(out)   :: dimlens(:)
    character(len=*),  optional, intent(out)   :: dimnames(:)
    integer,           optional, intent(out)   :: unlimDimID
    character(len=*),  optional, intent(in)    :: varname

    ! Local variables
    integer                                    :: ret     ! PIO return value
    integer                                    :: i
    integer                                    :: err_handling
    character(len=128)                         :: errsuff
    !-----------------------------------------------------------------------
    ! We will handle errors for this routine

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
    call esp_pio_handle_error(ret, 'ESP_PIO_VAR_INFO: Error with num dimensions')
    if (size(dimids) < ndims) then
      call shr_sys_abort('ESP_PIO_VAR_INFO: dimids too small'//trim(errsuff))
    end if
    ret = PIO_inq_vardimid(ncid, varid, dimids(1:ndims))
    call esp_pio_handle_error(ret, 'ESP_PIO_VAR_INFO: Error with inq dim ids'//trim(errsuff))
    if (size(dimlens) < ndims) then
      call shr_sys_abort('ESP_PIO_VAR_INFO: dimlens too small'//trim(errsuff))
    end if
    do i = 1, ndims
      ret = PIO_inq_dimlen(ncid, dimids(i), dimlens(i))
      call esp_pio_handle_error(ret, 'ESP_PIO_VAR_INFO: Error with inq dimlens')
      if (present(dimnames)) then
        ret = PIO_inq_dimname(ncid, dimids(i), dimnames(i))
        call esp_pio_handle_error(ret, 'ESP_PIO_VAR_INFO: Error with inq dimnames')
      end if
    end do
    if (present(unlimDimID)) then
      ret = PIO_inquire(ncid, unlimitedDimID=unlimDimID)
      call esp_pio_handle_error(ret, 'ESP_PIO_VAR_INFO: Error with inquire')
    end if
    call PIO_seterrorhandling(ncid, err_handling)

  end subroutine esp_pio_var_info

  subroutine esp_pio_find_var(ncid, varname, varid, found)
    use pio, only: file_desc_t, var_desc_t
    use pio, only: pio_inq_varid, pio_noerr
    use pio, only: PIO_seterrorhandling, PIO_BCAST_ERROR

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

    call PIO_seterrorhandling(ncid, PIO_BCAST_ERROR, err_handling)
    ret = PIO_inq_varid(ncid, trim(varname), varid)
    found = (ret == PIO_NOERR)
    call PIO_seterrorhandling(ncid, err_handling)

  end subroutine esp_pio_find_var

  subroutine esp_pio_newdecomp(iodesc, piosys, iotype, dims, dof, dtype)
    use pio, only: pio_initdecomp, pio_offset_kind, pio_iotype_pnetcdf
    use pio, only: io_desc_t, iosystem_desc_t, PIO_REARR_SUBSET, PIO_REARR_BOX

    type(io_desc_t)                           :: iodesc
    type(iosystem_desc_t), pointer            :: piosys
    integer,                       intent(in) :: iotype
    integer,                       intent(in) :: dims(:)
    integer(kind=PIO_OFFSET_KIND), intent(in) :: dof(:)
    integer,                       intent(in) :: dtype

    integer                                   :: pio_rearranger

    if(iotype == pio_iotype_pnetcdf) then
       pio_rearranger = PIO_REARR_SUBSET
    else
       pio_rearranger = PIO_REARR_BOX
    endif

    call pio_initdecomp(piosys, dtype, dims, dof, iodesc, rearr=pio_rearranger)

  end subroutine esp_pio_newdecomp

  subroutine esp_pio_modify_variable(id, comm, filename, varname, found)
    use mpi,         only: MPI_LOGICAL, MPI_LAND
    use shr_mpi_mod, only: shr_mpi_commsize, shr_mpi_commrank
    use shr_pio_mod, only: shr_pio_getiosys, shr_pio_getiotype
    use pio,         only: PIO_write, file_desc_t, pio_offset_kind
    use pio,         only: io_desc_t, var_desc_t, pio_freedecomp, PIO_DOUBLE
    use pio,         only: pio_read_darray, pio_write_darray, iosystem_desc_t

    integer,          intent(in)          :: id
    integer,          intent(in)          :: comm
    character(len=*), intent(in)          :: filename
    character(len=*), intent(in)          :: varname
    logical,          intent(out)         :: found

    type(file_desc_t)                     :: file
    integer                               :: ierr
    integer                               :: mode
    integer                               :: pio_iotype
    type(var_desc_t)                      :: varid
    integer                               :: ndims
    integer                               :: dimids(7)
    integer                               :: dimlens(7)
    integer                               :: totalsize, mysize
    integer                               :: i
    integer                               :: uid
    integer                               :: npes, iam
    logical                               :: equiv, all_equiv
    integer                               :: offset
    real(r8),                 allocatable :: varr(:)
    integer(pio_offset_kind), allocatable :: ldof(:)
    type(io_desc_t)                       :: iodesc
    type(iosystem_desc_t),    pointer     :: pio_subsystem
    character(len=*),         parameter   :: subname = 'ESP_PIO_MODIFY_VARIABLE: '

    pio_subsystem => shr_pio_getiosys(id)
    pio_iotype =  shr_pio_getiotype(id)
    mode = PIO_WRITE

    call esp_pio_openfile(file, filename, pio_subsystem, pio_iotype, mode)
    ! Find the variable
    call esp_pio_find_var(file, varname, varid, found)
    if (found) then
      ! Check dimensions
      call esp_pio_var_info(file, varid, ndims, dimids, dimlens,              &
           varname=varname, unlimDimID=uid)
      ! Skip the unlimited dimension if it is in varname
      ierr = 1
      do i = 1, ndims
        if (i > ierr) then
          dimids(ierr) = dimids(i)
          dimlens(ierr) = dimlens(i)
        end if
        if (dimids(i) /= uid) then
          ierr = ierr + 1
        end if
      end do
      ndims = ndims - COUNT(dimids(1:ndims) == uid)
      ! Calculate global and local array sizes
      totalsize = PRODUCT(dimlens(1:ndims))
      call shr_mpi_commsize(comm, npes)
      call shr_mpi_commrank(comm, iam)
      mysize = totalsize / npes
      if (iam < MOD(totalsize, npes)) then
        mysize = mysize + 1
      end if
      allocate(varr(mysize))
      allocate(ldof(mysize))
      offset = (iam * (totalsize / npes)) + MIN(iam, MOD(totalsize, npes))
      do i = 1, mysize
        ldof(i) = i + offset
      end do
      call esp_pio_newdecomp(iodesc, pio_subsystem, pio_iotype,               &
           dimlens(1:ndims), ldof, PIO_DOUBLE)
      call pio_read_darray(file, varid, iodesc, varr, ierr)
      call esp_pio_handle_error(ierr, subname//'Error reading variable '//trim(varname))
      ! See if we have a constant zero value
      equiv = ALL(varr == 0.0_r8)
      call mpi_allreduce(equiv, all_equiv, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
      ! Modify and write back to file
      call esp_modify_array(varr, 1.0e-12_r8, all_equiv)
      call pio_write_darray(file, varid, iodesc, varr, ierr)
      call esp_pio_handle_error(ierr, subname//'Error writing variable '//trim(varname))
      ! Cleanup
      call pio_freedecomp(pio_subsystem, iodesc)
    end if

    call esp_pio_closefile(file)
    if (allocated(varr)) then
      deallocate(varr)
    end if
    if (allocated(ldof)) then
      deallocate(ldof)
    end if

  end subroutine esp_pio_modify_variable

end module esp_utils
