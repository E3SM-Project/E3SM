!>
!! @file
!! @brief Module containing netcdf-specific PIO unit tests
!<

module ncdf_tests

  use pio
  use global_vars

  Implicit None
  private
  save

  public :: test_redef
  public :: test_enddef

Contains

  Subroutine test_redef(test_id, err_msg)
    ! test_redef():
    ! * Open file, enter define mode, add dimension / variable / attribute
    !   * Try to run PIO_redef from define mode, check for error
    ! * Leave define mode, close file
    !   * Try to run PIO_redef with closed file
    ! Routines used in test: PIO_initdecomp, PIO_openfile, PIO_redef, PIO_def_dim,
    !                        PIO_def_var, PIO_put_att, PIO_enddef,
    !                        PIO_write_darray, PIO_closefile, PIO_freedecomp

    ! Input / Output Vars
    integer,                intent(in)  :: test_id
    character(len=str_len), intent(out) :: err_msg

    ! Local Vars
    character(len=str_len) :: filename
    integer                :: iotype, ret_val

    ! Data used to test writing
    integer,          dimension(2) :: data_to_write, compdof
    integer,          dimension(1) :: dims
    type(io_desc_t)                :: iodesc_nCells
    integer                        :: pio_dim
    type(var_desc_t)               :: pio_var

    err_msg = "no_error"

    dims(1) = 2*ntasks
    compdof = 2*my_rank+(/1,2/)  ! Where in the global array each task writes
    data_to_write = 1+my_rank

    call PIO_initdecomp(pio_iosystem, PIO_int, dims, compdof, iodesc_nCells)

    filename = fnames(test_id)
    iotype   = iotypes(test_id)

    ! Open existing file, write data to it
    ret_val = PIO_openfile(pio_iosystem, pio_file, iotype, filename, PIO_write)
    if (ret_val .ne. PIO_NOERR) then
       ! Error in PIO_openfile
       err_msg = "Could not open " // trim(filename) // " in write mode"
       return
    end if

    ! Enter define mode
    ret_val = PIO_redef(pio_file)
    if (ret_val .ne. PIO_NOERR) then
       ! Error in PIO_redef
       err_msg = "Could not enter define mode"
       call PIO_closefile(pio_file)
       return
    end if

    ! Define a new dimension M (already has 'N' from previous tests)
    ret_val = PIO_def_dim(pio_file, 'M', int(2*ntasks,pio_offset_kind), pio_dim)
    if (ret_val .ne. PIO_NOERR) then
       err_msg = "Could not define dimension M"
       call PIO_closefile(pio_file)
       return
    end if

    ! Define a new variable foo2 (already has 'foo' from previous tests)
    ret_val = PIO_def_var(pio_file, 'foo2', PIO_int, &
         (/pio_dim/), pio_var)
    if (ret_val .ne. PIO_NOERR) then
       ! Error in PIO_def_var
       err_msg = "Could not define variable foo2"
       call PIO_closefile(pio_file)
       return
    end if

    ret_val = PIO_put_att(pio_file, pio_var, "max_val", ntasks)
    if (ret_val .ne. PIO_NOERR) then
       ! Error in PIO_put_att
       err_msg = "Could not define max_val attribute for foo2"
       call PIO_closefile(pio_file)
       return
    end if

    ret_val = PIO_put_att(pio_file, PIO_global, "created_by", "PIO unit tests")
    if (ret_val .ne. PIO_NOERR) then
       ! Error in PIO_put_att
       err_msg = "Could not define global attribute"
       call PIO_closefile(pio_file)
       return
    end if

    ! Try to enter define mode again
    if(master_task) write(*,"(6x,A,1x)") "trying to enter define mode in define mode, error expected ... "
    call mpi_barrier(MPI_COMM_WORLD,ret_val)

    ret_val = PIO_redef(pio_file)
    if (ret_val .eq. PIO_NOERR) then
       ! Error in PIO_redef
       err_msg = "Entered define mode from define mode"
       call PIO_closefile(pio_file)
       return
    end if

    ! Leave define mode
    ret_val = PIO_enddef(pio_file)
    if (ret_val .ne. PIO_NOERR) then
       ! Error in PIO_enddef
       print *,__FILE__,__LINE__,ret_val
       err_msg = "Could not end define mode"
       return
    end if

    ! Write foo2
    call PIO_write_darray(pio_file, pio_var, iodesc_nCells, data_to_write, ret_val)
    if (ret_val .ne. PIO_NOERR) then
       ! Error in PIO_write_darray
       err_msg = "Could not write data"
       return
    end if

    ! Close file
    call PIO_closefile(pio_file)

    ! Try to enter define mode again
    if(master_task) write(*,"(6x,A,1x)") "trying to enter define mode in closed file, error expected ... "
    call mpi_barrier(MPI_COMM_WORLD,ret_val)
    ret_val = PIO_redef(pio_file)
    if (ret_val .eq. PIO_NOERR) then
       ! Error in PIO_redef
       err_msg = "Entered define mode from a closed file"
       return
    end if

    ! Free decomp
    call PIO_freedecomp(pio_iosystem, iodesc_nCells)

  End Subroutine test_redef

  Subroutine test_enddef(test_id, err_msg)
    ! test_enddef():
    ! * Open file with PIO_nowrite, try to enter define mode, check for error
    ! * Open file with PIO_write, enter define mode, leave define mode
    ! * Try calling PIO_enddef from data mode, check for error
    ! * Close file
    !   * Try to run PIO_enddef with closed file
    ! Routines used in test: PIO_openfile, PIO_redef, PIO_enddef, PIO_closefile

    ! Input / Output Vars
    integer,                intent(in)  :: test_id
    character(len=str_len), intent(out) :: err_msg

    ! Local Vars
    character(len=str_len) :: filename
    integer                :: iotype, ret_val

    err_msg = "no_error"
    filename = fnames(test_id)
    iotype   = iotypes(test_id)

    ! Open existing file (read-only)
    ret_val = PIO_openfile(pio_iosystem, pio_file, iotype, filename, PIO_nowrite)
    if (ret_val .ne. PIO_NOERR) then
       ! Error in PIO_openfile
       err_msg = "Could not open " // trim(filename) // " in write mode"
       return
    end if

    ! Enter define mode
    if(master_task) write(*,"(6x,A,1x)") "trying to enter define mode in read only file, error expected ... "
    call mpi_barrier(MPI_COMM_WORLD,ret_val)
    ret_val = PIO_redef(pio_file)
    if (ret_val .eq. PIO_NOERR) then
       ! Error in PIO_redef
       err_msg = "Entered define mode in read-only file"
       call PIO_closefile(pio_file)
       return
    end if

    ! Close file
    call PIO_closefile(pio_file)

    ! Open existing file
    ret_val = PIO_openfile(pio_iosystem, pio_file, iotype, filename, PIO_write)
    if (ret_val .ne. PIO_NOERR) then
       ! Error in PIO_openfile
       err_msg = "Could not open " // trim(filename) // " in write mode"
       return
    end if

    ! Enter define mode
    ret_val = PIO_redef(pio_file)
    if (ret_val .ne. PIO_NOERR) then
       ! Error in PIO_redef
       err_msg = "Could not enter define mode"
       call PIO_closefile(pio_file)
       return
    end if

    ! End define mode
    ret_val = PIO_enddef(pio_file)
    if (ret_val .ne. PIO_NOERR) then
       ! Error in PIO_enddef
       print *,__FILE__,__LINE__,ret_val
       err_msg = "Could not end define mode"
       call PIO_closefile(pio_file)
       return
    end if

    ! Try to end define mode from data mode
    if(master_task) write(*,"(6x,A,1x)") "trying to end define mode in data mode, error expected ... "
    call mpi_barrier(MPI_COMM_WORLD,ret_val)
    ret_val = PIO_enddef(pio_file)
    if (ret_val .eq. PIO_NOERR) then
       ! Error in PIO_enddef
       err_msg = "Ended define mode while in data mode"
       call PIO_closefile(pio_file)
       return
    end if

    call mpi_barrier(MPI_COMM_WORLD,ret_val)
    ! Close file
    call PIO_closefile(pio_file)

    ! Try to end define mode in un-opened file
    if(master_task) write(*,"(6x,A,1x)") "trying to end define mode in closed file, error expected ... "
    call mpi_barrier(MPI_COMM_WORLD,ret_val)
    ret_val = PIO_enddef(pio_file)
    if (ret_val .eq. PIO_NOERR) then
       ! Error in PIO_enddef
       err_msg = "Ended define mode in a file that was already closed"
       return
    end if
    call mpi_barrier(MPI_COMM_WORLD,ret_val)

  End Subroutine test_enddef

end module ncdf_tests
