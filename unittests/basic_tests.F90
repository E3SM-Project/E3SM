!>
!! @file
!! @brief Module containing basic unit tests that are run for both
!!        binary and netcdf file types.
!<

module basic_tests

  use pio 
  use global_vars

  Implicit None
  private
  save

  public :: test_create
  public :: test_open

  Contains

    Subroutine test_create(test_id, err_msg)
    ! test_create():
    ! * Create an empty file, close it, test that it can be opened
    ! * For netcdf / pnetcdf:
    !   - Create same file with PIO_CLOBBER mode
    !   - Try to create same file with PIO_NOCLOBBER mode, check for error
    ! Routines used in test: PIO_createfile, PIO_openfile, PIO_closefile
    ! Also uses PIO_enddef for [p]netcdf tests

      ! Input / Output Vars
      integer,                intent(in)  :: test_id
      character(len=str_len), intent(out) :: err_msg

      ! Local Vars
      character(len=str_len) :: filename
      integer                :: iotype, ret_val, pio_dim

      err_msg = "no_error"

      filename = fnames(test_id)
      iotype   = iotypes(test_id)

      ! Delete file before initial create
!      if (master_task) call system("rm -f " // trim(filename))

      call PIO_deletefile(pio_iosystem, filename)

      ! Original file creation
      ret_val = PIO_createfile(pio_iosystem, pio_file, iotype, filename)
      if (ret_val.ne.0) then
        ! Error in PIO_createfile
        err_msg = "Could not create " // trim(filename)
        return
      end if
      

      ! netcdf files need to end define mode before closing
      if (is_netcdf(iotype)) then
        ret_val = PIO_enddef(pio_file)
        if (ret_val.ne.0) then
          ! Error in PIO_enddef
          err_msg = "Could not end define mode"
          call PIO_closefile(pio_file)
          return
        end if
      end if
      call PIO_closefile(pio_file)


      ! Test opening of file
      ret_val = PIO_openfile(pio_iosystem, pio_file, iotype, filename, PIO_nowrite)
      if (ret_val.ne.0) then
        ! Error in PIO_openfile
        err_msg = "Could not open " // trim(filename)
        return
      end if

      ! Close file
      call PIO_closefile(pio_file)

      ! Recreate file with CLOBBER (netcdf / pnetcdf only)
      if (is_netcdf(iotype)) then
        ret_val = PIO_createfile(pio_iosystem, pio_file, iotype, filename, PIO_CLOBBER)
        if (ret_val.ne.0) then
          ! Error in PIO_createfile
          err_msg = "Could not clobber " // trim(filename)
          return
        end if

        ! Leave define mode
        ret_val = PIO_enddef(pio_file)
        if (ret_val.ne.0) then
          ! Error in PIO_enddef
          err_msg = "Could not end define mode in clobbered file"
          call PIO_closefile(pio_file)
          return
        end if

        ! Close file
        call PIO_closefile(pio_file)
      end if

      ! Recreate file with NOCLOBBER
      if (is_netcdf(iotype)) then
        ret_val = PIO_createfile(pio_iosystem, pio_file, iotype, filename, PIO_NOCLOBBER)
        if (ret_val.eq.0) then
          ! Error in PIO_createfile
          err_msg = "Was able to clobber file despite PIO_NOCLOBBER"
          ret_val = PIO_enddef(pio_file)
          call PIO_closefile(pio_file)
          return
        end if
      end if

    End Subroutine test_create

    Subroutine test_open(test_id, err_msg)
    ! test_open():
    ! * Try to open file that doesn't exist, check for error
    ! * Open a file with PIO_write, write something, close
    ! * Open a file with PIO_nowrite, try to write, check for error
    ! * For netcdf / pnetcdf:
    !   - Try to open non-netcdf file, check for error
    ! Routines used in test: PIO_initdecomp, PIO_openfile, PIO_write_darray,
    !                        PIO_closefile, PIO_freedecomp
    ! Also uses PIO_createfile for binary tests
    !           PIO_redef, PIO_def_dim, PIO_def_var, PIO_enddef for [p]netcdf tests

      ! Input / Output Vars
      integer,                intent(in)  :: test_id
      character(len=str_len), intent(out) :: err_msg

      ! Local Vars
      character(len=str_len) :: filename
      integer                :: iotype, ret_val

      ! Data used to test writing
      integer,          dimension(3) :: data_to_write, compdof
      integer,          dimension(1) :: dims
      type(io_desc_t)                :: iodesc_nCells
      integer                        :: pio_dim
      type(var_desc_t)               :: pio_var

      err_msg = "no_error"
      dims(1) = 3*ntasks
      compdof = 3*my_rank+(/1,2,3/)  ! Where in the global array each task writes
      data_to_write = my_rank

      call PIO_initdecomp(pio_iosystem, PIO_int, dims, compdof, iodesc_nCells)

      filename = fnames(test_id)
      iotype   = iotypes(test_id)

      ! Open file that doesn't exist
      ret_val = PIO_openfile(pio_iosystem, pio_file, iotype, "FAKE.FILE", &
                             PIO_nowrite, CheckMPI=.false.)
      if (ret_val.eq.0) then
        ! Error in PIO_openfile
        err_msg = "Successfully opened file that doesn't exist"
        call PIO_closefile(pio_file)
        return
      end if

      ! Open existing file, write data to it (for binary file, need to create new file)
      if (is_netcdf(iotype)) then
        ret_val = PIO_openfile(pio_iosystem, pio_file, iotype, filename, PIO_write)
      else
        ret_val = PIO_createfile(pio_iosystem, pio_file, iotype, filename)
      end if
      if (ret_val.ne.0) then
        ! Error in PIO_openfile (or PIO_createfile)
        err_msg = "Could not open " // trim(filename) // " in write mode"
        return
      end if

      ! Enter define mode for netcdf files
      if (is_netcdf(iotype)) then
        ret_val = PIO_redef(pio_file)
        if (ret_val.ne.0) then
          ! Error in PIO_redef
          err_msg = "Could not enter redef mode"
          call PIO_closefile(pio_file)
          return
        end if

        ! Define a new dimension N
        ret_val = PIO_def_dim(pio_file, 'N', 3*ntasks, pio_dim)
        if (ret_val.ne.0) then
          ! Error in PIO_def_dim
          err_msg = "Could not define dimension N"
          call PIO_closefile(pio_file)
          return
        end if

        ! Define a new variable foo
        ret_val = PIO_def_var(pio_file, 'foo', PIO_int, &
                              (/pio_dim/), pio_var)
        if (ret_val.ne.0) then
          ! Error in PIO_def_var
          err_msg = "Could not define variable foo"
          call PIO_closefile(pio_file)
          return
        end if

        ! Leave define mode
        ret_val = PIO_enddef(pio_file)
        if (ret_val.ne.0) then
          ! Error in PIO_enddef
           print *,__FILE__,__LINE__,ret_val
          err_msg = "Could not end define mode"
          call PIO_closefile(pio_file)
          return
        end if
      end if

      ! Write foo
      call PIO_write_darray(pio_file, pio_var, iodesc_nCells, data_to_write, ret_val)
      if (ret_val.ne.0) then
        ! Error in PIO_write_darray
        err_msg = "Could not write data"
        call PIO_closefile(pio_file)
        return
      end if

      ! Close file
      call PIO_closefile(pio_file)

      ! Open existing file with PIO_nowrite, try to write (netcdf only)
      if (is_netcdf(iotype)) then
        ret_val = PIO_openfile(pio_iosystem, pio_file, iotype, filename, PIO_nowrite)
        if (ret_val.ne.0) then
          ! Error opening file
          err_msg = "Could not open file in NoWrite mode"
          return
        end if

        ! Try to write (should fail)
        call PIO_write_darray(pio_file, pio_var, iodesc_nCells, data_to_write, ret_val)
        if (ret_val.eq.0) then
          ! Error in PIO_write_darray
          err_msg = "Wrote to file opened in NoWrite mode"
          call PIO_closefile(pio_file)
          return
        end if
        ! Close file
        call PIO_closefile(pio_file)
      end if

      ! Try to open standard binary file as netcdf (if iotype = netcdf)
      if (is_netcdf(iotype)) then
        ret_val = PIO_openfile(pio_iosystem, pio_file, iotype, &
                               "not_netcdf.ieee", PIO_nowrite)
        if (ret_val.eq.0) then
          ! Error in PIO_openfile
          err_msg = "Opened a non-netcdf file as netcdf"
          call PIO_closefile(pio_file)
          return
        end if
      end if

      ! Free decomp
      call PIO_freedecomp(pio_iosystem, iodesc_nCells)

    End Subroutine test_open

end module basic_tests
