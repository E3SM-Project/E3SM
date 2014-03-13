!>
!! @file 
!! @brief User interface Module for PIO, this is the only file a user program should 'use'
!! 
!! $Revision$
!! $LastChangedDate$
!<

module pio
! Package all exposed variables and functions under one roof

! only pio_offset is intended for export from kinds

  use pio_kinds, only : pio_offset

  use piolib_mod, only : pio_initdecomp, &
       pio_openfile, pio_closefile, pio_createfile, pio_setdebuglevel, &
       pio_seterrorhandling, pio_setframe, pio_init, pio_get_local_array_size, &
       pio_freedecomp, pio_syncfile, &
       pio_finalize, pio_set_hint, pio_getnumiotasks, pio_file_is_open, &
       PIO_deletefile, PIO_get_numiotasks, PIO_get_iorank

  use pio_types, only : io_desc_t, file_desc_t, var_desc_t, iosystem_desc_t, &
	pio_int, pio_real, pio_double, pio_noerr, iotype_netcdf, &
	iotype_pnetcdf, iotype_binary, iotype_direct_pbinary, iotype_pbinary, &
        PIO_iotype_binary, PIO_iotype_direct_pbinary, PIO_iotype_pbinary, &
        pio_iotype_netcdf4p, pio_iotype_netcdf4c, pio_iotype_pnetcdf,pio_iotype_netcdf, &
	pio_global, pio_char, pio_write, pio_nowrite, pio_clobber, pio_noclobber, &
	pio_max_name, pio_max_var_dims, pio_rearr_none, &
#if defined(_NETCDF) || defined(_PNETCDF)
	pio_nofill, pio_unlimited, &
#endif
	pio_64bit_offset, pio_64bit_data, &
	pio_iotype_vdc2, &
        pio_rearr_box, pio_internal_error, pio_bcast_error, pio_return_error

  use piodarray, only : pio_read_darray, pio_write_darray, pio_set_buffer_size_limit  

  use pio_nf, only:        &
       PIO_enddef,            &
       PIO_inquire ,          &
       PIO_inq_attname ,  &
       PIO_inq_att ,          &
       PIO_inq_attlen ,    &
       PIO_inq_varid ,      &
       PIO_inq_varname ,  &
       PIO_inq_vartype ,  &
       PIO_inq_varndims ,&
       PIO_inq_vardimid ,&
       PIO_inq_varnatts ,&
       PIO_inq_dimid ,      &
       PIO_inq_dimname ,  &
       PIO_inq_dimlen ,    &
       PIO_def_dim   ,        &
       PIO_def_var   ,        &
       PIO_redef     ,          &
!       PIO_copy_att  ,       &   
       PIO_inquire_variable , &
       PIO_inquire_dimension 

  use pionfatt_mod, only : PIO_put_att   => put_att,        &
       PIO_get_att   => get_att
  use pionfput_mod, only : PIO_put_var   => put_var
  use pionfget_mod, only : PIO_get_var   => get_var
  use iso_c_binding

  implicit none
  public
contains
  subroutine pio_set_blocksize(blocksize)
    integer :: blocksize
    integer :: ierr
    interface
       integer(C_INT) function PIOc_set_blocksize(blocksize) &
            bind(C,name="PIOc_set_blocksize")
         use iso_c_binding
         integer(C_INT), intent(in), value :: blocksize
       end function PIOc_set_blocksize
    end interface
    ierr = PIOc_set_blocksize(blocksize)
  end subroutine pio_set_blocksize


  function pio_iam_iotask(iosystem) result(task)
    use iso_c_binding
    type(iosystem_desc_t), intent(in) :: iosystem
    logical :: task
    integer :: ierr
    logical(C_BOOL) :: ctask
    interface
       integer(C_INT) function PIOc_iam_iotask(iosysid, iotask) &
            bind(C,name="PIOc_iam_iotask")
         use iso_c_binding
         integer(C_INT), intent(in), value :: iosysid
         logical(C_BOOL), intent(out) :: iotask
       end function PIOc_iam_iotask
    end interface
    
    ierr = PIOc_iam_iotask(iosystem%iosysid, ctask)
    task = ctask
  end function pio_iam_iotask
  function pio_iotask_rank(iosystem) result(rank)
    type(iosystem_desc_t), intent(in) :: iosystem
    integer :: rank, ierr
    interface
       integer(C_INT) function PIOc_iotask_rank(iosysid, rank) &
            bind(C,name="PIOc_iotask_rank")
         use iso_c_binding
         integer(C_INT), intent(in), value :: iosysid
         integer(C_INT), intent(out) :: rank
       end function PIOc_iotask_rank
    end interface
    
    ierr = PIOc_iotask_rank(iosystem%iosysid, rank)
  end function pio_iotask_rank

  subroutine pio_iosystem_is_active(iosystem, active)
    use iso_c_binding
    type(iosystem_desc_t), intent(in) :: iosystem
    logical, intent(out) :: active
    logical(C_BOOL) :: lactive
    integer :: ierr
    interface
       integer(C_INT) function PIOc_iosystem_is_active(iosysid, active) &
            bind(C,name="PIOc_iosystem_is_active")
         use iso_c_binding
         integer(C_INT), intent(in), value :: iosysid
         logical(C_BOOL), intent(out) :: active
       end function PIOc_iosystem_is_active
    end interface

    ierr = PIOc_iosystem_is_active(iosystem%iosysid, lactive)
    active = lactive


  end subroutine pio_iosystem_is_active


end module pio

