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
       pio_freedecomp, pio_syncfile,pio_numtowrite,pio_numtoread,pio_setiotype, &
       pio_dupiodesc, pio_finalize, pio_set_hint, pio_getnumiotasks, pio_file_is_open, &
       pio_setnum_OST, pio_getnum_OST

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
#if defined(_PNETCDF)
	pio_64bit_offset, pio_64bit_data, &
#endif
	pio_iotype_vdc2, &
        pio_rearr_box, pio_internal_error, pio_bcast_error, pio_return_error

  use piodarray, only : pio_read_darray, pio_write_darray, pio_set_buffer_size_limit  

  use nf_mod, only:        &
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
       PIO_copy_att  ,       &
       PIO_inquire_variable , &
       PIO_inquire_dimension 

  use pionfatt_mod, only : PIO_put_att   => put_att,        &
       PIO_get_att   => get_att
  use pionfput_mod, only : PIO_put_var   => put_var
  use pionfget_mod, only : PIO_get_var   => get_var

  use calcdecomp, only : pio_set_blocksize
   


  implicit none
  public
contains
  function pio_iam_iotask(iosystem) result(task)
    type(iosystem_desc_t), intent(in) :: iosystem
    logical :: task
    task = iosystem%ioproc
  end function pio_iam_iotask
  function pio_iotask_rank(iosystem) result(rank)
    type(iosystem_desc_t), intent(in) :: iosystem
    integer :: rank
    rank = iosystem%io_rank
  end function pio_iotask_rank

end module pio

