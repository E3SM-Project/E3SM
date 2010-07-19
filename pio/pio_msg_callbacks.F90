subroutine create_file_handler(iosystem)
  use mpi !_EXTERNAL
  use pio, only : iosystem_desc_t, file_desc_t, pio_createfile
  use pio_kinds, only : char_len
  use pio_msg_mod, only : add_to_file_list
  use pio_support, only : debugAsync
  implicit none
  type(iosystem_desc_t) :: iosystem

  integer :: ierr
  integer :: iotype, amode
  
  character(len=char_len) :: fname
  type(file_desc_t), pointer :: file
  
  call mpi_bcast(fname, char_len, mpi_character, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(iotype, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  
  call mpi_bcast(amode, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)

  allocate(file)
  
  ierr= pio_createfile(iosystem, file, iotype, trim(fname), amode )

  call add_to_file_list(file)


end subroutine create_file_handler

subroutine open_file_handler(iosystem)
  use mpi !_EXTERNAL
  use pio
  use piolib_mod
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
  implicit none
  type(iosystem_desc_t) :: iosystem

  integer :: ierr
  integer :: iotype, amode
  
  character(len=char_len) :: fname
  type(file_desc_t), pointer :: file
  
  call mpi_bcast(fname, char_len, mpi_character, iosystem%compmaster, iosystem%intercomm, ierr)

  call mpi_bcast(iotype, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  
  call mpi_bcast(amode, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)

  allocate(file)
  
  ierr= pio_openfile(iosystem, file, iotype, trim(fname), amode)

  call add_to_file_list(file)

end subroutine open_file_handler

subroutine def_dim_handler(iosystem)
  use mpi !_EXTERNAL
  use pio
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
  implicit none

  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: ierr
  
  integer :: fh, len, clen, dimid
  character(len=PIO_MAX_NAME) :: dimname

  if(debugAsync) print *,__FILE__,__LINE__

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(debugAsync) print *,__FILE__,__LINE__,fh
  call mpi_bcast(len, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(debugAsync) print *,__FILE__,__LINE__,len
  call mpi_bcast(clen, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(debugAsync) print *,__FILE__,__LINE__,clen
  call mpi_bcast(dimname(1:clen), clen, mpi_character, iosystem%compmaster, iosystem%intercomm, ierr)
  if(debugAsync) print *,__FILE__,__LINE__,dimname(1:clen)

  file=> lookupfile(fh)

  if(debugAsync) print *, __FILE__,__LINE__,file%fh, fh

  
  ierr = pio_def_dim(file, dimname(1:clen), len, dimid)
  
  call mpi_bcast(dimid, 1, mpi_integer, iosystem%iomaster, iosystem%intercomm, ierr)

end subroutine def_dim_handler

subroutine def_var_handler(iosystem)
  use mpi !_EXTERNAL
  use pio
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
  implicit none

  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: ierr
  
  integer :: fh, type, clen, dimcnt, dimids(pio_max_var_dims)
  type(var_desc_t) :: vardesc
  
  character(len=PIO_MAX_NAME) :: varname


  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(type, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(clen, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(varname(1:clen), clen, mpi_character, iosystem%compmaster, iosystem%intercomm, ierr)

  call mpi_bcast(dimcnt, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(dimids(1:dimcnt), dimcnt, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)

  if(debugAsync) print *,__FILE__,__LINE__,varname(1:clen)
  
  file=> lookupfile(fh)

  if(debugAsync) print *, __FILE__,__LINE__,file%fh, fh

  
  ierr = pio_def_var(file, varname(1:clen),type, dimids(1:dimcnt), vardesc)
  

end subroutine def_var_handler

subroutine enddef_handler(iosystem)
  use pio
  use pio_msg_mod
  use mpi ! _EXTERNAL
  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file

  integer :: fh
  integer :: ierr

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)

  file=> lookupfile(fh)

  ierr = pio_enddef(file) 

end subroutine enddef_handler

subroutine initdecomp_dof_handler(iosystem)
  use mpi !_EXTERNAL
  use pio
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
  implicit none

  type(iosystem_desc_t) :: iosystem

  type(io_desc_t), pointer :: iodesc
  integer :: ierr
  integer(i4) :: basepiotype, dims_size, dims(PIO_MAX_VAR_DIMS), dof_size, sandc_size
  integer(i4) :: compdof(1)
  integer(kind=pio_offset), allocatable :: iostart(:), iocount(:)

  call mpi_bcast(basepiotype, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(dims_size, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(dims(1:dims_size), dims_size, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)  
  
  allocate(iodesc)

  compdof=0
  call add_to_iodesc_list(iodesc)

  call mpi_bcast(iodesc%async_id, 1, mpi_integer, iosystem%iomaster, iosystem%intercomm, ierr)  


  call pio_initdecomp(iosystem, basepiotype, dims(1:dims_size), compdof, iodesc)


end subroutine initdecomp_dof_handler

subroutine writedarray_handler(iosystem)
  use mpi !_EXTERNAL
  use pio
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
  implicit none

  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  type(var_desc_t) :: v
  type(io_desc_t), pointer :: iodesc

  integer :: ierr, type, fh, fillv, iod_id

  integer(i4) :: fillval_int, aint(1)
  real(r4) :: fillval_real, areal(1)
  real(r8) :: fillval_double, adouble(1)
  



  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(v%varid, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(v%rec, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(v%ndims, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(iod_id, 1, mpi_integer , iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(type, 1, mpi_integer , iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(fillv, 1, mpi_integer , iosystem%compmaster, iosystem%intercomm, ierr)

  file=> lookupfile(fh)
  iodesc => lookupiodesc(iod_id)

  select case(type)
  case(mpi_integer)
     if(fillv==1) then
        call mpi_bcast(fillval_int, 1, mpi_integer , iosystem%compmaster, iosystem%intercomm, ierr)
        call pio_write_darray(file, v, iodesc, aint, ierr, fillval_int)
     else
        call pio_write_darray(file, v, iodesc, aint, ierr)
     end if
  case(mpi_real4)
     if(fillv==1) then
        call mpi_bcast(fillval_real, 1, mpi_real4 , iosystem%compmaster, iosystem%intercomm, ierr)
        call pio_write_darray(file, v, iodesc, areal, ierr, fillval_real)
     else
        call pio_write_darray(file, v, iodesc, areal, ierr)
     end if
  case(mpi_real8)
     if(fillv==1) then
        call mpi_bcast(fillval_double, 1, mpi_real8, iosystem%compmaster, iosystem%intercomm, ierr)
        call pio_write_darray(file, v, iodesc, adouble, ierr, fillval_double)
     else
        call pio_write_darray(file, v, iodesc, adouble, ierr)
     end if
  end select


end subroutine writedarray_handler


subroutine readdarray_handler(iosystem)
  use mpi !_EXTERNAL
  use pio
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
  implicit none

  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  type(var_desc_t) :: v
  type(io_desc_t), pointer :: iodesc

  integer :: ierr, type, fh, iod_id

  integer(i4) :: aint(1)
  real(r4) :: areal(1)
  real(r8) :: adouble(1)
  


  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(v%varid, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(v%rec, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(iod_id, 1, mpi_integer , iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(type, 1, mpi_integer , iosystem%compmaster, iosystem%intercomm, ierr)

  file=> lookupfile(fh)
  iodesc => lookupiodesc(iod_id)

  select case(type)
  case(mpi_integer)
     call pio_read_darray(file, v, iodesc, aint, ierr)
  case(mpi_real4)
     call pio_read_darray(file, v, iodesc, areal, ierr)
  case(mpi_real8)
     call pio_read_darray(file, v, iodesc, adouble, ierr)
  end select


end subroutine readdarray_handler

subroutine close_file_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_closefile
  use pio_msg_mod, only: lookupfile, delete_from_file_list
  use mpi !_EXTERNAL
  use pio_support, only : debugAsync
  implicit none
  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, ierr

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  file=> lookupfile(fh)

  call delete_from_file_list(fh)
  call pio_closefile(file)
  deallocate(file)


end subroutine close_file_handler

subroutine inq_varndims_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_inq_varndims
  use mpi !_EXTERNAL
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, ierr, varid, ndims

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  call mpi_bcast(varid, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,varid
  file=> lookupfile(fh)

  ierr =  pio_inq_varndims(file, varid, ndims)



end subroutine inq_varndims_handler

subroutine inq_varid_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_inq_varid, pio_max_name
  use mpi !_EXTERNAL
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, ierr, nlen, varid
  character(len=PIO_MAX_NAME) :: name

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  call mpi_bcast(nlen, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(name(1:nlen), nlen, mpi_character, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,name(1:nlen)
  file=> lookupfile(fh)

  ierr =  pio_inq_varid(file, name(1:nlen), varid)



end subroutine inq_varid_handler

subroutine finalize_handler(iosystem)
  use pio, only : iosystem_desc_t, pio_finalize
  use pio_support, only : debugAsync
  implicit none
  type(iosystem_desc_t) :: iosystem
  integer :: ierr

  call pio_finalize(iosystem, ierr)

end subroutine finalize_handler
