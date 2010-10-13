#include "dtypes.h"

subroutine create_file_handler(iosystem)
  use pio, only : iosystem_desc_t, file_desc_t, pio_createfile
  use pio_kinds, only : char_len
  use pio_msg_mod, only : add_to_file_list
  use pio_support, only : debugAsync
  implicit none
  include 'mpif.h' !_EXTERNAL

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
  if(Debugasync) print *,__FILE__,__LINE__,file%fh

end subroutine create_file_handler

subroutine open_file_handler(iosystem)
  use pio
  use piolib_mod
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
  implicit none
  include 'mpif.h' !_EXTERNAL

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
  if(Debugasync) print *,__FILE__,__LINE__,file%fh

end subroutine open_file_handler

subroutine def_dim_handler(iosystem)
  use pio
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
  implicit none
  include 'mpif.h' !_EXTERNAL

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

  
end subroutine def_dim_handler

subroutine def_var_handler(iosystem)
  use pio
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
  implicit none
  include 'mpif.h' !_EXTERNAL

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
  implicit none
  include 'mpif.h' !_EXTERNAL

  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file

  integer :: fh
  integer :: ierr

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)

  file=> lookupfile(fh)

  ierr = pio_enddef(file) 

end subroutine enddef_handler

subroutine initdecomp_dof_handler(iosystem)

  use pio
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
  implicit none
  include 'mpif.h' !_EXTERNAL

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
  use pio
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
  implicit none
  include 'mpif.h' !_EXTERNAL

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
  use pio
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
  implicit none
  include 'mpif.h' !_EXTERNAL

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
  use pio_support, only : debugAsync
  implicit none
  include 'mpif.h' !_EXTERNAL


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
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL


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


subroutine sync_file_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_syncfile
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL


  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, ierr, varid, ndims

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  file=> lookupfile(fh)

  call pio_syncfile(file)



end subroutine sync_file_handler



subroutine inq_att_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_inq_att, pio_max_name
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL


  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, ierr, namelen, xtype, alen, varid
  character(len=PIO_MAX_NAME) :: name

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  call mpi_bcast(varid, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,varid
  call mpi_bcast(namelen, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,varid
  call mpi_bcast(name, namelen, mpi_character, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,varid
  file=> lookupfile(fh)

  ierr =  pio_inq_att(file, varid, name(1:namelen),xtype, alen)

end subroutine inq_att_handler

subroutine inq_attname_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_inq_attname, pio_max_name
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL


  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, ierr,varid, attnum
  character(len=PIO_MAX_NAME) :: name

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  call mpi_bcast(varid, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,varid
  call mpi_bcast(attnum, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,varid

  file=> lookupfile(fh)

  ierr =  pio_inq_attname(file, varid, attnum, name)

end subroutine inq_attname_handler

subroutine inq_attlen_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_inq_attlen, pio_max_name
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL


  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, ierr,varid, alen, namelen
  character(len=PIO_MAX_NAME) :: name

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  call mpi_bcast(varid, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,varid
  call mpi_bcast(namelen, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,varid
  call mpi_bcast(name, namelen, mpi_character, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,varid

  file=> lookupfile(fh)

  ierr =  pio_inq_attlen(file, varid, name(1:namelen), alen)

end subroutine inq_attlen_handler


subroutine inq_varnatts_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_inq_varnatts
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL


  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, ierr, varid, natts

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  call mpi_bcast(varid, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,varid
  file=> lookupfile(fh)

  ierr =  pio_inq_varnatts(file, varid, natts)



end subroutine inq_varnatts_handler


subroutine inq_vartype_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_inq_vartype
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL


  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, ierr, varid, type

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  call mpi_bcast(varid, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,varid
  file=> lookupfile(fh)

  ierr =  pio_inq_vartype(file, varid, type)


end subroutine inq_vartype_handler

subroutine inq_varid_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_inq_varid, pio_max_name
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL

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

subroutine inq_varname_internal(file, varid, namelen)
  use pio, only : file_desc_t, pio_inq_varname
  use pio_support, only : debugAsync

  type(file_desc_t) :: file
  integer, intent(in) ::  varid, namelen
  character(len=namelen) :: name  
  integer :: ierr

  if(Debugasync) print *,__FILE__,__LINE__,varid
  ierr =  pio_inq_varname(file, varid, name)

end subroutine inq_varname_internal

subroutine inq_varname_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL

  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, varid, namelen, ierr


  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  call mpi_bcast(varid, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(namelen, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)

  file=> lookupfile(fh)

  call inq_varname_internal(file, varid, namelen)


end subroutine inq_varname_handler

subroutine inq_dimname_internal(file, dimid, namelen)
  use pio, only : file_desc_t, pio_inq_dimname
  use pio_support, only : debugAsync

  type(file_desc_t) :: file
  integer, intent(in) ::  dimid, namelen
  character(len=namelen) :: name  
  integer :: ierr

  if(Debugasync) print *,__FILE__,__LINE__,dimid
  ierr =  pio_inq_dimname(file, dimid, name)

end subroutine inq_dimname_internal

subroutine inq_dimname_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL

  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, dimid, namelen, ierr


  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  call mpi_bcast(dimid, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(namelen, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)

  file=> lookupfile(fh)

  call inq_dimname_internal(file, dimid, namelen)


end subroutine inq_dimname_handler





subroutine inq_dimid_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_inq_dimid, pio_max_name
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL

  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, ierr, nlen, dimid
  character(len=PIO_MAX_NAME) :: name

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  call mpi_bcast(nlen, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(name(1:nlen), nlen, mpi_character, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,name(1:nlen)
  file=> lookupfile(fh)

  ierr =  pio_inq_dimid(file, name(1:nlen), dimid)



end subroutine inq_dimid_handler

subroutine inq_vardimid_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_inq_vardimid
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL

  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, ierr, size_dimids, varid
  integer, allocatable :: dimids(:)

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  call mpi_bcast(varid, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  call mpi_bcast(size_dimids, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  file=> lookupfile(fh)
  allocate(dimids(size_dimids))
  ierr =  pio_inq_vardimid(file, varid, dimids)
  deallocate(dimids)
end subroutine inq_vardimid_handler

subroutine inq_dimlen_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_inq_dimlen
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL

  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, ierr, dimlen, dimid

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  call mpi_bcast(dimid, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  file=> lookupfile(fh)

  ierr =  pio_inq_dimlen(file, dimid, dimlen)
end subroutine inq_dimlen_handler

subroutine inquire_handler(iosystem)
  use pio, only: iosystem_desc_t, file_desc_t, pio_inquire
  use pio_msg_mod, only : lookupfile
  use pio_support, only : debugAsync
  
  implicit none
  include 'mpif.h' !_EXTERNAL

  type(iosystem_desc_t) :: iosystem
  type(file_desc_t), pointer :: file
  integer :: fh, ierr, nDimensions,nVariables,nAttributes,unlimitedDimID

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__,fh
  file=> lookupfile(fh)

  ierr =  pio_inquire(file,nDimensions,nVariables,nAttributes,unlimitedDimID)
end subroutine inquire_handler

subroutine seterrorhandling_handler(ios)
  use pio, only : iosystem_desc_t, pio_seterrorhandling
  implicit none
  include 'mpif.h' !_EXTERNAL
  type(iosystem_desc_t), intent(inout) :: ios
  integer :: method, ierr

  call mpi_bcast(method, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  
  call pio_seterrorhandling(ios, method)

end subroutine seterrorhandling_handler

subroutine string_handler_for_att(file, varid, name, strlen, msg)
  use pio_msg_mod, only : pio_msg_getatt
  use pio, only : file_desc_t, pio_get_att, pio_put_att
  use pio_support, only : debugasync
  type(file_desc_t) :: file
  integer, intent(in) :: varid, strlen, msg
  character(len=*) :: name
  character(len=strlen) :: str
  if(msg==PIO_MSG_GETATT) then
     if(Debugasync) print *,__FILE__,__LINE__, varid, name
     ierr = pio_get_att(file, varid, name, str )  
     if(Debugasync) print *,__FILE__,__LINE__, str
  else
     ierr = pio_put_att(file, varid, name, str )  
  end if
end subroutine string_handler_for_att

subroutine att_handler(ios, msg)
  
  use pio, only : iosystem_desc_t, file_desc_t, pio_get_att, pio_max_name, pio_put_att
  use pio_kinds, only : i4, r4, r8
  use pio_msg_mod, only : lookupfile, pio_msg_putatt, pio_msg_getatt
  use pio_support, only : debugAsync, piodie
  implicit none
  include 'mpif.h' !_EXTERNAL
  integer, intent(in) :: msg
  type(iosystem_desc_t), intent(inout) :: ios
  type(file_desc_t), pointer :: file
  integer :: fh, varid, ierr, itype, strlen, nlen
  character(len=PIO_MAX_NAME) :: name

  real(r4) :: rvar
  real(r8) :: dvar
  integer(i4) :: ivar

  if(Debugasync) print *,__FILE__,__LINE__
  
  call mpi_bcast(fh, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  call mpi_bcast(varid, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  call mpi_bcast(itype, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  call mpi_bcast(nlen, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  call mpi_bcast(name(1:nlen), nlen, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  if(Debugasync) print *,__FILE__,__LINE__, itype,nlen

  file=> lookupfile(fh)
  
  select case(itype)
  case (TYPETEXT)
     call mpi_bcast(strlen, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
       if(Debugasync) print *,__FILE__,__LINE__, strlen,nlen
     call string_handler_for_att (file, varid, name(1:nlen), strlen, msg)
  case (TYPEREAL)
     if(msg==PIO_MSG_GETATT) then
        ierr = pio_get_att(file, varid, name(1:nlen), rvar)
     else
        ierr = pio_put_att(file, varid, name(1:nlen), rvar)
     end if
  case (TYPEDOUBLE)
     if(msg==PIO_MSG_GETATT) then
        ierr = pio_get_att(file, varid, name(1:nlen), dvar)
     else
        ierr = pio_put_att(file, varid, name(1:nlen), dvar)
     end if
  case (TYPEINT)
     if(msg==PIO_MSG_GETATT) then
        ierr = pio_get_att(file, varid, name(1:nlen), ivar)
     else
        ierr = pio_put_att(file, varid, name(1:nlen), ivar)
     end if
  end select
  
end subroutine att_handler



subroutine att_1d_handler(ios, msg)
  
  use pio, only : iosystem_desc_t, file_desc_t, pio_get_att, pio_max_name, pio_put_att
  use pio_kinds, only : i4, r4, r8
  use pio_msg_mod, only : lookupfile, pio_msg_getatt_1d, pio_msg_putatt_1d
  use pio_support, only : debugAsync, piodie
  implicit none
  include 'mpif.h' !_EXTERNAL
  type(iosystem_desc_t), intent(inout) :: ios
  integer, intent(in) :: msg
  type(file_desc_t), pointer :: file
  integer :: fh, varid, ierr, itype, strlen, nlen, clen
  character(len=PIO_MAX_NAME) :: name

  real(r4), allocatable :: rvar(:)
  real(r8), allocatable :: dvar(:)
  integer(i4), allocatable :: ivar(:)
  
  call mpi_bcast(fh, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  call mpi_bcast(varid, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  call mpi_bcast(itype, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  call mpi_bcast(nlen, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  call mpi_bcast(name(1:nlen), nlen, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  call MPI_BCAST(clen,1,MPI_INTEGER,ios%CompMaster, ios%intercomm , ierr)

  file=> lookupfile(fh)
  
  select case(itype)
  case (TYPEREAL)
     allocate(rvar(clen))
     if(msg==pio_msg_getatt_1d) then
        ierr = pio_get_att(file, varid, name(1:nlen), rvar)
     else
        ierr = pio_put_att(file, varid, name(1:nlen), rvar)
     end if
     deallocate(rvar)
  case (TYPEDOUBLE)
     allocate(dvar(clen))
     if(msg==pio_msg_getatt_1d) then
        ierr = pio_get_att(file, varid, name(1:nlen), dvar)
     else
        ierr = pio_put_att(file, varid, name(1:nlen), dvar)
     end if
     deallocate(dvar)
  case (TYPEINT)
     allocate(ivar(clen))
     if(msg==pio_msg_getatt_1d) then
        ierr = pio_get_att(file, varid, name(1:nlen), ivar)
     else
        ierr = pio_put_att(file, varid, name(1:nlen), ivar)
     end if
     deallocate(ivar)
  end select
  
end subroutine att_1d_handler


subroutine finalize_handler(iosystem)
  use pio, only : iosystem_desc_t, pio_finalize
  use pio_support, only : debugAsync
  implicit none
  type(iosystem_desc_t) :: iosystem
  integer :: ierr

  call pio_finalize(iosystem, ierr)

end subroutine finalize_handler
