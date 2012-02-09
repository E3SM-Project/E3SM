#include "dtypes.h"
#define __PIO_FILE__ "pio_msg_callbacks.F90"
subroutine pio_callback_handler(iosystem, msg)
  use pio
  use pio_msg_mod
  use pio_support, only : debugAsync, piodie
#ifndef NO_MPIMOD
  use mpi !_EXTERNAL
#endif
  implicit none
#ifdef NO_MPIMOD
  include 'mpif.h' !_EXTERNAL
#endif
  type(iosystem_desc_t) :: iosystem
  integer, intent(in) :: msg

  type(file_desc_t), pointer :: file
  integer fh, ierr


  integer :: len, id, dimids(PIO_MAX_VAR_DIMS), type, info(4)
  character(len=pio_max_name) :: name
  type(var_desc_t) :: vardesc

  call mpi_bcast(fh, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)
  file=> lookupfile(fh)


  select case(msg)
  case (PIO_MSG_CLOSE_FILE)
     call delete_from_file_list(fh)
     call pio_closefile(file)
     deallocate(file)
  case (PIO_MSG_DEF_DIM)
     ierr = pio_def_dim(file, name, len, id)   
  case (PIO_MSG_DEF_VAR)
     ierr = pio_def_var(file, name, type, dimids, vardesc)
  case (PIO_MSG_ENDDEF)
     ierr = pio_enddef(file)
  case (PIO_MSG_REDEF)
     ierr = pio_redef(file)
  case (PIO_MSG_INQ_VARNDIMS)
     ierr = pio_inq_varndims(file, id, len)
  case (PIO_MSG_INQ_VARNATTS)
     ierr = pio_inq_varnatts(file, id, len)
  case (PIO_MSG_INQ_VARDIMID)
     ierr = pio_inq_vardimid(file, id, dimids)
  case (PIO_MSG_INQ_VARID)
     ierr = pio_inq_varid(File, name, id)
  case (PIO_MSG_INQ_VARNAME)
     ierr = pio_inq_varname(file, id, name)
  case (PIO_MSG_INQ_VARTYPE)
     ierr = pio_inq_vartype(file, id, type)
  case (PIO_MSG_INQ_DIMID)
     ierr = pio_inq_dimid(file, name, id)
  case (PIO_MSG_INQ_DIMLEN)
     ierr = pio_inq_dimlen(file, id, len)
  case (PIO_MSG_INQ_DIMNAME)
     ierr = pio_inq_dimname(file, id, name)
  case(PIO_MSG_INQUIRE)
     ierr = pio_inquire(file, info(1), info(2), info(3), info(4))
  case(PIO_MSG_INQ_ATT)
     ierr = pio_inq_att(file, id, name, type, len)
  case(PIO_MSG_INQ_ATTNAME)
     ierr = pio_inq_attname(file, id, type, name)
  case(PIO_MSG_INQ_ATTLEN)
     ierr = pio_inq_attlen(file, id, name, len)
  case(PIO_MSG_SYNC_FILE)
     call pio_syncfile(file)
  case default
     print *, 'PIO Got unrecognized message ', msg, ierr
     call piodie(__PIO_FILE__,__LINE__)  
  end select

end subroutine pio_callback_handler

subroutine create_file_handler(iosystem)
  use pio, only : iosystem_desc_t, file_desc_t, pio_createfile
  use pio_kinds, only : char_len
  use pio_msg_mod, only : add_to_file_list
  use pio_support, only : debugAsync
#ifndef NO_MPIMOD
  use mpi !_EXTERNAL
#endif
  implicit none
#ifdef NO_MPIMOD
  include 'mpif.h' !_EXTERNAL
#endif
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
  if(Debugasync) print *,__PIO_FILE__,__LINE__,file%fh

end subroutine create_file_handler

subroutine open_file_handler(iosystem)
  use pio
  use piolib_mod
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
#ifndef NO_MPIMOD
  use mpi !_EXTERNAL
#endif
  implicit none
#ifdef NO_MPIMOD
  include 'mpif.h' !_EXTERNAL
#endif
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
  if(Debugasync) print *,__PIO_FILE__,__LINE__,file%fh

end subroutine open_file_handler

subroutine initdecomp_dof_handler(iosystem)

  use pio
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
#ifndef NO_MPIMOD
  use mpi !_EXTERNAL
#endif
  implicit none
#ifdef NO_MPIMOD
  include 'mpif.h' !_EXTERNAL
#endif
  type(iosystem_desc_t) :: iosystem

  type(io_desc_t), pointer :: iodesc
  integer :: ierr
  integer(i4) :: basepiotype, dims_size, dims(PIO_MAX_VAR_DIMS), dof_size, sandc_size
  integer(kind=pio_offset) :: compdof(1)
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
#ifndef NO_MPIMOD
  use mpi !_EXTERNAL
#endif
  implicit none
#ifdef NO_MPIMOD
  include 'mpif.h' !_EXTERNAL
#endif
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
#ifndef _MPISERIAL
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
#endif

end subroutine writedarray_handler


subroutine readdarray_handler(iosystem)
  use pio
  use pio_kinds
  use pio_msg_mod
  use pio_support, only : debugAsync
#ifndef NO_MPIMOD
  use mpi !_EXTERNAL
#endif
  implicit none
#ifdef NO_MPIMOD
  include 'mpif.h' !_EXTERNAL
#endif
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
#ifndef _MPISERIAL
  select case(type)
  case(mpi_integer)
     call pio_read_darray(file, v, iodesc, aint, ierr)
  case(mpi_real4)
     call pio_read_darray(file, v, iodesc, areal, ierr)
  case(mpi_real8)
     call pio_read_darray(file, v, iodesc, adouble, ierr)
  end select
#endif

end subroutine readdarray_handler

subroutine seterrorhandling_handler(ios)
  use pio, only : iosystem_desc_t, pio_seterrorhandling
#ifndef NO_MPIMOD
  use mpi !_EXTERNAL
#endif
  implicit none
#ifdef NO_MPIMOD
  include 'mpif.h' !_EXTERNAL
#endif 
  type(iosystem_desc_t), intent(inout) :: ios
  integer :: method, ierr

  call mpi_bcast(method, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  
  call pio_seterrorhandling(ios, method)

end subroutine seterrorhandling_handler

subroutine string_handler_for_att(file, varid, name, strlen, msg)
  use pio_msg_mod, only : pio_msg_getatt
  use pio, only : file_desc_t, pio_get_att, pio_put_att
  use pio_support, only : debugasync
  implicit none

  type(file_desc_t) :: file
  integer, intent(in) :: varid, strlen, msg
  character(len=*) :: name
  character(len=strlen) :: str
  integer :: ierr

  if(msg==PIO_MSG_GETATT) then
     if(Debugasync) print *,__PIO_FILE__,__LINE__, varid, name
     ierr = pio_get_att(file, varid, name, str )  
     if(Debugasync) print *,__PIO_FILE__,__LINE__, str
  else
     ierr = pio_put_att(file, varid, name, str )  
  end if
end subroutine string_handler_for_att

subroutine att_handler(ios, msg)
  
  use pio, only : iosystem_desc_t, file_desc_t, pio_get_att, pio_max_name, pio_put_att
  use pio_kinds, only : i4, r4, r8
  use pio_msg_mod, only : lookupfile, pio_msg_putatt, pio_msg_getatt
  use pio_support, only : debugAsync, piodie
#ifndef NO_MPIMOD
  use mpi !_EXTERNAL
#endif
  implicit none
#ifdef NO_MPIMOD
  include 'mpif.h' !_EXTERNAL
#endif
  integer, intent(in) :: msg
  type(iosystem_desc_t), intent(inout) :: ios
  type(file_desc_t), pointer :: file
  integer :: fh, varid, ierr, itype, strlen, nlen
  character(len=PIO_MAX_NAME) :: name

  real(r4) :: rvar
  real(r8) :: dvar
  integer(i4) :: ivar

  if(Debugasync) print *,__PIO_FILE__,__LINE__
  
  call mpi_bcast(fh, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  call mpi_bcast(varid, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  call mpi_bcast(itype, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  call mpi_bcast(nlen, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  call mpi_bcast(name(1:nlen), nlen, mpi_integer, ios%compmaster, ios%intercomm, ierr)
  if(Debugasync) print *,__PIO_FILE__,__LINE__, itype,nlen

  file=> lookupfile(fh)
  
  select case(itype)
  case (TYPETEXT)
     call mpi_bcast(strlen, 1, mpi_integer, ios%compmaster, ios%intercomm, ierr)
       if(Debugasync) print *,__PIO_FILE__,__LINE__, strlen,nlen
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
#ifndef NO_MPIMOD
  use mpi !_EXTERNAL
#endif
  implicit none
#ifdef NO_MPIMOD
  include 'mpif.h' !_EXTERNAL
#endif

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
