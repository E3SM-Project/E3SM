subroutine wrap_create (path, cmode, ncid)
  implicit none
  include 'netcdf.inc'
      
  character*(*) path
  integer cmode, ncid

  integer ret

  ret =      nf_create (path, cmode, ncid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_open (path, omode, ncid)
  implicit none
  include 'netcdf.inc'

  character*(*) path
  integer omode
  integer ncid
  integer ret
  
  ret =      nf_open (path, omode, ncid)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_OPEN: nf_open failed for file ',path
    call handle_error (ret)
  end if
end subroutine

subroutine wrap_inq_varid (nfid, varname, varid)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character*(*) varname

  integer ret

  ret = nf_inq_varid (nfid, varname, varid)
  if (ret.ne.NF_NOERR) then
    write(6,*)'nf_inq_varid: ',trim(varname),' not found'
    call handle_error (ret)
  end if
end subroutine wrap_inq_varid

subroutine wrap_inq_dimname (nfid, dimid, dimname)
  implicit none
  include 'netcdf.inc'
  
  integer nfid, dimid
  character*(*) dimname
  
  integer ret
  
  ret =      nf_inq_dimname (nfid, dimid, dimname)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_inq_dimlen (nfid, dimid, dimlen)
  implicit none
  include 'netcdf.inc'

  integer nfid, dimid, dimlen

  integer ret

  ret = nf_inq_dimlen (nfid, dimid, dimlen)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_inq_dimlen

subroutine wrap_inq_dimid (nfid, dimname, dimid)
  implicit none
  include 'netcdf.inc'

  integer nfid, dimid
  character*(*) dimname

  integer ret

  ret = nf_inq_dimid (nfid, dimname, dimid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_inq_dimid

subroutine wrap_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, xtype, ndims, dimids(nf_max_dims), natts
  character*(*) varname

  integer ret

  ret = nf_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_inq_var

subroutine wrap_def_dim (nfid, dimname, len, dimid)
  implicit none
  include 'netcdf.inc'

  integer nfid, len, dimid
  character*(*) dimname

  integer ret

  ret = nf_def_dim (nfid, dimname, len, dimid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_def_dim

subroutine wrap_def_var (nfid, name, xtype, nvdims, vdims, varid)
  implicit none
  include 'netcdf.inc'

  integer nfid, xtype, nvdims, varid
  integer vdims(nvdims)
  character*(*) name
  
  integer ret

  ret =      nf_def_var (nfid, name, xtype, nvdims, vdims, varid)
!      write(6,*)'WRAP_DEF_VAR: ',name,' has varid ',varid,nvdims,
!     $'dimensions of ids ',(vdims(i),i=1,nvdims)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_get_var_double (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  real*8 arr(*)

  integer ret

  ret = nf_get_var_double (nfid, varid, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_get_var_double

subroutine wrap_get_var_int (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer arr(*)

  integer ret

  ret = nf_get_var_int (nfid, varid, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_get_var_int

subroutine wrap_get_vara_double (nfid, varid, start, count, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(nf_max_dims), count(nf_max_dims)
  real*8 arr(*)

  integer ret

  ret = nf_get_vara_double (nfid, varid, start, count, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_get_vara_double

subroutine wrap_get_vara_int (nfid, varid, start, count, arr)
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(*), count(*)
  integer arr(*)

  integer ret

  ret =      nf_get_vara_int (nfid, varid, start, count, arr)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_GET_VARA8: error reading varid =', varid
    call handle_error (ret)
  end if
end subroutine

subroutine wrap_put_vara_text (nfid, varid, start, count, text)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), count(*)
  character*(*) text

  integer ret

  ret =      nf_put_vara_text (nfid, varid, start, count, text)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_PUT_VARA_TEXT: error writing varid =', varid
    call handle_error (ret)
  end if
end subroutine

!------------------------------------------------------------------------------

subroutine wrap_get_vara_text (nfid, varid, start, count, text)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), count(*)
  character*(*) text

  integer ret

  ret =      nf_get_vara_text (nfid, varid, start, count, text)
  if (ret.ne.NF_NOERR) then
    write(6,*)'WRAP_GET_VARA_TEXT: error writing varid =', varid
    call handle_error (ret)
  end if
end subroutine

subroutine wrap_put_var8 (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  real*8 arr(*)

  integer ret
  ret = nf_put_var_double (nfid, varid, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_put_var8

subroutine wrap_put_vara_double (nfid, varid, start, count, arr)
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), count(*)
  real(r8) arr(*)
  
  integer ret
  ret =      nf_put_vara_double (nfid, varid, start, count, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine handle_error (ret)
  implicit none
  include 'netcdf.inc'

  integer ret

  write(6,*) nf_strerror (ret)
  call abort
end subroutine handle_error

subroutine wrap_put_vara_int (nfid, varid, start, count, arr)
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(*), count(*)
  integer arr(*)

  integer ret
  ret =      nf_put_vara_int (nfid, varid, start, count, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_inq (nfid, ndims, nvars, ngatts, unlimdimid)
  implicit none
  include 'netcdf.inc'
  
  integer nfid, ndims, nvars, ngatts, unlimdimid
  
  integer ret

  ret =      nf_inq (nfid, ndims, nvars, ngatts, unlimdimid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_inq_dim (nfid, dimid, name, length)
  implicit none
  include 'netcdf.inc'
  
  integer nfid, dimid, length
  character*(*) name
  
  integer ret

  ret =      nf_inq_dim (nfid, dimid, name, length)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_inq_attname (nfid, varid, num, attname)
  implicit none
      
  include 'netcdf.inc'

  integer nfid, varid, num
  character*(*) attname

  integer ret

  ret =      nf_inq_attname (nfid, varid, num, attname)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_get_att_text (nfid, varid, attname, atttext)
  implicit none
      
  include 'netcdf.inc'

  integer nfid, varid
  character*(*) attname, atttext

  integer ret

  ret =      nf_get_att_text (nfid, varid, attname, atttext)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

subroutine wrap_copy_att (nfid, varid, attname, nfido, varido)
  implicit none
      
  include 'netcdf.inc'

  integer nfid, varid, nfido, varido
  character*(*) attname

  integer ret

  ret =      nf_copy_att (nfid, varid, attname, nfido, varido)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine

