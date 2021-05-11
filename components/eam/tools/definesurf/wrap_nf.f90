subroutine wrap_inq_varid (nfid, varname, varid)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character*(*) varname

  integer ret

  ret = nf_inq_varid (nfid, varname, varid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_inq_varid

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

subroutine wrap_get_var8 (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  real*8 arr(*)

  integer ret

  ret = nf_get_var_double (nfid, varid, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_get_var8

subroutine wrap_put_var8 (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  real*8 arr(*)

  integer ret
  ret = nf_put_var_double (nfid, varid, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_put_var8

subroutine wrap_get_vara8 (nfid, varid, start, count, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(*), count(*)
  real*8 arr(*)

  integer ret

  ret = nf_get_vara_double (nfid, varid, start, count, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_get_vara8

subroutine wrap_put_vara8 (nfid, varid, start, count, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), count(*)
  real*8 arr(*)

  integer ret
  ret = nf_put_vara_double (nfid, varid, start, count, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_put_vara8

subroutine wrap_put_att_text (nfid, varid, attname, atttext)
  implicit none
  include 'netcdf.inc'

  integer, intent(in):: nfid
  integer, intent(in):: varid
  character*(*), intent(in):: attname
  character*(*), intent(in):: atttext
  
  integer ret      ! NetCDF return code
  integer siz

  siz = len_trim(atttext)
  ret = nf_put_att_text (nfid, varid, attname, siz, atttext)
  if (ret/=NF_NOERR) call handle_error (ret)
end subroutine wrap_put_att_text

subroutine wrap_put_att_double (nfid, varid, name, xtype, len, dvals)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, xtype, len
  character*(*) name
  real*8 dvals

  integer ret

  ret = nf_put_att_double (nfid, varid, name, xtype, len, dvals)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_put_att_double

