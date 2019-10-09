!===============================================================================
!
! Wrapper functions for netcdf.  Print message and abort upon failure return.
!
!===============================================================================

subroutine wrap_nf_open (file, mode, nfid)
   implicit none
   include 'netcdf.inc'

   character(len=*) file
   integer mode, nfid

   integer ret

   ret = nf_open (file, mode, nfid)
   if (ret /= nf_noerr) then
      write(6,*)nf_strerror(ret)
      write(6,*)'Unable to open file ', file
      call abort
   end if
end subroutine wrap_nf_open

subroutine wrap_nf_close (nfid)
   implicit none
   include 'netcdf.inc'

   integer nfid

   integer ret

   ret = nf_close (nfid)
   if (ret /= nf_noerr) then
      write(6,*)nf_strerror(ret)
      write(6,*)'Unable to close nfid ', nfid
      call abort
   end if
end subroutine wrap_nf_close

subroutine wrap_nf_inq_varid (nfid, varname, varid)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character*(*) varname

  integer ret

  ret = nf_inq_varid (nfid, varname, varid)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_inq_varid

subroutine wrap_nf_inq_dimlen (nfid, dimid, dimlen)
  implicit none
  include 'netcdf.inc'

  integer nfid, dimid, dimlen

  integer ret

  ret = nf_inq_dimlen (nfid, dimid, dimlen)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_inq_dimlen

subroutine wrap_nf_inq_dimid (nfid, dimname, dimid)
  implicit none
  include 'netcdf.inc'

  integer nfid, dimid
  character*(*) dimname

  integer ret

  ret = nf_inq_dimid (nfid, dimname, dimid)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_inq_dimid

subroutine wrap_nf_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, xtype, ndims, dimids(nf_max_dims), natts
  character*(*) varname

  integer ret

  ret = nf_inq_var (nfid, varid, varname, xtype, ndims, dimids, natts)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_inq_var

subroutine wrap_nf_def_dim (nfid, dimname, len, dimid)
  implicit none
  include 'netcdf.inc'

  integer nfid, len, dimid
  character*(*) dimname

  integer ret

  ret = nf_def_dim (nfid, dimname, len, dimid)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_def_dim

subroutine wrap_nf_def_var (nfid, name, xtype, nvdims, vdims, varid)
   implicit none
   include 'netcdf.inc'

   integer, intent(in):: nfid
   integer, intent(in)::xtype
   integer, intent(in)::nvdims
   integer, intent(out)::varid
   integer, intent(in):: vdims(nvdims+1)
   character*(*), intent(in):: name

   integer ret      ! NetCDF return code

   ret = nf_def_var (nfid, name, xtype, nvdims, vdims, varid)
   if (ret/=NF_NOERR) then
      write(6,*) 'wrap_nf_def_var failed for ',name
      write(6,*) nf_strerror (ret)
      call abort
   end if
end subroutine wrap_nf_def_var

subroutine wrap_nf_get_var_double (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  double precision arr(*)

  integer ret

  ret = nf_get_var_double (nfid, varid, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_get_var_double

subroutine wrap_nf_get_var_int (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer arr(*)

  integer ret

  ret = nf_get_var_int (nfid, varid, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_get_var_int

subroutine wrap_nf_put_var_double (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  double precision arr(*)

  integer ret
  ret = nf_put_var_double (nfid, varid, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_put_var_double

subroutine wrap_nf_put_var_int (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer arr(*)

  integer ret
  ret = nf_put_var_int (nfid, varid, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_put_var_int

subroutine wrap_nf_put_vara_int (nfid, varid, start, count, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(*), count(*)
  integer arr(*)

  integer ret

  ret = nf_put_vara_int (nfid, varid, start, count, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_put_vara_int

subroutine wrap_nf_get_vara_int (nfid, varid, start, count, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(*), count(*)
  integer arr(*)

  integer ret

  ret = nf_get_vara_int (nfid, varid, start, count, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_get_vara_int

subroutine wrap_nf_get_vara_double (nfid, varid, start, count, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, start(*), count(*)
  double precision arr(*)

  integer ret

  ret = nf_get_vara_double (nfid, varid, start, count, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_get_vara_double

subroutine wrap_nf_put_vara_double (nfid, varid, start, count, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), count(*)
  double precision arr(*)

  integer ret
  ret = nf_put_vara_double (nfid, varid, start, count, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_put_vara_double

subroutine wrap_nf_put_vara_real (nfid, varid, start, count, arr)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  integer start(*), count(*)
  real arr(*)

  integer ret
  ret = nf_put_vara_real (nfid, varid, start, count, arr)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_put_vara_real

subroutine wrap_nf_get_att_text (nfid, varid, name, text)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character*(*) name, text

  integer ret

  ret = nf_get_att_text (nfid, varid, name, text)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret), ": ", trim(name)
     call abort
  end if
end subroutine wrap_nf_get_att_text

subroutine wrap_nf_put_att_text (nfid, varid, name, len, text)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, len
  character*(*) name, text

  integer ret

  ret = nf_put_att_text (nfid, varid, name, len, text)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_put_att_text

subroutine wrap_nf_get_att_double (nfid, varid, name, dvals)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid
  character*(*) name
  double precision dvals

  integer ret

  ret = nf_get_att_double (nfid, varid, name, dvals)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret), ": ", trim(name)
     call abort
  end if
end subroutine wrap_nf_get_att_double

subroutine wrap_nf_put_att_double (nfid, varid, name, xtype, len, dvals)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, xtype, len
  character*(*) name
  double precision dvals

  integer ret

  ret = nf_put_att_double (nfid, varid, name, xtype, len, dvals)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_put_att_double

subroutine wrap_nf_put_att_real (nfid, varid, name, xtype, len, dvals)
  implicit none
  include 'netcdf.inc'

  integer nfid, varid, xtype, len
  character*(*) name
  real dvals

  integer ret

  ret = nf_put_att_real (nfid, varid, name, xtype, len, dvals)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_put_att_real

subroutine wrap_nf_copy_att (nfid_in, varid_in, name, nfid_out, varid_out)
  implicit none
  include 'netcdf.inc'

  integer nfid_in, varid_in, nfid_out, varid_out
  character*(*) name

  integer ret

  ret = nf_copy_att (nfid_in, varid_in, name, nfid_out, varid_out)
  if (ret /= NF_NOERR) then
     write(6,*) nf_strerror (ret)
     call abort
  end if
end subroutine wrap_nf_copy_att

