!===============================================================================
!
! Wrapper functions for netcdf.  Print message and abort upon failure return.
! WARNING: If enabling 32-bit internal reals, this code makes the assumption that
! selected_real_kind(6) (see prec.f90) will return "4", i.e. the number of bytes 
! in a 32-bit real number
!
!===============================================================================

   subroutine wrap_nf_create (path, omode, ncid)
      implicit none
      include 'netcdf.inc'

      character*(*), intent(in):: path
      integer, intent(in):: omode
      integer, intent(out):: ncid

      integer ret      ! NetCDF return code

      ret = nf_create (path, omode, ncid)
      if (ret /= NF_NOERR) then
         write(6,*)'WRAP_NF_CREATE: nf_create failed for file ',path
         call handle_error (ret)
      end if
   end subroutine wrap_nf_create

!===============================================================================

   subroutine wrap_nf_open (path, omode, ncid)
      implicit none
      include 'netcdf.inc'

      character*(*), intent(in):: path
      integer, intent(in):: omode
      integer, intent(out):: ncid

      integer ret      ! NetCDF return code

      ret = nf_open (path, omode, ncid)
      if (ret /= NF_NOERR) then
         write(6,*)'WRAP_NF_OPEN: nf_open failed for file ',path
         call handle_error (ret)
      end if
   end subroutine wrap_nf_open

   subroutine wrap_nf_put_att_text (nfid, varid, attname, atttext)
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
      if (ret/=NF_NOERR) then
         write(6,*)'wrap_nf_put_att_text varid, attname, text=', varid, attname, ' ', atttext
         call handle_error (ret)
      end if
   end subroutine wrap_nf_put_att_text

!===============================================================================

   subroutine wrap_nf_def_dim (nfid, dimname, len, dimid)
      implicit none
      include 'netcdf.inc'

      integer, intent(in):: nfid
      integer, intent(in):: len
      integer, intent(out):: dimid
      character*(*), intent(in):: dimname
      
      integer ret      ! NetCDF return code

      ret = nf_def_dim (nfid, dimname, len, dimid)
      if (ret/=NF_NOERR) then
         write(6,*)'wrap_nf_def_dim: nfid, dimname, len=', nfid, dimname, len
         call handle_error (ret)
      end if
   end subroutine wrap_nf_def_dim

!===============================================================================

   subroutine wrap_nf_inq_dimid (nfid, name, dimid)
      implicit none
      include 'netcdf.inc'

      integer, intent(in)  :: nfid
      integer, intent(out) :: dimid
      character(len=*), intent(in) :: name
   
      integer ret      ! NetCDF return code

      ret = nf_inq_dimid (nfid, name, dimid)
      if (ret /= NF_NOERR) then
         write(6,*)'wrap_nf_inq_dimid: nfid,name=', nfid, name
         call handle_error (ret)
      end if
   end subroutine wrap_nf_inq_dimid

!===============================================================================

   subroutine wrap_nf_inq_dimlen (nfid, dimid, dimlen)
      implicit none
      include 'netcdf.inc'

      integer, intent(in)  :: nfid
      integer, intent(in)  :: dimid
      integer, intent(out) :: dimlen
   
      integer ret      ! NetCDF return code

      ret = nf_inq_dimlen (nfid, dimid, dimlen)
      if (ret /= NF_NOERR) then
         write(6,*)'wrap_nf_inq_dimlen: nfid,dimid=', nfid, dimid
         call handle_error (ret)
      end if
   end subroutine wrap_nf_inq_dimlen

!===============================================================================

   subroutine wrap_nf_inq_varid (nfid, name, varid)
      implicit none
      include 'netcdf.inc'

      integer, intent(in)  :: nfid
      integer, intent(out) :: varid
      character(len=*), intent(in) :: name
   
      integer ret      ! NetCDF return code

      ret = nf_inq_varid (nfid, name, varid)
      if (ret /= NF_NOERR) then
         write(6,*)'wrap_nf_inq_varid: nfid,name=', nfid, name
         call handle_error (ret)
      end if
   end subroutine wrap_nf_inq_varid

!===============================================================================

   subroutine wrap_nf_def_var (nfid, name, xtype, nvdims, vdims, varid)
      implicit none
      include 'netcdf.inc'

      integer, intent(in):: nfid
      integer, intent(in)::xtype
      integer, intent(in)::nvdims
      integer, intent(out)::varid
      integer, intent(in):: vdims(nvdims)
      character*(*), intent(in):: name
   
      integer ret      ! NetCDF return code

      ret = nf_def_var (nfid, name, xtype, nvdims, vdims, varid)
      if (ret/=NF_NOERR) then
         write(6,*) 'wrap_nf_def_var: nfid, varname, nvdims, vdims=', name, nvdims, vdims
         call handle_error (ret)
      end if
   end subroutine wrap_nf_def_var

!===============================================================================

   subroutine wrap_nf_get_vara_int (nfid, varid, start, count, arr)
      use prec
      implicit none

      include 'netcdf.inc'

      integer, intent(in):: nfid
      integer, intent(in)::varid
      integer, intent(in)::start(*)
      integer, intent(in)::count(*)
      integer, intent(out):: arr(*)

      integer ret      ! NetCDF return code

      ret = nf_get_vara_int (nfid, varid, start, count, arr)
      if (ret/=NF_NOERR) then
         write(6,*)'WRAP_NF_GET_VARA_INT: error reading varid =', varid
         call handle_error (ret)
      end if
   end subroutine wrap_nf_get_vara_int

!===============================================================================

   subroutine wrap_nf_get_var_int (nfid, varid, arr)
      implicit none
      include 'netcdf.inc'

      integer, intent(in):: nfid
      integer, intent(in):: varid
      integer, intent(out):: arr(*)

      integer ret      ! NetCDF return code

      ret = nf_get_var_int (nfid, varid, arr)
      if (ret/=NF_NOERR) then
         write(6,*)'WRAP_NF_GET_VAR_INT: error reading varid =', varid
         call handle_error (ret)
      end if
   end subroutine wrap_nf_get_var_int

!===============================================================================

   subroutine wrap_nf_get_vara_double (nfid, varid, start, count, arr)
      use prec
      implicit none

      include 'netcdf.inc'

      integer, intent(in):: nfid
      integer, intent(in)::varid
      integer, intent(in)::start(*)
      integer, intent(in)::count(*)
      real(r8), intent(out):: arr(*)

      integer ret      ! NetCDF return code

      if (selected_real_kind(6) /= 4) then
         write(6,*)'WRAP_NF_GET_VARA_DOUBLE: cannot determine r4 vs. r8'
         call abort ()
      end if
         
      if (r8 == 4) then
         ret = nf_get_vara_real (nfid, varid, start, count, arr)
      else
         ret = nf_get_vara_double (nfid, varid, start, count, arr)
      end if
      if (ret/=NF_NOERR) then
         write(6,*)'WRAP_NF_GET_VARA_DOUBLE: error reading varid =', varid
         call handle_error (ret)
      end if
   end subroutine wrap_nf_get_vara_double

!===============================================================================

   subroutine wrap_nf_get_var_double (nfid, varid, arr)
      use prec
      implicit none

      include 'netcdf.inc'

      integer, intent(in):: nfid
      integer, intent(in)::varid
      real(r8), intent(out):: arr(*)

      integer ret      ! NetCDF return code

      if (selected_real_kind(6) /= 4) then
         write(6,*)'WRAP_NF_GET_VAR_DOUBLE: cannot determine r4 vs. r8'
         call abort ()
      end if
         
      if (r8 == 4) then
         ret = nf_get_var_real (nfid, varid, arr)
      else
         ret = nf_get_var_double (nfid, varid, arr)
      end if
      if (ret/=NF_NOERR) then
         write(6,*)'WRAP_NF_GET_VAR_DOUBLE: error reading varid =', varid
         call handle_error (ret)
      end if
   end subroutine wrap_nf_get_var_double

!===============================================================================

   subroutine wrap_nf_close (ncid)
      implicit none
      include 'netcdf.inc'

      integer, intent(in):: ncid

      integer ret      ! NetCDF return code
   
      ret = nf_close (ncid)
      if (ret/=NF_NOERR) then
         write(6,*)'WRAP_NF_CLOSE: nf_close failed for id ',ncid
         call handle_error (ret)
      end if
   end subroutine wrap_nf_close

!===============================================================================

   subroutine wrap_nf_put_var_int (nfid, varid, arr)
      implicit none
      include 'netcdf.inc'

      integer, intent(in):: nfid
      integer, intent(in):: varid
      integer, intent(in):: arr(*)

      integer ret      ! NetCDF return code

      ret = nf_put_var_int (nfid, varid, arr)
      if (ret /= NF_NOERR) call handle_error (ret)
   end subroutine wrap_nf_put_var_int

!===============================================================================

   subroutine wrap_nf_put_vara_int (nfid, varid, start, count, arr)
      implicit none
      include 'netcdf.inc'

      integer, intent(in):: nfid
      integer, intent(in):: varid
      integer, intent(in):: start(*)
      integer, intent(in):: count(*)
      integer, intent(in):: arr(*)

      integer ret      ! NetCDF return code

      ret = nf_put_vara_int (nfid, varid, start, count, arr)
      if (ret /= NF_NOERR) call handle_error (ret)
   end subroutine wrap_nf_put_vara_int

!===============================================================================

   subroutine wrap_nf_put_var_double (nfid, varid, arr)
      use prec
      implicit none
      include 'netcdf.inc'
   
      integer, intent(in):: nfid
      integer, intent(in):: varid
      real(r8), intent(in):: arr(*)

      integer ret      ! NetCDF return code

      if (selected_real_kind(6) /= 4) then
         write(6,*)'WRAP_NF_PUT_VARA_DOUBLE: cannot determine r4 vs. r8'
         call abort ()
      end if
         
      if (r8 == 4) then
         ret = nf_put_var_real (nfid, varid, arr)
      else
         ret = nf_put_var_double (nfid, varid, arr)
      end if
      if (ret /= NF_NOERR) call handle_error (ret)
   end subroutine wrap_nf_put_var_double

!===============================================================================

   subroutine wrap_nf_put_vara_double (nfid, varid, start, count, arr)
      use prec
      implicit none
      include 'netcdf.inc'
   
      integer, intent(in):: nfid
      integer, intent(in):: varid
      integer, intent(in):: start(*)
      integer, intent(in):: count(*)
      real(r8), intent(in):: arr(*)

      integer ret      ! NetCDF return code

      if (selected_real_kind(6) /= 4) then
         write(6,*)'WRAP_NF_PUT_VARA_DOUBLE: cannot determine r4 vs. r8'
         call abort ()
      end if
         
      if (r8 == 4) then
         ret = nf_put_vara_real (nfid, varid, start, count, arr)
      else
         ret = nf_put_vara_double (nfid, varid, start, count, arr)
      end if
      if (ret /= NF_NOERR) call handle_error (ret)
   end subroutine wrap_nf_put_vara_double

!===============================================================================

   subroutine handle_error(ret)
      implicit none
      include 'netcdf.inc'
   
      integer, intent(in):: ret
      
      write(6,*)nf_strerror(ret)
      call abort ()
   end subroutine handle_error

!===============================================================================
