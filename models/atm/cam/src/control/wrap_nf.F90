!-------------------------------------------------------------------------------
!
! WARNING: USE OF THIS MODULE WITHIN CAM IS DEPRECATED. ALL
!          HANDLING OF NETCDF FILES SHOULD ULTIMATELY BE DONE
!          BY PIO, OR FOR PORTABLE CODE, SHOULD CONTACT THE
!          NETCDF F90 INTERFACE DIRECTLY.
!
!          DO NOT USE THIS MODULE.
!
!-------------------------------------------------------------------------------

module wrap_nf
  use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use abortutils,   only: endrun
  use cam_logfile,  only: iulog
  use netcdf

!-------------------------------------------------------------------------------
!
! Purpose:
!
! Wrapper routines for the netCDF library for input and output data.
!
! Author: Jim Rosinski
!
! $Id$
!
!-------------------------------------------------------------------------------


contains
!===============================================================================

!===============================================================================

   subroutine wrap_redef (nfid)
   implicit none
   
   integer, intent(in):: nfid

   integer ret      ! NetCDF return code

   ret = nf90_redef (nfid)
   if (ret/=NF90_NOERR) call handle_error (ret)
   
   end subroutine wrap_redef
!===============================================================================

   subroutine wrap_enddef (nfid)
   implicit none
   
   
   integer, intent(in):: nfid

   integer ret      ! NetCDF return code

   ret = nf90_enddef (nfid)
   if (ret/=NF90_NOERR) call handle_error (ret)

   end subroutine wrap_enddef

   subroutine wrap_create (path, cmode, ncid)
   implicit none
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Create a netCDF file for reading and/or writing
!
!-------------------------------------------------------------------------------
   
   character*(*), intent(in):: path
   integer, intent(in):: cmode
   integer, intent(out):: ncid

   integer ret      ! NetCDF return code

   ret = nf90_create (path, cmode, ncid)
   if (ret/=NF90_NOERR) call handle_error (ret)

   end subroutine wrap_create

!===============================================================================

  subroutine wrap_inq_unlimdim (nfid, dimid)
    implicit none
    !-------------------------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Get dimid for the unlimited dimension.
    !
    !-------------------------------------------------------------------------------
    integer, intent(in)::  nfid
    integer, intent(out):: dimid

    integer ret      ! NetCDF return code

    ret = nf90_inquire(nfid, unlimitedDimId=dimid)
    if (ret/=NF90_NOERR) call handle_error (ret)
  end subroutine wrap_inq_unlimdim

   subroutine wrap_inq_dim (nfid, dimid, dimname, dimlen)
   implicit none
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets dimension name for a given dimension id
!
!-------------------------------------------------------------------------------
   integer, intent(in)::  nfid
   integer, intent(in):: dimid
   integer, intent(out):: dimlen
   character*(*), intent(out):: dimname

   integer ret      ! NetCDF return code

   ret = nf90_inquire_dimension (nfid, dimid, dimname, dimlen)
   if (ret/=NF90_NOERR) call handle_error (ret)

   end subroutine wrap_inq_dim

   subroutine wrap_inq_nvars (nfid, nvars)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets number of variables in file
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in)::  nfid
   integer, intent(out):: nvars

   integer ret      ! NetCDF return code

   ret = nf90_inquire (nfid, nvars)
   if (ret/=NF90_NOERR) call handle_error (ret)

   end subroutine wrap_inq_nvars

   subroutine wrap_inq_ndims (nfid, ndims)
   implicit none
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets number of dimensions in file
!
!-------------------------------------------------------------------------------
   integer, intent(in)::  nfid
   integer, intent(out):: ndims

   integer ret      ! NetCDF return code

   ret = nf90_inquire(nfid, ndims)
   if (ret/=NF90_NOERR) call handle_error (ret)

 end subroutine wrap_inq_ndims

!===============================================================================

   subroutine wrap_inq_dimid (nfid, dimname, dimid)
   implicit none
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the dimension id
!
!-------------------------------------------------------------------------------
   
   integer, intent(in):: nfid
   integer, intent(out):: dimid
   character*(*), intent(in):: dimname

   integer ret      ! NetCDF return code

   ret = nf90_inq_dimid (nfid, dimname, dimid)
   if(ret==NF90_NOERR) return
   if (ret/=NF90_EBADDIM) call handle_error (ret)
   dimid=-1  ! do not exist on bad dim.  This allows the user to check for dims that may not
             ! be in the file  
   end subroutine wrap_inq_dimid

!===============================================================================

   subroutine wrap_inq_dimlen (nfid, dimid, dimlen)
   implicit none
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the dimension length for a given dimension
!
!-------------------------------------------------------------------------------
   
   integer, intent(in)::  nfid
   integer, intent(in)::  dimid 
   integer, intent(out):: dimlen
   
   integer ret      ! NetCDF return code

   ret = nf90_inquire_dimension (nfid, dimid, len=dimlen)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_inq_dimlen

!===============================================================================

   subroutine wrap_inq_vardimid (nfid, varid, dimids)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the dimension Id's from a variable
!
!-------------------------------------------------------------------------------
   implicit none
   
   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(out)::  dimids(:)
   
   integer ret      ! NetCDF return code

   ret = nf90_inquire_variable (nfid, varid, dimids=dimids)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_inq_vardimid

!===============================================================================

   subroutine wrap_inq_varndims (nfid, varid, ndims)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the dimension Id's from a variable
!
!-------------------------------------------------------------------------------
   implicit none
   
   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(out):: ndims
   
   integer ret      ! NetCDF return code

   ret = nf90_inquire_variable (nfid, varid, ndims=ndims)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_inq_varndims

!===============================================================================

   subroutine wrap_inq_varid (nfid, varname, varid, abort)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the variable ID
!
!-------------------------------------------------------------------------------
   implicit none
   
   integer, intent(in):: nfid
   integer, intent(out):: varid
   character*(*), intent(in):: varname
   logical, optional :: abort
 
   integer ret      ! NetCDF return code
   logical :: call_endrun

   ret = nf90_inq_varid (nfid, varname, varid)
   if (ret/=NF90_NOERR ) then
      call_endrun = .true.
      if ( present(abort) ) then
         call_endrun = abort
      endif

      if ( call_endrun ) then
         write(iulog,*)'wrap_inq_varid: id for ',trim(varname),' not found'
         call handle_error (ret)
      else
         varid = -1
      endif
   end if
   end subroutine wrap_inq_varid

!===============================================================================

   subroutine wrap_inq_var (nfid, varid, varname, xtype, ndims, &
                            dimids, natts)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the variable name, type, number of dimensions, dimension ID's, and number of attributes
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(out):: xtype
   integer, intent(out):: ndims
   integer, intent(out):: dimids(:)
   integer, intent(out):: natts
   character*(*), intent(out):: varname

   integer ret      ! NetCDF return code

   ret = nf90_inquire_variable (nfid, varid, varname, xtype, ndims, dimids, &
                     natts)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_inq_var

!===============================================================================

   subroutine wrap_inq_varname (nfid, varid, varname)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the variable name from the dimension ID
!
!-------------------------------------------------------------------------------
   implicit none

   integer ret      ! NetCDF return code

   integer, intent(in):: nfid
   integer, intent(in):: varid
   character*(*), intent(out):: varname

   ret = nf90_inquire_variable(nfid, varid, varname)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_inq_varname

!===============================================================================

   subroutine wrap_get_att_text (nfid, varid, attname, atttext)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the attribute text from the given variable ID and attribute name
!
!-------------------------------------------------------------------------------
   implicit none
   
   integer, intent(in):: nfid
   integer, intent(in):: varid
   character*(*), intent(in):: attname
   character*(*), intent(out):: atttext

   integer ret      ! NetCDF return code

   ret = nf90_get_att(nfid, varid, attname, atttext)
   if (ret/=NF90_NOERR) then
      write(iulog,*)'WRAP_GET_ATT_TEXT: error reading attribute '//trim(attname)
      call handle_error (ret)
   endif
   end subroutine wrap_get_att_text

!===============================================================================

   subroutine wrap_put_att_text (nfid, varid, attname, atttext)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Puts the given attribute text to variable ID.
!
! This routine violates the convetion that the wrapper codes take an identical
! set of arguments as the netcdf library code.  The length of the character
! argument is computed inside the wrapper.
!
!-------------------------------------------------------------------------------
   implicit none
   
   integer, intent(in):: nfid
   integer, intent(in):: varid
   character*(*), intent(in):: attname
   character*(*), intent(in):: atttext

   integer ret      ! NetCDF return code
   integer siz

   ret = nf90_put_att(nfid, varid, attname, atttext)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_put_att_text

!===============================================================================

   subroutine wrap_put_att_realx (nfid, varid, attname, xtype, len, &
                                  attval)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Puts the given real attribute to the variable id
!
!-------------------------------------------------------------------------------
   implicit none
   
   integer , intent(in):: nfid
   integer , intent(in):: varid
   integer , intent(in):: xtype
   integer , intent(in):: len
   character*(*) , intent(in):: attname
   real(r8) , intent(in):: attval(len)

   integer ret      ! NetCDF return code

   ret = nf90_put_att(nfid, varid, attname, attval)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_put_att_realx
!===============================================================================

   subroutine wrap_def_dim (nfid, dimname, len, dimid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Defines the input dimension
!
!-------------------------------------------------------------------------------
   implicit none
   integer, intent(in):: nfid
   integer, intent(in):: len
   integer, intent(out):: dimid
   character*(*), intent(in):: dimname
   
   integer ret      ! NetCDF return code

   ret = nf90_def_dim (nfid, dimname, len, dimid)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_def_dim

!===============================================================================

   subroutine wrap_def_var (nfid, name, xtype, nvdims, vdims, varid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Defines the given variable
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in)::xtype
   integer, intent(in)::nvdims
   integer, intent(out)::varid
   integer, intent(in):: vdims(nvdims+1)
   character*(*), intent(in):: name
   
   integer ret      ! NetCDF return code

   ret = nf90_def_var(nfid, name, xtype, vdims, varid)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_def_var

!===============================================================================

   subroutine wrap_get_var_realx (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the given real variable from a input file
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   real(r8), intent(out):: arr(:)

   integer ret      ! NetCDF return code

   ret = nf90_get_var (nfid, varid, arr)
   if (ret/=NF90_NOERR) then
     write(iulog,*)'WRAP_GET_VAR_REALX: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_var_realx

!===============================================================================

   subroutine wrap_get_var_real4 (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the given real variable from a input file
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   real(r4), intent(out):: arr(:)

   integer ret      ! NetCDF return code

   ret = nf90_get_var (nfid, varid, arr)
   if (ret/=NF90_NOERR) then
     write(iulog,*)'WRAP_GET_VAR_REAL4: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_var_real4

!===============================================================================

   subroutine wrap_get_scalar_realx (nfid, varid, x)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the given real variable from a input file
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   real(r8), intent(out):: x

   integer ret      ! NetCDF return code

   ret = nf90_get_var (nfid, varid, x)
   if (ret/=NF90_NOERR) then
     write(iulog,*)'WRAP_GET_SCALAR_REALX: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_scalar_realx

!===============================================================================

   subroutine wrap_get_var_int (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the given integer variable from a input file
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(out):: arr(:)

   integer ret      ! NetCDF return code

   ret = nf90_get_var (nfid, varid, arr)
   if (ret/=NF90_NOERR) then
     write(iulog,*)'WRAP_GET_VAR_INT: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_var_int

!===============================================================================

   subroutine wrap_get_scalar_int (nfid, varid, x)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the given integer variable from a input file
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(out):: x

   integer ret      ! NetCDF return code

   ret = nf90_get_var (nfid, varid, x)
   if (ret/=NF90_NOERR) then
     write(iulog,*)'WRAP_GET_SCALAR_INT: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_scalar_int

!===============================================================================

   subroutine wrap_get_vara_realx (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets a range of the given real variable from a input file
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in)::varid
   integer, intent(in)::start(:)
   integer, intent(in)::count(:)
   real(r8), intent(out):: arr(:)

   integer ret      ! NetCDF return code

   ret = nf90_get_var (nfid, varid, arr, start, count)
   if (ret/=NF90_NOERR) then
     write(iulog,*)'WRAP_GET_VARA_REALX: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_vara_realx

!===============================================================================

   subroutine wrap_get_vara_int (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets a range of the given integer variable from a input file.
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(:)
   integer, intent(in):: count(:)
   integer, intent(out):: arr(:)

   integer ret      ! NetCDF return code

   ret = nf90_get_var (nfid, varid, arr, start, count)
   if (ret/=NF90_NOERR) then
     write(iulog,*)'WRAP_GET_VARA_INT: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_vara_int

!===============================================================================

   subroutine wrap_get_vara_text (nfid, varid, start, count, text)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets a range of the given text variable to input file.
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(:)
   integer, intent(in):: count(:)
   character(len=*), intent(out):: text(:)

   integer ret      ! NetCDF return code

   ret = nf90_get_var (nfid, varid, text, start, count)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_get_vara_text

!===============================================================================

   subroutine wrap_open (path, omode, ncid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Open a netCDF file
!
!-------------------------------------------------------------------------------
   implicit none

   character*(*), intent(in):: path
   integer, intent(in):: omode
   integer, intent(out):: ncid

   integer ret      ! NetCDF return code

   ret = nf90_open (path, omode, ncid)
   if (ret/=NF90_NOERR) then
     write(iulog,*)'WRAP_OPEN: nf90_open failed for file ',path
     call handle_error (ret)
   end if
   end subroutine wrap_open

!===============================================================================

   subroutine wrap_close (ncid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Close netCDF file
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: ncid

   integer ret      ! NetCDF return code

   ret = nf90_close (ncid)
   if (ret/=NF90_NOERR) then
     write(iulog,*)'WRAP_CLOSE: nf90_close failed for id ',ncid
     call handle_error (ret)
   end if
   end subroutine wrap_close

!===============================================================================

   subroutine wrap_put_var_int (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put a integer variable on output file.
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: arr(:)

   integer ret      ! NetCDF return code

   ret = nf90_put_var (nfid, varid, arr)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_put_var_int

!===============================================================================

   subroutine wrap_put_var1_int (nfid, varid, index, ival)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put a variable on output file at a given index.
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: index(:)
   integer, intent(in):: ival

   integer ret      ! NetCDF return code

   ret = nf90_put_var (nfid, varid, ival, index)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_put_var1_int

!===============================================================================

   subroutine wrap_put_vara_int (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put a range of a integer variable on a output file.
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(:)
   integer, intent(in):: count(:)
   integer, intent(in):: arr(:)

   integer ret      ! NetCDF return code

   ret = nf90_put_var (nfid, varid, arr, start, count)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_put_vara_int

!===============================================================================

   subroutine wrap_put_vara_text (nfid, varid, start, count, text)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put a range of the given text variable to output file.
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(:)
   integer, intent(in):: count(:)
   character(len=*), intent(in):: text(:)

   integer ret      ! NetCDF return code

   ret = nf90_put_var (nfid, varid, text, start, count)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_put_vara_text

!===============================================================================

   subroutine wrap_put_var1_realx (nfid, varid, index, val)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put the given real variable to output file at given index.
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: index(:)
   real(r8), intent(in):: val

   integer ret      ! NetCDF return code

   ret = nf90_put_var (nfid, varid, val, index)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_put_var1_realx

!===============================================================================

   subroutine wrap_put_vara_realx (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Output the given portion of the real array.
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(:)
   integer, intent(in):: count(:)
   real(r8), intent(in):: arr(:)

   integer ret      ! NetCDF return code
   ret = nf90_put_var (nfid, varid, arr, start, count)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_put_vara_realx

!===============================================================================

   subroutine wrap_put_vara_real (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Output the given portion of the real array.
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(:)
   integer, intent(in):: count(:)
   real(r4), intent(in):: arr(:)

   integer ret      ! NetCDF return code
   ret = nf90_put_var (nfid, varid, arr, start, count)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_put_vara_real

!===============================================================================

   subroutine wrap_put_var_realx (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put the given real variable to output file.
!
!-------------------------------------------------------------------------------
   implicit none

   integer, intent(in):: nfid
   integer, intent(in):: varid
   real(r8), intent(in):: arr(:)

   integer ret      ! NetCDF return code

   ret = nf90_put_var (nfid, varid, arr)
   if (ret/=NF90_NOERR) call handle_error (ret)
   end subroutine wrap_put_var_realx

!===============================================================================

   subroutine handle_error(ret)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Handle netCDF errors.
!
!-------------------------------------------------------------------------------

   implicit none

   integer, intent(in):: ret
   
   write(iulog,*)nf90_strerror(ret)
   call endrun ('HANDLE_ERROR')
   end subroutine handle_error

!===============================================================================
 end module wrap_nf
