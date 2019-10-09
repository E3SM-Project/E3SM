module fill_positions

use shr_kind_mod, only: r8 => shr_kind_r8
use control,      only: verbose
   
implicit none
private

include 'netcdf.inc'

public :: varspec_t, fillvar

type varspec_t
   character*(nf_max_name) :: name
   character*8 :: vshape

   integer :: varid
   integer :: totsiz
   integer :: xtype
   integer :: nx, ny, nz
   integer :: tpos
   integer :: count(nf_max_var_dims)

   integer :: x_vid
   integer :: y_vid
   integer :: z_vid
end type varspec_t

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

subroutine fillvar(ncid, name, xtype, vshape, dimnames, &
                   varid, nvdims, vardids, varspec)

   ! arguments

   integer,                 intent(in) :: ncid
   character*(*),           intent(in) :: name
   integer,                 intent(in) :: xtype
   character*8,             intent(in) :: vshape
   character*(nf_max_name), intent(in) :: dimnames(nvdims)
   integer,                 intent(in) :: varid
   integer,                 intent(in) :: nvdims
   integer,                 intent(in) :: vardids(nvdims)

   type(varspec_t),         intent(out) :: varspec

   ! Local workspace

   integer :: nx, ny, nz
   integer :: dimlen
   integer :: n
   integer :: unlimdimid
   character(len=*), parameter :: sub = 'fillvar: '
   !-------------------------------------------------------------------------------

   varspec%name   = name
   varspec%xtype  = xtype
   varspec%vshape = vshape
   varspec%varid  = varid
 
   varspec%totsiz   = 1
   varspec%tpos     = 1
   varspec%count(:) = 1

   nx    = 1 
   ny    = 1 
   nz    = 1

   ! Get dimension lengths.

   ! Identify the dimension ID of the unlimited dimension if there is one.
   if (nf_inq_unlimdim(ncid, unlimdimid) /= nf_noerr) then
      write(6,*) sub//'ERROR return from nf_inq_unlimdim'
   end if

   do n = 1, nvdims

      ! If the dimension ID is the unlimited one, assume that it's time
      ! and set the count to 1 since will work with 1 record at a time.
      if (vardids(n) == unlimdimid) then
         varspec%tpos     = n
         varspec%count(n) = 1
      else

         ! If the dimension ID isn't the unlimited one, get the corresponding size and
         ! increment the total size of the variable.
         if (nf_inq_dimlen(ncid, vardids(n), dimlen) == nf_noerr) then
            varspec%count(n) = dimlen
            varspec%totsiz   = varspec%totsiz * dimlen
         else
            write(6,*) sub//'ERROR return from nf_inq_dimlen for dimension '//trim(dimnames(n))
            stop
         end if
      end if
   end do

   ! Now get the IDs of the variables that contain the coordinate information for
   ! each spatial dimension

   if (vshape == 'xy') then

      varspec%nx = varspec%count(1)
      varspec%ny = varspec%count(2)
      varspec%nz = 1

      call fill_xpos(ncid,  dimnames(1), varspec%x_vid)
      call fill_yzpos(ncid, dimnames(2), varspec%y_vid)

   else if (vshape == 'xyz') then

      varspec%nx = varspec%count(1)
      varspec%ny = varspec%count(2)
      varspec%nz = varspec%count(3)

      call fill_xpos(ncid,  dimnames(1), varspec%x_vid)
      call fill_yzpos(ncid, dimnames(2), varspec%y_vid)
      call fill_yzpos(ncid, dimnames(3), varspec%z_vid)

   else if (vshape == 'xzy') then

      varspec%nx = varspec%count(1)
      varspec%ny = varspec%count(3)
      varspec%nz = varspec%count(2)

      call fill_xpos(ncid,  dimnames(1), varspec%x_vid)
      call fill_yzpos(ncid, dimnames(3), varspec%y_vid)
      call fill_yzpos(ncid, dimnames(2), varspec%z_vid)

   else if (vshape == 'nz') then

      varspec%nx = varspec%count(1)
      varspec%ny = 1
      varspec%nz = varspec%count(2)

      call fill_xpos(ncid,  dimnames(1), varspec%x_vid)
      call fill_yzpos(ncid, dimnames(1), varspec%y_vid)
      call fill_yzpos(ncid, dimnames(2), varspec%z_vid)

   else if (vshape == 'n') then

      varspec%nx = varspec%count(1)
      varspec%ny = 1
      varspec%nz = 1

      call fill_xpos(ncid,  dimnames(1), varspec%x_vid)
      call fill_yzpos(ncid, dimnames(1), varspec%y_vid)

   end if

   if (verbose) then
      if (nvdims == 1) then
         print*,'fillvar: ',trim(name),'(',trim(dimnames(1)),'=',varspec%count(1),')'
      else if (nvdims == 2) then
         print*,'fillvar: ',trim(name),'(',trim(dimnames(1)),'=',varspec%count(1),trim(dimnames(2)),'=',varspec%count(2),')'
      else if (nvdims == 3) then
         print*,'fillvar: ',trim(name),'(',trim(dimnames(1)),'=',varspec%count(1),trim(dimnames(2)),'=',varspec%count(2),&
                                           trim(dimnames(3)),'=',varspec%count(3),')'
      else if (nvdims == 4) then
         print*,'fillvar: ',trim(name),'(',trim(dimnames(1)),'=',varspec%count(1),trim(dimnames(2)),'=',varspec%count(2),&
                                           trim(dimnames(3)),'=',varspec%count(3),trim(dimnames(4)),'=',varspec%count(4),')'
      end if
   end if

end subroutine fillvar

!-------------------------------------------------------------------------------

subroutine fill_xpos(ncid, dimname, x_vid)

   ! arguments
   integer,                    intent(in) :: ncid
   character(len=nf_max_name), intent(in) :: dimname

   integer,                    intent(out) :: x_vid
   !-------------------------------------------------------------

   if (dimname == 'ncol') then
      call wrap_inq_varid(ncid, 'lon', x_vid)
   else
      call wrap_inq_varid(ncid, dimname, x_vid)
   end if

end subroutine fill_xpos

!-------------------------------------------------------------------------------

subroutine fill_yzpos(ncid, dimname, yz_vid)

   ! arguments
   integer,                    intent(in) :: ncid
   character(len=nf_max_name), intent(in) :: dimname

   integer,                    intent(out) :: yz_vid
   !-------------------------------------------------------------

   if (dimname == 'ncol') then
      call wrap_inq_varid(ncid, 'lat', yz_vid)
   else
      call wrap_inq_varid(ncid, dimname, yz_vid)
   end if

end subroutine fill_yzpos

!-------------------------------------------------------------------------------
end module fill_positions

