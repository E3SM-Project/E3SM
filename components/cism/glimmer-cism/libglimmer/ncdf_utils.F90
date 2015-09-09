!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   ncdf_utils.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

!> This code provides a simple interface to create and then add
!! to a netcdf file containing time-slices of a single 2D field,
!! for use in debugging.
module ncdf_utils
  
  use netcdf
  use glimmer_global, only: sp, dp

  implicit none

  type ncdf_utils_type
     integer :: id,varid,dimid1,dimid2,dimid3,d3id
     integer :: next=1
     character(100) :: fname
  end type ncdf_utils_type

contains

  !TODO - Change d1 and d2 to dp?
  !       Note: This subroutine currently is not called, as far as I can tell

  subroutine ncdf_utils_create(handle,fname,varname,d1name,d2name,d1,d2)

    type(ncdf_utils_type),intent(out) :: handle !< Netcdf file handles
    character(*),         intent(in) :: fname   !< File name
    character(*),         intent(in) :: varname !< Variable name
    character(*),         intent(in) :: d1name  !< Name of first dimension
    character(*),         intent(in) :: d2name  !< Name of second dimension
    real(sp),dimension(:),intent(in) :: d1      !< Dimension 1
    real(sp),dimension(:),intent(in) :: d2      !< Dimension 2

    integer :: ncerr,d1id,d2id

    ! Create file

    ncerr=nf90_create(fname,0,handle%id)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)
    handle%fname=fname

    ! Define dimensions

    ncerr=nf90_def_dim(handle%id,d1name,size(d1),handle%dimid1)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)
    ncerr=nf90_def_dim(handle%id,d2name,size(d2),handle%dimid2)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)
    ncerr=nf90_def_dim(handle%id,'time',NF90_UNLIMITED,handle%dimid3)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)

    ! Define dimension variables

    ncerr=nf90_def_var(handle%id,d1name,NF90_FLOAT,(/handle%dimid1/),d1id)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)
    ncerr=nf90_def_var(handle%id,d2name,NF90_FLOAT,(/handle%dimid2/),d2id)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)
    ncerr=nf90_def_var(handle%id,'time',NF90_FLOAT,(/handle%dimid3/),handle%d3id)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)

    ! Define 2D variable

    ncerr=nf90_def_var(handle%id,varname,NF90_DOUBLE, &
         (/handle%dimid1,handle%dimid2,handle%dimid3/),handle%varid)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)

    ! Exit define mode and save dimension variables

    ncerr=nf90_enddef(handle%id)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)
    ncerr=nf90_put_var(handle%id,d1id,d1)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)
    ncerr=nf90_put_var(handle%id,d2id,d2)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)

  end subroutine ncdf_utils_create

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine ncdf_utils_write(handle,var,time)

    type(ncdf_utils_type),  intent(inout) :: handle
    real(dp),dimension(:,:),intent(in)    :: var
    real(dp),               intent(in)    :: time

    integer :: ncerr

    ncerr=nf90_put_var(handle%id,handle%varid,real(var,dp),(/1,1,handle%next/))
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)
    ncerr=nf90_put_var(handle%id,handle%d3id,time,(/handle%next/))
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)
    ncerr=nf90_sync(handle%id)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)

    handle%next=handle%next+1

  end subroutine ncdf_utils_write

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine ncdf_utils_close(handle)

    type(ncdf_utils_type),   intent(in) :: handle

    integer :: ncerr

    ncerr=nf90_close(handle%id)    

  end subroutine ncdf_utils_close

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine ncdf_utils_read_slice(filename,varname,slice,array)

    character(*),           intent(in)  :: filename
    character(*),           intent(in)  :: varname
    integer,                intent(in)  :: slice
    real(dp),dimension(:,:),intent(out) :: array

    integer :: ncerr,fileid,varid

    ncerr=nf90_open(filename,0,fileid)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)
    ncerr=nf90_inq_varid(fileid,varname,varid)
    if (ncerr/=NF90_NOERR) call ncerr_handle(ncerr)
    ncerr=nf90_get_var(fileid,varid,array,(/1,1,slice/))

  end subroutine ncdf_utils_read_slice

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine ncerr_handle(ncerr)

    integer,intent(in) :: ncerr

    print*,nf90_strerror(ncerr)
    stop

  end subroutine ncerr_handle

end module ncdf_utils
