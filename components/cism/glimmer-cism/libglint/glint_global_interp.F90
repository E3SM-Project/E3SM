!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_global_interp.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!TODO - Is this module still needed?  It is not currently used by other Glint modules.

module glint_global_interp

  use glint_global_grid
  use glimmer_global, only: dp
  use glimmer_physcon, only: pi
  implicit none

contains

  subroutine global_interp (in_grid,a,out_grid,ao,in_mask,out_mask,missing,error)

    ! This subroutine does an area weighted average from one grid,
    ! on a spherical earth, to another.  Logical masks may be assigned
    ! for each grid, and only those grid boxes which are masked true
    ! on both grids will be used.  A value of amm will be assigned
    ! to all nodes of the new grid which are initially false or have
    ! no data from the old grid. The new mask will also be changed to
    ! false where no data is available.
    !
    ! Restrictions:  longitude must be the first dimension and it
    !                be monotonically increasing (west to east).
    !
    !                latitude must be the second dimension and it
    !                must be monotonic.
    !
    !                values for longitude and latitude must be in
    !                degrees.
    !
    !                arrays that wrap around must repeat longitudes
    !                with a 360 degree increment.  it will be assumed
    !                that values in the wrapped input and mask arrays
    !                will also repeat (wrapped values in these arrays
    !                will not be used).
    !
    ! input
    !
    ! integer   idl    first dimension of input a and mask.
    ! integer   il     number of grid boxes in longitude for a and mask.
    ! real      alon   longitude (deg) limits of grid boxes for a and mask.
    ! integer   jl     number of grid boxes in latitude for a and mask.
    ! real      alat   latitude (deg) limits of grid boxes for a and mask.
    ! real      a      array of input data.
    ! logical   mask   mask for input data (.false. to mask out data).
    !
    ! output
    !
    ! integer   idlo   first dimension of output ao and masko.
    ! integer   ilo    number of grid boxes in longitude for ao and masko.
    ! real      alono  longitude (deg) limits of grid boxes for ao and masko.
    ! integer   jlo    number of grid boxes in latitude for ao and masko.
    ! real      alato  latitude (deg) limits of grid boxes for ao and masko.
    ! real      ao     array of output data.
    ! logical   masko  mask for output data (.false. to mask out data).
    ! integer   ier    error indication:
    !                  (values may be summed for multiple errors)
    !                  0  no errors
    !                  1  input longitude dimension and/or length <=0.
    !                  2  output dimension and/or length <=0.
    !                  4  input latititude dimension <=0.
    !                  8  output latitude dimension <=0.
    !                 16  wrap-around on input longitude grid doesn't
    !                     repeat (+360).
    !                 32  wrap-around on output longitude grid doesn't
    !                     repeat (+360).
    !                 64  longitude of input is not monotonic increasing.
    !                128  longitude of output is not monotonic increasing.
    !                256  latitude of input is not monotonic.
    !                512  latitude of output is not monotonic.
    !               1024  input longitude wraps but doesn't repeat identically.
    !               2048  output longitude wraps but doesn't repeat identically.
    !                 -1  output mask is changed.
    !                 -2  output mask contains all false values.

    ! --------------------------------------------------------
    ! Subroutine arguments
    ! --------------------------------------------------------

    type(global_grid)                              :: in_grid
    real(dp),dimension(:,:),         intent(in)    :: a
    type(global_grid)                              :: out_grid
    real(dp),dimension(:,:),         intent(out)   :: ao
    logical, dimension(:,:),optional,intent(inout) :: in_mask
    logical, dimension(:,:),optional,intent(inout) :: out_mask
    real(dp),               optional,intent(in)    :: missing
    integer,                optional,intent(out)   :: error

    ! --------------------------------------------------------
    ! Automatic arrays
    ! --------------------------------------------------------

    logical, dimension(size(a,1) ,size(a,2))  :: mask
    logical, dimension(size(ao,1),size(ao,2)) :: masko
    real(dp),dimension(size(a,1)+1)  :: alon
    real(dp),dimension(size(a,2)+1)  :: alat
    real(dp),dimension(size(ao,1)+1) :: alono
    real(dp),dimension(size(ao,2)+1) :: alato

    ! --------------------------------------------------------
    ! Internal variables
    ! --------------------------------------------------------

    integer  :: idl,il,jl,idlo,ilo,jlo,ier
    real(dp) :: amm,almx,almn,sgn,al,dln,almxo,almno,dlno,amnlto
    real(dp) :: amxlto,amnlt,amxlt,amnlno,amxlno,amnln,amxln,wt,avg
    real(dp) :: slatmx,wlat,slatmn,slon,slonp,slonmx,slonmn,delon
    integer  :: i,j,iil,iilo,j1,j2,jj,i1,i2,k,ii,iii,iip

    ! --------------------------------------------------------
    ! Set up array sizes and check things match up.
    ! --------------------------------------------------------

    idl=size(a,1)
    il=size(alon)-1
    jl=size(alat)-1
    idlo=size(ao,1)
    ilo=size(alono)-1
    jlo=size(alato)-1

    alon=in_grid%lon_bound
    alat=in_grid%lat_bound
    alono=out_grid%lon_bound
    alato=out_grid%lat_bound

    ! Check array sizes --------------------------------------

    if (idl/=in_grid%nx.or. &
         jl/=in_grid%ny.or. &
         idlo/=out_grid%nx.or. &
         jlo/=out_grid%ny) then
       print*,'Array size mismatch in global_interp'
       stop
    end if

    ! Deal with optional mask input --------------------------

    if (present(in_mask)) then
       mask=in_mask
    else
       mask=.true.
    endif

    if (present(out_mask)) then
       masko=out_mask
    else
       masko=.true.
    endif

    ! Set up missing value -----------------------------------

    if (present(missing)) then
       amm=missing
    else
       amm=50.
    end if
    ier=0

    ! Check that the sizes of the arrays given are sensible --

    if (idl < il.or.il <= 0) ier=1
    if (idlo < ilo.or.ilo <= 0) ier=ier+2
    if (jl <= 0)  ier=ier+4
    if (jlo <= 0) ier=ier+8
    if (ier > 0) then
       if (present(error)) error=ier
       return
    end if

    ! Check monotonic increasing input longitudes ------------

    do i=2,il
       if (alon(i) <= alon(i-1)) then
          ier=ier+64
          exit
       endif
    end do

    ! Check monotonic increasing output longitudes -----------

    do i=2,ilo 
       if (alono(i) <= alono(i-1)) then
          ier=ier+128
          exit
       endif
    end do

    ! Check monotonicity of input latitudes ------------------

    sgn=(alat(2)-alat(1))
    do j=2,jl
       if (sgn < 0.0) then
          if (alat(j)-alat(j-1) >= 0) then
             ier=ier+256
             exit
          endif
       else if (sgn > 0.0) then
          if (alat(j)-alat(j-1) <= 0.0) then
             ier=ier+256
             exit
          endif
       else
          ier=ier+256
          exit
       endif
    end do

    ! Check monotonicity of output latitudes ------------------

    sgn=(alato(2)-alato(1))
    do j=2,jlo 
       if (sgn < 0.0) then
          if (alato(j)-alato(j-1) >= 0.0) then
             ier=ier+512
             exit
          endif
       else if (sgn > 0.0) then
          if (alato(j)-alato(j-1) <= 0.0) then
             ier=ier+512
             exit
          endif
       else
          ier=ier+512
          exit
       endif
    end do

    ! Find wrap around of input grid, if it exists ------------

    iil=il
    almx=alon(1)
    almn=alon(1)

    do i=2,il+1
       almx=max(almx,alon(i))
       almn=min(almn,alon(i))
       al=abs(alon(i)-alon(1))-360.0
       if (abs(al) <= 1.e-4) then
          iil=i-1
          exit
       else if (al > 0.0) then
          ier=ier+1024
          go to 12
       endif
    end do

    dln=0.0
    if (almn < 0.0) then
       dln=int(-almn/360.0+.001)*360.0
    else if (almn > 360.0) then
       dln=-int(almn/360.0+.001)*360.0
    endif
12  continue

    ! Find wrap around of output grid, if it exists -----------

    iilo=ilo
    almxo=alono(1)
    almno=alono(1)
    do i=2,ilo+1
       almxo=max(almxo,alono(i))
       almno=min(almno,alono(i))
       al=abs(alono(i)-alono(1))-360.0
       if (abs(al) <= 1.e-4) then
          iilo=i-1
          exit
       else if (al > 0.0) then
          ier=ier+2048
          go to 15
       endif
    end do

    dlno=0.0
    if (almno < 0.0) then
       dlno=int(-almno/360.0+.001)*360.0
    else if (almno > 360.0) then
       dlno=-int(almno/360.0+.001)*360.0
    endif
15  continue

    ! Test for errors.  return if any --------------------------

    if (ier /= 0) then
       if (present(error)) error=ier
       return
    end if

    ! The output grid needs to begin with or after the input grid.

    if (almno+dlno < almn+dln) dlno=dlno+360.0

    do j=1,jlo ! loop 200 - over output latitudes
       ! find index limits in latitude to cover the new grid.
       j1=jl+1
       j2=0
       amnlto=min(alato(j),alato(j+1))
       amxlto=max(alato(j),alato(j+1))

       ! search for index limits in j.

       do jj=1,jl
          amnlt=min(alat(jj),alat(jj+1))
          amxlt=max(alat(jj),alat(jj+1))
          ! find jj limits
          if (amxlt > amnlto.and.amnlt < amxlto) then
             j1=min(jj,j1)
             j2=max(jj,j2)
          endif
       end do

       ! if input grid doesn't at least partially cover the
       ! output grid box, no values will be assigned.  mask out
       ! all values for the latitude.

       if (j2 < j1) then
          do i=1,iilo 
             ao(i,j)=amm
             if (masko(i,j)) ier=-1
             masko(i,j)=.false.
          end do
          cycle
       endif

       do i=1,iilo ! loop 100

          ! no need to compute if it is masked out.
          if (.not.masko(i,j)) cycle

          ! find index limits in longitude to cover the new grid.
          i1=3*il+1
          i2=0
          amnlno=min(alono(i),alono(i+1))+dlno
          amxlno=max(alono(i),alono(i+1))+dlno

          ! search for index limits in i.
          ! because of wrap around it is necessary to
          ! look through the data twice.
          ! the output grid longitudes have been adjusted
          ! (using dlno) such that the first longitude in
          ! the output grid is greater than the first
          ! longitude on the input grid.

          do k=0,1
             do ii=1,iil 
                amnln=min(alon(ii),alon(ii+1))+dln+k*360.0
                amxln=max(alon(ii),alon(ii+1))+dln+k*360.0
                ! find ii limits
                if (amxln > amnlno.and.amnln < amxlno) then
                   i1=min(ii+k*il,i1)
                   i2=max(ii+k*il,i2)
                endif
             end do
          end do

          ! if input grid doesn't partially cover the output
          ! grid box, no values will be assigned.  mask out
          ! the grid box.

          if (i2 < i1) then
             ao(i,j)=amm
             if (masko(i,j)) ier=-1
             masko(i,j)=.false.
             cycle
          endif

          wt=0.0
          avg=0.0

          do jj=j1,j2
             slatmx=max(alat(jj),alat(jj+1))
             slatmn=min(alat(jj),alat(jj+1))
             wlat=max(sin(min(amxlto,slatmx)*pi/180.d0)-sin(max(amnlto,slatmn)*pi/180.d0),0.d0)
             if (wlat /= 0.0) then
                do iii=i1,i2
                   slon=dln
                   slonp=dln
                   if (iii > iil) then
                      slon=slon+360.
                      slonp=slonp+360.
                   endif
                   ii=mod(iii-1,iil)+1
                   iip=ii+1
                   if (mask(ii,jj)) then
                      slon=slon+alon(ii)
                      slonp=slonp+alon(iip)
                      slonmx=max(slon,slonp)
                      slonmn=min(slon,slonp)
                      delon=max(min(amxlno,slonmx)-max(amnlno,slonmn),0.d0)
                      wt=wt+wlat*delon
                      avg=avg+a(ii,jj)*wlat*delon
                   endif
                end do
             endif
          end do

          if (wt > 0.0) then
             ao(i,j)=avg/wt
          else
             ao(i,j)=amm
             if (masko(i,j)) ier=-1
             masko(i,j)=.false.
          endif
100       continue
       end do
200    continue
    end do

    ! Finish filling the output array from wrap-around.

    if (iilo < ilo) then
       do j=1,jlo 
          do i=iilo+1,ilo
             ao(i,j)=ao(i-iilo,j)
             masko(i,j)=masko(i-iilo,j)
          end do
       end do
    endif

    ! Check if output masko is all false.

    if (all(.not.masko)) ier=-2

    ! Copy outputs if necessary

    if (present(error)) error=ier
    if (present(in_mask)) in_mask=mask
    if (present(out_mask)) out_mask=masko

  end subroutine global_interp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine interp_error(err,line)

    use glimmer_log

    integer :: err,line
    character(50) :: message
    integer :: e

    write(message,'(A,I6)')'Interpolation errors at line ',line
    call write_log(message)

    e=err

    call err_check(e,2048,'output longitude wraps but doesn''t repeat identically')
    call err_check(e,1024,'input longitude wraps but doesn''t repeat identically')
    call err_check(e,512, 'latitude of output is not monotonic')
    call err_check(e,256, 'latitude of input is not monotonic')
    call err_check(e,128, 'longitude of output is not monotonic increasing')
    call err_check(e,64,  'longitude of input is not monotonic increasing')
    call err_check(e,32,  'wrap-around on output longitude grid doesn''t repeat (+360)')
    call err_check(e,16,  'wrap-around on input longitude grid doesn''t repeat (+360)')
    call err_check(e,8,   'output latitude dimension <=0')
    call err_check(e,4,   'input latitude dimension <=0')
    call err_check(e,2,   'output dimension and/or length <=0')
    call err_check(e,1,   'input longitude dimension and/or length <=0')
    stop

  end subroutine interp_error

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine err_check(e,en,out)

    use glimmer_log

    integer :: e,en
    character(*) :: out

    if (e>en) then
       call write_log(out)
       e=e-en
    end if

  end subroutine err_check

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glint_global_interp
