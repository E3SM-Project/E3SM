subroutine sm121 (a, plon, nlat, nlon)

  use shr_kind_mod, only: r8 => shr_kind_r8

!      
! perform 1-2-1 smoothing using data array a.  On reduced grid, linearly
! interpolate to a rectangular grid (nlon(j),3) before interpolating
!
!-----------------------------------------------------------------------
  implicit none
!-----------------------------Arguments---------------------------------

  integer plon             ! Input:  Lon dim
  integer nlat             ! Input:  Lat dim
  integer nlon(nlat)       ! Number of longitudes per latitude
  real(r8)  a(plon,nlat)     ! I/O:    Array to be smoothed

!--------------------------Local variables------------------------------

  integer i,j              ! Indices
  integer imin,imax        ! Indices
  integer jmax,jmin        ! Indices
!
! Dynamic
!
  real(r8) xin(plon,nlat)
  real(r8) xout(plon)
  real(r8) temp(plon,nlat)   ! Temp array
  real(r8) tempjmin(plon)    ! Temp array
  real(r8) tempjmax(plon)    ! Temp array
!
!-----------------------------------------------------------------------
!     
  temp(:,:) = a(:,:)
!
! first do the S and N boundaries.
!
  do i=1,nlon(1)
    imin = i - 1
    imax = i + 1
    if( imin .lt.   1 )    imin = imin + nlon(1)
    if( imax .gt. nlon(1)) imax = imax - nlon(1)
    a(i,1) =    (temp(imin,1)  + 3.*temp(i,1)   +temp(imax,1))/5.
  end do

  do i=1,nlon(nlat)
    imin = i - 1
    imax = i + 1
    if( imin .lt.   1 )       imin = imin + nlon(nlat)
    if( imax .gt. nlon(nlat)) imax = imax - nlon(nlat)
    a(i,nlat) = (temp(imin,nlat)+3.*temp(i,nlat)+temp(imax,nlat))/5.
  end do
!
! Define xin array for each latitude
!
  do j=1,nlat
    do i=1,nlon(j)
      xin(i,j) = (i-1)*360./nlon(j)
    end do
  end do
!
! Linearly interpolate data N and S of each target latitude to the longitudes
! of each target latitude before applying 1-2-1 filter
!
  do j=2,nlat-1
    jmin = j - 1
    jmax = j + 1
    xout(:) = xin(:,j)
    call lininterp (temp(1,jmin), nlon(jmin), 1, xin(1,jmin), &
                    tempjmin, nlon(j), 1, xout, .true.)
    call lininterp (temp(1,jmax), nlon(jmax), 1, xin(1,jmax), &
                    tempjmax, nlon(j), 1, xout, .true.)

    do i=1,nlon(j)
      imin = i - 1
      imax = i + 1
      if( imin .lt.   1 )    imin = imin + nlon(j)
      if( imax .gt. nlon(j)) imax = imax - nlon(j)
      a(i,j) =            (tempjmin(i)  + &
          temp(imin,j) + 4.*temp(i,j)  + temp(imax,j) + &
                           tempjmax(i)  ) / 8.
    enddo
  enddo
!
  return
end subroutine sm121
