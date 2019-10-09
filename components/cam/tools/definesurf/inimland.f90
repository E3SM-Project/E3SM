subroutine inimland (plon, nlat, nlon_reduced, mlatcnts, mloncnts, topofile, &
                     verbose, make_ross, landm_reduced) 
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
!
! Input arguments
!
  integer , intent(in) :: plon                ! number of longitudes
  integer , intent(in) :: nlat                ! number of latitudes
  integer , intent(in) :: nlon_reduced(nlat)  ! number of reduced latitudes
  real(r8), intent(in) :: mlatcnts(nlat)      ! latitude at center of grid cell
  real(r8), intent(in) :: mloncnts(plon,nlat) ! model cell ceneter longitudes
  character(len=*), intent(in) :: topofile    ! high res topo file
  logical, intent(in) :: verbose              ! verbose output  
  logical, intent(in) :: make_ross            ! flag to make Ross ice shelf
!
! Output arguments
!
  real(r8), intent(out) :: landm_reduced(plon,nlat) ! landm on reduced grid

! Local variables

  real(r8) landm(plon,nlat)       ! landm on full grid
  real(r8) clon(plon)
  real(r8) clon_reduced(plon,nlat)
  real(r8) cont(plon,nlat)
  real(r8) temp(plon,nlat)
  real(r8) dmax
  real(r8) arad
  real(r8) dist
  real(r8) sum
  real(r8) cs(nlat)
  real(r8) ss(nlat)
  real(r8) c1
  real(r8) s1
  real(r8) c2
  real(r8) s2
  real(r8) dx
  real(r8) dy
  real(r8) term
  real(r8) pi
  real(r8) sgh(plon,nlat)    ! required by SGHPHIS (unused locally)
  real(r8) phis(plon,nlat)   ! required by SGHPHIS (unused locally)
  real(r8) oro(plon,nlat)    ! land/ocean flag 
  real(r8) fland(plon,nlat)  ! land fraction output from SGHPHIS
  real(r8) mloncnts_full(plon,nlat) ! longitudes for rectangular grid

  integer i
  integer j
  integer ii
  integer jj
  integer iplm1
  integer jof
  integer iof
  integer itmp
  integer jmin, jmax
  integer nlon(nlat)
  integer latid

  pi = acos(-1.d0)
!
! Define longitudes for a rectangular grid: index nlat/2+1 will be a latitude
! closest to the equator, i.e. with the most points in a reduced grid.
!
  nlon(:) = plon
  do j=1,nlat
    mloncnts_full(:,j) = mloncnts(:,nlat/2+1)
  end do

  call sghphis (plon, nlat, nlon, mlatcnts, mloncnts_full, topofile, &
                verbose, sgh, phis, fland)
!
! Define land mask.  Set all non-land points to ocean (i.e. not sea ice).
!
  where (fland(:,:) >= 0.5)
    oro(:,:) = 1.
  elsewhere
    oro(:,:) = 0.
  endwhere
!
! Overwrite ORO flag as land for Ross ice shelf: note that the ORO field
! defined in this routine is only used locally.
!
  do j=1,nlat
    if (make_ross .and. mlatcnts(j) < -79.) then
      do i=1,plon
        oro(i,j) = 1.
      end do
    end if
  end do
!
! Code lifted directly from cldwat.F
!
  dmax = 2.e6 ! distance to carry the mask
  arad = 6.37e6
  do i = 1,plon
    clon(i) = 2.*(i-1)*pi/plon
  end do
!
! first isolate the contenents 
! as land points not surrounded by ocean or ice
!
  do j = 1,nlat
    cs(j) = cos(mlatcnts(j)*pi/180.)
    ss(J) = sin(mlatcnts(j)*pi/180.)
    do i = 1,plon
      cont(i,j) = 0.
      if (nint(oro(i,j)) .eq. 1) then
        cont(i,j) = 1.
      endif
    end do
    temp(1,j) = cont(1,j)
    temp(plon,j) = cont(plon,j)
  end do

  do i = 1,plon
    temp(i,1) = cont(i,1)
    temp(i,nlat) = cont(i,nlat)
  end do
!
! get rid of one and two point islands
!
  do j = 2,nlat-1
    do i = 2,plon-1
      sum =  cont(i  ,j+1) + cont(i  ,j-1) &
           + cont(i+1,j+1) + cont(i+1,j-1) &
           + cont(i-1,j+1) + cont(i-1,j-1) &
           + cont(i+1,j  ) + cont(i-1,j) &
           + cont(i  ,j  )
      if (sum.le.2.) then
        temp(i,j) = 0.
      else
        temp(i,j) = 1.
      endif
    enddo
  end do

  do j = 1,nlat
    do i = 1,plon
      cont(i,j) = temp(i,j)
    end do
  end do
!
! construct a function which is one over land, 
! zero over ocean points beyond dmax from land
!
  iplm1 = 2*plon - 1
  dy = pi*arad/nlat
  jof = dmax/dy + 1
!      write (6,*) ' lat bands to check ', 2*jof+1
  do j = 1,nlat
    c1 = cs(j)
    s1 = ss(j)
    dx = 2*pi*arad*cs(j)/plon
!
! if dx is too small, int(dmax/dx) may exceed the maximum size
! of an integer, especially on Suns, causing a core dump. Test
! to avoid that.
!
    if (dx .lt. 1. .and. dmax .gt. 10000.) then
      iof = plon
    else
      iof = min(int(dmax/dx) + 1, plon)
    end if
    do i = 1,plon
      temp(i,j) = 0.
      landm(i,j) = 0.
      jmin = max(1,j-jof)
      jmax = min(nlat,j+jof)
      do jj = jmin, jmax
        s2 = ss(jj)
        c2 = cs(jj)
        do itmp = -iof,iof
          ii = mod(i+itmp+iplm1,plon)+1
          term = s1*s2 + c1*c2*cos(clon(ii)-clon(i))
          if (term.gt.0.9999999) term = 1.
          dist = arad*acos(term)
          landm(i,j) = max(landm(i,j), (1.-dist/dmax)*cont(ii,jj))
!          if (dist.lt.dmax .and. cont(ii,jj).eq.1) then
!            landm(i,j) = max(landm(i,j), 1.-dist/dmax)
!          endif
        end do
      end do
    end do
  end do
!
! Interpolate to reduced grid.  Redefine clon in terms of degrees for interpolation
!
  do i = 1,plon
    clon(i) = (i-1)*360./plon
  end do
  do j=1,nlat
    do i=1,nlon_reduced(j)
      clon_reduced(i,j) = (i-1)*360./nlon_reduced(j)
    end do
  end do

  do j=1,nlat
    call lininterp (landm(1,j), plon, 1, clon, &
                    landm_reduced(1,j), nlon_reduced(j), 1, clon_reduced(1,j), .true.)
  end do

  return
  end
