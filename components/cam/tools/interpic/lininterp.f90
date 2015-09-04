subroutine lininterp (arrin, nxin, incin, xin, &
                      arrout, nxout, incout, xout, periodic)
!-----------------------------------------------------------------------
!
! Do a linear interpolation from input mesh defined by xin to output
! mesh defined by xout.  Where extrapolation is necessary, values will
! be copied from the extreme edge of the input grid.
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Rosinski
!
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
!-----------------------------------------------------------------------
!
! Arguments
!
  integer nxin, incin
  integer nxout, incout

  real(r8) xin(nxin), xout(nxout)
  real(r8) arrin(incin,nxin)
  real(r8) arrout(incout,nxout)

  logical periodic
!
! Local workspace
!
  integer i, ii                ! input grid indices
  integer im, ip, iiprev       ! input grid indices
  integer icount               ! number of values

  real(r8) extrap                ! percent grid non-overlap
  real(r8) dxinwrap              ! delta-x on input grid for 2-pi
  real(r8) avgdxin               ! avg input delta-x
  real(r8) ratio                 ! compare dxinwrap to avgdxin
!
! Dynamic
!
  integer iim(nxout)           ! interp. indices minus
  integer iip(nxout)           ! interp. indices plus

  real(r8) wgtm(nxout)           ! interp. weight minus
  real(r8) wgtp(nxout)           ! interp. weight plus
!
! Just copy the data and return if input dimensions are 1
!
  if (nxin.eq.1 .and. nxout.eq.1) then
    arrout(1,1) = arrin(1,1)
  else if (nxin.eq.1) then
    write(6,*)'LININTERP: Must have at least 2 input points'
    call abort
  end if
  icount = 0
  do i=1,nxin-1
    if (xin(i).gt.xin(i+1)) icount = icount + 1
  end do
  do i=1,nxout-1
    if (xout(i).gt.xout(i+1)) icount = icount + 1
  end do
  if (icount.gt.0) then
    write(6,*)'LININTERP: Non-monotonic coordinate array(s) found'
    call abort
  end if
!
! Initialize index arrays for later checking
!
  do i=1,nxout
    iim(i) = 0
    iip(i) = 0
  end do
  if (periodic) then
!
! Periodic case: for values which extend beyond boundaries, assume 
! periodicity and interpolate between endpoints.  First check for sane 
! periodicity assumption.
!
    if (xin(nxin).gt.360.) then
      write(6,*)'LININTERP: Periodic input x-grid must not be greater than 360'
      call abort
    end if
    if (xout(nxout).gt.360.) then
      write(6,*)'LININTERP: Output x-grid must not be greater than 360'
      call abort
    end if
    dxinwrap = xin(1) + 360. - xin(nxin)
    avgdxin = (xin(nxin)-xin(1))/(nxin-1.)
    ratio = dxinwrap/avgdxin
    if (ratio.lt.0.9 .or. ratio.gt.1.1) then
      write(6,*)'LININTERP: Insane dxinwrap value =',dxinwrap,' avg=', avgdxin
      call abort
    end if
    do im=1,nxout
      if (xout(im).gt.xin(1)) exit
      iim(im) = nxin
      iip(im) = 1
      wgtm(im) = (xin(1)        - xout(im)) /dxinwrap
      wgtp(im) = (xout(im)+360. - xin(nxin))/dxinwrap
    end do
    do ip=nxout,1,-1
      if (xout(ip).le.xin(nxin)) exit
      iim(ip) = nxin
      iip(ip) = 1
      wgtm(ip) = (xin(1)+360. - xout(ip)) /dxinwrap
      wgtp(ip) = (xout(ip)    - xin(nxin))/dxinwrap
    end do
  else
!
! Non-periodic case: for values which extend beyond boundaries, set weights
! such that values will just be copied.
!
    do im=1,nxout
      if (xout(im).gt.xin(1)) exit
      iim(im) = 1
      iip(im) = 1
      wgtm(im) = 1.
      wgtp(im) = 0.
    end do
    do ip=nxout,1,-1
      if (xout(ip).le.xin(nxin)) exit
      iim(ip) = nxin
      iip(ip) = nxin
      wgtm(ip) = 1.
      wgtp(ip) = 0.
    end do
  end if
!
! Loop though output indices finding input indices and weights
!
  iiprev = 1
  do i=im,ip
    do ii=iiprev,nxin-1
      if (xout(i).gt.xin(ii) .and. xout(i).le.xin(ii+1)) then
        iim(i) = ii
        iip(i) = ii + 1
        wgtm(i) = (xin(ii+1)-xout(i))/(xin(ii+1)-xin(ii))
        wgtp(i) = (xout(i)-xin(ii))/(xin(ii+1)-xin(ii))
        goto 30
      end if
    end do
    write(6,*)'LININTERP: Failed to find interp values'
30  iiprev = ii
  end do
!
! Check grid overlap
!
  extrap = 100.*((im - 1.) + (nxout - ip))/nxout
  if (extrap.gt.30.) then
    write(6,*)'********LININTERP WARNING:',extrap,' % of output', &
              ' grid will have to be extrapolated********'
  end if
!
! Check that interp/extrap points have been found for all outputs
!
  icount = 0
  do i=1,nxout
    if (iim(i).eq.0 .or. iip(i).eq.0) icount = icount + 1
  end do
  if (icount.gt.0) then
    write(6,*)'LININTERP: Point found without interp indices'
    call abort
  end if
!
! Do the interpolation
!
  do i=1,nxout
    arrout(1,i) = arrin(1,iim(i))*wgtm(i) + arrin(1,iip(i))*wgtp(i)
  end do
  return
end subroutine lininterp
 
