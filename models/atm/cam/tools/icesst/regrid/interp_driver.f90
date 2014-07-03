subroutine interp_driver (nxi, nyi, nzi, numxi, xposi, yposi, zposi, xiziyi, &
                          nxo, nyo, nzo, numxo, xposo, yposo, zposo, xozoyo)
!-----------------------------------------------------------------------
!
! Driver code for linear interpolation
!
!-----------------------------------------------------------------------
  use precision, only: r8

  implicit none
!
! Arguments
!
  integer, intent(in) :: nxi, nyi, nzi            ! dimensions of input grid
  integer, intent(in) :: numxi(nyi)               ! number of x points per y
  real(r8), intent(in) :: xposi(nxi,nyi)          ! x-positions of input grid
  real(r8), intent(in) :: yposi(nyi)              ! y-positions of input grid
  real(r8), intent(in) :: zposi(nzi)              ! z-positions of input grid
  real(r8), intent(in) :: xiziyi(nxi,nzi,nyi)     ! field to be interpolated

  integer, intent(in) :: nxo, nyo, nzo            ! dimensions of output grid
  integer, intent(in) :: numxo(nyo)               ! number of x points per y
  real(r8), intent(in) :: xposo(nxo,nyo)          ! x-positions of output grid
  real(r8), intent(in) :: yposo(nyo)              ! y-positions of output grid
  real(r8), intent(in) :: zposo(nzo)              ! z-positions of output grid
  real(r8), intent(out) :: xozoyo(nxo,nzo,nyo)    ! output field
!
! Local workspace
!
  integer :: i,j,k             ! spatial indices
  integer :: numxis, numxin    ! number of xpoints input to the south, north
  integer :: numxoj            ! numxo(j)
  integer :: jj, jjs, jjn      ! index in y-direction
  integer :: count             ! number of values found

  real(r8) :: wgts, wgtn       ! interpolation weights
!
! Intermediate interpolation arrays
!
  real(r8) :: xizoyi(nxi,nzo,nyi)
  real(r8) :: xtemp(nxo,2)
!
! Interpolate in z
!
  if (nzi == 1 .and. nzo == 1) then
    xizoyi(:,1,:) = xiziyi(:,1,:)
  else
    do j=1,nyi
      do i=1,nxi
        call lininterp (xiziyi(i,1,j), nzi, nxi, zposi, &
                        xizoyi(i,1,j), nzo, nxi, zposo, .false.)
      end do
    end do
  end if
!
! Check monotonicity of y-coordinate variable before interpolating.  z and
! x monotonicity is checked inside lininterp.
!
  count = 0
  do j=1,nyi-1
    if (yposi(j) > yposi(j+1)) count = count + 1
  end do
  do j=1,nyo-1
    if (yposo(j) > yposo(j+1)) count = count + 1
  end do

  if (count > 0) then
    call err_exit ('interp_driver: non-monotonic coordinate array(s) found')
  end if
!
! Interpolate in x and y
!
  do j=1,nyo
    numxoj = numxo(j)

    jjs = -1
    jjn = -1
    
    if (yposi(1) >= yposo(j)) then          ! extrapolate south
      jjs = 1
      jjn = 1
    else if (yposi(nyi) < yposo(j)) then    ! extrapolate north
      jjs = nyi
      jjn = nyi
    else                                    ! interpolate
      do jj=1,nyi-1
        if (yposi(jj  ) < yposo(j) .and. &
            yposi(jj+1) >= yposo(j)) then
          jjs = jj
          jjn = jj+1
          exit
        end if
      end do
    end if

    if (jjs < 0 .or. jjn < 0) then
      call err_exit ('interp_driver: bad index calculation')
    end if

    numxis = numxi(jjs)
    numxin = numxi(jjn)

    if (jjs /= jjn) then
      wgts = (yposi(jjn) - yposo(j)) / (yposi(jjn) - yposi(jjs))
      wgtn = (yposo(j) - yposi(jjs)) / (yposi(jjn) - yposi(jjs))
      if (abs ((wgts+wgtn)-1.) > 1.e-6) then
        call err_exit ('interp_driver: bad weight calculation')
      end if
    end if

    do k=1,nzo
!
! X interp
!
      call lininterp (xizoyi(1,k,jjs), numxis, 1, xposi(1,jjs), &
                      xtemp(1,1),      numxoj, 1, xposo(1,j), .true.)
      if (jjs == jjn) then
        xozoyo(:numxoj,k,j) = xtemp(:numxoj,1)
      else
        call lininterp (xizoyi(1,k,jjn), numxin, 1, xposi(1,jjn), &
                        xtemp(1,2),      numxoj, 1, xposo(1,j), .true.)
!
! Y interp
!
        xozoyo(:numxoj,k,j) = xtemp(:numxoj,1)*wgts + xtemp(:numxoj,2)*wgtn
      end if
    end do     ! k=1,nzo
  end do       ! j=1,nyo

  return
end subroutine interp_driver
