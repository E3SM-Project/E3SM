subroutine interp_driver (xiziyi, xozoyo, vari, varo, nxi, &
                          nzi, nyi, nxo, nzo, nyo)
  use shr_kind_mod, only: r8 => shr_kind_r8
  use varspecs_mod, only: varspecs

  implicit none
!
! Input arguments
!
  integer nxi, nzi, nyi
  integer nxo, nzo, nyo

  type(varspecs) :: vari, varo

  real(r8) :: xiziyi(nxi,nzi,nyi)
  real(r8) :: xozoyo(nxo,nzo,nyo)
!
! Local workspace
!
  integer :: i,j,k
  integer :: numxis, numxin
  integer :: numxo
  integer :: jj, jjs, jjn
  integer :: count

  real(r8) :: wgts, wgtn
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
        call lininterp (xiziyi(i,1,j), nzi, nxi, vari%zpos, &
                        xizoyi(i,1,j), nzo, nxi, varo%zpos, .false.)
      end do
    end do
  end if
!
! Check monotonicity of y-coordinate variable before interpolating.  z and
! x monotonicity is checked inside lininterp.
!
  count = 0
  do j=1,nyi-1
    if (vari%ypos(j) > vari%ypos(j+1)) count = count + 1
  end do
  do j=1,nyo-1
    if (varo%ypos(j) > varo%ypos(j+1)) count = count + 1
  end do

  if (count > 0) then
    call err_exit ('interp_driver: non-monotonic coordinate array(s) found')
  end if
!
! Interpolate in x and y
!
  do j=1,nyo
    numxo = varo%numx(j)

    jjs = -1
    jjn = -1
    
    if (vari%ypos(1) >= varo%ypos(j)) then
      jjs = 1
      jjn = 1
    else if (vari%ypos(nyi) < varo%ypos(j)) then
      jjs = nyi
      jjn = nyi
    else
      do jj=1,nyi-1
        if (vari%ypos(jj  ) < varo%ypos(j) .and. &
            vari%ypos(jj+1) >= varo%ypos(j)) then
          jjs = jj
          jjn = jj+1
          exit
        end if
      end do
    end if

    if (jjs < 0 .or. jjn < 0) then
      call err_exit ('interp_driver: bad index calculation')
    end if

    numxis = vari%numx(jjs)
    numxin = vari%numx(jjn)

    if (jjs /= jjn) then
      wgts = (vari%ypos(jjn) - varo%ypos(j)) / (vari%ypos(jjn) - vari%ypos(jjs))
      wgtn = (varo%ypos(j) - vari%ypos(jjs)) / (vari%ypos(jjn) - vari%ypos(jjs))
      if (abs ((wgts+wgtn)-1.) > 1.e-6) then
        call err_exit ('interp_driver: bad weight calculation')
      end if
    end if

    do k=1,nzo
!
! X interp
!
      call lininterp (xizoyi(1,k,jjs), numxis, 1, vari%xpos(1,jjs), &
                      xtemp(1,1),      numxo,  1, varo%xpos(1,j), .true.)
      if (jjs == jjn) then
        xozoyo(:numxo,k,j) = xtemp(:numxo,1)
      else
        call lininterp (xizoyi(1,k,jjn), numxin, 1, vari%xpos(1,jjn), &
                        xtemp(1,2),      numxo,  1, varo%xpos(1,j), .true.)
!
! Y interp
!
        xozoyo(:numxo,k,j) = xtemp(:numxo,1)*wgts + xtemp(:numxo,2)*wgtn
      end if
    end do
  end do

  return
end subroutine interp_driver
