module advect_um_lib

! function & subroutine library for the 5th ULTIMATE-MACHO scheme

  use grid
  use params, only: crm_rknd
  implicit none

  logical, parameter :: fct = .true.  ! apply FCT for monotone

  ! Used to judge if courant unuber needs to be updated.
  ! Courant number is same for same icycle.
  integer :: nstep_adv = 0
  logical, dimension(4) :: updated_cn = .false. ! maximum icycle is 4. See kurant.f90

  ! Courant number
  real(crm_rknd), dimension(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm) :: cu
  real(crm_rknd), dimension(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm) :: cv
  real(crm_rknd), dimension(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) :: cw

  ! Inverse of adz, adzw, rho and rhow*adz, common for same icycle
  real(crm_rknd), dimension(nzm) :: iadz, iadzw, irho, irhow

  ! f for advective form update, face values
  real(crm_rknd), dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) :: fadv
  real(crm_rknd), dimension(0:nxp2,          dimy1_s:dimy2_s, nzm) :: fx
  real(crm_rknd), dimension(dimx1_s:dimx2_s, 0:nyp2,          nzm) :: fy
  real(crm_rknd), dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nz ) :: fz

contains

  !--------------------------------------------------------------------------------------------------

  real(crm_rknd) function face_2nd( f_im1, f_i, cn )

    ! Returns face value at left side of i-th control volume, f_i

    implicit none

    ! input
    real(crm_rknd), intent(in) :: f_im1, f_i, cn

    ! Face value
    face_2nd = 0.5 * ( f_i + f_im1 - cn * ( f_i - f_im1 ) )

  end function face_2nd

  !--------------------------------------------------------------------------------------------------

  real(crm_rknd) function face_2nd_z( f_im1, f_i, cn, i )

    ! Returns face value for non-uniform grid. Only used for vertial grid

    implicit none

    ! input
    real(crm_rknd), intent(in) :: f_im1, f_i, cn
    integer, intent(in) :: i

    ! Face value
    face_2nd_z = 0.5 * ( f_i + f_im1 - cn * adz(icrm,i) * iadzw(i) * ( f_i - f_im1 ) )

  end function face_2nd_z

  !--------------------------------------------------------------------------------------------------

  real(crm_rknd) function face_3rd( f_im2, f_im1, f_i, f_ip1, cn )

    implicit none

    ! input
    real(crm_rknd), intent(in) :: f_im2, f_im1, f_i, f_ip1, cn

    ! local
    real(crm_rknd) :: difference2, difference3

    ! 2nd & 3rd difference
    difference2 = f_ip1 - f_i - f_im1 + f_im2
    difference3 = f_ip1 - 3.*f_i + 3.*f_im1 - f_im2

    ! Face value
    face_3rd = 0.5 * ( f_i + f_im1 - cn * ( f_i - f_im1 ) &
      + 1./6. * (cn*cn-1.) * ( difference2 - sign(1._crm_rknd,cn) * difference3 ) )

  end function face_3rd

  !--------------------------------------------------------------------------------------------------

  real(crm_rknd) function face_3rd_z( f_im2, f_im1, f_i, f_ip1, cn, i )

    implicit none

    ! input
    real(crm_rknd), intent(in) :: f_im2, f_im1, f_i, f_ip1, cn
    integer, intent(in) :: i

    ! local
    real(crm_rknd) :: positive_3rd, negative_3rd

    positive_3rd = 1./6. * ( cn * cn * adz(icrm,i) * adz(icrm,i) * iadzw(i) * iadzw(i-1) - 1. ) &
      * ( adzw(icrm,i-1) * iadz(i-1) * ( f_i - f_im1 ) - adzw(icrm,i) * iadz(i-1) * ( f_im1 - f_im2 ) )
    negative_3rd = 1./6. * ( cn * cn * adz(icrm,i) * adz(icrm,i) * iadzw(i+1) * adzw(icrm,i) - 1. ) &
      * ( adzw(icrm,i) * iadz(i) * ( f_ip1 - f_i ) - adzw(icrm,i+1) * iadz(i) * ( f_i - f_im1 ) )

    ! Face value
    face_3rd_z = 0.5 * ( f_i + f_im1 - cn * adz(icrm,i) * iadzw(i) * ( f_i - f_im1 ) &
      + positive_3rd + negative_3rd + sign(1._crm_rknd,cn) * ( positive_3rd - negative_3rd ) )

  end function face_3rd_z

  !--------------------------------------------------------------------------------------------------

  real(crm_rknd) function face_5th( f_im3, f_im2, f_im1, f_i, f_ip1, f_ip2, cn )

    implicit none

    ! input
    real(crm_rknd), intent(in) :: f_im3, f_im2, f_im1, f_i, f_ip1, f_ip2, cn

    ! local
    real(crm_rknd) :: difference2, difference3, difference4, difference5

    ! 2-5th difference
    difference2 = f_ip1 - f_i - f_im1 + f_im2
    difference3 = f_ip1 - 3.*f_i + 3.*f_im1 - f_im2
    difference4 = f_ip2 - 3.*f_ip1 + 2.*f_i + 2.*f_im1 - 3.*f_im2 + f_im3
    difference5 = f_ip2 - 5.*f_ip1 + 10.*f_i - 10.*f_im1 + 5.*f_im2 - f_im3

    ! Face value
    face_5th = 0.5 * ( f_i + f_im1 - cn * ( f_i - f_im1 ) &
      + 1./6. * (cn*cn-1.) * ( difference2 - 0.5 * cn * difference3 ) &
      + 1./120. * (cn*cn-1.) * (cn*cn-4.) * ( difference4 - sign(1._crm_rknd,cn) * difference5 ) )

  end function face_5th

  !--------------------------------------------------------------------------------------------------

  real(crm_rknd) function face_5th_z( f_im3, f_im2, f_im1, f_i, f_ip1, f_ip2, cn, i )

    implicit none

    ! input
    real(crm_rknd), intent(in) :: f_im3, f_im2, f_im1, f_i, f_ip1, f_ip2, cn
    integer, intent(in) :: i

    ! local
    real(crm_rknd) :: positive_5th, negative_5th

    positive_5th = 1./120. * ( cn * cn * adz(icrm,i) * adz(icrm,i) * iadzw(i) * iadzw(i-1) - 1. ) &
      * ( cn * cn * adz(icrm,i) * adz(icrm,i) * iadzw(i+1) * iadzw(i-2) - 4. ) &
      * ( adzw(icrm,i-1) * adzw(icrm,i-2) * iadz(i) * iadz(i-1) * ( f_ip1 - f_i ) &
        - adzw(icrm,i+1) * adzw(icrm,i-2) * (adzw(icrm,i-1) * adz(icrm,i-1) + adzw(icrm,i-1) * adz(icrm,i) + adzw(icrm,i) * adz(icrm,i)) &
          * iadzw(i) * iadz(i) * iadz(i-1) * iadz(i-1) * ( f_i - f_im1 ) &
        + adzw(icrm,i+1) * adzw(icrm,i-2) * (adzw(icrm,i-1) * adz(icrm,i-2) + adzw(icrm,i) * adz(icrm,i-2) + adzw(icrm,i) *adz(icrm,i-1))&
          * iadzw(i-1) * iadz(i-1) * iadz(i-1) * iadz(i-2) * ( f_im1 - f_im2 ) &
        - adzw(icrm,i+1) * adzw(icrm,i) * iadz(i-1) * iadz(i-2) * ( f_im2 - f_im3 ) )

    negative_5th = 1./120. * ( cn * cn * adz(icrm,i) * adz(icrm,i) * iadzw(i+1) * iadzw(i) - 1. ) &
      * ( cn  * cn * adz(icrm,i) * adz(icrm,i) * iadzw(i+2) * iadzw(i-1) - 4. ) &
      * ( adzw(icrm,i) * adzw(icrm,i-1) * iadz(i+1) * iadz(i) * ( f_ip2 - f_ip1 ) &
        - adzw(icrm,i+2) * adzw(icrm,i-1) * (adzw(icrm,i) * adz(icrm,i) + adzw(icrm,i) * adz(icrm,i+1) + adzw(icrm,i+1) * adz(icrm,i+1)) &
          * iadzw(i+1) * iadz(i+1) * iadz(i) * iadz(i) * ( f_ip1 - f_i ) &
        + adzw(icrm,i+2) * adzw(icrm,i-1) * (adzw(icrm,i) * adz(icrm,i-1) + adzw(icrm,i+1) * adz(icrm,i-1) + adzw(icrm,i+1) *adz(icrm,i))&
          * iadzw(i) * iadz(i) * iadz(i) * iadz(i-1) * ( f_i - f_im1 ) &
        - adzw(icrm,i+2) * adzw(icrm,i+1) * iadz(i) * iadz(i-1) * ( f_im1 - f_im2 ) )

    ! Face value
    face_5th_z = 0.5 * ( f_i + f_im1 - cn * adz(icrm,i) * iadzw(i) * ( f_i - f_im1 ) &
      + 1./3. * ( cn * cn * adz(icrm,i) * adz(icrm,i) * iadzw(i+1) * iadzw(i-1) - 1. ) &
        * ( adzw(icrm,i-1) / ( adz(icrm,i) + adz(icrm,i-1) ) * ( f_ip1 - f_i ) &
          - adzw(icrm,i+1) / ( adz(icrm,i) + adz(icrm,i-1) )  * ( f_im1 - f_im2 ) ) &
      - 1./12. * ( cn * cn * adz(icrm,i) * adz(icrm,i) * iadzw(i+1) * adzw(icrm,i-1) - 1. ) * cn*adz(icrm,i)*iadzw(i)&
        * ( adzw(icrm,i-1) * iadz(i) * ( f_ip1 - f_i ) &
          - adzw(icrm,i+1)*adzw(icrm,i-1)*(adz(icrm,i-1)+adz(icrm,i))*iadzw(i)*iadz(i)*iadz(i-1) * ( f_i - f_im1 ) &
          + adzw(icrm,i+1) * iadz(i-1) * ( f_im1 - f_im2 ) ) &
      + positive_5th + negative_5th + sign(1._crm_rknd,cn) * ( positive_5th - negative_5th ) )

  end function face_5th_z

  !--------------------------------------------------------------------------------------------------

  real(crm_rknd) function advective_cn( cn_left, cn_right )

    ! Returns advective courant number
    implicit none
    real(crm_rknd), intent(in) :: cn_left, cn_right

    ! original method to estimate advective velocity for advective update
    if ( (cn_right > 0.).and.(cn_left >= 0.) ) then
      advective_cn = cn_left
    else if ( (cn_right <= 0.).and.(cn_left < 0.) ) then
      advective_cn = cn_right
    else
      advective_cn = 0.
    endif

  end function advective_cn

  !--------------------------------------------------------------------------------------------------

  subroutine face_x_5th( x1, x2, y1, y2 )

    use grid
    implicit none

    ! input
    integer, intent(in) :: x1, x2, y1, y2

    ! local
    integer :: i, j, k

    do k = 1, nzm
      do j = y1, y2
        do i = x1, x2
          fx(i,j,k) = face_5th( fadv(i-3,j,k), fadv(i-2,j,k), fadv(i-1,j,k), fadv(i,j,k), &
                                fadv(i+1,j,k), fadv(i+2,j,k), cu(i,j,k) )
        enddo
      enddo
    enddo

  end subroutine face_x_5th

  !--------------------------------------------------------------------------------------------------

  subroutine face_y_5th( x1, x2, y1, y2 )

    use grid
    implicit none

    ! input
    integer, intent(in) :: x1, x2, y1, y2

    ! local
    integer :: i, j, k

    do k = 1, nzm
      do j = y1, y2
        do i = x1, x2
          fy(i,j,k) = face_5th( fadv(i,j-3,k), fadv(i,j-2,k), fadv(i,j-1,k), fadv(i,j,k), &
                                fadv(i,j+1,k), fadv(i,j+2,k), cv(i,j,k) )
        enddo
      enddo
    enddo

  end subroutine face_y_5th

  !--------------------------------------------------------------------------------------------------

  subroutine face_z_5th( x1, x2, y1, y2 )

    use grid
    implicit none

    ! input
    integer, intent(in) :: x1, x2, y1, y2

    ! local
    integer :: i, j, k

    do j = y1, y2
      do i = x1, x2
        fz(i,j,2) = face_2nd_z( fadv(i,j,1), fadv(i,j,2), cw(i,j,2), 2 )
        fz(i,j,3) = face_3rd_z( fadv(i,j,1), fadv(i,j,2), fadv(i,j,3), fadv(i,j,4), cw(i,j,3),3)
        fz(i,j,nzm-1) = face_3rd_z( fadv(i,j,nzm-3), fadv(i,j,nzm-2), fadv(i,j,nzm-1), &
                                  fadv(i,j,nzm), cw(i,j,nzm-1), nzm-1 )
        fz(i,j,nzm) = face_2nd_z( fadv(i,j,nzm-1), fadv(i,j,nzm), cw(i,j,nzm), nzm )
      enddo
    enddo
    do k = 4, nzm-2
      do j = y1, y2
        do i = x1, x2
          fz(i,j,k) = face_5th_z( fadv(i,j,k-3), fadv(i,j,k-2), fadv(i,j,k-1), fadv(i,j,k), &
                                  fadv(i,j,k+1), fadv(i,j,k+2), cw(i,j,k), k )
        enddo
      enddo
    enddo

  end subroutine face_z_5th

  !--------------------------------------------------------------------------------------------------

  subroutine adv_form_update_x( x1, x2, y1, y2 )

    use grid
    implicit none

    ! input
    integer, intent(in) :: x1, x2, y1, y2

    ! local
    integer :: i, j, k

    do k = 1, nzm
      do j = y1, y2
        do i = x1, x2
          fadv(i,j,k) = fadv(i,j,k) &
            + advective_cn( cu(i,j,k), cu(i+1,j,k) )  * ( fx(i,j,k) - fx(i+1,j,k) )
        enddo
      enddo
    enddo

  end subroutine adv_form_update_x

  !--------------------------------------------------------------------------------------------------

  subroutine adv_form_update_y( x1, x2, y1, y2 )

    use grid
    implicit none

    ! input
    integer, intent(in) :: x1, x2, y1, y2

    ! local
    integer :: i, j, k

    do k = 1, nzm
      do j = y1, y2
        do i = x1, x2
          fadv(i,j,k) = fadv(i,j,k) &
            + advective_cn( cv(i,j,k), cv(i,j+1,k) ) * ( fy(i,j,k) - fy(i,j+1,k) )
        enddo
      enddo
    enddo

  end subroutine adv_form_update_y

  !--------------------------------------------------------------------------------------------------

  subroutine adv_form_update_z( x1, x2, y1, y2 )

    use grid
    implicit none

    ! input
    integer, intent(in) :: x1, x2, y1, y2

    ! local
    integer :: i, j, k

    do k = 2, nzm-1
      do j = y1, y2
        do i = x1, x2
          fadv(i,j,k) = fadv(i,j,k) &
            + advective_cn( cw(i,j,k), cw(i,j,k+1) ) * ( fz(i,j,k) - fz(i,j,k+1) )
        enddo
      enddo
    enddo

  end subroutine adv_form_update_z

  !--------------------------------------------------------------------------------------------------

  subroutine fct3D( f, u, v, w, flux )

    ! Flux corrected transport to enforce monotonicity
    ! input: u, v, w: mass weighted courant number

    use grid
    implicit none

    ! input & output
    real(crm_rknd), dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm), intent(inout) :: f

    ! input
    real(crm_rknd), dimension(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm, ncrms), intent(in) :: u
    real(crm_rknd), dimension(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm, ncrms), intent(in) :: v
    real(crm_rknd), dimension(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz , ncrms), intent(in) :: w

    ! output
    real(crm_rknd), dimension(nz), intent(out) :: flux

    ! local
    real(crm_rknd), dimension(-1:nxp3, -1:nyp2, nzm) :: flx_x
    real(crm_rknd), dimension(-1:nxp2, -1:nyp3, nzm) :: flx_y
    real(crm_rknd), dimension(-1:nxp2, -1:nyp2, nz ) :: flx_z
    real(crm_rknd), dimension(0:nxp1, 0:nyp1, nzm) :: mn, mx
    integer :: i, j, k, km1, kp1
    real(crm_rknd), parameter :: eps = 1.e-10

    ! Set bottom and top vertical flux zero, also horizontal sum of vertical flux for output
    flx_z(:,:,1)  = 0.
    flx_z(:,:,nz) = 0.
    flux(:) = 0.

    ! min and max bound with f
    do k = 1, nzm
      km1 = max(1,k-1)
      kp1 = min(k+1,nzm)
      do j = 0, nyp1
        do i = 0, nxp1
          mn(i,j,k) = min( f(i,j,k), f(i-1,j,k), f(i+1,j,k), f(i,j-1,k), f(i,j+1,k), &
                           f(i,j,km1), f(i,j,kp1) )

          mx(i,j,k) = max( f(i,j,k), f(i-1,j,k), f(i+1,j,k), f(i,j-1,k), f(i,j+1,k), &
                           f(i,j,km1), f(i,j,kp1) )
        enddo
      enddo
    enddo

    ! 1st order upwind flux
    do k = 1, nzm
      do j = -1, nyp2
        do i = -1, nxp3
          flx_x(i,j,k) = f(i-1,j,k) * max( 0., u(icrm,i,j,k) ) + f(i,j,k) * min( 0., u(icrm,i,j,k) )
        enddo
      enddo

      do j = -1, nyp3
        do i = -1, nxp2
          flx_y(i,j,k) = f(i,j-1,k) * max( 0., v(icrm,i,j,k) ) + f(i,j,k) * min( 0., v(icrm,i,j,k) )
        enddo
      enddo
    enddo
    do k = 2, nzm
      do j = -1, nyp2
        do i = -1, nxp2
          flx_z(i,j,k) = f(i,j,k-1) * max( 0., w(icrm,i,j,k) ) + f(i,j,k) * min( 0., w(icrm,i,j,k) )
        enddo
      enddo
    enddo

    ! 1st order upwind update
    do k = 1, nzm
      do j = -1, nyp2
        do i = -1, nxp2
          f(i,j,k) = f(i,j,k) &
            + ( flx_x(i,j,k) - flx_x(i+1,j,k) + flx_y(i,j,k) - flx_y(i,j+1,k) &
            + ( flx_z(i,j,k) - flx_z(i,j,k+1) ) * iadz(k) ) * irho(k)
        enddo
      enddo
      do j = 1, ny
        do i = 1, nx
          flux(k) = flux(k) + flx_z(i,j,k)
        enddo
      enddo
    enddo

    ! Antidiffusive flux
    do k = 1, nzm
      do j = 0, nyp1
        do i = 0, nxp2
          flx_x(i,j,k) = u(icrm,i,j,k) * fx(i,j,k) - flx_x(i,j,k)
        enddo
      enddo

      do j = 0, nyp2
        do i = 0, nxp1
          flx_y(i,j,k) = v(icrm,i,j,k) * fy(i,j,k) - flx_y(i,j,k)
        enddo
      enddo
    enddo
    do k = 2, nzm
      do j = 0, nyp1
        do i = 0, nxp1
          flx_z(i,j,k) = w(icrm,i,j,k) * fz(i,j,k) - flx_z(i,j,k)
        enddo
      enddo
    enddo

    ! min and max bounds with upwind-updated f
    ! convert mn and mx to outflow and inflow fct scale factor
    do k = 1, nzm
      km1 = max(1,k-1)
      kp1 = min(k+1,nzm)
      do j = 0, nyp1
        do i = 0, nxp1

          mn(i,j,k) = min( f(i,j,k), f(i-1,j,k), f(i+1,j,k), f(i,j-1,k), f(i,j+1,k), &
                           f(i,j,km1), f(i,j,kp1), mn(i,j,k) )

          mx(i,j,k) = max( f(i,j,k), f(i-1,j,k), f(i+1,j,k), f(i,j-1,k), f(i,j+1,k), &
                           f(i,j,km1), f(i,j,kp1), mx(i,j,k) )

          ! total higher-order outflow flux
          ! outflow fct factor
          mn(i,j,k) = ( f(i,j,k) - mn(i,j,k) ) &
                / ( ( max(0.,flx_x(i+1,j,k)) - min(0.,flx_x(i,j,k))  &
                    + max(0.,flx_y(i,j+1,k)) - min(0.,flx_y(i,j,k))  &
                  + ( max(0.,flx_z(i,j,k+1)) - min(0.,flx_z(i,j,k)))*iadz(k) )*irho(k) + eps )

          ! total higher-order inflow flux
          ! inflow fct factor
          mx(i,j,k) = ( mx(i,j,k) - f(i,j,k) ) &
                / ( ( max(0.,flx_x(i,j,k)) - min(0.,flx_x(i+1,j,k)) &
                    + max(0.,flx_y(i,j,k)) - min(0.,flx_y(i,j+1,k)) &
                  + ( max(0.,flx_z(i,j,k)) - min(0.,flx_z(i,j,k+1)))*iadz(k) )*irho(k) + eps )
        enddo
      enddo
    enddo

    ! Limit
    do k = 1, nzm
      do j = 1, ny
        do i = 1, nxp1
          flx_x(i,j,k) = max( 0., flx_x(i,j,k) ) * min( 1., mn(i-1,j,k), mx(i,j,k) ) &
                       + min( 0., flx_x(i,j,k) ) * min( 1., mn(i,j,k), mx(i-1,j,k) )
        enddo
      enddo

      do j = 1, nyp1
        do i = 1, nx
          flx_y(i,j,k) = max( 0., flx_y(i,j,k) ) * min( 1., mn(i,j-1,k), mx(i,j,k) ) &
                       + min( 0., flx_y(i,j,k) ) * min( 1., mn(i,j,k), mx(i,j-1,k) )
        enddo
      enddo
    enddo
    do k = 2, nzm
      do j = 1, ny
        do i = 1, nx
          flx_z(i,j,k) = max( 0., flx_z(i,j,k) ) * min( 1., mn(i,j,k-1), mx(i,j,k) ) &
                       + min( 0., flx_z(i,j,k) ) * min( 1., mn(i,j,k), mx(i,j,k-1) )
        enddo
      enddo
    enddo

    ! Final update
    do k = 1, nzm
      do j = 1, ny
        do i = 1, nx
          f(i,j,k) = f(i,j,k) &
            + ( flx_x(i,j,k) - flx_x(i+1,j,k) + flx_y(i,j,k) - flx_y(i,j+1,k) &
            + ( flx_z(i,j,k) - flx_z(i,j,k+1) ) * iadz(k) ) * irho(k)
          flux(k) = flux(k) + flx_z(i,j,k)
        enddo
      enddo
    enddo

  end subroutine fct3D

  !--------------------------------------------------------------------------------------------------

  subroutine fct2D( f, u, w, flux )

    ! Flux corrected transport to enforce monotonicity
    ! input: u, w: mass weighted courant number

    use grid
    implicit none

    ! input & output
    real(crm_rknd), dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm), intent(inout) :: f

    ! input
    real(crm_rknd), dimension(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm, ncrms), intent(in) :: u
    real(crm_rknd), dimension(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz , ncrms), intent(in) :: w

    ! input & output
    real(crm_rknd), dimension(nz), intent(out) :: flux

    ! local
    real(crm_rknd), dimension(-1:nxp3, nzm) :: flx_x
    real(crm_rknd), dimension(-1:nxp2, nz ) :: flx_z
    real(crm_rknd), dimension(0:nxp1, nzm) :: mn, mx
    integer :: i, k, km1, kp1
    real(crm_rknd), parameter :: eps = 1.e-10

    ! Set bottom and top vertical flux zero, also horizontal sum of vertical flux for output
    flx_z(:,1)  = 0.
    flx_z(:,nz) = 0.
    flux(:) = 0.

    ! min and max bounds with f
    do k = 1, nzm
      km1 = max(1,k-1)
      kp1 = min(k+1,nzm)
      do i = 0, nxp1
        mn(i,k) = min( f(i,1,k), f(i-1,1,k), f(i+1,1,k), f(i,1,km1), f(i,1,kp1) )
        mx(i,k) = max( f(i,1,k), f(i-1,1,k), f(i+1,1,k), f(i,1,km1), f(i,1,kp1) )
      enddo
    enddo

    ! 1st order upwind face value and residual higher-order flux
    do k = 1, nzm
      do i = -1, nxp3
        flx_x(i,k) = f(i-1,1,k) * max( 0., u(icrm,i,1,k) ) + f(i,1,k) * min( 0., u(icrm,i,1,k) )
      enddo
    enddo
    do k = 2, nzm
      do i = -1, nxp2
        flx_z(i,k) = f(i,1,k-1) * max( 0., w(icrm,i,1,k) ) + f(i,1,k) * min( 0., w(icrm,i,1,k) )
      enddo
    enddo

    ! 1st order upwind update value
    do k = 1, nzm
      do i = -1, nxp2
        f(i,1,k) = f(i,1,k) + ( flx_x(i,k) - flx_x(i+1,k) &
                            + ( flx_z(i,k) - flx_z(i,k+1) ) * iadz(k) ) * irho(k)
      enddo
      do i = 1, nx
        flux(k) = flux(k) + flx_z(i,k)
      enddo
    enddo

    ! Antidiffusive flux
    do k = 1, nzm
      do i = 0, nxp2
        flx_x(i,k) = u(icrm,i,1,k) * fx(i,1,k) - flx_x(i,k)
      enddo
    enddo
    do k = 2, nzm
      do i = 0, nxp1
        flx_z(i,k) = w(icrm,i,1,k) * fz(i,1,k) - flx_z(i,k)
      enddo
    enddo

    ! min and max bounds
    ! Convert mn and mx to outflow and inflow fct scale factor
    do k = 1, nzm
      km1 = max(1,k-1)
      kp1 = min(k+1,nzm)
      do i = 0, nxp1
        mn(i,k) = min( f(i,1,k), f(i-1,1,k), f(i+1,1,k), f(i,1,km1), f(i,1,kp1), mn(i,k) )
        mx(i,k) = max( f(i,1,k), f(i-1,1,k), f(i+1,1,k), f(i,1,km1), f(i,1,kp1), mx(i,k) )

        ! total higher-order outflow flux
        ! outflow fct factor
        mn(i,k) = ( f(i,1,k) - mn(i,k) ) / ( ( max( 0., flx_x(i+1,k) ) - min( 0., flx_x(i,k) ) &
           + ( max( 0., flx_z(i,k+1) ) - min( 0., flx_z(i,k) ) ) * iadz(k) ) * irho(k) + eps )

        ! total higher-order inflow flux
        ! inflow fct factor
        mx(i,k) = ( mx(i,k) - f(i,1,k) ) / ( ( max( 0., flx_x(i,k) ) - min( 0., flx_x(i+1,k) ) &
           + ( max( 0., flx_z(i,k) ) - min( 0., flx_z(i,k+1) ) ) * iadz(k) ) * irho(k) + eps )
      enddo
    enddo

    ! Limit
    do k = 1, nzm
      do i = 1, nxp1
        flx_x(i,k) = max( 0., flx_x(i,k) ) * min( 1., mn(i-1,k), mx(i,k) ) &
                   + min( 0., flx_x(i,k) ) * min( 1., mn(i,k), mx(i-1,k) )
      enddo
    enddo
    do k = 2, nzm
      do i = 1, nx
        flx_z(i,k) = max( 0., flx_z(i,k) ) * min( 1., mn(i,k-1), mx(i,k) ) &
                   + min( 0., flx_z(i,k) ) * min( 1., mn(i,k), mx(i,k-1) )
      enddo
    enddo

    ! Final updatex
    do k = 1, nzm
      do i = 1, nx
        f(i,1,k) = f(i,1,k) + ( flx_x(i,k) - flx_x(i+1,k) &
                            + ( flx_z(i,k) - flx_z(i,k+1) ) * iadz(k) ) * irho(k)
        flux(k) = flux(k) + flx_z(i,k)
      enddo
    enddo

  end subroutine fct2D

  !--------------------------------------------------------------------------------------------------

end module advect_um_lib
