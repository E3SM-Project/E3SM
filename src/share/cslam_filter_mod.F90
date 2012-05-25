#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!MODULE CSLAM_FILTER_MOD------------------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
! This module contains everything  to do (ONLY) a CUBIC (3rd order) filtering       ! 
! which expect coefficient from the cubic reconstruciton                            !
!-----------------------------------------------------------------------------------!
module cslam_filter_mod

#ifndef MESH

  use kinds, only                  : int_kind, real_kind
  use dimensions_mod, only         : nc,nhc,nhe
  use coordinate_systems_mod, only : cartesian2D_t,cartesian3D_t
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  use parallel_mod, only : haltmp
  implicit none
  private
  public :: monotonic_gradient_cart
contains
! ----------------------------------------------------------------------------------!
!SUBROUTINE MONOTONIC_GRADIENT_CART-------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
! DESCRIPTION: apply a monotonic filter to the calculated gradient from the cubic   !
!              reconstruction                                                       !
!                                                                                   !
! CALLS: monotonic_interior, monotonic_exterior, monotonic_exteriorswap             !
! INPUT: fcube  ...  tracer values incl. the halo zone                              !
!        cslam  ...  structure incl. tracer values aso                              !
!        desc   ...  edge descriptor, needed for reordering                         !
! INPUT/OUTPUT:                                                                     !
!        recons ...  array of reconstructed coefficients, which will be scaled by   !
!                    the filter                                                     !
! REMARK:                                                                           !
!   This monotonizing scheme is based on the monotone scheme for unstructured       !
!   grids of 'The design and application of upwind schemes on unstructed meshes',   !
!   Barth and Jespersen, AIAA-89-0366, 27th Aerospace Sciences Meeting,             !
!   January 9-12, 1989.                                                             !
!-----------------------------------------------------------------------------------!
subroutine monotonic_gradient_cart(fcube,cslam,recons,desc)
  use cslam_control_volume_mod, only: cslam_struct
  use edge_mod, only : edgedescriptor_t
  
  implicit none
  real (kind=real_kind), dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc),intent(in)    :: fcube
  type (cslam_struct), intent(in)                                            :: cslam   
  real (kind=real_kind),dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe),intent(inout) :: recons  
  type (edgedescriptor_t),intent(in)                                         :: desc
  
  integer (kind=int_kind)          :: i, j, m, n,tmpi,tmpj
  real (kind=real_kind)            :: local_min, local_max, value, phi, min_phi
  real (kind=real_kind)            :: disc, mx, my
  
  call monotonic_int(fcube,recons,cslam%acartx,cslam%acarty,cslam%spherecentroid, &
                     cslam%jx_min,cslam%jx_max,cslam%jy_min,cslam%jy_max,cslam%cubeboundary)

  ! all elements, which have a cube edge
  if(cslam%cubeboundary>0)then
    if(cslam%swap1) then
      call monotonic_haloswap(fcube,recons,cslam%acartx1,cslam%acarty1,&
                                  cslam%spherecentroid,cslam%jx_min1,cslam%jx_max1,&
                                  cslam%jy_min1,cslam%jy_max1,cslam%cubeboundary,desc)
    else
      call monotonic_halo(fcube,recons,cslam%acartx1,cslam%acarty1,&
                              cslam%spherecentroid,cslam%jx_min1,cslam%jx_max1,&
                              cslam%jy_min1,cslam%jy_max1,cslam%cubeboundary,desc)
    endif
    ! only for corner elements 
    if(cslam%cubeboundary>4)then
      if(cslam%swap2) then
        call monotonic_haloswap(fcube,recons,cslam%acartx2,cslam%acarty2,&
                                    cslam%spherecentroid,cslam%jx_min2,cslam%jx_max2,&
                                    cslam%jy_min2,cslam%jy_max2,cslam%cubeboundary,desc)
      else
        call monotonic_halo(fcube,recons,cslam%acartx2,cslam%acarty2,&
                                cslam%spherecentroid,cslam%jx_min2,cslam%jx_max2,&
                                cslam%jy_min2,cslam%jy_max2,cslam%cubeboundary,desc)
      endif
    endif
  endif
end subroutine monotonic_gradient_cart
!END SUBROUTINE MONOTONIC_GRADIENT_CART---------------------------------CE-for CSLAM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE MONOTONIC_INT-----------------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
! DESCRIPTION: filtering for interior elements                                      !
!              this subroutine calls an optimized minmax search algorithm, and does !
!              not have to reorder the coefficients and centroids                   !
!                                                                                   !
! CALLS: minmax_patch, recons_val_cart, slopelimiter_val                            !
! INPUT: fcube  ...  tracer values incl. the halo zone                              !
!        cslam  ...  structure incl. tracer values aso                              !
!        acartx ...  x cartesian coordinats of the arrival grid on the cube         !
!        acarty ...  y cartesian coordinats of the arrival grid on the cube         !
!        centroid..  x,y,x^2,y^2,xy                                                 !
!        jx_min, jx_max, jy_min, jy_max                                             !
!               ... cell boundaries in x and y directions                           !
!        cubeboundary. 0-8, depending of an element is interior, west, south, aso   !
! INPUT/OUTPUT:                                                                     !
!        recons ...  array of reconstructed coefficients, which will be scaled by   !
!                    the filter                                                     !
!-----------------------------------------------------------------------------------!
subroutine monotonic_int(fcube,recons,acartx,acarty,centroid,jx_min,jx_max,jy_min,jy_max,&
                         cubeboundary)
  
  implicit none
  real (kind=real_kind), dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc),intent(in)    :: fcube  
  real (kind=real_kind),dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe),intent(inout) :: recons  
  real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in)   :: acartx, acarty
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(in)  :: centroid
  integer, intent(in)                                  :: cubeboundary
  
  integer, intent(in)                                  :: jx_min,jx_max,jy_min,jy_max
  
  integer (kind=int_kind)          :: i, j, m, n
  real (kind=real_kind)            :: local_min, local_max, value, phi, min_phi
  real (kind=real_kind)            :: disc, mx, my

  ! Apply monotone limiting, go through all coefficients
  do j = jy_min, jy_max-1
    do i = jx_min, jx_max-1
      call minmax_patch(fcube, i, j, local_min, local_max, cubeboundary)
      ! Initialize the limiter
      min_phi = 1.0D0
      ! Check the minima and maxima at the corner points of the element
      do m = i, i+1
        do n = j, j+1
          ! Evaluate the function at each corner point
          call recons_val_cart(fcube, acartx(m),acarty(n),centroid,recons,i,j,value)
          call slopelimiter_val(value, fcube(i,j), local_min, local_max, min_phi)
        enddo
      enddo
      ! For the third order method, the minima and maxima may occur along
      ! the line segments given by du/dx = 0 and du/dy = 0.  Also check
      ! for the presence of a maxima / minima of the quadratic within
      ! the domain.
      disc =  4.0D0 * recons(4,i,j) * recons(3,i,j) - recons(5,i,j)**2
      ! Check if the quadratic is minimized within the element
      ! Extrema in the interior of the element (there might be just one candidate)
      ! DO NOT NEED ABS here, if disc<0 we have a saddle point (no maximum or minimum)
      if (abs(disc) > 1.0D-12) then
        mx = recons(5,i,j) * recons(2,i,j)        &
             - 2.0D0 * recons(4,i,j) * recons(1,i,j)
        my = recons(5,i,j) * recons(1,i,j)        &
             - 2.0D0 * recons(3,i,j) * recons(2,i,j)

        mx = mx / disc + centroid(1,i,j)
        my = my / disc + centroid(2,i,j)

        if ((mx - acartx(i)   > -1.0D-12) .and. &
            (mx - acartx(i+1) <  1.0D-12) .and. &
            (my - acarty(j)   > -1.0D-12) .and. &
            (my - acarty(j+1) <  1.0D-12)) then       
          call recons_val_cart(fcube, mx, my, centroid, recons, i, j, value)
          call slopelimiter_val(value, fcube(i,j),local_min, local_max, min_phi)
        endif
      endif
      ! Check all potential minimizer points along element boundaries
      if (abs(recons(5,i,j)) > 1.0D-12) then
        ! Left/right edge, intercept with du/dx = 0
        do m = i, i+1
          my = - recons(1,i,j) - 2.0D0 * recons(3,i,j) * &
                 (acartx(m) - centroid(1,i,j))
          my = my / recons(5,i,j) + centroid(2,i,j)
          if ((my > acarty(j)) .and. (my < acarty(j+1))) then
            call recons_val_cart(fcube, acartx(m), my, centroid, recons, i, j, value)
            call slopelimiter_val(value, fcube(i,j),local_min, local_max, min_phi)
          endif
        enddo
        ! Top/bottom edge, intercept with du/dy = 0
        do n = j, j+1
          mx = - recons(2,i,j) - 2.0D0 * recons(4,i,j) * &
               (acarty(n) - centroid(2,i,j))
          mx = mx / recons(5,i,j) + centroid(1,i,j)
          if ((mx > acartx(i)) .and. (mx < acartx(i+1))) then
            call recons_val_cart(fcube, mx, acarty(n),centroid,recons, i, j, value)
            call slopelimiter_val(value, fcube(i,j),local_min, local_max, min_phi)
          endif
        enddo
      endif
      ! Top/bottom edge, y=const., du/dx=0
      if (abs(recons(3,i,j)) > 1.0D-12) then
        do n = j, j+1
          mx = - recons(1,i,j) - recons(5,i,j) * &
               (acarty(n) - centroid(2,i,j))

          mx = mx / (2.0D0 * recons(3,i,j)) + centroid(1,i,j)
          if ((mx > acartx(i)) .and. (mx < acartx(i+1))) then
            call recons_val_cart(fcube, mx, acarty(n), centroid, recons, i, j, value)
            call slopelimiter_val(value, fcube(i,j),local_min, local_max, min_phi)
          endif
        enddo
      endif
        ! Left/right edge, x=const., du/dy=0
      if (abs(recons(4,i,j)) > 1.0D-12) then
        do m = i, i+1
          my = - recons(2,i,j) - recons(5,i,j) * &
                 (acartx(m) - centroid(1,i,j))
          my = my / (2.0D0 * recons(4,i,j)) + centroid(2,i,j)
          if ((my > acarty(j)) .and. (my < acarty(j+1))) then
            call recons_val_cart(fcube, acartx(m), my, centroid, recons, i, j, value)
            call slopelimiter_val(value, fcube(i,j),local_min, local_max, min_phi)
          endif
        enddo
      endif
      
      if ((min_phi < -1.0D-12) .or. (min_phi > 1.0D0 + 1.0D-12)) then
        write(*,*) value, fcube(i,j),local_min, local_max, min_phi
        call haltmp("Fatal Error in monotonic_int: slope limiter out of range (not in [0,1])!")        
        stop
      endif
      ! Apply monotone limiter to all reconstruction coefficients
      recons(1,i,j) = min_phi * recons(1,i,j)
      recons(2,i,j) = min_phi * recons(2,i,j)
      recons(3,i,j) = min_phi * recons(3,i,j)
      recons(4,i,j) = min_phi * recons(4,i,j)
      recons(5,i,j) = min_phi * recons(5,i,j)
    enddo
  enddo
end subroutine monotonic_int
!END SUBROUTINE MONOTONIC_INT-------------------------------------------CE-for CSLAM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE MONOTONIC_HALO----------------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
! DESCRIPTION: filtering for elements with an edge on the cube boundary             !
!             for certain elements we have to reorder the coefficients and centroids!
!                                                                                   !
! CALLS: minmax_patch_halo, recons_val_cart, slopelimiter_val                       !
! INPUT: fcube  ...  tracer values incl. the halo zone                              !
!        cslam  ...  structure incl. tracer values aso                              !
!        acartx ...  x cartesian coordinats of the arrival grid on the cube         !
!        acarty ...  y cartesian coordinats of the arrival grid on the cube         !
!        centroid..  x,y,x^2,y^2,xy                                                 !
!        jx_min, jx_max, jy_min, jy_max                                             !
!               ... cell boundaries in x and y directions                           !
!        cubeboundary. 0-8, depending of an element is interior, west, south, aso   !
!        desc   ... edge descriptor, needed for reordering                          !
! INPUT/OUTPUT:                                                                     !
!        recons ...  array of reconstructed coefficients, which will be scaled by   !
!                    the filter                                                     !
!-----------------------------------------------------------------------------------!
subroutine monotonic_halo(fcube,recons,acartx,acarty,centroid,&
                              jx_min,jx_max,jy_min,jy_max,cubeboundary,desc)
  use edge_mod, only : edgedescriptor_t
  
  implicit none
  real (kind=real_kind), dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc),intent(in)    :: fcube  
  real (kind=real_kind),dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe),intent(inout) :: recons  
  real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in)   :: acartx, acarty
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(in)  :: centroid
  integer, intent(in)                                  :: jx_min,jx_max,jy_min,jy_max
  integer, intent(in)                                  :: cubeboundary
  type (edgedescriptor_t),intent(in)                   :: desc
  
  integer (kind=int_kind)          :: i, j, m, n, tmpi, tmpj
  real (kind=real_kind)            :: local_min, local_max, value, phi, min_phi
  real (kind=real_kind)            :: disc, mx, my

  do j = jy_min, jy_max-1
    do i = jx_min, jx_max-1
      if(desc%reverse(west) .or. desc%reverse(east)) then
        tmpj=jy_max-1+jy_min-j
      else
        tmpj=j
      endif
      if(desc%reverse(north) .or. desc%reverse(south)) then
        tmpi=jx_max-1+jx_min-i
      else
        tmpi=i
      endif
      call minmax_patch_halo(fcube, tmpi, tmpj, cubeboundary,local_min, local_max)
      ! Initialize the limiter
      min_phi = 1.0D0
      ! For the second-order calculation, the minima and maxima will occur
      ! at the corner points of the element
      do m = i, i+1
        do n = j, j+1
        ! Evaluate the function at each corner point
          call recons_val_cart(fcube, acartx(m),acarty(n),centroid,recons,tmpi,tmpj,value)
          call slopelimiter_val(value, fcube(tmpi,tmpj), local_min, local_max, min_phi) 
        enddo
      enddo 

      disc = 4.0D0 * recons(4,tmpi,tmpj) * recons(3,tmpi,tmpj) - recons(5,tmpi,tmpj)**2
      ! Check if the quadratic is minimized within the element
      if (abs(disc) > 1.0D-12) then
        mx =  recons(5,tmpi,tmpj) * recons(2,tmpi,tmpj)        &
             - 2.0D0 * recons(4,tmpi,tmpj) * recons(1,tmpi,tmpj)
        my =   recons(5,tmpi,tmpj) * recons(1,tmpi,tmpj)        &
             - 2.0D0 * recons(3,tmpi,tmpj) * recons(2,tmpi,tmpj)

        mx = mx / disc + centroid(1,tmpi,tmpj)
        my = my / disc + centroid(2,tmpi,tmpj)
        if ((mx - acartx(i)   > -1.0D-12) .and. &
            (mx - acartx(i+1) <  1.0D-12) .and. &
            (my - acarty(j)   > -1.0D-12) .and. &
            (my - acarty(j+1) <  1.0D-12)) then       
          call recons_val_cart(fcube, mx, my, centroid, recons, tmpi, tmpj, value)
          call slopelimiter_val(value, fcube(tmpi,tmpj),local_min, local_max, min_phi)
        endif
      endif
      ! Check all potential minimizer points along element boundaries
      if (abs(recons(5,tmpi,tmpj)) > 1.0D-12) then
        ! Left/right edge, intercept with du/dx = 0
        do m = i, i+1
          my = - recons(1,tmpi,tmpj) - 2.0D0 * recons(3,tmpi,tmpj) * &
                (acartx(m) - centroid(1,tmpi,tmpj))
          my = my / recons(5,tmpi,tmpj) + centroid(2,tmpi,tmpj)
          if ((my > acarty(j)) .and. (my < acarty(j+1))) then
            call recons_val_cart(fcube, acartx(m), my, centroid, recons, tmpi, tmpj, value)
            call slopelimiter_val(value, fcube(tmpi,tmpj),local_min, local_max, min_phi)
          endif
        enddo
        ! Top/bottom edge, intercept with du/dy = 0
        do n = j, j+1
          mx = - recons(2,tmpi,tmpj) - 2.0D0 * recons(4,tmpi,tmpj) * &
              (acarty(n) - centroid(2,tmpi,tmpj))
          mx = mx / recons(5,tmpi,tmpj) + centroid(1,tmpi,tmpj)
          if ((mx > acartx(i)) .and. (mx < acartx(i+1))) then
            call recons_val_cart(fcube, mx, acarty(n),centroid,recons, tmpi, tmpj, value)
            call slopelimiter_val(value, fcube(tmpi,tmpj),local_min, local_max, min_phi)
          endif
        enddo
      endif
      ! Top/bottom edge, intercept with du/dx = 0
      if (abs(recons(3,tmpi,tmpj)) > 1.0D-12) then
        do n = j, j+1
          mx = - recons(1,tmpi,tmpj) - recons(5,tmpi,tmpj) * &
              (acarty(n) - centroid(2,tmpi,tmpj))

          mx = mx / (2.0D0 * recons(3,tmpi,tmpj)) + centroid(1,tmpi,tmpj)
          if ((mx > acartx(i)) .and. (mx < acartx(i+1))) then
            call recons_val_cart(fcube, mx, acarty(n), centroid, recons, tmpi, tmpj, value)
            call slopelimiter_val(value, fcube(tmpi,tmpj),local_min, local_max, min_phi)
          endif
        enddo
      endif
       ! Left/right edge, intercept with du/dy = 0
      if (abs(recons(4,tmpi,tmpj)) > 1.0D-12) then
        do m = i, i+1
          my = - recons(2,tmpi,tmpj) - recons(5,tmpi,tmpj) * &
                (acartx(m) - centroid(1,tmpi,tmpj))
          my = my / (2.0D0 * recons(4,tmpi,tmpj)) + centroid(2,tmpi,tmpj)
          if ((my > acarty(j)) .and. (my < acarty(j+1))) then
            call recons_val_cart(fcube, acartx(m), my, centroid, recons, tmpi, tmpj, value)
            call slopelimiter_val(value, fcube(tmpi,tmpj),local_min, local_max, min_phi)
          endif
        enddo
      endif
      if ((min_phi < -1.0D-12) .OR. (min_phi > 1.0D0 + 1.0D-12)) then
        call haltmp("Fatal Error in monotonic_halo: slope limiter out of range (no in [0,1])!")        
        stop
      endif
      recons(1,tmpi,tmpj) = min_phi * recons(1,tmpi,tmpj)
      recons(2,tmpi,tmpj) = min_phi * recons(2,tmpi,tmpj)
      recons(3,tmpi,tmpj) = min_phi * recons(3,tmpi,tmpj)
      recons(4,tmpi,tmpj) = min_phi * recons(4,tmpi,tmpj)
      recons(5,tmpi,tmpj) = min_phi * recons(5,tmpi,tmpj)
    enddo
  enddo
end subroutine monotonic_halo
!END SUBROUTINE MONOTONIC_HALO------------------------------------------CE-for CSLAM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE MONOTONIC_HALOSWAP------------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
! DESCRIPTION: filtering for elements with an edge on the cube boundary             !
!              use this if we have to swap x and y coordinates                      !
!             for certain elements we have to reorder the coefficients and centroids!
!                                                                                   !
! CALLS: minmax_patch_halo, recons_val_cart, slopelimiter_val                       !
! INPUT: fcube  ...  tracer values incl. the halo zone                              !
!        cslam  ...  structure incl. tracer values aso                              !
!        acartx ...  x cartesian coordinats of the arrival grid on the cube         !
!        acarty ...  y cartesian coordinats of the arrival grid on the cube         !
!        centroid..  x,y,x^2,y^2,xy                                                 !
!        jx_min, jx_max, jy_min, jy_max                                             !
!               ... cell boundaries in x and y directions                           !
!        cubeboundary. 0-8, depending of an element is interior, west, south, aso   !
!        desc   ... edge descriptor, needed for reordering                          !
! INPUT/OUTPUT:                                                                     !
!        recons ...  array of reconstructed coefficients, which will be scaled by   !
!                    the filter                                                     !
!-----------------------------------------------------------------------------------!
subroutine monotonic_haloswap(fcube,recons,acartx,acarty,centroid,&
                              jx_min,jx_max,jy_min,jy_max,cubeboundary,desc)
  use edge_mod, only : edgedescriptor_t
  
  implicit none
  real (kind=real_kind), dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc),intent(in)    :: fcube  
  real (kind=real_kind),dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe),intent(inout) :: recons  
  real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in)   :: acartx, acarty
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(in)  :: centroid
  integer, intent(in)                                  :: jx_min,jx_max,jy_min,jy_max
  integer, intent(in)                                  :: cubeboundary
  type (edgedescriptor_t),intent(in)                   :: desc
  
  integer (kind=int_kind)          :: i, j, m, n, tmpi, tmpj
  real (kind=real_kind)            :: local_min, local_max, value, phi, min_phi
  real (kind=real_kind)            :: disc, mx, my

  do j = jy_min, jy_max-1
    do i = jx_min, jx_max-1
      if(desc%reverse(west) .or. desc%reverse(east)) then
        tmpj=jy_max-1+jy_min-j
      else
        tmpj=j
      endif
      if(desc%reverse(north) .or. desc%reverse(south)) then
        tmpi=jx_max-1+jx_min-i
      else
        tmpi=i
      endif
     
      call minmax_patch_halo(fcube, tmpi, tmpj, cubeboundary,local_min, local_max)
     ! Initialize the limiter
      min_phi = 1.0D0
     ! For the second-order calculation, the minima and maxima will occur
     ! at the corner points of the element
      do m = i, i+1
        do n = j, j+1
        ! Evaluate the function at each corner point
          call recons_val_cart(fcube, acarty(n),acartx(m),centroid,recons,tmpi,tmpj,value)
          call slopelimiter_val(value, fcube(tmpi,tmpj), local_min, local_max, min_phi) 
        enddo
      enddo 
      disc = 4.0D0 * recons(4,tmpi,tmpj) * recons(3,tmpi,tmpj) - recons(5,tmpi,tmpj)**2 
      ! Check if the quadratic is minimized within the element
      if (abs(disc) > 1.0D-12) then
        mx =   recons(5,tmpi,tmpj) * recons(2,tmpi,tmpj)        &
             - 2.0D0 * recons(4,tmpi,tmpj) * recons(1,tmpi,tmpj)
        my =   recons(5,tmpi,tmpj) * recons(1,tmpi,tmpj)        &
             - 2.0D0 * recons(3,tmpi,tmpj) * recons(2,tmpi,tmpj)
        mx = mx / disc + centroid(1,tmpi,tmpj)
        my = my / disc + centroid(2,tmpi,tmpj)

        if ((mx - acarty(j)   > -1.0D-12) .and. &
            (mx - acarty(j+1) <  1.0D-12) .and. &
            (my - acartx(i)   > -1.0D-12) .and. &
            (my - acartx(i+1) <  1.0D-12)) then       
          call recons_val_cart(fcube, mx, my, centroid, recons, tmpi, tmpj, value)
          call slopelimiter_val(value, fcube(tmpi,tmpj),local_min, local_max, min_phi)
        endif
      endif
      ! Check all potential minimizer points along element boundaries
      if (abs(recons(5,tmpi,tmpj)) > 1.0D-12) then
        ! Left/right edge, intercept with du/dx = 0
        do m = j, j+1
          my = - recons(1,tmpi,tmpj) - 2.0D0 * recons(3,tmpi,tmpj) * &
                (acarty(m) - centroid(1,tmpi,tmpj))
          my = my / recons(5,tmpi,tmpj) + centroid(2,tmpi,tmpj)
          if ((my > acartx(i)) .and. (my < acartx(i+1))) then
            call recons_val_cart(fcube, acarty(m), my, centroid, recons, tmpi, tmpj, value)
            call slopelimiter_val(value, fcube(tmpi,tmpj),local_min, local_max, min_phi)
          endif
        enddo
        ! Top/bottom edge, intercept with du/dy = 0
        do n = i, i+1
          mx = - recons(2,tmpi,tmpj) - 2.0D0 * recons(4,tmpi,tmpj) * &
                (acartx(n) - centroid(2,tmpi,tmpj))
          mx = mx / recons(5,tmpi,tmpj) + centroid(1,tmpi,tmpj)
          if ((mx > acarty(j)) .and. (mx < acarty(j+1))) then
            call recons_val_cart(fcube, mx, acartx(n),centroid,recons, tmpi, tmpj, value)
            call slopelimiter_val(value, fcube(tmpi,tmpj),local_min, local_max, min_phi)
          endif
        enddo
      endif
      ! Top/bottom edge, intercept with du/dx = 0
      if (abs(recons(3,tmpi,tmpj)) > 1.0D-12) then
        do n = i, i+1
          mx = - recons(1,tmpi,tmpj) - recons(5,tmpi,tmpj) * &
               (acartx(n) - centroid(2,tmpi,tmpj))

          mx = mx / (2.0D0 * recons(3,tmpi,tmpj)) + centroid(1,tmpi,tmpj)
          if ((mx > acarty(j)) .and. (mx < acarty(j+1))) then
            call recons_val_cart(fcube, mx, acartx(n), centroid, recons, tmpi, tmpj, value)
            call slopelimiter_val(value, fcube(tmpi,tmpj),local_min, local_max, min_phi)
          endif
        enddo
      endif
      ! Left/right edge, intercept with du/dy = 0
      if (abs(recons(4,tmpi,tmpj)) > 1.0D-12) then
        do m = j, j+1
          my = - recons(2,tmpi,tmpj) - recons(5,tmpi,tmpj) * &
                 (acarty(m) - centroid(1,tmpi,tmpj))
          my = my / (2.0D0 * recons(4,tmpi,tmpj)) + centroid(2,tmpi,tmpj)
          if ((my > acartx(i)) .and. (my < acartx(i+1))) then
            call recons_val_cart(fcube, acarty(m), my, centroid, recons, tmpi, tmpj, value)
            call slopelimiter_val(value, fcube(tmpi,tmpj),local_min, local_max, min_phi)
          endif
        enddo
      endif
      if ((min_phi < -1.0D-12) .or. (min_phi > 1.0D0 + 1.0D-12)) then
        call haltmp("Fatal Error in monotonic_haloswap: slope limiter out of range (not in [0,1])!")                
        stop
      endif
      recons(1,tmpi,tmpj) = min_phi * recons(1,tmpi,tmpj)
      recons(2,tmpi,tmpj) = min_phi * recons(2,tmpi,tmpj)
      recons(3,tmpi,tmpj) = min_phi * recons(3,tmpi,tmpj)
      recons(4,tmpi,tmpj) = min_phi * recons(4,tmpi,tmpj)
      recons(5,tmpi,tmpj) = min_phi * recons(5,tmpi,tmpj)
    enddo
  enddo
end subroutine monotonic_haloswap
!END SUBROUTINE MONOTONIC_HALOSWAP--------------------------------------CE-for CSLAM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE MINMAX_PATCH------------------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
! DESCRIPTION: search for the min/max in a patch (around a the element (a,b))       !
!              this is optimized and just valid for interior elements               !
!                                                                                   !
! INPUT: fcube  ...  tracer values incl. the halo zone                              !
!        a,b    ...  index of the considered cell                                   !
!        cubeboundary. 0-8, depending of an element is interior, west, south, aso   !
! OUTPUT:                                                                           !
!        minval ...  minmal value in the patch                                      !
!        maxval ...  maximal value in the patch                                     !
!-----------------------------------------------------------------------------------!
subroutine minmax_patch(fcube, a, b, min_val, max_val,cubeboundary)
  implicit none
  real (kind=real_kind), dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc), intent(in) :: fcube
  integer, intent(in)                                                      :: a, b
  real    (kind=real_kind), intent(out)                          :: min_val, max_val
  integer, intent(in)                                  :: cubeboundary


  integer               :: i, j

  min_val = fcube(a,b)
  max_val = fcube(a,b)
  do j = b-1,b+1
    do i = a-1,a+1
      !BLUEFIRE does not allow if statements with NaN
      if (cubeboundary == swest .AND. i<1 .AND. j<1) CYCLE
      if (cubeboundary == seast .AND. i>nc .AND. j<1) CYCLE
      if (cubeboundary == nwest .AND. i<1 .AND. j>nc) CYCLE
      if (cubeboundary == neast .AND. i>nc .AND. j>nc) CYCLE

!FROST does not like the min/max function 
      if(min_val > fcube(i,j)) then
        min_val=fcube(i,j)
      endif
      if(max_val < fcube(i,j)) then
        max_val=fcube(i,j)
      endif
    enddo
  enddo
end subroutine minmax_patch
!END SUBROUTINE MINMAX_PATCH--------------------------------------------CE-for CSLAM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE MINMAX_PATCH_HALO-------------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
! DESCRIPTION: search for the min/max in a patch (around a the element (a,b))       !
!              this is a geneneral subroutine, here only used for the halo zone     !
!                                                                                   !
! INPUT: fcube  ...  tracer values incl. the halo zone                              !
!        a,b    ...  index of the considered cell                                   !
!        cubeboundary. 0-8, depending of an element is interior, west, south, aso   !
! OUTPUT:                                                                           !
!        minval ...  minmal value in the patch                                      !
!        maxval ...  maximal value in the patch                                     !
!-----------------------------------------------------------------------------------!
subroutine minmax_patch_halo(fcube, a, b, cubeboundary, min_val, max_val)
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  implicit none

  real (kind=real_kind), dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc), intent(in) :: fcube

  integer, intent(in)                                            :: a, b
  integer, intent(in)                                            :: cubeboundary  
  
  real    (kind=real_kind), intent(out)                          :: min_val, max_val

  ! Local variables
  integer               :: i, j
  real (kind=real_kind) :: value

  min_val = fcube(a,b)
  max_val = fcube(a,b)
  value   = fcube(a,b)

  do j = b-1,b+1
    do i = a-1,a+1
      !Interior and elements without a corner
      if (cubeboundary <= 4) then   !interior, west, east, south or north element
         value = fcube(i,j)
      else !for corner elements, value might not exist or is on another halo
        !CASE WEST
        !Interior and elements without a corner
        select case (cubeboundary)
          case(swest)
            if ((a < 1) .and. (b==1) .and. (j==0))then
              value = fcube(1,i)
            elseif((a == 1) .and. (b<1) .and. (i==0))then
              value = fcube(j,1)
            elseif(i<1 .AND. j<1) then
              CYCLE 
            else
              value = fcube(i,j)    
            endif
          case(seast)
            if ((a > nc) .and. (b==1) .and. (j==0))then
              value = fcube(nc,nc+1-i)
            elseif((a == nc) .and. (b<1) .and. (i==nc+1))then
              value = fcube(nc+1-j,1)
            elseif(i>nc .AND. j<1) then
               CYCLE
            else
              value = fcube(i,j)    
            endif  
          case(neast)
            if ((a > nc) .and. (b==nc) .and. (j==nc+1))then
              value = fcube(nc,i)
            elseif((a == nc) .and. (b>nc) .and. (i==nc+1))then
              value = fcube(j,nc)
            elseif(i>nc .AND. j>nc) then
              CYCLE
            else
              value = fcube(i,j)    
            endif
          case(nwest)
            if ((a < 1) .and. (b==nc) .and. (j==nc+1))then
              value = fcube(1,nc+1-i)
            elseif((a == 1) .and. (b>nc) .and. (i==0))then
              value = fcube(nc+1-j,nc)
            elseif(i<1 .AND. j>nc) then
               CYCLE
            else
              value = fcube(i,j)    
            endif  
        end select
      end if
!       min_val = min(min_val, value)
!       max_val = max(max_val, value)
      !FROST does not like the min/max function 
      if(min_val>value) then
        min_val=value
      endif
      if(max_val<value) then
        max_val=value
      endif
    enddo
  enddo
end subroutine minmax_patch_halo
!END SUBROUTINE MINMAX_PATCH_HALO---------------------------------------CE-for CSLAM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE RECONS_VAL_CART---------------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
! DESCRIPTION: returns the value from the reconstruction (3rd order Taylor polynom) !
!              at the point (cartx,carty) -> in cube CARTESIAN coordinates          !
!                                                                                   !
! INPUT: fcube  ...  tracer values incl. the halo zone                              !
!        cartx ...  x cartesian coordinate of the evaluation point                  !
!        carty ...  y cartesian coordinate of the evaluation point                  !
!        centroid..  x,y,x^2,y^2,xy                                                 !
!        recons ...  array of reconstructed coefficients                            !
!        a,b   ...   index of the cell, we do the evaluation                        !
! OUTPUT: value ... evaluation at a given point                                     !
!-----------------------------------------------------------------------------------!
subroutine recons_val_cart(fcube, cartx, carty, centroid, recons, a, b, value)

  implicit none
  real (kind=real_kind), dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc), intent(in) :: fcube
  real (kind=real_kind), intent(in)                                   :: cartx, carty 
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe),intent(in) :: centroid
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(in):: recons
  integer (kind=int_kind), intent(in)                                      :: a, b
  real    (kind=real_kind), intent(out)                                    :: value

  ! Evaluate constant order terms
  value = fcube(a,b) + &
  ! Evaluate linear order terms
          recons(1,a,b) * (cartx - centroid(1,a,b)) + &
          recons(2,a,b) * (carty - centroid(2,a,b)) + &
  ! Evaluate second order terms
          recons(3,a,b) * (centroid(1,a,b)**2 - centroid(3,a,b)) + &
          recons(4,a,b) * (centroid(2,a,b)**2 - centroid(4,a,b)) + &
          recons(5,a,b) * (centroid(1,a,b) * centroid(2,a,b) - centroid(5,a,b)) + &

          recons(3,a,b) * (cartx - centroid(1,a,b))**2 + &
          recons(4,a,b) * (carty - centroid(2,a,b))**2 + &
          recons(5,a,b) * (cartx - centroid(1,a,b)) * (carty - centroid(2,a,b))
end subroutine recons_val_cart
!END SUBROUTINE RECONS_VAL_CART-----------------------------------------CE-for CSLAM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE SLOPELIMITER_VAL--------------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 30.November 2011                                         !
! DESCRIPTION: returns the value from the reconstruction (3rd order Taylor polynom) !
!              at the point (cartx,carty) -> in cube CARTESIAN coordinates          !
!                                                                                   !
! INPUT: value  ...  point value (calculated here by recons_val_cart)               !
!        cell_value ...  tracer value (in the cell center) of the cell              !
!        local_min ...  minmal value in the patch                                   !
!        local_max ...  maximal value in the patch                                  !
! INPUT/OUTPUT: min_phi ... slope limiter, inout because we go through any possible !
!                           extrema on the cell                                     !
!-----------------------------------------------------------------------------------!
subroutine slopelimiter_val(value, cell_value, local_min, local_max, min_phi)

  implicit none
  real (kind=real_kind), intent(in)    :: value, cell_value
  real (kind=real_kind), intent(in)    :: local_min, local_max
  real (kind=real_kind), intent(inout) :: min_phi

  real (kind=real_kind) :: phi 
  
  phi= 0.0D0
  ! Check against the minimum bound on the reconstruction
  if (value - cell_value > 1.0D-12 * value) then
    phi = (local_max - cell_value) / (value - cell_value)
    if (phi < min_phi) then
      min_phi = phi
    endif
  ! Check against the maximum bound on the reconstruction
  elseif (value - cell_value < -1.0D-12 * value) then
    phi = (local_min - cell_value) / (value - cell_value)    
    if(phi < min_phi) then
      min_phi = phi
    endif
  endif

end subroutine slopelimiter_val
!END SUBROUTINE SLOPELIMITER_VAL----------------------------------------CE-for CSLAM!

#endif

end module cslam_filter_mod
