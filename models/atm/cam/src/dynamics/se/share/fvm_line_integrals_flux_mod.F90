#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!MODULE FVM_LINE_INTEGRALS_MOD--------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! Given four vertices of a simply connected cell (x(i),y(i)), i=1,4 ordered in a    !
! counter clockwise manner, compute coefficients for line integral.                 !               
! This module contains everything  to do that and is base on the weights computation!
! from Peter Lauritzens code, here adapted for HOMME                                ! 
!                                                                                   !
!-----------------------------------------------------------------------------------!
module fvm_line_integrals_flux_mod

  use kinds, only               : int_kind, real_kind
  use dimensions_mod, only      : nc, nhe, ngpc
  use parallel_mod, only : abortmp

  implicit none
  private
  real (kind=real_kind),parameter, public   :: bignum = 1.0D20
  real (kind=real_kind),parameter, public   :: tiny   = 1.0D-12
  real (kind=real_kind),parameter           :: fuzzy_width = 10.0*tiny
  ! ALLGAUSS: calculate all line integrals with Gaussian quadrature
  logical                                   :: ALLGAUSS=.FALSE.
  logical :: debugon=.FALSE.
  
  public :: compute_weights_xflux, compute_weights_yflux
contains
! ----------------------------------------------------------------------------------!
!SUBROUTINE COMPUTE_WEIGHTS-----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: computes the coefficients for line integral, set up everything to use!
!              the existing subroutines for the extended (incl. halo zone) element  !
!                                                                                   !
! CALLS: compute_weights_cell, getdep_cellboundaries                                !
! INPUT:  fvm  ... structure, see fvm_control_volume_mod.F90                    ! 
!         nreconstruction ... choose the reconstruction: 1-linear, 3-quadratic      !
!                             6-cubic                                               !
! OUTPUT: weights_all ... value of the weights, number depends of the different     !
!                         intersection/overlap areas to the Eulerian cell           !
!         weights_eul_index_all...correspoding Eulerian cell index of the overlap   !
!         weights_lgr_index_all...which arrival cell is conntected with the         !
!                                 departure cell                                    !
!         jall...number of intersections of an element                              !
!-----------------------------------------------------------------------------------!
subroutine compute_weights_xflux(fvm,nreconstruction,weights_all,weights_eul_index_all, &
                                           weights_lgr_index_all,klev,jall, debug)  
  use fvm_control_volume_mod, only:  fvm_struct                                         
  use coordinate_systems_mod,  only :  cartesian2D_t, spherical_polar_t, &
                                       cart2cubedspherexy, spherical_to_cart
  use physical_constants, only : DD_PI
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest

  implicit none
  type (fvm_struct), intent(inout)                                :: fvm
  integer (kind=int_kind), intent(in)                            :: nreconstruction
  ! arrays for collecting cell data
  real (kind=real_kind),dimension(10*(nc+2*nhe)*(nc+2*nhe),nreconstruction), &
                                                intent(out) :: weights_all
  integer (kind=int_kind), dimension(10*(nc+2*nhe)*(nc+2*nhe),2), &
                                                intent(out) :: weights_eul_index_all
  integer (kind=int_kind), dimension(10*(nc+2*nhe)*(nc+2*nhe),2), &
                                                intent(out) :: weights_lgr_index_all
  integer (kind=int_kind), intent(in)                       :: klev
  integer (kind=int_kind), intent(out)                      :: jall
  logical, intent(in)                                     :: debug
  ! local workspace
  ! max number of line segments is:
  ! (number of longitudes)*(max average number of crossings per line segment = 3)*ncube*2
  !

  integer (kind=int_kind)                     :: jx,jy
  integer                                     :: jx_min, jx_max, jy_min, jy_max
  integer                                     :: jx_min1, jx_max1, jy_min1, jy_max1
  integer                                     :: jx_min2, jx_max2, jy_min2, jy_max2
  logical                                     :: swap1, swap2
  
  integer (kind=int_kind)                     :: i, jtmp
  
  type (cartesian2D_t)                        :: dcart(nc+1,nc+1)       ! Cartesian coordinates 
  
  real (kind=real_kind), dimension(0:5)       :: xcell,ycell
  real (kind=real_kind), dimension(0:5)       :: xcell2,ycell2
  logical                                     :: twofluxareas=.FALSE., zeroflux=.FALSE.
  integer (kind=int_kind)                     :: fluxsign1, fluxsign2 
  
  
  integer (kind=int_kind)                     :: inttmp
  real (kind=real_kind)                       :: tmp
  logical                                     :: swap
  ! for Gaussian quadrature
  real (kind=real_kind), dimension(ngpc)      :: gsweights, gspts
  ! weight-variables for individual cells
  integer (kind=int_kind) :: jmax_segments_cell
  real (kind=real_kind)   , dimension(nhe*50,nreconstruction)   :: weights_cell
  integer (kind=int_kind),  dimension(nhe*50,2)                 :: weights_eul_index_cell
  integer (kind=int_kind)                                       :: jcollect_cell


  jx_min=fvm%jx_min; jx_max=fvm%jx_max; 
  jy_min=fvm%jy_min; jy_max=fvm%jy_max;
  if (fvm%cubeboundary > 0) then
    jx_min1=fvm%jx_min1; jx_max1=fvm%jx_max1; 
    jy_min1=fvm%jy_min1; jy_max1=fvm%jy_max1;
    swap1=fvm%swap1
    if (fvm%cubeboundary > 4) then
      jx_min2=fvm%jx_min2; jx_max2=fvm%jx_max2;
      jy_min2=fvm%jy_min2; jy_max2=fvm%jy_max2;
      swap2=fvm%swap2
    endif
  endif

  jmax_segments_cell = nhe*50
  call gauss_points(ngpc,gsweights,gspts)
  tmp =0.0D0
  jall=1
  ! 
  ! FYI: THESE is CACULATED twice in (xlfux and yflux)
  ! calculate xy Cartesian on the cube of departure points on the corresponding face  
  do jy=1,nc+1
     do jx=1,nc+1               
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%faceno,dcart(jx,jy))              
     end do
  end do
  
  do jy=1, nc
    do jx=0, nc
      !call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)    
      ! xflux: 
      xcell(1)=fvm%acartx(jx+1); ycell(1)=fvm%acarty(jy)
      xcell(3)=fvm%acartx(jx+1); ycell(3)=fvm%acarty(jy+1)
      xcell(2)=dcart(jx+1,jy+1)%x; ycell(2)=dcart(jx+1,jy+1)%y
      xcell(0)=dcart(jx+1,jy)%x; ycell(0)=dcart(jx+1,jy)%y
      
     
      
      if ((debug) .and. (jx==2) .and. (jy==2)) then
        write(*,*) 'eul', fvm%acartx(jx), fvm%acarty(jy) 
        write(*,*) 'eul', fvm%acartx(jx), fvm%acarty(jy+1)
        write(*,*) 'eul', fvm%acartx(jx+1), fvm%acarty(jy+1)
        write(*,*) 'eul', fvm%acartx(jx+1), fvm%acarty(jy)
        write(*,*)
        write(*,*) 'dep', dcart(jx,jy)%x, dcart(jx,jy)%y
        write(*,*) 'dep', dcart(jx,jy+1)%x, dcart(jx,jy+1)%y
        write(*,*) 'dep', dcart(jx+1,jy+1)%x, dcart(jx+1,jy+1)%y
        write(*,*) 'dep', dcart(jx+1,jy)%x, dcart(jx+1,jy)%y
         debugon=.False.
      else
         debugon=.FALSE.
      endif
      
      call makefluxarea(xcell,ycell,xcell2,ycell2,.FALSE.,twofluxareas,zeroflux,fluxsign1, fluxsign2)
      
      if (.not. zeroflux) then
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
             fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
           
        if (jcollect_cell>0) then
          weights_all(jall:jall+jcollect_cell-1,:) = fluxsign1*weights_cell(1:jcollect_cell,:)
          weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                        weights_eul_index_cell(1:jcollect_cell,:)
          weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
          weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
          jall = jall+jcollect_cell          
        endif
        
        ! only if the flow gives us hourglass flux area
        if (twofluxareas) then
          call compute_weights_cell(xcell2,ycell2,jx,jy,nreconstruction,&
               fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
               tmp,ngpc,gsweights,gspts,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
           
          if (jcollect_cell>0) then
            weights_all(jall:jall+jcollect_cell-1,:) = fluxsign2*weights_cell(1:jcollect_cell,:)
            weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                          weights_eul_index_cell(1:jcollect_cell,:)
            weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
            weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
            jall = jall+jcollect_cell          
          endif
        endif
      end if  ! end zeroflux
      
    end do
  end do

  !WEST SIDE
  if (fvm%cubeboundary == west) then
    jx=1
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(west),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
           fvm%nbrsface(west),dcart(jx+1,jy))                 
    end do
    do jy=1,nc
!       call getdep_cellboundaries(xcell,ycell,jx,jy,fvm%nbrsface(west),fvm%dsphere) 
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
          
      if(swap1) then  !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else  
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)+nhe-1
        end do
      end if
      if (fvm%faceno==6) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,2)=jy_max1-jy_min1-weights_eul_index_cell(i,2)-nhe-nhe+1
          end do
        end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  endif
  !EAST SIDE
  if (fvm%cubeboundary == east) then
    jx=nc
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(east),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
           fvm%nbrsface(east),dcart(jx+1,jy))                 
    end do
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=jy_max1-jy_min1-weights_eul_index_cell(i,2)-nhe-nhe+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  endif  
  !NORTH SIDE 
  if (fvm%cubeboundary == north) then
    jy=nc
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(north),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
           fvm%nbrsface(north),dcart(jx,jy+1))                 
    end do
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  
      end if
      !I have to correct the number
      if (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==3) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1))-nhe-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  end if

  !SOUTH SIDE
  if (fvm%cubeboundary == south) then
    jy=1
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(south),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
           fvm%nbrsface(south),dcart(jx,jy+1))                 
    end do
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  
      end if
      !I have to correct the number   
      if  (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
        end do
      end if
      if  (fvm%faceno==3) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if  (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if  (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  end if

  !SOUTHWEST Corner
  if (fvm%cubeboundary == swest) then
    jy=1
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(south),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
           fvm%nbrsface(south),dcart(jx,jy+1))                 
    end do
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if ((fvm%faceno==3) .OR. (fvm%faceno==5)) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
          weights_eul_index_cell(i,2)=(weights_eul_index_cell(i,2)+nhe-1)
        end do
      end if
      if (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
        end do
      end if
      if (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do

    jx=1
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(west),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
           fvm%nbrsface(west),dcart(jx+1,jy))                 
    end do
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap2) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min2,jx_min2,nreconstruction,&
           fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)+nhe-1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2)+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell (1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                          weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  endif

  ! SOUTHEAST Corner
  if (fvm%cubeboundary == seast) then
    jy=1
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(south),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
           fvm%nbrsface(south),dcart(jx,jy+1))                 
    end do
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number   
      if  (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
        end do
      end if
      if  (fvm%faceno==3) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if  (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if  (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do

    jx=nc
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(east),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
           fvm%nbrsface(east),dcart(jx+1,jy))                 
    end do
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
      if(swap2) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min2,jx_min2,nreconstruction,&
           fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2)+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  endif
  
  !NORTHEAST Corner
  if (fvm%cubeboundary == neast) then
    jy=nc
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(north),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
           fvm%nbrsface(north),dcart(jx,jy+1))                 
    end do
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=(weights_eul_index_cell(i,2))-nhe+1
        end do
      end if
      if (fvm%faceno==3) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
    
    jx=nc
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(east),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
           fvm%nbrsface(east),dcart(jx+1,jy))                 
    end do
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap2) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min2,jx_min2,nreconstruction,&
           fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=jy_max2-jy_min2-weights_eul_index_cell(i,2)-nhe-nhe+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)-nhe+1
        end do
      end if
    
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  endif

  !NORTH WEST CORNER 
  if (fvm%cubeboundary == nwest) then
    jy=nc
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(north),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
           fvm%nbrsface(north),dcart(jx,jy+1))                 
    end do
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==3) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1))+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  
    jx=1
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(west),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
           fvm%nbrsface(west),dcart(jx+1,jy))                 
    end do
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap2) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min2,jx_min2,nreconstruction,&
           fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)+nhe-1
        end do
      end if
      if (fvm%faceno==6) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,2)=jy_max2-jy_min2-weights_eul_index_cell(i,2)-nhe-nhe+1
          end do
        end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  endif
  jall=jall-1
end subroutine compute_weights_xflux  



!!!!!!!!!!!!!!! all for make flux area BEGIN

 subroutine makefluxarea(xdep,ydep,xdep2,ydep2,l_yflux,twofluxareas,zeroflux,fluxsign1, fluxsign2)
   IMPLICIT NONE
   !
   ! orientation: South to North / West to East
   !
   !
   ! NOTE: xdep input only uses index (0:3) - xdep is used as a placeholder
   !       xdep output is the full departure area (0:5)
   !       similarly for ydep, xdep2, ydep2
   !
   !
   ! for xflux (l_yflux=.false.): flux-edge points are xydep(1) and xydep(3) where
   !                              xy(1) is the lower one
   !                              departure point for xydep(1) is xydep(0)
   !                              departure point for xydep(3) is xydep(2)
   !
   !     3
   !     |  
   !  2* |
   !     |
   !     |
   !     |
   ! 0*  1
   !
   !
   ! for yflux (l_yflux=.true.): flux-edge points are xydep(1) and xydep(3) where
   !                             xy(1) is the right-most one
   !                             departure point for xydep(1) is xydep(0)
   !                             departure point for xydep(3) is xydep(2)
   !
   !     3--------1
   ! 
   !             0*
   !     2*
   !
   real (kind=real_kind), dimension(0:5), intent(INOUT) :: xdep,ydep,xdep2,ydep2
   !
   ! l_yflux = .true.  for face with constant beta
   ! l_yflux = .false. for face with constant alpha
   !  
   logical, intent(in)  :: l_yflux
   !
   ! The upstream flux region is two triangles or non-convex quadrilateral
   ! 
   logical, intent(out)    :: twofluxareas,zeroflux  
   !
   ! +1 in-flow and -1 for out-flow
   !
   integer (kind=int_kind), intent(out) :: fluxsign1, fluxsign2 
   real (kind=real_kind)   :: xorigo,yorigo
   !
   ! flux-areas are defined as "departure cells"
   !
   real (kind=real_kind), dimension(0:3) :: x,y,xtmp,ytmp
   real (kind=real_kind) :: a,b,y_intersect !to define line and compute crossing
   real (kind=real_kind) :: area !just for DBG
   integer (kind=int_kind) :: i
   !
   real (kind=real_kind), dimension(0:5) :: xdep_tmp,ydep_tmp,xdep2_tmp,ydep2_tmp

   twofluxareas = .false.
   zeroflux = .false.

   if (l_yflux) then
      xorigo = xdep(1); yorigo = ydep(1)
      call translate_then_rotate(xdep(0:3),ydep(0:3),x,y,0,3,xorigo,yorigo)
   else
      x=xdep(0:3)
      y=ydep(0:3)
   end if


   if (abs(x(0)-x(1))<tiny.and.abs(x(2)-x(3))<tiny) then
      !
      ! zero flux area
      !
      zeroflux = .true.
      xdep  = bignum;ydep  = bignum
      xdep2 = bignum;ydep2 = bignum
      return
   end if


   if (abs(x(0)-x(1))<tiny) then
      if (x(2)<x(3)) then
         !WRITE(*,*) "epsilon case 1A"
         !
         ! epsilon cases
         !
         !                 3
         !           2*    |
         !                 |
         !                 |
         !                 |
         !                 |
         !                0/1
         !       
         fluxsign1 = -1
         !
         ! define flux-area
         !
         xdep_tmp(1) = x(0); ydep_tmp(1) = y(0)
         xdep_tmp(2) = x(2); ydep_tmp(2) = y(2)
         xdep_tmp(3) = x(3); ydep_tmp(3) = y(3)
         xdep_tmp(4) = x(3); ydep_tmp(4) = y(3)
      else
         !WRITE(*,*) "epsilon case 1B"
         !
         ! epsilon cases
         !
         !                 3   *2
         !                 |
         !                 |
         !                 |
         !                 |
         !                 |
         !                0/1
         !       
         fluxsign1 = 1
         !
         ! define flux-area
         !
         xdep_tmp(1) = x(0); ydep_tmp(1) = y(0)
         xdep_tmp(2) = x(3); ydep_tmp(2) = y(3)
         xdep_tmp(3) = x(2); ydep_tmp(3) = y(2)
         xdep_tmp(4) = x(2); ydep_tmp(4) = y(2)
      end if
   else if (abs(x(2)-x(3))<tiny) then
      if (x(0)<x(1)) then
         !WRITE(*,*) "epsilon case 2A"
         !
         ! epsilon cases
         !
         !                2/3
         !                 |
         !                 |
         !                 |
         !                 |
         !                 |
         !            0*   1
         !       
         fluxsign1 = -1
         !
         ! define flux-area
         !
         xdep_tmp(1) = x(0); ydep_tmp(1) = y(0)
         xdep_tmp(2) = x(2); ydep_tmp(2) = y(2)
         xdep_tmp(3) = x(1); ydep_tmp(3) = y(1)
         xdep_tmp(4) = x(1); ydep_tmp(4) = y(1)
      else
         !WRITE(*,*) "epsilon case 2B"
         !
         ! epsilon cases
         !
         !                2/3
         !                 |
         !                 |
         !                 |
         !                 |
         !                 |
         !                 1     *0
         !       
         fluxsign1 = 1
         !
         ! define flux-area
         !
         xdep_tmp(1) = x(1); ydep_tmp(1) = y(1)
         xdep_tmp(2) = x(2); ydep_tmp(2) = y(2)
         xdep_tmp(3) = x(0); ydep_tmp(3) = y(0)
         xdep_tmp(4) = x(0); ydep_tmp(4) = y(0)
      end if
   else
      !
      ! CASE ALPHA 1
      !
      !                 3
      !           2*    |
      !                 |
      !                 |
      !                 |
      !                 |
      !                 1
      !        0*
      !
      if ((x(0) < x(1)) .and. (x(2) < x(1))) then
         fluxsign1 = -1
         !
         ! define flux-area
         !
         xdep_tmp(1) = x(0); ydep_tmp(1) = y(0)
         xdep_tmp(2) = x(2); ydep_tmp(2) = y(2)
         xdep_tmp(3) = x(3); ydep_tmp(3) = y(3)
         xdep_tmp(4) = x(1); ydep_tmp(4) = y(1)
         !
         !WRITE(*,*) "CASE 1: simple quadrilateral flux; fluxsign1 ",fluxsign1
      else if ((x(0) > x(1)) .and. (x(2) > x(1))) then
         !
         ! CASE ALPHA 2
         !
         !     3
         !     |    2*
         !     |
         !     |
         !     |
         !     |      0*
         !     1
         !               
         fluxsign1 =  1       
         !
         ! define flux-area
         !
         xdep_tmp(1) = x(1); ydep_tmp(1) = y(1)
         xdep_tmp(2) = x(3); ydep_tmp(2) = y(3)
         xdep_tmp(3) = x(2); ydep_tmp(3) = y(2)
         xdep_tmp(4) = x(0); ydep_tmp(4) = y(0)
         !WRITE(*,*) "CASE 2: simple quadrilateral flux; fluxsign1 ",fluxsign1
      else if ((x(0) < x(1)) .and. (x(2) > x(1))) then
         !
         ! CASE ALPHA 3
         !
         ! "hour glass" flux area or non-convex quadrilateral with x(0) to the left and x(2) to the right
         !
         ! compute crossing
         !
         a = (y(2)-y(0))/(x(2)-x(0)) !slope
         b = y(0)-a*x(0) 
         y_intersect = a*x(1)+b                    
         if (abs(y_intersect-y(1))<tiny) y_intersect = y(1)
         if (abs(y_intersect-y(3))<tiny) y_intersect = y(3)
         !
         !
         ! CASE ALPHA 3A
         !
         !     3
         !     |    
         !     |    2*
         !     |    /  
         !     |   /
         !     |  /    
         !     1 /
         !      /
         !     /
         !    /
         !  0*
         !
         if (y_intersect < y(1)) then
            twofluxareas = .true.
            !
            ! break into two triangles since quadrilateral is non-convex
            ! (topo search in CAM supports non-convex shapes but not the current search in HOMME)
            !
            fluxsign1 = 1
            !
            ! define flux-area 1
            !
            xdep_tmp(1) = x(0); ydep_tmp(1) = y(0)
            xdep_tmp(2) = x(1); ydep_tmp(2) = y(1)
            xdep_tmp(3) = x(2); ydep_tmp(3) = y(2)
            xdep_tmp(4) = x(2); ydep_tmp(4) = y(2)
            !
            ! second triangle
            !
            fluxsign2 = 1
            !
            ! define flux-area 2
            !
            xdep2_tmp(1) = x(1); ydep2_tmp(1) = y(1)
            xdep2_tmp(2) = x(3); ydep2_tmp(2) = y(3)
            xdep2_tmp(3) = x(2); ydep2_tmp(3) = y(2)
            xdep2_tmp(4) = x(2); ydep2_tmp(4) = y(2)
            
            !WRITE(*,*) "CASE 3A: fluxsign1s ",fluxsign1,fluxsign2
         else if (y_intersect > y(3)) then
            !
            ! CASE ALPHA 3C
            !
            !      2*
            !      /
            !     /
            !    / 
            !   /
            !  /  3
            ! 0*  |    
            !     |
            !     |
            !     |
            !     |
            !     1
            !                  
            !
            ! break into two triangles since quadrilateral is non-convex
            !
            twofluxareas = .true.
            fluxsign1 = -1
            !
            ! define flux-area 1
            !
            xdep_tmp(1) = x(1); ydep_tmp(1) = y(1)
            xdep_tmp(2) = x(0); ydep_tmp(2) = y(0)
            xdep_tmp(3) = x(3); ydep_tmp(3) = y(3)
            xdep_tmp(4) = x(3); ydep_tmp(4) = y(3)
            !
            ! second triangle
            !
            fluxsign2 = -1
            !
            ! define flux-area 2
            !
            xdep2_tmp(1) = x(0); ydep2_tmp(1) = y(0)
            xdep2_tmp(2) = x(2); ydep2_tmp(2) = y(2)
            xdep2_tmp(3) = x(3); ydep2_tmp(3) = y(3)
            xdep2_tmp(4) = x(3); ydep2_tmp(4) = y(3)
            
            !WRITE(*,*) "CASE 3C: fluxsign1s ",fluxsign1,fluxsign2
         else
            !
            ! CASE ALPHA 3B
            !
            !
            !     3
            !     |    *2
            !     |   /
            !     |  /
            !     | /
            !     |/
            !     |
            !    /|             
            !   / 1    
            ! 0*
            !
            !
            ! break into two triangles since quadrilateral is non-convex
            !
            twofluxareas = .true.
            fluxsign1 = -1
            !
            ! define flux-area 1
            !
            xdep_tmp(1) = x(0); ydep_tmp(1) = y(0)
            xdep_tmp(2) = x(1); ydep_tmp(2) = y_intersect
            xdep_tmp(3) = x(1); ydep_tmp(3) = y(1)
            xdep_tmp(4) = x(1); ydep_tmp(4) = y(1)
            !
            ! second triangle
            !
            fluxsign2 = 1
            !
            ! define flux-area 2
            !
            xdep2_tmp(1) = x(1); ydep2_tmp(1) = y_intersect
            xdep2_tmp(2) = x(3); ydep2_tmp(2) = y(3)
            xdep2_tmp(3) = x(2); ydep2_tmp(3) = y(2)
            xdep2_tmp(4) = x(2); ydep2_tmp(4) = y(2)
            
            !WRITE(*,*) "CASE 3B: fluxsign1s ",fluxsign1,fluxsign2
         end if
      else 
         !
         ! entering here means that ((x(0) >= x(1)) .and. (x(2) <= x(1))) then
         !
         !
         ! CASE ALPHA 4
         !
         ! "hour glass" flux area or non-convex quadrilateral with x(0) 
         ! to the right and x(2) to the left
         !
         !
         ! compute crossing
         !
         a = (y(0)-y(2))/(x(0)-x(2)) !slope
         b = y(0)-a*x(0) 
         y_intersect = a*x(1)+b                    
         if (abs(y_intersect-y(1))<tiny) y_intersect = y(1)
         if (abs(y_intersect-y(3))<tiny) y_intersect = y(3)         
         !
         !WRITE(*,*) "y_intersect",y_intersect
         !
         if (Y_intersect < y(1)) then
            !
            !
            ! CASE ALPHA 4A
            !             
            !
            !     3
            !     |    
            !     |
            !2*   |
            ! \   |
            !  \  |
            !   \ 1 
            !     \
            !      \
            !      0*
            !
            !
            ! break into two triangles since quadrilateral is non-convex
            !
            twofluxareas = .true.
            fluxsign1 = -1
            !
            ! define flux-area 1
            !
            xdep_tmp(1) = x(0); ydep_tmp(1) = y(0)
            xdep_tmp(2) = x(2); ydep_tmp(2) = y(2)
            xdep_tmp(3) = x(1); ydep_tmp(3) = y(1)
            xdep_tmp(4) = x(1); ydep_tmp(4) = y(1)
            !
            ! second triangle
            !
            fluxsign2 = -1
            !
            ! define flux-area 2
            !
            xdep2_tmp(1) = x(2); ydep2_tmp(1) = y(2)
            xdep2_tmp(2) = x(3); ydep2_tmp(2) = y(3)
            xdep2_tmp(3) = x(1); ydep2_tmp(3) = y(1)
            xdep2_tmp(4) = x(1); ydep2_tmp(4) = y(1)
            
            !WRITE(*,*) "CASE 4A: fluxsign1s ",fluxsign1,fluxsign2
         else if (y_intersect > y(3)) then
            !
            ! CASE ALPHA 4C
            !
            !
            !  2*
            !    \
            !     \
            !      \
            !     3  
            !     |    *0
            !     |
            !     |
            !     1
            !      
            ! break into two triangles since quadrilateral is non-convex
            !
            twofluxareas = .true.
            fluxsign1 = 1
            !
            ! define flux-area 1
            !
            xdep_tmp(1) = x(1); ydep_tmp(1) = y(1)
            xdep_tmp(2) = x(3); ydep_tmp(2) = y(3)
            xdep_tmp(3) = x(0); ydep_tmp(3) = y(0)
            xdep_tmp(4) = x(0); ydep_tmp(4) = y(0)
            !
            ! second triangle
            !
            fluxsign2 = 1
            !
            ! define flux-area 2
            !
            xdep2_tmp(1) = x(3); ydep2_tmp(1) = y(3)
            xdep2_tmp(2) = x(2); ydep2_tmp(2) = y(2)
            xdep2_tmp(3) = x(0); ydep2_tmp(3) = y(0)
            xdep2_tmp(4) = x(0); ydep2_tmp(4) = y(0)
            !
            !WRITE(*,*) "CASE 4C: fluxsign1s ",fluxsign1,fluxsign2
         else
            !
            ! CASE ALPHA 4B
            !
            !     3
            ! *2  |
            !   \ |
            !    \|
            !     |
            !     |\
            !     | \
            !     |  *0           
            !     1    
            ! 
            !
            !
            ! break into two triangles since quadrilateral is non-convex
            !
            twofluxareas = .true.
            fluxsign1 = 1
            !
            ! define flux-area 1
            !
            xdep_tmp(1) = x(1); ydep_tmp(1) = y(1)
            xdep_tmp(2) = x(1); ydep_tmp(2) = y_intersect
            xdep_tmp(3) = x(0); ydep_tmp(3) = y(0)
            xdep_tmp(4) = x(0); ydep_tmp(4) = y(0)
            !
            ! second triangle
            !
            fluxsign2 = -1
            !
            ! define flux-area 2
            !
            xdep2_tmp(1) = x(1); ydep2_tmp(1) = y_intersect
            xdep2_tmp(2) = x(2); ydep2_tmp(2) = y(2)
            xdep2_tmp(3) = x(3); ydep2_tmp(3) = y(3)
            xdep2_tmp(4) = x(3); ydep2_tmp(4) = y(3)
            !
            !WRITE(*,*) "CASE 4B: fluxsign1s ",fluxsign1,fluxsign2
         end if
      end if
   end if
   !
   !************************************************************************
   !
   ! DONE DEFINING "DEPARTURE AREAS"
   !
   ! compute fluxes
   !
   !************************************************************************
   !

   if (l_yflux) then
      call rotate_then_translate(xdep_tmp(1:4),ydep_tmp(1:4),xdep(1:4),ydep(1:4),&
           1,4,xorigo,yorigo)
      !          
      !          xdep_tmp(1:4) = xdep(1:4); ydep_tmp(1:4) = ydep(1:4)
      !          !
      !          ! make sure xdep is lower left corner
      !          !  
      !          
      !          xdep(1) = xdep_tmp(2); ydep(1) = ydep_tmp(2); 
      !          xdep(2) = xdep_tmp(3); ydep(2) = ydep_tmp(3); 
      !          xdep(3) = xdep_tmp(4); ydep(3) = ydep_tmp(4); 
      !          xdep(4) = xdep_tmp(1); ydep(4) = ydep_tmp(1); 
   else
      xdep(1:4) = xdep_tmp(1:4); ydep(1:4) = ydep_tmp(1:4);
   end if
   xdep(5) = xdep(1); ydep(5) = ydep(1)
   xdep(0) = xdep(4); ydep(0) = ydep(4)      

   if (twofluxareas) then
      if (l_yflux) then
         call rotate_then_translate(xdep2_tmp(1:4),ydep2_tmp(1:4),&
              xdep2(1:4),ydep2(1:4),1,4,xorigo,yorigo)
      else
         xdep2(1:4) = xdep2_tmp(1:4); ydep2(1:4) = ydep2_tmp(1:4);
      end if
      

      xdep2(5) = xdep2(1); ydep2(5) = ydep2(1)
      xdep2(0) = xdep2(4); ydep2(0) = ydep2(4)      
   else
      xdep2 = bignum;ydep2 = bignum;
   end if
 end subroutine makefluxarea


subroutine translate_then_rotate(x_in,y_in,x,y,imin,imax,xorigo,yorigo)
  !
  ! translates x_in and y_in to (x_in(3),y_in(3)) and rotates data -90 degrees about xorigo,yorigo
  !
  IMPLICIT NONE
  integer (kind=int_kind)                    , intent(in)  :: imin,imax
  real (kind=real_kind), dimension(imin:imax), intent(in)  :: x_in,y_in
  real (kind=real_kind)                                    :: xorigo, yorigo
  real (kind=real_kind), dimension(imin:imax), intent(out) :: x,y
  !
  !
  !
  real (kind=real_kind), dimension(imin:imax) :: xtmp,ytmp
  integer (kind=int_kind) :: i
  !
  !
  ! move origo to xorigo,yorigo
  !
  do i=imin,imax
     xtmp(i) = x_in(i)-xorigo
     ytmp(i) = y_in(i)-yorigo
  end do
  !
  ! rotate +90 degrees
  !
  x =  ytmp
  y = -xtmp
end subroutine translate_then_rotate

subroutine rotate_then_translate(x_in,y_in,x,y,imin,imax,xorigo,yorigo)
  !
  ! reverse of translate_then_rotate subroutine
  !
  IMPLICIT NONE
  integer (kind=int_kind)                    , intent(in)  :: imin,imax
  real (kind=real_kind), dimension(imin:imax), intent(in)  :: x_in,y_in
  real (kind=real_kind)                                    :: xorigo, yorigo
  real (kind=real_kind), dimension(imin:imax), intent(out) :: x,y
  !
  !
  !
  real (kind=real_kind), dimension(imin:imax) :: xtmp,ytmp
  integer (kind=int_kind) :: i
  !
  !
  ! rotate -90 degrees
  !
  xtmp =  -y_in
  ytmp =   x_in
  !
  ! move origo
  !
  do i=imin,imax
     x(i) = xtmp(i)+xorigo
     y(i) = ytmp(i)+yorigo
  end do
end subroutine rotate_then_translate

!!!!!!!!!!!!!!! all for make flux area END

!!!! Y flux calculation
subroutine compute_weights_yflux(fvm,nreconstruction,weights_all,weights_eul_index_all, &
                                           weights_lgr_index_all,klev,jall,debug)  
  use fvm_control_volume_mod, only:  fvm_struct                                         
  use coordinate_systems_mod,  only :  cartesian2D_t, spherical_polar_t, &
                                       cart2cubedspherexy, spherical_to_cart
  use physical_constants, only : DD_PI
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest

  implicit none
  type (fvm_struct), intent(inout)                                :: fvm
  integer (kind=int_kind), intent(in)                            :: nreconstruction
  ! arrays for collecting cell data
  real (kind=real_kind),dimension(10*(nc+2*nhe)*(nc+2*nhe),nreconstruction), &
                                                intent(out) :: weights_all
  integer (kind=int_kind), dimension(10*(nc+2*nhe)*(nc+2*nhe),2), &
                                                intent(out) :: weights_eul_index_all
  integer (kind=int_kind), dimension(10*(nc+2*nhe)*(nc+2*nhe),2), &
                                                intent(out) :: weights_lgr_index_all
  integer (kind=int_kind), intent(in)                       :: klev
  integer (kind=int_kind), intent(out)                      :: jall
  logical, intent(in)                                     :: debug

  ! local workspace
  ! max number of line segments is:
  ! (number of longitudes)*(max average number of crossings per line segment = 3)*ncube*2
  !

  integer (kind=int_kind)                     :: jx,jy
  integer                                     :: jx_min, jx_max, jy_min, jy_max
  integer                                     :: jx_min1, jx_max1, jy_min1, jy_max1
  integer                                     :: jx_min2, jx_max2, jy_min2, jy_max2
  logical                                     :: swap1, swap2
  
  integer (kind=int_kind)                     :: i, jtmp
  
  type (cartesian2D_t)                        :: dcart(nc+1,nc+1)       ! Cartesian coordinates 
  
  real (kind=real_kind), dimension(0:5)       :: xcell,ycell
  real (kind=real_kind), dimension(0:5)       :: xcell2,ycell2
  logical                                     :: twofluxareas=.FALSE., zeroflux=.FALSE.
  integer (kind=int_kind)                     :: fluxsign1, fluxsign2 
  
  
  integer (kind=int_kind)                     :: inttmp
  real (kind=real_kind)                       :: tmp
  logical                                     :: swap
  ! for Gaussian quadrature
  real (kind=real_kind), dimension(ngpc)      :: gsweights, gspts
  ! weight-variables for individual cells
  integer (kind=int_kind) :: jmax_segments_cell
  real (kind=real_kind)   , dimension(nhe*50,nreconstruction)   :: weights_cell
  integer (kind=int_kind),  dimension(nhe*50,2)                 :: weights_eul_index_cell
  integer (kind=int_kind)                                       :: jcollect_cell

!   debugon=debug
  
  jx_min=fvm%jx_min; jx_max=fvm%jx_max; 
  jy_min=fvm%jy_min; jy_max=fvm%jy_max;
  if (fvm%cubeboundary > 0) then
    jx_min1=fvm%jx_min1; jx_max1=fvm%jx_max1; 
    jy_min1=fvm%jy_min1; jy_max1=fvm%jy_max1;
    swap1=fvm%swap1
    if (fvm%cubeboundary > 4) then
      jx_min2=fvm%jx_min2; jx_max2=fvm%jx_max2;
      jy_min2=fvm%jy_min2; jy_max2=fvm%jy_max2;
      swap2=fvm%swap2
    endif
  endif

  jmax_segments_cell = nhe*50
  call gauss_points(ngpc,gsweights,gspts)
  tmp =0.0D0
  jall=1
  ! 
  ! FYI: THESE is CACULATED twice in (xlfux and yflux)
  ! calculate xy Cartesian on the cube of departure points on the corresponding face  
  do jy=1,nc+1
     do jx=1,nc+1               
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%faceno,dcart(jx,jy))              
     end do
  end do
  
  do jy=0, nc
    do jx=1, nc
      !call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)    
      ! xflux: 
      xcell(1)=fvm%acartx(jx+1); ycell(1)=fvm%acarty(jy+1)
      xcell(3)=fvm%acartx(jx); ycell(3)=fvm%acarty(jy+1)
      xcell(2)=dcart(jx,jy+1)%x; ycell(2)=dcart(jx,jy+1)%y
      xcell(0)=dcart(jx+1,jy+1)%x; ycell(0)=dcart(jx+1,jy+1)%y
      
      if ((debug) .and. (jx==2) .and. (jy==2)) then
        write(*,*) 'eul', fvm%acartx(jx), fvm%acarty(jy) 
        write(*,*) 'eul', fvm%acartx(jx), fvm%acarty(jy+1)
        write(*,*) 'eul', fvm%acartx(jx+1), fvm%acarty(jy+1)
        write(*,*) 'eul', fvm%acartx(jx+1), fvm%acarty(jy)
        write(*,*)
        write(*,*) 'dep', dcart(jx,jy)%x, dcart(jx,jy)%y
        write(*,*) 'dep', dcart(jx,jy+1)%x, dcart(jx,jy+1)%y
        write(*,*) 'dep', dcart(jx+1,jy+1)%x, dcart(jx+1,jy+1)%y
        write(*,*) 'dep', dcart(jx+1,jy)%x, dcart(jx+1,jy)%y
         debugon=.TRUE.
      else
         debugon=.FALSE.
      endif
      
      
      call makefluxarea(xcell,ycell,xcell2,ycell2,.TRUE.,twofluxareas,zeroflux,fluxsign1, fluxsign2)
      
      if (.not. zeroflux) then
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
             fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
           
        if (jcollect_cell>0) then
          weights_all(jall:jall+jcollect_cell-1,:) = fluxsign1*weights_cell(1:jcollect_cell,:)
          weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                        weights_eul_index_cell(1:jcollect_cell,:)
          weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
          weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
          jall = jall+jcollect_cell          
        endif
        
        ! only if the flow gives us hourglass flux area
        if (twofluxareas) then
          call compute_weights_cell(xcell2,ycell2,jx,jy,nreconstruction,&
               fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
               tmp,ngpc,gsweights,gspts,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
           
          if (jcollect_cell>0) then
            weights_all(jall:jall+jcollect_cell-1,:) = fluxsign2*weights_cell(1:jcollect_cell,:)
            weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                          weights_eul_index_cell(1:jcollect_cell,:)
            weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
            weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
            jall = jall+jcollect_cell          
          endif
        endif
      end if  ! end zeroflux
      
    end do
  end do

  !WEST SIDE
  if (fvm%cubeboundary == west) then
    jx=1
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(west),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
           fvm%nbrsface(west),dcart(jx+1,jy))                 
    end do
    do jy=1,nc
!       call getdep_cellboundaries(xcell,ycell,jx,jy,fvm%nbrsface(west),fvm%dsphere) 
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
          
      if(swap1) then  !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else  
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)+nhe-1
        end do
      end if
      if (fvm%faceno==6) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,2)=jy_max1-jy_min1-weights_eul_index_cell(i,2)-nhe-nhe+1
          end do
        end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  endif
  !EAST SIDE
  if (fvm%cubeboundary == east) then
    jx=nc
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(east),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
           fvm%nbrsface(east),dcart(jx+1,jy))                 
    end do
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=jy_max1-jy_min1-weights_eul_index_cell(i,2)-nhe-nhe+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  endif  
  !NORTH SIDE 
  if (fvm%cubeboundary == north) then
    jy=nc
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(north),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
           fvm%nbrsface(north),dcart(jx,jy+1))                 
    end do
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  
      end if
      !I have to correct the number
      if (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==3) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1))-nhe-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  end if

  !SOUTH SIDE
  if (fvm%cubeboundary == south) then
    jy=1
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(south),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
           fvm%nbrsface(south),dcart(jx,jy+1))                 
    end do
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  
      end if
      !I have to correct the number   
      if  (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
        end do
      end if
      if  (fvm%faceno==3) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if  (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if  (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  end if

  !SOUTHWEST Corner
  if (fvm%cubeboundary == swest) then
    jy=1
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(south),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
           fvm%nbrsface(south),dcart(jx,jy+1))                 
    end do
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if ((fvm%faceno==3) .OR. (fvm%faceno==5)) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
          weights_eul_index_cell(i,2)=(weights_eul_index_cell(i,2)+nhe-1)
        end do
      end if
      if (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
        end do
      end if
      if (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do

    jx=1
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(west),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
           fvm%nbrsface(west),dcart(jx+1,jy))                 
    end do
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap2) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min2,jx_min2,nreconstruction,&
           fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)+nhe-1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2)+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell (1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                          weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  endif

  ! SOUTHEAST Corner
  if (fvm%cubeboundary == seast) then
    jy=1
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(south),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
           fvm%nbrsface(south),dcart(jx,jy+1))                 
    end do
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number   
      if  (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
        end do
      end if
      if  (fvm%faceno==3) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if  (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if  (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do

    jx=nc
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(east),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
           fvm%nbrsface(east),dcart(jx+1,jy))                 
    end do
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
      if(swap2) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min2,jx_min2,nreconstruction,&
           fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2)+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  endif
  
  !NORTHEAST Corner
  if (fvm%cubeboundary == neast) then
    jy=nc
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(north),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
           fvm%nbrsface(north),dcart(jx,jy+1))                 
    end do
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=(weights_eul_index_cell(i,2))-nhe+1
        end do
      end if
      if (fvm%faceno==3) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
    
    jx=nc
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(east),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
           fvm%nbrsface(east),dcart(jx+1,jy))                 
    end do
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap2) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min2,jx_min2,nreconstruction,&
           fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=jy_max2-jy_min2-weights_eul_index_cell(i,2)-nhe-nhe+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)-nhe+1
        end do
      end if
    
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  endif

  !NORTH WEST CORNER 
  if (fvm%cubeboundary == nwest) then
    jy=nc
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(north),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
           fvm%nbrsface(north),dcart(jx,jy+1))                 
    end do
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
           fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==3) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1))+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  
    jx=1
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=1,nc+1
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%nbrsface(west),dcart(jx,jy))  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
           fvm%nbrsface(west),dcart(jx+1,jy))                 
    end do
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap2) then !flip orientation
        call compute_weights_cell(xcell,ycell,jy_min2,jx_min2,nreconstruction,&
           fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)+nhe-1
        end do
      end if
      if (fvm%faceno==6) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,2)=jy_max2-jy_min2-weights_eul_index_cell(i,2)-nhe-nhe+1
          end do
        end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                         weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  endif
  jall=jall-1
end subroutine compute_weights_yflux  




!END SUBROUTINE COMPUTE_WEIGHTS-------------------------------------------CE-for FVM!

! ----------------------------------------------------------------------------------!
!SUBROUTINE GETDEP_CELLBOUNDARIESXYVEC------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 20.March 2012                                            !
! DESCRIPTION: needed to apply the search subroutines                               !
!                                                                                   !
! CALLS: cart2cubedspherexy                                                         !
! INPUT:  jx,jy   ... index of the cell                                             !
!         dcart ... Cartesian coordinates of the cell on the corresponding face     !
! OUTPUT: xcell, ycell ... x and y Cartesian coordinates on the cube of the cell    !
!                          dsphere                                                  !
!-----------------------------------------------------------------------------------!
subroutine getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)
use coordinate_systems_mod, only : cartesian2D_t
  implicit none
  real (kind=real_kind), dimension(0:5), intent(inout)     :: xcell,ycell
  integer (kind=int_kind), intent(in)                      :: jx, jy
  type (cartesian2D_t), intent(in)                         :: dcart(nc+1,nc+1)

   xcell(1) = dcart(jx  ,jy  )%x 
   ycell(1) = dcart(jx  ,jy  )%y
   xcell(2) = dcart(jx  ,jy+1)%x 
   ycell(2) = dcart(jx  ,jy+1)%y
   xcell(3) = dcart(jx+1,jy+1)%x 
   ycell(3) = dcart(jx+1,jy+1)%y
   xcell(4) = dcart(jx+1,jy  )%x 
   ycell(4) = dcart(jx+1,jy  )%y
   xcell(5) = xcell(1)        
   ycell(5) = ycell(1)          
   xcell(0) = xcell(4)        
   ycell(0) = ycell(4)
end subroutine getdep_cellboundariesxyvec
!END SUBROUTINE GETDEP_CELLBOUNDARIESXYVEC--------------------------------CE-for FVM!
    
! ----------------------------------------------------------------------------------!
!NEXT SUBROUTINES ARE TAKEN FROM------------------------------------------CE-for FVM!
! PETER LAURITZENs CODE adapted to use in HOMME                                    !
! ----------------------------------------------------------------------------------!

  subroutine compute_weights_cell(xcell_in,ycell_in,jx,jy,nreconstruction,xgno,ygno,&
       jx_min, jx_max, jy_min, jy_max,tmp,&
       ngauss,gauss_weights,abscissae,weights,weights_eul_index,jcollect,jmax_segments)

    implicit none
    integer (kind=int_kind)                  , intent(in):: nreconstruction, jx,jy,ngauss,jmax_segments
    real (kind=real_kind)   ,  dimension(0:5), intent(in):: xcell_in,ycell_in
    !
    integer (kind=int_kind), intent(in)               :: jx_min, jy_min, jx_max, jy_max
    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in) :: xgno
    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in) :: ygno
    !
    ! for Gaussian quadrature
    !
    real (kind=real_kind), dimension(ngauss), intent(in) :: gauss_weights, abscissae
    !
    ! boundaries of domain
    !
    real (kind=real_kind):: tmp
    !
    ! Number of Eulerian sub-cell integrals for the cell in question
    !
    integer (kind=int_kind), intent(out)       :: jcollect
    !
    ! local workspace
    !
    !
    ! max number of line segments is:
    !
    ! (number of longitudes)*(max average number of crossings per line segment = 3)*ncube*2
    !
    real (kind=real_kind)   ,  &
         dimension(jmax_segments,nreconstruction), intent(out) :: weights
    integer (kind=int_kind),  &
         dimension(jmax_segments,2), intent(out)      :: weights_eul_index
    
    real (kind=real_kind), dimension(0:3) :: x,y
    integer (kind=int_kind),dimension(0:5) :: jx_eul, jy_eul
    integer (kind=int_kind) :: jsegment,i
    !
    ! variables for registering crossings with Eulerian latitudes and longitudes
    !
    integer (kind=int_kind)  :: jcross_lat, iter
    !
    ! max. crossings per side is 2*nhe
    !
    real (kind=real_kind), &
         dimension(8*nhe,2) :: r_cross_lat
    integer (kind=int_kind), &
         dimension(8*nhe,2) :: cross_lat_eul_index
    real (kind=real_kind)   ,  dimension(4) :: xcell,ycell

    integer (kind=int_kind) :: nvertex
       

    nvertex = 4
    xcell = xcell_in(1:nvertex)
    ycell = ycell_in(1:nvertex)

    jsegment   = 0
    weights    = 0.0D0
    jcross_lat = 0
          
    call side_integral(xcell,ycell,4,jsegment,jmax_segments,&
         weights,weights_eul_index,nreconstruction,jx,jy,xgno,ygno,jx_min, jx_max, jy_min, jy_max,&
         ngauss,gauss_weights,abscissae,&
         jcross_lat,r_cross_lat,cross_lat_eul_index)
    !
    !**********************
    ! 
    ! Do inner integrals
    !
    !**********************
    !    
    call compute_inner_line_integrals_lat(r_cross_lat,cross_lat_eul_index,&
         jcross_lat,jsegment,jmax_segments,xgno,jx_min, jx_max, jy_min, jy_max,&
         weights,weights_eul_index,&
         nreconstruction,ngauss,gauss_weights,abscissae)
    
    
    IF (ABS((jcross_lat/2)-DBLE(jcross_lat)/2.0)>tiny) then
      WRITE(*,*) "number of latitude crossings are not even: ABORT",jcross_lat,jx,jy
      STOP
    END IF
    
    !
    ! collect line-segment that reside in the same Eulerian cell
    !
    if (jsegment>0) then
      call collect(weights,weights_eul_index,nreconstruction,jcollect,jsegment,jmax_segments)
      !
      ! DBG
      !
      tmp=0.0
      !       WRITE(*,*) "max area",maxval(weights(1:jcollect,1))
      !       WRITE(*,*) "min area",minval(weights(1:jcollect,1))
      !       stop
      do i=1,jcollect     
        IF (weights(i,1)<0.0) THEN
          IF (weights(i,1)<-1.0E-10) THEN
            WRITE(*,*) "negative cell area",weights(i,1)
            STOP
          END IF
          !           weights(i,2:nreconstruction) = 0.0
        END IF
        !       end do
        
        tmp=tmp+weights(i,1)
      enddo

      IF (abs(tmp)>0.04) THEN
        WRITE(*,*) "sum of weights too large",tmp
        stop
      END IF
      IF (tmp<-1.0E-9) THEN
        WRITE(*,*) "sum of weights is negative - negative area?",tmp,jx,jy
        stop
      END IF
    else
      jcollect = 0
    end if
  end subroutine compute_weights_cell


  !
  !****************************************************************************
  !
  ! organize data and store it
  !
  !****************************************************************************
  !
  subroutine collect(weights,weights_eul_index,nreconstruction,jcollect,jsegment,jmax_segments)
    implicit none
    integer (kind=int_kind),                                  INTENT(IN   ) :: jsegment,jmax_segments
    integer (kind=int_kind)                                 , intent(in)    :: nreconstruction
    real (kind=real_kind)   , dimension(jmax_segments,nreconstruction), intent(inout) :: weights
    integer (kind=int_kind), dimension(jmax_segments,2     ), intent(inout) :: weights_eul_index
    integer (kind=int_kind),                                  INTENT(OUT  ) :: jcollect
    !
    ! local workspace
    !
    integer (kind=int_kind) :: imin, imax, jmin, jmax, i,j,k,h
    logical                 :: ltmp

    real (kind=real_kind)   , dimension(jmax_segments,nreconstruction) :: weights_out
    integer (kind=int_kind), dimension(jmax_segments,2     ) :: weights_eul_index_out

    weights_out           = 0.0D0
    weights_eul_index_out = -100

    imin = MINVAL(weights_eul_index(1:jsegment,1))
    imax = MAXVAL(weights_eul_index(1:jsegment,1))
    jmin = MINVAL(weights_eul_index(1:jsegment,2))
    jmax = MAXVAL(weights_eul_index(1:jsegment,2))

    ltmp = .FALSE.

    jcollect = 1

    do j=jmin,jmax
       do i=imin,imax
          do k=1,jsegment
             if (weights_eul_index(k,1)==i.AND.weights_eul_index(k,2)==j) then
                weights_out(jcollect,1:nreconstruction) = &
                     weights_out(jcollect,1:nreconstruction) + weights(k,1:nreconstruction)
                ltmp = .TRUE.
                h = k
             endif
          enddo
          if (ltmp) then
             weights_eul_index_out(jcollect,:) = weights_eul_index(h,:)
             jcollect = jcollect+1
          endif
          ltmp = .FALSE.
       enddo
    enddo
    jcollect = jcollect-1
    weights           = weights_out
    weights_eul_index = weights_eul_index_out
  end subroutine collect

  !    
  ! compute crossings with Eulerian latitudes and longitudes
  !
  !

  !
  !*****************************************************************************************
  !
  ! 
  !
  !*****************************************************************************************
  !
  subroutine compute_inner_line_integrals_lat(r_cross_lat,cross_lat_eul_index,&
       jcross_lat,jsegment,jmax_segments,xgno,jx_min,jx_max,jy_min, jy_max,weights,weights_eul_index,&
       nreconstruction,ngauss,gauss_weights,abscissae)    
    implicit none
    !
    ! variables for registering crossings with Eulerian latitudes and longitudes
    !
    integer (kind=int_kind),         intent(in):: jcross_lat, jmax_segments,nreconstruction,ngauss
    integer (kind=int_kind),         intent(inout):: jsegment
    !
    ! for Gaussian quadrature
    !
    real (kind=real_kind), dimension(ngauss), intent(in) :: gauss_weights, abscissae
    !
    ! max. crossings per side is 2*nhe
    !
    real (kind=real_kind), &
         dimension(8*nhe,2), intent(in):: r_cross_lat
    integer (kind=int_kind), &
         dimension(8*nhe,2), intent(in):: cross_lat_eul_index
    integer (kind=int_kind), intent(in):: jx_min, jx_max, jy_min, jy_max

    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in)  :: xgno
    real (kind=real_kind)   ,  &
         dimension(jmax_segments,nreconstruction), intent(inout) :: weights
    integer (kind=int_kind),  &
         dimension(jmax_segments,2), intent(inout)               :: weights_eul_index
    real (kind=real_kind)   , dimension(nreconstruction)         :: weights_tmp

    integer (kind=int_kind) :: imin, imax, jmin, jmax, i,j,k, isgn, h, eul_jx, eul_jy
    integer (kind=int_kind) :: idx_start_y,idx_end_y
    logical                 :: ltmp,lcontinue
    real (kind=real_kind), dimension(2)  :: rstart,rend,rend_tmp
    real (kind=real_kind), dimension(2)  :: xseg, yseg
    5   FORMAT(10e14.6)
    lcontinue = .TRUE.
    if (jcross_lat>0) then
       do i=MINVAL(cross_lat_eul_index(1:jcross_lat,2)),MAXVAL(cross_lat_eul_index(1:jcross_lat,2))
          !
          ! find "first" crossing with Eulerian cell i
          !
          do k=1,jcross_lat
             if (cross_lat_eul_index(k,2)==i) exit
          enddo
          do j=k+1,jcross_lat
             !
             ! find "second" crossing with Eulerian cell i
             !
             if (cross_lat_eul_index(j,2)==i) then
                if (r_cross_lat(k,1)<r_cross_lat(j,1)) then
                   rstart = r_cross_lat(k,1:2)
                   rend   = r_cross_lat(j,1:2)
                   imin   = cross_lat_eul_index(k,1)
                   imax   = cross_lat_eul_index(j,1)
                else
                   rstart = r_cross_lat(j,1:2)
                   rend   = r_cross_lat(k,1:2)
                   imin   = cross_lat_eul_index(j,1)
                   imax   = cross_lat_eul_index(k,1)
                endif
                do h=imin,imax
                   if (h==imax) then
                      rend_tmp = rend
                   else
                      rend_tmp(1) = xgno(h+1)
                      rend_tmp(2) = r_cross_lat(k,2)
                   endif
                   xseg(1) = rstart(1)
                   xseg(2) = rend_tmp(1)
                   yseg(1) = rstart(2)
                   yseg(2) = rend_tmp(2)
                   call get_weights_exact(weights_tmp,xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)

!phl                   if (i.LE.nc.AND.i.GE.1.AND.h.LE.nc.AND.h.GE.1) then
                   if (i.LE.jy_max-1.AND.i.GE.jy_min.AND.h.LE.jx_max-1.AND.h.GE.jx_min) then
                      jsegment=jsegment+1
                      weights_eul_index(jsegment,1) = h 
                      weights_eul_index(jsegment,2) = i
                      weights(jsegment,1:nreconstruction) = -weights_tmp
                   endif
                   !
                   ! subtract the same weights on the west side of the line
                   !
                   if (i.LE.jy_max.AND.i.GE.jy_min+1.AND.h.LE.jx_max-1.AND.h.GE.jx_min) then
!phl                   if (i.GE.2.AND.i.LE.nc+1.AND.h.LE.nc.AND.h.GE.1) then
                      jsegment = jsegment+1
                      weights_eul_index(jsegment,1) = h 
                      weights_eul_index(jsegment,2) = i-1
                      weights(jsegment,1:nreconstruction) = weights_tmp
                   endif
                   !
                   ! prepare for next iteration
                   !
                   rstart = rend_tmp
                enddo
             endif
          enddo
       enddo
    endif
  end subroutine compute_inner_line_integrals_lat

  !
  ! line integral from (a1_in,a2_in) to (b1_in,b2_in)
  ! If line is coniciding with an Eulerian longitude or latitude the routine
  ! needs to know where an adjacent side is located to determine which
  ! reconstruction must be used. therefore (c1,c2) is passed to the routine
  !
  !   

  subroutine side_integral(&
       x_in,y_in,nvertex,jsegment,jmax_segments,&
       weights,weights_eul_index,nreconstruction,jx,jy,xgno,ygno,jx_min,jx_max,jy_min,jy_max,&
       ngauss,gauss_weights,abscissae,&!)!phl add jx_min etc.
       jcross_lat,r_cross_lat,cross_lat_eul_index)
    implicit none
    integer (kind=int_kind),            intent(in)    :: nreconstruction,jx,jy,jmax_segments,ngauss
    integer (kind=int_kind), intent(in)               :: nvertex
    !
    ! for Gaussian quadrature
    !
    real (kind=real_kind), dimension(ngauss), intent(in) :: gauss_weights, abscissae
    real (kind=real_kind), dimension(1:nvertex)        , intent(in)    :: x_in,y_in

    integer (kind=int_kind), intent(in)               :: jx_min, jy_min, jx_max, jy_max
    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in) :: xgno
    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in) :: ygno
    integer (kind=int_kind),            intent(inout) :: jsegment
!    integer (kind=int_kind),dimension(0:2),intent(in)    :: jx_eul_in, jy_eul_in
    real (kind=real_kind)   ,  &
         dimension(jmax_segments,nreconstruction), intent(out) :: weights
    integer (kind=int_kind),  &
         dimension(jmax_segments,2), intent(out) :: weights_eul_index


    !
    ! variables for registering crossings with Eulerian latitudes and longitudes
    !
    integer (kind=int_kind),         intent(inout):: jcross_lat
    !
    ! max. crossings per side is 2*nhe
    !
    real (kind=real_kind), &
         dimension(8*nhe,2), intent(inout):: r_cross_lat
    integer (kind=int_kind), &
         dimension(8*nhe,2), intent(inout):: cross_lat_eul_index
    !
    ! local variables
    !
    real (kind=real_kind) :: dist_lon,dist_lat, tmp_a1, tmp_a2, tmp_x(1), tmp_b2, a1, a2, b2
    real (kind=real_kind) :: dist
    real (kind=real_kind), dimension(2) :: xseg,yseg 
    real (kind=real_kind), dimension(0:3) :: x,y
    real (kind=real_kind)               :: lon_x,lat_y,lon_y,lat_x
    real (kind=real_kind)               :: xeul,yeul,xcross,ycross,slope
    integer (kind=int_kind) ::    jx_eul_tmp,jy_eul_tmp
    integer (kind=int_kind)            :: xsgn1,ysgn1,xsgn2,ysgn2
    integer (kind=int_kind) :: ifrom_left, iter,previous_jy_eul_cross
    logical :: lcontinue, lregister_cross, lsame_cell_x, lsame_cell_y

    integer (kind=int_kind) :: jx_eul, jy_eul, side_count,jdbg
    real (kind=real_kind), dimension(0:nvertex+2)  :: xcell,ycell
    real (kind=real_kind), dimension(0:nvertex+2)  :: xcell_tmp,ycell_tmp
    

5   FORMAT(10e14.6)
    !
    !***********************************************
    !
    ! find jx_eul and jy_eul for (x(1),y(1))
    !
    !***********************************************
    !
    jx_eul = jx; jy_eul = jy    
    xcell(1:nvertex)=x_in; ycell(1:nvertex)=y_in
    DO iter=1,nvertex
      CALL truncate_vertex(xcell(iter),jx_eul,xgno)
      CALL truncate_vertex(ycell(iter),jy_eul,ygno)
    END DO
    xcell(0) = xcell(nvertex); xcell(nvertex+1)=xcell(1); xcell(nvertex+2)=xcell(2);
    ycell(0) = ycell(nvertex); ycell(nvertex+1)=ycell(1); ycell(nvertex+2)=ycell(2);


    IF (MAXVAL(xcell).LE.xgno(jx_min).OR.MINVAL(xcell).GE.xgno(jx_max).OR.&
        MAXVAL(ycell).LE.ygno(jy_min).OR.MINVAL(ycell).GE.ygno(jy_max)) THEN
      !
      ! entire cell off panel
      !
    ELSE             
      jx_eul = MIN(MAX(jx,jx_min),jx_max)
      jy_eul = MIN(MAX(jy,jy_min),jy_max)
      CALL which_eul_cell(xcell(1:3),jx_eul,xgno)
      CALL which_eul_cell(ycell(1:3),jy_eul,ygno)
      
      side_count = 1
      DO WHILE (side_count<nvertex+1)
        jdbg = 0
        iter = 0
        lcontinue = .TRUE.
        x(0:3) = xcell(side_count-1:side_count+2); y(0:3) = ycell(side_count-1:side_count+2); 
        DO while (lcontinue)
          iter = iter+1
          IF (iter>10) THEN
            WRITE(*,*) "search not converging",iter
            STOP
          END IF
          lsame_cell_x = (x(2).GE.xgno(jx_eul).AND.x(2).LE.xgno(jx_eul+1))
          lsame_cell_y = (y(2).GE.ygno(jy_eul).AND.y(2).LE.ygno(jy_eul+1))
          IF (lsame_cell_x.AND.lsame_cell_y) THEN
            !
            !****************************
            !
            ! same cell integral
            !
            !****************************
            !
            xseg(1) = x(1); yseg(1) = y(1); xseg(2) = x(2); yseg(2) = y(2)
            jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul; 
            lcontinue = .FALSE.
            !
            ! prepare for next side if (x(2),y(2)) is on a grid line
            !
            IF (x(2).EQ.xgno(jx_eul+1).AND.x(3)>xgno(jx_eul+1)) THEN
              !
              ! cross longitude jx_eul+1
              !
              jx_eul=jx_eul+1
            ELSE IF (x(2).EQ.xgno(jx_eul  ).AND.x(3)<xgno(jx_eul)) THEN
              !
              ! cross longitude jx_eul
              !
              jx_eul=jx_eul-1
            END IF
            IF (y(2).EQ.ygno(jy_eul+1).AND.y(3)>ygno(jy_eul+1)) THEN
              !
              ! register crossing with latitude: line-segments point Northward
              !
              jcross_lat = jcross_lat + 1
              jy_eul     = jy_eul     + 1
              cross_lat_eul_index(jcross_lat,1) = jx_eul
              cross_lat_eul_index(jcross_lat,2) = jy_eul
              r_cross_lat(jcross_lat,1) = x(2)
              r_cross_lat(jcross_lat,2) = y(2)
            ELSE IF (y(2).EQ.ygno(jy_eul  ).AND.y(3)<ygno(jy_eul)) THEN
              !
              ! register crossing with latitude: line-segments point Southward
              !
              jcross_lat = jcross_lat+1
              cross_lat_eul_index(jcross_lat,1) = jx_eul
              cross_lat_eul_index(jcross_lat,2) = jy_eul
              r_cross_lat(jcross_lat,1) = x(2)
              r_cross_lat(jcross_lat,2) = y(2)
              !
              jy_eul=jy_eul-1
            END IF
            lcontinue=.FALSE.
          ELSE
            !
            !****************************
            !
            ! not same cell integral
            !
            !****************************
            !
            IF (lsame_cell_x) THEN
              ysgn1 = (1+INT(SIGN(1.0D0,y(2)-y(1))))/2 !"1" if y(2)>y(1) else "0"
              ysgn2 = INT(SIGN(1.0D0,y(2)-y(1)))       !"1" if y(2)>y(1) else "-1"
              !
              !*******************************************************************************
              !
              ! there is at least one crossing with latitudes but no crossing with longitudes
              !
              !*******************************************************************************
              !
              yeul   = ygno(jy_eul+ysgn1)
              IF (x(1).EQ.x(2)) THEN
                !
                ! line segment is parallel to longitude (infinite slope)
                !
                xcross = x(1)
              ELSE
                slope  = (y(2)-y(1))/(x(2)-x(1))
                xcross = x_cross_eul_lat(x(1),y(1),yeul,slope)
                !
                ! constrain crossing to be "physically" possible
                !
                xcross = MIN(MAX(xcross,xgno(jx_eul)),xgno(jx_eul+1))
                !
                ! debugging
                !
                IF (xcross.GT.xgno(jx_eul+1).OR.xcross.LT.xgno(jx_eul)) THEN
                  WRITE(*,*) "xcross is out of range",jx,jy
                  WRITE(*,*) "xcross-xgno(jx_eul+1), xcross-xgno(jx_eul))",&
                       xcross-xgno(jx_eul+1), xcross-ygno(jx_eul)
                  STOP
                END IF
              END IF
              xseg(1) = x(1); yseg(1) = y(1); xseg(2) = xcross; yseg(2) = yeul
              jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul; 
              !
              ! prepare for next iteration
              !
              x(0) = x(1); y(0) = y(1); x(1) = xcross; y(1) = yeul; jy_eul = jy_eul+ysgn2
              !
              ! register crossing with latitude
              !
              jcross_lat = jcross_lat+1
              cross_lat_eul_index(jcross_lat,1) = jx_eul
              if (ysgn2>0) then                
                cross_lat_eul_index(jcross_lat,2) = jy_eul
              else
                cross_lat_eul_index(jcross_lat,2) = jy_eul+1
              end if
              r_cross_lat(jcross_lat,1) = xcross
              r_cross_lat(jcross_lat,2) = yeul
            ELSE IF (lsame_cell_y) THEN
              !
              !*******************************************************************************
              !
              ! there is at least one crossing with longitudes but no crossing with latitudes
              !
              !*******************************************************************************
              !
              xsgn1 = (1+INT(SIGN(1.0D0,x(2)-x(1))))/2 !"1" if x(2)>x(1) else "0"
              xsgn2 = INT(SIGN(1.0D0,x(2)-x(1))) !"1" if x(2)>x(1) else "-1"
              xeul   = xgno(jx_eul+xsgn1)
              IF (ABS(x(2)-x(1))<fuzzy_width) THEN
                ! fuzzy crossing
                ycross = 0.5*(y(2)-y(1))
              ELSE
                slope  = (y(2)-y(1))/(x(2)-x(1))
                ycross = y_cross_eul_lon(x(1),y(1),xeul,slope)
              END IF
              !
              ! constrain crossing to be "physically" possible
              !
              ycross = MIN(MAX(ycross,ygno(jy_eul)),ygno(jy_eul+1))
              !
              ! debugging
              !
              IF (ycross.GT.ygno(jy_eul+1).OR.ycross.LT.ygno(jy_eul)) THEN
                WRITE(*,*) "ycross is out of range"
                WRITE(*,*) "jx,jy,jx_eul,jy_eul",jx,jy,jx_eul,jy_eul
                WRITE(*,*) "ycross-ygno(jy_eul+1), ycross-ygno(jy_eul))",&
                     ycross-ygno(jy_eul+1), ycross-ygno(jy_eul)
                STOP
              END IF
              xseg(1) = x(1); yseg(1) = y(1); xseg(2) = xeul; yseg(2) = ycross
              jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul; 
              !
              ! prepare for next iteration
              !
              x(0) = x(1); y(0) = y(1); x(1) = xeul; y(1) = ycross; jx_eul = jx_eul+xsgn2
            ELSE
              !
              !*******************************************************************************
              !
              ! there are crossings with longitude(s) and latitude(s)
              !
              !*******************************************************************************
              ! 
              xsgn1 = (1+INT(SIGN(1.0D0,x(2)-x(1))))/2 !"1" if x(2)>x(1) else "0"
              xsgn2 = (INT(SIGN(1.0D0,x(2)-x(1)))) !"1" if x(2)>x(1) else "0"
              xeul   = xgno(jx_eul+xsgn1) 
              ysgn1 = (1+INT(SIGN(1.0D0,y(2)-y(1))))/2 !"1" if y(2)>y(1) else "0"
              ysgn2 = INT(SIGN(1.0D0,y(2)-y(1)))       !"1" if y(2)>y(1) else "-1"
              yeul   = ygno(jy_eul+ysgn1)
              
              slope  = (y(2)-y(1))/(x(2)-x(1))
              IF (ABS(x(2)-x(1))<fuzzy_width) THEN
                ycross = 0.5*(y(2)-y(1))
              ELSE
                ycross = y_cross_eul_lon(x(1),y(1),xeul,slope)
              END IF
              xcross = x_cross_eul_lat(x(1),y(1),yeul,slope)

              
              IF ((xsgn2>0.AND.xcross.LE.xeul).OR.(xsgn2<0.AND.xcross.GE.xeul)) THEN
                !
                ! cross latitude
                !
                xseg(1) = x(1); yseg(1) = y(1); xseg(2) = xcross; yseg(2) = yeul
                jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul; 
                !
                ! prepare for next iteration
                !
                x(0) = x(1); y(0) = y(1); x(1) = xcross; y(1) = yeul; jy_eul = jy_eul+ysgn2
                !
                ! register crossing with latitude
                !
                jcross_lat = jcross_lat+1
                cross_lat_eul_index(jcross_lat,1) = jx_eul
                if (ysgn2>0) then                
                  cross_lat_eul_index(jcross_lat,2) = jy_eul
                else
                  cross_lat_eul_index(jcross_lat,2) = jy_eul+1
                end if
                r_cross_lat(jcross_lat,1) = xcross
                r_cross_lat(jcross_lat,2) = yeul
              ELSE
                !
                ! cross longitude
                !
                xseg(1) = x(1); yseg(1) = y(1); xseg(2) = xeul; yseg(2) = ycross
                jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul; 
                !
                ! prepare for next iteration
                !
                x(0) = x(1); y(0) = y(1); x(1) = xeul; y(1) = ycross; jx_eul = jx_eul+xsgn2
              END IF
              
            END IF
          END IF
          !
          ! register line-segment (don't register line-segment if outside of panel)
          !
          if (jx_eul_tmp>=jx_min.AND.jy_eul_tmp>=jy_min.AND.&
               jx_eul_tmp<=jx_max-1.AND.jy_eul_tmp<=jy_max-1) then
            jsegment=jsegment+1
            weights_eul_index(jsegment,1) = jx_eul_tmp
            weights_eul_index(jsegment,2) = jy_eul_tmp
            call get_weights_gauss(weights(jsegment,1:nreconstruction),&
                 xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
          ELSE
            !
            ! segment outside of panel
            !
          END IF
          
        END DO
        side_count = side_count+1
      END DO
    END IF
  end subroutine side_integral
 

  real (kind=real_kind) function compute_slope(x,y)
    implicit none
    real (kind=real_kind), dimension(2), intent(in) :: x,y
!    if (fuzzy(ABS(y(2)-y(1)),fuzzy_width)==0) then
!      compute_slope = 0.0
!    else
    if (fuzzy(ABS(x(2)-x(1)),fuzzy_width)>0) THEN
      compute_slope = (y(2)-y(1))/(x(2)-x(1))
    else
      compute_slope = bignum
    end if
!  end if
  end function compute_slope

  real (kind=real_kind) function y_cross_eul_lon(x,y,xeul,slope)
    implicit none
    real (kind=real_kind), intent(in) :: x,y
    real (kind=real_kind)              , intent(in) :: xeul,slope
    !    
    ! line: y=a*x+b
    !
    real (kind=real_kind) :: a,b
    
    b = y-slope*x 
    y_cross_eul_lon = slope*xeul+b
  end function y_cross_eul_lon

  real (kind=real_kind) function x_cross_eul_lat(x,y,yeul,slope)
    implicit none
    real (kind=real_kind), intent(in) :: x,y
    real (kind=real_kind)              , intent(in) :: yeul,slope

    if (fuzzy(ABS(slope),fuzzy_width)>0) THEN
        x_cross_eul_lat = x+(yeul-y)/slope
    ELSE
      x_cross_eul_lat = bignum
    END IF
  end function x_cross_eul_lat

  subroutine get_weights_exact(weights,xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
    use fvm_analytic_mod, only: I_00, I_10, I_01, I_20, I_02, I_11
  
    implicit none
    integer (kind=int_kind), intent(in) :: nreconstruction, ngauss
    real (kind=real_kind), dimension(nreconstruction), intent(out) :: weights
    real (kind=real_kind), dimension(ngauss), intent(in) :: gauss_weights, abscissae
    
    
    real (kind=real_kind), dimension(2     ), intent(in) :: xseg,yseg
    !
    ! compute weights
    !
    real (kind=real_kind) :: tmp,slope,b,integral,dx2,xc
    integer (kind=int_kind) :: i

    if(.not. ALLGAUSS) then
      weights(1) = ((I_00(xseg(2),yseg(2))-I_00(xseg(1),yseg(1))))
      if (ABS(weights(1))>1.0) THEN
        WRITE(*,*) "1 exact weights(jsegment)",weights(1),xseg,yseg
        stop
      end if
      if (nreconstruction>1) then
         weights(2) = ((I_10(xseg(2),yseg(2))-I_10(xseg(1),yseg(1))))
         weights(3) = ((I_01(xseg(2),yseg(2))-I_01(xseg(1),yseg(1))))
      endif
      if (nreconstruction>3) then
         weights(4) = ((I_20(xseg(2),yseg(2))-I_20(xseg(1),yseg(1))))
         weights(5) = ((I_02(xseg(2),yseg(2))-I_02(xseg(1),yseg(1))))
         weights(6) = ((I_11(xseg(2),yseg(2))-I_11(xseg(1),yseg(1))))
      endif
    else
      call get_weights_gauss(weights,xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
    endif
  end subroutine get_weights_exact



  subroutine get_weights_gauss(weights,xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
    use fvm_analytic_mod, only: I_00, I_10, I_01, I_20, I_02, I_11
    
    implicit none
    integer (kind=int_kind), intent(in) :: nreconstruction,ngauss
    real (kind=real_kind), dimension(nreconstruction), intent(out) :: weights
    real (kind=real_kind), dimension(2     ), intent(in) :: xseg,yseg
    real (kind=real_kind) :: slope
    !
    ! compute weights
    !
    !
    ! for Gaussian quadrature
    !
    real (kind=real_kind), dimension(ngauss), intent(in) :: gauss_weights, abscissae

    ! if line-segment parallel to x or y use exact formulaes else use qudrature
    !
    real (kind=real_kind) :: tmp,b,integral,dx2,xc,x,y
    integer (kind=int_kind) :: i

!    if (fuzzy(abs(xseg(1) -xseg(2)),fuzzy_width)==0)then
    if (xseg(1).EQ.xseg(2))then
      weights = 0.0D0
    else
      slope    = (yseg(2)-yseg(1))/(xseg(2)-xseg(1))
      b        = yseg(1)-slope*xseg(1)
      dx2      = 0.5D0*(xseg(2)-xseg(1))
      xc       = 0.5D0*(xseg(1)+xseg(2))
      integral = 0.0D0
      do i=1,ngauss
        x        = xc+abscissae(i)*dx2
        y        = slope*x+b
        integral = integral+gauss_weights(i)*F_00(x,y)
      enddo
      weights(1) = integral*dx2  
      if (nreconstruction>1) then
        integral = 0.0D0
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_10(x,y)
        enddo
        weights(2) = integral*dx2  
        integral = 0.0D0
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_01(x,y)
        enddo
        weights(3) = integral*dx2  
      endif
      if (nreconstruction>3) then
        integral = 0.0D0
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_20(x,y)
        enddo
        weights(4) = integral*dx2  
        integral = 0.0D0
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_02(x,y)
        enddo
        weights(5) = integral*dx2  
        integral = 0.0D0
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_11(x,y)
        enddo
        weights(6) = integral*dx2  
      endif
    end if
  end subroutine get_weights_gauss

  real (kind=real_kind) function F_00(x_in,y_in)
    implicit none
    real (kind=real_kind), intent(in) :: x_in,y_in
    real (kind=real_kind) :: x,y,tmp
    !
    x = x_in
    y = y_in
    F_00 =y/((1.0D0+x*x)*SQRT(1.0D0+x*x+y*y))
  end function F_00

  real (kind=real_kind) function F_10(x_in,y_in)
    implicit none
    real (kind=real_kind), intent(in) :: x_in,y_in
    real (kind=real_kind) :: x,y,tmp

    x = x_in
    y = y_in

    F_10 =x*y/((1.0D0+x*x)*SQRT(1.0D0+x*x+y*y))
  end function F_10

  real (kind=real_kind) function F_01(x_in,y_in)
    implicit none
    real (kind=real_kind), intent(in) :: x_in,y_in
    real (kind=real_kind) :: x,y,tmp

    x = x_in
    y = y_in

    F_01 =-1.0D0/(SQRT(1.0D0+x*x+y*y))
  end function F_01

  real (kind=real_kind) function F_20(x_in,y_in)
    implicit none
    real (kind=real_kind), intent(in) :: x_in,y_in
    real (kind=real_kind) :: x,y,tmp

    x = x_in
    y = y_in

    F_20 =x*x*y/((1.0D0+x*x)*SQRT(1.0D0+x*x+y*y))
  end function F_20

  real (kind=real_kind) function F_02(x_in,y_in)
    implicit none
    real (kind=real_kind), intent(in) :: x_in,y_in
    real (kind=real_kind) :: x,y,alpha, tmp

    x = x_in
    y = y_in

    alpha = ATAN(x)
!     F_02 =-y/SQRT(1.0D0+x*x+y*y)+ASINH(y*COS(alpha))
    tmp=y*COS(alpha)
    F_02 =-y/SQRT(1.0D0+x*x+y*y)+log(tmp+sqrt(tmp*tmp+1))
    
    !
    ! cos(alpha) = 1/sqrt(1+x*x)
    !
  end function F_02

  real (kind=real_kind) function F_11(x_in,y_in)
    implicit none
    real (kind=real_kind), intent(in) :: x_in,y_in
    real (kind=real_kind) :: x,y,tmp

    x = x_in
    y = y_in

    F_11 =-x/(SQRT(1.0D0+x*x+y*y))
  end function F_11

  subroutine which_eul_cell(x,j_eul,gno)
    implicit none
    integer (kind=int_kind)                               , intent(inout) :: j_eul
    real (kind=real_kind), dimension(3)                    , intent(in)    :: x
    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in)    :: gno !phl
!    real (kind=real_kind), intent(in)    :: eps
    
    real (kind=real_kind) :: d1,d2,d3,d1p1
    logical                 :: lcontinue
    integer :: iter
    
    lcontinue = .TRUE.
    iter = 0 

    DO WHILE (lcontinue)
      iter = iter+1 
      IF (x(1).GE.gno(j_eul).AND.x(1).LT.gno(j_eul+1)) THEN
        lcontinue = .FALSE.
        !
        ! special case when x(1) is on top of grid line
        !
        IF (x(1).EQ.gno(j_eul)) THEN
          !
          ! x(1) is on top of gno(J_eul)
          !
          IF (x(2).GT.gno(j_eul)) THEN
            j_eul = j_eul
          ELSE IF (x(2).LT.gno(j_eul)) THEN
            j_eul = j_eul-1
          ELSE
            !
            ! x(2) is on gno(j_eul) grid line; need x(3) to determine Eulerian cell 
            !
            IF (x(3).GT.gno(j_eul)) THEN
              !
              ! x(3) to the right
              !
              j_eul = j_eul
            ELSE IF (x(3).LT.gno(j_eul)) THEN
              !
              ! x(3) to the left
              !
              j_eul = j_eul-1
            ELSE
              WRITE(*,*) "inconsistent cell: x(1)=x(2)=x(3)"
              STOP
            END IF
          END IF
        END IF
      ELSE
        ! 
        ! searching - prepare for next iteration
        !
        IF (x(1).GE.gno(j_eul+1)) THEN
          j_eul = j_eul + 1
        ELSE
          !
          ! x(1).LT.gno(j_eul)
          !
          j_eul = j_eul - 1
        END IF
      END IF
      IF (iter>1000.OR.j_eul<-nhe.OR.j_eul>nc+2+nhe) THEN
        WRITE(*,*) "search is which_eul_cell not converging!", iter, j_eul,nhe,nc+2+nhe
        WRITE(*,*) "gno", gno(nc+2+nhe), gno(-nhe)
        write(*,*) gno
        STOP
      END IF
    END DO
  END subroutine which_eul_cell


  subroutine truncate_vertex(x,j_eul,gno)
    implicit none
    integer (kind=int_kind)                               , intent(inout) :: j_eul
    real (kind=real_kind)                    , intent(inout)    :: x
    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in)    :: gno !phl
!    real (kind=real_kind), intent(in)    :: eps
    
    logical                 :: lcontinue
    integer :: iter, xsgn
    real (kind=real_kind) :: dist,dist_new,tmp
    
    lcontinue = .TRUE.
    iter = 0 
    dist = bignum
    xsgn     = INT(SIGN(1.0D00,x-gno(j_eul)))
    
    DO WHILE (lcontinue)
      if ((j_eul<-nhe) .or. (j_eul>nc+2+nhe)) then
        write(*,*) 'somthing is wrong', j_eul, -nhe,nc+2+nhe, iter
        stop
      endif
      iter     = iter+1 
      tmp      = x-gno(j_eul)
      dist_new = ABS(tmp)
      IF (dist_new>dist) THEN
        lcontinue = .FALSE.
!      ELSE IF (ABS(tmp)<1.0E-11) THEN
      ELSE IF (ABS(tmp)<1.0E-9) THEN
!      ELSE IF (ABS(tmp)<1.0E-4) THEN
        x = gno(j_eul)
        lcontinue = .FALSE.
      ELSE
        j_eul = j_eul+xsgn
        dist = dist_new
      END IF
      IF (iter>100) THEN
        WRITE(*,*) "truncate vertex not converging"
        STOP
      END IF
    END DO
  END subroutine truncate_vertex




!********************************************************************************
!
! Gauss-Legendre quadrature
!
! Tabulated values
!
!********************************************************************************
subroutine gauss_points(n,weights,points)
  implicit none
  integer (kind=int_kind)           , intent(in ) :: n
  real (kind=real_kind), dimension(n), intent(out) :: weights, points
  
  select case (n)
!    CASE(1)
!       abscissae(1) = 0.0D0
!       weights(1)   = 2.0D0
  case(2)
     points(1)    = -sqrt(1.0D0/3.0D0)
     points(2)    =  sqrt(1.0D0/3.0D0)
     weights(1)   =  1.0D0
     weights(2)   =  1.0D0
  case(3)
     points(1)    = -0.774596669241483377035853079956D0
     points(2)    =  0.0D0
     points(3)    =  0.774596669241483377035853079956D0
     weights(1)   =  0.555555555555555555555555555556D0
     weights(2)   =  0.888888888888888888888888888889D0
     weights(3)   =  0.555555555555555555555555555556D0
  case(4)
     points(1)    = -0.861136311594052575223946488893D0
     points(2)    = -0.339981043584856264802665659103D0
     points(3)    =  0.339981043584856264802665659103D0
     points(4)    =  0.861136311594052575223946488893D0
     weights(1)   =  0.347854845137453857373063949222D0
     weights(2)   =  0.652145154862546142626936050778D0 
     weights(3)   =  0.652145154862546142626936050778D0 
     weights(4)   =  0.347854845137453857373063949222D0      
  case(5)
     points(1)    = -(1.0D0/3.0D0)*sqrt(5.0D0+2.0D0*sqrt(10.0D0/7.0D0))
     points(2)    = -(1.0D0/3.0D0)*sqrt(5.0D0-2.0D0*sqrt(10.0D0/7.0D0))
     points(3)    =  0.0D0
     points(4)    =  (1.0D0/3.0D0)*sqrt(5.0D0-2.0D0*sqrt(10.0D0/7.0D0))
     points(5)    =  (1.0D0/3.0D0)*sqrt(5.0D0+2.0D0*sqrt(10.0D0/7.0D0))
     weights(1)   = (322.0D0-13.0D0*sqrt(70.0D0))/900.0D0
     weights(2)   = (322.0D0+13.0D0*sqrt(70.0D0))/900.0D0
     weights(3)   = 128.0D0/225.0D0
     weights(4)   = (322.0D0+13.0D0*sqrt(70.0D0))/900.0D0
     weights(5)   = (322.0D0-13.0D0*sqrt(70.0D0))/900.0D0
  case default
     write(*,*) 'n out of range in glwp of module gll. n=',n
     write(*,*) '0<n<5'
     stop
  end select

end subroutine gauss_points

!------------------------------------------------------------------------------
! FUNCTION SIGNUM_FUZZY
!
! Description:
!   Gives the sign of the given real number, returning zero if x is within 
!     a small amount from zero.
!------------------------------------------------------------------------------
  function signum_fuzzy(x)
    implicit none

    real (kind=real_kind) :: signum_fuzzy
    real (kind=real_kind) :: x

    IF (x > fuzzy_width) THEN
      signum_fuzzy = 1.0D0
    ELSEIF (x < fuzzy_width) THEN
      signum_fuzzy = -1.0D0
    ELSE
      signum_fuzzy = 0.0D0
    ENDIF
  end function

  function fuzzy(x,epsilon)
    implicit none

    integer (kind=int_kind) :: fuzzy
    real (kind=real_kind), intent(in) :: epsilon
    real (kind=real_kind) :: x

    IF (ABS(x)<epsilon) THEN
      fuzzy = 0
    ELSE IF (x >epsilon) THEN
      fuzzy = 1
    ELSE !IF (x < fuzzy_width) THEN
      fuzzy = -1
    ENDIF
  end function

end module fvm_line_integrals_flux_mod
