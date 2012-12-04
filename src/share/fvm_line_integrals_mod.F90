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
module fvm_line_integrals_mod

  use kinds, only               : int_kind, real_kind
  use dimensions_mod, only      : nc, nhe, ngpc
  use parallel_mod, only : abortmp

  implicit none
  private
  real (kind=real_kind),parameter, public   :: bignum = 1.0D20
  real (kind=real_kind),parameter, public   :: tiny   = 1.0D-12
  real (kind=real_kind),parameter           :: fuzzy_width = 10.0*tiny
  real (kind=real_kind)                     :: toleul=1.0D-10
  
  logical :: ldbg=.FALSE.
  
  public :: compute_weights
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
subroutine compute_weights(fvm,nreconstruction,weights_all,weights_eul_index_all, &
                                           weights_lgr_index_all,klev,jall)  
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
  
  type (cartesian2D_t)                        :: dcart(-1:nc+3,-1:nc+3)       ! Cartesian coordinates 
  
  real (kind=real_kind), dimension(0:5)       :: xcell,ycell
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
  
  integer (kind=int_kind)                    :: jallactual, jallactual_eul, ja
  real (kind=real_kind)                      :: da_cslam(1-nhe:nc+nhe,1-nhe:nc+nhe), centroid_cslam(5,1-nhe:nc+nhe,1-nhe:nc+nhe)
  real (kind=real_kind)                      :: area

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
  jall = 1
  ! 
  ! calculate xy Cartesian on the cube of departure points on the corresponding face  
  do jy=-1,nc+3
     do jx=-1,nc+3  
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%faceno,dcart(jx,jy))  
     end do
  end do
  
  do jy=1, nc
    do jx=1, nc            
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
      call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
           fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
           tmp,ngpc,gsweights,gspts,&
           weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                      weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
  end do
  jallactual=jall
  
!!!!  WEIGHTS CORRECTION FOR THE interior element/cells
    do jy=jy_min-1, 0
      do jx=jx_min-1, jx_max  
        if ((fvm%cubeboundary == swest) .and. (jx<1) .and. (jy<1)) then
          cycle
        endif      
        if ((fvm%cubeboundary == seast) .and. (jx>nc) .and. (jy<1)) then
          cycle
        endif     
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)  
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
             fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
        if (jcollect_cell>0) then
          weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
          weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                        weights_eul_index_cell(1:jcollect_cell,:)
          weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
          weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
          jall = jall+jcollect_cell          
        endif
      end do
    end do
    
    do jy=nc+1, jy_max
      do jx=jx_min-1, jx_max      
        if ((fvm%cubeboundary == nwest) .and. (jx<1) .and. (jy>nc)) then
          cycle
        endif      
        if ((fvm%cubeboundary == neast) .and. (jx>nc) .and. (jy>nc)) then
          cycle
        endif        
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)  
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
             fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
        if (jcollect_cell>0) then
          weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
          weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                        weights_eul_index_cell(1:jcollect_cell,:)
          weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
          weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
          jall = jall+jcollect_cell          
        endif
      end do
    end do
    do jx=jx_min-1, 0
      do jy=1, nc             
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)  
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
             fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
        if (jcollect_cell>0) then
          weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
          weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                        weights_eul_index_cell(1:jcollect_cell,:)
          weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
          weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
          jall = jall+jcollect_cell          
        endif
      end do
    end do
    do jx=nc+1, jx_max
      do jy=1, nc             
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)  
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
             fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
        if (jcollect_cell>0) then
          weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
          weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                                        weights_eul_index_cell(1:jcollect_cell,:)
          weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
          weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
          jall = jall+jcollect_cell          
        endif
      end do
    end do
    da_cslam=0.0D0
    centroid_cslam=0.0D0
    do ja=1,jall-1
      jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
      da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
      centroid_cslam(:,jx,jy) = centroid_cslam(:,jx,jy)+weights_all(ja,2:6)
    end do   
  
  jall=jallactual
  jallactual_eul=jall
  
  !WEST SIDE
  if (fvm%cubeboundary == west) then
    ! calculate xy Cartesian on the cube of departure points on the corresponding face
    do jx=-1,2      
      do jy=-1,nc+3
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(west),dcart(jx,jy))                  
      end do
    end do
    jx=1
    do jy=1,nc
!       call getdep_cellboundaries(xcell,ycell,jx,jy,fvm%nbrsface(west),fvm%dsphere) 
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     

      if(swap1) then  !flip orientation
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    
    jallactual=jall
   ! for Eulerian Correction  
    do jx=-1,1      
      do jy=-1,nc+2
       if ((jx==1) .and. (jy>0) .and. (jy<nc+1)) then 
         cycle
       endif 
  !       call getdep_cellboundaries(xcell,ycell,jx,jy,fvm%nbrsface(west),fvm%dsphere) 
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     

        if(swap1) then  !flip orientation
          call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    end do
    
    do ja=jallactual_eul,jall-1
      jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
      da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
      centroid_cslam(:,jx,jy) = centroid_cslam(:,jx,jy)+weights_all(ja,2:6)
    end do
!     
    jall=jallactual
  endif
    
  !EAST SIDE
  if (fvm%cubeboundary == east) then
    ! calculate xy Cartesian on the cube of departure points on the corresponding face 
    do jx=nc,nc+3 
      do jy=-1,nc+3
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(east),dcart(jx,jy))                  
      end do
    end do
    jx=nc
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    
    jallactual=jall
    ! for Eulerian Correction  
    do jx=nc,nc+2      
      do jy=-1,nc+2
       if ((jx==nc) .and. (jy>0) .and. (jy<nc+1)) then 
         cycle
       endif 
         call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
          if(swap1) then !flip orientation
            call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    end do

    do ja=jallactual_eul,jall-1
      jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
      da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
      centroid_cslam(:,jx,jy) = centroid_cslam(:,jx,jy)+weights_all(ja,2:6)
    end do

    jall=jallactual
  endif  
  
  !NORTH SIDE 
  if (fvm%cubeboundary == north) then
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=nc,nc+3
      do jx=-1,nc+3
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(north),dcart(jx,jy))                  
      end do
    end do
    jy=nc
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
   
   jallactual=jall
   ! for Eulerian Correction      
    do jy=nc,nc+2      
      do jx=-1,nc+2
       if ((jy==nc) .and. (jx>0) .and. (jx<nc+1)) then 
         cycle
       endif 
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
        if(swap1) then !flip orientation
          call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    end do

    do ja=jallactual_eul,jall-1
      jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
      da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
      centroid_cslam(:,jx,jy) = centroid_cslam(:,jx,jy)+weights_all(ja,2:6)
    end do
!     
    jall=jallactual    
  end if

  !SOUTH SIDE
  if (fvm%cubeboundary == south) then
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=-1,2
      do jx=-1,nc+3
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(south),dcart(jx,jy))                   
      end do
    end do
    jy=1
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    
    jallactual=jall
    
!for Eulerian correction    
    do jy=-1,1      
      do jx=-1,nc+2
       if ((jy==1) .and. (jx>0) .and. (jx<nc+1)) then 
         cycle
       endif 
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
        if(swap1) then !flip orientation
          call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    end do

    do ja=jallactual_eul,jall-1
      jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
      da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
      centroid_cslam(:,jx,jy) = centroid_cslam(:,jx,jy)+weights_all(ja,2:6)
    end do

    jall=jallactual    
  end if

  !SOUTHWEST Corner
  if (fvm%cubeboundary == swest) then
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=-1,2
      do jx=-1,nc+3
        if ((jy<1) .and. (jx<1)) then   ! in the southwest corner are no values!!!
          cycle
        end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(south),dcart(jx,jy))                  
      end do
    end do
    jy=1
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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

    jallactual=jall
  
    do jy=-1,1      
      do jx=-1,nc+2
       if (((jy<1) .and. (jx<1)) .or. ((jy==1) .and. (jx>0) .and. (jx<nc+1))) then 
         cycle
       endif 
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
        if(swap1) then !flip orientation
          call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    end do

    do ja=jallactual_eul,jall-1
      jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
      da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
      centroid_cslam(:,jx,jy) = centroid_cslam(:,jx,jy)+weights_all(ja,2:6)
    end do
    
    
    jall=jallactual    
    jallactual_eul=jall  


    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=-1,2
      do jy=-1,nc+3
        if ((jy<1) .and. (jx<1)) then
          cycle
        end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(west),dcart(jx,jy))                  
      end do
    end do
    jx=1
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap2) then !flip orientation
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    
    jallactual=jall
    
    do jx=-1,1      
      do jy=-1,nc+2
       if (((jx<1) .and. (jy<1)) .or. ((jx==1) .and. (jy>0) .and. (jy<nc+1))) then 
         cycle
       endif
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
        if(swap2) then !flip orientation
          call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    end do

    do ja=jallactual_eul,jall-1
      jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
      da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
      centroid_cslam(:,jx,jy) = centroid_cslam(:,jx,jy)+weights_all(ja,2:6)
    end do
!     
    jall=jallactual      
    
  endif

  ! SOUTHEAST Corner
  if (fvm%cubeboundary == seast) then
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=-1,2
      do jx=-1,nc+3
        if ((jy<1) .and. (jx>nc+1)) then   ! in the southwest corner are no values!!!
          cycle
        end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(south),dcart(jx,jy))  
      end do
    end do
    jy=1
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    jallactual=jall
! Eulerian weight correction    
    do jy=-1,1
      do jx=-1,nc+2
        if (((jy<1) .and. (jx>nc)) .or. ((jy==1) .and. (jx>0) .and. (jx<nc+1))) then 
          cycle
        endif
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
        if(swap1) then !flip orientation
          call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    end do

    do ja=jallactual_eul,jall-1
      jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
      da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
      centroid_cslam(:,jx,jy) = centroid_cslam(:,jx,jy)+weights_all(ja,2:6)
    end do

    jall=jallactual    
    jallactual_eul=jall

    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=nc,nc+3
      do jy=-1,nc+3
        if ((jy<1) .and. (jx>nc+1)) then
          cycle
        end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(east),dcart(jx,jy))                   
      end do
    end do
    jx=nc
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
      if(swap2) then !flip orientation
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    
    jallactual=jall
    ! Eulerian weight correction        
    do jx=nc,nc+2
      do jy=-1,nc+2
        if (((jx>nc) .and. (jy<1)) .or. ((jx==nc) .and. (jy>0) .and. (jy<nc+1))) then 
          cycle
        endif        
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
        if(swap2) then !flip orientation
          call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    end do 
    
    do ja=jallactual_eul,jall-1
      jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
      da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
      centroid_cslam(:,jx,jy) = centroid_cslam(:,jx,jy)+weights_all(ja,2:6)
    end do
  !     
    jall=jallactual
    
  endif
  
  !NORTHEAST Corner
  if (fvm%cubeboundary == neast) then
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=nc,nc+3
      do jx=-1,nc+3
        if ((jy>nc+1) .and. (jx>nc+1)) then   ! in the southwest corner are no values!!!
          cycle
        end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(north),dcart(jx,jy))                
      end do
    end do
    jy=nc
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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


    jallactual=jall

    do jy=nc,nc+2      
      do jx=-1,nc+2
        if (((jy>nc) .and. (jx>nc)) .or. ((jy==nc) .and. (jx>0) .and. (jx<nc+1))) then 
          cycle
        endif
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
        if(swap1) then !flip orientation
          call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    end do

    do ja=jallactual_eul,jall-1
      jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
      da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
      centroid_cslam(:,jx,jy) = centroid_cslam(:,jx,jy)+weights_all(ja,2:6)
    end do


    jall=jallactual    
    jallactual_eul=jall  


    
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=nc,nc+3
      do jy=-1,nc+3
        if ((jy>nc+1) .and. (jx>nc+1)) then
          cycle
        end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(east),dcart(jx,jy))                  
      end do
    end do
    jx=nc
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap2) then !flip orientation
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    
    jallactual=jall
    
    do jx=nc, nc+2
      do jy=-1,nc+2
        if (((jx>nc) .and. (jy>nc)) .or. ((jx==nc) .and. (jy>0) .and. (jy<nc+1))) then 
          cycle
        endif
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
        if(swap2) then !flip orientation
          call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    end do
    
    do ja=jallactual_eul,jall-1
      jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
      da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
      centroid_cslam(:,jx,jy) = centroid_cslam(:,jx,jy)+weights_all(ja,2:6)
    end do
!     
    jall=jallactual    
    
  endif

  !NORTH WEST CORNER 
  if (fvm%cubeboundary == nwest) then
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=nc, nc+3
      do jx=-1,nc+3
        if ((jy>nc+1) .and. (jx<1)) then   ! in the southwest corner are no values!!!
          cycle
        end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(north),dcart(jx,jy))                   
      end do
    end do
    jy=nc
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
  
    jallactual=jall

    do jy=nc,nc+2
      do jx=-1,nc+2
        if (((jy>nc) .and. (jx<1)) .or. ((jy==nc) .and. (jx>0) .and. (jx<nc+1))) then 
          cycle
        endif
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
        if(swap1) then !flip orientation
          call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    end do
    do ja=jallactual_eul,jall-1
      jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
      da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
      centroid_cslam(:,jx,jy) = centroid_cslam(:,jx,jy)+weights_all(ja,2:6)
    end do
    jall=jallactual    
    jallactual_eul=jall  

    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=-1,2
      do jy=-1,nc+3
        if ((jy>nc+1) .and. (jx<1)) then
          cycle
        end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(west),dcart(jx,jy))                  
      end do
    end do
    jx=1
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap2) then !flip orientation
        call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    
    jallactual=jall
    do jx=-1,1
      do jy=-1,nc+2
        if (((jx<1) .and. (jy>nc)) .or. ((jx==1) .and. (jy>0) .and. (jy<nc+1))) then 
          cycle
        endif
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
        if(swap2) then !flip orientation
          call compute_weights_cell(xcell,ycell,jx,jy,nreconstruction,&
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
    end do
    
    do ja=jallactual_eul,jall-1
      jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
      da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
      centroid_cslam(:,jx,jy) = centroid_cslam(:,jx,jy)+weights_all(ja,2:6)
    end do
!     
    jall=jallactual       
  endif
  
  
  jall=jall-1
!   toleul=1.0D-11
!   do jx=0,nc+1
!     do jy=0,nc+1
! !         write(*,*) 'area', jx, jy, da_cslam(jx,jy)-fvm%area_sphere(jx,jy)
!       if (abs(da_cslam(jx,jy)-fvm%area_sphere(jx,jy))>toleul) then
!         write(*,*) 'area diff', jx, jy, da_cslam(jx,jy)-fvm%area_sphere(jx,jy)
!         write(*,*) 'cubeboundary', fvm%cubeboundary, 'faceno', fvm%faceno
! !           stop
!       endif
!       if (abs(centroid_cslam(1,jx,jy)-fvm%area_sphere(jx,jy)*fvm%spherecentroid(1,jx,jy))>toleul) then
!         write(*,*) 'centroid 1 diff', jx, jy, centroid_cslam(1,jx,jy)-fvm%area_sphere(jx,jy)*fvm%spherecentroid(1,jx,jy)
!         write(*,*) 'cubeboundary', fvm%cubeboundary, 'faceno', fvm%faceno
! !           stop
!       endif
!       if (abs(centroid_cslam(2,jx,jy)-fvm%area_sphere(jx,jy)*fvm%spherecentroid(2,jx,jy))>toleul) then
!         write(*,*) 'centroid 2 diff', jx, jy, centroid_cslam(2,jx,jy)-fvm%area_sphere(jx,jy)*fvm%spherecentroid(2,jx,jy)
!         write(*,*) 'cubeboundary', fvm%cubeboundary, 'faceno', fvm%faceno
! !           stop
!       endif
!       if (abs(centroid_cslam(3,jx,jy)-fvm%area_sphere(jx,jy)*fvm%spherecentroid(3,jx,jy))>toleul) then
!         write(*,*) 'centroid 3 diff', jx, jy, centroid_cslam(3,jx,jy)-fvm%area_sphere(jx,jy)*fvm%spherecentroid(3,jx,jy)
!         write(*,*) 'cubeboundary', fvm%cubeboundary, 'faceno', fvm%faceno
! !           stop
!       endif
!       if (abs(centroid_cslam(4,jx,jy)-fvm%area_sphere(jx,jy)*fvm%spherecentroid(4,jx,jy))>toleul) then
!         write(*,*) 'centroid 4 diff', jx, jy, centroid_cslam(4,jx,jy)-fvm%area_sphere(jx,jy)*fvm%spherecentroid(4,jx,jy)
!         write(*,*) 'cubeboundary', fvm%cubeboundary, 'faceno', fvm%faceno
! !           stop
!       endif
!       if (abs(centroid_cslam(5,jx,jy)-fvm%area_sphere(jx,jy)*fvm%spherecentroid(5,jx,jy))>toleul) then
!         write(*,*) 'centroid 5 diff', jx, jy, centroid_cslam(5,jx,jy)-fvm%area_sphere(jx,jy)*fvm%spherecentroid(5,jx,jy)
!         write(*,*) 'cubeboundary', fvm%cubeboundary, 'faceno', fvm%faceno
! !           stop
!       endif        
! !         write(*,*) 'centroid', jx, jy, centroid_cslam(:,jx,jy)-fvm%area_sphere(jx,jy)*fvm%spherecentroid(:,jx,jy)
!     end do
!   end do
 
! here is the correction, uncomment it if you want to run the scheme without weight correction!  
!   do ja=1,jall
!    jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
!    area=fvm%area_sphere(jx,jy)
! !    
!    weights_all(ja,1) = weights_all(ja,1)*abs(area/da_cslam(jx,jy))
!    weights_all(ja,2) = weights_all(ja,2)*abs(area*fvm%spherecentroid(1,jx,jy)/centroid_cslam(1,jx,jy))
!    weights_all(ja,3) = weights_all(ja,3)*abs(area*fvm%spherecentroid(2,jx,jy)/centroid_cslam(2,jx,jy))
!    weights_all(ja,4) = weights_all(ja,4)*abs(area*fvm%spherecentroid(3,jx,jy)/centroid_cslam(3,jx,jy))  
!    weights_all(ja,5) = weights_all(ja,5)*abs(area*fvm%spherecentroid(4,jx,jy)/centroid_cslam(4,jx,jy))
!    weights_all(ja,6) = weights_all(ja,6)*abs(area*fvm%spherecentroid(5,jx,jy)/centroid_cslam(5,jx,jy))
!   end do
  
 
end subroutine compute_weights  


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
  type (cartesian2D_t), intent(in)                         :: dcart(-1:nc+3,-1:nc+3)

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
!phl
    !
    ! ipanel is just for debugging
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

!    IF (ltest_search) THEN
!      IF (ipanel==1.AND.((jx==22.AND.jy==2).OR.(jx==22.AND.jy==3))) THEN
!        ldbg = .true. 
!      else
!        ldbg = .false.
!      end if
!    END IF
    
    

    if (ldbg) WRITE(*,*) "*********************************************"
    if (ldbg) WRITE(*,*) "looking at cell (Eulerian index):",jx,jy


    jsegment = 0
    weights  = 0.0D0
    jcross_lat = 0
          

    IF (ldbg) WRITE(*,*) xcell(1),ycell(1)
    IF (ldbg) WRITE(*,*) xcell(2),ycell(2)
    IF (ldbg) WRITE(*,*) xcell(3),ycell(3)
    IF (ldbg) WRITE(*,*) xcell(4),ycell(4)
    IF (ldbg) WRITE(*,*) xcell(1),ycell(1)

    IF (ldbg) WRITE(*,*) "vertex 1  displacement",xcell(1)-xgno(jx),ycell(1)-ygno(jy)
    IF (ldbg) WRITE(*,*) "vertex 2  displacement",xcell(2)-xgno(jx),ycell(2)-ygno(jy+1)
    IF (ldbg) WRITE(*,*) "vertex 3  displacement",xcell(3)-xgno(jx+1),ycell(3)-ygno(jy+1)
    IF (ldbg) WRITE(*,*) "vertex 4  displacement",xcell(4)-xgno(jx+1),ycell(4)-ygno(jy)
    IF (ldbg) WRITE(*,*) "---------------------"




!        do i=1,4
!          IF (ldbg) WRITE(*,*) "----------------------------------------------------"
!          IF (ldbg) WRITE(*,*) "starting side ",i
!          x(0) = xcell(i-1); x(1) = xcell(i); x(2) = xcell(i+1)
!          y(0) = ycell(i-1); y(1) = ycell(i); y(2) = ycell(i+1)
!    WRITE(*,*) jx,jy
    call side_integral(xcell,ycell,4,jsegment,jmax_segments,&
         weights,weights_eul_index,nreconstruction,jx,jy,xgno,ygno,jx_min, jx_max, jy_min, jy_max,&
         ngauss,gauss_weights,abscissae,&
         jcross_lat,r_cross_lat,cross_lat_eul_index)

!        enddo
        !
       !**********************
       ! 
       ! Do inner integrals
       !
       !**********************
       !

     if (ldbg) WRITE(*,*) "jcross_lat",jcross_lat
       
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
       !         jcollect = jsegment
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
       !          write(*,*) "tmp=",tmp
       !          if (ldbg) then
       IF (abs(tmp)>0.04) THEN
         WRITE(*,*) "sum of weights too large",tmp
         !              ldbg=.TRUE.
         stop
       END IF
       IF (tmp<-1.0E-9) THEN
         WRITE(*,*) "sum of weights is negative - negative area?",tmp,jx,jy
         !              ldbg=.TRUE.
         stop
       END IF
       !          end if
     else
       jcollect = 0
     end if
!   else
!     jcollect = 0
!     IF (ldbg) WRITE(*,*) "line segment outside of bounds/panel"
!   end if
 
   
   IF (ldbg) WRITE(*,*) "done cell"
   IF (ldbg) WRITE(*,*) "*********************************************"
   

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
    integer (kind=int_kind)                                 , intent(in)    :: nreconstruction
    real (kind=real_kind)   , dimension(jmax_segments,nreconstruction), intent(inout) :: weights
    integer (kind=int_kind), dimension(jmax_segments,2     ), intent(inout) :: weights_eul_index
    integer (kind=int_kind),                                  INTENT(OUT  ) :: jcollect
    integer (kind=int_kind),                                  INTENT(IN   ) :: jsegment,jmax_segments
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
!             weights_eul_index_out(jcollect,:) = weights_eul_index(jcollect,:)
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
       nreconstruction,ngauss,gauss_weights,abscissae)!phl add jx_min etc.
    
    implicit none
    !
    ! for Gaussian quadrature
    !
    real (kind=real_kind), dimension(ngauss), intent(in) :: gauss_weights, abscissae
    !
    ! variables for registering crossings with Eulerian latitudes and longitudes
    !
    integer (kind=int_kind),         intent(in):: jcross_lat, jmax_segments,nreconstruction,ngauss
    integer (kind=int_kind),         intent(inout):: jsegment
    !
    ! max. crossings per side is 2*nhe
    !
    real (kind=real_kind), &
         dimension(8*nhe,2), intent(in):: r_cross_lat
    integer (kind=int_kind), &
         dimension(8*nhe,2), intent(in):: cross_lat_eul_index
!phl    real (kind=real_kind), dimension(0:nc+2), intent(in)    :: xgno
    integer (kind=int_kind), intent(in)            ::jx_min, jx_max, jy_min, jy_max
    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in) :: xgno
    real (kind=real_kind)   ,  &
         dimension(jmax_segments,nreconstruction), intent(inout) :: weights
    integer (kind=int_kind),  &
         dimension(jmax_segments,2), intent(inout) :: weights_eul_index
    real (kind=real_kind)   , dimension(nreconstruction) :: weights_tmp

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
!                   if (abs(rend_tmp(1)-rend(1))<tiny) then
!                      EXIT !are we done already?
!                   else
                      rstart = rend_tmp
!                   endif
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
    !
    ! for Gaussian quadrature
    !
    real (kind=real_kind), dimension(ngauss), intent(in) :: gauss_weights, abscissae
    real (kind=real_kind), dimension(1:nvertex)        , intent(in)    :: x_in,y_in

    integer (kind=int_kind), intent(in)               :: jx_min, jy_min, jx_max, jy_max
    integer (kind=int_kind), intent(in)               :: nvertex
    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in) :: xgno
    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in) :: ygno
    integer (kind=int_kind),            intent(inout) :: jsegment
!    integer (kind=int_kind),dimension(0:2),intent(in)    :: jx_eul_in, jy_eul_in
    integer (kind=int_kind),            intent(in)    :: nreconstruction,jx,jy,jmax_segments,ngauss
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

    IF (ldbg) THEN
      WRITE(*,*) "from side_integral: cell vertices"
      DO iter=1,nvertex
        WRITE(*,*) "x(iter),y(iter)",iter, xcell(iter),ycell(iter)
      END DO
    END IF




    IF (MAXVAL(xcell).LE.xgno(jx_min).OR.MINVAL(xcell).GE.xgno(jx_max).OR.&
        MAXVAL(ycell).LE.ygno(jy_min).OR.MINVAL(ycell).GE.ygno(jy_max)) THEN
      IF (ldbg)  WRITE(*,*) "entire cell off panel"
    ELSE             
      jx_eul = MIN(MAX(jx,jx_min),jx_max)
      jy_eul = MIN(MAX(jy,jy_min),jy_max)
      CALL which_eul_cell(xcell(1:3),jx_eul,xgno)
      CALL which_eul_cell(ycell(1:3),jy_eul,ygno)
      IF (ldbg) WRITE(*,*) "x(1),y(1) in cell",jx_eul,jy_eul
      
      side_count = 1
      DO WHILE (side_count<nvertex+1)
        jdbg = 0
        iter = 0
        lcontinue = .TRUE.
        x(0:3) = xcell(side_count-1:side_count+2); y(0:3) = ycell(side_count-1:side_count+2); 
        IF (ldbg) WRITE(*,*) "+++++++++++++++++++++++++++++++++++++++"
        IF (ldbg) WRITE(*,*) "side",side_count
        DO while (lcontinue)
          IF (ldbg) WRITE(*,*) "iter",iter
          IF (ldbg) WRITE(*,*) "x,y(1)",x(1),y(1)
          IF (ldbg) WRITE(*,*) "x,y(2)",x(2),y(2)
          IF (ldbg) WRITE(*,*) "jx_eul,jy_eul",jx_eul,jy_eul
          IF (ldbg) WRITE(*,*) "xgno",xgno(jx_eul),xgno(jx_eul+1)
          IF (ldbg) WRITE(*,*) "ygno",ygno(jy_eul),ygno(jy_eul+1)

          iter = iter+1
          IF (iter>10) THEN
            WRITE(*,*) "search not converging",iter
            STOP
          END IF
          lsame_cell_x = (x(2).GE.xgno(jx_eul).AND.x(2).LE.xgno(jx_eul+1))
          lsame_cell_y = (y(2).GE.ygno(jy_eul).AND.y(2).LE.ygno(jy_eul+1))
          IF (ldbg) WRITE(*,*) "lsame_cell_x,lsame_cell_y=",lsame_cell_x,lsame_cell_y
          IF (lsame_cell_x.AND.lsame_cell_y) THEN
            !
            !****************************
            !
            ! same cell integral
            !
            !****************************
            !
            IF (ldbg) WRITE(*,*) "same cell integral",jx_eul,jy_eul
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
              IF (ldbg) WRITE(*,*) "cross longitude",jx_eul+1
              jx_eul=jx_eul+1
            ELSE IF (x(2).EQ.xgno(jx_eul  ).AND.x(3)<xgno(jx_eul)) THEN
              !
              ! cross longitude jx_eul
              !
              IF (ldbg) WRITE(*,*) "cross longitude",jx_eul
              jx_eul=jx_eul-1
            END IF
            IF (y(2).EQ.ygno(jy_eul+1).AND.y(3)>ygno(jy_eul+1)) THEN
              !
              ! register crossing with latitude: line-segments point Northward
              !
              jcross_lat = jcross_lat + 1
              jy_eul     = jy_eul     + 1
              IF (ldbg) WRITE(*,*) "cross latitude",jy_eul
              cross_lat_eul_index(jcross_lat,1) = jx_eul
              cross_lat_eul_index(jcross_lat,2) = jy_eul
              r_cross_lat(jcross_lat,1) = x(2)
              r_cross_lat(jcross_lat,2) = y(2)
            ELSE IF (y(2).EQ.ygno(jy_eul  ).AND.y(3)<ygno(jy_eul)) THEN
              !
              ! register crossing with latitude: line-segments point Southward
              !
              IF (ldbg) WRITE(*,*) "cross latitude",jy_eul
              jcross_lat = jcross_lat+1
              cross_lat_eul_index(jcross_lat,1) = jx_eul
              cross_lat_eul_index(jcross_lat,2) = jy_eul
              r_cross_lat(jcross_lat,1) = x(2)
              r_cross_lat(jcross_lat,2) = y(2)
              
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
              IF (ldbg) WRITE(*,*) "same cell x"
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
                IF (ldbg) WRITE(*,*) "line segment parallel to longitude"
                xcross = x(1)
              ELSE
                slope  = (y(2)-y(1))/(x(2)-x(1))
                xcross = x_cross_eul_lat(x(1),y(1),yeul,slope)
                !
                ! constrain crossing to be "physically" possible
                !
                xcross = MIN(MAX(xcross,xgno(jx_eul)),xgno(jx_eul+1))


                IF (ldbg) WRITE(*,*) "cross latitude"
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
              IF (ldbg) WRITE(*,*) "same cell y"
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
              IF (ldbg) WRITE(*,*) " crossing longitude",jx_eul+xsgn1
              IF (ABS(x(2)-x(1))<fuzzy_width) THEN
                ycross = 0.5*(y(2)-y(1))
                IF (ldbg) WRITE(*,*) "fuzzy crossing"
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
              IF (ldbg) WRITE(*,*) "not same cell x; not same cell y"
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
                IF (ldbg) WRITE(*,*) "crossing latitude",jy_eul+ysgn1
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
                IF (ldbg) WRITE(*,*) "crossing longitude",jx_eul+xsgn1
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

         !    jdbg=jdbg+1
         ! 
         !             wdbg(jx,jy,side_count,jdbg,1:nreconstruction) = weights(jsegment,:)
         ! 
         !             if (xseg(1).EQ.xseg(2))then
         !               slope = bignum
         !             else if (abs(yseg(1) -yseg(2))<fuzzy_width) then
         !               slope = 0.0
         !             else
         !               slope    = (yseg(2)-yseg(1))/(xseg(2)-xseg(1))
         !             end if
         ! 
         !             wdbg          (jx,jy,side_count,jdbg,nreconstruction+1) = slope
         !             wdbg          (jx,jy,side_count,jdbg,nreconstruction+2) = xseg(2)-xseg(1)
         !             wdbg          (jx,jy,side_count,jdbg,nreconstruction+3) = yseg(2)-yseg(1)
         ! 
         ! 
         !             wdbg_eul_index(jx,jy,side_count,jdbg,1) = jx_eul_tmp
         !             wdbg_eul_index(jx,jy,side_count,jdbg,2) = jy_eul_tmp
         !             wdbg2(jx,jy,side_count) = jdbg
          ELSE
            IF (ldbg) WRITE(*,*) "segment outside of panel"
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

    ! line: y=a*x+b
    real (kind=real_kind) :: a,b

!    a = (y(2)-y(1))/(x(2)-x(1))
!    WRITE(*,*) "haeldning",a,(x(2)-x(1)),(y(2)-y(1))
    
    b = y-slope*x 
!    WRITE(*,*) "b from y_cross_lon",b
    !
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
    real (kind=real_kind), dimension(nreconstruction) :: weightsgauss
    real (kind=real_kind), dimension(ngauss), intent(in) :: gauss_weights, abscissae
    
    
    real (kind=real_kind), dimension(2     ), intent(in) :: xseg,yseg
    !
    ! compute weights
    !
    real (kind=real_kind) :: tmp,slope,b,integral,dx2,xc
    integer (kind=int_kind) :: i
!    weights(:) = -half*(xseg(1)*yseg(2)-xseg(2)*yseg(1)) !dummy for testing

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

    call get_weights_gauss(weightsgauss,xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
!     do i=1, 6
!       if (abs(weightsgauss(i)-weights(i))>1.0D-12) then
!          write(*,*)'weightsdiff', i, abs(weightsgauss(i)-weights(i)), weightsgauss(i), weights(i) !/abs(weightsgauss(i))
!               stop
!       endif  
!     enddo
    weights=weightsgauss

    IF (ldbg) THEN
      WRITE(*,*) "from get weights exact"
!      WRITE(*,*) "weights in get_weights gauss",weights
!      WRITE(*,*) "xseg(1),yseg(1)",xseg(1),yseg(1)
!      WRITE(*,*) "xseg(2),yseg(2)",xseg(2),yseg(2)
!      WRITE(*,*) "slope,b",slope,b
!      WRITE(*,*) "dx2",dx2
!      write(*,*) I_00(xseg(2),yseg(2))-I_00(xseg(1),yseg(1))
      DO i=1,nreconstruction
        write(*,*) "i,weights",i,weights(i)
      end do
!      stop
    END IF
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
!      else if (abs(yseg(1) -yseg(2))<fuzzy_width) then
!         !
!         ! line segment parallel to latitude - compute weights exactly
!         !
!     !    if (ldbg) write(*,*) "line segment parallel to latitude - compute weights exactly"
!        weights(1) = ((I_00(xseg(2),yseg(2))-I_00(xseg(1),yseg(1))))
!        if (nreconstruction>1) then
!          weights(2) = ((I_10(xseg(2),yseg(2))-I_10(xseg(1),yseg(1))))
!          weights(3) = ((I_01(xseg(2),yseg(2))-I_01(xseg(1),yseg(1))))
!        endif
!        if (nreconstruction>3) then
!          weights(4) = ((I_20(xseg(2),yseg(2))-I_20(xseg(1),yseg(1))))
!          weights(5) = ((I_02(xseg(2),yseg(2))-I_02(xseg(1),yseg(1))))
!          weights(6) = ((I_11(xseg(2),yseg(2))-I_11(xseg(1),yseg(1))))
!        endif
    else
      
      
      slope    = (yseg(2)-yseg(1))/(xseg(2)-xseg(1))
      b        = yseg(1)-slope*xseg(1)
      dx2      = 0.5D0*(xseg(2)-xseg(1))
      if (ldbg) WRITE(*,*) "dx2 and slope in gauss weight",dx2,slope
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
!    IF (MAXVAL(ABS(weights(1:nreconstruction)))>0.1) THEN
    IF (ldbg) THEN
      WRITE(*,*) "from get weights gauss"
!      WRITE(*,*) "weights in get_weights gauss",weights
!      WRITE(*,*) "xseg(1),yseg(1)",xseg(1),yseg(1)
!      WRITE(*,*) "xseg(2),yseg(2)",xseg(2),yseg(2)
!      WRITE(*,*) "slope,b",slope,b
!      WRITE(*,*) "dx2",dx2
!      write(*,*) I_00(xseg(2),yseg(2))-I_00(xseg(1),yseg(1))
      DO i=1,nreconstruction
        write(*,*) "i,weights",i,weights(i)
      end do
!      stop
    END IF
!    weights = 0.0
!    weights(:) = -half*(xseg(1)*yseg(2)-xseg(2)*yseg(1)) !dummy for testing
  end subroutine get_weights_gauss

  real (kind=real_kind) function F_00(x_in,y_in)
    implicit none
    real (kind=real_kind), intent(in) :: x_in,y_in
    real (kind=real_kind) :: x,y,tmp

    x = x_in
    y = y_in
!    x = x_in
!    y = y_in

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
    IF (ldbg) WRITE(*,*) "from which_eul_cell",x(1),x(2),x(3)
    DO WHILE (lcontinue)
      iter = iter+1 
      IF (x(1).GE.gno(j_eul).AND.x(1).LT.gno(j_eul+1)) THEN
        lcontinue = .FALSE.
        !
        ! special case when x(1) is on top of grid line
        !
        IF (x(1).EQ.gno(j_eul)) THEN
          IF (ldbg) WRITE(*,*) "x(1) is on top of gno(J_eul)"
          IF (x(2).GT.gno(j_eul)) THEN
            j_eul = j_eul
          ELSE IF (x(2).LT.gno(j_eul)) THEN
            j_eul = j_eul-1
          ELSE
            IF (ldbg) WRITE(*,*) "x(2) is on top of gno(J_eul)"
            !
            ! x(2) is on gno(j_eul) grid line; need x(3) to determine Eulerian cell 
            !
            IF (x(3).GT.gno(j_eul)) THEN
              IF (ldbg) WRITE(*,*) "x(3) to the right"
              j_eul = j_eul
            ELSE IF (x(3).LT.gno(j_eul)) THEN
              IF (ldbg) WRITE(*,*) "x(3) to the left x(3)-x(2),x(3),x(2)",x(3)-x(2),x(3),x(2)
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
    integer :: iter
    real (kind=real_kind) :: xsgn,dist,dist_new,tmp
    
    lcontinue = .TRUE.
    iter = 0 
    dist = bignum
    xsgn     = INT(SIGN(1.0D00,x-gno(j_eul)))
    DO WHILE (lcontinue)
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
  real (kind=real_kind), dimension(n), intent(out) :: weights, points
  integer (kind=int_kind)           , intent(in ) :: n
  
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


end module fvm_line_integrals_mod
