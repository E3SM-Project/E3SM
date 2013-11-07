#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!MODULE FVM_LINE_INTEGRALS_MOD_FLUX-------------------------------------------------!
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
  logical                                   :: lexact_horizontal_line_integrals=.TRUE.
  public :: compute_weights_fluxform
contains
  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE COMPUTE_WEIGHTS_FLUXFORM------------------------------------------------!
  !                                                                                   !
  ! THIS CODE IS BASED ON COMPUTE_WEIGHTS CODED BY CHRISTOPH ERATH AND MODIFIED FOR   !
  ! FLUX-FORM CSLAM BY PETER HJORT LAURITZEN                                          !
  !                                                                                   !
  ! DESCRIPTION:                                                                      !
  !-----------------------------------------------------------------------------------!
  subroutine compute_weights_fluxform(fvm,nreconstruction,weights_all,weights_eul_index_all, &
       weights_lgr_index_all,klev,jall)  
    use fvm_control_volume_mod, only:  fvm_struct                                         
    use coordinate_systems_mod,  only :  cartesian2D_t, spherical_polar_t, &
         cart2cubedspherexy, spherical_to_cart
    use physical_constants, only : DD_PI
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    
    use fvm_line_integrals_mod, only: gauss_points, getdep_cellboundariesxyvec
    
    implicit none
    type (fvm_struct), intent(inout)                                :: fvm
    integer (kind=int_kind), intent(in)                            :: nreconstruction
    ! arrays for collecting cell data
    real (kind=real_kind),dimension(10*(nc+2*nhe)*(nc+2*nhe),nreconstruction,2), &
         intent(out) :: weights_all
    integer (kind=int_kind), dimension(10*(nc+2*nhe)*(nc+2*nhe),2,2), &
         intent(out) :: weights_eul_index_all
    integer (kind=int_kind), dimension(10*(nc+2*nhe)*(nc+2*nhe),2,2), &
         intent(out) :: weights_lgr_index_all
    integer (kind=int_kind), intent(in)                       :: klev
    integer (kind=int_kind), dimension(2), intent(out)        :: jall
    
    ! local workspace
    ! max number of line segments is:
    ! (number of longitudes)*(max average number of crossings per line segment = 3)*ncube*2
    !
    
    integer (kind=int_kind)                     :: jx,jy
    integer                                     :: jx_min, jx_max, jy_min, jy_max
    integer                                     :: jx_min1, jx_max1, jy_min1, jy_max1
    integer                                     :: jx_min2, jx_max2, jy_min2, jy_max2
    logical                                     :: swap1, swap2
    
    integer (kind=int_kind)                     :: i, j, jtmp
    
    type (cartesian2D_t)                        :: dcart(-1:nc+3,-1:nc+3)       ! Cartesian coordinates 
    
    real (kind=real_kind), dimension(0:5)       :: xcell,ycell
    integer (kind=int_kind)                     :: inttmp
    real (kind=real_kind)                       :: tmp
    logical                                     :: swap
    ! for Gaussian quadrature
    real (kind=real_kind), dimension(ngpc)      :: gsweights, gspts
    ! weight-variables for individual cells
    integer (kind=int_kind) :: jmax_segments_cell
    real (kind=real_kind)   , dimension(nhe*50,nreconstruction,2)   :: weights_flux_cell
    integer (kind=int_kind),  dimension(nhe*50,2,2)                 :: weights_eul_index_cell
    integer (kind=int_kind), dimension(2)                           :: jcollect_cell
    !
    !
    !    x--x--x--x--x
    !    |  |  |  |  |
    !    x--x--x--x--x
    !    |  |  |  |  |
    !    x--x--x--x--x
    !    |  |  |  |  |
    !    x--x--x--x--x
    !    |  |  |  |  |
    !    x--x--x--x--x
    !
    !
    !            =========
    !            |       |
    !            |       |
    !            |       |
    !            |       |
    !            |       |
    !    =================================
    !    |       |       |       |       |
    !    |       |       |       |       |
    !    |       |       |       |       |
    !    |       |       |       |       |
    !    |       |       |       |       |
    !    =================================
    !            |       |
    !            |       |
    !            |       |
    !            |       |
    !            |       |
    !            =========
    !
    jx_min=fvm%jx_min; jx_max=fvm%jx_max; 
    jy_min=fvm%jy_min; jy_max=fvm%jy_max;
    !
    ! fvm%cubeboundary=0 means interior element
    !
    if (fvm%cubeboundary > 0) then
       !
       ! element is on panel side
       !
!       write(*,*) "element is on panel side"!dbg
       jx_min1=fvm%jx_min1; jx_max1=fvm%jx_max1; 
       jy_min1=fvm%jy_min1; jy_max1=fvm%jy_max1;
       swap1=fvm%swap1
       if (fvm%cubeboundary > 4) then
          !
          ! element is on a panel corner
          !
!          write(*,*) "element is on panel corner"!dbg
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
    do jy=1,nc+1
       do jx=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%faceno,dcart(jx,jy))  
       end do
    end do
    
    do jy=1, nc
       do jx=1, nc            
          !
          ! define departure cell
          !fvm%cubeboundary
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
          call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jx,jy,nreconstruction,&
               fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
               tmp,ngpc,gsweights,gspts,&
               weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
    end do
    

    !WEST SIDE
    if (fvm%cubeboundary == west) then
       !
       !    x--x--x--x--W-  jy=nc+1
       !    |  |  |  |  |
       !    x--x--x--x--W-  jy=nc
       !    |  |  |  |  | 
       !    x--x--x--x--W-  ...
       !    |  |  |  |  |
       !    x--x--x--x--W-  jy=2
       !    |  |  |  |  |
       !    x--x--x--x--W-  jy=1
       !
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face
       !
       jx=1
       do jy=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(west),dcart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
               fvm%nbrsface(west),dcart(jx+1,jy))                  
       end do
       do jy=1,nc
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
          
          if(swap1) then  !flip orientation
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else  
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==5) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,1,j)+nhe-1
                end do
             end do
          end if
          if (fvm%faceno==6) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=jy_max1-jy_min1-weights_eul_index_cell(i,2,j)-nhe-nhe+1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,:) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,:)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
    endif

    !EAST SIDE
    if (fvm%cubeboundary == east) then
       !
       ! no "swapping case"
       !
       !
       ! jy=nc+1  -E--x--x--x--x
       !           |  |  |  |  |
       ! jy=nc    -E--x--x--x--x
       !           |  |  |  |  |
       ! ...      -E--x--x--x--x
       !           |  |  |  |  |
       ! jy=2     -E--x--x--x--x
       !           |  |  |  |  |
       ! jy=1     -E--x--x--x--x
       !
       !        
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face 
       jx=nc
       do jy=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(east),dcart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
               fvm%nbrsface(east),dcart(jx+1,jy))                  
       end do

       do jy=1,nc
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
          if(swap1) then !flip orientation
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          !I have to correct the number
          if (fvm%faceno==5) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=jy_max1-jy_min1-weights_eul_index_cell(i,2,j)-nhe-nhe+1
                end do
             end do
          end if
          if (fvm%faceno==6) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,1,j)-nhe+1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
       
    endif

    !NORTH SIDE 
    if (fvm%cubeboundary == north) then
       !
       ! no "swapping case"
       !
       !
       !  jx=1   ...  jx=nc+1
       !     jx=2  jx=nc
       !
       !    N--N--N--N--N
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face
       !
       jy=nc
       do jx=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(north),dcart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
               fvm%nbrsface(north),dcart(jx,jy+1))                  
       end do

       do jx=1,nc
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
          if(swap1) then !flip orientation
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  
          end if
          if (fvm%faceno==2) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                end do
             end do
          end if
          if (fvm%faceno==3) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j))-nhe-nhe+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                end do
             end do
          end if
          if (fvm%faceno==4) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)-nhe-nhe+1
                end do
             end do
          end if
          if (fvm%faceno==6) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)-nhe-nhe+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
             
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
       
    endif
    !SOUTH SIDE
    if (fvm%cubeboundary == south) then
       !
       ! no "swapping case"
       !
       !
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    S--S--S--S--S
       ! 
       !  jx=1   ...  jx=nc+1
       !     jx=2  jx=nc
       !
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face
       !
       jy=1
       do jx=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(south),dcart(jx,jy))                   
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
               fvm%nbrsface(south),dcart(jx,jy+1))                   
       end do

       do jx=1,nc
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
          if(swap1) then !flip orientation
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  
          end if
          if  (fvm%faceno==2) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j)+nhe)-nhe+1
                end do
             end do
          end if
          if  (fvm%faceno==3) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j)+nhe)-nhe+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          if  (fvm%faceno==4) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          if  (fvm%faceno==5) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j)+nhe)-nhe+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
       
    endif
    !SOUTHWEST Corner
    if (fvm%cubeboundary == swest) then
       !
       ! no "swapping case"
       !
       !
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !   SW-SW-SW-SW-SW--
       !                |
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face
       jy=1
       do jx=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(south),dcart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
               fvm%nbrsface(south),dcart(jx,jy+1))                  
       end do

       do jx=1,nc
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
          if(swap1) then !flip orientation
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if ((fvm%faceno==3) .OR. (fvm%faceno==5)) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)+1
                   weights_eul_index_cell(i,2,j)=(weights_eul_index_cell(i,2,j)+nhe-1)
                end do
             end do
          end if
          if (fvm%faceno==2) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)+1
                end do
             end do
          end if
          if (fvm%faceno==4) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
       
       
       ! calculate xy Cartesian on the cube of departure points on the corresponding face
       jx=1
       do jy=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(west),dcart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
               fvm%nbrsface(west),dcart(jx+1,jy))                  
       end do

       do jy=1,nc
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
          if(swap2) then !flip orientation
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jy_min2,jx_min2,nreconstruction,&
                  fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==5) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,1,j)+nhe-1
                end do
             end do
          end if
          if (fvm%faceno==6) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2,j)+1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell (1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
    endif
    
    ! SOUTHEAST Corner
    if (fvm%cubeboundary == seast) then
       !
       ! no "swapping case"
       !
       !
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !  -SE-SE-SE-SE-SE
       !    |            
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face 
       !
       jy=1
       do jx=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(south),dcart(jx,jy))  
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
               fvm%nbrsface(south),dcart(jx,jy+1))  
       end do
       do jx=1,nc
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
          if(swap1) then !flip orientation
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          !I have to correct the number   
          if  (fvm%faceno==2) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j)+nhe)-nhe+1
                end do                   
             end do
          end if
          if  (fvm%faceno==3) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j)+nhe)-nhe+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          if  (fvm%faceno==4) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          if  (fvm%faceno==5) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j)+nhe)-nhe+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
       
       ! calculate xy Cartesian on the cube of departure points on the corresponding face  
       jx=nc
       do jy=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(east),dcart(jx,jy))                   
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
               fvm%nbrsface(east),dcart(jx+1,jy))                   
       end do

       do jy=1,nc
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
          if(swap2) then !flip orientation
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jy_min2,jx_min2,nreconstruction,&
                  fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
          end do
          else
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==5) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2,j)+1
                end do
             end do
          end if
          if (fvm%faceno==6) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,1,j)-nhe+1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
    endif
    
    !NORTHEAST Corner
    if (fvm%cubeboundary == neast) then
       !
       !    |
       ! --NE-NE-NE-NE-NE
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face  
       jy=nc
       do jx=1,nc+1
             call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
                  fvm%nbrsface(north),dcart(jx,jy))                
             call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
                  fvm%nbrsface(north),dcart(jx,jy+1))                
       end do

       do jx=1,nc
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
          if(swap1) then !flip orientation
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==2) then
             do j=1,2
                do i=1,jcollect_cell(j)
                weights_eul_index_cell(i,2,1)=(weights_eul_index_cell(i,2,j))-nhe+1
             end do
          end do
          end if
          if (fvm%faceno==3) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)-nhe-nhe+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                end do
             end do
          end if
          if (fvm%faceno==6) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)-nhe-nhe+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                end do
             end do
          end if
          if (fvm%faceno==4) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)-nhe-nhe+1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
       
       ! calculate xy Cartesian on the cube of departure points on the corresponding face  
       jx=nc
       do jy=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(east),dcart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
               fvm%nbrsface(east),dcart(jx+1,jy))                  
       end do
       jx=nc
       do jy=1,nc
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
          if(swap2) then !flip orientation
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jy_min2,jx_min2,nreconstruction,&
                  fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==5) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=jy_max2-jy_min2-weights_eul_index_cell(i,2,j)-nhe-nhe+1
                end do
             end do
          end if
          if (fvm%faceno==6) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,1,j)-nhe+1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
    endif
    
    !NORTH WEST CORNER 
    if (fvm%cubeboundary == nwest) then
       !
       !
       !                |
       !   NW-NW-NW-NW-NW--
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !    |  |  |  |  |
       !    x--x--x--x--x
       !
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face
       jy=nc
       do jx=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(north),dcart(jx,jy))                   
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
               fvm%nbrsface(north),dcart(jx,jy+1))                   
       end do
       do jx=1,nc
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
          if(swap1) then !flip orientation
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==2) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                end do
             end do
          end if
          if (fvm%faceno==3) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j))+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                end do
             end do
          end if
          if (fvm%faceno==4) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)+1
                end do
             end do
          end if
          if (fvm%faceno==6) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
       
       ! calculate xy Cartesian on the cube of departure points on the corresponding face
       jx=1
       do jy=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(west),dcart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
               fvm%nbrsface(west),dcart(jx+1,jy))                  
       end do
       do jy=1,nc
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
          if(swap2) then !flip orientation
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jy_min2,jx_min2,nreconstruction,&
                  fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(.true.,.true.,xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==5) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,1,j)+nhe-1
                end do
             end do
          end if
          if (fvm%faceno==6) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=jy_max2-jy_min2-weights_eul_index_cell(i,2,j)-nhe-nhe+1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)                 
             endif
          end do
       end do
    endif
    jall=jall-1
    
  end subroutine compute_weights_fluxform


  subroutine compute_weights_flux_cell(ldo_xflux,ldo_yflux,xcell_in,ycell_in,jx,jy,nreconstruction,xgno,ygno,&
       jx_min, jx_max, jy_min, jy_max,tmp,&
       ngauss,gauss_weights,abscissae,weights,weights_eul_index,jcollect,jmax_segments)

    use fvm_line_integrals_mod, only : compute_weights_cell
    implicit none
    logical, intent(in) :: ldo_xflux, ldo_yflux
    integer (kind=int_kind)                  , intent(in):: nreconstruction, jx,jy,ngauss,jmax_segments
    real (kind=real_kind)   ,  dimension(0:5), intent(inout):: xcell_in,ycell_in
    !
    integer (kind=int_kind), intent(in)               :: jx_min, jy_min, jx_max, jy_max
    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(inout) :: xgno
    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(inout) :: ygno
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
    integer (kind=int_kind), dimension(2), intent(out)       :: jcollect
    !
    ! local workspace
    !
    !
    ! max number of line segments is:
    !
    ! (number of longitudes)*(max average number of crossings per line segment = 3)*ncube*2
    !
    real (kind=real_kind)   ,  &
         dimension(jmax_segments,nreconstruction,2), intent(out) :: weights
    integer (kind=int_kind),  &
         dimension(jmax_segments,2,2), intent(out)      :: weights_eul_index
    !
    ! local workspace
    !
    integer :: i,j,nvertex,jcollect1,jcollect2
    real (kind=real_kind)   ,  dimension(0:5) :: xcell_flux,ycell_flux
    real (kind=real_kind)   ,  dimension(0:5) :: xcell_flux2,ycell_flux2
    real (kind=real_kind)                     :: weight_sign,weight_sign2
    logical :: lzero_flux
    real (kind=real_kind)   ,  &
         dimension(jmax_segments,nreconstruction,2) :: weights2

    ! dbg - start
    nvertex=4
    jcollect = 0
    call compute_weights_cell(nvertex,lexact_horizontal_line_integrals,&
         xcell_in(1:nvertex),ycell_in(1:nvertex),jx,jy,nreconstruction,xgno,ygno,&
         jx_min, jx_max, jy_min, jy_max,tmp,&
         ngauss,gauss_weights,abscissae,weights(:,:,1),weights_eul_index(:,:,1),jcollect(1),jmax_segments)
    call compute_weights_cell(nvertex,lexact_horizontal_line_integrals,&
         xcell_in(1:nvertex),ycell_in(1:nvertex),jx,jy,nreconstruction,xgno,ygno,&
         jx_min, jx_max, jy_min, jy_max,tmp,&
         ngauss,gauss_weights,abscissae,weights(:,:,2),weights_eul_index(:,:,2),jcollect(2),jmax_segments)
    return
    ! dbg - end

    lzero_flux = .false.
    ! dbg - start
    do i=-nhe,nc+2+nhe
       xgno(i)=i*1.0
       ygno(i)=i*1.0
    end do
    !
    ! make unit test for all cases - also epsilon cases!
    !

    xcell_in(1) = jx  ; ycell_in(1) = jy+0.5;
    xcell_in(2) = jx  ; ycell_in(2) = jy+1;
    xcell_in(3) = jx+1; ycell_in(3) = jy+1;
    xcell_in(4) = jx+1; ycell_in(4) = jy+0.5;

    ! dbg - end

    !
    ! figure out how to loop over flux-edges (in calling routine)
    !
    !
    !               gno(jx ,jy+1)-----------------gno(jx+1,jy+1)
    !              /|                             |  
    !             / |                             |  
    !            /  |                             |  
    ! cell(jx,jy+1) |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !           |   gno(jx ,jy  )-----------------gno(jx+1,jy  )
    !           |  /                             /
    !           | /                             /
    !           |/                             /
    !       cell(jx,jy)-------------------cell(jx+1,jy)
    !
    !    
    jcollect=-1
    if (ldo_xflux) then
       !
       ! constuct xflux-cell
       !
       xcell_flux(1) = xcell_in  (1 ); ycell_flux(1) = ycell_in  (1   ); 
       xcell_flux(2) = xcell_in  (2 ); ycell_flux(2) = ycell_in  (2   ); 
       xcell_flux(3) = xgno      (jx); ycell_flux(3) = ygno      (jy+1);
       xcell_flux(4) = xgno      (jx); ycell_flux(4) = ygno      (jy  );
       xcell_flux(0) = xcell_flux(4 ); ycell_flux(0) = ycell_flux(4   );
       xcell_flux(5) = xcell_flux(1 ); ycell_flux(5) = ycell_flux(1   );
       
       call make_flux_area(jx,jy,xgno,ygno,xcell_flux,ycell_flux,xcell_flux2,ycell_flux2,&
            lzero_flux,nvertex,weight_sign,weight_sign2,.true.)
       
       if (lzero_flux) then
          jcollect(1)=0
       else
          call compute_weights_cell(nvertex,lexact_horizontal_line_integrals,&
               xcell_flux(0:nvertex+1),ycell_flux(0:nvertex+1),jx,jy,nreconstruction,xgno,ygno,&
               jx_min, jx_max, jy_min, jy_max,tmp,&
               ngauss,gauss_weights,abscissae,weights(:,:,1),weights_eul_index,jcollect1,jmax_segments)
          weights(1:jcollect1,:,1) = weight_sign*weights(1:jcollect1,:,1)
          if (weight_sign2>-2) then
             !
             ! hour-glass flow situation
             !
             call compute_weights_cell(nvertex,lexact_horizontal_line_integrals,&
                  xcell_flux(0:nvertex+1),ycell_flux(0:nvertex+1),jx,jy,nreconstruction,xgno,ygno,&
                  jx_min, jx_max, jy_min, jy_max,tmp,&
                  ngauss,gauss_weights,abscissae,weights2(:,:,1),weights_eul_index,jcollect2,jmax_segments)
             weights(jcollect1+1:jcollect1+jcollect2,:,1) = weight_sign2*weights2(1:jcollect2,:,1)
             jcollect1 = jcollect1+jcollect2
          end if
       end if
       jcollect(1) = jcollect1
       !dbg start
       write(*,*) "zero xflux?",lzero_flux
       if (.not.lzero_flux) then
          WRITE(*,*) "xflux area cell       ",jx,jy
          WRITE(*,*) "xflux area cell - sign",weight_sign
          do i=1,nvertex
             write(*,*) i,xcell_flux(i),ycell_flux(i)
          end do
          WRITE(*,*) ""
          if (weight_sign2>-2.0) then
             WRITE(*,*) "xflux area2 cell       ",jx,jy
             WRITE(*,*) "xflux area2 cell - sign",weight_sign2
             do i=1,nvertex
                write(*,*) i,xcell_flux2(i),ycell_flux2(i)
             end do
             WRITE(*,*) ""
          end if
          !dbg end
       end if
    end if

    if (ldo_yflux) then
       !
       ! xxx jcollect should be separate for x nd y fluxes
       ! xxx idem for weights_eul_index
       !
       !
       ! constuct yflux-cell
       !
       xcell_flux(1) = xcell_in  (4   ); ycell_flux(1) = ycell_in  (4 ); 
       xcell_flux(2) = xcell_in  (1   ); ycell_flux(2) = ycell_in  (1 ); 
       xcell_flux(3) = xgno      (jx  ); ycell_flux(3) = ygno      (jy);
       xcell_flux(4) = xgno      (jx+1); ycell_flux(4) = ygno      (jy);
       xcell_flux(0) = xcell_flux(4   ); ycell_flux(0) = ycell_flux(4 );
       xcell_flux(5) = xcell_flux(1   ); ycell_flux(5) = ycell_flux(1 );
       
       call make_flux_area(jx,jy,xgno,ygno,xcell_flux,ycell_flux,xcell_flux2,ycell_flux2,&
            lzero_flux,nvertex,weight_sign,weight_sign2,.false.)
       
       if (lzero_flux) then
          jcollect(2)=0
       else
          call compute_weights_cell(nvertex,lexact_horizontal_line_integrals,&
               xcell_flux(0:nvertex+1),ycell_flux(0:nvertex+1),jx,jy,nreconstruction,xgno,ygno,&
               jx_min, jx_max, jy_min, jy_max,tmp,&
               ngauss,gauss_weights,abscissae,weights(:,:,2),weights_eul_index,jcollect1,jmax_segments)
          weights(1:jcollect1,:,2) = weight_sign*weights(1:jcollect1,:,2)
          if (weight_sign2>-2) then
             !
             ! hour-glass flow situation
             !
             call compute_weights_cell(nvertex,lexact_horizontal_line_integrals,&
                  xcell_flux(0:nvertex+1),ycell_flux(0:nvertex+1),jx,jy,nreconstruction,xgno,ygno,&
                  jx_min, jx_max, jy_min, jy_max,tmp,&
                  ngauss,gauss_weights,abscissae,weights2(:,:,2),weights_eul_index,jcollect2,jmax_segments)
             weights(jcollect1+1:jcollect1+jcollect2,:,2) = weight_sign2*weights2(1:jcollect2,:,2)
             jcollect1 = jcollect1+jcollect2
          end if
          jcollect(2) = jcollect1
       end if
       !dbg start
       write(*,*) "zero yflux?",lzero_flux
       if (.not.lzero_flux) then
          WRITE(*,*) "yflux area cell       ",jx,jy
          WRITE(*,*) "yflux area cell - sign",weight_sign
          do i=1,nvertex
             write(*,*) i,xcell_flux(i),ycell_flux(i)
          end do
          WRITE(*,*) ""
          if (weight_sign2>-2.0) then
             WRITE(*,*) "yflux area2 cell       ",jx,jy
             WRITE(*,*) "yflux area2 cell - sign",weight_sign2
             do i=1,nvertex
                write(*,*) i,xcell_flux2(i),ycell_flux2(i)
             end do
             WRITE(*,*) ""
          end if
       end if
          !dbg end
    end if
       
    stop!dbg
    
 
    weights(:,:,2) = -99.9E9

  end subroutine compute_weights_flux_cell

  subroutine make_flux_area(jx,jy,xgno,ygno,xcell_flux,ycell_flux,xcell_flux2,ycell_flux2,&
       lzero_flux,nvertex,weight_sign,weight_sign2,lx)
    use fvm_line_integrals_mod, only: compute_slope, y_cross_eul_lon, x_cross_eul_lat
    implicit none
    integer (kind=int_kind)              , intent(in   ):: jx,jy
    real (kind=real_kind), dimension(0:5), intent(inout):: xcell_flux   ,ycell_flux
    real (kind=real_kind), dimension(0:5), intent(  out):: xcell_flux2  ,ycell_flux2
    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in) :: xgno
    real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in) :: ygno
    logical                , intent(out) :: lzero_flux
    integer (kind=int_kind), intent(out):: nvertex
    real (kind=real_kind), intent(out) :: weight_sign,weight_sign2
    logical, intent(in) :: lx !.true. if xflux otherwise yflux
    !
    ! local workspace
    !
    logical :: ltrajec1_zero,ltrajec2_zero
    real (kind=real_kind), dimension(0:5) :: xcell_flux_tmp, ycell_flux_tmp
    real (kind=real_kind)                 :: slope, xcross, ycross
    !
    !
    !          cell_flux(3)------------------------
    !              /|                             |  
    !             / |                             |  
    !            /  |                             |  
    ! cell_flux(2)  |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !           |cell_flux(4)---------------------|
    !           |  /                             /
    !           | /                             /
    !           |/                             /
    !       cell_flux(1)----------------------/
    !
    
    !
    ! xxx this is different for yflux
    !
    weight_sign  = -999999.99
    weight_sign2 = -999999.99
    ltrajec1_zero=.false.; ltrajec2_zero=.false.
    if (lx) then
       if (abs(xcell_flux(1)-xcell_flux(4))<tiny) ltrajec1_zero=.true.
       if (abs(xcell_flux(2)-xcell_flux(3))<tiny) ltrajec2_zero=.true.
    else
       if (abs(ycell_flux(1)-ycell_flux(4))<tiny) ltrajec1_zero=.true.
       if (abs(ycell_flux(2)-ycell_flux(3))<tiny) ltrajec2_zero=.true.
    end if

    lzero_flux = (ltrajec1_zero.and.ltrajec2_zero)
    if (lzero_flux) then
       xcell_flux=-999999.99; ycell_flux=-999999.99; nvertex=-1
    else
       xcell_flux_tmp = xcell_flux; ycell_flux_tmp = ycell_flux;
       if (ltrajec1_zero) then
          !
          ! flux-area a triangle
          !
          nvertex = 3          
          call orient(xcell_flux(0:nvertex+1),ycell_flux(0:nvertex+1),nvertex,weight_sign)
       else if (ltrajec2_zero) then
          nvertex = 3
          xcell_flux(1:2) = xcell_flux_tmp(1:2); ycell_flux(1:2) = ycell_flux_tmp(1:2);
          xcell_flux(3  ) = xcell_flux_tmp(4  ); ycell_flux(3  ) = ycell_flux_tmp(4  );
          call orient(xcell_flux(0:nvertex+1),ycell_flux(0:nvertex+1),nvertex,weight_sign)
       else
          slope = compute_slope(xcell_flux(1:2),ycell_flux(1:2))
          if (lx) then
             ycross = y_cross_eul_lon(xcell_flux(1),ycell_flux(1),xcell_flux(4),slope)
             if (ycross+tiny>ycell_flux(4).and.ycross-tiny<ycell_flux(3)) then
                !
                ! hour-glass flow situation
                !
                nvertex = 3
                xcell_flux(1) = xcell_flux_tmp(1); ycell_flux(1) = ycell_flux_tmp(1);
                xcell_flux(2) = xcell_flux_tmp(4); ycell_flux(2) = ycross;
                xcell_flux(3) = xcell_flux_tmp(4); ycell_flux(3) = ycell_flux_tmp(4);
                call orient(xcell_flux(0:nvertex+1),ycell_flux(0:nvertex+1),nvertex,weight_sign)             
                !
                xcell_flux2(1) = xcell_flux_tmp(4); ycell_flux2(1) = ycross;
                xcell_flux2(2) = xcell_flux_tmp(2); ycell_flux2(2) = ycell_flux_tmp(2);
                xcell_flux2(3) = xcell_flux_tmp(3); ycell_flux2(3) = ycell_flux_tmp(3);
                call orient(xcell_flux2(0:nvertex+1),ycell_flux2(0:nvertex+1),nvertex,weight_sign2)             
             else
                !
                ! flux-area is a quadrilateral
                !
                nvertex = 4
                call orient(xcell_flux(0:nvertex+1),ycell_flux(0:nvertex+1),nvertex,weight_sign)
             end if
          else
             xcross = x_cross_eul_lat(xcell_flux(1),ycell_flux(1),ycell_flux(4),slope)
             if (xcross-tiny<xcell_flux(4).and.xcross+tiny>xcell_flux(3)) then
                !
                ! hour-glass flow situation
                !
                nvertex = 3
                xcell_flux(1) = xcell_flux_tmp(1); ycell_flux(1) = ycell_flux_tmp(1);
                xcell_flux(2) = xcross           ; ycell_flux(2) = ycell_flux_tmp(4);
                xcell_flux(3) = xcell_flux_tmp(4); ycell_flux(3) = ycell_flux_tmp(4);
                call orient(xcell_flux(0:nvertex+1),ycell_flux(0:nvertex+1),nvertex,weight_sign)             
                !
                xcell_flux2(1) = xcross           ; ycell_flux2(1) = ycell_flux_tmp(4);
                xcell_flux2(2) = xcell_flux_tmp(2); ycell_flux2(2) = ycell_flux_tmp(2);
                xcell_flux2(3) = xcell_flux_tmp(3); ycell_flux2(3) = ycell_flux_tmp(3);
                call orient(xcell_flux2(0:nvertex+1),ycell_flux2(0:nvertex+1),nvertex,weight_sign2)             
             else
                !
                ! flux-area is a quadrilateral
                !
                nvertex = 4
                call orient(xcell_flux(0:nvertex+1),ycell_flux(0:nvertex+1),nvertex,weight_sign)
             end if
          end if
       end if
    end if
    
  end subroutine make_flux_area

  real (kind=real_kind) function area(xseg,yseg,nvertex)
    implicit none
    real (kind=real_kind)  , dimension(nvertex), intent(in):: xseg,yseg
    integer (kind=int_kind)                    , intent(in):: nvertex
    !
    integer (kind=int_kind):: i
    area = 0.0
    do i=1,nvertex-1
       area = area + (yseg(i+1)-yseg(i))/(xseg(i+1)+xseg(i))
    end do
    area = area + (yseg(1)-yseg(nvertex))/(xseg(1)+xseg(nvertex))
    if (abs(area)< tiny) area = 0.0
  end function area

  subroutine orient(x,y,nvertex,weight_sign)
    implicit none
    real (kind=real_kind), dimension(0:nvertex+1), intent(inout):: x,y
    real (kind=real_kind)                        , intent(out  ):: weight_sign
    !
    real (kind=real_kind), dimension(0:nvertex+1) :: xtmp,ytmp
    integer (kind=int_kind)           , intent(in):: nvertex

    if (area(x(1:nvertex),y(1:nvertex),nvertex)<0) then
       xtmp(1:nvertex)=x(1:nvertex); ytmp(1:nvertex)=y(1:nvertex)

       x(1:nvertex) = xtmp(nvertex:1:-1); y(1:nvertex) = ytmp(nvertex:1:-1);
       x(0        ) = xtmp(nvertex     ); y(0        ) = ytmp(nvertex     );
       x(nvertex+1) = xtmp(1           ); y(nvertex+1) = ytmp(1           );
       weight_sign  = -1.0
    else
       xtmp(1:nvertex)=x(1:nvertex); ytmp(1:nvertex)=y(1:nvertex)
       weight_sign =  1.0
    end if

  end subroutine orient

end module fvm_line_integrals_flux_mod
