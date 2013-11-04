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
    integer (kind=int_kind), dimension(2)       :: inttmp
    real (kind=real_kind)                       :: tmp
    logical                                     :: swap
    ! for Gaussian quadrature
    real (kind=real_kind), dimension(ngpc)      :: gsweights, gspts
    ! weight-variables for individual cells
    integer (kind=int_kind) :: jmax_segments_cell
    real (kind=real_kind)   , dimension(nhe*50,nreconstruction,2)   :: weights_flux_cell
    integer (kind=int_kind),  dimension(nhe*50,2,2)                 :: weights_eul_index_cell
    integer (kind=int_kind)                                         :: jcollect_cell
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
       write(*,*) "element is on panel side"!dbg
       jx_min1=fvm%jx_min1; jx_max1=fvm%jx_max1; 
       jy_min1=fvm%jy_min1; jy_max1=fvm%jy_max1;
       swap1=fvm%swap1
       if (fvm%cubeboundary > 4) then
          !
          ! element is on a panel corner
          !
          write(*,*) "element is on panel corner"!dbg
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
    
    write(*,*) "doing interior from/to [1:nc]x[1:nc]",1,nc!dbg
    do jy=1, nc
       do jx=1, nc            
          !
          ! define departure cell
          !
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
          call compute_weights_flux_cell(xcell,ycell,jx,jy,nreconstruction,&
               fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
               tmp,ngpc,gsweights,gspts,&
               weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
          if (jcollect_cell>0) then
             weights_all(jall:jall+jcollect_cell-1,:,:) = weights_flux_cell(1:jcollect_cell,:,:)
             weights_eul_index_all(jall:jall+jcollect_cell-1,:,:) = &
                  weights_eul_index_cell(1:jcollect_cell,:,:)
             weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
             weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
             jall = jall+jcollect_cell          
          endif
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
             call compute_weights_flux_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1,:)
                weights_eul_index_cell(i,1,:)=weights_eul_index_cell(i,2,:)
                weights_eul_index_cell(i,2,:)=inttmp
             end do
          else  
             call compute_weights_flux_cell(xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==5) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=weights_eul_index_cell(i,1,1)+nhe-1
                weights_eul_index_cell(i,1,2)=weights_eul_index_cell(i,1,2)+nhe-1
             end do
          end if
          if (fvm%faceno==6) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,2,1)=jy_max1-jy_min1-weights_eul_index_cell(i,2,1)-nhe-nhe+1
                weights_eul_index_cell(i,2,2)=jy_max1-jy_min1-weights_eul_index_cell(i,2,2)-nhe-nhe+1
             end do
          end if
          if (jcollect_cell>0) then
             weights_all(jall:jall+jcollect_cell-1,:,:) = weights_flux_cell(1:jcollect_cell,:,:)
             weights_eul_index_all(jall:jall+jcollect_cell-1,:,:) = &
                  weights_eul_index_cell(1:jcollect_cell,:,:)
             weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
             weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
             jall = jall+jcollect_cell          
          endif
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
             call compute_weights_flux_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1,:)
                weights_eul_index_cell(i,1,:)=weights_eul_index_cell(i,2,:)
                weights_eul_index_cell(i,2,:)=inttmp
             end do
          else
             call compute_weights_flux_cell(xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          !I have to correct the number
          if (fvm%faceno==5) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,2,1)=jy_max1-jy_min1-weights_eul_index_cell(i,2,1)-nhe-nhe+1
                weights_eul_index_cell(i,2,2)=jy_max1-jy_min1-weights_eul_index_cell(i,2,2)-nhe-nhe+1
             end do
          end if
          if (fvm%faceno==6) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=weights_eul_index_cell(i,1,1)-nhe+1
                weights_eul_index_cell(i,1,2)=weights_eul_index_cell(i,1,2)-nhe+1
             end do
          end if
          if (jcollect_cell>0) then
             weights_all(jall:jall+jcollect_cell-1,:,:) = weights_flux_cell(1:jcollect_cell,:,:)
             weights_eul_index_all(jall:jall+jcollect_cell-1,:,:) = &
                  weights_eul_index_cell(1:jcollect_cell,:,:)
             weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
             weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
             jall = jall+jcollect_cell          
          endif
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
             call compute_weights_flux_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1,:)
                weights_eul_index_cell(i,1,:)=weights_eul_index_cell(i,2,:)
                weights_eul_index_cell(i,2,:)=inttmp
             end do
          else
             call compute_weights_flux_cell(xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  
          end if
          if (fvm%faceno==2) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)-nhe+1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)-nhe+1
             end do
          end if
          if (fvm%faceno==3) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,1))-nhe-nhe+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,2))-nhe-nhe+1
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)-nhe+1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)-nhe+1
             end do
          end if
          if (fvm%faceno==4) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,1)-nhe-nhe+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,2)-nhe-nhe+1
             end do
          end if
          if (fvm%faceno==6) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,1)-nhe-nhe+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,2)-nhe-nhe+1
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)-nhe+1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)-nhe+1
             end do
          end if
          if (jcollect_cell>0) then
             weights_all(jall:jall+jcollect_cell-1,:,:) = weights_flux_cell(1:jcollect_cell,:,:)
             weights_eul_index_all(jall:jall+jcollect_cell-1,:,:) = &
                  weights_eul_index_cell(1:jcollect_cell,:,:)
             weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
             weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
             jall = jall+jcollect_cell          
          endif
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
             call compute_weights_flux_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1,:)
                weights_eul_index_cell(i,1,:)=weights_eul_index_cell(i,2,:)
                weights_eul_index_cell(i,2,:)=inttmp
             end do
          else
             call compute_weights_flux_cell(xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  
          end if
          if  (fvm%faceno==2) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,1)+nhe)-nhe+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,2)+nhe)-nhe+1
             end do
          end if
          if  (fvm%faceno==3) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,1)+nhe)-nhe+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,2)+nhe)-nhe+1
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)+nhe-1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)+nhe-1
             end do
          end if
          if  (fvm%faceno==4) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)+nhe-1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)+nhe-1
             end do
          end if
          if  (fvm%faceno==5) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,1)+nhe)-nhe+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,2)+nhe)-nhe+1
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)+nhe-1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)+nhe-1
             end do
          end if
          if (jcollect_cell>0) then
             weights_all(jall:jall+jcollect_cell-1,:,:) = weights_flux_cell(1:jcollect_cell,:,:)
             weights_eul_index_all(jall:jall+jcollect_cell-1,:,:) = &
                  weights_eul_index_cell(1:jcollect_cell,:,:)
             weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
             weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
             jall = jall+jcollect_cell          
          endif
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
             call compute_weights_flux_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1,:)
                weights_eul_index_cell(i,1,:)=weights_eul_index_cell(i,2,:)
                weights_eul_index_cell(i,2,:)=inttmp
             end do
          else
             call compute_weights_flux_cell(xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if ((fvm%faceno==3) .OR. (fvm%faceno==5)) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,1)+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,2)+1
                weights_eul_index_cell(i,2,1)=(weights_eul_index_cell(i,2,1)+nhe-1)
                weights_eul_index_cell(i,2,2)=(weights_eul_index_cell(i,2,2)+nhe-1)
             end do
          end if
          if (fvm%faceno==2) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,1)+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,2)+1
             end do
          end if
          if (fvm%faceno==4) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)+nhe-1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)+nhe-1
             end do
          end if
          if (jcollect_cell>0) then
             weights_all(jall:jall+jcollect_cell-1,:,:) = weights_flux_cell(1:jcollect_cell,:,:)
             weights_eul_index_all(jall:jall+jcollect_cell-1,:,:) = &
                  weights_eul_index_cell(1:jcollect_cell,:,:)
             weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
             weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
             jall = jall+jcollect_cell          
          endif
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
             call compute_weights_flux_cell(xcell,ycell,jy_min2,jx_min2,nreconstruction,&
                  fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1,:)
                weights_eul_index_cell(i,1,:)=weights_eul_index_cell(i,2,:)
                weights_eul_index_cell(i,2,:)=inttmp
             end do
          else
             call compute_weights_flux_cell(xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==5) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=weights_eul_index_cell(i,1,1)+nhe-1
                weights_eul_index_cell(i,1,2)=weights_eul_index_cell(i,1,2)+nhe-1
             end do
          end if
          if (fvm%faceno==6) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,2,1)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2,1)+1
                weights_eul_index_cell(i,2,2)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2,2)+1
             end do
          end if
          if (jcollect_cell>0) then
             weights_all(jall:jall+jcollect_cell-1,:,:) = weights_flux_cell (1:jcollect_cell,:,:)
             weights_eul_index_all(jall:jall+jcollect_cell-1,:,:) = &
                  weights_eul_index_cell(1:jcollect_cell,:,:)
             weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
             weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
             jall = jall+jcollect_cell          
          endif
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
             call compute_weights_flux_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1,:)
                weights_eul_index_cell(i,1,:)=weights_eul_index_cell(i,2,:)
                weights_eul_index_cell(i,2,:)=inttmp
             end do
          else
             call compute_weights_flux_cell(xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          !I have to correct the number   
          if  (fvm%faceno==2) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,1)+nhe)-nhe+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,2)+nhe)-nhe+1
             end do
          end if
          if  (fvm%faceno==3) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,1)+nhe)-nhe+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,2)+nhe)-nhe+1
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)+nhe-1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)+nhe-1
             end do
          end if
          if  (fvm%faceno==4) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)+nhe-1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)+nhe-1
             end do
          end if
          if  (fvm%faceno==5) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,1)+nhe)-nhe+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,2)+nhe)-nhe+1
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)+nhe-1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)+nhe-1
             end do
          end if
          if (jcollect_cell>0) then
             weights_all(jall:jall+jcollect_cell-1,:,:) = weights_flux_cell(1:jcollect_cell,:,:)
             weights_eul_index_all(jall:jall+jcollect_cell-1,:,:) = &
                  weights_eul_index_cell(1:jcollect_cell,:,:)
             weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
             weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
             jall = jall+jcollect_cell          
          endif
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
             call compute_weights_flux_cell(xcell,ycell,jy_min2,jx_min2,nreconstruction,&
                  fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1,:)
                weights_eul_index_cell(i,1,:)=weights_eul_index_cell(i,2,:)
                weights_eul_index_cell(i,2,:)=inttmp
             end do
          else
             call compute_weights_flux_cell(xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==5) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,2,1)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2,1)+1
                weights_eul_index_cell(i,2,2)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2,2)+1
             end do
          end if
          if (fvm%faceno==6) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=weights_eul_index_cell(i,1,1)-nhe+1
                weights_eul_index_cell(i,1,2)=weights_eul_index_cell(i,1,2)-nhe+1
             end do
          end if
          if (jcollect_cell>0) then
             weights_all(jall:jall+jcollect_cell-1,:,:) = weights_flux_cell(1:jcollect_cell,:,:)
             weights_eul_index_all(jall:jall+jcollect_cell-1,:,:) = &
                  weights_eul_index_cell(1:jcollect_cell,:,:)
             weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
             weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
             jall = jall+jcollect_cell          
          endif
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
             call compute_weights_flux_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1,:)
                weights_eul_index_cell(i,1,:)=weights_eul_index_cell(i,2,:)
                weights_eul_index_cell(i,2,:)=inttmp
             end do
          else
             call compute_weights_flux_cell(xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==2) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,2,1)=(weights_eul_index_cell(i,2,1))-nhe+1
                weights_eul_index_cell(i,2,2)=(weights_eul_index_cell(i,2,2))-nhe+1
             end do
          end if
          if (fvm%faceno==3) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,1)-nhe-nhe+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,2)-nhe-nhe+1
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)-nhe+1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)-nhe+1
             end do
          end if
          if (fvm%faceno==6) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,1)-nhe-nhe+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,2)-nhe-nhe+1
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)-nhe+1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)-nhe+1
             end do
          end if
          if (fvm%faceno==4) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,1)-nhe-nhe+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,2)-nhe-nhe+1
             end do
          end if
          if (jcollect_cell>0) then
             weights_all(jall:jall+jcollect_cell-1,:,:) = weights_flux_cell(1:jcollect_cell,:,:)
             weights_eul_index_all(jall:jall+jcollect_cell-1,:,:) = &
                  weights_eul_index_cell(1:jcollect_cell,:,:)
             weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
             weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
             jall = jall+jcollect_cell          
          endif
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
             call compute_weights_flux_cell(xcell,ycell,jy_min2,jx_min2,nreconstruction,&
                  fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1,:)
                weights_eul_index_cell(i,1,:)=weights_eul_index_cell(i,2,:)
                weights_eul_index_cell(i,2,:)=inttmp
             end do
          else
             call compute_weights_flux_cell(xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==5) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,2,1)=jy_max2-jy_min2-weights_eul_index_cell(i,2,1)-nhe-nhe+1
                weights_eul_index_cell(i,2,2)=jy_max2-jy_min2-weights_eul_index_cell(i,2,2)-nhe-nhe+1
             end do
          end if
          if (fvm%faceno==6) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=weights_eul_index_cell(i,1,1)-nhe+1
                weights_eul_index_cell(i,1,2)=weights_eul_index_cell(i,1,2)-nhe+1
             end do
          end if
          
          if (jcollect_cell>0) then
             weights_all(jall:jall+jcollect_cell-1,:,:) = weights_flux_cell(1:jcollect_cell,:,:)
             weights_eul_index_all(jall:jall+jcollect_cell-1,:,:) = &
                  weights_eul_index_cell(1:jcollect_cell,:,:)
             weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
             weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
             jall = jall+jcollect_cell          
          endif
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
             call compute_weights_flux_cell(xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1,:)
                weights_eul_index_cell(i,1,:)=weights_eul_index_cell(i,2,:)
                weights_eul_index_cell(i,2,:)=inttmp
             end do
          else
             call compute_weights_flux_cell(xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==2) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)-nhe+1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)-nhe+1
             end do
          end if
          if (fvm%faceno==3) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,1))+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,2))+1
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)-nhe+1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)-nhe+1
             end do
          end if
          if (fvm%faceno==4) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,1)+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,2)+1
             end do
          end if
          if (fvm%faceno==6) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,1)+1
                weights_eul_index_cell(i,1,2)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,2)+1
                weights_eul_index_cell(i,2,1)=weights_eul_index_cell(i,2,1)-nhe+1
                weights_eul_index_cell(i,2,2)=weights_eul_index_cell(i,2,2)-nhe+1
             end do
          end if
          if (jcollect_cell>0) then
             weights_all(jall:jall+jcollect_cell-1,:,:) = weights_flux_cell(1:jcollect_cell,:,:)
             weights_eul_index_all(jall:jall+jcollect_cell-1,:,:) = &
                  weights_eul_index_cell(1:jcollect_cell,:,:)
             weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
             weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
             jall = jall+jcollect_cell          
          endif
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
             call compute_weights_flux_cell(xcell,ycell,jy_min2,jx_min2,nreconstruction,&
                  fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1,:)
                weights_eul_index_cell(i,1,:)=weights_eul_index_cell(i,2,:)
                weights_eul_index_cell(i,2,:)=inttmp
             end do
          else
             call compute_weights_flux_cell(xcell,ycell,jx,jy,nreconstruction,&
                  fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          if (fvm%faceno==5) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,1,1)=weights_eul_index_cell(i,1,1)+nhe-1
                weights_eul_index_cell(i,1,2)=weights_eul_index_cell(i,1,2)+nhe-1
             end do
          end if
          if (fvm%faceno==6) then
             do i=1,jcollect_cell
                weights_eul_index_cell(i,2,1)=jy_max2-jy_min2-weights_eul_index_cell(i,2,1)-nhe-nhe+1
                weights_eul_index_cell(i,2,2)=jy_max2-jy_min2-weights_eul_index_cell(i,2,2)-nhe-nhe+1
             end do
          end if
          if (jcollect_cell>0) then
             weights_all(jall:jall+jcollect_cell-1,:,:) = weights_flux_cell(1:jcollect_cell,:,:)
             weights_eul_index_all(jall:jall+jcollect_cell-1,:,:) = &
                  weights_eul_index_cell(1:jcollect_cell,:,:)
             weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
             weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
             jall = jall+jcollect_cell          
          endif
       end do
    endif
    
    jall=jall-1
    
  end subroutine compute_weights_fluxform


  subroutine compute_weights_flux_cell(xcell_in,ycell_in,jx,jy,nreconstruction,xgno,ygno,&
       jx_min, jx_max, jy_min, jy_max,tmp,&
       ngauss,gauss_weights,abscissae,weights,weights_eul_index,jcollect,jmax_segments)

    use fvm_line_integrals_mod, only : compute_weights_cell
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
         dimension(jmax_segments,nreconstruction,2), intent(out) :: weights
    integer (kind=int_kind),  &
         dimension(jmax_segments,2,2), intent(out)      :: weights_eul_index
    integer :: i,j

    
!    write(*,*) "y isolines"
!    do j=jy_min,jy_max
!       do i=jx_min,jx_max-1
!          write(*,*) xgno(i  ),ygno(j)
!          write(*,*) xgno(i+1),ygno(j)
!          write(*,*) "  "
!       end do
!    end do
!    write(*,*) "x isolines"
!    do j=jy_min,jy_max-1
!       do i=jx_min,jx_max
!          write(*,*) xgno(i),ygno(j)
!          write(*,*) xgno(i),ygno(j+1)
!          write(*,*) "  "
!       end do
!    end do
!    stop

    !
    ! figure out how to loop over flux-edges (in calling routine)
    !
    ! pseudo-code:
    !
    ! if (zero_flux_area) return 
    !
    ! zero_flux should be determined from trajectories; can not use area since flux-area may be bowtie
    !
    ! if (triangle flux_area) nvertex=3
    !
    ! triangle flux-area determined from trajctories
    !
    ! if (non_self_intersecting)
    !
    ! self-intersection can be computed through line-integral area computations!
    !
    !
    !
    !


    call compute_weights_cell(4,lexact_horizontal_line_integrals,&
         xcell_in,ycell_in,jx,jy,nreconstruction,xgno,ygno,&
         jx_min, jx_max, jy_min, jy_max,tmp,&
         ngauss,gauss_weights,abscissae,weights(:,:,1),weights_eul_index,jcollect,jmax_segments)


    weights(:,:,2) = -99.9E9




  end subroutine compute_weights_flux_cell
end module fvm_line_integrals_flux_mod
