#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!MODULE FVM_CONTROL_VOLUME_MOD---------------------------------------------CE-for FVM
! AUTHOR: Christoph Erath, 11.June 2011                                             !
! This module contains everything to initialize the arrival. It also provides the   !
! interpolation points for the reconstruction (projection from one face to another  !
! when the element is on the cube edge)                                             !
! It also intialize the start values, see also fvm_analytic                         !
!-----------------------------------------------------------------------------------! 
module fvm_control_volume_mod
  ! ---------------------------------------------------------------------------------
  use kinds, only : real_kind, int_kind
  ! ---------------------------------------------------------------------------------
  use coordinate_systems_mod, only : spherical_polar_t, cartesian2D_t
  ! ---------------------------------------------------------------------------------
  use element_mod, only: timelevels, element_t
  ! ---------------------------------------------------------------------------------
  use dimensions_mod, only: nc, nhc, nhe, nlev, ntrac, ntrac_d, ne, np
  ! ---------------------------------------------------------------------------------
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
#ifndef MESH
  ! ---------------------------------------------------------------------------------
  use cube_mod, only     : cube_xstart, cube_xend, cube_ystart, cube_yend
#endif

  use parallel_mod, only : abortmp
  
  implicit none
  private
  
  type, public :: fvm_struct
    sequence
    ! fvm tracer mixing ratio: (kg/kg)
    real (kind=real_kind) :: c(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev,ntrac_d,timelevels) 
    real (kind=real_kind) :: cstart(1:nc,1:nc)
!-----------------------------------------------------------------------------------!
    ! define the arrival grid, which is build on the original HOMME elements
    type (spherical_polar_t) :: asphere(nc+1,nc+1) ! spherical coordinates
    ! save velocity at time t 
    real (kind=real_kind)    :: vn0(np,np,2,nlev)
!-----------------------------------------------------------------------------------!
    real (kind=real_kind)    :: dalpha, dbeta    ! grid step in alpha/beta coord.
!-----------------------------------------------------------------------------------!
    ! define the departure grid (depends on the examples, usually the velocity)
    ! they are NOT stored for the Levels
    type (spherical_polar_t) :: dsphere(nc+1,nc+1)     ! Spherical coordinates
    type (spherical_polar_t) :: centersphere(nc,nc)      ! Spherical coordinates of fvm grid
!-----------------------------------------------------------------------------------!  
    !area of the cell on the sphere 
    real (kind=real_kind)    :: area_sphere(1-nhe:nc+nhe,1-nhe:nc+nhe)    
    real (kind=real_kind)    :: area_spherenew(1-nhe:nc+nhe,1-nhe:nc+nhe)    
    
    real (kind=real_kind)    :: spherecentroid(5,1-nhe:nc+nhe,1-nhe:nc+nhe)    
    integer                  :: faceno
    integer                  :: nbrsface(8)    ! store the neighbours in north, south 
    ! number of south,....,swest and 0 for interior element 
    integer                  :: cubeboundary                                                 
!-----------------------------------------------------------------------------------!     
    ! next maybe used just for test reasons
    real (kind=real_kind)    :: elem_mass
    real (kind=real_kind)    :: maxcfl(2,nlev)
!-----------------------------------------------------------------------------------!     
    ! data for arrival grid   
    integer                  :: jx_min, jx_max, jy_min, jy_max
    integer                  :: jx_min1, jx_max1, jy_min1, jy_max1
    integer                  :: jx_min2, jx_max2, jy_min2, jy_max2
    !gnomonic coordinates, the real Cartesian coordinates on the cube
    !first and last entry of acartx* is bignum because of search algorithm
    real (kind=real_kind)    :: acartx(-nhe:nc+2+nhe), acarty(-nhe:nc+2+nhe)
    ! points for the first halo zone (only for west, east, north, south, swest,
    !                                 seast, neast, nwest)
    real (kind=real_kind)    :: acartx1(-nhe:nc+2+nhe), acarty1(-nhe:nc+2+nhe)
    ! points for the second halo zone (only for swest, seast, neast, nwest)
    real (kind=real_kind)    :: acartx2(-nhe:nc+2+nhe), acarty2(-nhe:nc+2+nhe)
    ! tells us if orientation is swapped, or if we have to reorder it
    logical                  :: swap1, swap2
    ! for the reconstruction, points can be in the wrong order, derivative has
    ! wrong sign
    integer                  :: invx1, invy1, invx2,invy2
!-----------------------------------------------------------------------------------!     
    ! provide fixed interpolation points with respect to the arrival grid for 
    ! reconstruction   
    real (kind=real_kind)    :: interp(-1:nc+2,2,2)
    real (kind=real_kind)    :: interphalo(-1:nc+2,2,2)
    real (kind=real_kind)    :: interphaloex(-1:nc+2,2,2)
    integer                  :: ibase(-1:nc+2,2,2)
    integer                  :: ibasehalo(-1:nc+2,2,2)
    integer                  :: ibasehaloex(-1:nc+2,2,2)     
  end type fvm_struct

  public :: fvm_mesh_ari
  
  real (kind=real_kind),parameter, public   :: bignum = 1.0D20
  
contains

! ----------------------------------------------------------------------------------!
!SUBROUTINE fvm_MESH_ARI--------------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 11.June 2011                                             !
! DESCRIPTION: Calculates the mesh for fvm, use the elem structure from HOMME       !
!                                                                                   !
! CALLS: create_ari, create_interpolation_points                                    !
! INPUT: elem     ...  element structure from HOMME                                 !
!        tl       ...  time level structure                                         !
! INTPUT/OUTPUT:                                                                    !
!        fvm   ...  structure                                                       !
!-----------------------------------------------------------------------------------!
subroutine fvm_mesh_ari(elem, fvm, tl)
  use time_mod, only : timelevel_t
  implicit none
  type (element_t), intent(in)      :: elem
  type (fvm_struct), intent(inout)   :: fvm
  type (timeLevel_t), intent(in)       :: tl              ! time level struct
  
  integer                              :: i,j
  logical                              :: corner

  fvm%faceno=elem%FaceNum
  ! write the neighbors in the structure
  fvm%cubeboundary=0
  corner=.FALSE.

! Jose Garcia: This code does not work with MESH
! yet we allow it so MESH and fvm can be compiled 
! in the same executable. Some execution paths that
! do not call this routine will work fine.
#ifndef MESH
  do j=1,8
    if (elem%vertex%nbrs(j)%used) then
      fvm%nbrsface(j)=elem%vertex%nbrs(j)%f
      ! note that if the element lies on a corner, it will be at j=5,6,7,8
      if ((fvm%nbrsface(j) /= fvm%faceno) .AND. (j<5)) then
        fvm%cubeboundary=j
      endif
    else   ! corner on the cube
      if (.NOT. corner) then
        fvm%nbrsface(j)=-1
        fvm%cubeboundary=j
        corner=.TRUE.
      else
        print *,'Error in fvm_CONTROL_VOLUME_MOD - Subroutine fvm_MESH_ARI: '
        call abortmp('Do not allow one element per face for fvm, please increase ne!')
      endif
    end if
  end do
#endif

  call create_ari(elem,fvm)
  call create_interpolation_points(elem,fvm)

end subroutine fvm_mesh_ari
!END SUBROUTINE fvm_MESH_ARI----------------------------------------------CE-for FVM

! ----------------------------------------------------------------------------------!
!SUBROUTINE CREATE_ARI----------------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 7.November 2011                                          !
! DESCRIPTION: Calculates the mesh for fvm, including for the halo zone             !
!                                                                                   !
! INPUT: elem     ...  element structure from HOMME                                 ! 
! INTPUT/OUTPUT:                                                                    !
!        fvm   ...  structure                                                       !
!-----------------------------------------------------------------------------------!
subroutine create_ari(elem, fvm)
  use coordinate_systems_mod, only : cart2spherical  
  implicit none
  type (element_t), intent(in)      :: elem
  type (fvm_struct), intent(inout)   :: fvm
  
  integer                              :: i,j, jx, jy, jtmp
  logical                              :: corner
  integer                              :: jx_min, jx_max, jy_min, jy_max
  integer                              :: jx_min1, jx_max1, jy_min1, jy_max1
  integer                              :: jx_min2, jx_max2, jy_min2, jy_max2
  integer                              :: invx1, invy1,invx2, invy2
  logical                              :: swap1, swap2
  real (kind=real_kind)                :: centerx,centery
  
  swap1=.FALSE.
  swap2=.FALSE.
  invx1=1 ! 1 order the same, -1 reorder in x direction
  invy1=1 ! 1 order the same, -1 reorder in y direction
  invx2=1 ! 1 order the same, -1 reorder in x direction
  invy2=1 ! 1 order the same, -1 reorder in y direction
! fvm mesh is (alpha,beta) equidistant
  fvm%dalpha=abs(elem%corners(1)%x-elem%corners(2)%x)/nc !in alpha 
  fvm%dbeta=abs(elem%corners(1)%y-elem%corners(4)%y)/nc  !in beta
  
  do i=1,nc+1
    fvm%acartx(i) = tan(elem%corners(1)%x+(i-1)*fvm%dalpha) 
    fvm%acarty(i) = tan(elem%corners(1)%y+(i-1)*fvm%dbeta)
  end do

  !
  do j=1,nc
    do i=1,nc
      centerx = tan(elem%corners(1)%x+(i-0.5)*fvm%dalpha)  
      centery = tan(elem%corners(1)%y+(j-0.5)*fvm%dbeta) 
      fvm%centersphere(i,j) = &
            cart2spherical(centerx,centery,fvm%faceno)
           
      fvm%asphere(i,j) = cart2spherical(fvm%acartx(i),fvm%acarty(j),fvm%faceno)
    enddo
    ! values for the nodes on the edges
    fvm%asphere(nc+1,j) = cart2spherical(fvm%acartx(nc+1),fvm%acarty(j),fvm%faceno)
    fvm%asphere(j,nc+1) = cart2spherical(fvm%acartx(j),fvm%acarty(nc+1),fvm%faceno)  
  enddo  
  ! value for the upper right corner
  fvm%asphere(nc+1,nc+1) = cart2spherical(fvm%acartx(nc+1),fvm%acarty(nc+1),fvm%faceno) 
  
  select case (fvm%cubeboundary)
  ! element is completly inside a face, note that this case will be the most one
  ! INTERIOR ELEMENT
    case (0) 
      jx_min=1-nhe
      jx_max=nc+1+nhe
      jy_min=1-nhe
      jy_max=nc+1+nhe
      !fill in halo zone
      do i=0,nhe-1
        fvm%acartx(-i) = tan(elem%corners(1)%x-(i+1)*fvm%dalpha) 
        fvm%acartx(nc+2+i) = tan(elem%corners(1)%x+(nc+1+i)*fvm%dalpha) 
        fvm%acarty(-i) = tan(elem%corners(1)%y-(i+1)*fvm%dbeta)
        fvm%acarty(nc+2+i) = tan(elem%corners(1)%y+(nc+1+i)*fvm%dbeta)
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
  !WEST EDGE ON THE CUBE SIDE
    case (west)  
      jx_min=1
      jx_max=nc+1+nhe
      jy_min=1-nhe
      jy_max=nc+1+nhe
      do i=0,nhe-1
        fvm%acartx(nc+2+i) = tan(elem%corners(1)%x+(nc+1+i)*fvm%dalpha) 
        fvm%acarty(-i) = tan(elem%corners(1)%y-(i+1)*fvm%dbeta)
        fvm%acarty(nc+2+i) = tan(elem%corners(1)%y+(nc+1+i)*fvm%dbeta)
      enddo
      jx_min1=1-nhe
      jx_max1=1
      jy_min1=1-nhe
      jy_max1=nc+1+nhe
      ! use symmetry of the cubed-sphere on the edge
      if (fvm%faceno<=4) then
        fvm%acarty1(jy_min1:jy_max1)=fvm%acarty(jy_min1:jy_max1)
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=-fvm%acartx(1+(jx_max1-jx))   
        end do
      end if
      if (fvm%faceno==5) then
        ! flip orientation
        swap1=.TRUE.
        invx1=-1
        fvm%acarty1(jy_min1:jy_max1)=fvm%acarty(jy_min1:jy_max1)
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=fvm%acartx(1+jx-jx_min1)   
        end do
      end if
      ! use symmetry of the cubed-sphere on the edge
      if (fvm%faceno==6) then
        ! flip orientation
        swap1=.TRUE.
        invy1=-1
        do jy=jy_min1,jy_max1
          jtmp=jy_max1+jy_min1-jy
          fvm%acarty1(jtmp)=-fvm%acarty(jy)
        end do
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=-fvm%acartx(1+jx_max1-jx)  
        end do 
      end if
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      
  !EAST EDGE ON THE CUBE SIDE
    case(east)
      jx_min=1-nhe
      jx_max=nc+1
      jy_min=1-nhe
      jy_max=nc+1+nhe 
      do i=0,nhe-1
        fvm%acartx(-i) = tan(elem%corners(1)%x-(i+1)*fvm%dalpha) 
        fvm%acarty(-i) = tan(elem%corners(1)%y-(i+1)*fvm%dbeta)
        fvm%acarty(nc+2+i) = tan(elem%corners(1)%y+(nc+1+i)*fvm%dbeta)
      enddo
      jx_min1=nc+1
      jx_max1=nc+1+nhe
      jy_min1=1-nhe
      jy_max1=nc+1+nhe 
      if (fvm%faceno<=4) then
        fvm%acarty1(jy_min1:jy_max1)=fvm%acarty(jy_min1:jy_max1)
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=-fvm%acartx(nc+1-(jx-jx_min1))   
        end do
        
      end if
      if (fvm%faceno==5) then
        ! flip orientation
        swap1=.TRUE.
        invy1=-1
        do jy=jy_min1,jy_max1
          jtmp=jy_max1+jy_min1-jy
          fvm%acarty1(jtmp)=-fvm%acarty(jy)
        end do
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=-fvm%acartx(nc+1-(jx-jx_min1))  
        end do 
      end if
      ! use symmetry of the cubed-sphere on the edge
      if (fvm%faceno==6) then
        ! flip orientation
        swap1=.TRUE.
        invx1=-1
        fvm%acarty1(jy_min1:jy_max1)=fvm%acarty(jy_min1:jy_max1)
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=fvm%acartx(nc+1-(jx_max1-jx))  
        end do 
      end if

  !NORTH EDGE ON THE CUBE SIDE        
    case(north)
      jx_min=1-nhe
      jx_max=nc+1+nhe
      jy_min=1-nhe
      jy_max=nc+1
      do i=0,nhe-1
        fvm%acartx(-i) = tan(elem%corners(1)%x-(i+1)*fvm%dalpha) 
        fvm%acartx(nc+2+i) = tan(elem%corners(1)%x+(nc+1+i)*fvm%dalpha) 
        fvm%acarty(-i) = tan(elem%corners(1)%y-(i+1)*fvm%dbeta)
      enddo
      jx_min1=1-nhe
      jx_max1=nc+1+nhe
      jy_min1=nc+1
      jy_max1=nc+1+nhe
      ! use symmetry of the cubed-sphere on the edge
      if ((fvm%faceno==1) .OR. (fvm%faceno==5)) then
        fvm%acartx1(jx_min1:jx_max1)=fvm%acartx(jx_min1:jx_max1)
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=-fvm%acarty(nc+1-(jy-jy_min1))   
        end do 
      end if
      ! use symmetry of the cubed-sphere on the edge
      if ((fvm%faceno==3) .OR. (fvm%faceno==6)) then
        invx1=-1
        invy1=-1
        do jx=jx_min1,jx_max1
          jtmp=jx_max1+jx_min1-jx
          fvm%acartx1(jtmp)=-fvm%acartx(jx)
        end do
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(nc+1-(jy_max1-jy))  
        end do 
      end if
      if (fvm%faceno==2) then
        ! flip orientation
        swap1=.TRUE.
        invy1=-1
        fvm%acartx1(jx_min1:jx_max1)=fvm%acartx(jx_min1:jx_max1)
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(nc+1-(jy_max1-jy))   
        end do 
      end if
      if (fvm%faceno==4) then
        ! flip orientation
        swap1=.TRUE.
        invx1=-1
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=-fvm%acarty(nc+1-(jy-jy_min1))  
        end do 
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=-fvm%acartx(jx_max1-(jx-jx_min1))
        end do
      end if
      
  !SOUTH EDGE ON THE CUBE SIDE    
    case(south)
      jx_min=1-nhe
      jx_max=nc+1+nhe
      jy_min=1
      jy_max=nc+1+nhe
      do i=0,nhe-1
        fvm%acartx(-i) = tan(elem%corners(1)%x-(i+1)*fvm%dalpha) 
        fvm%acartx(nc+2+i) = tan(elem%corners(1)%x+(nc+1+i)*fvm%dalpha) 
        fvm%acarty(nc+2+i) = tan(elem%corners(1)%y+(nc+1+i)*fvm%dbeta)
      enddo
      jx_min1=1-nhe
      jx_max1=nc+1+nhe
      jy_min1=1-nhe
      jy_max1=1
      ! use symmetry of the cubed-sphere on the edge
      if ((fvm%faceno==1) .OR. (fvm%faceno==6)) then
        fvm%acartx1(jx_min1:jx_max1)=fvm%acartx(jx_min1:jx_max1)
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=-fvm%acarty(1+(jy_max1-jy))   
        end do 
      end if
      if (fvm%faceno==2) then
        ! flip orientation
        swap1=.TRUE.
        invx1=-1
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=-fvm%acartx(jx_max1-(jx-jx_min1))
        end do
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=-fvm%acarty(1+(jy_max1-jy))   
        end do 
      end if
      ! use symmetry of the cubed-sphere on the edge
      if (fvm%faceno==3) then
        invx1=-1
        invy1=-1
        do jx=jx_min1,jx_max1
          jtmp=jx_max1+jx_min1-jx
          fvm%acartx1(jtmp)=-fvm%acartx(jx)
        end do
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(1+(jy-jy_min1))  
        end do 
      end if
      if (fvm%faceno==4) then
        ! flip orientation
        swap1=.TRUE.
        invy1=-1
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(1+(jy-jy_min1))   
        end do 
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=fvm%acartx(jx)
        end do
      end if
      ! use symmetry of the cubed-sphere on the edge
      if (fvm%faceno==5) then
        invx1=-1
        invy1=-1
        do jx=jx_min1,jx_max1
          jtmp=jx_max1+jx_min1-jx
          fvm%acartx1(jtmp)=-fvm%acartx(jx)
        end do
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(1+(jy-jy_min1))  
        end do 
      end if          
      
  !SOUTH AND WEST EDGE ON THE CUBE SIDE        
    case(swest)
      jx_min=1
      jx_max=nc+1+nhe
      jy_min=1
      jy_max=nc+1+nhe 
      do i=0,nhe-1
        fvm%acartx(nc+2+i) = tan(elem%corners(1)%x+(nc+1+i)*fvm%dalpha) 
        fvm%acarty(nc+2+i) = tan(elem%corners(1)%y+(nc+1+i)*fvm%dbeta)
      enddo  
      jx_min1=1
      jx_max1=nc+1+nhe
      jy_min1=1-nhe
      jy_max1=1
      ! use symmetry of the cubed-sphere on the edge
      if ((fvm%faceno==1) .OR. (fvm%faceno==6)) then
        fvm%acartx1(jx_min1:jx_max1)=fvm%acartx(jx_min1:jx_max1)
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=-fvm%acarty(1+(jy_max1-jy))   
        end do 
      end if
      if (fvm%faceno==2) then
        ! flip orientation
        swap1=.TRUE.
        invx1=-1
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=-fvm%acartx(jx_max1+1-jx)
        end do
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=-fvm%acarty(1+(jy_max1-jy))  
        end do 
      end if
      ! use symmetry of the cubed-sphere on the edge
      if (fvm%faceno==3) then
        invx1=-1
        invy1=-1
        do jx=jx_min1,jx_max1
          jtmp=jx_max1+jx_min1-jx
          fvm%acartx1(jtmp)=-fvm%acartx(jx)
        end do
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(1+(jy-jy_min1))  
        end do 
      end if
      if (fvm%faceno==4) then
        ! flip orientation
        swap1=.TRUE.
        invy1=-1
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(1+(jy-jy_min1))   
        end do 
        fvm%acartx1(jx_min1:jx_max1)=fvm%acartx(jx_min1:jx_max1)
      end if
      if (fvm%faceno==5) then
        invx1=-1
        invy1=-1
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=-fvm%acartx(jx_max1+1-jx)   
        end do 
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(1+(jy-jy_min1))
        end do
      end if 
       
      jx_min2=1-nhe
      jx_max2=1
      jy_min2=1
      jy_max2=nc+1+nhe
      if (fvm%faceno<=4) then
        fvm%acarty2(jy_min2:jy_max2)=fvm%acarty(jy_min2:jy_max2)
        do jx=jx_min2,jx_max2
          fvm%acartx2(jx)=-fvm%acartx(1+(jx_max2-jx))   
        end do 
      end if
      if (fvm%faceno==5) then
        ! flip orientation
        swap2=.TRUE.
        invx2=-1
        fvm%acarty2(jy_min2:jy_max2)=fvm%acarty(jy_min2:jy_max2)
        do jx=jx_min2,jx_max2
          fvm%acartx2(jx)=fvm%acartx(1+jx-jx_min2)   
        end do 
      end if
      ! use symmetry of the cubed-sphere on the edge
      if (fvm%faceno==6) then
        ! flip orientation
        swap2=.TRUE.
        invy2=-1
        do jy=jy_min2,jy_max2
          jtmp=jy_max2+jy_min2-jy
          fvm%acarty2(jtmp)=-fvm%acarty(jy)
        end do
        do jx=jx_min2,jx_max2
          fvm%acartx2(jx)=-fvm%acartx(1+jx_max2-jx)  
        end do 
      end if
      
  !SOUTH AND EAST EDGE ON THE CUBE SIDE
    case(seast)
      jx_min=1-nhe
      jx_max=nc+1
      jy_min=1
      jy_max=nc+1+nhe 
      do i=0,nhe-1
        fvm%acartx(-i) = tan(elem%corners(1)%x-(i+1)*fvm%dalpha) 
        fvm%acarty(nc+2+i) = tan(elem%corners(1)%y+(nc+1+i)*fvm%dbeta)
      enddo
      jx_min1=1-nhe
      jx_max1=nc+1
      jy_min1=1-nhe
      jy_max1=1
      ! use symmetry of the cubed-sphere on the edge
      if ((fvm%faceno==1) .OR. (fvm%faceno==6)) then
        fvm%acartx1(jx_min1:jx_max1)=fvm%acartx(jx_min1:jx_max1)
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=-fvm%acarty(1+(jy_max1-jy))   
        end do 
      end if
      if (fvm%faceno==2) then
        ! flip orientation
        swap1=.TRUE.
        invx1=-1
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=-fvm%acartx(nc+1-(jx-jx_min1))
        end do
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=-fvm%acarty(1+(jy_max1-jy))   
        end do 
      end if
      ! use symmetry of the cubed-sphere on the edge
      if (fvm%faceno==3) then
        invx1=-1
        invy1=-1
        do jx=jx_min1,jx_max1
          jtmp=jx_max1+jx_min1-jx
          fvm%acartx1(jtmp)=-fvm%acartx(jx)
        end do
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(1+(jy-jy_min1))  
        end do 
      end if
      if (fvm%faceno==4) then
        ! flip orientation
        swap1=.TRUE.
        invy1=-1
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(1+(jy-jy_min1))   
        end do 
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=fvm%acartx(jx)
        end do
      end if
      ! use symmetry of the cubed-sphere on the edge
      if (fvm%faceno==5) then
        invx1=-1
        invy1=-1
        do jx=jx_min1,jx_max1
          jtmp=jx_max1+jx_min1-jx
          fvm%acartx1(jtmp)=-fvm%acartx(jx)
        end do
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(1+(jy-jy_min1))  
        end do 
      end if

      jx_min2=nc+1
      jx_max2=nc+1+nhe
      jy_min2=1
      jy_max2=nc+1+nhe
      if (fvm%faceno<=4) then
        fvm%acarty2(jy_min2:jy_max2)=fvm%acarty(jy_min2:jy_max2)
        do jx=jx_min2,jx_max2
          fvm%acartx2(jx)=-fvm%acartx(nc+1-(jx-jx_min2))   
        end do 
      end if
      if (fvm%faceno==5) then
        ! flip orientation
        swap2=.TRUE.
        invy2=-1
        do jy=jy_min2,jy_max2
          jtmp=jy_max2+jy_min2-jy
          fvm%acarty2(jtmp)=-fvm%acarty(jy)
        end do
        do jx=jx_min2,jx_max2
          fvm%acartx2(jx)=-fvm%acartx(nc+1-(jx-jx_min2))   
        end do 
      end if
      ! use symmetry of the cubed-sphere on the edge
      if (fvm%faceno==6) then
        ! flip orientation
        swap2=.TRUE.
        invx2=-1
        fvm%acarty2(jy_min2:jy_max2)=fvm%acarty(jy_min2:jy_max2)
        do jx=jx_min2,jx_max2
          fvm%acartx2(jx)=fvm%acartx(nc+1-(jx_max2-jx))  
        end do 
      end if

  !NORTH AND EAST EDGE ON THE CUBE SIDE    
    case(neast) 
      jx_min=1-nhe
      jx_max=nc+1
      jy_min=1-nhe
      jy_max=nc+1  
      do i=0,nhe-1
        fvm%acartx(-i) = tan(elem%corners(1)%x-(i+1)*fvm%dalpha) 
        fvm%acarty(-i) = tan(elem%corners(1)%y-(i+1)*fvm%dbeta)
      enddo
      jx_min1=1-nhe
      jx_max1=nc+1
      jy_min1=nc+1
      jy_max1=nc+1+nhe
      ! use symmetry of the cubed-sphere on the edge
      if ((fvm%faceno==1) .OR. (fvm%faceno==5)) then
        fvm%acartx1(jx_min1:jx_max1)=fvm%acartx(jx_min1:jx_max1)
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=-fvm%acarty(nc+1-(jy-jy_min1))   
        end do 
      end if
      if (fvm%faceno==2) then
        ! flip orientation
        swap1=.TRUE.
        invy1=-1
        fvm%acartx1(jx_min1:jx_max1)=fvm%acartx(jx_min1:jx_max1)
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(nc+1-(jy_max1-jy))   
        end do 
      end if
      if (fvm%faceno==4) then
        ! flip orientation
        swap1=.TRUE.
        invx1=-1
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=-fvm%acarty(nc+1-(jy-jy_min1))   
        end do 
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=-fvm%acartx(jx_max-(jx-jx_min1))
        end do
      end if
      ! use symmetry of the cubed-sphere on the edge
      if ((fvm%faceno==3) .OR. (fvm%faceno==6)) then
        invx1=-1
        invy1=-1
        do jx=jx_min1,jx_max1
          jtmp=jx_max1+jx_min1-jx
          fvm%acartx1(jtmp)=-fvm%acartx(jx)
        end do
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(nc+1-(jy_max1-jy)) 
        end do 
      end if

      jx_min2=nc+1
      jx_max2=nc+1+nhe
      jy_min2=1-nhe
      jy_max2=nc+1
      if (fvm%faceno<=4) then
        fvm%acarty2(jy_min2:jy_max2)=fvm%acarty(jy_min2:jy_max2)
        do jx=jx_min2,jx_max2
          fvm%acartx2(jx)=-fvm%acartx(nc+1-(jx-jx_min2))   
        end do 
      end if
      if (fvm%faceno==5) then
        ! flip orientation
        swap2=.TRUE.
        invy2=-1
        do jy=jy_min2,jy_max2
          jtmp=jy_max2+jy_min2-jy
          fvm%acarty2(jtmp)=-fvm%acarty(jy)
        end do
        do jx=jx_min2,jx_max2
          fvm%acartx2(jx)=-fvm%acartx(nc+1-(jx-jx_min2))   
        end do 
      end if
      ! use symmetry of the cubed-sphere on the edge
      if (fvm%faceno==6) then
        ! flip orientation
        swap2=.TRUE.
        invx2=-1
        fvm%acarty2(jy_min2:jy_max2)=fvm%acarty(jy_min2:jy_max2)
        do jx=jx_min2,jx_max2
          fvm%acartx2(jx)=fvm%acartx(nc+1-(jx_max2-jx))  
        end do 
      end if
      
  !NORTH AND WEST EDGE ON THE CUBE SIDE      
    case(nwest)
      jx_min=1
      jx_max=nc+1+nhe
      jy_min=1-nhe
      jy_max=nc+1
      do i=0,nhe-1
        fvm%acartx(nc+2+i) = tan(elem%corners(1)%x+(nc+1+i)*fvm%dalpha) 
        fvm%acarty(-i) = tan(elem%corners(1)%y-(i+1)*fvm%dbeta)
      enddo
      jx_min1=1
      jx_max1=nc+1+nhe
      jy_min1=nc+1
      jy_max1=nc+1+nhe
      ! use symmetry of the cubed-sphere on the edge
      if ((fvm%faceno==1) .OR. (fvm%faceno==5)) then
        fvm%acartx1(jx_min1:jx_max1)=fvm%acartx(jx_min1:jx_max1)
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=-fvm%acarty(nc+1-(jy-jy_min1))   
        end do 
      end if
      ! use symmetry of the cubed-sphere on the edge
      if ((fvm%faceno==3) .OR. (fvm%faceno==6)) then
        invx1=-1
        invy1=-1
        do jx=jx_min1,jx_max1
          jtmp=jx_max1+jx_min1-jx
          fvm%acartx1(jtmp)=-fvm%acartx(jx)
        end do
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(nc+1-(jy_max1-jy))  
        end do 
      end if
      if (fvm%faceno==2) then
        ! flip orientation
        swap1=.TRUE.
        invy1=-1
        fvm%acartx1(jx_min1:jx_max1)=fvm%acartx(jx_min1:jx_max1)
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=fvm%acarty(nc+1-(jy_max1-jy))   
        end do 
      end if
      if (fvm%faceno==4) then
        ! flip orientation
        swap1=.TRUE.
        invx1=-1
        do jy=jy_min1,jy_max1
          fvm%acarty1(jy)=-fvm%acarty(nc+1-(jy-jy_min1))   
        end do 
        do jx=jx_min1,jx_max1
          fvm%acartx1(jx)=-fvm%acartx(jx_max1-(jx-jx_min1))
        end do
      end if
      
      jx_min2=1-nhe
      jx_max2=1
      jy_min2=1-nhe
      jy_max2=nc+1
      if (fvm%faceno<=4) then
        fvm%acarty2(jy_min2:jy_max2)=fvm%acarty(jy_min2:jy_max2)
        do jx=jx_min2,jx_max2
          fvm%acartx2(jx)=-fvm%acartx(1+(jx_max2-jx))   
        end do 
      end if
      if (fvm%faceno==5) then
        ! flip orientation
        swap2=.TRUE.
        invx2=-1
        fvm%acarty2(jy_min2:jy_max2)=fvm%acarty(jy_min2:jy_max2)
        do jx=jx_min2,jx_max2
          fvm%acartx2(jx)=fvm%acartx(1+jx-jx_min2)   
        end do 
      end if
      ! use symmetry of the cubed-sphere on the edge
      if (fvm%faceno==6) then
        ! flip orientation
        swap2=.TRUE.
        invy2=-1
        do jy=jy_min2,jy_max2
          jtmp=jy_max2+jy_min2-jy
          fvm%acarty2(jtmp)=-fvm%acarty(jy)
        end do
        do jx=jx_min2,jx_max2
          fvm%acartx2(jx)=-fvm%acartx(1+jx_max2-jx)  
        end do 
      end if      
      
  !THIS CASE SHOULD NOT HAPPEN, SINCE WE ASSUME NE>=2!     
    case default
      print *, 'Fatal Error in fvm_line_integrals_mod.F90.'
      call abortmp('Selected case for cubeboundary does not exists!')     
  end select
  
  fvm%acartx(jx_min-1)=-bignum 
  fvm%acartx(jx_max+1)=bignum
  fvm%acarty(jy_min-1)=-bignum 
  fvm%acarty(jy_max+1)=bignum
  fvm%jx_min=jx_min; fvm%jx_max=jx_max; 
  fvm%jy_min=jy_min; fvm%jy_max=jy_max;
  
  if(fvm%cubeboundary>0) then   !only for elements on the cube edge
    fvm%acartx1(jx_min1-1)=-bignum; fvm%acartx1(jx_max1+1)=bignum;
    fvm%acarty1(jy_min1-1)=-bignum; fvm%acarty1(jy_max1+1)=bignum;    
    fvm%jx_min1=jx_min1; fvm%jx_max1=jx_max1; 
    fvm%jy_min1=jy_min1; fvm%jy_max1=jy_max1;
    fvm%swap1=swap1;
    fvm%invx1=invx1;fvm%invy1=invy1
    if(fvm%cubeboundary>4) then !only for elements on the corner
      fvm%acartx2(jx_min2-1)=-bignum; fvm%acartx2(jx_max2+1)=bignum;
      fvm%acarty2(jy_min2-1)=-bignum; fvm%acarty2(jy_max2+1)=bignum;
      fvm%jx_min2=jx_min2; fvm%jx_max2=jx_max2;
      fvm%jy_min2=jy_min2; fvm%jy_max2=jy_max2;
      fvm%swap2=swap2
      fvm%invx2=invx2;fvm%invy2=invy2
    endif
  endif
end subroutine create_ari
!END SUBROUTINE fvm_MESH_ARI----------------------------------------------CE-for FVM!

! ----------------------------------------------------------------------------------!
!SUBROUTINE CREATE_INTERPOLATIION_POINTS----------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: for elements, which share a cube edge, we have to do some            !
!        interpolation on different cubic faces, also in the halo region:           !
!        because we also need the reconstruction coefficients in the halo zone,     !
!        which is basically calculated twice, on the original cell of an element    !
!        on face A and on a cell in the halo region of an element of face B         !
!        The crux is, that the interpolation has to be the same to ensure           !
!        conservation of the scheme                                                 !
!        SYMMETRY of the CUBE is used for calucaltion the interpolation_point       !
!                                                                                   !
! CALLS: interpolation_point                                                        !
! INPUT/OUTPUT:                                                                     !
!        elem     ...  element structure from HOMME                                 ! 
!        fvm   ...  structure                                                       !                 
!-----------------------------------------------------------------------------------!
subroutine create_interpolation_points(elem,fvm)
  implicit none
  type (element_t), intent(in)      :: elem
  type (fvm_struct), intent(inout)   :: fvm  
  
  real    (kind=real_kind), dimension(1-nhc:nc+nhc)  :: gnomxstart, gnomxend, &
                                                        gnomystart, gnomyend
  integer                                       :: i, halo, ida, ide, iref1, iref2
  type (cartesian2D_t)                          :: tmpgnom     

! Jose Garcia
! creating an empty routine when MESH is defined.
! We allow this so both MESH and fvm can be compiled together
! because some executions paths that do not use this routine are 
! important.

#ifndef MESH  
  ! element is not on a corner, but shares a cube edge (call of subroutine)
  if(fvm%cubeboundary <= 4) then
    gnomxstart(1-nhc)=elem%corners(1)%x-(nhc-0.5)*fvm%dalpha
    gnomystart(1-nhc)=elem%corners(1)%y-(nhc-0.5)*fvm%dbeta
    do i=2-nhc,nc+nhc
      gnomxstart(i)=gnomxstart(i-1)+fvm%dalpha
      gnomystart(i)=gnomystart(i-1)+fvm%dbeta
    end do
    select case (fvm%cubeboundary)
      !INTERIOR element
      case(0)
        ! nothing to do!
      !CASE WEST
      case(west) 
        ! west zone
        ida=-3
        ide=nc+4
        do halo=1,2
          iref1=ida
          tmpgnom%x=cube_xstart-(halo-0.5)*fvm%dalpha
          do i=-1,nc+2
            tmpgnom%y=gnomystart(i)
            call interpolation_point(tmpgnom,gnomystart,1,4,1,0, fvm%interp(i,halo,1),&
                                     ida,ide,iref1,fvm%ibase(i,halo,1))
          end do
        end do
        !westhalo: east zone (because of the reconstruction in the halo, conservation!)
        do halo=1,2
          iref1=ida
          tmpgnom%x=cube_xend+(halo-0.5)*fvm%dalpha
          do i=-1,nc+2
            tmpgnom%y=gnomystart(i)
            call interpolation_point(tmpgnom,gnomystart,4,1,1,0, fvm%interphalo(i,halo,1),&
                                     ida,ide,iref1,fvm%ibasehalo(i,halo,1))                                    
          end do
        end do
      !CASE EAST  
      case(east)
        ! east zone
        ida=-3
        ide=nc+4
        do halo=1,2
          iref1=ida
          tmpgnom%x=cube_xend+(halo-0.5)*fvm%dalpha
          do i=-1,nc+2
            tmpgnom%y=gnomystart(i)
            call interpolation_point(tmpgnom,gnomystart,1,2,1,0, fvm%interp(i,halo,1),&
                                     ida,ide,iref1,fvm%ibase(i,halo,1))                                                  
          end do 
        end do  
        !easthalo: west zone (because of the reconstruction in the halo, conservation!)
        do halo=1,2
          iref1=ida
          tmpgnom%x=cube_xstart-(halo-0.5)*fvm%dalpha
          do i=-1,nc+2
            tmpgnom%y=gnomystart(i)
            call interpolation_point(tmpgnom,gnomystart,1,4,1,0, fvm%interphalo(i,halo,1),&
                                     ida,ide,iref1,fvm%ibasehalo(i,halo,1))                                    
          end do
        end do
      !CASE NORTH    
      case(north)
        ! north zone
          ida=-3
          ide=nc+4
        do halo=1,2 
          tmpgnom%y=cube_yend+(halo-0.5)*fvm%dbeta
          iref1=ida 
          do i=-1,nc+2
            tmpgnom%x=gnomxstart(i)
            call interpolation_point(tmpgnom,gnomxstart,1,6,0,0, fvm%interp(i,halo,2),&
                                     ida,ide,iref1,fvm%ibase(i,halo,2))                                      
          end do
        end do
        !northhalo: south zone (because of the reconstruction in the halo, conservation!)
        do halo=1,2
          iref1=ida
          tmpgnom%y=cube_ystart-(halo-0.5)*fvm%dbeta
          do i=-1,nc+2
            tmpgnom%x=gnomxstart(i)
            call interpolation_point(tmpgnom,gnomxstart,1,5,0,0, fvm%interphalo(i,halo,2),&
                                     ida,ide,iref1,fvm%ibasehalo(i,halo,2))                                    
          end do
        end do            
      !CASE SOUTH   
      case(south)
       !south zone
       ida=-3
       ide=nc+4
       do halo=1,2
          iref1=ida
          tmpgnom%y=cube_ystart-(halo-0.5)*fvm%dbeta
          do i=-1,nc+2
            tmpgnom%x=gnomxstart(i)
            call interpolation_point(tmpgnom,gnomxstart,1,5,0,0, fvm%interp(i,halo,2),&
                                     ida,ide,iref1,fvm%ibase(i,halo,2))                                      
          end do
        end do      

        !southhalo: north zone (because of the reconstruction in the halo, conservation!) 
        do halo=1,2
          iref1=ida
          tmpgnom%y=cube_yend+(halo-0.5)*fvm%dbeta
          do i=-1,nc+2
            tmpgnom%x=gnomxstart(i)
            call interpolation_point(tmpgnom,gnomxstart,1,6,0,0, fvm%interphalo(i,halo,2),&
                                     ida,ide,iref1,fvm%ibasehalo(i,halo,2))                                    
          end do
        end do
      !THIS CASE SHOULD NOT HAPPEN!     
        case default
           print *,'Fatal Error in first select statement:'
           call abortmp('fvm_reconstruction_mod.F90 subroutine fillhalo_cubic!' )
    end select
  !CORNER TREATMENT
  else  
    gnomxstart(1-nhc)=cube_xstart-(nhc-0.5)*fvm%dalpha
    gnomxend(nc+nhc)=cube_xend+(nhc-0.5)*fvm%dalpha
    gnomystart(1-nhc)=cube_ystart-(nhc-0.5)*fvm%dbeta
    gnomyend(nc+nhc)=cube_yend+(nhc-0.5)*fvm%dbeta
    do i=2-nhc,nc+nhc
      gnomxstart(i)=gnomxstart(i-1)+fvm%dalpha
      gnomxend(nc+1-i)=gnomxend(nc+2-i)-fvm%dalpha
      gnomystart(i)=gnomystart(i-1)+fvm%dbeta
      gnomyend(nc+1-i)=gnomyend(nc+2-i)-fvm%dbeta
    end do
   
    select case (fvm%cubeboundary)
      !CASE SOUTH WEST
      case(swest) 
        ! west zone
        do halo=1,2
          tmpgnom%x=cube_xstart-(halo-0.5)*fvm%dalpha
          ida=1
          ide=nc+4
          iref1=ida
          do i=0,nc+2
            tmpgnom%y=gnomystart(i)
            call interpolation_point(tmpgnom,gnomystart,1,4,1,-1, fvm%interp(i,halo,1),&
                                     ida,ide,iref1,fvm%ibase(i,halo,1))                                      
          end do
        ! south zone
          tmpgnom%y=cube_ystart-(halo-0.5)*fvm%dbeta
          iref2=ida
          do i=0,nc+2
            tmpgnom%x=gnomxstart(i)
            call interpolation_point(tmpgnom,gnomxstart,1,5,0,-1, fvm%interp(i,halo,2),&
                                     ida,ide,iref2,fvm%ibase(i,halo,2))                                      
          end do
        end do
        !westhalo: east zone & south zone
        !(because of the reconstruction in the halo, conservation!) 
        do halo=1,2
          ida=1
          ide=nc+4
          tmpgnom%x=cube_xend+(halo-0.5)*fvm%dalpha
          iref1=ida
          do i=0,nc+2
            tmpgnom%y=gnomystart(i)
            call interpolation_point(tmpgnom,gnomystart,4,1,1,-1, fvm%interphalo(i,halo,1),&
                                     ida,ide,iref1,fvm%ibasehalo(i,halo,1))                                    
          end do
          tmpgnom%y=cube_ystart-(halo-0.5)*fvm%dbeta

          ida=nc-nhc+1
          ide=nc 
          iref2=ida
          do i=0,nhe+1
            tmpgnom%x=gnomxend(nc-nhe+i)
            call interpolation_point(tmpgnom,gnomxend,1,5,0,1, fvm%interphaloex(i,halo,1),&
                                     ida,ide,iref2,fvm%ibasehaloex(i,halo,1))                                          
          end do
        end do
        !southhalo: north zone & west zone
        !(because of the reconstruction in the halo, conservation!) 
        do halo=1,2
          ida=1
          ide=nc+4
          tmpgnom%y=cube_yend+(halo-0.5)*fvm%dbeta
          iref1=ida
          do i=0,nc+2
            tmpgnom%x=gnomxstart(i)
            call interpolation_point(tmpgnom,gnomxstart,1,6,0,-1, fvm%interphalo(i,halo,2),&
                                     ida,ide,iref1,fvm%ibasehalo(i,halo,2))                                    
          end do
          tmpgnom%x=cube_xstart-(halo-0.5)*fvm%dalpha
          ida=nc-nhc+1
          ide=nc
          iref2=ida
          do i=0,nhe+1
            tmpgnom%y=gnomyend(nc-nhe+i)
            call interpolation_point(tmpgnom,gnomyend,1,4,1,1, fvm%interphaloex(i,halo,2),&
                                     ida,ide,iref2,fvm%ibasehaloex(i,halo,2))                                          
          end do
        end do
      !CASE SOUTH EAST  
      case(seast)
        ! east zone
        do halo=1,2
          tmpgnom%x=cube_xend+(halo-0.5)*fvm%dalpha
          ida=1
          ide=nc+4
          iref1=ida
          do i=0,nc+2
            tmpgnom%y=gnomystart(i)
            call interpolation_point(tmpgnom,gnomystart,1,2,1,-1, fvm%interp(i,halo,1),&
                                     ida,ide,iref1,fvm%ibase(i,halo,1))                                      
          end do
        ! south zone
          tmpgnom%y=cube_ystart-(halo-0.5)*fvm%dbeta
          ida=-3
          ide=nc
          iref2=ida
          do i=-1,nc+1
            tmpgnom%x=gnomxend(i)
            call interpolation_point(tmpgnom,gnomxend,1,5,0,1, fvm%interp(i,halo,2),&
                                     ida,ide,iref2,fvm%ibase(i,halo,2))                                       
          end do
        end do      
        !easthalo: west zone & south zone
        !(because of the reconstruction in the halo, conservation!) 
        do halo=1,2
          tmpgnom%x=cube_xstart-(halo-0.5)*fvm%dalpha
          ida=1
          ide=nc+4
          iref1=ida
          do i=0,nc+2
            tmpgnom%y=gnomystart(i)
            call interpolation_point(tmpgnom,gnomystart,1,4,1,-1, fvm%interphalo(i,halo,1),&
                                     ida,ide,iref1,fvm%ibasehalo(i,halo,1))                                    
          end do
          tmpgnom%y=cube_ystart-(halo-0.5)*fvm%dbeta
          ida=1
          ide=nhc 
          iref2=ida
          do i=0,nhe+1
            tmpgnom%x=gnomxstart(i)
            call interpolation_point(tmpgnom,gnomystart,1,5,0,-1, fvm%interphaloex(i,halo,1),&
                                     ida,ide,iref2,fvm%ibasehaloex(i,halo,1))                                      
          end do
        end do
        !(because of the reconstruction in the halo, conservation!) 
        !southhalo: north zone & west zone
        do halo=1,2
          tmpgnom%y=cube_yend+(halo-0.5)*fvm%dbeta
          ida=-3
          ide=nc
          iref1=ida
          do i=-1,nc+1
            tmpgnom%x=gnomxend(i)
            call interpolation_point(tmpgnom,gnomxend,1,6,0,1, fvm%interphalo(i,halo,2),&
                                     ida,ide,iref1,fvm%ibasehalo(i,halo,2))                                    
          end do
          tmpgnom%x=cube_xend+(halo-0.5)*fvm%dalpha
          ida=nc-nhc+1
          ide=nc 
          iref2=ida
          do i=0,nhe+1
            tmpgnom%y=gnomyend(nc-nhe+i)
            call interpolation_point(tmpgnom,gnomyend,1,2,1,1, fvm%interphaloex(i,halo,2),&
                                     ida,ide,iref2,fvm%ibasehaloex(i,halo,2))                                          
          end do
        end do      
      !CASE NORTH EAST     
      case(neast)
        ! east zone
        do halo=1,2
          tmpgnom%x=cube_xend+(halo-0.5)*fvm%dalpha
          ida=-3
          ide=nc
          iref1=ida
          do i=-1,nc+1 
            tmpgnom%y=gnomyend(i)
            call interpolation_point(tmpgnom,gnomyend,1,2,1,1, fvm%interp(i,halo,1),&
                                     ida,ide,iref1,fvm%ibase(i,halo,1))                                      
          end do
        ! north zone
          tmpgnom%y=cube_yend+(halo-0.5)*fvm%dbeta
          ida=-3
          ide=nc
          iref2=ida
          do i=-1,nc+1
            tmpgnom%x=gnomxend(i)
            call interpolation_point(tmpgnom,gnomxend,1,6,0,1, fvm%interp(i,halo,2),&
                                     ida,ide,iref2,fvm%ibase(i,halo,2))                                      
          end do
        end do
        !easthalo: west zone & north zone
        !(because of the reconstruction in the halo, conservation!) 
        do halo=1,2
          ida=-3
          ide=nc
          iref1=ida
          tmpgnom%x=cube_xstart-(halo-0.5)*fvm%dalpha
          do i=-1,nc+1
            tmpgnom%y=gnomyend(i)
            call interpolation_point(tmpgnom,gnomyend,1,4,1,1, fvm%interphalo(i,halo,1),&
                                     ida,ide,iref1,fvm%ibasehalo(i,halo,1))                                    
          end do
          tmpgnom%y=cube_yend+(halo-0.5)*fvm%dbeta
          ida=1
          ide=nhc
          iref2=ida
          do i=0,nhe+1
            tmpgnom%x=gnomxstart(i)
            call interpolation_point(tmpgnom,gnomxstart,1,6,0,-1, fvm%interphaloex(i,halo,1),&
                                     ida,ide,iref2,fvm%ibasehaloex(i,halo,1))                                      
          end do
        end do
        !northhalo: south zone & east zone
        !(because of the reconstruction in the halo, conservation!)         
        do halo=1,2
          ida=-3
          ide=nc
          iref1=ida
          tmpgnom%y=cube_ystart-(halo-0.5)*fvm%dbeta
          do i=-1,nc+1
            tmpgnom%x=gnomxend(i)
            call interpolation_point(tmpgnom,gnomyend,1,5,0,1, fvm%interphalo(i,halo,2),&
                                     ida,ide,iref1,fvm%ibasehalo(i,halo,2))                                    
          end do
          tmpgnom%x=cube_xend+(halo-0.5)*fvm%dalpha
          ida=1
          ide=nhc
          iref2=ida
          do i=0,nhe+1
            tmpgnom%y=gnomystart(i)
            call interpolation_point(tmpgnom,gnomystart,1,2,1,-1, fvm%interphaloex(i,halo,2),&
                                     ida,ide,iref2,fvm%ibasehaloex(i,halo,2))                                      
          end do
        end do            
      !CASE NORTH WEST   
      case(nwest)
        ! west zone
        do halo=1,2
          tmpgnom%x=cube_xstart-(halo-0.5)*fvm%dalpha
          ida=-3
          ide=nc
          iref1=ida
          do i=-1,nc+1
            tmpgnom%y=gnomyend(i)
            call interpolation_point(tmpgnom,gnomyend,1,4,1,1, fvm%interp(i,halo,1),&
                                     ida,ide,iref1,fvm%ibase(i,halo,1))                                      
          end do
        ! north zone
          tmpgnom%y=cube_yend+(halo-0.5)*fvm%dbeta
          ida=1
          ide=nc+4
          iref2=ida
          do i=0,nc+2
            tmpgnom%x=gnomxstart(i)
            call interpolation_point(tmpgnom,gnomxstart,1,6,0,-1, fvm%interp(i,halo,2),&
                                     ida,ide,iref2,fvm%ibase(i,halo,2))                                      
          end do
        end do
        !westhalo: east zone & north zone
        !(because of the reconstruction in the halo, conservation!) 
        do halo=1,2
          ida=-3
          ide=nc
          iref1=ida
          tmpgnom%x=cube_xend+(halo-0.5)*fvm%dalpha
          do i=-1,nc+1
            tmpgnom%y=gnomyend(i)
            call interpolation_point(tmpgnom,gnomyend,4,1,1,1, fvm%interphalo(i,halo,1),&
                                     ida,ide,iref1,fvm%ibasehalo(i,halo,1))                                    
          end do
          tmpgnom%y=cube_yend+(halo-0.5)*fvm%dbeta
          ida=nc-nhc+1
          ide=nc
          iref2=ida
          do i=0,nhe+1
            tmpgnom%x=gnomxend(nc-nhe+i)
            call interpolation_point(tmpgnom,gnomxend,1,6,0,1, fvm%interphaloex(i,halo,1),&
                                     ida,ide,iref2,fvm%ibasehaloex(i,halo,1))                                          
          end do
        end do        
        !northhalo: south zone & west zone
        !(because of the reconstruction in the halo, conservation!) 
        do halo=1,2
          ida=1
          ide=nc+4
          iref1=ida
          tmpgnom%y=cube_ystart-(halo-0.5)*fvm%dbeta
          do i=0,nc+2
            tmpgnom%x=gnomxstart(i)
            call interpolation_point(tmpgnom,gnomxstart,1,5,0,-1, fvm%interphalo(i,halo,2),&
                                     ida,ide,iref1,fvm%ibasehalo(i,halo,2))                                    
          end do
          tmpgnom%x=cube_xstart-(halo-0.5)*fvm%dalpha
          ida=1
          ide=nhc
          iref2=ida
          do i=0,nhe+1
            tmpgnom%y=gnomystart(i)
            call interpolation_point(tmpgnom,gnomystart,1,4,1,-1, fvm%interphaloex(i,halo,2),&
                                     ida,ide,iref2,fvm%ibasehaloex(i,halo,2))                                      
          end do
        end do
        !THIS CASE SHOULD NOT HAPPEN!     
          case default
            print *,'Fatal Error in second select statement:'
            call abortmp('fvm_reconstruction_mod.F90 subroutine create_interpolationpoint!')
         end select
  endif

#endif 
!  ^
!  |
! endif for ifndef MESH

end subroutine create_interpolation_points
!END SUBROUTINE CREATE_INTERPOLATION_POINTS-------------------------------CE-for FVM!



! ----------------------------------------------------------------------------------!
!SUBROUTINE INTERPOLATION_POINT-------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 14.November 2011                                         !
! DESCRIPTION: calculates the interpolation point on from face 1 in face 2 in       !
!        alpha/beta coordinates, only 1D                                            !
!                                                                                   !
! CALLS: cubedsphere2cart, cart2cubedsphere                                         !
! INPUT: gnom... 1D coordinates                                                     !
!        gnom1d... 1d coordinates                                                   !
!        face1... orginal face                                                      !
!        face2... target face (where the interpolation has to be done)              !
!        xy ... 0 for alpha coordinate, any other for beta                          !
!        except.which type, interior, left edge (-1), right edge (1)                !
!        point... interpolation point                                               !
!        ida  ... begin of interpval                                                !
!        ide  ... end of interpval                                                  !


! INPUT/OUTPUT/RETURN:                                                              !
!        iref ... where we start the search, is also an OUTPUT, so we know for the  !
!                 next point where to start                                         !
!-----------------------------------------------------------------------------------!
                                      !
! DESCRIPTION: searchs where the interpolation point has to be (iref), two values   !
!        of interpval on the left and on the right, except if we are out of range   !
!        which is indicated through ia and ie, respectively                         !
!        It is a 1D interpolation, use alpha/beta coordinates!!!                    !
!                                                                                   !
! CALLS: cubic_equispace_interp                                                     !
! INPUT: iref ... where we start the search, is also an OUTPUT, so we know for the  !
!                 next point where to start                                         !
!        ibaseref ... startindex of the four tracer value for the reconstruction    !
!        point    ... provides the difference of the interpolation point to use it  !
!                     directly in CUBIC_EQUISPACE_INTERP                            !                
!-----------------------------------------------------------------------------------!
subroutine interpolation_point(gnom,gnom1d,face1,face2,xy,except, point,ida,ide,iref,ibaseref)
  
  use coordinate_systems_mod, only : cubedsphere2cart, cart2cubedsphere, &
                                     cartesian2D_t,cartesian3D_t

  implicit none
  type (cartesian2D_t), intent(in)                                :: gnom  
  real (kind=real_kind), dimension(1-nhc:nc+nhc), intent(in)      :: gnom1d  
  integer, intent(in)                                             :: face1, face2, xy
  integer,intent(in)                                              :: except,ida, ide
  integer,intent(inout)                                           :: iref,ibaseref
  real (kind=real_kind), intent(inout)                            :: point
  
  
  type(cartesian3D_t)                 :: tmpcart3d
  type (cartesian2D_t)                :: tmpgnom
  
  tmpcart3d=cubedsphere2cart(gnom,face1)
  tmpgnom=cart2cubedsphere(tmpcart3d,face2)
  if(xy==0) then
    point=tmpgnom%x
  else
    point=tmpgnom%y
  end if 
  
  do while ((iref .ne. nc+3) .AND. (point>gnom1d(iref) ))
    iref = iref + 1
  end do
  if ((iref<=ida+1).AND.(except==-1)) then
    ibaseref=ida
  elseif ((iref>=ide).AND.(except==1)) then
    ibaseref=ide-3
  else
    ibaseref=iref-2  
  end if

  point=point-gnom1d(ibaseref)
  
end subroutine interpolation_point
!END SUBROUTINE INTERPOLATION_POINT---------------------------------------CE-for FVM!

end module fvm_control_volume_mod
