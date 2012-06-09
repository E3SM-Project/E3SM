#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!MODULE FVM_ANALYTIC_MOD--------------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! This module contains all analytical terms for fvm                                 ! 
!-----------------------------------------------------------------------------------!
module fvm_analytic_mod

  use kinds, only : real_kind, int_kind
  use dimensions_mod, only: nc, nhe, ntrac
  
  implicit none
  private
  public :: computexytosphere_moments
  public :: I_00, I_10, I_01, I_20, I_02, I_11
contains
! ----------------------------------------------------------------------------------!
!SUBROUTINE COMPUTEXYTOSPHERE_MOMENTS-------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 20.July 2011                                             !
! DESCRIPTION: Compute area and centroids/moments via line integrals                !
!                                                                                   !
! INPUT:  fvm ... fvm structure                                                     !
!         desc... edge descriptor, needed for reordering                            !
!-----------------------------------------------------------------------------------!
subroutine computexytosphere_moments(fvm,desc)
  use fvm_control_volume_mod, only:  fvm_struct 
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  use edge_mod, only : edgedescriptor_t
  implicit none
  
  type (fvm_struct), intent(inout)                  :: fvm
  type (edgedescriptor_t),intent(in)                  :: desc
  
  real (kind=real_kind)                               :: tmpval(5,1-nhe:nc+nhe)
  real (kind=real_kind)                               :: tmpvalarea(1-nhe:nc+nhe)
  
  integer                                             :: i,j,tmpi,tmpj
  
! compute area and centroid for the interior and halo zone of interior elements
  call moment_onsphere(fvm%area_sphere,fvm%spherecentroid,fvm%acartx,fvm%acarty,&
                       fvm%jx_min,fvm%jx_max,fvm%jy_min,fvm%jy_max)
  
  if(fvm%cubeboundary>0) then   !only for elements on the cube edge
    if(fvm%swap1) then
      call moment_onsphereswap(fvm%area_sphere,fvm%spherecentroid, &
                           fvm%acartx1,fvm%acarty1,fvm%jx_min1,fvm%jx_max1, &
                           fvm%jy_min1,fvm%jy_max1)
    else             
      call moment_onsphere(fvm%area_sphere,fvm%spherecentroid,fvm%acartx1,&
                fvm%acarty1,fvm%jx_min1,fvm%jx_max1,fvm%jy_min1,fvm%jy_max1)
    endif
    !I have to reorder it, note that the corner is already included (southwest in 
    ! south and west), note, no element can have a west and a east or a north and a
    ! south, we do not allow just one element per face
    if(desc%reverse(west) .or. desc%reverse(east)) then
      do i=fvm%jx_min1,fvm%jx_max1-1
        tmpval=fvm%spherecentroid(:,i,:)
        tmpvalarea=fvm%area_sphere(i,:)
        do j=fvm%jy_min1,fvm%jy_max1-1
          tmpj=fvm%jy_max1-1+fvm%jy_min1-j
          fvm%spherecentroid(:,i,tmpj)=tmpval(:,j)
          fvm%area_sphere(i,tmpj)=tmpvalarea(j)
        enddo
      enddo
    endif
    if(desc%reverse(north) .or. desc%reverse(south)) then
      do j=fvm%jy_min1,fvm%jy_max1-1
        tmpval=fvm%spherecentroid(:,:,j)
        tmpvalarea=fvm%area_sphere(:,j)
        do i=fvm%jx_min1,fvm%jx_max1-1
          tmpi=fvm%jx_max1-1+fvm%jx_min1-i
          fvm%spherecentroid(:,tmpi,j)=tmpval(:,i)
          fvm%area_sphere(tmpi,j)=tmpvalarea(i)
        enddo
      enddo
    endif
    
    if(fvm%cubeboundary>4) then !only for elements on the corner
      if(fvm%swap2) then
        call moment_onsphereswap(fvm%area_sphere,fvm%spherecentroid, &
                           fvm%acartx2,fvm%acarty2,fvm%jx_min2,fvm%jx_max2, &
                           fvm%jy_min2,fvm%jy_max2)
      else
        call moment_onsphere(fvm%area_sphere,fvm%spherecentroid,fvm%acartx2,&
               fvm%acarty2,fvm%jx_min2,fvm%jx_max2,fvm%jy_min2,fvm%jy_max2)
      endif

      if(desc%reverse(west) .or. desc%reverse(east)) then
        do i=fvm%jx_min2,fvm%jx_max2-1
          tmpval=fvm%spherecentroid(:,i,:)
          tmpvalarea=fvm%area_sphere(i,:)
          do j=fvm%jy_min2,fvm%jy_max2-1
            tmpj=fvm%jy_max2-1+fvm%jy_min2-j
            fvm%spherecentroid(:,i,tmpj)=tmpval(:,j)
            fvm%area_sphere(i,tmpj)=tmpvalarea(j)
          enddo
        enddo
      endif

      if(desc%reverse(north) .or. desc%reverse(south)) then
        do j=fvm%jy_min2,fvm%jy_max2-1
          tmpval=fvm%spherecentroid(:,:,j)
          tmpvalarea=fvm%area_sphere(:,j)
          do i=fvm%jx_min2,fvm%jx_max2-1
            tmpi=fvm%jx_max2-1+fvm%jx_min2-i
            fvm%spherecentroid(:,tmpi,j)=tmpval(:,i)
            fvm%area_sphere(tmpi,j)=tmpvalarea(i)
          enddo
        enddo
      endif
    endif
  endif
!   
end subroutine computexytosphere_moments
!END SUBROUTINE COMPUTEXYTOSPHERE_MOMENTS---------------------------------CE-for FVM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE MOMENT_ONSPHERE-----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 20.July 2011                                             !
! DESCRIPTION: Compute area and centroids/moments via line integrals                !
!                                                                                   !
! INPUT:  x  ...  x cartesian coordinats of the arrival grid on the cube            !
!         y  ...  y cartesian coordinats of the arrival grid on the cube            !
!         jx_min, jx_max, jy_min, jy_max                                            !
!            ... cell boundaries in x and y directions                              !
! INPUT/OUTPUT:                                                                     !
!         area      ... area of cells on the sphere                                 ! 
!         centroid  ... x,y,x^2,y^2,xy                                              !
!-----------------------------------------------------------------------------------!
subroutine moment_onsphere(area,centroid,x,y,jx_min,jx_max,jy_min,jy_max)
  use coordinate_systems_mod, only : surfareaxy
  implicit none
  real (kind=real_kind), dimension(1-nhe:nc+nhe,1-nhe:nc+nhe), intent(inout) :: area
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(inout) :: centroid
  real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in)   :: x, y
  integer, intent(in)                                  :: jx_min,jx_max,jy_min,jy_max
  integer                                              :: i,j

  do j=jy_min,jy_max-1
    do i=jx_min,jx_max-1
!       area(i,j) = surfareaxy(x(i),x(i+1),y(j),y(j+1)) 
      area(i,j) = (I_00(x(i+1),y(j+1)) - I_00(x(i),y(j+1)) + &
                         I_00(x(i),y(j)) - I_00(x(i+1),y(j)))
      ! Compute centroids via line integrals
      centroid(1,i,j) = (I_10(x(i+1),y(j+1)) - I_10(x(i),y(j+1)) + &
                         I_10(x(i),y(j)) - I_10(x(i+1),y(j))) / area(i,j)
      centroid(2,i,j) = (I_01(x(i+1),y(j+1)) - I_01(x(i),y(j+1)) + &
                         I_01(x(i),y(j)) - I_01(x(i+1),y(j))) / area(i,j)
      !ADD PHL START
     ! TAN(alpha)^2 component
      centroid(3,i,j) = (I_20(x(i+1),y(j+1)) - I_20(x(i),y(j+1)) + &
                         I_20(x(i),y(j)) - I_20(x(i+1),y(j))) / area(i,j)
     ! TAN(beta)^2 component
      centroid(4,i,j) = (I_02(x(i+1),y(j+1)) - I_02(x(i),y(j+1)) + &
                         I_02(x(i),y(j)) - I_02(x(i+1),y(j))) / area(i,j)
     ! TAN(alpha) TAN(beta) component
      centroid(5,i,j) = (I_11(x(i+1),y(j+1)) - I_11(x(i),y(j+1)) + &
                         I_11(x(i),y(j)) - I_11(x(i+1),y(j))) / area(i,j)
     !ADD PHL END
    end do
  end do  
end subroutine moment_onsphere
!END SUBROUTINE MOMENT_ONSPHERE-------------------------------------------CE-for FVM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE MOMENT_ONSPHERESWAP-------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 16. November 2011                                        !
! DESCRIPTION: Compute area and centroids/moments via line integrals for            !
!              swap version (if x and y coordinates has to be swapped)              !
!                                                                                   !
! INPUT:  x  ...  x cartesian coordinats of the arrival grid on the cube            !
!         y  ...  y cartesian coordinats of the arrival grid on the cube            !
!         jx_min, jx_max, jy_min, jy_max                                            !
!            ... cell boundaries in x and y directions                              !
! INPUT/OUTPUT:                                                                     !
!         area      ... area of cells on the sphere                                 ! 
!         centroid  ... x,y,x^2,y^2,xy                                              !
!-----------------------------------------------------------------------------------!
subroutine moment_onsphereswap(area,centroid,x,y,jx_min,jx_max,jy_min,jy_max)
  use coordinate_systems_mod, only : surfareaxy
  implicit none
  real (kind=real_kind), dimension(1-nhe:nc+nhe,1-nhe:nc+nhe), intent(inout) :: area
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(inout) :: centroid
  real (kind=real_kind), dimension(-nhe:nc+2+nhe), intent(in)   :: x, y
  integer, intent(in)                                  :: jx_min,jx_max,jy_min,jy_max
  integer                                              :: i,j

  do j=jy_min,jy_max-1
    do i=jx_min,jx_max-1
!       area(i,j) = surfareaxy(y(j),y(j+1),x(i),x(i+1)) 
      area(i,j) = (I_00(y(j+1),x(i+1)) - I_00(y(j+1),x(i)) + &
                         I_00(y(j),x(i)) - I_00(y(j),x(i+1)))
      ! Compute centroids via line integrals
      centroid(1,i,j) = (I_10(y(j+1),x(i+1)) - I_10(y(j+1),x(i)) + &
                         I_10(y(j),x(i)) - I_10(y(j),x(i+1))) / area(i,j)
      centroid(2,i,j) = (I_01(y(j+1),x(i+1)) - I_01(y(j+1),x(i)) + &
                         I_01(y(j),x(i)) - I_01(y(j),x(i+1))) / area(i,j)
      !ADD PHL START
     ! TAN(alpha)^2 component
      centroid(3,i,j) = (I_20(y(j+1),x(i+1)) - I_20(y(j+1),x(i)) + &
                         I_20(y(j),x(i)) - I_20(y(j),x(i+1))) / area(i,j)
     ! TAN(beta)^2 component
      centroid(4,i,j) = (I_02(y(j+1),x(i+1)) - I_02(y(j+1),x(i)) + &
                         I_02(y(j),x(i)) - I_02(y(j),x(i+1))) / area(i,j)
     ! TAN(alpha) TAN(beta) component
      centroid(5,i,j) = (I_11(y(j+1),x(i+1)) - I_11(y(j+1),x(i)) + &
                         I_11(y(j),x(i)) - I_11(y(j),x(i+1))) / area(i,j)
     !ADD PHL END
    end do
  end do  
end subroutine moment_onsphereswap
!END SUBROUTINE MOMENT_ONSPHERESWAP---------------------------------------CE-for FVM!

! ----------------------------------------------------------------------------------!
!SUBROUTINES I_00, I_01, I_20, I_02, I11----------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: calculates the exact integrals                                       !
!                                                                                   !
! CALLS: none                                                                       !
! INPUT: x    ... x coordinate of the evaluation point (Cartesian on the cube)      !
!        y    ... y coordinate of the evaluation point (Cartesian on the cube)      !                                                     !
! OUTPUT: I_00, I_01, I_20, I_02, I11                                               !
!-----------------------------------------------------------------------------------!
function I_00(x,y)
  implicit none
  real (kind=real_kind)                 :: I_00
  real (kind=real_kind), intent(in)     :: x,y

  I_00 = ATAN(x*y/SQRT(1.0D0+x*x+y*y))
end function I_00

function I_10(x,y)
  implicit none
  real (kind=real_kind)                 :: I_10
  real (kind=real_kind), intent(in)     :: x,y
  real (kind=real_kind)                 :: tmp

!   tmp = ATAN(x)
!   I_10 = -ASINH(y*COS(tmp))
  tmp = y*COS(ATAN(x))
  I_10 = -log(tmp+sqrt(tmp*tmp+1))
end function I_10


function I_01(x,y)
  implicit none
  real (kind=real_kind)                 :: I_01
  real (kind=real_kind), intent(in)     :: x,y
  real (kind=real_kind)                 :: tmp

!   I_01 = -ASINH(x/SQRT(1+y*y))
  tmp=x/SQRT(1+y*y)
  I_01 = -log(tmp+sqrt(tmp*tmp+1)) 
end function I_01

function I_20(x,y)
  implicit none
  real (kind=real_kind)                 :: I_20
  real (kind=real_kind), intent(in)     :: x,y
  real (kind=real_kind)                 :: tmp,tmp1

  tmp = 1.0D0+y*y
  tmp1=x/SQRT(tmp)
  I_20 = y*log(tmp1+sqrt(tmp1*tmp1+1))+ACOS(x*y/(SQRT((1.0D0+x*x)*tmp)))
end function I_20

function I_02(x,y)
  implicit none
  real (kind=real_kind)                  :: I_02
  real (kind=real_kind), intent(in)      :: x,y
  real (kind=real_kind)                  :: tmp,tmp1

!   tmp=1.0D0+x*x
!   I_02 = x*ASINH(y/SQRT(tmp))+ACOS(x*y/SQRT(tmp*(1+y*y)))
  tmp=1.0D0+x*x
  tmp1=y/SQRT(tmp) 
  
  I_02 = x*log(tmp1+sqrt(tmp1*tmp1+1))+ACOS(x*y/SQRT(tmp*(1+y*y)))  
  
end function I_02

function I_11(x,y)
  implicit none
  real (kind=real_kind)                   :: I_11
  real (kind=real_kind), intent(in)       :: x,y

  I_11 = -SQRT(1+x*x+y*y)
end function I_11
!END SUBROUTINES I_00, I_01, I_20, I_02, I11------------------------------CE-for FVM!


end module fvm_analytic_mod
