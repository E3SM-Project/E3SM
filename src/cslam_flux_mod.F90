#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!CE-for CSLAM-----------------------------------------------------------------------!
! AUTHOR: Christoph Erath, 11.June 2011                                             !
! This modules contains the flux version of CSLAM                                   !
! AVAILABLE: only the simplified version is implemented (14.June 2011)              !
!-----------------------------------------------------------------------------------!
  
module cslam_flux_mod
  ! ---------------------------------------------------------------------------------
  use kinds, only : real_kind
  ! --------------------------------------------------------------------------------- 
  implicit none
  private
  
  public :: simplified_flux

contains

!CE-for CSLAM-----------------------------------------------------------------------!
! AUTHOR: Christoph Erath, 11.June 2011                                             !
!-----------------------------------------------------------------------------------!
subroutine simplified_flux(fcube,cslam, tl, nets, nete)
  use dimensions_mod, only: nc, nhc
  use cslam_control_volume_mod, only: cslam_struct
  use time_mod, only : timelevel_t
  ! --------------------------------------------------------------------------------- 
  implicit none
  real (kind=real_kind), dimension(:,:,:,:), intent(inout) :: fcube
  
  type (cslam_struct),intent(inout)     :: cslam(:)
  type (TimeLevel_t), intent(in)        :: tl              ! time level struct
  integer, intent(in)                   :: nets, nete
  
  real (kind=real_kind)                 :: tmpflux
  integer                               :: ie,i,j

  do ie=nets, nete 
    do j=1,nc
      do i=1,nc+1
        ! xfluxes
        call flux_areas_x(cslam(ie)%acart(i,j), cslam(ie)%acart(i,j+1), &
                          cslam(ie)%dcart(i,j), cslam(ie)%dcart(i,j+1), &
                          fcube(i-1,j,1,tl%n0),fcube(i,j,1,tl%n0),&
                          cslam(ie)%asphere(i,j), cslam(ie)%asphere(i,j+1), &
                          cslam(ie)%dsphere(i,j), cslam(ie)%dsphere(i,j+1),cslam(ie)%faceno,tmpflux)
        fcube(i-1,j,1,tl%np1)=fcube(i-1,j,1,tl%np1)-tmpflux
        fcube(i,j,1,tl%np1)=fcube(i,j,1,tl%np1)+tmpflux    
        !yfluxes         
        call flux_areas_y(cslam(ie)%acart(j,i), cslam(ie)%acart(j+1,i), &
                          cslam(ie)%dcart(j,i), cslam(ie)%dcart(j+1,i), &
                          fcube(j,i-1,1,tl%n0),fcube(j,i,1,tl%n0),&
                          cslam(ie)%asphere(j,i), cslam(ie)%asphere(j+1,i), &
                          cslam(ie)%dsphere(j,i), cslam(ie)%dsphere(j+1,i),cslam(ie)%faceno,tmpflux)
                
        fcube(j,i-1,1,tl%np1)=fcube(j,i-1,1,tl%np1)+tmpflux
        fcube(j,i,1,tl%np1)=fcube(j,i,1,tl%np1)-tmpflux 
       end do
    end do 
    ! finalize scheme
    do j=1,nc
      do i=1,nc        
        fcube(i,j,1,tl%np1)=fcube(i,j,1,tl%n0)+   &
                 (fcube(i,j,1,tl%np1))/cslam(ie)%area_sphere(i,j)  
      end do
    end do     
  end do
end subroutine simplified_flux
!-----------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------!
subroutine flux_areas_x(ap1,ap2,dp1,dp2,cvalleft,cvalright,&
                        asph1,asph2,dsph1,dsph2, faceno, value)
  use coordinate_systems_mod, only: cartesian2D_t,spherical_polar_t, unit_face_based_cube_to_unit_sphere
  
  implicit none
  logical                                :: lcross, lcrosstest,orient
  type(cartesian2D_t), intent (in)       :: ap1,ap2,dp1, dp2
  type (spherical_polar_t), intent (in)  :: asph1,asph2,dsph1,dsph2
  real (kind=real_kind),intent(in)       :: cvalleft, cvalright
  integer, intent(in)                    :: faceno
  real (kind=real_kind),intent(out)      :: value
  
  real (kind=real_kind)                  :: lat(4), lon(4), latt1(3), lont1(3), latt2(3), lont2(3)
  type(cartesian2D_t)                    :: interpointxy
  type(spherical_polar_t)                :: interpointsphere
  
  
  call check_lines_cross(ap1,dp1,ap2,dp2,lcrosstest,interpointxy)
  if (lcrosstest) then
    print *, "FATAL Error in flux_areas_x, lines between departure and arrival node cross!"
    STOP "Exit program!"
  endif
  
  call check_lines_cross(ap1,ap2,dp1,dp2,lcross,interpointxy)
  value=0
  if ((.not. lcross)) then  
    call point_orientation4(ap1,ap2,dp2,dp1,orient)
    if(orient) then
      lat(1)=asph1%lat
      lat(2)=asph2%lat
      lat(3)=dsph2%lat
      lat(4)=dsph1%lat   
      lon(1)=asph1%lon
      lon(2)=asph2%lon
      lon(3)=dsph2%lon
      lon(4)=dsph1%lon
      value=cvalleft*abs(spherical_area_rect(lat, lon))
    else
      lat(1)=asph1%lat
      lat(2)=dsph1%lat
      lat(3)=dsph2%lat
      lat(4)=asph2%lat   
      lon(1)=asph1%lon
      lon(2)=dsph1%lon
      lon(3)=dsph2%lon
      lon(4)=asph2%lon
      value=-cvalright*abs(spherical_area_rect(lat, lon))  
    endif
  else    
    interpointsphere=unit_face_based_cube_to_unit_sphere(interpointxy,faceno)
    latt1(1)=asph1%lat
    latt1(2)=interpointsphere%lat
    latt1(3)=dsph1%lat
    lont1(1)=asph1%lon
    lont1(2)=interpointsphere%lon
    lont1(3)=dsph1%lon
    
    latt2(1)=asph2%lat
    latt2(2)=interpointsphere%lat
    latt2(3)=dsph2%lat
    lont2(1)=asph2%lon
    lont2(2)=interpointsphere%lon
    lont2(3)=dsph2%lon
    call point_orientation3(ap1,interpointxy,dp1,orient) 
    if (orient) then
      value=cvalleft*abs(spherical_area_tria(latt1, lont1)) - &
            cvalright*abs(spherical_area_tria(latt2, lont2))
    else
      value=-cvalright*abs(spherical_area_tria(latt1, lont1)) + &
            cvalleft*abs(spherical_area_tria(latt2, lont2))
    endif
  endif
end subroutine flux_areas_x
!-----------------------------------------------------------------------------------!
subroutine flux_areas_y(ap1,ap2,dp1,dp2,cvallower,cvalupper, &
                        asph1,asph2,dsph1,dsph2,faceno, value)
  use coordinate_systems_mod, only: cartesian2D_t,spherical_polar_t, unit_face_based_cube_to_unit_sphere
  implicit none
  logical                                :: lcross, lcrosstest, orient
  type(cartesian2D_t), intent (in)       :: ap1,ap2,dp1, dp2
  type (spherical_polar_t), intent (in)  :: asph1,asph2,dsph1,dsph2
  real (kind=real_kind),intent(in)       :: cvallower, cvalupper
  integer, intent(in)                    :: faceno
  real (kind=real_kind),intent(out)      :: value
  
  real (kind=real_kind)                  :: lat(4), lon(4), latt1(3), lont1(3), latt2(3), lont2(3),test1, test2
  type(cartesian2D_t)                    :: interpointxy
  type(spherical_polar_t)                :: interpointsphere  
  
  call check_lines_cross(ap1,dp1,ap2,dp2,lcrosstest,interpointxy)
  if (lcrosstest) then
    print *, "FATAL Error in flux_areas_y, lines between departure and arrival node cross!"
    STOP "Exit program!"
  endif
  call check_lines_cross(ap1,ap2,dp1,dp2,lcross,interpointxy)
  value=0
  if ((.not. lcross)) then
    call point_orientation4(ap1,ap2,dp2,dp1,orient)
     if(orient) then
       lat(1)=asph1%lat
       lat(2)=asph2%lat
       lat(3)=dsph2%lat
       lat(4)=dsph1%lat   
       lon(1)=asph1%lon
       lon(2)=asph2%lon
       lon(3)=dsph2%lon
       lon(4)=dsph1%lon
       value=cvalupper*abs(spherical_area_rect(lat, lon))        
     else
       lat(1)=asph1%lat
       lat(2)=dsph1%lat
       lat(3)=dsph2%lat
       lat(4)=asph2%lat   
       lon(1)=asph1%lon
       lon(2)=dsph1%lon
       lon(3)=dsph2%lon
       lon(4)=asph2%lon
       value=-cvallower*abs(spherical_area_rect(lat, lon))       
     endif
   else
     interpointsphere=unit_face_based_cube_to_unit_sphere(interpointxy,faceno)
     latt1(1)=asph1%lat
     latt1(2)=interpointsphere%lat
     latt1(3)=dsph1%lat
     lont1(1)=asph1%lon
     lont1(2)=interpointsphere%lon
     lont1(3)=dsph1%lon

     latt2(1)=asph2%lat
     latt2(2)=interpointsphere%lat
     latt2(3)=dsph2%lat
     lont2(1)=asph2%lon
     lont2(2)=interpointsphere%lon
     lont2(3)=dsph2%lon
     call point_orientation3(ap1,interpointxy,dp1,orient) 
     if (orient) then
       value=cvalupper*abs(spherical_area_tria(latt1, lont1)) - &
             cvallower*abs(spherical_area_tria(latt2, lont2))
     else
       value=-cvallower*abs(spherical_area_tria(latt1, lont1)) + &
             cvalupper*abs(spherical_area_tria(latt2, lont2))
     endif
   endif 
end subroutine flux_areas_y
!-----------------------------------------------------------------------------------!
subroutine point_orientation4(p1,p2,p3,p4,orient)
  use coordinate_systems_mod, only: cartesian2D_t
  implicit none
  type(cartesian2D_t), intent (in)    :: p1,p2,p3,p4
  
  logical, intent(out)                :: orient
  real (kind=real_kind)               :: tmp
  ! 
  tmp=(p1%y+p2%y)*(p1%x-p2%x)+(p2%y+p3%y)*(p2%x-p3%x)+ &
             (p3%y+p4%y)*(p3%x-p4%x)+(p4%y+p1%y)*(p4%x-p1%x)
  orient=.FALSE.
  if (tmp>0.0D0) then
     orient=.TRUE.   !points are counter clockweise
   endif
end subroutine point_orientation4
!-----------------------------------------------------------------------------------!
subroutine point_orientation3(p1,p2,p3,orient)
  use coordinate_systems_mod, only: cartesian2D_t
  implicit none
  type(cartesian2D_t), intent (in)    :: p1, p2, p3 
  logical, intent(out)                :: orient
  real (kind=real_kind)               :: tmp
  ! 
  tmp=(p1%y+p2%y)*(p1%x-p2%x)+(p2%y+p3%y)*(p2%x-p3%x)+ &
             (p3%y+p1%y)*(p3%x-p1%x)
  orient=.FALSE.
  if (tmp>0.0D0) then
     orient=.TRUE.   !points are counter clockweise
   endif
end subroutine point_orientation3
!-----------------------------------------------------------------------------------!
! see, e.g., http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
subroutine check_lines_cross(p1,p2,q1,q2,lcross, point)
  use coordinate_systems_mod, only: cartesian2D_t
  implicit none
  type(cartesian2D_t), intent (in)  :: p1,p2,q1,q2
  type(cartesian2D_t), intent (out) :: point
  logical, intent(out)              :: lcross
  !
  ! local workspace
  !
  REAL (kind=real_kind)    :: dn,tp,tq, NaN=-1.0

  NaN=sqrt(NaN)
  dn = (q2%y-q1%y)*(p2%x-p1%x)-(q2%x-q1%x)*(p2%y-p1%y)
  point%x=NaN
  point%y=NaN
  if (abs(dn)<1.0D-12) then
    ! lines are parallel
    lcross = .false.
  else
    tp = ((q2%x-q1%x)*(p1%y-q1%y)-(q2%y-q1%y)*(p1%x-q1%x))/dn
    tq = ((p2%x-p1%x)*(p1%y-q1%y)-(p2%y-p1%y)*(p1%x-q1%x))/dn
    ! implement the next two lines faster!!!!!!
    if (tp>-1.0D-12 .and. tp<1.0D00+1.0D-12 .and. &
        tq>-1.0D-12 .and. tq<1.0D00+1.0D-12) then
      lcross = .true.
      point%x=p1%x+tp*(p2%x-p1%x)
      point%y=p1%y+tp*(p2%y-p1%y)
    else
      ! not parallel but not crossing (inside)
      lcross = .false.
    endif
  endif
end subroutine check_lines_cross

function spherical_area_rect(lat, lon) result(area)
  implicit none
  real (kind=real_kind), intent(in)    :: lat(:), lon(:)
  real (kind=real_kind)                :: area
  
  integer                              :: i
  real (kind=real_kind)          :: dlon, dlat, side(5), sum, exc1,exc2,tmp
  
  do i=1,3
    dlon=lon(i+1)-lon(i)
    dlat=lat(i+1)-lat(i)
    side(i)=sin(dlat/2)*sin(dlat/2)+cos(lat(i))*cos(lat(i+1))*sin(dlon/2)*sin(dlon/2)
    side(i)=2*atan2(sqrt(side(i)),sqrt(1-side(i)))
  end do
  dlon=lon(1)-lon(4)
  dlat=lat(1)-lat(4)
  side(4)=sin(dlat/2)*sin(dlat/2)+cos(lat(4))*cos(lat(1))*sin(dlon/2)*sin(dlon/2)
  side(4)=2*atan2(sqrt(side(4)),sqrt(1-side(4)))  
    
  dlon=lon(4)-lon(2)
  dlat=lat(4)-lat(2)
  side(5)=sin(dlat/2)*sin(dlat/2)+cos(lat(2))*cos(lat(4))*sin(dlon/2)*sin(dlon/2)
  side(5)=2*atan2(sqrt(side(5)),sqrt(1-side(5)))  
  
  sum=(side(1)+side(5)+side(4))/2
  tmp=tan(sum/2)*tan((sum-side(1))/2)*tan((sum-side(5))/2)*tan((sum-side(4))/2)
  tmp=sqrt(abs(tmp))  
  exc1=4*atan(tmp)
  sum=(side(2)+side(3)+side(5))/2
  tmp=tan(sum/2)*tan((sum-side(2))/2)*tan((sum-side(3))/2)*tan((sum-side(5))/2)
  tmp=sqrt(abs(tmp))
  exc2=4*atan(tmp)
  area=exc1+exc2  
end function spherical_area_rect

function spherical_area_tria(lat, lon) result(area)
  implicit none
  real (kind=real_kind), intent(in)    :: lat(:), lon(:)
  real (kind=real_kind)                :: area
  
  integer                              :: i
  real (kind=real_kind)          :: dlon, dlat, side(3), sum, tmp
  
  do i=1,2
    dlon=lon(i+1)-lon(i)
    dlat=lat(i+1)-lat(i)
    side(i)=sin(dlat/2)*sin(dlat/2)+cos(lat(i))*cos(lat(i+1))*sin(dlon/2)*sin(dlon/2)
    side(i)=2*atan2(sqrt(side(i)),sqrt(1-side(i)))
  end do
  dlon=lon(1)-lon(3)
  dlat=lat(1)-lat(3)
  side(3)=sin(dlat/2)*sin(dlat/2)+cos(lat(3))*cos(lat(1))*sin(dlon/2)*sin(dlon/2)
  side(3)=2*atan2(sqrt(side(3)),sqrt(1-side(3)))  
  
  sum=(side(1)+side(2)+side(3))/2
  tmp=tan(sum/2)*tan((sum-side(1))/2)*tan((sum-side(2))/2)*tan((sum-side(3))/2)
  tmp=sqrt(abs(tmp))
  area=4*atan(tmp)    
end function spherical_area_tria

end module cslam_flux_mod
