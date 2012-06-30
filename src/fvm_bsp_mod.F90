! MODULE FVM_BSP_MOD------------------------------------------------------CE-for FVM!
! AUTHOR: Christoph Erath, 11.June 2011                                             !
! This module contains everything to run testcases, including fvm_mesh_dep          !
! which computes the departure grid for the selected test case                      !
! SELECTION of TEST CASES just in a COMMENT and UNCOMMENT basis                     !
! Available tests: SOLID Body, Bommerang with slotted cylinders                     !
!-----------------------------------------------------------------------------------! 
module fvm_bsp_mod
  ! ---------------------------------------------------------------------------------
  use kinds, only                  : real_kind
  use dimensions_mod, only: np, nc, ntrac,nlev
  
  implicit none
  private
  public :: fvm_bsp, boomerang, set_boomerang_velocities_gll, get_boomerang_velocities_gll, &
            fvm_init_tracer

  public :: solidbody
  public :: analytical_function

contains


! ----------------------------------------------------------------------------------!
!SUBROUTINE SOLIDBODY-----------------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 11.June 2011                                             !
! DESCRIPTION: calculates the deparute grid for the solid body test case            !
!                                                                                   !
! CALLS: vortex_rotatedsphere, vortex_rotatedsphereback                             !
! INPUT:  asphere ... arrival grid in polar coordinates                             !
! OUTPUT: dsphere ...  departure grid in polar coordinates                          !
!-----------------------------------------------------------------------------------!
subroutine solidbody(asphere, dsphere)
  use kinds, only : real_kind
  use time_mod, only : tstep, nmax, Time_at
  use physical_constants, only : DD_PI
  use fvm_transformation_mod, only: vortex_rotatedsphere, vortex_rotatedsphereback
  use coordinate_systems_mod, only : spherical_polar_t

  implicit none
  type (spherical_polar_t),intent(in)   :: asphere
  type (spherical_polar_t),intent(out)  :: dsphere

! for test case solid-body rotation on the sphere with alpha
  real (kind=real_kind)                 :: alpha, omega
  real (kind=real_kind)                 :: lap,thp,lamrot, therot,tmplondep,tmplatdep
  real (kind=real_kind)                 ::tmpflux
  integer ie,i,j

  ! set values for solid-body rotation on the sphere with alpha, this should be 
  ! outside 
  alpha=1.3!0.78
  omega=2*DD_PI/Time_at(nmax)          ! angular velocity: around the earth
  omega=2*DD_PI/1036800
  lap=DD_PI
  thp=DD_PI/2-alpha

  if (alpha==0.0D0) then      ! move along the equator
    tmplondep=asphere%lon-omega*tstep
    if (tmplondep<0.0D0) then
      tmplondep=tmplondep+2*DD_PI
    endif
    tmplatdep=asphere%lat
  else            
    ! rotate sphere with the pole (lap,thp) with respect to the unrotated system
    call vortex_rotatedsphere(lap,thp,asphere%lon,asphere%lat,lamrot,therot)
    lamrot=lamrot-omega*tstep
!     if (lamrot<0.0D0) lamrot=lamrot+2*DD_PI        
    call vortex_rotatedsphereback(lap,thp,lamrot,therot,&
                                 tmplondep,tmplatdep)                                                      
  endif
  dsphere%lon=tmplondep
  dsphere%lat=tmplatdep
  dsphere%r=asphere%r
end subroutine solidbody
!END SUBROUTINE SOLIDBODY-------------------------------------------------CE-for FVM!

! ----------------------------------------------------------------------------------!
!SUBROUTINE ANALYTICAL_FUNCTION-------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 11.June 2011                                             !
! DESCRIPTION: evaluates the analytical function for the solid body test case       !
!                                                                                   !
! INPUT:  sphere ... coordinates in polar coordinates                               !
!         itr    ... number of the tracer                                           !
! INPUT/OUTPUT: value                                                               !
!-----------------------------------------------------------------------------------!
subroutine analytical_function(value,sphere,klev,itr)
  use coordinate_systems_mod, only: spherical_polar_t
  use physical_constants, only : DD_PI

  implicit none
  type (spherical_polar_t),intent(in)     :: sphere
  integer, intent(in)                     :: klev,itr
  real (kind=real_kind), intent(inout)    :: value 

  real (kind=real_kind)                   :: R0, h0
  real (kind=real_kind)                   :: lon, lat, Rg
  real (kind=real_kind)                   :: lon1, lat1, lon2, lat2, Rg1, Rg2

!   ! temporary: set all parameters for the Solid Body rotation test
!   Solid Body test ----------------------------------------------------------------!
  R0=sphere%r/3.0D0
  h0=1000.0D0
  lon1=4.0D0*DD_PI/5.0D0
  lat1=0.0D0
  Rg1 = acos(sin(lat1)*sin(sphere%lat)+cos(lat1)*cos(sphere%lat)*cos(sphere%lon-lon1))
  lon2=6.0D0*DD_PI/5.0D0
  lat2=0.0D0
  Rg2 = acos(sin(lat2)*sin(sphere%lat)+cos(lat2)*cos(sphere%lat)*cos(sphere%lon-lon2))
  if (Rg1 .le. R0) then
    value = 1.0D0+(h0/2.0D0)*(1.0D0+cos(DD_PI*Rg1/R0))
!     value=sin(sphere%lon)*cos(sphere%lat)
  elseif (Rg2 .le. R0) then
    value = 1.0D0+(h0/2.0D0)*(1.0D0+cos(DD_PI*Rg2/R0))
  else
    value = 1.0D0
  endif



  !Non-smooth scalar field (slotted cylinder) --------------------------------------!
!   R0=0.5D0
!   lon1=4.0D0*DD_PI/5.0D0
!   lat1=0.0D0
!   Rg1 = acos(sin(lat1)*sin(sphere%lat)+cos(lat1)*cos(sphere%lat)*cos(sphere%lon-lon1))
!   lon2=6.0D0*DD_PI/5.0D0
!   lat2=0.0D0
!   Rg2 = acos(sin(lat2)*sin(sphere%lat)+cos(lat2)*cos(sphere%lat)*cos(sphere%lon-lon2))   
! 
!   if ((Rg1 .le. R0) .AND. (abs(sphere%lon-lon1).ge. R0/6)) then
!     value = 1.0D0 + itr*1.0D0/ntrac + klev*1.0D0/nlev
!   elseif ((Rg2 .le. R0) .AND. (abs(sphere%lon-lon2).ge. R0/6)) then
!     value = 1.0D0 + itr*1.0D0/ntrac + klev*1.0D0/nlev
!   elseif ((Rg1 .le. R0) .AND. (abs(sphere%lon-lon1) < R0/6) &
!                         .AND. (sphere%lat-lat1 < -5.0D0*R0/12.0D0)) then
!     value = 1.0D0 + itr*1.0D0/ntrac + klev*1.0D0/nlev
!   elseif ((Rg2 .le. R0) .AND. (abs(sphere%lon-lon2) < R0/6) &
!                         .AND. (sphere%lat-lat2 > 5.0D0*R0/12.0D0)) then
!     value = 1.0D0 + itr*1.0D0/ntrac + klev*1.0D0/nlev
!   else
!     value = 0.1D0 + itr*1.0D0/ntrac + klev*1.0D0/nlev   
!   endif  
end subroutine analytical_function
!END SUBROUTINE ANALYTICAL_FUNCTION---------------------------------------CE-for FVM!



! ----------------------------------------------------------------------------------!
!SUBROUTINE FVM_BSP-------------------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 11.June 2011                                             !
! DESCRIPTION: controls the test cases                                              !
!                                                                                   !
! CALLS: analytical_function                                                        !
! INPUT: elem     ...  element structure from HOMME                                 !
!        tl       ...  time level structure                                         !
! INTPUT/OUTPUT:                                                                    !
!        fvm   ...  structure                                                       !
!-----------------------------------------------------------------------------------!
subroutine fvm_bsp(fvm, tl)
  use element_mod, only : element_t
  use time_mod, only : timelevel_t
  use fvm_control_volume_mod, only: fvm_struct
  implicit none
  type (fvm_struct), intent(inout)      :: fvm
  type (timeLevel_t), intent(in)       :: tl              ! time level struct
  
  integer                              :: i,j,k,itr
  
  do k=1, nlev
    fvm%c(:,:,k,1,tl%n0)=1.0D0    !density of the air
    do itr=2,ntrac
      do j=1,nc
        do i=1,nc               
          call analytical_function(fvm%c(i,j,k,itr,tl%n0),fvm%centersphere(i,j),k,itr)      
        end do
      end do
    end do 
  end do
end subroutine fvm_bsp
!END SUBROUTINE FVM_BSP---------------------------------------------------CE-for FVM!

function get_boomerang_velocities_gll(elem, time) result(vstar)
  use coordinate_systems_mod, only : cart2cubedspherexy, spherical_to_cart
  use element_mod, only : element_t
  use fvm_control_volume_mod, only: fvm_struct
  use physical_constants, only : DD_PI, rearth
  
  implicit none
  type (element_t), intent(in)   :: elem
  real (kind=real_kind),intent(in)             :: time  ! time of the arrival grid
  real (kind=real_kind)                        :: vstar(np,np,2)
  
  integer                              :: i,j
  real (kind=real_kind)       :: slat, clat, slon, clon, &
                                 slon2,clon2, u, v, tt, ck, omega, lon, lat, tmp_time
  tt = 5.0D0  !total time
  ck = 10.0D0/tt
  omega = DD_PI/tt

  !tmp_time=(nstep)*5.0D0/nmax
  tmp_time=5*time/(12*3600*24)  ! convert from dimensional time to dimensionless
  
  do i=1,np
    do j=1,np
      lon = elem%spherep(i,j)%lon
      lat = elem%spherep(i,j)%lat

      slat  = sin(lat)
      clat  = cos(lat)

      slon  = sin(lon-tmp_time*2*DD_PI/tt)!solid-body rotation added
      slon2 = sin(2*(lon-tmp_time*2*DD_PI/tt))!solid-body rotation added
      clon  = cos(lon-tmp_time*2*DD_PI/tt)!solid-body rotation added
      clon2 = cos(2*(lon-tmp_time*2*DD_PI/tt))!solid-body rotation added

      u =  ck*slon*slon*sin(2*lat)*cos(tmp_time*omega) + clat*2*DD_PI/(tt)
      v =  ck*slon2*clat*cos(tmp_time*omega)
      ! convert from radians per dimensionless time to 
      ! meters/sec 
      vstar(i,j,1)=u * Rearth /( 12*3600*24/5)
      vstar(i,j,2)=v * Rearth /( 12*3600*24/5)
    enddo
  enddo

end function get_boomerang_velocities_gll     

subroutine set_boomerang_velocities_gll(elem, time,klev)
  use coordinate_systems_mod, only : cart2cubedspherexy, spherical_to_cart
  use element_mod, only : element_t
  use fvm_control_volume_mod, only: fvm_struct
  use physical_constants, only : DD_PI, rearth
  
  implicit none
  type (element_t), intent(inout)   :: elem
  real (kind=real_kind),intent(in)             :: time  ! time of the arrival grid
  integer,intent(in)                   :: klev
  
  integer                              :: i,j
  real (kind=real_kind)       :: slat, clat, slon, clon, &
                                 slon2,clon2, u, v, tt, ck, omega, lon, lat, tmp_time
  tt = 5.0D0  !total time
  ck = 10.0D0/tt
  omega = DD_PI/tt

  !tmp_time=(nstep)*5.0D0/nmax
  tmp_time=5*time/(12*3600*24)  ! convert from dimensional time to dimensionless
  
  do i=1,np
    do j=1,np
      lon = elem%spherep(i,j)%lon
      lat = elem%spherep(i,j)%lat

      slat  = sin(lat)
      clat  = cos(lat)

      slon  = sin(lon-tmp_time*2*DD_PI/tt)!solid-body rotation added
      slon2 = sin(2*(lon-tmp_time*2*DD_PI/tt))!solid-body rotation added
      clon  = cos(lon-tmp_time*2*DD_PI/tt)!solid-body rotation added
      clon2 = cos(2*(lon-tmp_time*2*DD_PI/tt))!solid-body rotation added

      u =  ck*slon*slon*sin(2*lat)*cos(tmp_time*omega) + clat*2*DD_PI/(tt)
      v =  ck*slon2*clat*cos(tmp_time*omega)
      ! convert from radians per dimensionless time to 
      ! meters/sec 
      elem%derived%vstar(i,j,1,klev)=u * Rearth /( 12*3600*24/5)
      elem%derived%vstar(i,j,2,klev)=v * Rearth /( 12*3600*24/5)
    enddo
  enddo

end subroutine set_boomerang_velocities_gll



!END SUBROUTINE FVM_MESH_DEP----------------------------------------------CE-for FVM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE BOOMERANG-----------------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 31. October 2011                                         !
! DESCRIPTION: calculates the deparute grid for the boomerang flow                  !
!                                                                                   !
! CALLS: vortex_rotatedsphere, vortex_rotatedsphereback                             !
! INPUT:  asphere ... arrival grid in polar coordinates                             !
!         nstep   ... actual step                                                   !
! OUTPUT: dsphere ...  departure grid in polar coordinates                          !
!-----------------------------------------------------------------------------------!
subroutine boomerang(asphere,dsphere,nstep)
  use physical_constants, only : DD_PI
  use coordinate_systems_mod, only : spherical_polar_t
  use time_mod, only : tstep,nmax

  implicit none
  type (spherical_polar_t),intent(in)   :: asphere
  type (spherical_polar_t),intent(out)  :: dsphere
  integer, intent(in)                   :: nstep
  
  
  integer                     :: iorder, i, iteration
  real (kind=real_kind)       :: slat, clat, slon, clon, xg, yg, zg, ca, sa, co, so, &
                                 slon2,clon2,slonh,clonh
  real (kind=real_kind)       :: tmp_time,lon,lat,tmp_dt,tmp_lm,tmp_th, ck, tt, omega
  real (kind=real_kind)       :: dplm,dpth,trm1,trm2,sslm,cwt,swt
  real (kind=real_kind)       :: dt2,dt3,dt4,dtt,udc,uexact, vexact

  iteration=10
  iorder=5
  
  tmp_time=(nstep)*5.0D0/nmax
  tmp_dt = 5.0D0/nmax/iteration
  dt2 = tmp_dt*tmp_dt/2.0D0
  dt3 = dt2*tmp_dt/3.0D0
  dt4 = dt3*tmp_dt/4.0D0
  
  tmp_lm = asphere%lon
  tmp_th = asphere%lat
  tt = 5.0D0  !total time
  ck = 10.0D0/tt
  
  DO i=1,iteration

    omega = DD_PI/tt
    IF (i>1) THEN
      tmp_time=tmp_time-tmp_dt
    END IF 
     
    cwt = cos(omega*tmp_time)
    swt = sin(omega*tmp_time)

    lon   = tmp_lm
    lat   = tmp_th
    slat  = sin(lat)
    clat  = cos(lat)
    
    slon  = sin(lon-tmp_time*2*DD_PI/tt)!solid-body rotation added
    slon2 = sin(2*(lon-tmp_time*2*DD_PI/tt))!solid-body rotation added
    clon  = cos(lon-tmp_time*2*DD_PI/tt)!solid-body rotation added
    clon2 = cos(2*(lon-tmp_time*2*DD_PI/tt))!solid-body rotation added
    sslm  = (sin(0.5*(lon-tmp_time*2*DD_PI/tt)))**2

    uexact =  ck*slon*slon*sin(2*lat)*cos(tmp_time*omega) + clat*2*DD_PI/(tt)
    vexact =  ck*slon2*clat*cos(tmp_time*omega)
    udc    =  2*ck*slon*slon*slat*cos(tmp_time*omega) + 2*DD_PI/tt
    ! 2nd-order
    tmp_lm = lon - tmp_dt * udc 
    
    ! 3rd-order
    if (iorder>2) then
         tmp_lm = tmp_lm - 2 * dt2 * slon * (slon * swt * omega * slat - 2 * udc * cwt * clon * slat&
         - vexact * slon * cwt * clat) * ck 
    endif
    ! 4th-order
    if (iorder>3) then
         tmp_lm = tmp_lm - dt3 * (-4 * slon * swt * (2 * udc * clon * slat + vexact * slon * clat) * &
         ck * omega + 2 * cwt * (-omega ** 2 * slat + omega ** 2 * slat * clon ** 2 &
         + 4 * udc ** 2 * clon ** 2 * slat - 2 * udc ** 2 * slat + 4 * udc * vexact * &
         slon * clat * clon - vexact ** 2 * slat + vexact ** 2 * slat * clon ** 2) * ck)
    endif
    !
    ! 5th-order
    !
    if (iorder>4) then
         tmp_lm = tmp_lm + dt4 * (-2 * swt * (-slat * omega ** 2 + omega ** 2 * slat * clon ** 2 + &
         12 * udc ** 2 * clon ** 2 * slat - 6 * udc ** 2 * slat + 12 * udc * vexact * &
         slon * clon * clat - 3 * vexact ** 2 * slat + 3 * vexact ** 2 * slat * clon ** 2) &
         * ck * omega - 2 * cwt * (vexact ** 3 * clat - 3 * vexact * omega ** 2 * clat * &
         clon ** 2 - vexact ** 3 * clat * clon ** 2 - 8 * udc ** 2 * vexact * clon ** 2 * clat &
         + 6 * udc * vexact ** 2 * slat * slon * clon + 3 * vexact * omega ** 2 * clat + 6 * &
         udc ** 2 * vexact * clat + 8 * udc ** 3 * slon * clon * slat + 6 * udc * slon * &
         omega ** 2 * clon * slat) * ck)
    endif
    ! 2nd-order
    tmp_th = lat - tmp_dt * vexact 
    ! 3rd-order
    if (iorder>2) then
         tmp_th = tmp_th &
!- 2 * dt2 * ck * (clat * swt * omega * slon * clon - 2 * udc * clat * cwt * &
!             clon ** 2 + udc * clat * cwt + vexact * slat * cwt * slon * clon)
         + dt2 * (-ck * slon2 * clat * swt * omega + 2 * udc * ck * clon2 *&
         clat * cwt - vexact * ck * slon2 * slat * cwt) 
    endif
    ! 4th-order
    if (iorder>3) then    
        tmp_th = tmp_th - dt3 * (-4 * omega * ck * swt * (2 * udc * clat * clon ** 2 &
         - udc * clat - vexact * slat * slon * clon) - 2 * ck * cwt * (clat * omega ** 2 * &
         slon * clon + 4 * udc ** 2 * slon * clat * clon + 12 * udc * vexact * slat &
         * clon ** 2 - 2 * udc * vexact * slat + vexact ** 2 * slon * clat * clon))
    endif
    ! 5th-order
    if (iorder>4) then
      tmp_th = tmp_th + dt4 * (2 * omega * ck * swt * (clat * omega ** 2 * slon * clon - 6 * &
      udc * vexact * slat + 12 * udc ** 2 * slon * clon * clat + 12 * udc * vexact * &
      slat * clon ** 2 + 3 * vexact ** 2 * slon * clon * clat) - 2 * ck * cwt * &
      (-vexact ** 3 * slat * slon * clon - 3 * udc * clat * omega ** 2 + 6 * udc * clat &
      * omega ** 2 * clon ** 2 - 3 * vexact * slat * omega ** 2 * slon * clon + 6 * &
      udc ** 3 * clat * clon ** 2 - 12 * udc ** 2 * vexact * slat * slon * clon + 6 * &
      udc * vexact ** 2 * clon ** 2 * clat - 4 * udc * vexact ** 2 * clat - 4 * &
      udc ** 3 * clat))
    endif
  end do
  dsphere%lon=tmp_lm
  dsphere%lat=tmp_th
  dsphere%r=asphere%r
end subroutine boomerang
!END SUBROUTINE BOOMERANG-------------------------------------------------CE-for FVM!

! ----------------------------------------------------------------------------------!
!SUBROUTINE FVM_INIT_TRACER-----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 11.June 2011                                             !
! DESCRIPTION: controls the test cases                                              !
!                                                                                   !
! CALLS: analytical_function                                                        !
! INPUT: elem     ...  element structure from HOMME                                 !
!        tl       ...  time level structure                                         !
! INTPUT/OUTPUT:                                                                    !
!        fvm   ...  structure                                                       !
!-----------------------------------------------------------------------------------!
subroutine fvm_init_tracer(fvm, tl)
  use element_mod, only : element_t
  use time_mod, only : timelevel_t
  use fvm_control_volume_mod, only: fvm_struct
  implicit none
  type (fvm_struct), intent(inout)      :: fvm
  type (timeLevel_t), intent(in)       :: tl              ! time level struct
  
  integer                              :: i,j,k,itr
  
  do k=1, nlev
    fvm%c(:,:,k,1,tl%n0)=1.0D0    !density of the air
!     fvm%c(:,:,k,2,tl%n0)=10.0D0    !density of the air
    
    do itr=2,ntrac
      do j=1,nc
        do i=1,nc               
          call analytical_function(fvm%c(i,j,k,itr,tl%n0),fvm%centersphere(i,j),k,itr)      
        end do
      end do
    end do 
  end do
end subroutine fvm_init_tracer
!END SUBROUTINE FVM_BSP---------------------------------------------------CE-for FVM!

end module fvm_bsp_mod
