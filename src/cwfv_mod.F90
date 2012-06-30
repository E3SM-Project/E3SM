!-----------------------------------------------------------------------------------!
!MODULE FVM_MOD-----------------------------------------------------------CE-for FVM!
! FVM_MOD File for the fvm project in HOMME                                         !
! Author: Christoph Erath                                                           !
! Date: 25.January 2011                                                             !
! MAIN module to run fvm on HOMME                                                   !
! 14.November 2011: reorganisation done                                             !
! 7.Februar 2012: cslam_run and cslam_runair                                            !
!-----------------------------------------------------------------------------------!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module cwfv_mod

  use kinds, only : real_kind, int_kind, longdouble_kind
  use edge_mod, only : ghostbuffertr_t, initghostbuffer, freeghostbuffertr, &
                       ghostVpack, ghostVunpack,  edgebuffer_t, initEdgebuffer
  use bndry_mod, only: ghost_exchangeV                     
  use dimensions_mod, only: nelem, nelemd, nelemdmax, nlev, ne, nlev, ntrac, np, ntrac_d
  use time_mod, only : timelevel_t
  use coordinate_systems_mod, only : spherical_polar_t, cartesian2D_t
  use element_mod, only : element_t, timelevels
  use hybrid_mod, only : hybrid_t

  implicit none
  private
  type (ghostBuffertr_t)                      :: cellghostbuf
  type (EdgeBuffer_t)                         :: edgeveloc
  
    type, public :: cwfv_struct
      sequence
      ! fvm tracer mixing ratio: (kg/kg)
      real (kind=real_kind)    :: c(1-np:np+np,1-np:np+np,nlev,ntrac_d,timelevels) 
      real (kind=real_kind)    :: vn0(np,np,2,nlev)
  !-----------------------------------------------------------------------------------!
      ! define the departure grid (depends on the examples, usually the velocity)
      ! they are NOT stored for the Levels
      type (spherical_polar_t) :: dsphere(np,np)     ! Spherical coordinates
      type (cartesian2D_t)     :: dcart(np,np)
      
      integer                  :: xstart(np,np)
      integer                  :: ystart(np,np)
  !-----------------------------------------------------------------------------------!  
      !area of the cell on the sphere 
      integer                  :: faceno
      integer                  :: nbrsface(8)    ! store the neighbours in north, south 
      ! number of south,....,swest and 0 for interior element 
      integer                  :: cubeboundary                                                 
  !-----------------------------------------------------------------------------------!     
      ! next maybe used just for test reasons
      real (kind=real_kind)    :: elem_mass
      real (kind=real_kind)    :: maxcfl(2,nlev)   
    end type cwfv_struct
  
  public :: cellghostbuf, edgeveloc, fvm_init1,fvm_init2, fvm_mcgregor, fvm_mcgregordss
  public :: fvm_init3, cwfv_run
contains

subroutine cwfv_run(elem,fvm,hybrid,deriv,tstep,tl,nets,nete)
  ! ---------------------------------------------------------------------------------
  use derivative_mod, only : derivative_t
  ! ---------------------------------------------------------------------------------
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------
  use interpolate_mod, only : interpolate_ce
   
  implicit none
  type (element_t), intent(inout)                :: elem(:)
  type (cwfv_struct), intent(inout)             :: fvm(:)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)
  real (kind=real_kind)                       :: tstep
  
  integer                                     :: i,j,k,ie,itr, jx, jy, jdx, jdy, h
  type (TimeLevel_t)                          :: tl              ! time level struct
  type (derivative_t)                         :: deriv           ! derivative struct
  
  integer            :: xstart, ystart, xend, yend  
   
  do ie=nets, nete
    do k=1,nlev
      call fvm_mesh_dep_gll(elem(ie),deriv,fvm(ie),tstep,tl,k)
      !-Departure fvm Meshes, initialization done                                                                  
      do itr=1,ntrac  
         do j=1,np
            do i=1,np
               xstart=fvm(ie)%xstart(i,j)
               xend=xstart+np-1
               ystart=fvm(ie)%ystart(i,j)
               yend=ystart+np-1
               
               call interpolate_ce(fvm(ie)%dcart(i,j), &
                    fvm(ie)%c(xstart:xend,ystart:yend,k,itr,tl%n0), &
                    np,fvm(ie)%c(i,j,k,itr,tl%np1))
            end do
         end do
      enddo  !End Tracer
    end do  !End Level
    !note write tl%np1 in buffer
    call ghostVpack(cellghostbuf, fvm(ie)%c,np,np,nlev,ntrac,0,tl%np1,timelevels,elem(ie)%desc)
  end do
  call t_startf('FVM Communication')
  call ghost_exchangeV(hybrid,cellghostbuf,np,np)
  call t_stopf('FVM Communication')
  !-----------------------------------------------------------------------------------!
  do ie=nets,nete
     call ghostVunpack(cellghostbuf, fvm(ie)%c, np, np,nlev,ntrac, 0, tl%np1, timelevels,elem(ie)%desc)
  enddo
end subroutine cwfv_run


! initialize global buffers shared by all threads
subroutine fvm_init1(par)
  use parallel_mod, only : parallel_t, haltmp
  type (parallel_t) :: par
  
  if (ntrac>ntrac_d) then
     if (par%masterproc) print *,'ntrac,ntrac_d=',ntrac,ntrac_d
     call haltmp("PARAMTER ERROR for fvm: ntrac > ntrac_d")
  endif

  call initghostbuffer(cellghostbuf,nlev,ntrac,np,np) !+1 for the air_density, which comes from SE
  
  call initEdgebuffer(edgeveloc,2*nlev)
end subroutine fvm_init1



! initialization that can be done in threaded regions
subroutine fvm_init2(elem,fvm,hybrid,nets,nete,tl)
  use fvm_control_volume_mod, only: fvm_mesh_ari
  use fvm_analytic_mod, only: computexytosphere_moments
  use bndry_mod, only: compute_ghost_corner_orientation 
  use parallel_mod, only : abortmp

  type (timelevel_t)                    :: tl
  type (cwfv_struct)                    :: fvm(:)
  type (element_t)                      :: elem(:)
  type (hybrid_t)                       :: hybrid
  integer                               :: ie,j,nets,nete
  logical                               :: corner
  

  call compute_ghost_corner_orientation(hybrid,elem,nets,nete)
  ! run some tests:
  !   call test_ghost(hybrid,elem,nets,nete)

  do ie=nets,nete
    fvm(ie)%faceno=elem(ie)%FaceNum
    ! write the neighbors in the structure
    fvm(ie)%cubeboundary=0
    corner=.FALSE.

! Jose Garcia: This code does not work with MESH
! yet we allow it so MESH and fvm can be compiled 
! in the same executable. Some execution paths that
! do not call this routine will work fine.
#ifndef MESH
    do j=1,8
      if (elem(ie)%vertex%nbrs(j)%used) then
        fvm(ie)%nbrsface(j)=elem(ie)%vertex%nbrs(j)%f
        ! note that if the element lies on a corner, it will be at j=5,6,7,8
        if ((fvm(ie)%nbrsface(j) /= fvm(ie)%faceno) .AND. (j<5)) then
          fvm(ie)%cubeboundary=j
        endif
      else   ! corner on the cube
        if (.NOT. corner) then
          fvm(ie)%nbrsface(j)=-1
          fvm(ie)%cubeboundary=j
          corner=.TRUE.
        else
          print *,'Error in fvm_CONTROL_VOLUME_MOD - Subroutine fvm_MESH_ARI: '
          call abortmp('Do not allow one element per face for fvm, please increase ne!')
        endif
      end if
    end do
#endif
  end do
end subroutine fvm_init2

! first communciation of FVM tracers
subroutine fvm_init3(elem,fvm,hybrid,nets,nete,tnp0)
  
  type (element_t),intent(inout)            :: elem(:)                 
  type (cwfv_struct),intent(inout)           :: fvm(:)  
                  
  type (hybrid_t),intent(in)                :: hybrid                  
                                                                            
  integer,intent(in)                        :: nets,nete,tnp0       
                                                                            
  integer                                   :: ie                 
  
  ! do it only for FVM tracers, FIRST TRACER will be the AIR DENSITY   
  do ie=nets,nete 
    call ghostVpack(cellghostbuf, fvm(ie)%c,np,np,nlev,ntrac,0,tnp0,timelevels,elem(ie)%desc)
  end do
  !exchange values for the initial data
  call ghost_exchangeV(hybrid,cellghostbuf,np,np)
  !-----------------------------------------------------------------------------------!
  do ie=nets,nete
    call ghostVunpack(cellghostbuf, fvm(ie)%c, np, np,nlev,ntrac, 0, tnp0, timelevels,elem(ie)%desc)
  enddo

end subroutine fvm_init3
! ----------------------------------------------------------------------------------!
!SUBROUTINE FVM_MESH_DEP--------------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 30.October 2011                                          !
! DESCRIPTION: Calculates the departure mesh                                        !
!              Please select the test cases by an COMMENT/UNCOMMENT basis           !
!                                                                                   !
! CALLS: solidbody, cart2cubedspherexy, analytical_function                         !
! INTPUT/OUTPUT:                                                                    !
!        fvm   ...  structure                                                       !
! INPUT: nstep   ... actual step                                                    !
!-----------------------------------------------------------------------------------!
subroutine fvm_mesh_dep_gll(elem, deriv, fvm, dt, tl, klev)
  use coordinate_systems_mod, only : cartesian2D_t
  use element_mod, only : element_t
  use time_mod, only : timelevel_t, time_at
  use parallel_mod, only : haltmp
  use control_mod, only : test_cfldep

  use derivative_mod, only : derivative_t
  use interpolate_mod, only : cube_facepoint_ne 
  

#ifndef _PRIM
  use fvm_bsp_mod, only: boomerang, solidbody
#endif

  implicit none
  type (derivative_t)  , intent(in) :: deriv
  type (cwfv_struct), intent(inout)   :: fvm
  type (timelevel_t),intent(in)        :: tl
  integer,intent(in)                   :: klev
  
  type (element_t), intent(inout)      :: elem
  real (kind=real_kind)                :: time,dt, dx, dy, cflx, cfly, maxcflx, maxcfly  
  integer                              :: i,j
  
  type (cartesian2D_t)                        :: cart
  integer                                     :: number

! for the benchmark test, use more accurate departure point creation
#if 1
!CE: define new mesh for fvm fvm on an equal angular grid
! go from alpha,beta -> cartesian xy on cube -> lat lon on the sphere
  do j=1,np
     do i=1,np               
!                 call solidbody(fvm%asphere(i,j), fvm%dsphere(i,j))
        call boomerang(elem%spherep(i,j), fvm%dsphere(i,j),tl%nstep)
     end do
  end do
#else

! for given velocities in the element structure
  call fvm_dep_from_gll(elem, deriv, elem%spherep,fvm%dsphere,dt,tl,klev)
#endif

  do j=1,np
    do i=1,np
      call cube_facepoint_ne(fvm%dsphere(i,j),ne,fvm%dcart(i,j), number)
            
        if(number==elem%GlobalId) then          !same element
           fvm%xstart(i,j)=1
           fvm%ystart(i,j)=1
        end if   
        if(number==elem%vertex%nbrs(1)%n) then  !west
           fvm%xstart(i,j)=1-np
           fvm%ystart(i,j)=1
         end if   
           
        if(number==elem%vertex%nbrs(2)%n) then   !east
           fvm%xstart(i,j)=np+1
           fvm%ystart(i,j)=1
         end if   
           
        if(number==elem%vertex%nbrs(3)%n) then   !south
           fvm%xstart(i,j)=1
           fvm%ystart(i,j)=1-np
         end if   
           
        if(number==elem%vertex%nbrs(4)%n) then   !north
           fvm%xstart(i,j)=1
           fvm%ystart(i,j)=np+1
         end if   
           
        if(number==elem%vertex%nbrs(5)%n) then   !southwest
           fvm%xstart(i,j)=1-np
           fvm%ystart(i,j)=1-np
         end if   
           
        if(number==elem%vertex%nbrs(6)%n) then   !southeast
           fvm%xstart(i,j)=np+1
           fvm%ystart(i,j)=1-np
         end if   

        if(number==elem%vertex%nbrs(7)%n) then   !northwest
           fvm%xstart(i,j)=1-np
           fvm%ystart(i,j)=np+1
         end if   
           
        if(number==elem%vertex%nbrs(8)%n) then   !northeast
           fvm%xstart(i,j)=np+1
           fvm%ystart(i,j)=np+1
         end if   
   
    end do
  end do
  
   if (test_cfldep) then
     call check_departurecell(fvm,klev) 
   endif 
  
end subroutine fvm_mesh_dep_gll


! ----------------------------------------------------------------------------------!
!SUBROUTINE FVM_DEP_FROM_GLL----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR 14. December 2011                            !
! DESCRIPTION: calculates the deparute grid for fvm coming from the gll points      !
!                                                                                   !
! CALLS: 
! INPUT: 
!        
! OUTPUT: 
!-----------------------------------------------------------------------------------!
subroutine fvm_dep_from_gll(elem, deriv, asphere,dsphere,dt,tl,klev)
  use physical_constants, only : DD_PI, rearth
  use coordinate_systems_mod, only : spherical_polar_t, cartesian3D_t, change_coordinates
  use time_mod, only : timelevel_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, interpolate_gll2fvm_corners

  implicit none
  type (element_t), intent(in)          :: elem
  type (derivative_t)  , intent(in) :: deriv
  type (spherical_polar_t),intent(in)   :: asphere(np,np)
  type (spherical_polar_t),intent(out)  :: dsphere(np,np)
  type (timelevel_t),intent(in)        :: tl
  real (kind=real_kind),intent(in)      :: dt
  integer,intent(in)                   :: klev
  
  real (kind=real_kind)                 :: uxyz_gll(np,np,3)
  real (kind=real_kind)                 :: un,ue,ur,clon,clat,slon,slat
  type(cartesian3D_t)                   :: acart
  integer :: i,j,k
  
  real (kind=real_kind)  :: vstar(np,np,2)
  
   ! convert velocity from lat/lon to cartesian 3D
   vstar = elem%derived%vstar(:,:,:,klev)
    
   do i=1,np
     do j=1,np
        clon = cos(elem%spherep(i,j)%lon)
        slon = sin(elem%spherep(i,j)%lon)
        clat = cos(elem%spherep(i,j)%lat)
        slat = sin(elem%spherep(i,j)%lat)

        ur = 0
        ue = vstar(i,j,1) 
        un = vstar(i,j,2)

        uxyz_gll(i,j,1)= clon*clat*ur - clon*slat*un - slon*ue
        uxyz_gll(i,j,2)= slon*clon*ur - slon*slat*un + clon*ue
        uxyz_gll(i,j,3)= slat          *ur + clat          *un
    
     enddo
   enddo
   ! compute departure point 
   ! crude, 1st order accurate approximation.  to be improved
   do i=1,np
      do j=1,np
         ! note: asphere()%r=1, so we need to convert velocity to radians/sec:
         acart = change_coordinates(asphere(i,j))  
         acart%x = acart%x - dt*uxyz_gll(i,j,1)/rearth
         acart%y = acart%y - dt*uxyz_gll(i,j,2)/rearth
         acart%z = acart%z - dt*uxyz_gll(i,j,3)/rearth
         dsphere(i,j) = change_coordinates(acart)
         dsphere(i,j)%r = asphere(i,j)%r
      enddo
   enddo
  
end subroutine fvm_dep_from_gll
!END SUBROUTINE FVM_DEP_FROM_GLL------------------------------------------CE-for FVM!

!SUBROUTINE CHECK_DEPARTURECELL-------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: Check if the departure cell is degenerated or if the                 !
!              CFL number > 1.0D0                                                   !
!         IF THE CHECK FAILS, THE PROGRAM WILL STOP                                 !
!                                                                                   !
! INPUT:  fvm  ... fvm structur                                                 ! 
!         klev   ... Level (vertical)                                               !
!-----------------------------------------------------------------------------------!
subroutine check_departurecell(fvm,klev)
  use coordinate_systems_mod, only: cartesian2D_t, cart2cubedspherexy, spherical_to_cart
  implicit none

  type (cwfv_struct), intent(inout)   :: fvm
  integer,intent(in)                   :: klev
  
  type (cartesian2D_t)                 :: dcart(np,np)
  real (kind=real_kind)                :: cflx, cfly, maxcflx, maxcfly  
  integer                              :: i,j
  logical                              :: crossline
            
  ! calculate xy Cartesian on the cube of departure points on the corresponding face  
  ! the calculation of dcart is also done in fvm_line_integrals_mod.F90, one could 
  ! save the points here...
  ! BUT: this subroutine should only be used for test reasons
  do j=1,np
     do i=1,np               
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(i,j)),&
             fvm%faceno,dcart(i,j))              
     end do
  end do                        
  cflx=0.0D0
  cfly=0.0D0
  fvm%maxcfl(1,klev)=cflx
  fvm%maxcfl(2,klev)=cfly
  crossline=.FALSE.
  do j=1,np-1
     do i=1,np-1  
      ! equidistant mesh in alpha/beta coordinates
!       cflx=abs((atan(dcart(i,j)%x)-atan(fvm%acartx(i)))/fvm%dalpha)
!       cfly=abs((atan(dcart(i,j)%y)-atan(fvm%acarty(j)))/fvm%dbeta)
!       if(cflx>fvm%maxcfl(1,klev)) fvm%maxcfl(1,klev)=cflx
!       if(cfly>fvm%maxcfl(2,klev)) fvm%maxcfl(2,klev)=cfly
    
      !one could stop immediately here, if crossline=.TRUE., but want to calculate 
      ! all CFL first
      call check_lines_cross(dcart(i,j),dcart(i+1,j),dcart(i,j+1),dcart(i+1,j+1),crossline)  
      call check_lines_cross(dcart(i,j),dcart(i,j+1),dcart(i+1,j),dcart(i+1,j+1),crossline)
     end do
  end do  
  ! nodes on the north and east boundary
!   do i=1,nc+1  
!     cflx=abs((atan(dcart(i,nc+1)%x)-atan(fvm%acartx(i)))/fvm%dalpha)
!     cfly=abs((atan(dcart(nc+1,j)%y)-atan(fvm%acarty(j)))/fvm%dbeta)
!     if(cflx>fvm%maxcfl(1,klev)) fvm%maxcfl(1,klev)=cflx
!     if(cfly>fvm%maxcfl(2,klev)) fvm%maxcfl(2,klev)=cfly
!   end do
! 
!   if (fvm%maxcfl(1,klev) > 1.0D0 .OR. fvm%maxcfl(2,klev) > 1.0D0) then
!     write(*,*) "Error in fvm_mod.F90: CFL number too big"
!     write(*,*) "CFL has to be < 1.0D0"
!     write(*,*) "Choose a smaller time step!"
!     write(*,*) "max CFL in this element: maxcflx", fvm%maxcfl(1,klev), "maxcfly", fvm%maxcfl(2,klev) 
!     STOP "Exit program!"
!   endif

  if (crossline) then
    write(*,*) "FATAL Error in fvm_mod.F90: departure cell is degenerated!"
    write(*,*) "Choose a smaller time step! max CFL in this element: maxcflx", fvm%maxcfl(1,klev), "maxcfly", fvm%maxcfl(2,klev) 
    STOP "Exit program!"
  endif
end subroutine check_departurecell


! see, e.g., http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
subroutine check_lines_cross(p1,p2,q1,q2,crossline)
  use coordinate_systems_mod, only: cartesian2D_t
  implicit none
  type(cartesian2D_t), intent (in)  :: p1,p2,q1,q2
  logical            , intent(inout)  :: crossline
  !
  ! local workspace
  !
  REAL (kind=real_kind)    :: dn,tp,tq
  
  dn = (q2%y-q1%y)*(p2%x-p1%x)-(q2%x-q1%x)*(p2%y-p1%y)

  if (abs(dn)>1.0D-12) then
    tp = ((q2%x-q1%x)*(p1%y-q1%y)-(q2%y-q1%y)*(p1%x-q1%x))/dn
    tq = ((p2%x-p1%x)*(p1%y-q1%y)-(p2%y-p1%y)*(p1%x-q1%x))/dn
    ! implement the next two lines faster!!!!!!
    if (tp>=0.0D0 .and. tp<=1.0D00 .and. &
        tq>=0.0D0 .and. tq<=1.0D00) then
        crossline=.TRUE.
    endif
  endif
end subroutine check_lines_cross


! ----------------------------------------------------------------------------------!
!SUBROUTINE FVM_MCGREGOR--------------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 09. January 2012                                         !
! DESCRIPTION: ! using McGregor AMS 1993 scheme: Economical Determination of        !
!                Departure Points for Semi-Lagrangian Models                        !
!                                                                                   !
! CALLS: 
! INPUT: 
!        
! OUTPUT: 
!-----------------------------------------------------------------------------------!
subroutine fvm_mcgregor(elem, deriv, tstep, vhat, vstar,order)
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, gradient_sphere, ugradv_sphere, vorticity_sphere
  implicit none

  type (element_t), intent(in)                                :: elem
  type (derivative_t), intent(in)                             :: deriv      ! derivative struct
  real (kind=real_kind), intent(in)                           :: tstep
  real (kind=real_kind), dimension(np,np,2), intent(inout)    :: vstar
  real (kind=real_kind), dimension(np,np,2), intent(in)       :: vhat
  
  integer, intent(in)                                         :: order

  integer                                            :: i
  real (kind=real_kind), dimension(np,np,2)          :: ugradv
  real (kind=real_kind)                              :: timetaylor
    
  ugradv=vstar
  timetaylor=1
  do i=1,order
!     tmp = 0.5D0*(vgradv(:,:,1)**2 + vgradv(:,:,2)**2) 
! 
!     gradvstar = gradient_sphere(tmp,deriv,elem%Dinv)    ! scalar -> latlon vector
!     tmp = vorticity_sphere(vgradv,deriv,elem)                 ! latlon vector -> scalar 
! 
!     ! v \nabla v expressed through gradient of mean scalar and relative velocity
!     ! see (17) in Williamson et.al. 1992, JCP 102 (211-224)
!     vgradv(:,:,1)= gradvstar(:,:,1) - vstarold(:,:,2)*tmp(:,:)
!     vgradv(:,:,2)= gradvstar(:,:,2) + vstarold(:,:,1)*tmp(:,:)
!     
!     timetaylor=-timetaylor*tstep/(i+1)
!     vstar=vstar + timetaylor*vgradv
    
    ugradv=ugradv_sphere(vhat,ugradv,deriv,elem)
    timetaylor=-timetaylor*tstep/(i+1)
    
    vstar=vstar + timetaylor*ugradv  
  end do
end subroutine fvm_mcgregor
!END SUBROUTINE FVM_MCGREGOR----------------------------------------------CE-for FVM!
! ----------------------------------------------------------------------------------!
!SUBROUTINE FVM_MCGREGORDSS-----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 26. May 2012                                             !
! DESCRIPTION: ! using McGregor AMS 1993 scheme: Economical Determination of        !
!                Departure Points for Semi-Lagrangian Models                        !
!                McGegror version with DSS every ugradv                             !
! CALLS: 
! INPUT: 
!        
! OUTPUT: 
!-----------------------------------------------------------------------------------!
subroutine fvm_mcgregordss(elem,fvm,nets,nete, hybrid, deriv, tstep, ordertaylor)
  use derivative_mod, only : derivative_t, ugradv_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  
  implicit none

  type (element_t), intent(inout)                :: elem(:)
  type (cwfv_struct), intent(in)                :: fvm(:)

  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)

  type (derivative_t), intent(in)                             :: deriv      ! derivative struct
  real (kind=real_kind), intent(in)                           :: tstep
  integer, intent(in)                                         :: ordertaylor

  real (kind=real_kind), dimension(nets:nete,np,np,2,nlev)    :: ugradv
  real (kind=real_kind), dimension(nets:nete,np,np,2,nlev)    :: vhat
  integer                                                     :: ie, k, order
  real (kind=real_kind), dimension(np,np,2)                   :: ugradvtmp
  real (kind=real_kind)                                       :: timetaylor

    !------------------------------------------------------------------------------------
  timetaylor=1  
  do  order=1,ordertaylor
    timetaylor=-timetaylor*tstep/(order+1)  
    do ie=nets,nete
      if (order==1)then
        ugradv(ie,:,:,:,:)=elem(ie)%derived%vstar(:,:,:,:) 
        vhat(ie,:,:,:,:)=(fvm(ie)%vn0(:,:,:,:) + ugradv(ie,:,:,:,:))/2 
      endif
      do k=1,nlev
        ugradvtmp=ugradv_sphere(vhat(ie,:,:,:,k),ugradv(ie,:,:,:,k),deriv,elem(ie))
        ugradv(ie,:,:,1,k) = elem(ie)%spheremp(:,:)*ugradvtmp(:,:,1) 
        ugradv(ie,:,:,2,k) = elem(ie)%spheremp(:,:)*ugradvtmp(:,:,2) 
      enddo 
      call edgeVpack(edgeveloc,ugradv(ie,:,:,1,:),nlev,0,elem(ie)%desc)
      call edgeVpack(edgeveloc,ugradv(ie,:,:,2,:),nlev,nlev,elem(ie)%desc)
    enddo 
    call bndry_exchangeV(hybrid,edgeveloc)
    do ie=nets,nete
       call edgeVunpack(edgeveloc,ugradv(ie,:,:,1,:),nlev,0,elem(ie)%desc)
       call edgeVunpack(edgeveloc,ugradv(ie,:,:,2,:),nlev,nlev,elem(ie)%desc)
       do k=1, nlev  
         ugradv(ie,:,:,1,k)=ugradv(ie,:,:,1,k)*elem(ie)%rspheremp(:,:)
         ugradv(ie,:,:,2,k)=ugradv(ie,:,:,2,k)*elem(ie)%rspheremp(:,:)
         elem(ie)%derived%vstar(:,:,:,k)=elem(ie)%derived%vstar(:,:,:,k) + timetaylor*ugradv(ie,:,:,:,k)
       end do
    end do
  end do  

end subroutine fvm_mcgregordss
!END SUBROUTINE FVM_MCGREGORDSS-------------------------------------------CE-for FVM!

end module cwfv_mod
