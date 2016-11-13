#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module unit_tests_mod

implicit none

public ::  test_ibyp
public ::  test_subcell_dss_fluxes
public ::  test_subcell_div_fluxes
public ::  test_subcell_Laplace_fluxes
public ::  test_subcell_div_fluxes_again
public ::  test_subcell_dss_fluxes_again
public ::  test_subcell_Laplace_fluxes_again
public ::  test_sub_integration
public ::  test_edge_flux

contains


  subroutine test_ibyp(elem, hybrid,  nets,   nete)
!
! Note: vector test functions should be co-variant since u is contra-variant
!  PHIvec = PHIcov  (test function)
!  PHIcon = DtD PHIcov
!
! weak grad:
!  < PHIcov du/dt > = < PHIcon grad(p) >    (output of grad is covariant)
!  < PHIcov du/dt > = -< div(PHIcon) p >    (input of div is contra)
!  verify:
!    gradient_sphere_wk(p) = - <div(PHIcon) p >
!    gradient_sphere_wk(p) + MASS*grad(p) = b.c. (b.c. are covariant)
!
! weak div:
!   < PHI div(u) > = -< grad(PHI) dot u >     u=contra, output of grad is covariant
! verify:
!   divergence_sphere_wk(u) = -<grad(PHI) dot u>
!   divergence_sphere_wk(u) + MASS*div(u) = b.c.  (b.c. are scalars)
!
! weak curl:
!  < PHIcov du/dt > = < PHIcov curl( a ) >    (output of curl is contra)
!  < PHIcov du/dt > = < vor(PHIcov) a >       (input to vor is covariant)
! verify:
!    curl_sphere_wk(a) = < vor(PHIcov) a >
!    curl_sphere_wk(a) - MASS*curl(a) = b.c. (b.c. are contra)
!
    ! ---------------------
    use kinds, only : real_kind, iulog
    ! ---------------------
    use physical_constants, only : rearth 
    ! ---------------------
    use dimensions_mod, only : np, nlev
    ! ---------------------
    use element_mod, only : element_t
    ! ---------------------
    use hybrid_mod, only : hybrid_t
    ! ---------------------
    use derivative_mod, only : derivative_t, gradient_sphere, divergence_sphere,vorticity_sphere,&
                               divergence_sphere_wk, curl_sphere, derivinit
    use viscosity_mod, only : make_c0, make_c0_vector
    use global_norms_mod
    use coordinate_systems_mod, only : cartesian3D_t, spherical_to_cart

    implicit none

    type (element_t)     , intent(inout), target :: elem(:)

    type (hybrid_t)      , intent(in) :: hybrid

    integer              , intent(in) :: nets
    integer              , intent(in) :: nete

#undef CURLGRAD_TEST
#define IBYP_TEST
    ! =================
    ! Local
    ! =================
    ! pointer ...
    real (kind=real_kind), dimension(:,:), pointer :: rspheremv,spheremv

    ! Thread private working set ...

    real (kind=real_kind), dimension(np,np,nets:nete) :: ptens
    real (kind=real_kind), dimension(np,np,nets:nete) :: ptens2
    real (kind=real_kind), dimension(np,np,nets:nete) :: ptens3

    real (kind=real_kind), dimension(np,np,2,nets:nete)    :: pv      ! p*v lat-lon
    real (kind=real_kind), dimension(np,np,nets:nete)            :: E          ! kinetic energy term
    real (kind=real_kind), dimension(np,np,nets:nete)  :: divbig
    real (kind=real_kind), dimension(np,np,nets:nete)  :: wdivbig
    real (kind=real_kind), dimension(np,np,2,nets:nete)    :: gradbig

    real (kind=real_kind), dimension(np,np,2)      :: grade
    real (kind=real_kind), dimension(np,np,2)      :: grade2
    real (kind=real_kind), dimension(np,np)      :: vor
    real (kind=real_kind), dimension(np,np)      :: div  
    real (kind=real_kind), dimension(np,np)      :: wdiv  

    real (kind=real_kind) ::  v1,v2,v3,mx

    real*8                :: st,et, time_adv
    integer    :: i,j,k,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    type (cartesian3D_t)             :: cart
    type (derivative_t)          :: deriv
    call derivinit(deriv)

    ! ===================================
    ! construct test functions for pv and E
    ! ===================================
    do ie=nets,nete
       do j=1,np
          do i=1,np
             cart = spherical_to_cart(elem(ie)%spherep(i,j))
             E(i,j,ie)=cart%x**2 + cart%y + cart%z
             pv(i,j,1,ie) = cart%x + cart%y**2 + cart%z
             pv(i,j,2,ie) = cart%x + cart%y + cart%z**2
          end do
       end do
    enddo
    call make_C0(E,elem,hybrid,nets,nete)
    call make_C0_vector(pv,elem,hybrid,nets,nete) 

#ifdef CURLGRAD_TEST
    ! check curl(grad(E)) 
    do ie=nets,nete
       if ( maxval(abs(E(:,:,ie))) > 0 ) then
       !write(iulog,'(a,i4,2e20.10)') 'maxval: E =',ie,maxval(E(:,:,ie))
       grade=curl_sphere(E(:,:,ie),deriv,elem(ie))
       div=divergence_sphere(grade,deriv,elem(ie))
       vor=vorticity_sphere(grade,deriv,elem(ie))
       if (maxval(abs(div))*rearth**2 > .2e-11) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: div(curl),  vor(curl)=',ie,maxval(abs(div))*rearth**2,maxval(abs(vor))*rearth**2
       endif

       grade=gradient_sphere(E(:,:,ie),deriv,elem(ie)%Dinv)
       vor=vorticity_sphere(grade,deriv,elem(ie))
       div=divergence_sphere(grade,deriv,elem(ie))
       if (maxval(abs(vor))*rearth**2 > .2e-11) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: curl(grad), div(grad)=',ie,maxval(abs(vor))*rearth**2,maxval(abs(div))*rearth**2
       endif
       endif
    enddo

    ! check div(curl(E)) with DSS 
    do ie=nets,nete
       gradbig(:,:,:,ie)=curl_sphere(E(:,:,ie),deriv,elem(ie))
    enddo
    call make_C0_vector(gradbig,elem,hybrid,nets,nete)
    do ie=nets,nete
       divbig(:,:,ie)=divergence_sphere(gradbig(:,:,:,ie),deriv,elem(ie))
    enddo
    call make_C0(divbig,elem,hybrid,nets,nete)
    do ie=nets,nete
       if (maxval(abs(divbig(:,:,ie)))*rearth**2 > .8e-12) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: [div([curl])]=',ie,maxval(abs(divbig(:,:,ie)))*rearth**2
       endif
    enddo


    ! check curl(grad(E)) with DSS 
    do ie=nets,nete
       gradbig(:,:,:,ie)=gradient_sphere(E(:,:,ie),deriv,elem(ie)%Dinv)
    enddo
    call make_C0_vector(gradbig,elem,hybrid,nets,nete)
    do ie=nets,nete
       divbig(:,:,ie)=vorticity_sphere(gradbig(:,:,:,ie),deriv,elem(ie))
    enddo
    call make_C0(divbig,elem,hybrid,nets,nete)
    
    do ie=nets,nete
       if (maxval(abs(divbig(:,:,ie)))*rearth**2 > .8e-12) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: [curl([gradl])]=',ie,maxval(abs(divbig(:,:,ie)))*rearth**2
       endif
    enddo
#endif


#ifdef IBYP_TEST
    ! compare <grad(E) dot pv> and <E div(pv)>  < E weak_div(pv) >
    v2=0
    do ie=nets,nete
       spheremv     => elem(ie)%spheremp(:,:)

       div = divergence_sphere(pv(1,1,1,ie),deriv,elem(ie))      ! latlon vector -> scalar 
       grade = gradient_sphere(E(1,1,ie),deriv,elem(ie)%Dinv)
       wdiv = divergence_sphere_wk(pv(1,1,1,ie),deriv,elem(ie)) 

       do j=1,np
          do i=1,np
!             write(iulog,'(3i3,3e22.14)') ie,i,j,pv(i,j,1,ie),pv(i,j,2,ie),div(i,j)
             ! (grad(E) dot pv )
             ptens3(i,j,ie) = grade(i,j,1)*pv(i,j,1,ie) + grade(i,j,2)*pv(i,j,2,ie)
             v2=v2+wdiv(i,j)*E(i,j,ie)
          end do
       end do
       ptens(:,:,ie)=div(:,:)*E(:,:,ie)   ! < E div(pv) >
       ptens2(:,:,ie)=wdiv(:,:)*E(:,:,ie)/spheremv(:,:)   ! < wdiv E >
       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
!       divbig(:,:,ie)=div(:,:)*spheremv(:,:)
!       wdivbig(:,:,ie)=wdiv(:,:)
       divbig(:,:,ie)=div(:,:)
       wdivbig(:,:,ie)=wdiv(:,:)/spheremv(:,:)
    end do
    call make_C0(divbig,elem,hybrid,nets,nete)
    call make_C0(wdivbig,elem,hybrid,nets,nete)

    v1=global_integral(elem,ptens,hybrid,np,nets,nete)                                                 
    v2=global_integral(elem,ptens2,hybrid,np,nets,nete)                                                
    v3=global_integral(elem,ptens3,hybrid,np,nets,nete)                                                
    mx =max(abs(v1),abs(v2),abs(v3))
    if (hybrid%masterthread) then   
       print *,'< E div(pv) >   =',v1                                                                     
       print *,'< E div_wk(pv) >=',v2                                                                     
       print *,'-<grad(E),pv >  =',-v3                                                                    
       print *,'integration by parts rel error:',(v1-v2)/mx,(v1+v3)/mx
    endif
    if (( abs(v1-v2)/mx .gt. 1e-12 ) .or. ( abs(v1+v3)/mx .gt. 1e-12 )) then
       stop 'integration by parts error1 too large?'
    endif
  
#if 0  
    do ie=nets,nete
       div(:,:)=divbig(:,:,ie)
       wdiv(:,:)=wdivbig(:,:,ie)
       do j=1,np
       do i=1,np
          ! < div(pv) >   vs < div_wk(pv) >
          if ( abs(div(i,j)-wdiv(i,j)) > .15e-17) then
             write(iulog,'(3i3,4e22.14)') ie,i,j,div(i,j),wdiv(i,j),div(i,j)-wdiv(i,j),E(i,j,ie)
          endif
       end do
       end do
    end do
#endif
    ! test function is O(1). gradient is O(1/rearth).  lets require agreement to
    ! 1e-11/rearth
    mx = rearth*maxval(abs(divbig(:,:,:)-wdivbig(:,:,:)))
    if (hybrid%masterthread) then   
       write(iulog,'(a,2e20.10)') 'div vs. weak div, max error on masterthread: ',mx
    endif
    if (mx >  1e-11 ) then
       write(iulog,'(a,2e20.10)') 'max diff div-wdiv: ',mx,maxval(divbig(:,:,:))
       stop 'integration by parts error2 too large?'
    endif
#endif
  end subroutine test_ibyp


  subroutine test_subcell_dss_fluxes(elem,deriv,nets,nete)
    use dimensions_mod, only : np
    use derivative_mod, only : subcell_dss_fluxes
    use derivative_mod, only : subcell_integration
    use element_mod,    only : element_t
    use derivative_mod, only : derivative_t
    use kinds,          only : real_kind

    implicit none

    type (element_t)     , intent(in) :: elem(:)
    type (derivative_t)  , intent(in) :: deriv
    integer              , intent(in) :: nets,nete

    integer              , parameter :: intervals=6 

    real (kind=real_kind)              :: dss(np,np)
    real (kind=real_kind)              :: values(intervals,intervals)
    real (kind=real_kind)              :: fluxes(intervals,intervals,4)
    real (kind=real_kind)              :: test(intervals,intervals)
    real (kind=real_kind)              :: cflux(2,2,2)
    real (kind=real_kind)              :: p,t

    integer                            :: ie,i,j,k
    logical                            :: success
    success = .true.


    do ie=nets,nete
      call random_number(p)

      if (ie <= np*np) then
        dss = 0 
        dss(1+mod(ie,np), 1+mod(ie/np,np)) = 1
      else
        t = dss(1+mod((7*ie),np), 1+mod((13*ie)/np,np))
        dss(1+mod(ie,np), 1+mod(ie/np,np)) = t + 10*p
        dss(1+mod(INT(32147*p),np), 1+mod(INT(1123*p)/np,np)) = 0
      end if

      dss(2:np-1,2:np-1) = 0

      do i=1,2
      do j=1,2
        do k=1,2
          call random_number(cflux(i,j,k))
        end do
        cflux(i,j,1) = cflux(i,j,1)*dss(np*i-np+2-i,np*j-np+2-j)+cflux(i,j,2)
        cflux(i,j,2) = dss(np*i-np+2-i,np*j-np+2-j)-cflux(i,j,1)
      end do
      end do

      values = subcell_integration(dss, np, intervals, elem(ie)%metdet) 
      fluxes = subcell_dss_fluxes (dss, np, intervals, elem(ie)%metdet, cflux) 

      test = SUM(fluxes,3)

      do i=1,intervals
      do j=1,intervals
      if (.00001<ABS(test(i,j)-values(i,j))) then
        print *,__FILE__,__LINE__,ie,i,j,test(i,j),values(i,j)
        success = .false.
      end if
      end do
      end do
    end do

    if (success) then
      print *,__FILE__,__LINE__," test_subcell_dss_fluxes test passed."
    else
      print *,__FILE__,__LINE__," test_subcell_dss_fluxes test FAILED."
    end if

  end subroutine test_subcell_dss_fluxes

  subroutine test_subcell_div_fluxes(elem,deriv,nets,nete)
    use dimensions_mod, only : np
    use derivative_mod, only : subcell_div_fluxes
    use derivative_mod, only : subcell_integration
    use derivative_mod, only : divergence_sphere
    use element_mod,    only : element_t
    use derivative_mod, only : derivative_t
    use kinds,          only : real_kind
    use coordinate_systems_mod, only: spherical_polar_t


    implicit none

    type (element_t)     , intent(in) :: elem(:)
    type (derivative_t)  , intent(in) :: deriv
    integer              , intent(in) :: nets,nete

    integer              , parameter :: intervals=6 

    real (kind=real_kind)              :: u(np,np,2), v(np,np,2)
    real (kind=real_kind)              :: div(np,np)
    real (kind=real_kind)              :: values(intervals,intervals)
    real (kind=real_kind)              :: fluxes(intervals,intervals,4)
    real (kind=real_kind)              :: test(intervals,intervals)
    real (kind=real_kind)              :: p,t

    type(spherical_polar_t)            :: s
    integer                            :: ie,i,j
    logical                            :: success
    success = .true.


    do ie=nets,nete

      call random_number(p)

      if (ie <= np*np) then
        u = 0 
        u(1+mod(ie,np), 1+mod(ie/np,np),1) = 1
      else if (ie <= 2*np*np) then
        u = 0 
        u(1+mod(ie,np), 1+mod(ie/np,np),2) = 1
      else
        t = u(1+mod((7*ie),np), 1+mod((13*ie)/np,np),1)
        u(1+mod(ie,np), 1+mod(ie/np,np),1) = t + 10*p
        u(1+mod(ie,np), 1+mod(ie/np,np),2) = t/2 + 15*p
        u(1+mod(INT(32147*p),np), 1+mod(INT(1123*p)/np,np),1) = 0
        u(1+mod(INT(23583*p),np), 1+mod(INT(1317*p)/np,np),2) = 0
      end if

      div  = divergence_sphere(u, deriv, elem(ie))
   
      v(:,:,1) = elem(ie)%Dinv(:,:,1,1)*u(:,:,1) + elem(ie)%Dinv(:,:,1,2)*u(:,:,2)
      v(:,:,2) = elem(ie)%Dinv(:,:,2,1)*u(:,:,1) + elem(ie)%Dinv(:,:,2,2)*u(:,:,2)

      values = subcell_integration(div, np, intervals, elem(ie)%metdet) 
      fluxes = subcell_div_fluxes (v,   np, intervals, elem(ie)%metdet) 

      test = SUM(fluxes,3)

      do i=1,intervals
      do j=1,intervals
!        t = ABS(test(i,j)-values(i,j))/MAX(ABS(test(i,j)),ABS(values(i,j)))
        t = ABS(test(i,j)-values(i,j))/MAX(ABS(MAXVAL(test(:,:))),ABS(MAXVAL(values(:,:))))
        if (.0000001<t) then
          print *,__FILE__,__LINE__,ie,i,j,test(i,j),values(i,j),t
          success = .false.
        end if
      end do
      end do
    end do

    if (success) then
      print *,__FILE__,__LINE__," test_subcell_div_fluxes test passed."
    else
      print *,__FILE__,__LINE__," test_subcell_div_fluxes test FAILED."
    end if

  end subroutine test_subcell_div_fluxes


  subroutine test_subcell_div_fluxes_again(elem,deriv,nets,nete)
    use physical_constants, only : rearth
    use dimensions_mod, only : np
    use derivative_mod, only : subcell_div_fluxes
    use element_mod,    only : element_t
    use derivative_mod, only : derivative_t
    use kinds,          only : real_kind
    use quadrature_mod, only : gausslobatto, quadrature_t

    implicit none

    type (element_t)     , intent(in) :: elem(:)
    type (derivative_t)  , intent(in) :: deriv
    integer              , intent(in) :: nets,nete

    integer              , parameter :: intervals=6 

    real (kind=real_kind)              :: v(np,np,2)
    real (kind=real_kind)              :: fluxes(intervals,intervals,4)
    real (kind=real_kind)              :: t

    real (kind=real_kind)              :: metdet(np,np)
    type (quadrature_t)                :: gll
    integer                            :: ie,i,j
    logical                            :: success
    real (kind=real_kind),   parameter :: EPS=.0000001
    success = .true.

    gll = gausslobatto(np)

    metdet = 1
    do ie=nets,nete

      v = 0 
      if (ie <= np*np) then
        do i = 1,np
          v(i,:,2) = gll%points(:)
        end do
      else if (ie <= 2*np*np) then
        do j = 1,np
          v(:,j,1) = gll%points(:)
        end do
      else if (ie <= 3*np*np) then
        do i = 1,np
          v(i,:,2) = 1+gll%points(:)
        end do
      else if (ie <= 4*np*np) then
        do j = 1,np
          v(:,j,1) = 1+gll%points(:)
        end do
      end if

      fluxes = subcell_div_fluxes(v, np, intervals, metdet) 
      fluxes = rearth*fluxes

      if (ie <= np*np) then
        t = 4./(intervals*intervals)
        do i=1,intervals
        do j=1,intervals
        ! check for fluxes from the bottom to the top
          if (EPS < ABS(fluxes(i,j,1)+fluxes(i,j,3)-t)) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
            success = .false.
          end if
          if (EPS < MAX(ABS(fluxes(i,j,2)),ABS(fluxes(i,j,4)))) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
            success = .false.
          end if
          if (1.lt.i) then
            if (EPS < ABS(fluxes(i-1,j,1)-fluxes(i,j,1))) then
              print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
              success = .false.
            end if
            if (EPS < ABS(fluxes(i-1,j,3)-fluxes(i,j,3))) then
              print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
              success = .false.
            end if
          end if
        end do
        end do
      else if (ie <= 2*np*np) then
        t = 4./(intervals*intervals)
        do i=1,intervals
        do j=1,intervals
        ! check for fluxes from the left to the right
          if (EPS < ABS(fluxes(i,j,2)+fluxes(i,j,4)-t)) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
            success = .false.
          end if
          if (EPS < MAX(ABS(fluxes(i,j,1)),ABS(fluxes(i,j,3)))) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
            success = .false.
          end if
          if (1.lt.j) then
            if (EPS < ABS(fluxes(i,j-1,2)-fluxes(i,j,2))) then
              print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
              success = .false.
            end if
            if (EPS < ABS(fluxes(i,j-1,4)-fluxes(i,j,4))) then
              print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
              success = .false.
            end if
          end if
        end do
        end do
      else if (ie <= 3*np*np) then
        t = 4./(intervals*intervals)
        do i=1,intervals
        do j=1,intervals
        ! check for fluxes from the bottom to the top
          if (EPS < fluxes(i,j,1)) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:)
            success = .false.
          end if
          if (EPS < ABS(fluxes(i,j,1)+fluxes(i,j,3)-t)) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
            success = .false.
          end if
          if (EPS < MAX(ABS(fluxes(i,j,2)),ABS(fluxes(i,j,4)))) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
            success = .false.
          end if
        end do
        end do
      else if (ie <= 4*np*np) then
        t = 4./(intervals*intervals)
        do i=1,intervals
        do j=1,intervals
        ! check for fluxes from the left to the right
          if (EPS < fluxes(i,j,4)) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:)
            success = .false.
          end if
          if (EPS < ABS(fluxes(i,j,2)+fluxes(i,j,4)-t)) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
            success = .false.
          end if
          if (EPS < MAX(ABS(fluxes(i,j,1)),ABS(fluxes(i,j,3)))) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
            success = .false.
          end if
        end do
        end do
      endif

      if (.not.success) then
        print *,__FILE__,__LINE__," test_subcell_div_fluxes_again test FAILED."
      end if
    end do

    if (success) then
      print *,__FILE__,__LINE__," test_subcell_div_fluxes_again test passed."
    else
      print *,__FILE__,__LINE__," test_subcell_div_fluxes_again test FAILED."
    end if

  end subroutine test_subcell_div_fluxes_again

  subroutine test_subcell_dss_fluxes_again(elem,deriv,nets,nete)
    use physical_constants, only : rearth
    use dimensions_mod, only : np
    use derivative_mod, only : subcell_dss_fluxes
    use element_mod,    only : element_t
    use derivative_mod, only : derivative_t
    use kinds,          only : real_kind
    use quadrature_mod, only : gausslobatto, quadrature_t

    implicit none

    type (element_t)     , intent(in) :: elem(:)
    type (derivative_t)  , intent(in) :: deriv
    integer              , intent(in) :: nets,nete

    integer              , parameter :: intervals=6 

    real (kind=real_kind)              :: dss(np,np)
    real (kind=real_kind)              :: fluxes(intervals,intervals,4)
    real (kind=real_kind)              :: t

    type (quadrature_t)                :: gll
    real (kind=real_kind)              :: metdet(np,np)
    real (kind=real_kind)              :: cflux(2,2,2)
    integer                            :: ie,i,j
    logical                            :: success
    real (kind=real_kind),   parameter :: EPS=.0000001
    success = .true.

    gll = gausslobatto(np)

    metdet = 1
    do ie=nets,nete

      dss = 0 
      if (ie <= np*np) then
        do i=2,np-1
          dss(i,1) = 1
        end do
      else if (ie <= 2*np*np) then
        do j=2,np-1
          dss(1,j) = 1
        end do
      end if

      do i=1,2
      do j=1,2
        cflux(i,j,1) = dss(np*i-np+2-i,np*j-np+2-j)/2
        cflux(i,j,2) = dss(np*i-np+2-i,np*j-np+2-j)/2
      end do
      end do

      fluxes = subcell_dss_fluxes (dss, np, intervals, metdet, cflux)

      if (ie <= np*np) then
        do i=1,intervals
        do j=1,intervals
        ! check for fluxes from the bottom to the top
          if (ABS(fluxes(i,j,1)).lt.ABS(fluxes(i,j,3))) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:)
            success = .false.
          end if
          if (EPS < MAX(ABS(fluxes(i,j,2)),ABS(fluxes(i,j,4)))) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:)
            success = .false.
          end if
          if (intervals.eq.j) then
            if (EPS < ABS(fluxes(i,j,3))) then
              print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:)
              success = .false.
            end if
          end if
        end do
        end do
      else if (ie <= 2*np*np) then
        do i=1,intervals
        do j=1,intervals
        ! check for fluxes from the left to the right
          if (ABS(fluxes(i,j,4)).lt.ABS(fluxes(i,j,2))) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:)
            success = .false.
          end if
          if (EPS < MAX(ABS(fluxes(i,j,1)),ABS(fluxes(i,j,3)))) then
            print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
            success = .false.
          end if
          if (intervals.eq.i) then
            if (EPS < ABS(fluxes(i,j,2))) then
              print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:),t
              success = .false.
            end if
          end if
        end do
        end do
      endif

      if (.not.success) then
        print *,__FILE__,__LINE__," test_subcell_dss_fluxes_again test FAILED."
      end if
    end do

    if (success) then
      print *,__FILE__,__LINE__," test_subcell_dss_fluxes_again test passed."
    else
      print *,__FILE__,__LINE__," test_subcell_dss_fluxes_again test FAILED."
    end if

  end subroutine test_subcell_dss_fluxes_again


  subroutine test_subcell_Laplace_fluxes_again(elem,deriv,nets,nete)
    use physical_constants, only : rearth
    use dimensions_mod, only : np
    use derivative_mod, only : subcell_Laplace_fluxes
    use element_mod,    only : element_t
    use derivative_mod, only : derivative_t
    use kinds,          only : real_kind
    use quadrature_mod, only : gausslobatto, quadrature_t

    implicit none

    type (element_t)     , intent(in) :: elem(:)
    type (derivative_t)  , intent(in) :: deriv
    integer              , intent(in) :: nets,nete

    integer              , parameter :: intervals=6 

    real (kind=real_kind)              :: laplace(np,np)
    real (kind=real_kind)              :: fluxes(intervals,intervals,4)
    real (kind=real_kind)              :: t

    type (quadrature_t)                :: gll
    integer                            :: ie,i,j
    logical                            :: success
    real (kind=real_kind),   parameter :: EPS=.0000001
    success = .true.

    gll = gausslobatto(np)

    do ie=nets,nete

      laplace = 0 
      if (ie <= np*np) then
        do i=1,np
          laplace(i,1) = 1
        end do
      else if (ie <= 2*np*np) then
        do j=2,np-1
          laplace(1,j) = 1
        end do
      end if

      fluxes = subcell_Laplace_fluxes (laplace, deriv, elem(ie), np, intervals)

      do i=1,intervals
        j = 1
        if (fluxes(i,j,1).ne.0) then
          print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:)
          success = .false.
        end if
        j = intervals
        if (fluxes(i,j,3).ne.0) then
          print *,__FILE__,__LINE__,ie,i,j,fluxes(i,j,:)
          success = .false.
        end if

        j = 1
        if (fluxes(j,i,4).ne.0) then
          print *,__FILE__,__LINE__,ie,j,i,fluxes(j,i,:)
          success = .false.
        end if
        j = intervals
        if (fluxes(j,i,2).ne.0) then
          print *,__FILE__,__LINE__,ie,j,i,fluxes(j,i,:)
          success = .false.
        end if
      end do
    end do

    if (success) then
      print *,__FILE__,__LINE__," test_subcell_laplace_fluxes_again test passed."
    else
      print *,__FILE__,__LINE__," test_subcell_laplace_fluxes_again test FAILED."
    end if

  end subroutine test_subcell_Laplace_fluxes_again

  subroutine test_subcell_Laplace_fluxes(elem,deriv,nets,nete)
    use physical_constants, only : rearth
    use dimensions_mod, only : np
    use derivative_mod, only : subcell_Laplace_fluxes
    use derivative_mod, only : subcell_integration
    use derivative_mod, only : laplace_sphere_wk, gradient_sphere
    use derivative_mod, only : divergence_sphere, divergence_sphere_wk
    use element_mod,    only : element_t
    use derivative_mod, only : derivative_t
    use kinds,          only : real_kind
    use coordinate_systems_mod, only: spherical_polar_t

    implicit none

    type (element_t)     , intent(in) :: elem(:)
    type (derivative_t)  , intent(in) :: deriv
    integer              , intent(in) :: nets,nete

    integer              , parameter :: intervals=6 

    real (kind=real_kind)              :: u(np,np)
    real (kind=real_kind)              :: laplace(np,np)
    real (kind=real_kind)              :: laplace_values(intervals,intervals)
    real (kind=real_kind)              :: laplace_fluxes(intervals,intervals,4)
    real (kind=real_kind)              :: laplace_test(intervals,intervals)
    real (kind=real_kind)              :: p,t

    type(spherical_polar_t)            :: s
    integer                            :: ie,i,j
    logical                            :: success
    success = .true.


    do ie=nets,nete

      call random_number(p)

      if (ie <= np*np) then
        u = 0 
        u(1+mod(ie,np), 1+mod(ie/np,np)) = 1
      else
        t = u(1+mod((7*ie),np), 1+mod((13*ie)/np,np))
        u(1+mod(ie,np), 1+mod(ie/np,np)) = t + 10*p
        u(1+mod(INT(32147*p),np), 1+mod(INT(1123*p)/np,np)) = 0
      end if

      laplace = laplace_sphere_wk(u,deriv,elem(ie),.false.)
      laplace = laplace / elem(ie)%spheremp
      laplace_values = subcell_integration(laplace, np, intervals, elem(ie)%metdet) 

      laplace_fluxes = subcell_Laplace_fluxes(u, deriv, elem(ie), np, intervals) 

      laplace_test = SUM(laplace_fluxes,3)

      laplace_test   = rearth*rearth*laplace_test
      laplace_values = rearth*rearth*laplace_values
      do i=1,intervals
      do j=1,intervals
        t = ABS(laplace_test(i,j)-laplace_values(i,j))/ &
            MAX(ABS(laplace_test(i,j)),ABS(laplace_values(i,j)))
        if (.0000001<t) then
          print *,__FILE__,__LINE__,ie,i,j,laplace_test(i,j),laplace_values(i,j),t
          success = .false.
        end if
      end do
      end do
    end do

    if (success) then
      print *,__FILE__,__LINE__," test_subcell_Laplace_fluxes test passed."
    else
      print *,__FILE__,__LINE__," test_subcell_Laplace_fluxes test FAILED."
    end if

  end subroutine test_subcell_Laplace_fluxes

  subroutine test_sub_integration(elem,deriv,nets,nete)
    use dimensions_mod, only : np
    use derivative_mod, only : subcell_integration
    use element_mod,    only : element_t
    use derivative_mod, only : derivative_t
    use kinds,          only : real_kind

    implicit none

    type (element_t)     , intent(in) :: elem(:)
    type (derivative_t)  , intent(in) :: deriv
    integer              , intent(in) :: nets,nete

    integer              , parameter :: intervals=6 

    real (kind=real_kind)              :: values(intervals,intervals)
    real (kind=real_kind)              :: V(np,np)
    real (kind=real_kind)              :: t, p
    integer                            :: ie,i,j
    logical                            :: success

    V = 0
    success = .true.
    do ie=nets,nete
      call random_number(p)
      if (ie <= np*np) then
        V = 0 
        V(1+mod(ie,np), 1+mod(ie/np,np)) = 1
      else
        t = V(1+mod((7*ie),np), 1+mod((13*ie)/np,np))
        V(1+mod(ie,np), 1+mod(ie/np,np)) = t + 10*p
        V(1+mod(INT(32147*p),np), 1+mod(INT(1123*p)/np,np)) = 0
      end if

      values = subcell_integration(V, np, intervals, elem(ie)%metdet) 

      t = 0
      do i = 1,np
        t    = t + DOT_PRODUCT(V(:,i),elem(ie)%spheremp(:,i))
      end do

      if (.00001<ABS(t-SUM(values))) then
        print *,__FILE__,__LINE__,ie,t,SUM(values)
        success = .false.
      end if
    end do
    if (success) then
      print *,__FILE__,__LINE__," test_sub_integration test passed."
    else
      print *,__FILE__,__LINE__," test_sub_integration test FAILED."
    end if

  end subroutine test_sub_integration


  subroutine test_edge_flux(elem,hybrid,deriv,nets,nete)
!
!  check local element vector dentities:
!*****
!  1. div and weak div are adjoints: (for all scalar test functions)
!     integral[  p div(u) ] + integral[ grad(p) dot u ] = boundary_integral[ p u dot n]
!       PHI div(u) spheremp - div_wk(u)(i,j) = boundary_integral[ u PHI]
!       where PHI = the delta function at (i,j)
!
!*****
!  2. grad and weak grad are adjoints: 
!     weak gradient is defined with CONTRA vector test functions
!     i.e. it returns vector:   [  integral[ p div(PHIcontra_1) ]       
!                               [  integral[ p div(PHIcontra_2) ]       
!     
!   integral[  p div(u) ] + integral[ grad(p) dot u ] = boundary_integral[ p u dot n]
! take u = PHIcontra_1 = (1,0) vector delta funciton at (i,j):
!  -grad_wk(p)_1(i,j) + spheremp PHIcontra_1 dot grad(p) = boundary_integral[ PHIcontra_1 p]
! and then take u = PHIcontra_2 = (0,1) vector delta function at (i,j):
!  -grad_wk(p)_2(i,j) + spheremp PHIcontra_2 dot grad(p) = boundary_integral[ PHIcontra_2 p]
!
! which is an equation for each covariant component:
! -grad_wk(p)_cov1 + spheremp grad(p)_cov1 = boundary_integral[ PHIcontra_1 p dot n]
! -grad_wk(p)_cov2 + spheremp grad(p)_cov2 = boundary_integral[ PHIcontra_2 p dot n]
!
! HOMME-SE works in latlon, so convert cov->lat/lon:
!
! -grad_wk(p) + spheremp grad(p) = D^-t * B 
!
! with
!    B1 = boundary_integral[ PHIcontra_1 p] 
!    B2 = boundary_integral[ PHIcontra_2 p]
!
!*****
! 3.  weak grid with COVARIANT test functions! 
!   integral[  p div(u) ] + integral[ grad(p) dot u ] = boundary_integral[ p u dot n]
! take u = PHIcov_1 = (1,0) vector delta funciton at (i,j):
!  -grad_wk(p)_1(i,j) + spheremp PHIcov_1 dot grad(p) = boundary_integral[ PHIcov_1 p]
! and then take u = PHIcov_2 = (0,1) vector delta function at (i,j):
!  -grad_wk(p)_2(i,j) + spheremp PHIcov_2 dot grad(p) = boundary_integral[ PHIcov_2 p]
!
! which is an equation for each CONTRA component:
! -grad_wk(p)_contra1 + spheremp grad(p)_contra1 = B1
! -grad_wk(p)_contra2 + spheremp grad(p)_contra2 = B2
!
! HOMME-SE works in latlon, so convert contra ->lat/lon:
!
! -grad_wk(p) + spheremp grad(p) = D * B 
!
! with
!    B1 = boundary_integral[ PHIcov_1 p] 
!    B2 = boundary_integral[ PHIcov_2 p]
!
!*****
! 4.  weak curl with COVARIANT test functions! 
!  integral[ u dot curl(v)] - integral[v dot curl(u)] = boundary_integral[ v cross u dot n]
!  curl(p) = curl(p*khat) = horizontal vector
!  vor(U) =  s*khat       = (which we treat as a scalar)
!   integral[ p * vor(u)  ] - integral[ u dot curl(p) ] = boundary_integral[ u cross p*khat  dot n]
!
! take u = PHIcov_1 = (1,0) vector delta funciton at (i,j):
!   curl_wk(p)_1(i,j) - spheremp PHIcov_1 dot curl(p) = boundary_integral[ perp(PHIcov_1) p]
! and then take u = PHIcov_2 = (0,1) vector delta function at (i,j):
!   curl_wk(p)_2(i,j) - spheremp PHIcov_2 dot curl(p) = boundary_integral[ perp(PHIcov_2) p]
!
! which is an equation for each CONTRA component:
! curl_wk(p)_contra1 - spheremp curl(p)_contra1 = B1
! curl_wk(p)_contra2 - spheremp curl(p)_contra2 = B2
!
! HOMME-SE works in latlon, so convert contra ->lat/lon:
!
! curl_wk(p) + spheremp curl(p) = D * B 
!
! with
!    B1 = boundary_integral[ PHIcov_1 p] 
!    B2 = boundary_integral[ PHIcov_2 p]
!
  use dimensions_mod, only : np, np, nlev
  use element_mod, only    : element_t
  use derivative_mod, only  : derivative_t, divergence_sphere, divergence_sphere_wk, &
                             element_boundary_integral, gradient_sphere, &
                             gradient_sphere_wk_testcontra,gradient_sphere_wk_testcov, &
                             curl_sphere, curl_sphere_wk_testcov
  use physical_constants, only : rrearth
  use kinds,          only : real_kind
  use hybrid_mod, only : hybrid_t

  implicit none
  
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  type (hybrid_t)      , intent(in) :: hybrid
  integer :: nets,nete
  ! local 
  real (kind=real_kind), dimension(np,np,2) :: ucontra,ulatlon,gradp,gradp_wk,ucov
  real (kind=real_kind), dimension(np,np) :: phidivu,ugradphi,rhs,lhs,p
  real (kind=real_kind), dimension(np,np) :: rhs2,lhs2
  integer :: i,j,ie,count

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hybrid%masterthread) print *,'integration by parts identity: check div/weak div:'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! test integration by parts identity for each Cardinal function PHI:
  ! div(u)*spheremp - div_wk(u) = boundary integral phi u dot n
  count = 0
  do ie=nets,nete
     call random_number(ucontra)
     ! contra->latlon
     ulatlon(:,:,1)=(elem(ie)%D(:,:,1,1)*ucontra(:,:,1) + elem(ie)%D(:,:,1,2)*ucontra(:,:,2))
     ulatlon(:,:,2)=(elem(ie)%D(:,:,2,1)*ucontra(:,:,1) + elem(ie)%D(:,:,2,2)*ucontra(:,:,2))
     phidivu = elem(ie)%spheremp(:,:)*divergence_sphere(ulatlon,deriv,elem(ie))
     ugradphi = divergence_sphere_wk(ulatlon,deriv,elem(ie))
     lhs = phidivu - ugradphi
     
     rhs = element_boundary_integral(ulatlon,deriv,elem(ie))
     
     
     do j=1,np
     do i=1,np
        if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
           write(*,'(a)') 'ERROR: div/div_wk integration by parts failure!'
           write(*,'(a,2i3,a,3e12.5)') 'for test function (i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),lhs(i,j)-rhs(i,j)
           count=count+1
        endif
     enddo
     enddo
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hybrid%masterthread) print *,'check grad/weak grad (gradient_sphere_wk_testcontra)'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PHIVEC = contra cardinal function 
  !          check each contra component seperately

  do ie=nets,nete
     call random_number(p)
     
     ! grad(p)  (lat/lon vector)
     gradp = gradient_sphere(p,deriv,elem(ie)%Dinv)
     gradp(:,:,1)=gradp(:,:,1)*elem(ie)%spheremp(:,:)  
     gradp(:,:,2)=gradp(:,:,2)*elem(ie)%spheremp(:,:)
     gradp_wk = gradient_sphere_wk_testcontra(p,deriv,elem(ie))
     
     ucontra(:,:,1)=p  ! PHIvec_1 * p
     ucontra(:,:,2)=0
     ! contra->latlon
     ulatlon(:,:,1)=(elem(ie)%D(:,:,1,1)*ucontra(:,:,1) + elem(ie)%D(:,:,1,2)*ucontra(:,:,2))
     ulatlon(:,:,2)=(elem(ie)%D(:,:,2,1)*ucontra(:,:,1) + elem(ie)%D(:,:,2,2)*ucontra(:,:,2))

     rhs = element_boundary_integral(ulatlon,deriv,elem(ie))
     lhs = gradp(:,:,1)-gradp_wk(:,:,1)

     ucontra(:,:,1)=0  ! PHIvec_2 * p
     ucontra(:,:,2)=p
     ! contra->latlon
     ulatlon(:,:,1)=(elem(ie)%D(:,:,1,1)*ucontra(:,:,1) + elem(ie)%D(:,:,1,2)*ucontra(:,:,2))
     ulatlon(:,:,2)=(elem(ie)%D(:,:,2,1)*ucontra(:,:,1) + elem(ie)%D(:,:,2,2)*ucontra(:,:,2))
     rhs2 = element_boundary_integral(ulatlon,deriv,elem(ie))
     lhs2 = gradp(:,:,2)-gradp_wk(:,:,2)  


     ! boundary integral gives covariant components. (see above) convert to latlon:
     ! cov -> latlon
     gradp(:,:,1)=rhs
     gradp(:,:,2)=rhs2
     rhs(:,:)=elem(ie)%Dinv(:,:,1,1)*gradp(:,:,1) + elem(ie)%Dinv(:,:,2,1)*gradp(:,:,2)
     rhs2(:,:)=elem(ie)%Dinv(:,:,1,2)*gradp(:,:,1) + elem(ie)%Dinv(:,:,2,2)*gradp(:,:,2)


     do j=1,np
     do i=1,np
        if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
           write(*,'(a)') 'ERROR: grad/grad_wk CONTRA (1) integration by parts failure!'
           write(*,'(a,2i3,a,4e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),&
                lhs(i,j)-rhs(i,j),lhs(i,j)/rhs(i,j)
           count=count+1
        endif
     enddo
     enddo


     do j=1,np
     do i=1,np
        if ( abs(lhs2(i,j)-rhs2(i,j)) .gt. 1d-20) then
           write(*,'(a)') 'ERROR: grad/grad_wk CONTRA (2) integration by parts failure!'
           write(*,'(a,2i2,a,3e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs2(i,j),rhs2(i,j),lhs2(i,j)-rhs2(i,j)
           count=count+1
        endif
     enddo
     enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hybrid%masterthread) print *,'check grad/weak grad (gradient_sphere_wk_testcov)'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete
     call random_number(p)

     
     ! grad(p)  (lat/lon vector)
     gradp = gradient_sphere(p,deriv,elem(ie)%Dinv)
     gradp(:,:,1)=gradp(:,:,1)*elem(ie)%spheremp(:,:)  
     gradp(:,:,2)=gradp(:,:,2)*elem(ie)%spheremp(:,:)
     gradp_wk = gradient_sphere_wk_testcov(p,deriv,elem(ie))
     lhs = gradp(:,:,1)-gradp_wk(:,:,1)
     lhs2 = gradp(:,:,2)-gradp_wk(:,:,2)  
     
     ucov(:,:,1)=p  ! PHIvec_1 * p
     ucov(:,:,2)=0
     ! cov->latlon
     ulatlon(:,:,1)=(elem(ie)%Dinv(:,:,1,1)*ucov(:,:,1) + elem(ie)%Dinv(:,:,2,1)*ucov(:,:,2))
     ulatlon(:,:,2)=(elem(ie)%Dinv(:,:,1,2)*ucov(:,:,1) + elem(ie)%Dinv(:,:,2,2)*ucov(:,:,2))
     rhs = element_boundary_integral(ulatlon,deriv,elem(ie))

     ucov(:,:,1)=0  ! PHIvec_2 * p
     ucov(:,:,2)=p
     ! cov->latlon
     ulatlon(:,:,1)=(elem(ie)%Dinv(:,:,1,1)*ucov(:,:,1) + elem(ie)%Dinv(:,:,2,1)*ucov(:,:,2))
     ulatlon(:,:,2)=(elem(ie)%Dinv(:,:,1,2)*ucov(:,:,1) + elem(ie)%Dinv(:,:,2,2)*ucov(:,:,2))
     rhs2 = element_boundary_integral(ulatlon,deriv,elem(ie))


     ! boundary integral gives contra components. (see above) convert to latlon:
     ! contra -> latlon
     gradp(:,:,1)=rhs
     gradp(:,:,2)=rhs2
     rhs(:,:) =elem(ie)%D(:,:,1,1)*gradp(:,:,1) + elem(ie)%D(:,:,1,2)*gradp(:,:,2)
     rhs2(:,:)=elem(ie)%D(:,:,2,1)*gradp(:,:,1) + elem(ie)%D(:,:,2,2)*gradp(:,:,2)


     do j=1,np
     do i=1,np
        if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
           write(*,'(a)') 'ERROR: grad/grad_wk COV (1) integration by parts failure!'
           write(*,'(a,2i2,a,4e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),lhs(i,j)-rhs(i,j),lhs(i,j)/rhs(i,j)
           count=count+1
        endif
     enddo
     enddo

     do j=1,np
     do i=1,np
        if ( abs(lhs2(i,j)-rhs2(i,j)) .gt. 1d-20) then
           write(*,'(a)') 'ERROR: grad/grad_wk COV (2) integration by parts failure!'
           write(*,'(a,2i2,a,3e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs2(i,j),rhs2(i,j),lhs2(i,j)-rhs2(i,j)
           count=count+1
        endif
     enddo
     enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hybrid%masterthread) print *,'check curl/weak curl (curl_sphere_wk_testcov)'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete
     call random_number(p)
     
     ! grad(p)  (lat/lon vector)
     gradp = curl_sphere(p,deriv,elem(ie))
     gradp(:,:,1)=gradp(:,:,1)*elem(ie)%spheremp(:,:)  
     gradp(:,:,2)=gradp(:,:,2)*elem(ie)%spheremp(:,:)
     gradp_wk = curl_sphere_wk_testcov(p,deriv,elem(ie))
     lhs =  gradp_wk(:,:,1)-gradp(:,:,1)
     lhs2 = gradp_wk(:,:,2)-gradp(:,:,2)
     
     ucov(:,:,1)=p  ! PHIvec_1 * p
     ucov(:,:,2)=0
     ! cov->latlon, and then u cross khat:
     ulatlon(:,:,2)=-(elem(ie)%Dinv(:,:,1,1)*ucov(:,:,1) + elem(ie)%Dinv(:,:,2,1)*ucov(:,:,2))
     ulatlon(:,:,1)= (elem(ie)%Dinv(:,:,1,2)*ucov(:,:,1) + elem(ie)%Dinv(:,:,2,2)*ucov(:,:,2))
     rhs = element_boundary_integral(ulatlon,deriv,elem(ie))

     ucov(:,:,1)=0  ! PHIvec_2 * p
     ucov(:,:,2)=p
     ! cov->latlon, and u cross khat:
     ulatlon(:,:,2)=-(elem(ie)%Dinv(:,:,1,1)*ucov(:,:,1) + elem(ie)%Dinv(:,:,2,1)*ucov(:,:,2))
     ulatlon(:,:,1)= (elem(ie)%Dinv(:,:,1,2)*ucov(:,:,1) + elem(ie)%Dinv(:,:,2,2)*ucov(:,:,2))
     rhs2 = element_boundary_integral(ulatlon,deriv,elem(ie))


     ! boundary integral gives contra components. (see above) convert to latlon:
     ! contra -> latlon
     gradp(:,:,1)=rhs
     gradp(:,:,2)=rhs2
     rhs(:,:) =elem(ie)%D(:,:,1,1)*gradp(:,:,1) + elem(ie)%D(:,:,1,2)*gradp(:,:,2)
     rhs2(:,:)=elem(ie)%D(:,:,2,1)*gradp(:,:,1) + elem(ie)%D(:,:,2,2)*gradp(:,:,2)


     do j=1,np
     do i=1,np
        if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
           write(*,'(a)') 'ERROR: curl/curl_wk COV (1) integration by parts failure!'
           write(*,'(a,2i2,a,4e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),lhs(i,j)-rhs(i,j),lhs(i,j)/rhs(i,j)
           count=count+1
        endif
     enddo
     enddo

     do j=1,np
     do i=1,np
        if ( abs(lhs2(i,j)-rhs2(i,j)) .gt. 1d-20) then
           write(*,'(a)') 'ERROR: curl/curl_wk COV (2) integration by parts failure!'
           write(*,'(a,2i2,a,3e12.4)') '(i,j)=',i,j,' lhs,rhs=',lhs2(i,j),rhs2(i,j),lhs2(i,j)-rhs2(i,j)
           count=count+1
        endif
     enddo
     enddo
  enddo

  if (hybrid%masterthread) print *,'done integration by parts identity check'
  if (count>0) stop 'ERROR: at least one integration by parts failuure'
  end subroutine

end module
