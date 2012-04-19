#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module state_mod
  ! ------------------------------
  use kinds, only : real_kind
  ! ------------------------------
  use dimensions_mod, only: nlev, np, nelemd
  ! ------------------------------
  use hybrid_mod, only : hybrid_t
  ! ------------------------------
  use element_mod, only : element_t
  ! ------------------------------
  use reduction_mod, only : parallelmin, parallelmax
  ! ------------------------------
  use global_norms_mod, only: wrap_repro_sum
  ! ------------------------------
  use parallel_mod, only: global_shared_buf, global_shared_sum
  ! ------------------------------
implicit none
private

  public :: printstate_dg
  public :: printstate
  public :: printstate_init

 contains

  subroutine printstate_init()

  end subroutine printstate_init

  
!======================================================================================================!  
  subroutine printstate_dg(elem,pmean,g,n0,hybrid,nets,nete)
    type(element_t), intent(in)  :: elem(:)
    real (kind=real_kind)        :: pmean,g
    integer                      :: n0
    type (hybrid_t),intent(in)   :: hybrid
    integer, intent(in)          :: nets,nete
    real (kind=real_kind)  :: umin_local(nets:nete),umax_local(nets:nete),usum_local(nets:nete), & 
			      vmin_local(nets:nete),vmax_local(nets:nete),vsum_local(nets:nete), &
			      pmin_local(nets:nete),pmax_local(nets:nete),psum_local(nets:nete)
    integer :: ie
    real (kind=real_kind) :: umin, vmin, pmin
    real (kind=real_kind) :: umax, vmax, pmax
    real (kind=real_kind) :: usum, vsum, psum
#ifdef _PRIMDG
    real (kind=real_kind) :: tmin_local(nets:nete),tmax_local(nets:nete),tsum_local(nets:nete)
    real (kind=real_kind) :: tsum, tmin, tmax
    integer, parameter :: nglen=4
#else
    real (kind=real_kind), dimension(np,np)  :: v1, v2
    real (kind=real_kind), dimension(np,np,2,nlev) ::ulatlon
    integer :: lev
    integer, parameter :: nglen=3
#endif
!=======================================================================================================!
    do ie=nets,nete
#ifdef _SWDG
!=== For DG 2D Version (Nair) ==================================
     ! Convert covariant wind to latlon 
     !  do lev=1,nlev
     !      v1 = elem(ie)%state%couv(:,:,1,lev)
     !      v2 = elem(ie)%state%couv(:,:,2,lev)
     !      ulatlon(:,:,1,lev)=elem(ie)%Dinv(1,1,:,:)*v1 + elem(ie)%Dinv(2,1,:,:)*v2   ! co->latlon
     !      ulatlon(:,:,2,lev)=elem(ie)%Dinv(1,2,:,:)*v1 + elem(ie)%Dinv(2,2,:,:)*v2   ! co->latlon
     !  end do
     ! Convert cotravariant wind to latlon 
        do lev=1,nlev
            v1 = elem(ie)%state%v(:,:,1,lev,n0)
            v2 = elem(ie)%state%v(:,:,2,lev,n0)
            ulatlon(:,:,1,lev)=elem(ie)%D(1,1,:,:)*v1 + elem(ie)%D(1,2,:,:)*v2 
            ulatlon(:,:,2,lev)=elem(ie)%D(2,1,:,:)*v1 + elem(ie)%D(2,2,:,:)*v2  
        end do
     !     ulatlon(:,:,1,:)=elem(ie)%state%v(:,:,1,:,n0)   !if output is in spherical (u,v) 
     !     ulatlon(:,:,2,:)=elem(ie)%state%v(:,:,2,:,n0)
     !     ulatlon(:,:,1,:)=elem(ie)%state%couv(:,:,1,:)   !if output is in spherical (u,v) 
     !     ulatlon(:,:,2,:)=elem(ie)%state%couv(:,:,2,:)
       umax_local(ie) = MAXVAL(ulatlon(:,:,1,:)) 
       vmax_local(ie) = MAXVAL(ulatlon(:,:,2,:))
       pmax_local(ie) = MAXVAL(elem(ie)%state%ht(:,:,:)) 
!======================================================
       umin_local(ie) = MINVAL(ulatlon(:,:,1,:)) 
       vmin_local(ie) = MINVAL(ulatlon(:,:,2,:))
       pmin_local(ie) = MINVAL(elem(ie)%state%ht(:,:,:))
!======================================================
       usum_local(ie) = SUM(ulatlon(:,:,1,:)) 
       vsum_local(ie) = SUM(ulatlon(:,:,2,:))
       psum_local(ie) = SUM(elem(ie)%state%ht(:,:,:)) 
!======================================================       
#else
#ifdef _PRIMDG
!=== For DG 3D version (nair) =============================
       umax_local(ie) = MAXVAL(elem(ie)%state%uv(:,:,1,:))
       vmax_local(ie) = MAXVAL(elem(ie)%state%uv(:,:,2,:))
       pmax_local(ie) = MAXVAL(elem(ie)%state%pr3d(:,:,nlev+1))       
       tmax_local(ie) = MAXVAL(elem(ie)%state%pt3d(:,:,:))
!======================================================
       umin_local(ie) = MINVAL(elem(ie)%state%uv(:,:,1,:))
       vmin_local(ie) = MINVAL(elem(ie)%state%uv(:,:,2,:))
       pmin_local(ie) = MINVAL(elem(ie)%state%pr3d(:,:,nlev+1))
       tmin_local(ie) = MINVAL(elem(ie)%state%pt3d(:,:,:))
!======================================================
       usum_local(ie) = SUM(elem(ie)%state%uv(:,:,1,:))
       vsum_local(ie) = SUM(elem(ie)%state%uv(:,:,2,:))
       psum_local(ie) = SUM(elem(ie)%state%pr3d(:,:,nlev+1))
       tsum_local(ie) = SUM(elem(ie)%state%pt3d(:,:,:))
!======================================================
!======================================================
#else
      umax_local(ie) = MAXVAL(elem(ie)%state%v(:,:,1,:,2))
      vmax_local(ie) = MAXVAL(elem(ie)%state%v(:,:,2,:,2))
      umin_local(ie) = MINVAL(elem(ie)%state%v(:,:,1,:,2))
      vmin_local(ie) = MINVAL(elem(ie)%state%v(:,:,2,:,2))
      usum_local(ie) = SUM(elem(ie)%state%v(:,:,1,:,2))
      vsum_local(ie) = SUM(elem(ie)%state%v(:,:,2,:,2))

#endif   
#endif
      global_shared_buf(ie,1) = usum_local(ie)
      global_shared_buf(ie,2) = vsum_local(ie)
      global_shared_buf(ie,3) = psum_local(ie)
# ifdef _PRIMDG
      global_shared_buf(ie,4) = tsum_local(ie)
# endif
!======================================================   
    end do

    umin = ParallelMin(umin_local,hybrid)
    umax = ParallelMax(umax_local,hybrid)

    vmin = ParallelMin(vmin_local,hybrid)
    vmax = ParallelMax(vmax_local,hybrid)

    pmin = ParallelMin(pmin_local,hybrid)
    pmax = ParallelMax(pmax_local,hybrid)

    call wrap_repro_sum(nvars=nglen, comm=hybrid%par%comm)
    usum = global_shared_sum(1)
    vsum = global_shared_sum(2)
    psum = global_shared_sum(3)

#ifdef _PRIMDG
    tmin = ParallelMin(tmin_local,hybrid)
    tmax = ParallelMax(tmax_local,hybrid)
    tsum = global_shared_sum(4)
#endif
    
    if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
      write (*,100) "u= ",umin,umax,usum
      write (*,100) "v= ",vmin,vmax,vsum
!======================================================
#ifdef _SWDG
      write (*,100) "ht= ",pmin,pmax,psum
#endif
#ifdef _PRIMDG
      write (*,100) "ps= ",pmin,pmax,psum      
      write (*,100) "temp= ",tmin,tmax,tsum
#endif   
!====================================================== 
      print *
    endif
    100 format (A5,3(E24.15))
     
  end subroutine printstate_dg
!======================================================================================================! 
!======================================================================================================!     
  subroutine printstate(elem,pmean,g,n0,hybrid,nets,nete, kmass)
    type(element_t), intent(in)  :: elem(:)
    real (kind=real_kind)        :: pmean,g
    integer                      :: n0
    type (hybrid_t),intent(in)   :: hybrid
    integer, intent(in)          :: nets,nete
  
    real (kind=real_kind)  :: umin_local(nets:nete),umax_local(nets:nete),usum_local(nets:nete), & 
			      vmin_local(nets:nete),vmax_local(nets:nete),vsum_local(nets:nete), &
			      pmin_local(nets:nete),pmax_local(nets:nete),psum_local(nets:nete)
    integer :: ie,k, i, j, k1,k2, k6

    real (kind=real_kind) :: umin, vmin, pmin
    real (kind=real_kind) :: umax, vmax, pmax
    real (kind=real_kind) :: usum, vsum, psum
    real (kind=real_kind), dimension(np, np) :: v1, v2
    real (kind=real_kind), dimension(np, np,2) :: ulatlon
!OG
    integer, intent(in)   :: kmass

    real (kind=real_kind)  :: tracer1(np*np), tracer2(np*np), w(np*np)

    real (kind=real_kind)  :: rm(nets:nete),rpu(nets:nete),os(nets:nete),weights_s(nets:nete)

    real (kind=real_kind)  :: rm_sum,rpu_sum,os_sum,weights_sum

!    real (kind=real_kind), parameter  :: surf=12.566370614359172*

    do k=1,nlev
    do ie=nets,nete

        ! Convert wind to lat-lon
       v1     = elem(ie)%state%v(:,:,1,k,n0)   ! contra
       v2     = elem(ie)%state%v(:,:,2,k,n0)   ! contra 
       ulatlon(:,:,1)=elem(ie)%D(1,1,:,:)*v1 + elem(ie)%D(1,2,:,:)*v2   ! contra->latlon
       ulatlon(:,:,2)=elem(ie)%D(2,1,:,:)*v1 + elem(ie)%D(2,2,:,:)*v2   ! contra->latlon

     if ((kmass>0).AND.(k.ne.kmass)) then
!======================================================
       umax_local(ie) = MAXVAL(ulatlon(:,:,1))
       vmax_local(ie) = MAXVAL(ulatlon(:,:,2))
#ifndef _PRIMDG
       pmax_local(ie) = MAXVAL((elem(ie)%state%p(:,:,k,n0)/elem(ie)%state%p(:,:,kmass,n0)+pmean)/g) &
		      + MAXVAL(elem(ie)%state%ps(:,:))
#endif
!======================================================
       umin_local(ie) = MINVAL(ulatlon(:,:,1))
       vmin_local(ie) = MINVAL(ulatlon(:,:,2))
#ifndef _PRIMDG
       pmin_local(ie) = MINVAL((elem(ie)%state%p(:,:,k,n0)/elem(ie)%state%p(:,:,kmass,n0)+pmean)/g) &
		      + MINVAL(elem(ie)%state%ps(:,:))
#endif
!======================================================
       usum_local(ie) = SUM(ulatlon(:,:,1))
       vsum_local(ie) = SUM(ulatlon(:,:,2))
       ! lets output tracer mass, instead of meaningless sum of grid point values:
       psum_local(ie) = SUM((elem(ie)%state%p(:,:,k,n0)*elem(ie)%spheremp(:,:)))
!======================================================
    else

!======================================================
       umax_local(ie) = MAXVAL(ulatlon(:,:,1))
       vmax_local(ie) = MAXVAL(ulatlon(:,:,2))
#ifndef _PRIMDG
       pmax_local(ie) = MAXVAL((elem(ie)%state%p(:,:,k,n0)+pmean)/g) &
                      + MAXVAL(elem(ie)%state%ps(:,:))
#endif
!======================================================
       umin_local(ie) = MINVAL(ulatlon(:,:,1))
       vmin_local(ie) = MINVAL(ulatlon(:,:,2))
#ifndef _PRIMDG
       pmin_local(ie) = MINVAL((elem(ie)%state%p(:,:,k,n0)+pmean)/g) &
                      + MINVAL(elem(ie)%state%ps(:,:))
#endif
!======================================================
       usum_local(ie) = SUM(ulatlon(:,:,1))
       vsum_local(ie) = SUM(ulatlon(:,:,2))
#ifdef _PRIMDG
       pmax_local(ie) = MAXVAL((elem(ie)%state%p(:,:,k,n0)+pmean)/g) &
                      + MAXVAL(elem(ie)%state%phis(:,:))
       pmin_local(ie) = MINVAL((elem(ie)%state%p(:,:,k,n0)+pmean)/g) &
                      + MINVAL(elem(ie)%state%phis(:,:))
       psum_local(ie) = SUM((elem(ie)%state%p(:,:,k,n0)+pmean)/g) &
                      + SUM(elem(ie)%state%phis(:,:))
#else
       psum_local(ie) = SUM((elem(ie)%state%p(:,:,k,n0)+pmean)/g) &
                      + SUM(elem(ie)%state%ps(:,:))
#endif
!======================================================
     endif

     global_shared_buf(ie,1) = usum_local(ie)
     global_shared_buf(ie,2) = vsum_local(ie)
     global_shared_buf(ie,3) = psum_local(ie)


    end do

    umin = ParallelMin(umin_local,hybrid)
    umax = ParallelMax(umax_local,hybrid)

    vmin = ParallelMin(vmin_local,hybrid)
    vmax = ParallelMax(vmax_local,hybrid)

    pmin = ParallelMin(pmin_local,hybrid)
    pmax = ParallelMax(pmax_local,hybrid)

    call wrap_repro_sum(nvars=3, comm=hybrid%par%comm)
    usum = global_shared_sum(1)
    vsum = global_shared_sum(2)
    psum = global_shared_sum(3)

    if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
      if (k==1) then 
         write (*,100) "u= ",umin,umax,usum
         write (*,100) "v= ",vmin,vmax,vsum
      endif
      write (*,100) "p= ",pmin,pmax,psum
      if (k==nlev) print *,'sqrt(g(h0+max(h)) = ',sqrt(pmean+pmax)
      if (k==nlev) print *
    endif
    100 format (A5,3(E24.15))
    110 format (A57,3(E24.15))
    120 format (A20,3(E24.15))
    enddo

!!!!this code is to see if limiters preserve const sum of two tracers
!made for swirl only and right now is not set to work(?)
    if((nlev>=5).and.(kmass>0))then
      k1=5 ! compiler will complain if code access index 5 below and nlev<5:
      do ie=nets,nete
	pmax_local(ie) = MAXVAL(elem(ie)%state%p(:,:,1,n0)/elem(ie)%state%p(:,:,kmass,n0) &
			+ elem(ie)%state%p(:,:,k1,n0)/elem(ie)%state%p(:,:,kmass,n0))
	pmin_local(ie) = MINVAL(elem(ie)%state%p(:,:,1,n0)/elem(ie)%state%p(:,:,kmass,n0) &
			+ elem(ie)%state%p(:,:,k1,n0)/elem(ie)%state%p(:,:,kmass,n0))
      end do

      pmin = ParallelMin(pmin_local,hybrid)
      pmax = ParallelMax(pmax_local,hybrid)

      if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
	write (*,110) "min/max(TRACER1+TRACER5)(meaningless if not swirl)= ",pmin,pmax
      endif
    endif
!!!!this code is to see if limiters preserve sum of two tracers = third tracer
!made for swirl only and right now is not set to work(?),
!requires initializetion level 7 = level 2 + level 6
    if((nlev>=7).and.(kmass>0))then
      do ie=nets,nete
        k6=6
        k2=2
	pmax_local(ie) = MAXVAL(elem(ie)%state%p(:,:,k2,n0)/elem(ie)%state%p(:,:,kmass,n0) &
			+ elem(ie)%state%p(:,:,k6,n0)/elem(ie)%state%p(:,:,kmass,n0) &
			- elem(ie)%state%p(:,:,k6+1,n0)/elem(ie)%state%p(:,:,kmass,n0) )
	pmin_local(ie) = MINVAL(elem(ie)%state%p(:,:,k2,n0)/elem(ie)%state%p(:,:,kmass,n0) &
			+ elem(ie)%state%p(:,:,k6,n0)/elem(ie)%state%p(:,:,kmass,n0) &
			- elem(ie)%state%p(:,:,k6+1,n0)/elem(ie)%state%p(:,:,kmass,n0) )
      end do

      pmin = ParallelMin(pmin_local,hybrid)
      pmax = ParallelMax(pmax_local,hybrid)

      if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
	write (*,110) "min/max(TR2+TR6-TR7)(meaningless if not swirl)= ",pmin,pmax
      endif
    endif

!this is mixing diagnostics, makes sense only in swirl because
!it requires setup of layers 1 and 6
    if(nlev>5)then
      do ie=nets,nete
	if(kmass>0)then
          k1=0
          do i=1,np
            do j=1,np
              k1=k1+1
              k6=6
              tracer1(k1)=elem(ie)%state%p(i,j,1,n0)/elem(ie)%state%p(i,j,kmass,n0)
              tracer2(k1)=elem(ie)%state%p(i,j,k6,n0)/elem(ie)%state%p(i,j,kmass,n0)
              w(k1)=elem(ie)%spheremp(i,j)
            enddo
          enddo
	else
          k1=0
          do i=1,np
            do j=1,np
              k1=k1+1
              k6=6
              tracer1(k1)=elem(ie)%state%p(i,j,1,n0)
              tracer2(k1)=elem(ie)%state%p(i,j,k6,n0)
              w(k1)=elem(ie)%spheremp(i,j)
            enddo
          enddo
        endif
	call correlation_diag(tracer1,tracer2,w,np*np,rm(ie),rpu(ie),os(ie))
        weights_s(ie)=sum(w)
        global_shared_buf(ie,1) = rm(ie)
        global_shared_buf(ie,2) = rpu(ie)
        global_shared_buf(ie,3) = os(ie)
        global_shared_buf(ie,4) = weights_s(ie)
      end do

      call wrap_repro_sum(nvars=4, comm=hybrid%par%comm)
      rm_sum = global_shared_sum(1)
      rpu_sum = global_shared_sum(2)
      os_sum = global_shared_sum(3)
!redundant, weights_sum is a const always
      weights_sum = global_shared_sum(4)

      if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
	if(weights_sum>0)then
	  write (*,120) "real mixing = ",rm_sum/weights_sum
	  write (*,120) "range_pres_unmixing = ",rpu_sum/weights_sum
	  write (*,120) "overshooting = ",os_sum/weights_sum
	  write (*,120) "weightsum = ",weights_sum
	else
	  write (*,120) 'Warning: mixing diagnostics: by some reason area of sphere is not positive?'
	endif
      endif
    endif

  end subroutine printstate
!======================================================================================================!
!======================================================================================================!

!!!!!! here is PL's code for mixing ratios
!for workshop in 2011
  SUBROUTINE correlation_diag(f1,f2,w,K,real_mixing,range_pres_unmixing,overshooting)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: K
    real (kind=real_kind), DIMENSION(K)      , INTENT(IN) :: f1,f2,w
    !
    ! local workspace
    !
    real(KIND=real_kind) :: root, tol,q1,q2,c
    real (kind=real_kind), INTENT(OUT)     :: real_mixing,overshooting,range_pres_unmixing
    INTEGER  :: j
    real(KIND=real_kind) :: q1_min,q1_max,q2_min,q2_max, a, cc, d, xx, tiny
    
    q1_min = 0.1d0
    q1_max = 1.0d0
    q2_min = 0.1d0
    q2_max = 0.892d0
    
    real_mixing          = 0.0d0
    overshooting         = 0.0d0
    range_pres_unmixing  = 0.0d0

    tiny = 1.0E-10

    DO j=1,K

       q1 = f1(j)
       q2 = f2(j)

!THESE FORMULAS FOR ROOT IN UNSCALED DISTANCE DO NOT DELETE
!        c=DBLE(5400)*q1+DBLE(6)*SQRT(-DBLE(7986)+DBLE(87120)*q2-DBLE(316800)*(q2**2)+&
!             DBLE(384000)*(q2**3)+DBLE(810000)*(q1**2))
!        c=c**(DBLE(1)/DBLE(3))
!        c=c/DBLE(24)
!        root=c-((-DBLE(11)/DBLE(96)+DBLE(5)*q2/DBLE(12))/c)


!        c=5400.0q0*q1+6.0q0*SQRT(-7986.0q0+87120.0q0*q2-316800.0q0*(q2**2)+&
!             384000.0q0*(q2**3)+810000.0q0*(q1**2))
!        c=c**(1.0q0/3.0q0)
!        c=c/24.0q0
!        root=c-((-11.0q0/96.0q0+5.0q0*q2/12.0q0)/c)
!-------------------------------- DO NOT DELETE

     ! IF (dist_fct(q1,q1,q2)>tiny) THEN
        
       c=(DBLE(65340)*q1+DBLE(12)*SQRT(-DBLE(1687296)+DBLE(12168000)*q2&
            -DBLE(29250000)*q2**2+DBLE(23437500)*q2**3+DBLE(29648025)*q1**2))**(DBLE(1)/DBLE(3))
       c=c/(DBLE(60))
       root = c-(-(DBLE(13)/DBLE(75))+(DBLE(5)/DBLE(12))*q2)/c
  

      
!         c=(65340.0q0 *q1+12.0q0*SQRT(-1687296.0q0+12168000.0q0*q2&
!              -29250000.0q0*q2**2+23437500.0q0*q2**3+29648025.0q0*q1**2))**(1.0q0/3.0q0)
!         c=c/(60.0q0)
!         root = c-(-(13.0q0/75.0q0)+(5.0q0/12.0q0)*q2)/c

	root = MAX(0.1D0,root)
	root = MIN(1.0D0,root)

!   write(6,*) 'aaa0',q1, q2
!   write(6,*) 'aaa1',q2, corr_fct(q1)
!   write(6,*) 'aaa2',dist_fct(root,q1,q2), root, dist_fct(q1,q1,q2)
! 
!   write(6,*) 'aaa3',a,cc,d
! 
!   stop
	IF (q2<corr_fct(q1).AND.q2>line_fct(q2,q1_min,q1_max,q2_min,q2_max)) THEN
	    !
	    ! `real' mixing
	    ! 
	    real_mixing = real_mixing + dist_fct(root,q1,q2)*w(j)                 
	ELSE IF (q1<q1_max.AND.q1>q1_min.AND.q2<q2_max.AND.q2>q2_min.AND.&
	      q2>line_fct(q2,q1_min,q1_max,q2_min,q2_max)) THEN
	    !
	    ! range-preserving unmixing
	    ! 
	    range_pres_unmixing = range_pres_unmixing+dist_fct(root,q1,q2)*w(j)                   
	ELSE
	    !
	    ! overshooting
	    !
	    overshooting = overshooting + dist_fct(root,q1,q2)*w(j)   
	END IF
   !    endif
    END DO
  END SUBROUTINE correlation_diag

  !
  ! correlation function
  !
  real (kind=real_kind) FUNCTION corr_fct(x)
  IMPLICIT NONE
    real (kind=real_kind), INTENT(IN)  :: x
    corr_fct = -0.8d0*x**2+0.9d0
  END FUNCTION corr_fct
  !
  ! Eucledian distance function
  !
  real (kind=real_kind) FUNCTION dist_fct(x,x0,y0)
    IMPLICIT NONE
    real (kind=real_kind), INTENT(IN)  :: x,x0,y0
    dist_fct = SQRT((x-x0)*(x-x0)/(0.9d0**2)+(corr_fct(x)-y0)*(corr_fct(x)-y0)/(0.792d0**2))
  END FUNCTION dist_fct
  !
  ! straight line line function
  !
  real (kind=real_kind) FUNCTION line_fct(x,xmin,xmax,ymin,ymax)
    IMPLICIT NONE
    real (kind=real_kind), INTENT(IN)  :: x,xmin,xmax,ymin,ymax
    real (kind=real_kind) :: a,b
    !
    ! line: y=a*x+b
    ! 
    a = (ymax-ymin)/(xmax-xmin)
    b = ymin-xmin*a
    line_fct = a*x+b
  END FUNCTION line_fct







end module state_mod
