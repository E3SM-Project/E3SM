!-----------------------------------------------------------------------------------!
!MODULE SPELT_MOD-------------------------------------------------------CE-for SPELT!
! SPELT_MOD File for the spelt project in HOMME                                      !
! Author: Christoph Erath                                                           !
! Date: 24.August 2012                                                              !
! MAIN module to run SPELT on HOMME                                                 !
!-----------------------------------------------------------------------------------!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module spelt_bench_mod
  use kinds, only : real_kind, longdouble_kind, int_kind
contains

subroutine spelt_run_bench(elem,spelt,hybrid,nets,nete,tl)
  ! ---------------------------------------------------------------------------------
  use fvm_bsp_mod, only: analytical_function, get_boomerang_velocities_gll,get_solidbody_velocities_gll
  ! ---------------------------------------------------------------------------------
  use derivative_mod, only : derivative_t, derivinit, v2pinit, &
                             interpolate_gll2spelt_points
                             
  use coordinate_systems_mod, only : cartesian2D_t                          
  ! ---------------------------------------------------------------------------------
  use edge_mod, only :  ghostVpack2d, ghostVunpack2d, freeghostbuffertr, freeedgebuffer
  use dimensions_mod, only: ne, np, nlev, ntrac, nc, nhe, nip, nipm, nep
  use element_mod, only : element_t, timelevels
  use hybrid_mod, only : hybrid_t
  use spelt_mod, only: cellghostbuf, factorR, edgeveloc, spelt_struct
  use spelt_mod, only: spelt_mcgregordss, spelt_rkdss,spelt_run, spelt_runlimit, spelt_runpos, spelt_runair, spelt_init3
  use spelt_mod, only: cip_coeff, cip_interpolate, metric_term,  cell_search, qmsl_cell_filter, cell_minmax, cip_cell_avr
  ! ---------------------------------------------------------------------------------
  use bndry_mod, only: ghost_exchangeV
  ! ---------------------------------------------------------------------------------
  use time_mod, only : tstep, nmax, timelevel_t, time_at, timelevel_update
  ! ---------------------------------------------------------------------------------
  use parallel_mod, only: global_shared_buf, global_shared_sum
  ! ---------------------------------------------------------------------------------
  use global_norms_mod, only: wrap_repro_sum
  ! ---------------------------------------------------------------------------------
  use reduction_mod, only : parallelmax, parallelmin
  ! ---------------------------------------------------------------------------------  
  use physical_constants, only : g, rearth, DD_PI
  ! ------EXTERNAL----------------
  use perf_mod, only : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  ! -----------------------------------------------
  use quadrature_mod, only : quadrature_t, gausslobatto
  use control_mod, only : test_cfldep
  
  
#ifdef PIO_INTERP
    use interp_movie_mod, only : interp_movie_init, interp_movie_output, interp_movie_finish
#else
    use shal_movie_mod, only : shal_movie_init, shal_movie_output, shal_movie_finish
#endif
  
  implicit none
  type (element_t), intent(inout)             :: elem(:)
  type (spelt_struct), intent(inout)           :: spelt(:)
  
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)

  type (TimeLevel_t)                          :: tl              ! time level struct
  type (derivative_t)                         :: deriv           ! derivative struct

  real (kind=real_kind)                       :: massstart, mass, maxc, maxcstart,minc, mincstart
  real (kind=real_kind)                       :: tmp1(nets:nete), tmp2(nets:nete)
  real (kind=real_kind)                       :: l1,l2, lmax, tmpref, area
  
  integer                                     :: i,j,k,ie,itr

  real (kind=real_kind)                       :: tmp
  
  real (kind=real_kind)                       :: ff(nip,nip)
  real (kind=real_kind)                       :: cf(nip,nip,1-nhe:nc+nhe,1-nhe:nc+nhe)
  real (kind=real_kind)                       :: minmax(1-nhe:nc+nhe,1-nhe:nc+nhe,2)
  
  real (kind=longdouble_kind)                 :: refnep(1:nep)
  
  integer                                     :: icell, jcell
  real (kind=real_kind)                       :: sga, dx, dy
  type (cartesian2D_t)                        :: dref, alphabeta
  real (kind=real_kind)                       :: vstar(np,np,2), v1, v2, u(1:nep,1:nep), v(1:nep,1:nep)
  real (kind=real_kind)                       :: maxcflx=0.0D0, maxcfly=0.0D0,maxcflxstep, maxcflystep  
  
  type (quadrature_t) :: gp   ! Quadrature points and weights on pressure grid
  
  integer                                     :: face_nodep
  integer  choosetrac, chooselev   !for test reason the output
 !-----------------------------------------------------------------------------------!  
 choosetrac=1
 chooselev=1
 
  if(hybrid%masterthread) then 
    print *,"!-----------------------------------------------------------------------!"
    print *,"!  Test CASE for SPELT, Christoph Erath                                 !" 
    print *,"!-----------------------------------------------------------------------!" 
  endif
     
  ! Initialize derivative structure
  ! SPELT nodes are equally spaced in alpha/beta
  ! HOMME with equ-angular gnomonic projection maps alpha/beta space
  ! to the reference element via simple scale + translation
  
  call derivinit(deriv)
  
  gp=gausslobatto(np)
! ! for the np grid
!   do i=1,nc
!     tmp=abs(gp%points(i+1)-gp%points(i))/nipm
!     refnep(1+(i-1)*nipm)=gp%points(i)
!     refnep(2+(i-1)*nipm)=gp%points(i)+tmp
!   end do
!   refnep(nep)=gp%points(np)
!  for the nc grid  
  tmp=nep-1
  do i=1,nep
    refnep(i)= 2*(i-1)/tmp - 1
  end do

  call v2pinit(deriv%Sfvm,gp%points,refnep,np,nep)
  
!Initialize test example with first communication-----------------------------------!      
  do ie=nets,nete
    do k=1, nlev
      do itr=1,ntrac
        do j=1,nep
          do i=1,nep  
        ! define test example
            if ((itr==0)) then
              spelt(ie)%c(i,j,k,itr,tl%n0)=1.0D0    !density of the air
            else
              call analytical_function(spelt(ie)%c(i,j,k,itr,tl%n0),spelt(ie)%asphere(i,j),k,itr)
            endif
            if ((mod(i,2)==0) .and.(mod(j,2)==0) .and. (itr==choosetrac)) then
               spelt(ie)%cstart(i/2,j/2)=spelt(ie)%c(i,j,chooselev,choosetrac,tl%n0)
            endif 
          end do
        end do
      end do 
    end do
  end do
  
  !first exchange of the initial values
  call spelt_init3(elem,spelt,hybrid,nets,nete,tl%n0)
!-----------------------------------------------------------------------------------!
! ONLY for output on ps needed
  do ie=nets,nete  
    do k=1, nlev
      do itr=1,ntrac
        do j=1-nhe,nc+nhe
          do i=1-nhe,nc+nhe
            icell=1+(i-1)*nipm
            jcell=1+(j-1)*nipm
            ff=spelt(ie)%c(icell:icell+nipm,jcell:jcell+nipm,k,itr,tl%n0)*spelt(ie)%sga(icell:icell+nipm,jcell:jcell+nipm)
            ! FOR TEST OUTPUT on gll only needed
            minmax(i,j,:)=cell_minmax(ff)
            call cip_coeff(spelt(ie)%drefx(i,j),spelt(ie)%drefy(i,j),ff,ff(2,2),cf(:,:,i,j))
          enddo
        enddo
        if ((itr==choosetrac) .and. (k==chooselev)) then
          do j=1,np    ! FOR TEST OUTPUT on gll only needed
            do i=1,np
              call cell_search(elem(ie),spelt(ie), elem(ie)%spherep(i,j), icell, jcell,dref, alphabeta,face_nodep)
              tmp=cip_interpolate(cf(:,:,icell,jcell),dref%x,dref%y)   
              tmp=qmsl_cell_filter(icell,jcell,minmax,tmp) 
              sga=metric_term(alphabeta) 
!               sga= metric_termref(elem(ie),dref)
              elem(ie)%state%p(i,j,1,tl%n0)=g*tmp/sga     
            enddo
          enddo
        endif 
      enddo
    enddo
    elem(ie)%state%ps(:,:)=0.0D0
  enddo 
!-----------------------------------------------------------------------------------! 
#ifdef PIO_INTERP
  call interp_movie_init(elem,hybrid,nets,nete,tl=tl)    
  call interp_movie_output(elem,tl, hybrid, 0D0, nets, nete,spelt)
#else
  call shal_movie_init(elem,hybrid)
  call shal_movie_output(elem,tl, hybrid, 0D0, nets, nete,deriv)
#endif
! 
!-----------------------------------------------------------------------------------! 
! for Quantities as mass aso 
  do ie=nets,nete
   spelt(ie)%elem_mass=0.0D0
   tmp1(ie)=-1.0D0-20
   tmp2(ie)=1.0D20
   do j=1,nc
     do i=1,nc
       area=spelt(ie)%area_sphere(i,j)
       icell=2+(i-1)*nipm
       jcell=2+(j-1)*nipm
       if (choosetrac==1) then   ! mass of air, code is not optimal
         spelt(ie)%elem_mass=spelt(ie)%elem_mass + &
                        area*spelt(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0)  ! that is the real mass, we have cell average
       elseif(choosetrac==2) then
         spelt(ie)%elem_mass=spelt(ie)%elem_mass + &
                       area*spelt(ie)%c(icell,jcell,chooselev,1,tl%n0)*spelt(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0)
       endif
!          spelt(ie)%cstart(i,j)=spelt(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0)
      tmp1(ie)=max(tmp1(ie),spelt(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0))
      tmp2(ie)=min(tmp2(ie),spelt(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0))
     enddo
   enddo
   ! for the mass value
   global_shared_buf(ie,1)=0.0D0
   global_shared_buf(ie,1)=spelt(ie)%elem_mass
   ! for the max value on the sphere
!    tmp1(ie) = MAXVAL(spelt(ie)%c(:,:,chooselev,choosetrac,tl%n0))
!    tmp2(ie) = MINVAL(spelt(ie)%c(:,:,chooselev,choosetrac,tl%n0))   
  end do 

  !need the buffer cellghostbuf in the time loop
  ! for mass calculation
  call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
  massstart=global_shared_sum(1)
  
  maxcstart = parallelmax(tmp1,hybrid)
  mincstart = parallelmin(tmp2,hybrid) 
!-----------------------------------------------------------------------------------!
  call t_barrierf('spelt time loop', hybrid%par%comm)
  call t_startf('spelt')
  
DO WHILE(tl%nstep<nmax) !nmax)
! ! start mcgregordss
  do ie=nets,nete
    do k=1,nlev
      !simulate gll velocities, I get the velocities on the gll from the Spectral Element simulation
      elem(ie)%derived%vstar(:,:,:,k)=get_boomerang_velocities_gll(elem(ie), time_at(tl%nstep+1))
      spelt(ie)%vn0(:,:,:,k)=get_boomerang_velocities_gll(elem(ie),time_at(tl%nstep))
      spelt(ie)%vn12(:,:,:,k)=get_boomerang_velocities_gll(elem(ie),(time_at(tl%nstep)+time_at(tl%nstep+1))/2.0D0)

!       elem(ie)%derived%vstar(:,:,:,k)=get_solidbody_velocities_gll(elem(ie), time_at(tl%nstep+1))
!       spelt(ie)%vn0(:,:,:,k)=get_solidbody_velocities_gll(elem(ie),time_at(tl%nstep))
!       spelt(ie)%vn1(:,:,:,k)=get_solidbody_velocities_gll(elem(ie),(time_at(tl%nstep)+time_at(tl%nstep+1))/2)


!       vstar=elem(ie)%derived%vstar(:,:,:,k)/rearth
! !       do j=1,np
! !         do i=1,np
! !           v1 = spelt(ie)%Dinv(1,1,i,j)*vstar(i,j,1) + spelt(ie)%Dinv(1,2,i,j)*vstar(i,j,2)
! !           v2 = spelt(ie)%Dinv(2,1,i,j)*vstar(i,j,1) + spelt(ie)%Dinv(2,2,i,j)*vstar(i,j,2)
! !           vstar(i,j,1)=v1
! !           vstar(i,j,2)=v2
! !         enddo
! !       enddo
! !     spelt(ie)%contrau(:,:,k)=interpolate_gll2spelt_points(vstar(:,:,1),deriv)
! !     spelt(ie)%contrav(:,:,k)=interpolate_gll2spelt_points(vstar(:,:,2),deriv)
! 
!     u=interpolate_gll2spelt_points(vstar(:,:,1),deriv)
!     v=interpolate_gll2spelt_points(vstar(:,:,2),deriv)
!     do j=1,nep
!       do i=1,nep
!         spelt%contrau(i,j,k)=spelt%Ainv(1,1,i,j)*u(i,j)+spelt%Ainv(2,1,i,j)*v(i,j)      
!         spelt%contrav(i,j,k)=spelt%Ainv(1,2,i,j)*u(i,j)+spelt%Ainv(2,2,i,j)*v(i,j)
!       end do
!     end do
!     
!     vstar=spelt(ie)%vn12(:,:,:,k)/rearth
! !     do j=1,np
! !       do i=1,np
! !         v1 = spelt(ie)%Dinv(1,1,i,j)*vstar(i,j,1) + spelt(ie)%Dinv(1,2,i,j)*vstar(i,j,2)
! !         v2 = spelt(ie)%Dinv(2,1,i,j)*vstar(i,j,1) + spelt(ie)%Dinv(2,2,i,j)*vstar(i,j,2)
! !         vstar(i,j,1)=v1
! !         vstar(i,j,2)=v2
! !       enddo
! !     enddo
! !     spelt(ie)%contrau1=interpolate_gll2spelt_points(vstar(:,:,1),deriv)
! !     spelt(ie)%contrav1=interpolate_gll2spelt_points(vstar(:,:,2),deriv)
!     u=interpolate_gll2spelt_points(vstar(:,:,1),deriv)
!     v=interpolate_gll2spelt_points(vstar(:,:,2),deriv)
!     do j=1,nep
!       do i=1,nep
!         spelt%contrau1(i,j)=spelt%Ainv(1,1,i,j)*u(i,j)+spelt%Ainv(2,1,i,j)*v(i,j)      
!         spelt%contrav1(i,j)=spelt%Ainv(1,2,i,j)*u(i,j)+spelt%Ainv(2,2,i,j)*v(i,j)
!       end do
!     end do    
!     
!       vstar=spelt(ie)%vn0(:,:,:,k)/rearth
! !     do j=1,np
! !       do i=1,np
! !         v1 = spelt(ie)%Dinv(1,1,i,j)*vstar(i,j,1) + spelt(ie)%Dinv(1,2,i,j)*vstar(i,j,2)
! !         v2 = spelt(ie)%Dinv(2,1,i,j)*vstar(i,j,1) + spelt(ie)%Dinv(2,2,i,j)*vstar(i,j,2)
! !         vstar(i,j,1)=v1
! !         vstar(i,j,2)=v2
! !       enddo
! !     enddo
! !   spelt(ie)%contrau2=interpolate_gll2spelt_points(vstar(:,:,1),deriv)
! !   spelt(ie)%contrav2=interpolate_gll2spelt_points(vstar(:,:,2),deriv)
!     u=interpolate_gll2spelt_points(vstar(:,:,1),deriv)
!     v=interpolate_gll2spelt_points(vstar(:,:,2),deriv)
!     do j=1,nep
!       do i=1,nep
!         spelt%contrau2(i,j)=spelt%Ainv(1,1,i,j)*u(i,j)+spelt%Ainv(2,1,i,j)*v(i,j)      
!         spelt%contrav2(i,j)=spelt%Ainv(1,2,i,j)*u(i,j)+spelt%Ainv(2,2,i,j)*v(i,j)
!       end do
!     end do    
!     
! !       u=interpolate_gll2spelt_points(vstar(:,:,1),deriv)
! !       v=interpolate_gll2spelt_points(vstar(:,:,2),deriv)
! !       do j=1,nep
! !         do i=1,nep
! !           spelt%contrau(i,j,k)=spelt%Ainv(1,1,i,j)*u(i,j)+spelt%Ainv(2,1,i,j)*v(i,j)      
! !           spelt%contrav(i,j,k)=spelt%Ainv(1,2,i,j)*u(i,j)+spelt%Ainv(2,2,i,j)*v(i,j)
! !         end do
! !       end do
    end do
  end do
  
!   call spelt_mcgregordss(elem,spelt,nets,nete, hybrid, deriv, tstep, 3)
!   call spelt_rkdss(elem,spelt,nets,nete, hybrid, deriv, tstep, 3)
! ! end mcgregordss
!   call spelt_runlimit(elem,spelt,hybrid,deriv,tstep,tl,nets,nete)
!   call spelt_runpos(elem,spelt,hybrid,deriv,tstep,tl,nets,nete)
  call spelt_run(elem,spelt,hybrid,deriv,tstep,tl,nets,nete)
!   call spelt_runair(elem,spelt,hybrid,deriv,tstep,tl,nets,nete)

  call TimeLevel_update(tl,"forward") 
!-----------------------------------------------------------------------------------! 
! for Quantities as mass aso
if (mod(tl%nstep,1)==0) then
  do ie=nets,nete
    spelt(ie)%elem_mass=0.0D0
    tmp1(ie)=-1.0D0-20
    tmp2(ie)=1.0D20
    do j=1,nc
      do i=1,nc
        area=spelt(ie)%area_sphere(i,j)
        icell= 2+(i-1)*nipm
        jcell=2+(j-1)*nipm
        sga=spelt(ie)%sga(icell,jcell)
        if (choosetrac==1) then   ! mass of air, code is not optimal
         spelt(ie)%elem_mass=spelt(ie)%elem_mass + &
                        area*spelt(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0)
        elseif (choosetrac==2) then
          spelt(ie)%elem_mass=spelt(ie)%elem_mass + &
                         area*spelt(ie)%c(icell,jcell,chooselev,1,tl%n0)*spelt(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0)
        endif
        
        tmp1(ie)=max(tmp1(ie),spelt(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0))
        tmp2(ie)=min(tmp2(ie),spelt(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0))
      enddo
    enddo
    ! for the mass value
    global_shared_buf(ie,1)=0.0D0
    global_shared_buf(ie,1)=spelt(ie)%elem_mass
    ! for the max value on the sphere
    
!     tmp1(ie) = MAXVAL(spelt(ie)%c(:,:,chooselev,choosetrac,tl%n0))
!     tmp2(ie) = MINVAL(spelt(ie)%c(:,:,chooselev,choosetrac,tl%n0))   
  end do
! 
  !need the buffer cellghostbuf in the time loop
  ! for mass calculation
  call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
  mass=global_shared_sum(1)
  maxc = parallelmax(tmp1,hybrid)
  minc = parallelmin(tmp2,hybrid)

! !-----------------------------------------------------------------------------------!
!   
  if (test_cfldep) then
    maxcflxstep = parallelmax(spelt(:)%maxcfl(1,chooselev),hybrid)
    maxcflystep = parallelmax(spelt(:)%maxcfl(2,chooselev),hybrid)
    maxcflx=max(maxcflxstep,maxcflx)
    maxcfly=max(maxcflystep,maxcfly)
    if  (hybrid%masterthread) then
      write(*,*) "CFL: maxcflxstep=", maxcflxstep, "maxcflystep=", maxcflystep 
    endif
  endif 
  if  (hybrid%masterthread) then 
    write(*,*) 'time=', time_at(tl%nstep), 'timeatmax',Time_at(nmax)
    write(*,*) 'chooselev=', chooselev, 'choosetrac=', choosetrac
    write(*,*) 'STEP',tl%nstep,'MAXSTEP',nmax, &
                't0', tl%n0, 't1', tl%np1
    write(*,*) 'massbegin', massstart, 'massend', mass 
    write(*,*) 'rel', (mass-massstart)/massstart           
    write(*,*) 'maxvaluestart:', maxcstart, 'minvaluestart:', mincstart
    write(*,*) 'maxvalue:     ', maxc,       'minvalue:    ', minc
    print *
  endif 
endif
!-----------------------------------------------------------------------------------!   
  !!! OUTPUT FOR TESTREASONS because I use gll points
  do ie=nets,nete
    do j=1-nhe,nc+nhe
      do i=1-nhe,nc+nhe
        icell=1+(i-1)*nipm
        jcell=1+(j-1)*nipm
        ff=spelt(ie)%c(icell:icell+nipm,jcell:jcell+nipm,chooselev,choosetrac,tl%n0)*spelt(ie)%sga(icell:icell+nipm,jcell:jcell+nipm)
        minmax(i,j,:)=cell_minmax(ff)
        call cip_coeff(spelt(ie)%drefx(i,j),spelt(ie)%drefy(i,j),ff,ff(2,2),cf(:,:,i,j))
      enddo
    enddo
    do j=1,np
      do i=1,np
        call cell_search(elem(ie),spelt(ie), elem(ie)%spherep(i,j), icell, jcell,dref, alphabeta,face_nodep)
        tmp=cip_interpolate(cf(:,:,icell,jcell),dref%x,dref%y)      
        tmp=qmsl_cell_filter(icell,jcell,minmax,tmp)
        sga=metric_term(alphabeta) 
!         sga=metric_termref(elem(ie),dref)
        elem(ie)%state%p(i,j,1,tl%n0)=g*tmp/sga
      enddo
    enddo 
    elem(ie)%state%ps(:,:)=0.0D0
  end do  

#ifdef PIO_INTERP
  call interp_movie_output(elem,tl, hybrid, 0D0, nets, nete,spelt)
#else
  call shal_movie_output(elem,tl, hybrid, 0D0, nets, nete,deriv)
#endif  

!-----------------------------------------------------------------------------------!
      
ENDDO  ! END TIME LOOP
  call t_stopf('spelt')

  call freeghostbuffertr(cellghostbuf)
  call freeghostbuffertr(factorR)
  call freeedgebuffer(edgeveloc)
#ifdef PIO_INTERP
    call interp_movie_finish
#else
    call shal_movie_finish
#endif  
!-----------------------------------------------------------------------------------!  
! Error analysis/ complicated, but for a first try o.k.
  do ie=nets,nete
    tmp=0.0D0
    tmpref=0.0D0
    global_shared_buf(ie,:)=0.0D0
    do j=1,nc
      do i=1,nc
        area=spelt(ie)%area_sphere(i,j)
        icell= 2+(i-1)*nipm
        jcell=2+(j-1)*nipm
        global_shared_buf(ie,1)=global_shared_buf(ie,1)+area*abs(spelt(ie)%c(i*nipm,j*nipm,chooselev,choosetrac,tl%n0)-spelt(ie)%cstart(i,j))
        global_shared_buf(ie,2)=global_shared_buf(ie,2)+area*abs(spelt(ie)%cstart(i,j))

        global_shared_buf(ie,3)=global_shared_buf(ie,3)+area*(spelt(ie)%c(i*nipm,j*nipm,chooselev,choosetrac,tl%n0)-spelt(ie)%cstart(i,j))* &
                                          (spelt(ie)%c(i*nipm,j*nipm,chooselev,choosetrac,tl%n0)-spelt(ie)%cstart(i,j))
        global_shared_buf(ie,4)=global_shared_buf(ie,4)+area*(spelt(ie)%cstart(i,j))*(spelt(ie)%cstart(i,j))
        tmp=max(tmp,abs(spelt(ie)%c(i*nipm,j*nipm,chooselev,choosetrac,tl%n0)-spelt(ie)%cstart(i,j)))
        tmpref=max(tmpref,abs(spelt(ie)%cstart(i,j)))
      end do
    end do
    tmp1(ie)=tmp
    tmp2(ie)=tmpref
  end do
  call wrap_repro_sum(nvars=4, comm=hybrid%par%comm)
  l1=global_shared_sum(1)/global_shared_sum(2)
  l2=sqrt(global_shared_sum(3)/global_shared_sum(4))

  lmax = parallelmax(tmp1,hybrid)/parallelmax(tmp2,hybrid)

  do ie=nets,nete
    spelt(ie)%elem_mass=0.0D0
    tmp1(ie)=-1.0D0-20
    tmp2(ie)=1.0D20
    do j=1,nc
      do i=1,nc
        area=spelt(ie)%area_sphere(i,j)
        icell= 2+(i-1)*nipm
        jcell=2+(j-1)*nipm
        if (choosetrac==1) then   ! mass of air, code is not optimal
         spelt(ie)%elem_mass=spelt(ie)%elem_mass + &
                        area*spelt(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0)
        elseif (choosetrac==2) then
          spelt(ie)%elem_mass=spelt(ie)%elem_mass + &
                         area*(spelt(ie)%c(icell,jcell,chooselev,1,tl%n0))*(spelt(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0))
        endif        
        tmp1(ie)=max(tmp1(ie),spelt(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0))
        tmp2(ie)=min(tmp2(ie),spelt(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0))
      enddo
    enddo
    ! for the mass value
    global_shared_buf(ie,1)=0.0D0
    global_shared_buf(ie,1)=spelt(ie)%elem_mass      
  end do
! 
  !need the buffer cellghostbuf in the time loop
  ! for mass calculation
  call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
  mass=global_shared_sum(1)
  maxc = parallelmax(tmp1,hybrid)
  minc = parallelmin(tmp2,hybrid)
!SUMMARY

  if(hybrid%masterthread) then 
    print *
    print *,"!-----------------------------------------------------------------------!"
    print *,"!  Test CASE for SPELT, Christoph Erath                                 !" 
    print *,"!-----------------------------------------------------------------------!"
    print *, "Summary"  
    write(*,*) 'number of elements', 6*ne*ne*nc*nc

    print *
    write(*,*) 'chooselev=', chooselev, 'choosetrac=', choosetrac
    write(*,*) 'massbegin', massstart, 'massend', mass 
    write(*,*) 'rel', (mass-massstart)/massstart
    write(*,*) 'maxvaluestart:', maxcstart, 'minvaluestart:', mincstart
    write(*,*) 'maxvalue:     ', maxc,      'minvalue:     ', minc
    write(*,*) "l1 = ", l1, "l2 = ", l2, "lmax = ", lmax
    write(*,*) "ne*nc = ", ne*nc, "timestep = ", tstep
    write(*,*) "CFL: maxcflx=", maxcflx, "maxcfly=", maxcfly 
    print *
  endif  
  

end subroutine spelt_run_bench

end module spelt_bench_mod
