!-----------------------------------------------------------------------------------!
!MODULE CWFV_MOD-------------------------------------------------------CE-for CWFV!
! CWFV_MOD File for the cwfv project in HOMME                                      !
! Author: Christoph Erath                                                           !
! Date: 13.January 2012                                                             !
! MAIN module to run CWFV on HOMME                                                  !
!-----------------------------------------------------------------------------------!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module cwfv_bench_mod
  use kinds, only : real_kind, int_kind
contains

subroutine cwfv_run_bench(elem,fvm,hybrid,nets,nete,tl)
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
  use cwfv_mod, only: cellghostbuf, edgeveloc, cwfv_struct
  use cwfv_mod, only: fvm_mcgregordsscwfv, spelt_run 
  use cwfv_mod, only: cip_coeff, cip_interpolate, metric_term, metric_termref, cell_search, qmsl_cell_filter, cell_minmax, cip_cell_avr
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
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------
  use quadrature_mod, only : quadrature_t, gausslobatto
  
  
#ifdef PIO_INTERP
    use interp_movie_mod, only : interp_movie_init, interp_movie_output, interp_movie_finish
#else
    use shal_movie_mod, only : shal_movie_init, shal_movie_output, shal_movie_finish
#endif
  
  implicit none
  type (element_t), intent(inout)             :: elem(:)
  type (cwfv_struct), intent(inout)           :: fvm(:)
  
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
  
  real (kind=real_kind)                       :: refnep(1:nep)
  
  integer                                     :: icell, jcell
  real (kind=real_kind)                       :: sga
  type (cartesian2D_t)                        :: dref, alphabeta
  real (kind=real_kind)                       :: vstar(np,np,2), v1, v2
  
  type (quadrature_t) :: gp   ! Quadrature points and weights on pressure grid
  
  integer  choosetrac, chooselev   !for test reason the output
 !-----------------------------------------------------------------------------------!  
 choosetrac=1
 chooselev=1
 
  if(hybrid%masterthread) then 
    print *,"!-----------------------------------------------------------------------!"
    print *,"!  Test CASE for CWFV, Christoph Erath                                  !" 
    print *,"!-----------------------------------------------------------------------!" 
  endif
     
  ! Initialize derivative structure
  ! CWFV nodes are equally spaced in alpha/beta
  ! HOMME with equ-angular gnomonic projection maps alpha/beta space
  ! to the reference element via simple scale + translation
  ! thus, CWFV nodes in reference element [-1,1] are a tensor product of
  ! array 'fvm_nodes(:)' computed below:
  call derivinit(deriv)
  tmp=nep-1
  do i=1,nep
    refnep(i)= 2*(i-1)/tmp - 1
  end do
  gp=gausslobatto(np)
  call v2pinit(deriv%Sfvm,gp%points,refnep,np,nep)
  
!Initialize test example with first communication-----------------------------------!      
  do ie=nets,nete
    do k=1, nlev
!         fvm(ie)%c(:,:,k,1,tl%n0)=1.0D0    !density of the air
      do itr=1,ntrac
        do j=1,nep
          do i=1,nep  
        ! define test example
            call analytical_function(fvm(ie)%c(i,j,k,itr,tl%n0),fvm(ie)%asphere(i,j),k,itr)
            sga=fvm(ie)%sga(i,j)
            fvm(ie)%c(i,j,k,itr,tl%n0)=sga*fvm(ie)%c(i,j,k,itr,tl%n0)  
            if ((mod(i,2)==0) .and. (mod(j,2)==0) .and. (k==chooselev) .and. (itr==choosetrac)) then
              fvm(ie)%cstart(i/2,j/2)=fvm(ie)%c(i,j,k,itr,tl%n0)  
            endif  
          end do
        end do
      end do 
    end do
    !first exchange of the initial values
    call ghostVpack2d(cellghostbuf,fvm(ie)%c,nipm, nep,nlev,ntrac,0, tl%n0, timelevels,elem(ie)%desc)
  end do
!-----------------------------------------------------------------------------------!  
  call ghost_exchangeV(hybrid,cellghostbuf,nipm,nep)
!-----------------------------------------------------------------------------------!    
  do ie=nets,nete
    call ghostVunpack2d(cellghostbuf,fvm(ie)%c,nipm, nep,nlev,ntrac,0, tl%n0, timelevels,elem(ie)%desc)
  end do
!-----------------------------------------------------------------------------------!
! initialize test example
  do ie=nets,nete  
    do k=1, nlev
      do itr=1,ntrac
        do j=1-nhe,nc+nhe
          do i=1-nhe,nc+nhe
            icell=1+(i-1)*nipm
            jcell=1+(j-1)*nipm
            ff=fvm(ie)%c(icell:icell+nipm,jcell:jcell+nipm,k,itr,tl%n0)
            ff(2,2)=cip_cell_avr(ff) 
            fvm(ie)%c(icell+1,jcell+1,k,itr,tl%n0)=ff(2,2)
            ! FOR TEST OUTPUT on gll only needed
            minmax(i,j,:)=cell_minmax(ff)
            call cip_coeff(ff,ff(2,2),cf(:,:,i,j))
          enddo
        enddo
        if ((itr==choosetrac) .and. (k==chooselev)) then
          do j=1,np    ! FOR TEST OUTPUT on gll only needed
            do i=1,np
              call cell_search(elem(ie), elem(ie)%spherep(i,j), icell, jcell,dref, alphabeta)
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
  call interp_movie_output(elem,tl, hybrid, 0D0, deriv, nets, nete)
#else
  call shal_movie_init(elem,hybrid)
  call shal_movie_output(elem,tl, hybrid, 0D0, nets, nete,deriv)
#endif
!-----------------------------------------------------------------------------------! 
! for Quantities as mass aso 
  do ie=nets,nete
   fvm(ie)%elem_mass=0
   tmp=abs(elem(ie)%corners(1)%x-elem(ie)%corners(2)%x)/nc
   do j=1,nc
     do i=1,nc
       icell= 2+(i-1)*nipm
       jcell=2+(j-1)*nipm
       if (choosetrac==1) then   ! mass of air, code is not optimal
         fvm(ie)%elem_mass=fvm(ie)%elem_mass + &
                       tmp*tmp*fvm(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0)
       else
         fvm(ie)%elem_mass=fvm(ie)%elem_mass + &
                       tmp*tmp**fvm(ie)%c(icell,jcell,chooselev,1,tl%n0)*&
                                                  fvm(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0)
       endif
  !           fvm(ie)%cstart(i,j)=fvm(ie)%c(i,j,chooselev,choosetrac,tl%n0)
     enddo
   enddo
   ! for the mass value
   global_shared_buf(ie,1)=0D0
   global_shared_buf(ie,1)=fvm(ie)%elem_mass
   ! for the max value on the sphere
   tmp1(ie) = MAXVAL(fvm(ie)%c(:,:,chooselev,choosetrac,tl%n0))
   tmp2(ie) = MINVAL(fvm(ie)%c(:,:,chooselev,choosetrac,tl%n0))   
  end do 


  !need the buffer cellghostbuf in the time loop
  ! for mass calculation
  call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
  massstart=global_shared_sum(1)
  maxcstart = parallelmax(tmp1,hybrid)
  mincstart = parallelmin(tmp2,hybrid) 
!-----------------------------------------------------------------------------------!
 
DO WHILE(tl%nstep<nmax)
! ! start mcgregordss
  do ie=nets,nete
    do k=1,nlev
      elem(ie)%derived%vstar(:,:,:,k)=get_boomerang_velocities_gll(elem(ie), time_at(tl%nstep+1))
      fvm(ie)%vn0(:,:,:,k)=get_boomerang_velocities_gll(elem(ie),time_at(tl%nstep))
!       elem(ie)%derived%vstar(:,:,:,k)=get_solidbody_velocities_gll(elem(ie), time_at(tl%nstep+1))
!       fvm(ie)%vn0(:,:,:,k)=get_solidbody_velocities_gll(elem(ie),time_at(tl%nstep))
      vstar=elem(ie)%derived%vstar(:,:,:,k)/rearth
      do j=1,np
        do i=1,np
!           v1 = elem(ie)%Dinv(1,1,i,j)*vstar(i,j,1) + elem(ie)%Dinv(1,2,i,j)*vstar(i,j,2)
!           v2 = elem(ie)%Dinv(2,1,i,j)*vstar(i,j,1) + elem(ie)%Dinv(2,2,i,j)*vstar(i,j,2)
          v1 = fvm(ie)%Dinv(1,1,i,j)*vstar(i,j,1) + fvm(ie)%Dinv(1,2,i,j)*vstar(i,j,2)
          v2 = fvm(ie)%Dinv(2,1,i,j)*vstar(i,j,1) + fvm(ie)%Dinv(2,2,i,j)*vstar(i,j,2)
          vstar(i,j,1)=v1
          vstar(i,j,2)=v2
        enddo
      enddo
      fvm(ie)%contrau=interpolate_gll2spelt_points(vstar(:,:,1),deriv)
      fvm(ie)%contrav=interpolate_gll2spelt_points(vstar(:,:,2),deriv)
    end do
  end do
  call fvm_mcgregordsscwfv(elem,fvm,nets,nete, hybrid, deriv, tstep, 3)
! ! end mcgregordss

  call spelt_run(elem,fvm,hybrid,deriv,tstep,tl,nets,nete)

  call TimeLevel_update(tl,"forward") 
!-----------------------------------------------------------------------------------! 
! for Quantities as mass aso
  do ie=nets,nete
    fvm(ie)%elem_mass=0
    tmp=abs(elem(ie)%corners(1)%x-elem(ie)%corners(2)%x)/nc
    do j=1,nc
      do i=1,nc
        icell= 2+(i-1)*nipm
        jcell=2+(j-1)*nipm
        if (choosetrac==1) then   ! mass of air, code is not optimal
          fvm(ie)%elem_mass=fvm(ie)%elem_mass + &
                         tmp*tmp*fvm(ie)%c(icell,jcell,chooselev,choosetrac,tl%n0)
        else
          fvm(ie)%elem_mass=fvm(ie)%elem_mass + &
                         tmp*tmp**fvm(ie)%c(icell,jcell,chooselev,1,tl%n0)*&
                                                    fvm(ie)%c(i,j,chooselev,choosetrac,tl%n0)
        endif
 !           fvm(ie)%cstart(i,j)=fvm(ie)%c(i,j,chooselev,choosetrac,tl%n0)
      enddo
    enddo
    ! for the mass value
    global_shared_buf(ie,1)=0D0
    global_shared_buf(ie,1)=fvm(ie)%elem_mass
    ! for the max value on the sphere
    tmp1(ie) = MAXVAL(fvm(ie)%c(:,:,chooselev,choosetrac,tl%n0))
    tmp2(ie) = MINVAL(fvm(ie)%c(:,:,chooselev,choosetrac,tl%n0))   
  end do

  !need the buffer cellghostbuf in the time loop
  ! for mass calculation
  call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
  mass=global_shared_sum(1)
  maxc = parallelmax(tmp1,hybrid)
  minc = parallelmin(tmp2,hybrid)
!-----------------------------------------------------------------------------------!
   
  if  (hybrid%masterthread) then 
    write(*,*) 'time=', time_at(tl%nstep), 'timeatmax',Time_at(nmax)
    write(*,*) 'chooselev=', chooselev, 'choosetrac=', choosetrac
    write(*,*) 'STEP',tl%nstep,'MAXSTEP',nmax, &
                't0', tl%n0, 't1', tl%np1
    write(*,*) 'massbegin', massstart, 'massend', mass 
    write(*,*) 'rel', (mass-massstart)/massstart           
    write(*,*) 'maxvaluestart:', maxcstart, 'minvaluestart:', mincstart
    write(*,*) 'maxvalue:     ', maxc,       'minvalue:    ', minc
!      write(*,*) "CFL: maxcflx=", maxcflx, "maxcfly=", maxcfly 
    print *
  endif 
!-----------------------------------------------------------------------------------!   
  !!!! OUTPUT FOR TESTREASONS because I use gll points
  do ie=nets,nete
    do j=1-nhe,nc+nhe
      do i=1-nhe,nc+nhe
        icell=1+(i-1)*nipm
        jcell=1+(j-1)*nipm
        ff=fvm(ie)%c(icell:icell+nipm,jcell:jcell+nipm,chooselev,choosetrac,tl%n0)
        minmax(i,j,:)=cell_minmax(ff)
        call cip_coeff(ff,ff(2,2),cf(:,:,i,j))
      enddo
    enddo
    do j=1,np
      do i=1,np
        call cell_search(elem(ie), elem(ie)%spherep(i,j), icell, jcell,dref, alphabeta)
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
  call interp_movie_output(elem,tl, hybrid, 0D0, deriv, nets, nete)
#else
  call shal_movie_output(elem,tl, hybrid, 0D0, nets, nete,deriv)
#endif  
!-----------------------------------------------------------------------------------!
      
ENDDO  ! END TIME LOOP
  call freeghostbuffertr(cellghostbuf)
  call freeedgebuffer(edgeveloc)
  
!-----------------------------------------------------------------------------------!  
! Error analysis/ complicated, but for a first try o.k.
  do ie=nets,nete
    tmp=0.0D0
    tmpref=0.0D0
    global_shared_buf(ie,:)=0.0D0
    area=abs(elem(ie)%corners(1)%x-elem(ie)%corners(2)%x)/nc
    area=area*area
    do j=1,nc
      do i=1,nc
        global_shared_buf(ie,1)=global_shared_buf(ie,1)+area*abs(fvm(ie)%c(i*nipm,j*nipm,chooselev,choosetrac,tl%n0)-fvm(ie)%cstart(i,j))
        global_shared_buf(ie,2)=global_shared_buf(ie,2)+area*abs(fvm(ie)%cstart(i,j))

        global_shared_buf(ie,3)=global_shared_buf(ie,3)+area*(fvm(ie)%c(i*nipm,j*nipm,chooselev,choosetrac,tl%n0)-fvm(ie)%cstart(i,j))* &
                                          (fvm(ie)%c(i*nipm,j*nipm,chooselev,choosetrac,tl%n0)-fvm(ie)%cstart(i,j))
        global_shared_buf(ie,4)=global_shared_buf(ie,4)+area*(fvm(ie)%cstart(i,j))*(fvm(ie)%cstart(i,j))

        tmp=max(tmp,abs(fvm(ie)%c(i,j,chooselev,choosetrac,tl%n0)-fvm(ie)%cstart(i,j)))
        tmpref=max(tmpref,abs(fvm(ie)%cstart(i,j)))
      end do
    end do
    tmp1(ie)=tmp
    tmp2(ie)=tmpref
  end do
  call wrap_repro_sum(nvars=4, comm=hybrid%par%comm)
  l1=global_shared_sum(1)/global_shared_sum(2)
  l2=sqrt(global_shared_sum(3)/global_shared_sum(4))

  lmax = parallelmax(tmp1,hybrid)/parallelmax(tmp2,hybrid)


!SUMMARY
  if(hybrid%masterthread) then 
    print *
    print *,"!-----------------------------------------------------------------------!"
    print *,"!  Test CASE for FVM, Christoph Erath                                   !" 
    print *,"!-----------------------------------------------------------------------!"
    print *  
    write(*,*) 'number of elements', 6*ne*ne*nc*nc

    print *
    write(*,*) 'chooselev=', chooselev, 'choosetrac=', choosetrac
    write(*,*) 'massbegin', massstart, 'massend', mass 
    write(*,*) 'rel', (mass-massstart)/massstart
    write(*,*) 'maxvaluestart:', maxcstart, 'minvaluestart:', mincstart
    write(*,*) 'maxvalue:     ', maxc,      'minvalue:     ', minc
    write(*,*) "l1 = ", l1, "l2 = ", l2, "lmax = ", lmax
    write(*,*) "ne*nc = ", ne*nc, "timestep = ", tstep
  endif  
  

end subroutine cwfv_run_bench

end module cwfv_bench_mod
