#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module vertremap_mod
  use vertremap_base, only: remap1

  use kinds, only                  : real_kind,int_kind
  use dimensions_mod, only         : np,nlev,qsize,nlevp,npsq
  use hybvcoord_mod, only          : hvcoord_t
  use element_mod, only            : element_t
  use perf_mod, only               : t_startf, t_stopf  ! _EXTERNAL
  use parallel_mod, only           : abortmp, parallel_t
  use control_mod, only : vert_remap_q_alg

  implicit none
  private
  public :: vertical_remap

contains


  subroutine vertical_remap(hybrid,elem,hvcoord,dt,np1,np1_qdp,nets,nete)

  ! This routine is called at the end of the vertically Lagrangian
  ! dynamics step to compute the vertical flux needed to get back
  ! to reference eta levels
  !
  ! input:
  !     derived%dp()  delta p on levels at beginning of timestep
  !     state%dp3d(np1)  delta p on levels at end of timestep
  ! output:
  !     state%ps_v(np1)          surface pressure at time np1
  !     derived%eta_dot_dpdn()   vertical flux from final Lagrangian
  !                              levels to reference eta levels
  !
  use kinds,          only: real_kind
  use hybvcoord_mod,  only: hvcoord_t
  use control_mod,    only: rsplit
  use hybrid_mod,     only: hybrid_t
! HHLEE 20190521
#ifdef FIVE
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use physconst,      only: rair, gravit
  use ppgrid,         only: pver, pverp
  use five_intr,      only: pver_five, pverp_five, &
                            masswgt_vert_avg, linear_interp, &
                            hyai_five_toshare, hyam_five_toshare, &
                            hybi_five_toshare, hybm_five_toshare
#endif

  type (hybrid_t),  intent(in)    :: hybrid  ! distributed parallel structure (shared)
  type (element_t), intent(inout) :: elem(:)
  type (hvcoord_t)                :: hvcoord
  real (kind=real_kind)           :: dt

  integer :: ie,i,j,k,np1,nets,nete,np1_qdp
  integer :: q

  real (kind=real_kind), dimension(np,np,nlev)  :: dp,dp_star
  real (kind=real_kind), dimension(np,np,nlev,2)  :: ttmp
#ifdef FIVE
  real (kind=real_kind), dimension(np,np,pverp)  :: pint_host
  real (kind=real_kind), dimension(np,np,pverp)  :: pint_star_host
  real (kind=real_kind), dimension(np,np,pver)   :: pmid_host
  real (kind=real_kind), dimension(np,np,pver)   :: pmid_star_host
  real (kind=real_kind), dimension(np,np,pverp_five)  :: pint_five
  real (kind=real_kind), dimension(np,np,pver_five)   :: pmid_five
  real (kind=real_kind), dimension(np,np,pver_five)   :: dp_five
  real (kind=real_kind), dimension(np,np,pver_five)   :: dp_star_five
  real (kind=real_kind), dimension(np,np,pverp_five)  :: pint_star_five
  real (kind=real_kind), dimension(np,np,pver_five)   :: pmid_star_five
  real (kind=real_kind), dimension(np,np,pver_five)   :: dp3d_five

  real (kind=real_kind), dimension(np,np,pver_five)  :: ttmp_five
  real (kind=real_kind), dimension(np,np,pver_five)  :: utmp_five
  real (kind=real_kind), dimension(np,np,pver_five)  :: vtmp_five
  real (kind=real_kind), dimension(np,np,pver_five,qsize)  :: Qdptmp_five
  real (kind=real_kind), dimension(np,np,pver_five)  :: ttmp_five_old
  real (kind=real_kind), dimension(np,np,pver_five)  :: utmp_five_old
  real (kind=real_kind), dimension(np,np,pver_five)  :: vtmp_five_old
  real (kind=real_kind), dimension(np,np,pver_five,qsize)  :: Qdptmp_five_old
  real (kind=real_kind), dimension(np,np,pver_five)  :: ttmp_five_tend
  real (kind=real_kind), dimension(np,np,pver_five)  :: utmp_five_tend
  real (kind=real_kind), dimension(np,np,pver_five)  :: vtmp_five_tend
  real (kind=real_kind), dimension(np,np,pver_five,qsize)  :: Qdptmp_five_tend
  real (kind=real_kind), dimension(np,np,pver)       :: ttmp_host
  real (kind=real_kind), dimension(np,np,pver)       :: utmp_host
  real (kind=real_kind), dimension(np,np,pver)       :: vtmp_host
  real (kind=real_kind), dimension(np,np,pver,qsize) :: Qdptmp_host

  real (kind=real_kind), dimension(np,np,pver_five)  :: rho_five
  real (kind=real_kind), dimension(np,np,pver_five)  :: dz_five
  real (kind=real_kind), dimension(np,np,pver)       :: rho_host
  real (kind=real_kind), dimension(np,np,pver)       :: dz_host
#endif
  call t_startf('vertical_remap')

  ! reference levels:
  !   dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v(i,j)
  !   hybi(1)=0          pure pressure at top of atmosphere
  !   hyai(1)=ptop
  !   hyai(nlev+1) = 0   pure sigma at bottom
  !   hybi(nlev+1) = 1
  !
  ! sum over k=1,nlev
  !  sum(dp(k)) = (hyai(nlev+1)-hyai(1))*ps0 + (hybi(nlev+1)-hybi(1))*ps_v
  !             = -ps0 + ps_v
  !  ps_v =  ps0+sum(dp(k))
  !
  ! reference levels:
  !    dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v
  ! floating levels:
  !    dp_star(k) = dp(k) + dt_q*(eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  ! hence:
  !    (dp_star(k)-dp(k))/dt_q = (eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  !

#ifdef FIVE
  pint_host = 0.
  pmid_host = 0.
  pint_star_host = 0. 
  ttmp_five = 0.
  utmp_five = 0.
  vtmp_five = 0.
  Qdptmp_five = 0.
  ttmp_host = 0.
  utmp_host = 0.
  vtmp_host = 0.
  Qdptmp_host = 0.
#endif
 
   do ie=nets,nete
     ! update final ps_v
     elem(ie)%state%ps_v(:,:,np1) = hvcoord%hyai(1)*hvcoord%ps0 + &
          sum(elem(ie)%state%dp3d(:,:,:,np1),3)
     do k=1,nlev
        dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
        if (rsplit==0) then
           dp_star(:,:,k) = dp(:,:,k) + dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) -&
                elem(ie)%derived%eta_dot_dpdn(:,:,k))
        else
           dp_star(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)
        endif
#ifdef FIVE
        pint_host(:,:,k) = hvcoord%hyai(k)*hvcoord%ps0+hvcoord%hybi(k)*elem(ie)%state%ps_v(:,:,np1)
        pmid_host(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0+hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,np1)
        pint_star_host(:,:,1) = pint_host(:,:,1)
        pint_star_host(:,:,k+1) = pint_star_host(:,:,k) + dp_star(:,:,k)
        pmid_star_host(:,:,k) = (pint_star_host(:,:,k)+pint_star_host(:,:,k+1))/2.0
#endif
     enddo
#ifdef FIVE
        pint_host(:,:,pverp) = hvcoord%hyai(pverp)*hvcoord%ps0+hvcoord%hybi(pverp)*elem(ie)%state%ps_v(:,:,np1)
#endif
     if (minval(dp_star)<0) then
        do k=1,nlev
        do i=1,np
        do j=1,np
           if (dp_star(i,j,k ) < 0) then
              print *,'level = ',k
              print *,"column location lat,lon (radians):",elem(ie)%spherep(i,j)%lat,elem(ie)%spherep(i,j)%lon
           endif
        enddo
        enddo
        enddo
        call abortmp('negative layer thickness.  timestep or remap time too large')
     endif

#ifdef FIVE
     do k = 1, pverp_five
        pint_five(:,:,k) = hvcoord%ps0*hyai_five_toshare(k) + elem(ie)%state%ps_v(:,:,np1)*hybi_five_toshare(k)
     enddo
     do k=1,pver_five
        dp_five(:,:,k) = pint_five(:,:,k+1)-pint_five(:,:,k)
        pmid_five(:,:,k) = (pint_five(:,:,k)+pint_five(:,:,k+1))/2.0
     enddo

     do i=1,np
     do j=1,np
        call linear_interp(pint_host(i,j,:),pint_five(i,j,:),pint_star_host(i,j,1:pverp),pint_star_five(i,j,1:pverp_five),pverp,pverp_five)  
     enddo
     enddo

     do k=1,pver_five
        dp_star_five(:,:,k) = pint_star_five(:,:,k+1)-pint_star_five(:,:,k)
        pmid_star_five(:,:,k) = (pint_star_five(:,:,k)+pint_star_five(:,:,k+1))/2.0
     enddo
#endif

#ifdef FIVE
     if (rsplit>0) then

     do i=1,np
     do j=1,np
! get t, v, u, dp. In the future, they shoudl pass from pbuf, like T_five
        call linear_interp(pmid_star_host(i,j,:),pmid_star_five(i,j,:),elem(ie)%state%t(i,j,:,np1),ttmp_five(i,j,:),pver,pver_five)
        call linear_interp(pmid_star_host(i,j,:),pmid_star_five(i,j,:),elem(ie)%state%v(i,j,1,:,np1),utmp_five(i,j,:),pver,pver_five)
        call linear_interp(pmid_star_host(i,j,:),pmid_star_five(i,j,:),elem(ie)%state%v(i,j,2,:,np1),vtmp_five(i,j,:),pver,pver_five)
       do q = 1, qsize 
         call linear_interp(pmid_star_host(i,j,:),pmid_star_five(i,j,:),elem(ie)%state%Qdp(i,j,:,q,np1_qdp),Qdptmp_five(i,j,:,q),pver,pver_five)
         Qdptmp_five_old(i,j,:,q) = Qdptmp_five(i,j,:,q)
       enddo 
     enddo 
     enddo

        ttmp_five_old = ttmp_five 
        ttmp_five = ttmp_five*dp_star_five
        call t_startf('vertical_remap1_FIVE_1')
        call remap1(ttmp_five,np,1,dp_star_five,dp_five)
        call t_stopf('vertical_remap1_FIVE_1')
        ttmp_five = ttmp_five/dp_five
        ttmp_five_tend = (ttmp_five - ttmp_five_old)/dt

        utmp_five_old = utmp_five 
        utmp_five = utmp_five*dp_star_five
        call t_startf('vertical_remap1_FIVE_2')
        call remap1(utmp_five,np,1,dp_star_five,dp_five)
        call t_stopf('vertical_remap1_FIVE_2')
        utmp_five = utmp_five/dp_five
        utmp_five = utmp_five/dp_five
        utmp_five_tend = (utmp_five - utmp_five_old)/dt

        vtmp_five_old = vtmp_five 
        vtmp_five = vtmp_five*dp_star_five
        call t_startf('vertical_remap1_FIVE_3')
        call remap1(vtmp_five,np,1,dp_star_five,dp_five)
        call t_stopf('vertical_remap1_FIVE_3')
        vtmp_five = vtmp_five/dp_five
        vtmp_five = vtmp_five/dp_five
        vtmp_five_tend = (vtmp_five - vtmp_five_old)/dt

     if (qsize>0) then

       call t_startf('vertical_remap1_FIVE_4')
       call remap1(Qdptmp_five,np,qsize,dp_star_five,dp_five)
       call t_stopf('vertical_remap1_FIVE_4')

     endif
! ----- For rho_host esimate on the reference levels
        ttmp(:,:,:,1)=elem(ie)%state%t(:,:,:,np1)
        ttmp(:,:,:,1)=ttmp(:,:,:,1)*dp_star

        call t_startf('vertical_remap1_1')
        call remap1(ttmp,np,1,dp_star,dp)
        call t_stopf('vertical_remap1_1')

!-------------------
! should calculate rho_five before or after remap??
     rho_five = pmid_five/(rair*ttmp_five)
     rho_host = pmid_host/(rair*ttmp(:,:,:,1)/dp)
     dz_five = dp_five/(rho_five*gravit)
     dz_host = dp/(rho_host*gravit)

     do i=1,np
     do j=1,np
      ! Mass weighted vertical average for temperature

        call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),dz_host(i,j,:),dz_five(i,j,:),&
                            pint_host(i,j,:),pmid_five(i,j,:),pmid_host(i,j,:),&
                            ttmp_five_tend(i,j,:),ttmp_host(i,j,:))

        call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),dz_host(i,j,:),dz_five(i,j,:),&
                            pint_host(i,j,:),pmid_five(i,j,:),pmid_host(i,j,:),&
                            utmp_five_tend(i,j,:),utmp_host(i,j,:))

        call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),dz_host(i,j,:),dz_five(i,j,:),&
                            pint_host(i,j,:),pmid_five(i,j,:),pmid_host(i,j,:),&
                            vtmp_five_tend(i,j,:),vtmp_host(i,j,:))

       do q = 1, qsize 
        Qdptmp_five_tend(i,j,:,q) = (Qdptmp_five(i,j,:,q) - Qdptmp_five_old(i,j,:,q))/dt

        call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),dz_host(i,j,:),dz_five(i,j,:),&
                            pint_host(i,j,:),pmid_five(i,j,:),pmid_host(i,j,:),&
                            Qdptmp_five_tend(i,j,:,q),Qdptmp_host(i,j,:,q))
       enddo 
     enddo 
     enddo
! update t, u, v, qdp
     elem(ie)%state%t(:,:,:,np1)=elem(ie)%state%t(:,:,:,np1)+ttmp_host(:,:,:)*dt
     elem(ie)%state%v(:,:,1,:,np1)=elem(ie)%state%v(:,:,1,:,np1)+utmp_host(:,:,:)*dt
     elem(ie)%state%v(:,:,2,:,np1)=elem(ie)%state%v(:,:,2,:,np1)+vtmp_host(:,:,:)*dt

    end if
#else

     if (rsplit>0) then
        !  REMAP u,v,T from levels in dp3d() to REF levels
#undef REMAP_TE
#ifdef REMAP_TE
        ! remap u,v and cp*T + .5 u^2
        ttmp(:,:,:,1)=(elem(ie)%state%v(:,:,1,:,np1)**2 + &
             elem(ie)%state%v(:,:,2,:,np1)**2)/2 + &
             elem(ie)%state%t(:,:,:,np1)*cp
#else
        ttmp(:,:,:,1)=elem(ie)%state%t(:,:,:,np1)
#endif
        ttmp(:,:,:,1)=ttmp(:,:,:,1)*dp_star

        call t_startf('vertical_remap1_1')
        call remap1(ttmp,np,1,dp_star,dp)
        call t_stopf('vertical_remap1_1')

        elem(ie)%state%t(:,:,:,np1)=ttmp(:,:,:,1)/dp

        ttmp(:,:,:,1)=elem(ie)%state%v(:,:,1,:,np1)*dp_star
        ttmp(:,:,:,2)=elem(ie)%state%v(:,:,2,:,np1)*dp_star

        call t_startf('vertical_remap1_2')
        call remap1(ttmp,np,2,dp_star,dp)
        call t_stopf('vertical_remap1_2')

        elem(ie)%state%v(:,:,1,:,np1)=ttmp(:,:,:,1)/dp
        elem(ie)%state%v(:,:,2,:,np1)=ttmp(:,:,:,2)/dp

#ifdef REMAP_TE
        ! back out T from TE
        elem(ie)%state%t(:,:,:,np1) = &
             ( elem(ie)%state%t(:,:,:,np1) - ( (elem(ie)%state%v(:,:,1,:,np1)**2+ &
             elem(ie)%state%v(:,:,2,:,np1)**2)/2))/cp
#endif

     endif
#endif

     ! remap the gll tracers from lagrangian levels (dp_star)  to REF levels dp
     if (qsize>0) then
#ifdef FIVE
      do q = 1, qsize
         elem(ie)%state%Qdp(:,:,:,q,np1_qdp) = elem(ie)%state%Qdp(:,:,:,q,np1_qdp) + &
                   Qdptmp_host(:,:,:,q) *dt
      enddo
#else
       call t_startf('vertical_remap1_3')
       call remap1(elem(ie)%state%Qdp(:,:,:,:,np1_qdp),np,qsize,dp_star,dp)
       call t_stopf('vertical_remap1_3')
#endif

     endif

  enddo
  call t_stopf('vertical_remap')
  end subroutine vertical_remap


end module 




