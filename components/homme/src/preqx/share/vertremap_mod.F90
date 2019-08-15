#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module vertremap_mod
  use vertremap_base, only: remap1

  use kinds, only                  : real_kind,int_kind
  use dimensions_mod, only         : np,nlev,qsize,nlevp,npsq, nelemd
  use hybvcoord_mod, only          : hvcoord_t
  use element_mod, only            : element_t
  use perf_mod, only               : t_startf, t_stopf  ! _EXTERNAL
  use parallel_mod, only           : abortmp, parallel_t
  use control_mod, only : vert_remap_q_alg
#ifdef FIVE
  use physics_buffer, only: physics_buffer_desc
#endif
  implicit none
  private
  public :: vertical_remap

contains


  subroutine vertical_remap(hybrid,elem,hvcoord,dt,np1,np1_qdp,nets,nete,single_column &
#ifdef FIVE
             , t_five, q_five, u_five, v_five &
#endif
  )

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
  use ppgrid,         only: pver, pverp
! HHLEE 20190521
#ifdef FIVE
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use physconst,      only: rair, gravit
  use five_intr,      only: pver_five, pverp_five, &
                            masswgt_vert_avg, linear_interp, &
                            tendency_low_to_high, compute_five_heights,&
                            hyai_five_toshare, hyam_five_toshare, &
                            hybi_five_toshare, hybm_five_toshare
#endif

  logical,          intent(in)    :: single_column
  type (hybrid_t),  intent(in)    :: hybrid  ! distributed parallel structure (shared)
  type (element_t), intent(inout) :: elem(:)
  type (hvcoord_t)                :: hvcoord
  real (kind=real_kind)           :: dt

  integer :: ie,i,j,k,np1,nets,nete,np1_qdp
  integer :: q, npt

  real (kind=real_kind), dimension(np,np,nlev)  :: dp,dp_star
  real (kind=real_kind), dimension(np,np,nlev,2)  :: ttmp
  real (kind=real_kind), dimension(np,np,pverp)   :: eta_dot_dpdn
#ifdef FIVE
  real(r8),  intent(in) :: t_five(np,np,pver_five,nelemd)
  real(r8),  intent(in) :: u_five(np,np,pver_five,nelemd)
  real(r8),  intent(in) :: v_five(np,np,pver_five,nelemd)
  real(r8),  intent(in) :: q_five(np,np,pver_five,qsize,nelemd)
 
  real(r8) :: t_five_low(np,np,pver,nelemd)
  real(r8) :: u_five_low(np,np,pver,nelemd)
  real(r8) :: v_five_low(np,np,pver,nelemd)
  real(r8) :: q_five_low(np,np,pver,qsize,nelemd)

  real(r8) :: t_five_tend_low(np,np,pver,nelemd)
  real(r8) :: u_five_tend_low(np,np,pver,nelemd)
  real(r8) :: v_five_tend_low(np,np,pver,nelemd)
  real(r8) :: q_five_tend_low(np,np,pver,qsize,nelemd)

  real(r8) :: t_five_tend(np,np,pver_five,nelemd)
  real(r8) :: u_five_tend(np,np,pver_five,nelemd)
  real(r8) :: v_five_tend(np,np,pver_five,nelemd)
  real(r8) :: q_five_tend(np,np,pver_five,qsize,nelemd)

  real (kind=real_kind), dimension(np,np,pverp)  :: pint_host
  real (kind=real_kind), dimension(np,np,pverp)  :: pint_star
  real (kind=real_kind), dimension(np,np,pver)   :: pmid_host
  real (kind=real_kind), dimension(np,np,pver)   :: pdel_host
  real (kind=real_kind), dimension(np,np,pverp_five)  :: pint_five
  real (kind=real_kind), dimension(np,np,pver_five)   :: pmid_five
  real (kind=real_kind), dimension(np,np,pver_five)   :: pdel_five
  real (kind=real_kind), dimension(np,np,pver_five)   :: dp_five
  real (kind=real_kind), dimension(np,np,pver_five)   :: dp_star_five
  real (kind=real_kind), dimension(np,np,pverp_five)  :: pint_star_five
  real (kind=real_kind), dimension(np,np,pver,qsize)  :: state_Q
  real (kind=real_kind), dimension(np,np,pver,qsize)  :: state_Qnew

  real (kind=real_kind), dimension(np,np,pver_five)  :: ttmp_five
  real (kind=real_kind), dimension(np,np,pver_five)  :: utmp_five
  real (kind=real_kind), dimension(np,np,pver_five)  :: vtmp_five
  real (kind=real_kind), dimension(np,np,pver_five,qsize)  :: Qtmp_five
  real (kind=real_kind), dimension(np,np,pver_five)  :: ttmp_five_old
  real (kind=real_kind), dimension(np,np,pver_five)  :: utmp_five_old
  real (kind=real_kind), dimension(np,np,pver_five)  :: vtmp_five_old
  real (kind=real_kind), dimension(np,np,pver_five,qsize)  :: Qtmp_five_old
  real (kind=real_kind), dimension(np,np,pver_five)  :: ttmp_five_tend
  real (kind=real_kind), dimension(np,np,pver_five)  :: utmp_five_tend
  real (kind=real_kind), dimension(np,np,pver_five)  :: vtmp_five_tend
  real (kind=real_kind), dimension(np,np,pver_five,qsize)  :: Qtmp_five_tend
  real (kind=real_kind), dimension(np,np,pver)       :: ttmp_host
  real (kind=real_kind), dimension(np,np,pver)       :: utmp_host
  real (kind=real_kind), dimension(np,np,pver)       :: vtmp_host
  real (kind=real_kind), dimension(np,np,pver,qsize) :: Qtmp_host

  real (kind=real_kind), dimension(np,np,pver_five)  :: rho_five
  real (kind=real_kind), dimension(np,np,pver_five)  :: dz_five
  real (kind=real_kind), dimension(np,np,pver)       :: rho_host
  real (kind=real_kind), dimension(np,np,pver)       :: dz_host
  real (kind=real_kind), dimension(np,np,pver)       :: zm_host
  real (kind=real_kind), dimension(np,np,pverp)      :: zi_host
  real (kind=real_kind), dimension(np,np,pver_five)        :: zm_five
  real (kind=real_kind), dimension(np,np,pverp_five)       :: zi_five

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
  ttmp_five = 0.
  utmp_five = 0.
  vtmp_five = 0.
  Qtmp_five = 0.
  ttmp_host = 0.
  utmp_host = 0.
  vtmp_host = 0.
  Qtmp_host = 0.
#endif
 
   do ie = nets, nete
     ! update final ps_v
     elem(ie)%state%ps_v(:,:,np1) = hvcoord%hyai(1)*hvcoord%ps0 + &
          sum(elem(ie)%state%dp3d(:,:,:,np1),3)

     if (single_column) then
        do k=1,nlev
           eta_dot_dpdn(:,:,k)=elem(1)%derived%omega_p(1,1,k)
        enddo
           eta_dot_dpdn(:,:,nlev+1) = eta_dot_dpdn(:,:,nlev)
     endif


     do k=1,nlev
        dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)


        if (single_column) then
           dp_star(:,:,k) = dp(:,:,k) + dt*(eta_dot_dpdn(:,:,k+1) - eta_dot_dpdn(:,:,k))
        else
          if (rsplit==0) then
             dp_star(:,:,k) = dp(:,:,k) + dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) -&
                  elem(ie)%derived%eta_dot_dpdn(:,:,k))
          else
             dp_star(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)
          endif
        endif 
 
#ifdef FIVE
        pint_host(:,:,k) = hvcoord%hyai(k)*hvcoord%ps0+hvcoord%hybi(k)*elem(ie)%state%ps_v(:,:,np1)
        pmid_host(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0+hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,np1)
#endif

     enddo
#ifdef FIVE
        pint_host(:,:,pverp) = hvcoord%hyai(pverp)*hvcoord%ps0+hvcoord%hybi(pverp)*elem(ie)%state%ps_v(:,:,np1)

        pint_star(:,:,1) = pint_host(:,:,1)

      do k = 1, nlev
         pint_star(:,:,k+1) = pint_star(:,:,k) + dp_star(:,:,k)
         pdel_host(:,:,k) = pint_host(:,:,k+1)-pint_host(:,:,k)
      enddo
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
        pdel_five(:,:,k) = pint_five(:,:,k+1)-pint_five(:,:,k)
   enddo

     do i=1,np
     do j=1,np
        !call linear_interp(pint_host(i,j,:),pint_five(i,j,:),eta_dot_dpdn(i,j,1:pverp),&
        !             eta_dot_dpdn_five(i,j,1:pverp_five),pverp,pverp_five)  
        call linear_interp(pint_host(i,j,:),pint_five(i,j,:),pint_star(i,j,1:pverp),&
                     pint_star_five(i,j,1:pverp_five),pverp,pverp_five)  
     enddo
     enddo

     do k=1,pver_five
        !dp_star_five(:,:,k) = dp_five(:,:,k) + dt*(eta_dot_dpdn_five(:,:,k+1) - eta_dot_dpdn_five(:,:,k))
        dp_star_five(:,:,k) = pint_star_five(:,:,k+1) - pint_star_five(:,:,k)
     enddo

! ----- syncronize FIVE variables
     do i=1,np
     do j=1,np
        call linear_interp(pmid_host(i,j,:),pmid_five(i,j,:),elem(ie)%state%t(i,j,:,np1),ttmp_five(i,j,:),pver,pver_five)
     enddo
     enddo

     rho_five = pmid_five/(rair*ttmp_five)
     rho_host = pmid_host/(rair*elem(ie)%state%t(:,:,:,np1))
     dz_five = dp_five/(rho_five*gravit)
     dz_host = dp/(rho_host*gravit)

    do i=1,np
    do j=1,np
      ! Mass weighted vertical average for temperature
      call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),dz_host(i,j,:),dz_five(i,j,:),&
                            pint_host(i,j,:),pmid_five(i,j,:),pmid_host(i,j,:),&
                            t_five(i,j,:,ie),t_five_low(i,j,:,ie))
                        
      ! Mass weighted vertical average for u wind
      call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),dz_host(i,j,:),dz_five(i,j,:),&
                            pint_host(i,j,:),pmid_five(i,j,:),pmid_host(i,j,:),&
                            u_five(i,j,:,ie),u_five_low(i,j,:,ie))
                        
      ! Mass weighted vertical average for v wind
      call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),dz_host(i,j,:),dz_five(i,j,:),&
                            pint_host(i,j,:),pmid_five(i,j,:),pmid_host(i,j,:),&
                            v_five(i,j,:,ie),v_five_low(i,j,:,ie))

      ! Mass weighted vertical average for tracers 
      do q = 1, qsize 
        call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),dz_host(i,j,:),dz_five(i,j,:),&
                            pint_host(i,j,:),pmid_five(i,j,:),pmid_host(i,j,:),&
                            q_five(i,j,:,q,ie),q_five_low(i,j,:,q,ie))
      enddo

    enddo
    enddo

    ! Next compute the tendency of FIVE variables from the state, this is 
    !   done on the E3SM grid
    do k=1,pver
      do i=1,np
        do j=1,np
        t_five_tend_low(i,j,k,ie) = (elem(ie)%state%t(i,j,k,np1) - t_five_low(i,j,k,ie))/dt
        u_five_tend_low(i,j,k,ie) = (elem(ie)%state%v(i,j,1,k,np1) - u_five_low(i,j,k,ie))/dt
        v_five_tend_low(i,j,k,ie) = (elem(ie)%state%v(i,j,2,k,np1) - v_five_low(i,j,k,ie))/dt

        do q=1,qsize
          state_Q(i,j,k,q) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp) / dp_star(i,j,k) 
          q_five_tend_low(i,j,k,q,ie) = (state_Q(i,j,k,q) - q_five_low(i,j,k,q,ie))/dt
        enddo

        enddo
      enddo
    enddo

    ! Now interpolate this tendency to the higher resolution FIVE grid, 
    !   using the interpolation method of Sheng and Zwiers (1998), 
    !   as documented in Yamaguchi et al. (2017) Appendix B
    do i=1,np
    do j=1,np
      call compute_five_heights(pmid_host(i,j,1:pver),pint_host(i,j,1:pver+1), &
             elem(ie)%state%t(i,j,1:pver,np1),&
             state_Q(i,j,1:pver,1),state_Q(i,j,1:pver,2),pdel_host(i,j,1:pver),pver,hvcoord%ps0,&
             zm_host(i,j,1:pver),zi_host(i,j,1:pver+1))

      call compute_five_heights(pmid_five(i,j,:),pint_five(i,j,:),t_five(i,j,:,ie),&
             q_five(i,j,:,1,ie),q_five(i,j,:,2,ie),pdel_five(i,j,:),pver_five,hvcoord%ps0,&
             zm_five(i,j,:),zi_five(i,j,:))
    enddo
    enddo

    do i=1,np
    do j=1,np

      call tendency_low_to_high(zm_host(i,j,:),zi_host(i,j,:),zm_five(i,j,:),&
             rho_host(i,j,:),rho_five(i,j,:),t_five_tend_low(i,j,:,ie),t_five_tend(i,j,:,ie))
        
      call tendency_low_to_high(zm_host(i,j,:),zi_host(i,j,:),zm_five(i,j,:),&
             rho_host(i,j,:),rho_five(i,j,:),u_five_tend_low(i,j,:,ie),u_five_tend(i,j,:,ie))
        
      call tendency_low_to_high(zm_host(i,j,:),zi_host(i,j,:),zm_five(i,j,:),&
             rho_host(i,j,:),rho_five(i,j,:),v_five_tend_low(i,j,:,ie),v_five_tend(i,j,:,ie))
        
      do q=1,qsize
         call tendency_low_to_high(zm_host(i,j,:),zi_host(i,j,:),zm_five(i,j,:),&
             rho_host(i,j,:),rho_five(i,j,:),q_five_tend_low(i,j,:,q,ie),q_five_tend(i,j,:,q,ie))
      enddo                             

    enddo       
    enddo       

    ! Finally, update FIVE prognostic variables based on this tendency, so 
    !   complete syncronization with E3SM
    do k=1,pver_five
      do i=1,np
      do j=1,np
        ttmp_five_old(i,j,k) = t_five(i,j,k,ie) + dt * t_five_tend(i,j,k,ie)
        utmp_five_old(i,j,k) = u_five(i,j,k,ie) + dt * u_five_tend(i,j,k,ie)
        vtmp_five_old(i,j,k) = v_five(i,j,k,ie) + dt * v_five_tend(i,j,k,ie)

        do q=1,qsize
          Qtmp_five_old(i,j,k,q) = q_five(i,j,k,q,ie) + dt * q_five_tend(i,j,k,q,ie)
          Qtmp_five_old(i,j,k,q) = max(Qtmp_five_old(i,j,k,q),0._r8)
        enddo

      enddo
      enddo
    enddo

#endif

#ifdef FIVE
     if (rsplit>0) then
goto 9999
     do i=1,np
     do j=1,np
! get t, v, u, dp. In the future, they shoudl pass from pbuf, like T_five

        call linear_interp(pmid_host(i,j,:),pmid_five(i,j,:),elem(ie)%state%t(i,j,:,np1),ttmp_five(i,j,:),pver,pver_five)
        call linear_interp(pmid_host(i,j,:),pmid_five(i,j,:),elem(ie)%state%v(i,j,1,:,np1),utmp_five(i,j,:),pver,pver_five)
        call linear_interp(pmid_host(i,j,:),pmid_five(i,j,:),elem(ie)%state%v(i,j,2,:,np1),vtmp_five(i,j,:),pver,pver_five)
       do q = 1, qsize
         state_Q(i,j,:,q) = elem(ie)%state%Qdp(i,j,:,q,np1_qdp) / dp_star(i,j,:) 
         call linear_interp(pmid_host(i,j,:),pmid_five(i,j,:),state_Q(i,j,:,q),Qtmp_five(i,j,:,q),pver,pver_five)
         Qtmp_five_old(i,j,:,q) = Qtmp_five(i,j,:,q)
         Qtmp_five(i,j,:,q) = Qtmp_five(i,j,:,q)*dp_star_five(i,j,:)

       enddo 
     enddo 
     enddo
9999 continue

     if (single_column) then
     do i=1,np
     do j=1,np
        do k=1,pver_five
          dp_star_five(i,j,k)  = dp_star_five(1,1,k)
          ttmp_five_old(i,j,k) = ttmp_five_old(1,1,k)
          utmp_five_old(i,j,k) = utmp_five_old(1,1,k)
          vtmp_five_old(i,j,k) = vtmp_five_old(1,1,k)
          do q = 1, qsize 
             Qtmp_five_old(i,j,k,q) = Qtmp_five_old(1,1,k,q)
          enddo
        enddo
     enddo 
     enddo
     endif

!        ttmp_five_old = ttmp_five 
!        ttmp_five = ttmp_five*dp_star_five
        ttmp_five = ttmp_five_old*dp_star_five
        call t_startf('vertical_remap1_FIVE_1')
        call remap1(ttmp_five,np,pver_five,1,dp_star_five,dp_five)
        call t_stopf('vertical_remap1_FIVE_1')
        ttmp_five = ttmp_five/dp_five
        ttmp_five_tend = (ttmp_five - ttmp_five_old)/dt

!        utmp_five_old = utmp_five
!        utmp_five = utmp_five*dp_star_five
        utmp_five = utmp_five_old * dp_star_five
        call t_startf('vertical_remap1_FIVE_2')
        call remap1(utmp_five,np,pver_five,1,dp_star_five,dp_five)
        call t_stopf('vertical_remap1_FIVE_2')
        utmp_five = utmp_five/dp_five
        utmp_five_tend = (utmp_five - utmp_five_old)/dt

!        vtmp_five_old = vtmp_five
!        vtmp_five = vtmp_five*dp_star_five
        vtmp_five = vtmp_five_old * dp_star_five
        call t_startf('vertical_remap1_FIVE_3')
        call remap1(vtmp_five,np,pver_five,1,dp_star_five,dp_five)
        call t_stopf('vertical_remap1_FIVE_3')
        vtmp_five = vtmp_five/dp_five
        vtmp_five_tend = (vtmp_five - vtmp_five_old)/dt

     if (qsize>0) then
       
       do q = 1, qsize
         Qtmp_five(:,:,:,q) = Qtmp_five_old(:,:,:,q)*dp_star_five(:,:,:)
       end do
       call t_startf('vertical_remap1_FIVE_4')
       call remap1(Qtmp_five,np,pver_five,qsize,dp_star_five,dp_five)
       call t_stopf('vertical_remap1_FIVE_4')

     endif


! update to remap temperature 
     rho_five = pmid_five/(rair*ttmp_five)
     rho_host = pmid_host/(rair*elem(ie)%state%t(:,:,:,np1))
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
        Qtmp_five(i,j,:,q) = Qtmp_five(i,j,:,q)/dp_five(i,j,:)
        Qtmp_five_tend(i,j,:,q) = (Qtmp_five(i,j,:,q) - Qtmp_five_old(i,j,:,q))/dt

        call masswgt_vert_avg(rho_host(i,j,:),rho_five(i,j,:),dz_host(i,j,:),dz_five(i,j,:),&
                            pint_host(i,j,:),pmid_five(i,j,:),pmid_host(i,j,:),&
                            Qtmp_five_tend(i,j,:,q),Qtmp_host(i,j,:,q))
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
        call remap1(ttmp,np,pver,1,dp_star,dp)
        call t_stopf('vertical_remap1_1')

        elem(ie)%state%t(:,:,:,np1)=ttmp(:,:,:,1)/dp

        ttmp(:,:,:,1)=elem(ie)%state%v(:,:,1,:,np1)*dp_star
        ttmp(:,:,:,2)=elem(ie)%state%v(:,:,2,:,np1)*dp_star

        call t_startf('vertical_remap1_2')
        call remap1(ttmp,np,pver,2,dp_star,dp)
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
         state_Qnew(:,:,:,q) = state_Q(:,:,:,q)+Qtmp_host(:,:,:,q) * dt
         state_Qnew(:,:,:,q) = max(state_Qnew(:,:,:,q),0.)
         elem(ie)%state%Qdp(:,:,:,q,np1_qdp) = state_Qnew(:,:,:,q) * dp(:,:,:)
      enddo
#else
       call t_startf('vertical_remap1_3')
       call remap1(elem(ie)%state%Qdp(:,:,:,:,np1_qdp),np,pver,qsize,dp_star,dp)
       call t_stopf('vertical_remap1_3')
#endif

     endif
!     do i=1,np
!     do j=1,np
!     do k=1,pver
!     do q =1, qsize
!        if (state_Q(i,j,k,q) .ne. 0. .or. state_Qup(i,j,k,q) .ne. 0.) then
!        print *, 'HHLEE LS', elem(ie)%state%t(i,j,1,np1)
!        endif
!     end do
!     end do
!     end do
!     end do

  enddo
  call t_stopf('vertical_remap')
  end subroutine vertical_remap


end module 

