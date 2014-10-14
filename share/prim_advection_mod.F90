#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#if 0
SUBROUTINES:
   prim_advec_tracers_remap_rk2()
      SEM 2D RK2 + monotone remap + hyper viscosity
      SEM 2D RK2 can use sign-preserving or monotone reconstruction

Notes on Lagrange+REMAP advection
dynamics will compute mean fluxes, so that (i.e. for qsplit=3)

    dp(t+3)-dp(t) = -3dt div(Udp_sum/3) - 3dt d(eta_dot_dpdn_sum/3)  + 3dt D(dpdiss_sum/3)

Where the floating lagrangian component:
    dp_star(t+3) = dp(t)  -3dt div(Udp_sum/3)  + 3dt D(dpdiss_sum/3)
OR:
    dp_star(t+3) = dp(t+1) + 3dt d( eta_dot_dpdn_ave(t) )


For RK2 advection of Q:  (example of 2 stage RK for tracers):   dtq = qsplit*dt
For consistency, if Q=1
  dp1  = dp(t)- dtq div[ U1 dp(t)]
  dp2  = dp1  - dtq div[ U2 dp1  ]  + 2*dtq D( dpdiss_ave )
  dp*  = (dp(t) + dp2 )/2
       =  dp(t) - dtq  div[ U1 dp(t) + U2 dp1 ]/2   + dtq D( dpdiss_ave )

so we require:
  U1 = Udp_ave / dp(t)
  U2 = Udp_ave / dp1

For tracer advection:
  Qdp1  = Qdp(t)- dtq div[ U1 Qdp(t)]
  Qdp2  = Qdp1  - dtq div[ U2 Qdp1  ]  + 2*dtq D( Q dpdiss_ave )
  Qdp*  = (Qdp(t) + Qdp2 )/2
       =  Qdp(t) - dtq  div[ U1 Qdp(t) + U2 Qdp1 ]   + dtq D( Q dpdiss_ave )

Qdp1:  limit Q, with Q = Qdp1-before-DSS/(dp1-before-DSS)      with dp1 as computed above
Qdp2:  limit Q, with Q = Qdp2-before-DSS/(dp2-before-DSS)      with dp2 as computed above

For dissipation: Q = Qdp1-after-DSS / dp1-after-DSS


last step:
  remap Qdp* to Qdp(t+1)   [ dp_star(t+1) -> dp(t+1) ]



#endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Begin GPU remap module  !!
!! by Rick Archibald, 2010  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vertremap_mod

  !**************************************************************************************
  !
  !  Purpose:
  !        Construct sub-grid-scale polynomials using piecewise spline method with
  !        monotone filters.
  !
  !  References: PCM - Zerroukat et al., Q.J.R. Meteorol. Soc., 2005. (ZWS2005QJR)
  !              PSM - Zerroukat et al., Int. J. Numer. Meth. Fluids, 2005. (ZWS2005IJMF)
  !
  !**************************************************************************************

  use kinds, only                  : real_kind,int_kind
  use dimensions_mod, only         : np,nlev,qsize,nlevp,npsq,ntrac,nc
  use hybvcoord_mod, only          : hvcoord_t
  use element_mod, only            : element_t
  use fvm_control_volume_mod, only : fvm_struct
  use spelt_mod, only              : spelt_struct
  use perf_mod, only               : t_startf, t_stopf  ! _EXTERNAL
  use parallel_mod, only           : abortmp, parallel_t
  use control_mod, only : vert_remap_q_alg

  public remap1                  ! remap any field, splines, monotone
  public remap1_nofilter         ! remap any field, splines, no filter
! todo: tweak interface to match remap1 above, rename remap1_ppm:
  public remap_q_ppm             ! remap state%Q, PPM, monotone

  contains

!=======================================================================================================!

!remap_calc_grids computes the vertical pressures and pressure differences for one vertical column for the reference grid
!and for the deformed Lagrangian grid. This was pulled out of each routine since it was a repeated task.
subroutine remap_calc_grids( hvcoord , ps , dt , eta_dot_dpdn , p_lag , p_ref , dp_lag , dp_ref )
  implicit none
  type(hvcoord_t)      , intent(in   ) :: hvcoord               !Derived type to hold vertical sigma grid parameters
  real(kind=real_kind) , intent(in   ) :: ps                    !Surface pressure for this column
  real(kind=real_kind) , intent(in   ) :: dt                    !Time step
  real(kind=real_kind) , intent(in   ) :: eta_dot_dpdn(nlev+1)  !Looks like a vertical pressure flux
                                                                !to compute deformed grid spacing
  real(kind=real_kind) , intent(  out) :: p_lag(nlev+1)         !Pressures at interfaces of the Lagrangian deformed grid
  real(kind=real_kind) , intent(  out) :: p_ref(nlev+1)         !Pressures at interfaces of the reference grid
  real(kind=real_kind) , intent(  out) :: dp_lag(nlev)          !Pressure differences on Lagrangian deformed grid
  real(kind=real_kind) , intent(  out) :: dp_ref(nlev)          !Pressure differences on reference grid
  integer :: k                                                  !Iterator
  p_ref(1) = 0  !Both grids have a model top pressure of zero
  p_lag(1) = 0  !Both grids have a model top pressure of zero
  do k = 1 , nlev
    dp_ref(k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) ) * hvcoord%ps0 + &
         ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * ps  !Reference pressure difference
    ! Lagrangian pressure difference (flux in - flux out over the time step)
    dp_lag(k) = dp_ref(k) + dt * ( eta_dot_dpdn(k+1) - eta_dot_dpdn(k) )
    p_ref(k+1) = p_ref(k) + dp_ref(k) !Pressure at interfaces accumulated using difference over each cell
    p_lag(k+1) = p_lag(k) + dp_lag(k) !Pressure at interfaces accumulated using difference over each cell
  enddo
end subroutine remap_calc_grids

!=======================================================================================================!



subroutine remap1(Qdp,nx,qsize,dp1,dp2)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass, monotone on Q=Qdp/dp
  !
  implicit none
  integer, intent(in) :: nx,qsize
  real (kind=real_kind), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=real_kind), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! ========================
  ! Local Variables
  ! ========================

  real (kind=real_kind), dimension(nlev+1)    :: rhs,lower_diag,diag,upper_diag,q_diag,zgam,z1c,z2c,zv
  real (kind=real_kind), dimension(nlev)      :: h,Qcol,dy,za0,za1,za2,zarg,zhdp,dp_star,dp_np1
  real (kind=real_kind)  :: f_xm,level1,level2,level3,level4,level5, &
                            peaks_min,peaks_max,tmp_cal,xm,xm_d,zv1,zv2, &
                            zero = 0,one = 1,tiny = 1e-12,qmax = 1d50
  integer(kind=int_kind) :: zkr(nlev+1),filter_code(nlev),peaks,im1,im2,im3,ip1,ip2, &
                            lt1,lt2,lt3,t0,t1,t2,t3,t4,tm,tp,ie,i,ilev,j,jk,k,q
  logical :: abort=.false.
  call t_startf('remap1')

  if (vert_remap_q_alg == 1 .or. vert_remap_q_alg == 2) then
     call remap_Q_ppm(qdp,nx,qsize,dp1,dp2)
     call t_stopf('remap1')
     return
  endif

#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,i,j,z1c,z2c,zv,k,dp_np1,dp_star,Qcol,zkr,ilev) &
!$omp    private(jk,zgam,zhdp,h,zarg,rhs,lower_diag,diag,upper_diag,q_diag,tmp_cal,filter_code) &
!$omp    private(dy,im1,im2,im3,ip1,t1,t2,t3,za0,za1,za2,xm_d,xm,f_xm,t4,tm,tp,peaks,peaks_min) &
!$omp    private(peaks_max,ip2,level1,level2,level3,level4,level5,lt1,lt2,lt3,zv1,zv2)
#endif
  do q=1,qsize
  do i=1,nx
    do j=1,nx

      z1c(1)=0 ! source grid
      z2c(1)=0 ! target grid
      do k=1,nlev
         z1c(k+1)=z1c(k)+dp1(i,j,k)
         z2c(k+1)=z2c(k)+dp2(i,j,k)
      enddo

      zv(1)=0
      do k=1,nlev
        Qcol(k)=Qdp(i,j,k,q)!  *(z1c(k+1)-z1c(k)) input is mass
        zv(k+1) = zv(k)+Qcol(k)
      enddo

      if (ABS(z2c(nlev+1)-z1c(nlev+1)).GE.0.000001) then
        write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
        write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
        write(6,*) 'DATA FOR MODEL LEVELS'
        write(6,*) 'PLEVMODEL=',z2c(nlev+1)
        write(6,*) 'PLEV     =',z1c(nlev+1)
        write(6,*) 'DIFF     =',z2c(nlev+1)-z1c(nlev+1)
        abort=.true.
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! quadratic splies with UK met office monotonicity constraints  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zkr  = 99
      ilev = 2
      zkr(1) = 1
      zkr(nlev+1) = nlev
      kloop: do k = 2,nlev
        do jk = ilev,nlev+1
          if (z1c(jk).ge.z2c(k)) then
            ilev      = jk
            zkr(k)   = jk-1
            cycle kloop
          endif
        enddo
      enddo kloop

      zgam  = (z2c(1:nlev+1)-z1c(zkr)) / (z1c(zkr+1)-z1c(zkr))
      zgam(1)      = 0.0
      zgam(nlev+1) = 1.0
      zhdp = z1c(2:nlev+1)-z1c(1:nlev)


      h = 1/zhdp
      zarg = Qcol * h
      rhs = 0
      lower_diag = 0
      diag = 0
      upper_diag = 0

      rhs(1)=3*zarg(1)
      rhs(2:nlev) = 3*(zarg(2:nlev)*h(2:nlev) + zarg(1:nlev-1)*h(1:nlev-1))
      rhs(nlev+1)=3*zarg(nlev)

      lower_diag(1)=1
      lower_diag(2:nlev) = h(1:nlev-1)
      lower_diag(nlev+1)=1

      diag(1)=2
      diag(2:nlev) = 2*(h(2:nlev) + h(1:nlev-1))
      diag(nlev+1)=2

      upper_diag(1)=1
      upper_diag(2:nlev) = h(2:nlev)
      upper_diag(nlev+1)=0

      q_diag(1)=-upper_diag(1)/diag(1)
      rhs(1)= rhs(1)/diag(1)

      do k=2,nlev+1
        tmp_cal    =  1/(diag(k)+lower_diag(k)*q_diag(k-1))
        q_diag(k) = -upper_diag(k)*tmp_cal
        rhs(k) =  (rhs(k)-lower_diag(k)*rhs(k-1))*tmp_cal
      enddo
      do k=nlev,1,-1
        rhs(k)=rhs(k)+q_diag(k)*rhs(k+1)
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  monotonicity modifications  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      filter_code = 0
      dy(1:nlev-1) = zarg(2:nlev)-zarg(1:nlev-1)
      dy(nlev) = dy(nlev-1)

      dy = merge(zero, dy, abs(dy) < tiny )

      do k=1,nlev
        im1=MAX(1,k-1)
        im2=MAX(1,k-2)
        im3=MAX(1,k-3)
        ip1=MIN(nlev,k+1)
        t1 = merge(1,0,(zarg(k)-rhs(k))*(rhs(k)-zarg(im1)) >= 0)
        t2 = merge(1,0,dy(im2)*(rhs(k)-zarg(im1)) > 0 .AND. dy(im2)*dy(im3) > 0 &
             .AND. dy(k)*dy(ip1) > 0 .AND. dy(im2)*dy(k) < 0 )
        t3 = merge(1,0,ABS(rhs(k)-zarg(im1)) > ABS(rhs(k)-zarg(k)))

        filter_code(k) = merge(0,1,t1+t2 > 0)
        rhs(k) = (1-filter_code(k))*rhs(k)+filter_code(k)*(t3*zarg(k)+(1-t3)*zarg(im1))
        filter_code(im1) = MAX(filter_code(im1),filter_code(k))
      enddo

      rhs = merge(qmax,rhs,rhs > qmax)
      rhs = merge(zero,rhs,rhs < zero)

      za0 = rhs(1:nlev)
      za1 = -4*rhs(1:nlev) - 2*rhs(2:nlev+1) + 6*zarg
      za2 =  3*rhs(1:nlev) + 3*rhs(2:nlev+1) - 6*zarg

      dy(1:nlev) = rhs(2:nlev+1)-rhs(1:nlev)
      dy = merge(zero, dy, abs(dy) < tiny )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Compute the 3 quadratic spline coeffients {za0, za1, za2}				   !!
      !! knowing the quadratic spline parameters {rho_left,rho_right,zarg}		   !!
      !! Zerroukat et.al., Q.J.R. Meteorol. Soc., Vol. 128, pp. 2801-2820 (2002).   !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      h = rhs(2:nlev+1)

      do k=1,nlev
        xm_d = merge(one,2*za2(k),abs(za2(k)) < tiny)
        xm = merge(zero,-za1(k)/xm_d, abs(za2(k)) < tiny)
        f_xm = za0(k) + za1(k)*xm + za2(k)*xm**2

        t1 = merge(1,0,ABS(za2(k)) > tiny)
        t2 = merge(1,0,xm <= zero .OR. xm >= 1)
        t3 = merge(1,0,za2(k) > zero)
        t4 = merge(1,0,za2(k) < zero)
        tm = merge(1,0,t1*((1-t2)+t3) .EQ. 2)
        tp = merge(1,0,t1*((1-t2)+(1-t3)+t4) .EQ. 3)

        peaks=0
        peaks = merge(-1,peaks,tm .EQ. 1)
        peaks = merge(+1,peaks,tp .EQ. 1)
        peaks_min = merge(f_xm,MIN(za0(k),za0(k)+za1(k)+za2(k)),tm .EQ. 1)
        peaks_max = merge(f_xm,MAX(za0(k),za0(k)+za1(k)+za2(k)),tp .EQ. 1)

        im1=MAX(1,k-1)
        im2=MAX(1,k-2)
        ip1=MIN(nlev,k+1)
        ip2=MIN(nlev,k+2)

        t1 = merge(abs(peaks),0,(dy(im2)*dy(im1) <= tiny) .OR. &
             (dy(ip1)*dy(ip2) <= tiny) .OR. (dy(im1)*dy(ip1) >= tiny) .OR. &
             (dy(im1)*float(peaks) <= tiny))

        filter_code(k) = merge(1,t1+(1-t1)*filter_code(k),(rhs(k) >= qmax) .OR. &
             (rhs(k) <= zero) .OR. (peaks_max > qmax) .OR. (peaks_min < tiny))

        if (filter_code(k) > 0) then
          level1 = rhs(k)
          level2 = (2*rhs(k)+h(k))/3
          level3 = 0.5*(rhs(k)+h(k))
          level4 = (1/3d0)*rhs(k)+2*(1/3d0)*h(k)
          level5 = h(k)

          t1 = merge(1,0,h(k) >= rhs(k))
          t2 = merge(1,0,zarg(k) <= level1 .OR.  zarg(k) >= level5)
          t3 = merge(1,0,zarg(k) >  level1 .AND. zarg(k) <  level2)
          t4 = merge(1,0,zarg(k) >  level4 .AND. zarg(k) <  level5)

          lt1 = t1*t2
          lt2 = t1*(1-t2+t3)
          lt3 = t1*(1-t2+1-t3+t4)

          za0(k) = merge(zarg(k),za0(k),lt1 .EQ. 1)
          za1(k) = merge(zero,za1(k),lt1 .EQ. 1)
          za2(k) = merge(zero,za2(k),lt1 .EQ. 1)

          za0(k) = merge(rhs(k),za0(k),lt2 .EQ. 2)
          za1(k) = merge(zero,za1(k),lt2 .EQ. 2)
          za2(k) = merge(3*(zarg(k)-rhs(k)),za2(k),lt2 .EQ. 2)

          za0(k) = merge(-2*h(k)+3*zarg(k),za0(k),lt3 .EQ. 3)
          za1(k) = merge(+6*h(k)-6*zarg(k),za1(k),lt3 .EQ. 3)
          za2(k) = merge(-3*h(k)+3*zarg(k),za2(k),lt3 .EQ. 3)

          t2 = merge(1,0,zarg(k) >= level1 .OR.  zarg(k) <= level5)
          t3 = merge(1,0,zarg(k) <  level1 .AND. zarg(k) >  level2)
          t4 = merge(1,0,zarg(k) <  level4 .AND. zarg(k) >  level5)

          lt1 = (1-t1)*t2
          lt2 = (1-t1)*(1-t2+t3)
          lt3 = (1-t1)*(1-t2+1-t3+t4)

          za0(k) = merge(zarg(k),za0(k),lt1 .EQ. 1)
          za1(k) = merge(zero,za1(k),lt1 .EQ. 1)
          za2(k) = merge(zero,za2(k),lt1 .EQ. 1)

          za0(k) = merge(rhs(k),za0(k),lt2 .EQ. 2)
          za1(k) = merge(zero,za1(k),lt2 .EQ. 2)
          za2(k) = merge(3*(zarg(k)-rhs(k)),za2(k),lt2 .EQ. 2)

          za0(k) = merge(-2*h(k)+3*zarg(k),za0(k),lt3 .EQ. 3)
          za1(k) = merge(+6*h(k)-6*zarg(k),za1(k),lt3 .EQ. 3)
          za2(k) = merge(-3*h(k)+3*zarg(k),za2(k),lt3 .EQ. 3)
        endif
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! start iteration from top to bottom of atmosphere !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zv1 = 0
      do k=1,nlev
        if (zgam(k+1)>1d0) then
          WRITE(*,*) 'r not in [0:1]', zgam(k+1)
          abort=.true.
        endif
        zv2 = zv(zkr(k+1))+(za0(zkr(k+1))*zgam(k+1)+(za1(zkr(k+1))/2)*(zgam(k+1)**2)+ &
             (za2(zkr(k+1))/3)*(zgam(k+1)**3))*zhdp(zkr(k+1))
        Qdp(i,j,k,q) = (zv2 - zv1) ! / (z2c(k+1)-z2c(k) ) dont convert back to mixing ratio
        zv1 = zv2
      enddo
    enddo
  enddo
  enddo ! q loop
  if (abort) call abortmp('Bad levels in remap1.  usually CFL violatioin')
  call t_stopf('remap1')
end subroutine remap1

subroutine remap1_nofilter(Qdp,nx,qsize,dp1,dp2)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  implicit none
  integer, intent(in) :: nx,qsize
  real (kind=real_kind), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=real_kind), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! ========================
  ! Local Variables
  ! ========================

  real (kind=real_kind), dimension(nlev+1)    :: rhs,lower_diag,diag,upper_diag,q_diag,zgam,z1c,z2c,zv
  real (kind=real_kind), dimension(nlev)      :: h,Qcol,dy,za0,za1,za2,zarg,zhdp,dp_star,dp_np1
  real (kind=real_kind)  :: f_xm,level1,level2,level3,level4,level5, &
                            peaks_min,peaks_max,tmp_cal,xm,xm_d,zv1,zv2, &
                            zero = 0,one = 1,tiny = 1e-12,qmax = 1d50
  integer(kind=int_kind) :: zkr(nlev+1),filter_code(nlev),peaks,im1,im2,im3,ip1,ip2, &
                            lt1,lt2,lt3,t0,t1,t2,t3,t4,tm,tp,ie,i,ilev,j,jk,k,q
  logical :: abort=.false.
!   call t_startf('remap1_nofilter')

#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,i,j,z1c,z2c,zv,k,dp_np1,dp_star,Qcol,zkr,ilev) &
!$omp    private(jk,zgam,zhdp,h,zarg,rhs,lower_diag,diag,upper_diag,q_diag,tmp_cal,filter_code) &
!$omp    private(dy,im1,im2,im3,ip1,t1,t2,t3,za0,za1,za2,xm_d,xm,f_xm,t4,tm,tp,peaks,peaks_min) &
!$omp    private(peaks_max,ip2,level1,level2,level3,level4,level5,lt1,lt2,lt3,zv1,zv2)
#endif
  do q=1,qsize
  do i=1,nx
    do j=1,nx

      z1c(1)=0 ! source grid
      z2c(1)=0 ! target grid
      do k=1,nlev
         z1c(k+1)=z1c(k)+dp1(i,j,k)
         z2c(k+1)=z2c(k)+dp2(i,j,k)
      enddo

      zv(1)=0
      do k=1,nlev
        Qcol(k)=Qdp(i,j,k,q)!  *(z1c(k+1)-z1c(k)) input is mass
        zv(k+1) = zv(k)+Qcol(k)
      enddo

      if (ABS(z2c(nlev+1)-z1c(nlev+1)).GE.0.000001) then
        write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
        write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
        write(6,*) 'DATA FOR MODEL LEVELS'
        write(6,*) 'PLEVMODEL=',z2c(nlev+1)
        write(6,*) 'PLEV     =',z1c(nlev+1)
        write(6,*) 'DIFF     =',z2c(nlev+1)-z1c(nlev+1)
        abort=.true.
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! quadratic splies with UK met office monotonicity constraints  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zkr  = 99
      ilev = 2
      zkr(1) = 1
      zkr(nlev+1) = nlev
      kloop: do k = 2,nlev
        do jk = ilev,nlev+1
          if (z1c(jk).ge.z2c(k)) then
            ilev      = jk
            zkr(k)   = jk-1
            cycle kloop
          endif
        enddo
      enddo kloop

      zgam  = (z2c(1:nlev+1)-z1c(zkr)) / (z1c(zkr+1)-z1c(zkr))
      zgam(1)      = 0.0
      zgam(nlev+1) = 1.0
      zhdp = z1c(2:nlev+1)-z1c(1:nlev)


      h = 1/zhdp
      zarg = Qcol * h
      rhs = 0
      lower_diag = 0
      diag = 0
      upper_diag = 0

      rhs(1)=3*zarg(1)
      rhs(2:nlev) = 3*(zarg(2:nlev)*h(2:nlev) + zarg(1:nlev-1)*h(1:nlev-1))
      rhs(nlev+1)=3*zarg(nlev)

      lower_diag(1)=1
      lower_diag(2:nlev) = h(1:nlev-1)
      lower_diag(nlev+1)=1

      diag(1)=2
      diag(2:nlev) = 2*(h(2:nlev) + h(1:nlev-1))
      diag(nlev+1)=2

      upper_diag(1)=1
      upper_diag(2:nlev) = h(2:nlev)
      upper_diag(nlev+1)=0

      q_diag(1)=-upper_diag(1)/diag(1)
      rhs(1)= rhs(1)/diag(1)

      do k=2,nlev+1
        tmp_cal    =  1/(diag(k)+lower_diag(k)*q_diag(k-1))
        q_diag(k) = -upper_diag(k)*tmp_cal
        rhs(k) =  (rhs(k)-lower_diag(k)*rhs(k-1))*tmp_cal
      enddo
      do k=nlev,1,-1
        rhs(k)=rhs(k)+q_diag(k)*rhs(k+1)
      enddo

      za0 = rhs(1:nlev)
      za1 = -4*rhs(1:nlev) - 2*rhs(2:nlev+1) + 6*zarg
      za2 =  3*rhs(1:nlev) + 3*rhs(2:nlev+1) - 6*zarg


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! start iteration from top to bottom of atmosphere !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zv1 = 0
      do k=1,nlev
        if (zgam(k+1)>1d0) then
          WRITE(*,*) 'r not in [0:1]', zgam(k+1)
          abort=.true.
        endif
        zv2 = zv(zkr(k+1))+(za0(zkr(k+1))*zgam(k+1)+(za1(zkr(k+1))/2)*(zgam(k+1)**2)+ &
             (za2(zkr(k+1))/3)*(zgam(k+1)**3))*zhdp(zkr(k+1))
        Qdp(i,j,k,q) = (zv2 - zv1) ! / (z2c(k+1)-z2c(k) ) dont convert back to mixing ratio
        zv1 = zv2
      enddo
    enddo
  enddo
  enddo ! q loop
  if (abort) call abortmp('Bad levels in remap1_nofilter.  usually CFL violatioin')
!   call t_stopf('remap1_nofilter')
end subroutine remap1_nofilter

!=======================================================================================================!


!This uses the exact same model and reference grids and data as remap_Q, but it interpolates
!using PPM instead of splines.
subroutine remap_Q_ppm(Qdp,nx,qsize,dp1,dp2)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  use control_mod, only        : prescribed_wind, vert_remap_q_alg
  implicit none
  integer,intent(in) :: nx,qsize
  real (kind=real_kind), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=real_kind), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! Local Variables
  integer, parameter :: gs = 2                              !Number of cells to place in the ghost region
  real(kind=real_kind), dimension(       nlev+2 ) :: pio    !Pressure at interfaces for old grid
  real(kind=real_kind), dimension(       nlev+1 ) :: pin    !Pressure at interfaces for new grid
  real(kind=real_kind), dimension(       nlev+1 ) :: masso  !Accumulate mass up to each interface
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: ao     !Tracer value on old grid
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpo    !change in pressure over a cell for old grid
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpn    !change in pressure over a cell for old grid
  real(kind=real_kind), dimension(3,     nlev   ) :: coefs  !PPM coefficients within each cell
  real(kind=real_kind), dimension(       nlev   ) :: z1, z2
  real(kind=real_kind) :: ppmdx(10,0:nlev+1)  !grid spacings
  real(kind=real_kind) :: mymass, massn1, massn2
  integer :: i, j, k, q, kk, kid(nlev)

  call t_startf('remap_Q_ppm')
  do j = 1 , nx
    do i = 1 , nx

      pin(1)=0
      pio(1)=0
      do k=1,nlev
         dpn(k)=dp2(i,j,k)
         dpo(k)=dp1(i,j,k)
         pin(k+1)=pin(k)+dpn(k)
         pio(k+1)=pio(k)+dpo(k)
      enddo



      pio(nlev+2) = pio(nlev+1) + 1.  !This is here to allow an entire block of k threads to run in the remapping phase.
                                      !It makes sure there's an old interface value below the domain that is larger.
      pin(nlev+1) = pio(nlev+1)       !The total mass in a column does not change.
                                      !Therefore, the pressure of that mass cannot either.
      !Fill in the ghost regions with mirrored values. if vert_remap_q_alg is defined, this is of no consequence.
      do k = 1 , gs
        dpo(1   -k) = dpo(       k)
        dpo(nlev+k) = dpo(nlev+1-k)
      enddo

      !Compute remapping intervals once for all tracers. Find the old grid cell index in which the
      !k-th new cell interface resides. Then integrate from the bottom of that old cell to the new
      !interface location. In practice, the grid never deforms past one cell, so the search can be
      !simplified by this. Also, the interval of integration is usually of magnitude close to zero
      !or close to dpo because of minimial deformation.
      !Numerous tests confirmed that the bottom and top of the grids match to machine precision, so
      !I set them equal to each other.
      do k = 1 , nlev
        kk = k  !Keep from an order n^2 search operation by assuming the old cell index is close.
        !Find the index of the old grid cell in which this new cell's bottom interface resides.
        do while ( pio(kk) <= pin(k+1) )
          kk = kk + 1
        enddo
        kk = kk - 1                   !kk is now the cell index we're integrating over.
        if (kk == nlev+1) kk = nlev   !This is to keep the indices in bounds.
                                      !Top bounds match anyway, so doesn't matter what coefficients are used
        kid(k) = kk                   !Save for reuse
        z1(k) = -0.5D0                !This remapping assumes we're starting from the left interface of an old grid cell
                                      !In fact, we're usually integrating very little or almost all of the cell in question
        z2(k) = ( pin(k+1) - ( pio(kk) + pio(kk+1) ) * 0.5 ) / dpo(kk)  !PPM interpolants are normalized to an independent
                                                                        !coordinate domain [-0.5,0.5].
      enddo

      !This turned out a big optimization, remembering that only parts of the PPM algorithm depends on the data, namely the
      !limiting. So anything that depends only on the grid is pre-computed outside the tracer loop.
      ppmdx(:,:) = compute_ppm_grids( dpo )

      !From here, we loop over tracers for only those portions which depend on tracer data, which includes PPM limiting and
      !mass accumulation
      do q = 1 , qsize
        !Accumulate the old mass up to old grid cell interface locations to simplify integration
        !during remapping. Also, divide out the grid spacing so we're working with actual tracer
        !values and can conserve mass. The option for ifndef ZEROHORZ I believe is there to ensure
        !tracer consistency for an initially uniform field. I copied it from the old remap routine.
        masso(1) = 0.
        do k = 1 , nlev
          ao(k) = Qdp(i,j,k,q)
          masso(k+1) = masso(k) + ao(k) !Accumulate the old mass. This will simplify the remapping
          ao(k) = ao(k) / dpo(k)        !Divide out the old grid spacing because we want the tracer mixing ratio, not mass.
        enddo
        !Fill in ghost values. Ignored if vert_remap_q_alg == 2
        do k = 1 , gs
          ao(1   -k) = ao(       k)
          ao(nlev+k) = ao(nlev+1-k)
        enddo
        !Compute monotonic and conservative PPM reconstruction over every cell
        coefs(:,:) = compute_ppm( ao , ppmdx )
        !Compute tracer values on the new grid by integrating from the old cell bottom to the new
        !cell interface to form a new grid mass accumulation. Taking the difference between
        !accumulation at successive interfaces gives the mass inside each cell. Since Qdp is
        !supposed to hold the full mass this needs no normalization.
        massn1 = 0.
        do k = 1 , nlev
          kk = kid(k)
          massn2 = masso(kk) + integrate_parabola( coefs(:,kk) , z1(k) , z2(k) ) * dpo(kk)
          Qdp(i,j,k,q) = massn2 - massn1
          massn1 = massn2
        enddo
      enddo
    enddo
  enddo
  call t_stopf('remap_Q_ppm')
end subroutine remap_Q_ppm


!=======================================================================================================!


!THis compute grid-based coefficients from Collela & Woodward 1984.
function compute_ppm_grids( dx )   result(rslt)
  use control_mod, only: vert_remap_q_alg
  implicit none
  real(kind=real_kind), intent(in) :: dx(-1:nlev+2)  !grid spacings
  real(kind=real_kind)             :: rslt(10,0:nlev+1)  !grid spacings
  integer :: j
  integer :: indB, indE

  !Calculate grid-based coefficients for stage 1 of compute_ppm
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-1
  else
    indB = 0
    indE = nlev+1
  endif
  do j = indB , indE
    rslt( 1,j) = dx(j) / ( dx(j-1) + dx(j) + dx(j+1) )
    rslt( 2,j) = ( 2.*dx(j-1) + dx(j) ) / ( dx(j+1) + dx(j) )
    rslt( 3,j) = ( dx(j) + 2.*dx(j+1) ) / ( dx(j-1) + dx(j) )
  enddo

  !Caculate grid-based coefficients for stage 2 of compute_ppm
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  do j = indB , indE
    rslt( 4,j) = dx(j) / ( dx(j) + dx(j+1) )
    rslt( 5,j) = 1. / sum( dx(j-1:j+2) )
    rslt( 6,j) = ( 2. * dx(j+1) * dx(j) ) / ( dx(j) + dx(j+1 ) )
    rslt( 7,j) = ( dx(j-1) + dx(j  ) ) / ( 2. * dx(j  ) + dx(j+1) )
    rslt( 8,j) = ( dx(j+2) + dx(j+1) ) / ( 2. * dx(j+1) + dx(j  ) )
    rslt( 9,j) = dx(j  ) * ( dx(j-1) + dx(j  ) ) / ( 2.*dx(j  ) +    dx(j+1) )
    rslt(10,j) = dx(j+1) * ( dx(j+1) + dx(j+2) ) / (    dx(j  ) + 2.*dx(j+1) )
  enddo
end function compute_ppm_grids

!=======================================================================================================!



!This computes a limited parabolic interpolant using a net 5-cell stencil, but the stages of computation are broken up into 3 stages
function compute_ppm( a , dx )    result(coefs)
  use control_mod, only: vert_remap_q_alg
  implicit none
  real(kind=real_kind), intent(in) :: a    (    -1:nlev+2)  !Cell-mean values
  real(kind=real_kind), intent(in) :: dx   (10,  0:nlev+1)  !grid spacings
  real(kind=real_kind) ::             coefs(0:2,   nlev  )  !PPM coefficients (for parabola)
  real(kind=real_kind) :: ai (0:nlev  )                     !fourth-order accurate, then limited interface values
  real(kind=real_kind) :: dma(0:nlev+1)                     !An expression from Collela's '84 publication
  real(kind=real_kind) :: da                                !Ditto
  ! Hold expressions based on the grid (which are cumbersome).
  real(kind=real_kind) :: dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8, dx9, dx10
  real(kind=real_kind) :: al, ar                            !Left and right interface values for cell-local limiting
  integer :: j
  integer :: indB, indE

  ! Stage 1: Compute dma for each cell, allowing a 1-cell ghost stencil below and above the domain
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-1
  else
    indB = 0
    indE = nlev+1
  endif
  do j = indB , indE
    da = dx(1,j) * ( dx(2,j) * ( a(j+1) - a(j) ) + dx(3,j) * ( a(j) - a(j-1) ) )
    dma(j) = minval( (/ abs(da) , 2. * abs( a(j) - a(j-1) ) , 2. * abs( a(j+1) - a(j) ) /) ) * sign(1.D0,da)
    if ( ( a(j+1) - a(j) ) * ( a(j) - a(j-1) ) <= 0. ) dma(j) = 0.
  enddo

  ! Stage 2: Compute ai for each cell interface in the physical domain (dimension nlev+1)
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  do j = indB , indE
    ai(j) = a(j) + dx(4,j) * ( a(j+1) - a(j) ) + dx(5,j) * ( dx(6,j) * ( dx(7,j) - dx(8,j) ) &
         * ( a(j+1) - a(j) ) - dx(9,j) * dma(j+1) + dx(10,j) * dma(j) )
  enddo

  ! Stage 3: Compute limited PPM interpolant over each cell in the physical domain
  ! (dimension nlev) using ai on either side and ao within the cell.
  if (vert_remap_q_alg == 2) then
    indB = 3
    indE = nlev-2
  else
    indB = 1
    indE = nlev
  endif
  do j = indB , indE
    al = ai(j-1)
    ar = ai(j  )
    if ( (ar - a(j)) * (a(j) - al) <= 0. ) then
      al = a(j)
      ar = a(j)
    endif
    if ( (ar - al) * (a(j) - (al + ar)/2.) >  (ar - al)**2/6. ) al = 3.*a(j) - 2. * ar
    if ( (ar - al) * (a(j) - (al + ar)/2.) < -(ar - al)**2/6. ) ar = 3.*a(j) - 2. * al
    !Computed these coefficients from the edge values and cell mean in Maple. Assumes normalized coordinates: xi=(x-x0)/dx
    coefs(0,j) = 1.5 * a(j) - ( al + ar ) / 4.
    coefs(1,j) = ar - al
    coefs(2,j) = -6. * a(j) + 3. * ( al + ar )
  enddo

  !If we're not using a mirrored boundary condition, then make the two cells bordering the top and bottom
  !material boundaries piecewise constant. Zeroing out the first and second moments, and setting the zeroth
  !moment to the cell mean is sufficient to maintain conservation.
  if (vert_remap_q_alg == 2) then
    coefs(0,1:2) = a(1:2)
    coefs(1:2,1:2) = 0.
    coefs(0,nlev-1:nlev) = a(nlev-1:nlev)
    coefs(1:2,nlev-1:nlev) = 0.D0
  endif
end function compute_ppm

!=======================================================================================================!


!Simple function computes the definite integral of a parabola in normalized coordinates, xi=(x-x0)/dx,
!given two bounds. Make sure this gets inlined during compilation.
function integrate_parabola( a , x1 , x2 )    result(mass)
  implicit none
  real(kind=real_kind), intent(in) :: a(0:2)  !Coefficients of the parabola
  real(kind=real_kind), intent(in) :: x1      !lower domain bound for integration
  real(kind=real_kind), intent(in) :: x2      !upper domain bound for integration
  real(kind=real_kind)             :: mass
  mass = a(0) * (x2 - x1) + a(1) * (x2 ** 2 - x1 ** 2) / 0.2D1 + a(2) * (x2 ** 3 - x1 ** 3) / 0.3D1
end function integrate_parabola


!=============================================================================================!



end module vertremap_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End GPU remap module    !!
!! by Rick Archibald, 2010  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!=======================================================================================================!




module prim_advection_mod
!
! two formulations.  both are conservative
! u grad Q formulation:
!
!    d/dt[ Q] +  U grad Q  +  eta_dot dp/dn dQ/dp  = 0
!                            ( eta_dot dQ/dn )
!
!    d/dt[ dp/dn ] = div( dp/dn U ) + d/dn ( eta_dot dp/dn )
!
! total divergence formulation:
!    d/dt[dp/dn Q] +  div( U dp/dn Q ) + d/dn ( eta_dot dp/dn Q ) = 0
!
! for convience, rewrite this as dp Q:  (since dn does not depend on time or the horizonal):
! equation is now:
!    d/dt[dp Q] +  div( U dp Q ) + d( eta_dot_dpdn Q ) = 0
!
!
  use kinds, only              : real_kind
  use dimensions_mod, only     : nlev, nlevp, np, qsize, ntrac, nc, nep
  use physical_constants, only : rgas, Rwater_vapor, kappa, g, rearth, rrearth, cp
  use derivative_mod, only     : gradient, vorticity, gradient_wk, derivative_t, divergence, &
                                 gradient_sphere, divergence_sphere
  use element_mod, only        : element_t
  use fvm_control_volume_mod, only        : fvm_struct
  use spelt_mod, only          : spelt_struct
  use filter_mod, only         : filter_t, filter_P
  use hybvcoord_mod, only      : hvcoord_t
  use time_mod, only           : TimeLevel_t, smooth, TimeLevel_Qdp
  use prim_si_mod, only        : preq_pressure
  use diffusion_mod, only      : scalar_diffusion, diffusion_init
  use control_mod, only        : integration, test_case, filter_freq_advection,  hypervis_order, &
        statefreq, moisture, TRACERADV_TOTAL_DIVERGENCE, TRACERADV_UGRADQ, &
        prescribed_wind, nu_q, nu_p, limiter_option, hypervis_subcycle_q, rsplit
  use edge_mod, only           : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack, initedgebuffer, edgevunpackmin, ghostbuffer3D_t
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use viscosity_mod, only      : biharmonic_wk_scalar, biharmonic_wk_scalar_minmax, neighbor_minmax
  use perf_mod, only           : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only   : abortmp

  implicit none

  private
  save

  public :: Prim_Advec_Init1, Prim_Advec_Init2
  public :: Prim_Advec_Tracers_remap, Prim_Advec_Tracers_remap_rk2, Prim_Advec_Tracers_remap_ALE
  public :: prim_advec_tracers_fvm
#if defined(_SPELT)
  public :: prim_advec_tracers_spelt
#endif
  public :: vertical_remap

  type (EdgeBuffer_t) :: edgeAdv, edgeAdvQ3, edgeAdv_p1, edgeAdvQ2, edgeAdv1,  edgeveloc
  type (ghostBuffer3D_t)   :: ghostbuf_tr

  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  real(kind=real_kind), allocatable :: qmin(:,:,:), qmax(:,:,:)

  type (derivative_t), public, allocatable   :: deriv(:) ! derivative struct (nthreads)

contains

  subroutine Prim_Advec_Init1(par, n_domains)
    use dimensions_mod, only : nlev, qsize, nelemd
    use parallel_mod, only : parallel_t
    type(parallel_t) :: par
    integer, intent(in) :: n_domains

    ! Shared buffer pointers.
    ! Using "=> null()" in a subroutine is usually bad, because it makes
    ! the variable have an implicit "save", and therefore shared between
    ! threads. But in this case we want shared pointers.
    real(kind=real_kind), pointer :: buf_ptr(:) => null()
    real(kind=real_kind), pointer :: receive_ptr(:) => null()

    ! this might be called with qsize=0
    ! allocate largest one first
    ! Currently this is never freed. If it was, only this first one should
    ! be freed, as only it knows the true size of the buffer.
    call initEdgeBuffer(par,edgeAdvQ3,max(nlev,qsize*nlev*3), buf_ptr, receive_ptr)  ! Qtens,Qmin, Qmax

    ! remaining edge buffers can share %buf and %receive with edgeAdvQ3
    ! (This is done through the optional 1D pointer arguments.)
    call initEdgeBuffer(par,edgeAdv1,nlev,buf_ptr,receive_ptr)
    call initEdgeBuffer(par,edgeAdv,qsize*nlev,buf_ptr,receive_ptr)
    call initEdgeBuffer(par,edgeAdv_p1,qsize*nlev + nlev,buf_ptr,receive_ptr)
    call initEdgeBuffer(par,edgeAdvQ2,qsize*nlev*2,buf_ptr,receive_ptr)  ! Qtens,Qmin, Qmax
    call initEdgeBuffer(par,edgeveloc,2*nlev)

    ! Don't actually want these saved, if this is ever called twice.
    nullify(buf_ptr)
    nullify(receive_ptr)

    allocate(deriv(0:n_domains-1))

    ! this static array is shared by all threads, so dimension for all threads (nelemd), not nets:nete:
    allocate (qmin(nlev,qsize,nelemd))
    allocate (qmax(nlev,qsize,nelemd))

  end subroutine Prim_Advec_Init1

  subroutine Prim_Advec_Init2(hybrid,fvm_corners, fvm_points, spelt_refnep)
    use kinds,          only : longdouble_kind
    use dimensions_mod, only : nc, nep
    use derivative_mod, only : derivinit

    type (hybrid_t), intent(in) :: hybrid
    real(kind=longdouble_kind), intent(in) :: fvm_corners(nc+1)
    real(kind=longdouble_kind), intent(in) :: fvm_points(nc)
    real(kind=longdouble_kind), intent(in) :: spelt_refnep(1:nep)

    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call derivinit(deriv(hybrid%ithr),fvm_corners, fvm_points, spelt_refnep)
  end subroutine Prim_Advec_Init2

#if defined(_SPELT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SPELT driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Prim_Advec_Tracers_spelt(elem, spelt, deriv,hvcoord,hybrid,&
        dt,tl,nets,nete)
    use perf_mod, only : t_startf, t_stopf            ! _EXTERNAL
    use spelt_mod, only : spelt_run, spelt_runpos, spelt_runair, spelt_runlimit, edgeveloc, spelt_mcgregordss, spelt_rkdss
    use derivative_mod, only : interpolate_gll2spelt_points
    use vertremap_mod, only: remap1_nofilter ! _EXTERNAL (actually INTERNAL)

    implicit none
    type (element_t), intent(inout)   :: elem(:)
    type (spelt_struct), intent(inout)  :: spelt(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t),     intent(in):: hybrid
    type (TimeLevel_t)                :: tl

    real(kind=real_kind) , intent(in) :: dt
    integer,intent(in)                :: nets,nete

    real (kind=real_kind), dimension(np,np,nlev)    :: dp_star
    real (kind=real_kind), dimension(np,np,nlev)    :: dp

    integer :: np1,ie,k, i, j


    call t_barrierf('sync_prim_advec_tracers_spelt', hybrid%par%comm)
    call t_startf('prim_advec_tracers_spelt')
    np1 = tl%np1
    ! departure algorithm requires two velocities:
    !
    ! spelt%v0:      velocity at beginning of tracer timestep (time n0_qdp)
    !                this was saved before the (possibly many) dynamics steps
    ! elem%derived%vstar:    
    !                velocity at end of tracer timestep (time np1 = np1_qdp)
    !                for lagrangian dynamics, this is on lagrangian levels
    !                for eulerian dynamcis, this is on reference levels
    !                and it should be interpolated.
    !
    do ie=nets,nete
       elem(ie)%derived%vstar(:,:,:,:)=elem(ie)%state%v(:,:,:,:,np1)
    enddo

    if (rsplit==0) then
       ! interpolate t+1 velocity from reference levels to lagrangian levels
       ! For rsplit=0, we need to first compute lagrangian levels based on vertical velocity
       ! which requires we first DSS mean vertical velocity from dynamics
       do ie=nets,nete
          do k=1,nlev
             elem(ie)%derived%eta_dot_dpdn(:,:,k) = elem(ie)%spheremp(:,:)*elem(ie)%derived%eta_dot_dpdn(:,:,k)
          enddo
          call edgeVpack(edgeAdv1,elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev),nlev,0,elem(ie)%desc)
       enddo
       call bndry_exchangeV(hybrid,edgeAdv1)
       do ie=nets,nete
          call edgeVunpack(edgeAdv1,elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev),nlev,0,elem(ie)%desc)
          do k=1,nlev
             elem(ie)%derived%eta_dot_dpdn(:,:,k)=elem(ie)%derived%eta_dot_dpdn(:,:,k)*elem(ie)%rspheremp(:,:)
          enddo
          ! SET VERTICAL VELOCITY TO ZERO FOR DEBUGGING
          elem(ie)%derived%eta_dot_dpdn(:,:,:)=0

          ! elem%state%u(np1)  = velocity at time t+1 on reference levels
          ! elem%derived%vstar = velocity at t+1 on floating levels (computed below) using eta_dot_dpdn
!           call remap_UV_ref2lagrange(np1,dt,elem,hvcoord,ie)
          do k=1,nlev
             dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                  ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
             dp_star(:,:,k) = dp(:,:,k) + dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) -&
                  elem(ie)%derived%eta_dot_dpdn(:,:,k))
          enddo
          !elem(ie)%derived%vstar=elem(ie)%state%v(:,:,:,:,np1) ; done above
          call remap1_nofilter(elem(ie)%derived%vstar,np,1,dp,dp_star)
          !take the average on level, should be improved later, because we know the SE velocity at t+1/2
          spelt(ie)%vn12=(spelt(ie)%vn0+elem(ie)%derived%vstar)/2.0D0
       end do
    else
       ! for rsplit>0:  do nothing.
       ! dynamics is also vertically lagrangian, so we do not need to
       ! remap the velocities and we dont need eta_dot_dpdn
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2D advection step
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------

!     call t_startf('spelt_depalg')
!     call spelt_mcgregordss(elem,spelt,nets,nete, hybrid, deriv, dt, 3)
    call spelt_rkdss(elem,spelt,nets,nete, hybrid, deriv, dt, 3)
!     call t_stopf('spelt_depalg')

    ! ! end mcgregordss
    ! spelt departure calcluation should use vstar.
    ! from c(n0) compute c(np1):
!     call spelt_run(elem,spelt,hybrid,deriv,dt,tl,nets,nete)

    call spelt_runair(elem,spelt,hybrid,deriv,dt,tl,nets,nete)
!     call spelt_runpos(elem,spelt,hybrid,deriv,dt,tl,nets,nete)
!       call spelt_runlimit(elem,spelt,hybrid,deriv,dt,tl,nets,nete)
    call t_stopf('prim_advec_tracers_spelt')
  end subroutine Prim_Advec_Tracers_spelt
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fvm driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Prim_Advec_Tracers_fvm(elem, fvm, deriv,hvcoord,hybrid,&
        dt,tl,nets,nete)
    use perf_mod, only : t_startf, t_stopf            ! _EXTERNAL
    use vertremap_mod, only: remap1_nofilter  ! _EXTERNAL (actually INTERNAL)
    use fvm_mod, only : cslam_runairdensity, cslam_runflux, edgeveloc
    use fvm_mod, only: fvm_mcgregor, fvm_mcgregordss, fvm_rkdss
    use fvm_mod, only : fvm_ideal_test, IDEAL_TEST_OFF, IDEAL_TEST_ANALYTICAL_WINDS
    use fvm_mod, only : fvm_test_type, IDEAL_TEST_BOOMERANG, IDEAL_TEST_SOLIDBODY
    use fvm_bsp_mod, only: get_boomerang_velocities_gll, get_solidbody_velocities_gll
    use control_mod, only : tracer_transport_type
    use control_mod, only : TRACERTRANSPORT_LAGRANGIAN_FVM, TRACERTRANSPORT_FLUXFORM_FVM
    use control_mod, only : use_semi_lagrange_transport
    use time_mod,    only : time_at

    implicit none
    type (element_t), intent(inout)   :: elem(:)
    type (fvm_struct), intent(inout)   :: fvm(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t),     intent(in):: hybrid
    type (TimeLevel_t)                :: tl

    real(kind=real_kind) , intent(in) :: dt
    integer,intent(in)                :: nets,nete


    real (kind=real_kind), dimension(np,np,nlev)    :: dp_star
    real (kind=real_kind), dimension(np,np,nlev)    :: dp
    real (kind=real_kind)                           :: eta_dot_dpdn(np,np,nlevp) 
    integer :: np1,ie,k

    real (kind=real_kind)  :: vstar(np,np,2)
    real (kind=real_kind)  :: vhat(np,np,2)
    real (kind=real_kind), dimension(np, np) :: v1, v2


    call t_barrierf('sync_prim_advec_tracers_fvm', hybrid%par%comm)
    call t_startf('prim_advec_tracers_fvm')
    np1 = tl%np1

    ! departure algorithm requires two velocities:
    !
    ! fvm%v0:        velocity at beginning of tracer timestep (time n0_qdp)
    !                this was saved before the (possibly many) dynamics steps
    ! elem%derived%vstar:    
    !                velocity at end of tracer timestep (time np1 = np1_qdp)
    !                for lagrangian dynamics, this is on lagrangian levels
    !                for eulerian dynamcis, this is on reference levels
    !                and it should be interpolated.
    !
    do ie=nets,nete
       elem(ie)%derived%vstar(:,:,:,:)=elem(ie)%state%v(:,:,:,:,np1)
    enddo


    if (rsplit==0) then
       ! interpolate t+1 velocity from reference levels to lagrangian levels
       ! For rsplit=0, we need to first compute lagrangian levels based on vertical velocity
       ! which requires we first DSS mean vertical velocity from dynamics
       ! note: we introduce a local eta_dot_dpdn() variable instead of DSSing elem%eta_dot_dpdn
       ! so as to preserve BFB results in some HOMME regression tests
       do ie=nets,nete
          do k=1,nlevp
             eta_dot_dpdn(:,:,k) = elem(ie)%derived%eta_dot_dpdn(:,:,k)*elem(ie)%spheremp(:,:)
          enddo
          ! eta_dot_dpdn at nlevp is zero, so we dont boundary exchange it:
          call edgeVpack(edgeAdv1,eta_dot_dpdn(:,:,1:nlev),nlev,0,elem(ie)%desc)
       enddo
       call bndry_exchangeV(hybrid,edgeAdv1)
       do ie=nets,nete
          ! restor interior values.  we could avoid this if we created a global array for eta_dot_dpdn
          do k=1,nlevp
             eta_dot_dpdn(:,:,k) = elem(ie)%derived%eta_dot_dpdn(:,:,k)*elem(ie)%spheremp(:,:)
          enddo
          ! unpack DSSed edge data
          call edgeVunpack(edgeAdv1,eta_dot_dpdn(:,:,1:nlev),nlev,0,elem(ie)%desc)
          do k=1,nlevp
             eta_dot_dpdn(:,:,k) = eta_dot_dpdn(:,:,k)*elem(ie)%rspheremp(:,:)
          enddo


          ! SET VERTICAL VELOCITY TO ZERO FOR DEBUGGING
          !        elem(ie)%derived%eta_dot_dpdn(:,:,:)=0
          ! elem%state%u(np1)  = velocity at time t+1 on reference levels
          ! elem%derived%vstar = velocity at t+1 on floating levels (computed below)
!           call remap_UV_ref2lagrange(np1,dt,elem,hvcoord,ie)
          do k=1,nlev
             dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                  ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
             dp_star(:,:,k) = dp(:,:,k) + dt*(eta_dot_dpdn(:,:,k+1) -   eta_dot_dpdn(:,:,k))
             if (fvm_ideal_test == IDEAL_TEST_ANALYTICAL_WINDS) then
               if (fvm_test_type == IDEAL_TEST_BOOMERANG) then
                 elem(ie)%derived%vstar(:,:,:,k)=get_boomerang_velocities_gll(elem(ie), time_at(np1))
               else if (fvm_test_type == IDEAL_TEST_SOLIDBODY) then
                 elem(ie)%derived%vstar(:,:,:,k)=get_solidbody_velocities_gll(elem(ie), time_at(np1))
               else
                 call abortmp('Bad fvm_test_type in prim_step')
               end if
             end if
          enddo
          call remap1_nofilter(elem(ie)%derived%vstar,np,1,dp,dp_star)
       end do
    else
       ! do nothing
       ! for rsplit>0:  dynamics is also vertically lagrangian, so we do not need 
       ! to interpolate v(np1). 
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2D advection step
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------
!     call t_startf('fvm_depalg')

!     call fvm_mcgregordss(elem,fvm,nets,nete, hybrid, deriv, dt, 3)
    call fvm_rkdss(elem,fvm,nets,nete, hybrid, deriv, dt, 3)
!     call t_stopf('fvm_depalg')

!------------------------------------------------------------------------------------
    
    ! fvm departure calcluation should use vstar.
    ! from c(n0) compute c(np1):
    if (tracer_transport_type == TRACERTRANSPORT_FLUXFORM_FVM) then
      call cslam_runflux(elem,fvm,hybrid,deriv,dt,tl,nets,nete)
    else if (tracer_transport_type == TRACERTRANSPORT_LAGRANGIAN_FVM) then
      call cslam_runairdensity(elem,fvm,hybrid,deriv,dt,tl,nets,nete)
    else
      call abortmp('Bad tracer_transport_type in Prim_Advec_Tracers_fvm')
    end if

    call t_stopf('prim_advec_tracers_fvm')
  end subroutine Prim_Advec_Tracers_fvm



!=================================================================================================!

  subroutine Prim_Advec_Tracers_remap( elem , deriv , hvcoord , flt , hybrid , dt , tl , nets , nete )
    use control_mod   , only : use_semi_lagrange_transport
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (filter_t)      , intent(in   ) :: flt
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=real_kind) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete


  if (use_semi_lagrange_transport) then
    call Prim_Advec_Tracers_remap_ALE( elem , deriv ,                 hybrid , dt , tl , nets , nete )

  else
    call Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , flt , hybrid , dt , tl , nets , nete )
  end if
  end subroutine Prim_Advec_Tracers_remap



subroutine  Prim_Advec_Tracers_remap_ALE( elem , deriv , hybrid , dt , tl , nets , nete )
  use coordinate_systems_mod, only : cartesian3D_t, cartesian2D_t
  use dimensions_mod,         only : max_neigh_edges
  use edge_mod,               only : initghostbuffer3D, ghostVpack_unoriented, ghostVunpack_unoriented
  use bndry_mod,              only : ghost_exchangevfull
  use interpolate_mod,        only : interpolate_tracers, minmax_tracers
  use control_mod   ,         only : qsplit
  use global_norms_mod,       only: wrap_repro_sum
  use parallel_mod,           only: global_shared_buf, global_shared_sum



  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  type (derivative_t)  , intent(in   ) :: deriv
  type (hybrid_t)      , intent(in   ) :: hybrid
  real(kind=real_kind) , intent(in   ) :: dt
  type (TimeLevel_t)   , intent(in   ) :: tl
  integer              , intent(in   ) :: nets
  integer              , intent(in   ) :: nete


  type(cartesian3D_t)                           :: dep_points  (np,np)
  integer                                       :: elem_indexes(np,np)
  type(cartesian2D_t)                           :: para_coords (np,np)
  real(kind=real_kind)                          :: Que         (np,np,nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: Que_t       (np,np,nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: minq        (np,np,nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: maxq        (np,np,nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: f                      (qsize)
  real(kind=real_kind)                          :: g                      (qsize)
  real(kind=real_kind)                          :: mass              (nlev,qsize)
  real(kind=real_kind)                          :: rho         (np,np,nlev,      nets:nete)

  real(kind=real_kind)                          :: neigh_q     (np,np,qsize,max_neigh_edges+1)
  real(kind=real_kind)                          :: u           (np,np,qsize)

  integer                                       :: i,j,k,l,n,q,ie,kptr, n0_qdp, np1_qdp
  integer                                       :: num_neighbors
  logical,save :: firstcall = .true.

  call t_barrierf('Prim_Advec_Tracers_remap_ALE', hybrid%par%comm)
  call t_startf('Prim_Advec_Tracers_remap_ALE')

  if (firstcall) then
!OMP BARRIER
!OMP MASTER
     firstcall=.false.
     call initghostbuffer3D(ghostbuf_tr,nlev*qsize,np)
!OMP END MASTER
!OMP BARRIER
  endif

  call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)

  ! compute displacements for departure grid
  ! store in elem%derived%vstar
  call ALE_RKdss (elem, nets, nete, hybrid, deriv, dt, tl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  run ghost exchange to get global ID of all neighbors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  call t_startf('Prim_Advec_Tracers_remap_ALE_ghost_exchange')
  do ie=nets,nete
     kptr=0
     do k=1,nlev
     do q=1,qsize
        ! note: pack so that tracers per level are contiguous so we can unpack into
        ! array neigh_q()
        elem(ie)%state%Q(:,:,k,q) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp) / elem(ie)%derived%dp(:,:,k)
        call ghostVpack_unoriented(ghostbuf_tr, elem(ie)%state%Q(:,:,k,q),np,1,kptr,elem(ie)%desc)
        kptr=kptr+1
     enddo
     enddo
  end do
  call ghost_exchangeVfull(hybrid,ghostbuf_tr)
  do ie=nets,nete
     num_neighbors = elem(ie)%desc%actual_neigh_edges+1
     do k=1,nlev

        call t_startf('Prim_Advec_Tracers_remap_ALE_departure_from_gll')
        ! find departure points
        call ALE_departure_from_gll     (dep_points, elem(ie)%derived%vstar(:,:,:,k), elem(ie), dt)
        call t_stopf('Prim_Advec_Tracers_remap_ALE_departure_from_gll')

        call t_startf('Prim_Advec_Tracers_remap_ALE_elems_with_dep_points')
        ! find element containing departure point
        call ALE_elems_with_dep_points  (elem_indexes, dep_points, num_neighbors, elem(ie)%desc%neigh_corners)
        call t_stopf('Prim_Advec_Tracers_remap_ALE_elems_with_dep_points')

        call t_startf('Prim_Advec_Tracers_remap_ALE_parametric_coords')
        ! compute the parametric points
        call ALE_parametric_coords      (para_coords, elem_indexes, dep_points, num_neighbors, elem(ie)%desc%neigh_corners)
        call t_stopf('Prim_Advec_Tracers_remap_ALE_parametric_coords')

        call t_startf('Prim_Advec_Tracers_remap_ALE_ghostVunpack_unoriented')
        ! for each level k, unpack all tracer neighbor data on that level
        kptr=(k-1)*qsize
        neigh_q=0
        u(:,:,:) = elem(ie)%state%Q(:,:,k,:)
        call ghostVunpack_unoriented (ghostbuf_tr, neigh_q, np, qsize, kptr, elem(ie)%desc, elem(ie)%GlobalId, u)
        call t_stopf('Prim_Advec_Tracers_remap_ALE_ghostVunpack_unoriented')

        call t_startf('Prim_Advec_Tracers_remap_ALE_interpolate_tracers_minmax')
        do i=1,np
        do j=1,np
          ! interpolate tracers to deperature grid
          call interpolate_tracers     (para_coords(i,j), neigh_q(:,:,:,elem_indexes(i,j)),f)
          elem(ie)%state%Q(i,j,k,:) = f 

          call minmax_tracers          (para_coords(i,j), neigh_q(:,:,:,elem_indexes(i,j)),f,g)
          minq(i,j,k,:,ie) = f 
          maxq(i,j,k,:,ie) = g 
        enddo
        enddo
        call t_stopf('Prim_Advec_Tracers_remap_ALE_interpolate_tracers_minmax')
     end do
  end do

  call t_stopf('Prim_Advec_Tracers_remap_ALE_ghost_exchange')
  ! compute original mass, at tl_1%n0
  ! NOTE: global_shared_buf was allocated outside threaded region probably assuming n<=10.
  do ie=nets,nete
    n=0
    do k=1,nlev
    do q=1,qsize
      n=n+1
      global_shared_buf(ie,n) = 0
      do j=1,np
        global_shared_buf(ie,n) = global_shared_buf(ie,n) + DOT_PRODUCT(elem(ie)%state%Qdp(:,j,k,q,n0_qdp),elem(ie)%spheremp(:,j))
      end do
    end do
    end do
  end do
  call wrap_repro_sum(nvars=n, comm=hybrid%par%comm)
  n=0
  do k=1,nlev
  do q=1,qsize
    n=n+1
    mass(k,q) = global_shared_sum(n)
  enddo
  enddo

  do ie=nets,nete
  do k=1,nlev
    rho(:,:,k,ie) = elem(ie)%spheremp(:,:)*elem(ie)%state%dp3d(:,:,k,tl%np1)
  end do
  end do

  do ie=nets,nete
    Que_t(:,:,:,:,ie) = elem(ie)%state%Q(:,:,:,:)
  end do

  call t_startf('Prim_Advec_Tracers_remap_ALE_Cobra')
  call Cobra_SLBQP(Que, Que_t, rho, minq, maxq, mass, hybrid, nets, nete)
  call t_stopf('Prim_Advec_Tracers_remap_ALE_Cobra')

  do ie=nets,nete
    elem(ie)%state%Q(:,:,:,:) =  Que(:,:,:,:,ie)
  end do

  do ie=nets,nete
  do k=1,nlev
  do q=1,qsize
     ! note: pack so that tracers per level are contiguous so we can unpack into
     ! array neigh_q()
     elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%state%Q(:,:,k,q) * elem(ie)%state%dp3d(:,:,k,tl%np1)
  enddo
  enddo
  end do

! do ie=nets,nete
!    do q = 1 , qsize
!       do k = 1 , nlev    !  Potential loop inversion (AAM)
!          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%spheremp(:,:)* elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
!       enddo
!    enddo
!     call edgeVpack(edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
! enddo
! call bndry_exchangeV( hybrid , edgeAdv    )
! do ie = nets , nete
!    call edgeVunpack( edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
!    do q = 1 , qsize
!       do k = 1 , nlev    !  Potential loop inversion (AAM)
!          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
!       enddo
!    enddo
! enddo


  call t_stopf('Prim_Advec_Tracers_remap_ALE')


end subroutine Prim_Advec_Tracers_remap_ALE

subroutine VDOT(rp,Que,rho,mass,hybrid,nets,nete)
  use parallel_mod,        only: global_shared_buf, global_shared_sum
  use global_norms_mod,    only: wrap_repro_sum

  implicit none
  integer             , intent(in)              :: nets
  integer             , intent(in)              :: nete
  real(kind=real_kind), intent(out)             :: rp                (nlev,qsize)
  real(kind=real_kind), intent(in)              :: Que         (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: rho         (np*np,nlev,      nets:nete)
  real(kind=real_kind), intent(in)              :: mass              (nlev,qsize)
  type (hybrid_t)     , intent(in)              :: hybrid

  integer                                       :: k,n,q,ie

  global_shared_buf = 0 
  do ie=nets,nete
    n=0
    do q=1,qsize
    do k=1,nlev
      n=n+1
      global_shared_buf(ie,n) = global_shared_buf(ie,n) + DOT_PRODUCT(Que(:,k,q,ie), rho(:,k,ie))
    end do
    end do
  end do

  call wrap_repro_sum(nvars=n, comm=hybrid%par%comm)

  n=0
  do q=1,qsize
  do k=1,nlev
    n=n+1
    rp(k,q) = global_shared_sum(n) - mass(k,q)
  enddo
  enddo
  
end subroutine VDOT

subroutine Cobra_SLBQP(Que, Que_t, rho, minq, maxq, mass, hybrid, nets, nete) 

  use parallel_mod,        only: global_shared_buf, global_shared_sum
  use global_norms_mod,    only: wrap_repro_sum

  implicit none
  integer             , intent(in)              :: nets
  integer             , intent(in)              :: nete
  real(kind=real_kind), intent(out)             :: Que         (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: Que_t       (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: rho         (np*np,nlev,      nets:nete)
  real(kind=real_kind), intent(in)              :: minq        (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: maxq        (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: mass              (nlev,qsize)
  type (hybrid_t)     , intent(in)              :: hybrid

  integer,                            parameter :: max_clip = 10
  real(kind=real_kind),               parameter :: eta = 1D-08           
  real(kind=real_kind),               parameter :: hfd = 1D-10             
  real(kind=real_kind)                          :: lambda_p          (nlev,qsize)
  real(kind=real_kind)                          :: lambda_c          (nlev,qsize)
  real(kind=real_kind)                          :: rp                (nlev,qsize)
  real(kind=real_kind)                          :: rc                (nlev,qsize)
  real(kind=real_kind)                          :: rd                (nlev,qsize)
  real(kind=real_kind)                          :: alpha             (nlev,qsize)
  integer                                       :: j,k,n,q,ie
  integer                                       :: nclip

  nclip = 0

  Que(:,:,:,:) = Que_t(:,:,:,:)

  Que = MIN(MAX(Que,minq),maxq)

  call VDOT(rp,Que,rho,mass,hybrid,nets,nete)
  nclip = nclip + 1

  if (MAXVAL(ABS(rp)).lt.eta) return

  do ie=nets,nete
  do q=1,qsize
  do k=1,nlev
     Que(:,k,q,ie) = hfd * rho(:,k,ie) + Que_t(:,k,q,ie)
  enddo
  enddo
  enddo

  Que = MIN(MAX(Que,minq),maxq)

  call VDOT(rc,Que,rho,mass,hybrid,nets,nete)

  rd = rc-rp
  if (MAXVAL(ABS(rd)).eq.0) return 
  
  alpha = 0
  WHERE (rd.ne.0) alpha = hfd / rd 

  lambda_p = 0
  lambda_c =  -alpha*rp

  do while (MAXVAL(ABS(rc)).gt.eta .and. nclip.lt.max_clip)

    do ie=nets,nete
    do q=1,qsize
    do k=1,nlev
       Que(:,k,q,ie) = (lambda_c(k,q) + hfd) * rho(:,k,ie) + Que_t(:,k,q,ie)
    enddo
    enddo
    enddo
    Que = MIN(MAX(Que,minq),maxq)

    call VDOT(rc,Que,rho,mass,hybrid,nets,nete)
    nclip = nclip + 1

    rd = rp-rc

    if (MAXVAL(ABS(rd)).eq.0) exit

    alpha = 0
    WHERE (rd.ne.0) alpha = (lambda_p - lambda_c) / rd 

    rp       = rc
    lambda_p = lambda_c

    lambda_c = lambda_c -  alpha * rc

  enddo
end subroutine Cobra_SLBQP

! ----------------------------------------------------------------------------------!
!SUBROUTINE ALE_RKDSS-----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR, 06. December 2012
!
! DESCRIPTION: ! create a runge kutta taylor serios mixture to calculate the departure grid
!
! CALLS:
! INPUT:
!
! OUTPUT:
!-----------------------------------------------------------------------------------!

! this will calculate the velocity at time t+1/2  along the trajectory s(t) given the velocities
! at the GLL points at time t and t+1 using a second order time accurate formulation.


subroutine ALE_RKdss(elem, nets, nete, hy, deriv, dt, tl)

  use derivative_mod,         only : derivative_t, ugradv_sphere
  use edge_mod      ,         only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod     ,         only : bndry_exchangev
  use kinds,                  only : real_kind
  use hybrid_mod,             only : hybrid_t
  use element_mod,            only : element_t
  use dimensions_mod,         only : np, nlev

  implicit none

  type (element_t)     , intent(inout)             :: elem(:)
  integer              , intent(in   )             :: nets
  integer              , intent(in   )             :: nete
  type (hybrid_t)      , intent(in)                :: hy ! distributed parallel structure (shared)
  type (derivative_t)  , intent(in)                :: deriv ! derivative struct
  real (kind=real_kind), intent(in)                :: dt ! timestep
  type (TimeLevel_t)   , intent(in)                :: tl

  integer                                          :: ie, k
  real (kind=real_kind), dimension(np,np,2)        :: vtmp
  integer :: np1

  np1 = tl%np1


! RK-SSP 2 stage 2nd order:
!     x*(t+1) = x(t) + U(x(t),t) dt
!     x(t+1) = x(t) +  1/2 ( U(x*(t+1),t+1) + U(x(t),t) ) dt
! apply taylor series:
!  U(x*(t+1),t+1) = U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1)
!
! x(t+1) = x(t) +  1/2 ( U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1) + U(x(t),t) ) dt
! (x(t+1) - x(t)) / dt =  1/2 ( U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1) + U(x(t),t) )
! (x(t+1) - x(t)) / dt =  1/2 ( U(x(t),t+1) + U(x(t),t) + (x*(t+1)-x(t)) gradU(x(t),t+1) )
! (x(t+1) - x(t)) / dt =  1/2 ( U(x(t),t+1) + U(x(t),t) + U(x(t),t) dt  gradU(x(t),t+1) )
!
!
!  (x(t+1)-x(t))/dt =  1/2(U(x(t),t+1) + U(x(t),t) + dt U(x(t),t) gradU(x(t),t+1))
!
! suppose dt = -ts (we go backward)
!  (x(t-ts)-x(t))/-ts =  1/2( U(x(t),t-ts)+U(x(t),t)) - ts 1/2 U(x(t),t) gradU(x(t),t-ts)
!
!  x(t-ts) = x(t)) -ts * [ 1/2( U(x(t),t-ts)+U(x(t),t)) - ts 1/2 U(x(t),t) gradU(x(t),t-ts) ]
!
!    !------------------------------------------------------------------------------------

  do ie=nets,nete
     ! vstarn0 = U(x,t)
     ! vstar   = U(x,t+1)
    do k=1,nlev
      vtmp(:,:,:)=ugradv_sphere(elem(ie)%state%v(:,:,:,k,np1), elem(ie)%derived%vstar(:,:,:,k),deriv,elem(ie))

      elem(ie)%derived%vstar(:,:,:,k) = &
           (elem(ie)%state%v(:,:,:,k,np1) + elem(ie)%derived%vstar(:,:,:,k))/2 - dt*vtmp(:,:,:)/2

      elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%spheremp(:,:)
      elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%spheremp(:,:)
    enddo
    call edgeVpack(edgeveloc,elem(ie)%derived%vstar,2*nlev,0,elem(ie)%desc)
  enddo
  call bndry_exchangeV(hy,edgeveloc)
  do ie=nets,nete
    call edgeVunpack(edgeveloc,elem(ie)%derived%vstar,2*nlev,0,elem(ie)%desc)
    do k=1, nlev
      elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%rspheremp(:,:)
      elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%rspheremp(:,:)
    end do
  end do
end subroutine ALE_RKdss

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
subroutine ALE_departure_from_gll(acart, vstar, elem, dt)
  use physical_constants,     only : rearth
  use coordinate_systems_mod, only : spherical_polar_t, cartesian3D_t, change_coordinates
  use time_mod,               only : timelevel_t
  use element_mod,            only : element_t
  use kinds,                  only : real_kind
  use dimensions_mod,         only : np

  implicit none

  type(cartesian3D_t)     ,intent(out)  :: acart(np,np)
  real (kind=real_kind)   ,intent(in)   :: vstar(np,np,2)
  type (element_t)        ,intent(in)   :: elem
  real (kind=real_kind)   ,intent(in)   :: dt

  integer                               :: i,j

  real (kind=real_kind)                 :: uxyz (np,np,3)

   ! convert velocity from lat/lon to cartesian 3D

  do i=1,3
     ! Summing along the third dimension is a sum over components for each point.
     ! (This is just a faster way of doing a dot product for each grid point,
     ! since reindexing the inputs to use the intrinsic effectively would be
     ! just asking for trouble.)
     uxyz(:,:,i)=sum( elem%vec_sphere2cart(:,:,i,:)*vstar(:,:,:) ,3)
  end do
  ! interpolate velocity to fvm nodes
  ! compute departure point
  ! crude, 1st order accurate approximation.  to be improved
  do i=1,np
     do j=1,np
        acart(i,j) = change_coordinates(elem%spherep(i,j)) 
        acart(i,j)%x = acart(i,j)%x - dt*uxyz(i,j,1)/rearth
        acart(i,j)%y = acart(i,j)%y - dt*uxyz(i,j,2)/rearth
        acart(i,j)%z = acart(i,j)%z - dt*uxyz(i,j,3)/rearth
     enddo
  enddo

end subroutine ALE_departure_from_gll




subroutine ALE_elems_with_dep_points (elem_indexes, dep_points, num_neighbors, ngh_corners)

  use element_mod,            only : element_t
  use dimensions_mod,         only : np
  use coordinate_systems_mod, only : cartesian3D_t, change_coordinates
  use interpolate_mod,        only : point_inside_quad

  implicit none

  ! The ngh_corners array is a list of corners of both elem and all of it's
  ! neighor elements all sorted by global id.
  integer              , intent(in)                :: num_neighbors
  type(cartesian3D_t),intent(in)                   :: ngh_corners(4,num_neighbors)
  integer              , intent(out)               :: elem_indexes(np,np)
  type(cartesian3D_t)  , intent(in)                :: dep_points(np,np)

  integer                                          :: i,j,n
  logical                                          :: inside

  elem_indexes = -1
  do i=1,np
    do j=1,np
! Just itererate the neighbors in global id order to get the same result on every processor.
      do n = 1, num_neighbors
! Mark Taylor's handy dandy point_inside_gc check.
       inside = point_inside_quad (ngh_corners(:,n), dep_points(i,j))
       if (inside) then
         elem_indexes(i,j) = n
         exit
       end if
     end do
    end do
  end do

  if (MINVAL(elem_indexes(:,:))==-1) then
    write (*,*) __FILE__,__LINE__,"Aborting because point not found in neighbor list. Info:"
    do i=1,np
      do j=1,np
        if (elem_indexes(i,j)==-1) then
          write (*,*)   " departure point ",dep_points(i,j)
          do n = 1, num_neighbors
            write (*,*) " quad checked    ",ngh_corners(1,n)
            write (*,*) "                 ",ngh_corners(2,n)
            write (*,*) "                 ",ngh_corners(3,n)
            write (*,*) "                 ",ngh_corners(4,n)
            write (*,*)
          end do
          exit
        end if
      end do
    end do
    call abortmp("ERROR elems_with_dep_points: Can't find departure grid. Time step too long?")
  end if
end subroutine ALE_elems_with_dep_points

subroutine  shape_fcn_deriv(dNds, pc)
  real (kind=real_kind), intent(out) :: dNds(4,2)
  real (kind=real_kind), intent(in)  ::  pc(2)
 
  dNds(1, 1) = - 0.25 * (1.0 - pc(2))
  dNds(1, 2) = - 0.25 * (1.0 - pc(1))

  dNds(2, 1) =   0.25 * (1.0 - pc(2))
  dNds(2, 2) = - 0.25 * (1.0 + pc(1))

  dNds(3, 1) =   0.25 * (1.0 + pc(2))
  dNds(3, 2) =   0.25 * (1.0 + pc(1))

  dNds(4, 1) = - 0.25 * (1.0 + pc(2))
  dNds(4, 2) =   0.25 * (1.0 - pc(1))
end subroutine 

subroutine inv_2x2(A_inv, A)
  real (kind=real_kind), intent(out) :: A_inv(2,2)
  real (kind=real_kind), intent(in)  :: A    (2,2)
  real (kind=real_kind) :: det, denom

  det = A(1,1) * A(2,2) - A(2,1) * A(1,2)
  denom = 1/det
  ! inverse:
  A_inv(1,1) =  denom * A(2,2)  !  dxidx
  A_inv(2,1) = -denom * A(2,1)  !  detadx
  A_inv(1,2) = -denom * A(1,2)  !  dxidy
  A_inv(2,2) =  denom * A(1,1)  !  detady
end subroutine

subroutine formDSDX(dsdx, coord, dNds)

  real (kind=real_kind), intent(out) :: dsdx(2,3)
  real (kind=real_kind), intent(in)  :: coord(4,3), dNds(4,2)

  real (kind=real_kind)  ::     dxds(3,2)
  real (kind=real_kind)  ::      ata(2,2)
  real (kind=real_kind)  ::  ata_inv(2,2)

  dxds = MATMUL(TRANSPOSE(coord), dNds)

  !     dxds = | dxdxi   dxdeta |
  !            | dydxi   dydeta |
  !            | dzdxi   dzdeta |
  ata  = MATMUL(TRANSPOSE(dxds), dxds)
  call inv_2x2(ata_inv, ata)
  dsdx = MATMUL(ata_inv, TRANSPOSE(dxds))
  !     dsdx = |  dxidx   dxidy   dxidz |
  !            | detadx  detady  detadz |

end subroutine

subroutine shape_fcn(N, pc)
  real (kind=real_kind), intent(out) :: N(4)
  real (kind=real_kind), intent(in)  :: pc(2)

  ! shape function for each node evaluated at param_coords
  N(1) = 0.25 * (1.0 - pc(1)) * (1.0 - pc(2)) 
  N(2) = 0.25 * (1.0 + pc(1)) * (1.0 - pc(2)) 
  N(3) = 0.25 * (1.0 + pc(1)) * (1.0 + pc(2)) 
  N(4) = 0.25 * (1.0 - pc(1)) * (1.0 + pc(2)) 
end subroutine


function fcn(pc, coords) result(g)
  real (kind=real_kind), intent(in) :: pc(2), coords(4,3)

  real (kind=real_kind)            :: N(4), g(3)
  N = 0
  g = 0
  call shape_fcn(N,pc)
  g = MATMUL(TRANSPOSE(coords), N)
end function


function cartesian_parametric_coordinates(sphere, corners3D) result (ref)
  use coordinate_systems_mod, only : cartesian2d_t, cartesian3D_t, spherical_polar_t, spherical_to_cart
  implicit none
  type (spherical_polar_t), intent(in) :: sphere
  type (cartesian3D_t)    , intent(in) :: corners3D(4)  !x,y,z coords of element corners

  type (cartesian2D_t)                 :: ref

  integer,               parameter :: MAXIT = 20
  real (kind=real_kind), parameter :: TOL   = 1.0E-9
  integer,               parameter :: n     = 3

  type (cartesian3D_t)             :: cart
  real (kind=real_kind)            :: coords(4,3), dNds(4,2), dsdx(2,3)
  real (kind=real_kind)            :: p(3), pc(2), dx(3), g(3), ds(2)
  real (kind=real_kind)            :: dist, step                          
  
  integer                          :: i,j,k,iter
  do i=1,4                               
    coords(i,1) = corners3D(i)%x 
    coords(i,2) = corners3D(i)%y 
    coords(i,3) = corners3D(i)%z 
  end do

  pc = 0
  p  = 0
  cart = spherical_to_cart(sphere)

  p(1) = cart%x
  p(2) = cart%y
  p(3) = cart%z 

  dx   = 0
  ds   = 0
  dNds = 0
  dsdx = 0

  /*-------------------------------------------------------------------------*/

  ! Initial guess, center of element
  dist = 9999999.                         
  step = 9999999.                         
  iter = 0 

  do while  (TOL*TOL.lt.dist .and. iter.lt.MAXIT .and. TOL*TOL.lt.step)
    iter = iter + 1

    call shape_fcn_deriv  (dNds, pc)
    call formDSDX         (dsdx, coords, dNds)
    g = fcn(pc, coords)

    dx = g - p
    dist = DOT_PRODUCT(dx,dx)
    ds = MATMUL(dsdx, dx)
    pc = pc - ds
    step = DOT_PRODUCT(ds,ds)
  enddo

  ref%x = pc(1)
  ref%y = pc(2)
end function


subroutine  ALE_parametric_coords (parametric_coord, elem_indexes, dep_points, num_neighbors, ngh_corners)
  use coordinate_systems_mod, only : cartesian2d_t, cartesian3D_t, spherical_polar_t, change_coordinates, distance
  use interpolate_mod,        only : parametric_coordinates
  use dimensions_mod,         only : np

  implicit none

  type(cartesian2D_t)       , intent(out)       :: parametric_coord(np,np)
  type(cartesian3D_t)       , intent(in)        :: dep_points(np,np)
  integer                   , intent(in)        :: elem_indexes(np,np)
  integer                   , intent(in)        :: num_neighbors
  type(cartesian3D_t)       , intent(in)        :: ngh_corners(4,num_neighbors)

  type (spherical_polar_t)                      :: sphere(np,np)
  integer                                       :: i,j,n
! type(cartesian2D_t)                           :: parametric_test
  real(kind=real_kind)                          :: d

  do j=1,np
    sphere(:,j) = change_coordinates(dep_points(:,j))
  end do

  do i=1,np
    do j=1,np
      n = elem_indexes(i,j)
      ! Mark will fill in  parametric_coordinates for corners.
      parametric_coord(i,j) = cartesian_parametric_coordinates(sphere(i,j),ngh_corners(:,n))
!     parametric_test       = parametric_coordinates(sphere(i,j),ngh_corners(:,n))
    end do
  end do
end subroutine ALE_parametric_coords


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! forward-in-time 2 level vertically lagrangian step
!  this code takes a lagrangian step in the horizontal
! (complete with DSS), and then applies a vertical remap
!
! This routine may use dynamics fields at timelevel np1
! In addition, other fields are required, which have to be
! explicitly saved by the dynamics:  (in elem(ie)%derived struct)
!
! Fields required from dynamics: (in
!    omega_p   it will be DSS'd here, for later use by CAM physics
!              we DSS omega here because it can be done for "free"
!    Consistent mass/tracer-mass advection (used if subcycling turned on)
!       dp()   dp at timelevel n0
!       vn0()  mean flux  < U dp > going from n0 to np1
!
! 3 stage
!    Euler step from t     -> t+.5
!    Euler step from t+.5  -> t+1.0
!    Euler step from t+1.0 -> t+1.5
!    u(t) = u(t)/3 + u(t+2)*2/3
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , flt , hybrid , dt , tl , nets , nete )
    use perf_mod      , only : t_startf, t_stopf            ! _EXTERNAL
    use derivative_mod, only : divergence_sphere
    use control_mod   , only : vert_remap_q_alg, qsplit
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (filter_t)      , intent(in   ) :: flt
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=real_kind) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    real (kind=real_kind), dimension(np,np,2     ) :: gradQ
    real (kind=real_kind), dimension(np,np  ,nlev) :: dp_star
    real (kind=real_kind), dimension(np,np  ,nlev) :: dp_np1
    integer :: i,j,k,l,ie,q,nmin
    integer :: nfilt,rkstage,rhs_multiplier
    integer :: n0_qdp, np1_qdp

    call t_barrierf('sync_prim_advec_tracers_remap_k2', hybrid%par%comm)
    call t_startf('prim_advec_tracers_remap_rk2')
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp) !time levels for qdp are not the same
    rkstage = 3 !   3 stage RKSSP scheme, with optimal SSP CFL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! RK2 2D advection step
    ! note: stage 3 we take the oppertunity to DSS omega
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! use these for consistent advection (preserve Q=1)
    ! derived%vdp_ave        =  mean horiz. flux:   U*dp
    ! derived%eta_dot_dpdn    =  mean vertical velocity (used for remap)
    ! derived%omega_p         =  advection code will DSS this for the physics, but otherwise
    !                            it is not needed
    ! Also: save a copy of div(U dp) in derived%div(:,:,:,1), which will be DSS'd
    !       and a DSS'ed version stored in derived%div(:,:,:,2)
    do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k, gradQ)
#endif
      do k=1,nlev
        ! div( U dp Q),
        gradQ(:,:,1)=elem(ie)%derived%vn0(:,:,1,k)
        gradQ(:,:,2)=elem(ie)%derived%vn0(:,:,2,k)
        elem(ie)%derived%divdp(:,:,k) = divergence_sphere(gradQ,deriv,elem(ie))
      enddo
      elem(ie)%derived%divdp_proj(:,:,:) = elem(ie)%derived%divdp(:,:,:)
    enddo

    !rhs_multiplier is for obtaining dp_tracers at each stage:
    !dp_tracers(stage) = dp - rhs_multiplier*dt*divdp_proj
    rhs_multiplier = 0
    call euler_step( np1_qdp , n0_qdp  , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSdiv_vdp_ave , rhs_multiplier )

    rhs_multiplier = 1
    call euler_step( np1_qdp , np1_qdp , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSeta         , rhs_multiplier )

    rhs_multiplier = 2
    call euler_step( np1_qdp , np1_qdp , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSomega       , rhs_multiplier )

    !to finish the 2D advection step, we need to average the t and t+2 results to get a second order estimate for t+1.
    call qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Dissipation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( limiter_option == 8  ) then
      ! dissipation was applied in RHS.
    else
      call advance_hypervis_scalar(edgeadv,elem,hvcoord,hybrid,deriv,tl%np1,np1_qdp,nets,nete,dt)
    endif

    call t_stopf('prim_advec_tracers_remap_rk2')
  end subroutine prim_advec_tracers_remap_rk2

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )
#if USE_CUDA_FORTRAN
    use cuda_mod, only: qdp_time_avg_cuda
#endif
    implicit none
    type(element_t)     , intent(inout) :: elem(:)
    integer             , intent(in   ) :: rkstage , n0_qdp , np1_qdp , nets , nete , limiter_option
    real(kind=real_kind), intent(in   ) :: nu_p
    integer :: ie
#if USE_CUDA_FORTRAN
    call qdp_time_avg_cuda( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )
    return
#endif
    do ie=nets,nete
      elem(ie)%state%Qdp(:,:,:,1:qsize,np1_qdp) =               &
                   ( elem(ie)%state%Qdp(:,:,:,1:qsize,n0_qdp) + &
                     (rkstage-1)*elem(ie)%state%Qdp(:,:,:,1:qsize,np1_qdp) ) / rkstage
    enddo
  end subroutine qdp_time_avg

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  ! ===================================
  ! This routine is the basic foward
  ! euler component used to construct RK SSP methods
  !
  !           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! n0 can be the same as np1.
  !
  ! DSSopt = DSSeta or DSSomega:   also DSS eta_dot_dpdn or omega
  !
  ! ===================================
  use kinds          , only : real_kind
  use dimensions_mod , only : np, npdg, nlev
  use hybrid_mod     , only : hybrid_t
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod       , only : edgevpack, edgevunpack
  use bndry_mod      , only : bndry_exchangev
  use hybvcoord_mod  , only : hvcoord_t
#if USE_CUDA_FORTRAN
  use cuda_mod, only: euler_step_cuda
#endif
  implicit none
  integer              , intent(in   )         :: np1_qdp, n0_qdp
  real (kind=real_kind), intent(in   )         :: dt
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer              , intent(in   )         :: rhs_multiplier

  ! local
  real(kind=real_kind), dimension(np,np                       ) :: divdp, dpdiss
  real(kind=real_kind), dimension(np,np,2                     ) :: gradQ
  real(kind=real_kind), dimension(np,np,2,nlev                ) :: Vstar
  real(kind=real_kind), dimension(np,np  ,nlev                ) :: Qtens
  real(kind=real_kind), dimension(np,np  ,nlev                ) :: dp,dp_star
  real(kind=real_kind), dimension(np,np  ,nlev,qsize,nets:nete) :: Qtens_biharmonic
  real(kind=real_kind), pointer, dimension(:,:,:)               :: DSSvar
  real(kind=real_kind) :: dp0
  integer :: ie,q,i,j,k
  integer :: rhs_viss = 0

  if ( npdg > 0 ) then
    call euler_step_dg( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
    return
  endif
#if USE_CUDA_FORTRAN
  call euler_step_cuda( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  return
#endif
! call t_barrierf('sync_euler_step', hybrid%par%comm)
!   call t_startf('euler_step')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   compute Q min/max values for lim8
  !   compute biharmonic mixing term f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rhs_viss = 0
  if ( limiter_option == 8  ) then
    ! when running lim8, we also need to limit the biharmonic, so that term needs
    ! to be included in each euler step.  three possible algorithms here:
    ! 1) most expensive:
    !     compute biharmonic (which also computes qmin/qmax) during all 3 stages
    !     be sure to set rhs_viss=1
    !     cost:  3 biharmonic steps with 3 DSS
    !
    ! 2) cheapest:
    !     compute biharmonic (which also computes qmin/qmax) only on first stage
    !     be sure to set rhs_viss=3
    !     reuse qmin/qmax for all following stages (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps with 1 DSS
    !     main concern:  viscosity
    !
    ! 3)  compromise:
    !     compute biharmonic (which also computes qmin/qmax) only on last stage
    !     be sure to set rhs_viss=3
    !     compute qmin/qmax directly on first stage
    !     reuse qmin/qmax for 2nd stage stage (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps, 2 DSS
    !
    !  NOTE  when nu_p=0 (no dissipation applied in dynamics to dp equation), we should
    !        apply dissipation to Q (not Qdp) to preserve Q=1
    !        i.e.  laplace(Qdp) ~  dp0 laplace(Q)
    !        for nu_p=nu_q>0, we need to apply dissipation to Q * diffusion_dp
    !
    ! initialize dp, and compute Q from Qdp (and store Q in Qtens_biharmonic)
    do ie = nets , nete
      ! add hyperviscosity to RHS.  apply to Q at timelevel n0, Qdp(n0)/dp
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k, q)
#endif
      do k = 1 , nlev    !  Loop index added with implicit inversion (AAM)
        dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - rhs_multiplier*dt*elem(ie)%derived%divdp_proj(:,:,k)
        do q = 1 , qsize
          Qtens_biharmonic(:,:,k,q,ie) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp)/dp(:,:,k)
        enddo
      enddo
    enddo

    ! compute element qmin/qmax
    if ( rhs_multiplier == 0 ) then
      do ie = nets , nete
        do k = 1 , nlev
          do q = 1 , qsize
            qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
            qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
            qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
          enddo
        enddo
      enddo
      ! update qmin/qmax based on neighbor data for lim8
      call neighbor_minmax(elem,hybrid,edgeAdvQ2,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
    endif

    ! lets just reuse the old neighbor min/max, but update based on local data
    if ( rhs_multiplier == 1 ) then
      do ie = nets , nete
        do k = 1 , nlev    !  Loop index added with implicit inversion (AAM)
          do q = 1 , qsize
            qmin(k,q,ie)=min(qmin(k,q,ie),minval(Qtens_biharmonic(:,:,k,q,ie)))
            qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
            qmax(k,q,ie)=max(qmax(k,q,ie),maxval(Qtens_biharmonic(:,:,k,q,ie)))
          enddo
        enddo
      enddo
    endif

    ! get niew min/max values, and also compute biharmonic mixing term
    if ( rhs_multiplier == 2 ) then
      rhs_viss = 3
      ! compute element qmin/qmax
      do ie = nets , nete
        do k = 1  ,nlev
          do q = 1 , qsize
            qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
            qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
            qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
          enddo
        enddo
      enddo
      ! two scalings depending on nu_p:
      ! nu_p=0:    qtens_biharmonic *= dp0                   (apply viscsoity only to q)
      ! nu_p>0):   qtens_biharmonc *= elem()%psdiss_ave      (for consistency, if nu_p=nu_q)
      if ( nu_p > 0 ) then
        do ie = nets , nete
#if (defined COLUMN_OPENMP)
          !$omp parallel do private(k, q, dp0, dpdiss)
#endif
          do k = 1 , nlev
            dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                  ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
#if 0
            dpdiss(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%derived%psdiss_ave(:,:)
#else
            dpdiss(:,:) = elem(ie)%derived%dpdiss_ave(:,:,k)
#endif
            do q = 1 , qsize
              ! NOTE: divide by dp0 since we multiply by dp0 below
              Qtens_biharmonic(:,:,k,q,ie)=Qtens_biharmonic(:,:,k,q,ie)*dpdiss(:,:)/dp0
            enddo
          enddo
        enddo
      endif
      call biharmonic_wk_scalar_minmax( elem , qtens_biharmonic , deriv , edgeAdvQ3 , hybrid , &
           nets , nete , qmin(:,:,nets:nete) , qmax(:,:,nets:nete) )
      do ie = nets , nete
#if (defined COLUMN_OPENMP)
        !$omp parallel do private(k, q, dp0)
#endif
        do k = 1 , nlev    !  Loop inversion (AAM)
          dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
          do q = 1 , qsize
            ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
            qtens_biharmonic(:,:,k,q,ie) = &
                     -rhs_viss*dt*nu_q*dp0*Qtens_biharmonic(:,:,k,q,ie) / elem(ie)%spheremp(:,:)
          enddo
        enddo
      enddo
    endif
  endif  ! compute biharmonic mixing term and qmin/qmax


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie = nets , nete
    ! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
    ! all zero so we only have to DSS 1:nlev
    if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
    if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
    if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

    ! Compute velocity used to advance Qdp
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
    do k = 1 , nlev    !  Loop index added (AAM)
      ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
      ! but that's ok because rhs_multiplier=0 on the first stage:
      dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - rhs_multiplier * dt * elem(ie)%derived%divdp_proj(:,:,k)
      Vstar(:,:,1,k) = elem(ie)%derived%vn0(:,:,1,k) / dp(:,:,k)
      Vstar(:,:,2,k) = elem(ie)%derived%vn0(:,:,2,k) / dp(:,:,k)
    enddo

    ! advance Qdp
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k,gradQ,dp_star,qtens,dpdiss)
#endif
    do q = 1 , qsize
      do k = 1 , nlev  !  dp_star used as temporary instead of divdp (AAM)
        ! div( U dp Q),
        gradQ(:,:,1) = Vstar(:,:,1,k) * elem(ie)%state%Qdp(:,:,k,q,n0_qdp)
        gradQ(:,:,2) = Vstar(:,:,2,k) * elem(ie)%state%Qdp(:,:,k,q,n0_qdp)
        dp_star(:,:,k) = divergence_sphere( gradQ , deriv , elem(ie) )
        Qtens(:,:,k) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp) - dt * dp_star(:,:,k)
        ! optionally add in hyperviscosity computed above:
        if ( rhs_viss /= 0 ) Qtens(:,:,k) = Qtens(:,:,k) + Qtens_biharmonic(:,:,k,q,ie)
      enddo

      if ( limiter_option == 8) then
        do k = 1 , nlev  ! Loop index added (AAM)
          ! UN-DSS'ed dp at timelevel n0+1:
          dp_star(:,:,k) = dp(:,:,k) - dt * elem(ie)%derived%divdp(:,:,k)
          if ( nu_p > 0 .and. rhs_viss /= 0 ) then
            ! add contribution from UN-DSS'ed PS dissipation
!            dpdiss(:,:) = ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * elem(ie)%derived%psdiss_biharmonic(:,:)
            dpdiss(:,:) = elem(ie)%derived%dpdiss_biharmonic(:,:,k)
            dp_star(:,:,k) = dp_star(:,:,k) - rhs_viss * dt * nu_q * dpdiss(:,:) / elem(ie)%spheremp(:,:)
          endif
        enddo
        ! apply limiter to Q = Qtens / dp_star
        call limiter_optim_iter_full( Qtens(:,:,:) , elem(ie)%spheremp(:,:) , qmin(:,q,ie) , &
                                      qmax(:,q,ie) , dp_star(:,:,:) )
      endif


      if ( limiter_option==84) then
        do k = 1 , nlev  ! Loop index added (AAM)
           qmin(k,1,1)=minval(elem(ie)%state%Qdp(:,:,k,q,n0_qdp)/dp(:,:,k))
           qmin(k,1,1)=max(qmin(k,1,1),0d0)
           qmax(k,1,1)=maxval(elem(ie)%state%Qdp(:,:,k,q,n0_qdp)/dp(:,:,k))
           
           ! UN-DSS'ed dp at timelevel n0+1:
           dp_star(:,:,k) = dp(:,:,k) - dt * elem(ie)%derived%divdp(:,:,k)
        enddo
        call limiter2d_minmax(Qtens,dp_star,elem(ie)%spheremp,qmin(:,1,1),qmax(:,1,1))
      endif



      ! apply mass matrix, overwrite np1 with solution:
      ! dont do this earlier, since we allow np1_qdp == n0_qdp
      ! and we dont want to overwrite n0_qdp until we are done using it
      do k = 1 , nlev
        elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%spheremp(:,:) * Qtens(:,:,k)
      enddo

      if ( limiter_option == 4 ) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        ! sign-preserving limiter, applied after mass matrix
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        call limiter2d_zero( elem(ie)%state%Qdp(:,:,:,q,np1_qdp))
      endif
    enddo

    if ( DSSopt == DSSno_var ) then
      call edgeVpack(edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
    else
      call edgeVpack(edgeAdv_p1 , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
      ! also DSS extra field
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
      do k = 1 , nlev
        DSSvar(:,:,k) = elem(ie)%spheremp(:,:) * DSSvar(:,:,k)
      enddo
      call edgeVpack( edgeAdv_p1 , DSSvar(:,:,1:nlev) , nlev , nlev*qsize , elem(ie)%desc )
    endif
  enddo

  if ( DSSopt == DSSno_var ) then
    call bndry_exchangeV( hybrid , edgeAdv    )
  else
    call bndry_exchangeV( hybrid , edgeAdv_p1 )
  endif

  do ie = nets , nete
    if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
    if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
    if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

    if ( DSSopt == DSSno_var ) then
      call edgeVunpack( edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,q)
#endif
      do q = 1 , qsize
        do k = 1 , nlev    !  Potential loop inversion (AAM)
          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
        enddo
      enddo
    else
      call edgeVunpack( edgeAdv_p1 , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
#if (defined COLUMN_OPENMP)
      !$omp parallel do private(q,k)
#endif
      do q = 1 , qsize
        do k = 1 , nlev    !  Potential loop inversion (AAM)
          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
        enddo
      enddo
      call edgeVunpack( edgeAdv_p1 , DSSvar(:,:,1:nlev) , nlev , qsize*nlev , elem(ie)%desc )

      do k = 1 , nlev
        DSSvar(:,:,k) = DSSvar(:,:,k) * elem(ie)%rspheremp(:,:)
      enddo
    endif
  enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
!   call t_stopf('euler_step')
  end subroutine euler_step

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine euler_step_dg(np1_qdp, n0_qdp, dt,elem,hvcoord,hybrid,deriv,nets,nete,&
      DSSopt,rhs_multiplier)
  ! ===================================
  ! This routine is the basic foward
  ! euler component used to construct RK SSP methods
  !
  !           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! n0 can be the same as np1.
  !
  ! DSSopt = DSSeta or DSSomega:   also DSS eta_dot_dpdn or omega
  !
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : np, npdg, nlev
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, divergence_sphere_wk, edge_flux_u_cg, gll_to_dgmodal, dgmodal_to_gll
  use edge_mod, only : edgevpack, edgevunpack, edgedgvunpack
  use bndry_mod, only : bndry_exchangev
  use hybvcoord_mod, only : hvcoord_t

  implicit none
  integer :: np1_qdp, n0_qdp, nets, nete, DSSopt, rhs_multiplier
  real (kind=real_kind), intent(in)  :: dt

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv

  ! local
  real (kind=real_kind), dimension(np,np)    :: divdp
  real (kind=real_kind), dimension(npdg,npdg)    :: pshat
  real (kind=real_kind), dimension(0:np+1,0:np+1,nlev,qsize)    :: qedges
  real (kind=real_kind), dimension(np,np,2)    :: vtemp
  real(kind=real_kind), dimension(np,np,nlev) :: dp,dp_star
  real(kind=real_kind), dimension(np,np,2,nlev) :: Vstar
  real (kind=real_kind), pointer, dimension(:,:,:)   :: DSSvar
! nelemd

  real(kind=real_kind) :: dp0
  integer :: ie,q,i,j,k
  integer :: rhs_viss=0

  call t_barrierf('sync_euler_step_dg', hybrid%par%comm)
  call t_startf('euler_step_dg')


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   compute Q min/max values for lim8
  !   compute biharmonic mixing term f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rhs_viss=0
  if (limiter_option /= 4) then
     call abortmp('only limiter_opiton=4 supported for dg advection')
     ! todo:  we need to track a 'dg' mass, and use that to back out Q
     ! then compute Qmin/Qmax here
  endif



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete

     ! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
     ! all zero so we only have to DSS 1:nlev
     if ( DSSopt == DSSeta) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
     if ( DSSopt == DSSomega) DSSvar => elem(ie)%derived%omega_p(:,:,:)
     if ( DSSopt == DSSdiv_vdp_ave) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

     if(DSSopt==DSSno_var)then
	call edgeVpack(edgeAdv,elem(ie)%state%Qdp(:,:,:,:,n0_qdp),nlev*qsize,0,elem(ie)%desc)
     else
	call edgeVpack(edgeAdv_p1,elem(ie)%state%Qdp(:,:,:,:,n0_qdp),nlev*qsize,0,elem(ie)%desc)
	! also DSS extra field
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
	do k=1,nlev
	    DSSvar(:,:,k) = elem(ie)%spheremp(:,:)*DSSvar(:,:,k)
	enddo
	call edgeVpack(edgeAdv_p1,DSSvar(:,:,1:nlev),nlev,nlev*qsize,elem(ie)%desc)
     endif

  end do

  if(DSSopt==DSSno_var)then
     call bndry_exchangeV(hybrid,edgeAdv)
  else
     call bndry_exchangeV(hybrid,edgeAdv_p1)
  endif

  do ie=nets,nete

     if ( DSSopt == DSSeta) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
     if ( DSSopt == DSSomega) DSSvar => elem(ie)%derived%omega_p(:,:,:)
     if ( DSSopt == DSSdiv_vdp_ave) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

     if(DSSopt==DSSno_var)then
	call edgeDGVunpack(edgeAdv,qedges,nlev*qsize,0,elem(ie)%desc)
     else
	call edgeDGVunpack(edgeAdv_p1,qedges,nlev*qsize,0,elem(ie)%desc)
	call edgeVunpack(edgeAdv_p1,DSSvar(:,:,1:nlev),nlev,qsize*nlev,elem(ie)%desc)
	do k=1,nlev
	  DSSvar(:,:,k)=DSSvar(:,:,k)*elem(ie)%rspheremp(:,:)
	enddo
     endif

     ! compute flux and advection term
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
     do k=1,nlev
        dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - &
             rhs_multiplier*dt*elem(ie)%derived%divdp_proj(:,:,k)
        Vstar(:,:,1,k) = elem(ie)%derived%vn0(:,:,1,k)/dp(:,:,k)
        Vstar(:,:,2,k) = elem(ie)%derived%vn0(:,:,2,k)/dp(:,:,k)
     enddo

#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k,vtemp,divdp,pshat,j,i)
#endif
     do q=1,qsize
        do k=1,nlev
           vtemp(:,:,1)=elem(ie)%state%Qdp(:,:,k,q,n0_qdp)*Vstar(:,:,1,k)
           vtemp(:,:,2)=elem(ie)%state%Qdp(:,:,k,q,n0_qdp)*Vstar(:,:,2,k)

           divdp = divergence_sphere_wk(vtemp,deriv,elem(ie)) + &
                edge_flux_u_cg( Vstar(:,:,:,k), elem(ie)%state%Qdp(:,:,k,q,n0_qdp),qedges(:,:,k,q),&
                deriv, elem(ie), u_is_contra=.false.)

           ! advance in time. GLL quadrature, cardinal function basis, under-integrated.
           ! local mass matrix is diagonal, with entries elem(ie)%spheremp(),
           ! so we divide through by elem(ie)%spheremp().
           elem(ie)%state%Qdp(:,:,k,q,np1_qdp)=elem(ie)%state%Qdp(:,:,k,q,n0_qdp) - dt*divdp/elem(ie)%spheremp

           if (npdg<np) then
              ! modal timestep, with exact integration.  using prognostic variable: p*metdet
              ! local mass matrix is diagonal assuming npdg<np so that GLL quadrature is exact)
              ! (note: GLL/modal conversion comutes with time-stepping)

              ! compute modal coefficients of p*metdet
              ! (spherical inner-product of Legendre polynomial and p)
              pshat = gll_to_dgmodal(elem(ie)%state%Qdp(:,:,k,q,np1_qdp)*elem(ie)%metdet(:,:),deriv)

              ! modal based limiter goes here
              ! apply a little dissipation to last mode:
              do j=1,npdg
              do i=1,npdg
                 !if ( (i-1)+(j-1) == 4) pshat(i,j)=pshat(i,j)*.75
                 !if ( (i-1)+(j-1) == 3) pshat(i,j)=pshat(i,j)*.90
                 if ( i==npdg) pshat(i,j)=pshat(i,j)*.90
                 if ( j==npdg) pshat(i,j)=pshat(i,j)*.90
              enddo
              enddo


              ! evalute modal expanion of p*metdet on GLL points
              divdp=dgmodal_to_gll(pshat,deriv)

              ! convert from p*metdet back to p:
              elem(ie)%state%Qdp(:,:,k,q,np1_qdp)=divdp/elem(ie)%metdet(:,:)
           endif
        enddo
        if(limiter_option == 4)then
           ! reuse CG limiter, which wants Qdp*spheremp:
           do k=1,nlev
              elem(ie)%state%Qdp(:,:,k,q,np1_qdp)=elem(ie)%state%Qdp(:,:,k,q,np1_qdp)*elem(ie)%spheremp(:,:)
           enddo
           call limiter2d_zero(elem(ie)%state%Qdp(:,:,:,q,np1_qdp))
           do k=1,nlev
              elem(ie)%state%Qdp(:,:,k,q,np1_qdp)=elem(ie)%state%Qdp(:,:,k,q,np1_qdp)/elem(ie)%spheremp(:,:)
           enddo
        endif
     end do
  end do
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
  call t_stopf('euler_step_dg')

  end subroutine euler_step_dg

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter_optim_iter_full(ptens,sphweights,minp,maxp,dpmass)
    !THIS IS A NEW VERSION OF LIM8, POTENTIALLY FASTER BECAUSE INCORPORATES KNOWLEDGE FROM
    !PREVIOUS ITERATIONS

    !The idea here is the following: We need to find a grid field which is closest
    !to the initial field (in terms of weighted sum), but satisfies the min/max constraints.
    !So, first we find values which do not satisfy constraints and bring these values
    !to a closest constraint. This way we introduce some mass change (addmass),
    !so, we redistribute addmass in the way that l2 error is smallest.
    !This redistribution might violate constraints thus, we do a few iterations.
    use kinds         , only : real_kind
    use dimensions_mod, only : np, np, nlev
    real (kind=real_kind), dimension(np*np,nlev), intent(inout)            :: ptens
    real (kind=real_kind), dimension(np*np     ), intent(in   )            :: sphweights
    real (kind=real_kind), dimension(      nlev), intent(inout)            :: minp
    real (kind=real_kind), dimension(      nlev), intent(inout)            :: maxp
    real (kind=real_kind), dimension(np*np,nlev), intent(in   ), optional  :: dpmass

    real (kind=real_kind), dimension(np*np,nlev) :: weights
    integer  k1, k, i, j, iter, i1, i2
    integer :: whois_neg(np*np), whois_pos(np*np), neg_counter, pos_counter
    real (kind=real_kind) :: addmass, weightssum, mass
    real (kind=real_kind) :: x(np*np),c(np*np)
    real (kind=real_kind) :: al_neg(np*np), al_pos(np*np), howmuch
    real (kind=real_kind) :: tol_limiter = 1e-15
    integer, parameter :: maxiter = 5

    do k = 1 , nlev
      weights(:,k) = sphweights(:) * dpmass(:,k)
      ptens(:,k) = ptens(:,k) / dpmass(:,k)
    enddo

    do k = 1 , nlev
      c = weights(:,k)
      x = ptens(:,k)

      mass = sum(c*x)

      ! relax constraints to ensure limiter has a solution:
      ! This is only needed if runnign with the SSP CFL>1 or
      ! due to roundoff errors
      if( (mass / sum(c)) < minp(k) ) then
        minp(k) = mass / sum(c)
      endif
      if( (mass / sum(c)) > maxp(k) ) then
        maxp(k) = mass / sum(c)
      endif

      addmass = 0.0d0
      pos_counter = 0
      neg_counter = 0

      ! apply constraints, compute change in mass caused by constraints
      do k1 = 1 , np*np
        if ( ( x(k1) >= maxp(k) ) ) then
          addmass = addmass + ( x(k1) - maxp(k) ) * c(k1)
          x(k1) = maxp(k)
          whois_pos(k1) = -1
        else
          pos_counter = pos_counter+1
          whois_pos(pos_counter) = k1
        endif
        if ( ( x(k1) <= minp(k) ) ) then
          addmass = addmass - ( minp(k) - x(k1) ) * c(k1)
          x(k1) = minp(k)
          whois_neg(k1) = -1
        else
          neg_counter = neg_counter+1
          whois_neg(neg_counter) = k1
        endif
      enddo

      ! iterate to find field that satifies constraints and is l2-norm closest to original
      weightssum = 0.0d0
      if ( addmass > 0 ) then
        do i2 = 1 , maxIter
          weightssum = 0.0
          do k1 = 1 , pos_counter
            i1 = whois_pos(k1)
            weightssum = weightssum + c(i1)
            al_pos(i1) = maxp(k) - x(i1)
          enddo

          if( ( pos_counter > 0 ) .and. ( addmass > tol_limiter * abs(mass) ) ) then
            do k1 = 1 , pos_counter
              i1 = whois_pos(k1)
              howmuch = addmass / weightssum
              if ( howmuch > al_pos(i1) ) then
                howmuch = al_pos(i1)
                whois_pos(k1) = -1
              endif
              addmass = addmass - howmuch * c(i1)
              weightssum = weightssum - c(i1)
              x(i1) = x(i1) + howmuch
            enddo
            !now sort whois_pos and get a new number for pos_counter
            !here neg_counter and whois_neg serve as temp vars
            neg_counter = pos_counter
            whois_neg = whois_pos
            whois_pos = -1
            pos_counter = 0
            do k1 = 1 , neg_counter
              if ( whois_neg(k1) .ne. -1 ) then
                pos_counter = pos_counter+1
                whois_pos(pos_counter) = whois_neg(k1)
              endif
            enddo
          else
            exit
          endif
        enddo
      else
         do i2 = 1 , maxIter
           weightssum = 0.0
           do k1 = 1 , neg_counter
             i1 = whois_neg(k1)
             weightssum = weightssum + c(i1)
             al_neg(i1) = x(i1) - minp(k)
           enddo

           if ( ( neg_counter > 0 ) .and. ( (-addmass) > tol_limiter * abs(mass) ) ) then
             do k1 = 1 , neg_counter
               i1 = whois_neg(k1)
               howmuch = -addmass / weightssum
               if ( howmuch > al_neg(i1) ) then
                 howmuch = al_neg(i1)
                 whois_neg(k1) = -1
               endif
               addmass = addmass + howmuch * c(i1)
               weightssum = weightssum - c(i1)
               x(i1) = x(i1) - howmuch
             enddo
             !now sort whois_pos and get a new number for pos_counter
             !here pos_counter and whois_pos serve as temp vars
             pos_counter = neg_counter
             whois_pos = whois_neg
             whois_neg = -1
             neg_counter = 0
             do k1 = 1 , pos_counter
               if ( whois_pos(k1) .ne. -1 ) then
                 neg_counter = neg_counter+1
                 whois_neg(neg_counter) = whois_pos(k1)
               endif
             enddo
           else
             exit
           endif
         enddo
      endif

      ptens(:,k) = x
    enddo

    do k = 1 , nlev
      ptens(:,k) = ptens(:,k) * dpmass(:,k)
    enddo
  end subroutine limiter_optim_iter_full


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter2d_minmax(Qdp,dp,spheremp,qmin,qmax)
!
! mass conserving limiter (2D only).  to be called just before DSS
!
! in pure 2D advection, the element mass will not be negative before DSS
! this routine will redistribute to remove negative values (conservative)
!
! call with Qdp and assocated dp
!
!
  implicit none
  real (kind=real_kind), intent(inout) :: Qdp(np,np,nlev)
  real (kind=real_kind), intent(in   ) :: spheremp(np,np)
  real (kind=real_kind), intent(in   ) ::  dp(np,np,nlev)
  real (kind=real_kind), intent(in   ) :: qmin(nlev),qmax(nlev)

  ! local
  integer i,j,k
  real (kind=real_kind) :: mass,mass_new,area,mass2
  real (kind=real_kind) :: Q(np,np,nlev)

  Q(:,:,:)=Qdp(:,:,:)/dp(:,:,:)  ! convert to concentration


#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,mass,area,mass2,mass_new,i,j)
#endif
  do k = 1 , nlev
    mass = sum( Qdp(:,:,k)*spheremp(:,:) )
    area = sum( dp(:,:,k)*spheremp(:,:) )

!   if (mass>0) print *,k,mass/area,qmin(k),qmax(k)
    ! max limiter
    if ( maxval(Q(:,:,k)) > qmax(k) ) then
      Q(:,:,k) = qmax(k) - Q(:,:,k)      ! some of these will be negative
      mass2 = area * qmax(k) - mass

      if (mass2 < 0) Q(:,:,k) = -Q(:,:,k)
      mass_new = 0
      do j = 1 , np
        do i = 1 , np
          if ( Q(i,j,k) < 0 ) then
            Q(i,j,k) = 0
          else
            mass_new = mass_new + Q(i,j,k) * dp(i,j,k) * spheremp(i,j)
          endif
        enddo
      enddo

      ! now scale the all positive values to restore mass
      if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass2) / mass_new
      if ( mass2    < 0 ) Q(:,:,k) = -Q(:,:,k)

      Q(:,:,k) = qmax(k) - Q(:,:,k)
    endif

    ! min limiter
    if ( minval(Q(:,:,k)) < qmin(k) ) then
      Q(:,:,k) = Q(:,:,k) - qmin(k)
      mass2 = mass - area * qmin(k)
      ! negative mass.  so reduce all postive values to zero
      ! then increase negative values as much as possible
      if ( mass2 < 0 ) Q(:,:,k) = -Q(:,:,k)
      mass_new = 0
      do j = 1 , np
        do i = 1 , np
          if ( Q(i,j,k) < 0 ) then
            Q(i,j,k) = 0
          else
            mass_new = mass_new + Q(i,j,k) * dp(i,j,k) * spheremp(i,j)
          endif
        enddo
      enddo

      ! now scale the all positive values to restore mass
      if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass2) / mass_new
      if ( mass2    < 0 ) Q(:,:,k) = -Q(:,:,k)

      Q(:,:,k) = Q(:,:,k) + qmin(k)
    endif

    Qdp(:,:,k) = Q(:,:,k) * dp(:,:,k)
  enddo
  end subroutine limiter2d_minmax

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter2d_zero(Q)
  ! mass conserving zero limiter (2D only).  to be called just before DSS
  !
  ! this routine is called inside a DSS loop, and so Q had already
  ! been multiplied by the mass matrix.  Thus dont include the mass
  ! matrix when computing the mass = integral of Q over the element
  !
  ! ps is only used when advecting Q instead of Qdp
  ! so ps should be at one timelevel behind Q
  implicit none
  real (kind=real_kind), intent(inout) :: Q(np,np,nlev)

  ! local
  real (kind=real_kind) :: dp(np,np)
  real (kind=real_kind) :: mass,mass_new,ml
  integer i,j,k

  do k = nlev , 1 , -1
    mass = 0
    do j = 1 , np
      do i = 1 , np
        !ml = Q(i,j,k)*dp(i,j)*spheremp(i,j)  ! see above
        ml = Q(i,j,k)
        mass = mass + ml
      enddo
    enddo

    ! negative mass.  so reduce all postive values to zero
    ! then increase negative values as much as possible
    if ( mass < 0 ) Q(:,:,k) = -Q(:,:,k)
    mass_new = 0
    do j = 1 , np
      do i = 1 , np
        if ( Q(i,j,k) < 0 ) then
          Q(i,j,k) = 0
        else
          ml = Q(i,j,k)
          mass_new = mass_new + ml
        endif
      enddo
    enddo

    ! now scale the all positive values to restore mass
    if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass) / mass_new
    if ( mass     < 0 ) Q(:,:,k) = -Q(:,:,k)
  enddo
  end subroutine limiter2d_zero

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine advance_hypervis_scalar( edgeAdv , elem , hvcoord , hybrid , deriv , nt , nt_qdp , nets , nete , dt2 )
  !  hyperviscsoity operator for foward-in-time scheme
  !  take one timestep of:
  !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
#if USE_CUDA_FORTRAN
  use cuda_mod       , only : advance_hypervis_scalar_cuda
#endif
  use kinds          , only : real_kind
  use dimensions_mod , only : np, nlev
  use hybrid_mod     , only : hybrid_t
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t
  use edge_mod       , only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod      , only : bndry_exchangev
  use perf_mod       , only : t_startf, t_stopf                          ! _EXTERNAL
  implicit none
  type (EdgeBuffer_t)  , intent(inout)         :: edgeAdv
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nt
  integer              , intent(in   )         :: nt_qdp
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  real (kind=real_kind), intent(in   )         :: dt2

  ! local
  real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: Qtens
  real (kind=real_kind), dimension(np,np,nlev                ) :: dp
  real (kind=real_kind), dimension(      nlev,qsize,nets:nete) :: min_neigh
  real (kind=real_kind), dimension(      nlev,qsize,nets:nete) :: max_neigh
  integer :: k , kptr , i , j , ie , ic , q

! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
!       data is incorrect (offset by a few numbers actually)
!       removed for now.
!  real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  real (kind=real_kind), dimension(np,np) :: lap_p
  real (kind=real_kind) :: v1,v2,dt,dp0
  integer :: density_scaling = 0
  if ( nu_q           == 0 ) return
  if ( hypervis_order /= 2 ) return
#if USE_CUDA_FORTRAN
  call advance_hypervis_scalar_cuda( edgeAdv , elem , hvcoord , hybrid , deriv , nt , nt_qdp , nets , nete , dt2 )
  return
#endif
!   call t_barrierf('sync_advance_hypervis_scalar', hybrid%par%comm)
  call t_startf('advance_hypervis_scalar')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dt = dt2 / hypervis_subcycle_q

  do ic = 1 , hypervis_subcycle_q
    do ie = nets , nete
      ! Qtens = Q/dp   (apply hyperviscsoity to dp0 * Q, not Qdp)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,dp0,q)
#endif
      do k = 1 , nlev
         ! various options:
         !   1)  biharmonic( Qdp )
         !   2)  dp0 * biharmonic( Qdp/dp )
         !   3)  dpave * biharmonic(Q/dp)
         ! For trace mass / mass consistenciy, we use #2 when nu_p=0
         ! and #e when nu_p>0, where dpave is the mean mass flux from the nu_p
         ! contribution from dynamics.
         dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) ) * hvcoord%ps0 + &
              ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * hvcoord%ps0
         dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - dt2*elem(ie)%derived%divdp_proj(:,:,k)
         if (nu_p>0) then
            do q = 1 , qsize
               Qtens(:,:,k,q,ie) = elem(ie)%derived%dpdiss_ave(:,:,k)*&
                    elem(ie)%state%Qdp(:,:,k,q,nt_qdp) / dp(:,:,k)
            enddo
         else
            do q = 1 , qsize
               Qtens(:,:,k,q,ie) = dp0*elem(ie)%state%Qdp(:,:,k,q,nt_qdp) / dp(:,:,k)
            enddo
         endif
      enddo
   enddo

    ! compute biharmonic operator. Qtens = input and output
    call biharmonic_wk_scalar( elem , Qtens , deriv , edgeAdv , hybrid , nets , nete )
    do ie = nets , nete
      !spheremp     => elem(ie)%spheremp
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k,j,i,dp0)
#endif
      do q = 1 , qsize
        do k = 1 , nlev
          dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) ) * hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * hvcoord%ps0
          do j = 1 , np
            do i = 1 , np
              ! advection Qdp.  For mass advection consistency:
              ! DIFF( Qdp) ~   dp0 DIFF (Q)  =  dp0 DIFF ( Qdp/dp )
              elem(ie)%state%Qdp(i,j,k,q,nt_qdp) = elem(ie)%state%Qdp(i,j,k,q,nt_qdp) * elem(ie)%spheremp(i,j) &
                                                   - dt * nu_q * Qtens(i,j,k,q,ie)
            enddo
          enddo
        enddo

        ! smooth some of the negativities introduced by diffusion:
        call limiter2d_zero( elem(ie)%state%Qdp(:,:,:,q,nt_qdp))
      enddo
      call edgeVpack  ( edgeAdv , elem(ie)%state%Qdp(:,:,:,:,nt_qdp) , qsize*nlev , 0 , elem(ie)%desc )
    enddo

    call bndry_exchangeV( hybrid , edgeAdv )

    do ie = nets , nete
      call edgeVunpack( edgeAdv , elem(ie)%state%Qdp(:,:,:,:,nt_qdp) , qsize*nlev , 0 , elem(ie)%desc )
      !rspheremp     => elem(ie)%rspheremp
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k)
#endif
      do q = 1 , qsize
        ! apply inverse mass matrix
        do k = 1 , nlev
          elem(ie)%state%Qdp(:,:,k,q,nt_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,nt_qdp)
        enddo
      enddo
    enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
  enddo
  call t_stopf('advance_hypervis_scalar')
  end subroutine advance_hypervis_scalar





  subroutine vertical_remap(hybrid,elem,fvm,hvcoord,dt,np1,np1_qdp,nets,nete)
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
  use kinds, only : real_kind
  use hybvcoord_mod, only : hvcoord_t
  use vertremap_mod, only : remap1, remap1_nofilter, remap_q_ppm ! _EXTERNAL (actually INTERNAL)
  use control_mod, only :  rsplit, tracer_transport_type, TRACER_GRIDTYPE_FVM
  use parallel_mod, only : abortmp
  use hybrid_mod     , only : hybrid_t
  use derivative_mod, only : interpolate_gll2fvm_points
#if defined(_SPELT)
  use spelt_mod, only: spelt_struct
#else
  use fvm_control_volume_mod, only : fvm_struct
#endif
#if USE_CUDA_FORTRAN
  use cuda_mod, only: vertical_remap_cuda
#endif

  type (hybrid_t), intent(in) :: hybrid  ! distributed parallel structure (shared)
#if defined(_SPELT)
  type(spelt_struct), intent(inout) :: fvm(:)
  real (kind=real_kind) :: cdp(1:nep,1:nep,nlev,ntrac-1)
  real (kind=real_kind)  :: psc(nep,nep), dpc(nep,nep,nlev),dpc_star(nep,nep,nlev)
#else
  type(fvm_struct), intent(inout) :: fvm(:)
  real (kind=real_kind) :: cdp(1:nc,1:nc,nlev,ntrac)
  real (kind=real_kind)  :: psc(nc,nc), dpc(nc,nc,nlev),dpc_star(nc,nc,nlev)
#endif

  !    type (hybrid_t), intent(in)       :: hybrid  ! distributed parallel structure (shared)
  type (element_t), intent(inout)   :: elem(:)
  type (hvcoord_t)                  :: hvcoord
  real (kind=real_kind)             :: dt

  integer :: ie,i,j,k,np1,nets,nete,np1_qdp
  real (kind=real_kind), dimension(np,np,nlev)  :: dp,dp_star
  real (kind=real_kind), dimension(np,np,nlev,2)  :: ttmp

#if USE_CUDA_FORTRAN
  call vertical_remap_cuda(elem,fvm,hvcoord,dt,np1,np1_qdp,nets,nete)
  return
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

  do ie=nets,nete
!     ! SET VERTICAL VELOCITY TO ZERO FOR DEBUGGING
!     elem(ie)%derived%eta_dot_dpdn(:,:,:)=0
     if (rsplit==0) then
        ! compute dp_star from eta_dot_dpdn():
        do k=1,nlev
           dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
           dp_star(:,:,k) = dp(:,:,k) + dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) -&
                elem(ie)%derived%eta_dot_dpdn(:,:,k))
        enddo
        if (minval(dp_star)<0) call abortmp('negative layer thickness.  timestep or remap time too large')
     else
        !  REMAP u,v,T from levels in dp3d() to REF levels
        !
        ! update final ps_v
        elem(ie)%state%ps_v(:,:,np1) = hvcoord%hyai(1)*hvcoord%ps0 + &
             sum(elem(ie)%state%dp3d(:,:,:,np1),3)
        do k=1,nlev
           dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
           dp_star(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)
        enddo
        if (minval(dp_star)<0) call abortmp('negative layer thickness.  timestep or remap time too large')

        ! remap the dynamics:
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
        call remap1(ttmp,np,1,dp_star,dp)
        elem(ie)%state%t(:,:,:,np1)=ttmp(:,:,:,1)/dp

        ttmp(:,:,:,1)=elem(ie)%state%v(:,:,1,:,np1)*dp_star
        ttmp(:,:,:,2)=elem(ie)%state%v(:,:,2,:,np1)*dp_star
        call remap1(ttmp,np,2,dp_star,dp)
!        call remap1_nofilter(ttmp,np,2,dp_star,dp)
        elem(ie)%state%v(:,:,1,:,np1)=ttmp(:,:,:,1)/dp
        elem(ie)%state%v(:,:,2,:,np1)=ttmp(:,:,:,2)/dp
#ifdef REMAP_TE
        ! back out T from TE
        elem(ie)%state%t(:,:,:,np1) = &
             ( elem(ie)%state%t(:,:,:,np1) - ( (elem(ie)%state%v(:,:,1,:,np1)**2 + &
             elem(ie)%state%v(:,:,2,:,np1)**2)/2))/cp

#endif
     endif

     ! remap the tracers from lagrangian levels (dp_star)  to REF levels dp
     if (qsize>0) then
        call remap1(elem(ie)%state%Qdp(:,:,:,:,np1_qdp),np,qsize,dp_star,dp)
     endif


     if (ntrac>0) then
#if defined(_SPELT)
        do i=1,nep
          do j=1,nep
            ! 1. compute surface pressure, 'ps_c', from SPELT air density
            psc(i,j)=sum(fvm(ie)%c(i,j,:,1,np1)) +  hvcoord%hyai(1)*hvcoord%ps0
            ! 2. compute dp_np1 using CSLAM air density and eta coordinate formula
            ! get the dp now on the eta coordinates (reference level)
            do k=1,nlev
              dpc(i,j,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                              (hvcoord%hybi(k+1) - hvcoord%hybi(k))*psc(i,j)
              cdp(i,j,k,1:(ntrac-1))=fvm(ie)%c(i,j,k,2:ntrac,np1)*fvm(ie)%c(i,j,k,1,np1)
              dpc_star(i,j,k)=fvm(ie)%c(i,j,k,1,np1)
            end do
          end do
        end do
        call remap1(cdp,nep,ntrac-1,dpc_star,dpc)
        do i=1,nep
          do j=1,nep
            do k=1,nlev
              fvm(ie)%c(i,j,k,1,np1)=dpc(i,j,k)
              fvm(ie)%c(i,j,k,2:ntrac,np1)=cdp(i,j,k,1:(ntrac-1))/dpc(i,j,k)
            end do
          end do
        end do
!         call remap_velocityCspelt(np1,dt,elem,fvm,hvcoord,ie)
#else
        ! create local variable  cdp(1:nc,1:nc,nlev,ntrac)
        ! cdp(:,:,:,n) = fvm%c(:,:,:,n,np1)*fvm%dp_fvm(:,:,:,np1)
        ! dp(:,:,:) = reference level thicknesses

        ! call remap1(cdp,nc,ntrac-1,fvm%c(:,:,:,1,np1),dp)

        ! convert back to mass:
        ! fvm%dp_fvm(:,:,:,np1) = dp(:,:,:) ??XXgoldyXX??
        ! fvm%c(:,:,:,n,np1) = fvm%c(:,:,:,n,np1)/dp(:,:,:)
        if (tracer_transport_type == TRACER_GRIDTYPE_FVM) then
          do i=1,nc
            do j=1,nc
              ! 1. compute surface pressure, 'ps_c', from FVMair density

! phl dp_cslam
              psc(i,j)=sum(fvm(ie)%dp_fvm(i,j,:,np1)) +  hvcoord%hyai(1)*hvcoord%ps0
              ! 2. compute dp_np1 using FVM air density and eta coordinate formula
              ! get the dp now on the eta coordinates (reference level)
              do k=1,nlev
                dpc(i,j,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                             (hvcoord%hybi(k+1) - hvcoord%hybi(k))*psc(i,j)
                cdp(i,j,k,1:ntrac)=fvm(ie)%c(i,j,k,1:ntrac,np1)*fvm(ie)%dp_fvm(i,j,k,np1)
              end do
            end do
          end do
          dpc_star=fvm(ie)%dp_fvm(1:nc,1:nc,:,np1)
          call remap1(cdp,nc,ntrac,dpc_star,dpc)
! phl ntrac-1 -> ntrac       call remap1(cdp,nc,ntrac-1,dpc_star,dpc)
          do k=1,nlev
            do j=1,nc
              do i=1,nc
                fvm(ie)%dp_fvm(i,j,k,np1)=dpc(i,j,k) !!XXgoldyXX??
                fvm(ie)%c(i,j,k,:,np1)=cdp(i,j,k,:)/dpc(i,j,k)
              end do
            end do
          end do
          ! Recompute dp_fvm
          do k = 1, nlev
            fvm(ie)%dp_fvm(1:nc,1:nc,k,np1)=interpolate_gll2fvm_points(dp(:,:,k),deriv(hybrid%ithr))
          end do
        end if
!         call remap_velocityC(np1,dt,elem,fvm,hvcoord,ie)
#endif
     endif

  enddo
  call t_stopf('vertical_remap')
  end subroutine vertical_remap




end module prim_advection_mod
