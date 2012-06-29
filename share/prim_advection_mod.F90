#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef CAM
#else
! sometimes used for debugging REMAP
#undef ZEROVERT
#undef ZEROHORZ
#endif

#if 0
SUBROUTINES:
   prim_advec_tracers()
      Full Euler + hypervis
      oscillatory, non-conservative QNEG column fixer
   prim_advec_tracers_remap_rk2()
      SEM 2D RK2 + monotone remap + hyper viscosity
      SEM 2D RK2 can use sign-preserving or monotone reconstruction

Notes on Lagrange+REMAP advection
dynamics looks like (i.e. for qsplit=3)

    dp(t+1)-dp(t)   = -dt div(Udp1) - dt d(eta_dot_dpdn1)  + dt D(dpdiss1)
    dp(t+2)-dp(t+1) = -dt div(Udp2) - dt d(eta_dot_dpdn2)  + dt D(dpdiss2)
    dp(t+3)-dp(t+2) = -dt div(Udp3) - dt d(eta_dot_dpdn3)  + dt D(dpdiss3)
    ---------------
    dp(t+3)-dp(t)        = -3dt div(Udp_sum/3) - 3dt d(eta_dot_dpdn_sum/3)  + 3dt D(dpdiss_sum/3)
    dpstart(t+3) - dp(t) = -3dt div(Udp_sum/3)  + 3dt D(dpdiss_sum/3)
OR:
    dp_star(t+3) = dp(t+1) + 3dt d( eta_dot_dpdn_ave(t) ) 


For RK2 advection of Q:  (example of 2 stage RK for tracers):   dtq = qsplit*dt
For consistency, if Q=1
  dp1  = dp(t)- dtq div[ U1 dp(t)]     
  dp2  = dp1  - dtq div[ U2 dp1  ]  + 2*dtq D( dpdiss_ave )   
  dp*  = (dp(t) + dp2 )/2
       =  dp(t) - dtq  div[ U1 dp(t) + U2 dp1 ]   + dtq D( dpdiss_ave )

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

  use kinds, only          : real_kind,int_kind
  use dimensions_mod, only : np,nlev,qsize,nlevp,npsq,ntrac
  use hybvcoord_mod, only  : hvcoord_t
  use element_mod, only    : element_t
  use fvm_control_volume_mod, only    : fvm_struct
  use perf_mod, only	   : t_startf, t_stopf  ! _EXTERNAL
  use parallel_mod, only   : abortmp  	
  contains

  subroutine remap_velocityUV(np1,dt,elem,hvcoord,nets,nete)

    implicit none
    real (kind=real_kind),  intent(in)          :: dt
    type (element_t),    intent(inout), target  :: elem(:)
    type (hvcoord_t),    intent(in)             :: hvcoord
    
    integer :: nets,nete,np1
     
    ! ========================
    ! Local Variables
    ! ========================

    real(kind=real_kind), dimension(nlev) :: dp,dp_star,Ustar,Vstar
    real(kind=real_kind), dimension(nlev) :: Uold,Vold,Unew,Vnew
    real(kind=real_kind), dimension(nlevp):: z1cU,z2cU,z1cV,z2cV
    real(kind=real_kind), dimension(nlev+1)    :: rhsU,lower_diagU,diagU,upper_diagU,q_diagU,zgamU, & 
                             rhsV,lower_diagV,diagV,upper_diagV,q_diagV,zgamV
    real(kind=real_kind), dimension(nlev)    :: hU,rho_barU,za0U,za1U,za2U,zhdpU, & 
                             hV,rho_barV,za0V,za1V,za2V,zhdpV
    real(kind=real_kind)            :: tmp_calU,zaccintegerbU,zacctopU,zaccbotU, & 
                             tmp_calV,zaccintegerbV,zacctopV,zaccbotV
    
    integer(kind=int_kind) :: zkrU(nlev+1),zkrV(nlev+1),ie,i,j,k,jl,jk,ilevU,itopU,ibotU,ilevV,itopV,ibotV,jsubz,ij
    
    call t_startf('remap_velocityUV')
    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(ij,i,j,k,dp,dp_star,Ustar,Vstar,z1cU,z2cU,z1cV,z2cV,Uold,Vold) &
!$omp    private(zkrU,ilevU,zgamU,zhdpU,zkrV,ilevV,zgamV,zhdpV,jl,jk,hU,rhsU,hV,rhsV) &
!$omp    private(lower_diagU,lower_diagV,diagU,diagV,upper_diagU,upper_diagV,q_diagU) &
!$omp    private(q_diagV,tmp_calU,tmp_calV,za0U,za1U,za2U,za0V,za1V,za2V,zaccintegerbU) &
!$omp    private(itopU,zacctopU,zaccintegerbV,itopV,zacctopV,ibotU,jsubz,zaccbotU,ibotV) &
!$omp    private(zaccbotV,rho_barU,rho_barV)
#endif
       do ij = 1, npsq
          j = (ij-1)/np + 1
          i = ij - (j-1)*np
          do k=1,nlev
            dp(k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,np1)
            dp_star(k) = dp(k) + dt*(elem(ie)%derived%eta_dot_dpdn(i,j,k+1) - elem(ie)%derived%eta_dot_dpdn(i,j,k) ) 
            Ustar(k) = elem(ie)%state%v(i,j,1,k,np1)*dp(k)
            Vstar(k) = elem(ie)%state%v(i,j,2,k,np1)*dp(k)
          enddo
          

          z1cU(1)=0
          z2cU(1)=0
          z1cV(1)=0
          z2cV(1)=0
          do k=1,nlev
            Uold(k) = Ustar(k)
            z1cU(k+1) = z1cU(k)+dp(k)
            z2cU(k+1) = z2cU(k)+dp_star(k)
            Vold(k) = Vstar(k)
            z1cV(k+1) = z1cV(k)+dp(k)
            z2cV(k+1) = z2cV(k)+dp_star(k)
          enddo
		  
          if (ABS(z2cU(nlev+1)-z1cU(nlev+1)).GE.0.000001) then
             write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
             write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
             write(6,*) 'DATA FOR MODEL LEVELS'
             write(6,*) 'PLEVMODEL=',z2cU(nlev+1)
             write(6,*) 'PLEV     =',z1cU(nlev+1)
             write(6,*) 'DIFF     =',z2cU(nlev+1)-z1cU(nlev+1)
             ! call ABORT
          endif
          
          if (ABS(z2cV(nlev+1)-z1cV(nlev+1)).GE.0.000001) then
             write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
             write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
             write(6,*) 'DATA FOR MODEL LEVELS'
             write(6,*) 'PLEVMODEL=',z2cV(nlev+1)
             write(6,*) 'PLEV     =',z1cV(nlev+1)
             write(6,*) 'DIFF     =',z2cV(nlev+1)-z1cV(nlev+1)
             ! call ABORT
          endif
          
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!  calculate quadratic splies !!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          zkrU  = 99
          ilevU = 2
          zkrU(1)       = 1
          zgamU(1)      = 0.0
          zkrU(nlev+1)  = nlev
          zgamU(nlev+1) = 1.0
          zhdpU(1) = z1cU(2)-z1cU(1)
          
          zkrV  = 99
          ilevV = 2
          zkrV(1)       = 1
          zgamV(1)      = 0.0
          zkrV(nlev+1)  = nlev
          zgamV(nlev+1) = 1.0
          zhdpV(1) = z1cV(2)-z1cV(1)
          do jl = 2,nlev
            zhdpU(jl) = z1cU(jl+1)-z1cU(jl)
            jkloopU: do jk = ilevU,nlev+1
              if (z1cU(jk).ge.z2cU(jl)) then
                ilevU      = jk
                zkrU(jl)   = jk-1
                zgamU(jl)   = (z2cU(jl)-z1cU(jk-1))/(z1cU(jk)-z1cU(jk-1))
                exit jkloopU
              endif
            enddo jkloopU
            zhdpV(jl) = z1cV(jl+1)-z1cV(jl)
            jkloopV: do jk = ilevV,nlev+1
              if (z1cV(jk).ge.z2cV(jl)) then
                ilevV      = jk
                zkrV(jl)   = jk-1
                zgamV(jl)   = (z2cV(jl)-z1cV(jk-1))/(z1cV(jk)-z1cV(jk-1))
                exit jkloopV
              endif
            enddo jkloopV
          enddo 

          hU = 1/zhdpU
          rho_barU = Uold * hU
          rhsU = 0
          lower_diagU = 0
          diagU = 0
          upper_diagU = 0
          
          hV = 1/zhdpV 
          rho_barV = Vold * hV          
          rhsV = 0
          lower_diagV = 0
          diagV = 0
          upper_diagV = 0

          rhsU(1)=3*rho_barU(1)
          rhsU(2:nlev) = 3*(rho_barU(2:nlev)*hU(2:nlev) + rho_barU(1:nlev-1)*hU(1:nlev-1)) 
          rhsU(nlev+1)=3*rho_barU(nlev)
          
          rhsV(1)=3*rho_barV(1)
          rhsV(2:nlev) = 3*(rho_barV(2:nlev)*hV(2:nlev) + rho_barV(1:nlev-1)*hV(1:nlev-1)) 
          rhsV(nlev+1)=3*rho_barV(nlev)

          lower_diagU(1)=1
          lower_diagU(2:nlev) = hU(1:nlev-1)
          lower_diagU(nlev+1)=1
          
          lower_diagV(1)=1
          lower_diagV(2:nlev) = hV(1:nlev-1)
          lower_diagV(nlev+1)=1

          diagU(1)=2
          diagU(2:nlev) = 2*(hU(2:nlev) + hU(1:nlev-1))
          diagU(nlev+1)=2
          
          diagV(1)=2
          diagV(2:nlev) = 2*(hV(2:nlev) + hV(1:nlev-1))
          diagV(nlev+1)=2

          upper_diagU(1)=1
          upper_diagU(2:nlev) = hU(2:nlev)
          upper_diagU(nlev+1)=0
          
          upper_diagV(1)=1
          upper_diagV(2:nlev) = hV(2:nlev)
          upper_diagV(nlev+1)=0

          q_diagU(1)=-upper_diagU(1)/diagU(1)
          rhsU(1)= rhsU(1)/diagU(1)
          
          q_diagV(1)=-upper_diagV(1)/diagV(1)
          rhsV(1)= rhsV(1)/diagV(1)
          do jl=2,nlev+1
            tmp_calU    =  1/(diagU(jl)+lower_diagU(jl)*q_diagU(jl-1))
            q_diagU(jl) = -upper_diagU(jl)*tmp_calU
            rhsU(jl) =  (rhsU(jl)-lower_diagU(jl)*rhsU(jl-1))*tmp_calU
            
            tmp_calV    =  1/(diagV(jl)+lower_diagV(jl)*q_diagV(jl-1))
            q_diagV(jl) = -upper_diagV(jl)*tmp_calV
            rhsV(jl) =  (rhsV(jl)-lower_diagV(jl)*rhsV(jl-1))*tmp_calV
          enddo
          do jl=nlev,1,-1
            rhsU(jl)=rhsU(jl)+q_diagU(jl)*rhsU(jl+1)
            rhsV(jl)=rhsV(jl)+q_diagV(jl)*rhsV(jl+1)
          enddo        

          za0U = rhsU(1:nlev)
          za1U = -4*rhsU(1:nlev) - 2*rhsU(2:nlev+1) + 6*rho_barU
          za2U = +3*rhsU(1:nlev) + 3*rhsU(2:nlev+1) - 6*rho_barU
          
          za0V = rhsV(1:nlev)
          za1V = -4*rhsV(1:nlev) - 2*rhsV(2:nlev+1) + 6*rho_barV
          za2V = +3*rhsV(1:nlev) + 3*rhsV(2:nlev+1) - 6*rho_barV
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! start iteration from top to bottom of atmosphere !! 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 
          zaccintegerbU = 0
          itopU = 1
          zacctopU = 0.0
          
          zaccintegerbV = 0
          itopV = 1
          zacctopV = 0.0


          do jk = 1,nlev
            ibotU = zkrU(jk+1)
            if (zgamU(jk+1)>1d0) then
               WRITE(*,*) 'r not in [0:1]', zgamU(jk+1)
            endif
            do jsubz=itopU,ibotU-1,1
              zaccintegerbU = zaccintegerbU + Uold(jsubz)
            enddo
            zaccbotU = zaccintegerbU + (za0U(ibotU)*zgamU(jk+1)+(za1U(ibotU)/2)*(zgamU(jk+1)**2)+(za2U(ibotU)/3)*(zgamU(jk+1)**3))*zhdpU(ibotU)
            elem(ie)%derived%vstar(i,j,1,jk) = (zaccbotU-zacctopU)/dp_star(jk)
            zacctopU        = zaccbotU
            itopU           = ibotU
            
            ibotV = zkrV(jk+1)
            if (zgamV(jk+1)>1d0) then
               WRITE(*,*) 'r not in [0:1]', zgamV(jk+1)
            endif
            do jsubz=itopV,ibotV-1,1
              zaccintegerbV = zaccintegerbV + Vold(jsubz)
            enddo
            zaccbotV = zaccintegerbV + (za0V(ibotV)*zgamV(jk+1)+(za1V(ibotV)/2)*(zgamV(jk+1)**2)+(za2V(ibotV)/3)*(zgamV(jk+1)**3))*zhdpV(ibotV)
            elem(ie)%derived%vstar(i,j,2,jk) = (zaccbotV-zacctopV)/dp_star(jk)
            zacctopV        = zaccbotV
            itopV           = ibotV
          enddo
       enddo
    enddo
    call t_stopf('remap_velocityUV')
  end subroutine remap_velocityUV
  
  
  
  subroutine remap_velocityQ(n0,np1,dt,elem,hvcoord,nets,nete,compute_diagnostics,rkstage)
  
    use physical_constants, only : cp, cpwater_vapor
	
    implicit none
    real (kind=real_kind),  intent(in)        :: dt
    type (element_t),    intent(inout), target  :: elem(:)
    type (hvcoord_t),    intent(in)        :: hvcoord
    logical,        intent(in)              :: compute_diagnostics
    
    integer :: nets,nete,n0,np1,rkstage
    
    ! ========================
    ! Local Variables
    ! ========================

    real (kind=real_kind), dimension(nlev+1)    :: rhs,lower_diag,diag,upper_diag,q_diag,zgam,z1c,z2c,zv
    real (kind=real_kind), dimension(nlev)      :: h,Qcol,dy,za0,za1,za2,zarg,zhdp
    real (kind=real_kind)  :: dp_star,dp_np1,f_xm,level1,level2,level3,level4,level5, &
								peaks_min,peaks_max,Q_vadv,tmp_cal,xm,xm_d,zv1,zv2, &
								zero = 0,one = 1,tiny = 1e-12,qmax = 1d50
    
    integer(kind=int_kind) :: zkr(nlev+1),filter_code(nlev),peaks,im1,im2,im3,ip1,ip2, & 
									lt1,lt2,lt3,t0,t1,t2,t3,t4,tm,tp,ie,i,ilev,j,jk,k,q
    
    call t_startf('remap_velocityQ')

    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,i,j,z1c,z2c,zv,k,dp_np1,dp_star,Qcol,zkr,ilev) &
!$omp    private(jk,zgam,zhdp,h,zarg,rhs,lower_diag,diag,upper_diag,q_diag,tmp_cal,filter_code) &
!$omp    private(dy,im1,im2,im3,ip1,t1,t2,t3,za0,za1,za2,xm_d,xm,f_xm,t4,tm,tp,peaks,peaks_min) &
!$omp    private(peaks_max,ip2,level1,level2,level3,level4,level5,lt1,lt2,lt3,zv1,zv2,Q_vadv)
#endif
      do q=1,qsize
        do i=1,np
          do j=1,np
             z1c(1)=0
             z2c(1)=0
             zv(1)=0
             do k=1,nlev
                dp_np1 = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                     (hvcoord%hybi(k+1) - hvcoord%hybi(k))*elem(ie)%state%ps_v(i,j,np1)
                dp_star = dp_np1 + dt*(elem(ie)%derived%eta_dot_dpdn(i,j,k+1) & 
                     -elem(ie)%derived%eta_dot_dpdn(i,j,k)) 
#ifdef ZEROVERT                        
                ! ignore the vertical motion
                dp_star = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                     (hvcoord%hybi(k+1) - hvcoord%hybi(k))*elem(ie)%state%ps_v(i,j,np1)
#endif

                z1c(k+1) = z1c(k)+dp_star
                z2c(k+1) = z2c(k)+dp_np1
#ifdef ZEROHORZ			  
                Qcol(k)=elem(ie)%state%Qdp(i,j,k,q,n0)
#else		
                Qcol(k)=(elem(ie)%state%Qdp(i,j,k,q,n0)+&
                     (rkstage-1)*elem(ie)%state%Qdp(i,j,k,q,np1))/rkstage
#endif			
                zv(k+1) = zv(k)+Qcol(k)
             enddo
             
             if (ABS(z2c(nlev+1)-z1c(nlev+1)).GE.0.000001) then
                write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
                write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
                write(6,*) 'DATA FOR MODEL LEVELS'
                write(6,*) 'PLEVMODEL=',z2c(nlev+1)
                write(6,*) 'PLEV     =',z1c(nlev+1)
                write(6,*) 'DIFF     =',z2c(nlev+1)-z1c(nlev+1)
                ! call ABORT
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
				endif
                zv2 = zv(zkr(k+1))+(za0(zkr(k+1))*zgam(k+1)+(za1(zkr(k+1))/2)*(zgam(k+1)**2)+ &
                     (za2(zkr(k+1))/3)*(zgam(k+1)**3))*zhdp(zkr(k+1))
                Q_vadv = (elem(ie)%state%Qdp(i,j,k,q,np1) - (zv2 - zv1)) / dt
                elem(ie)%state%Qdp(i,j,k,q,np1) = (zv2 - zv1)	
                zv1 = zv2
#ifdef ENERGY_DIAGNOSTICS
                if (compute_diagnostics .and. q==1) then
                   elem(ie)%accum%IEvert1_wet(i,j) = elem(ie)%accum%IEvert1_wet(i,j) + (Cpwater_vapor-Cp)*elem(ie)%state%T(i,j,k,n0)*Q_vadv
                endif
#endif	
             enddo
          enddo
       enddo
    enddo
 enddo
 
 call t_stopf('remap_velocityQ')
 end subroutine remap_velocityQ
 


!This uses the exact same model and reference grids and data as remap_velocityQ, but it interpolates
!using PPM instead of splines.
subroutine remap_velocityQ_ppm(n0,np1,dt,elem,hybrid,hvcoord,nets,nete,compute_diagnostics,rkstage)
  use hybrid_mod, only: hybrid_t
  use physical_constants, only : cp, cpwater_vapor
  use control_mod, only        : compute_mean_flux, prescribed_wind, vert_remap_q_alg
  implicit none
  real (kind=real_kind), intent(in)            :: dt
  type (element_t)     , intent(inout), target :: elem(:)
    type (hybrid_t),     intent(in):: hybrid
  type (hvcoord_t)     , intent(in)            :: hvcoord
  logical              , intent(in)            :: compute_diagnostics
  integer              , intent(in)            :: nets,nete,n0,np1,rkstage
  ! Local Variables
  integer, parameter :: gs = 2                              !Number of cells to place in the ghost region
  real(kind=real_kind), dimension(       nlev+2 ) :: pio    !Pressure at interfaces for old grid
  real(kind=real_kind), dimension(       nlev+1 ) :: pin    !Pressure at interfaces for new grid
  real(kind=real_kind), dimension(       nlev+1 ) :: masso  !Accumulate mass up to each interface
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: ao     !Tracer value on old grid
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpo    !change in pressure over a cell for old grid
  real(kind=real_kind), dimension(3,     nlev   ) :: coefs  !PPM coefficients within each cell
  real(kind=real_kind), dimension(       nlev   ) :: z1, z2
  real(kind=real_kind) :: ppmdx(10,0:nlev+1)  !grid spacings
  real(kind=real_kind) :: mymass, massn1, massn2, dpn
  integer :: i, j, k, q, ie, kk, kid(nlev)

  call t_startf('remap_velocityQ_ppm')
  do ie = nets , nete
    do j = 1 , np
      do i = 1 , np
        !Form old and new grids for remapping. It counts from the top of the atmosphere to the bottom. Starting with no mass at the top.
        pio(1) = 0.
        pin(1) = 0.
        do k = 1 , nlev
          dpn = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + (hvcoord%hybi(k+1) - hvcoord%hybi(k))*elem(ie)%state%ps_v(i,j,np1)
          if ( compute_mean_flux == 1 .and. prescribed_wind == 0 ) then
#ifdef ZEROVERT
            dpo(k) = dpn  ! ignore the vertical motion
#else
            dpo(k) = elem(ie)%derived%dp(i,j,k) - dt * elem(ie)%derived%divdp_proj(i,j,k)
#endif
          else
#ifdef ZEROVERT
            dpo(k) = dpn  ! ignore the vertical motion
#else
            dpo(k) = dpn + dt * ( elem(ie)%derived%eta_dot_dpdn(i,j,k+1) - elem(ie)%derived%eta_dot_dpdn(i,j,k) ) 
#endif
          endif
          pio(k+1) = pio(k) + dpo(k)
          pin(k+1) = pin(k) + dpn
        enddo
        pio(nlev+2) = pio(nlev+1) + 1.  !This is here to allow an entire block of k threads to run in the remapping phase.
                                        !It makes sure there's an old interface value below the domain that is larger.
        pin  (nlev+1) = pio  (nlev+1)   !The total mass in a column does not change. Therefore, the pressure of that mass cannot either.
        !Fill in the ghost regions with mirrored values. if vert_remap_q_alg is defined, this is of no consequence.
        do k = 1 , gs
          dpo(1   -k) = dpo(       k)
          dpo(nlev+k) = dpo(nlev+1-k)
        enddo

        !Compute remapping intervals once for all tracers. Find the old grid cell index in which the k-th new cell interface resides. Then integrate
        !from the bottom of that old cell to the new interface location. In practice, the grid never deforms past one cell, so the search can be
        !simplified by this. Also, the interval of integration is usually of magnitude close to zero or close to dpo because of minimial deformation.
        !Numerous tests confirmed that the bottom and top of the grids match to machine precision, so I set them equal to each other.
        do k = 1 , nlev
          kk = k  !Keep from an order n^2 search operation by assuming the old cell index is close.
          !Find the index of the old grid cell in which this new cell's bottom interface resides.
          do while ( pio(kk) <= pin(k+1) )
            kk = kk + 1
          enddo
          kk = kk - 1                   !kk is now the cell index we're integrating over.
          if (kk == nlev+1) kk = nlev   !This is to keep the indices in bounds. Top bounds match anyway, so doesn't matter what coefficients are used
          kid(k) = kk                   !Save for reuse
          z1(k) = -0.5D0                !This remapping assumes we're starting from the left interface of an old grid cell
                                        !In fact, we're usually integrating very little or almost all of the cell in question
          z2(k) = ( pin(k+1) - ( pio(kk) + pio(kk+1) ) * 0.5 ) / dpo(kk)  !PPM interpolants are normalized to an independent coordinate domain [-0.5,0.5].
        enddo

        !This turned out a big optimization, remembering that only parts of the PPM algorithm depends on the data, namely the limiting. So anything that
        !depends only on the grid is pre-computed outside the tracer loop.
        ppmdx(:,:) = compute_ppm_grids( dpo )

        !From here, we loop over tracers for only those portions which depend on tracer data, which includes PPM limiting and mass accumulation
        do q = 1 , qsize
          !Accumulate the old mass up to old grid cell interface locations to simplify integration during remapping. Also, divide out the
          !grid spacing so we're working with actual tracer values and can conserve mass. The option for ifndef ZEROHORZ I believe is there
          !to ensure tracer consistency for an initially uniform field. I copied it from the old remap routine.
          masso(1) = 0.
          do k = 1 , nlev
#ifdef ZEROHORZ
            ao(k) = elem(ie)%state%Qdp(i,j,k,q,n0)
#else
            ao(k) = ( elem(ie)%state%Qdp(i,j,k,q,n0) + (rkstage-1) * elem(ie)%state%Qdp(i,j,k,q,np1) ) / rkstage
#endif
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
          !Compute tracer values on the new grid by integrating from the old cell bottom to the new cell interface to form a new grid mass accumulation
          !Taking the difference between accumulation at successive interfaces gives the mass inside each cell. Since Qdp is supposed to hold the full mass
          !this needs no normalization.
          massn1 = 0.
          do k = 1 , nlev
            kk = kid(k)
            massn2 = masso(kk) + integrate_parabola( coefs(:,kk) , z1(k) , z2(k) ) * dpo(kk)
            elem(ie)%state%Qdp(i,j,k,q,np1) = massn2 - massn1
            massn1 = massn2
          enddo
        enddo
      enddo
    enddo
  enddo
  call t_stopf('remap_velocityQ_ppm')
end subroutine remap_velocityQ_ppm




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
  real(kind=real_kind) :: dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8, dx9, dx10 !Holds expressions based on the grid which are cumbersome
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
    ai(j) = a(j) + dx(4,j) * ( a(j+1) - a(j) ) + dx(5,j) * ( dx(6,j) * ( dx(7,j) - dx(8,j) ) * ( a(j+1) - a(j) ) - dx(9,j) * dma(j+1) + dx(10,j) * dma(j) )
  enddo

  ! Stage 3: Compute limited PPM interpolant over each cell in the physical domain (dimension nlev) using ai on either side and ao within the cell.
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

  !If we're not using a mirrored boundary condition, then make the two cells bordering the top and bottom material boundaries piecewise constant.
  !Zeroing out the first and second moments, and setting the zeroth moment to the cell mean is sufficient to maintain conservation.
  if (vert_remap_q_alg == 2) then
    coefs(0,1:2) = a(1:2)
    coefs(1:2,1:2) = 0.
    coefs(0,nlev-1:nlev) = a(nlev-1:nlev)
    coefs(1:2,nlev-1:nlev) = 0.D0
  endif
end function compute_ppm




!Simple function computes the definite integral of a parabola in normalized coordinates, xi=(x-x0)/dx, given two bounds. Make sure this gets inlined during compilation.
function integrate_parabola( a , x1 , x2 )    result(mass)
  implicit none
  real(kind=real_kind), intent(in) :: a(0:2)  !Coefficients of the parabola
  real(kind=real_kind), intent(in) :: x1      !lower domain bound for integration
  real(kind=real_kind), intent(in) :: x2      !upper domain bound for integration
  real(kind=real_kind)             :: mass
  mass = a(0) * (x2 - x1) + a(1) * (x2 ** 2 - x1 ** 2) / 0.2D1 + a(2) * (x2 ** 3 - x1 ** 3) / 0.3D1
end function integrate_parabola


 




  subroutine remap_velocityC(n0,np1,dt,elem,fvm,hvcoord,nets,nete,compute_diagnostics)
  
    use physical_constants, only : cp, cpwater_vapor
	
    implicit none
    real (kind=real_kind),  intent(in)        :: dt
    type (element_t),    intent(inout), target  :: elem(:)
    type (fvm_struct),    intent(inout), target  :: fvm(:)
    type (hvcoord_t),    intent(in)        :: hvcoord
    logical,        intent(in)              :: compute_diagnostics
    
    integer :: nets,nete,n0,np1
    
    ! ========================
    ! Local Variables
    ! ========================

    real (kind=real_kind), dimension(nlev+1)    :: rhs,lower_diag,diag,upper_diag,q_diag,zgam,z1c,z2c,zv
    real (kind=real_kind), dimension(nlev)      :: h,Qcol,dy,za0,za1,za2,zarg,zhdp
    real (kind=real_kind)  :: dp_star,dp_np1,f_xm,level1,level2,level3,level4,level5, &
								peaks_min,peaks_max,Q_vadv,tmp_cal,xm,xm_d,zv1,zv2, &
								zero = 0,one = 1,tiny = 1e-12,qmax = 1d50
    
    integer(kind=int_kind) :: zkr(nlev+1),filter_code(nlev),peaks,im1,im2,im3,ip1,ip2, & 
									lt1,lt2,lt3,t0,t1,t2,t3,t4,tm,tp,ie,i,ilev,j,jk,k,q
    
    call t_startf('remap_velocityC')

    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,i,j,z1c,z2c,zv,k,dp_np1,dp_star,Qcol,zkr,ilev) &
!$omp    private(jk,zgam,zhdp,h,zarg,rhs,lower_diag,diag,upper_diag,q_diag,tmp_cal,filter_code) &
!$omp    private(dy,im1,im2,im3,ip1,t1,t2,t3,za0,za1,za2,xm_d,xm,f_xm,t4,tm,tp,peaks,peaks_min) &
!$omp    private(peaks_max,ip2,level1,level2,level3,level4,level5,lt1,lt2,lt3,zv1,zv2,Q_vadv)
#endif
      do q=1,ntrac
        do i=1,np
          do j=1,np
             z1c(1)=0
             z2c(1)=0
             zv(1)=0
             do k=1,nlev
                dp_np1 = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                     (hvcoord%hybi(k+1) - hvcoord%hybi(k))*elem(ie)%state%ps_v(i,j,np1)

#ifdef ZEROVERT                        
                dp_star = dp_np1  ! ignore the vertical motion
#else
                dp_star = dp_np1 + dt*(elem(ie)%derived%eta_dot_dpdn(i,j,k+1) & 
                     -elem(ie)%derived%eta_dot_dpdn(i,j,k)) 
#endif        
                z1c(k+1) = z1c(k)+dp_star
                z2c(k+1) = z2c(k)+dp_np1
#ifdef ZEROHORZ			  
                Qcol(k)=fvm(ie)%c(i,j,k,q,n0)  ! ignore horizontal motion
#else		
                Qcol(k)=fvm(ie)%c(i,j,k,q,np1)
#endif			
                zv(k+1) = zv(k)+Qcol(k)
             enddo
             
             if (ABS(z2c(nlev+1)-z1c(nlev+1)).GE.0.000001) then
                write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
                write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
                write(6,*) 'DATA FOR MODEL LEVELS'
                write(6,*) 'PLEVMODEL=',z2c(nlev+1)
                write(6,*) 'PLEV     =',z1c(nlev+1)
                write(6,*) 'DIFF     =',z2c(nlev+1)-z1c(nlev+1)
                ! call ABORT
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
				endif
                zv2 = zv(zkr(k+1))+(za0(zkr(k+1))*zgam(k+1)+(za1(zkr(k+1))/2)*(zgam(k+1)**2)+ &
                     (za2(zkr(k+1))/3)*(zgam(k+1)**3))*zhdp(zkr(k+1))
                Q_vadv = (fvm(ie)%c(i,j,k,q,np1) - (zv2 - zv1)) / dt
                fvm(ie)%c(i,j,k,q,np1) = (zv2 - zv1)	
                zv1 = zv2
#ifdef ENERGY_DIAGNOSTICS
                if (compute_diagnostics .and. q==1) then
                   elem(ie)%accum%IEvert1_wet(i,j) = elem(ie)%accum%IEvert1_wet(i,j) + (Cpwater_vapor-Cp)*elem(ie)%state%T(i,j,k,n0)*Q_vadv
                endif
#endif	
             enddo
          enddo
       enddo
    enddo
 enddo
 
 call t_stopf('remap_velocityC')
 end subroutine remap_velocityC
 
 end module vertremap_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End GPU remap module    !!
!! by Rick Archibald, 2010  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






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
  use dimensions_mod, only     : nlev, nlevp, np, qsize, ntrac
  use physical_constants, only : rgas, Rwater_vapor, kappa, g, rearth, rrearth, cp, cpwater_vapor
  use derivative_mod, only     : gradient, vorticity, gradient_wk, derivative_t, divergence, &
                                 gradient_sphere, divergence_sphere
  use element_mod, only        : element_t
  use fvm_control_volume_mod, only        : fvm_struct
  use filter_mod, only         : filter_t, filter_P
  use hybvcoord_mod, only      : hvcoord_t
  use time_mod, only           : TimeLevel_t, smooth
  use prim_si_mod, only        : preq_pressure
  use diffusion_mod, only      : scalar_diffusion, diffusion_init
  use control_mod, only        : integration, test_case, filter_freq_advection,  hypervis_order, &
        statefreq, moisture, TRACERADV_TOTAL_DIVERGENCE, TRACERADV_UGRADQ, &
        tracer_advection_formulation, prescribed_wind, nu_q, nu_p, limiter_option, hypervis_subcycle_q
  use edge_mod, only           : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack, initedgebuffer, edgevunpackmin
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use viscosity_mod, only      : biharmonic_wk_scalar, biharmonic_wk_scalar_minmax, neighbor_minmax
  use perf_mod, only           : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only   : abortmp  	

  implicit none
  
  private  

  public :: Prim_Advec_Init, Prim_Advec_Tracers_remap_rk2, Prim_Advec_Tracers_lf
  public :: prim_advec_tracers_fvm

  type (EdgeBuffer_t) :: edgeAdv, edgeAdvQ3, edgeAdv_p1, edgeAdvQ2, edgeAdv1

  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  real(kind=real_kind), allocatable :: qmin(:,:,:), qmax(:,:,:)


contains

  subroutine Prim_Advec_Init()
    use dimensions_mod, only : nlev, qsize, nelemd

    ! this might be called with qsize=0
    ! allocate largest one first
    call initEdgeBuffer(edgeAdvQ3,max(nlev,qsize*nlev*3))  ! Qtens,Qmin, Qmax
    ! remaining edge buffers can share %buf and %receive with edgeAdvQ3:
    call initEdgeBuffer(edgeAdv1,nlev,edgeAdvQ3%buf,edgeAdvQ3%receive)
    call initEdgeBuffer(edgeAdv,qsize*nlev,edgeAdvQ3%buf,edgeAdvQ3%receive)
    call initEdgeBuffer(edgeAdv_p1,qsize*nlev + nlev,edgeAdvQ3%buf,edgeAdvQ3%receive) 
    call initEdgeBuffer(edgeAdvQ2,qsize*nlev*2,edgeAdvQ3%buf,edgeAdvQ3%receive)  ! Qtens,Qmin, Qmax




    ! this static array is shared by all threads, so dimension for all threads (nelemd), not nets:nete:
    allocate (qmin(nlev,qsize,nelemd))
    allocate (qmax(nlev,qsize,nelemd))

  end subroutine Prim_Advec_Init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fvm driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Prim_Advec_Tracers_fvm(elem, fvm, deriv,hvcoord,hybrid,&
        dt,tl,nets,nete, compute_diagnostics)
    use perf_mod, only : t_startf, t_stopf            ! _EXTERNAL
    use vertremap_mod, only: remap_velocityC,remap_velocityUV  ! _EXTERNAL (actually INTERNAL)
    use fvm_mod, only : cslam_run, cslam_runairdensity, edgeveloc, fvm_mcgregor, fvm_mcgregordss
    
    implicit none
    type (element_t), intent(inout)   :: elem(:)
    type (fvm_struct), intent(inout)   :: fvm(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t),     intent(in):: hybrid
    type (TimeLevel_t)                :: tl

    real(kind=real_kind) , intent(in) :: dt
    integer,intent(in)                :: nets,nete

    logical,intent(in)                :: compute_diagnostics


    real (kind=real_kind), dimension(np,np,nlev)    :: dp_star
    real (kind=real_kind), dimension(np,np,nlev)    :: dp_np1
   
    integer :: n0,np1,ie,k
    
    real (kind=real_kind)  :: vstar(np,np,2)
    real (kind=real_kind)  :: vhat(np,np,2)
    real (kind=real_kind), dimension(np, np) :: v1, v2
    

    call t_barrierf('sync_prim_advec_tracers_fvm', hybrid%par%comm)
    call t_startf('prim_advec_tracers_fvm')
    n0  = tl%n0
    np1 = tl%np1

!     do ie=nets,nete
!       do k=1,nlev
!         do i=1,np
!           do j=1,np      
!           elem(ie)%derived%dp(i,j,k)=
!           enddo
!         enddo
!       enddo
!     enddo


! dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
!      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,n0)

    ! mean vertical velocity computed in dynamics needs to be DSS'd:
    ! this should be moved into one of the DSS operations performed by the
    ! dynamics to be more efficient, but only if ntrac>0
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
!        elem(ie)%derived%eta_dot_dpdn(:,:,:)=0
    end do
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2D advection step
    ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! elem%state%u(np1)  = velocity at time t+1 on reference levels
    ! elem%derived%vstar = velocity at t+1 on floating levels (computed below)
    call remap_velocityUV(np1,dt,elem,hvcoord,nets,nete)
!------------------------------------------------------------------------------------    
    call t_startf('fvm_mcgregor')
    ! using McGregor AMS 1993 scheme: Economical Determination of Departure Points for
    ! Semi-Lagrangian Models 
!     do ie=nets,nete
!       do k=1,nlev
!         vstar=elem(ie)%derived%vstar(:,:,:,k) 
!         vhat=(fvm(ie)%vn0(:,:,:,k) + elem(ie)%derived%vstar(:,:,:,k))/2
!         ! calculate high order approximation
!         call fvm_mcgregor(elem(ie), deriv, dt, vhat, vstar,1)
!         ! apply DSS to make vstar C0
!         elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%spheremp(:,:)*vstar(:,:,1) 
!         elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%spheremp(:,:)*vstar(:,:,2) 
!       enddo 
!       call edgeVpack(edgeveloc,elem(ie)%derived%vstar(:,:,1,:),nlev,0,elem(ie)%desc)
!       call edgeVpack(edgeveloc,elem(ie)%derived%vstar(:,:,2,:),nlev,nlev,elem(ie)%desc)
!     enddo 
!     call bndry_exchangeV(hybrid,edgeveloc)
!     do ie=nets,nete
!        call edgeVunpack(edgeveloc,elem(ie)%derived%vstar(:,:,1,:),nlev,0,elem(ie)%desc)
!        call edgeVunpack(edgeveloc,elem(ie)%derived%vstar(:,:,2,:),nlev,nlev,elem(ie)%desc)
!        do k=1, nlev  
!          elem(ie)%derived%vstar(:,:,1,k)=elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%rspheremp(:,:)
!          elem(ie)%derived%vstar(:,:,2,k)=elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%rspheremp(:,:)
!        end do
!     end do
!     do ie=nets,nete
!       do k=1,nlev
!          ! Convert wind to lat-lon
!         v1     = fvm(ie)%vn0(:,:,1,k)  ! contra
!         v2     = fvm(ie)%vn0(:,:,2,k)  ! contra 
!         fvm(ie)%vn0(:,:,1,k)=elem(ie)%D(1,1,:,:)*v1 + elem(ie)%D(1,2,:,:)*v2   ! contra->latlon
!         fvm(ie)%vn0(:,:,2,k)=elem(ie)%D(2,1,:,:)*v1 + elem(ie)%D(2,2,:,:)*v2   ! contra->latlon
!     
! !          ! Convert wind to lat-lon
! !         v1     = elem(ie)%derived%vstar(:,:,1,k)  ! contra
! !         v2     = elem(ie)%derived%vstar(:,:,2,k)   ! contra 
! !         elem(ie)%derived%vstar(:,:,1,k)=elem(ie)%D(1,1,:,:)*v1 + elem(ie)%D(1,2,:,:)*v2   ! contra->latlon
! !         elem(ie)%derived%vstar(:,:,2,k)=elem(ie)%D(2,1,:,:)*v1 + elem(ie)%D(2,2,:,:)*v2   ! contra->latlon
!     
!       enddo  
!     end do
!     call fvm_mcgregordss(elem,fvm,nets,nete, hybrid, deriv, dt, 3)
    call t_stopf('fvm_mcgregor')

!------------------------------------------------------------------------------------    
    
    

    ! fvm departure calcluation should use vstar.
    ! from c(n0) compute c(np1): 
    call cslam_runairdensity(elem,fvm,hybrid,deriv,dt,tl,nets,nete)

    ! apply vertical remap back to reference levels
    call remap_velocityC(n0,np1,dt,elem,fvm,hvcoord,nets,nete,compute_diagnostics)


    call t_stopf('prim_advec_tracers_fvm')
  end subroutine Prim_Advec_Tracers_fvm



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
!    Non-consistent scheme used with leapfrog dynamics, no subcycling
!       vn0()           U at timelevel n0 
!       eta_dot_dpdn()  mean vertical velocity
! 
!
! 3 stage
!    Euler step from t     -> t+.5
!    Euler step from t+.5  -> t+1.0
!    Euler step from t+1.0 -> t+1.5
!    u(t) = u(t)/3 + u(t+2)*2/3
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine Prim_Advec_Tracers_remap_rk2(elem, deriv,hvcoord,flt,hybrid,&
        dt,tl,nets,nete, compute_diagnostics)
    use perf_mod, only : t_startf, t_stopf            ! _EXTERNAL
    use derivative_mod, only : divergence_sphere
    use vertremap_mod, only: remap_velocityuv,remap_velocityq, remap_velocityq_ppm  ! _EXTERNAL (actually INTERNAL)
    use control_mod, only: vert_remap_q_alg

    implicit none
    type (element_t), intent(inout)   :: elem(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (filter_t)                   :: flt
    type (hybrid_t),     intent(in):: hybrid
    type (TimeLevel_t)                :: tl

    real(kind=real_kind) , intent(in) :: dt
    integer,intent(in)                :: nets,nete

    logical,intent(in)                :: compute_diagnostics
    real (kind=real_kind), dimension(np,np,2)    :: gradQ


    real (kind=real_kind), dimension(np,np,nlev) :: Q_vadv 
    real (kind=real_kind), dimension(np,np,nlev)    :: dp_star
    real (kind=real_kind), dimension(np,np,nlev)    :: dp_np1
   
    integer :: i,j,k,l,ie,q,nmin
    integer :: n0,np1,nfilt,rkstage,rhs_multiplier


    call t_barrierf('sync_prim_advec_tracers_remap_k2', hybrid%par%comm)
    call t_startf('prim_advec_tracers_remap_rk2')
    n0  = tl%n0
    np1 = tl%np1
    rkstage=3 !   3 stage RKSSP scheme, with optimal SSP CFL

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
#if (defined ELEMENT_OPENMP)
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
      call euler_step(np1,n0,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,&
           compute_diagnostics,DSSdiv_vdp_ave,rhs_multiplier)
      
      rhs_multiplier = 1
      call euler_step(np1,np1,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,&
           .false.,DSSeta,rhs_multiplier)
      
      rhs_multiplier = 2
      call euler_step(np1,np1,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,&
           .false.,DSSomega,rhs_multiplier)
      

    ! to finish the 2D advection step, we need to average the t and t+2 results
    ! to get a second order estimate for t+1.  We then apply the vertical
    ! remap.  These two steps have been merged into one for efficienty:
    if (vert_remap_q_alg == 0) then
      call remap_velocityQ(n0,np1,dt,elem,hvcoord,nets,nete,compute_diagnostics,rkstage)
    elseif (vert_remap_q_alg == 1 .or. vert_remap_q_alg == 2) then
      call remap_velocityQ_ppm(n0,np1,dt,elem,hybrid,hvcoord,nets,nete,compute_diagnostics,rkstage)
    else
      call abortmp('specification for vert_remap_q_alg must be 0, 1, or 2.')
    endif


#ifndef ZEROHORZ
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Dissipation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (limiter_option == 8 .or. nu_p>0) then
       ! dissipation was applied in RHS.  
    else
       call advance_hypervis_scalar(edgeadv,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt)
    endif
#endif

    ! update Q from Qdp
    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q)
#endif
       do k=1,nlev    !  Loop inversion (AAM)
          dp_np1(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
          do q=1,qsize
            elem(ie)%state%Q(:,:,k,q,np1)=elem(ie)%state%Qdp(:,:,k,q,np1)/dp_np1(:,:,k) 
          enddo
       enddo
    enddo


    call t_stopf('prim_advec_tracers_remap_rk2')
  end subroutine prim_advec_tracers_remap_rk2

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine Prim_Advec_Tracers_lf(elem, deriv,hvcoord,flt,hybrid,dt,tl,nets,nete, compute_diagnostics)
    use perf_mod, only : t_startf, t_stopf              ! _EXTERNAL
    implicit none
    type (element_t), intent(inout)   :: elem(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (filter_t)                   :: flt

    type (hybrid_t),     intent(in):: hybrid

    real(kind=real_kind) , intent(in) :: dt
    type (TimeLevel_t)                :: tl

    integer,intent(in)                :: nets,nete
    logical,intent(in)                :: compute_diagnostics

    ! local
    real(kind=real_kind) :: dp, pmid, dt2
    integer :: i,j,k,l,ie,q
    integer :: n0,nm1,np1,nfilt, nstep


    call t_barrierf('sync_prim_advec_tracers_lf', hybrid%par%comm)
    call t_startf('prim_advec_tracers_lf')

    n0  = tl%n0
    nm1 = tl%nm1
    np1 = tl%np1

    nfilt = n0
    nstep = tl%nstep
    dt2 = 2 * dt


    if (filter_freq_advection>0) then
    if(nstep>0 .and.  modulo(nstep,filter_freq_advection)==0) then
       do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q)
#endif
          do q=1,qsize	
             do k=1,nlev    !  Potential loop inversion (AAM)
                call filter_P(elem(ie)%state%Q(:,:,k,q,nfilt),flt)
                elem(ie)%state%Q(:,:,k,q,nfilt) = elem(ie)%mp(:,:)*elem(ie)%state%Q(:,:,k,q,nfilt)
             end do
          end do
          call edgeVpack(edgeadv,elem(ie)%state%Q(:,:,:,:,nfilt),nlev*qsize,0,elem(ie)%desc)
       end do
       call bndry_exchangeV(hybrid,edgeadv)
       do ie=nets,nete
          call edgeVunpack(edgeadv,elem(ie)%state%Q(:,:,:,:,nfilt),nlev*qsize,0,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q)    !  Potential loop inversion (AAM)
#endif
          do q=1,qsize	
             do k=1,nlev
                elem(ie)%state%Q(:,:,k,q,nfilt) = elem(ie)%rmp(:,:)*elem(ie)%state%Q(:,:,k,q,nfilt)
             enddo
          end do
       end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
    end if
    end if


    if (nstep==0) then
       ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))
       call compute_and_apply_rhs(np1,n0,n0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.)
       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))
       call compute_and_apply_rhs(np1,n0,np1,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics)
    else
       call compute_and_apply_rhs(np1,nm1,n0,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics)
    endif


    if (hypervis_order==0) then
       call scalar_diffusion(elem, nets,nete,np1,deriv,dt2,hybrid)
    else
       ! hypervis_order==1:  weak form laplacian  (only 1 DSS)
       ! hypervis_order==2:  weak form biharmonic (2 DSS's) 
       call advance_hypervis_scalar_lf(edgeadv,elem,hvcoord,hybrid,deriv,np1,n0,nets,nete,dt2)
    endif

    if (nu_q>0) then
    ! if nu_q=0, we are running an inviscid test, skip fixer
    !
    ! apply negative Q fixer
    !
    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q)
#endif
       do q=1,qsize
          do k=1,nlev
             elem(ie)%state%Q(:,:,k,q,np1) = elem(ie)%spheremp(:,:)*elem(ie)%state%Q(:,:,k,q,np1)
          enddo

          ! limiter3d_noncon: no negative values, even if mass is added
          call limiter3d_noncon(elem(ie)%state%Q(:,:,:,q,np1),&
              elem(ie)%state%ps_v(:,:,n0),&
              hvcoord,elem(ie)%accum%mass_added(q))

       end do
       call edgeVpack(edgeadv,elem(ie)%state%Q(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
    end do
    call bndry_exchangeV(hybrid,edgeadv)
    do ie=nets,nete
       call edgeVunpack(edgeadv,elem(ie)%state%Q(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q)
#endif
       do q=1,qsize	
          do k=1,nlev    !  Potential loop inversion (AAM)
             elem(ie)%state%Q(:,:,k,q,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%Q(:,:,k,q,np1)
          enddo
       end do
    end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
    endif

    call t_stopf('prim_advec_tracers_lf')

    if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
       do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q,j,i,dp)
#endif
          do q=1,qsize
             do k=1,nlev
                do j=1,np
                   do i=1,np
                      ! timestep was done in Q.  copy over to Qdp:                                              
                      elem(ie)%state%Qdp(i,j,k,q,np1)=elem(ie)%state%Q(i,j,k,q,np1)
                      ! recompute Q from dpQ for consistency                                                    
                      dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                           ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,n0)
                      elem(ie)%state%Q(i,j,k,q,np1) =elem(ie)%state%Qdp(i,j,k,q,np1)/dp
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
       

  end subroutine Prim_Advec_Tracers_lf

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine compute_and_apply_rhs(np1,nm1,n0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine is normally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  ! if  dt2<0, then the DSS'd RHS is returned in timelevel np1
  !
  ! Combining the RHS and DSS pack operation in one routine 
  ! allows us to fuse these two loops for more cache reuse
  !
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : np, nlev
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use hybvcoord_mod, only : hvcoord_t

  implicit none
  integer :: np1,nm1,n0,nets,nete
  real (kind=real_kind) :: dt2
  logical  :: compute_diagnostics

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv

  ! local
  real (kind=real_kind) ::  qtens(np,np,nlev)
  real (kind=real_kind) ::  Q_vadv(np,np,nlev)
  real (kind=real_kind) ::  rpdel(np,np,nlev)
  real (kind=real_kind) ::  gradQ(np,np,2)
  real (kind=real_kind) ::  divdp(np,np)
  real (kind=real_kind) ::  v1,v2
  integer :: i,j,k,ie,q

  logical ::  use_explicit_eta_dot=.true. ! reuse from dynamics or recompute?


  do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q,j,i,gradQ,divdp,qtens,v1,v2,Q_vadv,rpdel)
#endif
     do q=1,qsize
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   2D contribution
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
           do k=1,nlev
              ! div( U dp Q),                                                                                 
              gradQ(:,:,1)=elem(ie)%state%v(:,:,1,k,n0)*elem(ie)%state%Qdp(:,:,k,q,n0)
              gradQ(:,:,2)=elem(ie)%state%v(:,:,2,k,n0)*elem(ie)%state%Qdp(:,:,k,q,n0)
              divdp = divergence_sphere(gradQ,deriv,elem(ie))
              do j=1,np
                 do i=1,np
                    qtens(i,j,k)=-divdp(i,j)
                 enddo
              enddo
           enddo
        else
           !   UGRADQ formulation
           do k=1,nlev
              gradQ = gradient_sphere(elem(ie)%state%Q(:,:,k,q,n0),deriv,elem(ie)%Dinv)
              do j=1,np	
                 do i=1,np
                    v1    = elem(ie)%state%v(i,j,1,k,n0)
                    v2    = elem(ie)%state%v(i,j,2,k,n0)
                    Qtens(i,j,k) = -( v1*gradQ(i,j,1) + v2*gradQ(i,j,2)  )
                 enddo
              enddo
           enddo
        endif
        


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! vertical advection
        ! evaluate at np1 for time-split scheme
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
           call preq_vertadv_dpQ(elem(ie)%state%Q(:,:,:,q,n0),elem(ie)%derived%eta_dot_dpdn,Q_vadv)
           ! advance in time, into Q, apply mass matrix
           do k=1,nlev
              elem(ie)%state%Q(:,:,k,q,np1) = elem(ie)%spheremp(:,:)*&
                   ( elem(ie)%state%Qdp(:,:,k,q,nm1)  + &
                   dt2*(qtens(:,:,k)-Q_vadv(:,:,k)) )
           enddo
        else
           if (use_explicit_eta_dot) then
              ! vertical advection term, using eta_dot_dpdn from advection timestep
              do k=1,nlev
                 do j=1,np
                    do i=1,np
                       v1 = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%ps_v(i,j,n0)
                       v2 = hvcoord%hyai(k+1)*hvcoord%ps0 + hvcoord%hybi(k+1)*elem(ie)%state%ps_v(i,j,n0)
                       rpdel(i,j,k) = 1.0D0/(v2-v1)
                    end do
                 end do
              enddo
              call preq_vertadvQ(elem(ie)%state%Q(:,:,:,q,n0),elem(ie)%derived%eta_dot_dpdn,rpdel,Q_vadv)
           else
              ! recompute eta_dot_dpdn using velocity at level n0
              ! and elem%derived%grad_lnps
              call preq_impsysQ(elem(ie),hvcoord,np1,n0,nm1,elem(ie)%state%Q(:,:,:,q,n0),Q_vadv)
           endif
           ! advance in time, apply mass matrix
           do k=1,nlev
              elem(ie)%state%Q(:,:,k,q,np1) = elem(ie)%spheremp(:,:)*&
                   ( elem(ie)%state%Q(:,:,k,q,nm1)  + &
                   dt2*(qtens(:,:,k)-Q_vadv(:,:,k)) )
           enddo
        endif
        
        if (nu_q>0) then
           ! if nu_q=0, we are running an inviscid test, skip fixer
           call limiter2d_zero(elem(ie)%state%Q(:,:,:,q,np1),&
             elem(ie)%state%ps_v(:,:,n0),hvcoord)
        endif
        
     end do
     call edgeVpack(edgeadv,elem(ie)%state%Q(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
  end do
  call bndry_exchangeV(hybrid,edgeadv)
  
  do ie=nets,nete
     call edgeVunpack(edgeadv,elem(ie)%state%Q(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
     
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, q)
#endif
     do q=1,qsize	
        do k=1,nlev    !  Potential loop inversion (AAM)
           elem(ie)%state%Q(:,:,k,q,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%Q(:,:,k,q,np1)
        enddo
     end do
  end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
end subroutine compute_and_apply_rhs

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine euler_step(np1,n0,dt,elem,hvcoord,hybrid,deriv,nets,nete,&
      compute_diagnostics,DSSopt,rhs_multiplier)
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
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use hybvcoord_mod, only : hvcoord_t

  implicit none
  integer :: np1,nm1,n0,nets,nete,DSSopt,rhs_multiplier
  real (kind=real_kind), intent(in)  :: dt
  logical  :: compute_diagnostics

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv

  ! local
  real (kind=real_kind), dimension(np,np)    :: divdp, dpdiss
  real (kind=real_kind), dimension(np,np,2)    :: gradQ
  real(kind=real_kind), dimension(np,np,nlev) :: Qtens
  real(kind=real_kind), dimension(np,np,nlev) :: dp,dp_star
  real(kind=real_kind), dimension(np,np,2,nlev) :: Vstar
  real (kind=real_kind), pointer, dimension(:,:,:)   :: DSSvar
  real(kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: Qtens_biharmonic
! nelemd

  real(kind=real_kind) :: dp0
  integer :: ie,q,i,j,k
  integer :: rhs_viss=0

  if (npdg>0) then
      call euler_step_dg(np1,n0,dt,elem,hvcoord,hybrid,deriv,nets,nete,&
           compute_diagnostics,DSSopt,rhs_multiplier)
     return
  endif

!   call t_barrierf('sync_euler_step', hybrid%par%comm)
  call t_startf('euler_step')


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   compute Q min/max values for lim8
  !   compute biharmonic mixing term f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rhs_viss=0
  if (limiter_option == 8 .or. nu_p>0) then
     !
     ! for limiter=0,4 or 8 we will apply dissipation in the RHS,
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


     ! initialize dp, and compute Q from Qdp (and store Q in Qtens_biharmonic)
     do ie=nets,nete
        ! add hyperviscosity to RHS.  apply to Q at timelevel n0, Qdp(n0)/dp
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, q)
#endif
        do k=1,nlev    !  Loop index added with implicit inversion (AAM)
           dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - &
                rhs_multiplier*dt*elem(ie)%derived%divdp_proj(:,:,k) 
           do q=1,qsize
              Qtens_biharmonic(:,:,k,q,ie) = elem(ie)%state%Qdp(:,:,k,q,n0)/dp(:,:,k)
           enddo
        enddo
     enddo

     ! compute element qmin/qmax
     if (rhs_multiplier == 0) then
        do ie=nets,nete
           do k=1,nlev    
              do q=1,qsize
                 qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
                 qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
                 qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
              enddo
           enddo
        enddo
        ! update qmin/qmax based on neighbor data for lim8
        if (limiter_option==8) &
             call neighbor_minmax(elem,hybrid,edgeAdvQ2,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
     endif

     ! lets just reuse the old neighbor min/max, but update based on local data
     if (rhs_multiplier == 1 ) then
        do ie=nets,nete
           do k=1,nlev    !  Loop index added with implicit inversion (AAM)
              do q=1,qsize
                 qmin(k,q,ie)=min(qmin(k,q,ie),minval(Qtens_biharmonic(:,:,k,q,ie)))
                 qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
                 qmax(k,q,ie)=max(qmax(k,q,ie),maxval(Qtens_biharmonic(:,:,k,q,ie)))
              enddo
           enddo
        enddo
     endif

     ! get niew min/max values, and also compute biharmonic mixing term
     if (rhs_multiplier == 2) then
        rhs_viss=3
        ! compute element qmin/qmax  
        do ie=nets,nete
           do k=1,nlev    
              do q=1,qsize
                 qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
                 qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
                 qmin(k,q,ie)=max(qmin(k,q,ie),0d0)
              enddo
           enddo
        enddo
        ! two scalings depending on nu_p:
        ! nu_p=0:    qtens_biharmonic *= dp0                   (apply viscsoity only to q)
        ! nu_p>0):   qtens_biharmonc *= elem()%psdiss_ave      (for consistency, if nu_p=nu_q)
        if (nu_p>0) then
           do ie=nets,nete
#if (defined ELEMENT_OPENMP)
              !$omp parallel do private(k, q, dp0)
#endif
              do k=1,nlev    
                 dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
                 dpdiss(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%derived%psdiss_ave(:,:)
                 do q=1,qsize
                    ! NOTE: divide by dp0 since we multiply by dp0 below
                    Qtens_biharmonic(:,:,k,q,ie)=Qtens_biharmonic(:,:,k,q,ie)*dpdiss(:,:)/dp0
                 enddo
              enddo
           enddo
        endif
        if (limiter_option==8) then
           ! biharmonic and update neighbor min/max
           call biharmonic_wk_scalar_minmax(elem,qtens_biharmonic,deriv,edgeAdvQ3,hybrid,&
             nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
        else
           ! regular biharmonic, no need to updat emin/max
           call biharmonic_wk_scalar(elem,qtens_biharmonic,deriv,edgeAdv,hybrid,nets,nete)
        endif
        do ie=nets,nete
#if (defined ELEMENT_OPENMP)
        !$omp parallel do private(k, q, dp0)
#endif
           do k=1,nlev    !  Loop inversion (AAM)
              dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                   ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0

              do q=1,qsize
                 ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
                 qtens_biharmonic(:,:,k,q,ie) =  -rhs_viss*dt*nu_q*dp0*Qtens_biharmonic(:,:,k,q,ie)&
                      / elem(ie)%spheremp(:,:)
              enddo
           enddo
        enddo
     endif
  endif  ! compute biharmonic mixing term and qmin/qmax




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete

     ! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
     ! all zero so we only have to DSS 1:nlev
     if ( DSSopt == DSSeta) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
     if ( DSSopt == DSSomega) DSSvar => elem(ie)%derived%omega_p(:,:,:)
     if ( DSSopt == DSSdiv_vdp_ave) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)



     ! Compute velocity used to advance Qdp 
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
      do k=1,nlev    !  Loop index added (AAM)
         ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
         ! but that's ok because rhs_multiplier=0 on the first stage:
         dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - &
              rhs_multiplier*dt*elem(ie)%derived%divdp_proj(:,:,k) 
         Vstar(:,:,1,k) = elem(ie)%derived%vn0(:,:,1,k)/dp(:,:,k)
         Vstar(:,:,2,k) = elem(ie)%derived%vn0(:,:,2,k)/dp(:,:,k)
      enddo
     

     ! advance Qdp
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k,gradQ,dp_star,qtens)
#endif
     do q=1,qsize

        do k=1,nlev  !  dp_star used as temporary instead of divdp (AAM)
           ! div( U dp Q), 
           gradQ(:,:,1)=Vstar(:,:,1,k)*elem(ie)%state%Qdp(:,:,k,q,n0)
           gradQ(:,:,2)=Vstar(:,:,2,k)*elem(ie)%state%Qdp(:,:,k,q,n0)
           dp_star(:,:,k) = divergence_sphere(gradQ,deriv,elem(ie))
           Qtens(:,:,k)=elem(ie)%state%Qdp(:,:,k,q,n0) - dt*dp_star(:,:,k)
           ! optionally add in hyperviscosity computed above:
           if (rhs_viss/=0) Qtens(:,:,k) = Qtens(:,:,k) + Qtens_biharmonic(:,:,k,q,ie)
        enddo
           
#ifdef ENERGY_DIAGNOSTICS
        if (compute_diagnostics .and. q==1) then
          do k=1,nlev  !  dp_star used as temporary instead of divdp (AAM)
              ! IEvert1_wet():  (Cpv-Cp) T Qdp_vadv  (Q equation)
              ! IEhorz1_wet():  (Cpv-Cp) T Qdp_hadv  (Q equation)
              elem(ie)%accum%IEhorz1_wet(:,:) = elem(ie)%accum%IEhorz1_wet(:,:) +&
                   (Cpwater_vapor-Cp)*elem(ie)%state%T(:,:,k,n0)*dp_star(:,:,k)
          enddo
        endif
#endif

        if(limiter_option == 8)then

           do k=1,nlev  ! Loop index added (AAM)
              ! UN-DSS'ed dp at timelevel n0+1:  
              dp_star(:,:,k) = dp(:,:,k) - dt*elem(ie)%derived%divdp(:,:,k)  
              if (nu_p>0 .and.  rhs_viss/=0) then
                 ! add contribution from UN-DSS'ed PS dissipation
                 dpdiss(:,:) = &
                   ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%derived%psdiss_biharmonic(:,:)
                 dp_star(:,:,k) = dp_star(:,:,k) - rhs_viss*dt*nu_q*dpdiss(:,:)/elem(ie)%spheremp(:,:)
              endif
           enddo
           ! apply limiter to Q = Qtens / dp_star 
	   call limiter_optim_iter_full(Qtens(:,:,:),elem(ie)%spheremp(:,:),&
	      qmin(:,q,ie),qmax(:,q,ie),dp_star(:,:,:))

	endif

        ! apply mass matrix, overwrite np1 with solution:
        ! dont do this earlier, since we allow np1 to be the same as n0
        ! and we dont want to overwrite n0 until we are done using it
        do k=1,nlev
           elem(ie)%state%Qdp(:,:,k,q,np1) = elem(ie)%spheremp(:,:)*Qtens(:,:,k) 
        enddo

        if(limiter_option == 4)then
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	  ! sign-preserving limiter, applied after mass matrix
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	  call limiter2d_zero(elem(ie)%state%Qdp(:,:,:,q,np1),&
             elem(ie)%state%ps_v(:,:,np1),hvcoord) ! ps_v argument not used
        endif


     end do

     if(DSSopt==DSSno_var)then
	call edgeVpack(edgeAdv,elem(ie)%state%Qdp(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
     else
	call edgeVpack(edgeAdv_p1,elem(ie)%state%Qdp(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
	! also DSS extra field
#if (defined ELEMENT_OPENMP)
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
	call edgeVunpack(edgeAdv,elem(ie)%state%Qdp(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q)
#endif
	do q=1,qsize
           do k=1,nlev    !  Potential loop inversion (AAM)
	      elem(ie)%state%Qdp(:,:,k,q,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%Qdp(:,:,k,q,np1)
	   enddo
	end do
     else
	call edgeVunpack(edgeAdv_p1,elem(ie)%state%Qdp(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
        !$omp parallel do private(k,q)
#endif
        do q=1,qsize
           do k=1,nlev    !  Potential loop inversion (AAM)
	      elem(ie)%state%Qdp(:,:,k,q,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%Qdp(:,:,k,q,np1)
           enddo
	end do
	call edgeVunpack(edgeAdv_p1,DSSvar(:,:,1:nlev),nlev,qsize*nlev,elem(ie)%desc)
        
	do k=1,nlev
           DSSvar(:,:,k)=DSSvar(:,:,k)*elem(ie)%rspheremp(:,:)
	enddo
        
     endif
  end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
  call t_stopf('euler_step')
  
  end subroutine euler_step

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine euler_step_dg(np1,n0,dt,elem,hvcoord,hybrid,deriv,nets,nete,&
      compute_diagnostics,DSSopt,rhs_multiplier)
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
  integer :: np1,nm1,n0,nets,nete,DSSopt,rhs_multiplier
  real (kind=real_kind), intent(in)  :: dt
  logical  :: compute_diagnostics

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
  if (limiter_option == 8) then
     call abortmp('limiter_opiton=8 not supported for dg advection')
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
	call edgeVpack(edgeAdv,elem(ie)%state%Qdp(:,:,:,:,n0),nlev*qsize,0,elem(ie)%desc)
     else
	call edgeVpack(edgeAdv_p1,elem(ie)%state%Qdp(:,:,:,:,n0),nlev*qsize,0,elem(ie)%desc)
	! also DSS extra field
#if (defined ELEMENT_OPENMP)
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
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q)
#endif
	call edgeVunpack(edgeAdv_p1,DSSvar(:,:,1:nlev),nlev,qsize*nlev,elem(ie)%desc)
	do k=1,nlev
	  DSSvar(:,:,k)=DSSvar(:,:,k)*elem(ie)%rspheremp(:,:)
	enddo
     endif

     ! compute flux and advection term
     do k=1,nlev
        dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - &
             rhs_multiplier*dt*elem(ie)%derived%divdp_proj(:,:,k) 
        Vstar(:,:,1,k) = elem(ie)%derived%vn0(:,:,1,k)/dp(:,:,k)
        Vstar(:,:,2,k) = elem(ie)%derived%vn0(:,:,2,k)/dp(:,:,k)
     enddo

     do q=1,qsize
        do k=1,nlev
           vtemp(:,:,1)=elem(ie)%state%Qdp(:,:,k,q,n0)*Vstar(:,:,1,k)
           vtemp(:,:,2)=elem(ie)%state%Qdp(:,:,k,q,n0)*Vstar(:,:,2,k)
           
           divdp = divergence_sphere_wk(vtemp,deriv,elem(ie)) + &
                edge_flux_u_cg( Vstar(:,:,:,k), elem(ie)%state%Qdp(:,:,k,q,n0),qedges(:,:,k,q),&
                deriv, elem(ie), u_is_contra=.false.)

           ! advance in time. GLL quadrature, cardinal function basis, under-integrated.  
           ! local mass matrix is diagonal, with entries elem(ie)%spheremp(),
           ! so we divide through by elem(ie)%spheremp().
           elem(ie)%state%Qdp(:,:,k,q,np1)=elem(ie)%state%Qdp(:,:,k,q,n0) - dt*divdp/elem(ie)%spheremp
           
           if (npdg<np) then
              ! modal timestep, with exact integration.  using prognostic variable: p*metdet
              ! local mass matrix is diagonal assuming npdg<np so that GLL quadrature is exact)
              ! (note: GLL/modal conversion comutes with time-stepping)
              
              ! compute modal coefficients of p*metdet
              ! (spherical inner-product of Legendre polynomial and p)
              pshat = gll_to_dgmodal(elem(ie)%state%Qdp(:,:,k,q,np1)*elem(ie)%metdet(:,:),deriv)

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
              elem(ie)%state%Qdp(:,:,k,q,np1)=divdp/elem(ie)%metdet(:,:)
           endif
        enddo
        if(limiter_option == 4)then
           ! reuse CG limiter, which wants Qdp*spheremp:
           do k=1,nlev
              elem(ie)%state%Qdp(:,:,k,q,np1)=elem(ie)%state%Qdp(:,:,k,q,np1)*elem(ie)%spheremp(:,:)
           enddo
           call limiter2d_zero(elem(ie)%state%Qdp(:,:,:,q,np1),&
                elem(ie)%state%ps_v(:,:,np1),hvcoord) ! ps_v argument not used
           do k=1,nlev
              elem(ie)%state%Qdp(:,:,k,q,np1)=elem(ie)%state%Qdp(:,:,k,q,np1)/elem(ie)%spheremp(:,:)
           enddo
        endif
     end do
  end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
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

    use kinds, only : real_kind
    use dimensions_mod, only : np, np, nlev

    real (kind=real_kind), dimension(nlev), intent(inout)   :: minp, maxp
    real (kind=real_kind), dimension(np,np,nlev), intent(inout)   :: ptens
    real (kind=real_kind), dimension(np,np,nlev), intent(in), optional  :: dpmass
    real (kind=real_kind), dimension(np,np), intent(in)   :: sphweights
 
    real (kind=real_kind), dimension(np,np,nlev) :: weights
    real (kind=real_kind), dimension(np,np) :: ptens_mass
    integer  k1, k, i, j, iter, i1, i2
    integer :: whois_neg(np*np), whois_pos(np*np), neg_counter, pos_counter

    real (kind=real_kind) :: addmass, weightssum, mass
    real (kind=real_kind) :: x(np*np),c(np*np)

    real (kind=real_kind) :: al_neg(np*np), al_pos(np*np), howmuch
    real (kind=real_kind) :: tol_limiter=1e-15

    integer, parameter :: maxiter = 5

    do k=1,nlev
       weights(:,:,k)=sphweights(:,:)*dpmass(:,:,k)
       ptens(:,:,k)=ptens(:,:,k)/dpmass(:,:,k)
    enddo

    do k=1,nlev

      k1=1
      do i=1,np
	do j=1,np
	  c(k1)=weights(i,j,k)
	  x(k1)=ptens(i,j,k)
	  k1=k1+1
	enddo
      enddo

      mass=sum(c*x)

      ! relax constraints to ensure limiter has a solution:
      ! This is only needed if runnign with the SSP CFL>1 or 
      ! due to roundoff errors
      if( (mass/sum(c))<minp(k) ) then
         minp(k)=mass/sum(c)
      endif
      if( (mass/sum(c))>maxp(k) ) then
         maxp(k)=mass/sum(c)
      endif

      addmass=0.0d0
      pos_counter=0;
      neg_counter=0;
      
      ! apply constraints, compute change in mass caused by constraints 
      do k1=1,np*np
         if((x(k1)>=maxp(k)) ) then
	    addmass=addmass+(x(k1)-maxp(k))*c(k1)
	    x(k1)=maxp(k)
	    whois_pos(k1)=-1
         else
            pos_counter=pos_counter+1;
	    whois_pos(pos_counter)=k1;
         endif
         if((x(k1)<=minp(k)) ) then
	    addmass=addmass-(minp(k)-x(k1))*c(k1)
	    x(k1)=minp(k)
	    whois_neg(k1)=-1
         else
	    neg_counter=neg_counter+1;
	    whois_neg(neg_counter)=k1;
         endif
      enddo
      
      ! iterate to find field that satifies constraints and is l2-norm closest to original 
      weightssum=0.0d0
      if(addmass>0)then
         do i2=1,maxIter
            weightssum=0.0
            do k1=1,pos_counter
               i1=whois_pos(k1)
               weightssum=weightssum+c(i1)
               al_pos(i1)=maxp(k)-x(i1)
            enddo
            
            if((pos_counter>0).and.(addmass>tol_limiter*abs(mass)))then
               do k1=1,pos_counter
                  i1=whois_pos(k1)
                  howmuch=addmass/weightssum
                  if(howmuch>al_pos(i1))then
                     howmuch=al_pos(i1)
                     whois_pos(k1)=-1
                  endif
                  addmass=addmass-howmuch*c(i1)
                  weightssum=weightssum-c(i1)
                  x(i1)=x(i1)+howmuch
               enddo
               !now sort whois_pos and get a new number for pos_counter
               !here neg_counter and whois_neg serve as temp vars
               neg_counter=pos_counter
               whois_neg=whois_pos
               whois_pos=-1
               pos_counter=0
               do k1=1,neg_counter
                  if(whois_neg(k1).ne.-1)then
                     pos_counter=pos_counter+1
                     whois_pos(pos_counter)=whois_neg(k1)
                  endif
               enddo
            else
               exit
            endif
         enddo
      else
         do i2=1,maxIter
            weightssum=0.0
            do k1=1,neg_counter
               i1=whois_neg(k1)
               weightssum=weightssum+c(i1)
               al_neg(i1)=x(i1)-minp(k)
            enddo
            
            if((neg_counter>0).and.((-addmass)>tol_limiter*abs(mass)))then
               
               do k1=1,neg_counter
                  i1=whois_neg(k1)
                  howmuch=-addmass/weightssum
                  if(howmuch>al_neg(i1))then
                     howmuch=al_neg(i1)
                     whois_neg(k1)=-1
                  endif
                  addmass=addmass+howmuch*c(i1)
                  weightssum=weightssum-c(i1)
                  x(i1)=x(i1)-howmuch
               enddo
               !now sort whois_pos and get a new number for pos_counter
               !here pos_counter and whois_pos serve as temp vars
               pos_counter=neg_counter
               whois_pos=whois_neg
               whois_neg=-1
               neg_counter=0
               do k1=1,pos_counter
                  if(whois_pos(k1).ne.-1)then
                     neg_counter=neg_counter+1
                     whois_neg(neg_counter)=whois_pos(k1)
                  endif
               enddo
            else
               exit
            endif
            
         enddo
         
      endif
      
      k1=1
      do i=1,np
         do j=1,np
            ptens(i,j,k)=x(k1)
            k1=k1+1
         enddo
      enddo
      
   enddo
   
   do k=1,nlev
      ptens(:,:,k)=ptens(:,:,k)*dpmass(:,:,k)
   enddo

  end subroutine limiter_optim_iter_full


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter2d_minmax(Q,dp,hvcoord,spheremp,qmin,qmax)
!
! mass conserving limiter (2D only).  to be called just before DSS
!
! in pure 2D advection, the element mass will not be negative before DSS
! this routine will redistribute to remove negative values (conservative)
!
! if used in 3D, should be applied with 2D/vertical split advection
! 
! call with Qdp and assocated dp
!
!
  implicit none
  real (kind=real_kind), intent(inout)  :: Q(np,np,nlev)
  real (kind=real_kind), intent(in)  :: spheremp(np,np)
  real (kind=real_kind), intent(in) ::  dp(np,np,nlev)
  type (hvcoord_t)                 :: hvcoord

  ! local
  integer i,j,k
  real (kind=real_kind) :: mass,mass_new,area,qmin(nlev),qmax(nlev),mass2


#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,mass,area,mass2,mass_new,i,j)
#endif
  do k=1,nlev
     mass=sum( Q(:,:,k)*spheremp(:,:) )
     area=sum( dp(:,:,k)*spheremp(:,:) )

     Q(:,:,k)=Q(:,:,k)/dp(:,:,k)  ! % convert to concentration


!     if (mass>0) print *,k,mass/area,qmin(k),qmax(k)
     ! max limiter
     if ( maxval(Q(:,:,k)) > qmax(k) ) then
        
        Q(:,:,k)=qmax(k)-Q(:,:,k)      ! some of these will be negative
        mass2 = area*qmax(k) - mass

        if (mass2 < 0) Q(:,:,k)=-Q(:,:,k) 
        mass_new=0
        do j=1,np	
        do i=1,np
           if (Q(i,j,k)<0) then
              Q(i,j,k)=0
           else
              mass_new = mass_new + Q(i,j,k)*dp(i,j,k)*spheremp(i,j)
           endif
        enddo
        enddo
     
        ! now scale the all positive values to restore mass
        if (mass_new>0) Q(:,:,k) = Q(:,:,k)*abs(mass2)/mass_new
        if (mass2 < 0) Q(:,:,k)=-Q(:,:,k) 
        
        Q(:,:,k)=qmax(k)-Q(:,:,k)
     endif


     ! min limiter
     if ( minval(Q(:,:,k)) < qmin(k) ) then
        Q(:,:,k)=Q(:,:,k)-qmin(k)
        mass2 = mass - area*qmin(k)
        ! negative mass.  so reduce all postive values to zero 
        ! then increase negative values as much as possible
        if (mass2 < 0) Q(:,:,k)=-Q(:,:,k) 
        mass_new=0
        do j=1,np	
           do i=1,np
              if (Q(i,j,k)<0) then
                 Q(i,j,k)=0
              else
                 mass_new = mass_new + Q(i,j,k)*dp(i,j,k)*spheremp(i,j)
              endif
           enddo
        enddo
        
        ! now scale the all positive values to restore mass
        if (mass_new>0) Q(:,:,k) = Q(:,:,k)*abs(mass2)/mass_new
        if (mass2 < 0) Q(:,:,k)=-Q(:,:,k) 

        Q(:,:,k)=Q(:,:,k)+qmin(k)
     endif

     Q(:,:,k)=Q(:,:,k)*dp(:,:,k)

  enddo

end subroutine 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter2d_zero(Q,ps,hvcoord)
!
! mass conserving zero limiter (2D only).  to be called just before DSS
!
! this routine is called inside a DSS loop, and so Q had already
! been multiplied by the mass matrix.  Thus dont include the mass
! matrix when computing the mass = integral of Q over the element
!
! ps is only used when advecting Q instead of Qdp
! so ps should be at one timelevel behind Q
!
  implicit none
  real (kind=real_kind), intent(inout)  :: Q(np,np,nlev)
  real (kind=real_kind), intent(in)  :: ps(np,np)
  type (hvcoord_t)                 :: hvcoord

  ! local
  real (kind=real_kind) ::  dp(np,np)
  integer i,j,k
  real (kind=real_kind) :: mass,mass_new,ml


  do k=nlev,1,-1

     mass=0
     do j=1,np	
     do i=1,np
        !ml = Q(i,j,k)*dp(i,j)*spheremp(i,j)  ! see above
        if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
           ml = Q(i,j,k)
        else
           dp(i,j) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*ps(i,j)
           ml = Q(i,j,k)*dp(i,j)
        endif
        mass = mass + ml
     enddo
     enddo

     ! negative mass.  so reduce all postive values to zero 
     ! then increase negative values as much as possible
     if (mass < 0) Q(:,:,k)=-Q(:,:,k) 
     mass_new=0
     do j=1,np	
     do i=1,np
        if (Q(i,j,k)<0) then
           Q(i,j,k)=0
        else
           if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
              ml = Q(i,j,k)
           else
              ml = Q(i,j,k)*dp(i,j)
           endif
           mass_new = mass_new + ml
        endif
     enddo
     enddo

     ! now scale the all positive values to restore mass
     if ( mass_new > 0) Q(:,:,k) = Q(:,:,k)*abs(mass)/mass_new
     if (mass < 0) Q(:,:,k)=-Q(:,:,k) 
  enddo

end subroutine 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter3d_noncon(Q,ps,hvcoord,mass_added)
!
! local element negative mass "fixer"  
! Trucate negative values, then use CAM style mass fixer to restore
! original mass.  Mass fixer is local to the element.  Will break C0
! do DSS needs to be done afterwards.
! 
! this routine is called inside a DSS loop, and so Q had already
! been multiplied by the mass matrix.  Thus dont include the mass
! matrix when computing the mass = integral of Q over the element
!
! ps is only used when advecting Q instead of Qdp
! so ps should be at one timelevel behind Q
!
  implicit none
  real (kind=real_kind), intent(inout)  :: Q(np,np,nlev)
  real (kind=real_kind), intent(in)  :: ps(np,np)
  real (kind=real_kind), intent(out)  :: mass_added
  type (hvcoord_t)                 :: hvcoord

  ! local
  real (kind=real_kind) ::  dp(np,np)
  integer i,j,k
  real (kind=real_kind) :: mass,mass_new,ml

  mass_added=0
  mass=0
  mass_new=0


  if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
  do k=nlev,1,-1
     do j=1,np	
        do i=1,np
!           ml = Q(i,j,k)*dp(i,j)*spheremp(i,j)  ! see above
           mass = mass + Q(i,j,k)
           if (Q(i,j,k)<0) then
              Q(i,j,k)=0
           else
              mass_new = mass_new + Q(i,j,k)
           endif
        enddo
     enddo
  enddo
  else
  do k=nlev,1,-1
     dp(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*ps(:,:)
     do j=1,np	
        do i=1,np
!           ml = Q(i,j,k)*dp(i,j)*spheremp(i,j)  ! see above
           ml = Q(i,j,k)*dp(i,j)
           mass = mass + ml
           
           if (Q(i,j,k)<0) then
              Q(i,j,k)=0
           else
              mass_new = mass_new + ml
           endif
        enddo
     enddo
  enddo
  endif
 
  ! now scale the all positive values to restore mass
  if ( mass > 0) then
     ! rescale positive values to restore original mass
     Q(:,:,:) = Q(:,:,:)*mass/mass_new
  else
     ! mass was negative.  set all values to zero
     Q(:,:,:) = 0
     mass_added = mass_added -mass
  endif
  
  
  end subroutine 
  
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine preq_impsysQ(elemin,hvcoord,np1,n0,nm1,Q,Q_vadv)
    use perf_mod, only : t_startf, t_stopf                    ! _EXTERNAL
    implicit none

    type (element_t), intent(in)     :: elemin
    type (hvcoord_t)                 :: hvcoord
    integer, intent(in)              :: np1
    integer, intent(in)              :: n0
    integer, intent(in)              :: nm1

    real(kind=real_kind), intent(in),target, dimension(np,np,nlev)   :: Q         
    real(kind=real_kind), intent(inout), dimension(np,np,nlev)       :: Q_vadv    

    ! ==============================
    ! Local Variables
    ! ==============================
    
    real(kind=real_kind), dimension(nlevp)      :: hyai
    real(kind=real_kind), dimension(nlevp)      :: hybi
    real(kind=real_kind), dimension(nlev)       :: hyam
    real(kind=real_kind), dimension(nlev)       :: hybm
    real(kind=real_kind), dimension(nlev)       :: hybd 
    real(kind=real_kind), dimension(np,np,nlev) :: div

    real(kind=real_kind), dimension(np,np)      :: ps                ! surface pressure
    real(kind=real_kind), dimension(np,np)      :: rps               ! 1/ps
    real(kind=real_kind), dimension(np,np,nlev) :: rpdel             ! 1./pdel

    real(kind=real_kind), dimension(np,np,nlevp)  :: eta_dot_dp_deta       ! eta dot * dp/deta at time level n
    real(kind=real_kind), dimension(np,np,nlev)   :: vgrad_ps              ! ps*(v.grad(lnps))

    real(kind=real_kind) :: pint(np,np,nlevp)
    real(kind=real_kind) :: pdel(np,np,nlev)
    real(kind=real_kind) :: pmid(np,np,nlev)

    real(kind=real_kind) :: v1, v2, vcon1, vcon2

    integer :: i,j,k,l

    call t_startf('preq_impsysQ')




    hyai   = hvcoord%hyai
    hybi   = hvcoord%hybi
    hyam   = hvcoord%hyam
    hybm   = hvcoord%hybm
    hybd   = hvcoord%hybd

    ps(:,:)  = EXP(elemin%state%lnps(:,:,n0))
    rps(:,:) = 1.0_real_kind/ps(:,:)

    call preq_pressure(hvcoord%ps0,  ps,&
         hyai, hybi, hyam, hybm,&
         pint, pmid, pdel)

    rpdel = 1.0_real_kind/pdel

!   v should be contravariant for explicit time step

    do k=1,nlev
       do j=1,np
          do i=1,np
             v1 = elemin%state%v(i,j,1,k,n0)
             v2 = elemin%state%v(i,j,2,k,n0)

!            vcon1 = elemin%metinv(1,1,i,j)*v1 + elemin%metinv(1,2,i,j)*v2
!            vcon2 = elemin%metinv(2,1,i,j)*v1 + elemin%metinv(2,2,i,j)*v2

             vcon1 = elemin%Dinv(1,1,i,j)*v1 + elemin%Dinv(1,2,i,j)*v2
             vcon2 = elemin%Dinv(2,1,i,j)*v1 + elemin%Dinv(2,2,i,j)*v2

             vgrad_ps(i,j,k) = ps(i,j)*(vcon1*elemin%derived%grad_lnps(i,j,1) + vcon2*elemin%derived%grad_lnps(i,j,2))
          end do
       end do
    end do

    eta_dot_dp_deta(:,:,1) = 0.0_real_kind

    do k=1,nlev
       div(:,:,k) = elemin%derived%div(:,:,k,n0)
       do j=1,np
          do i=1,np
             eta_dot_dp_deta(i,j,k+1) = eta_dot_dp_deta(i,j,k) + vgrad_ps(i,j,k)*hybd(k) + div(i,j,k)*pdel(i,j,k)
          end do
       end do
    end do

    do k=1,nlev-1
       do j=1,np
          do i=1,np
             eta_dot_dp_deta(i,j,k+1) = hybi(k+1)*eta_dot_dp_deta(i,j,nlev+1) - eta_dot_dp_deta(i,j,k+1)
          end do
       end do
    end do

    eta_dot_dp_deta(:,:,nlev+1) = 0.0_real_kind

    call preq_vertadvQ(Q(:,:,:),eta_dot_dp_deta,rpdel,Q_vadv)
    call t_stopf('preq_impsysQ')

  end subroutine preq_impsysQ



  subroutine preq_vertadvQ(Q,eta_dot_dp_deta, rpdel, Q_vadv)
    use perf_mod, only : t_startf, t_stopf                   ! _EXTERNAL
    implicit none

    real (kind=real_kind), intent(in)  :: Q(np,np,nlev)
    real (kind=real_kind), intent(in)  :: eta_dot_dp_deta(np,np,nlevp)
    real (kind=real_kind), intent(in)  :: rpdel(np,np,nlev)

    real (kind=real_kind), intent(inout) :: Q_vadv(np,np,nlev)

    ! ========================
    ! Local Variables
    ! ========================

    integer                            :: i,j,k,l
    real (kind=real_kind)              :: facp, facm

    ! ===========================================================
    ! Compute vertical advection of Q from eq. (3.b.1)
    !
    ! k = 1 case:
    ! ===========================================================
    call t_startf('preq_vertadvQ')

    k=1
    do j=1,np
       do i=1,np 
          facp  = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k+1)             
          Q_vadv(i,j,k) = facp*(Q(i,j,k+1)- Q(i,j,k))
       end do
    end do

    ! ===========================================================
    ! vertical advection
    !
    ! 1 < k < nlev case:
    ! ===========================================================

    do k=2,nlev-1
       do j=1,np
          do i=1,np
             facp = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k+1)
             facm = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k)                
             Q_vadv(i,j,k) = facp*(Q(i,j,k+1)- Q(i,j,k)) + facm*(Q(i,j,k)- Q(i,j,k-1))
          enddo
       end do
    end do

    ! ===========================================================
    ! vertical advection
    !
    ! k = nlev case:
    ! ===========================================================

    k=nlev
    do j=1,np
       do i=1,np
          facm = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k)             
          Q_vadv(i,j,k) = facm*(Q(i,j,k)- Q(i,j,k-1))
       enddo
    end do
    call t_stopf('preq_vertadvQ')

  end subroutine preq_vertadvQ




  subroutine advance_hypervis_scalar_lf(edgeAdv,elem,hvcoord,hybrid,deriv,nt,n0,nets,nete,dt2)
  !
  !  hyperviscosity operator used by leapfrog scheme
  !  take one timestep of:  
  !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use kinds, only : real_kind
  use dimensions_mod, only : np, nlev
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use perf_mod, only: t_startf, t_stopf                  ! _EXTERNAL
  implicit none

  type (hvcoord_t), intent(in)      :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edgeAdv
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: dt2
  integer :: nets,nete,nt,n0

  
  ! local
  integer :: k,kptr,i,j,ie,ic,q
  real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: Qtens

! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
!       data is incorrect (offset by a few numbers actually)
!       removed for now.  
!  real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp

  real (kind=real_kind), dimension(np,np) :: lap_p
  real (kind=real_kind) :: v1,v2,dt,nu_scale,dp,dp0
  integer :: density_scaling = 0
  real(kind=real_kind), dimension(:,:), pointer :: viscosity => NULL()
  
  if (nu_q == 0) return;
  call t_barrierf('sync_advance_hypervisc_scalar_lf', hybrid%par%comm)
  call t_startf('advance_hypervis_scalar_lf')

  if (tracer_advection_formulation == TRACERADV_UGRADQ ) then
     ! conservative advection using non-conservative form of the equations
     ! in this case, we use a density scaled viscosity coefficieint: 
     density_scaling = 1
  endif
  
  
  dt=dt2/hypervis_subcycle_q
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     do ic=1,hypervis_subcycle_q
        do ie=nets,nete
           viscosity     => elem(ie)%variable_hyperviscosity
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q,lap_p)
#endif
	       do q=1,qsize
              do k=1,nlev    !  Potential loop inversion (AAM)
                 lap_p=laplace_sphere_wk(elem(ie)%state%Q(:,:,k,q,nt),deriv,elem(ie),viscosity)
                 ! advace in time.  (note: DSS commutes with time stepping, so we
                 ! can time advance and then DSS.
                 elem(ie)%state%Q(:,:,k,q,nt)=elem(ie)%state%Q(:,:,k,q,nt)*elem(ie)%spheremp(:,:)  +  dt*nu_q*lap_p(:,:) 
              enddo
           enddo
           call edgeVpack(edgeAdv, elem(ie)%state%Q(:,:,:,:,nt),nlev*qsize,0,elem(ie)%desc)
        enddo
           
        call bndry_exchangeV(hybrid,edgeAdv)
        
        do ie=nets,nete
           call edgeVunpack(edgeAdv, elem(ie)%state%Q(:,:,:,:,nt),nlev*qsize,0,elem(ie)%desc)
              !rspheremp     => elem(ie)%rspheremp
              ! apply inverse mass matrix
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q,i,j)
#endif
	       do q=1,qsize
              do k=1,nlev    !  Potential loop inversion (AAM)
              do j=1,np
              do i=1,np             
                 elem(ie)%state%Q(i,j,k,q,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%Q(i,j,k,q,nt)
              enddo
              enddo
              enddo
           enddo
        enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo  ! subcycle
  endif
        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle_q
        do ie=nets,nete
           if (density_scaling==1) then
              ! state%Q really is Q !
              Qtens(:,:,:,:,ie)=elem(ie)%state%Q(:,:,:,1:qsize,nt)
           else
              ! state%Q is really Qdp.  but we only apply diffusion on Q
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,dp)
#endif
              do k=1,nlev
              do j=1,np
              do i=1,np
                 ! note: use ps(t+1) to get exact consistency
                 dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,nt)
                 Qtens(i,j,k,:,ie)=elem(ie)%state%Q(i,j,k,1:qsize,nt)/dp
              enddo
              enddo
              enddo
           endif
        enddo
        ! compute biharmonic operator. Qtens = input and output 
        call biharmonic_wk_scalar(elem,Qtens,deriv,edgeAdv,hybrid,nets,nete)
        do ie=nets,nete
           !spheremp     => elem(ie)%spheremp
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k,i,j,dp0,dp,nu_scale)
#endif
           do q=1,qsize
           do k=1,nlev
           do j=1,np
           do i=1,np

              dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                   ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0

              if (density_scaling==1) then
                 ! advection Q.  For conservation:     dp0/dp DIFF(Q)
                 ! scale velosity by 1/rho (normalized to be O(1))
                 ! dp = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                 ! NOTE: use ps(t) to get exact mass conservation
                 dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,n0)
                 nu_scale = dp0/dp
                 elem(ie)%state%Q(i,j,k,q,nt)  =  elem(ie)%state%Q(i,j,k,q,nt)*elem(ie)%spheremp(i,j) &
                              -dt*nu_q*nu_scale*Qtens(i,j,k,q,ie)
              else
                 ! advection Qdp.  For mass advection consistency:
                 ! DIFF( Qdp) ~   dp0 DIFF (Q)  =  dp0 DIFF ( Qdp/dp )  
!                 elem(ie)%state%Q(i,j,k,q,nt)  =  elem(ie)%state%Q(i,j,k,q,nt)*elem(ie)%spheremp(i,j) &
!                        -dt*nu_q*Qtens(i,j,k,q,ie)
                 elem(ie)%state%Q(i,j,k,q,nt)  =  elem(ie)%state%Q(i,j,k,q,nt)*elem(ie)%spheremp(i,j) &
                        -dt*nu_q*dp0*Qtens(i,j,k,q,ie)
              endif
           enddo
           enddo
           enddo

           ! smooth some of the negativities introduced by diffusion:
           ! note: ps_v not used if advecting Qdp
           call limiter2d_zero(elem(ie)%state%Q(:,:,:,q,nt),&
                elem(ie)%state%ps_v(:,:,n0),hvcoord)

           enddo
           call edgeVpack(edgeAdv,elem(ie)%state%Q(:,:,:,:,nt),qsize*nlev,0,elem(ie)%desc)
        enddo

        call bndry_exchangeV(hybrid,edgeAdv)

        do ie=nets,nete
        call edgeVunpack(edgeAdv, elem(ie)%state%Q(:,:,:,:,nt), qsize*nlev, 0, elem(ie)%desc)
        !rspheremp     => elem(ie)%rspheremp
        ! apply inverse mass matrix
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q,i,j)
#endif
	       do q=1,qsize
           do k=1,nlev    !  Potential loop inversion (AAM)
           do j=1,np
           do i=1,np
              elem(ie)%state%Q(i,j,k,q,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%Q(i,j,k,q,nt)
           enddo
           enddo
           enddo
           enddo
        enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo
  endif

  call t_stopf('advance_hypervis_scalar_lf')
  
  end subroutine advance_hypervis_scalar_lf



  subroutine advance_hypervis_scalar(edgeAdv,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2)
  !
  !  hyperviscsoity operator for foward-in-time scheme
  !  take one timestep of:  
  !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use kinds, only : real_kind
  use dimensions_mod, only : np, nlev
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use perf_mod, only: t_startf, t_stopf                          ! _EXTERNAL
  implicit none

  type (hvcoord_t), intent(in)      :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edgeAdv
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: dt2
  integer :: nets,nete,nt

  
  ! local
  integer :: k,kptr,i,j,ie,ic,q
  real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: Qtens
  real (kind=real_kind), dimension(np,np,nlev) :: dp

  real (kind=real_kind), dimension(nlev,qsize,nets:nete) :: min_neigh
  real (kind=real_kind), dimension(nlev,qsize,nets:nete) :: max_neigh

! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
!       data is incorrect (offset by a few numbers actually)
!       removed for now.  
!  real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp

  real (kind=real_kind), dimension(np,np) :: lap_p
  real (kind=real_kind) :: v1,v2,dt,dp0
  integer :: density_scaling = 0
  
  if (nu_q == 0) return;
  if (hypervis_order /= 2) return
  call t_barrierf('sync_advance_hypervis_scalar', hybrid%par%comm)
  call t_startf('advance_hypervis_scalar')
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dt=dt2/hypervis_subcycle_q


  do ic=1,hypervis_subcycle_q
     do ie=nets,nete
        ! Qtens = Q/dp   (apply hyperviscsoity to dp0 * Q, not Qdp)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,q)
#endif
        do k=1,nlev
           dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,nt)
           do q=1,qsize
              Qtens(:,:,k,q,ie)=elem(ie)%state%Qdp(:,:,k,q,nt)/dp(:,:,k)
           enddo
        enddo
     enddo

     ! compute biharmonic operator. Qtens = input and output 
     call biharmonic_wk_scalar(elem,Qtens,deriv,edgeAdv,hybrid,nets,nete)
     do ie=nets,nete
        !spheremp     => elem(ie)%spheremp
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k,i,j)
#endif
        do q=1,qsize
           do k=1,nlev
              dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                   ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
           do j=1,np
           do i=1,np
              ! advection Qdp.  For mass advection consistency:
              ! DIFF( Qdp) ~   dp0 DIFF (Q)  =  dp0 DIFF ( Qdp/dp )  
              elem(ie)%state%Qdp(i,j,k,q,nt)  =  elem(ie)%state%Qdp(i,j,k,q,nt)*elem(ie)%spheremp(i,j) &
                   -dt*nu_q*dp0*Qtens(i,j,k,q,ie)
           enddo
           enddo
           enddo

           ! smooth some of the negativities introduced by diffusion:
           call limiter2d_zero(elem(ie)%state%Qdp(:,:,:,q,nt),&
                elem(ie)%state%ps_v(:,:,nt),hvcoord)
        enddo
        call edgeVpack(edgeAdv,elem(ie)%state%Qdp(:,:,:,:,nt),qsize*nlev,0,elem(ie)%desc)
     enddo

     call bndry_exchangeV(hybrid,edgeAdv)
     
     do ie=nets,nete
        call edgeVunpack(edgeAdv, elem(ie)%state%Qdp(:,:,:,:,nt), qsize*nlev, 0, elem(ie)%desc)
        !rspheremp     => elem(ie)%rspheremp
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k)
#endif
        do q=1,qsize    
           ! apply inverse mass matrix
           do k=1,nlev
              elem(ie)%state%Qdp(:,:,k,q,nt)=elem(ie)%rspheremp(:,:)*elem(ie)%state%Qdp(:,:,k,q,nt)
           enddo
        enddo
     enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
  enddo
  call t_stopf('advance_hypervis_scalar')
  end subroutine advance_hypervis_scalar









  subroutine preq_vertadv_dpQ(Q,eta_dot_dp_deta, Q_vadv)
    implicit none
    
    real (kind=real_kind), intent(in)  :: Q(np,np,nlev)
    real (kind=real_kind), intent(in)  :: eta_dot_dp_deta(np,np,nlevp)
    real (kind=real_kind), intent(inout) :: Q_vadv(np,np,nlev)

    ! ========================
    ! Local Variables
    ! ========================

    integer                            :: i,j,k,l
    real (kind=real_kind)              :: qp,qm
    !
    ! compute  d(eta_dot_dp_deta Q)
    ! 
    !
    qp=0
    qm=0
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,qp,qm)
#endif
    do k=1,nlev
       do i=1,np
       do j=1,np
       ! at k=1,nlev, eta_dot_dp_deta is zero, so value of q does not matter
#undef QADV_UPWIND
#ifdef QADV_UPWIND
       ! UPWIND, 1st order:
          if (k/=nlev) then
             if (eta_dot_dp_deta(i,j,k+1) > 0 ) then
                qp=Q(i,j,k)
             else
                qp=Q(i,j,k+1)
             endif
          endif
          if (k/=1) then
             if (eta_dot_dp_deta(i,j,k) > 0 ) then
                qm=Q(i,j,k-1)
             else
                qm=Q(i,j,k)
             endif
          endif
#else
          ! CENTERED, 2nd order:
          if (k/=nlev) then
             qp=.5*(Q(i,j,k)+Q(i,j,k+1))
          endif
          if (k/=1) then
             qm=.5*(Q(i,j,k)+Q(i,j,k-1))
          endif
#endif
          Q_vadv(i,j,k)= eta_dot_dp_deta(i,j,k+1)*qp - eta_dot_dp_deta(i,j,k)*qm
       enddo
      enddo
    enddo

end subroutine preq_vertadv_dpQ
end module prim_advection_mod









 
