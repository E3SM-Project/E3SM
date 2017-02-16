#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_state_mod

  use kinds,            only: real_kind, iulog
  use dimensions_mod,   only: nlev, np, qsize_d, qsize, nelemd
  use parallel_mod,     only:  iam, ordered, parallel_t, syncmp
  use parallel_mod,     only: global_shared_buf, global_shared_sum
  use global_norms_mod, only: wrap_repro_sum
  use hybrid_mod,       only: hybrid_t
  use time_mod,         only: tstep, secpday, timelevel_t, TimeLevel_Qdp, time_at
  use control_mod,      only: integration, test_case,  moisture, &
                              qsplit, ftype, use_cpstar, rsplit,&
                              theta_hydrostatic_mode
  use hybvcoord_mod,    only: hvcoord_t
  use global_norms_mod, only: global_integral, linf_snorm, l1_snorm, l2_snorm
  use element_mod,      only: element_t
  use element_ops,      only: get_field, get_pnh_and_exner
  use viscosity_mod,    only: compute_zeta_C0
  use reduction_mod,    only: parallelmax,parallelmin
  use perf_mod,         only: t_startf, t_stopf
  use physical_constants, only : p0,Cp,g,Rgas

implicit none
private
  character(len=*), private, parameter :: massfname = "mass.out"

  integer, public :: naccum


  public :: prim_printstate
  public :: prim_printstate_init
  public :: prim_energy_halftimes
  public :: prim_diag_scalars

contains
!=======================================================================================================! 


  subroutine prim_printstate_init(par)
    type (parallel_t) :: par

    real (kind=real_kind) :: time
    integer               :: c0

    if (par%masterproc) then
       time=0.0D0
       c0  =0.0D0

#ifndef CAM
#ifndef _BGL
       open(unit=10,file=massfname,form="formatted",status="replace")
!       write(10,*)time,c0
       close(10)
#endif
#endif
    end if

  end subroutine prim_printstate_init
!=======================================================================================================! 

  subroutine prim_printstate(elem, tl,hybrid,hvcoord,nets,nete)

    use physical_constants,     only: dd_pi

    type(element_t),            intent(in) :: elem(:)
    type(TimeLevel_t),target,   intent(in) :: tl
    type(hybrid_t),             intent(in) :: hybrid
    type(hvcoord_t),            intent(in) :: hvcoord
    integer,                    intent(in) :: nets,nete

    ! Local variables...
    integer :: i,j,k,ie
    integer,parameter  :: type=ORDERED

    real (kind=real_kind)  :: Mass2,Mass
    real (kind=real_kind)  :: TOTE(4),KEner(4),PEner(4),IEner(4)
    real (kind=real_kind)  :: Qvar(qsize_d,4),Qmass(qsize_d,4),Q1mass(qsize_d)
    real (kind=real_kind),save  :: time0
    real (kind=real_kind),save  :: TOTE0=0,Qmass0(qsize_d)=0
    real (kind=real_kind)  :: I_div(nlev)
    real (kind=real_kind)  :: I_vort(nlev)

    real (kind=real_kind)  :: tmp(np,np,nets:nete)
    real (kind=real_kind),allocatable  :: tmp3d(:,:,:,:)
    real (kind=real_kind)  :: tmp1(nets:nete)
    real (kind=real_kind)  :: ps(np,np)
    real (kind=real_kind)  :: dp(np,np)
    real (kind=real_kind)  :: dphi(np,np,nlev-1)
    real (kind=real_kind)  :: tdiag(np,np,nlev)
    !    real (kind=real_kind)  :: E(np,np)

    real (kind=real_kind) :: umin_local(nets:nete), umax_local(nets:nete), usum_local(nets:nete), &
         vmin_local(nets:nete), vmax_local(nets:nete), vsum_local(nets:nete), &
         tmin_local(nets:nete), tmax_local(nets:nete), tsum_local(nets:nete), &
         psmin_local(nets:nete),psmax_local(nets:nete),pssum_local(nets:nete), &
         fumin_local(nets:nete),fumax_local(nets:nete),fusum_local(nets:nete), &
         fvmin_local(nets:nete),fvmax_local(nets:nete),fvsum_local(nets:nete), &
         ftmin_local(nets:nete),ftmax_local(nets:nete),ftsum_local(nets:nete), &
         fqmin_local(nets:nete),fqmax_local(nets:nete),fqsum_local(nets:nete), &
         wmin_local(nets:nete),wmax_local(nets:nete),wsum_local(nets:nete),&
         phimin_local(nets:nete),phimax_local(nets:nete),phisum_local(nets:nete),&
         dpmin_local(nets:nete), dpmax_local(nets:nete), dpsum_local(nets:nete)


    real (kind=real_kind) :: umin_p, vmin_p, tmin_p, qvmin_p(qsize_d),&
         psmin_p, dpmin_p


    real (kind=real_kind) :: umax_p, vmax_p, tmax_p, qvmax_p(qsize_d),&
         psmax_p, dpmax_p

    real (kind=real_kind) :: usum_p, vsum_p, tsum_p, qvsum_p(qsize_d),&
         pssum_p, dpsum_p

    real (kind=real_kind) :: fusum_p, fvsum_p, ftsum_p, fqsum_p
    real (kind=real_kind) :: fumin_p, fvmin_p, ftmin_p, fqmin_p
    real (kind=real_kind) :: fumax_p, fvmax_p, ftmax_p, fqmax_p
    real (kind=real_kind) :: wmax_p, wmin_p, wsum_p
    real (kind=real_kind) :: phimax_p, phimin_p, phisum_p


    real(kind=real_kind) :: vsum_t(1), relvort
    real(kind=real_kind) :: v1, v2, vco(np,np,2,nlev)

    real (kind=real_kind) :: time, time2,time1, scale, dt, dt_split
    real (kind=real_kind) :: KEvert,IEvert,PEvert,T1,T2,S1,S2,P1,P2
    real (kind=real_kind) :: KEhorz,KEhorz2,PEhorz,IEhorz,KEH1,KEH2
    real (kind=real_kind) :: ddt_tot,ddt_diss
    integer               :: n0, nm1, np1, n0q
    integer               :: npts,n,q
    
    call t_startf('prim_printstate')
    if (hybrid%masterthread) then 
       if (Time_at(tl%nstep) <= 3600) then
          write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)," [s]"
       else if (Time_at(tl%nstep) <= 24*3600) then
          write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/3600," [h]"
       else
          write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
       endif
    end if

    TOTE     = 0
    KEner    = 0
    PEner    = 0
    IEner    = 0
    ! dynamics timelevels
    n0=tl%n0
    nm1=tl%nm1
    np1=tl%np1
    call TimeLevel_Qdp( tl, qsplit, n0q) !get n0 level into t2_qdp 


    dt=tstep*qsplit
    if (rsplit>0) dt = tstep*qsplit*rsplit  ! vertical REMAP timestep 
    !
    !   dynamics variables in n0 are at time =  'time' 
    !
    !   Diagnostics computed a 4 different levels during one compute REMAP step
    ! in RK code:
    !   E(:,:,2)-E(:,:,1)   change due to dynamics step  from time-dt to time
    !   E(:,:,3)-E(:,:,2)   change due to energy fixer   
    !   E(:,:,1)-E(:,:,4)   impact of forcing
    !
    ! Dissipation rates were computed during the first dynamics timstep, and represent
    ! the change going from 'time-dt' to 'time-dt+tstep' (one dynamics step)
    !
    time=tl%nstep*tstep
    time2 = time 
    time1 = time - dt


    vsum_t(1) = 0.0D0

    ! npts = np
    npts=SIZE(elem(1)%state%ps_v(:,:,n0),1)

    do q=1,qsize
       do ie=nets,nete
          tmp1(ie) = MINVAL(elem(ie)%state%Q(:,:,:,q))
       enddo
       qvmin_p(q) = ParallelMin(tmp1,hybrid)
       do ie=nets,nete
          tmp1(ie) = MAXVAL(elem(ie)%state%Q(:,:,:,q))
       enddo

       qvmax_p(q) = ParallelMax(tmp1,hybrid)

       do ie=nets,nete
          global_shared_buf(ie,1) = SUM(elem(ie)%state%Q(:,:,:,q))
       enddo
       call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
       qvsum_p(q) = global_shared_sum(1)
    enddo

    !
    do ie=nets,nete
       if (theta_hydrostatic_mode) then
          ! show min/max/num of temperature as a diagnostic
          call get_field(elem(ie),'temperature',tdiag,hvcoord,n0,n0q)
       else
          ! show min/max/num of dpnh / dp as a diagnostic
          call get_field(elem(ie),'dpnh_dp',tdiag,hvcoord,n0,n0q)
       endif

       tmp(:,:,ie)=elem(ie)%state%ps_v(:,:,n0)
       do k=1,nlev-1
          dphi(:,:,k)=elem(ie)%state%phi(:,:,k,n0)-elem(ie)%state%phi(:,:,k+1,n0)
       enddo
       !======================================================  
       umax_local(ie)    = MAXVAL(elem(ie)%state%v(:,:,1,:,n0))
       vmax_local(ie)    = MAXVAL(elem(ie)%state%v(:,:,2,:,n0))
       wmax_local(ie)    = MAXVAL(elem(ie)%state%w(:,:,:,n0))
       phimax_local(ie)  = MAXVAL(dphi(:,:,:))

       fumax_local(ie)    = MAXVAL(elem(ie)%derived%FM(:,:,1,:))
       fvmax_local(ie)    = MAXVAL(elem(ie)%derived%FM(:,:,2,:))

       tmax_local(ie)    = MAXVAL(tdiag)

       if (rsplit>0) &
            dpmax_local(ie)    = MAXVAL(elem(ie)%state%dp3d(:,:,:,n0))

       psmax_local(ie) = MAXVAL(tmp(:,:,ie))
       ftmax_local(ie)    = MAXVAL(elem(ie)%derived%FT(:,:,:))
       fqmax_local(ie)    = MAXVAL(elem(ie)%derived%FQ(:,:,:,1))
       !======================================================

       umin_local(ie)    = MINVAL(elem(ie)%state%v(:,:,1,:,n0))
       vmin_local(ie)    = MINVAL(elem(ie)%state%v(:,:,2,:,n0))
       Wmin_local(ie)    = MINVAL(elem(ie)%state%w(:,:,:,n0))
       phimin_local(ie)  = MINVAL(dphi)

       Fumin_local(ie)    = MINVAL(elem(ie)%derived%FM(:,:,1,:))
       Fvmin_local(ie)    = MINVAL(elem(ie)%derived%FM(:,:,2,:))

       tmin_local(ie)    = MINVAL(tdiag)

       if (rsplit>0) &
            dpmin_local(ie)    = MINVAL(elem(ie)%state%dp3d(:,:,:,n0))

       Ftmin_local(ie)    = MINVAL(elem(ie)%derived%FT(:,:,:))
       Fqmin_local(ie) = MINVAL(elem(ie)%derived%FQ(:,:,:,1))


       psmin_local(ie) = MINVAL(tmp(:,:,ie))
       !======================================================

       usum_local(ie)    = SUM(elem(ie)%state%v(:,:,1,:,n0))
       vsum_local(ie)    = SUM(elem(ie)%state%v(:,:,2,:,n0))
       Wsum_local(ie)    = SUM(elem(ie)%state%w(:,:,:,n0))
       phisum_local(ie)  = SUM(dphi)
       Fusum_local(ie)   = SUM(elem(ie)%derived%FM(:,:,1,:))
       Fvsum_local(ie)   = SUM(elem(ie)%derived%FM(:,:,2,:))

       tsum_local(ie)    = SUM(tdiag)
       if (rsplit>0) then
          dpsum_local(ie)    = SUM(elem(ie)%state%dp3d(:,:,:,n0))
       else
          ! Make sure to initialize this to prevent possible
          ! floating point exceptions.
          dpsum_local(ie)    = 0.0D0
       end if

       Ftsum_local(ie)    = SUM(elem(ie)%derived%FT(:,:,:))
       FQsum_local(ie) = SUM(elem(ie)%derived%FQ(:,:,:,1))

       pssum_local(ie) = SUM(tmp(:,:,ie))
       !======================================================

       global_shared_buf(ie,1) = usum_local(ie)
       global_shared_buf(ie,2) = vsum_local(ie)
       global_shared_buf(ie,3) = tsum_local(ie)
       global_shared_buf(ie,4) = pssum_local(ie)
       global_shared_buf(ie,5) = FUsum_local(ie)
       global_shared_buf(ie,6) = FVsum_local(ie)
       global_shared_buf(ie,7) = FTsum_local(ie)
       global_shared_buf(ie,8) = FQsum_local(ie)
       global_shared_buf(ie,9) = Wsum_local(ie)
       global_shared_buf(ie,10)= dpsum_local(ie)
       global_shared_buf(ie,11)= phisum_local(ie)
    end do

    !JMD This is a Thread Safe Reduction 
    umin_p = ParallelMin(umin_local,hybrid)
    umax_p = ParallelMax(umax_local,hybrid)

    vmin_p = ParallelMin(vmin_local,hybrid)
    vmax_p = ParallelMax(vmax_local,hybrid)

    tmin_p = ParallelMin(tmin_local,hybrid)
    tmax_p = ParallelMax(tmax_local,hybrid)
    
    if (rsplit>0)then
       dpmin_p = ParallelMin(dpmin_local,hybrid)
       dpmax_p = ParallelMax(dpmax_local,hybrid)
    endif

    psmin_p = ParallelMin(psmin_local,hybrid)
    psmax_p = ParallelMax(psmax_local,hybrid)

    FUmin_p = ParallelMin(FUmin_local,hybrid)
    FUmax_p = ParallelMax(FUmax_local,hybrid)

    FVmin_p = ParallelMin(FVmin_local,hybrid)
    FVmax_p = ParallelMax(FVmax_local,hybrid)

    FTmin_p = ParallelMin(FTmin_local,hybrid)
    FTmax_p = ParallelMax(FTmax_local,hybrid)

    FQmin_p = ParallelMin(FQmin_local,hybrid)
    FQmax_p = ParallelMax(FQmax_local,hybrid)

    Wmin_p = ParallelMin(Wmin_local,hybrid)
    Wmax_p = ParallelMax(Wmax_local,hybrid)

    phimin_p = ParallelMin(phimin_local,hybrid)
    phimax_p = ParallelMax(phimax_local,hybrid)

    call wrap_repro_sum(nvars=11, comm=hybrid%par%comm)
    usum_p = global_shared_sum(1)
    vsum_p = global_shared_sum(2)
    tsum_p = global_shared_sum(3)
    pssum_p = global_shared_sum(4)
    FUsum_p = global_shared_sum(5)
    FVsum_p = global_shared_sum(6)
    FTsum_p = global_shared_sum(7)
    FQsum_p = global_shared_sum(8)
    Wsum_p = global_shared_sum(9)
    dpsum_p = global_shared_sum(10)
    phisum_p = global_shared_sum(11)


    scale=1/g                                  ! assume code is using Pa
    if (hvcoord%ps0 <  2000 ) scale=100*scale  ! code is using mb
    ! after scaling, Energy is in J/m**2,  Mass kg/m**2
    ! so rate of changes are W/m**2


    !   mass = integral( ps-p(top) )
    do ie=nets,nete
       tmp(:,:,ie)=elem(ie)%state%ps_v(:,:,n0) 
    enddo
    Mass2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    do ie=nets,nete
       tmp(:,:,ie)=elem(ie)%state%ps_v(:,:,np1) 
    enddo

    !
    !   ptop =  hvcoord%hyai(1)*hvcoord%ps0)  + hvcoord%hybi(1)*ps(i,j)
    !   but we assume hybi(1) is zero at top of atmosphere (pure pressure coordinates)
    !    Mass = (Mass2-(hvcoord%hyai(1)*hvcoord%ps0) )*scale  ! this correction is a constant,
    !                                                         ! ~20 kg/m^2 (effects 4th digit of Mass)
    !   BUT: CAM EUL defines mass as integral( ps ), so to be consistent, ignore ptop contribution; 
    Mass = Mass2*scale

    if(hybrid%masterthread) then
       write(iulog,100) "u     = ",umin_p,umax_p,usum_p
       write(iulog,100) "v     = ",vmin_p,vmax_p,vsum_p
       write(iulog,100) "w     = ",wmin_p,wmax_p,wsum_p
       write(iulog,100) "tdiag = ",tmin_p,tmax_p,tsum_p
       write(iulog,100) "dz(m) = ",phimin_p/g,phimax_p/g,phisum_p/g
       if (rsplit>0) &
       write(iulog,100) "dp    = ",dpmin_p,dpmax_p,dpsum_p

       do q=1,qsize
          write(iulog,100) "qv= ",qvmin_p(q), qvmax_p(q), qvsum_p(q)
       enddo
       write(iulog,100) "ps= ",psmin_p,psmax_p,pssum_p
       write(iulog,'(a,E23.15,a,E23.15,a)') "      M = ",Mass,' kg/m^2',Mass2,' mb     '

       if(fumin_p.ne.fumax_p) write(iulog,100) "fu = ",fumin_p,fumax_p,fusum_p
       if(fvmin_p.ne.fvmax_p) write(iulog,100) "fv = ",fvmin_p,fvmax_p,fvsum_p
       if(ftmin_p.ne.ftmax_p) write(iulog,100) "ft = ",ftmin_p,ftmax_p,ftsum_p
       if(fqmin_p.ne.fqmax_p) write(iulog,100) "fq = ",fqmin_p, fqmax_p, fqsum_p
    end if
 

    if ( test_case(1:10) == "baroclinic" ) then
       ! zeta does not need to be made continious, but we  
       allocate(tmp3d(np,np,nlev,nets:nete))
       call compute_zeta_C0(tmp3d,elem,hybrid,nets,nete,n0)
       tmp=tmp3d(:,:,nlev,:)**2
       relvort = SQRT(global_integral(elem, tmp,hybrid,npts,nets,nete))
       deallocate(tmp3d)
    endif

    if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
       if ( test_case(1:10) == "baroclinic" ) then
          write(iulog,101) "2-norm relative vorticity = ",relvort
101       format (A30,E24.15)

          open(unit=10,file=massfname,form="formatted",position="append")
          write(10,'(99E24.15)')time/secpday,Mass,time2/secpday,KEner(2),IEner(2),PEner(2)
          close(10)

       else
#ifndef CAM
#ifndef _BGL
          open(unit=10,file=massfname,form="formatted",position="append")
          write(10,*)time/secpday,(Mass2-hvcoord%ps0)/hvcoord%ps0
          close(10)
#endif
#endif
       end if
    endif


    if (tl%nstep < tl%nstep0) then
       call t_stopf('prim_printstate')
       return
    endif

    ! ====================================================================
    !	
    !   Detailed diagnostics.  Computed with quantities computed by 
    !   prim_advance() and prim_advection()
    !
    !   only compute on full leapfrog timesteps (tl%nstep >= tl%nstep2)
    ! ====================================================================

!   Compute Energies at time1 and time2 (half levels between leapfrog steps)
    do n=1,4
       
       do ie=nets,nete
          tmp(:,:,ie)=elem(ie)%accum%IEner(:,:,n)
       enddo
       IEner(n) = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
       IEner(n) = IEner(n)*scale
       
       do ie=nets,nete
          tmp(:,:,ie)=elem(ie)%accum%KEner(:,:,n)
       enddo
       KEner(n) = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
       KEner(n) = KEner(n)*scale
       
       do ie=nets,nete
          tmp(:,:,ie)=elem(ie)%accum%PEner(:,:,n)
       enddo
       PEner(n) = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
       PEner(n) = PEner(n)*scale
       TOTE(n)=IEner(n)+PEner(n)+KEner(n)
       
       
       do q=1,qsize
          do ie=nets,nete
             tmp(:,:,ie)=elem(ie)%accum%Qvar(:,:,q,n)
          enddo
          Qvar(q,n) = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
          Qvar(q,n) = Qvar(q,n)*scale
          
          do ie=nets,nete
             tmp(:,:,ie)=elem(ie)%accum%Qmass(:,:,q,n)
          enddo
          Qmass(q,n) = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
          Qmass(q,n) = Qmass(q,n)*scale
          
          if (n==2) then
             do ie=nets,nete
                tmp(:,:,ie)=elem(ie)%accum%Q1mass(:,:,q)
             enddo
             Q1mass(q) = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
             Q1mass(q) = Q1mass(q)*scale
          endif
       enddo
       
    enddo
    


    !
    !   All of these transport terms are at time-tstep = (time1+time2)/2
    !   Vertical transport terms
#ifdef ENERGY_DIAGNOSTICS
    do ie=nets,nete
      tmp(:,:,ie) = elem(ie)%accum%KEhoriz1
    enddo
    KEhorz = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEhorz = KEhorz*scale

    do ie=nets,nete
      tmp(:,:,ie) = elem(ie)%accum%KE1
    enddo
    KEH1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEH1 = KEH1*scale

    do ie=nets,nete
      tmp(:,:,ie) = elem(ie)%accum%KE2
    enddo
    KEH2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEH2 = KEH2*scale



    do ie=nets,nete
      tmp(:,:,ie) = elem(ie)%accum%KEhoriz2
    enddo
    KEhorz2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEhorz2 = KEhorz2*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%KEvert1 + elem(ie)%accum%KEvert2
    enddo
    KEvert = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEvert = KEvert*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%IEvert1 
    enddo
    IEvert = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    IEvert = IEvert*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%PEhoriz1 
    enddo
    PEhorz = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    PEhorz = PEhorz*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%PEvert1
    enddo
    PEvert = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    PEvert = PEvert*scale

    !   KE->IE
    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%T1
    enddo
    T1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    T1 = T1*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%T2
    enddo
    T2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    T2 = T2*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%S1
    enddo
    S1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    S1 = S1*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%S2
    enddo
    S2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    S2 = S2*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%P1
    enddo
    P1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    P1 = P1*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%P2
    enddo
    P2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    P2 = P2*scale



#else
    T1=0; T2=0; S1=0; S2=0; P1=0; P2=0; KEvert=0; IEvert=0; KEhorz=0; IEhorz=0
    KEhorz2=0;
#endif


    if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
       if(moisture /= "dry")then
          if (qsize>=1) then
             write(iulog,'(a,E23.15,a,E23.15,a)') "    dry M = ",Mass-Q1mass(1),' kg/m^2'
          endif
       endif
       
       ! note: diagnostics not yet coded for semi-implicit 
       if (integration == "explicit") then
          
          write(iulog,'(3a25)') "**DYNAMICS**        J/m^2","   W/m^2","W/m^2    "
#ifdef ENERGY_DIAGNOSTICS
          ! terms computed during prim_advance, if ENERGY_DIAGNOSTICS is enabled
          write(iulog,'(a,2e22.14)')'horiz adv KE-u terms abs/rel, should = 0 :',KEhorz, abs(KEH1+KEH2)/dsqrt(KEH1**2+KEH2**2)
          write(iulog,'(a,2e22.14)')'horiz adv KE-u terms, should+to 0        :',KEH1,KEH2
          write(iulog,'(a,2e22.14)')'vert adv etadot KE-uw terms = 0          :',KEvert
          write(iulog,'(a,2e22.14)')'horiz adv KE-w terms, possibly nonzero   :',KEhorz2
          write(iulog,'(a,2e22.14)')'Tot IE advection vert =0                 :',IEvert
          write(iulog,'(a,2e22.14)')'Tot PE advection horiz, vert = 0         :',PEhorz,PEvert       
          write(iulog,'(a,2e22.14)')'(T1+S1 abs/rel = 0)                      :',S1+T1, abs(S1+T1)/dsqrt(S1**2+T1**2)
          write(iulog,'(a,2e22.14)')'(S1,T1)                                  :',S1,T1
          write(iulog,'(a,2e22.14)')'(T2+S2 = 0)                              :',T2+S2
          
          ddt_tot =  (KEner(2)-KEner(1))/(dt)
          ddt_diss = ddt_tot -(KEhorz+KEhorz2+KEvert+T1+T2+P1) 
          write(iulog,'(a,3E22.14)') "KE,d/dt,diss:",KEner(2),ddt_tot,ddt_diss
          
          ddt_tot =  (IEner(2)-IEner(1))/(dt)
          ddt_diss = ddt_tot - (S1+S2+IEVert)
          write(iulog,'(a,3E22.14)') "IE,d/dt,diss:",IEner(2),ddt_tot,ddt_diss
          
          ddt_tot = (PEner(2)-PEner(1))/(dt)
          ddt_diss = ddt_tot - (PEhorz+PEvert+P2)
          write(iulog,'(a,3E22.14)') "PE,d/dt,diss:",PEner(2),ddt_tot,ddt_diss
#else
          write(iulog,'(a,3E22.14)') "KE,d/dt      ",KEner(2),(KEner(2)-KEner(1))/(dt)
          write(iulog,'(a,3E22.14)') "IE,d/dt      ",IEner(2),(IEner(2)-IEner(1))/(dt)
          write(iulog,'(a,3E22.14)') "PE,d/dt      ",PEner(2),(PEner(2)-PEner(1))/(dt)
#endif
          ddt_tot = (TOTE(2)-TOTE(1))/(dt)
          write(iulog,'(a,3E22.14)') " E,dE/dt     ",TOTE(2),ddt_tot
          
          do q=1,qsize
             write(iulog,'(a,i3,a,E22.14,a,2E15.7)') "Q",q,",Q diss, dQ^2/dt:",Qmass(q,2)," kg/m^2",&
                  (Qmass(q,2)-Qmass(q,1))/dt,(Qvar(q,2)-Qvar(q,1))/dt
          enddo

          write(iulog,'(a)') 'Physics tendencies applied by dycore:'
          write(iulog,'(a,2e15.7)') 'dKE/dt(W/m^2): ',(KEner(1)-KEner(3))/dt
          write(iulog,'(a,2e15.7)') 'dIE/dt(W/m^2): ',(IEner(1)-IEner(3))/dt
          write(iulog,'(a,2e15.7)') 'dPE/dt(W/m^2): ',(PEner(1)-PEner(3))/dt
          q=1
          if (qsize>0) write(iulog,'(a,2e15.7)') 'dQ1/dt(kg/sm^2)',(Qmass(q,1)-Qmass(q,3))/dt
       endif

       ! Print change in energy and tracer mass from the start of the simulation
#ifdef CAM
       ! for cam, this information only makes sense for short test runs with
       ! physics disabled (ftype<0).  So disable this diagnostic if ftype>=0:
       if (ftype>=0) TOTE0=-1  
#endif       
       if (TOTE0>0) then
          write(iulog,100) "(E-E0)/E0    ",(TOTE(2)-TOTE0)/TOTE0
          do q=1,qsize
             if(Qmass0(q)>0.0D0) then
                write(iulog,'(a,E23.15,a,i1)') "(Q-Q0)/Q0 ",(Qmass(q,2)-Qmass0(q))/Qmass0(q),"   Q",q
             end if
          enddo
       endif
    endif
    
100 format (A10,3(E23.15))
    
    
    ! initialize "E0" for printout of E-E0/E0
    ! energy is computed in timestep loop.  Save the value 
    ! after the first real timestep, not the half-timesteps used
    ! to bootstrap leapfrog: 
    if (tl%nstep >= tl%nstep0 .and. TOTE0==0 ) then  
       TOTE0=TOTE(2)
       do q=1,qsize
          Qmass0(q)=Qmass(q,1)
       enddo
       time0=time1
    endif

  call t_stopf('prim_printstate')
  end subroutine prim_printstate
   
   

subroutine prim_energy_halftimes(elem,hvcoord,tl,n,t_before_advance,nets,nete)
! 
!  called at the end of a timestep, before timelevel update.  
!  Solution is known at two timelevels.  We originally tried to
!  compute the energy at the midpoint between these levels, but now
!  we compute the energy at:
!
!  t_before_advance==.true.     tl%n0    begining of timestep
!  t_before_advance==.false.    tl%np1   completed timestep, before tl update  
!
!  This routine is called 4 times:  for historical reasons they are out
!  of sequence.  
!
!    n=1:    after CAM forcing, before timestep
!    n=2:    after timestep, including remap, before tl update
!    n=3:    after CAM forcing, before timestep
!    n=4:    (not used)
!
! LF case we use staggered formulas:
!  compute the energies at time half way between timelevels t1 and t2,
!  using timelevels t1,t2.  store in location n of our elem()%accum%* variables
!    
!  in the Dry case, PE (and similar for IE) will be conserved exactly by Leapfrog, 
!  if it is defined as:
!
!   PE(n+.5)  = .5*(  T(n+1) dp(n) + T(n) dp(n+1) )
!
!  KE cannot be exactly conserved with Leapfrog.  but we still use:
!  
!   KE(n+.5) = .5*(  .5 u(n+1)^2 dp(n) +  .5 u(n)^2 dp(n+1) )
!  

    use kinds, only : real_kind
    use dimensions_mod, only : np, np, nlev
    use hybvcoord_mod, only : hvcoord_t
    use element_mod, only : element_t
    use physical_constants, only : Cp, cpwater_vapor
    use physics_mod, only : Virtual_Specific_Heat
    use prim_si_mod, only : preq_hydrostatic

    integer :: t1,n,nets,nete
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (TimeLevel_t), intent(in)       :: tl
    logical :: t_before_advance

    integer :: ie,k,i,j
    real (kind=real_kind), dimension(np,np,nlev)  :: dpt1  ! delta pressure
    real (kind=real_kind), dimension(np,np)  :: E
    real (kind=real_kind), dimension(np,np)  :: suml,suml2,v1,v2
    real (kind=real_kind), dimension(np,np,nlev)  :: sumlk, suml2k
    real (kind=real_kind) :: cp_star1,qval_t1,qval_t2
    real (kind=real_kind) :: pnh(np,np,nlev)   ! nh nonhyrdo pressure
    real (kind=real_kind) :: dpnh(np,np,nlev) ! nh nonhyrdo pressure interfaces
    real (kind=real_kind) :: exner(np,np,nlev)  ! exner nh pressure


    integer:: tmp, t1_qdp   ! the time pointers for Qdp are not the same

    if (t_before_advance) then
       t1=tl%n0
       call TimeLevel_Qdp( tl, qsplit, t1_qdp) !get n0 level into t2_qdp 
    else
       t1=tl%np1
       call TimeLevel_Qdp(tl, qsplit, tmp, t1_qdp) !get np1 into t2_qdp
    endif
    do ie=nets,nete
       do k=1,nlev
          dpt1(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,t1)
       enddo

       call get_pnh_and_exner(hvcoord,elem(ie)%state%theta(:,:,:,t1),dpt1,&
            elem(ie)%state%phi(:,:,:,t1), &
            elem(ie)%state%phis(:,:),elem(ie)%state%Qdp(:,:,:,1,t1_qdp),pnh,dpnh,exner)
   
 !   KE   .5 dp/dn U^2
       do k=1,nlev
          E = ( elem(ie)%state%v(:,:,1,k,t1)**2 +  &
                elem(ie)%state%v(:,:,2,k,t1)**2 + &
                elem(ie)%state%w(:,:,k,t1)**2 ) *0.5 
          sumlk(:,:,k) = E*dpt1(:,:,k)
       enddo
       suml=0
       do k=1,nlev
          suml(:,:) = suml(:,:) + sumlk(:,:,k)
       enddo
       elem(ie)%accum%KEner(:,:,n)=suml(:,:)


    
    !   PE   dp/dn PHIs
       suml=0
       do k=1,nlev
          suml = suml + elem(ie)%state%phi(:,:,k,t1)*dpt1(:,:,k)
       enddo
       elem(ie)%accum%PEner(:,:,n)=suml(:,:)
       

    !  IE = c_p^* dp/deta T - dp/ds phi  + psurf phisurf 
       suml=0
       suml2=0
       do k=1,nlev
          suml(:,:)=suml(:,:)+&
               Cp * dpt1(:,:,k) * elem(ie)%state%theta(:,:,k,t1)*exner(:,:,k) 
          suml2(:,:) = suml2(:,:)  -dpnh(:,:,k)*elem(ie)%state%phi(:,:,k,t1)
       enddo
       elem(ie)%accum%IEner(:,:,n)=suml(:,:) + suml2(:,:) +&
            elem(ie)%state%ps_v(:,:,t1) * elem(ie)%state%phis(:,:)

       if (theta_hydrostatic_mode) then
          ! hydrostatic case: combine IE and PE
          elem(ie)%accum%IEner(:,:,n)=elem(ie)%accum%IEner(:,:,n) + elem(ie)%accum%PEner(:,:,n)
          elem(ie)%accum%PEner(:,:,n) = 0
       endif
       enddo
    
end subroutine prim_energy_halftimes
    
!=======================================================================================================! 
  

subroutine prim_diag_scalars(elem,hvcoord,tl,n,t_before_advance,nets,nete)
! 
!  called at the end of a timestep, before timelevel update.  Solution known at
!  timelevel nm1,n0,np1.  
!
!  This routine is called twice:  n=1:    t1=nm1, t2=n0
!                                 n=2:    t1=n0, t2=np1
! 
!  in all cases:
!     Q1mass(n) = Qdp(t2)    (or:  Q(t2)*dp(t2) )
! 
!  We also compute 
!     Qmass(n)  = .5*(  Q(t2) dp(t1)  Q(t1) dp(t2) )       Leapfrog, concentration formulation
!     Qmass(n)  = .5*(  Qdp(t1) + Qdp(t2) )                Leapfrog, conservation formulation
!     Qmass(n)  =  Qdp(t2)                                 RK2, conservation formulation
!
!
    use kinds, only : real_kind
    use dimensions_mod, only : np, np, nlev
    use hybvcoord_mod, only : hvcoord_t
    use element_mod, only : element_t

    integer :: t1,n,nets,nete
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    
    integer :: ie,k,q
    real (kind=real_kind), dimension(np,np)  :: ps         ! pressure
    real (kind=real_kind), dimension(np,np)  :: dp         ! delta pressure
    real (kind=real_kind), dimension(np,np)  :: E
    real (kind=real_kind), dimension(np,np)  :: suml
    type (TimeLevel_t), intent(in)       :: tl
    logical :: t_before_advance
    integer:: t1_qdp, tmp   ! the time pointer for Qdp are not the same

    if (t_before_advance) then
       t1=tl%n0
       call TimeLevel_Qdp( tl, qsplit, t1_qdp) !get n0 level into t2_qdp 
    else
       t1=tl%np1
       call TimeLevel_Qdp(tl, qsplit, tmp, t1_qdp) !get np1 into t2_qdp (don't need tmp)
    endif


    do ie=nets,nete
       do q=1,qsize
          suml=0
          do k=1,nlev
             suml = suml + elem(ie)%state%Qdp(:,:,k,q,t1_qdp)*elem(ie)%state%Q(:,:,k,q)
          enddo
          elem(ie)%accum%Qvar(:,:,q,n)=suml(:,:)
       enddo
    enddo
    
    do ie=nets,nete
       do q=1,qsize
          suml=0
          do k=1,nlev
             suml = suml + elem(ie)%state%Qdp(:,:,k,q,t1_qdp)
          enddo
          elem(ie)%accum%Q1mass(:,:,q)=suml(:,:)
          elem(ie)%accum%Qmass(:,:,q,n)=suml(:,:)
       enddo
    enddo
 

end subroutine prim_diag_scalars
end module prim_state_mod
