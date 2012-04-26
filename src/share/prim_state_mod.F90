#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_state_mod
  ! ------------------------------
  use kinds, only : real_kind, iulog
  ! ------------------------------
  use dimensions_mod, only : nlev, np, qsize_d, qsize, nelemd, ntrac, ntrac_d
  ! ------------------------------
  use parallel_mod, only :  ordered, parallel_t, syncmp
  use parallel_mod, only: global_shared_buf, global_shared_sum
  ! ------------------------------
  use global_norms_mod, only: wrap_repro_sum
  ! ------------------------------
  use hybrid_mod, only : hybrid_t
  ! ------------------------------
  use physical_constants, only : p0,Cp,g
  ! ------------------------------
  use time_mod, only : tstep, secpday, timelevel_t 
  ! ------------------------------
  use control_mod, only : integration, test_case, runtype, moisture, tracer_advection_formulation,&
       TRACERADV_TOTAL_DIVERGENCE,TRACERADV_UGRADQ,tstep_type,energy_fixer, qsplit, ftype, use_cpstar
  ! ------------------------------
  use hybvcoord_mod, only : hvcoord_t 
  ! ------------------------------
  use global_norms_mod, only : global_integral 
  ! ------------------------------
  use element_mod, only : element_t
  use cslam_control_volume_mod, only : cslam_struct
  ! ------------------------------
  use viscosity_mod, only : compute_zeta_C0_2d
  ! ------------------------------
  use reduction_mod, only : parallelmax,parallelmin
  ! ------------------------------
implicit none
private
  character(len=*), private, parameter :: massfname = "mass.out"

  integer, public :: naccum


  public :: prim_printstate
  public :: prim_printstate_init
  public :: prim_energy_halftimes
  public :: prim_diag_scalars

contains


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

  subroutine prim_printstate(elem, tl,hybrid,hvcoord,nets,nete, cslam)
    type (element_t), intent(in) :: elem(:)
    type (cslam_struct), intent(in), optional :: cslam(:)
    type (TimeLevel_t), target, intent(in) :: tl
    type (hybrid_t),intent(in)     :: hybrid
    type (hvcoord_t), intent(in)   :: hvcoord
    integer,intent(in)             :: nets,nete



    ! Local variables...
    integer :: i,j,k,ie
    integer,parameter  :: type=ORDERED

    real (kind=real_kind)  :: Mass2,Mass,Mass_np1
    real (kind=real_kind)  :: TOTE(4),KEner(4),PEner(4),IEner(4),IEner_wet(4)
    real (kind=real_kind)  :: Qvar(qsize_d,4),Qmass(qsize_d,4),Q1mass(qsize_d)
    real (kind=real_kind)  :: Qmass_added(qsize_d)
    real (kind=real_kind),save  :: time0
    real (kind=real_kind),save  :: TOTE0=0,Qmass0(qsize_d)=0
    real (kind=real_kind)  :: I_div(nlev)
    real (kind=real_kind)  :: I_vort(nlev)

    real (kind=real_kind)  :: tmp(np,np,nets:nete)
    real (kind=real_kind)  :: tmp1(nets:nete)
    real (kind=real_kind)  :: ps(np,np)
    real (kind=real_kind)  :: dp(np,np)
    !    real (kind=real_kind)  :: E(np,np)


    real (kind=real_kind) :: umin_local(nets:nete), umax_local(nets:nete), usum_local(nets:nete), &
         vmin_local(nets:nete), vmax_local(nets:nete), vsum_local(nets:nete), &
         tmin_local(nets:nete), tmax_local(nets:nete), tsum_local(nets:nete), &
         psmin_local(nets:nete),psmax_local(nets:nete),pssum_local(nets:nete), &
         fumin_local(nets:nete),fumax_local(nets:nete),fusum_local(nets:nete), &
         fvmin_local(nets:nete),fvmax_local(nets:nete),fvsum_local(nets:nete), &
         ftmin_local(nets:nete),ftmax_local(nets:nete),ftsum_local(nets:nete), &
         fqmin_local(nets:nete),fqmax_local(nets:nete),fqsum_local(nets:nete), &
         omegamin_local(nets:nete),omegamax_local(nets:nete),omegasum_local(nets:nete)


    real (kind=real_kind) :: umin_p, &
         vmin_p, &
         tmin_p, &
         qvmin_p(qsize_d), cmin(ntrac_d),&
         psmin_p


    real (kind=real_kind) :: umax_p, &
         vmax_p, &
         tmax_p, &
         qvmax_p(qsize_d), cmax(ntrac_d),&
         psmax_p

    real (kind=real_kind) :: usum_p, &
         vsum_p, &
         tsum_p, &
         qvsum_p(qsize_d), csum(ntrac_d),&
         pssum_p

    real (kind=real_kind) :: fusum_p, fvsum_p, ftsum_p, fqsum_p
    real (kind=real_kind) :: fumin_p, fvmin_p, ftmin_p, fqmin_p
    real (kind=real_kind) :: fumax_p, fvmax_p, ftmax_p, fqmax_p
    real (kind=real_kind) :: omegamax_p, omegamin_p, omegasum_p


    real(kind=real_kind) :: vsum_t(1), relvort
    real(kind=real_kind) :: v1, v2, vco(np,np,2,nlev)

    real (kind=real_kind) :: time, time2,time1, scale, dt, dt_split
    real (kind=real_kind) :: KEvert,IEvert,T1,T2,T2_s,T2_m,S1,S2,S1_wet
    real (kind=real_kind) :: KEhorz,IEhorz,IEhorz_wet,IEvert_wet
    real (kind=real_kind) :: ddt_tot,ddt_diss
    integer               :: n0, nm1, pnm1, np1
    integer               :: npts,n,q

    if (.not. present(cslam) .and. ntrac>0) then
       print *,'ERROR: prim_state_mod.F90: optional cslam argument requred if ntrac>0'
    endif

    TOTE     = 0
    KEner    = 0
    PEner    = 0
    IEner    = 0
    IEner_wet = 0
    ! dynamics timelevels
    n0=tl%n0
    nm1=tl%nm1
    np1=tl%np1

    ! forcing timelevel
#ifdef CAM
    pnm1 = 1
#else
    pnm1=tl%nm1  
#endif

    dt = tstep*qsplit
    !
    !   state variables in n0 are at time =  'time' 
    !   state variables in nm1 are at time =  'time-dt' 
    !
    !   Other diagnostics were computed during the timestep, and thus:
    !   Energy variables *(:,:,2)       time-0.5*dt    "time2"
    !   Energy variables *(:,:,1)       time-1.5*dt    "time1"
    !   dE/dt variables                 time-1.0*dt 
    ! advecting Q
    !   Mass                            time (conserved w/o forcing, so never changes)
    !   conserved mass is at half levels
    ! 
    ! advecting Qdp with leapfrog (tstep_type=0)
    !   Note: Qmass conserved at whole timelevels, but we compute at half-time levels 
    !         so change in Qmass will exactly balance forcing term.  
    !         could also use time, time-2*dt, but then Qmass uses a different dt from all other terms
    !
    !   Mass                            time ( Mass-Q1mass(:,:) conserved )
    !   Q1Mass                          time  
    !   Qmass(:,:,4)                    time-0.5*dt  after robert
    !   Qmass(:,:,3)                    time-0.5*dt  after forcing
    !   Qmass(:,:,2)                    time-0.5*dt  after timestep to time
    !   Qmass(:,:,1)                    time-1.5*dt  
    !
    !   Qmass is computed at t+.5, t-.5 so that dQ/dt exactly balances FQ
    !   To check conservation of dry mass, use Qmass at time t ("Q1mass")
    !
    ! advecting Qdp with RK-SSP (tstep_type=1)
    !   Mass                            time ( Mass-Q1mass(:,:) conserved )
    !   Q1Mass                          time  
    !   Qmass(:,:,2)                    time after timestep to time
    !   Qmass(:,:,1)                    time-dt
    !   Qmass is computed at t+1, t so that dQ/dt exactly balances FQ
    !   Qmass(:,:,3) , Qmass(:,:,4) and Q1Mass are all the same since there is
    !   no robert filter
    !
    time=tl%nstep* dt
    time2 = time - dt/2
    time1 = time - 3*dt/2



    vsum_t(1) = 0.0D0

    ! npts = np
    npts=SIZE(elem(1)%state%lnps(:,:,n0),1)

    do q=1,qsize
       do ie=nets,nete
          tmp1(ie) = MINVAL(elem(ie)%state%Q(:,:,:,q,n0))
       enddo
       qvmin_p(q) = ParallelMin(tmp1,hybrid)
       do ie=nets,nete
          tmp1(ie) = MAXVAL(elem(ie)%state%Q(:,:,:,q,n0))
       enddo
       qvmax_p(q) = ParallelMax(tmp1,hybrid)
       do ie=nets,nete
          global_shared_buf(ie,1) = SUM(elem(ie)%state%Q(:,:,:,q,n0))
       enddo
       call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
       qvsum_p(q) = global_shared_sum(1)
    enddo
    do q=1,ntrac
       do ie=nets,nete
          tmp1(ie) = MINVAL(cslam(ie)%c(:,:,:,q,n0))
       enddo
       cmin(q) = ParallelMin(tmp1,hybrid)
       do ie=nets,nete
          tmp1(ie) = MAXVAL(cslam(ie)%c(:,:,:,q,n0))
       enddo
       cmax(q) = ParallelMax(tmp1,hybrid)
       do ie=nets,nete
          global_shared_buf(ie,1) = SUM(cslam(ie)%c(:,:,:,q,n0))
       enddo
       call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
       csum(q) = global_shared_sum(1)
    enddo

    do ie=nets,nete

       tmp(:,:,ie)=elem(ie)%state%ps_v(:,:,n0)
       !       tmp(:,:,ie)=EXP(elem(ie)%state%lnps(:,:,n0))


       !======================================================  
       umax_local(ie)    = MAXVAL(elem(ie)%state%v(:,:,1,:,n0))
       vmax_local(ie)    = MAXVAL(elem(ie)%state%v(:,:,2,:,n0))

       fumax_local(ie)    = MAXVAL(elem(ie)%derived%FM(:,:,1,:,pnm1))
       fvmax_local(ie)    = MAXVAL(elem(ie)%derived%FM(:,:,2,:,pnm1))

       tmax_local(ie)    = MAXVAL(elem(ie)%state%T(:,:,:,n0))

       psmax_local(ie) = MAXVAL(tmp(:,:,ie))
       ftmax_local(ie)    = MAXVAL(elem(ie)%derived%FT(:,:,:,pnm1))
       fqmax_local(ie)    = MAXVAL(elem(ie)%derived%FQ(:,:,:,1,pnm1))
       omegamax_local(ie)    = MAXVAL(elem(ie)%derived%Omega_p(:,:,:))
       !======================================================

       umin_local(ie)    = MINVAL(elem(ie)%state%v(:,:,1,:,n0))
       vmin_local(ie)    = MINVAL(elem(ie)%state%v(:,:,2,:,n0))

       Fumin_local(ie)    = MINVAL(elem(ie)%derived%FM(:,:,1,:,pnm1))
       Fvmin_local(ie)    = MINVAL(elem(ie)%derived%FM(:,:,2,:,pnm1))

       tmin_local(ie)    = MINVAL(elem(ie)%state%T(:,:,:,n0))


       Ftmin_local(ie)    = MINVAL(elem(ie)%derived%FT(:,:,:,pnm1))
       Fqmin_local(ie) = MINVAL(elem(ie)%derived%FQ(:,:,:,1,pnm1))
       Omegamin_local(ie) = MINVAL(elem(ie)%derived%Omega_p(:,:,:))


       psmin_local(ie) = MINVAL(tmp(:,:,ie))
       !======================================================

       usum_local(ie)    = SUM(elem(ie)%state%v(:,:,1,:,n0))
       vsum_local(ie)    = SUM(elem(ie)%state%v(:,:,2,:,n0))
       Fusum_local(ie)    = SUM(elem(ie)%derived%FM(:,:,1,:,pnm1))
       Fvsum_local(ie)    = SUM(elem(ie)%derived%FM(:,:,2,:,pnm1))

       tsum_local(ie)    = SUM(elem(ie)%state%t(:,:,:,n0))

       Ftsum_local(ie)    = SUM(elem(ie)%derived%FT(:,:,:,pnm1))
       FQsum_local(ie) = SUM(elem(ie)%derived%FQ(:,:,:,1,pnm1))
       Omegasum_local(ie) = SUM(elem(ie)%derived%Omega_p(:,:,:))

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
       global_shared_buf(ie,9) = Omegasum_local(ie)
    end do

    !JMD This is a Thread Safe Reduction 
    umin_p = ParallelMin(umin_local,hybrid)
    umax_p = ParallelMax(umax_local,hybrid)

    vmin_p = ParallelMin(vmin_local,hybrid)
    vmax_p = ParallelMax(vmax_local,hybrid)

    tmin_p = ParallelMin(tmin_local,hybrid)
    tmax_p = ParallelMax(tmax_local,hybrid)

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

    Omegamin_p = ParallelMin(Omegamin_local,hybrid)
    Omegamax_p = ParallelMax(Omegamax_local,hybrid)

    call wrap_repro_sum(nvars=9, comm=hybrid%par%comm)
    usum_p = global_shared_sum(1)
    vsum_p = global_shared_sum(2)
    tsum_p = global_shared_sum(3)
    pssum_p = global_shared_sum(4)
    FUsum_p = global_shared_sum(5)
    FVsum_p = global_shared_sum(6)
    FTsum_p = global_shared_sum(7)
    FQsum_p = global_shared_sum(8)
    Omegasum_p = global_shared_sum(9)

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
    Mass_np1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)

    !
    !   ptop =  hvcoord%hyai(1)*hvcoord%ps0)  + hvcoord%hybi(1)*ps(i,j)
    !   but we assume hybi(1) is zero at top of atmosphere (pure pressure coordinates)
    !    Mass = (Mass2-(hvcoord%hyai(1)*hvcoord%ps0) )*scale  ! this correction is a constant, ~20 kg/m^2 (effects 4th digit of Mass)
    !   BUT: CAM EUL defines mass as integral( ps ), so to be consistent, ignore ptop contribution; 
    Mass = Mass2*scale
    Mass_np1 = Mass_np1*scale


    !   sum the mass added by the mass fixer
    !   mass fixer computed mass within each element.  now just global sum:
    !   result in kg/m^2
    do q=1,qsize
       do ie=nets,nete
          tmp1(ie) = elem(ie)%accum%mass_added(q)
          global_shared_buf(ie,1) = tmp1(ie)
       enddo
       call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
       Qmass_added(q) = global_shared_sum(1)*scale
    enddo


    if(hybrid%masterthread) then
       write(iulog,100) "u     = ",umin_p,umax_p,usum_p
       write(iulog,100) "v     = ",vmin_p,vmax_p,vsum_p
       write(iulog,100) "omega = ",omegamin_p,omegamax_p,omegasum_p
       write(iulog,100) "t     = ",tmin_p,tmax_p,tsum_p
       do q=1,qsize
          write(iulog,100) "qv= ",qvmin_p(q), qvmax_p(q), qvsum_p(q)
       enddo
       do q=1,ntrac
          write(iulog,100) " c= ",cmin(q), cmax(q), csum(q)
       enddo
       write(iulog,100) "ps= ",psmin_p,psmax_p,pssum_p
       write(iulog,'(a,E23.15,a,E23.15,a)') "      M = ",Mass,' kg/m^2',Mass2,' mb     '

       if(fumin_p.ne.fumax_p) write(iulog,100) "fu = ",fumin_p,fumax_p,fusum_p
       if(fvmin_p.ne.fvmax_p) write(iulog,100) "fv = ",fvmin_p,fvmax_p,fvsum_p
       if(ftmin_p.ne.ftmax_p) write(iulog,100) "ft = ",ftmin_p,ftmax_p,ftsum_p
       if(fqmin_p.ne.fqmax_p) write(iulog,100) "fq = ",fqmin_p, fqmax_p, fqsum_p
       do q=1,qsize
          if (Qmass_added(q) /= 0) then
             write(iulog,'(a,i1,a,E23.15)') "Q",q," qnegfix mass added:: ",Qmass_added(q)
          endif
       enddo
    endif
 

    if ( test_case(1:10) == "baroclinic" ) then
       ! zeta does not need to be made continious, but we  
       call compute_zeta_C0_2d(tmp,elem,hybrid,nets,nete,n0,nlev)
       tmp=tmp**2
       relvort = SQRT(global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete))
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


    if (tl%nstep < tl%nstep0) return

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
          tmp(:,:,ie)=elem(ie)%accum%IEner_wet(:,:,n)
       enddo
       IEner_wet(n) = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
       IEner_wet(n) = IEner_wet(n)*scale
       
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
       tmp(:,:,ie) = elem(ie)%accum%KEvert1 + elem(ie)%accum%KEvert2
    enddo
    KEvert = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEvert = KEvert*scale
    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%IEvert1 + elem(ie)%accum%IEvert2
    enddo
    IEvert = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    IEvert = IEvert*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%KEhorz1 + elem(ie)%accum%KEhorz2
    enddo
    KEhorz = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEhorz = KEhorz*scale
    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%IEhorz1 + elem(ie)%accum%IEhorz2
    enddo
    IEhorz = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    IEhorz = IEhorz*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%IEhorz1_wet + elem(ie)%accum%IEhorz2_wet
    enddo
    IEhorz_wet = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    IEhorz_wet = IEhorz_wet*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%IEvert1_wet + elem(ie)%accum%IEvert2_wet
    enddo
    IEvert_wet = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    IEvert_wet = IEvert_wet*scale

    !   KE->PE
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
       tmp(:,:,ie) = elem(ie)%accum%T2_s
    enddo
    T2_s = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    T2_s = T2_s*scale
    T2_m = T2 - T2_s

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%S1
    enddo
    S1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    S1 = S1*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%S1_wet
    enddo
    S1_wet = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    S1_wet = S1_wet*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%S2
    enddo
    S2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    S2 = S2*scale


#else
    T1=0; T2=0; T2_s=0; T2_m=0; S1_wet = 0; S1=0; S2=0; KEvert=0; IEvert=0; KEhorz=0; IEhorz=0
#endif


    if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
       if(moisture /= "dry")then
          if (qsize>=1) then
             if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
                write(iulog,'(a,E23.15,a,E23.15,a)') "    dry M = ",Mass-Q1mass(1),' kg/m^2'
             else
                write(iulog,'(a,E23.15,a,E23.15,a)') "    dry M = ",.5*(Mass_np1+Mass)-Qmass(1,2),' kg/m^2'
             endif
          endif
       endif
       
       ! note: diagnostics not yet coded for semi-implicit 
       if (integration == "explicit") then
          
          
          write(iulog,'(3a25)') "**DYNAMICS**        J/m^2","   W/m^2","W/m^2    "
#ifdef ENERGY_DIAGNOSTICS
          ! terms computed during prim_advance, if ENERGY_DIAGNOSTICS is enabled
          write(iulog,'(a,2e22.14)') 'Tot KE advection horiz, vert:      ',-KEhorz,-KEvert
          write(iulog,'(a,2e22.14)') 'Tot IE dry advection horiz, vert:  ',-IEhorz,-IEvert
          if (use_cpstar==1) then
             write(iulog,'(a,2e22.14)') 'Tot IE wet advection horiz, vert:  ',-IEhorz_wet,-IEvert_wet
          endif
          
          write(iulog,'(a,2e22.14)') 'Transfer:   KE->IE, KE->PE:        ',-(T1+T2_m), -T2_s
          write(iulog,'(a,2e22.14)') 'Transfer:   IE->KE, PE->KE:        ',-S1,-S2
          if (use_cpstar==1) then
             write(iulog,'(a,2e22.14)') 'Transfer:   IE->KE(wet only):      ',-S1_wet
          endif
          
          ddt_tot =  (KEner(2)-KEner(1))/(dt)
          ddt_diss = ddt_tot -(T1+T2) 
          write(iulog,'(a,3E22.14)') "KE,d/dt,diss:",KEner(2),ddt_tot,ddt_diss
          
          if (use_cpstar==1) then
             ddt_tot =  ( (IEner(2)-IEner_wet(2)) - (IEner(1)-IEner_wet(1))  )/dt
             ddt_diss = ddt_tot - (S1 -S1_wet)
             write(iulog,'(a,3E22.14,a)') "IE,d/dt,diss:",IEner(2)-IEner_wet(2),ddt_tot,ddt_diss,' (dry)'
             ddt_tot =  ( IEner_wet(2) - IEner_wet(1) )/dt 
             ddt_diss =  ddt_tot - S1_wet
             write(iulog,'(a,3E22.14,a)') "IE,d/dt,diss:",IEner_wet(2),ddt_tot,ddt_diss,' (wet)'
          else
             ddt_tot =  (IEner(2)-IEner(1))/(dt)
             ddt_diss = ddt_tot - S1 
             write(iulog,'(a,3E22.14)') "IE,d/dt,diss:",IEner(2),ddt_tot,ddt_diss
          endif
          
          ddt_tot = (PEner(2)-PEner(1))/(dt)
          ddt_diss = ddt_tot - S2
          write(iulog,'(a,3E22.14)') "PE,d/dt,diss:",PEner(2),ddt_tot,ddt_diss
#else
          write(iulog,'(a,3E22.14)') "KE,d/dt      ",KEner(2),(KEner(2)-KEner(1))/(dt)
          write(iulog,'(a,3E22.14)') "IE,d/dt      ",IEner(2),(IEner(2)-IEner(1))/(dt)
          write(iulog,'(a,3E22.14)') "PE,d/dt      ",PEner(2),(PEner(2)-PEner(1))/(dt)
#endif
          ddt_tot = (TOTE(2)-TOTE(1))/(dt)
          write(iulog,'(a,3E22.14)') " E,dE/dt     ",TOTE(2),ddt_tot
          do q=1,qsize
             write(iulog,'(a,i1,a,E22.14,a,2E15.7)') "Q",q,",Q diss, dQ^2/dt:",Qmass(q,2)," kg/m^2",&
                  (Qmass(q,2)-Qmass(q,1))/dt,(Qvar(q,2)-Qvar(q,1))/dt
          enddo
          
          ! LF code diagnostics
          if (tstep_type==0) then  ! leapfrog
             write(iulog,'(a)') 'Robert filter, Physics (except adjustments):'
             write(iulog,'(a,2e15.7)') 'dKE/dt(W/m^2): ',(KEner(4)-KEner(3))/dt,(KEner(3)-KEner(2))/dt
             write(iulog,'(a,2e15.7)') 'dIE/dt(W/m^2): ',(IEner(4)-IEner(3))/dt,(IEner(3)-IEner(2))/dt
             write(iulog,'(a,2e15.7)') 'dPE/dt(W/m^2): ',(PEner(4)-PEner(3))/dt,(PEner(3)-PEner(2))/dt
          endif
          ! RK code diagnostics
          if (tstep_type==1) then  
             write(iulog,'(a)') 'Energy Fixer, Physics (except adjustments):'
             write(iulog,'(a,2e15.7)') 'dKE/dt(W/m^2): ',(KEner(3)-KEner(2))/dt,(KEner(1)-KEner(4))/dt
             write(iulog,'(a,2e15.7)') 'dIE/dt(W/m^2): ',(IEner(3)-IEner(2))/dt,(IEner(1)-IEner(4))/dt
             write(iulog,'(a,2e15.7)') 'dPE/dt(W/m^2): ',(PEner(3)-PEner(2))/dt,(PEner(1)-PEner(4))/dt
          endif
       endif

       ! Print change in energy and tracer mass from the start of the simulation
#ifdef CAM
       ! for cam, this information only makes sense for short test runs with
       ! physics disabled (ftype<0).  So disable this diagnostic if ftype>=0:
       if (ftype>=0) TOTE0=-1  
#endif       
       if (TOTE0>0) then
          write(iulog,100) "(E-E0)/E0    ",(TOTE(4)-TOTE0)/TOTE0
          do q=1,qsize
             if(Qmass0(q)>0.0) then
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
       TOTE0=TOTE(4)
       do q=1,qsize
          Qmass0(q)=Qmass(q,1)
       enddo
       time0=time1
    endif

  end subroutine prim_printstate
   
   




subroutine prim_energy_halftimes(elem,hvcoord,tl,n,t_before_advance,nets,nete,tQ)
! 
!  called at the end of a timestep, before timelevel update.  Solution known at
!  dynamics:     nm1,  n0,  np1.  
!
!
!  This routien is called 4 times:  n=1:    t1=nm1, t2=n0
!                                   n=2:    t1=n0, t2=np1
!                                   n=3:    t1=n0, t2=np1
!                                   n=4:    t1=n0, t2=np1
!
!  Solution is given at times u(t-1),u(t),u(t+1)
!  E(n=1) = energy at time t-.5
!  E(n=2) = energy at time t+.5
!  E(n=3) = energy at time t+.5  after Forcing applied
!  E(n=4) = energy at time t+.5  after Robert filter
!
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
!
    use kinds, only : real_kind
    use dimensions_mod, only : np, np, nlev
    use control_mod, only : use_cpstar
    use hybvcoord_mod, only : hvcoord_t
    use element_mod, only : element_t
    use physical_constants, only : Cp, cpwater_vapor
    use physics_mod, only : Virtual_Specific_Heat

    integer :: t1,t2,n,nets,nete
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (TimeLevel_t), intent(in)       :: tl
    integer, intent(in),optional    :: tQ
    logical :: t_before_advance

    integer :: ie,k,i,j,nm_f
    real (kind=real_kind), dimension(np,np,nlev)  :: dpt1,dpt2   ! delta pressure
    real (kind=real_kind), dimension(np,np)  :: E
    real (kind=real_kind), dimension(np,np)  :: suml,suml2,v1,v2
    real (kind=real_kind), dimension(np,np,nlev)  :: sumlk, suml2k
    real (kind=real_kind) :: cp_star1,cp_star2,qval_t1,qval_t2
    logical tstagger

    nm_f = 1
    if (t_before_advance) then
       t1=tl%nm1
       t2=tl%n0
    else
       t1=tl%n0
       t2=tl%np1
    endif

    tstagger = .true.

! energy_fixer
!     -2         disabled, but compute energy non-staggered in time
!     -1         disabled, but compute energy non-staggered in time
!     0          disabled, but compute energy staggered in time 
!     1          cp_star(t1)*dp(t1)*T(t2)
!     2          cp*dp(t1)*T(t2)
!     3          cp_star(t2)*dp(t2)*T(t2)
!     4          cp*dp(t2)*T(t2)
!
    if (energy_fixer==-2) tstagger = .false.
    if (energy_fixer==-1) tstagger = .false.
    if (energy_fixer==3) tstagger = .false.
    if (energy_fixer==4) tstagger = .false.


    !   IE   Cp*dpdn*T  + (Cpv-Cp) Qdpdn*T
    !        Cp*dpdn(n)*T(n+1) + (Cpv-Cp) Qdpdn(n)*T(n+1)
    !        [Cp + (Cpv-Cp) Q(n)] *dpdn(n)*T(n+1) 
    do ie=nets,nete

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
          dpt1(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,t1)
          dpt2(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,t2)
       enddo

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,i,j,cp_star1,cp_star2,qval_t1,qval_t2)
#endif
       do k=1,nlev
          sumlk(:,:,k)=0
          suml2k(:,:,k)=0
          do i=1,np
          do j=1,np
             if(use_cpstar == 1)  then
                ! Cp_star = cp + (Cpwater_vapor - cp)*qval
                if (present(tQ)) then
                   ! we should interpolate to t +/- dt_dynamics/2, but not worth the trouble
                   cp_star2= Virtual_Specific_Heat(elem(ie)%state%Q(i,j,k,1,tQ))
                   cp_star1= cp_star2  
                else
                   if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
                      qval_t1 = elem(ie)%state%Qdp(i,j,k,1,t1)/dpt1(i,j,k)
                      qval_t2 = elem(ie)%state%Qdp(i,j,k,1,t2)/dpt2(i,j,k)
                   else
                      qval_t1 = elem(ie)%state%Q(i,j,k,1,t1)
                      qval_t2 = elem(ie)%state%Q(i,j,k,1,t2)
                   endif
                   cp_star1= Virtual_Specific_Heat(qval_t1)
                   cp_star2= Virtual_Specific_Heat(qval_t2)
                endif
             else
                cp_star1=cp
                cp_star2=cp
             endif
             if (tstagger) then
                sumlk(i,j,k) = sumlk(i,j,k) + cp_star1*elem(ie)%state%T(i,j,k,t2) *dpt1(i,j,k)/2 
                sumlk(i,j,k) = sumlk(i,j,k) + Cp_star2*elem(ie)%state%T(i,j,k,t1) *dpt2(i,j,k)/2
                suml2k(i,j,k) = suml2k(i,j,k) + (cp_star1-cp)*elem(ie)%state%T(i,j,k,t2) *dpt1(i,j,k)/2 
                suml2k(i,j,k) = suml2k(i,j,k) + (cp_star2-cp)*elem(ie)%state%T(i,j,k,t1) *dpt2(i,j,k)/2
             else
                sumlk(i,j,k) = sumlk(i,j,k) + Cp_star2*elem(ie)%state%T(i,j,k,t2) *dpt2(i,j,k)
                suml2k(i,j,k) = suml2k(i,j,k) + (cp_star2-cp)*elem(ie)%state%T(i,j,k,t2) *dpt2(i,j,k)
             endif
          enddo
          enddo
       enddo
       suml=0
       suml2=0
       do k=1,nlev
          suml(:,:) = suml(:,:) + sumlk(:,:,k)
          suml2(:,:) = suml2(:,:) + suml2k(:,:,k)
       enddo
       elem(ie)%accum%IEner(:,:,n)=suml(:,:)
       elem(ie)%accum%IEner_wet(:,:,n)=suml2(:,:)

    
    !   KE   .5 dp/dn U^2
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,E)
#endif
       do k=1,nlev
          if (tstagger) then
             E = (elem(ie)%state%v(:,:,1,k,t2)**2 +  &
                  elem(ie)%state%v(:,:,2,k,t2)**2 ) / 2 
             sumlk(:,:,k) = E*dpt1(:,:,k)/2
             E = (elem(ie)%state%v(:,:,1,k,t1)**2 +  &
                  elem(ie)%state%v(:,:,2,k,t1)**2 ) / 2 
             sumlk(:,:,k) = sumlk(:,:,k) + E*dpt2(:,:,k)/2
          else
             E = (elem(ie)%state%v(:,:,1,k,t2)**2 +  &
                  elem(ie)%state%v(:,:,2,k,t2)**2 ) / 2 
             sumlk(:,:,k) = E*dpt2(:,:,k)
          endif
       enddo
       suml=0
       do k=1,nlev
          suml(:,:) = suml(:,:) + sumlk(:,:,k)
       enddo
       elem(ie)%accum%KEner(:,:,n)=suml(:,:)

    
    !   PE   dp/dn PHIs
       suml=0
       do k=1,nlev
          if (tstagger) then
             suml = suml + elem(ie)%state%phis(:,:)*dpt1(:,:,k)/2
             suml = suml + elem(ie)%state%phis(:,:)*dpt2(:,:,k)/2
          else
             suml = suml + elem(ie)%state%phis(:,:)*dpt2(:,:,k)
          endif
       enddo
       elem(ie)%accum%PEner(:,:,n)=suml(:,:)
    enddo
    
end subroutine prim_energy_halftimes
    
    

subroutine prim_diag_scalars(elem,hvcoord,tl,n,t_before_advance,nets,nete)
! 
!  called at the end of a timestep, before timelevel update.  Solution known at
!  timelevel nm1,n0,np1.  
!
!  This routien is called twice:  n=1:    t1=nm1, t2=n0
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

    integer :: t1,t2,n,nets,nete
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord

    integer :: ie,k,q,nm_f
    real (kind=real_kind), dimension(np,np)  :: ps         ! pressure
    real (kind=real_kind), dimension(np,np)  :: dp         ! delta pressure
    real (kind=real_kind), dimension(np,np)  :: E
    real (kind=real_kind), dimension(np,np)  :: suml
    type (TimeLevel_t), intent(in)       :: tl
    logical :: t_before_advance

    nm_f = 1
    if (t_before_advance) then
       t1=tl%nm1     
       t2=tl%n0
    else
       t1=tl%n0
       t2=tl%np1
    endif

    !
    !  leapfrog of Qdp.  compute at t1 and t2, take average
    !
    if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE .and. &
         tstep_type==0) then

    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k,suml)
#endif
    do q=1,qsize
       suml=0
       do k=1,nlev
          suml = suml + .5*( elem(ie)%state%Qdp(:,:,k,q,t1)*elem(ie)%state%Q(:,:,k,q,t1) + & 
                             elem(ie)%state%Qdp(:,:,k,q,t2)*elem(ie)%state%Q(:,:,k,q,t2) )
       enddo
       elem(ie)%accum%Qvar(:,:,q,n)=suml(:,:)
    enddo
    enddo

    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k,suml)
#endif
    do q=1,qsize
       suml=0
       do k=1,nlev
          suml = suml + .5*(elem(ie)%state%Qdp(:,:,k,q,t1) + elem(ie)%state%Qdp(:,:,k,q,t2) )
       enddo
       elem(ie)%accum%Qmass(:,:,q,n)=suml(:,:)

       suml=0
       do k=1,nlev
          suml = suml + elem(ie)%state%Qdp(:,:,k,q,t2)
       enddo
       elem(ie)%accum%Q1mass(:,:,q)=suml(:,:)
    enddo
    enddo
    endif


    !
    !  RK2 forward scheme.  compute everything at t2
    !
    if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE .and. &
         tstep_type==1) then

    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k,suml)
#endif
    do q=1,qsize
       suml=0
       do k=1,nlev
          suml = suml + elem(ie)%state%Qdp(:,:,k,q,t2)*elem(ie)%state%Q(:,:,k,q,t2)
       enddo
       elem(ie)%accum%Qvar(:,:,q,n)=suml(:,:)
    enddo
    enddo

    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k,suml)
#endif
    do q=1,qsize
       suml=0
       do k=1,nlev
          suml = suml + elem(ie)%state%Qdp(:,:,k,q,t2)
       enddo
       elem(ie)%accum%Q1mass(:,:,q)=suml(:,:)
       elem(ie)%accum%Qmass(:,:,q,n)=suml(:,:)
    enddo
    enddo
    endif


    !
    !  leapfrog of Q - use staggered in time formula
    !
    if (tracer_advection_formulation==TRACERADV_UGRADQ) then
    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k,suml,ps,dp)
#endif
    do q=1,qsize
       suml=0
       ps(:,:)=elem(ie)%state%ps_v(:,:,t1)
       do k=1,nlev
          dp(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*ps
          suml = suml + elem(ie)%state%Q(:,:,k,q,t2)**2 *dp(:,:) 
       enddo
       ps(:,:)=elem(ie)%state%ps_v(:,:,t2)
       do k=1,nlev
          dp(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*ps
          suml = suml + elem(ie)%state%Q(:,:,k,q,t1)**2 *dp(:,:)  
       enddo
       elem(ie)%accum%Qvar(:,:,q,n)=suml(:,:)/2
    enddo
    enddo
    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(q,k,suml,ps,dp)
#endif
    do q=1,qsize
       suml=0
       ps(:,:)=elem(ie)%state%ps_v(:,:,t1)
       do k=1,nlev
          dp(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*ps
          suml = suml + elem(ie)%state%Q(:,:,k,q,t2) *dp(:,:) 
       enddo
       ps(:,:)=elem(ie)%state%ps_v(:,:,t2)
       do k=1,nlev
          dp(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*ps
          suml = suml + elem(ie)%state%Q(:,:,k,q,t1) *dp(:,:)  
       enddo
       elem(ie)%accum%Qmass(:,:,q,n)=suml(:,:)/2
    enddo
    enddo



    endif


end subroutine prim_diag_scalars
end module prim_state_mod
