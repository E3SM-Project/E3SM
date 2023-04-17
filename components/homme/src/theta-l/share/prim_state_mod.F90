! July 2019 O. Guba Added extrema locations in output

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_state_mod

  use kinds,            only: real_kind, iulog
  use dimensions_mod,   only: nlev, np, qsize_d, qsize, nlevp
  use parallel_mod,     only:  iam, ordered, parallel_t, syncmp, abortmp
  use parallel_mod,     only: global_shared_buf, global_shared_sum
  use global_norms_mod, only: wrap_repro_sum
  use hybrid_mod,       only: hybrid_t
  use time_mod,         only: tstep, secpday, timelevel_t, TimeLevel_Qdp, time_at
  use control_mod,      only: integration, test_case,  use_moisture, &
                              qsplit, ftype, rsplit,&
                              theta_hydrostatic_mode
  use hybvcoord_mod,    only: hvcoord_t
  use global_norms_mod, only: global_integral, linf_snorm, l1_snorm, l2_snorm
  use element_mod,      only: element_t
  use element_ops,      only: get_field, get_phi, get_field_i
  use element_state,    only: max_itercnt,max_deltaerr,max_reserr,diagtimes
  use eos,              only: pnh_and_exner_from_eos
  use viscosity_mod,    only: compute_zeta_C0
  use reduction_mod,    only: parallelmax,parallelmin,parallelmaxwithindex,parallelminwithindex
  use perf_mod,         only: t_startf, t_stopf
  use physical_constants, only : p0,Cp,g,Rgas

implicit none
private
  character(len=*), private, parameter :: massfname = "mass.out"
  real (kind=real_kind)                :: BIGVAL=1e15

  integer, public :: naccum


  public :: prim_printstate
  public :: prim_printstate_init
  public :: prim_energy_halftimes
  public :: prim_diag_scalars

contains
!=======================================================================================================! 


  subroutine prim_printstate_init(par, elem)
    use dimensions_mod, only: nelemd

    type (parallel_t) :: par
    type (element_t), pointer :: elem(:)

    real (kind=real_kind) :: time
    integer               :: c0,ie

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

    do ie = 1,nelemd
       elem(ie)%accum%Qvar=0
       elem(ie)%accum%Qmass=0
       elem(ie)%accum%Q1mass=0
       elem(ie)%accum%KEner=0
       elem(ie)%accum%IEner=0
       elem(ie)%accum%PEner=0
    end do
  end subroutine prim_printstate_init
!=======================================================================================================! 

  subroutine prim_printstate(elem, tl,hybrid,hvcoord,nets,nete)

    use physical_constants,     only: dd_pi

    type(element_t),            intent(inout), target :: elem(:)
    type(TimeLevel_t),target,   intent(in) :: tl
    type(hybrid_t),             intent(in) :: hybrid
    type(hvcoord_t),            intent(in) :: hvcoord
    integer,                    intent(in) :: nets,nete

    ! Local variables...
    integer :: i,j,k,ie
    integer,parameter  :: type=ORDERED

    real (kind=real_kind)  :: Mass2,Mass
    real (kind=real_kind)  :: TOTE(diagtimes),KEner(diagtimes),PEner(diagtimes),IEner(diagtimes)
    real (kind=real_kind)  :: Qvar(qsize_d,diagtimes),Qmass(qsize_d,diagtimes),Q1mass(qsize_d)
    real (kind=real_kind),save  :: time0
    real (kind=real_kind),save  :: TOTE0=0,Qmass0(qsize_d)=0
    real (kind=real_kind)  :: I_div(nlev)
    real (kind=real_kind)  :: I_vort(nlev)

    real (kind=real_kind)  :: tmp(np,np,nets:nete)
    real (kind=real_kind),allocatable  :: tmp3d(:,:,:,:)
    real (kind=real_kind)  :: tmp1(nets:nete)
    real (kind=real_kind)  :: ps(np,np)
    real (kind=real_kind)  :: dp(np,np)
    real (kind=real_kind)  :: phi_i(np,np,nlevp)
    real (kind=real_kind)  :: dphi(np,np,nlev)
    real (kind=real_kind)  :: w_over_dz(np,np,nlev),w_over_dz_p
    real (kind=real_kind)  :: tdiag(np,np,nlev), muvalue(np,np,nlevp), tempext
    !    real (kind=real_kind)  :: E(np,np)
    integer                :: location(3)

    real (kind=real_kind) :: usum_local(nets:nete), vsum_local(nets:nete), tsum_local(nets:nete), thetasum_local(nets:nete),&
                             psmin_local(nets:nete),psmax_local(nets:nete),pssum_local(nets:nete), &
                             fusum_local(nets:nete), fvsum_local(nets:nete), ftsum_local(nets:nete), fqsum_local(nets:nete), &
                             wsum_local(nets:nete), phisum_local(nets:nete), dpsum_local(nets:nete)

    real (kind=real_kind) :: umin_p, vmin_p, tmin_p, qvmin_p(qsize_d),&
         psmin_p, dpmin_p, TSmin_p

    real (kind=real_kind) :: umax_p, vmax_p, tmax_p, qvmax_p(qsize_d),&
         psmax_p, dpmax_p, TSmax_p

    real (kind=real_kind) :: usum_p, vsum_p, tsum_p, qvsum_p(qsize_d),&
         pssum_p, dpsum_p, thetasum_p, wsum_p

    real (kind=real_kind) :: fusum_p, fvsum_p, ftsum_p, fqsum_p
    real (kind=real_kind) :: fumin_p, fvmin_p, ftmin_p, fqmin_p
    real (kind=real_kind) :: fumax_p, fvmax_p, ftmax_p, fqmax_p
    real (kind=real_kind) :: phimax_p, phimin_p, phisum_p

    real (kind=real_kind) :: umax_local(2), umin_local(2), vmax_local(2), vmin_local(2),&
                             wmax_local(2), wmin_local(2), thetamin_local(2), thetamax_local(2),&
                             tmin_local(2), tmax_local(2), phimin_local(2), phimax_local(2), &
                             fumin_local(2), fumax_local(2), fvmin_local(2), fvmax_local(2), &
                             ftmin_local(2), ftmax_local(2), fqmin_local(2), fqmax_local(2), &
                             dpmin_local(2), dpmax_local(2), mumax_local(2), mumin_local(2), &
                             w_over_dz_max_local(2)
    character(len=3)      :: which
 
    real(kind=real_kind) :: relvort
    real(kind=real_kind) :: v1, v2, vco(np,np,2,nlev)

    real (kind=real_kind) :: time, time2,time1, scale, dt, dt_f
    real (kind=real_kind) :: IEvert1,IEvert2,PEvert1,PEvert2
    real (kind=real_kind) :: T1,T2,S1,S2,P1,P2
    real (kind=real_kind) :: PEhorz1,PEhorz2
    real (kind=real_kind) :: KEH1,KEH2,KEV1,KEV2
    real (kind=real_kind) :: KEwH1,KEwH2,KEwH3,KEwV1,KEwV2
    real (kind=real_kind) :: ddt_tot,ddt_diss, ddt_diss_adj
    integer               :: n0, n0q
    integer               :: npts,n,q
    integer :: max_itercnt_g
    real (kind=real_kind) :: max_deltaerr_g, max_reserr_g

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
    muvalue  = 0
    ! dynamics timelevels
    n0=tl%n0
    call TimeLevel_Qdp(tl, qsplit, n0q) ! get n0 level into n0q


    dt=tstep*qsplit
    if (rsplit>0) dt = tstep*qsplit*rsplit  ! vertical REMAP timestep 
    dt_f=dt
    if (ftype==4) dt_f=tstep

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


    npts = np
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
          global_shared_buf(ie,1) = 0
          do k=1,nlev
             global_shared_buf(ie,1) = global_shared_buf(ie,1) + &
                  SUM(elem(ie)%spheremp*elem(ie)%state%Qdp(:,:,k,q,n0q))
          enddo
       enddo
       call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
       qvsum_p(q) = global_shared_sum(1)
    enddo


    ! compute TS min/max
    do ie=nets,nete
       call get_field(elem(ie),'temperature',tdiag,hvcoord,n0,n0q)
       psmin_local(ie) = MINVAL(tdiag(:,:,nlev))
       psmax_local(ie) = MAXVAL(tdiag(:,:,nlev))
    enddo
    TSmin_p = ParallelMin(psmin_local,hybrid)
    TSmax_p = ParallelMax(psmax_local,hybrid)




    !
    do ie=nets,nete
       call get_field(elem(ie),'temperature',tdiag,hvcoord,n0,n0q)
       call get_field_i(elem(ie),'mu_i',muvalue,hvcoord,n0)

       ! layer thickness
       call get_phi(elem(ie),dphi,phi_i,hvcoord,n0)
       do k=1,nlev
          dphi(:,:,k)=-(phi_i(:,:,k+1)-phi_i(:,:,k))
          w_over_dz(:,:,k)=&
             g*max(elem(ie)%state%w_i(:,:,k,n0)/dphi(:,:,k),&
             elem(ie)%state%w_i(:,:,k+1,n0)/dphi(:,:,k))
       enddo

       ! surface pressure
       tmp(:,:,ie)=hvcoord%hyai(1)*hvcoord%ps0 + sum(elem(ie)%state%dp3d(:,:,:,n0),3) 

       !======================================================  

       psmax_local(ie) = MAXVAL(tmp(:,:,ie))

       call extremumLevelHelper(tmax_local,tdiag,'max',logical(ie == nets),nlev)
       call extremumLevelHelper(mumax_local,muvalue,'max',logical(ie == nets),nlevp)
       call extremumLevelHelper(phimax_local,dphi,'max',logical(ie == nets),nlev)
       call extremumLevelHelper(w_over_dz_max_local,w_over_dz,'max',logical(ie == nets),nlev)

       !======================================================

       psmin_local(ie) = MINVAL(tmp(:,:,ie))

       call extremumLevelHelper(tmin_local,tdiag,'min',logical(ie == nets),nlev)
       call extremumLevelHelper(mumin_local,muvalue,'min',logical(ie == nets),nlevp)
       call extremumLevelHelper(phimin_local,dphi,'min',logical(ie == nets),nlev)

       !======================================================

       usum_local(ie)    = SUM(elem(ie)%state%v(:,:,1,:,n0))
       vsum_local(ie)    = SUM(elem(ie)%state%v(:,:,2,:,n0))
       Wsum_local(ie)    = SUM(elem(ie)%state%w_i(:,:,:,n0))
       thetasum_local(ie)= SUM(elem(ie)%state%vtheta_dp(:,:,:,n0))
       phisum_local(ie)  = SUM(dphi)
       Fusum_local(ie)   = SUM(elem(ie)%derived%FM(:,:,1,:))
       Fvsum_local(ie)   = SUM(elem(ie)%derived%FM(:,:,2,:))
       tsum_local(ie)    = SUM(tdiag)

       dpsum_local(ie)    = SUM(elem(ie)%state%dp3d(:,:,:,n0))

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
       global_shared_buf(ie,12)=thetasum_local(ie)


    end do

    which = 'u'
    call findExtremumWithLevel(elem,umax_local,which,'max',n0,hybrid,hvcoord,nets,nete,.false.)
    call findExtremumWithLevel(elem,umin_local,which,'min',n0,hybrid,hvcoord,nets,nete,.false.)

    which = 'v'
    call findExtremumWithLevel(elem,vmax_local,which,'max',n0,hybrid,hvcoord,nets,nete,.false.)
    call findExtremumWithLevel(elem,vmin_local,which,'min',n0,hybrid,hvcoord,nets,nete,.false.)

    which = 'vth'
    call findExtremumWithLevel(elem,thetamax_local,which,'max',n0,hybrid,hvcoord,nets,nete,.false.)
    call findExtremumWithLevel(elem,thetamin_local,which,'min',n0,hybrid,hvcoord,nets,nete,.false.)

    which = 't'
    call findExtremumWithLevel(elem,tmax_local,which,'max',n0,hybrid,hvcoord,nets,nete,.true.)
    call findExtremumWithLevel(elem,tmin_local,which,'min',n0,hybrid,hvcoord,nets,nete,.true.)

    which = 'mu'
    call findExtremumWithLevel(elem,mumax_local,which,'max',n0,hybrid,hvcoord,nets,nete,.true.)
    call findExtremumWithLevel(elem,mumin_local,which,'min',n0,hybrid,hvcoord,nets,nete,.true.)

    which = 'dp'
    call findExtremumWithLevel(elem,dpmax_local,which,'max',n0,hybrid,hvcoord,nets,nete,.false.)
    call findExtremumWithLevel(elem,dpmin_local,which,'min',n0,hybrid,hvcoord,nets,nete,.false.)

    which = 'dg'
    call findExtremumWithLevel(elem,phimax_local,which,'max',n0,hybrid,hvcoord,nets,nete,.true.)
    call findExtremumWithLevel(elem,phimin_local,which,'min',n0,hybrid,hvcoord,nets,nete,.true.)

    which = 'cfl'
    call findExtremumWithLevel(elem,w_over_dz_max_local,which,'max',n0,hybrid,hvcoord,nets,nete,.true.)

    which = 'fu'
    call findExtremumWithLevel(elem,fumax_local,which,'max',n0,hybrid,hvcoord,nets,nete,.false.)
    call findExtremumWithLevel(elem,fumin_local,which,'min',n0,hybrid,hvcoord,nets,nete,.false.)

    which = 'fv'
    call findExtremumWithLevel(elem,fvmax_local,which,'max',n0,hybrid,hvcoord,nets,nete,.false.)
    call findExtremumWithLevel(elem,fvmin_local,which,'min',n0,hybrid,hvcoord,nets,nete,.false.)

    which = 'ft'
    call findExtremumWithLevel(elem,ftmax_local,which,'max',n0,hybrid,hvcoord,nets,nete,.false.)
    call findExtremumWithLevel(elem,ftmin_local,which,'min',n0,hybrid,hvcoord,nets,nete,.false.)

    which = 'fq1'
    call findExtremumWithLevel(elem,fqmax_local,which,'max',n0,hybrid,hvcoord,nets,nete,.false.)
    call findExtremumWithLevel(elem,fqmin_local,which,'min',n0,hybrid,hvcoord,nets,nete,.false.)

    which = 'w_i'
    call findExtremumWithLevel(elem,wmax_local,which,'max',n0,hybrid,hvcoord,nets,nete,.false.)
    call findExtremumWithLevel(elem,wmin_local,which,'min',n0,hybrid,hvcoord,nets,nete,.false.)

    psmin_p = ParallelMin(psmin_local,hybrid)
    psmax_p = ParallelMax(psmax_local,hybrid)

    call wrap_repro_sum(nvars=12, comm=hybrid%par%comm)
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
    thetasum_p= global_shared_sum(12)


    scale=1/g                                  ! assume code is using Pa
    if (hvcoord%ps0 <  2000 ) scale=100*scale  ! code is using mb
    ! after scaling, Energy is in J/m**2,  Mass kg/m**2
    ! so rate of changes are W/m**2


    !   mass = integral( ps-p(top) ) = sum(dp3d)
    do ie=nets,nete
       tmp(:,:,ie)=sum(elem(ie)%state%dp3d(:,:,:,n0),3) 
#ifdef CAM
       ! CAM neglicts the small, constant p(top) term.  Add this
       ! term in to be consistent
       tmp(:,:,ie)=tmp(:,:,ie) + hvcoord%hyai(1)*hvcoord%ps0 
#endif
    enddo
    Mass2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    Mass = Mass2*scale

    if(hybrid%masterthread) then
       write(iulog,109) "u     = ",umin_local(1)," (",nint(umin_local(2)),")",&
                                   umax_local(1)," (",nint(umax_local(2)),")",usum_p
       write(iulog,109) "v     = ",vmin_local(1)," (",nint(vmin_local(2)),")",&
                                   vmax_local(1)," (",nint(vmax_local(2)),")",vsum_p
       if(.not. theta_hydrostatic_mode) then
         write(iulog,109) "w     = ",wmin_local(1)," (",nint(wmin_local(2)),")",&
                                     wmax_local(1)," (",nint(wmax_local(2)),")",wsum_p
       endif

       write(iulog,109) "vTh_dp= ",thetamin_local(1)," (",nint(thetamin_local(2)),")",&
                                   thetamax_local(1)," (",nint(thetamax_local(2)),")",thetasum_p

       write(iulog,109) "T     = ",tmin_local(1)," (",nint(tmin_local(2)),")",&
                                   tmax_local(1)," (",nint(tmax_local(2)),")",tsum_p
       write(iulog,109) "mu    = ",mumin_local(1)," (",nint(mumin_local(2)),")",&
                                   mumax_local(1)," (",nint(mumax_local(2)),")"

       write(iulog,109) "dz(m) = ",phimin_local(1)/g," (",nint(phimin_local(2)),")",&
                                   phimax_local(1)/g," (",nint(phimax_local(2)),")",phisum_p/g

       write(iulog,109) "dp    = ",dpmin_local(1)," (",nint(dpmin_local(2)),")",&
                                   dpmax_local(1)," (",nint(dpmax_local(2)),")",dpsum_p

       write(iulog,109) "fu  = ",fumin_local(1)," (",nint(fumin_local(2)),")",&
                                 fumax_local(1)," (",nint(fumax_local(2)),")",fusum_p
       write(iulog,109) "fv  = ",fvmin_local(1)," (",nint(fvmin_local(2)),")",&
                                 fvmax_local(1)," (",nint(fvmax_local(2)),")",fvsum_p
       write(iulog,109) "ft  = ",ftmin_local(1)," (",nint(ftmin_local(2)),")",&
                                 ftmax_local(1)," (",nint(ftmax_local(2)),")",ftsum_p
       write(iulog,109) "fq1 = ",fqmin_local(1)," (",nint(fqmin_local(2)),")",&
                                 fqmax_local(1)," (",nint(fqmax_local(2)),")",fqsum_p

       if (.not.theta_hydrostatic_mode) then
          write(iulog,110) "   min dz/w (CFL condition) = ",1/max(1d-12,w_over_dz_max_local(1)),&
               " (",nint(w_over_dz_max_local(2)),")"
       endif

       do q=1,qsize
          write(iulog,102) "qv(",q,")= ",qvmin_p(q), qvmax_p(q), qvsum_p(q)
       enddo
       write(iulog,100) "TBOT= ",TSmin_p,TSmax_p
       write(iulog,100) "ps= ",psmin_p,psmax_p,pssum_p
       write(iulog,'(a,E23.15,a,E23.15,a)') "      M = ",Mass,' kg/m^2',Mass2,'mb     '

    end if
 
100 format (A10,3(E23.15))
102 format (A4,I3,A3,3(E24.16))
108 format (A10,E23.15,A6,E23.15,A6,E23.15)
109 format (A10,E23.15,A2,I3,A1,E23.15,A2,I3,A1,E23.15)
110 format (A33,E23.15,A2,I3,A1)

    if ( test_case(1:10) == "baroclinic" ) then
       ! zeta does not need to be made continious, but we  
       allocate(tmp3d(np,np,nlev,nets:nete))
       call compute_zeta_C0(tmp3d,elem,hybrid,nets,nete,n0)
       tmp=tmp3d(:,:,nlev,:)**2
       relvort = SQRT(global_integral(elem, tmp,hybrid,npts,nets,nete))
       deallocate(tmp3d)
    endif

    if(hybrid%masterthread) then 
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
    do n=1,diagtimes
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
#if defined(ENERGY_DIAGNOSTICS) && !defined (HOMMEXX_ENABLE_GPU)
    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%KEu_horiz1
    enddo
    !if(hybrid%masterthread) print *,'KEH1'
    KEH1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEH1 = KEH1*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%KEu_horiz2
    enddo
    !if(hybrid%masterthread) print *,'KEH2'
    KEH2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEH2 = KEH2*scale
    
    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%KEu_vert1
    enddo
    !if(hybrid%masterthread) print *,'KEV1'
    KEV1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEV1 = KEV1*scale
    
    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%KEu_vert2
    enddo
    !if(hybrid%masterthread) print *,'KEV2'
    KEV2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEV2 = KEV2*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%KEw_horiz1
    enddo
    !if(hybrid%masterthread) print *,'KEwH1'
    KEwH1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEwH1 = KEwH1*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%KEw_horiz2
    enddo
    !if(hybrid%masterthread) print *,'KEwH2'
    KEwH2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEwH2 = KEwH2*scale
    
    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%KEw_horiz3
    enddo
    !if(hybrid%masterthread) print *,'KEwH3'
    KEwH3 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEwH3 = KEwH3*scale
    
    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%KEw_vert1
    enddo
    !if(hybrid%masterthread) print *,'KEwV1'
    KEwV1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEwV1 = KEwV1*scale

    do ie=nets,nete
      tmp(:,:,ie) = elem(ie)%accum%KEw_vert2
    enddo
    !if(hybrid%masterthread) print *,'KEwV2'
    KEwV2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    KEwV2 = KEwV2*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%IEvert1 
    enddo
    !if(hybrid%masterthread) print *,'IEvert1'
    IEvert1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    IEvert1 = IEvert1*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%IEvert2
    enddo
    !if(hybrid%masterthread) print *,'IEvert2'
    IEvert2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    IEvert2 = IEvert2*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%PEhoriz1 
    enddo
    !if(hybrid%masterthread) print *,'PEhoriz1'
    PEhorz1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    PEhorz1 = PEhorz1*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%PEhoriz2
    enddo
    !if(hybrid%masterthread) print *,'PEhoriz2'
    PEhorz2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    PEhorz2 = PEhorz2*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%PEvert1
    enddo
    !if(hybrid%masterthread) print *,'PEver1'
    PEvert1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    PEvert1 = PEvert1*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%PEvert2
    enddo
    !if(hybrid%masterthread) print *,'PEver2'
    PEvert2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    PEvert2 = PEvert2*scale
    
    !   KE->IE
    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%T01
    enddo
    !if(hybrid%masterthread) print *,'T01'
    T1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    T1 = T1*scale
    
    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%T2
    enddo
    !if(hybrid%masterthread) print *,'T2'
    T2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    T2 = T2*scale
    
    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%S1
    enddo
    !if(hybrid%masterthread) print *,'S1'
    S1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    S1 = S1*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%S2
    enddo
    !if(hybrid%masterthread) print *,'S2'
    S2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    S2 = S2*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%P1
    enddo
    !if(hybrid%masterthread) print *,'P1'
    P1 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    P1 = P1*scale

    do ie=nets,nete
       tmp(:,:,ie) = elem(ie)%accum%P2
    enddo
    !if(hybrid%masterthread) print *,'P2'
    P2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    P2 = P2*scale

#else
    T1=0; T2=0; S1=0; S2=0; P1=0; P2=0; 
    IEvert1=0; PEhorz1=0; PEhorz2=0; PEvert1=0; PEvert2=0; 
    KEH1=0; KEH2=0;  KEV1=0; KEV2=0;  KEwH1=0; KEwH2=0; KEwH3=0;  KEwV1=0; KEwV2=0; 
#endif
    max_itercnt_g=parallelmax(max_itercnt,hybrid)
    max_deltaerr_g=parallelmax(max_deltaerr,hybrid)
    max_reserr_g=parallelmax(max_reserr,hybrid)
    max_itercnt=0 ; max_deltaerr=0; max_reserr=0 ! reset max counters each statefreq output




    if(hybrid%masterthread) then 
       if(use_moisture)then
          if (qsize>=1) then
             write(iulog,'(a,E23.15,a,E23.15,a)') "    dry M = ",Mass-Q1mass(1),' kg/m^2'
          endif
       endif
       
       write(iulog,'(3a25)') "**DYNAMICS**        J/m^2","   W/m^2","W/m^2    "
       if (ftype==4) &
            write(iulog,*) "NOTE:ftype=4 so d/dt and diss diagnostics include effects of forcing"
#if defined(ENERGY_DIAGNOSTICS) && !defined (HOMMEXX_ENABLE_GPU)
       ! terms computed during prim_advance, if ENERGY_DIAGNOSTICS is enabled
       if (theta_hydrostatic_mode) then
          write(iulog,'(a,2e22.14)')'KEu h-adv,sum=0:',KEH1,KEH2
          write(iulog,'(a,2e22.14)')'KEu v-adv,sum=0:',KEV1,KEV2
          write(iulog,'(a,3e22.14)')'IE  v-adv,sum=0:',IEvert1,PEvert1
          write(iulog,'(a,2e22.14)')'KE->I+P,I+P->KE:',(T1+PEhorz2),(S1+PEhorz1)
          
          ddt_tot  =  (KEner(2)-KEner(1))/dt
          ddt_diss = ddt_tot -(T1+PEhorz2)
          write(iulog,'(a,3E22.14)') " KE,d/dt,diss:",KEner(2),ddt_tot,ddt_diss

          ddt_tot  =  (IEner(2)-IEner(1))/dt
          write(iulog,'(a,3E22.14)') " IE,d/dt     :",IEner(2),ddt_tot
          ddt_tot  =  (PEner(2)-PEner(1))/dt
          write(iulog,'(a,3E22.14)') " PE,d/dt     :",PEner(2),ddt_tot

          ddt_tot =  (PEner(2)+IEner(2)-PEner(1)-IEner(1))/dt
          ddt_diss = ddt_tot - (S1+PEhorz1)
          write(iulog,'(a,3E22.14)') "I+P,d/dt,diss:",PEner(2)+IEner(2),ddt_tot,ddt_diss
          ddt_tot = (TOTE(2)-TOTE(1))/dt
          write(iulog,'(a,3E22.14)') "  E,d/dt,diss:",TOTE(2),ddt_tot
       else
          write(iulog,'(a,2e22.14)')'KEu h-adv,sum=0:',KEH1,KEH2
          write(iulog,'(a,3e22.14)')'KEw h-adv,sum=0:',KEwH1+KEwH3,KEwH2
          write(iulog,'(a,2e22.14)')'KEu v-adv,sum=0:',KEV1,KEV2
          write(iulog,'(a,2e22.14)')'KEw v-adv,sum=0:',KEwV1,KEwV2
          write(iulog,'(a,2e22.14)')'PE h-adv, sum=0:',PEhorz1,PEhorz2
          write(iulog,'(a,2e22.14)')'PE v-adv, sum=0:',PEvert1,PEvert2
          write(iulog,'(a,3e22.14)')'IE v-adv, sum=0:',IEvert1,IEvert2
          write(iulog,'(a,2e22.14)')'KEw wvor, ~0:   ',KEwH3
          write(iulog,'(a,2e22.14)')'KE->PE, PE->KE :',P1,P2
          write(iulog,'(a,2e22.14)')'KE->IE, IE->KE :',T1+T2,S1+S2
          
          ddt_tot  =  (KEner(2)-KEner(1))/dt
          ddt_diss = ddt_tot -(T1+T2+P1) 
          write(iulog,'(a,3E22.14)') "KE,d/dt,diss:",KEner(2),ddt_tot,ddt_diss
          !ddt_diss_adj = ddt_tot -(T1+T2+P1+KEwH1+KEwH2)
          !write(iulog,'(a,3E22.14)') "KE diss(adj):",ddt_diss_adj
          
          ddt_tot =  (IEner(2)-IEner(1))/dt
          ddt_diss = ddt_tot - (S1+S2)
          write(iulog,'(a,3E22.14)') "IE,d/dt,diss:",IEner(2),ddt_tot,ddt_diss
          !ddt_diss_adj = ddt_tot - (S1+S2+IEvert1+IEvert2)
          !write(iulog,'(a,3E22.14)') "IE diss(adj):",ddt_diss_adj
          
          ddt_tot = (PEner(2)-PEner(1))/dt
          ddt_diss = ddt_tot - P2
          write(iulog,'(a,3E22.14)') "PE,d/dt,diss:",PEner(2),ddt_tot,ddt_diss
          ddt_tot = (TOTE(2)-TOTE(1))/dt
          !ddt_diss = ddt_tot - (KEwH1+KEwH2+IEvert1+IEvert2)
          write(iulog,'(a,3E22.14)') " E,d/dt,diss:",TOTE(2),ddt_tot!,ddt_diss
       endif
#else
       write(iulog,'(a,3E22.14)') "KE,d/dt      ",KEner(2),(KEner(2)-KEner(1))/dt
       write(iulog,'(a,3E22.14)') "IE,d/dt      ",IEner(2),(IEner(2)-IEner(1))/dt
       write(iulog,'(a,3E22.14)') "PE,d/dt      ",PEner(2),(PEner(2)-PEner(1))/dt
       ddt_tot = (TOTE(2)-TOTE(1))/dt
       write(iulog,'(a,3E22.14)') " E,dE/dt     ",TOTE(2),ddt_tot
#endif
       
       do q=1,qsize
          write(iulog,'(a,i3,a,E22.14,a,2E15.7)') "Q",q,",Q diss, dQ^2/dt:",Qmass(q,2)," kg/m^2",&
               (Qmass(q,2)-Qmass(q,1))/dt,(Qvar(q,2)-Qvar(q,1))/dt
       enddo


       ! changes due to viscosity were with tstep
       ! changes due to forcing depend on ftype
       write(iulog,'(a)') 'Change from dribbled phys tendencies, viscosity, remap:'
       write(iulog,'(a,3e15.7)') 'dKE/dt(W/m^2): ',(KEner(1)-KEner(3))/dt_f,&
            (KEner(6)-KEner(5))/tstep,(KEner(2)-KEner(4))/dt
       write(iulog,'(a,3e15.7)') 'dIE/dt(W/m^2): ',(IEner(1)-IEner(3))/dt_f,&
            (IEner(6)-IEner(5))/tstep,(IEner(2)-IEner(4))/dt
       write(iulog,'(a,3e15.7)') 'dPE/dt(W/m^2): ',(PEner(1)-PEner(3))/dt_f,&
            (PEner(6)-PEner(5))/tstep,(PEner(2)-PEner(4))/dt
       q=1
       if (qsize>0) write(iulog,'(a,2e15.7)') 'dQ1/dt(kg/sm^2)',(Qmass(q,1)-Qmass(q,3))/dt

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

      ! IMEX diagnostics
      write(iulog,'(a,I3,2E23.15)') 'IMEX max iterations, error:',max_itercnt_g,max_deltaerr_g,max_reserr_g
    endif
    
    
    
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
    use dimensions_mod, only : np, np, nlev,nlevp
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
    real (kind=real_kind) :: phi(np,np,nlev)
    real (kind=real_kind) :: phi_i(np,np,nlevp)  
    real (kind=real_kind) :: pnh(np,np,nlev)   ! nh nonhyrdo pressure
    real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp) 
    real (kind=real_kind) :: exner(np,np,nlev)  ! exner nh pressure
    real (kind=real_kind) :: pnh_i(np,np,nlevp)  ! pressure on intefaces


    integer:: tmp, t1_qdp   ! the time pointers for Qdp are not the same

    if (t_before_advance) then
       t1=tl%n0
       call TimeLevel_Qdp( tl, qsplit, t1_qdp) !get n0 level into t2_qdp 
    else
       t1=tl%np1
       call TimeLevel_Qdp(tl, qsplit, tmp, t1_qdp) !get np1 into t2_qdp
    endif
    do ie=nets,nete
       dpt1=elem(ie)%state%dp3d(:,:,:,t1)
       call pnh_and_exner_from_eos(hvcoord,elem(ie)%state%vtheta_dp(:,:,:,t1),dpt1,&
            elem(ie)%state%phinh_i(:,:,:,t1),pnh,exner,dpnh_dp_i,pnh_i,'prim_state_mod')
       call get_phi(elem(ie),phi,phi_i,hvcoord,t1)

   
       !   KE   .5 dp/dn U^2
       do k=1,nlev
          E = ( elem(ie)%state%v(:,:,1,k,t1)**2 +  &
               elem(ie)%state%v(:,:,2,k,t1)**2  )/2
          sumlk(:,:,k) = E*dpt1(:,:,k)

          E=(elem(ie)%state%w_i(:,:,k,t1)**2+elem(ie)%state%w_i(:,:,k+1,t1)**2)/4
          suml2k(:,:,k) = E*dpt1(:,:,k)
       enddo
       suml=0
       suml2=0
       do k=1,nlev
          suml(:,:) = suml(:,:) + sumlk(:,:,k)
          suml2(:,:) = suml2(:,:) + suml2k(:,:,k)
       enddo
       if ( theta_hydrostatic_mode) then
          elem(ie)%accum%KEner(:,:,n)=suml(:,:)
       else
          elem(ie)%accum%KEner(:,:,n)=suml(:,:)+ suml2(:,:)
       endif

    !   PE   dp/dn PHIs
       suml=0
       do k=1,nlev
          suml = suml + phi(:,:,k)*dpt1(:,:,k)
       enddo
       elem(ie)%accum%PEner(:,:,n)=suml(:,:)
       

    !  IE = c_p^* dp/deta T - pnh dphi/deta  + ptop phi_top
       suml=0
       suml2=0
       do k=1,nlev
          suml(:,:)=suml(:,:)+&
                Cp*elem(ie)%state%vtheta_dp(:,:,k,t1)*exner(:,:,k) 
          suml2(:,:) = suml2(:,:)+(phi_i(:,:,k+1)-phi_i(:,:,k))*pnh(:,:,k)
       enddo
       elem(ie)%accum%IEner(:,:,n)=suml(:,:) + suml2(:,:) +&
            pnh_i(:,:,1)* phi_i(:,:,1)

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


!helper routine to compute min/max for derived quantities
subroutine extremumLevelHelper(res,field,operation,first,klev)
   use kinds, only : real_kind
   use dimensions_mod, only : np, np, nlev
   implicit none
   real (kind=real_kind), intent(inout) :: res(1:2) ! extremum and level where it happened
   character(len=*),      intent(in)    :: operation
   logical,               intent(in)    :: first
   real (kind=real_kind), intent(in)    :: field(np,np,klev)
   integer,               intent(in)    :: klev 

   real (kind=real_kind)                :: val
   integer                              :: location(3)

   if((operation /= 'max').and.(operation /= 'min')) call abortmp('unknown operation in extremumLevelHelper()')

   if ( first ) then
      if( operation == 'max' ) then
         res(1) = MAXVAL(field(:,:,1:klev))
         location = MAXLOC(field(:,:,1:klev))
         res(2) = location(3)
      else
         res(1) = MINVAL(field(:,:,1:klev))
         location = MINLOC(field(:,:,1:klev))
         res(2) = location(3)
      endif
   else
      if( operation == 'max' ) then
         val = MAXVAL(field(:,:,1:klev))
         if ( val > res(1) ) then
            res(1) = val
            location = MAXLOC(field(:,:,1:klev))
            res(2) = location(3)
         endif
      else
         val = MINVAL(field(:,:,1:klev))
         if ( val < res(1) ) then
            res(1) = val
            location = MINLOC(field(:,:,1:klev))
            res(2) = location(3)
         endif
      endif
   endif
end subroutine extremumLevelHelper


!doing extrema with level for all elems
!if processed, res already contains each thread's min/max value with level
!if not processed, compute min/max and level
subroutine findExtremumWithLevel(elem,res,which,operation,n0,hybrid,hvcoord,nets,nete,processed)
    use kinds, only : real_kind
    use dimensions_mod, only : np, np, nlev, nlevp
    implicit none
    real (kind=real_kind), intent(inout) :: res(1:2) ! extremum and level where it happened
    character(len=*),      intent(in)    :: which, operation
    integer,               intent(in)    :: nets,nete,n0
    type (element_t),      intent(in), target :: elem(:)
    type(hybrid_t),        intent(in)    :: hybrid
    type(hvcoord_t),       intent(in)    :: hvcoord
    logical,               intent(in)    :: processed

    integer                              :: ie,ksize
    real (kind=real_kind)                :: val
    integer                              :: location(3)
    real (kind=real_kind)                :: field(np,np,nlevp)

    if((operation /= 'max').and.(operation /= 'min')) call abortmp('unknown operation in findExtremumWithLevel()')

    if ( .not. processed ) then

       res(2) = -1
       if( operation == 'min' ) res(1) =  BIGVAL
       if( operation == 'max' ) res(1) = -BIGVAL
          
       !decide size of the column
       ksize=-1
       if( which == 'u' .or.  which == 'v' .or. which == 'vth' .or. which == 'dp' ) ksize=nlev
       if( which == 'fu' .or.  which == 'fv' .or. which == 'ft' .or. which == 'fq1' ) ksize=nlev
       if( which == 'w_i' ) ksize=nlevp

       if( ksize < 1) call abortmp('unset ksize in routine findExtremumWithLevel()')
   
       ! find max or min
       do ie=nets,nete

          if( which == 'u' )then
             field(1:np,1:np,1:ksize) = elem(ie)%state%v(1:np,1:np,1,1:ksize,n0)
          elseif( which == 'v' )then
             field(1:np,1:np,1:ksize) = elem(ie)%state%v(1:np,1:np,2,1:ksize,n0)
          elseif( which == 'vth' )then
             field(1:np,1:np,1:ksize) = elem(ie)%state%vtheta_dp(1:np,1:np,1:ksize,n0)
          elseif( which == 'w_i' )then
             field(1:np,1:np,1:ksize) = elem(ie)%state%w_i(1:np,1:np,1:ksize,n0)
          elseif( which == 'dp' )then
             field(1:np,1:np,1:ksize) = elem(ie)%state%dp3d(1:np,1:np,1:ksize,n0)
          elseif( which == 'fu' )then
             field(1:np,1:np,1:ksize) = elem(ie)%derived%FM(1:np,1:np,1,1:ksize)
          elseif( which == 'fv' )then
             field(1:np,1:np,1:ksize) = elem(ie)%derived%FM(1:np,1:np,2,1:ksize)
          elseif( which == 'ft' )then
             field(1:np,1:np,1:ksize) = elem(ie)%derived%FT(1:np,1:np,1:ksize)
          elseif( which == 'fq1' )then
             field(1:np,1:np,1:ksize) = elem(ie)%derived%FQ(1:np,1:np,1:ksize,1)
          endif

          if( operation == 'min' )then
             val = MINVAL(field(:,:,1:ksize))
             if( val < res(1) ) then
                location = MINLOC(field(:,:,1:ksize))
                res(2)   = location(3)
                res(1)   = val
             endif 
          elseif ( operation == 'max' )then
             val = MAXVAL(field(:,:,1:ksize))
             if( val > res(1) ) then
                location = MAXLOC(field(:,:,1:ksize))
                res(2)   = location(3)
                res(1)   = val          
             endif
          endif
       enddo !ie 

    endif ! processed

    if( operation == 'max' )then
       call ParallelMaxWithIndex(res, hybrid)
    else
       call ParallelMinWithIndex(res, hybrid)
    endif

end subroutine findExtremumWithLevel

end module prim_state_mod
