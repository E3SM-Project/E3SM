#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_state_mod
  ! ------------------------------
  use kinds, only : real_kind, iulog
  ! ------------------------------
  use dimensions_mod, only : nlev, np, nc, qsize_d, qsize, nelemd, ntrac, ntrac_d
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
  use time_mod, only : tstep, secpday, timelevel_t, TimeLevel_Qdp, time_at
  ! ------------------------------
  use control_mod, only : integration, test_case, runtype, moisture, &
       tstep_type,energy_fixer, qsplit, ftype, use_cpstar, rsplit
  ! ------------------------------
  use hybvcoord_mod, only : hvcoord_t 
  ! ------------------------------
  use global_norms_mod, only : global_integral, linf_snorm, l1_snorm, l2_snorm
  ! ------------------------------
  use element_mod, only : element_t
  ! ------------------------------
  use fvm_control_volume_mod, only : fvm_struct
  use spelt_mod, only : spelt_struct
  ! ------------------------------
  use viscosity_mod, only : compute_zeta_C0
  ! ------------------------------
  use reduction_mod, only : parallelmax,parallelmin
  ! ------------------------------
#ifdef _REFSOLN
  use ref_state_mod, only : ref_state_read, ref_state_write
#endif

implicit none
private
  character(len=*), private, parameter :: massfname = "mass.out"

  integer, public :: naccum


  public :: prim_printstate
  public :: prim_printstate_par
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

  subroutine prim_printstate(elem, tl,hybrid,hvcoord,nets,nete, fvm)
    use physical_constants, only : dd_pi
    use control_mod, only : tracer_transport_type
    use control_mod, only : TRACERTRANSPORT_LAGRANGIAN_FVM, TRACERTRANSPORT_FLUXFORM_FVM, TRACERTRANSPORT_SE_GLL

    type (element_t), intent(in) :: elem(:)
    
#if defined(_SPELT)
      type(spelt_struct), optional, intent(in) :: fvm(:)
#else
      type(fvm_struct), optional, intent(in) :: fvm(:)
#endif
    type (TimeLevel_t), target, intent(in) :: tl
    type (hybrid_t),intent(in)     :: hybrid
    type (hvcoord_t), intent(in)   :: hvcoord
    integer,intent(in)             :: nets,nete

    ! Local variables...
    integer :: i,j,k,ie
    integer,parameter  :: type=ORDERED

    real (kind=real_kind)  :: Mass2,Mass
    real (kind=real_kind)  :: TOTE(4),KEner(4),PEner(4),IEner(4),IEner_wet(4)
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
    !    real (kind=real_kind)  :: E(np,np)

    real (kind=real_kind) :: umin_local(nets:nete), umax_local(nets:nete), usum_local(nets:nete), &
         vmin_local(nets:nete), vmax_local(nets:nete), vsum_local(nets:nete), &
         tmin_local(nets:nete), tmax_local(nets:nete), tsum_local(nets:nete), &
         psmin_local(nets:nete),psmax_local(nets:nete),pssum_local(nets:nete), &
         fumin_local(nets:nete),fumax_local(nets:nete),fusum_local(nets:nete), &
         fvmin_local(nets:nete),fvmax_local(nets:nete),fvsum_local(nets:nete), &
         ftmin_local(nets:nete),ftmax_local(nets:nete),ftsum_local(nets:nete), &
         fqmin_local(nets:nete),fqmax_local(nets:nete),fqsum_local(nets:nete), &
         omegamin_local(nets:nete),omegamax_local(nets:nete),omegasum_local(nets:nete),&
         dpmin_local(nets:nete), dpmax_local(nets:nete), dpsum_local(nets:nete)


    real (kind=real_kind) :: umin_p, vmin_p, tmin_p, qvmin_p(qsize_d), cmin(ntrac_d),&
         psmin_p, dpmin_p


    real (kind=real_kind) :: umax_p, vmax_p, tmax_p, qvmax_p(qsize_d), cmax(ntrac_d),&
         psmax_p, dpmax_p

    real (kind=real_kind) :: usum_p, vsum_p, tsum_p, qvsum_p(qsize_d), csum(ntrac_d),&
         pssum_p, dpsum_p

    !
    ! for fvm diagnostics
    !
    real (kind=real_kind) :: psc_mass, psc_min, psc_max,dp_fvm_mass, dp_fvm_min, dp_fvm_max 


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
    
    if (hybrid%masterthread) then 
       write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
    end if
    if (.not. present(fvm) .and. ntrac>0) then
       print *,'ERROR: prim_state_mod.F90: optional fvm argument required if ntrac>0'
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
    npts=SIZE(elem(1)%state%lnps(:,:,n0),1)

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

       tmp(:,:,ie)=elem(ie)%state%ps_v(:,:,n0)
       !       tmp(:,:,ie)=EXP(elem(ie)%state%lnps(:,:,n0))


       !======================================================  
       umax_local(ie)    = MAXVAL(elem(ie)%state%v(:,:,1,:,n0))
       vmax_local(ie)    = MAXVAL(elem(ie)%state%v(:,:,2,:,n0))

       fumax_local(ie)    = MAXVAL(elem(ie)%derived%FM(:,:,1,:,pnm1))
       fvmax_local(ie)    = MAXVAL(elem(ie)%derived%FM(:,:,2,:,pnm1))

       tmax_local(ie)    = MAXVAL(elem(ie)%state%T(:,:,:,n0))

       if (rsplit>0) &
            dpmax_local(ie)    = MAXVAL(elem(ie)%state%dp3d(:,:,:,n0))

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

       if (rsplit>0) &
            dpmin_local(ie)    = MINVAL(elem(ie)%state%dp3d(:,:,:,n0))

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
       if (rsplit>0) then
          dpsum_local(ie)    = SUM(elem(ie)%state%dp3d(:,:,:,n0))
       else
          ! Make sure to initialize this to prevent possible
          ! floating point exceptions.
          dpsum_local(ie)    = 0.0D0
       end if

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
       global_shared_buf(ie,10) = dpsum_local(ie)
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

    Omegamin_p = ParallelMin(Omegamin_local,hybrid)
    Omegamax_p = ParallelMax(Omegamax_local,hybrid)

    call wrap_repro_sum(nvars=10, comm=hybrid%par%comm)
    usum_p = global_shared_sum(1)
    vsum_p = global_shared_sum(2)
    tsum_p = global_shared_sum(3)
    pssum_p = global_shared_sum(4)
    FUsum_p = global_shared_sum(5)
    FVsum_p = global_shared_sum(6)
    FTsum_p = global_shared_sum(7)
    FQsum_p = global_shared_sum(8)
    Omegasum_p = global_shared_sum(9)
    dpsum_p = global_shared_sum(10)

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

    !
    ! fvm diagnostics
    !
    if (ntrac>0) then
    if (tracer_transport_type == TRACERTRANSPORT_FLUXFORM_FVM.or.&
        tracer_transport_type == TRACERTRANSPORT_LAGRANGIAN_FVM) then
       do q=1,ntrac
          do ie=nets,nete
             tmp1(ie) = MINVAL(fvm(ie)%c(1:nc,1:nc,:,q,n0))
          enddo
          cmin(q) = ParallelMin(tmp1,hybrid)
          do ie=nets,nete
             tmp1(ie) = MAXVAL(fvm(ie)%c(1:nc,1:nc,:,q,n0))
          enddo
          cmax(q) = ParallelMax(tmp1,hybrid)
          !
          ! compute total tracer mass
          !
          global_shared_buf(:,1) = 0.0D0
          do k=1,nlev
             do ie=nets,nete
                global_shared_buf(ie,1) = global_shared_buf(ie,1)+&
                     SUM(fvm(ie)%c(1:nc,1:nc,k,q,n0)*&
                         fvm(ie)%dp_fvm(1:nc,1:nc,k,n0)*&
                         fvm(ie)%area_sphere(1:nc,1:nc))
             end do
          enddo
          call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
          csum(q) = global_shared_sum(1)/(dble(nlev)*4.0D0*DD_PI)
       enddo
       !
       ! psC diagnostics
       !
       do ie=nets,nete
          tmp1(ie) = MINVAL(fvm(ie)%psc(1:nc,1:nc))
       enddo
       psc_min = ParallelMin(tmp1,hybrid)
       do ie=nets,nete
          tmp1(ie) = MAXVAL(fvm(ie)%psc(1:nc,1:nc))
       enddo
       !
       ! surface pressure mass implied by fvm
       !
       psc_max = ParallelMax(tmp1,hybrid)
       do ie=nets,nete
          global_shared_buf(ie,1) = SUM(fvm(ie)%psc(1:nc,1:nc)*fvm(ie)%area_sphere(1:nc,1:nc))
       enddo
       call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
       psc_mass = global_shared_sum(1)/(4.0D0*DD_PI)
       !
       ! dp_fvm
       !
       do ie=nets,nete
          tmp1(ie) = MINVAL(fvm(ie)%dp_fvm(1:nc,1:nc,:,n0))
       enddo
       dp_fvm_min = ParallelMin(tmp1,hybrid)
       do ie=nets,nete
          tmp1(ie) = MAXVAL(fvm(ie)%dp_fvm(1:nc,1:nc,:,n0))
       enddo
       dp_fvm_max = ParallelMax(tmp1,hybrid)
       
       global_shared_buf(:,1) = 0.0D0
       do k=1,nlev
          do ie=nets,nete
             global_shared_buf(ie,1) = global_shared_buf(ie,1)+&
                  SUM(fvm(ie)%dp_fvm(1:nc,1:nc,k,n0)*fvm(ie)%area_sphere(1:nc,1:nc))
          end do
       enddo
       call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
       dp_fvm_mass = global_shared_sum(1)/(4.0D0*DD_PI)
    end if
    end if


    if(hybrid%masterthread) then
       write(iulog,100) "u     = ",umin_p,umax_p,usum_p
       write(iulog,100) "v     = ",vmin_p,vmax_p,vsum_p
       write(iulog,100) "omega = ",omegamin_p,omegamax_p,omegasum_p
       write(iulog,100) "t     = ",tmin_p,tmax_p,tsum_p
       if (rsplit>0) &
       write(iulog,100) "dp    = ",dpmin_p,dpmax_p,dpsum_p

       if (tstep_type>0) then  !no longer support tracer advection with tstep_type = 0
          do q=1,qsize
             write(iulog,100) "qv= ",qvmin_p(q), qvmax_p(q), qvsum_p(q)
          enddo
       endif
       write(iulog,100) "ps= ",psmin_p,psmax_p,pssum_p
       write(iulog,'(a,E23.15,a,E23.15,a)') "      M = ",Mass,' kg/m^2',Mass2,' mb     '






       if(fumin_p.ne.fumax_p) write(iulog,100) "fu = ",fumin_p,fumax_p,fusum_p
       if(fvmin_p.ne.fvmax_p) write(iulog,100) "fv = ",fvmin_p,fvmax_p,fvsum_p
       if(ftmin_p.ne.ftmax_p) write(iulog,100) "ft = ",ftmin_p,ftmax_p,ftsum_p
       if(fqmin_p.ne.fqmax_p) write(iulog,100) "fq = ",fqmin_p, fqmax_p, fqsum_p
       !
       ! fvm diagnostics
       !
       if (ntrac>0) then
       if (tracer_transport_type == TRACERTRANSPORT_FLUXFORM_FVM.or.&
           tracer_transport_type == TRACERTRANSPORT_LAGRANGIAN_FVM) then
          write(iulog,'(A36)') "-----------------------------------"
          write(iulog,'(A36)') "fvm diagnostics                    "
          write(iulog,'(A36)') "-----------------------------------"
          do q=1,ntrac
             write(iulog,'(A36,I1,3(E23.15))')&
                  "#c,min(c  ), max(c  ), mass(c  ) = ",q,cmin(q), cmax(q), csum(q)
          enddo
          write(iulog,'(A37,3(E23.15))')&
                  "   min(dp_), max(dp_), mass(dp_) =  ",dp_fvm_min, dp_fvm_max, dp_fvm_mass
          write(iulog,'(A37,3(E23.15))')&
                  "   min(psC), max(psC), mass(psC) =  ",psc_min, psc_max, psC_mass          
          write(iulog,'(A36)') "                                   "

       end if
       end if
    endif
 

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
             write(iulog,'(a,E23.15,a,E23.15,a)') "    dry M = ",Mass-Q1mass(1),' kg/m^2'
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
          
          if (tstep_type>0) then  !no longer support tracer advection with tstep_type = 0
             do q=1,qsize
                write(iulog,'(a,i1,a,E22.14,a,2E15.7)') "Q",q,",Q diss, dQ^2/dt:",Qmass(q,2)," kg/m^2",&
                     (Qmass(q,2)-Qmass(q,1))/dt,(Qvar(q,2)-Qvar(q,1))/dt
             enddo
          endif

          ! LF code diagnostics
          if (tstep_type==0) then  ! leapfrog
             write(iulog,'(a)') 'Robert filter, Physics (except adjustments):'
             write(iulog,'(a,2e15.7)') 'dKE/dt(W/m^2): ',(KEner(4)-KEner(3))/dt,(KEner(3)-KEner(2))/dt
             write(iulog,'(a,2e15.7)') 'dIE/dt(W/m^2): ',(IEner(4)-IEner(3))/dt,(IEner(3)-IEner(2))/dt
             write(iulog,'(a,2e15.7)') 'dPE/dt(W/m^2): ',(PEner(4)-PEner(3))/dt,(PEner(3)-PEner(2))/dt
          else
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
          if (tstep_type>0) then  !no longer support tracer advection with tstep_type = 0
             do q=1,qsize
                if(Qmass0(q)>0.0D0) then
                   write(iulog,'(a,E23.15,a,i1)') "(Q-Q0)/Q0 ",(Qmass(q,2)-Qmass0(q))/Qmass0(q),"   Q",q
                end if
             enddo
          endif
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
   
   
  subroutine prim_printstate_par(elem, tl,hybrid,hvcoord,nets,nete, par)
    type (element_t), intent(in) :: elem(:)
    type (TimeLevel_t), target, intent(in) :: tl
    type (hybrid_t),intent(in)     :: hybrid
    type (hvcoord_t), intent(in)   :: hvcoord
    integer,intent(in)             :: nets,nete
    character(len=*), parameter    :: fstub = "state_norms"
    integer	                   :: simday
    type(parallel_t)               :: par

    real (kind=real_kind)  :: v(np,np,2,nlev,nets:nete)
    real (kind=real_kind)  :: t(np,np,nlev,nets:nete)
    real (kind=real_kind)  :: ps_v(np,np,nets:nete)
    real (kind=real_kind)  :: vp(np,np,2,nlev,nets:nete)
    real (kind=real_kind)  :: tp(np,np,nlev,nets:nete)
    real (kind=real_kind)  :: ps_vp(np,np,nets:nete)
    real (kind=real_kind) :: l1,l2,linf
    integer               :: n0,i,j,k,ie,npts

    npts=SIZE(elem(1)%state%lnps(:,:,n0),1)
    n0=tl%n0
    do ie=nets,nete
       v(:,:,:,:,ie)=elem(ie)%state%v(:,:,:,:,n0) 
       T(:,:,:,ie)=elem(ie)%state%T(:,:,:,n0) 
       ps_v(:,:,ie)=elem(ie)%state%ps_v(:,:,n0) 
    enddo
       simday = 0

#ifdef _REFSOLN
! parallel write file with state vector in unformatted blocks for later calculation of norms
!    call ref_state_write(v(:,:,:,:,nets:nete),T(:,:,:,nets:nete),ps_v(:,:,nets:nete), & 
!	fstub,simday,nets,nete,par)
!    do ie=nets,nete
!       vp(:,:,:,:,ie)=v(:,:,:,:,ie)
!       Tp(:,:,:,ie)=T(:,:,:,ie)
!       ps_vp(:,:,ie)=ps_v(:,:,ie)
!    end do
#endif

#ifdef _REFSOLN
! parallel read file with state vector in unformatted blocks as written above
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
!  Parallel version of ref_state, comment out if writing above
!    call ref_state_read(vp(:,:,:,:,nets:nete),Tp(:,:,:,nets:nete),ps_vp(:,:,nets:nete), & 
!	fstub,simday,nets,nete,par)
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif

    npts=np

    l1   = l1_snorm(elem,ps_v(:,:,nets:nete),  ps_vp(:,:,nets:nete),hybrid,npts,nets,nete)
    l2   = l2_snorm(elem,ps_v(:,:,nets:nete),  ps_vp(:,:,nets:nete),hybrid,npts,nets,nete)
    linf = linf_snorm(ps_v(:,:,nets:nete),ps_vp(:,:,nets:nete),hybrid,npts,nets,nete)

    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       print *,simday, "L1=",l1
       print *,simday, "L2=",l2
       print *,simday, "Linf=",linf
    end if
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
#endif


  end subroutine prim_printstate_par

!=======================================================================================================! 


subroutine prim_energy_halftimes(elem,hvcoord,tl,n,t_before_advance,nets,nete)
! 
!  called at the end of a timestep, before timelevel update.  Solution known at
!  dynamics:     nm1,  n0,  np1.  
!
!                                                           LF code            RK code
!  This routine is called 4 times:  n=1:    t1=nm1, t2=n0   start              after forcing, before timestep
!                                   n=2:    t1=n0, t2=np1   after timestep     after timestep, including remap
!                                   n=3:    t1=n0, t2=np1   after forcing      after fixer 
!                                   n=4:    t1=n0, t2=np1   after Robert       before forcing
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
!
    use kinds, only : real_kind
    use dimensions_mod, only : np, np, nlev
    use control_mod, only : use_cpstar
    use hybvcoord_mod, only : hvcoord_t
    use element_mod, only : element_t
    use physical_constants, only : Cp, cpwater_vapor
    use physics_mod, only : Virtual_Specific_Heat, Virtual_Temperature
    use prim_si_mod, only : preq_hydrostatic

    integer :: t1,t2,n,nets,nete
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (TimeLevel_t), intent(in)       :: tl
    logical :: t_before_advance

    integer :: ie,k,i,j,nm_f
    real (kind=real_kind), dimension(np,np,nlev)  :: dpt1,dpt2   ! delta pressure
    real (kind=real_kind), dimension(np,np)  :: E
    real (kind=real_kind), dimension(np,np)  :: suml,suml2,v1,v2
    real (kind=real_kind), dimension(np,np,nlev)  :: sumlk, suml2k
    real (kind=real_kind), dimension(np,np,nlev)  :: p,T_v,phi
    real (kind=real_kind) :: cp_star1,cp_star2,qval_t1,qval_t2
    real (kind=real_kind) :: Qt
    logical :: wet


    logical tstagger
    integer:: t2_qdp, t1_qdp   ! the time pointers for Qdp are not the same

    nm_f = 1
    if (t_before_advance) then
       t1=tl%nm1
       t2=tl%n0
       call TimeLevel_Qdp( tl, qsplit, t2_qdp, t1_qdp) !get n0 level into t2_qdp 
    else
       t1=tl%n0
       t2=tl%np1
       call TimeLevel_Qdp(tl, qsplit, t1_qdp, t2_qdp) !get np1 into t2_qdp
    endif



! energy_fixer
!     <0         disabled, but compute energy non-staggered in time
!     0          disabled, but compute energy staggered in time (for use with leapfrog code)
!    >0          Enabled.  energy fixer requires dry formulation (use_cpstar=0) and non-staggered in time
!
    tstagger = .false.
    if (energy_fixer==0) tstagger = .true.


    !   IE   Cp*dpdn*T  + (Cpv-Cp) Qdpdn*T
    !        Cp*dpdn(n)*T(n+1) + (Cpv-Cp) Qdpdn(n)*T(n+1)
    !        [Cp + (Cpv-Cp) Q(n)] *dpdn(n)*T(n+1) 
    do ie=nets,nete

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
          dpt1(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,t1)
          dpt2(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,t2)
       enddo

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,cp_star1,cp_star2,qval_t1,qval_t2)
#endif
       do k=1,nlev
          sumlk(:,:,k)=0
          suml2k(:,:,k)=0
          do i=1,np
          do j=1,np
             if(use_cpstar == 1)  then
                ! Cp_star = cp + (Cpwater_vapor - cp)*qval
                qval_t1 = elem(ie)%state%Qdp(i,j,k,1,t1_qdp)/dpt1(i,j,k)
                qval_t2 = elem(ie)%state%Qdp(i,j,k,1,t2_qdp)/dpt2(i,j,k)
                cp_star1= Virtual_Specific_Heat(qval_t1)
                cp_star2= Virtual_Specific_Heat(qval_t2)
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
#if (defined COLUMN_OPENMP)
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



!      compute alternate PE term which matches what is used in CAM physics
       wet =(moisture /= "dry")
       do k=1,nlev
          p(:,:,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,t2)
          do j=1,np
             do i=1,np
                if (wet) then
                   Qt = elem(ie)%state%Qdp(i,j,k,1,t2_qdp)/dpt2(i,j,k)
                   T_v(i,j,k) = Virtual_Temperature(elem(ie)%state%T(i,j,k,t2),Qt)
                else
                   T_v(i,j,k) = elem(ie)%state%T(i,j,k,t2)
                endif
             end do
          end do
       end do

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
    integer:: t2_qdp, tmp   ! the time pointer for Qdp are not the same

    nm_f = 1
    if (t_before_advance) then
       t1=tl%nm1     
       t2=tl%n0
       call TimeLevel_Qdp( tl, qsplit, t2_qdp) !get n0 level into t2_qdp 
    else
       t1=tl%n0
       t2=tl%np1
       call TimeLevel_Qdp(tl, qsplit, tmp, t2_qdp) !get np1 into t2_qdp (don't need tmp)
    endif


    !
    !  RK2 forward scheme.  compute everything at t2
    !  (used by CAM)
    !   Q has only one time dimension
    if (tstep_type>0) then

       do ie=nets,nete
#if (defined COLUMN_OPENMP)
          !$omp parallel do private(q,k,suml)
#endif
          do q=1,qsize
             suml=0
             do k=1,nlev
                suml = suml + elem(ie)%state%Qdp(:,:,k,q,t2_qdp)*elem(ie)%state%Q(:,:,k,q)
             enddo
             elem(ie)%accum%Qvar(:,:,q,n)=suml(:,:)
          enddo
       enddo
       
       do ie=nets,nete
#if (defined COLUMN_OPENMP)
          !$omp parallel do private(q,k,suml)
#endif
          do q=1,qsize
             suml=0
             do k=1,nlev
                suml = suml + elem(ie)%state%Qdp(:,:,k,q,t2_qdp)
             enddo
             elem(ie)%accum%Q1mass(:,:,q)=suml(:,:)
             elem(ie)%accum%Qmass(:,:,q,n)=suml(:,:)
          enddo
       enddo
    endif



end subroutine prim_diag_scalars
end module prim_state_mod
