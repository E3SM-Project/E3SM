#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module flops_mod
  ! ---------------------
  use kinds, only : real_kind, iulog
  ! ---------------------
  use dimensions_mod, only : nv, np, nelem, nelemd, nlev, ne
  ! ---------------------
  use control_mod, only : filter_counter, filter_freq, precon_method, integration
  ! ---------------------
  use time_mod, only : nmax
  ! ---------------------
  use hybrid_mod, only : hybrid_t
  ! ---------------------
  use parallel_mod, only : iam, ncompoints, npackpoints, global_shared_buf, global_shared_sum
  ! ---------------------
  use global_norms_mod, only: wrap_repro_sum
  ! ---------------------
  use reduction_mod,only :  red_max, red_flops, red_timer, pmax_mt, pmin_mt
  ! ---------------------
  implicit none
  private

  real (kind=real_kind),public, parameter :: flops_sgrad=2.0*(nv+np)*nv*(2.0*np) + 2*nv*nv
  real (kind=real_kind),public, parameter :: flops_vgrad=2.0*nv*nv*(2.0*nv) + 2*nv*nv
  real (kind=real_kind),public, parameter :: flops_vort=2.0*nv*nv*(2.0*nv) + 3*nv*nv

  real (kind=real_kind),public, parameter :: flops_avg3= 3*4*(nv+1)
  real (kind=real_kind),public, parameter :: flops_avg2= 2*4*(nv+1)

  real (kind=real_kind),public, parameter :: flops_div=2.0*(nv+np)*np*(2.0*nv) + 3*np*np
  real (kind=real_kind),public, parameter :: flops_v2p=np*(np+nv)*(2.0*nv)

  real (kind=real_kind),public, parameter :: flops_exp_adv = (10.+21.+5.+14.)*nv*nv+(4.+7.)*np*np
  real (kind=real_kind),public, parameter :: flops_filter = 4.*np*np*np + 8.*nv*nv*nv + 4.*nv*nv
  real (kind=real_kind),public, parameter :: flops_exp = flops_sgrad + flops_vgrad + &
       flops_vort  + flops_avg3  + &
       flops_div   + flops_v2p   + &
       flops_exp_adv 

  ! =============================
  ! solver iteration independent
  ! semi implicit Flops count
  ! =============================

  real (kind=real_kind),public, parameter :: flops_sgrad_si = 3*flops_sgrad
  real (kind=real_kind),public, parameter :: flops_div_si   = 2*flops_div
  real (kind=real_kind),public, parameter :: flops_si_adv = (10.+23.+6.+16.)*nv*nv + (6.+7.)*np*np
  real (kind=real_kind),public, parameter :: flops_si = flops_sgrad_si + flops_vgrad + &
       flops_vort  + flops_avg3  + &
       flops_div_si   + flops_v2p   + &
       flops_avg2  + flops_si_adv

  ! =============================================
  ! Solver flop counts (per layer-iteration)
  ! =============================================

  real (kind=real_kind),public, parameter :: flops_cg0_it1  = (6.+3.)*np*np
  real (kind=real_kind),public, parameter :: flops_cg0_itn  = (4.+8.)*np*np
  real (kind=real_kind),public, parameter :: flops_lap_adv  = (8.+2.)*nv*nv
  real (kind=real_kind),public, parameter :: flops_blkjac   = 2.*(np*np)*(np*np)  ! matvec flops
  real (kind=real_kind),public, parameter :: flops_lap      = flops_sgrad + flops_div + &
       flops_avg2  + flops_lap_adv
  real (kind=real_kind),public, parameter :: flops_helm     = flops_lap + (4.)*np*np

  real (kind=real_kind),public            :: total_comm,total_adv,time_helm
  real (kind=real_kind),public            :: time_tmp

  real*8,public,parameter       :: timer_granularity=25.0
  integer,public                :: numlevels_exchanged   = (3+2)*nlev

  integer,public               :: min_number,max_number
  integer,public               :: min_message,max_message
  integer,public               :: TnSend

  real (kind=real_kind),public :: avg_message,avg_volume,avg_number
  real (kind=real_kind),public :: min_volume,max_volume
  real (kind=real_kind),public :: LocalComVolume,TotalComVolume


  public :: flops_report
  public :: prim_flops_report

contains

  subroutine flops_report(timer,hybrid,tot_iter)

    type(timer_t), intent(in) :: timer
    type (hybrid_t), intent(in) :: hybrid
    real*8, intent(in)          :: tot_iter

! make gtimer static to work around bug in g95 (0.90).  Otherwise, 
! associated(gtimer%time_buf)=.true. even before it is allocated.
    type(timer_t),save :: gtimer

    real*8                      :: tcount(2)

    real (kind=real_kind) :: flops_solver, flops_cg0, flops_tot
    real (kind=real_kind) :: flops_tmp(MAX_TIMER_CNT)
    logical, parameter :: Debug = .FALSE.

    if(Debug) write(iulog,*)'flops_report: point #1'
    flops_tmp = real(timer%time_buf,kind=real_kind)

    if(Debug) write(iulog,*)'flops_report: point #2'

    call pmax_mt(red_timer,flops_tmp,MAX_TIMER_CNT,hybrid)
    if(Debug) write(iulog,*)'flops_report: point #3'
    call timer_initialize(gtimer)

    gtimer%time_buf = red_timer%buf 

#ifdef _HTRACE
    call TRACE_FINALIZE(tcount)
    flops_tot = 1.0E6*tcount(2)

    ! ============================================================
    !  Perform the global sum across all threads and MPI processes
    ! ============================================================
    global_shared_buf(1,1) = flops_tot
    call wrap_repro_sum(nvars=1, comm=hybrid%par%comm, nsize=1)
    red_flops%buf(1) = global_shared_sum(1)
    flops_tot = red_flops%buf(1)
    if(Debug) write(iulog,*)'flops_report: point #5'
#else
    tcount(2) = 1.0E-6*flops_tot
#endif
    if(Debug) write(iulog,*)'flops_report: point #6'
    if (hybrid%par%masterproc .and. hybrid%ithr==0) then
       if (integration=="semi_imp") then
          if(Debug) write(iulog,*)'flops_report: point #6.1'
          flops_cg0   = nmax*flops_cg0_it1 + (tot_iter-2*nmax)*flops_cg0_itn 
          if(Debug) write(iulog,*)'flops_report: point #6.2'
          if (precon_method == "block_jacobi") then
             flops_solver = flops_cg0     + (tot_iter-nmax)*(flops_helm + flops_blkjac)
             if(Debug) write(iulog,*)'flops_report: point #6.3'
          else if (precon_method == "identity") then
             flops_solver = flops_cg0     + (tot_iter-nmax)*(flops_helm)
             if(Debug) write(iulog,*)'flops_report: point #6.4'
          end if
          if(Debug) write(iulog,*)'flops_report: point #7'

          ! assumes all levels converge in same iteration count

          flops_tot    = (flops_solver + nmax*flops_si + nmax*flops_filter/real(filter_freq,kind=real_kind))*(nelem*nlev)
          flops_solver = flops_solver*(nelem*nlev)
          flops_cg0    = flops_cg0*(nelem*nlev)
          if(Debug) write(iulog,*)'flops_report: point #8'

          write(iulog,9) tot_iter/nmax         

          write(iulog,300) gtimer%time_buf(T_init)
          write(iulog,301) gtimer%time_buf(T_topology)
          write(iulog,302) gtimer%time_buf(T_partition)
          write(iulog,303) gtimer%time_buf(T_metagraph)
          write(iulog,304) gtimer%time_buf(T_schedule)
          write(iulog,305) gtimer%time_buf(T_mass)
          write(iulog,306) gtimer%time_buf(T_topo_init),gtimer%time_buf(T_topo_init)/nelemd
          write(iulog,307) gtimer%time_buf(T_checksum)
          write(iulog,308) gtimer%time_buf(T_initrestart)

          !JMD  I don't understand why I am not using the global timer gtimer here 
          !JMD total_adv = gtimer%time_buf(T_adv - timer%solver - gtimer%time_buf(T_boundary - timer%boundary2 - timer%filter_boundary

          total_adv=gtimer%time_buf(T_adv)-gtimer%time_buf(T_solver)-gtimer%time_buf(T_boundary)&
               -gtimer%time_buf(T_boundary2)-gtimer%time_buf(T_filter_boundary)
          write(iulog,10) total_adv/nmax,(nmax*flops_si*nelem*nlev)/total_adv
          if(gtimer%time_buf(T_div) > 0.0) write(iulog,35) gtimer%time_buf(T_div)/nmax,&
               (nmax*flops_div_si*nelem*nlev)/gtimer%time_buf(T_div)
          if(gtimer%time_buf(T_sgrad) > 0.0) write(iulog,36) gtimer%time_buf(T_sgrad)&
               /nmax,(nmax*flops_sgrad_si*nelem*nlev)/gtimer%time_buf(T_sgrad)
          if(gtimer%time_buf(T_vgrad) > 0.0) write(iulog,37) gtimer%time_buf(T_vgrad)&
               /nmax,(nmax*flops_vgrad*nelem*nlev)/gtimer%time_buf(T_vgrad)
          if(gtimer%time_buf(T_vort) > 0.0) write(iulog,38) gtimer%time_buf(T_vort)&
               /nmax,(nmax*flops_vort*nelem*nlev)/gtimer%time_buf(T_vort)
          if(gtimer%time_buf(T_v2p) > 0.0) write(iulog,39) gtimer%time_buf(T_v2p)&
               /nmax,(nmax*flops_v2p*nelem*nlev)/gtimer%time_buf(T_v2p)
          if(gtimer%time_buf(T_filter) > 0.0) write(iulog,41) gtimer%time_buf(T_filter)/nmax, &
               (nmax*flops_filter*nelem*nlev)/(real(filter_freq,kind=real_kind)*gtimer%time_buf(T_filter))
          if(gtimer%time_buf(T_prepackcomp) > 0.0D0) write(iulog,100) gtimer%time_buf(T_prepackcomp)/nmax
          if(gtimer%time_buf(T_pack) > 0.0D0)   write(iulog,110)gtimer%time_buf(T_pack)/nmax,&
               (dble(nmax*nPackPoints*3*nlev*8)/gtimer%time_buf(T_pack))
          if(gtimer%time_buf(T_unpack) > 0.0D0) write(iulog,140)gtimer%time_buf(T_unpack)/nmax,&
               (dble(nmax*nPackPoints*3*nlev*8)/gtimer%time_buf(T_unpack))
          if(gtimer%time_buf(T_rotate) > 0.0) write(iulog,42) gtimer%time_buf(T_rotate)/nmax
          if(gtimer%time_buf(T_postunpackcomp) > 0.0D0) write(iulog,150) &
               gtimer%time_buf(T_postunpackcomp)/nmax
          if(gtimer%time_buf(T_si_adv) > 0.0) write(iulog,43) gtimer%time_buf(T_si_adv)/nmax,&
               (nmax*flops_si_adv*nelem*nlev)/gtimer%time_buf(T_si_adv)
          if(gtimer%time_buf(T_pointers) > 0.0) write(iulog,44) gtimer%time_buf(T_pointers)/nmax
          if(Debug) write(iulog,*)'flops_report: point #9'

          if(gtimer%time_buf(T_prepackcomp2) > 0.0D0) write(iulog,160) &
               gtimer%time_buf(T_prepackcomp2)/nmax
          if(gtimer%time_buf(T_pack2) > 0.0D0) write(iulog,170) gtimer%time_buf(T_pack2)/nmax,&
               (dble(nmax*nPackPoints*2*nlev*8)/gtimer%time_buf(T_pack2))
          if(gtimer%time_buf(T_unpack2) > 0.0D0)  write(iulog,200) gtimer%time_buf(T_unpack2)/nmax, &
               (dble(nmax*nPackPoints*2*nlev*8)/gtimer%time_buf(T_unpack2))
          if(gtimer%time_buf(T_postunpackcomp2) > 0.0D0) write(iulog,210) &
               gtimer%time_buf(T_postunpackcomp2)/nmax

          write(iulog,20)(gtimer%time_buf(T_solver)-gtimer%time_buf(T_reduce) -&
               gtimer%time_buf(T_lapcomm) - gtimer%time_buf(T_REDUCEBARRIER))/nmax, &
               flops_solver/(gtimer%time_buf(T_solver)-gtimer%time_buf(T_reduce)-gtimer%time_buf(T_lapcomm)&
               - gtimer%time_buf(T_REDUCEBARRIER)), &
               1.0E-6*flops_solver


          time_tmp = gtimer%time_buf(T_cg)-gtimer%time_buf(T_reduce) - gtimer%time_buf(T_REDUCEBARRIER)
          if(time_tmp > 0.0D0) write(iulog,30)time_tmp/nmax,flops_cg0/(time_tmp)

          if(Debug) write(iulog,*)'flops_report: point #10'

          if (precon_method == "block_jacobi") then
             if(gtimer%time_buf(T_precon) > 0.0D0) write(iulog,240)gtimer%time_buf(T_precon)&
                  /nmax,(tot_iter*flops_blkjac*nelem*nlev)/gtimer%time_buf(T_precon)
          else if (precon_method == "identity") then
             if(gtimer%time_buf(T_precon) > 0.0D0) write(iulog,240)gtimer%time_buf(T_precon)/nmax,0./gtimer%time_buf(T_precon)
          end if

          time_helm = gtimer%time_buf(T_helm) - gtimer%time_buf(T_lapcomm) - gtimer%time_buf(T_precon)
          if(time_helm > 0.0D0) write(iulog,50)time_helm/nmax,(tot_iter*flops_helm*nelem*nlev)/time_helm

          if(Debug) write(iulog,*)'flops_report: point #11'
          total_comm = gtimer%time_buf(T_lapcomm) + gtimer%time_buf(T_reduce) +gtimer%time_buf(T_boundary) +  &
               gtimer%time_buf(T_boundary2) + gtimer%time_buf(T_filter_boundary) + gtimer%time_buf(T_REDUCEBARRIER)
          if(total_comm > 0.0D0) write(iulog,6)total_comm/nmax, &
               (dble(nComPoints*(tot_iter*2+nmax*(2+3))*nlev*8))/total_comm

          if(gtimer%time_buf(T_filter_boundary) > 0.0D0) write(iulog,220) gtimer%time_buf(T_filter_boundary)/nmax
          if(gtimer%time_buf(T_send_time) > 0.0) write(iulog,91)gtimer%time_buf(T_send_time)/nmax
          if(gtimer%time_buf(T_recv_time) > 0.0) write(iulog,92)gtimer%time_buf(T_recv_time)/nmax
          if(gtimer%time_buf(T_wait_time) > 0.0) write(iulog,93)gtimer%time_buf(T_wait_time)/nmax
          if(gtimer%time_buf(T_copy_time) > 0.0) write(iulog,94)gtimer%time_buf(T_copy_time)/nmax

          if(gtimer%time_buf(T_boundary) > 0.0D0) write(iulog,130) gtimer%time_buf(T_boundary)/nmax, &
               dble(nmax*nComPoints*3*nlev*8)/gtimer%time_buf(T_boundary)

          if(gtimer%time_buf(T_boundary2) > 0.0D0) write(iulog,190) gtimer%time_buf(T_boundary2)/nmax, &
               dble(nmax*nComPoints*2*nlev*8)/gtimer%time_buf(T_boundary2)
          if(gtimer%time_buf(T_lapcomm) > 0.0D0 ) write(iulog,55) gtimer%time_buf(T_lapcomm)/nmax, &
               dble(tot_iter*nComPoints*2*nlev*8)/gtimer%time_buf(T_lapcomm)
          if(gtimer%time_buf(T_reduce) > 0.0D0) then 
             write(iulog,131)gtimer%time_buf(T_REDUCEBARRIER)/nmax
             write(iulog,135)gtimer%time_buf(T_reduce)/nmax
             write(iulog,136)gtimer%time_buf(T_reduce_smp)/nmax
             write(iulog,137)gtimer%time_buf(T_reduce_mpi)/nmax
          endif
          if(Debug) write(iulog,*)'flops_report: point #12'
#ifdef _HTRACE
          call TRACE_FINALIZE(tcount)
#else
          tcount(2) = 1.0E-6*flops_tot
#endif
          write(iulog,60)hybrid%par%nprocs,hybrid%Nthreads,np,nv,ne,nlev,gtimer%time_buf(T_adv)/nmax,  &
               flops_tot/gtimer%time_buf(T_adv),1.0E-6*flops_tot/dble(nmax),tcount(2)/dble(nmax)
          write(iulog,*) 'nmax,filter_counter ',nmax,filter_counter
       else
          ! =======================================
          ! This is for the explicit time-stepping
          ! =======================================

#ifdef _HTRACE

#else
          flops_tot = nmax*(flops_exp+flops_filter/real(filter_freq,kind=real_kind))*(nelem*nlev)
#endif

          write(iulog,300) gtimer%time_buf(T_init)
          write(iulog,301) gtimer%time_buf(T_topology)
          write(iulog,302) gtimer%time_buf(T_partition)
          write(iulog,303) gtimer%time_buf(T_metagraph)
          write(iulog,304) gtimer%time_buf(T_schedule)
          write(iulog,305) gtimer%time_buf(T_mass)
          write(iulog,306) gtimer%time_buf(T_topo_init),gtimer%time_buf(T_topo_init)/nelemd
          write(iulog,307) gtimer%time_buf(T_checksum)
          write(iulog,308) gtimer%time_buf(T_initrestart)
          total_comm = gtimer%time_buf(T_boundary) + gtimer%time_buf(T_filter_boundary)
          total_adv = gtimer%time_buf(T_adv) - total_comm
          if(Debug) write(iulog,*)'flops_report: point #13'
          if(total_adv > 0.0D0) &
               write(iulog,4) (total_adv)/nmax,flops_tot/(total_adv)
          if(gtimer%time_buf(T_div) > 0.0) write(iulog,35) gtimer%time_buf(T_div)/nmax,&
               (nmax*flops_div*nelem*nlev)/gtimer%time_buf(T_div)
          if(gtimer%time_buf(T_sgrad) > 0.0) write(iulog,36) gtimer%time_buf(T_sgrad)/nmax,&
               (nmax*flops_sgrad*nelem*nlev)/gtimer%time_buf(T_sgrad)
          if(gtimer%time_buf(T_vgrad) > 0.0) write(iulog,37) gtimer%time_buf(T_vgrad)/nmax,&
               (nmax*flops_vgrad*nelem*nlev)/gtimer%time_buf(T_vgrad)
          if(gtimer%time_buf(T_vort) > 0.0) write(iulog,38) gtimer%time_buf(T_vort)/nmax,&
               (nmax*flops_vort*nelem*nlev)/gtimer%time_buf(T_vort)
          if(gtimer%time_buf(T_v2p) > 0.0) write(iulog,39) gtimer%time_buf(T_v2p)/nmax,&
               (nmax*flops_v2p*nelem*nlev)/gtimer%time_buf(T_v2p)
          if(gtimer%time_buf(T_exp_adv) > 0.0) write(iulog,40) gtimer%time_buf(T_exp_adv)/nmax, &
               (nmax*flops_exp_adv*nelem*nlev)/gtimer%time_buf(T_exp_adv)
          if(gtimer%time_buf(T_filter) > 0.0) write(iulog,41) gtimer%time_buf(T_filter)/nmax, &
               (nmax*flops_filter*nelem*nlev)/(real(filter_freq,kind=real_kind)*gtimer%time_buf(T_filter))
          if(gtimer%time_buf(T_rotate) > 0.0) write(iulog,42) gtimer%time_buf(T_rotate)/nmax

          if(gtimer%time_buf(T_pack) > 0.0D0)   write(iulog,110) gtimer%time_buf(T_pack)/nmax,&
               dble(nmax*nPackPoints*3*nlev*8)/gtimer%time_buf(T_pack)
          if(gtimer%time_buf(T_unpack) > 0.0D0) write(iulog,140) gtimer%time_buf(T_unpack)/nmax,&
               dble(nmax*nPackPoints*3*nlev*8)/gtimer%time_buf(T_unpack)

          if(Debug) write(iulog,*)'flops_report: point #14'
          if(total_comm > 0.0) write(iulog,66) (total_comm)/nmax

          if(gtimer%time_buf(T_boundary) > 0.0)  &
               write(iulog,130) gtimer%time_buf(T_boundary)/nmax,dble(nmax)*dble(TotalComVolume)/gtimer%time_buf(T_boundary)
          if(gtimer%time_buf(T_filter_boundary) > 0.0D0) write(iulog,220) gtimer%time_buf(T_filter_boundary)/nmax
          if(gtimer%time_buf(T_send_time) > 0.0) write(iulog,91)gtimer%time_buf(T_send_time)/nmax
          if(gtimer%time_buf(T_recv_time) > 0.0) write(iulog,92)gtimer%time_buf(T_recv_time)/nmax
          if(gtimer%time_buf(T_wait_time) > 0.0) write(iulog,93)gtimer%time_buf(T_wait_time)/nmax
          if(gtimer%time_buf(T_copy_time) > 0.0) write(iulog,94)gtimer%time_buf(T_copy_time)/nmax

          write(iulog,5)hybrid%par%nprocs,hybrid%Nthreads,np,nv,ne,nlev,gtimer%time_buf(T_adv)/nmax, &
               flops_tot/gtimer%time_buf(T_adv),1.0E-6*flops_tot/dble(nmax),tcount(2)/dble(nmax)
          if(Debug) write(iulog,*)'flops_report: point #15'

       end if
    end if

    call timer_finalize(gtimer)

300 format("INITALIZATION:         (usec)= ",f11.3)
301 format("   TOPOLOGY            (usec)= ",f11.3)  
302 format("   PARTITION           (usec)= ",f11.3)  
303 format("   METAGRAPH           (usec)= ",f11.3)  
304 format("   SCHEDULE            (usec)= ",f11.3)  
305 format("   MASS                (usec)= ",f11.3)  
306 format("   TOPO_INIT           (usec)= ",f11.3," (usec/element) =",f11.3)  
307 format("   CHECKSUM            (usec)= ",f11.3)  
308 format("   INITRESTART         (usec)= ",f11.3)  
4   format("EXP FLOPS              (usec/step)=",f11.3," Mflops     =",f16.3)
35  format("    DIVERGENCE_V2P     (usec/step)=",f11.3," MFlops     =",f16.3)
36  format("    GRADIENT_WK        (usec/step)=",f11.3," MFlops     =",f16.3)
37  format("    GRADIENT_V2V       (usec/step)=",f11.3," MFlops     =",f16.3)
38  format("    VORTICITY          (usec/step)=",f11.3," MFlops     =",f16.3)
39  format("    INTERPOLATE_V2P    (usec/step)=",f11.3," MFlops     =",f16.3)
40  format("    ADVANCE            (usec/step)=",f11.3," MFlops     =",f16.3)
41  format("    FILTER             (usec/step)=",f11.3," Mflops     =",f16.3)
42  format("    ROTATE             (usec/step)=",f11.3)
43  format("    ADVANCE_SI         (usec/step)=",f11.3," MFlops     =",f16.3)
44  format("    POINTERS           (usec/step)=",f11.3)
6   format("COMMUNICATION TOTAL    (usec/step)=",f11.3," Mbytes/sec =",f11.3)
66  format("COMMUNICATION TOTAL    (usec/step)=",f11.3)
5   format("ND=",i5," THR=",i2," NP=",i2," NV=",i2," NE=",i2," NLEV=",i2, &
         " EXP TOT  (usec/step)=",f11.3," Mflops = ",f16.3," Total (Mflop/step) = ",e16.8,e16.8)
91  format("              SEND     (usec/step)=",f11.3)
92  format("              RECV     (usec/step)=",f11.3)
93  format("              WAIT     (usec/step)=",f11.3)
94  format("              COPY     (usec/step)=",f11.3)

9   format("MEAN SOLVER ITERATIONS = ",f9.4)
10  format("SI STEP                (usec/step)=",f11.3," Mflops     =",f16.3)
100 format("   PRE PACK COMP    #1 (usec/step)=",f11.3)
110 format("   PACK             #1 (usec/step)=",f11.3," Mbytes/sec =",f11.3)
130 format("   ADV BNDRY EXCH   #1 (usec/step)=",f11.3," Mbytes/sec =",f11.3)
140 format("   UNPACK           #1 (usec/step)=",f11.3," Mbytes/sec =",f11.3)
150 format("   POST UNPACK COMP #1 (usec/step)=",f11.3)
160 format("   PRE PACK COMP    #2 (usec/step)=",f11.3)
170 format("   PACK             #2 (usec/step)=",f11.3," Mbytes/sec =",f11.3)
190 format("   ADV BNDRY EXCH   #2 (usec/step)=",f11.3," Mbytes/sec =",f11.3)
200 format("   UNPACK           #2 (usec/step)=",f11.3," Mbytes/sec =",f11.3)
210 format("   POST UNPACK COMP #2 (usec/step)=",f11.3)
220 format("   FILTER BNDRY EXCH   (usec/step)=",f11.3)

20  format("PCG SOLVER             (usec/step)=",f11.3," Mflops     =",f16.3," (Mflop/step) = ",e16.8)
30  format("   CONGRAD CORE        (usec/step)=",f11.3," Mflops     =",f16.3)
135 format("   CONGRAD REDUCE      (usec/step)=",f11.3)
131 format("   REDUCE BARRIER      (usec/step)=",f11.3)
136 format("           SMP         (usec/step)=",f11.3)
137 format("           MPI         (usec/step)=",f11.3)
240 format("   PRECONDITIONER      (usec/step)=",f11.3," Mflops     =",f16.3)
50  format("   HELMHOLTZ OPER      (usec/step)=",f11.3," Mflops     =",f16.3)
55  format("   HELM BNDRY EXCH     (usec/step)=",f11.3," Mbytes/sec =",f11.3)

60  format("ND=",i5," THR=",i2," NP=",i2," NV=",i2," NE=",i2," NLEV=",i2, &
         " SI  TOT  (usec/step)=",f11.3," Mflops =",f16.3," Total (Mflop/step) = ",e16.8,e16.8)
  end subroutine flops_report

  subroutine prim_flops_report(timer,hybrid,tot_iter)

    type (timer_t), intent(in)  :: timer
    !     real*8,          intent(in)  :: timer
    type (hybrid_t), intent(in) :: hybrid
    real*8, intent(in)          :: tot_iter
    real*8                      :: tcount(2)
    real*8                      :: tmax,tmin,prim_seam_min,prim_seam_adv_min,prim_seam_adv1_min
    real (kind=real_kind) :: flops_tot
    real (kind=real_kind) :: flops_exper
    real (kind=real_kind) :: flops_tmp(MAX_TIMER_CNT)
    real (kind=real_kind) :: avg_barrier_time,max_barrier_time,imbalance_time
    real (kind=real_kind) :: time_adv,time_boundary,time_filter_boundary
    real (kind=real_kind) :: total_adv,total_comm,total_comp,total_boundary, &
                             total_adv1,boundary_min,adv_per_step
    real (kind=real_kind) :: max_comm_size
! make gtimer static to work around bug in g95 (0.90).  Otherwise, 
! associated(gtimer%time_buf)=.true. even before it is allocated.
    type(timer_t),save         :: gtimer
    real (kind=real_kind) :: time_tmp,flop_rate

!#if 0
    !JMD     write(iulog,*)'ITHR: ',hybrid%ithr,' prim_flops_report (BEFORE):  ',timer(t_adv),timer(T_boundary)

    flops_tmp = real(timer%time_buf,kind=real_kind)

    call pmax_mt(red_timer,flops_tmp,MAX_TIMER_CNT,hybrid)

    call timer_initialize(gtimer)

    gtimer%time_buf = red_timer%buf

!TRV     flops_tmp(1)   = timer%time_buf(T_adv)
!     call pmax_mt(red_max,flops_tmp,1,hybrid)
!     time_adv       = red_max%buf(1)

!TRV     flops_tmp(1)   = timer%time_buf(T_boundary)
!     call pmax_mt(red_max,flops_tmp,1,hybrid)
!     time_boundary = red_max%buf(1)
    !JMD     write(iulog,*)'ITHR: ',hybrid%ithr,' prim_flops_report (AFTER):  ',time_adv,time_boundary

#ifdef _HTRACE
    call TRACE_FINALIZE(tcount)
    flops_tot = 1.0E6*tcount(2)

    ! ============================================================
    !  Perform the global sum across all threads and MPI processes
    ! ============================================================
    global_shared_buf(1,1) = flops_tot
    call wrap_repro_sum(nvars=1, comm=hybrid%par%comm, nsize=1)
    red_flops%buf(1) = global_shared_sum(1)
    flops_tot = red_flops%buf(1)
#else
    tcount(2) = 1.0E-6*flops_tot
#endif

!     flops_tmp(1)=timer%boundary2/nmax
!     call psum_mt(red_flops,flops_tmp,1,hybrid)
!     avg_barrier_time = red_flops%buf(1)/hybrid%par%nprocs

!     flops_tmp(1)=timer%boundary2/nmax
!     call pmax_mt(red_flops,flops_tmp,1,hybrid) 
!     max_barrier_time = red_flops%buf(1)

!TRV todo    if(avg_barrier_time > 0.0) imbalance_time = max_barrier_time/avg_barrier_time 


    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       if (integration=="semi_imp") then

#ifndef _HTRACE
          flops_exper = 0

          selectcase(ne)
          !----------------
          ! Resolution C28 
          !----------------
          case(18)     
             ! ----------------
             !  Levels 
             ! ----------------
             selectcase(nlev)
             ! L16 
             case(20)
                
                !---------------
                ! Filter frequency
                !---------------
                selectcase(filter_freq)
                case(1)
                   flops_exper = 5549.16
                case default 
                   call FlopCountError(ne,nlev,filter_freq)  
                end select  ! end filter frequency 

             case(40)

                selectcase(filter_freq)
                case(1)
                   flops_exper = 12554.31
                case default
                   call FlopCountError(ne,nlev,filter_freq)
                end select
             case default
                call FlopCountError(ne,nlev,filter_freq)  
             end select   ! end  NLEV 
          case(36)

             selectcase(nlev)
             case(40)
                selectcase(filter_freq)
                case(1)
                   flops_exper = 49394.19
                case default
                   call FlopCountError(ne,nlev,filter_freq)
                end select  ! end filter frequency
             case default
                call FlopCountError(ne,nlev,filter_freq)
             end select
          case default
             call FlopCountError(ne,nlev,filter_freq)  
          end select   ! end resolution 
          
          flops_tot = 1.0E6*nmax*flops_exper

#endif

          write(iulog,*)'semi-implicit timing output'
          ! timers in semi_imp:
          ! preq_filter(passed timer) T_FILTER,T_FILTER_BOUNDARY,T_PACK,T_ROTATE
          ! bndry_exchangev(passed timer) T_BOUNDARY,T_BOUNDARY2,T_SEND_TIME,T_SEND_TIME,
          !   T_RECV_TIME,T_WAIT_TIME,T_COPY_TIME

          ! --- timing calculations ---
          ! * computation timers are T_adv (sub-timers are T_filter,T_pack,T_unpack,
          !   T_rotate)
          ! * communication timers are T_boundary, T_boundary2, T_filter_boundary, 
          !   T_filter_boundary2

          time_boundary = gtimer%time_buf(T_boundary)
          time_filter_boundary = gtimer%time_buf(T_filter_boundary)/(nmax/filter_freq)
          total_comm = gtimer%time_buf(T_boundary) + gtimer%time_buf(T_filter_boundary) + &
               gtimer%time_buf(T_boundary2) + gtimer%time_buf(T_filter_boundary2)
          total_boundary = gtimer%time_buf(T_boundary2) + gtimer%time_buf(T_filter_boundary2)
          total_comp = gtimer%time_buf(T_adv) - total_comm
          adv_per_step = (gtimer%time_buf(T_adv) - total_comm)/nmax

          ! --- print the timing stats ---

          ! Initialization times
          write(iulog,300)gtimer%time_buf(T_init)
          write(iulog,301)gtimer%time_buf(T_topology)
          write(iulog,302)gtimer%time_buf(T_partition)
          write(iulog,303)gtimer%time_buf(T_metagraph)
          write(iulog,304)gtimer%time_buf(T_schedule)
          write(iulog,305)gtimer%time_buf(T_mass)
          write(iulog,306)gtimer%time_buf(T_topo_init),gtimer%time_buf(T_topo_init)/nelemd
          write(iulog,307)gtimer%time_buf(T_checksum)
          write(iulog,308)gtimer%time_buf(T_initrestart)

          write(*,309)gtimer%time_buf(T_main_loop)/(nmax-1)

          ! applycolumnmodel() time
          write(iulog,310)gtimer%time_buf(T_apply_column_model)/(nmax-1)

!           ! Total advance time
!           write(*,310)gtimer%time_buf(T_adv)/nmax, (flops_tot/(gtimer%time_buf(T_adv)))

          ! advance time 
          write(iulog,444)(gtimer%time_buf(T_adv))/nmax !, flops_tot/gtimer%time_buf(T_adv)

          ! prim tracer advection time
          write(iulog,44)gtimer%time_buf(T_ADVEC_TRACERS)/nmax

          ! timelevel_update()
          write(iulog,311)gtimer%time_buf(T_timelevel_update)/(nmax-1)

          ! movie output time
          write(iulog,313)gtimer%time_buf(T_movie_output)/(nmax-1)

          ! accum output
          write(iulog,314)gtimer%time_buf(T_accum_output)/(nmax-1)

          if(total_comm > 0.0)  & 
               write(iulog,6)(total_comm)/nmax,  (dble(nmax)*dble(TotalComVolume)/total_comm)
          ! advance and tracer advect boundary exchange time
          if(time_boundary > 0.0)  &
               write(iulog,131)(time_boundary)/nmax,(dble(nmax)*dble(TotalComVolume)/time_boundary)
          ! prim_filter boundary exchange time
          if(gtimer%time_buf(T_filter_boundary) > 0.0) &
               write(iulog,220)(gtimer%time_buf(T_filter_boundary)/nmax)
          if(timer%time_buf(T_send_time) > 0.0) write(iulog,91)timer%time_buf(T_send_time)/nmax
          if(timer%time_buf(T_recv_time) > 0.0) write(iulog,92)timer%time_buf(T_recv_time)/nmax
          if(timer%time_buf(T_wait_time) > 0.0) write(iulog,93)timer%time_buf(T_wait_time)/nmax
          if(timer%time_buf(T_copy_time) > 0.0) write(iulog,94)timer%time_buf(T_copy_time)/nmax

          write(iulog,51) hybrid%par%nprocs, hybrid%Nthreads, np, nv, ne, nlev, &
               gtimer%time_buf(T_main_loop)/(nmax-1), 0, 0

       ! ------------------------------
       ! --- Explicit timing output ---
       ! ------------------------------
       else

#ifndef _HTRACE

          selectcase(ne)
          !----------------
          ! Resolution C28 
          !----------------
          case(4)     
             ! ----------------
             !  Levels 
             ! ----------------
             selectcase(nlev)
             ! L16 
             case(16)
                
                !---------------
                ! Filter frequency
                !---------------
                selectcase(filter_freq)
                case(12)
                   flops_exper = 0.51287163E+02 
                case(-1)
                   flops_exper = 0.50795865E+02 
                case default 
                   call FlopCountError(ne,nlev,filter_freq)  
                end select  ! end filter frequency 

             case default
                call FlopCountError(ne,nlev,filter_freq)  
             end select   ! end  NLEV 
          !----------------
          ! Resolution C56 
          !----------------
          case(8)
             ! ----------------
             !  Levels 
             ! ----------------
             selectcase(nlev)
                
             ! L16
             case(16)
                !---------------
                ! Filter frequency
                !---------------
                selectcase(filter_freq)
                case(12)
                   ! flops_exper = 0.20497647E+03 
                   flops_exper = 0.18200433E+03
                case(6)
                   ! flops_exper = 0.24150208E+03 
                   flops_exper = 0.18469502E+03
                case(-1)
                   ! flops_exper = 0.34384912E+03
                   flops_exper = 0.17931515E+03
                case default 
                   call FlopCountError(ne,nlev,filter_freq)  
                end select  ! end filter frequency 
                
             case default
                call FlopCountError(ne,nlev,filter_freq)  
             end select  ! end NLEV
          !----------------
          ! Resolution C154
          !----------------
          case(22)
             ! ----------------
             !  Levels
             ! ----------------
             selectcase(nlev)
             ! L16
             case(16)
                !---------------
                ! Filter frequency
                !---------------
                selectcase(filter_freq)
                case(-1)
                   flops_exper = 0.158504E+04    ! Old value flops_exper = 0.17038378E+04
                case default
                   call FlopCountError(ne,nlev,filter_freq)
                end select  ! end filter frequency
             case default
                call FlopCountError(ne,nlev,filter_freq)
             end select   ! end NLEV
          case default
             call FlopCountError(ne,nlev,filter_freq)  
          end select   ! end resolution 
          
          flops_tot = 1.0E6*nmax*flops_exper
#endif

! --- timing calculations ---
! * computation timers are T_adv and T_advec_tracers 
!   * (sub-timers are T_filter,T_pack,T_unpack,T_rotate)
! * communication timers are T_boundary, T_boundary2, T_filter_boundary, T_filter_boundary2
!
          time_boundary = gtimer%time_buf(T_boundary)
          time_filter_boundary = gtimer%time_buf(T_filter_boundary)/(nmax/filter_freq)
          total_comm = gtimer%time_buf(T_boundary) + gtimer%time_buf(T_filter_boundary) + &
               gtimer%time_buf(T_boundary2) + gtimer%time_buf(T_filter_boundary2)
          total_boundary = gtimer%time_buf(T_boundary2) + gtimer%time_buf(T_filter_boundary2)
          total_comp = gtimer%time_buf(T_adv) + gtimer%time_buf(T_advec_tracers) - total_comm
          adv_per_step = (gtimer%time_buf(T_adv))/nmax

!print the timing stats

          ! Initialization times
          write(iulog,300)gtimer%time_buf(T_init)
          write(iulog,301)gtimer%time_buf(T_topology)
          write(iulog,302)gtimer%time_buf(T_partition)
          write(iulog,303)gtimer%time_buf(T_metagraph)
          write(iulog,304)gtimer%time_buf(T_schedule)
          write(iulog,305)gtimer%time_buf(T_mass)
          write(iulog,306)gtimer%time_buf(T_topo_init),gtimer%time_buf(T_topo_init)/nelemd
          write(iulog,307)gtimer%time_buf(T_checksum)
          write(iulog,308)gtimer%time_buf(T_initrestart)

          write(*,309)gtimer%time_buf(T_main_loop)/(nmax-1)
          ! Total main timestepping loop time (advance+tracer advection)
!           write(*,310)(gtimer%time_buf(T_adv)+gtimer%time_buf(T_advec_tracers))/nmax, &
!                flops_tot/(gtimer%time_buf(T_adv) + gtimer%time_buf(T_advec_tracers))

          !           (flops_tot/(gtimer%time_buf(T_adv) + gtimer%time_buf(T_adv1)))
               

          !       write(iulog,4)total_adv, flops_tot/(gtimer%time_buf(T_adv)-total_comm)

          ! applycolumnmodel() time
          write(iulog,310)gtimer%time_buf(T_apply_column_model)/(nmax-1)

          ! advance time 
          if(gtimer%time_buf(T_adv) > 0.0)  &
		write(iulog,4)(gtimer%time_buf(T_adv))/nmax, flops_tot/gtimer%time_buf(T_adv)

          ! advance sections
          if(TIMER_DETAIL(3,timer)) then 
             write(iulog,608)'prim_filter            (usec/step) = ',gtimer%time_buf(T_filter)/nmax
             write(iulog,608)'gradient               (usec/step) = ',gtimer%time_buf(T_gradient)/nmax
             write(iulog,608)'state.grad_lnps update (usec/step) = ',gtimer%time_buf(T_grad_lnps_update)/nmax
             write(iulog,608)'surface pressure       (usec/step) = ',gtimer%time_buf(T_surf_press)/nmax
             write(iulog,608)'compute p and delta p  (usec/step) = ',gtimer%time_buf(T_p_deltap)/nmax
             write(iulog,608)'compute vgrad_lnps     (usec/step) = ',gtimer%time_buf(T_vgrad_lnps)/nmax
             write(iulog,608)'divergence             (usec/step) = ',gtimer%time_buf(T_div)/nmax
             write(iulog,608)'vorticity              (usec/step) = ',gtimer%time_buf(T_vort)/nmax
             write(iulog,608)'hydrostatic            (usec/step) = ',gtimer%time_buf(T_hydrostatic)/nmax
             write(iulog,608)'omega_p eq             (usec/step) = ',gtimer%time_buf(T_omega_p)/nmax
             write(iulog,608)'lnps tendency          (usec/step) = ',gtimer%time_buf(T_lnps_tendency)/nmax
             write(iulog,608)'eta_dot_dp_deta        (usec/step) = ',gtimer%time_buf(T_eta_dot_dp_deta)/nmax
             write(iulog,608)'vertical advection (3) (usec/step) = ',gtimer%time_buf(T_vert_advect)/nmax
             write(iulog,608)'phi + kinetic energy   (usec/step) = ',gtimer%time_buf(T_phik)/nmax
             write(iulog,608)'compute gradp          (usec/step) = ',gtimer%time_buf(T_gradp)/nmax
             write(iulog,608)'interpolate vel. grid  (usec/step) = ',gtimer%time_buf(T_inter_vel)/nmax
          end if

          if(TIMER_DETAIL(2,timer)) then
             write(*,110)gtimer%time_buf(T_pack)/nmax
             write(*,140)gtimer%time_buf(T_unpack)/nmax
             write(*,42)gtimer%time_buf(T_rotate)/nmax
          endif

          ! prim tracer advection time
          write(iulog,44)gtimer%time_buf(T_ADVEC_TRACERS)/nmax
          if(TIMER_DETAIL(2,timer)) then
             write(iulog,604)'RkIntegrator           (usec/step) = ',gtimer%time_buf(T_rkintegrator)/nmax
             
             if(TIMER_DETAIL(3,timer)) then
                write(iulog,608)'filter_v                    (usec/step) = ',gtimer%time_buf(T_filter_v)/nmax
                write(iulog,608)'RkFlux                      (usec/step) = ',gtimer%time_buf(T_rkflux)/nmax

                write(iulog,604)'Prim_Condense()'
                write(iulog,608)'Saturation_Vapor_Pressure() (usec/step) = ',gtimer%time_buf(T_sat_vap_pressure)/nmax
                write(iulog,608)'Vapor_Pressure()            (usec/step) = ',gtimer%time_buf(T_vap_pressure)/nmax
                write(iulog,608)'Mixing_Ratio()              (usec/step) = ',gtimer%time_buf(T_mixing_ratio)/nmax
             end if

          end if

          ! timelevel_update()
          write(iulog,311)gtimer%time_buf(T_timelevel_update)/(nmax-1)
 
          ! movie output time
          write(iulog,313)gtimer%time_buf(T_movie_output)/(nmax-1)

          ! accum output
          write(iulog,314)gtimer%time_buf(T_accum_output)/(nmax-1)

          if(timer%detail > 1) write(iulog,3) iam, avg_barrier_time,max_barrier_time,imbalance_time
          if(total_comm > 0.0)  & 
               write(iulog,6)(total_comm)/nmax,  (dble(nmax)*dble(TotalComVolume)/total_comm)
          ! advance and tracer advect boundary exchange time
          if(time_boundary > 0.0)  &
               write(iulog,131)(time_boundary)/nmax,(dble(nmax)*dble(TotalComVolume)/time_boundary)
          ! prim_filter boundary exchange time
          if(gtimer%time_buf(T_filter_boundary) > 0.0) &
               write(iulog,220)(gtimer%time_buf(T_filter_boundary)/nmax)
          if(timer%time_buf(T_send_time) > 0.0) write(iulog,91)timer%time_buf(T_send_time)/nmax
          if(timer%time_buf(T_recv_time) > 0.0) write(iulog,92)timer%time_buf(T_recv_time)/nmax
          if(timer%time_buf(T_wait_time) > 0.0) write(iulog,93)timer%time_buf(T_wait_time)/nmax
          if(timer%time_buf(T_copy_time) > 0.0) write(iulog,94)timer%time_buf(T_copy_time)/nmax

!          write(iulog,5)hybrid%par%nprocs,hybrid%Nthreads,ne,nlev,time_adv/nmax, &
!               flops_tot/time_adv,1.0E-6*flops_tot/dble(nmax)

          time_tmp = gtimer%time_buf(T_adv)+gtimer%time_buf(T_advec_tracers)

          if(time_tmp > 0.0) then 
		flop_rate = flops_tot/time_tmp
	  else
		flop_rate = 0.0
          endif
          write(iulog,5)hybrid%par%nprocs,hybrid%Nthreads,np,nv,ne,nlev, &
               gtimer%time_buf(T_main_loop)/(nmax-1), flop_rate, &
               1.0E-6*flops_tot/dble(nmax)

604       format (t4,a,f13.3)
608       format (t8,a,f13.3)

       endif
    endif





300 format("INITIALIZATION:        (usec)= ",f11.3)
301 format("   TOPOLOGY            (usec)= ",f11.3)
302 format("   PARTITION           (usec)= ",f11.3)
303 format("   METAGRAPH           (usec)= ",f11.3)
304 format("   SCHEDULE            (usec)= ",f11.3)
305 format("   MASS                (usec)= ",f11.3)
306 format("   TOPO_INIT           (usec)= ",f11.3," (usec/element) =",f11.3)
307 format("   CHECKSUM            (usec)= ",f11.3)
308 format("   INITRESTART         (usec)= ",f11.3)
309 format("MAIN TSTEP LOOP        (usec/step)=",f11.3)
310 format("APPLYCOLUMNMODEL()     (usec/step)=",f11.3)
311 format("TIMELEVELUPDATE()      (usec/step)=",f11.3)
312 format("TOT (ADVN+ADVEC Tracr) (usec/step)=",f11.3," Mflops     =",f16.3)
313 format("MOVIE OUTPUT           (usec/step)=",f11.3)
314 format("ACCUM OUTPUT           (usec/step)=",f11.3)
3   format("IAM: ",i4," LOAD INBALANCE   (usec/step) {AVG,MAX,Ratio}=",f11.3,f11.3,f5.2)
4   format("ADVANCE EXP FLOPS      (usec/step)=",f11.3," Mflops     =",f16.3)
444 format("ADVANCE                (usec/step)=",f11.3)
335 format("    DIVERGENCE_V2V     (usec/step)=",f11.3," MFlops     =",f16.3)
337 format("    GRADIENT_V2V       (usec/step)=",f11.3," MFlops     =",f16.3)
338 format("    VORTICITY          (usec/step)=",f11.3," MFlops     =",f16.3)
 42 format("    ROTATE             (usec/step)=",f11.3)
44  format("ADVEC TRACERS          (usec/step)=",f11.3)
6   format("COMMUNICATION TOTAL    (usec/step)=",f11.3," Mbytes/sec =",f11.3)
91  format("              SEND     (usec/step)=",f11.3)
92  format("              RECV     (usec/step)=",f11.3)
93  format("              WAIT     (usec/step)=",f11.3)
94  format("              COPY     (usec/step)=",f11.3)

5   format("ND=",i5," THR=",i2," NP=",i2," NV=",i2," NE=",i2," NLEV=",i2, &
         " EXP TOT  (usec/step)=",f11.3," Mflops = ",f16.3," Total (Mflop/step) = ",e16.8)
51  format("ND=",i5," THR=",i2," NP=",i2," NV=",i2," NE=",i2," NLEV=",i2, &
         " SEMI-IMPLICIT TOT  (usec/step)=",f11.3," Mflops = ",f16.3," Total (Mflop/step) = ",e16.8)

9   format("MEAN SOLVER ITERATIONS = ",f9.4)
10  format("SI STEP                (usec/step)=",f11.3," Mflops     =",f16.3)
100 format("   PRE PACK COMP    #1 (usec/step)=",f11.3)
110 format("   PACK             #1 (usec/step)=",f11.3," Mbytes/sec =",f11.3)
130 format("   ADV BNDRY EXCH   #1 (usec/step)=",f11.3," Mbytes/sec =",f11.3)
131 format("   ADV BNDRY EXCH      (usec/step)=",f11.3," Mbytes/sec =",f11.3)
140 format("   UNPACK           #1 (usec/step)=",f11.3," Mbytes/sec =",f11.3)
150 format("   POST UNPACK COMP #1 (usec/step)=",f11.3)
160 format("   PRE PACK COMP    #2 (usec/step)=",f11.3)
170 format("   PACK             #2 (usec/step)=",f11.3," Mbytes/sec =",f11.3)
190 format("   ADV BNDRY EXCH   #2 (usec/step)=",f11.3," Mbytes/sec =",f11.3)
200 format("   UNPACK           #2 (usec/step)=",f11.3," Mbytes/sec =",f11.3)
210 format("   POST UNPACK COMP #2 (usec/step)=",f11.3)
220 format("   FILTER BNDRY EXCH   (usec/step)=",f11.3)

20  format("PCG SOLVER             (usec/step)=",f11.3," Mflops     =",f16.3," (Mflop/step) = ",e16.8)
30  format("   CONGRAD CORE        (usec/step)=",f11.3," Mflops     =",f16.3)
35  format("   CONGRAD REDUCE      (usec/step)=",f11.3)
31  format("   REDUCE BARRIER      (usec/step)=",f11.3)
36  format("           SMP         (usec/step)=",f11.3)
37  format("           MPI         (usec/step)=",f11.3)
40  format("   PRECONDITIONER      (usec/step)=",f11.3," Mflops     =",f16.3)
50  format("   HELMHOLTZ OPER      (usec/step)=",f11.3," Mflops     =",f16.3)
55  format("   HELM BNDRY EXCH     (usec/step)=",f11.3," Mbytes/sec =",f16.3)
60  format("ND=",i5," THR=",i2," NE=",i2," NLEV=",i2, &
         " SI  TOT  (usec/step)=",f11.3," Mflops =",f16.3," Total (Mflop/step) = ",e16.8,e16.8)
  end subroutine prim_flops_report

  subroutine FlopCountError(ne,nlev,filter_freq)

    integer, intent(in) :: ne,nlev,filter_freq

#if 0
    write(iulog,10) ne,nlev,filter_freq 
#endif

10  format("prim_flops_report:  Sorry I don't have flop counts for ","NE=",i3,"NLEV=",i3,"FILTER_FREQ=",i3)

  end subroutine FlopCountError

end module flops_mod
