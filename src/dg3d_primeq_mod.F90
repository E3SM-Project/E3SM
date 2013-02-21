#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


!=======================================================================================================! 
module dg3d_primeq_mod
!=======================================================================================================!
 implicit none
 integer, public:: tasknum
!=======================================================================================================! 
 contains
!=======================================================================================================!
subroutine primeq_dg(elem,edge1,edge2,edge3,red,par,ithr,nets,nete, my_hvcoord)
    !-----------------
    use dg3d_vertical_mod, only: eta_levels
    !-----------------
    use kinds, only : real_kind
    !-----------------
    use element_mod, only : element_t
    !-----------------
    use physical_constants, only : g
    !-----------------
    use parallel_mod, only : parallel_t, syncmp
    !-----------------
    use thread_mod, only : nthreads
    !-----------------
    use hybrid_mod, only : hybrid_t, hybrid_create
    !-----------------
    use hybvcoord_mod, only : hvcoord_t
    !-----------------
    use time_mod, only : timelevel_t , tstep, secpday, time_at, nmax, timelevel_update, timelevel_init
    !-----------------
    use derivative_mod, only : derivative_t, derivative_stag_t, derivinit, deriv_print
    !-----------------
    use dimensions_mod, only : np, nlev, npsq
    !-----------------
    use interp_movie_mod, only: interp_movie_init, interp_movie_output, interp_movie_finish
    !-----------------
    use global_norms_mod, only : test_global_integral
    !-----------------
    use quadrature_mod, only : quadrature_t, gauss,gausslobatto
    !-----------------
    use edge_mod, only : EdgeBuffer_t
    ! ----------------
    use reduction_mod, only : ReductionBuffer_ordered_1d_t
    !-----------------
    use state_mod, only : printstate_dg
    !-----------------
    use filter_mod, only : filter_t, taylor_filter_create, fm_filter_create, fm_transfer, bv_transfer
    !-----------------
    use solver_mod, only : blkjac_t, blkjac_init
    !-----------------
    use cg_mod, only : cg_t, cg_create
    !-----------------
    use restart_io_mod, only : readrestart, writerestart
    !-----------------
    use control_mod, only : integration, filter_mu, filter_type, transfer_type, debug_level, test_case, &
         restartfreq, statefreq, runtype, s_bv, p_bv, wght_fm, kcut_fm, accumfreq, topology
    !----------------- 
    use dg3d_prim_mod, only :  dg3d_advance
    !----------------- 
    use dg3d_errors_mod, only : jw_bcl_errors, cube2polar, jw_bcl_zonal   
    !----------------- 
!=======================================================================================================!   
    implicit none
#ifdef _HPM
#include "f_hpm.h"
#endif
!=======================================================================================================! 
    integer, parameter :: facs = 4            ! starting face number to print
    integer, parameter :: face = 4            ! ending  face number to print
    type (element_t), intent(inout)                :: elem(:)
    type (EdgeBuffer_t), intent(in)             :: edge1 ! edge buffer entity             (shared)
    type (EdgeBuffer_t), intent(in)             :: edge2 ! edge buffer entity             (shared)
    type (EdgeBuffer_t), intent(in)             :: edge3 ! edge buffer entity             (shared)
    type (ReductionBuffer_ordered_1d_t),intent(in):: red   ! reduction buffer               (shared)
    type (parallel_t), intent(in)               :: par   ! distributed parallel structure (shared)
    integer, intent(in)                         :: ithr  ! thread number                  (private)
    integer, intent(in)                         :: nets  ! starting thread element number (private)
    integer, intent(in)                         :: nete  ! ending thread element number   (private)
!=======================================================================================================!  
!	Local thread (private) memory
!=======================================================================================================!  
    type (hvcoord_t)            :: hybvc, my_hvcoord           ! hybrid vertical coordinate
    real (kind=real_kind)       :: dt              ! 'timestep dependent' timestep
    real (kind=real_kind)       :: pmean           ! mean geopotential
    type (derivative_t)         :: deriv           ! derivative struct
    type (TimeLevel_t)          :: tl              ! time level struct
    type (blkjac_t),allocatable :: blkjac(:)  
    type (cg_t)                 :: cg              ! conjugate gradient struct
    real (kind=real_kind)       :: lambdasq(nlev)  ! Helmholtz length scale
    type (hybrid_t)             :: hybrid
    type (quadrature_t)         :: gll,gs          ! gauss-lobatto and gauss wts and pts


    real (kind=real_kind) :: Tp(np)          ! transfer function (pressure and velocity grid)
    type (filter_t)       :: flt           ! Filter structure for both v and p grid
    type (quadrature_t)   :: gp           ! quadratures on velocity and pressure grids
    real (kind=real_kind) :: solver_wts(npsq,nete-nets+1)        ! solver weights array for nonstaggered grid

    integer  :: simday

    integer :: point
    integer :: i,j,iptr
    integer :: it,ie,k
    integer :: ntmp
    integer :: nm1,n0,np1
    integer :: nstep
    integer :: stage
    real*8  :: tot_iter
!=======================================================================================================!
#ifdef _HPM
    type (parallel_t):: taskid
    integer :: tasktmp= 0
#endif
!=======================================================================================================!    
    logical, parameter :: Debug      = .FALSE.    
    logical, parameter :: CPU_Report = .FALSE.
!=======================================================================================================!
    ! initialize vertical coordinate
    call eta_levels(hybvc%hyai, hybvc%hybi, hybvc%etam, hybvc%etai)

    if(Debug) print *,'seam: point #1'
    ! ==========================
    ! begin executable code
    ! ==========================
    hybrid = hybrid_create(par,ithr,NThreads)

    if (topology == 'cube') then
       call test_global_integral(elem,hybrid,nets,nete)  
    end if

    if(Debug) print *,'seam: point #2'
    ! ==================================
    ! Initialize derivative structure
    ! ==================================

    call derivinit(deriv)

    if (hybrid%par%masterproc .AND. ithr == 0) then
       call deriv_print(deriv)
    end if

    ! ========================================
    ! Initialize velocity and pressure grid
    ! quadrature points...
    ! ========================================

    gp =gausslobatto(np)

    if(Debug) print *,'seam: point #3'
    ! ==========================================
    ! Initialize pressure and velocity grid 
    ! filter matrix...
    ! ==========================================

    if (hybrid%par%masterproc .AND. ithr==0) then
       print *,'transfer function type in seam=',transfer_type
       print *,'filter type            in seam=',filter_type
    end if

    if (transfer_type == 'bv') then
       Tp= bv_transfer(p_bv,s_bv,np)
    else if (transfer_type == 'fm') then
       Tp= fm_transfer(kcut_fm,wght_fm,np)
    end if

    if (filter_type == 'taylor') then
       flt= taylor_filter_create(Tp, filter_mu, gp)
    else if (filter_type == 'fischer') then
       flt= fm_filter_create(Tp, filter_mu, gp)
    end if
    if(Debug) print *,'seam: point #4'

    if (hybrid%par%masterproc .AND. ithr == 0) then
#if 0
       print *,'Filter:'
       do j=1,np
          do k=1,np
             print *,'F(',k,',',j,')=',flt%FmatV(k,j)
          end do
       end do

       if (transfer_type=='bv') then
          call bvsigma_test(p_bv)
       end if
#endif
    end if

    do ie=nets,nete
       iptr=1
       do j=1,np
          do i=1,np
             solver_wts(iptr,ie-nets+1) = elem(ie)%solver_wts(i,j)
             iptr=iptr+1
          end do
       end do
    end do
!=======================================================================================================!
!	Sync-up to make sure timing is clean								!
!=======================================================================================================!    
 if (ithr==0) then
   call syncmp(par)
 endif
#ifdef _HPM
 taskid = hybrid%par
 tasktmp= ithr
if (tasknum==-1) then
 call f_hpminit(taskid,'dg3d')
else
 call f_hpminit(tasknum,'dg3d') 
endif
#endif 
!=======================================================================================================!
!	Begin threaded region...									!
!=======================================================================================================!
!=======================================================================================================!
!=======================================================================================================!
    call TimeLevel_init(tl)
!=======================================================================================================!
    if(Debug) print *,'seam: point #5'
    ! =================================================================
    ! Initialize geopotential and velocity for different test cases...
    ! =================================================================

    if (topology == 'cube') then
       if (runtype .eq. 1) then  	! runtype = 1

          if (hybrid%par%masterproc.and.ithr==0) then
             print *,'runtype: RESTART of DG Shallow Water equations'
          end if
          !============================
          ! Read in the restarted state 
          !============================
          call ReadRestart(elem,ithr,nete,nets,tl)
          !================================================
          ! Print out the state variables 
          !================================================          
	  if(hybrid%par%masterproc .AND. ithr == 0) then
             print *,tl%nstep,'time=',Time_at(tl%nstep)/secpday,' days'
          endif
          call printstate_dg(elem,pmean,g,tl%n0,hybrid,nets,nete)
 
       elseif (runtype .eq. 0) then 	! runtype = 0

          if (hybrid%par%masterproc.and.ithr==0) then
             print *,'runtype: INITIAL of DG Primitive equations'
          end if
          ! ==============================================
          ! Output initial picture of geopotential...
          ! ============================================== 
if (statefreq>0) then	
          call interp_movie_init(elem, hybrid, nets, nete, my_hvcoord, tl)
endif
          if(Debug) print *,'seam: point #6' 
	  
          ! =================================
          ! Call advance: Initialization
          ! =================================
          if (integration == 'explicit') then
             if(hybrid%par%masterproc .AND. ithr == 0) then 
                print *,'calling initial DG 3D advance...'
	     endif
             if (test_case(1:6) == 'jw_bcl' .or. test_case(1:10) == 'heldsuarez') then
                call dg3d_advance(elem,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
             endif
          endif
          if(Debug) print *,'seam: point #7'

          ! ===============================================================
          ! Print Min/Max/Sum of State vars before first timestep
          ! ===============================================================
       if (tl%nstep==0) then 
          if(hybrid%par%masterproc .AND. ithr == 0) then
             print *,tl%nstep,'time=',Time_at(tl%nstep)/secpday,' days'
          endif
          call printstate_dg(elem,pmean,g,tl%n0,hybrid,nets,nete)  
       endif
if (statefreq>0) then
          call interp_movie_output(elem, tl, hybrid, pmean, nets, nete, hvcoord=my_hvcoord )
          !call jw_bcl_errors(elem, 7, tl, hybrid, nets, nete)
	  !call jw_bcl_zonal(elem, tl, hybrid)    

endif

	  ! =================================
          ! update time level pointers ...
          ! =================================
          call TimeLevel_update(tl,'forward')
          if(Debug) print *,'seam: point #8'

       endif  ! if initial run 

    end if ! if topology == 'cube'
!=======================================================================================================!
!=======================================================================================================!
!	Sync-up to make sure timing is clean								!
!=======================================================================================================!    
 if (ithr==0) then
   call syncmp(par)
 endif
 dt = tstep
 tot_iter=0.0    
!=======================================================================================================!
!	Main timestepping loop										!
!=======================================================================================================!    
 do while(tl%nstep<=nmax)

       ! =================================
       ! Call advance
       ! =================================
!=======================================================================================================!
!	MPI TRACE START											!
!=======================================================================================================!
#ifdef _BGL
 call trace_start()
#endif        
#ifdef _HTRACE
 point= 1
 call EVENT_POINT(point)
#endif
#ifdef _HPM
 call f_hpmstart(5,'dg3d advance')
#endif
!=======================================================================================================!
       if(Debug) print *,'seam: point #9'

       if(integration == 'explicit') then
	if (test_case(1:6)=='jw_bcl' .or. test_case(1:10) == 'heldsuarez') then
           call dg3d_advance(elem,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
	endif
       endif	
       if(Debug) print *,'seam: point #10'
!=======================================================================================================!
!	MPI TRACE STOP											!
!=======================================================================================================!
#ifdef _BGL
 call trace_stop()
#endif
#ifdef _HTRACE
 point = 2
 call EVENT_POINT(point)
#endif
#ifdef _HPM
 call f_hpmstop(5)
#endif
!=======================================================================================================!  
       ! =================================
       ! update time level pointers
       ! =================================
!       call TimeLevel_update(tl,'leapfrog')
       
       ! ============================================================
       ! Instrumentation alley:
       !
       ! Shallow Water Test Case output files
       ! ============================================================
if (statefreq>0) then
       if (MODULO(tl%nstep,statefreq)==0 .or. MODULO(tl%nstep,nmax)==0 .or. tl%nstep==1) then 
          if(hybrid%par%masterproc .AND. ithr == 0) then
             print *,tl%nstep,'time=',Time_at(tl%nstep)/secpday,' days'
          endif
          call printstate_dg(elem,pmean,g,tl%n0,hybrid,nets,nete)          
          !call jw_bcl_errors(elem, 7, tl, hybrid, nets, nete) 	  
	  !call jw_bcl_zonal(elem, tl, hybrid)         
       end if        
       call interp_movie_output(elem, tl, hybrid, pmean, nets, nete, hvcoord=my_hvcoord )
else
       if (MODULO(tl%nstep,nmax)==0 .or. tl%nstep==1) then 
          if(hybrid%par%masterproc .AND. ithr == 0) then
             print *,tl%nstep,'time=',Time_at(tl%nstep)/secpday,' days'
          endif
          call printstate_dg(elem,pmean,g,tl%n0,hybrid,nets,nete)        
       end if   
endif        
       ! ============================================================
       ! Write restart files if required
       ! ============================================================

       if((restartfreq > 0) .and. (MODULO(tl%nstep,restartfreq) ==0)) then
          call WriteRestart(elem,ithr,nets,nete,tl)
       endif      
       ! =================================
       ! update time level pointers
       ! =================================
       call TimeLevel_update(tl,'forward')

 enddo
!=======================================================================================================!  
    if (ithr==0) then
       call syncmp(par)
    end if    
    if (hybrid%par%masterproc.and.ithr==0) then
	print *,'Finishing DG 3D Time Integration...'
    endif
!=======================================================================================================!
10  format(1x,'Time solver total time (sec)=',f10.2)
!=======================================================================================================!
#ifdef _HPM
 taskid = hybrid%par
 tasktmp= ithr
if (tasknum == -1) then
 call f_hpm_terminate(taskid)
else 
 call f_hpm_terminate(tasknum) 
endif    
#endif  
!=======================================================================================================!
    ! ======================================================
    ! compute and report times...
    ! ======================================================
    ! ==========================
    ! end of the hybrid program
    ! ==========================
if (statefreq>0) then    
   call interp_movie_finish
endif

end subroutine primeq_dg
!=======================================================================================================!
end module dg3d_primeq_mod

