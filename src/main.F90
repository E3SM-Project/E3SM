#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


program main
  use init_mod, only : init
  ! -----------------------------------------------
  use sweq_mod, only : sweq, sweq_rk, define_g, define_kmass
  ! -----------------------------------------------
  !      use kinds
  ! -----------------------------------------------
  use parallel_mod, only : parallel_t, initmp, syncmp, haltmp
  ! -----------------------------------------------
  use thread_mod, only : nthreads, omp_get_thread_num, omp_set_num_threads
  ! -----------------------------------------------
  !      use time_mod
  ! -----------------------------------------------
  use dimensions_mod, only : nelemd
  ! -----------------------------------------------
  use domain_mod, only : domain1d_t, decompose
  ! -----------------------------------------------
  use element_mod, only : element_t
  ! -----------------------------------------------
  use fvm_control_volume_mod, only : fvm_struct
  ! -----------------------------------------------
  !      use state_mod
  ! -----------------------------------------------
  use edge_mod, only : EdgeBuffer_t
  ! -----------------------------------------------
  use reduction_mod, only : ReductionBuffer_ordered_1d_t
  ! -----------------------------------------------
  use perf_mod, only : t_initf, t_prf, t_finalizef, t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------
  use control_mod, only : integration

  implicit none
  type (element_t), pointer :: elem(:)
  type (fvm_struct), pointer  :: fvm(:)
  
  type (EdgeBuffer_t)  :: edge1            ! 1 component edge buffer (1, 3d scalar field)
  type (EdgeBuffer_t)  :: edge2            ! 2 component edge buffer (1, 3d vector field)
  type (EdgeBuffer_t)  :: edge3            ! 3 component edge buffer (1, 3d vector + 1 3d scalar field)
  type (ReductionBuffer_ordered_1d_t)  :: red    ! reduction buffer for cg
  type (parallel_t)    :: par              ! parallel structure for distributed memory programming
  type (domain1d_t), allocatable :: dom_mt(:)

  integer nets,nete
  integer ithr
  integer ierr

  ! =====================================================
  ! Begin executable code set distributed memory world...
  ! =====================================================

  par=initmp()
  call t_initf('input.nl',LogPrint=par%masterproc, &
       Mpicom=par%comm, MasterTask=par%masterproc)
  call t_startf('Total')
  
#ifdef _FVM
  if(par%masterproc) then
    print *, '----------------------------'
    print *, 'RUN SWEQX with FVM TRACERS  '
    print *, '----------------------------'
  endif
#endif

#ifdef _FVM
    call init(elem,edge1,edge2,edge3,red,par,fvm)
#else    
    call init(elem,edge1,edge2,edge3,red,par)
#endif    
  ! =====================================================
  ! Allocate state variables
  ! =====================================================

  if(par%masterproc) print *,"allocating state variables..."
  !JMD allocate(state(nelemd))

  ! =====================================
  ! Set number of threads...
  ! =====================================

  if(par%masterproc) print *,"Main:NThreads=",NThreads

  call omp_set_num_threads(NThreads)

  allocate(dom_mt(0:NThreads-1))
  do ithr=0,NThreads-1
     dom_mt(ithr)=decompose(1,nelemd,NThreads,ithr)
  end do

  ! =====================================
  !  Sync-up to make sure timing is clean
  ! =====================================

  call syncmp(par)

  ! =====================================
  ! Begin threaded region...
  ! =====================================
#if (! defined ELEMENT_OPENMP)
  !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ithr,nets,nete)
#endif
  ithr=omp_get_thread_num()
  nets=dom_mt(ithr)%start
  nete=dom_mt(ithr)%end

  !
  ! ================================================
  ! Initialize thread decomposition
  ! ================================================
  !
  write(6,9) par%rank,ithr,nets,nete 
9 format("process: ",i2,1x,"thread: ",i2,1x,"element limits: ",i4," - ",i4)

  if(integration == "runge_kutta")then
!define_g is setting g=1 for swirl test case, for all other tcases g=g
!define_kmass sets kmass to the level with air density, if provided. right now it is provided 
!only for swirl test case. also, kmass is only for lim8 to work, and lim8 is only for RK
     call define_g
     call define_kmass
     call sweq_rk(elem,edge1,edge2,edge3,red,par,ithr,nets,nete)
  else
     call sweq(elem,fvm,edge1,edge2,edge3,red,par,ithr,nets,nete)
  endif

#if (! defined ELEMENT_OPENMP)
  !$OMP END PARALLEL
#endif
  ! ================================================
  ! End distributed memory region
  ! ================================================
  call t_stopf('Total')
  call t_prf('HommeSWTime',par%comm)
  call t_finalizef()
  call haltmp("exiting program...")
  deallocate(elem)
end program main
