#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

program prim_main

  ! Main routine for stand-alone HOMME atmospheric dyanmcial-core

  use hybvcoord_mod,  only: hvcoord_t, hvcoord_init
  use parallel_mod,   only: parallel_t, initmp, syncmp, haltmp, abortmp
  use hybrid_mod,     only: hybrid_t
  use time_mod,       only: tstep, nendstep, timelevel_t
  use dimensions_mod, only: nelemd
  use domain_mod,     only: domain1d_t
  use element_mod,    only: element_t
  use spelt_mod,      only: spelt_struct
  use common_io_mod,  only: output_dir
  use perf_mod,       only: t_initf,t_prf,t_finalizef,t_startf,t_stopf
  use restart_io_mod, only: writerestart
  use hybrid_mod,     only: hybrid_create
  use common_movie_mod, only: nextoutputstep
  use fvm_control_volume_mod, only: fvm_struct

  use control_mod, only: restartfreq, vfile_mid, vfile_int, &
  runtype, integration, statefreq, tstep_type

  use prim_driver_mod, only: prim_init1, prim_init2, prim_run, &
  prim_finalize, leapfrog_bootstrap, prim_run_subcycle

  use thread_mod, only: nthreads, vert_num_threads, omp_get_thread_num,&
  omp_set_num_threads, omp_get_nested

#ifdef _REFSOLN
  use prim_state_mod, only : prim_printstate_par
#endif

  ! Select data IO method (parallel-IO interpolate, or prim native)
#ifdef PIO_INTERP
  use interp_movie_mod, only: interp_movie_output, interp_movie_finish, interp_movie_init
  use interpolate_driver_mod, only : interpolate_driver
#else
  use prim_movie_mod, only : prim_movie_output, prim_movie_finish, prim_movie_init
#endif

  implicit none

  type (element_t),  pointer  :: elem(:)        ! pointer to element array
  type (domain1d_t), pointer  :: dom_mt(:)      ! element start,end for each thread
  type (parallel_t)           :: par            ! mpi data structure
  type (hybrid_t)             :: hybrid         ! mpi / omp data-structure
  type (TimeLevel_t)          :: tl             ! time level index data struct
  type (hvcoord_t)            :: hvcoord        ! hybrid vertical coordinates

  integer nets,nete                             ! start, end element indices
  integer ithr                                  ! thread index
  integer ierr                                  ! error state
  integer nstep                                 ! main loop timestep counter

  ! Select data structure for tracers
#if defined(_SPELT)
  type (spelt_struct), pointer :: fvm(:)
#else
  type (fvm_struct),   pointer :: fvm(:)
#endif

  ! Enable HTRACE if needed
#ifdef _HTRACE
  integer htype,pflag
#endif

  ! Initialize MPI
  par=initmp()

  ! Initialize timers
  call t_initf('input.nl',logprint=par%masterproc, mpicom=par%comm, mastertask=par%masterproc)
  call t_startf('Total')

  ! PRIM_INIT1: Initialize mesh and restart files. Partition domain
  call t_startf('prim_init1')
  call prim_init1(elem,  fvm, par,dom_mt,tl)
  call t_stopf('prim_init1')

  ! Initialize HTRACE if needed
#ifdef _HTRACE
  htype = 19
  pflag = 0
  call TRACE_INIT(htype,pflag)
#endif

  ! Enable thread-nesting for omp in both Horizontal and Vertical
#if (defined HORIZ_OPENMP && defined COLUMN_OPENMP)
  call omp_set_nested(.true.)
  if (omp_get_nested() == 0) call haltmp("Nested threading required but not available. Set OMP_NESTED=true")
#endif

  ! Set number of vertical threads for Horizontal omp
#if (defined HORIZ_OPENMP)
  !$OMP PARALLEL NUM_THREADS(nthreads), DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid)
  call omp_set_num_threads(vert_num_threads)
#endif

  ! Display element partition
  ithr = omp_get_thread_num()
  nets = dom_mt(ithr)%start
  nete = dom_mt(ithr)%end
#if (! defined ELEMENT_OPENMP)
  !$OMP CRITICAL
#endif
  if (par%rank<100) write(6,9) par%rank,ithr,nets,nete
9 format("process: ",i2,1x,"thread: ",i2,1x,"element limits: ",i5," - ",i5)
#if (defined HORIZ_OPENMP)
  !$OMP END CRITICAL
  !$OMP END PARALLEL
#endif

  ! Initialize vertical coordinates (CAM initializes hvcoord externally)
  ithr    = omp_get_thread_num()
  hybrid  = hybrid_create(par,ithr,1)
  nets    = 1
  nete    = nelemd
  hvcoord = hvcoord_init(vfile_mid, vfile_int, .true., hybrid%masterthread, ierr)
  if (ierr /= 0) call haltmp("error in hvcoord_init")

  ! If runtype<0, interpolate netcdf file and exit
#ifdef PIO_INTERP
  if(runtype<0) then
     call interpolate_driver(elem, hybrid)
     call haltmp('interpolation complete')
  end if
#endif

  ! PRIM_INIT2: Set initial conditions. Initialize filters. Compute max CFL.
  if(par%masterproc) print *,"Primitive Equation Initialization..."
  ithr   = omp_get_thread_num()
  hybrid = hybrid_create(par,ithr,nthreads)
  nets   = dom_mt(ithr)%start
  nete   = dom_mt(ithr)%end
  call t_startf('prim_init2')
  call prim_init2(elem, fvm,  hybrid,nets,nete,tl, hvcoord)
  call t_stopf('prim_init2')

#if (defined HORIZ_OPENMP)
  !$OMP END PARALLEL
#endif

  ! Ensure that the output directory for history files exists
  ithr   = omp_get_thread_num()
  hybrid = hybrid_create(par,ithr,1)
  if (hybrid%masterthread) then
     open(unit=447,file=trim(output_dir) // "/output_dir_test",iostat=ierr)
     if ( ierr==0 ) then
        print *,'Directory ',output_dir, ' does exist: initialing IO'
        close(447)
     else
        print *,'Error creating file in directory ',output_dir
        call haltmp("Please be sure the directory exist or specify 'output_dir' in the namelist.")
     end if
  endif

  ! Initialize history files
#ifdef PIO_INTERP
   call interp_movie_init( elem, hybrid, 1, nelemd, hvcoord, tl )
#else
  call prim_movie_init( elem, par, hvcoord, tl )
#endif

  ! Output initial state for NEW runs (not restarts or branch runs)
  if (runtype == 0 ) then
#ifdef PIO_INTERP
     call interp_movie_output(elem, tl, hybrid, 0d0, 1, nelemd,fvm=fvm, hvcoord=hvcoord)
#else
     call prim_movie_output(elem, tl, hvcoord, hybrid, 1,nelemd, fvm)
#endif
  endif

  ! Use leapfrog-boostrap to start semi-implicit runs
  if(integration == 'semi_imp') then
     if (runtype /= 1 ) then
        if(hybrid%masterthread) print *,"Leapfrog bootstrap initialization..."
        call leapfrog_bootstrap(elem, hybrid,1,nelemd,tstep,tl,hvcoord)
     endif
  endif

  if(par%masterproc) print *,"Entering main timestepping loop"

  ! Begin main timestepping loop
  do while(tl%nstep < nEndStep)

#if (defined HORIZ_OPENMP)
     !$OMP PARALLEL NUM_THREADS(nthreads), DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid)
     call omp_set_num_threads(vert_num_threads)
#endif
     ithr   = omp_get_thread_num()
     hybrid = hybrid_create(par,ithr,nthreads)
     nets   = dom_mt(ithr)%start
     nete   = dom_mt(ithr)%end

     ! Loop until next output step
     nstep = nextoutputstep(tl)
     do while(tl%nstep<nstep)

        ! Advanced solvers
        call t_startf('prim_run')
        if (tstep_type>0) then  ! forward in time subcycled methods
           call prim_run_subcycle(elem, fvm, hybrid,nets,nete, tstep, tl, hvcoord,1)
        else                    ! leapfrog
           call prim_run(elem, hybrid,nets,nete, tstep, tl, hvcoord, "leapfrog")
        endif
        call t_stopf('prim_run')

     end do

#if (defined HORIZ_OPENMP)
     !$OMP END PARALLEL
#endif

     ! Output data to history file
     ithr   = omp_get_thread_num()
     hybrid = hybrid_create(par,ithr,1)
#ifdef PIO_INTERP
     call interp_movie_output(elem, tl, hybrid, 0d0, 1, nelemd,fvm=fvm, hvcoord=hvcoord)
#else
     call prim_movie_output(elem, tl, hvcoord, hybrid, 1,nelemd, fvm)
#endif

    ! Print state to screen
#ifdef _REFSOLN
     call prim_printstate_par(elem, tl,hybrid,hvcoord,nets,nete, par)
#endif 

     ! Write restart files if required 
     if((restartfreq > 0) .and. (MODULO(tl%nstep,restartfreq) ==0)) then
        call WriteRestart(elem, ithr,1,nelemd,tl)
     endif

  end do

  ! Exit main time-stepping loop
  if(par%masterproc) print *,"Finished main timestepping loop",tl%nstep
  call prim_finalize(hybrid)
  if(par%masterproc) print *,"closing history files"
#ifdef PIO_INTERP
  call interp_movie_finish
#else
  call prim_movie_finish
#endif

  ! Stop all timers and write timing data
  call t_stopf('Total')
  if(par%masterproc) print *,"writing timing data"
  call t_prf('HommeTime', par%comm)
  if(par%masterproc) print *,"calling t_finalizef"
  call t_finalizef()

  ! Finalize MPI and exit
  call haltmp("exiting program...")

end program prim_main








