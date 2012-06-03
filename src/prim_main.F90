#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

program prim_main
  ! -----------------------------------------------
#ifdef _PRIM
  use prim_driver_mod, only : prim_init1, prim_init2, prim_run, prim_finalize,&
                              leapfrog_bootstrap, prim_run_subcycle
  use hybvcoord_mod, only : hvcoord_t, hvcoord_init
#else
! use these for dg3d

#endif

  ! -----------------------------------------------
  use parallel_mod, only : parallel_t, initmp, syncmp, haltmp, abortmp
  ! -----------------------------------------------
  use hybrid_mod, only : hybrid_t
  ! -----------------------------------------------
  use thread_mod, only : nthreads, omp_get_thread_num, omp_set_num_threads
  ! ----------------------------------------------- 
  use time_mod, only : tstep, nendstep, timelevel_t, TimeLevel_init
  ! -----------------------------------------------
  use dimensions_mod, only : nelemd, qsize, ntrac
  ! -----------------------------------------------
  use control_mod, only : restartfreq, vfile_mid, vfile_int, runtype, integration, statefreq, tstep_type
  ! -----------------------------------------------
  use domain_mod, only : domain1d_t, decompose
  ! -----------------------------------------------
  use element_mod, only : element_t
  !-----------------------------------------------
  use physics_mod, only : elem_physics_t
  !-----------------------------------------------
  use cslam_control_volume_mod, only : cslam_struct
  ! -----------------------------------------------
  use common_io_mod, only:  output_dir
  ! -----------------------------------------------

#ifdef PIO_INTERP
  use interp_movie_mod, only : interp_movie_output, interp_movie_finish, interp_movie_init
  use interpolate_driver_mod, only : interpolate_driver
#else
  use prim_movie_mod, only : prim_movie_output, prim_movie_finish,prim_movie_init
#endif
  use common_movie_mod, only : nextoutputstep
  use perf_mod, only : t_initf, t_prf, t_finalizef, t_startf, t_stopf ! _EXTERNAL

  !-----------------
  use restart_io_mod , only : restartheader_t, writerestart
  !-----------------
  use hybrid_mod, only : hybrid_create


	
  implicit none
  type (element_t), pointer  :: elem(:)
  type (elem_physics_t), pointer  :: elem_physics(:)
  type (cslam_struct), pointer  :: cslam(:)

  type (hybrid_t)       :: hybrid ! parallel structure for shared memory/distributed memory
  type (parallel_t)                    :: par  ! parallel structure for distributed memory programming
  type (domain1d_t), pointer :: dom_mt(:)
  type (RestartHeader_t)                  :: RestartHeader
  type (TimeLevel_t)    :: tl             ! Main time level struct
  type (hvcoord_t)     :: hvcoord         ! hybrid vertical coordinate struct

  real*8 timeit, et, st
  integer nets,nete
  integer ithr
  integer ierr
  integer nstep
  
  character (len=20)                          :: numproc_char
  character (len=20)                          :: numtrac_char
  
  logical :: dir_e ! boolean existence of directory where output netcdf goes
  
#ifdef _HTRACE
  integer htype,pflag
#endif
  ! =====================================================
  ! Begin executable code set distributed memory world...
  ! =====================================================
  par=initmp()



  ! =====================================
  ! Set number of threads...
  ! =====================================
  call t_initf('input.nl',LogPrint=par%masterproc, &
	Mpicom=par%comm, MasterTask=par%masterproc)
  call t_startf('Total')
  call t_startf('prim_init1')
  call prim_init1(elem, elem_physics,  cslam, par,dom_mt,tl)
  call t_stopf('prim_init1')



#ifdef _HTRACE
  htype = 19
  pflag = 0
  call TRACE_INIT(htype,pflag)
#endif

  ! =====================================
  ! Begin threaded region so each thread can print info
  ! =====================================
#if (! defined ELEMENT_OPENMP)
  !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid)
#endif
  ithr=omp_get_thread_num()
  nets=dom_mt(ithr)%start
  nete=dom_mt(ithr)%end

  ! ================================================
  ! Initialize thread decomposition
  ! ================================================
  if (par%rank<100) then 
     write(6,9) par%rank,ithr,nets,nete 
  endif
9 format("process: ",i2,1x,"thread: ",i2,1x,"element limits: ",i5," - ",i5)
#if (! defined ELEMENT_OPENMP)
  !$OMP END PARALLEL
#endif

  
! back to single threaded
  ithr=omp_get_thread_num()
  hybrid = hybrid_create(par,ithr,1)
  nets=1
  nete=nelemd


  ! ==================================
  ! Initialize the vertical coordinate  (cam initializes hvcoord externally)
  ! ==================================
  hvcoord = hvcoord_init(vfile_mid, vfile_int, .true., hybrid%masterthread, ierr)
  if (ierr /= 0) then
     call haltmp("error in hvcoord_init")
  end if



#ifdef PIO_INTERP
  if(runtype<0) then
     ! Interpolate a netcdf file from one grid to another
     call interpolate_driver(elem, hybrid)
     call haltmp('interpolation complete')
  end if
#endif

  if(par%masterproc) print *,"Primitive Equation Initialization..."
#if (! defined ELEMENT_OPENMP)
  !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid)
#endif
  ithr=omp_get_thread_num()
  hybrid = hybrid_create(par,ithr,nthreads)
  nets=dom_mt(ithr)%start
  nete=dom_mt(ithr)%end

  call t_startf('prim_init2')
  call prim_init2(elem, elem_physics, cslam,  hybrid,nets,nete,tl, hvcoord)
  call t_stopf('prim_init2')
#if (! defined ELEMENT_OPENMP)
  !$OMP END PARALLEL
#endif

  ithr=omp_get_thread_num()
  hybrid = hybrid_create(par,ithr,1)
  
  ! Here we get sure the directory specified
  ! in the input namelist file in the 
  ! variable 'output_dir' does exist.
  ! this avoids a abort deep within the PIO 
  ! library (SIGABRT:signal 6) which in most
  ! architectures produces a core dump.
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
#if 0
  this ALWAYS fails on lustre filesystems.  replaced with the check above
  inquire( file=output_dir, exist=dir_e )
  if ( dir_e ) then
     if(hybrid%masterthread) print *,'Directory ',output_dir, ' does exist: initialing IO'
  else
     if(hybrid%masterthread) print *,'Directory ',output_dir, ' does not exist: stopping'
     call haltmp("Please get sure the directory exist or specify one via output_dir in the namelist file.")
  end if
#endif
  

#ifdef PIO_INTERP
  ! initialize history files.  filename constructed with restart time
  ! so we have to do this after ReadRestart in prim_init2 above
  call interp_movie_init( elem, hybrid, 1, nelemd, hvcoord, tl )
#else
  call prim_movie_init( elem, par, hvcoord, tl )
#endif


  ! output initial state for NEW runs (not restarts or branch runs)
  if (runtype == 0 ) then
#ifdef PIO_INTERP
     call interp_movie_output(elem, tl, hvcoord, hybrid, 1, nelemd)
#else
     call prim_movie_output(elem, tl, hvcoord, hybrid, 1,nelemd, cslam)
#endif
  endif


  ! advance_si not yet upgraded to be self-starting.  use leapfrog bootstrap procedure:
  if(integration == 'semi_imp') then
     if (runtype /= 1 ) then
        if(hybrid%masterthread) print *,"Leapfrog bootstrap initialization..."
        call leapfrog_bootstrap(elem, elem_physics, hybrid,1,nelemd,tstep,tl,hvcoord)
     endif
  endif
  

  if(par%masterproc) print *,"Entering main timestepping loop"
  do while(tl%nstep < nEndStep)
#if (! defined ELEMENT_OPENMP)
     !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid)
#endif
     ithr=omp_get_thread_num()
     hybrid = hybrid_create(par,ithr,nthreads)
     nets=dom_mt(ithr)%start
     nete=dom_mt(ithr)%end
     
     nstep = nextoutputstep(tl)
     do while(tl%nstep<nstep)
        call t_startf('prim_run')
        if (tstep_type==1) then  ! foreward in time methods
           call prim_run_subcycle(elem, cslam, hybrid,nets,nete, tstep, tl, hvcoord)
        else  ! leapfrog
           call prim_run(elem, elem_physics, hybrid,nets,nete, tstep, tl, hvcoord, "leapfrog")
        endif
        call t_stopf('prim_run')
     end do
#if (! defined ELEMENT_OPENMP)
     !$OMP END PARALLEL
#endif

     ithr=omp_get_thread_num()
     hybrid = hybrid_create(par,ithr,1)
#ifdef PIO_INTERP
     call interp_movie_output(elem, tl,hvcoord, hybrid, 1,nelemd)
#else
     call prim_movie_output(elem, tl, hvcoord, hybrid, 1,nelemd, cslam)
#endif
     ! ============================================================
     ! Write restart files if required 
     ! ============================================================
     if((restartfreq > 0) .and. (MODULO(tl%nstep,restartfreq) ==0)) then 
        call WriteRestart(elem, ithr,1,nelemd,tl)
     endif
  end do
  if(par%masterproc) print *,"Finished main timestepping loop",tl%nstep
  call prim_finalize(hybrid)
  if(par%masterproc) print *,"closing history files"
#ifdef PIO_INTERP
  call interp_movie_finish
#else
  call prim_movie_finish
#endif


  call t_stopf('Total')
  if(par%masterproc) print *,"writing timing data"
!   write(numproc_char,*) par%nprocs
!   write(numtrac_char,*) ntrac
!   call system('mkdir -p '//'time/'//trim(adjustl(numproc_char))//'-'//trim(adjustl(numtrac_char))) 
!   call t_prf('time/HommeCSLAMTime-'//trim(adjustl(numproc_char))//'-'//trim(adjustl(numtrac_char)),par%comm)
  call t_prf('HommeTime', par%comm)
  if(par%masterproc) print *,"calling t_finalizef"
  call t_finalizef()
  call haltmp("exiting program...")
end program prim_main








