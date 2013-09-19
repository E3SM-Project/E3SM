#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if(!defined ELEMENT_OPENMP)
#define EOMP_OFF 1
#else
#define EOMP_OFF 0
#endif

program prim_main

  use parallel_mod,     only: parallel_t,initmp,haltmp
  use hybrid_mod,       only: hybrid_t
  use thread_mod,       only: nthreads,omp_get_thread_num
  use time_mod,         only: tstep,nendstep,timelevel_t
  use dimensions_mod,   only: nelemd
  use control_mod,      only: restartfreq,vfile_mid,vfile_int,runtype,integration,tstep_type
  use domain_mod,       only: domain1d_t
  use element_mod,      only: element_t
  use spelt_mod,        only: spelt_struct
  use common_io_mod,    only: output_dir
  use common_movie_mod, only: nextoutputstep
  use restart_io_mod,   only: restartheader_t, writerestart
  use hybrid_mod,       only: hybrid_create
  use perf_mod,         only: t_initf,t_prf,t_finalizef,t_startf,t_stopf ! _EXTERNAL
  use fvm_control_volume_mod, only: fvm_struct

#ifdef _PRIM
  ! Modules specific to preqx exectuable
  use prim_driver_mod, only: prim_init1,prim_init2,prim_run,prim_run_subcycle,prim_finalize,leapfrog_bootstrap
  use hybvcoord_mod,   only: hvcoord_t,hvcoord_init
#endif

#ifdef PIO_INTERP
  ! Modules for interpolated output
  use interp_movie_mod,       only: interp_movie_output,interp_movie_finish,interp_movie_init
  use interpolate_driver_mod, only: interpolate_driver
#else
  ! Modules for native grid output
  use prim_movie_mod, only: prim_movie_output,prim_movie_finish,prim_movie_init
#endif

  implicit none

  type (element_t),  pointer  :: elem(:)                              ! element  array pointer
  type (domain1d_t), pointer  :: dom(:)                               ! domain1d array pointer
  type (hybrid_t)             :: hybrid                               ! parallel struct for shared mem + distributed mem
  type (parallel_t)           :: par                                  ! parallel struct for distributed mem
  type (restartHeader_t)      :: restartHeader
  type (timeLevel_t)          :: tl                                   ! main time level struct
  type (hvcoord_t)            :: hvcoord                              ! hybrid vertical coordinate struct
  integer                     :: nets, nete                           ! start, end element indices
  integer                     :: ithr                                 ! thread number
  integer                     :: ierr                                 ! error code
  integer                     :: nstep                                ! step counter
  logical                     :: dir_e                                ! flag: netcdf directory exists
  logical, parameter          :: lprint = .true.                      ! print hvcoord levels

#if defined(_SPELT)
  ! Replace finite-volume structure if using SPELT transport scheme
  type (spelt_struct),pointer :: fvm(:)                               ! let fvm = spelt structure
#else
  type (fvm_struct),  pointer :: fvm(:)                               ! let fvm = finite-volume structure
#endif

#ifdef _HTRACE
  integer htype,pflag
#endif

  par = initmp()                                                      ! initialize mpi

  call t_initf  ('input.nl',LogPrint=par%masterproc,mpicom=par%comm,masterTask=par%masterproc) ! init timers
  call t_startf ('Total')                                             ! start performance timer

  call t_startf ('prim_init1')                                        
  call prim_init1(elem,  fvm, par,dom,tl)                             ! init mesh and element array
  call t_stopf  ('prim_init1')                                        

#ifdef _HTRACE
  htype = 19; pflag = 0; call TRACE_INIT(htype,pflag)
#endif

  call display_proc_and_thread_ids(par)                               ! display ids for first 100 procs

  call initialize_vertical_coords(hvcoord)                            ! set vertical coords (CAM inits hvcoord externally)

  call interpolate_and_exit(hybrid)                                   ! interpolate netcdf file, if runtype<0

  if(par%masterproc) print *,"Initializing Primitive Equation."

  !$OMP if(EOMP_OFF) PARALLEL DEFAULT(SHARED), PRIVATE(ithr,hybrid)

  ithr   = omp_get_thread_num()
  hybrid = hybrid_create(par,ithr,nthreads)
  nets   = dom(ithr)%start
  nete   = dom(ithr)%end

  call t_startf('prim_init2')                                         
  call prim_init2(elem,fvm,hybrid,nets,nete,tl,hvcoord)               ! set initial state
  call t_stopf('prim_init2')

  !$OMP if(EOMP_OFF) END PARALLEL

  ithr   = omp_get_thread_num()
  hybrid = hybrid_create(par,ithr,1)
  call ensure_output_dir_exists(hybrid)                               ! ensure output dir exists

  call initialize_output_files()                                      ! init output files. write first output

  if(integration=='semi_imp' .and. runtype/=1 ) then
    call leapfrog_bootstrap(elem,hybrid,1,nelemd,tstep,tl,hvcoord)    ! use bootstrap to start semi-implicit runs
  endif

  if(par%masterproc) print *,"Entering main timestepping loop"

  do while(tl%nstep < nEndStep)                                       ! loop until end step

     !$OMP if(EOMP_OFF) PARALLEL DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid)

     ithr   = omp_get_thread_num()
     hybrid = hybrid_create(par,ithr,nthreads)
     nets   = dom(ithr)%start
     nete   = dom(ithr)%end

     nstep  = nextoutputstep(tl)                                      ! get next output step

     do while(tl%nstep < nstep)                                       ! loop until next output step
        call t_startf('prim_run')                                     

        if (tstep_type>0) then                                        ! advance with    subcycling
           call prim_run_subcycle(elem,fvm,hybrid,nets,nete,tstep,tl,hvcoord,1)
        else                                                          ! advance without subcycling
           call prim_run(elem,hybrid,nets,nete,tstep,tl,hvcoord,"leapfrog")
        endif

        call t_stopf('prim_run')                                      
     end do

     !$OMP if(EOMP_OFF) END PARALLEL

     call write_data_to_file()                                        ! output interp or native data

     if((restartfreq>0) .and. (MODULO(tl%nstep,restartfreq)==0)) then 
       call writeRestart(elem,ithr,1,nelemd,tl)                       ! write restart files occasionally
     endif
  end do                                                              ! exit main loop

  if(par%masterproc) print *,"Finished main timestepping loop",tl%nstep

  call prim_finalize(hybrid)                                          ! finalize mpi

  call close_data_files()                                             ! close files, flush output

  call t_stopf('Total')                                               ! halt performance timer

  if(par%masterproc) print *,"writing timing data"
  call t_prf('HommeTime', par%comm)                                   ! write timer data to file

  if(par%masterproc) print *,"calling t_finalizef"
  call t_finalizef()

  call haltmp("exiting program...")                                   ! exit

  contains

!_____________________________________________________________________
! Interpolate data between grids if runtype <0 and PIO_INTERP is active

subroutine interpolate_and_exit(hybrid)
  type (hybrid_t) :: hybrid

  if(runtype>=0) return

#ifdef PIO_INTERP
  call interpolate_driver(elem, hybrid)
  call haltmp('interpolation complete')
#endif

end subroutine

!_____________________________________________________________________
! Display process and thread ids for first 100 procs

subroutine display_proc_and_thread_ids(par)

  type (parallel_t), intent(in) :: par

  integer :: ithr                                                     

  !$OMP if(EOMP_OFF) PARALLEL DEFAULT(SHARED), PRIVATE(ithr,hybrid)

  ! Display process and thread ids for first 100 procs
  if (par%rank<100) write(6,9) par%rank, omp_get_thread_num(), dom(ithr)%start, dom(ithr)%end 
9 format("process: ",i2,1x,"thread: ",i2,1x,"element limits: ",i5," - ",i5)

  !$OMP if(EOMP_OFF) END PARALLEL

end subroutine

!_____________________________________________________________________
! Ensure data output dir exists to avoid error in PIO

subroutine ensure_output_dir_exists(hybrid)

  type (hybrid_t), intent(in) :: hybrid

  if (hybrid%masterthread) then 
     open(unit=447,file=trim(output_dir) // "/output_dir_test",iostat=ierr)
     if ( ierr==0 ) then
        print *,'Directory ',output_dir, 'exists: initialing IO'
        close(447)
     else
        print *,'Error creating file in directory ',output_dir
        call haltmp("Please be sure the directory exist or specify 'output_dir' in the namelist.")
     end if
  endif

end subroutine

!_____________________________________________________________________
! Read vertical coordinate coefficients from interface and midpoint files

subroutine initialize_vertical_coords(hvcoord)

  type (hvcoord_t), intent(inout):: hvcoord

  integer :: ithr
  type (hybrid_t) :: hybrid                               

  ithr    = omp_get_thread_num()
  hybrid  = hybrid_create(par, ithr, 1)
  hvcoord = hvcoord_init(vfile_mid, vfile_int, lprint, hybrid%masterthread, ierr)
  if (ierr /= 0) call haltmp("error in hvcoord_init")

end subroutine

!_____________________________________________________________________
! Initialize data files & output initial state for NEW runs (not restarts or branch runs)

subroutine initialize_output_files()

#ifdef PIO_INTERP
  call interp_movie_init( elem, hybrid, 1, nelemd, hvcoord, tl )
  if (runtype==0) call interp_movie_output(elem, tl, hybrid, 0d0, 1, nelemd,fvm=fvm, hvcoord=hvcoord)
#else
  call prim_movie_init( elem, par, hvcoord, tl )
  if (runtype==0) call prim_movie_output(elem, tl, hvcoord, hybrid, 1,nelemd, fvm)
#endif

end subroutine

!_____________________________________________________________________
subroutine close_data_files()

if(par%masterproc) print *,"closing history files"
#ifdef PIO_INTERP
  call interp_movie_finish
#else
  call prim_movie_finish
#endif

end subroutine

!_____________________________________________________________________
subroutine write_data_to_file()

  ithr   = omp_get_thread_num()
  hybrid = hybrid_create(par,ithr,1)
#ifdef PIO_INTERP
  call interp_movie_output(elem,tl,hybrid,0d0,1,nelemd,fvm=fvm,hvcoord=hvcoord)
#else
  call prim_movie_output(elem,tl,hvcoord,hybrid,1,nelemd,fvm)
#endif
end subroutine

end program prim_main

