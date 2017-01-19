#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_driver_base

  use cg_mod,           only: cg_t
  use derivative_mod,   only: derivative_t, derivinit
  use dimensions_mod,   only: np, nlev, nlevp, nelem, nelemd, nelemdmax, GlobalUniqueCols, qsize
  use element_mod,      only: element_t, allocate_element_desc, setup_element_pointers
  use element_state,    only: timelevels
  use hybrid_mod,       only: hybrid_t
  use kinds,            only: real_kind, iulog
  use perf_mod,         only: t_startf, t_stopf
  use prim_si_ref_mod,  only: ref_state_t
  use quadrature_mod,   only: quadrature_t, test_gauss, test_gausslobatto, gausslobatto
  use reduction_mod,    only: reductionbuffer_ordered_1d_t, red_min, red_max, red_max_int, &
                              red_sum, red_sum_int, red_flops, initreductionbuffer
  use solver_mod,       only: blkjac_t
  use thread_mod,       only: nThreadsHoriz, omp_get_num_threads

#ifndef CAM
  use prim_restart_mod, only : initrestartfile
  use restart_io_mod ,  only : RestFile,readrestart
  use test_mod,         only: set_test_initial_conditions, compute_test_forcing
#endif

  implicit none

  private
  public :: prim_init1, prim_init2, prim_finalize, smooth_topo_datasets

  type (cg_t), allocatable  :: cg(:)              ! conjugate gradient struct (nthreads)
  type (quadrature_t)   :: gp                     ! element GLL points
  type (ReductionBuffer_ordered_1d_t), save :: red   ! reduction buffer               (shared)
  type (derivative_t), public  :: deriv1

contains

  subroutine prim_init1(elem, par, dom_mt, Tl)

    use thread_mod,      only : nthreads, omp_get_thread_num
    use control_mod,     only : runtype, restartfreq, integration, topology, partmethod, use_semi_lagrange_transport
    use prim_state_mod,  only : prim_printstate_init
    use namelist_mod,    only : readnl
    use mesh_mod,        only : MeshUseMeshFile
    use time_mod,        only : nmax, time_at, timelevel_init, timelevel_t
    use mass_matrix_mod, only : mass_matrix
    use cube_mod,        only : cubeedgecount , cubeelemcount, cubetopology
    use mesh_mod,        only : MeshSetCoordinates, MeshUseMeshFile, MeshCubeTopology, MeshCubeElemCount, MeshCubeEdgeCount
    use cube_mod,        only : cube_init_atomic, set_corner_coordinates, assign_node_numbers_to_elem
    use metagraph_mod,   only : metavertex_t, metaedge_t, localelemcount, initmetagraph, printmetavertex
    use gridgraph_mod,   only : gridvertex_t, gridedge_t, allocate_gridvertex_nbrs, deallocate_gridvertex_nbrs
    use schedtype_mod,   only : schedule
    use schedule_mod,    only : genEdgeSched,  PrintSchedule
    use prim_advection_mod, only: prim_advec_init1
    use prim_advance_mod,   only: prim_advance_init

#ifdef TRILINOS
    use prim_implicit_mod, only : prim_implicit_init
#endif
    use parallel_mod, only : iam, parallel_t, syncmp, abortmp, global_shared_buf, nrepro_vars
#ifdef _MPI
    use parallel_mod, only : mpiinteger_t, mpireal_t, mpi_max, mpi_sum, haltmp
#endif
    use metis_mod, only : genmetispart
    use spacecurve_mod, only : genspacepart
    use dof_mod, only : global_dof, CreateUniqueIndex, SetElemOffset
    use params_mod, only : SFCURVE
    use domain_mod, only : domain1d_t, decompose
    use physical_constants, only : dd_pi
    use bndry_mod, only : sort_neighbor_buffer_mapping
#ifndef CAM
    use repro_sum_mod,      only: repro_sum, repro_sum_defaultopts, repro_sum_setopts
#else
    use infnan,             only: nan, assignment(=)
    use shr_reprosum_mod,   only: repro_sum => shr_reprosum_calc
#endif

#ifdef TRILINOS
    use prim_implicit_mod,  only : prim_implicit_init
#endif

    implicit none

    type (element_t),   pointer     :: elem(:)
    type (parallel_t),  intent(in)  :: par
    type (domain1d_t),  pointer     :: dom_mt(:)
    type (timelevel_t), intent(out) :: Tl

    type (GridVertex_t), target,allocatable :: GridVertex(:)
    type (GridEdge_t),   target,allocatable :: Gridedge(:)
    type (MetaVertex_t), target,allocatable :: MetaVertex(:)

    integer :: ii,ie, ith
    integer :: nets, nete
    integer :: nelem_edge,nedge
    integer :: nstep
    integer :: nlyr
    integer :: iMv
    integer :: err, ierr, l, j
    logical, parameter :: Debug = .FALSE.

    real(kind=real_kind), allocatable :: aratio(:,:)
    real(kind=real_kind) :: area(1),xtmp

    integer  :: i
    integer,allocatable :: TailPartition(:)
    integer,allocatable :: HeadPartition(:)

    integer total_nelem
    real(kind=real_kind) :: approx_elements_per_task
    integer :: n_domains

#ifndef CAM
    logical :: repro_sum_use_ddpdd, repro_sum_recompute
    real(kind=real_kind) :: repro_sum_rel_diff_max
#endif

    ! =====================================
    ! Read in model control information
    ! =====================================
    ! cam readnl is called in spmd_dyn (needed prior to mpi_init)
#ifndef CAM
    call readnl(par)
    if (MeshUseMeshFile) then
       total_nelem = MeshCubeElemCount()
    else
       total_nelem = CubeElemCount()
    end if

    approx_elements_per_task = dble(total_nelem)/dble(par%nprocs)
    if  (approx_elements_per_task < 1.0D0) then
       if(par%masterproc) print *,"number of elements=", total_nelem
       if(par%masterproc) print *,"number of procs=", par%nprocs
       call abortmp('There is not enough parallelism in the job, that is, there is less than one elements per task.')
    end if
    ! ====================================
    ! initialize reproducible sum module
    ! ====================================
    call repro_sum_defaultopts(                           &
       repro_sum_use_ddpdd_out=repro_sum_use_ddpdd,       &
       repro_sum_rel_diff_max_out=repro_sum_rel_diff_max, &
       repro_sum_recompute_out=repro_sum_recompute        )
    call repro_sum_setopts(                              &
       repro_sum_use_ddpdd_in=repro_sum_use_ddpdd,       &
       repro_sum_rel_diff_max_in=repro_sum_rel_diff_max, &
       repro_sum_recompute_in=repro_sum_recompute,       &
       repro_sum_master=par%masterproc,                      &
       repro_sum_logunit=6                           )
       if(par%masterproc) print *, "Initialized repro_sum"
#endif

#ifndef CAM
    if (par%masterproc) then
       ! =============================================
       ! Compute total simulated time...
       ! =============================================
       write(iulog,*)""
       write(iulog,*)" total simulated time = ",Time_at(nmax)
       write(iulog,*)""
       ! =============================================
       ! Perform Gauss/Gauss Lobatto tests...
       ! =============================================
       call test_gauss(np)
       call test_gausslobatto(np)
    end if
#endif
    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call derivinit(deriv1)

    ! ===============================================================
    ! Allocate and initialize the graph (array of GridVertex_t types)
    ! ===============================================================

    if (topology=="cube") then

       if (par%masterproc) then
          write(iulog,*)"creating cube topology..."
       end if

       if (MeshUseMeshFile) then
           nelem = MeshCubeElemCount()
           nelem_edge = MeshCubeEdgeCount()
       else
           nelem      = CubeElemCount()
           nelem_edge = CubeEdgeCount()
       end if

       ! we want to exit elegantly when we are using too many processors.
       if (nelem < par%nprocs) then
          call abortmp('Error: too many MPI tasks. set dyn_npes <= nelem')
       end if

       allocate(GridVertex(nelem))
       allocate(GridEdge(nelem_edge))

       do j =1,nelem
          call allocate_gridvertex_nbrs(GridVertex(j))
       end do

       if (MeshUseMeshFile) then
           if (par%masterproc) then
               write(iulog,*) "Set up grid vertex from mesh..."
           end if
           call MeshCubeTopology(GridEdge, GridVertex)
       else
           call CubeTopology(GridEdge,GridVertex)
        end if

       if(par%masterproc)       write(iulog,*)"...done."
    end if
    if(par%masterproc) write(iulog,*)"total number of elements nelem = ",nelem

    !DBG if(par%masterproc) call PrintGridVertex(GridVertex)

    if(partmethod .eq. SFCURVE) then
       if(par%masterproc) write(iulog,*)"partitioning graph using SF Curve..."
       call genspacepart(GridEdge,GridVertex)
    else
        if(par%masterproc) write(iulog,*)"partitioning graph using Metis..."
       call genmetispart(GridEdge,GridVertex)
    endif

    ! ===========================================================
    ! given partition, count number of local element descriptors
    ! ===========================================================
    allocate(MetaVertex(1))
    allocate(Schedule(1))

    nelem_edge=SIZE(GridEdge)

    allocate(TailPartition(nelem_edge))
    allocate(HeadPartition(nelem_edge))
    do i=1,nelem_edge
       TailPartition(i)=GridEdge(i)%tail%processor_number
       HeadPartition(i)=GridEdge(i)%head%processor_number
    enddo

    ! ====================================================
    !  Generate the communication graph
    ! ====================================================
    call initMetaGraph(iam,MetaVertex(1),GridVertex,GridEdge)

    nelemd = LocalElemCount(MetaVertex(1))
    if(par%masterproc .and. Debug) then 
        call PrintMetaVertex(MetaVertex(1))
    endif

    if(nelemd .le. 0) then
       call abortmp('Not yet ready to handle nelemd = 0 yet' )
       stop
    endif
#ifdef _MPI
    call mpi_allreduce(nelemd,nelemdmax,1,MPIinteger_t,MPI_MAX,par%comm,ierr)
#else
    nelemdmax=nelemd
#endif

    if (nelemd>0) then
       allocate(elem(nelemd))
       call setup_element_pointers(elem)
       call allocate_element_desc(elem)
    endif

    ! ====================================================
    !  Generate the communication schedule
    ! ====================================================

    call genEdgeSched(elem,iam,Schedule(1),MetaVertex(1))

    allocate(global_shared_buf(nelemd,nrepro_vars))
    global_shared_buf=0.0_real_kind
    !  nlyr=edge3p1%nlyr
    !  call MessageStats(nlyr)
    !  call testchecksum(par,GridEdge)

    ! ========================================================
    ! load graph information into local element descriptors
    ! ========================================================

    !  do ii=1,nelemd
    !     elem(ii)%vertex = MetaVertex(iam)%members(ii)
    !  enddo

    call syncmp(par)

    ! =================================================================
    ! Set number of domains (for 'decompose') equal to number of threads
    !  for OpenMP across elements, equal to 1 for OpenMP within element
    ! =================================================================
    n_domains = min(Nthreads,nelemd)

    ! =================================================================
    ! Initialize shared boundary_exchange and reduction buffers
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'init shared boundary_exchange buffers'
    call InitReductionBuffer(red,3*nlev,n_domains)
    call InitReductionBuffer(red_sum,5)
    call InitReductionBuffer(red_sum_int,1)
    call InitReductionBuffer(red_max,1)
    call InitReductionBuffer(red_max_int,1)
    call InitReductionBuffer(red_min,1)
    call initReductionBuffer(red_flops,1)

    gp=gausslobatto(np)  ! GLL points

    if (topology=="cube") then
       if(par%masterproc) write(iulog,*) "initializing cube elements..."
       if (MeshUseMeshFile) then
           call MeshSetCoordinates(elem)
       else
           do ie=1,nelemd
               call set_corner_coordinates(elem(ie))
           end do
           call assign_node_numbers_to_elem(elem, GridVertex)
       end if
       do ie=1,nelemd
          call cube_init_atomic(elem(ie),gp%points)
       enddo
    end if

    ! =================================================================
    ! Initialize mass_matrix
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'running mass_matrix'
    call mass_matrix(par,elem)
    allocate(aratio(nelemd,1))

    if (topology=="cube") then
       area = 0
       do ie=1,nelemd
          aratio(ie,1) = sum(elem(ie)%mp(:,:)*elem(ie)%metdet(:,:))
       enddo
       call repro_sum(aratio, area, nelemd, nelemd, 1, commid=par%comm)
       area(1) = 4*dd_pi/area(1)  ! ratio correction
       deallocate(aratio)
       if (par%masterproc) &
            write(iulog,'(a,f20.17)') " re-initializing cube elements: area correction=",area(1)

       do ie=1,nelemd
          call cube_init_atomic(elem(ie),gp%points,area(1))
       enddo
    end if

    if(par%masterproc) write(iulog,*) 're-running mass_matrix'
    call mass_matrix(par,elem)

    ! =================================================================
    ! Determine the global degree of freedome for each gridpoint
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'running global_dof'
    call global_dof(par,elem)

    ! =================================================================
    ! Create Unique Indices
    ! =================================================================

    do ie=1,nelemd
       call CreateUniqueIndex(elem(ie)%GlobalId,elem(ie)%gdofP,elem(ie)%idxP)
    enddo

    call SetElemOffset(par,elem, GlobalUniqueCols)

    call prim_printstate_init(par)
    ! Initialize output fields for plotting...

    ! initialize flux terms to 0

    do ie=1,nelemd
       elem(ie)%derived%FM=0.0
       elem(ie)%derived%FQ=0.0
       elem(ie)%derived%FQps=0.0
       elem(ie)%derived%FT=0.0

       elem(ie)%accum%Qvar=0
       elem(ie)%accum%Qmass=0
       elem(ie)%accum%Q1mass=0

       elem(ie)%derived%Omega_p=0
       elem(ie)%state%dp3d=0

    enddo

    ! ==========================================================
    !  This routines initalizes a Restart file.  This involves:
    !      I)  Setting up the MPI datastructures
    ! ==========================================================
#ifndef CAM
    if(restartfreq > 0 .or. runtype>=1)  then
       call initRestartFile(elem(1)%state,par,RestFile)
    endif
#endif
    !DBG  write(iulog,*) 'prim_init: after call to initRestartFile'

    deallocate(GridEdge)
    do j =1,nelem
       call deallocate_gridvertex_nbrs(GridVertex(j))
    end do
    deallocate(GridVertex)

    do j = 1, MetaVertex(1)%nmembers
       call deallocate_gridvertex_nbrs(MetaVertex(1)%members(j))
    end do
    deallocate(MetaVertex)
    deallocate(TailPartition)
    deallocate(HeadPartition)

    n_domains = min(Nthreads,nelemd)
    nthreads = n_domains

    ! =====================================
    ! Set number of threads...
    ! =====================================
    if(par%masterproc) then
       write(iulog,*) "Main:NThreads=",NThreads
       write(iulog,*) "Main:n_domains = ",n_domains
    endif

    allocate(dom_mt(0:n_domains-1))
    do ith=0,n_domains-1
       dom_mt(ith)=decompose(1,nelemd,n_domains,ith)
    end do
    ith=0
    nets=1
    nete=nelemd
    ! set the actual number of threads which will be used in the horizontal
    nThreadsHoriz = n_domains

    allocate(cg(0:n_domains-1))
    call prim_advance_init(par,elem,integration)
#ifdef TRILINOS
    call prim_implicit_init(par, elem)
#endif
    call Prim_Advec_Init1(par, elem,n_domains)


    if ( use_semi_lagrange_transport) then
      call sort_neighbor_buffer_mapping(par, elem,1,nelemd)
    end if

    call TimeLevel_init(tl)
    if(par%masterproc) write(iulog,*) 'end of prim_init'

  end subroutine prim_init1

  !_____________________________________________________________________
  subroutine prim_init2(elem, hybrid, nets, nete, tl, hvcoord)

    use control_mod,          only: runtype, integration, test_case, &
                                    debug_level, vfile_int, vform, vfile_mid, &
                                    topology,columnpackage, rsplit, qsplit, rk_stage_user,&
                                    sub_case, limiter_option, nu, nu_q, nu_div, tstep_type, hypervis_subcycle, &
                                    hypervis_subcycle_q, moisture, use_moisture
    use global_norms_mod,     only: test_global_integral, print_cfl
    use hybvcoord_mod,        only: hvcoord_t
    use parallel_mod,         only: parallel_t, haltmp, syncmp, abortmp
    use prim_state_mod,       only: prim_printstate, prim_diag_scalars
    use prim_si_ref_mod,      only: prim_si_refstate_init, prim_set_mass
    use prim_advance_mod,     only: prim_advance_init2
    use prim_advection_mod,   only: prim_advec_init2
    use solver_init_mod,      only: solver_init2
    use time_mod,             only: timelevel_t, tstep, phys_tscale, timelevel_init, nendstep, smooth, nsplit, TimeLevel_Qdp
    use thread_mod,           only: nthreads

#ifndef CAM
    use control_mod,          only: pertlim
#endif

#ifdef TRILINOS
    use prim_derived_type_mod ,only : derived_type, initialize
    use, intrinsic :: iso_c_binding
#endif

    type (element_t),   intent(inout) :: elem(:)
    type (hybrid_t),    intent(in)    :: hybrid
    type (TimeLevel_t), intent(inout) :: tl       ! time level struct
    type (hvcoord_t),   intent(inout) :: hvcoord  ! hybrid vertical coordinate struct
    integer,            intent(in)    :: nets     ! starting thread element number (private)
    integer,            intent(in)    :: nete     ! ending thread element number   (private)

    ! ==================================
    ! Local variables
    ! ==================================

    real (kind=real_kind) :: dt                 ! timestep

    ! variables used to calculate CFL
    real (kind=real_kind) :: dtnu               ! timestep*viscosity parameter
    real (kind=real_kind) :: dt_dyn_vis         ! viscosity timestep used in dynamics
    real (kind=real_kind) :: dt_tracer_vis      ! viscosity timestep used in tracers

    real (kind=real_kind) :: dp
    real (kind=real_kind) :: ps(np,np)          ! surface pressure

    character(len=80)     :: fname
    character(len=8)      :: njusn
    character(len=4)      :: charnum

    real (kind=real_kind) :: Tp(np)             ! transfer function

    integer :: simday
    integer :: i,j,k,ie,iptr,t,q
    integer :: ierr
    integer :: nfrc
    integer :: n0_qdp

#ifdef TRILINOS
     integer :: lenx
    real (c_double) ,allocatable, dimension(:) :: xstate(:)
    ! state_object is a derived data type passed thru noxinit as a pointer
    type(derived_type) ,target         :: state_object
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object

    type(derived_type) ,target         :: pre_object
    type(derived_type) ,pointer        :: pptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_pre

    type(derived_type) ,target         :: jac_object
    type(derived_type) ,pointer        :: jptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_jac

!    type(element_t)                    :: pc_elem(size(elem))
!    type(element_t)                    :: jac_elem(size(elem))

    logical :: compute_diagnostics
    integer :: qn0
    real (kind=real_kind) :: eta_ave_w

  interface
    subroutine noxinit(vectorSize,vector,comm,v_container,p_container,j_container) &
        bind(C,name='noxinit')
    use ,intrinsic :: iso_c_binding
      integer(c_int)                :: vectorSize,comm
      real(c_double)  ,dimension(*) :: vector
      type(c_ptr)                   :: v_container
      type(c_ptr)                   :: p_container  !precon ptr
      type(c_ptr)                   :: j_container  !analytic jacobian ptr
    end subroutine noxinit

  end interface
#endif

    if (topology == "cube") then
       call test_global_integral(elem, hybrid,nets,nete)
    end if

    ! should we assume Q(:,:,:,1) has water vapor:
    use_moisture = ( moisture /= "dry") 
    if (qsize<1) use_moisture = .false.  


    ! compute most restrictive dt*nu for use by variable res viscosity:
    ! compute timestep seen by viscosity operator:
    dt_dyn_vis = tstep
    if (qsplit>1 .and. tstep_type == 1) then
       ! tstep_type==1: RK2 followed by LF.  internal LF stages apply viscosity at 2*dt
       dt_dyn_vis = 2*tstep
    endif
    dt_tracer_vis=tstep*qsplit
    
    ! compute most restrictive condition:
    ! note: dtnu ignores subcycling
    dtnu=max(dt_dyn_vis*max(nu,nu_div), dt_tracer_vis*nu_q)
    ! compute actual viscosity timesteps with subcycling
    dt_tracer_vis = dt_tracer_vis/hypervis_subcycle_q
    dt_dyn_vis = dt_dyn_vis/hypervis_subcycle

#ifdef TRILINOS

      lenx=(np*np*nlev*3 + np*np*1)*(nete-nets+1)  ! 3 3d vars plus 1 2d vars
      allocate(xstate(lenx))
      xstate(:) = 0d0
      compute_diagnostics = .false.
      qn0 = -1 ! dry case for testing right now
      eta_ave_w = 1d0 ! divide by qsplit for mean flux interpolation

      call initialize(state_object, lenx, elem, hvcoord, compute_diagnostics, &
        qn0, eta_ave_w, hybrid, deriv1, tstep, tl, nets, nete)

      call initialize(pre_object, lenx, elem, hvcoord, compute_diagnostics, &
        qn0, eta_ave_w, hybrid, deriv1, tstep, tl, nets, nete)

      call initialize(jac_object, lenx, elem, hvcoord, .false., &
        qn0, eta_ave_w, hybrid, deriv1, tstep, tl, nets, nete)

!      pc_elem = elem
!      jac_elem = elem

      fptr => state_object
      c_ptr_to_object =  c_loc(fptr)
      pptr => state_object
      c_ptr_to_pre =  c_loc(pptr)
      jptr => state_object
      c_ptr_to_jac =  c_loc(jptr)

      call noxinit(size(xstate), xstate, 1, c_ptr_to_object, c_ptr_to_pre, c_ptr_to_jac)

#endif


#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
    if (hybrid%ithr==0) then
       call syncmp(hybrid%par)
    end if
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif

    if (topology /= "cube") then
       call abortmp('Error: only cube topology supported for primaitve equations')
    endif

#ifndef CAM

    ! initialize dynamics solver before setting initial conditions
    call prim_advance_init2(elem, nets, nete, hybrid, hvcoord)

    ! =================================
    ! HOMME stand alone initialization
    ! =================================

    if(runtype >= 1) then

       ! ===========================================================
       ! runtype==1   Exact Restart
       ! runtype==2   Initial run, but take inital condition from Restart file
       ! ===========================================================

       if (hybrid%masterthread) then
          write(iulog,*) 'runtype: RESTART of primitive equations'
       end if

       call ReadRestart(elem,hybrid%ithr,nets,nete,tl)

       ! scale PS to achieve prescribed dry mass
       if (runtype /= 1) call prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)

       if (runtype==2) then
          ! copy prognostic variables: tl%n0 into tl%nm1
          do ie=nets,nete
             elem(ie)%state%v(:,:,:,:,tl%nm1) = elem(ie)%state%v(:,:,:,:,tl%n0)
             elem(ie)%state%T(:,:,:,tl%nm1)   = elem(ie)%state%T(:,:,:,tl%n0)
             elem(ie)%state%ps_v(:,:,tl%nm1)  = elem(ie)%state%ps_v(:,:,tl%n0)
          enddo
       endif ! runtype==2

    else

      ! runtype=0: initial run
      if (hybrid%masterthread) write(iulog,*) ' runtype: initial run'
      call set_test_initial_conditions(elem,deriv1,hybrid,hvcoord,tl,nets,nete)
      if (hybrid%masterthread) write(iulog,*) '...done'

      ! scale PS to achieve prescribed dry mass
      call prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)

      do ie=nets,nete
        ! set perlim in ctl_nl namelist for temperature field initial perturbation
        elem(ie)%state%T=elem(ie)%state%T * (1.0_real_kind + pertlim)
      enddo

    endif !runtype

#endif
!$OMP MASTER
    tl%nstep0=2                   ! compute diagnostics starting with step 2 if LEAPFROG
    tl%nstep0=1                   ! compute diagnostics starting with step 1 if RK
    if (runtype==1) then
       tl%nstep0=tl%nstep+1       ! compute diagnostics after 1st step, leapfrog or RK
    endif
    if (runtype==2) then
       ! branch run
       ! reset time counters to zero since timestep may have changed
       nEndStep = nEndStep-tl%nstep ! used by standalone HOMME.  restart code set this to nmax + tl%nstep
       tl%nstep=0
    endif
!$OMP END MASTER
!$OMP BARRIER


    ! For new runs, and branch runs, convert state variable to (Qdp)
    ! because initial conditon reads in Q, not Qdp
    ! restart runs will read dpQ from restart file
    ! need to check what CAM does on a branch run
    if (runtype==0 .or. runtype==2) then
       do ie=nets,nete
          elem(ie)%derived%omega_p(:,:,:) = 0D0
       end do
       do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k, t, q, i, j, dp)
#endif
          do k=1,nlev    !  Loop inversion (AAM)
             do q=1,qsize
                do i=1,np
                   do j=1,np
                      dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                           ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,tl%n0)
                      
                      elem(ie)%state%Qdp(i,j,k,q,1)=elem(ie)%state%Q(i,j,k,q)*dp
                      elem(ie)%state%Qdp(i,j,k,q,2)=elem(ie)%state%Q(i,j,k,q)*dp
                      
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif


    if (runtype==1) then
       call TimeLevel_Qdp( tl, qsplit, n0_qdp)
       do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k, t, q, i, j, dp)
#endif
          do k=1,nlev    !  Loop inversion (AAM)
             do t=tl%n0,tl%n0
                do q=1,qsize
                   do i=1,np
                      do j=1,np
                         dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                              ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,t)
                         elem(ie)%state%Q(i,j,k,q)=elem(ie)%state%Qdp(i,j,k,q, n0_qdp)/dp
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif


    ! timesteps to use for advective stability:  tstep*qsplit and tstep
    call print_cfl(elem,hybrid,nets,nete,dtnu)

    if (hybrid%masterthread) then
       ! CAM has set tstep based on dtime before calling prim_init2(),
       ! so only now does HOMME learn the timstep.  print them out:
       write(iulog,'(a,2f9.2)') "dt_remap: (0=disabled)   ",tstep*qsplit*rsplit
       if (qsize>0) then
          write(iulog,'(a,2f9.2)') "dt_tracer (SE), per RK stage: ",tstep*qsplit,(tstep*qsplit)/(rk_stage_user-1)
       end if
       write(iulog,'(a,2f9.2)')    "dt_dyn:                  ",tstep
       write(iulog,'(a,2f9.2)')    "dt_dyn (viscosity):      ",dt_dyn_vis
       write(iulog,'(a,2f9.2)')    "dt_tracer (viscosity):   ",dt_tracer_vis


#ifdef CAM
       if (phys_tscale/=0) then
          write(iulog,'(a,2f9.2)') "CAM physics timescale:       ",phys_tscale
       endif
       write(iulog,'(a,2f9.2)') "CAM dtime (dt_phys):         ",tstep*nsplit*qsplit*max(rsplit,1)
#endif
    end if


    if (hybrid%masterthread) write(iulog,*) "initial state:"
    call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)

    call solver_init2(elem(:), deriv1)
    call prim_advec_init2(elem(:), hvcoord, hybrid)

  end subroutine prim_init2

!=======================================================================================================!

  subroutine prim_finalize()

#ifdef TRILINOS
  interface
    subroutine noxfinish() bind(C,name='noxfinish')
    use ,intrinsic :: iso_c_binding
    end subroutine noxfinish
  end interface

  call noxfinish()

#endif

    ! ==========================
    ! end of the hybrid program
    ! ==========================
  end subroutine prim_finalize

   subroutine smooth_topo_datasets(phis,sghdyn,sgh30dyn,elem,hybrid,nets,nete)
    use control_mod, only : smooth_phis_numcycle,smooth_sgh_numcycle
    use hybrid_mod, only : hybrid_t
    use bndry_mod, only : bndry_exchangev
    use derivative_mod, only : derivative_t , laplace_sphere_wk
    use viscosity_mod, only : smooth_phis
    implicit none

    integer , intent(in) :: nets,nete
    real (kind=real_kind), intent(inout)   :: phis(np,np,nets:nete)
    real (kind=real_kind), intent(inout)   :: sghdyn(np,np,nets:nete)
    real (kind=real_kind), intent(inout)   :: sgh30dyn(np,np,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    ! local
    integer :: ie
    real (kind=real_kind) :: minf

    minf=-9e9
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to PHIS"
    call smooth_phis(phis,elem,hybrid,deriv1,nets,nete,minf,smooth_phis_numcycle)

    minf=0
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to SGH"
    call smooth_phis(sghdyn,elem,hybrid,deriv1,nets,nete,minf,smooth_sgh_numcycle)
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to SGH30"
    call smooth_phis(sgh30dyn,elem,hybrid,deriv1,nets,nete,minf,smooth_sgh_numcycle)

    end subroutine smooth_topo_datasets


end module prim_driver_base



