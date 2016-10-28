#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"file: ",__FILE__," line: ",__LINE__," ithr: ",hybrid%ithr
#define _DBG_
module prim_driver_mod

  use cg_mod,           only: cg_t
  use derivative_mod,   only: derivative_t
  use dimensions_mod,   only: np, nlev, nlevp, nelem, nelemd, nelemdmax, GlobalUniqueCols, ntrac, qsize, nc,nhc
  use element_mod,      only: element_t, timelevels,  allocate_element_desc
  use fvm_mod,          only: fvm_init1,fvm_init2, fvm_init3
  use fvm_control_volume_mod, only: fvm_struct
  use hybrid_mod,       only: hybrid_t
  use kinds,            only: real_kind, iulog, longdouble_kind
  use perf_mod,         only: t_startf, t_stopf
  use prim_si_ref_mod,  only: ref_state_t
  use quadrature_mod,   only: quadrature_t, test_gauss, test_gausslobatto, gausslobatto
  use reduction_mod,    only: reductionbuffer_ordered_1d_t, red_min, red_max, red_max_int, &
                              red_sum, red_sum_int, red_flops, initreductionbuffer
  use solver_mod,       only: blkjac_t
  use thread_mod,       only: nThreadsHoriz, omp_get_num_threads

#ifndef CAM
  use column_types_mod, only : ColumnModel_t
  use prim_restart_mod, only : initrestartfile
  use restart_io_mod ,  only : RestFile,readrestart
  use Manager
  use test_mod,         only: set_test_initial_conditions, apply_test_forcing
#endif

  implicit none

  private
  public :: prim_init1, prim_init2 , prim_run_subcycle, prim_finalize
  public :: smooth_topo_datasets

  type (cg_t), allocatable  :: cg(:)              ! conjugate gradient struct (nthreads)
  type (quadrature_t)   :: gp                     ! element GLL points
  real(kind=longdouble_kind)  :: fvm_corners(nc+1)     ! fvm cell corners on reference element
  real(kind=longdouble_kind)  :: fvm_points(nc)     ! fvm cell centers on reference element

#ifndef CAM
  type (ColumnModel_t), allocatable :: cm(:) ! (nthreads)
#endif
  type (ReductionBuffer_ordered_1d_t), save :: red   ! reduction buffer               (shared)

contains

  subroutine prim_init1(elem, fvm, par, dom_mt, Tl)

    ! --------------------------------
    use thread_mod, only : nthreads, omp_get_thread_num, vert_num_threads
    ! --------------------------------
    use control_mod, only : runtype, restartfreq, integration, topology, &
         partmethod, while_iter, use_semi_lagrange_transport
    ! --------------------------------
    use prim_state_mod, only : prim_printstate_init
    ! --------------------------------
    use namelist_mod, only : readnl
    ! --------------------------------
    use mesh_mod, only : MeshUseMeshFile
    ! --------------------------------
    use time_mod, only : nmax, time_at, timelevel_init, timelevel_t
    ! --------------------------------
    use element_mod, only : setup_element_pointers
    ! --------------------------------
    use mass_matrix_mod, only : mass_matrix
    ! --------------------------------
    use cube_mod,  only : cubeedgecount , cubeelemcount, cubetopology
    ! --------------------------------
    use mesh_mod, only : MeshSetCoordinates, MeshUseMeshFile, MeshCubeTopology, &
         MeshCubeElemCount, MeshCubeEdgeCount
    use cube_mod, only : cube_init_atomic, rotation_init_atomic, set_corner_coordinates, assign_node_numbers_to_elem
    ! --------------------------------
    use metagraph_mod, only : metavertex_t, metaedge_t, localelemcount, initmetagraph, printmetavertex
    ! --------------------------------
    use derivative_mod, only : allocate_subcell_integration_matrix
    ! --------------------------------
    use gridgraph_mod, only : gridvertex_t, gridedge_t, allocate_gridvertex_nbrs, deallocate_gridvertex_nbrs
    ! --------------------------------
    use schedtype_mod, only : schedule
    ! --------------------------------
    use schedule_mod, only : genEdgeSched,  PrintSchedule
    ! --------------------------------
    use prim_advection_mod, only: prim_advec_init1
    ! --------------------------------
    use prim_advance_mod, only: prim_advance_init
    ! --------------------------------
#ifdef TRILINOS
    use prim_implicit_mod, only : prim_implicit_init
#endif
    ! --------------------------------
    use parallel_mod, only : iam, parallel_t, syncmp, abortmp, global_shared_buf, nrepro_vars
#ifdef _MPI
    use parallel_mod, only : mpiinteger_t, mpireal_t, mpi_max, mpi_sum, haltmp
#endif
    ! --------------------------------
    use metis_mod, only : genmetispart
    ! --------------------------------
    use spacecurve_mod, only : genspacepart
    ! --------------------------------
    use dof_mod, only : global_dof, CreateUniqueIndex, SetElemOffset
    ! --------------------------------
    use params_mod, only : SFCURVE
    ! --------------------------------
    use domain_mod, only : domain1d_t, decompose
    ! --------------------------------
    use physical_constants, only : dd_pi
    ! --------------------------------
    use bndry_mod, only : sort_neighbor_buffer_mapping
    ! --------------------------------
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
    type (fvm_struct),  pointer     :: fvm(:)
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
    character(len=80) rot_type   ! cube edge rotation type

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
    ! ====================================
    ! Set cube edge rotation type for model
    ! unnecessary complication here: all should
    ! be on the same footing. RDL
    ! =====================================
    rot_type="contravariant"

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
#ifndef CAM
       call ManagerInit()
#endif
    endif

    if (ntrac>0) then
       allocate(fvm(nelemd))
    else
       ! Even if fvm not needed, still desirable to allocate it as empty
       ! so it can be passed as a (size zero) array rather than pointer.
       allocate(fvm(0))
    end if

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

    ! fvm nodes are equally spaced in alpha/beta
    ! HOMME with equ-angular gnomonic projection maps alpha/beta space
    ! to the reference element via simple scale + translation
    ! thus, fvm nodes in reference element [-1,1] are a tensor product of
    ! array 'fvm_corners(:)' computed below:
    xtmp=nc
    do i=1,nc+1
       fvm_corners(i)= 2*(i-1)/xtmp - 1  ! [-1,1] including end points
    end do
    do i=1,nc
       fvm_points(i)= ( fvm_corners(i)+fvm_corners(i+1) ) /2
    end do

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
          !call rotation_init_atomic(elem(ie),rot_type)
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

    do ie=1,nelemd
       elem(ie)%idxV=>elem(ie)%idxP
    end do

    !JMD call PrintDofP(elem)
    !JMD call PrintDofV(elem)

    call prim_printstate_init(par)
    ! Initialize output fields for plotting...

    while_iter = 0

    ! initialize flux terms to 0

    do ie=1,nelemd
       elem(ie)%derived%FM=0.0
       elem(ie)%derived%FQ=0.0
       elem(ie)%derived%FQps=0.0
       elem(ie)%derived%FT=0.0
       elem(ie)%derived%pecnd=0.0

       elem(ie)%accum%Qvar=0
       elem(ie)%accum%Qmass=0
       elem(ie)%accum%Q1mass=0

       elem(ie)%derived%Omega_p=0
       elem(ie)%state%dp3d=0

#ifdef CAM
       elem(ie)%derived%etadot_prescribed = nan
       elem(ie)%derived%u_met = nan
       elem(ie)%derived%v_met = nan
       elem(ie)%derived%dudt_met = nan
       elem(ie)%derived%dvdt_met = nan
       elem(ie)%derived%T_met = nan
       elem(ie)%derived%dTdt_met = nan
       elem(ie)%derived%ps_met = nan
       elem(ie)%derived%dpsdt_met = nan
       elem(ie)%derived%nudge_factor = nan

       elem(ie)%derived%Utnd=0.D0
       elem(ie)%derived%Vtnd=0.D0
       elem(ie)%derived%Ttnd=0.D0
#endif
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
#ifndef CAM
    allocate(cm(0:n_domains-1))
#endif
    allocate(cg(0:n_domains-1))
    call prim_advance_init(par,elem,integration)
#ifdef TRILINOS
    call prim_implicit_init(par, elem)
#endif
    call Prim_Advec_Init1(par, elem,n_domains)
    if (ntrac>0) then
      call fvm_init1(par,elem)
    endif

    ! =======================================================
    ! Allocate memory for subcell flux calculations.
    ! =======================================================
    call allocate_subcell_integration_matrix(np, nc)

    if ( use_semi_lagrange_transport) then
      call sort_neighbor_buffer_mapping(par, elem,1,nelemd)
    end if

    call TimeLevel_init(tl)
    if(par%masterproc) write(iulog,*) 'end of prim_init'

  end subroutine prim_init1

  !_____________________________________________________________________
  subroutine prim_init2(elem, fvm, hybrid, nets, nete, tl, hvcoord)

    use control_mod,          only: runtype, integration, test_case, &
                                    debug_level, vfile_int, vform, vfile_mid, &
                                    topology,columnpackage, moisture, rsplit, qsplit, rk_stage_user,&
                                    sub_case, limiter_option, nu, nu_q, nu_div, tstep_type, hypervis_subcycle, &
                                    hypervis_subcycle_q, tracer_transport_type
    use derivative_mod,       only: derivinit, interpolate_gll2fvm_points, v2pinit
    use filter_mod,           only: filter_t, fm_filter_create, taylor_filter_create, fm_transfer, bv_transfer
    use fvm_control_volume_mod, only: n0_fvm, np1_fvm,fvm_supercycling
    use global_norms_mod,     only: test_global_integral, print_cfl
    use hybvcoord_mod,        only: hvcoord_t
    use parallel_mod,         only: parallel_t, haltmp, syncmp, abortmp
    use prim_state_mod,       only: prim_printstate, prim_diag_scalars
    use prim_si_ref_mod,      only: prim_si_refstate_init, prim_set_mass
    use prim_advection_mod,   only: prim_advec_init2, prim_advec_init_deriv, deriv
    use solver_init_mod,      only: solver_init2
    use time_mod,             only: timelevel_t, tstep, phys_tscale, timelevel_init, nendstep, smooth, nsplit, TimeLevel_Qdp
    use thread_mod,           only: nthreads

#ifndef CAM
    use column_model_mod,     only: InitColumnModel
    use control_mod,          only: pertlim                     
#endif

#ifdef TRILINOS
    use prim_derived_type_mod ,only : derived_type, initialize
    use, intrinsic :: iso_c_binding
#endif

    type (element_t),   intent(inout) :: elem(:)
    type (fvm_struct),  intent(inout) :: fvm(:)
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
        qn0, eta_ave_w, hybrid, deriv(hybrid%ithr), tstep, tl, nets, nete)

      call initialize(pre_object, lenx, elem, hvcoord, compute_diagnostics, &
        qn0, eta_ave_w, hybrid, deriv(hybrid%ithr), tstep, tl, nets, nete)

      call initialize(jac_object, lenx, elem, hvcoord, .false., &
        qn0, eta_ave_w, hybrid, deriv(hybrid%ithr), tstep, tl, nets, nete)

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

    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call Prim_Advec_Init_deriv(hybrid, fvm_corners, fvm_points)

    ! ================================================
    ! fvm initialization
    ! ================================================
    if (ntrac>0) then
      call fvm_init2(elem,fvm,hybrid,nets,nete,tl)
    endif

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

    ! =================================
    ! HOMME stand alone initialization
    ! =================================

    call InitColumnModel(elem, cm(hybrid%ithr), hvcoord, hybrid, tl,nets,nete,runtype)

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
             elem(ie)%state%lnps(:,:,tl%nm1)  = elem(ie)%state%lnps(:,:,tl%n0)
          enddo
       endif ! runtype==2

    else

      ! runtype=0: initial run
      if (hybrid%masterthread) write(iulog,*) ' runtype: initial run'
      call set_test_initial_conditions(elem,deriv(hybrid%ithr),hybrid,hvcoord,tl,nets,nete,fvm)
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

    if (ntrac>0) then

      ! do it only for FVM tracers, dp_fvm field will be the AIR DENSITY
      ! should be optimize and combined with the above caculation
      do ie=nets,nete
        do k=1,nlev
          do i=1,np
            do j=1,np
              elem(ie)%derived%dp(i,j,k)=( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                                         ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,tl%n0)
            enddo
          enddo
          !write air density in dp_fvm field of FVM
          fvm(ie)%dp_fvm(1:nc,1:nc,k,n0_fvm)=interpolate_gll2fvm_points(elem(ie)%derived%dp(:,:,k),deriv(hybrid%ithr))
        enddo
      enddo
      call fvm_init3(elem,fvm,hybrid,nets,nete,n0_fvm) !boundary exchange
      do ie=nets,nete
        do i=1-nhc,nc+nhc
          do j=1-nhc,nc+nhc
            !phl is it necessary to compute psc here?
            fvm(ie)%psc(i,j) = sum(fvm(ie)%dp_fvm(i,j,:,n0_fvm)) +  hvcoord%hyai(1)*hvcoord%ps0
          enddo
        enddo
      enddo
      if (hybrid%masterthread) then
         write(iulog,*) 'FVM tracers initialized.'
      end if
    endif

    ! for restart runs, we read in Qdp for exact restart, and rederive Q
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

       if (ntrac>0) then
          write(iulog,'(a,2f9.2)') "dt_tracer (fvm)          ",tstep*qsplit*fvm_supercycling
       end if
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
    call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete, fvm)

    call solver_init2(elem(:), deriv(hybrid%ithr))
    call Prim_Advec_Init2(elem(:), hvcoord, hybrid)

  end subroutine prim_init2

!=======================================================================================================!



  subroutine prim_run_subcycle(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,nsubstep)

    !   advance dynamic variables and tracers (u,v,T,ps,Q,C) from time t to t + dt_q
    !
    !   input:
    !       tl%nm1   not used
    !       tl%n0    data at time t
    !       tl%np1   new values at t+dt_q
    !
    !   then we update timelevel pointers:
    !       tl%nm1 = tl%n0
    !       tl%n0  = tl%np1
    !   so that:
    !       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
    !       tl%n0    time t + dt_q

    use control_mod,        only: statefreq, energy_fixer, ftype, qsplit, rsplit, test_cfldep, disable_diagnostics
    use fvm_control_volume_mod, only: n0_fvm
    use hybvcoord_mod,      only: hvcoord_t
    use parallel_mod,       only: abortmp
    use prim_advance_mod,   only: applycamforcing, applycamforcing_dynamics
    use prim_state_mod,     only: prim_printstate, prim_diag_scalars, prim_energy_halftimes
    use prim_advection_mod, only: vertical_remap
    use reduction_mod,      only: parallelmax
    use time_mod,           only: TimeLevel_t, timelevel_update, timelevel_qdp, nsplit

#if USE_OPENACC
    use openacc_utils_mod,  only: copy_qdp_h2d, copy_qdp_d2h
#endif

    type (element_t) ,    intent(inout) :: elem(:)
    type (fvm_struct),    intent(inout) :: fvm(:)
    type (hybrid_t),      intent(in)    :: hybrid                       ! distributed parallel structure (shared)
    type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets                         ! starting thread element number (private)
    integer,              intent(in)    :: nete                         ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt                           ! "timestep dependent" timestep
    type (TimeLevel_t),   intent(inout) :: tl
    integer,              intent(in)    :: nsubstep                     ! nsubstep = 1 .. nsplit

    real(kind=real_kind) :: dp, dt_q, dt_remap
    real(kind=real_kind) :: dp_np1(np,np)
    integer :: ie,i,j,k,n,q,t
    integer :: n0_qdp,np1_qdp,r,nstep_end
    logical :: compute_diagnostics, compute_energy

    ! compute timesteps for tracer transport and vertical remap

    dt_q      = dt*qsplit
    dt_remap  = dt_q
    nstep_end = tl%nstep + qsplit
    if (rsplit>0) then
       dt_remap  = dt_q*rsplit   ! rsplit=0 means use eulerian code, not vert. lagrange
       nstep_end = tl%nstep + qsplit*rsplit  ! nstep at end of this routine
    endif

    ! activate energy diagnostics if using an energy fixer
    compute_energy = energy_fixer > 0

    ! activate diagnostics periodically for display to stdout
    compute_diagnostics   = .false.
    if (MODULO(nstep_end,statefreq)==0 .or. nstep_end==tl%nstep0) then
       compute_diagnostics= .true.
       compute_energy     = .true.
    endif
    if(disable_diagnostics) compute_diagnostics= .false.

    ! compute scalar diagnostics if currently active
    if (compute_diagnostics) then
      call t_startf("prim_diag_scalars")
      call prim_diag_scalars(elem,hvcoord,tl,4,.true.,nets,nete)
      call t_stopf("prim_diag_scalars")
    endif

#ifdef CAM

    ! Apply CAM Physics forcing

    !   ftype= 2: Q was adjusted by physics, but apply u,T forcing here
    !   ftype= 1: forcing was applied time-split in CAM coupling layer
    !   ftype= 0: apply all forcing here
    !   ftype=-1: do not apply forcing

    call TimeLevel_Qdp(tl, qsplit, n0_qdp)

    if (ftype==0) then
      call t_startf("ApplyCAMForcing")
      call ApplyCAMForcing(elem, fvm, hvcoord,tl%n0,n0_qdp, dt_remap,nets,nete)
      call t_stopf("ApplyCAMForcing")

    elseif (ftype==2) then
      call t_startf("ApplyCAMForcing_dynamics")
      call ApplyCAMForcing_dynamics(elem, hvcoord,tl%n0,dt_remap,nets,nete)
      call t_stopf("ApplyCAMForcing_dynamics")
    endif
#else
    ! Apply HOMME test case forcing
    call apply_test_forcing(elem,fvm,hybrid,hvcoord,tl%n0,n0_qdp,dt_remap,nets,nete)

#endif

    ! E(1) Energy after CAM forcing
    if (compute_energy) then
      call t_startf("prim_energy_halftimes")
      call prim_energy_halftimes(elem,hvcoord,tl,1,.true.,nets,nete)
      call t_stopf("prim_energy_halftimes")
    endif

    ! qmass and variance, using Q(n0),Qdp(n0)
    if (compute_diagnostics) then
      call t_startf("prim_diag_scalars")
      call prim_diag_scalars(elem,hvcoord,tl,1,.true.,nets,nete)
      call t_stopf("prim_diag_scalars")
    endif

    ! initialize dp3d from ps
    do ie=nets,nete
       do k=1,nlev
          elem(ie)%state%dp3d(:,:,k,tl%n0)=&
               ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
       enddo
    enddo


#if (USE_OPENACC)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    call t_startf("copy_qdp_h2d")
    call copy_qdp_h2d( elem , n0_qdp )
    call t_stopf("copy_qdp_h2d")
#endif

    ! Loop over rsplit vertically lagrangian timesteps
    call t_startf("prim_step_rX")
    call prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,compute_diagnostics,1)
    call t_stopf("prim_step_rX")

    do r=2,rsplit
       call TimeLevel_update(tl,"leapfrog")
       call t_startf("prim_step_rX")
       call prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,.false.,r)
       call t_stopf("prim_step_rX")
    enddo
    ! defer final timelevel update until after remap and diagnostics

#if (USE_OPENACC)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    call t_startf("copy_qdp_h2d")
    call copy_qdp_d2h( elem , np1_qdp )
    call t_stopf("copy_qdp_h2d")
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  apply vertical remap
    !  always for tracers
    !  if rsplit>0:  also remap dynamics and compute reference level ps_v
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !compute timelevels for tracers (no longer the same as dynamics)
    ! note: time level update for fvm tracers takes place in fvm_mod
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    call vertical_remap(hybrid,elem,fvm,hvcoord,dt_remap,tl%np1,np1_qdp,n0_fvm,nets,nete)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! time step is complete.  update some diagnostic variables:
    ! lnps (we should get rid of this)
    ! Q    (mixing ratio)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call t_startf("prim_run_subcyle_diags")
    do ie=nets,nete
       !dir$ simd
       elem(ie)%state%lnps(:,:,tl%np1)= LOG(elem(ie)%state%ps_v(:,:,tl%np1))
#if (defined COLUMN_OPENMP)
       !$omp parallel do default(shared), private(k,q,dp_np1)
#endif
       do k=1,nlev    !  Loop inversion (AAM)
          !dir$ simd
          dp_np1(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%np1)
          !dir$ simd
          do q=1,qsize
             elem(ie)%state%Q(:,:,k,q)=elem(ie)%state%Qdp(:,:,k,q,np1_qdp)/dp_np1(:,:)
          enddo
       enddo
    enddo
    call t_stopf("prim_run_subcyle_diags")

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - 2*dt
    !   u(n0)    dynamics at  t+dt_remap - dt
    !   u(np1)   dynamics at  t+dt_remap
    !
    !   Q(1)   Q at t+dt_remap
    if (compute_diagnostics) then
      call t_startf("prim_diag_scalars")
      call prim_diag_scalars(elem,hvcoord,tl,2,.false.,nets,nete)
      call t_stopf("prim_diag_scalars")
    endif
    if (compute_energy) then
      call t_startf("prim_energy_halftimes")
      call prim_energy_halftimes(elem,hvcoord,tl,2,.false.,nets,nete)
      call t_stopf("prim_energy_halftimes")
    endif

    if (energy_fixer > 0) then
      call t_startf("prim_energy_fixer")
      call prim_energy_fixer(elem,hvcoord,hybrid,tl,nets,nete,nsubstep)
      call t_stopf("prim_energy_fixer")
    endif

    if (compute_diagnostics) then
      call t_startf("prim_diag_scalars")
      call prim_diag_scalars(elem,hvcoord,tl,3,.false.,nets,nete)
      call t_stopf("prim_diag_scalars")

      call t_startf("prim_energy_halftimes")
      call prim_energy_halftimes(elem,hvcoord,tl,3,.false.,nets,nete)
      call t_stopf("prim_energy_halftimes")
    endif

    ! =================================
    ! update dynamics time level pointers
    ! =================================
    call TimeLevel_update(tl,"leapfrog")
    ! note: time level update for fvm tracers takes place in fvm_mod

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - dt       
    !   u(n0)    dynamics at  t+dt_remap
    !   u(np1)   undefined

    ! ============================================================
    ! Print some diagnostic information
    ! ============================================================
    if (compute_diagnostics) then
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete, fvm)
    end if
  end subroutine prim_run_subcycle



  subroutine prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord, compute_diagnostics,rstep)
  !
  !   Take qsplit dynamics steps and one tracer step
  !   for vertically lagrangian option, this subroutine does only the horizontal step
  !
  !   input:
  !       tl%nm1   not used
  !       tl%n0    data at time t
  !       tl%np1   new values at t+dt_q
  !
  !   then we update timelevel pointers:
  !       tl%nm1 = tl%n0
  !       tl%n0  = tl%np1
  !   so that:
  !       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
  !       tl%n0    time t + dt_q
  !
  !
    use control_mod,        only: statefreq, integration, ftype, qsplit, nu_p, test_cfldep, rsplit
    use control_mod,        only: use_semi_lagrange_transport, tracer_transport_type
    use control_mod,        only: tracer_grid_type, TRACER_GRIDTYPE_GLL
    use derivative_mod,     only: subcell_integration
    use fvm_bsp_mod,        only: get_boomerang_velocities_gll, get_solidbody_velocities_gll
    use fvm_control_volume_mod, only : fvm_supercycling
    use fvm_mod,            only: fvm_ideal_test, IDEAL_TEST_OFF, IDEAL_TEST_ANALYTICAL_WINDS
    use fvm_mod,            only: fvm_test_type, IDEAL_TEST_BOOMERANG, IDEAL_TEST_SOLIDBODY
    use hybvcoord_mod,      only : hvcoord_t
    use parallel_mod,       only: abortmp
    use prim_advance_mod,   only: prim_advance_exp
    use prim_advection_mod, only: prim_advec_tracers_fvm
    use prim_advection_mod, only: prim_advec_tracers_remap, deriv
    use reduction_mod,      only: parallelmax
    use time_mod,           only: time_at,TimeLevel_t, timelevel_update, nsplit

    type(element_t),      intent(inout) :: elem(:)
    type(fvm_struct),     intent(inout) :: fvm(:)
    type(hybrid_t),       intent(in)    :: hybrid   ! distributed parallel structure (shared)
    type(hvcoord_t),      intent(in)    :: hvcoord  ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets     ! starting thread element number (private)
    integer,              intent(in)    :: nete     ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt       ! "timestep dependent" timestep
    type(TimeLevel_t),    intent(inout) :: tl
    integer,              intent(in)    :: rstep    ! vertical remap subcycling step

    real(kind=real_kind) :: st, st1, dp, dt_q
    integer :: ie, t, q,k,i,j,n, n_Q

    real (kind=real_kind)                          :: maxcflx, maxcfly

    real (kind=real_kind) ::  tempdp3d(np,np), x
    real (kind=real_kind) ::  tempmass(nc,nc)
    real (kind=real_kind) ::  tempflux(nc,nc,4)

    real (kind=real_kind) :: dp_np1(np,np)
    logical :: compute_diagnostics

    call t_startf("prim_step_init")
    dt_q = dt*qsplit
    if (ntrac>0.and.rstep==1) then
       !
       ! save velocity at time t for fvm trajectory algorithm
       !       
       do ie=nets,nete
          fvm(ie)%vn0=elem(ie)%state%v(:,:,:,:,tl%n0)
          elem(ie)%sub_elem_mass_flux=0
       end do
    end if
 
    ! ===============
    ! initialize mean flux accumulation variables and save some variables at n0
    ! for use by advection
    ! ===============
    do ie=nets,nete
      elem(ie)%derived%eta_dot_dpdn=0     ! mean vertical mass flux
      elem(ie)%derived%vn0=0              ! mean horizontal mass flux
      elem(ie)%derived%omega_p=0
      if (nu_p>0) then
         elem(ie)%derived%dpdiss_ave=0
         elem(ie)%derived%dpdiss_biharmonic=0
      endif
      ! save velocity at time t for seme-legrangian transport
      !
      ! this code is broken!
      !
      if (fvm_ideal_test == IDEAL_TEST_ANALYTICAL_WINDS) then
         stop
        do k = 1, nlev
          if (fvm_test_type == IDEAL_TEST_BOOMERANG) then
            elem(ie)%derived%vstar(:,:,:,k)=get_boomerang_velocities_gll(elem(ie), time_at(tl%n0))
            stop
          else if (fvm_test_type == IDEAL_TEST_SOLIDBODY) then
            elem(ie)%derived%vstar(:,:,:,k)=get_solidbody_velocities_gll(elem(ie), time_at(tl%n0))
            stop
          else
            call abortmp('Bad fvm_test_type in prim_step')
          end if
        end do
      else if (use_semi_lagrange_transport) then
        elem(ie)%derived%vstar=elem(ie)%state%v(:,:,:,:,tl%n0)
      end if

      elem(ie)%derived%dp(:,:,:)=elem(ie)%state%dp3d(:,:,:,tl%n0)
    enddo
    call t_stopf("prim_step_init")

    ! ===============
    ! Dynamical Step
    ! ===============
    call t_startf("prim_step_dyn")
    n_Q = tl%n0  ! n_Q = timelevel of FV tracers at time t.  need to save this
                 ! FV tracers still carry 3 timelevels
                 ! SE tracers only carry 2 timelevels
    call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
         hybrid, dt, tl, nets, nete, compute_diagnostics)
    do n=2,qsplit
       call TimeLevel_update(tl,"leapfrog")
       call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
            hybrid, dt, tl, nets, nete, .false.)
       ! defer final timelevel update until after Q update.
    enddo
#ifdef HOMME_TEST_SUB_ELEMENT_MASS_FLUX
    if (0<ntrac.and.rstep==1) then
      do ie=nets,nete
      do k=1,nlev
        tempdp3d = elem(ie)%state%dp3d(:,:,k,tl%np1) - &
                   elem(ie)%derived%dp(:,:,k) 
        tempmass = subcell_integration(tempdp3d, np, nc, elem(ie)%metdet)
        tempflux = dt_q*elem(ie)%sub_elem_mass_flux(:,:,:,k)
        do i=1,nc
        do j=1,nc
          x = SUM(tempflux(i,j,:))
          if (ABS(tempmass(i,j)).lt.1e-11 .and. 1e-11.lt.ABS(x)) then
            print *,__FILE__,__LINE__,"**********",ie,k,i,j,tempmass(i,j),x
          elseif (1e-5.lt.ABS((tempmass(i,j)-x)/tempmass(i,j))) then
            print *,__FILE__,__LINE__,"**********",ie,k,i,j,tempmass(i,j),x,&
                   ABS((tempmass(i,j)-x)/tempmass(i,j))
          endif
        end do
        end do
      end do
      end do
    end if
#endif
    call t_stopf("prim_step_dyn")

    ! current dynamics state variables:
    !    derived%dp              =  dp at start of timestep
    !    derived%vstar           =  velocity at start of tracer timestep
    !    derived%vn0             =  mean horiz. flux:   U*dp
    !    state%dp3d(:,:,:,np1)   = dp3d
    ! rsplit=0
    !        state%v(:,:,:,np1)      = velocity on reference levels
    ! rsplit>0
    !        state%v(:,:,:,np1)      = velocity on lagrangian levels 
    !        
    !


    ! ===============
    ! Tracer Advection.  
    ! in addition, this routine will apply the DSS to:
    !        derived%eta_dot_dpdn    =  mean vertical velocity (used for remap below)
    !        derived%omega           =
    ! Tracers are always vertically lagrangian.  
    ! For rsplit=0: 
    !   if tracer scheme needs v on lagrangian levels it has to vertically interpolate
    !   if tracer scheme needs dp3d, it needs to derive it from ps_v
    ! ===============
    ! Advect tracers if their count is > 0.  
    ! special case in CAM: if CSLAM tracers are turned on , qsize=1 but this tracer should 
    ! not be advected.  This will be cleaned up when the physgrid is merged into CAM trunk
    ! Currently advecting all species
    call t_startf("prim_step_advec")
    if (qsize > 0) then
      call t_startf("PAT_remap")
      call Prim_Advec_Tracers_remap(elem, deriv(hybrid%ithr),hvcoord,hybrid,dt_q,tl,nets,nete)
      call t_stopf("PAT_remap")
    end if
    !
    ! only run fvm transport every fvm_supercycling rstep
    !
    if ((ntrac > 0) .and. (mod(rstep,fvm_supercycling) == 0)) then
       !
       ! FVM transport
       !

      if ( n_Q /= tl%n0 ) then
        ! make sure tl%n0 contains tracers at start of timestep
        do ie=nets,nete
          fvm(ie)%c     (:,:,:,1:ntrac,tl%n0)  = fvm(ie)%c     (:,:,:,1:ntrac,n_Q)
          fvm(ie)%dp_fvm(:,:,:,        tl%n0)  = fvm(ie)%dp_fvm(:,:,:,        n_Q)
        end do
      end if
       call t_startf("PAT_fvm")
       call Prim_Advec_Tracers_fvm(elem, fvm, deriv(hybrid%ithr),hvcoord,hybrid,&
            dt_q,tl,nets,nete)
       call t_stopf("PAT_fvm")
       if (rstep.ne.rsplit) then
          !
          ! save velocity for fvm trajecotry algorithm for next fvm time-level update
          !
          do ie=nets,nete
             fvm(ie)%vn0=elem(ie)%state%v(:,:,:,:,tl%np1)
          end do
       end if




       if(test_cfldep) then
         maxcflx=0.0D0
         maxcfly=0.0D0
         do k=1, nlev

!            maxcflx = parallelmax(fvm(:)%maxcfl(1,k),hybrid)
!            maxcfly = parallelmax(fvm(:)%maxcfl(2,k),hybrid)
           maxcflx = max(maxcflx,parallelmax(fvm(:)%maxcfl(1,k),hybrid))
           maxcfly = max(maxcfly,parallelmax(fvm(:)%maxcfl(2,k),hybrid))
          end do

           if(hybrid%masterthread) then
             write(*,*) "nstep",tl%nstep,"dt_fvm=", dt_q*fvm_supercycling, "maximum over all Level"
             write(*,*) "CFL: maxcflx=", maxcflx, "maxcfly=", maxcfly
             print *
           endif
       endif
    endif
    call t_stopf("prim_step_advec")

  end subroutine prim_step


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



!=======================================================================================================!
  subroutine prim_energy_fixer(elem,hvcoord,hybrid,tl,nets,nete,nsubstep)
!
! non-subcycle code:
!  Solution is given at times u(t-1),u(t),u(t+1)
!  E(n=1) = energy before dynamics
!  E(n=2) = energy after dynamics
!
!  fixer will add a constant to the temperature so E(n=2) = E(n=1)
!
    use parallel_mod, only: global_shared_buf, global_shared_sum
    use kinds, only : real_kind
    use hybvcoord_mod, only : hvcoord_t
    use physical_constants, only : Cp
    use time_mod, only : timelevel_t
    use control_mod, only : use_cpstar, energy_fixer
    use hybvcoord_mod, only : hvcoord_t
    use global_norms_mod, only: wrap_repro_sum
    use parallel_mod, only : abortmp
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)
    integer :: t2,n,nets,nete
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (TimeLevel_t), intent(inout)       :: tl
    integer, intent(in)                    :: nsubstep

    integer :: ie,k,i,j,nmax
    real (kind=real_kind), dimension(np,np,nlev)  :: dp   ! delta pressure
    real (kind=real_kind), dimension(np,np,nlev)  :: sumlk
    real (kind=real_kind), pointer  :: PEner(:,:,:)
    real (kind=real_kind), dimension(np,np)  :: suml
    real (kind=real_kind) :: psum(nets:nete,4),psum_g(4),beta

    ! when forcing is applied during dynamics timstep, actual forcing is
    ! slightly different (about 0.1 W/m^2) then expected by the physics
    ! since u & T are changing while FU and FT are held constant.
    ! to correct for this, save compute de_from_forcing at step 1
    ! and then adjust by:  de_from_forcing_step1 - de_from_forcing_stepN
    real (kind=real_kind),save :: de_from_forcing_step1
    real (kind=real_kind)      :: de_from_forcing

    t2=tl%np1    ! timelevel for T
    if (use_cpstar /= 0 ) then
       call abortmp('Energy fixer requires use_cpstar=0')
    endif


    psum = 0
    do ie=nets,nete

#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k,dp)
#endif
       do k=1,nlev
          dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,t2)
       enddo
       suml=0
#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k, i, j)
#endif
       do k=1,nlev
          do i=1,np
          do j=1,np
                sumlk(i,j,k) = cp*dp(i,j,k)
          enddo
          enddo
       enddo
       suml=0
       do k=1,nlev
          do i=1,np
          do j=1,np
             suml(i,j) = suml(i,j) + sumlk(i,j,k)
          enddo
          enddo
       enddo
       PEner => elem(ie)%accum%PEner(:,:,:)

       ! psum(:,4) = energy before forcing
       ! psum(:,1) = energy after forcing, before dynamics
       ! psum(:,2) = energy after dynamics
       ! psum(:,3) = cp*dp (internal energy added is beta*psum(:,3))
       psum(ie,3) = psum(ie,3) + SUM(suml(:,:)*elem(ie)%spheremp(:,:))
       do n=1,2
          psum(ie,n) = psum(ie,n) + SUM(  elem(ie)%spheremp(:,:)*&
               (PEner(:,:,n) + &
               elem(ie)%accum%IEner(:,:,n) + &
               elem(ie)%accum%KEner(:,:,n) ) )
       enddo
    enddo

    nmax=3

    do ie=nets,nete
       do n=1,nmax
          global_shared_buf(ie,n) = psum(ie,n)
       enddo
    enddo
    call wrap_repro_sum(nvars=nmax, comm=hybrid%par%comm)
    do n=1,nmax
       psum_g(n) = global_shared_sum(n)
    enddo

    beta = ( psum_g(1)-psum_g(2) )/psum_g(3)

    ! apply fixer
    do ie=nets,nete
       elem(ie)%state%T(:,:,:,t2) =  elem(ie)%state%T(:,:,:,t2) + beta
    enddo
    end subroutine prim_energy_fixer
!=======================================================================================================!



    subroutine smooth_topo_datasets(phis,sghdyn,sgh30dyn,elem,hybrid,nets,nete)
    use control_mod, only : smooth_phis_numcycle,smooth_sgh_numcycle
    use hybrid_mod, only : hybrid_t
    use bndry_mod, only : bndry_exchangev
    use derivative_mod, only : derivative_t , laplace_sphere_wk
    use viscosity_mod, only : biharmonic_wk
    use prim_advance_mod, only : smooth_phis
    use prim_advection_mod, only: deriv
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
    call smooth_phis(phis,elem,hybrid,deriv(hybrid%ithr),nets,nete,minf,smooth_phis_numcycle)

    minf=0
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to SGH"
    call smooth_phis(sghdyn,elem,hybrid,deriv(hybrid%ithr),nets,nete,minf,smooth_sgh_numcycle)
    if (hybrid%masterthread) &
       write(iulog,*) "Applying hyperviscosity smoother to SGH30"
    call smooth_phis(sgh30dyn,elem,hybrid,deriv(hybrid%ithr),nets,nete,minf,smooth_sgh_numcycle)

    end subroutine smooth_topo_datasets

end module prim_driver_mod



