#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"file: ",__FILE__," line: ",__LINE__," ithr: ",hybrid%ithr
#define _DBG_
module prim_driver_mod
  use kinds, only : real_kind, iulog, longdouble_kind
  use dimensions_mod, only : np, nlev, nlevp, nelem, nelemd, nelemdmax, GlobalUniqueCols, ntrac, qsize, nc,nhc, nep, nipm
  use cg_mod, only : cg_t
  use hybrid_mod, only : hybrid_t
  use quadrature_mod, only : quadrature_t, test_gauss, test_gausslobatto, gausslobatto
#ifndef CAM
  use column_types_mod, only : ColumnModel_t
  use prim_restart_mod, only : initrestartfile
  use restart_io_mod , only : RestFile,readrestart
  use Manager
#endif
  use prim_si_ref_mod, only : ref_state_t
  use solver_mod, only : blkjac_t
  use filter_mod, only : filter_t
  use derivative_mod, only : derivative_t
  use reduction_mod, only : reductionbuffer_ordered_1d_t, red_min, red_max, &
         red_sum, red_sum_int, red_flops, initreductionbuffer

  use fvm_mod, only : fvm_init1,fvm_init2, fvm_init3
  use fvm_control_volume_mod, only : fvm_struct
#if defined(_SPELT)
  use spelt_mod, only : spelt_struct, spelt_init1,spelt_init2, spelt_init3
#endif

  use element_mod, only : element_t, timelevels,  allocate_element_desc
  use thread_mod, only : omp_get_num_threads
  implicit none
  private
  public :: prim_init1, prim_init2 , prim_run, prim_run_subcycle, prim_finalize, leapfrog_bootstrap
  public :: smooth_topo_datasets

  type (cg_t), allocatable  :: cg(:)              ! conjugate gradient struct (nthreads)
  type (quadrature_t)   :: gp                     ! element GLL points
  real(kind=longdouble_kind)  :: fvm_corners(nc+1)     ! fvm cell corners on reference element
  real(kind=longdouble_kind)  :: fvm_points(nc)     ! fvm cell centers on reference element
  real (kind=longdouble_kind) :: spelt_refnep(1:nep)


#ifndef CAM
  type (ColumnModel_t), allocatable :: cm(:) ! (nthreads)
#endif
  type (ref_state_t)    :: refstate        ! semi-implicit vertical reference state
  type (blkjac_t),allocatable  :: blkjac(:)  ! (nets:nete)
  type (filter_t)       :: flt             ! Filter struct for v and p grid
  type (filter_t)       :: flt_advection   ! Filter struct for v grid for advection only
  real*8  :: tot_iter
  type (ReductionBuffer_ordered_1d_t), save :: red   ! reduction buffer               (shared)

contains

  subroutine prim_init1(elem, fvm, par, dom_mt, Tl)

    ! --------------------------------
    use thread_mod, only : nthreads, omp_get_thread_num, omp_set_num_threads, &
                           vert_num_threads
    ! --------------------------------
    use control_mod, only : runtype, restartfreq, filter_counter, integration, topology, &
         partmethod, while_iter
    ! --------------------------------
    use prim_state_mod, only : prim_printstate_init
    ! --------------------------------
    use namelist_mod, only : readnl
    ! --------------------------------
    use mesh_mod, only : MeshUseMeshFile
    ! --------------------------------
    use time_mod, only : nmax, time_at, timelevel_init, timelevel_t
    ! --------------------------------
    ! --------------------------------
    use mass_matrix_mod, only : mass_matrix
    ! --------------------------------
    use cube_mod,  only : cubeedgecount , cubeelemcount, cubetopology
    ! --------------------------------
    use mesh_mod, only : MeshSetCoordinates, MeshUseMeshFile, MeshCubeTopology, &
         MeshCubeElemCount, MeshCubeEdgeCount
    use cube_mod, only : cube_init_atomic, rotation_init_atomic, set_corner_coordinates, assign_node_numbers_to_elem
    ! --------------------------------
    use metagraph_mod, only : metavertex_t, metaedge_t, localelemcount, initmetagraph
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
    use diffusion_mod, only      : diffusion_init
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
#ifndef CAM
    use repro_sum_mod, only: repro_sum, repro_sum_defaultopts, &
         repro_sum_setopts
#else
    use shr_reprosum_mod, only: repro_sum => shr_reprosum_calc
#endif
    implicit none
    type (element_t), pointer :: elem(:)
#if defined(_SPELT)
    type (spelt_struct), pointer   :: fvm(:)
#else
     type (fvm_struct), pointer   :: fvm(:)
#endif
    type (parallel_t), intent(in) :: par
    type (domain1d_t), pointer :: dom_mt(:)
    type (timelevel_t), intent(out) :: Tl
    ! Local Variables

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

    !debug  call PrintGridVertex(GridVertex)


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

    xtmp=nep-1
    do i=1,nep
      spelt_refnep(i)= 2*(i-1)/xtmp - 1
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
          call rotation_init_atomic(elem(ie),rot_type)
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
    filter_counter = 0

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
    call omp_set_num_threads(n_domains)

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
#ifndef CAM
    allocate(cm(0:n_domains-1))
#endif
    allocate(cg(0:n_domains-1))
    call prim_advance_init(par,integration)
    call Prim_Advec_Init1(par, n_domains)
    call diffusion_init(par)
    if (ntrac>0) then
#if defined(_SPELT)
      call spelt_init1(par)
#else
      call fvm_init1(par)
#endif
    endif
    call TimeLevel_init(tl)
    if(par%masterproc) write(iulog,*) 'end of prim_init'
  end subroutine prim_init1
!=======================================================================================================!

  subroutine prim_init2(elem, fvm, hybrid, nets, nete, tl, hvcoord)

    use parallel_mod, only : parallel_t, haltmp, syncmp, abortmp
    use time_mod, only : timelevel_t, tstep, phys_tscale, timelevel_init, nendstep, smooth, nsplit, TimeLevel_Qdp
    use prim_state_mod, only : prim_printstate, prim_diag_scalars
    use filter_mod, only : filter_t, fm_filter_create, taylor_filter_create, &
         fm_transfer, bv_transfer
    use control_mod, only : runtype, integration, filter_mu, filter_mu_advection, test_case, &
         debug_level, vfile_int, filter_freq, filter_freq_advection, &
         transfer_type, vform, vfile_mid, filter_type, kcut_fm, wght_fm, p_bv, &
         s_bv, topology,columnpackage, moisture, precon_method, rsplit, qsplit, rk_stage_user,&
         sub_case, &
         limiter_option, nu, nu_q, nu_div, tstep_type, hypervis_subcycle, &
         hypervis_subcycle_q, use_semi_lagrange_transport
    use control_mod, only : tracer_transport_type
    use control_mod, only : TRACERTRANSPORT_LAGRANGIAN_FVM, TRACERTRANSPORT_FLUXFORM_FVM, TRACERTRANSPORT_SE_GLL
#ifndef CAM
    use control_mod, only : pertlim                     !used for homme temperature perturbations
#endif
    use prim_si_ref_mod, only: prim_si_refstate_init, prim_set_mass
    use bndry_mod, only : sort_neighbor_buffer_mapping
#ifdef TRILINOS
    use prim_derived_type_mod ,only : derived_type, initialize
    use, intrinsic :: iso_c_binding
#endif
    use thread_mod, only : nthreads
    use derivative_mod, only : derivinit, interpolate_gll2fvm_points, interpolate_gll2spelt_points, v2pinit
    use global_norms_mod, only : test_global_integral, print_cfl
    use hybvcoord_mod, only : hvcoord_t
    use prim_advection_mod, only: prim_advec_init2, deriv
#ifdef CAM
#else
    use column_model_mod, only : InitColumnModel
    use held_suarez_mod, only : hs0_init_state
    use baroclinic_inst_mod, only : binst_init_state, jw_baroclinic
    use asp_tests, only : asp_tracer, asp_baroclinic, asp_rossby, asp_mountain, asp_gravity_wave, dcmip2_schar
    use aquaplanet, only : aquaplanet_init_state
    use fvm_control_volume_mod, only : n0_fvm, np1_fvm
#endif
#if USE_CUDA_FORTRAN
    use cuda_mod, only: cuda_mod_init
#endif

    type (element_t), intent(inout) :: elem(:)
#if defined(_SPELT)
    type (spelt_struct), intent(inout)   :: fvm(:)
#else
     type (fvm_struct), intent(inout)    :: fvm(:)
#endif
    type (hybrid_t), intent(in) :: hybrid

    type (TimeLevel_t), intent(inout)    :: tl              ! time level struct
    type (hvcoord_t), intent(inout)      :: hvcoord         ! hybrid vertical coordinate struct

     integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)


    ! ==================================
    ! Local variables
    ! ==================================

    real (kind=real_kind) :: dt              ! "timestep dependent" timestep
!   variables used to calculate CFL
    real (kind=real_kind) :: dtnu            ! timestep*viscosity parameter
    real (kind=real_kind) :: dt_dyn_vis      ! viscosity timestep used in dynamics
    real (kind=real_kind) :: dt_tracer_vis      ! viscosity timestep used in tracers

    real (kind=real_kind) :: dp


    real (kind=real_kind) :: ps(np,np)       ! surface pressure

    character(len=80)     :: fname
    character(len=8)      :: njusn
    character(len=4)      :: charnum

    real (kind=real_kind) :: Tp(np)     ! transfer function

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
    type(derived_type) ,pointer         :: pptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_pre

    type(derived_type) ,target         :: jac_object
    type(derived_type) ,pointer         :: jptr=>NULL()
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

    ! ==========================
    ! begin executable code
    ! ==========================
    if (topology == "cube") then
       call test_global_integral(elem, hybrid,nets,nete)
    end if


    ! compute most restrictive dt*nu for use by variable res viscosity:
    if (tstep_type == 0) then
       ! LF case: no tracers, timestep seen by viscosity is 2*tstep
       dt_tracer_vis = 0
       dt_dyn_vis = 2*tstep
       dtnu = 2.0d0*tstep*max(nu,nu_div)
    else
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
    endif

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
    call Prim_Advec_Init2(hybrid, fvm_corners, fvm_points, spelt_refnep)

    ! ================================================
    ! fvm initialization
    ! ================================================
    if (ntrac>0) then
#if defined(_SPELT)
      call spelt_init2(elem,fvm,hybrid,nets,nete,tl)
#else
      call fvm_init2(elem,fvm,hybrid,nets,nete,tl)
#endif
    endif
    ! ====================================
    ! In the semi-implicit case:
    ! initialize vertical structure and
    ! related matrices..
    ! ====================================
#if (defined HORIZ_OPENMP)
!$OMP MASTER
#endif
    if (integration == "semi_imp") then
       refstate = prim_si_refstate_init(.false.,hybrid%masterthread,hvcoord)
    endif
#if (defined HORIZ_OPENMP)
!$OMP END MASTER
#endif
    ! ==========================================
    ! Initialize pressure and velocity grid
    ! filter matrix...
    ! ==========================================
    if (transfer_type == "bv") then
       Tp    = bv_transfer(p_bv,s_bv,np)
    else if (transfer_type == "fm") then
       Tp    = fm_transfer(kcut_fm,wght_fm,np)
    end if
    if (filter_type == "taylor") then
       flt           = taylor_filter_create(Tp, filter_mu,gp)
       flt_advection = taylor_filter_create(Tp, filter_mu_advection,gp)
    else if (filter_type == "fischer") then
       flt           = fm_filter_create(Tp, filter_mu, gp)
       flt_advection = fm_filter_create(Tp, filter_mu_advection, gp)
    end if



    if (hybrid%masterthread) then
       if (filter_freq>0 .or. filter_freq_advection>0) then
          write(iulog,*) "transfer function type in preq=",transfer_type
          write(iulog,*) "filter type            in preq=",filter_type
          write(*,'(a,99f10.6)') "dynamics: I-mu + mu*Tp(:) = ",&
               (1-filter_mu)+filter_mu*Tp(:)
          write(*,'(a,99f10.6)') "advection: I-mu + mu*Tp(:) = ",&
               (1-filter_mu_advection)+filter_mu_advection*Tp(:)
       endif
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

       if (test_case(1:10) == "aquaplanet") then
          if (hybrid%masterthread) then
             write(iulog,*)  'Initializing aqua planet with MJO-type forcing'
          end if
          if(moisture.eq."dry") then
             call binst_init_state(elem, hybrid,nets,nete,hvcoord)
          end if
          call aquaplanet_init_state(elem, hybrid,hvcoord,nets,nete,integration)
       end if

       call ReadRestart(elem,hybrid%ithr,nets,nete,tl)

       ! scale PS to achieve prescribed dry mass
       if (runtype /= 1) &
            call prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)

       if (runtype==2) then
          ! copy prognostic variables:  tl%n0 into tl%nm1
          do ie=nets,nete
             elem(ie)%state%v(:,:,:,:,tl%nm1)=elem(ie)%state%v(:,:,:,:,tl%n0)
             elem(ie)%state%T(:,:,:,tl%nm1)=elem(ie)%state%T(:,:,:,tl%n0) 
             elem(ie)%state%ps_v(:,:,tl%nm1)=elem(ie)%state%ps_v(:,:,tl%n0)
             elem(ie)%state%lnps(:,:,tl%nm1)=elem(ie)%state%lnps(:,:,tl%n0)
          enddo
       endif ! runtype==2
    else  ! initial run  RUNTYPE=0
       ! ===========================================================
       ! Initial Run  - compute initial condition
       ! ===========================================================
       if (hybrid%masterthread) then
          write(iulog,*) ' runtype: INITIAL primitive equations'
       endif
       ! ========================================================
       ! Initialize the test cases
       ! ========================================================
       if (test_case(1:10) == "baroclinic") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing Polvani-Scott-Thomas baroclinic instability test'
          end if

          call binst_init_state(elem, hybrid,nets,nete,hvcoord)
       else if (test_case(1:16) == "asp_gravity_wave") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing ASP gravity wave test'
          end if
          call asp_gravity_wave(elem, hybrid,hvcoord,nets,nete, sub_case)
       else if (test_case(1:12) == "asp_mountain") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing ASP mountain Rossby test'
          end if
          call asp_mountain(elem, hybrid,hvcoord,nets,nete)
       else if (test_case(1:10) == "asp_rossby") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing ASP Rossby Haurwitz test'
          end if
          call asp_rossby(elem, hybrid,hvcoord,nets,nete)
       else if (test_case(1:10) == "asp_tracer") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing pure tracer advection tests'
          end if
          call asp_tracer(elem, hybrid,hvcoord,nets,nete)
       else if (test_case(1:14) == "asp_baroclinic") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing Jablonowski and Williamson ASP baroclinic instability test'
          end if
          call asp_baroclinic(elem, hybrid,hvcoord,nets,nete,fvm)
       else if (test_case(1:13) == "jw_baroclinic") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing Jablonowski and Williamson baroclinic instability test V1'
          end if
          call jw_baroclinic(elem, hybrid,hvcoord,nets,nete)
       else if (test_case(1:12) == "held_suarez0") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing Held-Suarez primitive equations test'
          end if
          call hs0_init_state(elem, hvcoord,nets,nete,300.0_real_kind)
       else if (test_case(1:10) == "aquaplanet") then
          if (hybrid%masterthread) then
             write(iulog,*)  'Initializing aqua planet with MJO-type forcing'
          end if
          if(moisture.eq."dry") then
             call binst_init_state(elem, hybrid,nets,nete,hvcoord)
          end if
          call aquaplanet_init_state(elem, hybrid,hvcoord,nets,nete,integration)
       else if (test_case(1:12) == "dcmip2_schar") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing DCMIP2 test 2-0'
          end if
          call dcmip2_schar(elem, hybrid,hvcoord,nets,nete)
       else
          call abortmp('Error: unrecognized test case')
       endif

       if (hybrid%masterthread) then
          write(iulog,*) '...done'
       end if

       ! scale PS to achieve prescribed dry mass
       call prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)

       do ie=nets,nete

          elem(ie)%state%T=elem(ie)%state%T &
                * (1.0_real_kind + pertlim)  ! set perlim in ctl_nl namelist for
                                             ! temperature field initial
                                             ! perterbation
       enddo
 
       ! ========================================
       ! Print state and movie output
       ! ========================================
    end if  ! runtype

!$OMP MASTER
    tl%nstep0=2   ! This will be the first full leapfrog step
    if (runtype==1) then
       tl%nstep0=tl%nstep+1            ! restart run: first step = first first full leapfrog step
    endif
    if (runtype==2) then
       ! branch run
       ! reset time counters to zero since timestep may have changed
       nEndStep = nEndStep-tl%nstep ! restart set this to nmax + tl%nstep
       tl%nstep=0
    endif
!$OMP END MASTER
!$OMP BARRIER
#endif

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
#if defined(_SPELT)
      ! do it only for SPELT tracers, FIRST TRACER will be the AIR DENSITY
      ! should be optimize and combined with the above caculation
      do ie=nets,nete
        do k=1,nlev
	    do i=1,np
	      do j=1,np
		  elem(ie)%derived%dp(i,j,k)=( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
		       ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,tl%n0)
	      enddo
	    enddo
          !write air density in tracer 1 of FVM 
          fvm(ie)%c(1:nep,1:nep,k,1,tl%n0)=interpolate_gll2spelt_points(elem(ie)%derived%dp(:,:,k),deriv(hybrid%ithr))
        enddo
      enddo
      call spelt_init3(elem,fvm,hybrid,nets,nete,tl%n0)
      do ie=nets,nete
	    do i=1-nipm,nep+nipm
	      do j=1-nipm,nep+nipm
	        fvm(ie)%psc(i,j) = sum(fvm(ie)%c(i,j,:,1,tl%n0) +  hvcoord%hyai(1)*hvcoord%ps0)
	      enddo
	    enddo
      enddo

      if (hybrid%masterthread) then
         write(iulog,*) 'FVM (Spelt) tracers (incl. in halo zone) initialized. FIRST tracer has air density!'
      end if
#else
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
!phl      call fvm_init3(elem,fvm,hybrid,nets,nete,n0_fvm) !boundary exchange now takes place in transport subroutine (cslam_runflux)
      do ie=nets,nete
!	    do i=1-nhc,nc+nhc
!	      do j=1-nhc,nc+nhc
!phl is it necessary to compute psc here?
	    do i=1,nc
	      do j=1,nc
	        fvm(ie)%psc(i,j) = sum(fvm(ie)%dp_fvm(i,j,:,n0_fvm)) +  hvcoord%hyai(1)*hvcoord%ps0
	      enddo
	    enddo
      enddo
      if (hybrid%masterthread) then
         write(iulog,*) 'FVM tracers initialized.'
      end if
#endif
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

       if (tracer_transport_type == TRACERTRANSPORT_FLUXFORM_FVM.or.&
           tracer_transport_type == TRACERTRANSPORT_LAGRANGIAN_FVM) then
          write(iulog,'(a,2f9.2)') "dt_tracer (fvm)          ",tstep*qsplit
       else if (tracer_transport_type == TRACERTRANSPORT_SE_GLL) then
          write(iulog,'(a,2f9.2)') "dt_tracer, per RK stage: ",tstep*qsplit,(tstep*qsplit)/(rk_stage_user-1)
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


#if USE_CUDA_FORTRAN
    !Inside this routine, we enforce an OMP BARRIER and an OMP MASTER. It's left out of here because it's ugly
    call cuda_mod_init(elem,hybrid,deriv(hybrid%ithr),hvcoord)
#endif
    if (hybrid%masterthread) write(iulog,*) "initial state:"
    call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete, fvm)

    if ( use_semi_lagrange_transport) then
      call sort_neighbor_buffer_mapping(hybrid, elem,nets,nete)
    end if
  end subroutine prim_init2

!=======================================================================================================!



  subroutine leapfrog_bootstrap(elem, hybrid,nets,nete,tstep,tl,hvcoord)

  !
  ! leapfrog bootstrap code.
  !
  ! take the equivilent of one timestep, but do it with a
  ! dt/2 euler and a dt/2 leapfrog step
  !
  use hybvcoord_mod, only : hvcoord_t
  use time_mod, only : TimeLevel_t

  type (element_t) , intent(inout)        :: elem(:)
  type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)
  type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct
  integer, intent(in)                     :: nets  ! starting thread element number (private)
  integer, intent(in)                     :: nete  ! ending thread element number   (private)
  real(kind=real_kind), intent(in)        :: tstep          ! "timestep dependent" timestep
  type (TimeLevel_t), intent(inout)       :: tl


  ! local
  real(kind=real_kind) :: tstep_tmp,tstep_dyn
  integer :: i,ie


  tstep_dyn = tstep
  ! forward euler to get to tstep_dyn/2 (keep t=0 in nm1 timelevel)
  ! (note: leapfrog tstep_dyn/4 with nm1=n0 is Euler with tstep_dyn/2 )
  tstep_tmp=tstep_dyn/4

  call prim_run(elem, hybrid,nets,nete, tstep_tmp, tl, hvcoord, "forward")

  ! leapfrog with tstep_dyn/2 to get to tstep_dyn (keep t=0 in nm1 timelevel)
  tstep_tmp=tstep_dyn/2
  call prim_run(elem, hybrid,nets,nete, tstep_tmp, tl, hvcoord, "forward")


  tl%nstep=tl%nstep-1        ! count all of that as 1 timestep

  end subroutine leapfrog_bootstrap


!=======================================================================================================!


  subroutine prim_run(elem, hybrid,nets,nete, dt, tl, hvcoord, advance_name)
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t, timelevel_update, smooth
    use control_mod, only: statefreq, integration, ftype, qsplit, disable_diagnostics
    use prim_advance_mod, only : prim_advance_exp, prim_advance_si, preq_robert3
    use prim_state_mod, only : prim_printstate, prim_diag_scalars, prim_energy_halftimes
    use prim_advection_mod, only: deriv
    use parallel_mod, only : abortmp
#ifndef CAM
    use column_model_mod, only : ApplyColumnModel
#endif

    type (element_t) , intent(inout)        :: elem(:)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    real(kind=real_kind), intent(in)        :: dt              ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout)       :: tl
    character(len=*), intent(in) :: advance_name
    real(kind=real_kind) :: st, st1, dp
    integer :: ie, t, q,k,i,j


    logical :: compute_diagnostics

    ! ===================================
    ! Main timestepping loop
    ! ===================================

    ! compute diagnostics and energy for STDOUT
    ! compute energy if we are using an energy fixer

    if (MODULO(tl%nstep+1,statefreq)==0 .or. tl%nstep+1==tl%nstep0) then
       compute_diagnostics=.true.
    else
       compute_diagnostics=.false.
    endif

    if(disable_diagnostics) compute_diagnostics=.false.

    tot_iter=0.0


    ! Forcing options for testing CAM-HOMME energy balance:
    if (ftype == -1) then
       ! disable all forcing, but allow moisture:
       do ie=nets,nete
          elem(ie)%derived%FQ = 0
          elem(ie)%derived%FM = 0
          elem(ie)%derived%FT = 0
       enddo
    endif
    if (ftype == -2) then
       ! disable moisture, but allow dynamics forcing
       do ie=nets,nete
          elem(ie)%state%Q = 0
          elem(ie)%state%Qdp = 0
          elem(ie)%derived%FQ = 0
       enddo
    endif
    if (ftype == -3) then
       ! disable forcing & moisture
       do ie=nets,nete
          elem(ie)%state%Q = 0
          elem(ie)%state%Qdp = 0
          elem(ie)%derived%FQ = 0
          elem(ie)%derived%FM = 0
          elem(ie)%derived%FT = 0
       enddo
    endif

    ! =================================
    ! energy, dissipation rate diagnostics.  Uses data at t-1,t
    ! to compute diagnostics at t - 0.5.
    ! small error in the t+.5 terms because at this
    ! point only state variables at t-1 has been Robert filtered.
    ! =================================
    if (compute_diagnostics) then
       call prim_energy_halftimes(elem,hvcoord,tl,1,.true.,nets,nete)
       call prim_diag_scalars(elem,hvcoord,tl,1,.true.,nets,nete)
    endif

    ! ===============
    ! initialize mean flux accumulation variables
    ! ===============
    do ie=nets,nete
       elem(ie)%derived%eta_dot_dpdn=0
       elem(ie)%derived%omega_p=0
    enddo

    ! ===============
    ! Dynamical Step  uses Q at tl%n0
    ! ===============
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
    if (integration == "semi_imp") then
       call prim_advance_si(elem, nets, nete, cg(hybrid%ithr), blkjac, red, &
            refstate, hvcoord, deriv(hybrid%ithr), flt, hybrid, tl, dt)
       tot_iter=tot_iter+cg(hybrid%ithr)%iter
    else if (integration == "full_imp") then
       call abortmp('full_imp integration requires tstep_type > 0')
    else
       call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
            hybrid, dt, tl, nets, nete, compute_diagnostics)

       ! keep lnps up to date (we should get rid of this requirement)
       do ie=nets,nete
          elem(ie)%state%lnps(:,:,tl%np1)= LOG(elem(ie)%state%ps_v(:,:,tl%np1))
       enddo
    end if

    ! =================================
    ! energy, dissipation rate diagnostics.  Uses data at t and t+1
    ! to compute diagnostics at t + 0.5.
    ! =================================
    if (compute_diagnostics) then
       call prim_energy_halftimes(elem,hvcoord,tl,2,.false.,nets,nete)
       call prim_diag_scalars(elem,hvcoord,tl,2,.false.,nets,nete)
    endif

    ! ===================================
    ! Compute Forcing Tendencies from nm1 data (for PROCESS SPLIT)
    ! or np1 data (for TIMESPLIT) and add tendencies into soluiton at timelevel np1
    ! ===================================
#ifdef CAM
    call abortmp('CAM-HOMME-SE requires RK timestepping option turned on')
#else
    call ApplyColumnModel(elem, hybrid, hvcoord, cm(hybrid%ithr),dt)
#endif
    ! measure the effects of forcing
    if (compute_diagnostics) then
       call prim_energy_halftimes(elem,hvcoord,tl,3,.false.,nets,nete)
       call prim_diag_scalars(elem,hvcoord,tl,3,.false.,nets,nete)
    endif


    ! =================================
    ! timestep is complete.
    ! =================================


    !Now apply robert filter to all prognostic variables
    if (smooth/=0) &
       call preq_robert3(tl%nm1,tl%n0,tl%np1,elem,hvcoord,nets,nete)
    ! measure the effects of Robert filter
    if (compute_diagnostics) then
       call prim_energy_halftimes(elem,hvcoord,tl,4,.false.,nets,nete)
       call prim_diag_scalars(elem,hvcoord,tl,4,.false.,nets,nete)
    endif


    ! =================================
    ! update dynamics time level pointers
    ! =================================
    call TimeLevel_update(tl,advance_name)

  ! ============================================================
    ! Print some diagnostic information
    ! ============================================================

    if (compute_diagnostics) then
       if (hybrid%masterthread) then
          if (integration == "semi_imp") write(iulog,*) "cg its=",cg(0)%iter
       end if
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)
    end if
  end subroutine prim_run

!=======================================================================================================!


  subroutine prim_run_subcycle(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,nsubstep)
!
!   advance all variables (u,v,T,ps,Q,C) from time t to t + dt_q
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
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t, timelevel_update, timelevel_qdp, nsplit
    use control_mod, only: statefreq,&
           energy_fixer, ftype, qsplit, rsplit, test_cfldep, disable_diagnostics
    use prim_advance_mod, only : applycamforcing, &
                                 applycamforcing_dynamics
    use prim_state_mod, only : prim_printstate, prim_diag_scalars, prim_energy_halftimes
    use parallel_mod, only : abortmp
    use reduction_mod, only : parallelmax
    use prim_advection_mod, only : vertical_remap
    use fvm_control_volume_mod, only : n0_fvm
#if USE_CUDA_FORTRAN
    use cuda_mod, only: copy_qdp_h2d, copy_qdp_d2h
#endif


    type (element_t) , intent(inout)        :: elem(:)

#if defined(_SPELT)
      type(spelt_struct), intent(inout) :: fvm(:)
#else
      type(fvm_struct), intent(inout) :: fvm(:)
#endif
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    real(kind=real_kind), intent(in)        :: dt  ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout)       :: tl
    integer, intent(in)                     :: nsubstep  ! nsubstep = 1 .. nsplit
    real(kind=real_kind) :: st, st1, dp, dt_q, dt_remap
    integer :: ie, t, q,k,i,j,n, n_Q
    integer :: n0_qdp,np1_qdp,r, nstep_end

    real (kind=real_kind)                          :: maxcflx, maxcfly
    real (kind=real_kind) :: dp_np1(np,np)
    logical :: compute_diagnostics, compute_energy


    ! ===================================
    ! Main timestepping loop
    ! ===================================
    dt_q = dt*qsplit
    dt_remap = dt_q
    nstep_end = tl%nstep + qsplit
    if (rsplit>0) then
       dt_remap=dt_q*rsplit   ! rsplit=0 means use eulerian code, not vert. lagrange
       nstep_end = tl%nstep + qsplit*rsplit  ! nstep at end of this routine
    endif



    ! compute diagnostics and energy for STDOUT
    ! compute energy if we are using an energy fixer
    compute_diagnostics=.false.
    compute_energy=energy_fixer > 0

    if (MODULO(nstep_end,statefreq)==0 .or. nstep_end==tl%nstep0) then
       compute_diagnostics=.true.
       compute_energy = .true.
    endif

    if(disable_diagnostics) compute_diagnostics=.false.

    if (compute_diagnostics) &
       call prim_diag_scalars(elem,hvcoord,tl,4,.true.,nets,nete)



#ifdef CAM
    ! ftype=2  Q was adjusted by physics, but apply u,T forcing here
    ! ftype=1  forcing was applied time-split in CAM coupling layer
    ! ftype=0 means forcing apply here
    ! ftype=-1 do not apply forcing
    call TimeLevel_Qdp( tl, qsplit, n0_qdp)
    if (ftype==0) then
      call ApplyCAMForcing(elem, fvm, hvcoord,tl%n0,n0_qdp, dt_remap,nets,nete)
    end if
    if (ftype==2) call ApplyCAMForcing_dynamics(elem, hvcoord,tl%n0,dt_remap,nets,nete)
#endif

    ! E(1) Energy after CAM forcing
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,1,.true.,nets,nete)

    ! qmass and variance, using Q(n0),Qdp(n0)
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,1,.true.,nets,nete)

    ! initialize dp3d from ps
    if (rsplit>0) then
    do ie=nets,nete
       do k=1,nlev
          elem(ie)%state%dp3d(:,:,k,tl%n0)=&
               ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
       enddo
       ! DEBUGDP step: ps_v should not be used for rsplit>0 code during prim_step
       ! vertical_remap.  so to this for debugging:
       elem(ie)%state%ps_v(:,:,tl%n0)=-9e9
    enddo
    endif

#if USE_CUDA_FORTRAN
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    call copy_qdp_h2d( elem , n0_qdp )
#endif

    ! loop over rsplit vertically lagrangian timesteps
    call prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,compute_diagnostics)
    do r=2,rsplit
       call TimeLevel_update(tl,"leapfrog")
       call prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord,.false.)
    enddo
    ! defer final timelevel update until after remap and diagnostics


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  apply vertical remap
    !  always for tracers
    !  if rsplit>0:  also remap dynamics and compute reference level ps_v
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !compute timelevels for tracers (no longer the same as dynamics)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    ! note: time level update for fvm tracers takes place in fvm_mod
    call vertical_remap(hybrid,elem,fvm,hvcoord,dt_remap,tl%np1,np1_qdp,n0_fvm,nets,nete)

#if USE_CUDA_FORTRAN
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    call copy_qdp_d2h( elem , np1_qdp )
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! time step is complete.  update some diagnostic variables:
    ! lnps (we should get rid of this)
    ! Q    (mixing ratio)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       elem(ie)%state%lnps(:,:,tl%np1)= LOG(elem(ie)%state%ps_v(:,:,tl%np1))
#if (defined COLUMN_OPENMP)
       !$omp parallel do default(shared), private(k,q,dp_np1)
#endif
       do k=1,nlev    !  Loop inversion (AAM)
          !if (k == 1) then
           !write(*,*) "In prim run there are ", omp_get_num_threads(), " in the current team in parallel region"
          !endif
          dp_np1(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%np1)
          do q=1,qsize
             elem(ie)%state%Q(:,:,k,q)=elem(ie)%state%Qdp(:,:,k,q,np1_qdp)/dp_np1(:,:)
          enddo
       enddo
    enddo





    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - 2*dt
    !   u(n0)    dynamics at  t+dt_remap - dt
    !   u(np1)   dynamics at  t+dt_remap
    !
    !   Q(1)   Q at t+dt_remap
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,2,.false.,nets,nete)
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,2,.false.,nets,nete)

    if (energy_fixer > 0) then
       call prim_energy_fixer(elem,hvcoord,hybrid,tl,nets,nete,nsubstep)
    endif

    if (compute_diagnostics) then
       call prim_diag_scalars(elem,hvcoord,tl,3,.false.,nets,nete)
       call prim_energy_halftimes(elem,hvcoord,tl,3,.false.,nets,nete)
     endif

    ! =================================
    ! update dynamics time level pointers
    ! =================================
    call TimeLevel_update(tl,"leapfrog")
    ! note: time level update for fvm tracers takes place in fvm_mod

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - dt       (Robert-filtered)
    !   u(n0)    dynamics at  t+dt_remap
    !   u(np1)   undefined


    ! ============================================================
    ! Print some diagnostic information
    ! ============================================================
    if (compute_diagnostics) then
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete, fvm)
    end if
  end subroutine prim_run_subcycle






  subroutine prim_step(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord, compute_diagnostics)
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
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t, timelevel_update, nsplit
    use control_mod, only: statefreq, integration, ftype, qsplit, nu_p, test_cfldep, rsplit
    use control_mod, only : use_semi_lagrange_transport, tracer_transport_type
    use control_mod, only : tracer_grid_type, TRACER_GRIDTYPE_GLL
    use fvm_mod,     only : fvm_ideal_test, IDEAL_TEST_OFF, IDEAL_TEST_ANALYTICAL_WINDS
    use fvm_mod,     only : fvm_test_type, IDEAL_TEST_BOOMERANG, IDEAL_TEST_SOLIDBODY
    use fvm_bsp_mod, only : get_boomerang_velocities_gll, get_solidbody_velocities_gll
    use prim_advance_mod, only : prim_advance_exp, overwrite_SEdensity
    use prim_advection_mod, only : prim_advec_tracers_remap, prim_advec_tracers_fvm, deriv
#if defined(_SPELT)
    use prim_advection_mod, only : prim_advec_tracers_spelt
#endif
    use parallel_mod, only : abortmp
    use reduction_mod, only : parallelmax
    use derivative_mod, only : interpolate_gll2spelt_points
    use time_mod,    only : time_at

    type (element_t) , intent(inout)        :: elem(:)

#if defined(_SPELT)
      type(spelt_struct), intent(inout) :: fvm(:)
#else
      type(fvm_struct), intent(inout) :: fvm(:)
#endif
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    real(kind=real_kind), intent(in)        :: dt  ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout)       :: tl
    real(kind=real_kind) :: st, st1, dp, dt_q
    integer :: ie, t, q,k,i,j,n

    real (kind=real_kind)                          :: maxcflx, maxcfly

    real (kind=real_kind) :: dp_np1(np,np)
    logical :: compute_diagnostics

    dt_q = dt*qsplit

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
      if (ntrac>0) then
        ! save velocity at time t for fvm
        fvm(ie)%vn0=elem(ie)%state%v(:,:,:,:,tl%n0)
      end if

      ! save velocity at time t for seme-legrangian transport
      if (fvm_ideal_test == IDEAL_TEST_ANALYTICAL_WINDS) then
        do k = 1, nlev
          if (fvm_test_type == IDEAL_TEST_BOOMERANG) then
            elem(ie)%derived%vstar(:,:,:,k)=get_boomerang_velocities_gll(elem(ie), time_at(tl%n0))
          else if (fvm_test_type == IDEAL_TEST_SOLIDBODY) then
            elem(ie)%derived%vstar(:,:,:,k)=get_solidbody_velocities_gll(elem(ie), time_at(tl%n0))
          else
            call abortmp('Bad fvm_test_type in prim_step')
          end if
        end do
      else if (use_semi_lagrange_transport) then
        elem(ie)%derived%vstar=elem(ie)%state%v(:,:,:,:,tl%n0)
      end if

      if (rsplit==0) then
        ! save dp at time t for use in tracers
#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k)
#endif
         do k=1,nlev
            elem(ie)%derived%dp(:,:,k)=&
                 ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                 ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
         enddo
      else
         ! dp at time t:  use floating lagrangian levels:
         elem(ie)%derived%dp(:,:,:)=elem(ie)%state%dp3d(:,:,:,tl%n0)
      endif
    enddo

    ! ===============
    ! Dynamical Step
    ! ===============
    call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
         hybrid, dt, tl, nets, nete, compute_diagnostics)
    do n=2,qsplit
       call TimeLevel_update(tl,"leapfrog")
       call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
            hybrid, dt, tl, nets, nete, .false.)
       ! defer final timelevel update until after Q update.
    enddo
    ! current dynamics state variables:
    !    derived%dp              =  dp at start of timestep
    !    derived%vstar           =  velocity at start of tracer timestep
    !    derived%vn0             =  mean horiz. flux:   U*dp
    ! rsplit=0
    !        state%v(:,:,:,np1)      = velocity on reference levels
    !        state%ps_v(:,:,:,np1)   = ps
    ! rsplit>0
    !        state%v(:,:,:,np1)      = velocity on lagrangian levels 
    !        state%dp3d(:,:,:,np1)   = dp3d
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
    if (tracer_grid_type == TRACER_GRIDTYPE_GLL) then
      call Prim_Advec_Tracers_remap(elem, deriv(hybrid%ithr),hvcoord,flt_advection,hybrid,&
           dt_q,tl,nets,nete)
    else
       !
       ! FVM transport
       !
#if defined(_SPELT)
       !
      call Prim_Advec_Tracers_spelt(elem, fvm, deriv(hybrid%ithr),hvcoord,hybrid,&
           dt_q,tl,nets,nete)
        do ie=nets,nete
!           do k=1, nlev
!             fvm(ie)%c(1:nep,1:nep,k,1,tl%np1)=interpolate_gll2spelt_points(elem(ie)%derived%dp(:,:,k),deriv(hybrid%ithr))
!           end do
	    do i=1-nipm,nep+nipm
	      do j=1-nipm,nep+nipm
	        fvm(ie)%psc(i,j) = sum(fvm(ie)%c(i,j,:,1,tl%np1)) +  hvcoord%hyai(1)*hvcoord%ps0
	      enddo
	    enddo
        enddo
        if (test_cfldep) then
          maxcflx=0.0D0
          maxcfly=0.0D0
          do k=1, nlev
            maxcflx = max(maxcflx,parallelmax(fvm(:)%maxcfl(1,k),hybrid))
            maxcfly = max(maxcfly,parallelmax(fvm(:)%maxcfl(2,k),hybrid))
          end do

          if  (hybrid%masterthread) then
            write(*,*) "nstep",tl%nstep,"dt_q=", dt_q, "maximum over all Level"
            write(*,*) "CFL: maxcflx=", maxcflx, "maxcfly=", maxcfly
            print *
         endif
       endif
#else
      call Prim_Advec_Tracers_fvm(elem, fvm, deriv(hybrid%ithr),hvcoord,hybrid,&
           dt_q,tl,nets,nete)

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
             write(*,*) "nstep",tl%nstep,"dt_q=", dt_q, "maximum over all Level"
             write(*,*) "CFL: maxcflx=", maxcflx, "maxcfly=", maxcfly
             print *
           endif
       endif
       !overwrite SE density by fvm(ie)%psc
!        call overwrite_SEdensity(elem,fvm,dt_q,hybrid,nets,nete,tl%np1)
#endif
    endif

  end subroutine prim_step


!=======================================================================================================!


  subroutine prim_finalize(hybrid)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

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
!$omp parallel do default(shared), private(k)
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
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
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



