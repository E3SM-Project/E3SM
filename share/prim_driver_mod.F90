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
  use spelt_mod, only : spelt_struct, spelt_init1,spelt_init2, spelt_init3
  
  use element_mod, only : element_t, timelevels
  use time_mod, only           : TimeLevel_Qdp

  implicit none
  private
  public :: prim_init1, prim_init2 , prim_run, prim_run_subcycle, prim_finalize, leapfrog_bootstrap
  public :: smooth_topo_datasets

  type (cg_t), allocatable  :: cg(:)              ! conjugate gradient struct (nthreads)
  type (quadrature_t)   :: gp                     ! element GLL points
  real(kind=longdouble_kind)  :: fvm_corners(nc+1)     ! fvm cell corners on reference element
  real(kind=longdouble_kind)  :: fvm_points(nc)     ! fvm cell centers on reference element


#ifndef CAM
  type (ColumnModel_t), allocatable :: cm(:) ! (nthreads)
#endif
  type (ref_state_t)    :: refstate        ! semi-implicit vertical reference state
  type (blkjac_t),allocatable  :: blkjac(:)  ! (nets:nete)
  type (filter_t)       :: flt             ! Filter struct for v and p grid
  type (filter_t)       :: flt_advection   ! Filter struct for v grid for advection only
  type (derivative_t), allocatable   :: deriv(:) ! derivative struct (nthreads)
  real*8  :: tot_iter
  type (ReductionBuffer_ordered_1d_t), save :: red   ! reduction buffer               (shared)

contains

  subroutine prim_init1(elem, fvm, par, dom_mt, Tl)

    ! --------------------------------
    use thread_mod, only : nthreads, omp_get_thread_num, omp_set_num_threads
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
    use gridgraph_mod, only : gridvertex_t, gridedge_t
    ! --------------------------------
    use schedule_mod, only : schedule, genEdgeSched,  PrintSchedule
    ! --------------------------------
    use prim_advection_mod, only: prim_advec_init
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
#ifdef TESTGRID
    use checksum_mod, only : testchecksum
#endif
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
    integer :: err, ierr, l

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
#ifndef CAM
       call ManagerInit()
#endif
    endif

    if (ntrac>0) allocate(fvm(nelemd))

    ! ====================================================
    !  Generate the communication schedule
    ! ====================================================

    call genEdgeSched(elem,iam,Schedule(1),MetaVertex(1))


    allocate(global_shared_buf(nelemd,nrepro_vars))
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
    n_domains = Nthreads
#if (defined ELEMENT_OPENMP)
    n_domains = 1
#endif


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
    ! Run the checksum to verify communication schedule
    ! =================================================================
#ifdef TESTGRID 
    if(par%masterproc)     write(iulog,*) 'running testchecksum ',iam,nelem,nelemd
    call testchecksum(par,GridEdge)
#endif

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
       elem(ie)%accum%mass_added=0
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

#ifndef TESTGRID
    deallocate(GridEdge)
    deallocate(GridVertex)
    deallocate(MetaVertex)
    deallocate(TailPartition)
    deallocate(HeadPartition)
#else
    ! here we need to call a function in gridgraph_mod.F to deallocate
    ! all of GridVertex's EdgeIndex pointers
#endif


    ! =====================================
    ! Set number of threads...
    ! =====================================
    if(par%masterproc) then
       write(iulog,*) "Main:NThreads=",NThreads
       write(iulog,*) "Main:n_domains = ",n_domains
    endif
    call omp_set_num_threads(NThreads)

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
    allocate(deriv(0:n_domains-1))
    allocate(cg(0:n_domains-1))
    call prim_advance_init(integration)
    call Prim_Advec_Init()
    call diffusion_init()
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
    use time_mod, only : timelevel_t, tstep, phys_tscale, timelevel_init, time_at, nendstep, smooth, nsplit, TimeLevel_Qdp
    use prim_state_mod, only : prim_printstate, prim_diag_scalars
    use filter_mod, only : filter_t, fm_filter_create, taylor_filter_create, &
         fm_transfer, bv_transfer
    use control_mod, only : runtype, integration, filter_mu, filter_mu_advection, test_case, &
         debug_level, vfile_int, filter_freq, filter_freq_advection, &
         transfer_type, vform, vfile_mid, filter_type, kcut_fm, wght_fm, p_bv, &
         s_bv, topology,columnpackage, moisture, precon_method, qsplit, rk_stage_user,&
         TRACERADV_TOTAL_DIVERGENCE, TRACERADV_UGRADQ, sub_case, &
         use_cpstar, energy_fixer, limiter_option, nu, nu_q, nu_div, tstep_type, hypervis_subcycle, &
         hypervis_subcycle_q
    use prim_si_ref_mod, only: prim_si_refstate_init, prim_set_mass
    use thread_mod, only : nthreads
    use derivative_mod, only : derivinit, interpolate_gll2fvm_points, interpolate_gll2spelt_points, v2pinit
    use global_norms_mod, only : test_global_integral, print_cfl
    use hybvcoord_mod, only : hvcoord_t
#ifdef CAM
#else
    use column_model_mod, only : InitColumnModel
    use held_suarez_mod, only : hs0_init_state
    use baroclinic_inst_mod, only : binst_init_state, jw_baroclinic
    use asp_tests, only : asp_tracer, asp_baroclinic, asp_rossby, asp_mountain, asp_gravity_wave 
    use aquaplanet, only : aquaplanet_init_state
#endif
#ifdef _ACCEL
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
    real (kind=real_kind) :: dt_dyn          ! viscosity timestep used in dynamics
    real (kind=real_kind) :: dt_tracers      ! viscosity timestep used in tracers

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
    real (kind=longdouble_kind)                 :: refnep(1:nep), tmp
    

    ! ==========================
    ! begin executable code
    ! ==========================
    if (topology == "cube") then
       call test_global_integral(elem, hybrid,nets,nete)  
    end if


    ! compute most restrictive dt*nu: 
    if (tstep_type.eq.0) then
       ! LF case, dynamics & tracers both use leapfrog with same timestep
       ! because of LF, timestep seen by viscosity is 2*tstep
       !dtnu = 2.0d0*tstep*max(nu/hypervis_subcycle,nu_q/hypervis_subcycle_q)
       dtnu = 2.0d0*tstep*max(nu,nu_div,nu_q)
    endif
    if (tstep_type.eq.1) then
       ! compute timestep seen by viscosity operator:
       if (qsplit==1) then
          dt_dyn = tstep   ! dynamics uses pure RK2 method
       else
          dt_dyn = 2*tstep  ! dynamics uses RK2 + (qsplit-1) LF.  LF viscosity uses 2dt. 
       endif
       dt_tracers=tstep*qsplit

       ! compute most restrictive condition:
       dtnu=max(dt_dyn*max(nu,nu_div), dt_tracers*nu_q)
    endif
    ! timesteps to use for advective stability:  tstep*qsplit and tstep
    call print_cfl(elem,hybrid,nets,nete,dtnu,tstep*qsplit,tstep)



    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call derivinit(deriv(hybrid%ithr),fvm_corners, fvm_points)

    ! ================================================
    ! fvm initialization
    ! ================================================
    if (ntrac>0) then
#if defined(_SPELT)
      tmp=nep-1
      do i=1,nep
        refnep(i)= 2*(i-1)/tmp - 1
      end do
      call v2pinit(deriv(hybrid%ithr)%Sfvm,gp%points,refnep,np,nep)
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
#if (! defined ELEMENT_OPENMP)
!$OMP MASTER
#endif
    if (integration == "semi_imp") then
       refstate = prim_si_refstate_init(.false.,hybrid%masterthread,hvcoord)
       if (precon_method == "block_jacobi") then
          allocate(blkjac(nets:nete))
       endif
    endif
#if (! defined ELEMENT_OPENMP)
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

#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif
    if (hybrid%ithr==0) then
       call syncmp(hybrid%par)
    end if
#if (! defined ELEMENT_OPENMP)
    !$OMP BARRIER
#endif

    if (topology /= "cube") then
       call abortmp('Error: only cube topology supported for primaitve equations') 
    endif


#ifndef CAM
    ! =================================
    ! HOMME stand alone initialization
    ! =================================
    tl%nstep0=2   ! This will be the first full leapfrog step
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
       
       tl%nstep0=tl%nstep+1            ! restart run: first step = first first full leapfrog step
       
       if (runtype==2) then
          ! branch run
          ! reset time counters to zero since timestep may have changed
          nEndStep = nEndStep-tl%nstep ! restart set this to nmax + tl%nstep
          tl%nstep=0
          tl%nstep0=2   
          ! copy prognostic variables:  tl%n0 into tl%nm1
          do ie=nets,nete
             elem(ie)%state%v(:,:,:,:,tl%nm1)=elem(ie)%state%v(:,:,:,:,tl%n0)
             elem(ie)%state%T(:,:,:,tl%nm1)=elem(ie)%state%T(:,:,:,tl%n0)
             elem(ie)%state%ps_v(:,:,tl%nm1)=elem(ie)%state%ps_v(:,:,tl%n0)
             elem(ie)%state%lnps(:,:,tl%nm1)=elem(ie)%state%lnps(:,:,tl%n0)

          enddo
       endif ! runtype==2
       if (hybrid%masterthread) then 
          write(iulog,*) "initial state from restart file:"
          write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
       end if
       
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
       else
          call abortmp('Error: unrecognized test case') 
       endif
       
       if (hybrid%masterthread) then
          write(iulog,*) '...done'
       end if
       
       ! scale PS to achieve prescribed dry mass
       call prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)  
       
       ! ========================================
       ! Print state and movie output
       ! ========================================
       
       if (hybrid%masterthread) then 
          write(iulog,*) "initial state:"
          write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
       end if
    end if  ! runtype

#endif

    ! For new runs, and branch runs, convert state variable to (Qdp)
    ! because initial conditon reads in Q, not Qdp
    ! restart runs will read dpQ from restart file
    ! need to check what CAM does on a branch run
    if (runtype==0 .or. runtype==2) then
       do ie=nets,nete
          elem(ie)%derived%omega_p(:,:,:) = 0D0
       end do

       call TimeLevel_Qdp( tl, qsplit, n0_qdp)
       do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, j, i, t, q, dp)
#endif
          do k=1,nlev    !  Loop inversion (AAM)
             do t=1,3
                do q=1,qsize       
                   do i=1,np
                      do j=1,np          
                         dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                              ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,t)
                         
                         elem(ie)%state%Qdp(i,j,k,q,n0_qdp)=elem(ie)%state%Q(i,j,k,q)*dp  
                         
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
    
 ! do it only for FVM tracers, FIRST TRACER will be the AIR DENSITY   
 ! should be optimize and combined with the above caculation 
    if (ntrac>0) then 
#if defined(_SPELT)
      ! do it only for FVM tracers, FIRST TRACER will be the AIR DENSITY   
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
          fvm(ie)%c(1:nep,1:nep,k,2,tl%n0)=interpolate_gll2spelt_points(elem(ie)%derived%dp(:,:,k),deriv(hybrid%ithr))
          do i=1,nep
            do j=1,nep
              fvm(ie)%c(i,j,k,2,tl%n0)=fvm(ie)%sga(i,j)*fvm(ie)%c(i,j,k,2,tl%n0)
!            fvm(ie)%c(i,j,k,1,tl%n0)=fvm(ie)%sga(i,j)
            end do
          end do
        enddo
      enddo
      call spelt_init3(elem,fvm,hybrid,nets,nete,tl%n0)
      do ie=nets,nete 
   	    do i=1-nipm,nep+nipm
   	      do j=1-nipm,nep+nipm  
   	        fvm(ie)%psc(i,j) = fvm(ie)%sga(i,j)*(sum(fvm(ie)%c(i,j,:,2,tl%n0)/fvm(ie)%sga(i,j)) +  hvcoord%hyai(1)*hvcoord%ps0)
   	      enddo
   	    enddo
      enddo
      
      if (hybrid%masterthread) then
         write(iulog,*) 'FVM (Spelt) tracers (incl. in halo zone) initialized. FIRST tracer has air density!'
      end if
#else
      ! do it only for FVM tracers, FIRST TRACER will be the AIR DENSITY   
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
          fvm(ie)%c(1:nc,1:nc,k,1,tl%n0)=interpolate_gll2fvm_points(elem(ie)%derived%dp(:,:,k),deriv(hybrid%ithr))
!            fvm(ie)%c(:,:,k,1,tl%n0)=1.0D0
        enddo
      enddo
      call fvm_init3(elem,fvm,hybrid,nets,nete,tl%n0)
      do ie=nets,nete 
   	    do i=1-nhc,nc+nhc
   	      do j=1-nhc,nc+nhc  
   	        fvm(ie)%psc(i,j) = sum(fvm(ie)%c(i,j,:,1,tl%n0)) +  hvcoord%hyai(1)*hvcoord%ps0
   	      enddo
   	    enddo
      enddo
      if (hybrid%masterthread) then
         write(iulog,*) 'FVM tracers (incl. in halo zone) initialized. FIRST tracer has air density!'	       
      end if
#endif      
    endif  

    ! for restart runs, we read in Qdp for exact restart, and rederive Q
    if (runtype==1) then
       call TimeLevel_Qdp( tl, qsplit, n0_qdp)
       do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, j, i, t, q, dp)
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

#ifdef CAM
    if (hybrid%masterthread) then 
       ! CAM has set tstep based on dtime before calling prim_init2(), 
       ! so only now does HOMME learn the timstep.  print them out:
       if (phys_tscale/=0) then
          write(iulog,'(a,2f9.2)') "CAM physics timescale:       ",phys_tscale
       endif
       write(iulog,'(a,2f9.2)') "CAM dtime (dt_phys):         ",tstep*nsplit*qsplit
       write(iulog,'(a,2f9.2)') "CAM dt_tracer, per RK stage: ",tstep*qsplit,(tstep*qsplit)/(rk_stage_user-1)
       write(iulog,'(a,2f9.2)') "CAM dt_dyn, per RK stage:    ",tstep*qsplit,tstep
 

       write(iulog,*) "initial state from CAM:"
       write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
    end if
#endif

#ifdef _ACCEL
    !Inside this routine, we enforce an OMP BARRIER and an OMP MASTER. It's left out of here because it's ugly
    call cuda_mod_init(elem,deriv(hybrid%ithr),hvcoord)
#endif

    call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete, fvm)
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
  use time_mod, only : TimeLevel_t, time_at
  use control_mod, only : tstep_type

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
    use time_mod, only : TimeLevel_t, time_at, timelevel_update, smooth
    use control_mod, only: statefreq, integration, &
           TRACERADV_TOTAL_DIVERGENCE,TRACERADV_UGRADQ, ftype, tstep_type, nu_p, qsplit
    use prim_advance_mod, only : prim_advance_exp, prim_advance_si, preq_robert3
    use prim_state_mod, only : prim_printstate, prim_diag_scalars, prim_energy_halftimes
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
    integer:: np1_qdp, no_qdp


    logical :: compute_diagnostics

    ! ===================================
    ! Main timestepping loop
    ! ===================================

    ! compute diagnostics and energy for STDOUT 
    ! compute energy if we are using an energy fixer

    compute_diagnostics=.false.
    if (MODULO(tl%nstep+1,statefreq)==0 .or. tl%nstep+1==tl%nstep0) then
       compute_diagnostics=.true.  
    endif
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
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
    if (integration == "explicit") then
       call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
            flt , hybrid,             &
            dt, tl, nets, nete, compute_diagnostics)
    else if (integration == "semi_imp") then
       call prim_advance_si(elem, nets, nete, cg(hybrid%ithr), blkjac, red, &
            refstate, hvcoord, deriv(hybrid%ithr), flt, hybrid, tl, dt)
       tot_iter=tot_iter+cg(hybrid%ithr)%iter
    end if


    ! ===============
    ! Tracer Advection  Needs U,V at timelevel n0 and eta_dot_dpdn at timellevel n0
    ! and maybe timelevel np1 which was computed in dynamics step above.  
    ! ===============
    !!!!!! NOTE: no longer suppoting advecting tracersw/leapfrog code


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
    ! timestep is complete.  Now apply robert filter to all prognostic variables
    ! =================================
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
          write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
          if (integration == "semi_imp") write(iulog,*) "cg its=",cg(0)%iter
       end if
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)
    end if
  end subroutine prim_run

!=======================================================================================================! 


  subroutine prim_run_subcycle(elem, fvm, hybrid,nets,nete, dt, tl, hvcoord)
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
    use time_mod, only : TimeLevel_t, time_at, timelevel_update, nsplit
    use control_mod, only: statefreq,&
           energy_fixer, ftype, qsplit, rsplit, nu_p, test_cfldep
    use prim_advance_mod, only : prim_advance_exp, prim_advance_si, preq_robert3, applycamforcing, &
                                 applycamforcing_dynamics, prim_advance_exp
    use prim_state_mod, only : prim_printstate, prim_diag_scalars, prim_energy_halftimes
    use parallel_mod, only : abortmp
    use reduction_mod, only : parallelmax
    use prim_advection_mod, only : vertical_remap

    

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
    if (compute_diagnostics) then
       call prim_diag_scalars(elem,hvcoord,tl,4,.true.,nets,nete)
       call prim_energy_halftimes(elem,hvcoord,tl,4,.true.,nets,nete,tl%n0)
    endif


#ifdef CAM
    ! ftype=2  Q was adjusted by physics, but apply u,T forcing here
    ! ftype=1  forcing was applied time-split in CAM coupling layer
    ! ftype=0 means forcing apply here
    ! ftype=-1 do not apply forcing
    call TimeLevel_Qdp( tl, qsplit, n0_qdp)
    if (ftype==0) call ApplyCAMForcing(elem, hvcoord,tl%n0,n0_qdp, dt_remap,nets,nete)
    if (ftype==2) call ApplyCAMForcing_dynamics(elem, hvcoord,tl%n0,dt_remap,nets,nete)
#endif

    !We are only storing Q now (just n0 - not nm1 and np1)  
    ! E(1) Energy at start of timestep, diagnostics at t-dt/2  (using t-dt, t and Q(t))
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,1,.true.,nets,nete,1)

    ! qmass and variance, using Q(n0),Qdp(n0)
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,1,.true.,nets,nete)

    if (rsplit>0) then
    ! vertically lagrangian code: initialize prognostic variable dp3d for floating levels
    do ie=nets,nete
      do k=1,nlev
         elem(ie)%state%dp3d(:,:,k,tl%n0)=&
              ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
              ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)

      enddo
    enddo
    endif

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
    !  also for dynamics if rsplit>0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !compute timelevels for tracers (no longer the same as dynamics)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp) 
    call vertical_remap(elem,fvm,hvcoord,dt_remap,tl%np1,np1_qdp,nets,nete)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! time step is complete.  update some diagnostic variables:
    ! lnps (we should get rid of this)
    ! Q    (mixing ratio)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
       !$omp parallel do private(k,q)
#endif
       elem(ie)%state%lnps(:,:,tl%np1)= LOG(elem(ie)%state%ps_v(:,:,tl%np1))
       do k=1,nlev    !  Loop inversion (AAM)
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
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,2,.false.,nets,nete,1)

    if (energy_fixer > 0) then
       if ( .not. compute_energy) call abortmp("ERROR: energy fixer needs compute_energy=.true")
       call prim_energy_fixer(elem,hvcoord,hybrid,tl,nets,nete)
    endif

    if (compute_diagnostics) then
       call prim_diag_scalars(elem,hvcoord,tl,3,.false.,nets,nete)
       call prim_energy_halftimes(elem,hvcoord,tl,3,.false.,nets,nete,1)
     endif

    ! =================================
    ! update dynamics time level pointers 
    ! =================================
    call TimeLevel_update(tl,"leapfrog")

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - dt       (Robert-filtered)
    !   u(n0)    dynamics at  t+dt_remap 
    !   u(np1)   undefined
 

    ! ============================================================
    ! Print some diagnostic information 
    ! ============================================================
    if (compute_diagnostics) then
       if (hybrid%masterthread) then 
          write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
       end if
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
    use time_mod, only : TimeLevel_t, time_at, timelevel_update, nsplit
    use control_mod, only: statefreq,&
           energy_fixer, ftype, qsplit, nu_p, test_cfldep, rsplit
    use prim_advance_mod, only : prim_advance_exp, prim_advance_si, preq_robert3, applycamforcing, &
                                 applycamforcing_dynamics, prim_advance_exp, overwrite_SEdensity
    use prim_advection_mod, only : prim_advec_tracers_remap_rk2, prim_advec_tracers_fvm, &
         prim_advec_tracers_spelt
    use prim_state_mod, only : prim_printstate, prim_diag_scalars, prim_energy_halftimes
    use parallel_mod, only : abortmp
    use reduction_mod, only : parallelmax

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
    integer :: ie, t, q,k,i,j,n, n_Q

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
         elem(ie)%derived%psdiss_ave=0
         elem(ie)%derived%psdiss_biharmonic=0
      endif
      if (ntrac>0) then
        ! save velocity at time t for fvm
        fvm(ie)%vn0=elem(ie)%state%v(:,:,:,:,tl%n0)
      end if
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, j, i)
#endif
      ! save dp at time t for use in tracers
      do k=1,nlev
         if (rsplit==0) then
            ! dp at time t is on REF levels
            elem(ie)%derived%dp(:,:,k)=&
                 ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                 ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
         else
            ! dp at time t:  use floating lagrangian levels:
            elem(ie)%derived%dp(:,:,k)=elem(ie)%state%dp3d(:,:,k,tl%n0)
         endif
      enddo
    enddo

    ! ===============
    ! Dynamical Step 
    ! ===============
    n_Q = tl%n0  ! n_Q = timelevel of FV tracers at time t.  need to save this
    call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
         flt , hybrid, dt, tl, nets, nete, compute_diagnostics)
    do n=2,qsplit
       call TimeLevel_update(tl,"leapfrog")
       call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
            flt , hybrid, dt, tl, nets, nete, .false.)
       ! defer final timelevel update until after Q update.
    enddo



    ! ===============
    ! Tracer Advection.  SE advection uses mean flux variables:
    !        derived%dp              =  dp at start of timestep
    !        derived%vn0             =  mean horiz. flux:   U*dp
    ! in addition, this routine will apply the DSS to:
    !        derived%eta_dot_dpdn    =  mean vertical velocity (used for remap below)
    !        derived%omega           =  
    ! ===============
    if (qsize>0) call Prim_Advec_Tracers_remap_rk2(elem, deriv(hybrid%ithr),hvcoord,flt_advection,hybrid,&
         dt_q,tl,nets,nete)


    if (ntrac>0) then
      if ( n_Q /= tl%n0 ) then
        ! make sure tl%n0 contains tracers at start of timestep
        do ie=nets,nete
          fvm(ie)%c(:,:,:,1:ntrac,tl%n0)  = fvm(ie)%c(:,:,:,1:ntrac,n_Q)
        enddo
      endif 
#if defined(_SPELT) 
      call Prim_Advec_Tracers_spelt(elem, fvm, deriv(hybrid%ithr),hvcoord,hybrid,&
           dt_q,tl,nets,nete)
       do ie=nets,nete 
    	    do i=1-nipm,nep+nipm
    	      do j=1-nipm,nep+nipm  
    	        fvm(ie)%psc(i,j) = fvm(ie)%sga(i,j)*(sum(fvm(ie)%c(i,j,:,2,tl%np1)/fvm(ie)%sga(i,j)) +  hvcoord%hyai(1)*hvcoord%ps0)
    	      enddo
    	    enddo
       enddo
#else      
      call Prim_Advec_Tracers_fvm(elem, fvm, deriv(hybrid%ithr),hvcoord,hybrid,&
           dt_q,tl,nets,nete)
           ! values in the halo zone are only in np1 at this time
       do ie=nets,nete 
         do i=1-nhc,nc+nhc
           do j=1-nhc,nc+nhc  
             fvm(ie)%psc(i,j) = sum(fvm(ie)%c(i,j,:,1,tl%np1)) +  hvcoord%hyai(1)*hvcoord%ps0
           enddo
         enddo
       enddo

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
!        call overwrite_SEdensity(elem,fvm,hybrid,nets,nete,tl%np1) 
!elem(ie)%state%ps_v(i,j,np1)
#endif
       ! dynamics computed a predictor surface pressure, now correct with fvm result
       ! fvm has computed a new dp(:,:,k) on the fvm grid
       !
       ! step 1:  computes sum of fvm density on fvm grid
       ! step 2:  interpolate to GLL grid
       ! step 3:  store in  elem(ie)%state%ps_v(:,:,tl%np1)
       ! step 4:  apply DSS to make ps_v continuous 
       !          (or use continuous reconstruction)
    endif

  end subroutine prim_step


!=======================================================================================================! 


  subroutine prim_finalize(hybrid)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    ! ==========================
    ! end of the hybrid program
    ! ==========================
  end subroutine prim_finalize

!=======================================================================================================! 
  subroutine prim_energy_fixer(elem,hvcoord,hybrid,tl,nets,nete)
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
    use physical_constants, only : Cp, cpwater_vapor,g,dd_pi
    use physics_mod, only : Virtual_Specific_Heat
    use time_mod, only : timelevel_t, TimeLevel_Qdp
    use control_mod, only : moisture, traceradv_total_divergence,energy_fixer, use_cpstar, qsplit
    use hybvcoord_mod, only : hvcoord_t
    use global_norms_mod, only: wrap_repro_sum
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)
    integer :: t1,t2,t3,n,nets,nete,t_beta, t_beta_qdp, t1_qdp, t2_qdp
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (TimeLevel_t), intent(inout)       :: tl

    integer :: ie,k,i,j,nm_f
    real (kind=real_kind), dimension(np,np,nlev)  :: dp   ! delta pressure
    real (kind=real_kind), dimension(np,np,nlev)  :: sumlk
    real (kind=real_kind), dimension(np,np)  :: suml
    real (kind=real_kind) :: cp_star,qval

    real (kind=real_kind) :: psum(nets:nete,3),psum_g(3),beta,scale
    logical :: use_cp_star

! energy_fixer
!     1          cp_star(t1)*dp(t1)*T(t2)         (staggered)
!     2          cp*dp(t1)*T(t2)                  (staggered)
!     3          cp_star(t2)*dp(t2)*T(t2)
!     4          cp*dp(t2)*T(t2)
!
    t1=tl%n0     ! timelevel for cp_star dp 
    t2=tl%np1    ! timelevel for T
    call TimeLevel_Qdp( tl, qsplit, t1_qdp, t2_qdp) 
    
    use_cp_star = (use_cpstar.eq.1)  ! by default, use global value
    ! but override to false for certain fixers:
    if (energy_fixer==2) use_cp_star = .false.   
    if (energy_fixer==4) use_cp_star = .false.
    
    t_beta = t1
    t_beta_qdp = t1_qdp
    scale = 2.0
    if (energy_fixer==3 .or. energy_fixer==4) then
       t_beta=t2
       t_beta_qdp = t2_qdp
       scale = 1.0
    endif
    ! with staggered-in-time formula:
    !    compute <dp cp_star> at t1   
    !    E(Tnew) =(< dp(t1) cp_star(t1) (T(t2)+beta) > + < dp(t2) cp_star(t2) T(t1) >)/2 
    !    E(Tnew)-E(t2) = < dp(t1) cp_star(t1) beta >/2
    ! with everyting at t2: 
    !    compute <dp cp_star> at t2   
    !    E(Tnew)-E(t2) = < dp(t2) cp_star(t2) beta >
    !

    psum = 0
    do ie=nets,nete

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
          dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,t_beta)
       enddo
       suml=0
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, j, i, qval, cp_star)
#endif
       do k=1,nlev
          do i=1,np
          do j=1,np
             if(use_cp_star)  then
                qval = elem(ie)%state%Qdp(i,j,k,1,t_beta_qdp)/dp(i,j,k)
                cp_star= Virtual_Specific_Heat(qval)
                sumlk(i,j,k) = cp_star*dp(i,j,k) 
             else
                sumlk(i,j,k) = cp*dp(i,j,k) 
             endif

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
       psum(ie,3) = psum(ie,3) + SUM(suml(:,:)*elem(ie)%spheremp(:,:))
       do n=1,2
          if (use_cp_star ) then
             psum(ie,n) = psum(ie,n) + SUM(  elem(ie)%spheremp(:,:)*&
                                    (elem(ie)%accum%PEner(:,:,n) + &
                                    elem(ie)%accum%IEner(:,:,n) + &
                                    elem(ie)%accum%KEner(:,:,n) ) )
          else
             psum(ie,n) = psum(ie,n) + SUM(  elem(ie)%spheremp(:,:)*&
                                    (elem(ie)%accum%PEner(:,:,n) + &
              (elem(ie)%accum%IEner(:,:,n)-elem(ie)%accum%IEner_wet(:,:,n)) + &
                                    elem(ie)%accum%KEner(:,:,n) ) )
          endif
       enddo
    enddo
    do ie=nets,nete
       do n=1,3
          global_shared_buf(ie,n) = psum(ie,n)
       enddo
    enddo
    call wrap_repro_sum(nvars=3, comm=hybrid%par%comm)
    do n=1,3
       psum_g(n) = global_shared_sum(n)
    enddo
    beta = scale*( psum_g(1)-psum_g(2) )/psum_g(3)
    do ie=nets,nete
       elem(ie)%state%T(:,:,:,t2) =  elem(ie)%state%T(:,:,:,t2) + beta
    enddo

!=======================================================================================================! 
    end subroutine prim_energy_fixer


    subroutine smooth_topo_datasets(phis,sghdyn,sgh30dyn,elem,hybrid,nets,nete)
    use control_mod, only : smooth_phis_numcycle,smooth_sgh_numcycle
    use hybrid_mod, only : hybrid_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
    use bndry_mod, only : bndry_exchangev
    use derivative_mod, only : derivative_t , laplace_sphere_wk
    use viscosity_mod, only : biharmonic_wk
    use time_mod, only : TimeLevel_t
    use prim_advance_mod, only : smooth_phis
    implicit none
    
    real (kind=real_kind), intent(inout)   :: phis(np,np,nets:nete)
    real (kind=real_kind), intent(inout)   :: sghdyn(np,np,nets:nete)
    real (kind=real_kind), intent(inout)   :: sgh30dyn(np,np,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    integer , intent(in) :: nets,nete
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



