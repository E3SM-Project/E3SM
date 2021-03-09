! ------------------------------------------------------------------------------------------------
! prim_driver_mod: 
!
! Revisions:
! 08/2016: O. Guba Inserting code for "espilon bubble" reference element map
! 03/2018: M. Taylor  fix memory leak
! 06/2018: O. Guba  code for new ftypes
! 06/2019: M. Taylor remove ps_v
!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_driver_base

  use derivative_mod,   only: derivative_t, derivinit
  use dimensions_mod,   only: np, nlev, nlevp, nelem, nelemd, nelemdmax, GlobalUniqueCols, qsize
  use element_mod,      only: element_t, allocate_element_desc, setup_element_pointers
  use element_ops,      only: copy_state
  use gridgraph_mod,    only: GridVertex_t, GridEdge_t
  use hybrid_mod,       only: hybrid_t
  use kinds,            only: real_kind, iulog
  use metagraph_mod,    only: MetaVertex_t
  use perf_mod,         only: t_startf, t_stopf
  use quadrature_mod,   only: quadrature_t, gausslobatto
  use reduction_mod,    only: reductionbuffer_ordered_1d_t, red_min, red_max, red_max_int, &
                              red_sum, red_sum_int, red_flops, initreductionbuffer, &
                              red_max_index, red_min_index
#ifndef CAM
#ifndef SCREAM
  use prim_restart_mod, only : initrestartfile
  use restart_io_mod ,  only : readrestart
#endif
! For SCREAM, we might do unit tests, so enable setting of initial conditions
  use test_mod,         only: set_test_initial_conditions, compute_test_forcing
#endif

  implicit none

  private
  public :: prim_init1, prim_init2 , prim_run_subcycle, prim_finalize
  public :: prim_init1_geometry, prim_init1_elem_arrays, prim_init1_buffers, prim_init1_cleanup
#ifndef CAM
  public :: prim_init1_no_cam
#endif

  public :: smooth_topo_datasets, deriv1

  public :: applyCAMforcing_tracers

  ! Service variables used to partition the mesh.
  ! Note: GridEdge and MeshVertex are public, cause kokkos targets need to access them
  type (GridVertex_t), pointer :: GridVertex(:)
  type (GridEdge_t),   public, pointer :: GridEdge(:)
  type (MetaVertex_t), public :: MetaVertex
  logical :: can_scalably_init_grid

  type (quadrature_t)   :: gp                     ! element GLL points
  type (ReductionBuffer_ordered_1d_t), save :: red   ! reduction buffer               (shared)
  type (derivative_t)  :: deriv1

contains

  subroutine prim_init1(elem, par, dom_mt, Tl)
    use domain_mod,    only : domain1d_t
    use parallel_mod,  only : parallel_t
    use time_mod,      only : TimeLevel_t, TimeLevel_init
    !
    ! Inputs
    !
    type (element_t),   pointer     :: elem(:)
    type (parallel_t),  intent(in)  :: par
    type (domain1d_t),  pointer     :: dom_mt(:)
    type (timelevel_t), intent(out) :: Tl
    !
    ! Locals
    !

    ! Note: resist the temptation to collapse these routines into prim_init1!
    !       They have been split up to allow kokkos targets to work correctly.
    !       For instance, kokkos builds do not need to init communication buffers,
    !       since they have their own C++ based communication system. Also,
    !       the C++ communication graph is built from the Fortran one.
    !       So, the prim_init1 routine in the kokkos targets' prim_driver_mod
    !       will call *some* of the following pieces, but not all. Moreover,
    !       we pulled GridEdge and MetaVertex out of the prim_init1_geometry
    !       chunk, so that the kokkos targets can do the following:
    !         a) call prim_init1_geometry
    !         b) use the info in GridEdge and MetaVertex to build a C++
    !            communication structure
    !         c) call prim_init1_cleanup to deallocate GridEdge and MetaVertex
    !       The overall behavior is the same as it was before, when there was
    !       a single big prim_init1 subroutine. However, this split-up
    !       implementation allows kokkos targets to avoid unnecessary pieces,
    !       as well as to inject code in between pieces that is needed to
    !       properly setup the C++ structures.

#ifndef CAM
    ! Initialize a few things that CAM would take care of (e.g., parsing namelist)
    call prim_init1_no_cam (par)
#endif

    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call derivinit(deriv1)

    ! ==================================
    ! Initialize and partition the geometry
    ! ==================================
    call prim_init1_geometry(elem,par,dom_mt)

    ! ==================================
    ! Initialize element pointers (if any)
    ! ==================================
    ! TODO: this is for OPENACC only. preqx_acc should define its own prim_init1,
    !       which should do the same things as in this prim_init1, including
    !       the call to setup_element_pointers, which should be removed from
    !       the base version of prim_init1
    call setup_element_pointers(elem)

    ! ==================================
    ! Initialize element arrays (fluxes and state)
    ! ==================================
    call prim_init1_elem_arrays(elem,par)

    call prim_init1_compose(par,elem)

    ! Cleanup the tmp stuff used in prim_init1_geometry
    call prim_init1_cleanup()

    ! ==================================
    ! Initialize the buffers for exchanges
    ! ==================================
    call prim_init1_buffers(elem,par)

    ! Initialize the time levels
    call TimeLevel_init(tl)

    if(par%masterproc) write(iulog,*) 'end of prim_init1'
  end subroutine prim_init1


#ifndef CAM
  subroutine prim_init1_no_cam(par)
    use mesh_mod,       only : MeshUseMeshFile, MeshCubeElemCount
    use cube_mod,       only : CubeElemCount
    use planar_mod,     only : PlaneElemCount
    use parallel_mod,   only : parallel_t, abortmp
    use namelist_mod,   only : readnl
    use quadrature_mod, only : test_gauss, test_gausslobatto
    use repro_sum_mod,  only : repro_sum_defaultopts, repro_sum_setopts
    use time_mod,       only : nmax, time_at
    use control_mod, only : topology
#ifndef HOMME_WITHOUT_PIOLIBRARY
    use common_io_mod,  only : homme_pio_init
#endif
    !
    ! Inputs
    !
    type (parallel_t),  intent(in)  :: par
    !
    ! Locals
    !
    logical :: repro_sum_use_ddpdd, repro_sum_recompute
    real(kind=real_kind) :: repro_sum_rel_diff_max
    real(kind=real_kind) :: approx_elements_per_task
    integer :: total_nelem

    ! =====================================
    ! Read in model control information
    ! =====================================
    ! cam readnl is called in spmd_dyn (needed prior to mpi_init)
    call readnl(par)
    if (MeshUseMeshFile) then
       total_nelem = MeshCubeElemCount()
    else
      if (topology=="cube") then
       total_nelem = CubeElemCount()
     else if (topology=="plane") then
       total_nelem      = PlaneElemCount()
    end if
    end if

#ifndef HOMME_WITHOUT_PIOLIBRARY
    call homme_pio_init(par%rank,par%comm)
#endif

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
  end subroutine prim_init1_no_cam
#endif

  subroutine prim_init1_geometry(elem, par, dom_mt)

    ! --------------------------------
    use thread_mod, only : nthreads, hthreads, vthreads
    ! --------------------------------
    use control_mod, only : topology, geometry, partmethod, z2_map_method, cubed_sphere_map
    ! --------------------------------
    use prim_state_mod, only : prim_printstate_init
    ! --------------------------------
    use mass_matrix_mod, only : mass_matrix
    ! --------------------------------
    use cube_mod,  only : cubeedgecount , cubeelemcount, cubetopology, cube_init_atomic, &
                          set_corner_coordinates
    ! --------------------------------
    use geometry_mod, only: set_area_correction_map0, set_area_correction_map2
    ! --------------------------------
    use mesh_mod, only : MeshSetCoordinates, MeshUseMeshFile, MeshCubeTopology, &
                         MeshCubeElemCount, MeshCubeEdgeCount, MeshCubeTopologyCoords
    ! --------------------------------
    use planar_mod,  only : PlaneEdgeCount , PlaneElemCount, PlaneTopology
    ! --------------------------------
    use planar_mesh_mod, only :   PlaneMeshSetCoordinates,      &
                          MeshPlaneTopology, MeshPlaneTopologyCoords
    ! --------------------------------
    use planar_mod, only : plane_init_atomic, plane_set_corner_coordinates
    ! --------------------------------
    use metagraph_mod, only : localelemcount, initmetagraph, printmetavertex
    ! --------------------------------
    use gridgraph_mod, only : allocate_gridvertex_nbrs
    ! --------------------------------
    use schedtype_mod, only : schedule
    ! --------------------------------
    use schedule_mod, only : genEdgeSched,  PrintSchedule
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
    use scalable_grid_init_mod, only : sgi_init_grid
    ! --------------------------------
    use dof_mod, only : global_dof, CreateUniqueIndex, SetElemOffset
    ! --------------------------------
    use params_mod, only : SFCURVE
    ! --------------------------------
    use zoltan_mod, only: genzoltanpart, getfixmeshcoordinates, printMetrics, is_zoltan_partition, is_zoltan_task_mapping
    ! --------------------------------
    use domain_mod, only : domain1d_t, decompose
    ! --------------------------------
    use physical_constants, only : dd_pi
    ! --------------------------------

    implicit none
    !
    ! Locals
    !

    type (element_t),   pointer     :: elem(:)
    type (parallel_t),  intent(in)  :: par
    type (domain1d_t),  pointer     :: dom_mt(:)

    integer :: ii,ie, ith
    integer :: nelem_edge,nedge
    integer :: nstep
    integer :: nlyr
    integer :: iMv
    integer :: err, ierr, l, j
    logical, parameter :: Debug = .FALSE.

    integer  :: i

    type (quadrature_t)   :: gp                     ! element GLL points


    real (kind=real_kind) ,  allocatable :: coord_dim1(:)
    real (kind=real_kind) ,  allocatable :: coord_dim2(:)
    real (kind=real_kind) ,  allocatable :: coord_dim3(:)
    integer :: coord_dimension = 3

    ! ===============================================================
    ! Allocate and initialize the graph (array of GridVertex_t types)
    ! ===============================================================
    if (MeshUseMeshFile) then
       nelem = MeshCubeElemCount()
       nelem_edge = MeshCubeEdgeCount()
    else if (topology=="cube") then
       nelem      = CubeElemCount()
       nelem_edge = CubeEdgeCount()
     else if (topology=="plane") then
       nelem      = PlaneElemCount()
       nelem_edge = PlaneEdgeCount()
    end if

    ! we want to exit elegantly when we are using too many processors.
    if (nelem < par%nprocs) then
       call abortmp('Error: too many MPI tasks. set dyn_npes <= nelem')
    end if

    can_scalably_init_grid = &
         topology == "cube" .and. &
         .not. MeshUseMeshFile .and. &
         partmethod .eq. SFCURVE .and. &
         .not. (is_zoltan_partition(partmethod) .or. is_zoltan_task_mapping(z2_map_method))

    if (can_scalably_init_grid) then
       call sgi_init_grid(par, GridVertex, GridEdge, MetaVertex)
    end if

    if (.not. can_scalably_init_grid) then

       if (par%masterproc) then
          write(iulog,*)"creating topology..."
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
           if (topology=="cube") then
           call MeshCubeTopologyCoords(GridEdge, GridVertex, coord_dim1, coord_dim2, coord_dim3, coord_dimension)
           !MD:TODO: still need to do the coordinate transformation for this case.
         else if  (topology=="plane") then
           call MeshPlaneTopologyCoords(GridEdge,GridVertex, coord_dim1, coord_dim2, coord_dim3, coord_dimension)
         end if

       else
         if (topology=="cube") then
           call CubeTopology(GridEdge,GridVertex)
           if (is_zoltan_partition(partmethod) .or. is_zoltan_task_mapping(z2_map_method)) then
              call getfixmeshcoordinates(GridVertex, coord_dim1, coord_dim2, coord_dim3, coord_dimension)
           endif
        else if (topology=="plane") then
         call PlaneTopology(GridEdge,GridVertex)
       end if
        end if

       if(par%masterproc)       write(iulog,*)"...done."
    end if
    if(par%masterproc) write(iulog,*)"total number of elements nelem = ",nelem

    !DBG if(par%masterproc) call PrintGridVertex(GridVertex)

    call t_startf('PartitioningTime')

    if (.not. can_scalably_init_grid) then
       if(partmethod .eq. SFCURVE) then
          if(par%masterproc) write(iulog,*)"partitioning graph using SF Curve..."
          !if the partitioning method is space filling curves
          call genspacepart(GridEdge,GridVertex)
          if (is_zoltan_task_mapping(z2_map_method)) then
             if(par%masterproc) write(iulog,*)"mapping graph using zoltan2 task mapping on the result of SF Curve..."
             call genzoltanpart(GridEdge,GridVertex, par%comm, coord_dim1, coord_dim2, coord_dim3, coord_dimension)
          endif
          !if zoltan2 partitioning method is asked to run.
       elseif ( is_zoltan_partition(partmethod)) then
          if(par%masterproc) write(iulog,*)"partitioning graph using zoltan2 partitioning/task mapping..."
          call genzoltanpart(GridEdge,GridVertex, par%comm, coord_dim1, coord_dim2, coord_dim3, coord_dimension)
       else
          if(par%masterproc) write(iulog,*)"partitioning graph using Metis..."
          call genmetispart(GridEdge,GridVertex)
       endif
    endif ! .not. can_scalably_init_grid

    call t_stopf('PartitioningTime')

    !call t_startf('PrintMetricTime')
    !print partitioning and mapping metrics
    !call printMetrics(GridEdge,GridVertex, par%comm)
    !call t_stopf('PrintMetricTime')

    ! ===========================================================
    ! given partition, count number of local element descriptors
    ! ===========================================================
    allocate(Schedule(1))

    nelem_edge=SIZE(GridEdge)

    ! ====================================================
    !  Generate the communication graph
    ! ====================================================
    if (.not. can_scalably_init_grid) then
       call initMetaGraph(iam,MetaVertex,GridVertex,GridEdge)
    end if

    nelemd = LocalElemCount(MetaVertex)
    if(par%masterproc .and. Debug) then 
        call PrintMetaVertex(MetaVertex)
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
       call allocate_element_desc(elem)
    endif

    ! ====================================================
    !  Generate the communication schedule
    ! ====================================================

    call genEdgeSched(elem,iam,Schedule(1),MetaVertex)


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
    !
    ! At this point, we can assume: 
    ! nthreads was set by CAM driver, or in namelist and error checked
    ! if CAM is running w/o threads, nthreads=0 
    ! vthreads=1 or read from namelist and verified consistent with COLUMN_OPENMP
    ! 
    ! set hthreads, and check that vthreads was not set too large
    if (vthreads > max(nthreads,1) .or. vthreads < 1) &
         call abortmp('Error: vthreads<1 or vthreads > NTHRDS_ATM')
#if defined(HORIZ_OPENMP) && defined (COLUMN_OPENMP)
    if (vthreads>1) call omp_set_nested(.true.)
#endif
    hthreads = max(nthreads,1) / vthreads
    hthreads = min(max(nthreads,1),nelemd)
    if(par%masterproc) then
       write(iulog,*) "prim_init1: total threads: nthreads = ",nthreads
       write(iulog,*) "threads across elements    hthreads = ",hthreads
       write(iulog,*) "threading within elements  vthreads = ",vthreads
    endif

#ifndef COLUMN_OPENMP
    if (vthreads>1) call abortmp('Error: vthreads>1 requires -DCOLUMN_OPENMP')
#endif
#ifndef HORIZ_OPENMP
    if (hthreads>1) call abortmp('Error: hthreads>1 requires -DHORIZ_OPENMP')
#endif
    


    ! =================================================================
    ! Initialize shared boundary_exchange and reduction buffers
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'init shared boundary_exchange buffers'
    call InitReductionBuffer(red,3*nlev,hthreads)
    call InitReductionBuffer(red_sum,5)
    call InitReductionBuffer(red_sum_int,1)
    call InitReductionBuffer(red_max,1)
    call InitReductionBuffer(red_max_int,1)
    call InitReductionBuffer(red_min,1)
    call initReductionBuffer(red_flops,1)
    call initReductionBuffer(red_min_index,2)
    call initReductionBuffer(red_max_index,2)

    gp=gausslobatto(np)  ! GLL points

        if(par%masterproc) write(6,*)"initializing elements..."
    
         if (MeshUseMeshFile) then
          if (geometry=="sphere") then
            call MeshSetCoordinates(elem)
          else if  (geometry=="plane") then
            call PlaneMeshSetCoordinates(elem)
          end if
         else
          if (geometry=="sphere") then
            do ie=1,nelemd
               call set_corner_coordinates(elem(ie))
            enddo
          else if (geometry=="plane") then
           do ie=1,nelemd
              call plane_set_corner_coordinates(elem(ie))
           enddo
          end if
           !call assign_node_numbers_to_elem(elem, GridVertex)
         endif
    
      if (geometry=="sphere") then
         do ie=1,nelemd
            call cube_init_atomic(elem(ie),gp%points)
         enddo
      else if (geometry=="plane") then
        do ie=1,nelemd
           call plane_init_atomic(elem(ie),gp%points)
        enddo
      end if

    ! =================================================================
    ! Initialize mass_matrix
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'running mass_matrix'
    call mass_matrix(par,elem)

    ! global area correction (default for cubed-sphere meshes)
    if( cubed_sphere_map == 0 ) then
       call set_area_correction_map0(elem, nelemd, par, gp)
    endif

    ! Epsilon bubble correction (default for RRM meshes).
    if(( cubed_sphere_map == 2 ).AND.( np > 2 )) then
       call set_area_correction_map2(elem, nelemd, par, gp)
    endif

    deallocate(gp%points)
    deallocate(gp%weights)


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

    allocate(dom_mt(0:hthreads-1))
    do ith=0,hthreads-1
       dom_mt(ith)=decompose(1,nelemd,hthreads,ith)
    end do

  end subroutine prim_init1_geometry

  subroutine prim_init1_elem_arrays (elem,par)
    ! --------------------------------
    use prim_state_mod, only : prim_printstate_init
    use parallel_mod,   only : parallel_t
    use control_mod,    only : runtype, restartfreq, transport_alg
    use bndry_mod,      only : sort_neighbor_buffer_mapping
#if !defined(CAM) && !defined(SCREAM)
    use restart_io_mod, only : RestFile,readrestart
#endif

    implicit none
    !
    ! Inputs
    !
    type (element_t),   pointer     :: elem(:)
    type (parallel_t),  intent(in)  :: par
    !
    ! Locals
    !
    integer :: ie

    ! Initialize output fields for plotting...
    call prim_printstate_init(par, elem)

    ! initialize flux terms to 0
    do ie=1,nelemd
       elem(ie)%derived%FM=0.0
       elem(ie)%derived%FQ=0.0
       elem(ie)%derived%FQps=0.0
       elem(ie)%derived%FT=0.0

       elem(ie)%derived%Omega_p=0
    enddo

    ! ==========================================================
    !  This routines initalizes a Restart file.  This involves:
    !      I)  Setting up the MPI datastructures
    ! ==========================================================
#if !defined(CAM) && !defined(SCREAM)
    if(restartfreq > 0 .or. runtype>=1)  then
       call initRestartFile(elem(1)%state,par,RestFile)
    endif
#endif

    if (transport_alg > 0) then
      call sort_neighbor_buffer_mapping(par, elem,1,nelemd)
    end if

  end subroutine prim_init1_elem_arrays

  subroutine prim_init1_compose(par, elem)
    use parallel_mod, only : parallel_t, abortmp
    use control_mod,  only : transport_alg, semi_lagrange_cdr_alg
#ifdef HOMME_ENABLE_COMPOSE
    use compose_mod,  only : kokkos_init, compose_init, cedr_set_ie2gci, cedr_unittest
#endif

    type (parallel_t), intent(in) :: par
    type (element_t), pointer, intent(in) :: elem(:)
    integer :: ie, ierr

    if (transport_alg > 0) then
#ifdef HOMME_ENABLE_COMPOSE
       call compose_init(par, elem, GridVertex)
       do ie = 1, nelemd
          call cedr_set_ie2gci(ie, elem(ie)%vertex%number)
       end do
#else
       call abortmp('COMPOSE SL transport was requested, but HOMME was built without COMPOSE.')
#endif
    end if
  end subroutine prim_init1_compose

  subroutine prim_init1_cleanup ()
    use gridgraph_mod, only : deallocate_gridvertex_nbrs
    use metagraph_mod, only : destroyMetaGraph
    use scalable_grid_init_mod, only : sgi_finalize

    integer :: j

    if (can_scalably_init_grid) then
       call sgi_finalize()
    else
       deallocate(GridEdge)
       call destroyMetaGraph(MetaVertex)
       do j =1,nelem
          call deallocate_gridvertex_nbrs(GridVertex(j))
       end do
       deallocate(GridVertex)
    end if

  end subroutine prim_init1_cleanup

  subroutine prim_init1_buffers (elem,par)
    use control_mod,        only : integration, transport_alg
    use edge_mod,           only : initedgebuffer, edge_g
    use parallel_mod,       only : parallel_t
    use prim_advance_mod,   only : prim_advance_init1
    use prim_advection_mod, only : prim_advec_init1
    use thread_mod,         only : hthreads
#ifdef TRILINOS
    use prim_implicit_mod,  only : prim_implicit_init
#endif
#ifdef HOMME_ENABLE_COMPOSE
    use dimensions_mod,     only : max_corner_elem
    use compose_mod,        only : compose_query_bufsz, compose_set_bufs
#endif
    !
    ! Inputs
    !
    type (element_t),   pointer     :: elem(:)
    type (parallel_t),  intent(in)  :: par

    integer :: edgesz, sendsz, recvsz, n, den

    call prim_advance_init1(par,elem,integration)
#ifdef TRILINOS
    call prim_implicit_init(par, elem)
#endif
    call Prim_Advec_Init1(par, elem)

    ! single global edge buffer for all models:
    ! hydrostatic 4*nlev      NH:  6*nlev+1  
    ! SL tracers: (qsize+1)*nlev   e3sm:  (qsize+3)*nlev+2
    ! if this is too small, code will abort with an error message    
    edgesz = max((qsize+3)*nlev+2,6*nlev+1)

#ifdef HOMME_ENABLE_COMPOSE
    if (transport_alg > 0) then
       ! slmm_init_impl has already run, called from compose_init. Now
       ! find out how much memory SLMM wants.
       call compose_query_bufsz(sendsz, recvsz)
       ! from initEdgeBuffer nbuf calc
       den = 4*(np+max_corner_elem)*nelemd
       n = (max(sendsz, recvsz) + den - 1)/den
       edgesz = max(edgesz, n)
    end if
#endif

    call initEdgeBuffer(par,edge_g,elem,edgesz)

#ifdef HOMME_ENABLE_COMPOSE
    if (transport_alg > 0) call compose_set_bufs(edge_g%buf, edge_g%receive)
#endif

  end subroutine prim_init1_buffers

  !_____________________________________________________________________
  subroutine prim_init2(elem, hybrid, nets, nete, tl, hvcoord)

    use control_mod,          only: runtype, test_case, &
                                    debug_level, vfile_int, vform, vfile_mid, &
                                    topology, dt_remap_factor, dt_tracer_factor,&
                                    sub_case, limiter_option, nu, nu_q, nu_div, tstep_type, hypervis_subcycle, &
                                    hypervis_subcycle_q, hypervis_subcycle_tom
    use global_norms_mod,     only: test_global_integral, print_cfl
    use hybvcoord_mod,        only: hvcoord_t
    use parallel_mod,         only: parallel_t, haltmp, syncmp, abortmp
    use prim_state_mod,       only: prim_printstate, prim_diag_scalars
    use prim_advection_mod,   only: prim_advec_init2
    use model_init_mod,       only: model_init2
    use time_mod,             only: timelevel_t, tstep, timelevel_init, nendstep, smooth, nsplit, TimeLevel_Qdp
    use control_mod,          only: smooth_phis_numcycle

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

    if (topology == "cube" .OR. topology=='plane') then
       call test_global_integral(elem, hybrid,nets,nete)
    end if


    ! compute most restrictive dt*nu for use by variable res viscosity:
    ! compute timestep seen by viscosity operator:
    dt_dyn_vis = tstep
    if (dt_tracer_factor>1 .and. tstep_type == 1) then
       ! tstep_type==1: RK2 followed by LF.  internal LF stages apply viscosity at 2*dt
       dt_dyn_vis = 2*tstep
    endif
    dt_tracer_vis=tstep*dt_tracer_factor
    
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

#if !defined(CAM) && !defined(SCREAM)

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

       call set_test_initial_conditions(elem,deriv1,hybrid,hvcoord,tl,nets,nete)

       call ReadRestart(elem,hybrid%ithr,nets,nete,tl)

       if (runtype==2) then
          do ie=nets,nete
             call copy_state(elem(ie),tl%n0,tl%nm1)
          enddo
       endif ! runtype==2

    else

      ! runtype=0: initial run
      if (hybrid%masterthread) write(iulog,*) ' runtype: initial run'
      call set_test_initial_conditions(elem,deriv1,hybrid,hvcoord,tl,nets,nete)
      if (hybrid%masterthread) write(iulog,*) '...done'

!      do ie=nets,nete
!        ! set perlim in ctl_nl namelist for temperature field initial perturbation
!        elem(ie)%state%T=elem(ie)%state%T * (1.0_real_kind + pertlim)
!      enddo

    endif !runtype

#endif
!$OMP MASTER
    if (runtype==2) then
       ! branch run
       ! reset time counters to zero since timestep may have changed
       nEndStep = nEndStep-tl%nstep ! used by standalone HOMME.  restart code set this to nmax + tl%nstep
       tl%nstep=0
    endif
    tl%nstep0=tl%nstep+1       ! compute diagnostics after 1st step
!$OMP END MASTER
!$OMP BARRIER

#ifdef CAM
    ! initialize dp3d from ps_v.  CAM IC/restart code reads ps_v, doesn't
    ! have access to hvcoord to compute dp3d:
    do ie=nets,nete
       do k=1,nlev
          elem(ie)%state%dp3d(:,:,k,tl%n0)=&
               ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
       enddo
    end do
#endif


    ! For new runs, and branch runs, convert state variable Q to (Qdp)
    ! because initial conditon reads in Q, not Qdp
    ! restart runs will read dpQ from restart file
    if (runtype==0 .or. runtype==2) then
       do ie=nets,nete
          elem(ie)%derived%omega_p(:,:,:) = 0D0
       end do
       do ie=nets,nete
          do q=1,qsize
             elem(ie)%state%Qdp(:,:,:,q,1)=elem(ie)%state%Q(:,:,:,q)*elem(ie)%state%dp3d(:,:,:,tl%n0)
             elem(ie)%state%Qdp(:,:,:,q,2)=elem(ie)%state%Q(:,:,:,q)*elem(ie)%state%dp3d(:,:,:,tl%n0)
          enddo
       enddo
    endif

    if (runtype==1) then
       call TimeLevel_Qdp( tl, dt_tracer_factor, n0_qdp)
       do ie=nets,nete
          do q=1,qsize
             elem(ie)%state%Q(:,:,:,q)=elem(ie)%state%Qdp(:,:,:,q,n0_qdp)/elem(ie)%state%dp3d(:,:,:,tl%n0)
          enddo
       enddo
    endif

    call model_init2(elem(:), hybrid,deriv1,hvcoord,tl,nets,nete)

    ! advective and viscious CFL estimates
    ! may also adjust tensor coefficients based on CFL
    call print_cfl(elem,hybrid,nets,nete,dtnu)

    if (hybrid%masterthread) then
       ! CAM has set tstep based on dtime before calling prim_init2(),
       ! so only now does HOMME learn the timstep.  print them out:
       write(iulog,'(a,2f9.2)') "dt_remap: (0=disabled)   ",tstep*dt_remap_factor
       if (qsize>0) then
          write(iulog,'(a,2f9.2)') "dt_tracer (SE), per RK stage: ", &
               tstep*dt_tracer_factor,(tstep*dt_tracer_factor)/2
       end if
       write(iulog,'(a,2f9.2)')    "dt_dyn:                  ",tstep
       write(iulog,'(a,2f9.2)')    "dt_dyn (viscosity):      ",dt_dyn_vis
       write(iulog,'(a,2f9.2)')    "dt_tracer (viscosity):   ",dt_tracer_vis
       if (hypervis_subcycle_tom==0) then                                                     
          ! applied with hyperviscosity                                                       
          write(iulog,'(a,2f9.2)') "dt_vis_TOM:  ",dt_dyn_vis                                 
       else                                                                                   
          write(iulog,'(a,2f9.2)') "dt_vis_TOM:  ",tstep/hypervis_subcycle_tom               
       endif                                                                 


#ifdef CAM
       write(iulog,'(a,2f9.2)') "CAM dtime (dt_phys):         ",tstep*nsplit*max(dt_remap_factor, dt_tracer_factor)
#endif
    end if

    if (hybrid%masterthread) write(iulog,*) "initial state:"
    call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)
    call Prim_Advec_Init2(elem(:), hvcoord, hybrid)

  end subroutine prim_init2

!=======================================================================================================!



  subroutine prim_run_subcycle(elem, hybrid,nets,nete, dt, single_column, tl, hvcoord,nsubstep)

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

    use control_mod,        only: statefreq, qsplit, rsplit, disable_diagnostics, &
         dt_remap_factor, dt_tracer_factor, transport_alg
    use hybvcoord_mod,      only: hvcoord_t
    use parallel_mod,       only: abortmp
    use prim_state_mod,     only: prim_printstate
    use vertremap_mod,      only: vertical_remap
    use reduction_mod,      only: parallelmax
    use time_mod,           only: TimeLevel_t, timelevel_update, timelevel_qdp, nsplit, tstep
#if USE_OPENACC
    use openacc_utils_mod,  only: copy_qdp_h2d, copy_qdp_d2h
#endif

    implicit none

    type (element_t) ,    intent(inout) :: elem(:)
    type (hybrid_t),      intent(in)    :: hybrid                       ! distributed parallel structure (shared)
    type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets                         ! starting thread element number (private)
    integer,              intent(in)    :: nete                         ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt                           ! "timestep dependent" timestep
    logical,              intent(in)    :: single_column
    type (TimeLevel_t),   intent(inout) :: tl
    integer,              intent(in)    :: nsubstep                     ! nsubstep = 1 .. nsplit

    real(kind=real_kind) :: dp, dt_q, dt_remap
    real(kind=real_kind) :: dp_np1(np,np)
    integer :: ie,i,j,k,n,q,t,scm_dum
    integer :: n0_qdp,np1_qdp,r,nstep_end,nets_in,nete_in,step_factor
    logical :: compute_diagnostics, independent_time_steps

    ! Use the flexible time stepper if dt_remap_factor == 0 (vertically Eulerian
    ! dynamics) or dt_remap < dt_tracer. This applies to SL transport only.
    independent_time_steps = transport_alg > 1 .and. dt_remap_factor < dt_tracer_factor

    ! compute timesteps for tracer transport and vertical remap
    dt_q = dt*dt_tracer_factor
    if (dt_remap_factor == 0) then
       dt_remap  = dt_q
       nstep_end = tl%nstep + dt_tracer_factor
    else
       ! dt_remap_factor = 0 means use eulerian code, not vert. lagrange
       dt_remap  = dt*dt_remap_factor
       step_factor = max(dt_remap_factor, dt_tracer_factor)
       nstep_end = tl%nstep + step_factor ! nstep at end of this routine
    endif

    ! activate diagnostics periodically for display to stdout and on first 2 timesteps
    compute_diagnostics   = .false.
    if (MODULO(nstep_end,statefreq)==0 .or. (tl%nstep <= tl%nstep0+(nstep_end-tl%nstep) )) then
       compute_diagnostics= .true.
    endif
    if(disable_diagnostics) compute_diagnostics= .false.

    ! compute scalar diagnostics if currently active
    if (compute_diagnostics) call run_diagnostics(elem,hvcoord,tl,3,.true.,nets,nete)

    if (.not. independent_time_steps) then
       call TimeLevel_Qdp(tl, dt_tracer_factor, n0_qdp, np1_qdp)

#ifndef CAM
       ! compute HOMME test case forcing
       ! by calling it here, it mimics eam forcings computations in standalone
       ! homme.
       call compute_test_forcing(elem,hybrid,hvcoord,tl%n0,n0_qdp,dt_remap,nets,nete,tl)
#endif

       call applyCAMforcing_remap(elem,hvcoord,tl%n0,n0_qdp,dt_remap,nets,nete)

       ! E(1) Energy after CAM forcing
       if (compute_diagnostics) call run_diagnostics(elem,hvcoord,tl,1,.true.,nets,nete)

#if (USE_OPENACC)
       !    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
       call t_startf("copy_qdp_h2d")
       call copy_qdp_h2d( elem , n0_qdp )
       call t_stopf("copy_qdp_h2d")
#endif

      if (.not. single_column) then 

        ! Loop over rsplit vertically lagrangian timesiteps
        call prim_step(elem, hybrid, nets, nete, dt, tl, hvcoord, compute_diagnostics)

        do r=2,rsplit
          call TimeLevel_update(tl,"leapfrog")
          call prim_step(elem, hybrid, nets, nete, dt, tl, hvcoord, .false.)
        enddo

      else 

        ! Single Column Case
        ! Loop over rsplit vertically lagrangian timesiteps
        call prim_step_scm(elem, nets, nete, dt, tl, hvcoord)
        do r=2,rsplit
          call TimeLevel_update(tl,"leapfrog")
          call prim_step_scm(elem, nets, nete, dt, tl, hvcoord)
        enddo

      endif

      ! defer final timelevel update until after remap and diagnostics
      !compute timelevels for tracers (no longer the same as dynamics)
      call TimeLevel_Qdp( tl, dt_tracer_factor, n0_qdp, np1_qdp)

#if (USE_OPENACC)
      call t_startf("copy_qdp_h2d")
      call copy_qdp_d2h( elem , np1_qdp )
      call t_stopf("copy_qdp_h2d")
#endif

      if (compute_diagnostics) call run_diagnostics(elem,hvcoord,tl,4,.false.,nets,nete)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  apply vertical remap
      !  always for tracers
      !  if dt_remap_factor>0:  also remap dynamics back to reference levels.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (single_column) then
        nets_in=1
        nete_in=1
      else
        nets_in=nets
        nete_in=nete
      endif

      call vertical_remap(hybrid,elem,hvcoord,dt_remap,tl%np1,np1_qdp,nets_in,nete_in)
    else
      ! This time stepping routine permits the vertical remap time
      ! step to be shorter than the tracer transport time step.
      call prim_step_flexible(hybrid, elem, nets, nete, dt, tl, hvcoord, compute_diagnostics)
    end if ! independent_time_steps

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! time step is complete.  update some diagnostic variables:
    ! Q    (mixing ratio)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - 2*dt
    !   u(n0)    dynamics at  t+dt_remap - dt
    !   u(np1)   dynamics at  t+dt_remap
    !
    !   Q(1)   Q at t+dt_remap
    if (compute_diagnostics) call run_diagnostics(elem,hvcoord,tl,2,.false.,nets,nete)
    
    ! =================================
    ! update dynamics time level pointers
    ! =================================
    call TimeLevel_update(tl,"leapfrog")
    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - dt       
    !   u(n0)    dynamics at  t+dt_remap
    !   u(np1)   undefined

    ! ============================================================
    ! Print some diagnostic information
    ! ============================================================
    if (compute_diagnostics) then
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)
    end if
  end subroutine prim_run_subcycle



  subroutine prim_step(elem, hybrid,nets,nete, dt, tl, hvcoord, compute_diagnostics)
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
    use control_mod,        only: statefreq, integration, ftype, nu_p, dt_tracer_factor, dt_remap_factor
    use control_mod,        only: transport_alg
    use hybvcoord_mod,      only: hvcoord_t
    use parallel_mod,       only: abortmp
    use prim_advance_mod,   only: prim_advance_exp, applycamforcing_dynamics
    use prim_advection_mod, only: prim_advec_tracers_remap
    use reduction_mod,      only: parallelmax
    use time_mod,           only: time_at,TimeLevel_t, timelevel_update, nsplit
    use prim_state_mod,     only: prim_printstate

    type(element_t),      intent(inout) :: elem(:)
    type(hybrid_t),       intent(in)    :: hybrid   ! distributed parallel structure (shared)
    type(hvcoord_t),      intent(in)    :: hvcoord  ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets     ! starting thread element number (private)
    integer,              intent(in)    :: nete     ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt       ! "timestep dependent" timestep
    type(TimeLevel_t),    intent(inout) :: tl

    real(kind=real_kind) :: st, st1, dp, dt_q
    integer :: ie, t, q,k,i,j,n
    real (kind=real_kind)                          :: maxcflx, maxcfly
    real (kind=real_kind) :: dp_np1(np,np)
    logical :: compute_diagnostics

    dt_q = dt*dt_tracer_factor

    call set_tracer_transport_derived_values(elem, nets, nete, tl)
 
    ! ===============
    ! Dynamical Step
    ! for ftype==4, also apply dynamics tendencies from forcing
    ! for ftype==4, energy diagnostics will be incorrect
    ! ===============
    if (ftype==4) then
       call ApplyCAMforcing_dynamics(elem,hvcoord,tl%n0,dt,nets,nete)
       ! E(1) Energy after CAM forcing applied
       ! with ftype==4, need (E(1)-E(3))/dt_dyn instead (E(1)-E(3))/dt_tracer
       if (compute_diagnostics) call run_diagnostics(elem,hvcoord,tl,1,.true.,nets,nete)
    endif
       
    call prim_advance_exp(elem,deriv1,hvcoord,hybrid,dt,tl,nets,nete,compute_diagnostics)
    do n=2,dt_tracer_factor
       call TimeLevel_update(tl,"leapfrog")
       if (ftype==4) call ApplyCAMforcing_dynamics(elem,hvcoord,tl%n0,dt,nets,nete)
       call prim_advance_exp(elem, deriv1, hvcoord,hybrid, dt, tl, nets, nete, .false.)
       ! defer final timelevel update until after Q update.
    enddo


    ! current dynamics state variables:
    !    derived%dp              =  dp at start of timestep
    !    derived%vstar           =  velocity at start of tracer timestep
    !    derived%vn0             =  mean horiz. flux:   U*dp
    !    state%dp3d(:,:,:,np1)   = dp3d
    ! dt_remap_factor=0
    !        state%v(:,:,:,np1)      = velocity on reference levels
    ! dt_remap_factor>0
    !        state%v(:,:,:,np1)      = velocity on lagrangian levels 
    !        
    ! Tracer Advection.  
    ! in addition, this routine will apply the DSS to:
    !        derived%eta_dot_dpdn    =  mean vertical velocity (used for remap below)
    !        derived%omega           =
    ! Tracers are always vertically lagrangian.  
    ! For dt_remap_factor=0: 
    !   if tracer scheme needs v on lagrangian levels it has to vertically interpolate

    call t_startf("prim_step_advec")
    if (qsize > 0) then
      call t_startf("PAT_remap")
      call Prim_Advec_Tracers_remap(elem, deriv1,hvcoord,hybrid,dt_q,tl,nets,nete)
      call t_stopf("PAT_remap")
    end if
    call t_stopf("prim_step_advec")

  end subroutine prim_step

  subroutine prim_step_flexible(hybrid, elem, nets, nete, dt, tl, hvcoord, compute_diagnostics)
    use control_mod,        only: ftype, nu_p, dt_tracer_factor, dt_remap_factor, prescribed_wind, transport_alg
    use hybvcoord_mod,      only: hvcoord_t
    use parallel_mod,       only: abortmp
    use prim_advance_mod,   only: prim_advance_exp, applycamforcing_dynamics
    use prim_advection_mod, only: prim_advec_tracers_remap
    use reduction_mod,      only: parallelmax
    use time_mod,           only: TimeLevel_t, timelevel_update, timelevel_qdp
    use prim_state_mod,     only: prim_printstate
    use vertremap_mod,      only: vertical_remap
    use sl_advection,       only: sl_vertically_remap_tracers

    type(element_t),      intent(inout) :: elem(:)
    type(hybrid_t),       intent(in)    :: hybrid   ! distributed parallel structure (shared)
    type(hvcoord_t),      intent(in)    :: hvcoord  ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets     ! starting thread element number (private)
    integer,              intent(in)    :: nete     ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt       ! "timestep dependent" timestep
    type(TimeLevel_t),    intent(inout) :: tl
    logical,              intent(in)    :: compute_diagnostics

    real(kind=real_kind) :: dt_q, dt_remap, dp(np,np,nlev)
    integer :: ie, q, k, n, n0_qdp, np1_qdp
    logical :: compute_diagnostics_it, apply_forcing

    dt_q = dt*dt_tracer_factor
    if (dt_remap_factor == 0) then
       dt_remap = dt
    else
       dt_remap = dt*dt_remap_factor
    end if

    call TimeLevel_Qdp(tl, dt_tracer_factor, n0_qdp, np1_qdp)

#ifndef CAM
    ! Compute test forcing over tracer time step.
    call compute_test_forcing(elem,hybrid,hvcoord,tl%n0,n0_qdp,dt_q,nets,nete,tl)
#endif

#ifdef CAM
    apply_forcing = ftype == 0
#else
    apply_forcing = ftype == 0 .or. ftype == 2 .or. ftype == 4
#endif
    if (apply_forcing) then
       ! Apply tracer forcings over tracer time step.
       do ie = nets,nete
          call ApplyCAMForcing_tracers(elem(ie),hvcoord,tl%n0,n0_qdp,dt_q,.false.)
       enddo
    end if

    call set_tracer_transport_derived_values(elem, nets, nete, tl)

    call t_startf("prim_step_dyn")
    do n = 1, dt_tracer_factor
       compute_diagnostics_it = logical(compute_diagnostics .and. n == 1)

       if (n > 1) call TimeLevel_update(tl, "leapfrog")

       if (ftype == 4) then
          ! also apply dynamics tendencies from forcing; energy
          ! diagnostics will be incorrect
          call ApplyCAMforcing_dynamics(elem,hvcoord,tl%n0,dt,nets,nete)
          if (compute_diagnostics_it) call run_diagnostics(elem,hvcoord,tl,1,.true.,nets,nete)
       else if (ftype == 2 .or. ftype == 0) then
          ! Apply dynamics forcing over the dynamics (vertically Eulerian) or
          ! vertical remap time step if we're at reference levels.
          if (dt_remap_factor > 0) then
             apply_forcing = modulo(n-1, dt_remap_factor) == 0
          else
             apply_forcing = .true.
          end if
          if (apply_forcing) then
             call ApplyCAMforcing_dynamics(elem,hvcoord,tl%n0,dt_remap,nets,nete)
             if (compute_diagnostics_it) call run_diagnostics(elem,hvcoord,tl,1,.true.,nets,nete)
          end if
       end if

       call prim_advance_exp(elem, deriv1, hvcoord, hybrid, dt, tl, nets, nete, &
            compute_diagnostics_it)

       if (dt_remap_factor == 0) then
          ! Set np1_qdp to -1. Since dt_remap == 0, the only part of
          ! vertical_remap that is active is the updates to
          ! ps_v(:,:,np1) and dp3d(:,:,:,np1).
          call vertical_remap(hybrid, elem, hvcoord, dt_remap, tl%np1, -1, nets, nete)
       else
          if (modulo(n, dt_remap_factor) == 0) then
             if (compute_diagnostics) call run_diagnostics(elem,hvcoord,tl,4,.false.,nets,nete)
             if (prescribed_wind == 1) then
                ! Prescribed winds are evaluated on reference levels,
                ! not floating levels, so don't remap, just update dp3d.
                do ie = nets,nete
                   elem(ie)%state%ps_v(:,:,tl%np1) = hvcoord%hyai(1)*hvcoord%ps0 + &
                        sum(elem(ie)%state%dp3d(:,:,:,tl%np1),3)
                   do k=1,nlev
                      dp(:,:,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                                  (hvcoord%hybi(k+1) - hvcoord%hybi(k))*elem(ie)%state%ps_v(:,:,tl%np1)
                   end do
                   elem(ie)%state%dp3d(:,:,:,tl%np1) = dp
                end do
             else
                ! Set np1_qdp to -1 to remap dynamics variables but
                ! not tracers.
                call vertical_remap(hybrid, elem, hvcoord, dt_remap, tl%np1, -1, nets, nete)
             end if
          end if
       end if
       ! defer final timelevel update until after Q update.
    enddo
    call t_stopf("prim_step_dyn")

    if (qsize > 0) then
       call t_startf("PAT_remap")
       call Prim_Advec_Tracers_remap(elem, deriv1, hvcoord, hybrid, dt_q, tl, nets, nete)
       call t_stopf("PAT_remap")
    end if

    if (dt_remap_factor == 0 .and. compute_diagnostics) then
       call TimeLevel_Qdp(tl, dt_tracer_factor, n0_qdp, np1_qdp)
       call run_diagnostics(elem,hvcoord,tl,4,.false.,nets,nete)
    end if

    ! Remap tracers.
    if (qsize > 0) then
       call sl_vertically_remap_tracers(hybrid, elem, nets, nete, tl, dt_q)
    end if
  end subroutine prim_step_flexible

  subroutine run_diagnostics(elem, hvcoord, tl, n, t_before_advance, nets, nete)
    use time_mod,           only: TimeLevel_t
    use hybvcoord_mod,      only: hvcoord_t
    use prim_state_mod,     only: prim_diag_scalars, prim_energy_halftimes

    type(element_t),      intent(inout) :: elem(:)
    type(TimeLevel_t),    intent(in)    :: tl
    type (hvcoord_t),     intent(in)    :: hvcoord
    integer,              intent(in)    :: nets, nete, n
    logical,              intent(in)    :: t_before_advance

    call t_startf("prim_diag")
    call prim_diag_scalars(elem, hvcoord, tl, n, t_before_advance, nets, nete)
    call prim_energy_halftimes(elem, hvcoord, tl, n, t_before_advance, nets, nete)
    call t_stopf("prim_diag")
  end subroutine run_diagnostics

  subroutine set_tracer_transport_derived_values(elem, nets, nete, tl)
    use control_mod,        only: nu_p, transport_alg
    use time_mod,           only: TimeLevel_t

    type(element_t),      intent(inout) :: elem(:)
    integer,              intent(in)    :: nets, nete
    type(TimeLevel_t),    intent(in)    :: tl

    integer :: ie

    ! ===============
    ! initialize mean flux accumulation variables and save some variables at n0
    ! for use by advection
    ! ===============
    do ie = nets,nete
       elem(ie)%derived%eta_dot_dpdn = 0     ! mean vertical mass flux
       elem(ie)%derived%vn0 = 0              ! mean horizontal mass flux
       elem(ie)%derived%omega_p = 0
       if (nu_p > 0) then
          elem(ie)%derived%dpdiss_ave = 0
          elem(ie)%derived%dpdiss_biharmonic = 0
       endif
       if (transport_alg > 0) then
          elem(ie)%derived%vstar = elem(ie)%state%v(:,:,:,:,tl%n0)
       end if
       elem(ie)%derived%dp(:,:,:) = elem(ie)%state%dp3d(:,:,:,tl%n0)
    enddo
  end subroutine set_tracer_transport_derived_values

!---------------------------------------------------------------------------
!
! Apply all forcing terms that are applied with frequency dt_remap 
!
! Note on ftypes:
!   ftype= 4: Q was adjusted by physics, dynamics tendencies applied elsewhere
!   ftype= 2: Q was adjusted by physics, but apply u,T forcing here
!   ftype= 1: forcing was applied time-split in CAM coupling layer
!   ftype= 0: apply all forcing here
!   ftype=-1: do not apply forcing
  subroutine applyCAMforcing_remap(elem,hvcoord,n0,n0qdp,dt_remap,nets,nete)
  use control_mod,        only : ftype
  use hybvcoord_mod,      only : hvcoord_t
  use prim_advance_mod,   only : applycamforcing_dynamics
  implicit none
  type (element_t),       intent(inout) :: elem(:)
  real (kind=real_kind),  intent(in)    :: dt_remap
  type (hvcoord_t),       intent(in)    :: hvcoord
  integer,                intent(in)    :: n0,n0qdp,nets,nete
  integer                               :: ie

  call t_startf("ApplyCAMForcing_remap")
  if (ftype==-1) then
    !do nothing
  elseif (ftype==0) then
    do ie = nets,nete
       call applyCAMforcing_tracers (elem(ie),hvcoord,n0,n0qdp,dt_remap,.false.)
    enddo
    call applyCAMforcing_dynamics(elem,hvcoord,n0,dt_remap,nets,nete)
  elseif (ftype==1) then
    !do nothing
  elseif (ftype==2) then
    ! with CAM physics, tracers were adjusted in dp coupling layer
#ifndef CAM
    do ie = nets,nete
       call ApplyCAMForcing_tracers (elem(ie),hvcoord,n0,n0qdp,dt_remap,.false.)
    enddo
#endif
    call ApplyCAMForcing_dynamics(elem,hvcoord,n0,dt_remap,nets,nete)
 elseif (ftype==4) then
    ! with CAM physics, tracers were adjusted in dp coupling layer
#ifndef CAM
    do ie = nets,nete
       call ApplyCAMForcing_tracers (elem(ie),hvcoord,n0,n0qdp,dt_remap,.false.)
    enddo
#endif
  endif

  call t_stopf("ApplyCAMForcing_remap")
  end subroutine applyCAMforcing_remap



  subroutine applyCAMforcing_tracers(elem,hvcoord,np1,np1_qdp,dt,adjustment)
  !
  ! Apply forcing to tracers
  !    adjustment=1:  apply forcing as hard adjustment, assume qneg check already done
  !    adjustment=0:  apply tracer tendency
  ! in both cases, update PS to conserve mass
  !
  ! For theta model, convert temperature tendency to theta/phi tendency
  ! this conversion is done assuming constant pressure except for changes to hydrostatic
  ! pressure from the water vapor tendencies. It is thus recomputed whenever
  ! water vapor tendency is applied
  ! 
  ! theta model hydrostatic requires this constant pressure assumption due to 
  ! phi/density being diagnostic.  theta model NH could do the conversion constant 
  ! density which would simplify this routine
  !
  ! NOTE about ps_v/dp3d
  ! init:
  !   (both ps_v and dp3d are valid)
  ! do: 
  !    physics  (uses ps_v to compute pressure levels. doesn't change ps_v)
  !    applyCAMforcing_tracers  use ps_v for initial pressure.  
  !                             may adjust dp3d for mass conservation (if adjust_ps=.false.)
  !                             ps_v no longer valid
  !    dynamics                 should only use dp3d
  !    remap                    remap back to ref levels.  ps_v now valid
  !    write restart files      ps_v ok for restart
  !
  use control_mod,        only : use_moisture, dt_remap_factor
  use hybvcoord_mod,      only : hvcoord_t
#ifdef MODEL_THETA_L
  use control_mod,        only : theta_hydrostatic_mode
  use physical_constants, only : cp, g, kappa, Rgas, p0
  use element_ops,        only : get_temperature, get_r_star, get_hydro_pressure
  use eos,                only : pnh_and_exner_from_eos
#ifdef HOMMEXX_BFB_TESTING
  use bfb_mod,            only : bfb_pow
#endif
#endif
  implicit none
  type (element_t),       intent(inout) :: elem
  real (kind=real_kind),  intent(in)    :: dt
  type (hvcoord_t),       intent(in)    :: hvcoord
  integer,                intent(in)    :: np1,np1_qdp
  logical,                intent(in)    :: adjustment

  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind)  :: fq
  real (kind=real_kind)  :: dp(np,np,nlev), ps(np,np), dp_adj(np,np,nlev)
  real (kind=real_kind)  :: phydro(np,np,nlev)  ! hydrostatic pressure
  logical :: adjust_ps   ! adjust PS or DP3D to conserve dry mass
#ifdef MODEL_THETA_L
  real (kind=real_kind)  :: pprime(np,np,nlev)
  real (kind=real_kind)  :: vthn1(np,np,nlev)
  real (kind=real_kind)  :: tn1(np,np,nlev)
  real (kind=real_kind)  :: pnh(np,np,nlev)
  real (kind=real_kind)  :: phi_n1(np,np,nlevp)
  real (kind=real_kind)  :: rstarn1(np,np,nlev)
  real (kind=real_kind)  :: exner(np,np,nlev)
  real (kind=real_kind)  :: dpnh_dp_i(np,np,nlevp)
#endif

#ifdef MODEL_THETA_L
  if (dt_remap_factor==0) then
     adjust_ps=.true.   ! stay on reference levels for Eulerian case
  else
     adjust_ps=.true.   ! Lagrangian case can support adjusting dp3d or ps
  endif
#else
  adjust_ps=.true.      ! preqx requires forcing to stay on reference levels
#endif

  dp=elem%state%dp3d(:,:,:,np1)
  dp_adj=dp
  ps=elem%state%ps_v(:,:,np1)
  !ps=hvcoord%hyai(1)*hvcoord%ps0 + sum(dp(:,:,:),3) ! introduces roundoff

  ! after calling this routine, ps_v may not be valid and should not be used
  elem%state%ps_v(:,:,np1)=0


#ifdef MODEL_THETA_L
   !compute temperatue and NH perturbation pressure before Q tendency
   do k=1,nlev
      phydro(:,:,k)=hvcoord%ps0*hvcoord%hyam(k) + ps(:,:)*hvcoord%hybm(k)
   enddo

   !one can set pprime=0 to hydro regime but it is not done in master
   !compute pnh, here only pnh is needed
   call pnh_and_exner_from_eos(hvcoord,elem%state%vtheta_dp(:,:,:,np1),dp,&
        elem%state%phinh_i(:,:,:,np1),pnh,exner,dpnh_dp_i)
   do k=1,nlev
      pprime(:,:,k) = pnh(:,:,k)-phydro(:,:,k)
   enddo
   call get_R_star(rstarn1,elem%state%Q(:,:,:,1))
   tn1=exner* elem%state%vtheta_dp(:,:,:,np1)*(Rgas/rstarn1) / dp
#endif

   if (adjustment) then 
      ! hard adjust Q from physics.  negativity check done in physics
      do k=1,nlev
         do j=1,np
            do i=1,np
               do q=1,qsize
                  ! apply forcing to Qdp
                  ! dyn_in%elem(ie)%state%Qdp(i,j,k,q,tl_fQdp) = &
                  !        dyn_in%elem(ie)%state%Qdp(i,j,k,q,tl_fQdp) + fq 
                  elem%state%Qdp(i,j,k,q,np1_qdp) = &
                       dp(i,j,k)*elem%derived%FQ(i,j,k,q)
                  
                  if (q==1) then
                     fq = dp(i,j,k)*( elem%derived%FQ(i,j,k,q) -&
                          elem%state%Q(i,j,k,q))
                     ! force ps to conserve mass:  
                     ps(i,j)=ps(i,j) + fq
                     dp_adj(i,j,k)=dp_adj(i,j,k) + fq   !  ps =  ps0+sum(dp(k))
                  endif
               enddo
            end do
         end do
      end do
   else ! end of adjustment
      ! apply forcing to Qdp
      elem%derived%FQps(:,:)=0
      do q=1,qsize
         do k=1,nlev
            do j=1,np
               do i=1,np
                  fq = dt*elem%derived%FQ(i,j,k,q)
                  if (elem%state%Qdp(i,j,k,q,np1_qdp) + fq < 0 .and. fq<0) then
                     if (elem%state%Qdp(i,j,k,q,np1_qdp) < 0 ) then
                        fq=0  ! Q already negative, dont make it more so
                     else
                        fq = -elem%state%Qdp(i,j,k,q,np1_qdp)
                     endif
                  endif
                  elem%state%Qdp(i,j,k,q,np1_qdp) = elem%state%Qdp(i,j,k,q,np1_qdp)+fq
                  if (q==1) then
                     elem%derived%FQps(i,j)=elem%derived%FQps(i,j)+fq/dt
                     dp_adj(i,j,k)=dp_adj(i,j,k) + fq
                  endif
               enddo
            enddo
         enddo
      enddo

      ! to conserve dry mass in the precese of Q1 forcing:
      ps(:,:) = ps(:,:) + dt*elem%derived%FQps(:,:)
   endif ! if adjustment


   if (use_moisture) then
      ! compute water vapor adjusted dp3d:
      if (adjust_ps) then
         ! compute new dp3d from adjusted ps()
         do k=1,nlev
            dp_adj(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                 ( hvcoord%hybi(k+1) - hvcoord%hybi(k))*ps(:,:)
         enddo
      endif
      elem%state%dp3d(:,:,:,np1)=dp_adj(:,:,:)
   endif

   ! Qdp(np1) was updated by forcing - update Q(np1)
   do q=1,qsize
      elem%state%Q(:,:,:,q) = elem%state%Qdp(:,:,:,q,np1_qdp)/elem%state%dp3d(:,:,:,np1)
   enddo
   

#ifdef MODEL_THETA_L
   if (use_moisture) then
      ! compute updated pnh and exner
      if (adjust_ps) then
         ! recompute hydrostatic pressure from ps
         do k=1,nlev  
            phydro(:,:,k)=hvcoord%ps0*hvcoord%hyam(k) + ps(:,:)*hvcoord%hybm(k)
         enddo
      else
         ! recompute hydrostatic pressure from dp3d
         call get_hydro_pressure(phydro,elem%state%dp3d(:,:,:,np1),hvcoord)
      endif
      do k=1,nlev
         pnh(:,:,k)=phydro(:,:,k) + pprime(:,:,k)
#ifdef HOMMEXX_BFB_TESTING
         exner(:,:,k)=bfb_pow(pnh(:,:,k)/p0,Rgas/Cp)
#else
         exner(:,:,k)=(pnh(:,:,k)/p0)**(Rgas/Cp)
#endif
      enddo
   endif
   
   !update temperature
   call get_R_star(rstarn1,elem%state%Q(:,:,:,1))
   tn1(:,:,:) = tn1(:,:,:) + dt*elem%derived%FT(:,:,:)
   
   
   ! now we have tn1,dp,pnh - compute corresponding theta and phi:
   vthn1 =  (rstarn1(:,:,:)/Rgas)*tn1(:,:,:)*elem%state%dp3d(:,:,:,np1)/exner(:,:,:)
     
   phi_n1(:,:,nlevp)=elem%state%phinh_i(:,:,nlevp,np1)
   do k=nlev,1,-1
      phi_n1(:,:,k)=phi_n1(:,:,k+1) + Rgas*vthn1(:,:,k)*exner(:,:,k)/pnh(:,:,k)
   enddo
   
   !finally, compute difference for FVTheta
   ! this method is using new dp, new exner, new-new r*, new t
   elem%derived%FVTheta(:,:,:) = &
        (vthn1 - elem%state%vtheta_dp(:,:,:,np1))/dt
   
   elem%derived%FPHI(:,:,:) = &
        (phi_n1 - elem%state%phinh_i(:,:,:,np1))/dt
   
#endif

  end subroutine applyCAMforcing_tracers
  
  
  subroutine prim_step_scm(elem, nets,nete, dt, tl, hvcoord)
  !
  !   prim_step version for single column model (SCM)
  !   Here we simply want to compute the floating level tendency
  !    based on the prescribed large scale vertical velocity
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
    use control_mod,        only: statefreq, integration, ftype, nu_p, dt_tracer_factor, dt_remap_factor
    use control_mod,        only: transport_alg
    use hybvcoord_mod,      only : hvcoord_t
    use parallel_mod,       only: abortmp
    use prim_advance_mod,   only: prim_advance_exp, applyCAMforcing_dynamics
    use reduction_mod,      only: parallelmax
    use time_mod,           only: time_at,TimeLevel_t, timelevel_update, timelevel_qdp, nsplit
    use prim_state_mod,     only: prim_printstate, prim_diag_scalars, prim_energy_halftimes

    type(element_t),      intent(inout) :: elem(:)
    type(hvcoord_t),      intent(in)    :: hvcoord  ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets     ! starting thread element number (private)
    integer,              intent(in)    :: nete     ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt       ! "timestep dependent" timestep
    type(TimeLevel_t),    intent(inout) :: tl

    real(kind=real_kind) :: st, st1, dp, dt_q
    integer :: ie, t, q,k,i,j,n,qn0
    real (kind=real_kind)                          :: maxcflx, maxcfly
    real (kind=real_kind) :: dp_np1(np,np)
 
    ! ===============
    ! initialize mean flux accumulation variables and save some variables at n0
    ! for use by advection
    ! ===============
    do ie=nets,nete
      elem(ie)%derived%eta_dot_dpdn=0     ! mean vertical mass flux
      elem(ie)%derived%vn0=0              ! mean horizontal mass flux
      if (nu_p>0) then
         elem(ie)%derived%dpdiss_ave=0
         elem(ie)%derived%dpdiss_biharmonic=0
      endif
      if (transport_alg > 0) then
        elem(ie)%derived%vstar=elem(ie)%state%v(:,:,:,:,tl%n0)
      end if
      elem(ie)%derived%dp(:,:,:)=elem(ie)%state%dp3d(:,:,:,tl%n0)
    enddo

    ! ===============
    ! Dynamical Step
    ! ===============
    
    call TimeLevel_Qdp(tl, dt_tracer_factor, qn0)  ! compute current Qdp() timelevel 
    call set_prescribed_scm(elem,dt,tl)
    
    do n=2,dt_tracer_factor
 
      call TimeLevel_update(tl,"leapfrog")
      if (ftype==4) call ApplyCAMforcing_dynamics(elem,hvcoord,tl%n0,dt,nets,nete)       

      ! get timelevel for accessing tracer mass Qdp() to compute virtual temperature      
      call TimeLevel_Qdp(tl, dt_tracer_factor, qn0)  ! compute current Qdp() timelevel      
      
      ! call the single column forcing
      call set_prescribed_scm(elem,dt,tl)
      
    enddo

  end subroutine prim_step_scm


!=======================================================================================================!


  subroutine prim_finalize()
#ifdef TRILINOS
  interface
    subroutine noxfinish() bind(C,name='noxfinish')
    use ,intrinsic :: iso_c_binding
    end subroutine noxfinish
  end interface
#endif

#ifdef HOMME_ENABLE_COMPOSE
    use compose_mod, only: compose_finalize
    use control_mod, only: transport_alg
#endif
    implicit none

#ifdef TRILINOS
    call noxfinish()
#endif

#ifdef HOMME_ENABLE_COMPOSE
    if (transport_alg > 0) call compose_finalize()
#endif

    ! ==========================
    ! end of the hybrid program
    ! ==========================
  end subroutine prim_finalize



    subroutine smooth_topo_datasets(elem,hybrid,nets,nete)
    use control_mod, only : smooth_phis_numcycle, smooth_phis_nudt
    use hybrid_mod, only : hybrid_t
    use bndry_mod, only : bndry_exchangev
    use derivative_mod, only : derivative_t , laplace_sphere_wk
    use viscosity_mod, only : smooth_phis
    implicit none

    integer , intent(in) :: nets,nete
    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    ! local
    integer :: ie
    real (kind=real_kind) :: minf
    real (kind=real_kind) :: phis(np,np,nets:nete)

    do ie=nets,nete
       phis(:,:,ie)=elem(ie)%state%phis(:,:)
    enddo
    
    minf=-9e9
    if (hybrid%masterthread) then
       write(iulog,*) "Applying hyperviscosity smoother to PHIS"
       write(iulog,'(a,i10)')  " smooth_phis_numcycle =",smooth_phis_numcycle
       write(iulog,'(a,e13.5)')" smooth_phis_nudt =",smooth_phis_nudt
    endif
    call smooth_phis(phis,elem,hybrid,deriv1,nets,nete,minf,smooth_phis_numcycle)


    do ie=nets,nete
       elem(ie)%state%phis(:,:)=phis(:,:,ie)
    enddo

    end subroutine smooth_topo_datasets
    
  !_____________________________________________________________________
  subroutine set_prescribed_scm(elem,dt,tl)
  
    ! Update the floating levels based on the prescribed
    !  large scale vertical velocity for single column model

    use dimensions_mod, only: qsize
    use time_mod, only: timelevel_qdp
    use control_mod, only: dt_tracer_factor  
    use time_mod,       only: timelevel_t

    type (element_t),      intent(inout), target  :: elem(:) 
    real (kind=real_kind), intent(in)             :: dt
    type (TimeLevel_t)   , intent(in)             :: tl
    
    real (kind=real_kind) :: dp(np,np)! pressure thickness, vflux
    real(kind=real_kind)  :: eta_dot_dpdn(np,np,nlevp)
    
    integer :: ie,k,p,n0,np1,n0_qdp,np1_qdp

    n0    = tl%n0
    np1   = tl%np1

    call TimeLevel_Qdp(tl, dt_tracer_factor, n0_qdp, np1_qdp)
    
    do k=1,nlev
      eta_dot_dpdn(:,:,k)=elem(1)%derived%omega_p(1,1,k)
    enddo  
    eta_dot_dpdn(:,:,nlev+1) = eta_dot_dpdn(:,:,nlev)
    
    do k=1,nlev
      elem(1)%state%dp3d(:,:,k,np1) = elem(1)%state%dp3d(:,:,k,n0) &
        + dt*(eta_dot_dpdn(:,:,k+1) - eta_dot_dpdn(:,:,k))
    enddo

    ! Update temperature variable for the theta-l or preqx implementation
#ifdef MODEL_THETA_L
    do k=1,nlev
      elem(1)%state%vtheta_dp(:,:,k,np1) = (elem(1)%state%vtheta_dp(:,:,k,n0)/ &
                elem(1)%state%dp3d(:,:,k,n0))*elem(1)%state%dp3d(:,:,k,np1)
    enddo
#else
   do k=1,nlev
     elem(1)%state%T(:,:,k,np1) = elem(1)%state%T(:,:,k,n0)
   enddo
#endif

    do p=1,qsize
      do k=1,nlev
        elem(1)%state%Qdp(:,:,k,p,np1_qdp)=elem(1)%state%Q(:,:,k,p)*&
          elem(1)%state%dp3d(:,:,k,np1)
      enddo
    enddo
    
  end subroutine set_prescribed_scm
    
end module prim_driver_base



