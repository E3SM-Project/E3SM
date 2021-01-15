
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module init_mod
contains
  subroutine init(elem, edge1,edge2,edge3,red,par, dom_mt)
    use kinds, only : real_kind, longdouble_kind
    ! --------------------------------
    use thread_mod, only : nthreads, hthreads, omp_set_num_threads
    ! --------------------------------
    use control_mod, only : restartfreq, topology, geometry, partmethod, cubed_sphere_map
    ! --------------------------------
    use namelist_mod, only : readnl
    ! --------------------------------
    use dimensions_mod, only : np, nelem, nlev, nelemd, nelemdmax,  &
	GlobalUniqueCols
    ! -------------------------------- 
    use time_mod, only : time_at, nmax
    ! --------------------------------
    use quadrature_mod, only :  test_gauss, test_gausslobatto, quadrature_t, gausslobatto
    ! --------------------------------
    use element_mod, only : element_t, allocate_element_desc
    ! --------------------------------
    use mass_matrix_mod, only : mass_matrix
    ! --------------------------------

    use mesh_mod, only : MeshUseMeshFile
    use cube_mod,  only : CubeEdgeCount , CubeElemCount, CubeTopology
    ! --------------------------------
    use mesh_mod, only :   MeshSetCoordinates,      &
                           MeshCubeTopology,  &
                           MeshCubeElemCount, &
                           MeshCubeEdgeCount
    ! --------------------------------
    use cube_mod, only : cube_init_atomic, set_corner_coordinates
    ! --------------------------------
    use geometry_mod, only : set_area_correction_map0, set_area_correction_map2
    ! --------------------------------
  use planar_mod,  only : PlaneEdgeCount , PlaneElemCount, PlaneTopology
  ! --------------------------------
  use planar_mesh_mod, only :   PlaneMeshSetCoordinates,      &
                         MeshPlaneTopology
  ! --------------------------------
  use planar_mod, only : plane_init_atomic, plane_set_corner_coordinates

    ! --------------------------------
    use edge_mod, only : initedgebuffer, edge_g
    use edgetype_mod, only : EdgeDescriptor_t, edgebuffer_t
    ! --------------------------------
    use reduction_mod, only : reductionbuffer_ordered_1d_t, red_min, red_flops, red_timer, &
         red_sum, red_sum_int, red_max, red_max_int, InitReductionBuffer
    ! --------------------------------
    use metagraph_mod, only : metavertex_t, metaedge_t, LocalElemCount, &
         initmetagraph
    ! --------------------------------
    use gridgraph_mod, only : gridvertex_t, gridedge_t, PrintGridEdge, PrintGridVertex, allocate_gridvertex_nbrs
    ! --------------------------------
    use schedule_mod, only : genEdgeSched
    ! --------------------------------
    use schedtype_mod, only : schedule
    ! --------------------------------
    use state_mod, only : printstate_init
    ! --------------------------------
    !  use global_norms_mod
    ! --------------------------------
    use parallel_mod, only : iam, parallel_t, haltmp, mpiinteger_t, abortmp,global_shared_buf, nrepro_vars
    ! --------------------------------
    use metis_mod, only : genmetispart
    ! --------------------------------
    use checksum_mod, only : testchecksum
    ! --------------------------------
    use spacecurve_mod, only : genspacepart
    ! --------------------------------
    use dof_mod, only : global_dof, CreateUniqueIndex, SetElemOffset
    ! --------------------------------
    use restart_io_mod, only : restfile
    ! --------------------------------
    use restart_mod, only : initRestartFile
    ! --------------------------------
    use params_mod, only : SFCURVE
    ! ---------------------------------
    use domain_mod, only: domain1d_t, decompose 
    ! ---------------------------------
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
    ! --------------------------------
    use repro_sum_mod, only: repro_sum, repro_sum_defaultopts, repro_sum_setopts
    ! --------------------------------
    !use physical_constants, only : dd_pi
    ! -------------------------------
    !use coordinate_systems_mod, only : sphere_tri_area
    ! --------------------------------
    use common_io_mod, only : homme_pio_init
    ! --------------------------------


    implicit none
#ifdef _MPI
! _EXTERNAL
#include <mpif.h> 
#endif
!   G95  "pointer attribute conflicts with INTENT attribute":  
!    type (element_t), intent(inout), pointer :: elem(:)
    type (element_t), pointer :: elem(:)
        
    type (EdgeBuffer_t)           :: edge1
    type (EdgeBuffer_t)           :: edge2
    type (EdgeBuffer_t)           :: edge3
    type (ReductionBuffer_ordered_1d_t) :: red
    type (parallel_t), intent(in) :: par
    type (domain1d_t), pointer :: dom_mt(:)
    ! Local Variables

    type (GridVertex_t), target,allocatable :: GridVertex(:)
    type (GridEdge_t),   target,allocatable :: Gridedge(:)
    type (MetaVertex_t), target,allocatable :: MetaVertex(:)
    type (MetaEdge_t),   target,allocatable :: MetaEdge(:)

    integer :: ierr
    integer :: iptr,i,ie,j
    integer :: nelem_edge,nedge
    integer :: nete,nets,ithr,nstep
    integer :: pflag,htype
    integer :: nlyr
    integer :: iMv
    integer :: nxyp, istartP
    real(kind=real_kind) :: et,st

    real(kind=real_kind), allocatable       :: mass(:,:,:)

    logical, parameter :: Debug = .FALSE.

    integer,allocatable :: TailPartition(:)
    integer,allocatable :: HeadPartition(:)
    real(kind=real_kind) :: approx_elements_per_task, xtmp
    type (quadrature_t)   :: gp                     ! element GLL points

    logical :: repro_sum_use_ddpdd, repro_sum_recompute
    real(kind=real_kind) :: repro_sum_rel_diff_max

    ! Read in model control information
    ! =====================================
    call t_startf('init')

    call readnl(par)
    call homme_pio_init(par%rank,par%comm)

    if (par%masterproc) then

       ! =============================================
       ! Compute total simulated time...
       ! =============================================

       write(6,*)""
       write(6,*)" total simulated time = ",Time_at(nmax)
       write(6,*)""

       ! =============================================
       ! Perform Gauss/Gauss Lobatto tests...
       ! =============================================

       !print *,'init: before call to test_gauss'
       call test_gauss(np)
       !print *,'init: before call to test_gausslobatto'
       call test_gausslobatto(np)

    end if

    ! ===============================================================
    ! Allocate and initialize the graph (array of GridVertex_t types)
    ! ===============================================================



       if (par%masterproc) then
          write(6,*)"creating topology..."
       end if
       
       if (MeshUseMeshFile) then
          nelem      = MeshCubeElemCount()
          nelem_edge = MeshCubeEdgeCount()
       else
        if (topology=="cube") then
          nelem      = CubeElemCount()
          nelem_edge = CubeEdgeCount()
        else if (topology=="plane") then
          nelem      = PlaneElemCount()
          nelem_edge = PlaneEdgeCount()
        end if
      end if

        approx_elements_per_task = dble(nelem)/dble(par%nprocs)
        if  (approx_elements_per_task < 1.0D0) then
            if(par%masterproc) print *,"number of elements=", nelem
            if(par%masterproc) print *,"number of procs=", par%nprocs
            call abortmp('There is not enough paralellism in the job, that is, there is less than one elements per task.')
        end if

       allocate(GridVertex(nelem))
       allocate(GridEdge(nelem_edge))
       
       do j =1,nelem
          call allocate_gridvertex_nbrs(GridVertex(j))
       end do

       if (MeshUseMeshFile) then
          if (par%masterproc) then
             write(6,*)"Set up grid vertex from mesh..."
          end if
          if (topology=="cube") then
            call MeshCubeTopology(GridEdge,GridVertex)
          else if  (topology=="plane") then
            call MeshPlaneTopology(GridEdge,GridVertex)
          end if
      else
          if (topology=="cube") then
            call CubeTopology(GridEdge,GridVertex)
          else if (topology=="plane") then
           call PlaneTopology(GridEdge,GridVertex)
          end if
      end if

       if(par%masterproc) write(6,*)"...done."
       

    if(par%masterproc) write(6,*)"partitioning graph..."


    if(partmethod .eq. SFCURVE) then 
       call genspacepart(GridEdge,GridVertex)
    else
       call genmetispart(GridEdge,GridVertex)
    endif

    if(par%masterproc) write(6,*)"...done."

    ! ===========================================================
    ! given partition, count number of local element descriptors
    ! ===========================================================

    !DBG call PrintGridVertex(GridVertex)

    nelem_edge=SIZE(GridEdge)
    allocate(TailPartition(nelem_edge))
    allocate(HeadPartition(nelem_edge))
    do i=1,nelem_edge
       ! TailPartition(i)=GridEdge(i)%tail%partition
       ! HeadPartition(i)=GridEdge(i)%head%partition
       TailPartition(i)=GridEdge(i)%tail%processor_number
       HeadPartition(i)=GridEdge(i)%head%processor_number
    enddo
    allocate(MetaVertex(1))
    allocate(Schedule(1))

    ! ====================================================
    !  Generate the communication graph 
    ! ====================================================
    !JMD call initMetaGraph(iam,MetaVertex(1),TailPartition,HeadPartition,GridVertex,GridEdge)

    if(par%masterproc) write(6,*)"initialize MetaGraph..."

    call initMetaGraph(iam,MetaVertex(1),GridVertex,GridEdge)

    if(par%masterproc) write(6,*)"...done."

    nelemd = LocalElemCount(MetaVertex(1))
#ifdef _MPI
    call mpi_allreduce(nelemd,nelemdmax,1, MPIinteger_t,MPI_MAX,par%comm,ierr)
#else
    nelemdmax=nelemd
#endif

    allocate(elem(nelemd))
    call allocate_element_desc(elem)

    ! ====================================================
    !  Generate the communication schedule
    ! ====================================================
    call genEdgeSched(elem, iam,Schedule(1),MetaVertex(1))
    !call PrintSchedule(Schedule(1))
    
    deallocate(TailPartition)
    deallocate(HeadPartition)

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

    ! initial 1D grids used to form tensor product element grids:
    gp=gausslobatto(np)

       ! ========================================================================
       ! Note it is more expensive to initialize each individual spectral element 
       ! ========================================================================
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
       
    ! This routine does not check whether gp is init-ed.
    if(( cubed_sphere_map == 2 ).AND.( np > 2 )) then
       call set_area_correction_map2(elem, nelemd, par, gp)
    endif

    deallocate(gp%points)
    deallocate(gp%weights)

    ! =================================================================
    ! Run the checksum to verify communication schedule
    ! =================================================================

    !call testchecksum(elem,par,GridEdge)  ! broken

    ! =================================================================
    ! Initialize mass_matrix
    ! =================================================================

    call mass_matrix(par,elem)

    ! =================================================================
    ! Determine the global degree of freedome for each gridpoint 
    ! =================================================================

    call global_dof(par,elem)

    ! =================================================================
    ! Create Unique Indices 
    ! =================================================================
    
    do ie=1,nelemd
       call CreateUniqueIndex(elem(ie)%GlobalId,elem(ie)%gdofP,elem(ie)%idxP)
    enddo
    call SetElemOffset(par,elem,GlobalUniqueCols)

    !JMD call PrintDofV(elem)
    !JMD call PrintDofP(elem)

    hthreads = min(nthreads,nelemd)
    allocate(dom_mt(0:hthreads-1))
    do ithr=0,hthreads-1
       dom_mt(ithr)=decompose(1,nelemd,hthreads,ithr)
    end do

    ! =================================================================
    ! Initialize shared boundary_exchange and reduction buffers
    ! =================================================================
    call initEdgeBuffer(par,edge_g,elem,max(2,nlev))
    call initEdgeBuffer(par,edge1,elem,nlev)
    call initEdgeBuffer(par,edge2,elem,2*nlev)
    call initEdgeBuffer(par,edge3,elem,11*nlev)

    allocate(global_shared_buf(nelemd,nrepro_vars))

    call InitReductionBuffer(red,3*nlev,hthreads)
    call InitReductionBuffer(red_sum,1)
    call InitReductionBuffer(red_sum_int,1)
    call InitReductionBuffer(red_max,1)
    call InitReductionBuffer(red_min,1)
    call InitReductionBuffer(red_flops,1)

    call printstate_init()
    ! Initialize output fields for plotting...

    allocate(mass(np,np,nelemd))

    ! ==========================================================
    !  This routines initalizes a Restart file.  This involves:
    !      I)  Setting up the MPI datastructures
    ! ==========================================================
    if(restartfreq > 0) then
       call initRestartFile(elem(1)%state,par,RestFile)
    endif
    
    call t_stopf('init')
  end subroutine init


end module init_mod
