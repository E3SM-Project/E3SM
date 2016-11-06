
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module init_mod
contains
! not a nice way to integrate fvm but otherwise DG does not work anymore
  subroutine init(elem, edge1,edge2,edge3,red,par, dom_mt, fvm)
    use kinds, only : real_kind, longdouble_kind
    ! --------------------------------
    use thread_mod, only : nthreads, omp_set_num_threads
    ! --------------------------------
    use control_mod, only : restartfreq, topology, partmethod
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
    use cube_mod,  only : cubeedgecount , cubeelemcount, cubetopology
    ! --------------------------------
    use mesh_mod, only :   MeshSetCoordinates,      &
                           MeshCubeTopology,  &
                           MeshCubeElemCount, &
                           MeshCubeEdgeCount
    ! --------------------------------
    use cube_mod, only : cube_init_atomic, set_corner_coordinates, &
         assign_node_numbers_to_elem

    ! --------------------------------
    use edge_mod, only : initedgebuffer
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
    use thread_mod, only : nThreadsHoriz, omp_get_num_threads
    ! ---------------------------------
    use domain_mod, only: domain1d_t, decompose 
    ! ---------------------------------
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
    ! --------------------------------
    use fvm_mod, only : fvm_init1
    use fvm_control_volume_mod, only : fvm_struct
    use dimensions_mod, only : nc, ntrac

    implicit none
#ifdef _MPI
! _EXTERNAL
#include <mpif.h> 
#endif
!   G95  "pointer attribute conflicts with INTENT attribute":  
!    type (element_t), intent(inout), pointer :: elem(:)
    type (element_t), pointer :: elem(:)
    type (fvm_struct), pointer, optional :: fvm(:)
        
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
    integer :: n_domains
    real(kind=real_kind) :: et,st

    real(kind=real_kind), allocatable       :: mass(:,:,:)

    logical, parameter :: Debug = .FALSE.

    integer,allocatable :: TailPartition(:)
    integer,allocatable :: HeadPartition(:)
    real(kind=real_kind) :: approx_elements_per_task, xtmp
    type (quadrature_t)   :: gp                     ! element GLL points


    ! =====================================
    ! Read in model control information
    ! =====================================
    call t_startf('init')

    call readnl(par)

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


    if (topology=="cube") then

       if (par%masterproc) then
          write(6,*)"creating cube topology..."
       end if
       
       if (MeshUseMeshFile) then
          nelem      = MeshCubeElemCount()
          nelem_edge = MeshCubeEdgeCount() 
       else
          nelem      = CubeElemCount()
          nelem_edge = CubeEdgeCount() 
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
          call MeshCubeTopology(GridEdge,GridVertex)
       else 
          call CubeTopology(GridEdge,GridVertex)
       end if

       if(par%masterproc) write(6,*)"...done."
       
    end if

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

    call initMetaGraph(iam,MetaVertex(1),GridVertex,GridEdge)

    nelemd = LocalElemCount(MetaVertex(1))
#ifdef _MPI
    call mpi_allreduce(nelemd,nelemdmax,1, MPIinteger_t,MPI_MAX,par%comm,ierr)
#else
    nelemdmax=nelemd
#endif

    allocate(elem(nelemd))
    call allocate_element_desc(elem)
    if (present(fvm)) then
    if (ntrac>0) then
       allocate(fvm(nelemd))
    else
       allocate(fvm(0))  ! create mesh, even if no cslam tracers
    endif
    endif

    ! ====================================================
    !  Generate the communication schedule
    ! ====================================================
    call genEdgeSched(elem, iam,Schedule(1),MetaVertex(1))
    !call PrintSchedule(Schedule(1))
    
    deallocate(TailPartition)
    deallocate(HeadPartition)

    ! initial 1D grids used to form tensor product element grids:
    gp=gausslobatto(np)
    if (topology=="cube") then
       ! ========================================================================
       ! Note it is more expensive to initialize each individual spectral element 
       ! ========================================================================
       if(par%masterproc) write(6,*)"initializing cube elements..."

       if (MeshUseMeshFile) then
          call MeshSetCoordinates(elem)
       else
          do ie=1,nelemd
             call set_corner_coordinates(elem(ie))
          enddo
          call assign_node_numbers_to_elem(elem, GridVertex)
       endif

       do ie=1,nelemd
          call cube_init_atomic(elem(ie),gp%points)
       enddo
       if(par%masterproc) write(6,*)"...done."
    end if
    ! =================================================================
    ! Run the checksum to verify communication schedule
    ! =================================================================
    call testchecksum(elem,par,GridEdge)

    

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

    ! use a copy.  has to be done *after* SetElemOffset()
    do ie=1,nelemd
       elem(ie)%idxV=>elem(ie)%idxP
    end do

    !JMD call PrintDofV(elem)
    !JMD call PrintDofP(elem)

    n_domains = min(Nthreads,nelemd)
    call omp_set_num_threads(n_domains)
    allocate(dom_mt(0:n_domains-1))
    do ithr=0,NThreads-1
       dom_mt(ithr)=decompose(1,nelemd,NThreads,ithr)
    end do
    nThreadsHoriz = NThreads

    ! =================================================================
    ! Initialize shared boundary_exchange and reduction buffers
    ! =================================================================
    call initEdgeBuffer(par,edge1,elem,nlev)
    call initEdgeBuffer(par,edge2,elem,2*nlev)
    call initEdgeBuffer(par,edge3,elem,11*nlev)

    allocate(global_shared_buf(nelemd,nrepro_vars))

    call InitReductionBuffer(red,3*nlev,nthreads)
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
    if (ntrac>0) call fvm_init1(par,elem)
    
    call t_stopf('init')
  end subroutine init


end module init_mod
