#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module init_mod
contains
! not a nice way to integrate fvm but otherwise DG does not work anymore
#ifdef _FVM
  subroutine init(elem, edge1,edge2,edge3,red,par, fvm)
#else  
  subroutine init(elem, edge1,edge2,edge3,red,par)
#endif
    use kinds, only : real_kind, longdouble_kind
    ! --------------------------------
    use thread_mod, only : nthreads
    ! --------------------------------
    use control_mod, only : filter_counter, restartfreq, topology, &
          partmethod, while_iter
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
    use element_mod, only : element_t
    ! --------------------------------
    use mass_matrix_mod, only : mass_matrix
    ! --------------------------------
    use mesh_mod, only : MeshUseMeshFile
#ifndef MESH
     use cube_mod,  only : cubeedgecount , cubeelemcount, cubetopology
#else
    ! --------------------------------
    use mesh_mod, only :   MeshSetCoordinates,      &
                           MeshCubeTopology,  &
                           MeshCubeElemCount, &
                           MeshCubeEdgeCount
#endif 
    ! --------------------------------
    use cube_mod, only : cube_init_atomic, rotation_init_atomic, set_corner_coordinates, &
         assign_node_numbers_to_elem

    ! --------------------------------
    use edge_mod, only : edgebuffer_t, initedgebuffer   
    ! --------------------------------
    use reduction_mod, only : reductionbuffer_ordered_1d_t, red_min, red_flops, red_timer, &
         red_sum, red_sum_int, red_max, red_max_int, InitReductionBuffer
    ! --------------------------------
    use metagraph_mod, only : metavertex_t, metaedge_t, LocalElemCount, &
         initmetagraph
    ! --------------------------------
    use gridgraph_mod, only : gridvertex_t, gridedge_t
    ! --------------------------------
    use schedule_mod, only : schedule, genEdgeSched
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
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
    ! --------------------------------
#ifdef _FVM 
  use fvm_mod, only : fvm_init1
  use fvm_control_volume_mod, only : fvm_struct
  use dimensions_mod, only : nc
#endif

    implicit none
#ifdef _MPI
#include <mpif.h> ! _EXTERNAL
#endif
!   G95  "pointer attribute conflicts with INTENT attribute":  
!    type (element_t), intent(inout), pointer :: elem(:)
    type (element_t), pointer :: elem(:)
#ifdef _FVM    
    type (fvm_struct), pointer, optional :: fvm(:)
#endif
        
    type (EdgeBuffer_t)           :: edge1
    type (EdgeBuffer_t)           :: edge2
    type (EdgeBuffer_t)           :: edge3
    type (ReductionBuffer_ordered_1d_t) :: red
    type (parallel_t), intent(in) :: par
    ! Local Variables

    type (GridVertex_t), target,allocatable :: GridVertex(:)
    type (GridEdge_t),   target,allocatable :: Gridedge(:)
    type (MetaVertex_t), target,allocatable :: MetaVertex(:)
    type (MetaEdge_t),   target,allocatable :: MetaEdge(:)

    integer :: ierr
    integer :: iptr,i,ie
    integer :: nelem_edge,nedge
    integer :: nete,nets,ithr,nstep
    integer :: pflag,htype
    integer :: nlyr
    integer :: iMv
    integer :: nxyp, istartP
    real(kind=real_kind) :: et,st

    character(len=80) rot_type   ! cube edge rotation type
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

    ! =====================================
    ! Set cube edge rotation type for model
    ! =====================================

    rot_type="contravariant"

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
#ifdef MESH
          nelem      = MeshCubeElemCount()
          nelem_edge = MeshCubeEdgeCount() 
#else
          call abortmp('Input file requires compilation with CPP macro MESH, but mesh support was not built in. Aborting.')
#endif

       else

#ifndef MESH 
          nelem      = CubeElemCount()
          nelem_edge = CubeEdgeCount() 
#else
          call abortmp('Input file does not require an external mesh file, yet the standard cube topology was not built in. Aborting.')
#endif
       end if


        approx_elements_per_task = dble(nelem)/dble(par%nprocs)
        if  (approx_elements_per_task < 1.0D0) then
            if(par%masterproc) print *,"number of elements=", nelem
            if(par%masterproc) print *,"number of procs=", par%nprocs
            call abortmp('There is not enough paralellism in the job, that is, there is less than one elements per task.')
        end if
       allocate(GridVertex(nelem))
       allocate(GridEdge(nelem_edge))

       if (MeshUseMeshFile) then
#ifdef MESH
          if (par%masterproc) then
             write(6,*)"Set up grid vertex from mesh..."
          end if
          call MeshCubeTopology(GridEdge,GridVertex)
#else
          call abortmp('Input file requires compilation with CPP macro MESH, but mesh support was not built in. Aborting.')
#endif 
       else 
#ifndef MESH    
          call CubeTopology(GridEdge,GridVertex)
#else
          call abortmp('Input file does not require an external mesh file, yet the standard cube topology was not built in. Aborting.')
#endif
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
#ifdef _PREDICT
    allocate(MetaVertex(npart))
    allocate(Schedule(npart))
#ifdef _MPI
    call haltmp("init: PREDICT code branch not supported under MPI")
#else
    nelemd = nelem
    print *,'init: before allocation of elem nelemd:',nelemd
    !JMDallocate(elem(nelemd))
    print *,'init: after allocation of elem'
#endif
    do iMv = 1,npart
       ! ====================================================
       !  Generate the communication graph 
       ! ====================================================
       !JMDcall initMetaGraph(iMv,MetaVertex(iMv),TailPartition,HeadPartition,GridVertex,GridEdge)
       call initMetaGraph(iMv,MetaVertex(iMv),GridVertex,GridEdge)

       !JMD nelemd = LocalElemCount(MetaVertex(iMv))

       ! ====================================================
       !  Generate the communication schedule
       ! ====================================================
       if (MODULO(iMv,100) ==0) print *,'init: before call to genEdgeSched iMv :',iMv
       call genEdgeSched(iMv,Schedule(iMv),MetaVertex(iMv))

    enddo

    print *,'init: before MessageStats'
    call MessageStats(nlyr)
    print *,'init: after MessageStats'
#else
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

#ifdef _FVM
    allocate(fvm(nelemd))
#endif

    ! ====================================================
    !  Generate the communication schedule
    ! ====================================================
    call genEdgeSched(elem, iam,Schedule(1),MetaVertex(1))

#endif
    deallocate(TailPartition)
    deallocate(HeadPartition)

    ! initial 1D grids used to form tensor product element grids:
    gp=gausslobatto(np)

    if (topology=="cube") then
       ! ========================================================================
       ! Note it is more expensive to initialize each individual spectral element 
       ! ========================================================================
       if(par%masterproc) write(6,*)"initializing cube elements..."
#ifdef MESH
       if (MeshUseMeshFile) then
          call MeshSetCoordinates(elem)
       else
          do ie=1,nelemd
             call set_corner_coordinates(elem(ie))
          enddo
          call assign_node_numbers_to_elem(elem, GridVertex)
       endif
#else
       do ie=1,nelemd
          call set_corner_coordinates(elem(ie))
       enddo
       call assign_node_numbers_to_elem(elem, GridVertex)
#endif 
       do ie=1,nelemd
          call cube_init_atomic(elem(ie),gp%points)
          call rotation_init_atomic(elem(ie),rot_type)
       enddo
       if(par%masterproc) write(6,*)"...done."
    end if
    ! =================================================================
    ! Run the checksum to verify communication schedule
    ! =================================================================
!    call testchecksum(par,GridEdge)
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

    ! =================================================================
    ! Initialize shared boundary_exchange and reduction buffers
    ! =================================================================

    call initEdgeBuffer(edge1,nlev)
#ifdef _PRIMDG
    call initEdgeBuffer(edge2,4*nlev)
    call initEdgeBuffer(edge3,11*nlev)
#else
    call initEdgeBuffer(edge2,2*nlev)
    call initEdgeBuffer(edge3,11*nlev)
#endif
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

#ifdef _HTRACE
    htype = 19
    pflag = 1
    call TRACE_INIT(htype,pflag)
#endif
    while_iter = 0
    filter_counter = 0

    ! ==========================================================
    !  This routines initalizes a Restart file.  This involves:
    !      I)  Setting up the MPI datastructures
    ! ==========================================================
    if(restartfreq > 0) then
       call initRestartFile(elem(1)%state,par,RestFile)
    endif
    !DBG  print *,'prim_init: after call to initRestartFile'
  
#ifdef _FVM
  call fvm_init1(par)
#endif    
    
    call t_stopf('init')
  end subroutine init


end module init_mod
