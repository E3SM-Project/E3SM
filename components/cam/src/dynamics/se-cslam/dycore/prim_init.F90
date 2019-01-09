module prim_init

  use shr_kind_mod,   only: r8=>shr_kind_r8
  use dimensions_mod, only: nc
  use reduction_mod,  only: reductionbuffer_ordered_1d_t
  use quadrature_mod, only: quadrature_t, gausslobatto

  implicit none
  private
  save

  public :: prim_init1

  real(r8), public :: fvm_corners(nc+1) ! fvm cell corners on reference element
  real(r8), public :: fvm_points(nc)    ! fvm cell centers on reference element

  type (quadrature_t), public               :: gp    ! element GLL points
  type (ReductionBuffer_ordered_1d_t)       :: red   ! reduction buffer (shared)

contains
  subroutine prim_init1(elem, fvm, par, Tl)
    use cam_logfile,            only: iulog
    use shr_sys_mod,            only: shr_sys_flush
    use thread_mod,             only: max_num_threads
    use dimensions_mod,         only: np, nlev, nelem, nelemd, nelemdmax
    use dimensions_mod,         only: GlobalUniqueCols, fv_nphys,irecons_tracer
    use control_mod,            only: topology, partmethod
    use element_mod,            only: element_t, allocate_element_desc
    use fvm_mod,                only: fvm_init1
    use mesh_mod,               only: MeshUseMeshFile
    use time_mod,               only: timelevel_init, timelevel_t
    use mass_matrix_mod,        only: mass_matrix
    use derivative_mod,         only: allocate_subcell_integration_matrix_cslam
    use derivative_mod,         only: allocate_subcell_integration_matrix_physgrid
    use cube_mod,               only: cubeedgecount , cubeelemcount, cubetopology
    use cube_mod,               only: cube_init_atomic, rotation_init_atomic, set_corner_coordinates
    use cube_mod,               only: assign_node_numbers_to_elem
    use mesh_mod,               only: MeshSetCoordinates, MeshUseMeshFile, MeshCubeTopology
    use mesh_mod,               only: MeshCubeElemCount, MeshCubeEdgeCount
    use metagraph_mod,          only: metavertex_t, localelemcount, initmetagraph, printmetavertex
    use gridgraph_mod,          only: gridvertex_t, gridedge_t
    use gridgraph_mod,          only: allocate_gridvertex_nbrs, deallocate_gridvertex_nbrs
    use schedtype_mod,          only: schedule
    use schedule_mod,           only: genEdgeSched
    use prim_advection_mod,     only: prim_advec_init1
    use cam_abortutils,         only: endrun
    use spmd_utils,             only: mpi_integer, mpi_max
    use parallel_mod,           only: parallel_t, syncmp, global_shared_buf, nrepro_vars
    use spacecurve_mod,         only: genspacepart
    use dof_mod,                only: global_dof, CreateUniqueIndex, SetElemOffset
    use params_mod,             only: SFCURVE
    use physconst,              only: pi
    use reduction_mod,          only: red_min, red_max, red_max_int, red_flops
    use reduction_mod,          only: red_sum, red_sum_int, initreductionbuffer
    use infnan,                 only: nan, assignment(=)
    use shr_reprosum_mod,       only: repro_sum => shr_reprosum_calc
    use fvm_analytic_mod,       only: compute_basic_coordinate_vars
    use fvm_control_volume_mod, only: fvm_struct, allocate_physgrid_vars

    type(element_t),  pointer        :: elem(:)
    type(fvm_struct), pointer        :: fvm(:)
    type(parallel_t),  intent(inout) :: par
    type(timelevel_t), intent(out)   :: Tl

    ! Local Variables
    type (GridVertex_t), target,allocatable :: GridVertex(:)
    type (GridEdge_t),   target,allocatable :: Gridedge(:)
    type (MetaVertex_t), target,allocatable :: MetaVertex(:)

    integer                        :: ie
    integer                        :: nets, nete
    integer                        :: nelem_edge
    integer                        :: ierr, j
    logical,           parameter   :: Debug = .FALSE.

    real(r8),          allocatable :: aratio(:,:)
    real(r8)                       :: area(1), xtmp
    character(len=80)              :: rot_type ! cube edge rotation type

    integer                        :: i

    character(len=128)             :: errmsg
    character(len=*),  parameter   :: subname = 'PRIM_INIT1: '

    ! ====================================
    ! Set cube edge rotation type for model
    ! unnecessary complication here: all should
    ! be on the same footing. RDL
    ! =====================================
    rot_type = "contravariant"

    ! ===============================================================
    ! Allocate and initialize the graph (array of GridVertex_t types)
    ! ===============================================================

    if (topology=="cube") then

      if (par%masterproc) then
        write(iulog,*) subname, "creating cube topology..."
        call shr_sys_flush(iulog)
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

      do j = 1, nelem
        call allocate_gridvertex_nbrs(GridVertex(j))
      end do

      if (MeshUseMeshFile) then
        if (par%masterproc) then
          write(iulog,*) subname, "Set up grid vertex from mesh..."
        end if
        call MeshCubeTopology(GridEdge, GridVertex)
      else
        call CubeTopology(GridEdge,GridVertex)
      end if

      if (par%masterproc) then
        write(iulog,*)"...done."
      end if
    end if
    if(par%masterproc) then
      write(iulog,*) subname, "total number of elements nelem = ",nelem
    end if

    if(partmethod == SFCURVE) then
      if(par%masterproc) then
        write(iulog,*) subname, "partitioning graph using SF Curve..."
      end if
      call genspacepart(GridVertex)
    else
      write(errmsg, *) 'Unsupported partition method, ',partmethod
      call endrun(subname//trim(errmsg))
    end if

    ! ===========================================================
    ! given partition, count number of local element descriptors
    ! ===========================================================
    allocate(MetaVertex(1))
    allocate(Schedule(1))

    nelem_edge = SIZE(GridEdge)

    ! ====================================================
    !  Generate the communication graph
    ! ====================================================
    call initMetaGraph(par%rank+1,MetaVertex(1),GridVertex,GridEdge)

    nelemd = LocalElemCount(MetaVertex(1))
    if (par%masterproc .and. Debug) then
      call PrintMetaVertex(MetaVertex(1))
    endif

    if(nelemd <= 0) then
      call endrun(subname//'Not yet ready to handle nelemd = 0 yet' )
    end if
    call mpi_allreduce(nelemd, nelemdmax, 1, MPI_INTEGER, MPI_MAX, par%comm, ierr)

    if (nelemd > 0) then
      allocate(elem(nelemd))
      call allocate_element_desc(elem)
    end if

    if (fv_nphys > 0) then
      allocate(fvm(nelemd))
      call allocate_physgrid_vars(fvm,par)
    else
      ! Even if fvm not needed, still desirable to allocate it as empty
      ! so it can be passed as a (size zero) array rather than pointer.
      allocate(fvm(0))
    end if

    ! ====================================================
    !  Generate the communication schedule
    ! ====================================================

    call genEdgeSched(par, elem, par%rank+1, Schedule(1), MetaVertex(1))

    allocate(global_shared_buf(nelemd, nrepro_vars))
    global_shared_buf = 0.0_r8

    call syncmp(par)

    ! =================================================================
    ! Set number of domains (for 'decompose') equal to number of threads
    !  for OpenMP across elements, equal to 1 for OpenMP within element
    ! =================================================================

    ! =================================================================
    ! Initialize shared boundary_exchange and reduction buffers
    ! =================================================================
    if(par%masterproc) then
      write(iulog,*) subname, 'init shared boundary_exchange buffers'
      call shr_sys_flush(iulog)
    end if
    call InitReductionBuffer(red,3*nlev,max_num_threads)
    call InitReductionBuffer(red_sum,5)
    call InitReductionBuffer(red_sum_int,1)
    call InitReductionBuffer(red_max,1)
    call InitReductionBuffer(red_max_int,1)
    call InitReductionBuffer(red_min,1)
    call initReductionBuffer(red_flops,1)

    gp = gausslobatto(np)  ! GLL points

    ! fvm nodes are equally spaced in alpha/beta
    ! HOMME with equ-angular gnomonic projection maps alpha/beta space
    ! to the reference element via simple scale + translation
    ! thus, fvm nodes in reference element [-1,1] are a tensor product of
    ! array 'fvm_corners(:)' computed below:
    xtmp = nc
    do i = 1, nc+1
      fvm_corners(i)= 2*(i-1)/xtmp - 1  ! [-1,1] including end points
    end do
    do i = 1, nc
      fvm_points(i)= ( fvm_corners(i)+fvm_corners(i+1) ) /2
    end do

    if (topology == "cube") then
      if(par%masterproc) then
        write(iulog,*) subname, "initializing cube elements..."
        call shr_sys_flush(iulog)
      end if
      if (MeshUseMeshFile) then
        call MeshSetCoordinates(elem)
      else
        do ie = 1, nelemd
          call set_corner_coordinates(elem(ie))
        end do
        call assign_node_numbers_to_elem(elem, GridVertex)
      end if
      do ie = 1, nelemd
        call cube_init_atomic(elem(ie),gp%points)
      end do
    end if

    ! =================================================================
    ! Initialize mass_matrix
    ! =================================================================
    if(par%masterproc) then
      write(iulog,*) subname, 'running mass_matrix'
      call shr_sys_flush(iulog)
    end if
    call mass_matrix(par, elem)
    allocate(aratio(nelemd,1))

    if (topology == "cube") then
      area = 0
      do ie = 1, nelemd
        aratio(ie,1) = sum(elem(ie)%mp(:,:)*elem(ie)%metdet(:,:))
      end do
      call repro_sum(aratio, area, nelemd, nelemd, 1, commid=par%comm)
      area(1) = 4.0_r8*pi/area(1)  ! ratio correction
      deallocate(aratio)
      if (par%masterproc) then
        write(iulog,'(2a,f20.17)') subname, "re-initializing cube elements: area correction=", area(1)
        call shr_sys_flush(iulog)
      end if

      do ie = 1, nelemd
        call cube_init_atomic(elem(ie),gp%points,area(1))
        call rotation_init_atomic(elem(ie),rot_type)
      end do
    end if

    if(par%masterproc) then
      write(iulog,*) subname, 're-running mass_matrix'
      call shr_sys_flush(iulog)
    end if
    call mass_matrix(par, elem)

    ! =================================================================
    ! Determine the global degree of freedome for each gridpoint
    ! =================================================================
    if(par%masterproc) then
      write(iulog,*) subname, 'running global_dof'
      call shr_sys_flush(iulog)
    end if
    call global_dof(par, elem)

    ! =================================================================
    ! Create Unique Indices
    ! =================================================================

    do ie = 1, nelemd
      call CreateUniqueIndex(elem(ie)%GlobalId,elem(ie)%gdofP,elem(ie)%idxP)
    end do

    call SetElemOffset(par,elem, GlobalUniqueCols)

    do ie = 1, nelemd
      elem(ie)%idxV=>elem(ie)%idxP
    end do

    ! initialize flux terms to 0
    do ie = 1, nelemd
      elem(ie)%derived%FM=0.0_r8
      elem(ie)%derived%FQ=0.0_r8
      elem(ie)%derived%FT=0.0_r8
      elem(ie)%derived%pecnd=0.0_r8

      elem(ie)%derived%Omega=0
      elem(ie)%state%dp3d=0

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

      elem(ie)%derived%Utnd=0._r8
      elem(ie)%derived%Vtnd=0._r8
      elem(ie)%derived%Ttnd=0._r8
    end do

    ! ==========================================================
    !  This routines initalizes a Restart file.  This involves:
    !      I)  Setting up the MPI datastructures
    ! ==========================================================
    deallocate(GridEdge)
    do j = 1, nelem
      call deallocate_gridvertex_nbrs(GridVertex(j))
    end do
    deallocate(GridVertex)

    do j = 1, MetaVertex(1)%nmembers
      call deallocate_gridvertex_nbrs(MetaVertex(1)%members(j))
    end do
    deallocate(MetaVertex)

    ! =====================================
    ! Set number of threads...
    ! =====================================
    if(par%masterproc) then
      write(iulog,*) subname, "max_num_threads=",max_num_threads
      call shr_sys_flush(iulog)
    end if

    nets = 1
    nete = nelemd
    call Prim_Advec_Init1(par, elem)
    if (fv_nphys > 0) then
      call fvm_init1(par,elem)
    end if

    ! =======================================================
    ! Allocate memory for subcell flux calculations.
    ! =======================================================
    call allocate_subcell_integration_matrix_cslam(np, nc)
    if (fv_nphys > 0) then
      call allocate_subcell_integration_matrix_physgrid(np, fv_nphys)
    end if

    call TimeLevel_init(tl)

    if (fv_nphys > 0) then
      if(par%masterproc) then
        write(iulog,*) subname, 'initialize basic fvm coordinate variables'
        call shr_sys_flush(iulog)
      end if
      do ie = 1, nelemd
        call compute_basic_coordinate_vars(elem(ie), nc, irecons_tracer,      &
             fvm(ie)%dalpha, fvm(ie)%dbeta, fvm(ie)%vtx_cart(:,:,1:nc,1:nc),  &
             fvm(ie)%center_cart(1:nc,1:nc), fvm(ie)%area_sphere(1:nc,1:nc),  &
             fvm(ie)%spherecentroid(:,1:nc,1:nc))
        call compute_basic_coordinate_vars(elem(ie), fv_nphys, irecons_tracer,&
             fvm(ie)%dalpha_physgrid, fvm(ie)%dbeta_physgrid,                 &
             fvm(ie)%vtx_cart_physgrid   (:,:,1:fv_nphys,1:fv_nphys),         &
             fvm(ie)%center_cart_physgrid(1:fv_nphys,1:fv_nphys),             &
             fvm(ie)%area_sphere_physgrid(1:fv_nphys,1:fv_nphys),             &
             fvm(ie)%spherecentroid_physgrid(:,1:fv_nphys,1:fv_nphys))
      end do
    end if

    if(par%masterproc) then
      write(iulog,*) subname, 'end of prim_init'
      call shr_sys_flush(iulog)
    end if
  end subroutine prim_init1
end module prim_init
