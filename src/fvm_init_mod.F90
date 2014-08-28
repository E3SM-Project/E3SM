#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module fvm_init_mod
contains
 subroutine fvm_init(elem,red,par)
    use kinds, only : real_kind, longdouble_kind
    ! --------------------------------
    use thread_mod, only : nthreads
    ! --------------------------------
    use control_mod, only : filter_counter, restartfreq, topology, &
          partmethod, while_iter
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
    use cube_mod, only : cubeelemcount, cubeedgecount, &
        cube_init_atomic, rotation_init_atomic, set_corner_coordinates, &
        assign_node_numbers_to_elem
    use cube_mod, only : cubetopology
    ! --------------------------------
    use edge_mod, only : edgebuffer_t, initedgebuffer   
    ! --------------------------------
    use reduction_mod, only : reductionbuffer_ordered_1d_t, red_min, red_flops, red_timer, &
         red_sum, red_sum_int, red_max, red_max_int, InitReductionBuffer
    ! --------------------------------
    use metagraph_mod, only : metavertex_t, metaedge_t, LocalElemCount, &
         initmetagraph
    ! --------------------------------
    use gridgraph_mod, only : gridvertex_t, gridedge_t, allocate_gridvertex_nbrs
    ! --------------------------------
    use schedule_mod,  only : genEdgeSched
    use schedtype_mod, only : Schedule_t, schedule
    ! --------------------------------
    use state_mod, only : printstate_init
    ! --------------------------------
    !  use global_norms_mod
    ! --------------------------------
    use parallel_mod, only : iam, parallel_t, haltmp, mpiinteger_t, abortmp, global_shared_buf, nrepro_vars
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

    implicit none
     ! _EXTERNAL
#ifdef _MPI
#include <mpif.h>
#endif
    type (element_t), pointer :: elem(:)
    
    type (ReductionBuffer_ordered_1d_t) :: red
    type (parallel_t), intent(in) :: par
    ! Local Variables

    type (GridVertex_t), target,allocatable :: GridVertex(:)
    type (GridEdge_t),   target,allocatable :: Gridedge(:)
    type (MetaVertex_t), target,allocatable :: MetaVertex(:)
    type (MetaEdge_t),   target,allocatable :: MetaEdge(:)

    integer :: ierr
    integer :: iptr,ii,ie,j
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
    integer :: i
    real(kind=real_kind) :: approx_elements_per_task
    type (quadrature_t)   :: gp                     ! element GLL points

    ! =====================================
    ! Read in model control information
    ! =====================================
    call t_startf('fvm_init')

    call fvm_readnl(par)

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

    end if

    ! ===============================================================
    ! Allocate and initialize the graph (array of GridVertex_t types)
    ! ===============================================================


     if (par%masterproc) then
        write(6,*)"creating cube topology..."
     end if


      nelem      = CubeElemCount()
      nelem_edge = CubeEdgeCount() 
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

     call CubeTopology(GridEdge,GridVertex)
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
    call allocate_element_desc(elem)


    ! ====================================================
    !  Generate the communication schedule
    ! ====================================================
    call genEdgeSched(elem, iam,Schedule(1),MetaVertex(1))

#endif
    deallocate(TailPartition)
    deallocate(HeadPartition)

     ! ========================================================================
     ! Note it is more expensive to initialize each individual spectral element 
     ! ========================================================================
     if(par%masterproc) write(6,*)"initializing cube elements..."

        do ie=1,nelemd
           call set_corner_coordinates(elem(ie))
        enddo
        call assign_node_numbers_to_elem(elem, GridVertex)

     ! initialize the GLL grid 
     gp=gausslobatto(np)
     do ie=1,nelemd
        call cube_init_atomic(elem(ie),gp%points)
        call rotation_init_atomic(elem(ie),rot_type)
     enddo
     if(par%masterproc) write(6,*)"...done."
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

    call t_stopf('fvm_init')
  end subroutine fvm_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CE: my own readnl for fvm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fvm_readnl(par)
  !-----------------
  use kinds, only : real_kind, iulog
  !-----------------
  use params_mod, only : recursive, sfcurve
  !-----------------
  use control_mod, only : &
       MAX_STRING_LEN,& 
       partmethod,    &       ! Mesh partitioning method (METIS)
       topology,      &       ! Mesh topology
       test_case,     &       ! test case
       uselapi,       &
       multilevel,    &    
       numnodes,      & 
       sub_case,      &
       tasknum,	      &       ! used dg model in AIX machine
       remapfreq,     &       ! number of steps per remapping call
#ifdef _PRIMDG
       remap_type,    &       ! selected remapping option
#endif
       statefreq,     &       ! number of steps per printstate call
       accumfreq,     &       ! frequency in steps of accumulation
       accumstart,    &       ! model day to start accumulating state variables
       accumstop,     &       ! model day to stop  accumulating state variables
       restartfreq,   &       
       restartfile,   &       ! name of the restart file for INPUT
       restartdir,    &       ! name of the restart directory for OUTPUT
       runtype,       &
       integration,   &       ! integration method
       tracer_advection_formulation, &   ! conservation or non-conservation formulaton
       cubed_sphere_map, &
       compute_mean_flux, &
       qsplit, &
       rk_stage_user, &
       LFTfreq,       &
       TRACERADV_TOTAL_DIVERGENCE, &
       TRACERADV_UGRADQ, &
       prescribed_wind, &
       ftype,        &
       energy_fixer,        &
       limiter_option, &
       nu,            &
       nu_s,          &
       nu_q,          &
       nu_p,          &
       nu_top,        &
       psurf_vis,    &  
       hypervis_order,    &  
       hypervis_power,    &  
       hypervis_subcycle, &
       hypervis_subcycle_q, &
       smooth_phis_numcycle, &
       smooth_sgh_numcycle, &
       smooth_phis_nudt, &
       initial_total_mass, &  ! set > 0 to set the initial_total_mass
       u_perturb,     &          ! J&W bareclinic test perturbation size
       columnpackage, &
       moisture,      &
       filter_type,   &
       transfer_type, &
       filter_freq,   &
       filter_mu,     &
       filter_freq_advection,   &
       filter_mu_advection,     &
       p_bv,          &
       s_bv,          &
       wght_fm,       &
       kcut_fm,       &
       vform,           &
       vfile_mid,       &
       vfile_int,       &    
       precon_method, &
       maxits,        &
       tol,           &
       debug_level,   &
       tracer_transport_type,           &
       TRACERTRANSPORT_SE_GLL,          &
       TRACERTRANSPORT_SEMILAGRANG_GLL, &
       TRACERTRANSPORT_LAGRANGIAN_FVM,  &
       TRACERTRANSPORT_FLUXFORM_FVM,    &
       TRACERTRANSPORT_SPELT_FVM,       &
       tracer_grid_type,                &
       TRACER_GRIDTYPE_GLL,             &
       TRACER_GRIDTYPE_FVM,             &
       test_cfldep                          
  !-----------------
  use thread_mod, only : nthreads
  !-----------------
  use dimensions_mod, only : ne, np,nc,ntrac, ntrac_d, nnodes, nmpi_per_node, npart, qsize, qsize_d
  !-----------------
#ifdef CAM
  use time_mod, only : smooth, phys_tscale
#else
  use time_mod, only : tstep, ndays,nmax, nendstep,secpday, smooth, secphr, phys_tscale
#endif
  !-----------------
  use parallel_mod, only : parallel_t,  iam, abortmp, &
       partitionfornodes, useframes, mpireal_t, mpilogical_t, mpiinteger_t, mpichar_t
  !-----------------
  use cg_mod, only : cg_no_debug
  !-----------------

  use interpolate_mod, only : vector_uvars, vector_vvars, max_vecvars, interpolate_analysis
#ifndef CAM
  use common_io_mod, only : &
       output_prefix,       &
       output_type,         &
       output_timeunits,    &
       output_start_time,   &
       output_end_time,     &
       output_frequency,    &
       output_dir,          &
       output_varnames1,    &
       output_varnames2,    &
       output_varnames3,    &
       output_varnames4,    &
       output_varnames5,    &
       max_output_variables,&
       max_output_streams,  &
       num_io_procs,       &
       io_stride,          &
       varname_len,         &
       infilenames,         &
       MAX_INFILES

#ifdef _PRIM
  use aquaplanet, only : cool_ampl, &
       cool_min,      &
       cool_max,      &
       qv_flag,       &
       qv_pert_flag,  &
       qv_pert_ampl,  &
       qv_pert_zmin,  &
       qv_pert_zmax,  &
       isrf_forc,     &
       h_dis,         &
       Cdrag,         &
       wstar,         &
       tsurf,         &
       qsurf,         &
       u0,            &
       zabsampl,      &
       zabsmid,       &
       zabsmin,       &
       noisef       
#endif
  use common_movie_mod, only : setvarnames


#endif
  use interpolate_mod, only : set_interp_parameter, get_interp_parameter
  use fvm_mod, only: fvm_ideal_test, ideal_test_off, ideal_test_analytical_departure, ideal_test_analytical_winds
  use fvm_mod, only: fvm_test_type, ideal_test_boomerang, ideal_test_solidbody

!=======================================================================================================!                                 
!	Adding for SW DG										!
!=======================================================================================================!
#ifdef _SWDG  
  ! ------------------------
  use dg_flux_mod, only: riemanntype
  ! ------------------------
  use dg_tests_mod, only : alpha_dg, alphatype    
  ! ------------------------  
  use dg_sweq_mod, only: stage_rk	  
  ! ------------------------
  use physical_constants, only: dd_pi  
  ! ------------------------	
#endif

  character(len=32) :: tracer_transport_method = 'cslam_fvm'
  character(len=32) :: cslam_ideal_test = 'off'
  character(len=32) :: cslam_test_type = 'boomerang'
  type (parallel_t), intent(in) ::  par
  ! ============================================
  ! Namelists
  ! ============================================
  namelist /ctl_nl/ partmethod,          &       ! Mesh partitioning method (METIS)
                    qsize,               &       ! number of tracers
                    ntrac,               &       ! number of tracers
                    nthreads,            &       ! Number of threads per process
                    ne,                  &       ! element resolution factor
                    test_case,           &       ! test case
                    nmax,                &       ! number of steps
                    tstep,               &       ! time step
                    ndays,               &       ! number of days to simulate
                    test_cfldep,         &       ! test shape of departure grid cell and cfl number
                    tracer_transport_method,&
                    cslam_ideal_test,       &
                    cslam_test_type
                   
  namelist /analysis_nl/                 &
                    output_prefix,       &
                    output_timeunits,    &
                    output_start_time,   &
                    output_end_time,     &
                    output_frequency,    &
                    output_dir,          &
                    output_type,         &
                    output_varnames1
  
  integer           :: interp_nlat, interp_nlon    
!   integer :: interp_nlat, interp_nlon, interp_gridtype, interp_type
  namelist /analysis_nl/                 &
         interp_nlat,         &
         interp_nlon      
                
  integer           :: i, ierr
!=======================================================================================================!
  ! ==========================
  ! Set the default partmethod
  ! ==========================

  partmethod    = RECURSIVE 
  npart    = 1
  ndays         = 0
  nmax          = 12
  nthreads = 1

  ! =======================
  ! Read namelist variables
  ! =======================

  if (par%masterproc) then
    write(iulog,*)"reading ctl namelist..."
#if defined(OSF1) || defined(_BGL) || defined(_NAMELIST_FROM_FILE)
    open(unit=7,file="input.nl",status="OLD")
    read(unit=7,nml=ctl_nl)
#else
    read(*,nml=ctl_nl)
#endif

#ifndef _USEMETIS
    !=================================
    ! override the selected partition
    ! method and set it to SFCURVE 
    !=================================
    partmethod = sfcurve
#endif

    ! ================================================
    ! if ndays is defined in the namelist, then 
    ! moviefreq and restartfreq are interpreted to be in units of days. 
    ! Both must be converted to numbers of steps.
    ! ================================================
    if (tstep <= 0) then
      call abortmp('tstep must be > 0')
    end if
    if (ndays .gt. 0) then
      nmax = ndays * (secpday/tstep)
    end if
    nEndStep = nmax 

    output_prefix = ""
    output_start_time=0
    output_end_time =0
    output_frequency=0
    output_timeunits=0
! Default is to write all variables at time 0 and time nendstep.
    output_end_time(1) = -1
    output_frequency(1) = nendstep
    output_dir = "./results/"
    output_varnames1='     '
    io_stride=0
    num_io_procs=0
    output_type = 'pnetcdf'

    write(iulog,*)"reading analysis namelist..."
       
#if defined(OSF1) || defined(_BGL) || defined(_NAMELIST_FROM_FILE)
    read(unit=7,nml=analysis_nl)
#else
    read(*,nml=analysis_nl)
#endif
    if (io_stride .eq.0 .and. num_io_procs .eq.0) then
       ! user did not set anything
       io_stride=1
       num_io_procs = par%nprocs
    else if (io_stride.eq.0) then  ! user set num_io_procs
       io_stride = max(1,par%nprocs/num_io_procs)
    else if (num_io_procs .eq. 0) then   ! user set io_stride
       num_io_procs=max(1,par%nprocs/io_stride)
    else ! user set both parameters
    endif
    ! sanity check
    if(num_io_procs*io_stride>par%nprocs) then
       if (io_stride > par%nprocs) io_stride=par%nprocs
       num_io_procs=par%nprocs/io_stride
    end if
    if(output_varnames1(1).eq.'     ') then
      if( runtype>=0) 	then
        call setvarnames(output_varnames1)
      else ! interpolation mode
        output_varnames1='all'
      end if
    end if
    do i=1,max_output_streams
      if(output_frequency(i)>0 .and. output_end_time(i)==0) output_end_time(i)=-1
      if(output_timeunits(i).eq.1) then  ! per_day
        output_frequency(i) = output_frequency(i)*(secpday/tstep)
        output_start_time(i)= output_start_time(i)*(secpday/tstep)
        output_end_time(i)  = output_end_time(i)*(secpday/tstep)
      else if(output_timeunits(i).eq.2) then  ! per_hour
        output_frequency(i) = output_frequency(i)*(secphr/tstep)
        output_start_time(i)= output_start_time(i)*(secphr/tstep)
        output_end_time(i)  = output_end_time(i)*(secphr/tstep)
      end if
      if(output_end_time(i)<0) then
        output_end_time(i)=nEndStep
      endif
    end do
#if defined(OSF1) || defined(_BGL) || defined(_NAMELIST_FROM_FILE)
    close(unit=7)
#endif
  end if  
!=======================================================================================================!
    	
  call MPI_barrier(par%comm,ierr)
  npart  = par%nprocs 
  ! =====================================
  !  Spread the namelist stuff around 
  ! =====================================

  call MPI_bcast(PARTMETHOD ,1,MPIinteger_t,par%root,par%comm,ierr) 
  call MPI_bcast(test_case    ,MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
  call MPI_bcast(tasknum ,1,MPIinteger_t,par%root,par%comm,ierr)

  call MPI_bcast( NE        ,1,MPIinteger_t,par%root,par%comm,ierr)
  
  call MPI_bcast(qsize     ,1,MPIinteger_t,par%root,par%comm,ierr)
  
  call MPI_bcast(ntrac     ,1,MPIinteger_t,par%root,par%comm,ierr)

  call MPI_bcast(tstep     ,1,MPIreal_t   ,par%root,par%comm,ierr) 
  call MPI_bcast(nmax      ,1,MPIinteger_t,par%root,par%comm,ierr)
  
  call MPI_bcast(NTHREADS  ,1,MPIinteger_t,par%root,par%comm,ierr)
  call MPI_bcast(ndays     ,1,MPIinteger_t,par%root,par%comm,ierr)
  call MPI_bcast(test_cfldep,1,MPIlogical_t,par%root,par%comm,ierr)
  
  nEndStep = nmax


  call MPI_bcast(output_prefix,MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
  call MPI_bcast(output_timeunits ,max_output_streams,MPIinteger_t,par%root,par%comm,ierr)
  call MPI_bcast(output_start_time ,max_output_streams,MPIinteger_t,par%root,par%comm,ierr)
  call MPI_bcast(output_end_time ,max_output_streams,MPIinteger_t,par%root,par%comm,ierr)
  call MPI_bcast(output_frequency ,max_output_streams,MPIinteger_t,par%root,par%comm,ierr)

  call MPI_bcast(output_dir   ,MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
  call MPI_bcast(output_varnames1 ,varname_len*max_output_variables,MPIChar_t,par%root,par%comm,ierr)
  call MPI_bcast(output_type , 9,MPIChar_t,par%root,par%comm,ierr)
  call MPI_bcast(cubed_sphere_map,1,MPIinteger_t ,par%root,par%comm,ierr)

  call MPI_bcast(interpolate_analysis, 7,MPIlogical_t,par%root,par%comm,ierr)
  call MPI_bcast(interp_nlat , 1,MPIinteger_t,par%root,par%comm,ierr)
  call MPI_bcast(interp_nlon , 1,MPIinteger_t,par%root,par%comm,ierr)
  call MPI_bcast(io_stride , 1,MPIinteger_t,par%root,par%comm,ierr)
  call MPI_bcast(num_io_procs , 1,MPIinteger_t,par%root,par%comm,ierr)

! Set and broadcast tracer transport type
    if (trim(tracer_transport_method) == 'se_gll') then
      tracer_transport_type = TRACERTRANSPORT_SE_GLL
      tracer_grid_type = TRACER_GRIDTYPE_GLL
      if (ntrac>0) then
         call abortmp('user specified ntrac should only be > 0 when tracer_transport_type is fvm')
      end if
    else if (trim(tracer_transport_method) == 'cslam_fvm') then
      tracer_transport_type = TRACERTRANSPORT_LAGRANGIAN_FVM
      tracer_grid_type = TRACER_GRIDTYPE_FVM
      if (qsize>0) then
         call abortmp('user specified qsize should only be > 0 when tracer_transport_type is se_gll')
      end if
    else if (trim(tracer_transport_method) == 'flux_form_cslam_fvm') then
      tracer_transport_type = TRACERTRANSPORT_FLUXFORM_FVM
      tracer_grid_type = TRACER_GRIDTYPE_FVM
      if (qsize>0) then
         call abortmp('user specified qsize should only be > 0 when tracer_transport_type is se_gll')
      end if
    else
      call abortmp('Unknown tracer transport method: '//trim(tracer_transport_method))
    end if
    call MPI_bcast(tracer_transport_type,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(tracer_grid_type,1,MPIinteger_t,par%root,par%comm,ierr)
! Set and broadcast CSLAM test options
    if (trim(cslam_ideal_test) == 'off') then
      fvm_ideal_test = IDEAL_TEST_OFF
    else if (trim(cslam_ideal_test) == 'analytical_departure') then
      fvm_ideal_test = IDEAL_TEST_ANALYTICAL_DEPARTURE
    else if (trim(cslam_ideal_test) == 'analytical_winds') then
      fvm_ideal_test = IDEAL_TEST_ANALYTICAL_WINDS
    else
      call abortmp('Unknown ideal_cslam_test: '//trim(cslam_ideal_test))
    end if
    if (trim(cslam_test_type) == 'boomerang') then
      fvm_test_type = IDEAL_TEST_BOOMERANG
    else if (trim(cslam_test_type) == 'solidbody') then
      fvm_test_type = IDEAL_TEST_SOLIDBODY
    else
      call abortmp('Unknown cslam test type: '//trim(cslam_test_type))
    end if
    call MPI_bcast(fvm_ideal_test,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(fvm_test_type,1,MPIinteger_t,par%root,par%comm,ierr)

  ! set map
  if (cubed_sphere_map<0) then
     cubed_sphere_map=0  ! default is equi-angle gnomonic
     if (ne.eq.0) cubed_sphere_map=2  ! element_local for var-res grids
  endif
  if (par%masterproc) write (iulog,*) "Reference element projection: cubed_sphere_map=",cubed_sphere_map


  if(any(interpolate_analysis)) then
     if (interp_nlon==0 .or. interp_nlat==0) then
        ! compute interpolation grid based on number of points around equator
        call set_interp_parameter('auto',4*ne*(np-1))
        
        interp_nlat = get_interp_parameter('nlat')
        interp_nlon = get_interp_parameter('nlon')
        
        call MPI_bcast(interp_nlat , 1,MPIinteger_t,par%root,par%comm,ierr)
        call MPI_bcast(interp_nlon , 1,MPIinteger_t,par%root,par%comm,ierr)
     else
        call set_interp_parameter('nlat',interp_nlat)
        call set_interp_parameter('nlon',interp_nlon)
     endif
  endif
  
  if (par%masterproc) then
     write(iulog,*)"done reading namelist..."
  
     write(iulog,*)"readnl: test_case     = ",TRIM(test_case)
     write(iulog,*)"readnl: ndays         = ",ndays
     write(iulog,*)"readnl: nmax          = ",nmax

     write(iulog,*)"readnl: qsize,qsize_d = ",qsize,qsize_d
     if (qsize>qsize_d) then
        call abortmp('user specified qsize > qsize_d parameter in dimensions_mod.F90')
     endif
     write(iulog,*)"readnl: NThreads      = ",NTHREADS
     
     write(iulog,*)"readnl: ne,np,nc      = ",NE,np,nc
     write(iulog,*)"readnl: ntrac, ntrac_d         = ",ntrac, ntrac_d
     
     ! Write CSLAM namelist values
     select case (tracer_transport_type)
     case (TRACERTRANSPORT_SE_GLL)
        write(iulog, *) 'Eulerian tracer advection on GLL grid'
     case (TRACERTRANSPORT_SEMILAGRANG_GLL)
        write(iulog, *) 'Classic semi-Lagrangian tracer advection on GLL grid'
     case (TRACERTRANSPORT_LAGRANGIAN_FVM)
        write(iulog, *) 'CSLAM tracer advection on FVM grid'
     case (TRACERTRANSPORT_FLUXFORM_FVM)
        write(iulog, *) 'Flux-form CSLAM tracer advection on FVM grid'
     case (TRACERTRANSPORT_SPELT_FVM)
        write(iulog, *) 'Spelt tracer advection on FVM grid'
     end select
     if (fvm_ideal_test /= IDEAL_TEST_OFF) then
        select case (fvm_test_type)
        case (IDEAL_TEST_BOOMERANG)
           write(iulog, *) 'Running boomerang CSLAM test'
        case (IDEAL_TEST_SOLIDBODY)
           write(iulog, *) 'Running solid body CSLAM test'
        end select
        select case (fvm_ideal_test)
        case (IDEAL_TEST_ANALYTICAL_DEPARTURE)
           write(iulog, *) 'Using analytical departure points for CSLAM test'
        case (IDEAL_TEST_ANALYTICAL_WINDS)
           write(iulog, *) 'Using analytical winds for CSLAM test'
        end select
     end if
     
     
     write(iulog,*)"readnl: partmethod    = ",PARTMETHOD
     write(iulog,*)'readnl: nmpi_per_node = ',nmpi_per_node
     write(iulog,*)'readnl: multilevel    = ',multilevel
     write(iulog,*)'readnl: useframes     = ',useframes
     write(iulog,*)'readnl: nnodes        = ',nnodes
     write(iulog,*)'readnl: npart         = ',npart
     

     write(iulog,*)"readnl: tstep         = ",tstep

     write(iulog,*)"readnl: test_cfldep         = ",test_cfldep


     write(iulog,*)"  analysis: output_prefix = ",TRIM(output_prefix)
     write(iulog,*)"  analysis: io_stride = ",io_stride
     write(iulog,*)"  analysis: num_io_procs = ",num_io_procs
     
     write(iulog,*)" analysis interpolation = ", interpolate_analysis
     if(any(interpolate_analysis)) then
        write(iulog,*)" analysis interp nlat = ",interp_nlat
        write(iulog,*)" analysis interp nlon = ",interp_nlon
     end if
  end if

end subroutine fvm_readnl


end module fvm_init_mod
