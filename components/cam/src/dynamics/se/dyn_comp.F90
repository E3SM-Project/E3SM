Module dyn_comp

  use shr_kind_mod, only: r8 => shr_kind_r8
  use domain_mod, only : domain1d_t
  use element_mod, only : element_t
  use time_mod, only : TimeLevel_t, se_nsplit=>nsplit
  use hybvcoord_mod, only : hvcoord_t
  use hybrid_mod, only : hybrid_t
  use thread_mod, only: nthreads, vthreads, omp_get_max_threads, omp_get_thread_num
  use perf_mod, only: t_startf, t_stopf
  use cam_logfile, only : iulog
  use time_manager, only: is_first_step
  use spmd_utils,  only : iam, npes_cam => npes
  use pio,         only: file_desc_t

  implicit none
  private
  save


  ! PUBLIC MEMBER FUNCTIONS:
  public dyn_init1, dyn_init2, dyn_run

  ! PUBLIC DATA MEMBERS:
  public dyn_import_t, dyn_export_t


  type (TimeLevel_t)   , public :: TimeLevel     ! main time level struct (used by tracers)

  type dyn_import_t
     type (element_t), pointer :: elem(:) => null()
  end type dyn_import_t

  type dyn_export_t
     type (element_t), pointer :: elem(:) => null()
  end type dyn_export_t
  type (hvcoord_t), public  :: hvcoord
  integer, parameter  ::  DYN_RUN_SUCCESS           = 0
  integer, parameter  ::  DYN_RUN_FAILURE           = -1

  ! !DESCRIPTION: This module implements the SE Dynamical Core as
  !               an ESMF gridded component.  It is specific to SE
  !               and does not use ESMF.
  !
  ! \paragraph{Overview}
  !
  !   This module contains an ESMF wrapper for the SE
  !   Dynamical Core used in the Community Atmospheric Model. 
  !
  ! !REVISION HISTORY:
  !
  !  JPE  06.05.31:  created
  !
  !----------------------------------------------------------------------

  ! Enumeration of DYNAMICS_IN_COUPLINGS


  logical, parameter         :: DEBUG = .true.

  real(r8), parameter        :: ONE    = 1.0_r8

  character(*), parameter, public :: MODULE_NAME = "dyn_comp"
  character(*), parameter, public :: VERSION     = "$Id$" 
  type (domain1d_t), pointer, public :: dom_mt(:) => null()

  ! Frontogenesis indices
  integer, public :: frontgf_idx = -1
  integer, public :: frontga_idx = -1

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dyn_init1(fh, NLFileName, dyn_in, dyn_out)

  ! Initialize the dynamical core

    use pio,                 only: file_desc_t
    use hycoef,              only: hycoef_init
    use ref_pres,            only: ref_pres_init

    use pmgrid,              only: dyndecomp_set
    use dyn_grid,            only: dyn_grid_init, elem, get_dyn_grid_parm,&
                                   set_horiz_grid_cnt_d, define_cam_grids
    use rgrid,               only: fullgrid
    use spmd_utils,          only: mpi_integer, mpicom, mpi_logical
    use spmd_dyn,            only: spmd_readnl
    use native_mapping,      only: create_native_mapping_files, native_mapping_readnl
    use time_manager,        only: get_nstep, dtime

    use dimensions_mod,   only: globaluniquecols, nelem, nelemd, nelemdmax
    use prim_driver_mod,  only: prim_init1
    use parallel_mod,     only: par, initmp
    use namelist_mod,     only: readnl
    use control_mod,      only: runtype, qsplit, rsplit
    use time_mod,         only: tstep
    use phys_control,     only: use_gw_front
    use physics_buffer,   only: pbuf_add_field, dtype_r8
    use ppgrid,           only: pcols, pver
    use cam_abortutils,   only : endrun

    ! PARAMETERS:
    type(file_desc_t),   intent(in)  :: fh       ! PIO file handle for initial or restart file
    character(len=*),    intent(in)  :: NLFileName
    type (dyn_import_t), intent(OUT) :: dyn_in
    type (dyn_export_t), intent(OUT) :: dyn_out

    integer :: neltmp(3)
    integer :: npes_se

    !----------------------------------------------------------------------

    if (use_gw_front) then
       call pbuf_add_field("FRONTGF", "global", dtype_r8, (/pcols,pver/), &
            frontgf_idx)
       call pbuf_add_field("FRONTGA", "global", dtype_r8, (/pcols,pver/), &
            frontga_idx)
    end if

    ! Initialize dynamics grid
    call dyn_grid_init()

    ! Read in the number of tasks to be assigned to SE (needed by initmp)
    call spmd_readnl(NLFileName, npes_se)
    ! Initialize the SE structure that holds the MPI decomposition information
    par=initmp(npes_se)

    ! Read the SE specific part of the namelist
    call readnl(par, NLFileName)

    ! override the setting in the SE namelist, it's redundant anyway
    if (.not. is_first_step()) runtype = 1

    ! Initialize hybrid coordinate arrays.
    call hycoef_init(fh)

    ! Initialize physics grid reference pressures (needed by initialize_radbuffer)
    call ref_pres_init()

    ! legacy reduced grid code -- should be removed
    fullgrid=.true.

#ifdef _OPENMP    
!   Set by driver
    nthreads = omp_get_max_threads()
    if(par%masterproc) then
       write(iulog,*) " "
       write(iulog,*) "dyn_init1: number of OpenMP threads = ", nthreads
       write(iulog,*) " "
    endif
#ifdef COLUMN_OPENMP
    call omp_set_nested(.true.)
    if (vthreads > nthreads .or. vthreads < 1) &
         call endrun('Error: vthreads<1 or vthreads > NTHRDS_ATM')
    nthreads = nthreads / vthreads
    if(par%masterproc) then
       write(iulog,*) " "
       write(iulog,*) "dyn_init1: using OpenMP across and within elements"
       write(iulog,*) "dyn_init1: nthreads=",nthreads,"vthreads=",vthreads
       write(iulog,*) " "
    endif
#else
    if (vthreads>1) &
         call endrun('Error: vthreads>1 requires -DCOLUMN_OPENMP')
#endif
#else
    nthreads = 1
    if(par%masterproc) then
       write(iulog,*) " "
       write(iulog,*) "dyn_init1: openmp not activated"
       write(iulog,*) " "
    endif
#endif
    if(iam < par%nprocs) then
       call prim_init1(elem,par,dom_mt,TimeLevel)

       dyn_in%elem => elem
       dyn_out%elem => elem
    
       call set_horiz_grid_cnt_d(GlobalUniqueCols)

       neltmp(1) = nelemdmax
       neltmp(2) = nelem
       neltmp(3) = get_dyn_grid_parm('plon')
    else
       nelemd = 0
       neltmp(1) = 0
       neltmp(2) = 0
       neltmp(3) = 0
    endif

    dyndecomp_set = .true.

    if (par%nprocs .lt. npes_cam) then
! Broadcast quantities to auxiliary processes
#ifdef SPMD
       call mpibcast(neltmp, 3, mpi_integer, 0, mpicom)
#endif
       if (iam .ge. par%nprocs) then
          nelemdmax = neltmp(1)
          nelem     = neltmp(2)
          call set_horiz_grid_cnt_d(neltmp(3))
       endif
    endif


    !
    ! This subroutine creates mapping files using SE basis functions if requested
    !
    call native_mapping_readnl(NLFileName)
    call create_native_mapping_files( par, elem,'native')
    call create_native_mapping_files( par, elem,'bilin')

    ! Dynamics timestep
    !
    !  Note: dtime = progress made in one timestep.  value in namelist
    !        dtime = the frequency at which physics is called
    !        tstep = the dynamics timestep:  
    !

    if (rsplit==0) then
       ! non-lagrangian code
       tstep = dtime/real(se_nsplit*qsplit,r8)
       TimeLevel%nstep = get_nstep()*se_nsplit*qsplit
   else
      ! lagrangian code
       tstep = dtime/real(se_nsplit*qsplit*rsplit,r8)
       TimeLevel%nstep = get_nstep()*se_nsplit*qsplit*rsplit
    endif

    ! Define the CAM grids (this has to be after dycore spinup).
    ! Physics-grid will be defined later by phys_grid_init
    call define_cam_grids()

  end subroutine dyn_init1


  subroutine dyn_init2(dyn_in)
    use dimensions_mod,   only: nlev, nelemd
    use prim_driver_mod,  only: prim_init2
    use prim_si_ref_mod,  only: prim_set_mass
    use hybrid_mod,       only: hybrid_create
    use hycoef,           only: hyam, hybm, hyai, hybi, ps0
    use parallel_mod,     only: par
    use time_mod,         only: time_at
    use control_mod,      only: moisture, runtype
    use cam_control_mod,  only: aqua_planet, ideal_phys, adiabatic
    use comsrf,           only: landm, sgh, sgh30
    use nctopo_util_mod,  only: nctopo_util_driver
    use cam_instance,     only: inst_index

    type (dyn_import_t), intent(inout) :: dyn_in

    type(element_t),    pointer :: elem(:)

    integer :: ithr, nets, nete, ie, k
    real(r8), parameter :: Tinit=300.0_r8
    real(r8) :: dyn_ps0
    type(hybrid_t) :: hybrid

    elem  => dyn_in%elem

    dyn_ps0=ps0
    hvcoord%hyam=hyam
    hvcoord%hyai=hyai
    hvcoord%hybm=hybm
    hvcoord%hybi=hybi
    hvcoord%ps0=dyn_ps0  
    do k=1,nlev
       hvcoord%hybd(k) = hvcoord%hybi(k+1) - hvcoord%hybi(k)
    end do
    if(iam < par%nprocs) then

#ifdef HORIZ_OPENMP
       if (iam==0) write (iulog,*) "dyn_init2: nthreads=",nthreads,&
                                   "max_threads=",omp_get_max_threads()
       !$OMP PARALLEL NUM_THREADS(nthreads), DEFAULT(SHARED), PRIVATE(ie,ithr,nets,nete,hybrid)
#endif
#ifdef COLUMN_OPENMP
       call omp_set_num_threads(vthreads)
#endif
       ithr=omp_get_thread_num()
       nets=dom_mt(ithr)%start
       nete=dom_mt(ithr)%end
       hybrid = hybrid_create(par,ithr,NThreads)

       moisture='moist'

       if(adiabatic) then
          moisture='dry'
          if(runtype == 0) then
             do ie=nets,nete
                elem(ie)%state%q(:,:,:,:)=0.0_r8
                elem(ie)%derived%fq(:,:,:,:,:)=0.0_r8
             end do
          end if
       else if(ideal_phys) then
          moisture='dry'
          if(runtype == 0) then
             do ie=nets,nete
                elem(ie)%state%ps_v(:,:,:) =dyn_ps0

                elem(ie)%state%phis(:,:)=0.0_r8

                elem(ie)%state%T(:,:,:,:) =Tinit

                elem(ie)%state%v(:,:,:,:,:) =0.0_r8

                elem(ie)%state%q(:,:,:,:)=0.0_r8

             end do
          end if
       else if(aqua_planet .and. runtype==0)  then
          do ie=nets,nete
             elem(ie)%state%phis(:,:)=0.0_r8
          end do
          if(allocated(landm)) landm=0.0_r8
          if(allocated(sgh)) sgh=0.0_r8
          if(allocated(sgh30)) sgh30=0.0_r8
       end if

       do ie=nets,nete
          elem(ie)%derived%FM=0.0_r8
          elem(ie)%derived%FT=0.0_r8
          elem(ie)%derived%FQ=0.0_r8
       end do

       ! scale PS to achieve prescribed dry mass
       if (runtype == 0) then
          ! new run, scale mass to value given in namelist, if needed
          call prim_set_mass(elem, TimeLevel,hybrid,hvcoord,nets,nete)
       endif
       call prim_init2(elem,hybrid,nets,nete, TimeLevel, hvcoord)
       !
       ! This subroutine is used to create nc_topo files, if requested
       ! 
       call nctopo_util_driver(elem,hybrid,nets,nete)
#ifdef HORIZ_OPENMP
       !$OMP END PARALLEL 
#endif
    end if

    if (inst_index == 1) then
       call write_grid_mapping(par, elem)
    end if

  end subroutine dyn_init2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-----------------------------------------------------------------------
  !BOP
  ! !ROUTINE:  RUN --- Driver for the 
  !
  ! !INTERFACE:
  subroutine dyn_run( dyn_state, rc )

    ! !USES:
    use parallel_mod,     only : par
    use prim_driver_mod,  only: prim_run_subcycle
    use dimensions_mod,   only : nlev
    use time_mod,         only: tstep
    use hybrid_mod,       only: hybrid_create
!    use perf_mod, only : t_startf, t_stopf
    implicit none


    type (dyn_export_t), intent(inout)       :: dyn_state   !  container
    type(hybrid_t) :: hybrid

    integer, intent(out)               :: rc      ! Return code
    integer ::  n
    integer :: nets, nete, ithr
    integer :: ie

    ! !DESCRIPTION:
    !
    if(iam < par%nprocs) then
#ifdef HORIZ_OPENMP
       !if (iam==0) write (iulog,*) "dyn_run: nthreads=",nthreads,&
       !                            "max_threads=",omp_get_max_threads()
       !$OMP PARALLEL NUM_THREADS(nthreads), DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid,n)
#endif
#ifdef COLUMN_OPENMP
       ! nested threads
       call omp_set_num_threads(vthreads)
#endif
       ithr=omp_get_thread_num()
       nets=dom_mt(ithr)%start
       nete=dom_mt(ithr)%end
       hybrid = hybrid_create(par,ithr,NThreads)

       do n=1,se_nsplit
          ! forward-in-time RK, with subcycling
          call t_startf("prim_run_sybcycle")
          call prim_run_subcycle(dyn_state%elem,hybrid,nets,nete,&
               tstep, TimeLevel, hvcoord, n)
          call t_stopf("prim_run_sybcycle")
       end do


#ifdef HORIZ_OPENMP
       !$OMP END PARALLEL
#endif
    end if
    rc = DYN_RUN_SUCCESS

    !EOC
  end subroutine dyn_run
  !-----------------------------------------------------------------------



  subroutine write_grid_mapping(par, elem)
    use parallel_mod,     only: parallel_t
    use element_mod, only : element_t
    use cam_pio_utils, only : cam_pio_createfile, pio_subsystem
    use pio, only : file_desc_t, pio_def_dim, var_desc_t, pio_int, pio_def_var, &
         pio_enddef, pio_closefile, pio_initdecomp, io_desc_t, pio_write_darray, &
         pio_freedecomp, pio_setdebuglevel
    use dimensions_mod, only : np, nelem, nelemd
    use dof_mod, only : createmetadata

    type(parallel_t) :: par
    type(element_t) :: elem(:)
    type(file_desc_t) :: nc
    type(var_desc_t) :: vid
    type(io_desc_t) :: iodesc
    integer :: dim1, dim2, ierr, i, j, ie, cc, base, ii, jj
    integer, parameter :: npm12 = (np-1)*(np-1)
    integer :: subelement_corners(npm12*nelemd,4)
    integer :: dof(npm12*nelemd*4)


    ! Create a CS grid mapping file for postprocessing tools

       ! write meta data for physics on GLL nodes
       call cam_pio_createfile(nc, 'SEMapping.nc', 0)
   
       ierr = pio_def_dim(nc, 'ncenters', npm12*nelem, dim1)
       ierr = pio_def_dim(nc, 'ncorners', 4, dim2)
       ierr = pio_def_var(nc, 'element_corners', PIO_INT, (/dim1,dim2/),vid)
    
       ierr = pio_enddef(nc)
       call createmetadata(par, elem, subelement_corners)

       jj=0
       do cc=0,3
          do ie=1,nelemd
             base = ((elem(ie)%globalid-1)+cc*nelem)*npm12
             ii=0
             do j=1,np-1
                do i=1,np-1
                   ii=ii+1
                   jj=jj+1
                   dof(jj) = base+ii
                end do
             end do
          end do
       end do

       call pio_initdecomp(pio_subsystem, pio_int, (/nelem*npm12,4/), dof, iodesc)

       call pio_write_darray(nc, vid, iodesc, reshape(subelement_corners,(/nelemd*npm12*4/)), ierr)
       
       call pio_freedecomp(nc, iodesc)
       
       call pio_closefile(nc)

  end subroutine write_grid_mapping

end module dyn_comp



