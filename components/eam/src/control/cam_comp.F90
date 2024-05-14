module cam_comp
!-----------------------------------------------------------------------
!
! Purpose:  The CAM Community Atmosphere Model component. Interfaces with 
!           a merged surface field that is provided outside of this module. 
!           This is the atmosphere only component. It can interface with a 
!           host of surface components.
!
!-----------------------------------------------------------------------
   use shr_kind_mod,      only: r8 => SHR_KIND_R8, cl=>SHR_KIND_CL, cs=>SHR_KIND_CS
   use pmgrid,            only: plat, plev
   use spmd_utils,        only: masterproc
   use cam_abortutils,        only: endrun
   use camsrfexch,        only: cam_out_t, cam_in_t     
   use shr_sys_mod,       only: shr_sys_flush
   use physics_types,     only: physics_state, physics_tend
   use cam_control_mod,   only: nsrest, print_step_cost, obliqr, lambm0, mvelpp, eccen
   use dyn_comp,          only: dyn_import_t, dyn_export_t
   use ppgrid,            only: begchunk, endchunk
   use perf_mod
   use cam_logfile,       only: iulog
   use physics_buffer,            only: physics_buffer_desc

   implicit none
   private
   save
   !
   ! Public access methods
   !
   public cam_init      ! First phase of CAM initialization
   public cam_run1      ! CAM run method phase 1
   public cam_run2      ! CAM run method phase 2
   public cam_run3      ! CAM run method phase 3
   public cam_run4      ! CAM run method phase 4
   public cam_final     ! CAM Finalization
   !
   ! Private module data
   !
#if ( defined SPMD )
   real(r8) :: cam_time_beg              ! Cam init begin timestamp
   real(r8) :: cam_time_end              ! Cam finalize end timestamp
   real(r8) :: stepon_time_beg = -1.0_r8 ! Stepon (run1) begin timestamp
   real(r8) :: stepon_time_end = -1.0_r8 ! Stepon (run4) end timestamp
   integer  :: nstep_beg = -1            ! nstep at beginning of run
#else
   integer  :: mpicom = 0
#endif

  real(r8) :: gw(plat)           ! Gaussian weights
  real(r8) :: etamid(plev)       ! vertical coords at midpoints
  real(r8) :: dtime              ! Time step for either physics or dynamics (set in dynamics init)

  type(dyn_import_t) :: dyn_in   ! Dynamics import container
  type(dyn_export_t) :: dyn_out  ! Dynamics export container

  type(physics_state), pointer :: phys_state(:) => null()
  type(physics_tend ), pointer :: phys_tend(:) => null()
  type(physics_buffer_desc), pointer :: pbuf2d(:,:) => null()

  real(r8) :: wcstart, wcend     ! wallclock timestamp at start, end of timestep
  real(r8) :: usrstart, usrend   ! user timestamp at start, end of timestep
  real(r8) :: sysstart, sysend   ! sys timestamp at start, end of timestep

!-----------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!
subroutine cam_init( cam_out, cam_in, mpicom_atm, &
                     start_ymd, start_tod, ref_ymd, ref_tod, stop_ymd, stop_tod, &
                     perpetual_run, perpetual_ymd, calendar)

   !-----------------------------------------------------------------------
   !
   ! Purpose:  CAM initialization.
   !
   !-----------------------------------------------------------------------

   use infnan,           only: nan, assignment(=)
   use history_defaults, only: bldfld
   use cam_initfiles,    only: cam_initfiles_open
   use inital,           only: cam_initial
   use cam_restart,      only: cam_read_restart
   use stepon,           only: stepon_init
   use physpkg,          only: phys_init
   
   use dycore,           only: dycore_is
#if (defined E3SM_SCM_REPLAY)
   use history_defaults, only: initialize_iop_history
#endif
!   use shr_orb_mod,      only: shr_orb_params
   use camsrfexch,       only: hub2atm_alloc, atm2hub_alloc
   use cam_history,      only: intht, init_masterlinkedlist
   use history_scam,     only: scm_intht
   use scamMod,          only: single_column
   use cam_pio_utils,    only: init_pio_subsystem
   use cam_instance,     only: inst_suffix
#if defined(CLDERA_PROFILING)
   use iso_c_binding, only: c_loc
   use cldera_interface_mod, only: cldera_add_partitioned_field, max_str_len, &
                                   cldera_set_field_part_extent, &
                                   cldera_set_field_part_data, &
                                   cldera_commit_all_fields,   &
                                   cldera_commit_field
   use physics_buffer,   only: physics_buffer_desc, col_type_grid, pbuf_get_index, &
                               pbuf_get_field_rank, pbuf_get_field_dims, pbuf_get_field, &
                               pbuf_get_field_name, pbuf_has_field, pbuf_get_field_persistence, &
                               persistence_global
   use ppgrid,           only: begchunk, endchunk, pcols, pver
   use phys_grid,        only: get_ncols_p, get_gcol_all_p, get_area_all_p
   use constituents,     only: pcnst, cnst_name
#endif

#if ( defined SPMD )   
   real(r8) :: mpi_wtime  ! External
#endif
   !-----------------------------------------------------------------------
   !
   ! Arguments
   !
   type(cam_out_t), pointer        :: cam_out(:)       ! Output from CAM to surface
   type(cam_in_t) , pointer        :: cam_in(:)        ! Merged input state to CAM
   integer            , intent(in) :: mpicom_atm       ! CAM MPI communicator
   integer            , intent(in) :: start_ymd        ! Start date (YYYYMMDD)
   integer            , intent(in) :: start_tod        ! Start time of day (sec)
   integer            , intent(in) :: ref_ymd          ! Reference date (YYYYMMDD)
   integer            , intent(in) :: ref_tod          ! Reference time of day (sec)
   integer            , intent(in) :: stop_ymd         ! Stop date (YYYYMMDD)
   integer            , intent(in) :: stop_tod         ! Stop time of day (sec)
   logical            , intent(in) :: perpetual_run    ! If in perpetual mode or not
   integer            , intent(in) :: perpetual_ymd    ! Perpetual date (YYYYMMDD)
   character(len=cs)  , intent(in) :: calendar         ! Calendar type
   !
   ! Local variables
   !
   integer :: dtime_cam        ! Time-step
   logical :: log_print        ! Flag to print out log information or not
   character(len=cs) :: filein ! Input namelist filename
#if defined(CLDERA_PROFILING)
   character(len=max_str_len) :: fname
   integer :: c, nfields, idx, rank, icmp, nparts, part_dim, ipart, fsize, ncols, icall, tag_loop
   integer :: nlcols,irank,part_alloc_size
   integer :: dims(3)
   integer, allocatable :: cols_gids(:)
   real(r8), allocatable :: cols_area(:)
   character(len=max_str_len) :: dimnames(3)
   logical :: in_pbuf, in_q
   real(r8), pointer :: field1d(:), field2d(:,:), field3d(:,:,:)
   type(physics_buffer_desc), pointer :: field_desc
   character(len=5) :: int_str
   character(len=4) :: diag(0:2) = (/'    ','_d1 ','_d2 '/)
   character(len=2) :: tagged_suffix(3) = (/'01', '02', '03'/)
#endif
   !-----------------------------------------------------------------------
   etamid = nan
   !
   ! Initialize CAM MPI communicator, number of processors, master processors 
   !	
#if ( defined SPMD )
   cam_time_beg = mpi_wtime()
#endif
   !
   ! Initialization needed for cam_history
   ! 
   call init_masterlinkedlist()
   !
   ! Set up spectral arrays
   !
   call trunc()
   !
   ! Determine input namelist filename
   !
   filein = "atm_in" // trim(inst_suffix)
   !
   ! Do appropriate dynamics and history initialization depending on whether initial, restart, or 
   ! branch.  On restart run intht need not be called because all the info is on restart dataset.
   !
   call init_pio_subsystem(filein)

   if ( nsrest == 0 )then

      call t_startf('cam_initfiles_open')
      call cam_initfiles_open()
      call t_stopf('cam_initfiles_open')

      call t_startf('cam_initial')
      call cam_initial(dyn_in, dyn_out, NLFileName=filein)
      call t_stopf('cam_initial')

      ! Allocate and setup surface exchange data
      call atm2hub_alloc(cam_out)
      call hub2atm_alloc(cam_in)

   else

      call t_startf('cam_read_restart')
      call cam_read_restart ( cam_in, cam_out, dyn_in, dyn_out, pbuf2d, stop_ymd, stop_tod, NLFileName=filein )
      call t_stopf('cam_read_restart')

     ! Commented out the hub2atm_alloc call as it overwrite cam_in, which is undesirable. The fields in cam_in are necessary for getting BFB restarts
	 ! There are no side effects of commenting out this call as this call allocates cam_in and cam_in allocation has already been done in cam_init
     !call hub2atm_alloc( cam_in )
#if (defined E3SM_SCM_REPLAY)
      call initialize_iop_history()
#endif
   end if

   call t_startf('phys_init')
   call phys_init( phys_state, phys_tend, pbuf2d,  cam_out )
   call t_stopf('phys_init')

   call bldfld ()       ! master field list (if branch, only does hash tables)

   !
   ! Setup the characteristics of the orbit
   ! (Based on the namelist parameters)
   !
   if (masterproc) then
      log_print = .true.
   else
      log_print = .false.
   end if

   call stepon_init( dyn_in, dyn_out ) ! dyn_out necessary?

   if (single_column) call scm_intht()
   call intht()

#if defined(CLDERA_PROFILING)
   call t_startf('cldera_add_fields')
   nparts = endchunk - begchunk + 1

   ! All fields are partitioned over cols index, which is the first
   part_dim = 1
   nlcols = 0
   do c = begchunk,endchunk
     nlcols = nlcols +  get_ncols_p(c)
   enddo
   dimnames(1) = "ncol"

   ! Col GIDs and area
   ! NOTE: .false. is to declare the field as a Copy of input data, rather than a view
   dims(1) = nlcols
   allocate(cols_gids(pcols))
   allocate(cols_area(pcols))
   part_alloc_size = pcols
   call cldera_add_partitioned_field ("col_gids",1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.,"int")
   call cldera_add_partitioned_field ("area",1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   do ipart = 1,nparts
     c = begchunk+ipart-1
     ncols = get_ncols_p(c)

     call cldera_set_field_part_extent("col_gids",ipart,ncols)
     call cldera_set_field_part_extent("area",ipart,ncols)
   enddo
   call cldera_commit_field("col_gids")
   call cldera_commit_field("area")
   do ipart = 1,nparts
     c = begchunk+ipart-1
     call get_gcol_all_p (c,pcols,cols_gids)
     call cldera_set_field_part_data("col_gids",ipart,cols_gids)
     call get_area_all_p(c,pcols,cols_area)
     call cldera_set_field_part_data("area",ipart,cols_area)
   enddo

   ! PBUF fields
   nfields = size(pbuf2d,1)
   do idx=1,nfields
     ! retrieve fields
     fname = pbuf_get_field_name(idx)

     ! We can only track persistent pbuf quantities
     if (pbuf_get_field_persistence(fname) .ne. persistence_global) then
       cycle
     endif

     ! retrieve field specs
     field_desc => pbuf2d(idx,begchunk)

     rank = pbuf_get_field_rank(idx)
     dims(:) = pbuf_get_field_dims(idx)
     dims(1) = nlcols

     do irank=2,rank
       if (dims(irank) .eq. plev) then
         dimnames(irank) = "lev"
       elseif (dims(irank) .eq. (plev+1)) then
         dimnames(irank) = "ilev"
       else
         write (int_str,"(I5)") dims(irank)
         dimnames(irank) = "dim"//trim(adjustl(int_str))
       endif
     enddo

     call cldera_add_partitioned_field(fname,2,dims,dimnames,nparts,part_dim,part_alloc_size)
     do ipart = 1,nparts
       c = begchunk+ipart-1 ! Chunk
       ncols = get_ncols_p(c)
       call cldera_set_field_part_extent(fname,ipart,ncols)
       if (rank .eq. 1) then
         call pbuf_get_field(pbuf2d, c, idx, field1d)
         call cldera_set_field_part_data(fname,ipart,field1d)
       elseif (rank .eq. 2) then
         call pbuf_get_field(pbuf2d, c, idx, field2d)
         call cldera_set_field_part_data(fname,ipart,field2d)
       else 
         call pbuf_get_field(pbuf2d, c, idx, field3d)
         call cldera_set_field_part_data(fname,ipart,field3d)
       endif
     enddo
   enddo

   ! TRACERS fields
   dims(2) = plev
   dimnames(2) = "lev"
   do idx=1,pcnst
     fname = cnst_name(idx)

     call cldera_add_partitioned_field(fname,2,dims,dimnames,nparts,part_dim,part_alloc_size)
     do ipart = 1,nparts
       c = begchunk+ipart-1 ! Chunk
       ncols = phys_state(c)%ncol
       field2d=>phys_state(c)%q(:,:,idx)
       call cldera_set_field_part_extent(fname,ipart,ncols)
       call cldera_set_field_part_data(fname,ipart,field2d)
     enddo
   enddo

   ! STATE fields
   dims(1) = nlcols
   dimnames(1) = 'ncol'

   !1d, horizontal
   call cldera_add_partitioned_field("lat",1,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("lon",1,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("ps",1,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("psdry",1,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("phis",1,dims,dimnames,nparts,part_dim,part_alloc_size)

   ! Last arg is view=false, since these fields are *not* views of EAM persistent data.
   call cldera_add_partitioned_field("AEROD_v", 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   call cldera_add_partitioned_field("AODALL", 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   call cldera_add_partitioned_field("ABSORB", 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   call cldera_add_partitioned_field("AODVIS", 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   call cldera_add_partitioned_field("AODABS", 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   call cldera_add_partitioned_field("aod"    , 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   call cldera_add_partitioned_field("aod_so2", 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   call cldera_add_partitioned_field("aod_ash", 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   call cldera_add_partitioned_field("aod_sulf",1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   do icall = 2,0,-1 ! profile climate calculation & two diags for now
      call cldera_add_partitioned_field("SOLIN"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSDS"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSNIRTOA"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSNRTOAC"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSNRTOAS"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSNT"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSNS"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSNTC"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSNSC"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSDSC"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSNTOA"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSUTOA"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSNTOAC"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSUTOAC"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("SOLS"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("SOLL"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("SOLSD"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("SOLLD"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSN200"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FSN200C"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("SWCF"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FLNT"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FLUT"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FLUTC"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FLNTC"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FLNS"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FLDSC"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FLNSC"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("LWCF"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FLN200"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FLN200C"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("FLDS"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("CLDTOT_ISCCP"//diag(icall), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)

   end do
   call cldera_add_partitioned_field("AODSO4", 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   call cldera_add_partitioned_field("BURDENSO4", 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   do tag_loop = 1,3 ! only three tags needed for now
      call cldera_add_partitioned_field("AODSO4"//tagged_suffix(tag_loop), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("BURDENSO4"//tagged_suffix(tag_loop), 1,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   end do

   ! 2d, mid points
   dims(2) = pver
   dimnames(2) = 'lev'
   call cldera_add_partitioned_field("T",2,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("u",2,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("v",2,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("s",2,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("omega",2,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("pmid",2,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("pmid_dry",2,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("pdel",2,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("pdel_dry",2,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("exner",2,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("zm",2,dims,dimnames,nparts,part_dim,part_alloc_size)

   call cldera_add_partitioned_field("NIHF"//diag(icall), 2,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   call cldera_add_partitioned_field("NIIMM"//diag(icall), 2,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   call cldera_add_partitioned_field("NIDEP"//diag(icall), 2,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   call cldera_add_partitioned_field("NIMEY"//diag(icall), 2,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)


   ! 2d, mid points (copy)
   do icall = 2,0,-1 ! profile climate calculation & two diags for now
      call cldera_add_partitioned_field('QRS'//diag(icall), 2,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field('QRSC'//diag(icall), 2,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("QRL"//diag(icall), 2,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
      call cldera_add_partitioned_field("QRLC"//diag(icall), 2,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   end do
   call cldera_add_partitioned_field("Mass_so4",2,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   do tag_loop = 1,3 ! only three tags needed for now
      call cldera_add_partitioned_field("Mass_so4"//tagged_suffix(tag_loop),2,dims,dimnames,nparts,part_dim,part_alloc_size,.false.)
   end do

   ! 2d, interfaces
   dims(2) = pver+1
   dimnames(2) = "ilev"
   call cldera_add_partitioned_field("pint",2,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("pint_dry",2,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("zi",2,dims,dimnames,nparts,part_dim,part_alloc_size)

   ! 1d, vertically integrated
   call cldera_add_partitioned_field("te",1,dims,dimnames,nparts,part_dim,part_alloc_size) ! total energy
   call cldera_add_partitioned_field("tw",1,dims,dimnames,nparts,part_dim,part_alloc_size) ! total water

   ! cam_in fields
   call cldera_add_partitioned_field("TREFHT", 1,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("QREFHT", 1,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("TS", 1,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("QFLX", 1,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("SHFLX", 1,dims,dimnames,nparts,part_dim,part_alloc_size)
   call cldera_add_partitioned_field("LHFLX", 1,dims,dimnames,nparts,part_dim,part_alloc_size)

   ! Set fields data
   do ipart = 1,nparts
     c = begchunk+ipart-1 ! Chunk
     ncols = phys_state(c)%ncol

     ! 1d (horiz) fields
     field1d => phys_state(c)%lat(:)
     call cldera_set_field_part_extent("lat",ipart,ncols)
     call cldera_set_field_part_data("lat",ipart,field1d)

     field1d => phys_state(c)%lon(:)
     call cldera_set_field_part_extent("lon",ipart,ncols)
     call cldera_set_field_part_data("lon",ipart,field1d)

     field1d => phys_state(c)%ps(:)
     call cldera_set_field_part_extent("ps",ipart,ncols)
     call cldera_set_field_part_data("ps",ipart,field1d)

     field1d => phys_state(c)%psdry(:)
     call cldera_set_field_part_extent("psdry",ipart,ncols)
     call cldera_set_field_part_data("psdry",ipart,field1d)

     field1d => phys_state(c)%phis(:)
     call cldera_set_field_part_extent("phis",ipart,ncols)
     call cldera_set_field_part_data("phis",ipart,field1d)

     ! 2d mid
     field2d => phys_state(c)%t(:,:)
     call cldera_set_field_part_extent("T",ipart,ncols)
     call cldera_set_field_part_data("T",ipart,field2d)

     field2d => phys_state(c)%u(:,:)
     call cldera_set_field_part_extent("u",ipart,ncols)
     call cldera_set_field_part_data("u",ipart,field2d)

     field2d => phys_state(c)%v(:,:)
     call cldera_set_field_part_extent("v",ipart,ncols)
     call cldera_set_field_part_data("v",ipart,field2d)

     field2d => phys_state(c)%s(:,:)
     call cldera_set_field_part_extent("s",ipart,ncols)
     call cldera_set_field_part_data("s",ipart,field2d)

     field2d => phys_state(c)%omega(:,:)
     call cldera_set_field_part_extent("omega",ipart,ncols)
     call cldera_set_field_part_data("omega",ipart,field2d)

     field2d => phys_state(c)%pmid(:,:)
     call cldera_set_field_part_extent("pmid",ipart,ncols)
     call cldera_set_field_part_data("pmid",ipart,field2d)

     field2d => phys_state(c)%pmiddry(:,:)
     call cldera_set_field_part_extent("pmid_dry",ipart,ncols)
     call cldera_set_field_part_data("pmid_dry",ipart,field2d)

     field2d => phys_state(c)%pdel(:,:)
     call cldera_set_field_part_extent("pdel",ipart,ncols)
     call cldera_set_field_part_data("pdel",ipart,field2d)

     field2d => phys_state(c)%pdeldry(:,:)
     call cldera_set_field_part_extent("pdel_dry",ipart,ncols)
     call cldera_set_field_part_data("pdel_dry",ipart,field2d)

     field2d => phys_state(c)%exner(:,:)
     call cldera_set_field_part_extent("exner",ipart,ncols)
     call cldera_set_field_part_data("exner",ipart,field2d)

     field2d => phys_state(c)%zm(:,:)
     call cldera_set_field_part_extent("zm",ipart,ncols)
     call cldera_set_field_part_data("zm",ipart,field2d)

     ! 2d int
     field2d => phys_state(c)%pint(:,:)
     call cldera_set_field_part_extent("pint",ipart,ncols)
     call cldera_set_field_part_data("pint",ipart,field2d)
     field2d => phys_state(c)%pintdry(:,:)
     call cldera_set_field_part_extent("pint_dry",ipart,ncols)
     call cldera_set_field_part_data("pint_dry",ipart,field2d)
     field2d => phys_state(c)%zi(:,:)
     call cldera_set_field_part_extent("zi",ipart,ncols)
     call cldera_set_field_part_data("zi",ipart,field2d)

     ! 1d, vertically integrated
     field1d => phys_state(c)%te_cur(:)
     call cldera_set_field_part_extent("te",ipart,ncols)
     call cldera_set_field_part_data("te",ipart,field1d)
     field1d => phys_state(c)%tw_cur(:)
     call cldera_set_field_part_extent("tw",ipart,ncols)
     call cldera_set_field_part_data("tw",ipart,field1d)

     ! cam_in fields
     field1d => cam_in(c)%tref(:)
     call cldera_set_field_part_extent("TREFHT",ipart,ncols)
     call cldera_set_field_part_data("TREFHT",ipart,field1d)
     field1d => cam_in(c)%qref(:)
     call cldera_set_field_part_extent("QREFHT",ipart,ncols)
     call cldera_set_field_part_data("QREFHT",ipart,field1d)
     field1d => cam_in(c)%ts(:)
     call cldera_set_field_part_extent("TS",ipart,ncols)
     call cldera_set_field_part_data("TS",ipart,field1d)
     field1d => cam_in(c)%cflx(:,1)
     call cldera_set_field_part_extent("QFLX",ipart,ncols)
     call cldera_set_field_part_data("QFLX",ipart,field1d)
     field1d => cam_in(c)%shf(:)
     call cldera_set_field_part_extent("SHFLX",ipart,ncols)
     call cldera_set_field_part_data("SHFLX",ipart,field1d)
     field1d => cam_in(c)%lhf(:)
     call cldera_set_field_part_extent("LHFLX",ipart,ncols)
     call cldera_set_field_part_data("LHFLX",ipart,field1d)

     ! Copied field
     call cldera_set_field_part_extent("AEROD_v", ipart,ncols)
     call cldera_set_field_part_extent("AODALL", ipart,ncols)
     call cldera_set_field_part_extent("ABSORB", ipart,ncols)
     call cldera_set_field_part_extent("AODVIS", ipart,ncols)
     call cldera_set_field_part_extent("AODABS", ipart,ncols)
     call cldera_set_field_part_extent("aod"    , ipart,ncols)
     call cldera_set_field_part_extent("aod_so2", ipart,ncols)
     call cldera_set_field_part_extent("aod_ash", ipart,ncols)
     call cldera_set_field_part_extent("aod_sulf",ipart,ncols)
     do icall = 2,0,-1 ! profile climate calculation & two diags for now
       call cldera_set_field_part_extent("SOLIN"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSDS"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSNIRTOA"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSNRTOAC"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSNRTOAS"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSNT"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSNS"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSNTC"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSNSC"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSDSC"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSNTOA"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSUTOA"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSNTOAC"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSUTOAC"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("SOLS"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("SOLL"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("SOLSD"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("SOLLD"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSN200"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FSN200C"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("SWCF"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FLNT"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FLUT"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FLUTC"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FLNTC"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FLNS"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FLDSC"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FLNSC"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("LWCF"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FLN200"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FLN200C"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("FLDS"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("QRS"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("QRSC"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("QRL"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("QRLC"//diag(icall), ipart,ncols)
       call cldera_set_field_part_extent("CLDTOT_ISCCP"//diag(icall), ipart,ncols)

     end do
     call cldera_set_field_part_extent("AODSO4", ipart,ncols)
     call cldera_set_field_part_extent("BURDENSO4", ipart,ncols)
     call cldera_set_field_part_extent("Mass_so4", ipart,ncols)
     call cldera_set_field_part_extent("NIHF", ipart,ncols)
     call cldera_set_field_part_extent("NIIMM", ipart,ncols)
     call cldera_set_field_part_extent("NIDEP", ipart,ncols)
     call cldera_set_field_part_extent("NIMEY", ipart,ncols)
     do tag_loop = 1,3 ! only three tags needed for now
       call cldera_set_field_part_extent("AODSO4"//tagged_suffix(tag_loop), ipart,ncols)
       call cldera_set_field_part_extent("BURDENSO4"//tagged_suffix(tag_loop), ipart,ncols)
       call cldera_set_field_part_extent("Mass_so4"//tagged_suffix(tag_loop), ipart,ncols)
     end do
   enddo

   call cldera_commit_all_fields()
   call t_stopf('cldera_add_fields')
#endif


end subroutine cam_init

!
!-----------------------------------------------------------------------
!
subroutine cam_run1(cam_in, cam_out)
!-----------------------------------------------------------------------
!
! Purpose:   First phase of atmosphere model run method.
!            Runs first phase of dynamics and first phase of
!            physics (before surface model updates).
!
!-----------------------------------------------------------------------
   
   use physpkg,          only: phys_run1
   use stepon,           only: stepon_run1
#if ( defined SPMD )
   use mpishorthand,     only: mpicom
#endif
   use time_manager,     only: get_nstep
   use scamMod,          only: single_column

   type(cam_in_t)  :: cam_in(begchunk:endchunk)
   type(cam_out_t) :: cam_out(begchunk:endchunk)

#if ( defined SPMD )
   real(r8) :: mpi_wtime
#endif
!-----------------------------------------------------------------------

#if ( defined SPMD )
   if (stepon_time_beg == -1.0_r8) stepon_time_beg = mpi_wtime()
   if (nstep_beg == -1) nstep_beg = get_nstep()
#endif
   if (masterproc .and. print_step_cost) then
      call t_stampf (wcstart, usrstart, sysstart)
   end if
   !----------------------------------------------------------
   ! First phase of dynamics (at least couple from dynamics to physics)
   ! Return time-step for physics from dynamics.
   !----------------------------------------------------------
   call t_barrierf ('sync_stepon_run1', mpicom)
   call t_startf ('stepon_run1')
   call stepon_run1( dtime, phys_state, phys_tend, pbuf2d, dyn_in, dyn_out )
   call t_stopf  ('stepon_run1')

   if (single_column) then
     call scam_use_iop_srf( cam_in)
   endif


   !
   !----------------------------------------------------------
   ! PHYS_RUN Call the Physics package
   !----------------------------------------------------------
   !
   call t_barrierf ('sync_phys_run1', mpicom)
   call t_startf ('phys_run1')
   call phys_run1(phys_state, dtime, phys_tend, pbuf2d,  cam_in, cam_out)
   call t_stopf  ('phys_run1')

end subroutine cam_run1

!
!-----------------------------------------------------------------------
!

subroutine cam_run2( cam_out, cam_in )
!-----------------------------------------------------------------------
!
! Purpose:   Second phase of atmosphere model run method.
!            Run the second phase physics, run methods that
!            require the surface model updates.  And run the
!            second phase of dynamics that at least couples
!            between physics to dynamics.
!
!-----------------------------------------------------------------------
   
   use physpkg,          only: phys_run2
   use stepon,           only: stepon_run2
   use time_manager,     only: is_first_step, is_first_restart_step
#if ( defined SPMD )
   use mpishorthand,     only: mpicom
#endif

   type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
   type(cam_in_t),  intent(inout) :: cam_in(begchunk:endchunk)

   !
   ! Second phase of physics (after surface model update)
   !
   call t_barrierf ('sync_phys_run2', mpicom)
   call t_startf ('phys_run2')
   call phys_run2(phys_state, dtime, phys_tend, pbuf2d,  cam_out, cam_in )
   call t_stopf  ('phys_run2')

   !
   ! Second phase of dynamics (at least couple from physics to dynamics)
   !
   call t_barrierf ('sync_stepon_run2', mpicom)
   call t_startf ('stepon_run2')
   call stepon_run2( phys_state, phys_tend, dyn_in, dyn_out )

   call t_stopf  ('stepon_run2')

   if (is_first_step() .or. is_first_restart_step()) then
      call t_startf ('cam_run2_memusage')
      call t_stopf  ('cam_run2_memusage')
   end if
end subroutine cam_run2

!
!-----------------------------------------------------------------------
!

subroutine cam_run3( cam_out )
!-----------------------------------------------------------------------
!
! Purpose:  Third phase of atmosphere model run method. This consists
!           of the third phase of the dynamics. For some dycores
!           this will be the actual dynamics run, for others the
!           dynamics happens before physics in phase 1.
!
!-----------------------------------------------------------------------
   use stepon,           only: stepon_run3
   use time_manager,     only: is_first_step, is_first_restart_step
#if ( defined SPMD )
   use mpishorthand,     only: mpicom
#endif

   type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
!-----------------------------------------------------------------------
   !
   ! Third phase of dynamics
   !
   call t_barrierf ('sync_stepon_run3', mpicom)
   call t_startf ('stepon_run3')
   call stepon_run3( dtime, cam_out, phys_state, dyn_in, dyn_out )

   call t_stopf  ('stepon_run3')

   if (is_first_step() .or. is_first_restart_step()) then
      call t_startf ('cam_run3_memusage')
      call t_stopf  ('cam_run3_memusage')
   end if
end subroutine cam_run3

!
!-----------------------------------------------------------------------
!

subroutine cam_run4( cam_out, cam_in, rstwr, nlend, &
                     yr_spec, mon_spec, day_spec, sec_spec )

!-----------------------------------------------------------------------
!
! Purpose:  Final phase of atmosphere model run method. This consists
!           of all the restart output, history writes, and other
!           file output.
!
!-----------------------------------------------------------------------
   use cam_history,      only: wshist, wrapup
   use cam_restart,      only: cam_write_restart
   use dycore,           only: dycore_is
#if ( defined SPMD )
   use mpishorthand,     only: mpicom
#endif

   type(cam_out_t), intent(inout)        :: cam_out(begchunk:endchunk)
   type(cam_in_t) , intent(inout)        :: cam_in(begchunk:endchunk)
   logical            , intent(in)           :: rstwr           ! true => write restart file
   logical            , intent(in)           :: nlend           ! true => this is final timestep
   integer            , intent(in), optional :: yr_spec         ! Simulation year
   integer            , intent(in), optional :: mon_spec        ! Simulation month
   integer            , intent(in), optional :: day_spec        ! Simulation day
   integer            , intent(in), optional :: sec_spec        ! Seconds into current simulation day

#if ( defined SPMD )
   real(r8) :: mpi_wtime
#endif
!-----------------------------------------------------------------------
! print_step_cost

   !
   !----------------------------------------------------------
   ! History and restart logic: Write and/or dispose history tapes if required
   !----------------------------------------------------------
   !
   call t_barrierf ('sync_wshist', mpicom)
   call t_startf ('wshist')
   call wshist ()
   call t_stopf  ('wshist')

#if ( defined SPMD )
   stepon_time_end = mpi_wtime()
#endif
   !
   ! Write restart files
   !
   if (rstwr) then
      call t_startf ('cam_write_restart')
      if (present(yr_spec).and.present(mon_spec).and.present(day_spec).and.present(sec_spec)) then
         call cam_write_restart( cam_in, cam_out, dyn_out, pbuf2d, &
              yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
      else
         call cam_write_restart( cam_in, cam_out, dyn_out, pbuf2d )
      end if
      call t_stopf  ('cam_write_restart')
   end if

   call t_startf ('cam_run4_wrapup')
   call wrapup(rstwr, nlend)
   call t_stopf  ('cam_run4_wrapup')

   if (masterproc .and. print_step_cost) then
      call t_startf ('cam_run4_print')
      call t_stampf (wcend, usrend, sysend)
      write(iulog,'(a,3f8.3,a)')'Prv timestep wallclock, usr, sys=', &
                            wcend-wcstart, usrend-usrstart, sysend-sysstart, &
                            ' seconds'
      call t_stopf  ('cam_run4_print')
   end if

#ifndef UNICOSMP
   call t_startf ('cam_run4_flush')
   call shr_sys_flush(iulog)
   call t_stopf  ('cam_run4_flush')
#endif

end subroutine cam_run4

!
!-----------------------------------------------------------------------
!

subroutine cam_final( cam_out, cam_in )
!-----------------------------------------------------------------------
!
! Purpose:  CAM finalization.
!
!-----------------------------------------------------------------------
   use units,            only: getunit
   use time_manager,     only: get_nstep, get_step_size
#if ( defined SPMD )
   use mpishorthand,     only: mpicom, mpiint, &
                               nsend, nrecv, nwsend, nwrecv
   use spmd_utils,       only: iam, npes
#endif
   use stepon,           only: stepon_final
   use physpkg,          only: phys_final
   use cam_initfiles,    only: cam_initfiles_close
   use camsrfexch,       only: atm2hub_deallocate, hub2atm_deallocate
   !
   ! Arguments
   !
   type(cam_out_t), pointer :: cam_out(:) ! Output from CAM to surface
   type(cam_in_t), pointer :: cam_in(:)   ! Input from merged surface to CAM
!-----------------------------------------------------------------------
   !
   ! Local variables
   !
   integer :: nstep           ! Current timestep number.

#if ( defined SPMD )   
   integer :: iu              ! SPMD Statistics output unit number
   character*24 :: filenam    ! SPMD Stats output filename
   integer :: signal          ! MPI message buffer
!------------------------------Externals--------------------------------
   real(r8) :: mpi_wtime
#endif

   call t_startf ('phys_final')
   call phys_final( phys_state, phys_tend , pbuf2d)
   call t_stopf ('phys_final')

   call t_startf ('stepon_final')
   call stepon_final(dyn_in, dyn_out)
   call t_stopf ('stepon_final')

   if(nsrest==0) then
      call t_startf ('cam_initfiles_close')
      call cam_initfiles_close()
      call t_stopf ('cam_initfiles_close')
   end if

   call hub2atm_deallocate(cam_in)
   call atm2hub_deallocate(cam_out)

#if ( defined SPMD )
   if (.false.) then
      write(iulog,*)'The following stats are exclusive of initialization/boundary datasets'
      write(iulog,*)'Number of messages sent by proc ',iam,' is ',nsend
      write(iulog,*)'Number of messages recv by proc ',iam,' is ',nrecv
   end if
#endif

   ! This flush attempts to ensure that asynchronous diagnostic prints from all 
   ! processes do not get mixed up with the "END OF MODEL RUN" message printed 
   ! by masterproc below.  The test-model script searches for this message in the 
   ! output log to figure out if CAM completed successfully.  This problem has 
   ! only been observed with the Linux Lahey compiler (lf95) which does not 
   ! support line-buffered output.  
#ifndef UNICOSMP
   call shr_sys_flush( 0 )   ! Flush all output to standard error
   call shr_sys_flush( iulog )   ! Flush all output to standard output
#endif

   if (masterproc) then
#if ( defined SPMD )
      cam_time_end = mpi_wtime()
#endif
      nstep = get_nstep()
      write(iulog,9300) nstep-1,nstep
9300  format (//'Number of completed timesteps:',i6,/,'Time step ',i6, &
                ' partially done to provide convectively adjusted and ', &
                'time filtered values for history tape.')
      write(iulog,*)'------------------------------------------------------------'
#if ( defined SPMD )
      write(iulog,*)
      write(iulog,*)' Total run time (sec) : ', cam_time_end-cam_time_beg
      write(iulog,*)' Time Step Loop run time(sec) : ', stepon_time_end-stepon_time_beg
      if (((nstep-1)-nstep_beg) > 0) then
         write(iulog,*)' SYPD : ',  &
         236.55_r8/((86400._r8/(dtime*((nstep-1)-nstep_beg)))*(stepon_time_end-stepon_time_beg))
      endif
      write(iulog,*)
#endif
      write(iulog,*)'******* END OF MODEL RUN *******'
   end if


#if ( defined SPMDSTATS )
   if (t_single_filef()) then
      write(filenam,'(a17)') 'spmdstats_cam.all'
      iu = getunit ()
      if (iam .eq. 0) then
         open (unit=iu, file=filenam, form='formatted', status='replace')
         signal = 1
      else
         call mpirecv(signal, 1, mpiint, iam-1, iam, mpicom) 
         open (unit=iu, file=filenam, form='formatted', status='old', position='append')
      endif
      write (iu,*)'************ PROCESS ',iam,' ************'
      write (iu,*)'iam ',iam,' msgs  sent =',nsend
      write (iu,*)'iam ',iam,' msgs  recvd=',nrecv
      write (iu,*)'iam ',iam,' words sent =',nwsend
      write (iu,*)'iam ',iam,' words recvd=',nwrecv
      write (iu,*)
      close(iu)
      if (iam+1 < npes) then
         call mpisend(signal, 1, mpiint, iam+1, iam+1, mpicom)
      endif
   else
      iu = getunit ()
      write(filenam,'(a14,i5.5)') 'spmdstats_cam.', iam
      open (unit=iu, file=filenam, form='formatted', status='replace')
      write (iu,*)'************ PROCESS ',iam,' ************'
      write (iu,*)'iam ',iam,' msgs  sent =',nsend
      write (iu,*)'iam ',iam,' msgs  recvd=',nrecv
      write (iu,*)'iam ',iam,' words sent =',nwsend
      write (iu,*)'iam ',iam,' words recvd=',nwrecv
      close(iu)
   endif
#endif

end subroutine cam_final

!
!-----------------------------------------------------------------------
!

end module cam_comp
