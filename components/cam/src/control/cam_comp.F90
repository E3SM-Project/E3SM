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
   use physpkg,          only: phys_init, phys_register
   
   use dycore,           only: dycore_is
#if (defined BFB_CAM_SCAM_IOP)
   use history_defaults, only: initialize_iop_history
#endif
!   use shr_orb_mod,      only: shr_orb_params
   use camsrfexch,       only: hub2atm_alloc, atm2hub_alloc
   use cam_history,      only: intht, init_masterlinkedlist
   use history_scam,     only: scm_intht
   use scamMod,          only: single_column
   use cam_pio_utils,    only: init_pio_subsystem
   use cam_instance,     only: inst_suffix

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
   ! Initialize index values for advected and non-advected tracers
   !
   call phys_register ()
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

      call cam_initfiles_open()
      call cam_initial(dyn_in, dyn_out, NLFileName=filein)

      ! Allocate and setup surface exchange data
      call atm2hub_alloc(cam_out)
      call hub2atm_alloc(cam_in)

   else

      call cam_read_restart ( cam_in, cam_out, dyn_in, dyn_out, pbuf2d, stop_ymd, stop_tod, NLFileName=filein )

     ! Commented out the hub2atm_alloc call as it overwrite cam_in, which is undesirable. The fields in cam_in are necessary for getting BFB restarts
	 ! There are no side effects of commenting out this call as this call allocates cam_in and cam_in allocation has already been done in cam_init
     !call hub2atm_alloc( cam_in )
#if (defined BFB_CAM_SCAM_IOP)
      call initialize_iop_history()
#endif
   end if


   call phys_init( phys_state, phys_tend, pbuf2d,  cam_out )

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

   call phys_final( phys_state, phys_tend , pbuf2d)
   call stepon_final(dyn_in, dyn_out)

   if(nsrest==0) then
      call cam_initfiles_close()
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
