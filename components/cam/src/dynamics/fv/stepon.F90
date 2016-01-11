!BOP
!
! !MODULE: stepon -- FV Dynamics specific time-stepping
!
! !INTERFACE:
module stepon

! !USES:
   use shr_kind_mod,       only: r8 => shr_kind_r8
   use shr_sys_mod,        only: shr_sys_flush
   use pmgrid,             only: plev, plevp, plat
   use spmd_utils,         only: iam, masterproc
   use dyn_internal_state, only: get_dyn_state, get_dyn_state_grid
   use cam_abortutils,         only: endrun
   use ppgrid,             only: begchunk, endchunk
   use physconst,          only: zvir, cappa
   use physics_types,      only: physics_state, physics_tend
   use dyn_comp,           only: dyn_import_t, dyn_export_t
   use dynamics_vars,      only: T_FVDYCORE_STATE, T_FVDYCORE_GRID
   use cam_control_mod,    only: nsrest, moist_physics
#if defined ( SPMD )
   use mpishorthand,       only: mpicom
#endif
   use perf_mod
   use cam_logfile,        only: iulog

   implicit none

   private

!
! !PUBLIC MEMBER FUNCTIONS:
!
  public stepon_init   ! Initialization
  public stepon_run1   ! run method phase 1
  public stepon_run2   ! run method phase 2
  public stepon_run3   ! run method phase 3
  public stepon_final  ! Finalization

!----------------------------------------------------------------------
!
! !DESCRIPTION: Module for FV dynamics specific time-stepping.
!
! !REVISION HISTORY:
!
! 2005.06.10  Sawyer    Adapted from FVdycore_GridCompMod
! 2005.09.16  Kluzek    Creation from stepon subroutine. 
! 2005.09.23  Kluzek    Create multiple run methods.
! 2005.11.10  Sawyer    Now using dyn_import/export_t containers
! 2006.04.13  Sawyer    Removed dependencies on prognostics
! 2006.06.29  Sawyer    Changed t3 to IJK; removed use_eta option
! 2006.07.01  Sawyer    Transitioned q3 to T_TRACERS
!
!EOP
!----------------------------------------------------------------------
!BOC
!
! !PRIVATE DATA MEMBERS:
!
  save

!-----------------------------------------------------------------------

! Magic numbers used in this module
   real(r8), parameter ::  D0_0                    =  0.0_r8
   real(r8), parameter ::  D1_0                    =  1.0_r8
   real(r8), parameter ::  D1E5                    =  1.0e5_r8

   integer :: pdt       ! Physics time step
   real(r8) :: dtime    ! Physics time step
   real(r8) :: te0            ! Total energy before dynamics

! for fv_out
   integer :: freq_diag ! Output frequency in seconds
   logical fv_monitor         ! Monitor Mean/Max/Min fields every time step
   data freq_diag  / 21600 /  ! time interval (sec) for calling fv_out
   data fv_monitor / .true. / ! This is CPU-time comsuming; set it to false for
                              ! production runs
!
! Pointers to variables in dyn_state%grid (for convenience)
   integer            :: ks
   real (r8)          :: ptop
   real (r8), pointer :: ak(:)
   real (r8), pointer :: bk(:)

   integer            :: im, jm, km, mq
   integer            :: jfirst, kfirst, jlast, klast, klastp
   integer            :: ifirstxy, ilastxy, jfirstxy, jlastxy


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_init --- Time stepping initialization
!
! !INTERFACE:
subroutine stepon_init( dyn_in, dyn_out )
! !USES:
   use constituents,       only: pcnst, cnst_get_type_byind
   use time_manager,       only: get_step_size
   use pmgrid,             only: myid_z, npr_z, twod_decomp
   use hycoef,             only: hyai, hybi, hyam, hybm
   use commap,             only: w
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! !OUTPUT PARAMETERS
!
  type (dyn_import_t)   :: dyn_in             ! Dynamics import container
  type (dyn_export_t)   :: dyn_out            ! Dynamics export container
! !DESCRIPTION:
!
! Allocate data, initialize values, setup grid locations and other
! work needed to prepare the FV dynamics to run. Return weights and 
! vertical coords to atmosphere calling program.
!
!EOP
!-----------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
   type (T_FVDYCORE_GRID), pointer :: grid
   integer i,k,j,m             ! longitude, level, latitude and tracer indices
   logical :: nlres = .false.  ! true => restart or branch run
!delta pressure dry
   real(r8), allocatable :: delpdryxy(:,:,:)
!-----------------------------------------------------------------------

   if (nsrest/=0) nlres=.true.

   grid => get_dyn_state_grid()
   im      =  grid%im
   jm      =  grid%jm
   km      =  grid%km

   jfirst    =  grid%jfirst
   jlast    =  grid%jlast
   kfirst    =  grid%kfirst
   klast    =  grid%klast
   klastp   =  grid%klastp

   ifirstxy  =  grid%ifirstxy
   ilastxy  =  grid%ilastxy
   jfirstxy  =  grid%jfirstxy
   jlastxy  =  grid%jlastxy

   ks     =  grid%ks
   ptop   =  grid%ptop
   ak     => grid%ak
   bk     => grid%bk

   !-----------------------------------------
   ! Use ak and bk as specified by CAM IC
   !-----------------------------------------
   do k = 1, km+1
      ak(k) = hyai(k) * D1E5
      bk(k) = hybi(k)
      if( bk(k) == D0_0 ) ks = k-1
   end do
   ptop = ak(1)
   if ( iam == 0 ) then
      write(iulog,*) 'Using hyai & hybi from IC:', 'KS=',ks,' PTOP=',ptop
   endif
   grid%ks   = ks
   grid%ptop = ptop

   !----------------------------------------------------------
   ! Lin-Rood dynamical core initialization
   !----------------------------------------------------------

   pdt = get_step_size()    ! Physics time step
   dtime = pdt

#if (!defined STAGGERED)
   write(iulog,*) "STEPON: pre-processor variable STAGGERED must be set"
   write(iulog,*) "Then recompile CAM. Quitting."
   call endrun
#endif

   do j = jfirstxy, jlastxy
      do i=ifirstxy, ilastxy
         dyn_in%pe(i,1,j) = ptop
      enddo
   enddo

   if ( nlres) then ! restart or branch run
      !
      ! read_restart_dynamics delivers phis, ps, u3s, v3s, delp, pt
      ! in XY decomposition

      !
      ! Do not recalculate delta pressure (delp) if this is a restart run.
      ! Re. SJ Lin: The variable "delp" (pressure thikness for a Lagrangian
      ! layer) must be in the restart file. This is because delp will be
      ! modified "after" the physics update (to account for changes in water
      ! vapor), and it can not be reproduced by surface pressure and the
      ! ETA coordinate's a's and b's.

!$omp parallel do private(i,j,k)
      do j = jfirstxy, jlastxy
        do k=1, km
          do i=ifirstxy, ilastxy
            dyn_in%pe(i,k+1,j) = dyn_in%pe(i,k,j) + dyn_in%delp(i,j,k)
          enddo
        enddo
      enddo
   else
 
      ! Initial run --> generate pe and delp from the surface pressure
 
!$omp parallel do private(i,j,k)
         do j = jfirstxy, jlastxy
            do k=1,km+1
               do i=ifirstxy, ilastxy
                  dyn_in%pe(i,k,j) = ak(k) + bk(k) * dyn_in%ps(i,j)
               enddo
            enddo
         enddo

!$omp parallel do private(i,j,k)
         do k = 1, km
            do j = jfirstxy, jlastxy
               do i= ifirstxy, ilastxy
                  dyn_in%delp(i,j,k) = dyn_in%pe(i,k+1,j) - dyn_in%pe(i,k,j)
               enddo
            enddo
         enddo
   endif

   !----------------------------------------------------------
   ! Check total dry air mass; set to 982.22 mb if initial run
   ! Print out diagnostic message if restart run
   !----------------------------------------------------------

   if ( moist_physics ) then
      call dryairm( grid, .true., dyn_in%ps, dyn_in%tracer,  &
                    dyn_in%delp, dyn_in%pe, nlres )
   endif

  ! Initialize pk, edge pressure to the cappa power.

!$omp parallel do private(i,j,k)
      do k = 1, km+1
         do j = jfirstxy, jlastxy
            do i = ifirstxy, ilastxy
               dyn_in%pk(i,j,k) = dyn_in%pe(i,k,j)**cappa
            enddo
         enddo
      enddo

   ! Generate pkz, the conversion factor betw pt and t3

   call pkez(1,      im,   km,       jfirstxy,  jlastxy,              &
             1,      km,   ifirstxy, ilastxy,    dyn_in%pe,    &
             dyn_in%pk, cappa,  ks, dyn_out%peln, dyn_out%pkz,  .false. )

   if ( .not. nlres ) then

      ! Compute pt for initial run: scaled virtual potential temperature
      ! defined as (virtual temp deg K)/pkz. pt will be written to restart (SJL)

!$omp parallel do private(i,j,k)
      do k = 1, km
         do j = jfirstxy, jlastxy
            do i = ifirstxy, ilastxy
               dyn_in%pt(i,j,k) =  dyn_in%t3(i,j,k)*            &
                (D1_0+zvir*dyn_in%tracer(i,j,k,1))    &
                /dyn_in%pkz(i,j,k) 
            enddo
         enddo
      enddo
   endif

   !----------------------------------------------------------------
   ! Convert mixing ratios initialized as dry to moist for dynamics
   !----------------------------------------------------------------

   if ( .not. nlres ) then
      ! on initial time step, dry mixing ratio advected constituents have been
      !    initialized to dry mixing ratios. dynpkg expects moist m.r. so convert here.

      ! first calculate delpdry. The set_pdel_state subroutine
      !   is called after the dynamics in d_p_coupling to set more variables.
      !   This is not in tracers.F90 because it is only used by LR dynamics.
      allocate (delpdryxy(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km))
      do k = 1, km
         do j = jfirstxy, jlastxy
            do i = ifirstxy, ilastxy
               delpdryxy(i,j,k) = dyn_in%delp(i,j,k)*          &
                    (D1_0-dyn_in%tracer(i,j,k,1))
            enddo
         enddo
      enddo
      do m = 1,pcnst
         if (cnst_get_type_byind(m).eq.'dry') then
            do k=1, km
               do j = jfirstxy, jlastxy
                  do i = ifirstxy, ilastxy
                     dyn_in%tracer(i,j,k,m) =               &
                             dyn_in%tracer(i,j,k,m)*        &
                             delpdryxy(i,j,k)/dyn_in%delp(i,j,k)
                  end do
               end do
            end do
         end if
      end do
      deallocate (delpdryxy)
      
   end if  ! .not. nlres

!EOC
end subroutine stepon_init

!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_run1 -- Phase 1 of dynamics run method.
!
! !INTERFACE:
subroutine stepon_run1( dtime_out, phys_state, phys_tend, pbuf2d,        &
                        dyn_in, dyn_out )
!-----------------------------------------------------------------------
!
!  ATTENTION *** ATTENTION *** ATTENTION *** ATTENTION *** ATTENTION
!
!
!     A 2D xy decomposition is used for handling the Lagrangian surface
!     remapping, the ideal physics, and (optionally) the geopotential
!     calculation.
!
!     The transpose from yz to xy decomposition takes place within dynpkg.
!     The xy decomposed variables are then transposed directly to the
!     physics decomposition within d_p_coupling.
!
!     The xy decomposed variables have names corresponding to the
!     yz decomposed variables: simply append "xy". Thus, "uxy" is the
!     xy decomposed version of "u".
!
!     To assure that the latitudinal decomposition operates
!     as efficiently as before, a separate parameter "twod_decomp" has
!     been defined; a value of 1 refers to the multi-2D decomposition with
!     transposes; a value of 0 means that the decomposition is effectively
!     one-dimensional, thereby enabling the transpose logic to be skipped;
!     there is an option to force computation of transposes even for case
!     where decomposition is effectively 1-D.
!
!     For questions/comments, contact Art Mirin, mirin@llnl.gov
!
!-----------------------------------------------------------------------
! !USES:
   use dp_coupling,       only: d_p_coupling
   use dyn_comp,          only: dyn_run
   
   use physics_buffer,    only: physics_buffer_desc
   use pmgrid,            only: twod_decomp
   use advect_tend,       only: compute_adv_tends_xyz
   use fv_control_mod,    only: nsplit, nspltrac
!-----------------------------------------------------------------------
! !OUTPUT PARAMETERS:
!
   real(r8), intent(out) :: dtime_out   ! Time-step
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   type(dyn_import_t)                 :: dyn_in  ! Dynamics import container
   type(dyn_export_t)                 :: dyn_out ! Dynamics export container

   type(T_FVDYCORE_STATE), pointer :: dyn_state

! !DESCRIPTION:
!
!  Phase 1 run of FV dynamics. Run the dynamics, and couple to physics.
!
!EOP
!-----------------------------------------------------------------------
!BOC
   
   integer  :: rc 
#if (! defined SPMD)
   integer  :: mpicom = 0
#endif

   dtime_out = dtime
   dyn_state => get_dyn_state()

   ! Dump state variables to IC file
   call t_barrierf('sync_diag_dynvar_ic', mpicom)
   call t_startf ('diag_dynvar_ic')
   call diag_dynvar_ic (dyn_state%grid, dyn_out%phis, dyn_out%ps,             &
                        dyn_out%t3, dyn_out%u3s, dyn_out%v3s, dyn_out%tracer  )
   call t_stopf  ('diag_dynvar_ic')

   call t_startf ('comp_adv_tends1')
   call compute_adv_tends_xyz(dyn_state%grid, dyn_in%tracer )
   call t_stopf  ('comp_adv_tends1')
   !
   !--------------------------------------------------------------------------
   ! Perform finite-volume dynamics -- this dynamical core contains some 
   ! yet to be published algorithms. Its use in the CAM is
   ! for software development purposes only. 
   ! Please contact S.-J. Lin (Shian-Jiann.Lin@noaa.gov)
   ! if you plan to use this mudule for scientific purposes. Contact S.-J. Lin
   ! or Will Sawyer (sawyer@gmao.gsfc.nasa.gov) if you plan to modify the
   ! software.
   !--------------------------------------------------------------------------

   !----------------------------------------------------------
   ! For 2-D decomposition, phisxy is input to dynpkg, and the other
   ! xy variables are output. Some are computed through direct
   ! transposes, and others are derived.
   !----------------------------------------------------------
   call t_barrierf('sync_dyn_run', mpicom)
   call t_startf ('dyn_run')
   call dyn_run(ptop,      pdt,     te0,         &
                dyn_state, dyn_in,  dyn_out,  rc )
   if ( rc /= 0 ) then
     write(iulog,*) "STEPON_RUN: dyn_run returned bad error code", rc
     write(iulog,*) "Quitting."
     call endrun
   endif 
   call t_stopf  ('dyn_run')

   call t_startf ('comp_adv_tends2')
   call compute_adv_tends_xyz(dyn_state%grid, dyn_out%tracer )
   call t_stopf  ('comp_adv_tends2')

   !----------------------------------------------------------
   ! Move data into phys_state structure.
   !----------------------------------------------------------
   call t_barrierf('sync_d_p_coupling', mpicom)
   call t_startf('d_p_coupling')
   call d_p_coupling(dyn_state%grid, phys_state, phys_tend,  pbuf2d, dyn_out)
   call t_stopf('d_p_coupling')

!EOC
end subroutine stepon_run1

!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_run2 -- second phase run method
!
! !INTERFACE:
subroutine stepon_run2( phys_state, phys_tend, dyn_in, dyn_out )
! !USES:
   use dp_coupling,      only: p_d_coupling
!
! !INPUT/OUTPUT PARAMETERS:
!
   type(physics_state), intent(inout):: phys_state(begchunk:endchunk)
   type(physics_tend), intent(inout):: phys_tend(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
   type (T_FVDYCORE_GRID), pointer :: grid
!
! !DESCRIPTION:
!
! Second phase run method. Couple from physics to dynamics.
!
!EOP
!-----------------------------------------------------------------------
!BOC
#if (! defined SPMD)
   integer  :: mpicom = 0
#endif

!-----------------------------------------------------------------------

   !----------------------------------------------------------
   ! Update dynamics variables using phys_state & phys_tend.
   ! 2-D decomposition: Compute ptxy and q3xy; for ideal
   !   physics, scale ptxy by (old) pkzxy; then transpose to yz variables
   ! 1-D decomposition: Compute dudt, dvdt, pt and q3; for ideal physics,
   !   scale pt by old pkz.
   ! Call uv3s_update to update u3s and v3s from dudt and dvdt.
   ! Call p_d_adjust to update pt, q3, pe, delp, ps, piln, pkz and pk.
   ! For adiabatic case, transpose to yz variables.
   !----------------------------------------------------------
   grid => get_dyn_state_grid()

   call t_barrierf('sync_p_d_coupling', mpicom)
   call t_startf ('p_d_coupling')
   call p_d_coupling(grid, phys_state, phys_tend, &
                     dyn_in, dtime, zvir, cappa, ptop)
   call t_stopf  ('p_d_coupling')

!EOC
end subroutine stepon_run2

!-----------------------------------------------------------------------

subroutine stepon_run3( dtime, cam_out, phys_state, dyn_in, dyn_out )
! !USES:
   use time_manager,     only: get_curr_date
   use fv_prints,        only: fv_out
   use camsrfexch,       only: cam_out_t    
   use pmgrid,           only: twod_decomp
!
! !INPUT PARAMETERS:
!
   type(physics_state), intent(in):: phys_state(begchunk:endchunk)
   real(r8), intent(in) :: dtime            ! Time-step
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
!
! !INPUT/OUTPUT PARAMETERS:
!
   type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
!
! !DESCRIPTION:
!
!	Final run phase of dynamics. Some printout and time index updates.
!
! !HISTORY:
!   2005.09.16  Kluzek     Creation
!   2006.04.13  Sawyer     Removed shift_time_indices (not needed in FV)
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   type (T_FVDYCORE_GRID), pointer :: grid
   integer :: ncdate            ! current date in integer format [yyyymmdd]
   integer :: ncsec             ! time of day relative to current date [seconds]
   integer :: yr, mon, day      ! year, month, day components of a date
   integer :: ncsecp
#if (! defined SPMD)
   integer  :: mpicom = 0
#endif

   !----------------------------------------------------------
   ! Monitor max/min/mean of selected fields
   !
   !  SEE BELOW  ****  SEE BELOW  ****  SEE BELOW
   
   ! Beware that fv_out uses both dynamics and physics instantiations.
   ! However, I think that they are used independently, so that the
   ! answers are correct. Still, this violates the notion that the
   ! physics state is no longer active after p_d_coupling.
   !----------------------------------------------------------
   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day
   ncsecp = pdt + ncsec      !  step complete, but nstep not incremented yet

   if ( fv_monitor .and. mod(ncsecp, freq_diag) == 0 ) then
      grid => get_dyn_state_grid()

      call t_barrierf('sync_fv_out', mpicom)
      call t_startf('fv_out')
      call fv_out(grid, dyn_out%pk, dyn_out%pt,         &
                  ptop, dyn_out%ps, dyn_out%tracer,     &
                  dyn_out%delp, dyn_out%pe, cam_out,    &
                   phys_state, ncdate, ncsecp, moist_physics)
      call t_stopf('fv_out')
   endif

!EOC
end subroutine stepon_run3

!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_final --- Dynamics finalization
!
! !INTERFACE:
subroutine stepon_final(dyn_in, dyn_out)

! !PARAMETERS:
  type (dyn_import_t), intent(out) :: dyn_in  ! Dynamics import container
  type (dyn_export_t), intent(out) :: dyn_out ! Dynamics export container
!
! !DESCRIPTION:
!
! Deallocate data needed for dynamics. Finalize any dynamics specific
! files or subroutines.
!
!EOP
!-----------------------------------------------------------------------
!BOC

!!! Not yet ready for the call to dyn_final
!!! call dyn_final( RESTART_FILE, dyn_state, dyn_in, dyn_out )
!EOC
end subroutine stepon_final

!-----------------------------------------------------------------------

end module stepon
