module q3d_cam
   use shr_kind_mod, only: r8 => SHR_KIND_R8, CS => SHR_KIND_CS
   use ppgrid,       only: pcols, pver

   implicit none
   private
   save

   ! Public interfaces
   public :: q3d_cam_register
   public :: q3d_cam_init
   public :: q3d_cam_tend
   public :: q3d_cam_initialize_tendency
   public :: q3d_cam_phys_to_vGCM
   public :: q3d_cam_vGCM_to_phys

#ifdef JUNG_TEST
   public :: q3d_cam_write
#endif

   ! Private data
   integer :: ixcldliq  = -1 ! cloud liquid amount index
   integer :: ixcldice  = -1 ! cloud ice amount index
   integer :: ixrain    = -1 ! rain amount index
   integer :: ixsnow    = -1 ! snow amount index
   integer :: ixgraupel = -1 ! graupel amount index

!=======================================================================
CONTAINS
!=======================================================================

   subroutine q3d_cam_register
      use constituents,   only: cnst_add
      use physconst,      only: mwh2o, cpair
      use cam_history,    only: addfld
      use q3d_runtime,    only: q3d_ncrm
      use dimensions_mod, only: ne, fv_nphys
      use parmsld,        only: parmsld_set_sizes

      ! Initialize CRM
      call cnst_add('CLDLIQ', mwh2o, cpair, 0._r8, ixcldliq,  &
           longname='Grid box averaged cloud liquid amount', is_convtran1=.true.)
      call cnst_add('CLDICE', mwh2o, cpair, 0._r8, ixcldice,  &
           longname='Grid box averaged cloud ice amount', is_convtran1=.true.)
      call cnst_add('RAINQM', mwh2o, cpair, 0._r8, ixrain,    &
           longname='Grid box averaged rain amount', is_convtran1=.true.)
      call cnst_add('SNOWQM', mwh2o, cpair, 0._r8, ixsnow,    &
           longname='Grid box averaged snow amount', is_convtran1=.true.)
      call cnst_add('GRAUPELQM', mwh2o, cpair, 0._r8, ixgraupel, &
           longname='Grid box averaged graupel amount', is_convtran1=.true.)

      ! Set VVM parameters
      call parmsld_set_sizes(ne*fv_nphys, q3d_ncrm)

   end subroutine q3d_cam_register

   subroutine q3d_cam_init()
      use time_manager,         only: get_nstep
      use cam_history,          only: addfld, add_default
      use constituents,         only: cnst_longname, cnst_name, apcnst, bpcnst
      use q3d_comm,             only: q3d_comm_init
      use q3d_interface_module, only: crm_init, crm_init_diag
      use q3d_runtime,          only: q3d_begchan, q3d_endchan, define_crm_grids
      use vvm_data_types,       only: channel => channel_data
      use vGCM_data_types,      only: vGCM => channel_vGCM

      integer                          :: step
      integer                          :: i

      ! Add diagnostic fields
      call addfld(cnst_name(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(ixcldliq))
      call addfld(apcnst(ixcldliq),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' after physics')
      call addfld(bpcnst(ixcldliq),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' before physics')
      call addfld(cnst_name(ixcldice), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(ixcldice))
      call addfld(apcnst(ixcldice),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' after physics')
      call addfld(bpcnst(ixcldice),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' before physics')
      call addfld(cnst_name(ixrain), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(ixrain))
      call addfld(apcnst(ixrain),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixrain))//' after physics')
      call addfld(bpcnst(ixrain),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixrain))//' before physics')
      call addfld(cnst_name(ixsnow), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(ixsnow))
      call addfld(apcnst(ixsnow),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixsnow))//' after physics')
      call addfld(bpcnst(ixsnow),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixsnow))//' before physics')
      call addfld(cnst_name(ixgraupel), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(ixgraupel))
      call addfld(apcnst(ixgraupel),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixgraupel))//' after physics')
      call addfld(bpcnst(ixgraupel),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixgraupel))//' before physics')

      step = get_nstep()

      ! Setup vGCM data structures and communications infrastructure
      call q3d_comm_init()
      ! Create CRM channel structures and decomposition
      call crm_init(q3d_begchan, q3d_endchan, channel, vGCM)

      ! Create coords and grids for CRM and vGCM output
      call define_crm_grids()

      ! Initialize CRM and vGCM diagnostics
      call crm_init_diag()

   end subroutine q3d_cam_init

   ! Clear vGCM tendencies for beginning of time step
   subroutine q3d_cam_initialize_tendency()
      use vGCM_data_types, only: channel_vGCM
      use q3d_runtime,     only: q3d_begchan, q3d_endchan

      integer :: chan, seg

      do chan = q3d_begchan, q3d_endchan
         do seg = 1, 4
            channel_vGCM(chan)%vGCM_tend(seg)%dT(:,:) = 0.0_r8
            channel_vGCM(chan)%vGCM_tend(seg)%dQV(:,:) = 0.0_r8
            channel_vGCM(chan)%vGCM_tend(seg)%dQC(:,:) = 0.0_r8
            channel_vGCM(chan)%vGCM_tend(seg)%dQI(:,:) = 0.0_r8
            channel_vGCM(chan)%vGCM_tend(seg)%dQR(:,:) = 0.0_r8
            channel_vGCM(chan)%vGCM_tend(seg)%dQS(:,:) = 0.0_r8
            channel_vGCM(chan)%vGCM_tend(seg)%dQG(:,:) = 0.0_r8
            channel_vGCM(chan)%vGCM_tend(seg)%dQT(:,:,:) = 0.0_r8
            channel_vGCM(chan)%vGCM_tend(seg)%dU(:,:) = 0.0_r8
            channel_vGCM(chan)%vGCM_tend(seg)%dV(:,:) = 0.0_r8
         end do
      end do

   end subroutine q3d_cam_initialize_tendency

   ! Derivatives needed for U,V,\omega,T,P,Q (in \lambda,\theta coords)
   subroutine q3d_cam_phys_to_vGCM(phys_state, phys_deriv, surf_in)
      use physics_types,   only: physics_state, physics_tend
      use ppgrid,          only: begchunk, endchunk
      use ppgrid,          only: pcols, pver, pverp
      use camsrfexch,      only: cam_in_t
      use q3d_comm,        only: q3d_comm_state_to_vGCM

      type(physics_state), intent(in)    :: phys_state(begchunk:endchunk) ! Physics state
      type(physics_tend),  intent(in)    :: phys_deriv(begchunk:endchunk) ! Physics tendency
      type(cam_in_t),      intent(in)    :: surf_in(begchunk:endchunk) ! Surface flux on CRM grid

      call q3d_comm_state_to_vGCM(phys_state)

   end subroutine q3d_cam_phys_to_vGCM

!#ifdef JUNG_TEST
!   subroutine q3d_cam_vGCM_to_phys(ztodt, phys_state, phys_deriv, surf_out, phys_dout)
!      use physics_types,   only: physics_state, physics_tend, physics_dout
!#endif

   subroutine q3d_cam_vGCM_to_phys(ztodt, phys_state, phys_deriv, surf_out)
      use physics_types,   only: physics_state, physics_tend

      use ppgrid,          only: begchunk, endchunk, pcols, pver, pverp
      use camsrfexch,      only: cam_out_t
      use q3d_comm,        only: q3d_comm_vGCM_to_tend

      real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0

      type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
      type(physics_tend),  intent(inout) :: phys_deriv(begchunk:endchunk)
      type(cam_out_t),     intent(inout) :: surf_out(begchunk:endchunk) ! Surface out on CRM grid

#ifdef JUNG_TEST
      type(physics_dout),  intent(inout) :: phys_dout(begchunk:endchunk)
#endif

      !    call physics_ptend_init(ptend, state%psetcols, 'q3d', ls=.true., lu=.true., lv=.true.)

!#ifdef JUNG_TEST
!      call q3d_comm_vGCM_to_tend(phys_state, ztodt, phys_deriv, phys_dout)
!      call q3d_comm_vGCM_to_dout(phys_state, phys_dout)
!#endif

      call q3d_comm_vGCM_to_tend(phys_state, ztodt, phys_deriv)

   end subroutine q3d_cam_vGCM_to_phys

   subroutine q3d_cam_tend(ztodt)
      !-----------------------------------------------------------------------
      !
      ! Purpose:
      !
      !-----------------------------------------------------------------------
      use shr_cal_mod,          only: shr_cal_ymd2julian
      use time_manager,         only: get_nstep, timemgr_get_calendar_cf
      use time_manager,         only: get_curr_date, get_prev_date
      use time_manager,         only: is_first_step, is_first_restart_step
      use q3d_interface_module, only: crm_run
      use vvm_data_types,       only: channel => channel_data
      use vGCM_data_types,      only: vGCM => channel_vGCM
      use q3d_runtime,          only: q3d_begchan, q3d_endchan
      !
      ! Input arguments
      !
      real(r8), intent(in) :: ztodt ! Two times model timestep (2 delta-t)
      !
      ! Output argument
      !
      !
      !---------------------------Local workspace-----------------------------

      integer              :: nstep                  ! step number
      integer              :: year, month, day, tod
      integer              :: year_c, month_c, day_c, tod_c
      character(len=CS)    :: calendar
      real(r8)             :: jday_beg, jday_end
      !-----------------------------------------------------------------------
      !


      ! Call Q3D timestep
      nstep = get_nstep()
      call get_prev_date(year, month, day, tod)
      calendar = timemgr_get_calendar_cf()
      call shr_cal_ymd2julian(year, month, day, tod, jday_beg, calendar)
      call get_curr_date(year_c, month_c, day_c, tod_c)
      call shr_cal_ymd2julian(year_c, month_c, day_c, tod_c, jday_end, calendar)

      call crm_run (is_first_restart_step(),year,jday_beg,jday_end, &
                    q3d_begchan,q3d_endchan,channel,vGCM)

   end subroutine q3d_cam_tend

#ifdef JUNG_TEST

   subroutine q3d_cam_write(phys_dout)
      !-----------------------------------------------------------------------
      !
      ! Purpose:
      !
      !-----------------------------------------------------------------------
      use ppgrid,          only: begchunk, endchunk
      use physics_types,   only: physics_dout
      !
      ! Input arguments
      type(physics_dout),  intent(in) :: phys_dout(begchunk:endchunk)
      !
      ! How to write a variable on "phys_grid""
      ! Sample here: phys_dout(:)%T, for example.

   end subroutine q3d_cam_write

#endif

end module q3d_cam
