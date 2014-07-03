!===============================================================================
!BOP ===========================================================================
!
! !MODULE: seq_hist_mod -- cpl7 history writing routines
!
! !DESCRIPTION:
!
!    Creates cpl7 history files, instantanious, time-avg, and auxilliary
!
! !REVISION HISTORY:
!     2009-Sep-25 - B. Kauffman - move from cpl7 main program into hist module
!     2009-mmm-dd - T. Craig - initial versions
!
! !INTERFACE: ------------------------------------------------------------------

module seq_hist_mod

! !USES:

   use shr_kind_mod,      only: R8 => SHR_KIND_R8, IN => SHR_KIND_IN
   use shr_kind_mod,      only: CL => SHR_KIND_CL, CS => SHR_KIND_CS
   use shr_sys_mod,       only: shr_sys_abort, shr_sys_flush
   use shr_cal_mod,       only: shr_cal_date2ymd
   use mct_mod           ! adds mct_ prefix to mct lib
   use ESMF

   use seq_infodata_mod  ! "infodata" gathers various control flags into one datatype
   use seq_timemgr_mod   ! clock & alarm routines
   use seq_io_mod        ! lower level io routines

   use seq_comm_mct     , only: seq_comm_getdata=>seq_comm_setptrs
   use seq_comm_mct,      only: seq_comm_setnthreads, seq_comm_iamin
   use seq_comm_mct,      only: CPLID, GLOID, logunit, loglevel
   use seq_comm_mct,      only: num_inst_atm, num_inst_lnd, num_inst_ocn
   use seq_comm_mct,      only: num_inst_ice, num_inst_glc, num_inst_wav
   use seq_comm_mct,      only: num_inst_rof

   use prep_ocn_mod,      only: prep_ocn_get_r2x_ox
   use prep_ocn_mod,      only: prep_ocn_get_x2oacc_ox
   use prep_ocn_mod,      only: prep_ocn_get_x2oacc_ox_cnt  
   use prep_atm_mod,      only: prep_atm_get_o2x_ax
   use prep_aoflux_mod,   only: prep_aoflux_get_xao_ox
   use prep_aoflux_mod,   only: prep_aoflux_get_xao_ax

   use component_type_mod

   implicit none

   private

! !PUBLIC TYPES:
  
   ! no public types

! !PUBLIC MEMBER FUNCTIONS

   public :: seq_hist_write     ! write instantaneous hist file
   public :: seq_hist_writeavg  ! write time-avg      hist file
   public :: seq_hist_writeaux  ! write auxiliary     hist files
   public :: seq_hist_spewav    ! write avs to history file for debugging

! !PUBLIC DATA MEMBERS:

   ! no public data

!EOP

   !----------------------------------------------------------------------------
   ! local/module data
   !----------------------------------------------------------------------------

   logical     :: iamin_CPLID            ! pe associated with CPLID
   integer(IN) :: mpicom_GLOID           ! MPI global communicator
   integer(IN) :: mpicom_CPLID           ! MPI cpl communicator

   integer(IN) :: nthreads_GLOID         ! OMP global number of threads
   integer(IN) :: nthreads_CPLID         ! OMP cpl number of threads
   logical     :: drv_threading          ! driver threading control

   logical     :: atm_present            ! .true.  => atm is present
   logical     :: lnd_present            ! .true.  => land is present
   logical     :: ice_present            ! .true.  => ice is present
   logical     :: ocn_present            ! .true.  => ocn is present
   logical     :: rof_present            ! .true.  => land runoff is present
   logical     :: glc_present            ! .true.  => glc is present
   logical     :: wav_present            ! .true.  => wav is present
   
   logical     :: atm_prognostic         ! .true.  => atm comp expects input
   logical     :: lnd_prognostic         ! .true.  => lnd comp expects input
   logical     :: ice_prognostic         ! .true.  => ice comp expects input
   logical     :: ocn_prognostic         ! .true.  => ocn comp expects input
   logical     :: ocnrof_prognostic      ! .true.  => ocn comp expects runoff input
   logical     :: rof_prognostic         ! .true.  => rof comp expects input
   logical     :: glc_prognostic         ! .true.  => glc comp expects input
   logical     :: wav_prognostic         ! .true.  => wav comp expects input

   logical     :: cdf64                  ! true => use 64 bit addressing in netCDF files

   !--- domain equivalent 2d grid size ---
   integer(IN) :: atm_nx, atm_ny         ! nx,ny of 2d grid, if known
   integer(IN) :: lnd_nx, lnd_ny         ! nx,ny of 2d grid, if known
   integer(IN) :: ice_nx, ice_ny         ! nx,ny of 2d grid, if known
   integer(IN) :: ocn_nx, ocn_ny         ! nx,ny of 2d grid, if known
   integer(IN) :: rof_nx, rof_ny         ! nx,ny of 2d grid, if known
   integer(IN) :: glc_nx, glc_ny         ! nx,ny of 2d grid, if known
   integer(IN) :: wav_nx, wav_ny         ! nx,ny of 2d grid, if known

   integer(IN) :: info_debug = 0         ! local info_debug level

   !--- temporary pointers ---
   type(mct_aVect), pointer :: r2x_ox(:)
   type(mct_aVect), pointer :: x2oacc_ox(:)
   integer        , pointer :: x2oacc_ox_cnt
   type(mct_aVect), pointer :: xao_ox(:)
   type(mct_aVect), pointer :: xao_ax(:)
   type(mct_aVect), pointer :: o2x_ax(:)

!===============================================================================
contains
!===============================================================================

subroutine seq_hist_write(infodata, EClock_d, &
     atm, lnd, ice, ocn, rof, glc, wav, &
     fractions_ax, fractions_lx, fractions_ix, fractions_ox, fractions_rx,  &
     fractions_gx, fractions_wx)

   implicit none
   !
   ! Arguments
   type(seq_infodata_type) , intent(in)    :: infodata
   type (ESMF_Clock)       , intent(in)    :: EClock_d   ! driver clock
   type (component_type)   , intent(inout) :: atm(:)
   type (component_type)   , intent(inout) :: lnd(:)
   type (component_type)   , intent(inout) :: ice(:)
   type (component_type)   , intent(inout) :: ocn(:)
   type (component_type)   , intent(inout) :: rof(:)
   type (component_type)   , intent(inout) :: glc(:)
   type (component_type)   , intent(inout) :: wav(:)
   type(mct_aVect)         , intent(inout) :: fractions_ax(:) ! Fractions on atm grid/decomp
   type(mct_aVect)         , intent(inout) :: fractions_lx(:) ! Fractions on lnd grid/decomp
   type(mct_aVect)         , intent(inout) :: fractions_ix(:) ! Fractions on ice grid/decomp
   type(mct_aVect)         , intent(inout) :: fractions_ox(:) ! Fractions on ocn grid/decomp
   type(mct_aVect)         , intent(inout) :: fractions_rx(:) ! Fractions on rof grid/decomp
   type(mct_aVect)         , intent(inout) :: fractions_gx(:) ! Fractions on glc grid/decomp
   type(mct_aVect)         , intent(inout) :: fractions_wx(:) ! Fractions on wav grid/decomp
   !
   ! Local Variables
   integer(IN)   :: curr_ymd     ! Current date YYYYMMDD
   integer(IN)   :: curr_tod     ! Current time-of-day (s)
   integer(IN)   :: start_ymd    ! Starting date YYYYMMDD
   integer(IN)   :: start_tod    ! Starting time-of-day (s)
   real(r8)      :: curr_time    ! Time interval since reference time
   integer(IN)   :: yy,mm,dd     ! year, month, day
   integer(IN)   :: fk           ! index
   character(CL) :: time_units   ! units of time variable
   character(CL) :: calendar     ! calendar type
   character(CL) :: case_name    ! case name
   character(CL) :: hist_file    ! Local path to history filename
   integer(IN)   :: lsize        ! local size of an aVect
   real(r8)      :: tbnds(2)     ! CF1.0 time bounds
   logical       :: whead,wdata  ! for writing restart/history cdf files

   type(mct_gsMap), pointer :: gsmap
   type(mct_gGrid), pointer :: dom    ! comp domain on cpl pes
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! get required infodata
   !----------------------------------------------------------------------------
   iamin_CPLID  = seq_comm_iamin(CPLID)

   call seq_comm_getdata(GLOID,mpicom=mpicom_GLOID,nthreads=nthreads_GLOID)
   call seq_comm_getdata(CPLID,mpicom=mpicom_CPLID,nthreads=nthreads_CPLID)

   call seq_infodata_getData(infodata,       &
        drv_threading=drv_threading,         &
        atm_present=atm_present,             &
        lnd_present=lnd_present,             &
        rof_present=rof_present,             &
        ice_present=ice_present,             &
        ocn_present=ocn_present,             &
        glc_present=glc_present,             &
        wav_present=wav_present,             &
        atm_prognostic=atm_prognostic,       &
        lnd_prognostic=lnd_prognostic,       &
        ice_prognostic=ice_prognostic,       &
        ocn_prognostic=ocn_prognostic,       &
        ocnrof_prognostic=ocnrof_prognostic, &
        rof_prognostic=rof_prognostic,       &
        glc_prognostic=glc_prognostic,       &
        wav_prognostic=wav_prognostic,       &
        atm_nx=atm_nx, atm_ny=atm_ny,        &
        lnd_nx=lnd_nx, lnd_ny=lnd_ny,        &
        rof_nx=rof_nx, rof_ny=rof_ny,        &
        ice_nx=ice_nx, ice_ny=ice_ny,        &
        glc_nx=glc_nx, glc_ny=glc_ny,        &
        wav_nx=wav_nx, wav_ny=wav_ny,        &
        ocn_nx=ocn_nx, ocn_ny=ocn_ny,        &
        cpl_cdf64=cdf64,                     & 
        case_name=case_name)

   !--- Get current date from clock needed to label the history pointer file ---

   call seq_timemgr_EClockGetData( EClock_d, curr_ymd=curr_ymd, curr_tod=curr_tod, &
        start_ymd=start_ymd, start_tod=start_tod, curr_time=curr_time, &
        calendar=calendar)
   call shr_cal_date2ymd(curr_ymd,yy,mm,dd)
   write(hist_file,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
      trim(case_name), '.cpl.hi.', yy,'-',mm,'-',dd,'-',curr_tod,'.nc'

   time_units = 'days since ' &
        // seq_io_date2yyyymmdd(start_ymd) // ' ' // seq_io_sec2hms(start_tod)

   if (iamin_CPLID) then

      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
      call seq_io_wopen(hist_file,clobber=.true.,cdf64=cdf64)

      ! loop twice, first time write header, second time write data for perf

      do fk = 1,2
         if (fk == 1) then
            whead = .true.
            wdata = .false.
         elseif (fk == 2) then
            whead = .false.
            wdata = .true.
            call seq_io_enddef(hist_file)
         else
            call shr_sys_abort('seq_hist_write fk illegal')
         end if

         tbnds = curr_time
         !------- tcx nov 2011 tbnds of same values causes problems in ferret
         if (tbnds(1) >= tbnds(2)) then
            call seq_io_write(hist_file,&
                 time_units=time_units, time_cal=calendar, time_val=curr_time, &
                 whead=whead, wdata=wdata)
         else
            call seq_io_write(hist_file, &
                 time_units=time_units, time_cal=calendar, time_val=curr_time, &
                 whead=whead, wdata=wdata, tbnds=tbnds)
         endif

         if (atm_present) then
            gsmap => component_get_gsmap_cx(atm(1))
            dom   => component_get_dom_cx(atm(1))
            call seq_io_write(hist_file, gsmap, dom%data, 'dom_ax',  &
                 nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata, pre='doma')
            call seq_io_write(hist_file, gsmap, fractions_ax, 'fractions_ax',  &
                 nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata, pre='fraca')
            call seq_io_write(hist_file, atm, 'x2c', 'x2a_ax', &
                 nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata, pre='x2a')
            call seq_io_write(hist_file, atm, 'c2x', 'a2x_ax', &
                 nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata, pre='a2x')
            !call seq_io_write(hist_file, gsmap, l2x_ax, 'l2x_ax',  &
            !    nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata, pre='l2x_ax')
            !call seq_io_write(hist_file, gsmap, o2x_ax, 'o2x_ax',  &
            !    nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata, pre='o2x_ax')
            !call seq_io_write(hist_file, gsmap, i2x_ax, 'i2x_ax',  &
            !    nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata, pre='i2x_ax')
         endif

         if (lnd_present) then
            gsmap => component_get_gsmap_cx(lnd(1))
            dom   => component_get_dom_cx(lnd(1))
            call seq_io_write(hist_file, gsmap, dom%data, 'dom_lx',  &
                 nx=lnd_nx, ny=lnd_ny, nt=1, whead=whead, wdata=wdata, pre='doml')
            call seq_io_write(hist_file, gsmap, fractions_lx, 'fractions_lx',  &
                 nx=lnd_nx, ny=lnd_ny, nt=1, whead=whead, wdata=wdata, pre='fracl')
            call seq_io_write(hist_file, lnd, 'c2x', 'l2x_lx', &
                 nx=lnd_nx, ny=lnd_ny, nt=1, whead=whead, wdata=wdata, pre='l2x')
            call seq_io_write(hist_file, lnd, 'x2c', 'x2l_lx',&
                 nx=lnd_nx, ny=lnd_ny, nt=1, whead=whead, wdata=wdata, pre='x2l')
         endif

         if (rof_present) then
            gsmap => component_get_gsmap_cx(rof(1))
            dom   => component_get_dom_cx(rof(1))
            call seq_io_write(hist_file, gsmap, dom%data, 'dom_rx',  &
                 nx=rof_nx, ny=rof_ny, nt=1, whead=whead, wdata=wdata, pre='domr')
            call seq_io_write(hist_file, gsmap, fractions_rx, 'fractions_rx',  &
                 nx=rof_nx, ny=rof_ny, nt=1, whead=whead, wdata=wdata, pre='fracr')
            call seq_io_write(hist_file, rof, 'c2x', 'r2x_rx', &
                 nx=rof_nx, ny=rof_ny, nt=1, whead=whead, wdata=wdata, pre='r2x')
            call seq_io_write(hist_file, rof, 'x2c', 'x2r_rx', &
                 nx=rof_nx, ny=rof_ny, nt=1, whead=whead, wdata=wdata, pre='x2r')
         endif

         if (rof_present .and. ocnrof_prognostic) then
            gsmap  => component_get_gsmap_cx(ocn(1))
            r2x_ox => prep_ocn_get_r2x_ox()
            call seq_io_write(hist_file, gsmap, r2x_ox, 'r2x_ox',  &
                 nx=ocn_nx, ny=ocn_ny, nt=1, whead=whead, wdata=wdata, pre='r2xo')
         endif

         if (ocn_present) then
            gsmap => component_get_gsmap_cx(ocn(1))
            dom   => component_get_dom_cx(ocn(1))
            call seq_io_write(hist_file, gsmap, dom%data, 'dom_ox',  &
                              nx=ocn_nx, ny=ocn_ny, nt=1, whead=whead, wdata=wdata, pre='domo')
            call seq_io_write(hist_file, gsmap, fractions_ox, 'fractions_ox',  &
                              nx=ocn_nx, ny=ocn_ny, nt=1, whead=whead, wdata=wdata, pre='fraco')
            call seq_io_write(hist_file, ocn, 'c2x', 'o2x_ox', &
                              nx=ocn_nx, ny=ocn_ny, nt=1, whead=whead, wdata=wdata, pre='o2x')
            !call seq_io_write(hist_file, ocn, 'x2c', 'x2o_ox', &
            !                  nx=ocn_nx, ny=ocn_ny, nt=1, whead=whead, wdata=wdata, pre='x2o')

            gsmap     => component_get_gsmap_cx(ocn(1))
            x2oacc_ox => prep_ocn_get_x2oacc_ox()
            call seq_io_write(hist_file, gsmap, x2oacc_ox, 'x2oacc_ox',  &
                 nx=ocn_nx, ny=ocn_ny, nt=1, whead=whead, wdata=wdata, pre='x2oacc')

            gsmap         => component_get_gsmap_cx(ocn(1))
            x2oacc_ox_cnt => prep_ocn_get_x2oacc_ox_cnt()
            call seq_io_write(hist_file, x2oacc_ox_cnt, 'x2oacc_ox_cnt',  &
                 whead=whead, wdata=wdata)
            gsmap => component_get_gsmap_cx(ocn(1))
            xao_ox => prep_aoflux_get_xao_ox()
            call seq_io_write(hist_file, gsmap, xao_ox, 'xao_ox',  &
                 nx=ocn_nx, ny=ocn_ny, nt=1, whead=whead, wdata=wdata, pre='xaoo')

            gsmap  => component_get_gsmap_cx(atm(1))
            o2x_ax => prep_atm_get_o2x_ax()
            call seq_io_write(hist_file, gsmap, o2x_ax, 'o2x_ax',  &
                 nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata, pre='o2xa')

            gsmap  => component_get_gsmap_cx(atm(1))
            xao_ax => prep_aoflux_get_xao_ax()
            call seq_io_write(hist_file, gsmap, xao_ax, 'xao_ax',  &
                 nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata, pre='xaoa')
         endif

         if (ice_present) then
            gsmap => component_get_gsmap_cx(ice(1))
            dom   => component_get_dom_cx(ice(1))
            call seq_io_write(hist_file, gsmap, dom%data, 'dom_ix',  &
                 nx=ice_nx, ny=ice_ny, nt=1, whead=whead, wdata=wdata, pre='domi')
            call seq_io_write(hist_file, gsmap, fractions_ix, 'fractions_ix',  &
                 nx=ice_nx, ny=ice_ny, nt=1, whead=whead, wdata=wdata, pre='fraci')
            call seq_io_write(hist_file, ice, 'c2x', 'i2x_ix', &
                 nx=ice_nx, ny=ice_ny, nt=1, whead=whead, wdata=wdata, pre='i2x')
            call seq_io_write(hist_file, ice, 'x2c', 'x2i_ix', &
                 nx=ice_nx, ny=ice_ny, nt=1, whead=whead, wdata=wdata, pre='x2i')
         endif

         if (glc_present) then
            gsmap => component_get_gsmap_cx(glc(1))
            dom   => component_get_dom_cx(glc(1))
            call seq_io_write(hist_file, gsmap, dom%data, 'dom_gx',  &
                 nx=glc_nx, ny=glc_ny, nt=1, whead=whead, wdata=wdata, pre='domg')
            call seq_io_write(hist_file, gsmap, fractions_gx, 'fractions_gx',  &
                 nx=glc_nx, ny=glc_ny, nt=1, whead=whead, wdata=wdata, pre='fracg')
            call seq_io_write(hist_file, glc, 'c2x', 'g2x_gx', &
                 nx=glc_nx, ny=glc_ny, nt=1, whead=whead, wdata=wdata, pre='g2x')
            call seq_io_write(hist_file, glc, 'x2c', 'x2g_gx', &
                 nx=glc_nx, ny=glc_ny, nt=1, whead=whead, wdata=wdata, pre='x2g')
         endif

         if (wav_present) then
            gsmap => component_get_gsmap_cx(wav(1))
            dom   => component_get_dom_cx(wav(1))
            call seq_io_write(hist_file, gsmap, dom%data, 'dom_wx',  &
                 nx=wav_nx, ny=wav_ny, nt=1, whead=whead, wdata=wdata, pre='domw')
            call seq_io_write(hist_file, gsmap, fractions_wx, 'fractions_wx',  &
                 nx=wav_nx, ny=wav_ny, nt=1, whead=whead, wdata=wdata, pre='fracw')
            call seq_io_write(hist_file, wav, 'c2x', 'w2x_wx', &
                 nx=wav_nx, ny=wav_ny, nt=1, whead=whead, wdata=wdata, pre='w2x')
            call seq_io_write(hist_file, wav, 'x2c', 'x2w_wx', &
                 nx=wav_nx, ny=wav_ny, nt=1, whead=whead, wdata=wdata, pre='x2w')
         endif
      enddo

      call seq_io_close(hist_file)
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
   endif

end subroutine seq_hist_write

!===============================================================================

subroutine seq_hist_writeavg(infodata, EClock_d, &
     atm, lnd, ice, ocn, rof, glc, wav, write_now)

   implicit none

   type(seq_infodata_type) ,  intent(in) :: infodata
   type (ESMF_Clock)       ,  intent(in) :: EClock_d   ! driver clock
   type (component_type)   ,  intent(in) :: atm(:)
   type (component_type)   ,  intent(in) :: lnd(:)
   type (component_type)   ,  intent(in) :: ice(:)
   type (component_type)   ,  intent(in) :: ocn(:)
   type (component_type)   ,  intent(in) :: rof(:)
   type (component_type)   ,  intent(in) :: glc(:)
   type (component_type)   ,  intent(in) :: wav(:)
   logical                 ,  intent(in) :: write_now  ! write or accumulate

   integer(IN)           :: curr_ymd     ! Current date YYYYMMDD
   integer(IN)           :: curr_tod     ! Current time-of-day (s)
   integer(IN)           :: prev_ymd     ! Previous date YYYYMMDD
   integer(IN)           :: prev_tod     ! Previous time-of-day (s)
   integer(IN)           :: start_ymd    ! Starting date YYYYMMDD
   integer(IN)           :: start_tod    ! Starting time-of-day (s)
   real(r8)              :: curr_time    ! Time interval since reference time
   real(r8)              :: prev_time    ! Time interval since reference time
   real(r8)              :: avg_time     ! Average time of tavg
   integer(IN)           :: yy, mm, dd     ! year,  month,  day
   integer(IN)           :: fk           ! index
   character(CL)         :: time_units   ! units of time variable
   character(CL)         :: calendar     ! calendar type
   integer(IN)           :: lsize        ! local size of an aVect
   character(CL)         :: case_name    ! case name
   character(CL)         :: hist_file    ! Local path to history filename
   logical               :: whead, wdata  ! flags write header vs. data
   integer(IN)           :: iidx ! component instance counter

   type(mct_aVect), save  :: a2x_ax_avg(num_inst_atm)   ! tavg aVect/bundle
   type(mct_aVect), save  :: x2a_ax_avg(num_inst_atm)
   type(mct_aVect), save  :: l2x_lx_avg(num_inst_lnd)
   type(mct_aVect), save  :: x2l_lx_avg(num_inst_lnd)
   type(mct_aVect), save  :: r2x_rx_avg(num_inst_rof)
   type(mct_aVect), save  :: x2r_rx_avg(num_inst_rof)
   type(mct_aVect), save  :: o2x_ox_avg(num_inst_ocn)
   type(mct_aVect), save  :: x2o_ox_avg(num_inst_ocn)
   type(mct_aVect), save  :: i2x_ix_avg(num_inst_ice)
   type(mct_aVect), save  :: x2i_ix_avg(num_inst_ice)
   type(mct_aVect), save  :: g2x_gx_avg(num_inst_glc)
   type(mct_aVect), save  :: x2g_gx_avg(num_inst_glc)
   type(mct_aVect), save  :: w2x_wx_avg(num_inst_wav)
   type(mct_aVect), save  :: x2w_wx_avg(num_inst_wav)

   integer(IN)    , save  :: cnt                 ! counts samples in tavg
   real(r8)       , save  :: tbnds(2)            ! CF1.0 time bounds

   logical        , save  :: first_call = .true. ! flags 1st call of this routine

   type(mct_gsMap),  pointer :: gsmap  ! component decomp on cpl pes
   type(mct_gGrid),  pointer :: dom    ! component domain on cpl pes
   type(mct_avect),  pointer :: c2x    ! component->coupler avs on cpl pes
   type(mct_avect),  pointer :: x2c    ! coupler->component avs on cpl pes
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! get required infodata
   !----------------------------------------------------------------------------
   iamin_CPLID  = seq_comm_iamin(CPLID)
   call seq_comm_getdata(GLOID, &
        mpicom=mpicom_GLOID, nthreads=nthreads_GLOID)
   call seq_comm_getdata(CPLID, &
        mpicom=mpicom_CPLID, nthreads=nthreads_CPLID)

   call seq_infodata_getData(infodata,       &
        drv_threading=drv_threading,         &
        atm_present=atm_present,             &
        lnd_present=lnd_present,             &
        rof_present=rof_present,             &
        ice_present=ice_present,             &
        ocn_present=ocn_present,             &
        glc_present=glc_present,             &
        wav_present=wav_present,             & 
        atm_prognostic=atm_prognostic,       &
        lnd_prognostic=lnd_prognostic,       &
        ice_prognostic=ice_prognostic,       &
        ocn_prognostic=ocn_prognostic,       &
        ocnrof_prognostic=ocnrof_prognostic, &
        glc_prognostic=glc_prognostic,       &
        wav_prognostic=wav_prognostic,       &
        atm_nx=atm_nx,  atm_ny=atm_ny,       &
        lnd_nx=lnd_nx,  lnd_ny=lnd_ny,       &
        rof_nx=rof_nx,  rof_ny=rof_ny,       &
        ice_nx=ice_nx,  ice_ny=ice_ny,       &
        glc_nx=glc_nx,  glc_ny=glc_ny,       &
        wav_nx=wav_nx,  wav_ny=wav_ny,       &
        ocn_nx=ocn_nx,  ocn_ny=ocn_ny,       &
        cpl_cdf64=cdf64 )

   ! Get current date from clock needed to label the histavg pointer file

   call seq_timemgr_EClockGetData( EClock_d,  curr_ymd=curr_ymd,  curr_tod=curr_tod,  &
        start_ymd=start_ymd,  start_tod=start_tod,  curr_time=curr_time,  prev_time=prev_time,  &
        calendar=calendar)

   if (first_call) then
      if (atm_present) then
         do iidx = 1,  num_inst_atm
            c2x => component_get_c2x_cx(atm(iidx))
            lsize = mct_aVect_lsize(c2x)
            call mct_aVect_init(a2x_ax_avg(iidx), c2x, lsize)
            call mct_aVect_zero(a2x_ax_avg(iidx))

            x2c => component_get_x2c_cx(atm(iidx))
            lsize = mct_aVect_lsize(x2c)
            call mct_aVect_init(x2a_ax_avg(iidx), x2c, lsize)
            call mct_aVect_zero(x2a_ax_avg(iidx))
         enddo
      endif
      if (lnd_present) then
         do iidx = 1,  num_inst_lnd
            c2x => component_get_c2x_cx(lnd(iidx))
            lsize = mct_aVect_lsize(c2x)
            call mct_aVect_init(l2x_lx_avg(iidx), c2x, lsize)
            call mct_aVect_zero(l2x_lx_avg(iidx))

            x2c => component_get_x2c_cx(lnd(iidx))
            lsize = mct_aVect_lsize(x2c)
            call mct_aVect_init(x2l_lx_avg(iidx), x2c, lsize)
            call mct_aVect_zero(x2l_lx_avg(iidx))
         enddo
      endif
      if (rof_present .and. ocnrof_prognostic) then
         do iidx = 1,  num_inst_rof
            c2x => component_get_c2x_cx(rof(iidx))
            lsize = mct_aVect_lsize(c2x)
            call mct_aVect_init(r2x_rx_avg(iidx), c2x, lsize)
            call mct_aVect_zero(r2x_rx_avg(iidx))

            x2c => component_get_x2c_cx(rof(iidx))
            lsize = mct_aVect_lsize(x2c)
            call mct_aVect_init(x2r_rx_avg(iidx), x2c, lsize)
            call mct_aVect_zero(x2r_rx_avg(iidx))
         enddo
      endif
      if (ocn_present) then
         do iidx = 1,  num_inst_ocn
            c2x => component_get_c2x_cx(ocn(iidx))
            lsize = mct_aVect_lsize(c2x)
            call mct_aVect_init(o2x_ox_avg(iidx), c2x, lsize)
            call mct_aVect_zero(o2x_ox_avg(iidx))

            x2c => component_get_x2c_cx(ocn(iidx))
            lsize = mct_aVect_lsize(x2c)
            call mct_aVect_init(x2o_ox_avg(iidx), x2c, lsize)
            call mct_aVect_zero(x2o_ox_avg(iidx))
         enddo
      endif
      if (ice_present) then
         do iidx = 1,  num_inst_ice
            c2x => component_get_c2x_cx(ice(iidx))
            lsize = mct_aVect_lsize(c2x)
            call mct_aVect_init(i2x_ix_avg(iidx), c2x, lsize)
            call mct_aVect_zero(i2x_ix_avg(iidx))

            x2c => component_get_x2c_cx(ice(iidx))
            lsize = mct_aVect_lsize(x2c)
            call mct_aVect_init(x2i_ix_avg(iidx), x2c, lsize)
            call mct_aVect_zero(x2i_ix_avg(iidx))
         enddo
      endif
      if (glc_present) then
         do iidx = 1,  num_inst_glc
            c2x => component_get_c2x_cx(glc(iidx))
            lsize = mct_aVect_lsize(c2x)
            call mct_aVect_init(g2x_gx_avg(iidx), c2x, lsize)
            call mct_aVect_zero(g2x_gx_avg(iidx))

            x2c => component_get_x2c_cx(glc(iidx))
            lsize = mct_aVect_lsize(x2c)
            call mct_aVect_init(x2g_gx_avg(iidx), x2c, lsize)
            call mct_aVect_zero(x2g_gx_avg(iidx))
         enddo
      endif
      if (wav_present) then
         do iidx = 1,  num_inst_wav
            c2x => component_get_c2x_cx(wav(iidx))
            lsize = mct_aVect_lsize(c2x)
            call mct_aVect_init(w2x_wx_avg(iidx), c2x, lsize)
            call mct_aVect_zero(w2x_wx_avg(iidx))

            x2c => component_get_x2c_cx(wav(iidx))
            lsize = mct_aVect_lsize(x2c)
            call mct_aVect_init(x2w_wx_avg(iidx), x2c, lsize)
            call mct_aVect_zero(x2w_wx_avg(iidx))
         enddo
      endif
      cnt = 0
      tbnds(1) = prev_time
      first_call = .false.
   endif

   if (.not.write_now) then
      cnt = cnt + 1
      if (atm_present) then
         do iidx = 1,  num_inst_atm
            c2x => component_get_c2x_cx(atm(iidx))
            x2c => component_get_x2c_cx(atm(iidx))
            a2x_ax_avg(iidx)%rAttr = a2x_ax_avg(iidx)%rAttr + c2x%rAttr
            x2a_ax_avg(iidx)%rAttr = x2a_ax_avg(iidx)%rAttr + x2c%rAttr
         enddo
      endif
      if (lnd_present) then
         do iidx = 1,  num_inst_lnd
            c2x => component_get_c2x_cx(lnd(iidx))
            x2c => component_get_x2c_cx(lnd(iidx))
            l2x_lx_avg(iidx)%rAttr = l2x_lx_avg(iidx)%rAttr + c2x%rAttr
            x2l_lx_avg(iidx)%rAttr = x2l_lx_avg(iidx)%rAttr + x2c%rAttr
         enddo
      endif
      if (rof_present .and. ocnrof_prognostic) then
         do iidx = 1,  num_inst_rof
            c2x => component_get_c2x_cx(rof(iidx))
            x2c => component_get_x2c_cx(rof(iidx))
            r2x_rx_avg(iidx)%rAttr = r2x_rx_avg(iidx)%rAttr + c2x%rAttr
            x2r_rx_avg(iidx)%rAttr = x2r_rx_avg(iidx)%rAttr + x2c%rAttr
         enddo
      endif
      if (ocn_present) then
         do iidx = 1,  num_inst_ocn
            c2x => component_get_c2x_cx(ocn(iidx))
            x2c => component_get_x2c_cx(ocn(iidx))
            o2x_ox_avg(iidx)%rAttr = o2x_ox_avg(iidx)%rAttr + c2x%rAttr
            x2o_ox_avg(iidx)%rAttr = x2o_ox_avg(iidx)%rAttr + x2c%rAttr
         enddo
      endif
      if (ice_present) then
         do iidx = 1,  num_inst_ice
            c2x => component_get_c2x_cx(ice(iidx))
            x2c => component_get_x2c_cx(ice(iidx))
            i2x_ix_avg(iidx)%rAttr = i2x_ix_avg(iidx)%rAttr + c2x%rAttr
            x2i_ix_avg(iidx)%rAttr = x2i_ix_avg(iidx)%rAttr + x2c%rAttr
         enddo
      endif
      if (glc_present) then
         do iidx = 1,  num_inst_glc
            c2x => component_get_c2x_cx(glc(iidx))
            x2c => component_get_x2c_cx(glc(iidx))
            g2x_gx_avg(iidx)%rAttr = g2x_gx_avg(iidx)%rAttr + c2x%rAttr
            x2g_gx_avg(iidx)%rAttr = x2g_gx_avg(iidx)%rAttr + x2c%rAttr
         enddo
      endif
      if (wav_present) then
         do iidx = 1,  num_inst_wav
            c2x => component_get_c2x_cx(wav(iidx))
            x2c => component_get_x2c_cx(wav(iidx))
            w2x_wx_avg(iidx)%rAttr = w2x_wx_avg(iidx)%rAttr + c2x%rAttr
            x2w_wx_avg(iidx)%rAttr = x2w_wx_avg(iidx)%rAttr + x2c%rAttr
         enddo
      endif

   else

      cnt = cnt + 1
      tbnds(2) = curr_time
      if (atm_present) then
         do iidx = 1,  num_inst_atm
            c2x => component_get_c2x_cx(atm(iidx))
            x2c => component_get_x2c_cx(atm(iidx))
            a2x_ax_avg(iidx)%rAttr = (a2x_ax_avg(iidx)%rAttr + c2x%rAttr) / (cnt * 1.0_r8)
            x2a_ax_avg(iidx)%rAttr = (x2a_ax_avg(iidx)%rAttr + x2c%rAttr) / (cnt * 1.0_r8)
         enddo
      endif
      if (lnd_present) then
         do iidx = 1,  num_inst_lnd
            c2x => component_get_c2x_cx(lnd(iidx))
            x2c => component_get_x2c_cx(lnd(iidx))
            l2x_lx_avg(iidx)%rAttr = (l2x_lx_avg(iidx)%rAttr + c2x%rAttr) / (cnt * 1.0_r8)
            x2l_lx_avg(iidx)%rAttr = (x2l_lx_avg(iidx)%rAttr + x2c%rAttr) / (cnt * 1.0_r8)
         enddo
      endif
      if (rof_present .and. ocnrof_prognostic) then
         do iidx = 1,  num_inst_rof
            c2x => component_get_c2x_cx(rof(iidx))
            x2c => component_get_x2c_cx(rof(iidx))
            r2x_rx_avg(iidx)%rAttr = (r2x_rx_avg(iidx)%rAttr + c2x%rAttr) / (cnt * 1.0_r8)
            x2r_rx_avg(iidx)%rAttr = (x2r_rx_avg(iidx)%rAttr + x2c%rAttr) / (cnt * 1.0_r8)
         enddo
      endif
      if (ocn_present) then
         do iidx = 1,  num_inst_ocn
            c2x => component_get_c2x_cx(ocn(iidx))
            x2c => component_get_x2c_cx(ocn(iidx))
            o2x_ox_avg(iidx)%rAttr = (o2x_ox_avg(iidx)%rAttr + c2x%rAttr) / (cnt * 1.0_r8)
            x2o_ox_avg(iidx)%rAttr = (x2o_ox_avg(iidx)%rAttr + x2c%rAttr) / (cnt * 1.0_r8)
         enddo
      endif
      if (ice_present) then
         do iidx = 1,  num_inst_ice
            c2x => component_get_c2x_cx(ice(iidx))
            x2c => component_get_x2c_cx(ice(iidx))
            i2x_ix_avg(iidx)%rAttr = (i2x_ix_avg(iidx)%rAttr + c2x%rAttr) / (cnt * 1.0_r8)
            x2i_ix_avg(iidx)%rAttr = (x2i_ix_avg(iidx)%rAttr + x2c%rAttr) / (cnt * 1.0_r8)
         enddo
      endif
      if (glc_present) then
         do iidx = 1,  num_inst_glc
            c2x => component_get_c2x_cx(glc(iidx))
            x2c => component_get_x2c_cx(glc(iidx))
            g2x_gx_avg(iidx)%rAttr = (g2x_gx_avg(iidx)%rAttr + c2x%rAttr) / (cnt * 1.0_r8)
            x2g_gx_avg(iidx)%rAttr = (x2g_gx_avg(iidx)%rAttr + x2c%rAttr) / (cnt * 1.0_r8)
         enddo
      endif
      if (wav_present) then
         do iidx = 1,  num_inst_wav
            c2x => component_get_c2x_cx(wav(iidx))
            x2c => component_get_x2c_cx(wav(iidx))
            w2x_wx_avg(iidx)%rAttr = (w2x_wx_avg(iidx)%rAttr + c2x%rAttr) / (cnt * 1.0_r8)
            x2w_wx_avg(iidx)%rAttr = (x2w_wx_avg(iidx)%rAttr + x2c%rAttr) / (cnt * 1.0_r8)
         enddo
      endif

      call seq_infodata_GetData( infodata,  case_name=case_name)
      call seq_timemgr_EClockGetData( EClock_d,  prev_ymd=prev_ymd,  prev_tod=prev_tod)
      if (seq_timemgr_histavg_type == seq_timemgr_type_nyear) then
         call shr_cal_date2ymd(prev_ymd, yy, mm, dd)
         write(hist_file, "(2a, i4.4, a)") &
            trim(case_name),  '.cpl.ha.',  yy, '.nc'
      elseif (seq_timemgr_histavg_type == seq_timemgr_type_nmonth) then
         call shr_cal_date2ymd(prev_ymd, yy, mm, dd)
         write(hist_file, "(2a, i4.4, a, i2.2, a)") &
            trim(case_name),  '.cpl.ha.',  yy, '-', mm, '.nc'
      elseif (seq_timemgr_histavg_type == seq_timemgr_type_nday) then
         call shr_cal_date2ymd(prev_ymd, yy, mm, dd)
         write(hist_file, "(2a, i4.4, a, i2.2, a, i2.2, a)") &
            trim(case_name),  '.cpl.ha.',  yy, '-', mm, '-', dd, '.nc'
      else
         call shr_cal_date2ymd(curr_ymd, yy, mm, dd)
         write(hist_file, "(2a, i4.4, a, i2.2, a, i2.2, a, i5.5, a)") &
            trim(case_name),  '.cpl.ha.',  yy, '-', mm, '-', dd, '-', curr_tod, '.nc'
      endif

      time_units = 'days since ' &
           // seq_io_date2yyyymmdd(start_ymd) // ' ' // seq_io_sec2hms(start_tod)

      if (iamin_CPLID) then

         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         call seq_io_wopen(hist_file, clobber=.true., cdf64=cdf64)

         ! loop twice,  first time write header,  second time write data for perf

         do fk = 1, 2
            if (fk == 1) then
               whead = .true.
               wdata = .false.
            elseif (fk == 2) then
               whead = .false.
               wdata = .true.
               call seq_io_enddef(hist_file)
            else
               call shr_sys_abort('seq_hist_writeavg fk illegal')
            end if

            avg_time = 0.5_r8 * (tbnds(1) + tbnds(2))
            !---------- tcx nov 2011 tbnds of same values causes problems in ferret
            if (tbnds(1) >= tbnds(2)) then
               call seq_io_write(hist_file, &
                    time_units=time_units, time_cal=calendar, time_val=avg_time, &
                    whead=whead, wdata=wdata)
            else
               call seq_io_write(hist_file, &
                    time_units=time_units, time_cal=calendar, time_val=avg_time, &
                    whead=whead, wdata=wdata, tbnds=tbnds)
            endif
            if (atm_present) then
               gsmap => component_get_gsmap_cx(atm(1))
               dom   => component_get_dom_cx(atm(1))
               call seq_io_write(hist_file, gsmap, dom%data, 'dom_ax',  &
                    nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata, pre='doma')

               call seq_io_write(hist_file, gsmap, x2a_ax_avg, 'x2a_ax',  & 
                    nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata,  &
                    pre='x2aavg', tavg=.true.)
               call seq_io_write(hist_file, gsmap, a2x_ax_avg, 'a2x_ax',  & 
                    nx=atm_nx, ny=atm_ny, nt=1, whead=whead, wdata=wdata,  &
                    pre='a2xavg', tavg=.true.)
            endif
            if (lnd_present) then
               gsmap => component_get_gsmap_cx(lnd(1))
               dom   => component_get_dom_cx(lnd(1))
               call seq_io_write(hist_file, gsmap, dom%data, 'dom_lx',  &
                    nx=lnd_nx, ny=lnd_ny, nt=1, whead=whead, wdata=wdata, pre='doml')
               call seq_io_write(hist_file, gsmap, l2x_lx_avg, 'l2x_lx',  &
                    nx=lnd_nx, ny=lnd_ny, nt=1, whead=whead, wdata=wdata,  &
                    pre='l2xavg', tavg=.true.)
               call seq_io_write(hist_file, gsmap, x2l_lx_avg, 'x2l_lx',  &
                    nx=lnd_nx, ny=lnd_ny, nt=1, whead=whead, wdata=wdata,  &
                    pre='x2lavg', tavg=.true.)
            endif

            if (rof_present .and. ocnrof_prognostic) then
               gsmap => component_get_gsmap_cx(rof(1))
               dom   => component_get_dom_cx(rof(1))
               call seq_io_write(hist_file, gsmap, dom%data, 'dom_rx',  &
                    nx=rof_nx, ny=rof_ny, nt=1, whead=whead, wdata=wdata, pre='domr')
               call seq_io_write(hist_file, gsmap, r2x_rx_avg, 'r2x_rx',  &
                    nx=rof_nx, ny=rof_ny, nt=1, whead=whead, wdata=wdata,  &
                    pre='r2xavg', tavg=.true.)
               call seq_io_write(hist_file, gsmap, x2r_rx_avg, 'x2r_rx',  &
                    nx=rof_nx, ny=rof_ny, nt=1, whead=whead, wdata=wdata,  &
                    pre='x2ravg', tavg=.true.)
            endif
            if (ocn_present) then
               gsmap => component_get_gsmap_cx(ocn(1))
               dom   => component_get_dom_cx(ocn(1))
               call seq_io_write(hist_file, gsmap, dom%data, 'dom_ox',  &
                    nx=ocn_nx, ny=ocn_ny, nt=1, whead=whead, wdata=wdata, pre='domo')
               call seq_io_write(hist_file, gsmap, o2x_ox_avg, 'o2x_ox',  &
                    nx=ocn_nx, ny=ocn_ny, nt=1, whead=whead, wdata=wdata,  &
                    pre='o2xavg', tavg=.true.)
               call seq_io_write(hist_file, gsmap, x2o_ox_avg, 'x2o_ox',  &
                    nx=ocn_nx, ny=ocn_ny, nt=1, whead=whead, wdata=wdata,  &
                    pre='x2oavg', tavg=.true.)
            endif
            if (ice_present) then
               gsmap => component_get_gsmap_cx(ice(1))
               dom   => component_get_dom_cx(ice(1))
               call seq_io_write(hist_file, gsmap, dom%data, 'dom_ix',  &
                                 nx=ice_nx, ny=ice_ny, nt=1, whead=whead, wdata=wdata, pre='domi')
               call seq_io_write(hist_file, gsmap, i2x_ix_avg, 'i2x_ix',  &
                                 nx=ice_nx, ny=ice_ny, nt=1, whead=whead, wdata=wdata,  &
                                 pre='i2xavg', tavg=.true.)
               call seq_io_write(hist_file, gsmap, x2i_ix_avg, 'x2i_ix',  &
                                 nx=ice_nx, ny=ice_ny, nt=1, whead=whead, wdata=wdata,  &
                                 pre='x2iavg', tavg=.true.)
            endif
            if (glc_present) then
               gsmap => component_get_gsmap_cx(glc(1))
               dom   => component_get_dom_cx(glc(1))
               call seq_io_write(hist_file, gsmap, dom%data, 'dom_gx',  &
                                 nx=glc_nx, ny=glc_ny, nt=1, whead=whead, wdata=wdata, pre='domg')
               call seq_io_write(hist_file, gsmap, g2x_gx_avg, 'g2x_gx',  &
                                 nx=glc_nx, ny=glc_ny, nt=1, whead=whead, wdata=wdata,  &
                                 pre='g2xavg', tavg=.true.)
               call seq_io_write(hist_file, gsmap, x2g_gx_avg, 'x2g_gx',  &
                                 nx=glc_nx, ny=glc_ny, nt=1, whead=whead, wdata=wdata,  &
                                 pre='x2gavg', tavg=.true.)
            endif
            if (wav_present) then
               gsmap => component_get_gsmap_cx(wav(1))
               dom   => component_get_dom_cx(wav(1))
               call seq_io_write(hist_file, gsmap, dom%data, 'dom_wx',  &
                                 nx=wav_nx, ny=wav_ny, nt=1, whead=whead, wdata=wdata, pre='domw')
               call seq_io_write(hist_file, gsmap, w2x_wx_avg, 'w2x_wx',  &
                                 nx=wav_nx, ny=wav_ny, nt=1, whead=whead, wdata=wdata,  &
                                 pre='w2xavg', tavg=.true.)
               call seq_io_write(hist_file, gsmap, x2w_wx_avg, 'x2w_wx',  &
                                 nx=wav_nx, ny=wav_ny, nt=1, whead=whead, wdata=wdata,  &
                                 pre='x2wavg', tavg=.true.)
            endif
         enddo

         call seq_io_close(hist_file)
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

         if (atm_present) then
            do iidx = 1,  num_inst_atm
               call mct_aVect_zero(a2x_ax_avg(iidx))
               call mct_aVect_zero(x2a_ax_avg(iidx))
            enddo
         endif
         if (lnd_present) then
            do iidx = 1,  num_inst_lnd
               call mct_aVect_zero(l2x_lx_avg(iidx))
               call mct_aVect_zero(x2l_lx_avg(iidx))
            enddo
         endif
         if (rof_present .and. ocnrof_prognostic) then
            do iidx = 1,  num_inst_rof
               call mct_aVect_zero(r2x_rx_avg(iidx))
               call mct_aVect_zero(x2r_rx_avg(iidx))
            enddo
         endif
         if (ocn_present) then
            do iidx = 1,  num_inst_ocn
               call mct_aVect_zero(o2x_ox_avg(iidx))
               call mct_aVect_zero(x2o_ox_avg(iidx))
            enddo
         endif
         if (ice_present) then
            do iidx = 1,  num_inst_ice
               call mct_aVect_zero(i2x_ix_avg(iidx))
               call mct_aVect_zero(x2i_ix_avg(iidx))
            enddo
         endif
         if (glc_present) then
            do iidx = 1,  num_inst_glc
               call mct_aVect_zero(g2x_gx_avg(iidx))
               call mct_aVect_zero(x2g_gx_avg(iidx))
            enddo
         endif 
         if (wav_present) then
            do iidx = 1,  num_inst_wav
               call mct_aVect_zero(w2x_wx_avg(iidx))
               call mct_aVect_zero(x2w_wx_avg(iidx))
            enddo
         endif
         cnt = 0
         tbnds(1) = curr_time

      endif
   endif

end subroutine seq_hist_writeavg

!===============================================================================

subroutine seq_hist_writeaux(infodata, EClock_d, comp, flow, aname, dname, &
     nx, ny, nt, write_now, flds, yr_offset)

   implicit none

   !--- arguments ---
   type (seq_infodata_type) , intent(inout)        :: infodata
   type(ESMF_Clock)         , intent(in)           :: EClock_d  ! driver clock
   type(component_type)     , intent(in)           :: comp      ! component instance
   character(len=3)         , intent(in)           :: flow      ! 'x2c' or 'c2x'
   character(*)             , intent(in)           :: aname     ! avect name for hist file
   character(*)             , intent(in)           :: dname     ! domain name for hist file
   integer(IN)              , intent(in)           :: nx        ! 2d global size nx
   integer(IN)              , intent(in)           :: ny        ! 2d global size ny
   integer(IN)              , intent(in)           :: nt        ! number of time samples per file
   logical                  , optional, intent(in) :: write_now ! write a sample now, if not used, write every call
   character(*)             , optional, intent(in) :: flds      ! list of fields to write
   integer                  , optional, intent(in) :: yr_offset ! offset to apply to current year when generating file name

   !--- local ---
   type(mct_gGrid), pointer :: dom
   type(mct_avect), pointer :: av
   type(mct_gsMap), pointer :: gsmap
   character(CL)            :: case_name  ! case name
   integer(IN)              :: curr_ymd   ! Current date YYYYMMDD
   integer(IN)              :: curr_tod   ! Current time-of-day (s)
   integer(IN)              :: start_ymd  ! Starting date YYYYMMDD
   integer(IN)              :: start_tod  ! Starting time-of-day (s)
   real(r8)                 :: curr_time  ! Time interval since reference time
   real(r8)                 :: prev_time  ! Time interval since reference time
   real(r8)                 :: avg_time   ! Average time for time average
   integer(IN)              :: yy, mm, dd ! year,  month,  day
   integer(IN)              :: n, fk, fk1 ! index
   character(CL)            :: time_units ! units of time variable
   character(CL)            :: calendar   ! calendar type
   integer(IN)              :: samples_per_file
   integer(IN)              :: lsize      ! local size of an aVect
   logical                  :: first_call
   integer(IN)              :: found = -10
   logical                  :: useavg
   logical                  :: lwrite_now     
   logical                  :: whead, wdata  ! for writing restart/history cdf files
   real(r8)                 :: tbnds(2)

   integer(IN), parameter   :: maxout = 20
   integer(IN)       , save :: ntout = 0
   character(CS)     , save :: tname(maxout) = 'x1y2z3'
   integer(IN)       , save :: ncnt(maxout)  = -10
   character(CL)     , save :: hist_file(maxout)       ! local path to history filename
   type(mct_aVect)   , save :: avavg(maxout)           ! av accumulator if needed
   integer(IN)       , save :: avcnt(maxout) = 0       ! accumulator counter
   logical           , save :: fwrite(maxout) = .true. ! first write
   real(r8)          , save :: tbnds1(maxout)          ! first time_bnds
   real(r8)          , save :: tbnds2(maxout)          ! second time_bnds

   type(mct_aVect)         :: avflds                  ! non-avg av for a subset of fields

   real(r8), parameter :: c0 = 0.0_r8 ! zero

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! get required infodata
   !----------------------------------------------------------------------------
   iamin_CPLID  = seq_comm_iamin(CPLID)

   call seq_comm_getdata(GLOID, &
        mpicom=mpicom_GLOID, nthreads=nthreads_GLOID)
   call seq_comm_getdata(CPLID, &
        mpicom=mpicom_CPLID, nthreads=nthreads_CPLID)

   call seq_infodata_getData(infodata, &
        drv_threading=drv_threading,   &
        atm_present=atm_present,       &
        lnd_present=lnd_present,       &
        rof_present=rof_present,       &
        ice_present=ice_present,       &
        ocn_present=ocn_present,       &
        glc_present=glc_present,       &
        wav_present=wav_present,       &
        cpl_cdf64=cdf64 )

   lwrite_now = .true.
   useavg = .false.
   if (present(write_now)) then
      useavg = .true.
      lwrite_now = write_now
   endif
 
   call seq_timemgr_EClockGetData( EClock_d, &
        curr_ymd=curr_ymd,                   &
        curr_tod=curr_tod,                   &
        start_ymd=start_ymd,                 &
        start_tod=start_tod,                 &
        curr_time=curr_time,                 &
        prev_time=prev_time,                 &
        calendar=calendar)

   first_call = .true.
   do n = 1, ntout
      if (trim(tname(n)) == trim(aname)) then
         first_call = .false.
         found = n
      endif
   enddo

   if (first_call) then
      ntout = ntout + 1
      if (ntout > maxout) then
         write(logunit, *) 'write_history_writeaux maxout exceeded', ntout, maxout
         call shr_sys_abort()
      endif
      tname(ntout) = trim(aname)
      ncnt(ntout) = -10
      if (iamin_CPLID .and. useavg) then
         lsize = mct_aVect_lsize(av)
         call mct_aVect_init(avavg(ntout), av, lsize)
         call mct_aVect_zero(avavg(ntout))
         avcnt(ntout) = 0
      endif
      tbnds1(ntout) = prev_time
      found = ntout
   endif

   if (iamin_CPLID) then !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      dom   => component_get_dom_cx(comp)
      gsmap => component_get_gsmap_cx(comp)

      samples_per_file = nt      

      if (useavg) then
         if (lwrite_now) then
            avcnt(found) = avcnt(found) + 1
            avavg(found)%rAttr = (avavg(found)%rAttr + av%rAttr) / (avcnt(found) * 1.0_r8)
         else
            avcnt(found) = avcnt(found) + 1
            avavg(found)%rAttr = avavg(found)%rAttr + av%rAttr
         endif
      endif
      
      if (lwrite_now) then
         
         ncnt(found) = ncnt(found) + 1
         if (ncnt(found) < 1 .or. ncnt(found) > samples_per_file) ncnt(found) = 1
         
         time_units = 'days since ' &
              // seq_io_date2yyyymmdd(start_ymd) // ' ' // seq_io_sec2hms(start_tod)
         tbnds2(found) = curr_time
         
         if (ncnt(found) == 1) then
            fk1 = 1
            call seq_infodata_GetData( infodata,  case_name=case_name)
            call shr_cal_date2ymd(curr_ymd, yy, mm, dd)
            
            ! Adjust yyyy in file name by yr_offset,  if present
            ! For example,  for a field written once a year,  this will make it so the file
            ! with fields from year 1 has time stamp 0001-01-01 rather than 0002-01-01, 
            ! which can simplify later reading by a data model
            if (present(yr_offset)) then
               yy = yy + yr_offset
            end if
            
            write(hist_file(found), "(a, i4.4, a, i2.2, a, i2.2, a)") &
                 trim(case_name)//'.cpl.h'//trim(aname)//'.',  yy, '-', mm, '-', dd, '.nc'
         else
            fk1 = 2
         endif
         
      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
      if (fk1 == 1) then
         call seq_io_wopen(hist_file(found), clobber=.true., cdf64=cdf64)
      else
         call seq_io_wopen(hist_file(found), clobber=.false., cdf64=cdf64)
      endif

      ! loop twice,  first time write header,  second time write data for perf

      tbnds(1) = tbnds1(found)
      tbnds(2) = tbnds2(found)

      do fk = fk1, 2
         if (fk == 1) then
            whead = .true.
            wdata = .false.
         elseif (fk == 2) then
            whead = .false.
            wdata = .true.
         else
            call shr_sys_abort('seq_hist_writeaux fk illegal')
         end if

         if (present(flds)) then
            if (fk == fk1) then
               lsize = mct_aVect_lsize(av)
               call mct_aVect_init(avflds,  rList=flds,  lsize=lsize)
               call mct_aVect_zero(avflds)
            end if
         end if

         avg_time = 0.5_r8 * (tbnds(1) + tbnds(2))
         !------- tcx nov 2011 tbnds of same values causes problems in ferret
         if (tbnds(1) >= tbnds(2)) then
            call seq_io_write(hist_file(found), &
                           time_units=time_units, time_cal=calendar, time_val=avg_time, &
                           nt=ncnt(found), whead=whead, wdata=wdata)
         else
            call seq_io_write(hist_file(found), &
                           time_units=time_units, time_cal=calendar, time_val=avg_time, &
                           nt=ncnt(found), whead=whead, wdata=wdata, tbnds=tbnds)
         endif

         if (fwrite(found)) then
            call seq_io_write(hist_file(found), gsmap, dom%data, trim(dname),  &
                              nx=nx, ny=ny, whead=whead, wdata=wdata, fillval=c0, pre=trim(dname))
         endif

         if (useavg) then
            if (present(flds)) then
               call mct_aVect_copy(aVin=avavg(found),  aVout=avflds)
               call seq_io_write(hist_file(found),  gsmap,  avflds,  trim(aname),  &
                                 nx=nx,  ny=ny,  nt=ncnt(found),  whead=whead,  wdata=wdata,  &
                                 pre=trim(aname), tavg=.true., use_float=.true.)
            else
               call seq_io_write(hist_file(found),  gsmap,  avavg(found),  trim(aname),  &
                                 nx=nx,  ny=ny,  nt=ncnt(found),  whead=whead,  wdata=wdata,  &
                                 pre=trim(aname), tavg=.true.,  use_float=.true.)
            end if
         else if (present(flds)) then
            call mct_aVect_copy(aVin=av,  aVout=avflds)
            call seq_io_write(hist_file(found),  gsmap,  avflds,  trim(aname),  &
                              nx=nx, ny=ny, nt=ncnt(found), whead=whead, wdata=wdata, pre=trim(aname), &
                              use_float=.true.)
         else
            call seq_io_write(hist_file(found),  gsmap,  av,  trim(aname),  &
                              nx=nx, ny=ny, nt=ncnt(found), whead=whead, wdata=wdata, pre=trim(aname), &
                              use_float=.true.)
         endif
   
         if (present(flds)) then
            if (fk == 2) then
               call mct_aVect_clean(avflds)
            end if
         end if

         if (fk == 1) call seq_io_enddef(hist_file(found))
         if (fk == 2) then
            fwrite(found) = .false.
            if (useavg) then
               call mct_aVect_zero(avavg(found))
               avcnt(found) = 0
            endif
            tbnds1(found) = curr_time
         endif
      enddo

      call seq_io_close(hist_file(found))
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

   endif   ! lwrite_now

   endif   ! iamin_CPLID <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end subroutine seq_hist_writeaux

!===============================================================================

subroutine seq_hist_spewav(infodata, aname, gsmap, av, nx, ny, nt, write_now, flds)

   implicit none

   type(seq_infodata_type) , intent(in)    :: infodata
   character(*)    , intent(in)           :: aname     ! avect name for hist file
   type(mct_gsmap) , intent(in)           :: gsmap     ! gsmap
   type(mct_aVect) , intent(in)           :: av        ! avect
   integer(IN)     , intent(in)           :: nx        ! 2d global size nx
   integer(IN)     , intent(in)           :: ny        ! 2d global size ny
   integer(IN)     , intent(in)           :: nt        ! number of time samples per file
   logical         , intent(in), optional :: write_now ! write a sample now,  if not used,  write every call
   character(*)    , intent(in), optional :: flds      ! list of fields to write

   !--- local ---
   character(CL)           :: case_name         ! case name
   integer(IN)             :: n,fk,fk1          ! index
   integer(IN)             :: samples_per_file
   integer(IN)             :: lsize             ! local size of an aVect
   logical                 :: first_call
   integer(IN)             :: found = -10
   logical                 :: useavg
   logical                 :: lwrite_now     
   logical                 :: whead,wdata  ! for writing restart/history cdf files
   real(r8)                :: tbnds(2)

   integer(IN),parameter   :: maxout = 20
   integer(IN)       ,save :: ntout = 0
   character(CS)     ,save :: tname(maxout) = 'x1y2z3'
   integer(IN)       ,save :: ncnt(maxout)  = -10
   integer(IN)       ,save :: nfiles(maxout) = 0
   character(CL)     ,save :: hist_file(maxout)       ! local path to history filename
   type(mct_aVect)   ,save :: avavg(maxout)           ! av accumulator if needed
   integer(IN)       ,save :: avcnt(maxout) = 0       ! accumulator counter
   logical           ,save :: fwrite(maxout) = .true. ! first write

   type(mct_aVect)         :: avflds                  ! non-avg av for a subset of fields

   real(r8),parameter :: c0 = 0.0_r8 ! zero

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! get required infodata
   !----------------------------------------------------------------------------
   iamin_CPLID  = seq_comm_iamin(CPLID)

   call seq_comm_getdata(GLOID, mpicom=mpicom_GLOID, nthreads=nthreads_GLOID)
   call seq_comm_getdata(CPLID, mpicom=mpicom_CPLID, nthreads=nthreads_CPLID)

   call seq_infodata_getData(infodata, &
        drv_threading=drv_threading,   &
        atm_present=atm_present,       &
        lnd_present=lnd_present,       &
        rof_present=rof_present,       &
        ice_present=ice_present,       &
        ocn_present=ocn_present,       &
        glc_present=glc_present,       &
        wav_present=wav_present,       &
        cpl_cdf64=cdf64 )

   lwrite_now = .true.
   useavg = .false.
   if (present(write_now)) then
      useavg = .true.
      lwrite_now = write_now
   endif
 
   first_call = .true.
   do n = 1, ntout
      if (trim(tname(n)) == trim(aname)) then
         first_call = .false.
         found = n
      endif
   enddo

   if (first_call) then
      ntout = ntout + 1
      if (ntout > maxout) then
         write(logunit, *) 'write_history_spewAV maxout exceeded', ntout, maxout
         call shr_sys_abort()
      endif
      tname(ntout) = trim(aname)
      ncnt(ntout) = -10
      nfiles(ntout) = 0
      if (iamin_CPLID .and. useavg) then
         lsize = mct_aVect_lsize(av)
         call mct_aVect_init(avavg(ntout), av, lsize)
         call mct_aVect_zero(avavg(ntout))
         avcnt(ntout) = 0
      endif
      found = ntout
   endif

!  if (.not. iamin_CPLID) return
   if (iamin_CPLID) then !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   samples_per_file = nt      

   if (useavg) then
      if (lwrite_now) then
         avcnt(found) = avcnt(found) + 1
         avavg(found)%rAttr = (avavg(found)%rAttr + av%rAttr) / (avcnt(found) * 1.0_r8)
      else
         avcnt(found) = avcnt(found) + 1
         avavg(found)%rAttr = avavg(found)%rAttr + av%rAttr
      endif
   endif

   if (lwrite_now) then

      ncnt(found) = ncnt(found) + 1
      if (ncnt(found) < 1 .or. ncnt(found) > samples_per_file) then
         ncnt(found) = 1
         nfiles(found) = nfiles(found) + 1
      endif

      if (ncnt(found) == 1) then
         fk1 = 1
         call seq_infodata_GetData( infodata,  case_name=case_name)
         write(hist_file(found), "(a, i4.4, a)") &
            trim(case_name)//'.cpl.h'//trim(aname)//'.', nfiles(found), '.nc'
      else
         fk1 = 2
      endif

      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
      if (fk1 == 1) then
         call seq_io_wopen(hist_file(found), clobber=.true., cdf64=cdf64)
      else
         call seq_io_wopen(hist_file(found), clobber=.false., cdf64=cdf64)
      endif

      ! loop twice,  first time write header,  second time write data for perf

      do fk = fk1, 2
         if (fk == 1) then
            whead = .true.
            wdata = .false.
         elseif (fk == 2) then
            whead = .false.
            wdata = .true.
         else
            call shr_sys_abort('seq_hist_spewav fk illegal')
         end if

         if (present(flds)) then
            if (fk == fk1) then
               lsize = mct_aVect_lsize(av)
               call mct_aVect_init(avflds,  rList=flds,  lsize=lsize)
               call mct_aVect_zero(avflds)
            end if
         end if

         tbnds = real(ncnt(found), r8)
!------- tcx nov 2011 tbnds of same values causes problems in ferret
         if (tbnds(1) >= tbnds(2)) then
           call seq_io_write(hist_file(found), &
                             time_units='nstep', time_cal='nstep', time_val=real(ncnt(found), r8), &
                             nt=ncnt(found), whead=whead, wdata=wdata)
         else
           call seq_io_write(hist_file(found), &
                             time_units='nstep', time_cal='nstep', time_val=real(ncnt(found), r8), &
                             nt=ncnt(found), whead=whead, wdata=wdata, tbnds=tbnds)
         endif

         if (useavg) then
            if (present(flds)) then
               call mct_aVect_copy(aVin=avavg(found),  aVout=avflds)
               call seq_io_write(hist_file(found),  gsmap,  avflds,  trim(aname),  &
                                 nx=nx,  ny=ny,  nt=ncnt(found),  whead=whead,  wdata=wdata,  &
                                 pre=trim(aname), tavg=.true., use_float=.true.)
            else
               call seq_io_write(hist_file(found),  gsmap,  avavg(found),  trim(aname),  &
                                 nx=nx,  ny=ny,  nt=ncnt(found),  whead=whead,  wdata=wdata,  &
                                 pre=trim(aname), tavg=.true.,  use_float=.true.)
            end if
         else if (present(flds)) then
            call mct_aVect_copy(aVin=av,  aVout=avflds)
            call seq_io_write(hist_file(found),  gsmap,  avflds,  trim(aname),  &
                              nx=nx, ny=ny, nt=ncnt(found), whead=whead, wdata=wdata, pre=trim(aname), &
                              use_float=.true.)
         else
            call seq_io_write(hist_file(found),  gsmap,  av,  trim(aname),  &
                              nx=nx, ny=ny, nt=ncnt(found), whead=whead, wdata=wdata, pre=trim(aname), &
                              use_float=.true.)
         endif
   
         if (present(flds)) then
            if (fk == 2) then
               call mct_aVect_clean(avflds)
            end if
         end if

         if (fk == 1) call seq_io_enddef(hist_file(found))
         if (fk == 2) then
            fwrite(found) = .false.
            if (useavg) then
               call mct_aVect_zero(avavg(found))
               avcnt(found) = 0
            endif
         endif
      enddo

      call seq_io_close(hist_file(found))
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

   endif   ! lwrite_now

   endif   ! iamin_CPLID <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end subroutine seq_hist_spewav

!===============================================================================

end module seq_hist_mod
