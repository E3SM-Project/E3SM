module ml_training

   ! TODO:
   ! - add physics state
   ! - figure out how to handles prescribed aerosol, ozone, etc. (see phys_timestep_init)
   ! - screen CRM variables from pbuf data
   ! - add calls in phys_run1() 
   ! - write "output" routines

   use shr_kind_mod,       only: r8 => shr_kind_r8
   use spmd_utils,         only: masterproc
   use constituents,       only: pcnst
   use iofileMod
   use cam_abortutils,     only: endrun
   use camsrfexch,         only: cam_in_t, cam_out_t
   use cam_logfile,        only: iulog
   use shr_sys_mod,        only: shr_sys_flush
   use pio,                only: file_desc_t, io_desc_t, var_desc_t, &
                                 pio_double, pio_int, pio_noerr, &
                                 pio_seterrorhandling, pio_bcast_error, &
                                 ! pio_inq_varid, pio_get_var,  &
                                 pio_def_var, pio_def_dim, pio_enddef, &
                                 pio_put_var, pio_write_darray
   use perf_mod,           only: t_startf, t_stopf

   implicit none
   private
   save
   
   ! Public interfaces
   public :: get_ml_input_filename
   public :: write_ml_training_input   ! write physics input data for ML training
   ! public :: write_ml_training_output  ! write physics output data for ML verification

   ! filename specifiers for master restart filename
   ! (%c = caseid, $y = year, $m = month, $d = day, $s = seconds in day, %t = number)
   integer, parameter :: filename_len = 256
   character(len=filename_len),public :: mli_filename_spec = '%c.eam.mli.%y-%m-%d-%s.nc'
   character(len=filename_len),public :: mlo_filename_spec = '%c.eam.mlo.%y-%m-%d-%s.nc'
   
   ! Private module data
   character(len=8)     :: num
   type(var_desc_t)     :: trefmxav_desc
   type(var_desc_t)     :: trefmnav_desc
   type(var_desc_t)     :: flwds_desc
   type(var_desc_t)     :: solld_desc
   type(var_desc_t)     :: sols_desc
   type(var_desc_t)     :: soll_desc
   type(var_desc_t)     :: solsd_desc

   ! type(var_desc_t)     :: bcphidry_desc
   ! type(var_desc_t)     :: bcphodry_desc
   ! type(var_desc_t)     :: ocphidry_desc
   ! type(var_desc_t)     :: ocphodry_desc
   ! type(var_desc_t)     :: dstdry1_desc
   ! type(var_desc_t)     :: dstdry2_desc
   ! type(var_desc_t)     :: dstdry3_desc
   ! type(var_desc_t)     :: dstdry4_desc

   type(var_desc_t)     :: cflx_desc(pcnst)
   type(var_desc_t)     :: lhf_desc
   type(var_desc_t)     :: shf_desc

CONTAINS
   !------------------------------------------------------------------------------------------------
   ! character(len=filename_len) function get_ml_input_filename() result(fname)
   function get_ml_input_filename(yr,mn,dy,sec) result(fname)
      use seq_timemgr_mod, only: seq_timemgr_EClockGetData
      use filenames,       only: interpret_filename_spec
      character(len=filename_len) :: fname  ! surface restart filename
      integer, intent(in) :: yr         ! current year
      integer, intent(in) :: mn         ! current month
      integer, intent(in) :: dy         ! current day
      integer, intent(in) :: sec        ! current time of day (sec)
      ! integer :: yr         ! current year
      ! integer :: mon        ! current month
      ! integer :: day        ! current day
      ! integer :: sec        ! current time of day (sec)

      ! call seq_timemgr_EClockGetData( EClock, curr_yr=yr, curr_mon=mon, &
      !                                 curr_day=day,curr_tod=sec)

      fname = interpret_filename_spec( mli_filename_spec, &
                                       yr_spec=yr, mon_spec=mn, &
                                       day_spec=dy, sec_spec=sec )

      ! return fname
   end function get_ml_input_filename
   !------------------------------------------------------------------------------------------------
   ! subroutine write_ml_training_input(file, cam_in, cam_out, pbuf2d)
   subroutine write_ml_training_input( pbuf2d, cam_in, cam_out, yr, mn, dy, sec )
#ifdef MMF_ML_TRAINING

      use physics_buffer,      only: pbuf_init_restart, physics_buffer_desc, pbuf_write_restart
      use time_manager,        only: timemgr_init_restart, timemgr_write_restart
      use ppgrid,              only: pver, pverp, pcols, begchunk, endchunk
      use phys_grid,           only: phys_decomp
      use chemistry,           only: chem_init_restart, chem_write_restart
      ! use prescribed_ozone,    only: init_prescribed_ozone_restart, write_prescribed_ozone_restart
      ! use prescribed_ghg,      only: init_prescribed_ghg_restart, write_prescribed_ghg_restart
      ! use prescribed_aero,     only: init_prescribed_aero_restart, write_prescribed_aero_restart
      ! use prescribed_volcaero, only: init_prescribed_volcaero_restart, write_prescribed_volcaero_restart
      use cam_grid_support,    only: cam_grid_id, cam_grid_header_info_t
      use cam_grid_support,    only: cam_grid_get_decomp, cam_grid_dimensions
      use cam_grid_support,    only: cam_grid_write_attr, cam_grid_write_var
      use cam_history_support, only: fillvalue
      use cam_pio_utils,       only: cam_pio_createfile, cam_pio_def_dim, cam_pio_closefile
      use phys_control,        only: phys_getopts
      use hycoef,              only: init_restart_hycoef
      ! use dimensions_mod,      only: np, ne, nelem
      use dyn_grid,            only: get_horiz_grid_d
      !-------------------------------------------------------------------------
      ! Input arguments
      ! type(file_desc_t), intent(inout) :: file
      type(physics_buffer_desc), pointer    :: pbuf2d(:,:)
      type(cam_in_t),            intent(in) :: cam_in(begchunk:endchunk)
      type(cam_out_t),           intent(in) :: cam_out(begchunk:endchunk)
      integer,                   intent(in) :: yr   ! Simulation year
      integer,                   intent(in) :: mn   ! Simulation month
      integer,                   intent(in) :: dy   ! Simulation day
      integer,                   intent(in) :: sec  ! Seconds into current simulation day
      !-------------------------------------------------------------------------
      ! Local workspace
      type(file_desc_t)            :: file
      integer                      :: hdim1_d,hdim2_d,ngcols ! number of physics grid columns
      integer                      :: hdimcnt
      integer                      :: dimids(4)
      integer, allocatable         :: hdimids(:)
      integer                      :: ndims, pver_id, pverp_id
      integer                      :: ncol_dimid
      ! integer                      :: lev_dimid(2)
      integer                      :: grid_id
      type(cam_grid_header_info_t) :: header_info ! A structure to hold the horz dims and coord info

      type(io_desc_t), pointer :: iodesc
      real(r8):: tmpfield(pcols, begchunk:endchunk)
      integer :: ierr, i, m
      integer :: physgrid
      integer :: gdims(3)
      integer :: nhdims
      !-------------------------------------------------------------------------
      ! Initialize file

      call cam_pio_createfile(file, trim(get_ml_input_filename(yr,mn,dy,sec)))

      grid_id = cam_grid_id('physgrid')
      call cam_grid_dimensions(grid_id, gdims(1:2), nhdims)
      call cam_grid_get_decomp(grid_id, (/pcols,endchunk-begchunk+1/), gdims(1:nhdims), pio_double, iodesc)

      write(iulog,222) grid_id, gdims(1), gdims(2), nhdims
      call shr_sys_flush(iulog)
222 format('WHDEBUG - grid_id:',i6,'  gdims(1):',i6,'  gdims(2):',i6,'  nhdims:',i6)

      ! ' lat:',f6.2,' lon:',f6.2,' k=',i2,'  t_vt_tend: ',E12.6,'  t_vt_in: ',E12.6,'  t_vt: ',E12.6)

      call cam_grid_write_attr(file, grid_id, header_info)
      call cam_grid_write_var(file, grid_id)
      
      hdimcnt = header_info%num_hdims()
      do i = 1, hdimcnt
         dimids(i) = header_info%get_hdimid(i)
      end do
      allocate(hdimids(hdimcnt))
      hdimids(1:hdimcnt) = dimids(1:hdimcnt)

      ! ndims = hdimcnt
      ! ndims = hdimcnt+1

      ! call cam_pio_def_dim(file, 'pcnst', pcnst, dimids(hdimcnt+1), existOK=.true.)
      !-------------------------------------------------------------------------

      ierr = pio_def_var(file, 'SHF',  pio_double, hdimids, shf_desc)
      ! ierr = pio_def_var(file, 'LHF',  pio_double, hdimids, lhf_desc)

      ierr = pio_enddef(file)
      !-------------------------------------------------------------------------

      ! call timemgr_write_restart(File)

      ! ! Physics buffer
      ! call pbuf_write_restart(file, pbuf2d)

      ! ! data for chemistry
      ! call chem_write_restart(file)
      ! call write_prescribed_ozone_restart(file)
      ! call write_prescribed_ghg_restart(file)
      ! call write_prescribed_aero_restart(file)
      ! call write_prescribed_volcaero_restart(file)

      ! Set missing portions of tmpfield - only do this once since cam_in/out vars all same shape
      do i = begchunk, endchunk
         if (cam_out(i)%ncol < pcols) tmpfield(cam_out(i)%ncol+1:, i) = fillvalue
      end do

      !-------------------------------------------------------------------------
      ! surface rad fluxes

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%flwds(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(file, flwds_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%sols(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(file, sols_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%soll(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(file, soll_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%solsd(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(file, solsd_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%solld(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(file, solld_desc, iodesc, tmpfield, ierr)

      !-------------------------------------------------------------------------
      ! surface constiuent fluxes

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%bcphidry(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(file, bcphidry_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%bcphodry(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(file, bcphodry_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%ocphidry(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(file, ocphidry_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%ocphodry(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(file, ocphodry_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%dstdry1(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(file, dstdry1_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%dstdry2(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(file, dstdry2_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%dstdry3(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(file, dstdry3_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%dstdry4(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(file, dstdry4_desc, iodesc, tmpfield, ierr)

      !-------------------------------------------------------------------------
      ! cam_in components

      ! do m = 1, pcnst
      !    do i = begchunk, endchunk
      !       tmpfield(:cam_in(i)%ncol, i) = cam_in(i)%cflx(:cam_in(i)%ncol, m)
      !    end do
      !    call pio_write_darray(file, cflx_desc(m), iodesc, tmpfield, ierr)
      ! end do

      do i = begchunk, endchunk
         tmpfield(:cam_in(i)%ncol, i) = cam_in(i)%shf(:cam_in(i)%ncol)
      end do
      call pio_write_darray(file, shf_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_in(i)%ncol, i) = cam_in(i)%lhf(:cam_in(i)%ncol)
      ! end do
      ! call pio_write_darray(file, lhf_desc, iodesc, tmpfield, ierr)

      !-------------------------------------------------------------------------
      ! close the file and clean-up

      if (allocated(hdimids)) deallocate(hdimids)

      call cam_pio_closefile(file)
      
#endif /* MMF_ML_TRAINING */
   end subroutine write_ml_training_input
   !------------------------------------------------------------------------------------------------
end module ml_training
