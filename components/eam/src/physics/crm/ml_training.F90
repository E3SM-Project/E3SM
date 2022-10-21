module ml_training

   ! TODO:
   ! - add physics state
   ! - figure out how to handles prescribed aerosol, ozone, etc.
   ! - screen CRM variables from pbuf data
   ! - add calls in phys_run1() 
   ! - write "output" routines

   use shr_kind_mod,       only: r8 => shr_kind_r8
   use spmd_utils,         only: masterproc
   use constituents,       only: pcnst
   use ioFileMod
   use cam_abortutils,     only: endrun
   use camsrfexch,         only: cam_in_t, cam_out_t
   use cam_logfile,        only: iulog
   use pio,                only: file_desc_t, io_desc_t, var_desc_t, &
                                 pio_double, pio_int, pio_noerr, &
                                 pio_seterrorhandling, pio_bcast_error, &
                                 pio_inq_varid, pio_def_var, pio_def_dim, &
                                 pio_put_var, pio_get_var
   use perf_mod,           only: t_startf, t_stopf

   implicit none
   private
   save
   
   ! Public interfaces
   public :: init_ml_training_input
   public :: write_ml_training_input   ! write physics input data for ML training
   ! public :: init_ml_training_output
   ! public :: write_ml_training_output  ! write physics output data for ML verification
   
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
   subroutine init_ml_training_input ( File, pbuf2d)
      use physics_buffer,      only: pbuf_init_restart, physics_buffer_desc
      use ppgrid,              only: pver, pverp, pcols
      use chemistry,           only: chem_init_restart
      use prescribed_ozone,    only: init_prescribed_ozone_restart
      use prescribed_ghg,      only: init_prescribed_ghg_restart
      use prescribed_aero,     only: init_prescribed_aero_restart
      use prescribed_volcaero, only: init_prescribed_volcaero_restart
      use cam_grid_support,    only: cam_grid_write_attr, cam_grid_id
      use cam_grid_support,    only: cam_grid_header_info_t
      use cam_pio_utils,       only: cam_pio_def_dim
      use phys_control,        only: phys_getopts
      !-------------------------------------------------------------------------
      ! Input arguments
      type(file_desc_t), intent(inout) :: file
      type(physics_buffer_desc), pointer :: pbuf2d(:,:)
      !-------------------------------------------------------------------------
      ! Local workspace
      integer                      :: hdimcnt, ierr, i, vsize
      integer                      :: dimids(4)
      integer, allocatable         :: hdimids(:)
      integer                      :: ndims, pver_id, pverp_id
      integer                      :: kiss_seed_dim
      type(cam_grid_header_info_t) :: info
      !-------------------------------------------------------------------------
      call pio_seterrorhandling(File, PIO_BCAST_ERROR)
      
      ! Probably should have the grid write this out (?)
      call cam_grid_write_attr(File, cam_grid_id('physgrid'), info)
      hdimcnt = info%num_hdims()

      do i = 1, hdimcnt
         dimids(i) = info%get_hdimid(i)
      end do
      allocate(hdimids(hdimcnt))
      hdimids(1:hdimcnt) = dimids(1:hdimcnt)

      call cam_pio_def_dim(File, 'lev', pver, pver_id, existOK=.true.)
      call cam_pio_def_dim(File, 'ilev', pverp, pverp_id, existOK=.true.)

      ndims=hdimcnt

      ndims=hdimcnt+1

      call pbuf_init_restart(File, pbuf2d)

      call chem_init_restart(File)

      call init_prescribed_ozone_restart(File)
      call init_prescribed_ghg_restart(File)
      call init_prescribed_aero_restart(File)
      call init_prescribed_volcaero_restart(File)

      call cam_pio_def_dim(File, 'pcnst', pcnst, dimids(hdimcnt+1), existOK=.true.)

      ierr = pio_def_var(File, 'FLWDS', pio_double, hdimids, flwds_desc)
      ierr = pio_def_var(File, 'SOLS', pio_double, hdimids, sols_desc)
      ierr = pio_def_var(File, 'SOLL', pio_double, hdimids, soll_desc)
      ierr = pio_def_var(File, 'SOLSD', pio_double, hdimids, solsd_desc)
      ierr = pio_def_var(File, 'SOLLD', pio_double, hdimids, solld_desc)

      ! ierr = pio_def_var(File, 'BCPHIDRY', pio_double, hdimids, bcphidry_desc)
      ! ierr = pio_def_var(File, 'BCPHODRY', pio_double, hdimids, bcphodry_desc)
      ! ierr = pio_def_var(File, 'OCPHIDRY', pio_double, hdimids, ocphidry_desc)
      ! ierr = pio_def_var(File, 'OCPHODRY', pio_double, hdimids, ocphodry_desc)
      ! ierr = pio_def_var(File, 'DSTDRY1',  pio_double, hdimids, dstdry1_desc)
      ! ierr = pio_def_var(File, 'DSTDRY2',  pio_double, hdimids, dstdry2_desc)
      ! ierr = pio_def_var(File, 'DSTDRY3',  pio_double, hdimids, dstdry3_desc)
      ! ierr = pio_def_var(File, 'DSTDRY4',  pio_double, hdimids, dstdry4_desc)

      ! cam_import variables -- write the constituent surface fluxes as individual 2D arrays
      ! rather than as a single variable with a pcnst dimension.  Note that the cflx components
      ! are only needed for those constituents that are not passed to the coupler.  The restart
      ! for constituents passed through the coupler are handled by the .rs. restart file.  But
      ! we don't currently have a mechanism to know whether the constituent is handled by the
      ! coupler or not, so we write all of cflx to the CAM restart file.
      do i = 1, pcnst
       write(num,'(i4.4)') i
       ierr = pio_def_var(File, 'CFLX'//num,  pio_double, hdimids, cflx_desc(i))
      end do
      ! Add LHF and SHF to restart file to fix non-BFB restart issue due to qneg4 correction at the restart time step
      ierr = pio_def_var(File, 'SHF',  pio_double, hdimids, shf_desc)
      ! ierr = pio_def_var(File, 'LHF',  pio_double, hdimids, lhf_desc)

   end subroutine init_ml_training_input
   !------------------------------------------------------------------------------------------------
   subroutine write_ml_training_input (File, cam_in, cam_out, pbuf2d)
      use physics_buffer,      only: physics_buffer_desc, pbuf_write_restart
      use phys_grid,           only: phys_decomp
      use ppgrid,              only: begchunk, endchunk, pcols, pverp
      use chemistry,           only: chem_write_restart
      use prescribed_ozone,    only: write_prescribed_ozone_restart
      use prescribed_ghg,      only: write_prescribed_ghg_restart
      use prescribed_aero,     only: write_prescribed_aero_restart
      use prescribed_volcaero, only: write_prescribed_volcaero_restart
      use cam_history_support, only: fillvalue
      use cam_grid_support,    only: cam_grid_id
      use cam_grid_support,    only: cam_grid_get_decomp, cam_grid_dimensions
      use cam_grid_support,    only: cam_grid_write_var
      use pio,                 only: pio_write_darray
      !-------------------------------------------------------------------------
      ! Input arguments
      type(file_desc_t), intent(inout) :: File
      type(cam_in_t),    intent(in)    :: cam_in(begchunk:endchunk)
      type(cam_out_t),   intent(in)    :: cam_out(begchunk:endchunk)
      type(physics_buffer_desc), pointer        :: pbuf2d(:,:)
      !-------------------------------------------------------------------------
      ! Local workspace
      type(io_desc_t), pointer :: iodesc
      real(r8):: tmpfield(pcols, begchunk:endchunk)
      integer :: i, m
      integer :: ierr
      integer :: physgrid
      integer :: gdims(3)
      integer :: nhdims
      !-------------------------------------------------------------------------
      ! Write grid vars
      call cam_grid_write_var(File, phys_decomp)

      ! Physics buffer
      call pbuf_write_restart(File, pbuf2d)

      physgrid = cam_grid_id('physgrid')
      call cam_grid_dimensions(physgrid, gdims(1:2), nhdims)

      ! data for chemistry
      call chem_write_restart(File)
      call write_prescribed_ozone_restart(File)
      call write_prescribed_ghg_restart(File)
      call write_prescribed_aero_restart(File)
      call write_prescribed_volcaero_restart(File)

      ! Get PIO decomposition
      call cam_grid_get_decomp(physgrid, (/pcols,endchunk-begchunk+1/), gdims(1:nhdims), pio_double, iodesc)

      ! Set missing portions of tmpfield - only do this once since cam_in/out vars all same shape
      do i = begchunk, endchunk
         if (cam_out(i)%ncol < pcols) tmpfield(cam_out(i)%ncol+1:, i) = fillvalue
      end do

      !-------------------------------------------------------------------------
      ! surface rad fluxes

      do i = begchunk, endchunk
         tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%flwds(:cam_out(i)%ncol)
      end do
      call pio_write_darray(File, flwds_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
         tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%sols(:cam_out(i)%ncol)
      end do
      call pio_write_darray(File, sols_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
         tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%soll(:cam_out(i)%ncol)
      end do
      call pio_write_darray(File, soll_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
         tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%solsd(:cam_out(i)%ncol)
      end do
      call pio_write_darray(File, solsd_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
         tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%solld(:cam_out(i)%ncol)
      end do
      call pio_write_darray(File, solld_desc, iodesc, tmpfield, ierr)

      !-------------------------------------------------------------------------
      ! surface constiuent fluxes

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%bcphidry(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(File, bcphidry_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%bcphodry(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(File, bcphodry_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%ocphidry(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(File, ocphidry_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%ocphodry(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(File, ocphodry_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%dstdry1(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(File, dstdry1_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%dstdry2(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(File, dstdry2_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%dstdry3(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(File, dstdry3_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_out(i)%ncol, i) = cam_out(i)%dstdry4(:cam_out(i)%ncol)
      ! end do
      ! call pio_write_darray(File, dstdry4_desc, iodesc, tmpfield, ierr)

      !-------------------------------------------------------------------------
      ! cam_in components

      do m = 1, pcnst
         do i = begchunk, endchunk
            tmpfield(:cam_in(i)%ncol, i) = cam_in(i)%cflx(:cam_in(i)%ncol, m)
         end do
         call pio_write_darray(File, cflx_desc(m), iodesc, tmpfield, ierr)
      end do

      do i = begchunk, endchunk
         tmpfield(:cam_in(i)%ncol, i) = cam_in(i)%shf(:cam_in(i)%ncol)
      end do
      call pio_write_darray(File, shf_desc, iodesc, tmpfield, ierr)

      ! do i = begchunk, endchunk
      !    tmpfield(:cam_in(i)%ncol, i) = cam_in(i)%lhf(:cam_in(i)%ncol)
      ! end do
      ! call pio_write_darray(File, lhf_desc, iodesc, tmpfield, ierr)
      
   end subroutine write_ml_training_input
   !------------------------------------------------------------------------------------------------
end module ml_training
