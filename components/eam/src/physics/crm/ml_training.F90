module ml_training

   ! TODO:
   ! - figure out how to handle prescribed aerosol, ozone, etc. (see phys_timestep_init)
   ! - allow data to be written at specific interval (i.e. not every time step)

   use shr_kind_mod,       only: r8 => shr_kind_r8
   use spmd_utils,         only: masterproc
   use constituents,       only: pcnst
   use ppgrid,             only: pver, pverp, pcols, begchunk, endchunk
   use cam_abortutils,     only: endrun
   use cam_history_support,only: fillvalue
   use camsrfexch,         only: cam_in_t, cam_out_t
   use cam_logfile,        only: iulog
   use shr_sys_mod,        only: shr_sys_flush
   use physics_types,      only: physics_state, physics_tend
   use physics_buffer,     only: physics_buffer_desc
   use pio,                only: file_desc_t, io_desc_t, var_desc_t, &
                                 pio_double, pio_int, pio_noerr, &
                                 pio_seterrorhandling, pio_bcast_error, &
                                 pio_def_var, pio_def_dim, pio_enddef, &
                                 pio_put_var, pio_write_darray, pio_closefile
   use perf_mod,           only: t_startf, t_stopf

   implicit none
   private
   save
   
   ! Public interfaces
   public :: write_ml_training

   ! internal private routines
   private :: get_ml_filename

   integer, parameter :: filename_len = 256

CONTAINS
   !------------------------------------------------------------------------------------------------
   function get_ml_filename(fspec_in,yr,mn,dy,sec) result(fname)
      ! return filename based on date/time info and fspec_in specifiers
      use seq_timemgr_mod, only: seq_timemgr_EClockGetData
      use filenames,       only: interpret_filename_spec
      character(len=4), intent(in) :: fspec_in
      integer,          intent(in) :: yr,mn,dy,sec   ! current year, month, day, and time of day (sec)
      character(len=filename_len) :: fspec_loc
      character(len=filename_len) :: fname  ! full file name to return
      ! (%c=caseid, $y=year, $m=month, $d=day, $s=sec in day, %t=number)
      fspec_loc = '%c.eam.'//trim(fspec_in)//'.%y-%m-%d-%s.nc'
      fname = interpret_filename_spec( fspec_loc, yr_spec=yr, mon_spec=mn, day_spec=dy, sec_spec=sec )
   end function get_ml_filename
   !------------------------------------------------------------------------------------------------
   subroutine write_ml_training( pbuf2d, phys_state, phys_tend, cam_in, cam_out, yr, mn, dy, sec, mode )
      use phys_grid,           only: phys_decomp
      use physics_buffer,      only: pbuf_init_restart_alt, pbuf_write_restart_alt
      use time_manager,        only: timemgr_init_restart, timemgr_write_restart
      use chemistry,           only: chem_init_restart, chem_write_restart
      use cam_grid_support,    only: cam_grid_id, cam_grid_header_info_t
      use cam_grid_support,    only: cam_grid_get_decomp, cam_grid_dimensions
      use cam_grid_support,    only: cam_grid_write_attr, cam_grid_write_var
      use cam_pio_utils,       only: cam_pio_createfile, cam_pio_def_dim, cam_pio_closefile
      use phys_control,        only: phys_getopts
      use hycoef,              only: init_restart_hycoef
      use dyn_grid,            only: get_horiz_grid_d
      !-------------------------------------------------------------------------
      ! Input arguments
      type(physics_buffer_desc), pointer    :: pbuf2d(:,:)
      type(physics_state),       intent(in) :: phys_state(begchunk:endchunk)
      type(physics_tend ),       intent(in) :: phys_tend (begchunk:endchunk)
      type(cam_in_t),            intent(in) :: cam_in    (begchunk:endchunk)
      type(cam_out_t),           intent(in) :: cam_out   (begchunk:endchunk)
      integer,                   intent(in) :: yr   ! Simulation year
      integer,                   intent(in) :: mn   ! Simulation month
      integer,                   intent(in) :: dy   ! Simulation day
      integer,                   intent(in) :: sec  ! Seconds into current simulation day
      integer,                   intent(in) :: mode ! used to select input / output modes
      !-------------------------------------------------------------------------
      ! Local workspace
      type(file_desc_t)            :: file
      integer                      :: pver_id, pverp_id, pcnst_id
      integer                      :: ncol_dimid
      integer                      :: grid_id
      type(cam_grid_header_info_t) :: header_info ! A structure to hold the horz dims and coord info

      integer :: ncol(begchunk:endchunk)                ! ncol value per chunk
      type(io_desc_t), pointer :: iodesc2d              ! 
      type(io_desc_t), pointer :: iodesc3d              ! 
      real(r8):: tmp2D(pcols, begchunk:endchunk)        ! temp variable for derived type data
      real(r8):: tmp3D(pcols, pver, begchunk:endchunk)  ! temp variable for derived type data

      integer :: ierr, i, m
      integer :: physgrid
      integer :: gdims(3)
      integer :: nhdims

      integer, parameter, dimension(1) :: dimids_hrz = (/1/)     ! horz dims only
      integer, parameter, dimension(2) :: dimids_3D1 = (/1,2/)   ! horz + pver
      integer, parameter, dimension(2) :: dimids_3D2 = (/1,3/)   ! horz + pverp

      character(len=8)     :: num      ! used for writing numeric charaters (i.e. constituent index)
      character(len=4)     :: fspec    ! string used after ".eam." in file name 
      logical              :: add_pbuf       = .false.
      logical              :: add_phys_state = .false.
      logical              :: add_phys_tend  = .false.
      logical              :: add_cam_in     = .false.
      logical              :: add_cam_out    = .false.

      ! file variable descriptions
      type(var_desc_t)     :: desc_flwds
      type(var_desc_t)     :: desc_solld
      type(var_desc_t)     :: desc_sols
      type(var_desc_t)     :: desc_soll
      type(var_desc_t)     :: desc_solsd

      type(var_desc_t)     :: state_desc_t
      type(var_desc_t)     :: state_desc_u
      type(var_desc_t)     :: state_desc_v
      type(var_desc_t)     :: state_desc_s
      type(var_desc_t)     :: state_desc_omega
      type(var_desc_t)     :: state_desc_pmid
      type(var_desc_t)     :: state_desc_pmiddry
      type(var_desc_t)     :: state_desc_pdel
      type(var_desc_t)     :: state_desc_pdeldry
      type(var_desc_t)     :: state_desc_rpdel
      type(var_desc_t)     :: state_desc_rpdeldry
      type(var_desc_t)     :: state_desc_lnpmid
      type(var_desc_t)     :: state_desc_lnpmiddry
      type(var_desc_t)     :: state_desc_exner
      type(var_desc_t)     :: state_desc_zm
      type(var_desc_t)     :: state_desc_q(pcnst)
      type(var_desc_t)     :: state_desc_pint
      type(var_desc_t)     :: state_desc_pintdry
      type(var_desc_t)     :: state_desc_lnpint
      type(var_desc_t)     :: state_desc_lnpintdry
      type(var_desc_t)     :: state_desc_zi

      type(var_desc_t)     :: tend_desc_dtdt
      type(var_desc_t)     :: tend_desc_dudt
      type(var_desc_t)     :: tend_desc_dvdt
      type(var_desc_t)     :: tend_desc_flx_net
      type(var_desc_t)     :: tend_desc_te_tnd
      type(var_desc_t)     :: tend_desc_tw_tnd

      type(var_desc_t)     :: desc_cflx(pcnst)
      type(var_desc_t)     :: desc_lhf
      type(var_desc_t)     :: desc_shf
      !-------------------------------------------------------------------------
      ! Initialize stuff
      !-------------------------------------------------------------------------
      if (mode==1) then
         fspec = 'mli'
         add_pbuf        = .true.
         add_phys_state  = .true.
         add_cam_in      = .true.
      end if
      if (mode==2) then
         fspec = 'mlo'
         add_pbuf       = .true.
         add_phys_state = .true.
         add_phys_tend  = .true.
         add_cam_out    = .true.
      end if

      do i=begchunk,endchunk
         ncol(i) = phys_state(i)%ncol
      end do

      grid_id = cam_grid_id('physgrid')
      call cam_grid_dimensions(grid_id, gdims(1:2), nhdims)
      if (nhdims==1) then
         call cam_grid_get_decomp(grid_id, (/pcols,endchunk-begchunk+1/), (/gdims(1)/),      pio_double, iodesc2d)
         call cam_grid_get_decomp(grid_id, (/pcols,endchunk-begchunk+1/), (/gdims(1),pver/), pio_double, iodesc3d)
      end if
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Initialize file and define variables
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      call cam_pio_createfile(file, trim(get_ml_filename(fspec,yr,mn,dy,sec)))      
      call cam_grid_write_attr(file, grid_id, header_info)

      call cam_pio_def_dim(file, 'lev',   pver,  pver_id,  existOK=.true.)
      call cam_pio_def_dim(file, 'ilev',  pverp, pverp_id, existOK=.true.)
      call cam_pio_def_dim(file, 'pcnst', pcnst, pcnst_id, existOK=.true.)

      call timemgr_init_restart(File)

      if (add_pbuf) call pbuf_init_restart_alt(file, pbuf2d)

      ! data for prognostic chemistry (does nothing for "chem none" option)
      call chem_init_restart(file)
      
      !-------------------------------------------------------------------------
      ! define physics state variables
      if (add_phys_state) then
         do m=1,pcnst
            write(num,'(i4.4)') m
            ierr = pio_def_var(file, 'state_q'//num, pio_double, dimids_3D1, state_desc_q(m))
         end do
         ierr = pio_def_var(file, 'state_t',         pio_double, dimids_3D1, state_desc_t)
         ierr = pio_def_var(file, 'state_u',         pio_double, dimids_3D1, state_desc_u)
         ierr = pio_def_var(file, 'state_v',         pio_double, dimids_3D1, state_desc_v)
         ierr = pio_def_var(file, 'state_s',         pio_double, dimids_3D1, state_desc_s)
         ierr = pio_def_var(file, 'state_omega',     pio_double, dimids_3D1, state_desc_omega)
         ierr = pio_def_var(file, 'state_pmid',      pio_double, dimids_3D1, state_desc_pmid)
         ierr = pio_def_var(file, 'state_pmiddry',   pio_double, dimids_3D1, state_desc_pmiddry)
         ierr = pio_def_var(file, 'state_pdel',      pio_double, dimids_3D1, state_desc_pdel)
         ierr = pio_def_var(file, 'state_pdeldry',   pio_double, dimids_3D1, state_desc_pdeldry)
         ierr = pio_def_var(file, 'state_rpdel',     pio_double, dimids_3D1, state_desc_rpdel)
         ierr = pio_def_var(file, 'state_rpdeldry',  pio_double, dimids_3D1, state_desc_rpdeldry)
         ierr = pio_def_var(file, 'state_lnpmid',    pio_double, dimids_3D1, state_desc_lnpmid)
         ierr = pio_def_var(file, 'state_lnpmiddry', pio_double, dimids_3D1, state_desc_lnpmiddry)
         ierr = pio_def_var(file, 'state_exner',     pio_double, dimids_3D1, state_desc_exner)
         ierr = pio_def_var(file, 'state_zm',        pio_double, dimids_3D1, state_desc_zm)
         ierr = pio_def_var(file, 'state_pint',      pio_double, dimids_3D2, state_desc_pint)
         ierr = pio_def_var(file, 'state_pintdry',   pio_double, dimids_3D2, state_desc_pintdry)
         ierr = pio_def_var(file, 'state_lnpint',    pio_double, dimids_3D2, state_desc_lnpint)
         ierr = pio_def_var(file, 'state_lnpintdry', pio_double, dimids_3D2, state_desc_lnpintdry)
         ierr = pio_def_var(file, 'state_zi',        pio_double, dimids_3D2, state_desc_zi)
      end if

      !-------------------------------------------------------------------------
      ! define physics tendency variables
      if (add_phys_tend) then
         ierr = pio_def_var(file, 'tend_dtdt',   pio_double, dimids_3D1, tend_desc_dtdt )
         ierr = pio_def_var(file, 'tend_dudt',   pio_double, dimids_3D1, tend_desc_dudt )
         ierr = pio_def_var(file, 'tend_dvdt',   pio_double, dimids_3D1, tend_desc_dvdt)
         ierr = pio_def_var(file, 'tend_flx_net',pio_double, dimids_hrz, tend_desc_flx_net)
         ierr = pio_def_var(file, 'tend_te_tnd', pio_double, dimids_hrz, tend_desc_te_tnd) ! cumulative boundary flux of total energy
         ierr = pio_def_var(file, 'tend_tw_tnd', pio_double, dimids_hrz, tend_desc_tw_tnd) ! cumulative boundary flux of total water
      end if

      !-------------------------------------------------------------------------
      ! define cam_out variables
      if (add_cam_out) then
         ierr = pio_def_var(file, 'cam_out_FLWDS', pio_double, dimids_hrz, desc_flwds )
         ierr = pio_def_var(file, 'cam_out_SOLS',  pio_double, dimids_hrz, desc_sols )
         ierr = pio_def_var(file, 'cam_out_SOLL',  pio_double, dimids_hrz, desc_soll )
         ierr = pio_def_var(file, 'cam_out_SOLSD', pio_double, dimids_hrz, desc_solsd )
         ierr = pio_def_var(file, 'cam_out_SOLLD', pio_double, dimids_hrz, desc_solld )
      end if

      !-------------------------------------------------------------------------
      ! define cam_in variables
      if (add_cam_in) then
         ierr = pio_def_var(file, 'cam_in_SHF',  pio_double, dimids_hrz, desc_shf)
         ierr = pio_def_var(file, 'cam_in_LHF',  pio_double, dimids_hrz, desc_lhf)
         do m=1,pcnst
            write(num,'(i4.4)') m
            ierr = pio_def_var(file, 'cam_in_CFLX'//num,  pio_double, dimids_hrz, desc_cflx(m))
         end do
      end if

      !-------------------------------------------------------------------------
      ! End variable definitions
      ierr = pio_enddef(file)

      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! write data to file
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      
      ! Set missing portions of tmp variables
      ! (only do this once since all chunked vars have the same shape)
      do i=begchunk,endchunk
         if (ncol(i) < pcols) then
            tmp2D(ncol(i)+1:,i)   = fillvalue
            tmp3D(ncol(i)+1:,:,i) = fillvalue
         end if
      end do

      call cam_grid_write_var(file, grid_id)

      call timemgr_write_restart(file)

      if (add_pbuf) call pbuf_write_restart_alt(file, pbuf2d)

      ! data for prognostic chemistry (does nothing for "chem none" option)
      call chem_write_restart(file)

      !-------------------------------------------------------------------------
      ! write physics state variables
      if (add_phys_state) then

         do m=1,pcnst
            do i=begchunk,endchunk
               tmp3D(:ncol(i),:,i) = phys_state(i)%q(:ncol(i),:,m) 
            end do
            call pio_write_darray(file, state_desc_q(m), iodesc3d, tmp3D, ierr)
         end do

         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%t(:ncol(i),:)
         end do
         call pio_write_darray(file, state_desc_t, iodesc3d, tmp3D, ierr)

         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%u(:ncol(i),:)
         end do
         call pio_write_darray(file, state_desc_u, iodesc3d, tmp3D, ierr)


         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%v(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_v, iodesc3d, tmp3D, ierr)

         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%s(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_s, iodesc3d, tmp3D, ierr)

         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%omega(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_omega, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%pmid(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_pmid, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%pmiddry(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_pmiddry, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%pdel(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_pdel, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%pdeldry(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_pdeldry, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%rpdel(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_rpdel, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%rpdeldry(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_rpdeldry, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%lnpmid(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_lnpmid, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%lnpmiddry(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_lnpmiddry, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%exner(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_exner, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%zm(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_zm, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%pint(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_pint, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%pintdry(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_pintdry, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%lnpint(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_lnpint, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%lnpintdry(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_lnpintdry, iodesc3d, tmp3D, ierr)
         
         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%zi(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_zi, iodesc3d, tmp3D, ierr)

      end if

      !-------------------------------------------------------------------------
      ! write physics state variables
      if (add_phys_tend) then

         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_tend(i)%dtdt(:ncol(i),:)
         end do
         call pio_write_darray(file, tend_desc_dtdt, iodesc3d, tmp3D, ierr)

         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_tend(i)%dudt(:ncol(i),:)
         end do
         call pio_write_darray(file, tend_desc_dudt, iodesc3d, tmp3D, ierr)

         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_tend(i)%dvdt(:ncol(i),:)
         end do 
         call pio_write_darray(file, tend_desc_dvdt, iodesc3d, tmp3D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i),i) = phys_tend(i)%flx_net(:ncol(i))
         end do
         call pio_write_darray(file, tend_desc_flx_net, iodesc2d, tmp2D, ierr)
         do i=begchunk,endchunk
            tmp2D(:ncol(i),i) = phys_tend(i)%te_tnd(:ncol(i))
         end do
         call pio_write_darray(file, tend_desc_te_tnd, iodesc2d, tmp2D, ierr)
         do i=begchunk,endchunk
            tmp2D(:ncol(i),i) = phys_tend(i)%tw_tnd(:ncol(i))
         end do
         call pio_write_darray(file, tend_desc_tw_tnd, iodesc2d, tmp2D, ierr)

      end if

      !-------------------------------------------------------------------------
      ! Write cam_in components

      if (add_cam_in) then

         do m=1,pcnst
            do i=begchunk,endchunk
               tmp2D(:ncol(i), i) = cam_in(i)%cflx(:ncol(i), m)
            end do
            call pio_write_darray(file, desc_cflx(m), iodesc2d, tmp2D, ierr)
         end do

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_in(i)%shf(:ncol(i))
         end do
         call pio_write_darray(file, desc_shf, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_in(i)%lhf(:ncol(i))
         end do
         call pio_write_darray(file, desc_lhf, iodesc2d, tmp2D, ierr)

      end if

      !-------------------------------------------------------------------------
      ! Write cam_out components

      if (add_cam_out) then

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_out(i)%flwds(:ncol(i))
         end do
         call pio_write_darray(file, desc_flwds, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_out(i)%sols(:ncol(i))
         end do
         call pio_write_darray(file, desc_sols, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_out(i)%soll(:ncol(i))
         end do
         call pio_write_darray(file, desc_soll, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_out(i)%solsd(:ncol(i))
         end do
         call pio_write_darray(file, desc_solsd, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_out(i)%solld(:ncol(i))
         end do
         call pio_write_darray(file, desc_solld, iodesc2d, tmp2D, ierr)

      end if

      !-------------------------------------------------------------------------
      ! close the file
      call pio_closefile(file)
      ! call cam_pio_closefile(file)
      
   end subroutine write_ml_training
   !------------------------------------------------------------------------------------------------
end module ml_training
