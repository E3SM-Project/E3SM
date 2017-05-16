module dyn_comp
!----------------------------------------------------------------------- 
! 
! Eulerian dycore interface module
!
!-----------------------------------------------------------------------

use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
use spmd_utils,   only: masterproc
use constituents, only: pcnst, cnst_name, cnst_longname
use constituents, only: sflxnam, tendnam, fixcnam, tottnam, hadvnam, vadvnam, cnst_get_ind
use pmgrid,       only: plev, plevp, dyndecomp_set
use hycoef,       only: hycoef_init
use cam_history,  only: addfld, add_default, horiz_only
use phys_control, only: phys_getopts
use eul_control_mod, only: dyn_eul_readnl, eul_nsplit
use cam_logfile,  only: iulog
use pio,          only: file_desc_t

#if (defined SPMD)
use spmd_dyn,     only: spmd_readnl, spmdinit_dyn
#endif

implicit none
private

public :: dyn_init, dyn_import_t, dyn_export_t

! these structures are not used in this dycore, but are included
! for source code compatibility.  
type dyn_import_t
   integer :: placeholder
end type dyn_import_t

type dyn_export_t
   integer :: placeholder
end type dyn_export_t

! Frontogenesis indices
integer, public :: frontgf_idx = -1
integer, public :: frontga_idx = -1

! Index into physics buffer for zonal mean zonal wind.
integer, public :: uzm_idx = -1

!#######################################################################
CONTAINS
!#######################################################################

subroutine dyn_init(file, nlfilename)
   use phys_control, only : use_gw_front
   use dyn_grid,     only: define_cam_grids, initgrid
   use scamMod,      only: single_column
   use physics_buffer,  only : pbuf_add_field, dtype_r8
   use ppgrid,          only : pcols, pver

   ! ARGUMENTS:
   type(file_desc_t), intent(in) :: file       ! PIO file handle for initial or restart file
   character(len=*),  intent(in) :: nlfilename

   ! Local workspace
   integer m                     ! Index
   integer :: ixcldice, ixcldliq ! constituent indices for cloud liquid and ice water.
   logical :: history_amwg       ! output for AMWG diagnostics
   logical :: history_budget     ! output tendencies and state variables for CAM4
                                 ! temperature, water vapor, cloud ice and cloud
                                 ! liquid budgets.
   integer :: history_budget_histfile_num  ! output history file number for budget fields
   !----------------------------------------------------------------------------

   dyndecomp_set=.true.

   call dyn_eul_readnl(nlfilename)
   if (masterproc) write(iulog,*) 'EUL subcycling - eul_nsplit = ', eul_nsplit

#if (defined SPMD)
   call spmd_readnl(nlfilename)
   call spmdinit_dyn()
#endif 

  !----------------------------------------------------------------------

  if (use_gw_front) then
     call pbuf_add_field("FRONTGF", "global", dtype_r8, (/pcols,pver/), &
          frontgf_idx)
     call pbuf_add_field("FRONTGA", "global", dtype_r8, (/pcols,pver/), &
          frontga_idx)
  end if

   ! Initialize hybrid coordinate arrays
   call hycoef_init(file)

   ! Run initgrid (the old initcom) which sets up coordinates and weights
   call initgrid()

   ! Define the CAM grids (must be before addfld calls)
   call define_cam_grids()

   call addfld ('ETADOT',(/ 'ilev' /),'A', '1/s','Vertical (eta) velocity',             gridname='gauss_grid')
   call addfld ('U&IC',  (/ 'lev' /), 'I', 'm/s','Zonal wind',                          gridname='gauss_grid' )
   call addfld ('V&IC',  (/ 'lev' /), 'I', 'm/s','Meridional wind',                     gridname='gauss_grid' )
   call add_default ('U&IC',0, 'I')
   call add_default ('V&IC',0, 'I')

   call addfld ('PS&IC',horiz_only,'I',    'Pa','Surface pressure',                     gridname='gauss_grid' )
   call addfld ('T&IC',(/ 'lev' /),'I', 'K','Temperature',                              gridname='gauss_grid' )
   call add_default ('PS&IC',0, 'I')
   call add_default ('T&IC',0, 'I')

   do m = 1, pcnst
      call addfld (trim(cnst_name(m))//'&IC',(/ 'lev' /),'I', 'kg/kg',cnst_longname(m), gridname='gauss_grid' )
      call add_default(trim(cnst_name(m))//'&IC',0, 'I')
      call addfld (hadvnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' horizontal advection tendency',  &
           gridname='gauss_grid')
      call addfld (vadvnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' vertical advection tendency',    &
           gridname='gauss_grid')
      call addfld (tendnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' total tendency',                 &
           gridname='gauss_grid')
      call addfld (tottnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' horz + vert + fixer tendency',   &
           gridname='gauss_grid')
      call addfld (fixcnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' tendency due to slt fixer',      &
           gridname='gauss_grid')
   end do

   call addfld ('DUH     ',(/ 'lev' /),'A', 'K/s     ','U horizontal diffusive heating',              gridname='gauss_grid')
   call addfld ('DVH     ',(/ 'lev' /),'A', 'K/s     ','V horizontal diffusive heating',              gridname='gauss_grid')
   call addfld ('DTH     ',(/ 'lev' /),'A', 'K/s     ','T horizontal diffusive heating',              gridname='gauss_grid')

   call addfld ('ENGYCORR',(/ 'lev' /),'A', 'W/m2    ','Energy correction for over-all conservation', gridname='gauss_grid')
   call addfld ('TFIX    ',horiz_only ,'A', 'K/s     ','T fixer (T equivalent of Energy correction)', gridname='gauss_grid')

   call addfld ('FU      ',(/ 'lev' /),'A', 'm/s2    ','Zonal wind forcing term',                     gridname='gauss_grid')
   call addfld ('FV      ',(/ 'lev' /),'A', 'm/s2    ','Meridional wind forcing term',                gridname='gauss_grid')
   call addfld ('UTEND   ',(/ 'lev' /),'A', 'm/s2    ','U tendency',                                  gridname='gauss_grid')
   call addfld ('VTEND   ',(/ 'lev' /),'A', 'm/s2    ','V tendency',                                  gridname='gauss_grid')
   call addfld ('TTEND   ',(/ 'lev' /),'A', 'K/s     ','T tendency',                                  gridname='gauss_grid')
   call addfld ('LPSTEN  ',horiz_only ,'A', 'Pa/s    ','Surface pressure tendency',                   gridname='gauss_grid')
   call addfld ('VAT     ',(/ 'lev' /),'A', 'K/s     ','Vertical advective tendency of T',            gridname='gauss_grid')
   call addfld ('KTOOP   ',(/ 'lev' /),'A', 'K/s     ','(Kappa*T)*(omega/P)',                         gridname='gauss_grid')

   call phys_getopts(history_amwg_out=history_amwg, &
        history_budget_out = history_budget, &
        history_budget_histfile_num_out = history_budget_histfile_num)

   if (history_amwg) then
      call add_default ('DTH     ', 1, ' ')
   end if

   if ( history_budget ) then
      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('CLDICE', ixcldice)
      ! The following variables are not defined for single column
      if (.not. single_column) then 
         call add_default(hadvnam(       1), history_budget_histfile_num, ' ')
         call add_default(hadvnam(ixcldliq), history_budget_histfile_num, ' ')
         call add_default(hadvnam(ixcldice), history_budget_histfile_num, ' ')
         call add_default(vadvnam(       1), history_budget_histfile_num, ' ')
         call add_default(vadvnam(ixcldliq), history_budget_histfile_num, ' ')
         call add_default(vadvnam(ixcldice), history_budget_histfile_num, ' ')
      end if
      call add_default(fixcnam(       1), history_budget_histfile_num, ' ')
      call add_default(fixcnam(ixcldliq), history_budget_histfile_num, ' ')
      call add_default(fixcnam(ixcldice), history_budget_histfile_num, ' ')
      call add_default(tottnam(       1), history_budget_histfile_num, ' ')
      call add_default(tottnam(ixcldliq), history_budget_histfile_num, ' ')
      call add_default(tottnam(ixcldice), history_budget_histfile_num, ' ')
      call add_default(tendnam(       1), history_budget_histfile_num, ' ')
      call add_default(tendnam(ixcldliq), history_budget_histfile_num, ' ')
      call add_default(tendnam(ixcldice), history_budget_histfile_num, ' ')
      call add_default('TTEND',           history_budget_histfile_num, ' ')
      call add_default('TFIX',            history_budget_histfile_num, ' ')
      call add_default('KTOOP',           history_budget_histfile_num, ' ')
      call add_default('VAT',             history_budget_histfile_num, ' ')
      call add_default('DTH',             history_budget_histfile_num, ' ')
   end if

end subroutine dyn_init

end module dyn_comp
