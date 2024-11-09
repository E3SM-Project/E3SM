module startup_initialconds
!----------------------------------------------------------------------- 
! 
! Wrapper for calls to initialize buffers and read initial/topo files
! 
!-----------------------------------------------------------------------

use pio,          only: file_desc_t

implicit none
private
save

public :: initial_conds ! Read in initial conditions (dycore dependent)
!added for orographic drag
public topo_OD_file_get_id
public setup_initial_OD
public close_initial_file_OD
type(file_desc_t), pointer :: ncid_topo_OD

!======================================================================= 
contains
!======================================================================= 

function topo_OD_file_get_id()
        type(file_desc_t), pointer :: topo_OD_file_get_id
        topo_OD_file_get_id => ncid_topo_OD
end function topo_OD_file_get_id

subroutine initial_conds(dyn_in)

   ! This routine does some initializing of buffers that should move to a
   ! non-dycore dependent location.  It also calls the dycore dependent
   ! routine read_inidat.

   use comsrf,        only: initialize_comsrf
#if (defined E3SM_SCM_REPLAY )
   use history_defaults, only: initialize_iop_history
#endif
   
   use pio,           only: file_desc_t
   use cam_initfiles, only: initial_file_get_id, topo_file_get_id
   use inidat,        only: read_inidat
   use dyn_comp,      only: dyn_import_t
   use cam_pio_utils, only: clean_iodesc_list
   use perf_mod,      only: t_startf, t_stopf
   
   type(dyn_import_t),  intent(inout) :: dyn_in

   ! Local variables
   type(file_desc_t), pointer :: fh_ini, fh_topo
   !-----------------------------------------------------------------------

   ! Initialize buffer, comsrf, and radbuffer variables 
   ! (which must occur after the call to phys_grid_init)
   call initialize_comsrf

#if (defined E3SM_SCM_REPLAY )
   call initialize_iop_history
#endif

   ! Read in initial data
   call t_startf('read_inidat')
   fh_ini  => initial_file_get_id()
   fh_topo => topo_file_get_id()
   call read_inidat(fh_ini, fh_topo, dyn_in)
   call t_stopf('read_inidat')

   call t_startf('clean_iodesc_list')
   call clean_iodesc_list()
   call t_stopf('clean_iodesc_list')

end subroutine initial_conds

!======================================================================= 

subroutine setup_initial_OD()
   use filenames,        only: bnd_topo
   use ioFileMod,        only: getfil
   use cam_pio_utils,    only: cam_pio_openfile
   use pio,              only: pio_nowrite
!
! Input arguments
!
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!-----------------------------------------------------------------------
   character(len=256) :: bnd_topo_loc   ! filepath of topo file on local disk
      allocate(ncid_topo_OD)
      call getfil(bnd_topo, bnd_topo_loc)
      call cam_pio_openfile(ncid_topo_OD, bnd_topo_loc, PIO_NOWRITE)
end subroutine setup_initial_OD

subroutine close_initial_file_OD
  use pio,          only: pio_closefile
        call pio_closefile(ncid_topo_OD)
        deallocate(ncid_topo_OD)
        nullify(ncid_topo_OD)
end subroutine close_initial_file_OD
!======================================================================= 





end module startup_initialconds
