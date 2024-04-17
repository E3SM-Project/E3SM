module inital

! Dynamics initialization

implicit none
private

public :: cam_initial

!=========================================================================
contains
!=========================================================================

subroutine cam_initial(dyn_in, dyn_out, NLFileName)

   use dyn_comp,             only: dyn_init1, dyn_init2, dyn_import_t, &
                                   dyn_export_t, frontgf_idx, frontga_idx
   use phys_grid,            only: phys_grid_init
   use physpkg,              only: phys_register
   use phys_control,         only: use_gw_front
   use physics_buffer,       only: pbuf_add_field, dtype_r8
   use ppgrid,               only: pcols, pver
   use chem_surfvals,        only: chem_surfvals_init
   use cam_initfiles,        only: initial_file_get_id
   use startup_initialconds, only: initial_conds
   use cam_logfile,          only: iulog
   use perf_mod,             only: t_startf, t_stopf
   use phys_grid_ctem,       only: phys_grid_ctem_reg
   
   ! modules from SE
   use parallel_mod, only : par

   type(dyn_import_t), intent(out) :: dyn_in
   type(dyn_export_t), intent(out) :: dyn_out
   character(len=*),   intent(in)  :: NLFileName
   !----------------------------------------------------------------------

   call t_startf('dyn_init1')
   call dyn_init1(initial_file_get_id(), NLFileName, dyn_in, dyn_out)
   call t_stopf('dyn_init1')

   ! Define physics data structures
   if(par%masterproc  ) write(iulog,*) 'Running phys_grid_init()'
   call phys_grid_init()

   ! Register zonal average grid for phys TEM diagnostics
   call phys_grid_ctem_reg()

   ! Initialize index values for advected and non-advected tracers
   call phys_register()

   ! Frontogenesis indices
   if (use_gw_front) then
      call pbuf_add_field("FRONTGF", "global", dtype_r8, (/pcols,pver/), &
           frontgf_idx)
      call pbuf_add_field("FRONTGA", "global", dtype_r8, (/pcols,pver/), &
           frontga_idx)
   end if

   ! Initialize ghg surface values before default initial distributions
   ! are set in inidat.
   call chem_surfvals_init()

   if(par%masterproc  ) write(iulog,*) 'Reading initial data'

   call t_startf('initial_conds')
   call initial_conds(dyn_in)
   call t_stopf('initial_conds')

   if(par%masterproc  ) write(iulog,*) 'Done Reading initial data'

   call t_startf('dyn_init2')
   call dyn_init2(dyn_in)
   call t_stopf('dyn_init2')

end subroutine cam_initial

end module inital
