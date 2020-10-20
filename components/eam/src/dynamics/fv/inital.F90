module inital

! Dynamics initialization

implicit none
private

public cam_initial

!=========================================================================
contains
!=========================================================================

subroutine cam_initial( dyn_in, dyn_out, NLFileName )

   use dyn_comp,             only : dyn_import_t, dyn_export_t, &
                                    dyn_init, frontgf_idx, frontga_idx, &
                                    uzm_idx
   use phys_grid,            only : phys_grid_init
   use physpkg,              only : phys_register
   use phys_control,         only : use_gw_front
   use physics_buffer,       only : pbuf_add_field, dtype_r8
   use ppgrid,               only : pcols, pver
   use qbo,                  only : qbo_use_forcing
   use chem_surfvals,        only : chem_surfvals_init
   use cam_initfiles,        only : initial_file_get_id
   use startup_initialconds, only : initial_conds
   use ref_pres,             only : ref_pres_init
   use dynamics_vars,        only : T_FVDYCORE_STATE   
   use dyn_internal_state,   only : get_dyn_state

   ! Arguments
   type(dyn_import_t), intent(out) :: dyn_in
   type(dyn_export_t), intent(out) :: dyn_out
   character(len=*),   intent(in)  :: NLFileName

   ! LOCAL VARIABLES:
   type (T_FVDYCORE_STATE), pointer :: dyn_state
   !-----------------------------------------------------------------------

   dyn_state => get_dyn_state()

   ! Initialize dynamics grid
   ! This has to be before dyn_init because it defines the grid coords (needed
   !      by define_cam_grids to define coords and grids).
   call initcom

   call dyn_init(initial_file_get_id(), dyn_state, dyn_in, dyn_out, NLFileName )

   ! Define physics data structures
   call phys_grid_init()

   ! Initialize index values for advected and non-advected tracers
   call phys_register()

   ! Frontogenesis indices
   if (use_gw_front) then
      call pbuf_add_field("FRONTGF", "global", dtype_r8, (/pcols,pver/), &
           frontgf_idx)
      call pbuf_add_field("FRONTGA", "global", dtype_r8, (/pcols,pver/), &
           frontga_idx)
   end if

   if (qbo_use_forcing) then
      call pbuf_add_field("UZM", "global", dtype_r8, (/pcols,pver/), &
           uzm_idx)
   end if

   ! Initialize ghg surface values before default initial distributions
   ! are set in inidat.
   call chem_surfvals_init()

   ! Initialize physics grid reference pressures (needed by initialize_radbuffer)
   call ref_pres_init()

   call initial_conds( dyn_in )

end subroutine cam_initial

!=========================================================================

end module inital
