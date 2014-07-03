module inital

! Dynamics initialization

implicit none
private

public :: cam_initial

!=========================================================================
contains
!=========================================================================

subroutine cam_initial(dyn_in, dyn_out, nlfilename)

   use dyn_comp,             only: dyn_import_t, dyn_export_t, dyn_init
   use prognostics,          only: initialize_prognostics
   use phys_grid,            only: phys_grid_init
   use chem_surfvals,        only: chem_surfvals_init
   use scanslt,              only: slt_alloc
   use cam_initfiles,        only: initial_file_get_id
   use startup_initialconds, only: initial_conds
   use ref_pres,             only: ref_pres_init
#if (defined SPMD)
   use spmd_dyn,             only: spmdbuf
#endif

   ! Arguments are not used in this dycore, included for compatibility
   type(dyn_import_t) :: dyn_in
   type(dyn_export_t) :: dyn_out
   character(len=*), intent(in) :: nlfilename
   !-----------------------------------------------------------------------

   ! Initialize ghg surface values before default initial distributions
   ! are set in inidat.
   call chem_surfvals_init()

   call dyn_init(initial_file_get_id(), nlfilename)

   ! Initialize prognostics variables
   call initialize_prognostics
   call slt_alloc()

   call initcom

   ! Define physics data structures
   call phys_grid_init

#if (defined SPMD)
   ! Allocate communication buffers for
   ! collective communications in realloc
   ! routines and in dp_coupling
   call spmdbuf ()
#endif

   ! Initialize physics grid reference pressures (needed by initialize_radbuffer)
   call ref_pres_init()

   call initial_conds(dyn_in)

end subroutine cam_initial

!=========================================================================

end module inital
