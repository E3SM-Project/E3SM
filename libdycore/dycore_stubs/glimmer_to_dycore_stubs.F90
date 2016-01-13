! The glimmer_to_dycore stubs module contains stubs for the Fortran side of the Glimmer-DyCore
! interface.  It uses the routines in dycore_to_glim_extern.cpp to create one
! or more instances of a dynamic core ice sheet model.  The dycore_model_index is
! the only parameter needed by glimmer_to_dycore subroutines to interact with a 
! specific instance of a dynamic core model. DMR--5/24/10   

module glimmer_to_dycore
  !*FD glimmer_to_dycore contains Fortran routines to couple Glimmer to a
  ! dynamic core model.
  use glide_types
  use simple_forcing

  contains

  subroutine gtd_init_dycore_interface(model,dycore_type,dycore_model_index)
    type(glide_global_type) :: model
    integer*4 dycore_model_index, error_code
    integer*4 dycore_type  ! 0=BISICLES, 1=Ymir

!    call dycore_init_registry()
!    call dycore_init_model(dycore_type,dycore_model_index,error_code)
!    call gtd_set_geometry_vars(model,dycore_model_index)

!    print *,"In init_dycore_interface, dycore_type = ",dycore_type
!    print *,"In init_dycore_interface, dycore1 = ",dycore_model_index
  end subroutine gtd_init_dycore_interface

  subroutine gtd_run_dycore(dycore_model_index)
    integer*4 dycore_model_index

!    call dycore_run_model(dycore_model_index)
  end subroutine gtd_run_dycore

  subroutine gtd_set_geometry_vars(model,dycore_model_index)
    type(glide_global_type) :: model
    integer*4 dycore_model_index

  end subroutine gtd_set_geometry_vars 


  subroutine gtd_set_velocity_vars(model,dycore_model_index)
    type(glide_global_type) :: model
    integer*4 dycore_model_index
 
  end subroutine gtd_set_velocity_vars  



end module glimmer_to_dycore
