module physics_update_mod 
  !============================================================================
  !An interface for physics update call:
  !For perturbation growth test, we need to output certain variables after every
  !physics_update call (when the state variables are updated). Placing outfld
  !directly in the existing physics_update call resulted in hard to resolve
  !circular dependencies. The alternative was to create this interface which calls 
  !physics_update and outfld subroutines one after another.

  !The circulare dependencies, if we add outfld call to physics_types module, were:
  ! [NOTE: '<-' can be read as "depends upon" as the example below can be read as 
  !"physics_type 'depends upon' cam_history" (due to the outfld call in physics_type )"]

  !1. physics_type<-cam_history<-subcol_utils<-physics_types
  !2. physics_types<-cam_history<-chem_surfvals<-mo_flbc<-phys_gmean<-physics_type

  !This module helps breaking these circular dependencies but it requires changes 
  !to other parts of the code where we need to replace physics_update call with
  !physics_update_intr call and change the associated use statement. 
  !============================================================================


  implicit none
  private
  public  :: physics_update_intr_init, physics_update_intr

  logical :: pergro_test_active
  
contains 

  subroutine physics_update_intr_init()
    use phys_control, only: phys_getopts

    call phys_getopts(pergro_test_active_out = pergro_test_active)

  end subroutine physics_update_intr_init


  subroutine physics_update_intr (state, ptend, dt, tend)

    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols
    use physics_types, only: physics_update, physics_ptend, physics_state, physics_tend
    use cam_history,   only: outfld, fieldname_len  
    
    !Arguments
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
    type(physics_state), intent(inout)  :: state   ! Physics state variables
    real(r8),            intent(in)     :: dt      ! time step
    
    !optional arguments
    type(physics_tend ), intent(inout), optional  :: tend  ! Physics tendencies over timestep
    
    !local vars
    character(len=fieldname_len)   :: vname
    
    call physics_update (state, ptend, dt, tend)
    if (pergro_test_active) then
          
       !BSINGH - output fields
       vname = 'T_'//trim(adjustl(ptend%name))
       call outfld( trim(adjustl(vname)), state%t, pcols, state%lchnk )
       
       vname = 'S_'//trim(adjustl(ptend%name))
       call outfld( trim(adjustl(vname)), state%s, pcols, state%lchnk )
       
       !vname = 'St_'//trim(adjustl(ptend%name))
       !call outfld( trim(adjustl(vname)), ptend%s, pcols, state%lchnk )
       
       vname = 'QV_'//trim(adjustl(ptend%name))
       call outfld( trim(adjustl(vname)), state%q(:,:,1), pcols, state%lchnk )
       
       !vname = 'QVt_'//trim(adjustl(ptend%name))
       !call outfld( trim(adjustl(vname)), ptend%q(:,:,1), pcols, state%lchnk )
    endif
  end subroutine physics_update_intr
      
end module physics_update_mod
