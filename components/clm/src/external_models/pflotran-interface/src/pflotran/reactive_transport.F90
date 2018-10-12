module Reactive_Transport_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Transport_module
  use Reaction_module

  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  public :: RTTimeCut, &
            RTSetup, &
            RTMaxChange, &
            RTUpdateEquilibriumState, &
            RTUpdateKineticState, &
            RTUpdateMassBalance, &
            RTResidual, &
            RTJacobian, &
            RTInitializeTimestep, &
            RTUpdateAuxVars, &
            RTComputeMassBalance, &
            RTDestroy, &
            RTUpdateTransportCoefs, &
            RTUpdateActivityCoefficients, &
            RTUpdateRHSCoefs, &
            RTCalculateRHS_t0, &
            RTCalculateRHS_t1, &
            RTCalculateTransportMatrix, &
            RTReact, &
            RTExplicitAdvection, &
            RTClearActivityCoefficients
  
contains

! ************************************************************************** !

subroutine RTTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/15/08
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Global_module
  use Secondary_Continuum_module, only : SecondaryRTTimeCut
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  
  PetscErrorCode :: ierr

  field => realization%field
  option => realization%option
 
  ! copy previous solution back to current solution
  call VecCopy(field%tran_yy,field%tran_xx,ierr);CHKERRQ(ierr)
  
  ! set densities and saturations to t+dt
  if (realization%option%nflowdof > 0) then
    call GlobalWeightAuxVars(realization, &
                             realization%option%transport%tran_weight_t1)
  endif

  if (option%use_mc) then
    call SecondaryRTTimeCut(realization)
  endif
 
end subroutine RTTimeCut

! ************************************************************************** !

subroutine RTSetup(realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Region_module
  use Coupler_module
  use Condition_module
  use Connection_module
  use Transport_Constraint_module
  use Fluid_module
  use Material_module
  use Material_Aux_class

  !geh: please leave the "only" clauses for Secondary_Continuum_XXX as this
  !      resolves a bug in the Intel Visual Fortran compiler.
  use Secondary_Continuum_Aux_module, only : sec_transport_type, &
                                             SecondaryAuxRTCreate
  use Secondary_Continuum_module, only : SecondaryRTAuxVarInit
  use Output_Aux_module
  use Generic_module
  use String_module
 
  implicit none

  type(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(output_variable_list_type), pointer :: list
  type(reaction_type), pointer :: reaction
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(fluid_property_type), pointer :: cur_fluid_property
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  type(coupler_type), pointer :: initial_condition
  type(tran_constraint_type), pointer :: sec_tran_constraint
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_property_type), pointer :: cur_material_property
  type(reactive_transport_param_type), pointer :: rt_parameter
  type(generic_parameter_type), pointer :: cur_generic_parameter

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: ghosted_id, iconn, sum_connection
  PetscInt :: iphase, local_id, i
  PetscBool :: error_found
  PetscInt :: flag(10)  
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  reaction => realization%reaction
  sec_tran_constraint => realization%sec_transport_constraint

  patch%aux%RT => RTAuxCreate(reaction%naqcomp,option%transport%nphase)
  rt_parameter => patch%aux%RT%rt_parameter
  ! rt_parameter %naqcomp and %nphase set in RTAuxCreate()
  rt_parameter%ncomp = reaction%ncomp
  rt_parameter%offset_aqueous = reaction%offset_aqueous
  rt_parameter%nimcomp = reaction%immobile%nimmobile
  rt_parameter%ngas = reaction%gas%nactive_gas
  rt_parameter%offset_immobile = reaction%offset_immobile
  if (reaction%gas%nactive_gas > 0 .and. option%transport%nphase == 1) then
    option%io_buffer = 'INCLUDE_GAS_PHASE must be included in SUBSURFACE_&
      &TRANSPORT block when simulating active gases.'
    call printErrMsg(option)
  endif

  if (reaction%ncollcomp > 0) then
    rt_parameter%ncoll = reaction%ncoll
    rt_parameter%offset_colloid  = reaction%offset_colloid 
    rt_parameter%ncollcomp = reaction%ncollcomp
    rt_parameter%offset_collcomp = reaction%offset_collcomp
    allocate(rt_parameter%pri_spec_to_coll_spec(reaction%naqcomp))
    rt_parameter%pri_spec_to_coll_spec = &
      reaction%pri_spec_to_coll_spec
    allocate(rt_parameter%coll_spec_to_pri_spec(reaction%ncollcomp))
    rt_parameter%coll_spec_to_pri_spec = &
      reaction%coll_spec_to_pri_spec
  endif
  if (reaction%immobile%nimmobile > 0) then
    rt_parameter%nimcomp = reaction%immobile%nimmobile
    rt_parameter%offset_immobile = reaction%offset_immobile
  endif
  ! loop over material properties and determine if any transverse 
  ! dispersivities are defined.
  cur_material_property => realization%material_properties                            
  do                                      
    if (.not.associated(cur_material_property)) exit
    if (maxval(cur_material_property%dispersivity(2:3)) > 0.d0) then
      rt_parameter%calculate_transverse_dispersion = PETSC_TRUE
      exit
    endif
    cur_material_property => cur_material_property%next
  enddo
  
  material_auxvars => patch%aux%Material%auxvars
  flag = 0
  !TODO(geh): change to looping over ghosted ids once the legacy code is 
  !           history and the communicator can be passed down.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)

    ! Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    
    if (material_auxvars(ghosted_id)%volume < 0.d0 .and. flag(1) == 0) then
      flag(1) = 1
      option%io_buffer = 'Non-initialized cell volume.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%porosity < 0.d0 .and. flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'Non-initialized porosity.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%tortuosity < 0.d0 .and. flag(3) == 0) then
      flag(3) = 1
      option%io_buffer = 'Non-initialized tortuosity.'
      call printMsg(option)
    endif
    if (reaction%neqkdrxn > 0) then
      if (material_auxvars(ghosted_id)%soil_particle_density < 0.d0 .and. &
          flag(4) == 0) then
        flag(4) = 1
        option%io_buffer = 'Non-initialized soil particle density.'
        call printMsg(option)
      endif
    endif
  enddo  
 
  if (maxval(flag) > 0) then
    option%io_buffer = &
      'Material property errors found in RTSetup (reactive transport).'
    call printErrMsg(option)
  endif  
  
!============== Create secondary continuum variables - SK 2/5/13 ===============

  
  if (option%use_mc) then
    patch%aux%SC_RT => SecondaryAuxRTCreate(option)
    initial_condition => patch%initial_condition_list%first
    allocate(rt_sec_transport_vars(grid%nlmax))  
    do local_id = 1, grid%nlmax
    ! Assuming the same secondary continuum type for all regions
      call SecondaryRTAuxVarInit(patch%material_property_array(1)%ptr, &
                                 rt_sec_transport_vars(local_id), &
                                 reaction,initial_condition, &
                                 sec_tran_constraint,option)
    enddo      
    patch%aux%SC_RT%sec_transport_vars => rt_sec_transport_vars      
  endif

!===============================================================================   

    
  ! allocate auxvar data structures for all grid cells
#ifdef COMPUTE_INTERNAL_MASS_FLUX
  option%iflag = 1 ! allocate mass_balance array
#else  
  option%iflag = 0 ! be sure not to allocate mass_balance array
#endif
  allocate(patch%aux%RT%auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call RTAuxVarInit(patch%aux%RT%auxvars(ghosted_id),reaction,option)
  enddo
  patch%aux%RT%num_aux = grid%ngmax
  
  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array 
    allocate(patch%aux%RT%auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call RTAuxVarInit(patch%aux%RT%auxvars_bc(iconn),reaction,option)
    enddo
  endif
  patch%aux%RT%num_aux_bc = sum_connection

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array 
    allocate(patch%aux%RT%auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call RTAuxVarInit(patch%aux%RT%auxvars_ss(iconn),reaction,option)
    enddo
  endif
  patch%aux%RT%num_aux_ss = sum_connection
  option%iflag = 0

  ! ensure that active gas species only have one aqueous counterpart
  do i = 1, reaction%gas%nactive_gas
    if (reaction%gas%acteqspecid(0,i) > 1) then
      write(word,*) reaction%gas%acteqspecid(0,i)
      option%io_buffer = 'Acdtive gas species "' // &
        trim(reaction%gas%active_names(i)) // '" is associated with ' // &
        trim(adjustl(word)) // ' aqueous species when it may only be &
        &associated with 1.'
    endif
  enddo

  ! initialize parameters
  cur_fluid_property => realization%fluid_properties
  do 
    if (.not.associated(cur_fluid_property)) exit
    iphase = cur_fluid_property%phase_id
    ! setting of phase diffusion coefficients must come before individual
    ! species below
    rt_parameter%diffusion_coefficient(:,iphase) = &
      cur_fluid_property%diffusion_coefficient
    rt_parameter%diffusion_activation_energy(:) = &
      cur_fluid_property%diffusion_activation_energy
    cur_fluid_property => cur_fluid_property%next
  enddo

  ! overwrite diffusion coefficients with species specific values
  ! aqueous diffusion
  iphase = option%liquid_phase
  cur_generic_parameter => reaction%aq_diffusion_coefficients
  do
    if (.not.associated(cur_generic_parameter)) exit
    i = GetPrimarySpeciesIDFromName(cur_generic_parameter%name, &
                                    reaction,PETSC_FALSE,option)
    if (Uninitialized(i)) then
      option%io_buffer = 'Species "' // trim(cur_generic_parameter%name) // &
        '" listed in aqueous diffusion coefficient list not found among &
        &aqueous species.'
      call printErrMsg(option)
    endif
    rt_parameter%diffusion_coefficient(i,iphase) = &
        cur_generic_parameter%rvalue
    cur_generic_parameter => cur_generic_parameter%next
  enddo
  if (associated(reaction%gas_diffusion_coefficients)) then
    if (rt_parameter%nphase <= 1) then
      option%io_buffer = 'GAS_DIFFUSION_COEFFICIENTS may not be set when &
        &INCLUDE_GAS_PHASE is not included in the SIMULATION,&
        &SUBSURFACE_TRANSPORT block.'
      call printErrMsg(option)
    endif
    ! gas diffusion
    iphase = option%gas_phase
    cur_generic_parameter => reaction%gas_diffusion_coefficients
    do
      if (.not.associated(cur_generic_parameter)) exit
      do i = 1, reaction%gas%nactive_gas
        if (StringCompare(cur_generic_parameter%name, &
                          reaction%gas%active_names(i))) then
          ! maps gas to aqueous primary species, for which there should only
          ! be one per gas
          rt_parameter%diffusion_coefficient(reaction%gas%acteqspecid(1,i), &
                                             iphase) = &
            cur_generic_parameter%rvalue
          exit
        endif
      enddo
      if (i > reaction%gas%nactive_gas) then
        option%io_buffer = 'Species "' // trim(cur_generic_parameter%name) // &
          '" listed in gas diffusion coefficient list has no aqueous &
          &counterpart'
        call printErrMsg(option)
      endif
      cur_generic_parameter => cur_generic_parameter%next
    enddo
  endif
 
  list => realization%output_option%output_snap_variable_list
  call RTSetPlotVariables(list,reaction,option, &
                          realization%output_option%tunit)
  if (.not.associated(realization%output_option%output_snap_variable_list, &
                 realization%output_option%output_obs_variable_list)) then
    list => realization%output_option%output_obs_variable_list
    call RTSetPlotVariables(list,reaction,option, &
                            realization%output_option%tunit)
  endif
  
end subroutine RTSetup

! ************************************************************************** !

subroutine RTComputeMassBalance(realization,mass_balance)
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/23/08
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module

  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%ntrandof, &
                            realization%option%transport%nphase)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(reaction_type), pointer :: reaction

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iphase
  PetscInt :: i, icomp, imnrl, ncomp, irate, irxn, naqcomp
  PetscInt :: nphase

  iphase = 1
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  reaction => realization%reaction

  rt_auxvars => patch%aux%RT%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  mass_balance = 0.d0
  naqcomp = reaction%naqcomp
  nphase = option%transport%nphase

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    do iphase = 1, nphase
    
!     mass_balance(:,iphase) = mass_balance(:,iphase) + &
      mass_balance(1:naqcomp,1) = &
        mass_balance(1:naqcomp,1) + &
        rt_auxvars(ghosted_id)%total(:,iphase) * &
        global_auxvars(ghosted_id)%sat(iphase) * &
        material_auxvars(ghosted_id)%porosity * &
        material_auxvars(ghosted_id)%volume*1000.d0
        
      if (iphase == 1) then
        ! add contribution of equilibrium sorption
        if (reaction%neqsorb > 0) then
          mass_balance(1:naqcomp,iphase) = mass_balance(1:naqcomp,iphase) + &
            rt_auxvars(ghosted_id)%total_sorb_eq(:) * &
            material_auxvars(ghosted_id)%volume
        endif

        ! add contribution from mineral volume fractions
        do imnrl = 1, reaction%mineral%nkinmnrl
          ncomp = reaction%mineral%kinmnrlspecid(0,imnrl)
          do i = 1, ncomp
            icomp = reaction%mineral%kinmnrlspecid(i,imnrl)
            mass_balance(icomp,iphase) = mass_balance(icomp,iphase) &
            + reaction%mineral%kinmnrlstoich(i,imnrl) &
            * rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl) &
            * material_auxvars(ghosted_id)%volume &
            / reaction%mineral%kinmnrl_molar_vol(imnrl)
          enddo
        enddo

        ! add contribution of immobile mass (still considered aqueous phase)
        do i = 1, reaction%immobile%nimmobile
          mass_balance(reaction%offset_immobile+i,iphase) = &
            mass_balance(reaction%offset_immobile+i,iphase) + &
            rt_auxvars(ghosted_id)%immobile(i) * &
            material_auxvars(ghosted_id)%volume
        enddo
      endif
    enddo
  enddo

end subroutine RTComputeMassBalance

! ************************************************************************** !

subroutine RTZeroMassBalanceDelta(realization)
  ! 
  ! Zeros mass balance delta array
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/08
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars_bc(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  rt_auxvars_bc => patch%aux%RT%auxvars_bc
  rt_auxvars_ss => patch%aux%RT%auxvars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%RT%num_aux
    patch%aux%RT%auxvars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  do iconn = 1, patch%aux%RT%num_aux_bc
    rt_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo

  do iconn = 1, patch%aux%RT%num_aux_ss
    rt_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine RTZeroMassBalanceDelta

! ************************************************************************** !

subroutine RTUpdateMassBalance(realization)
  ! 
  ! Updates mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/08
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars_bc(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  rt_auxvars_bc => patch%aux%RT%auxvars_bc
  rt_auxvars_ss => patch%aux%RT%auxvars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%RT%num_aux
    patch%aux%RT%auxvars(iconn)%mass_balance = &
      patch%aux%RT%auxvars(iconn)%mass_balance + &
      patch%aux%RT%auxvars(iconn)%mass_balance_delta*option%tran_dt
  enddo
#endif

  do iconn = 1, patch%aux%RT%num_aux_bc
    rt_auxvars_bc(iconn)%mass_balance = &
      rt_auxvars_bc(iconn)%mass_balance + &
      rt_auxvars_bc(iconn)%mass_balance_delta*option%tran_dt
  enddo

  do iconn = 1, patch%aux%RT%num_aux_ss
    rt_auxvars_ss(iconn)%mass_balance = &
      rt_auxvars_ss(iconn)%mass_balance + &
      rt_auxvars_ss(iconn)%mass_balance_delta*option%tran_dt
  enddo

end subroutine RTUpdateMassBalance

! ************************************************************************** !

subroutine RTInitializeTimestep(realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Realization_Subsurface_class

  type(realization_subsurface_type) :: realization
  
  call RTUpdateFixedAccumulation(realization)
  ! geh: never use transport coefs evaluated at time k
!  call RTUpdateTransportCoefs(realization)

end subroutine RTInitializeTimestep

! ************************************************************************** !

subroutine RTUpdateEquilibriumState(realization)
  ! 
  ! Updates equilibrium state variables after a
  ! successful time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/04/08
  ! 

  use Realization_Subsurface_class
  use Discretization_module
  use Patch_module
  use Option_module
  use Grid_module
  use Reaction_module
  !geh: please leave the "only" clauses for Secondary_Continuum_XXX as this
  !      resolves a bug in the Intel Visual Fortran compiler.
  use Secondary_Continuum_Aux_module, only : sec_transport_type
  use Secondary_Continuum_module, only : SecondaryRTUpdateEquilState
 
  implicit none

  type(realization_subsurface_type) :: realization

  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscInt :: ghosted_id, local_id
  PetscReal :: conc, max_conc, min_conc
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  reaction => realization%reaction
  grid => patch%grid

  call VecCopy(realization%field%tran_xx,realization%field%tran_yy, &
               ierr);CHKERRQ(ierr)
  call DiscretizationGlobalToLocal(realization%discretization, &
                                   realization%field%tran_xx, &
                                   realization%field%tran_xx_loc,NTRANDOF)
  
  rt_auxvars => patch%aux%RT%auxvars
  global_auxvars => patch%aux%Global%auxvars

  ! update:                        cells      bcs         act coefs
  call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_FALSE,PETSC_FALSE)

  ! update secondary continuum variables
  if (option%use_mc) then
    rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
        call SecondaryRTUpdateEquilState(rt_sec_transport_vars(local_id), &
                                          global_auxvars(ghosted_id), &
                                          reaction,option)                     
    enddo
  endif
  
end subroutine RTUpdateEquilibriumState

! ************************************************************************** !

subroutine RTUpdateKineticState(realization)
  ! 
  ! Updates kinetic state variables for reactive
  ! transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/27/13
  ! 

  use Realization_Subsurface_class
  use Discretization_module
  use Patch_module
  use Option_module
  use Grid_module
  use Reaction_module
  !geh: please leave the "only" clauses for Secondary_Continuum_XXX as this
  !      resolves a bug in the Intel Visual Fortran compiler.
  use Secondary_Continuum_Aux_module, only : sec_transport_type
  use Secondary_Continuum_module, only : SecondaryRTUpdateKineticState
 
  implicit none

  type(realization_subsurface_type) :: realization

  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscInt :: ghosted_id, local_id
  PetscReal :: conc, max_conc, min_conc
  PetscErrorCode :: ierr
  PetscReal :: sec_porosity
  
  option => realization%option
  patch => realization%patch
  reaction => realization%reaction
  grid => patch%grid

  rt_auxvars => patch%aux%RT%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  ! update mineral volume fractions, multirate sorption concentrations, 
  ! kinetic sorption concentration etc.  These updates must take place
  ! within reaction so that auxiliary variables are updated when only
  ! run in reaction mode.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

    if (.not.option%use_isothermal) then
      call RUpdateTempDependentCoefs(global_auxvars(ghosted_id),reaction, &
                                    PETSC_FALSE,option)
    endif

    call RUpdateKineticState(rt_auxvars(ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             reaction,option)
  enddo
  
  ! update secondary continuum variables
  if (option%use_mc) then
    rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
        sec_porosity = patch%material_property_array(1)%ptr% &
                        secondary_continuum_porosity

        call SecondaryRTUpdateKineticState(rt_sec_transport_vars(local_id), &
                                           global_auxvars(ghosted_id), &
                                           reaction,sec_porosity,option)                     
    enddo
  endif
  
end subroutine RTUpdateKineticState

! ************************************************************************** !

subroutine RTUpdateFixedAccumulation(realization)
  ! 
  ! Computes derivative of accumulation term in
  ! residual function
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/15/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Reactive_Transport_Aux_module
  use Option_module
  use Field_module  
  use Grid_module
  use Secondary_Continuum_Aux_module  

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reaction_type), pointer :: reaction
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscReal, pointer :: xx_p(:), accum_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: dof_offset, istart, iendaq, iendall
  PetscInt :: istartim, iendim
  PetscInt :: istartcoll, iendcoll
  PetscErrorCode :: ierr
  PetscReal :: vol_frac_prim
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  rt_auxvars => patch%aux%RT%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  
  grid => patch%grid
  reaction => realization%reaction
  if (option%use_mc) then
    rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars
  endif

  ! cannot use tran_xx_loc vector here as it has not yet been updated.
  call VecGetArrayReadF90(field%tran_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%tran_accum, accum_p, ierr);CHKERRQ(ierr)
  
  vol_frac_prim = 1.d0
  
! Do not use RTUpdateAuxVars() as it loops over ghosted ids

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    
    ! compute offset in solution vector for first dof in grid cell
    dof_offset = (local_id-1)*reaction%ncomp
    
    ! calculate range of aqueous species
    istart = dof_offset + 1
    iendaq = dof_offset + reaction%naqcomp
    iendall = dof_offset + reaction%ncomp

    ! copy primary aqueous species
    rt_auxvars(ghosted_id)%pri_molal = xx_p(istart:iendaq)
    
    if (reaction%ncoll > 0) then
      istartcoll = dof_offset + reaction%offset_colloid + 1
      iendcoll = dof_offset + reaction%offset_colloid + reaction%ncoll
      rt_auxvars(ghosted_id)%colloid%conc_mob = xx_p(istartcoll:iendcoll)* &
        global_auxvars(ghosted_id)%den_kg(1)*1.d-3
    endif
    
    if (reaction%immobile%nimmobile > 0) then
      istartim = dof_offset + reaction%offset_immobile + 1
      iendim = dof_offset + reaction%offset_immobile + reaction%immobile%nimmobile
      rt_auxvars(ghosted_id)%immobile = xx_p(istartim:iendim)
    endif

    if (.not.option%use_isothermal) then
      call RUpdateTempDependentCoefs(global_auxvars(ghosted_id),reaction, &
                                     PETSC_FALSE,option)
    endif
    
    ! DO NOT RECOMPUTE THE ACTIVITY COEFFICIENTS BEFORE COMPUTING THE
    ! FIXED PORTION OF THE ACCUMULATION TERM - geh
    call RTAuxVarCompute(rt_auxvars(ghosted_id), &
                         global_auxvars(ghosted_id), &
                         material_auxvars(ghosted_id), &
                         reaction,option)
    call RTAccumulation(rt_auxvars(ghosted_id), &
                        global_auxvars(ghosted_id), &
                        material_auxvars(ghosted_id), &
                        reaction,option, &
                        accum_p(istart:iendall)) 
    if (reaction%neqsorb > 0) then
      call RAccumulationSorb(rt_auxvars(ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             reaction,option,accum_p(istart:iendall))
    endif
        
    if (option%use_mc) then
      vol_frac_prim = rt_sec_transport_vars(local_id)%epsilon
      accum_p(istart:iendall) = accum_p(istart:iendall)*vol_frac_prim
    endif
    
  enddo

  call VecRestoreArrayReadF90(field%tran_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%tran_accum, accum_p, ierr);CHKERRQ(ierr)

end subroutine RTUpdateFixedAccumulation

! ************************************************************************** !

subroutine RTUpdateTransportCoefs(realization)
  ! 
  ! Calculates coefficients for transport matrix
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/24/10
  ! 

  use Realization_Subsurface_class
  use Discretization_module
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reactive_transport_param_type), pointer :: rt_parameter
  PetscInt :: local_id, ghosted_id, ghosted_face_id, id
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set  
  PetscInt :: sum_connection, iconn, num_connections
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal, allocatable :: cell_centered_Darcy_velocities(:,:)
  PetscReal, allocatable :: cell_centered_Darcy_velocities_ghosted(:,:,:)
  PetscReal ::  local_Darcy_velocities_up(3,2)
  PetscReal ::  local_Darcy_velocities_dn(3,2)
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: i
  PetscInt :: iphase
  PetscInt :: nphase
  PetscErrorCode :: ierr
    
  option => realization%option
  field => realization%field
  patch => realization%patch
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  nphase = rt_parameter%nphase
  
  local_Darcy_velocities_up = UNINITIALIZED_DOUBLE
  local_Darcy_velocities_dn = UNINITIALIZED_DOUBLE

  if (rt_parameter%calculate_transverse_dispersion) then
    allocate(cell_centered_Darcy_velocities_ghosted(3,nphase, &
                                                    patch%grid%ngmax))
    cell_centered_Darcy_velocities_ghosted = 0.d0
    allocate(cell_centered_Darcy_velocities(3,patch%grid%nlmax))
    do iphase = 1, nphase
      call PatchGetCellCenteredVelocities(patch,iphase, &
                                          cell_centered_Darcy_velocities)
      ! at this point, velocities are at local cell centers; we need 
      ! ghosted too.
      do i=1,3
        call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
        vec_ptr(:) = cell_centered_Darcy_velocities(i,:)
        call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
        call DiscretizationGlobalToLocal(realization%discretization, &
                                         field%work, &
                                         field%work_loc,ONEDOF)
        call VecGetArrayF90(field%work_loc,vec_ptr,ierr);CHKERRQ(ierr)
        cell_centered_Darcy_velocities_ghosted(i,iphase,:) = vec_ptr(:)
        call VecRestoreArrayF90(field%work_loc,vec_ptr,ierr);CHKERRQ(ierr)
      enddo
    enddo
    deallocate(cell_centered_Darcy_velocities)
  endif
  
  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle
          
      ! have to use temporary array since unallocated arrays cannot be
      ! indexed in call to subroutine.
      if (allocated(cell_centered_Darcy_velocities_ghosted)) then
        local_Darcy_velocities_up(:,1:nphase) = &
          cell_centered_Darcy_velocities_ghosted(:,1:nphase,ghosted_id_up)
        local_Darcy_velocities_dn(:,1:nphase) = &
          cell_centered_Darcy_velocities_ghosted(:,1:nphase,ghosted_id_dn)
      endif

      call TDispersion(global_auxvars(ghosted_id_up), &
                      material_auxvars(ghosted_id_up), &
                      local_Darcy_velocities_up, &
                      patch%material_property_array(patch%imat(ghosted_id_up))% &
                        ptr%dispersivity, &
                      global_auxvars(ghosted_id_dn), &
                      material_auxvars(ghosted_id_dn), &
                      local_Darcy_velocities_dn, &
                      patch%material_property_array(patch%imat(ghosted_id_dn))% &
                        ptr%dispersivity, &
                      cur_connection_set%dist(:,iconn), &
                      rt_parameter,option, &
                      patch%internal_velocities(:,sum_connection), &
                      patch%internal_tran_coefs(:,:,sum_connection))
                                           
    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
  
! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
 
    cur_connection_set => boundary_condition%connection_set
    num_connections = cur_connection_set%num_connections
    do iconn = 1, num_connections
      sum_connection = sum_connection + 1
  
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      if (allocated(cell_centered_Darcy_velocities_ghosted)) then
        local_Darcy_velocities_up(:,1:nphase) = &
          cell_centered_Darcy_velocities_ghosted(:,1:nphase,ghosted_id)
      endif
      
      call TDispersionBC(boundary_condition%tran_condition%itype, &
                        global_auxvars_bc(sum_connection), &
                        global_auxvars(ghosted_id), &
                        material_auxvars(ghosted_id), &
                        local_Darcy_velocities_up, &
                        patch%material_property_array(patch%imat(ghosted_id))% &
                          ptr%dispersivity, &
                        cur_connection_set%dist(:,iconn), &
                        rt_parameter,option, &
                        patch%boundary_velocities(:,sum_connection), &
                        patch%boundary_tran_coefs(:,:,sum_connection))
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
  if (allocated(cell_centered_Darcy_velocities_ghosted)) &
    deallocate(cell_centered_Darcy_velocities_ghosted)

end subroutine RTUpdateTransportCoefs

! ************************************************************************** !

subroutine RTUpdateRHSCoefs(realization)
  ! 
  ! Updates coefficients for the right hand side of
  ! linear transport equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/25/10
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  PetscReal, pointer :: rhs_coef_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: iphase
  PetscErrorCode :: ierr
    
  option => realization%option
  field => realization%field
  patch => realization%patch
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  grid => patch%grid

  ! Get vectors
  call VecGetArrayF90(field%tran_rhs_coef,rhs_coef_p,ierr);CHKERRQ(ierr)

  iphase = 1
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    rhs_coef_p(local_id) = material_auxvars(ghosted_id)%porosity* &
                           global_auxvars(ghosted_id)%sat(iphase)* &
! total already has den_kg within 
!                           global_auxvars(ghosted_id)%den_kg(iphase)* &
                           1000.d0* &
                           material_auxvars(ghosted_id)%volume/option%tran_dt
  enddo

  ! Restore vectors
  call VecRestoreArrayF90(field%tran_rhs_coef,rhs_coef_p,ierr);CHKERRQ(ierr)

end subroutine RTUpdateRHSCoefs

! ************************************************************************** !

subroutine RTCalculateRHS_t0(realization)
  ! 
  ! Calculate porition of RHS of transport system
  ! at time t0 or time level k
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/25/10
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reaction_type), pointer :: reaction
  PetscReal, pointer :: rhs_coef_p(:)
  PetscReal, pointer :: rhs_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: iphase
  PetscInt :: istartaq, iendaq

  PetscErrorCode :: ierr
    
  option => realization%option
  field => realization%field
  patch => realization%patch
  rt_auxvars => patch%aux%RT%auxvars
  grid => patch%grid
  reaction => realization%reaction

  ! Get vectors
  call VecGetArrayReadF90(field%tran_rhs_coef,rhs_coef_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%tran_rhs,rhs_p,ierr);CHKERRQ(ierr)

  iphase = 1
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle    
    iendaq = local_id*reaction%naqcomp
    istartaq = iendaq-reaction%naqcomp+1
    rhs_p(istartaq:iendaq) = rt_auxvars(ghosted_id)%total(:,iphase)* &
                             rhs_coef_p(local_id)
  enddo

  ! Restore vectors
  call VecRestoreArrayReadF90(field%tran_rhs_coef,rhs_coef_p, &
                              ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%tran_rhs,rhs_p,ierr);CHKERRQ(ierr)

end subroutine RTCalculateRHS_t0

! ************************************************************************** !

subroutine RTCalculateRHS_t1(realization)
  ! 
  ! Calculate porition of RHS of transport system
  ! at time level k+1
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/25/10
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reaction_type), pointer :: reaction
  PetscReal, pointer :: rhs_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: iphase
  PetscReal :: coef_up(1), coef_dn(1)
  PetscReal :: msrc(2)
  PetscReal :: Res(realization%reaction%naqcomp)
  PetscInt :: istartaq, iendaq

  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: source_sink
  type(reactive_transport_param_type), pointer :: rt_parameter
  PetscInt :: sum_connection, iconn  
  PetscReal :: qsrc(2)
  PetscInt :: offset, istartcoll, iendcoll, istartall, iendall, icomp, ieqgas
  PetscBool :: volumetric
  PetscInt :: flow_src_sink_type
  PetscReal :: coef_in(2), coef_out(2)
  PetscInt :: nphase
  PetscErrorCode :: ierr
    
  option => realization%option
  field => realization%field
  patch => realization%patch
  rt_auxvars => patch%aux%RT%auxvars
  rt_auxvars_bc => patch%aux%RT%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  grid => patch%grid
  reaction => realization%reaction
  rt_parameter => patch%aux%RT%rt_parameter
  nphase = rt_parameter%nphase
  
  iphase = 1
  option%io_buffer = 'RTCalculateRHS_t1 must be refactored'
  call printErrMsg(option)

  ! Get vectors
  call VecGetArrayF90(field%tran_rhs,rhs_p,ierr);CHKERRQ(ierr)

  ! add in inflowing boundary conditions
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
  
    cur_connection_set => boundary_condition%connection_set
  
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
  
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      call TFluxCoef(rt_parameter,option,cur_connection_set%area(iconn), &
                     patch%boundary_velocities(:,sum_connection), &
                     patch%boundary_tran_coefs(:,:,sum_connection), &
                     0.5d0, & ! fraction upwind (0.d0 upwind, 0.5 central)
                     coef_up,coef_dn)

      ! coef_dn not needed 
      iendaq = local_id*reaction%naqcomp
      istartaq = iendaq-reaction%naqcomp+1
      
      rhs_p(istartaq:iendaq) = rhs_p(istartaq:iendaq) + &
        coef_up(iphase)*rt_auxvars_bc(sum_connection)%total(:,iphase)

    enddo
    boundary_condition => boundary_condition%next
  enddo  

  ! add in inflowing sources
#if 1
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    qsrc = 0.d0
    flow_src_sink_type = 0
    if (associated(source_sink%flow_condition) .and. &
        associated(source_sink%flow_condition%rate)) then
      qsrc = source_sink%flow_condition%rate%dataset%rarray(1)
      flow_src_sink_type = source_sink%flow_condition%rate%itype
    endif
    
    ! only handle injection on rhs
    if (qsrc(1) < 0.d0) then
      source_sink => source_sink%next
      cycle
    endif
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      offset = (local_id-1)*reaction%ncomp

      if (patch%imat(ghosted_id) <= 0) cycle
      
      istartaq = reaction%offset_aqueous + 1
      iendaq = reaction%offset_aqueous + reaction%naqcomp
      
      if (reaction%ncoll > 0) then
        istartcoll = reaction%offset_colloid + 1
        iendcoll = reaction%offset_colloid + reaction%ncoll
      endif

      if (associated(patch%ss_flow_vol_fluxes)) then
        qsrc(1:nphase) = patch%ss_flow_vol_fluxes(1:nphase,sum_connection)
      endif
      call TSrcSinkCoef(rt_parameter,qsrc,source_sink%tran_condition%itype, &
                        coef_in,coef_out)      

      Res(istartaq:iendaq) = & !coef_in*rt_auxvars(ghosted_id)%total(:,iphase) + &
                             coef_out*source_sink%tran_condition%cur_constraint_coupler% &
                                        rt_auxvar%total(:,iphase)
      if (reaction%ncoll > 0) then
        Res(istartcoll:iendcoll) = & !coef_in*rt_auxvars(ghosted_id)%colloid%conc_mob(:) !+ &
                                   coef_out*source_sink%tran_condition%cur_constraint_coupler% &
                                              rt_auxvar%colloid%conc_mob(:)
      endif
      istartall = offset + 1
      iendall = offset + reaction%ncomp
      ! subtract since the contribution is on the rhs
      rhs_p(istartall:iendall) = rhs_p(istartall:iendall) - Res(1:reaction%ncomp)                                  
    enddo
    source_sink => source_sink%next
  enddo

#endif

  ! Restore vectors
  call VecRestoreArrayF90(field%tran_rhs,rhs_p,ierr);CHKERRQ(ierr)

  ! Mass Transfer
  if (field%tran_mass_transfer /= PETSC_NULL_VEC) then
    ! scale by -1.d0 for contribution to residual.  A negative contribution
    ! indicates mass being added to system.
    call VecAXPY(field%tran_rhs,-1.d0,field%tran_mass_transfer, &
                 ierr);CHKERRQ(ierr)
  endif
  
end subroutine RTCalculateRHS_t1

! ************************************************************************** !

subroutine RTCalculateTransportMatrix(realization,T)
  ! 
  ! Calculate transport matrix
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/25/10
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Grid_module
  use Patch_module
  use Field_module
  use Coupler_module
  use Connection_module
  use Debug_module

  implicit none
      
  type(realization_subsurface_type) :: realization
  Mat :: T
  
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: iphase
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: source_sink
  type(reactive_transport_param_type), pointer :: rt_parameter
  PetscInt :: sum_connection, iconn
  PetscReal :: coef
  PetscReal :: coef_up(2), coef_dn(2)
  PetscReal :: qsrc(2)
  PetscBool :: volumetric  
  PetscInt :: flow_pc
  PetscInt :: flow_src_sink_type
  PetscReal :: coef_in(2), coef_out(2)
  PetscViewer :: viewer
  PetscInt :: nphase
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: string

  option => realization%option
  field => realization%field
  patch => realization%patch
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  grid => patch%grid  
  rt_parameter => patch%aux%RT%rt_parameter

  nphase = rt_parameter%nphase
    
  call MatZeroEntries(T,ierr);CHKERRQ(ierr)
  
  ! Get vectors

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle

      call TFluxCoef(rt_parameter,option,cur_connection_set%area(iconn), &
                     patch%internal_velocities(:,sum_connection), &
                     patch%internal_tran_coefs(:,:,sum_connection), &
                     cur_connection_set%dist(-1,iconn), &
                     coef_up,coef_dn)

!      coef_up = coef_up*global_auxvars(ghosted_id_up)%den_kg*1.d-3
!      coef_dn = coef_dn*global_auxvars(ghosted_id_dn)%den_kg*1.d-3
                     
      if (local_id_up > 0) then
        call MatSetValuesLocal(T,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                               coef_up,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesLocal(T,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                               coef_dn,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
      if (local_id_dn > 0) then
        coef_up = -coef_up
        coef_dn = -coef_dn
        call MatSetValuesLocal(T,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                               coef_dn,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesLocal(T,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                               coef_up,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
                     
    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
  
  ! add in outflowing boundary conditions
  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
  
    cur_connection_set => boundary_condition%connection_set
  
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
  
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      call TFluxCoef(rt_parameter,option,cur_connection_set%area(iconn), &
                     patch%boundary_velocities(:,sum_connection), &
                     patch%boundary_tran_coefs(:,:,sum_connection), &
                     0.5d0, & ! fraction upwind (0.d0 upwind, 0.5 central)
                     coef_up,coef_dn)

 !     coef_dn = coef_dn*global_auxvars(ghosted_id)%den_kg*1.d-3

      !Jup not needed 
      coef_dn = -coef_dn
      call MatSetValuesLocal(T,1,ghosted_id-1,1,ghosted_id-1,coef_dn, &
                             ADD_VALUES,ierr);CHKERRQ(ierr)
    
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
  ! Accumulation term
  iphase = 1
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle    
    coef = material_auxvars(ghosted_id)%porosity* &
           global_auxvars(ghosted_id)%sat(iphase)* &
!geh           global_auxvars(ghosted_id)%den_kg(iphase)* &
           1000.d0* &
           material_auxvars(ghosted_id)%volume/option%tran_dt
    call MatSetValuesLocal(T,1,ghosted_id-1,1,ghosted_id-1,coef, &
                           ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo
                        
  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    qsrc = 0.d0
    flow_src_sink_type = 0
    if (associated(source_sink%flow_condition) .and. &
        associated(source_sink%flow_condition%rate)) then
      qsrc = source_sink%flow_condition%rate%dataset%rarray(1)
      flow_src_sink_type = source_sink%flow_condition%rate%itype
    endif
      
    ! only handle extraction on lhs
    if (qsrc(1) > 0.d0) then
      source_sink => source_sink%next
      cycle
    endif
      
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      if (associated(patch%ss_flow_vol_fluxes)) then
        qsrc(1:nphase) = patch%ss_flow_vol_fluxes(1:nphase,sum_connection)
      endif
      call TSrcSinkCoef(rt_parameter,qsrc,source_sink%tran_condition%itype, &
                        coef_in,coef_out)

      coef_dn = coef_in
      !geh: do not remove this conditional as otherwise MatSetValuesLocal() 
      !     will be called for injection too (wasted calls)
      if (coef_dn(1) > 0.d0) then
        call MatSetValuesLocal(T,1,ghosted_id-1,1,ghosted_id-1,coef_dn, &
                               ADD_VALUES,ierr);CHKERRQ(ierr)
      endif 

    enddo
    source_sink => source_sink%next
  enddo

  ! All CO2 source/sinks are handled on the RHS for now

  call MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  if (patch%aux%RT%inactive_cells_exist) then
    coef = 1.d0
    call MatZeroRowsLocal(T,patch%aux%RT%n_zero_rows, &
                          patch%aux%RT%zero_rows_local_ghosted,coef, &
                          PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
  endif

  if (realization%debug%matview_Jacobian) then
    string = 'Tmatrix'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(T,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif  
  
end subroutine RTCalculateTransportMatrix

! ************************************************************************** !

subroutine RTReact(realization)
  ! 
  ! Calculate reaction
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/03/10
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  
  use Secondary_Continuum_Aux_module
  use Logging_module
  
!$ use omp_lib
     
  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reaction_type), pointer :: reaction
  type(option_type), pointer :: option
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend, iendaq
  PetscInt :: iphase
  PetscInt :: ithread, vector_length
  PetscReal, pointer :: tran_xx_p(:)
  PetscReal, pointer :: mask_p(:)
  PetscInt :: num_iterations
#ifdef OS_STATISTICS
  PetscInt :: sum_iterations
  PetscInt :: max_iterations
#endif
  PetscInt :: icount
  PetscErrorCode :: ierr

#ifdef OS_STATISTICS
  PetscInt :: call_count
  PetscInt :: sum_newton_iterations
  PetscReal :: ave_newton_iterations_in_a_cell
  PetscInt :: max_newton_iterations_in_a_cell
  PetscInt :: max_newton_iterations_on_a_core
  PetscInt :: min_newton_iterations_on_a_core
  PetscInt :: temp_int_in(3)
  PetscInt :: temp_int_out(3)
#endif
  
  call PetscLogEventBegin(logging%event_rt_react,ierr);CHKERRQ(ierr)
                          
#ifdef OS_STATISTICS
  call_count = 0
  sum_newton_iterations = 0
  max_newton_iterations_in_a_cell = -99999999
  max_newton_iterations_on_a_core = -99999999
  min_newton_iterations_on_a_core = 99999999
#endif   
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  global_auxvars => patch%aux%Global%auxvars
  rt_auxvars => patch%aux%RT%auxvars
  material_auxvars => patch%aux%Material%auxvars
  grid => patch%grid
  reaction => realization%reaction

  ! need up update aux vars based on current density/saturation,
  ! but NOT activity coefficients
  call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_FALSE,PETSC_FALSE)

  ! Get vectors
  call VecGetArrayReadF90(field%tran_xx,tran_xx_p,ierr);CHKERRQ(ierr)
      
  iphase = 1
  ithread = 1
#ifdef OS_STATISTICS
  sum_iterations = 0
  max_iterations = 0
  icount = 0
#endif

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    
    istart = (local_id-1)*reaction%ncomp+1
    iend = istart + reaction%ncomp - 1
    iendaq = istart + reaction%naqcomp - 1
    
    
    call RReact(rt_auxvars(ghosted_id),global_auxvars(ghosted_id), &
                material_auxvars(ghosted_id), &
                tran_xx_p(istart:iend), &
                num_iterations,reaction,option)
    ! set primary dependent var back to free-ion molality
    tran_xx_p(istart:iendaq) = rt_auxvars(ghosted_id)%pri_molal
    if (reaction%immobile%nimmobile > 0) then
      tran_xx_p(reaction%offset_immobile: &
                reaction%offset_immobile + reaction%immobile%nimmobile) = &
        rt_auxvars(ghosted_id)%immobile
    endif
#ifdef OS_STATISTICS
    if (num_iterations > max_iterations) then
      max_iterations = num_iterations
    endif
    sum_iterations = sum_iterations + num_iterations
    icount = icount + 1
#endif
  enddo
  
#ifdef OS_STATISTICS
  patch%aux%RT%rt_parameter%newton_call_count = icount
  patch%aux%RT%rt_parameter%sum_newton_call_count = &
    patch%aux%RT%rt_parameter%sum_newton_call_count + dble(icount)
  patch%aux%RT%rt_parameter%newton_iterations = sum_iterations
  patch%aux%RT%rt_parameter%sum_newton_iterations = &
    patch%aux%RT%rt_parameter%sum_newton_iterations + dble(sum_iterations)
  patch%aux%RT%rt_parameter%max_newton_iterations = max_iterations
#endif
  
  ! Restore vectors
  call VecRestoreArrayReadF90(field%tran_xx,tran_xx_p,ierr);CHKERRQ(ierr)

  if (option%compute_mass_balance_new) then
    call RTZeroMassBalanceDelta(realization)
    call RTComputeBCMassBalanceOS(realization)
  endif
  
#ifdef OS_STATISTICS
  call_count = call_count + cur_patch%aux%RT%rt_parameter%newton_call_count
  sum_newton_iterations = sum_newton_iterations + &
    cur_patch%aux%RT%rt_parameter%newton_iterations
  if (cur_patch%aux%RT%rt_parameter%max_newton_iterations > &
      max_newton_iterations_in_a_cell) then
    max_newton_iterations_in_a_cell = &
      cur_patch%aux%RT%rt_parameter%max_newton_iterations
  endif
  if (cur_patch%aux%RT%rt_parameter%max_newton_iterations > &
      cur_patch%aux%RT%rt_parameter%overall_max_newton_iterations) then
    cur_patch%aux%RT%rt_parameter%overall_max_newton_iterations = &
      cur_patch%aux%RT%rt_parameter%max_newton_iterations
  endif
#endif 

  ! Logging must come before statistics since the global reductions
  ! will synchonize the cores
  call PetscLogEventEnd(logging%event_rt_react,ierr);CHKERRQ(ierr)
                        
#ifdef OS_STATISTICS
  temp_int_in(1) = call_count
  temp_int_in(2) = sum_newton_iterations
  call MPI_Allreduce(temp_int_in,temp_int_out,TWO_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  ave_newton_iterations_in_a_cell = float(temp_int_out(2)) / temp_int_out(1)

  temp_int_in(1) = max_newton_iterations_in_a_cell
  temp_int_in(2) = sum_newton_iterations ! to calc max # iteration on a core
  temp_int_in(3) = -sum_newton_iterations ! to calc min # iteration on a core
  call MPI_Allreduce(temp_int_in,temp_int_out,THREE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  max_newton_iterations_in_a_cell = temp_int_out(1)
  max_newton_iterations_on_a_core = temp_int_out(2)
  min_newton_iterations_on_a_core = -temp_int_out(3)
  
  if (option%print_screen_flag) then
    write(*, '(" OS Reaction Statistics: ",/, &
             & "   Ave Newton Its / Cell: ",1pe12.4,/, &
             & "   Max Newton Its / Cell: ",i4,/, &
             & "   Max Newton Its / Core: ",i6,/, &
             & "   Min Newton Its / Core: ",i6)') &
               ave_newton_iterations_in_a_cell, &
               max_newton_iterations_in_a_cell, &
               max_newton_iterations_on_a_core, &
               min_newton_iterations_on_a_core
  endif

  if (option%print_file_flag) then
    write(option%fid_out, '(" OS Reaction Statistics: ",/, &
             & "   Ave Newton Its / Cell: ",1pe12.4,/, &
             & "   Max Newton Its / Cell: ",i4,/, &
             & "   Max Newton Its / Core: ",i6,/, &
             & "   Min Newton Its / Core: ",i6)') &
               ave_newton_iterations_in_a_cell, &
               max_newton_iterations_in_a_cell, &
               max_newton_iterations_on_a_core, &
               min_newton_iterations_on_a_core
  endif

#endif   

end subroutine RTReact

! ************************************************************************** !

subroutine RTComputeBCMassBalanceOS(realization)
  ! 
  ! Calculates mass balance at boundary
  ! conditions for operator split mode
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/04/10
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Transport_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  
  implicit none

  type(realization_subsurface_type) :: realization  

  PetscInt :: local_id, ghosted_id
  PetscInt, parameter :: iphase = 1
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(reactive_transport_param_type), pointer :: rt_parameter
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars_bc(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:) 
  type(global_auxvar_type), pointer :: global_auxvars_ss(:) 
  PetscReal :: Res(realization%reaction%ncomp)
  
  PetscReal, pointer :: face_fluxes_p(:)

  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: flow_src_sink_type
  PetscReal :: qsrc(2)
  
  PetscReal :: coef_up(realization%option%transport%nphase)
  PetscReal :: coef_dn(realization%option%transport%nphase)
  PetscReal :: coef_in(2), coef_out(2)
  PetscInt :: nphase
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  patch => realization%patch
  reaction => realization%reaction
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  rt_auxvars => patch%aux%RT%auxvars
  rt_auxvars_bc => patch%aux%RT%auxvars_bc
  rt_auxvars_ss => patch%aux%RT%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  nphase = rt_parameter%nphase
  
! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
  
    cur_connection_set => boundary_condition%connection_set
  
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
  
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      ! TFluxCoef accomplishes the same as what TBCCoef would
      call TFluxCoef(rt_parameter,option,cur_connection_set%area(iconn), &
                     patch%boundary_velocities(:,sum_connection), &
                     patch%boundary_tran_coefs(:,:,sum_connection), &
                     0.5d0, &
                     coef_up,coef_dn)
      ! TFlux accomplishes the same as what TBCFlux would
      call TFlux(rt_parameter, &
                 rt_auxvars_bc(sum_connection), &
                 global_auxvars_bc(sum_connection), &
                 rt_auxvars(ghosted_id), &
                 global_auxvars(ghosted_id), &
                 coef_up,coef_dn,option,Res)

    ! contribution to boundary 
      rt_auxvars_bc(sum_connection)%mass_balance_delta(:,iphase) = &
        rt_auxvars_bc(sum_connection)%mass_balance_delta(:,iphase) - Res
!        ! contribution to internal 
!        rt_auxvars(ghosted_id)%mass_balance_delta(:,iphase) = &
!          rt_auxvars(ghosted_id)%mass_balance_delta(:,iphase) + Res
    
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set

    flow_src_sink_type = 0
    if (associated(source_sink%flow_condition) .and. &
        associated(source_sink%flow_condition%rate)) then
      qsrc = source_sink%flow_condition%rate%dataset%rarray(1)
      flow_src_sink_type = source_sink%flow_condition%rate%itype
    endif
      
    do iconn = 1, cur_connection_set%num_connections 
      sum_connection = sum_connection + 1     
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      if (associated(patch%ss_flow_vol_fluxes)) then
        qsrc(1:nphase) = patch%ss_flow_vol_fluxes(1:nphase,sum_connection)
      endif
      call TSrcSinkCoef(rt_parameter,qsrc,source_sink%tran_condition%itype, &
                        coef_in,coef_out)
      
      Res = coef_in*rt_auxvars(ghosted_id)%total(:,iphase) + &
            coef_out*source_sink%tran_condition%cur_constraint_coupler% &
            rt_auxvar%total(:,iphase)
      if (option%compute_mass_balance_new) then
        ! contribution to boundary 
        rt_auxvars_ss(sum_connection)%mass_balance_delta(:,iphase) = &
          rt_auxvars_ss(sum_connection)%mass_balance_delta(:,iphase) + Res
        ! contribution to internal 
!        rt_auxvars(ghosted_id)%mass_balance_delta(:,iphase) = &
!          rt_auxvars(ghosted_id)%mass_balance_delta(:,iphase) - Res
        endif
    enddo
    source_sink => source_sink%next
  enddo

end subroutine RTComputeBCMassBalanceOS

! ************************************************************************** !

subroutine RTNumericalJacobianTest(realization)
  ! 
  ! Computes the a test numerical jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/20/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Field_module

  implicit none

  Vec :: xx
  type(realization_subsurface_type) :: realization

  Vec :: xx_pert
  Vec :: res
  Vec :: res_pert
  Mat :: A
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  PetscReal :: derivative, perturbation
  
  PetscReal, pointer :: vec_p(:), vec2_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  
  PetscInt :: idof, idof2, icell

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  call VecDuplicate(field%tran_xx,xx_pert,ierr);CHKERRQ(ierr)
  call VecDuplicate(field%tran_xx,res,ierr);CHKERRQ(ierr)
  call VecDuplicate(field%tran_xx,res_pert,ierr);CHKERRQ(ierr)
  
  call MatCreate(option%mycomm,A,ierr);CHKERRQ(ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE, &
                   grid%nlmax*option%ntrandof, &
                   grid%nlmax*option%ntrandof,ierr);CHKERRQ(ierr)
  call MatSetType(A,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetFromOptions(A,ierr);CHKERRQ(ierr)
    
  call RTResidual(PETSC_NULL_SNES,field%tran_xx,res,realization,ierr)
  call VecGetArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)
  do idof = 1,grid%nlmax*option%ntrandof
    icell = (idof-1)/option%ntrandof+1
    if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    call VecCopy(field%tran_xx,xx_pert,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
    perturbation = vec_p(idof)*perturbation_tolerance
    vec_p(idof) = vec_p(idof)+perturbation
    call vecrestorearrayf90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
    call RTResidual(PETSC_NULL_SNES,xx_pert,res_pert,realization,ierr)
    call vecgetarrayf90(res_pert,vec_p,ierr);CHKERRQ(ierr)
    do idof2 = 1, grid%nlmax*option%ntrandof
      derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
      if (dabs(derivative) > 1.d-30) then
        call matsetvalue(a,idof2-1,idof-1,derivative,insert_values, &
                         ierr);CHKERRQ(ierr)
      endif
    enddo
    call VecRestoreArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
  enddo
  call VecRestoreArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call PetscViewerASCIIOpen(option%mycomm,'RTnumerical_jacobian.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call MatView(A,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  call MatDestroy(A,ierr);CHKERRQ(ierr)
  
  call VecDestroy(xx_pert,ierr);CHKERRQ(ierr)
  call VecDestroy(res,ierr);CHKERRQ(ierr)
  call VecDestroy(res_pert,ierr);CHKERRQ(ierr)
  
end subroutine RTNumericalJacobianTest

! ************************************************************************** !

subroutine RTResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Logging_module
  use Debug_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  PetscReal, pointer :: xx_p(:), log_xx_p(:)
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  PetscViewer :: viewer  
  
  character(len=MAXSTRINGLENGTH) :: string

  call PetscLogEventBegin(logging%event_rt_residual,ierr);CHKERRQ(ierr)

  patch => realization%patch
  field => realization%field
  discretization => realization%discretization
  option => realization%option

  ! Communication -----------------------------------------
  if (realization%reaction%use_log_formulation) then
    ! have to convert the log concentration to non-log form
    call VecGetArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(xx,log_xx_p,ierr);CHKERRQ(ierr)
    xx_p(:) = exp(log_xx_p(:))
    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(xx,log_xx_p,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                     field%tran_xx_loc,NTRANDOF)
  else
    call DiscretizationGlobalToLocal(discretization,xx,field%tran_xx_loc, &
                                     NTRANDOF)
  endif

  ! pass #1 for internal and boundary flux terms
  call RTResidualFlux(snes,xx,r,realization,ierr)

  ! pass #2 for everything else
  call RTResidualNonFlux(snes,xx,r,realization,ierr)

  if (realization%debug%vecview_residual) then
    string = 'RTresidual'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%vecview_solution) then
    string = 'RTxx'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(field%tran_xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  
  call PetscLogEventEnd(logging%event_rt_residual,ierr);CHKERRQ(ierr)

end subroutine RTResidual

! ************************************************************************** !

subroutine RTResidualFlux(snes,xx,r,realization,ierr)
  ! 
  ! Computes the flux terms in the residual function for
  ! reactive transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Transport_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  use Secondary_Continuum_Aux_module
  
  implicit none

  type :: flux_ptrs
    PetscReal, dimension(:), pointer :: flux_p 
  end type

  type (flux_ptrs), dimension(0:2) :: fluxes
  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(inout) :: r
  type(realization_subsurface_type) :: realization  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: r_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt, parameter :: iphase = 1
  PetscInt :: i, istart, iend                        
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(reactive_transport_param_type), pointer :: rt_parameter
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:), rt_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 
  
  PetscReal, pointer :: face_fluxes_p(:)

  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn
  
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscReal :: vol_frac_prim

  PetscInt :: ii
  character(len=MAXSTRINGLENGTH) :: string
  
#ifdef CENTRAL_DIFFERENCE  
  PetscReal :: T_11(realization%option%transport%nphase)
  PetscReal :: T_12(realization%option%transport%nphase)
  PetscReal :: T_21(realization%option%transport%nphase)
  PetscReal :: T_22(realization%option%transport%nphase)
  PetscReal :: Res_1(realization%reaction%ncomp)
  PetscReal :: Res_2(realization%reaction%ncomp)
#else
  PetscReal :: coef_up(realization%patch%aux%RT%rt_parameter%naqcomp, &
                       realization%option%transport%nphase)
  PetscReal :: coef_dn(realization%patch%aux%RT%rt_parameter%naqcomp, &
                       realization%option%transport%nphase)
  PetscReal :: Res(realization%reaction%ncomp)
#endif

  ! CO2-specific
  PetscReal :: msrc(1:realization%option%nflowspec)
  PetscInt :: icomp, ieqgas

  ! zero out xy-direction transport
  PetscReal :: unitvec_xyz(3)

  option => realization%option
  field => realization%field
  patch => realization%patch
  reaction => realization%reaction
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  rt_auxvars => patch%aux%RT%auxvars
  rt_auxvars_bc => patch%aux%RT%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  if (option%use_mc) then
    rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars
  endif

  if (reaction%act_coef_update_frequency == &
      ACT_COEF_FREQUENCY_NEWTON_ITER) then
    ! update:                        cells      bcs        act. coefs.
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)
  else
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE,PETSC_FALSE)
  endif
  
  if (option%compute_mass_balance_new) then
    call RTZeroMassBalanceDelta(realization)
  endif
  
  ! Get pointer to Vector data
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
 
  r_p = 0.d0
  vol_frac_prim = 1.d0

  ! Interior Flux Terms -----------------------------------
  unitvec_xyz(1:2) = 0.d0
  unitvec_xyz(3)   = 1.d0

  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle

      !----------------
      if ( option%flow%only_vertical_flow .or. option%transport%only_vertical_tran ) then
         if(abs(dot_product(cur_connection_set%dist(1:3,iconn),unitvec_xyz)) < 1.d-20) cycle
      end if
      !----------------

      ! TFluxCoef will eventually be moved to another routine where it should be
      ! called only once per flux interface at the beginning of a transport
      ! time step.
      
      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(local_id_up)%epsilon
      endif  
      
      
#ifndef CENTRAL_DIFFERENCE        
      call TFluxCoef(rt_parameter,option,cur_connection_set%area(iconn), &
                patch%internal_velocities(:,sum_connection), &
                patch%internal_tran_coefs(:,:,sum_connection)*vol_frac_prim, &
                cur_connection_set%dist(-1,iconn), &
                coef_up,coef_dn)
                      
      call TFlux(rt_parameter, &
                  rt_auxvars(ghosted_id_up), &
                  global_auxvars(ghosted_id_up), &
                  rt_auxvars(ghosted_id_dn), &
                  global_auxvars(ghosted_id_dn), &
                  coef_up,coef_dn,option,Res)

      ! checking NaN or INF
      do ii=1,reaction%ncomp
        if(Res(ii) /= Res(ii) .or. &
          abs(Res(ii))>huge(Res(ii)) ) then
          write(string,*) 'local_id: ', local_id_up, 'Res: ', ii, Res(ii)
          option%io_buffer = ' NaN or INF of Residuals @ reactive_transport.F90: RTResidualFlux - Interior Flux ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo


#ifdef COMPUTE_INTERNAL_MASS_FLUX
      rt_auxvars(local_id_up)%mass_balance_delta(:,iphase) = &
        rt_auxvars(local_id_up)%mass_balance_delta(:,iphase) - Res        
#endif
      if (local_id_up>0) then
        iend = local_id_up*reaction%ncomp
        istart = iend-reaction%ncomp+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:reaction%ncomp)
      endif
      
      if (local_id_dn>0) then
        iend = local_id_dn*reaction%ncomp
        istart = iend-reaction%ncomp+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:reaction%ncomp)
      endif

      if (associated(patch%internal_tran_fluxes)) then
        patch%internal_tran_fluxes(1:reaction%ncomp,iconn) = &
            Res(1:reaction%ncomp)
      endif

#else
      call TFluxCoef_CD(option,cur_connection_set%area(iconn), &
                 patch%internal_velocities(:,sum_connection), &
                 patch%internal_tran_coefs(:,:,sum_connection)*vol_frac_prim, &
                 cur_connection_set%dist(-1,iconn), &
                 T_11,T_12,T_21,T_22)
      call TFlux_CD(rt_parameter, &
                  rt_auxvars(ghosted_id_up), &
                  global_auxvars(ghosted_id_up), &
                  rt_auxvars(ghosted_id_dn), &
                  global_auxvars(ghosted_id_dn), &
                  T_11,T_12,T_21,T_22,option,Res_1,Res_2)

      ! checking NaN or INF
      do ii=1,reaction%ncomp
        if(Res_1(ii) /= Res_1(ii) .or. &
          abs(Res_1(ii))>huge(Res_1(ii)) ) then
          write(string,*) 'local_id: ', local_id_up, 'Res: ', ii, Res_1(ii)
          option%io_buffer = ' NaN or INF of Residuals @ reactive_transport.F90: RTResidualFlux - Interior Upwind ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo
      do ii=1,reaction%ncomp
        if(Res_2(ii) /= Res_2(ii) .or. &
          abs(Res_2(ii))>huge(Res_2(ii)) ) then
          write(string,*) 'local_id: ', local_id_dn, 'Res: ', ii, Res_2(ii)
          option%io_buffer = ' NaN or INF of Residuals @ reactive_transport.F90: RTResidualFlux - Interior Downwind ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo

      if (local_id_up>0) then
        iend = local_id_up*reaction%ncomp
        istart = iend-reaction%ncomp+1
        r_p(istart:iend) = r_p(istart:iend) + Res_1(1:reaction%ncomp)
      endif
      
      if (local_id_dn>0) then
        iend = local_id_dn*reaction%ncomp
        istart = iend-reaction%ncomp+1
        r_p(istart:iend) = r_p(istart:iend) + Res_2(1:reaction%ncomp)
      endif

      if (associated(patch%internal_tran_fluxes)) then
        patch%internal_tran_fluxes(1:reaction%ncomp,iconn) = &
            Res_1(1:reaction%ncomp) + Res_2(1:reaction%ncomp)
      endif

#endif


    enddo
    cur_connection_set => cur_connection_set%next
  enddo
    
! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(local_id)%epsilon
      endif  
      
#ifndef CENTRAL_DIFFERENCE
      ! TFluxCoef accomplishes the same as what TBCCoef would
      call TFluxCoef(rt_parameter,option,cur_connection_set%area(iconn), &
                  patch%boundary_velocities(:,sum_connection), &
                  patch%boundary_tran_coefs(:,:,sum_connection)*vol_frac_prim, &
                  0.5d0, &
                  coef_up,coef_dn)
      ! TFlux accomplishes the same as what TBCFlux would
      call TFlux(rt_parameter, &
                  rt_auxvars_bc(sum_connection), &
                  global_auxvars_bc(sum_connection), &
                  rt_auxvars(ghosted_id), &
                  global_auxvars(ghosted_id), &
                  coef_up,coef_dn,option,Res)
      iend = local_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      r_p(istart:iend)= r_p(istart:iend) - Res(1:reaction%ncomp)

      if (option%compute_mass_balance_new) then
      ! contribution to boundary 
        rt_auxvars_bc(sum_connection)%mass_balance_delta(:,iphase) = &
          rt_auxvars_bc(sum_connection)%mass_balance_delta(:,iphase) - Res
!        ! contribution to internal 
!        rt_auxvars(ghosted_id)%mass_balance_delta(:,iphase) = &
!          rt_auxvars(ghosted_id)%mass_balance_delta(:,iphase) + Res
        endif  

      if (associated(patch%boundary_tran_fluxes)) then
        patch%boundary_tran_fluxes(1:reaction%ncomp,sum_connection) = &
            -Res(1:reaction%ncomp)

      ! checking NaN or INF
      do ii=1,reaction%ncomp
        if(Res(ii) /= Res(ii) .or. &
          abs(Res(ii))>huge(Res(ii)) ) then
          write(string,*) 'local_id: ', local_id, 'Res: ', ii, Res(ii)
          option%io_buffer = ' NaN or INF of Residuals @ reactive_transport.F90: RTResidualFlux - BC flux ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo

#else
      call TFluxCoef_CD(option,cur_connection_set%area(iconn), &
                patch%boundary_velocities(:,sum_connection), &
                patch%boundary_tran_coefs(:,:,sum_connection)*vol_frac_prim, &
                0.5d0, & ! fraction upwind (0.d0 upwind, 0.5 central)
                T_11,T_12,T_21,T_22)
      call TFlux_CD(rt_parameter, &
                  rt_auxvars_bc(sum_connection), &
                  global_auxvars_bc(sum_connection), &
                  rt_auxvars(ghosted_id), &
                  global_auxvars(ghosted_id), &
                  T_11,T_12,T_21,T_22,option,Res_1,Res_2)

      iend = local_id*reaction%ncomp
      istart = iend-reaction%ncomp+1
      r_p(istart:iend)= r_p(istart:iend) + Res_2(1:reaction%ncomp)

      if (option%compute_mass_balance_new) then
      ! contribution to boundary 
        rt_auxvars_bc(sum_connection)%mass_balance_delta(:,iphase) = &
          rt_auxvars_bc(sum_connection)%mass_balance_delta(:,iphase) - Res_2
!        ! contribution to internal 
!        rt_auxvars(ghosted_id)%mass_balance_delta(:,iphase) = &
!          rt_auxvars(ghosted_id)%mass_balance_delta(:,iphase) + Res
        endif  
      
      if (associated(patch%boundary_tran_fluxes)) then
        patch%boundary_tran_fluxes(1:reaction%ncomp,sum_connection) = &
            -Res_2(1:reaction%ncomp)

      ! checking NaN or INF
      do ii=1,option%ncomp
        if(Res_1(ii) /= Res_1(ii) .or. &
          abs(Res_1(ii))>huge(Res_1(ii)) ) then
          write(string,*) 'local_id: ', local_id, 'Res: ', ii, Res_1(ii)
          option%io_buffer = ' NaN or INF of Residuals @ reactive_transport.F90: RTResidualFlux - BC Upwind ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo
      do ii=1,option%ncomp
        if(Res_2(ii) /= Res_2(ii) .or. &
          abs(Res_2(ii))>huge(Res_2(ii)) ) then
          write(string,*) 'local_id: ', local_id, 'Res: ', ii, Res_2(ii)
          option%io_buffer = ' NaN or INF of Residuals @ reactive_transport.F90: RTResidualFlux - BC Downwind ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo

#endif

      endif
    enddo
    boundary_condition => boundary_condition%next

  enddo


  ! Restore vectors
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
 
end subroutine RTResidualFlux

! ************************************************************************** !

subroutine RTResidualNonFlux(snes,xx,r,realization,ierr)
  ! 
  ! Computes the non-flux terms in the residual function for
  ! reactive transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Transport_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module
  use Debug_module
  use Logging_module
  !geh: please leave the "only" clauses for Secondary_Continuum_XXX as this
  !      resolves a bug in the Intel Visual Fortran compiler.
  use Secondary_Continuum_Aux_module, only : sec_transport_type
  use Secondary_Continuum_module, only : SecondaryRTResJacMulti
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(inout) :: r
  type(realization_subsurface_type) :: realization  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: r_p(:), accum_p(:), vec_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: iphase
  PetscInt :: i
  PetscInt :: istartaq, iendaq
  PetscInt :: istartcoll, iendcoll
  PetscInt :: istartall, iendall
  PetscInt :: idof
  PetscInt :: offset
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(reactive_transport_param_type), pointer :: rt_parameter
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:) 
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal :: Res(realization%reaction%ncomp)
  
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscReal :: qsrc(2), molality
  PetscInt :: flow_src_sink_type
  PetscReal :: scale, coef_in(2), coef_out(2)
  PetscReal :: Jup(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscBool :: volumetric
  PetscInt :: sum_connection
  PetscInt :: nphase

  ! CO2-specific
  PetscReal :: msrc(1:realization%option%nflowspec)
  PetscInt :: icomp, ieqgas

#ifdef CLM_PFLOTRAN
  ! temporarily changing option%iflag to pass 'ghosted_id' from CLM to PF RT bgc
  PetscInt :: option_iflag
#endif

  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscReal :: vol_frac_prim
  PetscReal :: sec_diffusion_coefficient
  PetscReal :: sec_porosity
  PetscReal :: res_sec_transport(realization%reaction%ncomp)

  PetscInt  :: ii
  character(len=MAXSTRINGLENGTH) :: string

  option => realization%option
  field => realization%field
  patch => realization%patch
  reaction => realization%reaction
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  rt_auxvars => patch%aux%RT%auxvars
  rt_auxvars_ss => patch%aux%RT%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  if (option%use_mc) then
    rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars
  endif
  nphase = rt_parameter%nphase
  
  ! Get pointer to Vector data
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
 
  vol_frac_prim = 1.d0

  if (.not.option%steady_state) then
#if 1
    call VecGetArrayF90(field%tran_accum, accum_p, ierr);CHKERRQ(ierr)
    r_p = r_p - accum_p / option%tran_dt
    call VecRestoreArrayF90(field%tran_accum, accum_p, ierr);CHKERRQ(ierr)
    ! Accumulation terms ------------------------------------
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle

      offset = (local_id-1)*reaction%ncomp
      istartall = offset + 1
      iendall = offset + reaction%ncomp

      call RTAccumulation(rt_auxvars(ghosted_id), &
                          global_auxvars(ghosted_id), &
                          material_auxvars(ghosted_id), &
                          reaction,option,Res)
      if (reaction%neqsorb > 0) then
        call RAccumulationSorb(rt_auxvars(ghosted_id), &
                               global_auxvars(ghosted_id), &
                               material_auxvars(ghosted_id), &
                               reaction,option,Res)
      endif
      Res = Res / option%tran_dt

      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(local_id)%epsilon
        Res = Res*vol_frac_prim
      endif        

      ! checking NaN or INF
      do ii=1,reaction%ncomp
        if(Res(ii) /= Res(ii) .or. abs(Res(ii))>huge(Res(ii)) ) then
          write(string,*) 'local_id: ', local_id, 'Res: ', ii, Res(ii)
          option%io_buffer = ' NaN or INF of Residuals @ reactive_transport.F90: RTResidualNonFlux - Accumulation of ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo
      
      r_p(istartall:iendall) = r_p(istartall:iendall) + Res(1:reaction%ncomp)
      
      ! Secondary continuum formation not implemented for Age equation
      if (reaction%calculate_water_age) then 
        call RAge(rt_auxvars(ghosted_id),global_auxvars(ghosted_id), &
                  material_auxvars(ghosted_id),option,reaction,Res)
        r_p(istartall:iendall) = r_p(istartall:iendall) + &
          Res(1:reaction%ncomp)
      endif
      if (reaction%calculate_tracer_age) then 
        call RAge(rt_auxvars(ghosted_id),global_auxvars(ghosted_id), &
                  material_auxvars(ghosted_id),option,reaction,Res)
        r_p(istartall:iendall) = r_p(istartall:iendall) + &
          Res(1:reaction%ncomp)
      endif
    enddo
  endif
#endif

#if 1

! ========== Secondary continuum transport source terms -- MULTICOMPONENT ======
  if (option%use_mc) then
  ! Secondary continuum contribution (SK 1/31/2013)
  ! only one secondary continuum for now for each primary continuum node
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      
      offset = (local_id-1)*reaction%ncomp
      istartall = offset + 1
      iendall = offset + reaction%ncomp
         
      sec_diffusion_coefficient = patch% &
                                  material_property_array(1)%ptr% &
                                  secondary_continuum_diff_coeff
      sec_porosity = patch%material_property_array(1)%ptr% &
                     secondary_continuum_porosity

      call SecondaryRTResJacMulti(rt_sec_transport_vars(local_id), &
                                  rt_auxvars(ghosted_id), &
                                  global_auxvars(ghosted_id), &
                                  material_auxvars(ghosted_id)%volume, &
                                  reaction, &
                                  sec_diffusion_coefficient, &
                                  sec_porosity, &
                                  option,res_sec_transport)

      ! checking NaN or INF
      do ii=1,reaction%ncomp
        if(res_sec_transport(ii) /= res_sec_transport(ii) .or. &
          abs(res_sec_transport(ii))>huge(res_sec_transport(ii)) ) then
          write(string,*) 'local_id: ', local_id, 'Res: ', ii, res_sec_transport(ii)
          option%io_buffer = ' NaN or INF of Residuals @ reactive_transport.F90: RTResidualNonFlux - Secondary continuum of ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo

      r_p(istartall:iendall) = r_p(istartall:iendall) - &
                               res_sec_transport(1:reaction%ncomp) ! in mol/s
                               
    enddo   
  endif
! ============== end secondary continuum coupling terms ========================

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set

    do iconn = 1, cur_connection_set%num_connections 
      sum_connection = sum_connection + 1     
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      offset = (local_id-1)*reaction%ncomp

      if (patch%imat(ghosted_id) <= 0) cycle
      
      istartaq = reaction%offset_aqueous + 1
      iendaq = reaction%offset_aqueous + reaction%naqcomp
      
      if (reaction%ncoll > 0) then
        istartcoll = reaction%offset_colloid + 1
        iendcoll = reaction%offset_colloid + reaction%ncoll
      endif

      if (associated(patch%ss_flow_vol_fluxes)) then
        qsrc(1:nphase) = patch%ss_flow_vol_fluxes(1:nphase,sum_connection)
      endif
      call TSrcSinkCoef(rt_parameter,qsrc,source_sink%tran_condition%itype, &
                        coef_in,coef_out)

      Res = 0.d0
      do iphase = 1, nphase
        Res(istartaq:iendaq) = &
          coef_in(iphase)*rt_auxvars(ghosted_id)%total(:,iphase) + &
          coef_out(iphase)*source_sink%tran_condition%cur_constraint_coupler% &
                                        rt_auxvar%total(:,iphase)
      enddo
      if (reaction%ncoll > 0) then
        Res(istartcoll:iendcoll) = coef_in*rt_auxvars(ghosted_id)%colloid%conc_mob(:) + &
                                   coef_out*source_sink%tran_condition%cur_constraint_coupler% &
                                              rt_auxvar%colloid%conc_mob(:)
      endif

      ! checking NaN or INF
      do ii=1,reaction%ncomp
        if(Res(ii) /= Res(ii) .or. &
          abs(Res(ii))>huge(Res(ii)) ) then
          write(string,*) 'local_id: ', local_id, 'Res: ', ii, Res(ii)
          option%io_buffer = ' NaN or INF of Residuals @ reactive_transport.F90: RTResidualNonFlux - SrcSink of ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo

      istartall = offset + 1
      iendall = offset + reaction%ncomp
      r_p(istartall:iendall) = r_p(istartall:iendall) + Res(1:reaction%ncomp)
      if (associated(patch%ss_tran_fluxes)) then
        patch%ss_tran_fluxes(:,sum_connection) = Res(:)
      endif
      if (option%compute_mass_balance_new) then
        ! contribution to boundary 
        iphase = 1
        rt_auxvars_ss(sum_connection)%mass_balance_delta(:,iphase) = &
          rt_auxvars_ss(sum_connection)%mass_balance_delta(:,iphase) + Res
        ! contribution to internal 
!        rt_auxvars(ghosted_id)%mass_balance_delta(:,iphase) = &
!          rt_auxvars(ghosted_id)%mass_balance_delta(:,iphase) - Res
        endif
    enddo
    source_sink => source_sink%next
  enddo

#endif

#if 1  
! Reactions
  if (associated(reaction)) then
  
    call PetscLogEventBegin(logging%event_rt_res_reaction,ierr);CHKERRQ(ierr)
    
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      offset = (local_id-1)*reaction%ncomp
      istartall = offset + 1
      iendall = offset + reaction%ncomp
      Res = 0.d0
      Jup = 0.d0
      if (.not.option%use_isothermal) then
        call RUpdateTempDependentCoefs(global_auxvars(ghosted_id),reaction, &
                                       PETSC_FALSE,option)
      endif      

!F.-M. YUAN: option%iflag IS used here as indexing of cell-id for passing data from
! clm_pf_idata%??? to PFLOTRAN for driving reaction sandboxes
! note: 'local_id' is used in those sandboxes, but after checking when in parallel mode,
! it should be 'ghosted_id', because in 'clm_pf_idata%???', those are defined as PETSC seq. vecs.
#ifdef CLM_PFLOTRAN
    option_iflag = option%iflag
    option%iflag = ghosted_id
#endif

      call RReaction(Res,Jup,PETSC_FALSE,rt_auxvars(ghosted_id), &
                     global_auxvars(ghosted_id), &
                     material_auxvars(ghosted_id), &
                     reaction,option)

#ifdef CLM_PFLOTRAN
    option%iflag = option_iflag
#endif

      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(local_id)%epsilon
        Res = Res*vol_frac_prim
      endif 

      ! checking NaN or INF
      do ii=1,reaction%ncomp
        if(Res(ii) /= Res(ii) .or. &
          abs(Res(ii))>huge(Res(ii)) ) then
          write(string,*) 'local_id: ', local_id, 'Res: ', ii, Res(ii)
          option%io_buffer = ' NaN or INF of Residuals @ reactive_transport.F90: RTResidualNonFlux - Reaction of ' // &
            trim(string)
          call printErrMsg(option)
        endif

      enddo

      r_p(istartall:iendall) = r_p(istartall:iendall) + Res(1:reaction%ncomp)                    

    enddo

    call PetscLogEventEnd(logging%event_rt_res_reaction,ierr);CHKERRQ(ierr)
  endif
#endif

  if (patch%aux%RT%inactive_cells_exist) then
    do i=1,patch%aux%RT%n_zero_rows
      r_p(patch%aux%RT%zero_rows_local(i)) = 0.d0
    enddo
  endif

  ! Restore vectors
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
 
  ! Mass Transfer
  if (field%tran_mass_transfer /= PETSC_NULL_VEC) then
    ! scale by -1.d0 for contribution to residual.  A negative contribution
    ! indicates mass being added to system.
    call VecGetArrayF90(field%tran_mass_transfer,vec_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%tran_mass_transfer,vec_p,ierr);CHKERRQ(ierr)
    call VecAXPY(r,-1.d0,field%tran_mass_transfer,ierr);CHKERRQ(ierr)
  endif

end subroutine RTResidualNonFlux

! ************************************************************************** !

subroutine RTJacobian(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  use Logging_module
  use Debug_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer  
  type(grid_type),  pointer :: grid
  character(len=MAXSTRINGLENGTH) :: string

  call PetscLogEventBegin(logging%event_rt_jacobian,ierr);CHKERRQ(ierr)

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif
    
  call MatZeroEntries(J,ierr);CHKERRQ(ierr)


  call PetscLogEventBegin(logging%event_rt_jacobian1,ierr);CHKERRQ(ierr)


  ! pass #1 for internal and boundary flux terms  
  call RTJacobianFlux(snes,xx,J,J,realization,ierr)

  call PetscLogEventEnd(logging%event_rt_jacobian1,ierr);CHKERRQ(ierr)
  call PetscLogEventBegin(logging%event_rt_jacobian2,ierr);CHKERRQ(ierr)
  
  ! pass #2 for everything else
  call RTJacobianNonFlux(snes,xx,J,J,realization,ierr)

  call PetscLogEventEnd(logging%event_rt_jacobian2,ierr);CHKERRQ(ierr)
    
  if (realization%debug%matview_Jacobian) then
    string = 'RTjacobian'
    call DebugCreateViewer(realization%debug,string,realization%option,viewer)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  if (realization%reaction%use_log_formulation) then
    call MatDiagonalScaleLocal(J,realization%field%tran_work_loc, &
                               ierr);CHKERRQ(ierr)

    if (realization%debug%matview_Jacobian) then
      string = 'RTjacobianLog'
      call DebugCreateViewer(realization%debug,string,realization%option,viewer)
      call MatView(J,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif
    
  endif

  call PetscLogEventEnd(logging%event_rt_jacobian,ierr);CHKERRQ(ierr)
  
end subroutine RTJacobian

! ************************************************************************** !

subroutine RTJacobianFlux(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the flux term entries in the Jacobian for
  ! reactive transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Transport_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  use Logging_module  
  use Secondary_Continuum_Aux_module
  
  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: r_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart, iend                        
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reactive_transport_param_type), pointer :: rt_parameter
  type(reaction_type), pointer :: reaction
      
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:), rt_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn, rdum
  
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscReal :: vol_frac_prim

#ifdef CENTRAL_DIFFERENCE
  PetscReal :: T_11(realization%option%transport%nphase)
  PetscReal :: T_12(realization%option%transport%nphase)
  PetscReal :: T_21(realization%option%transport%nphase)
  PetscReal :: T_22(realization%option%transport%nphase)
  PetscReal :: J_11(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: J_12(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: J_21(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: J_22(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Res(realization%reaction%ncomp)  
#else
  PetscReal :: coef_up(realization%patch%aux%RT%rt_parameter%naqcomp, &
                       realization%option%transport%nphase)
  PetscReal :: coef_dn(realization%patch%aux%RT%rt_parameter%naqcomp, &
                       realization%option%transport%nphase)
  PetscReal :: Jup(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Jdn(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Res(realization%reaction%ncomp)  
#endif

  ! zero out xy-direction transport
  PetscReal :: unitvec_xyz(3)

  option => realization%option
  field => realization%field
  patch => realization%patch  
  grid => patch%grid
  reaction => realization%reaction
  rt_parameter => patch%aux%RT%rt_parameter
  rt_auxvars => patch%aux%RT%auxvars
  rt_auxvars_bc => patch%aux%RT%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  if (option%use_mc) then
    rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars
  endif


  vol_frac_prim = 1.d0

  ! Interior Flux Terms -----------------------------------
  unitvec_xyz(1:2) = 0.d0
  unitvec_xyz(3)   = 1.d0

  ! must zero out Jacobian blocks

  call PetscLogEventBegin(logging%event_rt_jacobian_flux,ierr);CHKERRQ(ierr)

  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle

      !----------------
      if ( option%flow%only_vertical_flow .or. option%transport%only_vertical_tran ) then
         if(abs(dot_product(cur_connection_set%dist(1:3,iconn),unitvec_xyz)) < 1.d-20) cycle
      end if
      !----------------

      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(local_id_up)%epsilon
      endif 

#ifndef CENTRAL_DIFFERENCE
      call TFluxCoef(rt_parameter,option,cur_connection_set%area(iconn), &
                patch%internal_velocities(:,sum_connection), &
                patch%internal_tran_coefs(:,:,sum_connection)*vol_frac_prim, &
                cur_connection_set%dist(-1,iconn), &
                coef_up,coef_dn)
      call TFluxDerivative(rt_parameter, &
                           rt_auxvars(ghosted_id_up), &
                           global_auxvars(ghosted_id_up), &
                           rt_auxvars(ghosted_id_dn), &
                           global_auxvars(ghosted_id_dn), &
                           coef_up,coef_dn,option,Jup,Jdn)
      if (local_id_up>0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
   
      if (local_id_dn>0) then
        Jup = -Jup
        Jdn = -Jdn
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif

#else
      call TFluxCoef_CD(option,cur_connection_set%area(iconn), &
                patch%internal_velocities(:,sum_connection), &
                patch%internal_tran_coefs(:,:,sum_connection)*vol_frac_prim, &
                cur_connection_set%dist(-1,iconn), &
                T_11,T_12,T_21,T_22)
      call TFluxDerivative_CD(rt_parameter, &
                           rt_auxvars(ghosted_id_up), &
                           global_auxvars(ghosted_id_up), &
                           rt_auxvars(ghosted_id_dn), &
                           global_auxvars(ghosted_id_dn), &
                           T_11,T_12,T_21,T_22,option, &
                           J_11,J_12,J_21,J_22)
      if (local_id_up>0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      J_11,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      J_12,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
   
      if (local_id_dn>0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      J_22,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      J_21,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
#endif


    enddo
    cur_connection_set => cur_connection_set%next
  enddo    

  call PetscLogEventEnd(logging%event_rt_jacobian_flux,ierr);CHKERRQ(ierr)
  
  ! Boundary Flux Terms -----------------------------------
  ! must zero out Jacobian block

  call PetscLogEventBegin(logging%event_rt_jacobian_fluxbc,ierr);CHKERRQ(ierr)

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle
    
      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(local_id)%epsilon
      endif 

#ifndef CENTRAL_DIFFERENCE
      ! TFluxCoef accomplishes the same as what TBCCoef would
      call TFluxCoef(rt_parameter,option,cur_connection_set%area(iconn), &
                patch%boundary_velocities(:,sum_connection), &
                patch%boundary_tran_coefs(:,:,sum_connection)*vol_frac_prim, &
                0.5d0, & ! fraction upwind (0.d0 upwind, 0.5 central)
                coef_up,coef_dn)
      ! TFluxDerivative accomplishes the same as what TBCFluxDerivative would
      call TFluxDerivative(rt_parameter, &
                           rt_auxvars_bc(sum_connection), &
                           global_auxvars_bc(sum_connection), &
                           rt_auxvars(ghosted_id), &
                           global_auxvars(ghosted_id), &
                           coef_up,coef_dn,option,Jup,Jdn)

      !Jup not needed 
      Jdn = -Jdn
      
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn,ADD_VALUES, &
                                    ierr);CHKERRQ(ierr)
 
#else
      call TFluxCoef_CD(option,cur_connection_set%area(iconn), &
                 patch%boundary_velocities(:,sum_connection), &
                 patch%boundary_tran_coefs(:,:,sum_connection)*vol_frac_prim, &
                 0.5d0, & ! fraction upwind (0.d0 upwind, 0.5 central)
                 T_11,T_12,T_21,T_22)
      call TFluxDerivative_CD(rt_parameter, &
                           rt_auxvars_bc(sum_connection), &
                           global_auxvars_bc(sum_connection), &
                           rt_auxvars(ghosted_id), &
                           global_auxvars(ghosted_id), &
                           T_11,T_12,T_21,T_22,option, &
                           J_11,J_12,J_21,J_22)
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,J_22,ADD_VALUES, &
                                    ierr);CHKERRQ(ierr)
#endif
 
    enddo
    boundary_condition => boundary_condition%next
  enddo
  call PetscLogEventEnd(logging%event_rt_jacobian_fluxbc,ierr);CHKERRQ(ierr)

end subroutine RTJacobianFlux

! ************************************************************************** !

subroutine RTJacobianNonFlux(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes non-flux term entries in the Jacobian for
  ! reactive transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Transport_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  use Logging_module
  use Secondary_Continuum_Aux_module

  
  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: work_loc_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: istartaq, iendaq
  PetscInt :: istart, iend
  PetscInt :: offset, idof                  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(reactive_transport_param_type), pointer :: rt_parameter
  PetscInt :: tran_pc
    
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:), rt_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal :: Jup(realization%reaction%ncomp,realization%reaction%ncomp)
  PetscReal :: Res(realization%reaction%ncomp)    
  
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection
  PetscReal :: qsrc(2), rdum
  PetscBool :: volumetric
  PetscInt :: flow_src_sink_type
  PetscReal :: coef_in(2), coef_out(2)
  PetscReal :: scale
  
  ! secondary continuum variables
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  PetscReal :: vol_frac_prim
  PetscReal :: sec_diffusion_coefficient
  PetscReal :: sec_porosity
  PetscReal :: jac_transport(realization%reaction%naqcomp,realization%reaction%naqcomp)
  PetscInt :: ncomp
  PetscInt :: nphase
  PetscInt :: iphase

#ifdef CLM_PFLOTRAN
  ! temporarily changing option%iflag to pass 'ghosted_id' from CLM to PF RT bgc
  PetscInt :: option_iflag
#endif
  
  option => realization%option
  field => realization%field
  patch => realization%patch  
  reaction => realization%reaction
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  rt_auxvars => patch%aux%RT%auxvars
  rt_auxvars_bc => patch%aux%RT%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
  if (option%use_mc) then
    rt_sec_transport_vars => patch%aux%SC_RT%sec_transport_vars
  endif
  nphase = rt_parameter%nphase

  vol_frac_prim = 1.d0
  
  if (.not.option%steady_state) then
  call PetscLogEventBegin(logging%event_rt_jacobian_accum,ierr);CHKERRQ(ierr)
#if 1  
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      
      call RTAccumulationDerivative(rt_auxvars(ghosted_id), &
                                    global_auxvars(ghosted_id), &
                                    material_auxvars(ghosted_id), &
                                    reaction,option,Jup) 
                                    
      if (reaction%neqsorb > 0) then
        call RAccumulationSorbDerivative(rt_auxvars(ghosted_id), &
                                         global_auxvars(ghosted_id), &
                                         material_auxvars(ghosted_id), &
                                         reaction,option,Jup)
      endif
      
      if (option%use_mc) then
      
        vol_frac_prim = rt_sec_transport_vars(local_id)%epsilon
        Jup = Jup*vol_frac_prim

        sec_diffusion_coefficient = patch%material_property_array(1)% &
                                    ptr%secondary_continuum_diff_coeff
        sec_porosity = patch%material_property_array(1)%ptr% &
                       secondary_continuum_porosity
                        
        if (realization%reaction%ncomp /= realization%reaction%naqcomp) then
          option%io_buffer = 'Current multicomponent implementation is for '// &
                             'aqueous reactions only'
          call printErrMsg(option)
        endif
        
        if (rt_sec_transport_vars(local_id)%sec_jac_update) then
          jac_transport = rt_sec_transport_vars(local_id)%sec_jac
        else
          option%io_buffer = 'RT secondary continuum term in primary '// &
                             'jacobian not updated'
          call printErrMsg(option)
        endif
         
        Jup = Jup - jac_transport                                                                   
                                                                                
      endif

      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup,ADD_VALUES, &
                                    ierr);CHKERRQ(ierr)
    enddo
#endif
  call PetscLogEventEnd(logging%event_rt_jacobian_accum,ierr);CHKERRQ(ierr)
  endif
#if 1
  ! Source/Sink terms -------------------------------------
  call PetscLogEventBegin(logging%event_rt_jacobian_ss,ierr);CHKERRQ(ierr)
  source_sink => patch%source_sink_list%first 
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    if (reaction%ncoll > 0) then
      option%io_buffer = 'Source/sink not yet implemented for colloids'
      call printErrMsg(option)
    endif

    cur_connection_set => source_sink%connection_set

    do iconn = 1, cur_connection_set%num_connections      
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      istartaq = reaction%offset_aqueous + 1
      iendaq = reaction%offset_aqueous + reaction%naqcomp
      
      if (associated(patch%ss_flow_vol_fluxes)) then
        qsrc(1:nphase) = patch%ss_flow_vol_fluxes(1:nphase,sum_connection)
      endif
      call TSrcSinkCoef(rt_parameter,qsrc,source_sink%tran_condition%itype, &
                        coef_in,coef_out)
      Jup = 0.d0
      ! coef_in is non-zero
      do iphase = 1, nphase
        Jup(istartaq:iendaq,istartaq:iendaq) = &
          Jup(istartaq:iendaq,istartaq:iendaq) + coef_in(iphase)* &
            rt_auxvars(ghosted_id)%aqueous%dtotal(:,:,iphase)
      enddo                       
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
    enddo
    source_sink => source_sink%next
  enddo
  
  call PetscLogEventEnd(logging%event_rt_jacobian_ss,ierr);CHKERRQ(ierr)
#endif


#if 1  
! Reactions
  if (associated(reaction)) then

    call PetscLogEventBegin(logging%event_rt_jac_reaction,ierr);CHKERRQ(ierr)
                              
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      Res = 0.d0
      Jup = 0.d0
      if (.not.option%use_isothermal) then
        call RUpdateTempDependentCoefs(global_auxvars(ghosted_id),reaction, &
                                       PETSC_FALSE,option)
      endif      

!F.-M. YUAN: option%iflag IS used here as indexing of cell-id for passing data from
! clm_pf_idata%??? to PFLOTRAN for driving reaction sandboxes
! note: 'local_id' is used in those sandboxes, but after checking when in parallel mode,
! it should be 'ghosted_id', because in 'clm_pf_idata%???', those are defined as PETSC seq. vecs.
#ifdef CLM_PFLOTRAN
    option_iflag = option%iflag
    option%iflag = ghosted_id
#endif

      call RReactionDerivative(Res,Jup,rt_auxvars(ghosted_id), &
                               global_auxvars(ghosted_id), &
                               material_auxvars(ghosted_id), &
                               reaction,option)

#ifdef CLM_PFLOTRAN
    option%iflag = option_iflag
#endif

      if (option%use_mc) then
        vol_frac_prim = rt_sec_transport_vars(local_id)%epsilon
        Jup = Jup*vol_frac_prim
      endif
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1, &
                                    Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
    enddo
    
    call PetscLogEventEnd(logging%event_rt_jac_reaction,ierr);CHKERRQ(ierr)
    
  endif
#endif
  
  ! Mass Transfer - since the current implementation of mass transfer has
  ! mass transfer being fixed.  Nothing to do here as the contribution to
  ! the derivatives is zero.
!  if (field%tran_mass_transfer /= 0) then
!  endif
 
  if (reaction%use_log_formulation) then
    call PetscLogEventBegin(logging%event_rt_jacobian_zero_calc, &
                            ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%tran_work_loc, work_loc_p, ierr);CHKERRQ(ierr)
    do ghosted_id = 1, grid%ngmax  ! For each local node do...
      offset = (ghosted_id-1)*reaction%ncomp
      if (patch%imat(ghosted_id) <= 0) then
        istart = offset + 1
        iend = offset + reaction%ncomp
        work_loc_p(istart:iend) = 1.d0
      else
        istartaq = offset + reaction%offset_aqueous + 1
        iendaq = offset + reaction%offset_aqueous + reaction%naqcomp
        work_loc_p(istartaq:iendaq) = rt_auxvars(ghosted_id)%pri_molal(:)
        if (reaction%immobile%nimmobile > 0) then
          istart = offset + reaction%offset_immobile + 1
          iend = offset + reaction%offset_immobile + reaction%immobile%nimmobile
          work_loc_p(istart:iend) = &
            rt_auxvars(ghosted_id)%immobile(:)
        endif
        if (reaction%ncoll > 0) then
          istart = offset + reaction%offset_colloid + 1
          iend = offset + reaction%offset_colloid + reaction%ncoll
          work_loc_p(istart:iend) = &
            rt_auxvars(ghosted_id)%colloid%conc_mob(:)
        endif
      endif
    enddo
    call VecRestoreArrayF90(field%tran_work_loc, work_loc_p,  &
                            ierr);CHKERRQ(ierr)
    call PetscLogEventEnd(logging%event_rt_jacobian_zero_calc, &
                          ierr);CHKERRQ(ierr)
  endif

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  
  if (patch%aux%RT%inactive_cells_exist) then
    call PetscLogEventBegin(logging%event_rt_jacobian_zero,ierr);CHKERRQ(ierr)
    rdum = 1.d0
    call MatZeroRowsLocal(A,patch%aux%RT%n_zero_rows, &
                          patch%aux%RT%zero_rows_local_ghosted,rdum, &
                          PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
    call PetscLogEventEnd(logging%event_rt_jacobian_zero,ierr);CHKERRQ(ierr)
  endif

end subroutine RTJacobianNonFlux

! ************************************************************************** !

subroutine RTJacobianEquilibrateCO2(J,realization)
  ! 
  ! Adds CO2 saturation constraint to Jacobian for
  ! reactive transport
  ! 
  ! Author: Glenn Hammond/Peter Lichtner
  ! Date: 12/12/14
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module

  implicit none

  Mat :: J
  type(realization_subsurface_type) :: realization  
  
  PetscInt :: local_id, ghosted_id
  PetscInt :: idof                  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  PetscErrorCode :: ierr
    
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscInt :: zero_rows(realization%patch%grid%nlmax * realization%option%ntrandof)
  PetscInt :: ghosted_rows(realization%patch%grid%nlmax)
  PetscInt :: zero_count
  PetscInt :: i, jco2
  PetscReal :: jacobian_entry
  PetscReal :: eps = 1.d-6

  option => realization%option
  field => realization%field
  patch => realization%patch  
  reaction => realization%reaction
  grid => patch%grid
  rt_auxvars => patch%aux%RT%auxvars
  global_auxvars => patch%aux%Global%auxvars

  ! loop over cells twice.  the first time to zero (all rows to be zeroed have
  ! to be zeroed in a single call by passing in a list).  the second loop to 
  ! add the equilibration

  jacobian_entry = 1.d0
  jco2 = reaction%species_idx%co2_aq_id
  zero_count = 0
  zero_rows = 0
  ghosted_rows = 0
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (global_auxvars(ghosted_id)%sat(GAS_PHASE) > eps .and. &
      global_auxvars(ghosted_id)%sat(GAS_PHASE) < 1.d0-eps) then
      zero_count = zero_count + 1
      zero_rows(zero_count) = jco2+(ghosted_id-1)*reaction%ncomp-1
      ghosted_rows(zero_count) = ghosted_id
    endif
  enddo

  call MatZeroRowsLocal(J,zero_count,zero_rows(1:zero_count),jacobian_entry, &
                        PETSC_NULL_VEC,PETSC_NULL_VEC, &
                        ierr);CHKERRQ(ierr)

  do i = 1, zero_count
    ghosted_id = ghosted_rows(i) ! zero indexing back to 1-based
    if (patch%imat(ghosted_id) <= 0) cycle
    if (reaction%use_log_formulation) then
      jacobian_entry = rt_auxvars(ghosted_id)%pri_molal(jco2)
    else
      jacobian_entry = 1.d0
    endif

    idof = (ghosted_id-1)*option%ntrandof + jco2
    call MatSetValuesLocal(J,1,idof-1,1,idof-1,jacobian_entry,INSERT_VALUES, &
                           ierr);CHKERRQ(ierr)
  enddo

  call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

end subroutine RTJacobianEquilibrateCO2

! ************************************************************************** !

subroutine RTUpdateActivityCoefficients(realization,update_cells,update_bcs)
  ! 
  ! Updates activity coeffficients for cell and boundary auxvars
  ! reactive transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/06/16
  ! 
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Option_module
  
  implicit none

  type(realization_subsurface_type) :: realization
  PetscBool :: update_bcs
  PetscBool :: update_cells
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: ghosted_id, local_id, sum_connection, iconn
  
  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  reaction => realization%reaction

  if (update_cells) then
    do ghosted_id = 1, grid%ngmax
      if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      call RActivityCoefficients(patch%aux%RT%auxvars(ghosted_id), &
                                 patch%aux%Global%auxvars(ghosted_id), &
                                 reaction,option)
    enddo
  endif

  if (update_bcs) then
    boundary_condition => patch%boundary_condition_list%first
    sum_connection = 0    
    do 
      if (.not.associated(boundary_condition)) exit
      cur_connection_set => boundary_condition%connection_set
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        
        if (patch%imat(ghosted_id) <= 0) cycle

        call RActivityCoefficients(patch%aux%RT%auxvars_bc(sum_connection), &
                                   patch%aux%Global% &
                                     auxvars_bc(sum_connection), &
                                   reaction,option)
      enddo ! iconn
      boundary_condition => boundary_condition%next
    enddo
  endif 

end subroutine RTUpdateActivityCoefficients

! ************************************************************************** !

subroutine RTUpdateAuxVars(realization,update_cells,update_bcs, &
                           update_activity_coefs)
  ! 
  ! Updates the auxiliary variables associated with
  ! reactive transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/15/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Option_module
  use Field_module
  use Logging_module

  implicit none

  type(realization_subsurface_type) :: realization
  PetscBool :: update_bcs
  PetscBool :: update_cells
  PetscBool :: update_activity_coefs
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn
  PetscInt :: istartaq, iendaq 
  PetscInt :: istartcoll, iendcoll
  PetscInt :: istartaq_loc, iendaq_loc
  PetscInt :: istartcoll_loc, iendcoll_loc
  PetscInt :: istartim, iendim
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%reaction%ncomp)
  PetscReal, pointer :: basis_molarity_p(:)
  PetscReal, pointer :: basis_coll_conc_p(:)
  PetscReal :: weight
  PetscInt, parameter :: iphase = 1
  PetscInt :: offset
  PetscErrorCode :: ierr
  PetscBool :: skip_equilibrate_constraint
  PetscInt, save :: icall
  
  data icall/0/

  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  field => realization%field
  reaction => realization%reaction


  call VecGetArrayReadF90(field%tran_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  if (update_cells) then

    call PetscLogEventBegin(logging%event_rt_auxvars,ierr);CHKERRQ(ierr)
  
    do ghosted_id = 1, grid%ngmax
      if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
      !geh - Ignore inactive cells with inactive materials

      if (patch%imat(ghosted_id) <= 0) cycle

      offset = (ghosted_id-1)*reaction%ncomp
      istartaq = offset + reaction%offset_aqueous + 1
      iendaq = offset + reaction%offset_aqueous + reaction%naqcomp
      
      patch%aux%RT%auxvars(ghosted_id)%pri_molal = xx_loc_p(istartaq:iendaq)
      if (reaction%immobile%nimmobile > 0) then
        istartim = offset + reaction%offset_immobile + 1
        iendim = offset + reaction%offset_immobile + reaction%immobile%nimmobile
        patch%aux%RT%auxvars(ghosted_id)%immobile = xx_loc_p(istartim:iendim)
      endif
      if (reaction%ncoll > 0) then
        istartcoll = offset + reaction%offset_colloid + 1
        iendcoll = offset + reaction%offset_colloid + reaction%ncoll
        patch%aux%RT%auxvars(ghosted_id)%colloid%conc_mob = &
          xx_loc_p(istartcoll:iendcoll)* &
          patch%aux%Global%auxvars(ghosted_id)%den_kg(1)*1.d-3
      endif
      if (.not.option%use_isothermal) then
        call RUpdateTempDependentCoefs(patch%aux%Global%auxvars(ghosted_id), &
                                       reaction,PETSC_FALSE, &
                                       option)
      endif
      if (update_activity_coefs) then
        call RActivityCoefficients(patch%aux%RT%auxvars(ghosted_id), &
                                   patch%aux%Global%auxvars(ghosted_id), &
                                   reaction,option)
      endif
      call RTAuxVarCompute(patch%aux%RT%auxvars(ghosted_id), &
                           patch%aux%Global%auxvars(ghosted_id), &
                           patch%aux%Material%auxvars(ghosted_id), &
                           reaction,option)
    enddo

    call PetscLogEventEnd(logging%event_rt_auxvars,ierr);CHKERRQ(ierr)
  endif

  if (update_bcs) then

    call PetscLogEventBegin(logging%event_rt_auxvars_bc,ierr);CHKERRQ(ierr)

    boundary_condition => patch%boundary_condition_list%first
    sum_connection = 0    
    do 
      if (.not.associated(boundary_condition)) exit
      cur_connection_set => boundary_condition%connection_set

      basis_molarity_p => boundary_condition%tran_condition% &
        cur_constraint_coupler%aqueous_species%basis_molarity

      if (reaction%ncoll > 0) then
        basis_coll_conc_p => boundary_condition%tran_condition% &
                             cur_constraint_coupler%colloids%basis_conc_mob
      endif

      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        
        if (patch%imat(ghosted_id) <= 0) cycle

        offset = (ghosted_id-1)*reaction%ncomp
        istartaq_loc = reaction%offset_aqueous + 1
        iendaq_loc = reaction%offset_aqueous + reaction%naqcomp
        istartaq = offset + istartaq_loc
        iendaq = offset + iendaq_loc
    
        if (reaction%ncoll > 0) then
          istartcoll_loc = reaction%offset_colloid + 1
          iendcoll_loc = reaction%offset_colloid + reaction%ncoll
          istartcoll = offset + istartcoll_loc
          iendcoll = offset + iendcoll_loc
        endif

          skip_equilibrate_constraint = PETSC_FALSE
        ! Chuan needs to fill this in.
          select case(boundary_condition%tran_condition%itype)
            case(CONCENTRATION_SS,DIRICHLET_BC,NEUMANN_BC)
              ! don't need to do anything as the constraint below provides all
              ! the concentrations, etc.
              
              !geh: terrible kludge, but should work for now.
              !geh: the problem is that ...%pri_molal() on first call is 
              !      zero and PETSC_TRUE is passed into 
              !      ReactionEquilibrateConstraint() below for 
              !      use_prev_soln_as_guess.  If the previous solution is 
              !      zero, the code will crash.
              if (patch%aux%RT%auxvars_bc(sum_connection)%pri_molal(1) < &
                  1.d-200) then
!               patch%aux%RT%auxvars_bc(sum_connection)%pri_molal = 1.d-9
                patch%aux%RT%auxvars_bc(sum_connection)%pri_molal = &
                    xx_loc_p(istartaq:iendaq)
              endif
            case(DIRICHLET_ZERO_GRADIENT_BC)
              if (patch%boundary_velocities(iphase,sum_connection) >= 0.d0) then
                  ! don't need to do anything as the constraint below 
                  ! provides all the concentrations, etc.
                  
                if (patch%aux%RT%auxvars_bc(sum_connection)%pri_molal(1) < &
                    1.d-200) then
!                 patch%aux%RT%auxvars_bc(sum_connection)%pri_molal = 1.d-9
                  patch%aux%RT%auxvars_bc(sum_connection)%pri_molal = &
                    xx_loc_p(istartaq:iendaq)
                endif
              else
                ! same as zero_gradient below
                skip_equilibrate_constraint = PETSC_TRUE
                patch%aux%RT%auxvars_bc(sum_connection)%pri_molal = &
                  xx_loc_p(istartaq:iendaq)
                if (reaction%ncoll > 0) then
                  patch%aux%RT%auxvars_bc(sum_connection)%colloid%conc_mob = &
                    xx_loc_p(istartcoll:iendcoll)* &
                    patch%aux%Global%auxvars_bc(sum_connection)%den_kg(1)*1.d-3
                endif
              endif
            case(ZERO_GRADIENT_BC)
              skip_equilibrate_constraint = PETSC_TRUE
              patch%aux%RT%auxvars_bc(sum_connection)%pri_molal = &
                xx_loc_p(istartaq:iendaq)
              if (reaction%ncoll > 0) then
                patch%aux%RT%auxvars_bc(sum_connection)%colloid%conc_mob = &
                  xx_loc_p(istartcoll:iendcoll)* &
                  patch%aux%Global%auxvars_bc(sum_connection)%den_kg(1)*1.d-3
              endif                
          end select
          ! no need to update boundary fluid density since it is already set
          if (.not.skip_equilibrate_constraint) then
            ! print *,'RT redo constrain on BCs: 1: ', sum_connection
            call ReactionEquilibrateConstraint( &
              patch%aux%RT%auxvars_bc(sum_connection), &
              patch%aux%Global%auxvars_bc(sum_connection), &
              patch%aux%Material%auxvars(ghosted_id),reaction, &
              boundary_condition%tran_condition%cur_constraint_coupler% &
                constraint_name, &
              boundary_condition%tran_condition%cur_constraint_coupler% &
                aqueous_species, &
              boundary_condition%tran_condition%cur_constraint_coupler% &
                free_ion_guess, &
              boundary_condition%tran_condition%cur_constraint_coupler% &
                minerals, &
              boundary_condition%tran_condition%cur_constraint_coupler% &
                colloids, &
              boundary_condition%tran_condition%cur_constraint_coupler% &
                immobile_species, &
              boundary_condition%tran_condition%cur_constraint_coupler% &
                num_iterations, &
              PETSC_TRUE,option)
            ! print *,'RT redo constrain on BCs: 2: ', sum_connection  
          endif         


      enddo ! iconn
      boundary_condition => boundary_condition%next
    enddo

    call PetscLogEventEnd(logging%event_rt_auxvars_bc,ierr);CHKERRQ(ierr)

  endif 

  call VecRestoreArrayReadF90(field%tran_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  icall = icall+ 1
  
end subroutine RTUpdateAuxVars

! ************************************************************************** !

subroutine RTMaxChange(realization,dcmax,dvfmax)
  ! 
  ! Computes the maximum change in the solution vector
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/15/08
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Patch_module
  use Grid_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  PetscReal :: dcmax
  PetscReal :: dvfmax
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field 
  type(reaction_type), pointer :: reaction
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscReal, pointer :: dxx_ptr(:), xx_ptr(:), yy_ptr(:)
  PetscInt :: local_id, ghosted_id, imnrl
  PetscReal :: delta_volfrac
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  reaction => realization%reaction
  patch => realization%patch
  grid => patch%grid
  rt_auxvars => patch%aux%RT%auxvars  

  dcmax = 0.d0
  dvfmax = 0.d0
  
  call VecWAXPY(field%tran_dxx,-1.d0,field%tran_xx,field%tran_yy, &
                ierr);CHKERRQ(ierr)
  
  call VecStrideNorm(field%tran_dxx,ZERO_INTEGER,NORM_INFINITY,dcmax, &
                     ierr);CHKERRQ(ierr)
                     
#if 1
  ! update mineral volume fractions
  if (reaction%mineral%nkinmnrl > 0) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      do imnrl = 1, reaction%mineral%nkinmnrl
        delta_volfrac = rt_auxvars(ghosted_id)%mnrl_rate(imnrl)* &
                        reaction%mineral%kinmnrl_molar_vol(imnrl)* &
                        option%tran_dt
        dvfmax = max(dabs(delta_volfrac),dvfmax)
      enddo
    enddo
  endif 
#endif
      
end subroutine RTMaxChange


! ************************************************************************** !

subroutine RTExplicitAdvection(realization)
  ! 
  ! Updates advective transport explicitly
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/12
  ! 

  use Realization_Subsurface_class

  use Discretization_module
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  
  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(reaction_type), pointer :: reaction
  type(discretization_type), pointer :: discretization
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  type(reactive_transport_param_type), pointer :: rt_parameter
  class(material_auxvar_type), pointer :: material_auxvars(:)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscInt :: iphase
  PetscInt :: id_up2, id_dn2
  PetscInt :: local_start, local_end, istart, iend
  PetscInt :: ntvddof
  PetscReal :: qsrc(2), coef_in(2), coef_out(2)
  PetscReal :: velocity, area, psv_t  
  PetscReal :: flux(realization%reaction%ncomp)
  PetscInt :: nphase
  
  PetscReal :: sum_flux(realization%reaction%ncomp,realization%patch%grid%ngmax)
  
  PetscReal, pointer :: tran_xx_p(:)
  PetscReal, pointer :: tvd_ghosts_p(:)
  PetscReal, pointer :: rhs_coef_p(:)
  PetscReal, pointer :: total_up2(:,:), total_dn2(:,:)
  PetscErrorCode :: ierr
  PetscViewer :: viewer

  ! zero out xy-direction transport
  PetscReal :: unitvec_xyz(3)

  procedure (TFluxLimiterDummy), pointer :: TFluxLimitPtr
  
  select case(realization%option%transport%tvd_flux_limiter)
    case(TVD_LIMITER_UPWIND)
      TFluxLimitPtr => TFluxLimitUpwind
    case(TVD_LIMITER_MC)
      TFluxLimitPtr => TFluxLimitMC
    case(TVD_LIMITER_MINMOD)
      TFluxLimitPtr => TFluxLimitMinmod
    case(TVD_LIMITER_SUPERBEE)
      TFluxLimitPtr => TFluxLimitSuperBee
    case(TVD_LIMITER_VAN_LEER)
      TFluxLimitPtr => TFluxLimitVanLeer
    case default
      TFluxLimitPtr => TFluxLimiter
  end select

  option => realization%option
  field => realization%field
  patch => realization%patch
  discretization => realization%discretization
  reaction => realization%reaction
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  rt_auxvars => patch%aux%RT%auxvars
  rt_auxvars_bc => patch%aux%RT%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
  
  ntvddof = patch%aux%RT%rt_parameter%naqcomp
  nphase = patch%aux%RT%rt_parameter%nphase
  
  if (realization%option%transport%tvd_flux_limiter /= TVD_LIMITER_UPWIND) then
    allocate(total_up2(nphase,ntvddof))
    allocate(total_dn2(nphase,ntvddof))
  else
    ! these must be nullifed so that the explicit scheme ignores them
    nullify(total_up2)
    nullify(total_up2)
  endif

  ! load total component concentrations into tran_xx_p.  it will be used
  ! as local storage here and eventually be overwritten upon leaving 
  ! this routine
  call VecGetArrayF90(field%tran_xx,tran_xx_p,ierr);CHKERRQ(ierr)
  tran_xx_p = 0.d0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    local_end = local_id*ntvddof
    local_start = local_end-ntvddof+1
    do iphase = 1, nphase
      tran_xx_p(local_start:local_end) = &
        rt_auxvars(ghosted_id)%total(:,iphase)
    enddo
  enddo
  call VecRestoreArrayF90(field%tran_xx,tran_xx_p,ierr);CHKERRQ(ierr)
  call VecScatterBegin(discretization%tvd_ghost_scatter,field%tran_xx, &
                       field%tvd_ghosts,INSERT_VALUES,SCATTER_FORWARD, &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd(discretization%tvd_ghost_scatter,field%tran_xx, &
                     field%tvd_ghosts,INSERT_VALUES,SCATTER_FORWARD, &
                     ierr);CHKERRQ(ierr)

! Update Boundary Concentrations------------------------------
  call VecGetArrayF90(field%tvd_ghosts,tvd_ghosts_p,ierr);CHKERRQ(ierr)
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    if (associated(cur_connection_set%id_dn2)) then
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        id_dn2 = cur_connection_set%id_dn2(iconn)
        if (id_dn2 < 0) then
          iend = abs(id_dn2)*ntvddof
          istart = iend-ntvddof+1
          tvd_ghosts_p(istart:iend) = rt_auxvars_bc(sum_connection)%total(1,:)
        endif
      enddo
    endif
    boundary_condition => boundary_condition%next
  enddo  
  call VecRestoreArrayF90(field%tvd_ghosts,tvd_ghosts_p,ierr);CHKERRQ(ierr)
#if TVD_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'tvd_ghosts.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call VecView(field%tvd_ghosts,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  sum_flux = 0.d0
  
  if (reaction%ncoll > 0) then
    option%io_buffer = &
      'Need to add colloidal source/sinks to RTExplicitAdvection()'
    call printErrMsg(option)
  endif
  if (patch%aux%RT%rt_parameter%nphase > 1) then
    option%io_buffer = &
      'Need to add multiphase source/sinks to RTExplicitAdvection()'
    call printErrMsg(option)
  endif
  if (reaction%ncomp /= reaction%naqcomp) then
    option%io_buffer = &
      'Need to account for non-aqueous species to RTExplicitAdvection()'
    call printErrMsg(option)
  endif
  if (option%compute_mass_balance_new) then  
    option%io_buffer = &
      'Mass balance not yet supported in RTExplicitAdvection()'
    call printErrMsg(option)
  endif
  
! Interior Flux Terms -----------------------------------
  unitvec_xyz(1:2) = 0.d0
  unitvec_xyz(3)   = 1.d0

  call VecGetArrayF90(field%tvd_ghosts,tvd_ghosts_p,ierr);CHKERRQ(ierr)
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle

      !----------------
      if ( option%flow%only_vertical_flow .or. option%transport%only_vertical_tran ) then
         if(abs(dot_product(cur_connection_set%dist(1:3,iconn),unitvec_xyz)) < 1.d-20) cycle
      end if
      !----------------
        
      if (associated(cur_connection_set%id_dn2)) then
        id_up2 = cur_connection_set%id_up2(iconn)
        if (id_up2 > 0) then
          total_up2 = rt_auxvars(id_up2)%total
        else
          iend = abs(id_up2)*ntvddof
          istart = iend-ntvddof+1
          total_up2(1,:) = tvd_ghosts_p(istart:iend)
        endif
        id_dn2 = cur_connection_set%id_dn2(iconn)
        if (id_dn2 > 0) then
          total_dn2 = rt_auxvars(id_dn2)%total
        else
          iend = abs(id_dn2)*ntvddof
          istart = iend-ntvddof+1
          total_dn2(1,:) = tvd_ghosts_p(istart:iend)
        endif
      endif
      call TFluxTVD(patch%aux%RT%rt_parameter, &
                    patch%internal_velocities(:,sum_connection), &
                    cur_connection_set%area(iconn), &
                    cur_connection_set%dist(:,iconn), &
                    total_up2, &
                    rt_auxvars(ghosted_id_up), &
                    rt_auxvars(ghosted_id_dn), &
                    total_dn2, &
                    TFluxLimitPtr, &
                    option,flux)
          
      ! contribution upwind
      sum_flux(:,ghosted_id_up) = sum_flux(:,ghosted_id_up) - flux
        
      ! contribution downwind
      sum_flux(:,ghosted_id_dn) = sum_flux(:,ghosted_id_dn) + flux
          
    enddo ! iconn
    cur_connection_set => cur_connection_set%next
  enddo
  call VecRestoreArrayF90(field%tvd_ghosts,tvd_ghosts_p,ierr);CHKERRQ(ierr)
    
! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      if (associated(cur_connection_set%id_dn2)) then
        total_up2 = rt_auxvars_bc(sum_connection)%total
        id_dn2 = cur_connection_set%id_dn2(iconn)
        if (id_dn2 > 0) then
          total_dn2 = rt_auxvars(id_dn2)%total
        else
          iend = abs(id_dn2)*ntvddof
          istart = iend-ntvddof+1
          total_dn2(1,:) = tvd_ghosts_p(istart:iend)
        endif
      endif
      call TFluxTVD(patch%aux%RT%rt_parameter, &
                    patch%boundary_velocities(:,sum_connection), &
                    cur_connection_set%area(iconn), &
                    cur_connection_set%dist(:,iconn), &
                    total_up2, &
                    rt_auxvars_bc(sum_connection), &
                    rt_auxvars(ghosted_id), &
                    total_dn2, &
                    TFluxLimitPtr, &
                    option,flux)

      ! contribution downwind
      sum_flux(:,ghosted_id) = sum_flux(:,ghosted_id) + flux
     
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  qsrc = 0.d0
  do 
    if (.not.associated(source_sink)) exit
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections 
      sum_connection = sum_connection + 1     
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      if (associated(patch%ss_flow_vol_fluxes)) then
        qsrc(1:nphase) = patch%ss_flow_vol_fluxes(1:nphase,sum_connection)
      endif
      call TSrcSinkCoef(rt_parameter,qsrc,source_sink%tran_condition%itype, &
                        coef_in,coef_out)
      do iphase = 1, nphase
        flux = coef_in*rt_auxvars(ghosted_id)%total(:,iphase) + &
               coef_out*source_sink%tran_condition%cur_constraint_coupler% &
                                          rt_auxvar%total(:,iphase)
        !geh: TSrcSinkCoef() unit are in L/s.
         sum_flux(:,ghosted_id) = sum_flux(:,ghosted_id) + flux
      enddo
    enddo
    source_sink => source_sink%next
  enddo
  
  call VecGetArrayF90(field%tran_xx,tran_xx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%tran_rhs_coef,rhs_coef_p,ierr);CHKERRQ(ierr)

  
  ! update concentration
  iphase = 1
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    local_end = local_id*ntvddof
    local_start = local_end-ntvddof+1
      ! psv_t must have same units [mol/sec] and be consistent with rhs_coef_p
      ! in RTUpdateRHSCoefs()
      psv_t = material_auxvars(ghosted_id)%porosity* &
              global_auxvars(ghosted_id)%sat(iphase)* &
              1000.d0* &
              material_auxvars(ghosted_id)%volume/option%tran_dt
      !geh: clearly dangerous that I reload into total, but I am going to do it!
      tran_xx_p(local_start:local_end) = &
        ((rhs_coef_p(local_id)*rt_auxvars(ghosted_id)%total(:,iphase)) + &
         sum_flux(:,ghosted_id)) / psv_t
!    enddo
  enddo
  
  if (associated(total_up2)) then
    deallocate(total_up2)
    nullify(total_up2)
    deallocate(total_dn2)
    nullify(total_dn2)
  endif
  
  ! Restore vectors
  call VecRestoreArrayF90(field%tran_xx,tran_xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%tran_rhs_coef,rhs_coef_p, &
                              ierr);CHKERRQ(ierr)
  
end subroutine RTExplicitAdvection

! ************************************************************************** !

subroutine RTWriteToHeader(fid,variable_string,cell_string,icolumn)
  ! 
  ! Appends formatted strings to header string
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/27/11
  ! 

  PetscInt :: fid
  character(len=*) :: variable_string
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXSTRINGLENGTH) :: variable_string_adj
  character(len=MAXWORDLENGTH) :: column_string
  PetscInt :: icolumn

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: len_cell_string

  variable_string_adj = variable_string
  !geh: Shift to left.  Cannot perform on same string since len=*
  variable_string_adj = adjustl(variable_string_adj)
  
  if (icolumn > 0) then
    icolumn = icolumn + 1
    write(column_string,'(i4,''-'')') icolumn
    column_string = trim(adjustl(column_string))
  else
    column_string = ''
  endif

  !geh: this is all to remove the lousy spaces
  len_cell_string = len_trim(cell_string) 

  if (len_cell_string > 0) then
    write(string,'('',"'',a,a,'' '',a,''"'')') trim(column_string), &
          trim(variable_string_adj), trim(cell_string)
  else
    write(string,'('',"'',a,a,''"'')') trim(column_string), &
          trim(variable_string_adj)
  endif
  write(fid,'(a)',advance="no") trim(string)

end subroutine RTWriteToHeader

! ************************************************************************** !

subroutine RTClearActivityCoefficients(realization)
  ! 
  ! Sets activity coefficients back to 1.
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/11/14
  ! 

  use Realization_Subsurface_class
  use Reactive_Transport_Aux_module
  use Option_module
  use Field_module  
  use Grid_module
  use Secondary_Continuum_Aux_module  

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscInt :: ghosted_id
  
  rt_auxvars => realization%patch%aux%RT%auxvars
  
  do ghosted_id = 1, realization%patch%grid%ngmax
    rt_auxvars(ghosted_id)%pri_act_coef = 1.d0
    if (associated(rt_auxvars(ghosted_id)%sec_act_coef)) then
      rt_auxvars(ghosted_id)%sec_act_coef = 1.d0
    endif
  enddo

end subroutine RTClearActivityCoefficients

! ************************************************************************** !

subroutine RTDestroy(realization)
  ! 
  ! Deallocates variables associated with Reactive Transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/09
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module

  type(realization_subsurface_type) :: realization
  
#ifdef OS_STATISTICS
  type(option_type), pointer :: option
  PetscErrorCode :: ierr
  
  PetscReal :: temp_real_in(3), temp_real_out(3)
  PetscReal :: call_count
  PetscReal :: sum_newton_iterations
  PetscReal :: ave_newton_iterations_in_a_cell
  PetscInt :: max_newton_iterations_in_a_cell
  PetscReal :: max_newton_iterations_on_a_core
  PetscReal :: min_newton_iterations_on_a_core
  
  PetscReal :: sum, ave, var, value
  PetscInt :: irank
  PetscReal, allocatable :: tot_newton_iterations(:)
  
  option => realization%option
  call_count = 0.d0
  sum_newton_iterations = 0.d0
  max_newton_iterations_in_a_cell = -99999999
  max_newton_iterations_on_a_core = -99999999.d0
  min_newton_iterations_on_a_core = 99999999.d0
  
#endif  

#ifdef OS_STATISTICS
  call_count = call_count + &
    cur_patch%aux%RT%rt_parameter%sum_newton_call_count
  sum_newton_iterations = sum_newton_iterations + &
    cur_patch%aux%RT%rt_parameter%sum_newton_iterations
  if (cur_patch%aux%RT%rt_parameter%overall_max_newton_iterations > &
      max_newton_iterations_in_a_cell) then
    max_newton_iterations_in_a_cell = &
      cur_patch%aux%RT%rt_parameter%overall_max_newton_iterations
  endif
#endif

#ifdef OS_STATISTICS
  if (option%reactive_transport_coupling == OPERATOR_SPLIT) then
    temp_real_in(1) = call_count
    temp_real_in(2) = sum_newton_iterations
    call MPI_Allreduce(temp_real_in,temp_real_out,TWO_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
    ave_newton_iterations_in_a_cell = temp_real_out(2)/temp_real_out(1)

    temp_real_in(1) = dble(max_newton_iterations_in_a_cell)
    temp_real_in(2) = sum_newton_iterations
    temp_real_in(3) = -sum_newton_iterations
    call MPI_Allreduce(temp_real_in,temp_real_out,THREE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
    max_newton_iterations_in_a_cell = int(temp_real_out(1)+1.d-4)
    max_newton_iterations_on_a_core = temp_real_out(2)
    min_newton_iterations_on_a_core = -temp_real_out(3)

    ! Now let's compute the variance!
    call OptionMeanVariance(sum_newton_iterations,ave,var,PETSC_TRUE,option)
  
    if (option%print_screen_flag) then
      write(*, '(/,/" OS Reaction Statistics (Overall): ",/, &
               & "       Ave Newton Its / Cell: ",1pe12.4,/, &
               & "       Max Newton Its / Cell: ",i4,/, &
               & "       Max Newton Its / Core: ",1pe12.4,/, &
               & "       Min Newton Its / Core: ",1pe12.4,/, &
               & "       Ave Newton Its / Core: ",1pe12.4,/, &
               & "   Std Dev Newton Its / Core: ",1pe12.4,/)') &
                 ave_newton_iterations_in_a_cell, &
                 max_newton_iterations_in_a_cell, &
                 max_newton_iterations_on_a_core, &
                 min_newton_iterations_on_a_core, &
                 ave, &
                 sqrt(var)

    endif

    if (option%print_file_flag) then
      write(option%fid_out, '(/,/" OS Reaction Statistics (Overall): ",/, &
               & "       Ave Newton Its / Cell: ",1pe12.4,/, &
               & "       Max Newton Its / Cell: ",i4,/, &
               & "       Max Newton Its / Core: ",1pe12.4,/, &
               & "       Min Newton Its / Core: ",1pe12.4,/, &
               & "       Ave Newton Its / Core: ",1pe12.4,/, &
               & "   Std Dev Newton Its / Core: ",1pe12.4,/)') &
                 ave_newton_iterations_in_a_cell, &
                 max_newton_iterations_in_a_cell, &
                 max_newton_iterations_on_a_core, &
                 min_newton_iterations_on_a_core, &
                 ave, &
                 sqrt(var)
    endif
  endif

#endif 


end subroutine RTDestroy

end module Reactive_Transport_module
