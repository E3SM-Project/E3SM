module TH_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use TH_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none
  
  private 

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: unit_z(3) = [0.d0,0.d0,1.d0]

  public THResidual,THJacobian, &
         THUpdateFixedAccumulation,THTimeCut,&
         THSetup, THNumericalJacobianTest, &
         THMaxChange, THUpdateSolution, &
         THGetTecplotHeader, THInitializeTimestep, &
         THComputeMassBalance, THResidualToMass, &
         THSecondaryHeat, THSecondaryHeatJacobian, & 
         THUpdateAuxVars, THDestroy, &
         THUpdateSurfaceBC, THAccumulation
         
  PetscInt, parameter :: jh2o = 1

contains

! ************************************************************************** !

subroutine THTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Field_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
 
  call THInitializeTimestep(realization)
 
end subroutine THTimeCut

! ************************************************************************** !

subroutine THSetup(realization)
  ! 
  ! Author: ???
  ! Date: 02/22/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Output_Aux_module

  type(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  type(output_variable_list_type), pointer :: list
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THSetupPatch(realization)
    cur_patch => cur_patch%next
  enddo

  list => realization%output_option%output_snap_variable_list
  call THSetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call THSetPlotVariables(realization,list)

end subroutine THSetup

! ************************************************************************** !

subroutine THSetupPatch(realization)
  ! 
  ! Creates arrays for auxiliary variables
  ! 
  ! Author: ???
  ! Date: 02/22/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Region_module
  use Coupler_module
  use Connection_module
  use Fluid_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
  use Characteristic_Curves_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(TH_auxvar_type), pointer :: TH_auxvars(:), TH_auxvars_bc(:)
  type(TH_auxvar_type), pointer :: TH_auxvars_ss(:)
  type(fluid_property_type), pointer :: cur_fluid_property
  type(sec_heat_type), pointer :: TH_sec_heat_vars(:)
  type(coupler_type), pointer :: initial_condition
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: area_per_vol

  PetscInt :: ghosted_id, iconn, sum_connection
  PetscInt :: i, iphase, local_id, material_id
  PetscBool :: error_found
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
    
  patch%aux%TH => THAuxCreate(option)
  patch%aux%SC_heat => SecondaryAuxHeatCreate(option)

! option%io_buffer = 'Before TH can be run, the TH_parameter object ' // &
!                    'must be initialized with the proper variables ' // &
!                    'THAuxCreate() is called anywhere.'
! call printErrMsg(option)

  allocate(patch%aux%TH%TH_parameter%sir(option%nphase, &
                                  size(patch%characteristic_curves_array)))

  !Jitu, 08/04/2010: Check these allocations. Currently assumes only single value in the array	<modified pcl 1-13-11>
  allocate(patch%aux%TH%TH_parameter%dencpr(size(patch%material_property_array)))
  allocate(patch%aux%TH%TH_parameter%ckwet(size(patch%material_property_array)))
  allocate(patch%aux%TH%TH_parameter%ckdry(size(patch%material_property_array)))
  allocate(patch%aux%TH%TH_parameter%alpha(size(patch%material_property_array)))
  !if (option%use_th_freezing) then
     allocate(patch%aux%TH%TH_parameter%ckfrozen(size(patch%material_property_array)))
     allocate(patch%aux%TH%TH_parameter%alpha_fr(size(patch%material_property_array)))
  !endif

  !Copy the values in the TH_parameter from the global realization 
  error_found = PETSC_FALSE
  do i = 1, size(patch%material_property_array)
    word = patch%material_property_array(i)%ptr%name 
    if (Uninitialized(patch%material_property_array(i)%ptr%specific_heat)) then
      option%io_buffer = 'ERROR: Non-initialized HEAT_CAPACITY in material ' &
                         // trim(word)
      call printMsg(option)
      error_found = PETSC_TRUE
    endif
    if (Uninitialized(patch%material_property_array(i)%ptr% &
                      thermal_conductivity_wet)) then
      option%io_buffer = 'ERROR: Non-initialized THERMAL_CONDUCTIVITY_WET in &
                         &material ' // &
                         trim(word)
      call printMsg(option)
      error_found = PETSC_TRUE
    endif
    if (Uninitialized(patch%material_property_array(i)%ptr% &
                      thermal_conductivity_dry)) then
      option%io_buffer = 'ERROR: Non-initialized THERMAL_CONDUCTIVITY_DRY in &
                         &material ' // &
                         trim(word)
      call printMsg(option)
      error_found = PETSC_TRUE
    endif
    if (option%use_th_freezing) then
      if (Uninitialized(patch%material_property_array(i)%ptr% &
                        thermal_conductivity_frozen)) then
        option%io_buffer = 'ERROR: Non-initialized THERMAL_CONDUCTIVITY_&
                           &FROZEN in material ' // trim(word)
        call printMsg(option)
        error_found = PETSC_TRUE
      endif
      if (Uninitialized(patch%material_property_array(i)%ptr% &
                        alpha_fr)) then
        option%io_buffer = 'ERROR: Non-initialized THERMAL_COND_EXPONENT&
                           &_FROZEN in material ' // trim(word)
        call printMsg(option)
        error_found = PETSC_TRUE
      endif
    endif
    material_id = abs(patch%material_property_array(i)%ptr%internal_id)
    ! kg rock/m^3 rock * J/kg rock-K * 1.e-6 MJ/J = MJ/m^3-K
    patch%aux%TH%TH_parameter%dencpr(material_id) = &
      patch%material_property_array(i)%ptr%rock_density*option%scale* &
        patch%material_property_array(i)%ptr%specific_heat
    patch%aux%TH%TH_parameter%ckwet(material_id) = &
      patch%material_property_array(i)%ptr%thermal_conductivity_wet* &
      option%scale  
    patch%aux%TH%TH_parameter%ckdry(material_id) = &
      patch%material_property_array(i)%ptr%thermal_conductivity_dry* &
      option%scale
    patch%aux%TH%TH_parameter%alpha(material_id) = &
      patch%material_property_array(i)%ptr%alpha
    if (option%use_th_freezing) then
      patch%aux%TH%TH_parameter%ckfrozen(material_id) = &
        patch%material_property_array(i)%ptr%thermal_conductivity_frozen* &
        option%scale
      patch%aux%TH%TH_parameter%alpha_fr(material_id) = &
        patch%material_property_array(i)%ptr%alpha_fr
    endif

  enddo 

  if (error_found) then
    option%io_buffer = 'Material property errors found in THSetup.'
    call printErrMsg(option)
  endif

  do i = 1, size(patch%characteristic_curves_array)
    patch%aux%TH%TH_parameter%sir(:,i) = &
        CharCurvesGetGetResidualSats(patch%characteristic_curves_array(i)%ptr,option)
  enddo

  ! allocate auxvar data structures for all grid cells
  allocate(TH_auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call THAuxVarInit(TH_auxvars(ghosted_id),option)
  enddo

 
  if (option%use_mc) then
    initial_condition => patch%initial_condition_list%first
    allocate(TH_sec_heat_vars(grid%nlmax))
  
    do local_id = 1, grid%nlmax
  
    ! Assuming the same secondary continuum for all regions (need to make it an array)
    ! S. Karra 07/18/12
      call SecondaryContinuumSetProperties( &
        TH_sec_heat_vars(local_id)%sec_continuum, &
        patch%material_property_array(1)%ptr%secondary_continuum_name, &
        patch%material_property_array(1)%ptr%secondary_continuum_length, &
        patch%material_property_array(1)%ptr%secondary_continuum_matrix_block_size, &
        patch%material_property_array(1)%ptr%secondary_continuum_fracture_spacing, &
        patch%material_property_array(1)%ptr%secondary_continuum_radius, &
        patch%material_property_array(1)%ptr%secondary_continuum_area, &
        option)
        
      TH_sec_heat_vars(local_id)%ncells = &
        patch%material_property_array(1)%ptr%secondary_continuum_ncells
      TH_sec_heat_vars(local_id)%aperture = &
        patch%material_property_array(1)%ptr%secondary_continuum_aperture
      TH_sec_heat_vars(local_id)%epsilon = &
        patch%material_property_array(1)%ptr%secondary_continuum_epsilon
      TH_sec_heat_vars(local_id)%log_spacing = &
        patch%material_property_array(1)%ptr%secondary_continuum_log_spacing
      TH_sec_heat_vars(local_id)%outer_spacing = &
        patch%material_property_array(1)%ptr%secondary_continuum_outer_spacing
                
      allocate(TH_sec_heat_vars(local_id)%area(TH_sec_heat_vars(local_id)%ncells))
      allocate(TH_sec_heat_vars(local_id)%vol(TH_sec_heat_vars(local_id)%ncells))
      allocate(TH_sec_heat_vars(local_id)%dm_minus(TH_sec_heat_vars(local_id)%ncells))
      allocate(TH_sec_heat_vars(local_id)%dm_plus(TH_sec_heat_vars(local_id)%ncells))
      allocate(TH_sec_heat_vars(local_id)%sec_continuum% &
             distance(TH_sec_heat_vars(local_id)%ncells))
    
      call SecondaryContinuumType(TH_sec_heat_vars(local_id)%sec_continuum, &
                                  TH_sec_heat_vars(local_id)%ncells, &
                                  TH_sec_heat_vars(local_id)%area, &
                                  TH_sec_heat_vars(local_id)%vol, &
                                  TH_sec_heat_vars(local_id)%dm_minus, &
                                  TH_sec_heat_vars(local_id)%dm_plus, &
                                  TH_sec_heat_vars(local_id)%aperture, &
                                  TH_sec_heat_vars(local_id)%epsilon, &
                                  TH_sec_heat_vars(local_id)%log_spacing, &
                                  TH_sec_heat_vars(local_id)%outer_spacing, &
                                  area_per_vol,option)
                                
      TH_sec_heat_vars(local_id)%interfacial_area = area_per_vol* &
        (1.d0 - TH_sec_heat_vars(local_id)%epsilon)* &
        patch%material_property_array(1)%ptr% &
        secondary_continuum_area_scaling

    ! Setting the initial values of all secondary node temperatures same as primary node 
    ! temperatures (with initial dirichlet BC only) -- sk 06/26/12
      allocate(TH_sec_heat_vars(local_id)%sec_temp(TH_sec_heat_vars(local_id)%ncells))
      
      if (option%set_secondary_init_temp) then
        TH_sec_heat_vars(local_id)%sec_temp = &
          patch%material_property_array(1)%ptr%secondary_continuum_init_temp
      else
        TH_sec_heat_vars(local_id)%sec_temp = &
        initial_condition%flow_condition%temperature%dataset%rarray(1)
      endif
          
    enddo
      
    patch%aux%SC_heat%sec_heat_vars => TH_sec_heat_vars    

  endif

  
  patch%aux%TH%auxvars => TH_auxvars
  patch%aux%TH%num_aux = grid%ngmax
  
  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  boundary_condition => patch%boundary_condition_list%first

  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo

  if (sum_connection > 0) then 
    allocate(TH_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call THAuxVarInit(TH_auxvars_bc(iconn),option)
    enddo
    patch%aux%TH%auxvars_bc => TH_auxvars_bc
  endif
  patch%aux%TH%num_aux_bc = sum_connection

  ! Create aux vars for source/sink
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(TH_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call THAuxVarInit(TH_auxvars_ss(iconn),option)
    enddo
    patch%aux%TH%auxvars_ss => TH_auxvars_ss
  endif
  patch%aux%TH%num_aux_ss = sum_connection
  
  ! initialize parameters
  cur_fluid_property => realization%fluid_properties
  do 
    if (.not.associated(cur_fluid_property)) exit
    iphase = cur_fluid_property%phase_id
    cur_fluid_property => cur_fluid_property%next
  enddo

end subroutine THSetupPatch

! ************************************************************************** !

subroutine THComputeMassBalance(realization, mass_balance)
  ! 
  ! THomputeMassBalance:
  ! Adapted from RichardsComputeMassBalance: need to be checked
  ! 
  ! Author: Jitendra Kumar
  ! Date: 07/21/2010
  ! 

  use Realization_Subsurface_class
  use Patch_module

  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)
   
  type(patch_type), pointer :: cur_patch

  mass_balance = 0.d0

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THComputeMassBalancePatch(realization, mass_balance)
    cur_patch => cur_patch%next
  enddo

end subroutine THComputeMassBalance    

! ************************************************************************** !

subroutine THComputeMassBalancePatch(realization,mass_balance)
  ! 
  ! THomputeMassBalancePatch:
  ! Adapted from RichardsComputeMassBalancePatch: need to be checked
  ! 
  ! Author: Jitendra Kumar
  ! Date: 07/21/2010
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_class, only : material_auxvar_type, &
                                 soil_compressibility_index, &
                                 MaterialCompressSoil
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(TH_auxvar_type),pointer :: TH_auxvars(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscReal :: compressed_porosity
  PetscReal :: por
  PetscReal :: dum1

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  TH_auxvars => patch%aux%TH%auxvars

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    ! mass = volume*saturation*density

    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvars(ghosted_id), &
                                global_auxvars(ghosted_id)%pres(1), &
                                compressed_porosity,dum1)
      por = compressed_porosity
    else
      por = material_auxvars(ghosted_id)%porosity
    endif

    mass_balance = mass_balance + &
      global_auxvars(ghosted_id)%den_kg* &
      global_auxvars(ghosted_id)%sat* &
      por* &
      material_auxvars(ghosted_id)%volume

    if (option%use_th_freezing) then
      ! mass = volume*saturation_ice*density_ice
      mass_balance = mass_balance + &
        TH_auxvars(ghosted_id)%ice%den_ice*FMWH2O* &
        TH_auxvars(ghosted_id)%ice%sat_ice* &
        por* &
        material_auxvars(ghosted_id)%volume
    endif

  enddo

end subroutine THComputeMassBalancePatch

! ************************************************************************** !

subroutine THZeroMassBalDeltaPatch(realization)
  ! 
  ! Zeros mass balance delta array
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 12/13/11
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%TH%num_aux
    patch%aux%Global%auxvars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%TH%num_aux_bc > 0) then
    do iconn = 1, patch%aux%TH%num_aux_bc
      global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
  if (patch%aux%TH%num_aux_ss > 0) then
    do iconn = 1, patch%aux%TH%num_aux_ss
      global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
 
end subroutine THZeroMassBalDeltaPatch

! ************************************************************************** !

subroutine THUpdateMassBalancePatch(realization)
  ! 
  ! Updates mass balance
  ! 
  ! Author: ???
  ! Date: 12/13/11
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%TH%num_aux
    patch%aux%Global%auxvars(iconn)%mass_balance = &
      patch%aux%Global%auxvars(iconn)%mass_balance + &
      patch%aux%Global%auxvars(iconn)%mass_balance_delta*FMWH2O* &
      option%flow_dt
  enddo
#endif

  if (patch%aux%TH%num_aux_bc > 0) then
    do iconn = 1, patch%aux%TH%num_aux_bc
      global_auxvars_bc(iconn)%mass_balance = &
        global_auxvars_bc(iconn)%mass_balance + &
        global_auxvars_bc(iconn)%mass_balance_delta*FMWH2O*option%flow_dt
    enddo
  endif
  if (patch%aux%TH%num_aux_ss > 0) then
    do iconn = 1, patch%aux%TH%num_aux_ss
      global_auxvars_ss(iconn)%mass_balance = &
        global_auxvars_ss(iconn)%mass_balance + &
        global_auxvars_ss(iconn)%mass_balance_delta*FMWH2O*option%flow_dt
    enddo
  endif


end subroutine THUpdateMassBalancePatch

! ************************************************************************** !

subroutine THUpdateAuxVars(realization)
  ! 
  ! Updates the auxiliary variables associated with
  ! the TH problem
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module

  type(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THUpdateAuxVarsPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine THUpdateAuxVars

! ************************************************************************** !

subroutine THUpdateAuxVarsPatch(realization)
  ! 
  ! Updates the auxiliary variables associated with
  ! the TH problem
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
   
  implicit none

  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(TH_auxvar_type), pointer :: TH_auxvars(:)
  type(TH_auxvar_type), pointer :: TH_auxvars_bc(:)
  type(TH_auxvar_type), pointer :: TH_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  !class(material_auxvar_type), pointer :: material_auxvars(:)  ! this unknownly would likely mess up %soil_properties(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(TH_parameter_type), pointer :: TH_parameter

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscReal, pointer :: xx(:)
  PetscReal :: tsrc1
  PetscErrorCode :: ierr
  PetscInt :: ithrm
  PetscReal, pointer :: ithrm_loc_p(:)

  !!
!  PetscReal, allocatable :: gradient(:,:)
  !!
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  
  !!
!  allocate(gradient(grid%ngmax,3))
!  gradient = 0.d0
  !!
  
  TH_auxvars => patch%aux%TH%auxvars
  TH_auxvars_bc => patch%aux%TH%auxvars_bc
  TH_auxvars_ss => patch%aux%TH%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  TH_parameter => patch%aux%TH%TH_parameter

  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    if (patch%imat(ghosted_id) <= 0) cycle
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1
    iphase = int(iphase_loc_p(ghosted_id))
    ithrm = int(ithrm_loc_p(ghosted_id))

    if (option%use_th_freezing) then
       call THAuxVarComputeFreezing(xx_loc_p(istart:iend), &
            TH_auxvars(ghosted_id),global_auxvars(ghosted_id), &
            material_auxvars(ghosted_id), &
            iphase, &
            patch%characteristic_curves_array(int(icap_loc_p(ghosted_id)))%ptr, &
            TH_parameter, ithrm, &
            option)
    else
       call THAuxVarComputeNoFreezing(xx_loc_p(istart:iend), &
            TH_auxvars(ghosted_id),global_auxvars(ghosted_id), &
            material_auxvars(ghosted_id), &
            iphase, &
            patch%characteristic_curves_array(int(icap_loc_p(ghosted_id)))%ptr, &
            TH_parameter, ithrm, &
            option)
    endif

    iphase_loc_p(ghosted_id) = iphase
  enddo

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
      ithrm = int(ithrm_loc_p(ghosted_id))

      do idof=1,option%nflowdof
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC, &
               HET_SURF_SEEPAGE_BC, &
               HET_DIRICHLET_BC,HET_SEEPAGE_BC,HET_CONDUCTANCE_BC)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)

            ! a special situation of flow-in, in which NO heat-conduction/diffusion,
            ! but with bulk-heat flux defined by fluid itself rather than boundary-cell.
            if (boundary_condition%flow_condition%itype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC .and. &
                boundary_condition%flow_aux_real_var(TH_PRESSURE_DOF,iconn) > floweps) then

              if(boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) /= UNINITIALIZED_DOUBLE) &
              xxbc(TH_TEMPERATURE_DOF) = boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn)
            endif

        end select
      enddo
      
      select case(boundary_condition%flow_condition%itype(TH_PRESSURE_DOF))
        case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC, &
             HET_DIRICHLET_BC,HET_SEEPAGE_BC,HET_CONDUCTANCE_BC)
          iphasebc = boundary_condition%flow_aux_int_var(1,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          iphasebc=int(iphase_loc_p(ghosted_id))                               
      end select

      if (option%use_th_freezing) then
         call THAuxVarComputeFreezing(xxbc,TH_auxvars_bc(sum_connection), &
              global_auxvars_bc(sum_connection), &
              material_auxvars(ghosted_id), &
              iphasebc, &
              patch%characteristic_curves_array(int(icap_loc_p(ghosted_id)))%ptr, &
              TH_parameter, ithrm, &
              option)
      else
         call THAuxVarComputeNoFreezing(xxbc,TH_auxvars_bc(sum_connection), &
              global_auxvars_bc(sum_connection), &
              material_auxvars(ghosted_id), &
              iphasebc, &
              patch%characteristic_curves_array(int(icap_loc_p(ghosted_id)))%ptr, &
              TH_parameter, ithrm, &
              option)
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! source/sinks
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  allocate(xx(option%nflowdof))
  do
    if (.not.associated(source_sink)) exit
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      ithrm = int(ithrm_loc_p(ghosted_id))

      iend = ghosted_id*option%nflowdof
      istart = iend-option%nflowdof+1
      iphase = int(iphase_loc_p(ghosted_id))

      select case(source_sink%flow_condition%itype(TH_TEMPERATURE_DOF))
        case (HET_DIRICHLET_BC)
          tsrc1 = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
        case (DIRICHLET_BC)
          tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
        case (ENERGY_RATE_SS,SCALED_ENERGY_RATE_SS,HET_ENERGY_RATE_SS, &
              ZERO_GRADIENT_BC)
          tsrc1 = xx_loc_p((ghosted_id-1)*option%nflowdof+2)
        case default
          option%io_buffer='Unsupported temperature flow condtion for ' // &
            'a source-sink in TH mode: ' // trim(source_sink%name)
          call printErrMsg(option)
      end select

      xx(1) = xx_loc_p(istart)
      xx(2) = tsrc1

      if (option%use_th_freezing) then
         call THAuxVarComputeFreezing(xx, &
              TH_auxvars_ss(sum_connection),global_auxvars_ss(sum_connection), &
              material_auxvars(ghosted_id), &
              iphase, &
              patch%characteristic_curves_array(int(icap_loc_p(ghosted_id)))%ptr, &
              TH_parameter, ithrm, &
              option)
      else
         call THAuxVarComputeNoFreezing(xx, &
              TH_auxvars_ss(sum_connection),global_auxvars_ss(sum_connection), &
              material_auxvars(ghosted_id), &
              iphase, &
              patch%characteristic_curves_array(int(icap_loc_p(ghosted_id)))%ptr, &
              TH_parameter, ithrm, &
              option)
      endif
    enddo
    source_sink => source_sink%next
  enddo
  deallocate(xx)

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)

  patch%aux%TH%auxvars_up_to_date = PETSC_TRUE

  ! Update a flag marking presence or absence of standing water for BC grid cells
  if (option%surf_flow_on) call THUpdateSurfaceWaterFlag(realization)

end subroutine THUpdateAuxVarsPatch

! ************************************************************************** !

subroutine THInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: ???
  ! Date: 02/20/08
  ! 

  use Realization_Subsurface_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  call THUpdateFixedAccumulation(realization)

end subroutine THInitializeTimestep

! ************************************************************************** !

subroutine THUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time step
  ! 
  ! Author: ???
  ! Date: 02/13/08
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(field_type), pointer :: field
  type(patch_type), pointer :: cur_patch
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  
  field => realization%field
    
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THUpdateSolutionPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine THUpdateSolution

! ************************************************************************** !

subroutine THUpdateSolutionPatch(realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 12/13/11, 02/28/14
  ! 


  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
    
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(TH_parameter_type), pointer :: TH_parameter
  type(TH_auxvar_type), pointer :: auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(sec_heat_type), pointer :: TH_sec_heat_vars(:)

  PetscInt :: istart, iend  
  PetscInt :: local_id, ghosted_id
  ! secondary continuum variables
  PetscReal :: sec_dencpr
  PetscErrorCode :: ierr
  PetscReal, pointer :: ithrm_loc_p(:)

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  TH_parameter => patch%aux%TH%TH_parameter
  auxvars => patch%aux%TH%auxvars
  global_auxvars => patch%aux%Global%auxvars
  
  if (option%use_mc) then
    TH_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars
  endif

  if (realization%option%compute_mass_balance_new) then
    call THUpdateMassBalancePatch(realization)
  endif

  if (option%use_mc) then
    call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
    
      sec_dencpr = TH_parameter%dencpr(int(ithrm_loc_p(local_id))) ! secondary rho*c_p same as primary for now
        
      call THSecHeatAuxVarCompute(TH_sec_heat_vars(local_id), &
                            global_auxvars(ghosted_id), &
                            TH_parameter%ckwet(int(ithrm_loc_p(local_id))), &
                            sec_dencpr, &
                            option)
                            
    enddo 
    call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
  endif


end subroutine THUpdateSolutionPatch

! ************************************************************************** !

subroutine THUpdateFixedAccumulation(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module

  type(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THUpdateFixedAccumPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine THUpdateFixedAccumulation

! ************************************************************************** !

subroutine THUpdateFixedAccumPatch(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Secondary_Continuum_Aux_module


  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(TH_auxvar_type), pointer :: TH_auxvars(:)
  type(TH_parameter_type), pointer :: TH_parameter
  type(sec_heat_type), pointer :: TH_sec_heat_vars(:)
  !class(material_auxvar_type), pointer :: material_auxvars(:)  ! this unknownly would likely mess up %soil_properties(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, istart, iend, iphase
  PetscReal, pointer :: xx_p(:), icap_loc_p(:), iphase_loc_p(:)
  PetscReal, pointer :: ithrm_loc_p(:), accum_p(:)
  PetscReal :: vol_frac_prim
  PetscInt :: ithrm
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  TH_parameter => patch%aux%TH%TH_parameter
  TH_auxvars => patch%aux%TH%auxvars
  global_auxvars => patch%aux%Global%auxvars
  TH_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)


  vol_frac_prim = 1.d0
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    iphase = int(iphase_loc_p(ghosted_id))
    ithrm = int(ithrm_loc_p(ghosted_id))

    if (option%use_th_freezing) then
       call THAuxVarComputeFreezing(xx_p(istart:iend), &
            TH_auxvars(ghosted_id),global_auxvars(ghosted_id), &
            material_auxvars(ghosted_id), &
            iphase, &
            patch%characteristic_curves_array(int(icap_loc_p(ghosted_id)))%ptr, &
            TH_parameter, ithrm, &
            option)
    else
       call THAuxVarComputeNoFreezing(xx_p(istart:iend), &
            TH_auxvars(ghosted_id),global_auxvars(ghosted_id), &
            material_auxvars(ghosted_id), &
            iphase, &
            patch%characteristic_curves_array(int(icap_loc_p(ghosted_id)))%ptr, &
            TH_parameter, ithrm, &
            option)
    endif


    if (option%use_mc) then
      vol_frac_prim = TH_sec_heat_vars(local_id)%epsilon
    endif
    
    iphase_loc_p(ghosted_id) = iphase
    call THAccumulation(TH_auxvars(ghosted_id),global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              TH_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                              option,vol_frac_prim,accum_p(istart:iend)) 
  enddo

  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc,icap_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

#if 0
   call THNumericalJacobianTest(field%flow_xx,realization)
#endif

end subroutine THUpdateFixedAccumPatch

! ************************************************************************** !

subroutine THNumericalJacobianTest(xx,realization)
  ! 
  ! Computes the a test numerical jacobian
  ! 
  ! Author: ???
  ! Date: 12/13/07
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
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  
  PetscInt :: idof, idof2, icell

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  
  call VecDuplicate(xx,xx_pert,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res_pert,ierr);CHKERRQ(ierr)
  
  call MatCreate(option%mycomm,A,ierr);CHKERRQ(ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,grid%nlmax*option%nflowdof,grid%nlmax*option%nflowdof, &
                   ierr);CHKERRQ(ierr)
  call MatSetType(A,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetFromOptions(A,ierr);CHKERRQ(ierr)
    
  call THResidual(PETSC_NULL_SNES,xx,res,realization,ierr)
  call VecGetArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)
  do icell = 1,grid%nlmax
    if (patch%imat(icell) <= 0) cycle
    do idof = (icell-1)*option%nflowdof+1,icell*option%nflowdof 
      call VecCopy(xx,xx_pert,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      perturbation = vec_p(idof)*perturbation_tolerance
      vec_p(idof) = vec_p(idof)+perturbation
      call VecRestoreArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      call THResidual(PETSC_NULL_SNES,xx_pert,res_pert,realization,ierr)
      call VecGetArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
      do idof2 = 1, grid%nlmax*option%nflowdof
        derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
        if (dabs(derivative) > 1.d-30) then
          call matsetvalue(a,idof2-1,idof-1,derivative,insert_values, &
                           ierr);CHKERRQ(ierr)
        endif
      enddo
      call VecRestoreArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
    enddo
  enddo
  call VecRestoreArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call PetscViewerASCIIOpen(option%mycomm,'numerical_jacobian.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call MatView(A,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  call MatDestroy(A,ierr);CHKERRQ(ierr)
  
  call VecDestroy(xx_pert,ierr);CHKERRQ(ierr)
  call VecDestroy(res,ierr);CHKERRQ(ierr)
  call VecDestroy(res_pert,ierr);CHKERRQ(ierr)
  
end subroutine THNumericalJacobianTest

! ************************************************************************** !

subroutine THAccumDerivative(TH_auxvar,global_auxvar, &
                             material_auxvar, &
                             rock_dencpr, &
                             th_parameter, &
                             ithrm, &
                             option,sat_func, &
                             vol_frac_prim,J)
  ! 
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 

  use Option_module
  use Characteristic_Curves_module
  use Material_Aux_class, only : material_auxvar_type, &
                                 soil_compressibility_index, &
                                 MaterialCompressSoil
  use EOS_Water_module
  
  implicit none

  type(TH_auxvar_type) :: TH_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: vol,por,rock_dencpr
  type(TH_parameter_type) :: th_parameter
  PetscInt :: ithrm
  class(Characteristic_Curves_type), pointer :: sat_func
  PetscReal :: J(option%nflowdof,option%nflowdof)
     
  PetscInt :: ispec 
  PetscReal :: porXvol, mol(option%nflowspec), eng, por1

  PetscInt :: iphase, ideriv
  type(TH_auxvar_type) :: TH_auxvar_pert
  type(global_auxvar_type) :: global_auxvar_pert
  ! leave as type
  type(material_auxvar_type) :: material_auxvar_pert
  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), pert
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: J_pert(option%nflowdof,option%nflowdof)
  PetscReal :: vol_frac_prim, tempreal
  PetscReal :: compressed_porosity, dcompressed_porosity_dp
  
  PetscReal :: pres, temp
  PetscReal :: sat, dsat_dp, dsat_dt
  PetscReal :: den, dden_dp, dden_dt
  PetscReal :: u, du_dp, du_dt

  ! ice variables
  PetscReal :: sat_g, p_g, den_g, p_sat, mol_g, u_g, C_g
  PetscReal :: dpsat_dt
  PetscReal :: ddeng_dt, dmolg_dt, dsatg_dt, dug_dt
  PetscReal :: ddeng_dp, dmolg_dp, dsatg_dp, dug_dp
  PetscReal :: sat_i, den_i, u_i
  PetscReal :: dsati_dp, dsati_dt
  PetscReal :: ddeni_dp, ddeni_dt
  PetscReal :: dui_dp, dui_dt
  PetscErrorCode :: ierr


  J(:,:) = 0.d0
  
  ! X = {p, T}; R = {R_p, R_T}
  
  vol = material_auxvar%volume
  pres = global_auxvar%pres(1)
  temp = global_auxvar%temp
  sat = global_auxvar%sat(1)
  den = global_auxvar%den(1)
  dden_dp = TH_auxvar%dden_dp
  dden_dt = TH_auxvar%dden_dt
  dsat_dp = TH_auxvar%dsat_dp
  dsat_dt = TH_auxvar%dsat_dt
  u = TH_auxvar%u
  du_dt = TH_auxvar%du_dt
  du_dp = TH_auxvar%du_dp
  
  if (soil_compressibility_index > 0) then
    tempreal = sat*den
    call MaterialCompressSoil(material_auxvar,pres, &
                              compressed_porosity,dcompressed_porosity_dp)
    por = compressed_porosity
  else
    por = material_auxvar%porosity_base
    dcompressed_porosity_dp = 0.d0
  endif

  porXvol = por*vol

  ! d(por*sat*den)/dP * vol
  J(TH_PRESSURE_DOF,TH_PRESSURE_DOF) = (sat*dden_dp + dsat_dp*den)*porXvol + &
    dcompressed_porosity_dp*sat*den*vol

  J(TH_PRESSURE_DOF,TH_TEMPERATURE_DOF) = (sat*dden_dt+dsat_dt*den)*porXvol

  !eng = sat*den*u*porXvol + (1.-por)*vol*rock_dencpr*temp
  J(TH_TEMPERATURE_DOF,TH_PRESSURE_DOF) = (dsat_dp*den*u + &
                                           sat*dden_dp*u + &
                                           sat*den*du_dp)*porXvol + &
                        (den*sat*u - rock_dencpr*temp)*vol*dcompressed_porosity_dp

  !eng = sat*den*u*porXvol + (1.-por)*vol*rock_dencpr*temp
  J(TH_TEMPERATURE_DOF,TH_TEMPERATURE_DOF) = (dsat_dt*den*u + &
                                              sat*dden_dt*u + &
                                              sat*den*du_dt )*porXvol +  &
                                             (1.d0 - por)*vol*rock_dencpr


  if (option%use_th_freezing) then
     ! SK, 11/17/11
     sat_g    = TH_auxvar%ice%sat_gas
     sat_i    = TH_auxvar%ice%sat_ice
     dsati_dt = TH_auxvar%ice%dsat_ice_dt
     dsati_dp = TH_auxvar%ice%dsat_ice_dp
     dsatg_dp = TH_auxvar%ice%dsat_gas_dp
     dsatg_dt = TH_auxvar%ice%dsat_gas_dt

     u_g      = TH_auxvar%ice%u_gas
     u_i      = TH_auxvar%ice%u_ice
     dug_dt   = TH_auxvar%ice%du_gas_dt
     dug_dp   = TH_auxvar%ice%du_gas_dp
     dui_dt   = TH_auxvar%ice%du_ice_dt
     dui_dp   = TH_auxvar%ice%du_ice_dp

     den_i    = TH_auxvar%ice%den_ice
     den_g    = TH_auxvar%ice%den_gas
     mol_g    = TH_auxvar%ice%mol_gas
     ddeni_dt = TH_auxvar%ice%dden_ice_dt
     ddeni_dp = TH_auxvar%ice%dden_ice_dp
     ddeng_dt = TH_auxvar%ice%dden_gas_dt
     ddeng_dp = TH_auxvar%ice%dden_gas_dp
     dmolg_dt = TH_auxvar%ice%dmol_gas_dt
     dmolg_dp = TH_auxvar%ice%dmol_gas_dp

     !mol(1) = mol(1) + (sat_g*den_g*mol_g + sat_i*den_i)*porXvol
     J(TH_PRESSURE_DOF,TH_PRESSURE_DOF) = J(TH_PRESSURE_DOF,TH_PRESSURE_DOF) + &
                                          (dsatg_dp*den_g*mol_g     + &
                                           sat_g   *den_g*dmolg_dp  + &
                                           sat_g   *mol_g*ddeng_dp  + &
                                           dsati_dp*den_i       + &
                                           sat_i   *ddeni_dp     )*porXvol + &
                                          (sat_g   *den_g*mol_g + &
                                           sat_i   *den_i        )*dcompressed_porosity_dp*vol

     J(TH_PRESSURE_DOF,TH_TEMPERATURE_DOF) = J(TH_PRESSURE_DOF,TH_TEMPERATURE_DOF) + &
                            (dsatg_dt * den_g    * mol_g            + &
                             sat_g    * ddeng_dt * mol_g            + &
                             sat_g    * den_g    * dmolg_dt         + &
                             dsati_dt * den_i                       + &
                             sat_i    * ddeni_dt                    )*porXvol

     !eng = eng + (sat_g*den_g*u_g + sat_i*den_i*u_i)*porXvol
     J(TH_TEMPERATURE_DOF,TH_PRESSURE_DOF) = J(TH_TEMPERATURE_DOF,TH_PRESSURE_DOF) + &
                     (dsatg_dp * den_g    * u_g                + &
                      sat_g    *(ddeng_dp*u_g+den_g*dug_dp)    + &
                      dsati_dp * den_i    * u_i                + &
                      sat_i    *(ddeni_dp*u_i+dui_dp*den_i) ) *porXvol   + &
                     (sat_g    * den_g    * u_g                + &
                      sat_i    * den_i    * u_i )*dcompressed_porosity_dp*vol

     J(TH_TEMPERATURE_DOF,TH_TEMPERATURE_DOF) = J(TH_TEMPERATURE_DOF,TH_TEMPERATURE_DOF) + &
                (dsatg_dt * den_g   * u_g                           + &
                 sat_g  * (ddeng_dt * u_g+ den_g *dug_dt)           + &
                  dsati_dt * den_i    * u_i                         + &
                  sat_i    * ddeni_dt * u_i                         + &
                  sat_i    * den_i    * dui_dt                      )*porXvol

  endif

  J = J/option%flow_dt
  J(option%nflowdof,:) = vol_frac_prim*J(option%nflowdof,:)

  ! If only solving the energy equation,
  !  - Set jacobian term corresponding to mass-equation to zero, and
  !  - Set off-diagonal jacobian terms to zero.
  if (option%flow%only_energy_eq) then
    J(1,1) = 1.d0
    J(1,2) = 0.d0
    J(2,1) = 0.d0
  endif

  if (option%flow%numerical_derivatives) then
    call GlobalAuxVarInit(global_auxvar_pert,option)  
    call MaterialAuxVarInit(material_auxvar_pert,option)  

    call THAuxVarCopy(TH_auxvar,TH_auxvar_pert,option)
    call GlobalAuxVarCopy(global_auxvar,global_auxvar_pert,option)
    call MaterialAuxVarCopy(material_auxvar,material_auxvar_pert,option)

    x(1) = global_auxvar%pres(1)
    x(2) = global_auxvar%temp
    
    call THAccumulation(TH_auxvar,global_auxvar,material_auxvar, &
                         rock_dencpr,option, &
                         vol_frac_prim,res)
    
    do ideriv = 1,option%nflowdof
      pert = x(ideriv)*perturbation_tolerance
      x_pert = x
      if (option%use_th_freezing) then
         if (ideriv == 1) then
            if (x_pert(ideriv) < option%reference_pressure) then
               pert = - pert
            endif
            x_pert(ideriv) = x_pert(ideriv) + pert
         endif
      
         if (ideriv == 2) then
            if (x_pert(ideriv) < 0.d0) then
               pert = - 1.d-8
            else
               pert =  1.d-8
            endif
            x_pert(ideriv) = x_pert(ideriv) + pert
         endif
      else
         x_pert(ideriv) = x_pert(ideriv) + pert
      endif

      if (option%use_th_freezing) then
         call THAuxVarComputeFreezing(x_pert,TH_auxvar_pert, &
                                 global_auxvar_pert,material_auxvar_pert, &
                                 iphase,sat_func, &
                                 TH_parameter, ithrm, &
                                 option)
      else
         call THAuxVarComputeNoFreezing(x_pert,TH_auxvar_pert,&
                              global_auxvar_pert,material_auxvar_pert,&
                              iphase,sat_func, &
                              TH_parameter,ithrm, &
                              option)
      endif

      call THAccumulation(TH_auxvar_pert,global_auxvar_pert, material_auxvar_pert, &
                           rock_dencpr,option,vol_frac_prim, &
                           res_pert)
      J_pert(:,ideriv) = (res_pert(:)-res(:))/pert
    enddo

    J = J_pert
    call GlobalAuxVarStrip(global_auxvar_pert)  
  endif
   
end subroutine THAccumDerivative

! ************************************************************************** !

subroutine THAccumulation(auxvar,global_auxvar, &
                          material_auxvar, &
                          rock_dencpr,option,vol_frac_prim,Res)
  ! 
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 

  use Option_module
  use Material_Aux_class, only : material_auxvar_type, &
                                 soil_compressibility_index, &
                                 MaterialCompressSoil
  use EOS_Water_module
  
  implicit none

  type(TH_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal ::rock_dencpr,por1

  PetscInt :: ispec
  PetscReal :: vol,por,uh
  PetscReal :: porXvol, mol(option%nflowspec), eng
  PetscReal :: vol_frac_prim
  PetscReal :: compressed_porosity, dcompressed_porosity_dp

  ! ice variables
  PetscReal :: sat_g, p_g, den_g, p_sat, mol_g, u_g, C_g
  PetscReal :: sat_i, den_i, u_i
  PetscErrorCode :: ierr
  
  vol = material_auxvar%volume
  
  if (soil_compressibility_index > 0) then
    call MaterialCompressSoil(material_auxvar,global_auxvar%pres(1), &
                              compressed_porosity,dcompressed_porosity_dp)
    material_auxvar%porosity = compressed_porosity
    por = compressed_porosity
  else
    por = material_auxvar%porosity_base
  endif
  auxvar%transient_por = por

  ! TechNotes, TH Mode: First term of Equation 8
  porXvol = por*vol

  mol(1) = global_auxvar%sat(1)*global_auxvar%den(1)*porXvol
  uh     = auxvar%u

! TechNotes, TH Mode: First term of Equation 9
  ! rock_dencpr [MJ/m^3 rock-K]
  eng = global_auxvar%sat(1) * &
        global_auxvar%den(1) * &
        uh * porXvol +         &
        (1.d0 - por) * vol * rock_dencpr * global_auxvar%temp

  if (option%use_th_freezing) then
     ! SK, 11/17/11
     sat_g = auxvar%ice%sat_gas
     sat_i = auxvar%ice%sat_ice
     u_i   = auxvar%ice%u_ice
     den_i = auxvar%ice%den_ice
     den_g = auxvar%ice%den_gas
     mol_g = auxvar%ice%mol_gas
     u_g = auxvar%ice%u_gas
     mol(1) = mol(1) + (sat_g*den_g*mol_g + sat_i*den_i)*porXvol
     eng = eng + (sat_g*den_g*u_g + sat_i*den_i*u_i)*porXvol
  endif

  Res(1:option%nflowdof-1) = mol(:)/option%flow_dt
  Res(option%nflowdof) = vol_frac_prim*eng/option%flow_dt

end subroutine THAccumulation

! ************************************************************************** !

subroutine THFluxDerivative(auxvar_up,global_auxvar_up, &
                            material_auxvar_up, &
                            sir_up, &
                            Dk_up, &
                            ithrm_up, &
                            auxvar_dn,global_auxvar_dn, &
                            material_auxvar_dn, &
                            sir_dn, &
                            Dk_dn, &
                            ithrm_dn, &
                            area, &
                            dist, upweight, &
                            option,sat_func_up,sat_func_dn, &
                            Dk_dry_up,Dk_dry_dn, &
                            Dk_ice_up,Dk_ice_dn, &
                            alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                            th_parameter, &
                            Jup,Jdn)

  !
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 
                             
  use Option_module 
  use Characteristic_Curves_module
  use Connection_module
  use EOS_Water_module
  use Utility_module
  
  implicit none
  
  type(TH_auxvar_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: Dk_dry_up, Dk_dry_dn
  PetscReal :: Dk_ice_up, Dk_ice_dn
  PetscReal :: alpha_up, alpha_dn
  PetscReal :: alpha_fr_up, alpha_fr_dn
  PetscInt :: ithrm_up, ithrm_dn
  PetscReal :: v_darcy, area
  PetscReal :: dist(-1:3)
  class(Characteristic_Curves_type), pointer :: sat_func_up, sat_func_dn
  type(TH_parameter_type) :: th_parameter
  PetscReal :: Jup(option%nflowdof,option%nflowdof), Jdn(option%nflowdof,option%nflowdof)
     
  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscInt :: ispec
  PetscReal :: fluxm,fluxe,q
  PetscReal :: uh,ukvr,DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
  
  PetscReal :: dden_ave_dp_up, dden_ave_dp_dn, dden_ave_dt_up, dden_ave_dt_dn
  PetscReal :: dgravity_dden_up, dgravity_dden_dn
  PetscReal :: dphi_dp_up, dphi_dp_dn, dphi_dt_up, dphi_dt_dn
  PetscReal :: dukvr_dp_up, dukvr_dp_dn, dukvr_dt_up, dukvr_dt_dn
  PetscReal :: duh_dp_up, duh_dp_dn, duh_dt_up, duh_dt_dn
  PetscReal :: dq_dp_up, dq_dp_dn, dq_dt_up, dq_dt_dn
  
  PetscReal :: Dk_eff_up, Dk_eff_dn
  PetscReal, parameter :: epsilon = 1.d-6
  PetscReal :: dKe_dt_up, dKe_dp_up
  PetscReal :: dKe_dt_dn, dKe_dp_dn
  PetscReal :: dDk_dt_up, dDk_dt_dn
  PetscReal :: dDk_dp_up, dDk_dp_dn

  PetscInt :: iphase, ideriv
  type(TH_auxvar_type) :: auxvar_pert_up, auxvar_pert_dn
  type(global_auxvar_type) :: global_auxvar_pert_up, global_auxvar_pert_dn
  PetscReal :: x_up(option%nflowdof), x_dn(option%nflowdof)
  PetscReal :: x_pert_up(option%nflowdof), x_pert_dn(option%nflowdof)
  PetscReal :: pert_up, pert_dn
  PetscReal :: res(option%nflowdof)
  PetscReal :: res_pert_up(option%nflowdof)
  PetscReal :: res_pert_dn(option%nflowdof)
  PetscReal :: J_pert_up(option%nflowdof,option%nflowdof)
  PetscReal :: J_pert_dn(option%nflowdof,option%nflowdof)
  class(material_auxvar_type), allocatable :: material_auxvar_pert_dn, &
                                              material_auxvar_pert_up

  ! ice variables
  PetscReal :: Ddiffgas_avg, Ddiffgas_up, Ddiffgas_dn
  PetscReal :: p_g
  PetscReal :: deng_up, deng_dn
  PetscReal :: psat_up, psat_dn
  PetscReal :: molg_up, molg_dn
  PetscReal :: satg_up, satg_dn
  PetscReal :: Diffg_up, Diffg_dn
  PetscReal :: ddeng_dt_up, ddeng_dt_dn, ddeng_dp_up, ddeng_dp_dn
  PetscReal :: dpsat_dt_up, dpsat_dt_dn
  PetscReal :: dmolg_dt_up, dmolg_dt_dn
  PetscReal :: dDiffg_dt_up, dDiffg_dt_dn
  PetscReal :: dDiffg_dp_up, dDiffg_dp_dn
  PetscReal :: dsatg_dp_up, dsatg_dp_dn, dsatg_dt_up, dsatg_dt_dn
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscErrorCode :: ierr
  PetscReal :: Ke_fr_up,Ke_fr_dn   ! frozen soil Kersten numbers
  PetscReal :: dKe_fr_dt_up, dKe_fr_dt_dn
  PetscReal :: dKe_fr_dp_up, dKe_fr_dp_dn
  PetscReal :: fv_up, fv_dn
  PetscReal :: dfv_dt_up, dfv_dt_dn
  PetscReal :: dfv_dp_up, dfv_dp_dn
  PetscReal :: dmolg_dp_up, dmolg_dp_dn
  PetscReal :: fv_up_pert
  PetscReal :: ugas_ave, dugas_ave_dt, dugas_ave_dp, fdiffgas, fdiffgas_dx

  !
  call ConnectionCalculateDistances(dist,option%gravity,dd_up,dd_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  por_up = material_auxvar_up%porosity
  por_dn = material_auxvar_dn%porosity

  tor_up = material_auxvar_up%tortuosity
  tor_dn = material_auxvar_dn%tortuosity

  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)

  fluxm = 0.D0
  fluxe = 0.D0
  v_darcy = 0.D0 
  
  Jup = 0.d0
  Jdn = 0.d0 
  
  dden_ave_dp_up = 0.d0
  dden_ave_dt_up = 0.d0
  dden_ave_dp_dn = 0.d0
  dden_ave_dt_dn = 0.d0
  dgravity_dden_up = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_up = 0.d0
  dphi_dp_dn = 0.d0
  dphi_dt_up = 0.d0
  dphi_dt_dn = 0.d0
  dukvr_dp_up = 0.d0
  dukvr_dp_dn = 0.d0
  dukvr_dt_up = 0.d0
  dukvr_dt_dn = 0.d0
  duh_dp_up = 0.d0
  duh_dp_dn = 0.d0
  duh_dt_up = 0.d0
  duh_dt_dn = 0.d0
  dq_dp_up = 0.d0
  dq_dp_dn = 0.d0
  dq_dt_up = 0.d0
  dq_dt_dn = 0.d0
  dDk_dt_up = 0.d0
  dDk_dt_dn = 0.d0
  dDk_dp_up = 0.d0
  dDk_dp_dn = 0.d0
  
  dfv_dt_up = 0.d0
  dfv_dt_dn = 0.d0
  dfv_dp_up = 0.d0
  dfv_dp_dn = 0.d0
  dmolg_dp_up = 0.d0
  dmolg_dp_dn = 0.d0
  dmolg_dt_up = 0.d0
  dmolg_dt_dn = 0.d0

  ! Flow term
  if (global_auxvar_up%sat(1) > sir_up .or.  &
      global_auxvar_dn%sat(1) > sir_dn) then
    if (global_auxvar_up%sat(1) <eps) then
      upweight=0.d0
    else if (global_auxvar_dn%sat(1) <eps) then 
      upweight=1.d0
    endif
    density_ave = upweight*global_auxvar_up%den(1)+ &
                  (1.D0-upweight)*global_auxvar_dn%den(1)
    dden_ave_dp_up = upweight*auxvar_up%dden_dp
    dden_ave_dp_dn = (1.D0-upweight)*auxvar_dn%dden_dp
    dden_ave_dt_up = upweight*auxvar_up%dden_dt
    dden_ave_dt_dn = (1.D0-upweight)*auxvar_dn%dden_dt

    gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
              (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
              * dist_gravity
    dgravity_dden_up = upweight*auxvar_up%avgmw*dist_gravity
    dgravity_dden_dn = (1.d0-upweight)*auxvar_dn%avgmw*dist_gravity

    if (option%ice_model /= DALL_AMICO &
        .and.option%ice_model /= PAINTER_KARRA_EXPLICIT_SMOOTH) then
      dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
      dphi_dp_up = 1.d0 + dgravity_dden_up*auxvar_up%dden_dp
      dphi_dp_dn = -1.d0 + dgravity_dden_dn*auxvar_dn%dden_dp
      dphi_dt_up = dgravity_dden_up*auxvar_up%dden_dt
      dphi_dt_dn = dgravity_dden_dn*auxvar_dn%dden_dt

    else
      dphi = auxvar_up%ice%pres_fh2o - auxvar_dn%ice%pres_fh2o + gravity
      dphi_dp_up =  auxvar_up%ice%dpres_fh2o_dp + &
                    dgravity_dden_up*auxvar_up%dden_dp
      dphi_dp_dn = -auxvar_dn%ice%dpres_fh2o_dp + &
                    dgravity_dden_dn*auxvar_dn%dden_dp
      dphi_dt_up =  auxvar_up%ice%dpres_fh2o_dt + &
                    dgravity_dden_up*auxvar_up%dden_dt
      dphi_dt_dn = -auxvar_dn%ice%dpres_fh2o_dt + &
                    dgravity_dden_dn*auxvar_dn%dden_dt
    endif

    if (dphi>=0.D0) then
      ukvr = auxvar_up%kvr
      dukvr_dp_up = auxvar_up%dkvr_dp
      dukvr_dt_up = auxvar_up%dkvr_dt
      
      uh = auxvar_up%h
      duh_dp_up = auxvar_up%dh_dp
      duh_dt_up = auxvar_up%dh_dt
    else
      ukvr = auxvar_dn%kvr
      dukvr_dp_dn = auxvar_dn%dkvr_dp
      dukvr_dt_dn = auxvar_dn%dkvr_dt
      
      uh = auxvar_dn%h
      duh_dp_dn = auxvar_dn%dh_dp
      duh_dt_dn = auxvar_dn%dh_dt
    endif      

    call InterfaceApprox(auxvar_up%kvr, auxvar_dn%kvr, &
                         auxvar_up%dkvr_dp, auxvar_dn%dkvr_dp, &
                         dphi, &
                         option%rel_perm_aveg, &
                         ukvr, dukvr_dp_up, dukvr_dp_dn)
    call InterfaceApprox(auxvar_up%kvr, auxvar_dn%kvr, &
                         auxvar_up%dkvr_dt, auxvar_dn%dkvr_dt, &
                         dphi, &
                         option%rel_perm_aveg, &
                         ukvr, dukvr_dt_up, dukvr_dt_dn)

    !call InterfaceApprox(auxvar_up%h, auxvar_dn%h, &
    !                     auxvar_up%dh_dp, auxvar_dn%dh_dp, &
    !                     dphi, &
    !                     option%rel_perm_aveg, &
    !                     uh, duh_dp_up, duh_dp_dn)
    !call InterfaceApprox(auxvar_up%h, auxvar_dn%h, &
    !                     auxvar_up%dh_dt, auxvar_dn%dh_dt, &
    !                     dphi, &
    !                     option%rel_perm_aveg, &
    !                     uh, duh_dt_up, duh_dp_dn)

    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi
   
      q = v_darcy * area
      dq_dp_up = Dq*(dukvr_dp_up*dphi+ukvr*dphi_dp_up)*area
      dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
      
      dq_dt_up = Dq*(dukvr_dt_up*dphi+ukvr*dphi_dt_up)*area
      dq_dt_dn = Dq*(dukvr_dt_dn*dphi+ukvr*dphi_dt_dn)*area
        

      ! If only solving the energy equation, ensure Jup(2,2) & Jdn(2,2)
      ! have no contribution from the mass equation
      if (option%flow%only_energy_eq) then
         v_darcy = 0.d0
         q = 0.d0
         dq_dt_up = 0.d0
         dq_dt_dn = 0.d0
      endif

      Jup(1,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)
      Jup(1,2) = (dq_dt_up*density_ave+q*dden_ave_dt_up)

      Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)
      Jdn(1,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)

      ! based on flux = q*density_ave*uh
      Jup(option%nflowdof,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)*uh+ &
                               q*density_ave*duh_dp_up
      Jup(option%nflowdof,2) = (dq_dt_up*density_ave+q*dden_ave_dt_up)*uh+ &
                               q*density_ave*duh_dt_up
      Jdn(option%nflowdof,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uh+ &
                               q*density_ave*duh_dp_dn
      Jdn(option%nflowdof,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uh+ &
                               q*density_ave*duh_dt_dn

    endif
  endif 

  if (option%use_th_freezing) then
    ! Added by Satish Karra, updated 11/11/11
    satg_up = auxvar_up%ice%sat_gas
    satg_dn = auxvar_dn%ice%sat_gas
    if ((satg_up > eps) .and. (satg_dn > eps)) then

      p_g = option%reference_pressure  ! set to reference pressure
      deng_up = auxvar_up%ice%den_gas
      deng_dn = auxvar_dn%ice%den_gas

      Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
      p_ref = 1.01325d5   ! in Pa
      T_ref = 25.d0       ! in deg C

      Diffg_up = Diffg_ref*(p_ref/p_g)*((max(-50d0, global_auxvar_up%temp) + 273.15d0) &
                 /(T_ref + 273.15d0))**(1.8d0)
      Diffg_dn = Diffg_ref*(p_ref/p_g)*((max(-50d0, global_auxvar_dn%temp) + 273.15d0) &
                 /(T_ref + 273.15d0))**(1.8d0)

      Ddiffgas_up = por_up*tor_up*satg_up*Diffg_up
      Ddiffgas_dn = por_dn*tor_dn*satg_dn*Diffg_dn

#if 0
      ! vapor pressure lowering due to capillary pressure
      fv_up = exp(-auxvar_up%pc/(global_auxvar_up%den(1)* &
           IDEAL_GAS_CONSTANT*(global_auxvar_up%temp + 273.15d0)))
      fv_dn = exp(-auxvar_dn%pc/(global_auxvar_dn%den(1)* &
           IDEAL_GAS_CONSTANT*(global_auxvar_dn%temp + 273.15d0)))

      dfv_dt_up = fv_up*(auxvar_up%pc/IDEAL_GAS_CONSTANT/(global_auxvar_up%den(1)* &
           (global_auxvar_up%temp + 273.15d0))**2)* &
           (auxvar_up%dden_dt*(global_auxvar_up%temp + 273.15d0) &
           + global_auxvar_up%den(1))
      dfv_dt_dn = fv_dn*(auxvar_dn%pc/IDEAL_GAS_CONSTANT/(global_auxvar_dn%den(1)* &
           (global_auxvar_dn%temp + 273.15d0))**2)* &
           (auxvar_dn%dden_dt*(global_auxvar_dn%temp + 273.15d0) &
           + global_auxvar_dn%den(1))

      dfv_dp_up = fv_up*(auxvar_up%pc/IDEAL_GAS_CONSTANT/(global_auxvar_up%den(1))**2/ &
           (global_auxvar_up%temp + 273.15d0)*auxvar_up%dden_dp &
           + 1.d0/IDEAL_GAS_CONSTANT/global_auxvar_up%den(1)/ &
           (global_auxvar_up%temp + 273.15d0))
      dfv_dp_dn = fv_dn*(auxvar_dn%pc/IDEAL_GAS_CONSTANT/(global_auxvar_dn%den(1))**2/ &
           (global_auxvar_dn%temp + 273.15d0)*auxvar_dn%dden_dp &
           + 1.d0/IDEAL_GAS_CONSTANT/global_auxvar_dn%den(1)/ &
           (global_auxvar_dn%temp + 273.15d0))

      molg_up = auxvar_up%ice%mol_gas*fv_up
      molg_dn = auxvar_dn%ice%mol_gas*fv_dn
      dmolg_dt_up = auxvar_up%ice%dmol_gas_dt*fv_up + auxvar_up%ice%mol_gas*dfv_dt_up
      dmolg_dt_dn = auxvar_dn%ice%dmol_gas_dt*fv_dn + auxvar_dn%ice%mol_gas*dfv_dt_dn
      dmolg_dp_up = auxvar_up%ice%dmol_gas_dp*fv_up + auxvar_up%ice%mol_gas*dfv_dp_up
      dmolg_dp_dn = auxvar_dn%ice%dmol_gas_dp*fv_dn + auxvar_dn%ice%mol_gas*dfv_dp_dn

#else
! seems having issues with vapor-only diffusion as above (not sure why)
      fv_up = 1.d0
      fv_dn = 1.d0
      dfv_dt_up = 0.d0
      dfv_dt_dn = 0.d0
      dfv_dp_up = 0.d0
      dfv_dp_dn = 0.d0

      molg_up = auxvar_up%ice%mol_gas
      molg_dn = auxvar_dn%ice%mol_gas
      dmolg_dt_up = auxvar_up%ice%dmol_gas_dt
      dmolg_dt_dn = auxvar_dn%ice%dmol_gas_dt
      dmolg_dp_up = auxvar_up%ice%dmol_gas_dp
      dmolg_dp_dn = auxvar_dn%ice%dmol_gas_dp

#endif

      ddeng_dt_up = auxvar_up%ice%dden_gas_dt
      ddeng_dt_dn = auxvar_dn%ice%dden_gas_dt
      ddeng_dp_up = auxvar_up%ice%dden_gas_dp
      ddeng_dp_dn = auxvar_dn%ice%dden_gas_dp

      dDiffg_dt_up = 1.8d0*Diffg_up/(max(-50.0d0,global_auxvar_up%temp) + 273.15d0)
      dDiffg_dt_dn = 1.8d0*Diffg_dn/(max(-50.0d0,global_auxvar_dn%temp) + 273.15d0)
      dDiffg_dp_up = 0.d0
      dDiffg_dp_dn = 0.d0

      dsatg_dp_up = auxvar_up%ice%dsat_gas_dp
      dsatg_dp_dn = auxvar_dn%ice%dsat_gas_dp
      dsatg_dt_up = auxvar_up%ice%dsat_gas_dt
      dsatg_dt_dn = auxvar_dn%ice%dsat_gas_dt
 
      if (deng_up*molg_up > deng_dn*molg_dn) then
      ! fmyuan: 'molg' is mole fraction of vapor in air, NOT vapor density
        upweight = 0.d0
        ugas_ave = auxvar_up%ice%u_gas
        dugas_ave_dt = auxvar_up%ice%du_gas_dt
        dugas_ave_dp = auxvar_up%ice%du_gas_dp
      else
        upweight = 1.d0
        ugas_ave = auxvar_dn%ice%u_gas
        dugas_ave_dt = auxvar_dn%ice%du_gas_dt
        dugas_ave_dp = auxvar_dn%ice%du_gas_dp
      endif

      Ddiffgas_avg = upweight*Ddiffgas_up + (1.D0 - upweight)*Ddiffgas_dn

      ! derivaives: Ddiffgas_avg*area*(deng_up*molg_up - deng_dn*molg_dn)/(dd_up + dd_dn)
      fdiffgas = area/(dd_up + dd_dn)*Ddiffgas_avg*(deng_up*molg_up - deng_dn*molg_dn)

      fdiffgas_dx = area/(dd_up+dd_dn) * &                                             !_up, _dp
                 ( Ddiffgas_avg * (ddeng_dp_up*molg_up + deng_up*dmolg_dp_up) + &
                   (deng_up*molg_up - deng_dn*molg_dn) * upweight * por_up*tor_up* &
                    (satg_up*dDiffg_dp_up + dsatg_dp_up*Diffg_up) )
      Jup(1,1) = Jup(1,1) + fdiffgas_dx
      Jup(2,1) = Jup(2,1) + (fdiffgas*dugas_ave_dp + ugas_ave*fdiffgas_dx)


      fdiffgas_dx = area/(dd_up+dd_dn) * &                                             !_up, _dt
                 ( Ddiffgas_avg * (ddeng_dt_up*molg_up + deng_up*dmolg_dt_up) + &
                   (deng_up*molg_up - deng_dn*molg_dn) * upweight * por_up*tor_up* &
                    (satg_up*dDiffg_dt_up + dsatg_dt_up*Diffg_up) )
      Jup(1,2) = Jup(1,2) + fdiffgas_dx
      Jup(2,2) = Jup(2,2) + (fdiffgas*dugas_ave_dt + ugas_ave*fdiffgas_dx)


      fdiffgas_dx = area/(dd_up+dd_dn) * &                                             !_dn, _dp
                 ( Ddiffgas_avg * (-ddeng_dp_dn*molg_dn - deng_dn*dmolg_dp_dn) + &
                   (deng_up*molg_up - deng_dn*molg_dn) * (1.d0-upweight) * por_dn*tor_dn* &
                    (satg_dn*dDiffg_dp_dn + dsatg_dp_dn*Diffg_dn) )
      Jdn(1,1) = Jdn(1,1) + fdiffgas_dx
      Jdn(2,1) = Jdn(2,1) + (fdiffgas*dugas_ave_dp + ugas_ave*fdiffgas_dx)


      fdiffgas_dx = area/(dd_up+dd_dn) * &                                             !_dn, _dt
                 ( Ddiffgas_avg * (-ddeng_dt_dn*molg_dn - deng_dn*dmolg_dt_dn) + &
                   (deng_up*molg_up - deng_dn*molg_dn) * (1.d0-upweight) * por_dn*tor_dn* &
                    (satg_dn*dDiffg_dt_dn + dsatg_dt_dn*Diffg_dn) )
      Jdn(1,2) = Jdn(1,2) + fdiffgas_dx
      Jdn(2,2) = Jdn(2,2) + (fdiffgas*dugas_ave_dt + ugas_ave*fdiffgas_dx)


   endif

  endif ! if (use_th_freezing)

        
  dKe_dp_up = auxvar_up%dKe_dp
  dKe_dp_dn = auxvar_dn%dKe_dp

  dKe_dt_up = auxvar_up%dKe_dt
  dKe_dt_dn = auxvar_dn%dKe_dt

  if (option%use_th_freezing) then
            
    Dk_eff_up = auxvar_up%Dk_eff
    Dk_eff_dn = auxvar_dn%Dk_eff

    Ke_fr_up = auxvar_up%ice%Ke_fr
    Ke_fr_dn = auxvar_dn%ice%Ke_fr

    dKe_fr_dt_up = auxvar_up%ice%dKe_fr_dt
    dKe_fr_dt_dn = auxvar_dn%ice%dKe_fr_dt

    dKe_fr_dp_up = auxvar_up%ice%dKe_fr_dp
    dKe_fr_dp_dn = auxvar_dn%ice%dKe_fr_dp

  else

    Dk_eff_up = auxvar_up%Dk_eff
    Dk_eff_dn = auxvar_dn%Dk_eff

  endif
 
  Dk = (Dk_eff_up * Dk_eff_dn) / (dd_dn*Dk_eff_up + dd_up*Dk_eff_dn)
  
  if (option%use_th_freezing) then

    ! 1/Dk = dd_dn/Dk_eff_dn + dd_up/Dk_eff_up
    dKe_dp_up = auxvar_up%dDk_eff_dp  ! effective one
    dKe_dp_dn = auxvar_dn%dDk_eff_dp
    dKe_dt_up = auxvar_up%dDk_eff_dt
    dKe_dt_dn = auxvar_dn%dDk_eff_dt

    dDk_dp_up = Dk*Dk * (dd_up/Dk_eff_up/Dk_eff_up*dKe_dp_up)
    dDk_dp_dn = Dk*Dk * (dd_dn/Dk_eff_dn/Dk_eff_dn*dKe_dp_dn)

    dDk_dt_up = Dk*Dk * (dd_up/Dk_eff_up/Dk_eff_up*dKe_dt_up)
    dDk_dt_dn = Dk*Dk * (dd_dn/Dk_eff_dn/Dk_eff_dn*dKe_dt_dn)

  else
  
    dDk_dt_up = Dk**2/Dk_eff_up**2*dd_up*(Dk_up - Dk_dry_up)*dKe_dt_up
    dDk_dt_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn - Dk_dry_dn)*dKe_dt_dn

    dDk_dp_up = Dk**2/Dk_eff_up**2*dd_up*(Dk_up - Dk_dry_up)*dKe_dp_up
    dDk_dp_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn - Dk_dry_dn)*dKe_dp_dn

  endif
    
  !  cond = Dk*area*(global_auxvar_up%temp-global_auxvar_dn%temp) 
  Jup(option%nflowdof,1) = Jup(option%nflowdof,1) + &
                           area*(global_auxvar_up%temp - &
                           global_auxvar_dn%temp)*dDk_dp_up
  Jdn(option%nflowdof,1) = Jdn(option%nflowdof,1) + &
                           area*(global_auxvar_up%temp - &
                           global_auxvar_dn%temp)*dDk_dp_dn
                           
  Jup(option%nflowdof,2) = Jup(option%nflowdof,2) + Dk*area + &
                           area*(global_auxvar_up%temp - & 
                           global_auxvar_dn%temp)*dDk_dt_up 
  Jdn(option%nflowdof,2) = Jdn(option%nflowdof,2) + Dk*area*(-1.d0) + &
                           area*(global_auxvar_up%temp - & 
                           global_auxvar_dn%temp)*dDk_dt_dn 

  ! If only solving the energy equation,
  !  - Set jacobian term corresponding to mass-equation to zero, and
  !  - Set off-diagonal jacobian terms to zero.
  if (option%flow%only_energy_eq) then
    Jup(1,1) = 0.d0
    Jup(1,2) = 0.d0
    Jup(option%nflowdof,1) = 0.d0

    Jdn(1,1) = 0.d0
    Jdn(1,2) = 0.d0
    Jdn(option%nflowdof,1) = 0.d0

  endif

  ! note: Res is the flux contribution, for node up J = J + Jup
  !                                              dn J = J - Jdn  

  if (option%flow%numerical_derivatives) then
    call THAuxVarCopy(auxvar_up,auxvar_pert_up,option)
    call THAuxVarCopy(auxvar_dn,auxvar_pert_dn,option)

    call GlobalAuxVarInit(global_auxvar_pert_up,option)
    call GlobalAuxVarInit(global_auxvar_pert_dn,option)  
    call GlobalAuxVarCopy(global_auxvar_up,global_auxvar_pert_up,option)
    call GlobalAuxVarCopy(global_auxvar_dn,global_auxvar_pert_dn,option)

    allocate(material_auxvar_pert_up,material_auxvar_pert_dn)
    call MaterialAuxVarInit(material_auxvar_pert_up,option)
    call MaterialAuxVarInit(material_auxvar_pert_dn,option)  
    call MaterialAuxVarCopy(material_auxvar_up,material_auxvar_pert_up,option)
    call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_dn,option)

    x_up(1) = global_auxvar_up%pres(1)
    x_up(2) = global_auxvar_up%temp
    x_dn(1) = global_auxvar_dn%pres(1)
    x_dn(2) = global_auxvar_dn%temp

    call THFlux( &
      auxvar_up,global_auxvar_up, &
      material_auxvar_up, &
      sir_up, &
      Dk_up, &
      auxvar_dn,global_auxvar_dn, &
      material_auxvar_dn, &
      sir_dn, &
      Dk_dn, &
      area, &
      dist, upweight, &
      option,v_darcy,Dk_dry_up,Dk_dry_dn, &
      Dk_ice_up,Dk_ice_dn, &
      alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
      res)

    do ideriv = 1,option%nflowdof
      pert_up = x_up(ideriv)*perturbation_tolerance
      pert_dn = x_dn(ideriv)*perturbation_tolerance
      x_pert_up = x_up
      x_pert_dn = x_dn

      if (option%use_th_freezing) then
        if (ideriv == 1) then
          if (x_pert_up(ideriv) < option%reference_pressure) then
            pert_up = - pert_up
          endif
          x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
          
          if (x_pert_dn(ideriv) < option%reference_pressure) then
            pert_dn = - pert_dn
          endif
          x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
        endif

        if (ideriv == 2) then
          if (x_pert_up(ideriv) < 0.d0) then
            pert_up = - 1.d-5
          else
            pert_up = 1.d-5
          endif
          x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up

          if (x_pert_dn(ideriv) < 0.d0) then
            pert_dn = - 1.d-5
          else
            pert_dn = 1.d-5
          endif
          x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
        endif

      else
         x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
         x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn

      endif

      if (option%use_th_freezing) then
        call THAuxVarComputeFreezing(x_pert_up,auxvar_pert_up, &
             global_auxvar_pert_up, material_auxvar_pert_up, &
             iphase,sat_func_up, &
             TH_parameter,ithrm_up, &
             option)
        call THAuxVarComputeFreezing(x_pert_dn,auxvar_pert_dn, &
             global_auxvar_pert_dn, material_auxvar_pert_up, &
             iphase,sat_func_dn, &
             TH_parameter,ithrm_up, &
             option)
      else
        call THAuxVarComputeNoFreezing(x_pert_up,auxvar_pert_up, &
             global_auxvar_pert_up, material_auxvar_pert_up, &
             iphase,sat_func_up, &
             TH_parameter,ithrm_up, &
             option)
        call THAuxVarComputeNoFreezing(x_pert_dn,auxvar_pert_dn, &
             global_auxvar_pert_dn,material_auxvar_pert_dn, &
             iphase,sat_func_dn, &
             TH_parameter,ithrm_dn, &
             option)
      endif

      call THFlux(auxvar_pert_up,global_auxvar_pert_up, &
                   material_auxvar_pert_up, &
                   sir_up, &
                   Dk_up, &
                   auxvar_dn,global_auxvar_dn, &
                   material_auxvar_pert_dn, &
                   sir_dn, &
                   Dk_dn, &
                   area, &
                   dist, upweight, &
                   option,v_darcy,Dk_dry_up, &
                   Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                   alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                   res_pert_up)
      call THFlux(auxvar_up,global_auxvar_up, &
                   material_auxvar_pert_up, &
                   sir_up,&
                   Dk_up, &
                   auxvar_pert_dn,global_auxvar_pert_dn, &
                   material_auxvar_pert_dn, &
                   sir_dn, &
                   Dk_dn, &
                   area, &
                   dist, upweight, &
                   option,v_darcy,Dk_dry_up, &
                   Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                   alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                   res_pert_dn)
                                             
      J_pert_up(:,ideriv) = (res_pert_up(:)-res(:))/pert_up
      J_pert_dn(:,ideriv) = (res_pert_dn(:)-res(:))/pert_dn
    enddo
    
    Jup = J_pert_up
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_auxvar_pert_up)
    call GlobalAuxVarStrip(global_auxvar_pert_dn)    
    call MaterialAuxVarStrip(material_auxvar_pert_up)
    call MaterialAuxVarStrip(material_auxvar_pert_dn)    
  endif

end subroutine THFluxDerivative

! ************************************************************************** !
subroutine THFlux(auxvar_up,global_auxvar_up, &
                  material_auxvar_up, &
                  sir_up, &
                  Dk_up, &
                  auxvar_dn,global_auxvar_dn, &
                  material_auxvar_dn, &
                  sir_dn, &
                  Dk_dn, &
                  area, &
                  dist, &
                  upweight, &
                  option,v_darcy,Dk_dry_up, &
                  Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                  alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                  Res)
  ! 
  ! Computes the internal flux terms for the residual
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 
                  
  use Option_module                              
  use Connection_module
  use EOS_Water_module
  use Utility_module

  implicit none
  
  type(TH_auxvar_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: Dk_dry_up, Dk_dry_dn
  PetscReal :: Dk_ice_up, Dk_ice_dn
  PetscReal :: alpha_up, alpha_dn
  PetscReal :: alpha_fr_up, alpha_fr_dn
  PetscReal :: Dk_eff_up, Dk_eff_dn
  PetscReal :: v_darcy,area
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist(-1:3)
  PetscInt :: ispec
  PetscReal :: fluxm,fluxe,q
  PetscReal :: uh,ukvr,DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi
  PetscReal, parameter :: epsilon = 1.d-6

  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: perm_up, perm_dn

  ! ice variables
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: Ddiffgas_avg, Ddiffgas_up, Ddiffgas_dn
  PetscReal :: p_g
  PetscReal :: deng_up, deng_dn
  PetscReal :: psat_up, psat_dn
  PetscReal :: molg_up, molg_dn
  PetscReal :: satg_up, satg_dn
  PetscReal :: Diffg_up, Diffg_dn
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscErrorCode :: ierr
  PetscReal :: Ke_fr_up,Ke_fr_dn   ! frozen soil Kersten numbers
  PetscReal :: fv_up, fv_dn
  PetscReal :: ugas_ave, fdiffgas

  !
  call ConnectionCalculateDistances(dist,option%gravity,dd_up,dd_dn, &
                                    dist_gravity,upweight)
  call material_auxvar_up%PermeabilityTensorToScalar(dist,perm_up)
  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

  por_up = material_auxvar_up%porosity
  por_dn = material_auxvar_dn%porosity

  tor_up = material_auxvar_up%tortuosity
  tor_dn = material_auxvar_dn%tortuosity

  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)

  fluxm = 0.D0
  fluxe = 0.D0
  v_darcy = 0.D0  

  ! Flow term
  if (global_auxvar_up%sat(1) > sir_up .or. global_auxvar_dn%sat(1) > sir_dn) then
    if (global_auxvar_up%sat(1) < eps) then 
      upweight=0.d0
    else if (global_auxvar_dn%sat(1) < eps) then 
      upweight=1.d0
    endif
    density_ave = upweight*global_auxvar_up%den(1)+(1.D0-upweight)*global_auxvar_dn%den(1) 

    gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
              (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
              * dist_gravity

    if (option%ice_model /= DALL_AMICO &
        .and.option%ice_model /= PAINTER_KARRA_EXPLICIT_SMOOTH) then
      dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
    else
      dphi = auxvar_up%ice%pres_fh2o - auxvar_dn%ice%pres_fh2o + gravity
    endif

    if (dphi >= 0.D0) then
      ukvr = auxvar_up%kvr
      uh = auxvar_up%u
    else
      ukvr = auxvar_dn%kvr
      uh = auxvar_dn%u

    endif

    call InterfaceApprox(auxvar_up%kvr, auxvar_dn%kvr, dphi, &
                         option%rel_perm_aveg, ukvr)

    if (ukvr > floweps) then
      v_darcy = Dq * ukvr * dphi
   
      ! If only solving the energy equation, ensure Res(2) has no
      ! contribution from mass equation by setting darcy velocity
      ! to be zero
      if (option%flow%only_energy_eq) v_darcy = 0.d0

      q = v_darcy * area
        
      fluxm = fluxm + q*density_ave
      fluxe = fluxe + q*density_ave*uh
    endif
  endif 

  
  if (option%use_th_freezing) then
    ! Added by Satish Karra, 10/24/11
    satg_up = auxvar_up%ice%sat_gas
    satg_dn = auxvar_dn%ice%sat_gas
    if ((satg_up > eps) .and. (satg_dn > eps)) then

      p_g = option%reference_pressure ! set to reference pressure
      deng_up = auxvar_up%ice%den_gas
      deng_dn = auxvar_dn%ice%den_gas

      Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
      p_ref = 1.01325d5 ! in Pa
      T_ref = 25.d0 ! in deg C

      Diffg_up = Diffg_ref*(p_ref/p_g)*((max(-50d0, global_auxvar_up%temp) + 273.15d0) &
                /(T_ref + 273.15d0)) **(1.8d0)
      Diffg_dn = Diffg_ref*(p_ref/p_g)*((max(-50d0, global_auxvar_dn%temp) + 273.15d0) &
                /(T_ref + 273.15d0)) **(1.8d0)
           
      Ddiffgas_up = por_up*tor_up*satg_up*Diffg_up
      Ddiffgas_dn = por_dn*tor_dn*satg_dn*Diffg_dn

#if 0
      ! vapor pressure lowering due to capillary pressure
      fv_up = exp(-auxvar_up%pc/(global_auxvar_up%den(1)* &
           IDEAL_GAS_CONSTANT*(global_auxvar_up%temp + 273.15d0)))
      fv_dn = exp(-auxvar_dn%pc/(global_auxvar_dn%den(1)* &
           IDEAL_GAS_CONSTANT*(global_auxvar_dn%temp + 273.15d0)))
      molg_up = auxvar_up%ice%mol_gas*fv_up !psat_up*fv_up/p_g
      molg_dn = auxvar_dn%ice%mol_gas*fv_dn !psat_dn*fv_dn/p_g
#else
      ! seems having issues with vapor-only diffusion as above (not sure why)
      fv_up = 1.d0
      fv_dn = 1.d0
      molg_up = auxvar_up%ice%mol_gas
      molg_dn = auxvar_dn%ice%mol_gas
#endif

        
      if (deng_up*molg_up > deng_dn*molg_dn) then
        ! fmyuan: 'molg' is mole fraction of vapor in air, NOT vapor density
        upweight = 0.d0
        ugas_ave = auxvar_up%ice%u_gas
      else 
        upweight = 1.d0
        ugas_ave = auxvar_dn%ice%u_gas
      endif

      Ddiffgas_avg = upweight*Ddiffgas_up + (1.D0 - upweight)*Ddiffgas_dn 

      fdiffgas = area/(dd_up + dd_dn) * &
                 Ddiffgas_avg*(deng_up*molg_up - deng_dn*molg_dn)

      fluxm = fluxm + fdiffgas
      fluxe = fluxe + fdiffgas*ugas_ave

    endif

  endif ! if (use_th_freezing)

  if (option%use_th_freezing) then

    Dk_eff_up = auxvar_up%Dk_eff
    Dk_eff_dn = auxvar_dn%Dk_eff
  else

    Dk_eff_up = auxvar_up%Dk_eff
    Dk_eff_dn = auxvar_dn%Dk_eff

  endif
 
  Dk = (Dk_eff_up * Dk_eff_dn) / (dd_dn*Dk_eff_up + dd_up*Dk_eff_dn)
  cond = Dk*area*(global_auxvar_up%temp - global_auxvar_dn%temp)

  fluxe = fluxe + cond

  ! If only solving the energy equation, ensure Res(1) is zero
  if (option%flow%only_energy_eq) fluxm = 0.d0

  Res(1:option%nflowdof-1) = fluxm
  Res(option%nflowdof) = fluxe
  
 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL  
 
  
end subroutine THFlux

! ************************************************************************** !

subroutine THBCFluxDerivative(ibndtype,auxvars, &
                              auxvar_up,global_auxvar_up, &
                              auxvar_dn,global_auxvar_dn, &
                              material_auxvar_dn, &
                              sir_dn, &
                              Dk_dn, &
                              area, &
                              dist, &
                              option, &
                              sat_func_dn,&
                              Dk_dry_dn, &
                              Dk_ice_dn, &
                              Jdn)
  ! 
  ! Computes the derivatives of the boundary flux
  ! terms for the Jacobian
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 
  use Option_module
  use Characteristic_Curves_module
  use Connection_module
  use EOS_Water_module
  use Utility_module

  implicit none
  
  PetscInt :: ibndtype(:)
  type(TH_auxvar_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_dn
  PetscReal :: auxvars(:) ! from aux_real_var array in boundary condition
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn
  PetscReal :: area
  class(Characteristic_Curves_type), pointer :: sat_func_dn
  PetscReal :: Dk_dry_dn
  PetscReal :: Dk_ice_dn
  PetscReal :: alpha_dn
  PetscReal :: alpha_fr_dn
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscReal :: dist(-1:3)
  
  PetscReal :: dist_gravity  ! distance along gravity vector
          
  PetscReal :: dd_dn
  PetscInt :: ispec
  PetscReal :: v_darcy
  PetscReal :: fluxm,fluxe,q,density_ave
  PetscReal :: uh,ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi

  PetscReal :: ddiff_dp_dn, ddiff_dt_dn
  PetscReal :: dden_ave_dp_dn, dden_ave_dt_dn
  PetscReal :: dgravity_dden_dn
  PetscReal :: dphi_dp_dn, dphi_dt_dn
  PetscReal :: dukvr_dp_dn, dukvr_dt_dn
  PetscReal :: duh_dp_dn, duh_dt_dn
  PetscReal :: dq_dp_dn, dq_dt_dn
  PetscReal :: Dk_eff_dn
  PetscReal :: dDk_dt_dn, dDk_dp_dn
  PetscReal :: dKe_dt_dn, dKe_dp_dn
  PetscReal :: dKe_fr_dt_dn, dKe_fr_dp_dn

  PetscInt :: iphase, ideriv
  type(TH_auxvar_type) :: auxvar_pert_dn, auxvar_pert_up
  type(global_auxvar_type) :: global_auxvar_pert_dn, global_auxvar_pert_up
  class(material_auxvar_type), allocatable :: material_auxvar_pert_dn, &
                                              material_auxvar_pert_up

  PetscReal :: perturbation
  PetscReal :: x_up(option%nflowdof), x_dn(option%nflowdof)
  PetscReal :: x_pert_up(option%nflowdof), x_pert_dn(option%nflowdof)
  PetscReal :: pert_up, pert_dn
  PetscReal :: res(option%nflowdof)
  PetscReal :: res_pert_up(option%nflowdof)
  PetscReal :: res_pert_dn(option%nflowdof)
  PetscReal :: J_pert_dn(option%nflowdof,option%nflowdof)

  PetscBool :: hw_present

  ! ice variables
  PetscReal :: Ddiffgas_avg, Ddiffgas_up, Ddiffgas_dn
  PetscReal :: p_g
  PetscReal :: deng_up, deng_dn
  PetscReal :: psat_up, psat_dn
  PetscReal :: molg_up, molg_dn
  PetscReal :: satg_up, satg_dn
  PetscReal :: Diffg_up, Diffg_dn
  PetscReal :: ddeng_dt_dn, ddeng_dp_dn
  PetscReal :: dpsat_dt_dn
  PetscReal :: dmolg_dt_dn, dmolg_dp_dn
  PetscReal :: dDiffg_dt_dn
  PetscReal :: dDiffg_dp_dn
  PetscReal :: dsatg_dt_dn, dsatg_dp_dn
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscErrorCode :: ierr
  PetscReal :: v_darcy_allowable
  PetscReal :: temp_real
  PetscReal :: dum1
  PetscReal :: T_th,fct,fctT,dfctT_dT
  PetscReal :: rho
  PetscReal :: dq_lin,dP_lin
  PetscReal :: q_approx,dq_approx
  PetscBool :: skip_thermal_conduction
  PetscReal :: ugas_ave, dugas_ave_dt, dugas_ave_dp, fdiffgas, fdiffgas_dx
  PetscReal :: deltaTf

  PetscBool :: skip_mass_flow

!----------------------------------------------------------------------------------------
  skip_thermal_conduction = PETSC_FALSE
  skip_mass_flow = PETSC_FALSE
  T_th  = 0.5d0

  fluxm = 0.d0
  fluxe = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

  Jdn = 0.d0 
  
  dden_ave_dp_dn = 0.d0
  dden_ave_dt_dn = 0.d0
  ddiff_dp_dn = 0.d0
  ddiff_dt_dn = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_dn = 0.d0
  dphi_dt_dn = 0.d0
  dukvr_dp_dn = 0.d0
  dukvr_dt_dn = 0.d0
  duh_dp_dn = 0.d0
  duh_dt_dn = 0.d0
  dq_dp_dn = 0.d0
  dq_dt_dn = 0.d0

  hw_present = PETSC_FALSE
  if (associated(auxvar_dn%surface)) then
    hw_present = auxvar_dn%surface%surf_wat
  endif
        
  dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
  dd_dn = dist(0)

  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  por_dn = material_auxvar_dn%porosity
  tor_dn = material_auxvar_dn%tortuosity

  ! Flow
  diffdp = por_dn*tor_dn/dd_dn*area
  select case(ibndtype(TH_PRESSURE_DOF))
    ! figure out the direction of flow
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC, &
         HET_DIRICHLET_BC,HET_SEEPAGE_BC,HET_CONDUCTANCE_BC)
      if (ibndtype(TH_PRESSURE_DOF) == CONDUCTANCE_BC .or. &
          ibndtype(TH_PRESSURE_DOF) == HET_CONDUCTANCE_BC) then
        Dq = auxvars(TH_CONDUCTANCE_DOF)
      else
        Dq = perm_dn / dd_dn
      endif
      ! Flow term
      if (global_auxvar_up%sat(1) > sir_dn .or. &
          global_auxvar_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        
        density_ave = upweight*global_auxvar_up%den(1)+ &
                      (1.D0-upweight)*global_auxvar_dn%den(1)
        dden_ave_dp_dn = (1.D0-upweight)*auxvar_dn%dden_dp
        dden_ave_dt_dn = (1.D0-upweight)*auxvar_dn%dden_dt

        if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
          dden_ave_dt_dn = dden_ave_dt_dn + upweight*auxvar_up%dden_dt
        endif
        
        gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
                  (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
                  * dist_gravity
        dgravity_dden_dn = (1.d0-upweight)*auxvar_dn%avgmw*dist_gravity

        if (option%ice_model /= DALL_AMICO &
            .and.option%ice_model /= PAINTER_KARRA_EXPLICIT_SMOOTH) then
          dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
          dphi_dp_dn = -1.d0 + dgravity_dden_dn*auxvar_dn%dden_dp
          dphi_dt_dn = dgravity_dden_dn*auxvar_dn%dden_dt

        else
          if (ibndtype(TH_PRESSURE_DOF) == SEEPAGE_BC .or. &
              ibndtype(TH_PRESSURE_DOF) == HET_SEEPAGE_BC) then
            ! F.-M. Yuan: SEEPAGE_BC outlet 'pres(1)' is for water table,
            ! while %pres_fh2o may be incorrect due to its %temp usually as first soil layer or external input
            dphi = global_auxvar_up%pres(1) - auxvar_dn%ice%pres_fh2o + gravity
            dphi_dp_dn = -auxvar_dn%ice%dpres_fh2o_dp + &
                        dgravity_dden_dn*auxvar_dn%dden_dp
            dphi_dt_dn = -auxvar_dn%ice%dpres_fh2o_dt + &
                        dgravity_dden_dn*auxvar_dn%dden_dt

          else
            dphi = auxvar_up%ice%pres_fh2o - auxvar_dn%ice%pres_fh2o + gravity
            dphi_dp_dn = -auxvar_dn%ice%dpres_fh2o_dp + &
                        dgravity_dden_dn*auxvar_dn%dden_dp
            dphi_dt_dn = -auxvar_dn%ice%dpres_fh2o_dt + &
                        dgravity_dden_dn*auxvar_dn%dden_dt
          endif
        endif

        if (ibndtype(TH_PRESSURE_DOF) == SEEPAGE_BC .or. &
            ibndtype(TH_PRESSURE_DOF) == CONDUCTANCE_BC .or. &
            ibndtype(TH_PRESSURE_DOF) == HET_SEEPAGE_BC .or. &
            ibndtype(TH_PRESSURE_DOF) == HET_CONDUCTANCE_BC &
           ) then
          ! boundary cell is <= pref 
          ! commenting out the following 'if ... endif' block
          ! so that always shut-off 'flow-in' because boundary-cells CAN be under water-table
          !if (global_auxvar_up%pres(1)-option%reference_pressure < eps) then
            ! skip thermal conduction whenever water table is lower than cell
            skip_thermal_conduction = PETSC_TRUE
            ! flow inward
            if (dphi >= 0.d0) then
              dphi = 0.d0
              dphi_dp_dn = 0.d0
              dphi_dt_dn = 0.d0
              skip_mass_flow = PETSC_TRUE   ! also shut-off other (e.g. gas) mass flow, if any
            endif

          !endif
        endif

        if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
                                   !( dgravity_dden_up                   ) (dden_dt_up)
          dphi_dt_dn = dphi_dt_dn + upweight*auxvar_up%avgmw*dist_gravity*auxvar_up%dden_dt
        endif
        
        if (dphi>=0.D0) then
          ukvr = auxvar_up%kvr
          dukvr_dp_dn = auxvar_up%dkvr_dp
          dukvr_dt_dn = auxvar_up%dkvr_dt

#if 0
          ! the following may be incorrect - better to calc. 'auxvar' with appropriate cell/BC's P/T
          ! (NOT make sense that using UP's state but with DN's derivatives)
          if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dukvr_dt_dn = auxvar_dn%dkvr_dt
          endif
#endif

        else
          ukvr = auxvar_dn%kvr
          dukvr_dp_dn = auxvar_dn%dkvr_dp
          dukvr_dt_dn = auxvar_dn%dkvr_dt
        endif      

        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
          dq_dt_dn = Dq*(dukvr_dt_dn*dphi+ukvr*dphi_dt_dn)*area

          ! for SEEPAGE_BC: if conduction OFF, Q W.r.t T should be off as well.
          ! i.e. no permissivity dependence on T.
          if(skip_thermal_conduction) then
            dq_dt_dn = 0.d0
          endif

        endif

#ifdef CLM_PFLOTRAN
        ! when coupled with CLM, 'pressure-type' BC due to water-column formed on the BC-contacted cell(s)
        ! should not be producing 'darcy-flow' over standing water availability
        ! NOTE: this IS not appropriate for 'injection' caused pressure-type BC.
        temp_real = global_auxvar_up%pres(1) -  &
                    max(global_auxvar_dn%pres(1), option%reference_pressure)
        if (ibndtype(TH_PRESSURE_DOF) == DIRICHLET_BC .and. temp_real >= floweps) then
          call EOSWaterdensity(option%reference_temperature, &
                                option%reference_pressure,   &
                                rho,dum1,ierr)

          v_darcy_allowable = temp_real/option%flow_dt/(-option%gravity(3))/rho

          if (v_darcy > v_darcy_allowable) then
            v_darcy = v_darcy_allowable
            q = v_darcy * area
            ukvr = v_darcy/Dq/dphi
            dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
            dq_dt_dn = Dq*(dukvr_dt_dn*dphi+ukvr*dphi_dt_dn)*area
          endif

        endif     !if (ibndtype(TH_PRESSURE_DOF) == DIRICHLET_BC)
#endif


      endif

    case(HET_SURF_SEEPAGE_BC)
      Dq = perm_dn / dd_dn
      ! Flow term
      if (global_auxvar_up%sat(1) > sir_dn .or. global_auxvar_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        
        density_ave = upweight*global_auxvar_up%den(1)+(1.D0-upweight)*global_auxvar_dn%den(1)
        dden_ave_dp_dn = (1.D0-upweight)*auxvar_dn%dden_dp
        dden_ave_dt_dn = (1.D0-upweight)*auxvar_dn%dden_dt

        if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
          dden_ave_dt_dn = dden_ave_dt_dn + upweight*auxvar_up%dden_dt
        endif
        
        gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
                  (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
                  * dist_gravity
        dgravity_dden_dn = (1.d0-upweight)*auxvar_dn%avgmw*dist_gravity

        if (option%ice_model /= DALL_AMICO &
           .and.option%ice_model /= PAINTER_KARRA_EXPLICIT_SMOOTH) then
          dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
          dphi_dp_dn = -1.d0 + dgravity_dden_dn*auxvar_dn%dden_dp
          dphi_dt_dn = dgravity_dden_dn*auxvar_dn%dden_dt

        else
          dphi = auxvar_up%ice%pres_fh2o - auxvar_dn%ice%pres_fh2o + gravity
          dphi_dp_dn = -auxvar_dn%ice%dpres_fh2o_dp + dgravity_dden_dn*auxvar_dn%dden_dp
          dphi_dt_dn = -auxvar_dn%ice%dpres_fh2o_dt + dgravity_dden_dn*auxvar_dn%dden_dt
        endif

        ! flow in         ! boundary cell is <= pref
        if (dphi > 0.d0 .and. global_auxvar_up%pres(1)-option%reference_pressure < eps) then
          dphi = 0.d0
          dphi_dp_dn = 0.d0
          dphi_dt_dn = 0.d0
        endif


        if (option%surf_flow_on) then
          ! ---------------------------
          ! Surface-subsurface simulation
          ! ---------------------------
          
          ! If surface-water is frozen, zero out the darcy velocity
          if (global_auxvar_up%temp < 0.d0) then
            dphi = 0.d0
            dphi_dp_dn = 0.d0
            dphi_dt_dn = 0.d0
          endif
        endif

        if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
                                   !( dgravity_dden_up                   ) (dden_dt_up)
          dphi_dt_dn = dphi_dt_dn + upweight*auxvar_up%avgmw*dist_gravity*auxvar_up%dden_dt
        endif
        
        if (dphi>=0.D0) then
          ukvr = auxvar_up%kvr
          dukvr_dp_dn = auxvar_up%dkvr_dp
          dukvr_dt_dn = auxvar_up%dkvr_dt

#if 0
          ! the following may be incorrect - better to calc. 'auxvar' with appropriate cell/BC's P/T
          ! (NOT make sense that using UP's state but with DN's derivatives)
          if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dukvr_dt_dn = auxvar_dn%dkvr_dt
          endif
#endif

        else
          ukvr = auxvar_dn%kvr
          dukvr_dp_dn = auxvar_dn%dkvr_dp
          dukvr_dt_dn = auxvar_dn%dkvr_dt
        endif      

        ! (fmyuan: the following was commented out, but NOT in FLUX alone calculation - so better turn it on as well)
        call InterfaceApprox(auxvar_up%kvr, auxvar_dn%kvr, &
                             auxvar_up%dkvr_dp, auxvar_dn%dkvr_dp, &
                             dphi, &
                             option%rel_perm_aveg, &
                             ukvr, dum1, dukvr_dp_dn)

        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
          dq_dt_dn = Dq*(dukvr_dt_dn*dphi+ukvr*dphi_dt_dn)*area

          if (option%surf_flow_on .and. &
              option%subsurf_surf_coupling /= DECOUPLED) then

            ! ---------------------------
            ! Surface-subsurface simulation
            ! ---------------------------

            ! Temperature-smoothing
            fctT = 1.d0
            dfctT_dT = 0.d0
            if (global_auxvar_up%temp < 0.d0) then
              ! surface water is frozen, so no flow can occur
              fctT = 0.d0
              dfctT_dT = 0.d0
            else
              ! if subsurface is close to frozen, smoothly throttle down the flow
              if (global_auxvar_dn%temp < 0.d0) then
                fctT      = 0.d0
                dfctT_dT  = 0.d0
              else if (global_auxvar_dn%temp > T_th) then
                fctT      = 1.d0
                dfctT_dt  = 0.d0
              else
                fct      = 1.d0-(global_auxvar_dn%temp/T_th)**2.d0
                fctT     = 1.d0-fct**2.d0
                dfctT_dT = 4.d0*global_auxvar_dn%temp/(T_th*T_th)*fct
              endif
            endif

            ! Pressure-smoothing
            if (.not. auxvar_dn%surface%bcflux_default_scheme) then
              if (global_auxvar_dn%pres(1) <= auxvar_dn%surface%P_min) then

                ! Linear approximation
                call Interpolate(auxvar_dn%surface%range_for_linear_approx(2), &
                                 auxvar_dn%surface%range_for_linear_approx(1), &
                                 global_auxvar_dn%pres(1), &
                                 auxvar_dn%surface%range_for_linear_approx(4), &
                                 auxvar_dn%surface%range_for_linear_approx(3), &
                                 q_approx)
                v_darcy = q_approx/area

                dP_lin = auxvar_dn%surface%range_for_linear_approx(2) - &
                         auxvar_dn%surface%range_for_linear_approx(1)
                dq_lin = auxvar_dn%surface%range_for_linear_approx(4) - &
                         auxvar_dn%surface%range_for_linear_approx(3)
                if (abs(dP_lin) < 1.d-10) dP_lin = 1.d-10

                dq_dp_dn = dq_lin/dP_lin

                ! Approximation:
                ! q_approx = q_min + slope*(P_dn - P_min)
                !
                ! Derivative of approximation w.r.t T_dn:
                ! d(q_approx)/dT_dn = d(q_min)/dT_dn +
                !                     d(slope)/dT_dn*(P_dn - P_min) +
                !                     slope*(0 - d(P_min)/dT_dn)
                ! Note: 
                !   d(q_min)/dT_dn = 0
                !   d(P_min)/dT_dn = 0
                !                     slope*(0 - d(P_min)/dT_dn)
                dq_dt_dn = auxvar_dn%surface%dlinear_slope_dT* &
                           (global_auxvar_dn%pres(1) - &
                            auxvar_dn%surface%range_for_linear_approx(1))

              else if (global_auxvar_dn%pres(1) <= auxvar_dn%surface%P_max) then

                ! Cubic approximation
                call CubicPolynomialEvaluate(auxvar_dn%surface%coeff_for_cubic_approx, &
                                             global_auxvar_dn%pres(1) - option%reference_pressure, &
                                             q_approx, dq_approx)
                v_darcy = q_approx/area
                dq_dp_dn = dq_approx

                call CubicPolynomialEvaluate(auxvar_dn%surface%coeff_for_deriv_cubic_approx, &
                                             global_auxvar_dn%pres(1) - option%reference_pressure, &
                                             dq_dt_dn, dum1)

              endif
            endif

            ! Apply temperature smoothing
            v_darcy  = v_darcy*fctT
            q        = v_darcy*area
            dq_dp_dn = dq_dp_dn*fctT
            dq_dt_dn = dq_dt_dn*fctT + q*dfctT_dT

          endif
        endif
      endif
      
    case(NEUMANN_BC)
      if (dabs(auxvars(TH_PRESSURE_DOF)) > floweps) then
        v_darcy = auxvars(TH_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = global_auxvar_up%den(1)
          dden_ave_dp_dn = auxvar_up%dden_dp
          dden_ave_dt_dn = auxvar_up%dden_dt
#if 0
          ! the following may be incorrect - better to calc. 'auxvar' with appropriate cell/BC's P/T
          ! (NOT make sense that using UP's state but with DN's derivatives)
          if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dden_ave_dt_dn = auxvar_dn%dden_dt
          endif
#endif

        else 
          density_ave = global_auxvar_dn%den(1)
          dden_ave_dp_dn = auxvar_dn%dden_dp
          dden_ave_dt_dn = auxvar_dn%dden_dt
        endif 
        q = v_darcy * area

      endif

      ! if NO flow at all, turn-off any mass-flow (later on, e.g. for vapor-diffusion)
      if(abs(auxvars(TH_PRESSURE_DOF)) <= floweps) skip_mass_flow = PETSC_TRUE

    case(ZERO_GRADIENT_BC)
      ! do nothing, but by-passing default

      ! if flux-type BC for T, the fluid is totally energy form without mass (may not be needed, but just in case)
      if(ibndtype(TH_TEMPERATURE_DOF) == NEUMANN_BC) skip_mass_flow = PETSC_TRUE

  end select

  if (v_darcy > 0.D0) then
    uh = auxvar_up%u
    duh_dp_dn = auxvar_up%du_dp
    duh_dt_dn = auxvar_up%du_dt

#if 0
    ! the following may be incorrect - better to calc. 'auxvar' with appropriate cell/BC's P/T
    ! (NOT make sense that using UP's state but with DN's derivatives)
    if (ibndtype(TH_PRESSURE_DOF) == ZERO_GRADIENT_BC) then
      duh_dp_dn = auxvar_dn%du_dp
    endif
    if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
      duh_dt_dn = auxvar_dn%du_dt
    endif
#endif

  else
    uh = auxvar_dn%u
    duh_dp_dn = auxvar_dn%du_dp
    duh_dt_dn = auxvar_dn%du_dt
  endif

      
  ! If only solving the energy equation, ensure Jdn(2,2) has no
  ! contribution from mass equation
  if (option%flow%only_energy_eq) then
    q = 0.d0
    dq_dt_dn = 0.d0
  endif

  Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)
  Jdn(1,2) = (dq_dt_dn*density_ave+q*dden_ave_dt_dn)

  ! based on flux = q*density_ave*uh
  Jdn(option%nflowdof,1) =  &
     ((dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uh+q*density_ave*duh_dp_dn)
  Jdn(option%nflowdof,2) =  &
     ((dq_dt_dn*density_ave+q*dden_ave_dt_dn)*uh+q*density_ave*duh_dt_dn)

  ! Conduction term
  select case(ibndtype(TH_TEMPERATURE_DOF))
    case(DIRICHLET_BC,HET_DIRICHLET_BC)
      Dk =  auxvar_dn%Dk_eff / dd_dn
      !cond = Dk*area*(global_auxvar_up%temp-global_auxvar_dn%temp)

      if (skip_thermal_conduction) then
        ! skip thermal conducton when the boundary pressure is below the
        ! reference pressure (e.g. river stage is below cell center).

        Dk_eff_dn    = 0.d0
        dKe_dp_dn    = 0.d0
        dKe_dt_dn    = 0.d0
        dKe_fr_dt_dn = 0.d0
        dKe_fr_dp_dn = 0.d0

        dDk_dt_dn = 0.d0
        dDk_dp_dn = 0.d0
        Dk = 0.d0
      else
        if (option%use_th_freezing) then
          Dk_eff_dn    = auxvar_dn%Dk_eff
          dKe_dp_dn    = auxvar_dn%dKe_dp
          dKe_dt_dn    = auxvar_dn%dKe_dt
          dKe_fr_dt_dn = auxvar_dn%ice%dKe_fr_dt
          dKe_fr_dp_dn = auxvar_dn%ice%dKe_fr_dp
          Dk           = Dk_eff_dn/dd_dn

         ! Dk = Dk_eff_dn/dd_dn
         dDk_dp_dn = auxvar_dn%dDk_eff_dp/dd_dn  ! effective one
         dDk_dt_dn = auxvar_dn%dDk_eff_dt/dd_dn

        else
  
          Dk_eff_dn = auxvar_dn%Dk_eff
          dKe_dp_dn = auxvar_dn%dKe_dp
          dKe_dt_dn = auxvar_dn%dKe_dt
          Dk        = Dk_eff_dn/dd_dn
  
          dDk_dt_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn - Dk_dry_dn)*dKe_dt_dn
          dDk_dp_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn - Dk_dry_dn)*dKe_dp_dn

        endif
      endif

      if (.not. option%surf_flow_on) then
        ! ---------------------------
        ! Subsurface only simulation
        ! ---------------------------

        Jdn(option%nflowdof,1) = Jdn(option%nflowdof,1) + &
           area*(global_auxvar_up%temp - global_auxvar_dn%temp)*dDk_dp_dn

        Jdn(option%nflowdof,2) = Jdn(option%nflowdof,2) + Dk*area*(-1.d0) + &
           area*(global_auxvar_up%temp - global_auxvar_dn%temp)*dDk_dt_dn

      else
        ! ---------------------------
        ! Surface-subsurface simulation
        ! ---------------------------
        if (ibndtype(TH_PRESSURE_DOF) /= HET_SURF_SEEPAGE_BC) then
          if (.not.(hw_present)) then
            Jdn(option%nflowdof,1) = Jdn(option%nflowdof,1) + &
                area*(global_auxvar_up%temp - global_auxvar_dn%temp)*dDk_dp_dn

            Jdn(option%nflowdof,2) = Jdn(option%nflowdof,2) + Dk*area*(-1.d0) + &
                area*(global_auxvar_up%temp - global_auxvar_dn%temp)*dDk_dt_dn
          else
            Jdn = 0.d0
          endif
        else
          ! Only add contribution to Jacboian term for heat equation if
          ! standing water is present
          if (hw_present) then
            Jdn(option%nflowdof,1) = Jdn(option%nflowdof,1) + &
                area*(global_auxvar_up%temp - global_auxvar_dn%temp)*dDk_dp_dn

            Jdn(option%nflowdof,2) = Jdn(option%nflowdof,2) + Dk*area*(-1.d0) + &
                area*(global_auxvar_up%temp - global_auxvar_dn%temp)*dDk_dt_dn
          endif
        endif
      endif
      if (option%use_th_freezing .and. (.not.skip_mass_flow)) then
        ! Added by Satish Karra, 11/21/11
        satg_up = auxvar_up%ice%sat_gas
        satg_dn = auxvar_dn%ice%sat_gas
        if ((satg_up > eps) .and. (satg_dn > eps)) then
          p_g = option%reference_pressure  ! set to reference pressure
          deng_up = auxvar_up%ice%den_gas
          deng_dn = auxvar_dn%ice%den_gas
        
          Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
          p_ref = 1.01325d5   ! in Pa
          T_ref = 25.d0       ! in deg C

          Diffg_up = Diffg_ref*(p_ref/p_g)*((max(-50d0, global_auxvar_up%temp) + &
                 273.15d0)/(T_ref + 273.15d0)) **(1.8)
          Diffg_dn = Diffg_ref*(p_ref/p_g)*((max(-50d0, global_auxvar_dn%temp) + &
                 273.15d0)/(T_ref + 273.15d0)) **(1.8)

          dDiffg_dt_dn = 1.8d0*Diffg_dn/(global_auxvar_dn%temp + 273.15d0)
          if (global_auxvar_dn%temp<-50d0) dDiffg_dt_dn = 0.d0
          dDiffg_dp_dn = 0.d0

          Ddiffgas_up = satg_up*Diffg_up
          Ddiffgas_dn = satg_dn*Diffg_dn
          molg_up = auxvar_up%ice%mol_gas
          molg_dn = auxvar_dn%ice%mol_gas

          dmolg_dt_dn = auxvar_dn%ice%dmol_gas_dt  ! no need for '_up'
          dmolg_dp_dn = auxvar_dn%ice%dmol_gas_dp
          ddeng_dt_dn = auxvar_dn%ice%dden_gas_dt
          ddeng_dp_dn = auxvar_dn%ice%dden_gas_dp
          dsatg_dp_dn = auxvar_dn%ice%dsat_gas_dp
          dsatg_dt_dn = auxvar_dn%ice%dsat_gas_dt

          if (deng_up*molg_up > deng_dn*molg_dn) then
            ! 'molg' is the mole fraction of vapor in air (or partial-pressure), NOT density
            upweight = 0.d0
            ugas_ave = auxvar_up%ice%u_gas
            dugas_ave_dt = auxvar_up%ice%du_gas_dt
            dugas_ave_dp = auxvar_up%ice%du_gas_dp
          else
            upweight = 1.d0
            ugas_ave = auxvar_dn%ice%u_gas
            dugas_ave_dt = auxvar_dn%ice%du_gas_dt
            dugas_ave_dp = auxvar_dn%ice%du_gas_dp
          endif
        
          Ddiffgas_avg = upweight*Ddiffgas_up+(1.D0 - upweight)*Ddiffgas_dn
    
          ! derivaives: por_dn*tor_dn*Ddiffgas_avg*area/dd_dn*(deng_up*molg_up - deng_dn*molg_dn)
          fdiffgas = por_dn*tor_dn*area/dd_dn * &
                     Ddiffgas_avg*(deng_up*molg_up - deng_dn*molg_dn)

          fdiffgas_dx = por_dn*tor_dn*area/dd_dn * &
                 ( Ddiffgas_avg * (-ddeng_dp_dn*molg_dn - deng_dn*dmolg_dp_dn) + &
                   (deng_up*molg_up - deng_dn*molg_dn) * (1.d0-upweight) * por_dn*tor_dn* &
                    (satg_dn*dDiffg_dp_dn + dsatg_dp_dn*Diffg_dn) )
          Jdn(1,1) = Jdn(1,1) + fdiffgas_dx
          Jdn(2,1) = Jdn(2,1) + fdiffgas*dugas_ave_dp + ugas_ave * fdiffgas_dx  ! fdiffgas * ugas_ave


          fdiffgas_dx = por_dn*tor_dn*area/dd_dn * &
                 ( Ddiffgas_avg * (-ddeng_dt_dn*molg_dn - deng_dn*dmolg_dt_dn) + &
                   (deng_up*molg_up - deng_dn*molg_dn) * (1.d0-upweight) * por_dn*tor_dn* &
                    (satg_dn*dDiffg_dt_dn + dsatg_dt_dn*Diffg_dn) )
          Jdn(1,2) = Jdn(1,2) + fdiffgas_dx
          Jdn(2,2) = Jdn(2,2) + fdiffgas*dugas_ave_dt + ugas_ave * fdiffgas_dx  ! fdiffgas * ugas_ave

         endif
      endif ! if (use_th_freezing)

  end select

  ! If only solving the energy equation,
  !  - Set jacobian term corresponding to mass-equation to zero, and
  !  - Set off-diagonal jacobian terms to zero.
  if (option%flow%only_energy_eq) then
    Jdn(1,1) = 0.d0
    Jdn(1,2) = 0.d0
    Jdn(option%nflowdof,1) = 0.d0
  endif

#if 0
  if (option%flow%numerical_derivatives) then
    allocate(material_auxvar_pert_up,material_auxvar_pert_dn)

    call MaterialAuxVarInit(material_auxvar_pert_up,option)
    call MaterialAuxVarInit(material_auxvar_pert_dn,option)  
    
    call GlobalAuxVarInit(global_auxvar_pert_up,option)
    call GlobalAuxVarInit(global_auxvar_pert_dn,option)  
    call THAuxVarCopy(auxvar_up,auxvar_pert_up,option)
    call THAuxVarCopy(auxvar_dn,auxvar_pert_dn,option)
    call GlobalAuxVarCopy(global_auxvar_up,global_auxvar_pert_up,option)
    call GlobalAuxVarCopy(global_auxvar_dn,global_auxvar_pert_dn,option)
    
    call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_up, &
                            option)
    call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_dn, &
                            option)

    x_up(1) = global_auxvar_up%pres(1)
    x_up(2) = global_auxvar_up%temp
    x_dn(1) = global_auxvar_dn%pres(1)
    x_dn(2) = global_auxvar_dn%temp
    do ideriv = 1,3
      if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
        x_up(ideriv) = x_dn(ideriv)
      endif
    enddo
    if (option%use_th_freezing) then
       call THAuxVarComputeFreezing(x_dn,auxvar_dn, &
            global_auxvar_dn, &
            material_auxvar_dn, &
            iphase,sat_func_dn, &
            TH_parameter,ithrm_up, &
            option)
       call THAuxVarComputeFreezing(x_up,auxvar_up, &
            global_auxvar_up, &
            material_auxvar_up, &
            iphase,sat_func_dn, &
            TH_parameter,ithrm_up, &
            option)
    else
       call THAuxVarComputeNoFreezing(x_dn,auxvar_dn, &
            global_auxvar_dn, &
            material_auxvar_dn, &
            iphase,sat_func_dn, &
            option)
       call THAuxVarComputeNoFreezing(x_up,auxvar_up, &
            global_auxvar_up, &
            material_auxvar_up, &
            iphase,sat_func_dn, &
            option)
    endif
    
    call THBCFlux(ibndtype,auxvars,auxvar_up,global_auxvar_up, &
                  material_auxvar_up, &
                  auxvar_dn,global_auxvar_dn, &
                  material_auxvar_dn, &
                  sir_dn, &
                  Dk_dn, &
                  area,dist_gravity,option,v_darcy, &
                  fluxe_bulk, fluxe_cond, &
                  res)
    if (ibndtype(TH_PRESSURE_DOF) == ZERO_GRADIENT_BC .or. &
        ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC ) then
      x_pert_up = x_up
    endif

    do ideriv = 1,option%nflowdof
      pert_dn = x_dn(ideriv)*perturbation_tolerance    
      x_pert_dn = x_dn
     
      if (option%use_th_freezing) then
      
         if (ideriv == 1) then
            if (x_pert_dn(ideriv) < option%reference_pressure) then
               pert_dn = - pert_dn
            endif
            x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
         endif
      
         if (ideriv == 2) then
            if (x_pert_dn(ideriv) < 0.d0) then
               pert_dn = - 1.d-5
            else
               pert_dn = 1.d-5
            endif
            x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
         endif
      else
         x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      endif
        
      x_pert_up = x_up
      if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
        x_pert_up(ideriv) = x_pert_dn(ideriv)
      endif   

      if (option%use_th_freezing) then
         call THAuxVarComputeFreezing(x_pert_dn,auxvar_pert_dn, &
              global_auxvar_pert_dn, &
              material_auxvar_pert_dn, &
              iphase,sat_func_dn, &
              option)
         call THAuxVarComputeFreezing(x_pert_up,auxvar_pert_up, &
              global_auxvar_pert_up, &
              material_auxvar_pert_up, &
              iphase,sat_func_dn, &
              option)
      else
         call THAuxVarComputeNoFreezing(x_pert_dn,auxvar_pert_dn, &
              global_auxvar_pert_dn, &
              material_auxvar_pert_dn, &
              iphase,sat_func_dn, &
              option)
         call THAuxVarComputeNoFreezing(x_pert_up,auxvar_pert_up, &
              global_auxvar_pert_up, &
              material_auxvar_pert_up, &
              iphase,sat_func_dn, &
              option)
      endif

      call THBCFlux(ibndtype,auxvars,auxvar_pert_up,global_auxvar_pert_up, &
                    material_auxvar_pert_up, &
                    auxvar_pert_dn,global_auxvar_pert_dn, &
                    material_auxvar_pert_dn, &
                    sir_dn, &
                    Dk_dn, &
                    area,dist_gravity,option,v_darcy, &
                    fluxe_bulk, fluxe_cond, &
                    res_pert_dn)
      J_pert_dn(:,ideriv) = (res_pert_dn(:)-res(:))/pert_dn
    enddo
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_auxvar_pert_up)
    call GlobalAuxVarStrip(global_auxvar_pert_dn)      
  endif
#endif

end subroutine THBCFluxDerivative

! ************************************************************************** !

subroutine THBCFlux(ibndtype,auxvars,auxvar_up,global_auxvar_up, &
                    auxvar_dn,global_auxvar_dn, &
                    material_auxvar_dn, &
                    sir_dn, &
                    Dk_dn, &
                    area, &
                    dist, &
                    option,v_darcy, &
                    fluxe_bulk, fluxe_cond, &
                    Res)
  !
  ! Computes the  boundary flux terms for the residual
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 
  use Option_module
  use Connection_module
  use EOS_Water_module
  use Condition_module
  use Utility_module
 
  implicit none
  
  PetscInt :: ibndtype(:)
  type(TH_auxvar_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_dn
  PetscReal :: auxvars(:) ! from aux_real_var array
  PetscReal :: Dk_dn
  PetscReal :: v_darcy, area
  PetscReal :: Res(1:option%nflowdof) 
  PetscReal :: dist(-1:3)
  PetscReal, intent(out) :: fluxe_bulk, fluxe_cond
  
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dd_dn
          
  PetscReal :: por_dn,perm_dn,tor_dn
  PetscInt :: ispec
  PetscReal :: fluxm,fluxe,q,density_ave
  PetscReal :: uh,ukvr,diff,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi
  PetscReal :: dphi_orig
  PetscBool :: hw_present
  
  ! ice variables
  PetscReal :: Ddiffgas_avg, Ddiffgas_dn, Ddiffgas_up
  PetscReal :: p_g
  PetscReal :: deng_dn, deng_up
  PetscReal :: psat_dn, psat_up
  PetscReal :: molg_dn, molg_up
  PetscReal :: satg_dn, satg_up
  PetscReal :: Diffg_dn, Diffg_up
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscErrorCode :: ierr
  PetscReal :: fv_up, fv_dn
  PetscReal :: v_darcy_allowable, temp_real
  PetscReal :: rho,dum1
  PetscReal :: q_approx, dq_approx
  PetscReal :: T_th,fctT,fct
  PetscBool :: skip_thermal_conduction,  skip_mass_flow
  PetscReal :: fdiffgas, ugas_ave
  PetscReal :: deltaTf

  skip_thermal_conduction = PETSC_FALSE
  skip_mass_flow = PETSC_FALSE
  T_th  = 0.5d0

  fluxm = 0.d0
  fluxe = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0
  fctT = 0.d0
  fluxe_bulk = 0.d0
  fluxe_cond = 0.d0

  hw_present = PETSC_FALSE
  if (associated(auxvar_dn%surface)) then
    hw_present = auxvar_dn%surface%surf_wat
  endif

  dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
  dd_dn = dist(0)

  call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)
  por_dn = material_auxvar_dn%porosity
  tor_dn = material_auxvar_dn%tortuosity

  ! Flow
  diffdp = por_dn*tor_dn/dd_dn*area
  select case(ibndtype(TH_PRESSURE_DOF))
    case(DIRICHLET_BC,HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC, &
         HET_DIRICHLET_BC,HET_SEEPAGE_BC,HET_CONDUCTANCE_BC)
      if (ibndtype(TH_PRESSURE_DOF) == CONDUCTANCE_BC) then
        Dq = auxvars(TH_CONDUCTANCE_DOF)
      else
        Dq = perm_dn / dd_dn
      endif

      ! Flow term
      if (global_auxvar_up%sat(1) > sir_dn .or. &
          global_auxvar_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        density_ave = upweight*global_auxvar_up%den(1)+(1.D0-upweight)*global_auxvar_dn%den(1)
   
        gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
                  (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
                  * dist_gravity

        if (option%ice_model /= DALL_AMICO &
            .and.option%ice_model /= PAINTER_KARRA_EXPLICIT_SMOOTH) then
          dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity

        else
          if (ibndtype(TH_PRESSURE_DOF) == SEEPAGE_BC .or. &
              ibndtype(TH_PRESSURE_DOF) == HET_SEEPAGE_BC) then
          ! F.-M. Yuan: SEEPAGE_BC outlet 'pres(1)' is for water table,
          ! while %pres_fh2o may be incorrect due to its %temp usually as first soil layer or external input
            dphi = global_auxvar_up%pres(1) - auxvar_dn%ice%pres_fh2o + gravity

          else
            dphi = auxvar_up%ice%pres_fh2o - auxvar_dn%ice%pres_fh2o + gravity
          endif

        endif

        if (ibndtype(TH_PRESSURE_DOF) == SEEPAGE_BC .or. &
            ibndtype(TH_PRESSURE_DOF) == CONDUCTANCE_BC .or. &
            ibndtype(TH_PRESSURE_DOF) == HET_SEEPAGE_BC .or. &
            ibndtype(TH_PRESSURE_DOF) == HET_CONDUCTANCE_BC &
           ) then
          ! boundary cell is <= pref 

          ! commenting out the following 'if ... endif' block
          ! so that always shut-off 'flow-in' because boundary-cells CAN be under water-table
          !if (global_auxvar_up%pres(1)-option%reference_pressure < eps) then
            ! skip thermal conduction whenever water table is lower than cell
            skip_thermal_conduction = PETSC_TRUE
            ! flow inward
            if (dphi > 0.d0) then
              dphi = 0.d0
              skip_mass_flow = PETSC_TRUE    ! also shut-off other (e.g. gas) flow, if any
            endif
          !endif
        endif
        
        if (dphi>=0.D0) then
          ukvr = auxvar_up%kvr
        else
          ukvr = auxvar_dn%kvr
        endif      

        call InterfaceApprox(auxvar_up%kvr, auxvar_dn%kvr, &
                             dphi, &
                             option%rel_perm_aveg, &
                             ukvr)

        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
        endif


#ifdef CLM_PFLOTRAN
        ! when coupled with CLM, 'pressure-type' BC due to water-column formed on the BC-contacted cell(s)
        ! should not be producing 'darcy-flow' over standing water availability
        ! NOTE: this IS not appropriate for 'injection' caused pressure-type BC.
        temp_real = global_auxvar_up%pres(1) -  &
                    max(global_auxvar_dn%pres(1), option%reference_pressure)
        if (ibndtype(TH_PRESSURE_DOF) == DIRICHLET_BC .and. temp_real >= floweps) then
          call EOSWaterdensity(option%reference_temperature, &
                                option%reference_pressure,   &
                                rho,dum1,ierr)

          v_darcy_allowable = temp_real/option%flow_dt/(-option%gravity(3))/rho

          if (v_darcy > v_darcy_allowable) then
            v_darcy = v_darcy_allowable
            q = v_darcy * area
          endif

        endif     !if (ibndtype(TH_PRESSURE_DOF) == DIRICHLET_BC)
#endif


      endif 

    case(HET_SURF_SEEPAGE_BC)
      Dq = perm_dn / dd_dn
      ! Flow term
      if (global_auxvar_up%sat(1) > sir_dn .or. global_auxvar_dn%sat(1) > sir_dn) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then 
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then 
          upweight=1.d0
        endif
        density_ave = upweight*global_auxvar_up%den(1)+(1.D0-upweight)*global_auxvar_dn%den(1)
        
        gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
             (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
             * dist_gravity
        
        if (option%ice_model /= DALL_AMICO &
            .and.option%ice_model /= PAINTER_KARRA_EXPLICIT_SMOOTH) then
          dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
        else
          dphi = auxvar_up%ice%pres_fh2o - auxvar_dn%ice%pres_fh2o + gravity
        endif
        
        if (dphi > 0.d0 .and. global_auxvar_up%pres(1) - option%reference_pressure < eps) then
          dphi = 0.d0
        endif
        
        if (option%surf_flow_on) then
          ! If surface-water is frozen, zero out the darcy velocity
          if (global_auxvar_up%temp < 0.d0) then
            dphi = 0.d0
          endif
        endif
        
        if (dphi>=0.D0) then
          ukvr = auxvar_up%kvr
        else
          ukvr = auxvar_dn%kvr
        endif
        
        call InterfaceApprox(auxvar_up%kvr, auxvar_dn%kvr, &
             dphi, &
             option%rel_perm_aveg, &
             ukvr)
        
        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          
          if (option%surf_flow_on .and. &
               option%subsurf_surf_coupling /= DECOUPLED) then
            
            ! ---------------------------
            ! Surface-subsurface simulation
            ! ---------------------------
            
            ! Temperature-smoothing
            fctT = 1.d0
            if (global_auxvar_up%temp < 0.d0) then
              ! surface water is frozen, so no flow can occur
              fctT = 0.d0
            else 
              ! if subsurface is close to frozen, smoothly throttle down the flow
              if (global_auxvar_dn%temp < 0.d0) then
                fctT = 0.d0
              else if (global_auxvar_dn%temp > T_th) then
                fctT = 1.d0
              else
                fct  = 1.d0-(global_auxvar_dn%temp/T_th)**2.d0
                fctT = 1.d0-fct**2.d0
              endif
            endif

            ! If needed, apply pressure-smoothing
            if (.not. auxvar_dn%surface%bcflux_default_scheme) then
              if (global_auxvar_dn%pres(1) <= auxvar_dn%surface%P_min) then
                
                ! Linear approximation
                call Interpolate(auxvar_dn%surface%range_for_linear_approx(2), &
                     auxvar_dn%surface%range_for_linear_approx(1), &
                     global_auxvar_dn%pres(1), &
                     auxvar_dn%surface%range_for_linear_approx(4), &
                     auxvar_dn%surface%range_for_linear_approx(3), &
                     q_approx)
                v_darcy = q_approx/area
                
              else if (global_auxvar_dn%pres(1) <= auxvar_dn%surface%P_max) then
                
                ! Cubic approximation
                call CubicPolynomialEvaluate(auxvar_dn%surface%coeff_for_cubic_approx, &
                     global_auxvar_dn%pres(1)-option%reference_pressure, &
                     !global_auxvar_dn%pres(1), &
                     q_approx, dq_approx)
                v_darcy = q_approx/area
              endif
            endif
            
            ! Now apply temperature-smoothing
            v_darcy = v_darcy*fctT
            
          endif
          
        endif
      endif
      v_darcy = min(v_darcy,option%max_infiltration_velocity)

    case(NEUMANN_BC)
      if (dabs(auxvars(TH_PRESSURE_DOF)) > floweps) then
        v_darcy = auxvars(TH_PRESSURE_DOF)
        if (v_darcy > 0.d0) then 
          density_ave = global_auxvar_up%den(1)
        else 
          density_ave = global_auxvar_dn%den(1)
        endif 

      endif

      ! if NO flow at all, turn-off any mass-flow (later on, e.g. for vapor-diffusion)
      if(abs(auxvars(TH_PRESSURE_DOF)) <= floweps) skip_mass_flow = PETSC_TRUE

    case(ZERO_GRADIENT_BC)

      ! if flux-type BC for T, the fluid is totally energy form without mass
      ! (in this case, any 'aux_var%' calculation involving PRESSURE_DOF will virtually use boundary-cell's.)
      if(ibndtype(TH_TEMPERATURE_DOF) == NEUMANN_BC) skip_mass_flow = PETSC_TRUE

    case default
      option%io_buffer = 'BC type for H: "' // trim(GetSubConditionName(ibndtype(TH_PRESSURE_DOF))) // &
        '" not implemented in TH mode.'
      call printErrMsg(option)

  end select

  q = v_darcy * area
  ! If only solving the energy equation, ensure Res(2) has no
  ! contribution from mass equation by setting darcy velocity
  ! to be zero
  if (option%flow%only_energy_eq) q = 0.d0

  if (v_darcy > 0.D0) then
    uh = auxvar_up%u
  else
    uh = auxvar_dn%u
  endif


  fluxm = fluxm + q*density_ave
  fluxe = fluxe + q*density_ave*uh
  fluxe_bulk = q*density_ave*uh

  ! Conduction term
  select case(ibndtype(TH_TEMPERATURE_DOF))
    case(DIRICHLET_BC,HET_DIRICHLET_BC)
      Dk =  auxvar_dn%Dk_eff / dd_dn
      cond = Dk*area*(global_auxvar_up%temp-global_auxvar_dn%temp)

      if (skip_thermal_conduction) then
        ! skip thermal conducton when the boundary pressure is below the
        ! reference pressure (e.g. river stage is below cell center).
        cond = 0.d0
      endif

      if (option%surf_flow_on) then

        ! ---------------------------
        ! Surface-subsurface simulation
        ! ---------------------------

        ! Check if the pressure BC is associated with surface-flow model and
        ! there is no standing water, set heat conduction to be zero.
        if (ibndtype(TH_PRESSURE_DOF) == HET_SURF_SEEPAGE_BC .and. &
            .not.(hw_present)) then
          cond = 0.d0
        endif

        if (ibndtype(TH_PRESSURE_DOF) /= HET_SURF_SEEPAGE_BC .and. &
            (hw_present)) then
          cond = 0.d0
        endif

      else  ! subsurface only simulation

      endif

      fluxe = fluxe + cond
      fluxe_cond = cond

      if (option%use_th_freezing .and. (.not.skip_mass_flow)) then
        ! Added by Satish Karra,
        satg_up = auxvar_up%ice%sat_gas
        satg_dn = auxvar_dn%ice%sat_gas
        if ((satg_up > eps) .and. (satg_dn > eps)) then
          p_g = option%reference_pressure ! set to reference pressure
          deng_up = auxvar_up%ice%den_gas
          deng_dn = auxvar_dn%ice%den_gas
  
          Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
          p_ref = 1.01325d5   ! in Pa
          T_ref = 25.d0       ! in deg C

          Diffg_up = Diffg_ref*(p_ref/p_g)*((max(-50d0, global_auxvar_up%temp) + &
                 273.15d0)/(T_ref + 273.15d0)) **(1.8)
          Diffg_dn = Diffg_ref*(p_ref/p_g)*((max(-50d0, global_auxvar_dn%temp) + &
                 273.15d0)/(T_ref + 273.15d0)) **(1.8)

          Ddiffgas_up = satg_up*Diffg_up
          Ddiffgas_dn = satg_dn*Diffg_dn

          molg_up = auxvar_up%ice%mol_gas
          molg_dn = auxvar_dn%ice%mol_gas

          if (deng_up*molg_up > deng_dn*molg_dn) then   ! fmyuan: 'molg' is mole fraction of vapor in air, NOT vapor density
            upweight = 0.d0
            ugas_ave = auxvar_up%ice%u_gas
          else
            upweight = 1.d0
            ugas_ave = auxvar_dn%ice%u_gas
          endif

          if (ibndtype(TH_PRESSURE_DOF) == SEEPAGE_BC .or. &
              ibndtype(TH_PRESSURE_DOF) == HET_SEEPAGE_BC) then
            deltaTf = 1.0d-10              ! half-width of smoothing zone of Freezing-Thawing (by default, nearly NO smoothing)
            if(option%frzthw_halfwidth /= UNINITIALIZED_DOUBLE) deltaTf = max(deltaTf,option%frzthw_halfwidth)
            if(global_auxvar_dn%temp<=deltaTF) then   ! when iced-water exists, shut-off heat bulk outlet
              ugas_ave = 0.d0
            endif
          endif

          Ddiffgas_avg = upweight*Ddiffgas_up + (1.D0 - upweight)*Ddiffgas_dn

          fdiffgas = por_dn*tor_dn*Ddiffgas_avg*area/dd_dn * &
                  (deng_up*molg_up - deng_dn*molg_dn)
          fluxm = fluxm + fdiffgas
          fluxe = fluxe + fdiffgas*ugas_ave

         endif
      endif ! if (use_th_freezing)

    case(NEUMANN_BC)
      !geh: default internal energy units are MJ (option%scale = 1.d-6 is for J->MJ)
      fluxe = fluxe + auxvars(TH_TEMPERATURE_DOF)*area*(1.d6*option%scale) ! added by SK 10/18/11
      fluxe_cond = auxvars(TH_TEMPERATURE_DOF)*area*(1.d6*option%scale)

    case(ZERO_GRADIENT_BC)
      ! No change in fluxe
      fluxe_cond = 0.d0    ! need this for output

    case default
      option%io_buffer = 'BC type for T: "' // trim(GetSubConditionName(ibndtype(TH_TEMPERATURE_DOF))) // &
        '" not implemented in TH mode.'
      call printErrMsg(option)
  end select

  ! If only solving the energy equation, set Res(1) is 0.d0
  if (option%flow%only_energy_eq) fluxm = 0.d0

  Res(1:option%nflowspec) = fluxm
  Res(option%nflowdof) = fluxe

end subroutine THBCFlux

! ************************************************************************** !

subroutine THResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module
  use Variables_module
  use Material_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: cur_patch
  type(option_type), pointer :: option

  field => realization%field
  discretization => realization%discretization
  option => realization%option
  
 ! check initial guess -----------------------------------------------
  ierr = THInitGuessCheck(xx,option)
  if (ierr<0) then
    call SNESSetFunctionDomainError(snes,ierr);CHKERRQ(ierr)
    return
  endif

  ! Communication -----------------------------------------
  ! These 3 must be called before THUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc,field%iphas_loc,ONEDOF)
  call DiscretizationLocalToLocal(discretization,field%icap_loc,field%icap_loc,ONEDOF)

  call DiscretizationLocalToLocal(discretization,field%ithrm_loc,field%ithrm_loc,ONEDOF)
  
  call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_X,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_X,ZERO_INTEGER)

  call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Y,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Y,ZERO_INTEGER)

  call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Z,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Z,ZERO_INTEGER)

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THResidualPatch(snes,xx,r,realization,ierr)
    cur_patch => cur_patch%next
  enddo

end subroutine THResidual

! ************************************************************************** !

subroutine THResidualPatch(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation at patch level
  ! 
  ! Author: ???
  ! Date: 12/10/07
  !

  

  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module  
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
  
  implicit none

  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(inout) :: r
  type(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: i, jn
  PetscInt :: ip1, ip2
  PetscInt :: local_id, ghosted_id, local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: accum_p(:)

  PetscReal, pointer :: r_p(:), xx_loc_p(:), xx_p(:), yy_p(:)
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)

  PetscInt :: iphase
  PetscInt :: icap_up, icap_dn, ithrm_up, ithrm_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: dd, f_up, f_dn, ff
  PetscReal :: perm_up, perm_dn
  PetscReal :: D_up, D_dn  ! thermal conductivity wet constants at upstream, downstream faces.
  PetscReal :: Dk_dry_up, Dk_dry_dn ! dry thermal conductivities
  PetscReal :: Dk_ice_up, Dk_ice_dn ! frozen soil thermal conductivities
  PetscReal :: alpha_up, alpha_dn
  PetscReal :: alpha_fr_up, alpha_fr_dn
  PetscReal :: dw_kg, dw_mol
  PetscReal :: qsrc1, csrc1, enth_src_h2o, enth_src_co2 , esrc1

  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof)
  PetscReal :: Res_src(realization%option%nflowdof)
  PetscViewer :: viewer

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(TH_parameter_type), pointer :: TH_parameter
  type(TH_auxvar_type), pointer :: auxvars(:), auxvars_bc(:)
  type(TH_auxvar_type), pointer :: auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  !class(material_auxvar_type), pointer :: material_auxvars(:)  ! this unknownly would likely mess up %soil_properties(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set  
  type(sec_heat_type), pointer :: TH_sec_heat_vars(:)
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal, pointer :: mmsrc(:)
  PetscReal :: well_status
  PetscReal :: well_factor
  PetscReal :: pressure_bh
  PetscReal :: pressure_max
  PetscReal :: pressure_min
  PetscReal :: well_inj_water
  PetscReal :: Dq, dphi, v_darcy, ukvr

  PetscInt :: iconn, idof, istart, iend, ii
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  PetscReal :: vol_frac_prim
  PetscReal :: fluxe_bulk, fluxe_cond

  ! secondary continuum variables
  PetscReal :: sec_density
  PetscReal :: sec_dencpr
  PetscReal :: res_sec_heat

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  TH_parameter => patch%aux%TH%TH_parameter
  auxvars => patch%aux%TH%auxvars
  auxvars_bc => patch%aux%TH%auxvars_bc
  auxvars_ss => patch%aux%TH%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  TH_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars
  
  call THUpdateAuxVarsPatch(realization)
  ! override flags since they will soon be out of date  
  patch%aux%TH%auxvars_up_to_date = PETSC_FALSE

  if (option%compute_mass_balance_new) then
    call THZeroMassBalDeltaPatch(realization)
  endif


! now assign access pointer to local variables
  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90( r, r_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
 
  call VecGetArrayF90(field%flow_yy,yy_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr);CHKERRQ(ierr)
  !print *,' Finished scattering non deriv'
  
  if (option%surf_flow_on) call THComputeCoeffsForSurfFlux(realization)
  
  ! Calculating volume fractions for primary and secondary continua

  vol_frac_prim = 1.d0
  r_p = 0.d0

  ! Accumulation terms ------------------------------------
  r_p = - accum_p

  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1


    if (option%use_mc) then
      vol_frac_prim = TH_sec_heat_vars(local_id)%epsilon
    endif

    call THAccumulation(auxvars(ghosted_id),global_auxvars(ghosted_id), &
                        material_auxvars(ghosted_id), &
                        TH_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                        option,vol_frac_prim,Res)
    r_p(istart:iend) = r_p(istart:iend) + Res


    do ii=1,option%nflowdof
      if(Res(ii) /= Res(ii) .or. abs(Res(ii))>huge(Res(ii)) ) then
        write(string,*) 'local_id: ', local_id, 'Res: ', ii, Res(ii)
        option%io_buffer = ' NaN or INF of Residuals @ th.F90: THResidualPatch - Accumulation of ' // &
          trim(string)
        call printErrMsg(option)
      endif
    enddo

  enddo


  ! ================== Secondary continuum heat source terms =====================
  if (option%use_mc) then
  ! Secondary continuum contribution (Added by SK 06/02/2012)
  ! only one secondary continuum for now for each primary continuum node
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
    
      sec_dencpr = TH_parameter%dencpr(int(ithrm_loc_p(local_id))) ! secondary rho*c_p same as primary for now
        
      call THSecondaryHeat(TH_sec_heat_vars(local_id), &
                          global_auxvars(ghosted_id), &
!                         TH_parameter%ckdry(int(ithrm_loc_p(local_id))), &
                          TH_parameter%ckwet(int(ithrm_loc_p(local_id))), &
                          sec_dencpr, &
                          option,res_sec_heat)

      r_p(iend) = r_p(iend) - res_sec_heat*material_auxvars(ghosted_id)%volume
    enddo   
  endif
  ! ============== end secondary continuum heat source ===========================

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first 
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit

    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      iend = local_id * option%nflowdof
      istart = iend - option%nflowdof + 1
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
       
      if (source_sink%flow_condition%rate%itype /= HET_MASS_RATE_SS) &
        qsrc1 = source_sink%flow_condition%rate%dataset%rarray(1)
      
      Res_src = 0.d0
      select case (source_sink%flow_condition%rate%itype)
        case(MASS_RATE_SS)
          qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
        case(SCALED_MASS_RATE_SS)
          qsrc1 = qsrc1 / FMWH2O * & 
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn) ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc1 = qsrc1*global_auxvars(ghosted_id)%den(1) ! den = kmol/m^3 
        case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc1 = qsrc1*global_auxvars(ghosted_id)%den(1)* & ! den = kmol/m^3
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(HET_MASS_RATE_SS)
          qsrc1 = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)/FMWH2O

        case default
          write(string,*) source_sink%flow_condition%rate%itype
          option%io_buffer='TH mode source_sink%flow_condition%rate%itype = ' // &
          trim(adjustl(string)) // ', not implemented.'
      end select

      Res_src(TH_PRESSURE_DOF) = qsrc1

      esrc1 = 0.d0
      select case(source_sink%flow_condition%itype(TH_TEMPERATURE_DOF))
        case (ENERGY_RATE_SS)
          esrc1 = source_sink%flow_condition%energy_rate%dataset%rarray(1)
        case (SCALED_ENERGY_RATE_SS)
          esrc1 = source_sink%flow_condition%energy_rate%dataset%rarray(1) * &
                  source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case (HET_ENERGY_RATE_SS)
          esrc1 = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
      end select
      ! convert J/s --> MJ/s
      !geh: default internal energy units are MJ (option%scale = 1.d-6 is for 
      !     J->MJ)
      Res_src(TH_TEMPERATURE_DOF) = esrc1*1.d6*option%scale

      ! Update residual term associated with T
      if (qsrc1 > 0.d0) then ! injection
        Res_src(TH_TEMPERATURE_DOF) = Res_src(TH_TEMPERATURE_DOF) + &
          qsrc1*auxvars_ss(sum_connection)%u
      else
        ! extraction
        Res_src(TH_TEMPERATURE_DOF) = Res_src(TH_TEMPERATURE_DOF) + &
          qsrc1*auxvars(ghosted_id)%u
      endif

      r_p(istart:iend) = r_p(istart:iend) - Res_src

      if (option%compute_mass_balance_new) then
        global_auxvars_ss(sum_connection)%mass_balance_delta(1,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(1,1) - Res_src(1)
      endif
      if (associated(patch%ss_flow_vol_fluxes)) then
        ! fluid flux [m^3/sec] = qsrc_mol [kmol/sec] / den [kmol/m^3]
        patch%ss_flow_vol_fluxes(1,sum_connection) = qsrc1 / &
                                           global_auxvars(ghosted_id)%den(1)
      endif
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(1,sum_connection) = qsrc1
      endif

      Res = Res_src
      do ii=1,option%nflowdof
        if(Res(ii) /= Res(ii) .or. abs(Res(ii))>huge(Res(ii)) ) then
          write(string, *) ' name -', source_sink%name, ' @local_id -', local_id, 'with Res -', ii, Res(ii)
          option%io_buffer = ' NaN or INF of Residuals @ th.F90: THResidualPatch - source_sink of ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo

    enddo
    source_sink => source_sink%next
  enddo

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

      if (option%flow%only_vertical_flow) then
        !geh: place second conditional within first to avoid excessive
        !     dot products when .not. option%flow%only_vertical_flow
        if (dot_product(cur_connection_set%dist(1:3,iconn),unit_z) < &
            1.d-10) cycle
      endif

#if 0
      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)


! better to use the following to avoid inconsistence (fmyuan)
#else
      call ConnectionCalculateDistances(cur_connection_set%dist(:,iconn), &
                                    option%gravity,dd_up,dd_dn, &
                                    distance_gravity,upweight)
#endif
        
      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))
   
      D_up = TH_parameter%ckwet(ithrm_up)
      D_dn = TH_parameter%ckwet(ithrm_dn)
      
      Dk_dry_up = TH_parameter%ckdry(ithrm_up)
      Dk_dry_dn = TH_parameter%ckdry(ithrm_dn)
      
      alpha_up = TH_parameter%alpha(ithrm_up)
      alpha_dn = TH_parameter%alpha(ithrm_dn)

      if (option%use_th_freezing) then
         Dk_ice_up = TH_parameter%ckfrozen(ithrm_up)
         DK_ice_dn = TH_parameter%ckfrozen(ithrm_dn)
      
         alpha_fr_up = TH_parameter%alpha_fr(ithrm_up)
         alpha_fr_dn = TH_parameter%alpha_fr(ithrm_dn)
      else
         Dk_ice_up = Dk_dry_up
         Dk_ice_dn = Dk_dry_dn
      
         alpha_fr_up = alpha_up
         alpha_fr_dn = alpha_dn
      endif

      call THFlux(auxvars(ghosted_id_up),global_auxvars(ghosted_id_up), &
                  material_auxvars(ghosted_id_up), &
                  TH_parameter%sir(1,icap_up), &
                  D_up, &
                  auxvars(ghosted_id_dn),global_auxvars(ghosted_id_dn), &
                  material_auxvars(ghosted_id_dn), &
                  TH_parameter%sir(1,icap_dn), &
                  D_dn, &
                  cur_connection_set%area(iconn), &
                  cur_connection_set%dist(:,iconn), &
                  upweight,option,v_darcy,Dk_dry_up, &
                  Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                  alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                  Res)

      patch%internal_velocities(1,sum_connection) = v_darcy
      patch%internal_flow_fluxes(:,sum_connection) = Res(:)

      if (local_id_up>0) then
        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
      endif
   
      if (local_id_dn>0) then
        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:option%nflowdof)
      endif

      do ii=1,option%nflowdof
        if(Res(ii) /= Res(ii) .or. abs(Res(ii))>huge(Res(ii)) ) then
          write(string, *) ' @local_id -', local_id_up, local_id_dn, ' with Res -', ii, Res(ii)
          option%io_buffer = ' NaN or INF of Residuals @ th.F90: THResidualPatch - interior flux between ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo

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

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ithrm_dn = int(ithrm_loc_p(ghosted_id))
      D_dn = TH_parameter%ckwet(ithrm_dn)

      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))

      icap_dn = int(icap_loc_p(ghosted_id))
  
      call THBCFlux(boundary_condition%flow_condition%itype, &
                                boundary_condition%flow_aux_real_var(:,iconn), &
                                auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                auxvars(ghosted_id), &
                                global_auxvars(ghosted_id), &
                                material_auxvars(ghosted_id), &
                                TH_parameter%sir(1,icap_dn), &
                                D_dn, &
                                cur_connection_set%area(iconn), &
                                cur_connection_set%dist(-1:3,iconn), &
                                option, &
                                v_darcy, &
                                fluxe_bulk, fluxe_cond, &
                                Res)

      patch%boundary_velocities(1,sum_connection) = v_darcy
      patch%boundary_flow_fluxes(:,sum_connection) = Res(:)
      patch%boundary_energy_flux(1,sum_connection) = fluxe_bulk
      patch%boundary_energy_flux(2,sum_connection) = fluxe_cond

      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_bc(sum_connection)%mass_balance_delta(1,1) = &
          global_auxvars_bc(sum_connection)%mass_balance_delta(1,1) - Res(1)
      endif

      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      r_p(istart:iend)= r_p(istart:iend) - Res(1:option%nflowdof)

      do ii=1,option%nflowdof
        if(Res(ii) /= Res(ii) .or. abs(Res(ii))>huge(Res(ii)) ) then
          write(string, *) ' name -', boundary_condition%name, ' @local_id -', local_id, 'with Res -', ii, Res(ii)
          option%io_buffer = ' NaN or INF of Residuals @ th.F90: THResidualPatch - boundary_condition of ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! scale the residual by the volume
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    r_p (istart:iend)= r_p(istart:iend)/material_auxvars(ghosted_id)%volume
  enddo

  if (option%use_isothermal) then
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      istart = TWO_INTEGER + (local_id-1)*option%nflowdof
      r_p(istart)=xx_loc_p(2 + (ghosted_id-1)*option%nflowdof)-yy_p(istart-1)
    enddo
  endif

  if (patch%aux%TH%inactive_cells_exist) then
    do i=1,patch%aux%TH%n_zero_rows
      r_p(patch%aux%TH%zero_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_yy, yy_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr);CHKERRQ(ierr)

  if (realization%debug%vecview_residual) then
    string = 'THresidual'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%vecview_solution) then
    string = 'THxx'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

end subroutine THResidualPatch

! ************************************************************************** !

subroutine THJacobian(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
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
  type(patch_type), pointer :: cur_patch
  type(grid_type),  pointer :: grid
  type(option_type),  pointer :: option
  PetscReal :: norm
  
  character(len=MAXSTRINGLENGTH) :: string

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THJacobianPatch(snes,xx,J,J,realization,ierr)
    cur_patch => cur_patch%next
  enddo
  
  if (realization%debug%matview_Jacobian) then
    string = 'THjacobian'
    call DebugCreateViewer(realization%debug,string,realization%option,viewer)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    option => realization%option
    call MatNorm(J,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option)
    call MatNorm(J,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option)
    call MatNorm(J,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option)
  endif
  
end subroutine THJacobian

! ************************************************************************** !

subroutine THJacobianPatch(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: ???
  ! Date: 12/13/07
  ! 
       
  

  use Connection_module
  use Option_module
  use Grid_module
  use Realization_Subsurface_class
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: nvar,neq,nr
  PetscInt :: ithrm_up, ithrm_dn, i
  PetscInt :: ip1, ip2 

  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: iphase_loc_p(:), icap_loc_p(:), ithrm_loc_p(:)
  PetscInt :: icap,iphas,icap_up,icap_dn
  PetscInt :: ii, jj
  PetscReal :: dw_kg,dw_mol,enth_src_co2,enth_src_h2o,rho
  PetscReal :: qsrc1,csrc1
  PetscReal :: dd_up, dd_dn, dd, f_up, f_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: D_up, D_dn  ! thermal conductivity wet constants at upstream, downstream faces.
  PetscReal :: Dk_dry_up, Dk_dry_dn ! dry thermal conductivities
  PetscReal :: Dk_ice_up, Dk_ice_dn ! frozen soil thermal conductivities
  PetscReal :: alpha_up, alpha_dn
  PetscReal :: alpha_fr_up, alpha_fr_dn
  PetscReal :: zero, norm
  PetscReal :: upweight
  PetscReal :: max_dev  
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
            Jdn(realization%option%nflowdof,realization%option%nflowdof), &
            Jsrc(realization%option%nflowdof,realization%option%nflowdof)
  
  PetscInt :: istart, iend
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, idof
  PetscInt :: sum_connection  
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity 
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(TH_parameter_type), pointer :: TH_parameter
  type(TH_auxvar_type), pointer :: auxvars(:), auxvars_bc(:),auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 
  !class(material_auxvar_type), pointer :: material_auxvars(:)  ! this unknownly would likely mess up %soil_properties(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  type(sec_heat_type), pointer :: sec_heat_vars(:)
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: ithrm

  PetscViewer :: viewer
  Vec :: debug_vec
  PetscReal :: vol_frac_prim
  
  ! secondary continuum variables
  PetscReal :: area_prim_sec
  PetscReal :: jac_sec_heat

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  TH_parameter => patch%aux%TH%TH_parameter
  auxvars => patch%aux%TH%auxvars
  auxvars_bc => patch%aux%TH%auxvars_bc
  auxvars_ss => patch%aux%TH%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars

  sec_heat_vars => patch%aux%SC_heat%sec_heat_vars
  
#if 0
   call THNumericalJacobianTest(xx,realization)
#endif

  call VecGetArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%iphas_loc, iphase_loc_p, ierr);CHKERRQ(ierr)
  
  vol_frac_prim = 1.d0

  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    ! Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    icap = int(icap_loc_p(ghosted_id))
    
    if (option%use_mc) then    
      vol_frac_prim = sec_heat_vars(local_id)%epsilon
    endif

    ithrm = int(ithrm_loc_p(ghosted_id))
    call THAccumDerivative(auxvars(ghosted_id),global_auxvars(ghosted_id), &
                            material_auxvars(ghosted_id), &
                            TH_parameter%dencpr(ithrm), &
                            TH_parameter, ithrm, option, &
                            patch%characteristic_curves_array(icap)%ptr, &
                            vol_frac_prim,Jup) 

    if (option%use_mc) then
      call THSecondaryHeatJacobian(sec_heat_vars(local_id), &
                        TH_parameter%ckwet(int(ithrm_loc_p(local_id))), &
                        TH_parameter%dencpr(int(ithrm_loc_p(local_id))), &
                        option,jac_sec_heat)
                        
      Jup(option%nflowdof,2) = Jup(option%nflowdof,2) - &
                               jac_sec_heat*material_auxvars(ghosted_id)%volume
    endif

    do ii=1,option%nflowdof
      do jj=1,option%nflowdof
        if(Jup(ii,jj) /= Jup(ii,jj) &
           .or. abs(Jup(ii,jj))>huge(Jup(ii,jj)) ) then
          write(string, *) ' @local_id -', local_id, 'with Jacobin -', ii,jj,Jup(ii,jj)
          option%io_buffer = ' NaN or INF of Jacobians @ th.F90: THJacobinPatch - Accumulation ' // &
            trim(string)
          call printErrMsg(option)
        endif
      enddo
    enddo

                            
    ! scale by the volume of the cell
    Jup = Jup/material_auxvars(ghosted_id)%volume

    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo


  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_accum'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

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

      if (patch%imat(ghosted_id) <= 0) cycle

      if (source_sink%flow_condition%rate%itype /= HET_MASS_RATE_SS) &
        qsrc1 = source_sink%flow_condition%rate%dataset%rarray(1)
      
      select case (source_sink%flow_condition%rate%itype)
        case(MASS_RATE_SS)
          qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
        case(SCALED_MASS_RATE_SS)
          qsrc1 = qsrc1 / FMWH2O * & 
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn) ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc1 = qsrc1*global_auxvars(ghosted_id)%den(1) ! den = kmol/m^3 
        case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc1 = qsrc1*global_auxvars(ghosted_id)%den(1)* & ! den = kmol/m^3
                   source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(HET_MASS_RATE_SS)
          qsrc1 = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)/FMWH2O
        case default
          write(string,*) source_sink%flow_condition%rate%itype
          option%io_buffer='TH mode source_sink%flow_condition%rate%itype = ' // &
          trim(adjustl(string)) // ', not implemented.'
      end select

      Jsrc = 0.d0

      if (qsrc1 > 0.d0) then ! injection
        Jsrc(TH_TEMPERATURE_DOF,TH_PRESSURE_DOF) = & 
          -qsrc1*auxvars_ss(sum_connection)%du_dp
        ! dresT_dt = -qsrc1*hw_dt ! since tsrc1 is prescribed, there is no derivative
        istart = ghosted_id*option%nflowdof
      else
        ! extraction
        Jsrc(TH_TEMPERATURE_DOF,TH_PRESSURE_DOF) = &
          -qsrc1*auxvars(ghosted_id)%du_dp
        Jsrc(TH_TEMPERATURE_DOF,TH_TEMPERATURE_DOF) = &
          -qsrc1*auxvars(ghosted_id)%du_dt
        istart = ghosted_id*option%nflowdof
      endif
      
      do ii=1,option%nflowdof
        do jj=1,option%nflowdof
          if(Jsrc(ii,jj) /= Jsrc(ii,jj) &
             .or. abs(Jsrc(ii,jj))>huge(Jsrc(ii,jj)) ) then
            write(string, *) ' name -', source_sink%name, ' @local_id -', local_id, 'with Jacobin -', ii,jj, Jsrc
            option%io_buffer = ' NaN or INF of Jacobians @ th.F90: THJacobinPatch - Source_Sink of ' // &
              trim(string)
            call printErrMsg(option)
          endif
        enddo
      enddo

      ! scale by the volume of the cell
      Jsrc = Jsrc/material_auxvars(ghosted_id)%volume
         
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jsrc, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)

    
    enddo
    source_sink => source_sink%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_srcsink'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
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

      if (patch%imat(ghosted_id_up) <= 0 .or. &
          patch%imat(ghosted_id_dn) <= 0) cycle

      if (option%flow%only_vertical_flow) then
        !geh: place second conditional within first to avoid excessive
        !     dot products when .not. option%flow%only_vertical_flow
        if (dot_product(cur_connection_set%dist(1:3,iconn),unit_z) < &
            1.d-10) cycle
      endif

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

#if 0
      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)


! better to use the following to avoid inconsistence
#else
      call ConnectionCalculateDistances(cur_connection_set%dist(:,iconn), &
                                    option%gravity,dd_up,dd_dn, &
                                    distance_gravity,upweight)
#endif
    
      ithrm_up = int(ithrm_loc_p(ghosted_id_up))
      ithrm_dn = int(ithrm_loc_p(ghosted_id_dn))
      
      D_up = TH_parameter%ckwet(ithrm_up)
      D_dn = TH_parameter%ckwet(ithrm_dn)
    
      Dk_dry_up = TH_parameter%ckdry(ithrm_up)
      Dk_dry_dn = TH_parameter%ckdry(ithrm_dn)
      
      alpha_up = TH_parameter%alpha(ithrm_up)
      alpha_dn = TH_parameter%alpha(ithrm_dn)

      if (option%use_th_freezing) then
         Dk_ice_up = TH_parameter%ckfrozen(ithrm_up)
         DK_ice_dn = TH_parameter%ckfrozen(ithrm_dn)
      
         alpha_fr_up = TH_parameter%alpha_fr(ithrm_up)
         alpha_fr_dn = TH_parameter%alpha_fr(ithrm_dn)
      else
         Dk_ice_up = Dk_dry_up
         Dk_ice_dn = Dk_dry_dn
      
         alpha_fr_up = alpha_up
         alpha_fr_dn = alpha_dn
      endif

      icap_up = int(icap_loc_p(ghosted_id_up))
      icap_dn = int(icap_loc_p(ghosted_id_dn))

      call THFluxDerivative(auxvars(ghosted_id_up),global_auxvars(ghosted_id_up), &
                             material_auxvars(ghosted_id_up), &
                             TH_parameter%sir(1,icap_up), &
                             D_up, &
                             ithrm_up, &
                             auxvars(ghosted_id_dn),global_auxvars(ghosted_id_dn), &
                             material_auxvars(ghosted_id_dn), &
                             TH_parameter%sir(1,icap_dn), &
                             D_dn, &
                             ithrm_dn, &
                             cur_connection_set%area(iconn), &
                             cur_connection_set%dist(-1:3,iconn), &
                             upweight,option, &
                             patch%characteristic_curves_array(icap_up)%ptr, &
                             patch%characteristic_curves_array(icap_dn)%ptr, &
                             Dk_dry_up,Dk_dry_dn, &
                             Dk_ice_up,Dk_ice_dn, &
                             alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                             TH_parameter, &
                             Jup,Jdn)
      
      do ii=1,option%nflowdof
        do jj=1,option%nflowdof
          if(Jup(ii,jj) /= Jup(ii,jj) .or. Jdn(ii,jj) /= Jdn(ii,jj) &
             .or. abs(Jup(ii,jj))>huge(Jup(ii,jj)) .or. abs(Jdn(ii,jj))>huge(Jdn(ii,jj)) ) then
            write(string, *) ' between local_id up/dn -', local_id_up, local_id_dn, &
                'with Jacobin -', ii,jj, Jup(ii,jj), Jdn(ii,jj)
            option%io_buffer = ' NaN or INF of Jacobians @ th.F90: THJacobinPatch - Interior flux ' // &
                trim(string)
            call printErrMsg(option)
          endif
        enddo
      enddo


!  scale by the volume of the cell                      
      
      if (local_id_up > 0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup/material_auxvars(ghosted_id_up)%volume,ADD_VALUES, &
                                      ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn/material_auxvars(ghosted_id_up)%volume,ADD_VALUES, &
                                      ierr);CHKERRQ(ierr)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      Jdn/material_auxvars(ghosted_id_dn)%volume,ADD_VALUES, &
                                      ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      Jup/material_auxvars(ghosted_id_dn)%volume,ADD_VALUES, &
                                      ierr);CHKERRQ(ierr)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_flux'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

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

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      ithrm_dn  = int(ithrm_loc_p(ghosted_id))
      D_dn      = TH_parameter%ckwet(ithrm_dn)
      Dk_dry_dn = TH_parameter%ckdry(ithrm_dn)
      alpha_dn  = TH_parameter%alpha(ithrm_dn)

      icap_dn = int(icap_loc_p(ghosted_id))

      if (option%use_th_freezing) then
         DK_ice_dn = TH_parameter%ckfrozen(ithrm_dn)
         alpha_fr_dn = TH_parameter%alpha_fr(ithrm_dn)
      else
         Dk_ice_dn = Dk_dry_dn
         alpha_fr_dn = alpha_dn
      endif

      call THBCFluxDerivative(boundary_condition%flow_condition%itype, &
                              boundary_condition%flow_aux_real_var(:,iconn), &
                              auxvars_bc(sum_connection), &
                              global_auxvars_bc(sum_connection), &
                              auxvars(ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              TH_parameter%sir(1,icap_dn), &
                              D_dn, &
                              cur_connection_set%area(iconn), &
                              cur_connection_set%dist(-1:3,iconn), &
                              option, &
                              patch%characteristic_curves_array(icap_dn)%ptr,&
                              Dk_dry_dn,Dk_ice_dn, &
                              Jdn)
      Jdn = -Jdn

      do ii=1,option%nflowdof
        do jj=1,option%nflowdof
          if(Jdn(ii,jj) /= Jdn(ii,jj) .or. &
             abs(Jdn(ii,jj))>huge(Jdn(ii,jj)) ) then
            write(string, *) ' name -', boundary_condition%name, ' @local_id -', local_id, 'with Jacobin -', ii,jj, Jdn(ii,jj)
            option%io_buffer = ' NaN or INF of Jacobians @ th.F90: THJacobinPatch - Boundary_Condition of ' // &
                trim(string)
            call printErrMsg(option)
          endif
        enddo
      enddo

      !  scale by the volume of the cell
      Jdn = Jdn/material_auxvars(ghosted_id)%volume
      
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn,ADD_VALUES, &
                                    ierr);CHKERRQ(ierr)
 
    enddo
    boundary_condition => boundary_condition%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_bcflux'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  
  call VecRestoreArrayF90(field%flow_xx_loc, xx_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%ithrm_loc, ithrm_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%icap_loc, icap_loc_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%iphas_loc, iphase_loc_p, ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

! zero out isothermal and inactive cells
#ifdef ISOTHERMAL_MODE_DOES_NOT_WORK
  zero = 0.d0
  call MatZeroRowsLocal(A,n_zero_rows,zero_rows_local_ghosted,zero, &
                        PETSC_NULL_VEC,PETSC_NULL_VEC, &
                        ierr);CHKERRQ(ierr)
  do i=1, n_zero_rows
    ii = mod(zero_rows_local(i),option%nflowdof)
    ip1 = zero_rows_local_ghosted(i)
    if (ii == 0) then
      ip2 = ip1-1
    else if (ii == option%nflowdof-1) then
      ip2 = ip1+1
    else
      ip2 = ip1
    endif
    call MatSetValuesLocal(A,1,ip1,1,ip2,1.d0,INSERT_VALUES, &
                           ierr);CHKERRQ(ierr)
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
#else
  if (patch%aux%TH%inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(A,patch%aux%TH%n_zero_rows, &
                          patch%aux%TH%zero_rows_local_ghosted,f_up, &
                          PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
  endif
#endif

end subroutine THJacobianPatch

! ************************************************************************** !

subroutine THMaxChange(realization,dpmax,dtmpmax)
  ! 
  ! Computes the maximum change in the solution vector
  ! 
  ! Author: ???
  ! Date: 01/15/08
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Field_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  
  PetscReal :: dpmax, dtmpmax
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field

  dpmax = 0.d0
  dtmpmax = 0.d0
  
  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy, &
                ierr);CHKERRQ(ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,dpmax, &
                     ierr);CHKERRQ(ierr)
  call VecStrideNorm(field%flow_dxx,ONE_INTEGER,NORM_INFINITY,dtmpmax, &
                     ierr);CHKERRQ(ierr)
    
end subroutine THMaxChange

! ************************************************************************** !

subroutine THResidualToMass(realization)
  ! 
  ! Computes mass balance from residual equation
  ! 
  ! Author: ???
  ! Date: 12/10/07
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module
  use Grid_module

  implicit none

  Vec :: ts_mass_balance
  type(realization_subsurface_type) :: realization
  
  type(field_type), pointer :: field
  type(patch_type), pointer :: cur_patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  PetscReal, pointer :: mass_balance_p(:)
  type(TH_auxvar_type), pointer :: auxvars(:) 
  type(global_auxvar_type), pointer :: global_auxvars(:) 
  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart
  
  option => realization%option
  field => realization%field

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    grid => cur_patch%grid
    auxvars => cur_patch%aux%TH%auxvars

    call VecGetArrayF90(field%flow_ts_mass_balance,mass_balance_p,  &
                        ierr);CHKERRQ(ierr)
  
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (cur_patch%imat(ghosted_id) <= 0) cycle
        
      istart = (ghosted_id-1)*option%nflowdof+1
      mass_balance_p(istart) = mass_balance_p(istart)/ &
                                global_auxvars(ghosted_id)%den(1)* &
                                global_auxvars(ghosted_id)%den_kg(1)
    enddo

    call VecRestoreArrayF90(field%flow_ts_mass_balance,mass_balance_p,  &
                            ierr);CHKERRQ(ierr)

    cur_patch => cur_patch%next
  enddo

end subroutine THResidualToMass

! ************************************************************************** !

function THGetTecplotHeader(realization,icolumn)
  ! 
  ! THLiteGetTecplotHeader: Returns TH contribution to
  ! Tecplot file header
  ! 
  ! Author: ???
  ! Date: 02/13/08
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Field_module

  implicit none
  
  character(len=MAXSTRINGLENGTH) :: THGetTecplotHeader
  type(realization_subsurface_type) :: realization
  PetscInt :: icolumn
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  PetscInt :: i
  
  option => realization%option
  field => realization%field
  
  string = ''

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-T [C]"'')') icolumn
  else
    write(string2,'('',"T [C]"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-P [Pa]"'')') icolumn
  else
    write(string2,'('',"P [Pa]"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Sl"'')') icolumn
  else
    write(string2,'('',"Sl"'')')
  endif
  string = trim(string) // trim(string2)

  if (option%use_th_freezing) then
     if (icolumn > -1) then
        icolumn = icolumn + 1
        write(string2,'('',"'',i2,''-Sg"'')') icolumn
     else
        write(string2,'('',"Sg"'')')
     endif
     string = trim(string) // trim(string2)

     if (icolumn > -1) then
        icolumn = icolumn + 1
        write(string2,'('',"'',i2,''-Si"'')') icolumn
     else
        write(string2,'('',"Si"'')')
     endif
     string = trim(string) // trim(string2)

     if (icolumn > -1) then
        icolumn = icolumn + 1
        write(string2,'('',"'',i2,''-deni"'')') icolumn
     else
        write(string2,'('',"deni"'')')
     endif
     string = trim(string) // trim(string2)
  endif

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-denl"'')') icolumn
  else
    write(string2,'('',"denl"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Ul"'')') icolumn
  else
    write(string2,'('',"Ul"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-visl"'')') icolumn
  else
    write(string2,'('',"visl"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-mobilityl"'')') icolumn
  else
    write(string2,'('',"mobilityl"'')')
  endif
  string = trim(string) // trim(string2)

  do i=1,option%nflowspec
    if (icolumn > -1) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-Xl('',i2,'')"'')') icolumn,i
    else
      write(string2,'('',"Xl('',i2,'')"'')') i
    endif
    string = trim(string) // trim(string2)
  enddo
  
  THGetTecplotHeader = string

end function THGetTecplotHeader

! ************************************************************************** !

subroutine THSetPlotVariables(realization,list)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  ! 
  
  use Realization_Subsurface_class
  use Output_Aux_module
  use Variables_module
  use Material_Aux_class
  use Option_module

  implicit none

  type(realization_subsurface_type) :: realization
  type(output_variable_list_type), pointer :: list

  type(output_variable_type) :: output_variable
  character(len=MAXWORDLENGTH) :: name, units
  
  if (associated(list%first)) then
    return
  endif
  
  if (list%flow_vars) then
  
    name = 'Liquid Pressure'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                                LIQUID_PRESSURE)

    name = 'Liquid Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                LIQUID_SATURATION)
                                
    name = 'Liquid Density'
    units = 'kg/m^3'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_DENSITY)

    name = 'Liquid Energy'
    units = 'kJ/mol'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_ENERGY)

    name = 'Liquid Viscosity'
    units = 'Pa.s'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_VISCOSITY)

    name = 'Liquid Mobility'
    units = '1/Pa.s'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_MOBILITY)
  
  endif
  
  if (list%energy_vars) then

    name = 'Temperature'
    units = 'C'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                TEMPERATURE)
  
  endif

  if (realization%option%use_th_freezing) then
  
    if (realization%option%ice_model /= DALL_AMICO) then
      name = 'Gas Saturation'
      units = ''
      call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
          GAS_SATURATION)
    endif

    name = 'Ice Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
        ICE_SATURATION)

    name = 'Ice Density'
    units = 'kg/m^3'
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
        ICE_DENSITY)
        
  endif

  if (soil_compressibility_index > 0) then
  
    name = 'Transient Porosity'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                 EFFECTIVE_POROSITY)
                                 
  endif
! name = 'Phase'
! units = ''
! output_variable%iformat = 1 ! integer
! call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
!                              PHASE)

end subroutine THSetPlotVariables

! ************************************************************************** !

subroutine THComputeGradient(grid, global_auxvars, ghosted_id, gradient, &
                              option) 
  ! 
  ! Computes the gradient of temperature (for now) using
  ! least square fit of values from neighboring cells
  ! See:I. Bijelonja, I. Demirdzic, S. Muzaferija -- A finite volume method
  ! for incompressible linear elasticity, CMAME
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 2/20/12
  ! 

#include "petsc/finclude/petscdmda.h"
  use petscdmda
  use Grid_module
  use Global_Aux_module
  use Option_module
  use Utility_module

  implicit none

  type(option_type) :: option
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)

  
  PetscInt :: ghosted_neighbors_size, ghosted_id
  PetscInt :: ghosted_neighbors(26)
  PetscReal :: gradient(3), disp_vec(3,1), disp_mat(3,3)
  PetscReal :: temp_weighted(3,1)
  PetscInt :: i
  
  PetscInt :: INDX(3)
  PetscInt :: D
   
  call GridGetGhostedNeighborsWithCorners(grid,ghosted_id, &
                                         DMDA_STENCIL_STAR, &
                                         ONE_INTEGER,ONE_INTEGER,ONE_INTEGER, &
                                         ghosted_neighbors_size, &
                                         ghosted_neighbors, &
                                         option)   

  disp_vec = 0.d0
  disp_mat = 0.d0
  temp_weighted = 0.d0
  do i = 1, ghosted_neighbors_size
    disp_vec(1,1) = grid%x(ghosted_neighbors(i)) - grid%x(ghosted_id)
    disp_vec(2,1) = grid%y(ghosted_neighbors(i)) - grid%y(ghosted_id)
    disp_vec(3,1) = grid%z(ghosted_neighbors(i)) - grid%z(ghosted_id)
    disp_mat = disp_mat + matmul(disp_vec,transpose(disp_vec))
    temp_weighted = temp_weighted + disp_vec* &
                    (global_auxvars(ghosted_neighbors(i))%temp - &
                     global_auxvars(ghosted_id)%temp)
  enddo

  call ludcmp(disp_mat,THREE_INTEGER,INDX,D)
  call lubksb(disp_mat,THREE_INTEGER,INDX,temp_weighted)
  
  gradient(:) = temp_weighted(:,1)
  
end subroutine THComputeGradient

! ************************************************************************** !

subroutine THSecondaryHeat(sec_heat_vars,global_auxvar, &
                            therm_conductivity,dencpr, &
                            option,res_heat)
  ! 
  ! Calculates the source term contribution due to secondary
  ! continuum in the primary continuum residual
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/2/12
  ! 
                            
  use Option_module 
  use Global_Aux_module
  use Secondary_Continuum_Aux_module
  
  implicit none
  
  type(sec_heat_type) :: sec_heat_vars
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: temp_primary_node
  PetscReal :: m
  PetscReal :: temp_current_N
  PetscReal :: res_heat
  
  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  dm_minus = sec_heat_vars%dm_minus
  area_fm = sec_heat_vars%interfacial_area
  temp_primary_node = global_auxvar%temp
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  
  alpha = option%flow_dt*therm_conductivity/dencpr

  
! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0
                        
  rhs = sec_heat_vars%sec_temp  ! secondary continuum values from previous time step
  rhs(ngcells) = rhs(ngcells) + & 
                 alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells))* &
                 temp_primary_node
                
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    rhs(i) = rhs(i) - m*rhs(i-1)
  enddo

  ! Back substitution
  ! We only need the temperature at the outer-most node (closest to primary node)
  temp_current_N = rhs(ngcells)/coeff_diag(ngcells)
  
  ! Calculate the coupling term
  res_heat = area_fm*therm_conductivity*(temp_current_N - temp_primary_node)/ &
             dm_plus(ngcells)
                          
end subroutine THSecondaryHeat

! ************************************************************************** !

subroutine THSecondaryHeatJacobian(sec_heat_vars, &
                                    therm_conductivity, &
                                    dencpr, &
                                    option,jac_heat)
  ! 
  ! Calculates the source term jacobian contribution
  ! due to secondary continuum in the primary continuum residual
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/6/12
  ! 
                                    
  use Option_module 
  use Global_Aux_module
  use Secondary_Continuum_Aux_module
  
  implicit none
  
  type(sec_heat_type) :: sec_heat_vars
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: m
  PetscReal :: Dtemp_N_Dtemp_prim
  PetscReal :: jac_heat
  
  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  area_fm = sec_heat_vars%interfacial_area
  dm_minus = sec_heat_vars%dm_minus
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  
  alpha = option%flow_dt*therm_conductivity/dencpr

! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0
                                        
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    ! We do not have to calculate rhs terms
  enddo

  ! We need the temperature derivative at the outer-most node (closest to primary node)
  Dtemp_N_Dtemp_prim = 1.d0/coeff_diag(ngcells)*alpha*area(ngcells)/ &
                       (dm_plus(ngcells)*vol(ngcells))
  
  ! Calculate the jacobian term
  jac_heat = area_fm*therm_conductivity*(Dtemp_N_Dtemp_prim - 1.d0)/ &
             dm_plus(ngcells)
                            
              
end subroutine THSecondaryHeatJacobian                                  


! ************************************************************************** !
function THInitGuessCheck(xx, option)
  !
  ! Checks if the initial guess is valid.
  ! Note: Only implemented for DALL_AMICO formulation.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/04/2014
  !
  use Option_module

  Vec :: xx
  type(option_type), pointer :: option

  PetscInt :: THInitGuessCheck
  PetscInt :: idx
  PetscReal :: pres_min, pres_max
  PetscReal :: temp_min, temp_max
  PetscInt :: ipass, ipass0
  PetscErrorCode :: ierr

  ipass = 1

  if (option%ice_model /= DALL_AMICO) then
    THInitGuessCheck = ipass
    return
  endif

  call VecStrideMin(xx,ZERO_INTEGER,idx,pres_min,ierr)
  call VecStrideMin(xx,ONE_INTEGER ,idx,temp_min,ierr)
  call VecStrideMax(xx,ZERO_INTEGER,idx,pres_max,ierr)
  call VecStrideMax(xx,ONE_INTEGER ,idx,temp_max,ierr)

  if (pres_min < -1.d10 .or. pres_min > 1.d10 .or. &
      temp_min < -100.d0 .or. temp_max > 100.d0) then
      ipass = -1
  endif

   call MPI_Barrier(option%mycomm,ierr)
   if (option%mycommsize>1)then
      call MPI_Allreduce(ipass,ipass0,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                         option%mycomm,ierr)
      if (ipass0 < option%mycommsize) ipass=-1
   endif
   THInitGuessCheck = ipass

end function THInitGuessCheck

! ************************************************************************** !

subroutine EnergyToTemperatureBisection(T,TL,TR,h,energy,Cwi,Pr,option)
  ! 
  ! Solves the following nonlinear equation using the bisection method
  !
  ! R(T) = rho(T) Cwi hw T - energy = 0
  !
  ! Author: Nathan Collier, ORNL
  ! Date: 11/2014
  ! 
  use EOS_Water_module
  use Option_module

  implicit none

  PetscReal :: T,TL,TR,h,energy,Cwi,Pr
  type(option_type), pointer :: option

  PetscReal :: Tp,rho,rho_t,f,fR,fL,rtol
  PetscInt :: iter,niter
  PetscBool :: found
  PetscErrorCode :: ierr

  call EOSWaterdensity(TR,Pr,rho,rho_T,ierr)
  fR = rho*Cwi*h*(TR+273.15d0) - energy
  call EOSWaterdensity(TL,Pr,rho,rho_T,ierr)
  fL = rho*Cwi*h*(TL+273.15d0) - energy

  if (fL*fR > 0.d0) then
     print *,"[TL,TR] = ",TL,TR
     print *,"[fL,fR] = ",fL,fR
     write(option%io_buffer,'("th.F90: EnergyToTemperatureBisection --> root is not bracketed")')
     call printErrMsg(option)
  endif

  T = 0.5d0*(TL+TR)
  call EOSWaterdensity(T,Pr,rho,rho_T,ierr)
  f = rho*Cwi*h*(T+273.15d0) - energy

  found = PETSC_FALSE
  niter = 200
  rtol  = 1.d-6
  do iter = 1,niter
     Tp = T
     if (fL*f < 0.d0) then
        TR = T
     else 
        TL = T
     endif

     T = 0.5d0*(TL+TR)

     call EOSWaterdensity(T,Pr,rho,rho_T,ierr)
     f = rho*Cwi*h*(T+273.15d0) - energy

     if (abs((T-Tp)/(T+273.15d0)) < rtol) then
        found = PETSC_TRUE
        exit
     endif
  enddo

  if (found .eqv. PETSC_FALSE) then
     print *,"[TL,T,TR] = ",TL,T,TR
     write(option%io_buffer,'("th.F90: EnergyToTemperatureBisection --> root not found!")')
     call printErrMsg(option)
  endif

end subroutine EnergyToTemperatureBisection

! ************************************************************************** !

subroutine THUpdateSurfaceBC(realization)
  ! 
  ! Updates pressure and temperature BC associated with surface-flow
  ! 
  ! Author: Gautam Bisht
  ! Date: 10/23/13
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Region_module
  use Coupler_module
  use Connection_module
  use Fluid_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
  use String_module
  use EOS_Water_module
  use PFLOTRAN_Constants_module, only : DUMMY_VALUE,UNINITIALIZED_DOUBLE

  implicit none

  type(realization_subsurface_type) :: realization

  PetscInt :: ghosted_id
  PetscInt :: local_id
  PetscInt :: sum_connection
  PetscInt :: iconn
  PetscInt :: iter
  PetscInt :: niter
  PetscReal :: eflux
  PetscReal :: eflux_bulk
  PetscReal :: eflux_cond
  PetscReal :: area
  PetscReal :: den
  PetscReal :: den_surf_at_Told
  PetscReal :: den_subsurf
  PetscReal :: den_aveg
  PetscReal :: dum1
  PetscReal :: head_old
  PetscReal :: head_new
  PetscReal :: dhead
  PetscReal :: surfpress_old
  PetscReal :: surfpress_new
  PetscReal :: eng_per_unitvol_old
  PetscReal :: eng_per_unitvol_new
  PetscReal :: eng_times_ht_per_unitvol_old
  PetscReal :: eng_times_ht_per_unitvol_new
  PetscReal :: deng_times_ht_per_unitvol
  PetscReal :: enthalpy
  PetscReal :: surftemp_old
  PetscReal :: Temp_upwind
  PetscReal :: surftemp_new,psurftemp_new,rtol
  PetscReal :: Cwi,TL,TR,one
  PetscBool :: found
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  type(TH_auxvar_type), pointer :: auxvars(:), auxvars_bc(:)

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  global_auxvars => realization%patch%aux%Global%auxvars
  auxvars => patch%aux%TH%auxvars
  auxvars_bc => patch%aux%TH%auxvars_bc

  ! GB: Should probably add this as a member of option
  Cwi = 4.188d3 ! [J/kg/K]
  one = 1.d0

  ! Maximum no. of iterations to compute updated temperature of surface-flow
  niter = 20
  rtol  = 1.d-12
  
  eflux_bulk = 0.d0
  eflux_cond = 0.d0

  ! boundary conditions
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    if (StringCompare(boundary_condition%name,'from_surface_bc')) then

      if (boundary_condition%flow_condition%itype(TH_PRESSURE_DOF) /= &
         HET_SURF_SEEPAGE_BC) then
        call printErrMsg(option,'from_surface_bc is not of type ' // &
                        'HET_SURF_SEEPAGE_BC')
      endif

      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)

        eflux      = patch%boundary_flow_fluxes(TH_TEMPERATURE_DOF,sum_connection) ! [MJ/s]
        eflux_bulk = patch%boundary_energy_flux(1,sum_connection) ! [MJ/s]
        eflux_cond = patch%boundary_energy_flux(2,sum_connection) ! [MJ/s]

        ! [MJ/s] to [J/s]
        !geh: default internal energy units are MJ (option%scale = 1.d-6 is for J->MJ)
        eflux      = eflux/(1.d6*option%scale)
        area = cur_connection_set%area(iconn) ! [m^2]

        surfpress_old = &
          boundary_condition%flow_aux_real_var(TH_PRESSURE_DOF,iconn)
        surftemp_old = &
          boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn)
        call EOSWaterdensity(surftemp_old,option%reference_pressure,den,dum1,ierr)

        head_old = (surfpress_old - option%reference_pressure)/den/abs(option%gravity(3)) ! [m]
        dhead    = patch%boundary_velocities(1,sum_connection)*option%flow_dt ! [m]
        head_new = head_old - dhead ! [m]
        surftemp_new = UNINITIALIZED_DOUBLE ! to ensure we end up setting this

        if (head_new <= MIN_SURFACE_WATER_HEIGHT) then
          surfpress_new = option%reference_pressure
          surftemp_new = DUMMY_VALUE
        else
          
          if (head_old <= MIN_SURFACE_WATER_HEIGHT) then

            ! Surface water was absent prior to subsurface step and exfiltration
            ! occured during the subsurface step.

            surftemp_new = global_auxvars(ghosted_id)%temp
            call EOSWaterdensity(surftemp_new,option%reference_pressure,den,dum1,ierr)
            surfpress_new = head_new*(abs(option%gravity(3)))*den + &
              option%reference_pressure
          else

            ! Surface water was present prior to subsurface step

            ! Compute the new and old energy states based on the energy flux
            call EOSWaterdensity(surftemp_old,option%reference_pressure,den_surf_at_Told,dum1,ierr)

            !noc: surftemp_new is uninitialized at this point, so moving
            !these two lines below (1). Alternatively we could
            !evaluate this density at global_auxvars(ghosted_id)%temp.
            !
            !call EOSWaterdensity(surftemp_new,option%reference_pressure,den_subsurf     ,dum1,ierr)
            !den_aveg = 0.5d0*(den_surf_at_Told + den_subsurf)

            ! 1) Find new surface-temperature due to heat transfer via conduction.
            den = den_surf_at_Told
            eng_per_unitvol_old = den*Cwi*(surftemp_old + 273.15d0)
            !geh: default internal energy units are MJ (option%scale = 1.d-6 is for J->MJ)
            eng_per_unitvol_new = eng_per_unitvol_old - eflux_cond/(1.d6*option%scale)*option%flow_dt/area

            TL = -100.d0
            TR =  100.d0
            if (den*Cwi*(TL+273.15d0) < eng_per_unitvol_new) then
              surftemp_new = surftemp_old
            else
              call EnergyToTemperatureBisection(surftemp_new,TL,TR, &
                                                one, &
                                                eng_per_unitvol_new, &
                                                Cwi, &
                                                option%reference_pressure, &
                                                option)
            endif

            ! 2) Find new surface-temperature due to heat transfer via bulk-movement
            !    water transport
            if (patch%boundary_velocities(1,sum_connection) < 0) then
               Temp_upwind = global_auxvars(ghosted_id)%temp
            else
               Temp_upwind = surftemp_old
            endif

            call EOSWaterdensity(surftemp_new,option%reference_pressure,den_subsurf,dum1,ierr)
            den_aveg = 0.5d0*(den_surf_at_Told + den_subsurf)

            ! In THBCFlux():
            ! fluxe_bulk = (rho*q*H)               [kmol/m^3 * m^3/s * MJ/kmol]     = [MJ/s]
            !            = (rho*v_darcy*area*H)    [kmol/m^3 * m/s * m^2 * MJ/kmol]

            ! Retrieve H in units of [J/kg] from fluxe_bulk
            !        = [MJ/s     * J/MJ       * m^3/kg * m^{-2} * s/m]
            if (abs(patch%boundary_velocities(1,sum_connection))<1.d-14) then ! avoid division by zero
              enthalpy = 0.d0 
            else
              !geh: default internal energy units are MJ (option%scale = 1.d-6 is for J->MJ)
              enthalpy = eflux_bulk/(1.d6*option%scale)/den_aveg/area/patch%boundary_velocities(1,sum_connection)
            endif

            surftemp_old = surftemp_new
            eng_times_ht_per_unitvol_old = den     *(Cwi*surftemp_old + Cwi*273.15d0)*head_old
            deng_times_ht_per_unitvol    = den_aveg*(enthalpy         + Cwi*273.15d0)*dhead
            eng_times_ht_per_unitvol_new = eng_times_ht_per_unitvol_old - deng_times_ht_per_unitvol

            TL = -100.d0
            TR =  100.d0
            if (den*Cwi*head_new*(TL+273.15d0) < eng_times_ht_per_unitvol_new) then
              surftemp_new = surftemp_old
            else
              call EnergyToTemperatureBisection(surftemp_new,TL,TR, &
                                                head_new, &
                                                eng_times_ht_per_unitvol_new, &
                                                Cwi, &
                                                option%reference_pressure, &
                                                option)
            endif

            call EOSWaterdensity(surftemp_new,option%reference_pressure,den,dum1,ierr)
            surfpress_new = head_new*(abs(option%gravity(3)))*den + &
              option%reference_pressure
          endif
        endif

        boundary_condition%flow_aux_real_var(TH_PRESSURE_DOF,iconn) = &
          surfpress_new
        boundary_condition%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
          surftemp_new
      enddo

    else
      sum_connection = sum_connection + cur_connection_set%num_connections
    endif

    boundary_condition => boundary_condition%next

  enddo

end subroutine THUpdateSurfaceBC

! ************************************************************************** !

subroutine THUpdateSurfaceWaterFlag(realization)
  !
  ! For BC cells, set the flag for presence or absence of standing water
  !
  ! Author: Gautam Bisht
  ! Date: 04/17/14
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use String_module

  implicit none

  type(realization_subsurface_type) :: realization

  type(coupler_type), pointer :: boundary_condition
  type(TH_auxvar_type), pointer :: TH_auxvars_bc(:)
  type(TH_auxvar_type), pointer :: TH_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: ghosted_id
  PetscInt :: local_id
  PetscInt :: sum_connection
  PetscInt :: iconn

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  TH_auxvars_bc => patch%aux%TH%auxvars_bc
  TH_auxvars => patch%aux%TH%auxvars

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit

    cur_connection_set => boundary_condition%connection_set

    if (StringCompare(boundary_condition%name,'from_surface_bc')) then

      sum_connection = sum_connection + 1
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        if (patch%imat(ghosted_id) <= 0) cycle

        if (global_auxvars_bc(sum_connection)%pres(1) - option%reference_pressure < eps) then
          TH_auxvars_bc(sum_connection)%surface%surf_wat = PETSC_FALSE
          TH_auxvars(ghosted_id)%surface%surf_wat = PETSC_FALSE
        else
          TH_auxvars_bc(sum_connection)%surface%surf_wat = PETSC_TRUE
          TH_auxvars(ghosted_id)%surface%surf_wat = PETSC_TRUE
        endif
      enddo

    else

      sum_connection = sum_connection + cur_connection_set%num_connections

    endif

    boundary_condition => boundary_condition%next
  enddo

end subroutine THUpdateSurfaceWaterFlag

! ************************************************************************** !

subroutine THComputeCoeffsForSurfFlux(realization)
  !
  ! This routine computes coefficients for approximation boundary darcy
  ! flux between surface and subsurface domains.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 05/21/14
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Logging_module
  use String_module
  use EOS_Water_module
  use Material_Aux_class
  use Utility_module

  implicit none

  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(th_auxvar_type), pointer :: th_auxvars_bc(:)
  type(th_auxvar_type), pointer :: th_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(th_parameter_type), pointer :: th_parameter
  type(th_auxvar_type) :: th_auxvar_max
  type(global_auxvar_type) :: global_auxvar_max
  type(th_auxvar_type),pointer :: th_auxvar_up, th_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  !class(material_auxvar_type), pointer :: material_auxvars(:)  ! this unknownly would likely mess up %soil_properties(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvar_dn

  PetscInt :: pressure_bc_type
  PetscInt :: ghosted_id
  PetscInt :: local_id
  PetscInt :: sum_connection
  PetscInt :: iconn
  PetscInt :: icap_dn
  PetscInt :: iphase

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist(-1:3)
  PetscReal :: gravity
  PetscReal :: Dq
  PetscReal :: sir_dn
  PetscReal :: P_max_pert,P_min_pert,temp_pert
  PetscReal :: perm_dn
  PetscReal :: area
  PetscReal, pointer :: iphase_loc_p(:)
  PetscReal :: coeff_for_cubic_approx_pert(4)
  PetscReal :: range_for_linear_approx_pert(4)
  PetscReal :: slope_1, slope_2
  PetscReal :: num, den
  PetscErrorCode :: ierr
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: ithrm_loc_p(:)
  PetscInt :: ithrm_up
  PetscInt :: ithrm_dn

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  th_parameter => patch%aux%TH%th_parameter
  material_auxvars => patch%aux%Material%auxvars

  th_auxvars => patch%aux%TH%auxvars
  th_auxvars_bc => patch%aux%TH%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc

  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_yy, xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)

  ! boundary conditions
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    if (StringCompare(boundary_condition%name,'from_surface_bc')) then

      pressure_bc_type = boundary_condition%flow_condition%itype(TH_PRESSURE_DOF)

      if (pressure_bc_type /= HET_SURF_SEEPAGE_BC) then
        call printErrMsg(option,'from_surface_bc is not of type ' // &
                        'HET_SURF_SEEPAGE_BC')
      endif

      do iconn = 1, cur_connection_set%num_connections

        sum_connection = sum_connection + 1
        local_id       = cur_connection_set%id_dn(iconn)
        ghosted_id     = grid%nL2G(local_id)
        iphase         = int(iphase_loc_p(ghosted_id))

        ! Step-1: Find P_max/P_min for cubic polynomial approximation

        global_auxvar_up = global_auxvars_bc(sum_connection)
        global_auxvar_dn = global_auxvars(ghosted_id)

        th_auxvar_up => th_auxvars_bc(sum_connection)
        th_auxvar_dn => th_auxvars(ghosted_id)
        material_auxvar_dn => material_auxvars(ghosted_id)

        if (xx_p(ghosted_id*option%nflowdof-1) > 100000.d0) then
          th_auxvar_dn%surface%bcflux_default_scheme = PETSC_TRUE
        else
          th_auxvar_dn%surface%bcflux_default_scheme = PETSC_FALSE
        endif

        th_auxvar_dn%surface%coeff_for_cubic_approx(:) = -99999.d0

        dist = cur_connection_set%dist(:,iconn)

        call material_auxvar_dn%PermeabilityTensorToScalar(dist,perm_dn)

        dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
        Dq = perm_dn / dist(0)
        area = cur_connection_set%area(iconn)

        icap_dn = patch%sat_func_id(ghosted_id)
        sir_dn = th_parameter%sir(1,icap_dn)

        ithrm_up = int(ithrm_loc_p(ghosted_id))
        ithrm_dn = int(ithrm_loc_p(ghosted_id))

        ! Compute coeff
        call ComputeCoeffsForApprox(global_auxvar_up%pres(1), &
                                    global_auxvar_up%temp, &
                                    ithrm_up, &
                                    global_auxvar_dn%pres(1), &
                                    global_auxvar_dn%temp, &
                                    ithrm_dn, &
                                    material_auxvars(ghosted_id), &
                                    TH_parameter, &
                                    iphase, &
                                    patch%characteristic_curves_array(icap_dn)%ptr, &
                                    dist_gravity, &
                                    area, &
                                    Dq, &
                                    sir_dn, &
                                    option, &
                                    th_auxvar_dn%surface%P_min, &
                                    th_auxvar_dn%surface%P_max, &
                                    th_auxvar_dn%surface%coeff_for_cubic_approx, &
                                    th_auxvar_dn%surface%range_for_linear_approx)

        temp_pert = global_auxvar_dn%temp*perturbation_tolerance

        call ComputeCoeffsForApprox(global_auxvar_up%pres(1), &
                                   global_auxvar_up%temp, &
                                   ithrm_up, &
                                   global_auxvar_dn%pres(1), &
                                   global_auxvar_dn%temp + temp_pert, &
                                   ithrm_dn, &
                                   material_auxvars(ghosted_id), &
                                   TH_parameter, &
                                   iphase, &
                                   patch%characteristic_curves_array(icap_dn)%ptr, &
                                   dist_gravity, &
                                   area, &
                                   Dq, &
                                   sir_dn, &
                                   option, &
                                   P_min_pert, &
                                   P_max_pert, &
                                   coeff_for_cubic_approx_pert, &
                                   range_for_linear_approx_pert)

        th_auxvar_dn%surface%coeff_for_deriv_cubic_approx(1) = &
                            (coeff_for_cubic_approx_pert(1) - &
                             th_auxvar_dn%surface%coeff_for_cubic_approx(1))/temp_pert

        th_auxvar_dn%surface%coeff_for_deriv_cubic_approx(2) = &
                            (coeff_for_cubic_approx_pert(2) - &
                             th_auxvar_dn%surface%coeff_for_cubic_approx(2))/temp_pert

        th_auxvar_dn%surface%coeff_for_deriv_cubic_approx(3) = &
                            (coeff_for_cubic_approx_pert(3) - &
                             th_auxvar_dn%surface%coeff_for_cubic_approx(3))/temp_pert

        th_auxvar_dn%surface%coeff_for_deriv_cubic_approx(4) = &
                            (coeff_for_cubic_approx_pert(4) - &
                             th_auxvar_dn%surface%coeff_for_cubic_approx(4))/temp_pert

        num = (th_auxvar_dn%surface%range_for_linear_approx(4) - &
               th_auxvar_dn%surface%range_for_linear_approx(3))
        den = (th_auxvar_dn%surface%range_for_linear_approx(2) - &
               th_auxvar_dn%surface%range_for_linear_approx(1))
        if (abs(den) < 1.d-10) den = 1.d-10
        slope_1 = num/den

        num = (range_for_linear_approx_pert(4) - &
               range_for_linear_approx_pert(3))
        den = (range_for_linear_approx_pert(2) - &
               range_for_linear_approx_pert(1))
        if (abs(den) < 1.d-10) den = 1.d-10
        slope_2 = num/den

        th_auxvar_dn%surface%dlinear_slope_dT = (slope_2 - slope_1)/temp_pert

      enddo

    else

      sum_connection = sum_connection + cur_connection_set%num_connections

    endif

    boundary_condition => boundary_condition%next

  enddo

  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_yy, xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)

end subroutine THComputeCoeffsForSurfFlux


! ************************************************************************** !

subroutine ComputeCoeffsForApprox(P_up, T_up, ithrm_up, &
                                  P_dn, T_dn, ithrm_dn, &
                                  material_auxvar, &
                                  th_parameter, &
                                  iphase, &
                                  saturation_function, &
                                  dist_gravity, &
                                  area, &
                                  Dq, &
                                  sir_dn, &
                                  option, &
                                  P_min, P_max, &
                                  coeff_for_cubic_approx, &
                                  range_for_linear_approx)
  !
  ! To smoothly approximation boundary darcy flux, this routine computes 
  !  (i) coefficients for polynomial approximation and
  !  (ii) range for linear approximation
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 05/30/14
  !

  use EOS_Water_module
  use Field_module
  use Material_Aux_class
  use Option_module
  use Characteristic_Curves_module
  use String_module
  use Utility_module

  implicit none

  PetscReal :: P_up, T_up
  PetscInt :: ithrm_up
  PetscReal :: P_dn, T_dn
  PetscInt :: ithrm_dn
  class(material_auxvar_type) :: material_auxvar
  type(TH_parameter_type) :: th_parameter
  PetscInt :: iphase
  class(characteristic_curves_type), pointer :: saturation_function
  PetscReal :: dist_gravity
  PetscReal :: area
  PetscReal :: Dq
  PetscReal :: sir_dn
  type(option_type) :: option
  PetscReal, intent(out) :: P_min
  PetscReal, intent(out) :: P_max
  PetscReal, intent(out) :: coeff_for_cubic_approx(4)
  PetscReal, intent(out) :: range_for_linear_approx(4)

  type(global_auxvar_type) :: global_auxvar_up
  type(global_auxvar_type) :: global_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_max
  type(th_auxvar_type) :: th_auxvar_up
  type(th_auxvar_type) :: th_auxvar_dn
  type(th_auxvar_type) :: th_auxvar_max

  PetscReal :: xx(option%nflowdof)
  PetscReal :: den
  PetscReal :: dum1
  PetscReal :: upweight,gravity,dphi
  PetscReal :: ukvr
  PetscReal :: P_allowable
  PetscReal :: v_darcy_allowable,v_darcy
  PetscReal :: q_allowable,q
  PetscReal :: dq_dp_dn
  PetscReal :: dP
  PetscReal :: density_ave
  PetscReal :: dgravity_dden_dn
  PetscReal :: dukvr_dp_dn
  PetscReal :: dphi_dp_dn
  PetscReal :: perm_dn
  PetscReal :: slope
  PetscErrorCode :: ierr

  ! Distance away from allowable pressure at which cubic approximation begins
  dP = 10 ! [Pa]

  call GlobalAuxVarInit(global_auxvar_up,option)
  call GlobalAuxVarInit(global_auxvar_dn,option)
  call GlobalAuxVarInit(global_auxvar_max,option)

  call THAuxVarInit(th_auxvar_up,option)
  call THAuxVarInit(th_auxvar_dn,option)
  call THAuxVarInit(th_auxvar_max,option)

  ! Step-1: Set auxvars (global and th) for up/dn
  if (option%use_th_freezing) then

    xx(1) = P_up
    xx(2) = T_up
    call THAuxVarComputeFreezing(xx, &
                                 th_auxvar_up, &
                                 global_auxvar_up, &
                                 material_auxvar, &
                                 iphase, &
                                 saturation_function, &
                                 th_parameter, &
                                 ithrm_up, &
                                 option)

    xx(1) = P_dn
    xx(2) = T_dn
    call THAuxVarComputeFreezing(xx, &
                                 th_auxvar_dn, &
                                 global_auxvar_dn, &
                                 material_auxvar, &
                                 iphase, &
                                 saturation_function, &
                                 th_parameter, &
                                 ithrm_up, &
                                 option)
  else

    xx(1) = P_up
    xx(2) = T_up
    call THAuxVarComputeNoFreezing(xx, &
                                   th_auxvar_up, &
                                   global_auxvar_up, &
                                   material_auxvar, &
                                   iphase, &
                                   saturation_function, &
                                   th_parameter, &
                                   ithrm_up, &
                                   option)

    xx(1) = P_dn
    xx(2) = T_dn
    call THAuxVarComputeNoFreezing(xx, &
                                   th_auxvar_dn, &
                                   global_auxvar_dn, &
                                   material_auxvar, &
                                   iphase, &
                                   saturation_function, &
                                   th_parameter, &
                                   ithrm_dn, &
                                   option)
  endif

  ! Step-2: Find P_max/P_min for cubic polynomial approximation
  call EOSWaterdensity(global_auxvar_up%temp,option%reference_pressure,den,dum1,ierr)

  gravity = den * dist_gravity

  dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity

  v_darcy_allowable = (global_auxvar_up%pres(1) - option%reference_pressure)/ &
                      option%flow_dt/(-option%gravity(3))/den
  q_allowable = v_darcy_allowable*area

  if (dphi>=0.D0) then
    ukvr = th_auxvar_up%kvr
  else
    ukvr = th_auxvar_dn%kvr
  endif

  P_allowable = global_auxvar_up%pres(1) + gravity - v_darcy_allowable/Dq/ukvr

  P_max       = P_allowable + dP
  P_min       = P_allowable

  ! Step-3: Find derivative at P_max

  if (option%use_th_freezing) then

    xx(1) = P_max
    xx(2) = T_up
    call THAuxVarComputeFreezing(xx, &
                                 th_auxvar_max, &
                                 global_auxvar_max, &
                                 material_auxvar, &
                                 iphase, &
                                 saturation_function, &
                                 th_parameter, &
                                 ithrm_up, &
                                 option)
  else

    xx(1) = P_max
    xx(2) = T_up
    call THAuxVarComputeNoFreezing(xx, &
                                   th_auxvar_max, &
                                   global_auxvar_max, &
                                   material_auxvar, &
                                   iphase, &
                                   saturation_function, &
                                   th_parameter, &
                                   ithrm_up, &
                                   option)

  endif

  if (global_auxvar_up%sat(1) > sir_dn .or. global_auxvar_max%sat(1) > sir_dn) then

    upweight=1.D0
    if (global_auxvar_up%sat(1) < eps) then
      upweight=0.d0
    else if (global_auxvar_max%sat(1) < eps) then
      upweight=1.d0
    endif

    density_ave = upweight*global_auxvar_up%den(1)+(1.D0-upweight)*global_auxvar_max%den(1)

    gravity = (upweight*       global_auxvar_up%den(1) + &
               (1.D0-upweight)*global_auxvar_max%den(1)) &
              * FMWH2O * dist_gravity
    dgravity_dden_dn = (1.d0-upweight)*FMWH2O*dist_gravity

    if (option%ice_model /= DALL_AMICO) then
      dphi = global_auxvar_up%pres(1) - global_auxvar_max%pres(1) + gravity
      dphi_dp_dn = -1.d0 + dgravity_dden_dn*th_auxvar_max%dden_dp
    else
      dphi = th_auxvar_up%ice%pres_fh2o - th_auxvar_max%ice%pres_fh2o + gravity
      dphi_dp_dn = -th_auxvar_max%ice%dpres_fh2o_dp + dgravity_dden_dn*th_auxvar_max%dden_dp
    endif

    ! flow in         ! boundary cell is <= pref
    if (dphi > 0.d0 .and. global_auxvar_up%pres(1) - option%reference_pressure < eps) then
      dphi = 0.d0
      dphi_dp_dn = 0.d0
    endif

    if (dphi>=0.D0) then
      ukvr = th_auxvar_up%kvr
      dukvr_dp_dn = 0.d0
    else
      ukvr = th_auxvar_max%kvr
      dukvr_dp_dn = th_auxvar_max%dkvr_dp
    endif

    call InterfaceApprox(th_auxvar_up%kvr, th_auxvar_max%kvr, &
                         th_auxvar_up%dkvr_dp, th_auxvar_max%dkvr_dp, &
                         dphi, &
                         option%rel_perm_aveg, &
                         ukvr, dum1, dukvr_dp_dn)

    if (ukvr*Dq>floweps) then

      v_darcy = Dq * ukvr * dphi
      q = v_darcy*area

      dq_dp_dn = Dq*(dukvr_dp_dn*dphi + ukvr*dphi_dp_dn)*area

      ! Step-3: Find coefficients of cubic polynomial curve

      ! Values of function at min/max
      coeff_for_cubic_approx(1) = 0.99d0*q_allowable
      coeff_for_cubic_approx(2) = q

      ! Values of function derivatives at min/max
      slope = min(-0.01d0*q_allowable/P_min, -1.d-8)
      slope = -0.01d0*q_allowable/P_min
      coeff_for_cubic_approx(3) = slope
      coeff_for_cubic_approx(4) = dq_dp_dn

      call CubicPolynomialSetup(P_min - option%reference_pressure, &
                                P_max - option%reference_pressure, &
                                coeff_for_cubic_approx)

      ! Step-4: Save values for linear approximation
      if (Equal(q_allowable,0.d0)) then
        range_for_linear_approx(1) = 0.d0
      else
        range_for_linear_approx(1) = P_min + 0.01d0*q_allowable/slope ! - option%reference_pressure
      endif
      range_for_linear_approx(2) = P_min
      range_for_linear_approx(3) = q_allowable
      range_for_linear_approx(4) = 0.99d0*q_allowable

    endif

  endif
  
  call GlobalAuxVarStrip(global_auxvar_up)
  call GlobalAuxVarStrip(global_auxvar_dn)
  call GlobalAuxVarStrip(global_auxvar_max)

  call THAuxVarDestroy(th_auxvar_up)
  call THAuxVarDestroy(th_auxvar_dn)
  call THAuxVarDestroy(th_auxvar_max)

end subroutine ComputeCoeffsForApprox

! ************************************************************************** !

subroutine THDestroy(patch)
  ! 
  ! Deallocates variables associated with Richard
  ! 
  ! Author: ???
  ! Date: 02/14/08
  ! 

  use Patch_module

  implicit none
  
  type(patch_type) :: patch
  
  ! need to free array in aux vars
  call THAuxDestroy(patch%aux%TH)

end subroutine THDestroy

end module TH_module
