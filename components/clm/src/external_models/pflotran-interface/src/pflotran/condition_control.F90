module Condition_Control_module
  
  ! This module store routines that operate on conditions from a level above 
  ! that of the realization_module.  This is necessary to access capability
  ! such as HDF5 which is unavailable from within the realization object
  ! and below.  Routines in this module will loop over realization, levels,
  ! and patches without calling underlying level/patch versions of the
  ! subroutines, which is common in realization.F90 - GEH
#include "petsc/finclude/petscvec.h"
  use petscvec

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: CondControlAssignFlowInitCond, &
            CondControlAssignTranInitCond, &
            CondControlAssignFlowInitCondSurface, &
            CondControlScaleSourceSink
 
contains

! ************************************************************************** !

subroutine CondControlAssignFlowInitCond(realization)
  ! 
  ! Assigns flow initial conditions to model
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07, 10/18/11
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Subsurface_class
  use Discretization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Dataset_Base_class
  use Dataset_Gridded_HDF5_class
  use Dataset_Common_HDF5_class
  use Grid_module
  use Patch_module
  use EOS_Water_module

  use Global_module
  use Global_Aux_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  
  PetscInt :: icell, iconn, idof, iface
  PetscInt :: local_id, ghosted_id, iend, ibegin
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
  PetscErrorCode :: ierr
  
  character(len=MAXSTRINGLENGTH) :: string
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: initial_condition
  type(patch_type), pointer :: cur_patch
  class(dataset_base_type), pointer :: dataset
  type(global_auxvar_type) :: global_aux
  PetscBool :: use_dataset
  PetscBool :: dataset_flag(realization%option%nflowdof)
  PetscInt :: num_connections
  PetscInt, pointer :: conn_id_ptr(:)
  PetscInt :: offset, istate
  PetscReal :: x(realization%option%nflowdof)
  PetscReal :: temperature, p_sat
  PetscReal :: tempreal,pru,sou,sgu,pbu
  PetscInt :: saturated_state
  type(global_auxvar_type), pointer :: global_auxvars(:)

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch

  ! to catch uninitialized grid cells.  see VecMin check at bottom.
  call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr);CHKERRQ(ierr)
  iphase_loc_p = UNINITIALIZED_DOUBLE
  call VecRestoreArrayF90(field%iphas_loc,iphase_loc_p,ierr);CHKERRQ(ierr)

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    grid => cur_patch%grid

    select case(option%iflowmode)

      case default
        ! assign initial conditions values to domain
        call VecGetArrayF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
        call VecGetArrayF90(field%iphas_loc,iphase_loc_p,ierr);CHKERRQ(ierr)
      
        xx_p = UNINITIALIZED_DOUBLE
      
        initial_condition => cur_patch%initial_condition_list%first
        do
      
          if (.not.associated(initial_condition)) exit

          use_dataset = PETSC_FALSE
          dataset_flag = PETSC_FALSE
          do idof = 1, option%nflowdof
            dataset =>  initial_condition%flow_condition% &
                              sub_condition_ptr(idof)%ptr%dataset
            select type(dataset_ptr => dataset)
              class is(dataset_gridded_hdf5_type)
                ! already mapped to flow_aux_real_var
              class is(dataset_common_hdf5_type)
                use_dataset = PETSC_TRUE
                dataset_flag(idof) = PETSC_TRUE
                call ConditionControlMapDatasetToVec(realization, &
                        initial_condition%flow_condition% &
                          sub_condition_ptr(idof)%ptr%dataset,idof, &
                        field%flow_xx,GLOBAL)
              class default
            end select
          enddo            
          if (.not.associated(initial_condition%flow_aux_real_var) .and. &
              .not.associated(initial_condition%flow_condition)) then
            option%io_buffer = 'Flow condition is NULL in initial condition'
            call printErrMsg(option)
          endif
          if (associated(initial_condition%flow_aux_real_var)) then
            num_connections = &
              initial_condition%connection_set%num_connections
            conn_id_ptr => initial_condition%connection_set%id_dn
          else
            num_connections = initial_condition%region%num_cells
            conn_id_ptr => initial_condition%region%cell_ids
          endif
          do iconn=1, num_connections
            local_id = conn_id_ptr(iconn)
            ghosted_id = grid%nL2G(local_id)
            iend = local_id*option%nflowdof
            ibegin = iend-option%nflowdof+1
            if (cur_patch%imat(ghosted_id) <= 0) then
              xx_p(ibegin:iend) = 0.d0
              iphase_loc_p(ghosted_id) = 0
              cycle
            endif
            if (associated(initial_condition%flow_aux_real_var)) then
              do idof = 1, option%nflowdof
                if (.not.dataset_flag(idof)) then
                  xx_p(ibegin+idof-1) =  &
                    initial_condition%flow_aux_real_var(idof,iconn)
                endif
              enddo
            else
              do idof = 1, option%nflowdof
                if (.not.dataset_flag(idof)) then
                  xx_p(ibegin+idof-1) = &
                    initial_condition%flow_condition% &
                      sub_condition_ptr(idof)%ptr%dataset%rarray(1)
                endif
              enddo
            endif
            ! TODO(geh): phase out field%iphas_loc
            iphase_loc_p(ghosted_id) = &
              initial_condition%flow_condition%iphase
          enddo
          initial_condition => initial_condition%next
        enddo
        call VecRestoreArrayF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)

    end select 
   
    cur_patch => cur_patch%next
  enddo

  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
                                   field%flow_xx_loc,NFLOWDOF)  

  call VecCopy(field%flow_xx, field%flow_yy, ierr);CHKERRQ(ierr)
  call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                  field%iphas_loc,ONEDOF)  
  call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                  field%iphas_old_loc,ONEDOF)

  ! cannot perform VecMin on local vector as the ghosted corner values are not
  ! updated during the local to local update.
  call DiscretizationLocalToGlobal(discretization,field%iphas_loc,field%work, &
                                   ONEDOF)
  call VecMin(field%work,PETSC_NULL_INTEGER,tempreal,ierr);CHKERRQ(ierr)
  if (tempreal < 0.d0) then
!    print *, tempreal
    option%io_buffer = 'Uninitialized cells in domain.'
    call printErrMsg(option)
  endif

end subroutine CondControlAssignFlowInitCond

! ************************************************************************** !

subroutine CondControlAssignTranInitCond(realization)
  ! 
  ! Assigns transport initial conditions to model
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07, 10/18/11
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec

  use Realization_Subsurface_class
  use Discretization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Condition_module
  use Transport_Constraint_module
  use Grid_module
  use Dataset_Base_class
  use Patch_module
  use Reactive_Transport_module, only : RTUpdateAuxVars, &
                                        RTUpdateActivityCoefficients
  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Reaction_module
  use HDF5_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  PetscInt :: icell, iconn, idof, isub_condition, temp_int, iimmobile
  PetscInt :: local_id, ghosted_id, iend, ibegin
  PetscInt :: irxn, isite, imnrl, ikinrxn
  PetscReal, pointer :: xx_p(:), xx_loc_p(:), vec_p(:), vec_p2(:)
  Vec :: vec1_loc
  Vec :: vec2_loc
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: initial_condition
  type(patch_type), pointer :: cur_patch
  type(reaction_type), pointer :: reaction
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(tran_constraint_coupler_type), pointer :: constraint_coupler
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: iphase
  PetscInt :: offset
  PetscBool :: re_equilibrate_at_each_cell
  character(len=MAXSTRINGLENGTH) :: string, string2
  class(dataset_base_type), pointer :: dataset
  PetscInt :: aq_dataset_to_idof(realization%reaction%naqcomp)
  PetscInt :: iaqdataset, num_aq_datasets
  PetscBool :: use_aq_dataset
  PetscReal :: ave_num_iterations
  PetscReal :: tempreal
  PetscInt :: prev_equilibrated_ghosted_id
  PetscReal, pointer :: iphase_loc_p(:)
  PetscReal, pointer :: flow_xx_p(:)
  PetscLogDouble :: tstart, tend
  
  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  reaction => realization%reaction
  
  iphase = 1
  vec1_loc = PETSC_NULL_VEC
  vec2_loc = PETSC_NULL_VEC
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    grid => cur_patch%grid
    rt_auxvars => cur_patch%aux%RT%auxvars
    global_auxvars => cur_patch%aux%Global%auxvars
    material_auxvars => cur_patch%aux%Material%auxvars

    ! assign initial conditions values to domain
    call VecGetArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
      
    xx_p = UNINITIALIZED_DOUBLE
      
    initial_condition => cur_patch%initial_condition_list%first
    do
      
      if (.not.associated(initial_condition)) exit
        
      constraint_coupler => initial_condition%tran_condition%cur_constraint_coupler

      re_equilibrate_at_each_cell = PETSC_FALSE
      use_aq_dataset = PETSC_FALSE
      num_aq_datasets = 0
      aq_dataset_to_idof = 0
      do idof = 1, reaction%naqcomp ! primary aqueous concentrations
        if (constraint_coupler%aqueous_species%external_dataset(idof)) then
          num_aq_datasets = num_aq_datasets + 1
          aq_dataset_to_idof(num_aq_datasets) = idof
          re_equilibrate_at_each_cell = PETSC_TRUE
          use_aq_dataset = PETSC_TRUE
          string = 'constraint ' // trim(constraint_coupler%constraint_name)
          dataset => DatasetBaseGetPointer(realization%datasets, &
                        constraint_coupler%aqueous_species%constraint_aux_string(idof), &
                        string,option)
          call ConditionControlMapDatasetToVec(realization,dataset,idof, &
                                                field%tran_xx_loc,LOCAL)
        endif
      enddo

      ! read in heterogeneous mineral volume fractions
      if (associated(constraint_coupler%minerals)) then
        do imnrl = 1, reaction%mineral%nkinmnrl
          if (constraint_coupler%minerals%external_vol_frac_dataset(imnrl)) then
            re_equilibrate_at_each_cell = PETSC_TRUE
            string = 'constraint ' // trim(constraint_coupler%constraint_name)
            dataset => DatasetBaseGetPointer(realization%datasets, &
                          constraint_coupler%minerals% &
                            constraint_vol_frac_string(imnrl), &
                          string,option)
            if (vec1_loc == PETSC_NULL_VEC) then
              ! cannot use field%work_loc as it is used within ConditionCo...
              call VecDuplicate(field%work_loc,vec1_loc, &
                                ierr);CHKERRQ(ierr)
            endif 
            idof = ONE_INTEGER
            call ConditionControlMapDatasetToVec(realization,dataset,idof, &
                                                 vec1_loc,LOCAL)
            call VecGetArrayF90(vec1_loc,vec_p,ierr);CHKERRQ(ierr)
            do icell=1,initial_condition%region%num_cells
              local_id = initial_condition%region%cell_ids(icell)
              ghosted_id = grid%nL2G(local_id)
              rt_auxvars(ghosted_id)%mnrl_volfrac0(imnrl) = vec_p(ghosted_id)
              rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl) = vec_p(ghosted_id)
            enddo
            call VecRestoreArrayF90(vec1_loc,vec_p,ierr);CHKERRQ(ierr)
          endif
        enddo
      endif
          
      ! read in heterogeneous mineral surface area
      if (associated(constraint_coupler%minerals)) then
        do imnrl = 1, reaction%mineral%nkinmnrl
          if (constraint_coupler%minerals%external_area_dataset(imnrl)) then
            re_equilibrate_at_each_cell = PETSC_TRUE
            string = 'constraint ' // trim(constraint_coupler%constraint_name)
            dataset => DatasetBaseGetPointer(realization%datasets, &
                          constraint_coupler%minerals% &
                          constraint_area_string(imnrl), &
                          string,option)
            if (vec1_loc == PETSC_NULL_VEC) then
              ! cannot use field%work_loc as it is used within ConditionCo...
              call VecDuplicate(field%work_loc,vec1_loc, &
                                ierr);CHKERRQ(ierr)
            endif 
            idof = ONE_INTEGER
            call ConditionControlMapDatasetToVec(realization,dataset,idof, &
                                                 vec1_loc,LOCAL)
            call VecScale(vec1_loc,constraint_coupler%minerals% &
                            constraint_area_conv_factor(imnrl), &
                          ierr);CHKERRQ(ierr)
            if (constraint_coupler%minerals%area_per_unit_mass(imnrl)) then
              if (constraint_coupler%minerals% &
                    external_vol_frac_dataset(imnrl)) then
                dataset => DatasetBaseGetPointer(realization%datasets, &
                              constraint_coupler%minerals% &
                                constraint_vol_frac_string(imnrl), &
                              string,option)
                if (vec2_loc == PETSC_NULL_VEC) then
                  call VecDuplicate(vec1_loc,vec2_loc, &
                                    ierr);CHKERRQ(ierr)
                endif 
                idof = ONE_INTEGER
                call ConditionControlMapDatasetToVec(realization,dataset, &
                                                     idof,vec2_loc,LOCAL)
                call VecPointwiseMult(vec1_loc,vec1_loc, &
                                      vec2_loc,ierr);CHKERRQ(ierr)
              else
                call VecScale(vec1_loc, &
                              constraint_coupler%minerals% &
                                constraint_vol_frac(imnrl), &
                              ierr);CHKERRQ(ierr)
              endif
            endif
            call VecGetArrayF90(vec1_loc,vec_p,ierr);CHKERRQ(ierr)
            do icell=1,initial_condition%region%num_cells
              local_id = initial_condition%region%cell_ids(icell)
              ghosted_id = grid%nL2G(local_id)
              rt_auxvars(ghosted_id)%mnrl_area0(imnrl) = vec_p(ghosted_id)
              rt_auxvars(ghosted_id)%mnrl_area(imnrl) = vec_p(ghosted_id)
            enddo
            call VecRestoreArrayF90(vec1_loc,vec_p,ierr);CHKERRQ(ierr)
          endif
        enddo
      endif
          
      ! read in heterogeneous immobile
      if (associated(constraint_coupler%immobile_species)) then
        do iimmobile = 1, reaction%immobile%nimmobile
          if (constraint_coupler%immobile_species%external_dataset(iimmobile)) then
            ! no need to requilibrate at each cell
            string = 'constraint ' // trim(constraint_coupler%constraint_name)
            dataset => DatasetBaseGetPointer(realization%datasets, &
                constraint_coupler%immobile_species%constraint_aux_string(iimmobile), &
                string,option)
            if (vec1_loc == PETSC_NULL_VEC) then
              ! cannot use field%work_loc as it is used within ConditionCo...
              call VecDuplicate(field%work_loc,vec1_loc, &
                                ierr);CHKERRQ(ierr)
            endif 
            idof = ONE_INTEGER
            call ConditionControlMapDatasetToVec(realization,dataset,idof, &
                                                 vec1_loc,LOCAL)
            call VecGetArrayF90(vec1_loc,vec_p,ierr);CHKERRQ(ierr)
            do icell=1,initial_condition%region%num_cells
              local_id = initial_condition%region%cell_ids(icell)
              ghosted_id = grid%nL2G(local_id)
              rt_auxvars(ghosted_id)%immobile(iimmobile) = vec_p(ghosted_id)
            enddo
            call VecRestoreArrayF90(vec1_loc,vec_p,ierr);CHKERRQ(ierr)
          endif
        enddo
      endif
          
      if (.not.option%use_isothermal) then
        re_equilibrate_at_each_cell = PETSC_TRUE
      endif
        
      if (use_aq_dataset) then
        call VecGetArrayF90(field%tran_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
        call PetscTime(tstart,ierr);CHKERRQ(ierr)
      endif
        
      ave_num_iterations = 0.d0
      prev_equilibrated_ghosted_id = 0
      do icell=1,initial_condition%region%num_cells
        local_id = initial_condition%region%cell_ids(icell)
        ghosted_id = grid%nL2G(local_id)
        iend = local_id*option%ntrandof
        ibegin = iend-option%ntrandof+1
        if (cur_patch%imat(ghosted_id) <= 0) then
          xx_p(ibegin:iend) = 1.d-200
          cycle
        endif
        if (re_equilibrate_at_each_cell) then
          if (use_aq_dataset) then
            offset = (ghosted_id-1)*option%ntrandof
            do iaqdataset = 1, num_aq_datasets
              ! remember that xx_loc_p holds the data set values that were read in
              temp_int = aq_dataset_to_idof(iaqdataset)
              constraint_coupler%aqueous_species%constraint_conc(temp_int) = &
                xx_loc_p(offset+temp_int)
            enddo
          endif
          option%iflag = grid%nG2A(grid%nL2G(local_id))
          if (prev_equilibrated_ghosted_id == 0) then
            call ReactionEquilibrateConstraint(rt_auxvars(ghosted_id), &
              global_auxvars(ghosted_id),material_auxvars(ghosted_id), &
              reaction, &
              constraint_coupler%constraint_name, &
              constraint_coupler%aqueous_species, &
              constraint_coupler%free_ion_guess, &
              constraint_coupler%minerals, &
              constraint_coupler%colloids, &
              constraint_coupler%immobile_species, &
              constraint_coupler%num_iterations, &
              PETSC_FALSE,option)
          else
            ! copy molalities from previous equilibrated auxvar as initial guess
            rt_auxvars(ghosted_id)%pri_molal = &
              rt_auxvars(prev_equilibrated_ghosted_id)%pri_molal
            call ReactionEquilibrateConstraint(rt_auxvars(ghosted_id), &
              global_auxvars(ghosted_id),material_auxvars(ghosted_id), &
              reaction, &
              constraint_coupler%constraint_name, &
              constraint_coupler%aqueous_species, &
              constraint_coupler%free_ion_guess, &
              constraint_coupler%minerals, &
              constraint_coupler%colloids, &
              constraint_coupler%immobile_species, &
              constraint_coupler%num_iterations, &
              PETSC_TRUE,option)
          endif
          option%iflag = 0
          ave_num_iterations = ave_num_iterations + &
            constraint_coupler%num_iterations
          ! update CO2 mole fraction for CO2 modes

          ! prev_eq_ghosted_id is only updated to the prev active cell
          prev_equilibrated_ghosted_id = ghosted_id
        endif
        ! ibegin is the local non-ghosted offset: (local_id-1)*option%ntrandof+1
        offset = ibegin + reaction%offset_aqueous - 1
        ! primary aqueous concentrations
        do idof = 1, reaction%naqcomp 
          xx_p(offset+idof) = &
            constraint_coupler%aqueous_species%basis_molarity(idof) / &
            global_auxvars(ghosted_id)%den_kg(iphase)*1000.d0 ! convert molarity -> molality
        enddo
        ! mineral volume fractions
        if (associated(constraint_coupler%minerals)) then
          do imnrl = 1, reaction%mineral%nkinmnrl
            ! if read from a dataset, the vol frac was set above.  Don't want to
            ! overwrite
            if (.not.constraint_coupler%minerals% &
                  external_vol_frac_dataset(imnrl)) then
              rt_auxvars(ghosted_id)%mnrl_volfrac0(imnrl) = &
                constraint_coupler%minerals%constraint_vol_frac(imnrl)
              rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl) = &
                constraint_coupler%minerals%constraint_vol_frac(imnrl)
            endif
            if (.not.constraint_coupler%minerals% &
                  external_area_dataset(imnrl)) then
              rt_auxvars(ghosted_id)%mnrl_area0(imnrl) = &
                constraint_coupler%minerals%constraint_area(imnrl)
              rt_auxvars(ghosted_id)%mnrl_area(imnrl) = &
                constraint_coupler%minerals%constraint_area(imnrl)
            endif
          enddo
        endif
        ! colloids fractions
        if (associated(constraint_coupler%colloids)) then
          offset = ibegin + reaction%offset_colloid - 1
          do idof = 1, reaction%ncoll ! primary aqueous concentrations
            xx_p(offset+idof) = &
              constraint_coupler%colloids%basis_conc_mob(idof) / &
              global_auxvars(ghosted_id)%den_kg(iphase)*1000.d0 ! convert molarity -> molality
            rt_auxvars(ghosted_id)%colloid%conc_imb(idof) = &
              constraint_coupler%colloids%basis_conc_imb(idof)
          enddo
        endif
        ! immobile
        if (associated(constraint_coupler%immobile_species)) then
          offset = ibegin + reaction%offset_immobile - 1
          do iimmobile = 1, reaction%immobile%nimmobile
            if (constraint_coupler%immobile_species%external_dataset(iimmobile)) then
              ! already read into rt_auxvars above.
              xx_p(offset+iimmobile) = &
                rt_auxvars(ghosted_id)%immobile(iimmobile)
            else
              xx_p(offset+iimmobile) = &
                constraint_coupler%immobile_species%constraint_conc(iimmobile)
              rt_auxvars(ghosted_id)%immobile(iimmobile) = &
                constraint_coupler%immobile_species%constraint_conc(iimmobile)
            endif
          enddo
        endif
      enddo ! icell=1,initial_condition%region%num_cells
      if (use_aq_dataset) then
        call PetscTime(tend,ierr);CHKERRQ(ierr)
        call VecRestoreArrayF90(field%tran_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
        ave_num_iterations = ave_num_iterations / &
          initial_condition%region%num_cells
        write(option%io_buffer,&
              '("Average number of iterations in ReactionEquilibrateConstraint():", &
              & f5.1)') ave_num_iterations
        call printMsg(option)
        write(option%io_buffer,'(f10.2," Seconds to equilibrate constraints")') &
          tend-tstart
        call printMsg(option)
      endif
      initial_condition => initial_condition%next
    enddo
      
    call VecRestoreArrayF90(field%tran_xx,xx_p, ierr);CHKERRQ(ierr)

    cur_patch => cur_patch%next
  enddo
  
  ! check to ensure that minimum concentration is not less than or equal
  ! to zero
  call VecMin(field%tran_xx,PETSC_NULL_INTEGER,tempreal,ierr);CHKERRQ(ierr)
  if (tempreal <= 0.d0) then
    option%io_buffer = 'ERROR: Zero concentrations found in initial ' // &
      'transport solution.'
    call printMsg(option)
    ! now figure out which species have zero concentrations
    do idof = 1, option%ntrandof
      call VecStrideMin(field%tran_xx,idof-1,offset,tempreal, &
                        ierr);CHKERRQ(ierr)
      if (tempreal <= 0.d0) then
        write(string,*) tempreal
        if (idof <= reaction%naqcomp) then
          string2 = '  Aqueous species "' // &
            trim(reaction%primary_species_names(idof))
        else
          string2 = '  Immobile species "' // &
            trim(reaction%immobile%names(idof-reaction%offset_immobile))
        endif
          string2 = trim(string2) // &
            '" has zero concentration (' // &
            trim(adjustl(string)) // ').'
        call printMsg(option,string2)
      endif
    enddo
    option%io_buffer = ''
    call printMsg(option)
    option%io_buffer = '*** Begin Note'
    call printMsg(option)
    option%io_buffer = 'If concentrations = -999., they have not ' // &
              'been initialized properly.'
    call printMsg(option)
    option%io_buffer = '*** End Note'
    call printMsg(option)
    option%io_buffer = 'Free ion concentations must be positive.  Try ' // &
      'using a small value such as 1.e-20 or 1.e-40 instead of zero.'
    call printErrMsg(option)
  endif
  
  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)  
  call VecCopy(field%tran_xx, field%tran_yy, ierr);CHKERRQ(ierr)

  ! override initial conditions if they are to be read from a file
  if (len_trim(option%initialize_transport_filename) > 1) then
    call CondControlReadTransportIC(realization, &
                                    option%initialize_transport_filename)
  endif
  call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_FALSE,PETSC_FALSE)
  ! at this point the auxvars have been computed with activity coef = 1.d0
  ! to use intitial condition with activity coefs /= 1.d0, must update
  ! activity coefs and recompute auxvars
  if (realization%reaction%act_coef_update_frequency /= &
      ACT_COEF_FREQUENCY_OFF) then
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_FALSE,PETSC_TRUE)
    !geh: you may ask, why call this twice....  We need to iterate at least
    !     once to ensure that the activity coefficients are more accurate.
    !     Otherwise, the total component concentrations can be quite
    !     different from what is defined in the input file.
    call RTUpdateAuxVars(realization,PETSC_TRUE,PETSC_FALSE,PETSC_TRUE)
  endif

  if (vec1_loc /= PETSC_NULL_VEC) then
    call VecDestroy(vec1_loc,ierr);CHKERRQ(ierr)
  endif
  if (vec2_loc /= PETSC_NULL_VEC) then
    call VecDestroy(vec2_loc,ierr);CHKERRQ(ierr)
  endif

end subroutine CondControlAssignTranInitCond

! ************************************************************************** !

subroutine ConditionControlMapDatasetToVec(realization,dataset,idof, &
                                           mdof_vec,vec_type)
  ! 
  ! maps an external dataset to a PETSc vec
  ! representing values at each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/23/12
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Dataset_Common_HDF5_class
  use Dataset_Base_class
  use HDF5_module
  use Discretization_module

  implicit none
  

  class(realization_subsurface_type) :: realization
  class(dataset_base_type), pointer :: dataset
  PetscInt :: idof
  Vec :: mdof_vec
  PetscInt :: vec_type

  type(field_type), pointer :: field
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr

  field => realization%field
  option => realization%option
  
  call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
  if (associated(dataset)) then
    select type(dataset)
      class is (dataset_common_hdf5_type)
        string = '' ! group name
        ! have to copy to string2 due to mismatch in string size
        string2 = dataset%hdf5_dataset_name
        call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                          dataset%filename, &
                                          string,string2, &
                                          dataset%realization_dependent)
        if (vec_type == GLOBAL) then
          call VecStrideScatter(field%work,idof-1,mdof_vec, &
                                INSERT_VALUES,ierr);CHKERRQ(ierr)
        else
          call DiscretizationGlobalToLocal(realization%discretization, &
                                           field%work, &
                                           field%work_loc,ONEDOF)
          call VecStrideScatter(field%work_loc,idof-1,mdof_vec, &
                                INSERT_VALUES,ierr);CHKERRQ(ierr)
        endif
      class default
        option%io_buffer = 'Dataset "' // trim(dataset%name) // &
          '" not supported in ConditionControlMapDatasetToVec.'
        call printErrMsg(option)
    end select
  endif

end subroutine ConditionControlMapDatasetToVec

! ************************************************************************** !

subroutine CondControlScaleSourceSink(realization)
  ! 
  ! Scales select source/sinks based on perms
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/03/08, 10/18/11
  ! 
#include "petsc/finclude/petscdmda.h"
  use petscdmda
      
  use Realization_Subsurface_class
  use Discretization_module
  use Region_module
  use Option_module
  use Field_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Grid_module
  use Patch_module
  use Material_Aux_class
  use Variables_module, only : PERMEABILITY_X

  implicit none
  
  class(realization_subsurface_type) :: realization
  
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: cur_source_sink
  type(connection_set_type), pointer :: cur_connection_set
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(patch_type), pointer :: cur_patch
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id, neighbor_ghosted_id
  PetscInt :: iconn
  PetscReal :: scale, sum
  PetscInt :: icount
  PetscInt :: x_count, y_count, z_count
  PetscInt, parameter :: x_width = 1, y_width = 1, z_width = 0
  
  PetscInt :: ghosted_neighbors(0:27)
  
  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  material_auxvars => realization%patch%aux%Material%auxvars
  
  ! GB: grid was uninitialized
  grid => patch%grid

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    ! BIG-TIME warning here.  I assume that all source/sink cells are within 
    ! a single patch - geh

    grid => cur_patch%grid

    cur_source_sink => cur_patch%source_sink_list%first
    do
      if (.not.associated(cur_source_sink)) exit

      call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)

      cur_connection_set => cur_source_sink%connection_set
    
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)

        select case(option%iflowmode)
          case(TH_MODE)
              call GridGetGhostedNeighbors(grid,ghosted_id,DMDA_STENCIL_STAR, &
                                          x_width,y_width,z_width, &
                                          x_count,y_count,z_count, &
                                          ghosted_neighbors,option)
              ! ghosted neighbors is ordered first in x, then, y, then z
              icount = 0
              sum = 0.d0
              ! x-direction
              do while (icount < x_count)
                icount = icount + 1
                neighbor_ghosted_id = ghosted_neighbors(icount)
                sum = sum + MaterialAuxVarGetValue(material_auxvars( &
                              neighbor_ghosted_id),PERMEABILITY_X) * &
                            grid%structured_grid%dy(neighbor_ghosted_id)* &
                            grid%structured_grid%dz(neighbor_ghosted_id)
                 
              enddo
              ! y-direction
              do while (icount < x_count + y_count)
                icount = icount + 1
                neighbor_ghosted_id = ghosted_neighbors(icount)                 
                sum = sum + MaterialAuxVarGetValue(material_auxvars( &
                              neighbor_ghosted_id),PERMEABILITY_X) * &
                            grid%structured_grid%dx(neighbor_ghosted_id)* &
                            grid%structured_grid%dz(neighbor_ghosted_id)
                 
              enddo
              ! z-direction
              do while (icount < x_count + y_count + z_count)
                icount = icount + 1
                neighbor_ghosted_id = ghosted_neighbors(icount)                 
                sum = sum + MaterialAuxVarGetValue(material_auxvars( &
                              neighbor_ghosted_id),PERMEABILITY_X) * &
                            grid%structured_grid%dx(neighbor_ghosted_id)* &
                            grid%structured_grid%dy(neighbor_ghosted_id)
              enddo
              vec_ptr(local_id) = vec_ptr(local_id) + sum
         end select

      enddo
        
      call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
      call VecNorm(field%work,NORM_1,scale,ierr);CHKERRQ(ierr)
      scale = 1.d0/scale
      call VecScale(field%work,scale,ierr);CHKERRQ(ierr)

      call VecGetArrayF90(field%work,vec_ptr, ierr);CHKERRQ(ierr)
      do iconn = 1, cur_connection_set%num_connections      
        local_id = cur_connection_set%id_dn(iconn)
        select case(option%iflowmode)
          case(TH_MODE)
            cur_source_sink%flow_aux_real_var(ONE_INTEGER,iconn) = &
              vec_ptr(local_id)
        end select

      enddo
      call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
        
      cur_source_sink => cur_source_sink%next
    enddo
    cur_patch => cur_patch%next
  enddo

end subroutine CondControlScaleSourceSink

! ************************************************************************** !

subroutine CondControlReadTransportIC(realization,filename)
  ! 
  ! Assigns transport initial condition from
  ! HDF5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/05/10
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use Patch_module
  use Reactive_Transport_module
  use Reaction_Aux_module
  use Discretization_module
  use HDF5_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  character(len=MAXSTRINGLENGTH) :: filename
  
  PetscInt :: local_id, idx, offset, idof
  PetscReal, pointer :: xx_p(:)
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscReal, pointer :: vec_p(:)  
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(patch_type), pointer :: cur_patch
  type(reaction_type), pointer :: reaction

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  reaction => realization%reaction

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    grid => cur_patch%grid

      ! assign initial conditions values to domain
    call VecGetArrayF90(field%tran_xx,xx_p, ierr);CHKERRQ(ierr)

    ! Primary species concentrations for all modes 
    do idof = 1, option%ntrandof ! primary aqueous concentrations
      offset = idof
      group_name = ''
      dataset_name = reaction%primary_species_names(idof)
      call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                        filename,group_name, &
                                        dataset_name,option%id>0)
      call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
      do local_id=1, grid%nlmax
        if (cur_patch%imat(grid%nL2G(local_id)) <= 0) cycle
        if (vec_p(local_id) < 1.d-40) then
          print *,  option%myrank, grid%nG2A(grid%nL2G(local_id)), &
            ': Zero free-ion concentration in Initial Condition read from file.'
        endif
        idx = (local_id-1)*option%ntrandof + offset
        xx_p(idx) = vec_p(local_id)
      enddo
      call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
     
    enddo     

    call VecRestoreArrayF90(field%tran_xx,xx_p, ierr);CHKERRQ(ierr)
        
    cur_patch => cur_patch%next
  enddo
   
  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                   field%tran_xx_loc,NTRANDOF)  
  call VecCopy(field%tran_xx, field%tran_yy, ierr);CHKERRQ(ierr)
  
end subroutine CondControlReadTransportIC

! ************************************************************************** !

subroutine CondControlAssignFlowInitCondSurface(surf_realization)

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Surface_class
  use Discretization_module
  use Region_module
  use Option_module
  use Surface_Field_module
  use Coupler_module
  use Condition_module
  use Grid_module
  use Patch_module
  use EOS_Water_module
  use Surface_TH_Aux_module
  use Surface_Global_Aux_module
  
  implicit none
  
  class(realization_surface_type) :: surf_realization
  
  PetscInt :: icell, iconn, idof, iface
  PetscInt :: local_id, ghosted_id, iend, ibegin
  PetscReal, pointer :: xx_p(:)!, iphase_loc_p(:)
  PetscErrorCode :: ierr
  
  PetscReal :: temperature, p_sat
  PetscReal :: pw, dw_kg, dw_mol, hw
  PetscReal :: temp
  PetscReal :: dpsat_dt
  character(len=MAXSTRINGLENGTH) :: string
  
  type(option_type), pointer :: option
  type(surface_field_type), pointer :: surf_field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(coupler_type), pointer :: initial_condition
  type(patch_type), pointer :: cur_patch
  type(Surface_TH_auxvar_type), pointer :: surf_th_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)

  option => surf_realization%option
  discretization => surf_realization%discretization
  surf_field => surf_realization%surf_field
  patch => surf_realization%patch

  if (option%iflowmode == TH_MODE) then
    surf_th_auxvars => patch%surf_aux%SurfaceTH%auxvars
    surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  endif

  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    grid => cur_patch%grid

    select case(option%iflowmode)
      
      case (TH_MODE)
        ! assign initial conditions values to domain
        call VecGetArrayF90(surf_field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
    
        xx_p = UNINITIALIZED_DOUBLE
      
        initial_condition => cur_patch%initial_condition_list%first
        do
      
          if (.not.associated(initial_condition)) exit

            if (.not.associated(initial_condition%flow_aux_real_var)) then
              if (.not.associated(initial_condition%flow_condition)) then
                option%io_buffer = 'Flow condition is NULL in initial condition'
                call printErrMsg(option)
              endif
              do icell=1,initial_condition%region%num_cells
                local_id = initial_condition%region%cell_ids(icell)
                ghosted_id = grid%nL2G(local_id)
                iend = local_id*option%nflowdof
                ibegin = iend-option%nflowdof+1
                if (cur_patch%imat(ghosted_id) <= 0) then
                  xx_p(ibegin:iend) = 0.d0
                  cycle
                endif
                do idof = 1, option%nflowdof
                  select case (idof)
                    case (ONE_INTEGER)
                      xx_p(ibegin+idof-1) = &
                        initial_condition%flow_condition% &
                        sub_condition_ptr(idof)%ptr%dataset%rarray(1)
                    case (TWO_INTEGER)
                      temp = &
                        initial_condition%flow_condition% &
                        sub_condition_ptr(idof)%ptr%dataset%rarray(1)
                      pw = option%reference_pressure
                        
                      call EOSWaterDensity(temp,pw,dw_kg,dw_mol,ierr)
                      ! [rho*h*T*Cw]
                      xx_p(ibegin+idof-1) = dw_kg*xx_p(ibegin)* &
                                            (temp + 273.15d0)* &
                                            surf_th_auxvars(ghosted_id)%Cw
                      surf_global_auxvars(ghosted_id)%den_kg(1) = dw_kg
                      surf_global_auxvars(ghosted_id)%temp = temp
                  end select
                enddo
                xx_p(ibegin:iend) = 0.0d0
              enddo
            else
              do iconn=1,initial_condition%connection_set%num_connections
                local_id = initial_condition%connection_set%id_dn(iconn)
                ghosted_id = grid%nL2G(local_id)
                iend = local_id*option%nflowdof
                ibegin = iend-option%nflowdof+1
                if (cur_patch%imat(ghosted_id) <= 0) then
                  xx_p(ibegin:iend) = 0.d0
                  cycle
                endif
                xx_p(ibegin:iend) = &
                      initial_condition%flow_aux_real_var(1:option%nflowdof,iconn)
                xx_p(ibegin:iend) = 0.0d0
              enddo
            endif
          initial_condition => initial_condition%next
        enddo
     
        call VecRestoreArrayF90(surf_field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
      case default
        option%io_buffer = 'CondControlAssignFlowInitCondSurface not ' // &
          'for this mode'
        call printErrMsg(option)
    end select 
   
    cur_patch => cur_patch%next
  enddo
   
  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization, surf_field%flow_xx, &
                                   surf_field%flow_xx_loc, NFLOWDOF)

end subroutine CondControlAssignFlowInitCondSurface

end module Condition_Control_module
