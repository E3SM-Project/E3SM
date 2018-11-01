module PMC_Subsurface_class

  use PMC_Base_class
  use Realization_Subsurface_class

  use PFLOTRAN_Constants_module

  implicit none

#include "petsc/finclude/petscsys.h"
  
  private

  type, public, extends(pmc_base_type) :: pmc_subsurface_type
    class(realization_subsurface_type), pointer :: realization
  contains
    procedure, public :: Init => PMCSubsurfaceInit
    procedure, public :: SetupSolvers => PMCSubsurfaceSetupSolvers
    procedure, public :: GetAuxData => PMCSubsurfaceGetAuxData
    procedure, public :: SetAuxData => PMCSubsurfaceSetAuxData
    procedure, public :: Destroy => PMCSubsurfaceDestroy
  end type pmc_subsurface_type
  
  public :: PMCSubsurfaceCreate
  
contains

! ************************************************************************** !

function PMCSubsurfaceCreate()
  ! 
  ! Allocates and initializes a new process_model_coupler
  ! object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pmc_subsurface_type), pointer :: PMCSubsurfaceCreate
  
  class(pmc_subsurface_type), pointer :: pmc

#ifdef DEBUG
  print *, 'PMCSubsurface%Create()'
#endif
  
  allocate(pmc)
  call pmc%Init()
  
  PMCSubsurfaceCreate => pmc  
  
end function PMCSubsurfaceCreate

! ************************************************************************** !

subroutine PMCSubsurfaceInit(this)
  ! 
  ! Initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  ! 
  ! for some reason, Intel with VS want this explicitly specified.
  use PMC_Base_class, only : PMCBaseInit
  
  implicit none
  
  class(pmc_subsurface_type) :: this
  
#ifdef DEBUG
  print *, 'PMCSubsurface%Init()'
#endif
  
  call PMCBaseInit(this)
  this%name = 'PMCSubsurface'
  nullify(this%realization)

end subroutine PMCSubsurfaceInit

! ************************************************************************** !

subroutine PMCSubsurfaceSetupSolvers(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Convergence_module
  use Discretization_module
  use Option_module
  use PMC_Base_class
  use PM_Base_Pointer_module
  use PM_Base_class
  use PM_Subsurface_Flow_class

  use PM_TH_class
  use PM_RT_class
  use Secondary_Continuum_module, only : SecondaryRTUpdateIterate  
  use Solver_module
  use Timestepper_Base_class
  use Timestepper_BE_class

  implicit none

  class(pmc_subsurface_type) :: this

  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  SNESLineSearch :: linesearch  
  character(len=MAXSTRINGLENGTH) :: string  
  PetscErrorCode :: ierr

#ifdef DEBUG
  call printMsg(this%option,'PMCSubsurface%SetupSolvers()')
#endif

  option => this%option
  
  if (associated(this%timestepper)) then
    select type(ts => this%timestepper)
      class is(timestepper_BE_type)
        solver => ts%solver
        call SolverCreateSNES(solver,option%mycomm)
        call SNESGetLineSearch(ts%solver%snes,linesearch, &
                               ierr);CHKERRQ(ierr)
        ! set solver pointer within pm for convergence purposes
        call this%pm_ptr%pm%SetSolver(solver)
        select type(pm => this%pm_ptr%pm)
  ! ----- subsurface flow
          class is(pm_subsurface_flow_type)
            call printMsg(option,"  Beginning setup of FLOW SNES ")
            if (solver%J_mat_type == MATAIJ .and. &
                option%iflowmode /= TH_MODE) then
              option%io_buffer = 'AIJ matrix not supported for current &
                &mode: '// option%flowmode
              call printErrMsg(option)
            endif
            if (OptionPrintToScreen(option)) then
              write(*,'(" number of dofs = ",i3,", number of &
                        &phases = ",i3,i2)') option%nflowdof,option%nphase
              select case(option%iflowmode)
                case(TH_MODE)
                  write(*,'(" mode = TH: p, T")')
              end select
            endif

            call SNESSetOptionsPrefix(solver%snes, "flow_",ierr);CHKERRQ(ierr)
            call SolverCheckCommandLine(solver)

            if (solver%Jpre_mat_type == '') then
              if (solver%J_mat_type /= MATMFFD) then
                solver%Jpre_mat_type = solver%J_mat_type
              else
                solver%Jpre_mat_type = MATBAIJ
              endif
            endif

           if (associated(solver%cprstash)) then
              call CPRWorkersCreate(pm, solver, option)
            endif

            call DiscretizationCreateJacobian(pm%realization%discretization, &
                                              NFLOWDOF, &
                                              solver%Jpre_mat_type, &
                                              solver%Jpre, &
                                              option)

            call MatSetOptionsPrefix(solver%Jpre,"flow_",ierr);CHKERRQ(ierr)

            if (solver%J_mat_type /= MATMFFD) then
              solver%J = solver%Jpre
            endif

            if (solver%use_galerkin_mg) then
              call DiscretizationCreateInterpolation( &
                             pm%realization%discretization,NFLOWDOF, &
                             solver%interpolation, &
                             solver%galerkin_mg_levels_x, &
                             solver%galerkin_mg_levels_y, &
                             solver%galerkin_mg_levels_z, &
                             option)
            endif
    
            if (solver%J_mat_type == MATMFFD) then
              call MatCreateSNESMF(solver%snes,solver%J,ierr);CHKERRQ(ierr)
            endif

            ! by default turn off line search
            call SNESGetLineSearch(solver%snes, linesearch, ierr);CHKERRQ(ierr)
            call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC,  &
                                        ierr);CHKERRQ(ierr)
            ! Have PETSc do a SNES_View() at the end of each solve if 
            ! verbosity > 0.
            if (option%verbosity >= 2) then
              string = '-flow_snes_view'
              call PetscOptionsInsertString(PETSC_NULL_OPTIONS,string, &
                                            ierr);CHKERRQ(ierr)
            endif

            ! If we are using a structured grid, set the corresponding flow 
            ! DA as the DA for the PCEXOTIC preconditioner, in case we 
            ! choose to use it. The PCSetDA() call is ignored if the 
            ! PCEXOTIC preconditioner is no used.  We need to put this call 
            ! after SolverCreateSNES() so that KSPSetFromOptions() will 
            ! already have been called.  I also note that this 
            ! preconditioner is intended only for the flow 
            ! solver.  --RTM
            if (pm%realization%discretization%itype == STRUCTURED_GRID) then
              call PCSetDM(solver%pc, &
                           pm%realization%discretization%dm_nflowdof%dm, &
                           ierr);CHKERRQ(ierr)
            endif

            call SNESSetConvergenceTest(solver%snes, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                                        PMCheckConvergence, &
                                        this%pm_ptr%pm, &
#else
                                        PMCheckConvergencePtr, &
                                        this%pm_ptr, &
#endif
                                        PETSC_NULL_FUNCTION,ierr);CHKERRQ(ierr)
            if (pm%check_post_convergence) then
              call SNESLineSearchSetPostCheck(linesearch, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                                              PMCheckUpdatePost, &
                                              this%pm_ptr%pm, &
#else
                                              PMCheckUpdatePostPtr, &
                                              this%pm_ptr, &
#endif
                                              ierr);CHKERRQ(ierr)
              !geh: it is possible that the other side has not been set
              pm%check_post_convergence = PETSC_TRUE
            endif
            select type(pm)
            !-------------------------------------
              !TODO(geh): Can't these be consolidated?
            !-------------------------------------
              class is(pm_th_type)
                if (Initialized(pm%pressure_dampening_factor) .or. &
                    Initialized(pm%pressure_change_limit) .or. &
                    Initialized(pm%temperature_change_limit)) then
                  call SNESLineSearchSetPreCheck(linesearch, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                                                 PMCheckUpdatePre, &
                                                 this%pm_ptr%pm, &
#else
                                                 PMCheckUpdatePrePtr, &
                                                 this%pm_ptr, &
#endif
                                                 ierr);CHKERRQ(ierr)
                endif
            end select
            call printMsg(option,"  Finished setting up FLOW SNES ")
  ! ----- subsurface reactive transport                
          class is(pm_rt_type)
            call printMsg(option,"  Beginning setup of TRAN SNES ")
            call SNESSetOptionsPrefix(solver%snes, "tran_",ierr);CHKERRQ(ierr)
            call SolverCheckCommandLine(solver)
    
            if (option%transport%reactive_transport_coupling == &
                GLOBAL_IMPLICIT) then
              if (solver%Jpre_mat_type == '') then
                if (solver%J_mat_type /= MATMFFD) then
                  solver%Jpre_mat_type = solver%J_mat_type
                else
                  solver%Jpre_mat_type = MATBAIJ
                endif
              endif
              call DiscretizationCreateJacobian(pm%realization%discretization, &
                                                NTRANDOF, &
                                                solver%Jpre_mat_type, &
                                                solver%Jpre,option)
            else
              solver%J_mat_type = MATAIJ
              solver%Jpre_mat_type = MATAIJ

              call DiscretizationCreateJacobian(pm%realization%discretization, &
                                                ONEDOF, &
                                                solver%Jpre_mat_type, &
                                                solver%Jpre,option)
            endif

            if (solver%J_mat_type /= MATMFFD) then
              solver%J = solver%Jpre
            endif
    
            call MatSetOptionsPrefix(solver%Jpre,"tran_",ierr);CHKERRQ(ierr)
    
            if (solver%use_galerkin_mg) then
              call DiscretizationCreateInterpolation( &
                             pm%realization%discretization,NTRANDOF, &
                             solver%interpolation, &
                             solver%galerkin_mg_levels_x, &
                             solver%galerkin_mg_levels_y, &
                             solver%galerkin_mg_levels_z, &
                             option)
            endif

            if (option%transport%reactive_transport_coupling == &
                GLOBAL_IMPLICIT) then

              if (solver%J_mat_type == MATMFFD) then
                call MatCreateSNESMF(solver%snes,solver%J, &
                                      ierr);CHKERRQ(ierr)
              endif
      
              ! this could be changed in the future if there is a way to 
              ! ensure that the linesearch update does not perturb 
              ! concentrations negative.
              call SNESGetLineSearch(solver%snes, linesearch, &
                                     ierr);CHKERRQ(ierr)
              call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC,  &
                                          ierr);CHKERRQ(ierr)
      
              if (option%use_mc) then
                call SNESLineSearchSetPostCheck(linesearch, &
                                            SecondaryRTUpdateIterate, &
                                            pm%realization,ierr);CHKERRQ(ierr)
              endif
      
              ! Have PETSc do a SNES_View() at the end of each solve if 
              ! verbosity > 0.
              if (option%verbosity >= 2) then
                string = '-tran_snes_view'
                call PetscOptionsInsertString(PETSC_NULL_OPTIONS, &
                                              string, ierr);CHKERRQ(ierr)
              endif

            endif

            if (option%transport%reactive_transport_coupling == &
                GLOBAL_IMPLICIT) then
              call SNESSetConvergenceTest(solver%snes, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                                        PMCheckConvergence, &
                                        this%pm_ptr%pm, &
#else
                                        PMCheckConvergencePtr, &
                                        this%pm_ptr, &
#endif
                                        PETSC_NULL_FUNCTION,ierr);CHKERRQ(ierr)
            endif
            if (pm%print_EKG .or. option%use_mc .or. &
                pm%check_post_convergence) then
              call SNESLineSearchSetPostCheck(linesearch, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                                              PMCheckUpdatePost, &
                                              this%pm_ptr%pm, &
#else
                                              PMCheckUpdatePostPtr, &
                                              this%pm_ptr, &
#endif
                                              ierr);CHKERRQ(ierr)
              if (pm%print_EKG) then
                pm%check_post_convergence = PETSC_TRUE
              endif
            endif
            if (pm%realization%reaction%check_update) then
              call SNESLineSearchSetPreCheck(linesearch, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                                             PMCheckUpdatePre, &
                                             this%pm_ptr%pm, &
#else
                                             PMCheckUpdatePrePtr, &
                                             this%pm_ptr, &
#endif
                                             ierr);CHKERRQ(ierr)
            endif          
            call printMsg(option,"  Finished setting up TRAN SNES ")        
        end select
        call SNESSetFunction(ts%solver%snes, &
                             this%pm_ptr%pm%residual_vec, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                             PMResidual, &
                             this%pm_ptr%pm, &
#else
                             PMResidualPtr, &
                             this%pm_ptr, &
#endif
                             ierr);CHKERRQ(ierr)
        call SNESSetJacobian(ts%solver%snes, &
                             solver%J, &
                             solver%Jpre, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                             PMJacobian, &
                             this%pm_ptr%pm, &
#else
                             PMJacobianPtr, &
                             this%pm_ptr, &
#endif
                             ierr);CHKERRQ(ierr)
        call SolverSetSNESOptions(solver,option)
        option%io_buffer = 'Solver: ' // trim(solver%ksp_type)
        call printMsg(option)
        option%io_buffer = 'Preconditioner: ' // trim(solver%pc_type)
        call printMsg(option)
    end select
  endif ! associated(pmc%timestepper)        

end subroutine PMCSubsurfaceSetupSolvers

! ************************************************************************** !

subroutine PMCSubsurfaceGetAuxData(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 10/24/13
  ! 

  implicit none

  class(pmc_subsurface_type) :: this

  if (this%option%surf_flow_on) call PMCSubsurfaceGetAuxDataFromSurf(this)

end subroutine PMCSubsurfaceGetAuxData

! ************************************************************************** !

subroutine PMCSubsurfaceSetAuxData(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 10/24/13
  ! 

  implicit none

  class(pmc_subsurface_type) :: this

  if (this%option%surf_flow_on) call PMCSubsurfaceSetAuxDataForSurf(this)

end subroutine PMCSubsurfaceSetAuxData

! ************************************************************************** !

subroutine PMCSubsurfaceGetAuxDataFromSurf(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/22/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Connection_module
  use Coupler_module
  use Field_module
  use Grid_module
  use Option_module
  use Patch_module
!  use Realization_Base_class
  use Realization_Subsurface_class
  use String_module
  use EOS_Water_module

  implicit none

  class(pmc_subsurface_type) :: this
  
  class(realization_subsurface_type), pointer :: realization
  type (patch_type),pointer :: patch
  type (grid_type),pointer :: grid
  type (coupler_list_type), pointer :: coupler_list
  type (coupler_type), pointer :: coupler
  type (option_type), pointer :: option
  type (field_type),pointer :: field
  type (connection_set_type), pointer :: cur_connection_set
  PetscBool :: coupler_found
  PetscInt :: iconn
  PetscReal :: den
  PetscReal :: dt
  PetscReal :: surfpress
  PetscReal :: dum1
  PetscReal, pointer :: mflux_p(:)
  PetscReal, pointer :: hflux_p(:)
  PetscReal, pointer :: head_p(:)
  PetscReal, pointer :: temp_p(:)
  PetscErrorCode :: ierr

#ifdef DEBUG
  print *, 'PMCSubsurfaceGetAuxData()'
#endif

  dt = this%option%surf_subsurf_coupling_flow_dt

  if (associated(this%sim_aux)) then

    select type (pmc => this)
      class is (pmc_subsurface_type)

      if (this%sim_aux%subsurf_mflux_exchange_with_surf /= PETSC_NULL_VEC) then
        ! PETSc Vector to store relevant mass-flux data between
        ! surface-subsurface model exists

        patch      => pmc%realization%patch
        grid       => pmc%realization%discretization%grid
        field      => pmc%realization%field
        option     => pmc%realization%option

        select case(this%option%iflowmode)

          case (TH_MODE)
            call VecScatterBegin(pmc%sim_aux%surf_to_subsurf, &
                                 pmc%sim_aux%surf_head, &
                                 pmc%sim_aux%subsurf_pres_top_bc, &
                                 INSERT_VALUES,SCATTER_FORWARD, &
                                 ierr);CHKERRQ(ierr)
            call VecScatterEnd(pmc%sim_aux%surf_to_subsurf, &
                               pmc%sim_aux%surf_head, &
                               pmc%sim_aux%subsurf_pres_top_bc, &
                               INSERT_VALUES,SCATTER_FORWARD, &
                               ierr);CHKERRQ(ierr)

            call VecScatterBegin(pmc%sim_aux%surf_to_subsurf, &
                                 pmc%sim_aux%surf_temp, &
                                 pmc%sim_aux%subsurf_temp_top_bc, &
                                 INSERT_VALUES,SCATTER_FORWARD, &
                                 ierr);CHKERRQ(ierr)
            call VecScatterEnd(pmc%sim_aux%surf_to_subsurf, &
                               pmc%sim_aux%surf_temp, &
                               pmc%sim_aux%subsurf_temp_top_bc, &
                               INSERT_VALUES,SCATTER_FORWARD, &
                               ierr);CHKERRQ(ierr)

            call VecScatterBegin(pmc%sim_aux%surf_to_subsurf, &
                                 pmc%sim_aux%surf_hflux_exchange_with_subsurf, &
                                 pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                                 INSERT_VALUES,SCATTER_FORWARD, &
                                 ierr);CHKERRQ(ierr)
            call VecScatterEnd(pmc%sim_aux%surf_to_subsurf, &
                               pmc%sim_aux%surf_hflux_exchange_with_subsurf, &
                               pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                               INSERT_VALUES,SCATTER_FORWARD, &
                               ierr);CHKERRQ(ierr)

            coupler_list => patch%boundary_condition_list
            coupler => coupler_list%first
            do
              if (.not.associated(coupler)) exit

              ! FLOW
              if (associated(coupler%flow_aux_real_var)) then
                ! Find the BC from the list of BCs
                if (StringCompare(coupler%name,'from_surface_bc')) then
                  coupler_found = PETSC_TRUE

                  call VecGetArrayF90(pmc%sim_aux%subsurf_pres_top_bc, &
                                      head_p,ierr);CHKERRQ(ierr)
                  call VecGetArrayF90(pmc%sim_aux%subsurf_temp_top_bc, &
                                      temp_p,ierr);CHKERRQ(ierr)

                  do iconn = 1,coupler%connection_set%num_connections

                    ! The pressure value needed to computed density should
                    ! be surf_press and not reference_pressure. But,
                    ! surf_pressure depends on density.
                    call EOSWaterdensity(temp_p(iconn), option%reference_pressure, &
                                         den,dum1,ierr)

                    surfpress = head_p(iconn)*(abs(option%gravity(3)))*den + &
                                option%reference_pressure
                    coupler%flow_aux_real_var(TH_PRESSURE_DOF,iconn) = &
                      surfpress
                    coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                      temp_p(iconn)
                  enddo

                  call VecRestoreArrayF90(pmc%sim_aux%subsurf_pres_top_bc, &
                                          head_p,ierr);CHKERRQ(ierr)
                  call VecRestoreArrayF90(pmc%sim_aux%subsurf_temp_top_bc, &
                                      temp_p,ierr);CHKERRQ(ierr)
                endif
              endif

              if (StringCompare(coupler%name,'from_atm_subsurface_bc')) then
                coupler_found = PETSC_TRUE

                call VecGetArrayF90(pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                                    mflux_p,ierr);CHKERRQ(ierr)

                do iconn = 1,coupler%connection_set%num_connections
                  coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                    mflux_p(iconn)
                enddo

                call VecRestoreArrayF90(pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                                    mflux_p,ierr);CHKERRQ(ierr)
              endif

              coupler => coupler%next
            enddo

          case default
            this%option%io_buffer='PMCSubsurfaceGetAuxData() not supported for this mode.'
            call printErrMsg(this%option)

        end select

        if ( .not. coupler_found) then
          option%io_buffer = 'Coupler not found in PMCSubsurfaceGetAuxData()'
          call printErrMsg(option)
        endif
      endif

    end select

  endif ! if (associated(this%sim_aux))

end subroutine PMCSubsurfaceGetAuxDataFromSurf

! ************************************************************************** !

subroutine PMCSubsurfaceSetAuxDataForSurf(this)
  ! 
  ! This routine sets auxiliary to be exchanged between process-models.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/21/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Grid_module
  use String_module
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Coupler_module
  use Field_module
  use Connection_module
  use Realization_Base_class
  use EOS_Water_module

  implicit none

  class(pmc_subsurface_type) :: this
  
  class(realization_subsurface_type), pointer :: realization
  type (patch_type),pointer :: patch
  type (grid_type),pointer :: grid
  type (coupler_list_type), pointer :: coupler_list
  type (coupler_type), pointer :: coupler
  type (option_type), pointer :: option
  type (field_type),pointer :: field
  type (connection_set_type), pointer :: cur_connection_set
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iconn
  PetscInt :: istart
  PetscInt :: iend
  PetscReal :: den
  PetscReal :: dum1
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: pres_top_bc_p(:)
  PetscReal, pointer :: temp_top_bc_p(:)
  PetscReal, pointer :: head_p(:)
  PetscErrorCode :: ierr

#ifdef DEBUG
  print *, 'PMCSubsurfaceSetAuxData()'
#endif

  if (associated(this%sim_aux)) then

    select type (pmc => this)
      class is (pmc_subsurface_type)

        if (this%sim_aux%subsurf_pres_top_bc/= PETSC_NULL_VEC) then
          ! PETSc Vector to store relevant subsurface-flow data for
          ! surface-flow model exists

          patch      => pmc%realization%patch
          grid       => pmc%realization%discretization%grid
          field      => pmc%realization%field
          option     => pmc%realization%option

          call EOSWaterdensity(option%reference_temperature, option%reference_pressure, &
                               den,dum1,ierr)
          coupler_list => patch%boundary_condition_list
          coupler => coupler_list%first
          do
            if (.not.associated(coupler)) exit

            ! FLOW
            if (associated(coupler%flow_aux_real_var)) then

              ! Find the BC from the list of BCs
              if (StringCompare(coupler%name,'from_surface_bc')) then
                select case(this%option%iflowmode)
                  case (TH_MODE)
                    call VecGetArrayF90(this%sim_aux%subsurf_pres_top_bc, &
                                        pres_top_bc_p,ierr);CHKERRQ(ierr)
                    call VecGetArrayF90(this%sim_aux%subsurf_temp_top_bc, &
                                        temp_top_bc_p,ierr);CHKERRQ(ierr)

                    do iconn = 1,coupler%connection_set%num_connections
                      pres_top_bc_p(iconn) = &
                        coupler%flow_aux_real_var(TH_PRESSURE_DOF,iconn)
                      temp_top_bc_p(iconn) = &
                        coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn)
                    enddo

                    call VecRestoreArrayF90(this%sim_aux%subsurf_pres_top_bc, &
                                            pres_top_bc_p,ierr);CHKERRQ(ierr)
                    call VecRestoreArrayF90(this%sim_aux%subsurf_temp_top_bc, &
                                            temp_top_bc_p,ierr);CHKERRQ(ierr)
                    case default
                      option%io_buffer = 'PMCSubsurfaceGetAuxData() not ' // &
                        'supported in this FLOW_MODE'
                      call printErrMsg(option)
                end select
              endif
            endif

            coupler => coupler%next
          enddo

        endif
    end select

  endif

end subroutine PMCSubsurfaceSetAuxDataForSurf

! ************************************************************************** !
!
! PMCSubsurfaceFinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMCSubsurfaceFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Option_module
  
  implicit none
  
  class(pmc_subsurface_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'PMCSubsurface%FinalizeRun()')
#endif
  
  nullify(this%realization)
  
end subroutine PMCSubsurfaceFinalizeRun

! ************************************************************************** !

subroutine CPRWorkersCreate(pm, solver, option)
  ! 
  ! create all the worker/storage matrices/vectors that will be needed for the
  ! cpr preconditioner
  !
  ! Author: Daniel Stone
  ! Date: Oct 2017 - March 2018
  ! 

  use PM_Subsurface_Flow_class
  use Solver_module
  use Option_module
  use Discretization_module
  
  implicit none

  class(pm_subsurface_flow_type) :: pm
  class(solver_type) :: solver
  class(option_type) :: option
  MatType :: cpr_ap_mat_type



  cpr_ap_mat_type =  MATAIJ
  call DiscretizationCreateJacobian(pm%realization%discretization, &
                                    ONEDOF, &
                                    cpr_ap_mat_type, &
                                    solver%cprstash%Ap, &
                                    option)

  call DiscretizationCreateVector(pm%realization%discretization, &
                                  NFLOWDOF, solver%cprstash%T1r, &
                                  GLOBAL, option)
  call DiscretizationCreateVector(pm%realization%discretization, &
                                  NFLOWDOF, solver%cprstash%r2, &
                                  GLOBAL, option)

  call DiscretizationCreateVector(pm%realization%discretization, &
                                  ONEDOF, solver%cprstash%s,  &
                                  GLOBAL, option)
  call DiscretizationCreateVector(pm%realization%discretization, &
                                  ONEDOF, solver%cprstash%z, &
                                  GLOBAL, option)

  call DiscretizationCreateVector(pm%realization%discretization, &
                                  NFLOWDOF, solver%cprstash%factors1vec, &
                                  GLOBAL, option)
  call DiscretizationCreateVector(pm%realization%discretization, &
                                  NFLOWDOF, solver%cprstash%factors2vec, &
                                  GLOBAL, option)
end subroutine CPRWorkersCreate

! ************************************************************************** !

subroutine PMCSubsurfaceStrip(this)
  !
  ! Deallocates members of PMC Subsurface.
  !
  ! Author: Glenn Hammond
  ! Date: 01/13/14
  
  implicit none
  
  class(pmc_subsurface_type) :: this

  call PMCBaseStrip(this)
  nullify(this%realization)

end subroutine PMCSubsurfaceStrip

! ************************************************************************** !

recursive subroutine PMCSubsurfaceDestroy(this)
  ! 
  ! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Option_module

  implicit none
  
  class(pmc_subsurface_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'PMCSubsurface%Destroy()')
#endif

  if (associated(this%child)) then
    call this%child%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%child)
    nullify(this%child)
  endif 
  
  if (associated(this%peer)) then
    call this%peer%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%peer)
    nullify(this%peer)
  endif
  
  call PMCSubsurfaceStrip(this)
  
end subroutine PMCSubsurfaceDestroy

! ************************************************************************** !
  
end module PMC_Subsurface_class
