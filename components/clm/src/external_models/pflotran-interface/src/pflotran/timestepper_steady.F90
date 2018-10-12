module Timestepper_Steady_class

  use Timestepper_BE_class
  use Timestepper_Base_class
  use Convergence_module
  use Solver_module
  use Waypoint_module
  use PFLOTRAN_Constants_module

  implicit none

#include "petsc/finclude/petscsys.h"
 
  type, public, extends(timestepper_BE_type) :: timestepper_steady_type

  contains

 !   procedure, public :: Init => TimestepperSteadyInit
    procedure, public :: StepDT => TimestepperSteadyStepDT
    procedure, public :: UpdateDT => TimestepperSteadyUpdateDT
    procedure, public :: InputRecord => TimestepperSteadyInputRecord

  end type timestepper_steady_type

  public :: TimestepperSteadyCreate, &
            TimestepperSteadyCreateFromBE

contains

! ************************************************************************** !

function TimestepperSteadyCreate()
  ! 
  ! This routine creates timestepper for steady solve
  ! 
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 04/07/2015
  ! 

  implicit none

  class(timestepper_steady_type), pointer :: TimestepperSteadyCreate

  class(timestepper_steady_type), pointer :: stepper

  allocate(stepper)
  call stepper%Init()

  stepper%solver => SolverCreate()

  TimestepperSteadyCreate => stepper

end function TimestepperSteadyCreate

! ************************************************************************** !

subroutine TimestepperSteadyCreateFromBE(timestepper_BE)
  ! 
  ! This routine creates timestepper for steady solve
  ! 
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 04/07/2015
  ! 

  implicit none

  class(timestepper_BE_type), pointer :: timestepper_BE

  class(timestepper_steady_type), pointer :: stepper

  allocate(stepper)
  
  stepper%name = timestepper_BE%name
  stepper%steps = timestepper_BE%steps
  stepper%num_constant_time_steps = timestepper_BE%num_constant_time_steps

  stepper%max_time_step = timestepper_BE%max_time_step
  stepper%max_time_step_cuts = timestepper_BE%max_time_step_cuts
  stepper%constant_time_step_threshold = timestepper_BE%constant_time_step_threshold

  stepper%cumulative_time_step_cuts = timestepper_BE%cumulative_time_step_cuts    
  stepper%cumulative_solver_time = timestepper_BE%cumulative_solver_time

  stepper%start_time = timestepper_BE%start_time
  stepper%start_time_step = timestepper_BE%start_time_step
  stepper%time_step_tolerance = timestepper_BE%time_step_tolerance
  stepper%target_time = timestepper_BE%target_time
  
  stepper%prev_dt = timestepper_BE%prev_dt
  stepper%dt = timestepper_BE%dt
  stepper%dt_init = timestepper_BE%dt_init
  stepper%dt_max = timestepper_BE%dt_max
  
  stepper%time_step_cut_flag = timestepper_BE%time_step_cut_flag
  
  stepper%init_to_steady_state = timestepper_BE%init_to_steady_state
  stepper%steady_state_rel_tol = timestepper_BE%steady_state_rel_tol
  stepper%run_as_steady_state = timestepper_BE%run_as_steady_state
  
  stepper%cur_waypoint => timestepper_BE%cur_waypoint
  stepper%prev_waypoint => timestepper_BE%prev_waypoint
  stepper%revert_dt = timestepper_BE%revert_dt
  stepper%num_contig_revert_due_to_sync = timestepper_BE%num_contig_revert_due_to_sync
      
  stepper%num_newton_iterations = timestepper_BE%num_newton_iterations
  stepper%num_linear_iterations = timestepper_BE%num_linear_iterations

  stepper%cumulative_newton_iterations = timestepper_BE%cumulative_newton_iterations
  stepper%cumulative_linear_iterations = timestepper_BE%cumulative_linear_iterations

  stepper%iaccel = timestepper_BE%iaccel
  stepper%ntfac = timestepper_BE%ntfac
  allocate(stepper%tfac(13))
  stepper%tfac(1) = timestepper_BE%tfac(1)  
  stepper%tfac(2) = timestepper_BE%tfac(2)  
  stepper%tfac(3) = timestepper_BE%tfac(3)  
  stepper%tfac(4) = timestepper_BE%tfac(4)  
  stepper%tfac(5) = timestepper_BE%tfac(5)  
  stepper%tfac(6) = timestepper_BE%tfac(6)  
  stepper%tfac(7) = timestepper_BE%tfac(7)  
  stepper%tfac(8) = timestepper_BE%tfac(8)  
  stepper%tfac(9) = timestepper_BE%tfac(9)  
  stepper%tfac(10) = timestepper_BE%tfac(10)  
  stepper%tfac(11) = timestepper_BE%tfac(11)  
  stepper%tfac(12) = timestepper_BE%tfac(12)  
  stepper%tfac(13) = timestepper_BE%tfac(13)  
  
  stepper%solver => timestepper_BE%solver

  timestepper_BE => stepper

end subroutine TimestepperSteadyCreateFromBE

! ************************************************************************** !

subroutine TimestepperSteadyInit(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 04/07/2015
  ! 

  implicit none
  
  class(timestepper_steady_type) :: this

  call TimestepperBaseInit(this)

end subroutine TimestepperSteadyInit


! ************************************************************************** !

subroutine TimestepperSteadyUpdateDT(this,process_model)
  ! 
  ! Updates time step
  ! 
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 04/07/2015
  ! 

  use PM_Base_class
  use Option_module
  
  implicit none

  class(timestepper_steady_type) :: this
  class(pm_base_type) :: process_model
  

end subroutine TimestepperSteadyUpdateDT

! ************************************************************************** !

subroutine TimestepperSteadyStepDT(this, process_model, stop_flag)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 04/07/2015
  ! 

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use Option_module
  use Output_module, only : Output

  use Solver_module

  implicit none

  class(timestepper_steady_type) :: this
  class(pm_base_type) :: process_model
  PetscInt :: stop_flag

  PetscBool :: failure
  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_end_time
  PetscErrorCode :: ierr
  PetscInt :: sum_newton_iterations, sum_linear_iterations
  PetscInt :: num_newton_iterations, num_linear_iterations
  PetscInt :: snes_reason
  PetscInt :: icut
  PetscReal :: fnorm
  PetscReal :: inorm
  PetscReal :: scaled_fnorm
  Vec :: residual_vec

  type(option_type), pointer :: option
  type(solver_type), pointer :: solver

  solver => this%solver
  option => process_model%option

  sum_newton_iterations = 0
  sum_linear_iterations = 0
  icut = 0

  call process_model%InitializeTimestep()

  call process_model%PreSolve()

  call PetscTime(log_start_time, ierr);CHKERRQ(ierr)

  call SNESSolve(solver%snes, PETSC_NULL_VEC, &
                 process_model%solution_vec, ierr);CHKERRQ(ierr)

  call PetscTime(log_end_time, ierr);CHKERRQ(ierr)

  this%cumulative_solver_time = &
      this%cumulative_solver_time + &
      (log_end_time - log_start_time)
     
  call SNESGetIterationNumber(solver%snes, num_newton_iterations,  &
                              ierr);CHKERRQ(ierr)
  call SNESGetLinearSolveIterations(solver%snes, num_linear_iterations,  &
                                    ierr);CHKERRQ(ierr)
  call SNESGetConvergedReason(solver%snes, snes_reason, ierr);CHKERRQ(ierr)

  sum_newton_iterations = sum_newton_iterations + num_newton_iterations
  sum_linear_iterations = sum_linear_iterations + num_linear_iterations

  if (snes_reason <= 0) then
    if (option%print_screen_flag) then
      print *, 'Newton solver failed to converge in steady-solve, reason: ', &
                snes_reason
    endif
    failure = PETSC_TRUE
    stop_flag = TS_STOP_END_SIMULATION
    return
  endif
  
  this%steps = this%steps + 1
  this%cumulative_newton_iterations = &
    this%cumulative_newton_iterations + sum_newton_iterations
  this%cumulative_linear_iterations = &
    this%cumulative_linear_iterations + sum_linear_iterations
  this%cumulative_time_step_cuts = &
    this%cumulative_time_step_cuts + icut

  this%num_newton_iterations = num_newton_iterations
  this%num_linear_iterations = num_linear_iterations  

  ! print screen output
  call SNESGetFunction(solver%snes,residual_vec,PETSC_NULL_FUNCTION, &
                       PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call VecNorm(residual_vec,NORM_2,fnorm,ierr);CHKERRQ(ierr)
  call VecNorm(residual_vec,NORM_INFINITY,inorm,ierr);CHKERRQ(ierr)
  if (option%print_screen_flag) then
    select type(pm => process_model)
    end select
    write(*,*) ''
    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'(" --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif
  
  if (option%print_screen_flag) print *, ""
  
  if (option%print_file_flag) then
    write(option%fid_out, '(" STEADY-SOLVE ",i6," snes_conv_reason: ",i4,/, &
      &"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]")') &
      this%steps, &
      snes_reason,sum_newton_iterations, &
      this%cumulative_newton_iterations,sum_linear_iterations, &
      this%cumulative_linear_iterations
  endif  

  option%time = this%target_time
  call process_model%FinalizeTimestep()
  
  if (option%print_screen_flag) print *, ""
  ! check if the steady state option is selected in input deck
  ! if yes then end simulation
  if (option%steady_state) stop_flag = TS_STOP_END_SIMULATION

end subroutine TimestepperSteadyStepDT

! ************************************************************************** !

subroutine TimestepperSteadyInputRecord(this)
  ! 
  ! Prints information about the time stepper to the input record.
  ! To get a## format, must match that in simulation types.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  ! 
  
  implicit none
  
  class(timestepper_steady_type) :: this

  PetscInt :: id
  character(len=MAXWORDLENGTH) :: word
   
  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pmc timestepper: '
  write(id,'(a)') this%name

end subroutine TimestepperSteadyInputRecord

end module Timestepper_Steady_class
