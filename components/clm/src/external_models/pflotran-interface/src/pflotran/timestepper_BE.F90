module Timestepper_BE_class
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Solver_module
  use Convergence_module
  use Timestepper_Base_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "petsc/finclude/petscsys.h"
 
  type, public, extends(timestepper_base_type) :: timestepper_BE_type
  
    PetscInt :: num_newton_iterations ! number of Newton iterations in a time step
    PetscInt :: num_linear_iterations ! number of linear solver iterations in a time step
    PetscInt :: cumulative_newton_iterations       ! Total number of Newton iterations
    PetscInt :: cumulative_linear_iterations     ! Total number of linear iterations
    PetscInt :: cumulative_wasted_linear_iterations

    PetscInt :: iaccel        ! Accelerator index
    ! An array of multiplicative factors that specify how to increase time step.
    PetscReal, pointer :: tfac(:)
    PetscInt :: ntfac             ! size of tfac
            
    type(solver_type), pointer :: solver
  
  contains
    
    procedure, public :: ReadInput => TimestepperBERead
    procedure, public :: Init => TimestepperBEInit
!    procedure, public :: SetTargetTime => TimestepperBaseSetTargetTime
    procedure, public :: StepDT => TimestepperBEStepDT
    procedure, public :: UpdateDT => TimestepperBEUpdateDT
    procedure, public :: CheckpointBinary => TimestepperBECheckpointBinary
    procedure, public :: RestartBinary => TimestepperBERestartBinary
#if defined(PETSC_HAVE_HDF5)
    procedure, public :: CheckpointHDF5 => TimestepperBECheckpointHDF5
    procedure, public :: RestartHDF5 => TimestepperBERestartHDF5
#endif
    procedure, public :: Reset => TimestepperBEReset
    procedure, public :: PrintInfo => TimestepperBEPrintInfo
    procedure, public :: InputRecord => TimestepperBEInputRecord
    procedure, public :: FinalizeRun => TimestepperBEFinalizeRun
    procedure, public :: Strip => TimestepperBEStrip
    procedure, public :: Destroy => TimestepperBEDestroy
    
  end type timestepper_BE_type
  
  ! For checkpointing
  type, public, extends(stepper_base_header_type) :: stepper_BE_header_type
    PetscInt :: cumulative_newton_iterations
    PetscInt :: cumulative_linear_iterations
    PetscInt :: num_newton_iterations
  end type stepper_BE_header_type

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
      import :: stepper_BE_header_type
      implicit none
#include "petsc/finclude/petscbag.h"      
      PetscBag :: bag
      class(stepper_BE_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData  

  public :: TimestepperBECreate, TimestepperBEPrintInfo, &
            TimestepperBEInit

contains

! ************************************************************************** !

function TimestepperBECreate()
  ! 
  ! Allocates and initializes a new Timestepper object
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  ! 

  implicit none
  
  class(timestepper_BE_type), pointer :: TimestepperBECreate
  
  class(timestepper_BE_type), pointer :: stepper
  
  allocate(stepper)
  call stepper%Init()
  
  stepper%solver => SolverCreate()
  
  TimestepperBECreate => stepper
  
end function TimestepperBECreate

! ************************************************************************** !

subroutine TimestepperBEInit(this)
  ! 
  ! Allocates and initializes a new Timestepper object
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  ! 

  implicit none
  
  class(timestepper_BE_type) :: this
  
  call TimestepperBaseInit(this)
  
  this%num_newton_iterations = 0
  this%num_linear_iterations = 0

  this%cumulative_newton_iterations = 0
  this%cumulative_linear_iterations = 0
  this%cumulative_wasted_linear_iterations = 0

  this%iaccel = 5
  this%ntfac = 13
  allocate(this%tfac(13))
  this%tfac(1)  = 2.0d0; this%tfac(2)  = 2.0d0
  this%tfac(3)  = 2.0d0; this%tfac(4)  = 2.0d0
  this%tfac(5)  = 2.0d0; this%tfac(6)  = 1.8d0
  this%tfac(7)  = 1.6d0; this%tfac(8)  = 1.4d0
  this%tfac(9)  = 1.2d0; this%tfac(10) = 1.0d0
  this%tfac(11) = 1.0d0; this%tfac(12) = 1.0d0
  this%tfac(13) = 1.0d0
  
  nullify(this%solver)
  
end subroutine TimestepperBEInit

! ************************************************************************** !

subroutine TimestepperBERead(this,input,option)
  ! 
  ! Reads parameters associated with time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none

  class(timestepper_BE_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: string

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','TIMESTEPPER_BE')
    call StringToUpper(keyword)   

    select case(trim(keyword))
  
      case('TS_ACCELERATION')
        call InputReadInt(input,option,this%iaccel)
        call InputDefaultMsg(input,option,'iaccel')

      case('DT_FACTOR')
        string='time_step_factor'
        call UtilityReadArray(this%tfac,NEG_ONE_INTEGER,string,input, &
            option)
        this%ntfac = size(this%tfac)

      case default
        call TimestepperBaseProcessKeyword(this,input,option,keyword)
    end select 
  
  enddo
  
  this%solver%print_ekg = this%print_ekg

end subroutine TimestepperBERead

! ************************************************************************** !

subroutine TimestepperBEUpdateDT(this,process_model)
  ! 
  ! Updates time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  ! 

  use PM_Base_class
  
  implicit none

  class(timestepper_BE_type) :: this
  class(pm_base_type) :: process_model
  
  PetscBool :: update_time_step
  
  update_time_step = PETSC_TRUE

  if (this%time_step_cut_flag) then
    this%num_constant_time_steps = 1
  else if (this%num_constant_time_steps > 0) then
    ! otherwise, only increment if the constant time step counter was
    ! initialized to 1
    this%num_constant_time_steps = &
      this%num_constant_time_steps + 1
  endif

  ! num_constant_time_steps = 0: normal time stepping with growing steps
  ! num_constant_time_steps > 0: restriction of constant time steps until
  !                              constant_time_step_threshold is met
  if (this%num_constant_time_steps > &
      this%constant_time_step_threshold) then
    this%num_constant_time_steps = 0
  else if (this%num_constant_time_steps > 0) then
    ! do not increase time step size
    update_time_step = PETSC_FALSE
  endif
    
  if (update_time_step .and. this%iaccel /= 0) then
      
    call process_model%UpdateTimestep(this%dt, &
                                      this%dt_min, &
                                      this%dt_max, &
                                      this%iaccel, &
                                      this%num_newton_iterations, &
                                      this%tfac, &
                                      this%time_step_max_growth_factor)
    
  endif

end subroutine TimestepperBEUpdateDT

! ************************************************************************** !

subroutine TimestepperBEStepDT(this,process_model,stop_flag)
  ! 
  ! Steps forward one step in time
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  ! 

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use Option_module
  use Output_module, only : Output, OutputFindNaNOrInfInVec
  use Output_EKG_module, only : IUNIT_EKG

  ! for checking grid cell index in which error occurs (F.-M. Yuan, 2017-02-24)
  use PM_Subsurface_Flow_class, only : pm_subsurface_flow_type
  use grid_module, only : grid_type
  use grid_structured_module, only : StructGridGetIJKFromLocalID
  
  implicit none

  class(timestepper_BE_type) :: this
  class(pm_base_type) :: process_model
  PetscInt :: stop_flag
  
  SNESConvergedReason :: snes_reason
  PetscInt :: icut
  
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  
  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_end_time
  PetscInt :: num_newton_iterations
  PetscInt :: num_linear_iterations
  PetscInt :: num_linear_iterations2
  PetscInt :: sum_newton_iterations
  PetscInt :: sum_linear_iterations
  PetscInt :: sum_wasted_linear_iterations
  character(len=MAXWORDLENGTH) :: tunit
  PetscReal :: tconv
  PetscReal :: fnorm, inorm, scaled_fnorm
  PetscBool :: snapshot_plot_flag, observation_plot_flag, massbal_plot_flag
  Vec :: residual_vec
  PetscErrorCode :: ierr

!fmy: for printing vecs if program stops
  PetscScalar, pointer :: solution_p(:)
  PetscScalar, pointer :: residual_p(:)
  PetscInt :: fileid_info
  PetscInt :: vecsize, i
  PetscInt :: cell_index, cell_i, cell_j, cell_k
  PetscInt :: cell_offset, cell_offset_x, cell_offset_y, cell_offset_z
  PetscReal:: max_res(2), max_res_global(2), max_res_value, max_res_solu
  PetscInt :: max_res_index, max_res_rank
  PetscErrorCode :: ierr2
  type(grid_type), pointer :: grid
!fmy: for printing vecs if program stops

  ! GNU -O3 can fail below in SNESGetFunction() as the compiler can set the
  ! initial value to -1, which CHKFORTRANNULLOBJECT() interprets as NULL.
  residual_vec = tVec(0)
  
  solver => this%solver
  option => process_model%option
  
!geh: for debugging
#ifdef DEBUG
  write(process_model%option%io_buffer,'(es12.5)') this%dt
  process_model%option%io_buffer = 'StepperStepDT(' // &
    trim(adjustl(process_model%option%io_buffer)) // ')'
  call printMsg(process_model%option)  
#endif

  tconv = process_model%output_option%tconv
  tunit = process_model%output_option%tunit
  sum_linear_iterations = 0
  sum_wasted_linear_iterations = 0
  sum_newton_iterations = 0
  icut = 0
  
  option%dt = this%dt
  option%time = this%target_time-this%dt

  call process_model%InitializeTimestep()
    
  do
      
    call process_model%PreSolve()
    
    call PetscTime(log_start_time, ierr);CHKERRQ(ierr)

    call SNESSolve(solver%snes,PETSC_NULL_VEC, &
                   process_model%solution_vec,ierr);

    if (ierr .ne. 0) then

      print *, ' <-- SNES Solver ERROR @TimeStepperBEStepDT --> '
      print *, ' Time (s): ', option%time, ' log_start_time: ', log_start_time
      print *, ' dT (s): ', option%dt
      print *, ' Linear Iterations: ', sum_linear_iterations
      print *, ' Newton Iterations: ', sum_newton_iterations
      print *, 'PETSC error id: ', ierr

      if (option%print_file_flag) then

        !------fmy - checking the solution to see what's going on
        call VecGetLocalSize(process_model%residual_vec, vecsize, ierr2)

        call VecGetArrayF90(process_model%residual_vec, residual_p, ierr2)
        call VecGetArrayF90(process_model%solution_vec, solution_p, ierr2)
        max_res_value= abs(residual_p(1))
        max_res_index= 1
        do i=2, vecsize
          if (max_res(1) < abs(residual_p(i))) then
            max_res_value = abs(residual_p(i))
            max_res_index = i                  ! local at this point
            max_res_solu  = solution_p(i)
          endif
        enddo

        ! a note here: because only operating on 'ierr .ne. 0' process, NO need to do MPI_Allreduce operation (TODO -checking)

        ! NOTE: the following will print info to an un-initialized txt file: fort.'max_res_rank+?'.
        fileid_info = option%myrank+1   ! when @ rank = 0, info goes to screen (not sure why), while +1 will work around

        write(fileid_info, *) ' <-- SNES Solver ERROR @TimeStepperBEStepDT -->'

        write(fileid_info, *) 'Time(s): ', option%time, 'rank: ', option%myrank, 'Petsc ErrCode: ', ierr
        ! not yet figure out how to get 'reaction_aux%ncomp' in
        ! and either the 2 vecs are different for flow-transport model and reaction model
        ! BE CAUTIOUS!
        if (option%nflowdof>0 .or. option%ntrandof>0) then
          cell_index = floor(real((max_res_index-1)/(option%nflowdof+option%ntrandof)))+1

          select type (pm => process_model)
            class is (pm_subsurface_flow_type)
              grid => pm%realization%patch%grid

              select case(grid%itype)
                case(STRUCTURED_GRID)
                  call StructGridGetIJKFromLocalID(grid%structured_grid, cell_index, &
                                                   cell_i, cell_j, cell_k)
                  cell_offset   = grid%structured_grid%global_offset
                  cell_offset_x = grid%structured_grid%lxs  ! global corner starting index (0-based)
                  cell_offset_y = grid%structured_grid%lys
                  cell_offset_z = grid%structured_grid%lzs
                  ! NOTE: cell_offset+cell_index ~ global cell index (1-based)
                  !       cell_offset_x/y/z + cell_i/j/k ~ global cell x/y/z index (0-based)
                  !  this info will help to locate which cell has very likely caused error due to its largest RES.

                case default
                  option%io_buffer = " only STRUCTURED_GRID can know cell_IJK."
                  call printMsg(option)
              end select
            !
          end select

          write(fileid_info, *) ' <--- cell offset ---- cell index --- vec. no. --elem. no. -- solution_vec with max. res ----> '
          do i=(cell_index-1)*(option%nflowdof+option%ntrandof)+1, &
                cell_index*(option%nflowdof+option%ntrandof)
            write(fileid_info, *) cell_offset,'(', cell_offset_x,cell_offset_y,cell_offset_z,')', &
                cell_index, '(',cell_i,cell_j,cell_k, ')', &
                i, i-(cell_index-1)*(option%nflowdof+option%ntrandof), solution_p(i)
          enddo

          write(fileid_info, *) '  '
          write(fileid_info, *) ' <--- cell offset ---- cell index --- vec no. --elem. no. -- max. residual_vec ----> '
          do i=(cell_index-1)*(option%nflowdof+option%ntrandof)+1, &
                cell_index*(option%nflowdof+option%ntrandof)
            write(fileid_info, *) cell_offset,'(', cell_offset_x,cell_offset_y,cell_offset_z,')', &
                cell_index, '(',cell_i,cell_j,cell_k, ')', &
                i, i-(cell_index-1)*(option%nflowdof+option%ntrandof), residual_p(i)
          enddo

        endif

        call VecRestoreArrayF90(process_model%residual_vec, residual_p, ierr2)
        call VecRestoreArrayF90(process_model%solution_vec, solution_p, ierr2)

        write(fileid_info, *) '  '
        write(fileid_info, *) ' Stop Executing! '

      endif

      print *, ' Stop Executing!'
      CHKERRQ(ierr)
    endif
!fmy: checking SNESSolver error and stop excuting/output messages if error occurs

    CHKERRQ(ierr)

    call PetscTime(log_end_time, ierr);CHKERRQ(ierr)

    this%cumulative_solver_time = &
      this%cumulative_solver_time + &
      (log_end_time - log_start_time)

    call SNESGetIterationNumber(solver%snes,num_newton_iterations, &
                                ierr);CHKERRQ(ierr)
    call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations, &
                                      ierr);CHKERRQ(ierr)
    call SNESGetConvergedReason(solver%snes,snes_reason,ierr);CHKERRQ(ierr)

    sum_newton_iterations = sum_newton_iterations + num_newton_iterations
    sum_linear_iterations = sum_linear_iterations + num_linear_iterations
  
    if (snes_reason <= 0 .or. .not. process_model%AcceptSolution()) then
      sum_wasted_linear_iterations = sum_wasted_linear_iterations + &
        num_linear_iterations
      ! The Newton solver diverged, so try reducing the time step.
      icut = icut + 1
      this%time_step_cut_flag = PETSC_TRUE
      ! if a cut occurs on the last time step, the stop_flag will have been
      ! set to TS_STOP_END_SIMULATION.  Set back to TS_CONTINUE to prevent
      ! premature ending of simulation.
      if (stop_flag /= TS_STOP_MAX_TIME_STEP) stop_flag = TS_CONTINUE

      if (icut > this%max_time_step_cuts .or. this%dt < this%dt_min) then

        !------fmy - checking the solution to see what's going on
        call VecGetLocalSize(process_model%residual_vec, vecsize, ierr2)

        call VecGetArrayF90(process_model%residual_vec, residual_p, ierr2)
        call VecGetArrayF90(process_model%solution_vec, solution_p, ierr2)
        max_res(1)   = abs(residual_p(1))
        max_res_index= 1
        do i=2, vecsize
          if (max_res(1) < abs(residual_p(i))) then
            max_res(1)    = abs(residual_p(i))
            max_res_index = i                  ! local at this point
            max_res_solu  = solution_p(i)
          endif
        enddo

        max_res(2)   = real(option%myrank)
        call MPI_Allreduce(max_res,max_res_global,ONE_INTEGER_MPI, &
                 MPI_2DOUBLE_PRECISION, MPI_MAXLOC, option%mycomm,ierr)
        max_res_value= max_res_global(1)
        max_res_rank = int(max_res_global(2))   ! global rank in which 'max_res(1)' located

        max_res(2)    = real(max_res_index)     ! local at this point
        call MPI_Allreduce(max_res,max_res_global,ONE_INTEGER_MPI, &
                 MPI_2DOUBLE_PRECISION, MPI_MAXLOC, option%mycomm,ierr)
        max_res_index = int(max_res_global(2))  ! global 'index' in which 'max_res(1)' located

        max_res(2)    = max_res_solu            ! local at this point
        call MPI_Allreduce(max_res,max_res_global,ONE_INTEGER_MPI, &
                 MPI_2DOUBLE_PRECISION, MPI_MAXLOC, option%mycomm,ierr)
        max_res_solu  = max_res_global(2)       ! global 'resolution' in which 'max_res(1)' located

        ! NOTE: the following will print info to an un-initialized txt file: fort.'max_res_rank+?'.
        fileid_info = option%myrank+1   ! when @ rank = 0, info goes to screen (not sure why), while +1 will work around
        if (option%myrank == max_res_rank) then

          write(fileid_info, *) ' <-- SNES Solver checking @TimeStepperBEStepDT -->, @ rank, @time', &
            max_res_rank, ' max_res: ', max_res_value, option%time
          if (icut > this%max_time_step_cuts) then
            write(fileid_info, *) ' -- MAX. CUTS reached --', icut
          elseif (this%dt < this%dt_min) then
            write(fileid_info, *) ' -- MIN. TIME-STEPS (sec.) reached --', this%dt
          endif


          ! not yet figure out how to get 'reaction_aux%ncomp' in
          ! and neither the 2 vecs are different for flow-transport model and reaction model
          ! BE CAUTIOUS!
          if (option%nflowdof>0 .or. option%ntrandof>0) then
            cell_index = floor(real((max_res_index-1))/(option%nflowdof+option%ntrandof))+1

            select type (pm => process_model)
              class is (pm_subsurface_flow_type)
                grid => pm%realization%patch%grid

                select case(grid%itype)
                  case(STRUCTURED_GRID)
                    call StructGridGetIJKFromLocalID(grid%structured_grid, cell_index, &
                                                   cell_i, cell_j, cell_k)
                    cell_offset   = grid%structured_grid%global_offset   ! total cells prior to this rank
                    cell_offset_x = grid%structured_grid%lxs
                    cell_offset_y = grid%structured_grid%lys
                    cell_offset_z = grid%structured_grid%lzs
                    ! NOTE: cell_offset+cell_index ~ global cell index (1-based)
                    !       cell_offset_x/y/z + cell_i/j/k ~ global cell x/y/z index (0-based)
                    !  this info will help to locate which cell has very likely caused error due to its largest RES.
                  case default
                    option%io_buffer = " only STRUCTURED_GRID can know cell_IJK."
                    call printMsg(option)
                end select
            end select

            write(fileid_info, *) ' <--- cell offset ---- cell index --- vec. no. --elem. no. -- solution_vec with max. res ----> '
            do i=(cell_index-1)*(option%nflowdof+option%ntrandof)+1, &
                 cell_index*(option%nflowdof+option%ntrandof)
              write(fileid_info, *) cell_offset,'(', cell_offset_x,cell_offset_y,cell_offset_z,')', &
                cell_index, '(',cell_i,cell_j,cell_k, ')', &
                i, i-(cell_index-1)*(option%nflowdof+option%ntrandof), solution_p(i)
            enddo

            write(fileid_info, *) '  '
            write(fileid_info, *) ' <--- cell offset ---- cell index --- vec no. --elem. no. -- max. residual_vec ----> '
            do i=(cell_index-1)*(option%nflowdof+option%ntrandof)+1, &
                  cell_index*(option%nflowdof+option%ntrandof)
              write(fileid_info, *) cell_offset,'(', cell_offset_x,cell_offset_y,cell_offset_z,')', &
                cell_index, '(',cell_i,cell_j,cell_k, ')', &
                i, i-(cell_index-1)*(option%nflowdof+option%ntrandof), residual_p(i)
            enddo
          endif

          call VecRestoreArrayF90(process_model%solution_vec, solution_p, ierr2)
          call VecRestoreArrayF90(process_model%residual_vec, residual_p, ierr2)

          write(fileid_info, *) '  '
          write(fileid_info, *) ' Stop Executing! @ myrank: ', option%myrank

        endif
        !------fmy

        write(option%io_buffer,'(" Stopping: Time step cut criteria exceeded!")')
        call printMsg(option)
        write(option%io_buffer,'("    icut =",i3,", max_time_step_cuts=",i3)') &
             icut,this%max_time_step_cuts
        call printMsg(option)
        write(option%io_buffer,'("    dt   =",es15.7,", dt_min=",es15.7)') &
             this%dt/tconv,this%dt_min/tconv
        call printMsg(option)
        
        process_model%output_option%plot_name = 'flow_cut_to_failure'
        snapshot_plot_flag = PETSC_TRUE
        observation_plot_flag = PETSC_FALSE
        massbal_plot_flag = PETSC_FALSE
        call Output(process_model%realization_base,snapshot_plot_flag, &
                    observation_plot_flag,massbal_plot_flag)
        stop_flag = TS_STOP_FAILURE
        return
      endif
 
      this%target_time = this%target_time - this%dt

      this%dt = this%time_step_reduction_factor * this%dt

#ifndef CLM_PFLOTRAN
      write(option%io_buffer,'('' -> Cut time step: snes='',i3, &
           &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.5, '' dt= '', &
           &   1pe12.5)')  snes_reason,icut,this%cumulative_time_step_cuts, &
           option%time/tconv, &
           this%dt/tconv
      if (option%print_screen_flag) &
      call printMsg(option)
#endif
      if (snes_reason < SNES_CONVERGED_ITERATING) then
        call SolverNewtonPrintFailedReason(solver,option)
        if (solver%verbose_error_msg) then
          select case(snes_reason)
            case(SNES_DIVERGED_FNORM_NAN)
              ! attempt to find cells with NaNs.
              call SNESGetFunction(solver%snes,residual_vec, &
                                   PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER, &
                                 ierr);CHKERRQ(ierr)
              call OutputFindNaNOrInfInVec(residual_vec, &
                                           process_model%realization_base% &
                                             discretization%grid,option)
          end select
        endif
        call KSPGetIterationNumber(solver%ksp,num_linear_iterations2, &
                                   ierr);CHKERRQ(ierr)
        sum_wasted_linear_iterations = sum_wasted_linear_iterations + &
          num_linear_iterations2
        sum_linear_iterations = sum_linear_iterations + num_linear_iterations2
      endif

      this%target_time = this%target_time + this%dt
      option%dt = this%dt
      call process_model%TimeCut()
  
    else
      ! The Newton solver converged, so we can exit.
      exit
    endif
  enddo

  this%steps = this%steps + 1      
  this%cumulative_newton_iterations = &
    this%cumulative_newton_iterations + sum_newton_iterations
  this%cumulative_linear_iterations = &
    this%cumulative_linear_iterations + sum_linear_iterations
  this%cumulative_wasted_linear_iterations = &
    this%cumulative_wasted_linear_iterations + sum_wasted_linear_iterations
  this%cumulative_time_step_cuts = &
    this%cumulative_time_step_cuts + icut

  this%num_newton_iterations = num_newton_iterations
  this%num_linear_iterations = num_linear_iterations  
  
  ! print screen output
  ! (TODO -checking, fmyuan) after updating TH mode, the following 'residual_vec' is NULL
  !            when surf_subsurface simulation is on
  !call SNESGetFunction(solver%snes,residual_vec,PETSC_NULL_FUNCTION, &
  !                     PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  !call VecNorm(residual_vec,NORM_2,fnorm,ierr);CHKERRQ(ierr)
  !call VecNorm(residual_vec,NORM_INFINITY,inorm,ierr);CHKERRQ(ierr)
  call VecNorm(process_model%residual_vec,NORM_2,fnorm,ierr);CHKERRQ(ierr)
  call VecNorm(process_model%residual_vec,NORM_INFINITY,inorm,ierr);CHKERRQ(ierr)
  if (option%print_screen_flag) then
      write(*, '(/," Step ",i6," Time= ",1pe12.5," Dt= ",1pe12.5, &
           & " [",a,"]", " snes_conv_reason: ",i4,/,"  newton = ",i3, &
           & " [",i8,"]", " linear = ",i5," [",i10,"]"," cuts = ",i2, &
           & " [",i4,"]")') &
           this%steps, &
           this%target_time/tconv, &
           this%dt/tconv, &
           trim(tunit),snes_reason,sum_newton_iterations, &
           this%cumulative_newton_iterations,sum_linear_iterations, &
           this%cumulative_linear_iterations,icut, &
           this%cumulative_time_step_cuts


    if (associated(process_model%realization_base%discretization%grid)) then
       scaled_fnorm = fnorm/process_model%realization_base% &
                        discretization%grid%nmax 
    else
       scaled_fnorm = fnorm
    endif

    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'("  --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif

!fmy: begining
#ifndef CLM_PFLOTRAN
! the following output produces a large ascii file if coupled with CLM
  if (option%print_file_flag) then
    write(option%fid_out, '(" Step ",i6," Time= ",1pe12.5," Dt= ",1pe12.5, &
      & " [",a,"]"," snes_conv_reason: ",i4,/,"  newton = ",i3, &
      & " [",i8,"]", " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      this%steps, &
      this%target_time/tconv, &
      this%dt/tconv, &
      trim(tunit),snes_reason,sum_newton_iterations, &
      this%cumulative_newton_iterations,sum_linear_iterations, &
      this%cumulative_linear_iterations,icut, &
      this%cumulative_time_step_cuts
  endif  
#endif
!fmy: ending

  
  option%time = this%target_time
  call process_model%FinalizeTimestep()
  
  if (this%print_ekg .and. OptionPrintToFile(option)) then
100 format(a32," TIMESTEP ",i10,2es16.8,a,i3,i5,i3,i5,i5,i10)
    write(IUNIT_EKG,100) trim(this%name), this%steps, this%target_time/tconv, &
      this%dt/tconv, trim(tunit), &
      icut, this%cumulative_time_step_cuts, &
      sum_newton_iterations, this%cumulative_newton_iterations, &
      sum_linear_iterations, this%cumulative_linear_iterations
  endif
  
  if (option%print_screen_flag) print *, ""  
  
end subroutine TimestepperBEStepDT

! ************************************************************************** !

subroutine TimestepperBECheckpointBinary(this,viewer,option)
  ! 
  ! Checkpoints parameters/variables associated with
  ! a time stepper.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/25/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"
  
  class(timestepper_BE_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option
  
  class(stepper_BE_header_type), pointer :: header
  type(stepper_BE_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  PetscErrorCode :: ierr

  bagsize = size(transfer(dummy_header,dummy_char))

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call TimestepperBERegisterHeader(this,bag,header)
  call TimestepperBESetHeader(this,bag,header)
  call PetscBagView(bag,viewer,ierr);CHKERRQ(ierr)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

end subroutine TimestepperBECheckpointBinary

! ************************************************************************** !

subroutine TimestepperBERegisterHeader(this,bag,header)
  ! 
  ! Register header entries.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/30/13
  ! 

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_BE_type) :: this
  class(stepper_BE_header_type) :: header
  PetscBag :: bag
  
  PetscErrorCode :: ierr
  
  call PetscBagRegisterInt(bag,header%cumulative_newton_iterations,0, &
                           "cumulative_newton_iterations","", &
                           ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%cumulative_linear_iterations,0, &
                           "cumulative_linear_iterations","", &
                           ierr);CHKERRQ(ierr)
! need to add cumulative wasted linear iterations
  call PetscBagRegisterInt(bag,header%num_newton_iterations,0, &
                           "num_newton_iterations","",ierr);CHKERRQ(ierr)

  call TimestepperBaseRegisterHeader(this,bag,header)
  
end subroutine TimestepperBERegisterHeader

! ************************************************************************** !

subroutine TimestepperBESetHeader(this,bag,header)
  ! 
  ! Sets values in checkpoint header.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/25/13
  ! 

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_BE_type) :: this
  class(stepper_BE_header_type) :: header
  PetscBag :: bag
  
  PetscErrorCode :: ierr
  
  header%cumulative_newton_iterations = this%cumulative_newton_iterations
  header%cumulative_linear_iterations = this%cumulative_linear_iterations
  header%num_newton_iterations = this%num_newton_iterations

  call TimestepperBaseSetHeader(this,bag,header)
  
end subroutine TimestepperBESetHeader

! ************************************************************************** !

subroutine TimestepperBERestartBinary(this,viewer,option)
  ! 
  ! Checkpoints parameters/variables associated with
  ! a time stepper.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/25/13
  ! 

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_BE_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option
  
  class(stepper_BE_header_type), pointer :: header
  type(stepper_BE_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  PetscErrorCode :: ierr

  bagsize = size(transfer(dummy_header,dummy_char))
  
  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call TimestepperBERegisterHeader(this,bag,header)
  call PetscBagLoad(viewer,bag,ierr);CHKERRQ(ierr)
  call TimestepperBEGetHeader(this,header)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

end subroutine TimestepperBERestartBinary

! ************************************************************************** !

#if defined(PETSC_HAVE_HDF5)
subroutine TimestepperBECheckpointHDF5(this, chk_grp_id, option)
  !
  ! Checkpoints parameters/variables associated with
  ! a time stepper.
  !
  ! Author: Gautam Bisht
  ! Date: 07/30/15
  !
  use Option_module
  use hdf5
  use Checkpoint_module, only : CheckPointWriteIntDatasetHDF5
  use Checkpoint_module, only : CheckPointWriteRealDatasetHDF5

  implicit none
  
  class(timestepper_BE_type) :: this
  PetscInt :: chk_grp_id
  type(option_type) :: option

#if defined(SCORPIO_WRITE)
  integer :: h5_chk_grp_id
  integer, pointer :: dims(:)
  integer, pointer :: start(:)
  integer, pointer :: stride(:)
  integer, pointer :: length(:)
  integer :: timestepper_grp_id
#else
  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
  integer(HID_T) :: timestepper_grp_id
  integer(HID_T) :: h5_chk_grp_id
#endif

  PetscMPIInt :: dataset_rank
  character(len=MAXSTRINGLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: string
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)
  PetscReal, pointer :: real_array(:)
  PetscMPIInt :: hdf5_err

  string = "Timestepper"
  h5_chk_grp_id = chk_grp_id
  call h5gcreate_f(h5_chk_grp_id, string, timestepper_grp_id, &
                   hdf5_err, OBJECT_NAMELEN_DEFAULT_F)

  allocate(start(1))
  allocate(dims(1))
  allocate(length(1))
  allocate(stride(1))
  allocate(int_array(1))
  allocate(real_array(1))

  dataset_rank = 1
  dims(1) = ONE_INTEGER
  start(1) = 0
  length(1) = ONE_INTEGER
  stride(1) = ONE_INTEGER

  dataset_name = "Cumulative_newton_iterations" // CHAR(0)
  int_array(1) = this%cumulative_newton_iterations
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Cumulative_linear_iterations" // CHAR(0)
  int_array(1) = this%cumulative_linear_iterations
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Num_newton_iterations" // CHAR(0)
  int_array(1) = this%num_newton_iterations
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Time" // CHAR(0)
  real_array(1) = this%target_time
  call CheckPointWriteRealDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, real_array, option)

  dataset_name = "Dt" // CHAR(0)
  real_array(1) = this%dt
  call CheckPointWriteRealDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, real_array, option)

  dataset_name = "Prev_dt" // CHAR(0)
  real_array(1) = this%prev_dt
  call CheckPointWriteRealDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, real_array, option)

  dataset_name = "Num_steps" // CHAR(0)
  int_array(1) = this%steps
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Cumulative_time_step_cuts" // CHAR(0)
  int_array(1) = this%cumulative_time_step_cuts
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Num_constant_time_steps" // CHAR(0)
  int_array(1) = this%num_constant_time_steps
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Num_contig_revert_due_to_sync" // CHAR(0)
  int_array(1) = this%num_contig_revert_due_to_sync
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Revert_dt" // CHAR(0)
  int_array(1) = ZERO_INTEGER
  if (this%revert_dt) int_array(1) = ONE_INTEGER
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  call h5gclose_f(timestepper_grp_id, hdf5_err)

  deallocate(start)
  deallocate(dims)
  deallocate(length)
  deallocate(stride)
  deallocate(int_array)
  deallocate(real_array)

end subroutine TimestepperBECheckpointHDF5

! ************************************************************************** !

subroutine TimestepperBERestartHDF5(this, chk_grp_id, option)
  !
  ! Restarts parameters/variables associated with
  ! a time stepper.
  !
  ! Author: Gautam Bisht
  ! Date: 08/16/15
  !
  use Option_module
  use hdf5
  use Checkpoint_module, only : CheckPointReadIntDatasetHDF5
  use Checkpoint_module, only : CheckPointReadRealDatasetHDF5

  implicit none
  
  class(timestepper_BE_type) :: this
  PetscInt :: chk_grp_id
  type(option_type) :: option

#if defined(SCORPIO_WRITE)
  integer :: h5_chk_grp_id
  integer, pointer :: dims(:)
  integer, pointer :: start(:)
  integer, pointer :: stride(:)
  integer, pointer :: length(:)
  integer :: timestepper_grp_id
#else
  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
  integer(HID_T) :: timestepper_grp_id
  integer(HID_T) :: h5_chk_grp_id
#endif

  PetscMPIInt :: dataset_rank
  character(len=MAXSTRINGLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: string
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)
  PetscReal, pointer :: real_array(:)
  PetscMPIInt :: hdf5_err

  string = "Timestepper"
  h5_chk_grp_id = chk_grp_id
  call h5gopen_f(h5_chk_grp_id, string, timestepper_grp_id, hdf5_err)

  allocate(start(1))
  allocate(dims(1))
  allocate(length(1))
  allocate(stride(1))
  allocate(int_array(1))
  allocate(real_array(1))

  dataset_rank = 1
  dims(1) = ONE_INTEGER
  start(1) = 0
  length(1) = ONE_INTEGER
  stride(1) = ONE_INTEGER

  dataset_name = "Cumulative_newton_iterations" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                    dataset_rank, dims, start, length, &
                                    stride, int_array, option)
  this%cumulative_newton_iterations = int_array(1)

  dataset_name = "Cumulative_linear_iterations" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                    dataset_rank, dims, start, length, &
                                    stride, int_array, option)
  this%cumulative_linear_iterations = int_array(1)

  dataset_name = "Num_newton_iterations" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                    dataset_rank, dims, start, length, &
                                    stride, int_array, option)
  this%num_newton_iterations = int_array(1)

  dataset_name = "Time" // CHAR(0)
  call CheckPointReadRealDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, real_array, option)
  this%target_time = real_array(1)

  dataset_name = "Dt" // CHAR(0)
  call CheckPointReadRealDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, real_array, option)
  this%dt = real_array(1)

  dataset_name = "Prev_dt" // CHAR(0)
  call CheckPointReadRealDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, real_array, option)
  this%prev_dt = real_array(1)

  dataset_name = "Num_steps" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)
  this%steps = int_array(1)

  dataset_name = "Cumulative_time_step_cuts" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)
  this%cumulative_time_step_cuts = int_array(1)

  dataset_name = "Num_constant_time_steps" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)
  this%num_constant_time_steps = int_array(1)

  dataset_name = "Num_contig_revert_due_to_sync" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)
  this%num_contig_revert_due_to_sync = int_array(1)

  dataset_name = "Revert_dt" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)
  this%revert_dt = (int_array(1) == ONE_INTEGER)

  call h5gclose_f(timestepper_grp_id, hdf5_err)

  deallocate(start)
  deallocate(dims)
  deallocate(length)
  deallocate(stride)
  deallocate(int_array)
  deallocate(real_array)

end subroutine TimestepperBERestartHDF5
#endif

! ************************************************************************** !

subroutine TimestepperBEGetHeader(this,header)
  ! 
  ! Gets values in checkpoint header.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/25/13
  ! 

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_BE_type) :: this
  class(stepper_BE_header_type) :: header
  
  this%cumulative_newton_iterations = header%cumulative_newton_iterations
  this%cumulative_linear_iterations = header%cumulative_linear_iterations
  this%num_newton_iterations = header%num_newton_iterations

  call TimestepperBaseGetHeader(this,header)
  
end subroutine TimestepperBEGetHeader

! ************************************************************************** !

subroutine TimestepperBEReset(this)
  ! 
  ! Zeros timestepper object members.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/20/14
  ! 

  implicit none

  class(timestepper_BE_type) :: this
  
  this%cumulative_newton_iterations = 0
  this%cumulative_linear_iterations = 0
  this%num_newton_iterations = 0

  call TimestepperBaseReset(this)
  
end subroutine TimestepperBEReset

! ************************************************************************** !

subroutine TimestepperBEPrintInfo(this,option)
  ! 
  ! Prints settings for base timestepper.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Option_module

  implicit none

  class(timestepper_BE_type) :: this
  type(option_type) :: option
  
  call TimestepperBasePrintInfo(this,option)
  call SolverPrintNewtonInfo(this%solver,this%name,option)
  call SolverPrintLinearInfo(this%solver,this%name,option)
  
end subroutine TimestepperBEPrintInfo

! ************************************************************************** !

subroutine TimestepperBEInputRecord(this)
  ! 
  ! Prints information about the time stepper to the input record.
  ! To get a## format, must match that in simulation types.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  ! 
  
  implicit none
  
  class(timestepper_BE_type) :: this

  PetscInt :: id
  character(len=MAXWORDLENGTH) :: word
   
  id = INPUT_RECORD_UNIT
  
  write(id,'(a29)',advance='no') 'pmc timestepper: '
  write(id,'(a)') this%name

  write(id,'(a29)',advance='no') 'initial timestep size: '
  write(word,*) this%dt_init
  write(id,'(a)') trim(adjustl(word)) // ' sec'

end subroutine TimestepperBEInputRecord

! ************************************************************************** !

recursive subroutine TimestepperBEFinalizeRun(this,option)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  ! 

  use Option_module
  
  implicit none
  
  class(timestepper_BE_type) :: this
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
#ifdef DEBUG
  call printMsg(option,'TimestepperBEFinalizeRun()')
#endif
  
  if (OptionPrintToScreen(option)) then
    write(*,'(/,a," TS BE steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            trim(this%name), &
            this%steps, &
            this%cumulative_newton_iterations, &
            this%cumulative_linear_iterations, &
            this%cumulative_time_step_cuts
    write(string,'(i12)') this%cumulative_wasted_linear_iterations
    write(*,'(a)') trim(this%name) // ' TS BE Wasted Linear Iterations = ' // &
      trim(adjustl(string))
    write(string,'(f12.1)') this%cumulative_solver_time
    write(*,'(a)') trim(this%name) // ' TS BE SNES time = ' // &
      trim(adjustl(string)) // ' seconds'
  endif
  
end subroutine TimestepperBEFinalizeRun

! ************************************************************************** !

subroutine TimestepperBEStrip(this)
  ! 
  ! Deallocates members of a time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  ! 

  implicit none
  
  class(timestepper_BE_type) :: this
  
  call TimestepperBaseStrip(this)
  call SolverDestroy(this%solver)

  if (associated(this%tfac)) deallocate(this%tfac)
  nullify(this%tfac)
  
end subroutine TimestepperBEStrip

! ************************************************************************** !

subroutine TimestepperBEDestroy(this)
  ! 
  ! Deallocates a time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  ! 

  implicit none
  
  class(timestepper_BE_type) :: this
  
  call TimestepperBEStrip(this)
  
!  deallocate(this)
!  nullify(this)
  
end subroutine TimestepperBEDestroy

end module Timestepper_BE_class
