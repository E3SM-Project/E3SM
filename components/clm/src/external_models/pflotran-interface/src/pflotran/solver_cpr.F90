module Solver_CPR_module

#include "petsc/finclude/petscts.h"
  use petscts
  use CPR_Preconditioner_module
  use Option_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: SolverCPRInit, &
            SolverCPRRead, & 
            SolverCPRInitializeStorage

contains

subroutine SolverCPRRead(stash, input, option, ierr)
  ! 
  ! Reads the CPR_OPTIONS card
  !
  ! Author: Daniel Stone
  ! Date: March 2018
  ! 
  use Input_Aux_module
  use String_module

  implicit none
  type(cpr_pc_type) :: stash
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscErrorCode :: ierr

  character(len=MAXWORDLENGTH) :: keyword, word, word2
  character(len=MAXSTRINGLENGTH) :: string

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CPR OPTIONS')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
      case('CPRT2_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'cprT2_type','CPR OPTIONS')
        call StringToUpper(word)
        select case(trim(word))
          case('SAILS')
            stash%T2_type = 'SAILS'
          case('PBJ')
            stash%T2_type = 'PBJ'
          case('NONE')
            stash%T2_type = 'NONE'
          case('PCASM')
            stash%T2_type = 'PCASM'
          case('PCGASM')
            stash%T2_type = 'PCGASM'
          case('PILUT')
            stash%T2_type = 'PILUT'
          case('EUCLID')
            stash%T2_type = 'EUCLID'
          case('ILU')
            ! this will cause crash if using more than 1 proc
            stash%T2_type = 'ILU'
          case default
            option%io_buffer  = 'CPR T2 Preconditioner type: ' // trim(word) // &
                                ' unknown.'
            call printErrMsg(option)
        end select

      case('CPRT1_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'cprT1_type','CPR OPTIONS')
        call StringToUpper(word)
        select case(trim(word))
          case('RICHARDSON')
            !option%CPRT1_type = 'RICHARDSON'
            stash%T1_type = 'RICHARDSON'
          case('FGMRES')
            !option%CPRT1_type = 'FGMRES'
            stash%T1_type = 'FGMRES'
          case default
            option%io_buffer  = 'CPR T1 KSP type: ' // trim(word) // &
                                ' unknown.'
            call printErrMsg(option)
        end select

      case('CPR_EXTRACTION_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'cpr_extraction_type','CPR OPTIONS')
        call StringToUpper(word)
        select case(trim(word))
          case('QIMPES_VARIABLE')
            stash%T1_type = 'QIMPES_VARIABLE'
          case('QIMPES_VARIABLE_FORCE')
            stash%T1_type = 'QIMPES_VARIABLE_FORCE'
          case default
            option%io_buffer  = 'CPR Extraction type: ' // trim(word) // &
                                ' unknown.'
            call printErrMsg(option)
        end select

      case('CPR_EX_OFFSET')
        call InputReadInt(input,option, stash%exrow_offset)

      case('T2FILLIN')
        call InputReadInt(input,option, stash%t2fillin)

      case('T2ASMOVERLAP')
        call InputReadInt(input,option, stash%asmoverlap)

      case('T2SHIFTINBLOCKS')
        stash%t2shiftinblocks= PETSC_TRUE

      case('ASMINPLACE')
        stash%asmfactorinplace = PETSC_TRUE

      case('CPROUTPUTDEBUG')
        stash%mv_output = PETSC_TRUE

      case('T2_ZEROING')
        stash%t2_zeroing= PETSC_TRUE

      case('CPRZEROING')
        stash%zeroing = PETSC_TRUE

      case('CPRSKIPT1')
        stash%skip_T1= PETSC_TRUE

      case('CPRAMGREPORT')
        stash%amg_report = PETSC_TRUE

      case('CPR_MANUAL_AMG_CONFIG')
        stash%amg_manual = PETSC_TRUE

      case('CPRGAMG')
        stash%useGAMG = PETSC_TRUE

      ! here is a sub card for setting boomeramg options for within
      ! the CPR PC, ONLY.
      ! Note the lack of a flow/transport prefix.
      case('CPR_HYPRE_OPTIONS')
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,keyword,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'CPR OPTIONS, HYPRE options')
          call StringToUpper(keyword)
          select case(trim(keyword))
            case('BOOMERAMG_CYCLE_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG cycle type', &
                                 'CPR OPTIONS, HYPRE options')
              call StringToLower(word)
              string = '-pc_hypre_boomeramg_cycle_type'
              select case(trim(word))
                case('V')
                  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                            trim(string),'1', &
                                            ierr);CHKERRQ(ierr)
                case('W')
                  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                            trim(string),'2', &
                                            ierr);CHKERRQ(ierr)
                case default
                  option%io_buffer  = 'HYPRE BoomerAMG cycle type: ' &
                                      // trim(word) // ' unknown.'
                  call printErrMsg(option)
              end select
            case('BOOMERAMG_MAX_LEVELS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG maximum levels', &
                                 'CPR OPTIONS, HYPRE options')
              string =  '-pc_hypre_boomeramg_max_levels'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_MAX_ITER')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG maximum iterations', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_max_iter'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_TOL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG convergence tolerance', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_tol'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_TRUNCFACTOR')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG interpolation truncation factor', &
                                 'CPR OPTIONS, HYPRE options')
              string =  '-pc_hypre_boomeramg_truncfactor'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_AGG_NL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG # levels aggressive coarsening', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_agg_nl'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_AGG_NUM_PATHS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                'BoomerAMG # paths for aggressive coarsening', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_agg_num_paths'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_STRONG_THRESHOLD')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                'BoomerAMG threshold for strong connectivity', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_strong_threshold'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_GRID_SWEEPS_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                         'BoomerAMG number of grid sweeps up and down cycles', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_grid_sweeps_all'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_GRID_SWEEPS_DOWN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                'BoomerAMG number of grid sweeps down cycles', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_grid_sweeps_down'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_GRID_SWEEPS_UP')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                  'BoomerAMG number of grid sweeps up cycles', &
                                 'CPR OPTIONS, HYPRE options')
              !string = trim(prefix) // 'pc_hypre_boomeramg_grid_sweeps_up'
              string = '-pc_hypre_boomeramg_grid_sweeps_up'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_GRID_SWEEPS_COARSE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                           'BoomerAMG number of grid sweeps for coarse level', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_grid_sweeps_coarse'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_TYPE_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                           'BoomerAMG relaxation type for up and down cycles', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_relax_type_all'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_TYPE_DOWN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                  'BoomerAMG relaxation type for down cycles', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_relax_type_down'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_TYPE_UP')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG relaxation type for up cycles', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_relax_type_up'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_TYPE_COARSE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG relaxation type for coarse grids', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_relax_type_coarse'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_WEIGHT_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG relaxation weight for all levels', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_relax_weight_all'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_WEIGHT_LEVEL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputReadWord(input,option,word2,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG relaxation weight for a level', &
                                 'CPR OPTIONS, HYPRE options')
              word = trim(word) // ' ' // trim(word2)
              string = '-pc_hypre_boomeramg_relax_weight_level'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_OUTER_RELAX_WEIGHT_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                           'BoomerAMG outer relaxation weight for all levels', &
                                 'CPR OPTIONS, HYPRE options')
              string =  '-pc_hypre_boomeramg_outer_relax_weight_all'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_OUTER_RELAX_WEIGHT_LEVEL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputReadWord(input,option,word2,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                              'BoomerAMG outer relaxation weight for a level', &
                                 'CPR OPTIONS, HYPRE options')
              word = trim(word) // ' ' // trim(word2)
              string = '-pc_hypre_boomeramg_outer_relax_weight_level'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_NO_CF')
              !string = trim(prefix) // 'pc_hypre_boomeramg_no_CF'
              string = '-pc_hypre_boomeramg_no_CF'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),'',ierr);CHKERRQ(ierr)
            case('BOOMERAMG_MEASURE_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG measure type', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_measure_type'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_COARSEN_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG coarsen type', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_coarsen_type'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_INTERPOLATION_TYPE','BOOMERAMG_INTERP_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG interpolation type', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_interp_type'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_NODAL_COARSEN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG set nodal coarsening', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_nodal_coarsen'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),'',ierr);CHKERRQ(ierr)
            case('BOOMERAMG_NODAL_RELAXATION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG nodal relaxation via Schwarz', &
                                 'CPR OPTIONS, HYPRE options')
              string = '-pc_hypre_boomeramg_nodal_relaxation'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),'',ierr);CHKERRQ(ierr)
            case default
              option%io_buffer  = 'HYPRE option: ' // trim(keyword) // &
                                  ' unknown.'
              call printErrMsg(option)
          end select
        enddo
    case default
      option%io_buffer  = 'CPR preconditioner option: ' // trim(keyword) // &
                          ' unknown.'
      call printErrMsg(option)
    end select
  enddo

end subroutine SolverCPRRead

! ************************************************************************** !

subroutine SolverCPRInit(J, stash, pcin, ierr, option)
  ! 
  ! Do needed set up and call routine to build a CPR 
  ! preconditioner.
  !
  ! Author: Daniel Stone
  ! Date: Oct 2017 - March 2018
  ! 
  implicit none
  !type(solver_type) :: solver
  Mat :: J
  type(cpr_pc_type) :: stash
  PC :: pcin
  MPI_Comm :: C
  PetscErrorCode :: ierr
  type(option_type) :: option

  call PetscObjectGetComm(pcin, C, ierr); CHKERRQ(ierr)

  call CPRmake(pcin, stash, C, ierr, option)

  ! set the A matrix in the stash to J:
  stash%A = J

end subroutine SolverCPRInit

! ************************************************************************** !

subroutine SolverCPRInitializeStorage(ctx)
  ! MUST CALL THIS before doing anything with a cpr_pc_type

  ! 
  ! initialize all the noncomplicated members of an
  ! cpr_pc_type object
  !
  ! Author:  Daniel Stone 
  ! Date: Oct 2017 - March 2018
  ! 

  implicit none

  type(cpr_pc_type) :: ctx

  ! ensure that first run flags are set correctly
  ctx%firstT1Call = PETSC_TRUE
  ctx%firstT2Call = PETSC_TRUE

  ctx%T1_type = "NONE"
  ctx%T2_type = "Jacobi"
  ctx%extract_type = "QIMPES_VARIABLE"

  ctx%asmfactorinplace = PETSC_FALSE
  ctx%t2shiftinblocks = PETSC_FALSE
  ctx%mv_output= PETSC_FALSE
  ctx%t2_zeroing= PETSC_FALSE
  ctx%skip_T1= PETSC_FALSE
  ctx%amg_report = PETSC_FALSE
  ctx%amg_manual = PETSC_FALSE

  ctx%zeroing = PETSC_FALSE
  ctx%useGAMG= PETSC_FALSE


  ctx%t2fillin = 0
  ctx%timescalled = 0

  ctx%asmoverlap = 0

  ctx%exrow_offset = 0

  nullify(ctx%option)


end subroutine SolverCPRInitializeStorage 

! ************************************************************************** !

end module Solver_CPR_module
