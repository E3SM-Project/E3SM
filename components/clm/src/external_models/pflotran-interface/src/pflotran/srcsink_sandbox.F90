module SrcSink_Sandbox_module

  use SrcSink_Sandbox_Base_class
  use SrcSink_Sandbox_Mass_Rate_class
  use SrcSink_Sandbox_Downreg_class
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  class(srcsink_sandbox_base_type), pointer, public :: ss_sandbox_list
  PetscBool :: print_mass_balance

  interface SSSandboxRead
    module procedure SSSandboxRead1
    module procedure SSSandboxRead2
  end interface
  
  interface SSSandboxDestroyList
    module procedure SSSandboxDestroyList1
    module procedure SSSandboxDestroyList2
  end interface
  
  public :: SSSandboxInit, &
            SSSandboxRead, &
            SSSandboxSetup, &
            SSSandboxUpdate, &
            SSSandbox, &
            SSSandboxDestroyList

contains

! ************************************************************************** !

subroutine SSSandboxInit(option)
  ! 
  ! Initializes the sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  implicit none
  type(option_type) :: option

  if (associated(ss_sandbox_list)) then
    call SSSandboxDestroyList()
  endif
  nullify(ss_sandbox_list)
  print_mass_balance = PETSC_FALSE

end subroutine SSSandboxInit

! ************************************************************************** !

subroutine SSSandboxRead1(input,option)
  ! 
  ! Reads input deck for source/sink sandbox parameters
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  type(input_type), pointer :: input
  type(option_type) :: option

  call SSSandboxRead(ss_sandbox_list,input,option)

end subroutine SSSandboxRead1

! ************************************************************************** !

subroutine SSSandboxRead2(local_sandbox_list,input,option)
  ! 
  ! Reads input deck for src/sink sandbox parameters
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  class(srcsink_sandbox_base_type), pointer :: local_sandbox_list  
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  class(srcsink_sandbox_base_type), pointer :: new_sandbox, cur_sandbox
  
  ! Ensure that transport is not being simulated as we have no way for 
  ! introducting solutes.
  if (option%ntrandof > 0) then
    option%io_buffer = 'Reactive transport may not be simulated when a &
      &SOURCE_SINK_SANDBOX exists in the input file since no source/sink &
      &capability exists in the source/sink sandbox for solute mass.'
    call printErrMsg(option)
  endif

  nullify(new_sandbox)
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','SOURCE_SINK_SANDBOX')
    call StringToUpper(word)   

    select case(trim(word))
      case('MASS_RATE')
        new_sandbox => MassRateCreate()
      case('MASS_RATE_DOWNREGULATED')
        new_sandbox => DownregCreate()
      case('MASS_BALANCE')
        print_mass_balance = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(word,'SRCSINK_SANDBOX',option)
    end select
    
    if (associated(new_sandbox)) then
      call new_sandbox%ReadInput(input,option)
      if (.not.associated(local_sandbox_list)) then
        local_sandbox_list => new_sandbox
      else
        cur_sandbox => local_sandbox_list
        do
          if (.not.associated(cur_sandbox%next)) exit
          cur_sandbox => cur_sandbox%next
        enddo
        cur_sandbox%next => new_sandbox
      endif
    endif
    nullify(new_sandbox)
  enddo
  
end subroutine SSSandboxRead2


! ************************************************************************** !

subroutine SSSandboxSetup(grid,option,output_option)
  ! 
  ! Calls all the initialization routines for all source/sinks in
  ! the sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module
  use Output_Aux_module
  use Grid_module
  
  implicit none
  
  type(grid_type) :: grid
  type(option_type) :: option
  type(output_option_type) :: output_option
  
  class(srcsink_sandbox_base_type), pointer :: cur_sandbox  
  class(srcsink_sandbox_base_type), pointer :: prev_sandbox  
  class(srcsink_sandbox_base_type), pointer :: next_sandbox  
  PetscBool :: exists

  ! sandbox source/sinks
  cur_sandbox => ss_sandbox_list
  nullify(prev_sandbox)
  do
    if (.not.associated(cur_sandbox)) exit
    next_sandbox => cur_sandbox%next
    call cur_sandbox%Setup(grid,option)
    ! destory if not on process
    if (.not.Initialized(cur_sandbox%local_cell_id)) then
      if (associated(prev_sandbox)) then
        prev_sandbox%next => next_sandbox
      else
        ss_sandbox_list => next_sandbox
      endif
      nullify(cur_sandbox%next)
      call SSSandboxDestroy(cur_sandbox)
    else
      if (print_mass_balance) then
        allocate(cur_sandbox%instantaneous_mass_rate(option%nflowdof))
        cur_sandbox%instantaneous_mass_rate = 0.d0
        allocate(cur_sandbox%cumulative_mass(option%nflowdof))
        cur_sandbox%cumulative_mass = 0.d0
      endif    
    endif
    if (associated(cur_sandbox)) prev_sandbox => cur_sandbox
    cur_sandbox => next_sandbox
  enddo
  
  if (print_mass_balance) then
    call SSSandboxOutputHeader(ss_sandbox_list,grid,option,output_option)
  endif

end subroutine SSSandboxSetup

! ************************************************************************** !

subroutine SSSandbox(residual,Jacobian,compute_derivative, &
                     grid,material_auxvars,option)
  ! 
  ! Evaluates source/sink term storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module
  use Grid_module
  use Material_Aux_class, only: material_auxvar_type
  
  implicit none

  PetscBool :: compute_derivative
  Vec :: residual
  Mat :: Jacobian
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(grid_type) :: grid
  type(option_type) :: option
  
  PetscReal, pointer :: r_p(:)
  PetscReal :: res(option%nflowdof)
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  class(srcsink_sandbox_base_type), pointer :: cur_srcsink
  PetscInt :: i, local_id, ghosted_id, istart, iend
  PetscReal :: aux_real(0)
  PetscErrorCode :: ierr
  
  if (.not.compute_derivative) then
    call VecGetArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif
  
  cur_srcsink => ss_sandbox_list
  do
    if (.not.associated(cur_srcsink)) exit
    local_id = cur_srcsink%local_cell_id
    ghosted_id = grid%nL2G(local_id)
    res = 0.d0
    Jac = 0.d0
    call cur_srcsink%Evaluate(res,Jac,compute_derivative, &
                              material_auxvars(ghosted_id), &
                              aux_real,option)
    if (compute_derivative) then
      call MatSetValuesBlockedLocal(Jacobian,1,ghosted_id-1,1, &
                                    ghosted_id-1,Jac,ADD_VALUES, &
                                    ierr);CHKERRQ(ierr)
    else
      iend = local_id*option%nflowdof
      istart = iend - option%nflowdof + 1
      r_p(istart:iend) = r_p(istart:iend) + res
    endif
    cur_srcsink => cur_srcsink%next
  enddo
  
  if (.not.compute_derivative) then
    call VecRestoreArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine SSSandbox

! ************************************************************************** !

subroutine SSSandboxUpdate(sandbox_list,option,output_option)
  ! 
  ! Updates datasets associated with a sandbox, if they exist
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/22/15
  ! 
  use Option_module
  use Output_Aux_module

  implicit none

  class(srcsink_sandbox_base_type), pointer :: sandbox_list
  type(option_type) :: option
  type(output_option_type) :: output_option

  class(srcsink_sandbox_base_type), pointer :: cur_sandbox
  
  cur_sandbox => sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    call cur_sandbox%Update(option)
    cur_sandbox => cur_sandbox%next
  enddo 
  
  if (print_mass_balance) then
    call SSSandboxOutput(sandbox_list,option,output_option)
  endif

end subroutine SSSandboxUpdate

! ************************************************************************** !

function SSSandboxOutputFilename(option)
  ! 
  ! Generates filename for source/sink sandbox output
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/16

  use Option_module

  implicit none
  
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: SSSandboxOutputFilename
  character(len=MAXWORDLENGTH) :: word

  write(word,'(i6)') option%myrank
  SSSandboxOutputFilename = trim(option%global_prefix) // &
                            trim(option%group_prefix) // &
                            '-ss_mass-' // trim(adjustl(word)) // '.dat'
  
end function SSSandboxOutputFilename  

! ************************************************************************** !

subroutine SSSandboxOutputHeader(sandbox_list,grid,option,output_option)
  ! 
  ! Writes header for source/sink sandbox output
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/16

  use Option_module
  use Output_Aux_module
  use Grid_module
  use Utility_module, only : BestFloat
  
  implicit none
  
  class(srcsink_sandbox_base_type), pointer :: sandbox_list
  type(grid_type) :: grid
  type(option_type) :: option
  type(output_option_type) :: output_option

  class(srcsink_sandbox_base_type), pointer :: cur_srcsink
  character(len=MAXSTRINGLENGTH) :: cell_string
  character(len=MAXWORDLENGTH) :: x_string, y_string, z_string
  character(len=MAXWORDLENGTH) :: units_string, variable_string
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid
  PetscInt :: icolumn, i
  PetscInt :: local_id, ghosted_id
  
  filename = SSSandboxOutputFilename(option)
  open(unit=IUNIT_TEMP,file=filename,action="write",status="replace")  
  
  if (output_option%print_column_ids) then
    icolumn = 1
  else
    icolumn = -1
  endif 
  
  write(IUNIT_TEMP,'(a)',advance="no") ' "Time [' // &
    trim(output_option%tunit) // ']"'

  cur_srcsink => sandbox_list
  do
    if (.not.associated(cur_srcsink)) exit
    local_id = cur_srcsink%local_cell_id
    ghosted_id = grid%nL2G(local_id)

    ! cell natural id
    write(cell_string,*) grid%nG2A(ghosted_id)
    cell_string = ' (' // trim(adjustl(cell_string)) // ')'
    ! coordinate of cell
    x_string = BestFloat(grid%x(ghosted_id),1.d4,1.d-2)
    y_string = BestFloat(grid%y(ghosted_id),1.d4,1.d-2)
    z_string = BestFloat(grid%z(ghosted_id),1.d4,1.d-2)
    cell_string = trim(cell_string) // &
             ' (' // trim(adjustl(x_string)) // &
             ' ' // trim(adjustl(y_string)) // &
             ' ' // trim(adjustl(z_string)) // ')'
    select case(option%iflowmode)
      case(TH_MODE)
        variable_string = ' Water'
        ! cumulative
        units_string = 'kg'
        call OutputWriteToHeader(IUNIT_TEMP,variable_string,units_string, &
                                 cell_string,icolumn)
        ! instantaneous
        units_string = 'kg/' // trim(adjustl(output_option%tunit))
        call OutputWriteToHeader(IUNIT_TEMP,variable_string,units_string, &
                                 cell_string,icolumn)
    end select
    select case(option%iflowmode)
      case(TH_MODE)
        variable_string = ' Gas Component'
        ! cumulative
        units_string = 'kg'
        call OutputWriteToHeader(IUNIT_TEMP,variable_string,units_string, &
                                 cell_string,icolumn)
        ! instantaneous
        units_string = 'kg/' // trim(adjustl(output_option%tunit))
        call OutputWriteToHeader(IUNIT_TEMP,variable_string,units_string, &
                                 cell_string,icolumn)
        variable_string = ' Energy'
        ! cumulative
        units_string = 'MJ'
        call OutputWriteToHeader(IUNIT_TEMP,variable_string,units_string, &
                                 cell_string,icolumn)
        ! instantaneous
        units_string = 'MJ/' // trim(adjustl(output_option%tunit))
        call OutputWriteToHeader(IUNIT_TEMP,variable_string,units_string, &
                                 cell_string,icolumn)
    end select
    cur_srcsink => cur_srcsink%next
  enddo
  
  close(IUNIT_TEMP)
  
end subroutine SSSandboxOutputHeader

! ************************************************************************** !

subroutine SSSandboxOutput(sandbox_list,option,output_option)
  ! 
  ! Writes output for for source/sink sandbox
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/16

  use Option_module
  use Output_Aux_module

  implicit none
  
  class(srcsink_sandbox_base_type), pointer :: sandbox_list
  type(option_type) :: option
  type(output_option_type) :: output_option
  
  class(srcsink_sandbox_base_type), pointer :: cur_srcsink
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: i
  PetscReal :: flow_dof_scale(3)  
  
  if (.not.associated(sandbox_list)) return

  flow_dof_scale = 1.d0
  select case(option%iflowmode)
    case(TH_MODE)
      flow_dof_scale(1) = FMWH2O
  end select  
  
100 format(100es16.8)

  filename = SSSandboxOutputFilename(option)
  open(unit=IUNIT_TEMP,file=filename,action="write",status="old", &
       position="append")

  ! this time is set at the end of the reactive transport step
  write(IUNIT_TEMP,100,advance="no") option%time / output_option%tconv
  
  cur_srcsink => sandbox_list
  do
    if (.not.associated(cur_srcsink)) exit
    do i = 1, option%nflowdof
      write(IUNIT_TEMP,100,advance="no") &
        cur_srcsink%cumulative_mass(i)*flow_dof_scale(i), &
        cur_srcsink%instantaneous_mass_rate(i)*flow_dof_scale(i)* &
          output_option%tconv
    enddo
    cur_srcsink => cur_srcsink%next
  enddo
  close(IUNIT_TEMP)
  
end subroutine SSSandboxOutput

! ************************************************************************** !

subroutine SSSandboxDestroy(sandbox)
  ! 
  ! Destroys arbitrary sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  implicit none

  class(srcsink_sandbox_base_type), pointer :: sandbox

  if (.not.associated(sandbox)) return

  call sandbox%Destroy()
  deallocate(sandbox)
  nullify(sandbox)

end subroutine SSSandboxDestroy


! ************************************************************************** !

subroutine SSSandboxDestroyList1()
  ! 
  ! Destroys master sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  implicit none

  call SSSandboxDestroyList(ss_sandbox_list)
  
end subroutine SSSandboxDestroyList1

! ************************************************************************** !

subroutine SSSandboxDestroyList2(local_sandbox_list)
  ! 
  ! Destroys arbitrary sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  implicit none

  class(srcsink_sandbox_base_type), pointer :: local_sandbox_list

  class(srcsink_sandbox_base_type), pointer :: cur_sandbox, prev_sandbox
  
  ! sandbox source/sinks
  cur_sandbox => local_sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    prev_sandbox => cur_sandbox%next
    call SSSandboxDestroy(cur_sandbox)
    cur_sandbox => prev_sandbox
  enddo  

end subroutine SSSandboxDestroyList2

end module SrcSink_Sandbox_module
