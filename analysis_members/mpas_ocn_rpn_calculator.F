! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ocn_rpn_calculator
!
!> \brief MPAS ocean analysis core member: rpn_calculator
!> \author Jon Woodring
!> \date   March 21, 2016
!> \details
!>  Flexible vector RPN calculator of MPAS fields for up to 2D fields.
!-----------------------------------------------------------------------
module ocn_rpn_calculator
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_dmpar
  use mpas_timekeeping
  use mpas_stream_manager

  use ocn_constants
  use ocn_diagnostics_routines

  implicit none
  private
  save

  ! Public parameters
  !--------------------------------------------------------------------

  ! Public member functions
  !--------------------------------------------------------------------
  public :: &
         ocn_init_rpn_calculator, &
         ocn_compute_rpn_calculator, &
         ocn_restart_rpn_calculator, &
         ocn_finalize_rpn_calculator

  ! Private module variables
  !--------------------------------------------------------------------

  type rpn_stack_value_type
    integer :: symbol_type
    integer :: number_of_dims

    type (field0DReal), pointer :: d0
    type (field1DReal), pointer :: d1
    type (field2DReal), pointer :: d2
  end type rpn_stack_value_type

  integer, parameter :: SYMBOL_NOT_FOUND = 0

  integer, parameter :: IS_OPERATOR = 10
  integer, parameter :: IS_VARIABLE = 100
  integer, parameter :: IS_TEMPORARY = 1000

  integer, parameter :: MAX_STACK_SIZE = StrKIND / 2

  character (len=1), dimension(8), parameter :: variable_names = &
    (/ 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h' /)
  character (len=3), dimension(4) :: operator_names = &
    (/ '*  ' , '+  ', '-  ', '/  ' /)
  integer, parameter :: MUL_OP = 1
  integer, parameter :: PLUS_OP = 2
  integer, parameter :: MINUS_OP = 3
  integer, parameter :: DIV_OP = 4
  ! TODO FIXME
  ! integer, parameter :: SUM_OP = 5

  character (len=1), dimension(4) :: expression_names = &
    (/ '1', '2', '3', '4' /)

  character (len=StrKIND), parameter :: VARIABLE_PREFIX = &
    'config_AM_rpnCalculator_variable_'
  character (len=StrKIND), parameter :: EXPRESSION_PREFIX = &
    'config_AM_rpnCalculator_expression_'
  character (len=StrKIND), parameter :: OUTPUT_PREFIX = &
    'config_AM_rpnCalculator_output_name_'

  character (len=StrKIND), parameter :: OUTPUT_STREAM_CONFIG = &
    'config_AM_rpnCalculator_output_stream'

  character (len=StrKIND), parameter :: NONE_TOKEN = 'none'

  character (len=StrKIND), parameter :: MPAS_CORE_NAME = 'MPAS-Ocean'

!***********************************************************************
contains



!***********************************************************************
! routine ocn_init_rpn_calculator
!
!> \brief Initialize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_init_rpn_calculator(domain, err)!{{{
  ! input variables

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: i, last, stack_pointer
  character (len=StrKIND) :: config, field_name
  character (len=StrKIND), pointer :: config_result
  type (rpn_stack_value_type) :: output_value
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE) :: stack

  ! start procedure
  err = 0

  ! typecheck all the expressions
  last = size(expression_names)
  do i = 1, last
    config = trim(EXPRESSION_PREFIX) // trim(expression_names(i))
    call mpas_pool_get_config(domain % configs, config, config_result)

    if (trim(config_result) /= trim(NONE_TOKEN)) then
      stack_pointer = -1 ! typecheck with an empty stack
      call eval_expression(domain, config_result, i, stack, stack_pointer)

      ! check the stack size
      if (stack_pointer /= 1) then
        call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
          'expression #' // trim(expression_names(i)) // &
          ' in the RPN calculator AM ' // &
          'resulted in the stack size not being equal to 1: ' // &
          'i.e., the return result of the expression should be the only ' // &
          'value on the stack after evaluation')
      end if

      ! check that it's a new value
      if (stack(stack_pointer) % symbol_type /= IS_TEMPORARY) then
        call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
          'expression #' // trim(expression_names(i)) // &
          ' in the RPN calculator AM did not calculate anything, ' // &
          ' i.e., it only pushed a variable onto the stack')
      end if

      ! rename the stack field and put in allFields pool
      config = trim(OUTPUT_PREFIX) // trim(expression_names(i))
      call mpas_pool_get_config(domain % configs, config, config_result)

      if (trim(config_result) == (NONE_TOKEN)) then
        call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
          'expression #' // trim(expression_names(i)) // &
          ' in the RPN calculator AM was set, but the output field name ' // &
          'for that expression was set to "none"')
      end if

      if (stack(1) % number_of_dims == 0) then
        stack(1) % d0 % fieldName = config_result
        call mpas_pool_add_field(domain % blocklist % allFields, &
          config_result, stack(1) % d0)
      else if (stack(1) % number_of_dims == 1) then
        stack(1) % d1 % fieldName = config_result
        call mpas_pool_add_field(domain % blocklist % allFields, &
          config_result, stack(1) % d1)
      else if (stack(1) % number_of_dims == 2) then
        stack(1) % d2 % fieldName = config_result
        call mpas_pool_add_field(domain % blocklist % allFields, &
          config_result, stack(1) % d2)
      else
        call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
          'the impossible happened, the dimensions of the result on the ' // &
          'stack, for expression #' // trim(expression_names(i)) // &
          ' was not between 0 and 2 in the RPN calculator AM')
      end if
      field_name = config_result

      ! put them in the stream if necessary
      call mpas_pool_get_config(domain % configs, &
        OUTPUT_STREAM_CONFIG, config_result)

      if (trim(config_result) /= trim(NONE_TOKEN)) then
        call mpas_stream_mgr_add_field(domain % streamManager, &
          config_result, field_name, ierr=err)
      end if
    end if
  end do

end subroutine ocn_init_rpn_calculator!}}}



!***********************************************************************
! routine ocn_compute_rpn_calculator
!
!> \brief Compute MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_compute_rpn_calculator(domain, timeLevel, err)!{{{
  ! input variables
  integer, intent(in) :: timeLevel

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: i, stack_pointer, last
  character (len=StrKIND) :: config
  character (len=StrKIND), pointer :: config_result
  type (rpn_stack_value_type) :: output_value
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE) :: stack
  type (field0DReal), pointer :: d0
  type (field1DReal), pointer :: d1, t1
  type (field2DReal), pointer :: d2, t2
  real (kind=RKIND), dimension(:), pointer :: s1
  real (kind=RKIND), dimension(:,:), pointer :: s2

  ! start procedure
  err = 0

  ! do all the expressions
  last = size(expression_names)
  do i = 1, last
    config = trim(EXPRESSION_PREFIX) // trim(expression_names(i))
    call mpas_pool_get_config(domain % configs, config, config_result)

    if (trim(config_result) /= trim(NONE_TOKEN)) then
      stack_pointer = 0 ! evaluate with an empty stack
      call eval_expression(domain, config_result, i, stack, stack_pointer)

      ! lookup the field and reassign pointers - then deallocate stack
      config = trim(OUTPUT_PREFIX) // trim(expression_names(i))
      call mpas_pool_get_config(domain % configs, config, config_result)
      if (stack(1) % number_of_dims == 0) then
        call mpas_pool_get_field(domain % blocklist % allFields, &
          config_result, d0, 1)
        d0 % scalar = stack(1) % d0 % scalar
        call mpas_deallocate_field(stack(1) % d0)
      else if (stack(1) % number_of_dims == 1) then
        call mpas_pool_get_field(domain % blocklist % allFields, &
          config_result, d1, 1)
        t1 => stack(1) % d1
        do while (associated(d1))
          s1 => d1 % array
          d1 % array => t1 % array
          t1 % array => s1

          d1 => d1 % next
          t1 => t1 % next
        end do
        call mpas_deallocate_field(stack(1) % d1)
      else
        call mpas_pool_get_field(domain % blocklist % allFields, &
          config_result, d2, 1)
        t2 => stack(1) % d2
        do while (associated(d2))
          s2 => d2 % array
          d2 % array => t2 % array
          t2 % array => s2

          d2 => d2 % next
          t2 => t2 % next
        end do
        call mpas_deallocate_field(stack(1) % d2)
      end if

    end if
  end do

end subroutine ocn_compute_rpn_calculator!}}}



!***********************************************************************
! routine ocn_restart_rpn_calculator
!
!> \brief Save restart for MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_restart_rpn_calculator(domain, err)!{{{
  ! input variables

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables

  ! start procedure
  err = 0

end subroutine ocn_restart_rpn_calculator!}}}



!***********************************************************************
! routine ocn_finalize_rpn_calculator
!
!> \brief Finalize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_finalize_rpn_calculator(domain, err)!{{{
  ! input variables

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables

  ! start procedure
  err = 0

end subroutine ocn_finalize_rpn_calculator!}}}

!
! local subroutines
!

!***********************************************************************
! routine eval_expression
!
!> \brief Given a character string, evaluate the stack expression
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a character string, evaluate the stack expression
!> and copy the bottom (top) of a 1-length stack into the target MPAS field.
!-----------------------------------------------------------------------
subroutine eval_expression (domain, expression, exp_number, &
    stack, stack_pointer)!{{{
  ! input variables
  character (len=StrKIND), intent(in) :: expression
  integer, intent(in) :: exp_number

  ! input/output variables
  type (domain_type), intent(inout) :: domain
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  ! output variables

  ! local variables
  integer :: symbol_type
  logical :: eol, typechecking
  character (len=StrKIND) :: symbol, remainder

  ! start procedure
  if(stack_pointer < 0) then
    typechecking = .true.
    stack_pointer = -1 - stack_pointer
  else
    typechecking = .false.
  end if

  eol = .false.
  remainder = expression

  ! get the first symbol
  call stack_token(symbol, remainder, eol)

  ! iterate over symbols
  do while(.not. eol)
    symbol_type = symbol_table(symbol)

    ! operator
    if ((symbol_type > IS_OPERATOR) .and. (symbol_type < IS_VARIABLE)) then
      call eval_operator(exp_number, &
        symbol_type - IS_OPERATOR, stack, stack_pointer, typechecking)
    ! variable
    else &
    if ((symbol_type > IS_VARIABLE) .and. (symbol_type < IS_TEMPORARY)) then
      call eval_variable(domain, exp_number, &
        symbol_type - IS_VARIABLE, stack, stack_pointer, typechecking)
    ! symbol not found
    else
      call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: "' // &
        trim(symbol) // '" found in expression #' // &
        trim(expression_names(exp_number)) // &
        ' in the RPN calculator AM was not found.')
    end if

    ! get the next symbol
    call stack_token(symbol, remainder, eol)
  end do

end subroutine eval_expression!}}}



!***********************************************************************
! routine stack_token
!
!> \brief Get the next stack token given a character string
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Parses a character string to get the next stack token to eval.
!-----------------------------------------------------------------------
subroutine stack_token(substr, next, eol)!{{{
  ! input variables

  ! input/output variables
  character (len=StrKIND), intent(inout) :: next

  ! output variables
  character (len=StrKIND), intent(out) :: substr
  logical, intent(out) :: eol

  ! local variables
  integer :: i
  character (len=StrKIND) :: copy

  ! make a copy
  copy = trim(next)

  ! if there's anything in it other than whitespace, pass through
  i = verify(copy, ' ')
  eol = i < 1
  if (eol) then
    return
  end if
  copy = trim(next(i:))

  ! find the first whitespace and split
  i = scan(copy, ' ')

  ! return that substring and the remainder
  if (i > 0) then
    substr = trim(copy(1:i-1))
    next = trim(copy(i+1:))
  else
    substr = trim(copy)
    next = ''
  end if

end subroutine stack_token!}}}



!***********************************************************************
! function symbol_table
!
!> \brief Tries to find the symbol in the symbol table and its value
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Will attempt to find a symbol in the symbol table. The value
!> results are dependent on the return code of the symbol table lookup.
!-----------------------------------------------------------------------
integer function symbol_table (symbol)!{{{
  ! input variables
  character (len=StrKIND), intent(in) :: symbol

  ! input/output variables

  ! local variables
  integer :: i, last

  ! start procedure

  ! check the operations
  last = size(variable_names)
  do i = 1, last
    if (trim(symbol) == trim(variable_names(i))) then
      symbol_table = IS_VARIABLE + i
      return
    end if
  end do

  ! check the variables
  last = size(operator_names)
  do i = 1, last
    if (trim(symbol) == trim(operator_names(i))) then
      symbol_table = IS_OPERATOR + i
      return
    end if
  end do

  ! else not found
  symbol_table = SYMBOL_NOT_FOUND

end function symbol_table!}}}



!***********************************************************************
! routine eval_operator
!
!> \brief Given a operator index number, put it on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a operator index number, put the result on the top of
!> the stack. It will combine whatever is on the stack to be able to
!> generate results and push them into the stack.
!-----------------------------------------------------------------------
subroutine eval_operator (exp_number, &
    op_index, stack, stack_pointer, type_checking)!{{{
  ! input variables
  integer, intent(in) :: exp_number, op_index
  logical, intent(in) :: type_checking

  ! input/output variables
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  ! output variables

  ! local variables

  ! start procedure
  if (op_index == MUL_OP) then
    call mul_operator(exp_number, stack, stack_pointer, type_checking)
  else if (op_index == PLUS_OP) then
    call plus_operator(exp_number, stack, stack_pointer, type_checking)
  else if (op_index == MINUS_OP) then
    call minus_operator(exp_number, stack, stack_pointer, type_checking)
  else if (op_index == DIV_OP) then
    call div_operator(exp_number, stack, stack_pointer, type_checking)
  ! TODO FIXME
  ! sum (and other reduces) needs to be fixed,
  ! because it is using (:) over decomposed dimensions, which is wrong
  !
  ! else if (op_index == SUM_OP) then
  !  call sum_operator(exp_number, stack, stack_pointer, type_checking)
  else
    call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
      'the impossible happened, tried to apply an unknown operator ' // &
      'in expression #' // trim(expression_names(exp_number)) // &
      ' in the RPN calculator AM')
  end if

end subroutine eval_operator!}}}



!***********************************************************************
! routine eval_variable
!
!> \brief Given a variable index number, put it on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a variable index number, put the result on the top of
!> the stack. This will look up the field names in the variable and look
!> it up from the framework to push the pointer onto the stack.
!-----------------------------------------------------------------------
subroutine eval_variable (domain, exp_number, &
    var_index, stack, stack_pointer, type_checking)!{{{
  ! input variables
  integer, intent(in) :: exp_number, var_index
  logical, intent(in) :: type_checking

  ! input/output variables
  type (domain_type), intent(inout) :: domain
  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  ! output variables

  ! local variables
  character (len=StrKIND) :: config
  character (len=StrKIND), pointer :: config_result
  type (mpas_pool_field_info_type) :: info

  ! start procedure
  config = trim(VARIABLE_PREFIX) // trim(variable_names(var_index))
  call mpas_pool_get_config(domain % configs, config, config_result)

  if (type_checking) then
    if (trim(config_result) == trim(NONE_TOKEN)) then
      call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'the MPAS field assigned to variable ' // &
        trim(variable_names(var_index)) // ' was evaluated, but it is ' // &
        'currently set to "none"')
    end if
  end if

  call mpas_pool_get_field_info &
    (domain % blocklist % allFields, config_result, info)

  ! check if it's real
  if (type_checking) then
    if (info % fieldType /= MPAS_POOL_REAL) then
      call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'the MPAS field "' // trim(config_result) // &
        '"assigned to variable ' // &
        trim(variable_names(var_index)) // ' in the RPN calculator AM is ' // &
        'not a real field')
    end if

    ! check if it's 0D-2D
    if (info % nDims > 2) then
      call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'the MPAS field "' // trim(config_result) // &
        '"assigned to variable ' // &
        trim(variable_names(var_index)) // ' in the RPN calculator AM is ' // &
        'not a 0D, 1D, or 2D field')
    end if
  end if

  ! increment the stack and put it on the stack
  stack_pointer = stack_pointer + 1
  stack(stack_pointer) % number_of_dims = info % nDims
  stack(stack_pointer) % symbol_type = IS_VARIABLE

  ! get the dimension name if it is 1D
  if (info % nDims == 0) then
    call mpas_pool_get_field(domain % blocklist % allFields, &
      config_result, stack(stack_pointer) % d0, 1)
  else if (info % nDims == 1) then
    call mpas_pool_get_field(domain % blocklist % allFields, &
      config_result, stack(stack_pointer) % d1, 1)
  else
    call mpas_pool_get_field(domain % blocklist % allFields, &
      config_result, stack(stack_pointer) % d2, 1)
  end if
end subroutine eval_variable!}}}



!***********************************************************************
! routine create_2d_field_from_1ds
!
!> \brief Generates a new 2D field from 1D fields
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details This will take two 1D fields (second and top) and
!> generate a new 2D field (head) with second's dimension as its
!> first dimension, and top's as it's second dimension. If top is
!> decomposed, head will be decomposed as well. If second has
!> constituent names, head will have constituent names as well.
!> Both fields need to be active, otherwise head will be inactive.
!-----------------------------------------------------------------------
subroutine create_2d_field_from_1ds(second, top_head, head)!{{{
#include "rpn_calc_inc/field_2d_from_1ds.inc"
end subroutine create_2d_field_from_1ds!}}}



!***********************************************************************
! routine create_1d_field_from_2d
!
!> \brief Generates a new 1D field from a 2D
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details This will take a 2D field (top_head) and
!> generate a new 1D field (head) with top's dimension as its
!> dimension. If top is decomposed, head will be decomposed as well.
!-----------------------------------------------------------------------
subroutine create_1d_field_from_2d(top_head, head)!{{{
#include "rpn_calc_inc/field_1d_from_2d.inc"
end subroutine create_1d_field_from_2d!}}}



!***********************************************************************
! routine create_0d_field_from_1d
!
!> \brief Generates a new 1D field from a 2D
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details This will take a 1D field (top) and
!> generate a new 0D field (dst).
!-----------------------------------------------------------------------
subroutine create_0d_field_from_1d(top, dst)!{{{
#include "rpn_calc_inc/field_0d_from_1d.inc"
end subroutine create_0d_field_from_1d!}}}



!***********************************************************************
! routine mul_operator
!
!> \brief Do mul on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a stack, take two arguments off the stack and
!> multiply them together, pushing the result back to the stack.
!-----------------------------------------------------------------------
subroutine mul_operator ( &
  exp_number, stack, stack_pointer, type_checking)!{{{
#include "rpn_calc_inc/binary_op_dispatch_start.inc"
  op_name = '*'
#include "rpn_calc_inc/binary_op_dispatch_0d_0d.inc"
  call mul_op_0d_0d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_0d_1d.inc"
  call mul_op_0d_1d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_0d_2d.inc"
  call mul_op_0d_2d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_0d.inc"
  call mul_op_1d_0d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_1d_same.inc"
  call mul_op_1d_1d_same(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_1d_diff.inc"
  call mul_op_1d_1d_diff(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_2d_first.inc"
  call mul_op_1d_2d_first(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_2d_second.inc"
  call mul_op_1d_2d_second(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_0d.inc"
  call mul_op_2d_0d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_1d_first.inc"
  call mul_op_2d_1d_first(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_1d_second.inc"
  call mul_op_2d_1d_second(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_2d.inc"
  call mul_op_2d_2d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_end.inc"
end subroutine mul_operator!}}}

subroutine mul_op_0d_0d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_0d_0d_1.inc"
  second * top
#include "rpn_calc_inc/binary_op_0d_0d_2.inc"
end subroutine mul_op_0d_0d!}}}

subroutine mul_op_0d_1d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_0d_1d_1.inc"
  second * top
#include "rpn_calc_inc/binary_op_0d_1d_2.inc"
end subroutine mul_op_0d_1d!}}}

subroutine mul_op_0d_2d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_0d_2d_1.inc"
  second * top
#include "rpn_calc_inc/binary_op_0d_2d_2.inc"
end subroutine mul_op_0d_2d!}}}

subroutine mul_op_1d_0d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_0d_1.inc"
  second * top
#include "rpn_calc_inc/binary_op_1d_0d_2.inc"
end subroutine mul_op_1d_0d!}}}

subroutine mul_op_1d_1d_same (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_1d_same_1.inc"
  second * top
#include "rpn_calc_inc/binary_op_1d_1d_same_2.inc"
end subroutine mul_op_1d_1d_same!}}}

subroutine mul_op_1d_1d_diff (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_1d_diff_1.inc"
  second(i) * top(j)
#include "rpn_calc_inc/binary_op_1d_1d_diff_2.inc"
end subroutine mul_op_1d_1d_diff!}}}

subroutine mul_op_1d_2d_first (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_2d_first_1.inc"
  second * top(:,j)
#include "rpn_calc_inc/binary_op_1d_2d_first_2.inc"
end subroutine mul_op_1d_2d_first!}}}

subroutine mul_op_1d_2d_second (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_2d_second_1.inc"
  second * top(i,:)
#include "rpn_calc_inc/binary_op_1d_2d_second_2.inc"
end subroutine mul_op_1d_2d_second!}}}

subroutine mul_op_2d_0d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_0d_1.inc"
  second * top
#include "rpn_calc_inc/binary_op_2d_0d_2.inc"
end subroutine mul_op_2d_0d!}}}

subroutine mul_op_2d_1d_first (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_1d_first_1.inc"
  second(:,j) * top
#include "rpn_calc_inc/binary_op_2d_1d_first_2.inc"
end subroutine mul_op_2d_1d_first!}}}

subroutine mul_op_2d_1d_second (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_1d_second_1.inc"
  second(i,:) * top
#include "rpn_calc_inc/binary_op_2d_1d_second_2.inc"
end subroutine mul_op_2d_1d_second!}}}

subroutine mul_op_2d_2d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_2d_1.inc"
  second * top
#include "rpn_calc_inc/binary_op_2d_2d_2.inc"
end subroutine mul_op_2d_2d!}}}



!***********************************************************************
! routine plus_operator
!
!> \brief Do plus on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a stack, take two arguments off the stack and
!> add them together, pushing the result back to the stack.
!-----------------------------------------------------------------------
subroutine plus_operator ( &
  exp_number, stack, stack_pointer, type_checking)!{{{
#include "rpn_calc_inc/binary_op_dispatch_start.inc"
  op_name = '+'
#include "rpn_calc_inc/binary_op_dispatch_0d_0d.inc"
  call plus_op_0d_0d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_0d_1d.inc"
  call plus_op_0d_1d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_0d_2d.inc"
  call plus_op_0d_2d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_0d.inc"
  call plus_op_1d_0d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_1d_same.inc"
  call plus_op_1d_1d_same(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_1d_diff.inc"
  call plus_op_1d_1d_diff(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_2d_first.inc"
  call plus_op_1d_2d_first(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_2d_second.inc"
  call plus_op_1d_2d_second(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_0d.inc"
  call plus_op_2d_0d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_1d_first.inc"
  call plus_op_2d_1d_first(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_1d_second.inc"
  call plus_op_2d_1d_second(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_2d.inc"
  call plus_op_2d_2d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_end.inc"
end subroutine plus_operator!}}}

subroutine plus_op_0d_0d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_0d_0d_1.inc"
  second + top
#include "rpn_calc_inc/binary_op_0d_0d_2.inc"
end subroutine plus_op_0d_0d!}}}

subroutine plus_op_0d_1d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_0d_1d_1.inc"
  second + top
#include "rpn_calc_inc/binary_op_0d_1d_2.inc"
end subroutine plus_op_0d_1d!}}}

subroutine plus_op_0d_2d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_0d_2d_1.inc"
  second + top
#include "rpn_calc_inc/binary_op_0d_2d_2.inc"
end subroutine plus_op_0d_2d!}}}

subroutine plus_op_1d_0d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_0d_1.inc"
  second + top
#include "rpn_calc_inc/binary_op_1d_0d_2.inc"
end subroutine plus_op_1d_0d!}}}

subroutine plus_op_1d_1d_same (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_1d_same_1.inc"
  second + top
#include "rpn_calc_inc/binary_op_1d_1d_same_2.inc"
end subroutine plus_op_1d_1d_same!}}}

subroutine plus_op_1d_1d_diff (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_1d_diff_1.inc"
  second(i) + top(j)
#include "rpn_calc_inc/binary_op_1d_1d_diff_2.inc"
end subroutine plus_op_1d_1d_diff!}}}

subroutine plus_op_1d_2d_first (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_2d_first_1.inc"
  second + top(:,j)
#include "rpn_calc_inc/binary_op_1d_2d_first_2.inc"
end subroutine plus_op_1d_2d_first!}}}

subroutine plus_op_1d_2d_second (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_2d_second_1.inc"
  second + top(i,:)
#include "rpn_calc_inc/binary_op_1d_2d_second_2.inc"
end subroutine plus_op_1d_2d_second!}}}

subroutine plus_op_2d_0d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_0d_1.inc"
  second + top
#include "rpn_calc_inc/binary_op_2d_0d_2.inc"
end subroutine plus_op_2d_0d!}}}

subroutine plus_op_2d_1d_first (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_1d_first_1.inc"
  second(:,j) + top
#include "rpn_calc_inc/binary_op_2d_1d_first_2.inc"
end subroutine plus_op_2d_1d_first!}}}

subroutine plus_op_2d_1d_second (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_1d_second_1.inc"
  second(i,:) + top
#include "rpn_calc_inc/binary_op_2d_1d_second_2.inc"
end subroutine plus_op_2d_1d_second!}}}

subroutine plus_op_2d_2d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_2d_1.inc"
  second + top
#include "rpn_calc_inc/binary_op_2d_2d_2.inc"
end subroutine plus_op_2d_2d!}}}



!***********************************************************************
! routine minus_operator
!
!> \brief Do minus on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a stack, take two arguments off the stack and
!> subtract top from second, pushing the result back to the stack.
!-----------------------------------------------------------------------
subroutine minus_operator ( &
  exp_number, stack, stack_pointer, type_checking)!{{{
#include "rpn_calc_inc/binary_op_dispatch_start.inc"
  op_name = '-'
#include "rpn_calc_inc/binary_op_dispatch_0d_0d.inc"
  call minus_op_0d_0d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_0d_1d.inc"
  if (type_checking) then
    call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
      'unable to subtract a 1d from a 0d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM')
  end if
#include "rpn_calc_inc/binary_op_dispatch_0d_2d.inc"
  if (type_checking) then
    call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
      'unable to subtract a 2d from a 0d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM')
  end if
#include "rpn_calc_inc/binary_op_dispatch_1d_0d.inc"
  call minus_op_1d_0d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_1d_same.inc"
  call minus_op_1d_1d_same(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_1d_diff.inc"
  call minus_op_1d_1d_diff(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_2d_first.inc"
  if (type_checking) then
    call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
      'unable to subtract a 2d from a 1d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM')
  end if
#include "rpn_calc_inc/binary_op_dispatch_1d_2d_second.inc"
  if (type_checking) then
    call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
      'unable to subtract a 2d from a 1d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM')
  end if
#include "rpn_calc_inc/binary_op_dispatch_2d_0d.inc"
  call minus_op_2d_0d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_1d_first.inc"
  call minus_op_2d_1d_first(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_1d_second.inc"
  call minus_op_2d_1d_second(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_2d.inc"
  call minus_op_2d_2d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_end.inc"
end subroutine minus_operator!}}}

subroutine minus_op_0d_0d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_0d_0d_1.inc"
  second - top
#include "rpn_calc_inc/binary_op_0d_0d_2.inc"
end subroutine minus_op_0d_0d!}}}

subroutine minus_op_1d_0d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_0d_1.inc"
  second - top
#include "rpn_calc_inc/binary_op_1d_0d_2.inc"
end subroutine minus_op_1d_0d!}}}

subroutine minus_op_1d_1d_same (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_1d_same_1.inc"
  second - top
#include "rpn_calc_inc/binary_op_1d_1d_same_2.inc"
end subroutine minus_op_1d_1d_same!}}}

subroutine minus_op_1d_1d_diff (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_1d_diff_1.inc"
  second(i) - top(j)
#include "rpn_calc_inc/binary_op_1d_1d_diff_2.inc"
end subroutine minus_op_1d_1d_diff!}}}

subroutine minus_op_2d_0d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_0d_1.inc"
  second - top
#include "rpn_calc_inc/binary_op_2d_0d_2.inc"
end subroutine minus_op_2d_0d!}}}

subroutine minus_op_2d_1d_first (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_1d_first_1.inc"
  second(:,j) - top
#include "rpn_calc_inc/binary_op_2d_1d_first_2.inc"
end subroutine minus_op_2d_1d_first!}}}

subroutine minus_op_2d_1d_second (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_1d_second_1.inc"
  second(i,:) - top
#include "rpn_calc_inc/binary_op_2d_1d_second_2.inc"
end subroutine minus_op_2d_1d_second!}}}

subroutine minus_op_2d_2d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_2d_1.inc"
  second - top
#include "rpn_calc_inc/binary_op_2d_2d_2.inc"
end subroutine minus_op_2d_2d!}}}



!***********************************************************************
! routine div_operator
!
!> \brief Do div on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a stack, take two arguments off the stack and
!> divide the second by the top, pushing the result back to the stack.
!-----------------------------------------------------------------------
subroutine div_operator ( &
  exp_number, stack, stack_pointer, type_checking)!{{{
#include "rpn_calc_inc/binary_op_dispatch_start.inc"
  op_name = '/'
#include "rpn_calc_inc/binary_op_dispatch_0d_0d.inc"
  call div_op_0d_0d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_0d_1d.inc"
  if (type_checking) then
    call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
      'unable to divide a 0d by a 1d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM')
  end if
#include "rpn_calc_inc/binary_op_dispatch_0d_2d.inc"
  if (type_checking) then
    call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
      'unable to divide a 0d by a 2d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM')
  end if
#include "rpn_calc_inc/binary_op_dispatch_1d_0d.inc"
  call div_op_1d_0d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_1d_same.inc"
  call div_op_1d_1d_same(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_1d_diff.inc"
  call div_op_1d_1d_diff(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_1d_2d_first.inc"
  if (type_checking) then
    call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
      'unable to divide a 1d by a 2d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM')
  end if
#include "rpn_calc_inc/binary_op_dispatch_1d_2d_second.inc"
  if (type_checking) then
    call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
      'unable to divide a 1d by a 2d in expression #' // &
      trim(expression_names(exp_number)) // ' in the RPN calculator AM')
  end if
#include "rpn_calc_inc/binary_op_dispatch_2d_0d.inc"
  call div_op_2d_0d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_1d_first.inc"
  call div_op_2d_1d_first(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_1d_second.inc"
  call div_op_2d_1d_second(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_2d_2d.inc"
  call div_op_2d_2d(stack, stack_pointer)
#include "rpn_calc_inc/binary_op_dispatch_end.inc"
end subroutine div_operator!}}}

function safe_divide_0d_0d(second, top)
  implicit none
  real (kind=RKIND), intent(in) :: second
  real (kind=RKIND), intent(in) :: top

  real (kind=RKIND) :: safe_divide_0d_0d

  if (abs(top) > 0.0_RKIND) then
    safe_divide_0d_0d = second / top
  else
    safe_divide_0d_0d = huge(second)
  end if
end function safe_divide_0d_0d

function safe_divide_1d_0d(second, top)
  implicit none
  real (kind=RKIND), dimension(:), intent(in) :: second
  real (kind=RKIND), intent(in) :: top

  real (kind=RKIND), dimension(size(second)) :: safe_divide_1d_0d

  if (abs(top) > 0.0_RKIND) then
    safe_divide_1d_0d = second / top
  else
    safe_divide_1d_0d = huge(second)
  end if
end function safe_divide_1d_0d

function safe_divide_2d_0d(second, top)
  implicit none
  real (kind=RKIND), dimension(:, :), intent(in) :: second
  real (kind=RKIND), intent(in) :: top

  real (kind=RKIND), dimension(size(second, 1), size(second, 2)) :: &
    safe_divide_2d_0d

  if (abs(top) > 0.0_RKIND) then
    safe_divide_2d_0d = second / top
  else
    safe_divide_2d_0d = huge(second)
  end if
end function safe_divide_2d_0d

function safe_divide_1d_1d(second, top)
  implicit none
  real (kind=RKIND), dimension(:), intent(in) :: second
  real (kind=RKIND), dimension(:), intent(in) :: top

  real (kind=RKIND), dimension(size(second)) :: safe_divide_1d_1d

  where (abs(top) > 0.0_RKIND)
    safe_divide_1d_1d = second / top
  elsewhere
    safe_divide_1d_1d = huge(second)
  end where
end function safe_divide_1d_1d

function safe_divide_2d_2d(second, top)
  implicit none
  real (kind=RKIND), dimension(:, :), intent(in) :: second
  real (kind=RKIND), dimension(:, :), intent(in) :: top

  real (kind=RKIND), dimension(size(second, 1), size(second, 2)) :: &
    safe_divide_2d_2d

  where (abs(top) > 0.0_RKIND)
    safe_divide_2d_2d = second / top
  elsewhere
    safe_divide_2d_2d = huge(second)
  end where
end function safe_divide_2d_2d

subroutine div_op_0d_0d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_0d_0d_1.inc"
  safe_divide_0d_0d(second, top)
#include "rpn_calc_inc/binary_op_0d_0d_2.inc"
end subroutine div_op_0d_0d!}}}

subroutine div_op_1d_0d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_0d_1.inc"
  safe_divide_1d_0d(second, top)
#include "rpn_calc_inc/binary_op_1d_0d_2.inc"
end subroutine div_op_1d_0d!}}}

subroutine div_op_1d_1d_same (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_1d_same_1.inc"
  safe_divide_1d_1d(second, top)
#include "rpn_calc_inc/binary_op_1d_1d_same_2.inc"
end subroutine div_op_1d_1d_same!}}}

subroutine div_op_1d_1d_diff (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_1d_1d_diff_1.inc"
  safe_divide_0d_0d(second(i), top(j))
#include "rpn_calc_inc/binary_op_1d_1d_diff_2.inc"
end subroutine div_op_1d_1d_diff!}}}

subroutine div_op_2d_0d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_0d_1.inc"
  safe_divide_2d_0d(second, top)
#include "rpn_calc_inc/binary_op_2d_0d_2.inc"
end subroutine div_op_2d_0d!}}}

subroutine div_op_2d_1d_first (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_1d_first_1.inc"
  safe_divide_1d_1d(second(:,j), top)
#include "rpn_calc_inc/binary_op_2d_1d_first_2.inc"
end subroutine div_op_2d_1d_first!}}}

subroutine div_op_2d_1d_second (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_1d_second_1.inc"
  safe_divide_1d_1d(second(i,:), top)
#include "rpn_calc_inc/binary_op_2d_1d_second_2.inc"
end subroutine div_op_2d_1d_second!}}}

subroutine div_op_2d_2d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/binary_op_2d_2d_1.inc"
  safe_divide_2d_2d(second, top)
#include "rpn_calc_inc/binary_op_2d_2d_2.inc"
end subroutine div_op_2d_2d!}}}

!***********************************************************************
! routine sum_operator
!
!> \brief Do sum on the stack
!> \author  Jon Woodring
!> \date    March 21, 2016
!> \details Given a stack, take sum argument off the stack and
!> sum along the first dimension, pushing the result back to the stack.
!-----------------------------------------------------------------------
subroutine sum_operator ( &
  exp_number, stack, stack_pointer, type_checking)!{{{
#include "rpn_calc_inc/reduce_op_dispatch_start.inc"
  op_name = 'sum'
#include "rpn_calc_inc/reduce_op_dispatch_1d.inc"
  call sum_op_1d(stack, stack_pointer)
#include "rpn_calc_inc/reduce_op_dispatch_2d.inc"
  call sum_op_2d(stack, stack_pointer)
#include "rpn_calc_inc/reduce_op_dispatch_end.inc"
end subroutine sum_operator

subroutine sum_op_1d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/reduce_op_1d_1.inc"
  0
#include "rpn_calc_inc/reduce_op_1d_2.inc"
  reduced + sum(top)
#include "rpn_calc_inc/reduce_op_1d_3.inc"
end subroutine sum_op_1d!}}}

subroutine sum_op_2d (stack, stack_pointer)!{{{
#include "rpn_calc_inc/reduce_op_2d_1.inc"
  0
#include "rpn_calc_inc/reduce_op_2d_2.inc"
  reduced(j) + sum(top(:,j))
#include "rpn_calc_inc/reduce_op_2d_3.inc"
end subroutine sum_op_2d!}}}

end module ocn_rpn_calculator
! vim: foldmethod=marker
