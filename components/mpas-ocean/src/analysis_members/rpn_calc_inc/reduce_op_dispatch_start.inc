  integer, intent(in) :: exp_number
  logical, intent(in) :: type_checking

  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  integer :: a_nd
  character (len=StrKIND) ::  op_name

! start -> 1d
