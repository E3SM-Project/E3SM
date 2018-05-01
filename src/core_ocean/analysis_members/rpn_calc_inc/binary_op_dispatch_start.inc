  integer, intent(in) :: exp_number
  logical, intent(in) :: type_checking

  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  integer :: a_nd, b_nd
  character (len=StrKIND) ::  op_name

! 0d 0d
