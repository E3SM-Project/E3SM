module mpas_stack

   implicit none

   private

   ! Public Subroutines and Structures
   public :: mpas_stack_is_empty
   public :: mpas_stack_push
   public :: mpas_stack_pop
   public :: mpas_stack_free

   public :: mpas_stack_type, mpas_stack_payload_type

   type mpas_stack_payload_type
   end type mpas_stack_payload_type

   type mpas_stack_type
      type (mpas_stack_type), pointer :: next => null()
      class (mpas_stack_payload_type), pointer :: payload => null()
   end type mpas_stack_type

   !***********************************************************************
   !
   !  module mpas_stack
   !
   !> \brief   MPAS Stack module
   !> \author  Miles A. Curry
   !> \date    04/04/19
   !> \details
   !>
   !> Introduction
   !> ==============
   !> The MPAS stack is a simple, extensible data stack data structure for use
   !> within the MPAS atmospheric model. It functions as a wrapper around a
   !> polymorphic data structure to provide usage in different areas.
   !>
   !>
   !> Creating a Stack
   !> ==================
   !> The stack data structure (`type (mpas_stack_type)`) is defined by a single 
   !> `next` pointer > and a pointer to a `type (mpas_stack_payload_type)`, which 
   !> is defined as a empty derived type.
   !>
   !> To use the stack, create a derived type that extends the `mpas_stack_payload_type` 
   !> type. Define your extended derived type with members that meets your application.
   !>
   !> For instance:
   !> ```
   !> type, extends(mpas_stack_payload_type) :: my_payload_name
   !>    ! Define the members of your type as you wish
   !> end type my_payload_name
   !>
   !> class (my_payload_name), pointer :: item1 => null(), item2 => null()
   !> ```
   !>
   !> The extended mpas_stack_payload_type will enable a user defined type to be 
   !> associated with a stack item. The stack stores references of a payload, thus 
   !> a single payload can be used in multiple push operations.
   !>
   !> You will then need to create a stack (or multiple stacks if you desire) as
   !> the following:
   !>
   !> ```
   !> type (mpas_stack_type), pointer :: stack1 => null(), stack2 => null()
   !> ```
   !>
   !>  Pushing onto a Stack
   !>  ====================
   !>  You can push your items onto a stack as:
   !>
   !> ```
   !> allocate(item1)
   !> stack1 => mpas_stack_push(stack1, item1)
   !> allocate(item2)
   !> stack1 => mpas_stack_push(stack1, item2)
   !> ```
   !>
   !> Popping an item off of the stack
   !> ================================
   !> Popping an item off of the stack will require a bit more work than pushing.
   !> Because the payload is a polymorphic class , we will need to use the select 
   !> case to get our type (or multiple types) back into a usable object:
   !> ```
   !> ! The item to pop items into
   !> class (mpas_stack_payload_type), pointer :: top
   !> type (my_payload_name), pointer :: my_item
   !>
   !> top => mpas_stack_pop(stack1)
   !> select type(top)
   !>    type is(my_payload_name)
   !>       my_item => top
   !> end  select
   !> ```
   !>
   !> Note: It is recommended to create your own `pop` function so you can reduce 
   !> the amount of coded needed. An example is provided at the bottom of
   !> this module as the function `user_pop(..)`
   !
   !-----------------------------------------------------------------------

   contains

   !***********************************************************************
   !
   !  routine mpas_stack_is_empty
   !
   !> \brief   Returns .true. if the stack is empty, otherwise .false.
   !> \author  Miles A. Curry
   !> \date    01/28/20
   !> Returns .true. If the stack is empty and/or if the stack is unassociated.
   !
   !-----------------------------------------------------------------------
   function mpas_stack_is_empty(stack) result(is_empty)

      implicit none
      type (mpas_stack_type), intent(in), pointer :: stack
      logical :: is_empty

      is_empty = .true.
      if (associated(stack)) then
         is_empty = .false.
         return
      endif

   end function mpas_stack_is_empty

   !***********************************************************************
   !
   !  routine mpas_stack_push
   !
   !> \brief   Push an item onto stack
   !> \author  Miles A. Curry
   !> \date    01/28/20
   !> \details
   !>
   !> Push a mpas_stack_payload_type type, onto `stack` and return the new stack. If
   !> `payload` is the first item to be pushed onto the stack, then `stack`
   !> should be unassociated.
   !
   !-----------------------------------------------------------------------
   function mpas_stack_push(stack, payload) result(new_stack)
      
      implicit none

      type(mpas_stack_type), intent(inout), pointer :: stack
      class(mpas_stack_payload_type), intent(inout), target :: payload

      type(mpas_stack_type), pointer :: new_stack

      allocate(new_stack)
      new_stack % payload => payload
      new_stack % next => stack

      return

   end function mpas_stack_push

   !***********************************************************************
   !
   !  function mpas_stack_pop
   !
   !> \brief   Pop off the last item added from a stack
   !> \author  Miles A. Curry
   !> \date    01/28/20
   !> \details
   !> Pop off and return the top item of the stack as a `class mpas_stack_payload_type`.
   !> If the stack is empty (or unassociated), then a null `class mpas_stack_payload_type`
   !> pointer will be returned. `select type` will need to be used to retrieve
   !> any extended members.
   !
   !-----------------------------------------------------------------------
   function mpas_stack_pop(stack) result(top)

      implicit none

      type (mpas_stack_type), intent(inout), pointer :: stack
      type (mpas_stack_type), pointer :: next => null()
      class(mpas_stack_payload_type), pointer :: top

      if ( .not. associated(stack)) then
         top => null()
         return
      endif

      top => stack % payload
      next => stack % next
      deallocate(stack)
      stack => next
      return

   end function mpas_stack_pop

   !***********************************************************************
   !
   !  function mpas_stack_free
   !
   !> \brief   Deallocate the entire stack. Optionally deallocate payloads
   !> \author  Miles A. Curry 
   !> \date    01/28/20
   !> \details
   !>  Deallocate the entire stack. If free_payload is set to `.true.` or if
   !>  absent then the payload will be deallocated. If not, then the payload will not
   !>  be deallocated. Upon success, the stack will be unassociated.
   !  
   !-----------------------------------------------------------------------
   subroutine mpas_stack_free(stack, free_payload)

      implicit none

      type(mpas_stack_type), intent(inout), pointer :: stack
      logical, intent(in), optional :: free_payload
      logical :: fpl

      type(mpas_stack_type), pointer :: cur

      if (present(free_payload)) then
         fpl = free_payload
      else
         fpl = .true.
      endif

      cur => stack
      do while(associated(stack))
         stack => stack % next
         if ( fpl ) then
            deallocate(cur % payload)
         endif
         deallocate(cur)
         cur => stack
      enddo

   end subroutine mpas_stack_free


   !***********************************************************************
   !
   !  Example user-defined pop function
   !
   !> \brief   Pop off the last item added from a stack and return it as our
   !>          defined type
   !> \author  Miles A. Curry
   !> \date    01/28/20
   !
   !-----------------------------------------------------------------------
   ! function user_pop(stack) result(item)
   !
   !    use mpas_stack, only : mpas_stack_type, mpas_stack_payload_type, mpas_stack_pop
   !
   !    implicit none
   !
   !    type(mpas_stack_type), intent(inout), pointer :: stack
   !
   !    type(my_item), pointer :: item    ! Our user defined mpas_stack_type
   !
   !    ! We will need to use the mpas_stack_payload_type type to use mpas_stack_pop(...)
   !    class(mpas_stack_payload_type), pointer :: top
   !
   !    !
   !    ! Handle a pop on an empty stack if we want to here
   !    ! Note the stack will return null if it is empty.
   !    !
   !    if (mpas_stack_is_empty(stack)) then
   !       item => null()
   !       return
   !    endif
   ! 
   !    top => mpas_stack_pop(stack)
   !    
   !    select type(top)
   !       type is(my_item)
   !          item => top
   !       class default
   !          write(0,*) "We got an Error and we should handle it if we need to!!"
   !          stop
   !    end select
   !
   ! end function user_pop

end module mpas_stack
