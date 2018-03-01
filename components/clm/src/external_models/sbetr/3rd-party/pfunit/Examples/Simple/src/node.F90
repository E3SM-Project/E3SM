
module node_mod
  implicit none
  private

  public :: node

  type node
     real, pointer :: xPtr
   contains
     procedure :: getPtr
     procedure :: getX
  end type node

  interface node
     module procedure new_node
  end interface node

contains

  function new_node(ptr) result(newNode)
    type(node) :: newNode
    real, pointer, intent(in) :: ptr
    newNode%xPtr => ptr
  end function new_node

  function getPtr(this) result(ptr)
    class(node), intent(inout) :: this
    real, pointer :: ptr
    ptr => this%xPtr
  end function getPtr

  function getX(this) result(x)
    class(node), intent(inout) :: this
    real :: x
    x = this%xPtr
  end function getX
    
end module node_mod
