module GetGlobalValuesMod

  ! Mock of GetGlobalValuesMod, which satisfies routine signatures with minimal
  ! dependencies

  implicit none
  
  public :: GetGlobalWrite

contains

  subroutine GetGlobalWrite(decomp_index, clmlevel)
    integer, intent(in) :: decomp_index
    character(len=*), intent(in) :: clmlevel

    ! do nothing
  end subroutine GetGlobalWrite

end module GetGlobalValuesMod
