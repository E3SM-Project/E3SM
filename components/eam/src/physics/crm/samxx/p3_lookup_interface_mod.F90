module p3_lookup_interface_mod
  use iso_c_binding
  implicit none
  public :: initialize_p3_lookup
  interface
    
    subroutine initialize_p3_lookup() bind(C,name="initialize_p3_lookup")
    end subroutine initialize_p3_lookup

    subroutine finalize_p3_lookup() bind(C,name="finalize_p3_lookup")
    end subroutine finalize_p3_lookup

  end interface
end module p3_lookup_interface_mod
