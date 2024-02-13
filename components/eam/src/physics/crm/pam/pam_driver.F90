module pam_driver_mod
  use iso_c_binding
  implicit none
  interface
    subroutine pam_driver() bind(C,name="pam_driver")
    end subroutine
    subroutine pam_finalize() bind(C,name="pam_finalize")
    end subroutine
  end interface
end module pam_driver_mod
