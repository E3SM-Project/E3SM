module da_ml_wrapper
  interface
  integer function ml_init() bind(C, name="ml_init_c")
    use iso_c_binding
    implicit none
  end function
  end interface
end module da_ml_wrapper
 
