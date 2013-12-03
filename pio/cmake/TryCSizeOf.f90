       program trycsizeof
         use iso_c_binding, only : c_sizeof
         integer :: b
         character(len=2) :: a(5)
         b = c_sizeof(a(1))
      end program trycsizeof
