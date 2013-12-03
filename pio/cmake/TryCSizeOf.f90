       program trycsizeof
         use iso_c_binding, only : c_sizeof
         integer :: a,b
         b = c_sizeof(a)
      end program trycsizeof
