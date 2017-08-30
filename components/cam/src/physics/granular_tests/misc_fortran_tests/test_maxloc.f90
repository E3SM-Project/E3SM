program main

  real :: array_1d(3)
  real :: array_2d(3,2)
  integer :: icol_1d(1)
  integer :: icol_2d(2)

  array_1d = (/1.,2.,3./)
  array_2d(:,1) = (/1.,2.,3./)
  array_2d(:,2) = (/4.,5.,6./)

  icol_1d = maxloc( abs(array_1d) )
  icol_2d = maxloc( abs(array_2d) )

  print*,  icol_1d, array_1d(icol_1d(1))
  print*,  icol_2d, array_2d(icol_2d(1),icol_2d(2))
end
