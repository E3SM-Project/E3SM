program main

   implicit none
   real :: a
   character(len=8) :: char 

   namelist /my_nl/ a, char

   write(*,*) '--------------------'
   write(*,*) 'before namelist read'
   write(*,*) '--------------------'
   write(*,nml=my_nl)
   write(*,*)

   read(5,nml=my_nl) 

   write(*,*) '--------------------'
   write(*,*) 'after namelist read'
   write(*,*) '--------------------'
   write(*,nml=my_nl)

end program
