module test_mod
   use iso_fortran_env, only : real64
  implicit none
  private

  integer, parameter  :: r8 = real64

  type, public :: my_type
     real(r8), pointer :: member1(:) => null()
     real(r8), pointer :: member2(:) => null()
  end type my_type

  type(my_type), public, target :: my_inst
contains

end module test_mod

program test_factory_registry

end program test_factory_registry

