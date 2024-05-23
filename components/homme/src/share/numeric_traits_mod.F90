module numeric_traits_mod
  use kinds, only : real_kind

  implicit none
  public

#if HOMME_SINGLE_PRECISION
  real(kind=real_kind), parameter :: large_real = 1d30
#else
  real(kind=real_kind), parameter :: large_real = 1d99
#endif

end module numeric_traits_mod
