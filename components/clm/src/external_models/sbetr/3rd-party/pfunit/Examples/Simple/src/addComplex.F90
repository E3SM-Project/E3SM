
module addComplex_mod

  implicit none

contains

  function add(z0,z1) result (z2)
    complex, intent(in), dimension(:) :: z0,z1
    complex, dimension(size(z0)) :: z2
    ! Assume size(z0) = size(z1) for now...
    z2 = z0 + z1
  end function add

end module addComplex_mod
