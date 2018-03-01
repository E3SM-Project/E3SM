module adgaquad_types_mod
  use carma_precision_mod

  integer, public, parameter      :: nf=50  !! Number of factorials in fact table.


  !! The the functions that are being integrated may need some extra
  !! data. In the F77, these were stored in common blocks. To make the
  !! code thread safe, we need to move them into passed parameters. For
  !! convenience, we put all of these variables into one structure and
  !! pass the entire structure to all functions that could be integrated
  !! by these routines.

  type, public :: adgaquad_vars_type
  
    !   alpha         packing coefficient
    !   nb            number of monomers
    !   a             monomer size
    !   df            fractal dimension
    !   k             absolute value of wavevector = 2*pi/wavelength

    real(kind=f)                                :: fact(0:nf)
    integer                                     :: u1
    integer                                     :: u2
    integer                                     :: u3
    integer                                     :: u4
    integer                                     :: u5
    integer                                     :: u6
    integer                                     :: pbes
    real(kind=f)                                :: kbes
    real(kind=f)                                :: alpha
    real(kind=f)                                :: nb
    real(kind=f)                                :: a
    real(kind=f)                                :: df
    real(kind=f)                                :: k
    real(kind=f)                                :: zed
    real(kind=f)                                :: coeff
  end type adgaquad_vars_type
end module