program homme_nvector_driver

  use HommeNVector,   only: NVec_t, MakeHommeNVector
  use element_mod,    only: element_t
  use hybvcoord_mod,  only: hvcoord_t
  use hybrid_mod,     only: hybrid_t
  use derivative_mod, only: derivative_t
  use dimensions_mod, only: np, nlev
  use, intrinsic :: iso_c_binding

  implicit none
  type (element_t), allocatable :: elemX(:), elemY(:), elemZ(:)
  type (hvcoord_t)              :: hvcoord
  type (hybrid_t)               :: hybrid
  type (derivative_t)           :: deriv
  type(NVec_t), target          :: x, y, z
  type(c_ptr)                   :: x_C, y_C, z_C
  integer                       :: nets, nete, qn0, tl_idx, ier

  integer                       :: numElem

  !=======Internals ============

  ! Set number of testing elements and time level index
  numElem = 1
  nets = 1
  nete = 1

  ! Allocate elem object and initialize state variables for vector x
  allocate(elemX(numElem))
  tl_idx = 1
  elemX(1)%state%v(1:np,1:np,1:2,1:nlev,tl_idx) = 1.d0
  elemX(1)%state%w(1:np,1:np,1:nlev,tl_idx) = 1.d0
  elemX(1)%state%phi(1:np,1:np,1:nlev,tl_idx) = 1.d0
  elemX(1)%state%theta_dp_cp(1:np,1:np,1:nlev,tl_idx) = 1.d0
  elemX(1)%state%dp3d(1:np,1:np,1:nlev,tl_idx) = 1.d0

  ! Createt HommeNVector object and corresponding c pointer
  call MakeHommeNVector(elemX, hvcoord, hybrid, deriv, nets, nete, qn0, tl_idx, x, ier)
  if (ier == 1) then
    stop "Error making Homme NVector for x"
  end if
  x_C = c_loc(x)

  ! Allocate elem object and initialize state variables for vector y
  allocate(elemY(numElem))
  tl_idx = 2
  elemY(1)%state%v(1:np,1:np,1:2,1:nlev,tl_idx) = 2.d0
  elemY(1)%state%w(1:np,1:np,1:nlev,tl_idx) = 2.d0
  elemY(1)%state%phi(1:np,1:np,1:nlev,tl_idx) = 2.d0
  elemY(1)%state%theta_dp_cp(1:np,1:np,1:nlev,tl_idx) = 2.d0
  elemY(1)%state%dp3d(1:np,1:np,1:nlev,tl_idx) = 2.d0

  ! Create HommeNVector object and corresponding c pointer
  call MakeHommeNVector(elemY, hvcoord, hybrid, deriv, nets, nete, qn0, tl_idx, y, ier)
  if (ier == 1) then
    stop "Error making Homme NVector for y"
  end if
  y_C = c_loc(y)

  ! Allocate elem object for vector z
  allocate(elemZ(numElem))
  tl_idx = 3

  ! Create HommeNVector object and corresponding c pointer
  call MakeHommeNVector(elemY, hvcoord, hybrid, deriv, nets, nete, qn0, tl_idx, z, ier)
  if (ier == 1) then
    stop "Error making Homme NVector for z"
  end if
  z_C = c_loc(z)


  ! Call function to test and print out result
  call FNVExtLinearSum(1.d0, x_C, 1.d0, y_C, z_C)
  call FNVExtPrint(z_C)

end program
