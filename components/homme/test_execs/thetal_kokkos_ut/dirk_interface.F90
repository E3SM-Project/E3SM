module dirk_interface

  use kinds, only: real_kind

  implicit none

contains
  
  subroutine init_dirk_f90(hyai, hybi, hyam, hybm, ps0) bind(c)
    use iso_c_binding,          only: c_int
    use hybvcoord_mod,          only: set_layer_locations
    use dimensions_mod,         only: nlev, nlevp, np
    use thetal_test_interface,  only: hvcoord

    real (kind=real_kind), intent(in)  :: hyai(nlevp), hybi(nlevp), hyam(nlev), hybm(nlev)
    real (kind=real_kind), intent(in)  :: ps0

    integer :: k

    hvcoord%hyai = hyai
    hvcoord%hybi = hybi
    hvcoord%hyam = hyam
    hvcoord%hybm = hybm
    hvcoord%ps0 = ps0
    do k = 1,nlev
      hvcoord%dp0(k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*ps0 + &
                       (hvcoord%hybi(k+1) - hvcoord%hybi(k))*ps0
    enddo
    call set_layer_locations (hvcoord,.false.,.false.)
  end subroutine init_dirk_f90

  subroutine compute_gwphis_f90(gwh_i,dp3d,v,gradphis) bind(c)
    use dimensions_mod,        only: nlev, nlevp, np
    use thetal_test_interface, only: hvcoord
    use imex_mod,              only: compute_gwphis

    real (kind=real_kind) :: gwh_i(np,np,nlevp)
    real (kind=real_kind) :: dp3d(np,np,nlev)
    real (kind=real_kind) :: v(np,np,2,nlev)
    real (kind=real_kind) :: gradphis(np,np,2)

    call compute_gwphis(gwh_i,dp3d,v,gradphis,hvcoord)
  end subroutine compute_gwphis_f90

  subroutine get_dirk_jacobian_f90(JacL,JacD,JacU,dt2,dp3d,dphi,pnh) bind(c)
    use dimensions_mod,        only: nlev, nlevp, np
    use imex_mod,              only: get_dirk_jacobian

    real (kind=real_kind), intent(in)    :: dp3d(np,np,nlev)
    real (kind=real_kind), intent(in)    :: dphi(np,np,nlev)
    real (kind=real_kind), value, intent(in) :: dt2
    real (kind=real_kind), intent(inout) :: pnh(np,np,nlev)
    real (kind=real_kind), intent(out)   :: JacD(nlev,np,np)
    real (kind=real_kind), intent(out)   :: JacL(nlev-1,np,np),JacU(nlev-1,np,np)

    call get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,dphi,pnh,1)
  end subroutine get_dirk_jacobian_f90

end module dirk_interface
