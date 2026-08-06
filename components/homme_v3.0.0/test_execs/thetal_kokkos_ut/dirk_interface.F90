module dirk_interface

  use kinds, only: real_kind

  implicit none

contains
  
  subroutine init_dirk_f90(ne, hyai, hybi, hyam, hybm, ps0) bind(c)
    use iso_c_binding,          only: c_int
    use hybvcoord_mod,          only: set_layer_locations
    use dimensions_mod,         only: nlev, nlevp, np
    use thetal_test_interface,  only: init_f90
    use edge_mod_base,          only: initEdgeBuffer, edge_g
    use geometry_interface_mod, only: par, elem

    real (kind=real_kind), intent(in) :: hyai(nlevp), hybi(nlevp), hyam(nlev), hybm(nlev)
    integer (kind=c_int), value, intent(in) :: ne
    real (kind=real_kind), value, intent(in) :: ps0

    real (kind=real_kind) :: mp(np,np), dvv(np,np)

    call init_f90(ne, hyai, hybi, hyam, hybm, dvv, mp, ps0)
    call initEdgeBuffer(par, edge_g, elem, 6*nlev+1)
  end subroutine init_dirk_f90

  subroutine pnh_and_exner_from_eos_f90(vtheta_dp,dp3d,dphi,pnh,exner,dpnh_dp_i) bind(c)
    use dimensions_mod,        only: nlev, nlevp, np
    use thetal_test_interface, only: hvcoord
    use eos,                   only: pnh_and_exner_from_eos2

    real (kind=real_kind), intent(in) :: vtheta_dp(np,np,nlev)   
    real (kind=real_kind), intent(in) :: dp3d(np,np,nlev)   
    real (kind=real_kind), intent(in) :: dphi(np,np,nlev)
    real (kind=real_kind), intent(out) :: pnh(np,np,nlev)        ! nh nonhyrdo pressure
    real (kind=real_kind), intent(out) :: dpnh_dp_i(np,np,nlevp) ! d(pnh) / d(pi)
    real (kind=real_kind), intent(out) :: exner(np,np,nlev)      ! exner nh pressure

    call pnh_and_exner_from_eos2(hvcoord, vtheta_dp, dp3d, dphi, pnh, exner, dpnh_dp_i, &
         'dirk_interface')
  end subroutine pnh_and_exner_from_eos_f90

  subroutine phi_from_eos_f90(phis,vtheta_dp,dp,phi_i) bind(c)
    use dimensions_mod,        only: nlev, nlevp, np
    use thetal_test_interface, only: hvcoord
    use eos,                   only: phi_from_eos

    real (kind=real_kind), intent(in) :: phis(np,np), vtheta_dp(np,np,nlev), dp(np,np,nlev)
    real (kind=real_kind), intent(out) :: phi_i(np,np,nlevp)

    call phi_from_eos(hvcoord, phis, vtheta_dp, dp, phi_i)
  end subroutine phi_from_eos_f90

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

  subroutine c2f_f90(n, cnlev, cnlevp, dp3d, w_i, v, vtheta_dp, phinh_i, gradphis, phis) bind(c)
    use iso_c_binding,         only: c_int
    use dimensions_mod,        only: nlev, nlevp, np
    use geometry_interface_mod,only: elem

    integer, parameter :: ntl = 3

    integer (kind=c_int), value, intent(in) :: n, cnlev, cnlevp
    real (kind=real_kind), intent(in) :: &
         dp3d(cnlev,np,np,ntl,n), w_i(cnlevp,np,np,ntl,n), v(cnlev,np,np,2,ntl,n), &
         vtheta_dp(cnlev,np,np,ntl,n), phinh_i(cnlevp,np,np,ntl,n), gradphis(np,np,2,n), &
         phis(np,np,n)

    integer :: ie, t, k, j, i

    if (n /= size(elem)) print *, 'n /= size(nelem)'
    if (cnlev < nlev) print *, 'nlev < cnlev'
    if (cnlevp < nlevp) print *, 'nlevp < cnlevp'

    do ie = 1,n
       do t = 1,ntl
          do k = 1,nlevp
             do j = 1,np
                do i = 1,np
                   elem(ie)%state%w_i(i,j,k,t) = w_i(k,j,i,t,ie)
                   elem(ie)%state%phinh_i(i,j,k,t) = phinh_i(k,j,i,t,ie)
                   if (k /= nlevp) then
                      elem(ie)%state%dp3d(i,j,k,t) = dp3d(k,j,i,t,ie)
                      elem(ie)%state%v(i,j,1,k,t) = v(k,j,i,1,t,ie)
                      elem(ie)%state%v(i,j,2,k,t) = v(k,j,i,2,t,ie)
                      elem(ie)%state%vtheta_dp(i,j,k,t) = vtheta_dp(k,j,i,t,ie)
                   end if
                end do
             end do
          end do
       end do
       do j = 1,np
          do i = 1,np
             elem(ie)%state%phis(i,j) = phis(j,i,ie)
             elem(ie)%derived%gradphis(i,j,1) = gradphis(j,i,1,ie)
             elem(ie)%derived%gradphis(i,j,2) = gradphis(j,i,2,ie)
          end do
       end do
    end do
  end subroutine c2f_f90

  subroutine f2c_f90(n, cnlev, cnlevp, dp3d, w_i, v, vtheta_dp, phinh_i, gradphis) bind(c)
    use iso_c_binding,         only: c_int
    use dimensions_mod,        only: nlev, nlevp, np
    use geometry_interface_mod,only: elem

    integer, parameter :: ntl = 3

    integer (kind=c_int), value, intent(in) :: n, cnlev, cnlevp
    real (kind=real_kind), intent(out) :: &
         dp3d(cnlev,np,np,ntl,n), w_i(cnlevp,np,np,ntl,n), v(cnlev,np,np,2,ntl,n), &
         vtheta_dp(cnlev,np,np,ntl,n), phinh_i(cnlevp,np,np,ntl,n), gradphis(np,np,2,n)

    integer :: ie, t, k, j, i

    if (n /= size(elem)) print *, 'n /= size(nelem)'
    if (cnlev < nlev) print *, 'nlev < cnlev'
    if (cnlevp < nlevp) print *, 'nlevp < cnlevp'

    do ie = 1,n
       do t = 1,ntl
          do k = 1,nlevp
             do j = 1,np
                do i = 1,np
                   w_i(k,j,i,t,ie) = elem(ie)%state%w_i(i,j,k,t)
                   phinh_i(k,j,i,t,ie) = elem(ie)%state%phinh_i(i,j,k,t)
                   if (k /= nlevp) then
                      dp3d(k,j,i,t,ie) = elem(ie)%state%dp3d(i,j,k,t)
                      v(k,j,i,1,t,ie) = elem(ie)%state%v(i,j,1,k,t)
                      v(k,j,i,2,t,ie) = elem(ie)%state%v(i,j,2,k,t)
                      vtheta_dp(k,j,i,t,ie) = elem(ie)%state%vtheta_dp(i,j,k,t)
                   end if
                end do
             end do
          end do
       end do
       do j = 1,np
          do i = 1,np
             gradphis(j,i,1,ie) = elem(ie)%derived%gradphis(i,j,1)
             gradphis(j,i,2,ie) = elem(ie)%derived%gradphis(i,j,2)
          end do
       end do
    end do
  end subroutine f2c_f90

  subroutine compute_stage_value_dirk_f90(nm1,alphadt_nm1,n0,alphadt_n0,np1,dt2) bind(c)
    use iso_c_binding,         only: c_int
    use hybrid_mod,            only: hybrid_t
    use derivative_mod,        only: derivative_t
    use geometry_interface_mod,only: elem
    use thetal_test_interface, only: hvcoord
    use imex_mod,              only: compute_stage_value_dirk

    integer (kind=c_int), value, intent(in)  :: nm1, n0, np1
    real (kind=real_kind), value, intent(in) :: alphadt_nm1, alphadt_n0, dt2

    ! compute_stage_value_dirk doesn't actually use these.
    type (hybrid_t) :: hybrid
    type (derivative_t) :: deriv
    integer :: itercount, nets, nete
    real (kind=real_kind) :: itererr

    nets = 1
    nete = size(elem)
    call compute_stage_value_dirk(nm1, alphadt_nm1, n0, alphadt_n0, np1, dt2, &
         elem, hvcoord, hybrid, deriv, nets, nete, itercount, itererr, &
         verbosity_in=0)
  end subroutine compute_stage_value_dirk_f90

end module dirk_interface
