#ifndef CAM
#include "config.h"

module planar_transport_tests
  use kinds, only: rl=>real_kind, iulog
  use element_mod, only: element_t
  use parallel_mod, only: parallel_t, abortmp
  use hybrid_mod, only: hybrid_t
  use hybvcoord_mod, only: hvcoord_t, set_layer_locations
  use derivative_mod, only: derivative_t, gradient_sphere
  use element_ops, only: set_state, set_state_i
  ! Planar geometry parameters.
  use physical_constants, only: Lx, Ly, dd_pi, &
       &                        Rgas, g, cp, pi => dd_pi, p0
  use dimensions_mod, only: ne_x, ne_y, qsize, qsize_d, nlev, nlevp, np, nelemd
  ! Test problem tools.
  use dcmip12_wrapper, only: get_evenly_spaced_z, set_hybrid_coefficients, &
       &                     pressure_thickness

  implicit none
  private

  public :: test_conv_planar_advection, print_conv_planar_advection_results

  real(rl), parameter :: &
       tau     = 12.d0 * 86400.d0, & ! period of motion 12 days
       T0      = 300.d0,           & ! temperature (K)
       ztop    = 12000.d0,         & ! model top (m)
       H       = Rgas * T0 / g       ! scale height

  character :: tc_major, tc_minor
  real(rl) :: zi(nlevp), zm(nlev)

contains  

  subroutine test_conv_planar_advection( &
       test_case, elem, hybrid, hvcoord, deriv, nets, nete, time, n0, n1)
    character(len=*), intent(in):: test_case
    type (element_t), intent(inout), target :: elem(:)
    type (hybrid_t), intent(in):: hybrid
    type (hvcoord_t), intent(inout) :: hvcoord
    type (derivative_t), intent(in):: deriv
    integer, intent(in):: nets, nete, n0, n1
    real(rl), intent(in):: time

    integer :: ie, k, j, i, qi
    real(rl) :: x, y, u, v, w, T, ps, phis, p, dp, z, q(qsize)
    real(rl):: grad_p(np,np,2), p_i(np,np), u_i(np,np), v_i(np,np)

    if (time <= 0.d0) then
       call init(test_case, hybrid, hvcoord)
    end if

    w = 0
    do ie = nets,nete
       do k = 1,nlevp
          do j = 1,np
             do i = 1,np
                x = elem(ie)%spherep(i,j)%lon
                y = elem(ie)%spherep(i,j)%lat
                if (k < nlevp) then
                   call get_values(time, x, y, hvcoord%hyam(k), hvcoord%hybm(k), &
                        &          ps, phis, p, z, T, u, v, q)
                   dp = pressure_thickness(ps,k,hvcoord)
                   if (time <= 0.d0) then
                      do qi = 1, qsize
                         elem(ie)%state%Q(i,j,k,qi) = q(qi)
                         elem(ie)%state%Qdp(i,j,k,qi,:) = elem(ie)%state%Q(i,j,k,qi) * dp
                      end do
                   end if
                   call set_state(u, v, w, T, ps, phis, p, dp, z, g, i, j, k, elem(ie), n0, n1)
                end if
                call get_values(time, x, y, hvcoord%hyai(k), hvcoord%hybi(k), &
                     &          ps, phis, p, z, T, u, v, q)
                call set_state_i(u, v, w, T, ps, phis, p, z, g, i, j, k, elem(ie), n0, n1)
                p_i(i,j) = p
                u_i(i,j) = u
                v_i(i,j) = v
             end do
          end do
          ! Get vertical mass flux. This is not used in the case of interest: 3D
          ! transport. But include it for completeness.
          grad_p = gradient_sphere(p_i,deriv,elem(ie)%Dinv)
          elem(ie)%derived%eta_dot_dpdn_prescribed(:,:,k) = &
               -u_i*grad_p(:,:,1) - v_i*grad_p(:,:,2)
       end do
       elem(ie)%derived%eta_dot_dpdn_prescribed(:,:,1)     = 0
       elem(ie)%derived%eta_dot_dpdn_prescribed(:,:,nlevp) = 0       
    end do
  end subroutine test_conv_planar_advection

  subroutine get_values(time, x, y, hya, hyb, &
       &                ps, phis, p, z, T, u, v, q)
    real(rl), intent(in) :: time, x, y, hya, hyb
    real(rl), intent(out) :: ps, phis, p, z, T, u, v, q(qsize)

    real(rl), parameter :: &
         xm = 0.25, &
         distm = 3.d0/8.d0, &
         h0 = 2000.d0, &
         ztop_t = 2000.d0, &
         z1_h = ztop_t + 1000.d0, &
         z2_h = ztop_t + 6000.d0, &
         z0_h = (z1_h + z2_h)/2

    real(rl) :: u_topo_fac, xm_t, dist, zs, ztaper

    zs = 0
    xm_t = xm*Lx
    ! The mountain center moves in time.
    u_topo_fac = -Lx/tau/2.d0
    xm_t = xm_t + sin(pi*time/tau)*(tau/pi)*u_topo_fac
    ! Mountain shape.
    dist = min(min(abs(x - xm_t), &
         &         abs(x - (xm_t - Lx))), &
         &         abs(x - (xm_t + Lx)))
    dist = dist/Lx
    if (dist < distm) then
       zs = h0*(1.d0 + cos(pi*(dist/distm)))*cos(pi*(6.d0*dist))**2.d0
    end if

    zs = -zs
    phis = g*zs
    ps = p0 * exp(-zs/H)
    p = hya*p0 + hyb*ps
    z = H * log(p0/p)

    if (z <= 0) then
       ztaper = 0
    elseif (z >= ztop_t) then
       ztaper = 1
    else
       ztaper = (1 + cos(pi*(1 + z/ztop_t)))/2
    end if

    ! Simple translation in x direction.
    u = (Lx/tau)*ztaper
    ! Account for moving ps.
    u = u + cos(pi*time/tau)*u_topo_fac*(1 - ztaper)

    v = 0
    T = T0

    if (time > 0.d0) return

    if (qsize == 0) return
    q = 0
    if (z < z2_h .and. z > z1_h) then
       q(1) = 0.5d0 * (1 + cos(2.d0*pi*(z-z0_h)/(z2_h-z1_h)))
    end if
    if (qsize > 2) q(3) = q(1)
    q(1) = q(1)*(1 + cos(6.d0*pi*(x/Lx)))
    q(2) = q(1)
    if (qsize > 3) q(3:qsize) = q(1)
  end subroutine get_values
  
  subroutine print_conv_planar_advection_results(test_case, elem, tl, hvcoord, par)
    use time_mod, only: timelevel_t
    use parallel_mod, only: global_shared_buf, global_shared_sum, pmax_1d
    use global_norms_mod, only: wrap_repro_sum

    character(len=*), intent(in) :: test_case
    type(element_t), intent(in) :: elem(:)
    type(timelevel_t), intent(in) :: tl
    type(hvcoord_t), intent(in) :: hvcoord
    type(parallel_t), intent(in) :: par

    integer :: ie, k, j, i, iq
    real(rl) :: time, x, y, ps, phis, p, z, T, u, v, q(np,np,qsize), &
         &      reldif, linf_num(qsize), linf_den(qsize), a, b
    
    ! Set time to 0 to get the initial conditions.
    time = 0

    linf_num = 0
    linf_den = 0
    do ie = 1,nelemd
       global_shared_buf(ie,:2*qsize) = 0
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                x = elem(ie)%spherep(i,j)%lon
                y = elem(ie)%spherep(i,j)%lat
                call get_values(time, x, y, hvcoord%hyam(k), hvcoord%hybm(k), &
                     &          ps, phis, p, z, T, u, v, q(i,j,:))
             end do
          end do
          do iq = 1,qsize
             global_shared_buf(ie,2*iq-1) = global_shared_buf(ie,2*iq-1) + &
                  sum(elem(ie)%spheremp*(elem(ie)%state%Q(:,:,k,iq) - q(:,:,iq))**2)
             global_shared_buf(ie,2*iq) = global_shared_buf(ie,2*iq) + &
                  sum(elem(ie)%spheremp*q(:,:,iq)**2)
             linf_num(iq) = max(linf_num(iq), &
                  maxval(abs(elem(ie)%state%Q(:,:,k,iq) - q(:,:,iq))))
             linf_den(iq) = max(linf_den(iq), &
                  maxval(abs(q(:,:,iq))))
          end do
       end do
    end do

    call wrap_repro_sum(nvars=2*qsize, comm=par%comm)
    do iq = 1, qsize
       linf_num(iq) = pmax_1d(linf_num(iq:iq), par)
       linf_den(iq) = pmax_1d(linf_den(iq:iq), par)
    end do
    
    if (par%masterproc) then
       write(iulog, '(a)') &
            'planar_conv>                          l2                    linf'
       do iq = 1,qsize
          a = global_shared_sum(2*iq-1)
          b = global_shared_sum(2*iq)
          reldif = sqrt(a/b)
          write(iulog, '(a,i2,es24.16,es24.16)') &
               'planar_conv> Q', iq, reldif, linf_num(iq)/linf_den(iq)
       end do
    end if
  end subroutine print_conv_planar_advection_results

  subroutine init(test_case, hybrid, hvcoord)
    character(len=*), intent(in):: test_case
    type (hybrid_t), intent(in):: hybrid
    type (hvcoord_t), intent(inout) :: hvcoord

    !$omp barrier
    !$omp master
    ! Major and minor test case codes.
    tc_major = test_case(17:17)
    tc_minor = test_case(18:18) ! currently unused
    ! Vertical dimension.
    call get_evenly_spaced_z(zi, zm, 0.d0, ztop)
    hvcoord%etai = exp(-zi/H)
    call set_hybrid_coefficients(hvcoord, hybrid, hvcoord%etai(1), 1.d0)
    call set_layer_locations(hvcoord, .true., hybrid%masterthread)
    !$omp end master
    !$omp barrier
  end subroutine init
  
end module planar_transport_tests

#endif
