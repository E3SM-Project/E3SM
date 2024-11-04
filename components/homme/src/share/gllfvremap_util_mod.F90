#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module gllfvremap_util_mod
  ! Utilities and extended tests for high-order, mass-conserving, optionally
  ! shape-preserving
  !     FV physics <-> GLL dynamics
  ! remap.
  !
  ! Tests operate at the level of gllfvremap_mod's API.
  !
  ! Utilities are to support homme_tool. So far these focus on creating
  ! topography files.
  !   A topography file contains FV data -- PHIS, SGH*, etc -- plus ncol_d and
  ! PHIS_d that are the original GLL data. The FV PHIS data are consistent with
  ! PHIS_d in the sense that an integral of either one over a finite volume
  ! subcell has the same value.
  !
  ! AMB 2019/07-2020/06 Initial

  use hybrid_mod, only: hybrid_t
  use kinds, only: real_kind
  use dimensions_mod, only: nelemd, np, nlev, nlevp, qsize
  use element_mod, only: element_t

  implicit none

  private

  real(kind=real_kind), parameter :: &
       zero = 0.0_real_kind, half = 0.5_real_kind, &
       one = 1.0_real_kind, two = 2.0_real_kind, &
       eps = epsilon(1.0_real_kind)

  type :: PhysgridData_t
     integer :: nphys
     real(kind=real_kind), allocatable :: ps(:,:), zs(:,:), T(:,:,:), uv(:,:,:,:), &
          omega_p(:,:,:), q(:,:,:,:)
  end type PhysgridData_t

  type :: State_t
     real(kind=real_kind), dimension(np,np,nlev) :: u, v, w, T, p, dp, rho, z
     real(kind=real_kind), dimension(np,np,nlevp) :: zi, wi
     real(kind=real_kind), dimension(np,np) :: ps, phis
  end type State_t

  type (PhysgridData_t), private :: pg_data

  logical :: is_sphere

  public :: &
       ! Test gllfvremap's main API.
       gfr_check_api, &
       ! Convert a topography file from pure GLL to physgrid format, as a
       ! convenience to avoid going through the full tool chain; see next
       ! function.
       gfr_convert_topo, &
       ! One part of the full physgrid from-scratch topography tool chain.
       gfr_pgn_to_smoothed_topo

contains
  
  function sphere2cart(sphere) result(cart)
    use coordinate_systems_mod, only: spherical_polar_t, cartesian3D_t, change_coordinates

    type (spherical_polar_t), intent(in) :: sphere
    type (cartesian3D_t) :: cart

    if (is_sphere) then
       cart = change_coordinates(sphere)
    else
       ! See conventions established in planar_mod::coordinates_atomic.
       cart%x = sphere%lon
       cart%y = sphere%lat
       cart%z = 0
    end if
  end function sphere2cart
  
  subroutine init(nphys)
    ! Init pg_data.

    integer, intent(in) :: nphys

    integer :: ncol

    ncol = nphys*nphys
    pg_data%nphys = nphys
    allocate(pg_data%ps(ncol,nelemd), pg_data%zs(ncol,nelemd), pg_data%T(ncol,nlev,nelemd), &
         pg_data%omega_p(ncol,nlev,nelemd), pg_data%uv(ncol,2,nlev,nelemd), &
         pg_data%q(ncol,nlev,qsize,nelemd))
  end subroutine init

  subroutine finish()
    ! Clean up pg_data.

    deallocate(pg_data%ps, pg_data%zs, pg_data%T, pg_data%uv, pg_data%omega_p, pg_data%q)
  end subroutine finish

  subroutine set_state(s, nt1, nt2, ntq, elem)
    ! Convenience wrapper to set_elem_state.

    use physical_constants, only: g
    use element_ops, only: set_elem_state

    type (State_t), intent(in) :: s
    integer, intent(in) :: nt1, nt2, ntq
    type (element_t), intent(inout) :: elem

    call set_elem_state(s%u, s%v, s%w, s%wi, s%T, s%ps, s%phis, s%p, s%dp, s%z, s%zi, &
         g, elem, nt1, nt2, ntq)
  end subroutine set_state

  subroutine set_gll_state(hvcoord, elem, nt1, nt2)
    ! Set all quantities used in gfr_dyn_to_fv_phys and gfr_fv_phys_to_dyn to
    ! C^inf functions for convergence testing (as well as property preservation
    ! testing, but the C^inf part is needed only for convergence testing).
    !   The C^inf functions are various sinusoidal 3D functions from which the
    ! values on the 2D sphere are extracted. This is an easy way to get a C^inf
    ! function on the sphere.

    use dimensions_mod, only: nlev, qsize
    use physical_constants, only: g, dd_pi, Lx, Ly
    use coordinate_systems_mod, only: cartesian3D_t
    use hybvcoord_mod, only: hvcoord_t
    use element_ops, only: get_field

    type (hvcoord_t) , intent(in) :: hvcoord
    type (element_t), intent(inout) :: elem
    integer, intent(in) :: nt1, nt2

    type (State_t) :: s1
    type (cartesian3D_t) :: p
    real(kind=real_kind) :: wr(np,np,nlev,2), fx, fy
    integer :: i, j, k, q, tl

    if (.not. is_sphere) then
       fx = 2*dd_pi/Lx
       fy = 2*dd_pi/Ly
    end if

    elem%state%Q(:,:,:,1) = zero ! moisture tracer is 0
    do j = 1,np
       do i = 1,np
          p = sphere2cart(elem%spherep(i,j))
          if (is_sphere) then
             do k = 1,nlev
                do q = 2,qsize
                   elem%state%Q(i,j,k,q) = one + &
                        half*sin((half + modulo(q,2))*p%x)* &
                        sin((half + modulo(q,3))*1.5d0*p%y)* &
                        sin((-2.3d0 + modulo(q,5))*(0.7d0 + p%z))
                end do
             end do
             s1%ps(i,j) = 1.0d3*(one + 0.05d0*sin(two*p%x+half)*sin(p%y+1.5d0)*sin(3*p%z+2.5d0))
             s1%phis(i,j) = one + half*sin(p%x-half)*sin(half*p%y+2.5d0)*sin(2*p%z-2.5d0)
             do k = 1,nlev
                ! u, v have to be set carefully because they are converted to
                ! contravariant velocity. Thus, at the poles, we need u, v to make
                ! sense to measure OOA correctly.
                wr(i,j,k,1) = sin(p%x)*sin(1.5*p%y)*cos(1.7*p%z)
                wr(i,j,k,2) = sin(half*p%x)*sin(1.5*p%y)*cos(half*dd_pi*p%z)
                elem%derived%omega_p(i,j,k) = wr(i,j,k,1)
             end do
             do k = 1,nlev
                s1%T(i,j,k) = one + half*sin(p%x+1.5d0)*sin(1.5d0*p%y+half)*sin(two*p%z-half)
             end do
          else
             do k = 1,nlev
                do q = 2,qsize
                   elem%state%Q(i,j,k,q) = one + half*sin(fx*p%x)*sin(fy*p%y)
                end do
             end do
             s1%ps(i,j) = 1.0d3*(one + 0.05d0*sin(two*fx*p%x+half)*sin(fy*p%y+1.5d0))
             s1%phis(i,j) = one + half*sin(fx*p%x-half)*sin(fy*p%y+2.5d0)
             do k = 1,nlev
                wr(i,j,k,1) = sin(fx*p%x)*sin(2*fy*fy*p%y)
                wr(i,j,k,2) = sin(fx*p%x)*sin(2*fy*p%y)
                elem%derived%omega_p(i,j,k) = wr(i,j,k,1)
             end do
             do k = 1,nlev
                s1%T(i,j,k) = one + half*sin(fx*p%x+1.5d0)*sin(fy*p%y+half)
             end do
          end if
       end do
    end do
    s1%u = wr(:,:,:,1)
    s1%v = wr(:,:,:,2)
    do k = 1,nlev
       s1%p(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*s1%ps
       s1%dp(:,:,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                      (hvcoord%hybi(k+1) - hvcoord%hybi(k))*s1%ps
    end do
    s1%z = zero
    s1%zi = zero
    ! a bit of a kludge
    call set_state(s1, nt1, nt2, nt1, elem)
    call get_field(elem, 'rho', wr(:,:,:,1), hvcoord, nt1, nt1)
    s1%w = -elem%derived%omega_p/(wr(:,:,:,1)*g)
    s1%wi(:,:,:nlev) = s1%w
    s1%wi(:,:,nlevp) = s1%w(:,:,nlev)
    call set_state(s1, nt1, nt2, nt1, elem)
    do q = 1,qsize
       do tl = nt1,nt2
          elem%state%Qdp(:,:,:,q,tl) = &
               elem%state%Q(:,:,:,q)*elem%state%dp3d(:,:,:,nt1)
       end do
    end do
  end subroutine set_gll_state

  function make_tendency(p) result(t)
    use coordinate_systems_mod, only: cartesian3D_t
    use physical_constants, only: dd_pi, Lx, Ly

    type (cartesian3D_t), intent(in) :: p
    real(kind=real_kind) :: t, r

    if (is_sphere) then
       r = sqrt(p%x*p%x + p%y*p%y + p%z*p%z)
       t = 0.25_real_kind*sin(3.2*p%x)*sin(4.2*p%y)*sin(0.7 + 2.3*p%z)
    else
       t = 0.25_real_kind*sin(2*dd_pi*p%x/Lx)*sin(2*dd_pi*p%y/Ly)
    end if
  end function make_tendency

  function run(hybrid, hvcoord, elem, nets, nete, nphys, tendency) result(nerr)
    ! Run 3 convergence and property-preservation whole-mesh
    ! tests. See below for descriptions of the three tests.

    use kinds, only: iulog
    use hybvcoord_mod, only: hvcoord_t
    use dimensions_mod, only: nlev, qsize
    use coordinate_systems_mod, only: cartesian3D_t
    use element_ops, only: get_temperature, get_field
    use prim_driver_base, only: applyCAMforcing_tracers
    use prim_advance_mod, only: applyCAMforcing_dynamics
    use parallel_mod, only: global_shared_buf, global_shared_sum
    use global_norms_mod, only: wrap_repro_sum
    use reduction_mod, only: ParallelMin, ParallelMax
    use physical_constants, only: g, p0, kappa
    use edge_mod, only: edgevpack_nlyr, edgevunpack_nlyr, edge_g
    use bndry_mod, only: bndry_exchangev
    use control_mod, only: ftype
    use gllfvremap_mod

    type (hybrid_t), intent(in) :: hybrid
    type (hvcoord_t) , intent(in) :: hvcoord
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete, nphys
    logical, intent(in) :: tendency

    character(32) :: msg
    type (cartesian3D_t) :: p
    real(kind=real_kind) :: wg(np,np,nlev), tend(np,np,nlev), f, a, b, rd, &
         qmin1(qsize+3), qmax1(qsize+3), qmin2, qmax2, mass1, mass2, &
         wg1(np,np,nlev), wg2(np,np,nlev), dt, pressure(np,np,nlev), &
         p_fv(np*np,nlev), wf1(np*np,nlev), wf2(np*np,nlev), wf3(np*np)
    integer :: nf, nf2, nt1, nt2, ie, i, j, k, q, qi, col, nerr
    logical :: domass

    nerr = 0
    nf = nphys
    nf2 = nf*nf
    nt1 = 1
    nt2 = 2
    dt = 0.42_real_kind

    !! Test 1.
    ! Test physgrid API.
    !   Test that if tendency is 0, then the original field is
    ! recovered with error eps ~ machine precision.
    !   OOA for pgN should min(N, 2). The OOA is limited to 2 because
    ! of the way the tendencies are created. Instead of setting the FV
    ! value to the average over the subcell, the sampled value at the
    ! cell midpoint (ref midpoint mapped to sphere) is used. This
    ! creates a 2nd-order error. Tests 2 and 3 test whether the remaps
    ! achieve design OOA, in particular OOA 3 for pg3 u, v, T, phis
    ! fields.

    ! Set analytical GLL values.
    do ie = nets,nete
       call set_gll_state(hvcoord, elem(ie), nt1, nt2)
    end do

    ! GLL -> FV.
    call gfr_dyn_to_fv_phys(hybrid, nt2, hvcoord, elem, nets, nete, &
         pg_data%ps, pg_data%zs, pg_data%T, pg_data%uv, pg_data%omega_p, pg_data%q)
    
    ! Set FV tendencies.
    if (tendency) then
       do ie = nets,nete
          do j = 1,nf
             do i = 1,nf
                col = nf*(j-1) + i
                call gfr_f_get_cartesian3d(ie, i, j, p)
                f = make_tendency(p)
                pg_data%uv(col,:,:,ie) = f/dt
                pg_data%T(col,:,ie) = f/dt
                ! no moisture adjustment => no dp3d adjustment
                pg_data%q(col,:,1,ie) = zero
                pg_data%q(col,:,2:qsize,ie) = pg_data%q(col,:,2:qsize,ie) + f
             end do
          end do
       end do
    else
       ! Test that if tendencies are 0, then the original fields are unchanged.
       pg_data%T = zero
       pg_data%uv = zero
    end if

    ! FV -> GLL.
    call gfr_fv_phys_to_dyn(hybrid, nt2, hvcoord, elem, nets, nete, &
         pg_data%T, pg_data%uv, pg_data%q)
    call gfr_f2g_dss(hybrid, elem, nets, nete)
    call gfr_pg1_reconstruct(hybrid, nt2, hvcoord, elem, nets, nete)

    ! Apply the tendencies.
    do ie = nets,nete
       if (ftype == 0) then
          do q = 1, qsize
             elem(ie)%derived%FQ(:,:,:,q) = elem(ie)%state%dp3d(:,:,:,nt1)* &
                  (elem(ie)%derived%FQ(:,:,:,q) - elem(ie)%state%Q(:,:,:,q))/dt
          end do
       end if
       call applyCAMforcing_tracers(elem(ie), hvcoord, nt2, nt2, dt, logical(ftype /= 0))
    end do
    call applyCAMforcing_dynamics(elem, hvcoord, nt2, dt, nets, nete)

    ! Test GLL state nt2 vs the original state nt1.
    if (hybrid%masterthread) write(iulog, '(a,l2)') 'gfrt> tendency', tendency
    tend = zero
    mass1 = zero; mass2 = zero
    qmin1 = one; qmax1 = -one
    qmin2 = one; qmax2 = -one
    do q = 2, qsize+4
       do ie = nets,nete
          do k = 1,nlev
             wg(:,:,k) = elem(ie)%spheremp
          end do
          if (tendency .and. q > 1) then
             do j = 1,np
                do i = 1,np
                   p = sphere2cart(elem(ie)%spherep(i,j))
                   tend(i,j,:) = make_tendency(p)
                end do
             end do
          end if
          if (q > qsize) then
             qi = q - qsize
             if (qi < 3) then
                global_shared_buf(ie,1) = &
                     sum(wg*( &
                     elem(ie)%state%v(:,:,qi,:,nt2) - &
                     (elem(ie)%state%v(:,:,qi,:,nt1) + tend))**2)
                global_shared_buf(ie,2) = &
                     sum(wg*(elem(ie)%state%v(:,:,qi,:,nt1) + tend)**2)
             elseif (qi == 3) then
                call get_temperature(elem(ie), wg1, hvcoord, nt1)
                call get_temperature(elem(ie), wg2, hvcoord, nt2)
                global_shared_buf(ie,1) = sum(wg*(wg2 - (wg1 + tend))**2)
                global_shared_buf(ie,2) = sum(wg*(wg1 + tend)**2)
             else
                ! Test omega_p, phis, ps. These were remapped to FV
                ! but don't get remapped to GLL. Make sure they all
                ! were remapped: the following should hold to nearly
                ! machine precision.
                !  omega_p
                !  True omega on GLL and FV grids.
                call get_field(elem(ie), 'omega', wg1, hvcoord, nt1, -1)
                call gfr_g2f_scalar(ie, elem(ie)%metdet, wg1, wf1)
                !  Convert omega_p on FV grid to omega for preqx
                call get_field(elem(ie), 'p', pressure, hvcoord, nt1, -1)
                call gfr_g2f_scalar(ie, elem(ie)%metdet, pressure, p_fv)
#ifdef MODEL_THETA_L
                wf2(:nf2,:) = pg_data%omega_p(:,:,ie)
#else                
                wf2(:nf2,:) = pg_data%omega_p(:,:,ie)*p_fv(:nf2,:)
#endif                
                !  Compare.
                global_shared_buf(ie,1) = sum((wf2(:nf2,:) - wf1(:nf2,:))**2)
                global_shared_buf(ie,2) = sum(wf1(:nf2,:)**2)
                !  phis
                call gfr_dyn_to_fv_phys_topo_elem(elem, ie, wf3)
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + &
                     sum((wf3(:nf2) - pg_data%zs(:,ie))**2)
                global_shared_buf(ie,2) = global_shared_buf(ie,2) + &
                     sum(pg_data%zs(:,ie)**2)
                !  ps
                wg(:,:,1) = elem(ie)%state%ps_v(:,:,nt1)
                wf1(:nf2,2) = pg_data%ps(:,ie)
                call gfr_g2f_scalar(ie, elem(ie)%metdet, wg(:,:,1:1), wf1(:,3:3))
                global_shared_buf(ie,1) = global_shared_buf(ie,1) + &
                     sum((wf1(:nf2,2) - wf1(:nf2,3))**2)
                global_shared_buf(ie,2) = global_shared_buf(ie,2) + &
                     sum(wf1(:nf2,2)**2)
             end if
          else
             global_shared_buf(ie,1) = &
                  sum(wg*( &
                  elem(ie)%state%Q(:,:,:,q) - &
                  (elem(ie)%state%Qdp(:,:,:,q,nt1)/elem(ie)%state%dp3d(:,:,:,nt1) + tend))**2)
             global_shared_buf(ie,2) = &
                  sum(wg*( &
                  elem(ie)%state%Qdp(:,:,:,q,nt1)/elem(ie)%state%dp3d(:,:,:,nt1) + tend)**2)
          end if
       end do
       call wrap_repro_sum(nvars=2, comm=hybrid%par%comm)
       rd = sqrt(global_shared_sum(1)/global_shared_sum(2))
       msg = ''
       if (.not. tendency .and. rd > 5*eps) then
          nerr = nerr + 1
          msg = ' ERROR'
       end if
       if (hybrid%masterthread) then
          write(iulog, '(a,i3,a,i3,es12.4,a8)') 'gfrt> test1 q l2', q, ' of', qsize, rd, msg
       end if
    end do

    !! Test 2.
    ! Test topo routines. For pgN, phis should have OOA min(N,2)
    ! because the field is limited.
  
    ! Stash initial phis for later comparison.
    do ie = nets,nete
       elem(ie)%derived%vstar(:,:,1,1) = elem(ie)%state%phis
    end do
    call gfr_dyn_to_fv_phys_topo(hybrid, elem, nets, nete, pg_data%zs)
    call gfr_fv_phys_to_dyn_topo(hybrid, elem, nets, nete, pg_data%zs)
    ! Do the DSS w/o (r)spheremp, as in inidata.F90. gfr_fv_phys_to_dyn_topo has
    ! prepped phis for this.
    do ie = nets, nete
       call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis, 1, 0, 1)
    end do
    call bndry_exchangeV(hybrid, edge_g)
    do ie = nets, nete
       call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis, 1, 0, 1)
    end do
    call gfr_pg1_reconstruct_topo(hybrid, elem, nets, nete)
    ! Compare GLL phis1 against GLL phis0.
    do ie = nets,nete
       global_shared_buf(ie,1) = sum(elem(ie)%spheremp*(elem(ie)%state%phis - &
            elem(ie)%derived%vstar(:,:,1,1))**2)
       global_shared_buf(ie,2) = sum(elem(ie)%spheremp*elem(ie)%derived%vstar(:,:,1,1)**2)
    end do
    call wrap_repro_sum(nvars=2, comm=hybrid%par%comm)
    if (hybrid%masterthread) then
       rd = sqrt(global_shared_sum(1)/global_shared_sum(2))
       write(iulog, '(a,es12.4)') 'gfrt> test2 topo l2', rd
    end if

    if (.not. tendency) return

    !! Test 3.
    ! Test FV fields that were not covered in the previous tests. This
    ! is done by copying them to look like tendencies.
    !   For pg4, u,v,T should have l2 errors that are eps because pg4
    ! reconstructions the fields exactly, even with DSS.
    !   For pg4, q should have l2 errors that converge at OOA >= 2. 2
    ! is the formal OOA b/c of the limiter. The limiter acts on the
    ! tendency, which in this case is exactly what is being examined.
    !   For pgN, N=1,2,3, u,v,T should have OOA N.
    !   For pgN, q should have OOA min(N,2).

    do ie = nets,nete
       call set_gll_state(hvcoord, elem(ie), nt1, nt2)
    end do
    call gfr_dyn_to_fv_phys(hybrid, nt2, hvcoord, elem, nets, nete, &
         pg_data%ps, pg_data%zs, pg_data%T, pg_data%uv, pg_data%omega_p, pg_data%q)
    ! Leave T, uv as they are. They will be mapped back as
    ! tendencies. Double Q so that this new value minus the
    ! original is Q.
    qmin1 = one; qmax1 = -one
    do ie = nets,nete
       pg_data%q(:nf2,:,:,ie) = two*pg_data%q(:nf2,:,:,ie)
       do q = 2,qsize
          qmin1(q) = min(qmin1(q), minval(elem(ie)%state%Q(:,:,1,q)))
          qmax1(q) = max(qmax1(q), maxval(elem(ie)%state%Q(:,:,1,q)))
          qmin1(q) = min(qmin1(q), minval(pg_data%q(:nf2,1,q,ie)))
          qmax1(q) = max(qmax1(q), maxval(pg_data%q(:nf2,1,q,ie)))
       end do
    end do
    call gfr_fv_phys_to_dyn(hybrid, nt2, hvcoord, elem, nets, nete, &
         pg_data%T, pg_data%uv, pg_data%q)
    call gfr_f2g_dss(hybrid, elem, nets, nete)
    call gfr_pg1_reconstruct(hybrid, nt2, hvcoord, elem, nets, nete)
    ! Don't apply forcings; rather, the forcing fields now have the
    ! remapped quantities we want to compare against the original.
    do q = 2, qsize+3
       domass = .true.
       mass1 = zero; mass2 = zero
       qmin2 = one; qmax2 = -one
       do ie = nets,nete
          do k = 1,nlev
             wg(:,:,k) = elem(ie)%spheremp
          end do
          if (q > qsize) then
             qi = q - qsize
             if (qi < 3) then
                ! With contravariant-velocity approach, we don't
                ! expect mass conservation.
                domass = .false.
                global_shared_buf(ie,3) = sum(wg(:,:,1)*elem(ie)%state%dp3d(:,:,1,nt1)* &
                     elem(ie)%derived%FM(:,:,qi,1))
                global_shared_buf(ie,4) = sum(wg(:,:,1)*elem(ie)%state%dp3d(:,:,1,nt1)* &
                     elem(ie)%state%v(:,:,qi,1,nt1))
                global_shared_buf(ie,1) = &
                     sum(wg*( &
                     elem(ie)%derived%FM(:,:,qi,:) - &
                     elem(ie)%state%v(:,:,qi,:,nt1))**2)
                global_shared_buf(ie,2) = &
                     sum(wg*elem(ie)%state%v(:,:,qi,:,nt1)**2)
             else
                call get_temperature(elem(ie), wg1, hvcoord, nt1)
                call get_field(elem(ie), 'p', wg2, hvcoord, nt1, -1)
                global_shared_buf(ie,1) = sum(wg*(elem(ie)%derived%FT - wg1)**2)
                global_shared_buf(ie,2) = sum(wg*wg1**2)
                wg1 = wg1*(p0/wg2)**kappa
                global_shared_buf(ie,4) = sum(wg(:,:,1)*elem(ie)%state%dp3d(:,:,1,nt1)*wg1(:,:,1))
                wg1 = elem(ie)%derived%FT
                wg1 = wg1*(p0/wg2)**kappa
                global_shared_buf(ie,3) = sum(wg(:,:,1)*elem(ie)%state%dp3d(:,:,1,nt1)*wg1(:,:,1))
             end if
          else
             ! Extrema in level 1.
             qmin2 = min(qmin2, minval(elem(ie)%derived%FQ(:,:,1,q)))
             qmax2 = max(qmax2, maxval(elem(ie)%derived%FQ(:,:,1,q)))
             ! Mass in level 1.
             global_shared_buf(ie,3) = sum(wg(:,:,1)*elem(ie)%state%dp3d(:,:,1,nt1)* &
                  elem(ie)%derived%FQ(:,:,1,q))
             global_shared_buf(ie,4) = sum(wg(:,:,1)*elem(ie)%state%dp3d(:,:,1,nt1)* &
                  two*elem(ie)%state%Q(:,:,1,q))
             ! l2 error in volume.
             global_shared_buf(ie,1) = &
                  sum(wg*( &
                  elem(ie)%derived%FQ(:,:,:,q) - &
                  two*elem(ie)%state%Q(:,:,:,q))**2)
             global_shared_buf(ie,2) = &
                  sum(wg*elem(ie)%state%Q(:,:,:,q)**2)
          end if
       end do
       call wrap_repro_sum(nvars=4, comm=hybrid%par%comm)
       qmin1(q) = ParallelMin(qmin1(q), hybrid)
       qmax1(q) = ParallelMax(qmax1(q), hybrid)
       qmin2 = ParallelMin(qmin2, hybrid)
       qmax2 = ParallelMax(qmax2, hybrid)
       rd = sqrt(global_shared_sum(1)/global_shared_sum(2))
       if (hybrid%masterthread) &
            write(iulog, '(a,i3,a,i3,es12.4)') 'gfrt> test3 q l2', q, ' of', qsize, rd
       b = max(abs(qmin1(q)), abs(qmax1(q)))
       if (q <= qsize .and. qmin2 < qmin1(q) - 5*eps*b .or. &
            qmax2 > qmax1(q) + 5*eps*b) then
          nerr = nerr + 1
          if (hybrid%masterthread) then
             write(iulog, '(a,i3,es12.4,es12.4,es12.4,es12.4,a)') 'gfrt> test3 q extrema', &
                  q, qmin1(q), qmin2-qmin1(q), qmax2-qmax1(q), qmax1(q), ' ERROR'
          end if
       end if
       if (domass) then
          a = global_shared_sum(3)
          b = global_shared_sum(4)
          f = 5
#ifdef HOMMEXX_BFB_TESTING
          ! Errors in bfb_pow mean mass conservation holds to only ~1e-6.
          if (q == qsize+3) f = 1.0e11_real_kind
#endif
          if (abs(b - a) > f*eps*abs(a)) then
             nerr = nerr + 1
             if (hybrid%masterthread) then
                write(iulog, '(a,i3,es12.4,es12.4,es12.4,a)') 'gfrt> test3 q mass', &
                     q, a, b, abs(b - a)/abs(a), ' ERROR'
             end if
          end if
       end if
    end do
  end function run

  function gfr_check_api(hybrid, nets, nete, hvcoord, elem) result(nerr)
    ! Drive run. Check nphys 1 through 4, ftypes 0 and 2, and pg1 with
    ! and without the OOA boost.

    use kinds, only: iulog
    use hybvcoord_mod, only: hvcoord_t
    use control_mod, only: ftype, geometry
    use gllfvremap_mod

    type (hybrid_t), intent(in) :: hybrid
    type (element_t), intent(inout) :: elem(:)
    type (hvcoord_t) , intent(in) :: hvcoord
    integer, intent(in) :: nets, nete

    integer :: nphys, ftype_in, ftype_idx, boost_idx, nerr
    logical :: boost_pg1

    is_sphere = trim(geometry) /= 'plane'

    ftype_in = ftype

    nerr = 0
    do nphys = 1, np
       do ftype_idx = 1,2
          do boost_idx = 1,2
             if (nphys > 1 .and. boost_idx > 1) exit
             boost_pg1 = boost_idx == 2

             ! This is meant to be called before threading starts.
             if (hybrid%ithr == 0) then
                ftype = 2
                if (ftype_idx == 2) ftype = 0
                ! check=2 means that the remap routines due
                ! element-level verification of properties and output
                ! messages if a property fails. check >= 1 means that
                ! global properties are checked.
                call gfr_init(hybrid%par, elem, nphys, 2, boost_pg1)
                call init(nphys)
             end if
             !$omp barrier

             if (hybrid%masterthread) write(iulog, '(a,i2)') 'gfrt> ftype', ftype

             nerr = nerr + run(hybrid, hvcoord, elem, nets, nete, nphys, .false.)
             nerr = nerr + run(hybrid, hvcoord, elem, nets, nete, nphys, .true.)

             ! This is meant to be called after threading ends.
             !$omp barrier
             if (hybrid%ithr == 0) then
                call gfr_finish()
                call finish()
             end if
             !$omp barrier
          end do
       end do
    end do

    !$omp barrier
    if (hybrid%ithr == 0) ftype = ftype_in
  end function gfr_check_api

  subroutine gfr_convert_topo(par, elem, nphys, intopofn, outtopoprefix)
    ! Read a pure-GLL topography file. Remap all fields to physgrid. Write a new

#if !defined(CAM) && !defined(SCREAM)
    use common_io_mod, only: varname_len
    use gllfvremap_mod, only: gfr_init, gfr_finish, gfr_dyn_to_fv_phys_topo_data, gfr_f_get_latlon
    use interpolate_driver_mod, only: read_gll_topo_file, write_physgrid_topo_file
    use physical_constants, only: dd_pi
#endif
    use parallel_mod, only: parallel_t

    type (parallel_t), intent(in) :: par
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nphys
    character(*), intent(in) :: intopofn, outtopoprefix

#if !defined(CAM) && !defined(SCREAM)
    real(real_kind), allocatable :: gll_fields(:,:,:,:), pg_fields(:,:,:), latlon(:,:,:)
    integer :: nf2, vari, ie, i, j, k
    logical :: square, augment
    character(len=varname_len) :: fieldnames(5)

    nf2 = nphys*nphys
    call gfr_init(par, elem, nphys)

    allocate(gll_fields(np,np,nelemd,5), pg_fields(nf2,nelemd,5), latlon(nf2,nelemd,2))

    call read_gll_topo_file(intopofn, elem, par, gll_fields, fieldnames)

    do vari = 1,size(fieldnames)
       if (trim(fieldnames(vari)) == 'PHIS') then
          do ie = 1,nelemd
             elem(ie)%state%phis = gll_fields(:,:,ie,vari)
          end do
          exit
       end if
    end do

    do vari = 1,size(fieldnames)
       square = fieldnames(vari)(1:3) == 'SGH'
       augment = trim(fieldnames(vari)) == 'SGH'
       call gfr_dyn_to_fv_phys_topo_data(par, elem, 1, nelemd, &
            gll_fields(:,:,:,vari), np*np*nelemd, pg_fields(:,:,vari), nf2*nelemd, &
            square, augment)
    end do

    do ie = 1,nelemd
       do j = 1,nphys
          do i = 1,nphys
             k = nphys*(j-1) + i
             call gfr_f_get_latlon(ie, i, j, latlon(k,ie,1), latlon(k,ie,2))
          end do
       end do
    end do
    ! Convert to degrees.
    latlon = latlon*(180.0_real_kind/dd_pi)

    call gfr_finish()

    call write_physgrid_topo_file(intopofn, outtopoprefix, elem, par, &
         gll_fields, pg_fields, latlon, fieldnames, nphys, &
         'Converted from '// trim(intopofn) // ' by HOMME gfr_convert_topo')

    deallocate(gll_fields, pg_fields, latlon)
#endif
  end subroutine gfr_convert_topo

  function gfr_pgn_to_smoothed_topo(par, elem, output_nphys, intopofn, outtopoprefix) result(stat)
#if !defined(CAM) && !defined(SCREAM)
    use common_io_mod, only: varname_len
    use gllfvremap_mod, only: gfr_init, gfr_finish, gfr_fv_phys_to_dyn_topo, &
         gfr_dyn_to_fv_phys_topo, gfr_f_get_latlon
    use interpolate_driver_mod, only: read_physgrid_topo_file, write_physgrid_smoothed_phis_file, &
       pio_read_phis
    use physical_constants, only: dd_pi
    use edge_mod, only: edgevpack_nlyr, edgevunpack_nlyr, edge_g
    use bndry_mod, only: bndry_exchangev
    use prim_driver_base, only: smooth_topo_datasets
    use hybrid_mod, only: hybrid_t, hybrid_create
#endif
    use parallel_mod, only: parallel_t

    type (parallel_t), intent(in) :: par
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: output_nphys
    character(*), intent(in) :: intopofn, outtopoprefix
    integer :: stat

#if !defined(CAM) && !defined(SCREAM)
    real(real_kind), allocatable :: gll_fields(:,:,:,:), pg_fields(:,:,:)
    integer :: intopo_nphys, ie, i, j, k, nvar, nf2
    character(len=varname_len) :: fieldnames(1)
    real(real_kind) :: rad2deg
    logical :: write_latlon
    type(hybrid_t) :: hybrid

    write_latlon = .false.

    nvar = 1
    if (write_latlon) nvar = 3

    allocate(gll_fields(np,np,nelemd,nvar), pg_fields(np*np,nelemd,nvar))
    hybrid = hybrid_create(par, 0, 1)


    ! Read the unsmoothed physgrid topo data from cube_to_target's first
    ! run. Here, pg4 will give the best quality.
    fieldnames(1) = 'PHIS'
    if (hybrid%masterthread) print *,'Attempting to read PG4 PHIS data for smoothing:'
    call read_physgrid_topo_file(intopofn, elem, par, fieldnames, intopo_nphys, pg_fields, stat)
    if (stat == 0) then
       ! Map this topo field to GLL.
       call gfr_init(par, elem, intopo_nphys)
       call gfr_fv_phys_to_dyn_topo(par, elem, pg_fields(:,:,1))
       do ie = 1,nelemd
          call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis, 1, 0, 1)
       end do
       call bndry_exchangeV(par, edge_g)
       do ie = 1,nelemd
          call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis, 1, 0, 1)
       end do
       call gfr_finish()
    else
       ! error trying to read PG input data. try GLL:
       if (hybrid%masterthread) print *,'Attempting to read GLL PHIS data for smoothing:'
       call pio_read_phis(elem,hybrid%par,'PHIS')
    endif

    call smooth_topo_datasets(elem, hybrid, 1, nelemd)
    do ie = 1,nelemd
       gll_fields(:,:,ie,1) = elem(ie)%state%phis
    end do

    ! Map the GLL data to the target physgrid, e.g., pg2.
    call gfr_init(par, elem, output_nphys)
    call gfr_dyn_to_fv_phys_topo(par, elem, pg_fields(:,:,1))
    if (write_latlon) then
       rad2deg = 180.0_real_kind/dd_pi
       do ie = 1,nelemd
          do j = 1,output_nphys
             do i = 1,output_nphys
                k = output_nphys*(j-1) + i
                call gfr_f_get_latlon(ie, i, j, pg_fields(k,ie,2), pg_fields(k,ie,3))
             end do
          end do
       end do
       nf2 = output_nphys*output_nphys
       pg_fields(:nf2,:,2:3) = pg_fields(:nf2,:,2:3)*rad2deg
       do ie = 1,nelemd
          do j = 1,np
             do i = 1,np
                gll_fields(i,j,ie,2) = elem(ie)%spherep(i,j)%lat*rad2deg
                gll_fields(i,j,ie,3) = elem(ie)%spherep(i,j)%lon*rad2deg
             end do
          end do
       end do
    end if
    call gfr_finish()

    ! Write the netcdf file that will be used as the --smoothed-topography file
    ! in the second run of cube_to_target. Only ncol, PHIS are needed for
    ! this. We include PHIS_d for use in the final assembled GLL-physgrid topo
    ! file. Optionally, we also include physgrid and GLL lat-lon data for easy
    ! visualization using just this file.
    call write_physgrid_smoothed_phis_file(outtopoprefix, elem, par, &
         gll_fields, pg_fields, output_nphys, &
         'Created from '// trim(intopofn) // ' by HOMME gfr_pgn_to_smoothed_topo', &
         write_latlon)

    deallocate(gll_fields, pg_fields)
#endif
    stat = 0
  end function gfr_pgn_to_smoothed_topo
end module gllfvremap_util_mod
