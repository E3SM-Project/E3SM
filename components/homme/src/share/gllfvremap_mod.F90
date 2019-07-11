#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module gllfvremap_mod
  ! High-order, mass-conserving, optionally shape-preserving
  !     FV physics <-> GLL dynamics
  ! remap.
  !
  ! AMB 2019/07 Initial

  use hybrid_mod, only: hybrid_t
  use kinds, only: real_kind
  use dimensions_mod, only: np, npsq, qsize, nelemd
  use element_mod, only: element_t
  use coordinate_systems_mod, only: spherical_polar_t

  implicit none

  private

  integer, parameter :: nphys_max = np

  real(kind=real_kind), parameter :: &
       zero = 0.0_real_kind, half = 0.5_real_kind, &
       one = 1.0_real_kind, two = 2.0_real_kind, &
       four = 4.0_real_kind, eps = epsilon(1.0_real_kind)

  ! Data type and functions for high-order, shape-preserving FV <-> GLL remap.
  type, public :: GllFvRemap_t
     integer :: nphys, npi
     logical :: check, have_fv_topo_file_phis
     real(kind=real_kind) :: tolfac ! for checking
     real(kind=real_kind) :: &
          ! Node or cell weights
          w_gg(np,np), &   ! GLL np
          w_ff(nphys_max*nphys_max), & ! FV nphys
          w_sgsg(np,np), & ! GLL npi
          ! Mixed mass matrices
          M_sgf(np,np,nphys_max,nphys_max), & ! GLL npi, FV nphys
          M_gf(np,np,nphys_max,nphys_max), &  ! GLL np,  FV nphys
          ! Interpolate from GLL npi to GLL np
          interp(np,np,np,np), &
          ! Remap FV nphys <-> GLL np
          g2f_remapd(np,np,nphys_max*nphys_max), &
          f2g_remapd(nphys_max*nphys_max,np,np)
     ! FV subcell areas; FV analogue of GLL elem(ie)%metdet arrays
     real(kind=real_kind), allocatable :: &
          fv_metdet(:,:), & ! (nphys*nphys,nelemd)
          ! Vector on ref elem -> vector on sphere
          D_f(:,:,:,:,:), &   ! (nphys,nphys,2,2,nelemd)
          ! Inverse of D_f
          Dinv_f(:,:,:,:,:), &
          qmin(:,:,:), qmax(:,:,:), &
          phis(:,:)
     type (spherical_polar_t), allocatable :: &
          spherep_f(:,:,:) ! (nphys,nphys,nelemd)
  end type GllFvRemap_t

  type (GllFvRemap_t), private :: gfr

  ! Main API.
  public :: &
       gfr_init, &
       gfr_finish, &
       gfr_dyn_to_fv_phys, &
       gfr_fv_phys_to_dyn, &
       gfr_dyn_to_fv_phys_topo, &
       gfr_fv_phys_to_dyn_topo, &
       gfr_f2g_dss

  ! For testing.
  public :: &
       gfr_test, &
       gfr_g2f_scalar, gfr_f2g_scalar, gfr_f_get_latlon, gfr_f_get_cartesian3d, &
       gfr_g_make_nonnegative, gfr_dyn_to_fv_phys_topo_elem

  interface gfr_dyn_to_fv_phys
     module procedure gfr_dyn_to_fv_phys_hybrid
     module procedure gfr_dyn_to_fv_phys_dom_mt
  end interface gfr_dyn_to_fv_phys

  interface gfr_fv_phys_to_dyn
     module procedure gfr_fv_phys_to_dyn_hybrid
     module procedure gfr_fv_phys_to_dyn_dom_mt
  end interface gfr_fv_phys_to_dyn

  interface gfr_fv_phys_to_dyn_topo
     module procedure gfr_fv_phys_to_dyn_topo_hybrid
     module procedure gfr_fv_phys_to_dyn_topo_dom_mt
  end interface gfr_fv_phys_to_dyn_topo

contains

  subroutine gfr_init(par, elem, nphys, check)
    use kinds, only: iulog
    use dimensions_mod, only: nlev
    use parallel_mod, only: parallel_t, abortmp
    use quadrature_mod, only : gausslobatto, quadrature_t
    use control_mod, only: ftype

    type (parallel_t), intent(in) :: par
    type (element_t), intent(in) :: elem(:)
    integer, intent(in) :: nphys
    logical, intent(in), optional :: check

    real(real_kind) :: R(npsq,nphys_max*nphys_max), tau(npsq)
    integer :: nphys2

    gfr%check = .false.
    if (present(check)) gfr%check = check

    gfr%tolfac = one
    if (par%masterproc) &
         write(iulog,  '(a,i3,a,l2,a,i2)') 'gfr> init nphys', nphys, ' check', gfr%check, ' ftype', ftype

    if (nphys > np) then
       ! The FV -> GLL map is defined only if nphys <= np. If we ever are
       ! interested in the case of nphys > np, we will need to write a different
       ! map. See "!assume" annotations for mathematical assumptions in
       ! particular routines.
       call abortmp('gllfvremap_mod: nphys must be <= np')
    end if
    if (qsize == 0) then
       call abortmp('gllfvremap_mod: qsize must be >= 1')
    end if

    gfr%have_fv_topo_file_phis = .false.
    gfr%nphys = nphys
    gfr%npi = max(2, nphys)
    nphys2 = nphys*nphys

    call gfr_init_w_gg(np, gfr%w_gg)
    call gfr_init_w_gg(gfr%npi, gfr%w_sgsg)
    call gfr_init_w_ff(nphys, gfr%w_ff)
    call gfr_init_M_gf(np, nphys, gfr%M_gf, .true.)
    gfr%g2f_remapd(:,:,:nphys2) = reshape(gfr%M_gf(:,:,:nphys,:nphys), (/np,np,nphys2/))
    call gfr_init_M_gf(gfr%npi, nphys, gfr%M_sgf, .false.)
    call gfr_init_R(gfr%npi, nphys, gfr%w_sgsg, gfr%M_sgf, R, tau)
    call gfr_init_interp_matrix(gfr%npi, gfr%interp)
    call gfr_init_f2g_remapd(gfr, R, tau)

    allocate(gfr%fv_metdet(nphys2,nelemd), &
         gfr%D_f(nphys,nphys,2,2,nelemd), gfr%Dinv_f(nphys,nphys,2,2,nelemd), &
         gfr%qmin(nlev,max(1,qsize),nelemd), gfr%qmax(nlev,max(1,qsize),nelemd), &
         gfr%phis(nphys2,nelemd), gfr%spherep_f(nphys,nphys,nelemd))
    call gfr_init_fv_metdet(elem, gfr)
    call gfr_init_geometry(elem, gfr)
  end subroutine gfr_init

  subroutine gfr_finish()
    if (.not. allocated(gfr%fv_metdet)) return
    deallocate(gfr%fv_metdet, gfr%D_f, gfr%Dinv_f, gfr%qmin, gfr%qmax, gfr%phis, &
         gfr%spherep_f)
  end subroutine gfr_finish

  subroutine gfr_dyn_to_fv_phys_hybrid(hybrid, nt, hvcoord, elem, nets, nete, &
       ps, phis, T, uv, omega_p, q)
    use dimensions_mod, only: nlev
    use hybvcoord_mod, only: hvcoord_t
    use element_ops, only: get_temperature
    use physical_constants, only: p0, kappa

    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: nt
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(in) :: elem(:)
    integer, intent(in) :: nets, nete
    real(kind=real_kind), intent(inout) :: ps(:,:), phis(:,:), T(:,:,:), &
         uv(:,:,:,:), omega_p(:,:,:), q(:,:,:,:)

    real(kind=real_kind), dimension(np,np,nlev) :: dp, dp_fv, wr1, wr2, p, p_fv
    real(kind=real_kind) :: qmin, qmax, ones(np,np)
    integer :: ie, nf, nf2, ncol, qi, qsize

    ones = one
    nf = gfr%nphys
    nf2 = nf*nf
    ncol = nf*nf

    qsize = size(q,3)
    
    do ie = nets,nete
       call gfr_g2f_scalar(ie, elem(ie)%metdet, elem(ie)%state%ps_v(:,:,nt:nt), &
            wr1(:,:,:1))
       ps(:ncol,ie) = reshape(wr1(:nf,:nf,1), (/ncol/))

       if (gfr%have_fv_topo_file_phis) then
          phis(:ncol,ie) = gfr%phis(:,ie)
       else
          call gfr_dyn_to_fv_phys_topo_elem(elem, ie, phis(:,ie))
       end if

       call calc_dp(hvcoord, elem(ie)%state%ps_v(:,:,nt), dp)
       call gfr_g2f_scalar(ie, elem(ie)%metdet, dp, dp_fv)

       call get_temperature(elem(ie), wr2, hvcoord, nt)
       call calc_p(hvcoord, elem(ie)%state%ps_v(:,:,nt), p)
       call gfr_g2f_scalar(ie, elem(ie)%metdet, p, p_fv)
       wr2 = wr2*(p/p0)**kappa
       call gfr_g2f_scalar_dp(gfr, ie, elem(ie)%metdet, dp, dp_fv, wr2, wr1)
       wr1(:nf,:nf,:) = wr1(:nf,:nf,:)/(p_fv(:nf,:nf,:)/p0)**kappa
       T(:ncol,:,ie) = reshape(wr1(:nf,:nf,:), (/ncol,nlev/))

       call gfr_g2f_vector(gfr, ie, elem, &
            elem(ie)%state%v(:,:,1,:,nt), elem(ie)%state%v(:,:,2,:,nt), &
            wr1, wr2)
       uv(:ncol,1,:,ie) = reshape(wr1(:nf,:nf,:), (/ncol,nlev/))
       uv(:ncol,2,:,ie) = reshape(wr2(:nf,:nf,:), (/ncol,nlev/))

       call gfr_g2f_scalar(ie, elem(ie)%metdet, elem(ie)%derived%omega_p, wr1)
       omega_p(:ncol,:,ie) = reshape(wr1(:nf,:nf,:), (/ncol,nlev/))

       do qi = 1,qsize
          call gfr_g2f_mixing_ratio(gfr, ie, elem(ie)%metdet, dp, dp_fv, &
               dp*elem(ie)%state%Q(:,:,:,qi), wr1)
          q(:ncol,:,qi,ie) = reshape(wr1(:nf,:nf,:), (/ncol,nlev/))
          if (gfr%check) then
             call check_g2f_mixing_ratio(gfr, hybrid, ie, qi, elem, dp, dp_fv, &
                  elem(ie)%state%Q(:,:,:,qi), wr1)
          end if
       end do
    end do
  end subroutine gfr_dyn_to_fv_phys_hybrid

  subroutine gfr_fv_phys_to_dyn_hybrid(hybrid, nt, dt, hvcoord, elem, nets, nete, T, uv, q)
    ! If ftype is in 1:4, then q is the full mixing ratio state and dt
    ! is not used; if it is not, then q is Qdp tendency, and dt must
    ! be the correct physics time step.

    use dimensions_mod, only: nlev
    use hybvcoord_mod, only: hvcoord_t
    use physical_constants, only: p0, kappa
    use control_mod, only: ftype

    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: nt
    real(kind=real_kind), intent(in) :: dt
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in), optional :: nets, nete
    real(kind=real_kind), intent(in) :: T(:,:,:), uv(:,:,:,:), q(:,:,:,:)

    real(kind=real_kind), dimension(np,np,nlev) :: dp, dp_fv, wr1, wr2, p, p_fv
    real(kind=real_kind) :: qmin, qmax
    integer :: ie, nf, ncol, k, qsize, qi
    logical :: q_adjustment

    nf = gfr%nphys
    ncol = nf*nf
    q_adjustment = ftype >= 1 .and. ftype <= 4
    qsize = size(q,3)

    do ie = nets,nete
       call calc_dp(hvcoord, elem(ie)%state%ps_v(:,:,nt), dp)
       call gfr_g2f_scalar(ie, elem(ie)%metdet, dp, dp_fv)

       wr1(:nf,:nf,:) = reshape(uv(:ncol,1,:,ie), (/nf,nf,nlev/))
       wr2(:nf,:nf,:) = reshape(uv(:ncol,2,:,ie), (/nf,nf,nlev/))
       call gfr_f2g_vector(gfr, ie, elem, &
            wr1, wr2, elem(ie)%derived%FM(:,:,1,:), elem(ie)%derived%FM(:,:,2,:))

       call calc_p(hvcoord, elem(ie)%state%ps_v(:,:,nt), p)
       call gfr_g2f_scalar(ie, elem(ie)%metdet, p, p_fv)
       wr1(:nf,:nf,:) = reshape(T(:ncol,:,ie), (/nf,nf,nlev/))
       wr1(:nf,:nf,:) = wr1(:nf,:nf,:)*(p_fv(:nf,:nf,:)/p0)**kappa
       call gfr_f2g_scalar_dp(gfr, ie, elem(ie)%metdet, dp_fv, dp, wr1, elem(ie)%derived%FT)
       elem(ie)%derived%FT = elem(ie)%derived%FT/(p/p0)**kappa

       do qi = 1,qsize
          if (q_adjustment) then
             ! FV Q_ten
             !   GLL Q0 -> FV Q0
             call gfr_g2f_mixing_ratio(gfr, ie, elem(ie)%metdet, dp, dp_fv, &
                  dp*elem(ie)%state%Q(:,:,:,qi), wr1)
             !   FV Q_ten = FV Q1 - FV Q0
             wr1(:nf,:nf,:) = reshape(q(:ncol,:,qi,ie), (/nf,nf,nlev/)) - wr1(:nf,:nf,:)
             ! GLL Q_ten
             call gfr_f2g_scalar_dp(gfr, ie, elem(ie)%metdet, dp_fv, dp, wr1, wr2)
             ! GLL Q1
             elem(ie)%derived%FQ(:,:,:,qi) = elem(ie)%state%Q(:,:,:,qi) + wr2
             ! Get limiter bounds.
             do k = 1,nlev
                gfr%qmin(k,qi,ie) = minval(q(:ncol,k,qi,ie))
                gfr%qmax(k,qi,ie) = maxval(q(:ncol,k,qi,ie))
             end do
          else
             ! FV Q_ten
             wr1(:nf,:nf,:) = reshape(q(:ncol,:,qi,ie), (/nf,nf,nlev/))
             wr1(:nf,:nf,:) = dt*wr1(:nf,:nf,:)/dp_fv(:nf,:nf,:)
             ! GLL Q_ten
             call gfr_f2g_scalar_dp(gfr, ie, elem(ie)%metdet, dp_fv, dp, wr1, wr2)
             ! GLL Q1
             elem(ie)%derived%FQ(:,:,:,qi) = elem(ie)%state%Q(:,:,:,qi) + wr2
             ! GLL Q0 -> FV Q0
             call gfr_g2f_mixing_ratio(gfr, ie, elem(ie)%metdet, dp, dp_fv, &
                  dp*elem(ie)%state%Q(:,:,:,qi), wr2)
             ! FV Q1
             wr2(:nf,:nf,:) = wr2(:nf,:nf,:) + wr1(:nf,:nf,:)
             ! Get limiter bounds.
             do k = 1,nlev
                gfr%qmin(k,qi,ie) = minval(wr2(:nf,:nf,k))
                gfr%qmax(k,qi,ie) = maxval(wr2(:nf,:nf,k))
             end do
          end if
       end do
    end do

    ! Halo exchange limiter bounds.
    call gfr_f2g_mixing_ratios_he(hybrid, nets, nete, gfr%qmin(:,:,nets:nete), &
         gfr%qmax(:,:,nets:nete))

    do ie = nets,nete
       call calc_dp(hvcoord, elem(ie)%state%ps_v(:,:,nt), dp)
       do qi = 1,qsize
          ! Limit GLL Q1.
          if (gfr%check) wr1 = elem(ie)%derived%FQ(:,:,:,qi)
          do k = 1,nlev
             ! Augment bounds with GLL Q0 bounds. This assures that if
             ! the tendency is 0, GLL Q1 = GLL Q0.
             gfr%qmin(k,qi,ie) = min(minval(elem(ie)%state%Q(:,:,k,qi)), gfr%qmin(k,qi,ie))
             gfr%qmax(k,qi,ie) = max(maxval(elem(ie)%state%Q(:,:,k,qi)), gfr%qmax(k,qi,ie))
             ! Final GLL Q1, except for DSS, which is not done in this routine.
             call limiter_clip_and_sum(np, elem(ie)%spheremp, gfr%qmin(k,qi,ie), &
                  gfr%qmax(k,qi,ie), dp(:,:,k), elem(ie)%derived%FQ(:,:,k,qi))
          end do
          if (gfr%check) then
             call check_f2g_mixing_ratio(gfr, hybrid, ie, qi, elem, gfr%qmin(:,qi,ie), &
                  gfr%qmax(:,qi,ie), dp, wr1, elem(ie)%derived%FQ(:,:,:,qi))
          end if
          if (.not. q_adjustment) then
             ! Convert to a tendency.
             elem(ie)%derived%FQ(:,:,:,qi) = &
                  dp*(elem(ie)%derived%FQ(:,:,:,qi) - elem(ie)%state%Q(:,:,:,qi))/dt
          end if
       end do
    end do
  end subroutine gfr_fv_phys_to_dyn_hybrid

  subroutine gfr_dyn_to_fv_phys_topo(hybrid, elem, nets, nete, phis)
    type (hybrid_t), intent(in) :: hybrid
    type (element_t), intent(in) :: elem(:)
    integer, intent(in), optional :: nets, nete
    real(kind=real_kind), intent(out) :: phis(:,:)

    integer :: ie

    do ie = nets,nete
       call gfr_dyn_to_fv_phys_topo_elem(elem, ie, phis(:,ie))
    end do
  end subroutine gfr_dyn_to_fv_phys_topo

  subroutine gfr_dyn_to_fv_phys_topo_elem(elem, ie, phis)
    type (element_t), intent(in) :: elem(:)
    integer, intent(in), optional :: ie
    real(kind=real_kind), intent(out) :: phis(:)

    real(kind=real_kind) :: wr(np,np,2), ones(np,np), qmin, qmax
    integer :: nf, nf2, ncol

    ones = one
    nf = gfr%nphys
    nf2 = nf*nf
    ncol = nf*nf

    wr(:,:,1) = elem(ie)%state%phis(:,:)
    call gfr_g2f_scalar(ie, elem(ie)%metdet, wr(:,:,1:1), wr(:,:,2:2))
    qmin = minval(elem(ie)%state%phis)
    qmax = maxval(elem(ie)%state%phis)
    wr(:nf,:nf,1) = reshape(gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie), (/nf,nf/))
    call limiter_clip_and_sum(gfr%nphys, wr(:,:,1), qmin, qmax, ones, wr(:nf,:nf,2))
    phis(:ncol) = reshape(wr(:nf,:nf,2), (/ncol/))
  end subroutine gfr_dyn_to_fv_phys_topo_elem

  subroutine gfr_fv_phys_to_dyn_topo_hybrid(hybrid, elem, nets, nete, phis)
    use kinds, only: iulog
    use edge_mod, only: edgeVpack_nlyr, edgeVunpack_nlyr, edge_g
    use bndry_mod, only: bndry_exchangeV

    type (hybrid_t), intent(in) :: hybrid
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in), optional :: nets, nete
    real(kind=real_kind), intent(in) :: phis(:,:)

    real(kind=real_kind) :: wr(np,np,2), ones(np,np,1)
    integer :: ie, nf, ncol

    ones = one
    nf = gfr%nphys
    ncol = nf*nf
    ! For now, map GLL topo back to FV, as ne30pg2 runs otherwise show grid
    ! imprint in the PHIS diagnostic. We may revise this choice later.
    !gfr%have_fv_topo_file_phis = .true.

    do ie = nets,nete
       gfr%phis(:,ie) = phis(:ncol,ie)
       wr(:nf,:nf,1) = reshape(phis(:ncol,ie), (/nf,nf/))
       gfr%qmin(:,:,ie) = minval(wr(:nf,:nf,1))
       gfr%qmax(:,:,ie) = maxval(wr(:nf,:nf,1))
       call gfr_f2g_scalar(ie, elem(ie)%metdet, wr(:,:,1:1), wr(:,:,2:2))
       elem(ie)%state%phis = wr(:,:,2)
    end do

    if (hybrid%par%dynproc) then
       call gfr_f2g_mixing_ratios_he(hybrid, nets, nete, gfr%qmin(:,:,nets:nete), &
            gfr%qmax(:,:,nets:nete))
    end if

    do ie = nets,nete
       if (gfr%check) wr(:,:,1) = elem(ie)%state%phis
       call limiter_clip_and_sum(np, elem(ie)%spheremp, gfr%qmin(1,1,ie), &
            gfr%qmax(1,1,ie), ones(:,:,1), elem(ie)%state%phis)
       if (gfr%check) then
          if (gfr%qmin(1,1,ie) < zero) then
             write(iulog,*) 'gfr> topo min:', hybrid%par%rank, hybrid%ithr, ie, gfr%qmin(1,1,ie)
          end if
          wr(:,:,2) = elem(ie)%state%phis
          call check_f2g_mixing_ratio(gfr, hybrid, ie, 1, elem, gfr%qmin(:1,1,ie), &
               gfr%qmax(:1,1,ie), ones, wr(:,:,:1), wr(:,:,2:))
       end if
    end do

    do ie = nets,nete
       elem(ie)%state%phis = elem(ie)%state%phis*elem(ie)%spheremp*elem(ie)%rspheremp
    end do
  end subroutine gfr_fv_phys_to_dyn_topo_hybrid

  subroutine gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)
    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use hybrid_mod, only: hybrid_t, hybrid_create
    use thread_mod, only: omp_get_thread_num, hthreads

    type (parallel_t), intent(in) :: par
    type (domain1d_t), intent(in) :: dom_mt(:)
    type (hybrid_t), intent(out) :: hybrid
    integer, intent(out) :: nets, nete

    integer :: ithr

    ithr = omp_get_thread_num()
    nets = dom_mt(ithr+1)%start
    nete = dom_mt(ithr+1)%end
    hybrid = hybrid_create(par, ithr, hthreads)
  end subroutine gfr_hybrid_create

  subroutine gfr_dyn_to_fv_phys_dom_mt(par, dom_mt, nt, hvcoord, elem, &
       ps, phis, T, uv, omega_p, q)
    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use hybvcoord_mod, only: hvcoord_t
    use thread_mod, only: hthreads

    type (parallel_t), intent(in) :: par
    type (domain1d_t), intent(in) :: dom_mt(:)
    integer, intent(in) :: nt
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(inout) :: ps(:,:), phis(:,:), T(:,:,:), &
         uv(:,:,:,:), omega_p(:,:,:), q(:,:,:,:)

    type (hybrid_t) :: hybrid
    integer :: nets, nete

#ifdef HORIZ_OPENMP
    !$omp parallel num_threads(hthreads), default(shared), private(nets,nete,hybrid)
#endif
    call gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)
    call gfr_dyn_to_fv_phys_hybrid(hybrid, nt, hvcoord, elem, nets, nete, &
         ps, phis, T, uv, omega_p, q)
#ifdef HORIZ_OPENMP
    !$omp end parallel
#endif
  end subroutine gfr_dyn_to_fv_phys_dom_mt

  subroutine gfr_fv_phys_to_dyn_dom_mt(par, dom_mt, nt, dt, hvcoord, elem, T, uv, q)
    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use hybvcoord_mod, only: hvcoord_t
    use thread_mod, only: hthreads

    type (parallel_t), intent(in) :: par
    type (domain1d_t), intent(in) :: dom_mt(:)
    integer, intent(in) :: nt
    real(kind=real_kind), intent(in) :: dt
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(inout) :: elem(:)
    real(kind=real_kind), intent(inout) :: T(:,:,:), uv(:,:,:,:), q(:,:,:,:)

    type (hybrid_t) :: hybrid
    integer :: nets, nete

#ifdef HORIZ_OPENMP
    !$omp parallel num_threads(hthreads), default(shared), private(nets,nete,hybrid)
#endif
    call gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)
    call gfr_fv_phys_to_dyn_hybrid(hybrid, nt, dt, hvcoord, elem, nets, nete, T, uv, q)
#ifdef HORIZ_OPENMP
    !$omp end parallel
#endif
  end subroutine gfr_fv_phys_to_dyn_dom_mt

  subroutine gfr_fv_phys_to_dyn_topo_dom_mt(par, dom_mt, elem, phis)
    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use thread_mod, only: hthreads

    type (parallel_t), intent(in) :: par
    type (domain1d_t), intent(in) :: dom_mt(:)
    type (element_t), intent(inout) :: elem(:)
    real(kind=real_kind), intent(in) :: phis(:,:)

    type (hybrid_t) :: hybrid
    integer :: nets, nete

#ifdef HORIZ_OPENMP
    !$omp parallel num_threads(hthreads), default(shared), private(nets,nete,hybrid)
#endif
    call gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)
    call gfr_fv_phys_to_dyn_topo_hybrid(hybrid, elem, nets, nete, phis)
#ifdef HORIZ_OPENMP
    !$omp end parallel
#endif
  end subroutine gfr_fv_phys_to_dyn_topo_dom_mt

  subroutine gfr_init_w_gg(np, w_gg)
    use quadrature_mod, only : gausslobatto, quadrature_t
    
    integer, intent(in) :: np
    real(kind=real_kind), intent(out) :: w_gg(:,:)

    type (quadrature_t) :: gll
    integer :: i,j

    gll = gausslobatto(np)

    do j = 1,np
       do i = 1,np
          w_gg(i,j) = gll%weights(i)*gll%weights(j)
       end do
    end do
    
    call gll_cleanup(gll)
  end subroutine gfr_init_w_gg

  subroutine gfr_init_w_ff(nphys, w_ff)
    integer, intent(in) :: nphys
    real(kind=real_kind), intent(out) :: w_ff(:)

    w_ff(:nphys*nphys) = two*two/real(nphys*nphys, real_kind)
  end subroutine gfr_init_w_ff

  subroutine gll_cleanup(gll)
    use quadrature_mod, only : quadrature_t

    type (quadrature_t), intent(inout) :: gll

    deallocate(gll%points, gll%weights)
  end subroutine gll_cleanup

  subroutine calc_dp(hvcoord, ps, dp)
    use hybvcoord_mod, only: hvcoord_t
    use dimensions_mod, only: nlev

    type (hvcoord_t), intent(in) :: hvcoord
    real(kind=real_kind), intent(in) :: ps(:,:)
    real(kind=real_kind), intent(out) :: dp(:,:,:)

    integer :: k

    do k = 1,nlev
       dp(:,:,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
            (hvcoord%hybi(k+1) - hvcoord%hybi(k))*ps
    end do
  end subroutine calc_dp

  subroutine calc_p(hvcoord, ps, p)
    use hybvcoord_mod, only: hvcoord_t
    use dimensions_mod, only: nlev

    type (hvcoord_t), intent(in) :: hvcoord
    real(kind=real_kind), intent(in) :: ps(:,:)
    real(kind=real_kind), intent(out) :: p(:,:,:)

    integer :: k

    do k = 1,nlev
       p(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps
    end do
  end subroutine calc_p

  subroutine eval_lagrange_bases(gll, np, x, y)
    ! Evaluate the GLL basis functions at x in [-1,1], writing the values to
    ! y(1:np).
    use quadrature_mod, only : quadrature_t
    
    type (quadrature_t), intent(in) :: gll
    integer, intent(in) :: np
    real(kind=real_kind), intent(in) :: x ! in [-1,1]
    real(kind=real_kind), intent(out) :: y(:)

    integer :: i, j
    real(kind=real_kind) :: f

    do i = 1,np
       f = one
       do j = 1,np
          if (j /= i) then
             f = f*((x - gll%points(j))/(gll%points(i) - gll%points(j)))
          end if
       end do
       y(i) = f
    end do
  end subroutine eval_lagrange_bases

  subroutine gfr_init_M_gf(np, nphys, M_gf, scale)
    use quadrature_mod, only : gausslobatto, quadrature_t

    integer, intent(in) :: np, nphys
    real(kind=real_kind), intent(out) :: M_gf(:,:,:,:)
    logical, intent(in) :: scale

    type (quadrature_t) :: gll
    integer :: gi, gj, fi, fj, qi, qj
    real(kind=real_kind) :: xs, xe, ys, ye, ref, bi(np), bj(np)

    gll = gausslobatto(np)

    M_gf = zero

    do fj = 1,nphys
       ! The subcell is [xs,xe]x[ys,ye].
       xs = two*real(fj-1, real_kind)/real(nphys, real_kind) - one
       xe = two*real(fj, real_kind)/real(nphys, real_kind) - one
       do fi = 1,nphys
          ys = two*real(fi-1, real_kind)/real(nphys, real_kind) - one
          ye = two*real(fi, real_kind)/real(nphys, real_kind) - one
          ! Use GLL quadrature within this subcell.
          do qj = 1,np
             ! (xref,yref) are w.r.t. the [-1,1]^2 reference domain mapped to
             ! the subcell.
             ref = xs + half*(xe - xs)*(1 + gll%points(qj))
             call eval_lagrange_bases(gll, np, ref, bj)
             do qi = 1,np
                ref = ys + half*(ye - ys)*(1 + gll%points(qi))
                call eval_lagrange_bases(gll, np, ref, bi)
                do gj = 1,np
                   do gi = 1,np
                      ! Accumulate each GLL basis's contribution to this
                      ! subcell.
                      M_gf(gi,gj,fi,fj) = M_gf(gi,gj,fi,fj) + &
                           gll%weights(qi)*gll%weights(qj)*bi(gi)*bj(gj)
                   end do
                end do
             end do
          end do
       end do
    end do

    M_gf = M_gf/real(nphys*nphys, real_kind)

    if (scale) then
       ! Scale so the sum over FV subcells gives the GLL weights to machine
       ! precision.
       do gj = 1,np
          do gi = 1,np
             M_gf(gi,gj,:nphys,:nphys) = M_gf(gi,gj,:nphys,:nphys)* &
                  ((gll%weights(gi)*gll%weights(gj))/ &
                  sum(M_gf(gi,gj,:nphys,:nphys)))
          end do
       end do
    end if

    call gll_cleanup(gll)
  end subroutine gfr_init_M_gf

  subroutine gfr_init_R(np, nphys, w_gg, M_gf, R, tau)
    ! We want to solve
    !     min_g 1/2 g'M_gg g - g' M_gf f
    !      st   M_gf' g = M_ff f,
    ! which gives
    !     [M_gg -M_gf] [g] = [M_gf f]
    !     [M_gf'  0  ] [y]   [M_ff f].
    ! Recall M_gg, M_ff are diag. Let
    !     S = M_gf' inv(M_gg) M_gf.
    ! Then
    !     g = inv(M_gg) M_gf inv(S) M_ff f.
    ! Compute the QR factorization sqrt(inv(M_gg)) M_gf = Q R so that S =
    ! R'R. In this module, we can take M_gg = diag(w_gg) and M_ff = diag(w_ff)
    ! with no loss of accuracy.
    !   If nphys = np, then the problem reduces to
    !     M_gf' g = M_ff f.
    ! M_gf is symmetric. We could use the same computations as above
    ! or, to gain a bit more accuracy, compute the simpler
    !     M_gf = Q R
    ! and later solve
    !     R'Q' g = M_ff f.
    !   In either case, all of this is one-time initialization; during
    ! time stepping, just a matvec is computed.
    !
    !assume nphys <= np

    integer, intent(in) :: np, nphys
    real(kind=real_kind), intent(in) :: w_gg(:,:), M_gf(:,:,:,:)
    real(kind=real_kind), intent(out) :: R(:,:), tau(:)

    real(kind=real_kind) :: wrk(np*np*nphys*nphys), v
    integer :: gi, gj, fi, fj, npsq, info

    do fj = 1,nphys
       do fi = 1,nphys
          do gi = 1,np
             do gj = 1,np
                v = M_gf(gi,gj,fi,fj)
                if (nphys < np) v = v/sqrt(w_gg(gi,gj))
                R(np*(gi-1) + gj, nphys*(fi-1) + fj) = v
             end do
          end do
       end do
    end do
    call dgeqrf(np*np, nphys*nphys, R, size(R,1), tau, wrk, np*np*nphys*nphys, info)
  end subroutine gfr_init_R

  subroutine gfr_init_interp_matrix(npsrc, interp)
    use quadrature_mod, only : gausslobatto, quadrature_t

    integer, intent(in) :: npsrc
    real(kind=real_kind), intent(out) :: interp(:,:,:,:)

    type (quadrature_t) :: glls, gllt
    integer :: si, sj, ti, tj
    real(kind=real_kind) :: bi(npsrc), bj(npsrc)

    glls = gausslobatto(npsrc)
    gllt = gausslobatto(np)

    do tj = 1,np
       call eval_lagrange_bases(glls, npsrc, real(gllt%points(tj), real_kind), bj)
       do ti = 1,np
          call eval_lagrange_bases(glls, npsrc, real(gllt%points(ti), real_kind), bi)
          do sj = 1,npsrc
             do si = 1,npsrc
                interp(si,sj,ti,tj) = bi(si)*bj(sj)
             end do
          end do
       end do
    end do

    call gll_cleanup(glls)
    call gll_cleanup(gllt)
  end subroutine gfr_init_interp_matrix

  subroutine gfr_init_f2g_remapd(gfr, R, tau)
    !assume nphys <= np

    type (GllFvRemap_t), intent(inout) :: gfr
    real(kind=real_kind), intent(in) :: R(:,:), tau(:)

    integer :: nf, fi, fj, gi, gj
    real(kind=real_kind) :: f(np,np), g(np,np)

    ! Apply gfr_init_f2g_remapd_op to the Id matrix to get the remap operator's
    ! matrix representation.
    gfr%f2g_remapd = zero
    f = zero
    nf = gfr%nphys
    do fi = 1,nf
       do fj = 1,nf
          f(fi,fj) = one
          call gfr_f2g_remapd_op(gfr, R, tau, f, g)
          gfr%f2g_remapd(fi + (fj-1)*nf,:,:) = g
          f(fi,fj) = zero
       end do
    end do
  end subroutine gfr_init_f2g_remapd

  subroutine gfr_f2g_remapd_op(gfr, R, tau, f, g)
    !assume nphys <= np

    type (GllFvRemap_t), intent(in) :: gfr
    real(kind=real_kind), intent(in) :: R(:,:), tau(:), f(:,:)
    real(kind=real_kind), intent(out) :: g(:,:)

    integer :: nf, nf2, npi, np2, gi, gj, fi, fj, info
    real(kind=real_kind) :: accum, wrk(gfr%nphys,gfr%nphys), x(np,np), wr(np*np)

    nf = gfr%nphys
    nf2 = nf*nf
    npi = gfr%npi
    np2 = np*np

    ! Solve the constrained projection described in gfr_init_R:
    !     g = inv(M_sgsg) M_sgf inv(S) M_ff f
    wrk = reshape(gfr%w_ff(:nf2), (/nf,nf/))*f(:nf,:nf)
    if (nf == npi) then
       call dtrsm('l', 'u', 't', 'n', nf2, 1, one, R, size(R,1), wrk, nf2)
       call dormqr('l', 'n', nf2, 1, nf2, R, size(R,1), tau, wrk, nf2, wr, np2, info)
       g(:npi,:npi) =  wrk
    else
       call dtrtrs('u', 't', 'n', nf2, 1, R, size(R,1), wrk, nf2, info)
       call dtrtrs('u', 'n', 'n', nf2, 1, R, size(R,1), wrk, nf2, info)
       g(:npi,:npi) = zero
       do fj = 1,nf
          do fi = 1,nf
             do gj = 1,npi
                do gi = 1,npi
                   g(gi,gj) = g(gi,gj) + gfr%M_sgf(gi,gj,fi,fj)*wrk(fi,fj)
                end do
             end do
          end do
       end do
    end if
    if (npi < np) then
       if (nf == npi) then
          x(:nf,:nf) = g(:nf,:nf)
       else
          ! Finish the projection:
          !     wrk = inv(M_sgsg) g
          do gj = 1,npi
             do gi = 1,npi
                x(gi,gj) = g(gi,gj)/gfr%w_sgsg(gi,gj)
             end do
          end do
       end if
       ! Interpolate from npi to np; if npi = np, this is just the Id matrix.
       do fj = 1,np
          do fi = 1,np
             accum = zero
             do gj = 1,npi
                do gi = 1,npi
                   accum = accum + gfr%interp(gi,gj,fi,fj)*x(gi,gj)
                end do
             end do
             g(fi,fj) = accum
          end do
       end do
    elseif (nf < npi) then
       ! Finish the projection.
       do gj = 1,np
          do gi = 1,np
             g(gi,gj) = g(gi,gj)/gfr%w_gg(gi,gj)
          end do
       end do
    end if
  end subroutine gfr_f2g_remapd_op

  subroutine gfr_init_fv_metdet(elem, gfr)
    type (element_t), intent(in) :: elem(:)
    type (GllFvRemap_t), intent(inout) :: gfr

    real(kind=real_kind) :: ones(np*np), ones2(np,np), wrk(np,np)
    integer :: nf, nf2, ie

    nf = gfr%nphys
    nf2 = nf*nf
    ones = one
    ones2 = one
    do ie = 1,nelemd
       call gfr_g2f_remapd(gfr, elem(ie)%metdet, ones, ones2, wrk)
       gfr%fv_metdet(:,ie) = reshape(wrk(:nf,:nf), (/nf2/))
    end do
  end subroutine gfr_init_fv_metdet

  subroutine gfr_f_ref_center(nphys, i, a)
    ! FV subcell center in ref [-1,1]^2 coord.

    integer, intent(in) :: nphys, i
    real(kind=real_kind), intent(out) :: a

    a = two*((real(i-1, real_kind) + half)/real(nphys, real_kind)) - one
  end subroutine gfr_f_ref_center

  subroutine gfr_f_ref_edges(nphys, i_fv, a)
    ! FV subcell edges in ref [-1,1]^2 coord.

    integer, intent(in) :: nphys, i_fv
    real(kind=real_kind), intent(out) :: a(2)

    integer :: i

    do i = 0,1
       a(i+1) = two*(real(i_fv+i-1, real_kind)/real(nphys, real_kind)) - one
    end do
  end subroutine gfr_f_ref_edges

  subroutine gfr_init_geometry(elem, gfr)
    use coordinate_systems_mod, only: cartesian3D_t, change_coordinates, sphere_tri_area
    use cube_mod, only: Dmap, ref2sphere
    use control_mod, only: cubed_sphere_map

    type (element_t), intent(in) :: elem(:)
    type (GllFvRemap_t), intent(inout) :: gfr

    type (spherical_polar_t) :: p_sphere
    type (cartesian3D_t) :: fv_corners_xyz(2,2)
    real(kind=real_kind) :: wrk(2,2), det, a, b, ae(2), be(2), tmp, spherical_area
    integer :: ie, nf, nf2, i, j, k, ai, bi

    nf = gfr%nphys
    nf2 = nf*nf
    
    if (cubed_sphere_map == 2) then
       ! Reset fv_metdet to match flux coupler's definition of FV
       ! subcell area, which is spherical area of the FV subcell
       ! spherical polygon.
       do ie = 1,nelemd
          do j = 1,nf
             do i = 1,nf
                k = i+(j-1)*nf

                call gfr_f_ref_edges(nf, i, ae)
                call gfr_f_ref_edges(nf, j, be)
                do ai = 1,2
                   do bi = 1,2
                      p_sphere = ref2sphere(ae(ai), be(bi), elem(ie)%corners3D, cubed_sphere_map, &
                           elem(ie)%corners, elem(ie)%facenum)
                      fv_corners_xyz(ai,bi) = change_coordinates(p_sphere)
                   end do
                end do
                call sphere_tri_area(fv_corners_xyz(1,1), fv_corners_xyz(2,1), fv_corners_xyz(2,2), &
                     spherical_area)
                call sphere_tri_area(fv_corners_xyz(1,1), fv_corners_xyz(2,2), fv_corners_xyz(1,2), &
                     tmp)
                spherical_area = spherical_area + tmp
                gfr%fv_metdet(k,ie) = spherical_area/gfr%w_ff(k)
             end do
          end do
       end do
    end if

    ! Make the spherical area of the element according to FV and GLL agree to
    ! machine precision.
    do ie = 1,nelemd
       gfr%fv_metdet(:nf2,ie) = gfr%fv_metdet(:nf2,ie)* &
            (sum(elem(ie)%spheremp)/sum(gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie)))
    end do

    do ie = 1,nelemd
       do j = 1,nf
          call gfr_f_ref_center(nf, j, b)
          do i = 1,nf
             call gfr_f_ref_center(nf, i, a)

             gfr%spherep_f(i,j,ie) = ref2sphere(a, b, elem(ie)%corners3D, cubed_sphere_map, &
                  elem(ie)%corners, elem(ie)%facenum)

             call Dmap(wrk, a, b, elem(ie)%corners3D, cubed_sphere_map, elem(ie)%cartp, &
                  elem(ie)%facenum)

             det = wrk(1,1)*wrk(2,2) - wrk(1,2)*wrk(2,1)

             ! fv_metdet was obtained by remapping metdet. Make det(D)
             ! = fv_metdet.
             k = i+(j-1)*nf
             wrk = wrk*sqrt(gfr%fv_metdet(k,ie)/abs(det))
             det = gfr%fv_metdet(k,ie)

             gfr%D_f(i,j,:,:,ie) = wrk

             gfr%Dinv_f(i,j,1,1,ie) =  wrk(2,2)/det
             gfr%Dinv_f(i,j,1,2,ie) = -wrk(1,2)/det
             gfr%Dinv_f(i,j,2,1,ie) = -wrk(2,1)/det
             gfr%Dinv_f(i,j,2,2,ie) =  wrk(1,1)/det
          end do
       end do
    end do
  end subroutine gfr_init_geometry

  subroutine gfr_g2f_scalar(ie, gll_metdet, g, f) ! no gfr b/c public for testing
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), g(:,:,:)
    real(kind=real_kind), intent(out) :: f(:,:,:)

    integer :: k

    do k = 1, size(g,3)
       call gfr_g2f_remapd(gfr, gll_metdet, gfr%fv_metdet(:,ie), g(:,:,k), f(:,:,k))
    end do
  end subroutine gfr_g2f_scalar

  subroutine gfr_g2f_scalar_dp(gfr, ie, gll_metdet, dp_g, dp_f, g, f)
    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), dp_g(:,:,:), dp_f(:,:,:), g(:,:,:)
    real(kind=real_kind), intent(out) :: f(:,:,:)

    integer :: nf

    nf = gfr%nphys
    call gfr_g2f_scalar(ie, gll_metdet, dp_g*g, f)
    f = f(:nf,:nf,:)/dp_f(:nf,:nf,:)
  end subroutine gfr_g2f_scalar_dp

  subroutine gfr_g2f_vector(gfr, ie, elem, u_g, v_g, u_f, v_f)
    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(in) :: u_g(:,:,:), v_g(:,:,:)
    real(kind=real_kind), intent(out) :: u_f(:,:,:), v_f(:,:,:)

    real(kind=real_kind) :: wg(np,np,2), wf(np,np,2), ones(np*np), ones2(np,np)
    integer :: k, d, nf, nlev

    nf = gfr%nphys
    ones = one
    ones2 = one

    nlev = size(u_g,3)
    do k = 1, nlev
       ! sphere -> GLL ref
       do d = 1,2
          wg(:,:,d) = elem(ie)%Dinv(:,:,d,1)*u_g(:,:,k) + elem(ie)%Dinv(:,:,d,2)*v_g(:,:,k)
       end do
       do d = 1,2
          call gfr_g2f_remapd(gfr, ones2, ones, wg(:,:,d), wf(:,:,d))
       end do
       ! FV ref -> sphere
       do d = 1,2
          wg(:nf,:nf,d) = &
               gfr%D_f(:nf,:nf,d,1,ie)*wf(:nf,:nf,1) + &
               gfr%D_f(:nf,:nf,d,2,ie)*wf(:nf,:nf,2)
       end do
       u_f(:nf,:nf,k) = wg(:nf,:nf,1)
       v_f(:nf,:nf,k) = wg(:nf,:nf,2)
    end do
  end subroutine gfr_g2f_vector

  subroutine gfr_g2f_vector_dp(gfr, ie, elem, dp_g, dp_f, u_g, v_g, u_f, v_f)
    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(in) :: dp_g(:,:,:), dp_f(:,:,:), u_g(:,:,:), v_g(:,:,:)
    real(kind=real_kind), intent(out) :: u_f(:,:,:), v_f(:,:,:)

    integer :: nf

    nf = gfr%nphys
    call gfr_g2f_vector(gfr, ie, elem, dp_g*u_g, dp_g*v_g, u_f, v_f)
    u_f(:nf,:nf,:) = u_f(:nf,:nf,:)/dp_f(:nf,:nf,:)
    v_f(:nf,:nf,:) = v_f(:nf,:nf,:)/dp_f(:nf,:nf,:)
  end subroutine gfr_g2f_vector_dp

  subroutine gfr_g2f_mixing_ratio(gfr, ie, gll_metdet, dp_g, dp_f, qdp_g, q_f)
    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), dp_g(:,:,:), dp_f(:,:,:), &
         qdp_g(:,:,:)
    real(kind=real_kind), intent(out) :: q_f(:,:,:)

    real(kind=real_kind) :: qmin, qmax, wrk(np,np)
    integer :: q, k, nf, nf2

    nf = gfr%nphys
    nf2 = nf*nf
    do k = 1, size(qdp_g,3)
       call gfr_g2f_remapd(gfr, gll_metdet, gfr%fv_metdet(:,ie), qdp_g(:,:,k), q_f(:,:,k))
       q_f(:nf,:nf,k) = q_f(:nf,:nf,k)/dp_f(:nf,:nf,k)
       wrk = qdp_g(:,:,k)/dp_g(:,:,k)
       qmin = minval(wrk)
       qmax = maxval(wrk)
       wrk(:nf,:nf) = reshape(gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie), (/nf,nf/))
       call limiter_clip_and_sum(gfr%nphys, wrk, qmin, qmax, dp_f(:,:,k), q_f(:,:,k))
    end do
  end subroutine gfr_g2f_mixing_ratio

  subroutine gfr_g2f_mixing_ratios(gfr, ie, gll_metdet, dp_g, dp_f, qdp_g, q_f)
    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), dp_g(:,:,:), dp_f(:,:,:), &
         qdp_g(:,:,:,:)
    real(kind=real_kind), intent(out) :: q_f(:,:,:,:)

    real(kind=real_kind) :: qmin, qmax, wrk(np,np)
    integer :: q, k, nf

    nf = gfr%nphys
    do q = 1, size(qdp_g,4)
       call gfr_g2f_mixing_ratio(gfr, ie, gll_metdet, dp_g, dp_f, qdp_g(:,:,:,q), q_f(:,:,:,q))
    end do
  end subroutine gfr_g2f_mixing_ratios

  subroutine gfr_f2g_scalar(ie, gll_metdet, f, g) ! no gfr b/c public for testing
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), f(:,:,:)
    real(kind=real_kind), intent(out) :: g(:,:,:)

    integer :: k

    do k = 1, size(g,3)
       call gfr_f2g_remapd(gfr, gll_metdet, gfr%fv_metdet(:,ie), f(:,:,k), g(:,:,k))
    end do
  end subroutine gfr_f2g_scalar

  subroutine gfr_f2g_scalar_dp(gfr, ie, gll_metdet, dp_f, dp_g, f, g)
    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), dp_f(:,:,:), dp_g(:,:,:), f(:,:,:)
    real(kind=real_kind), intent(out) :: g(:,:,:)

    integer :: nf

    nf = gfr%nphys
    call gfr_f2g_scalar(ie, gll_metdet, dp_f(:nf,:nf,:)*f(:nf,:nf,:), g)
    g = g/dp_g
  end subroutine gfr_f2g_scalar_dp

  subroutine gfr_f2g_vector(gfr, ie, elem, u_f, v_f, u_g, v_g)
    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(in) :: u_f(:,:,:), v_f(:,:,:)
    real(kind=real_kind), intent(out) :: u_g(:,:,:), v_g(:,:,:)

    real(kind=real_kind) :: wg(np,np,2), wf(np,np,2), ones(np*np), ones2(np,np)
    integer :: k, d, nf, nlev

    nf = gfr%nphys
    ones = one
    ones2 = one

    nlev = size(u_g,3)
    do k = 1, nlev
       ! sphere -> FV ref
       do d = 1,2
          wf(:nf,:nf,d) = &
               gfr%Dinv_f(:nf,:nf,d,1,ie)*u_f(:nf,:nf,k) + &
               gfr%Dinv_f(:nf,:nf,d,2,ie)*v_f(:nf,:nf,k)
       end do
       do d = 1,2
          call gfr_f2g_remapd(gfr, ones2, ones, wf(:,:,d), wg(:,:,d))
       end do
       ! GLL ref -> sphere
       do d = 1,2
          wf(:,:,d) = elem(ie)%D(:,:,d,1)*wg(:,:,1) + elem(ie)%D(:,:,d,2)*wg(:,:,2)
       end do
       u_g(:,:,k) = wf(:,:,1)
       v_g(:,:,k) = wf(:,:,2)
    end do
  end subroutine gfr_f2g_vector

  subroutine gfr_f2g_vector_dp(gfr, ie, elem, dp_f, dp_g, u_f, v_f, u_g, v_g)
    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(in) :: dp_f(:,:,:), dp_g(:,:,:), u_f(:,:,:), v_f(:,:,:)
    real(kind=real_kind), intent(out) :: u_g(:,:,:), v_g(:,:,:)

    integer :: nf

    nf = gfr%nphys
    call gfr_f2g_vector(gfr, ie, elem, dp_f(:nf,:nf,:)*u_f(:nf,:nf,:), &
         dp_f(:nf,:nf,:)*v_f(:nf,:nf,:), u_g, v_g)
    u_g = u_g/dp_g
    v_g = v_g/dp_g
  end subroutine gfr_f2g_vector_dp

  subroutine gfr_f2g_mixing_ratios_he(hybrid, nets, nete, qmin, qmax)
    use viscosity_mod, only: neighbor_minmax
    use prim_advection_base, only: edgeAdvQminmax

    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: nets, nete
    real(kind=real_kind), intent(inout) :: qmin(:,:,:), qmax(:,:,:)

    call neighbor_minmax(hybrid, edgeAdvQminmax, nets, nete, qmin, qmax)
  end subroutine gfr_f2g_mixing_ratios_he

  subroutine gfr_f2g_dss(hybrid, elem, nets, nete)
    use dimensions_mod, only: nlev, qsize
    use edge_mod, only: edgevpack_nlyr, edgevunpack_nlyr, edge_g
    use bndry_mod, only: bndry_exchangev

    type (hybrid_t), intent(in) :: hybrid
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete

    integer :: ie, q, k

    do ie = nets, nete
       do q = 1,qsize
          do k = 1,nlev
             elem(ie)%derived%FQ(:,:,k,q) = elem(ie)%derived%FQ(:,:,k,q)*elem(ie)%spheremp(:,:)
          end do
       end do
       do q = 1,2
          do k = 1,nlev
             elem(ie)%derived%FM(:,:,q,k) = elem(ie)%derived%FM(:,:,q,k)*elem(ie)%spheremp(:,:)
          end do
       end do
       do k = 1,nlev
          elem(ie)%derived%FT(:,:,k) = elem(ie)%derived%FT(:,:,k)*elem(ie)%spheremp(:,:)
       end do
    end do
    do ie = nets, nete
       call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FQ, qsize*nlev, 0, qsize*nlev)
    end do
    call bndry_exchangeV(hybrid, edge_g)
    do ie = nets, nete
       call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FQ, qsize*nlev, 0, qsize*nlev)
    end do
    do ie = nets, nete
       call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FM, 2*nlev, 0, 3*nlev)
       call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FT, nlev, 2*nlev, 3*nlev)
    end do
    call bndry_exchangeV(hybrid, edge_g)
    do ie = nets, nete
       call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FM, 2*nlev, 0, 3*nlev)
       call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FT, nlev, 2*nlev, 3*nlev)
    end do
    do ie = nets, nete
       do q = 1,qsize
          do k = 1,nlev
             elem(ie)%derived%FQ(:,:,k,q) = elem(ie)%derived%FQ(:,:,k,q)*elem(ie)%rspheremp(:,:)
          end do
       end do
       do q = 1,2
          do k = 1,nlev
             elem(ie)%derived%FM(:,:,q,k) = elem(ie)%derived%FM(:,:,q,k)*elem(ie)%rspheremp(:,:)
          end do
       end do
       do k = 1,nlev
          elem(ie)%derived%FT(:,:,k) = elem(ie)%derived%FT(:,:,k)*elem(ie)%rspheremp(:,:)
       end do
    end do
  end subroutine gfr_f2g_dss

  subroutine gfr_g_make_nonnegative(gll_metdet, g)
    real(kind=real_kind), intent(in) :: gll_metdet(:,:)
    real(kind=real_kind), intent(inout) :: g(:,:,:)

    integer :: k, i, j
    real(kind=real_kind) :: nmass, spheremp(np,np), w(np,np)

    spheremp = gfr%w_gg*gll_metdet
    do k = 1, size(g,3)
       nmass = zero
       do j = 1,np
          do i = 1,np
             if (g(i,j,k) < zero) then
                nmass = nmass + spheremp(i,j)*g(i,j,k)
                g(i,j,k) = zero
                w(i,j) = zero
             else
                w(i,j) = spheremp(i,j)*g(i,j,k)
             end if
          end do
       end do
       if (nmass == zero) cycle
       w = (w/sum(w))/spheremp
       g(:,:,k) = g(:,:,k) + w*nmass
    end do
  end subroutine gfr_g_make_nonnegative

  subroutine gfr_f_get_latlon(ie, i, j, lat, lon)
    integer, intent(in) :: ie, i, j
    real(kind=real_kind), intent(out) :: lat, lon

    lat = gfr%spherep_f(i,j,ie)%lat
    lon = gfr%spherep_f(i,j,ie)%lon
  end subroutine gfr_f_get_latlon

  subroutine gfr_f_get_cartesian3d(ie, i, j, p)
    use coordinate_systems_mod, only: cartesian3D_t, change_coordinates

    integer, intent(in) :: ie, i, j
    type (cartesian3D_t), intent(out) :: p

    p = change_coordinates(gfr%spherep_f(i,j,ie))
  end subroutine gfr_f_get_cartesian3d

  ! d suffix means the inputs, outputs are densities.
  subroutine gfr_g2f_remapd(gfr, gll_metdet, fv_metdet, g, f)
    type (GllFvRemap_t), intent(in) :: gfr
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), fv_metdet(:), g(:,:)
    real(kind=real_kind), intent(out) :: f(:,:)

    integer :: nf, nf2, gi, gj, fi, fj, k
    real(kind=real_kind) :: gw(np,np)

    nf = gfr%nphys
    nf2 = nf*nf
    gw = g*gll_metdet
    do fj = 1,nf
       do fi = 1,nf
          k = fi + (fj-1)*nf
          f(fi,fj) = sum(gfr%g2f_remapd(:,:,k)*gw)/ &
               (gfr%w_ff(k)*fv_metdet(k))
       end do
    end do
  end subroutine gfr_g2f_remapd

  subroutine gfr_f2g_remapd(gfr, gll_metdet, fv_metdet, f, g)
    type (GllFvRemap_t), intent(in) :: gfr
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), fv_metdet(:), f(:,:)
    real(kind=real_kind), intent(out) :: g(:,:)

    integer :: nf, nf2, gi, gj, fi, fj
    real(kind=real_kind) :: wrk(np*np)

    nf = gfr%nphys
    nf2 = nf*nf
    wrk(:nf2) = reshape(f(:nf,:nf), (/nf2/))*fv_metdet(:nf2)
    do gj = 1,np
       do gi = 1,np
          g(gi,gj) = sum(gfr%f2g_remapd(:nf2,gi,gj)*wrk(:nf2))/ &
               gll_metdet(gi,gj)
       end do
    end do
  end subroutine gfr_f2g_remapd

  subroutine limiter_clip_and_sum(n, spheremp, qmin, qmax, dp, q)
    integer, intent(in) :: n
    real (kind=real_kind), intent(in) :: spheremp(:,:), dp(:,:)
    real (kind=real_kind), intent(inout) :: qmin, qmax, q(:,:)

    integer :: k1, i, j
    logical :: modified
    real(kind=real_kind) :: addmass, mass, sumc, den
    real(kind=real_kind) :: c(n*n), v(n*n), x(n*n)

    x = reshape(q(:n,:n), (/n*n/))
    c = reshape(spheremp(:n,:n)*dp(:n,:n), (/n*n/))

    sumc = sum(c)
    mass = sum(c*x)
    ! In the case of an infeasible problem, prefer to conserve mass
    ! and violate a bound.
    if (mass < qmin*sumc) qmin = mass / sumc
    if (mass > qmax*sumc) qmax = mass / sumc

    addmass = zero

    ! Clip.
    modified = .false.
    do k1 = 1, n*n
       if (x(k1) > qmax) then
          modified = .true.
          addmass = addmass + (x(k1) - qmax)*c(k1)
          x(k1) = qmax
       elseif (x(k1) < qmin) then
          modified = .true.
          addmass = addmass + (x(k1) - qmin)*c(k1)
          x(k1) = qmin
       end if
    end do
    if (.not. modified) return

    if (addmass /= zero) then
       ! Determine weights.
       if (addmass > zero) then
          v = qmax - x
       else
          v = x - qmin
       end if
       den = sum(v*c)
       if (den > zero) x = x + addmass*(v/den)
    end if

    q(:n,:n) = reshape(x, (/n,n/))
  end subroutine limiter_clip_and_sum

  ! --- testing ---
  ! Everything below is for testing.

  subroutine set_ps_Q(elem, nets, nete, timeidx, qidx, nlev)
    use coordinate_systems_mod, only: cartesian3D_t, change_coordinates

    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete, timeidx, qidx, nlev

    integer :: ie, i, j, k
    type (cartesian3D_t) :: p
    real(kind=real_kind) :: q

    do ie = nets, nete
       do j = 1,np
          do i = 1,np
             p = change_coordinates(elem(ie)%spherep(i,j))
             elem(ie)%state%ps_v(i,j,timeidx) = &
                  1.0d3*(1 + 0.05*sin(2*p%x+0.5)*sin(p%y+1.5)*sin(3*p%z+2.5))
             q = 0.5*(1 + sin(3*p%x)*sin(3*p%y)*sin(4*p%z))
             do k = 1,nlev
                elem(ie)%state%Q(i,j,k,qidx) = q
             end do
          end do
       end do
    end do
  end subroutine set_ps_Q

  subroutine check_g2f_mixing_ratio(gfr, hybrid, ie, qi, elem, dp, dp_fv, q_g, q_f)
    use kinds, only: iulog

    type (GllFvRemap_t), intent(in) :: gfr
    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: ie, qi
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(in) :: dp(:,:,:), dp_fv(:,:,:), q_g(:,:,:), q_f(:,:,:)

    real(kind=real_kind) :: qmin_f, qmin_g, qmax_f, qmax_g, mass_f, mass_g, den
    integer :: q, k, nf, nf2

    nf = gfr%nphys
    nf2 = nf*nf
    do k = 1,size(dp,3)
       qmin_f = minval(q_f(:nf,:nf,k))
       qmax_f = maxval(q_f(:nf,:nf,k))
       qmin_g = minval(elem(ie)%state%Q(:,:,k,qi))
       qmax_g = maxval(elem(ie)%state%Q(:,:,k,qi))
       den = gfr%tolfac*max(1e-10, maxval(abs(elem(ie)%state%Q(:,:,k,qi))))
       mass_f = sum(reshape(gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie), (/nf,nf/))* &
            dp_fv(:nf,:nf,k)*q_f(:nf,:nf,k))
       mass_g = sum(elem(ie)%spheremp*dp(:,:,k)*q_g(:,:,k))
       if (qmin_f < qmin_g - 10*eps*den .or. qmax_f > qmax_g + 10*eps*den) then
          write(iulog,*) 'gfr> g2f mixing ratio limits:', hybrid%par%rank, hybrid%ithr, ie, qi, k, &
               qmin_g, qmin_f-qmin_g, qmax_f-qmax_g, qmax_g, mass_f, mass_g
       end if
       if (abs(mass_f - mass_g) > gfr%tolfac*20*eps*max(mass_f, mass_g)) then
          write(iulog,*) 'gfr> g2f mixing ratio mass:', hybrid%par%rank, hybrid%ithr, ie, qi, k, &
               qmin_g, qmax_g, mass_f, mass_g
       end if
    end do
  end subroutine check_g2f_mixing_ratio

  subroutine check_f2g_mixing_ratio(gfr, hybrid, ie, qi, elem, qmin, qmax, dp, q0_g, q1_g)
    use kinds, only: iulog

    type (GllFvRemap_t), intent(in) :: gfr
    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: ie, qi
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(in) :: qmin(:), qmax(:), dp(:,:,:), q0_g(:,:,:), q1_g(:,:,:)

    real(kind=real_kind) :: qmin_f, qmin_g, qmax_f, qmax_g, mass_f, mass0, mass1, den, &
         wr(np,np)
    integer :: q, k

    do k = 1,size(dp,3)
       qmin_f = qmin(k)
       qmax_f = qmax(k)
       qmin_g = minval(q1_g(:,:,k))
       qmax_g = maxval(q1_g(:,:,k))
       den = gfr%tolfac*max(1e-10, maxval(abs(q0_g(:,:,k))))
       if (qmin_g < qmin_f - 50*eps*den .or. qmax_g > qmax_f + 50*eps*den) then
          write(iulog,*) 'gfr> f2g mixing ratio limits:', hybrid%par%rank, hybrid%ithr, ie, qi, k, &
               qmin_f, qmin_g-qmin_f, qmax_g-qmax_f, qmax_f, mass0, mass1
       end if
       mass0 = sum(elem(ie)%spheremp*dp(:,:,k)*q0_g(:,:,k))
       mass1 = sum(elem(ie)%spheremp*dp(:,:,k)*q1_g(:,:,k))
       den = sum(elem(ie)%spheremp*dp(:,:,k)*maxval(abs(q0_g(:,:,k))))
       if (abs(mass1 - mass0) > gfr%tolfac*20*eps*den) then
          write(iulog,*) 'gfr> f2g mixing ratio mass:', hybrid%par%rank, hybrid%ithr, ie, qi, k, &
               qmin_f, qmin_g, qmax_g, qmax_f, mass0, mass1
       end if
    end do
  end subroutine check_f2g_mixing_ratio
  
  subroutine check_nonnegative(elem, nets, nete)
    ! Check gfr_g_make_nonnegative.

    use kinds, only: iulog

    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete

    real(kind=real_kind) :: wrk3(np,np,1), mass0, mass1, rd
    integer :: ie, i, j, sign

    do ie = nets,nete
       sign = 1
       do j = 1,np
          do i = 1,np
             wrk3(i,j,1) = one + sign*(one + cos(real(i,real_kind)))*j
             sign = -sign
          end do
       end do
       mass0 = sum(elem(ie)%spheremp*wrk3(:,:,1))
       call gfr_g_make_nonnegative(elem(ie)%metdet, wrk3)
       mass1 = sum(elem(ie)%spheremp*wrk3(:,:,1))
       rd = (mass1 - mass0)/mass0
       if (rd /= rd .or. rd > 20*eps .or. any(wrk3(:,:,1) < zero)) then
          write(iulog,*) 'gfr> nonnegative', ie, rd, mass0, mass1, wrk3(:,:,1), 'ERROR'
       end if
    end do
  end subroutine check_nonnegative

  subroutine check(par, dom_mt, gfr, elem, verbose)
    use kinds, only: iulog
    use parallel_mod, only: parallel_t
    use dimensions_mod, only: nlev, qsize
    use domain_mod, only: domain1d_t
    use edge_mod, only: edge_g, edgevpack_nlyr, edgevunpack_nlyr
    use bndry_mod, only: bndry_exchangev
    use viscosity_mod, only: neighbor_minmax
    use parallel_mod, only: global_shared_buf, global_shared_sum
    use global_norms_mod, only: wrap_repro_sum
    use reduction_mod, only: ParallelMin, ParallelMax
    use prim_advection_base, only: edgeAdvQminmax

    type (parallel_t), intent(in) :: par
    type (domain1d_t), intent(in) :: dom_mt(:)
    type (GllFvRemap_t), intent(in) :: gfr
    type (element_t), intent(inout) :: elem(:)
    logical, intent(in) :: verbose

    real(kind=real_kind) :: a, b, rd, x, y, f0(np,np), f1(np,np), g(np,np), &
         wrk(np,np), qmin, qmax, qmin1, qmax1
    integer :: nf, nf2, ie, i, j, iremap, info, ilimit
    real(kind=real_kind), allocatable :: Qdp_fv(:,:,:), ps_v_fv(:,:,:), &
         qmins(:,:,:), qmaxs(:,:,:)
    logical :: limit
    character(32) :: msg

    ! Purposely construct our own hybrid object to test gfr_hybrid_create.
    type (hybrid_t) :: hybrid
    integer :: nets, nete

    nf = gfr%nphys
    nf2 = nf*nf

    call gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)

    if (hybrid%masterthread) then
       write(iulog,  '(a,i3,a,i3)'), 'gfr> npi', gfr%npi, ' nphys', nf
       if (verbose) then
          write(iulog,*) 'gfr> w_ff', nf, gfr%w_ff(:nf2)
          write(iulog,*) 'gfr> w_gg', np, gfr%w_gg(:np, :np)
          write(iulog,*) 'gfr> w_sgsg', gfr%npi, gfr%w_sgsg(:gfr%npi, :gfr%npi)
          write(iulog,*) 'gfr> M_gf', np, nf, gfr%M_gf(:np, :np, :nf, :nf)
          write(iulog,*) 'gfr> M_sgf', gfr%npi, nf, gfr%M_sgf(:gfr%npi, :gfr%npi, :nf, :nf)
          write(iulog,*) 'gfr> interp', gfr%npi, np, gfr%interp(:gfr%npi, :gfr%npi, :np, :np)
          write(iulog,*) 'gfr> f2g_remapd', np, nf, gfr%f2g_remapd(:nf*nf,:,:)
       end if
    end if

    ! Cell-local correctness checks
    do ie = nets, nete
       ! Check that areas match.
       a = sum(elem(ie)%metdet * gfr%w_gg)
       b = sum(gfr%fv_metdet(:nf2,ie) * gfr%w_ff(:nf2))
       rd = abs(b - a)/abs(a)
       if (rd /= rd .or. rd > 10*eps) write(iulog,*) 'gfr> area', ie, a, b, rd

       ! Check FV geometry.
       f0(:nf,:nf) = gfr%D_f(:,:,1,1,ie)*gfr%D_f(:,:,2,2,ie) - &
            gfr%D_f(:,:,1,2,ie)*gfr%D_f(:,:,2,1,ie)
       rd = maxval(reshape(abs(f0(:nf,:nf)), (/nf2/)) - gfr%fv_metdet(:nf2,ie))/ &
            maxval(gfr%fv_metdet(:nf2,ie))
       if (rd > 10*eps) write(iulog,*) 'gfr> D', ie, rd
       f0(:nf,:nf) = gfr%Dinv_f(:,:,1,1,ie)*gfr%Dinv_f(:,:,2,2,ie) - &
            gfr%Dinv_f(:,:,1,2,ie)*gfr%Dinv_f(:,:,2,1,ie)
       rd = maxval(reshape(abs(f0(:nf,:nf)), (/nf2/)) - one/gfr%fv_metdet(:nf2,ie))/ &
            maxval(one/gfr%fv_metdet(:nf2,ie))
       if (rd > 10*eps) write(iulog,*) 'gfr> Dinv', ie, rd

       ! Check that FV -> GLL -> FV recovers the original FV values exactly
       ! (with no DSS and no limiter).
       do j = 1,nf
          x = real(j-1, real_kind)/real(nf, real_kind)
          do i = 1,nf
             y = real(i-1, real_kind)/real(nf, real_kind)
             f0(i,j) = real(ie)/nelemd + x*x + ie*x + cos(ie + 4.2*y)
          end do
       end do
       call gfr_f2g_remapd(gfr, elem(ie)%metdet, gfr%fv_metdet(:,ie), f0, g)
       call gfr_g2f_remapd(gfr, elem(ie)%metdet, gfr%fv_metdet(:,ie), g, f1)
       wrk(:nf,:nf) = reshape(gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie), (/nf,nf/))
       a = sum(wrk(:nf,:nf)*abs(f1(:nf,:nf) - f0(:nf,:nf)))
       b = sum(wrk(:nf,:nf)*abs(f0(:nf,:nf)))
       rd = a/b
       if (rd /= rd .or. rd > 10*eps) &
            write(iulog,*) 'gfr> recover', ie, a, b, rd, gfr%fv_metdet(:nf2,ie)
    end do
    call check_nonnegative(elem, nets, nete)

    ! For convergence testing. Run this testing routine with a sequence of ne
    ! values and plot log l2 error vs log ne.
    allocate(Qdp_fv(gfr%nphys, gfr%nphys, nets:nete), ps_v_fv(gfr%nphys, gfr%nphys, nets:nete))
    allocate(qmins(nlev,qsize,nets:nete), qmaxs(nlev,qsize,nets:nete))
    do ilimit = 0,1
       limit = ilimit > 0
       ! 0. Create synthetic q and ps_v.
       call set_ps_Q(elem, nets, nete, 1, 1, nlev)
       call set_ps_Q(elem, nets, nete, 2, 2, nlev)
       do iremap = 1,1
          ! 1. GLL -> FV
          do ie = nets, nete
             call gfr_g2f_remapd(gfr, elem(ie)%metdet, gfr%fv_metdet(:,ie), &
                  elem(ie)%state%ps_v(:,:,1)*elem(ie)%state%Q(:,:,1,1), Qdp_fv(:,:,ie))
             if (limit) then
                call gfr_g2f_remapd(gfr, elem(ie)%metdet, gfr%fv_metdet(:,ie), &
                     elem(ie)%state%ps_v(:,:,1), ps_v_fv(:,:,ie))
                qmin = minval(elem(ie)%state%Q(:,:,1,1))
                qmax = maxval(elem(ie)%state%Q(:,:,1,1))
                wrk(:nf,:nf) = Qdp_fv(:nf,:nf,ie)/ps_v_fv(:nf,:nf,ie)
                f0(:nf,:nf) = reshape(gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie), (/nf,nf/))
                call limiter_clip_and_sum(nf, f0, qmin, qmax, ps_v_fv(:,:,ie), wrk)
                Qdp_fv(:nf,:nf,ie) = wrk(:nf,:nf)*ps_v_fv(:nf,:nf,ie)
             end if
          end do
          ! 2. FV -> GLL
          if (limit) then
             ! 2a. Get q bounds
             do ie = nets, nete
                wrk(:nf,:nf) = Qdp_fv(:nf,:nf,ie)/ps_v_fv(:nf,:nf,ie)
                qmins(:,:,ie) = minval(wrk(:nf,:nf))
                qmaxs(:,:,ie) = maxval(wrk(:nf,:nf))
             end do
             ! 2b. Halo exchange q bounds.
             call neighbor_minmax(hybrid, edgeAdvQminmax, nets, nete, qmins, qmaxs)
          endif
          ! 2c. Remap
          do ie = nets, nete
             wrk = elem(ie)%state%Q(:,:,1,1)
             call gfr_f2g_remapd(gfr, elem(ie)%metdet, gfr%fv_metdet(:,ie), &
                  Qdp_fv(:,:,ie), elem(ie)%state%Q(:,:,1,1))
             elem(ie)%state%Q(:,:,1,1) = elem(ie)%state%Q(:,:,1,1)/elem(ie)%state%ps_v(:,:,1)
             if (limit) then
                qmins(1,1,ie) = min(qmins(1,1,ie), minval(wrk))
                qmaxs(1,1,ie) = max(qmaxs(1,1,ie), maxval(wrk))
                call limiter_clip_and_sum(np, elem(ie)%spheremp, & ! same as w_gg*gll_metdet
                     qmins(1,1,ie), qmaxs(1,1,ie), elem(ie)%state%ps_v(:,:,1), &
                     elem(ie)%state%Q(:,:,1,1))
             end if
          end do
          ! 3. DSS
          do ie = nets, nete
             elem(ie)%state%Q(:,:,1,1) = &
                  elem(ie)%state%ps_v(:,:,1)*elem(ie)%state%Q(:,:,1,1)*elem(ie)%spheremp(:,:)
             call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Q(:,:,1,1), 1, 0, 1)
          end do
          call bndry_exchangeV(hybrid, edge_g)
          do ie = nets, nete
             call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Q(:,:,1,1), 1, 0, 1)
             elem(ie)%state%Q(:,:,1,1) = &
                  (elem(ie)%state%Q(:,:,1,1)*elem(ie)%rspheremp(:,:))/elem(ie)%state%ps_v(:,:,1)
          end do
       end do
       ! 5. Compute error.
       qmin = two
       qmax = -two
       qmin1 = two
       qmax1 = -two
       do ie = nets, nete
          wrk = gfr%w_gg(:,:)*elem(ie)%metdet(:,:)
          ! L2 on q. Might switch to q*ps_v.
          global_shared_buf(ie,1) = &
               sum(wrk*(elem(ie)%state%Q(:,:,1,1) - elem(ie)%state%Q(:,:,1,2))**2)
          global_shared_buf(ie,2) = &
               sum(wrk*elem(ie)%state%Q(:,:,1,2)**2)
          ! Mass conservation.
          wrk = wrk*elem(ie)%state%ps_v(:,:,1)
          global_shared_buf(ie,3) = sum(wrk*elem(ie)%state%Q(:,:,1,2))
          global_shared_buf(ie,4) = sum(wrk*elem(ie)%state%Q(:,:,1,1))
          qmin = min(qmin, minval(elem(ie)%state%Q(:,:,1,1)))
          qmin1 = min(qmin1, minval(elem(ie)%state%Q(:,:,1,2)))
          qmax = max(qmax, maxval(elem(ie)%state%Q(:,:,1,1)))
          qmax1 = max(qmax1, maxval(elem(ie)%state%Q(:,:,1,2)))
       end do
       call wrap_repro_sum(nvars=4, comm=hybrid%par%comm)
       qmin = ParallelMin(qmin, hybrid)
       qmax = ParallelMax(qmax, hybrid)
       qmin1 = ParallelMin(qmin1, hybrid)
       qmax1 = ParallelMax(qmax1, hybrid)
       if (hybrid%masterthread) then
          write(iulog, '(a,i3)') 'gfr> limiter', ilimit
          rd = sqrt(global_shared_sum(1)/global_shared_sum(2))
          write(iulog, '(a,es12.4)') 'gfr> l2  ', rd
          rd = abs(global_shared_sum(4) - global_shared_sum(3))/global_shared_sum(3)
          msg = ''
          if (rd > 10*eps) msg = ' ERROR'
          write(iulog, '(a,es11.3,a8)') 'gfr> mass', rd, msg
          msg = ''
          if (limit .and. (qmin < qmin1 - 5*eps .or. qmax > qmax1 + 5*eps)) msg = ' ERROR'
          write(iulog, '(a,es11.3,es11.3,a8)') 'gfr> limit', min(zero, qmin - qmin1), &
               max(zero, qmax - qmax1), msg
       end if
    end do
    deallocate(Qdp_fv, ps_v_fv, qmins, qmaxs)
  end subroutine check

  subroutine gfr_test(hybrid, dom_mt, hvcoord, deriv, elem)
    use domain_mod, only: domain1d_t
    use derivative_mod, only: derivative_t
    use hybvcoord_mod, only: hvcoord_t

    type (hybrid_t), intent(in) :: hybrid
    type (domain1d_t), intent(in) :: dom_mt(:)
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(inout) :: elem(:)
    type (hvcoord_t) , intent(in) :: hvcoord

    integer :: nphys

    do nphys = 1, np
       ! This is meant to be called before threading starts.
       if (hybrid%ithr == 0) call gfr_init(hybrid%par, elem, nphys)
       !$omp barrier

       call check(hybrid%par, dom_mt, gfr, elem, .false.)

       ! This is meant to be called after threading ends.
       !$omp barrier
       if (hybrid%ithr == 0) call gfr_finish()
       !$omp barrier
    end do
  end subroutine gfr_test
end module gllfvremap_mod
