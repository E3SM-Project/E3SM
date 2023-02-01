#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module gllfvremap_mod
  ! High-order, mass-conserving, optionally shape-preserving
  !     FV physics <-> GLL dynamics
  ! remap.
  !
  ! This module implements algorithms to remap state variables and
  ! tendencies between the spectral element (SE) dynamics grid and the
  ! finite volume (FV) physics grid. The two grids are the same at the
  ! element level. The dynamics grid has np x np GLL grid points
  ! inside each element. The physics grid has nphys x nphys FV
  ! subcells inside each element. We frequently shorten nphys to nf or
  ! just N in the pgN physics grid specification. Thus, the dynamics
  ! grid specification may be written ne30np4 and the physics grid one
  ! as ne30pg2.
  !   The remap algorithms have the following properties:
  !     * Conserve mass.
  !     * Have order of accuracy (OOA) N for N in pgN without a limiter.
  !     * Optionally prevent new mixing ratio (Q) extrema by a
  !       nonlinear limiter. OOA is then min(N,2). The extrema are
  !       limited in the full state, not just the tendency, so Q state
  !       or tendency returned from the physics is assured to be
  !       physically valid.
  !     * pg1's OOA is boosted to ~1.6.
  !   This module works with cubed_sphere_map 0 and 2, but only 2 will
  ! work in the fully coupled E3SM.
  !   This module supports all ftype values 0 to 4, as the calculations
  ! in this module are independent of ftype.
  !   We find in practice that pg1 is too coarse. To see this, run the
  ! Homme-standalone test dcmip2016_test1_pg1 and compare results with
  ! dcmip2016_test1_pg2 and dcmip2016_test1 (np4). pg2 and np4 fields are
  ! nearly identical out to day 30, whereas pg1 fields differ visibly.
  !
  ! AMB 2019/07-2020/06 Initial

  use hybrid_mod, only: hybrid_t
  use kinds, only: real_kind
  use dimensions_mod, only: np, npsq, qsize, nelemd
  use element_mod, only: element_t
  use coordinate_systems_mod, only: cartesian3D_t
#ifdef HOMMEXX_BFB_TESTING
  use bfb_mod, only: bfb_pow
#endif

  implicit none

  private

  ! Main API.
  public :: &
       ! Initialize this module.
       gfr_init, &
       ! Clean up this module.
       gfr_finish, &
       ! Remap GLL state to FV grid.
       gfr_dyn_to_fv_phys, &
       ! Remap FV tendencies or state to GLL grid.
       gfr_fv_phys_to_dyn, &
       ! Remap phis.
       gfr_dyn_to_fv_phys_topo, &
       gfr_fv_phys_to_dyn_topo, &
       gfr_dyn_to_fv_phys_topo_data, &
       ! If nphys == 1, reconstruct the field to boost the OOA. If
       ! nphys > 1, returns immediately.
       gfr_pg1_reconstruct_topo, & ! call after the gfr_fv_phys_to_dyn_topo and the DSS
       gfr_pg1_reconstruct         ! call after gfr_fv_phys_to_dyn and the DSS

  integer, parameter :: nphys_max = np

  real(kind=real_kind), parameter :: &
       zero = 0.0_real_kind, half = 0.5_real_kind, &
       one = 1.0_real_kind, two = 2.0_real_kind, &
       four = 4.0_real_kind, eps = epsilon(1.0_real_kind)

  ! Type for special case of pg1.
  type, private :: Pg1SolverData_t
     real(kind=real_kind) :: Achol(np*np,np*np), B(4,np*np), s(np*np), sts
  end type Pg1SolverData_t

  ! Top-level data type and functions for high-order, shape-preserving FV <->
  ! GLL remap.
  type, private :: GllFvRemap_t
     integer :: nphys, npi, check
     logical :: have_fv_topo_file_phis, boost_pg1, check_ok
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
     real(kind=real_kind), allocatable :: &
          ! FV subcell areas; FV analogue of GLL elem(ie)%metdet arrays
          fv_metdet(:,:), & ! (nphys*nphys,nelemd)
          ! Vector on ref elem -> vector on sphere
          D_f(:,:,:,:), & ! (nphys*nphys,2,2,nelemd)
          ! Inverse of D_f
          Dinv_f(:,:,:,:), &
          qmin(:,:,:), qmax(:,:,:), &
          phis(:,:), &
          ! For 'check', when it's on
          check_areas(:,:)
     type (cartesian3D_t), allocatable :: &
          center_f(:,:,:), & ! (nphys,nphys,nelemd)
          corners_f(:,:,:,:) ! (4,nphys,nphys,nelemd)
     type (Pg1SolverData_t) :: pg1sd
  end type GllFvRemap_t

  type (GllFvRemap_t), private :: gfr

  ! For testing in gllfvremap_util_mod and gllfvremap_ut.
  public :: &
       gfr_test, &
       gfr_g2f_scalar, gfr_f2g_scalar, gfr_g2f_vector, &
       gfr_f_get_area, gfr_f_get_latlon, gfr_f_get_corner_latlon, gfr_f_get_cartesian3d, &
       gfr_g_make_nonnegative, gfr_dyn_to_fv_phys_topo_elem, gfr_f2g_dss
  ! For testing in gllfvremap_ut.
  public :: &
       limiter1_clip_and_sum, calc_dp_fv, gfr_get_nphys
  ! For C++ dycore.
  public :: &
       gfr_init_hxx

  ! Interfaces to support calling inside or outside a horizontally
  ! threaded region.
  interface gfr_dyn_to_fv_phys
     module procedure gfr_dyn_to_fv_phys_hybrid
     module procedure gfr_dyn_to_fv_phys_dom_mt
  end interface gfr_dyn_to_fv_phys

  interface gfr_fv_phys_to_dyn
     module procedure gfr_fv_phys_to_dyn_hybrid
     module procedure gfr_fv_phys_to_dyn_dom_mt
  end interface gfr_fv_phys_to_dyn

  interface gfr_dyn_to_fv_phys_topo
     module procedure gfr_dyn_to_fv_phys_topo_hybrid
     module procedure gfr_dyn_to_fv_phys_topo_dom_mt
     module procedure gfr_dyn_to_fv_phys_topo_mpi_only
  end interface gfr_dyn_to_fv_phys_topo

  interface gfr_fv_phys_to_dyn_topo
     module procedure gfr_fv_phys_to_dyn_topo_hybrid
     module procedure gfr_fv_phys_to_dyn_topo_dom_mt
     module procedure gfr_fv_phys_to_dyn_topo_mpi_only
  end interface gfr_fv_phys_to_dyn_topo

  interface gfr_pg1_reconstruct
     module procedure gfr_pg1_reconstruct_hybrid
     module procedure gfr_pg1_reconstruct_dom_mt
  end interface gfr_pg1_reconstruct

  interface gfr_pg1_reconstruct_topo
     module procedure gfr_pg1_reconstruct_topo_hybrid
     module procedure gfr_pg1_reconstruct_topo_dom_mt
  end interface gfr_pg1_reconstruct_topo

contains

  ! ----------------------------------------------------------------------
  ! Public API.

  subroutine gfr_init(par, elem, nphys, check, boost_pg1)
    ! Initialize the gfr internal data structure.
    !   nphys is N in pgN.
    !   check is optional and defaults to 0, no checking. It will produce very
    ! verbose output if something goes wrong. It is also expensive. It is
    ! intended to be used in unit testing and (if ever needed) for
    ! debugging. There are three levels: 0, no checking; 1, global properties
    ! only; 2, also element-local properties.

    use kinds, only: iulog
    use dimensions_mod, only: nlev
    use parallel_mod, only: parallel_t, abortmp
    use quadrature_mod, only : gausslobatto, quadrature_t

    type (parallel_t), intent(in) :: par
    type (element_t), intent(in) :: elem(:)
    integer, intent(in) :: nphys
    integer, intent(in), optional :: check
    logical, intent(in), optional :: boost_pg1

    real(real_kind) :: R(npsq,nphys_max*nphys_max), tau(npsq)
    integer :: nphys2

    gfr%check = 0
    if (present(check)) gfr%check = check
    gfr%check_ok = .true.

    gfr%boost_pg1 = .false.
    if (nphys == 1 .and. present(boost_pg1)) gfr%boost_pg1 = boost_pg1    

    gfr%tolfac = one
    if (par%masterproc) then
       write(iulog,*) 'gfr> Running with dynamics and physics on separate grids (physgrid).'
       write(iulog, '(a,i3,a,i2,a,l2)') 'gfr> init nphys', nphys, ' check', gfr%check, &
            ' boost_pg1', gfr%boost_pg1
       if (nphys == 1) then
          ! Document state of pg1. dcmip2016_test1 shows it is too coarse. For
          ! boost_pg1 = true, stepon's DSS loop needs to be separated from its
          ! tendency application loop.
          write(iulog,*) 'gfr> Warning: pg1 is too coarse; see comments at top of gllfvremap_mod.F90'
          if (.not. gfr%boost_pg1) then
             write(iulog,*) 'gfr> Warning: If you want to try pg1, use the boosted-accuracy &
                  &boost_pg1 option and call gfr_pg1_reconstruct(_topo).'
          end if
       end if
    end if

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
    ! npi is the internal GLL np parameter. The high-order remap operator remaps
    ! from FV to npi-GLL grids, then interpolates from npi-GLL to np-GLL
    ! grids. In the case of nphys=1, npi must be 2 for GLL to make sense.
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
         gfr%D_f(nphys2,2,2,nelemd), gfr%Dinv_f(nphys2,2,2,nelemd), &
         gfr%qmin(nlev,max(1,qsize),nelemd), gfr%qmax(nlev,max(1,qsize),nelemd), &
         gfr%phis(nphys2,nelemd), gfr%center_f(nphys,nphys,nelemd), &
         gfr%corners_f(4,nphys,nphys,nelemd))
    call gfr_init_geometry(elem, gfr)
    call gfr_init_Dmap(elem, gfr)

    if (nphys == 1 .and. gfr%boost_pg1) call gfr_pg1_init(gfr)

    if (gfr%check > 0) call check_areas(par, gfr, elem, 1, nelemd)
  end subroutine gfr_init

  subroutine gfr_init_hxx() bind(c)
#if KOKKOS_TARGET
    use control_mod, only: theta_hydrostatic_mode
    use iso_c_binding, only: c_bool
    interface
       subroutine init_gllfvremap_c(nelemd, np, nf, nf_max, theta_hydrostatic_mode, &
            fv_metdet, g2f_remapd, f2g_remapd, D_f, Dinv_f) bind(c)
         use iso_c_binding, only: c_bool, c_int, c_double
         integer (c_int), value, intent(in) :: nelemd, np, nf, nf_max
         logical (c_bool), value, intent(in) :: theta_hydrostatic_mode
         real (c_double), dimension(nf*nf,nelemd), intent(in) :: fv_metdet
         real (c_double), dimension(np,np,nf_max*nf_max), intent(in) :: g2f_remapd
         real (c_double), dimension(nf_max*nf_max,np,np), intent(in) :: f2g_remapd
         real (c_double), dimension(nf*nf,2,2,nelemd), intent(in) :: D_f, Dinv_f
       end subroutine init_gllfvremap_c
    end interface
    logical (c_bool) :: thm
    thm = theta_hydrostatic_mode
    call init_gllfvremap_c(nelemd, np, gfr%nphys, nphys_max, thm, &
         gfr%fv_metdet, gfr%g2f_remapd, gfr%f2g_remapd, gfr%D_f, gfr%Dinv_f)
#endif
  end subroutine gfr_init_hxx

  function gfr_get_nphys() result(nf)
    integer :: nf
    nf = gfr%nphys
  end function gfr_get_nphys

  subroutine gfr_finish()
    ! Deallocate the internal gfr structure.

    if (.not. allocated(gfr%fv_metdet)) return
    deallocate(gfr%fv_metdet, gfr%D_f, gfr%Dinv_f, gfr%qmin, gfr%qmax, gfr%phis, &
         gfr%center_f, gfr%corners_f)
  end subroutine gfr_finish

  subroutine gfr_dyn_to_fv_phys_hybrid(hybrid, nt, hvcoord, elem, nets, nete, &
       ps, phis, T, uv, omega_p, q)
    ! Remap ps, phis, T, uv, omega_p, q from GLL to FV grids.

    use dimensions_mod, only: nlev
    use hybvcoord_mod, only: hvcoord_t
    use element_ops, only: get_temperature, get_field
    use physical_constants, only: p0, kappa

    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: nt
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(in) :: elem(:)
    integer, intent(in) :: nets, nete
    real(kind=real_kind), intent(out) :: ps(:,:), phis(:,:), T(:,:,:), &
         uv(:,:,:,:), omega_p(:,:,:), q(:,:,:,:)

    real(kind=real_kind), dimension(np,np,nlev) :: wg1, dp, p
    real(kind=real_kind), dimension(np*np,nlev) :: wf1, dp_fv, p_fv
    real(kind=real_kind) :: qmin, qmax, ones(np,np)
    integer :: ie, nf, nf2, qi, qsize, k, nerr

    ones = one
    nf = gfr%nphys
    nf2 = nf*nf

    qsize = size(q,3)

    do ie = nets,nete
       call gfr_g2f_scalar(ie, elem(ie)%metdet, elem(ie)%state%ps_v(:,:,nt:nt), wf1(:,:1))
       ps(:nf2,ie) = wf1(:nf2,1)
       dp = elem(ie)%state%dp3d(:,:,:,nt)
       call calc_dp_fv(nf, hvcoord, ps(:,ie), dp_fv)

       if (gfr%have_fv_topo_file_phis) then
          phis(:nf2,ie) = gfr%phis(:,ie)
       else
          call gfr_g2f_scalar_and_limit(gfr, ie, elem(ie)%metdet, elem(ie)%state%phis, &
               phis(:,ie))
       end if

       call get_temperature(elem(ie), wg1, hvcoord, nt)
       call get_field(elem(ie), 'p', p, hvcoord, nt, -1)
       call gfr_g2f_scalar(ie, elem(ie)%metdet, p, p_fv)
#ifndef HOMMEXX_BFB_TESTING
       wg1 = wg1*(p0/p)**kappa
       call gfr_g2f_scalar_dp(gfr, ie, elem(ie)%metdet, dp, dp_fv, wg1, wf1)
       T(:nf2,:,ie) = wf1(:nf2,:)*(p_fv(:nf2,:)/p0)**kappa
#else
       wg1 = wg1*bfb_pow(p0/p, kappa)
       call gfr_g2f_scalar_dp(gfr, ie, elem(ie)%metdet, dp, dp_fv, wg1, wf1)
       T(:nf2,:,ie) = wf1(:nf2,:)*bfb_pow(p_fv(:nf2,:)/p0, kappa)
#endif
       
       call gfr_g2f_vector(ie, elem, &
            elem(ie)%state%v(:,:,1,:,nt), elem(ie)%state%v(:,:,2,:,nt), &
            uv(:,1,:,ie), uv(:,2,:,ie))

       call get_field(elem(ie), 'omega', wg1, hvcoord, nt, -1)
       call gfr_g2f_scalar(ie, elem(ie)%metdet, wg1, wf1)
#ifdef MODEL_THETA_L
       omega_p(:nf2,:,ie) = wf1(:nf2,:)
#else
       ! for preqx, omega_p = omega/p
       omega_p(:nf2,:,ie) = wf1(:nf2,:)/p_fv(:nf2,:)
#endif

       do qi = 1,qsize
          call gfr_g2f_mixing_ratio(gfr, ie, elem(ie)%metdet, dp, dp_fv, &
               elem(ie)%state%Q(:,:,:,qi), wf1)
          q(:nf2,:,qi,ie) = wf1(:nf2,:)
          if (gfr%check > 1) then
             nerr = check_g2f_mixing_ratio(gfr, hybrid, ie, qi, elem, dp, dp_fv, &
                  elem(ie)%state%Q(:,:,:,qi), wf1)
             if (nerr > 0) gfr%check_ok = .false.
          end if
       end do
    end do

    if (gfr%check > 0) then
       nerr = check_global_properties(gfr, hybrid, hvcoord, elem, nt, nets, nete, .true., q)
       if (nerr > 0) gfr%check_ok = .false.
    end if
  end subroutine gfr_dyn_to_fv_phys_hybrid

  subroutine gfr_fv_phys_to_dyn_hybrid(hybrid, nt, hvcoord, elem, nets, nete, T, uv, q)
    ! Remap T, uv, and q tendencies and from FV to GLL grids.
    !   On input, T and uv are tendencies and q is state.
    !   On output, FM and FT are tendencies and FQ is state.

#ifdef __INTEL_COMPILER
# if __INTEL_COMPILER >= 1700 && __INTEL_COMPILER < 1800
    ! On Anvil with Intel 17, and reproduced on one other platform
    ! using specfically
    !   icpc version 17.0.0 (gcc version 4.7.4 compatibility)
    ! -O3 causes buggy code to be emitted for the limiter block near
    ! the end of this routine. Work around this by asking the compiler
    ! to compile this routine no higher than -O2.
    !DIR$ OPTIMIZE:2
# endif
#endif

    use element_ops, only: get_field
    use dimensions_mod, only: nlev
    use hybvcoord_mod, only: hvcoord_t
    use physical_constants, only: p0, kappa

    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: nt
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete
    real(kind=real_kind), intent(in) :: T(:,:,:), uv(:,:,:,:), q(:,:,:,:)

    real(kind=real_kind), dimension(np,np,nlev) :: dp, wg1, p
    real(kind=real_kind), dimension(np*np,nlev) :: wf1, wf2, dp_fv, p_fv
    real(kind=real_kind) :: qmin, qmax
    integer :: ie, nf, nf2, k, qsize, qi, nerr

    nf = gfr%nphys
    nf2 = nf*nf
    qsize = size(q,3)

    do ie = nets,nete
       dp = elem(ie)%state%dp3d(:,:,:,nt)
       call gfr_g2f_scalar(ie, elem(ie)%metdet, elem(ie)%state%ps_v(:,:,nt:nt), wf1(:,:1))
       call calc_dp_fv(nf, hvcoord, wf1(:,1), dp_fv)

       call gfr_f2g_vector(gfr, ie, elem, &
            uv(:nf2,1,:,ie), uv(:nf2,2,:,ie), &
            elem(ie)%derived%FM(:,:,1,:), elem(ie)%derived%FM(:,:,2,:))

       call get_field(elem(ie), 'p', p, hvcoord, nt, -1)
       call gfr_g2f_scalar(ie, elem(ie)%metdet, p, p_fv)
#ifndef HOMMEXX_BFB_TESTING
       wf1(:nf2,:) = T(:nf2,:,ie)*(p0/p_fv(:nf2,:))**kappa
       call gfr_f2g_scalar_dp(gfr, ie, elem(ie)%metdet, dp_fv, dp, wf1, &
            elem(ie)%derived%FT)
       elem(ie)%derived%FT = elem(ie)%derived%FT*(p/p0)**kappa
#else
       wf1(:nf2,:) = T(:nf2,:,ie)*bfb_pow(p0/p_fv(:nf2,:), kappa)
       call gfr_f2g_scalar_dp(gfr, ie, elem(ie)%metdet, dp_fv, dp, wf1, &
            elem(ie)%derived%FT)
       elem(ie)%derived%FT = elem(ie)%derived%FT*bfb_pow(p/p0, kappa)
#endif

       do qi = 1,qsize
          ! FV Q_ten
          !   GLL Q0 -> FV Q0
          call gfr_g2f_mixing_ratio(gfr, ie, elem(ie)%metdet, dp, dp_fv, &
               elem(ie)%state%Q(:,:,:,qi), wf1)
          !   FV Q_ten = FV Q1 - FV Q0
          wf1(:nf2,:) = q(:nf2,:,qi,ie) - wf1(:nf2,:)
          if (nf > 1 .or. .not. gfr%boost_pg1) then
             ! GLL Q_ten
             call gfr_f2g_scalar_dp(gfr, ie, elem(ie)%metdet, dp_fv, dp, wf1, wg1)
             ! GLL Q1
             elem(ie)%derived%FQ(:,:,:,qi) = elem(ie)%state%Q(:,:,:,qi) + wg1
          else
             ! GLL Q_ten
             do k = 1,nlev
                elem(ie)%derived%FQ(:,:,k,qi) = wf1(1,k)
             end do
          end if
          ! Get limiter bounds.
          do k = 1,nlev
             gfr%qmin(k,qi,ie) = minval(q(:nf2,k,qi,ie))
             gfr%qmax(k,qi,ie) = maxval(q(:nf2,k,qi,ie))
          end do
       end do
    end do

    ! Halo exchange limiter bounds.
    call gfr_f2g_mixing_ratios_he(hybrid, nets, nete, gfr%qmin(:,:,nets:nete), &
         gfr%qmax(:,:,nets:nete))

    if (nf == 1 .and. gfr%boost_pg1) return

    do ie = nets,nete
       dp = elem(ie)%state%dp3d(:,:,:,nt)
       do qi = 1,qsize
          ! Limit GLL Q1.
          if (gfr%check > 1) wg1 = elem(ie)%derived%FQ(:,:,:,qi)
          do k = 1,nlev
             ! Augment bounds with GLL Q0 bounds. This assures that if
             ! the tendency is 0, GLL Q1 = GLL Q0.
             gfr%qmin(k,qi,ie) = min(minval(elem(ie)%state%Q(:,:,k,qi)), gfr%qmin(k,qi,ie))
             gfr%qmax(k,qi,ie) = max(maxval(elem(ie)%state%Q(:,:,k,qi)), gfr%qmax(k,qi,ie))
             ! Final GLL Q1, except for DSS, which is not done in this routine.
             call limiter_clip_and_sum(elem(ie)%spheremp, gfr%qmin(k,qi,ie), &
                  gfr%qmax(k,qi,ie), dp(:,:,k), elem(ie)%derived%FQ(:,:,k,qi))
          end do
          if (gfr%check > 1) then
             nerr = check_f2g_mixing_ratio(gfr, hybrid, ie, qi, elem, gfr%qmin(:,qi,ie), &
                  gfr%qmax(:,qi,ie), dp, wg1, elem(ie)%derived%FQ(:,:,:,qi))
             if (nerr > 0) gfr%check_ok = .false.
          end if
       end do
    end do

    if (gfr%check > 0) then
       nerr = check_global_properties(gfr, hybrid, hvcoord, elem, nt, nets, nete, .false., q)
       if (nerr > 0) gfr%check_ok = .false.
    end if
  end subroutine gfr_fv_phys_to_dyn_hybrid

  subroutine gfr_dyn_to_fv_phys_topo_hybrid(hybrid, elem, nets, nete, phis)
    ! If needed, remap topography data defined on the GLL grid to the
    ! FV grid. The intended EAM configuration is to use topography
    ! data on the FV grid, so this routine is unlikely to be used.

    type (hybrid_t), intent(in) :: hybrid
    type (element_t), intent(in) :: elem(:)
    integer, intent(in) :: nets, nete
    real(kind=real_kind), intent(out) :: phis(:,:)

    integer :: ie

    do ie = nets,nete
       call gfr_g2f_scalar_and_limit(gfr, ie, elem(ie)%metdet, elem(ie)%state%phis, phis(:,ie))
    end do
  end subroutine gfr_dyn_to_fv_phys_topo_hybrid

  subroutine gfr_dyn_to_fv_phys_topo_data(par, elem, nets, nete, g, gsz, p, psz, square, augment)
    ! Remap SGH, SGH30, phis, landm_coslat, landfrac from GLL to FV grids. For
    ! SGH* fields, square should be true. Then variance is remapped and so
    ! conserved. For SGH, but not SGH30, set augment to true to add the variance
    ! increases due to the additional smoothness of remapping from GLL to
    ! (first-order) FV bases. SGH30 does not require additional variance to be
    ! added since it's the variance due to truncation of wave numbers at a
    ! separation length scale; this is independent of grid.

    use parallel_mod, only: parallel_t, abortmp
    use hybrid_mod, only: hybrid_create

    type (parallel_t), intent(in) :: par
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete, gsz, psz
    real(kind=real_kind), intent(in) :: g(gsz)
    real(kind=real_kind), intent(out) :: p(psz)
    logical, intent(in), optional :: square, augment

    type (hybrid_t) :: hybrid
    integer :: ie, nf2
    logical :: augment_in

    augment_in = .false.
    if (present(augment)) augment_in = augment

    nf2 = gfr%nphys*gfr%nphys

    do ie = nets,nete
       call gfr_dyn_to_fv_phys_topo_data_elem(ie, elem, square, augment_in, &
            g(npsq*(ie-nets)+1 : npsq*(ie-nets+1)), &
            p(nf2*(ie-nets)+1 : nf2*(ie-nets+1)))
    end do
  end subroutine gfr_dyn_to_fv_phys_topo_data

  subroutine gfr_dyn_to_fv_phys_topo_data_elem(ie, elem, square, augment_variance, g, p)
    ! Element-level impl of gfr_dyn_to_fv_phys_topo_data.

    use physical_constants, only: grav => g

    integer, intent(in) :: ie
    type (element_t), intent(in) :: elem(:)
    logical, intent(in) :: square, augment_variance
    real(kind=real_kind), intent(in) :: g(:)
    real(kind=real_kind), intent(out) :: p(:)

    real(kind=real_kind) :: wg(np,np,1), wf(np*np,2), ones(np,np), qmin, qmax, phispg(npsq)
    integer :: nf, nf2, i, j, k

    ones = one
    nf = gfr%nphys
    nf2 = nf*nf

    if (augment_variance) then
       ! Estimate additional variance due to remapping from GLL to FV
       ! bases.
       call gfr_g2f_scalar_and_limit(gfr, ie, elem(ie)%metdet, elem(ie)%state%phis, phispg)
       do k = 1,nf2
          ! Integrate (phis_gll - phis_fv)^2 over FV subcell (i,j). Do this
          ! using gfr_g2f_scalar; thus, only one entry out of nf^2 is used.
          wg(:,:,1) = ((elem(ie)%state%phis - phispg(k))/grav)**2
          call gfr_g2f_scalar(ie, elem(ie)%metdet, wg(:,:,:1), wf(:,1:1))
          ! Use just entry (i,j).
          wf(k,2) = max(zero, wf(k,1))
       end do

       ! Original SGH. augment_variance implies we need to square and sqrt
       ! quantities.
       wg(:,:,1) = reshape(g(:npsq)**2, (/np,np/))
       call gfr_g2f_scalar_and_limit(gfr, ie, elem(ie)%metdet, wg(:,:,1), p(:nf2))

       ! Combine the two sources of variance.
       p(:nf2) = sqrt(p(:nf2) + wf(:nf2,2))
    else
       wg(:,:,1) = reshape(g(:npsq), (/np,np/))
       if (square) wg(:,:,1) = wg(:,:,1)**2
       call gfr_g2f_scalar_and_limit(gfr, ie, elem(ie)%metdet, wg(:,:,1), p(:nf2))
       if (square) p(:nf2) = sqrt(p(:nf2))
    end if
  end subroutine gfr_dyn_to_fv_phys_topo_data_elem
  
  subroutine gfr_fv_phys_to_dyn_topo_hybrid(hybrid, elem, nets, nete, phis)
    ! Remap FV topography data to the GLL grid. Prevent new
    ! extrema. Conserve the integral of height.

    use kinds, only: iulog
    use edge_mod, only: edgeVpack_nlyr, edgeVunpack_nlyr, edge_g
    use bndry_mod, only: bndry_exchangeV

    type (hybrid_t), intent(in) :: hybrid
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete
    real(kind=real_kind), intent(in) :: phis(:,:)

    real(kind=real_kind) :: wg(np,np,2), ones(np,np,1)
    integer :: ie, nf, nf2, nerr

    ones = one
    nf = gfr%nphys
    nf2 = nf*nf
    ! For now, map GLL topo back to FV, as ne30pg2 runs otherwise show grid
    ! imprint in the PHIS diagnostic. We may revise this choice later.
    !gfr%have_fv_topo_file_phis = .true.

    do ie = nets,nete
       gfr%phis(:nf2,ie) = phis(:nf2,ie)
       gfr%qmin(:,:,ie) = minval(phis(:nf2,ie))
       gfr%qmax(:,:,ie) = maxval(phis(:nf2,ie))
       if (nf > 1) then
          call gfr_f2g_scalar(ie, elem(ie)%metdet, phis(:nf2,ie:ie), wg(:,:,2:2))
          elem(ie)%state%phis = wg(:,:,2)
       else
          elem(ie)%state%phis = phis(1,ie)
       end if
    end do

    call gfr_f2g_mixing_ratios_he(hybrid, nets, nete, gfr%qmin(:,:,nets:nete), &
         gfr%qmax(:,:,nets:nete))

    if (nf > 1 .or. .not. gfr%boost_pg1) then
       do ie = nets,nete
          if (gfr%check > 1) wg(:,:,1) = elem(ie)%state%phis
          call limiter_clip_and_sum(elem(ie)%spheremp, gfr%qmin(1,1,ie), &
               gfr%qmax(1,1,ie), ones(:,:,1), elem(ie)%state%phis)
          if (gfr%check > 1) then
             if (gfr%qmin(1,1,ie) < zero) then
                write(iulog,*) 'gfr> topo min:', hybrid%par%rank, hybrid%ithr, ie, &
                     gfr%qmin(1,1,ie), 'ERROR'
                gfr%check_ok = .false.
             end if
             wg(:,:,2) = elem(ie)%state%phis
             nerr = check_f2g_mixing_ratio(gfr, hybrid, ie, 1, elem, gfr%qmin(:1,1,ie), &
                  gfr%qmax(:1,1,ie), ones, wg(:,:,:1), wg(:,:,2:))
             if (nerr > 0) gfr%check_ok = .false.
          end if
       end do
    end if

    do ie = nets,nete
       elem(ie)%state%phis = elem(ie)%state%phis*elem(ie)%spheremp*elem(ie)%rspheremp
    end do
  end subroutine gfr_fv_phys_to_dyn_topo_hybrid

  subroutine gfr_dyn_to_fv_phys_dom_mt(par, dom_mt, nt, hvcoord, elem, &
       ps, phis, T, uv, omega_p, q)
    ! Wrapper to the hybrid-threading main routine.
    
    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use hybvcoord_mod, only: hvcoord_t
    use thread_mod, only: hthreads

    type (parallel_t), intent(in) :: par
    type (domain1d_t), intent(in) :: dom_mt(:)
    integer, intent(in) :: nt
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(out) :: ps(:,:), phis(:,:), T(:,:,:), &
         uv(:,:,:,:), omega_p(:,:,:), q(:,:,:,:)

    type (hybrid_t) :: hybrid
    integer :: nets, nete

    if (.not. par%dynproc) return
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

  subroutine gfr_fv_phys_to_dyn_dom_mt(par, dom_mt, nt, hvcoord, elem, T, uv, q)
    ! Wrapper to the hybrid-threading main routine.

    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use hybvcoord_mod, only: hvcoord_t
    use thread_mod, only: hthreads

    type (parallel_t), intent(in) :: par
    type (domain1d_t), intent(in) :: dom_mt(:)
    integer, intent(in) :: nt
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(inout) :: elem(:)
    real(kind=real_kind), intent(inout) :: T(:,:,:), uv(:,:,:,:), q(:,:,:,:)

    type (hybrid_t) :: hybrid
    integer :: nets, nete

    if (.not. par%dynproc) return
#ifdef HORIZ_OPENMP
    !$omp parallel num_threads(hthreads), default(shared), private(nets,nete,hybrid)
#endif
    call gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)
    call gfr_fv_phys_to_dyn_hybrid(hybrid, nt, hvcoord, elem, nets, nete, T, uv, q)
#ifdef HORIZ_OPENMP
    !$omp end parallel
#endif
  end subroutine gfr_fv_phys_to_dyn_dom_mt

  subroutine gfr_dyn_to_fv_phys_topo_dom_mt(par, dom_mt, elem, phis)
    ! Wrapper to the hybrid-threading main routine.

    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use thread_mod, only: hthreads

    type (parallel_t), intent(in) :: par
    type (domain1d_t), intent(in) :: dom_mt(:)
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(out) :: phis(:,:)

    type (hybrid_t) :: hybrid
    integer :: nets, nete

    if (.not. par%dynproc) return
#ifdef HORIZ_OPENMP
    !$omp parallel num_threads(hthreads), default(shared), private(nets,nete,hybrid)
#endif
    call gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)
    call gfr_dyn_to_fv_phys_topo_hybrid(hybrid, elem, nets, nete, phis)
#ifdef HORIZ_OPENMP
    !$omp end parallel
#endif
  end subroutine gfr_dyn_to_fv_phys_topo_dom_mt

  subroutine gfr_dyn_to_fv_phys_topo_mpi_only(par, elem, phis)
    ! Wrapper to the hybrid-threading main routine.

    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use hybrid_mod, only: hybrid_create

    type (parallel_t), intent(in) :: par
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(out) :: phis(:,:)

    type (hybrid_t) :: hybrid
    integer :: nets, nete

    if (.not. par%dynproc) return
    hybrid = hybrid_create(par, 0, 1)
    call gfr_dyn_to_fv_phys_topo_hybrid(hybrid, elem, 1, nelemd, phis)
  end subroutine gfr_dyn_to_fv_phys_topo_mpi_only

  subroutine gfr_fv_phys_to_dyn_topo_dom_mt(par, dom_mt, elem, phis)
    ! Wrapper to the hybrid-threading main routine.

    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use thread_mod, only: hthreads

    type (parallel_t), intent(in) :: par
    type (domain1d_t), intent(in) :: dom_mt(:)
    type (element_t), intent(inout) :: elem(:)
    real(kind=real_kind), intent(in) :: phis(:,:)

    type (hybrid_t) :: hybrid
    integer :: nets, nete

    if (.not. par%dynproc) return
#ifdef HORIZ_OPENMP
    !$omp parallel num_threads(hthreads), default(shared), private(nets,nete,hybrid)
#endif
    call gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)
    call gfr_fv_phys_to_dyn_topo_hybrid(hybrid, elem, nets, nete, phis)
#ifdef HORIZ_OPENMP
    !$omp end parallel
#endif
  end subroutine gfr_fv_phys_to_dyn_topo_dom_mt

  subroutine gfr_fv_phys_to_dyn_topo_mpi_only(par, elem, phis)
    ! Wrapper to the hybrid-threading main routine.

    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use hybrid_mod, only: hybrid_create

    type (parallel_t), intent(in) :: par
    type (element_t), intent(inout) :: elem(:)
    real(kind=real_kind), intent(in) :: phis(:,:)

    type (hybrid_t) :: hybrid
    integer :: nets, nete

    if (.not. par%dynproc) return
    hybrid = hybrid_create(par, 0, 1)
    call gfr_fv_phys_to_dyn_topo_hybrid(hybrid, elem, 1, nelemd, phis)
  end subroutine gfr_fv_phys_to_dyn_topo_mpi_only

  ! ----------------------------------------------------------------------
  ! Internal initialization routines.

  subroutine gfr_init_w_gg(np, w_gg)
    ! Init GLL w(i)*w(j) values on the reference element.

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
    ! Init FV w(i)*w(j) values on the reference element.
    
    integer, intent(in) :: nphys
    real(kind=real_kind), intent(out) :: w_ff(:)

    w_ff(:nphys*nphys) = four/real(nphys*nphys, real_kind)
  end subroutine gfr_init_w_ff

  subroutine gll_cleanup(gll)
    ! Deallocate the GLL object's data.

    use quadrature_mod, only : quadrature_t

    type (quadrature_t), intent(inout) :: gll

    deallocate(gll%points, gll%weights)
  end subroutine gll_cleanup

  subroutine eval_lagrange_bases(gll, np, x, y)
    ! Evaluate the GLL basis functions at x in [-1,1], writing the
    ! values to y(1:np). This implements the Lagrange interpolant.

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
    ! Compute the mixed mass matrix with range the FV subcells and
    ! domain the GLL nodes.

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
             ref = xs + half*(xe - xs)*(one + gll%points(qj))
             call eval_lagrange_bases(gll, np, ref, bj)
             do qi = 1,np
                ref = ys + half*(ye - ys)*(one + gll%points(qi))
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
    ! Compute the matrix that interpolates from the npi-GLL nodes to
    ! the np-GLL nodes.

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
    ! Apply gfr_init_f2g_remapd_op to the Id matrix to get the remap operator's
    ! matrix representation.

    !assume nphys <= np

    type (GllFvRemap_t), intent(inout) :: gfr
    real(kind=real_kind), intent(in) :: R(:,:), tau(:)

    integer :: nf, fi, fj, gi, gj
    real(kind=real_kind) :: f(np,np), g(np,np)

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
    ! This operator implements the linear operator that solves the
    ! problem described in gfr_init_R.

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
       call dtrsm('L', 'U', 'T', 'N', nf2, 1, one, R, size(R,1), wrk, nf2)
       call dormqr('L', 'N', nf2, 1, nf2, R, size(R,1), tau, wrk, nf2, wr, np2, info)
       g(:npi,:npi) =  wrk
    else
       call dtrtrs('U', 'T', 'N', nf2, 1, R, size(R,1), wrk, nf2, info)
       call dtrtrs('U', 'N', 'N', nf2, 1, R, size(R,1), wrk, nf2, info)
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
       call apply_interp(gfr%interp, np, npi, x, g)
    elseif (nf < npi) then
       ! Finish the projection.
       do gj = 1,np
          do gi = 1,np
             g(gi,gj) = g(gi,gj)/gfr%w_gg(gi,gj)
          end do
       end do
    end if
  end subroutine gfr_f2g_remapd_op

  subroutine gfr_init_geometry(elem, gfr)
    use kinds, only: iulog
    use control_mod, only: cubed_sphere_map
    use coordinate_systems_mod, only: cartesian3D_t, spherical_polar_t, &
         sphere_tri_area, change_coordinates
    use cube_mod, only: ref2sphere

    type (element_t), intent(in) :: elem(:)
    type (GllFvRemap_t), intent(inout) :: gfr

    type (spherical_polar_t) :: p_sphere
    type (cartesian3D_t) :: fv_corners_xyz(2,2), ctr
    real(kind=real_kind) :: ones(np*np), ones2(np,np), ae(2), be(2), &
         spherical_area, tmp, ac, bc
    integer :: nf, nf2, ie, i, j, k, ai, bi, idx

    nf = gfr%nphys
    nf2 = nf*nf
    do ie = 1,nelemd
       do j = 1,nf
          call gfr_f_ref_edges(nf, j, be)
          call gfr_f_ref_center(nf, j, bc)
          do i = 1,nf
             call gfr_f_ref_edges(nf, i, ae)
             call gfr_f_ref_center(nf, i, ac)
             k = i+(j-1)*nf
             ctr%x = zero; ctr%y = zero; ctr%z = zero
             do bi = 1,2
                do ai = 1,2
                   if ( (i == 1  .and. ai == 1 .and. j == 1  .and. bi == 1) .or. &
                        (i == 1  .and. ai == 1 .and. j == nf .and. bi == 2) .or. &
                        (i == nf .and. ai == 2 .and. j == 1  .and. bi == 1) .or. &
                        (i == nf .and. ai == 2 .and. j == nf .and. bi == 2)) then
                      ! Use the element corner if we are at it.
                      idx = 2*(bi-1)
                      if (bi == 1) then
                         idx = idx + ai
                      else
                         idx = idx + 3 - ai
                      end if
                      fv_corners_xyz(ai,bi) = elem(ie)%corners3D(idx)
                   else
                      ! p_sphere is unused. fv_corners_xyz(ai,bi) contains
                      ! the cartesian point before it's converted to lat-lon.
                      p_sphere = ref2sphere(ae(ai), be(bi), elem(ie)%corners3D, &
                           cubed_sphere_map, elem(ie)%corners, elem(ie)%facenum, &
                           fv_corners_xyz(ai,bi))
                      if (cubed_sphere_map == 0) then
                         ! In this case, fv_corners_xyz above is not set.
                         fv_corners_xyz(ai,bi) = change_coordinates(p_sphere)
                      end if
                   end if
                   ctr%x = ctr%x + fv_corners_xyz(ai,bi)%x
                   ctr%y = ctr%y + fv_corners_xyz(ai,bi)%y
                   ctr%z = ctr%z + fv_corners_xyz(ai,bi)%z
                end do
             end do

             if (cubed_sphere_map == 2) then
                call sphere_tri_area(fv_corners_xyz(1,1), fv_corners_xyz(2,1), &
                     fv_corners_xyz(2,2), spherical_area)
                call sphere_tri_area(fv_corners_xyz(1,1), fv_corners_xyz(2,2), &
                     fv_corners_xyz(1,2), tmp)
                spherical_area = spherical_area + tmp
                gfr%fv_metdet(k,ie) = spherical_area/gfr%w_ff(k)

                ! Center is average of 4 corner points projected to sphere.
                ctr%x = ctr%x/four; ctr%y = ctr%y/four; ctr%z = ctr%z/four
                gfr%center_f(i,j,ie) = ctr
             end if

             do bi = 1,2
                do ai = 1,2
                   ! CCW with ai the fast direction.
                   idx = 2*(bi-1)
                   if (bi == 1) then
                      idx = idx + ai
                   else
                      idx = idx + 3 - ai
                   end if
                   gfr%corners_f(idx,i,j,ie) = fv_corners_xyz(ai,bi)
                end do
             end do
          end do
       end do
    end do

    if (cubed_sphere_map == 0) then
       ! For cubed_sphere_map == 0, we set the center so that it maps to the ref
       ! element center and set fv_metdet so that it corresponds to the integral
       ! of metdet over the FV subcell. TempestRemap establishes the
       ! cubed_sphere_map == 2 convention, but for cubed_sphere_map == 0 there
       ! is no external convention.
       ones = one
       ones2 = one
       do ie = 1,nelemd
          call gfr_g2f_remapd(gfr, elem(ie)%metdet, ones, ones2, gfr%fv_metdet(:nf2,ie))
          do j = 1,nf
             call gfr_f_ref_center(nf, j, bc)
             do i = 1,nf
                call gfr_f_ref_center(nf, i, ac)
                p_sphere = ref2sphere(ac, bc, elem(ie)%corners3D, cubed_sphere_map, &
                     elem(ie)%corners, elem(ie)%facenum)
                gfr%center_f(i,j,ie) = change_coordinates(p_sphere)
             end do
          end do
       end do
    end if

    ! Make the spherical area of the element according to FV and GLL agree to
    ! machine precision.
    if (gfr%check > 0) allocate(gfr%check_areas(1,nelemd))
    do ie = 1,nelemd
       if (gfr%check > 0) gfr%check_areas(1,ie) = sum(gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie))
       gfr%fv_metdet(:nf2,ie) = gfr%fv_metdet(:nf2,ie)* &
            (sum(elem(ie)%spheremp)/sum(gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie)))
    end do
  end subroutine gfr_init_geometry

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

  subroutine gfr_init_Dmap(elem, gfr)
    use control_mod, only: cubed_sphere_map
    use cube_mod, only: Dmap, ref2sphere
    use coordinate_systems_mod, only: cartesian3D_t, change_coordinates

    type (element_t), intent(in) :: elem(:)
    type (GllFvRemap_t), intent(inout) :: gfr

    type (cartesian3D_t) :: sphere
    real(kind=real_kind) :: wrk(2,2), det, a, b
    integer :: ie, nf, nf2, i, j, k

    nf = gfr%nphys
    nf2 = nf*nf
    
    ! Jacobian matrices to map a vector between reference element and sphere.
    do ie = 1,nelemd
       do j = 1,nf
          do i = 1,nf
             if (cubed_sphere_map == 2) then
                call gfr_f_get_cartesian3d(ie, i, j, sphere)
                call sphere2ref(elem(ie)%corners3D, sphere, a, b)
             else
                call gfr_f_ref_center(nf, i, a)
                call gfr_f_ref_center(nf, j, b)
             end if

             call Dmap(wrk, a, b, elem(ie)%corners3D, cubed_sphere_map, elem(ie)%cartp, &
                  elem(ie)%facenum)

             det = wrk(1,1)*wrk(2,2) - wrk(1,2)*wrk(2,1)

             ! Make det(D) = fv_metdet. The two should be equal, and fv_metdet
             ! must be consistent with spherep. Thus, D must be adjusted by a
             ! scalar.
             k = i + (j-1)*nf
             wrk = wrk*sqrt(gfr%fv_metdet(k,ie)/abs(det))
             det = gfr%fv_metdet(k,ie)

             gfr%D_f(k,:,:,ie) = wrk

             gfr%Dinv_f(k,1,1,ie) =  wrk(2,2)/det
             gfr%Dinv_f(k,1,2,ie) = -wrk(1,2)/det
             gfr%Dinv_f(k,2,1,ie) = -wrk(2,1)/det
             gfr%Dinv_f(k,2,2,ie) =  wrk(1,1)/det
          end do
       end do
    end do
  end subroutine gfr_init_Dmap

  ! ----------------------------------------------------------------------
  ! Time stepping routines called by the GLL <-> FV remap API routines.

  ! GLL -> FV (g2f)

  subroutine gfr_g2f_remapd(gfr, gll_metdet, fv_metdet, g, f)
    ! Core remap operator. Conservative remap on the reference
    ! element.

    type (GllFvRemap_t), intent(in) :: gfr
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), fv_metdet(:), g(:,:)
    real(kind=real_kind), intent(out) :: f(:)

    integer :: nf, nf2, gi, gj, k
    real(kind=real_kind) :: gw(np,np)

    nf = gfr%nphys
    nf2 = nf*nf
    gw = g*gll_metdet
    do k = 1,nf2
       f(k) = sum(gfr%g2f_remapd(:,:,k)*gw)/(gfr%w_ff(k)*fv_metdet(k))
    end do
  end subroutine gfr_g2f_remapd

  subroutine gfr_g2f_scalar(ie, gll_metdet, g, f) ! no gfr b/c public
    ! Wrapper to remapd, where g and f are densities.

    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), g(:,:,:)
    real(kind=real_kind), intent(out) :: f(:,:)

    integer :: nlev, k

    nlev = size(g,3)
    do k = 1, nlev
       call gfr_g2f_remapd(gfr, gll_metdet, gfr%fv_metdet(:,ie), g(:,:,k), f(:,k))
    end do
  end subroutine gfr_g2f_scalar

  subroutine gfr_g2f_scalar_dp(gfr, ie, gll_metdet, dp_g, dp_f, g, f)
    ! Wrapper to remapd, where g and f are mixing ratios.

    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), dp_g(:,:,:), dp_f(:,:), g(:,:,:)
    real(kind=real_kind), intent(out) :: f(:,:)

    integer :: nf2

    nf2 = gfr%nphys*gfr%nphys
    call gfr_g2f_scalar(ie, gll_metdet, dp_g*g, f)
    f(:nf2,:) = f(:nf2,:)/dp_f(:nf2,:)
  end subroutine gfr_g2f_scalar_dp

  subroutine gfr_g2f_vector(ie, elem, u_g, v_g, u_f, v_f) ! no gfr b/c public
    ! Remap a vector on the sphere by doing the actual remap on the
    ! reference element, thus avoiding steep gradients at the poles.

    integer, intent(in) :: ie
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(in) :: u_g(:,:,:), v_g(:,:,:)
    real(kind=real_kind), intent(out) :: u_f(:,:), v_f(:,:)

    real(kind=real_kind) :: wg(np,np,2), wf(np*np,2), ones(np*np), ones2(np,np)
    integer :: k, d, nf, nf2, nlev

    nf = gfr%nphys
    nf2 = nf*nf
    ones = one
    ones2 = one

    nlev = size(u_g,3)
    do k = 1, nlev
       ! sphere -> GLL ref
       do d = 1,2
          wg(:,:,d) = elem(ie)%Dinv(:,:,d,1)*u_g(:,:,k) + elem(ie)%Dinv(:,:,d,2)*v_g(:,:,k)
       end do
       do d = 1,2
          call gfr_g2f_remapd(gfr, ones2, ones, wg(:,:,d), wf(:,d))
       end do
       ! FV ref -> sphere
       u_f(:nf2,k) = gfr%D_f(:nf2,1,1,ie)*wf(:nf2,1) + gfr%D_f(:nf2,1,2,ie)*wf(:nf2,2)
       v_f(:nf2,k) = gfr%D_f(:nf2,2,1,ie)*wf(:nf2,1) + gfr%D_f(:nf2,2,2,ie)*wf(:nf2,2)
    end do
  end subroutine gfr_g2f_vector

  subroutine gfr_g2f_mixing_ratio(gfr, ie, gll_metdet, dp_g, dp_f, q_g, q_f)
    ! Remap a mixing ratio conservatively and preventing new extrema.

    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), dp_g(:,:,:), dp_f(:,:), q_g(:,:,:)
    real(kind=real_kind), intent(out) :: q_f(:,:)

    real(kind=real_kind) :: qmin, qmax, wg(np,np), wf1(np*np), wf2(np*np)
    integer :: q, k, nf, nf2, nlev

    nf = gfr%nphys
    nf2 = nf*nf
    nlev = size(q_g,3)
    do k = 1,nlev
       wg = q_g(:,:,k)
       qmin = minval(wg)
       qmax = maxval(wg)
       wg = dp_g(:,:,k)*wg
       call gfr_g2f_remapd(gfr, gll_metdet, gfr%fv_metdet(:,ie), wg, wf1)
       wf1(:nf2) = wf1(:nf2)/dp_f(:nf2,k)
       wf2(:nf2) = gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie)
       call limiter1_clip_and_sum(nf2, wf2, qmin, qmax, dp_f(:,k), wf1)
       q_f(:nf2,k) = wf1(:nf2)
    end do
  end subroutine gfr_g2f_mixing_ratio

  subroutine gfr_g2f_scalar_and_limit(gfr, ie, gll_metdet, g, f)
    ! After remap, limit using extremal values from g.

    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), g(:,:)
    real(kind=real_kind), intent(out) :: f(:)

    real(kind=real_kind) :: wg(np,np,1), wf(np*np,1), ones(np*np), qmin, qmax
    integer :: nf, nf2

    ones = one
    nf = gfr%nphys
    nf2 = nf*nf

    qmin = minval(g(:np,:np))
    qmax = maxval(g(:np,:np))
    wg(:np,:np,1) = g
    call gfr_g2f_scalar(ie, gll_metdet, wg(:,:,:1), wf(:,:1))
    f(:nf2) = wf(:nf2,1)
    wf(:nf2,1) = gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie)
    call limiter1_clip_and_sum(nf2, wf(:,1), qmin, qmax, ones, f)
  end subroutine gfr_g2f_scalar_and_limit

  ! FV -> GLL (f2g)

  subroutine gfr_f2g_scalar(ie, gll_metdet, f, g) ! no gfr b/c public for testing
    ! Wrapper to remapd, where g and f are densities.

    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), f(:,:)
    real(kind=real_kind), intent(out) :: g(:,:,:)

    integer :: k

    do k = 1, size(g,3)
       call gfr_f2g_remapd(gfr, gll_metdet, gfr%fv_metdet(:,ie), f(:,k), g(:,:,k))
    end do
  end subroutine gfr_f2g_scalar

  subroutine gfr_f2g_scalar_dp(gfr, ie, gll_metdet, dp_f, dp_g, f, g)
    ! Wrapper to remapd, where g and f are mixing ratios.

    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), dp_f(:,:), dp_g(:,:,:), f(:,:)
    real(kind=real_kind), intent(out) :: g(:,:,:)

    integer :: nf2

    nf2 = gfr%nphys*gfr%nphys
    call gfr_f2g_scalar(ie, gll_metdet, dp_f(:nf2,:)*f(:nf2,:), g)
    g = g/dp_g
  end subroutine gfr_f2g_scalar_dp

  subroutine gfr_f2g_vector(gfr, ie, elem, u_f, v_f, u_g, v_g)
    ! Remap a vector on the sphere by doing the actual remap on the
    ! reference element, thus avoiding steep gradients at the poles.

    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(in) :: u_f(:,:), v_f(:,:)
    real(kind=real_kind), intent(out) :: u_g(:,:,:), v_g(:,:,:)

    real(kind=real_kind) :: wg(np,np,2), wf(np*np,2), ones(np*np), ones2(np,np)
    integer :: k, d, nf, nf2, nlev

    nf = gfr%nphys
    nf2 = nf*nf
    ones = one
    ones2 = one

    nlev = size(u_g,3)
    do k = 1, nlev
       ! sphere -> FV ref
       do d = 1,2
          wf(:nf2,d) = &
               gfr%Dinv_f(:nf2,d,1,ie)*u_f(:nf2,k) + &
               gfr%Dinv_f(:nf2,d,2,ie)*v_f(:nf2,k)
       end do
       do d = 1,2
          call gfr_f2g_remapd(gfr, ones2, ones, wf(:,d), wg(:,:,d))
       end do
       ! GLL ref -> sphere
       u_g(:,:,k) = elem(ie)%D(:,:,1,1)*wg(:,:,1) + elem(ie)%D(:,:,1,2)*wg(:,:,2)
       v_g(:,:,k) = elem(ie)%D(:,:,2,1)*wg(:,:,1) + elem(ie)%D(:,:,2,2)*wg(:,:,2)
    end do
  end subroutine gfr_f2g_vector

  subroutine gfr_f2g_mixing_ratios_he(hybrid, nets, nete, qmin, qmax)
    ! Exchange qmin/qmax among element neighbors.

    use viscosity_mod, only: neighbor_minmax
    use prim_advection_base, only: edgeAdvQminmax

    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: nets, nete
    real(kind=real_kind), intent(inout) :: qmin(:,:,:), qmax(:,:,:)

    if (hybrid%par%dynproc) call neighbor_minmax(hybrid, edgeAdvQminmax, nets, nete, qmin, qmax)
  end subroutine gfr_f2g_mixing_ratios_he

  subroutine gfr_f2g_dss(hybrid, elem, nets, nete)
    ! DSS FQ, FM, FT.

    use dimensions_mod, only: nlev, qsize
    use edge_mod, only: edgevpack_nlyr, edgevunpack_nlyr, edge_g
    use bndry_mod, only: bndry_exchangev

    type (hybrid_t), intent(in) :: hybrid
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete

    real(kind=real_kind) :: tmp(np,np,nlev)
    integer :: ie, q, k, npack

    npack = (qsize + 3)*nlev
    do ie = nets, nete
       do q = 1,qsize
          do k = 1,nlev
             elem(ie)%derived%FQ(:,:,k,q) = elem(ie)%derived%FQ(:,:,k,q)*elem(ie)%spheremp(:,:)
          end do
       end do
       call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FQ, qsize*nlev, 0, npack)
       do q = 1,2
          do k = 1,nlev
             elem(ie)%derived%FM(:,:,q,k) = elem(ie)%derived%FM(:,:,q,k)*elem(ie)%spheremp(:,:)
          end do
       end do
       do q = 1,2
          tmp = elem(ie)%derived%FM(:,:,q,:)
          call edgeVpack_nlyr(edge_g, elem(ie)%desc, tmp, nlev, (qsize+q-1)*nlev, npack)
       end do
       do k = 1,nlev
          elem(ie)%derived%FT(:,:,k) = elem(ie)%derived%FT(:,:,k)*elem(ie)%spheremp(:,:)
       end do
       call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FT, nlev, (qsize+2)*nlev, npack)
    end do
    if (hybrid%par%dynproc) call bndry_exchangeV(hybrid, edge_g)
    do ie = nets, nete
       call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FQ, qsize*nlev, 0, npack)
       do q = 1,qsize
          do k = 1,nlev
             elem(ie)%derived%FQ(:,:,k,q) = elem(ie)%derived%FQ(:,:,k,q)*elem(ie)%rspheremp(:,:)
          end do
       end do
       do q = 1,2
          tmp = elem(ie)%derived%FM(:,:,q,:)
          call edgeVunpack_nlyr(edge_g, elem(ie)%desc, tmp, nlev, (qsize+q-1)*nlev, npack)
          elem(ie)%derived%FM(:,:,q,:) = tmp
       end do
       do q = 1,2
          do k = 1,nlev
             elem(ie)%derived%FM(:,:,q,k) = elem(ie)%derived%FM(:,:,q,k)*elem(ie)%rspheremp(:,:)
          end do
       end do
       call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FT, nlev, (qsize+2)*nlev, npack)
       do k = 1,nlev
          elem(ie)%derived%FT(:,:,k) = elem(ie)%derived%FT(:,:,k)*elem(ie)%rspheremp(:,:)
       end do
    end do
  end subroutine gfr_f2g_dss

  subroutine gfr_f2g_remapd(gfr, gll_metdet, fv_metdet, f, g)
    ! Core remap operator. Conservative remap on the reference
    ! element.

    type (GllFvRemap_t), intent(in) :: gfr
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), fv_metdet(:), f(:)
    real(kind=real_kind), intent(out) :: g(:,:)

    integer :: nf, nf2, gi, gj, fi, fj
    real(kind=real_kind) :: wrk(np*np)

    nf = gfr%nphys
    nf2 = nf*nf
    wrk(:nf2) = f(:nf2)*fv_metdet(:nf2)
    do gj = 1,np
       do gi = 1,np
          g(gi,gj) = sum(gfr%f2g_remapd(:nf2,gi,gj)*wrk(:nf2))/ &
               gll_metdet(gi,gj)
       end do
    end do
  end subroutine gfr_f2g_remapd

  ! ----------------------------------------------------------------------
  ! All routines to boost pg1 OOA from 1 to ~1.6.

  subroutine gfr_pg1_init(gfr)
    ! Initializes pg1sd data structure; see gfr_pg1_solve doc for details.

    type (GllFvRemap_t), intent(inout) :: gfr

    real(kind=real_kind) :: Mnp2(np*np,4), wr(np*np)
    integer :: i, j, k, info, n

    if (gfr%nphys /= 1) return

    call make_mass_matrix_2d(np, np, gfr%pg1sd%Achol)
    call make_mass_matrix_2d(np, 2, Mnp2)

    n = np*np

    call dpotrf('U', n, gfr%pg1sd%Achol, size(gfr%pg1sd%Achol,1), info)
    if (info /= 0) print *, 'gfr ERROR> dpotrf returned', info

    do i = 1,n
       gfr%pg1sd%B(:,i) = Mnp2(i,:)
    end do

    ! Constraint vector c is just w_gg(gfr%pg1sd%inner), so don't store it explicitly.
    gfr%pg1sd%s = reshape(gfr%w_gg(:np,:np), (/np*np/))

    ! Form R's = c
    call dtrtrs('U', 'T', 'N', n, 1, gfr%pg1sd%Achol, size(gfr%pg1sd%Achol,1), &
         gfr%pg1sd%s, np*np, info)
    if (info /= 0) print *, 'gfr ERROR> dtrtrs returned', info
    gfr%pg1sd%sts = sum(gfr%pg1sd%s*gfr%pg1sd%s)
  end subroutine gfr_pg1_init

  subroutine gfr_pg1_solve(gfr, s, g)
    ! Solve
    !   min_g* M44 g* = M42 g(corners)  ! project a bilinear function onto 4-GLL basis
    !    st   spheremp'g* = spheremp g  ! subject to conserving mass
    ! where g(corners) extracts the corner values from the 4-GLL
    ! coefficient vector g. (In this explantation, for concreteness I
    ! assume np=4, but it doesn't have to be.) This problem
    ! reconstructs the field in the element using just the corner
    ! values and while conserving mass. This boosts the OOA of pg1.

    type (GllFvRemap_t), intent(in) :: gfr
    type (Pg1SolverData_t), intent(in) :: s
    real(kind=real_kind), intent(inout) :: g(:,:)

    real(kind=real_kind) :: x(np*np), mass, wr(4)
    integer :: np2, i, n, info

    np2 = np*np
    n = np2

    ! Form RHS M42 g(corners).
    wr(1) = g(1,1); wr(2) = g(np,1); wr(3) = g(1,np); wr(4) = g(np,np)
    do i = 1,n
       x(i) = sum(s%B(:4,i)*wr(:4))
    end do
    mass = sum(gfr%w_gg*g)

    ! Solve R'z = b.
    call dtrtrs('U', 'T', 'N', n, 1, s%Achol, size(s%Achol,1), x, np*np, info)
    ! Assemble z + (d - s'z)/(s's) s.
    x(:n) = x(:n) + ((mass - sum(s%s(:n)*x(:n)))/s%sts)*s%s(:n)
    ! Solve R x = z + (d - s'z)/(s's) s.
    call dtrtrs('U', 'N', 'N', n, 1, s%Achol, size(s%Achol,1), x, np*np, info)

    ! Extract g(I).
    g = reshape(x(:n), (/np,np/))
  end subroutine gfr_pg1_solve

  subroutine make_mass_matrix_2d(np1, np2, M, npq_in)
    ! Full mass matrix for a 2D element.

    use quadrature_mod, only : gausslobatto, quadrature_t

    integer, intent(in) :: np1, np2
    real(kind=real_kind), intent(out) :: M(:,:)
    integer, intent(in), optional :: npq_in

    type (quadrature_t) :: gll1, gll2, quad
    real(kind=real_kind) :: iv1(np1), jv1(np1), iv2(np2), jv2(np2), ir, jr
    integer :: np1sq, np2sq, npq, i1, j1, k1, i2, j2, k2, iq, jq

    npq = (np1 + np2 + 2)/2
    if (present(npq_in)) npq = npq_in

    np1sq = np1*np1
    np2sq = np2*np2

    gll1 = gausslobatto(np1)
    gll2 = gausslobatto(np2)
    quad = gausslobatto(npq)

    M(:np1sq,:np2sq) = zero

    do jq = 1,npq
       jr = quad%points(jq)
       call eval_lagrange_bases(gll1, np1, jr, jv1)
       call eval_lagrange_bases(gll2, np2, jr, jv2)
       do iq = 1,npq
          ir = quad%points(iq)
          call eval_lagrange_bases(gll1, np1, ir, iv1)
          call eval_lagrange_bases(gll2, np2, ir, iv2)
          do j2 = 1,np2
             do i2 = 1,np2
                k2 = (j2-1)*np2 + i2
                do j1 = 1,np1
                   do i1 = 1,np1
                      k1 = (j1-1)*np1 + i1
                      M(k1,k2) = M(k1,k2) + &
                           quad%weights(iq)*quad%weights(jq)* &
                           iv1(i1)*iv2(i2)*jv1(j1)*jv2(j2)
                   end do
                end do
             end do
          end do
       end do
    end do
    
    call gll_cleanup(quad)
    call gll_cleanup(gll2)
    call gll_cleanup(gll1)
  end subroutine make_mass_matrix_2d

  subroutine gfr_pg1_reconstruct_topo_hybrid(hybrid, elem, nets, nete)
    ! pg1 reconstruction routine for topography.

    use kinds, only: iulog
    use edge_mod, only: edgevpack_nlyr, edgevunpack_nlyr, edge_g
    use bndry_mod, only: bndry_exchangev

    type (hybrid_t), intent(in) :: hybrid
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete

    real(kind=real_kind) :: wr(np,np,2), ones(np,np,1)
    integer :: ie, nf, nerr

    if (gfr%nphys /= 1 .or. .not. gfr%boost_pg1) return

    ones = one

    do ie = nets,nete
       wr(:,:,1) = elem(ie)%state%phis
       call gfr_pg1_g_reconstruct_scalar(gfr, ie, elem(ie)%metdet, wr(:,:,:1))
       elem(ie)%state%phis = wr(:,:,1)
       call limiter_clip_and_sum(elem(ie)%spheremp, gfr%qmin(1,1,ie), &
            gfr%qmax(1,1,ie), ones(:,:,1), elem(ie)%state%phis)
       if (gfr%check > 1) then
          if (gfr%qmin(1,1,ie) < zero) then
             write(iulog,*) 'gfr> topo min:', hybrid%par%rank, hybrid%ithr, ie, &
                  gfr%qmin(1,1,ie), 'ERROR'
             gfr%check_ok = .false.
          end if
          wr(:,:,2) = elem(ie)%state%phis
          nerr = check_f2g_mixing_ratio(gfr, hybrid, ie, 1, elem, gfr%qmin(:1,1,ie), &
               gfr%qmax(:1,1,ie), ones, wr(:,:,1:1), wr(:,:,2:2))
          if (nerr > 0) gfr%check_ok = .false.
       end if
    end do

    if (hybrid%par%dynproc) then
       do ie = nets, nete
          elem(ie)%state%phis = elem(ie)%state%phis*elem(ie)%spheremp(:,:)
          call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis, 1, 0, 1)
       end do
       call bndry_exchangeV(hybrid, edge_g)
       do ie = nets, nete
          call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%phis, 1, 0, 1)
          elem(ie)%state%phis = elem(ie)%state%phis*elem(ie)%rspheremp(:,:)
       end do
    end if
  end subroutine gfr_pg1_reconstruct_topo_hybrid

  subroutine gfr_pg1_reconstruct_hybrid(hybrid, nt, hvcoord, elem, nets, nete)
    ! pg1 reconstruction routine for tendencies and states.

    use element_ops, only: get_field
    use dimensions_mod, only: nlev, qsize
    use hybvcoord_mod, only: hvcoord_t
    use physical_constants, only: p0, kappa

    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: nt
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete

    real(kind=real_kind), dimension(np,np,nlev) :: dp, p, wr1
    real(kind=real_kind) :: qmin, qmax
    integer :: ie, k, qi, nerr

    if (gfr%nphys /= 1 .or. .not. gfr%boost_pg1) return

    do ie = nets,nete
       dp = elem(ie)%state%dp3d(:,:,:,nt)

       call get_field(elem(ie), 'p', p, hvcoord, nt, -1)
       wr1 = (p0/p)**kappa
       elem(ie)%derived%FT = elem(ie)%derived%FT*wr1
       call gfr_pg1_g_reconstruct_scalar_dp(gfr, ie, elem(ie)%metdet, dp, &
            elem(ie)%derived%FT)
       elem(ie)%derived%FT = elem(ie)%derived%FT/wr1

       call gfr_pg1_g_reconstruct_vector(gfr, ie, elem, elem(ie)%derived%FM)

       do qi = 1,qsize
          ! Reconstruct Q_ten.
          call gfr_pg1_g_reconstruct_scalar_dp(gfr, ie, elem(ie)%metdet, dp, &
               elem(ie)%derived%FQ(:,:,:,qi))
          ! GLL Q.
          elem(ie)%derived%FQ(:,:,:,qi) = elem(ie)%derived%FQ(:,:,:,qi) + &
               elem(ie)%state%Q(:,:,:,qi)
          if (gfr%check > 1) wr1 = elem(ie)%derived%FQ(:,:,:,qi)
          do k = 1,nlev
             ! Augment bounds with GLL Q0 bounds.
             gfr%qmin(k,qi,ie) = min(minval(elem(ie)%state%Q(:,:,k,qi)), gfr%qmin(k,qi,ie))
             gfr%qmax(k,qi,ie) = max(maxval(elem(ie)%state%Q(:,:,k,qi)), gfr%qmax(k,qi,ie))
             call limiter_clip_and_sum(elem(ie)%spheremp, gfr%qmin(k,qi,ie), &
                  gfr%qmax(k,qi,ie), dp(:,:,k), elem(ie)%derived%FQ(:,:,k,qi))
          end do
          if (gfr%check > 1) then
             nerr = check_f2g_mixing_ratio(gfr, hybrid, ie, qi, elem, gfr%qmin(:,qi,ie), &
                  gfr%qmax(:,qi,ie), dp, wr1, elem(ie)%derived%FQ(:,:,:,qi))
             if (nerr > 0) gfr%check_ok = .false.
          end if
       end do
    end do

    call gfr_f2g_dss(hybrid, elem, nets, nete)
  end subroutine gfr_pg1_reconstruct_hybrid

  subroutine gfr_pg1_g_reconstruct_scalar(gfr, ie, gll_metdet, g)
    ! Wrapper to core solve routine. g is a density.

    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:)
    real(kind=real_kind), intent(inout) :: g(:,:,:)

    integer :: nlev, k

    nlev = size(g,3)
    do k = 1, nlev
       g(:,:,k) = g(:,:,k)*gll_metdet
       call gfr_pg1_solve(gfr, gfr%pg1sd, g(:,:,k))
       g(:,:,k) = g(:,:,k)/gll_metdet
    end do
  end subroutine gfr_pg1_g_reconstruct_scalar

  subroutine gfr_pg1_g_reconstruct_scalar_dp(gfr, ie, gll_metdet, dp, g)
    ! Wrapper to core solve routine. g is a mixing ratio.

    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), dp(:,:,:)
    real(kind=real_kind), intent(inout) :: g(:,:,:)

    real(kind=real_kind) wr(np*np)

    g = dp*g
    call gfr_pg1_g_reconstruct_scalar(gfr, ie, gll_metdet, g)
    g = g/dp
  end subroutine gfr_pg1_g_reconstruct_scalar_dp

  subroutine gfr_pg1_g_reconstruct_vector(gfr, ie, elem, v)
    ! Wrapper to core solve routine. v is a vector.

    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(inout) :: v(:,:,:,:)

    real(kind=real_kind) :: wr(np,np,2)
    integer :: nlev, k, d

    nlev = size(v,4)
    do k = 1,nlev
       ! sphere -> ref
       do d = 1,2
          wr(:,:,d) = elem(ie)%Dinv(:,:,d,1)*v(:,:,1,k) + elem(ie)%Dinv(:,:,d,2)*v(:,:,2,k)
       end do
       do d = 1,2
          call gfr_pg1_g_reconstruct_scalar(gfr, ie, elem(ie)%metdet, wr(:,:,d:d))
       end do
       ! ref -> sphere
       do d = 1,2
          v(:,:,d,k) = elem(ie)%D(:,:,d,1)*wr(:,:,1) + elem(ie)%D(:,:,d,2)*wr(:,:,2)
       end do
    end do
  end subroutine gfr_pg1_g_reconstruct_vector

  subroutine gfr_pg1_reconstruct_dom_mt(par, dom_mt, nt, hvcoord, elem)
    ! Wrapper to the hybrid-threading main routine.

    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use thread_mod, only: hthreads
    use hybvcoord_mod, only: hvcoord_t

    type (parallel_t), intent(in) :: par
    type (domain1d_t), intent(in) :: dom_mt(:)
    integer, intent(in) :: nt
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(inout) :: elem(:)

    type (hybrid_t) :: hybrid
    integer :: nets, nete
    
    if (.not. par%dynproc) return
#ifdef HORIZ_OPENMP
    !$omp parallel num_threads(hthreads), default(shared), private(nets,nete,hybrid)
#endif
    call gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)
    call gfr_pg1_reconstruct_hybrid(hybrid, nt, hvcoord, elem, nets, nete)
#ifdef HORIZ_OPENMP
    !$omp end parallel
#endif
  end subroutine gfr_pg1_reconstruct_dom_mt

  subroutine gfr_pg1_reconstruct_topo_dom_mt(par, dom_mt, elem)
    ! Wrapper to the hybrid-threading main routine.

    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use thread_mod, only: hthreads

    type (parallel_t), intent(in) :: par
    type (domain1d_t), intent(in) :: dom_mt(:)
    type (element_t), intent(inout) :: elem(:)

    type (hybrid_t) :: hybrid
    integer :: nets, nete
    
    if (.not. par%dynproc) return
#ifdef HORIZ_OPENMP
    !$omp parallel num_threads(hthreads), default(shared), private(nets,nete,hybrid)
#endif
    call gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)
    call gfr_pg1_reconstruct_topo_hybrid(hybrid, elem, nets, nete)
#ifdef HORIZ_OPENMP
    !$omp end parallel
#endif
  end subroutine gfr_pg1_reconstruct_topo_dom_mt

  ! ----------------------------------------------------------------------
  ! Internal helpers.

  subroutine gfr_dyn_to_fv_phys_topo_elem(elem, ie, phis)
    ! Remap topo within an element.

    type (element_t), intent(in) :: elem(:)
    integer, intent(in) :: ie
    real(kind=real_kind), intent(out) :: phis(:)

    call gfr_g2f_scalar_and_limit(gfr, ie, elem(ie)%metdet, elem(ie)%state%phis, phis)
  end subroutine gfr_dyn_to_fv_phys_topo_elem

  subroutine gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)
    ! Create a hybrid_t given a dom_mt.

    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use hybrid_mod, only: hybrid_create
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

  subroutine apply_interp(interp, np, npi, x, y)
    ! Apply the npi -> np interpolation matrix 'interp' to x to get y.

    real(kind=real_kind), intent(in) :: interp(:,:,:,:), x(:,:)
    integer, intent(in) :: np, npi
    real(kind=real_kind), intent(out) :: y(:,:)

    integer :: gi, gj, fi, fj, info
    real(kind=real_kind) :: accum

    do fj = 1,np
       do fi = 1,np
          accum = zero
          do gj = 1,npi
             do gi = 1,npi
                accum = accum + gfr%interp(gi,gj,fi,fj)*x(gi,gj)
             end do
          end do
          y(fi,fj) = accum
       end do
    end do
  end subroutine apply_interp

  subroutine gfr_g_make_nonnegative(gll_metdet, g)
    ! Move mass around as needed to make g nonnegative, where g is a
    ! density.

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

  function gfr_f_get_area(ie, i, j) result(area)
    ! Get (lat,lon) of FV point i,j.

    integer, intent(in) :: ie, i, j
    real(kind=real_kind) :: area
    
    integer :: k

    k = gfr%nphys*(j-1) + i
    area = gfr%w_ff(k)*gfr%fv_metdet(k,ie)
  end function gfr_f_get_area

  subroutine gfr_f_get_latlon(ie, i, j, lat, lon)
    ! Get (lat,lon) of FV point i,j.

    use coordinate_systems_mod, only: spherical_polar_t, change_coordinates

    integer, intent(in) :: ie, i, j
    real(kind=real_kind), intent(out) :: lat, lon

    type (spherical_polar_t) :: p

    p = change_coordinates(gfr%center_f(i,j,ie))
    lat = p%lat
    lon = p%lon
  end subroutine gfr_f_get_latlon

  subroutine gfr_f_get_corner_latlon(ie, i, j, c, lat, lon)
    ! Get (lat,lon) of FV point i,j.

    use coordinate_systems_mod, only: spherical_polar_t, change_coordinates

    integer, intent(in) :: ie, i, j, c
    real(kind=real_kind), intent(out) :: lat, lon

    type (spherical_polar_t) :: p

    p = change_coordinates(gfr%corners_f(c,i,j,ie))
    lat = p%lat
    lon = p%lon
  end subroutine gfr_f_get_corner_latlon

  subroutine gfr_f_get_cartesian3d(ie, i, j, p)
    ! Get (x,y,z) of FV point i,j.

    integer, intent(in) :: ie, i, j
    type (cartesian3D_t), intent(out) :: p

    p = gfr%center_f(i,j,ie)
  end subroutine gfr_f_get_cartesian3d

  subroutine calc_dp_fv(nf, hvcoord, ps, dp_fv)
    ! Compute pressure level increments on the FV grid given ps on the FV
    ! grid. Directly projecting dp_gll to dp_fv disagrees numerically with the
    ! loop in this subroutine. This loop is essentially how CAM computes pdel,
    ! so we must use it, too.

    use hybvcoord_mod, only: hvcoord_t
    use dimensions_mod, only: nlev
    
    integer, intent(in) :: nf
    type (hvcoord_t), intent(in) :: hvcoord
    real (kind=real_kind), intent(in) :: ps(:)
    real (kind=real_kind), intent(out) :: dp_fv(:,:)

    integer :: k, nf2

    nf2 = nf*nf
    do k = 1,nlev
       dp_fv(:nf2,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
                       (hvcoord%hybi(k+1) - hvcoord%hybi(k))*ps(:nf2)
    end do
  end subroutine calc_dp_fv

  subroutine limiter1_clip_and_sum(n, spheremp, qmin, qmax, dp, q)
    ! CAAS as described in Alg 3.1 of doi:10.1137/18M1165414. q is a
    ! mixing ratio. Solve
    !    min_q* norm(dp q - dp q*, 1)
    !     st    spheremp'(dp q*) = spheremp'(dp q)
    !           qmin < q* < qmax

    integer, intent(in) :: n
    real (kind=real_kind), intent(in) :: spheremp(n), dp(n)
    real (kind=real_kind), intent(inout) :: qmin, qmax, q(n)

    integer :: n2, k1, i, j
    logical :: modified
    real(kind=real_kind) :: addmass, mass, sumc, den
    real(kind=real_kind) :: x(n), c(n), v(n)

    x = q(:n)
    c = spheremp(:n)*dp(:n)

    sumc = sum(c)
    mass = sum(c*x)
    ! In the case of an infeasible problem, prefer to conserve mass
    ! and violate a bound.
    if (mass < qmin*sumc) qmin = mass / sumc
    if (mass > qmax*sumc) qmax = mass / sumc

    addmass = zero

    ! Clip.
    modified = .false.
    do k1 = 1, n
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

    q(:n) = x
  end subroutine limiter1_clip_and_sum

  subroutine limiter_clip_and_sum(spheremp, qmin, qmax, dp, q)
    real (kind=real_kind), intent(in) :: spheremp(np,np), dp(np,np)
    real (kind=real_kind), intent(inout) :: qmin, qmax, q(np,np)

    call limiter1_clip_and_sum(np*np, spheremp, qmin, qmax, dp, q)
  end subroutine limiter_clip_and_sum

  subroutine ref2spherea_deriv(c, a, b, s_ab, s)
    ! For cubed_sphere_map = 2.

    real(real_kind), intent(in) :: c(4,3), a, b
    real(real_kind), intent(out) :: s_ab(3,2), s(3)

    real(real_kind) :: q(4), q_ab(4,2), r2
    integer :: i, j

    q(1) = (1-a)*(1-b); q(2) = (1+a)*(1-b); q(3) = (1+a)*(1+b); q(4) = (1-a)*(1+b)
    q = q/four
    s = zero
    do i = 1,3
       s(i) = sum(c(:,i)*q)
    end do
    r2 = sum(s**2)
    q_ab(1,1) = -(1-b); q_ab(2,1) =  (1-b); q_ab(3,1) = (1+b); q_ab(4,1) = -(1+b)
    q_ab(1,2) = -(1-a); q_ab(2,2) = -(1+a); q_ab(3,2) = (1+a); q_ab(4,2) =  (1-a)
    q_ab = q_ab/four
    s_ab = zero
    do j = 1,2
       do i = 1,3
          s_ab(i,j) = sum(c(:,i)*q_ab(:,j))
       end do
       s_ab(:,j) = s_ab(:,j)/sqrt(r2) - (sum(s*s_ab(:,j))/sqrt(r2**3))*s
    end do
    s = s/sqrt(r2)
  end subroutine ref2spherea_deriv

  subroutine sphere2ref(corners, sphere, a, b, tol_in, maxit_in)
    ! For cubed_sphere_map = 2.

    use coordinate_systems_mod, only: cartesian3D_t

    type (cartesian3D_t), intent(in) :: corners(4), sphere
    real(real_kind), intent(out) :: a, b
    real(real_kind), intent(in), optional :: tol_in
    integer, intent(in), optional :: maxit_in

    real(real_kind) :: tol, c(4,3), s_in(3), s_ab(3,2), s(3), r(3), fac(3), x(2)
    integer :: maxit, i, it

    tol = eps
    maxit = 10
    if (present(tol_in)) tol = tol_in
    if (present(maxit_in)) maxit = maxit_in
    tol = tol**2

    do i = 1,4
       c(i,1) = corners(i)%x; c(i,2) = corners(i)%y; c(i,3) = corners(i)%z
    end do
    s_in(1) = sphere%x; s_in(2) = sphere%y; s_in(3) = sphere%z

    a = zero; b = zero
    do it = 1,maxit
       call ref2spherea_deriv(c, a, b, s_ab, s)
       r = s - s_in
       if (sum(r**2) <= tol) exit
       !! QR for s_ab d_ab = r
       ! Q
       fac(1) = sqrt(sum(s_ab(:,1)**2))
       s_ab(:,1) = s_ab(:,1)/fac(1)
       fac(2) = sum(s_ab(:,2)*s_ab(:,1))
       s_ab(:,2) = s_ab(:,2) - fac(2)*s_ab(:,1)
       fac(3) = sqrt(sum(s_ab(:,2)**2))
       s_ab(:,2) = s_ab(:,2)/fac(3)
       ! x = Q'r
       x(1) = sum(s_ab(:,1)*r)
       x(2) = sum(s_ab(:,2)*r)
       ! x = R \ x
       x(2) = x(2) / fac(3)
       x(1) = (x(1) - fac(2)*x(2)) / fac(1)
       !! Newton update
       a = a - x(1)
       b = b - x(2)
    end do
  end subroutine sphere2ref

  ! ----------------------------------------------------------------------
  ! Everything below is for internal unit testing of this module. For
  ! integration-level testing, see gllfvremap_util_mod and
  ! dcmip2016_test1_pg_forcing.

  subroutine set_ps_Q(elem, nets, nete, timeidx, qidx, nlev)
    ! Make up a test function for use in unit tests.

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

  function check_global_properties(gfr, hybrid, hvcoord, elem, nt, nets, nete, &
       use_state_Q, q_f) result(nerr)

    ! Compare global mass on dynamics and physics grids.

    use parallel_mod, only: global_shared_buf, global_shared_sum, nrepro_vars
    use global_norms_mod, only: wrap_repro_sum
    use kinds, only: iulog
    use dimensions_mod, only: nlev, qsize_d, qsize
    use hybvcoord_mod, only: hvcoord_t

    type (GllFvRemap_t), intent(in) :: gfr
    type (hybrid_t), intent(in) :: hybrid
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(in) :: elem(:)
    integer, intent(in) :: nt, nets, nete
    logical, intent(in) :: use_state_Q
    real (kind=real_kind), intent(in) :: q_f(:,:,:,:)

    integer :: nf, nf2, ie, k, qi, ic, nchunk, qi0, nq, qic, b1, b2, cnt, nerr
    real (kind=real_kind) :: dp(np,np,nlev), dp_fv(np*np,nlev), wg(np,np), &
         wf(np*np,1), mass(2,qsize), tol

    nerr = 0
    nf = gfr%nphys
    nf2 = nf*nf
    nq = nrepro_vars/2
    nchunk = (qsize*2 + 2*nq - 1)/(2*nq)
    do ic = 1,nchunk
       qi0 = nq*(ic-1)
       do ie = nets,nete
          dp = elem(ie)%state%dp3d(:,:,:,nt)
          call gfr_g2f_scalar(ie, elem(ie)%metdet, elem(ie)%state%ps_v(:,:,nt:nt), &
               wf(:,:1))
          call calc_dp_fv(nf, hvcoord, wf(:,1), dp_fv)
          global_shared_buf(ie,:) = 0
          cnt = 0
          do qic = 1,nq
             qi = qi0 + qic
             if (qi > qsize) exit
             cnt = cnt + 2
             b1 = 2*(qic-1) + 1
             b2 = 2*(qic-1) + 2
             do k = 1,nlev
                if (use_state_Q) then
                   wg = elem(ie)%state%Q(:,:,k,qi)
                else
                   wg = elem(ie)%derived%FQ(:,:,k,qi)
                end if
                global_shared_buf(ie,b1) = global_shared_buf(ie,b1) + &
                     sum(elem(ie)%spheremp(:,:)*dp(:,:,k)*wg)
                global_shared_buf(ie,b2) = global_shared_buf(ie,b2) + &
                     sum(gfr%fv_metdet(:nf2,ie)*gfr%w_ff(:nf2)* &
                     dp_fv(:nf2,k)*q_f(:nf2,k,qi,ie))
             end do
          end do
       end do
       call wrap_repro_sum(nvars=cnt, comm=hybrid%par%comm)
       do qic = 1,nq
          qi = qi0 + qic
          if (qi > qsize) exit
          b1 = 2*(qic-1) + 1
          b2 = 2*(qic-1) + 2
          mass(1,qi) = global_shared_sum(b1)
          mass(2,qi) = global_shared_sum(b2)
       end do
    end do
    tol = 10*eps
    if (hybrid%masterthread) then
       do qi = 1,qsize
          if (abs(mass(2,qi) - mass(1,qi)) > tol*abs(mass(1,qi))) then
             nerr = nerr + 1
             write (iulog,'(a,l2,i3,es24.16,es12.4)') 'gfr> mass err', &
                  use_state_Q, qi, mass(1,qi), &
                  abs(mass(2,qi) - mass(1,qi))/maxval(abs(mass(1:2,qi)))
          end if
       end do
    end if
  end function check_global_properties

  function check_g2f_mixing_ratio(gfr, hybrid, ie, qi, elem, dp, dp_fv, q_g, q_f) result(nerr)
    ! Check that gfr_g2f_mixing_ratio found a property-preserving
    ! solution.

    use kinds, only: iulog

    type (GllFvRemap_t), intent(in) :: gfr
    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: ie, qi
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(in) :: dp(:,:,:), dp_fv(:,:), q_g(:,:,:), q_f(:,:)

    real(kind=real_kind) :: qmin_f, qmin_g, qmax_f, qmax_g, mass_f, mass_g, den
    integer :: q, k, nf, nf2, nerr

    nerr = 0
    nf = gfr%nphys
    nf2 = nf*nf
    do k = 1,size(dp,3)
       qmin_f = minval(q_f(:nf2,k))
       qmax_f = maxval(q_f(:nf2,k))
       qmin_g = minval(elem(ie)%state%Q(:,:,k,qi))
       qmax_g = maxval(elem(ie)%state%Q(:,:,k,qi))
       den = gfr%tolfac*max(1e-10_real_kind, maxval(abs(elem(ie)%state%Q(:,:,k,qi))))
       mass_f = sum((gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie))*dp_fv(:nf2,k)*q_f(:nf2,k))
       mass_g = sum(elem(ie)%spheremp*dp(:,:,k)*q_g(:,:,k))
       if (qmin_f < qmin_g - 10*eps*den .or. qmax_f > qmax_g + 10*eps*den) then
          write(iulog,*) 'gfr> g2f mixing ratio limits:', hybrid%par%rank, hybrid%ithr, ie, qi, k, &
               qmin_g, qmin_f-qmin_g, qmax_f-qmax_g, qmax_g, mass_f, mass_g, 'ERROR'
          nerr = nerr + 1
       end if
       if (abs(mass_f - mass_g) > gfr%tolfac*20*eps*max(mass_f, mass_g)) then
          write(iulog,*) 'gfr> g2f mixing ratio mass:', hybrid%par%rank, hybrid%ithr, ie, qi, k, &
               qmin_g, qmax_g, mass_f, mass_g, 'ERROR'
          nerr = nerr + 1
       end if
    end do
  end function check_g2f_mixing_ratio

  function check_f2g_mixing_ratio(gfr, hybrid, ie, qi, elem, qmin, qmax, dp, q0_g, q1_g) result(nerr)
    ! Check that a property-preserving solution was found in the FV ->
    ! GLL direction.

    use kinds, only: iulog

    type (GllFvRemap_t), intent(in) :: gfr
    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: ie, qi
    type (element_t), intent(in) :: elem(:)
    real(kind=real_kind), intent(in) :: qmin(:), qmax(:), dp(:,:,:), q0_g(:,:,:), q1_g(:,:,:)

    real(kind=real_kind) :: qmin_f, qmin_g, qmax_f, qmax_g, mass_f, mass0, mass1, den, &
         wr(np,np)
    integer :: q, k, nerr

    nerr = 0
    do k = 1,size(dp,3)
       qmin_f = qmin(k)
       qmax_f = qmax(k)
       qmin_g = minval(q1_g(:,:,k))
       qmax_g = maxval(q1_g(:,:,k))
       den = gfr%tolfac*max(1e-10_real_kind, maxval(abs(q0_g(:,:,k))))
       mass0 = sum(elem(ie)%spheremp*dp(:,:,k)*q0_g(:,:,k))
       mass1 = sum(elem(ie)%spheremp*dp(:,:,k)*q1_g(:,:,k))
       if (qmin_g < qmin_f - 50*eps*den .or. qmax_g > qmax_f + 50*eps*den) then
          write(iulog,*) 'gfr> f2g mixing ratio limits:', hybrid%par%rank, hybrid%ithr, ie, qi, k, &
               qmin_f, qmin_g-qmin_f, qmax_g-qmax_f, qmax_f, mass0, mass1, 'ERROR'
          nerr = nerr + 1
       end if
       den = sum(elem(ie)%spheremp*dp(:,:,k)*maxval(abs(q0_g(:,:,k))))
       if (abs(mass1 - mass0) > gfr%tolfac*20*eps*den) then
          write(iulog,*) 'gfr> f2g mixing ratio mass:', hybrid%par%rank, hybrid%ithr, ie, qi, k, &
               qmin_f, qmin_g, qmax_g, qmax_f, mass0, mass1, 'ERROR'
          nerr = nerr + 1
       end if
    end do
  end function check_f2g_mixing_ratio
  
  function check_nonnegative(elem, nets, nete) result(nerr)
    ! Check gfr_g_make_nonnegative.

    use kinds, only: iulog

    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete

    real(kind=real_kind) :: wrk3(np,np,1), mass0, mass1, rd
    integer :: ie, i, j, sign, nerr

    nerr = 0
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
          nerr = nerr + 1
       end if
    end do
  end function check_nonnegative

  subroutine check_areas(par, gfr, elem, nets, nete)
    ! Check global area

    use kinds, only: iulog
    use parallel_mod, only: parallel_t, global_shared_buf, global_shared_sum
    use physical_constants, only: dd_pi
    ! Can't use wrap_repro_sum because this routine needs to support
    ! unit testing in an already threaded region.
#ifdef CAM
    use shr_reprosum_mod, only: repro_sum => shr_reprosum_calc
#else
    use repro_sum_mod, only: repro_sum
#endif

    type (parallel_t), intent(in) :: par
    type (GllFvRemap_t), intent(inout) :: gfr
    type (element_t), intent(in) :: elem(:)
    integer, intent(in) :: nets, nete

    integer :: ie, i, j, nf
    real(kind=real_kind) :: area, sphere_area, re

    nf = gfr%nphys
    do ie = nets,nete
       global_shared_buf(ie,1) = gfr%check_areas(1,ie)
       area = zero
       do j = 1,nf
          do i = 1,nf
             area = area + gfr_f_get_area(ie, i, j)
          end do
       end do
       global_shared_buf(ie,2) = area
       global_shared_buf(ie,3) = sum(elem(ie)%spheremp)
    end do
    call repro_sum(global_shared_buf, global_shared_sum, nelemd, nelemd, 3, commid=par%comm)
    sphere_area = 4*dd_pi
    if (par%masterproc) then
       write(iulog,*) 'gfr> area fv raw', global_shared_sum(1), &
            abs(global_shared_sum(1) - sphere_area)/sphere_area
       ! fv vs gll
       re = abs(global_shared_sum(2) - global_shared_sum(3))/global_shared_sum(3)
       write(iulog,*) 'gfr> area fv adj', &
            abs(global_shared_sum(2) - sphere_area)/sphere_area, re
       if (re > 2*eps) then
          write(iulog,*) 'gfr> check_areas ERROR'
          gfr%check_ok = .false.
       end if
       write(iulog,*) 'gfr> area gll   ', &
            abs(global_shared_sum(3) - sphere_area)/sphere_area
    end if
    deallocate(gfr%check_areas)
  end subroutine check_areas

  function check(par, dom_mt, gfr, elem, verbose) result(nerr)
    ! Run a bunch of unit tests.

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

    real(kind=real_kind) :: a, b, rd, x, y, f0(np*np), f1(np*np), g(np,np), &
         wf(np*np), wg(np,np), qmin, qmax, qmin1, qmax1
    integer :: nf, nf2, ie, i, j, k, iremap, info, ilimit, it
    real(kind=real_kind), allocatable :: Qdp_fv(:,:), ps_v_fv(:,:), &
         qmins(:,:,:), qmaxs(:,:,:)
    logical :: limit
    character(32) :: msg

    ! Purposely construct our own hybrid object to test gfr_hybrid_create.
    type (hybrid_t) :: hybrid
    integer :: nets, nete, nerr, ic

    nerr = 0
    if (.not. par%dynproc) return

    nf = gfr%nphys
    nf2 = nf*nf

    call gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)

    if (hybrid%masterthread) then
       write(iulog,  '(a,i3,a,i3)') 'gfr> npi', gfr%npi, ' nphys', nf
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
       if (rd /= rd .or. rd > 10*eps) then
          nerr = nerr + 1
          write(iulog,*) 'gfr> area', ie, a, b, rd
       end if

       ! Check FV geometry.
       f0(:nf2) = gfr%D_f(:,1,1,ie)*gfr%D_f(:,2,2,ie) - &
            gfr%D_f(:,1,2,ie)*gfr%D_f(:,2,1,ie)
       rd = maxval(abs(f0(:nf2)) - gfr%fv_metdet(:nf2,ie))/ &
            maxval(gfr%fv_metdet(:nf2,ie))
       if (rd > 10*eps) then
          nerr = nerr + 1
          write(iulog,*) 'gfr> D', ie, rd
       end if
       f0(:nf2) = gfr%Dinv_f(:,1,1,ie)*gfr%Dinv_f(:,2,2,ie) - &
            gfr%Dinv_f(:,1,2,ie)*gfr%Dinv_f(:,2,1,ie)
       rd = maxval(abs(f0(:nf2)) - one/gfr%fv_metdet(:nf2,ie))/ &
            maxval(one/gfr%fv_metdet(:nf2,ie))
       if (rd > 10*eps) then
          nerr = nerr + 1
          write(iulog,*) 'gfr> Dinv', ie, rd
       end if

       ! Check that FV -> GLL -> FV recovers the original FV values exactly
       ! (with no DSS and no limiter).
       do j = 1,nf
          x = real(j-1, real_kind)/real(nf, real_kind)
          do i = 1,nf
             y = real(i-1, real_kind)/real(nf, real_kind)
             k = i + (j-1)*nf
             f0(k) = real(ie)/nelemd + x*x + ie*x + cos(ie + 4.2*y)
          end do
       end do
       call gfr_f2g_remapd(gfr, elem(ie)%metdet, gfr%fv_metdet(:,ie), f0, g)
       call gfr_g2f_remapd(gfr, elem(ie)%metdet, gfr%fv_metdet(:,ie), g, f1)
       wf(:nf2) = gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie)
       a = sum(wf(:nf2)*abs(f1(:nf2) - f0(:nf2)))
       b = sum(wf(:nf2)*abs(f0(:nf2)))
       rd = a/b
       if (rd /= rd .or. rd > 10*eps) then
          nerr = nerr + 1
          write(iulog,*) 'gfr> recover', ie, a, b, rd, gfr%fv_metdet(:nf2,ie)
       end if
    end do
    nerr = nerr + check_nonnegative(elem, nets, nete)

    ! For convergence testing. Run this testing routine with a sequence of ne
    ! values and plot log l2 error vs log ne.
    allocate(Qdp_fv(nf2, nets:nete), ps_v_fv(nf2, nets:nete))
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
                  elem(ie)%state%ps_v(:,:,1)*elem(ie)%state%Q(:,:,1,1), Qdp_fv(:,ie))
             call gfr_g2f_remapd(gfr, elem(ie)%metdet, gfr%fv_metdet(:,ie), &
                  elem(ie)%state%ps_v(:,:,1), ps_v_fv(:,ie))
             if (limit) then
                qmin = minval(elem(ie)%state%Q(:,:,1,1))
                qmax = maxval(elem(ie)%state%Q(:,:,1,1))
                wf(:nf2) = Qdp_fv(:nf2,ie)/ps_v_fv(:nf2,ie)
                f0(:nf2) = gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie)
                call limiter1_clip_and_sum(nf2, f0, qmin, qmax, ps_v_fv(:,ie), wf)
                Qdp_fv(:nf2,ie) = wf(:nf2)*ps_v_fv(:nf2,ie)
             end if
          end do
          ! 2. FV -> GLL
          if (limit) then
             ! 2a. Get q bounds
             do ie = nets, nete
                wf(:nf2) = Qdp_fv(:nf2,ie)/ps_v_fv(:nf2,ie)
                qmins(:,:,ie) = minval(wf(:nf2))
                qmaxs(:,:,ie) = maxval(wf(:nf2))
             end do
             ! 2b. Halo exchange q bounds.
             call neighbor_minmax(hybrid, edgeAdvQminmax, nets, nete, qmins, qmaxs)
             ! 2c. Augment bounds with current values.
             do ie = nets, nete
                wg = elem(ie)%state%Q(:,:,1,1)
                qmins(1,1,ie) = min(qmins(1,1,ie), minval(wg))
                qmaxs(1,1,ie) = max(qmaxs(1,1,ie), maxval(wg))                
             end do
          endif
          ! 2d. Remap
          if (nf == 1 .and. gfr%boost_pg1) then
             do ie = nets, nete
                elem(ie)%state%Q(:,:,1,1) = Qdp_fv(1,ie)/ps_v_fv(1,ie)
             end do
          else
             do ie = nets, nete
                call gfr_f2g_remapd(gfr, elem(ie)%metdet, gfr%fv_metdet(:,ie), &
                     Qdp_fv(:,ie), elem(ie)%state%Q(:,:,1,1))
                elem(ie)%state%Q(:,:,1,1) = elem(ie)%state%Q(:,:,1,1)/elem(ie)%state%ps_v(:,:,1)
                if (limit) then
                   call limiter_clip_and_sum(elem(ie)%spheremp, & ! same as w_gg*gll_metdet
                        qmins(1,1,ie), qmaxs(1,1,ie), elem(ie)%state%ps_v(:,:,1), &
                        elem(ie)%state%Q(:,:,1,1))
                end if
             end do
          end if
          do it = 1,2
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
             if (it == 2 .or. nf > 1) exit
             ! 4. pg1 OOA boost.
             if (nf == 1 .and. gfr%boost_pg1) then
                do ie = nets, nete
                   if (limit) then
                      qmins(1,1,ie) = min(minval(elem(ie)%state%Q(:,:,1,1)), qmins(1,1,ie))
                      qmaxs(1,1,ie) = max(maxval(elem(ie)%state%Q(:,:,1,1)), qmaxs(1,1,ie))
                   end if
                   call gfr_pg1_g_reconstruct_scalar_dp(gfr, ie, elem(ie)%metdet, &
                        elem(ie)%state%ps_v(:,:,:1), elem(ie)%state%Q(:,:,:1,1))
                   if (limit) then
                      call limiter_clip_and_sum(gfr%w_gg*elem(ie)%metdet, qmins(1,1,ie), &
                           qmaxs(1,1,ie), elem(ie)%state%ps_v(:,:,1), elem(ie)%state%Q(:,:,1,1))
                   end if
                end do
             end if
          end do
       end do
       ! 5. Compute error.
       qmin = two
       qmax = -two
       qmin1 = two
       qmax1 = -two
       ic = 2
       do ie = nets, nete
          wg = gfr%w_gg(:,:)*elem(ie)%metdet(:,:)
          ! L2 on q. Might switch to q*ps_v.
          global_shared_buf(ie,1) = &
               sum(wg*(elem(ie)%state%Q(:,:,1,1) - elem(ie)%state%Q(:,:,1,ic))**2)
          global_shared_buf(ie,2) = &
               sum(wg*elem(ie)%state%Q(:,:,1,ic)**2)
          ! Mass conservation.
          wg = wg*elem(ie)%state%ps_v(:,:,1)
          global_shared_buf(ie,3) = sum(wg*elem(ie)%state%Q(:,:,1,ic))
          global_shared_buf(ie,4) = sum(wg*elem(ie)%state%Q(:,:,1,1))
          qmin = min(qmin, minval(elem(ie) %state%Q(:,:,1,1)))
          qmin1 = min(qmin1, minval(elem(ie)%state%Q(:,:,1,ic)))
          qmax = max(qmax, maxval(elem(ie)%state%Q(:,:,1,1)))
          qmax1 = max(qmax1, maxval(elem(ie)%state%Q(:,:,1,ic)))
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
          if (rd > 10*eps) then
             nerr = nerr + 1
             msg = ' ERROR'
          end if
          write(iulog, '(a,es11.3,a8)') 'gfr> mass', rd, msg
          msg = ''
          if (limit .and. (qmin < qmin1 - 5*eps .or. qmax > qmax1 + 5*eps)) then
             nerr = nerr + 1
             msg = ' ERROR'
          end if
          write(iulog, '(a,es11.3,es11.3,a8)') 'gfr> limit', min(zero, qmin - qmin1), &
               max(zero, qmax - qmax1), msg
       end if
    end do
    deallocate(Qdp_fv, ps_v_fv, qmins, qmaxs)
    if (.not. gfr%check_ok) nerr = nerr + 1
  end function check

  function test_sphere2ref() result(nerr)
    use coordinate_systems_mod, only: cartesian3D_t
    use kinds, only: iulog

    type (cartesian3D_t) :: corners(4), sphere
    real (real_kind) :: refin(2), refout(2), err
    integer :: i, j, n, nerr

    nerr = 0

    corners(1)%x =  0.24; corners(1)%y = -0.7; corners(1)%z = 0.3; call normalizecart(corners(1))
    corners(2)%x =  0.44; corners(2)%y =  0.5; corners(2)%z = 0.4; call normalizecart(corners(2))
    corners(3)%x = -0.34; corners(3)%y =  0.6; corners(3)%z = 0.1; call normalizecart(corners(3))
    corners(4)%x = -0.14; corners(4)%y = -0.5; corners(4)%z = 0.2; call normalizecart(corners(4))

    n = 77
    nerr = 0
    do i = 1,n
       refin(1) = -1 + (1.0_real_kind/(n-1))*(i-1)
       do j = 1,n
          refin(2) = -1 + (1.0_real_kind/(n-1))*(j-1)
          call ref2sphere(corners, refin(1), refin(2), sphere)
          call sphere2ref(corners, sphere, refout(1), refout(2))
          err = abs(refin(1) - refout(1)) + abs(refin(2) - refout(2))
          if (err > 15*eps .or. &
               maxval(abs(refout)) > 1 + 5*eps .or. &
               any(refout /= refout)) then
             write(iulog,*) refin(1), refin(2)
             write(iulog,*) refout(1), refout(2)
             write(iulog,*) err
             nerr = nerr + 1
          end if
       end do
    end do
    if (nerr /= 0) write(iulog,*) 'test_sphere2ref FAILED'

  contains
    subroutine normalizecart(sphere)
      type (cartesian3D_t), intent(inout) :: sphere
      real(real_kind) :: r
      r = sqrt(sphere%x**2 + sphere%y**2 + sphere%z**2)
      sphere%x = sphere%x/r; sphere%y = sphere%y/r; sphere%z = sphere%z/r
    end subroutine normalizecart

    subroutine ref2sphere(corners, a, b, sphere)
      type (cartesian3D_t), intent(in) :: corners(4)
      real(real_kind), intent(in) :: a, b
      type (cartesian3D_t), intent(out) :: sphere

      real(real_kind) :: c(4,3), s(3)
      integer :: i

      do i = 1,4
         c(i,1) = corners(i)%x; c(i,2) = corners(i)%y; c(i,3) = corners(i)%z
      end do
      call ref2spherea(c, a, b, s)
      sphere%x = s(1); sphere%y = s(2); sphere%z = s(3)
    end subroutine ref2sphere

    subroutine ref2spherea(c, a, b, s)
      real(real_kind), intent(in) :: c(4,3), a, b
      real(real_kind), intent(out) :: s(3)

      real(real_kind) :: q(4)
      integer :: i

      q(1) = (1-a)*(1-b); q(2) = (1+a)*(1-b); q(3) = (1+a)*(1+b); q(4) = (1-a)*(1+b)
      q = q/four
      s = zero
      do i = 1,3
         s(i) = sum(c(:,i)*q)
      end do
      s = s/sqrt(sum(s**2))
    end subroutine ref2spherea
  end function test_sphere2ref

  function gfr_test(hybrid, dom_mt, hvcoord, deriv, elem) result(nerr)
    ! Driver for check subroutine.

    use domain_mod, only: domain1d_t
    use derivative_mod, only: derivative_t
    use hybvcoord_mod, only: hvcoord_t

    type (hybrid_t), intent(in) :: hybrid
    type (domain1d_t), intent(in) :: dom_mt(:)
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(inout) :: elem(:)
    type (hvcoord_t) , intent(in) :: hvcoord

    integer :: nphys, bi, nerr
    logical :: boost_pg1

    nerr = 0    
    if (hybrid%masterthread) nerr = nerr + test_sphere2ref()
    do nphys = 1, np
       do bi = 1,2
          if (nphys > 1 .and. bi > 1) exit
          boost_pg1 = bi == 2

          ! This is meant to be called before threading starts.
          if (hybrid%ithr == 0) call gfr_init(hybrid%par, elem, nphys, 2, boost_pg1)
          !$omp barrier

          nerr = nerr + check(hybrid%par, dom_mt, gfr, elem, .false.)

          ! This is meant to be called after threading ends.
          !$omp barrier
          if (hybrid%ithr == 0) call gfr_finish()
          !$omp barrier
       end do
    end do
  end function gfr_test
end module gllfvremap_mod
