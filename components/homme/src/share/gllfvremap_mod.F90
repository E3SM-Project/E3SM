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
  !   This module supports all ftype values 0 to 4.
  !   We find in practice that pg1 is too coarse. To see this, run the
  ! Homme-standalone test dcmip2016_test1_pg1 and compare results with
  ! dcmip2016_test1_pg2 and dcmip2016_test1 (np4). pg2 and np4 fields are nearly
  ! identical out to day 30, whereas pg1 fields differ visibly.
  !
  ! AMB 2019/07 Initial

  use hybrid_mod, only: hybrid_t
  use kinds, only: real_kind
  use dimensions_mod, only: np, npsq, qsize, nelemd
  use element_mod, only: element_t
  use coordinate_systems_mod, only: spherical_polar_t

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
     integer :: nphys, npi
     logical :: check, have_fv_topo_file_phis, boost_pg1
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
     type (Pg1SolverData_t) :: pg1sd
  end type GllFvRemap_t

  type (GllFvRemap_t), private :: gfr

  ! For testing.
  public :: &
       gfr_test, &
       gfr_g2f_scalar, gfr_f2g_scalar, gfr_f_get_latlon, gfr_f_get_cartesian3d, &
       gfr_g_make_nonnegative, gfr_dyn_to_fv_phys_topo_elem, gfr_f2g_dss

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
    !   check is optional and defaults to false. It will produce very
    ! verbose output if something goes wrong. It is also expensive. It
    ! is intended to be used in unit testing and (if ever needed) for
    ! debuggin.

    use kinds, only: iulog
    use dimensions_mod, only: nlev
    use parallel_mod, only: parallel_t, abortmp
    use quadrature_mod, only : gausslobatto, quadrature_t
    use control_mod, only: ftype

    type (parallel_t), intent(in) :: par
    type (element_t), intent(in) :: elem(:)
    integer, intent(in) :: nphys
    logical, intent(in), optional :: check, boost_pg1

    real(real_kind) :: R(npsq,nphys_max*nphys_max), tau(npsq)
    integer :: nphys2

    gfr%check = .false.
    if (present(check)) gfr%check = check

    gfr%boost_pg1 = .false.
    if (present(boost_pg1)) gfr%boost_pg1 = boost_pg1    

    gfr%tolfac = one
    if (par%masterproc) then
       write(iulog, '(a,i3,a,l2,a,i2,a,l2)') 'gfr> init nphys', nphys, ' check', gfr%check, &
            ' ftype', ftype, ' boost_pg1', gfr%boost_pg1
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
         gfr%D_f(nphys,nphys,2,2,nelemd), gfr%Dinv_f(nphys,nphys,2,2,nelemd), &
         gfr%qmin(nlev,max(1,qsize),nelemd), gfr%qmax(nlev,max(1,qsize),nelemd), &
         gfr%phis(nphys2,nelemd), gfr%spherep_f(nphys,nphys,nelemd))
    call gfr_init_fv_metdet(elem, gfr)
    call gfr_init_geometry(elem, gfr)

    if (nphys == 1 .and. gfr%boost_pg1) call gfr_pg1_init(gfr)
  end subroutine gfr_init

  subroutine gfr_finish()
    ! Deallocate the internal gfr structure.

    if (.not. allocated(gfr%fv_metdet)) return
    deallocate(gfr%fv_metdet, gfr%D_f, gfr%Dinv_f, gfr%qmin, gfr%qmax, gfr%phis, &
         gfr%spherep_f)
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
          call gfr_g2f_scalar_and_limit(gfr, ie, elem(ie)%metdet, elem(ie)%state%phis, phis(:,ie))
       end if

       dp = elem(ie)%state%dp3d(:,:,:,nt)
       call gfr_g2f_scalar(ie, elem(ie)%metdet, dp, dp_fv)

       call get_temperature(elem(ie), wr2, hvcoord, nt)
       call get_field(elem(ie), 'p', p, hvcoord, nt, -1)
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

       call get_field(elem(ie), 'omega', wr2, hvcoord, nt, -1)
       call gfr_g2f_scalar(ie, elem(ie)%metdet, wr2, wr1)
#ifdef MODEL_THETA_L
       omega_p(:ncol,:,ie) = reshape(wr1(:nf,:nf,:), (/ncol,nlev/))
#else
       ! for preqx, omega_p = omega/p
       omega_p(:ncol,:,ie) = reshape(wr1(:nf,:nf,:)/p_fv(:nf,:nf,:), (/ncol,nlev/))
#endif
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
    ! Remap T, uv, q states or tendencies from FV to GLL grids.
    !   If ftype is in 1:4, then q is the full mixing ratio state and
    ! dt is not used; if it is not, then q is Qdp tendency, and dt
    ! must be the correct physics time step.

    use element_ops, only: get_field
    use dimensions_mod, only: nlev
    use hybvcoord_mod, only: hvcoord_t
    use physical_constants, only: p0, kappa
    use control_mod, only: ftype

    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: nt
    real(kind=real_kind), intent(in) :: dt
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete
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
       dp = elem(ie)%state%dp3d(:,:,:,nt)
       call gfr_g2f_scalar(ie, elem(ie)%metdet, dp, dp_fv)

       wr1(:nf,:nf,:) = reshape(uv(:ncol,1,:,ie), (/nf,nf,nlev/))
       wr2(:nf,:nf,:) = reshape(uv(:ncol,2,:,ie), (/nf,nf,nlev/))
       call gfr_f2g_vector(gfr, ie, elem, &
            wr1, wr2, elem(ie)%derived%FM(:,:,1,:), elem(ie)%derived%FM(:,:,2,:))

       call get_field(elem(ie), 'p', p, hvcoord, nt, -1)
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
             if (nf > 1 .or. .not. gfr%boost_pg1) then
                ! GLL Q_ten
                call gfr_f2g_scalar_dp(gfr, ie, elem(ie)%metdet, dp_fv, dp, wr1, wr2)
                ! GLL Q1
                elem(ie)%derived%FQ(:,:,:,qi) = elem(ie)%state%Q(:,:,:,qi) + wr2
             else
                ! GLL Q_ten
                do k = 1,nlev
                   elem(ie)%derived%FQ(:,:,k,qi) = wr1(1,1,k)
                end do
             end if
             ! Get limiter bounds.
             do k = 1,nlev
                gfr%qmin(k,qi,ie) = minval(q(:ncol,k,qi,ie))
                gfr%qmax(k,qi,ie) = maxval(q(:ncol,k,qi,ie))
             end do
          else
             ! FV Q_ten
             wr1(:nf,:nf,:) = reshape(q(:ncol,:,qi,ie), (/nf,nf,nlev/))
             wr1(:nf,:nf,:) = dt*wr1(:nf,:nf,:)/dp_fv(:nf,:nf,:)
             if (nf > 1 .or. .not. gfr%boost_pg1) then
                ! GLL Q_ten
                call gfr_f2g_scalar_dp(gfr, ie, elem(ie)%metdet, dp_fv, dp, wr1, wr2)
                ! GLL Q1
                elem(ie)%derived%FQ(:,:,:,qi) = elem(ie)%state%Q(:,:,:,qi) + wr2
             else
                ! GLL Q_ten
                do k = 1,nlev
                   elem(ie)%derived%FQ(:,:,k,qi) = wr1(1,1,k)
                end do
             end if
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

    if (nf == 1 .and. gfr%boost_pg1) return

    do ie = nets,nete
       dp = elem(ie)%state%dp3d(:,:,:,nt)
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
    integer :: ie, ncol
    logical :: augment_in

    augment_in = .false.
    if (present(augment)) augment_in = augment

    ncol = gfr%nphys*gfr%nphys

    do ie = nets,nete
       call gfr_dyn_to_fv_phys_topo_data_elem(ie, elem, square, augment_in, &
            g(npsq*(ie-nets)+1 : npsq*(ie-nets+1)), &
            p(ncol*(ie-nets)+1 : ncol*(ie-nets+1)))
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

    real(kind=real_kind) :: wr(np,np,3), ones(np,np), qmin, qmax, phispg(npsq)
    integer :: nf, ncol, i, j, k

    ones = one
    nf = gfr%nphys
    ncol = nf*nf

    if (augment_variance) then
       ! Estimate additional variance due to remapping from GLL to FV
       ! bases.
       call gfr_g2f_scalar_and_limit(gfr, ie, elem(ie)%metdet, elem(ie)%state%phis, phispg)
       do j = 1,nf
          do i = 1,nf
             ! Integrate (phis_gll - phis_fv)^2 over FV subcell (i,j). Do this
             ! using gfr_g2f_scalar; thus, only one entry out of nf^2 is used.
             k = nf*(j-1) + i
             wr(:,:,2) = ((elem(ie)%state%phis - phispg(k))/grav)**2
             call gfr_g2f_scalar(ie, elem(ie)%metdet, wr(:,:,2:2), wr(:,:,1:1))
             ! Use just entry (i,j).
             wr(i,j,3) = max(zero, wr(i,j,1))
          end do
       end do

       ! Original SGH. augment_variance implies we need to square and sqrt
       ! quantities.
       wr(:,:,1) = reshape(g(:npsq)**2, (/np,np/))
       call gfr_g2f_scalar_and_limit(gfr, ie, elem(ie)%metdet, wr(:,:,1), p(:ncol))

       ! Combine the two sources of variance.
       wr(:nf,:nf,2) = sqrt(reshape(p(:ncol), (/nf,nf/)) + wr(:nf,:nf,3))
       p(:ncol) = reshape(wr(:nf,:nf,2), (/ncol/))
    else
       wr(:,:,1) = reshape(g(:npsq), (/np,np/))
       if (square) wr(:,:,1) = wr(:,:,1)**2
       call gfr_g2f_scalar_and_limit(gfr, ie, elem(ie)%metdet, wr(:,:,1), p(:ncol))
       if (square) p(:ncol) = sqrt(p(:ncol))
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
       if (nf > 1) then
          call gfr_f2g_scalar(ie, elem(ie)%metdet, wr(:,:,1:1), wr(:,:,2:2))
          elem(ie)%state%phis = wr(:,:,2)
       else
          elem(ie)%state%phis = wr(1,1,1)
       end if
    end do

    call gfr_f2g_mixing_ratios_he(hybrid, nets, nete, gfr%qmin(:,:,nets:nete), &
         gfr%qmax(:,:,nets:nete))

    if (nf > 1 .or. .not. gfr%boost_pg1) then
       do ie = nets,nete
          if (gfr%check) wr(:,:,1) = elem(ie)%state%phis
          call limiter_clip_and_sum(np, elem(ie)%spheremp, gfr%qmin(1,1,ie), &
               gfr%qmax(1,1,ie), ones(:,:,1), elem(ie)%state%phis)
          if (gfr%check) then
             if (gfr%qmin(1,1,ie) < zero) then
                write(iulog,*) 'gfr> topo min:', hybrid%par%rank, hybrid%ithr, ie, &
                     gfr%qmin(1,1,ie), 'ERROR'
             end if
             wr(:,:,2) = elem(ie)%state%phis
             call check_f2g_mixing_ratio(gfr, hybrid, ie, 1, elem, gfr%qmin(:1,1,ie), &
                  gfr%qmax(:1,1,ie), ones, wr(:,:,:1), wr(:,:,2:))
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
    real(kind=real_kind), intent(inout) :: ps(:,:), phis(:,:), T(:,:,:), &
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

  subroutine gfr_fv_phys_to_dyn_dom_mt(par, dom_mt, nt, dt, hvcoord, elem, T, uv, q)
    ! Wrapper to the hybrid-threading main routine.

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

    if (.not. par%dynproc) return
#ifdef HORIZ_OPENMP
    !$omp parallel num_threads(hthreads), default(shared), private(nets,nete,hybrid)
#endif
    call gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)
    call gfr_fv_phys_to_dyn_hybrid(hybrid, nt, dt, hvcoord, elem, nets, nete, T, uv, q)
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

    w_ff(:nphys*nphys) = two*two/real(nphys*nphys, real_kind)
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

  subroutine gfr_init_fv_metdet(elem, gfr)
    ! Compute the reference-element-to-sphere Jacobian (metdet) for FV
    ! subcells consistent with those for GLL nodes.

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
    ! Compute various geometric quantities, described below.

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

    ! Jacobian matrices to map a vector between reference element and
    ! sphere.
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

  ! ----------------------------------------------------------------------
  ! Time stepping routines called by the GLL <-> FV remap API routines.

  ! GLL -> FV (g2f)

  subroutine gfr_g2f_remapd(gfr, gll_metdet, fv_metdet, g, f)
    ! Core remap operator. Conservative remap on the reference
    ! element.

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

  subroutine gfr_g2f_scalar(ie, gll_metdet, g, f) ! no gfr b/c public for testing
    ! Wrapper to remapd, where g and f are densities.

    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), g(:,:,:)
    real(kind=real_kind), intent(out) :: f(:,:,:)

    integer :: nlev, k

    nlev = size(g,3)
    do k = 1, nlev
       call gfr_g2f_remapd(gfr, gll_metdet, gfr%fv_metdet(:,ie), g(:,:,k), f(:,:,k))
    end do
  end subroutine gfr_g2f_scalar

  subroutine gfr_g2f_scalar_dp(gfr, ie, gll_metdet, dp_g, dp_f, g, f)
    ! Wrapper to remapd, where g and f are mixing ratios.

    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), dp_g(:,:,:), dp_f(:,:,:), g(:,:,:)
    real(kind=real_kind), intent(out) :: f(:,:,:)

    integer :: nf

    nf = gfr%nphys
    call gfr_g2f_scalar(ie, gll_metdet, dp_g*g, f)
    f(:nf,:nf,:) = f(:nf,:nf,:)/dp_f(:nf,:nf,:)
  end subroutine gfr_g2f_scalar_dp

  subroutine gfr_g2f_vector(gfr, ie, elem, u_g, v_g, u_f, v_f)
    ! Remap a vector on the sphere by doing the actual remap on the
    ! reference element, thus avoiding steep gradients at the poles.

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
    ! Remap dp_g*(u_g,v_g).

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
    ! Remap a mixing ratio conservatively and preventing new extrema.

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
    ! Remap mixing ratios conservatively and preventing new extrema.

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

  subroutine gfr_g2f_scalar_and_limit(gfr, ie, gll_metdet, g, f)
    ! After remap, limit using extremal values from g.

    type (GllFvRemap_t), intent(in) :: gfr
    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), g(:,:)
    real(kind=real_kind), intent(out) :: f(:)

    real(kind=real_kind) :: wr(np,np,2), ones(np,np), qmin, qmax
    integer :: nf, ncol

    ones = one
    nf = gfr%nphys
    ncol = nf*nf

    wr(:np,:np,1) = g
    call gfr_g2f_scalar(ie, gll_metdet, wr(:,:,1:1), wr(:,:,2:2))
    qmin = minval(g(:np,:np))
    qmax = maxval(g(:np,:np))
    wr(:nf,:nf,1) = reshape(gfr%w_ff(:ncol)*gfr%fv_metdet(:ncol,ie), (/nf,nf/))
    call limiter_clip_and_sum(gfr%nphys, wr(:,:,1), qmin, qmax, ones, wr(:nf,:nf,2))
    f(:ncol) = reshape(wr(:nf,:nf,2), (/ncol/))
  end subroutine gfr_g2f_scalar_and_limit

  ! FV -> GLL (f2g)

  subroutine gfr_f2g_scalar(ie, gll_metdet, f, g) ! no gfr b/c public for testing
    ! Wrapper to remapd, where g and f are densities.

    integer, intent(in) :: ie
    real(kind=real_kind), intent(in) :: gll_metdet(:,:), f(:,:,:)
    real(kind=real_kind), intent(out) :: g(:,:,:)

    integer :: k

    do k = 1, size(g,3)
       call gfr_f2g_remapd(gfr, gll_metdet, gfr%fv_metdet(:,ie), f(:,:,k), g(:,:,k))
    end do
  end subroutine gfr_f2g_scalar

  subroutine gfr_f2g_scalar_dp(gfr, ie, gll_metdet, dp_f, dp_g, f, g)
    ! Wrapper to remapd, where g and f are mixing ratios.

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
    ! Remap a vector on the sphere by doing the actual remap on the
    ! reference element, thus avoiding steep gradients at the poles.

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
    ! Remap dp_f*(u_f,v_f).

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
       call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FM, 2*nlev, qsize*nlev, npack)
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
       call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%FM, 2*nlev, qsize*nlev, npack)
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

    call dpotrf('u', n, gfr%pg1sd%Achol, size(gfr%pg1sd%Achol,1), info)
    if (info /= 0) print *, 'gfr ERROR> dpotrf returned', info

    do i = 1,n
       gfr%pg1sd%B(:,i) = Mnp2(i,:)
    end do

    ! Constraint vector c is just w_gg(gfr%pg1sd%inner), so don't store it explicitly.
    gfr%pg1sd%s = reshape(gfr%w_gg(:np,:np), (/np*np/))

    ! Form R's = c
    call dtrtrs('u', 't', 'n', n, 1, gfr%pg1sd%Achol, size(gfr%pg1sd%Achol,1), &
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
    call dtrtrs('u', 't', 'n', n, 1, s%Achol, size(s%Achol,1), x, np*np, info)
    ! Assemble z + (d - s'z)/(s's) s.
    x(:n) = x(:n) + ((mass - sum(s%s(:n)*x(:n)))/s%sts)*s%s(:n)
    ! Solve R x = z + (d - s'z)/(s's) s.
    call dtrtrs('u', 'n', 'n', n, 1, s%Achol, size(s%Achol,1), x, np*np, info)

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
    integer :: ie, nf, ncol

    if (gfr%nphys /= 1 .or. .not. gfr%boost_pg1) return

    ones = one

    do ie = nets,nete
       wr(:,:,1) = elem(ie)%state%phis
       call gfr_pg1_g_reconstruct_scalar(gfr, ie, elem(ie)%metdet, wr(:,:,:1))
       elem(ie)%state%phis = wr(:,:,1)
       call limiter_clip_and_sum(np, elem(ie)%spheremp, gfr%qmin(1,1,ie), &
            gfr%qmax(1,1,ie), ones(:,:,1), elem(ie)%state%phis)
       if (gfr%check) then
          if (gfr%qmin(1,1,ie) < zero) then
             write(iulog,*) 'gfr> topo min:', hybrid%par%rank, hybrid%ithr, ie, &
                  gfr%qmin(1,1,ie), 'ERROR'
          end if
          wr(:,:,2) = elem(ie)%state%phis
          call check_f2g_mixing_ratio(gfr, hybrid, ie, 1, elem, gfr%qmin(:1,1,ie), &
               gfr%qmax(:1,1,ie), ones, wr(:,:,1:1), wr(:,:,2:2))
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

  subroutine gfr_pg1_reconstruct_hybrid(hybrid, nt, dt, hvcoord, elem, nets, nete)
    ! pg1 reconstruction routine for tendencies and states.

    use element_ops, only: get_field
    use dimensions_mod, only: nlev, qsize
    use hybvcoord_mod, only: hvcoord_t
    use physical_constants, only: p0, kappa
    use control_mod, only: ftype

    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: nt
    real(kind=real_kind), intent(in) :: dt
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete

    real(kind=real_kind), dimension(np,np,nlev) :: dp, p, wr1
    real(kind=real_kind) :: qmin, qmax
    integer :: ie, k, qi
    logical :: q_adjustment

    if (gfr%nphys /= 1 .or. .not. gfr%boost_pg1) return

    q_adjustment = ftype >= 1 .and. ftype <= 4

    do ie = nets,nete
       dp = elem(ie)%state%dp3d(:,:,:,nt)

       call get_field(elem(ie), 'p', p, hvcoord, nt, -1)
       wr1 = (p/p0)**kappa
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
          if (gfr%check) wr1 = elem(ie)%derived%FQ(:,:,:,qi)
          do k = 1,nlev
             ! Augment bounds with GLL Q0 bounds.
             gfr%qmin(k,qi,ie) = min(minval(elem(ie)%state%Q(:,:,k,qi)), gfr%qmin(k,qi,ie))
             gfr%qmax(k,qi,ie) = max(maxval(elem(ie)%state%Q(:,:,k,qi)), gfr%qmax(k,qi,ie))
             call limiter_clip_and_sum(np, gfr%w_gg*elem(ie)%metdet, gfr%qmin(k,qi,ie), &
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
    integer :: nlev, k, edgeidx

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

  subroutine gfr_pg1_reconstruct_dom_mt(par, dom_mt, nt, dt, hvcoord, elem)
    ! Wrapper to the hybrid-threading main routine.

    use parallel_mod, only: parallel_t
    use domain_mod, only: domain1d_t
    use thread_mod, only: hthreads
    use hybvcoord_mod, only: hvcoord_t

    type (parallel_t), intent(in) :: par
    type (domain1d_t), intent(in) :: dom_mt(:)
    integer, intent(in) :: nt
    real(kind=real_kind), intent(in) :: dt
    type (hvcoord_t), intent(in) :: hvcoord
    type (element_t), intent(inout) :: elem(:)

    type (hybrid_t) :: hybrid
    integer :: nets, nete
    
    if (.not. par%dynproc) return
#ifdef HORIZ_OPENMP
    !$omp parallel num_threads(hthreads), default(shared), private(nets,nete,hybrid)
#endif
    call gfr_hybrid_create(par, dom_mt, hybrid, nets, nete)
    call gfr_pg1_reconstruct_hybrid(hybrid, nt, dt, hvcoord, elem, nets, nete)
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

  subroutine gfr_f_get_latlon(ie, i, j, lat, lon)
    ! Get (lat,lon) of FV point i,j.

    integer, intent(in) :: ie, i, j
    real(kind=real_kind), intent(out) :: lat, lon

    lat = gfr%spherep_f(i,j,ie)%lat
    lon = gfr%spherep_f(i,j,ie)%lon
  end subroutine gfr_f_get_latlon

  subroutine gfr_f_get_cartesian3d(ie, i, j, p)
    ! Get (x,y,z) of FV point i,j.

    use coordinate_systems_mod, only: cartesian3D_t, change_coordinates

    integer, intent(in) :: ie, i, j
    type (cartesian3D_t), intent(out) :: p

    p = change_coordinates(gfr%spherep_f(i,j,ie))
  end subroutine gfr_f_get_cartesian3d

  subroutine limiter_clip_and_sum(n, spheremp, qmin, qmax, dp, q)
    ! CAAS as described in Alg 3.1 of doi:10.1137/18M1165414. q is a
    ! mixing ratio. Solve
    !    min_q* norm(dp q - dp q*, 1)
    !     st    spheremp'(dp q*) = spheremp'(dp q)
    !           qmin < q* < qmax

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

  subroutine check_g2f_mixing_ratio(gfr, hybrid, ie, qi, elem, dp, dp_fv, q_g, q_f)
    ! Check that gfr_g2f_mixing_ratio found a property-preserving
    ! solution.

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
       den = gfr%tolfac*max(1e-10_real_kind, maxval(abs(elem(ie)%state%Q(:,:,k,qi))))
       mass_f = sum(reshape(gfr%w_ff(:nf2)*gfr%fv_metdet(:nf2,ie), (/nf,nf/))* &
            dp_fv(:nf,:nf,k)*q_f(:nf,:nf,k))
       mass_g = sum(elem(ie)%spheremp*dp(:,:,k)*q_g(:,:,k))
       if (qmin_f < qmin_g - 10*eps*den .or. qmax_f > qmax_g + 10*eps*den) then
          write(iulog,*) 'gfr> g2f mixing ratio limits:', hybrid%par%rank, hybrid%ithr, ie, qi, k, &
               qmin_g, qmin_f-qmin_g, qmax_f-qmax_g, qmax_g, mass_f, mass_g, 'ERROR'
       end if
       if (abs(mass_f - mass_g) > gfr%tolfac*20*eps*max(mass_f, mass_g)) then
          write(iulog,*) 'gfr> g2f mixing ratio mass:', hybrid%par%rank, hybrid%ithr, ie, qi, k, &
               qmin_g, qmax_g, mass_f, mass_g, 'ERROR'
       end if
    end do
  end subroutine check_g2f_mixing_ratio

  subroutine check_f2g_mixing_ratio(gfr, hybrid, ie, qi, elem, qmin, qmax, dp, q0_g, q1_g)
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
    integer :: q, k

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
       end if
       den = sum(elem(ie)%spheremp*dp(:,:,k)*maxval(abs(q0_g(:,:,k))))
       if (abs(mass1 - mass0) > gfr%tolfac*20*eps*den) then
          write(iulog,*) 'gfr> f2g mixing ratio mass:', hybrid%par%rank, hybrid%ithr, ie, qi, k, &
               qmin_f, qmin_g, qmax_g, qmax_f, mass0, mass1, 'ERROR'
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

    real(kind=real_kind) :: a, b, rd, x, y, f0(np,np), f1(np,np), g(np,np), &
         wrk(np,np), qmin, qmax, qmin1, qmax1, wr1(np,np,1)
    integer :: nf, nf2, ie, i, j, iremap, info, ilimit, it
    real(kind=real_kind), allocatable :: Qdp_fv(:,:,:), ps_v_fv(:,:,:), &
         qmins(:,:,:), qmaxs(:,:,:)
    logical :: limit
    character(32) :: msg

    ! Purposely construct our own hybrid object to test gfr_hybrid_create.
    type (hybrid_t) :: hybrid
    integer :: nets, nete

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
             call gfr_g2f_remapd(gfr, elem(ie)%metdet, gfr%fv_metdet(:,ie), &
                  elem(ie)%state%ps_v(:,:,1), ps_v_fv(:,:,ie))
             if (limit) then
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
             ! 2c. Augment bounds with current values.
             do ie = nets, nete
                wrk = elem(ie)%state%Q(:,:,1,1)
                qmins(1,1,ie) = min(qmins(1,1,ie), minval(wrk))
                qmaxs(1,1,ie) = max(qmaxs(1,1,ie), maxval(wrk))                
             end do
          endif
          ! 2d. Remap
          if (nf == 1 .and. gfr%boost_pg1) then
             do ie = nets, nete
                elem(ie)%state%Q(:,:,1,1) = Qdp_fv(1,1,ie)/ps_v_fv(1,1,ie)
             end do
          else
             do ie = nets, nete
                call gfr_f2g_remapd(gfr, elem(ie)%metdet, gfr%fv_metdet(:,ie), &
                     Qdp_fv(:,:,ie), elem(ie)%state%Q(:,:,1,1))
                elem(ie)%state%Q(:,:,1,1) = elem(ie)%state%Q(:,:,1,1)/elem(ie)%state%ps_v(:,:,1)
                if (limit) then
                   call limiter_clip_and_sum(np, elem(ie)%spheremp, & ! same as w_gg*gll_metdet
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
                      call limiter_clip_and_sum(np, gfr%w_gg*elem(ie)%metdet, qmins(1,1,ie), &
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
    ! Driver for check subroutine.

    use domain_mod, only: domain1d_t
    use derivative_mod, only: derivative_t
    use hybvcoord_mod, only: hvcoord_t

    type (hybrid_t), intent(in) :: hybrid
    type (domain1d_t), intent(in) :: dom_mt(:)
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(inout) :: elem(:)
    type (hvcoord_t) , intent(in) :: hvcoord

    integer :: nphys, bi
    logical :: boost_pg1

    do nphys = 1, np
       do bi = 1,2
          if (nphys > 1 .and. bi > 1) exit
          boost_pg1 = bi == 2

          ! This is meant to be called before threading starts.
          if (hybrid%ithr == 0) call gfr_init(hybrid%par, elem, nphys, .true., boost_pg1)
          !$omp barrier

          call check(hybrid%par, dom_mt, gfr, elem, .false.)

          ! This is meant to be called after threading ends.
          !$omp barrier
          if (hybrid%ithr == 0) call gfr_finish()
          !$omp barrier
       end do
    end do
  end subroutine gfr_test
end module gllfvremap_mod
