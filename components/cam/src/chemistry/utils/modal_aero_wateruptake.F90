module modal_aero_wateruptake

!   RCE 07.04.13:  Adapted from MIRAGE2 code

use shr_kind_mod,     only: r8 => shr_kind_r8
use physconst,        only: pi, rhoh2o, r_universal
use ppgrid,           only: pcols, pver
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field

use wv_saturation,    only: qsat_water
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_props, rad_cnst_get_mode_num
use cam_history,      only: addfld, add_default, outfld
use cam_logfile,      only: iulog
use ref_pres,         only: top_lev => clim_modal_aero_top_lev
use phys_control,     only: phys_getopts
use cam_abortutils,       only: endrun

implicit none
private
save

public :: &
   modal_aero_wateruptake_init, &
   modal_aero_wateruptake_dr

public :: modal_aero_wateruptake_reg

real(r8), parameter :: third = 1._r8/3._r8
real(r8), parameter :: pi43  = pi*4.0_r8/3.0_r8
real(r8), parameter :: mwh20 = 18.0_r8 ! molecular weight of water (g/mole)
real(r8), parameter :: surften = 76._r8 ! surface tension of water/air interface at 273K (mN/m)
integer, parameter :: max_iterations = 200


! Physics buffer indices
integer :: cld_idx        = 0
integer :: dgnum_idx      = 0
integer :: dgnumwet_idx   = 0
integer :: wetdens_ap_idx = 0
integer :: qaerwat_idx    = 0

logical :: pergro_mods         = .false.
!===============================================================================
contains
!===============================================================================

subroutine modal_aero_wateruptake_reg()

  use physics_buffer,   only: pbuf_add_field, dtype_r8
  use rad_constituents, only: rad_cnst_get_info

   integer :: nmodes

   call rad_cnst_get_info(0, nmodes=nmodes)
   call pbuf_add_field('DGNUMWET',   'global',  dtype_r8, (/pcols, pver, nmodes/), dgnumwet_idx)
   call pbuf_add_field('WETDENS_AP', 'physpkg', dtype_r8, (/pcols, pver, nmodes/), wetdens_ap_idx)

   ! 1st order rate for direct conversion of strat. cloud water to precip (1/s)
   call pbuf_add_field('QAERWAT',    'physpkg', dtype_r8, (/pcols, pver, nmodes/), qaerwat_idx)

end subroutine modal_aero_wateruptake_reg

!===============================================================================
!===============================================================================

subroutine modal_aero_wateruptake_init(pbuf2d)
   use time_manager,  only: is_first_step
   use physics_buffer,only: pbuf_set_field

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   integer :: m, nmodes
   logical :: history_aerosol      ! Output the MAM aerosol variables and tendencies
   logical :: history_verbose      ! produce verbose history output

   character(len=3) :: trnum       ! used to hold mode number (as characters)
   !----------------------------------------------------------------------------

   cld_idx        = pbuf_get_index('CLD')
   dgnum_idx      = pbuf_get_index('DGNUM')

   ! assume for now that will compute wateruptake for climate list modes only

   call rad_cnst_get_info(0, nmodes=nmodes)

   ! determine default variables
   call phys_getopts(history_aerosol_out = history_aerosol, &
                     history_verbose_out = history_verbose, &
                     pergro_mods_out = pergro_mods)

   do m = 1, nmodes
      write(trnum, '(i3.3)') m
      call addfld('dgnd_a'//trnum(2:3), (/ 'lev' /), 'A', 'm', &
         'dry dgnum, interstitial, mode '//trnum(2:3))
      call addfld('dgnw_a'//trnum(2:3), (/ 'lev' /), 'A', 'm', &
         'wet dgnum, interstitial, mode '//trnum(2:3))
      call addfld('wat_a'//trnum(3:3), (/ 'lev' /), 'A', 'm', &
         'aerosol water, interstitial, mode '//trnum(2:3))

      if (history_aerosol) then
         if (history_verbose) then
            call add_default('dgnd_a'//trnum(2:3), 1, ' ')
            call add_default('dgnw_a'//trnum(2:3), 1, ' ')
            call add_default('wat_a'//trnum(3:3),  1, ' ')
         endif
      endif

   end do

   ! Add total aerosol water
   if (history_aerosol .and. .not. history_verbose) then
      call addfld('aero_water', (/ 'lev' /), 'A', 'm', &
         'sum of aerosol water of interstitial modes wat_a1+wat_a2+wat_a3+wat_a4' )
      call add_default( 'aero_water',  1, ' ')
   endif

   if (is_first_step()) then
      ! initialize fields in physics buffer
      call pbuf_set_field(pbuf2d, dgnumwet_idx, 0.0_r8)
   endif

end subroutine modal_aero_wateruptake_init

!===============================================================================


subroutine modal_aero_wateruptake_dr(state, pbuf, list_idx_in, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m)
!-----------------------------------------------------------------------
!
! CAM specific driver for modal aerosol water uptake code.
!
! *** N.B. *** The calculation has been enabled for diagnostic mode lists
!              via optional arguments.  If the list_idx arg is present then
!              all the optional args must be present.
!
!-----------------------------------------------------------------------

   ! Arguments
   type(physics_state), target, intent(in)    :: state          ! Physics state variables
   type(physics_buffer_desc),   pointer       :: pbuf(:)        ! physics buffer

   integer,  optional,          intent(in)    :: list_idx_in
   real(r8), optional,          pointer       :: dgnumdry_m(:,:,:)
   real(r8), optional,          pointer       :: dgnumwet_m(:,:,:)
   real(r8), optional,          pointer       :: qaerwat_m(:,:,:)
   real(r8), optional,          pointer       :: wetdens_m(:,:,:)

   ! local variables

   integer  :: lchnk              ! chunk index
   integer  :: ncol               ! number of columns
   integer  :: list_idx           ! radiative constituents list index
   integer  :: stat

   integer :: i, k, l, m
   integer :: itim_old
   integer :: nmodes
   integer :: nspec

   real(r8), pointer :: h2ommr(:,:) ! specific humidity
   real(r8), pointer :: t(:,:)      ! temperatures (K)
   real(r8), pointer :: pmid(:,:)   ! layer pressure (Pa)
   real(r8), pointer :: raer(:,:)   ! aerosol species MRs (kg/kg and #/kg)

   real(r8), pointer :: cldn(:,:)      ! layer cloud fraction (0-1)
   real(r8), pointer :: dgncur_a(:,:,:)
   real(r8), pointer :: dgncur_awet(:,:,:)
   real(r8), pointer :: wetdens(:,:,:)
   real(r8), pointer :: qaerwat(:,:,:)

   real(r8), allocatable :: maer(:,:,:)      ! aerosol wet mass MR (including water) (kg/kg-air)
   real(r8), allocatable :: hygro(:,:,:)     ! volume-weighted mean hygroscopicity (--)
   real(r8), allocatable :: naer(:,:,:)      ! aerosol number MR (bounded!) (#/kg-air)
   real(r8), allocatable :: dryvol(:,:,:)    ! single-particle-mean dry volume (m3)
   real(r8), allocatable :: drymass(:,:,:)   ! single-particle-mean dry mass  (kg)
   real(r8), allocatable :: dryrad(:,:,:)    ! dry volume mean radius of aerosol (m)

   real(r8), allocatable :: wetrad(:,:,:)    ! wet radius of aerosol (m)
   real(r8), allocatable :: wetvol(:,:,:)    ! single-particle-mean wet volume (m3)
   real(r8), allocatable :: wtrvol(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)

   real(r8), allocatable :: rhcrystal(:)
   real(r8), allocatable :: rhdeliques(:)
   real(r8), allocatable :: specdens_1(:)

   real(r8) :: dryvolmr(pcols,pver)          ! volume MR for aerosol mode (m3/kg)
   real(r8) :: specdens
   real(r8) :: spechygro, spechygro_1
   real(r8) :: duma, dumb
   real(r8) :: sigmag
   real(r8) :: alnsg
   real(r8) :: v2ncur_a
   real(r8) :: drydens               ! dry particle density  (kg/m^3)
   real(r8) :: rh(pcols,pver)        ! relative humidity (0-1)

   real(r8) :: es(pcols)             ! saturation vapor pressure
   real(r8) :: qs(pcols)             ! saturation specific humidity
   real(r8) :: cldn_thresh
   real(r8) :: aerosol_water(pcols,pver) !sum of aerosol water (wat_a1 + wat_a2 + wat_a3 + wat_a4)
   logical :: history_aerosol      ! Output the MAM aerosol variables and tendencies
   logical :: history_verbose      ! produce verbose history output

   character(len=3) :: trnum       ! used to hold mode number (as characters)
   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol = state%ncol

   ! determine default variables
   call phys_getopts(history_aerosol_out = history_aerosol, &
                     history_verbose_out = history_verbose)

   list_idx = 0
   if (present(list_idx_in)) then
      list_idx = list_idx_in

      ! check that all optional args are present
      if (.not. present(dgnumdry_m) .or. .not. present(dgnumwet_m) .or. &
          .not. present(qaerwat_m)  .or. .not. present(wetdens_m)) then
         call endrun('modal_aero_wateruptake_dr called for'// &
                     'diagnostic list but required args not present')
      end if

      ! arrays for diagnostic calculations must be associated
      if (.not. associated(dgnumdry_m) .or. .not. associated(dgnumwet_m) .or. &
          .not. associated(qaerwat_m)  .or. .not. associated(wetdens_m)) then
         call endrun('modal_aero_wateruptake_dr called for'// &
                     'diagnostic list but required args not associated')
      end if
   end if

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   allocate( &
      maer(pcols,pver,nmodes),     &
      hygro(pcols,pver,nmodes),    &
      naer(pcols,pver,nmodes),     &
      dryvol(pcols,pver,nmodes),   &
      drymass(pcols,pver,nmodes),  &
      dryrad(pcols,pver,nmodes),   &
      wetrad(pcols,pver,nmodes),   &
      wetvol(pcols,pver,nmodes),   &
      wtrvol(pcols,pver,nmodes),   &
      rhcrystal(nmodes),           &
      rhdeliques(nmodes),          &
      specdens_1(nmodes)           )

   maer(:,:,:)     = 0._r8
   hygro(:,:,:)    = 0._r8


   if (list_idx == 0) then
      call pbuf_get_field(pbuf, dgnum_idx,      dgncur_a )
      call pbuf_get_field(pbuf, dgnumwet_idx,   dgncur_awet )
      call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens)
      call pbuf_get_field(pbuf, qaerwat_idx,    qaerwat)
   else
      dgncur_a    => dgnumdry_m
      dgncur_awet => dgnumwet_m
      qaerwat     => qaerwat_m
      wetdens     => wetdens_m
   end if

   do m = 1, nmodes

      dryvolmr(:,:) = 0._r8

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigmag,  &
         rhcrystal=rhcrystal(m), rhdeliques=rhdeliques(m))

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      do l = 1, nspec

         ! get species interstitial mixing ratio ('a')
         call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, raer)
         call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, hygro_aer=spechygro)

         if (l == 1) then
            ! save off these values to be used as defaults
            specdens_1(m)  = specdens
            spechygro_1    = spechygro
         end if

         do k = top_lev, pver
            do i = 1, ncol
               duma          = raer(i,k)
               maer(i,k,m)   = maer(i,k,m) + duma
               dumb          = duma/specdens
               dryvolmr(i,k) = dryvolmr(i,k) + dumb
               hygro(i,k,m)  = hygro(i,k,m) + dumb*spechygro
            end do
         end do
      end do

      alnsg = log(sigmag)

      do k = top_lev, pver
         do i = 1, ncol

            if (dryvolmr(i,k) > 1.0e-30_r8) then
               hygro(i,k,m) = hygro(i,k,m)/dryvolmr(i,k)
            else
               hygro(i,k,m) = spechygro_1
            end if

            ! dry aerosol properties

            v2ncur_a = 1._r8 / ( (pi/6._r8)*(dgncur_a(i,k,m)**3._r8)*exp(4.5_r8*alnsg**2._r8) )
            ! naer = aerosol number (#/kg)
            naer(i,k,m) = dryvolmr(i,k)*v2ncur_a

            ! compute mean (1 particle) dry volume and mass for each mode
            ! old coding is replaced because the new (1/v2ncur_a) is equal to
            ! the mean particle volume
            ! also moletomass forces maer >= 1.0e-30, so (maer/dryvolmr)
            ! should never cause problems (but check for maer < 1.0e-31 anyway)
            if (maer(i,k,m) .gt. 1.0e-31_r8) then
               drydens = maer(i,k,m)/dryvolmr(i,k)
            else
               drydens = 1.0_r8
            end if
            dryvol(i,k,m)   = 1.0_r8/v2ncur_a
            drymass(i,k,m)  = drydens*dryvol(i,k,m)
            dryrad(i,k,m)   = (dryvol(i,k,m)/pi43)**third

         end do
      end do

   end do    ! modes

   ! relative humidity calc

   h2ommr => state%q(:,:,1)
   t      => state%t
   pmid   => state%pmid

   itim_old    =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx, cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   do k = top_lev, pver
      call qsat_water(t(:ncol,k), pmid(:ncol,k), es(:ncol), qs(:ncol))
      do i = 1, ncol
         if (qs(i) > h2ommr(i,k)) then
            rh(i,k) = h2ommr(i,k)/qs(i)
         else
            rh(i,k) = 0.98_r8
         endif
         rh(i,k) = max(rh(i,k), 0.0_r8)
         rh(i,k) = min(rh(i,k), 0.98_r8)
         if(pergro_mods) then
            cldn_thresh = 0.9998_r8
         else
            cldn_thresh = 1.0_r8 !original code
         endif
         if (cldn(i,k) .lt. cldn_thresh) then !BSINGH
            rh(i,k) = (rh(i,k) - cldn(i,k)) / (1.0_r8 - cldn(i,k))  ! clear portion
         end if
         rh(i,k) = max(rh(i,k), 0.0_r8)
      end do
   end do

   call modal_aero_wateruptake_sub( &
      ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
      hygro, rh, dryvol, wetrad, wetvol,           &
      wtrvol)

   do m = 1, nmodes

      do k = top_lev, pver
         do i = 1, ncol

            dgncur_awet(i,k,m) = dgncur_a(i,k,m) * (wetrad(i,k,m)/dryrad(i,k,m))
            qaerwat(i,k,m)     = rhoh2o*naer(i,k,m)*wtrvol(i,k,m)

            ! compute aerosol wet density (kg/m3)
            if (wetvol(i,k,m) > 1.0e-30_r8) then
               wetdens(i,k,m) = (drymass(i,k,m) + rhoh2o*wtrvol(i,k,m))/wetvol(i,k,m)
            else
               wetdens(i,k,m) = specdens_1(m)
            end if
         end do
      end do

   end do    ! modes

   if (list_idx == 0) then

      aerosol_water(:ncol,:) = 0._r8
      do m = 1, nmodes
         ! output to history
         write( trnum, '(i3.3)' ) m
         call outfld( 'wat_a'//trnum(3:3),  qaerwat(:,:,m),     pcols, lchnk)
         call outfld( 'dgnd_a'//trnum(2:3), dgncur_a(:,:,m),    pcols, lchnk)
         call outfld( 'dgnw_a'//trnum(2:3), dgncur_awet(:,:,m), pcols, lchnk)
         if (history_aerosol .and. .not. history_verbose) &
         aerosol_water(:ncol,:) = aerosol_water(:ncol,:) + qaerwat(:ncol,:,m)
      end do

      if (history_aerosol .and. .not. history_verbose) &
         call outfld( 'aero_water',  aerosol_water(:ncol,:),    ncol, lchnk)

   end if

   deallocate( &
      maer,   hygro,  naer,   dryvol,    drymass,    dryrad, &
      wetrad, wetvol, wtrvol, rhcrystal, rhdeliques, specdens_1)

end subroutine modal_aero_wateruptake_dr

!===============================================================================

subroutine modal_aero_wateruptake_sub( &
   ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
   hygro, rh, dryvol, wetrad, wetvol,           &
   wtrvol)

!-----------------------------------------------------------------------
!
! Purpose: Compute aerosol wet radius
!
! Method:  Kohler theory
!
! Author:  S. Ghan
!
!-----------------------------------------------------------------------

   ! Arguments
   integer, intent(in)  :: ncol                    ! number of columns
   integer, intent(in)  :: nmodes

   real(r8), intent(in) :: rhcrystal(:)
   real(r8), intent(in) :: rhdeliques(:)
   real(r8), intent(in) :: dryrad(:,:,:)         ! dry volume mean radius of aerosol (m)
   real(r8), intent(in) :: hygro(:,:,:)          ! volume-weighted mean hygroscopicity (--)
   real(r8), intent(in) :: rh(:,:)               ! relative humidity (0-1)
   real(r8), intent(in) :: dryvol(:,:,:)

   real(r8), intent(out) :: wetrad(:,:,:)        ! wet radius of aerosol (m)
   real(r8), intent(out) :: wetvol(:,:,:)        ! single-particle-mean wet volume (m3)
   real(r8), intent(out) :: wtrvol(:,:,:)        ! single-particle-mean water volume in wet aerosol (m3)

   ! local variables

   integer :: i, k, m

   real(r8) :: hystfac                ! working variable for hysteresis
   !-----------------------------------------------------------------------


   ! loop over all aerosol modes
   do m = 1, nmodes

      hystfac = 1.0_r8 / max(1.0e-5_r8, (rhdeliques(m) - rhcrystal(m)))

      do k = top_lev, pver
         do i = 1, ncol

            ! compute wet radius for each mode
            call modal_aero_kohler(dryrad(i:i,k,m), hygro(i:i,k,m), rh(i:i,k), wetrad(i:i,k,m), 1)

            wetrad(i,k,m) = max(wetrad(i,k,m), dryrad(i,k,m))
            wetvol(i,k,m) = pi43*wetrad(i,k,m)**3
            wetvol(i,k,m) = max(wetvol(i,k,m), dryvol(i,k,m))
            wtrvol(i,k,m) = wetvol(i,k,m) - dryvol(i,k,m)
            wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0_r8)

            ! apply simple treatment of deliquesence/crystallization hysteresis
            ! for rhcrystal < rh < rhdeliques, aerosol water is a fraction of
            ! the "upper curve" value, and the fraction is a linear function of rh
            if (rh(i,k) < rhcrystal(m)) then
               wetrad(i,k,m) = dryrad(i,k,m)
               wetvol(i,k,m) = dryvol(i,k,m)
               wtrvol(i,k,m) = 0.0_r8
            else if (rh(i,k) < rhdeliques(m)) then
               wtrvol(i,k,m) = wtrvol(i,k,m)*hystfac*(rh(i,k) - rhcrystal(m))
               wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0_r8)
               wetvol(i,k,m) = dryvol(i,k,m) + wtrvol(i,k,m)
               wetrad(i,k,m) = (wetvol(i,k,m)/pi43)**third
            end if

         end do  ! columns
      end do     ! levels

   end do ! modes

end subroutine modal_aero_wateruptake_sub

! Computes the value of Kohler polynomial (Pruppacher & Klett 2nd ed. eqn 6-37)
pure function kohler_polynomial(rwet, rdry, hyg, s, a)
  real(r8) :: kohler_polynomial
  real(r8), intent(in) :: rwet ! wet radius (microns)
  real(r8), intent(in) :: rdry ! dry radius (microns)
  real(r8), intent(in) :: hyg ! hygroscopicity
  real(r8), intent(in) :: s ! relative humidity
  real(r8), intent(in) :: a ! constant kelvin parameter for surface tension effects
  !
  real(r8) :: rdcube, rwcube, logrh

  rdcube = rdry**3
  rwcube = rwet**3
  logrh = log(s)

  kohler_polynomial = logrh*rwet*rwcube - a*rwcube + rdcube*(hyg-logrh)*rwet + a*rdcube
end function

! Computes the derivative of Kohler polynomial (Pruppacher & Klett 2nd ed. eqn 6-37) w.r.t. wet radius
pure function kohler_derivative(rwet, rdry, hyg, s, a)
  real(r8) :: kohler_derivative
  real(r8), intent(in) :: rwet ! wet radius (microns)
  real(r8), intent(in) :: rdry ! dry radius (microns)
  real(r8), intent(in) :: hyg ! hygroscopicity
  real(r8), intent(in) :: s ! relative humidity
  real(r8), intent(in) :: a ! constant kelvin parameter for surface tension effects
  !
  real(r8) :: rwsq, rdcube, logrh

  rdcube = rdry**3
  rwsq = rwet*rwet
  logrh = log(s)

  kohler_derivative = 4.0_r8*logrh*rwet*rwsq - 3.0_r8*a*rwsq + rdcube*(hyg-logrh)
end function

subroutine newton_solve_kohler(rwet, rdry, hyg, s, a, tol, niter)
  real(r8), intent(inout) :: rwet ! wet radius (microns)
  real(r8), intent(in) :: rdry ! dry radius (microns)
  real(r8), intent(in) :: hyg ! hygroscopicity
  real(r8), intent(in) :: s ! relative humidity
  real(r8), intent(in) :: a ! constant kelvin parameter for surface tension effects
  real(r8), intent(in) :: tol ! convergence tolerance
  integer, intent(out) :: niter ! iteration counter
  !
  real(r8) :: dr, k, kprime, new_rwet

  dr = huge(1.0_r8)
  niter = 0
  rwet = 700
  do while (dr > tol .and. niter < max_iterations)
    niter = niter + 1
    k = kohler_polynomial(rwet, rdry, hyg, s, a)
    kprime = kohler_derivative(rwet, rdry, hyg, s, a)
    new_rwet = rwet - k/kprime
    dr = abs(new_rwet - rwet)
    rwet = new_rwet
  enddo
end subroutine

subroutine bisection_solve_kohler(rwet, rdry, hyg, s, a, tol, niter)
  real(r8), intent(inout) :: rwet ! wet radius (microns)
  real(r8), intent(in) :: rdry ! dry radius (microns)
  real(r8), intent(in) :: hyg ! hygroscopicity
  real(r8), intent(in) :: s ! relative humidity
  real(r8), intent(in) :: a ! constant kelvin parameter for surface tension effects
  real(r8), intent(in) :: tol ! convergence tolerance
  integer, intent(out) :: niter ! iteration counter
  !
  real(r8) :: left, right
  real(r8) :: kleft, kmid

  left = 0.9_r8*rdry
  right = 750.0_r8
  niter = 0
  kleft = kohler_polynomial(left, rdry, hyg, s, a)
  rwet = 0.5_r8*(left + right)
  do while (right-left > tol .and. niter < max_iterations)
    niter = niter + 1
    kmid = kohler_polynomial(rwet, rdry, hyg, s, a)
    if (kmid*kleft < 0) then
      left = left
      right = rwet
      kleft = kleft
    else
      left = rwet
      right = right
      kleft = kmid
    endif
    rwet = 0.5*(left + right)
  enddo
end subroutine

!-----------------------------------------------------------------------
      subroutine modal_aero_kohler(   &
          rdry_in, hygro, s, rwet_out, im )

! calculates equlibrium radius r of haze droplets as function of
! dry particle mass and relative humidity s using kohler solution
! given in pruppacher and klett (eqn 6-35) or (6-37, 2nd ed.)

! for multiple aerosol types, assumes an internal mixture of aerosols

      implicit none

! arguments
      integer :: im         ! number of grid points to be processed
      real(r8) :: rdry_in(:)    ! aerosol dry radius (m)
      real(r8) :: hygro(:)      ! aerosol volume-mean hygroscopicity (--)
      real(r8) :: s(:)          ! relative humidity (1 = saturated)
      real(r8) :: rwet_out(:)   ! aerosol wet radius (m)

! local variables
      integer, parameter :: imax=200
      integer :: i
      real(r8) :: r(im)         ! wet radius (microns)
      real(r8) :: rdry(im)      ! dry radius (microns)
      real(r8), parameter :: ugascon = r_universal * 1.0e-3_r8 ! convert gas constant to per mole (not per kmole) units
      real(r8), parameter :: tair = 273.15K ! temperature of air / water interface for surface tension
      real(r8), parameter :: a = 2.0_r8 * mwh20 * surften /(ugascon * tair * rhoh2o)

      real(r8), parameter :: convergence_tol = 1.0e-5_r8
      integer :: ctr

      do i=1,im
        rdry(i) = rdry_in(i)*1.0e6_r8   ! convert (m) to (microns)
        call newton_solve_kohler(r(i), rdry(i), hygro(i), s(i), a, convergence_tol, ctr)
        if (r(i) < 0) then
          ! Newton's method found the wrong root.  Revert to bisection method (guaranteed to converge)
          call bisection_solve_kohler(r(i), rdry(i), hygro(i), s(i), a, convergence_tol, ctr)
          write(iulog,*) "modal_aero_kohler: Newton's method found the wrong root."
          write(iulog,*) " reverting to bisection method (guaranteed to converge)."
        endif
        if (ctr >= max_iterations) then
          write(iulog,*) "modal_aero_kohler: no solution found (failed to converge)."
          write(iulog,*) "(hyg, rel_humid) = ", [hygro(i), s(i)]
          write(iulog,*) " setting radius to dry radius=", rdry(i)
          r(i) = rdry(i)
        endif
      end do

! bound and convert from microns to m
      do i=1,im
         r(i) = min(r(i),30._r8) ! upper bound based on 1 day lifetime
         rwet_out(i) = r(i)*1.e-6_r8
      end do

      return
      end subroutine modal_aero_kohler



   end module modal_aero_wateruptake


