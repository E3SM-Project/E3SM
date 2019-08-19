module modal_aero_wateruptake

!   RCE 07.04.13:  Adapted from MIRAGE2 code

use shr_kind_mod,     only: r8 => shr_kind_r8
use physconst,        only: pi, rhoh2o
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
   modal_aero_wateruptake_dr,   &
   modal_aero_wateruptake_sub

public :: modal_aero_wateruptake_reg

! JS added on 08-05-2019
public :: aw_scheme, aw_rh_max, aw_mam_rh_min, aw_m7_rh_min, &
          aw_l_hysteresis, aw_l_rh_clearsky, aw_l_m7_rh_min

real(r8), parameter :: third = 1._r8/3._r8
real(r8), parameter :: pi43  = pi*4.0_r8/3.0_r8

! JS added on 08-05-2019: Initialize the new added variables
integer  :: aw_scheme        = 0          ! 0 - default MAM scheme; 1 - k-kohler ; 2 - M7 (Jacobson 1996); 3 - M7 (Jacobson 2005)
real(r8) :: aw_mam_rh_min    = 0.35_r8    ! the minimum rh below which wet radius is equal to dry radius, used for MAM
real(r8) :: aw_m7_rh_min     = 0.45_r8    ! the minimum rh below which wet radius is equal to dry radius, used for M7
real(r8) :: aw_rh_max        = 0.98_r8    ! the maximum rh to solve for water uptake process
logical  :: aw_l_hysteresis  = .TRUE.     ! TRUE if hysteresis effect is considered for MAM
logical  :: aw_l_rh_clearsky = .TRUE.     ! TRUE if clear-sky RH is used for water uptake process
logical  :: aw_l_m7_rh_min   = .FALSE.    ! FALSE if RHmin and highest molality are not considered

! Physics buffer indices
integer :: cld_idx        = 0
integer :: dgnum_idx      = 0
integer :: dgnumwet_idx   = 0
integer :: wetdens_ap_idx = 0
integer :: qaerwat_idx    = 0

!!== KZ_INSITU 
integer :: dgnumwet_rh40_idx = 0
integer :: dgnumwet_rh55_idx = 0
integer :: dgnumwet_rh65_idx = 0
integer :: dgnumwet_rh75_idx = 0
integer :: dgnumwet_rh85_idx = 0
integer :: qaerwat_rh40_idx  = 0
integer :: qaerwat_rh55_idx  = 0
integer :: qaerwat_rh65_idx  = 0
integer :: qaerwat_rh75_idx  = 0
integer :: qaerwat_rh85_idx  = 0
!!== KZ_INSITU 

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
!!== KZ_INSITU 
   call pbuf_add_field('DGNUMWET40', 'global',  dtype_r8, (/pcols, pver, nmodes/), dgnumwet_rh40_idx)
   call pbuf_add_field('DGNUMWET55', 'global',  dtype_r8, (/pcols, pver, nmodes/), dgnumwet_rh55_idx)
   call pbuf_add_field('DGNUMWET65', 'global',  dtype_r8, (/pcols, pver, nmodes/), dgnumwet_rh65_idx)
   call pbuf_add_field('DGNUMWET75', 'global',  dtype_r8, (/pcols, pver, nmodes/), dgnumwet_rh75_idx)
   call pbuf_add_field('DGNUMWET85', 'global',  dtype_r8, (/pcols, pver, nmodes/), dgnumwet_rh85_idx)
   call pbuf_add_field('QAERWAT40',  'physpkg', dtype_r8, (/pcols, pver, nmodes/), qaerwat_rh40_idx)
   call pbuf_add_field('QAERWAT55',  'physpkg', dtype_r8, (/pcols, pver, nmodes/), qaerwat_rh55_idx)
   call pbuf_add_field('QAERWAT65',  'physpkg', dtype_r8, (/pcols, pver, nmodes/), qaerwat_rh65_idx)
   call pbuf_add_field('QAERWAT75',  'physpkg', dtype_r8, (/pcols, pver, nmodes/), qaerwat_rh75_idx)
   call pbuf_add_field('QAERWAT85',  'physpkg', dtype_r8, (/pcols, pver, nmodes/), qaerwat_rh85_idx)
!!== KZ_INSITU 

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

   use modal_aero_data, only: lspectype_amode, specname_amode

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

   integer :: i, j, k, l, m
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
!!== KZ_INSITU
   real(r8), pointer :: dgncur_awet40(:,:,:)
   real(r8), pointer :: dgncur_awet55(:,:,:)
   real(r8), pointer :: dgncur_awet65(:,:,:)
   real(r8), pointer :: dgncur_awet75(:,:,:)
   real(r8), pointer :: dgncur_awet85(:,:,:)
   real(r8), pointer :: qaerwat40(:,:,:)
   real(r8), pointer :: qaerwat55(:,:,:)
   real(r8), pointer :: qaerwat65(:,:,:)
   real(r8), pointer :: qaerwat75(:,:,:)
   real(r8), pointer :: qaerwat85(:,:,:)
!!== KZ_INSITU

   real(r8), allocatable :: maer(:,:,:)      ! aerosol wet mass MR (including water) (kg/kg-air)
   real(r8), allocatable :: hygro(:,:,:)     ! volume-weighted mean hygroscopicity (--)
   real(r8), allocatable :: naer(:,:,:)      ! aerosol number MR (bounded!) (#/kg-air)
   real(r8), allocatable :: dryvol(:,:,:)    ! single-particle-mean dry volume (m3)
   real(r8), allocatable :: drymass(:,:,:)   ! single-particle-mean dry mass  (kg)
   real(r8), allocatable :: dryrad(:,:,:)    ! dry volume mean radius of aerosol (m)

   real(r8), allocatable :: wetrad(:,:,:)    ! wet radius of aerosol (m)
   real(r8), allocatable :: wetvol(:,:,:)    ! single-particle-mean wet volume (m3)
   real(r8), allocatable :: wtrvol(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)

!!== KZ_INSITU
   real(r8), allocatable :: wetrad40(:,:,:)    ! wet radius of aerosol (m)
   real(r8), allocatable :: wetvol40(:,:,:)    ! single-particle-mean wet volume (m3)
   real(r8), allocatable :: wtrvol40(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)

   real(r8), allocatable :: wetrad55(:,:,:)    ! wet radius of aerosol (m)
   real(r8), allocatable :: wetvol55(:,:,:)    ! single-particle-mean wet volume (m3)
   real(r8), allocatable :: wtrvol55(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)

   real(r8), allocatable :: wetrad65(:,:,:)    ! wet radius of aerosol (m)
   real(r8), allocatable :: wetvol65(:,:,:)    ! single-particle-mean wet volume (m3)
   real(r8), allocatable :: wtrvol65(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)

   real(r8), allocatable :: wetrad75(:,:,:)    ! wet radius of aerosol (m)
   real(r8), allocatable :: wetvol75(:,:,:)    ! single-particle-mean wet volume (m3)
   real(r8), allocatable :: wtrvol75(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)

   real(r8), allocatable :: wetrad85(:,:,:)    ! wet radius of aerosol (m)
   real(r8), allocatable :: wetvol85(:,:,:)    ! single-particle-mean wet volume (m3)
   real(r8), allocatable :: wtrvol85(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)

   real(r8) :: rhfix(pcols,pver)        ! relative humidity (0-1)
!!== KZ_INSITU

   real(r8), allocatable :: rhcrystal(:)
   real(r8), allocatable :: rhdeliques(:)
   real(r8), allocatable :: specdens_1(:)

   real(r8), allocatable :: so4_aer(:,:,:)   ! mass MR of sulfate (kg/kg-air)
   real(r8), allocatable :: ss_aer(:,:,:)    ! mass MR of sea salt (kg/kg-air) 

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
!!== KZ_INSITU
      wetrad40(pcols,pver,nmodes),   &
      wetvol40(pcols,pver,nmodes),   &
      wtrvol40(pcols,pver,nmodes),   &
      wetrad55(pcols,pver,nmodes),   &
      wetvol55(pcols,pver,nmodes),   &
      wtrvol55(pcols,pver,nmodes),   &
      wetrad65(pcols,pver,nmodes),   &
      wetvol65(pcols,pver,nmodes),   &
      wtrvol65(pcols,pver,nmodes),   &
      wetrad75(pcols,pver,nmodes),   &
      wetvol75(pcols,pver,nmodes),   &
      wtrvol75(pcols,pver,nmodes),   &
      wetrad85(pcols,pver,nmodes),   &
      wetvol85(pcols,pver,nmodes),   &
      wtrvol85(pcols,pver,nmodes),   &
!!== KZ_INSITU
      rhcrystal(nmodes),           &
      rhdeliques(nmodes),          &
      specdens_1(nmodes)           )
   if ( aw_scheme .gt. 1 ) then
      allocate(so4_aer(pcols,pver,nmodes),ss_aer(pcols,pver,nmodes))
   end if

   maer(:,:,:)     = 0._r8
   hygro(:,:,:)    = 0._r8


   if (list_idx == 0) then
      call pbuf_get_field(pbuf, dgnum_idx,      dgncur_a )
      call pbuf_get_field(pbuf, dgnumwet_idx,   dgncur_awet )
      call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens)
      call pbuf_get_field(pbuf, qaerwat_idx,    qaerwat)
!!== KZ_INSITU
      call pbuf_get_field(pbuf, dgnumwet_rh40_idx,   dgncur_awet40)
      call pbuf_get_field(pbuf, dgnumwet_rh55_idx,   dgncur_awet55)
      call pbuf_get_field(pbuf, dgnumwet_rh65_idx,   dgncur_awet65)
      call pbuf_get_field(pbuf, dgnumwet_rh75_idx,   dgncur_awet75)
      call pbuf_get_field(pbuf, dgnumwet_rh85_idx,   dgncur_awet85)
      call pbuf_get_field(pbuf, qaerwat_rh40_idx,    qaerwat40)
      call pbuf_get_field(pbuf, qaerwat_rh55_idx,    qaerwat55)
      call pbuf_get_field(pbuf, qaerwat_rh65_idx,    qaerwat65)
      call pbuf_get_field(pbuf, qaerwat_rh75_idx,    qaerwat75)
      call pbuf_get_field(pbuf, qaerwat_rh85_idx,    qaerwat85)
!!== KZ_INSITU
   else
      dgncur_a    => dgnumdry_m
      dgncur_awet => dgnumwet_m
      qaerwat     => qaerwat_m
      wetdens     => wetdens_m
!!== KZ_INSITU
      allocate(dgncur_awet40(pcols,pver,nmodes))
      allocate(dgncur_awet55(pcols,pver,nmodes))
      allocate(dgncur_awet65(pcols,pver,nmodes))
      allocate(dgncur_awet75(pcols,pver,nmodes))
      allocate(dgncur_awet85(pcols,pver,nmodes))
      allocate(qaerwat40(pcols,pver,nmodes))
      allocate(qaerwat55(pcols,pver,nmodes))
      allocate(qaerwat65(pcols,pver,nmodes))
      allocate(qaerwat75(pcols,pver,nmodes))
      allocate(qaerwat85(pcols,pver,nmodes))
!!== KZ_INSITU
   end if

   do m = 1, nmodes

      dryvolmr(:,:) = 0._r8

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigmag,  &
         rhcrystal=rhcrystal(m), rhdeliques=rhdeliques(m))

      ! JS added on 08-05-2019, reset rhcrystal to a user-defined lower bound
      rhcrystal(m) = aw_mam_rh_min

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
               if ( aw_scheme .gt. 1 ) then          ! JS added on 08-05-2019, store sulfate and sea salt mass mixing ratio for M7 scheme
                  j = lspectype_amode(l,m)
                  if (trim(specname_amode(j)) .eq. 'sulfate') then 
                     so4_aer(i,k,m) = duma
                  else if (trim(specname_amode(j)) .eq. 'seasalt') then
                     ss_aer(i,k,m) = duma
                  else
                     continue
                  end if
               end if
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
!            rh(i,k) = 0.98_r8                 ! original code
            rh(i,k) = aw_rh_max                ! JS added on 08-05-2019
         endif
         rh(i,k) = max(rh(i,k), 0.0_r8)
!         rh(i,k) = min(rh(i,k), 0.98_r8)      ! original code
         rh(i,k) = min(rh(i,k), aw_rh_max)     ! JS added on 08-05-2019
         if(pergro_mods) then
            cldn_thresh = 0.9998_r8
         else			
            cldn_thresh = 1.0_r8 !original code
         endif
         if (cldn(i,k) .lt. cldn_thresh) then !BSINGH 
            if (aw_l_rh_clearsky) then        ! JS added on 08-05-2019
            rh(i,k) = (rh(i,k) - cldn(i,k)) / (1.0_r8 - cldn(i,k))  ! clear portion
         end if
         end if
         rh(i,k) = max(rh(i,k), 0.0_r8)
      end do
   end do

   ! When reading in namelist, we should check that aw_scheme is >=0 and <=3
   if (aw_scheme .lt. 2) then
         ! solve kohler or kappa-kohler equation
         call modal_aero_wateruptake_sub( ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
                                          hygro, rh, dryvol, wetrad, wetvol, wtrvol )

!!== KZ_INSITU
   rhfix(:,:) = 0.40_r8
   call modal_aero_wateruptake_sub( &
      ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
      hygro, rhfix, dryvol, wetrad40, wetvol40,           &
      wtrvol40)

   rhfix(:,:) = 0.55_r8
   call modal_aero_wateruptake_sub( &
      ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
      hygro, rhfix, dryvol, wetrad55, wetvol55,           &
      wtrvol55)

   rhfix(:,:) = 0.65_r8
   call modal_aero_wateruptake_sub( &
      ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
      hygro, rhfix, dryvol, wetrad65, wetvol65,           &
      wtrvol65)

   rhfix(:,:) = 0.75_r8
   call modal_aero_wateruptake_sub( &
      ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
      hygro, rhfix, dryvol, wetrad75, wetvol75,           &
      wtrvol75)

   rhfix(:,:) = 0.85_r8
   call modal_aero_wateruptake_sub( &
      ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
      hygro, rhfix, dryvol, wetrad85, wetvol85,           &
      wtrvol85)
!!== KZ_INSITU

   else
         ! use M7 schemes: ZSR method for sulfate-sea salt mixing (Jacobson, 1996, 2005)
         !                 A generalised form of the Kelvin equation for sulfate-non sea salt mixing (Stier et al., 2005), used for Aitken mode only
         call m7_wateruptake_scheme( ncol, nmodes, rhcrystal, rhdeliques, rh, t, pmid, maer, &
                                     naer, so4_aer, ss_aer, dryrad, wetrad, wetvol, wtrvol, dryvol )
   end if

   do m = 1, nmodes

      do k = top_lev, pver
         do i = 1, ncol

            dgncur_awet(i,k,m) = dgncur_a(i,k,m) * (wetrad(i,k,m)/dryrad(i,k,m))
            qaerwat(i,k,m)     = rhoh2o*naer(i,k,m)*wtrvol(i,k,m)
!!== KZ_INSITU
            dgncur_awet40(i,k,m) = dgncur_a(i,k,m) * (wetrad40(i,k,m)/dryrad(i,k,m))
            dgncur_awet55(i,k,m) = dgncur_a(i,k,m) * (wetrad55(i,k,m)/dryrad(i,k,m))
            dgncur_awet65(i,k,m) = dgncur_a(i,k,m) * (wetrad65(i,k,m)/dryrad(i,k,m))
            dgncur_awet75(i,k,m) = dgncur_a(i,k,m) * (wetrad75(i,k,m)/dryrad(i,k,m))
            dgncur_awet85(i,k,m) = dgncur_a(i,k,m) * (wetrad85(i,k,m)/dryrad(i,k,m))
            qaerwat40(i,k,m)     = rhoh2o*naer(i,k,m)*wtrvol40(i,k,m)
            qaerwat55(i,k,m)     = rhoh2o*naer(i,k,m)*wtrvol55(i,k,m)
            qaerwat65(i,k,m)     = rhoh2o*naer(i,k,m)*wtrvol65(i,k,m)
            qaerwat75(i,k,m)     = rhoh2o*naer(i,k,m)*wtrvol75(i,k,m)
            qaerwat85(i,k,m)     = rhoh2o*naer(i,k,m)*wtrvol85(i,k,m)
!!== KZ_INSITU

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

!!== KZ_INSITU
   deallocate( wetrad40, wetvol40, wtrvol40) 
   deallocate( wetrad55, wetvol55, wtrvol55) 
   deallocate( wetrad65, wetvol65, wtrvol65) 
   deallocate( wetrad75, wetvol75, wtrvol75) 
   deallocate( wetrad85, wetvol85, wtrvol85) 

   if (list_idx == 0) then
   else
      deallocate(dgncur_awet40)
      deallocate(dgncur_awet55)
      deallocate(dgncur_awet65)
      deallocate(dgncur_awet75)
      deallocate(dgncur_awet85)
      deallocate(qaerwat40)
      deallocate(qaerwat55)
      deallocate(qaerwat65)
      deallocate(qaerwat75)
      deallocate(qaerwat85)
   end if


!!== KZ_INSITU

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
            ! aw_scheme = 1 - solve kappa-kohler equation based on Petters and Kreidenweis, 2007, ACP
            !                 use bisection method to do the root-finding
            !           = 0 - olve kohler equation in the original code
            call modal_aero_kohler(dryrad(i:i,k,m), hygro(i:i,k,m), rh(i:i,k), wetrad(i:i,k,m), 1, aw_scheme)

            wetrad(i,k,m) = max(wetrad(i,k,m), dryrad(i,k,m))
            wetvol(i,k,m) = pi43*wetrad(i,k,m)**3
            wetvol(i,k,m) = max(wetvol(i,k,m), dryvol(i,k,m))
            wtrvol(i,k,m) = wetvol(i,k,m) - dryvol(i,k,m)
            wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0_r8)

            ! apply simple treatment of deliquesence/crystallization hysteresis
            ! for rhcrystal < rh < rhdeliques, aerosol water is a fraction of
            ! the "upper curve" value, and the fraction is a linear function of rh
            if ( aw_scheme .lt. 2 ) then            ! apply hyst to MAM and k-kohler only
            if (rh(i,k) < rhcrystal(m)) then
               wetrad(i,k,m) = dryrad(i,k,m)
               wetvol(i,k,m) = dryvol(i,k,m)
               wtrvol(i,k,m) = 0.0_r8
            else if (rh(i,k) < rhdeliques(m)) then
                  if ( aw_l_hysteresis ) then       ! JS added on 08-05-2019
               wtrvol(i,k,m) = wtrvol(i,k,m)*hystfac*(rh(i,k) - rhcrystal(m))
               wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0_r8)
               wetvol(i,k,m) = dryvol(i,k,m) + wtrvol(i,k,m)
               wetrad(i,k,m) = (wetvol(i,k,m)/pi43)**third
            end if
               end if
            end if

         end do  ! columns
      end do     ! levels

   end do ! modes

end subroutine modal_aero_wateruptake_sub

!-----------------------------------------------------------------------
      subroutine modal_aero_kohler(   &
          rdry_in, hygro, s, rwet_out, im, aw_scheme )

! calculates equlibrium radius r of haze droplets as function of
! dry particle mass and relative humidity s using kohler solution
! given in pruppacher and klett (eqn 6-35)

! for multiple aerosol types, assumes an internal mixture of aerosols

      implicit none

! arguments
      integer :: im         ! number of grid points to be processed
      real(r8) :: rdry_in(:)    ! aerosol dry radius (m)
      real(r8) :: hygro(:)      ! aerosol volume-mean hygroscopicity (--)
      real(r8) :: s(:)          ! relative humidity (1 = saturated)
      real(r8) :: rwet_out(:)   ! aerosol wet radius (m)
      integer  :: aw_scheme     ! 1 - k-kohler equation ; other - kohler equation (MAM default)

! local variables
      integer, parameter :: imax=200
      integer :: i, n, nsol

      real(r8) :: a, b
      real(r8) :: p40(imax),p41(imax),p42(imax),p43(imax) ! coefficients of polynomial
      real(r8) :: p30(imax),p31(imax),p32(imax) ! coefficients of polynomial
      real(r8) :: p
      real(r8) :: r3, r4
      real(r8) :: r(im)         ! wet radius (microns)
      real(r8) :: rdry(imax)    ! radius of dry particle (microns)
      real(r8) :: ss            ! relative humidity (1 = saturated)
      real(r8) :: slog(imax)    ! log relative humidity
      real(r8) :: vol(imax)     ! total volume of particle (microns**3)
      real(r8) :: xi, xr

      complex(r8) :: cx4(4,imax),cx3(3,imax)

      real(r8), parameter :: eps = 1.e-4_r8
      real(r8), parameter :: mw = 18._r8
      real(r8), parameter :: pi = 3.14159_r8
      real(r8), parameter :: rhow = 1._r8
      real(r8), parameter :: third = 1._r8/3._r8
      real(r8), parameter :: ugascon = 8.3e7_r8
      real(r8)            :: surften, tair

      if  (aw_scheme .eq. 1) then
          surften = 72._r8
          tair = 298.15_r8
      else
          surften = 76._r8
          tair = 273._r8
      end if

!     effect of organics on surface tension is neglected
      a=2.e4_r8*mw*surften/(ugascon*tair*rhow)

      do i=1,im
           rdry(i) = rdry_in(i)*1.0e6_r8   ! convert (m) to (microns)
           vol(i) = rdry(i)**3          ! vol is r**3, not volume
           b = vol(i)*hygro(i)

!          quartic
           ss=min(s(i),1._r8-eps)
           ss=max(ss,1.e-10_r8)
           slog(i)=log(ss)
           p43(i)=-a/slog(i)
           p42(i)=0._r8
           p41(i)=b/slog(i)-vol(i)
           p40(i)=a*vol(i)/slog(i)
!          cubic for rh=1
           p32(i)=0._r8
           p31(i)=-b/a
           p30(i)=-vol(i)
      end do


       do 100 i=1,im

!       if(vol(i).le.1.e-20)then
        if(vol(i).le.1.e-12_r8)then
           r(i)=rdry(i)
           go to 100
        endif

        p=abs(p31(i))/(rdry(i)*rdry(i))
        if(p.lt.eps)then
!          approximate solution for small particles
           r(i)=rdry(i)*(1._r8+p*third/(1._r8-slog(i)*rdry(i)/a))
        else
           if ( aw_scheme .eq. 1 ) then     ! JS added on 08-05-2019
              call solve_k_kohler( rdry(i), s(i), a, hygro(i), r(i) )
           else
           call makoh_quartic(cx4(1,i),p43(i),p42(i),p41(i),p40(i),1)
!          find smallest real(r8) solution
           r(i)=1000._r8*rdry(i)
           nsol=0
           do n=1,4
              xr=real(cx4(n,i))
              xi=aimag(cx4(n,i))
              if(abs(xi).gt.abs(xr)*eps) cycle  
              if(xr.gt.r(i)) cycle  
              if(xr.lt.rdry(i)*(1._r8-eps)) cycle  
              if(xr.ne.xr) cycle  
              r(i)=xr
              nsol=n
           end do  
           if(nsol.eq.0)then
              write(iulog,*)   &
               'ccm kohlerc - no real(r8) solution found (quartic)'
              write(iulog,*)'roots =', (cx4(n,i),n=1,4)
              write(iulog,*)'p0-p3 =', p40(i), p41(i), p42(i), p43(i)
              write(iulog,*)'rh=',s(i)
              write(iulog,*)'setting radius to dry radius=',rdry(i)
              r(i)=rdry(i)
!             stop
           endif
        endif
        endif

! JS commented on 08-5-2019: this section is called only when rh is very close to 1.0
!                            should never be called unless the rh_ceiling value is changed 
        if(s(i).gt.1._r8-eps)then
!          save quartic solution at s=1-eps
           r4=r(i)
!          cubic for rh=1
           p=abs(p31(i))/(rdry(i)*rdry(i))
           if(p.lt.eps)then
              r(i)=rdry(i)*(1._r8+p*third)
           else
              if ( aw_scheme .eq. 1 ) then     ! JS added on 08-05-2019
                 call solve_k_kohler( rdry(i), s(i), a, hygro(i), r(i) )
              else
              call makoh_cubic(cx3,p32,p31,p30,im)
!             find smallest real(r8) solution
              r(i)=1000._r8*rdry(i)
              nsol=0
              do n=1,3
                 xr=real(cx3(n,i))
                 xi=aimag(cx3(n,i))
                 if(abs(xi).gt.abs(xr)*eps) cycle  
                 if(xr.gt.r(i)) cycle  
                 if(xr.lt.rdry(i)*(1._r8-eps)) cycle  
                 if(xr.ne.xr) cycle  
                 r(i)=xr
                 nsol=n
              end do  
              if(nsol.eq.0)then
                 write(iulog,*)   &
                  'ccm kohlerc - no real(r8) solution found (cubic)'
                 write(iulog,*)'roots =', (cx3(n,i),n=1,3)
                 write(iulog,*)'p0-p2 =', p30(i), p31(i), p32(i)
                 write(iulog,*)'rh=',s(i)
                 write(iulog,*)'setting radius to dry radius=',rdry(i)
                 r(i)=rdry(i)
!                stop
              endif
           endif
           endif
           r3=r(i)
!          now interpolate between quartic, cubic solutions
           r(i)=(r4*(1._r8-s(i))+r3*(s(i)-1._r8+eps))/eps
        endif

  100 continue

! bound and convert from microns to m
      do i=1,im
         r(i) = min(r(i),30._r8) ! upper bound based on 1 day lifetime
         rwet_out(i) = r(i)*1.e-6_r8
      end do

      return
      end subroutine modal_aero_kohler


!-----------------------------------------------------------------------
      subroutine makoh_cubic( cx, p2, p1, p0, im )
!
!     solves  x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2 are real
!
      integer, parameter :: imx=200
      integer :: im
      real(r8) :: p0(imx), p1(imx), p2(imx)
      complex(r8) :: cx(3,imx)

      integer :: i
      real(r8) :: eps, q(imx), r(imx), sqrt3, third
      complex(r8) :: ci, cq, crad(imx), cw, cwsq, cy(imx), cz(imx)

      save eps
      data eps/1.e-20_r8/

      third=1._r8/3._r8
      ci=cmplx(0._r8,1._r8,r8)
      sqrt3=sqrt(3._r8)
      cw=0.5_r8*(-1+ci*sqrt3)
      cwsq=0.5_r8*(-1-ci*sqrt3)

      do i=1,im
      if(p1(i).eq.0._r8)then
!        completely insoluble particle
         cx(1,i)=(-p0(i))**third
         cx(2,i)=cx(1,i)
         cx(3,i)=cx(1,i)
      else
         q(i)=p1(i)/3._r8
         r(i)=p0(i)/2._r8
         crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
         crad(i)=sqrt(crad(i))

         cy(i)=r(i)-crad(i)
         if (abs(cy(i)).gt.eps) cy(i)=cy(i)**third
         cq=q(i)
         cz(i)=-cq/cy(i)

         cx(1,i)=-cy(i)-cz(i)
         cx(2,i)=-cw*cy(i)-cwsq*cz(i)
         cx(3,i)=-cwsq*cy(i)-cw*cz(i)
      endif
      enddo

      return
      end subroutine makoh_cubic


!-----------------------------------------------------------------------
      subroutine makoh_quartic( cx, p3, p2, p1, p0, im )

!     solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2, p3 are real
!
      integer, parameter :: imx=200
      integer :: im
      real(r8) :: p0(imx), p1(imx), p2(imx), p3(imx)
      complex(r8) :: cx(4,imx)

      integer :: i
      real(r8) :: third, q(imx), r(imx)
      complex(r8) :: cb(imx), cb0(imx), cb1(imx),   &
                     crad(imx), cy(imx), czero


      czero=cmplx(0.0_r8,0.0_r8,r8)
      third=1._r8/3._r8

      do 10 i=1,im

      q(i)=-p2(i)*p2(i)/36._r8+(p3(i)*p1(i)-4*p0(i))/12._r8
      r(i)=-(p2(i)/6)**3+p2(i)*(p3(i)*p1(i)-4*p0(i))/48._r8   &
       +(4*p0(i)*p2(i)-p0(i)*p3(i)*p3(i)-p1(i)*p1(i))/16

      crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
      crad(i)=sqrt(crad(i))

      cb(i)=r(i)-crad(i)
      if(cb(i).eq.czero)then
!        insoluble particle
         cx(1,i)=(-p1(i))**third
         cx(2,i)=cx(1,i)
         cx(3,i)=cx(1,i)
         cx(4,i)=cx(1,i)
      else
         cb(i)=cb(i)**third

         cy(i)=-cb(i)+q(i)/cb(i)+p2(i)/6

         cb0(i)=sqrt(cy(i)*cy(i)-p0(i))
         cb1(i)=(p3(i)*cy(i)-p1(i))/(2*cb0(i))

         cb(i)=p3(i)/2+cb1(i)
         crad(i)=cb(i)*cb(i)-4*(cy(i)+cb0(i))
         crad(i)=sqrt(crad(i))
         cx(1,i)=(-cb(i)+crad(i))/2._r8
         cx(2,i)=(-cb(i)-crad(i))/2._r8

         cb(i)=p3(i)/2-cb1(i)
         crad(i)=cb(i)*cb(i)-4*(cy(i)-cb0(i))
         crad(i)=sqrt(crad(i))
         cx(3,i)=(-cb(i)+crad(i))/2._r8
         cx(4,i)=(-cb(i)-crad(i))/2._r8
      endif
   10 continue

      return
      end subroutine makoh_quartic

!----------------------------------------------------------------------

      ! this subroutine uses a bisection method to solve the k-kohler equation
      subroutine solve_k_kohler( rd, rh, kelvin_A, hygro, rw )
      real(r8), intent(in)  :: rd,       &  ! dry radius (m)
                               rh,       &  ! relative humidity
                               kelvin_A, &  ! constant in the Kelvin term
                               hygro        ! volume-mean hygroscopicity
      real(r8), intent(out) :: rw           ! wet radius (m)

      ! Local variable
      real(r8) :: left, right, middle, f_left, f_right, f_middle
      integer  :: niter
      integer, parameter  :: niter_max = 1000
      real(r8), parameter :: toler = epsilon(1._r8)  ! tolerance for bisection 

      left    = rd
      right   = rd * 1.e3_r8
      f_left  = k_kohler(rh, kelvin_A, hygro, rd, left)
      f_right = k_kohler(rh, kelvin_A, hygro, rd, right)

      ! start the iteration
      niter       = 0
      do
         middle   = (left + right) / 2.0_r8
         f_middle = k_kohler(rh, kelvin_A, hygro, rd, middle)

         if (abs(left-right)/middle < toler) then
            rw = middle
            exit
         else if (f_middle .eq. 0._r8) then
            rw = middle
            exit
         else if (f_middle*f_left < 0.0_r8) then
            right   = middle
            f_right = f_middle
         else
            left    = middle
            f_left  = f_middle
         end if
         niter = niter + 1
         if (niter > niter_max) then
            write(iulog,*) 'Too many iterations for solving k-kohler equation...'
            write(iulog,*) 'Reset wet radius = dry radius...'
            rw = rd
            exit
         end if
      end do

      end subroutine solve_k_kohler

!----------------------------------------------------------------------

      ! this subroutine is used to evaluate the k-kohler equation based on Petters and Kreidenweis, 2007, ACP 
      real(r8) function k_kohler( rh, A, B, rd, r )
      implicit none
      real(r8), intent(in) :: rh, A, B, rd, r

      ! Evaluate the function F = ln(RH)-[ln(r^3-rd^3)-ln(r^3-(1-k)rd^3)+kelvin_term]

      k_kohler = rh - (r**3-rd**3)/(r**3-(1._r8-B)*(rd**3))*exp(A/r)

      end function k_kohler

!----------------------------------------------------------------------

      ! this subroutine is heavily borrowed from M7 codes to calculate the wet radius of aerosols due to water uptake
      subroutine m7_wateruptake_scheme( ncol, nmodes, rhcrystal, rhdeliques, &
                                        rh, t, pmid, maer, naer, so4_aer, ss_aer,  &
                                        dryrad, wetrad, wetvol, wtrvol, dryvol )

      use modal_aero_data, only: modeptr_pcarbon
      use physconst,       only: avogad, r_universal, mwdry
      use modal_aero_data, only: specmw_so4_amode, specmw_seasalt_amode, &
                                 specdens_so4_amode, specdens_seasalt_amode

      implicit none
      integer,  intent(in)  :: ncol, nmodes
      real(r8), intent(in)  :: rhcrystal(nmodes), rhdeliques(nmodes)
      real(r8), pointer     :: t(:,:)                     ! temperatures (K)
      real(r8), pointer     :: pmid(:,:)                  ! layer pressure (Pa)
      real(r8), intent(in)  :: rh(pcols,pver)             ! relative humidity  
      real(r8), intent(in)  :: maer(pcols,pver,nmodes)    ! aerosol mass MR (kg/kg-air)
      real(r8), intent(in)  :: naer(pcols,pver,nmodes)    ! aerosol number MR (#/kg-air)
      real(r8), intent(in)  :: so4_aer(pcols,pver,nmodes) ! aerosol mass MR (kg/kg-air)
      real(r8), intent(in)  :: ss_aer(pcols,pver,nmodes)  ! aerosol number MR (#/kg-air)
      real(r8), intent(in)  :: dryvol(pcols,pver,nmodes), &
                               dryrad(pcols,pver,nmodes)
      real(r8), intent(out) :: wetvol(pcols,pver,nmodes), &
                               wtrvol(pcols,pver,nmodes), &
                               wetrad(pcols,pver,nmodes)

      ! local vairables
      integer               :: i, k, n
      real(r8)              :: ztk, zrh, zpa, zlnm, zlnm2, zss2,    &
                               ztk2, zln3, zln32, zmso4h2o, zwso4,  &
                               zdso4h2o, h2omass, hystfac, zdryvol, &
                               zdryvol_mean, rdry, zwetvol, zaw,    &
                               zwtrvol, rwet, gf, pttn
      real(r8)              :: zmm(5), zmmr(5), zmol(5), zmo(4), zmmt(4)
      real(r8)              :: dso4, dss, wso4, wss, dair

      !
      ! From M7 module
      !
      real(r8), parameter   :: dna2so4 = 2.68E3_r8,                  &  ! Density of Na2SO4   [kg/m3]
                               dnahso4 = 2.435E3_r8,                 &  ! Density of NaHSO4   [kg/m3]
                               dh2so4  = 1.841E3_r8,                 &  ! Density of H2SO4    [kg/m3]
                               wh2so4  = 98.0734_r8,                 &  ! Molecular weight of H2SO4    [g/mol]
                               wna2so4 = 142.0376_r8,                &  ! Molecular weight of Na2SO4   [g/mol]
                               wnahso4 = 120.0555_r8                    ! Molecular weight of NaHSO4   [g/mol]

      real(r8), parameter   :: wvb(17) = &
                (/   95.80188_r8,      -28.5257_r8,    -1.082153_r8, &
                    0.1466501_r8,     -20627.51_r8,    0.0461242_r8, &
                    -0.003935_r8,      -3.36115_r8,  -0.00024137_r8, &
                  0.067938345_r8, 0.00000649899_r8,  8616124.373_r8, &
                  1.168155578_r8,  -0.021317481_r8,  0.000270358_r8, &
                -1353332314.0_r8,  -0.002403805_r8  /)

      real(r8), parameter   :: gmb(9) = &
                (/ 1.036391467_r8,   0.00728531_r8, -0.011013887_r8,  &
                  -0.068887407_r8,  0.001047842_r8,  0.001049607_r8,  &
                   0.000740534_r8, -1.081202685_r8, -0.0000029113_r8 /)

      dso4 = specdens_so4_amode                                     ! Density of sulfate  [kg/m3]
      dss  = specdens_seasalt_amode                                 ! Density of sea salt [kg/m3] 
      wso4 = specmw_so4_amode                                       ! Molecular weight of sulfate  [g/mol]
      wss  = specmw_seasalt_amode                                   ! Molecular weight of sea salt [g/mol]

mode_loop: do n = 1, nmodes

         if (n .eq. modeptr_pcarbon) then
            ! M7 scheme is not applicable for primary carbon mode, because of no sea salt or sulfate 
            ! Assume no water uptake occurs in this mode
            wetrad(:,:,n) = dryrad(:,:,n)
            wtrvol(:,:,n) = 0._r8
            cycle mode_loop
         end if

         hystfac = 1.0_r8 / max(1.0e-5_r8, (rhdeliques(n) - rhcrystal(n)))

         do k = top_lev, pver
            do i = 1, ncol
               dair = pmid(i,k) / (t(i,k) * r_universal / mwdry)                                  ! dry density of air [kg/m3]
               if ( so4_aer(i,k,n) .gt. 0._r8 .and. ss_aer(i,k,n)*1E9*dair .lt. 1.e-15 ) then     ! M7 uses [ug/m3] for sea salt
                  ! The calculations of the ambient particle properties are based on parameterisations of the mass of 
                  ! sulfate and density derived by Julian Wilson from a regression analysis of results of solving
                  ! the generalised Kelvin equation using (F. J. Zeleznik, J. Phys. Chem. Ref. Data 20, 1157, 1991), 
                  ! for an H2SO4-H2O mixture
                  pttn  = so4_aer(i,k,n)*avogad / (naer(i,k,n)*wso4)      !molecules of sulfate per particle
                  ztk   = t(i,k)
                  ztk   = MAX(ztk, 240._r8)
                  zrh   = MAX(rh(i,k), 0.05_r8)
                  zrh   = MIN(rh(i,k), 0.90_r8)
                  zpa   = pmid(i,k)
                  zlnm  = LOG(pttn)
                  zlnm2 = zlnm*zlnm
                  zss2  = zrh**2
                  ztk2  = ztk*ztk
                  zln3  = zlnm/3.0_r8
                  zln32 = zln3*zln3
                  ! Percentage by mass of sulfate in a H2O-H2SO4 particle [%] 
                  ! (Here we ignore any insoluble mass.)
                  zwso4 = wvb(1) + wvb(2)*zlnm + wvb(3)*zrh*zlnm + &
                          wvb(4)*ztk*zlnm + wvb(5)*zrh/ztk + &
                          wvb(6)*zlnm2*zrh + wvb(7)*zlnm2*ztk + &
                          wvb(8)*zlnm*zss2 + wvb(9)*zlnm*ztk2 + &
                          wvb(10)*zlnm2*zss2 + wvb(11)*zlnm2*ztk2 + &
                          wvb(12)*zss2/ztk2 + wvb(13)*zlnm2 + &
                          wvb(14)*zlnm2*zlnm + wvb(15)*zlnm2*zlnm2 + &
                          wvb(16)*zss2*zrh/(ztk2*ztk) + &
                          wvb(17)*LOG(zrh*ztk/zpa)
                  ! Mass of sulfate + water in an average particle, [kg]
                  zmso4h2o      = so4_aer(i,k,n) / (naer(i,k,n)*zwso4/100.0_r8)
                  ! zdso4h2o    = density of sulfate-h2o fraction of a particle [g/cm3]
                  zdso4h2o      = gmb(1) + gmb(2)*zwso4 + gmb(3)*zln3 + &
                                  gmb(4)*zrh + gmb(5)*ztk + gmb(6)*zln32 + &
                                  gmb(7)*zln3/zrh + gmb(8)*zln3/ztk + &
                                  gmb(9)*ztk2
                  ! Limits for zdso4h2o: H2O and pure H2SO4:
                  zdso4h2o      = max(zdso4h2o,rhoh2o/1E3)
                  zdso4h2o      = min(zdso4h2o,dh2so4/1E3)
                  ! Dry volume of a particle without sulfate [m3]
                  zdryvol_mean  = dryvol(i,k,n) - so4_aer(i,k,n) / naer(i,k,n) / dso4 
                  ! Wet volume of a particle with sulfate-water fraction [m3]
                  wetvol(i,k,n) = zdryvol_mean + zmso4h2o * 1.E-3 / zdso4h2o 
                  wetvol(i,k,n) = max(wetvol(i,k,n), zdryvol_mean)
                  ! Wet radius  [m]
                  wetrad(i,k,n) = (wetvol(i,k,n)/pi43)**third
                  wetrad(i,k,n) = min(wetrad(i,k,n), 30.e-6_r8)
                  wetrad(i,k,n) = max(wetrad(i,k,n), dryrad(i,k,n))
                  ! Volume of absorbed aerosol water [m3]
                  wtrvol(i,k,n) = pi43*(wetrad(i,k,n)**3._r8) - dryvol(i,k,n)
                  wtrvol(i,k,n) = max(wtrvol(i,k,n), 0.0_r8)

               else
                 
                  if ( rh(i,k) .lt. aw_m7_rh_min ) then
                     ! If relative humidity is lower than critical RH, the aerosol particles stay dry
                     wetrad(i,k,n) = dryrad(i,k,n)
                     wetvol(i,k,n) = pi43*wetrad(i,k,n)**3
                     wtrvol(i,k,n) = 0.0_r8
                  else
                     ! Apply ZSR method for sulfate-sea salt mixing, ignore the impact of other species 
                     !--- 1.1) Calculate initial concentrations of the compounds
                     zmm(1)  = ss_aer(i,k,n) * 1.E3 / wss     ! n(Na) [mole/kg-air]
                     zmm(2)  = zmm(1)                         ! n(Cl) [mole/kg-air]
                     zmm(3)  = so4_aer(i,k,n) * 1.E3 / wso4   ! n(H2SO4) [mole/kg-air]
                     zmm(4)  = 0._r8                          ! n(HSO4)  [mole/kg-air]
                     !--- 1.2) Calculation of the concentration of the different species:
                     !         The ions are supposed to be arranged such that sodium is associated to sulphate 
                     !         first in Na2SO4, the remaining sodium is associated to Cl in form of NaCl.
                     zmmt(1) = zmm(1)                         ! n(Na)
                     zmmt(3) = zmm(3)                         ! n(SO4)
                     zmmr(3) = min(zmmt(1)/2._r8 , zmmt(3))   ! n(Na2SO4)
                     zmmt(1) = zmmt(1)-2._r8*zmmr(3)          ! Remaining n(Na) after association 
                                                              ! with Na2SO4: n(Na)=n(Na)-2*n(Na2SO4)
                     zmmr(1) = min(zmm(2), zmmt(1))           ! n(NaCl) 
                     zmm(2)  = zmmr(1)                        ! n(Cl) bound in NaCl, the rest is
                                                              ! assumed to evaporate in presence of SO4
                     zmmr(2) = 0._r8                          ! n(NaHSO4)
                     zmmr(4) = zmm(3) - zmmr(2) - zmmr(3)     ! n(H2-SO4)(t)=n(H2SO4)(t0)-n(NaHSO4)-n(Na2SO4)
                                                              ! as n(H2SO4)(t0)=n(SO4--)(t0)
                     !--- 1.3) Total aerosol dry volume [m3/kg-air]:
                     zdryvol = dryvol(i,k,n) * naer(i,k,n)   - &
                               so4_aer(i,k,n) / dso4         - &
                               ss_aer(i,k,n)  / dss          + &
                               zmmr(1)*wss*1.E-3/dss         + &
                               zmmr(2)*wnahso4*1.E-3/dnahso4 + &
                               zmmr(3)*wna2so4*1.E-3/dna2so4 + &
                               zmmr(4)*wh2so4*1.E-3/dh2so4 
                     !--- 1.4) Mean aerosol dry volume [m3]:
                     zdryvol_mean = zdryvol / naer(i,k,n)
                     !--- 1.5) Dry radius [m]:
                     rdry         = (zdryvol_mean/pi43)**third
                     !--- 2.1) Molality as function of the water activity:
                     !         Currently sulfate is assumed to be fully dissociated, i.e. zmmr(2)=0. and zmo(2) is not calculated.
                     !         Assume water activity  = relative humidity
                     zaw        = rh(i,k)
                     if (aw_scheme .eq. 2) then
                         ! old M7 coefficients from Jacobson 1996, Eq. 29
                         if ( .not. aw_l_m7_rh_min ) then
                            !--- NaCl:
                            zmo(1) = (  - 1.918004E2_r8            + 2.001540E3_r8*zaw            &
                                        - 8.557205E3_r8*zaw**2     + 1.987670E4_r8*zaw**3         &
                                        - 2.717192E4_r8*zaw**4     + 2.187103E4_r8*zaw**5         &
                                        - 9.591577E3_r8*zaw**6     + 1.763672E3_r8*zaw**7      )**2
                            !--- NaHSO4:
                            zmo(2) = (  + 4.662777E0_r8            - 1.128472E1_r8*zaw            &
                                        + 7.049464E1_r8*zaw**2     - 2.788050E2_r8*zaw**3         &
                                        + 6.103105E2_r8*zaw**4     - 7.409417E2_r8*zaw**5         &
                                        + 4.614577E2_r8*zaw**6     - 1.150735E2_r8*zaw**7      )**2
                            !--- Na2SO4:
                            zmo(3) = (  - 3.295311E3_r8            + 3.188349E4_r8*zaw            &
                                        - 1.305168E5_r8*zaw**2     + 2.935608E5_r8*zaw**3         &
                                        - 3.920423E5_r8*zaw**4     + 3.109519E5_r8*zaw**5         &
                                        - 1.356439E5_r8*zaw**6     + 2.510249E4_r8*zaw**7      )**2
                            !--- H2SO4
                            zmo(4) = (  + 5.611895_r8              - 1.387446E1_r8*zaw            &
                                        + 1.750682E1_r8*zaw**2     + 7.138146E1_r8*zaw**3         &
                                        - 3.109173E2_r8*zaw**4     + 4.662288E2_r8*zaw**5         &
                                        - 3.128612E2_r8*zaw**6     + 7.76097E1_r8*zaw**7       )**2
                         else
                            ! apply RHmin and highest molality constraints
                            !--- NaCl:
!                            zaw    = rh(i,k)
!                            zaw    = max(zaw, 0.43_r8)
                            zmo(1) = (  - 1.918004E2_r8            + 2.001540E3_r8*zaw            &
                                        - 8.557205E3_r8*zaw**2     + 1.987670E4_r8*zaw**3         &
                                        - 2.717192E4_r8*zaw**4     + 2.187103E4_r8*zaw**5         &
                                        - 9.591577E3_r8*zaw**6     + 1.763672E3_r8*zaw**7      )**2
!                            zmo(1) = min(zmo(1), 13.6_r8)
                            !--- NaHSO4:
!                            zaw    = rh(i,k)
!                            zaw    = max(zaw, 0._r8)
                            zmo(2) = (  + 4.662777E0_r8            - 1.128472E1_r8*zaw            &
                                        + 7.049464E1_r8*zaw**2     - 2.788050E2_r8*zaw**3         &
                                        + 6.103105E2_r8*zaw**4     - 7.409417E2_r8*zaw**5         &
                                        + 4.614577E2_r8*zaw**6     - 1.150735E2_r8*zaw**7      )**2
!                            zmo(2) = min(zmo(2), 20.5_r8)
                            !--- Na2SO4:
!                            zaw    = rh(i,k)
!                            zaw    = max(zaw, 0.51_r8)
                            zmo(3) = (  - 3.295311E3_r8            + 3.188349E4_r8*zaw            &
                                        - 1.305168E5_r8*zaw**2     + 2.935608E5_r8*zaw**3         &
                                        - 3.920423E5_r8*zaw**4     + 3.109519E5_r8*zaw**5         &
                                        - 1.356439E5_r8*zaw**6     + 2.510249E4_r8*zaw**7      )**2
!                            zmo(3) = min(zmo(3), 12.8_r8)
                            !--- H2SO4  
!                            zaw    = rh(i,k)
!                            zaw    = max(zaw, 0._r8)
                            zmo(4) = (  + 5.611895_r8              - 1.387446E1_r8*zaw            &
                                        + 1.750682E1_r8*zaw**2     + 7.138146E1_r8*zaw**3         &
                                        - 3.109173E2_r8*zaw**4     + 4.662288E2_r8*zaw**5         &
                                        - 3.128612E2_r8*zaw**6     + 7.76097E1_r8*zaw**7       )**2
!                            zmo(4) = min(zmo(4), 29.1_r8)
                         end if
                     else
                         ! new coefficients and formula from Jacobson's book, 2005, Eq. 17.66 and Table B.10
                         if ( .not. aw_l_m7_rh_min ) then
                            !--- NaCl:
                            zmo(1) = (     5.875248E1_r8           - 1.8781997E2_r8*zaw           &
                                        + 2.7211377E2_r8*zaw**2    - 1.8458287E2_r8*zaw**3        &
                                        +  4.153689E1_r8*zaw**4                                   )
                            !--- NaHSO4:
                            zmo(2) = (    1.8457001681E2_r8        - 1.6147765817E3_r8*zaw        &
                                        +  8.444076586E3_r8*zaw**2 - 2.6813441936E4_r8*zaw**3     &
                                        + 5.0821277356E4_r8*zaw**4 - 5.5964847603E4_r8*zaw**5     &
                                        + 3.2945298603E4_r8*zaw**6 -  8.002609678E3_r8*zaw**7     )
                            !--- Na2SO4:
                            zmo(3) = (     5.5983158E2_r8          - 2.56942664E3_r8*zaw          &
                                        + 4.47450201E3_r8*zaw**2   - 3.45021842E3_r8*zaw**3       &
                                        +  9.8527913E2_r8*zaw**4                                  )
                            !--- H2SO4
                            zmo(4) = (    3.0391387536E1_r8        - 1.8995058929E2_r8*zaw        &
                                        + 9.7428231047E2_r8*zaw**2 - 3.1680155761E3_r8*zaw**3     &
                                        + 6.1400925314E3_r8*zaw**4 - 6.9116348199E3_r8*zaw**5     &
                                        + 4.1631475226E3_r8*zaw**6 - 1.0383424491E3_r8*zaw**7     )
                         else
                            ! apply RHmin and highest molality constraints
                            !--- NaCl:
!                            zaw    = rh(i,k)
!                            zaw    = max(zaw, 0.47_r8)
                            zmo(1) = (     5.875248E1_r8           - 1.8781997E2_r8*zaw           &
                                        + 2.7211377E2_r8*zaw**2    - 1.8458287E2_r8*zaw**3        &
                                        +  4.153689E1_r8*zaw**4                                   )
!                            zmo(1) = min(zmo(1), 13.5_r8)
                            !--- NaHSO4:
!                            zaw    = rh(i,k)
!                            zaw    = max(zaw, 0.019_r8)
                            zmo(2) = (    1.8457001681E2_r8        - 1.6147765817E3_r8*zaw        &
                                        +  8.444076586E3_r8*zaw**2 - 2.6813441936E4_r8*zaw**3     &
                                        + 5.0821277356E4_r8*zaw**4 - 5.5964847603E4_r8*zaw**5     &
                                        + 3.2945298603E4_r8*zaw**6 -  8.002609678E3_r8*zaw**7     )
!                            zmo(2) = min(zmo(2), 158._r8)
                            !--- Na2SO4:
!                            zaw    = rh(i,k)
!                            zaw    = max(zaw, 0.58_r8)
                            zmo(3) = (     5.5983158E2_r8          - 2.56942664E3_r8*zaw          &
                                        + 4.47450201E3_r8*zaw**2   - 3.45021842E3_r8*zaw**3       &
                                        +  9.8527913E2_r8*zaw**4                                  )
!                            zmo(3) = min(zmo(3), 13.1_r8)
                            !--- H2SO4
!                            zaw    = rh(i,k)
!                            zaw    = max(zaw, 0._r8)
                            zmo(4) = (    3.0391387536E1_r8        - 1.8995058929E2_r8*zaw        &
                                        + 9.7428231047E2_r8*zaw**2 - 3.1680155761E3_r8*zaw**3     &
                                        + 6.1400925314E3_r8*zaw**4 - 6.9116348199E3_r8*zaw**5     &
                                        + 4.1631475226E3_r8*zaw**6 - 1.0383424491E3_r8*zaw**7     )
!                            zmo(4) = min(zmo(4), 30.4_r8)
                         end if
                     end if
                     !--- 2.2) Calculation of the water content in kg-water/kg-air:
                     h2omass       = zmmr(1)/zmo(1) + zmmr(2)/zmo(2) + zmmr(3)/zmo(3) + zmmr(4)/zmo(4)
                     !--- 2.3) Mean water volume [m3]
                     zwtrvol       = h2omass / ( rhoh2o * naer(i,k,n) )
                     zwtrvol       = max(zwtrvol, 0.0_r8)
                     !--- 2.4) Mean wet volume [m3]
                     zwetvol       = zdryvol_mean + zwtrvol
                     zwetvol       = max(zwetvol, zdryvol_mean)
                     !--- 2.5) Wet radius [m]:
                     rwet          = (zwetvol/pi43)**third
                     rwet          = min(rwet, 30.e-6_r8)          ! upper bound based on 1 day lifetime
                     rwet          = max(rwet, rdry)
                     gf            = rwet / rdry                   ! growth factor of a particle
                     !--- 2.6) Assume the gf is the same and we apply it to MAM's dry/wet radius [m]:
                     wetrad(i,k,n) = gf * dryrad(i,k,n)            
                     wetrad(i,k,n) = max(wetrad(i,k,n), dryrad(i,k,n))
                     wetvol(i,k,n) = pi43*wetrad(i,k,n)**3
                     wetvol(i,k,n) = max(wetvol(i,k,n), dryvol(i,k,n))
                     wtrvol(i,k,n) = wetvol(i,k,n) - dryvol(i,k,n)
                     wtrvol(i,k,n) = max(wtrvol(i,k,n), 0.0_r8)
                  end if ! if for relative humidity < crh

               end if    ! if for sulfate + non sea-salt species

               ! apply simple treatment of deliquesence/crystallization hysteresis
               ! for rhcrystal < rh < rhdeliques, aerosol water is a fraction of
               ! the "upper curve" value, and the fraction is a linear function of rh
!               if (rh(i,k) < rhcrystal(n)) then
!                  wetrad(i,k,n) = dryrad(i,k,n)
!                  wetvol(i,k,n) = dryvol(i,k,n)
!                  wtrvol(i,k,n) = 0.0_r8
!               else if (rh(i,k) < rhdeliques(n)) then
!                  if ( aw_l_hysteresis ) then       ! JS added on 08-05-2019
!                     wtrvol(i,k,n) = wtrvol(i,k,n)*hystfac*(rh(i,k) - rhcrystal(n))
!                     wtrvol(i,k,n) = max(wtrvol(i,k,n), 0.0_r8)
!                     wetvol(i,k,n) = dryvol(i,k,n) + wtrvol(i,k,n)
!                     wetrad(i,k,n) = (wetvol(i,k,n)/pi43)**third
!                  end if
!               end if

            end do ! i, ncol loop
         end do    ! k, pver loop

      end do mode_loop  ! n, mode loop

      end subroutine m7_wateruptake_scheme 

!----------------------------------------------------------------------

   end module modal_aero_wateruptake


