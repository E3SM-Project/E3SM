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
#if ( defined MOSAIC_SPECIES )
use time_manager,              only: get_nstep, get_step_size
use module_mosaic_box_aerchem, only: MOSAIC_aerosol_water_only
use module_data_mosaic_aero,   only: naer_mosaic => naer, &
                                     inh4_a, ilim2_a, iso4_a, ina_a, icl_a, ibc_a, ioin_a, ioc_a, imom_a, &
                                     ino3_a, icl_a, ica_a, ico3_a, &
                                     ilim2_g, ih2so4_g, inh3_g, ihno3_g, ihcl_g, &
                                     jhyst_up, jtotal, &
                                     nbin_a, nbin_a_max, ngas_volatile, nmax_astem, nmax_mesa, no_aerosol, nsalt, &
                                     mosaic_vars_aa_type, mw_aer_mac, mw_gas
use constituents,              only: pcnst
use physconst,                 only: r_universal, rair
use modal_aero_data,           only: dgnumlo_amode, numptr_amode, ntot_amode, nsoa, npoa, nbc, &
                                     lptr_so4_a_amode, lptr_no3_a_amode, lptr_cl_a_amode, &
                                     lptr_msa_a_amode, lptr_nh4_a_amode, lptr_nacl_a_amode, &
                                     lptr_dust_a_amode, lptr_ca_a_amode, lptr_co3_a_amode, &
                                     lptr_h2so4_g_amode, lptr_hno3_g_amode, lptr_hcl_g_amode, &
                                     lptr_nh3_g_amode, lptr_soa_a_amode, lptr_bc_a_amode, &
                                     lptr_pom_a_amode, lptr_mom_a_amode, lptr2_soa_g_amode
use modal_aero_amicphys,       only: hygro_bc, hygro_pom, hygro_mom, hygro_soa, hygro_dst
#endif

implicit none
private
save

public :: &
   modal_aero_wateruptake_init, &
   modal_aero_wateruptake_dr

public :: modal_aero_wateruptake_reg

real(r8), parameter :: third = 1._r8/3._r8
real(r8), parameter :: pi43  = pi*4.0_r8/3.0_r8


! Physics buffer indices
integer :: cld_idx        = 0
integer :: dgnum_idx      = 0
integer :: dgnumwet_idx   = 0
integer :: wetdens_ap_idx = 0
integer :: qaerwat_idx    = 0

logical :: pergro_mods         = .false.

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
!$OMP THREADPRIVATE(maer, hygro, naer, dryvol, drymass, dryrad, wetrad, wetvol, wtrvol, rhcrystal, rhdeliques, specdens_1)

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
   use time_manager,   only: is_first_step
   use physics_buffer, only: pbuf_set_field
   use shr_log_mod ,   only: errmsg => shr_log_errmsg

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   integer :: m, nmodes, istat
   logical :: history_aerosol      ! Output the MAM aerosol variables and tendencies
   logical :: history_verbose      ! produce verbose history output

   character(len=3) :: trnum       ! used to hold mode number (as characters)
   !----------------------------------------------------------------------------

   cld_idx        = pbuf_get_index('CLD')    
   dgnum_idx      = pbuf_get_index('DGNUM')    

   ! assume for now that will compute wateruptake for climate list modes only

   call rad_cnst_get_info(0, nmodes=nmodes)

   !$OMP PARALLEL
   allocate(maer(pcols,pver,nmodes),stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate maer:       "//errmsg(__FILE__,__LINE__) )
   allocate(hygro(pcols,pver,nmodes),   stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate hygro:      "//errmsg(__FILE__,__LINE__) )
   allocate(naer(pcols,pver,nmodes),    stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate naer:       "//errmsg(__FILE__,__LINE__) )
   allocate(dryvol(pcols,pver,nmodes),  stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate dryvol:     "//errmsg(__FILE__,__LINE__) )
   allocate(drymass(pcols,pver,nmodes), stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate drymass:    "//errmsg(__FILE__,__LINE__) )
   allocate(dryrad(pcols,pver,nmodes),  stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate dryradr:    "//errmsg(__FILE__,__LINE__) )
   allocate(wetrad(pcols,pver,nmodes),  stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate wetrad:     "//errmsg(__FILE__,__LINE__) )
   allocate(wetvol(pcols,pver,nmodes),  stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate wetvol:     "//errmsg(__FILE__,__LINE__) )
   allocate(wtrvol(pcols,pver,nmodes),  stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate wtrvol:     "//errmsg(__FILE__,__LINE__) )
   allocate(rhcrystal(nmodes),          stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate rhcrystal:  "//errmsg(__FILE__,__LINE__) )
   allocate(rhdeliques(nmodes),         stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate rhdeliques: "//errmsg(__FILE__,__LINE__) )
   allocate(specdens_1(nmodes),         stat=istat)
   if (istat .ne. 0) call endrun("Unable to allocate specdens_1: "//errmsg(__FILE__,__LINE__) )
   !$OMP END PARALLEL


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
                                     qaerwat_m, wetdens_m, clear_rh_in)

  use shr_log_mod ,   only: errmsg => shr_log_errmsg
   !----------------------------------------------------------------------------
   !
   ! CAM specific driver for modal aerosol water uptake code.
   !
   ! *** N.B. *** The calculation has been enabled for diagnostic mode lists
   !              via optional arguments.  If the list_idx arg is present then
   !              all the optional args must be present.
   !
   !----------------------------------------------------------------------------

   ! Arguments
   type(physics_state), target, intent(in)    :: state          ! Physics state variables
   type(physics_buffer_desc),   pointer       :: pbuf(:)        ! physics buffer
   ! Optional inputs for diagnostic mode
   integer,  optional,          intent(in)    :: list_idx_in
   real(r8), optional, allocatable, target, intent(in)    :: dgnumdry_m(:,:,:)
   real(r8), optional, allocatable, target, intent(inout)   :: dgnumwet_m(:,:,:)
   real(r8), optional, allocatable, target, intent(inout)   :: qaerwat_m(:,:,:)
   real(r8), optional, allocatable, target, intent(inout)   :: wetdens_m(:,:,:)
   ! optional input relative humidty (overrides clearsky RH estimate below)
   real(r8), optional,          intent(in)    :: clear_rh_in(pcols,pver)

   !----------------------------------------------------------------------------
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
   logical :: compute_wetdens
#if ( defined MOSAIC_SPECIES )
   integer, parameter :: bigint = huge(1)
   integer :: ierr, lq, nstep
   integer, dimension(nbin_a_max) :: jaerosolstate     !Aerosol state (solid, liquid, gas)
   integer, dimension(nbin_a_max) :: jhyst_leg

   real(r8), parameter ::  oneatminv = 1.0_r8/1.01325e5_r8
   real(r8), pointer :: q(:,:,:)
   real(r8), dimension(ngas_volatile) :: gas
   real(r8), dimension(naer_mosaic) :: kappa_nonelectro
   real(r8), dimension(naer_mosaic,3,nbin_a_max) :: aer
   real(r8) :: mwuse_soa(nsoa), mwuse_poa(npoa)
   real(r8) :: RH_pc             !relative humidity (%)
   real(r8) :: P_atm             !Pressure in atm units
   real(r8) :: airdens
   real(r8) :: aircon            !air molar density (kmol/m3)
   real(r8) :: cair_mol_m3       !Air molar density (mol/m3)
   real(r8) :: dtime
   real(r8) :: factaermass
   real(r8) :: factaernumb
   real(r8) :: factaerwatr
   real(r8) :: factgasmass
   real(r8) :: tmpa, tmp_rh, tmp_hy
   real(r8), dimension(nbin_a_max) :: tmp_num_a
   real(r8), dimension(nbin_a_max) :: dens_dry_a
   real(r8), dimension(nbin_a_max) :: mass_dry_a
   real(r8), dimension(nbin_a_max) :: num_a, Dp_dry_a
   real(r8), dimension(nbin_a_max) :: Dp_wet_a
   real(r8), dimension(nbin_a_max) :: water_a

   type (mosaic_vars_aa_type)  :: mosaic_vars_aa
#endif

   character(len=3) :: trnum       ! used to hold mode number (as characters)
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol = state%ncol

   ! determine default variables
   call phys_getopts(history_aerosol_out = history_aerosol, &
                     history_verbose_out = history_verbose)

   list_idx = 0
   if (present(list_idx_in)) then
      list_idx = list_idx_in

      ! check that all optional args for diagnostic mode are present
      if (.not. present(dgnumdry_m) .or. .not. present(dgnumwet_m) .or. &
          .not. present(qaerwat_m)) then
         call endrun('modal_aero_wateruptake_dr called '// &
              'with list_idx_in but required args not present '//errmsg(__FILE__,__LINE__))
      end if

      ! arrays for diagnostic calculations must be allocated
      if (.not. allocated(dgnumdry_m) .or. .not. allocated(dgnumwet_m) .or. &
          .not. allocated(qaerwat_m)) then
         call endrun('modal_aero_wateruptake_dr called '// &
              'with list_idx_in but required args not allocated '//errmsg(__FILE__,__LINE__))
      end if
   end if ! if present(list_idx_in)

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   !initialize to an invalid value
   naer(:,:,:)    = huge(1.0_r8)
   dryvol(:,:,:)  = huge(1.0_r8)
   drymass(:,:,:) = huge(1.0_r8)
   dryrad(:,:,:)  = huge(1.0_r8)
   wetrad(:,:,:)  = huge(1.0_r8)
   wetvol(:,:,:)  = huge(1.0_r8)
   wtrvol(:,:,:)  = huge(1.0_r8)

   rhcrystal(:)   = huge(1.0_r8)
   rhdeliques(:)  = huge(1.0_r8)
   specdens_1(:)  = huge(1.0_r8)

   maer(:,:,:)     = 0._r8
   hygro(:,:,:)    = 0._r8

   !by default set compute_wetdens to be true
   compute_wetdens = .true.
   if (.not. present(list_idx_in)) then
      call pbuf_get_field(pbuf, dgnum_idx,      dgncur_a )
      call pbuf_get_field(pbuf, dgnumwet_idx,   dgncur_awet )
      call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens)
      call pbuf_get_field(pbuf, qaerwat_idx,    qaerwat)
   else
      dgncur_a    => dgnumdry_m
      dgncur_awet => dgnumwet_m
      qaerwat     => qaerwat_m
      if(present(wetdens_m)) then
         if (.not. allocated(wetdens_m)) then
            call endrun('modal_aero_wateruptake_dr called '// &
                 'with list_idx_in but wetdens_m is not allocated '//errmsg(__FILE__,__LINE__))
         endif
         wetdens     => wetdens_m
      else
         !set compute_wetdens to flase if wetdens is not present
         compute_wetdens = .false.
      endif
   end if

   !----------------------------------------------------------------------------
   ! retreive aerosol properties

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
            end do ! i = 1, ncol
         end do ! k = top_lev, pver
      end do ! l = 1, nspec

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

         end do ! i = 1, ncol
      end do ! k = top_lev, pver

   end do    ! modes

   !----------------------------------------------------------------------------
   ! specify clear air relative humidity

   if (present(clear_rh_in)) then

      ! use input relative humidity
      rh(1:ncol,1:pver) = clear_rh_in(1:ncol,1:pver)

      ! check that values are reasonable and apply upper limit
      do k = top_lev, pver
         do i = 1, ncol
            if ( rh(i,k)<0 ) then
               write(iulog,*) 'modal_aero_wateruptake_dr: clear_rh_in is negative - rh:',rh(i,k),' k=',k
               call endrun('modal_aero_wateruptake_dr: clear_rh_in cannot be negative')
            end if
            ! limit RH to 98% to be consistent with behavior when clear_rh_in is not provided
            rh(i,k) = min(rh(i,k), 0.98_r8)
         end do ! i
      end do ! k

   else

      ! estimate clear air relative humidity using cloud fraction
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
            if (cldn(i,k) .lt. cldn_thresh) then
               rh(i,k) = (rh(i,k) - cldn(i,k)) / (1.0_r8 - cldn(i,k))  ! RH of clear portion
            end if
            rh(i,k) = max(rh(i,k), 0.0_r8)
         end do ! i = 1, ncol
      end do ! k = top_lev, pver

   end if ! if present(clear_rh_in)


   !----------------------------------------------------------------------------
   ! compute aerosol wet radius and aerosol water

#if ( defined MOSAIC_SPECIES )
   ! call mam_amicphys to get aerosol water
   q => state%q

   ! allocate( gas(ngas_volatile) )
   ! allocate the allocatable parts of mosaic_vars_aa
   allocate( mosaic_vars_aa%iter_mesa(nbin_a_max), stat=ierr )
   if (ierr /= 0) call endrun( &
      '*** subr mosaic_gasaerexch_1subarea_intr - allocate error for mosaic_vars_aa%iter_mesa' )

   ! set kappa values for non-electrolyte species
   ! reason for doing this here is that if cam eventually has multiple varieties of dust and/or pom, 
   ! then the dust hygroscopicity may vary spatially and temporally,
   ! and the kappa values cannot be constants
   kappa_nonelectro(:)       = 0.0_r8
   kappa_nonelectro(ibc_a  ) = hygro_bc
   kappa_nonelectro(ioc_a  ) = hygro_pom
   kappa_nonelectro(imom_a ) = hygro_mom
   kappa_nonelectro(ilim2_a) = hygro_soa
   kappa_nonelectro(ioin_a ) = hygro_dst

   dtime = get_step_size()
   nstep = get_nstep()

   do k = top_lev, pver
      do i = 1, ncol
         RH_pc = 100._r8*rh(i,k)
         P_atm = pmid(i,k) * oneatminv      ! Pressure (atm)
         airdens = pmid(i,k)/(rair*t(i,k))
         factaermass = airdens * 1.0e12_r8  ! converts (kg-aer/kg-air) to (ng-aer/m3) 
         factaernumb = airdens * 1.0e-6_r8  ! converts (#/kg-air) to (#/cm3)
         factaerwatr = airdens              ! converts (kg-h2o/kg-air) to (kg-h2o/m3) 
         factgasmass = airdens * 1.0e12_r8  ! converts (kg-gas/kg-air) to (ng-gas/m3) 

         ! Populate aersols
         aer(:,:,:) = 0.0_r8 ! initialized to zero at every grid point for safety
         num_a(:)   = 0.0_r8 ! initialized to zero
         water_a(:) = 0.0_r8 ! initialized to zero
         
         do m = 1, ntot_amode
            ! Notes:
            ! 1. NCL(sea salt) of CAM is mapped in NA and CL of MOSAIC
            ! 2. SOA of CAM is lumped into LIM2 species of MOSAIC !BALLI *ASK RAHUL and Dick
            ! 3. Species NO3, MSA, CO3, Ca do not exist in CAM therefore not mapped here
            ! 4. Species ARO1, ARO2, ALK1, OLE1, API1, API2, LIM1 are SOA species in MOSAIC
            !    which are not used in CAM-MOSAIC framework as of now
            ! 5. CAM units are (kg/kg of air) which are converted to Mosaic units (nano mol/m3).

            ! Units conversion:  aer [nmol/m3] = q [kg/kg-air] * factaermass /
            ! specmw
            lq = lptr_nh4_a_amode(m)
            if (1 <= lq .and. lq <= pcnst) &
               aer(inh4_a,  jtotal, m) = q(i,k,lq)*factaermass/mw_aer_mac(inh4_a)

            lq = lptr_soa_a_amode(m)
            if (1 <= lq .and. lq <= pcnst) &
               aer(ilim2_a, jtotal, m) = q(i,k,lq)*factaermass/mw_aer_mac(ilim2_a)

            lq = lptr_so4_a_amode(m)
            if (1 <= lq .and. lq <= pcnst) &
               aer(iso4_a,  jtotal, m) = q(i,k,lq)*factaermass/mw_aer_mac(iso4_a)

            lq = lptr_nacl_a_amode(m)
            if (1 <= lq .and. lq <= pcnst) &
               aer(ina_a,   jtotal, m) = q(i,k,lq)*factaermass/mw_aer_mac(ina_a)

            lq = lptr_cl_a_amode(m)
            if (1 <= lq .and. lq <= pcnst) &
               aer(icl_a,   jtotal, m) = q(i,k,lq)*factaermass/mw_aer_mac(icl_a)

            lq = lptr_no3_a_amode(m)
            if (1 <= lq .and. lq <= pcnst) &
               aer(ino3_a,  jtotal, m) = q(i,k,lq)*factaermass/mw_aer_mac(ino3_a)

            lq = lptr_ca_a_amode(m)
            if (1 <= lq .and. lq <= pcnst) &
               aer(ica_a,   jtotal, m) = q(i,k,lq)*factaermass/mw_aer_mac(ica_a)

            lq = lptr_co3_a_amode(m)
            if (1 <= lq .and. lq <= pcnst) &
               aer(ico3_a,  jtotal, m) = q(i,k,lq)*factaermass/mw_aer_mac(ico3_a)

            ! Units conversion:  aer [ng/m3] = q [kg/kg-air] * factaermass
            lq = lptr_bc_a_amode(m)
            if (1 <= lq .and. lq <= pcnst) &
               aer(ibc_a,   jtotal, m) = q(i,k,lq)*factaermass

            lq = lptr_pom_a_amode(m)
            if (1 <= lq .and. lq <= pcnst) &
               aer(ioc_a,   jtotal, m) = q(i,k,lq)*factaermass

            lq = lptr_mom_a_amode(m)
            if (1 <= lq .and. lq <= pcnst) &
               aer(imom_a,  jtotal, m) = q(i,k,lq)*factaermass

            lq = lptr_dust_a_amode(m)
            if (1 <= lq .and. lq <= pcnst) &
               aer(ioin_a,  jtotal, m) = q(i,k,lq)*factaermass
             
            ! Populate aerosol number and water species
            ! Units conversion:  num_a [#/cm3] = q [#/kg-air] * factaernumb
            ! lq = numptr_amode(m)
            ! if (1 <= lq .and. lq <= pcnst) &
            !    num_a(m) = q(i,k,lq)*factaernumb

            ! steve -- try using this bounded number concentration to get a reasonable size in MOSAIC
            num_a(m) = naer(i,k,m)*factaernumb
            tmp_num_a(m) = num_a(m)

            ! (1/v2ncur_a) = volume of 1 particle whose size (dgn) is bounded
            ! naer = dryvolmr/v2ncur_a = number mixing ratio of size-bounded particles
            !
            ! when q(i,k,numptr_amode(n)) is much too small for some reason, 
            ! so that the unbounded particle size is much too big,
            ! and num_a(m) is set from q(i,k,numptr_amode(n)), then
            !   > the dp_dry_a and dp_wet_a from MOSAIC will be much too big
            !   > qaerwat calculated from naer and [(dp_wet^3 - dp_dry^3)*pi/6] will be much too big
            !   > the water_a from MOSAIC may be OK, and qaerwat calculated from water_a may be OK
            !
            ! if we want to calculate qaerwat from naer and [(dp_wet^3 - dp_dry^3)*pi/6],
            !    similar to the kohler approach, then it follows that num_a should be set using naer
            ! this is also need to get reasonable dgnumwet values
         end do ! m

         ! Populate gases
         gas(:) = 0.0_r8 ! Initialized to zero

         ! Units conversion:  gas [nmol/m3] = q [kg-gas/kg-air] * factgasmass / specmw
         lq = lptr_h2so4_g_amode
         if (1 <= lq .and. lq <= pcnst) &
            gas(ih2so4_g) = q(i,k,lq) * factaermass / mw_gas(ih2so4_g)

         lq = lptr_hno3_g_amode
         if (1 <= lq .and. lq <= pcnst) &
            gas(ihno3_g) = q(i,k,lq) * factaermass / mw_gas(ihno3_g)

         lq = lptr_hcl_g_amode
         if (1 <= lq .and. lq <= pcnst) &
            gas(ihcl_g) = q(i,k,lq) * factaermass / mw_gas(ihcl_g)

         lq = lptr_nh3_g_amode
         if (1 <= lq .and. lq <= pcnst) &
            gas(inh3_g) = q(i,k,lq) * factaermass / mw_gas(inh3_g)

         lq = lptr2_soa_g_amode(1)
         if (1 <= lq .and. lq <= pcnst) &
            gas(ilim2_g) = q(i,k,lq) * factaermass / mw_gas(ilim2_g)

         ! initialize mosaic_vars_aa
         mosaic_vars_aa%jastem_fail = 0
         mosaic_vars_aa%it_mosaic = nstep
         mosaic_vars_aa%iter_mesa(1:nbin_a_max) = 0
         mosaic_vars_aa%zero_water_flag = .false.
         mosaic_vars_aa%flag_itr_kel = .false.
         mosaic_vars_aa%hostgridinfo(1)   = i
         mosaic_vars_aa%hostgridinfo(2)   = k
         mosaic_vars_aa%hostgridinfo(3)   = lchnk
         mosaic_vars_aa%hostgridinfo(4:6) = bigint
         mosaic_vars_aa%it_host = 0
         mosaic_vars_aa%f_mos_fail = -1
         mosaic_vars_aa%isteps_astem = 0
         mosaic_vars_aa%isteps_astem_max = 0
         mosaic_vars_aa%jastem_call = 0
         mosaic_vars_aa%jmesa_call = 0
         mosaic_vars_aa%jmesa_fail = 0
         mosaic_vars_aa%niter_mesa_max = 0
         mosaic_vars_aa%nmax_astem = nmax_astem
         mosaic_vars_aa%nmax_mesa = nmax_mesa
         mosaic_vars_aa%cumul_steps_astem = 0.0_r8
         mosaic_vars_aa%niter_mesa = 0.0_r8
         mosaic_vars_aa%xnerr_astem_negative(:,:) = 0.0_r8
         mosaic_vars_aa%fix_astem_negative = 0

         ! these initializations may not be needed but are done for safety
         jaerosolstate(:) = no_aerosol
         jhyst_leg(:) = jhyst_up
         Dp_dry_a(:) = 0.0_r8
         Dp_wet_a(:) = 0.0_r8
         mass_dry_a(:) = 0.0_r8
         dens_dry_a(:) = 0.0_r8

         call MOSAIC_aerosol_water_only(rh(i,k),  t(i,k), &  ! Intent-ins
            P_atm,             RH_pc,      dtime,         &
            kappa_nonelectro,                             &
            jaerosolstate,     jhyst_leg,                 &  ! Intent-inouts
            aer,               num_a,      water_a,  gas, &
            Dp_dry_a,          Dp_wet_a,                  &
            mosaic_vars_aa,                               &
            mass_dry_a,        dens_dry_a)                   ! Intent-outs:489

         do m = 1, ntot_amode
            if (jaerosolstate(m) == no_aerosol) then
               ! mosaic did not calculate water because aerosol mass mixing ratio is very small 
               wetrad(i,k,m) = dryrad(i,k,m)
            else
               dryrad(i,k,m) = 0.5e-2_r8*Dp_dry_a(m)   ! convert from cm to m
               wetrad(i,k,m) = max( 0.5e-2_r8*Dp_wet_a(m), dryrad(i,k,m) )
            endif

            if ( (wetrad(i,k,m) > dryrad(i,k,m)*10.0_r8   ) .or. &
                 (dryrad(i,k,m) < dgnumlo_amode(m)*0.05_r8) ) then
               ! wetrad unreasonably larger than dryrad, or dryrad unreasonably small
               ! (may want to deactivate this error message in long runs)
               write(iulog,'(2a,2i10,4i5,1p,e17.9,7e12.4)') 'm_a_wateruptake error - ', &
                  'nstep,lchnk,i,k,m,state, RH_pc,hygro,dryrad,wetrad=', &
                  nstep, lchnk, i, k, m, jaerosolstate(m), &
                  RH_pc, hygro(i,k,m), dryrad(i,k,m), wetrad(i,k,m), &
                  naer(i,k,m), maer(i,k,m), drymass(i,k,m), dryvol(i,k,m)

               dryrad(i,k,m) = max( dryrad(i,k,m), dgnumlo_amode(m)*0.05_r8 )
               if (rh(i,k) <= rhcrystal(m)) then
                  wetrad(i,k,m) = dryrad(i,k,m)  ! no water
               else
                  ! use the following which ignores kelvin effect
                  ! (1/rh) = 1 + hygro*dryvol/wtrvol
                  ! wtrvol = dryvol*[hygro/((1/rh) - 1)]
                  ! wetvol = dryvol + wtrvol = dryvol*[1 + hygro/((1/rh) - 1)]
                  tmp_hy = max( 0.00_r8,    min( 2.00_r8, hygro(i,k,m) ) )
                  tmp_rh = max( 1.0e-20_r8, min( 0.99_r8, rh(i,k)    ) )
                  tmpa = ( 1.0_r8 + tmp_hy/((1.0_r8/tmp_rh) - 1.0_r8) )**third
                  wetrad(i,k,m) = dryrad(i,k,m)*max( tmpa, 1.0_r8 )
               endif
            endif 
            
            dryvol(i,k,m) = pi43*dryrad(i,k,m)**3
            wetvol(i,k,m) = pi43*wetrad(i,k,m)**3
            wtrvol(i,k,m) = max( wetvol(i,k,m) - dryvol(i,k,m), 0.0_r8 )
            ! if (dryrad(i,k,m).le.1.e-20) then
            !    write(iulog,*)'mode,dryrad,wetrad=',m,dryrad(i,k,m),wetrad(i,k,m)
            !    if(lptr_pom_a_amode(m)>=1) write(iulog,*) 'q for pom  =', q(i, k, lptr_pom_a_amode(m))
            !    if(lptr_bc_a_amode(m)>=1)  write(iulog,*) 'q for bc   =', q(i, k, lptr_bc_a_amode(m))
            !    if(lptr_so4_a_amode(m)>=1) write(iulog,*) 'q for so4  =', q(i, k, lptr_so4_a_amode(m))
            !    if(lptr_soa_a_amode(m)>=1) write(iulog,*) 'q for soa  =', q(i, k, lptr_soa_a_amode(m))
            !    if(lptr_nacl_a_amode(m)>=1)write(iulog,*) 'q for nacl =', q(i, k, lptr_nacl_a_amode(m))
            !    write(iulog,*) ’num=‘, num_a(m))
            !    call endrun('dryrad too small')
            ! endif
            ! if ( (m.eq.6) .and. (k.eq.pver) .and. (RH_pc>95.) .and. (dryvol(i,k,m).gt.1.e-30)) then
            !    if (wtrvol(i,k,m)/dryvol(i,k,m).lt.1.)then
            !       write(iulog,*) 'wtrvol/dryvol = ', wtrvol(i,k,m)/dryvol(i,k,m)
            !    endif 
            !    if (wtrvol(i,k,m)/dryvol(i,k,m).gt.1.e10)then
            !       write(iulog,*) 'RH_pc,wtrvol/dryvol = ', RH_pc, wtrvol(i,k,m)/dryvol(i,k,m)
            !       write(iulog,*) 'lchnk,i,k,m = ',lchnk,i,k,m
            !       call endrun('excessive aerosol water')
            !    endif
            ! endif
            
            naer(i,k,m) = num_a(m)/factaernumb  ! do this in case MOSAIC adjusted num_a
                                                ! this will not affect state%q(:,:,numptr_amode(:))
         end do ! m
      end do ! i
   end do ! k
#else
   call modal_aero_wateruptake_sub( &
      ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
      hygro, rh, dryvol, wetrad, wetvol,           &
      wtrvol)
#endif

   do m = 1, nmodes

      do k = top_lev, pver
         do i = 1, ncol

            dgncur_awet(i,k,m) = dgncur_a(i,k,m) * (wetrad(i,k,m)/dryrad(i,k,m))
            qaerwat(i,k,m)     = rhoh2o*naer(i,k,m)*wtrvol(i,k,m)

            ! compute aerosol wet density (kg/m3)
            if(compute_wetdens) then
               if (wetvol(i,k,m) > 1.0e-30_r8) then
                  wetdens(i,k,m) = (drymass(i,k,m) + rhoh2o*wtrvol(i,k,m))/wetvol(i,k,m)
               else
                  wetdens(i,k,m) = specdens_1(m)
               end if
            endif
         end do ! i = 1, ncol
      end do ! k = top_lev, pver

   end do ! m = 1, nmodes

   !----------------------------------------------------------------------------
   ! write history output if not in diagnostic mode

   if (.not.present(list_idx_in)) then

      aerosol_water(:ncol,:) = 0._r8
      do m = 1, nmodes
         ! output to history
         write( trnum, '(i3.3)' ) m
         call outfld( 'wat_a'//trnum(3:3),  qaerwat(:,:,m),     pcols, lchnk)
         call outfld( 'dgnd_a'//trnum(2:3), dgncur_a(:,:,m),    pcols, lchnk)
         call outfld( 'dgnw_a'//trnum(2:3), dgncur_awet(:,:,m), pcols, lchnk)
         if (history_aerosol .and. .not. history_verbose) &
         aerosol_water(:ncol,:) = aerosol_water(:ncol,:) + qaerwat(:ncol,:,m)
      end do ! m = 1, nmodes

      if (history_aerosol .and. .not. history_verbose) &
         call outfld( 'aero_water',  aerosol_water(:ncol,:),    ncol, lchnk)

   end if
#if ( defined MOSAIC_SPECIES )
   ! deallocate the allocatable parts of mosaic_vars_aa
   deallocate( mosaic_vars_aa%iter_mesa, stat=ierr )
   if (ierr /= 0) call endrun( &
      '*** subr modal_aero_wateruptake_dr - deallocate error for mosaic_vars_aa%iter_mesa' )
#endif

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

!-----------------------------------------------------------------------
      subroutine modal_aero_kohler(   &
          rdry_in, hygro, s, rwet_out, im )

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
      real(r8), parameter :: surften = 76._r8
      real(r8), parameter :: tair = 273._r8
      real(r8), parameter :: third = 1._r8/3._r8
      real(r8), parameter :: ugascon = 8.3e7_r8


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

        if(s(i).gt.1._r8-eps)then
!          save quartic solution at s=1-eps
           r4=r(i)
!          cubic for rh=1
           p=abs(p31(i))/(rdry(i)*rdry(i))
           if(p.lt.eps)then
              r(i)=rdry(i)*(1._r8+p*third)
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

   end module modal_aero_wateruptake


