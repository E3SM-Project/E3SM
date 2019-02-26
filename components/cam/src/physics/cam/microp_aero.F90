
module microp_aero

!---------------------------------------------------------------------------------
! Purpose:
!   CAM driver layer for aerosol activation processes.
!
! ***N.B.*** This module is currently hardcoded to recognize only the aerosols/modes that
!            affect the climate calculation.  This is implemented by using list
!            index 0 in all the calls to rad_constituent interfaces.
!
! Author: Andrew Gettelman
! Based on code from: Hugh Morrison, Xiaohong Liu and Steve Ghan
! May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)         
! for questions contact Andrew Gettelman  (andrew@ucar.edu)
! Modifications: A. Gettelman Nov 2010  - changed to support separation of 
!                  microphysics and macrophysics and concentrate aerosol information here
!                B. Eaton, Sep 2014 - Refactored to move CAM interface code into the CAM
!                  interface modules and preserve just the driver layer functionality here.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8=>shr_kind_r8
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols, pver, pverp
use ref_pres,         only: top_lev => trop_cloud_top_lev
use physconst,        only: rair, gravit, pi
use constituents,     only: cnst_get_ind
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,     only: phys_getopts, cam_chempkg_is, use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num

use nucleate_ice_cam, only: use_preexisting_ice, nucleate_ice_cam_readnl, nucleate_ice_cam_register, &
                            nucleate_ice_cam_init, nucleate_ice_cam_calc

use ndrop,            only: ndrop_init, dropmixnuc
use ndrop_bam,        only: ndrop_bam_init, ndrop_bam_run, ndrop_bam_ccn

use hetfrz_classnuc_cam, only: hetfrz_classnuc_cam_readnl, hetfrz_classnuc_cam_register, hetfrz_classnuc_cam_init, &
                               hetfrz_classnuc_cam_save_cbaero, hetfrz_classnuc_cam_calc

use cam_history,      only: addfld, add_default, outfld
use cam_logfile,      only: iulog
use cam_abortutils,       only: endrun
use perf_mod,         only: t_startf, t_stopf

implicit none
private
save

public :: microp_aero_init, microp_aero_run, microp_aero_readnl, microp_aero_register

! Private module data

character(len=16)   :: eddy_scheme
logical             :: micro_do_icesupersat

!!   icenul_wsub_scheme = 1 : f(TKE) as default
!!                        2 : Mean updraft calculated from Gausssian PDF, with stddev=f(TKE)
integer             :: icenul_wsub_scheme = 1
 
! contact freezing due to dust
! dust number mean radius (m), Zender et al JGR 2003 assuming number mode radius of 0.6 micron, sigma=2
real(r8), parameter :: rn_dst1 = 0.258e-6_r8
real(r8), parameter :: rn_dst2 = 0.717e-6_r8
real(r8), parameter :: rn_dst3 = 1.576e-6_r8
real(r8), parameter :: rn_dst4 = 3.026e-6_r8

real(r8) :: bulk_scale    ! prescribed aerosol bulk sulfur scale factor

! smallest mixing ratio considered in microphysics
real(r8), parameter :: qsmall = 1.e-18_r8

! minimum allowed cloud fraction
real(r8), parameter :: mincld = 0.0001_r8

! indices in state%q and pbuf structures
integer :: cldliq_idx = -1
integer :: cldice_idx = -1
integer :: numliq_idx = -1
integer :: numice_idx = -1
integer :: kvh_idx = -1
integer :: tke_idx = -1
integer :: wp2_idx = -1
integer :: ast_idx = -1
integer :: alst_idx = -1
integer :: aist_idx = -1

integer :: cldo_idx = -1
integer :: dgnumwet_idx = -1

! Bulk aerosols
character(len=20), allocatable :: aername(:)
real(r8), allocatable :: num_to_mass_aer(:)

integer :: naer_all      ! number of aerosols affecting climate
integer :: idxsul   = -1 ! index in aerosol list for sulfate
integer :: idxdst2  = -1 ! index in aerosol list for dust2
integer :: idxdst3  = -1 ! index in aerosol list for dust3
integer :: idxdst4  = -1 ! index in aerosol list for dust4

! modal aerosols
logical :: clim_modal_aero

integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
integer :: mode_coarse_idx = -1  ! index of coarse mode
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode
integer :: mode_fine_dst_idx = -1   ! index of dust in fine dust mode
integer :: mode_pcarbon_idx  = -1  ! index of dust in accum mode
integer :: accum_dust_idx    = -1  ! index of dust in accum mode
logical :: dem_in            = .false.           ! use DeMott IN

integer :: npccn_idx, rndst_idx, nacon_idx

logical  :: separate_dust = .false.
logical  :: liqcf_fix
real(r8), parameter :: unset_r8   = huge(1.0_r8)
real(r8) :: wsubmin = unset_r8 !PMA sets a much lower lower bound


contains
!=========================================================================================

subroutine microp_aero_register
   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Register pbuf fields for aerosols needed by microphysics
   ! 
   ! Author: Cheryl Craig October 2012
   ! 
   !-----------------------------------------------------------------------
   use ppgrid,         only: pcols
   use physics_buffer, only: pbuf_add_field, dtype_r8

   call pbuf_add_field('NPCCN',      'physpkg',dtype_r8,(/pcols,pver/), npccn_idx)

   call pbuf_add_field('RNDST',      'physpkg',dtype_r8,(/pcols,pver,4/), rndst_idx)
   call pbuf_add_field('NACON',      'physpkg',dtype_r8,(/pcols,pver,4/), nacon_idx)
 
   call nucleate_ice_cam_register()
   call hetfrz_classnuc_cam_register()

end subroutine microp_aero_register

!=========================================================================================

subroutine microp_aero_init

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Initialize constants for aerosols needed by microphysics
   ! 
   ! Author: Andrew Gettelman May 2010
   ! 
   !-----------------------------------------------------------------------

   ! local variables
   integer  :: iaer, ierr
   integer  :: m, n, nmodes, nspec

   character(len=32) :: str32
   character(len=*), parameter :: routine = 'microp_aero_init'
   logical :: history_amwg
   !-----------------------------------------------------------------------

   ! Query the PBL eddy scheme
   call phys_getopts(eddy_scheme_out          = eddy_scheme, &
        history_amwg_out = history_amwg, &
        micro_do_icesupersat_out = micro_do_icesupersat, &
        liqcf_fix_out    = liqcf_fix,    & 
        demott_ice_nuc_out = dem_in      ) 
   
   if(masterproc)write(iulog,*)'DEMOTT is:', dem_in 

   ! Access the physical properties of the aerosols that are affecting the climate
   ! by using routines from the rad_constituents module.

   ! get indices into state and pbuf structures
   call cnst_get_ind('CLDLIQ', cldliq_idx)
   call cnst_get_ind('CLDICE', cldice_idx)
   call cnst_get_ind('NUMLIQ', numliq_idx)
   call cnst_get_ind('NUMICE', numice_idx)

   select case(trim(eddy_scheme))
   case ('diag_TKE', 'SHOC_SGS')
      tke_idx      = pbuf_get_index('tke')   
   case ('CLUBB_SGS')
      wp2_idx = pbuf_get_index('WP2_nadv')
   case default
      kvh_idx      = pbuf_get_index('kvh')
   end select

   ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
   ! The modal aerosols can be either prognostic or prescribed.
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   ast_idx      = pbuf_get_index('AST')
   if(liqcf_fix) then
      alst_idx      = pbuf_get_index('ALST')
      aist_idx      = pbuf_get_index('AIST')
   endif

   if (clim_modal_aero) then

      cldo_idx     = pbuf_get_index('CLDO')
      dgnumwet_idx = pbuf_get_index('DGNUMWET')

      call ndrop_init()

      ! Init indices for specific modes/species

      ! mode index for specified mode types
      do m = 1, nmodes
         call rad_cnst_get_info(0, m, mode_type=str32)
         select case (trim(str32))
         case ('accum')
            mode_accum_idx = m
         case ('aitken')
            mode_aitken_idx = m
         case ('coarse')
            mode_coarse_idx = m
         case ('coarse_dust')
            mode_coarse_dst_idx = m
         case ('coarse_seasalt')
            mode_coarse_slt_idx = m
         case ('fine_dust')
            mode_fine_dst_idx = m
         case ('primary_carbon')
            mode_pcarbon_idx  = m            
         end select
      end do

      ! check if coarse dust is in separate mode
      separate_dust = mode_coarse_dst_idx > 0

      ! for 3-mode 
      if ( mode_coarse_dst_idx<0 ) mode_coarse_dst_idx = mode_coarse_idx
      if ( mode_coarse_slt_idx<0 ) mode_coarse_slt_idx = mode_coarse_idx

      ! Check that required mode types were found
      if (mode_accum_idx == -1 .or. mode_aitken_idx == -1 .or. &
          mode_coarse_dst_idx == -1.or. mode_coarse_slt_idx == -1) then
         write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
            mode_accum_idx, mode_aitken_idx, mode_coarse_dst_idx, mode_coarse_slt_idx
         call endrun(routine//': ERROR required mode type not found')
      end if

      ! species indices for specified types
      ! find indices for the dust and seasalt species in the coarse mode
      call rad_cnst_get_info(0, mode_coarse_dst_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_dst_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('dust')
            coarse_dust_idx = n
         end select
      end do
      call rad_cnst_get_info(0, mode_coarse_slt_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_slt_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('seasalt')
            coarse_nacl_idx = n
         end select
      end do
      
      if(dem_in) then
         call rad_cnst_get_info(0, mode_accum_idx, nspec=nspec)
         do n = 1, nspec
            call rad_cnst_get_info(0, mode_accum_idx, n, spec_type=str32)
            select case (trim(str32))
            case ('dust')
               accum_dust_idx = n
            end select
         end do
      endif

      ! Check that required mode specie types were found
      if ( coarse_dust_idx == -1 .or. coarse_nacl_idx == -1) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_dust_idx, coarse_nacl_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if

   else

      ! Props needed for BAM number concentration calcs.

      call rad_cnst_get_info(0, naero=naer_all)
      allocate( &
         aername(naer_all),        &
         num_to_mass_aer(naer_all) )

      do iaer = 1, naer_all
         call rad_cnst_get_aer_props(0, iaer, &
            aername         = aername(iaer), &
            num_to_mass_aer = num_to_mass_aer(iaer) )

         ! Look for sulfate, dust, and soot in this list (Bulk aerosol only)
         if (trim(aername(iaer)) == 'SULFATE') idxsul = iaer
         if (trim(aername(iaer)) == 'DUST2') idxdst2 = iaer
         if (trim(aername(iaer)) == 'DUST3') idxdst3 = iaer
         if (trim(aername(iaer)) == 'DUST4') idxdst4 = iaer
      end do

      call ndrop_bam_init()

   end if

   call addfld('LCLOUD', (/ 'lev' /), 'A', ' ', 'Liquid cloud fraction used in stratus activation')

   call addfld('WSUB',  (/ 'lev' /), 'A', 'm/s', 'Diagnostic sub-grid vertical velocity'                   )
   call addfld('WSUBI', (/ 'lev' /), 'A', 'm/s', 'Diagnostic sub-grid vertical velocity for ice'           )

   call addfld('WLARGE',(/ 'lev' /), 'A', 'm/s', 'Large-scale vertical velocity'                           )
   call addfld('WSIG',  (/ 'lev' /), 'A', 'm/s', 'Subgrid standard deviation of vertical velocity'         )
   call addfld('WSUBI2',(/ 'lev' /), 'A', 'm/s', 'Mean updraft, with stddev=f(TKE)'                        )
   call addfld('RHICE', (/ 'lev' /), 'A', '0-1', 'RHi for ice nucleation'                                  )

   if (history_amwg) then
      call add_default ('WSUB     ', 1, ' ')
   end if

   call nucleate_ice_cam_init(mincld, bulk_scale)
   call hetfrz_classnuc_cam_init(mincld)

end subroutine microp_aero_init

!=========================================================================================

subroutine microp_aero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Namelist variables
   real(r8) :: microp_aero_bulk_scale = 2._r8  ! prescribed aerosol bulk sulfur scale factor
   real(r8) :: microp_aero_wsubmin  = 0.2_r8
   integer  :: microp_aero_wsub_scheme = 1     ! updraft velocity parameterization option for ice nucleation
 
   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'microp_aero_readnl'

   namelist /microp_aero_nl/ microp_aero_bulk_scale, microp_aero_wsub_scheme,microp_aero_wsubmin

   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'microp_aero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, microp_aero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variable
   call mpibcast(microp_aero_bulk_scale, 1, mpir8, 0, mpicom)
   call mpibcast(microp_aero_wsub_scheme, 1, mpiint, 0, mpicom)
   call mpibcast(microp_aero_wsubmin,     1, mpir8, 0, mpicom)
#endif

   ! set local variables
   bulk_scale = microp_aero_bulk_scale
   icenul_wsub_scheme = microp_aero_wsub_scheme
   wsubmin = microp_aero_wsubmin

   call nucleate_ice_cam_readnl(nlfile)
   call hetfrz_classnuc_cam_readnl(nlfile)

end subroutine microp_aero_readnl

!=========================================================================================

subroutine microp_aero_run ( &
   state, ptend, deltatin, pbuf, liqcldfo )

   ! input arguments
   type(physics_state), target, intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend
   real(r8),                    intent(in)    :: deltatin     ! time step (s)
   real(r8),                    intent(in)    :: liqcldfo(pcols,pver)  ! old liquid cloud fraction
   type(physics_buffer_desc),   pointer       :: pbuf(:)

   ! local workspace
   ! all units mks unless otherwise stated

   integer :: i, k, m
   integer :: itim_old
   integer :: nmodes
   real(r8):: dst1_num_to_mass 

   real(r8), pointer :: ast(:,:)        
   real(r8), pointer :: alst(:,:)        
   real(r8), pointer :: aist(:,:)        

   real(r8), pointer :: npccn(:,:)      ! number of CCN (liquid activated)

   real(r8), pointer :: rndst(:,:,:)    ! radius of 4 dust bins for contact freezing
   real(r8), pointer :: nacon(:,:,:)    ! number in 4 dust bins for contact freezing

   real(r8), pointer :: num_coarse(:,:) ! number m.r. of coarse mode
   real(r8), pointer :: coarse_dust(:,:) ! mass m.r. of coarse dust
   real(r8), pointer :: coarse_nacl(:,:) ! mass m.r. of coarse nacl

   real(r8), pointer :: accum_dust(:,:) ! mass m.r. of accum dust
   real(r8), pointer :: num_fine(:,:)   ! number m.r. of fine dust
   real(r8), pointer :: num_pcarbon(:,:)! number m.r. of primary carbon

   real(r8), pointer :: kvh(:,:)        ! vertical eddy diff coef (m2 s-1)
   real(r8), pointer :: tke(:,:)        ! TKE from the UW PBL scheme (m2 s-2)
   real(r8), pointer :: wp2(:,:)        ! CLUBB vertical velocity variance

   real(r8), pointer :: cldn(:,:)       ! cloud fraction
   real(r8), pointer :: cldo(:,:)       ! old cloud fraction

   real(r8), pointer :: dgnumwet(:,:,:) ! aerosol mode diameter

   real(r8), pointer :: aer_mmr(:,:)    ! aerosol mass mixing ratio

   real(r8)          :: icecldf(pcols,pver)    ! ice cloud fraction   
   real(r8)          :: liqcldf(pcols,pver)    ! liquid cloud fraction

   real(r8) :: rho(pcols,pver)     ! air density (kg m-3)

   real(r8) :: lcldm(pcols,pver)   ! liq cloud fraction

   real(r8) :: lcldn(pcols,pver)   ! fractional coverage of new liquid cloud
   real(r8) :: lcldo(pcols,pver)   ! fractional coverage of old liquid cloud
   real(r8) :: qcld                ! total cloud water
   real(r8) :: nctend_mixnuc(pcols,pver)
   real(r8) :: dum, dum2           ! temporary dummy variable
   real(r8) :: dmc, ssmc           ! variables for modal scheme.

   real(r8) :: so4_num                               ! so4 aerosol number (#/cm^3)
   real(r8) :: soot_num                              ! soot (hydrophilic) aerosol number (#/cm^3)
   real(r8) :: dst1_num,dst2_num,dst3_num,dst4_num   ! dust aerosol number (#/cm^3)
   real(r8) :: organic_num                           
   real(r8) :: dst_num                               ! total dust aerosol number (#/cm^3)

   real(r8) :: qs(pcols)            ! liquid-ice weighted sat mixing rat (kg/kg)
   real(r8) :: es(pcols)            ! liquid-ice weighted sat vapor press (pa)
   real(r8) :: gammas(pcols)        ! parameter for cond/evap of cloud water

   ! bulk aerosol variables
   real(r8), allocatable :: naer2(:,:,:)    ! bulk aerosol number concentration (1/m3)
   real(r8), allocatable :: maerosol(:,:,:) ! bulk aerosol mass conc (kg/m3)

   real(r8) :: wsub(pcols,pver)    ! diagnosed sub-grid vertical velocity st. dev. (m/s)
   real(r8) :: wsubi(pcols,pver)   ! diagnosed sub-grid vertical velocity ice (m/s)
   real(r8) :: wsubice(pcols,pver) ! final updraft velocity for ice nucleation (m/s)
   real(r8) :: wsig(pcols,pver)    ! diagnosed standard deviation of vertical velocity ~ f(TKE)
   real(r8) :: nucboast

   real(r8) :: w0(pcols,pver)      ! large scale velocity (m/s) 
   real(r8) :: w2(pcols,pver)      ! subgrid mean updraft velocity, Gaussian PDF, stddev=f(tke)

   real(r8) :: wght

   real(r8), allocatable :: factnum(:,:,:) ! activation fraction for aerosol number
   !-------------------------------------------------------------------------------

   associate( &
      lchnk => state%lchnk,             &
      ncol  => state%ncol,              &
      t     => state%t,                 &
      qc    => state%q(:pcols,:pver,cldliq_idx), &
      qi    => state%q(:pcols,:pver,cldice_idx), &
      nc    => state%q(:pcols,:pver,numliq_idx), &
      omega => state%omega,             &
      pmid  => state%pmid               )


   call t_startf('microp_aero_run_init')

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx,      ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   if(liqcf_fix) then
      call pbuf_get_field(pbuf, alst_idx,     alst, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, aist_idx,     aist, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   endif

   if(liqcf_fix) then !BSINGH(09/22/2014): It mimics compl treatment
      liqcldf(:ncol,:pver) = alst(:ncol,:pver) 
      icecldf(:ncol,:pver) = aist(:ncol,:pver)
   else
      liqcldf(:ncol,:pver) = ast(:ncol,:pver)
      icecldf(:ncol,:pver) = ast(:ncol,:pver)
   endif

   call pbuf_get_field(pbuf, npccn_idx, npccn)

   call pbuf_get_field(pbuf, nacon_idx, nacon)
   call pbuf_get_field(pbuf, rndst_idx, rndst)

   if (clim_modal_aero) then

      itim_old = pbuf_old_tim_idx()
      
      if (micro_do_icesupersat) then
        call pbuf_get_field(pbuf, cldo_idx, cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))        
      else
        call pbuf_get_field(pbuf, ast_idx,  cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
      endif

      call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )


      call rad_cnst_get_info(0, nmodes=nmodes)
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )

      allocate(factnum(pcols,pver,nmodes))

   end if

   ! initialize output
   npccn(1:ncol,1:pver)    = 0._r8  

   nacon(1:ncol,1:pver,:)  = 0._r8

   ! set default or fixed dust bins for contact freezing
   rndst(1:ncol,1:pver,1) = rn_dst1
   rndst(1:ncol,1:pver,2) = rn_dst2
   rndst(1:ncol,1:pver,3) = rn_dst3
   rndst(1:ncol,1:pver,4) = rn_dst4

   ! save copy of cloud borne aerosols for use in heterogeneous freezing
   if (use_hetfrz_classnuc) then
      call hetfrz_classnuc_cam_save_cbaero(state, pbuf)
   end if

   ! initialize time-varying parameters
   do k = top_lev, pver
      do i = 1, ncol
         rho(i,k) = pmid(i,k)/(rair*t(i,k))
      end do
   end do

   if (clim_modal_aero) then
      ! mode number mixing ratios
      call rad_cnst_get_mode_num(0, mode_coarse_dst_idx, 'a', state, pbuf, num_coarse)
      if(dem_in)then
         if(mode_fine_dst_idx > 0)call rad_cnst_get_mode_num(0, mode_fine_dst_idx, 'a', state, pbuf, num_fine)
         if(mode_pcarbon_idx  > 0)call rad_cnst_get_mode_num(0, mode_pcarbon_idx, 'a', state, pbuf, num_pcarbon)      
      endif
      ! mode specie mass m.r.
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust_idx, 'a', state, pbuf, coarse_dust)
      call rad_cnst_get_aer_mmr(0, mode_coarse_slt_idx, coarse_nacl_idx, 'a', state, pbuf, coarse_nacl)
      if(dem_in)call rad_cnst_get_aer_mmr(0, mode_accum_idx, accum_dust_idx, 'a', state, pbuf, accum_dust) 

   else
      ! init number/mass arrays for bulk aerosols
      allocate( &
         naer2(pcols,pver,naer_all), &
         maerosol(pcols,pver,naer_all))

      do m = 1, naer_all
         call rad_cnst_get_aer_mmr(0, m, state, pbuf, aer_mmr)
         maerosol(:ncol,:,m) = aer_mmr(:ncol,:)*rho(:ncol,:)
         
         if (m .eq. idxsul) then
            naer2(:ncol,:,m) = maerosol(:ncol,:,m)*num_to_mass_aer(m)*bulk_scale
         else
            naer2(:ncol,:,m) = maerosol(:ncol,:,m)*num_to_mass_aer(m)
         end if
      end do
   end if

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! More refined computation of sub-grid vertical velocity 
   ! Set to be zero at the surface by initialization.

   select case (trim(eddy_scheme))
   case ('diag_TKE', 'SHOC_SGS')
      call pbuf_get_field(pbuf, tke_idx, tke)
   case ('CLUBB_SGS')
      itim_old = pbuf_old_tim_idx()
      call pbuf_get_field(pbuf, wp2_idx, wp2, start=(/1,1,itim_old/),kount=(/pcols,pverp,1/))
      allocate(tke(pcols,pverp))
      tke(:ncol,:) = (3._r8/2._r8)*wp2(:ncol,:)

   case default
      call pbuf_get_field(pbuf, kvh_idx, kvh)
   end select

   ! Set minimum values above top_lev.

   !PMA no longer needs the minimum value that is designed for CAM5-UW scheme which 
   !produces very low values

   wsub(:ncol,:top_lev-1)  = wsubmin
   wsubi(:ncol,:top_lev-1) = 0.001_r8
   wsig(:ncol,:top_lev-1)  = 0.001_r8

   do k = top_lev, pver
      do i = 1, ncol

         select case (trim(eddy_scheme))
         case ('diag_TKE', 'CLUBB_SGS', 'SHOC_SGS')
            wsub(i,k) = sqrt(0.5_r8*(tke(i,k) + tke(i,k+1))*(2._r8/3._r8))
!            write(*,*) 'WSUB ', wsub(i,k), tke(i,k)
            wsub(i,k) = min(wsub(i,k),10._r8)
            wsig(i,k) = max(0.001_r8, wsub(i,k))
         case default 
            ! get sub-grid vertical velocity from diff coef.
            ! following morrison et al. 2005, JAS
            ! assume mixing length of 30 m
            dum = (kvh(i,k) + kvh(i,k+1))/2._r8/30._r8
            ! use maximum sub-grid vertical vel of 10 m/s
            dum = min(dum, 10._r8)
            ! set wsub to value at current vertical level
            wsub(i,k)  = dum
         end select

         if (eddy_scheme == 'CLUBB_SGS') then
            wsubi(i,k) = max(0.2_r8, wsub(i,k))
            wsubi(i,k) = min(wsubi(i,k), 10.0_r8)
         else
            wsubi(i,k) = max(0.001_r8, wsub(i,k))
            if (.not. use_preexisting_ice) then
               wsubi(i,k) = min(wsubi(i,k), 0.2_r8)
            endif
         endif

         wsub(i,k)  = max(wsubmin, wsub(i,k))

      end do
   end do

   !!.......................................................... 
   !! Initialization
   !!.......................................................... 

   w0(1:ncol,1:pver) = 0._r8
   w2(1:ncol,1:pver) = 0._r8
   wsubice(1:ncol,1:pver) = 0._r8

   !!.......................................................... 
   !!  Convert from omega to w 
   !!  Negative omega means rising motion
   !!.......................................................... 

   do k = top_lev, pver
      do i = 1, ncol
         w0(i,k) = -1._r8*omega(i,k)/(rho(i,k)*gravit)
      enddo
   enddo

   call t_stopf('microp_aero_run_init')

   !!.......................................................... 
   !! icenul_wsub_scheme = 2 : Mean updraft calculated from Gausssian PDF, with
   !stddev=f(TKE)    
   !!.......................................................... 

   call t_startf('subgrid_mean_updraft')
   call subgrid_mean_updraft(ncol, w0, wsig, w2)
   call t_stopf('subgrid_mean_updraft')


   select case (icenul_wsub_scheme)

   case(1)
         wsubice(1:ncol,1:pver) = wsubi(1:ncol,1:pver)
   case(2)
         wsubice(1:ncol,1:pver) = w2(1:ncol,1:pver)
   case default
         call endrun('nucleate_ice_cam_calc : icenul_wsub_scheme not set')
   end select

   call outfld('WSUB',   wsub, pcols, lchnk)
   call outfld('WSUBI',  wsubice, pcols, lchnk)
   call outfld('WSIG',   wsig, pcols, lchnk)
   call outfld('WLARGE', w0, pcols, lchnk)
   call outfld('WSUBI2', w2, pcols, lchnk)



   if (trim(eddy_scheme) == 'CLUBB_SGS') deallocate(tke)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !ICE Nucleation

   call t_startf('nucleate_ice_cam_calc')
   call nucleate_ice_cam_calc(state, wsubice, pbuf)
   call t_stopf('nucleate_ice_cam_calc')

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! get liquid cloud fraction, check for minimum

   do k = top_lev, pver
      do i = 1, ncol
         lcldm(i,k) = max(ast(i,k), mincld)
      end do
   end do

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Droplet Activation

   if (clim_modal_aero) then

      ! for modal aerosol

      ! partition cloud fraction into liquid water part
      lcldn = 0._r8
      lcldo = 0._r8
      do k = top_lev, pver
         do i = 1, ncol
            qcld = qc(i,k) + qi(i,k)
            if (qcld > qsmall) then
               if(liqcf_fix) then
                  lcldn(i,k)=liqcldf(i,k)
                  lcldo(i,k)=liqcldfo(i,k)
               else
                  lcldn(i,k) = cldn(i,k)*qc(i,k)/qcld
                  lcldo(i,k) = cldo(i,k)*qc(i,k)/qcld
               endif
            end if
         end do
      end do

      call outfld('LCLOUD', lcldn, pcols, lchnk)

      call t_startf('dropmixnuc')
      call dropmixnuc( &
         state, ptend, deltatin, pbuf, wsub, &
         lcldn, lcldo, nctend_mixnuc, factnum)
      call t_stopf('dropmixnuc')

      npccn(:ncol,:) = nctend_mixnuc(:ncol,:)

   else

      ! for bulk aerosol

      ! no tendencies returned from ndrop_bam_run, so just init ptend here
      call physics_ptend_init(ptend, state%psetcols, 'mic_aero_run')

      call t_startf('droplet_act_bulk_aero')
      do k = top_lev, pver
         do i = 1, ncol

            if (qc(i,k) >= qsmall) then

               ! get droplet activation rate

               call ndrop_bam_run( &
                  wsub(i,k), t(i,k), rho(i,k), naer2(i,k,:), naer_all, &
                  naer_all, maerosol(i,k,:),  &
                  dum2)
               dum = dum2
            else
               dum = 0._r8
            end if

            npccn(i,k) = (dum*lcldm(i,k) - nc(i,k))/deltatin
         end do
      end do
      call t_stopf('droplet_act_bulk_aero')

   end if


   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Contact freezing  (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
   ! estimate rndst and nanco for 4 dust bins here to pass to MG microphysics

   do k = top_lev, pver
      do i = 1, ncol

         if (t(i,k) < 269.15_r8) then

            if (clim_modal_aero) then

               ! For modal aerosols:
               !  use size '3' for dust coarse mode...
               !  scale by dust fraction in coarse mode
               
               dmc  = coarse_dust(i,k)
               ssmc = coarse_nacl(i,k)

               if ( separate_dust ) then
                  ! 7-mode -- has separate dust and seasalt mode types and no need for weighting 
                  wght = 1._r8
               else
                  ! 3-mode -- needs weighting for dust since dust and seasalt are combined in the "coarse" mode type
                  wght = dmc/(ssmc + dmc)
               endif

               if (dmc > 0.0_r8) then
                  nacon(i,k,3) = wght*num_coarse(i,k)*rho(i,k)
               else
                  nacon(i,k,3) = 0._r8
               end if

               !also redefine parameters based on size...

               rndst(i,k,3) = 0.5_r8*dgnumwet(i,k,mode_coarse_dst_idx)
               if (rndst(i,k,3) <= 0._r8) then 
                  rndst(i,k,3) = rn_dst3
               end if

            else

               !For Bulk Aerosols: set equal to aerosol number for dust for bins 2-4 (bin 1=0)

               if (idxdst2 > 0) then 
                  nacon(i,k,2) = naer2(i,k,idxdst2)
               end if
               if (idxdst3 > 0) then 
                  nacon(i,k,3) = naer2(i,k,idxdst3)
               end if
               if (idxdst4 > 0) then 
                  nacon(i,k,4) = naer2(i,k,idxdst4)
               end if
            end if

         end if
      end do
   end do

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !bulk aerosol ccn concentration (modal does it in ndrop, from dropmixnuc)

   if (.not. clim_modal_aero) then

      ! ccn concentration as diagnostic
      call t_startf('ndrop_bam_ccn')
      call ndrop_bam_ccn(lchnk, ncol, maerosol, naer2)
      call t_stopf('ndrop_bam_ccn')

      deallocate( &
         naer2,    &
         maerosol)

   end if

   ! heterogeneous freezing
   if (use_hetfrz_classnuc) then

      call t_startf('hetfrz_classnuc_cam_calc')
      call hetfrz_classnuc_cam_calc(state, deltatin, factnum, pbuf)
      call t_stopf('hetfrz_classnuc_cam_calc')

   end if

   if (clim_modal_aero) then
      deallocate(factnum)
   end if

   end associate

end subroutine microp_aero_run

!=========================================================================================

subroutine subgrid_mean_updraft(ncol, w0, wsig, ww)

!---------------------------------------------------------------------------------
! Purpose: Calculate the mean updraft velocity inside a GCM grid assuming the 
!          vertical velocity distribution is Gaussian and peaks at the 
!          GCM resolved large-scale vertical velocity. 
!          When icenul_wsub_scheme = 2, the model uses the mean updraft velocity as the 
!          characteristic updraft velocity to calculate the ice nucleation rate. 
! Author:  Kai Zhang (kai.zhang@pnnl.gov) 
! Last Modified: Oct, 2015 
!---------------------------------------------------------------------------------

   !! interface 

   integer,  intent(in) :: ncol              ! number of cols 
   real(r8), intent(in) :: wsig(pcols,pver ) ! standard deviation (m/s)
   real(r8), intent(in) :: w0(pcols,pver ) ! large scale vertical velocity (m/s) 
   real(r8), intent(out):: ww(pcols,pver) ! mean updraft velocity(m/s) -> characteristic w*

   !! local 
   integer, parameter :: nbin = 50

   real(r8) :: wlarge,sigma
   real(r8) :: xx, yy 
   real(r8) :: zz(nbin) 
   real(r8) :: wa(nbin) 
   integer  :: kp(nbin) 
   integer  :: i, k
   integer  :: ibin

   !! program begins 

   do k = 1, pver
   do i = 1, ncol

      sigma  = max(0.001_r8, wsig(i,k))
      wlarge = w0(i,k)

      xx = 6._r8 * sigma / nbin

      do ibin = 1, nbin
         yy = wlarge - 3._r8*sigma + 0.5*xx
         yy = yy + (ibin-1)*xx
         !! wbar = integrator < w * f(w) * dw > 
         zz(ibin) = yy * exp(-1.*(yy-wlarge)**2/(2*sigma**2))/(sigma*sqrt(2*pi))*xx
      end do 

      kp(:) = 0 
      wa(:) = 0._r8 
 
      where(zz.gt.0._r8) 
         kp = 1 
         wa = zz
      elsewhere 
         kp = 0 
         wa = 0._r8 
      end where 

      if(sum(kp).gt.0) then 
         !! wbar = integrator < w * f(w) * dw > 
         ww(i,k) = sum(wa)
      else 
         ww(i,k) = 0.001_r8
      end if 

      !!write(6,*) 'i, k, w0, wsig, ww : ', i, k, w0(i,k), wsig(i,k), ww(i,k) 

  end do
  end do

end subroutine subgrid_mean_updraft
!================================================================================================

end module microp_aero
