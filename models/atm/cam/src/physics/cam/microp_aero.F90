module microp_aero

!---------------------------------------------------------------------------------
! Purpose:
!   CAM Interface for aerosol activation 
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
!                microphysics and macrophysics and concentrate aerosol information here
!
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8=>shr_kind_r8
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols, pver, pverp
use physconst,        only: rair, tmelt
use constituents,     only: cnst_get_ind, pcnst
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,     only: phys_getopts
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num, rad_cnst_get_mode_props
use shr_spfn_mod,     only: erf => shr_spfn_erf, &
                            erfc => shr_spfn_erfc
use wv_saturation,    only: qsat_water
use nucleate_ice,     only: nucleati
use ndrop,            only: ndrop_init, dropmixnuc
use ndrop_bam,        only: ndrop_bam_init, ndrop_bam_run, ndrop_bam_ccn
use cam_history,      only: addfld, phys_decomp, add_default, outfld
use cam_logfile,      only: iulog
use abortutils,       only: endrun
use ref_pres,         only: top_lev => trop_cloud_top_lev

implicit none
private
save

public :: microp_aero_init, microp_aero_run, microp_aero_readnl, microp_aero_register

! Private module data

character(len=16)   :: eddy_scheme  ! eddy scheme

! contact freezing due to dust
! dust number mean radius (m), Zender et al JGR 2003 assuming number mode radius of 0.6 micron, sigma=2
real(r8), parameter :: rn_dst1 = 0.258e-6_r8
real(r8), parameter :: rn_dst2 = 0.717e-6_r8
real(r8), parameter :: rn_dst3 = 1.576e-6_r8
real(r8), parameter :: rn_dst4 = 3.026e-6_r8

real(r8), public :: bulk_scale    ! prescribed aerosol bulk sulfur scale factor

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
integer :: cldo_idx = -1
integer :: dgnum_idx    = -1
integer :: dgnumwet_idx = -1

! Bulk aerosols
character(len=20), allocatable :: aername(:)
real(r8), allocatable :: num_to_mass_aer(:)

integer :: naer_all      ! number of aerosols affecting climate
integer :: idxsul   = -1 ! index in aerosol list for sulfate
integer :: idxdst1  = -1 ! index in aerosol list for dust1
integer :: idxdst2  = -1 ! index in aerosol list for dust2
integer :: idxdst3  = -1 ! index in aerosol list for dust3
integer :: idxdst4  = -1 ! index in aerosol list for dust4
integer :: idxbcphi = -1 ! index in aerosol list for Soot (BCPHIL)

! modal aerosols
logical :: prog_modal_aero
logical :: clim_modal_aero

integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
integer :: mode_coarse_idx = -1  ! index of coarse mode
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode

integer :: naai_idx, naai_hom_idx, npccn_idx, rndst_idx, nacon_idx

real(r8) :: sigmag_aitken
logical  :: separate_dust = .false.

!===============================================================================
contains
!===============================================================================
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

   call pbuf_add_field('NAAI',       'physpkg',dtype_r8,(/pcols,pver/), naai_idx)
   call pbuf_add_field('NAAI_HOM',   'physpkg',dtype_r8,(/pcols,pver/), naai_hom_idx)
   call pbuf_add_field('NPCCN',      'physpkg',dtype_r8,(/pcols,pver/), npccn_idx)
   call pbuf_add_field('RNDST',      'physpkg',dtype_r8,(/pcols,pver,4/), rndst_idx)
   call pbuf_add_field('NACON',      'physpkg',dtype_r8,(/pcols,pver,4/), nacon_idx)
 

end subroutine microp_aero_register

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
   integer  :: iaer
   integer  :: m, n, nmodes, nspec

   character(len=32) :: str32
   character(len=*), parameter :: routine = 'microp_aero_init'
   logical :: history_amwg
   !-----------------------------------------------------------------------

   ! Query the PBL eddy scheme
   call phys_getopts(eddy_scheme_out          = eddy_scheme, &
                     history_amwg_out = history_amwg)

   ! Access the physical properties of the aerosols that are affecting the climate
   ! by using routines from the rad_constituents module.

   ! get indices into state and pbuf structures
   call cnst_get_ind('CLDLIQ', cldliq_idx)
   call cnst_get_ind('CLDICE', cldice_idx)
   call cnst_get_ind('NUMLIQ', numliq_idx)
   call cnst_get_ind('NUMICE', numice_idx)

   select case(trim(eddy_scheme))
   case ('diag_TKE')
      tke_idx      = pbuf_get_index('tke')
   case ('CLUBB_SGS')
      wp2_idx = pbuf_get_index('WP2')
   case default
      kvh_idx      = pbuf_get_index('kvh')
   end select

   ! prog_modal_aero determines whether prognostic modal aerosols are present in the run.
   call phys_getopts(prog_modal_aero_out=prog_modal_aero)

   ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
   ! The modal aerosols can be either prognostic or prescribed.
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   ast_idx      = pbuf_get_index('AST')

   if (clim_modal_aero) then

      cldo_idx     = pbuf_get_index('CLDO')
      dgnum_idx    = pbuf_get_index('DGNUM' )
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

      ! Check that required mode specie types were found
      if ( coarse_dust_idx == -1 .or. coarse_nacl_idx == -1) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_dust_idx, coarse_nacl_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if

      ! get specific mode properties
      call rad_cnst_get_mode_props(0, mode_aitken_idx, sigmag=sigmag_aitken)

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
         if (trim(aername(iaer)) == 'DUST1') idxdst1 = iaer
         if (trim(aername(iaer)) == 'DUST2') idxdst2 = iaer
         if (trim(aername(iaer)) == 'DUST3') idxdst3 = iaer
         if (trim(aername(iaer)) == 'DUST4') idxdst4 = iaer
         if (trim(aername(iaer)) == 'BCPHIL') idxbcphi = iaer
      end do

      call ndrop_bam_init()

   end if

   call addfld('LCLOUD', ' ', pver, 'A', 'Liquid cloud fraction used in stratus activation', phys_decomp)

   call addfld('WSUB     ', 'm/s     ', pver, 'A', 'Diagnostic sub-grid vertical velocity'                   ,phys_decomp)
   call addfld('WSUBI    ', 'm/s     ', pver, 'A', 'Diagnostic sub-grid vertical velocity for ice'           ,phys_decomp)
   call addfld('NIHF',  '1/m3', pver, 'A', 'Activated Ice Number Concentation due to homogenous freezing',  phys_decomp)
   call addfld('NIDEP', '1/m3', pver, 'A', 'Activated Ice Number Concentation due to deposition nucleation',phys_decomp)
   call addfld('NIIMM', '1/m3', pver, 'A', 'Activated Ice Number Concentation due to immersion freezing',   phys_decomp)
   call addfld('NIMEY', '1/m3', pver, 'A', 'Activated Ice Number Concentation due to meyers deposition',    phys_decomp)

   if (history_amwg) then
      call add_default ('WSUB     ', 1, ' ')
   end if

end subroutine microp_aero_init

!===============================================================================

subroutine microp_aero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Namelist variables
   real(r8) :: microp_aero_bulk_scale = 2._r8  ! prescribed aerosol bulk sulfur scale factor
 
   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'microp_aero_readnl'

   namelist /microp_aero_nl/ microp_aero_bulk_scale
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
#endif

   ! set local variables
   bulk_scale = microp_aero_bulk_scale

end subroutine microp_aero_readnl

!===============================================================================

subroutine microp_aero_run ( &
   state, ptend, deltatin, pbuf)

   ! input arguments
   type(physics_state), target, intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend
   real(r8),                    intent(in)    :: deltatin     ! time step (s)
   type(physics_buffer_desc),   pointer       :: pbuf(:)




   ! local workspace
   ! all units mks unless otherwise stated

   integer :: i, k, m
   integer :: itim_old
   integer :: lchnk
   integer :: ncol
   integer :: nmodes
   integer :: nucboast

   real(r8), pointer :: ast(:,:)        

   real(r8)          :: icecldf(pcols,pver)    ! ice cloud fraction   
   real(r8)          :: liqcldf(pcols,pver)    ! liquid cloud fraction

   real(r8), pointer :: naai(:,:)       ! number of activated aerosol for ice nucleation 
   real(r8), pointer :: naai_hom(:,:)   ! number of activated aerosol for ice nucleation (homogeneous freezing only)
   real(r8), pointer :: npccn(:,:)      ! number of CCN (liquid activated)
   real(r8), pointer :: rndst(:,:,:)    ! radius of 4 dust bins for contact freezing
   real(r8), pointer :: nacon(:,:,:)    ! number in 4 dust bins for contact freezing

   real(r8), pointer :: t(:,:)          ! input temperature (K)
   real(r8), pointer :: qn(:,:)         ! input water vapor mixing ratio (kg/kg)
   ! note: all input cloud variables are grid-averaged
   real(r8), pointer :: qc(:,:)         ! cloud water mixing ratio (kg/kg)
   real(r8), pointer :: qi(:,:)         ! cloud ice mixing ratio (kg/kg)
   real(r8), pointer :: nc(:,:)         ! cloud water number conc (1/kg)
   real(r8), pointer :: ni(:,:)         ! cloud ice number conc (1/kg)
   real(r8), pointer :: pmid(:,:)       ! pressure at layer midpoints (pa)
   real(r8), pointer :: pdel(:,:)       ! pressure difference across level (pa)
   real(r8), pointer :: pint(:,:)       ! air pressure layer interfaces (pa)
   real(r8), pointer :: rpdel(:,:)      ! inverse pressure difference across level (pa)
   real(r8), pointer :: zm(:,:)         ! geopotential height of model levels (m)
   real(r8), pointer :: omega(:,:)      ! vertical velocity (Pa/s)
   real(r8), pointer :: num_accum(:,:)  ! number m.r. of accumulation mode
   real(r8), pointer :: num_aitken(:,:) ! number m.r. of aitken mode
   real(r8), pointer :: num_coarse(:,:) ! number m.r. of coarse mode
   real(r8), pointer :: coarse_dust(:,:) ! mass m.r. of coarse dust
   real(r8), pointer :: coarse_nacl(:,:) ! mass m.r. of coarse nacl

   real(r8), pointer :: kvh(:,:)        ! vertical eddy diff coef (m2 s-1)
   real(r8), pointer :: tke(:,:)        ! TKE from the UW PBL scheme (m2 s-2)
   real(r8), pointer :: wp2(:,:)        ! CLUBB vertical velocity variance

   real(r8), pointer :: cldn(:,:)       ! cloud fraction
   real(r8), pointer :: cldo(:,:)       ! old cloud fraction

   real(r8), pointer :: dgnum(:,:,:)    ! aerosol mode dry diameter
   real(r8), pointer :: dgnumwet(:,:,:) ! aerosol mode diameter

   real(r8), pointer :: aer_mmr(:,:)    ! aerosol mass mixing ratio

   real(r8) :: rho(pcols,pver)     ! air density (kg m-3)
   real(r8) :: relhum(pcols,pver)  ! relative humidity
   real(r8) :: icldm(pcols,pver)   ! ice cloud fraction
   real(r8) :: lcldm(pcols,pver)   ! liq cloud fraction
   real(r8) :: nfice(pcols,pver)   ! fice variable
   real(r8) :: dumfice             ! dummy var in fice calc
   real(r8) :: lcldn(pcols,pver)   ! fractional coverage of new liquid cloud
   real(r8) :: lcldo(pcols,pver)   ! fractional coverage of old liquid cloud
   real(r8) :: qcld                ! total cloud water
   real(r8) :: nctend_mixnuc(pcols,pver)
   real(r8) :: dum, dum2           ! temporary dummy variable
   real(r8) :: dmc, ssmc           ! variables for modal scheme.

   real(r8) :: so4_num                               ! so4 aerosol number (#/cm^3)
   real(r8) :: soot_num                              ! soot (hydrophilic) aerosol number (#/cm^3)
   real(r8) :: dst1_num,dst2_num,dst3_num,dst4_num   ! dust aerosol number (#/cm^3)
   real(r8) :: dst_num                               ! total dust aerosol number (#/cm^3)

   real(r8) :: qs(pcols)            ! liquid-ice weighted sat mixing rat (kg/kg)
   real(r8) :: es(pcols)            ! liquid-ice weighted sat vapor press (pa)
   real(r8) :: gammas(pcols)        ! parameter for cond/evap of cloud water

   ! bulk aerosol variables
   real(r8), allocatable :: naer2(:,:,:)    ! bulk aerosol number concentration (1/m3)
   real(r8), allocatable :: maerosol(:,:,:) ! bulk aerosol mass conc (kg/m3)

   real(r8) :: wsub(pcols,pver)    ! diagnosed sub-grid vertical velocity st. dev. (m/s)
   real(r8) :: wsubi(pcols,pver)   ! diagnosed sub-grid vertical velocity ice (m/s)

   ! history output for ice nucleation
   real(r8) :: nihf(pcols,pver)  !output number conc of ice nuclei due to heterogenous freezing (1/m3)
   real(r8) :: niimm(pcols,pver) !output number conc of ice nuclei due to immersion freezing (hetero nuc) (1/m3)
   real(r8) :: nidep(pcols,pver) !output number conc of ice nuclei due to deoposion nucleation (hetero nuc) (1/m3)
   real(r8) :: nimey(pcols,pver) !output number conc of ice nuclei due to meyers deposition (1/m3)

   real(r8) :: wght

   !-------------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   t     => state%t
   qn    => state%q(:,:,1)
   qc    => state%q(:,:,cldliq_idx)
   qi    => state%q(:,:,cldice_idx)
   nc    => state%q(:,:,numliq_idx)
   ni    => state%q(:,:,numice_idx)
   pmid  => state%pmid
   pdel  => state%pdel
   pint  => state%pint
   rpdel => state%rpdel
   zm    => state%zm
   omega => state%omega

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx,      ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   liqcldf(:ncol,:pver) = ast(:ncol,:pver)
   icecldf(:ncol,:pver) = ast(:ncol,:pver)

   call pbuf_get_field(pbuf, naai_idx, naai)
   call pbuf_get_field(pbuf, naai_hom_idx, naai_hom)
   call pbuf_get_field(pbuf, npccn_idx, npccn)
   call pbuf_get_field(pbuf, nacon_idx, nacon)
   call pbuf_get_field(pbuf, rndst_idx, rndst)

   if (clim_modal_aero) then

      itim_old = pbuf_old_tim_idx()
      call pbuf_get_field(pbuf, ast_idx,  cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
      call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

      call rad_cnst_get_info(0, nmodes=nmodes)
      call pbuf_get_field(pbuf, dgnum_idx,    dgnum,    start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
   end if

   ! initialize output
   naai(1:ncol,1:pver)     = 0._r8  
   naai_hom(1:ncol,1:pver) = 0._r8  
   npccn(1:ncol,1:pver)    = 0._r8  
   nacon(1:ncol,1:pver,:)  = 0._r8

   ! set default or fixed dust bins for contact freezing
   rndst(1:ncol,1:pver,1) = rn_dst1
   rndst(1:ncol,1:pver,2) = rn_dst2
   rndst(1:ncol,1:pver,3) = rn_dst3
   rndst(1:ncol,1:pver,4) = rn_dst4

   ! initialize history output fields for ice nucleation
   nihf(1:ncol,1:pver)  = 0._r8  
   niimm(1:ncol,1:pver) = 0._r8  
   nidep(1:ncol,1:pver) = 0._r8 
   nimey(1:ncol,1:pver) = 0._r8 

   ! initialize time-varying parameters
   do k = top_lev, pver
      do i = 1, ncol
         rho(i,k) = pmid(i,k)/(rair*t(i,k))
      end do
   end do

   if (clim_modal_aero) then
      ! mode number mixing ratios
      call rad_cnst_get_mode_num(0, mode_accum_idx,  'a', state, pbuf, num_accum)
      call rad_cnst_get_mode_num(0, mode_aitken_idx, 'a', state, pbuf, num_aitken)
      call rad_cnst_get_mode_num(0, mode_coarse_dst_idx, 'a', state, pbuf, num_coarse)

      ! mode specie mass m.r.
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust_idx, 'a', state, pbuf, coarse_dust)
      call rad_cnst_get_aer_mmr(0, mode_coarse_slt_idx, coarse_nacl_idx, 'a', state, pbuf, coarse_nacl)

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
   case ('diag_TKE')
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
   wsub(:ncol,:top_lev-1)  = 0.20_r8
   wsubi(:ncol,:top_lev-1) = 0.001_r8

   do k = top_lev, pver
      do i = 1, ncol

         select case (trim(eddy_scheme))
         case ('diag_TKE', 'CLUBB_SGS')
               wsub(i,k) = sqrt(0.5_r8*(tke(i,k) + tke(i,k+1))*(2._r8/3._r8))
               wsub(i,k) = min(wsub(i,k),10._r8)
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

         wsubi(i,k) = max(0.001_r8, wsub(i,k))
         wsubi(i,k) = min(wsubi(i,k), 0.2_r8)
	 
#ifdef CLUBB_SGS
	 if (wsubi(i,k) .le. 0.04_r8) then
           nucboast=100._r8
	   wsubi(i,k)=nucboast*wsubi(i,k)  ! boost ice SGS vertical velocity in CAM-CLUBB
	   				   ! to force nucleation in upper-level stratiform 
					   ! clouds.  Temporary fix until cloud-top radiative
					   ! cooling parameterization is added to CLUBB similar
					   ! to the one of appendix C of Bretherton and Park (2009).  
	 endif
#endif
	 
         wsub(i,k)  = max(0.20_r8, wsub(i,k))
      end do
   end do
   call outfld( 'WSUB'       , wsub,      pcols, lchnk )
   call outfld( 'WSUBI'      , wsubi,     pcols, lchnk )

   if (trim(eddy_scheme) == 'CLUBB_SGS') deallocate(tke)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !Get humidity and saturation vapor pressures

   ! find wet bulk temperature and saturation value for provisional t and q without
   ! condensation

   do k = top_lev, pver

      call qsat_water(t(:ncol,k), pmid(:ncol,k), &
           es(:ncol), qs(:ncol), gam=gammas(:ncol))

      do i = 1, ncol

         relhum(i,k) = qn(i,k)/qs(i)

         ! get cloud fraction, check for minimum
         icldm(i,k) = max(icecldf(i,k), mincld)
         lcldm(i,k) = max(liqcldf(i,k), mincld)

         ! calculate nfice based on liquid and ice mmr (no rain and snow mmr available yet)
         nfice(i,k) = 0._r8
         dumfice    = qc(i,k) + qi(i,k)
         if (dumfice > qsmall .and. qi(i,k) > qsmall) then
            nfice(i,k) = qi(i,k)/dumfice
         end if
      end do
   end do

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !ICE Nucleation

   do k = top_lev, pver
      do i = 1, ncol

         if (t(i,k).lt.tmelt - 5._r8) then


            ! compute aerosol number for so4, soot, and dust with units #/cm^3
            so4_num  = 0._r8
            soot_num = 0._r8
            dst1_num = 0._r8
            dst2_num = 0._r8
            dst3_num = 0._r8
            dst4_num = 0._r8
            dst_num  = 0._r8

            if (clim_modal_aero) then
               !For modal aerosols, assume for the upper troposphere:
               ! soot = accumulation mode
               ! sulfate = aiken mode
               ! dust = coarse mode
               ! since modal has internal mixtures.
               soot_num = num_accum(i,k)*rho(i,k)*1.0e-6_r8
               dmc  = coarse_dust(i,k)*rho(i,k)
               ssmc = coarse_nacl(i,k)*rho(i,k)

               if ( separate_dust ) then
                  ! 7-mode -- has separate dust and seasalt mode types and no need for weighting 
                  wght = 1._r8
               else
                  ! 3-mode -- needs weighting for dust since dust and seasalt are combined in the "coarse" mode type
                  wght = dmc/(ssmc + dmc)
               endif

               if (dmc > 0._r8) then
                  dst_num = wght * num_coarse(i,k)*rho(i,k)*1.0e-6_r8
               else 
                  dst_num = 0.0_r8
               end if

               if (dgnum(i,k,mode_aitken_idx) > 0._r8) then
                  ! only allow so4 with D>0.1 um in ice nucleation
                  so4_num  = num_aitken(i,k)*rho(i,k)*1.0e-6_r8 &
                     * (0.5_r8 - 0.5_r8*erf(log(0.1e-6_r8/dgnum(i,k,mode_aitken_idx))/  &
                     (2._r8**0.5_r8*log(sigmag_aitken))))
               else 
                  so4_num = 0.0_r8 
               end if
               so4_num = max(0.0_r8, so4_num)

            else

               if (idxsul > 0) then 
                  so4_num = naer2(i,k,idxsul)/25._r8 *1.0e-6_r8
               end if
               if (idxbcphi > 0) then 
                  soot_num = naer2(i,k,idxbcphi)/25._r8 *1.0e-6_r8
               end if
               if (idxdst1 > 0) then 
                  dst1_num = naer2(i,k,idxdst1)/25._r8 *1.0e-6_r8
               end if
               if (idxdst2 > 0) then 
                  dst2_num = naer2(i,k,idxdst2)/25._r8 *1.0e-6_r8
               end if
               if (idxdst3 > 0) then 
                  dst3_num = naer2(i,k,idxdst3)/25._r8 *1.0e-6_r8
               end if
               if (idxdst4 > 0) then 
                  dst4_num = naer2(i,k,idxdst4)/25._r8 *1.0e-6_r8
               end if
               dst_num = dst1_num + dst2_num + dst3_num + dst4_num

            end if

            ! *** Turn off soot nucleation ***
            soot_num = 0.0_r8

            call nucleati( &
               wsubi(i,k), t(i,k), relhum(i,k), icldm(i,k), qc(i,k), &
               nfice(i,k), rho(i,k), so4_num, dst_num, soot_num,     &
               naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k))

            naai_hom(i,k) = nihf(i,k)

            ! output activated ice (convert from #/kg -> #/m3)
            nihf(i,k)     = nihf(i,k) *rho(i,k)
            niimm(i,k)    = niimm(i,k)*rho(i,k)
            nidep(i,k)    = nidep(i,k)*rho(i,k)
            nimey(i,k)    = nimey(i,k)*rho(i,k)
         end if
      end do
   end do

   call outfld('NIHF',   nihf, pcols, lchnk)
   call outfld('NIIMM', niimm, pcols, lchnk)
   call outfld('NIDEP', nidep, pcols, lchnk)
   call outfld('NIMEY', nimey, pcols, lchnk)


   if (clim_modal_aero) then

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !droplet activation for modal aerosol

      ! partition cloud fraction into liquid water part
      lcldn = 0._r8
      lcldo = 0._r8
      do k = top_lev, pver
         do i = 1, ncol
            qcld = qc(i,k) + qi(i,k)
            if (qcld > qsmall) then
               lcldn(i,k) = cldn(i,k)*qc(i,k)/qcld
               lcldo(i,k) = cldo(i,k)*qc(i,k)/qcld
            end if
         end do
      end do

      call outfld('LCLOUD', lcldn, pcols, lchnk)

      call dropmixnuc( &
         state, ptend, deltatin, pbuf, wsub, &
         lcldn, lcldo, nctend_mixnuc)

      npccn(:ncol,:) = nctend_mixnuc(:ncol,:)

   else

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !droplet activation for bulk aerosol

      ! no tendencies returned from ndrop_bam_run, so just init ptend here
      call physics_ptend_init(ptend, state%psetcols, 'none')

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

            ! note: deltatin/2.  accounts for sub step in microphysics
            ! ***** This assumes two sub-steps in microphysics.  It's dangerous to 
            ! ***** make that assumption here.  Should move all coding related to 
            ! ***** microphysics substepping into the microphysics.
            npccn(i,k) = (dum - nc(i,k)/lcldm(i,k))/(deltatin/2._r8)*lcldm(i,k)
         end do
      end do

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
      call ndrop_bam_ccn(lchnk, ncol, maerosol, naer2)

      deallocate( &
         naer2,    &
         maerosol)

   end if

end subroutine microp_aero_run

!===============================================================================


!===============================================================================

end module microp_aero

