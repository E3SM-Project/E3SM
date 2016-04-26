module nucleate_ice_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for nucleate_ice module.
!
!  B. Eaton - Sept 2014
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver
use physconst,      only: pi, rair, tmelt
use constituents,   only: cnst_get_ind
use physics_types,  only: physics_state
use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,   only: use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num, rad_cnst_get_mode_props

use physics_buffer, only: pbuf_add_field, dtype_r8, pbuf_old_tim_idx, &
                          pbuf_get_index, pbuf_get_field
use cam_history,    only: addfld, add_default, outfld

use ref_pres,       only: top_lev => trop_cloud_top_lev
use wv_saturation,  only: qsat_water
use shr_spfn_mod,   only: erf => shr_spfn_erf

use cam_logfile,    only: iulog
use cam_abortutils, only: endrun

use nucleate_ice,   only: nucleati_init, nucleati


implicit none
private
save

public :: &
   nucleate_ice_cam_readnl,   &
   nucleate_ice_cam_register, &
   nucleate_ice_cam_init,     &
   nucleate_ice_cam_calc
   

! Namelist variables
logical, public, protected :: use_preexisting_ice = .false.
logical                    :: hist_preexisting_ice = .false.
real(r8)                   :: nucleate_ice_subgrid
real(r8)                   :: so4_sz_thresh_icenuc = huge(1.0_r8) !ice nucleation SO2 size threshold for aitken mode

! Vars set via init method.
real(r8) :: mincld      ! minimum allowed cloud fraction
real(r8) :: bulk_scale  ! prescribed aerosol bulk sulfur scale factor

! constituent indices
integer :: &
   cldliq_idx = -1, &
   cldice_idx = -1, &
   numice_idx = -1

integer :: &
   naai_idx,     &
   naai_hom_idx

integer :: &
   ast_idx   = -1, &
   dgnum_idx = -1

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
logical :: clim_modal_aero

integer :: nmodes = -1
integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
integer :: mode_coarse_idx = -1  ! index of coarse mode
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode

integer :: coarse_so4_idx = -1  ! index of so4 in coarse mode

#if (defined MODAL_AERO_4MODE_MOM)
integer :: coarse_mom_idx = -1  ! index of mom in coarse mode
#endif

#if (defined RAIN_EVAP_TO_COARSE_AERO) 
integer :: coarse_bc_idx = -1  ! index of bc in coarse mode
integer :: coarse_pom_idx = -1  ! index of pom in coarse mode
integer :: coarse_soa_idx = -1  ! index of soa in coarse mode
#endif

integer :: mode_fine_dst_idx = -1   ! index of dust in fine dust mode
integer :: mode_pcarbon_idx  = -1  ! index of dust in accum mode
integer :: accum_dust_idx    = -1  ! index of dust in accum mode
logical :: dem_in            = .false.           ! use DeMott IN

logical  :: separate_dust = .false.
real(r8) :: sigmag_aitken

!===============================================================================
contains
!===============================================================================

subroutine nucleate_ice_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'nucleate_ice_cam_readnl'

  namelist /nucleate_ice_nl/ use_preexisting_ice, hist_preexisting_ice, &
       nucleate_ice_subgrid,so4_sz_thresh_icenuc

  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'nucleate_ice_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, nucleate_ice_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

  end if

#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(use_preexisting_ice,  1, mpilog, 0, mpicom)
  call mpibcast(hist_preexisting_ice, 1, mpilog, 0, mpicom)
  call mpibcast(nucleate_ice_subgrid, 1, mpir8, 0, mpicom)
  call mpibcast(so4_sz_thresh_icenuc, 1, mpir8, 0, mpicom)
#endif

end subroutine nucleate_ice_cam_readnl

!================================================================================================

subroutine nucleate_ice_cam_register()

   call pbuf_add_field('NAAI',     'physpkg', dtype_r8, (/pcols,pver/), naai_idx)
   call pbuf_add_field('NAAI_HOM', 'physpkg', dtype_r8, (/pcols,pver/), naai_hom_idx)

end subroutine nucleate_ice_cam_register

!================================================================================================

subroutine nucleate_ice_cam_init(mincld_in, bulk_scale_in)

   real(r8), intent(in) :: mincld_in
   real(r8), intent(in) :: bulk_scale_in

   ! local variables
   integer  :: iaer
   integer  :: m, n, nspec

   character(len=32) :: str32
   character(len=*), parameter :: routine = 'nucleate_ice_cam_init'
   !--------------------------------------------------------------------------------------------

   mincld     = mincld_in
   bulk_scale = bulk_scale_in

   call cnst_get_ind('CLDLIQ', cldliq_idx)
   call cnst_get_ind('CLDICE', cldice_idx)
   call cnst_get_ind('NUMICE', numice_idx)

   call addfld('NIHF', (/ 'lev' /), 'A',  '1/m3', 'Activated Ice Number Concentation due to homogenous freezing')
   call addfld('NIDEP', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to deposition nucleation')
   call addfld('NIIMM', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to immersion freezing')
   call addfld('NIMEY', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to meyers deposition')

   if (use_preexisting_ice) then
      call addfld('fhom', (/ 'lev' /), 'A', 'fraction', 'Fraction of cirrus where homogeneous freezing occur'   ) 
      call addfld ('WICE', (/ 'lev' /), 'A', 'm/s','Vertical velocity Reduction caused by preexisting ice'  )
      call addfld ('WEFF', (/ 'lev' /), 'A', 'm/s','Effective Vertical velocity for ice nucleation' )
      call addfld ('INnso4', (/ 'lev' /), 'A','1/m3','Number Concentation so4 used for ice_nucleation')
      call addfld ('INnbc', (/ 'lev' /), 'A','1/m3','Number Concentation bc  used for ice_nucleation')
      call addfld ('INndust', (/ 'lev' /), 'A','1/m3','Number Concentation dustused for ice_nucleation')
      call addfld ('INhet', (/ 'lev' /), 'A','1/m3', &
                'contribution for in-cloud ice number density increase by het nucleation in ice cloud')
      call addfld ('INhom', (/ 'lev' /), 'A','1/m3', &
                'contribution for in-cloud ice number density increase by hom nucleation in ice cloud')
      call addfld ('INFrehom',(/ 'lev' /),'A','frequency','hom IN frequency ice cloud')
      call addfld ('INFreIN',(/ 'lev' /),'A','frequency','frequency of ice nucleation occur')

      if (hist_preexisting_ice) then
         call add_default ('WSUBI   ', 1, ' ')  ! addfld/outfld calls are in microp_aero

         call add_default ('fhom    ', 1, ' ') 
         call add_default ('WICE    ', 1, ' ')
         call add_default ('WEFF    ', 1, ' ')
         call add_default ('INnso4  ', 1, ' ')
         call add_default ('INnbc   ', 1, ' ')
         call add_default ('INndust ', 1, ' ')
         call add_default ('INhet   ', 1, ' ')
         call add_default ('INhom   ', 1, ' ')
         call add_default ('INFrehom', 1, ' ')
         call add_default ('INFreIN ', 1, ' ')
      end if
   end if

   ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
   ! The modal aerosols can be either prognostic or prescribed.
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   if (clim_modal_aero) then

      dgnum_idx    = pbuf_get_index('DGNUM' )

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
      if (mode_coarse_dst_idx < 0) mode_coarse_dst_idx = mode_coarse_idx
      if (mode_coarse_slt_idx < 0) mode_coarse_slt_idx = mode_coarse_idx

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


      call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('sulfate')
            coarse_so4_idx = n
         end select
      end do

      ! Check that required mode specie types were found
      if ( coarse_so4_idx == -1) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_so4_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if

#if (defined MODAL_AERO_4MODE_MOM)
      call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('m-organic')
            coarse_mom_idx = n
         end select
      end do

      ! Check that required mode specie types were found
      if ( coarse_mom_idx == -1) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_mom_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if
#endif

#if (defined RAIN_EVAP_TO_COARSE_AERO )
      call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('black-c')
            coarse_bc_idx = n
         end select
      end do

      call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('p-organic')
            coarse_pom_idx = n
         end select
      end do

      call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('s-organic')
            coarse_soa_idx = n
         end select
      end do

      ! Check that required mode specie types were found
      if ( coarse_bc_idx == -1 .or. coarse_pom_idx == -1 .or. coarse_soa_idx == -1 ) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_bc_idx, coarse_pom_idx, coarse_soa_idx 
         call endrun(routine//': ERROR required mode-species type not found')
      end if
#endif

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
   end if


   call nucleati_init(use_preexisting_ice, use_hetfrz_classnuc, iulog, pi, &
        mincld, nucleate_ice_subgrid)

   ! get indices for fields in the physics buffer
   ast_idx      = pbuf_get_index('AST')

end subroutine nucleate_ice_cam_init

!================================================================================================

subroutine nucleate_ice_cam_calc( &
   state, wsubi, pbuf)

   ! arguments
   type(physics_state), target, intent(in)    :: state
   real(r8),                    intent(in)    :: wsubi(:,:)
   type(physics_buffer_desc),   pointer       :: pbuf(:)
 
   ! local workspace

   ! naai and naai_hom are the outputs shared with the microphysics
   real(r8), pointer :: naai(:,:)       ! number of activated aerosol for ice nucleation 
   real(r8), pointer :: naai_hom(:,:)   ! number of activated aerosol for ice nucleation (homogeneous freezing only)

   integer :: lchnk, ncol
   integer :: itim_old
   integer :: i, k, m

   real(r8), pointer :: t(:,:)          ! input temperature (K)
   real(r8), pointer :: qn(:,:)         ! input water vapor mixing ratio (kg/kg)
   real(r8), pointer :: qc(:,:)         ! cloud water mixing ratio (kg/kg)
   real(r8), pointer :: qi(:,:)         ! cloud ice mixing ratio (kg/kg)
   real(r8), pointer :: ni(:,:)         ! cloud ice number conc (1/kg)
   real(r8), pointer :: pmid(:,:)       ! pressure at layer midpoints (pa)

   real(r8), pointer :: num_accum(:,:)  ! number m.r. of accumulation mode
   real(r8), pointer :: num_aitken(:,:) ! number m.r. of aitken mode
   real(r8), pointer :: num_coarse(:,:) ! number m.r. of coarse mode
   real(r8), pointer :: coarse_dust(:,:) ! mass m.r. of coarse dust
   real(r8), pointer :: coarse_nacl(:,:) ! mass m.r. of coarse nacl

   real(r8), pointer :: coarse_so4(:,:) ! mass m.r. of coarse so4

#if (defined MODAL_AERO_4MODE_MOM)
   real(r8), pointer :: coarse_mom(:,:) ! mass m.r. of coarse mom
#endif

#if (defined RAIN_EVAP_TO_COARSE_AERO) 
   real(r8), pointer :: coarse_bc(:,:) ! mass m.r. of coarse bc
   real(r8), pointer :: coarse_pom(:,:) ! mass m.r. of coarse pom
   real(r8), pointer :: coarse_soa(:,:) ! mass m.r. of coarse soa 
#endif

   real(r8), pointer :: aer_mmr(:,:)    ! aerosol mass mixing ratio
   real(r8), pointer :: dgnum(:,:,:)    ! mode dry radius

   real(r8), pointer :: ast(:,:)
   real(r8) :: icecldf(pcols,pver)  ! ice cloud fraction

   real(r8) :: rho(pcols,pver)      ! air density (kg m-3)

   real(r8), allocatable :: naer2(:,:,:)    ! bulk aerosol number concentration (1/m3)
   real(r8), allocatable :: maerosol(:,:,:) ! bulk aerosol mass conc (kg/m3)

   real(r8) :: qs(pcols)            ! liquid-ice weighted sat mixing rat (kg/kg)
   real(r8) :: es(pcols)            ! liquid-ice weighted sat vapor press (pa)
   real(r8) :: gammas(pcols)        ! parameter for cond/evap of cloud water

   real(r8) :: relhum(pcols,pver)  ! relative humidity
   real(r8) :: icldm(pcols,pver)   ! ice cloud fraction

   real(r8) :: so4_num                               ! so4 aerosol number (#/cm^3)
   real(r8) :: soot_num                              ! soot (hydrophilic) aerosol number (#/cm^3)
   real(r8) :: dst1_num,dst2_num,dst3_num,dst4_num   ! dust aerosol number (#/cm^3)
   real(r8) :: organic_num
   real(r8) :: dst_num                               ! total dust aerosol number (#/cm^3)
   real(r8) :: wght
   real(r8) :: dmc
   real(r8) :: ssmc
   real(r8) :: so4mc
   real(r8) :: mommc
   real(r8) :: bcmc
   real(r8) :: pommc
   real(r8) :: soamc

   ! For pre-existing ice
   real(r8) :: fhom(pcols,pver)    ! how much fraction of cloud can reach Shom
   real(r8) :: wice(pcols,pver)    ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at Shom 
   real(r8) :: weff(pcols,pver)    ! effective Vertical velocity for ice nucleation (m/s); weff=wsubi-wice 
   real(r8) :: INnso4(pcols,pver)   ! #/m3, so4 aerosol number used for ice nucleation
   real(r8) :: INnbc(pcols,pver)    ! #/m3, bc aerosol number used for ice nucleation
   real(r8) :: INndust(pcols,pver)  ! #/m3, dust aerosol number used for ice nucleation
   real(r8) :: INhet(pcols,pver)    ! #/m3, ice number from het freezing
   real(r8) :: INhom(pcols,pver)    ! #/m3, ice number from hom freezing
   real(r8) :: INFrehom(pcols,pver) !  hom freezing occurence frequency.  1 occur, 0 not occur.
   real(r8) :: INFreIN(pcols,pver)  !  ice nucleation occerence frequency.   1 occur, 0 not occur.

   ! history output for ice nucleation
   real(r8) :: nihf(pcols,pver)  !output number conc of ice nuclei due to heterogenous freezing (1/m3)
   real(r8) :: niimm(pcols,pver) !output number conc of ice nuclei due to immersion freezing (hetero nuc) (1/m3)
   real(r8) :: nidep(pcols,pver) !output number conc of ice nuclei due to deoposion nucleation (hetero nuc) (1/m3)
   real(r8) :: nimey(pcols,pver) !output number conc of ice nuclei due to meyers deposition (1/m3)


   !-------------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   t     => state%t
   qn    => state%q(:,:,1)
   qc    => state%q(:,:,cldliq_idx)
   qi    => state%q(:,:,cldice_idx)
   ni    => state%q(:,:,numice_idx)
   pmid  => state%pmid

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

      call rad_cnst_get_aer_mmr(0, mode_coarse_idx, coarse_so4_idx, 'a', state, pbuf, coarse_so4)

#if (defined MODAL_AERO_4MODE_MOM)
      call rad_cnst_get_aer_mmr(0, mode_coarse_idx, coarse_mom_idx, 'a', state, pbuf, coarse_mom)
#endif

#if (defined RAIN_EVAP_TO_COARSE_AERO) 
      call rad_cnst_get_aer_mmr(0, mode_coarse_idx, coarse_bc_idx, 'a', state, pbuf, coarse_bc)
      call rad_cnst_get_aer_mmr(0, mode_coarse_idx, coarse_pom_idx, 'a', state, pbuf, coarse_pom)
      call rad_cnst_get_aer_mmr(0, mode_coarse_idx, coarse_soa_idx, 'a', state, pbuf, coarse_soa)
#endif

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

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx, ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   icecldf(:ncol,:pver) = ast(:ncol,:pver)

   if (clim_modal_aero) then
      call pbuf_get_field(pbuf, dgnum_idx, dgnum)
   end if

   ! naai and naai_hom are the outputs from this parameterization
   call pbuf_get_field(pbuf, naai_idx, naai)
   call pbuf_get_field(pbuf, naai_hom_idx, naai_hom)
   naai(1:ncol,1:pver)     = 0._r8  
   naai_hom(1:ncol,1:pver) = 0._r8  

   ! initialize history output fields for ice nucleation
   nihf(1:ncol,1:pver)  = 0._r8  
   niimm(1:ncol,1:pver) = 0._r8  
   nidep(1:ncol,1:pver) = 0._r8 
   nimey(1:ncol,1:pver) = 0._r8 

   if (use_preexisting_ice) then
      fhom(:,:)     = 0.0_r8
      wice(:,:)     = 0.0_r8
      weff(:,:)     = 0.0_r8
      INnso4(:,:)   = 0.0_r8
      INnbc(:,:)    = 0.0_r8
      INndust(:,:)  = 0.0_r8
      INhet(:,:)    = 0.0_r8
      INhom(:,:)    = 0.0_r8
      INFrehom(:,:) = 0.0_r8
      INFreIN(:,:)  = 0.0_r8
   endif

   do k = top_lev, pver

      ! Get humidity and saturation vapor pressures
      call qsat_water(t(:ncol,k), pmid(:ncol,k), &
           es(:ncol), qs(:ncol), gam=gammas(:ncol))

      do i = 1, ncol

         relhum(i,k) = qn(i,k)/qs(i)

         ! get cloud fraction, check for minimum
         icldm(i,k) = max(icecldf(i,k), mincld)

      end do
   end do


   do k = top_lev, pver
      do i = 1, ncol

         if (t(i,k) < tmelt - 5._r8) then

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
               so4mc  = coarse_so4(i,k)*rho(i,k)

#if (defined MODAL_AERO_4MODE_MOM)
               mommc  = coarse_mom(i,k)*rho(i,k)
#endif

#if (defined RAIN_EVAP_TO_COARSE_AERO) 
               bcmc  = coarse_bc(i,k)*rho(i,k)
               pommc  = coarse_pom(i,k)*rho(i,k)
               soamc  = coarse_soa(i,k)*rho(i,k)
#endif

               if (dmc > 0._r8) then
                  if ( separate_dust ) then
                     ! 7-mode -- has separate dust and seasalt mode types and
                     !           no need for weighting 
                     wght = 1._r8
                  else
                     ! 3-mode -- needs weighting for dust since dust and seasalt
                     !           are combined in the "coarse" mode type
#if (defined MODAL_AERO_4MODE_MOM && defined RAIN_EVAP_TO_COARSE_AERO )
                     wght = dmc/(ssmc + dmc + so4mc + bcmc + pommc + soamc + mommc)
#elif (defined MODAL_AERO_4MODE_MOM)
                     wght = dmc/(ssmc + dmc + so4mc + mommc)
#elif (defined RAIN_EVAP_TO_COARSE_AERO) 
                     wght = dmc/(ssmc + dmc + so4mc + bcmc + pommc + soamc)
#else
                     wght = dmc/(ssmc + dmc) ! to keep FC5 binary identical 
#endif

                  endif
                  dst_num = wght * num_coarse(i,k)*rho(i,k)*1.0e-6_r8
               else 
                  dst_num = 0.0_r8
               end if

               if (dgnum(i,k,mode_aitken_idx) > 0._r8) then
                  if (.not. use_preexisting_ice) then
                     ! only allow so4 with D>0.1 um in ice nucleation
                     so4_num  = num_aitken(i,k)*rho(i,k)*1.0e-6_r8 &
                        * (0.5_r8 - 0.5_r8*erf(log(so4_sz_thresh_icenuc/dgnum(i,k,mode_aitken_idx))/  &
                        (2._r8**0.5_r8*log(sigmag_aitken))))
                  else
                     ! all so4 from aitken
                     so4_num  = num_aitken(i,k)*rho(i,k)*1.0e-6_r8
                  end if
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
               wsubi(i,k), t(i,k), pmid(i,k), relhum(i,k), icldm(i,k),   &
               qc(i,k), qi(i,k), ni(i,k), rho(i,k),                      &
               so4_num, dst_num, soot_num,                               &
               naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
               wice(i,k), weff(i,k), fhom(i,k),                          &
               dst1_num,dst2_num,dst3_num,dst4_num,organic_num,          &
               dem_in, clim_modal_aero)

            naai_hom(i,k) = nihf(i,k)

            ! output activated ice (convert from #/kg -> #/m3)
            nihf(i,k)     = nihf(i,k) *rho(i,k)
            niimm(i,k)    = niimm(i,k)*rho(i,k)
            nidep(i,k)    = nidep(i,k)*rho(i,k)
            nimey(i,k)    = nimey(i,k)*rho(i,k)

            if (use_preexisting_ice) then
               INnso4(i,k) =so4_num*1e6_r8  ! (convert from #/cm3 -> #/m3)
               INnbc(i,k)  =soot_num*1e6_r8
               INndust(i,k)=dst_num*1e6_r8
               INFreIN(i,k)=1.0_r8          ! 1,ice nucleation occur
               INhet(i,k) = niimm(i,k) + nidep(i,k)   ! #/m3, nimey not in cirrus
               INhom(i,k) = nihf(i,k)                 ! #/m3
               if (INhom(i,k).gt.1e3_r8)   then ! > 1/L
                  INFrehom(i,k)=1.0_r8       ! 1, hom freezing occur
               endif

               ! exclude  no ice nucleaton 
               if ((INFrehom(i,k) < 0.5_r8) .and. (INhet(i,k) < 1.0_r8))   then   
                  INnso4(i,k) =0.0_r8
                  INnbc(i,k)  =0.0_r8
                  INndust(i,k)=0.0_r8
                  INFreIN(i,k)=0.0_r8
                  INhet(i,k) = 0.0_r8
                  INhom(i,k) = 0.0_r8
                  INFrehom(i,k)=0.0_r8    
                  wice(i,k) = 0.0_r8
                  weff(i,k) = 0.0_r8 
                  fhom(i,k) = 0.0_r8
               endif
            end if

         end if
      end do
   end do

   if (.not. clim_modal_aero) then

      deallocate( &
         naer2,    &
         maerosol)

   end if

   call outfld('NIHF',   nihf, pcols, lchnk)
   call outfld('NIIMM', niimm, pcols, lchnk)
   call outfld('NIDEP', nidep, pcols, lchnk)
   call outfld('NIMEY', nimey, pcols, lchnk)

   if (use_preexisting_ice) then
      call outfld( 'fhom' , fhom, pcols, lchnk)
      call outfld( 'WICE' , wice, pcols, lchnk)
      call outfld( 'WEFF' , weff, pcols, lchnk)
      call outfld('INnso4  ',INnso4 , pcols,lchnk)
      call outfld('INnbc   ',INnbc  , pcols,lchnk)
      call outfld('INndust ',INndust, pcols,lchnk)
      call outfld('INhet   ',INhet  , pcols,lchnk)
      call outfld('INhom   ',INhom  , pcols,lchnk)
      call outfld('INFrehom',INFrehom,pcols,lchnk)
      call outfld('INFreIN ',INFreIN, pcols,lchnk)
   end if

end subroutine nucleate_ice_cam_calc

!================================================================================================

end module nucleate_ice_cam
