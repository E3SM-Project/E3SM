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
!kzm ++
use constituents,   only: pcnst
!use physics_types,  only: physics_state
use rad_constituents, only: rad_cnst_get_mode_num_idx, &
                            rad_cnst_get_mam_mmr_idx
use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
use phys_control,     only: phys_getopts ! namelist variable
!kzm --
use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,   only: use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num, rad_cnst_get_mode_props,rad_cnst_get_mode_num_idx, &
                            rad_cnst_get_mam_mmr_idx
use physics_buffer, only: pbuf_set_field
use wv_saturation,  only: svp_water, svp_ice

use physics_buffer, only: pbuf_add_field, dtype_r8, pbuf_old_tim_idx, &
                          pbuf_get_index, pbuf_get_field
use cam_history,    only: addfld, add_default, outfld

use ref_pres,       only: top_lev => trop_cloud_top_lev
use wv_saturation,  only: qsat_water
#ifndef HAVE_ERF_INTRINSICS
use shr_spfn_mod,   only: erf => shr_spfn_erf
#endif
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
logical, public, protected :: use_nie_nucleate = .false.
logical, public, protected :: use_dem_nucleate = .false.
!real(r8)                   :: nucleate_ice_subgrid
real(r8)                   :: so4_sz_thresh_icenuc = huge(1.0_r8) !ice nucleation SO2 size threshold for aitken mode
logical                    :: nucleate_ice_use_troplev = .false.  !kzm ++
logical                    :: aero_nucleation_removal = .false. !kzm
logical                    :: nucleate_ice_incloud = .false. !kzm ++
real(r8)                   :: nucleate_ice_subgrid = -1._r8 !kzm ++
real(r8)                   :: nucleate_ice_subgrid_strat = -1._r8  !kzm ++
real(r8)                   :: nucleate_ice_strat = 1.0_r8      !kzm ++
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
!kzm ++
logical :: prog_modal_aero
real(r8) :: sigmag_accum
!logical :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
integer :: cnum_coarse_idx, ccoarse_dst_idx, ccoarse_so4_idx
integer :: cnum_strat_coarse_idx, cstrat_coarse_so4_idx!kzm
integer :: accum_so4_idx,accum_pom_idx,accum_soa_idx,accum_bc_idx,accum_dst_idx,accum_ncl_idx,accum_mom_idx
!kzm --
integer :: nmodes = -1
integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
integer :: mode_coarse_idx = -1  ! index of coarse mode
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode

integer :: coarse_so4_idx = -1  ! index of so4 in coarse mode
! kzm ++
! for so4 only aitken mode so4 to get so4_num for nucleation
! if defined MAM7S add stratosphere aitken to troposphere aitken
integer :: mode_strat_sulfate1_idx = -1 
integer :: mode_strat_coarse_idx = -1 
integer :: strat_coarse_so4_idx = -1
real(r8) :: sigmag_coarse, sigmag_strat_coarse
real(r8) :: sigmag_aitken
! kzm --
#if (defined MODAL_AERO_4MODE_MOM  || defined MODAL_AERO_5MODE)
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
integer :: fine_dust_idx    = -1   ! index of dust in fine mode

logical  :: separate_dust = .false.


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
                             use_nie_nucleate, use_dem_nucleate,        &
                             nucleate_ice_subgrid, so4_sz_thresh_icenuc 

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
  call mpibcast(use_nie_nucleate, 1, mpilog, 0, mpicom)
  call mpibcast(use_dem_nucleate, 1, mpilog, 0, mpicom)
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
   logical :: modal_strat_sulfate_ice_nucleation !SICE
   logical :: modal_strat_sulfate_wet_removal 
   !--------------------------------------------------------------------------------------------

   mincld     = mincld_in
   bulk_scale = bulk_scale_in
!kzm ++
      ! Initialize naai.
    ! determine default variables
   call phys_getopts(modal_strat_sulfate_ice_nucleation_out &
                     = modal_strat_sulfate_ice_nucleation,  &
                     modal_strat_sulfate_wet_removal_out    &
                     = modal_strat_sulfate_wet_removal  )             
   nucleate_ice_use_troplev = modal_strat_sulfate_ice_nucleation
   aero_nucleation_removal = modal_strat_sulfate_wet_removal
   prog_modal_aero = aero_nucleation_removal
   if (nucleate_ice_use_troplev) then
        nucleate_ice_strat = 1.0_r8 ! turn on strat_ice_nucleation
   else
        nucleate_ice_strat = -1.0_r8
   end if
!kzm --   
   if( masterproc ) then
      write(iulog,*) 'nucleate_ice parameters:'
      write(iulog,*) '  mincld                     = ', mincld_in !set as 1.0E-04 in WACCM6
      write(iulog,*) '  bulk_scale                 = ', bulk_scale_in !set as 2.0 in WACCM6
      write(iulog,*) '  use_preexisiting_ice       = ', use_preexisting_ice ! T 
      write(iulog,*) '  hist_preexisiting_ice      = ', hist_preexisting_ice ! F
      write(iulog,*) '  nucleate_ice_subgrid       = ', nucleate_ice_subgrid !set as 1.2 in waccm6, modulate the RH staturation value, 
!      write(iulog,*) '  nucleate_ice_subgrid_strat = ', nucleate_ice_subgrid_strat ! this is the same as nucleate_ice_subgrid in
!      WACCM6
      write(iulog,*) '  nucleate_ice_strat         = ', nucleate_ice_strat !set as 1.0 in WACCM6
!      write(iulog,*) '  nucleate_ice_incloud       = ', nucleate_ice_incloud !set as F in WACCM6
      write(iulog,*) '  nucleate_ice_use_troplev   = ', nucleate_ice_use_troplev !set as T in WACCM6
      write(iulog,*) '  modal_strat_sulfate_wet_removal   = ', modal_strat_sulfate_wet_removal 
      write(iulog,*) '  aero_nucleation_removal    = ', aero_nucleation_removal ! remove through nucleation process
   end if

!kzm --


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
         !kzm ++
         case ('strat_sulfate1')
            mode_strat_sulfate1_idx = m
         case ('strat_coarse')	
            mode_strat_coarse_idx = m !kzm MAM5 cse
	     write(iulog,*)'kzm_MAM5_strat_coarse_mode_shown'

         !kzm --   
         end select
      end do

      if (use_nie_nucleate .or. use_dem_nucleate) then
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
            end select
         end do
      end if

      ! check if coarse dust is in separate mode
      separate_dust = mode_coarse_dst_idx > 0

      ! for 3-mode 
      if (mode_coarse_dst_idx < 0) mode_coarse_dst_idx = mode_coarse_idx
      if (mode_coarse_slt_idx < 0) mode_coarse_slt_idx = mode_coarse_idx
      if (mode_fine_dst_idx < 0 .and. (use_nie_nucleate .or. use_dem_nucleate)) mode_fine_dst_idx = mode_accum_idx

      ! Check that required mode types were found
      if (use_nie_nucleate .or. use_dem_nucleate) then
         if (mode_accum_idx == -1 .or. mode_aitken_idx == -1 .or. &
             mode_coarse_dst_idx == -1.or. mode_coarse_slt_idx == -1 .or. &
             mode_fine_dst_idx == -1) then
            write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
               mode_accum_idx, mode_aitken_idx, mode_coarse_dst_idx, mode_coarse_slt_idx, mode_fine_dst_idx
            call endrun(routine//': ERROR required mode type not found')
         end if
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

      if (use_nie_nucleate .or. use_dem_nucleate) then
         call rad_cnst_get_info(0, mode_fine_dst_idx, nspec=nspec)
         do n = 1, nspec
            call rad_cnst_get_info(0, mode_fine_dst_idx, n, spec_type=str32)
            select case (trim(str32))
            case ('dust')
               fine_dust_idx = n
            end select
         end do
      end if

      call rad_cnst_get_info(0, mode_coarse_slt_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_slt_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('seasalt')
            coarse_nacl_idx = n
         end select
      end do

      ! Check that required mode specie types were found
      if (use_nie_nucleate .or. use_dem_nucleate) then
         if ( coarse_dust_idx == -1 .or. coarse_nacl_idx == -1 .or. &
              fine_dust_idx == -1) then
            write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
               coarse_dust_idx, coarse_nacl_idx, fine_dust_idx
            call endrun(routine//': ERROR required mode-species type not found')
         end if
      end if

      if ( coarse_dust_idx == -1 .or. coarse_nacl_idx == -1) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_dust_idx, coarse_nacl_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if
!kzm ++
      if (mode_strat_coarse_idx>0) then
         call rad_cnst_get_info(0, mode_strat_coarse_idx, nspec=nspec)
         do n = 1, nspec
            call rad_cnst_get_info(0, mode_strat_coarse_idx, n, spec_type=str32)
            select case (trim(str32))
            case ('sulfate')
               strat_coarse_so4_idx = n    ! find sulfate
            end select
         end do
      endif
      !for accumulation mode
!      write(iulog,*)'kzm_mode_accum_idx', mode_accum_idx
      if (mode_accum_idx>0) then
         call rad_cnst_get_info(0, mode_accum_idx, nspec=nspec)
         do n = 1, nspec
            call rad_cnst_get_info(0, mode_accum_idx, n, spec_type=str32)
            write(iulog,*) trim(str32)
            select case (trim(str32))
            case ('sulfate')
               accum_so4_idx = n    ! find sulfate
            case ('p-organic')
               accum_pom_idx = n    ! find pom
            case ('s-organic')
               accum_soa_idx = n    ! find soa
            case ('black-c')
               accum_bc_idx = n    ! find bc
            case ('dust')
               accum_dst_idx = n    ! find dust
            case ('seasalt')
               accum_ncl_idx = n    ! find dust 
            case ('m-organic')
               accum_mom_idx = n    ! find dust   
            end select
         end do
      endif
      !write(iulog,*)'kzm_nuc_ini_point1'
!kzm --       

      if (mode_coarse_idx > 0) then
         call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
         do n = 1, nspec
            call rad_cnst_get_info(0, mode_coarse_idx, n, spec_type=str32)
            select case (trim(str32))
            case ('sulfate')
               coarse_so4_idx = n
            end select
         end do
      end if

      ! Check that required mode specie types were found
      if (mode_coarse_idx > 0) then
         if ( coarse_so4_idx == -1) then
            write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
               coarse_so4_idx
            call endrun(routine//': ERROR required mode-species type not found')
         end if
      end if

#if (defined MODAL_AERO_4MODE_MOM  || defined MODAL_AERO_5MODE )
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
      if (use_nie_nucleate .or. use_dem_nucleate) then
         call rad_cnst_get_mode_props(0, mode_coarse_dst_idx, sigmag=sigmag_coarse)
      end if
!kzm ++
      !call rad_cnst_get_mode_props(0, mode_accum_idx, sigmag=sigmag_accum) !kzm ++
!      write(iulog,*)'kzm_nuc_ini_point2' 
      if (prog_modal_aero .and. aero_nucleation_removal) then
         !call rad_cnst_get_mode_num_idx(mode_coarse_dst_idx, cnum_coarse_idx)
         !call rad_cnst_get_mam_mmr_idx(mode_coarse_dst_idx, coarse_dust_idx, ccoarse_dst_idx) !question
         !if (mode_coarse_idx>0) then
         !   call rad_cnst_get_mam_mmr_idx(mode_coarse_idx, coarse_so4_idx,ccoarse_so4_idx) !question
         !end if
         !lq(cnum_coarse_idx) = .true.
         !lq(ccoarse_dst_idx) = .true.
         !lq(ccoarse_so4_idx) = .true.
         !kzm to remove strat coarse sulfate
         !if (mode_strat_coarse_idx>0) then
         !   call rad_cnst_get_mode_num_idx(mode_strat_coarse_idx, cnum_strat_coarse_idx) !strat_coarse mode number
         !   call rad_cnst_get_mam_mmr_idx(mode_strat_coarse_idx, strat_coarse_so4_idx, cstrat_coarse_so4_idx) !strat_coarse mode mass
         !endif
         !lq(cnum_strat_coarse_idx) = .true.
         !lq(cstrat_coarse_so4_idx) = .true.
      endif
!kzm --

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

   call nucleati_init(use_preexisting_ice, use_hetfrz_classnuc,  &
                      use_nie_nucleate, use_dem_nucleate,        &
                      iulog, pi, mincld, nucleate_ice_subgrid)

   ! get indices for fields in the physics buffer
   ast_idx      = pbuf_get_index('AST')

end subroutine nucleate_ice_cam_init

!================================================================================================

subroutine nucleate_ice_cam_calc( &
   state, wsubi, pbuf)
!kzm ++
   use tropopause,          only : tropopause_find, TROP_ALG_HYBSTOB, TROP_ALG_CLIMATE
   use time_manager,   only: get_nstep  
!kzm --   
   ! arguments
   type(physics_state), target, intent(in)    :: state
   real(r8),                    intent(in)    :: wsubi(:,:)
   type(physics_buffer_desc),   pointer       :: pbuf(:)
   integer  :: troplev(pcols)       ! tropopause level 
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
   real(r8), pointer :: fine_dust(:,:)   ! mass m.r. of fine dust
   real(r8), pointer :: coarse_nacl(:,:) ! mass m.r. of coarse nacl

   real(r8), pointer :: coarse_so4(:,:) ! mass m.r. of coarse so4
!kzm ++
   real(r8), pointer :: num_strcrs(:,:)  ! number m.r. of strat. coarse mode !kzm
   real(r8), pointer :: cld_num_coarse(:,:) ! number m.r. of coarse mode
   real(r8), pointer :: cld_num_strat_coarse(:,:) ! number m.r. of strat coarse mode
   real(r8), pointer :: cld_coarse_dust(:,:) ! mass m.r. of coarse dust
   real(r8), pointer :: cld_coarse_so4(:,:) ! mass m.r. of coarse so4
   real(r8), pointer :: cld_strat_coarse_so4(:,:) ! mass m.r. of strat coarse dust
   real(r8), pointer :: strat_coarse_so4(:,:) ! mass m.r. of strat coarse so4
   !accum_so4,accum_pom,accum_soa,accum_bc,accum_dst,accum_ncl,accum_mom
   real(r8), pointer :: accum_so4(:,:)
   real(r8), pointer :: accum_pom(:,:)
   real(r8), pointer :: accum_soa(:,:)
   real(r8), pointer :: accum_bc(:,:)
   real(r8), pointer :: accum_dst(:,:)
   real(r8), pointer :: accum_ncl(:,:)
   real(r8), pointer :: accum_mom(:,:)
   !cld_num_accum, cld_accum_so4
   real(r8), pointer :: cld_num_accum(:,:)
   real(r8), pointer :: cld_accum_so4(:,:)
   real(r8), pointer :: cld_accum_bc(:,:)
!kzm --   

!   real(r8), pointer :: num_strcrs(:,:)  ! number m.r. of strat. coarse mode !kzm 
#if (defined MODAL_AERO_4MODE_MOM  || defined MODAL_AERO_5MODE)
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

!kzm ++
   real(r8) :: oso4_num
   real(r8) :: odst_num
   real(r8) :: osoot_num
   real(r8) :: dso4_num
   real(r8) :: so4_num_ac ! adjusted accum nucleated number
   real(r8) :: so4_num_at!kzm so4 num in aitken																		
   real(r8) :: so4_num_cr,bc_num_accum
   real(r8) :: so4_num_st_cr,fso4_m3,fso4_m5
   real(r8) :: so4mc_accum,pommc_accum,soamc_accum,bcmc_accum,dstmc_accum,nclmc_accum,mommc_accum,so4_num_accum             
   real(r8) :: ramp

  !real(r8) :: subgrid(pcols,pver) !kzm note: not add in for now
   real(r8) :: trop_pd(pcols,pver)
   integer :: nstep                             ! current timestep number
   real(r8) :: regm(pcols,pver)  !output temperature thershold for nucleation regime
   real(r8) :: cld_num_strat_coarse_tend2, num_strcrs_tend2,cld_mass_strat_coarse_so4_tend2,mass_strat_coarse_so4_tend2
   real(r8) :: num_accum_tend2,cld_num_accum_tend2,num_coarse_tend2,cld_num_coarse_tend2
   real(r8) :: cld_mass_accum_so4_tend2,mass_accum_so4_tend2,cld_mass_coarse_so4_tend2, mass_coarse_so4_tend2   
   real(r8) :: cld_num_coarse_tend1, num_coarse_tend1,cld_mass_coarse_dust_tend1, mass_coarse_dust_tend1
   real(r8) :: cld_num_accum_tend1, num_accum_tend1, cld_mass_accum_bc_tend1, mass_accum_bc_tend1
!kzm --


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

   real(r8) :: soot_num_to_mass, dst1_num_to_mass
   real(r8) :: soot_sfc_to_mass, dst1_sfc_to_mass, alnsg
   real(r8) :: dst1_sfc_to_num,  dst3_sfc_to_num
   real(r8) :: soot_sfc, organic_sfc, dst_sfc, &
               dst1_sfc, dst2_sfc, dst3_sfc, dst4_sfc     ! aerosol surface area (m2/cm^3)

   !-------------------------------------------------------------------------------
   nstep = get_nstep() !kzm ++
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
      !kzm ++
      if (mode_strat_coarse_idx > 0) then
          call rad_cnst_get_mode_num(0, mode_strat_coarse_idx,  'a', state, pbuf, num_strcrs) !kzm
      endif
      ! mode specie mass m.r.
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust_idx, 'a', state, pbuf, coarse_dust)
      call rad_cnst_get_aer_mmr(0, mode_coarse_slt_idx, coarse_nacl_idx, 'a', state, pbuf, coarse_nacl)
      !kzm --    
      if (mode_coarse_idx > 0) then
         call rad_cnst_get_aer_mmr(0, mode_coarse_idx, coarse_so4_idx, 'a', state, pbuf, coarse_so4)
      end if

      if (use_nie_nucleate .or. use_dem_nucleate) then
         call rad_cnst_get_aer_mmr(0, mode_fine_dst_idx, fine_dust_idx, 'a', state, pbuf, fine_dust)
      end if

#if (defined MODAL_AERO_4MODE_MOM || defined MODAL_AERO_5MODE)
      call rad_cnst_get_aer_mmr(0, mode_coarse_idx, coarse_mom_idx, 'a', state, pbuf, coarse_mom)
#endif

#if (defined RAIN_EVAP_TO_COARSE_AERO) 
      call rad_cnst_get_aer_mmr(0, mode_coarse_idx, coarse_bc_idx, 'a', state, pbuf, coarse_bc)
      call rad_cnst_get_aer_mmr(0, mode_coarse_idx, coarse_pom_idx, 'a', state, pbuf, coarse_pom)
      call rad_cnst_get_aer_mmr(0, mode_coarse_idx, coarse_soa_idx, 'a', state, pbuf, coarse_soa)
#endif
!kzm ++
     ! call rad_cnst_get_mode_num(0, mode_coarse_dst_idx, 'c', state, pbuf, cld_num_coarse)
     ! call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust_idx, 'c', state, pbuf, cld_coarse_dust)
     ! call physics_ptend_init(ptend, state%psetcols, 'nucleatei', lq=lq) !kzm not now
     !get accumulation mode aerosols
     !accum_so4,accum_pom,accum_soa,accum_bc,accum_dst,accum_ncl,accum_mom
      call rad_cnst_get_aer_mmr(0, mode_accum_idx, accum_so4_idx, 'a', state, pbuf, accum_so4)
      call rad_cnst_get_aer_mmr(0, mode_accum_idx, accum_pom_idx, 'a', state, pbuf, accum_pom)
      call rad_cnst_get_aer_mmr(0, mode_accum_idx, accum_soa_idx, 'a', state, pbuf, accum_soa)
      call rad_cnst_get_aer_mmr(0, mode_accum_idx, accum_bc_idx, 'a', state, pbuf, accum_bc)
      call rad_cnst_get_aer_mmr(0, mode_accum_idx, accum_dst_idx, 'a', state, pbuf, accum_dst)
      call rad_cnst_get_aer_mmr(0, mode_accum_idx, accum_ncl_idx, 'a', state, pbuf, accum_ncl)
      call rad_cnst_get_aer_mmr(0, mode_accum_idx, accum_mom_idx, 'a', state, pbuf, accum_mom)
      !write(iulog,*)'kzm_nuc_cal_point3'
      !write(iulog,*)'kzm_pass_point1_nuc_calc'
      if (aero_nucleation_removal) then
      if (mode_strat_coarse_idx > 0) then
      ! strat coarse mode so4 to removal        
      call rad_cnst_get_aer_mmr(0, mode_strat_coarse_idx, strat_coarse_so4_idx, 'a', state, pbuf, strat_coarse_so4)
      call rad_cnst_get_aer_mmr(0, mode_strat_coarse_idx, strat_coarse_so4_idx, 'c', state, pbuf, cld_strat_coarse_so4)
      call rad_cnst_get_mode_num(0, mode_strat_coarse_idx, 'c', state, pbuf, cld_num_strat_coarse)
      !/ CLD_NUM_STRAT_COARSE
      endif
      ! coarse mode dst and so4 to removal
      call rad_cnst_get_mode_num(0, mode_coarse_dst_idx, 'c', state, pbuf, cld_num_coarse)
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust_idx, 'c', state, pbuf, cld_coarse_dust)
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_so4_idx, 'c', state, pbuf, cld_coarse_so4)
      ! accum mode so4 to removal cld_num_accum, cld_accum_so4
      call rad_cnst_get_mode_num(0, mode_accum_idx, 'c', state, pbuf, cld_num_accum)
      call rad_cnst_get_aer_mmr(0, mode_accum_idx, accum_so4_idx, 'c', state, pbuf, cld_accum_so4)
      call rad_cnst_get_aer_mmr(0, mode_accum_idx, accum_bc_idx, 'c', state, pbuf, cld_accum_bc)
      !call physics_ptend_init(ptend, state%psetcols, 'nucleatei', lq=lq) !kzm
      endif
!kzm --
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
!kzm ++
   ! Use the same criteria that is used in chemistry and in CLUBB (for cloud fraction)
   ! to determine whether to use tropospheric or stratospheric settings. Include the
   ! tropopause level so that the cold point tropopause will use the stratospheric values.
   !call tropopause_findChemTrop(state, troplev)
   call tropopause_find(state, tropLev, primary=TROP_ALG_HYBSTOB, backup=TROP_ALG_CLIMATE)

!   if ((nucleate_ice_subgrid .eq. -1._r8) .or. (nucleate_ice_subgrid_strat .eq. -1._r8)) then
!      call pbuf_get_field(pbuf, qsatfac_idx, qsatfac)
!   end if

 !  trop_pd(:,:) = 0._r8

!   do k = top_lev, pver
!      do i = 1, ncol
!         trop_pd(i, troplev(i)) = 1._r8
!
!         if (k <= troplev(i)) then
!            if (nucleate_ice_subgrid_strat .eq. -1._r8) then
!               subgrid(i, k) = 1._r8 / qsatfac(i, k)
!            else
!               subgrid(i, k) = nucleate_ice_subgrid_strat
!            end if
!         else
!            if (nucleate_ice_subgrid .eq. -1._r8) then
!               subgrid(i, k) = 1._r8 / qsatfac(i, k)
!            else
!               subgrid(i, k) = nucleate_ice_subgrid
!            end if
!         end if
!      end do
!   end do
!kzm --


   ! initialize history output fields for ice nucleation
   nihf(1:ncol,1:pver)  = 0._r8  
   niimm(1:ncol,1:pver) = 0._r8  
   nidep(1:ncol,1:pver) = 0._r8 
   nimey(1:ncol,1:pver) = 0._r8 

   soot_num_to_mass = 4.751e+16_r8           ! #/kg
   dst1_num_to_mass = 3.484e+15_r8           ! #/kg, for dust in accumulation mode
   soot_sfc_to_mass = 9213.4_r8              ! m2/kg, Clarke et al., 1997; 2004; 2007: rg=0.1um, sig=1.6
   dst1_sfc_to_mass = 3464.0_r8              ! m2/kg
   dst1_sfc_to_num  = 9.943e-13_r8           ! m2/#, individual particle sfc
   dst3_sfc_to_num  = 0.0_r8

   soot_sfc    = 0.0_r8
   dst_sfc     = 0.0_r8
   organic_sfc = 0.0_r8
   dst1_sfc    = 0.0_r8
   dst2_sfc    = 0.0_r8
   dst3_sfc    = 0.0_r8
   dst4_sfc    = 0.0_r8

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

               if (mode_coarse_idx > 0) then
                  so4mc  = coarse_so4(i,k)*rho(i,k)
               endif

#if (defined MODAL_AERO_4MODE_MOM || defined MODAL_AERO_5MODE)
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
!kzm ++
#elif (defined MODAL_AERO_5MODE && defined RAIN_EVAP_TO_COARSE_AERO )
                     wght = dmc/(ssmc + dmc + so4mc + bcmc + pommc + soamc + mommc)
#elif (defined MODAL_AERO_4MODE_MOM)
                     wght = dmc/(ssmc + dmc + so4mc + mommc)
#elif (defined RAIN_EVAP_TO_COARSE_AERO) 
                     wght = dmc/(ssmc + dmc + so4mc + bcmc + pommc + soamc)
#else
                     wght = dmc/(ssmc + dmc) ! to keep FC5 binary identical 
#endif

                  endif
                  dst3_num = wght * num_coarse(i,k)*rho(i,k)*1.0e-6_r8
               else 
                  dst3_num = 0.0_r8
               end if

               if (use_nie_nucleate .or. use_dem_nucleate) then
                  dst1_num = fine_dust(i,k) * rho(i,k) * dst1_num_to_mass * 1.0e-6_r8

                  alnsg = log(sigmag_coarse)
                  dst3_sfc_to_num = pi*dgnum(i,k,mode_coarse_idx)**2.0_r8*exp(2.0_r8*alnsg**2.0_r8) ! m2/#, individual particle sfc 

                  dst_num = dst1_num + dst3_num
               else
                  dst_num = dst3_num
               end if
!kzm ++
! this part is to calculate the so4 num in the coarse mode (mode 3)
               if ( separate_dust ) then
                  ! 7-mode -- the 7 mode scheme does not support
                  ! stratospheric sulfates, and the sulfates are mixed in
                  ! with the separate soot and dust modes, so just ignore
                  ! for now.
                  so4_num_cr = 0.0_r8
               else
                  ! 3-mode -- needs weighting for dust since dust, seasalt,
                  !           and sulfate are combined in the "coarse" mode
                  !           type
                  so4mc    = coarse_so4(i,k)*rho(i,k)

                  if (so4mc > 0._r8) then
                    wght = so4mc/(ssmc + dmc + so4mc)
                    so4_num_cr = wght * num_coarse(i,k)*rho(i,k)*1.0e-6_r8
                  else
                    so4_num_cr = 0.0_r8
                  end if
               endif

               !accum_so4,accum_pom,accum_soa,accum_bc,accum_dst,accum_ncl,accum_mom
               !so4mc_accum,pommc_accum,soamc_accum,bcmc_accum,dstmc_accum,nclmc_accum,mommc_accum,so4_num_accum
               so4mc_accum = accum_so4(i,k)
               pommc_accum = accum_pom(i,k)
               soamc_accum = accum_soa(i,k)
               bcmc_accum = accum_bc(i,k)
               dstmc_accum = accum_dst(i,k)
               nclmc_accum = accum_ncl(i,k)
               mommc_accum = accum_mom(i,k)
               wght = so4mc_accum/(so4mc_accum+pommc_accum+soamc_accum+bcmc_accum+dstmc_accum+nclmc_accum+mommc_accum)
               so4_num_accum = wght*num_accum(i,k)
               so4_num_accum = so4_num_accum*rho(i,k)*1.0e-6_r8 
               wght = bcmc_accum/(so4mc_accum+pommc_accum+soamc_accum+bcmc_accum+dstmc_accum+nclmc_accum+mommc_accum)
               bc_num_accum = wght*num_accum(i,k)
               bc_num_accum = bc_num_accum*rho(i,k)*1.0e-6_r8
!kzm --


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
            organic_num = 0.0_r8

            call nucleati( &
               wsubi(i,k), t(i,k), pmid(i,k), relhum(i,k), icldm(i,k),   &
               qc(i,k), qi(i,k), ni(i,k), rho(i,k),                      &
               so4_num, dst_num, soot_num,                               &
               dst1_sfc_to_num, dst3_sfc_to_num,                         &
               naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
               wice(i,k), weff(i,k), fhom(i,k),                          &
               dst1_num,dst2_num,dst3_num,dst4_num,organic_num,          &
               clim_modal_aero,                                          &
               oso4_num, odst_num, osoot_num)   !kzm ++
!kzm ++
            ! Move aerosol used for nucleation from interstial to cloudborne,
            ! otherwise the same coarse mode aerosols will be available again
            ! in the next timestep and will supress homogeneous freezing.
            !if (odst_num > 0.0_r8) then
               ! write(iulog,*)'kzm_odst_num', odst_num
               ! write(iulog,*)'kzm_osoot_num', osoot_num
            !endif
            if (prog_modal_aero .and. aero_nucleation_removal) then
            if ( aero_nucleation_removal) then
               if (separate_dust) then
                  call endrun('nucleate_ice_cam: use_preexisting_ice is not supported in separate_dust mode (MAM7)')
               endif
               !write(iulog,*)'kzm_dust_removal_in_cirrus_cloud'
               ! removal coarse dst num
               !if (aero_nucleation_removal)then
                  !cld_num_coarse_tend1, num_coarse_tend1,cld_mass_coarse_dust_tend1,mass_coarse_dust_tend1
                  !cld_num_accum_tend1,num_accum_tend1,cld_mass_accum_bc_tend1,mass_accum_bc_tend1     
                  cld_num_coarse_tend1 = 0.0_r8  
                  num_coarse_tend1 = 0.0_r8

                  cld_mass_coarse_dust_tend1 = 0.0_r8
                  mass_coarse_dust_tend1 = 0.0_r8

                  cld_num_accum_tend1 = 0.0_r8
                  num_accum_tend1 = 0.0_r8

                  cld_mass_accum_bc_tend1 = 0.0_r8
                  mass_accum_bc_tend1 = 0.0_r8



               if (odst_num > 0.0_r8) then
                  ! transfer number from interstial to cloudborne
                  !cld_num_coarse(i,k) = max((cld_num_coarse(i,k) + (odst_num * icldm(i,k))/rho(i,k)/1e-6_r8), 0.0_r8)
                  !num_coarse(i,k) = max((num_coarse(i,k) - (odst_num * icldm(i,k))/rho(i,k)/1e-6_r8), 0.0_r8)
                  cld_num_coarse_tend1 = max( (odst_num * icldm(i,k))/rho(i,k)/1e-6_r8, 0.0_r8)
                  num_coarse_tend1 = -max( (odst_num * icldm(i,k))/rho(i,k)/1e-6_r8, 0.0_r8)

                  ! transfer mass from interstial to cloudborne
                  !cld_coarse_dust(i,k) = max((cld_coarse_dust(i,k) + odst_num / dst_num *icldm(i,k) * coarse_dust(i,k)), 0.0_r8)
                  !coarse_dust(i,k) = max((coarse_dust(i,k) - odst_num / dst_num *icldm(i,k) * coarse_dust(i,k)), 0.0_r8)
                  cld_mass_coarse_dust_tend1 = max((odst_num / dst_num *icldm(i,k) * coarse_dust(i,k)), 0.0_r8)
                  mass_coarse_dust_tend1 = -max((odst_num / dst_num *icldm(i,k) * coarse_dust(i,k)), 0.0_r8)


                  !cnum_coarse_idx is the inter coarse number
                  !this ptend is to remove the interstatial number      
                  !ptend%q(i,k,cnum_coarse_idx) = -(odst_num * icldm(i,k))/rho(i,k)/1e-6_r8/dtime
                  !cld_num_coarse(i,k)   = cld_num_coarse(i,k) + (odst_num * icldm(i,k))/rho(i,k)/1e-6_r8
               ! removal coarse dst mass
                  !ptend%q(i,k,ccoarse_dst_idx) = - odst_num / dst_num * icldm(i,k) * coarse_dust(i,k) / dtime!
                  !cld_coarse_dust(i,k) = cld_coarse_dust(i,k) + odst_num / dst_num *icldm(i,k) * coarse_dust(i,k)
               end if
               if (1 > 2) then
               if (osoot_num > 0.0_r8 ) then
                  ! transfer number from interstial to cloudborne
                  !cld_num_accum(i,k) = max((cld_num_accum(i,k) + (osoot_num * icldm(i,k))/rho(i,k)/1e-6_r8), 0.0_r8)
                  !num_accum(i,k) = max((num_accum(i,k) - (osoot_num * icldm(i,k))/rho(i,k)/1e-6_r8), 0.0_r8)
                  cld_num_accum_tend1 = max(osoot_num * icldm(i,k)/rho(i,k)/1e-6_r8, 0.0_r8)
                  num_accum_tend1 = -max(osoot_num * icldm(i,k)/rho(i,k)/1e-6_r8, 0.0_r8)
                  ! transfer mass from interstial to cloudborne
                  !cld_accum_bc(i,k) = max((cld_accum_bc(i,k) + osoot_num/soot_num *icldm(i,k) *accum_bc(i,k)), 0.0_r8)
                  !accum_bc(i,k) = max((accum_bc(i,k) - osoot_num/soot_num *icldm(i,k) *accum_bc(i,k)), 0.0_r8)
                  cld_mass_accum_bc_tend1 = max((osoot_num/soot_num *icldm(i,k) *accum_bc(i,k)), 0.0_r8)
                  mass_accum_bc_tend1 = -max((osoot_num/soot_num *icldm(i,k) *accum_bc(i,k)), 0.0_r8)
                  !write(iulog,*)'kzm_osoot_num', osoot_num     
               endif
               endif
            endif
            endif


            ! Liu&Penner does not generate enough nucleation in the polar winter
            ! stratosphere, which affects surface area density, dehydration and
            ! ozone chemistry. Part of this is that there are a larger number of
            ! particles in the accumulation mode than in the Aitken mode. In volcanic
            ! periods, the coarse mode may also be important. As a short
            ! term work around, include the accumulation and coarse mode particles
            ! and assume a larger fraction of the sulfates nucleate in the polar
            ! stratosphere.
            !
            ! Do not include the tropopause level, as stratospheric aerosols
            ! only exist above the tropopause level.
            !
            ! NOTE: This may still not represent the proper particles that
            ! participate in nucleation, because it doesn't include STS and NAT
            ! particles. It may not represent the proper saturation threshold for
            ! nucleation, and wsubi from CLUBB is probably not representative of
            ! wave driven varaibility in the polar stratosphere.
            if ( clim_modal_aero) then
               !cld_num_strat_coarse_tend2, num_strcrs_tend2,cld_mass_strat_coarse_so4_tend2,mass_strat_coarse_so4_tend2 
               !num_accum_tend2,cld_num_accum_tend2,num_coarse_tend2,cld_num_coarse_tend2
                     !cld_mass_accum_so4_tend2,mass_accum_so4_tend2,cld_mass_coarse_so4_tend2,mass_coarse_so4_tend2
               cld_num_strat_coarse_tend2 = 0._r8
               num_strcrs_tend2 = 0._r8
               cld_mass_strat_coarse_so4_tend2 = 0._r8
               mass_strat_coarse_so4_tend2 = 0._r8
               !num_accum_tend2,cld_num_accum_tend2,num_coarse_tend2,cld_num_coarse_tend2
                     !cld_mass_accum_so4_tend2,mass_accum_so4_tend2,cld_mass_coarse_so4_tend2,mass_coarse_so4_tend2
               num_accum_tend2 = 0._r8
               cld_num_accum_tend2 = 0._r8
               num_coarse_tend2 = 0._r8
               cld_num_coarse_tend2 = 0._r8
              cld_mass_accum_so4_tend2 = 0._r8
              mass_accum_so4_tend2 = 0._r8
              cld_mass_coarse_so4_tend2 = 0._r8
              mass_coarse_so4_tend2 = 0._r8      
              !num_coarse_tend2,num_accum_tend2

              !if ((k < troplev(i)) .and. (nucleate_ice_strat > 0._r8)) then
              if ( (nucleate_ice_strat > 0._r8) ) then
                 if (oso4_num > 0._r8) then
                    !oso4_num #/cm3     
                    !so4_num_ac = num_accum(i,k)*rho(i,k)*1.0e-6_r8
                    if (mode_strat_coarse_idx > 0._r8) then
                        so4_num_st_cr = num_strcrs(i,k)*rho(i,k)*1.0e-6_r8 !kzm add in stratosphere coarse
                        !so4_num_ac = num_accum(i,k)*rho(i,k)*1.0e-6_r8 ! over write weighted so4_num_ac
                        ! use 10 times aitken nucleated ice to adjust (Barahona and Nenes 2008)
                        !so4_num_ac = min(oso4_num*rho(i,k)*1.0e-6_r8*10.0_r8, so4_num_accum*0.1_r8) 
                        if (k < troplev(i)) then 
                           so4_num_ac = max(0.0_r8, so4_num_accum*0.25_r8) 
                        else
                           so4_num_ac = max(0.0_r8, so4_num_accum*0.05_r8)
                        endif     
                        !write(iulog,*)'kzm_oso4_num', oso4_num, so4_num_accum, so4_num_ac, so4_num_st_cr
                        dso4_num = max(0._r8, (nucleate_ice_strat * (so4_num_cr + so4_num_st_cr + so4_num_ac )) &  !kzm change only include coarse
                                   ) * 1e6_r8 / rho(i,k) !kzm
                    else
                        dso4_num = max(0._r8, (nucleate_ice_strat * (so4_num_cr )) ) * 1e6_r8 / rho(i,k)
                    endif
                    if (2<1) then
                       if (mode_coarse_idx > 0._r8  .and. mode_strat_coarse_idx > 0._r8 ) then
                          write(iulog,*)'kzm_oso4_num', oso4_num, so4_num_accum, so4_num_cr, so4_num_st_cr
                          write(iulog,*)'kzm_nuc_mass_removal', accum_so4(i,k)*icldm(i,k), coarse_so4(i,k)*icldm(i,k), strat_coarse_so4(i,k)*icldm(i,k)
                       else if  (mode_coarse_idx > 0._r8) then
                          write(iulog,*)'kzm_oso4_num', oso4_num, so4_num_accum, so4_num_cr
                          write(iulog,*)'kzm_nuc_mass_removal', accum_so4(i,k)*icldm(i,k), coarse_so4(i,k)*icldm(i,k)
                       endif      
                    endif

                    if (aero_nucleation_removal .and. dso4_num > 0.0_r8) then !if we need to remove it
                    if ( mode_coarse_idx > 0._r8 ) then !MAM4 case
                          !remove coarse mode so4
                      
                           ! transfer number from interstial to cloudborne
                           !accumulation mode
                           !cld_num_accum(i,k) = max((cld_num_accum(i,k) + (so4_num_accum* 1e6_r8 / rho(i,k) * icldm(i,k))), 0.0_r8)
                           !num_accum(i,k) = max((num_coarse(i,k) - ( so4_num_accum* 1e6_r8 / rho(i,k)* icldm(i,k))), 0.0_r8) 
                           cld_num_accum_tend2 = max(((so4_num_accum* 1e6_r8 / rho(i,k) * icldm(i,k))), 0.0_r8)!cld accumulation number
                           num_accum_tend2 = -max(((so4_num_accum* 1e6_r8 / rho(i,k) * icldm(i,k))), 0.0_r8)!accum number
                           !coarse mode
                           !cld_num_coarse(i,k) = max((cld_num_coarse(i,k) + (so4_num_cr* 1e6_r8 / rho(i,k) * icldm(i,k))), 0.0_r8)
                           !num_coarse(i,k) = max((num_coarse(i,k) - (so4_num_cr* 1e6_r8 / rho(i,k) * icldm(i,k))), 0.0_r8)
                           cld_num_coarse_tend2 = max(((so4_num_cr* 1e6_r8 / rho(i,k) * icldm(i,k))), 0.0_r8)!cld coarse num
                           num_coarse_tend2 = -max(((so4_num_cr* 1e6_r8 / rho(i,k) * icldm(i,k))), 0.0_r8) 

                           ! transfer mass from interstial to cloudborne
                           !accumulation mode
                           !cld_accum_so4(i,k) = max((cld_accum_so4(i,k) &
                           !        +  accum_so4(i,k) *icldm(i,k) ), 0.0_r8)
                           !accum_so4(i,k) = max((accum_so4(i,k) &
                           !        -  accum_so4(i,k) *icldm(i,k) ), 0.0_r8)
                           cld_mass_accum_so4_tend2 = max(accum_so4(i,k) *icldm(i,k), 0.0_r8) 
                           mass_accum_so4_tend2 = -max(accum_so4(i,k) *icldm(i,k),0.0_r8) 
                           !coarse mode
                           !cld_coarse_so4(i,k) = max((cld_coarse_so4(i,k) &
                           !        + coarse_so4(i,k) *icldm(i,k)), 0.0_r8)
                           !coarse_so4(i,k) = max((coarse_so4(i,k) &
                           !        - coarse_so4(i,k) *icldm(i,k) ), 0.0_r8)
                           cld_mass_coarse_so4_tend2 = max(coarse_so4(i,k) *icldm(i,k), 0.0_r8)
                           mass_coarse_so4_tend2 = -max(coarse_so4(i,k) *icldm(i,k), 0.0_r8)
                          !ptend%q(i,k,cnum_coarse_idx) = -(dso4_num * icldm(i,k))/rho(i,k)/1e-6_r8/dtime
                          !cld_num_coarse(i,k)   = cld_num_coarse(i,k) + (dso4_num * icldm(i,k))/rho(i,k)/1e-6_r8
                          ! removal coarse so4 mass
                          !ptend%q(i,k,ccoarse_so4_idx) = - dso4_num / so4_num_cr * icldm(i,k) * coarse_so4(i,k) / dtime!
                          !cld_coarse_so4(i,k) = cld_coarse_dust(i,k) + dso4_num / so4_num_cr *icldm(i,k) * coarse_so4(i,k)
                    endif! end coarse mode removal  
                    if (mode_strat_coarse_idx > 0._r8) then !MAM5 case
                          !remove as fraction of so4_num_st_cr to so4_num_cr
                          !fso4_m3 = so4_num_cr / (so4_num_st_cr + so4_num_cr +) !fraction of mode 3 so4
                          !fso4_m5 = so4_num_st_cr / (so4_num_st_cr + so4_num_cr) !fraction of mode 5 so4
                          !strat coarse mode num
                           !cld_num_strat_coarse(i,k) = max((cld_num_strat_coarse(i,k) + (so4_num_st_cr* 1e6_r8 / rho(i,k) * icldm(i,k))), 0.0_r8)
                           !num_strcrs(i,k) = max((num_strcrs(i,k) - (so4_num_st_cr* 1e6_r8 / rho(i,k) * icldm(i,k))), 0.0_r8)
                           cld_num_strat_coarse_tend2 = max(((so4_num_st_cr* 1e6_r8 / rho(i,k) * icldm(i,k))), 0.0_r8)
                           num_strcrs_tend2 = -max(((so4_num_st_cr* 1e6_r8 / rho(i,k) * icldm(i,k))), 0.0_r8)
                          !strat coarse mass 
                          !cld_strat_coarse_so4(i,k) = max((cld_strat_coarse_so4(i,k) &
                          !         + strat_coarse_so4(i,k) *icldm(i,k)), 0.0_r8)
                          !strat_coarse_so4(i,k) = max((strat_coarse_so4(i,k) &
                          !         - strat_coarse_so4(i,k) *icldm(i,k) ), 0.0_r8)
                          cld_mass_strat_coarse_so4_tend2 = max(strat_coarse_so4(i,k) *icldm(i,k), 0.0_r8)
                          mass_strat_coarse_so4_tend2 = -max(strat_coarse_so4(i,k) *icldm(i,k), 0.0_r8)


                          !remove coarse mode so4
                          !write(iulog,*) 'kzm_cnum_coarse_idx', cnum_coarse_idx
                          !write(iulog,*) 'kzm_q_cnum_coarse_idx ', -(dso4_num * fso4_m3 * icldm(i,k))/rho(i,k)/1e-6_r8/dtime
                    endif!mode5 removal
                    endif ! aero_nucleation_removal
!kzm --
                    ! end of removal
                    naai(i,k) = naai(i,k) + dso4_num
                    nihf(i,k) = nihf(i,k) + dso4_num
                 end if !if oso4_num > 0
               end if !if nucleate_ice_strat > 0
               if (aero_nucleation_removal) then
                  if (1>2) then      
                  if (mode_coarse_idx > 0._r8) then 
                     !num_accum_tend2,cld_num_accum_tend2,num_coarse_tend2,cld_num_coarse_tend2
                     !cld_mass_accum_so4_tend2,mass_accum_so4_tend2,cld_mass_coarse_so4_tend2,mass_coarse_so4_tend2     
                     num_accum(i,k) = num_accum(i,k) + num_accum_tend1 + num_accum_tend2
                     cld_num_accum(i,k) = cld_num_accum(i,k) + cld_num_accum_tend1 + cld_num_accum_tend2

                     num_coarse(i,k) = num_coarse(i,k) + num_coarse_tend1 + num_coarse_tend2
                     cld_num_coarse(i,k) = cld_num_coarse(i,k) + cld_num_coarse_tend1 + cld_num_coarse_tend2
                     
                     !remove bc in accumulation
                     accum_bc(i,k) = accum_bc(i,k) + mass_accum_bc_tend1
                     cld_accum_bc(i,k) = cld_accum_bc(i,k) + cld_mass_accum_bc_tend1 

                     !remove dst in coarse
                     cld_coarse_dust(i,k) = cld_coarse_dust(i,k) + cld_mass_coarse_dust_tend1
                     coarse_dust(i,k) = coarse_dust(i,k) + mass_coarse_dust_tend1

                     !remove so4 in accum
                     cld_accum_so4(i,k) = cld_accum_so4(i,k) + cld_mass_accum_so4_tend2
                     accum_so4(i,k) = accum_so4(i,k) + mass_accum_so4_tend2
                           
                     !remove so4 in coarse
                     cld_coarse_so4(i,k) = cld_coarse_so4(i,k) + cld_mass_coarse_so4_tend2
                     coarse_so4(i,k) = coarse_so4(i,k) + mass_coarse_so4_tend2                   
                  endif
                  endif
                  if (mode_strat_coarse_idx > 0._r8 ) then
                     !cld_num_strat_coarse_tend2, num_strcrs_tend2,cld_mass_strat_coarse_so4_tend2,mass_strat_coarse_so4_tend2       
                     cld_num_strat_coarse(i,k) = cld_num_strat_coarse(i,k) + cld_num_strat_coarse_tend2
                     num_strcrs(i,k) = num_strcrs(i,k) + num_strcrs_tend2
                     !remove so4 in strat coarse
                     cld_strat_coarse_so4(i,k) = cld_strat_coarse_so4(i,k) + cld_mass_strat_coarse_so4_tend2
                     strat_coarse_so4(i,k) = strat_coarse_so4(i,k) + mass_strat_coarse_so4_tend2    
                  endif   
               endif
          

                   
            else
                                  ! This maintains backwards compatibility with the previous version.
              if (pmid(i,k) <= 12500._r8 .and. pmid(i,k) > 100._r8 .and. abs(state%lat(i)) >= 60._r8 * pi / 180._r8) then
                 ramp = 1._r8 - min(1._r8, max(0._r8, (pmid(i,k) - 10000._r8) / 2500._r8))

                 if (oso4_num > 0._r8) then
                    dso4_num = (max(oso4_num, ramp * nucleate_ice_strat * so4_num) - oso4_num) * 1e6_r8 / rho(i,k)
                    naai(i,k) = naai(i,k) + dso4_num
                    nihf(i,k) = nihf(i,k) + dso4_num
                 end if
              end if
            end if


!kzm --


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
