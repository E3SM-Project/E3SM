module hetfrz_classnuc_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for hetfrz_classnuc module.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver, begchunk, endchunk
use physconst,      only: rair, cpair, rh2o, rhoh2o, mwh2o, tmelt, pi
use constituents,   only: cnst_get_ind
use physics_types,  only: physics_state
use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,   only: phys_getopts, use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_idx, rad_cnst_get_spec_idx, &
                            rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num, rad_cnst_get_mode_props

use physics_buffer, only: pbuf_add_field, dtype_r8, pbuf_old_tim_idx, &
                          pbuf_get_index, pbuf_get_field
use cam_history,    only: addfld, phys_decomp, add_default, outfld

use ref_pres,       only: top_lev => trop_cloud_top_lev
use wv_saturation,  only: svp_water, svp_ice

use cam_logfile,    only: iulog
use error_messages, only: handle_errmsg, alloc_err
use cam_abortutils, only: endrun

use hetfrz_classnuc,   only: hetfrz_classnuc_init, hetfrz_classnuc_calc

implicit none
private
save

public :: &
   hetfrz_classnuc_cam_readnl,   &
   hetfrz_classnuc_cam_register, &
   hetfrz_classnuc_cam_init,     &
   hetfrz_classnuc_cam_calc,     &
   hetfrz_classnuc_cam_save_cbaero

! Namelist variables
logical :: hist_hetfrz_classnuc = .false.

! Vars set via init method.
real(r8) :: mincld      ! minimum allowed cloud fraction

! constituent indices
integer :: &
   cldliq_idx = -1, &
   cldice_idx = -1, &
   numliq_idx = -1, &
   numice_idx = -1

! pbuf indices for fields provided by heterogeneous freezing
integer :: &
   frzimm_idx, &
   frzcnt_idx, &
   frzdep_idx

! pbuf indices for fields needed by heterogeneous freezing
integer :: &
   ast_idx = -1

! modal aerosols
integer, parameter :: MAM3_nmodes = 3
integer, parameter :: MAM7_nmodes = 7
integer :: nmodes = -1             ! number of aerosol modes

! mode indices
integer :: mode_accum_idx  = -1    ! accumulation mode
integer :: mode_coarse_idx = -1    ! coarse mode
integer :: mode_finedust_idx = -1  ! fine dust mode
integer :: mode_coardust_idx = -1  ! coarse dust mode
integer :: mode_pcarbon_idx = -1   ! primary carbon mode

! mode properties
real(r8) :: alnsg_mode_accum
real(r8) :: alnsg_mode_coarse
real(r8) :: alnsg_mode_finedust
real(r8) :: alnsg_mode_coardust
real(r8) :: alnsg_mode_pcarbon

! specie properties
real(r8) :: specdens_dust
real(r8) :: specdens_so4
real(r8) :: specdens_bc
real(r8) :: specdens_soa
real(r8) :: specdens_pom

! List all species
integer :: ncnst = 0     ! Total number of constituents (mass and number) needed
                         ! by the parameterization (depends on aerosol model used)

integer :: so4_accum     ! sulfate in accumulation mode
integer :: bc_accum      ! black-c in accumulation mode
integer :: pom_accum     ! p-organic in accumulation mode
integer :: soa_accum     ! s-organic in accumulation mode
integer :: dst_accum     ! dust in accumulation mode
integer :: ncl_accum     ! seasalt in accumulation mode
integer :: num_accum     ! number in accumulation mode

integer :: dst_coarse    ! dust in coarse mode
integer :: ncl_coarse    ! seasalt in coarse mode
integer :: so4_coarse    ! sulfate in coarse mode
integer :: num_coarse    ! number in coarse mode

integer :: dst_finedust  ! dust in finedust mode
integer :: so4_finedust  ! sulfate in finedust mode
integer :: num_finedust  ! number in finedust mode

integer :: dst_coardust  ! dust in coardust mode
integer :: so4_coardust  ! sulfate in coardust mode
integer :: num_coardust  ! number in coardust mode

integer :: bc_pcarbon    ! black-c in primary carbon mode
integer :: pom_pcarbon   ! p-organic in primary carbon mode
integer :: num_pcarbon   ! number in primary carbon mode

! Index arrays for looping over all constituents
integer, allocatable :: mode_idx(:)
integer, allocatable :: spec_idx(:)

! Copy of cloud borne aerosols before modification by droplet nucleation
! The basis is converted from mass to volume.
real(r8), allocatable :: aer_cb(:,:,:,:)

! Copy of interstitial aerosols with basis converted from mass to volume.
real(r8), allocatable :: aer(:,:,:,:)

!===============================================================================
contains
!===============================================================================

subroutine hetfrz_classnuc_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'hetfrz_classnuc_cam_readnl'

  namelist /hetfrz_classnuc_nl/ hist_hetfrz_classnuc

  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'hetfrz_classnuc_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, hetfrz_classnuc_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

  end if

#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(hist_hetfrz_classnuc, 1, mpilog, 0, mpicom)
#endif

end subroutine hetfrz_classnuc_cam_readnl

!================================================================================================

subroutine hetfrz_classnuc_cam_register()

   if (.not. use_hetfrz_classnuc) return

   ! pbuf fields provided by hetfrz_classnuc
   call pbuf_add_field('FRZIMM', 'physpkg', dtype_r8, (/pcols,pver/), frzimm_idx)
   call pbuf_add_field('FRZCNT', 'physpkg', dtype_r8, (/pcols,pver/), frzcnt_idx)
   call pbuf_add_field('FRZDEP', 'physpkg', dtype_r8, (/pcols,pver/), frzdep_idx)

end subroutine hetfrz_classnuc_cam_register

!================================================================================================

subroutine hetfrz_classnuc_cam_init(mincld_in)

   real(r8), intent(in) :: mincld_in

   ! local variables
   logical  :: prog_modal_aero
   integer  :: m, n, nspec
   integer  :: istat

   real(r8) :: sigma_logr_aer

   character(len=32) :: str32
   character(len=*), parameter :: routine = 'hetfrz_classnuc_cam_init'
   !--------------------------------------------------------------------------------------------

   if (.not. use_hetfrz_classnuc) return

   ! This parameterization currently assumes that prognostic modal aerosols are on.  Check...
   call phys_getopts(prog_modal_aero_out=prog_modal_aero)
   if (.not. prog_modal_aero) call endrun(routine//': cannot use hetfrz_classnuc without prognostic modal aerosols')

   mincld = mincld_in

   call cnst_get_ind('CLDLIQ', cldliq_idx)
   call cnst_get_ind('CLDICE', cldice_idx)
   call cnst_get_ind('NUMLIQ', numliq_idx)
   call cnst_get_ind('NUMICE', numice_idx)

   ! pbuf fields used by hetfrz_classnuc
   ast_idx      = pbuf_get_index('AST')

   call addfld('bc_num', '#/cm3', pver, 'A', 'total bc number', phys_decomp)
   call addfld('dst1_num', '#/cm3', pver, 'A', 'total dst1 number', phys_decomp)
   call addfld('dst3_num', '#/cm3', pver, 'A', 'total dst3 number', phys_decomp)
   call addfld('bcc_num', '#/cm3', pver, 'A', 'coated bc number', phys_decomp)
   call addfld('dst1c_num', '#/cm3', pver, 'A', 'coated dst1 number', phys_decomp)
   call addfld('dst3c_num', '#/cm3', pver, 'A', 'coated dst3 number', phys_decomp)
   call addfld('bcuc_num', '#/cm3', pver, 'A', 'uncoated bc number', phys_decomp)
   call addfld('dst1uc_num', '#/cm3', pver, 'A', 'uncoated dst1 number', phys_decomp)
   call addfld('dst3uc_num', '#/cm3', pver, 'A', 'uncoated dst3 number', phys_decomp)

   call addfld('bc_a1_num', '#/cm3', pver, 'A', 'interstitial bc number', phys_decomp)
   call addfld('dst_a1_num', '#/cm3', pver, 'A', 'interstitial dst1 number', phys_decomp)
   call addfld('dst_a3_num', '#/cm3', pver, 'A', 'interstitial dst3 number', phys_decomp)
   call addfld('bc_c1_num', '#/cm3', pver, 'A', 'cloud borne bc number', phys_decomp)
   call addfld('dst_c1_num', '#/cm3', pver, 'A', 'cloud borne dst1 number', phys_decomp)
   call addfld('dst_c3_num', '#/cm3', pver, 'A', 'cloud borne dst3 number', phys_decomp)
    
   call addfld('fn_bc_c1_num', '#/cm3', pver, 'A', 'cloud borne bc number derived from fn', phys_decomp)
   call addfld('fn_dst_c1_num', '#/cm3', pver, 'A', 'cloud borne dst1 number derived from fn', phys_decomp)
   call addfld('fn_dst_c3_num', '#/cm3', pver, 'A', 'cloud borne dst3 number derived from fn', phys_decomp)

   call addfld('na500', '#/cm3', pver, 'A', 'interstitial aerosol number with D>500 nm', phys_decomp)
   call addfld('totna500', '#/cm3', pver, 'A', 'total aerosol number with D>500 nm', phys_decomp)

   call addfld('FREQIMM', 'fraction', pver, 'A', 'Fractional occurance of immersion  freezing', phys_decomp)
   call addfld('FREQCNT', 'fraction', pver, 'A', 'Fractional occurance of contact    freezing', phys_decomp)
   call addfld('FREQDEP', 'fraction', pver, 'A', 'Fractional occurance of deposition freezing', phys_decomp)
   call addfld('FREQMIX', 'fraction', pver, 'A', 'Fractional occurance of mixed-phase clouds' , phys_decomp)

   call addfld('DSTFREZIMM', 'm-3s-1', pver, 'A', 'dust immersion  freezing rate', phys_decomp)
   call addfld('DSTFREZCNT', 'm-3s-1', pver, 'A', 'dust contact    freezing rate', phys_decomp)
   call addfld('DSTFREZDEP', 'm-3s-1', pver, 'A', 'dust deposition freezing rate', phys_decomp)

   call addfld('BCFREZIMM', 'm-3s-1', pver, 'A', 'bc immersion  freezing rate', phys_decomp)
   call addfld('BCFREZCNT', 'm-3s-1', pver, 'A', 'bc contact    freezing rate', phys_decomp)
   call addfld('BCFREZDEP', 'm-3s-1', pver, 'A', 'bc deposition freezing rate', phys_decomp)

   call addfld('NIMIX_IMM', '#/m3', pver, 'A', &
               'Activated Ice Number Concentration due to het immersion freezing in Mixed Clouds', phys_decomp)
   call addfld('NIMIX_CNT', '#/m3', pver, 'A', &
               'Activated Ice Number Concentration due to het contact freezing in Mixed Clouds', phys_decomp)
   call addfld('NIMIX_DEP', '#/m3', pver, 'A', &
               'Activated Ice Number Concentration due to het deposition freezing in Mixed Clouds', phys_decomp)

   call addfld('DSTNIDEP', '#/m3', pver, 'A', &
               'Activated Ice Number Concentration due to dst dep freezing in Mixed Clouds', phys_decomp)
   call addfld('DSTNICNT', '#/m3', pver, 'A', &
               'Activated Ice Number Concentration due to dst cnt freezing in Mixed Clouds', phys_decomp)
   call addfld('DSTNIIMM', '#/m3', pver, 'A', &
               'Activated Ice Number Concentration due to dst imm freezing in Mixed Clouds', phys_decomp)

   call addfld('BCNIDEP', '#/m3', pver, 'A', &
               'Activated Ice Number Concentration due to bc dep freezing in Mixed Clouds', phys_decomp)
   call addfld('BCNICNT', '#/m3', pver, 'A', &
               'Activated Ice Number Concentration due to bc cnt freezing in Mixed Clouds', phys_decomp)
   call addfld('BCNIIMM', '#/m3', pver, 'A', &
               'Activated Ice Number Concentration due to bc imm freezing in Mixed Clouds', phys_decomp)

   call addfld('NUMICE10s', '#/m3', pver, 'A', &
               'Ice Number Concentration due to het freezing in Mixed Clouds during 10-s period', phys_decomp)
   call addfld('NUMIMM10sDST', '#/m3', pver, 'A', &
               'Ice Number Concentration due to imm freezing by dst in Mixed Clouds during 10-s period', phys_decomp)
   call addfld('NUMIMM10sBC', '#/m3', pver, 'A', &
               'Ice Number Concentration due to imm freezing by bc in Mixed Clouds during 10-s period', phys_decomp)

   if (hist_hetfrz_classnuc) then

      call add_default('bc_num', 1, ' ')
      call add_default('dst1_num', 1, ' ')
      call add_default('dst3_num', 1, ' ')
      call add_default('bcc_num', 1, ' ')
      call add_default('dst1c_num', 1, ' ')
      call add_default('dst3c_num', 1, ' ')
      call add_default('bcuc_num', 1, ' ')
      call add_default('dst1uc_num', 1, ' ')
      call add_default('dst3uc_num', 1, ' ')

      call add_default('bc_a1_num', 1, ' ')
      call add_default('dst_a1_num', 1, ' ')
      call add_default('dst_a3_num', 1, ' ')
      call add_default('bc_c1_num', 1, ' ')
      call add_default('dst_c1_num', 1, ' ')
      call add_default('dst_c3_num', 1, ' ')
    
      call add_default('fn_bc_c1_num', 1, ' ')
      call add_default('fn_dst_c1_num', 1, ' ')
      call add_default('fn_dst_c3_num', 1, ' ')

      call add_default('na500', 1, ' ')
      call add_default('totna500', 1, ' ')

      call add_default('FREQIMM', 1, ' ')
      call add_default('FREQCNT', 1, ' ')
      call add_default('FREQDEP', 1, ' ')
      call add_default('FREQMIX', 1, ' ')

      call add_default('DSTFREZIMM', 1, ' ')
      call add_default('DSTFREZCNT', 1, ' ')
      call add_default('DSTFREZDEP', 1, ' ')

      call add_default('BCFREZIMM', 1, ' ')
      call add_default('BCFREZCNT', 1, ' ')
      call add_default('BCFREZDEP', 1, ' ')

      call add_default('NIMIX_IMM', 1, ' ')
      call add_default('NIMIX_CNT', 1, ' ')  
      call add_default('NIMIX_DEP', 1, ' ')

      call add_default('DSTNIDEP', 1, ' ')
      call add_default('DSTNICNT', 1, ' ')
      call add_default('DSTNIIMM', 1, ' ')

      call add_default('BCNIDEP', 1, ' ')
      call add_default('BCNICNT', 1, ' ')
      call add_default('BCNIIMM', 1, ' ')

      call add_default('NUMICE10s', 1, ' ')
      call add_default('NUMIMM10sDST', 1, ' ')
      call add_default('NUMIMM10sBC', 1, ' ')

   end if

   ! The following code sets indices of the mode specific species used
   ! in the module.  Having a list of the species needed allows us to
   ! allocate temporary space for just those species rather than for all the
   ! CAM species (pcnst) which may be considerably more than needed.
   !
   ! The indices set below are for use with the CAM rad_constituents
   ! interfaces.  Using the rad_constituents interfaces isolates the physics
   ! parameterization which requires constituent information from the chemistry
   ! code which provides that information.

   ! nmodes is the total number of modes
   call rad_cnst_get_info(0, nmodes=nmodes)

   ! Determine mode indices for all modes referenced in this module.
   mode_accum_idx    = rad_cnst_get_mode_idx(0, 'accum')
   mode_coarse_idx   = rad_cnst_get_mode_idx(0, 'coarse')
   mode_finedust_idx = rad_cnst_get_mode_idx(0, 'fine_dust')
   mode_coardust_idx = rad_cnst_get_mode_idx(0, 'coarse_dust')
   mode_pcarbon_idx  = rad_cnst_get_mode_idx(0, 'primary_carbon')

   ! Check that required mode types were found
   if (nmodes == MAM3_nmodes) then
      if (mode_accum_idx == -1 .or. mode_coarse_idx == -1) then
         write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
            mode_accum_idx, mode_coarse_idx
         call endrun(routine//': ERROR required mode type not found')
      end if

   else if (nmodes == MAM7_nmodes) then
      if (mode_coardust_idx == -1 .or. mode_finedust_idx == -1 .or. mode_pcarbon_idx == -1) then
         write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
            mode_coardust_idx, mode_finedust_idx, mode_pcarbon_idx
         call endrun(routine//': ERROR required mode type not found')
      end if
   end if

   ! Set some mode properties

   call rad_cnst_get_mode_props(0, mode_accum_idx, sigmag=sigma_logr_aer)
   alnsg_mode_accum = log(sigma_logr_aer)
    
   if (nmodes == MAM3_nmodes) then
      call rad_cnst_get_mode_props(0, mode_coarse_idx, sigmag=sigma_logr_aer)
      alnsg_mode_coarse = log(sigma_logr_aer)

   else if (nmodes == MAM7_nmodes) then
      call rad_cnst_get_mode_props(0, mode_finedust_idx, sigmag=sigma_logr_aer)
      alnsg_mode_finedust = log(sigma_logr_aer)

      call rad_cnst_get_mode_props(0, mode_coardust_idx, sigmag=sigma_logr_aer)
      alnsg_mode_coardust = log(sigma_logr_aer)

      call rad_cnst_get_mode_props(0, mode_pcarbon_idx, sigmag=sigma_logr_aer)
      alnsg_mode_pcarbon = log(sigma_logr_aer)
   end if

   ! Set list indices for all constituents (mass and number) used in this module.
   ! The list is specific to the aerosol model used.  Note that the order of the 
   ! constituents in these lists is arbitrary.

   if (nmodes == MAM3_nmodes) then
      ncnst = 11
      so4_accum  =  1
      bc_accum   =  2
      pom_accum  =  3
      soa_accum  =  4
      dst_accum  =  5
      ncl_accum  =  6
      num_accum  =  7
      dst_coarse =  8
      ncl_coarse =  9
      so4_coarse = 10
      num_coarse = 11
   else if (nmodes == MAM7_nmodes) then
      ncnst = 15
      so4_accum    =  1
      bc_accum     =  2
      pom_accum    =  3
      soa_accum    =  4
      ncl_accum    =  6
      num_accum    =  7
      dst_finedust =  8
      so4_finedust =  9
      num_finedust = 10
      dst_coardust = 11
      so4_coardust = 12
      num_coardust = 13
      bc_pcarbon   = 5
      pom_pcarbon  = 14
      num_pcarbon  = 15
   end if

   ! Allocate arrays to hold specie and mode indices for all constitutents (mass and number) 
   ! needed in this module.
   allocate(mode_idx(ncnst), spec_idx(ncnst), stat=istat)
   call alloc_err(istat, routine, 'mode_idx, spec_idx', ncnst)
   mode_idx = -1
   spec_idx = -1

   ! Allocate space for copy of cloud borne aerosols before modification by droplet nucleation.
   allocate(aer_cb(pcols,pver,ncnst,begchunk:endchunk), stat=istat)
   call alloc_err(istat, routine, 'aer_cb', pcols*pver*ncnst*(endchunk-begchunk+1))

   ! Allocate space for copy of interstitial aerosols with modified basis
   allocate(aer(pcols,pver,ncnst,begchunk:endchunk), stat=istat)
   call alloc_err(istat, routine, 'aer', pcols*pver*ncnst*(endchunk-begchunk+1))

   ! The following code sets the species and mode indices for each constituent
   ! in the list.  The indices are identical in the interstitial and the cloud
   ! borne phases.
   ! Specie index 0 is used to indicate the mode number mixing ratio

   ! Indices for species in accumulation mode (so4, bc, pom, soa, nacl, dust)
   spec_idx(num_accum) = 0
   mode_idx(num_accum) = mode_accum_idx
   spec_idx(so4_accum) = rad_cnst_get_spec_idx(0, mode_accum_idx, 'sulfate')
   mode_idx(so4_accum) = mode_accum_idx
   spec_idx(bc_accum)  = rad_cnst_get_spec_idx(0, mode_accum_idx, 'black-c')
   mode_idx(bc_accum)  = mode_accum_idx
   spec_idx(pom_accum) = rad_cnst_get_spec_idx(0, mode_accum_idx, 'p-organic')
   mode_idx(pom_accum) = mode_accum_idx
   spec_idx(soa_accum) = rad_cnst_get_spec_idx(0, mode_accum_idx, 's-organic')
   mode_idx(soa_accum) = mode_accum_idx
   spec_idx(ncl_accum) = rad_cnst_get_spec_idx(0, mode_accum_idx, 'seasalt')
   mode_idx(ncl_accum) = mode_accum_idx
   if (nmodes == MAM3_nmodes) then
      spec_idx(dst_accum) = rad_cnst_get_spec_idx(0, mode_accum_idx, 'dust')
      mode_idx(dst_accum) = mode_accum_idx
   end if

   ! Indices for species in coarse mode (dust, nacl, so4)
   if (mode_coarse_idx > 0) then
      spec_idx(num_coarse) = 0
      mode_idx(num_coarse) = mode_coarse_idx
      spec_idx(ncl_coarse) = rad_cnst_get_spec_idx(0, mode_coarse_idx, 'seasalt')
      mode_idx(ncl_coarse) = mode_coarse_idx
      spec_idx(dst_coarse) = rad_cnst_get_spec_idx(0, mode_coarse_idx, 'dust')
      mode_idx(dst_coarse) = mode_coarse_idx
      spec_idx(so4_coarse) = rad_cnst_get_spec_idx(0, mode_coarse_idx, 'sulfate')
      mode_idx(so4_coarse) = mode_coarse_idx
   end if

   ! Indices for species in fine dust mode (dust, so4)
   if (mode_finedust_idx > 0) then
      spec_idx(num_finedust) = 0
      mode_idx(num_finedust) = mode_finedust_idx
      spec_idx(dst_finedust) = rad_cnst_get_spec_idx(0, mode_finedust_idx, 'dust')
      mode_idx(dst_finedust) = mode_finedust_idx
      spec_idx(so4_finedust) = rad_cnst_get_spec_idx(0, mode_finedust_idx, 'sulfate')
      mode_idx(so4_finedust) = mode_finedust_idx
   end if

   ! Indices for species in coarse dust mode (dust, so4)
   if (mode_coardust_idx > 0) then
      spec_idx(num_coardust) = 0
      mode_idx(num_coardust) = mode_coardust_idx
      spec_idx(dst_coardust) = rad_cnst_get_spec_idx(0, mode_coardust_idx, 'dust')
      mode_idx(dst_coardust) = mode_coardust_idx
      spec_idx(so4_coardust) = rad_cnst_get_spec_idx(0, mode_coardust_idx, 'sulfate')
      mode_idx(so4_coardust) = mode_coardust_idx
   end if

    ! Indices for species in primary carbon mode (bc, pom)
   if (mode_pcarbon_idx > 0) then
      spec_idx(num_pcarbon) = 0
      mode_idx(num_pcarbon) = mode_pcarbon_idx
      spec_idx(bc_pcarbon)  = rad_cnst_get_spec_idx(0, mode_pcarbon_idx, 'black-c')
      mode_idx(bc_pcarbon)  = mode_pcarbon_idx
      spec_idx(pom_pcarbon) = rad_cnst_get_spec_idx(0, mode_pcarbon_idx, 'p-organic')
      mode_idx(pom_pcarbon) = mode_pcarbon_idx
   end if
 
   ! Check that all required specie types were found
   if (any(spec_idx == -1)) then
      write(iulog,*) routine//': ERROR required species type not found - indicies:', spec_idx
      call endrun(routine//': ERROR required species type not found')
   end if

   ! Get some specie specific properties.
   if (nmodes == MAM3_nmodes) then
      call rad_cnst_get_aer_props(0, mode_idx(dst_accum), spec_idx(dst_accum), density_aer=specdens_dust)
   else if (nmodes == MAM7_nmodes) then
      call rad_cnst_get_aer_props(0, mode_idx(dst_finedust), spec_idx(dst_finedust), density_aer=specdens_dust)
   end if
   call rad_cnst_get_aer_props(0, mode_idx(so4_accum), spec_idx(so4_accum), density_aer=specdens_so4)
   call rad_cnst_get_aer_props(0, mode_idx(bc_accum),  spec_idx(bc_accum),  density_aer=specdens_bc)
   call rad_cnst_get_aer_props(0, mode_idx(soa_accum), spec_idx(soa_accum), density_aer=specdens_soa)
   call rad_cnst_get_aer_props(0, mode_idx(pom_accum), spec_idx(pom_accum), density_aer=specdens_pom)

   call hetfrz_classnuc_init( &
      rair, cpair, rh2o, rhoh2o, mwh2o, &
      tmelt, pi, iulog)

end subroutine hetfrz_classnuc_cam_init

!================================================================================================

subroutine hetfrz_classnuc_cam_calc( &
   state, deltatin, factnum, pbuf)

   ! arguments
   type(physics_state), target, intent(in)    :: state
   real(r8),                    intent(in)    :: deltatin       ! time step (s)
   real(r8),                    intent(in)    :: factnum(:,:,:) ! activation fraction for aerosol number
   type(physics_buffer_desc),   pointer       :: pbuf(:)
 
   ! local workspace

   ! outputs shared with the microphysics via the pbuf
   real(r8), pointer :: frzimm(:,:)
   real(r8), pointer :: frzcnt(:,:)
   real(r8), pointer :: frzdep(:,:)

   integer :: itim_old
   integer :: i, k

   real(r8) :: rho(pcols,pver)          ! air density (kg m-3)

   real(r8), pointer :: ast(:,:)        

   real(r8) :: lcldm(pcols,pver)

   real(r8), pointer :: ptr2d(:,:)

   real(r8) :: fn(3)
   real(r8) :: awcam(pcols,pver,3)
   real(r8) :: awfacm(pcols,pver,3)
   real(r8) :: hetraer(pcols,pver,3)
   real(r8) :: dstcoat(pcols,pver,3)
   real(r8) :: total_interstitial_aer_num(pcols,pver,3)
   real(r8) :: total_cloudborne_aer_num(pcols,pver,3)
   real(r8) :: total_aer_num(pcols,pver,3)
   real(r8) :: coated_aer_num(pcols,pver,3)
   real(r8) :: uncoated_aer_num(pcols,pver,3)

   real(r8) :: fn_cloudborne_aer_num(pcols,pver,3)


   real(r8) :: con1, r3lx, supersatice

   real(r8) :: qcic
   real(r8) :: ncic

   real(r8) :: frzbcimm(pcols,pver), frzduimm(pcols,pver)
   real(r8) :: frzbccnt(pcols,pver), frzducnt(pcols,pver)
   real(r8) :: frzbcdep(pcols,pver), frzdudep(pcols,pver)

   real(r8) :: freqimm(pcols,pver), freqcnt(pcols,pver), freqdep(pcols,pver), freqmix(pcols,pver)
   real(r8) :: nnuccc_bc(pcols,pver), nnucct_bc(pcols,pver), nnudep_bc(pcols,pver)
   real(r8) :: nnuccc_dst(pcols,pver), nnucct_dst(pcols,pver), nnudep_dst(pcols,pver)
   real(r8) :: niimm_bc(pcols,pver), nicnt_bc(pcols,pver), nidep_bc(pcols,pver)
   real(r8) :: niimm_dst(pcols,pver), nicnt_dst(pcols,pver), nidep_dst(pcols,pver)
   real(r8) :: numice10s(pcols,pver)
   real(r8) :: numice10s_imm_dst(pcols,pver)
   real(r8) :: numice10s_imm_bc(pcols,pver)

   real(r8) :: na500(pcols,pver)
   real(r8) :: tot_na500(pcols,pver)

   character(128) :: errstring   ! Error status
   !-------------------------------------------------------------------------------

   associate( &
      lchnk => state%lchnk,             &
      ncol  => state%ncol,              &
      t     => state%t,                 &
      qc    => state%q(:,:,cldliq_idx), &
      nc    => state%q(:,:,numliq_idx), &
      pmid  => state%pmid               )

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx, ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   do k = top_lev, pver
      do i = 1, ncol
         rho(i,k) = pmid(i,k)/(rair*t(i,k))
      end do
   end do

   do k = top_lev, pver
      do i = 1, ncol
         lcldm(i,k) = max(ast(i,k), mincld)
      end do
   end do

   ! Convert interstitial and cloud borne aerosols from a mass to a volume basis before
   ! being used in get_aer_num
   do i = 1, ncnst
      aer_cb(:ncol,:,i,lchnk) = aer_cb(:ncol,:,i,lchnk) * rho(:ncol,:)

      ! Check whether constituent is a mass or number mixing ratio
      if (spec_idx(i) == 0) then
         call rad_cnst_get_mode_num(0, mode_idx(i), 'a', state, pbuf, ptr2d)
      else
         call rad_cnst_get_aer_mmr(0, mode_idx(i), spec_idx(i), 'a', state, pbuf, ptr2d)
      end if
      aer(:ncol,:,i,lchnk) = ptr2d(:ncol,:) * rho(:ncol,:)
   end do

   ! Init top levels of outputs of get_aer_num
   total_aer_num              = 0._r8
   coated_aer_num             = 0._r8
   uncoated_aer_num           = 0._r8
   total_interstitial_aer_num = 0._r8
   total_cloudborne_aer_num   = 0._r8
   hetraer                    = 0._r8
   awcam                      = 0._r8
   awfacm                     = 0._r8
   dstcoat                    = 0._r8
   na500                      = 0._r8
   tot_na500                  = 0._r8

   ! output aerosols as reference information for heterogeneous freezing
   do i = 1, ncol
      do k = top_lev, pver
         call get_aer_num(i, k, ncnst, aer(:,:,:,lchnk), aer_cb(:,:,:,lchnk), rho(i,k), &
            total_aer_num(i,k,:), coated_aer_num(i,k,:), uncoated_aer_num(i,k,:),       &
            total_interstitial_aer_num(i,k,:), total_cloudborne_aer_num(i,k,:),         &
            hetraer(i,k,:), awcam(i,k,:), awfacm(i,k,:), dstcoat(i,k,:),                &
            na500(i,k), tot_na500(i,k))

         fn_cloudborne_aer_num(i,k,1) = total_aer_num(i,k,1)*factnum(i,k,mode_accum_idx)  ! bc
         if (nmodes == MAM3_nmodes) then
            fn_cloudborne_aer_num(i,k,2) = total_aer_num(i,k,2)*factnum(i,k,mode_accum_idx)  ! dst_a1
            fn_cloudborne_aer_num(i,k,3) = total_aer_num(i,k,3)*factnum(i,k,mode_coarse_idx) ! dst_a3
         else if (nmodes == MAM7_nmodes) then
            fn_cloudborne_aer_num(i,k,2) = total_aer_num(i,k,2)*factnum(i,k,mode_finedust_idx) 
            fn_cloudborne_aer_num(i,k,3) = total_aer_num(i,k,3)*factnum(i,k,mode_coardust_idx) 
         end if
      end do
   end do

   call outfld('bc_num',        total_aer_num(:,:,1),    pcols, lchnk)
   call outfld('dst1_num',      total_aer_num(:,:,2),    pcols, lchnk)
   call outfld('dst3_num',      total_aer_num(:,:,3),    pcols, lchnk)

   call outfld('bcc_num',       coated_aer_num(:,:,1),   pcols, lchnk)
   call outfld('dst1c_num',     coated_aer_num(:,:,2),   pcols, lchnk)
   call outfld('dst3c_num',     coated_aer_num(:,:,3),   pcols, lchnk)

   call outfld('bcuc_num',      uncoated_aer_num(:,:,1), pcols, lchnk)
   call outfld('dst1uc_num',    uncoated_aer_num(:,:,2), pcols, lchnk)
   call outfld('dst3uc_num',    uncoated_aer_num(:,:,3), pcols, lchnk)

   call outfld('bc_a1_num',     total_interstitial_aer_num(:,:,1), pcols, lchnk)
   call outfld('dst_a1_num',    total_interstitial_aer_num(:,:,2), pcols, lchnk)
   call outfld('dst_a3_num',    total_interstitial_aer_num(:,:,3), pcols, lchnk)

   call outfld('bc_c1_num',     total_cloudborne_aer_num(:,:,1),   pcols, lchnk)
   call outfld('dst_c1_num',    total_cloudborne_aer_num(:,:,2),   pcols, lchnk)
   call outfld('dst_c3_num',    total_cloudborne_aer_num(:,:,3),   pcols, lchnk)

   call outfld('fn_bc_c1_num',  fn_cloudborne_aer_num(:,:,1),      pcols, lchnk)
   call outfld('fn_dst_c1_num', fn_cloudborne_aer_num(:,:,2),      pcols, lchnk)
   call outfld('fn_dst_c3_num', fn_cloudborne_aer_num(:,:,3),      pcols, lchnk)
        
   call outfld('na500',         na500,     pcols, lchnk)
   call outfld('totna500',      tot_na500, pcols, lchnk)

   ! frzimm, frzcnt, frzdep are the outputs of this parameterization used by the microphysics
   call pbuf_get_field(pbuf, frzimm_idx, frzimm)
   call pbuf_get_field(pbuf, frzcnt_idx, frzcnt)
   call pbuf_get_field(pbuf, frzdep_idx, frzdep)
    
   frzimm(:ncol,:) = 0._r8
   frzcnt(:ncol,:) = 0._r8
   frzdep(:ncol,:) = 0._r8

   frzbcimm(:ncol,:) = 0._r8
   frzduimm(:ncol,:) = 0._r8
   frzbccnt(:ncol,:) = 0._r8
   frzducnt(:ncol,:) = 0._r8
   frzbcdep(:ncol,:) = 0._r8
   frzdudep(:ncol,:) = 0._r8

   freqimm(:ncol,:) = 0._r8
   freqcnt(:ncol,:) = 0._r8
   freqdep(:ncol,:) = 0._r8
   freqmix(:ncol,:) = 0._r8

   numice10s(:ncol,:)         = 0._r8
   numice10s_imm_dst(:ncol,:) = 0._r8
   numice10s_imm_bc(:ncol,:)  = 0._r8

   do i = 1, ncol
      do k = top_lev, pver

         if (t(i,k) > 235.15_r8 .and. t(i,k) < 269.15_r8) then
            qcic = min(qc(i,k)/lcldm(i,k), 5.e-3_r8)
            ncic = max(nc(i,k)/lcldm(i,k), 0._r8)

            con1 = 1._r8/(1.333_r8*pi)**0.333_r8
            r3lx = con1*(rho(i,k)*qcic/(rhoh2o*max(ncic*rho(i,k), 1.0e6_r8)))**0.333_r8 ! in m
            r3lx = max(4.e-6_r8, r3lx)
            supersatice = svp_water(t(i,k))/svp_ice(t(i,k))

            fn(1) = factnum(i,k,mode_accum_idx)  ! bc accumulation mode
            if (nmodes == MAM3_nmodes) then
                fn(2) = factnum(i,k,mode_accum_idx)  ! dust_a1 accumulation mode
                fn(3) = factnum(i,k,mode_coarse_idx) ! dust_a3 coarse mode
            else if (nmodes == MAM7_nmodes) then
                fn(2) = factnum(i,k,mode_finedust_idx)  
                fn(3) = factnum(i,k,mode_coardust_idx)   
            end if

            call hetfrz_classnuc_calc( &
               deltatin,  t(i,k),  pmid(i,k),  supersatice,   &
               fn,  r3lx,  ncic*rho(i,k)*1.0e-6_r8,  frzbcimm(i,k),  frzduimm(i,k),   &
               frzbccnt(i,k),  frzducnt(i,k),  frzbcdep(i,k),  frzdudep(i,k),  hetraer(i,k,:), &
               awcam(i,k,:), awfacm(i,k,:), dstcoat(i,k,:), total_aer_num(i,k,:),  &
               coated_aer_num(i,k,:), uncoated_aer_num(i,k,:), total_interstitial_aer_num(i,k,:), &
               total_cloudborne_aer_num(i,k,:), errstring)

            call handle_errmsg(errstring, subname="hetfrz_classnuc_calc")

            frzimm(i,k) = frzbcimm(i,k) + frzduimm(i,k)
            frzcnt(i,k) = frzbccnt(i,k) + frzducnt(i,k)
            frzdep(i,k) = frzbcdep(i,k) + frzdudep(i,k)

            if (frzimm(i,k) > 0._r8) freqimm(i,k) = 1._r8
            if (frzcnt(i,k) > 0._r8) freqcnt(i,k) = 1._r8
            if (frzdep(i,k) > 0._r8) freqdep(i,k) = 1._r8
            if ((frzimm(i,k) + frzcnt(i,k) + frzdep(i,k)) > 0._r8) freqmix(i,k) = 1._r8
         else
            frzimm(i,k) = 0._r8
            frzcnt(i,k) = 0._r8
            frzdep(i,k) = 0._r8
         end if

         nnuccc_bc(i,k) = frzbcimm(i,k)*1.0e6_r8*ast(i,k)
         nnucct_bc(i,k) = frzbccnt(i,k)*1.0e6_r8*ast(i,k)
         nnudep_bc(i,k) = frzbcdep(i,k)*1.0e6_r8*ast(i,k)

         nnuccc_dst(i,k) = frzduimm(i,k)*1.0e6_r8*ast(i,k)
         nnucct_dst(i,k) = frzducnt(i,k)*1.0e6_r8*ast(i,k)     
         nnudep_dst(i,k) = frzdudep(i,k)*1.0e6_r8*ast(i,k)

         niimm_bc(i,k) = frzbcimm(i,k)*1.0e6_r8*deltatin
         nicnt_bc(i,k) = frzbccnt(i,k)*1.0e6_r8*deltatin
         nidep_bc(i,k) = frzbcdep(i,k)*1.0e6_r8*deltatin

         niimm_dst(i,k) = frzduimm(i,k)*1.0e6_r8*deltatin
         nicnt_dst(i,k) = frzducnt(i,k)*1.0e6_r8*deltatin
         nidep_dst(i,k) = frzdudep(i,k)*1.0e6_r8*deltatin

         numice10s(i,k) = (frzimm(i,k)+frzcnt(i,k)+frzdep(i,k))*1.0e6_r8*deltatin*(10._r8/deltatin)
         numice10s_imm_dst(i,k) = frzduimm(i,k)*1.0e6_r8*deltatin*(10._r8/deltatin)
         numice10s_imm_bc(i,k) = frzbcimm(i,k)*1.0e6_r8*deltatin*(10._r8/deltatin)
      end do
   end do

   call outfld('FREQIMM', freqimm, pcols, lchnk)
   call outfld('FREQCNT', freqcnt, pcols, lchnk)
   call outfld('FREQDEP', freqdep, pcols, lchnk)
   call outfld('FREQMIX', freqmix, pcols, lchnk)

   call outfld('DSTFREZIMM', nnuccc_dst, pcols, lchnk)
   call outfld('DSTFREZCNT', nnucct_dst, pcols, lchnk)
   call outfld('DSTFREZDEP', nnudep_dst, pcols, lchnk)

   call outfld('BCFREZIMM', nnuccc_bc, pcols, lchnk)
   call outfld('BCFREZCNT', nnucct_bc, pcols, lchnk)
   call outfld('BCFREZDEP', nnudep_bc, pcols, lchnk)

   call outfld('NIMIX_IMM', niimm_bc+niimm_dst, pcols, lchnk)
   call outfld('NIMIX_CNT', nicnt_bc+nicnt_dst, pcols, lchnk)   
   call outfld('NIMIX_DEP', nidep_bc+nidep_dst, pcols, lchnk)

   call outfld('DSTNICNT', nicnt_dst, pcols, lchnk)
   call outfld('DSTNIDEP', nidep_dst, pcols, lchnk)
   call outfld('DSTNIIMM', niimm_dst, pcols, lchnk)

   call outfld('BCNICNT', nicnt_bc, pcols, lchnk)
   call outfld('BCNIDEP', nidep_bc, pcols, lchnk)
   call outfld('BCNIIMM', niimm_bc, pcols, lchnk)

   call outfld('NUMICE10s', numice10s, pcols, lchnk)
   call outfld('NUMIMM10sDST', numice10s_imm_dst, pcols, lchnk)
   call outfld('NUMIMM10sBC', numice10s_imm_bc, pcols, lchnk)

   end associate

end subroutine hetfrz_classnuc_cam_calc

!====================================================================================================

subroutine hetfrz_classnuc_cam_save_cbaero(state, pbuf)

   ! Save the required cloud borne aerosol constituents.
   type(physics_state),         intent(in)    :: state
   type(physics_buffer_desc),   pointer       :: pbuf(:)

   ! local variables
   integer :: i, lchnk
   real(r8), pointer :: ptr2d(:,:)
   !-------------------------------------------------------------------------------

   lchnk = state%lchnk

   ! loop over the cloud borne constituents required by this module and save
   ! a local copy

   do i = 1, ncnst

      ! Check whether constituent is a mass or number mixing ratio
      if (spec_idx(i) == 0) then
         call rad_cnst_get_mode_num(0, mode_idx(i), 'c', state, pbuf, ptr2d)
      else
         call rad_cnst_get_aer_mmr(0, mode_idx(i), spec_idx(i), 'c', state, pbuf, ptr2d)
      end if
      aer_cb(:,:,i,lchnk) = ptr2d
   end do

end subroutine hetfrz_classnuc_cam_save_cbaero

!====================================================================================================

subroutine get_aer_num(ii, kk, ncnst, aer, aer_cb, rhoair,&
                       total_aer_num,                     &
                       coated_aer_num,                    &
                       uncoated_aer_num,                  &
                       total_interstial_aer_num,          &
                       total_cloudborne_aer_num,          &
                       hetraer, awcam, awfacm, dstcoat,   &
                       na500, tot_na500)
    
   !*****************************************************************************
   ! Purpose: Calculate BC and Dust number, including total number(interstitial+
   !          cloud borne), one monolayer coated number, and uncoated number
   !
   ! Author: Yong Wang and Xiaohong Liu, UWyo, 12/2012
   !*****************************************************************************

   ! input
   integer,  intent(in) :: ii, kk, ncnst
   real(r8), intent(in) :: aer(pcols,pver,ncnst)    ! interstitial aerosols, volume basis
   real(r8), intent(in) :: aer_cb(pcols,pver,ncnst) ! cloud borne aerosols, volume basis
   real(r8), intent(in) :: rhoair                   ! air density (kg/m3)

   ! The interstitial and cloud borne aerosol concentrations are accessed from
   ! module variables local to this module.

   ! output
   real(r8), intent(out) :: total_aer_num(3)            ! #/cm^3
   real(r8), intent(out) :: total_interstial_aer_num(3) ! #/cm^3
   real(r8), intent(out) :: total_cloudborne_aer_num(3) ! #/cm^3
   real(r8), intent(out) :: coated_aer_num(3)           ! #/cm^3 
   real(r8), intent(out) :: uncoated_aer_num(3)         ! #/cm^3
   real(r8), intent(out) :: hetraer(3)                  ! BC and Dust mass mean radius [m]
   real(r8), intent(out) :: awcam(3)                    ! modal added mass [mug m-3]
   real(r8), intent(out) :: awfacm(3)                   ! (OC+BC)/(OC+BC+SO4)
   real(r8), intent(out) :: dstcoat(3)                  ! coated fraction
   real(r8), intent(out) :: na500                  ! #/cm^3 interstitial aerosol number with D>500 nm (#/cm^3)
   real(r8), intent(out) :: tot_na500              ! #/cm^3 total aerosol number with D>500 nm (#/cm^3)


   !local variables
   !------------coated variables--------------------
   real(r8), parameter :: n_so4_monolayers_dust = 1.0_r8 ! number of so4(+nh4) monolayers needed to coat a dust particle
   real(r8), parameter :: dr_so4_monolayers_dust = n_so4_monolayers_dust * 4.76e-10_r8
   real(r8), parameter :: spechygro_so4 = 0.507_r8          ! Sulfate hygroscopicity
   real(r8), parameter :: spechygro_soa = 0.14_r8           ! SOA hygroscopicity
   real(r8), parameter :: spechygro_pom = 0.1_r8            ! POM hygroscopicity
   real(r8), parameter :: soa_equivso4_factor = spechygro_soa/spechygro_so4
   real(r8), parameter :: pom_equivso4_factor = spechygro_pom/spechygro_so4
   real(r8) :: vol_shell(3)
   real(r8) :: vol_core(3) 
   real(r8) :: fac_volsfc_dust_a1, fac_volsfc_dust_a3, fac_volsfc_bc
   real(r8) :: tmp1, tmp2
   real(r8) :: bc_num                     ! bc number in accumulation mode for MAM3
                                          ! bc number in accumulation and primary carbon mode for MAM7
   real(r8) :: dst1_num, dst3_num         ! dust number in accumulation and corase mode for MAM3
                                          ! dust number in fine dust and corase dust mode for MAM7
   logical :: num_to_mass_in = .true.
   real(r8), parameter :: bc_num_to_mass = 4.669152e+17_r8      ! #/kg from emission
   real(r8), parameter :: dst1_num_to_mass = 3.484e+15_r8       ! #/kg for dust in accumulation mode

   real(r8) :: dmc, ssmc

   real(r8) :: as_so4, as_du, as_soa
   real(r8) :: dst1_num_imm, dst3_num_imm, bc_num_imm
   real(r8) :: dmc_imm, ssmc_imm
   real(r8) :: as_bc, as_pom, as_ss

   real(r8) :: r_bc                         ! model radii of BC modes [m]   
   real(r8) :: r_dust_a1, r_dust_a3         ! model radii of dust modes [m]   

   integer :: i
   real(r8) :: dst1_scale
   !-------------------------------------------------------------------------------

   ! init output vars
   total_aer_num             = 0._r8
   total_interstial_aer_num  = 0._r8
   total_cloudborne_aer_num  = 0._r8
   coated_aer_num            = 0._r8
   uncoated_aer_num          = 0._r8
   hetraer                   = 0._r8
   awcam                     = 0._r8
   awfacm                    = 0._r8
   dstcoat                   = 0._r8
   na500                     = 0._r8
   tot_na500                 = 0._r8

   !*****************************************************************************
   !                calculate intersitial aerosol 
   !*****************************************************************************

   if (nmodes == MAM3_nmodes) then

      if (.not. num_to_mass_in) then

         as_so4 = aer(ii,kk,so4_accum)
         as_bc  = aer(ii,kk,bc_accum)
         as_pom = aer(ii,kk,pom_accum)
         as_soa = aer(ii,kk,soa_accum)
         as_ss  = aer(ii,kk,ncl_accum)
         as_du  = aer(ii,kk,dst_accum)

         if (as_du > 0._r8) then
            dst1_num = as_du/(as_so4+as_bc+as_pom+as_soa+as_ss+as_du)  &
                       * aer(ii,kk,num_accum)*1.0e-6_r8 ! #/cm^3
         else
            dst1_num = 0.0_r8
         end if

         if (as_bc > 0._r8) then
            bc_num = as_bc/(as_so4+as_bc+as_pom+as_soa+as_ss+as_du)  &
                     * aer(ii,kk,num_accum)*1.0e-6_r8 ! #/cm^3
         else
            bc_num = 0.0_r8
         end if

      else

         dst1_num = aer(ii,kk,dst_accum) * dst1_num_to_mass*1.0e-6_r8 ! #/cm^3, dust # in accumulation mode
         bc_num = aer(ii,kk,bc_accum) * bc_num_to_mass*1.0e-6_r8      ! #/cm^3
      end if
      dmc = aer(ii,kk,dst_coarse)
      ssmc = aer(ii,kk,ncl_coarse)
      
      if (dmc > 0._r8 ) then
         dst3_num = dmc/(ssmc+dmc) * aer(ii,kk,num_coarse)*1.0e-6_r8 ! #/cm^3
      else
         dst3_num = 0.0_r8
      end if

   else if (nmodes == MAM7_nmodes) then
      bc_num = (aer(ii,kk,bc_accum)+aer(ii,kk,bc_pcarbon)) * bc_num_to_mass*1.0e-6_r8 ! #/cm^3
      dst1_num = aer(ii,kk,num_finedust)*1.0e-6_r8  ! #/cm^3
      dst3_num = aer(ii,kk,num_coardust)*1.0e-6_r8 ! #/cm^3    
   end if

   !*****************************************************************************
   !                calculate cloud borne aerosol 
   !*****************************************************************************    

   if (nmodes == MAM3_nmodes) then

      as_so4 = aer_cb(ii,kk,so4_accum)
      as_bc  = aer_cb(ii,kk,bc_accum)
      as_pom = aer_cb(ii,kk,pom_accum)
      as_soa = aer_cb(ii,kk,soa_accum)
      as_ss  = aer_cb(ii,kk,ncl_accum)
      as_du  = aer_cb(ii,kk,dst_accum)

      if (as_du > 0._r8) then
         dst1_num_imm = as_du/(as_so4+as_bc+as_pom+as_soa+as_ss+as_du)  &
                       * aer_cb(ii,kk,num_accum)*1.0e-6_r8 ! #/cm^3
      else
         dst1_num_imm = 0.0_r8
      end if

      if (as_bc > 0._r8) then
         bc_num_imm = as_bc/(as_so4+as_bc+as_pom+as_soa+as_ss+as_du)  &
                    * aer_cb(ii,kk,num_accum)*1.0e-6_r8 ! #/cm^3
      else
         bc_num_imm = 0.0_r8
      end if
       
      dmc_imm = aer_cb(ii,kk,dst_coarse)
      ssmc_imm = aer_cb(ii,kk,ncl_coarse)

      if (dmc_imm > 0._r8) then
         dst3_num_imm = dmc_imm/(ssmc_imm+dmc_imm) * aer_cb(ii,kk,num_coarse)*1.0e-6_r8 ! #/cm^3
      else
         dst3_num_imm = 0.0_r8
      end if

   else if (nmodes == MAM7_nmodes) then
      ! primary carbon mode is insoluble and thus don't consider its cloud-borne state
      as_so4 = aer_cb(ii,kk,so4_accum)
      as_bc  = aer_cb(ii,kk,bc_accum)
      as_pom = aer_cb(ii,kk,pom_accum)
      as_soa = aer_cb(ii,kk,soa_accum)
      as_ss  = aer_cb(ii,kk,ncl_accum)
      if (as_bc > 0._r8) then
         bc_num_imm = as_bc/(as_so4+as_bc+as_pom+as_soa+as_ss)  &
                    * aer_cb(ii,kk,num_accum)*1.0e-6_r8 ! #/cm^3
      else
         bc_num_imm = 0.0_r8
      end if
      dst1_num_imm = aer_cb(ii,kk,num_finedust)*1.0e-6_r8 ! #/cm^3
      dst3_num_imm = aer_cb(ii,kk,num_coardust)*1.0e-6_r8 ! #/cm^3
   end if

   total_interstial_aer_num(1) = bc_num
   total_interstial_aer_num(2) = dst1_num
   total_interstial_aer_num(3) = dst3_num 

   total_cloudborne_aer_num(1) = bc_num_imm 
   total_cloudborne_aer_num(2) = dst1_num_imm 
   total_cloudborne_aer_num(3) = dst3_num_imm
  
   !*****************************************************************************    
   ! calculate mass mean radius
   !*****************************************************************************    

   if (nmodes == MAM3_nmodes) then

      if (aer(ii,kk,bc_accum)*1.0e-3_r8 > 1.0e-30_r8 .and. bc_num > 1.0e-3_r8) then
         r_bc = ( 3._r8/(4*pi*specdens_bc)*aer(ii,kk,bc_accum)/(bc_num*1.0e6_r8) )**(1._r8/3._r8)
      else
         r_bc = 0.04e-6_r8
      end if

      if (aer(ii,kk,dst_accum)*1.0e-3_r8 > 1.0e-30_r8 .and. dst1_num > 1.0e-3_r8) then
         r_dust_a1 = ( 3._r8/(4*pi*specdens_dust)*aer(ii,kk,dst_accum)/(dst1_num*1.0e6_r8) )**(1._r8/3._r8)
      else
         r_dust_a1 = 0.258e-6_r8
      end if

      if (aer(ii,kk,dst_coarse)*1.0e-3_r8 > 1.0e-30_r8 .and. dst3_num > 1.0e-3_r8) then
         r_dust_a3 = ( 3._r8/(4*pi*specdens_dust)*aer(ii,kk,dst_coarse)/(dst3_num*1.0e6_r8) )**(1._r8/3._r8)
      else
         r_dust_a3 = 1.576e-6_r8
      end if

   else if (nmodes == MAM7_nmodes) then

      if ((aer(ii,kk,bc_accum)+aer(ii,kk,bc_pcarbon))*1.0e-3_r8 > 1.0e-30_r8 &
          .and. bc_num > 1.0e-3_r8) then
          r_bc = ( 3._r8/(4*pi*specdens_bc)*(aer(ii,kk,bc_accum)+aer(ii,kk,bc_pcarbon))/ &
                  (bc_num*1.0e6_r8) )**(1._r8/3._r8)
      else
          r_bc = 0.067e-6_r8 ! from emission size
      end if
      
      if (aer(ii,kk,dst_finedust)*1.0e-3_r8 > 1.0e-30_r8 .and. dst1_num > 1.0e-3_r8) then
         r_dust_a1 = ( 3._r8/(4*pi*specdens_dust)*aer(ii,kk,dst_finedust)/(dst1_num*1.0e6_r8) )**(1._r8/3._r8)
      else
         r_dust_a1 = 0.258e-6_r8
      end if

      if (aer(ii,kk,dst_coardust)*1.0e-3_r8 > 1.0e-30_r8 .and. dst3_num > 1.0e-3_r8) then
         r_dust_a3 = ( 3._r8/(4*pi*specdens_dust)*aer(ii,kk,dst_coardust)/(dst3_num*1.0e6_r8) )**(1._r8/3._r8)
      else
         r_dust_a3 = 1.576e-6_r8
      end if
   end if    

   hetraer(1) = r_bc
   hetraer(2) = r_dust_a1
   hetraer(3) = r_dust_a3

   !*****************************************************************************
   !                calculate coated fraction 
   !*****************************************************************************

   if (nmodes == MAM3_nmodes) then

      fac_volsfc_bc      = exp(2.5_r8*alnsg_mode_accum**2)
      fac_volsfc_dust_a1 = exp(2.5_r8*alnsg_mode_accum**2)
      fac_volsfc_dust_a3 = exp(2.5_r8*alnsg_mode_coarse**2)

      vol_shell(2) = ( aer(ii,kk,so4_accum)/specdens_so4 + &
                       aer(ii,kk,pom_accum)*pom_equivso4_factor/specdens_pom + &
                       aer(ii,kk,soa_accum)*soa_equivso4_factor/specdens_soa )/rhoair

      vol_core(2) = aer(ii,kk,dst_accum)/(specdens_dust*rhoair)

      !   ratio1 = vol_shell/vol_core = 
      !      actual hygroscopic-shell-volume/dust-core-volume
      !   ratio2 = 6.0_r8*dr_so4_monolayers_pcage/(dgncur_a*fac_volsfc_dust)
      !      = (shell-volume corresponding to n_so4_monolayers_pcage)/core-volume 
      !      The 6.0/(dgncur_a*fac_volsfc_dust) = (mode-surface-area/mode-volume)
      !   Note that vol_shell includes both so4, pom, AND soa as "equivalent so4",
      !      The soa_equivso4_factor accounts for the lower hygroscopicity of soa.
      !
      !   Define xferfrac_pcage = min( 1.0, ratio1/ratio2)
      !   But ratio1/ratio2 == tmp1/tmp2, and coding below avoids possible overflow 

      ! bc
      vol_shell(1) = vol_shell(2)
      vol_core(1) = aer(ii,kk,bc_accum)/(specdens_bc*rhoair)
      tmp1 = vol_shell(1)*(r_bc*2._r8)*fac_volsfc_bc
      tmp2 = max(6.0_r8*dr_so4_monolayers_dust*vol_core(1), 0.0_r8)
      dstcoat(1) = tmp1/tmp2

      ! dust_a1
      tmp1 = vol_shell(2)*(r_dust_a1*2._r8)*fac_volsfc_dust_a1
      tmp2 = max(6.0_r8*dr_so4_monolayers_dust*vol_core(2), 0.0_r8)
      dstcoat(2) = tmp1/tmp2

      ! dust_a3
      vol_shell(3) = aer(ii,kk,so4_coarse)/(specdens_so4*rhoair)
      vol_core(3)  = aer(ii,kk,dst_coarse)/(specdens_dust*rhoair)
      tmp1 = vol_shell(3)*(r_dust_a3*2._r8)*fac_volsfc_dust_a3
      tmp2 = max(6.0_r8*dr_so4_monolayers_dust*vol_core(3), 0.0_r8)
      dstcoat(3) = tmp1/tmp2

   else if (nmodes == MAM7_nmodes) then

      ! for BC, only consider primary carbon mode,
      ! because most of particles in this mode are uncoated
      ! and nearly all particles in accumulation mode are coated
      fac_volsfc_bc      = exp(2.5_r8*alnsg_mode_pcarbon**2)

      vol_shell(1) = ( aer(ii,kk,pom_pcarbon)*pom_equivso4_factor/specdens_pom )/rhoair
      vol_core(1)  = aer(ii,kk,bc_pcarbon)/(specdens_bc*rhoair)
      tmp1 = vol_shell(1)*(r_bc*2._r8)*fac_volsfc_bc
      tmp2 = max(6.0_r8*dr_so4_monolayers_dust*vol_core(1), 0.0_r8)
      dstcoat(1) = tmp1/tmp2

      fac_volsfc_dust_a1 = exp(2.5_r8*alnsg_mode_finedust**2)
      fac_volsfc_dust_a3 = exp(2.5_r8*alnsg_mode_coardust**2)

      vol_shell(2) = aer(ii,kk,so4_finedust)/(specdens_so4*rhoair)
      vol_core(2)  = aer(ii,kk,dst_finedust)/(specdens_dust*rhoair)

      tmp1 = vol_shell(2)*(r_dust_a1*2._r8)*fac_volsfc_dust_a1
      tmp2 = max(6.0_r8*dr_so4_monolayers_dust*vol_core(2), 0.0_r8)
      dstcoat(2) = tmp1/tmp2

      vol_shell(3) = aer(ii,kk,so4_coardust)/(specdens_so4*rhoair)
      vol_core(3)  = aer(ii,kk,dst_coardust)/(specdens_dust*rhoair)
      tmp1 = vol_shell(3)*(r_dust_a3*2._r8)*fac_volsfc_dust_a3 
      tmp2 = max(6.0_r8*dr_so4_monolayers_dust*vol_core(3), 0.0_r8)
      dstcoat(3) = tmp1/tmp2

   end if

   if (dstcoat(1) > 1._r8)    dstcoat(1) = 1._r8
   if (dstcoat(1) < 0.001_r8) dstcoat(1) = 0.001_r8
   if (dstcoat(2) > 1._r8)    dstcoat(2) = 1._r8
   if (dstcoat(2) < 0.001_r8) dstcoat(2) = 0.001_r8
   if (dstcoat(3) > 1._r8)    dstcoat(3) = 1._r8
   if (dstcoat(3) < 0.001_r8) dstcoat(3) = 0.001_r8
   
   do i = 1, 3
      total_aer_num(i)    = total_interstial_aer_num(i) + total_cloudborne_aer_num(i)
      coated_aer_num(i)   = total_interstial_aer_num(i)*dstcoat(i)
      uncoated_aer_num(i) = total_interstial_aer_num(i)*(1._r8-dstcoat(i))
   end do

   if (nmodes == MAM7_nmodes) then
      coated_aer_num(1)   = (aer(ii,kk,bc_pcarbon)*bc_num_to_mass*1.0e-6_r8)*dstcoat(1)+ &
                            (aer(ii,kk,bc_accum)*bc_num_to_mass*1.0e-6_r8)
      uncoated_aer_num(1) = (aer(ii,kk,bc_pcarbon)*bc_num_to_mass*1.0e-6_r8)*(1._r8-dstcoat(1))
   end if

   if (nmodes == MAM3_nmodes) then
      dst1_scale = 0.488_r8    ! scaled for D>0.5-1 um from 0.1-1 um
   else if (nmodes == MAM7_nmodes) then
      dst1_scale = 0.566_r8    ! scaled for D>0.5-2 um from 0.1-2 um
   end if

   tot_na500 = total_aer_num(1)*0.0256_r8          & ! scaled for D>0.5 um using Clarke et al., 1997; 2004; 2007: rg=0.1um, sig=1.6
      + total_aer_num(2)*dst1_scale + total_aer_num(3)

   na500 = total_interstial_aer_num(1)*0.0256_r8   & ! scaled for D>0.5 um using Clarke et al., 1997; 2004; 2007: rg=0.1um, sig=1.6
      + total_interstial_aer_num(2)*dst1_scale + total_interstial_aer_num(3)

   !*****************************************************************************
   !                prepare some variables for water activity 
   !*****************************************************************************

   if (nmodes == MAM3_nmodes) then

      ! accumulation mode for dust_a1 
      if (aer(ii,kk,num_accum) > 0._r8) then 
         awcam(2) = (dst1_num*1.0e6_r8)/aer(ii,kk,num_accum)* &
            ( aer(ii,kk,so4_accum) + aer(ii,kk,soa_accum) + &
            aer(ii,kk,pom_accum) + aer(ii,kk,bc_accum) )*1.0e9_r8 ! [mug m-3]
      else
         awcam(2) = 0._r8
      end if

      if (awcam(2) > 0._r8) then   
         awfacm(2) = ( aer(ii,kk,bc_accum) + aer(ii,kk,soa_accum) + aer(ii,kk,pom_accum) )/ &
            ( aer(ii,kk,soa_accum) + aer(ii,kk,pom_accum) + aer(ii,kk,so4_accum) + aer(ii,kk,bc_accum) )
      else
         awfacm(2) = 0._r8
      end if

      ! accumulation mode for bc
      if (aer(ii,kk,num_accum) > 0._r8) then
         awcam(1) = (bc_num*1.0e6_r8)/aer(ii,kk,num_accum)* &
            ( aer(ii,kk,so4_accum) + aer(ii,kk,soa_accum) + aer(ii,kk,pom_accum) + aer(ii,kk,bc_accum) )*1.0e9_r8 ! [mug m-3]
      else
         awcam(1) = 0._r8
      end if
      awfacm(1) = awfacm(2)

      ! coarse mode for dust_a3
      if (aer(ii,kk,num_coarse) > 0._r8) then
         awcam(3) = (dst3_num*1.0e6_r8)/aer(ii,kk,num_coarse)* aer(ii,kk,so4_coarse)*1.0e9_r8
      else
         awcam(3) = 0._r8
      end if
      awfacm(3) = 0._r8

   else if (nmodes == MAM7_nmodes) then

      ! accumulation mode for bc (primary carbon mode is insoluble)
      if (aer(ii,kk,num_accum) > 0._r8) then
         awcam(1) = (bc_num*1.0e6_r8)/aer(ii,kk,num_accum)* &
            ( aer(ii,kk,so4_accum) + aer(ii,kk,soa_accum) + aer(ii,kk,pom_accum) + aer(ii,kk,bc_accum) )*1.0e9_r8 ! [mug m-3]
      else
         awcam(1) = 0._r8
      end if

      if (awcam(1) > 0._r8) then
         awfacm(1) = ( aer(ii,kk,bc_accum) + aer(ii,kk,soa_accum) + aer(ii,kk,pom_accum) )/ &
            ( aer(ii,kk,soa_accum) + aer(ii,kk,pom_accum) + aer(ii,kk,so4_accum) + aer(ii,kk,bc_accum) )
      else
         awfacm(1) = 0._r8
      end if

      if (aer(ii,kk,num_finedust) > 0._r8) then
         awcam(2) = (dst1_num*1.0e6_r8)/aer(ii,kk,num_finedust)* aer(ii,kk,so4_finedust)*1.0e9_r8
      else
         awcam(2) = 0._r8
      end if
      awfacm(2) = 0._r8

      if (aer(ii,kk,num_coardust) > 0._r8) then
         awcam(3) = (dst3_num*1.0e6_r8)/aer(ii,kk,num_coardust)* aer(ii,kk,so4_coardust)*1.0e9_r8
      else 
         awcam(3) = 0._r8
      end if
      awfacm(3) = 0._r8

   end if
  
end subroutine get_aer_num

!====================================================================================================

end module hetfrz_classnuc_cam
