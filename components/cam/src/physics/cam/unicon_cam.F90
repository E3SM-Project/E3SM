!===========================================================================
! CAM interface to the UNIFIED CONVECTION SCHEME (UNICON)
!
! The USE_UNICON macro converts this module to a stub interface which allows
! CAM to be built without the unicon and unicon_utils modules.
!
!===========================================================================

module unicon_cam

use shr_kind_mod,     only: r8 => shr_kind_r8, i4 => shr_kind_i4
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols, pver, pverp, begchunk, endchunk
use physconst,        only: rair, cpair, gravit, latvap, latice, zvir, mwdry

use constituents,     only: pcnst, cnst_add, qmin, cnst_get_type_byind, cnst_get_ind, cnst_name
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_num_idx, rad_cnst_get_mam_mmr_idx
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
use camsrfexch,       only: cam_in_t
use physics_buffer,   only: pbuf_add_field, dtype_r8, dyn_time_lvls, pbuf_old_tim_idx, &
                            physics_buffer_desc, pbuf_get_index, pbuf_get_field, pbuf_set_field

use phys_control,     only: phys_getopts
use cam_history,      only: outfld, addfld, horiz_only, add_default

use time_manager,     only: is_first_step
use cam_abortutils,   only: endrun

#ifdef USE_UNICON
use unicon,           only: unicon_init, compute_unicon
use unicon_utils,     only: unicon_utils_init, positive_moisture, positive_tracer
#endif

implicit none
private
save

public :: &
   unicon_cam_readnl,      &
   unicon_cam_register,    &
   unicon_cam_init,        &
   unicon_implements_cnst, &
   unicon_init_cnst,       &
   unicon_out_t,           &
   unicon_cam_tend,        &
   unicon_cam_org_diags

! namelist variables
logical          :: unicon_offline_dat_out   = .false.
integer          :: unicon_offline_dat_hfile = 2

! properties
real(r8) :: xlv    ! Latent heat of vaporization
real(r8) :: xlf    ! Latent heat of fusion
real(r8) :: xls    ! Latent heat of sublimation
real(r8) :: cp     ! Specific heat of dry air

integer, parameter :: &
   nseg = 1,      &! Number of updraft segments [ # ]
   mix  = pcols,  &! Maximum number of columns
   mkx  = pver,   &! Number of vertical layers
   ncnst = pcnst   ! Number of advected constituents

! For advecting organization-related variables
integer, parameter :: n_org = 5                      ! Number of constituents
character(len=8), dimension(n_org), parameter :: &   ! Constituent names
   cnst_names = (/'ORGawk  ','ORGthl  ','ORGqto  ','ORGuoo  ','ORGvoo  '/)

integer :: awk_cnst_ind, thl_cnst_ind, qt_cnst_ind, u_cnst_ind, v_cnst_ind 

! fields added to physics buffer by this module
integer :: &
   cushavg_idx, &
   cuorg_idx, &
   awk_PBL_idx, &
   delta_thl_PBL_idx, &
   delta_qt_PBL_idx, &
   delta_u_PBL_idx, &
   delta_v_PBL_idx, &
   delta_tr_PBL_idx, &
   cu_cmfr_idx, &
   cu_thlr_idx, &
   cu_qtr_idx, &
   cu_ur_idx, &
   cu_vr_idx, &
   cu_qlr_idx, &
   cu_qir_idx, &
   cu_trr_idx, &
   cmfr_det_idx, &
   qlr_det_idx, &
   qir_det_idx,&
   rqcr_l_idx, &
   rqcr_i_idx, &
   rncr_l_idx, &
   rncr_i_idx, &
   rice2_idx

! fields expected to be in the physics buffer
integer :: &
   ast_idx   = -1,   &
   tke_idx   = -1,   &
   bprod_idx = -1, &
   kpblh_idx  = -1, &
   pblh_idx   = -1, &
   went_idx   = -1,   &
   cush_idx   = -1, &
   shfrc_idx  = -1, &
   icwmrsh_idx = -1, &
   rprdsh_idx  = -1, &
   prec_sh_idx = -1, &
   snow_sh_idx = -1, &
   nevapr_shcu_idx = -1, &
   am_evp_st_idx = -1, &    !  Evaporation area of stratiform precipitation [fraction]
   evprain_st_idx = -1, &   !  Grid-mean evaporation rate of stratiform rain [kg/kg/s] >= 0.
   evpsnow_st_idx = -1      !  Grid-mean evaporation rate of stratiform snow [kg/kg/s] >= 0.
   
! constituent indices
integer :: ixcldliq, ixcldice, ixnumliq, ixnumice

! unicon output fields
type unicon_out_t
   real(r8) :: cmfmc(mix,mkx+1)   ! Upward     convective mass flux at the interface [ kg / s / m2 ]
   real(r8) :: slflx(mix,mkx+1)   ! Net upward convective flux of liquid static energy [ J / s / m2 ]
   real(r8) :: qtflx(mix,mkx+1)   ! Net upward convective flux of total specific humidity [ kg / s / m2 ]
   real(r8) :: rqc(mix,mkx)       ! Prod of suspended LWC+IWC by expel of excessive in-cumulus condensate [ kg / kg / s ] > 0
   real(r8) :: rliq(mix)          ! Vertical integral of 'rqc' in flux unit [ m / s ]
   real(r8) :: cnt(mix)           ! Cloud top  interface index ( ki = kpen )
   real(r8) :: cnb(mix)           ! Cloud base interface index ( ki = krel-1 )
end type unicon_out_t

! logical array to identify constituents that are mode number concentrations
logical :: cnst_is_mam_num(ncnst)
! logical array to identify constituents that are mode specie mass mixing ratios
logical :: cnst_is_mam_mmr(ncnst)

!==================================================================================================
contains
!==================================================================================================
  
!> \brief Read namelist group unicon_nl
!!
!! \param[in] nlfile  ! filepath for file containing namelist input

subroutine unicon_cam_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'unicon_cam_readnl'
   
   namelist /unicon_nl/ unicon_offline_dat_out, unicon_offline_dat_hfile

   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'unicon_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, unicon_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast (unicon_offline_dat_out,   1,  mpilog, 0, mpicom)
   call mpibcast (unicon_offline_dat_hfile, 1,  mpiint, 0, mpicom)
#endif

end subroutine unicon_cam_readnl

!================================================================================================

subroutine unicon_cam_register

! Register fields in the constituent array and the physics buffer.


   ! Jun.02.2012. Sungsu for advecting organization-related horizontal heterogeneity
   !              within PBL. 
   !              For the time being, in order to save computation time, advection of aerosol perturbations 
   !              are simply neglected.

#ifdef USE_UNICON

   call cnst_add(cnst_names(1), mwdry, cpair,    0._r8, awk_cnst_ind, & 
      'Wake area within PBL associated with organization', readiv=.false., mixtype = 'dry')
   call cnst_add(cnst_names(2), mwdry, cpair,    0._r8, thl_cnst_ind, & 
      'Perturbation of  thl associated with organization', readiv=.false., mixtype = 'dry')
   call cnst_add(cnst_names(3), mwdry, cpair,    0._r8,  qt_cnst_ind, & 
      'Perturbation of  qt  associated with organization', readiv=.false., mixtype = 'dry')
   call cnst_add(cnst_names(4), mwdry, cpair,    0._r8,   u_cnst_ind, & 
      'Perturbation of  u   associated with organization', readiv=.false., mixtype = 'dry')
   call cnst_add(cnst_names(5), mwdry, cpair,    0._r8,   v_cnst_ind, & 
      'Perturbation of  v   associated with organization', readiv=.false., mixtype = 'dry')


   call pbuf_add_field('cushavg',       'global', dtype_r8, (/pcols,dyn_time_lvls/),            cushavg_idx)
   call pbuf_add_field('cuorg',         'global', dtype_r8, (/pcols,dyn_time_lvls/),            cuorg_idx)
   call pbuf_add_field('awk_PBL',       'global', dtype_r8, (/pcols,dyn_time_lvls/),            awk_PBL_idx)
   call pbuf_add_field('delta_thl_PBL', 'global', dtype_r8, (/pcols,dyn_time_lvls/),            delta_thl_PBL_idx)
   call pbuf_add_field('delta_qt_PBL',  'global', dtype_r8, (/pcols,dyn_time_lvls/),            delta_qt_PBL_idx)
   call pbuf_add_field('delta_u_PBL',   'global', dtype_r8, (/pcols,dyn_time_lvls/),            delta_u_PBL_idx)
   call pbuf_add_field('delta_v_PBL',   'global', dtype_r8, (/pcols,dyn_time_lvls/),            delta_v_PBL_idx)
   call pbuf_add_field('delta_tr_PBL',  'global', dtype_r8, (/pcols,pcnst,dyn_time_lvls/),      delta_tr_PBL_idx)
   call pbuf_add_field('cu_cmfr',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),       cu_cmfr_idx)
   call pbuf_add_field('cu_thlr',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),       cu_thlr_idx)
   call pbuf_add_field('cu_qtr',        'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),       cu_qtr_idx)
   call pbuf_add_field('cu_ur',         'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),       cu_ur_idx)
   call pbuf_add_field('cu_vr',         'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),       cu_vr_idx)
   call pbuf_add_field('cu_qlr',        'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),       cu_qlr_idx)
   call pbuf_add_field('cu_qir',        'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),       cu_qir_idx)
   call pbuf_add_field('cu_trr',        'global', dtype_r8, (/pcols,pver,pcnst,dyn_time_lvls/), cu_trr_idx)
   call pbuf_add_field('cmfr_det',      'global', dtype_r8, (/pcols,pver/),                     cmfr_det_idx)
   call pbuf_add_field('qlr_det',       'global', dtype_r8, (/pcols,pver/),                     qlr_det_idx)
   call pbuf_add_field('qir_det',       'global', dtype_r8, (/pcols,pver/),                     qir_det_idx)

   call pbuf_add_field('rqcr_l',        'global', dtype_r8, (/pcols,pver/),                      rqcr_l_idx)
   call pbuf_add_field('rqcr_i',        'global', dtype_r8, (/pcols,pver/),                      rqcr_i_idx)
   call pbuf_add_field('rncr_l',        'global', dtype_r8, (/pcols,pver/),                      rncr_l_idx)
   call pbuf_add_field('rncr_i',        'global', dtype_r8, (/pcols,pver/),                      rncr_i_idx)

   call pbuf_add_field('rice2',         'global', dtype_r8, (/pcols/),                            rice2_idx)

#endif

end subroutine unicon_cam_register

!==================================================================================================

subroutine unicon_cam_init(pbuf2d)

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local variables
   character(len=*), parameter :: sub='unicon_cam_init: '
   integer :: i, icnst, j, l, m, nmodes, nspec

   character(len=8)  :: units
   character(len=30) :: varname
   character(len=60) :: surname
   character(len=2)  :: numcha
   integer           :: msfc
   ! ------------------------------------------------------------------------------------------- !

#ifdef USE_UNICON

   ! constants
   xlv   = latvap
   xlf   = latice
   xls   = xlv + xlf
   cp    = cpair

   call unicon_init(latvap, cpair, latice, zvir, rair, gravit)

   ! save some constituent indices
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('NUMLIQ', ixnumliq)
   call cnst_get_ind('NUMICE', ixnumice)

   ! save some physics buffer indices
   ast_idx      = pbuf_get_index('AST')
   tke_idx      = pbuf_get_index('tke')
   bprod_idx    = pbuf_get_index('bprod')
   kpblh_idx    = pbuf_get_index('kpblh')
   pblh_idx    = pbuf_get_index('pblh')
   went_idx    = pbuf_get_index('went')
   cush_idx    = pbuf_get_index('cush')
   shfrc_idx   = pbuf_get_index('shfrc')
   icwmrsh_idx = pbuf_get_index('ICWMRSH')
   rprdsh_idx  = pbuf_get_index('RPRDSH')
   prec_sh_idx = pbuf_get_index('PREC_SH')
   snow_sh_idx = pbuf_get_index('SNOW_SH')
   nevapr_shcu_idx = pbuf_get_index('NEVAPR_SHCU')
   am_evp_st_idx  = pbuf_get_index('am_evp_st')
   evprain_st_idx = pbuf_get_index('evprain_st')
   evpsnow_st_idx = pbuf_get_index('evpsnow_st')

   ! physics buffer fields that need initializers -- these are only
   ! fields that have been added to pbuf by this module, or by the 
   ! convection driver layer.
   if (is_first_step()) then

      if (cush_idx > 0) then
         call pbuf_set_field(pbuf2d, cush_idx,          1.e3_r8)
      else
         call endrun(sub//'cush not in pbuf')
      end if
      if (cushavg_idx > 0) then
         call pbuf_set_field(pbuf2d, cushavg_idx,       1.e3_r8)
      else
         call endrun(sub//'cushavg not in pbuf')
      end if
      if (cuorg_idx > 0) then
         call pbuf_set_field(pbuf2d, cuorg_idx,         0.0_r8)
      else
         call endrun(sub//'cuorg not in pbuf')
      end if
      if (awk_pbl_idx > 0) then
         call pbuf_set_field(pbuf2d, awk_pbl_idx,       0.0_r8)
      else
         call endrun(sub//'awk_PBL not in pbuf')
      end if
      if (delta_thl_PBL_idx > 0) then
         call pbuf_set_field(pbuf2d, delta_thl_PBL_idx, 0.0_r8)
      else
         call endrun(sub//'delta_thl_PBL not in pbuf')
      end if
      if (delta_qt_PBL_idx > 0) then
         call pbuf_set_field(pbuf2d, delta_qt_PBL_idx,  0.0_r8)
      else
         call endrun(sub//'delta_qt_PBL not in pbuf')
      end if
      if (delta_u_PBL_idx > 0) then
         call pbuf_set_field(pbuf2d, delta_u_PBL_idx,   0.0_r8)
      else
         call endrun(sub//'delta_u_PBL not in pbuf')
      end if
      if (delta_v_PBL_idx > 0) then
         call pbuf_set_field(pbuf2d, delta_v_PBL_idx,   0.0_r8)
      else
         call endrun(sub//'delta_v_PBL not in pbuf')
      end if
      if (delta_tr_PBL_idx > 0) then
         call pbuf_set_field(pbuf2d, delta_tr_PBL_idx,  0.0_r8)
      else
         call endrun(sub//'delta_tr_PBL not in pbuf')
      end if
      if (cu_cmfr_idx > 0) then
         call pbuf_set_field(pbuf2d, cu_cmfr_idx,       0.0_r8)
      else
         call endrun(sub//'cu_cmfr not in pbuf')
      end if
      if (cu_thlr_idx > 0) then
         call pbuf_set_field(pbuf2d, cu_thlr_idx,       0.0_r8)
      else
         call endrun(sub//'cu_thlr not in pbuf')
      end if
      if (cu_qtr_idx > 0) then
         call pbuf_set_field(pbuf2d, cu_qtr_idx,        0.0_r8)
      else
         call endrun(sub//'cu_qtr not in pbuf')
      end if
      if (cu_ur_idx > 0) then
         call pbuf_set_field(pbuf2d, cu_ur_idx,         0.0_r8)
      else
         call endrun(sub//'cu_ur not in pbuf')
      end if
      if (cu_vr_idx > 0) then
         call pbuf_set_field(pbuf2d, cu_vr_idx,         0.0_r8)
      else
         call endrun(sub//'cu_vr not in pbuf')
      end if
      if (cu_qlr_idx > 0) then
         call pbuf_set_field(pbuf2d, cu_qlr_idx,        0.0_r8)
      else
         call endrun(sub//'cu_qlr not in pbuf')
      end if
      if (cu_qir_idx > 0) then
         call pbuf_set_field(pbuf2d, cu_qir_idx,        0.0_r8)
      else
         call endrun(sub//'cu_qir not in pbuf')
      end if
      if (cu_trr_idx > 0) then
         call pbuf_set_field(pbuf2d, cu_trr_idx,        0.0_r8)
      else
         call endrun(sub//'cu_trr not in pbuf')
      end if
      if (cmfr_det_idx > 0) then
         call pbuf_set_field(pbuf2d, cmfr_det_idx,      0.0_r8)
      else
         call endrun(sub//'cmfr_det not in pbuf')
      end if
      if (qlr_det_idx > 0) then
         call pbuf_set_field(pbuf2d, qlr_det_idx,      0.0_r8)
      else
         call endrun(sub//'qlr_det not in pbuf')
      end if
      if (qir_det_idx > 0) then
         call pbuf_set_field(pbuf2d, qir_det_idx,      0.0_r8)
      else
         call endrun(sub//'qir_det not in pbuf')
      end if
      if (rqcr_l_idx > 0) then
         call pbuf_set_field(pbuf2d, rqcr_l_idx,      0.0_r8)
      else
         call endrun(sub//'rqcr_l not in pbuf')
      end if
      if (rqcr_i_idx > 0) then
         call pbuf_set_field(pbuf2d, rqcr_i_idx,      0.0_r8)
      else
         call endrun(sub//'rqcr_i not in pbuf')
      end if
      if (rncr_l_idx > 0) then
         call pbuf_set_field(pbuf2d, rncr_l_idx,      0.0_r8)
      else
         call endrun(sub//'rncr_l not in pbuf')
      end if
      if (rncr_i_idx > 0) then
         call pbuf_set_field(pbuf2d, rncr_i_idx,      0.0_r8)
      else
         call endrun(sub//'rncr_i not in pbuf')
      end if
      if (rice2_idx > 0) then
         call pbuf_set_field(pbuf2d, rice2_idx,      0.0_r8)
      else
         call endrun(sub//'rice2 not in pbuf')
      end if

   end if

   ! Set arrays to identify the modal aerosol constituents
   cnst_is_mam_num = .false.
   cnst_is_mam_mmr = .false.
   call rad_cnst_get_info(0, nmodes=nmodes)
   do i = 1, nmodes
      call rad_cnst_get_mode_num_idx(i, icnst)
      cnst_is_mam_num(icnst) = .true.
      call rad_cnst_get_info(0, i, nspec=nspec)
      do j = 1, nspec
         call rad_cnst_get_mam_mmr_idx(i, j, icnst)
         cnst_is_mam_mmr(icnst) = .true.
      end do
   end do

   ! ------------------------- !
   ! Internal Output Variables !
   ! ------------------------- !

   ! Sungsu for advection of convective organization

   call addfld('ORGawk', (/ 'lev' /),  'A', 'no', 'Convective Organization - Wake Area' )
   call addfld('ORGthl', (/ 'lev' /),  'A', 'K', 'Convective Organization - Perturbation of thl in the non-wake area' )
   call addfld('ORGqto', (/ 'lev' /),  'A', 'kg/kg', 'Convective Organization - Perturbation of qt  in the non-wake area' )
   call addfld('ORGuoo', (/ 'lev' /),  'A', 'm/s', 'Convective Organization - Perturbation of u  in the non-wake area' )
   call addfld('ORGvoo', (/ 'lev' /),  'A', 'm/s', 'Convective Organization - Perturbation of v  in the non-wake area' )

   ! From the main unified convection scheme

   call addfld( 'flxrain_SP' , (/ 'ilev' /), 'A', 'kg/m/s2', 'Convective net rain flux' )
   call addfld( 'flxsnow_SP' , (/ 'ilev' /), 'A', 'kg/m/s2', 'Convective net snow flux' )

   call addfld('cmf_SP', (/ 'ilev' /), 'A',   'kg/m2/s', 'Convective net mass flux')
   call addfld('qtflx_SP', (/ 'ilev' /), 'A', 'kg/m2/s', 'Convective net qt flux')
   call addfld('slflx_SP', (/ 'ilev' /), 'A', 'J/m2/s' , 'Convective net sl flux')
   call addfld('uflx_SP', (/ 'ilev' /), 'A',  'kg/m/s2', 'Convective net u flux')
   call addfld('vflx_SP', (/ 'ilev' /), 'A',  'kg/m/s2', 'Convective net v flux')

   call addfld('cmf_u_SP', (/ 'ilev' /), 'A',   'kg/m2/s', 'Convective updraft mass flux')
   call addfld('qtflx_u_SP', (/ 'ilev' /), 'A', 'kg/m2/s', 'Convective updraft qt flux')
   call addfld('slflx_u_SP', (/ 'ilev' /), 'A', 'J/m2/s' , 'Convective updraft sl flux')
   call addfld('uflx_u_SP', (/ 'ilev' /), 'A',  'kg/m/s2', 'Convective updraft u flux')
   call addfld('vflx_u_SP', (/ 'ilev' /), 'A',  'kg/m/s2', 'Convective updraft v flux')

   call addfld('cmf_d_SP', (/ 'ilev' /), 'A',   'kg/m2/s', 'Convective downdraft mass flux')
   call addfld('qtflx_d_SP', (/ 'ilev' /), 'A', 'kg/m2/s', 'Convective downdraft qt flux')
   call addfld('slflx_d_SP', (/ 'ilev' /), 'A', 'J/m2/s' , 'Convective downdraft sl flux')
   call addfld('uflx_d_SP', (/ 'ilev' /), 'A',  'kg/m/s2', 'Convective downdraft u flux')
   call addfld('vflx_d_SP', (/ 'ilev' /), 'A',  'kg/m/s2', 'Convective downdraft v flux')

   call addfld('qtten_u_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qt by updraft')
   call addfld('slten_u_SP',  (/ 'lev' /), 'A', 'J/kg/s', 'Convective tendency sl by updraft')
   call addfld('uten_u_SP',   (/ 'lev' /), 'A',  'm/s/s', 'Convective tendency u  by updraft')
   call addfld('vten_u_SP',  (/ 'lev' /), 'A',  'm/s/s' , 'Convective tendency v  by updraft')
   call addfld('sten_u_SP',  (/ 'lev' /), 'A',  'J/kg/s', 'Convective tendency s  by updraft')
   call addfld('qvten_u_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qv by updraft')
   call addfld('qlten_u_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency ql by updraft')
   call addfld('qiten_u_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qi by updraft')
   call addfld('nlten_u_SP',  (/ 'lev' /), 'A', '1/kg/s', 'Convective tendency nl by updraft')
   call addfld('niten_u_SP',  (/ 'lev' /), 'A', '1/kg/s', 'Convective tendency ni by updraft')

   call addfld('qtten_d_SP', (/ 'lev' /), 'A',  'kg/kg/s', 'Convective tendency qt by downdraft')
   call addfld('slten_d_SP',  (/ 'lev' /), 'A',  'J/kg/s', 'Convective tendency sl by downdraft')
   call addfld('uten_d_SP',   (/ 'lev' /), 'A',   'm/s/s', 'Convective tendency u  by downdraft')
   call addfld('vten_d_SP',   (/ 'lev' /), 'A',   'm/s/s', 'Convective tendency v  by downdraft')
   call addfld('sten_d_SP',  (/ 'lev' /), 'A',   'J/kg/s', 'Convective tendency s  by downdraft')
   call addfld('qvten_d_SP', (/ 'lev' /), 'A',  'kg/kg/s', 'Convective tendency qv by downdraft')
   call addfld('qlten_d_SP', (/ 'lev' /), 'A',  'kg/kg/s', 'Convective tendency ql by downdraft')
   call addfld('qiten_d_SP', (/ 'lev' /), 'A',  'kg/kg/s', 'Convective tendency qi by downdraft')
   call addfld('nlten_d_SP',  (/ 'lev' /), 'A',  '1/kg/s', 'Convective tendency nl by downdraft')
   call addfld('niten_d_SP',  (/ 'lev' /), 'A',  '1/kg/s', 'Convective tendency ni by downdraft')

   call addfld('qtten_evp_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qt by evap of precip within environment')
   call addfld('slten_evp_SP',  (/ 'lev' /), 'A', 'J/kg/s', 'Convective tendency sl by evap of precip within environment')
   call addfld('uten_evp_SP',   (/ 'lev' /), 'A',  'm/s/s', 'Convective tendency u  by evap of precip within environment')
   call addfld('vten_evp_SP',   (/ 'lev' /), 'A',  'm/s/s', 'Convective tendency v  by evap of precip within environment')
   call addfld('sten_evp_SP',  (/ 'lev' /), 'A',  'J/kg/s', 'Convective tendency s  by evap of precip within environment')
   call addfld('qvten_evp_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qv by evap of precip within environment')
   call addfld('qlten_evp_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency ql by evap of precip within environment')
   call addfld('qiten_evp_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qi by evap of precip within environment')
   call addfld('nlten_evp_SP',  (/ 'lev' /), 'A', '#/kg/s', 'Convective tendency nl by evap of precip within environment')
   call addfld('niten_evp_SP',  (/ 'lev' /), 'A', '#/kg/s', 'Convective tendency ni by evap of precip within environment')

   call addfld('qtten_dis_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qt by dissipative heating')
   call addfld('slten_dis_SP',  (/ 'lev' /), 'A', 'J/kg/s', 'Convective tendency sl by dissipative heating')
   call addfld('uten_dis_SP',   (/ 'lev' /), 'A',  'm/s/s', 'Convective tendency u  by dissipative heating')
   call addfld('vten_dis_SP',   (/ 'lev' /), 'A',  'm/s/s', 'Convective tendency v  by dissipative heating')
   call addfld('sten_dis_SP',  (/ 'lev' /), 'A',  'J/kg/s', 'Convective tendency s  by dissipative heating')
   call addfld('qvten_dis_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qv by dissipative heating')
   call addfld('qlten_dis_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency ql by dissipative heating')
   call addfld('qiten_dis_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qi by dissipative heating')
   call addfld('nlten_dis_SP',  (/ 'lev' /), 'A', '1/kg/s', 'Convective tendency nl by dissipative heating')
   call addfld('niten_dis_SP',  (/ 'lev' /), 'A', '1/kg/s', 'Convective tendency ni by dissipative heating')

   call addfld('qtten_pos_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qt by positive tracer constraints')
   call addfld('slten_pos_SP',  (/ 'lev' /), 'A', 'J/kg/s', 'Convective tendency sl by positive tracer constraints')
   call addfld('uten_pos_SP',   (/ 'lev' /), 'A',  'm/s/s', 'Convective tendency u  by positive tracer constraints')
   call addfld('vten_pos_SP',   (/ 'lev' /), 'A',  'm/s/s', 'Convective tendency v  by positive tracer constraints')
   call addfld('sten_pos_SP',  (/ 'lev' /), 'A',  'J/kg/s', 'Convective tendency s  by positive tracer constraints')
   call addfld('qvten_pos_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qv by positive tracer constraints')
   call addfld('qlten_pos_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency ql by positive tracer constraints')
   call addfld('qiten_pos_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qi by positive tracer constraints')
   call addfld('nlten_pos_SP',  (/ 'lev' /), 'A', '1/kg/s', 'Convective tendency nl by positive tracer constraints')
   call addfld('niten_pos_SP',  (/ 'lev' /), 'A', '1/kg/s', 'Convective tendency ni by positive tracer constraints')

   call addfld('qlten_sub_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency ql by compensating subsidence')
   call addfld('qiten_sub_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qi by compensating subsidence')

   call addfld('qlten_det_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency ql by detrainment')
   call addfld('qiten_det_SP', (/ 'lev' /), 'A', 'kg/kg/s', 'Convective tendency qi by detrainment')

   !f  call addfld('exit_Cu_SP', '1', 1, 'A', 'Exit identifier of UNICON')
   call addfld('cush_SP', horiz_only, 'A',    'm', 'Cumulus top height from UNICON')
   call addfld('cushavg_SP', horiz_only, 'A', 'm', 'Mean cumulus top height from UNICON')
   call addfld('cuorg_SP', horiz_only, 'A',   '1', 'Cumulus organization parameter from UNICON')
   call addfld('Radius_SP', horiz_only, 'A',  'm', 'Cumulus plume radius at surface from UNICON')
   !d  call addfld('orgforce_SP', 'kg/m/s^2', 1, 'A', 'Various organization forcing of UNICON')
   !d  call addfld('tau_org_SP',  's',        1, 'A', 'Damping time scale of convective organization of UNICON')
   !d  call addfld('tau_TKE_SP',  's',        1, 'A', 'Damping time scale of meso-scale TKE of UNICON')
   call addfld('sgh_SP', horiz_only, 'A',   'm', 'Standard deviation of subgrid topography')
   call addfld('sgh30_SP', horiz_only, 'A', 'm', 'Standard deviation of subgrid topography at 30 sec')

   call addfld('CMFR_DET', (/ 'lev' /), 'A', 'kg/m2/s', 'Detrained convective mass flux')
   call addfld('QLR_DET',   (/ 'lev' /), 'A',  'kg/kg', 'Detrained convective LWC')
   call addfld('QIR_DET',   (/ 'lev' /), 'A',  'kg/kg', 'Detrained convective IWC')

   call addfld('kw_SP', horiz_only, 'A', 'no', 'Internally computed kw from SPARKCONV' )

   call addfld('sigma_w_SP',  horiz_only, 'A',   'm/s', 'Standard deviation of updraft w at surface from UNICON')
   call addfld('sigma_thl_SP',    horiz_only, 'A', 'K', 'Standard deviation of updraft thl at surface from UNICON')
   call addfld('sigma_qt_SP', horiz_only, 'A',  'g/kg', 'Standard deviation of updraft qt at surface from UNICON')
   call addfld('sigma_u_SP',  horiz_only, 'A',   'm/s', 'Standard deviation of updraft u at surface from UNICON')
   call addfld('sigma_v_SP',  horiz_only, 'A',   'm/s', 'Standard deviation of updraft v at surface from UNICON')

   call addfld('w_org_SP', horiz_only, 'A',   'm2/s2', 'Organization-generated additional w at surface from UNICON')
   call addfld('thl_org_SP',     horiz_only, 'A', 'K', 'Organization-generated additional thl at surface from UNICON')
   call addfld('qt_org_SP',  horiz_only, 'A',  'g/kg', 'Organization-generated additional qt at surface from UNICON')
   call addfld('u_org_SP',   horiz_only, 'A',   'm/s', 'Organization-generated additional u at surface from UNICON')
   call addfld('v_org_SP',   horiz_only, 'A',   'm/s', 'Organization-generated additional v at surface from UNICON')

   call addfld('tkes_SP', horiz_only, 'A',     'm2/s2', 'TKE at surface within UNICON')
   call addfld('went_SP',   horiz_only, 'A',     'm/s', 'Entrainment rate at the PBL top interface from UW PBL')
   call addfld('went_eff_SP',   horiz_only, 'A', 'm/s', 'Effective entrainment rate at the PBL top interface in UNICON')

   call addfld('am_u_SP', (/ 'lev' /), 'A', '1', 'Convective updraft fractional area at mid-layer')
   call addfld('am_d_SP', (/ 'lev' /), 'A', '1', 'Convective downdraft fractional area at mid-layer')

   call addfld('qlm_u_SP', (/ 'lev' /), 'A', 'kg/kg', 'Area-weighted in-cumulus LWC condensate of updraft at mid-layer')
   call addfld('qlm_d_SP', (/ 'lev' /), 'A', 'kg/kg', 'Area-weighted in-cumulus IWC condensate of downdraft at mid-layer')

   call addfld('qim_u_SP', (/ 'lev' /), 'A', 'kg/kg', 'Area-weighted in-cumulus LWC condensate of updraft at mid-layer')
   call addfld('qim_d_SP', (/ 'lev' /), 'A', 'kg/kg', 'Area-weighted in-cumulus IWC condensate of downdraft at mid-layer')

   call addfld('thl_u_SP',     (/ 'ilev' /), 'A',  'K', 'Mass-flux weighted updraft mean thl')
   call addfld('qt_u_SP', (/ 'ilev' /), 'A',   'kg/kg', 'Mass-flux weighted updraft mean qt')
   call addfld('u_u_SP',   (/ 'ilev' /), 'A',    'm/s', 'Mass-flux weighted updraft mean u')
   call addfld('v_u_SP',   (/ 'ilev' /), 'A',    'm/s', 'Mass-flux weighted updraft mean v')
   call addfld('w_u_SP',   (/ 'ilev' /), 'A',    'm/s', 'Mass-flux weighted updraft mean w')
   call addfld('ql_u_SP', (/ 'ilev' /), 'A',   'kg/kg', 'Mass-flux weighted updraft mean ql')
   call addfld('qi_u_SP', (/ 'ilev' /), 'A',   'kg/kg', 'Mass-flux weighted updraft mean qi')
   call addfld('wa_u_SP',   (/ 'ilev' /), 'A',   'm/s', 'Area weighted updraft mean w')
   call addfld('qla_u_SP', (/ 'ilev' /), 'A',  'kg/kg', 'Area weighted updraft mean ql')
   call addfld('qia_u_SP', (/ 'ilev' /), 'A',  'kg/kg', 'Area weighted updraft mean qi')
   call addfld('a_u_SP',     (/ 'ilev' /), 'A',    '1', 'Convective updraft fractional area')
   call addfld('rad_u_SP',     (/ 'ilev' /), 'A',  'm', 'Number-weighted effective radius of updraft plumes')
   call addfld('num_u_SP', (/ 'ilev' /), 'A',  '1/m^2', 'Number concentration of updraft plumes')
   call addfld('gamw_u_SP', (/ 'ilev' /), 'A', 'ratio', 'The ratio of w_u to wa_u')
   call addfld('nl_u_SP',  (/ 'ilev' /), 'A',   '1/kg', 'Mass-flux weighted updraft mean nl')
   call addfld('ni_u_SP',  (/ 'ilev' /), 'A',   '1/kg', 'Mass-flux weighted updraft mean ni')
   call addfld('thva_u_SP',     (/ 'ilev' /), 'A', 'K', 'Area weighted updraft mean thv')

   call addfld('thl_d_SP',     (/ 'ilev' /), 'A', 'K', 'Mass-flux weighted downdraft mean thl')
   call addfld('qt_d_SP', (/ 'ilev' /), 'A',  'kg/kg', 'Mass-flux weighted downdraft mean qt')
   call addfld('u_d_SP',   (/ 'ilev' /), 'A',   'm/s', 'Mass-flux weighted downdraft mean u')
   call addfld('v_d_SP',   (/ 'ilev' /), 'A',   'm/s', 'Mass-flux weighted downdraft mean v')
   call addfld('w_d_SP',   (/ 'ilev' /), 'A',   'm/s', 'Mass-flux weighted downdraft mean w')
   call addfld('ql_d_SP', (/ 'ilev' /), 'A',  'kg/kg', 'Mass-flux weighted downdraft mean ql')
   call addfld('qi_d_SP', (/ 'ilev' /), 'A',  'kg/kg', 'Mass-flux weighted downdraft mean qi')
   call addfld('wa_d_SP' ,   (/ 'ilev' /), 'A', 'm/s', 'Area weighted downdraft mean w')
   call addfld('qla_d_SP', (/ 'ilev' /), 'A', 'kg/kg', 'Area weighted downdraft mean ql')
   call addfld('qia_d_SP', (/ 'ilev' /), 'A', 'kg/kg', 'Area weighted downdraft mean qi')
   call addfld('a_d_SP',     (/ 'ilev' /), 'A',   '1', 'Convective downdraft fractional area')
   call addfld('nl_d_SP',  (/ 'ilev' /), 'A',  '1/kg', 'Mass-flux weighted downdraft mean nl')
   call addfld('ni_d_SP',  (/ 'ilev' /), 'A',  '1/kg', 'Mass-flux weighted downdraft mean ni')

   call addfld('thv_b_SP', (/ 'ilev' /), 'A',   'K', 'thv_b : Environmental buoyancy at the base interface')
   call addfld('thv_t_SP', (/ 'ilev' /), 'A',   'K', 'thv_t : Environmental buoyancy at the top interface')
   call addfld('thv_mt_SP', (/ 'ilev' /), 'A',  'K', 'thv_mt : Environmental buoyancy at the top interface of lower layer')
   call addfld('thv_min_SP', (/ 'ilev' /), 'A', 'K', 'thv_min : Minimum environmental buoyancy for downdraft sorting')

   !a  call addfld('CFL_SP', '1', pver, 'A', 'Numerical stability parameter of UNICON')

   call addfld('cu_cmfr_SP', (/ 'lev' /), 'A', 'kg/m2/s', 'Mass flux of mixing environmental airs')
   call addfld('cu_thlr_SP',       (/ 'lev' /), 'A', 'K', 'thl       of mixing environmental airs')
   call addfld('cu_qtr_SP',   (/ 'lev' /), 'A',  'kg/kg', 'qt        of mixing environmental airs')
   call addfld('cu_qlr_SP',   (/ 'lev' /), 'A',  'kg/kg', 'ql        of mixing environmental airs')
   call addfld('cu_qir_SP',   (/ 'lev' /), 'A',  'kg/kg', 'qi        of mixing environmental airs')
   call addfld('cu_ur_SP',     (/ 'lev' /), 'A',   'm/s', 'u         of mixing environmental airs')
   call addfld('cu_vr_SP',     (/ 'lev' /), 'A',   'm/s', 'v         of mixing environmental airs')
   call addfld('cu_thvr_SP',       (/ 'lev' /), 'A', 'K', 'thv       of mixing environmental airs')
   call addfld('cu_rhr_SP',       (/ 'lev' /), 'A',  '1', 'rh        of mixing environmental airs')
   call addfld('cu_nlr_SP',    (/ 'lev' /), 'A',  '1/kg', 'nl        of mixing environmental airs')
   call addfld('cu_nir_SP',    (/ 'lev' /), 'A',  '1/kg', 'ni        of mixing environmental airs')

   call addfld( 'a_p_SP'    , (/ 'ilev' /), 'A', 'fraction', 'Convective precipitation area' &
    )
   call addfld( 'am_evp_SP' , (/ 'lev' /),  'A', 'fraction', 'Convective evaporation area'   &
    )
   call addfld( 'am_pu_SP'  , (/ 'lev' /),  'A', 'fraction', 'Overlapping area between conv precipitation and sat updraft area'    &
    )
   call addfld( 'x_um_SP'   ,        (/ 'lev' /),  'A', 'm', 'Zonal displacement of the updraft area from the surface'             &
    )
   call addfld( 'y_um_SP'   ,        (/ 'lev' /),  'A', 'm', 'Meridional displacement of the updraft area from the surface'        &
    )
   call addfld( 'x_p_SP'    ,        (/ 'ilev' /), 'A', 'm', 'Zonal displacement of the precipitation area from the surface'       &
    )
   call addfld( 'y_p_SP'    ,        (/ 'ilev' /), 'A', 'm', 'Meridional displacement of the precipitation area from the surface'  &
    )

   do l = 1, pcnst

      if (cnst_is_mam_num(l) .or. cnst_is_mam_mmr(l)) then

         units = '1/kg/s'
         if (cnst_is_mam_mmr(l)) units = 'kg/kg/s'

         varname = trim(cnst_name(l))//'_u_SP'
         surname = trim(cnst_name(l))//' tendency by convective updraft from UNICON'
         call addfld(trim(varname), (/ 'lev' /), 'A', units, trim(surname))

         varname = trim(cnst_name(l))//'_d_SP'
         surname = trim(cnst_name(l))//' tendency by convective downdraft from UNICON'
         call addfld(trim(varname), (/ 'lev' /), 'A', units, trim(surname))

         varname = trim(cnst_name(l))//'_evp_SP'
         surname = trim(cnst_name(l))//' tendency by evap. of precip in env. from UNICON'
         call addfld(trim(varname), (/ 'lev' /), 'A', units, trim(surname))

         varname = trim(cnst_name(l))//'_dis_SP'
         surname = trim(cnst_name(l))//' tendency by dissipative heating from UNICON'
         call addfld(trim(varname), (/ 'lev' /), 'A', units, trim(surname))

         varname = trim(cnst_name(l))//'_pos_SP'
         surname = trim(cnst_name(l))//' tendency by positive moisture from UNICON'
         call addfld(trim(varname), (/ 'lev' /), 'A', units, trim(surname))

      end if
   end do

   ! Nov.15.2012. Below output corresponding to individual updraft segment is designed to write out individual 
   !              segment values for writing UNICON_II paper.

   do msfc = 1, nseg
      write(numcha,'(i2.2)') msfc

      ! The properties of individual updraft segment       

      call addfld('thl_u'//numcha//'_SP',        (/ 'ilev' /), 'A', 'K', numcha//' updraft segment : updraft thl'            )
      call addfld('qt_u'//numcha//'_SP',    (/ 'ilev' /), 'A',  'kg/kg', numcha//' updraft segment : updraft qt'             )
      call addfld('u_u'//numcha//'_SP',      (/ 'ilev' /), 'A',   'm/s', numcha//' updraft segment : updraft u'              )
      call addfld('v_u'//numcha//'_SP',      (/ 'ilev' /), 'A',   'm/s', numcha//' updraft segment : updraft v'              )
      call addfld('w_u'//numcha//'_SP',      (/ 'ilev' /), 'A',   'm/s', numcha//' updraft segment : updraft w'              )
      call addfld('ql_u'//numcha//'_SP',    (/ 'ilev' /), 'A',  'kg/kg', numcha//' updraft segment : updraft ql'             )
      call addfld('qi_u'//numcha//'_SP',    (/ 'ilev' /), 'A',  'kg/kg', numcha//' updraft segment : updraft qi'             )
      call addfld('cmf_u'//numcha//'_SP', (/ 'ilev' /), 'A', 'kg/s/m^2', numcha//' updraft segment : updraft cmf'            )
      call addfld('a_u'//numcha//'_SP',        (/ 'ilev' /), 'A',   '1', numcha//' updraft segment : updraft fractional area')
      call addfld('num_u'//numcha//'_SP',    (/ 'ilev' /), 'A', '1/m^2', numcha//' updraft segment : updraft number density' )
      call addfld('rad_u'//numcha//'_SP',        (/ 'ilev' /), 'A', 'm', numcha//' updraft segment : updraft plume radius'   )
      call addfld('nl_u'//numcha//'_SP',     (/ 'ilev' /), 'A',  '1/kg', numcha//' updraft segment : updraft nl'             )
      call addfld('ni_u'//numcha//'_SP',     (/ 'ilev' /), 'A',  '1/kg', numcha//' updraft segment : updraft ni'             )

      call addfld('eps0_u'//numcha//'_SP',     (/ 'ilev' /), 'A',     '1/Pa', numcha//' updraft segment : updraft eps0'    )
      call addfld('eps_u'//numcha//'_SP' ,     (/ 'ilev' /), 'A',     '1/Pa', numcha//' updraft segment : updraft eps'     )
      call addfld('del_u'//numcha//'_SP' ,     (/ 'ilev' /), 'A',     '1/Pa', numcha//' updraft segment : updraft del'     )
      call addfld('eeps_u'//numcha//'_SP',        (/ 'ilev' /), 'A',     '1', numcha//' updraft segment : updraft eeps'    )
      call addfld('ddel_u'//numcha//'_SP',        (/ 'ilev' /), 'A',     '1', numcha//' updraft segment : updraft ddel'    )
      call addfld('xc_u'//numcha//'_SP',        (/ 'ilev' /), 'A',       '1', numcha//' updraft segment : updraft xc'      )
      call addfld('xs_u'//numcha//'_SP',        (/ 'ilev' /), 'A',       '1', numcha//' updraft segment : updraft xs'      )
      call addfld('xemin_u'//numcha//'_SP',        (/ 'ilev' /), 'A',    '1', numcha//' updraft segment : updraft xemin'   )
      call addfld('xemax_u'//numcha//'_SP',        (/ 'ilev' /), 'A',    '1', numcha//' updraft segment : updraft xemax'   )
      call addfld('cridis_u'//numcha//'_SP',        (/ 'ilev' /), 'A',   'm', numcha//' updraft segment : updraft cridis'  )
      call addfld('thvcuenv_u'//numcha//'_SP',        (/ 'ilev' /), 'A', 'K', numcha//' updraft segment : updraft thvcuenv')
      call addfld('thvegenv_u'//numcha//'_SP',        (/ 'ilev' /), 'A', 'K', numcha//' updraft segment : updraft thvegenv')
      call addfld('thvxsenv_u'//numcha//'_SP',        (/ 'ilev' /), 'A', 'K', numcha//' updraft segment : updraft thvxsenv')
      call addfld('fmix_u'//numcha//'_SP',        (/ 'ilev' /), 'A',     '1', numcha//' updraft segment : updraft fmix'    )
      call addfld('cmfumix_u'//numcha//'_SP', (/ 'ilev' /), 'A',  'kg/s/m^2', numcha//' updraft segment : updraft cmfumix' )

      ! call addfld('ktop'//numcha//'_SP', '1', 1, 'A', numcha//' updraft segment : top layer index')
      call addfld('ptop'//numcha//'_SP', horiz_only, 'A', 'Pa', numcha//' updraft segment : updraft top pressure')
      call addfld('ztop'//numcha//'_SP',  horiz_only, 'A', 'm', numcha//' updraft segment : updraft top height')

      ! The properties of mass flux weighted ( or area-weighted or net=sum ) downdraft properties for individual updraft segment       

      call addfld('thl_d'//numcha//'_SP',         (/ 'ilev' /), 'A', 'K',&
         'Mass-flux weighted  mean downdraft thl for '// numcha//' updraft segment')
      call addfld('qt_d'//numcha//'_SP' ,     (/ 'ilev' /), 'A', 'kg/kg',&
         'Mass-flux weighted  mean downdraft  qt for '// numcha//' updraft segment')
      call addfld('u_d'//numcha//'_SP'  ,       (/ 'ilev' /), 'A', 'm/s',&
         'Mass-flux weighted  mean downdraft   u for '// numcha//' updraft segment')
      call addfld('v_d'//numcha//'_SP'  ,       (/ 'ilev' /), 'A', 'm/s',&
         'Mass-flux weighted  mean downdraft   v for '// numcha//' updraft segment')
      call addfld('w_d'//numcha//'_SP'  ,       (/ 'ilev' /), 'A', 'm/s',&
         'Mass-flux weighted  mean downdraft   w for '// numcha//' updraft segment')
      call addfld('ql_d'//numcha//'_SP' ,     (/ 'ilev' /), 'A', 'kg/kg',&
         'Mass-flux weighted  mean downdraft  ql for '// numcha//' updraft segment')
      call addfld('qi_d'//numcha//'_SP' ,     (/ 'ilev' /), 'A', 'kg/kg',&
         'Mass-flux weighted  mean downdraft  qi for '// numcha//' updraft segment')
      call addfld('wa_d'//numcha//'_SP' ,       (/ 'ilev' /), 'A', 'm/s',&
         'Area      weighted  mean downdraft   w for '// numcha//' updraft segment')
      call addfld('qla_d'//numcha//'_SP',     (/ 'ilev' /), 'A', 'kg/kg',&
         'Area      weighted  mean downdraft  ql for '// numcha//' updraft segment')
      call addfld('qia_d'//numcha//'_SP',     (/ 'ilev' /), 'A', 'kg/kg',&
         'Area      weighted  mean downdraft  qi for '// numcha//' updraft segment')
      call addfld('cmf_d'//numcha//'_SP',  (/ 'ilev' /), 'A', 'kg/s/m^2',&
         'Net                      downdraft cmf for '// numcha//' updraft segment')
      call addfld('a_d'//numcha//'_SP'  ,  (/ 'ilev' /), 'A', 'fraction',&
         'Net                      downdraft   a for '// numcha//' updraft segment')
      call addfld('nl_d'//numcha//'_SP' ,      (/ 'ilev' /), 'A', '#/kg',&
         'Mass-flux weighted  mean downdraft  nl for '// numcha//' updraft segment')
      call addfld('ni_d'//numcha//'_SP' ,      (/ 'ilev' /), 'A', '#/kg',&
         'Mass-flux weighted  mean downdraft  ni for '// numcha//' updraft segment')

   enddo

   ! Nov.16.2012. Additional detailed diagnostic output

   call addfld('thl_orgfce_SP',      horiz_only, 'A', 'K/s', &
      'Total   organization forcing generating thl difference between non-wake and grid-mean areas')
   call addfld('qt_orgfce_SP',  horiz_only, 'A',  'kg/kg/s', &
      'Total   organization forcing generating qt  difference between non-wake and grid-mean areas')
   call addfld('u_orgfce_SP',    horiz_only, 'A',   'm/s/s', &
      'Total   organization forcing generating u   difference between non-wake and grid-mean areas')
   call addfld('v_orgfce_SP',    horiz_only, 'A',   'm/s/s', &
      'Total   organization forcing generating v   difference between non-wake and grid-mean areas')
   call addfld('awk_orgfce_SP',      horiz_only, 'A', '1/s', &
      'Total   organization forcing generating wake area')

   call addfld('thl_orgfce_f_SP',      horiz_only, 'A', 'K/s', &
      'PBL top flux-related forcing generating thl difference between non-wake and grid-mean areas')
   call addfld('qt_orgfce_f_SP',  horiz_only, 'A',  'kg/kg/s', &
      'PBL top flux-related forcing generating qt  difference between non-wake and grid-mean areas')
   call addfld('u_orgfce_f_SP',    horiz_only, 'A',   'm/s/s', &
      'PBL top flux-related forcing generating u   difference between non-wake and grid-mean areas')
   call addfld('v_orgfce_f_SP',    horiz_only, 'A',   'm/s/s', &
      'PBL top flux-related forcing generating v   difference between non-wake and grid-mean areas')
   call addfld('awk_orgfce_f_SP',      horiz_only, 'A', '1/s', &
      'PBL top flux-related forcing generating wake area')

   call addfld('thl_orgfce_u_SP',     horiz_only, 'A', 'K/s', &
      'Up-and-Down diabatic forcing generating thl difference between non-wake and grid-mean areas')
   call addfld('qt_orgfce_u_SP', horiz_only, 'A',  'kg/kg/s', &
      'Up-and-Down diabatic forcing generating qt  difference between non-wake and grid-mean areas')
   call addfld('u_orgfce_u_SP',   horiz_only, 'A',   'm/s/s', &
      'Up-and-Down diabatic forcing generating u   difference between non-wake and grid-mean areas')
   call addfld('v_orgfce_u_SP',   horiz_only, 'A',   'm/s/s', &
      'Up-and-Down diabatic forcing generating v   difference between non-wake and grid-mean areas')
   call addfld('awk_orgfce_m_SP',     horiz_only, 'A', '1/s', &
      'Lateral-Mixing       forcing for wake area')

   call addfld('thl_orgfce_e_SP',      horiz_only, 'A', 'K/s', &
      'Environment diabatic forcing generating thl difference between non-wake and grid-mean areas')
   call addfld('qt_orgfce_e_SP',  horiz_only, 'A',  'kg/kg/s', &
      'Environment diabatic forcing generating qt  difference between non-wake and grid-mean areas')
   call addfld('u_orgfce_e_SP',    horiz_only, 'A',   'm/s/s', &
      'Environment diabatic forcing generating u   difference between non-wake and grid-mean areas')
   call addfld('v_orgfce_e_SP',    horiz_only, 'A',   'm/s/s', &
      'Environment diabatic forcing generating v   difference between non-wake and grid-mean areas')
   call addfld('cmf_d_orgh_SP', horiz_only, 'A',   'kg/m^2/s', &
      'Organization-inducing downdraft mass flux at the PBL top interface')

   call addfld('taui_thl_SP', horiz_only, 'A', '1/s', &
      'Inverse of damping time scale of the difference between off-wake and grid-mean thl')
   call addfld('taui_qt_SP', horiz_only, 'A',  '1/s', &
      'Inverse of damping time scale of the difference between off-wake and grid-mean qt')
   call addfld('taui_u_SP', horiz_only, 'A',   '1/s', &
      'Inverse of damping time scale of the difference between off-wake and grid-mean u')
   call addfld('taui_v_SP', horiz_only, 'A',   '1/s', &
      'Inverse of damping time scale of the difference between off-wake and grid-mean v')
   call addfld('taui_awk_SP', horiz_only, 'A', '1/s', &
      'Inverse of damping time scale of the wake area')

   call addfld('del_org_SP', horiz_only, 'A',  '1/s', &
      'Detrainment rate of the cold pool from UNICON')
   call addfld('del0_org_SP', horiz_only, 'A', '1/s', &
      'Effective detrainment rate of the cold pool from UNICON')

#endif

end subroutine unicon_cam_init

!==================================================================================================

function unicon_implements_cnst(name)

   ! Return true if specified constituent is implemented by this package

   character(len=*), intent(in) :: name  ! constituent name
   logical :: unicon_implements_cnst     ! return value

   integer :: m
   !-----------------------------------------------------------------------

   unicon_implements_cnst = .false.

#ifdef USE_UNICON

   do m = 1, n_org
      if (name == cnst_names(m)) then
         unicon_implements_cnst = .true.
         return
      end if
   end do

#endif

end function unicon_implements_cnst

!==================================================================================================

subroutine unicon_init_cnst(name, q, gcid)

   ! Initialize constituents if they are not read from the initial file

   character(len=*), intent(in)  :: name     ! constituent name
   real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
   integer,          intent(in)  :: gcid(:)  ! global column id
   !-----------------------------------------------------------------------

#ifdef USE_UNICON

   if ( name == 'ORGawk' ) then
      q = 0.0_r8
      return
! JHYoon
!  else if ( name == 'ORGthl' ) then
!     q = 100.0_r8
!     return
!  else if ( name == 'ORGqto' ) then
!     q = 100.0_r8
!     return
!  else if ( name == 'ORGuoo' ) then
!     q = 100.0_r8
!     return
!  else if ( name == 'ORGvoo' ) then
!     q = 100.0_r8
!     return
   else if ( name == 'ORGthl' ) then
      q = 10.0_r8
      return
   else if ( name == 'ORGqto' ) then
      q = 0.01_r8
      return
   else if ( name == 'ORGuoo' ) then
      q = 10.0_r8
      return
   else if ( name == 'ORGvoo' ) then
      q = 10.0_r8
      return
! JHYoon
   end if

#endif

end subroutine unicon_init_cnst

!==================================================================================================

subroutine unicon_cam_tend(dt, state, cam_in, sgh30, &
                           pbuf, ptend, out)


   ! ---------------------- !
   ! Input-output Arguments !
   ! ---------------------- !

   real(r8),                  intent(in)  :: dt         ! Time step [s]
   type(physics_state),       intent(in)  :: state      ! Physics state variables
   type(cam_in_t),            intent(in)  :: cam_in     ! import state
   real(r8),                  intent(in)  :: sgh30(mix) ! Std dev of subgrid topographic height at 30 s horizontal area [ meter ]
   type(physics_buffer_desc), pointer     :: pbuf(:)    ! physics buffer
   type(physics_ptend),       intent(out) :: ptend      ! parameterization tendencies
   type(unicon_out_t),        intent(out) :: out        ! parameterization outputs

   ! -------------------------------------------------------- !
   ! Internal output and local variables for positive tracers !
   ! -------------------------------------------------------- !
   
   real(r8) :: sten_ori(mix,mkx)         !  Tendency of dry static energy [ J / kg / s ]
   real(r8) :: qvten_ori(mix,mkx)        !  Tendency of water vapor specific humidity [ kg / kg / s ]
   real(r8) :: qlten_ori(mix,mkx)        !  Tendency of liquid water mixing ratio [ kg / kg / s ]
   real(r8) :: qiten_ori(mix,mkx)        !  Tendency of ice mixing ratio [ kg / kg / s ]
   real(r8) :: trten_ori(mix,mkx,ncnst)  !  Tendency of tracers [ # / kg / s, kg / kg / s ]

   real(r8) :: slten_pos_inv(mix,mkx)    ! 
   real(r8) :: qtten_pos_inv(mix,mkx)    ! 
   real(r8) :: uten_pos_inv(mix,mkx)     ! 
   real(r8) :: vten_pos_inv(mix,mkx)     ! 
   real(r8) :: sten_pos_inv(mix,mkx)     ! 
   real(r8) :: qvten_pos_inv(mix,mkx)    ! 
   real(r8) :: qlten_pos_inv(mix,mkx)    ! 
   real(r8) :: qiten_pos_inv(mix,mkx)    ! 
   real(r8) :: trten_pos_inv(mix,mkx,ncnst)

   ! --------------- !
   ! Local variables !
   ! --------------- !

   integer :: iend
   integer :: lchnk     
   integer :: itim

   ! fields in physics buffer
   real(r8), pointer, dimension(:,:) :: & ! (mix,mkx)
      ast0_inv    ! Physical stratus fraction [ fraction ]

   real(r8), pointer, dimension(:,:) :: & ! (mix,mkx+1)
      tke0_inv,    &! TKE at the interface [ m2/s2 ]
      bprod0_inv    ! Buoyancy production at the interface [ m2/s3 ]

   integer(i4), pointer, dimension(:) :: & ! (mix)
      kpblh_inv       ! Layer index with PBL top in it or at the base interface

   real(r8), pointer, dimension(:) :: & ! (mix)
      pblh,          &! PBL top height [m]
      went,          &! Entrainment rate at the PBL top interface directly from UW PBL scheme [ m / s ]. went = 0 for STL.
      cush,          &! Cumulus top height [ m ]
      cushavg,       &! Mean cumulus top height weighted by updraft mass flux at surface [ m ]
      cuorg,         &! Convective organization parameter [ 0-1 ]

      awk_PBL,       &! Wake area within PBL [ 0 - 1 ]
      delta_thl_PBL, &! Diff of thl between off-wake region and grid-mean value averaged over the PBL [ K ]
      delta_qt_PBL,  &! Diff of qt  between off-wake region and grid-mean value averaged over the PBL [ kg/kg ]
      delta_u_PBL,   &! Diff of u   between off-wake region and grid-mean value averaged over the PBL [ m/s ]
      delta_v_PBL     ! Diff of v   between off-wake region and grid-mean value averaged over the PBL [ m/s ]

   real(r8), pointer, dimension(:,:) :: & ! (mix,ncnst)
      delta_tr_PBL ! Diff of tr  between off-wake region and grid-mean value avg over the PBL [ kg/kg, #/kg ]

   real(r8), dimension(mix,mkx) :: & ! (mix,mkx)
      cu_cmfum  ! The mass involved in the updraft buoyancy sorting at the previous time step [ kg/s/m2 ]

   real(r8), pointer, dimension(:,:) :: & ! (mix,mkx)
      cu_cmfr, &! The detrained mass from convective up and downdraft at the previous time step [ kg/s/m2 ]
      cu_thlr, &! Mass-flux wghted mean 'thl' of detrained mass from conv up and downdraft at prev step [ K ]
      cu_qtr,  &! Mass-flux wghted mean 'qt'  of detrained mass from conv up and downdraft at prev step [ kg/kg ]
      cu_ur,   &! Mass-flux wghted mean 'u'   of detrained mass from conv up and downdraft at prev step [ m/s ]
      cu_vr,   &! Mass-flux wghted mean 'v'   of detrained mass from conv up and downdraft at prev step [ m/s ]
      cu_qlr,  &! Mass-flux wghted mean 'ql'  of detrained mass from conv up and downdraft at prev step [ kg/kg ]
      cu_qir,  &! Mass-flux wghted mean 'qi'  of detrained mass from conv up and downdraft at prev step [ kg/kg ]
      cmfr_det,&! The detrained mass from convective up and downdraft at the previous time step [ kg/s/m2 ]
      qlr_det, &! Mass-flux wghted mean 'ql'  of detrained mass from conv up and downdraft at prev step [ kg/kg ]
      qir_det, &! Mass-flux wghted mean 'qi'  of detrained mass from conv up and downdraft at prev step [ kg/kg ]
      rqcr_l,  &! Grid-mean production rate of the mass   of detrained convective liquid condensate [ #/kg/s ] >= 0.
      rqcr_i,  &! Grid-mean production rate of the mass   of detrained convective ice    condensate [ #/kg/s ] >= 0.
      rncr_l,  &! Grid-mean production rate of the number of detrained convective liquid condensate [ #/kg/s ] >= 0.
      rncr_i    ! Grid-mean production rate of the number of detrained convective ice    condensate [ #/kg/s ] >= 0.

   real(r8), pointer, dimension(:,:,:) :: & ! (mix,mkx,ncnst)
      cu_trr    ! Mass-flux wghted mean 'tr'  of detrained mass from conv up and downdraft at prev step [ kg/kg ]

   real(r8), dimension(mix,mkx) :: & ! (mix,mkx)
      cu_cmfrd,&! The amount of detrained mass from convective downdraft at the previous time step [ kg/s/m2 ]
      cu_thlrd,&! Mass-flux wghted mean 'thl' of detrained mass from conv downdraft at prev step [ K ]
      cu_qtrd, &! Mass-flux wghted mean 'qt'  of detrained mass from conv downdraft at prev step [ kg/kg ]
      cu_urd,  &! Mass-flux wghted mean 'u'   of detrained mass from conv downdraft at prev step [ m/s ]
      cu_vrd,  &! Mass-flux wghted mean 'v'   of detrained mass from conv downdraft at prev step [ m/s ]
      cu_qlrd, &! Mass-flux wghted mean 'ql'  of detrained mass from conv downdraft at prev step [ kg/kg ]
      cu_qird   ! Mass-flux wghted mean 'qi'  of detrained mass from conv downdraft at prev step [ kg/kg ]

   real(r8), dimension(mix,mkx,ncnst) :: & ! (mix,mkx,ncnst)
      cu_trrd   ! Mass-flux wghted mean 'tr'  of detrained mass from conv downdraft at prev step [ kg/kg ]

   real(r8), pointer, dimension(:,:) :: & ! (mix,mkx)
      shfrc,     &! Convective updraft fractional area
      icwmrsh,   &! In-cloud LWC+IWC within convective updraft
      rprdsh,    &! Prod of rain+snow by lateral expels of cumulus condensate [ kg / kg / s ]
      evapc_inv, &! Evaporation rate of convective precipitation within environment [ kg/kg/s ]
      am_evp_st_inv,   &!  Evaporation area of stratiform precipitation [fraction]
      evprain_st_inv,  &!  Grid-mean evaporation rate of stratiform rain [kg/kg/s] >= 0.
      evpsnow_st_inv    !  Grid-mean evaporation rate of stratiform snow [kg/kg/s] >= 0.

   real(r8), pointer, dimension(:) :: & ! (mix)
      precip,    &!  Precipitation flux at surface in flux unit [ m / s ]
      snow,      &!  Snow flux at surface in flux unit [ m / s ]
      rice2       ! Vertical integral of reserved detrained convective ice condensate converted in liquid unit [m/s] >= 0.

   real(r8) :: ps0(mix,0:mkx)            !  Environmental pressure at full sigma levels
   real(r8) :: zs0(mix,0:mkx)            !  Environmental height at full sigma levels
   real(r8) :: p0(mix,mkx)               !  Environmental pressure at half sigma levels
   real(r8) :: z0(mix,mkx)               !  Environmental height at half sigma levels
   real(r8) :: dp0(mix,mkx)              !  Environmental layer pressure thickness
   real(r8) :: dpdry0(mix,mkx)           !  Environmental layer dry pressure thickness
   real(r8) :: u0(mix,mkx)               !  Environmental zonal wind
   real(r8) :: v0(mix,mkx)               !  Environmental meridional wind
   real(r8) :: qv0(mix,mkx)              !  Environmental specific humidity
   real(r8) :: ql0(mix,mkx)              !  Environmental liquid water mixing ratio
   real(r8) :: qi0(mix,mkx)              !  Environmental ice mixing ratio
   real(r8) :: tr0(mix,mkx,ncnst)        !  Environmental tracers [ #/kg, kg/kg ]
   real(r8) :: t0(mix,mkx)               !  Environmental temperature
   real(r8) :: s0(mix,mkx)               !  Environmental dry static energy
   real(r8) :: ast0(mix,mkx)             !  Physical stratiform cloud fraction [ fraction ]
   real(r8) :: tke0(mix,0:mkx)           !  TKE [ m2/s2 ]
   real(r8) :: bprod0(mix,0:mkx)         !  Buoyancy production [ m2/s3 ]
   real(r8) :: am_evp_st(mix,mkx)        !  Evaporation area of stratiform precipitation [fraction]
   real(r8) :: evprain_st(mix,mkx)       !  Grid-mean evaporation rate of stratiform rain [kg/kg/s] >= 0.
   real(r8) :: evpsnow_st(mix,mkx)       !  Grid-mean evaporation rate of stratiform snow [kg/kg/s] >= 0.

   integer(i4) :: kpblh(mix)             !  Layer index with PBL top in it or at the base interface


   real(r8) :: am_u(mix,mkx)             !  Convective updraft fractional area
   real(r8) :: qlm_u(mix,mkx)            !  In-cloud LWC within convective updraft [ kg / kg ]
   real(r8) :: qim_u(mix,mkx)            !  In-cloud IWC within convective updraft [ kg / kg ]

   real(r8) :: am_d(mix,mkx)             !  Convective downdraft fractional area
   real(r8) :: qlm_d(mix,mkx)            !  In-cloud LWC within downdraft updraft [ kg / kg ]
   real(r8) :: qim_d(mix,mkx)            !  In-cloud IWC within downdraft updraft [ kg / kg ]

   real(r8) :: cmf_u(mix,0:mkx)          !  Upward     convective mass flux at the interface [ kg / s / m2 ]
   real(r8) :: slflx(mix,0:mkx)          !  Net upward convective flux of liquid static energy [ J / s / m2 ]
   real(r8) :: qtflx(mix,0:mkx)          !  Net upward convective flux of total specific humidity [ kg / s / m2 ]

   real(r8) :: qvten(mix,mkx)            !  Tendency of water vapor specific humidity [ kg / kg / s ]
   real(r8) :: qlten(mix,mkx)            !  Tendency of liquid water mixing ratio [ kg / kg / s ]
   real(r8) :: qiten(mix,mkx)            !  Tendency of ice mixing ratio [ kg / kg / s ]
   real(r8) :: trten(mix,mkx,ncnst)      !  Tendency of tracers [ # / kg / s, kg / kg / s ]

   real(r8) :: sten(mix,mkx)             !  Tendency of dry static energy [ J / kg / s ]
   real(r8) :: uten(mix,mkx)             !  Tendency of zonal wind [ m / s / s ]
   real(r8) :: vten(mix,mkx)             !  Tendency of meridional wind [ m / s / s ]

   real(r8) :: qrten(mix,mkx)            !  Production rate of rain by lateral expels of cumulus condensate [ kg / kg / s ]
   real(r8) :: qsten(mix,mkx)            !  Production rate of snow by lateral expels of cumulus condensate [ kg / kg / s ]

   real(r8) :: evapc(mix,mkx)            !  Evaporation rate of convective precipitation within environment [ kg/kg/s ]

   real(r8) :: rqc(mix,mkx)              !  Production rate of detrained LWC+IWC [ kg / kg / s ] > 0
   real(r8) :: rqc_l(mix,mkx)            !  Production rate of detrained LWC     [ kg / kg / s ] > 0
   real(r8) :: rqc_i(mix,mkx)            !  Production rate of detrained IWC     [ kg / kg / s ] > 0
   real(r8) :: rnc_l(mix,mkx)            !  Production rate of detrained droplet number of cloud liquid droplets [ # / kg / s ] > 0
   real(r8) :: rnc_i(mix,mkx)            !  Production rate of detrained droplet number of cloud    ice droplets [ # / kg / s ] > 0

   real(r8) :: cnt(mix)                  !  Cloud top  interface index ( ki = kpen )
   real(r8) :: cnb(mix)                  !  Cloud base interface index ( ki = krel-1 )

   real(r8) :: pdel0(mix,mkx)            !  Environmental pressure thickness ( either dry or moist ) [ Pa ]
   real(r8) :: trmin                     !  Minimum concentration of tracer [ # / kg ]

   real(r8) :: cmf_det(mix,mkx)          ! Det mass flux from convective updraft (not from environment) and downdraft [kg/s/m^2]
   real(r8) :: ql_det(mix,mkx)           ! Det conv LWC from convective updraft (not from environment) and downdraft [kg/s/m^2]
   real(r8) :: qi_det(mix,mkx)           ! Det conv IWC from convective updraft (not from environment) and downdraft [kg/s/m^2]

   ! For prognostically updated state variables

   real(r8) :: qv0_c(mix,mkx)            !  Environmental specific humidity
   real(r8) :: ql0_c(mix,mkx)            !  Environmental liquid water mixing ratio
   real(r8) :: qi0_c(mix,mkx)            !  Environmental ice mixing ratio
   real(r8) :: t0_c(mix,mkx)             !  Environmental temperature
   real(r8) :: s0_c(mix,mkx)             !  Environmental dry static energy
   real(r8) :: tr0_c(mix,mkx,ncnst)      !  Environmental tracers [ # / kg, kg / kg ]

   ! Layer index variables
   integer  :: k                         !  Vertical index for local fields 
   integer  :: k_inv                     !  Vertical index for incoming fields
   integer  :: mt                        !  Tracer index [ no ]
   integer  :: m

   ! For aerosol tendency output
   character(len=30) :: varname

   logical :: lq(ncnst)

   ! --------- !
   ! Main body !
   ! --------- !

#ifdef USE_UNICON

   iend  = state%ncol
   lchnk = state%lchnk
   
   ! Associate pointers with physics buffer fields

   itim = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx,         ast0_inv,    start=(/1,1,itim/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, tke_idx,         tke0_inv)
   call pbuf_get_field(pbuf, bprod_idx,       bprod0_inv)
   call pbuf_get_field(pbuf, kpblh_idx,       kpblh_inv)
   call pbuf_get_field(pbuf, pblh_idx,        pblh)
   call pbuf_get_field(pbuf, went_idx,        went)
   call pbuf_get_field(pbuf, cush_idx,        cush,       start=(/1,itim/), kount=(/pcols,1/))
   call pbuf_get_field(pbuf, cushavg_idx,     cushavg,    start=(/1,itim/), kount=(/pcols,1/))
   call pbuf_get_field(pbuf, cuorg_idx,       cuorg,      start=(/1,itim/), kount=(/pcols,1/))

   call pbuf_get_field(pbuf, awk_PBL_idx,     awk_PBL,    start=(/1,itim/), kount=(/pcols,1/))
   call pbuf_get_field(pbuf, delta_thl_PBL_idx, delta_thl_PBL, start=(/1,itim/), kount=(/pcols,1/))
   call pbuf_get_field(pbuf, delta_qt_PBL_idx,  delta_qt_PBL,  start=(/1,itim/), kount=(/pcols,1/))
   call pbuf_get_field(pbuf, delta_u_PBL_idx,   delta_u_PBL,   start=(/1,itim/), kount=(/pcols,1/))
   call pbuf_get_field(pbuf, delta_v_PBL_idx,   delta_v_PBL,   start=(/1,itim/), kount=(/pcols,1/))
   call pbuf_get_field(pbuf, delta_tr_PBL_idx,  delta_tr_PBL,  start=(/1,1,itim/), kount=(/pcols,pcnst,1/))

   call pbuf_get_field(pbuf, cu_cmfr_idx,  cu_cmfr,  start=(/1,1,itim/),   kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, cu_thlr_idx,  cu_thlr,  start=(/1,1,itim/),   kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, cu_qtr_idx,   cu_qtr,   start=(/1,1,itim/),   kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, cu_ur_idx,    cu_ur,    start=(/1,1,itim/),   kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, cu_vr_idx,    cu_vr,    start=(/1,1,itim/),   kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, cu_qlr_idx,   cu_qlr,   start=(/1,1,itim/),   kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, cu_qir_idx,   cu_qir,   start=(/1,1,itim/),   kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, cu_trr_idx,   cu_trr,   start=(/1,1,1,itim/), kount=(/pcols,pver,pcnst,1/))

   call pbuf_get_field(pbuf, shfrc_idx,        shfrc)
   call pbuf_get_field(pbuf, icwmrsh_idx,     icwmrsh)
   call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh)
   call pbuf_get_field(pbuf, prec_sh_idx,     precip)
   call pbuf_get_field(pbuf, snow_sh_idx,     snow)
   call pbuf_get_field(pbuf, nevapr_shcu_idx, evapc_inv)
   call pbuf_get_field(pbuf, am_evp_st_idx,   am_evp_st_inv)
   call pbuf_get_field(pbuf, evprain_st_idx,  evprain_st_inv)
   call pbuf_get_field(pbuf, evpsnow_st_idx,  evpsnow_st_inv)

   call pbuf_get_field(pbuf, cmfr_det_idx,    cmfr_det)
   call pbuf_get_field(pbuf, qlr_det_idx,      qlr_det)
   call pbuf_get_field(pbuf, qir_det_idx,      qir_det)

   call pbuf_get_field(pbuf, rqcr_l_idx,        rqcr_l)
   call pbuf_get_field(pbuf, rqcr_i_idx,        rqcr_i)
   call pbuf_get_field(pbuf, rncr_l_idx,        rncr_l)
   call pbuf_get_field(pbuf, rncr_i_idx,        rncr_i)

   call pbuf_get_field(pbuf, rice2_idx,          rice2)

   ! Reverse variables defined at the layer mid-point

   do k = 1, mkx
      k_inv               = mkx + 1 - k
      p0(:iend,k)         = state%pmid(:iend,k_inv)
      u0(:iend,k)         = state%u(:iend,k_inv)
      v0(:iend,k)         = state%v(:iend,k_inv)
      z0(:iend,k)         = state%zm(:iend,k_inv)
      dp0(:iend,k)        = state%pdel(:iend,k_inv)
      dpdry0(:iend,k)     = state%pdeldry(:iend,k_inv)
      qv0(:iend,k)        = state%q(:iend,k_inv,1)
      ql0(:iend,k)        = state%q(:iend,k_inv,ixcldliq)
      qi0(:iend,k)        = state%q(:iend,k_inv,ixcldice)
      t0(:iend,k)         = state%t(:iend,k_inv)
      s0(:iend,k)         = state%s(:iend,k_inv)
      ast0(:iend,k)       = ast0_inv(:iend,k_inv)
      am_evp_st(:iend,k)  = am_evp_st_inv(:iend,k_inv)
      evprain_st(:iend,k) = evprain_st_inv(:iend,k_inv)
      evpsnow_st(:iend,k) = evpsnow_st_inv(:iend,k_inv)
      do mt = 1, ncnst
         tr0(:iend,k,mt)  = state%q(:iend,k_inv,mt)
      end do
   end do

   ! Reverse variables defined at the interfaces
    
   do k = 0, mkx
      k_inv               = mkx + 1 - k
      ps0(:iend,k)        = state%pint(:iend,k_inv)
      zs0(:iend,k)        = state%zi(:iend,k_inv)
      tke0(:iend,k)       = tke0_inv(:iend,k_inv)
      bprod0(:iend,k)     = bprod0_inv(:iend,k_inv)
   end do

   ! Reverse the index of ambiguous layer

   kpblh(:iend) = mkx + 1 - kpblh_inv(:iend)

   
   call compute_unicon( mix       , mkx        , iend      , ncnst    , dt      ,           &
                        ps0       , zs0        , p0        , z0       , dp0     , dpdry0  , &
                        t0        , qv0        , ql0       , qi0      , tr0     ,           & 
                        u0        , v0         , ast0      , tke0     , bprod0  ,           &
                        kpblh     , pblh       , went      ,                                & 
                        cam_in%cflx(:,1), cam_in%shf, cam_in%wsx, cam_in%wsy, cam_in%cflx , &
                        cam_in%landfrac, sgh30   ,                                          &
                        am_evp_st , evprain_st , evpsnow_st,                                &
                        cush      , cushavg    , cuorg     ,                                &
                        awk_PBL                , delta_thl_PBL        , delta_qt_PBL      , & 
                        delta_u_PBL            , delta_v_PBL          , delta_tr_PBL      , &
                        cu_cmfum  , cu_cmfr    , cu_thlr   , cu_qtr   , cu_ur   , cu_vr   , &
                        cu_qlr    , cu_qir     , cu_trr    ,                                & 
                        cu_cmfrd  , cu_thlrd   , cu_qtrd   , cu_urd   , cu_vrd  ,           &
                        cu_qlrd   , cu_qird    , cu_trrd   ,                                & 
                        am_u      , qlm_u      , qim_u     ,                                &
                        am_d      , qlm_d      , qim_d     ,                                &
                        cmf_u     , slflx      , qtflx     ,                                & 
                        qvten     , qlten      , qiten     , trten    ,                     &
                        sten      , uten       , vten      ,                                &
                        qrten     , qsten      ,                                            &
                        rqc_l     , rqc_i      , rqc       , rnc_l    , rnc_i   ,           &
                        out%rliq  , precip     , snow      , evapc    ,                     &
                        cnt       , cnb        , cmf_det   , ql_det   , qi_det  ,           &
                        lchnk )

      
   ! Initialize output ptend
   lq(:) = .true.
   call physics_ptend_init(ptend, state%psetcols, 'unicon', ls=.true., lu=.true., lv=.true., lq=lq)

   ! ----------------------------------------------------- !
   ! Treatment of reserved liquid-ice water for CAM        !
   ! All the reserved condensate is converted into liquid. !
   ! Jan.04.2012. Relocated from below to here to prevent  !
   !              energy and water conservation error.     !
   !              Also add corresponding tendency of cloud !
   !              droplet number concentration.            !
   ! Note that since 'positive_moisture' will convert      !
   ! negative 'ql,qi,nl,ni' into positive, below may cause !
   ! larger 'ql,qi' and smaller 'qv' than it should be if  !
   ! 'ql0,qi0' further below becomes negative.             !
   ! This feature is inevitable at the current stage.      !
   ! Note that in future, I can impose consistency between !
   ! positive_moisture and positive_tracer for nl,ni.      !     
   ! ----------------------------------------------------- !  

   qlten(:iend,:mkx) = qlten(:iend,:mkx) - rqc_l(:iend,:mkx)
   qiten(:iend,:mkx) = qiten(:iend,:mkx) - rqc_i(:iend,:mkx)
   sten(:iend,:mkx)  =  sten(:iend,:mkx) - ( xls - xlv ) * rqc_i(:iend,:mkx)
   trten(:iend,:mkx,ixnumliq) = trten(:iend,:mkx,ixnumliq) - rnc_l(:iend,:mkx)
   trten(:iend,:mkx,ixnumice) = trten(:iend,:mkx,ixnumice) - rnc_i(:iend,:mkx)

   ! --------------------------------------------- !
   ! Prevent negative cloud condensate and tracers !
   ! --------------------------------------------- !

   qv0_c(:iend,:mkx)  = qv0(:iend,:mkx) + qvten(:iend,:mkx)*dt
   ql0_c(:iend,:mkx)  = ql0(:iend,:mkx) + qlten(:iend,:mkx)*dt
   qi0_c(:iend,:mkx)  = qi0(:iend,:mkx) + qiten(:iend,:mkx)*dt
   t0_c(:iend,:mkx)   =  t0(:iend,:mkx) + (1._r8/cp)*sten(:iend,:mkx)*dt
   s0_c(:iend,:mkx)   =  s0(:iend,:mkx) +  sten(:iend,:mkx)*dt
   do mt = 1, ncnst
      tr0_c(:iend,:mkx,mt)  = tr0(:iend,:mkx,mt) + trten(:iend,:mkx,mt)*dt
   enddo

   ! Note : Since 'positive_moisture' will only perform condensation not the evaporation, 
   !        we don't need to modify corresponding 'nl,ni'.
   !        Thus, current version is completely OK.

   sten_ori(:iend,:mkx)  =  sten(:iend,:mkx)
   qvten_ori(:iend,:mkx) = qvten(:iend,:mkx)
   qlten_ori(:iend,:mkx) = qlten(:iend,:mkx)
   qiten_ori(:iend,:mkx) = qiten(:iend,:mkx)
   do mt = 1, ncnst
      trten_ori(:iend,:mkx,mt) = trten(:iend,:mkx,mt)
   enddo

   call positive_moisture( &
      cp, xlv, xls, mix, iend, &
      mkx, dt, qmin(1), qmin(ixcldliq), qmin(ixcldice), &
      dp0, qv0_c, ql0_c, qi0_c, t0_c, &
      s0_c, qvten, qlten, qiten, sten )

   do mt = 1, ncnst

      if( cnst_get_type_byind(mt) .eq. 'wet' ) then
         pdel0(:iend,:mkx) = dp0(:iend,:mkx)
      else
         pdel0(:iend,:mkx) = dpdry0(:iend,:mkx)
      endif

      if (cnst_is_mam_num(mt)) then
         trmin = 1.e-5_r8
      else
         trmin = qmin(mt)
      end if

      call positive_tracer( mix, iend, mkx, dt, trmin, pdel0, tr0_c(:,:,mt), trten(:,:,mt) )

   enddo

   ! ----------------------------------------------------- !
   ! Treatment of reserved liquid-ice water for CAM        !
   ! All the reserved condensate is converted into liquid. !
   ! ----------------------------------------------------- !

   ! Important : While below treatment looks different from what is being done in the uwshcu.F90,
   !             below is actually the same as what is being done in the uwshcu.F90. The process
   !             used below is (1) convert detrained ice into liquid, which involves melting cooling,
   !             and then (2) subtract this total detrained liquid + ice from the grid-mean qlten.
   ! Sep.08.2010. In the new scheme, it will be 'rqc_l = rqc_i = 0'. Thus, below block does nothing.
   !              But it is no harm to keep below block.
   ! Jan.04.2012. Below block will be used again for many reasons : (1) inctease TGCLDLWP, (2) reduce PREH20,
   !              (2) reduce too strong SWCF over the far eastern equatorial Pacific. 
   !              Currently, below produces conservation error of both energy and water.
   ! Jan.04.2012. Conservation errors of energy and moisture dissappears if I locate below block
   !              before applying 'positive_moisture'. Thus, I relocated below block above.

   ! qlten(:iend,:mkx) = qlten(:iend,:mkx) + rqc_i(:iend,:mkx)
   ! qiten(:iend,:mkx) = qiten(:iend,:mkx) - rqc_i(:iend,:mkx)
   ! sten(:iend,:mkx)  =  sten(:iend,:mkx) - ( xls - xlv ) * rqc_i(:iend,:mkx)
   ! qlten(:iend,:mkx) = qlten(:iend,:mkx) - rqc_l(:iend,:mkx) - rqc_i(:iend,:mkx)
   ! trten(:iend,:mkx,ixnumliq) = trten(:iend,:mkx,ixnumliq) - rnc_l(:iend,:mkx)
   ! trten(:iend,:mkx,ixnumice) = trten(:iend,:mkx,ixnumice) - rnc_i(:iend,:mkx)

   ! --------------------------- !
   ! Reverse vertical coordinate !
   ! --------------------------- !

   ! Reverse cloud top/base interface indices

   out%cnt(:iend) = real(mkx,r8) + 1._r8 - cnt(:iend)
   out%cnb(:iend) = real(mkx,r8) + 1._r8 - cnb(:iend)

   ! Reverse variables defined at the interfaces
   
   do k = 0, mkx
      k_inv                     = mkx + 1 - k
      out%cmfmc(:iend,k_inv)    = cmf_u(:iend,k)       !  Updraft mass flux at top of layer
      out%slflx(:iend,k_inv)    = slflx(:iend,k)       !  Net liquid static energy flux
      out%qtflx(:iend,k_inv)    = qtflx(:iend,k)       !  Net total water flux
   end do

   ! Reverse variables defined at the layer mid-point

   do k = 1, mkx
      k_inv = mkx + 1 - k
      ptend%q(:iend,k_inv,1)        = qvten(:iend,k)       ! Convective tendency of specific humidity
      ptend%q(:iend,k_inv,ixcldliq) = qlten(:iend,k)       ! Convective tendency of liquid water mixing ratio
      ptend%q(:iend,k_inv,ixcldice) = qiten(:iend,k)       ! Convective tendency of ice mixing ratio
      ptend%s(:iend,k_inv)      = sten(:iend,k)        ! Convective tendency of static energy
      ptend%u(:iend,k_inv)      = uten(:iend,k)        ! Convective tendency of zonal wind
      ptend%v(:iend,k_inv)      = vten(:iend,k)        ! Convective tendency of meridional wind

      out%rqc(:iend,k_inv)      = rqc(:iend,k)         ! Convective tendency of reserved (suspended) liquid + ice water

      ! These quantities are being put into physics buffer fields that were meant for shallow convection.
      ! This should be fixed.
      shfrc(:iend,k_inv)   = am_u(:iend,k)
      icwmrsh(:iend,k_inv) = qlm_u(:iend,k) + qim_u(:iend,k)
      rprdsh(:iend,k_inv)  = qrten(:iend,k) + qsten(:iend,k)
      evapc_inv(:iend,k_inv) = evapc(:iend,k)       ! Evaporation rate of convective precipitation within environment

      do mt = 2, ncnst
         if (mt /= ixcldliq .and. mt /= ixcldice) then
            ptend%q(:iend,k_inv,mt) = trten(:iend,k,mt)
         end if
      enddo

      ! Additional diagnostic output associated with positive tracer

      sten_pos_inv(:iend,k_inv)   = sten(:iend,k)  -  sten_ori(:iend,k)
      qvten_pos_inv(:iend,k_inv)  = qvten(:iend,k) - qvten_ori(:iend,k)
      qlten_pos_inv(:iend,k_inv)  = qlten(:iend,k) - qlten_ori(:iend,k)
      qiten_pos_inv(:iend,k_inv)  = qiten(:iend,k) - qiten_ori(:iend,k)
      qtten_pos_inv(:iend,k_inv)  = qvten_pos_inv(:iend,k_inv) + qlten_pos_inv(:iend,k_inv) &
                                                               + qiten_pos_inv(:iend,k_inv)
      slten_pos_inv(:iend,k_inv)  = sten_pos_inv(:iend,k_inv) - xlv * qlten_pos_inv(:iend,k_inv) &
                                                              - xls * qiten_pos_inv(:iend,k_inv)
      uten_pos_inv(:iend,k_inv)   = 0._r8
      vten_pos_inv(:iend,k_inv)   = 0._r8
      do mt = 1, ncnst
         trten_pos_inv(:iend,k_inv,mt)  = trten(:iend,k,mt) - trten_ori(:iend,k,mt)
      enddo

      ! Save detrainment variables to pbuf for use in the cloud macrophysics
      cmfr_det(:iend,k_inv) =  cmf_det(:iend,k)
      qlr_det(:iend,k_inv)  =   ql_det(:iend,k)
      qir_det(:iend,k_inv)  =   qi_det(:iend,k)

      rqcr_l(:iend,k_inv)   =    rqc_l(:iend,k)
      rqcr_i(:iend,k_inv)   =    rqc_i(:iend,k)
      rncr_l(:iend,k_inv)   =    rnc_l(:iend,k)
      rncr_i(:iend,k_inv)   =    rnc_i(:iend,k)

   end do

   call outfld('slten_pos_SP' , slten_pos_inv, mix, lchnk) 
   call outfld('qtten_pos_SP' , qtten_pos_inv, mix, lchnk) 
   call outfld('uten_pos_SP'  ,  uten_pos_inv, mix, lchnk) 
   call outfld('vten_pos_SP'  ,  vten_pos_inv, mix, lchnk)
   call outfld('sten_pos_SP'  ,  sten_pos_inv, mix, lchnk) 
   call outfld('qvten_pos_SP' , qvten_pos_inv, mix, lchnk) 
   call outfld('qlten_pos_SP' , qlten_pos_inv, mix, lchnk) 
   call outfld('qiten_pos_SP' , qiten_pos_inv, mix, lchnk) 
   call outfld('nlten_pos_SP' , trten_pos_inv(:,:,ixnumliq), mix, lchnk) 
   call outfld('niten_pos_SP' , trten_pos_inv(:,:,ixnumice), mix, lchnk) 
   call outfld('CMFR_DET'     , cmfr_det,      pcols, lchnk)
   call outfld('QLR_DET'      , qlr_det,       pcols, lchnk)
   call outfld('QIR_DET'      , qir_det,       pcols, lchnk)

   do m = 1, ncnst
      if (cnst_is_mam_num(m) .or. cnst_is_mam_mmr(m)) then
         varname = trim(cnst_name(m))//'_pos_SP'
         call outfld(trim(varname), trten_pos_inv(:,:,m), mix, lchnk)
      end if
   end do

#endif

end subroutine unicon_cam_tend

subroutine unicon_cam_org_diags(state, pbuf)

   ! ------------------------------------------------------------------------
   ! Insert the organization-related heterogeneities computed inside the
   ! UNICON into the tracer arrays here before performing advection.
   ! This is necessary to prevent any modifications of organization-related
   ! heterogeneities by non convection-advection process, such as
   ! dry and wet deposition of aerosols, MAM, etc.
   ! Again, note that only UNICON and advection schemes are allowed to
   ! changes to organization at this stage, although we can include the
   ! effects of other physical processes in future.
   ! ------------------------------------------------------------------------

   ! Arguments
   type(physics_state), intent(inout) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)

   ! Local variables
   real(r8), pointer, dimension(:  ) :: awk_PBL
   real(r8), pointer, dimension(:  ) :: delta_thl_PBL
   real(r8), pointer, dimension(:  ) :: delta_qt_PBL
   real(r8), pointer, dimension(:  ) :: delta_u_PBL
   real(r8), pointer, dimension(:  ) :: delta_v_PBL

   integer :: i, itim, ncol
   ! ------------------------------------------------------------------------

#ifdef USE_UNICON

   ncol = state%ncol

   itim = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, awk_PBL_idx,       awk_PBL,       start=(/1,itim/), kount=(/pcols,1/))
   call pbuf_get_field(pbuf, delta_thl_PBL_idx, delta_thl_PBL, start=(/1,itim/), kount=(/pcols,1/))
   call pbuf_get_field(pbuf, delta_qt_PBL_idx,  delta_qt_PBL,  start=(/1,itim/), kount=(/pcols,1/))
   call pbuf_get_field(pbuf, delta_u_PBL_idx,   delta_u_PBL,   start=(/1,itim/), kount=(/pcols,1/))
   call pbuf_get_field(pbuf, delta_v_PBL_idx,   delta_v_PBL,   start=(/1,itim/), kount=(/pcols,1/))

   do i = 1, ncol

      state%q(i,:,awk_cnst_ind) = awk_PBL(i) 

      ! Add a constant offset of '100._r8' to 'delta_xxx' variables to be
      ! consistent with the reading sentence of unicon.F90 and so to prevent
      ! negative tracers.
! JHYoon
!     state%q(i,:,thl_cnst_ind) = max( 0._r8, delta_thl_PBL(i) + 100._r8 )
!     state%q(i,:,qt_cnst_ind)  = max( 0._r8,  delta_qt_PBL(i) + 100._r8 )
!     state%q(i,:,u_cnst_ind)   = max( 0._r8,   delta_u_PBL(i) + 100._r8 )
!     state%q(i,:,v_cnst_ind)   = max( 0._r8,   delta_v_PBL(i) + 100._r8 )
      state%q(i,:,thl_cnst_ind) = max( 0._r8, delta_thl_PBL(i) + 10._r8 )
      state%q(i,:,qt_cnst_ind)  = max( 0._r8,  delta_qt_PBL(i) + 0.01_r8 )
      state%q(i,:,u_cnst_ind)   = max( 0._r8,   delta_u_PBL(i) + 10._r8 )
      state%q(i,:,v_cnst_ind)   = max( 0._r8,   delta_v_PBL(i) + 10._r8 )
! JHYoon
   end do

#endif

end subroutine unicon_cam_org_diags

!==================================================================================================
end module unicon_cam
