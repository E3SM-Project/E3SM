module cloud_rad_props

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use ppgrid,           only: pcols, pver, pverp
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field, pbuf_old_tim_idx
use radconstants,     only: nswbands, nlwbands, idx_sw_diag, ot_length, idx_lw_diag
use cam_abortutils,       only: endrun
use rad_constituents, only: iceopticsfile, liqopticsfile
use oldcloud,         only: oldcloud_lw, old_liq_get_rad_props_lw, old_ice_get_rad_props_lw, oldcloud_init

use ebert_curry,      only: scalefactor
use cam_logfile,      only: iulog

use interpolate_data, only: interp_type, lininterp_init, lininterp, &
     extrap_method_bndry, lininterp_finish

implicit none
private
save

public :: &
   cloud_rad_props_init,          &
   cloud_rad_props_get_sw,        & ! return SW optical props of total bulk aerosols
   cloud_rad_props_get_lw,        & ! return LW optical props of total bulk aerosols
   get_ice_optics_sw,             & ! return Mitchell SW ice radiative properties
   ice_cloud_get_rad_props_lw,    & ! Mitchell LW ice rad props
   mc6_ice_get_rad_props_lw,      & ! MC6 ice optics scheme, added by UM team on Dec.18
   get_liquid_optics_sw,          & ! return Conley SW rad props
   liquid_cloud_get_rad_props_lw, & ! return Conley LW rad props
   snow_cloud_get_rad_props_lw,   &
   get_snow_optics_sw

integer :: nmu, nlambda
real(r8), allocatable :: g_mu(:)           ! mu samples on grid
real(r8), allocatable :: g_lambda(:,:)     ! lambda scale samples on grid
real(r8), allocatable :: ext_sw_liq(:,:,:)
real(r8), allocatable :: ssa_sw_liq(:,:,:)
real(r8), allocatable :: asm_sw_liq(:,:,:)
real(r8), allocatable :: abs_lw_liq(:,:,:)

integer :: n_g_d
real(r8), allocatable :: g_d_eff(:)        ! radiative effective diameter samples on grid
real(r8), allocatable :: ext_sw_ice(:,:)
real(r8), allocatable :: ssa_sw_ice(:,:)
real(r8), allocatable :: asm_sw_ice(:,:)
real(r8), allocatable :: abs_lw_ice(:,:)

! 
! indexes into pbuf for optical parameters of MG clouds
! 
   integer :: i_dei, i_mu, i_lambda, i_iciwp, i_iclwp, i_des, i_icswp
   integer :: i_rei, i_cld ! added by UM team on Dec.18

! indexes into constituents for old optics
   integer :: &
        ixcldice,           & ! cloud ice water index
        ixcldliq              ! cloud liquid water index


!==============================================================================
contains
!==============================================================================

subroutine cloud_rad_props_init()

   use netcdf
   use spmd_utils,     only: masterproc
   use ioFileMod,      only: getfil
   use error_messages, only: handle_ncerr
#if ( defined SPMD )
   use mpishorthand
#endif
   use constituents,   only: cnst_get_ind
   use slingo,         only: slingo_rad_props_init
   use ebert_curry,    only: ec_rad_props_init, scalefactor

   character(len=256) :: liquidfile 
   character(len=256) :: icefile 
   character(len=256) :: locfn

   integer :: ncid, dimid, f_nlwbands, f_nswbands, ierr
   integer :: vdimids(NF90_MAX_VAR_DIMS), ndims, templen
   ! liquid clouds
   integer :: mudimid, lambdadimid
   integer :: mu_id, lambda_id, ext_sw_liq_id, ssa_sw_liq_id, asm_sw_liq_id, abs_lw_liq_id

   ! ice clouds
   integer :: d_dimid ! diameters
   integer :: d_id, ext_sw_ice_id, ssa_sw_ice_id, asm_sw_ice_id, abs_lw_ice_id

   integer :: err

   liquidfile = liqopticsfile 
   icefile = iceopticsfile

   call slingo_rad_props_init
   call ec_rad_props_init
   call oldcloud_init

   i_dei    = pbuf_get_index('DEI',errcode=err)
   i_mu     = pbuf_get_index('MU',errcode=err)
   i_lambda = pbuf_get_index('LAMBDAC',errcode=err)
   i_iciwp  = pbuf_get_index('ICIWP',errcode=err)
   i_iclwp  = pbuf_get_index('ICLWP',errcode=err)
   i_des    = pbuf_get_index('DES',errcode=err)
   i_icswp  = pbuf_get_index('ICSWP',errcode=err)
   i_rei    = pbuf_get_index('REI',errcode=err)  ! added by UM team on Dec.18
   i_cld    = pbuf_get_index('CLD',errcode=err)  ! added by UM team on Dec.18

   ! old optics
   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('CLDLIQ', ixcldliq)

   ! read liquid cloud optics
   if(masterproc) then
   call getfil( trim(liquidfile), locfn, 0)
   call handle_ncerr( nf90_open(locfn, NF90_NOWRITE, ncid), 'liquid optics file missing')
   write(iulog,*)' reading liquid cloud optics from file ',locfn

   call handle_ncerr(nf90_inq_dimid( ncid, 'lw_band', dimid), 'getting lw_band dim')
   call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=f_nlwbands), 'getting n lw bands')
   if (f_nlwbands /= nlwbands) call endrun('number of lw bands does not match')

   call handle_ncerr(nf90_inq_dimid( ncid, 'sw_band', dimid), 'getting sw_band_dim')
   call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=f_nswbands), 'getting n sw bands')
   if (f_nswbands /= nswbands) call endrun('number of sw bands does not match')

   call handle_ncerr(nf90_inq_dimid( ncid, 'mu', mudimid), 'getting mu dim')
   call handle_ncerr(nf90_inquire_dimension( ncid, mudimid, len=nmu), 'getting n mu samples')

   call handle_ncerr(nf90_inq_dimid( ncid, 'lambda_scale', lambdadimid), 'getting lambda dim')
   call handle_ncerr(nf90_inquire_dimension( ncid, lambdadimid, len=nlambda), 'getting n lambda samples')
   endif ! if (masterproc)

#if ( defined SPMD )
   call mpibcast(nmu, 1, mpiint, 0, mpicom, ierr)
   call mpibcast(nlambda, 1, mpiint, 0, mpicom, ierr)
#endif

   allocate(g_mu(nmu))
   allocate(g_lambda(nmu,nlambda))
   allocate(ext_sw_liq(nmu,nlambda,nswbands) )
   allocate(ssa_sw_liq(nmu,nlambda,nswbands))
   allocate(asm_sw_liq(nmu,nlambda,nswbands))
   allocate(abs_lw_liq(nmu,nlambda,nlwbands))

   if(masterproc) then
   call handle_ncerr( nf90_inq_varid(ncid, 'mu', mu_id),&
      'cloud optics mu get')
   call handle_ncerr( nf90_get_var(ncid, mu_id, g_mu),&
      'read cloud optics mu values')

   call handle_ncerr( nf90_inq_varid(ncid, 'lambda', lambda_id),&
      'cloud optics lambda get')
   call handle_ncerr( nf90_get_var(ncid, lambda_id, g_lambda),&
      'read cloud optics lambda values')

   call handle_ncerr( nf90_inq_varid(ncid, 'k_ext_sw', ext_sw_liq_id),&
      'cloud optics ext_sw_liq get')
   call handle_ncerr( nf90_get_var(ncid, ext_sw_liq_id, ext_sw_liq),&
      'read cloud optics ext_sw_liq values')

   call handle_ncerr( nf90_inq_varid(ncid, 'ssa_sw', ssa_sw_liq_id),&
      'cloud optics ssa_sw_liq get')
   call handle_ncerr( nf90_get_var(ncid, ssa_sw_liq_id, ssa_sw_liq),&
      'read cloud optics ssa_sw_liq values')

   call handle_ncerr( nf90_inq_varid(ncid, 'asm_sw', asm_sw_liq_id),&
      'cloud optics asm_sw_liq get')
   call handle_ncerr( nf90_get_var(ncid, asm_sw_liq_id, asm_sw_liq),&
      'read cloud optics asm_sw_liq values')

   call handle_ncerr( nf90_inq_varid(ncid, 'k_abs_lw', abs_lw_liq_id),&
      'cloud optics abs_lw_liq get')
   call handle_ncerr( nf90_get_var(ncid, abs_lw_liq_id, abs_lw_liq),&
      'read cloud optics abs_lw_liq values')

   call handle_ncerr( nf90_close(ncid), 'liquid optics file missing')
   endif ! if masterproc

#if ( defined SPMD )
    call mpibcast(g_mu, nmu, mpir8, 0, mpicom, ierr)
    call mpibcast(g_lambda, nmu*nlambda, mpir8, 0, mpicom, ierr)
    call mpibcast(ext_sw_liq, nmu*nlambda*nswbands, mpir8, 0, mpicom, ierr)
    call mpibcast(ssa_sw_liq, nmu*nlambda*nswbands, mpir8, 0, mpicom, ierr)
    call mpibcast(asm_sw_liq, nmu*nlambda*nswbands, mpir8, 0, mpicom, ierr)
    call mpibcast(abs_lw_liq, nmu*nlambda*nlwbands, mpir8, 0, mpicom, ierr)
#endif
   ! I forgot to convert kext from m^2/Volume to m^2/Kg
   ext_sw_liq = ext_sw_liq / 0.9970449e3_r8 
   abs_lw_liq = abs_lw_liq / 0.9970449e3_r8 

   ! read ice cloud optics
   if(masterproc) then
   call getfil( trim(icefile), locfn, 0)
   call handle_ncerr( nf90_open(locfn, NF90_NOWRITE, ncid), 'ice optics file missing')
   write(iulog,*)' reading ice cloud optics from file ',locfn

   call handle_ncerr(nf90_inq_dimid( ncid, 'lw_band', dimid), 'getting lw_band dim')
   call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=f_nlwbands), 'getting n lw bands')
   if (f_nlwbands /= nlwbands) call endrun('number of lw bands does not match')

   call handle_ncerr(nf90_inq_dimid( ncid, 'sw_band', dimid), 'getting sw_band_dim')
   call handle_ncerr(nf90_inquire_dimension( ncid, dimid, len=f_nswbands), 'getting n sw bands')
   if (f_nswbands /= nswbands) call endrun('number of sw bands does not match')

   call handle_ncerr(nf90_inq_dimid( ncid, 'd_eff', d_dimid), 'getting deff dim')
   call handle_ncerr(nf90_inquire_dimension( ncid, d_dimid, len=n_g_d), 'getting n deff samples')

   endif ! if (masterproc)

#if ( defined SPMD )
   call mpibcast(n_g_d, 1, mpiint, 0, mpicom, ierr)
!   call mpibcast(nswbands, 1, mpiint, 0, mpicom, ierr)
!   call mpibcast(nlwbands, 1, mpiint, 0, mpicom, ierr)
#endif

   allocate(g_d_eff(n_g_d))
   allocate(ext_sw_ice(n_g_d,nswbands))
   allocate(ssa_sw_ice(n_g_d,nswbands))
   allocate(asm_sw_ice(n_g_d,nswbands))
   allocate(abs_lw_ice(n_g_d,nlwbands))

   if(masterproc) then
   call handle_ncerr( nf90_inq_varid(ncid, 'd_eff', d_id),&
      'cloud optics deff get')
   call handle_ncerr( nf90_get_var(ncid, d_id, g_d_eff),&
      'read cloud optics deff values')

   call handle_ncerr( nf90_inq_varid(ncid, 'sw_ext', ext_sw_ice_id),&
      'cloud optics ext_sw_ice get')
   call handle_ncerr(nf90_inquire_variable ( ncid, ext_sw_ice_id, ndims=ndims, dimids=vdimids),&
       'checking dimensions of ext_sw_ice')
   call handle_ncerr(nf90_inquire_dimension( ncid, vdimids(1), len=templen),&
       'getting first dimension sw_ext')
   !write(iulog,*) 'expected length',n_g_d,'actual len',templen
   call handle_ncerr(nf90_inquire_dimension( ncid, vdimids(2), len=templen),&
       'getting first dimension sw_ext')
   !write(iulog,*) 'expected length',nswbands,'actual len',templen
   call handle_ncerr( nf90_get_var(ncid, ext_sw_ice_id, ext_sw_ice),&
      'read cloud optics ext_sw_ice values')

   call handle_ncerr( nf90_inq_varid(ncid, 'sw_ssa', ssa_sw_ice_id),&
      'cloud optics ssa_sw_ice get')
   call handle_ncerr( nf90_get_var(ncid, ssa_sw_ice_id, ssa_sw_ice),&
      'read cloud optics ssa_sw_ice values')

   call handle_ncerr( nf90_inq_varid(ncid, 'sw_asm', asm_sw_ice_id),&
      'cloud optics asm_sw_ice get')
   call handle_ncerr( nf90_get_var(ncid, asm_sw_ice_id, asm_sw_ice),&
      'read cloud optics asm_sw_ice values')

   call handle_ncerr( nf90_inq_varid(ncid, 'lw_abs', abs_lw_ice_id),&
      'cloud optics abs_lw_ice get')
   call handle_ncerr( nf90_get_var(ncid, abs_lw_ice_id, abs_lw_ice),&
      'read cloud optics abs_lw_ice values')

   call handle_ncerr( nf90_close(ncid), 'ice optics file missing')

   endif ! if masterproc
#if ( defined SPMD )
    call mpibcast(g_d_eff, n_g_d, mpir8, 0, mpicom, ierr)
    call mpibcast(ext_sw_ice, n_g_d*nswbands, mpir8, 0, mpicom, ierr)
    call mpibcast(ssa_sw_ice, n_g_d*nswbands, mpir8, 0, mpicom, ierr)
    call mpibcast(asm_sw_ice, n_g_d*nswbands, mpir8, 0, mpicom, ierr)
    call mpibcast(abs_lw_ice, n_g_d*nlwbands, mpir8, 0, mpicom, ierr)
#endif

    return

end subroutine cloud_rad_props_init

!==============================================================================

subroutine cloud_rad_props_get_sw(state, pbuf, &
                                  tau, tau_w, tau_w_g, tau_w_f,&
                                  diagnosticindex, oldliq, oldice)

! return totaled (across all species) layer tau, omega, g, f 
! for all spectral interval for aerosols affecting the climate

   ! Arguments
   type(physics_state), intent(in)   :: state
   type(physics_buffer_desc),pointer :: pbuf(:)
   integer, optional,   intent(in) :: diagnosticindex      ! index (if present) to radiation diagnostic information

   real(r8), intent(out) :: tau    (nswbands,pcols,pver) ! aerosol extinction optical depth
   real(r8), intent(out) :: tau_w  (nswbands,pcols,pver) ! aerosol single scattering albedo * tau
   real(r8), intent(out) :: tau_w_g(nswbands,pcols,pver) ! aerosol assymetry parameter * tau * w
   real(r8), intent(out) :: tau_w_f(nswbands,pcols,pver) ! aerosol forward scattered fraction * tau * w

   logical, optional, intent(in) :: oldliq,oldice

   ! Local variables

   integer :: ncol
   integer :: lchnk
   integer :: k, i    ! lev and daycolumn indices
   integer :: iswband ! sw band indices

   ! optical props for each aerosol
   real(r8), pointer :: h_ext(:,:)
   real(r8), pointer :: h_ssa(:,:)
   real(r8), pointer :: h_asm(:,:)
   real(r8), pointer :: n_ext(:)
   real(r8), pointer :: n_ssa(:)
   real(r8), pointer :: n_asm(:)

   ! rad properties for liquid clouds
   real(r8) :: liq_tau    (nswbands,pcols,pver) ! aerosol extinction optical depth
   real(r8) :: liq_tau_w  (nswbands,pcols,pver) ! aerosol single scattering albedo * tau
   real(r8) :: liq_tau_w_g(nswbands,pcols,pver) ! aerosol assymetry parameter * tau * w
   real(r8) :: liq_tau_w_f(nswbands,pcols,pver) ! aerosol forward scattered fraction * tau * w

   ! rad properties for ice clouds
   real(r8) :: ice_tau    (nswbands,pcols,pver) ! aerosol extinction optical depth
   real(r8) :: ice_tau_w  (nswbands,pcols,pver) ! aerosol single scattering albedo * tau
   real(r8) :: ice_tau_w_g(nswbands,pcols,pver) ! aerosol assymetry parameter * tau * w
   real(r8) :: ice_tau_w_f(nswbands,pcols,pver) ! aerosol forward scattered fraction * tau * w

   !-----------------------------------------------------------------------------

   ncol  = state%ncol
   lchnk = state%lchnk

   ! initialize to conditions that would cause failure
   tau     (:,:,:) = -100._r8
   tau_w   (:,:,:) = -100._r8
   tau_w_g (:,:,:) = -100._r8
   tau_w_f (:,:,:) = -100._r8

   ! initialize layers to accumulate od's
   tau    (:,1:ncol,:) = 0._r8
   tau_w  (:,1:ncol,:) = 0._r8
   tau_w_g(:,1:ncol,:) = 0._r8
   tau_w_f(:,1:ncol,:) = 0._r8


   call get_liquid_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f)

   call get_ice_optics_sw   (state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f)

   tau    (:,1:ncol,:) =  liq_tau    (:,1:ncol,:) + ice_tau    (:,1:ncol,:)
   tau_w  (:,1:ncol,:) =  liq_tau_w  (:,1:ncol,:) + ice_tau_w  (:,1:ncol,:)
   tau_w_g(:,1:ncol,:) =  liq_tau_w_g(:,1:ncol,:) + ice_tau_w_g(:,1:ncol,:)
   tau_w_f(:,1:ncol,:) =  liq_tau_w_f(:,1:ncol,:) + ice_tau_w_f(:,1:ncol,:)

end subroutine cloud_rad_props_get_sw
!==============================================================================

subroutine cloud_rad_props_get_lw(state, pbuf, cld_abs_od, diagnosticindex, oldliq, oldice, oldcloud)

! Purpose: Compute cloud longwave absorption optical depth
!    cloud_rad_props_get_lw() is called by radlw() 

   ! Arguments
   type(physics_state), intent(in)  :: state
   type(physics_buffer_desc),pointer:: pbuf(:)
   real(r8),            intent(out) :: cld_abs_od(nlwbands,pcols,pver) ! [fraction] absorption optical depth, per layer
   integer, optional,   intent(in)  :: diagnosticindex
   logical, optional,   intent(in)  :: oldliq  ! use old liquid optics
   logical, optional,   intent(in)  :: oldice  ! use old ice optics
   logical, optional,   intent(in)  :: oldcloud  ! use old optics for both (b4b)

   ! Local variables

   integer :: bnd_idx     ! LW band index
   integer :: i           ! column index
   integer :: k           ! lev index
   integer :: ncol        ! number of columns
   integer :: lchnk

   ! rad properties for liquid clouds
   real(r8) :: liq_tau_abs_od(nlwbands,pcols,pver) ! liquid cloud absorption optical depth

   ! rad properties for ice clouds
   real(r8) :: ice_tau_abs_od(nlwbands,pcols,pver) ! ice cloud absorption optical depth

   !-----------------------------------------------------------------------------

   ncol = state%ncol
   lchnk = state%lchnk

   ! compute optical depths cld_absod 
   cld_abs_od = 0._r8

   if(present(oldcloud))then
      if(oldcloud) then
         ! make diagnostic calls to these first to output ice and liq OD's
         !call old_liq_get_rad_props_lw(state, pbuf, liq_tau_abs_od, oldliqwp=.false.)
         !call old_ice_get_rad_props_lw(state, pbuf, ice_tau_abs_od, oldicewp=.false.)
         ! This affects climate (cld_abs_od)
         call oldcloud_lw(state,pbuf,cld_abs_od,oldwp=.false.)
         return
      endif
   endif

   if(present(oldliq))then
      if(oldliq) then
         call old_liq_get_rad_props_lw(state, pbuf, liq_tau_abs_od, oldliqwp=.false.)
      else
         call liquid_cloud_get_rad_props_lw(state, pbuf, liq_tau_abs_od)
      endif
   else
      call liquid_cloud_get_rad_props_lw(state, pbuf, liq_tau_abs_od)
   endif

   if(present(oldice))then
      if(oldice) then
         call old_ice_get_rad_props_lw(state, pbuf, ice_tau_abs_od, oldicewp=.false.)
      else
         call ice_cloud_get_rad_props_lw(state, pbuf, ice_tau_abs_od)
      endif
   else
      call ice_cloud_get_rad_props_lw(state, pbuf, ice_tau_abs_od)
   endif
      
   cld_abs_od(:,1:ncol,:) = liq_tau_abs_od(:,1:ncol,:) + ice_tau_abs_od(:,1:ncol,:) 

end subroutine cloud_rad_props_get_lw

!==============================================================================

subroutine get_snow_optics_sw(state, pbuf, tau, tau_w, tau_w_g, tau_w_f)
   type(physics_state), intent(in)   :: state
   type(physics_buffer_desc),pointer :: pbuf(:)

   real(r8),intent(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
   real(r8),intent(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
   real(r8),intent(out) :: tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
   real(r8),intent(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

   real(r8), pointer :: icswpth(:,:), des(:,:)

   ! This does the same thing as get_ice_optics_sw, except with a different
   ! water path and effective diameter.
   call pbuf_get_field(pbuf, i_icswp, icswpth)
   call pbuf_get_field(pbuf, i_des,   des)

   call interpolate_ice_optics_sw(state%ncol, icswpth, des, tau, tau_w, &
        tau_w_g, tau_w_f)

end subroutine get_snow_optics_sw   

!==============================================================================
! Private methods
!==============================================================================

subroutine get_ice_optics_sw(state, pbuf, tau, tau_w, tau_w_g, tau_w_f)
   type(physics_state), intent(in)   :: state
   type(physics_buffer_desc),pointer :: pbuf(:)

   real(r8),intent(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
   real(r8),intent(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
   real(r8),intent(out) :: tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
   real(r8),intent(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

   real(r8), pointer :: iciwpth(:,:), dei(:,:)

   ! Get relevant pbuf fields, and interpolate optical properties from
   ! the lookup tables.
   call pbuf_get_field(pbuf, i_iciwp, iciwpth)
   call pbuf_get_field(pbuf, i_dei,   dei)

   call interpolate_ice_optics_sw(state%ncol, iciwpth, dei, tau, tau_w, &
        tau_w_g, tau_w_f)

end subroutine get_ice_optics_sw

!==============================================================================

subroutine interpolate_ice_optics_sw(ncol, iciwpth, dei, tau, tau_w, &
     tau_w_g, tau_w_f)

  integer, intent(in) :: ncol
  real(r8), intent(in) :: iciwpth(pcols,pver)
  real(r8), intent(in) :: dei(pcols,pver)

  real(r8),intent(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
  real(r8),intent(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
  real(r8),intent(out) :: tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
  real(r8),intent(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

  type(interp_type) :: dei_wgts

  integer :: i, k, swband
  real(r8) :: ext(nswbands), ssa(nswbands), asm(nswbands)

  do k = 1,pver
     do i = 1,ncol
        if( iciwpth(i,k) < 1.e-80_r8 .or. dei(i,k) == 0._r8) then
           ! if ice water path is too small, OD := 0
           tau    (:,i,k) = 0._r8
           tau_w  (:,i,k) = 0._r8
           tau_w_g(:,i,k) = 0._r8
           tau_w_f(:,i,k) = 0._r8
        else
           ! for each cell interpolate to find weights in g_d_eff grid.
           call lininterp_init(g_d_eff, n_g_d, dei(i:i,k), 1, &
                extrap_method_bndry, dei_wgts)
           ! interpolate into grid and extract radiative properties
           do swband = 1, nswbands
              call lininterp(ext_sw_ice(:,swband), n_g_d, &
                   ext(swband:swband), 1, dei_wgts)
              call lininterp(ssa_sw_ice(:,swband), n_g_d, &
                   ssa(swband:swband), 1, dei_wgts)
              call lininterp(asm_sw_ice(:,swband), n_g_d, &
                   asm(swband:swband), 1, dei_wgts)
           end do
           tau    (:,i,k) = iciwpth(i,k) * ext
           tau_w  (:,i,k) = tau(:,i,k) * ssa
           tau_w_g(:,i,k) = tau_w(:,i,k) * asm
           tau_w_f(:,i,k) = tau_w_g(:,i,k) * asm
           call lininterp_finish(dei_wgts)
        endif
     enddo
  enddo

end subroutine interpolate_ice_optics_sw

!==============================================================================

subroutine get_liquid_optics_sw(state, pbuf, tau, tau_w, tau_w_g, tau_w_f)
   type(physics_state), intent(in)   :: state
   type(physics_buffer_desc),pointer :: pbuf(:)

   real(r8),intent(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
   real(r8),intent(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
   real(r8),intent(out) :: tau_w_g(nswbands,pcols,pver) ! asymetry parameter * tau * w
   real(r8),intent(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

   real(r8), pointer, dimension(:,:) :: lamc, pgam, iclwpth
   real(r8), dimension(pcols,pver) :: kext
   integer i,k,swband,lchnk,ncol

   lchnk = state%lchnk
   ncol = state%ncol


   call pbuf_get_field(pbuf, i_lambda,  lamc)
   call pbuf_get_field(pbuf, i_mu,      pgam)
   call pbuf_get_field(pbuf, i_iclwp,   iclwpth)
   
   do k = 1,pver
      do i = 1,ncol
         if(lamc(i,k) > 0._r8) then ! This seems to be clue from microphysics of no cloud
            call gam_liquid_sw(iclwpth(i,k), lamc(i,k), pgam(i,k), &
                tau(1:nswbands,i,k), tau_w(1:nswbands,i,k), tau_w_g(1:nswbands,i,k), tau_w_f(1:nswbands,i,k))
         else
            tau(1:nswbands,i,k) = 0._r8
            tau_w(1:nswbands,i,k) = 0._r8
            tau_w_g(1:nswbands,i,k) = 0._r8
            tau_w_f(1:nswbands,i,k) = 0._r8
         endif
      enddo
   enddo

end subroutine get_liquid_optics_sw

!==============================================================================

subroutine liquid_cloud_get_rad_props_lw(state, pbuf, abs_od)
   type(physics_state), intent(in)    :: state
   type(physics_buffer_desc),pointer  :: pbuf(:)
   real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)

   integer :: lchnk, ncol
   real(r8), pointer, dimension(:,:) :: lamc, pgam, iclwpth

   integer lwband, i, k

   abs_od = 0._r8

   lchnk = state%lchnk
   ncol = state%ncol

   call pbuf_get_field(pbuf, i_lambda,  lamc)
   call pbuf_get_field(pbuf, i_mu,      pgam)
   call pbuf_get_field(pbuf, i_iclwp,   iclwpth)

   do k = 1,pver
      do i = 1,ncol
         if(lamc(i,k) > 0._r8) then ! This seems to be the clue for no cloud from microphysics formulation
            call gam_liquid_lw(iclwpth(i,k), lamc(i,k), pgam(i,k), abs_od(1:nlwbands,i,k))
         else
            abs_od(1:nlwbands,i,k) = 0._r8
         endif
      enddo
   enddo

end subroutine liquid_cloud_get_rad_props_lw
!==============================================================================

subroutine snow_cloud_get_rad_props_lw(state, pbuf, abs_od)
   type(physics_state), intent(in)    :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)

   real(r8), pointer :: icswpth(:,:), des(:,:)

   ! This does the same thing as ice_cloud_get_rad_props_lw, except with a
   ! different water path and effective diameter.
   call pbuf_get_field(pbuf, i_icswp, icswpth)
   call pbuf_get_field(pbuf, i_des,   des)

   call interpolate_ice_optics_lw(state%ncol,icswpth, des, abs_od)

end subroutine snow_cloud_get_rad_props_lw

!==============================================================================

subroutine ice_cloud_get_rad_props_lw(state, pbuf, abs_od)
   type(physics_state), intent(in)     :: state
   type(physics_buffer_desc), pointer  :: pbuf(:)
   real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)

   real(r8), pointer :: iciwpth(:,:), dei(:,:)

   ! Get relevant pbuf fields, and interpolate optical properties from
   ! the lookup tables.
   call pbuf_get_field(pbuf, i_iciwp, iciwpth)
   call pbuf_get_field(pbuf, i_dei,   dei)

   call interpolate_ice_optics_lw(state%ncol,iciwpth, dei, abs_od)

end subroutine ice_cloud_get_rad_props_lw

!==============================================================================

subroutine interpolate_ice_optics_lw(ncol, iciwpth, dei, abs_od)

  integer, intent(in) :: ncol
  real(r8), intent(in) :: iciwpth(pcols,pver)
  real(r8), intent(in) :: dei(pcols,pver)

  real(r8),intent(out) :: abs_od(nlwbands,pcols,pver)

  type(interp_type) :: dei_wgts

  integer :: i, k, lwband
  real(r8) :: absor(nlwbands)

  do k = 1,pver
     do i = 1,ncol
        ! if ice water path is too small, OD := 0
        if( iciwpth(i,k) < 1.e-80_r8 .or. dei(i,k) == 0._r8) then
           abs_od (:,i,k) = 0._r8
        else
           ! for each cell interpolate to find weights in g_d_eff grid.
           call lininterp_init(g_d_eff, n_g_d, dei(i:i,k), 1, &
                extrap_method_bndry, dei_wgts)
           ! interpolate into grid and extract radiative properties
           do lwband = 1, nlwbands
              call lininterp(abs_lw_ice(:,lwband), n_g_d, &
                   absor(lwband:lwband), 1, dei_wgts)
           enddo
           abs_od(:,i,k) = iciwpth(i,k) * absor
!!== KZ_BUGFIX_20161011 
!!         Limit the snow/ice optical depths fix the case where there are large masses 
!!         of snow at high altitudes
           where(abs_od(:,i,k) > 50.0_r8) abs_od(:,i,k) = 50.0_r8
!!== KZ_BUGFIX_20161011 
           call lininterp_finish(dei_wgts)
        endif
     enddo
  enddo

end subroutine interpolate_ice_optics_lw

!==============================================================================

subroutine gam_liquid_lw(clwptn, lamc, pgam, abs_od)
  real(r8), intent(in) :: clwptn ! cloud water liquid path new (in cloud) (in g/m^2)?
  real(r8), intent(in) :: lamc   ! prognosed value of lambda for cloud
  real(r8), intent(in) :: pgam   ! prognosed value of mu for cloud
  real(r8), intent(out) :: abs_od(1:nlwbands)

  integer :: lwband ! sw band index

  type(interp_type) :: mu_wgts
  type(interp_type) :: lambda_wgts

  if (clwptn < 1.e-80_r8) then
    abs_od = 0._r8
    return
  endif

  call get_mu_lambda_weights(lamc, pgam, mu_wgts, lambda_wgts)

  do lwband = 1, nlwbands
     call lininterp(abs_lw_liq(:,:,lwband), nmu, nlambda, &
          abs_od(lwband:lwband), 1, mu_wgts, lambda_wgts)
  enddo

  abs_od = clwptn * abs_od

  call lininterp_finish(mu_wgts)
  call lininterp_finish(lambda_wgts)

end subroutine gam_liquid_lw

!==============================================================================

subroutine gam_liquid_sw(clwptn, lamc, pgam, tau, tau_w, tau_w_g, tau_w_f)
  real(r8), intent(in) :: clwptn ! cloud water liquid path new (in cloud) (in g/m^2)?
  real(r8), intent(in) :: lamc   ! prognosed value of lambda for cloud
  real(r8), intent(in) :: pgam   ! prognosed value of mu for cloud
  real(r8), intent(out) :: tau(1:nswbands), tau_w(1:nswbands), tau_w_f(1:nswbands), tau_w_g(1:nswbands)

  integer :: swband ! sw band index

  real(r8) :: ext(nswbands), ssa(nswbands), asm(nswbands)

  type(interp_type) :: mu_wgts
  type(interp_type) :: lambda_wgts

  if (clwptn < 1.e-80_r8) then
    tau = 0._r8
    tau_w = 0._r8
    tau_w_g = 0._r8
    tau_w_f = 0._r8
    return
  endif

  call get_mu_lambda_weights(lamc, pgam, mu_wgts, lambda_wgts)

  do swband = 1, nswbands
     call lininterp(ext_sw_liq(:,:,swband), nmu, nlambda, &
          ext(swband:swband), 1, mu_wgts, lambda_wgts)
     call lininterp(ssa_sw_liq(:,:,swband), nmu, nlambda, &
          ssa(swband:swband), 1, mu_wgts, lambda_wgts)
     call lininterp(asm_sw_liq(:,:,swband), nmu, nlambda, &
          asm(swband:swband), 1, mu_wgts, lambda_wgts)
  enddo

  ! compute radiative properties
  tau = clwptn * ext
  tau_w = tau * ssa
  tau_w_g = tau_w * asm
  tau_w_f = tau_w_g * asm

  call lininterp_finish(mu_wgts)
  call lininterp_finish(lambda_wgts)

end subroutine gam_liquid_sw

!==============================================================================

subroutine get_mu_lambda_weights(lamc, pgam, mu_wgts, lambda_wgts)
  real(r8), intent(in) :: lamc   ! prognosed value of lambda for cloud
  real(r8), intent(in) :: pgam   ! prognosed value of mu for cloud
  ! Output interpolation weights. Caller is responsible for freeing these.
  type(interp_type), intent(out) :: mu_wgts
  type(interp_type), intent(out) :: lambda_wgts

  integer :: ilambda
  real(r8) :: g_lambda_interp(nlambda)

  ! Make interpolation weights for mu.
  ! (Put pgam in a temporary array for this purpose.)
  call lininterp_init(g_mu, nmu, [pgam], 1, extrap_method_bndry, mu_wgts)

  ! Use mu weights to interpolate to a row in the lambda table.
  do ilambda = 1, nlambda
     call lininterp(g_lambda(:,ilambda), nmu, &
          g_lambda_interp(ilambda:ilambda), 1, mu_wgts)
  end do

  ! Make interpolation weights for lambda.
  call lininterp_init(g_lambda_interp, nlambda, [lamc], 1, &
       extrap_method_bndry, lambda_wgts)

end subroutine get_mu_lambda_weights

!============================
subroutine mc6_ice_get_rad_props_lw(state, pbuf, ext_od, abs_od, ssa_od, xmomc_od, oldicewp)
   use physconst,      only: gravit

!-----------------------------------------------------
! added by UM team on Dec.18
!
! Description:
!   MC6 ice cloud optics scheme. Ref: Kuo et al. (2020, JQSRT, https://doi.org/10.1016/j.jqsrt.2019.106683)
!
! History:
!   2016/05/20  ver 1.0 , ADD
!   2017/08/03  ver 1.1 , ADD extinction optical depth and absorption optical depth
!   2017/09/07  ver 2.1 , DEBUG use REI, ice effective radius, directly from CESM becasue it's already times a factor of 1.5 
!   2017/09/19  ver 2.2 , DEBUG iciwp, previous version use inconsistent iciwp with the Mitchell one, and the unit is also wrong
!
! Author:
!   Yi-Hsuan Chen
!   yihsuan@umich.edu
!-----------------------------------------------------

   use rrlw_cld, only: abscld1, absliq0, absliq1, &
                       absice0, absice1, absice2, absice3, &
                       absice4, extice4, ssaice4, asyice4, &! CPKuo@TAMU
                       absice5, extice5, ssaice5, asyice5, &! CPKuo@TAMU
                       absice6, extice6, ssaice6, asyice6   ! CPKuo@TAMU

!--- Input arguments ---!
   type(physics_state), intent(in)     :: state
   type(physics_buffer_desc), pointer  :: pbuf(:)
   logical , intent(in) :: oldicewp

   real(r8), pointer, dimension(:,:) :: rei     ! CAM ice effective radius from mg scheme (microns)
   real(r8), pointer, dimension(:,:) :: dei     ! CAM ice effective diameter. 2020/02/10, Xianwen Jing.
   real(r8), pointer, dimension(:,:) :: iciwpth ! CAM in-cloud ice water path (kg/m2)
   real(r8), pointer, dimension(:,:) :: cldn    ! CAM cloud fraction

!--- Output arguments ---!
   real(r8), intent(out) :: ext_od  (nlwbands,pcols,pver)        ! cloud ice extinction optical depth, i.e. including absorption and scattering
   real(r8), intent(out) :: abs_od  (nlwbands,pcols,pver)        ! cloud ice absorption optical depth
   real(r8), intent(out) :: ssa_od  (nlwbands,pcols,pver)        ! single scattering albedo
   !real(r8), intent(out) :: xmomc_od(0:16,nlwbands,pcols,pver)   ! phase function
   real(r8), intent(out) :: xmomc_od(0:1,nlwbands,pcols,pver)   ! phase function, use (0:16) edison will have segmentation fault

!--- Local ---!
   real(r8) :: diaice                        ! cloud ice effective diameter (microns)
   real(r8) :: invrad 
   real(r8) :: extcoice(nlwbands)            ! ice mass-extinction coefficients of TAMU scheme (m2/g)
   real(r8) :: ssacoice(nlwbands)            ! ice single scattering albedo of TAMU scheme (unitless)
   real(r8) :: asycoice(nlwbands)            ! ice asymmetric factor of TAMU scheme (unitless)
   real(r8) :: ext_tamu(nlwbands,pcols,pver) ! ice cloud extinction optical thickness (unitless)
   real(r8) :: abs_tamu(nlwbands,pcols,pver) ! ice cloud absorption optical thickness (unitless)
   real(r8) :: ssa_tamu(nlwbands,pcols,pver) ! ice cloud single scattering albedo (unitless)
   real(r8) :: xmomc_tamu(0:1,nlwbands,pcols,pver) ! ice cloud asymmetric factor (unitless), use (0:16) edison will have segmentation fault
   real(r8) :: iciwp(pcols,pver)             ! work array of in-cloud ice water path (kg/m2)
   real(r8) :: taucloud, ssacloud, asycloud, xmomcloud(0:1), fpeak

   real(r8) :: gicewp(pcols,pver)            ! work array of grid-box ice water path (kg/m2)
   real(r8) :: cicewp(pcols,pver)            ! work array of in-cloud ice water path (kg/m2)    

   integer  :: icb(nlwbands,0:2)
   data icb /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, &
             1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5, &
             1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/

   logical :: option_pbuf

   integer :: rei_idx, dei_idx, iciwp_idx, err, itim !add dei_idx. 2020/02/10, Xianwen Jing.
   integer,parameter :: iceflag=3
       !  iceflag= 1 : MC6 ice cloud optical properties (MC6 ice crystal shape and 
       !                  0.1 variance gamma PSD, absorption and scattering)
       !           2 : THM ice cloud optical properties (Two-Habit model ice crystal 
       !                  shape and 0.1 variance gamma PSD, absorption and scattering)
       !           3 : TAMU MC6 ice cloud optical properties (MC6 ice crystal shape 
       !                  and Bryan Baum in-situ PSD, absorption and scattering)   
   real(r8) :: diaice_max        !uplimit of DEI for ice rad properties, Xianwen Jing.
   integer :: ib,i,k,ncol,lchnk,iceind

!----------------------------
! initialize return arrays
!----------------------------
   ext_tamu   = 0._r8
   abs_tamu   = 0._r8
   ssa_tamu   = 0._r8
   xmomc_tamu = 0._r8

   ext_od   = 0._r8
   abs_od   = 0._r8
   ssa_od   = 0._r8
   xmomc_od = 0._r8

   ncol = state%ncol
   lchnk = state%lchnk

!-----------------
! get CAM fields
!-----------------

!*** ice water path ***
!
!note - UM team Dec.18, 2019  
! if ICIWP is not in cam buffer. So, pbuf_get_field(pbuf,iciwp_idx, iciwpth) would fail. 
! the error message is 
!   ENDRUN:index out of range
!   /home/UM team/model/cesm1_1_1/models/atm/cam/src/physics/cam/physics_buffer.F90
!            1171           -1
!  
   if(oldicewp) then   ! if use ebertcurry, ICIWP was not found in pbuf. So, ICIWP is calculated here.
     itim = pbuf_old_tim_idx()
     call pbuf_get_field(pbuf, i_cld,   cldn, start=(/1,1,itim/),kount=(/pcols,pver,1/))

     do k=1,pver
         do i = 1,ncol
            gicewp(i,k) = state%q(i,k,ixcldice)*state%pdel(i,k)/gravit ! Grid box ice water path.
            cicewp(i,k) = gicewp(i,k) / max(0.01_r8,cldn(i,k))         ! In-cloud ice water path.
         end do
     end do

     iciwp(:,:) = cicewp(:,:)

     call pbuf_get_field(pbuf, i_rei, rei)

   else

     ! in-cloud ice water path (kg/m2)
     call pbuf_get_field(pbuf, i_iciwp, iciwpth)
     iciwp(:,:) = iciwpth(:,:)

     !*** ice effective radius ***
     !rei_idx    = pbuf_get_index('REI',errcode=err)  ! added by UM team on Dec.18
     !call pbuf_get_field(pbuf, rei_idx  , rei)

     !*** ice effective diameter ***
     dei_idx    = pbuf_get_index('DEI',errcode=err)  ! added by UM team on Feb. 10, 2020.
     call pbuf_get_field(pbuf, dei_idx  , dei)

   end if

!--------------------------------------------------------
! compute absorption coefficient and cloud optical depth
!--------------------------------------------------------
   do i = 1,ncol
      do k = 1,pver
        !diaice = 2._r8*rei(i,k) ! effective diameter is defined by that particle-volume 
                                 ! divided by particle-projected-area (i.e. effective diameter) times 1.5.
                                 ! Chia-Pang Kuo @ TAMU, personal communication
                                 ! because the REI in CESM is already times 1.5, so don't need to times 1.5 here.
                                 ! (ref: eq. 4.154 in CAM5 scientific description)
        diaice = dei(i,k)  ! effective diameter defined by microphysics (Xianwen)
        extcoice(:) = 0._r8
        ssacoice(:) = 0._r8
        asycoice(:) = 0._r8

        !*** if there is no ice cloud ***
        if( iciwp(i,k) < 1.e-80_r8 .or. diaice .eq. 0._r8) then
          ext_tamu (:,i,k) = 0._r8
          abs_tamu (:,i,k) = 0._r8
          ssa_tamu (:,i,k) = 0._r8
          xmomc_tamu (:,:,i,k) = 0._r8
          iceind = 0

        !*** if there is an ice cloud layer ***
        else

          !diaice = min(diaice, 370._r8)  ! upper bound of ice effective diameter of TAMU scheme is 370 microns
                                         ! because mass absorption coefficient larger than 370um vary slightly with diameter,
                                         ! apply mass absorption coefficient of 370um to those larger ones is not a bad approximation
                                         ! Chia-Pang Kuo @ TAMU, personal communication

          ! if ice effective diameter in approciate range (3-370 microns)

        !-- set upper bound according to ice optical option -->
        !   iceflag=1. MC6 scattering + 0.1 variance Gamma PSD
        !           2. THM scattering + 0.1 variance Gamma PSD
        !           3. TAMU MC6 scattering + Bryan Baum in-situ PSD
          if (iceflag==1) then    
             diaice_max=500._r8
          else if (iceflag==2) then
             diaice_max=500._r8
          else if (iceflag==3) then
             diaice_max=370._r8
          else 
            call endrun('ERROR in UMRad: current iceflag is not available.')
          end if

          diaice = min(diaice, diaice_max)

          if (diaice .ge. 3.0_r8 .and. diaice .le. diaice_max) then
             if (iceflag==1) then
                invrad = 1.0/diaice - 0.05_r8
                do ib=1,nlwbands           
                   if (diaice .ge. 20.0_r8) then
                      extcoice(ib) = ((extice4(1,1,ib) * invrad + &
                                       extice4(1,2,ib)) * invrad + &
                                       extice4(1,3,ib)) * invrad + &
                                       extice4(1,4,ib)
                      ssacoice(ib) = ((ssaice4(1,1,ib) * invrad + &
                                       ssaice4(1,2,ib)) * invrad + &
                                       ssaice4(1,3,ib)) * invrad + &
                                       ssaice4(1,4,ib)
                      asycoice(ib) = ((asyice4(1,1,ib) * invrad + &
                                       asyice4(1,2,ib)) * invrad + &
                                       asyice4(1,3,ib)) * invrad + &
                                       asyice4(1,4,ib)
                   else
                      extcoice(ib) = (((((extice4(2,1,ib) * invrad + &
                                           extice4(2,2,ib)) * invrad + &
                                           extice4(2,3,ib)) * invrad + &
                                           extice4(2,4,ib)) * invrad + &
                                           extice4(2,5,ib)) * invrad + &
                                           extice4(2,6,ib)) * invrad + &
                                           extice4(2,7,ib)
                      ssacoice(ib) = (((((ssaice4(2,1,ib) * invrad + &
                                           ssaice4(2,2,ib)) * invrad + &
                                           ssaice4(2,3,ib)) * invrad + &
                                           ssaice4(2,4,ib)) * invrad + &
                                           ssaice4(2,5,ib)) * invrad + &
                                           ssaice4(2,6,ib)) * invrad + &
                                           ssaice4(2,7,ib)
                      asycoice(ib) = (((((asyice4(2,1,ib) * invrad + &
                                           asyice4(2,2,ib)) * invrad + &
                                           asyice4(2,3,ib)) * invrad + &
                                           asyice4(2,4,ib)) * invrad + &
                                           asyice4(2,5,ib)) * invrad + &
                                           asyice4(2,6,ib)) * invrad + &
                                           asyice4(2,7,ib)
                   endif  ! end if of diaice .ge. 25.
                 enddo    ! end do of lwbands for cloud radiative coefficients
             else if (iceflag==2) then
                 invrad = 1.0/diaice - 0.05_r8
                 do ib=1,nlwbands           
                   if (diaice .ge. 20.0_r8) then
                      extcoice(ib) = ((extice5(1,1,ib) * invrad + &
                                       extice5(1,2,ib)) * invrad + &
                                       extice5(1,3,ib)) * invrad + &
                                       extice5(1,4,ib)
                      ssacoice(ib) = ((ssaice5(1,1,ib) * invrad + &
                                       ssaice5(1,2,ib)) * invrad + &
                                       ssaice5(1,3,ib)) * invrad + &
                                       ssaice5(1,4,ib)
                      asycoice(ib) = ((asyice5(1,1,ib) * invrad + &
                                       asyice5(1,2,ib)) * invrad + &
                                       asyice5(1,3,ib)) * invrad + &
                                       asyice5(1,4,ib)

                   else
                      extcoice(ib) = (((((extice5(2,1,ib) * invrad + &
                                           extice5(2,2,ib)) * invrad + &
                                           extice5(2,3,ib)) * invrad + &
                                           extice5(2,4,ib)) * invrad + &
                                           extice5(2,5,ib)) * invrad + &
                                           extice5(2,6,ib)) * invrad + &
                                           extice5(2,7,ib)
                      ssacoice(ib) = (((((ssaice5(2,1,ib) * invrad + &
                                           ssaice5(2,2,ib)) * invrad + &
                                           ssaice5(2,3,ib)) * invrad + &
                                           ssaice5(2,4,ib)) * invrad + &
                                           ssaice5(2,5,ib)) * invrad + &
                                           ssaice5(2,6,ib)) * invrad + &
                                           ssaice5(2,7,ib)
                      asycoice(ib) = (((((asyice5(2,1,ib) * invrad + &
                                           asyice5(2,2,ib)) * invrad + &
                                           asyice5(2,3,ib)) * invrad + &
                                           asyice5(2,4,ib)) * invrad + &
                                           asyice5(2,5,ib)) * invrad + &
                                           asyice5(2,6,ib)) * invrad + &
                                           asyice5(2,7,ib)
                   endif  ! end if of diaice .ge. 20.
                 enddo    ! end do of lwbands for cloud radiative coefficients
             else if (iceflag==3) then
                 invrad = 1.0/diaice - 0.04_r8
                 do ib=1,nlwbands           
                   if (diaice .ge. 25.0_r8) then
                      extcoice(ib) = ((extice6(1,1,ib) * invrad + &
                                       extice6(1,2,ib)) * invrad + &
                                       extice6(1,3,ib)) * invrad + &
                                       extice6(1,4,ib)
                      ssacoice(ib) = ((ssaice6(1,1,ib) * invrad + &
                                       ssaice6(1,2,ib)) * invrad + &
                                       ssaice6(1,3,ib)) * invrad + &
                                       ssaice6(1,4,ib)
                      asycoice(ib) = ((asyice6(1,1,ib) * invrad + &
                                       asyice6(1,2,ib)) * invrad + &
                                       asyice6(1,3,ib)) * invrad + &
                                       asyice6(1,4,ib)

                   else
                      extcoice(ib) = ((((((extice6(2,1,ib) * invrad + &
                                           extice6(2,2,ib)) * invrad + &
                                           extice6(2,3,ib)) * invrad + &
                                           extice6(2,4,ib)) * invrad + &
                                           extice6(2,5,ib)) * invrad + &
                                           extice6(2,6,ib)) * invrad + &
                                           extice6(2,7,ib)) * invrad + &
                                           extice6(2,8,ib)
                      ssacoice(ib) = ((((((ssaice6(2,1,ib) * invrad + &
                                           ssaice6(2,2,ib)) * invrad + &
                                           ssaice6(2,3,ib)) * invrad + &
                                           ssaice6(2,4,ib)) * invrad + &
                                           ssaice6(2,5,ib)) * invrad + &
                                           ssaice6(2,6,ib)) * invrad + &
                                           ssaice6(2,7,ib)) * invrad + &
                                           ssaice6(2,8,ib)
                      asycoice(ib) = ((((((asyice6(2,1,ib) * invrad + &
                                           asyice6(2,2,ib)) * invrad + &
                                           asyice6(2,3,ib)) * invrad + &
                                           asyice6(2,4,ib)) * invrad + &
                                           asyice6(2,5,ib)) * invrad + &
                                           asyice6(2,6,ib)) * invrad + &
                                           asyice6(2,7,ib)) * invrad + &
                                           asyice6(2,8,ib)
                   endif  ! end if of diaice .ge. 25.
                 enddo    ! end do of lwbands for cloud radiative coefficients
             else 
                 call endrun('ERROR in UMRad: current iceflag is not available.')
             end if  ! iceflag=1
             iceind = 2

             do ib=1,nlwbands
                taucloud = iciwp(i,k) * 1000._r8 * extcoice(icb(ib,iceind))  ! change iwp from kg/m2 to g/m2 then compute ice cloud optical thickness
                ssacloud = ssacoice(icb(ib,iceind))
                ! delta-scaling (Liou, 2002, 313p)
                asycloud = asycoice(icb(ib,iceind))
                fpeak = asycloud * asycloud

                abs_tamu  (ib,i,k)    = taucloud * (1._r8-ssacloud)  ! absorption optical depth

                taucloud = (1._r8-ssacloud*fpeak) * taucloud   ! delta-scaling technique, ref: Joseph, Wiscombe and Weinman (1976, JAS)
                ssacloud = (1._r8-fpeak)*ssacloud / &
                                   (1._r8-ssacloud*fpeak)
                asycloud = (asycloud-fpeak) / (1.0_r8-fpeak)
                xmomcloud(0) = 1._r8
                xmomcloud(1) = asycloud

                ext_tamu  (ib,i,k)    = taucloud       ! cloud extinction optical depth, i.e. including absorption and scattering
                ssa_tamu  (ib,i,k)    = ssacloud       ! single scattering albedo
                xmomc_tamu(0:,ib,i,k) = xmomcloud(0:)  ! asymmetric factor
             enddo    ! end do of lwbands for cloud radiative coefficients
  
          ! if ice effective diameter is less than 3 microns
          else
            call endrun('ice effective diameter is less than lower boundl 3 microns) in TAMU cloud ice optics')
  
          endif ! end if of ice effective diameter range
        endif   ! end if of ice cloud layer
      enddo     ! end do of k
   enddo        ! end do of i


!---------
! return
!---------
   ext_od  (:,:,:) = ext_tamu  (:,:,:) 
   abs_od  (:,:,:) = abs_tamu  (:,:,:) ! return abs_tamu to return variable, abs_od
   ssa_od  (:,:,:) = ssa_tamu  (:,:,:)
   xmomc_od(:,:,:,:) = xmomc_tamu(:,:,:,:)

  return
end subroutine mc6_ice_get_rad_props_lw

!==============================================================================

end module cloud_rad_props
