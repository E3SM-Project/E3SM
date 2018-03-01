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

!==============================================================================

end module cloud_rad_props
