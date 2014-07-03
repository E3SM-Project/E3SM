module cloud_rad_props

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use ppgrid,           only: pcols, pver, pverp
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field, pbuf_old_tim_idx
use radconstants,     only: nswbands, nlwbands, idx_sw_diag, ot_length, idx_lw_diag
use abortutils,       only: endrun
use rad_constituents, only: iceopticsfile, liqopticsfile
use oldcloud,         only: oldcloud_lw, old_liq_get_rad_props_lw, old_ice_get_rad_props_lw, oldcloud_init

use ebert_curry,      only: scalefactor
use cam_logfile,      only: iulog

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

real(r8) :: dmin = 2. *  13. / scalefactor
real(r8) :: dmax = 2. * 130. / scalefactor

! Minimum cloud amount (as a fraction of the grid-box area) to 
! distinguish from clear sky
! 
   real(r8) cldmin
   parameter (cldmin = 1.0e-80_r8)
!
! Decimal precision of cloud amount (0 -> preserve full resolution;
! 10^-n -> preserve n digits of cloud amount)
! 
   real(r8) cldeps
   parameter (cldeps = 0.0_r8)

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

subroutine get_snow_optics_sw   (state, pbuf, tau, tau_w, tau_w_g, tau_w_f)
   type(physics_state), intent(in)   :: state
   type(physics_buffer_desc),pointer :: pbuf(:)

   real(r8),intent(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
   real(r8),intent(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
   real(r8),intent(out) :: tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
   real(r8),intent(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

   real(r8), pointer, dimension(:,:) :: dei, iciwpth
   real(r8) :: dlimited ! d limited to dmin,dmax range

   integer :: i,k,i_d_grid,k_d_eff,i_swband,ncol,lchnk
   real :: wd, onemwd, ext, ssa, asm

   ncol = state%ncol
   lchnk = state%lchnk

! temporary code to support diagnostics of snow radiation
   call pbuf_get_field(pbuf, i_icswp, iciwpth)
   call pbuf_get_field(pbuf, i_des,   dei)
! temporary code to support diagnostics of snow radiation

   do i = 1,ncol
      do k = 1,pver
         dlimited = dei(i,k) ! min(dmax,max(dei(i,k),dmin))
         if( iciwpth(i,k) < 1.e-80_r8 .or. dlimited .eq. 0._r8) then
         ! if ice water path is too small, OD := 0
            tau    (:,i,k) = 0._r8
            tau_w  (:,i,k) = 0._r8
            tau_w_g(:,i,k) = 0._r8
            tau_w_f(:,i,k) = 0._r8
         else 
            !if (dlimited < g_d_eff(1) .or. dlimited > g_d_eff(n_g_d)) then
               !write(iulog,*) 'dei from prognostic cldwat2m',dei(i,k)
               !write(iulog,*) 'grid values of deff ice from optics file',g_d_eff
               !call endrun ('deff of ice exceeds limits')
            !endif
            ! for each cell interpolate to find weights and indices in g_d_eff grid.
            if (dlimited <= g_d_eff(1)) then
               k_d_eff = 2
               wd = 1._r8
               onemwd = 0._r8
            elseif (dlimited >= g_d_eff(n_g_d)) then
               k_d_eff = n_g_d
               wd = 0._r8
               onemwd = 1._r8 
            else
               do i_d_grid = 1, n_g_d
                  k_d_eff = i_d_grid
                  if(g_d_eff(i_d_grid) > dlimited) exit
               enddo
               wd = (g_d_eff(k_d_eff) - dlimited)/(g_d_eff(k_d_eff) - g_d_eff(k_d_eff-1))
               onemwd = 1._r8 - wd
            endif
            ! interpolate into grid and extract radiative properties
            do i_swband = 1, nswbands
               ext = wd*ext_sw_ice(k_d_eff-1,i_swband) + &
                 onemwd*ext_sw_ice(k_d_eff  ,i_swband) 
               ssa = wd*ssa_sw_ice(k_d_eff-1,i_swband) + &
                 onemwd*ssa_sw_ice(k_d_eff  ,i_swband) 
               asm = wd*asm_sw_ice(k_d_eff-1,i_swband) + &
                 onemwd*asm_sw_ice(k_d_eff  ,i_swband) 
               tau    (i_swband,i,k)=iciwpth(i,k) * ext
               tau_w  (i_swband,i,k)=iciwpth(i,k) * ext * ssa
               tau_w_g(i_swband,i,k)=iciwpth(i,k) * ext * ssa * asm
               tau_w_f(i_swband,i,k)=iciwpth(i,k) * ext * ssa * asm * asm
            enddo
         endif
      enddo
   enddo

   return
end subroutine get_snow_optics_sw   

!==============================================================================
! Private methods
!==============================================================================

subroutine get_ice_optics_sw   (state, pbuf, tau, tau_w, tau_w_g, tau_w_f)
   type(physics_state), intent(in)   :: state
   type(physics_buffer_desc),pointer :: pbuf(:)

   real(r8),intent(out) :: tau    (nswbands,pcols,pver) ! extinction optical depth
   real(r8),intent(out) :: tau_w  (nswbands,pcols,pver) ! single scattering albedo * tau
   real(r8),intent(out) :: tau_w_g(nswbands,pcols,pver) ! assymetry parameter * tau * w
   real(r8),intent(out) :: tau_w_f(nswbands,pcols,pver) ! forward scattered fraction * tau * w

   real(r8), pointer, dimension(:,:) :: dei, iciwpth
   real(r8) :: dlimited ! d limited to dmin,dmax range

   integer :: i,k,i_d_grid,k_d_eff,i_swband,ncol,lchnk
   real :: wd, onemwd, ext, ssa, asm

   ncol = state%ncol
   lchnk = state%lchnk

   call pbuf_get_field(pbuf, i_iciwp, iciwpth)
   call pbuf_get_field(pbuf, i_dei,   dei)

   do i = 1,ncol
      do k = 1,pver
         dlimited = dei(i,k) ! min(dmax,max(dei(i,k),dmin))
         if( iciwpth(i,k) < 1.e-80_r8 .or. dlimited .eq. 0._r8) then
         ! if ice water path is too small, OD := 0
            tau    (:,i,k) = 0._r8
            tau_w  (:,i,k) = 0._r8
            tau_w_g(:,i,k) = 0._r8
            tau_w_f(:,i,k) = 0._r8
         else 
            !if (dlimited < g_d_eff(1) .or. dlimited > g_d_eff(n_g_d)) then
               !write(iulog,*) 'dei from prognostic cldwat2m',dei(i,k)
               !write(iulog,*) 'grid values of deff ice from optics file',g_d_eff
               !call endrun ('deff of ice exceeds limits')
            !endif
            ! for each cell interpolate to find weights and indices in g_d_eff grid.
            if (dlimited <= g_d_eff(1)) then
               k_d_eff = 2
               wd = 1._r8
               onemwd = 0._r8
            elseif (dlimited >= g_d_eff(n_g_d)) then
               k_d_eff = n_g_d
               wd = 0._r8
               onemwd = 1._r8 
            else
               do i_d_grid = 1, n_g_d
                  k_d_eff = i_d_grid
                  if(g_d_eff(i_d_grid) > dlimited) exit
               enddo
               wd = (g_d_eff(k_d_eff) - dlimited)/(g_d_eff(k_d_eff) - g_d_eff(k_d_eff-1))
               onemwd = 1._r8 - wd
            endif
            ! interpolate into grid and extract radiative properties
            do i_swband = 1, nswbands
               ext = wd*ext_sw_ice(k_d_eff-1,i_swband) + &
                 onemwd*ext_sw_ice(k_d_eff  ,i_swband) 
               ssa = wd*ssa_sw_ice(k_d_eff-1,i_swband) + &
                 onemwd*ssa_sw_ice(k_d_eff  ,i_swband) 
               asm = wd*asm_sw_ice(k_d_eff-1,i_swband) + &
                 onemwd*asm_sw_ice(k_d_eff  ,i_swband) 
               tau    (i_swband,i,k)=iciwpth(i,k) * ext
               tau_w  (i_swband,i,k)=iciwpth(i,k) * ext * ssa
               tau_w_g(i_swband,i,k)=iciwpth(i,k) * ext * ssa * asm
               tau_w_f(i_swband,i,k)=iciwpth(i,k) * ext * ssa * asm * asm
            enddo
         endif
      enddo
   enddo

   return
end subroutine get_ice_optics_sw   

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

   real(r8), pointer, dimension(:,:) :: dei, iciwpth
   real(r8) :: dlimited ! d limited to range dmin,dmax

   integer :: i,k,i_d_grid,k_d_eff,i_lwband,ncol,lchnk
   real :: wd, onemwd, absor

   abs_od = 0._r8

   ncol = state%ncol
   lchnk = state%lchnk

! note that this code makes the "ice path" point to the "snow path from CAM"
   
   call pbuf_get_field(pbuf, i_icswp, iciwpth)
   call pbuf_get_field(pbuf, i_des,   dei)

! note that this code makes the "ice path" point to the "snow path from CAM"

   do i = 1,ncol
      do k = 1,pver
         dlimited = dei(i,k) ! min(dmax,max(dei(i,k),dmin))
         ! if ice water path is too small, OD := 0
         if( iciwpth(i,k) < 1.e-80_r8 .or. dlimited .eq. 0._r8) then
            abs_od (:,i,k) = 0._r8
         !else if (dlimited < g_d_eff(1) .or. dlimited > g_d_eff(n_g_d)) then
         !   write(iulog,*) 'dlimited prognostic cldwat2m',dlimited
         !   write(iulog,*) 'grid values of deff ice from optics file',g_d_eff(1),' -> ',g_d_eff(n_g_d)
         !   !call endrun ('deff of ice exceeds limits')
         else
            ! for each cell interpolate to find weights and indices in g_d_eff grid.
            if (dlimited <= g_d_eff(1)) then
               k_d_eff = 2
               wd = 1._r8
               onemwd = 0._r8
            elseif (dlimited >= g_d_eff(n_g_d)) then
               k_d_eff = n_g_d
               wd = 0._r8
               onemwd = 1._r8 
            else
               do i_d_grid = 2, n_g_d
                  k_d_eff = i_d_grid
                  if(g_d_eff(i_d_grid) > dlimited) exit
               enddo
               wd = (g_d_eff(k_d_eff) - dlimited)/(g_d_eff(k_d_eff) - g_d_eff(k_d_eff-1))
               onemwd = 1._r8 - wd
            endif
            ! interpolate into grid and extract radiative properties
            do i_lwband = 1, nlwbands
               absor = wd*abs_lw_ice(k_d_eff-1,i_lwband) + &
                   onemwd*abs_lw_ice(k_d_eff  ,i_lwband)
               abs_od (i_lwband,i,k)=  iciwpth(i,k) * absor 
            enddo
         endif
      enddo
   enddo

end subroutine snow_cloud_get_rad_props_lw

!==============================================================================

subroutine ice_cloud_get_rad_props_lw(state, pbuf, abs_od)
   type(physics_state), intent(in)     :: state
   type(physics_buffer_desc), pointer  :: pbuf(:)
   real(r8), intent(out) :: abs_od(nlwbands,pcols,pver)

   real(r8), pointer, dimension(:,:) :: dei, iciwpth
   real(r8) :: dlimited ! d limited to range dmin,dmax

   integer :: i,k,i_d_grid,k_d_eff,i_lwband,ncol,lchnk
   real :: wd, onemwd, absor

   abs_od = 0._r8

   ncol = state%ncol
   lchnk = state%lchnk

   call pbuf_get_field(pbuf, i_iciwp, iciwpth)
   call pbuf_get_field(pbuf, i_dei,   dei)

   do i = 1,ncol
      do k = 1,pver
         dlimited = dei(i,k) ! min(dmax,max(dei(i,k),dmin))
         ! if ice water path is too small, OD := 0
         if( iciwpth(i,k) < 1.e-80_r8 .or. dlimited .eq. 0._r8) then
            abs_od (:,i,k) = 0._r8
         !else if (dlimited < g_d_eff(1) .or. dlimited > g_d_eff(n_g_d)) then
         !   write(iulog,*) 'dlimited prognostic cldwat2m',dlimited
         !   write(iulog,*) 'grid values of deff ice from optics file',g_d_eff(1),' -> ',g_d_eff(n_g_d)
         !   !call endrun ('deff of ice exceeds limits')
         else
            ! for each cell interpolate to find weights and indices in g_d_eff grid.
            if (dlimited <= g_d_eff(1)) then
               k_d_eff = 2
               wd = 1._r8
               onemwd = 0._r8
            elseif (dlimited >= g_d_eff(n_g_d)) then
               k_d_eff = n_g_d
               wd = 0._r8
               onemwd = 1._r8 
            else
               do i_d_grid = 2, n_g_d
                  k_d_eff = i_d_grid
                  if(g_d_eff(i_d_grid) > dlimited) exit
               enddo
               wd = (g_d_eff(k_d_eff) - dlimited)/(g_d_eff(k_d_eff) - g_d_eff(k_d_eff-1))
               onemwd = 1._r8 - wd
            endif
            ! interpolate into grid and extract radiative properties
            do i_lwband = 1, nlwbands
               absor = wd*abs_lw_ice(k_d_eff-1,i_lwband) + &
                   onemwd*abs_lw_ice(k_d_eff  ,i_lwband)
               abs_od (i_lwband,i,k)=  iciwpth(i,k) * absor 
            enddo
         endif
      enddo
   enddo

end subroutine ice_cloud_get_rad_props_lw

!==============================================================================

subroutine gam_liquid_lw(clwptn, lamc, pgam, abs_od)
  real(r8), intent(in) :: clwptn ! cloud water liquid path new (in cloud) (in g/m^2)?
  real(r8), intent(in) :: lamc   ! prognosed value of lambda for cloud
  real(r8), intent(in) :: pgam   ! prognosed value of mu for cloud
  real(r8), intent(out) :: abs_od(1:nlwbands)
  ! for interpolating into mu/lambda
  integer :: imu, kmu, wmu, onemwmu
  integer :: ilambda, klambda, wlambda, onemwlambda, lambdaplus, lambdaminus
  integer :: lwband ! sw band index
  real(r8) :: absc

  if (clwptn < 1.e-80_r8) then
    abs_od = 0._r8
    return
  endif

  if (pgam < g_mu(1) .or. pgam > g_mu(nmu)) then
    write(iulog,*)'pgam from prognostic cldwat2m',pgam
    write(iulog,*)'g_mu from file',g_mu
    call endrun ('pgam exceeds limits')
  endif
  do imu = 1, nmu
    kmu = imu
    if (g_mu(kmu) > pgam) exit
  enddo
  wmu = (g_mu(kmu) - pgam)/(g_mu(kmu) - g_mu(kmu-1))
  onemwmu = 1._r8 - wmu

  do ilambda = 1, nlambda
    klambda = ilambda
    if (wmu*g_lambda(kmu-1,ilambda) + onemwmu*g_lambda(kmu,ilambda) < lamc) exit
  enddo
  if (klambda <= 1 .or. klambda > nlambda) call endrun('lamc  exceeds limits')
  lambdaplus = wmu*g_lambda(kmu-1,klambda  ) + onemwmu*g_lambda(kmu,klambda  )
  lambdaminus= wmu*g_lambda(kmu-1,klambda-1) + onemwmu*g_lambda(kmu,klambda-1)
  wlambda = (lambdaplus - lamc) / (lambdaplus - lambdaminus)
  onemwlambda = 1._r8 - wlambda

  do lwband = 1, nlwbands
     absc=     wlambda*    wmu*abs_lw_liq(kmu-1,klambda-1,lwband) + &
           onemwlambda*    wmu*abs_lw_liq(kmu-1,klambda  ,lwband) + &
               wlambda*onemwmu*abs_lw_liq(kmu  ,klambda-1,lwband) + &
           onemwlambda*onemwmu*abs_lw_liq(kmu  ,klambda  ,lwband)

     abs_od(lwband) = clwptn * absc
  enddo

  return
end subroutine gam_liquid_lw

!==============================================================================

subroutine gam_liquid_sw(clwptn, lamc, pgam, tau, tau_w, tau_w_g, tau_w_f)
  real(r8), intent(in) :: clwptn ! cloud water liquid path new (in cloud) (in g/m^2)?
  real(r8), intent(in) :: lamc   ! prognosed value of lambda for cloud
  real(r8), intent(in) :: pgam   ! prognosed value of mu for cloud
  real(r8), intent(out) :: tau(1:nswbands), tau_w(1:nswbands), tau_w_f(1:nswbands), tau_w_g(1:nswbands)
  ! for interpolating into mu/lambda
  integer :: imu, kmu, wmu, onemwmu
  integer :: ilambda, klambda, wlambda, onemwlambda, lambdaplus, lambdaminus
  integer :: swband ! sw band index
  real(r8) :: ext, ssa, asm

  if (clwptn < 1.e-80_r8) then
    tau = 0._r8
    tau_w = 0._r8
    tau_w_g = 0._r8
    tau_w_f = 0._r8
    return
  endif

  if (pgam < g_mu(1) .or. pgam > g_mu(nmu)) then
    write(iulog,*)'pgam from prognostic cldwat2m',pgam
    write(iulog,*)'g_mu from file',g_mu
    call endrun ('pgam exceeds limits')
  endif
  do imu = 1, nmu
    kmu = imu
    if (g_mu(kmu) > pgam) exit
  enddo
  wmu = (g_mu(kmu) - pgam)/(g_mu(kmu) - g_mu(kmu-1))
  onemwmu = 1._r8 - wmu

  do ilambda = 1, nlambda
     klambda = ilambda
     if (wmu*g_lambda(kmu-1,ilambda) + onemwmu*g_lambda(kmu,ilambda) < lamc) exit
  enddo
  if (klambda <= 1 .or. klambda > nlambda) call endrun('lamc  exceeds limits')
     lambdaplus = wmu*g_lambda(kmu-1,klambda  ) + onemwmu*g_lambda(kmu,klambda  )
     lambdaminus= wmu*g_lambda(kmu-1,klambda-1) + onemwmu*g_lambda(kmu,klambda-1)
     wlambda = (lambdaplus - lamc) / (lambdaplus - lambdaminus)
  onemwlambda = 1._r8 - wlambda

  do swband = 1, nswbands
     ext =     wlambda*    wmu*ext_sw_liq(kmu-1,klambda-1,swband) + &
           onemwlambda*    wmu*ext_sw_liq(kmu-1,klambda  ,swband) + &
               wlambda*onemwmu*ext_sw_liq(kmu  ,klambda-1,swband) + &
           onemwlambda*onemwmu*ext_sw_liq(kmu  ,klambda  ,swband)
     ! probably should interpolate ext*ssa
     ssa =     wlambda*    wmu*ssa_sw_liq(kmu-1,klambda-1,swband) + &
           onemwlambda*    wmu*ssa_sw_liq(kmu-1,klambda  ,swband) + &
               wlambda*onemwmu*ssa_sw_liq(kmu  ,klambda-1,swband) + &
           onemwlambda*onemwmu*ssa_sw_liq(kmu  ,klambda  ,swband)
     ! probably should interpolate ext*ssa*asm
     asm =     wlambda*    wmu*asm_sw_liq(kmu-1,klambda-1,swband) + &
           onemwlambda*    wmu*asm_sw_liq(kmu-1,klambda  ,swband) + &
               wlambda*onemwmu*asm_sw_liq(kmu  ,klambda-1,swband) + &
           onemwlambda*onemwmu*asm_sw_liq(kmu  ,klambda  ,swband)
     ! compute radiative properties
     tau(swband) = clwptn * ext
     tau_w(swband) = clwptn * ext * ssa
     tau_w_g(swband) = clwptn * ext * ssa * asm
     tau_w_f(swband) = clwptn * ext * ssa * asm * asm
  enddo

  return
end subroutine gam_liquid_sw

!==============================================================================

end module cloud_rad_props
