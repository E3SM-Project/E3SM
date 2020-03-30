module cloud_rad_props

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use radconstants,     only: nswbands, nlwbands
use rad_constituents, only: iceopticsfile, liqopticsfile
use interpolate_data, only: interp_type, lininterp_init, lininterp, &
                            extrap_method_bndry, lininterp_finish
use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun

implicit none
private
save

public :: &
   cloud_rad_props_init,    &
   mitchell_ice_optics_sw,  & ! return Mitchell SW ice radiative properties
   mitchell_ice_optics_lw,  & ! Mitchell LW ice rad props
   gammadist_liq_optics_sw, & ! return Conley SW rad props
   gammadist_liq_optics_lw    ! return Conley LW rad props

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

   if (.not.allocated(g_mu)) allocate(g_mu(nmu))
   if (.not.allocated(g_lambda)) allocate(g_lambda(nmu,nlambda))
   if (.not.allocated(ext_sw_liq)) allocate(ext_sw_liq(nmu,nlambda,nswbands) )
   if (.not.allocated(ssa_sw_liq)) allocate(ssa_sw_liq(nmu,nlambda,nswbands))
   if (.not.allocated(asm_sw_liq)) allocate(asm_sw_liq(nmu,nlambda,nswbands))
   if (.not.allocated(abs_lw_liq)) allocate(abs_lw_liq(nmu,nlambda,nlwbands))

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

   if (.not.allocated(g_d_eff)) allocate(g_d_eff(n_g_d))
   if (.not.allocated(ext_sw_ice)) allocate(ext_sw_ice(n_g_d,nswbands))
   if (.not.allocated(ssa_sw_ice)) allocate(ssa_sw_ice(n_g_d,nswbands))
   if (.not.allocated(asm_sw_ice)) allocate(asm_sw_ice(n_g_d,nswbands))
   if (.not.allocated(abs_lw_ice)) allocate(abs_lw_ice(n_g_d,nlwbands))

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

subroutine gammadist_liq_optics_sw(ncol, nlev, iclwpth, lamc, pgam, tau, tau_w, tau_w_g, tau_w_f)
   integer, intent(in) :: ncol, nlev

   ! Inputs have dimension ncol,nlev
   real(r8), intent(in), dimension(:,:) :: lamc, pgam, iclwpth

   ! Outputs have dimension nbnd,ncol,nlev
   real(r8),intent(out) :: tau    (:,:,:) ! extinction optical depth
   real(r8),intent(out) :: tau_w  (:,:,:) ! single scattering albedo * tau
   real(r8),intent(out) :: tau_w_g(:,:,:) ! asymetry parameter * tau * w
   real(r8),intent(out) :: tau_w_f(:,:,:) ! forward scattered fraction * tau * w
   integer i,k,swband

   do k = 1,nlev
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

end subroutine gammadist_liq_optics_sw

!==============================================================================

subroutine gammadist_liq_optics_lw(ncol, nlev, iclwpth, lamc, pgam, abs_od)
   integer, intent(in) :: ncol, nlev
   real(r8), intent(in), dimension(:,:) :: lamc, pgam, iclwpth
   real(r8), intent(out) :: abs_od(:,:,:)
   integer lwband, i, k

   abs_od = 0._r8
   do k = 1,nlev
      do i = 1,ncol
         if(lamc(i,k) > 0._r8) then ! This seems to be the clue for no cloud from microphysics formulation
            call gam_liquid_lw(iclwpth(i,k), lamc(i,k), pgam(i,k), abs_od(1:nlwbands,i,k))
         else
            abs_od(1:nlwbands,i,k) = 0._r8
         endif
      enddo
   enddo
end subroutine gammadist_liq_optics_lw

!==============================================================================

subroutine mitchell_ice_optics_sw(ncol, nlev, iciwpth, dei, tau, tau_w, &
     tau_w_g, tau_w_f)

  ! Inputs have dimension ncol,nlev
  integer, intent(in) :: ncol, nlev
  real(r8), intent(in) :: iciwpth(:,:)
  real(r8), intent(in) :: dei(:,:)

  ! Outputs have dimension nswbands,ncol,nlev
  real(r8),intent(out) :: tau    (:,:,:) ! extinction optical depth
  real(r8),intent(out) :: tau_w  (:,:,:) ! single scattering albedo * tau
  real(r8),intent(out) :: tau_w_g(:,:,:) ! assymetry parameter * tau * w
  real(r8),intent(out) :: tau_w_f(:,:,:) ! forward scattered fraction * tau * w

  type(interp_type) :: dei_wgts

  integer :: i, k, swband
  real(r8) :: ext(nswbands), ssa(nswbands), asm(nswbands)

  do k = 1,nlev
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

end subroutine mitchell_ice_optics_sw

!==============================================================================

subroutine mitchell_ice_optics_lw(ncol, nlev, iciwpth, dei, abs_od)

  integer, intent(in) :: ncol, nlev
  real(r8), intent(in) :: iciwpth(:,:)
  real(r8), intent(in) :: dei(:,:)

  real(r8),intent(out) :: abs_od(:,:,:)

  type(interp_type) :: dei_wgts

  integer :: i, k, lwband
  real(r8) :: absor(nlwbands)

  do k = 1,nlev
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
           call lininterp_finish(dei_wgts)
        endif
     enddo
  enddo

end subroutine mitchell_ice_optics_lw

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
