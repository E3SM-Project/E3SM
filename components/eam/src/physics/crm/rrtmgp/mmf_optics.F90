module mmf_optics

   use shr_kind_mod, only: r8=>shr_kind_r8
   use assertions, only: assert, assert_valid, assert_range
   use radconstants, only: nswbands, nlwbands
   use cam_abortutils, only: endrun

   implicit none
   private

   public get_cloud_optics_sw, &
          get_cloud_optics_lw, &
          sample_cloud_optics_sw, &
          sample_cloud_optics_lw

contains

   !-------------------------------------------------------------------------------

   subroutine get_cloud_optics_sw(ncol, nlev, nbnd, cld, iclwp, iciwp, rel, rei, tau_out, ssa_out, asm_out)

      use cloud_rad_props, only: mitchell_ice_optics_sw, gammadist_liq_optics_sw
      use ebert_curry, only: ec_ice_optics_sw
      use slingo, only: slingo_liq_optics_sw

      integer, intent(in) :: ncol, nlev, nbnd
      real(r8), intent(in), dimension(:,:) :: cld, iclwp, iciwp, rel, rei

      ! Outputs are shortwave cloud optical properties *by band*. Dimensions should
      ! be nswbands,ncol,nlev. Generally, this should be able to handle cases were
      ! ncol might be something like nday, and nlev could be arbitrary so long as
      ! corresponding fields were defined for all indices of nlev. This
      ! isn't the case right now I don't think, as cloud_rad_props makes explicit
      ! assumptions about array sizes.
      real(r8), intent(out), dimension(:,:,:) :: tau_out, ssa_out, asm_out

      ! Temporary variables to hold cloud optical properties before combining into output arrays.
      real(r8), dimension(nbnd,ncol,nlev) :: &
            liq_tau, liq_tau_ssa, liq_tau_ssa_g, liq_tau_ssa_f, &
            ice_tau, ice_tau_ssa, ice_tau_ssa_g, ice_tau_ssa_f, &
            cld_tau, cld_tau_ssa, cld_tau_ssa_g, cld_tau_ssa_f

      integer :: iband, ilev, icol

      ! Initialize outputs
      tau_out = 0._r8
      ssa_out = 0._r8
      asm_out = 0._r8

      ! Initialize local variables
      ice_tau = 0._r8
      ice_tau_ssa = 0._r8
      ice_tau_ssa_g = 0._r8
      ice_tau_ssa_f = 0._r8
      liq_tau = 0._r8
      liq_tau_ssa = 0._r8
      liq_tau_ssa_g = 0._r8
      liq_tau_ssa_f = 0._r8

      ! Get ice cloud optics
      call ec_ice_optics_sw(ncol, nlev, cld, iciwp, rei, ice_tau, ice_tau_ssa, ice_tau_ssa_g, ice_tau_ssa_f)
      
      ! Get liquid cloud optics
      call slingo_liq_optics_sw(ncol, nlev, cld, iclwp, rel, liq_tau, liq_tau_ssa, liq_tau_ssa_g, liq_tau_ssa_f)

      ! Combine all cloud optics from CAM routines
      cld_tau = ice_tau + liq_tau
      cld_tau_ssa = ice_tau_ssa + liq_tau_ssa
      cld_tau_ssa_g = ice_tau_ssa_g + liq_tau_ssa_g
      
      ! Copy to output arrays, converting to optical depth, single scattering
      ! albedo, and assymmetry parameter from the products that the CAM routines
      ! return. Make sure we do not try to divide by zero...
      do iband = 1,nbnd
         tau_out(:ncol,:nlev,iband) = cld_tau(iband,:ncol,:nlev)
         where (cld_tau(iband,:ncol,:nlev) > 0)
            ssa_out(:ncol,:nlev,iband) &
               = cld_tau_ssa(iband,:ncol,:nlev) / cld_tau(iband,:ncol,:nlev)
         elsewhere
            ssa_out(:ncol,:nlev,iband) = 1._r8
         endwhere
         where (cld_tau_ssa(iband,:ncol,:nlev) > 0)
            asm_out(:ncol,:nlev,iband) &
               = cld_tau_ssa_g(iband,:ncol,:nlev) / cld_tau_ssa(iband,:ncol,:nlev)
         elsewhere
            asm_out(:ncol,:nlev,iband) = 0._r8
         end where
      end do
   end subroutine get_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine get_cloud_optics_lw(ncol, nlev, nbnd, cld, iclwp, iciwp, rei, tau_out)

      use ebert_curry, only: ec_ice_optics_lw
      use slingo, only: slingo_liq_optics_lw
      use radconstants, only: nlwbands

      integer, intent(in) :: ncol, nlev, nbnd
      real(r8), intent(in), dimension(:,:) :: cld, iclwp, iciwp, rei
      real(r8), intent(out), dimension(:,:,:) :: tau_out

      ! Temporary variables to hold absorption optical depth
      real(r8), dimension(nbnd,ncol,nlev) :: ice_tau, liq_tau, cld_tau

      integer :: iband

      ! Initialize outputs
      tau_out = 0._r8

      ! Initialize local variables
      ice_tau(:,:,:) = 0._r8
      liq_tau(:,:,:) = 0._r8
      cld_tau(:,:,:) = 0._r8

      ! Get ice optics
      call ec_ice_optics_lw(ncol, nlev, cld, iclwp, iciwp, rei, ice_tau)

      ! Get liquid optics
      call slingo_liq_optics_lw(ncol, nlev, cld, iclwp, iciwp, liq_tau)

      ! Combined cloud optics
      cld_tau = liq_tau + ice_tau

      ! Set output optics
      do iband = 1,nbnd
         tau_out(1:ncol,1:nlev,iband) = cld_tau(iband,1:ncol,1:nlev)
      end do

   end subroutine get_cloud_optics_lw

   !----------------------------------------------------------------------------

   ! Do MCICA sampling of optics here. This will map bands to gpoints,
   ! while doing stochastic sampling of cloud state
   subroutine sample_cloud_optics_sw( &
         ncol, nlev, ngpt, &
         gpt2bnd, pmid, cld, &
         tau_bnd, ssa_bnd, asm_bnd, &
         tau_gpt, ssa_gpt, asm_gpt)
      use mcica_subcol_gen, only: mcica_subcol_mask
      integer, intent(in) :: ncol, nlev, ngpt
      integer, intent(in), dimension(:) :: gpt2bnd
      real(r8), intent(in), dimension(:,:) :: pmid, cld
      real(r8), intent(in), dimension(:,:,:) :: tau_bnd, ssa_bnd, asm_bnd
      real(r8), intent(out), dimension(:,:,:) :: tau_gpt, ssa_gpt, asm_gpt
      logical, dimension(ngpt,ncol,nlev) :: iscloudy

      ! For MCICA sampling routine, how many times to permute random seed
      integer, parameter :: changeseed = 1

      ! Loop variables
      integer :: icol, ilev, igpt

      ! Get stochastic subcolumn cloud mask
      call mcica_subcol_mask(ngpt, ncol, nlev, changeseed, &
                             pmid(1:ncol,1:nlev), &
                             cld(1:ncol,1:nlev), &
                             iscloudy(1:ngpt,1:ncol,1:nlev))
      ! Generate subcolumns for homogeneous clouds
      do igpt = 1,ngpt
         do ilev = 1,nlev
            do icol = 1,ncol
               if (iscloudy(igpt,icol,ilev) .and. cld(icol,ilev) > 0._r8) then
                  tau_gpt(icol,ilev,igpt) = tau_bnd(icol,ilev,gpt2bnd(igpt))
                  ssa_gpt(icol,ilev,igpt) = ssa_bnd(icol,ilev,gpt2bnd(igpt))
                  asm_gpt(icol,ilev,igpt) = asm_bnd(icol,ilev,gpt2bnd(igpt))
               else
                  tau_gpt(icol,ilev,igpt) = 0._r8
                  ssa_gpt(icol,ilev,igpt) = 1._r8
                  asm_gpt(icol,ilev,igpt) = 0._r8
               end if
            end do
         end do
      end do
   end subroutine sample_cloud_optics_sw

   !----------------------------------------------------------------------------

   ! Do MCICA sampling of optics here. This will map bands to gpoints,
   ! while doing stochastic sampling of cloud state
   subroutine sample_cloud_optics_lw(&
         ncol, nlev, ngpt, &
         gpt2bnd, pmid, cld, &
         tau_bnd, tau_gpt)
      use mcica_subcol_gen, only: mcica_subcol_mask

      integer, intent(in) :: ncol, nlev, ngpt
      integer, intent(in), dimension(:) :: gpt2bnd
      real(r8), intent(in), dimension(:,:) :: pmid, cld
      real(r8), intent(in), dimension(:,:,:) :: tau_bnd
      real(r8), intent(out), dimension(:,:,:) :: tau_gpt
      logical, dimension(ngpt,ncol,nlev) :: iscloudy

      ! For MCICA sampling routine, how many times to permute random seed
      integer, parameter :: changeseed = 1

      ! Loop variables
      integer :: icol, ilev, igpt

      ! Get the stochastic subcolumn cloudy mask
      call mcica_subcol_mask(ngpt, ncol, nlev, changeseed, &
                             pmid(1:ncol,1:nlev), &
                             cld(1:ncol,1:nlev), &
                             iscloudy(1:ngpt,1:ncol,1:nlev))

      ! Map optics to g-points, selecting a single subcolumn for each
      ! g-point. This implementation generates homogeneous clouds, but it would be
      ! straightforward to extend this to handle horizontally heterogeneous clouds
      ! as well.
      do igpt = 1,ngpt
         do ilev = 1,nlev
            do icol = 1,ncol
               if (iscloudy(igpt,icol,ilev) .and. cld(icol,ilev) > 0._r8) then
                  tau_gpt(icol,ilev,igpt) = tau_bnd(icol,ilev,gpt2bnd(igpt))
               else
                  tau_gpt(icol,ilev,igpt) = 0._r8
               end if
            end do
         end do
      end do
   end subroutine sample_cloud_optics_lw

   !----------------------------------------------------------------------------

end module mmf_optics
