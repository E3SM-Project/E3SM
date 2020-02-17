module slingo

   !------------------------------------------------------------------------------------------------
   !  Implements Slingo Optics for MG/RRTMG for liquid clouds and
   !  a copy of the old cloud routine for reference
   !------------------------------------------------------------------------------------------------

   use shr_kind_mod,     only: r8 => shr_kind_r8
   use radconstants,     only: nswbands, nlwbands, get_sw_spectral_boundaries

   implicit none
   private
   save

   public :: &
      slingo_rad_props_init,        &
      slingo_liq_optics_lw, &
      slingo_liq_optics_sw

contains

   !===========================================================================

   subroutine slingo_rad_props_init()

      return

   end subroutine slingo_rad_props_init

   !===========================================================================

   subroutine slingo_liq_optics_sw(ncol, nlev, cldn, cliqwp, rel, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f)

      integer, intent(in) :: ncol, nlev

      ! Inputs have dimension ncol,nlev
      real(r8), intent(in), dimension(:,:) :: rel
      real(r8), intent(in), dimension(:,:) :: cldn
      real(r8), intent(in), dimension(:,:) :: cliqwp 

      ! Outputs have dimension nbnd,ncol,nlev
      real(r8),intent(out) :: liq_tau    (:,:,:) ! extinction optical depth
      real(r8),intent(out) :: liq_tau_w  (:,:,:) ! single scattering albedo * tau
      real(r8),intent(out) :: liq_tau_w_g(:,:,:) ! assymetry parameter * tau * w
      real(r8),intent(out) :: liq_tau_w_f(:,:,:) ! forward scattered fraction * tau * w

      real(r8), dimension(nswbands)     :: wavmin
      real(r8), dimension(nswbands)     :: wavmax

      ! Minimum cloud amount (as a fraction of the grid-box area) to 
      ! distinguish from clear sky
      real(r8), parameter :: cldmin = 1.0e-80_r8

      ! Decimal precision of cloud amount (0 -> preserve full resolution;
      ! 10^-n -> preserve n digits of cloud amount)
      real(r8), parameter :: cldeps = 0.0_r8

      ! A. Slingo's data for cloud particle radiative properties (from 'A GCM
      ! Parameterization for the Shortwave Properties of Water Clouds' JAS
      ! vol. 46 may 1989 pp 1419-1427)
      real(r8) :: abarl(4) = &  ! A coefficient for extinction optical depth
         (/ 2.817e-02_r8, 2.682e-02_r8,2.264e-02_r8,1.281e-02_r8/)
      real(r8) :: bbarl(4) = &  ! B coefficient for extinction optical depth
         (/ 1.305_r8    , 1.346_r8    ,1.454_r8    ,1.641_r8    /)
      real(r8) :: cbarl(4) = &  ! C coefficient for single scat albedo
         (/-5.62e-08_r8 ,-6.94e-06_r8 ,4.64e-04_r8 ,0.201_r8    /)
      real(r8) :: dbarl(4) = &  ! D coefficient for single  scat albedo
         (/ 1.63e-07_r8 , 2.35e-05_r8 ,1.24e-03_r8 ,7.56e-03_r8 /)
      real(r8) :: ebarl(4) = &  ! E coefficient for asymmetry parameter
         (/ 0.829_r8    , 0.794_r8    ,0.754_r8    ,0.826_r8    /)
      real(r8) :: fbarl(4) = &  ! F coefficient for asymmetry parameter
         (/ 2.482e-03_r8, 4.226e-03_r8,6.560e-03_r8,4.353e-03_r8/)

      real(r8) :: abarli        ! A coefficient for current spectral band
      real(r8) :: bbarli        ! B coefficient for current spectral band
      real(r8) :: cbarli        ! C coefficient for current spectral band
      real(r8) :: dbarli        ! D coefficient for current spectral band
      real(r8) :: ebarli        ! E coefficient for current spectral band
      real(r8) :: fbarli        ! F coefficient for current spectral band

      ! Caution... A. Slingo recommends no less than 4.0 micro-meters nor
      ! greater than 20 micro-meters. Here we set effective radius limits
      ! for liquid to the range 4.2 < rel < 16 micron (Slingo 89)
      real(r8), parameter :: rel_min = 4.2_r8
      real(r8), parameter :: rel_max = 16._r8

      integer :: ns, i, k, indxsl
      real(r8) :: tmp1l, tmp2l, tmp3l, g

      call get_sw_spectral_boundaries(wavmin,wavmax,'microns')
      do ns = 1, nswbands
         ! Set index for cloud particle properties based on the wavelength,
         ! according to A. Slingo (1989) equations 1-3:
         ! Use index 1 (0.25 to 0.69 micrometers) for visible
         ! Use index 2 (0.69 - 1.19 micrometers) for near-infrared
         ! Use index 3 (1.19 to 2.38 micrometers) for near-infrared
         ! Use index 4 (2.38 to 4.00 micrometers) for near-infrared
         if(wavmax(ns) <= 0.7_r8) then
            indxsl = 1
         else if(wavmax(ns) <= 1.25_r8) then
            indxsl = 2
         else if(wavmax(ns) <= 2.38_r8) then
            indxsl = 3
         else if(wavmax(ns) > 2.38_r8) then
            indxsl = 4
         end if

         ! Set cloud extinction optical depth, single scatter albedo,
         ! asymmetry parameter, and forward scattered fraction:
         abarli = abarl(indxsl)
         bbarli = bbarl(indxsl)
         cbarli = cbarl(indxsl)
         dbarli = dbarl(indxsl)
         ebarli = ebarl(indxsl)
         fbarli = fbarl(indxsl)

         do k=1,nlev
            do i=1,ncol

               ! note that optical properties for liquid valid only
               ! in range of 4.2 > rel > 16 micron (Slingo 89)
               if (cldn(i,k) >= cldmin .and. cldn(i,k) >= cldeps) then
                  tmp1l = abarli + bbarli/min(max(rel_min,rel(i,k)),rel_max)
                  liq_tau(ns,i,k) = 1000._r8*cliqwp(i,k)*tmp1l
               else
                  liq_tau(ns,i,k) = 0.0_r8
               endif

               tmp2l = 1._r8 - cbarli - dbarli*min(max(rel_min,rel(i,k)),rel_max)
               tmp3l = fbarli*min(max(rel_min,rel(i,k)),rel_max)
               ! Do not let single scatter albedo be 1.  Delta-eddington solution
               ! for non-conservative case has different analytic form from solution
               ! for conservative case, and raddedmx is written for non-conservative case.
               liq_tau_w(ns,i,k) = liq_tau(ns,i,k) * min(tmp2l,.999999_r8)
               g = ebarli + tmp3l
               liq_tau_w_g(ns,i,k) = liq_tau_w(ns,i,k) * g
               liq_tau_w_f(ns,i,k) = liq_tau_w(ns,i,k) * g * g

            end do ! End do i=1,ncol
         end do ! End do k=1,nlev
      end do ! nswbands

   end subroutine slingo_liq_optics_sw

   !===========================================================================

   subroutine slingo_liq_optics_lw(ncol, nlev, cldn, iclwpth, iciwpth, abs_od)

      integer, intent(in) :: ncol, nlev
      real(r8), intent(in), dimension(:,:) :: cldn, iclwpth, iciwpth
      real(r8), intent(out) :: abs_od(:,:,:)

      real(r8) :: ficemr(ncol,nlev)
      real(r8) :: cwp(ncol,nlev)
      real(r8) :: cldtau(ncol,nlev)

      integer :: lwband, i, k

      real(r8) :: kabs, kabsi

      ! longwave liquid absorption coeff (m**2/g)
      real(r8), parameter :: kabsl = 0.090361_r8

      do k=1,nlev
         do i = 1,ncol
            cwp   (i,k) = 1000.0_r8 * iclwpth(i,k) + 1000.0_r8 * iciwpth(i, k)
            ficemr(i,k) = 1000.0_r8 * iciwpth(i,k)/(max(1.e-18_r8, cwp(i,k)))
         end do
      end do

      do k=1,nlev
         do i=1,ncol
            ! Note from Andrew Conley:
            !  Optics for RK no longer supported, This is constructed to get
            !  close to bit for bit.  Otherwise we could simply use liquid water path
            !note that optical properties for ice valid only
            !in range of 13 > rei > 130 micron (Ebert and Curry 92)
            kabs = kabsl*(1._r8-ficemr(i,k))
            cldtau(i,k) = kabs*cwp(i,k)
         end do
      end do

      do lwband = 1,nlwbands
         abs_od(lwband,1:ncol,1:nlev)=cldtau(1:ncol,1:nlev)
      enddo

   end subroutine slingo_liq_optics_lw

   !===========================================================================

end module slingo
