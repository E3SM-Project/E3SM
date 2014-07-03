
      module mo_setair

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setair

      contains

      subroutine setair( z, nw, wc, airlev, dtrl, &
                         cz, o2top )
!-----------------------------------------------------------------------------
!   purpose:
!   set up an altitude profile of air molecules.  subroutine includes a
!   shape-conserving scaling method that allows scaling of the entire
!   profile to a given sea-level pressure.
!-----------------------------------------------------------------------------
!   parameters:
!   pmbnew  - real(r8), sea-level pressure (mb) to which profile should be
!             scaled.  if pmbnew < 0, no scaling is done
!   nz      - integer, number of specified altitude levels in the working (i)
!             grid
!   z       - real(r8), specified altitude working grid (km)                  (i)
!   nw      - integer, number of specified intervals + 1 in working       (i)
!             wavelength grid
!   wl      - real(r8), vector of lower limits of wavelength intervals in     (i)
!             working wavelength grid
!   airlev  - real(r8), air density (molec/cc) at each specified altitude     (o)
!   dtrl    - real(r8), rayleigh optical depth at each specified altitude     (o)
!             and each specified wavelength
!   cz      - real(r8), number of air molecules per cm^2 at each specified    (o)
!             altitude layer
!-----------------------------------------------------------------------------

      use mo_params, only : kw
      use ppgrid,    only : pver, pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)  :: nw
      real(r8), intent(in) :: o2top
      real(r8), intent(in) :: wc(kw)
      real(r8), intent(in) :: z(pverp)
      real(r8), intent(in) :: airlev(pverp)
!-----------------------------------------------------------------------------
! 	... air density (molec cm-3) at each grid level
!           rayleigh optical depths
!-----------------------------------------------------------------------------
      real(r8), intent(out) :: dtrl(pver,nw)
      real(r8), intent(out) :: cz(pverp)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer  :: k, kp1, wn
      real(r8) :: hscale
      real(r8) :: srayl
      real(r8) :: deltaz
      real(r8) :: wmicrn, xx 

!-----------------------------------------------------
! 	... compute column increments (logarithmic integrals)
!-----------------------------------------------------
      do k = 1,pver
         kp1 = k + 1
         deltaz = 1.e5_r8 * (z(kp1) - z(k)) 
         cz(k)  =  (airlev(kp1) - airlev(k)) / log( airlev(kp1)/airlev(k) ) * deltaz
      end do

!-----------------------------------------------------
! 	... include exponential tail integral from infinity to 50 km,
!           fold tail integral into top layer
!           specify scale height near top of data.  (scale height at 40km????)
!-----------------------------------------------------
!     hscale = 8.05e5
      cz(pverp) = o2top/.2095_r8

!-----------------------------------------------------
! 	... compute rayleigh cross sections and depths:
!-----------------------------------------------------
      do wn = 1,nw
!-----------------------------------------------------
! 	... rayleigh scattering cross section from wmo 1985 (originally from
!           nicolet, m., on the molecular scattering in the terrestrial atmosphere:
!           an empirical formula for its calculation in the homoshpere, planet.
!           space sci., 32, 1467-1468, 1984.
!-----------------------------------------------------
         wmicrn =  1.e-3_r8*wc(wn)
         if( wmicrn <= .55_r8 ) then
            xx = 3.6772_r8 + 0.389_r8*wmicrn + 0.09426_r8/wmicrn
         else
            xx = 4.04_r8
         end if
         srayl             = 4.02e-28_r8/(wmicrn)**xx
         dtrl(1:pver-1,wn) = cz(1:pver-1)*srayl
         dtrl(pver,wn)     = (cz(pver) + cz(pverp))*srayl
      end do

      end subroutine setair

      end module mo_setair
