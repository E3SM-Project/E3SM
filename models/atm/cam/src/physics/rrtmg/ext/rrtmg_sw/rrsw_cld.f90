      module rrsw_cld

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_sw cloud property coefficients
!
! Initial: J.-J. Morcrette, ECMWF, oct1999
! Revised: J. Delamere/MJIacono, AER, aug2005
! Revised: MJIacono, AER, nov2005
! Revised: MJIacono, AER, jul2006
!------------------------------------------------------------------
!
!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! xxxliq1 : real   : optical properties (extinction coefficient, single 
!                    scattering albedo, assymetry factor) from
!                    Hu & Stamnes, j. clim., 6, 728-742, 1993.  
! xxxice2 : real   : optical properties (extinction coefficient, single 
!                    scattering albedo, assymetry factor) from streamer v3.0,
!                    Key, streamer user's guide, cooperative institude 
!                    for meteorological studies, 95 pp., 2001.
! xxxice3 : real   : optical properties (extinction coefficient, single 
!                    scattering albedo, assymetry factor) from
!                    Fu, j. clim., 9, 1996.
! xbari   : real   : optical property coefficients for five spectral 
!                    intervals (2857-4000, 4000-5263, 5263-7692, 7692-14285,
!                    and 14285-40000 wavenumbers) following 
!                    Ebert and Curry, jgr, 97, 3831-3836, 1992.
!------------------------------------------------------------------

      real(kind=r8) :: extliq1(58,16:29), ssaliq1(58,16:29), asyliq1(58,16:29)
      real(kind=r8) :: extice2(43,16:29), ssaice2(43,16:29), asyice2(43,16:29)
      real(kind=r8) :: extice3(46,16:29), ssaice3(46,16:29), asyice3(46,16:29)
      real(kind=r8) :: fdlice3(46,16:29)
      real(kind=r8) :: abari(5),bbari(5),cbari(5),dbari(5),ebari(5),fbari(5)

      end module rrsw_cld

