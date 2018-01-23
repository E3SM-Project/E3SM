      module rrsw_cld

      use parkind, only : im => kind_im, rb => kind_rb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_sw cloud property coefficients
!
! Initial: J.-J. Morcrette, ECMWF, oct1999
! Revised: J. Delamere/MJIacono, AER, aug2005
! Revised: MJIacono, AER, nov2005
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------
!
!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! xxxliq1 : real   : optical properties (extinction coefficient, single 
!                    scattering albedo, assymetry factor) from
!                    Hu & Stamnes, j. clim., 6, 728-742, 1993.  
! xxxice2 : real   : optical properties (extinction coefficient, single 
!                    scattering albedo, assymetry factor) from streamer v3.0,
!                    Key, streamer user guide, cooperative institude 
!                    for meteorological studies, 95 pp., 2001.
! xxxice3 : real   : optical properties (extinction coefficient, single 
!                    scattering albedo, assymetry factor) from
!                    Fu, j. clim., 9, 1996.
! xbari   : real   : optical property coefficients for five spectral 
!                    intervals (2857-4000, 4000-5263, 5263-7692, 7692-14285,
!                    and 14285-40000 wavenumbers) following 
!                    Ebert and Curry, jgr, 97, 3831-3836, 1992.
!------------------------------------------------------------------

      real(kind=rb) :: extliq1(58,16:29), ssaliq1(58,16:29), asyliq1(58,16:29)
      real(kind=rb) :: extice2(43,16:29), ssaice2(43,16:29), asyice2(43,16:29)
      real(kind=rb) :: extice3(46,16:29), ssaice3(46,16:29), asyice3(46,16:29)
      real(kind=rb) :: fdlice3(46,16:29)
      real(kind=rb) :: abari(5),bbari(5),cbari(5),dbari(5),ebari(5),fbari(5)

      end module rrsw_cld

