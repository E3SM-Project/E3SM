      module rrsw_kg26

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb
      use parrrsw, only : ng26

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 26
! band 26: 22650-29000 cm-1 (low - nothing; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!sfluxrefo: real     
! raylo   : real     
!-----------------------------------------------------------------

      integer, parameter :: no26 = 16

      real(kind=r8) :: sfluxrefo(no26)
      real(kind=r8) :: raylo(no26)

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 26
! band 26: 22650-29000 cm-1 (low - nothing; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! sfluxref: real     
! rayl    : real     
!-----------------------------------------------------------------

      real(kind=r8) :: sfluxref(ng26)
      real(kind=r8) :: rayl(ng26)

      end module rrsw_kg26

