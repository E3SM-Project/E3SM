
      module mo_tracname
!-----------------------------------------------------------
! 	... List of advected and non-advected trace species, and
!           surface fluxes for the advected species.
!-----------------------------------------------------------

      use chem_mods, only : grpcnt, gas_pcnst

      implicit none

      character(len=16) :: solsym(gas_pcnst)   ! species names

      end module mo_tracname
