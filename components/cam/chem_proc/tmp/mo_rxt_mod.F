
      module mo_rxt_mod

      private
      public :: rxt_mod

      contains

      subroutine rxt_mod( rate, het_rates, grp_ratios, chnkpnts )

      use chem_mods, only : rxntot, hetcnt, grpcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!---------------------------------------------------------------------------
!       ... dummy arguments
!---------------------------------------------------------------------------
      integer, intent(in) ::  plnplv
      real(r8), intent(inout) ::  rate(:,:)
      real(r8), intent(inout) ::  het_rates(:,:)
      real(r8), intent(in)    ::  grp_ratios(:,:)


      end subroutine rxt_mod

      end module mo_rxt_mod
