
      module mo_make_grp_vmr

      private
      public :: mak_grp_vmr

      contains

      subroutine mak_grp_vmr( vmr, group_ratios, group_vmrs, plonl )

      use chem_mods, only : grpcnt, gas_pcnst
      use ppgrid, only : pver
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------------------------------------
!        ... dummy arguments
!----------------------------------------------------------------------------
      integer, intent(in) :: plonl
      real(r8), intent(in)    :: vmr(:,:,:)
      real(r8), intent(in)    :: group_ratios(:,:,:)
      real(r8), intent(out)   :: group_vmrs(:,:,:)

!----------------------------------------------------------------------------
!        ... local variables
!----------------------------------------------------------------------------
      integer ::  k

      end subroutine mak_grp_vmr

      end module mo_make_grp_vmr
