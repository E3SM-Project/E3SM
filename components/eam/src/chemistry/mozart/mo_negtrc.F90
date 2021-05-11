
      module mo_negtrc

      private
      public :: negtrc

      contains

      subroutine negtrc( header, fld, ncol )
!-----------------------------------------------------------------------
!  	... Check for negative constituent values and
!	    replace with zero value
!-----------------------------------------------------------------------

      use shr_kind_mod, only: r8 => shr_kind_r8
      use chem_mods,   only : gas_pcnst
      use ppgrid,      only : pver

      implicit none

!-----------------------------------------------------------------------
!  	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)          :: ncol
      character(len=*), intent(in) :: header
      real(r8), intent(inout)          :: fld(ncol,pver,gas_pcnst) ! field to check

!-----------------------------------------------------------------------
!  	... Local variables
!-----------------------------------------------------------------------
      integer :: m
      integer :: nneg                       ! flag counter

      do m  = 1,gas_pcnst
         nneg = count( fld(:,:,m) < 0._r8 )
	 if( nneg > 0 ) then
            where( fld(:,:,m) < 0._r8 )
	       fld(:,:,m) = 0._r8
	    endwhere
	 end if
      end do

      end subroutine negtrc

      end module mo_negtrc
