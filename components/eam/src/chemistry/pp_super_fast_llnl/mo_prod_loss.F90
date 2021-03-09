






      module mo_prod_loss

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: exp_prod_loss
      public :: imp_prod_loss

      contains

      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid, only : pver

      implicit none

!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:,:,:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(:,:,:)
      real(r8), intent(in) :: rxt(:,:,:)
      real(r8), intent(in) :: het_rates(:,:,:)



!--------------------------------------------------------------------
! ... loss and production for Explicit method
!--------------------------------------------------------------------


         loss(:,:,1) = (rxt(:,:,16)* y(:,:,2) + het_rates(:,:,8))* y(:,:,8)
         prod(:,:,1) = 0._r8

      end subroutine exp_prod_loss

      subroutine imp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid, only : pver

      implicit none

!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(:)
      real(r8), intent(in) :: rxt(:)
      real(r8), intent(in) :: het_rates(:)



!--------------------------------------------------------------------
! ... loss and production for Implicit method
!--------------------------------------------------------------------


         loss(10) = (rxt(7)* y(2) +rxt(8)* y(3) +rxt(12)* y(5) +rxt(28)* y(13) &
                  +rxt(30)* y(15) + rxt(1) + het_rates(1))* y(1)
         prod(10) =rxt(3)*y(6)
         loss(13) = (rxt(7)* y(1) +rxt(9)* y(3) +rxt(11)* y(4) +rxt(14)* y(6) +rxt(16) &
                 * y(8) +rxt(17)* y(9) +rxt(19)* y(11) + (rxt(24) +rxt(25))* y(12) &
                  +rxt(26)* y(13) +rxt(29)* y(15) + rxt(15) + het_rates(2))* y(2)
         prod(13) = (2.000_r8*rxt(1) +rxt(8)*y(3))*y(1) +rxt(13)*y(5)*y(3) &
                  +2.000_r8*rxt(2)*y(4) +rxt(6)*y(11)
         loss(11) = (rxt(8)* y(1) +rxt(9)* y(2) + 2._r8*rxt(10)* y(3) +rxt(13)* y(5) &
                  +rxt(18)* y(10) + het_rates(3))* y(3)
         prod(11) = (rxt(16)*y(8) +rxt(7)*y(1) +rxt(11)*y(4) +rxt(17)*y(9))*y(2) &
                  + (rxt(21)*y(5) +.800_r8*rxt(22)*y(10))*y(10) +.060_r8*rxt(30)*y(15) &
                 *y(1) +2.000_r8*rxt(4)*y(9) +rxt(6)*y(11)
         loss(4) = (rxt(11)* y(2) +rxt(27)* y(13) + rxt(2) + het_rates(4))* y(4)
         prod(4) =rxt(10)*y(3)*y(3)
         loss(12) = (rxt(12)* y(1) +rxt(13)* y(3) +rxt(21)* y(10) + het_rates(5)) &
                 * y(5)
         prod(12) =rxt(3)*y(6)
         loss(9) = (rxt(14)* y(2) + rxt(3) + rxt(23) + het_rates(6))* y(6)
         prod(9) = (rxt(12)*y(1) +rxt(13)*y(3) +rxt(21)*y(10))*y(5)
         loss(1) = ( + het_rates(7))* y(7)
         prod(1) = (.500_r8*rxt(23) +rxt(14)*y(2))*y(6)
         loss(7) = (rxt(17)* y(2) + rxt(4) + rxt(5) + het_rates(9))* y(9)
         prod(7) = (rxt(21)*y(5) +2.000_r8*rxt(22)*y(10))*y(10) + (rxt(6) + &
                 rxt(20)*y(2))*y(11) +.870_r8*rxt(30)*y(15)*y(1)
         loss(14) = (rxt(18)* y(3) +rxt(21)* y(5) + 2._r8*rxt(22)* y(10) &
                  + het_rates(10))* y(10)
         prod(14) = (rxt(15) +rxt(19)*y(11))*y(2) +1.860_r8*rxt(30)*y(15)*y(1)
         loss(8) = ((rxt(19) +rxt(20))* y(2) + rxt(6) + het_rates(11))* y(11)
         prod(8) =rxt(18)*y(10)*y(3)
         loss(3) = ((rxt(24) +rxt(25))* y(2) + het_rates(12))* y(12)
         prod(3) = 0._r8
         loss(5) = (rxt(28)* y(1) +rxt(26)* y(2) +rxt(27)* y(4) + het_rates(13)) &
                 * y(13)
         prod(5) = (rxt(24)*y(12) +.750_r8*rxt(25)*y(12))*y(2)
         loss(2) = ( + het_rates(14))* y(14)
         prod(2) = (rxt(26)*y(2) +rxt(27)*y(4) +rxt(28)*y(1))*y(13)
         loss(6) = (rxt(30)* y(1) +rxt(29)* y(2) + het_rates(15))* y(15)
         prod(6) = 0._r8

      end subroutine imp_prod_loss

      end module mo_prod_loss
