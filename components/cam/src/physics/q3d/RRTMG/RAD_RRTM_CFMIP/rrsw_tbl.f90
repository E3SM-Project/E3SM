      module rrsw_tbl

      use parkind, only : im => kind_im, rb => kind_rb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_sw lookup table arrays

! Initial version: MJIacono, AER, may2007
! Revised: MJIacono, AER, aug2007
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ntbl   :  integer: Lookup table dimension
! tblint :  real   : Lookup table conversion factor
! tau_tbl:  real   : Clear-sky optical depth 
! exp_tbl:  real   : Exponential lookup table for transmittance
! od_lo  :  real   : Value of tau below which expansion is used
!                  : in place of lookup table
! pade   :  real   : Pade approximation constant
! bpade  :  real   : Inverse of Pade constant
!------------------------------------------------------------------

      integer(kind=im), parameter :: ntbl = 10000

      real(kind=rb), parameter :: tblint = 10000.0_rb

      real(kind=rb), parameter :: od_lo = 0.06_rb

      real(kind=rb) :: tau_tbl
      real(kind=rb) , dimension(0:ntbl) :: exp_tbl

      real(kind=rb), parameter :: pade = 0.278_rb
      real(kind=rb) :: bpade

      end module rrsw_tbl

