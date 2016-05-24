      module rrlw_tbl

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw exponential lookup table arrays

! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, Jun 2006
! Revised: MJIacono, AER, Aug 2007
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ntbl   :  integer: Lookup table dimension
! tblint :  real   : Lookup table conversion factor
! tau_tbl:  real   : Clear-sky optical depth (used in cloudy radiative
!                    transfer)
! exp_tbl:  real   : Transmittance lookup table
! tfn_tbl:  real   : Tau transition function; i.e. the transition of
!                    the Planck function from that for the mean layer
!                    temperature to that for the layer boundary
!                    temperature as a function of optical depth.
!                    The "linear in tau" method is used to make 
!                    the table.
! pade   :  real   : Pade constant   
! bpade  :  real   : Inverse of Pade constant   
!------------------------------------------------------------------

      integer, parameter :: ntbl = 10000

      real(kind=r8), parameter :: tblint = 10000.0_r8

      real(kind=r8) , dimension(0:ntbl) :: tau_tbl
      real(kind=r8) , dimension(0:ntbl) :: exp_tbl
      real(kind=r8) , dimension(0:ntbl) :: tfn_tbl

      real(kind=r8), parameter :: pade = 0.278_r8
      real(kind=r8) :: bpade

      end module rrlw_tbl

